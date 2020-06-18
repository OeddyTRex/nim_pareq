c-----------------------------------------------------------------------
c     file nimplot_mgt.f:  contains management routines for finite
c     element computations performed for diagnostic calculations
c     in nimplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  jfromb.
c     2.  energies.
c     3.  divb.
c     4.  parallel_current.
c     5.  e_dot_j.
c     6.  heat_flux.
c     7.  divb_element.
c     8.  j_element.
c     9.  dyn_element.
c     10. mach.
c-----------------------------------------------------------------------
c     subprogram 1. jfromb.
c     find total current density as a bilinear quantity from magnetic
c     field.  a lumped mass matrix is used to find j as a vertex
c     quantity, and this may distort to the local curl of B, which
c     is what's usually used throughout nimrod.
c-----------------------------------------------------------------------
      SUBROUTINE jfromb(mode_st,mode_en)
      USE local
      USE input
      USE physdat
      USE nimplot_ints
      USE rblock
      USE tblock
      USE computation_pointers
      USE global
      USE fields
      USE edge
      USE iter_cg
      USE matrix_storage_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: mode_st,mode_en

      INTEGER(i4) :: ibl,iq,im,its
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     copy B into the work1 array for the curl integrand.
c     always copy the entire array--this unmolested B is used for
c     finding the poloidal flux function in xy vs. pol flux slices.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        work1(ibl)=be(ibl)
        ja(ibl)=0
      ENDDO
c-----------------------------------------------------------------------
c     find integral(alpha*curl(B)).
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
        CALL rblock_get_rhs(rb(ibl),crhs(ibl),curl,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         3_i4,nmodes)
        CALL tblock_get_rhs(tb(ibl),crhs(ibl),curl,3_i4,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     network seams that connect block borders.
c-----------------------------------------------------------------------
      CALL edge_network(3_i4,nmodes,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix for each component.
c-----------------------------------------------------------------------
      DO im=mode_st,mode_en
        DO iq=1,3
          DO ibl=1,nbl
            csln(ibl)=0
            cvectr(ibl)%arr(1,:,:)=crhs(ibl)%arr(iq,:,:,im)
            IF (ASSOCIATED(csln(ibl)%arrh).AND.ibl<=nrbl) THEN
              cvectr(ibl)%arrh(1,:,:,:)=crhs(ibl)%arrh(iq,:,:,:,im)
              cvectr(ibl)%arrv(1,:,:,:)=crhs(ibl)%arrv(iq,:,:,:,im)
              cvectr(ibl)%arri(1,:,:,:)=crhs(ibl)%arri(iq,:,:,:,im)
            ENDIF
          ENDDO
          CALL iter_cg_2d_solve(cmass_mat,cmass_fac,csln,cvectr,1_i4,
     $                          tol,maxit,solver,err,its,seed)
          IF (err>tol) THEN
            WRITE(msg,'(a,i4,a,es10.3,3a)') 'Jfromb: no convergence: ',
     $        its,' its ',err,' err ',seed,' seed'
            CALL nim_stop(msg)
          ENDIF
          DO ibl=1,nbl
            ja(ibl)%arr(iq,:,:,im)=csln(ibl)%arr(1,:,:)/mu0
            IF (ASSOCIATED(csln(ibl)%arrh).AND.ibl<=nrbl) THEN
              ja(ibl)%arrh(iq,:,:,:,im)=csln(ibl)%arrh(1,:,:,:)/mu0
              ja(ibl)%arrv(iq,:,:,:,im)=csln(ibl)%arrv(1,:,:,:)/mu0
              ja(ibl)%arri(iq,:,:,:,im)=csln(ibl)%arri(1,:,:,:)/mu0
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jfromb
c-----------------------------------------------------------------------
c     subprogram 2. energies.
c     modified from the nimrod history.f file.
c-----------------------------------------------------------------------
      SUBROUTINE energies(first_data)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data

      INTEGER :: read_stat,num_chars
      INTEGER(i4) :: ibl,ierror,imode,iq
      REAL(r8), DIMENSION(:), ALLOCATABLE :: kinetic_energy,
     $          magnetic_energy
      REAL(r8) :: total_energy,internal_energy
      CHARACTER(6) :: hist_pos,ctmp
      CHARACTER(6), SAVE :: ansl
c-----------------------------------------------------------------------
c     allocate arrays if this is the first call to energies.
c-----------------------------------------------------------------------
      ALLOCATE(magnetic_energy(nmodes))
      ALLOCATE(kinetic_energy(nmodes))
      magnetic_energy=0
      internal_energy=0
      kinetic_energy=0
c-----------------------------------------------------------------------
c     accumulate integral[b**2/2mu0] and integral[ro*v**2/2] for each
c     mode and integral[p/(gamma-1)] for n=0.  for nonlinear runs, this
c     includes the equilibrium fields, but for linear runs it doesn't.
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          CALL vector_type_alloc(cell_crhs(ibl),0_i4,rb(ibl)%mx,
     $                           rb(ibl)%my,4_i4,nmodes)
          CALL rblock_get_rhs(rb(ibl),cell_crhs(ibl),energy_density,
     $                        4_i4,nmodes)
        ELSE
          CALL vector_type_alloc(cell_crhs(ibl),0_i4,tb(ibl)%mcell,1_i4,
     $                           4_i4,nmodes)
          CALL tblock_get_rhs(tb(ibl),cell_crhs(ibl),energy_density,
     $                        4_i4,nmodes)
        ENDIF
        DO imode=1,nmodes
          magnetic_energy(imode)=magnetic_energy(imode)
     $       +SUM(cell_crhs(ibl)%arri(1,1,:,:,imode))
          kinetic_energy(imode)=kinetic_energy(imode)
     $       +SUM(cell_crhs(ibl)%arri(2,1,:,:,imode))
          internal_energy=internal_energy
     $       +SUM(cell_crhs(ibl)%arri(3,1,:,:,imode))
     $       +SUM(cell_crhs(ibl)%arri(4,1,:,:,imode))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     total energy:
c-----------------------------------------------------------------------
      total_energy=SUM(magnetic_energy+kinetic_energy)+internal_energy
c-----------------------------------------------------------------------
c     geometric factor:
c-----------------------------------------------------------------------
      IF (geom=='lin') THEN
        magnetic_energy=magnetic_energy*per_length
        kinetic_energy=kinetic_energy*per_length
        internal_energy=internal_energy*per_length
        total_energy=total_energy*per_length
      ELSE
        magnetic_energy=magnetic_energy*twopi
        kinetic_energy=kinetic_energy*twopi
        internal_energy=internal_energy*twopi
        total_energy=total_energy*twopi
      ENDIF
c-----------------------------------------------------------------------
c     plot log if desired.
c-----------------------------------------------------------------------
      IF (first_data) THEN
        ansl='lin'
        WRITE(nim_wr,'(/5a)') '>>> Hit return for linear plots, enter ',
     $                   "'",'ln',"'",' for'
        WRITE(nim_wr,'(5a)') '>>> natural logs, or enter ',
     $                  "'",'log',"'",' for log base 10.'
        WRITE(nim_wr,'(a,$)') '>>>? '
        READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $       SIZE=num_chars) ctmp
        IF (num_chars>0) ansl=ctmp
      ENDIF
      IF (ansl=='ln') THEN
        magnetic_energy=LOG(MAX(magnetic_energy,TINY(dtm)))
        kinetic_energy=LOG(MAX(kinetic_energy,TINY(dtm)))
      ELSE IF (ansl=='log') THEN
        magnetic_energy=LOG10(MAX(magnetic_energy,TINY(dtm)))
        kinetic_energy=LOG10(MAX(kinetic_energy,TINY(dtm)))
      ENDIF
c-----------------------------------------------------------------------
c     output:
c-----------------------------------------------------------------------
      IF (first_data) THEN
        hist_pos='REWIND'
      ELSE
        hist_pos='APPEND'
      ENDIF
      CALL open_bin(en_unit,'energy.bin','UNKNOWN',hist_pos,32_i4)
c-----------------------------------------------------------------------
c     array constructors are needed to avoid a write error on the c90
c     and j90.
c-----------------------------------------------------------------------
      DO iq=1,nmodes
        WRITE(en_unit) (/ REAL(istep,4),REAL(t,4),REAL(iq,4),
     $                    REAL(keff(iq),4),
     $                    REAL(magnetic_energy(iq),4),
     $                    REAL(kinetic_energy (iq),4) /)
      ENDDO
      WRITE(en_unit)
      CALL close_bin(en_unit,'energy.bin' )
      DEALLOCATE(magnetic_energy,kinetic_energy)
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(cell_crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE energies
c-----------------------------------------------------------------------
c     subprogram 3. divb.
c     find divb as a vertex quantity.
c-----------------------------------------------------------------------
      SUBROUTINE divb(first_data,datflag,mode_st,mode_en,phii,ndcon)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en,ndcon
      REAL(r8), INTENT(IN) :: phii
      TYPE(vector_type), DIMENSION(nbl) :: divbn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     find integral[alpha*div(b)]
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
        CALL rblock_get_rhs(rb(ibl),crhs(ibl),divb_int,1_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         1_i4,nmodes)
        CALL tblock_get_rhs(tb(ibl),crhs(ibl),divb_int,1_i4,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     network block seams.
c-----------------------------------------------------------------------
      CALL edge_network(1_i4,nmodes,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix to find div(b) at the vertices.
c-----------------------------------------------------------------------
      DO im=1,nmodes
        DO ibl=1,nbl
          csln(ibl)=0
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),im,1_i4)
        ENDDO
        CALL iter_cg_2d_solve(cmass_mat,cmass_fac,csln,cvectr,1_i4,
     $                        tol,maxit,solver,err,its,seed)
        IF (err>tol) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Divb: no convergence: ',
     $      its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        DO ibl=1,nbl
          CALL cvector_assign_cvec2(work3(ibl),csln(ibl),im)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     initialize contour plot if needed, otherwise just open the file.
c     option c plots div(b) at a specified angle in configuration space.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        nq=2*(mode_en-mode_st+1)
        iqst=mode_st
        iqen=mode_en
      ELSE
        nq=2
      ENDIF
      IF (first_data) THEN
        CALL contour_init(nq,ndcon,"divb.bin")
      ELSE
        CALL open_bin(con_unit,"divb.bin","OLD","APPEND",32_i4)
      ENDIF
c-----------------------------------------------------------------------
c     write the appropriate data.  save and restore the n=0 part if
c     plotting at angles.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        DO ibl=1,nrbl
          CALL contour_laq_write(rb(ibl)%work3,iqst,iqen,ndcon)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL contour_tl_write(tb(ibl)%work3,iqst,iqen)
        ENDDO
      ELSE
        jphi=0
        DO ibl=1,nbl
          mxb=SIZE(work3(ibl)%arr,2)-1
          myb=SIZE(work3(ibl)%arr,3)-1
          CALL vector_type_alloc(divbn0(ibl),poly_degree,mxb,myb,1_i4)
          CALL vector_assign_cvec(divbn0(ibl),work3(ibl),'real',1_i4)
        ENDDO
        phi_loop: DO
          IF (datflag=='c'.AND.jphi/=1) THEN
            jphi=jphi+1
            IF (jphi==1) CYCLE phi_loop
            EXIT phi_loop
          ENDIF
          CALL eval_at_angle(work3,MODULO(MIN(jphi*phii,1._r8),1._r8))
          DO ibl=1,nrbl
            CALL contour_laq_write(rb(ibl)%work3,1_i4,1_i4,ndcon)
            CALL cvector_assign_vec(work3(ibl),divbn0(ibl),'real',1_i4)
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL contour_tl_write(tb(ibl)%work3,1_i4,1_i4)
            CALL cvector_assign_vec(work3(ibl),divbn0(ibl),'real',1_i4)
          ENDDO
          IF (jphi*phii>=1) EXIT phi_loop
          jphi=jphi+1
        ENDDO phi_loop
        DO ibl=1,nbl
          CALL vector_type_dealloc(divbn0(ibl))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close the plot file.
c-----------------------------------------------------------------------
      CALL close_bin(con_unit,"divb.bin")
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE divb
c-----------------------------------------------------------------------
c     subprogram 4. parallel_current.
c     find mu0*J.B/B**2 as a vertex quantity.
c-----------------------------------------------------------------------
      SUBROUTINE parallel_current(first_data,datflag,mode_st,mode_en,
     $                            phii,ndcon)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en,ndcon
      REAL(r8), INTENT(IN) :: phii
      TYPE(vector_type), DIMENSION(nbl) :: lamn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     find integral[alpha*mu0*J.B/B**2]
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
        CALL rblock_get_rhs(rb(ibl),crhs(ibl),jpar_int,1_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         1_i4,nmodes)
        CALL tblock_get_rhs(tb(ibl),crhs(ibl),jpar_int,1_i4,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     network block seams.
c-----------------------------------------------------------------------
      CALL edge_network(1_i4,nmodes,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix to find J/B at the vertices.
c-----------------------------------------------------------------------
      DO iq=1,nmodes
        DO ibl=1,nbl
          csln(ibl)=0
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),iq,1_i4)
        ENDDO
        CALL iter_cg_2d_solve(cmass_mat,cmass_fac,csln,cvectr,1_i4,
     $                        tol,maxit,solver,err,its,seed)
        IF (err>tol) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Parallel_current: no '
     $      //'convergence: ',its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        DO ibl=1,nbl
          CALL cvector_assign_cvec2(work3(ibl),csln(ibl),iq)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     initialize contour plot if needed, otherwise just open the file.
c     option c plots Jpar/B at a specified angle in configuration space.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        nq=2*(mode_en-mode_st+1)
        iqst=mode_st
        iqen=mode_en
      ELSE
        nq=2
      ENDIF
      IF (first_data) THEN
        CALL contour_init(nq,ndcon,"parallel_current.bin")
      ELSE
        CALL open_bin(con_unit,"parallel_current.bin","OLD","APPEND",
     $                32_i4)
      ENDIF
c-----------------------------------------------------------------------
c     write the appropriate data.  save and restore the n=0 part if
c     plotting at angles.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        DO ibl=1,nrbl
          CALL contour_laq_write(rb(ibl)%work3,iqst,iqen,ndcon)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL contour_tl_write(tb(ibl)%work3,iqst,iqen)
        ENDDO
      ELSE
        jphi=0
        DO ibl=1,nbl
          mxb=SIZE(work3(ibl)%arr,2)-1
          myb=SIZE(work3(ibl)%arr,3)-1
          CALL vector_type_alloc(lamn0(ibl),poly_degree,mxb,myb,1_i4)
          CALL vector_assign_cvec(lamn0(ibl),work3(ibl),'real',1_i4)
        ENDDO
        phi_loop: DO
          IF (datflag=='c'.AND.jphi/=1) THEN
            jphi=jphi+1
            IF (jphi==1) CYCLE phi_loop
            EXIT phi_loop
          ENDIF
          CALL eval_at_angle(work3,MODULO(MIN(jphi*phii,1._r8),1._r8))
          DO ibl=1,nrbl
            CALL contour_laq_write(rb(ibl)%work3,1_i4,1_i4,ndcon)
            CALL cvector_assign_vec(work3(ibl),lamn0(ibl),'real',1_i4)
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL contour_tl_write(tb(ibl)%work3,1_i4,1_i4)
            CALL cvector_assign_vec(work3(ibl),lamn0(ibl),'real',1_i4)
          ENDDO
          IF (jphi*phii>=1) EXIT phi_loop
          jphi=jphi+1
        ENDDO phi_loop
        DO ibl=1,nbl
          CALL vector_type_dealloc(lamn0(ibl))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close the plot file.
c-----------------------------------------------------------------------
      CALL close_bin(con_unit,"parallel_current.bin")
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE parallel_current
c-----------------------------------------------------------------------
c     subprogram 5. e_dot_j.
c     find <E>.<J>, the electromagnetic power transferred from/to the
c     toroidally averaged current density.  <E> is broken into different
c     physical constituents.
c-----------------------------------------------------------------------
      SUBROUTINE e_dot_j(first_data,flag,mode,ndcon)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: flag
      INTEGER(i4), INTENT(IN) :: mode,ndcon

      INTEGER(i4) :: ibl,iq,jq,its
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     find integral[alpha*<-vXb>.<J>], and integral[alpha*eta*J**2]
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,4_i4*nmodes+3_i4)
        CALL rblock_get_rhs(rb(ibl),rhs(ibl),edotj_int,
     $                      4_i4*nmodes+3_i4)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         4_i4*nmodes+3_i4)
        CALL tblock_get_rhs(tb(ibl),rhs(ibl),edotj_int,
     $                      4_i4*nmodes+3_i4)
      ENDDO
c-----------------------------------------------------------------------
c     network block seams.
c-----------------------------------------------------------------------
      CALL edge_network(4_i4*nmodes+3_i4,0_i4,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix to find all E.J products at the vertices.
c-----------------------------------------------------------------------
      DO iq=1,4*nmodes+3
        DO ibl=1,nbl
          sln(ibl)=0
          vectr(ibl)%arr(1,:,:)=rhs(ibl)%arr(iq,:,:)
          IF (ASSOCIATED(rhs(ibl)%arrv).AND.ibl<=nrbl) THEN
            vectr(ibl)%arrh(1,:,:,:)=rhs(ibl)%arrh(iq,:,:,:)
            vectr(ibl)%arrv(1,:,:,:)=rhs(ibl)%arrv(iq,:,:,:)
            vectr(ibl)%arri(1,:,:,:)=rhs(ibl)%arri(iq,:,:,:)
          ENDIF
        ENDDO
        CALL iter_cg_2d_solve(mass_mat,mass_fac,sln,vectr,1_i4,
     $                        tol,maxit,solver,err,its,seed)
        IF (err>tol) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'E_dot_j: no convergence: ',
     $      its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        DO ibl=1,nbl
          rhs(ibl)%arr(iq,:,:)=sln(ibl)%arr(1,:,:)
          IF (ASSOCIATED(rhs(ibl)%arrh).AND.ibl<=nrbl) THEN
            rhs(ibl)%arrh(iq,:,:,:)=sln(ibl)%arrh(1,:,:,:)
            rhs(ibl)%arrv(iq,:,:,:)=sln(ibl)%arrv(1,:,:,:)
            rhs(ibl)%arri(iq,:,:,:)=sln(ibl)%arri(1,:,:,:)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     select the appropriate data.
c     option a plots the correlations as the sum over all Fourier
c     components.
c-----------------------------------------------------------------------
      IF (flag=='a') THEN
        jq=1
        DO ibl=1,nbl
          DO iq=2,nmodes
            rhs(ibl)%arr(1,:,:)=rhs(ibl)%arr(1 ,:,:)
     $                         +rhs(ibl)%arr(iq,:,:)
            rhs(ibl)%arr(  nmodes+1,:,:)=rhs(ibl)%arr(  nmodes+1 ,:,:)
     $                                  +rhs(ibl)%arr(  nmodes+iq,:,:)
            rhs(ibl)%arr(2*nmodes+1,:,:)=rhs(ibl)%arr(2*nmodes+1 ,:,:)
     $                                  +rhs(ibl)%arr(2*nmodes+iq,:,:)
            rhs(ibl)%arr(3*nmodes+1,:,:)=rhs(ibl)%arr(3*nmodes+1 ,:,:)
     $                                  +rhs(ibl)%arr(3*nmodes+iq,:,:)
          ENDDO
          IF (ASSOCIATED(rhs(ibl)%arrh)) THEN
            DO iq=2,nmodes
              rhs(ibl)%arrh(1,:,:,:)=rhs(ibl)%arrh(1 ,:,:,:)
     $                              +rhs(ibl)%arrh(iq,:,:,:)
              rhs(ibl)%arrh(  nmodes+1,:,:,:)=
     $           rhs(ibl)%arrh(  nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrh(  nmodes+iq,:,:,:)
              rhs(ibl)%arrh(2*nmodes+1,:,:,:)=
     $           rhs(ibl)%arrh(2*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrh(2*nmodes+iq,:,:,:)
              rhs(ibl)%arrh(3*nmodes+1,:,:,:)=
     $           rhs(ibl)%arrh(3*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrh(3*nmodes+iq,:,:,:)
            ENDDO
          ENDIF
          IF (ASSOCIATED(rhs(ibl)%arrv)) THEN
            DO iq=2,nmodes
              rhs(ibl)%arrv(1,:,:,:)=rhs(ibl)%arrv(1 ,:,:,:)
     $                              +rhs(ibl)%arrv(iq,:,:,:)
              rhs(ibl)%arrv(  nmodes+1,:,:,:)=
     $           rhs(ibl)%arrv(  nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrv(  nmodes+iq,:,:,:)
              rhs(ibl)%arrv(2*nmodes+1,:,:,:)=
     $           rhs(ibl)%arrv(2*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrv(2*nmodes+iq,:,:,:)
              rhs(ibl)%arrv(3*nmodes+1,:,:,:)=
     $           rhs(ibl)%arrv(3*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arrv(3*nmodes+iq,:,:,:)
            ENDDO
          ENDIF
          IF (ASSOCIATED(rhs(ibl)%arri)) THEN
            DO iq=2,nmodes
              rhs(ibl)%arri(1,:,:,:)=rhs(ibl)%arri(1 ,:,:,:)
     $                              +rhs(ibl)%arri(iq,:,:,:)
              rhs(ibl)%arri(  nmodes+1,:,:,:)=
     $           rhs(ibl)%arri(  nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arri(  nmodes+iq,:,:,:)
              rhs(ibl)%arri(2*nmodes+1,:,:,:)=
     $           rhs(ibl)%arri(2*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arri(2*nmodes+iq,:,:,:)
              rhs(ibl)%arri(3*nmodes+1,:,:,:)=
     $           rhs(ibl)%arri(3*nmodes+1 ,:,:,:)
     $          +rhs(ibl)%arri(3*nmodes+iq,:,:,:)
            ENDDO
          ENDIF
        ENDDO
      ELSE
        jq=mode
      ENDIF
c-----------------------------------------------------------------------
c     initialize contour plot if needed, otherwise just open the file.
c-----------------------------------------------------------------------
      IF (first_data) THEN
        CALL contour_init(10_i4,ndcon,"e_dot_j.bin")
      ELSE
        CALL open_bin(con_unit,"e_dot_j.bin","OLD","APPEND",32_i4)
      ENDIF
c-----------------------------------------------------------------------
c     transfer the data to work2, so standard contour write routines
c     may be used.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        work2(ibl)%arr(1,:,:,1)=rhs(ibl)%arr(         jq,:,:)
     $                         +rhs(ibl)%arr(  nmodes+jq,:,:)
        work2(ibl)%arr(2,:,:,1)=rhs(ibl)%arr(         jq,:,:)
        work2(ibl)%arr(3,:,:,1)=rhs(ibl)%arr(  nmodes+jq,:,:)
        work2(ibl)%arr(4,:,:,1)=rhs(ibl)%arr(2*nmodes+jq,:,:)
     $                         +rhs(ibl)%arr(3*nmodes+jq,:,:)
        work2(ibl)%arr(5,:,:,1)=rhs(ibl)%arr(2*nmodes+jq,:,:)
        work2(ibl)%arr(1,:,:,1)=work2(ibl)%arr(1,:,:,1)+
     $                    (0,1)*rhs(ibl)%arr(3*nmodes+jq,:,:)
        work2(ibl)%arr(2,:,:,1)=work2(ibl)%arr(2,:,:,1)+
     $                    (0,1)*(rhs(ibl)%arr(4*nmodes+1 ,:,:)
     $                          +rhs(ibl)%arr(4*nmodes+2 ,:,:))
        work2(ibl)%arr(3,:,:,1)=work2(ibl)%arr(3,:,:,1)+
     $                    (0,1)*rhs(ibl)%arr(4*nmodes+1 ,:,:)
        work2(ibl)%arr(4,:,:,1)=work2(ibl)%arr(4,:,:,1)+
     $                    (0,1)*rhs(ibl)%arr(4*nmodes+2 ,:,:)
        work2(ibl)%arr(5,:,:,1)=work2(ibl)%arr(5,:,:,1)+
     $                    (0,1)*rhs(ibl)%arr(4*nmodes+3 ,:,:)
        IF (ASSOCIATED(rhs(ibl)%arrh).AND.ibl<=nrbl) THEN
          work2(ibl)%arrh(1,:,:,:,1)=rhs(ibl)%arrh(         jq,:,:,:)
     $                              +rhs(ibl)%arrh(  nmodes+jq,:,:,:)
          work2(ibl)%arrh(2,:,:,:,1)=rhs(ibl)%arrh(         jq,:,:,:)
          work2(ibl)%arrh(3,:,:,:,1)=rhs(ibl)%arrh(  nmodes+jq,:,:,:)
          work2(ibl)%arrh(4,:,:,:,1)=rhs(ibl)%arrh(2*nmodes+jq,:,:,:)
     $                              +rhs(ibl)%arrh(3*nmodes+jq,:,:,:)
          work2(ibl)%arrh(5,:,:,:,1)=rhs(ibl)%arrh(2*nmodes+jq,:,:,:)
          work2(ibl)%arrh(1,:,:,:,1)=work2(ibl)%arrh(1,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrh(3*nmodes+jq,:,:,:)
          work2(ibl)%arrh(2,:,:,:,1)=work2(ibl)%arrh(2,:,:,:,1)+
     $                         (0,1)*(rhs(ibl)%arrh(4*nmodes+1,:,:,:)
     $                               +rhs(ibl)%arrh(4*nmodes+2,:,:,:))
          work2(ibl)%arrh(3,:,:,:,1)=work2(ibl)%arrh(3,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrh(4*nmodes+1,:,:,:)
          work2(ibl)%arrh(4,:,:,:,1)=work2(ibl)%arrh(4,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrh(4*nmodes+2,:,:,:)
          work2(ibl)%arrh(5,:,:,:,1)=work2(ibl)%arrh(5,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrh(4*nmodes+3,:,:,:)
          work2(ibl)%arrv(1,:,:,:,1)=rhs(ibl)%arrv(         jq,:,:,:)
     $                              +rhs(ibl)%arrv(  nmodes+jq,:,:,:)
          work2(ibl)%arrv(2,:,:,:,1)=rhs(ibl)%arrv(         jq,:,:,:)
          work2(ibl)%arrv(3,:,:,:,1)=rhs(ibl)%arrv(  nmodes+jq,:,:,:)
          work2(ibl)%arrv(4,:,:,:,1)=rhs(ibl)%arrv(2*nmodes+jq,:,:,:)
     $                              +rhs(ibl)%arrv(3*nmodes+jq,:,:,:)
          work2(ibl)%arrv(5,:,:,:,1)=rhs(ibl)%arrv(2*nmodes+jq,:,:,:)
          work2(ibl)%arrv(1,:,:,:,1)=work2(ibl)%arrv(1,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrv(3*nmodes+jq,:,:,:)
          work2(ibl)%arrv(2,:,:,:,1)=work2(ibl)%arrv(2,:,:,:,1)+
     $                         (0,1)*(rhs(ibl)%arrv(4*nmodes+1,:,:,:)
     $                               +rhs(ibl)%arrv(4*nmodes+2,:,:,:))
          work2(ibl)%arrv(3,:,:,:,1)=work2(ibl)%arrv(3,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrv(4*nmodes+1,:,:,:)
          work2(ibl)%arrv(4,:,:,:,1)=work2(ibl)%arrv(4,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrv(4*nmodes+2,:,:,:)
          work2(ibl)%arrv(5,:,:,:,1)=work2(ibl)%arrv(5,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arrv(4*nmodes+3,:,:,:)
          work2(ibl)%arri(1,:,:,:,1)=rhs(ibl)%arri(         jq,:,:,:)
     $                              +rhs(ibl)%arri(  nmodes+jq,:,:,:)
          work2(ibl)%arri(2,:,:,:,1)=rhs(ibl)%arri(         jq,:,:,:)
          work2(ibl)%arri(3,:,:,:,1)=rhs(ibl)%arri(  nmodes+jq,:,:,:)
          work2(ibl)%arri(4,:,:,:,1)=rhs(ibl)%arri(2*nmodes+jq,:,:,:)
     $                              +rhs(ibl)%arri(3*nmodes+jq,:,:,:)
          work2(ibl)%arri(5,:,:,:,1)=rhs(ibl)%arri(2*nmodes+jq,:,:,:)
          work2(ibl)%arri(1,:,:,:,1)=work2(ibl)%arri(1,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arri(3*nmodes+jq,:,:,:)
          work2(ibl)%arri(2,:,:,:,1)=work2(ibl)%arri(2,:,:,:,1)+
     $                         (0,1)*(rhs(ibl)%arri(4*nmodes+1,:,:,:)
     $                               +rhs(ibl)%arri(4*nmodes+2,:,:,:))
          work2(ibl)%arri(3,:,:,:,1)=work2(ibl)%arri(3,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arri(4*nmodes+1,:,:,:)
          work2(ibl)%arri(4,:,:,:,1)=work2(ibl)%arri(4,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arri(4*nmodes+2,:,:,:)
          work2(ibl)%arri(5,:,:,:,1)=work2(ibl)%arri(5,:,:,:,1)+
     $                         (0,1)*rhs(ibl)%arri(4*nmodes+3,:,:,:)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     plot the data.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL contour_laq_write(rb(ibl)%work2,1_i4,1_i4,ndcon)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL contour_tl_write(tb(ibl)%work2,1_i4,1_i4)
      ENDDO
      CALL close_bin(con_unit,"e_dot_j.bin")
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE e_dot_j
c-----------------------------------------------------------------------
c     subprogram 6. heat_flux.
c     compute and plot the conductive heat flux.
c-----------------------------------------------------------------------
      SUBROUTINE heat_flux(first_data,flag,mode,ndcon)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: flag
      INTEGER(i4), INTENT(IN) :: mode,ndcon

      INTEGER(i4) :: ibl,iq,jq,its
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     if using all Fourier components, set jmode to -1, otherwise set
c     it to the selected component to communicate this to the integrand
c     routine.
c-----------------------------------------------------------------------
      IF (flag=='a') THEN
        jmode=-1
      ELSE
        jmode=mode
      ENDIF
c-----------------------------------------------------------------------
c     find integral[alpha*q_parallel], and integral[alpha*q_perp] using
c     the k_pll and k_perp diffusivities, regardless of what has been
c     used in the simulation.
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,10_i4)
        CALL rblock_get_rhs(rb(ibl),rhs(ibl),heat_flux_int,10_i4)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,10_i4)
        CALL tblock_get_rhs(tb(ibl),rhs(ibl),heat_flux_int,10_i4)
      ENDDO
c-----------------------------------------------------------------------
c     network block seams.
c-----------------------------------------------------------------------
      CALL edge_network(10_i4,0_i4,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix to find heat flux vectors at the nodes.
c     load results into the real and imaginary parts of work2 to
c     facilitate plotting.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        work2(ibl)=0
      ENDDO
      DO iq=1,10
        DO ibl=1,nbl
          sln(ibl)=0
          vectr(ibl)%arr(1,:,:)=rhs(ibl)%arr(iq,:,:)
          IF (ASSOCIATED(rhs(ibl)%arrv).AND.ibl<=nrbl) THEN
            vectr(ibl)%arrh(1,:,:,:)=rhs(ibl)%arrh(iq,:,:,:)
            vectr(ibl)%arrv(1,:,:,:)=rhs(ibl)%arrv(iq,:,:,:)
            vectr(ibl)%arri(1,:,:,:)=rhs(ibl)%arri(iq,:,:,:)
          ENDIF
        ENDDO
        CALL iter_cg_2d_solve(mass_mat,mass_fac,sln,vectr,1_i4,
     $                        tol,maxit,solver,err,its,seed)
        IF (err>tol) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Heat_flux: no convergence: ',
     $      its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        IF (iq<=5) THEN
          DO ibl=1,nbl
            work2(ibl)%arr(iq,:,:,1)=sln(ibl)%arr(1,:,:)
            IF (ASSOCIATED(rhs(ibl)%arrh).AND.ibl<=nrbl) THEN
              work2(ibl)%arrh(iq,:,:,:,1)=sln(ibl)%arrh(1,:,:,:)
              work2(ibl)%arrv(iq,:,:,:,1)=sln(ibl)%arrv(1,:,:,:)
              work2(ibl)%arri(iq,:,:,:,1)=sln(ibl)%arri(1,:,:,:)
            ENDIF
          ENDDO
        ELSE
          jq=iq-5
          DO ibl=1,nbl
            work2(ibl)%arr(jq,:,:,1)=work2(ibl)%arr(jq,:,:,1)
     $                          +(0,1)*sln(ibl)%arr(1,:,:)
            IF (ASSOCIATED(rhs(ibl)%arrh).AND.ibl<=nrbl) THEN
              work2(ibl)%arrh(jq,:,:,:,1)=work2(ibl)%arrh(jq,:,:,:,1)
     $                               +(0,1)*sln(ibl)%arrh(1,:,:,:)
              work2(ibl)%arrv(jq,:,:,:,1)=work2(ibl)%arrv(jq,:,:,:,1)
     $                               +(0,1)*sln(ibl)%arrv(1,:,:,:)
              work2(ibl)%arri(jq,:,:,:,1)=work2(ibl)%arri(jq,:,:,:,1)
     $                               +(0,1)*sln(ibl)%arri(1,:,:,:)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     initialize contour plot if needed, otherwise just open the file.
c-----------------------------------------------------------------------
      IF (first_data) THEN
        CALL contour_init(10_i4,ndcon,"heat_flux.bin")
      ELSE
        CALL open_bin(con_unit,"heat_flux.bin","OLD","APPEND",32_i4)
      ENDIF
c-----------------------------------------------------------------------
c     plot the data.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL contour_laq_write(rb(ibl)%work2,1_i4,1_i4,ndcon,
     $                         nqlim=5_i4)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL contour_tl_write(tb(ibl)%work2,1_i4,1_i4,nqlim=5_i4)
      ENDDO
      CALL close_bin(con_unit,"heat_flux.bin")
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE heat_flux
c-----------------------------------------------------------------------
c     subprogram 7. divb_element.
c     find divb within elements and create a file that can
c     be used to create contour plots for fields that are discontinuous
c     at element boundaries.
c-----------------------------------------------------------------------
      SUBROUTINE divb_element(datflag,mode_st,mode_en)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      USE eldata_type_mod
      USE tecplot_mod
      USE vtk_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en
      TYPE(vector_type), DIMENSION(nbl) :: divbn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb, 
     $               ix,iy,ib,iby,ibx,n1,n2,n3,n4,ipt,iel,idat,ifc
      REAL(r8) :: err,dx,dy,x,y,eps=1.e-6
      REAL(r8) :: drdx,drdy,dzdx,dzdy
      REAL(r8) :: jac,dxdr,dxdz,dydr,dydz
      COMPLEX(r8) :: dbrdx,dbrdy,dbzdx,dbzdy,bp,br,divb
      CHARACTER(1) :: c1,c2
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg,fname
      CHARACTER(8), DIMENSION(mode_st:mode_en) :: fcomps

      TYPE(r4eldata_type), DIMENSION(nbl) :: eldivb
      LOGICAL :: renumel
c-----------------------------------------------------------------------
c     create a name for the file, and load a character
c     array with the Fourier component indices.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(a)') '>>>> Enter the root name for your divb'//
     $                    ' data file.'
      WRITE(nim_wr,'(a)',ADVANCE='NO') '>>>>? '
      READ(nim_rd,*) fname
      IF (elout_format=='vtk') THEN
        ipt=0
        renumel=.false.
        fname=TRIM(fname)//'.vtk'
      ELSE  !  'tecplot'
        ipt=1
        renumel=.true.
        fname=TRIM(fname)//'.plt'
      ENDIF
c-----------------------------------------------------------------------
c     create component labels.
c-----------------------------------------------------------------------
      DO im=mode_st,mode_en
        IF (geom=='tor') THEN
          n1=NINT(keff(im))
        ELSE
          n1=NINT(per_length*keff(im)/twopi)
        ENDIF
        IF (n1>=10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE IF (n1>=0) THEN
          WRITE(fcomps(im),'(i1)') n1
        ELSE IF (n1>-10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE 
          WRITE(fcomps(im),'(i3)') n1
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     prepare the generic element-data structure.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL eldata_alloc(eldivb(ibl),"quad",poly_degree,
     $                    2_i4*(mode_en-mode_st+1_i4),
     $                    mxb*myb,ipt,renumel)
        IF (geom=='tor') THEN
          eldivb(ibl)%coordlabels='RZ'
        ELSE
          eldivb(ibl)%coordlabels='xy'
        ENDIF
        ifc=1
        DO im=mode_st,mode_en
          eldivb(ibl)%datalabels(2*ifc-1)='Re(divb) n='//
     $      TRIM(fcomps(im))
          eldivb(ibl)%datalabels(2*ifc  )='Im(divb) n='//
     $      TRIM(fcomps(im))
          ifc=ifc+1
        ENDDO
c-----------------------------------------------------------------------
c       elements may be in any order, because points are only defined
c       within elements.  however, data within elements matches the
c       connectivity defined in eldata_type_mod.
c-----------------------------------------------------------------------
        iel=1
        DO iy=1,myb
          DO ix=1,mxb
            idat=1
            DO iby=0,poly_degree
              dy=eps+(1-2*eps)*REAL(iby)/poly_degree
              y=iy-1+dy
              DO ibx=0,poly_degree
                dx=eps+(1-2*eps)*REAL(ibx)/poly_degree
                x=ix-1+dx
c-----------------------------------------------------------------------
c               evaluate metric information.
c-----------------------------------------------------------------------
                CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
                CALL lagr_quad_eval(rb(ibl)%be,x,y,1_i4)
                drdx=rb(ibl)%rz%fx(1)
                drdy=rb(ibl)%rz%fy(1)
                dzdx=rb(ibl)%rz%fx(2)
                dzdy=rb(ibl)%rz%fy(2)
                jac=drdx*dzdy-drdy*dzdx
                dxdr= dzdy/jac
                dxdz=-drdy/jac
                dydr=-dzdx/jac
                dydz= drdx/jac
c-----------------------------------------------------------------------
c               save the physical coordinates.
c-----------------------------------------------------------------------
                eldivb(ibl)%coords(:,idat,iel)=rb(ibl)%rz%f
c-----------------------------------------------------------------------
c               compute div(b) for each Fourier component.
c-----------------------------------------------------------------------
                ifc=1
                DO im=mode_st,mode_en
                  dbrdx=rb(ibl)%be%fx(1,im)
                  dbrdy=rb(ibl)%be%fy(1,im)
                  dbzdx=rb(ibl)%be%fx(2,im)
                  dbzdy=rb(ibl)%be%fy(2,im)
                  bp=rb(ibl)%be%f(3,im)
                  br=rb(ibl)%be%f(1,im)
                  IF (geom=='tor') THEN
                    divb=dbrdx*dxdr+dbrdy*dydr+dbzdx*dxdz+dbzdy*dydz
     $                  +(0,1)*keff(im)*bp/rb(ibl)%rz%f(1)
     $                  +br/rb(ibl)%rz%f(1)
                  ELSE
                    divb=dbrdx*dxdr+dbrdy*dydr+dbzdx*dxdz+dbzdy*dydz
     $                  +(0,1)*keff(im)*bp
                  ENDIF
                  eldivb(ibl)%data(2*ifc-1,idat,iel)=REAL(divb,4)
                  eldivb(ibl)%data(2*ifc,idat,iel)=REAL(AIMAG(divb),4)
                  ifc=ifc+1
                ENDDO
                idat=idat+1
              ENDDO
            ENDDO
            iel=iel+1
          ENDDO
        ENDDO       
      ENDDO
c-----------------------------------------------------------------------
c     write the data to the specified file.
c-----------------------------------------------------------------------
      SELECT CASE (elout_format)
      CASE('vtk')
        CALL vtk_eldata_write(eldivb,TRIM(fname),60_i4)
      CASE('tecplot')
        CALL tec_eldata_write(eldivb,TRIM(fname),60_i4)
      END SELECT

      DO ib=1,nrbl
        CALL eldata_dealloc(eldivb(ib))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE divb_element
c-----------------------------------------------------------------------
c     subroutine 8.  j_element
c     find J within elements and create a file that can
c     be used to create contour plots for fields that are discontinuous
c     at element boundaries.
c-----------------------------------------------------------------------
      SUBROUTINE j_element(datflag,mode_st,mode_en)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      USE eldata_type_mod
      USE tecplot_mod
      USE vtk_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en
      TYPE(vector_type), DIMENSION(nbl) :: divbn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb, 
     $               ix,iy,ib,iby,ibx,n1,n2,n3,n4,ipt,iel,idat,ifc
      REAL(r8) :: err,dx,dy,x,y,eps=1.e-6
      REAL(r8) :: drdx,drdy,dzdx,dzdy
      REAL(r8) :: jac,dxdr,dxdz,dydr,dydz
      COMPLEX(r8) :: dbrdx,dbrdy,dbzdx,dbzdy,dbpdx,dbpdy,
     $               br,bz,bp,jr,jz,jp
      CHARACTER(1) :: c1,c2,c3
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg,fname
      CHARACTER(8), DIMENSION(mode_st:mode_en) :: fcomps

      TYPE(r4eldata_type), DIMENSION(nbl) :: elj
      LOGICAL :: renumel
c-----------------------------------------------------------------------
c     create a name for the file, and load a character
c     array with the Fourier component indices.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(a)') '>>>> Enter the root name for your j'//
     $                    ' data file.'
      WRITE(nim_wr,'(a)',ADVANCE='NO') '>>>>? '
      READ(nim_rd,*) fname
      IF (elout_format=='vtk') THEN
        ipt=0
        renumel=.false.
        fname=TRIM(fname)//'.vtk'
      ELSE  !  'tecplot'
        ipt=1
        renumel=.true.
        fname=TRIM(fname)//'.plt'
      ENDIF
c-----------------------------------------------------------------------
c     create component labels.
c-----------------------------------------------------------------------
      DO im=mode_st,mode_en
        IF (geom=='tor') THEN
          n1=NINT(keff(im))
        ELSE
          n1=NINT(per_length*keff(im)/twopi)
        ENDIF
        IF (n1>=10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE IF (n1>=0) THEN
          WRITE(fcomps(im),'(i1)') n1
        ELSE IF (n1>-10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE 
          WRITE(fcomps(im),'(i3)') n1
        ENDIF
      ENDDO
      IF (geom=='tor') THEN
        c1='R'
        c2='Z'
        c3='P'
      ELSE
        c1='x'
        c2='y'
        c3='z'
      ENDIF
c-----------------------------------------------------------------------
c     prepare the generic element-data structure.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL eldata_alloc(elj(ibl),"quad",poly_degree,
     $                    6_i4*(mode_en-mode_st+1_i4),
     $                    mxb*myb,ipt,renumel)
        elj(ibl)%coordlabels=c1//c2
        ifc=1
        DO im=mode_st,mode_en
          elj(ibl)%datalabels(6*ifc-5)=
     $      'Re(J'//c1//') n='//TRIM(fcomps(im))
          elj(ibl)%datalabels(6*ifc-4)=
     $      'Re(J'//c2//') n='//TRIM(fcomps(im))
          elj(ibl)%datalabels(6*ifc-3)=
     $      'Re(J'//c3//') n='//TRIM(fcomps(im))
          elj(ibl)%datalabels(6*ifc-2)=
     $      'Im(J'//c1//') n='//TRIM(fcomps(im))
          elj(ibl)%datalabels(6*ifc-1)=
     $      'Im(J'//c2//') n='//TRIM(fcomps(im))
          elj(ibl)%datalabels(6*ifc  )=
     $      'Im(J'//c3//') n='//TRIM(fcomps(im))
          ifc=ifc+1
        ENDDO
c-----------------------------------------------------------------------
c       elements may be in any order, because points are only defined
c       within elements.  however, data within elements matches the
c       connectivity defined in eldata_type_mod.
c-----------------------------------------------------------------------
        iel=1
        DO iy=1,myb
          DO ix=1,mxb
            idat=1
            DO iby=0,poly_degree
              dy=eps+(1-2*eps)*REAL(iby)/poly_degree
              y=iy-1+dy
              DO ibx=0,poly_degree
                dx=eps+(1-2*eps)*REAL(ibx)/poly_degree
                x=ix-1+dx
c-----------------------------------------------------------------------
c               evaluate metric information.
c-----------------------------------------------------------------------
                CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
                CALL lagr_quad_eval(rb(ibl)%be,x,y,1_i4)
                drdx=rb(ibl)%rz%fx(1)
                drdy=rb(ibl)%rz%fy(1)
                dzdx=rb(ibl)%rz%fx(2)
                dzdy=rb(ibl)%rz%fy(2)
                jac=drdx*dzdy-drdy*dzdx
                dxdr= dzdy/jac
                dxdz=-drdy/jac
                dydr=-dzdx/jac
                dydz= drdx/jac
c-----------------------------------------------------------------------
c               save the physical coordinates.
c-----------------------------------------------------------------------
                elj(ibl)%coords(:,idat,iel)=rb(ibl)%rz%f
c-----------------------------------------------------------------------
c               compute curl(b) for each Fourier component.
c-----------------------------------------------------------------------
                ifc=1
                DO im=mode_st,mode_en
                  dbrdx=rb(ibl)%be%fx(1,im)
                  dbrdy=rb(ibl)%be%fy(1,im)
                  dbzdx=rb(ibl)%be%fx(2,im)
                  dbzdy=rb(ibl)%be%fy(2,im)
                  dbpdx=rb(ibl)%be%fx(3,im)
                  dbpdy=rb(ibl)%be%fy(3,im)
                  br=rb(ibl)%be%f(1,im)
                  bz=rb(ibl)%be%f(2,im)
                  bp=rb(ibl)%be%f(3,im)
                  IF (geom=='tor') THEN
                    jr=(dbpdx*dxdz+dbpdy*dydz
     $                 -(0,1)*keff(im)*bz/rb(ibl)%rz%f(1))/mu0
                    jz=(((0,1)*keff(im)*br-bp)/rb(ibl)%rz%f(1)
     $                 -dbpdx*dxdr-dbpdy*dydr)/mu0
                  ELSE
                    jr=(dbpdx*dxdz+dbpdy*dydz-(0,1)*keff(im)*bz)/mu0
                    jz=((0,1)*keff(im)*br-dbpdx*dxdr-dbpdy*dydr)/mu0
                  ENDIF
                  jp=(dbzdx*dxdr+dbzdy*dydr-dbrdx*dxdz-dbrdy*dydz)/mu0
                  elj(ibl)%data(6*ifc-5:6*ifc-3,idat,iel)=
     $              REAL((/jr,jz,jp/),4)
                  elj(ibl)%data(6*ifc-2:6*ifc  ,idat,iel)=
     $              REAL(AIMAG((/jr,jz,jp/)),4)
                  ifc=ifc+1
                ENDDO
                idat=idat+1
              ENDDO
            ENDDO
            iel=iel+1
          ENDDO
        ENDDO       
      ENDDO
c-----------------------------------------------------------------------
c     write the data to the specified file.
c-----------------------------------------------------------------------
      SELECT CASE (elout_format)
      CASE('vtk')
        CALL vtk_eldata_write(elj,TRIM(fname),60_i4)
      CASE('tecplot')
        CALL tec_eldata_write(elj,TRIM(fname),60_i4)
      END SELECT

      DO ib=1,nrbl
        CALL eldata_dealloc(elj(ib))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE j_element
c-----------------------------------------------------------------------
c     subroutine 9.  dyn_element
c     find the local MHD and Hall dynamo fields parallel to the
c     average B within elements and create a file that can
c     be used to create contour plots for fields that are discontinuous
c     at element boundaries.
c-----------------------------------------------------------------------
      SUBROUTINE dyn_element(datflag,mode_st,mode_en)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      USE physdat
      USE eldata_type_mod
      USE tecplot_mod
      USE vtk_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en
      TYPE(vector_type), DIMENSION(nbl) :: divbn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb, 
     $               ix,iy,ib,iby,ibx,n1,n2,n3,n4,ipt,iel,idat,ifc
      REAL(r8) :: err,dx,dy,x,y,eps=1.e-6
      REAL(r8) :: drdx,drdy,dzdx,dzdy
      REAL(r8) :: jac,dxdr,dxdz,dydr,dydz
      REAL(r8) :: bh1,bh2,bh3,bxv,jxb,nda,bmag
      COMPLEX(r8) :: dbrdx,dbrdy,dbzdx,dbzdy,dbpdx,dbpdy,
     $               br,bz,bp,jr,jz,jp,vr,vz,vp
      CHARACTER(1) :: c1,c2
      CHARACTER(64) :: msg,fname
      CHARACTER(8), DIMENSION(mode_st:mode_en) :: fcomps

      TYPE(r4eldata_type), DIMENSION(nbl) :: eldyn
      LOGICAL :: renumel
c-----------------------------------------------------------------------
c     create a name for the file, and load a character
c     array with the Fourier component indices.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(a)') '>>>> Enter the root name for your dynamo'//
     $                    ' data file.'
      WRITE(nim_wr,'(a)',ADVANCE='NO') '>>>>? '
      READ(nim_rd,*) fname
      IF (elout_format=='vtk') THEN
        ipt=0
        renumel=.false.
        fname=TRIM(fname)//'.vtk'
      ELSE  !  'tecplot'
        ipt=1
        renumel=.true.
        fname=TRIM(fname)//'.plt'
      ENDIF
c-----------------------------------------------------------------------
c     create component labels.
c-----------------------------------------------------------------------
      DO im=mode_st,mode_en
        IF (geom=='tor') THEN
          n1=NINT(keff(im))
        ELSE
          n1=NINT(per_length*keff(im)/twopi)
        ENDIF
        IF (n1>=10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE IF (n1>=0) THEN
          WRITE(fcomps(im),'(i1)') n1
        ELSE IF (n1>-10) THEN
          WRITE(fcomps(im),'(i2)') n1
        ELSE
          WRITE(fcomps(im),'(i3)') n1
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     prepare the generic element-data structure.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL eldata_alloc(eldyn(ibl),"quad",poly_degree,
     $                    3_i4*(mode_en-mode_st+1_i4),
     $                    mxb*myb,ipt,renumel)
        IF (geom=='tor') THEN
          eldyn(ibl)%coordlabels='RZ'
        ELSE
          eldyn(ibl)%coordlabels='xy'
        ENDIF
        ifc=1
        DO im=mode_st,mode_en
          eldyn(ibl)%datalabels(3*ifc-2)='eta*J_z n='//TRIM(fcomps(im))
          eldyn(ibl)%datalabels(3*ifc-1)='-<vXb>_z n='//TRIM(fcomps(im))
          eldyn(ibl)%datalabels(3*ifc)=
     $      '<jxb>_z/<n>e n='//TRIM(fcomps(im))
          ifc=ifc+1
        ENDDO
c-----------------------------------------------------------------------
c       elements may be in any order, because points are only defined
c       within elements.  however, data within elements matches the
c       connectivity defined in eldata_type_mod.
c-----------------------------------------------------------------------
        iel=1
        DO iy=1,myb
          DO ix=1,mxb
            idat=1
            DO iby=0,poly_degree
              dy=eps+(1-2*eps)*REAL(iby)/poly_degree
              y=iy-1+dy
              DO ibx=0,poly_degree
                dx=eps+(1-2*eps)*REAL(ibx)/poly_degree
                x=ix-1+dx
c-----------------------------------------------------------------------
c               evaluate metric information.
c-----------------------------------------------------------------------
                CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
                CALL lagr_quad_eval(rb(ibl)%be,x,y,1_i4)
                CALL lagr_quad_eval(rb(ibl)%ve,x,y,0_i4)
                CALL lagr_quad_eval(rb(ibl)%nd,x,y,0_i4)
                CALL lagr_quad_eval(rb(ibl)%be_eq,x,y,1_i4)
                CALL lagr_quad_eval(rb(ibl)%nd_eq,x,y,0_i4)
                drdx=rb(ibl)%rz%fx(1)
                drdy=rb(ibl)%rz%fy(1)
                dzdx=rb(ibl)%rz%fx(2)
                dzdy=rb(ibl)%rz%fy(2)
                jac=drdx*dzdy-drdy*dzdx
                dxdr= dzdy/jac
                dxdz=-drdy/jac
                dydr=-dzdx/jac
                dydz= drdx/jac
c-----------------------------------------------------------------------
c               save the physical coordinates.
c-----------------------------------------------------------------------
                eldyn(ibl)%coords(:,idat,iel)=rb(ibl)%rz%f
c-----------------------------------------------------------------------
c               first find the direction of the average magnetic field.
c-----------------------------------------------------------------------
                bh1=rb(ibl)%be_eq%f(1)
                bh2=rb(ibl)%be_eq%f(2)
                bh3=rb(ibl)%be_eq%f(3)
                nda=rb(ibl)%nd_eq%f(1)
                IF (nonlinear) THEN
                  bh1=bh1+rb(ibl)%be%f(1,1)
                  bh2=bh2+rb(ibl)%be%f(2,1)
                  bh3=bh3+rb(ibl)%be%f(3,1)
                  nda=nda+rb(ibl)%nd%f(1,1)
                ENDIF
                bmag=SQRT(bh1**2+bh2**2+bh3**2)
c-----------------------------------------------------------------------
c               two-dimensional computations just consider perp-comp.
c-----------------------------------------------------------------------
                IF (lphi<=1) THEN
                  bh1=0._r8
                  bh2=0._r8
                  bh3=1._r8
                ELSE
                  bh1=bh1/bmag
                  bh2=bh2/bmag
                  bh3=bh3/bmag
                ENDIF
c-----------------------------------------------------------------------
c               start a loop over the selected range of Fourier
c               components.
c-----------------------------------------------------------------------
                ifc=1
                DO im=mode_st,mode_en
                  dbrdx=rb(ibl)%be%fx(1,im)
                  dbrdy=rb(ibl)%be%fy(1,im)
                  dbzdx=rb(ibl)%be%fx(2,im)
                  dbzdy=rb(ibl)%be%fy(2,im)
                  dbpdx=rb(ibl)%be%fx(3,im)
                  dbpdy=rb(ibl)%be%fy(3,im)
                  br=rb(ibl)%be%f(1,im)
                  bz=rb(ibl)%be%f(2,im)
                  bp=rb(ibl)%be%f(3,im)
                  IF (keff(im)==0) THEN
                    br=br+rb(ibl)%be_eq%f(1)
                    bz=bz+rb(ibl)%be_eq%f(2)
                    bp=bp+rb(ibl)%be_eq%f(3)
                    dbrdx=dbrdx+rb(ibl)%be_eq%fx(1)
                    dbrdy=dbrdy+rb(ibl)%be_eq%fy(1)
                    dbzdx=dbzdx+rb(ibl)%be_eq%fx(2)
                    dbzdy=dbzdy+rb(ibl)%be_eq%fy(2)
                    dbpdx=dbpdx+rb(ibl)%be_eq%fx(3)
                    dbpdy=dbpdy+rb(ibl)%be_eq%fy(3)
                  ENDIF
                  vr=rb(ibl)%ve%f(1,im)
                  vz=rb(ibl)%ve%f(2,im)
                  vp=rb(ibl)%ve%f(3,im)
c-----------------------------------------------------------------------
c                 compute current density from the local curl of b.
c-----------------------------------------------------------------------
                  IF (geom=='tor') THEN
                    jr=(dbpdx*dxdz+dbpdy*dydz
     $                 -(0,1)*keff(im)*bz/rb(ibl)%rz%f(1))/mu0
                    jz=(((0,1)*keff(im)*br-bp)/rb(ibl)%rz%f(1)
     $                 -dbpdx*dxdr-dbpdy*dydr)/mu0
                  ELSE
                    jr=(dbpdx*dxdz+dbpdy*dydz-(0,1)*keff(im)*bz)/mu0
                    jz=((0,1)*keff(im)*br-dbpdx*dxdr-dbpdy*dydr)/mu0
                  ENDIF
                  jp=(dbzdx*dxdr+dbzdy*dydr-dbrdx*dxdz-dbrdy*dydz)/mu0
c-----------------------------------------------------------------------
c                 find the parallel contribution to the MHD dynamo.
c-----------------------------------------------------------------------
                  bxv=bh1*(bz*CONJG(vp)-bp*CONJG(vz))+
     $                bh2*(bp*CONJG(vr)-br*CONJG(vp))+
     $                bh3*(br*CONJG(vz)-bz*CONJG(vr))
                  IF (keff(im)/=0) THEN
                    bxv=bxv+bh1*(CONJG(bz)*vp-CONJG(bp)*vz)+
     $                      bh2*(CONJG(bp)*vr-CONJG(br)*vp)+
     $                      bh3*(CONJG(br)*vz-CONJG(bz)*vr)
                  ENDIF
c-----------------------------------------------------------------------
c                 find the parallel contribution to the Hall dynamo.
c-----------------------------------------------------------------------
                  jxb=bh1*(jz*CONJG(bp)-jp*CONJG(bz))+
     $                bh2*(jp*CONJG(br)-jr*CONJG(bp))+
     $                bh3*(jr*CONJG(bz)-jz*CONJG(br))
                  IF (keff(im)/=0) THEN
                    jxb=jxb+bh1*(CONJG(jz)*bp-CONJG(jp)*bz)+
     $                      bh2*(CONJG(jp)*br-CONJG(jr)*bp)+
     $                      bh3*(CONJG(jr)*bz-CONJG(jz)*br)
                  ENDIF
                  jxb=jxb/(nda*elementary_q)
c-----------------------------------------------------------------------
c                 save the electric-field values.
c-----------------------------------------------------------------------
                  eldyn(ibl)%data(3*ifc-2:3*ifc,idat,iel)=
     $              (/REAL(elecd*mu0*jp,4),REAL(bxv,4),REAL(jxb,4)/)
                  ifc=ifc+1
                ENDDO
                idat=idat+1
              ENDDO
            ENDDO
            iel=iel+1
          ENDDO
        ENDDO       
      ENDDO
c-----------------------------------------------------------------------
c     write the data to the specified file.
c-----------------------------------------------------------------------
      SELECT CASE (elout_format)
      CASE('vtk')
        CALL vtk_eldata_write(eldyn,TRIM(fname),60_i4)
      CASE('tecplot')
        CALL tec_eldata_write(eldyn,TRIM(fname),60_i4)
      END SELECT

      DO ib=1,nrbl
        CALL eldata_dealloc(eldyn(ib))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dyn_element
c-----------------------------------------------------------------------
c     subprogram 10. mach.
c     find ion acoustic speed, Mach number and parallel Mach number at
c     the nodes for xdraw plots.
c-----------------------------------------------------------------------
      SUBROUTINE mach(first_data,datflag,mode_st,mode_en,phii,ndcon)
      USE local
      USE fields
      USE input
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimplot_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en,ndcon
      REAL(r8), INTENT(IN) :: phii
      TYPE(vector_type), DIMENSION(nbl) :: machn0

      INTEGER(i4) :: ibl,iq,nq,ir,ii,iqst,iqen,jphi,im,its,mxb,myb
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     find integral[alpha*v/ion_aspd], integral[alpha*ion_aspd]
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes)
        CALL rblock_get_rhs(rb(ibl),crhs(ibl),mach_int,
     $                        3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes)
        CALL rblock_get_rhs(rb(ibl),crhs(ibl),mach_int,
     $                        3_i4,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     network block seams.
c-----------------------------------------------------------------------
      CALL edge_network(3_i4,nmodes,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c     invert the mass matrix for each component and each computed
c     quantity.
c-----------------------------------------------------------------------
      DO im=mode_st,mode_en
        DO iq=1,3
          DO ibl=1,nbl
            csln(ibl)=0
            cvectr(ibl)%arr(1,:,:)=crhs(ibl)%arr(iq,:,:,im)
            IF (ASSOCIATED(csln(ibl)%arrh).AND.ibl<=nrbl) THEN
              cvectr(ibl)%arrh(1,:,:,:)=crhs(ibl)%arrh(iq,:,:,:,im)
              cvectr(ibl)%arrv(1,:,:,:)=crhs(ibl)%arrv(iq,:,:,:,im)
              cvectr(ibl)%arri(1,:,:,:)=crhs(ibl)%arri(iq,:,:,:,im)
            ENDIF
          ENDDO
          CALL iter_cg_2d_solve(cmass_mat,cmass_fac,csln,cvectr,1_i4,
     $                          tol,maxit,solver,err,its,seed)
          IF (err>tol) THEN
            WRITE(msg,'(a,i4,a,es10.3,3a)') 'mach: no convergence: ',
     $        its,' its ',err,' err ',seed,' seed'
            CALL nim_stop(msg)
          ENDIF
          DO ibl=1,nbl
            work1(ibl)%arr(iq,:,:,im)=csln(ibl)%arr(1,:,:)
            IF (ASSOCIATED(csln(ibl)%arrh).AND.ibl<=nrbl) THEN
              work1(ibl)%arrh(iq,:,:,:,im)=csln(ibl)%arrh(1,:,:,:)
              work1(ibl)%arrv(iq,:,:,:,im)=csln(ibl)%arrv(1,:,:,:)
              work1(ibl)%arri(iq,:,:,:,im)=csln(ibl)%arri(1,:,:,:)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     initialize contour plot if needed, otherwise just open the file.
c     option c plots machn at a specified angle in configuration space.
c     there are 6 real (3 complex) fields per Fourier component.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        nq=6*(mode_en-mode_st+1)
        iqst=mode_st
        iqen=mode_en
      ELSE
        nq=6
      ENDIF
      IF (first_data) THEN
        CALL contour_init(nq,ndcon,"mach.bin")
      ELSE
        CALL open_bin(con_unit,"mach.bin","OLD","APPEND",32_i4)
      ENDIF
c-----------------------------------------------------------------------
c     write the appropriate data.  save and restore the n=0 part if
c     plotting at angles.
c-----------------------------------------------------------------------
      IF (datflag=='a'.OR.datflag=='o') THEN
        DO ibl=1,nrbl
          CALL contour_laq_write(rb(ibl)%work1,iqst,iqen,ndcon)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL contour_tl_write(tb(ibl)%work1,iqst,iqen)
        ENDDO
      ELSE
        jphi=0
        DO ibl=1,nbl
          mxb=SIZE(work1(ibl)%arr,2)-1
          myb=SIZE(work1(ibl)%arr,3)-1
          CALL vector_type_alloc(machn0(ibl),poly_degree,mxb,myb,3_i4)
          CALL vector_assign_cvec(machn0(ibl),work1(ibl),'real',1_i4)
        ENDDO
        phi_loop: DO
          IF (datflag=='c'.AND.jphi/=1) THEN
            jphi=jphi+1
            IF (jphi==1) CYCLE phi_loop
            EXIT phi_loop
          ENDIF
          CALL eval_at_angle(work1,MODULO(MIN(jphi*phii,1._r8),1._r8))
          DO ibl=1,nrbl
            CALL contour_laq_write(rb(ibl)%work1,1_i4,1_i4,ndcon)
            CALL cvector_assign_vec(work1(ibl),machn0(ibl),'real',1_i4)
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL contour_tl_write(tb(ibl)%work1,1_i4,1_i4)
            CALL cvector_assign_vec(work1(ibl),machn0(ibl),'real',1_i4)
          ENDDO
          IF (jphi*phii>=1) EXIT phi_loop
          jphi=jphi+1
        ENDDO phi_loop
        DO ibl=1,nbl
          CALL vector_type_dealloc(machn0(ibl))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close the plot file.
c-----------------------------------------------------------------------
      CALL close_bin(con_unit,"mach.bin")
c-----------------------------------------------------------------------
c     deallocate rhs storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mach
