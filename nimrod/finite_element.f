c-----------------------------------------------------------------------
c     file finite_element.f
c     module containing routines that call the appropriate finite
c     element computation routines for each block.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. finite_element_mod.
c     1. real_matrix_create.
c     2. comp_matrix_create.
c     3. get_rhs.
c     4. real_fe_postsolve.
c     5. comp_fe_postsolve.
c-----------------------------------------------------------------------
c     0. module declarations for finite_element_mod.
c-----------------------------------------------------------------------
      MODULE finite_element_mod
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE surface
      USE matrix_mod
      USE time
      IMPLICIT NONE

      INTERFACE matrix_create
        MODULE PROCEDURE real_matrix_create,comp_matrix_create
      END INTERFACE

      INTERFACE fe_postsolve
        MODULE PROCEDURE real_fe_postsolve,comp_fe_postsolve
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. real_matrix_create.
c     calls the block-dependent finite element computation for a
c     real matrix specified in the parameter list.
c-----------------------------------------------------------------------
      SUBROUTINE real_matrix_create(matrix_structure,factor_structure,
     $                              integrand,bc_routine,bc_flag,ident,
     $                              int_elim,solver,mass_type)
      USE iter_cg
      USE pardata
      USE regularity
      USE mpi_nim

      TYPE(global_matrix_type), DIMENSION(:), POINTER ::
     $                          matrix_structure
      TYPE(matrix_factor_type), DIMENSION(:), POINTER ::
     $                          factor_structure
      LOGICAL, INTENT(IN) :: int_elim
      CHARACTER(*), INTENT(IN) :: bc_flag,ident,solver
      CHARACTER(*), INTENT(IN), OPTIONAL :: mass_type

      INTEGER(i4) :: ibl,nn,iv,nq,nqdis,ix,myb,iq,ierror
      REAL(r8), DIMENSION(2) :: df
      REAL(r8) :: dscale,tmp
      LOGICAL :: do_mass,do_fac
c-----------------------------------------------------------------------
c     interface block for the integrand routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE integrand(int,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the boundary condition routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE bc_routine(mat,flag,scale,symm)
        USE local
        USE matrix_type_mod
        TYPE(global_matrix_type), INTENT(INOUT) :: mat
        CHARACTER(*), INTENT(IN) :: flag
        REAL(r8), INTENT(IN) :: scale
        CHARACTER(*), INTENT(IN), OPTIONAL :: symm
        END SUBROUTINE bc_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     write messages to the output file.
c-----------------------------------------------------------------------
      IF (node==0) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF
        WRITE(out_unit,'(/a,i6)')
     $   ' Matrix_create is computing a new '//TRIM(ident)//
     $   ' matrix at cycle ',istep
      ENDIF
c-----------------------------------------------------------------------
c     check if certain operations are needed.
c-----------------------------------------------------------------------
      do_mass=.true.
      IF (PRESENT(mass_type)) THEN
        IF (mass_type=='none') do_mass=.false.
      ENDIF
      do_fac=ASSOCIATED(factor_structure)
c-----------------------------------------------------------------------
c     begin the loop over modes--jmode is in global and is used in the
c     integrand routines to determine wavenumber.
c-----------------------------------------------------------------------
      mode_loop: DO jmode=1,nmodes
c-----------------------------------------------------------------------
c       determine the matrix vector component index range for continuous
c       (nq) and any discontinuous (nqdis) fields.
c-----------------------------------------------------------------------
        nq=matrix_structure(jmode)%nqty
        nqdis=matrix_structure(jmode)%nqdis
c-----------------------------------------------------------------------
c       evaluate the operator from the finite element integrations.
c       also set the global-module variables ncontb and ndiscb to the
c       numbers of continuous and discontinuous bases in each block.
c-----------------------------------------------------------------------
        DO ibl=1,nrbl
          ncontb=matrix_structure(jmode)%rbl_mat(ibl)%nbasis_cont
          ndiscb=matrix_structure(jmode)%rbl_mat(ibl)%nbasis_disc
          CALL rblock_make_matrix
     $      (rb(ibl),matrix_structure(jmode)%rbl_mat(ibl),integrand,
     $       MAX(nq,nqdis))
        ENDDO
        DO ibl=nrbl+1,nbl
          ncontb=3
          ndiscb=0
          CALL tblock_make_matrix
     $      (tb(ibl),matrix_structure(jmode)%tbl_mat(ibl)%lmat,
     $       integrand,nq)
        ENDDO
c-----------------------------------------------------------------------
c       add mass matrix, then determine a scaling factor for diagonal
c       matrix elements that is based on grid-vertex entries.
c
c       the scaling factor, regularity, and boundary operations are
c       not called if there are no continuous fields.
c-----------------------------------------------------------------------
        IF (do_mass) CALL add_mass(matrix_structure(jmode),nq)
        IF (nq>0) THEN
          dscale=0._r8
          DO ibl=1,nrbl
            DO iq=1,nq
              dscale=MAX(dscale,MAXVAL(ABS(matrix_structure(jmode)%
     $                                     rbl_mat(ibl)%mat(1,1)%
     $                                     arr(iq,0,0,iq,:,:))))
            ENDDO
          ENDDO
          DO ibl=nrbl+1,nbl
            DO iv=0,SIZE(matrix_structure(jmode)%tbl_mat(ibl)%lmat)-1
              DO iq=1,nq
                dscale=MAX(dscale,ABS(matrix_structure(jmode)%
     $                                tbl_mat(ibl)%lmat(iv)%
     $                                element(iq,iq,0)))
              ENDDO
            ENDDO
          ENDDO
          IF (nprocs_layer>1) THEN
            CALL mpi_allreduce(dscale,tmp,1,mpi_nim_real,mpi_max,
     $           comm_layer,ierror)
            dscale=tmp
          ENDIF
          matrix_structure(jmode)%diag_scale=dscale
c-----------------------------------------------------------------------
c         call the regularity condition routine for matrices then
c         the boundary condition routine.
c-----------------------------------------------------------------------
          CALL regular_op(matrix_structure(jmode),dscale)
          CALL bc_routine(matrix_structure(jmode),bc_flag,dscale)
        ENDIF
c-----------------------------------------------------------------------
c       eliminate connections to cell-centered data (basis functions
c       of degree greater than linear) and save the inverse of the
c       interior block.
c-----------------------------------------------------------------------
        IF (int_elim) THEN
          CALL timer(timestart)
          CALL matelim_real_inv_int(matrix_structure(jmode),nq)
          CALL timer(timeend)
          time_stcon=time_stcon+timeend-timestart
        ENDIF
c-----------------------------------------------------------------------
c       sum off-diagonal elements along degenerate boundaries into the
c       adjacent diagonal elements.  this does not change the matrix
c       computed by a matrix-vector multiply, but it helps when forming
c       the preconditioners.  note: direct-solve preconditioning does
c       not use all connections, and the collection tends to create
c       singular matrices when nybl=1.  do not apply the collection for
c       these cases.
c-----------------------------------------------------------------------
        deg_bl: DO ibl=1,nrbl
          IF (rb(ibl)%degenerate.AND.nq>0) THEN
            IF (solver(1:8)=='bl_drect') THEN
              myb=rb(ibl)%my
              DO ix=1,rb(ibl)%mx
                df=ABS(rb(ibl)%rz%fs(:,ix,0)-rb(ibl)%rz%fs(:,ix,myb))
                IF (df(1)<1.e-12.AND.df(2)<1.e-12) CYCLE deg_bl
              ENDDO
            ENDIF
            CALL matrix_degen_collect_real
     $        (matrix_structure(jmode)%rbl_mat(ibl),nq,
     $         matrix_structure(jmode)%symmetric)
          ENDIF
        ENDDO deg_bl
c-----------------------------------------------------------------------
c       call matrix factorization routine.
c-----------------------------------------------------------------------
        IF (do_fac) CALL iter_factor(matrix_structure(jmode),
     $                   factor_structure(jmode),nq,solver,off_diag_fac)
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE real_matrix_create
c-----------------------------------------------------------------------
c     subprogram 2. comp_matrix_create.
c     calls the block-dependent finite element computation for a
c     complex matrix specified in the parameter list.
c-----------------------------------------------------------------------
      SUBROUTINE comp_matrix_create(matrix_structure,factor_structure,
     $                              integrand,bc_routine,bc_flag,ident,
     $                              int_elim,solver,mass_type)
      USE iter_cg
      USE pardata
      USE regularity
      USE mpi_nim

      TYPE(complex_matrix_type), DIMENSION(:), POINTER ::
     $                           matrix_structure
      TYPE(complex_factor_type), DIMENSION(:), POINTER ::
     $                           factor_structure
      LOGICAL, INTENT(IN) :: int_elim
      CHARACTER(*), INTENT(IN) :: bc_flag,ident,solver
      CHARACTER(*), INTENT(IN), OPTIONAL :: mass_type

      INTEGER(i4) :: ibl,nn,iv,nq,nqdis,ix,myb,iq,ierror
      REAL(r8), DIMENSION(2) :: df
      REAL(r8) :: dscale,tmp
      LOGICAL :: do_mass,do_fac
c-----------------------------------------------------------------------
c     interface block for the integrand routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE integrand(int,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the boundary condition routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE bc_routine(mat,flag,scale)
        USE local
        USE matrix_type_mod
        TYPE(complex_matrix_type), INTENT(INOUT) :: mat
        CHARACTER(*), INTENT(IN) :: flag
        REAL(r8), INTENT(IN) :: scale
        END SUBROUTINE bc_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     write messages to the output file.
c-----------------------------------------------------------------------
      IF (node==0) THEN
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF
        WRITE(out_unit,'(/a,i6)')
     $   ' Matrix_create is computing a new '//TRIM(ident)//
     $   ' matrix at cycle ',istep
      ENDIF
c-----------------------------------------------------------------------
c     check if certain operations are needed.
c-----------------------------------------------------------------------
      do_mass=.true.
      IF (PRESENT(mass_type)) THEN
        IF (mass_type=='none') do_mass=.false.
      ENDIF
      do_fac=ASSOCIATED(factor_structure)
c-----------------------------------------------------------------------
c     begin the loop over modes--jmode is in global and is used in the
c     integrand routines to determine wavenumber.
c-----------------------------------------------------------------------
      mode_loop: DO jmode=1,nmodes
c-----------------------------------------------------------------------
c       determine the matrix vector component index range for continuous
c       (nq) and any discontinuous (nqdis) fields.
c-----------------------------------------------------------------------
        nq=matrix_structure(jmode)%nqty
        nqdis=matrix_structure(jmode)%nqdis
c-----------------------------------------------------------------------
c       evaluate the operator from the finite element integrations.
c       also set the global-module variables ncontb and ndiscb to the
c       numbers of continuous and discontinuous bases in each block.
c-----------------------------------------------------------------------
        DO ibl=1,nrbl
          ncontb=matrix_structure(jmode)%rbl_mat(ibl)%nbasis_cont
          ndiscb=matrix_structure(jmode)%rbl_mat(ibl)%nbasis_disc
          CALL rblock_make_matrix
     $      (rb(ibl),matrix_structure(jmode)%rbl_mat(ibl),integrand,
     $       MAX(nq,nqdis))
        ENDDO
        DO ibl=nrbl+1,nbl
          ncontb=3
          ndiscb=0
          CALL tblock_make_matrix
     $      (tb(ibl),matrix_structure(jmode)%tbl_mat(ibl)%lmat,
     $       integrand,nq)
        ENDDO
c-----------------------------------------------------------------------
c       if this is a nonlinear computation without dealiasing, the
c       largest Fourier component must have a real equation.
c-----------------------------------------------------------------------
        IF (nonlinear.AND..NOT.dealiase.AND.nindex(jmode)==nphi/2) THEN
          DO ibl=1,nrbl
            CALL matrix_rbl_make_real(matrix_structure(jmode)%
     $                                rbl_mat(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL matrix_tbl_make_real(matrix_structure(jmode)%
     $                                tbl_mat(ibl))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       add mass matrix, then determine a scaling factor for diagonal
c       matrix elements that is based on grid-vertex entries.
c
c       the scaling factor, regularity, and boundary operations are
c       not called if there are no continuous fields.
c-----------------------------------------------------------------------
        IF (do_mass) CALL add_mass_comp(matrix_structure(jmode),nq)
        IF (nq>0) THEN
          dscale=0._r8
          DO ibl=1,nrbl
            DO iq=1,nq
              dscale=MAX(dscale,MAXVAL(ABS(matrix_structure(jmode)%
     $                                     rbl_mat(ibl)%mat(1,1)%
     $                                     arr(iq,0,0,iq,:,:))))
            ENDDO
          ENDDO
          DO ibl=nrbl+1,nbl
            DO iv=0,SIZE(matrix_structure(jmode)%tbl_mat(ibl)%lmat)-1
              DO iq=1,nq
                dscale=MAX(dscale,ABS(matrix_structure(jmode)%
     $                                tbl_mat(ibl)%lmat(iv)%
     $                                element(iq,iq,0)))
              ENDDO
            ENDDO
          ENDDO
          IF (nprocs_layer>1) THEN
            CALL mpi_allreduce(dscale,tmp,1,mpi_nim_real,mpi_max,
     $           comm_layer,ierror)
            dscale=tmp
          ENDIF
          matrix_structure(jmode)%diag_scale=dscale
c-----------------------------------------------------------------------
c         call the regularity condition routine for matrices then
c         the boundary condition routine.
c-----------------------------------------------------------------------
          CALL regular_comp_op(matrix_structure(jmode),dscale)
          CALL bc_routine(matrix_structure(jmode),bc_flag,dscale)
        ENDIF
c-----------------------------------------------------------------------
c       eliminate connections to cell-centered data (basis functions
c       of degree greater than linear) and save the inverse of the
c       interior block.
c-----------------------------------------------------------------------
        IF (int_elim) THEN
          CALL timer(timestart)
          CALL matelim_comp_inv_int(matrix_structure(jmode),nq)
          CALL timer(timeend)
          time_stcon=time_stcon+timeend-timestart
        ENDIF
c-----------------------------------------------------------------------
c       sum off-diagonal elements along degenerate boundaries into the
c       adjacent diagonal elements.  this does not change the matrix
c       computed by a matrix-vector multiply, but it helps when forming
c       the preconditioners.  note: direct-solve preconditioning does
c       not use all connections, and the collection tends to create
c       singular matrices when nybl=1.  do not apply the collection for
c       these cases.
c-----------------------------------------------------------------------
        deg_bl: DO ibl=1,nrbl
          IF (rb(ibl)%degenerate.AND.nq>0) THEN
            IF (solver(1:8)=='bl_drect') THEN
              myb=rb(ibl)%my
              DO ix=1,rb(ibl)%mx
                df=ABS(rb(ibl)%rz%fs(:,ix,0)-rb(ibl)%rz%fs(:,ix,myb))
                IF (df(1)<1.e-12.AND.df(2)<1.e-12) CYCLE deg_bl
              ENDDO
            ENDIF
            CALL matrix_degen_collect_comp
     $        (matrix_structure(jmode)%rbl_mat(ibl),nq,
     $         matrix_structure(jmode)%hermitian)
          ENDIF
        ENDDO deg_bl
c-----------------------------------------------------------------------
c       call matrix factorization routine.
c-----------------------------------------------------------------------
        IF (do_fac) CALL iter_factor(matrix_structure(jmode),
     $                   factor_structure(jmode),nq,solver,off_diag_fac)
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE comp_matrix_create
c-----------------------------------------------------------------------
c     subprogram 3. get_rhs.
c     calls the block-dependent finite element rhs computation for the
c     complex (only) integrand specified in the parameter list.
c
c     this routine is suitable for equations with fields that are
c     continuous across element borders, those with only discontinuous
c     fields, and those with a combination.  the passed data structure
c     rhsdum, and in particular which of its arrays are allocated, is
c     used to determine the necessary steps.  if there are no continuous
c     fields, only rhsdum%arri is allocated and used.  if there are both
c     continuous and discontinuous fields, rhsdum%arrtmp is used for the
c     discontinuous fields.
c
c     if there are no continuous fields, the presolve call fully solves
c     the element-local system.
c-----------------------------------------------------------------------
      SUBROUTINE get_rhs(integrand,rhsdum,bc_routine,bc_flag,r0_flag,
     $                   do_surf,surf_int,rmat_elim,cmat_elim)
      USE seam_storage_mod
      USE edge
      USE regularity

      CHARACTER(*), INTENT(IN) :: bc_flag,r0_flag
      LOGICAL, INTENT(IN) :: do_surf
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: rhsdum
      TYPE(global_matrix_type), DIMENSION(:), POINTER, OPTIONAL :: 
     $                          rmat_elim
      TYPE(complex_matrix_type), DIMENSION(:), POINTER, OPTIONAL :: 
     $                            cmat_elim

      INTEGER(i4) :: ibl,ibe,nside,imode,i_ri,iv,ivp,nqty,nfour,poly_deg
      INTEGER(i4) :: n_int,ndisc,ngrid,nqd,ncomb,iq,ib,ix,iy
      TYPE(vector_type), DIMENSION(:), POINTER :: vtmp1,vtmp2,vtmp3
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: ctmp1,ctmp2,ctmp3
      CHARACTER(8) :: flag
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE integrand(int,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the boundary condition routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE bc_routine(rhs,seam,flag,nqty,symm,svess)
        USE local
        USE edge_type_mod
        USE vector_type_mod
        INTEGER(i4), INTENT(IN) :: nqty
        TYPE(cvector_type), INTENT(INOUT) :: rhs
        TYPE(edge_type), INTENT(IN) :: seam
        CHARACTER(*), INTENT(IN) :: flag
        CHARACTER(*), INTENT(IN), OPTIONAL :: symm,svess
        END SUBROUTINE bc_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the surface integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE surf_int(sint,rb,tb,x,y,bigr,norm,
     $                      ijcell,alpha,dxdr,dydr,dxdz,dydz)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: sint 
        TYPE(rblock_type), INTENT(INOUT) :: rb
        TYPE(tblock_type), INTENT(INOUT) :: tb
        REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
        REAL(r8), DIMENSION(2), INTENT(IN) :: norm
        REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
        INTEGER(i4), DIMENSION(2) :: ijcell
        
        END SUBROUTINE surf_int
      END INTERFACE
c-----------------------------------------------------------------------
c     find the quantity and Fourier component dimensions.
c
c     if there are continuous basis functions, poly_deg is set according
c     to coefficient arrays for the continuous expansion.
c-----------------------------------------------------------------------
      ngrid=0
      poly_deg=0
      IF (ASSOCIATED(rhsdum(1)%arr)) THEN
        nqty=SIZE(rhsdum(1)%arr,1)
        nfour=SIZE(rhsdum(1)%arr,4)
        ngrid=1
        poly_deg=1
      ELSE IF (ASSOCIATED(rhsdum(1)%arrh)) THEN
        nqty=SIZE(rhsdum(1)%arrh,1)
        nfour=SIZE(rhsdum(1)%arrh,5)
      ELSE IF (ASSOCIATED(rhsdum(1)%arri)) THEN
        nqty=SIZE(rhsdum(1)%arri,1)
        nfour=SIZE(rhsdum(1)%arri,5)
      ELSE
        nqty=0
        nfour=0
      ENDIF

      IF (ASSOCIATED(rhsdum(1)%arrh)) THEN
        nside=SIZE(rhsdum(1)%arrh,2)
        poly_deg=nside+1
      ELSE
        nside=0
      ENDIF
c-----------------------------------------------------------------------
c     check if an auxiliary discontinuous field is used together with
c     the continuous field.  if so, its coefficients are placed in
c     arrtmp of rhsdum.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhsdum(1)%arrtmp)) THEN
        nqd=SIZE(rhsdum(1)%arrtmp,1)
        ndisc=SIZE(rhsdum(1)%arrtmp,2)
        IF (nfour==0) nfour=SIZE(rhsdum(1)%arrtmp,5)
      ELSE
        nqd=0
        ndisc=0
      ENDIF
c-----------------------------------------------------------------------
c     if there are no continuous fields, the discontinuous field
c     coefficients are placed in arri of rhstmp for compatibility with
c     the permanent storage for discontinuous fields.  reset poly_deg
c     according to the coefficient array for the discontinuous
c     expansion and multiply by -1 as a flag for the vector allocates
c     that appear below.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhsdum(1)%arri)) THEN
        n_int=SIZE(rhsdum(1)%arri,2)
        IF (poly_deg==0) poly_deg=-(NINT(SQRT(REAL(n_int)))-1_i4)
      ELSE
        n_int=0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for the appropriate integrand.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mpseudo=mpsq_block(ibl)
        ipseust=ipqst_block(ibl)
        ipseuen=ipqen_block(ibl)
        ncontb=4*ngrid+4*nside+nside**2
        IF (poly_deg>0) THEN
          ndiscb=ndisc
        ELSE
          ndiscb=n_int
        ENDIF
        CALL rblock_get_comp_rhs_q(rb(ibl),rhsdum(ibl),integrand)
      ENDDO
      DO ibl=nrbl+1,nbl
        mpseudo=mpsq_block(ibl)
        ipseust=ipqst_block(ibl)
        ipseuen=ipqen_block(ibl)
        ncontb=3
        ndiscb=0
        CALL tblock_get_comp_rhs_q(tb(ibl),rhsdum(ibl),integrand)
      ENDDO
c-----------------------------------------------------------------------
c     contributions from surface integrals for continuous fields,
c     where needed.
c-----------------------------------------------------------------------
      IF (do_surf) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          DO iv=1,seam(ibe)%nvert
            ivp=iv-1
            IF (ivp==0) ivp=seam(ibe)%nvert
            IF (seam(ibe)%expoint(iv).AND.seam(ibe)%expoint(ivp)) THEN
              IF (ibe<=nrbl) THEN
                CALL surface_rbl_comp_rhs(rb(ibe),rhsdum(ibe),surf_int,
     $                                    nqty,nfour,
     $                                    seam(ibe)%segment(iv)%intxys,
     $                                    seam(ibe)%segment(iv)%intxyp,
     $                                    seam(ibe)%segment(iv)%intxyn,
     $                                    seam(ibe)%segment(iv)%h_side,
     $                                    met_spl,nside+1_i4,geom)
              ELSE
                CALL surface_tbl_comp_rhs(tb(ibe),rhsdum(ibe),surf_int,
     $                                    nqty,nfour,
     $                                    seam(ibe)%segment(iv)%intxys,
     $                                    seam(ibe)%segment(iv)%intxyp,
     $                                    seam(ibe)%segment(iv)%intxyn,
     $                                    nside+1_i4,geom)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if this is a nonlinear computation without dealiasing, the
c     largest Fourier component must have a real equation.
c-----------------------------------------------------------------------
      IF (nonlinear.AND..NOT.dealiase) THEN
        DO imode=1,nmodes
          IF (nindex(imode)==nphi/2) THEN
            DO ibl=1,nbl
              CALL cvector_real_comp(rhsdum(ibl),imode)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if a matrix has been passed as an optional argument, it contains
c     information to eliminate the cell-centered data.  this must be
c     done before seaming.
c
c     note that if nqty=3, it is assumed that the real matrix is either
c     for the real r and z components and minus the imaginary phi
c     component, or it is for the imaginary r and z and real phi.
c     these cases must not have nqty present in the assignment calls.
c-PRE
c     no discontinuous fields for real solves at this point.
c-----------------------------------------------------------------------
      IF (PRESENT(rmat_elim).AND.n_int>0) THEN
        ALLOCATE(vtmp1(nrbl),vtmp2(nrbl))
        DO ibl=1,nrbl
          CALL vector_type_alloc(vtmp1(ibl),poly_deg,rb(ibl)%mx,
     $                           rb(ibl)%my,nqty)
          CALL vector_type_alloc(vtmp2(ibl),poly_deg,rb(ibl)%mx,
     $                           rb(ibl)%my,nqty)
        ENDDO
        mode_loop: DO imode=1,nfour
          DO i_ri=0,1
            IF (nqty==3) THEN
              IF (i_ri==0) THEN
                flag='r12mi3'
              ELSE
                flag='i12r3'
              ENDIF
            ELSE
              IF (i_ri==0) THEN
                flag='real'
              ELSE
                flag='imag'
              ENDIF
            ENDIF
            IF (keff(imode)==0) THEN
              IF (i_ri>0) CYCLE mode_loop
              flag='real'
            ENDIF
            DO ibl=1,nrbl
              CALL vector_assign_cvec(vtmp1(ibl),rhsdum(ibl),flag,imode)
            ENDDO

            CALL timer(timestart)
            CALL matelim_presolve(rmat_elim(imode),vtmp1,vtmp2,nqty)
            CALL timer(timeend)
            time_stcon=time_stcon+timeend-timestart

            DO ibl=1,nrbl
              CALL cvector_assign_vec(rhsdum(ibl),vtmp2(ibl),flag,imode)
            ENDDO
          ENDDO
        ENDDO mode_loop
        DO ibl=1,nrbl
          CALL vector_type_dealloc(vtmp1(ibl))
          CALL vector_type_dealloc(vtmp2(ibl))
        ENDDO
        DEALLOCATE(vtmp1,vtmp2)
c-----------------------------------------------------------------------
c     elimination for complex matrices--both real and imaginary parts
c     are done simultaneously.
c
c     if there is an auxiliary discontinuous field, the interior part
c     of the continuous field and the discontinuous field need to be
c     combined prior to the presolve call.  here, ctmp3 is used to swap 
c     memory space of the correct size.
c
c     note that the allocted ctmp1 and ctmp2 do not have arrtmp space.
c     this is checked by the vector-type assign operations, which skip
c     arrtmp operations when either vector does not have arrtmp space.
c-----------------------------------------------------------------------
      ELSE IF (PRESENT(cmat_elim).AND.(n_int>0.OR.ndisc>0)) THEN
        ALLOCATE(ctmp1(nrbl),ctmp2(nrbl))
        DO ibl=1,nrbl
          CALL vector_type_alloc(ctmp1(ibl),poly_deg,rb(ibl)%mx,
     $                           rb(ibl)%my,nqty)
          CALL vector_type_alloc(ctmp2(ibl),poly_deg,rb(ibl)%mx,
     $                           rb(ibl)%my,nqty)
        ENDDO

        IF (nqd>0) THEN
          ALLOCATE(ctmp3(nrbl))
          ncomb=nqty*n_int+nqd*ndisc
        ENDIF
        DO imode=1,nfour
          IF (nqty>0) THEN
            DO ibl=1,nrbl
              CALL cvector_2D_assign_cvec(ctmp1(ibl),rhsdum(ibl),imode)
            ENDDO
          ENDIF
          IF (nqd>0) THEN   !   swapping part
            DO ibl=1,nrbl
              ctmp3(ibl)%arrtmp=>ctmp1(ibl)%arri
              ctmp3(ibl)%arri=>ctmp2(ibl)%arri
              ALLOCATE(ctmp3(ibl)%arrh(ncomb,1,rb(ibl)%mx,rb(ibl)%my))
              ALLOCATE(ctmp3(ibl)%arrv(ncomb,1,rb(ibl)%mx,rb(ibl)%my))
              ctmp1(ibl)%arri=>ctmp3(ibl)%arrh
              ctmp2(ibl)%arri=>ctmp3(ibl)%arrv
              CALL cvector_2D_pack_cvec(rhsdum(ibl),ctmp1(ibl),
     $                                  nqty,n_int,nqd,ndisc,imode)
            ENDDO
          ENDIF

          CALL timer(timestart)
          CALL matelim_presolve(cmat_elim(imode),ctmp1,ctmp2,nqty)
          CALL timer(timeend)
          time_stcon=time_stcon+timeend-timestart

          IF (nqd>0) THEN   !   undo the swapping
            DO ibl=1,nrbl
              CALL cvector_2D_unpack_cvec(ctmp2(ibl),ctmp3(ibl),
     $                                    rhsdum(ibl),nqty,n_int,
     $                                    nqd,ndisc,imode)
              ! ctmp3%arri now holds int and rhsdum%arrtmp holds disc
              DEALLOCATE(ctmp1(ibl)%arri,ctmp2(ibl)%arri)
              ctmp1(ibl)%arri=>ctmp3(ibl)%arrtmp
              ctmp2(ibl)%arri=>ctmp3(ibl)%arri
              NULLIFY(ctmp3(ibl)%arrh,ctmp3(ibl)%arrv,
     $                ctmp3(ibl)%arri,ctmp3(ibl)%arrtmp)
            ENDDO
          ENDIF
          IF (nqty>0) THEN
            DO ibl=1,nrbl
              CALL cvector_assign_cvec2(rhsdum(ibl),ctmp2(ibl),imode)
            ENDDO
          ENDIF
        ENDDO
        DO ibl=1,nrbl
          CALL vector_type_dealloc(ctmp1(ibl))
          CALL vector_type_dealloc(ctmp2(ibl))
        ENDDO
        DEALLOCATE(ctmp1,ctmp2)
        IF (nqd>0) THEN
          DEALLOCATE(ctmp3)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     network to seams and impose external boundary conditions.
c     if blocks touch an R=0 border, apply regularity conditions.
c
c     return if there are no continuous fields.
c-----------------------------------------------------------------------
      IF (.NOT.ASSOCIATED(rhsdum(1)%arr)) RETURN
      DO ibl=1,nbl
        CALL edge_load_carr(rhsdum(ibl),nqty,1_i4,nfour,nside,seam(ibl))
      ENDDO
      CALL edge_network(nqty,nfour,nside,.false.)
      DO ibl=1,nbl
        CALL edge_unload_carr(rhsdum(ibl),nqty,1_i4,nfour,nside,
     $                        seam(ibl))
      ENDDO
      DO ibl=1,SIZE(r0block_list)
        ibe=r0block_list(ibl)
        CALL regular_vec(rhsdum(ibe),seam(ibe),r0_flag,nqty,nmodes,
     $                   nindex)
      ENDDO
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL bc_routine(rhsdum(ibe),seam(ibe),bc_flag,nqty)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_rhs
c-----------------------------------------------------------------------
c     subprogram 4. real_fe_postsolve.
c     when cell-interior data has been eliminated for a matrix solve,
c     the new interior data is generated through this routine.
c
c     rhs and sln are the right hand side and solution for a matrix
c     solve for single Fourier component.  old_cvec holds the solution
c     vector for all Fourier components at the start of this time
c     split.  crhs has the rhs for all Fourier components, where those
c     located in the cell interior have been multiplied by the inverse
c     of the interior to interior matrix.
c-----------------------------------------------------------------------
      SUBROUTINE real_fe_postsolve(mat_str,rhs,sln,old_cvec,crhs,
     $                             nq,imode,flag)

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: rhs,sln
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: old_cvec,crhs
      INTEGER(i4), INTENT(IN) :: nq,imode
      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     since the matrix equations are set-up for the change of a 
c     quantity, and the iterative solve is for the final value,
c     get the change from the iterative solution and the old vector.
c
c     A_ii**-1.b_i is stored in crhs%arri
c-----------------------------------------------------------------------
c-PRE
      DO ibl=1,nrbl
        CALL vector_assign_cvec(rhs(ibl),old_cvec(ibl),flag,imode)
        CALL vector_add(rhs(ibl),sln(ibl),v1fac=-1._r8)
        SELECT CASE(flag)
        CASE ('real','REAL')
          sln(ibl)%arri=crhs(ibl)%arri(1:nq,:,:,:,imode)
        CASE ('imag','IMAG')
          sln(ibl)%arri=AIMAG(crhs(ibl)%arri(1:nq,:,:,:,imode))
        CASE ('r12mi3')
          sln(ibl)%arri(1:2,:,:,:)=crhs(ibl)%arri(1:2,:,:,:,imode)
          sln(ibl)%arri(3,:,:,:)=-AIMAG(crhs(ibl)%arri(3,:,:,:,imode))
        CASE ('i12r3')
          sln(ibl)%arri(1:2,:,:,:)=
     $      AIMAG(crhs(ibl)%arri(1:2,:,:,:,imode))
          sln(ibl)%arri(3,:,:,:)=crhs(ibl)%arri(3,:,:,:,imode)
        CASE DEFAULT
          CALL nim_stop
     $      ('Real_fe_postsolve: flag '//flag//' not recognized.')
        END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     find x_i = A_ii**-1.b - A_ii**-1.A_io.x_o
c     where x_o is the change in the solution at grid vertices and
c     line segments, and x_i is the change in the interior solution.
c-----------------------------------------------------------------------
      CALL timer(timestart)
      CALL matelim_real_postsolve(mat_str,rhs,sln,nq)
      CALL timer(timeend)
      time_stcon=time_stcon+timeend-timestart
c-----------------------------------------------------------------------
c     for the interior nodes, add the old solution to the change to
c     find the new solution.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        SELECT CASE(flag)
        CASE ('real','REAL')
          sln(ibl)%arri=sln(ibl)%arri
     $                 +old_cvec(ibl)%arri(:,:,:,:,imode)
        CASE ('imag','IMAG')
          sln(ibl)%arri=sln(ibl)%arri
     $                 +AIMAG(old_cvec(ibl)%arri(:,:,:,:,imode))
        CASE ('r12mi3')
          sln(ibl)%arri(1:2,:,:,:)=sln(ibl)%arri(1:2,:,:,:)
     $      +old_cvec(ibl)%arri(1:2,:,:,:,imode)
          sln(ibl)%arri(3,:,:,:)=sln(ibl)%arri(3,:,:,:)
     $      -AIMAG(old_cvec(ibl)%arri(3,:,:,:,imode))
        CASE ('i12r3')
          sln(ibl)%arri(1:2,:,:,:)=sln(ibl)%arri(1:2,:,:,:)
     $      +AIMAG(old_cvec(ibl)%arri(1:2,:,:,:,imode))
          sln(ibl)%arri(3,:,:,:)=sln(ibl)%arri(3,:,:,:)
     $      +old_cvec(ibl)%arri(3,:,:,:,imode)
        END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE real_fe_postsolve
c-----------------------------------------------------------------------
c     subprogram 5. comp_fe_postsolve.
c     when cell-interior data has been eliminated for a matrix solve,
c     the new interior data is generated through this routine.
c
c     rhs and sln are the right hand side and solution for a matrix
c     solve for single Fourier component.  old_cvec holds the solution
c     vector for all Fourier components at the start of this time
c     split.  crhs has the rhs for all Fourier components, where those
c     located in the cell interior have been multiplied by the inverse
c     of the interior to interior matrix.
c-----------------------------------------------------------------------
      SUBROUTINE comp_fe_postsolve(mat_str,rhs,sln,old_cvec,crhs,
     $                             nq,imode,dsc)

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(INOUT) :: rhs,sln
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: old_cvec,crhs
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT), OPTIONAL :: dsc
      INTEGER(i4), INTENT(IN) :: nq,imode

      INTEGER(i4) :: ibl,nqd,ndisc,n_int,ncomb,iv,ib,iq,ix,iy,mxb,myb
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: ctmp
c-----------------------------------------------------------------------
c     first check if the eliminated data includes coefficients of 
c     fields that are discontinuous across element borders.  if so,
c     use the ctmp data structure to combine them with interior data.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(crhs(1)%arrtmp)) THEN
        ALLOCATE(ctmp(nrbl))
        nqd=SIZE(crhs(1)%arrtmp,1)
        ndisc=SIZE(crhs(1)%arrtmp,2)
        n_int=SIZE(crhs(1)%arri,2)
        ncomb=nq*n_int+nqd*ndisc
      ELSE
        nqd=0
        ndisc=0
      ENDIF
c-----------------------------------------------------------------------
c     since the matrix equations are set-up for the change of a 
c     quantity, and the iterative solve is for the final value,
c     get the change from the iterative solution and the old vector.
c
c     A_ii**-1.b_i is stored in crhs%arri
c
c     like get_rhs, pointer swapping is used to combine interior data
c     and discontinuous-field data before the matelim call, if
c     necessary.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL cvector_2D_assign_cvec(rhs(ibl),old_cvec(ibl),imode)
        CALL vector_add(rhs(ibl),sln(ibl),v1fac=-1._r8)
        IF (nqd==0) THEN
          sln(ibl)%arri=crhs(ibl)%arri(1:nq,:,:,:,imode)
        ELSE
          mxb=SIZE(rhs(ibl)%arri,3)
          myb=SIZE(rhs(ibl)%arri,4)
          ctmp(ibl)%arri=>sln(ibl)%arri
          ctmp(ibl)%arrtmp=>rhs(ibl)%arri
          ALLOCATE(ctmp(ibl)%arrh(ncomb,1,mxb,myb))
          ALLOCATE(ctmp(ibl)%arrv(ncomb,1,mxb,myb))
          sln(ibl)%arri=>ctmp(ibl)%arrh
          rhs(ibl)%arri=>ctmp(ibl)%arrv
          CALL cvector_2D_pack_cvec(crhs(ibl),sln(ibl),
     $                              nq,n_int,nqd,ndisc,imode)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     find x_i = A_ii**-1.b - A_ii**-1.A_io.x_o
c     where x_o is the change in the solution at grid vertices and
c     line segments, and x_i is the change in the interior solution.
c-----------------------------------------------------------------------
      CALL timer(timestart)
      CALL matelim_comp_postsolve(mat_str,rhs,sln,nq)
      CALL timer(timeend)
      time_stcon=time_stcon+timeend-timestart
c-----------------------------------------------------------------------
c     more data rearranging if needed for discontinuous fields.
c     discontinuous fields are added to the dsc structure.
c-----------------------------------------------------------------------
      IF (nqd>0) THEN
        DO ibl=1,nrbl
          IF (PRESENT(dsc)) THEN  !  save the discontinuous part in dsc
            CALL cvector_2D_unpack_add_cvec(sln(ibl),ctmp(ibl),dsc(ibl),
     $                                      nq,n_int,nqd,ndisc,imode)
          ELSE
            CALL cvector_2D_unpack_cvec2(sln(ibl),ctmp(ibl),
     $                                   nq,n_int,nqd,ndisc,.false.)
            ! disc part is not retained
          ENDIF
          DEALLOCATE(sln(ibl)%arri,rhs(ibl)%arri)
          sln(ibl)%arri=>ctmp(ibl)%arri
          rhs(ibl)%arri=>ctmp(ibl)%arrtmp
          NULLIFY(ctmp(ibl)%arrh,ctmp(ibl)%arrv,
     $            ctmp(ibl)%arri,ctmp(ibl)%arrtmp)
        ENDDO
        DEALLOCATE(ctmp)
      ENDIF
c-----------------------------------------------------------------------
c     for the interior nodes, add the old solution to the change to
c     find the new solution.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        sln(ibl)%arri=sln(ibl)%arri
     $               +old_cvec(ibl)%arri(:,:,:,:,imode)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE comp_fe_postsolve
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE finite_element_mod
