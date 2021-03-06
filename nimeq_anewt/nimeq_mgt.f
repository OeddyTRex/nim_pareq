c-----------------------------------------------------------------------
c     file nimeq_mgt.f:  contains management routines for finite
c     element computations performed for diagnostic calculations
c     in nimplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  gs_inter_block
c     1.  gssolve
c     2.  gsfree
c     3.  delstar
c     4.  gs_coil_read
c     5.  gs_surface_flux
c     6.  gs_work_init
c     7.  gs_find_nl
c     8.  calc_pres
c     9.  p_func
c     10. calc_f
c     11. f_func
c     12. f_adjust
c     13. pflux_range
c     14. jeq_poloidal
c     15. normalize_psi
c     16. calc_blen
c     17. plasma_cur
c     18. par_cont_init
c     19. equil_project
c     20. equil_plot
c     21. pflux_block
c     22. update_qfp
c     23. advect_btop
c     24. advect_solve
c     25. pfsmth_init
c-----------------------------------------------------------------------
c     0. module gs_inter_block
c     this module provides the interfaces for external programs that
c     are in this file.
c-----------------------------------------------------------------------
      MODULE gs_inter_block

c-----------------------------------------------------------------------
c     interface block for the calc_pres routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_pres(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_pres
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the calc_f routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_f(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_f
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the f_adjust subroutine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE f_adjust(plcar,r2ar,change_f1,dipdf)
        USE local
        IMPLICIT NONE
        REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plcar,r2ar
        REAL(r8), INTENT(OUT) :: dipdf
        LOGICAL, INTENT(IN) :: change_f1
        END SUBROUTINE f_adjust
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the gs_coil_read subroutine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE gs_coil_read(ncurf,curf)
        USE local
        IMPLICIT NONE
        INTEGER(i4), INTENT(OUT) :: ncurf
        REAL(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: curf
        END SUBROUTINE gs_coil_read
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the plasma_cur subroutine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE plasma_cur(muz,plcur,plcar,r2ar)
        USE local
        USE fields
        IMPLICIT NONE
        REAL(r8), INTENT(IN) :: muz
        REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plcar,r2ar
        REAL(r8), INTENT(OUT) :: plcur
        END SUBROUTINE plasma_cur
      END INTERFACE

      END MODULE gs_inter_block
c-----------------------------------------------------------------------
c     subprogram 1. gssolve
c     solve the Grad-Shafranov equation to find poloidal flux where the
c     flux along the surface is prescribed.
c
c     in this routine, the work arrays are used in the following way:
c       rwork1 -- psi/R**2 at domain surface, and 0 in interior
c                 (computed from an externally supplied B-normal)
c                 during nonlinear iterations.  normalized GS residual
c                 afterward.
c       rwork2 -- psi/R**2 in the interior and 0 at the surface
c                 (provided by the solution to GS with homogeneous
c                  Dirichlet boundary conditions).
c		  after the nonlinear iteration, rwork2 holds the
c		  total psi/R**2.
c       rwork3 -- holds FF' (F=R*B_phi) during the iterations.
c		  elsewhere, it just has F.
c       be_n0  -- holds the fields needed to determine Bpol for
c		  field-line tracing when btop_check is not none.  the
c		  4-vector is used to hold (lam,bigR,R,Z) and not B,
c		  itself.  if geom=="tor", bigR=R, otherwise bigR=1.
c       fllen  -- used to store field-line length and referenced by
c                 calc_pres and calc_f.
c
c     the poloidal flux function storage is pflux.  it is defined as the
c	physical poloidal flux divided by 2*pi in toroidal geometry and
c	physical flux divided by the periodic length in linear geometry.
c-----------------------------------------------------------------------
      SUBROUTINE gssolve(first_data,datflag,phii,ndcon)
      USE local
      USE fields
      USE input
      USE input_eq
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimeq_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_gmres_r2d
      USE matrix_storage_mod
      USE plot_data
      USE boundary
      USE nimeq_mod
      USE mpi_nim
      USE gs_inter_block
      USE iter_utils
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      REAL(r8), INTENT(IN) :: phii
      INTEGER(i4), INTENT(IN) :: ndcon

      INTEGER(i4) :: ibl,ibe,jphi,its,iq,nmodes_save
      INTEGER(i4) :: iters,maxiters
      REAL(r8) :: err,err2,er0
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
      LOGICAL :: converged,apply_blen,set_jdir,first_delst

      INTEGER(i4) :: ibasis,ix,iy,ipr,ibloc,ncurfil,ierror
      REAL(r8) :: dx,dy,b2,lam,tollin,p0,ref_len,oold,rtmp,afldl,gsuse
      REAL(r8) :: tol_incr=2._r8,tol_max=1.e-5_r8
      REAL(r8), DIMENSION(3) :: vec
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: curfil
      CHARACTER(1) :: bcsymm="n"
c-----------------------------------------------------------------------
c     interface block for eqdot_delst
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE eqdot_delst(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE eqdot_delst
      END INTERFACE
c-----------------------------------------------------------------------
c     input variables are now included in the input_eq module.  local
c     variables are initialized here.
c-----------------------------------------------------------------------
      iters=0
      converged=.FALSE.
      set_jdir=.TRUE.
      maxiters=eq_iters     
      smallnum=SQRT(TINY(smallnum))
      er0=1._r8
      err2=1._r8
      relerr=10._r8
      gsuse=gscenter
      tollin=gsh_tol
c-----------------------------------------------------------------------
c     rlencl indicates what fraction of the initial range of the
c     scalar-field is to be considered closed.  the advective algorithm
c     works better with something closer to 0 than to 1.  also set
c     parameters for edge-profile smoothing.
c-----------------------------------------------------------------------
      IF (btop_check=="passive adv") rlencl=0.1_r8
      IF (psin_smthpp<1._r8.OR.psin_smthfp<1._r8) CALL pfsmth_init
c-----------------------------------------------------------------------
c     set the symmetry flag for boundary conditions according to
c     the input parameter symm_region.  if this is the top region of
c     a vertically symmetric equilibrium, bcsymm is "b" for bottom
c     boundary.
c-----------------------------------------------------------------------
      IF (symm_region=="top") THEN
        bcsymm="b"
      ELSE IF (symm_region=="bottom") THEN
        bcsymm="t"
      ENDIF
c-----------------------------------------------------------------------
c     initialize work arrays.
c-----------------------------------------------------------------------
      CALL gs_work_init(ref_len)
c-----------------------------------------------------------------------
c     read coil information from a binary file.
c-----------------------------------------------------------------------
      CALL gs_coil_read(ncurfil,curfil)
c-----------------------------------------------------------------------
c     create arrays that contain all coil information.
c-----------------------------------------------------------------------
      ncoil_tot=neqcoil+ncurfil
      ALLOCATE(allcoil_r(ncoil_tot),allcoil_z(ncoil_tot),
     $         allcoil_i(ncoil_tot))
      allcoil_r=(/eqcoil_r(1:neqcoil),curfil(1,1:ncurfil)/)
      allcoil_z=(/eqcoil_z(1:neqcoil),curfil(2,1:ncurfil)/)
      allcoil_i=(/eqcoil_i(1:neqcoil),curfil(3,1:ncurfil)/)
c-----------------------------------------------------------------------
c     initialize the poloidal flux by inverting the del-star operator
c     using current density and boundary flux from the dump file.
c-----------------------------------------------------------------------
      IF (use_delst_init) THEN
        CALL gs_surface_flux(psio_min,psio_max,.false.)
        oflux=psio_max-psio_min
        first_delst=.true.
        CALL delstar(first_delst,datflag,phii,ndcon)

        DO ibl=1,nbl
          CALL vector_add(rwork2(ibl),rwork1(ibl))   !  set lambda
          CALL vector_assignq_vec(be_n0(ibl),rwork2(ibl),nqty=1_i4,
     $                            nstart1=1_i4)
          CALL vector_add(rwork2(ibl),rwork1(ibl),v2fac=-1._r8)
          IF (btop_check/="none") THEN
            fllen(ibl)=0._r8 ! compute lengths in first pass
          ELSE
            fllen(ibl)=1._r8
          ENDIF

          pres_eq(ibl)=0._r8
          CALL vector_assignq_vec(rwork3(ibl),ja_eq(ibl),1_i4,1_i4,3_i4)
          CALL vector_mult(rwork3(ibl),-mu0)
          rwork3(ibl)%arr(1,:,:)=rwork3(ibl)%arr(1,:,:)*
     $      be_n0(ibl)%arr(2,:,:)**2
          IF (ASSOCIATED(rwork3(ibl)%arrh).AND.ibl<=nrbl) THEN
            rwork3(ibl)%arrh(1,:,:,:)=rwork3(ibl)%arrh(1,:,:,:)*
     $        be_n0(ibl)%arrh(2,:,:,:)**2
            rwork3(ibl)%arrv(1,:,:,:)=rwork3(ibl)%arrv(1,:,:,:)* 
     $        be_n0(ibl)%arrv(2,:,:,:)**2
            rwork3(ibl)%arri(1,:,:,:)=rwork3(ibl)%arri(1,:,:,:)* 
     $        be_n0(ibl)%arri(2,:,:,:)**2
          ENDIF
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%rwork3,rb(ibl)%qrwork1,
     $                            rb(ibl))
            rb(ibl)%qpres_eq%qpf=0._r8
          ELSE
            CALL tblock_qp_update(tb(ibl)%rwork3,tb(ibl)%qrwork1,
     $                            tb(ibl))
            tb(ibl)%qpres_eq%qpf=0._r8
          ENDIF
        ENDDO
        SELECT CASE(btop_check)
        CASE ("beq-trace")
          CALL calc_blen(ref_len,stiff_solver)
        CASE ("passive adv")
          integrand_flag="forward"
          CALL advect_solve(bcsymm,ref_len)
        END SELECT
        apply_blen=.true.
      ELSE
        DO ibl=1,nbl
          fllen(ibl)=1._r8 ! treats all field lines as closed to start
        ENDDO
        apply_blen=.false.
        psic_min=0._r8
        psic_max=0._r8
        sep_flux=0._r8
        cflux=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     nimeq solves work with only 1 Fourier component (n=0).
c-----------------------------------------------------------------------
      nmodes_save=nmodes
      nmodes=1
c-----------------------------------------------------------------------
c     allocate vectors needed for computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,1_i4)
        sln(ibl)=0._r8
        crhs(ibl)=0._r8  !  needed for cvector_assign_vec below
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         1_i4,1_i4) 
        sln(ibl)=0._r8
        crhs(ibl)=0._r8  !  needed for cvector_assign_vec below
      ENDDO
c-----------------------------------------------------------------------
c     determine the inhomogeneous part of the poloidal flux along the
c     border and save it (divided by 2*pi*R**2 in toroidal geometry
c     or by per_length in linear) in rwork1.
c-----------------------------------------------------------------------
      CALL gs_surface_flux(psio_min,psio_max,.true.)
      oflux=psio_max-psio_min
c-----------------------------------------------------------------------
c     start the nonlinear iteration loop, which proceeds by:
c     1) finding the residual
c        a. update poloidal flux (and F parameter if needed)
c        b. check field-line length if needed
c        c. compute FF' and P'
c        d. adjust F coefficients for I if needed
c        e. find residual vector through FE operations
c     2) check if convergence criterion is met
c     3) solve the linearized system
c-----------------------------------------------------------------------
      non_it: DO
c-----------------------------------------------------------------------
c       1) a. if needed, update the F parameter, then update the
c             poloidal flux function from the solution for psi/R**2
c-----------------------------------------------------------------------
        IF (dellam>0._r8.AND.ip_equil/=0._r8.AND.relerr<1._r8.AND.
     $      ip_weight>0._r8) THEN
          afldl=0._r8
          DO ibl=1,nbl
            CALL vector_assignq_vec(rwork3(ibl),diperr(ibl),1_i4,
     $                              nstart1=1_i4,nstart2=2_i4)
            CALL iter_dot(afldl,rwork3(ibl),sln(ibl),seam(ibl),ibl,nrbl,
     $                    poly_degree,.true.)
            IF (nprocs>1) THEN
              CALL mpi_allreduce(afldl,rtmp,1,mpi_nim_real,mpi_sum,
     $             comm_layer,ierror)
              afldl=rtmp
            ENDIF
          ENDDO

          SELECT CASE (f_model)
          CASE("quad_closed")
            f1_closed=f1_closed-(afldl+iperr)/(ip_weight*dipdf)
          CASE("cubic_closed")
            f_axis=f_axis-(afldl+iperr)/(ip_weight*dipdf)
          CASE("linear")
            dfdpsi=dfdpsi-(afldl+iperr)/(ip_weight*dipdf)
          END SELECT
        ENDIF

        DO ibl=1,nbl
          CALL vector_add(sln(ibl),rwork2(ibl),v1fac=gsuse,v2fac=1._r8)
          rwork2(ibl)=sln(ibl)
          CALL vector_add(sln(ibl),rwork1(ibl))
          CALL pflux_block(sln(ibl),ibl)
        ENDDO
c-----------------------------------------------------------------------
c       1) b. update data for the field tracing, then call the
c             field-line computation.
c-----------------------------------------------------------------------
        IF (iters==3.AND.btop_check/="none".AND.
     $      .NOT.apply_blen) THEN
          apply_blen=.true.
          DO ibl=1,nbl
            fllen(ibl)=0._r8
          ENDDO
        ENDIF
        IF (btop_check/="none".AND.apply_blen) THEN
          DO ibl=1,nbl
            CALL vector_assignq_vec(be_n0(ibl),sln(ibl),nqty=1_i4,
     $                              nstart1=1_i4)
          ENDDO
          SELECT CASE(btop_check)
          CASE ("beq-trace")
            CALL calc_blen(ref_len,stiff_solver)
          CASE ("passive adv")
            integrand_flag="forward"
            CALL advect_solve(bcsymm,ref_len)
          END SELECT
        ENDIF
c-----------------------------------------------------------------------
c       1) c. if dellam>0, evaluate a centered difference of the
c             nonlinear terms about the current value of lambda for use
c             in the linear solves.  note that field-line lengths are
c             not varied.
c
c       the pflux_block calls should be part of the linearization, but
c       they seem to slow convergence when using ip_equil.
c-----------------------------------------------------------------------
        IF (dellam>0._r8) THEN
          IF (relerr<1._r8) THEN
            DO ibl=1,nbl
              crhs(ibl)=0.5_r8*dellam
              CALL dirichlet_rhs(crhs(ibl),seam(ibl),'all',1_i4,bcsymm)
              CALL vector_assign_cvec(rhs(ibl),crhs(ibl),'REAL',1_i4)
              CALL vector_add(sln(ibl),rhs(ibl))
              IF (ip_equil==0._r8) CALL pflux_block(sln(ibl),ibl)
            ENDDO
            CALL gs_findnl("no f1 update",iters,set_jdir)
            DO ibl=1,nbl
              nd_n0(ibl)=rwork3(ibl)      !  F*F' at incremented lambda
              pres_n0(ibl)=pres_eq(ibl)   !  P'   at incremented lambda
              CALL vector_add(sln(ibl),rhs(ibl),v2fac=-2._r8)
              IF (ip_equil==0._r8) CALL pflux_block(sln(ibl),ibl)
            ENDDO
            CALL gs_findnl("no f1 update",iters,set_jdir)
            DO ibl=1,nbl  !  form delta(F*F') and delta(P')
              CALL vector_add(nd_n0(ibl),rwork3(ibl),v2fac=-1._r8)
              CALL vector_add(pres_n0(ibl),pres_eq(ibl),v2fac=-1._r8)
              IF (ibl<=nrbl) THEN
                CALL rblock_qp_update(rb(ibl)%nd_n0,
     $                                rb(ibl)%qnd_n0,rb(ibl))
                CALL rblock_qp_update(rb(ibl)%pres_n0,
     $                                rb(ibl)%qpres_n0,rb(ibl))
              ELSE
                CALL tblock_qp_update(tb(ibl)%nd_n0,
     $                                tb(ibl)%qnd_n0,tb(ibl))
                CALL tblock_qp_update(tb(ibl)%pres_n0,
     $                                tb(ibl)%qpres_n0,tb(ibl))
              ENDIF
              CALL vector_add(sln(ibl),rhs(ibl)) ! reset original lambda
              IF (ip_equil==0._r8) CALL pflux_block(sln(ibl),ibl)
            ENDDO
            CALL gs_findnl("dFFdf only",iters,set_jdir)
c-----------------------------------------------------------------------
c           find the matrix elements that couple the F parameter and
c           lambda at the nodes to linearize the Ip convergence.  these
c           matrix elements are labeled A_lam_f and A_f_lam (diperr
c           quantity indices 1 & 2, respectively) in comments.
c-----------------------------------------------------------------------
            IF (ip_equil/=0._r8.AND.ip_weight>0._r8) THEN
              DO ibl=1,nrbl
                CALL rblock_get_rhs(rb(ibl),diperr(ibl),err_ip,2_i4)
                CALL edge_load_arr(diperr(ibl),2_i4,
     $                             poly_degree-1_i4,seam(ibl))
              ENDDO 
              DO ibl=nrbl+1,nbl
                CALL tblock_get_rhs(tb(ibl),diperr(ibl),err_ip,2_i4)
                CALL edge_load_arr(diperr(ibl),2_i4,0_i4,seam(ibl))
              ENDDO 
              CALL edge_network(2_i4,0_i4,poly_degree-1_i4,.false.)
              DO ibl=1,nbl
                CALL edge_unload_arr(diperr(ibl),2_i4,poly_degree-1_i4,
     $                               seam(ibl))
              ENDDO
              DO ibl=1,SIZE(exblock_list)
                ibe=exblock_list(ibl)
                crhs(ibe)=0._r8
                CALL cvector_assign_vec(crhs(ibe),diperr(ibe),'REAL',
     $                                  1_i4,nqty=1_i4)
                CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4,
     $                             bcsymm)
                CALL vector_assign_cvec(diperr(ibe),crhs(ibe),'REAL',
     $                                  1_i4,nqty=1_i4)
              ENDDO
            ENDIF
c-----------------------------------------------------------------------
c         just zero-out the arrays for the linearized FF' and P' terms
c         used by the integrand.
c-----------------------------------------------------------------------
          ELSE
            DO ibl=1,nrbl
              rb(ibl)%qnd_n0%qpf=0._r8
              rb(ibl)%qpres_n0%qpf=0._r8
              diperr(ibl)=0._r8
            ENDDO
            DO ibl=nrbl+1,nbl
              tb(ibl)%qnd_n0%qpf=0._r8
              tb(ibl)%qpres_n0%qpf=0._r8
              diperr(ibl)=0._r8
            ENDDO
            CALL gs_findnl("update f1",iters,set_jdir)
          ENDIF
c-----------------------------------------------------------------------
c       otherwise, find and save the coefficients for the
c       nonlinear terms with updates for matching ip_equil.
c-----------------------------------------------------------------------
        ELSE
          CALL gs_findnl("update f1",iters,set_jdir)
        ENDIF
c-----------------------------------------------------------------------
c       1) d. find the plasma current contained in the domain.  also set
c             the scaled error in plasma current.
c-----------------------------------------------------------------------
        CALL plasma_cur(mu0,pli,curarr,r2arr)
        iperr=ip_weight*(pli-ip_equil)
c-----------------------------------------------------------------------
c       1) e. project -(delstar(psi_tot)+FF'+mu0*R**2*P') onto
c             the basis functions as the residual vector, then
c             evaluate just the delstar part for normalization.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            integrand_flag="all terms"
            CALL rblock_get_rhs(rb(ibl),rhs(ibl),err_gs,1_i4)
            integrand_flag="delstar only"
            CALL rblock_get_rhs(rb(ibl),sln(ibl),err_gs,1_i4)
          ELSE
            integrand_flag="all terms"
            CALL tblock_get_rhs(tb(ibl),rhs(ibl),err_gs,1_i4)
            integrand_flag="delstar only"
            CALL tblock_get_rhs(tb(ibl),sln(ibl),err_gs,1_i4)
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       network block seams
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'REAL',1_i4)
          CALL cvector_assign_vec(crhs(ibl),sln(ibl),'IMAG',1_i4)
        ENDDO
        CALL edge_network(1_i4,1_i4,poly_degree-1_i4,.TRUE.)
c-----------------------------------------------------------------------
c       with fixed-boundary solves, the residual is 0 along the
c       exterior boundary.  [regularity routines do nothing  for n=0
c       scalars, so skip them.]
c-----------------------------------------------------------------------
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4,bcsymm)
        ENDDO
        DO ibl=1,nbl
          CALL vector_assign_cvec(rhs(ibl),crhs(ibl),'REAL',1_i4)
          CALL vector_assign_cvec(sln(ibl),crhs(ibl),'IMAG',1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       find the relative error.
c-----------------------------------------------------------------------
        IF (iters>0) THEN
          CALL iter_2_norm(sln,er0,1.0_r8,nbl,nrbl,poly_degree,.TRUE.)
          CALL iter_2_norm(rhs,err2,1._r8,nbl,nrbl,poly_degree,.TRUE.)
          relerr=err2/(er0*linffp)
          IF (relerr<1._r8.AND.dellam>0._r8) THEN
            gsuse=1._r8
          ELSE
            gsuse=gscenter
          ENDIF

          IF (ip_equil/=0._r8) err2=SQRT(err2**2+iperr**2)
          IF (write_iters.AND.node==0) THEN
            WRITE(*,'(a,i4,2(a,es11.4),a,i4,a,es11.4)') " GS it ",iters,
     $        " ||delst(psi)||=",er0," ||resid||=",err2,
     $        " CG its=",its ," Ip=",pli
          ENDIF
          IF (err2/er0<gsh_tol) converged=.true.
          IF (f_model=="quad_closed".AND.iters<2) converged=.false.
        ENDIF
c-----------------------------------------------------------------------
c       2) exit loop if convergence criterion is met.  before doing so,
c          combine rwork1 and rwork2 to have the full psi/R**2
c          distribution in one place.  also save the residual in rwork1.
c-----------------------------------------------------------------------
        IF (converged.OR.iters==maxiters) THEN
          DO ibl=1,nbl
            CALL vector_add(rwork2(ibl),rwork1(ibl))
            rwork1(ibl)=rhs(ibl)
            CALL vector_mult(rwork1(ibl),1._r8/er0)
          ENDDO

          IF (ip_equil/=0._r8.AND.node==0) THEN
            IF (f_model=="quad_closed") THEN
              WRITE(*,'(a,es12.5)') " Final f1_closed=",f1_closed
            ELSE IF (f_model=="cubic_closed") THEN
              WRITE(*,'(a,es12.5)') " Final f_axis=",f_axis
            ENDIF
          ENDIF

          IF (iters==maxiters) WRITE(*,*) "no convergence"

          EXIT non_it
        ENDIF
        iters=iters+1
c-----------------------------------------------------------------------
c       alter the rhs array as part of the matrix partitioning when
c       performing newton iteration with a targeted plasma current.
c-----------------------------------------------------------------------
        IF (dellam>0._r8.AND.ip_equil/=0._r8.AND.relerr<1._r8.AND.
     $      ip_weight>0._r8) THEN
          DO ibl=1,nbl
            CALL vector_assignq_vec(sln(ibl),diperr(ibl),nqty=1_i4,
     $                              nstart1=1_i4,nstart2=1_i4)
            CALL vector_add(rhs(ibl),sln(ibl),
     $                      v2fac=-iperr/(ip_weight*dipdf))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       3) invert the approximate Jacobian matrix to find
c          delta(lambda)=delta(psi/R**2) in the interior.
c
c          symmetry conditions are used for vertically unstable
c          equilibria where the solves are sensitive.
c
c          the norm-adjusting linear tolerance avoids progressively
c          restrictive absolute tolerance as the nonlinear resid, hence
c          the rhs for this system, shrinks.  see the team meeting
c          presentation from 4/2/16.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          sln(ibl)=0._r8
          CALL vector_mult(rhs(ibl),-1._r8)
        ENDDO

        CALL iter_gmr_r2d_mfsolve(eqdot_delst,delstar_mat,delstar_fac,
     $                            sln,rhs,1_i4,tollin,linmaxit,
     $                            nimeq_solver,err,its,seed)
        IF (err>tollin) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Delstar: no convergence: ',
     $          its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        tollin=MIN(tol_max,tollin*tol_incr)
      ENDDO non_it

      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     create the expansions of the equilibrium fields that are used
c     in nimrod.  also write output for plotting.
c-----------------------------------------------------------------------
      CALL equil_project
      CALL equil_plot(first_data,ndcon)
      nmodes=nmodes_save
c-----------------------------------------------------------------------
c     transfer equilibrium fields to n=0 if requested.
c-----------------------------------------------------------------------
      IF (nimeq_tx.AND.nonlinear) THEN
        CALL eq_swap(nmodes,keff,geom,eq_flow,ndeq_notx,rbph_notx)
        IF (.NOT.coil_tx.AND.geom=='tor'.AND.ncoil_tot>0)
     $    CALL replace_coil(nmodes,keff,ncoil_tot,allcoil_r,
     $                      allcoil_z,allcoil_i)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gssolve

c-----------------------------------------------------------------------
c     subprogram 2. gsfree
c     solve the Grad-Shafranov equation to find poloidal flux for
c     'free-boundary' computations, where flux along the surface of the
c     domain is consistent with external vacuum and coils.  this version
c     does not use a separate surface component.  it solves throughout
c     the domain.
c
c     in this routine, the work arrays are used in the following way:
c
c       rwork1 -- psi/R**2 at domain surface, and 0 in interior
c                 (computed from an externally supplied B-normal)
c                 during initialization.  just a work array later.
c       rwork2 -- psi/R**2 in the interior and along the surface.
c       rwork3 -- holds FF' (F=R*B_phi) during the iterations.
c		  elsewhere, it just has F.
c       be_n0  -- holds the fields needed to determine Bpol for
c		  field-line tracing when btop_check is not none.  the
c		  4-vector is used to hold (lam,bigR,R,Z) and not B,
c		  itself.  if geom=="tor", bigR=R, otherwise bigR=1.
c       fllen  -- used to store field-line length and referenced by
c                 calc_pres and calc_f.
c
c     the poloidal flux function storage is pflux.  it is defined as the
c	physical poloidal flux divided by 2*pi in toroidal geometry and
c	physical flux divided by the periodic length in linear geometry.
c
c     this routine uses the "matrix-free" solver, so static condensation
c       is not applied prior to sending the linear system to the solver.
c-----------------------------------------------------------------------
      SUBROUTINE gsfree(first_data,datflag,phii,ndcon)
      USE local
      USE fields
      USE input
      USE input_eq
      USE global
      USE time
      USE rblock
      USE tblock
      USE surface
      USE nimeq_ints
      USE nimeq_srf_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_gmres_r2d
      USE matrix_storage_mod
      USE plot_data
      USE boundary
      USE nimeq_mod
      USE nimeq_free
      USE mpi_nim
      USE gs_inter_block
      USE iter_utils
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      REAL(r8), INTENT(IN) :: phii
      INTEGER(i4), INTENT(IN) :: ndcon

      INTEGER(i4) :: ibl,ibe,jphi,its,iq,nmodes_save
      INTEGER(i4) :: iters,maxiters
      REAL(r8) :: err,err2,er0
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
      LOGICAL :: converged,apply_blen,set_jdir,first_delst

      INTEGER(i4) :: ibasis,ix,iy,ipr,ibloc,ncurfil,ierror,iv,ivp
      REAL(r8) :: dx,dy,b2,lam,p0,ref_len,oold,rtmp,exmn,exmx,tollin,
     $            afldl,gsuse
      REAL(r8) :: tol_incr=2._r8,tol_max=1.e-5_r8
      REAL(r8), DIMENSION(3) :: vec
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: curfil
      CHARACTER(1) :: bcsymm="n"
    
      TYPE(vector_type), DIMENSION(:), POINTER :: sl2
      TYPE(vector_type), DIMENSION(:), ALLOCATABLE :: rt2
c-----------------------------------------------------------------------
c     interface block for eqdot_delstlj
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE eqdot_delstlj(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE eqdot_delstlj
      END INTERFACE
c-----------------------------------------------------------------------
c     input variables are now included in the input_eq module.  local
c     variables are initialized here.
c-----------------------------------------------------------------------
      iters=0
      converged=.FALSE.
      set_jdir=.TRUE.
      maxiters=eq_iters     
      smallnum=SQRT(TINY(smallnum))
      bweight=delstlj_mat%diag_scale  !  used for border-lambda weight
      er0=1._r8
      err2=1._r8
      relerr=10._r8
      gsuse=gscenter
      tollin=gsh_tol
c-----------------------------------------------------------------------
c     rlencl indicates what fraction of the initial range of the
c     scalar-field is to be considered closed.  the advective algorithm
c     works better with something closer to 0 than to 1.  also set
c     parameters for edge-profile smoothing.
c-----------------------------------------------------------------------
      IF (btop_check=="passive adv") rlencl=0.1_r8
      IF (psin_smthpp<1._r8.OR.psin_smthfp<1._r8) CALL pfsmth_init
c-----------------------------------------------------------------------
c     set the symmetry flag for boundary conditions according to
c     the input parameter symm_region.  if this is the top region of
c     a vertically symmetric equilibrium, bcsymm is "b" for bottom
c     boundary.
c
c-PRE further development is needed for using symm_region with free-
c     boundary computations.
c-----------------------------------------------------------------------
      IF (symm_region=="top") THEN
        bcsymm="b"
      ELSE IF (symm_region=="bottom") THEN
        bcsymm="t"
      ENDIF
c-----------------------------------------------------------------------
c     initialize work arrays.
c-----------------------------------------------------------------------
      CALL gs_work_init(ref_len)
c-----------------------------------------------------------------------
c     read coil information from a binary file.
c-----------------------------------------------------------------------
      CALL gs_coil_read(ncurfil,curfil)
c-----------------------------------------------------------------------
c     create arrays that contain all coil information, and initialize
c     the surface-flux computations for free-boundary equilibria.
c-----------------------------------------------------------------------
      ncoil_tot=neqcoil
      IF (nrfbcoil>0) THEN
        fst_rfbc=ncoil_tot+1
        ncoil_tot=ncoil_tot+nrfbcoil
      ENDIF
      IF (nvfbcoil>0) THEN
        fst_vfbc=ncoil_tot+1
        ncoil_tot=ncoil_tot+nvfbcoil
      ENDIF
      ncoil_tot=ncoil_tot+ncurfil
      ALLOCATE(allcoil_r(ncoil_tot),allcoil_z(ncoil_tot),
     $         allcoil_i(ncoil_tot))
      allcoil_r=(/eqcoil_r(1:neqcoil),rfbcoil_r(1:nrfbcoil),
     $            vfbcoil_r(1:nvfbcoil),curfil(1,1:ncurfil)/)
      allcoil_z=(/eqcoil_z(1:neqcoil),rfbcoil_z(1:nrfbcoil),
     $            vfbcoil_z(1:nvfbcoil),curfil(2,1:ncurfil)/)
      allcoil_i=(/eqcoil_i(1:neqcoil),rfbcoil_i(1:nrfbcoil)*rfbc_curr,
     $            vfbcoil_i(1:nvfbcoil)*vfbc_curr,curfil(3,1:ncurfil)/)
      IF (.NOT.ALLOCATED(psi_from_j)) CALL nimeq_free_init
c-----------------------------------------------------------------------
c     load mu0*ja_eq into the quadrature-point storage for computing
c     surface flux.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_mult(ja_eq(ibl),mu0)
        CALL rblock_qp_update(rb(ibl)%ja_eq,rb(ibl)%qja_eq,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_mult(ja_eq(ibl),mu0)
        CALL tblock_qp_update(tb(ibl)%ja_eq,tb(ibl)%qja_eq,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     initialize the poloidal flux by inverting the del-star operator
c     using current density and boundary flux from the dump file.
c     here the interior and suface fluxes are separate (rwork2 and
c     rwork1, respectively).  use nimeq_free_eval to set the surface
c     flux values before calling delstar.
c-----------------------------------------------------------------------
      IF (use_delst_init) THEN
        CALL nimeq_free_eval(rwork1,1_i4,poly_degree,1._r8,
     $                       psio_max,psio_min)
        oflux=psio_max-psio_min
        first_delst=.true.

        DO ibl=1,nbl
          CALL vector_mult(ja_eq(ibl),1._r8/mu0) ! delstar needs J
        ENDDO
        CALL delstar(first_delst,'surf set',phii,ndcon)
        DO ibl=1,nbl
          CALL vector_mult(ja_eq(ibl),mu0) ! gsfree needs mu0*J
        ENDDO

        DO ibl=1,nbl
          IF (geom=='tor') THEN                      !  set flux values
            CALL vector_mult(pflux(ibl),1._r8/twopi)
          ELSE
            CALL vector_mult(pflux(ibl),1._r8/per_length)
          ENDIF
          CALL vector_add(rwork2(ibl),rwork1(ibl))   !  set lambda
          CALL vector_assignq_vec(be_n0(ibl),rwork2(ibl),nqty=1_i4,
     $                            nstart1=1_i4)
          IF (btop_check/="none") THEN
            fllen(ibl)=0._r8 ! compute lengths in first pass
          ELSE
            fllen(ibl)=1._r8
          ENDIF

          pres_eq(ibl)=0._r8
          CALL vector_assignq_vec(rwork3(ibl),ja_eq(ibl),1_i4,1_i4,3_i4)
          CALL vector_mult(rwork3(ibl),-1._r8)
          rwork3(ibl)%arr(1,:,:)=rwork3(ibl)%arr(1,:,:)*
     $      be_n0(ibl)%arr(2,:,:)**2
          IF (ASSOCIATED(rwork3(ibl)%arrh).AND.ibl<=nrbl) THEN
            rwork3(ibl)%arrh(1,:,:,:)=rwork3(ibl)%arrh(1,:,:,:)*
     $        be_n0(ibl)%arrh(2,:,:,:)**2
            rwork3(ibl)%arrv(1,:,:,:)=rwork3(ibl)%arrv(1,:,:,:)* 
     $        be_n0(ibl)%arrv(2,:,:,:)**2
            rwork3(ibl)%arri(1,:,:,:)=rwork3(ibl)%arri(1,:,:,:)* 
     $        be_n0(ibl)%arri(2,:,:,:)**2
          ENDIF
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%rwork3,rb(ibl)%qrwork1,
     $                            rb(ibl))
            rb(ibl)%qpres_eq%qpf=0._r8
          ELSE
            CALL tblock_qp_update(tb(ibl)%rwork3,tb(ibl)%qrwork1,
     $                            tb(ibl))
            tb(ibl)%qpres_eq%qpf=0._r8
          ENDIF
        ENDDO
c-PRE check that delst_init still works
c       SELECT CASE(btop_check)
c       CASE ("beq-trace")
c         CALL calc_blen(ref_len,stiff_solver)
c       CASE ("passive adv")
c         integrand_flag="forward"
c         CALL advect_solve(bcsymm,ref_len)
c       END SELECT
c       CALL pflux_range(btop_check,set_jdir)
        apply_blen=.true.
c-----------------------------------------------------------------------
c     if poloidal flux is not computed from existing J_phi, find the
c     range of flux values from external coils only.
c-----------------------------------------------------------------------
      ELSE
        apply_blen=.false.
        psic_min=0._r8
        psic_max=0._r8
        sep_flux=0._r8
        cflux=0._r8
        CALL gs_surface_flux(psio_min,psio_max,.true.)
        oflux=psio_max-psio_min
        DO ibl=1,nbl
          rwork2(ibl)=rwork1(ibl) ! copies surface-flux values
          fllen(ibl)=1._r8 ! treats all field lines as closed to start
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     nimeq solves work with only 1 Fourier component (n=0).
c-----------------------------------------------------------------------
      nmodes_save=nmodes
      nmodes=1
c-----------------------------------------------------------------------
c     allocate vectors needed for computation.  2-vectors are used
c     to find expansions of lambda and mu0*J_phi/R simultaneously during
c     free-boundary computations.
c-----------------------------------------------------------------------
      ALLOCATE(sl2(nbl),rt2(nbl))
      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,2_i4)
        CALL vector_type_alloc(sl2(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,2_i4)
        CALL vector_type_alloc(rt2(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,2_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,2_i4,1_i4)
        sl2(ibl)=0._r8
        rt2(ibl)=0._r8
        crhs(ibl)=0._r8  !  needed for cvector_assign_vec below
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,2_i4)
        CALL vector_type_alloc(sl2(ibl),1_i4,tb(ibl)%mvert,0_i4,2_i4)
        CALL vector_type_alloc(rt2(ibl),1_i4,tb(ibl)%mvert,0_i4,2_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         2_i4,1_i4) 
        sl2(ibl)=0._r8
        rt2(ibl)=0._r8
        crhs(ibl)=0._r8  !  needed for cvector_assign_vec below
      ENDDO
c-----------------------------------------------------------------------
c     start the nonlinear iteration loop, which proceeds by:
c     1) finding the residual
c        a. check field-line length if needed
c        b. compute FF' and P' and adjust F coefficients for I if needed
c        c. find the plasma current
c        d. find residual vector through FE operations
c     2) check if convergence criterion is met
c     3) find the linearized coefficients for approximate-Newton
c     4) solve the linearized system
c     5) update poloidal flux (and F parameter if needed)
c-----------------------------------------------------------------------
      non_it: DO
        IF (iters==0.AND..NOT.use_delst_init) set_jdir=.true.
c-----------------------------------------------------------------------
c       1) a. update the data for the field tracing, then call the
c             field-line computation.
c-----------------------------------------------------------------------
        IF (iters==3.AND.btop_check/="none".AND.
     $      .NOT.apply_blen) THEN
          apply_blen=.true.
          DO ibl=1,nbl
            fllen(ibl)=0._r8
          ENDDO
        ENDIF
        IF (btop_check/="none".AND.apply_blen) THEN
          DO ibl=1,nbl
            CALL vector_assignq_vec(be_n0(ibl),sl2(ibl),nqty=1_i4,
     $                              nstart1=1_i4,nstart2=1_i4)
          ENDDO
          SELECT CASE(btop_check)
          CASE ("beq-trace")
            CALL calc_blen(ref_len,stiff_solver)
          CASE ("passive adv")
            integrand_flag="forward"
            CALL advect_solve(bcsymm,ref_len)
          END SELECT
        ENDIF
c-----------------------------------------------------------------------
c       1) b. update the surface-flux, the FF' and P' flux functions,
c             and adjust the f1 parameter if targeting an Ip-value and
c             the approximate-Newton scheme is not used for Ip.
c-----------------------------------------------------------------------
        IF (dellam>0._r8.AND.ip_weight>0._r8.AND.relerr<1._r8) THEN
          CALL gs_findnl("no f1 update",iters,set_jdir)
        ELSE
          CALL gs_findnl("update f1",iters,set_jdir)
        ENDIF
c-----------------------------------------------------------------------
c       1) c. find the plasma current contained in the domain.  also set
c             the scaled error in plasma current.
c-----------------------------------------------------------------------
        CALL plasma_cur(mu0,pli,curarr,r2arr)
        iperr=ip_weight*(pli-ip_equil)
c-----------------------------------------------------------------------
c       1) d. project -(delstar(psi_tot)+FF'+mu0*R**2*P') and
c             ( (mu0*J_phi/R)*R**2-FF'-mu0*R**2*P' ) onto
c             the basis functions as the residual vector.
c-----------------------------------------------------------------------
        integrand_flag="all terms"
        DO ibl=1,nrbl
          CALL rblock_get_rhs(rb(ibl),rhs(ibl),err_gs,2_i4)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_get_rhs(tb(ibl),rhs(ibl),err_gs,2_i4)
        ENDDO
        CALL edge_network(2_i4,0_i4,poly_degree-1_i4,.TRUE.)
c-----------------------------------------------------------------------
c       form the surface residual which is the difference between
c       the lambda along the surface and the sum of all contributions
c       from interior and exterior currents.  here, sln is
c       used as temporary space.  the result needs to replace the
c       surface-node values in the residual vector (rhs), so send
c       the rhs data into the boundary-condition routines, then add
c       the surface residual.
c
c       [the second component, J_phi/R is like an n=0 scalar,
c       regularity-wise, so nothing needs to be applied along R=0.]
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          crhs(ibl)=0._r8
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'REAL',1_i4)
          CALL cvector_assign_vec(crhs(ibl),rwork2(ibl),'IMAG',1_i4,
     $                            nqty=1_i4)
        ENDDO
     
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4,bcsymm)
        ENDDO

        DO ibl=1,nbl
          CALL vector_assign_cvec(rhs(ibl),crhs(ibl),'REAL',1_i4)
          CALL vector_assign_cvec(sln(ibl),crhs(ibl),'IMAG',1_i4,
     $                            nqty=1_i4)
          CALL vector_add(sln(ibl),rwork2(ibl),v1fac=-1._r8,v2fac=1._r8)
          CALL vector_add(sln(ibl),rwork1(ibl),         ! weighted edge
     $                    v1fac=bweight,v2fac=-bweight) ! residual
          sl2(ibl)=0._r8
          CALL vector_assignq_vec(sl2(ibl),sln(ibl),nqty=1_i4,
     $                            nstart1=1_i4,nstart2=1_i4)
          CALL vector_add(rhs(ibl),sl2(ibl))
        ENDDO
c-----------------------------------------------------------------------
c       now that sl2 is available, project -delstar(psi_tot)
c       onto the basis functions for normalization.
c-----------------------------------------------------------------------
        integrand_flag="delstar only"
        DO ibl=1,nrbl
          CALL rblock_get_rhs(rb(ibl),sl2(ibl),err_gs,2_i4)
          CALL edge_load_arr(sl2(ibl),2_i4,poly_degree-1_i4,seam(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_get_rhs(tb(ibl),sl2(ibl),err_gs,2_i4)
          CALL edge_load_arr(sl2(ibl),2_i4,0_i4,seam(ibl))
        ENDDO
        CALL edge_network(2_i4,0_i4,poly_degree-1_i4,.FALSE.)
        DO ibl=1,nrbl
          CALL edge_unload_arr(sl2(ibl),2_i4,poly_degree-1_i4,seam(ibl))
        ENDDO
        DO ibl=1,nbl
          CALL edge_unload_arr(sl2(ibl),2_i4,0_i4,seam(ibl))
        ENDDO
c-----------------------------------------------------------------------
c       find the relative error.
c-----------------------------------------------------------------------
        CALL iter_2_norm(sl2,er0,1.0_r8,nbl,nrbl,poly_degree,.TRUE.)
        CALL iter_2_norm(rhs,err2,1._r8,nbl,nrbl,poly_degree,.TRUE.)
        relerr=err2/(er0*linffp)
        IF (relerr<1._r8.AND.dellam>0._r8) THEN
          gsuse=1._r8
        ELSE
          gsuse=gscenter
        ENDIF

        IF (ip_equil/=0._r8) err2=SQRT(err2**2+iperr**2)
        IF (iters>0.OR.use_delst_init) THEN
          IF (write_iters.AND.node==0) THEN
            WRITE(*,'(a,i4,2(a,es11.4),a,i4,a,es11.4)') " GS it ",iters,
     $        " ||delst(psi)||=",er0," ||resid||=",err2,
     $        " CG its=",its ," Ip=",pli
            IF (nvfbcoil>0.AND.node==0)
     $        WRITE(*,'(a,es12.5)') "   Vert. FB current=",vfbc_curr
            IF (nrfbcoil>0.AND.node==0)
     $        WRITE(*,'(a,es12.5)') "   Horiz. FB current=",rfbc_curr
          ENDIF
          IF (err2/er0<gsh_tol) converged=.true.
          IF (f_model=="quad_closed".AND.iters<2.AND.
     $        .NOT.use_delst_init) converged=.false.
        ENDIF
c-----------------------------------------------------------------------
c       2) exit loop if convergence criterion is met.  before doing so,
c          save the residual in rwork1.
c-----------------------------------------------------------------------
        IF (converged.OR.iters==maxiters) THEN
          DO ibl=1,nbl
            CALL vector_assignq_vec(rwork1(ibl),rhs(ibl),nqty=1_i4,
     $                              nstart1=1_i4,nstart2=1_i4)
            CALL vector_mult(rwork1(ibl),1._r8/er0)
          ENDDO
          IF (ip_equil/=0._r8.AND.node==0) THEN
            IF (f_model=="quad_closed") THEN
              WRITE(*,'(a,es12.5)') " Final f1_closed=",f1_closed
            ELSE IF (f_model=="cubic_closed") THEN
              WRITE(*,'(a,es12.5)') " Final f_axis=",f_axis
            ENDIF
          ENDIF
          IF (nvfbcoil>0.AND.node==0)
     $      WRITE(*,'(a,es12.5)') " Final Vert. FB current=",vfbc_curr
          IF (nrfbcoil>0.AND.node==0)
     $       WRITE(*,'(a,es12.5)') " Final Horiz. FB current=",rfbc_curr

          IF (iters==maxiters) WRITE(*,*) "no convergence"

          EXIT non_it
        ENDIF
        iters=iters+1
c-----------------------------------------------------------------------
c       3) if dellam>0, evaluate a centered difference of the nonlinear
c          terms about the current value of lambda for use in the linear
c          solves.  note that field-line lengths are not varied.
c-----------------------------------------------------------------------
        IF (dellam>0._r8) THEN
          IF (relerr<1._r8) THEN
            DO ibl=1,nbl
              sln(ibl)=0.5_r8*dellam
              rt2(ibl)=0._r8
              CALL vector_assignq_vec(rt2(ibl),sln(ibl),nqty=1_i4)
              CALL vector_add(sl2(ibl),rt2(ibl))
              CALL pflux_block(sl2(ibl),ibl)
            ENDDO
            CALL gs_findnl("no f1 update",iters,set_jdir)
            DO ibl=1,nbl
              nd_n0(ibl)=rwork3(ibl)      !  F*F' at incremented lambda
              pres_n0(ibl)=pres_eq(ibl)   !  P'   at incremented lambda
              CALL vector_add(sl2(ibl),rt2(ibl),v2fac=-2._r8)
              CALL pflux_block(sl2(ibl),ibl)
            ENDDO
            CALL gs_findnl("no f1 update",iters,set_jdir)
            DO ibl=1,nbl  !  form delta(F*F') and delta(P')
              CALL vector_add(nd_n0(ibl),rwork3(ibl),v2fac=-1._r8)
              CALL vector_add(pres_n0(ibl),pres_eq(ibl),v2fac=-1._r8)
              IF (ibl<=nrbl) THEN
                CALL rblock_qp_update(rb(ibl)%nd_n0,
     $                                rb(ibl)%qnd_n0,rb(ibl))
                CALL rblock_qp_update(rb(ibl)%pres_n0,
     $                                rb(ibl)%qpres_n0,rb(ibl))
              ELSE
                CALL tblock_qp_update(tb(ibl)%nd_n0,
     $                                tb(ibl)%qnd_n0,tb(ibl))
                CALL tblock_qp_update(tb(ibl)%pres_n0,
     $                                tb(ibl)%qpres_n0,tb(ibl))
              ENDIF
              CALL vector_add(sl2(ibl),rt2(ibl)) ! reset original lambda
              CALL pflux_block(sl2(ibl),ibl)
            ENDDO
c-----------------------------------------------------------------------
c           find the matrix elements that couple the F parameter and
c           lambda at the nodes to linearize the Ip convergence.  these
c           matrix elements are labeled A_lam_f and A_f_lam (diperr
c           quantity indices 1 & 2, respectively) in comments.
c-----------------------------------------------------------------------
            IF (ip_equil/=0._r8.AND.ip_weight>0._r8) THEN
              CALL gs_findnl("dFFdf only",iters,set_jdir)
              DO ibl=1,nrbl
                CALL rblock_get_rhs(rb(ibl),diperr(ibl),err_ip,2_i4)
                CALL edge_load_arr(diperr(ibl),2_i4,
     $                             poly_degree-1_i4,seam(ibl))
              ENDDO 
              DO ibl=nrbl+1,nbl
                CALL tblock_get_rhs(tb(ibl),diperr(ibl),err_ip,2_i4)
                CALL edge_load_arr(diperr(ibl),2_i4,0_i4,seam(ibl))
              ENDDO 
              CALL edge_network(2_i4,0_i4,poly_degree-1_i4,.false.)
              DO ibl=1,nbl
                CALL edge_unload_arr(diperr(ibl),2_i4,poly_degree-1_i4,
     $                               seam(ibl))
              ENDDO
            ELSE
              CALL gs_findnl("no f1 update",iters,set_jdir)
            ENDIF
c-----------------------------------------------------------------------
c         just zero-out the arrays for the linearized FF' and P' terms
c         used by the integrand.
c-----------------------------------------------------------------------
          ELSE
            DO ibl=1,nrbl
              rb(ibl)%qnd_n0%qpf=0._r8
              rb(ibl)%qpres_n0%qpf=0._r8
              rb(ibl)%qrwork4%qpf=0._r8
              diperr(ibl)=0._r8
            ENDDO
            DO ibl=nrbl+1,nbl
              tb(ibl)%qnd_n0%qpf=0._r8
              tb(ibl)%qpres_n0%qpf=0._r8
              tb(ibl)%qrwork4%qpf=0._r8
              diperr(ibl)=0._r8
            ENDDO
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       alter the rhs array as part of the matrix partitioning when
c       performing newton iteration with a targeted plasma current.
c-----------------------------------------------------------------------
        IF (dellam>0._r8.AND.ip_equil/=0._r8.AND.relerr<1._r8.AND.
     $      ip_weight>0._r8) THEN
          DO ibl=1,nbl
            rt2(ibl)=0._r8
            CALL vector_assignq_vec(rt2(ibl),diperr(ibl),nqty=1_i4)
            CALL vector_add(rhs(ibl),rt2(ibl),
     $                      v2fac=-iperr/(ip_weight*dipdf))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       4) invert the approximate Jacobian matrix to find the changes in
c          psi/R**2 and in mu0*J_phi/R.  the rhs vector is the residual
c          vector from above, multiplied by -1, so that the solution is
c          the change in lambda and mu0*J_phi/R.
c
c       symmetry conditions are used for vertically unstable equilibria
c       where the solves are sensitive.
c
c       the norm-adjusting linear tolerance avoids progressively
c       restrictive absolute tolerance as the nonlinear resid, hence
c       the rhs for this system, shrinks.  see the team meeting
c       presentation from 4/2/16.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          sl2(ibl)=0._r8
          CALL vector_mult(rhs(ibl),-1._r8)
        ENDDO

        CALL iter_gmr_r2d_mfsolve(eqdot_delstlj,delstlj_mat,delstlj_fac,
     $                            sl2,rhs,2_i4,tollin,linmaxit,
     $                            nimeq_solver,err,its,seed)

        IF (err>tollin) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Delstar: no convergence: ',
     $          its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        tollin=MIN(tol_max,tollin*tol_incr)
c-----------------------------------------------------------------------
c       5) if needed, update the F parameter, then update the
c          poloidal flux function from the solution for psi/R**2
c-----------------------------------------------------------------------
        IF (dellam>0._r8.AND.ip_equil/=0._r8.AND.relerr<1._r8.AND.
     $      ip_weight>0._r8) THEN
          afldl=0._r8
          DO ibl=1,nbl
            CALL vector_assignq_vec(rwork3(ibl),diperr(ibl),1_i4,
     $                              nstart1=1_i4,nstart2=2_i4)
            CALL vector_assignq_vec(sln(ibl),sl2(ibl),1_i4)
            CALL iter_dot(afldl,rwork3(ibl),sln(ibl),seam(ibl),ibl,nrbl,
     $                    poly_degree,.true.)
            IF (nprocs>1) THEN
              CALL mpi_allreduce(afldl,rtmp,1,mpi_nim_real,mpi_sum,
     $             comm_layer,ierror)
              afldl=rtmp
            ENDIF
          ENDDO

          SELECT CASE (f_model)
          CASE("quad_closed")
            f1_closed=f1_closed-(afldl+iperr)/(ip_weight*dipdf)
          CASE("cubic_closed")
            f_axis=f_axis-(afldl+iperr)/(ip_weight*dipdf)
          CASE("linear")
            dfdpsi=dfdpsi-(afldl+iperr)/(ip_weight*dipdf)
          END SELECT
        ENDIF
c-----------------------------------------------------------------------
c       update the poloidal flux function with the solution for
c       psi/R**2, and store mu*J_phi/R.  before the following
c       block, sl2 holds the change in (lam,mu0*J_phi/R).
c       afterwards, it holds their updated values.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL vector_assignq_vec(rt2(ibl),rwork2(ibl),nqty=1_i4,
     $                            nstart1=1_i4,nstart2=1_i4)
          CALL vector_assignq_vec(rt2(ibl),ja_eq(ibl),nqty=1_i4,
     $                            nstart1=2_i4,nstart2=3_i4)
          CALL vector_add(sl2(ibl),rt2(ibl),v1fac=gsuse,v2fac=1._r8)
          CALL vector_assignq_vec(rwork2(ibl),sl2(ibl),nqty=1_i4,
     $                            nstart1=1_i4,nstart2=1_i4)
          CALL vector_assignq_vec(ja_eq(ibl),sl2(ibl),nqty=1_i4,
     $                            nstart1=3_i4,nstart2=2_i4)
          CALL pflux_block(sl2(ibl),ibl)
        ENDDO

        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%ja_eq,rb(ibl)%qja_eq,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%ja_eq,tb(ibl)%qja_eq,tb(ibl))
        ENDDO
c-----------------------------------------------------------------------
c-PRE   the following 4 lines result from merging the new feedback
c       scheme.  how it would fit with the approximate newton iteration
c       needs more thought.
c-----------------------------------------------------------------------
c       CALL nimeq_free_eval(rwork1,1_i4,poly_degree,1._r8,
c    $                       psio_max,psio_min,iters>=3_i4)
c       oflux=psio_max-psio_min
c       CALL pflux_range(btop_check,set_jdir)
      ENDDO non_it
c-----------------------------------------------------------------------
c     deallocate arrays used for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
        CALL vector_type_dealloc(sl2(ibl))
        CALL vector_type_dealloc(rt2(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
      DEALLOCATE(sl2,rt2)
c-----------------------------------------------------------------------
c     create the expansions of the equilibrium fields that are used
c     in nimrod.  also write output for plotting.
c-----------------------------------------------------------------------
      CALL equil_project
      CALL equil_plot(first_data,ndcon)
      nmodes=nmodes_save
c-----------------------------------------------------------------------
c     transfer equilibrium fields to n=0 if requested.
c-----------------------------------------------------------------------
      IF (nimeq_tx.AND.nonlinear) THEN
        CALL eq_swap(nmodes,keff,geom,eq_flow,ndeq_notx,rbph_notx)
        IF (.NOT.coil_tx.AND.geom=='tor'.AND.ncoil_tot>0)
     $    CALL replace_coil(nmodes,keff,ncoil_tot,allcoil_r,
     $                      allcoil_z,allcoil_i)
      ENDIF
c-----------------------------------------------------------------------
c     write the distribution of internal current density if this is
c     a free-boundary computation.
c-----------------------------------------------------------------------
      IF (gs_type=='free') CALL nimeq_free_iwrite
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gsfree

c-----------------------------------------------------------------------
c     subprogram 3. delstar
c     calculates poloidal flux from toroidal current
c-----------------------------------------------------------------------
      SUBROUTINE delstar(first_data,datflag,phii,ndcon)
      USE local
      USE fields
      USE input
      USE input_eq
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimeq_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      USE boundary
      USE regularity
      USE nimeq_mod
      USE mpi_nim
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      CHARACTER(*), INTENT(IN) :: datflag
      REAL(r8), INTENT(IN) :: phii
      INTEGER(i4) :: ibl,jphi,its,nmodes_save
      INTEGER(i4) :: ndcon,ibe,ibloc,ierror
      LOGICAL :: first_point,first_r0
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
      TYPE(vector_type), DIMENSION(:), POINTER :: vtmp1,vtmp2
c-----------------------------------------------------------------------
c     determine the inhomogeneous part of the poloidal flux along the
c     border and save it (divided by 2*pi*R**2 in toroidal geometry
c     or by per_length in linear) in rwork1.
c-----------------------------------------------------------------------
      IF (datflag/='surf set') THEN
        DO ibl=1,nbl
          rwork1(ibl)=0._r8
        ENDDO
        CALL gs_surface_flux(psio_min,psio_max,.false.)
      ENDIF
c-----------------------------------------------------------------------
c     find the integral -[alpha*u0*jphi]
c     -[gradalpha*grad(Psi/r^2)*r^2)] along the surface
c
c     rhs structures are now allocated and deallocated on the fly to
c     save memory. Uses a dummie crhs for use with bc routines
c-----------------------------------------------------------------------
      ALLOCATE(vtmp1(nrbl),vtmp2(nrbl))

      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,1_i4)       
        CALL vector_type_alloc(vtmp1(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4)
        CALL vector_type_alloc(vtmp2(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4)
        CALL rblock_get_rhs(rb(ibl),rhs(ibl),delstar_rhs,1_i4)
        crhs(ibl)=0._r8
      END DO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         1_i4,1_i4) 
        CALL tblock_get_rhs(tb(ibl),rhs(ibl),delstar_rhs,1_i4)
        crhs(ibl)=0._r8
      END DO   
c-----------------------------------------------------------------------
c     eliminate the element interiors (static condensation) before
c     calling the iterative solver.     
c-----------------------------------------------------------------------
      IF (poly_degree>1) THEN
        DO ibl=1,nrbl
          vtmp1(ibl)=rhs(ibl)
        ENDDO
        CALL matelim_presolve(delstar_mat,vtmp1,vtmp2,1_i4)
        DO ibl=1,nrbl
          CALL cvector_assign_vec(crhs(ibl),vtmp2(ibl),'real',
     $                            1_i4,1_i4)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'real',1_i4)
        ENDDO
      ELSE
        DO ibl=1,nbl
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'real',1_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     network block seams and apply boundary conditions
c-----------------------------------------------------------------------
      CALL edge_network(1_i4,1_i4,poly_degree-1_i4,.true.)

      nmodes_save=nmodes
      nmodes=1
      DO ibl=1,SIZE(r0block_list)
         ibe=r0block_list(ibl)
         CALL regular_vec(crhs(ibe),seam(ibe),'all',1_i4,
     $                    nmodes,nindex)
      ENDDO
      nmodes=nmodes_save

      DO ibl=1,SIZE(exblock_list)
         ibe=exblock_list(ibl)
         CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4)
      ENDDO
c-----------------------------------------------------------------------
c     transfer data for the system solve.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
         CALL vector_assign_cvec(rhs(ibl),crhs(ibl),'real',1_i4,1_i4)
         sln(ibl)=0        
      ENDDO 
c-----------------------------------------------------------------------
c     invert the delstar matrix to find psi
c-----------------------------------------------------------------------
      CALL iter_cg_2d_solve(delstar_mat,delstar_fac,sln,rhs,1_i4,tol,
     $                      linmaxit,nimeq_solver,err,its,seed)
      IF (err>tol) THEN
         WRITE(msg,'(a,i4,a,es10.3,3a)') 'Delstar: no convergence: ',
     $         its,' its ',err,' err ',seed,' seed'
         CALL nim_stop(msg)
      ENDIF
c-----------------------------------------------------------------------
c     complete the solution at element interiors.
c-----------------------------------------------------------------------
      IF (poly_degree>1) THEN
        CALL matelim_real_postsolve(delstar_mat,sln,rhs,1_i4)
        DO ibl=1,nrbl
          sln(ibl)%arri=rhs(ibl)%arri
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     copy the solution to rwork2 to help initialize a subsequent GS
c     computation, then combine it with the wall flux and set pflux.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        rwork2(ibl)=sln(ibl)
        CALL vector_add(sln(ibl),rwork1(ibl))
        CALL pflux_block(sln(ibl),ibl)
        IF (geom=='tor') THEN
          CALL vector_mult(pflux(ibl),twopi)
        ELSE
          CALL vector_mult(pflux(ibl),per_length)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     the following write scheme assumes that all processors have
c     access to the same disk space.
c-----------------------------------------------------------------------
      IF (first_data) CALL par_cont_init("pflux.bin",1_i4,ndcon)
      DO ibl=1,nrbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,"pflux.bin","OLD","APPEND",32_i4)
          CALL contour_laq2_write(rb(ibloc)%pflux,ndcon)
          CALL close_bin(con_unit,"pflux.bin")
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
      DO ibl=nrbl_total+1,nbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,"pflux.bin","OLD","APPEND",32_i4)
          CALL contour_tl2_write(tb(ibloc)%pflux)
          CALL close_bin(con_unit,"pflux.bin")
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
c-----------------------------------------------------------------------
c     deallocate rhs storage
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
        CALL vector_type_dealloc(crhs(ibl))
        IF (ibl<=nrbl) THEN
          CALL vector_type_dealloc(vtmp1(ibl))
          CALL vector_type_dealloc(vtmp2(ibl))
        ENDIF
      ENDDO
      DEALLOCATE(vtmp1,vtmp2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstar

c-----------------------------------------------------------------------
c     subprogram 4. gs_coil_read
c     this subroutine reads boundary flux values and locations from a
c     binary file.
c-----------------------------------------------------------------------
      SUBROUTINE gs_coil_read(ncurfil,curfil)
      USE local
      USE io
      USE input_eq
      USE pardata
      USE mpi_nim
      IMPLICIT NONE

      INTEGER(i4), INTENT(OUT) :: ncurfil
      REAL(r8), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: curfil

      INTEGER(i4), PARAMETER :: nrecmx=2500000
      INTEGER(i4) :: ird,ntmp,ierror
      INTEGER :: read_stat
      REAL(r8), DIMENSION(3) :: data
c-----------------------------------------------------------------------
c     first determine the number of records in the file.
c-----------------------------------------------------------------------
      IF (node==0) THEN
        IF (TRIM(curfil_file)=="none") THEN
          ncurfil=0
        ELSE
          CALL open_bin(binary_unit,TRIM(curfil_file),'OLD',
     $                  'REWIND',64_i4)
          DO ncurfil=1,nrecmx+1
            READ(binary_unit,IOSTAT=read_stat) data
            IF (read_stat/=0) EXIT
          ENDDO
          CALL close_bin(binary_unit,TRIM(curfil_file))
          IF (ncurfil>nrecmx) THEN
            CALL nim_stop('Too many records; increase nrecmx.')
          ELSE IF (ncurfil==1) THEN
            CALL nim_stop
     $        ('File '//TRIM(curfil_file)//' is empty?? (ncurfil=0)')
          ENDIF
          ncurfil=ncurfil-1
        ENDIF
      ENDIF
      IF (nprocs>1) CALL mpi_bcast(ncurfil,1,mpi_nim_int,0,
     $                             mpi_comm_world,ierror)
c-----------------------------------------------------------------------
c     read the data into the curfil array.
c-----------------------------------------------------------------------
      ALLOCATE(curfil(3,ncurfil))
      IF (node==0.AND.ncurfil>0) THEN
        CALL open_bin(binary_unit,TRIM(curfil_file),'OLD',
     $                'REWIND',64_i4)
        DO ird=1,ncurfil
          READ(binary_unit) curfil(:,ird)
        ENDDO
        CALL close_bin(binary_unit,TRIM(curfil_file))
      ENDIF
      IF (nprocs>1.AND.ncurfil>0)
     $  CALL mpi_bcast(curfil(1,1),3*ncurfil,mpi_nim_real,
     $                 0,mpi_comm_world,ierror)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gs_coil_read

c-----------------------------------------------------------------------
c     subprogram 5. gs_surface_flux
c     use the normal magnetic field along the surface of the domain
c     to determine the inhomogenous part of the flux boundary
c     conditions.
c
c     modified to treat top and bottom of up/down symmetric pieces
c     symmetrically.
c-----------------------------------------------------------------------
      SUBROUTINE gs_surface_flux(psi_min,psi_max,from_gs)
      USE local
      USE fields
      USE input
      USE input_eq
      USE computation_pointers
      USE seam_storage_mod
      USE physdat
      USE nimeq_mod
      USE pardata
      USE mpi_nim
      IMPLICIT NONE

      REAL(r8), INTENT(OUT) :: psi_min,psi_max
      LOGICAL, INTENT(IN) :: from_gs

      INTEGER(i4) :: ix,iy,iv0,nv0,ip0,np0,iv,nv,ig,ib,ip,np,ivp,
     $               iv0p,np0p,ip0p,dix,diy,is,ibloc,itmp,ierror
      INTEGER(i4), SAVE :: ng
      REAL(r8) :: x,y,dx,dy,psis,per_psi,psir0,per_fac,psi_coil,psicst,
     $            rtmp
      REAL(r8), DIMENSION(2) :: bsurf,dsurf,rzsurf
      REAL(r8), DIMENSION(:), ALLOCATABLE, SAVE :: xnode,xg1d,wg1d
      REAL(r8), EXTERNAL :: aph_eval
      LOGICAL :: connected,r0side,add_coil

      REAL(r8), PARAMETER :: postol=1.e-8
      INTEGER(i4) :: iv00,iv0dir,iv0st,iv0en,iss,ivseg
      INTEGER(i4), DIMENSION(nprocs) :: ipr0,iptmp
c-----------------------------------------------------------------------
c     find the 1D node distribution along element edges, and determine
c     weights and locations for Gaussian quadrature.
c-----------------------------------------------------------------------
      IF (.NOT.ALLOCATED(xnode)) THEN
        ALLOCATE(xnode(0:poly_degree))
        CALL poly_nodes(poly_degree,xnode)
        ng=MAX(2_i4,poly_degree/2_i4+1_i4)
        ALLOCATE(xg1d(ng))
        ALLOCATE(wg1d(ng))
        CALL gauleg(0._r8,1._r8,xg1d,wg1d,ng)
      ENDIF
c-----------------------------------------------------------------------
c     initialize the flux for all images of the last location in the
c     seam0 structure.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        per_fac=1._r8/twopi
      ELSE
        per_fac=1._r8/per_length
      ENDIF
      psi_max=-HUGE(psi_max)
      psi_min= HUGE(psi_max)
      psir0=0._r8
      psis=0._r8
      nv0=seam0%nvert
      add_coil=from_gs.AND.geom=='tor'.AND.ncoil_tot>0
c-----------------------------------------------------------------------
c     for the top and bottom pieces, start along the inboard surface
c     at the symmetry plane, and integrate in the direction away from
c     the symmetry plane.
c-----------------------------------------------------------------------
      iv0st=0
      IF (symm_region=="neither") THEN
        iv0dir=1
        iv0en=nv0
      ELSE
        DO iv0=1,nv0
          ib=seam0%vertex(iv0)%ptr(1,1)
          iv=seam0%vertex(iv0)%ptr(2,1)
          IF (block2proc(ib)==node) THEN
            ibloc=global2local(ib)
            ix=seam(ibloc)%vertex(iv)%intxy(1)
            iy=seam(ibloc)%vertex(iv)%intxy(2)
            IF (ibloc<=nrbl) THEN
              rzsurf=rb(ibloc)%rz%fs(:,ix,iy)
            ELSE
              rzsurf=(/tb(ibloc)%tgeom%xs(ix),tb(ibloc)%tgeom%ys(ix)/)
            ENDIF
            IF (ABS(rzsurf(1)-xmin)<=postol.AND.
     $          ABS(rzsurf(2))<=postol) THEN
              iv0st=iv0
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(iv0st,itmp,1,mpi_nim_int,mpi_max,
     $         mpi_comm_world,ierror)
          iv0st=itmp
        ENDIF
        IF (iv0st==0) CALL nim_stop("Start location not found.")

        IF (symm_region=="top") THEN
          iv0dir=-1
          iv0en=iv0st+1
          IF (iv0en>nv0) iv0en=1
        ELSE
          iv0dir=1
          iv0en=iv0st-1
          IF (iv0en<1) iv0en=nv0
        ENDIF
      ENDIF

      np0=SIZE(seam0%vertex(iv0en)%ptr,2)
      DO ip0=1,np0
        ib=seam0%vertex(iv0en)%ptr(1,ip0)
        IF (block2proc(ib)==node) THEN
          ibloc=global2local(ib)
          iv=seam0%vertex(iv0en)%ptr(2,ip0)
          ix=seam(ibloc)%vertex(iv)%intxy(1)
          iy=seam(ibloc)%vertex(iv)%intxy(2)
          rwork1(ibloc)%arr(1,ix,iy)=psis
        ENDIF
      ENDDO
      iv0p=iv0en
      np0p=np0
      psicst=0._r8
c-----------------------------------------------------------------------
c     loop over seam0 to locate adjacent mesh vertices along the border
c     of the domain.  the block seam structures are defined such that
c     the seam segment lies between the previous point and the next
c     point.
c
c     the first part identifies the current block for places along
c     seam0 where two grid blocks join.
c
c     in parallel, this works 1 node at a time, so it is a possible
c     bottleneck, but this loop is expected to be a tiny fraction of
c     the cpu time.
c-----------------------------------------------------------------------
      iv0=iv0en
      DO iv00=1,nv0
        iv0=iv0+iv0dir
        IF (iv0<1) iv0=nv0
        IF (iv0>nv0) iv0=1
        np0=SIZE(seam0%vertex(iv0)%ptr,2)
        CALL gs_surface_search(connected,ib,iv)
c-----------------------------------------------------------------------
c       if there is a break in the seam0 vertices, it is from a
c       periodic gap.  increment psis by the specified amount of flux.
c-----------------------------------------------------------------------
        IF (.NOT.connected) THEN
          IF (symm_region=="neither") THEN
            psis=psis+per_fac*per_flux
            psi_max=MAX(psi_max,psis)
            psi_min=MIN(psi_min,psis)
            per_fac=-1._r8*per_fac
          ELSE
            CALL nim_stop("Error finding next vertex.")
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       determine the block-indices and their increments based on
c       seam-segment data.
c-----------------------------------------------------------------------
        IF (ib>0) THEN     !    block is on this processor
          IF (iv0dir>0) THEN
            ix=seam(ib)%segment(iv)%intxyp(1)
            iy=seam(ib)%segment(iv)%intxyp(2)
            dix=seam(ib)%segment(iv)%intxyn(1)-ix
            diy=seam(ib)%segment(iv)%intxyn(2)-iy
          ELSE
            ivseg=iv+1
            IF (ivseg>seam(ib)%nvert) ivseg=1
            ix=seam(ib)%segment(ivseg)%intxyp(1)
            iy=seam(ib)%segment(ivseg)%intxyp(2)
            dix=seam(ib)%segment(ivseg)%intxyn(1)-ix
            diy=seam(ib)%segment(ivseg)%intxyn(2)-iy
          ENDIF
c-----------------------------------------------------------------------
c         loop over nodes along the side of the element and quadrature
c         points for each increment.
c
c         the surface-b integration is used for fixed-equilibrium
c         computations.  for both fixed and free, subtract R*A_phi
c         contributions from external coils directly.
c-----------------------------------------------------------------------
          IF (iv0dir>0) THEN
            is=1
          ELSE
            is=poly_degree-1
          ENDIF

          DO iss=1,poly_degree
            dx=dix*(xnode(is)-xnode(is-iv0dir))
            dy=diy*(xnode(is)-xnode(is-iv0dir))
            IF (gs_type=="fixed".OR..NOT.from_gs) THEN
              DO ig=1,ng
                x=ix+xg1d(ig)*dx+dix*xnode(is-iv0dir)
                y=iy+xg1d(ig)*dy+diy*xnode(is-iv0dir)
                CALL gs_surface_b(bsurf,ib,x,y,.false.)
                CALL gs_surface_ds(dsurf,ib,x,y,dx,dy)
                psis=psis-SUM(bsurf*dsurf)*wg1d(ig)
              ENDDO
            ENDIF
            IF (add_coil) THEN
              IF (ib<=nrbl) THEN
                x=ix+dx+dix*xnode(is-iv0dir)
                y=iy+dy+diy*xnode(is-iv0dir)
                CALL lagr_quad_eval(rb(ib)%rz,x,y,0_i4)
                rzsurf=rb(ib)%rz%f(:)
              ELSE
                rzsurf=(/tb(ib)%tgeom%xs(ix),tb(ib)%tgeom%ys(ix)/)
              ENDIF
              psi_coil=psicst-aph_eval(rzsurf(1),rzsurf(2),ncoil_tot,
     $                                allcoil_r,allcoil_z,allcoil_i,mu0)
            ELSE
              psi_coil=0._r8
            ENDIF

            psi_max=MAX(psi_max,psis+psi_coil)
            psi_min=MIN(psi_min,psis+psi_coil)
            IF (iss<poly_degree) THEN
              IF (dix>0) THEN
                rwork1(ib)%arrh(1,is,ix+dix,iy)=psis+psi_coil
              ELSE IF (diy>0) THEN
                rwork1(ib)%arrv(1,is,ix,iy+diy)=psis+psi_coil
              ELSE IF (dix<0) THEN
                rwork1(ib)%arrh(1,poly_degree-is,ix,iy)=psis+psi_coil
              ELSE
                rwork1(ib)%arrv(1,poly_degree-is,ix,iy)=psis+psi_coil
              ENDIF
            ENDIF
            is=is+iv0dir
          ENDDO
        ELSE
          psis=0._r8
          psi_coil=0._r8
        ENDIF
        IF (nprocs>1) THEN
          CALL mpi_allreduce(psis,rtmp,1,mpi_nim_real,mpi_sum,
     $         mpi_comm_world,ierror)
          psis=rtmp
          CALL mpi_allreduce(psi_coil,rtmp,1,mpi_nim_real,mpi_sum,
     $         mpi_comm_world,ierror)
          psi_coil=rtmp
        ENDIF
        DO ip0=1,np0
          ib=seam0%vertex(iv0)%ptr(1,ip0)
          IF (block2proc(ib)==node) THEN
            ibloc=global2local(ib)
            iv=seam0%vertex(iv0)%ptr(2,ip0)
            ix=seam(ibloc)%vertex(iv)%intxy(1)
            iy=seam(ibloc)%vertex(iv)%intxy(2)
            rwork1(ibloc)%arr(1,ix,iy)=psis+psi_coil
          ENDIF
        ENDDO
        iv0p=iv0
        np0p=np0
      ENDDO
      IF (nprocs>1) THEN
        CALL mpi_allreduce(psi_max,rtmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        psi_max=rtmp
        CALL mpi_allreduce(psi_min,rtmp,1,mpi_nim_real,mpi_min,
     $       mpi_comm_world,ierror)
        psi_min=rtmp
      ENDIF
c-----------------------------------------------------------------------
c     if the geometry is toroidal, two more adjustments are needed.
c     if R=0 is part of the domain, shift the flux values to make it
c     zero on axis.  then, divide by R**2.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        ipr0=0
        DO iv0=1,nv0
          ib=seam0%vertex(iv0)%ptr(1,1)
          iv=seam0%vertex(iv0)%ptr(2,1)
          IF (block2proc(ib)==node) THEN
            ibloc=global2local(ib)
            IF (seam(ibloc)%r0point(iv)) THEN
              ipr0(node+1)=1
              ix=seam(ibloc)%vertex(iv)%intxy(1)
              iy=seam(ibloc)%vertex(iv)%intxy(2)
              psir0=rwork1(ibloc)%arr(1,ix,iy)
            ENDIF
          ENDIF
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(ipr0,iptmp,nprocs,mpi_nim_int,mpi_max,
     $         mpi_comm_world,ierror)
          ipr0=iptmp
          IF (MAXVAL(ipr0)>0) THEN
            DO itmp=1,nprocs
              IF (ipr0(itmp)>0) EXIT
            ENDDO
            CALL mpi_bcast(psir0,1,mpi_nim_real,itmp-1_i4,
     $                     mpi_comm_world,ierror)
          ENDIF
        ENDIF
        psi_min=psi_min-psir0
        psi_max=psi_max-psir0

        iv0p=nv0
        np0p=SIZE(seam0%vertex(nv0)%ptr,2)
        DO iv0=1,nv0
          np0=SIZE(seam0%vertex(iv0)%ptr,2)
          CALL gs_surface_search(connected,ib,iv)
          IF (ib>0) THEN
            ix=seam(ib)%segment(iv)%intxyp(1)
            iy=seam(ib)%segment(iv)%intxyp(2)
            dix=seam(ib)%segment(iv)%intxyn(1)-ix
            diy=seam(ib)%segment(iv)%intxyn(2)-iy
            ivp=iv-1
            IF (ivp==0) ivp=seam(ib)%nvert
            r0side=seam(ib)%r0point(ivp).AND.seam(ib)%r0point(iv)
            DO is=1,poly_degree-1
              IF (dix>0) THEN
                IF (r0side) THEN
                  rwork1(ib)%arrh(1,is,ix+dix,iy)=0._r8
                ELSE
                  rwork1(ib)%arrh(1,is,ix+dix,iy)=
     $              (rwork1(ib)%arrh(1,is,ix+dix,iy)-psir0)/
     $               rb(ib)%rz%fsh(1,is,ix+dix,iy)**2
                ENDIF
              ELSE IF (diy>0) THEN
                IF (r0side) THEN
                  rwork1(ib)%arrv(1,is,ix,iy+diy)=0._r8
                ELSE
                  rwork1(ib)%arrv(1,is,ix,iy+diy)=
     $              (rwork1(ib)%arrv(1,is,ix,iy+diy)-psir0)/
     $               rb(ib)%rz%fsv(1,is,ix,iy+diy)**2
                ENDIF
              ELSE IF (dix<0) THEN
                IF (r0side) THEN
                  rwork1(ib)%arrh(1,poly_degree-is,ix,iy)=0._r8
                ELSE
                  rwork1(ib)%arrh(1,poly_degree-is,ix,iy)=
     $              (rwork1(ib)%arrh(1,poly_degree-is,ix,iy)-psir0)/
     $               rb(ib)%rz%fsh(1,poly_degree-is,ix,iy)**2
                ENDIF
              ELSE
                IF (r0side) THEN
                  rwork1(ib)%arrv(1,poly_degree-is,ix,iy)=0._r8
                ELSE
                  rwork1(ib)%arrv(1,poly_degree-is,ix,iy)=
     $             (rwork1(ib)%arrv(1,poly_degree-is,ix,iy)-psir0)/
     $              rb(ib)%rz%fsv(1,poly_degree-is,ix,iy)**2
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          DO ip0=1,np0
            ib=seam0%vertex(iv0)%ptr(1,ip0)
            IF (block2proc(ib)==node) THEN
              ibloc=global2local(ib)
              iv=seam0%vertex(iv0)%ptr(2,ip0)
              ix=seam(ibloc)%vertex(iv)%intxy(1)
              iy=seam(ibloc)%vertex(iv)%intxy(2)
              IF (seam(ibloc)%r0point(iv)) THEN
                IF (seam(ibloc)%expoint(iv)) THEN
                  CALL gs_surface_b(bsurf,ibloc,REAL(ix,r8),REAL(iy,r8),
     $                              add_coil)
                  rwork1(ibloc)%arr(1,ix,iy)=0.5_r8*bsurf(2)
                ELSE
                  rwork1(ibloc)%arr(1,ix,iy)=0._r8
                ENDIF
              ELSE
                IF (ibloc<=nrbl) THEN
                  rwork1(ibloc)%arr(1,ix,iy)=
     $              (rwork1(ibloc)%arr(1,ix,iy)-psir0)/
     $               rb(ibloc)%rz%fs(1,ix,iy)**2
                ELSE
                  rwork1(ibloc)%arr(1,ix,iy)=
     $              (rwork1(ibloc)%arr(1,ix,iy)-psir0)/
     $               tb(ibloc)%tgeom%xs(ix)**2
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          iv0p=iv0
          np0p=np0
        ENDDO
      ENDIF

      RETURN
c-----------------------------------------------------------------------
c     internal subroutines for gs_surface_flux:
c-----------------------------------------------------------------------
      CONTAINS

c-----------------------------------------------------------------------
c       gs_surface_search determines which block and vertex continues
c       from the previous point in seam0.  for parallel computation,
c       set ibl<0 to indicate that this processor does not own the
c       block.
c-----------------------------------------------------------------------
        SUBROUTINE gs_surface_search(connect,ibl,ivert)

        LOGICAL, INTENT(OUT) :: connect
        INTEGER(i4), INTENT(OUT) :: ibl,ivert

        INTEGER(i4) :: ivp,ibt,ibgl,ibgp,icon,itmp

        connect=.false.
        icon=0
        ibl=-1
        p0_search: DO ip0=1,np0
          ibgl=seam0%vertex(iv0)%ptr(1,ip0)
          ivert=seam0%vertex(iv0)%ptr(2,ip0)
          IF (block2proc(ibgl)==node) THEN
            ibt=global2local(ibgl)
            ivp=ivert-iv0dir
            IF (ivp==0) ivp=seam(ibt)%nvert
            IF (ivp>seam(ibt)%nvert) ivp=1
            DO ip0p=1,np0p
              ibgp=seam0%vertex(iv0p)%ptr(1,ip0p)
              IF (block2proc(ibgp)==node.AND.
     $            global2local(ibgp)==ibt.AND.
     $            seam0%vertex(iv0p)%ptr(2,ip0p)==ivp) THEN
                icon=1
                ibl=ibt
                EXIT p0_search
              ENDIF
            ENDDO
          ENDIF
        ENDDO p0_search

        IF (nprocs>1) THEN
          CALL mpi_allreduce(icon,itmp,1,mpi_nim_int,mpi_max,
     $         mpi_comm_world,ierror)
          icon=itmp
        ENDIF
        IF (icon>0) connect=.true.
c-----------------------------------------------------------------------
c       if the previous point is not found, assume that it is due to a
c       periodic break.  set the block to the block of the current
c       seam0 connection, and set the vertex to the one previous to
c       the current seam0 connection.  leave connect false to signify
c       the break.
c-----------------------------------------------------------------------
        IF (.NOT.connect) THEN
          ibgl=seam0%vertex(iv0)%ptr(1,1)
          ivert=seam0%vertex(iv0)%ptr(2,1)
          IF (block2proc(ibgl)==node) THEN
            ibl=global2local(ibgl)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       end of internal subroutine.
c-----------------------------------------------------------------------
        RETURN
        END SUBROUTINE gs_surface_search

c-----------------------------------------------------------------------
c       gs_surface_b finds the magnetic field at the position specified
c       in the calling parameters.
c-----------------------------------------------------------------------
        SUBROUTINE gs_surface_b(bb,ibl,xpos,ypos,incl_coil)
        USE physdat
        USE global

        REAL(r8), DIMENSION(2), INTENT(OUT) :: bb
        REAL(r8), INTENT(IN) :: xpos,ypos
        INTEGER(i4), INTENT(IN) :: ibl
        LOGICAL, INTENT(IN) :: incl_coil

        REAL(r8), DIMENSION(2) :: xxyy
        INTEGER(i4) :: im,icoil
c-----------------------------------------------------------------------
c       the evaluation routine depends on the block type.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          CALL lagr_quad_eval(rb(ibl)%be_eq,xpos,ypos,0_i4)
          bb=rb(ibl)%be_eq%f(1:2)
          IF (nonlinear) THEN
            CALL lagr_quad_eval(rb(ibl)%be,xpos,ypos,0_i4)
            DO im=1,SIZE(keff)
              IF (keff(im)==0.) THEN
                bb=bb+rb(ibl)%be%f(1:2,im)
                EXIT
              ENDIF
            ENDDO
          ENDIF
c-PRE
        ELSE   !   triangles
        ENDIF
c-----------------------------------------------------------------------
c       add magnetic field from external coils if this is a GS
c       solve in toroidal geometry.  adding it to B here is only
c       used when resetting psi for regularity.
c-----------------------------------------------------------------------
        IF (ncoil_tot>0.AND.geom=='tor'.AND.incl_coil.AND.ibl<=nrbl)
     $    THEN
          CALL lagr_quad_eval(rb(ibl)%rz,xpos,ypos,1_i4)
          xxyy=rb(ibl)%rz%f
          IF (xxyy(1)<1.e-8*MAXVAL(ABS(rb(ibl)%rz%fs))) THEN
            DO icoil=1,ncoil_tot
              bb(2)=bb(2)-
     $          0.5_r8*mu0*allcoil_i(icoil)*allcoil_r(icoil)**2/
     $          (allcoil_r(icoil)**2
     $          +(xxyy(2)-allcoil_z(icoil))**2)**1.5_r8
            ENDDO
          ELSE
            CALL nim_stop("Gs_surface_b: coil selection error.")
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       end of internal subroutine.
c-----------------------------------------------------------------------
        RETURN
        END SUBROUTINE gs_surface_b

c-----------------------------------------------------------------------
c       gs_surface_ds finds outward-point differential surface vector
c       at the position specified in the calling parameters.
c-----------------------------------------------------------------------
        SUBROUTINE gs_surface_ds(ds,ibl,xpos,ypos,dxpos,dypos)

        REAL(r8), DIMENSION(2), INTENT(OUT) :: ds
        REAL(r8), INTENT(IN) :: xpos,ypos,dxpos,dypos
        INTEGER(i4), INTENT(IN) :: ibl
c-----------------------------------------------------------------------
c       the evaluation routine depends on the block type.
c       the outward pointing differential surface is proportional
c       to dL x phi-hat, where dL is along the surface (dxpos,dypos),
c       and phi-hat is the unit direction vector perpendicular to the
c       plane.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          CALL lagr_quad_eval(rb(ibl)%rz,xpos,ypos,1_i4)
          ds(1)= rb(ibl)%rz%fx(2)*dxpos+rb(ibl)%rz%fy(2)*dypos
          ds(2)=-rb(ibl)%rz%fx(1)*dxpos-rb(ibl)%rz%fy(1)*dypos
          IF (geom=='tor') ds=ds*rb(ibl)%rz%f(1)
c-PRE
        ELSE   !   triangles
        ENDIF
c-----------------------------------------------------------------------
c       end of internal subroutine.
c-----------------------------------------------------------------------
        RETURN
        END SUBROUTINE gs_surface_ds

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE gs_surface_flux

c-----------------------------------------------------------------------
c     subprogram 6. gs_work_init
c     this collects the initialization operations that are common to
c     to gssolve and gsfree.
c-----------------------------------------------------------------------
      SUBROUTINE gs_work_init(refl)
      USE local
      USE fields
      USE pardata
      USE mpi_nim
      USE input
      USE global
      IMPLICIT NONE

      REAL(r8), INTENT(OUT) :: refl

      REAL(r8) :: rtmp
      INTEGER(i4) :: ibl,ix,ierror
c-----------------------------------------------------------------------
c     most of the initialization is for the data used in field-line
c     tracing (temporarily stored in be_n0).
c-----------------------------------------------------------------------
      refl=0._r8
      DO ibl=1,nbl
        rwork1(ibl)=0._r8
        rwork2(ibl)=0._r8
        rwork3(ibl)=0._r8
        pflux(ibl)=0._r8
        be_n0(ibl)=1._r8   !   non-zero initial vals for trace except R:
        IF (ibl<=nrbl.AND.geom=="tor") THEN
          be_n0(ibl)%arr(2,:,:)=MAX(smallnum,rb(ibl)%rz%fs(1,:,:))
          IF (poly_degree>1) THEN
            be_n0(ibl)%arrh(2,:,:,:)=
     $        MAX(smallnum,rb(ibl)%rz%fsh(1,:,:,:))
            be_n0(ibl)%arrv(2,:,:,:)=
     $        MAX(smallnum,rb(ibl)%rz%fsv(1,:,:,:))
            be_n0(ibl)%arri(2,:,:,:)=
     $        MAX(smallnum,rb(ibl)%rz%fsi(1,:,:,:))
          ENDIF
        ELSE IF (geom=="tor") THEN
          DO ix=0,tb(ibl)%mvert
            be_n0(ibl)%arr(2,ix,0_i4)=MAX(smallnum,tb(ibl)%tgeom%xs(ix))
          ENDDO
        ENDIF
        IF (ibl<=nrbl) THEN
          be_n0(ibl)%arr(3:4,:,:)=rb(ibl)%rz%fs(1:2,:,:)
          IF (poly_degree>1) THEN
            be_n0(ibl)%arrh(3:4,:,:,:)=rb(ibl)%rz%fsh(1:2,:,:,:)
            be_n0(ibl)%arrv(3:4,:,:,:)=rb(ibl)%rz%fsv(1:2,:,:,:)
            be_n0(ibl)%arri(3:4,:,:,:)=rb(ibl)%rz%fsi(1:2,:,:,:)
          ENDIF
        ELSE
          DO ix=0,tb(ibl)%mvert
            be_n0(ibl)%arr(3,ix,0_i4)=tb(ibl)%tgeom%xs(ix)
            be_n0(ibl)%arr(4,ix,0_i4)=tb(ibl)%tgeom%ys(ix)
          ENDDO
        ENDIF
        refl=MAX(refl,MAXVAL(be_n0(ibl)%arr(3,:,:)))
      ENDDO 
      IF (nprocs>1) THEN
        CALL mpi_allreduce(refl,rtmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        refl=rtmp
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gs_work_init

c-----------------------------------------------------------------------
c     subprogram 7. gs_findnl
c     evaluates the nonlinear terms of the GS equation using stored
c     values for the poloidal-flux function and the field-line length
c     distributions.
c-----------------------------------------------------------------------
      SUBROUTINE gs_findnl(f_instruct,iters,setj)
      USE local
      USE nimeq_mod
      USE input_eq
      USE input
      USE fields
      USE gs_inter_block
      USE rblock
      USE tblock
      USE nimeq_free
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: f_instruct
      INTEGER(i4), INTENT(IN) :: iters
      LOGICAL, INTENT(INOUT) :: setj

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     update surface-flux computation for free-boundary solves.
c-----------------------------------------------------------------------
      IF (gs_type=="free") THEN
        CALL nimeq_free_eval(rwork1,1_i4,poly_degree,1._r8,
     $                       psio_max,psio_min)
        oflux=psio_max-psio_min
      ENDIF
c-----------------------------------------------------------------------
c     determine the range of closed poloidal-flux values and the
c     separatrix-flux value.
c-----------------------------------------------------------------------
      IF (iters>0.OR.use_delst_init.OR.gs_type=="free")
     $  CALL pflux_range(btop_check,setj)
c-----------------------------------------------------------------------
c     find the flux functions for the current flux distribution.
c-----------------------------------------------------------------------
      CALL calc_f(rwork3,"ffprime")
      CALL calc_pres(pres_eq,"pprime")
      CALL update_qfp
c-----------------------------------------------------------------------
c     if a target plasma current is specified, adjust coefficients for
c     the F profile to converge on the target.
c-----------------------------------------------------------------------
      IF (ip_equil/=0._r8.AND.(f_model=="quad_closed".OR.
     $    f_model=="cubic_closed".OR.f_model=="linear").AND.iters>0.AND.
     $    (f_instruct=="update f1".OR.f_instruct=="dFFdf only")) THEN
        IF (f_instruct=="update f1") THEN
          CALL f_adjust(curarr,r2arr,.true.,dipdf)
        ELSE IF (f_instruct=="dFFdf only") THEN
          CALL f_adjust(curarr,r2arr,.false.,dipdf)
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%rwork4,
     $                            rb(ibl)%qrwork4,rb(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%rwork4,
     $                            tb(ibl)%qrwork4,tb(ibl))
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gs_findnl
c-----------------------------------------------------------------------
c     subprogram 8. calc_pres
c     if the input flag is "pprime" then the subroutine calculates the
c     value of mu0*dP/dpsi at each of the nodes.  otherwise, just
c     compute mu*P.
c-----------------------------------------------------------------------
      SUBROUTINE calc_pres(vect,flag)
      USE local
      USE fields
      USE vector_type_mod
      IMPLICIT NONE

      TYPE(vector_type), DIMENSION(:), POINTER :: vect
      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: ibl,mxb,myb,ix,iy,ibase,dord
      REAL(r8), EXTERNAL :: p_func
c-----------------------------------------------------------------------
c     evaluate coefficients according to flag.
c     the fllen data structure holds field-line length relative to
c     the computed maximum which is used to indicated topology.
c-----------------------------------------------------------------------
      IF (flag=='pprime') THEN
        dord=1
      ELSE
        dord=0
      ENDIF
      DO ibl=1,nbl
        mxb=SIZE(vect(ibl)%arr,2)-1
        myb=SIZE(vect(ibl)%arr,3)-1
        DO iy=0,myb
          DO ix=0,mxb
            vect(ibl)%arr(1,ix,iy)=
     $        p_func(pflux(ibl)%arr(1,ix,iy),dord,
     $        fllen(ibl)%arr(1,ix,iy),.false.)
          ENDDO
        ENDDO
        IF (ASSOCIATED(vect(ibl)%arrh)) THEN
          DO iy=0,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(vect(ibl)%arrh,2)
                vect(ibl)%arrh(1,ibase,ix,iy)=
     $            p_func(pflux(ibl)%arrh(1,ibase,ix,iy),dord,
     $                   fllen(ibl)%arrh(1,ibase,ix,iy),.false.)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(vect(ibl)%arrv)) THEN
          DO iy=1,myb
            DO ix=0,mxb
              DO ibase=1,SIZE(vect(ibl)%arrv,2)
                vect(ibl)%arrv(1,ibase,ix,iy)=
     $            p_func(pflux(ibl)%arrv(1,ibase,ix,iy),dord,
     $                   fllen(ibl)%arrv(1,ibase,ix,iy),.false.)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(vect(ibl)%arri)) THEN
          DO iy=1,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(vect(ibl)%arri,2)
                vect(ibl)%arri(1,ibase,ix,iy)=
     $            p_func(pflux(ibl)%arri(1,ibase,ix,iy),dord,
     $                   fllen(ibl)%arri(1,ibase,ix,iy),.false.)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE calc_pres

c-----------------------------------------------------------------------
c     subprogram 9. p_func
c     this function returns the value or derivative of mu0*P(psi), where
c     psi is the poloidal flux.  use of this external function
c     consolidates evaluation of user-specified input.
c
c     rlen is the relative length to indicate magnetic topology, where
c     rlen<1 is considered open.  pisnormed indicates whether the
c     input psi value is normalized to the range of closed flux.
c-----------------------------------------------------------------------
      FUNCTION p_func(psi,d_order,rlen,psinormed) RESULT(pval)
      USE local
      USE input_eq
      USE input
      USE nimeq_mod
      USE spline
      IMPLICIT NONE

      REAL(r8) :: pval
      REAL(r8), INTENT(IN) :: psi,rlen
      INTEGER(i4), INTENT(IN) :: d_order
      LOGICAL, INTENT(IN) :: psinormed

      REAL(r8) :: psi_norm,dpsn,pr0
      REAL(r8) :: small=1.e-4
      INTEGER(i4) :: isq
c-----------------------------------------------------------------------
c     evaluate mu0*pressure or its derivative with respect to psi
c     according to d_order and pres_model.
c-----------------------------------------------------------------------
      SELECT CASE (pres_model)
      CASE ("constant")
        IF (d_order==0) THEN
          pval=pressure0
        ELSE
          pval=0._r8
        ENDIF
        RETURN
c-----------------------------------------------------------------------
c     the quad_closed model uses p1_closed to specify the difference
c     between p at the magnetic axis and the open p; p2_closed is a
c     quadratic shaping coefficient that does not affect the axis
c     and open pressures.
c-----------------------------------------------------------------------
      CASE ("quad_closed")
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            pval=p_open+p1_closed
          ELSE
            pval=-p1_closed/MAX(oflux,small)
          ENDIF
        ELSE IF (d_order==0) THEN
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            pval=p_open+p1_closed*(1._r8-psi_norm)+
     $            4._r8*p2_closed*psi_norm*(psi_norm-1._r8)
          ELSE
            pval=p_open
          ENDIF
        ELSE
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            pval=
     $        (4._r8*p2_closed*(2._r8*psi_norm-1._r8)-p1_closed)*dpsn
          ELSE
            pval=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     the cubic_closed model uses one of the Hermite cubic basis
c     functions so that derivatives wrt psi are zero at the axis and
c     at the separatrix.
c-----------------------------------------------------------------------
      CASE ("cubic_closed")
        pr0=p_axis-p_open
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            pval=p_axis
          ELSE
            pval=-pr0/MAX(oflux,small)
          ENDIF
        ELSE IF (d_order==0) THEN
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            pval=pr0*(1._r8-3._r8*psi_norm**2+2._r8*psi_norm**3)+p_open
          ELSE
            pval=p_open
          ENDIF
        ELSE
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            pval=6._r8*pr0*psi_norm*(psi_norm-1._r8)*dpsn
          ELSE
            pval=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     use data from fluxgrid's 1d.bin file.
c-----------------------------------------------------------------------
      CASE ("read1d")
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            pval=neqsq%fs(0,3)
          ELSE
            pval=-neqsq%fs(0,3)
          ENDIF
        ELSE IF (psi_norm>1._r8.OR.rlen<=rlencl) THEN
          IF (d_order==0) THEN
            pval=neqsq%fs(neqsq%nodes,3)+p_open
          ELSE
            pval=0._r8
          ENDIF
        ELSE
          CALL spline_eval(neqsq,psi_norm,1_i4)
          IF (d_order==0) THEN
            pval=neqsq%f(3)+p_open
          ELSE
            pval=neqsq%f1(3)*dpsn
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     exit if specification is not recognized.
c-----------------------------------------------------------------------
      CASE DEFAULT
        CALL nim_stop("The pres_model value is not recognized.")
      END SELECT
c-----------------------------------------------------------------------
c     if specified, use the Hermite cubic expansion to smooth p near
c     the last closed flux surface.
c-----------------------------------------------------------------------
      IF (psin_smthpp<1._r8.AND.psi_norm>psin_smthpp.AND.
     $    psi_norm<1._r8.AND.rlen>rlencl.AND.psmthset) THEN
        psi_norm=(psi_norm-psin_smthpp)
        IF (d_order==0) THEN
          pval=pvalin* (1._r8-3._r8*(psi_norm/pdpsi)**2+
     $                  2._r8*(psi_norm/pdpsi)**3)+
     $         pvalout*(3._r8*(psi_norm/pdpsi)**2-
     $                  2._r8*(psi_norm/pdpsi)**3)+
     $         pdrvin* (psi_norm-2._r8*psi_norm**2/pdpsi+
     $                  psi_norm**3/pdpsi**2)
        ELSE
          pval=(pvalin* (-6._r8*psi_norm/pdpsi**2+
     $                    6._r8*psi_norm**2/pdpsi**3)+
     $          pvalout*(6._r8*psi_norm/pdpsi**2-
     $                   6._r8*psi_norm**2/pdpsi**3)+
     $          pdrvin* (1._r8-4._r8*psi_norm/pdpsi+
     $                   3._r8*(psi_norm/pdpsi)**2))*dpsn
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION p_func

c-----------------------------------------------------------------------
c     subprogram 10. calc_f
c     if the input flag is "ffprime" then the subroutine calculates the
c     value of FF' at each of the nodes.  if it is "fprime", just
c     return F', and if it is "f only", return F.
c-----------------------------------------------------------------------
      SUBROUTINE calc_f(vect,flag)
      USE local
      USE fields
      USE vector_type_mod
      IMPLICIT NONE

      TYPE(vector_type), DIMENSION(:), POINTER :: vect
      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: ibl,mxb,myb,ix,iy,ibase,iord
      REAL(r8) :: val
      REAL(r8), EXTERNAL :: f_func
c-----------------------------------------------------------------------
c     evaluate coefficients according to flag.
c     the fllen data structure holds field-line length relative to
c     the computed maximum which is used to indicated topology.
c-----------------------------------------------------------------------
      IF (flag=='fprime') THEN
        iord=1
      ELSE
        iord=0
      ENDIF

      DO ibl=1,nbl
        mxb=SIZE(vect(ibl)%arr,2)-1
        myb=SIZE(vect(ibl)%arr,3)-1
        DO iy=0,myb
          DO ix=0,mxb
            val=f_func(pflux(ibl)%arr(1,ix,iy),iord,
     $                 fllen(ibl)%arr(1,ix,iy),.false.)
            IF (flag=='ffprime')
     $        val=val*f_func(pflux(ibl)%arr(1,ix,iy),1_i4,
     $                       fllen(ibl)%arr(1,ix,iy),.false.)
            vect(ibl)%arr(1,ix,iy)=val
          ENDDO
        ENDDO
        IF (ASSOCIATED(vect(ibl)%arrh)) THEN
          DO iy=0,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(vect(ibl)%arrh,2)
                val=f_func(pflux(ibl)%arrh(1,ibase,ix,iy),iord,
     $                     fllen(ibl)%arrh(1,ibase,ix,iy),.false.)
                IF (flag=='ffprime')
     $            val=val*f_func(pflux(ibl)%arrh(1,ibase,ix,iy),1_i4,
     $                           fllen(ibl)%arrh(1,ibase,ix,iy),.false.)
                vect(ibl)%arrh(1,ibase,ix,iy)=val
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(vect(ibl)%arrv)) THEN
          DO iy=1,myb
            DO ix=0,mxb
              DO ibase=1,SIZE(vect(ibl)%arrv,2)
                val=f_func(pflux(ibl)%arrv(1,ibase,ix,iy),iord,
     $                     fllen(ibl)%arrv(1,ibase,ix,iy),.false.)
                IF (flag=='ffprime')
     $            val=val*f_func(pflux(ibl)%arrv(1,ibase,ix,iy),1_i4,
     $                           fllen(ibl)%arrv(1,ibase,ix,iy),.false.)
                vect(ibl)%arrv(1,ibase,ix,iy)=val
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(vect(ibl)%arri)) THEN
          DO iy=1,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(vect(ibl)%arri,2)
                val=f_func(pflux(ibl)%arri(1,ibase,ix,iy),iord,
     $                     fllen(ibl)%arri(1,ibase,ix,iy),.false.)
                IF (flag=='ffprime')
     $            val=val*f_func(pflux(ibl)%arri(1,ibase,ix,iy),1_i4,
     $                           fllen(ibl)%arri(1,ibase,ix,iy),.false.)
                vect(ibl)%arri(1,ibase,ix,iy)=val
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE calc_f

c-----------------------------------------------------------------------
c     subprogram 11. f_func
c     this function returns the value or derivative of f(psi), where
c     psi is the poloidal flux, and f is R*B_phi.  use of this
c     external function consolidates evaluation of user-specified
c     input.
c
c     rlen is the relative length to indicate magnetic topology, where
c     rlen<1 is considered open.  pisnormed indicates whether the
c     input psi value is normalized to the range of closed flux.
c-----------------------------------------------------------------------
      FUNCTION f_func(psi,d_order,rlen,psinormed) RESULT(fval)
      USE local
      USE input_eq
      USE nimeq_mod
      IMPLICIT NONE

      REAL(r8) :: fval
      REAL(r8), INTENT(IN) :: psi,rlen
      INTEGER(i4), INTENT(IN) :: d_order
      LOGICAL, INTENT(IN) :: psinormed

      REAL(r8) :: psi_norm,dpsn
      REAL(r8) :: small=1.e-4
      INTEGER(i4) :: isq
c-----------------------------------------------------------------------
c     evaluate f or its derivative according to d_order and f_model.
c
c     the "linear" model uses the un-normalized flux value.
c-----------------------------------------------------------------------
      SELECT CASE (f_model)
      CASE ("linear")
        IF (d_order==0) THEN
          fval=f0+dfdpsi*psi
        ELSE
          fval=dfdpsi
        ENDIF 
        RETURN
c-----------------------------------------------------------------------
c     quadratic dependence on the normalized closed flux.
c-----------------------------------------------------------------------
      CASE ("quad_closed")
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            fval=f_open+f1_closed
          ELSE
            fval=-f1_closed/MAX(oflux,small)
          ENDIF
        ELSE IF (d_order==0) THEN
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            fval=f_open+f1_closed*(1._r8-psi_norm)+
     $            4._r8*f2_closed*psi_norm*(psi_norm-1._r8)
          ELSE
            fval=f_open
          ENDIF
        ELSE
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            fval=
     $        (4._r8*f2_closed*(2._r8*psi_norm-1._r8)-f1_closed)*dpsn
          ELSE
            fval=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     cubic dependence on the normalized closed flux.
c-----------------------------------------------------------------------
      CASE ("cubic_closed")
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            fval=f_axis
          ELSE
            fval=(f_open-f_axis)/MAX(oflux,small)
          ENDIF
        ELSE IF (d_order==0) THEN
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            fval=f_open+(f_axis-f_open)*
     $                  (1._r8-3._r8*psi_norm**2+2._r8*psi_norm**3)
          ELSE
            fval=f_open
          ENDIF
        ELSE
          IF (psi_norm<=1._r8.AND.rlen>rlencl) THEN
            fval=6._r8*(f_axis-f_open)*psi_norm*(psi_norm-1._r8)*dpsn
          ELSE
            fval=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     use data from fluxgrid's 1d.bin file.
c-----------------------------------------------------------------------
      CASE ("read1d")
        CALL normalize_psi(psi,psi_norm,dpsn,psinormed)
        IF (dpsn==0.) THEN	! no closed flux -- rough guess to get
				! started assuming d(psi_norm)/d(psi)>0
          IF (d_order==0) THEN
            fval=neqsq%fs(0,2)
          ELSE
            fval=0
          ENDIF
        ELSE IF (psi_norm>1._r8.OR.rlen<=rlencl) THEN
          IF (d_order==0) THEN
            fval=neqsq%fs(neqsq%nodes,2)
          ELSE
            fval=0._r8
          ENDIF
        ELSE
          CALL spline_eval(neqsq,psi_norm,1_i4)
          IF (d_order==0) THEN
            fval=neqsq%f(2)
          ELSE
            fval=neqsq%f1(2)*dpsn
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     exit if specification is not recognized.
c-----------------------------------------------------------------------
      CASE DEFAULT
        CALL nim_stop("The f_model value is not recognized.")
      END SELECT
c-----------------------------------------------------------------------
c     if specified, use the Hermite cubic expansion to smooth f near
c     the last closed flux surface.
c-----------------------------------------------------------------------
      IF (psin_smthfp<1._r8.AND.psi_norm>psin_smthfp.AND.
     $    psi_norm<1._r8.AND.rlen>rlencl.AND.fsmthset) THEN
        psi_norm=(psi_norm-psin_smthfp)
        IF (d_order==0) THEN
          fval=fvalin* (1._r8-3._r8*(psi_norm/fdpsi)**2+
     $                  2._r8*(psi_norm/fdpsi)**3)+
     $         fvalout*(3._r8*(psi_norm/fdpsi)**2-
     $                  2._r8*(psi_norm/fdpsi)**3)+
     $         fdrvin* (psi_norm-2._r8*psi_norm**2/fdpsi+
     $                  psi_norm**3/fdpsi**2)
        ELSE
          fval=(fvalin* (-6._r8*psi_norm/fdpsi**2+
     $                    6._r8*psi_norm**2/fdpsi**3)+
     $          fvalout*(6._r8*psi_norm/fdpsi**2-
     $                   6._r8*psi_norm**2/fdpsi**3)+
     $          fdrvin* (1._r8-4._r8*psi_norm/fdpsi+
     $                   3._r8*(psi_norm/fdpsi)**2))*dpsn
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION f_func

c-----------------------------------------------------------------------
c     subprogram 12. f_adjust
c     this subroutine adjusts f1_closed or f1_axis to achieve
c     a desired plasma current.
c
c     if change_f1 is false, find the partial derivatives with respect
c     to changes in f1_closed or f1_axis, but do not update their
c     values.  partial(FF')/partial(f1) at the finite-element nodes is
c     saved in rwork4 so that it can be used elsewhere.  note that
c     rwork3 is not over-written when change_f1 is false.
c-----------------------------------------------------------------------
      SUBROUTINE f_adjust(plcarr,r2arr,change_f1,dipdf)
      USE local
      USE physdat
      USE input
      USE input_eq
      USE fields
      USE rblock
      USE tblock
      USE pardata
      IMPLICIT NONE

      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plcarr,r2arr
      REAL(r8), INTENT(OUT) :: dipdf
      LOGICAL, INTENT(IN) :: change_f1

      REAL(r8) :: fa_ip,f1old,f2save,fosave,fsmth_save
      REAL(r8), PARAMETER :: deli_lim=0.75_r8
      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     interface block for the calc_pres routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_pres(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_pres
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the calc_f routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_f(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_f
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the plasma_cur subroutine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE plasma_cur(muz,plcur,plcar,r2ar)
        USE local
        USE fields
        IMPLICIT NONE
        REAL(r8), INTENT(IN) :: muz
        REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plcar,r2ar
        REAL(r8), INTENT(OUT) :: plcur
        END SUBROUTINE plasma_cur
      END INTERFACE
c-----------------------------------------------------------------------
c     find the current current.
c-----------------------------------------------------------------------
      CALL plasma_cur(mu0,fa_ip,plcarr,r2arr)
c-----------------------------------------------------------------------
c     vary f1_closed or f1_axis and recompute current to find
c     d(Ip)/df_parameter.  four f_func evaluations for F, F', dF and
c     dF' are needed to create F*dF'+dF*F'.  prese_eq, tele_eq, and
c     tion_eq are used for temporary storage.
c-----------------------------------------------------------------------
      CALL calc_f(rwork4,"f only")
      CALL calc_f(prese_eq,"fprime")
      SELECT CASE (f_model)
      CASE("quad_closed")
        f1old=f1_closed
        f2save=f2_closed
        fosave=f_open
        f1_closed=1._r8
        f2_closed=0._r8
        f_open=0._r8
      CASE("cubic_closed")
        f1old=f_axis
        fosave=f_open
        f_axis=1._r8
        f_open=0._r8
      CASE("linear")
        f1old=dfdpsi
        fosave=f0
        dfdpsi=1._r8
        f0=0._r8
      END SELECT
      fsmth_save=psin_smthfp
      psin_smthfp=10. ! preclude edge smoothing
      CALL calc_f(tion_eq,"fprime")
      CALL calc_f(tele_eq,"f only")
      psin_smthfp=fsmth_save

      DO ibl=1,nbl
        rwork4(ibl)%arr=rwork4(ibl)%arr*tion_eq(ibl)%arr+
     $                prese_eq(ibl)%arr*tele_eq(ibl)%arr
        IF (poly_degree>1.AND.ibl<=nrbl) THEN
          rwork4(ibl)%arrh=rwork4(ibl)%arrh*tion_eq(ibl)%arrh+
     $                   prese_eq(ibl)%arrh*tele_eq(ibl)%arrh
          rwork4(ibl)%arrv=rwork4(ibl)%arrv*tion_eq(ibl)%arrv+
     $                   prese_eq(ibl)%arrv*tele_eq(ibl)%arrv
          rwork4(ibl)%arri=rwork4(ibl)%arri*tion_eq(ibl)%arri+
     $                   prese_eq(ibl)%arri*tele_eq(ibl)%arri
        ENDIF
        IF (ibl<=nrbl) THEN
          CALL rblock_qp_update(rb(ibl)%rwork4,
     $                          rb(ibl)%qrwork1,rb(ibl))
          rb(ibl)%qpres_eq%qpf=0._r8
        ELSE
          CALL tblock_qp_update(tb(ibl)%rwork4,
     $                          tb(ibl)%qrwork1,tb(ibl))
          tb(ibl)%qpres_eq%qpf=0._r8
        ENDIF
      ENDDO
      CALL plasma_cur(mu0,dipdf,plcarr,r2arr)
c-----------------------------------------------------------------------
c     adjust the f parameter to achieve the desired current, and
c     update or reset the F*F' distribution and quad-point storage.
c-----------------------------------------------------------------------
      IF (change_f1) THEN
        IF (gs_type=='free') THEN
          f1old=f1old+deli_lim*(ip_equil-fa_ip)/dipdf
        ELSE
          f1old=f1old+(ip_equil-fa_ip)/dipdf
        ENDIF
      ENDIF

      SELECT CASE (f_model)
      CASE("quad_closed")
        f1_closed=f1old
        f2_closed=f2save
        f_open=fosave
      CASE("cubic_closed")
        f_axis=f1old
        f_open=fosave
      CASE("linear")
        dfdpsi=f1old
        f0=fosave
      END SELECT
      CALL calc_f(rwork3,"ffprime")
      CALL update_qfp  !  always needed regardless of change_f1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE f_adjust

c-----------------------------------------------------------------------
c     subprogram 13. pflux_range
c     this subroutine determines the range of closed-flux values based
c     on the current solution and the min and max open-flux values.
c-----------------------------------------------------------------------
      SUBROUTINE pflux_range(btop_ch,setj)
      USE local
      USE fields
      USE nimeq_mod
      USE mpi_nim
      USE pardata
      USE input_eq
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: btop_ch
      LOGICAL, INTENT(INOUT) :: setj

      INTEGER(i4) :: ibl,ierror
      INTEGER(i4), SAVE :: jdir
      REAL(r8) :: pfmin,pfmax,rtmp
      REAL(r8) :: small=1.e-6
c-----------------------------------------------------------------------
c     loop over blocks to get min and max of the current solution.
c-----------------------------------------------------------------------
      pfmin=HUGE(0._r8)
      pfmax=-pfmin
      DO ibl=1,nbl
        pfmin=MIN(pfmin,MINVAL(pflux(ibl)%arr))
        pfmax=MAX(pfmax,MAXVAL(pflux(ibl)%arr))
        IF (ASSOCIATED(pflux(ibl)%arrh)) THEN
          pfmin=MIN(pfmin,MINVAL(pflux(ibl)%arrh))
          pfmax=MAX(pfmax,MAXVAL(pflux(ibl)%arrh))
        ENDIF
        IF (ASSOCIATED(pflux(ibl)%arrv)) THEN
          pfmin=MIN(pfmin,MINVAL(pflux(ibl)%arrv))
          pfmax=MAX(pfmax,MAXVAL(pflux(ibl)%arrv))
        ENDIF
        IF (ASSOCIATED(pflux(ibl)%arri)) THEN
          pfmin=MIN(pfmin,MINVAL(pflux(ibl)%arri))
          pfmax=MAX(pfmax,MAXVAL(pflux(ibl)%arri))
        ENDIF
      ENDDO
      IF (nprocs>1) THEN
        CALL mpi_allreduce(pfmin,rtmp,1,mpi_nim_real,mpi_min,
     $       mpi_comm_world,ierror)
        pfmin=rtmp
        CALL mpi_allreduce(pfmax,rtmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        pfmax=rtmp
      ENDIF
c-----------------------------------------------------------------------
c     determine the values corresponding to closed flux.
c-----------------------------------------------------------------------
      IF (oflux<=small*(pfmax-pfmin)) THEN  !  no open flux
        psic_min=pfmin
        psic_max=pfmax
        IF (0.5*(psic_max+psic_min)<psio_min) THEN
          sep_flux=psic_max
        ELSE
          sep_flux=psic_min
        ENDIF
        open_flux=.false.
      ELSE
        open_flux=.true.
c-----------------------------------------------------------------------
c       if it is available, use the relative length of field lines to
c       determine which points to consider for the range of poloidal-
c       flux in the closed flux region.
c-----------------------------------------------------------------------
        IF (btop_ch/="none") THEN
          psic_min=HUGE(0._r8)
          psic_max=-psic_min
          DO ibl=1,nbl
            psic_min=MIN(psic_min,MINVAL(pflux(ibl)%arr,
     $                   mask=fllen(ibl)%arr>rlencl))
            psic_max=MAX(psic_max,MAXVAL(pflux(ibl)%arr,
     $                   mask=fllen(ibl)%arr>rlencl))
            IF (ASSOCIATED(pflux(ibl)%arrh)) THEN
              psic_min=MIN(psic_min,MINVAL(pflux(ibl)%arrh,
     $                     mask=fllen(ibl)%arrh>rlencl))
              psic_max=MAX(psic_max,MAXVAL(pflux(ibl)%arrh,
     $                     mask=fllen(ibl)%arrh>rlencl))
            ENDIF
            IF (ASSOCIATED(pflux(ibl)%arrv)) THEN
              psic_min=MIN(psic_min,MINVAL(pflux(ibl)%arrv,
     $                     mask=fllen(ibl)%arrv>rlencl))
              psic_max=MAX(psic_max,MAXVAL(pflux(ibl)%arrv,
     $                     mask=fllen(ibl)%arrv>rlencl))
            ENDIF
            IF (ASSOCIATED(pflux(ibl)%arri)) THEN
              psic_min=MIN(psic_min,MINVAL(pflux(ibl)%arri,
     $                     mask=fllen(ibl)%arri>rlencl))
              psic_max=MAX(psic_max,MAXVAL(pflux(ibl)%arri,
     $                     mask=fllen(ibl)%arri>rlencl))
            ENDIF
          ENDDO
          IF (nprocs>1) THEN
            CALL mpi_allreduce(psic_min,rtmp,1,mpi_nim_real,mpi_min,
     $           mpi_comm_world,ierror)
            psic_min=rtmp
            CALL mpi_allreduce(psic_max,rtmp,1,mpi_nim_real,mpi_max,
     $           mpi_comm_world,ierror)
            psic_max=rtmp
          ENDIF
          IF (setj) THEN
            IF (ip_equil/=0._r8) THEN
              jdir=NINT(ip_equil/ABS(ip_equil))
            ELSE IF ((pfmax-psio_max)>small*oflux) THEN
              jdir=-1
            ELSE
              jdir=1
            ENDIF
            setj=.false.
          ENDIF
          sep_flux=0.5_r8*((1_i4-jdir)*psic_min+(1_i4+jdir)*psic_max)
c-----------------------------------------------------------------------
c       otherwise, use the range of flux values, but note that the
c       range of closed flux values and the separatrix flux will be
c       incorrect if there is private flux.
c-----------------------------------------------------------------------
        ELSE
          IF ((pfmax-psio_max)>small*oflux) THEN
            psic_min=psio_max
            psic_max=pfmax
            sep_flux=psio_max
          ELSE IF ((psio_min-pfmin)>small*oflux) THEN
            psic_min=pfmin
            psic_max=psio_min
            sep_flux=psio_min
          ELSE
            psic_min=psio_min
            psic_max=psio_min
            sep_flux=psio_min
          ENDIF
        ENDIF
      ENDIF
      cflux=psic_max-psic_min
      IF (cflux<=small*(pfmax-pfmin)) THEN  !  no closed flux
        closed_flux=.false.
      ELSE
        closed_flux=.true.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pflux_range

c-----------------------------------------------------------------------
c     subprogram 14. jeq_poloidal
c     find the poloidal components of mu0*J_eq from B_eq and F' using
c     mu0*J_poloidal=-(dF/dpsi)*B_poloidal.
c-----------------------------------------------------------------------
      SUBROUTINE jeq_poloidal
      USE local
      USE fields
      IMPLICIT NONE

      INTEGER(i4) :: ibl,mxb,myb,ix,iy,ibase
      REAL(r8) :: coef
      REAL(r8), EXTERNAL :: f_func
c-----------------------------------------------------------------------
c     evaluate coefficients according to flag.
c     the fllen data structure holds field-line length relative to
c     the computed maximum which is used to indicated topology.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        mxb=SIZE(ja_eq(ibl)%arr,2)-1
        myb=SIZE(ja_eq(ibl)%arr,3)-1
        DO iy=0,myb
          DO ix=0,mxb
            coef=-f_func(pflux(ibl)%arr(1,ix,iy),1_i4,
     $                   fllen(ibl)%arr(1,ix,iy),.false.)
            ja_eq(ibl)%arr(1:2,ix,iy)=coef*be_eq(ibl)%arr(1:2,ix,iy)
          ENDDO
        ENDDO
        IF (ASSOCIATED(ja_eq(ibl)%arrh)) THEN
          DO iy=0,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(ja_eq(ibl)%arrh,2)
                coef=-f_func(pflux(ibl)%arrh(1,ibase,ix,iy),1_i4,
     $                       fllen(ibl)%arrh(1,ibase,ix,iy),.false.)
                ja_eq(ibl)%arrh(1:2,ibase,ix,iy)=
     $             coef*be_eq(ibl)%arrh(1:2,ibase,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(ja_eq(ibl)%arrv)) THEN
          DO iy=1,myb
            DO ix=0,mxb
              DO ibase=1,SIZE(ja_eq(ibl)%arrv,2)
                coef=-f_func(pflux(ibl)%arrv(1,ibase,ix,iy),1_i4,
     $                       fllen(ibl)%arrv(1,ibase,ix,iy),.false.)
                ja_eq(ibl)%arrv(1:2,ibase,ix,iy)=
     $             coef*be_eq(ibl)%arrv(1:2,ibase,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        IF (ASSOCIATED(ja_eq(ibl)%arri)) THEN
          DO iy=1,myb
            DO ix=1,mxb
              DO ibase=1,SIZE(ja_eq(ibl)%arri,2)
                coef=-f_func(pflux(ibl)%arri(1,ibase,ix,iy),1_i4,
     $                       fllen(ibl)%arri(1,ibase,ix,iy),.false.)
                ja_eq(ibl)%arri(1:2,ibase,ix,iy)=
     $             coef*be_eq(ibl)%arri(1:2,ibase,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jeq_poloidal

c-----------------------------------------------------------------------
c     subprogram 15. normalize_psi.
c     use this subroutine to uniquely define the normalization of
c     poloidal flux for the equilibrium profile options.
c-----------------------------------------------------------------------
      SUBROUTINE normalize_psi(psi,psi_norm,dpsn,normdone)
      USE local
      USE input_eq
      USE nimeq_mod
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: psi
      REAL(r8), INTENT(OUT) :: psi_norm,dpsn
      LOGICAL, INTENT(IN) :: normdone
c-----------------------------------------------------------------------
c     it is convenient to allow the calling routine an option for
c     accepting already normalized values.
c-----------------------------------------------------------------------
      IF (normdone) THEN
        psi_norm=psi
        dpsn=1._r8
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     find the normalized value of the flux function based on the limits
c     of the closed flux values.  also find d(psi_norm)/d(psi).
c-----------------------------------------------------------------------
      IF (psic_max-sep_flux>0.5_r8*cflux) THEN
        psi_norm=1._r8-(psi-sep_flux)/cflux
        dpsn=-1._r8/cflux
      ELSE IF (sep_flux-psic_min>0.5*cflux) THEN
        psi_norm=1._r8-(sep_flux-psi)/cflux
        dpsn=1._r8/cflux
      ELSE  !  no closed flux -- set dpsn=0 as a flag.
        psi_norm=psi
        dpsn=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE normalize_psi

c-----------------------------------------------------------------------
c     subprogram 16. calc_blen
c     find the length of equilibrium field-lines, starting from each
c     node of the finite-element expansion.
c
c     the input refln is the maximum length for the field-line
c     integration and will be used to test whether each line is open.
c-----------------------------------------------------------------------
      SUBROUTINE calc_blen(refln,use_stiff)
      USE local
      USE fields
      USE vector_type_mod
      USE edge
      USE nimeq_btr
      USE input
      USE input_eq
      USE nimeq_mod
      USE pardata
      USE nimeq_all
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: refln
      LOGICAL, INTENT(IN) :: use_stiff

      INTEGER(i4) :: ibl,mxb,myb,ix,iy,ibase,iv,ivp,ixp,iyp,ixs,iys
      REAL(r8) :: xb,yb,len,psi_norm,dpsn,pfl,pflp,rr,zz
      REAL(r8), PARAMETER :: offset=1.e-2_r8  !  fraction of 1 cell
      REAL(r8), PARAMETER :: pcheck=0.80_r8   !  flux-value limit
      REAL(r8), PARAMETER :: ptola=1.e-9_r8   !  abs max/min tolerance
c-----------------------------------------------------------------------
c     parallel computations need to update the flux values in the
c     assembled blocks (rball, tball) for field-line tracing.
c-----------------------------------------------------------------------
      IF (nprocs>1) CALL nimeq_all_update
c-----------------------------------------------------------------------
c     also check for variation in the flux value along the domain
c     boundary, signifying open field there.  rwork3 is used for
c     temporary storage without conflicting the FF' storage elsewhere
c     in the loop.
c-----------------------------------------------------------------------
c     DO ibl=1,nbl
c       rwork3(ibl)=0._r8
c       ivp=seam(ibl)%nvert
c       ixp=seam(ibl)%vertex(ivp)%intxy(1)
c       iyp=seam(ibl)%vertex(ivp)%intxy(2)
c       DO iv=1,seam(ibl)%nvert
c         ix=seam(ibl)%vertex(iv)%intxy(1)
c         iy=seam(ibl)%vertex(iv)%intxy(2)
c         IF (seam(ibl)%expoint(iv).AND.
c    $        seam(ibl)%expoint(ivp)) THEN
c           pfl=be_n0(ibl)%arr(1,ix,iy)*be_n0(ibl)%arr(2,ix,iy)**2
c           pflp=be_n0(ibl)%arr(1,ixp,iyp)*be_n0(ibl)%arr(2,ixp,iyp)**2
c           IF (ABS(pfl-pflp)>ptola*MAX(ptola,ABS(pfl),ABS(pflp))) THEN
c             rwork3(ibl)%arr(1,ix,iy)=1
c             rwork3(ibl)%arr(1,ixp,iyp)=1
c             IF (ibl<=nrbl.AND.poly_degree>1) THEN
c               ixs=seam(ibl)%segment(iv)%intxys(1)
c               iys=seam(ibl)%segment(iv)%intxys(2)
c               IF (seam(ibl)%segment(iv)%h_side) THEN
c                 rwork3(ibl)%arrh(1,:,ixs,iys)=1
c               ELSE
c                 rwork3(ibl)%arrv(1,:,ixs,iys)=1
c               ENDIF
c             ENDIF
c           ELSE
c             rwork3(ibl)%arr(1,ix,iy)=-100
c             rwork3(ibl)%arr(1,ixp,iyp)=-100
c             IF (ibl<=nrbl.AND.poly_degree>1) THEN
c               ixs=seam(ibl)%segment(iv)%intxys(1)
c               iys=seam(ibl)%segment(iv)%intxys(2)
c               IF (seam(ibl)%segment(iv)%h_side) THEN
c                 rwork3(ibl)%arrh(1,:,ixs,iys)=-100
c               ELSE
c                 rwork3(ibl)%arrv(1,:,ixs,iys)=-100
c               ENDIF
c             ENDIF
c           ENDIF
c         ENDIF
c         ivp=iv
c         ixp=ix
c         iyp=iy
c       ENDDO
c     ENDDO
c-----------------------------------------------------------------------
c     share nearby flux values across block borders.
c-----------------------------------------------------------------------
c     DO ibl=1,nbl
c       CALL edge_load_arr(rwork3(ibl),1_i4,poly_degree-1_i4,seam(ibl))
c     ENDDO
c     CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.false.)
c     DO ibl=1,nbl
c       CALL edge_unload_arr(rwork3(ibl),1_i4,poly_degree-1_i4,
c    $                       seam(ibl))
c     ENDDO
c-----------------------------------------------------------------------
c     loop over all blocks and call the tracing function from the
c     nimeq_btr module for each node.  save the results in the fllen
c     data structure.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb=SIZE(fllen(ibl)%arr,2)-1
        myb=SIZE(fllen(ibl)%arr,3)-1
        DO ibase=1,SIZE(rb(ibl)%rz%ix0)
          DO iy=rb(ibl)%rz%iy0(ibase),myb
            DO ix=rb(ibl)%rz%ix0(ibase),mxb
              xb=ix-rb(ibl)%rz%ix0(ibase)+rb(ibl)%rz%dx(ibase)
              IF (rb(ibl)%degenerate) xb=MAX(offset,xb)
              yb=iy-rb(ibl)%rz%iy0(ibase)+rb(ibl)%rz%dy(ibase)
              CALL lagr_quad_eval(rb(ibl)%rz,xb,yb,0_i4)
              CALL lagr_quad_eval(rb(ibl)%fllen,xb,yb,0_i4)
              CALL lagr_quad_eval(rb(ibl)%pflux,xb,yb,0_i4)
              CALL normalize_psi(rb(ibl)%pflux%f(1),psi_norm,dpsn,
     $                           .false.)
              IF (psi_norm>pcheck.OR.rb(ibl)%fllen%f(1)<rlencl) THEN
                CALL lagr_quad_eval(rb(ibl)%rwork3,xb,yb,0_i4)
c-TEST
c               IF (rb(ibl)%rwork3%f(1)>ptola) THEN
c                 len=0._r8  !  open-field on boundary
c               ELSE
                  rr=rb(ibl)%rz%f(1)  !  passing local variables avoids
                  zz=rb(ibl)%rz%f(2)  !  strange ptr error in nimeq_tr
                  len=beq_tr(rr,zz,refln,xb,yb,ibl,use_stiff)/
     $                (refln*blen_check)
c               ENDIF
                IF (ibase==1) THEN
                  fllen(ibl)%arr(1,ix,iy)=len
                ELSE IF (ibase<=poly_degree) THEN
                  fllen(ibl)%arrh(1,ibase-1,ix,iy)=len
                ELSE IF (ibase<2*poly_degree) THEN
                  fllen(ibl)%arrv(1,ibase-poly_degree,ix,iy)=len
                ELSE
                  fllen(ibl)%arri(1,ibase-2*poly_degree+1,ix,iy)=len
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        rwork3(ibl)=0._r8
      ENDDO

      DO ibl=nrbl+1,nbl
        mxb=SIZE(fllen(ibl)%arr,2)-1
        DO ix=0,mxb
          rr=tb(ibl)%tgeom%xs(ix)
          zz=tb(ibl)%tgeom%ys(ix)
          fllen(ibl)%arr(1,ix,0_i4)=
     $      beq_tr(rr,zz,refln,0._r8,0._r8,ibl,use_stiff)/
     $      (refln*blen_check)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     average along borders between blocks to ensure a unique value
c     at each node.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL edge_load_arr(fllen(ibl),1_i4,poly_degree-1_i4,
     $                     seam(ibl))
        DO ix=1,seam(ibl)%nvert
          seam(ibl)%vertex(ix)%seam_in(1)=
     $      seam(ibl)%vertex(ix)%seam_in(1)*
     $      seam(ibl)%vertex(ix)%ave_factor
          IF (ibl<=nrbl) THEN
            DO iy=1,poly_degree-1
              seam(ibl)%segment(ix)%seam_in(iy)=
     $          seam(ibl)%segment(ix)%seam_in(iy)*
     $          seam(ibl)%segment(ix)%ave_factor
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.false.)
      DO ibl=1,nbl
        CALL edge_unload_arr(fllen(ibl),1_i4,poly_degree-1_i4,
     $                       seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE calc_blen
c-----------------------------------------------------------------------
c     subprogram 17. plasma_cur
c     this subroutine computes the net toroidal current from current
c     density at the quadrature points of the numerical integration.
c     it allows two contributions: one from part of -mu0*R*J_phi at
c     quadrature points, saved in the qrwork1 data structure in
c     each block (from F*F'), and one from the rest of -mu0*J_phi/R in
c     qpres_eq (from P').  this is convenient for the GS iteration,
c     although this routine can be used elsewhere.
c-----------------------------------------------------------------------
      SUBROUTINE plasma_cur(muz,plcur,plcarr,r2arr)
      USE local
      USE fields
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: muz
      REAL(r8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: plcarr,r2arr
      REAL(r8), INTENT(OUT) :: plcur

      INTEGER(i4) :: njphi,jb,jcol,jelm,jx,jy,jg,ierror
      REAL(r8) :: rtmp
c-----------------------------------------------------------------------
c     if the geometric arrays have not been created, construct them.
c     the operations are based on those for the psi_from_j matrix in
c     nimeq_free.
c-----------------------------------------------------------------------
      IF (.NOT.ALLOCATED(plcarr)) THEN
        njphi=0
        DO jb=1,nrbl
          njphi=njphi+SIZE(rb(jb)%wjac)
        ENDDO
        DO jb=nrbl+1,nbl
          njphi=njphi+tb(jb)%mcell*tb(jb)%ng
        ENDDO
        ALLOCATE(plcarr(njphi),r2arr(njphi))

        jcol=0
        DO jb=1,nrbl
          jelm=1
          DO jy=0,rb(jb)%my-1
            DO jx=0,rb(jb)%mx-1
              DO jg=1,rb(jb)%ng
                jcol=jcol+1
                CALL lagr_quad_eval(rb(jb)%rz,jx+rb(jb)%xg(jg),
     $                              jy+rb(jb)%yg(jg),0_i4)
                plcarr(jcol)=
     $            rb(jb)%wjac(jg,jelm)/rb(jb)%bigr(jg,jelm)**2
                r2arr(jcol)=rb(jb)%bigr(jg,jelm)**2
              ENDDO
              jelm=jelm+1
            ENDDO
          ENDDO
        ENDDO

c-PRE
        DO jb=nrbl+1,nbl
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     sum the contributions to net current.
c-----------------------------------------------------------------------
      plcur=0._r8
      jcol=0
      DO jb=1,nrbl
        jelm=1
        DO jy=0,rb(jb)%my-1
          DO jx=0,rb(jb)%mx-1
            DO jg=1,rb(jb)%ng
              jcol=jcol+1
              plcur=plcur-plcarr(jcol)*
     $              (rb(jb)%qrwork1% qpf(1,jg,jelm)+
     $               rb(jb)%qpres_eq%qpf(1,jg,jelm)*r2arr(jcol))
            ENDDO
            jelm=jelm+1
          ENDDO
        ENDDO
      ENDDO
c-PRE
      DO jb=nrbl+1,nbl
      ENDDO
c-----------------------------------------------------------------------
c     sum contributions over processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(plcur,rtmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        plcur=rtmp
      ENDIF
      plcur=plcur/muz
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plasma_cur
c-----------------------------------------------------------------------
c     subprogram 18. par_cont_init
c     this subroutine initializes binary files for contour plots,
c     as used by nimeq when running in parallel.
c-----------------------------------------------------------------------
      SUBROUTINE par_cont_init(file,nq,ninterp)
      USE local
      USE fields
      USE mpi_nim
      USE pardata
      USE contour_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: file
      INTEGER(i4), INTENT(IN) :: nq,ninterp

      INTEGER(i4) :: ibl,ibloc,ierror
c-----------------------------------------------------------------------
c     follow the write sequence in contour_init with one processor
c     writing its information at a time.
c-----------------------------------------------------------------------
      IF (node==0) THEN
        CALL open_bin(con_unit,TRIM(file),"UNKNOWN","REWIND",32_i4)
        WRITE(con_unit) INT(nrbl_total,4),INT(nbl_total-nrbl_total,4),
     $                  INT(nq,4)
        CALL close_bin(con_unit,TRIM(file))
      ENDIF
      CALL mpi_barrier(mpi_comm_world,ierror)
      DO ibl=1,nrbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,TRIM(file),"OLD","APPEND",32_i4)
          WRITE(con_unit) INT(ninterp*rb(ibloc)%mx,4),
     $      INT(ninterp*rb(ibloc)%my,4)
          CALL contour_laq2_write(rb(ibloc)%rz,ninterp,
     $                            one_record=.true.)
          CALL close_bin(con_unit,TRIM(file))
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
      DO ibl=nrbl_total+1,nbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,TRIM(file),"OLD","APPEND",32_i4)
          WRITE(con_unit) INT(tb(ibloc)%mvert,4),
     $      INT(tb(ibloc)%mcell,4)
          WRITE(con_unit) REAL(tb(ibloc)%tgeom%xs,4),
     $      REAL(tb(ibloc)%tgeom%ys,4)
          WRITE(con_unit) INT(tb(ibloc)%tgeom%vertex,4)
          CALL close_bin(con_unit,TRIM(file))
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE par_cont_init
c-----------------------------------------------------------------------
c     subprogram 19. equil_project
c     this subroutine projects flux-related computations to construct
c     the Lagrange-type expansions of the equilibrium fields.  these
c     operations had been part of gssolve. 
c-----------------------------------------------------------------------
      SUBROUTINE equil_project
      USE local
      USE fields
      USE input
      USE input_eq
      USE global
      USE time
      USE rblock
      USE tblock
      USE nimeq_ints
      USE computation_pointers
      USE pardata
      USE edge
      USE contour_mod
      USE iter_cg
      USE matrix_storage_mod
      USE plot_data
      USE boundary
      USE regularity
      USE nimeq_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ibe,jphi,its,iq,nmodes_save
      REAL(r8) :: err,err2,er0
      CHARACTER(8) :: seed
      CHARACTER(64) :: msg
      LOGICAL :: converged

      INTEGER(i4) :: ibasis,ix,iy,ipr,ibloc
      REAL(r8) :: dx,dy,b2,lam,p0
      REAL(r8), DIMENSION(3) :: vec
    
      TYPE(vector_type), DIMENSION(:), POINTER :: vtmp1,vtmp2
c-TEST
      REAL(r8) :: psinorm,dpsn
c-----------------------------------------------------------------------
c     interface block for the calc_pres routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_pres(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_pres
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the calc_f routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE calc_f(vect,flag)
        USE local
        USE fields
        USE vector_type_mod
        IMPLICIT NONE
        TYPE(vector_type), DIMENSION(:), POINTER :: vect
        CHARACTER(*), INTENT(IN) :: flag
        END SUBROUTINE calc_f
      END INTERFACE
c-----------------------------------------------------------------------
c     nimeq solves work with only 1 Fourier component (n=0).
c-----------------------------------------------------------------------
      nmodes_save=nmodes
      nmodes=1
c-----------------------------------------------------------------------
c     calculate B-poloidal from grad(phi)Xgrad(psi) and RB_phi
c     directly from F.
c-----------------------------------------------------------------------
      CALL calc_f(rwork3,"f only")

      ALLOCATE(vtmp1(nrbl),vtmp2(nrbl))
      DO ibl=1,nrbl
        CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,1_i4)
        CALL vector_type_alloc(vtmp1(ibl),poly_degree,rb(ibl)%mx,
     $                          rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vtmp2(ibl),poly_degree,rb(ibl)%mx,
     $                          rb(ibl)%my,3_i4)
        CALL rblock_get_rhs(rb(ibl),rhs(ibl),get_B,3_i4) 
        crhs(ibl)=0._r8
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                         3_i4,1_i4) 
        CALL vector_type_alloc(vtmp1(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vtmp2(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL tblock_get_rhs(tb(ibl),rhs(ibl),get_B,3_i4)
        crhs(ibl)=0._r8
      ENDDO
c-----------------------------------------------------------------------
c     eliminate cell centered data
c-----------------------------------------------------------------------
      IF (poly_degree>1) THEN
        DO ibl=1,nrbl
          vtmp1(ibl)=rhs(ibl)
        ENDDO
        CALL matelim_presolve(bfield_mat,vtmp1,vtmp2,3_i4)
        DO ibl=1,nrbl
          CALL cvector_assign_vec(crhs(ibl),vtmp2(ibl),'real',1_i4,3_i4)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'real',1_i4)
        ENDDO
      ELSE
        DO ibl=1,nbl
          CALL cvector_assign_vec(crhs(ibl),rhs(ibl),'real',1_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate rhs contributions across block borders and solve.
c-----------------------------------------------------------------------
      CALL edge_network(3_i4,1_i4,poly_degree-1_i4,.true.)

      DO ibl=1,SIZE(r0block_list)
        ibe=r0block_list(ibl)
        CALL regular_vec(crhs(ibe),seam(ibe),'cyl_vec',3_i4,
     $                   nmodes,nindex)
      ENDDO

      DO ibl=1,nbl
        CALL vector_assign_cvec(rhs(ibl),crhs(ibl),'REAL',1_i4)
        be_eq(ibl)=0
      END DO

      CALL iter_cg_2d_solve(bfield_mat,bfield_fac,be_eq,rhs,3_i4,tol,
     $                      linmaxit,nimeq_solver,err,its,seed)
      IF (err>tol) THEN
         WRITE(msg,'(a,i4,a,es10.3,3a)') 'Get B: no convergence: ',
     $         its,' its ',err,' err ',seed,' seed'
        CALL nim_stop(msg)
      ENDIF
      IF (poly_degree>1) THEN
        CALL matelim_real_postsolve(bfield_mat,be_eq,rhs,3_i4)
        DO ibl=1,nrbl
          be_eq(ibl)%arri=rhs(ibl)%arri
        ENDDO
      ENDIF

      DO ibl=1,nbl
        CALL vector_type_dealloc(rhs(ibl))
        CALL vector_type_dealloc(crhs(ibl))
        CALL vector_type_dealloc(vtmp1(ibl))
        CALL vector_type_dealloc(vtmp2(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     calculate J from B and psi.  note that
c
c       mu0*J=grad(F)xgrad(phi)+delstar(psi)*grad(phi)
c
c     so the poloidal part is just -(dF/dpsi)*B_eq.
c-----------------------------------------------------------------------
      CALL jeq_poloidal
c-----------------------------------------------------------------------
c     the toroidal part, J_phi/R is evaluated from -FF'/mu0*R**2 - P'.
c     free-boundary computations already have mu0*J_phi/R.
c-----------------------------------------------------------------------
      IF (gs_type/='free') THEN
        DO ibl=1,nrbl
          CALL vector_type_alloc(rhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,1_i4)
          CALL rblock_get_rhs(rb(ibl),rhs(ibl),jphi_rhs,1_i4) 
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_type_alloc(rhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
          CALL tblock_get_rhs(tb(ibl),rhs(ibl),jphi_rhs,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       network seams that connect block borders.
c-----------------------------------------------------------------------
        CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.true.)
c-----------------------------------------------------------------------
c       invert the mass matrix to project onto the finite element bases.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          sln(ibl)=0._r8
        ENDDO
        CALL iter_cg_2d_solve(mass_mat,mass_fac,sln,rhs,1_i4,tol,
     $                        linmaxit,nimeq_solver,err,its,seed)
        IF (err>tol) THEN
          WRITE(msg,'(a,i4,a,es10.3,3a)') 'Jphi: no convergence: ',
     $          its,' its ',err,' err ',seed,' seed'
          CALL nim_stop(msg)
        ENDIF
        DO ibl=1,nbl
          CALL vector_type_dealloc(rhs(ibl))
          CALL vector_assignq_vec(ja_eq(ibl),sln(ibl),nqty=1_i4,
     $                            nstart1=3_i4,nstart2=1_i4)
        ENDDO
      ENDIF
      DO ibl=1,nbl
        CALL vector_mult(ja_eq(ibl),1._r8/mu0)
      ENDDO
c-----------------------------------------------------------------------
c     convert stream function to poloidal flux, and evaluate pressure.
c-----------------------------------------------------------------------
      CALL calc_pres(pres_eq,"p only")
      DO ibl=1,nbl
        CALL vector_mult(pres_eq(ibl),1._r8/mu0)
        prese_eq(ibl)=pres_eq(ibl)
        CALL vector_mult(prese_eq(ibl),pe_frac)
      ENDDO

      DO ibl=1,nbl
        IF (geom=='tor') THEN
          CALL vector_mult(pflux(ibl),twopi)
        ELSE
          CALL vector_mult(pflux(ibl),per_length)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     alter the number density profile so that
c     p/p0=(n/ndens)**gamma_nimeq.
c-----------------------------------------------------------------------
      IF (gamma_nimeq>0) THEN
        IF (pres_model=="read1d") THEN
          p0=neqsq%fs(0,3)
        ELSE IF (pres_model=="quad_closed") THEN
          p0=p_open+p1_closed
        ELSE IF (pres_model=="cubic_closed") THEN
          p0=p_axis
        ELSE
          p0=pressure0
        ENDIF
        DO ibl=1,nbl
          nd_eq(ibl)%arr=ndens*
     $      (mu0*pres_eq(ibl)%arr/p0)**(1._r8/gamma_nimeq)
          IF (poly_degree==1) CYCLE
          nd_eq(ibl)%arrh=ndens*
     $      (mu0*pres_eq(ibl)%arrh/p0)**(1._r8/gamma_nimeq)
          nd_eq(ibl)%arrv=ndens*
     $      (mu0*pres_eq(ibl)%arrv/p0)**(1._r8/gamma_nimeq)
          nd_eq(ibl)%arri=ndens*
     $      (mu0*pres_eq(ibl)%arri/p0)**(1._r8/gamma_nimeq)
        ENDDO
c-TEST base n and diff_shape on the fllen function
      ELSE IF (gamma_nimeq<0.AND.btop_check/='none') THEN
        DO ibl=1,nbl
          DO iy=0,rb(ibl)%my
            DO ix=0,rb(ibl)%mx
              psinorm=rb(ibl)%fllen%fs(1,ix,iy)-rlencl
              rb(ibl)%diff_shape%fs(1,ix,iy)=1._r8+
     $          0.5_r8*dvac*(1._r8-TANH(dexp*psinorm))
              rb(ibl)%nd_eq%fs(1,ix,iy)=nedge+0.5_r8*(ndens-nedge)*
     $          (1._r8+TANH(dexp*psinorm))
            ENDDO
          ENDDO
          IF (poly_degree==1) CYCLE
          DO iy=0,rb(ibl)%my
            DO ix=1,rb(ibl)%mx
              DO ipr=1,SIZE(fllen(ibl)%arrh,2)
                psinorm=rb(ibl)%fllen%fsh(1,ipr,ix,iy)-rlencl
                rb(ibl)%diff_shape%fsh(1,ipr,ix,iy)=1._r8+
     $            0.5_r8*dvac*(1._r8-TANH(dexp*psinorm))
                rb(ibl)%nd_eq%fsh(1,ipr,ix,iy)=nedge+0.5_r8*
     $            (ndens-nedge)*(1._r8+TANH(dexp*psinorm))
              ENDDO
            ENDDO
          ENDDO
          DO iy=1,rb(ibl)%my
            DO ix=0,rb(ibl)%mx
              DO ipr=1,SIZE(fllen(ibl)%arrv,2)
                psinorm=rb(ibl)%fllen%fsv(1,ipr,ix,iy)-rlencl
                rb(ibl)%diff_shape%fsv(1,ipr,ix,iy)=1._r8+
     $            0.5_r8*dvac*(1._r8-TANH(dexp*psinorm))
                rb(ibl)%nd_eq%fsv(1,ipr,ix,iy)=nedge+0.5_r8*
     $            (ndens-nedge)*(1._r8+TANH(dexp*psinorm))
              ENDDO
            ENDDO
          ENDDO
          DO iy=1,rb(ibl)%my
            DO ix=1,rb(ibl)%mx
              DO ipr=1,SIZE(fllen(ibl)%arri,2)
                psinorm=rb(ibl)%fllen%fsi(1,ipr,ix,iy)-rlencl
                rb(ibl)%diff_shape%fsi(1,ipr,ix,iy)=1._r8+
     $            0.5_r8*dvac*(1._r8-TANH(dexp*psinorm))
                rb(ibl)%nd_eq%fsi(1,ipr,ix,iy)=nedge+0.5_r8*
     $            (ndens-nedge)*(1._r8+TANH(dexp*psinorm))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     ion diamagnetic flow as (1/(1+meomi)-pe_frac)*J_perp/(n_e*e)
c-----------------------------------------------------------------------
      IF (eq_flow=='diamagnetic') THEN
        DO ibl=1,nbl
          DO ibasis=1,SIZE(rb(ibl)%ve_eq%ix0)
            DO iy=rb(ibl)%ve_eq%iy0(ibasis),rb(ibl)%ve_eq%my
              DO ix=rb(ibl)%ve_eq%ix0(ibasis),rb(ibl)%ve_eq%mx
                dx=ix-rb(ibl)%ve_eq%ix0(ibasis)+rb(ibl)%ve_eq%dx(ibasis)
                dy=iy-rb(ibl)%ve_eq%iy0(ibasis)+rb(ibl)%ve_eq%dy(ibasis)
                CALL lagr_quad_eval(rb(ibl)%be_eq,dx,dy,0_i4)
                CALL lagr_quad_eval(rb(ibl)%ja_eq,dx,dy,0_i4)
                CALL lagr_quad_eval(rb(ibl)%nd_eq,dx,dy,0_i4)
                IF (geom=='tor') THEN
                  CALL lagr_quad_eval(rb(ibl)%rz,dx,dy,0_i4)
                  rb(ibl)%be_eq%f(3)=rb(ibl)%be_eq%f(3)/rb(ibl)%rz%f(1)
                  rb(ibl)%ja_eq%f(3)=rb(ibl)%ja_eq%f(3)*rb(ibl)%rz%f(1)
                ENDIF
                b2=SUM(rb(ibl)%be_eq%f**2)
                lam=SUM(rb(ibl)%be_eq%f*rb(ibl)%ja_eq%f)/b2
                vec=(1._r8/(1._r8+meomi)-pe_frac)*
     $              (rb(ibl)%ja_eq%f-lam*rb(ibl)%be_eq%f)/
     $              (rb(ibl)%nd_eq%f(1)*elementary_q)
                CALL lagr_quad_basis_assign_loc(
     $               rb(ibl)%ve_eq,vec,ibasis,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     deallocate temporary storage.
c-----------------------------------------------------------------------
      DEALLOCATE(vtmp1,vtmp2)
      nmodes=nmodes_save
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_project

c-----------------------------------------------------------------------
c     subprogram 20. equil_project
c     write the binary files for xdraw plots of the new equilibrium.
c     this had been part of gssolve. 
c-----------------------------------------------------------------------
      SUBROUTINE equil_plot(first_data,ndcon)
      USE local
      USE io
      USE fields
      USE pardata
      USE contour_mod
      USE mpi_nim
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: first_data
      INTEGER(i4), INTENT(IN) :: ndcon

      INTEGER(i4) :: ibl,ibloc,ierror
c-----------------------------------------------------------------------
c     the following write scheme assumes that all processors have
c     access to the same disk space.
c-----------------------------------------------------------------------
      IF (first_data) CALL par_cont_init("grad.bin",10_i4,ndcon)
      DO ibl=1,nrbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,"grad.bin","OLD","APPEND",32_i4)
          CALL contour_laq2_write(rb(ibloc)%be_eq,ndcon)
          CALL contour_laq2_write(rb(ibloc)%ja_eq,ndcon)
          CALL contour_laq2_write(rb(ibloc)%pflux,ndcon)
          CALL contour_laq2_write(rb(ibloc)%pres_eq,ndcon)
          CALL contour_laq2_write(rb(ibloc)%rwork1,ndcon)
          CALL contour_laq2_write(rb(ibloc)%fllen,ndcon)
          CALL close_bin(con_unit,"grad.bin")
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
      DO ibl=nrbl_total+1,nbl_total
        IF (block2proc(ibl)==node) THEN
          ibloc=global2local(ibl)
          CALL open_bin(con_unit,"grad.bin","OLD","APPEND",32_i4)
          CALL contour_tl2_write(tb(ibloc)%be_eq)
          CALL contour_tl2_write(tb(ibloc)%ja_eq)
          CALL contour_tl2_write(tb(ibloc)%pflux)
          CALL contour_tl2_write(tb(ibloc)%pres_eq)
          CALL contour_tl2_write(tb(ibloc)%rwork1)
          CALL contour_tl2_write(tb(ibloc)%fllen)
          CALL close_bin(con_unit,"grad.bin")
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE equil_plot

c-----------------------------------------------------------------------
c     subprogram 21. pflux_block
c     this external subroutine sets the pflux arrays based on a passed
c     data structure for lambda=psi/R**2.
c-----------------------------------------------------------------------
      SUBROUTINE pflux_block(lam,ibl)
      USE local
      USE fields
      USE input
      IMPLICIT NONE

      TYPE(vector_type), INTENT(IN) :: lam
      INTEGER(i4), INTENT(IN) :: ibl
c-----------------------------------------------------------------------
c     multiply by R**2 for toroidal geometry computations.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        pflux(ibl)%arr(1,:,:)=lam%arr(1,:,:)*
     $    rb(ibl)%rz%fs(1,:,:)**2
        IF (ASSOCIATED(lam%arrh).AND.ibl<=nrbl) THEN
          pflux(ibl)%arrh(1,:,:,:)=lam%arrh(1,:,:,:)*
     $      rb(ibl)%rz%fsh(1,:,:,:)**2
          pflux(ibl)%arrv(1,:,:,:)=lam%arrv(1,:,:,:)* 
     $      rb(ibl)%rz%fsv(1,:,:,:)**2
          pflux(ibl)%arri(1,:,:,:)=lam%arri(1,:,:,:)* 
     $      rb(ibl)%rz%fsi(1,:,:,:)**2
        ENDIF
      ELSE
        CALL vector_assignq_vec(pflux(ibl),lam,nqty=1_i4,
     $                          nstart1=1_i4,nstart2=1_i4)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pflux_block

c-----------------------------------------------------------------------
c     subprogram 22. update_qfp
c     update the qrwork1 and qpres_eq quad-point storage of F or F'
c     and P or P'.
c-----------------------------------------------------------------------
      SUBROUTINE update_qfp
      USE local
      USE fields
      USE rblock
      USE tblock
      IMPLICIT NONE

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     interpolate the flux functions to the quadrature-point
c     storage for reuse.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL rblock_qp_update(rb(ibl)%rwork3,rb(ibl)%qrwork1,rb(ibl))
        CALL rblock_qp_update(rb(ibl)%pres_eq,rb(ibl)%qpres_eq,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_qp_update(tb(ibl)%rwork3,tb(ibl)%qrwork1,tb(ibl))
        CALL tblock_qp_update(tb(ibl)%pres_eq,tb(ibl)%qpres_eq,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE update_qfp
c-----------------------------------------------------------------------
c     subprogram 24. advect_btop.
c     this routine calls the passive-advection PDE solver twice, with
c     the "flow" run in opposite directions, then averages the results.
c
c-PRE
c     this symmetrization step should help improve the accuracy of
c     the solution, but it doesn't, so it is not used at present.
c-----------------------------------------------------------------------
      SUBROUTINE advect_btop(bcsym,refln)
      USE local
      USE computation_pointers
      USE global
      USE fields
      USE input
      IMPLICIT NONE

      CHARACTER(1), INTENT(IN) :: bcsym
      REAL(r8), INTENT(IN) :: refln

      INTEGER(i4) :: ibl
      TYPE(vector_type), DIMENSION(:), POINTER, SAVE :: flfwd
      LOGICAL, SAVE :: first_call=.true.
c-----------------------------------------------------------------------
c     create space to save the forward computation.
c-----------------------------------------------------------------------
      IF (first_call) THEN
        ALLOCATE(flfwd(nbl))
        DO ibl=1,nrbl
          CALL vector_type_alloc(flfwd(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,1_i4)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_type_alloc(flfwd(ibl),1_i4,tb(ibl)%mvert,
     $                           0_i4,1_i4)
        ENDDO
        first_call=.false.
      ENDIF
c-----------------------------------------------------------------------
c     call advect_solve advecting parallel to B_pol, then save the
c     result in rhs, which is not used by advect_solve.
c-----------------------------------------------------------------------
      integrand_flag="forward"
      CALL advect_solve(bcsym,refln)
      DO ibl=1,nbl
        flfwd(ibl)=fllen(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     call advect_solve advecting antiparallel to B_pol, then add the
c     forward result.
c-----------------------------------------------------------------------
      integrand_flag="reverse"
      CALL advect_solve(bcsym,refln)
      DO ibl=1,nbl
        CALL vector_add(fllen(ibl),flfwd(ibl),v1fac=0.5_r8,v2fac=0.5_r8)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advect_btop
c-----------------------------------------------------------------------
c     subprogram 24. advect_solve.
c     solve a passive advection-like equation using least-squares
c     projection to distinguish open and closed flux.
c-----------------------------------------------------------------------
      SUBROUTINE advect_solve(bcsym,refln)
      USE local
      USE input
      USE physdat
      USE nimeq_ints
      USE nimeq_mod
      USE boundary
      USE matrix_storage_mod
      USE computation_pointers
      USE matrix_mod
      USE rblock
      USE tblock
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE pardata
      IMPLICIT NONE

      CHARACTER(1), INTENT(IN) :: bcsym
      REAL(r8), INTENT(IN) :: refln

      INTEGER(i4) :: ibl,its,iv,ibe
      REAL(r8) :: err,dscale
      CHARACTER(8) :: seed
      CHARACTER(128) :: msg
      TYPE(vector_type), DIMENSION(:), POINTER, SAVE :: vtmp1,vtmp2
      LOGICAL, SAVE :: first_call=.true.
c-----------------------------------------------------------------------
c     estimate the magnitude of B_pol based on flux values and the
c     reference length, and make it accessible to the integrand
c     routines via the global module.
c-----------------------------------------------------------------------
      flux_n0=(MAX(psic_max,psio_max)-MIN(psic_min,psio_min))/refln**2
c-----------------------------------------------------------------------
c     create the advection operator.  with full matrix storage,
c     this needs to be updated each step.
c-PRE matrix-free could avoid some of this.
c-----------------------------------------------------------------------
      dscale=0._r8
      DO ibl=1,nrbl
        CALL rblock_make_matrix(rb(ibl),fladv_mat%rbl_mat(ibl),
     $                          scal_adv_op,1_i4)
        dscale=MAX(dscale,MAXVAL(ABS(fladv_mat%rbl_mat(ibl)%mat(1,1)%
     $                               arr(1,0,0,1,:,:))))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_make_matrix(tb(ibl),fladv_mat%tbl_mat(ibl)%lmat,
     $                          scal_adv_op,1_i4)
        DO iv=0,SIZE(fladv_mat%tbl_mat(ibl)%lmat)-1
          dscale=MAX(dscale,ABS(fladv_mat%tbl_mat(ibl)%lmat(iv)%
     $                          element(1,1,0)))
        ENDDO
      ENDDO
      fladv_mat%diag_scale=dscale

      CALL dirichlet_op(fladv_mat,"all",dscale)
      CALL matelim_real_inv_int(fladv_mat,1_i4)
      DO ibl=1,nrbl
        IF (rb(ibl)%degenerate)
     $    CALL matrix_degen_collect_real
     $        (fladv_mat%rbl_mat(ibl),1_i4,
     $         fladv_mat%symmetric)
      ENDDO
      CALL iter_factor(fladv_mat,fladv_fac,1_i4,nimeq_solver,
     $                 off_diag_fac)
c-----------------------------------------------------------------------
c     create space for algebra if not already done.
c-----------------------------------------------------------------------
      IF (first_call) THEN
        ALLOCATE(vtmp1(nbl),vtmp2(nrbl))
        DO ibl=1,nrbl
          CALL vector_type_alloc(vtmp1(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,1_i4)
          CALL vector_type_alloc(vtmp2(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,1_i4)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_type_alloc(vtmp1(ibl),1_i4,tb(ibl)%mvert,
     $                           0_i4,1_i4)
        ENDDO
        first_call=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL rblock_get_rhs(rb(ibl),vtmp1(ibl),scal_adv_rhs,1_i4)
        vtmp2(ibl)=vtmp1(ibl)
        fllen(ibl)=0._r8
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_get_rhs(tb(ibl),vtmp1(ibl),scal_adv_rhs,1_i4)
        fllen(ibl)=0._r8
      ENDDO
      IF (poly_degree>1) THEN
        CALL matelim_presolve(fladv_mat,vtmp2,vtmp1,1_i4)
      ENDIF
      DO ibl=1,nbl
        CALL edge_load_arr(vtmp1(ibl),1_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     network block seams and apply boundary conditions
c-----------------------------------------------------------------------
      CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.false.)
      DO ibl=1,nbl
        CALL edge_unload_arr(vtmp1(ibl),1_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     apply homogeneous essential conditions along the exterior
c     border, so that advection zeroes-out fllen in the open field.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        crhs(ibe)=0._r8
        CALL cvector_assign_vec(crhs(ibe),vtmp1(ibe),'REAL',
     $                          1_i4,nqty=1_i4)
        CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4,bcsym)
        CALL vector_assign_cvec(vtmp1(ibe),crhs(ibe),'REAL',
     $                          1_i4,nqty=1_i4)
      ENDDO
c-----------------------------------------------------------------------
c     solve the algebraic system.  the presumably weaker gsh_tol should
c     be sufficient.
c-----------------------------------------------------------------------
      CALL iter_cg_2d_solve(fladv_mat,fladv_fac,fllen,vtmp1,1_i4,
     $                      gsh_tol,linmaxit,nimeq_solver,err,its,seed)
      IF (err>gsh_tol) THEN
        WRITE(msg,'(2a,i4,a,es10.3,3a)') 'Advect_btop: ',
     $    'no convergence: ',its,' its ',err,' err ',seed,' seed'
        CALL nim_stop(msg)
      ENDIF
c-----------------------------------------------------------------------
c     complete the solution at element interiors.
c-----------------------------------------------------------------------
      IF (poly_degree>1) THEN
        CALL matelim_real_postsolve(fladv_mat,fllen,vtmp1,1_i4)
        DO ibl=1,nrbl
          fllen(ibl)%arri=vtmp1(ibl)%arri
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advect_solve
c-----------------------------------------------------------------------
c     subprogram 25. pfsmth_init.
c     this subroutine sets the coefficients for the Hermite cubic
c     expansions used to smooth P and F near the last closed flux
c     surface.
c-----------------------------------------------------------------------
      SUBROUTINE pfsmth_init
      USE local
      USE nimeq_mod
      USE input_eq
      IMPLICIT NONE

      REAL(r8), EXTERNAL :: p_func,f_func
c-----------------------------------------------------------------------
c     evaluate the P and F functions at the start of the respective
c     smoothing regions and in the open flux.  also evaluate their
c     derivatives with respect to normalized flux at the start of the
c     regions.
c-----------------------------------------------------------------------
      pvalin=p_func(psin_smthpp,0_i4,2._r8*rlencl,.true.)
      pdrvin=p_func(psin_smthpp,1_i4,2._r8*rlencl,.true.)
      pvalout=p_func(1.1_r8,0_i4,0.1_r8*rlencl,.true.)
      pdpsi=1._r8-psin_smthpp
      psmthset=.true.

      fvalin=f_func(psin_smthfp,0_i4,2._r8*rlencl,.true.)
      fdrvin=f_func(psin_smthfp,1_i4,2._r8*rlencl,.true.)
      fvalout=f_func(1.1_r8,0_i4,0.1_r8*rlencl,.true.)
      fdpsi=1._r8-psin_smthfp
      fsmthset=.true.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pfsmth_init
