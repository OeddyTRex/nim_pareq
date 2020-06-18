c-----------------------------------------------------------------------
c     file iter_gmres_r2d.f
c
c     this module contains routines for GMRES solves for real non-
c     symmetric systems.  if the mfsolve routine is used, the matrix-
c     vector products are computed with an externally supplied
c     subroutine, analogous to the complex 3d solves.
c
c     this module is based on iter_dir_nonsym, which is used for complex
c     2d solves.
c
c     The GMRES algorithm is from "Iterative Methods for Sparse Linear
c     Systems," 2nd ed., Y. Saad, SIAM (2003).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_gmres_r2d
c     1. iter_gmr_r2d_solve.
c     2. iter_gmr_r2d_mfsolve.
c     3. iter_gmr_r2d_init.
c     4. iter_gmr_r2d_dealloc.
c     5. iter_r2d_ortho.
c     6. iter_r2d_grot.
c     7. iter_r2d_min.
c     8. iter_gmr_r2d_mfpre.
c-----------------------------------------------------------------------
c     module containing the conjugate gradient precon.
c-----------------------------------------------------------------------
      MODULE iter_gmres_r2d
      USE local
      USE edge_type_mod
      USE matrix_type_mod
      USE vector_type_mod
      USE matrix_mod
      USE edge
      USE seam_storage_mod
      USE iter_utils
      USE iter_real_direct
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     module-specific definitions.
c-----------------------------------------------------------------------
      INTEGER(i4), PRIVATE :: nrb,ntotb,poly_deg,poly_d2
      REAL(r8), PRIVATE, PARAMETER :: one=1._r8
      LOGICAL, PRIVATE :: use_int
      LOGICAL, PARAMETER, PRIVATE :: writecheck=.false.

      TYPE :: gm_basis_type
        TYPE(vector_type), POINTER, DIMENSION(:) :: base
      END TYPE gm_basis_type

      INTEGER(i4) :: n_gm
      REAL(r8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: hess,givn
      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: errvec,minvec

      TYPE(gm_basis_type), DIMENSION(:), POINTER, PRIVATE :: ep_d,pinv

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iter_gmr_r2d_solve.
c     manages the gmres computation using the supplied matrix data
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_r2d_solve(mat_str,fac_str,sln,rhs,nqty,tol,
     $                              maxit,precon,err,its,seed)
      USE matrix_mod
      USE factor_type_mod
      USE mpi_nim
      USE pardata
      USE time

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(matrix_factor_type), INTENT(IN) :: fac_str
      TYPE(vector_type), DIMENSION(:), POINTER :: sln,rhs
      INTEGER(i4), INTENT(IN) :: nqty,maxit
      REAL(r8), INTENT(IN) :: tol
      CHARACTER(*), INTENT(IN) :: precon
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(vector_type), DIMENSION(:), POINTER :: res
      INTEGER(i4) :: ibl,it
      INTEGER(i4) :: ix,iy,iv,mv,iq
      REAL(r8) :: er0
      CHARACTER(64) :: msg

      REAL(r8) :: tmp,timestart_dir,timeend_dir
      INTEGER(i4) :: ierror,maxloop
      LOGICAL, SAVE :: first_call=.true.
      LOGICAL, PARAMETER :: residcheck=.true.
      LOGICAL, PARAMETER :: everycheck=.false.
c-----------------------------------------------------------------------
c     interface block for 2D preconditioner driver.
c-----------------------------------------------------------------------
      INCLUDE "iter_precon_intf.inc"
c-----------------------------------------------------------------------
c     verify that preconditioner choice is available.
c-----------------------------------------------------------------------
      IF (.NOT.(precon=='seq_slu'.OR.precon(1:5)=='slu_d'.OR.
     $          precon=='diagonal'.OR.precon=='no prec'))
     $  CALL nim_stop
     $  ('iter_gmr_r2d_solve: nonsymmetric precon '//
     $   'choice is not available.')
c-----------------------------------------------------------------------
c     allocate storage; determine best guess and initial error.
c-----------------------------------------------------------------------
      CALL timer(timestart_dir)
      IF (first_call) THEN
        ALLOCATE(hess(maxit+1,maxit))
        ALLOCATE(errvec(maxit+1))
        ALLOCATE(minvec(maxit))
        ALLOCATE(givn(2,maxit+1))
        ALLOCATE(ep_d(maxit+1))
        ALLOCATE(pinv(SIZE(ep_d)))
        first_call=.false.
      ENDIF
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      CALL iter_gmr_r2d_init(mat_str,sln,rhs,nqty,er0,
     $                       tol,res,.NOT.mat_str%eliminated)
c-----------------------------------------------------------------------
c     pre-step for GMRES.
c-----------------------------------------------------------------------
      CALL iter_resid(mat_str,sln,res,rhs,nqty,ntotb,nrb,poly_deg)
      CALL iter_2_norm(res,err,1._r8,ntotb,nrb,poly_deg,use_int)
      seed="guess"
      IF (err>er0) THEN
        DO ibl=1,ntotb
          res(ibl)=rhs(ibl)
          sln(ibl)=0._r8
        ENDDO
        err=er0
        seed="0 vec"
      ENDIF
      errvec(1)=err
      IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,es15.8)')
     $   '   Initial GMR is ',err

      IF (err/er0>0._r8) THEN
        DO ibl=1,ntotb
          CALL vector_mult(res(ibl),1._r8/err)
          ep_d(1)%base(ibl)=res(ibl)
        ENDDO
      ENDIF
      IF (ABS(errvec(1))/er0<=tol) THEN
        n_gm=0
        maxloop=0
      ELSE
        maxloop=maxit
      ENDIF
c-----------------------------------------------------------------------
c     begin iterations.
c-----------------------------------------------------------------------
      itloop: DO it=1,maxloop
c-----------------------------------------------------------------------
c       apply right preconditioning.  note that res is just
c       used as temporary storage.  also save P^-1 dotted into the
c       current basis vector for more reliable finite-precision
c       math.
c-----------------------------------------------------------------------
        CALL iter_pre_real(mat_str,fac_str,ep_d(it)%base,pinv(it)%base,
     $                     res,nqty,nrb,ntotb,precon)
c-----------------------------------------------------------------------
c       compute the inner product with the matrix.
c-----------------------------------------------------------------------
        CALL matvec(mat_str,pinv(it)%base,ep_d(it+1)%base,nqty)
c-----------------------------------------------------------------------
c       orthogonalize and add another column to the Hessenberg matrix.
c-----------------------------------------------------------------------
        CALL iter_r2d_ortho(it)
c-----------------------------------------------------------------------
c       find the Givens rotation for the last Hessenberg entry
c-----------------------------------------------------------------------
        CALL iter_r2d_grot(it)
c-----------------------------------------------------------------------
c       check tolerance.
c-----------------------------------------------------------------------
        n_gm=it
        IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,i4,a,es15.8)')
     $     '   GMR at step',it,' is ',ABS(errvec(it+1))
        IF (ABS(errvec(it+1))/er0<=tol) EXIT
      ENDDO itloop
c-----------------------------------------------------------------------
c     determine the minimizing vector of coefficients, and 
c     construct the solution with the saved (P^-1 . basis) vectors.
c-----------------------------------------------------------------------
      IF (n_gm>=1) THEN
        CALL iter_r2d_min
        DO it=1,n_gm
          DO ibl=1,ntotb
            CALL vector_add(sln(ibl),pinv(it)%base(ibl),
     $                      one,minvec(it))
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     check tolerance.
c-----------------------------------------------------------------------
      IF (residcheck) THEN
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,*) 'Final iterate  residual=',ABS(errvec(n_gm+1))
        CALL iter_resid(mat_str,sln,res,rhs,nqty,ntotb,nrb,poly_deg)
        CALL iter_2_norm(res,err,1._r8,ntotb,nrb,poly_deg,use_int)
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,*) 'Final separate residual=',err
      ENDIF

      err=err/er0
      its=n_gm

      CALL iter_gmr_r2d_dealloc(res)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(comm_layer,ierror)
      CALL timer(timeend_dir)
      time_iter = time_iter + timeend_dir-timestart_dir

      RETURN
      END SUBROUTINE iter_gmr_r2d_solve
c-----------------------------------------------------------------------
c     subprogram 2. iter_gmr_r2d_mfsolve.
c     manages the "matrix-free" gmres computation using the external
c     subroutine to compute matrix-vector products.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_r2d_mfsolve(dot_routine,mat_str,fac_str,sln,
     $                                rhs,nqty,tol,maxit,precon,err,its,
     $                                seed)
      USE matrix_mod
      USE factor_type_mod
      USE mpi_nim
      USE pardata
      USE time

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(matrix_factor_type), INTENT(IN) :: fac_str
      TYPE(vector_type), DIMENSION(:), POINTER :: sln,rhs
      INTEGER(i4), INTENT(IN) :: nqty,maxit
      REAL(r8), INTENT(IN) :: tol
      CHARACTER(*), INTENT(IN) :: precon
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(vector_type), DIMENSION(:), POINTER :: res,rtmp
      INTEGER(i4) :: ibl,it
      INTEGER(i4) :: ix,iy,iv,mv,iq
      REAL(r8) :: er0
      CHARACTER(64) :: msg

      REAL(r8) :: tmp,timestart_dir,timeend_dir
      INTEGER(i4) :: ierror,maxloop
      LOGICAL, SAVE :: first_call=.true.
      LOGICAL, PARAMETER :: residcheck=.true.
      LOGICAL, PARAMETER :: everycheck=.false.
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.  the format is
c     the same as those used for 3d solves, except that the data
c     structures are for real data.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dot_routine(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE dot_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     verify that preconditioner choice is available.
c-----------------------------------------------------------------------
      IF (.NOT.(precon=='seq_slu'.OR.precon(1:5)=='slu_d'.OR.
     $          precon=='diagonal'.OR.precon=='no prec'))
     $  CALL nim_stop('iter_gmr_r2d_mfsolve: '//
     $    'nonsymmetric precon choice is not available.')
c-----------------------------------------------------------------------
c     allocate storage; determine best guess and initial error.
c     when using matrix-free computation, interior nodes are used
c     in the GMRES iteration, hence the .true. specification for
c     the init routine's intinfo.
c-----------------------------------------------------------------------
      CALL timer(timestart_dir)
      IF (first_call) THEN
        ALLOCATE(hess(maxit+1,maxit))
        ALLOCATE(errvec(maxit+1))
        ALLOCATE(minvec(maxit))
        ALLOCATE(givn(2,maxit+1))
        ALLOCATE(ep_d(maxit+1))
        ALLOCATE(pinv(SIZE(ep_d)))
        first_call=.false.
      ENDIF
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      CALL iter_gmr_r2d_init(mat_str,sln,rhs,nqty,er0,
     $                       tol,res,.true.,rtmp)
c-----------------------------------------------------------------------
c     compute the residual and its norm.
c-----------------------------------------------------------------------
      CALL dot_routine(sln,res,.true.)
      DO ibl=1,ntotb
        CALL vector_add(res(ibl),rhs(ibl),v1fac=-1._r8)
      ENDDO
      CALL iter_2_norm(res,err,1._r8,ntotb,nrb,poly_deg,use_int)
c-----------------------------------------------------------------------
c     pre-step for GMRES.
c-----------------------------------------------------------------------
      seed="guess"
      IF (err>er0) THEN
        DO ibl=1,ntotb
          res(ibl)=rhs(ibl)
          sln(ibl)=0._r8
        ENDDO
        err=er0
        seed="0 vec"
      ENDIF
      errvec(1)=err
      IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,es15.8)')
     $   '   Initial GMR is ',err

      IF (err/er0>0._r8) THEN
        DO ibl=1,ntotb
          CALL vector_mult(res(ibl),1._r8/err)
          ep_d(1)%base(ibl)=res(ibl)
        ENDDO
      ENDIF
      IF (ABS(errvec(1))/er0<=tol) THEN
        n_gm=0
        maxloop=0
      ELSE
        maxloop=maxit
      ENDIF
c-----------------------------------------------------------------------
c     begin iterations.
c-----------------------------------------------------------------------
      itloop: DO it=1,maxloop
c-----------------------------------------------------------------------
c       apply right preconditioning.  note that res is just
c       used as temporary storage.  also save P^-1 dotted into the
c       current basis vector for more reliable finite-precision
c       math.
c-----------------------------------------------------------------------
        CALL iter_gmr_r2d_mfpre(mat_str,fac_str,it,res,rtmp,nqty,nrb,
     $                          ntotb,precon)
c-----------------------------------------------------------------------
c       compute the inner product with the matrix.
c-----------------------------------------------------------------------
        CALL dot_routine(pinv(it)%base,ep_d(it+1)%base,.false.)
c-----------------------------------------------------------------------
c       orthogonalize and add another column to the Hessenberg matrix.
c-----------------------------------------------------------------------
        CALL iter_r2d_ortho(it)
c-----------------------------------------------------------------------
c       find the Givens rotation for the last Hessenberg entry
c-----------------------------------------------------------------------
        CALL iter_r2d_grot(it)
c-----------------------------------------------------------------------
c       check tolerance.
c-----------------------------------------------------------------------
        n_gm=it
        IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,i4,a,es15.8)')
     $     '   GMR at step',it,' is ',ABS(errvec(it+1))
        IF (ABS(errvec(it+1))/er0<=tol) EXIT
      ENDDO itloop
c-----------------------------------------------------------------------
c     determine the minimizing vector of coefficients, and 
c     construct the solution with the saved (P^-1 . basis) vectors.
c-----------------------------------------------------------------------
      IF (n_gm>=1) THEN
        CALL iter_r2d_min
        DO it=1,n_gm
          DO ibl=1,ntotb
            CALL vector_add(sln(ibl),pinv(it)%base(ibl),
     $                      one,minvec(it))
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     check tolerance.
c-----------------------------------------------------------------------
      IF (residcheck) THEN
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,*) 'Final iterate  residual=',ABS(errvec(n_gm+1))
        CALL dot_routine(sln,res,.false.)
        DO ibl=1,ntotb
          CALL vector_add(res(ibl),rhs(ibl),v1fac=-1._r8)
        ENDDO
        CALL iter_2_norm(res,err,1._r8,ntotb,nrb,poly_deg,use_int)
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,*) 'Final separate residual=',err
      ENDIF

      err=err/er0
      its=n_gm

      CALL iter_gmr_r2d_dealloc(res,rtmp)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(comm_layer,ierror)
      CALL timer(timeend_dir)
      time_iter = time_iter + timeend_dir-timestart_dir

      RETURN
      END SUBROUTINE iter_gmr_r2d_mfsolve
c-----------------------------------------------------------------------
c     subprogram 3. iter_gmr_r2d_init.
c     allocates storage used during the iterative solve.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_r2d_init(mat,sln,rhs,nqty,er0,tol,
     $                             res,intinfo,rextra)
      USE mpi_nim
      USE pardata

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(vector_type), DIMENSION(:), POINTER :: sln,rhs
      TYPE(vector_type), DIMENSION(:), POINTER :: res
      TYPE(vector_type), DIMENSION(:), POINTER, OPTIONAL :: rextra
      INTEGER(i4), INTENT(IN) :: nqty
      REAL(r8), INTENT(OUT) :: er0
      REAL(r8), INTENT(IN) :: tol
      LOGICAL, INTENT(IN) :: intinfo

      INTEGER(i4) :: ibl,mx,my,ierror,igm,iv,ix,iy,iq,ibase
      REAL(r8) :: er_0v
c-----------------------------------------------------------------------
c     determine the number of basis functions.
c-----------------------------------------------------------------------
      IF (nrb>0) THEN
        IF (mat%rbl_mat(1)%nbtype>1) THEN
          poly_deg=SUM(mat%rbl_mat(1)%nb_type(1:2))
        ELSE
          poly_deg=1
        ENDIF
      ELSE
c-PRE
        poly_deg=1
      ENDIF
      poly_d2=poly_deg**2
c-----------------------------------------------------------------------
c     set a module flag that indicates whether the matrix includes
c     data at cell interiors, as relevant to the calling solver routine.
c-----------------------------------------------------------------------
      use_int=intinfo
c-----------------------------------------------------------------------
c     ensure that the guess has exactly the same value in different
c     images of the same node along block borders.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL edge_load_arr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_in(1:nqty)=
     $      seam(ibl)%vertex(iv)%seam_in(1:nqty)*
     $      seam(ibl)%vertex(iv)%ave_factor
          IF (poly_deg>1) THEN
            seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))=
     $        seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))*
     $        seam(ibl)%segment(iv)%ave_factor
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(nqty,0_i4,poly_deg-1_i4,.false.)
      DO ibl=1,ntotb
        CALL edge_unload_arr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     allocate storage for the gmres operations.
c-----------------------------------------------------------------------
      ALLOCATE(res(ntotb))
      IF (PRESENT(rextra)) ALLOCATE(rextra(ntotb))
      DO ibl=1,ntotb
        mx=SIZE(rhs(ibl)%arr,2)-1
        my=SIZE(rhs(ibl)%arr,3)-1
        IF (SIZE(rhs(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_r2d_init: rhs dimension inconsistent with nqty')
        IF (SIZE(sln(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_r2d_init: sln dimension inconsistent with nqty')
        CALL vector_type_alloc(res(ibl),poly_deg,mx,my,nqty)
        IF (PRESENT(rextra))
     $    CALL vector_type_alloc(rextra(ibl),poly_deg,mx,my,nqty)
      ENDDO
      DO igm=1,SIZE(ep_d)
        ALLOCATE(ep_d(igm)%base(ntotb))
        ALLOCATE(pinv(igm)%base(ntotb))
        DO ibl=1,ntotb
          mx=SIZE(rhs(ibl)%arr,2)-1
          my=SIZE(rhs(ibl)%arr,3)-1
          CALL vector_type_alloc(ep_d(igm)%base(ibl),
     $                           poly_deg,mx,my,nqty)
          CALL vector_type_alloc(pinv(igm)%base(ibl),
     $                           poly_deg,mx,my,nqty)
        ENDDO
      ENDDO
      CALL iter_2_norm(rhs,er_0v,1._r8,ntotb,nrb,poly_deg,use_int)
      er0=MAX(er_0v,SQRT(TINY(er0)))
      hess=0._r8
      errvec=0._r8
      minvec=0._r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gmr_r2d_init
c-----------------------------------------------------------------------
c     subprogram 4. iter_gmr_r2d_dealloc.
c     deallocates storage used during the solve.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_r2d_dealloc(res,rextra)

      INTEGER(i4) :: ibl,igm
      TYPE(vector_type), DIMENSION(:), POINTER :: res
      TYPE(vector_type), DIMENSION(:), POINTER, OPTIONAL :: rextra
c-----------------------------------------------------------------------
c     deallocate storage.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL vector_type_dealloc(res(ibl))
        IF (PRESENT(rextra)) CALL vector_type_dealloc(rextra(ibl))
      ENDDO
      DEALLOCATE(res)
      IF (PRESENT(rextra)) DEALLOCATE(rextra)
      DO igm=1,SIZE(ep_d)
        DO ibl=1,ntotb
          CALL vector_type_dealloc(ep_d(igm)%base(ibl))
          CALL vector_type_dealloc(pinv(igm)%base(ibl))
        ENDDO
        DEALLOCATE(ep_d(igm)%base)
        DEALLOCATE(pinv(igm)%base)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gmr_r2d_dealloc
c-----------------------------------------------------------------------
c     subprogram 5. iter_r2d_ortho.
c     orthogonalize with modified Gram-Schmidt.
c-----------------------------------------------------------------------
      SUBROUTINE iter_r2d_ortho(jout)
      USE pardata

      INTEGER(i4), INTENT(IN) :: jout

      INTEGER(i4) :: ibl,ii,ierror
      REAL(r8) :: rdot,rsc
      REAL(r8) :: rmag,rmag0
      REAL(r8), PARAMETER :: small=1.e-14,orthtol=1.e-2
c-----------------------------------------------------------------------
c     compute the norm of A.v_jout prior to orthogonalization.
c-----------------------------------------------------------------------
      CALL iter_2_norm(ep_d(jout+1)%base,rmag0,1._r8,ntotb,nrb,
     $                 poly_deg,use_int)
c-----------------------------------------------------------------------
c     loop over existing basis vectors and save the inner products
c     in the hess array.
c-----------------------------------------------------------------------
      DO ii=1,jout
        rdot=0._r8
        DO ibl=1,ntotb
          CALL iter_dot(rdot,ep_d(ii)%base(ibl),ep_d(jout+1)%base(ibl),
     $                  seam(ibl),ibl,nrb,poly_deg,use_int)
        ENDDO
        IF (nprocs_layer>1) THEN
          CALL mpi_allreduce(rdot,hess(ii,jout),1,mpi_nim_real,
     $         mpi_sum,comm_layer,ierror)
        ELSE
          hess(ii,jout)=rdot
        ENDIF
        DO ibl=1,ntotb
          CALL vector_add(ep_d(jout+1)%base(ibl),ep_d(ii)%base(ibl),
     $                    one,-hess(ii,jout))
        ENDDO
      ENDDO
      CALL iter_2_norm(ep_d(jout+1)%base,rmag,1._r8,ntotb,nrb,
     $                 poly_deg,use_int)
c-----------------------------------------------------------------------
c     repeat orthogonalization of this vector if its magnitude is
c     reduced by more than a factor of orthtol.
c-----------------------------------------------------------------------
      IF (rmag/rmag0<orthtol) THEN
        DO ii=1,jout
          rdot=0._r8
          DO ibl=1,ntotb
            CALL iter_dot(rdot,ep_d(ii)%base(ibl),
     $                    ep_d(jout+1)%base(ibl),
     $                    seam(ibl),ibl,nrb,poly_deg,use_int)
          ENDDO
          IF (nprocs_layer>1) THEN
            CALL mpi_allreduce(rdot,rsc,1,mpi_nim_real,
     $           mpi_sum,comm_layer,ierror)
            rdot=rsc
          ENDIF
          hess(ii,jout)=hess(ii,jout)+rdot
          DO ibl=1,ntotb
            CALL vector_add(ep_d(jout+1)%base(ibl),ep_d(ii)%base(ibl),
     $                      one,-rdot)
          ENDDO
        ENDDO
        CALL iter_2_norm(ep_d(jout+1)%base,rmag,1._r8,ntotb,nrb,
     $                   poly_deg,use_int)
      ENDIF
c-----------------------------------------------------------------------
c     normalize.
c-----------------------------------------------------------------------
      IF (rmag<small*MAXVAL(ABS(hess(1:jout,jout)))) THEN
        hess(jout+1,jout)=0._r8
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,i4)') '   No new vector at step ',jout
      ELSE
        hess(jout+1,jout)=rmag
        DO ibl=1,ntotb
          CALL vector_mult(ep_d(jout+1)%base(ibl),1._r8/rmag)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_r2d_ortho
c-----------------------------------------------------------------------
c     subprogram 6. iter_r2d_grot.
c     find and apply the Givens rotation to remove below-diagonal
c     entries from the Hessenberg matrix and update the error vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_r2d_grot(jout)

      INTEGER(i4), INTENT(IN) :: jout

      INTEGER(i4) :: ii,jj
      REAL(r8) :: gden
      REAL(r8) :: gs,gc,rtmp
c-----------------------------------------------------------------------
c     Apply existing rotations to the new column.
c-----------------------------------------------------------------------
      DO ii=1,jout-1
        rtmp=givn(1,ii)*hess(ii,jout)+
     $       givn(2,ii)*hess(ii+1,jout)
        hess(ii+1,jout)=givn(1,ii)*hess(ii+1,jout)-
     $                  givn(2,ii)*hess(ii,jout)
        hess(ii,jout)=rtmp
      ENDDO
c-----------------------------------------------------------------------
c     Compute the new Givens sine and cosine coefficients.
c-----------------------------------------------------------------------
      gden=hess(jout,jout)*hess(jout,jout)+
     $     hess(jout+1,jout)*hess(jout+1,jout)
      gden=SQRT(gden)
      givn(1,jout)=hess(jout,jout)/gden
      givn(2,jout)=hess(jout+1,jout)/gden
c-----------------------------------------------------------------------
c     Apply the new rotation to the new Hessenberg column and the
c     error vector.
c-----------------------------------------------------------------------
      gc=givn(1,jout)
      gs=givn(2,jout)
      rtmp=gc*hess(jout,jout)+gs*hess(jout+1,jout)
      hess(jout+1,jout)=givn(1,jout)*hess(jout+1,jout)-
     $                  givn(2,jout)*hess(jout,jout)
      hess(jout,jout)=rtmp
      errvec(jout+1)=-givn(2,jout)*errvec(jout)
      errvec(jout)=gc*errvec(jout)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_r2d_grot
c-----------------------------------------------------------------------
c     subprogram 7. iter_r2d_min.
c     find the minimizing vector of coefficients for the orthogonal
c     basis vectors.
c-----------------------------------------------------------------------
      SUBROUTINE iter_r2d_min

      INTEGER(i4) :: ii,jj

      DO ii=n_gm,1,-1
        minvec(ii)=errvec(ii)
        DO jj=ii+1,n_gm
          minvec(ii)=minvec(ii)-hess(ii,jj)*minvec(jj)
        ENDDO
        minvec(ii)=minvec(ii)/hess(ii,ii)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_r2d_min
c-----------------------------------------------------------------------
c     subprogram 8. iter_gmr_r2d_mfpre.
c     apply preconditioning for the matrix-free solves, which may need
c     a static-condensation step prior to calling the standard
c     methods.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_r2d_mfpre(mat,fac,indx,rptr,rpt2,nqty,nrb,
     $                              ntotb,precon)
      USE time
      USE factor_type_mod

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(matrix_factor_type), INTENT(IN) :: fac
      TYPE(vector_type), DIMENSION(:), POINTER :: rptr,rpt2 ! work space
      INTEGER(i4), INTENT(IN) :: indx,nqty,nrb,ntotb
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: ibl,iv,ix,iy
c-----------------------------------------------------------------------
c     interface block for 2D preconditioner driver.
c-----------------------------------------------------------------------
      INCLUDE "iter_precon_intf.inc"
c-----------------------------------------------------------------------
c     if no preconditioning is used, just copy the current residual.
c-----------------------------------------------------------------------
      IF (precon=='no prec') THEN
        DO ibl=1,ntotb
          pinv(indx)%base(ibl)=ep_d(indx)%base(ibl)
        ENDDO
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     if interior coefficients of the preconditioning matrix are
c     eliminated prior to factorization, find the rhs vector for the
c     modified system.  first multiply block-border coefficients
c     by their weights to mimic finite_elements operations.  also
c     save original values along block borders.
c-----------------------------------------------------------------------
      IF (mat%eliminated.AND.poly_deg>1) THEN
        DO ibl=1,ntotb
          rpt2(ibl)=ep_d(indx)%base(ibl)
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            rpt2(ibl)%arr(:,ix,iy)=rpt2(ibl)%arr(:,ix,iy)*
     $        seam(ibl)%vertex(iv)%ave_factor
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            IF (seam(ibl)%segment(iv)%h_side) THEN
              rpt2(ibl)%arrh(:,:,ix,iy)=rpt2(ibl)%arrh(:,:,ix,iy)*
     $          seam(ibl)%segment(iv)%ave_factor
            ELSE
              rpt2(ibl)%arrv(:,:,ix,iy)=rpt2(ibl)%arrv(:,:,ix,iy)*
     $          seam(ibl)%segment(iv)%ave_factor
            ENDIF
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       eliminate interior coefficients.
c-----------------------------------------------------------------------
        CALL timer(timestart)
        CALL matelim_presolve(mat,rpt2,rptr,nqty)
        CALL timer(timeend)
        time_stcon=time_stcon+timeend-timestart
c-----------------------------------------------------------------------
c       combine contributions across block borders.
c-----------------------------------------------------------------------
        DO ibl=1,ntotb
          CALL edge_load_arr(rptr(ibl),nqty,poly_deg-1_i4,seam(ibl))
        ENDDO
        CALL edge_network(nqty,0_i4,poly_deg-1_i4,.false.)
        DO ibl=1,ntotb
          CALL edge_unload_arr(rptr(ibl),nqty,poly_deg-1_i4,seam(ibl))
        ENDDO
c-----------------------------------------------------------------------
c       call the preconditioning routine with the reduced vector.
c-----------------------------------------------------------------------
        CALL iter_pre_real(mat,fac,rptr,pinv(indx)%base,
     $                     rpt2,nqty,nrb,ntotb,precon)
      ELSE
c-----------------------------------------------------------------------
c       call the preconditioning routine with the full vector.
c-----------------------------------------------------------------------
        CALL iter_pre_real(mat,fac,ep_d(indx)%base,pinv(indx)%base,
     $                     rpt2,nqty,nrb,ntotb,precon)
      ENDIF
c-----------------------------------------------------------------------
c     reconstruct the interior coefficients.
c-----------------------------------------------------------------------
      IF (mat%eliminated.AND.poly_deg>1) THEN
        CALL timer(timestart)
        CALL matelim_postsolve(mat,pinv(indx)%base,rptr,nqty)
        CALL timer(timeend)
        time_stcon=time_stcon+timeend-timestart

        DO ibl=1,ntotb
          pinv(indx)%base(ibl)%arri(:,:,:,:)=rptr(ibl)%arri(:,:,:,:)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gmr_r2d_mfpre

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_gmres_r2d
