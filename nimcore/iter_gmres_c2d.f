c-----------------------------------------------------------------------
c     file iter_gmres_c2d.f (formerly iter_dir_nonsym.f)
c
c     This module contains a GMRES solver for non-Hermitian, complex
c     2D systems.  Preconditioning steps uses the methods that are
c     available in the iter_precon_comp module.
c 
c     The GMRES algorithm is from "Iterative Methods for Sparse Linear
c     Systems," 2nd ed., Y. Saad, SIAM (2003).
c
c     Note that the storage and factorization are identical to what
c     is used for iter_cg_comp with a direct solve as the
c     preconditioner, so iter_fac_alloc_comp and iter_factor_comp
c     should be used for matrices that are solved with iter_dir_solve.
c
c     The iter_dir_init routine now performs a seaming operation
c     to ensure that the provided guess has absolutely no discrepancies
c     among different images along block borders.  This prevents an
c     unusual error that otherwise appears in some computations.
c     11/26/07, C. Sovinec
c
c     Two changes have been made to improve accuracy with ill-
c     conditioned systems.  First, the set of basis vectors with
c     the preconditioner applied are saved separately.  Second, a
c     reorthogonalization step may be used in iter_dir_ortho.  More
c     details are provided in the change log for iter_3d_cg.f.
c     1/4/08, C. Sovinec
c
c     Development for the vblock_from344 branch made a couple of updates
c     that have been ported back to the nimuw trunk.  This includes
c     passing the rhs and sln data structures into dummy structures that
c     are not pointers.
c     12/29/18, C. Sovinec
c
c     As part of the iterative-solver clean-up, this file and module
c     are being renamed from iter_dir_nonsym (it has been a misnomer
c     for a long time) to iter_gmres_c2d.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_gmres_c2d.
c     1. iter_gmr_c2d_solve.
c     2. iter_gmr_c2d_init.
c     3. iter_gmr_c2d_dealloc.
c     4. iter_c2d_ortho.
c     5. iter_c2d_grot.
c     6. iter_c2d_min.
c-----------------------------------------------------------------------
c     module containing the conjugate gradient precon.
c-----------------------------------------------------------------------
      MODULE iter_gmres_c2d
      USE local
      USE edge_type_mod
      USE matrix_type_mod
      USE vector_type_mod
      USE matrix_mod
      USE edge
      USE seam_storage_mod
      USE iter_utils
      USE iter_comp_direct
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     module-specific definitions.
c-----------------------------------------------------------------------
      INTEGER(i4), PRIVATE :: nrb,ntotb,poly_deg,poly_d2
      COMPLEX(r8), PRIVATE, PARAMETER :: c_one=1._r8
      LOGICAL, PRIVATE :: use_int
      LOGICAL, PARAMETER, PRIVATE :: writecheck=.false.

      TYPE :: gm_basis_type
        TYPE(cvector_2D_type), POINTER, DIMENSION(:) :: base
      END TYPE gm_basis_type

      INTEGER(i4) :: n_gm
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: hess,givn
      COMPLEX(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: errvec,minvec

      TYPE(gm_basis_type), DIMENSION(:), POINTER, PRIVATE :: ep_d,pinv

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iter_gmr_c2d_solve.
c     manages the preconditioned GMRES solve for a single Fourier
c     component.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_c2d_solve(mat_str,fac_str,sln,rhs,nqty,tol,
     $                              maxit,precon,err,its,seed)
      USE matrix_mod
      USE factor_type_mod
      USE mpi_nim
      USE pardata
      USE time

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      TYPE(complex_factor_type), INTENT(IN) :: fac_str
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: rhs
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(INOUT) :: sln
      INTEGER(i4), INTENT(IN) :: nqty,maxit
      REAL(r8), INTENT(IN) :: tol
      CHARACTER(*), INTENT(IN) :: precon
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res
      INTEGER(i4) :: ibl,it
      INTEGER(i4) :: ix,iy,iv,mv,iq
      REAL(r8) :: er0
      CHARACTER(64) :: msg

      REAL(r8) :: tmp,timestart_dir,timeend_dir
      INTEGER(i4) :: ierror,maxloop
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
     $    ('iter_gmr_c2d_solve: nonsymmetric precon choice '//
     $     'is not available.')
c-----------------------------------------------------------------------
c     allocate storage; determine best guess and initial error.
c-----------------------------------------------------------------------
      CALL timer(timestart_dir)
      IF (.NOT.ALLOCATED(hess)) THEN
        ALLOCATE(hess(maxit+1,maxit))
        ALLOCATE(errvec(maxit+1))
        ALLOCATE(minvec(maxit))
        ALLOCATE(givn(2,maxit+1))
        ALLOCATE(ep_d(maxit+1))
        ALLOCATE(pinv(SIZE(ep_d)))
      ENDIF
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      CALL iter_gmr_c2d_init(mat_str,sln,rhs,nqty,precon,er0,tol,res)
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
        CALL iter_pre_comp(mat_str,fac_str,ep_d(it)%base,pinv(it)%base,
     $                     res,nqty,nrb,ntotb,precon)
c-----------------------------------------------------------------------
c       compute the inner product with the matrix.
c-----------------------------------------------------------------------
        CALL matvec(mat_str,pinv(it)%base,ep_d(it+1)%base,nqty)
c-----------------------------------------------------------------------
c       orthogonalize and add another column to the Hessenberg matrix.
c-----------------------------------------------------------------------
        CALL iter_c2d_ortho(it)
c-----------------------------------------------------------------------
c       find the Givens rotation for the last Hessenberg entry
c-----------------------------------------------------------------------
        CALL iter_c2d_grot(it)
c-----------------------------------------------------------------------
c       check tolerance.
c-----------------------------------------------------------------------
        n_gm=it
        IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,i4,a,es15.8)')
     $     '   GMR at step',it,' is ',ABS(errvec(it+1))
        IF (ABS(errvec(it+1))/er0<=tol) THEN
          err=errvec(it+1)
          EXIT
        ENDIF
      ENDDO itloop
c-----------------------------------------------------------------------
c     determine the minimizing vector of coefficients, and 
c     construct the solution with the saved (P^-1 . basis) vectors.
c-----------------------------------------------------------------------
      IF (n_gm>=1) THEN
        CALL iter_c2d_min
        DO it=1,n_gm
          DO ibl=1,ntotb
            CALL vector_add(sln(ibl),pinv(it)%base(ibl),
     $                      c_one,minvec(it))
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

      CALL iter_gmr_c2d_dealloc(res)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(comm_layer,ierror)
      CALL timer(timeend_dir)
      time_iter = time_iter + timeend_dir-timestart_dir

      RETURN
      END SUBROUTINE iter_gmr_c2d_solve
c-----------------------------------------------------------------------
c     subprogram 2. iter_gmr_c2d_init.
c     allocates storage used during the iterative solve.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_c2d_init(mat,sln,rhs,nqty,precon,er0,tol,res)
      USE mpi_nim
      USE pardata

      TYPE(complex_matrix_type), INTENT(IN) :: mat
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: rhs
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(INOUT) :: sln
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res
      INTEGER(i4), INTENT(IN) :: nqty
      CHARACTER(*), INTENT(IN) :: precon
      REAL(r8), INTENT(OUT) :: er0
      REAL(r8), INTENT(IN) :: tol

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
c     data at cell interiors.
c-----------------------------------------------------------------------
      use_int=.NOT.mat%eliminated
c-----------------------------------------------------------------------
c     ensure that the guess has exactly the same value in different
c     images of the same node along block borders.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL edge_load_2D_carr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_cin(1:nqty)=
     $      seam(ibl)%vertex(iv)%seam_cin(1:nqty)*
     $      seam(ibl)%vertex(iv)%ave_factor
          IF (poly_deg>1) THEN
            seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1))=
     $        seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1))*
     $        seam(ibl)%segment(iv)%ave_factor
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(nqty,1_i4,poly_deg-1_i4,.false.)
      DO ibl=1,ntotb
        CALL edge_unload_2D_carr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     allocate storage for the gmres operations.
c-----------------------------------------------------------------------
      ALLOCATE(res(ntotb))
      DO ibl=1,ntotb
        mx=SIZE(rhs(ibl)%arr,2)-1
        my=SIZE(rhs(ibl)%arr,3)-1
        IF (SIZE(rhs(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_c2d_init: rhs dimension inconsistent with nqty')
        IF (SIZE(sln(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_c2d_init: sln dimension inconsistent with nqty')
        CALL vector_type_alloc(res(ibl),poly_deg,mx,my,nqty,
     $                         alloc_int=use_int)
      ENDDO
      DO igm=1,SIZE(ep_d)
        ALLOCATE(ep_d(igm)%base(ntotb))
        ALLOCATE(pinv(igm)%base(ntotb))
        DO ibl=1,ntotb
          mx=SIZE(rhs(ibl)%arr,2)-1
          my=SIZE(rhs(ibl)%arr,3)-1
          CALL vector_type_alloc(ep_d(igm)%base(ibl),
     $                           poly_deg,mx,my,nqty,alloc_int=use_int)
          CALL vector_type_alloc(pinv(igm)%base(ibl),
     $                           poly_deg,mx,my,nqty,alloc_int=use_int)
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
      END SUBROUTINE iter_gmr_c2d_init
c-----------------------------------------------------------------------
c     subprogram 3. iter_gmr_c2d_dealloc.
c     deallocates storage used during the solve.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_c2d_dealloc(res)

      INTEGER(i4) :: ibl,igm
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res
c-----------------------------------------------------------------------
c     deallocate storage.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL vector_type_dealloc(res(ibl))
      ENDDO
      DEALLOCATE(res)
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
      END SUBROUTINE iter_gmr_c2d_dealloc
c-----------------------------------------------------------------------
c     subprogram 4. iter_c2d_ortho.
c     orthogonalize with modified Gram-Schmidt.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c2d_ortho(jout)
      USE pardata

      INTEGER(i4), INTENT(IN) :: jout

      INTEGER(i4) :: ibl,ii,ierror
      COMPLEX(r8) :: cdot,csc
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
        cdot=0._r8
        DO ibl=1,ntotb
          CALL iter_dot(cdot,ep_d(ii)%base(ibl),ep_d(jout+1)%base(ibl),
     $                  seam(ibl),ibl,nrb,poly_deg,use_int)
        ENDDO
        IF (nprocs_layer>1) THEN
          CALL mpi_allreduce(cdot,hess(ii,jout),1,mpi_nim_comp,
     $         mpi_sum,comm_layer,ierror)
        ELSE
          hess(ii,jout)=cdot
        ENDIF
        DO ibl=1,ntotb
          CALL vector_add(ep_d(jout+1)%base(ibl),ep_d(ii)%base(ibl),
     $                    c_one,-hess(ii,jout))
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
          cdot=0._r8
          DO ibl=1,ntotb
            CALL iter_dot(cdot,ep_d(ii)%base(ibl),
     $                    ep_d(jout+1)%base(ibl),
     $                    seam(ibl),ibl,nrb,poly_deg,use_int)
          ENDDO
          IF (nprocs_layer>1) THEN
            CALL mpi_allreduce(cdot,csc,1,mpi_nim_comp,
     $           mpi_sum,comm_layer,ierror)
            cdot=csc
          ENDIF
          hess(ii,jout)=hess(ii,jout)+cdot
          DO ibl=1,ntotb
            CALL vector_add(ep_d(jout+1)%base(ibl),ep_d(ii)%base(ibl),
     $                      c_one,-cdot)
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
      END SUBROUTINE iter_c2d_ortho
c-----------------------------------------------------------------------
c     subprogram 5. iter_c2d_grot.
c     find and apply the Givens rotation to remove below-diagonal
c     entries from the Hessenberg matrix and update the error vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c2d_grot(jout)

      INTEGER(i4), INTENT(IN) :: jout

      INTEGER(i4) :: ii,jj
      REAL(r8) :: gden
      COMPLEX(r8) :: gs,gc,ctmp
c-----------------------------------------------------------------------
c     Apply existing rotations to the new column.
c-----------------------------------------------------------------------
      DO ii=1,jout-1
        ctmp=CONJG(givn(1,ii))*hess(ii,jout)+
     $       CONJG(givn(2,ii))*hess(ii+1,jout)
        hess(ii+1,jout)=givn(1,ii)*hess(ii+1,jout)-
     $                  givn(2,ii)*hess(ii,jout)
        hess(ii,jout)=ctmp
      ENDDO
c-----------------------------------------------------------------------
c     Compute the new Givens sine and cosine coefficients.
c-----------------------------------------------------------------------
      gden=hess(jout,jout)*CONJG(hess(jout,jout))+
     $     hess(jout+1,jout)*hess(jout+1,jout)
      gden=SQRT(gden)
      givn(1,jout)=hess(jout,jout)/gden
      givn(2,jout)=hess(jout+1,jout)/gden
c-----------------------------------------------------------------------
c     Apply the new rotation to the new Hessenberg column and the
c     error vector.
c-----------------------------------------------------------------------
      gc=CONJG(givn(1,jout))
      gs=CONJG(givn(2,jout))
      ctmp=gc*hess(jout,jout)+gs*hess(jout+1,jout)
      hess(jout+1,jout)=givn(1,jout)*hess(jout+1,jout)-
     $                  givn(2,jout)*hess(jout,jout)
      hess(jout,jout)=ctmp
      errvec(jout+1)=-givn(2,jout)*errvec(jout)
      errvec(jout)=gc*errvec(jout)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c2d_grot
c-----------------------------------------------------------------------
c     subprogram 6. iter_c2d_min.
c     find the minimizing vector of coefficients for the orthogonal
c     basis vectors.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c2d_min

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
      END SUBROUTINE iter_c2d_min

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_gmres_c2d
