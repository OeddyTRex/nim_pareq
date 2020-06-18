c-----------------------------------------------------------------------
c     file iter_krylov_c3d.f (formerly iter_3d_cg.f)
c     module that contains routines for a matrix-free 3D Krylov-space
c     solvers.  originally, this only had CG but it now allows GMRES
c     for non-Hermitian systems.  preconditioning is performed with
c     the 2D routines available in iter_precon_comp.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_ky_c3d_mod
c     1.  iter_ky_c3d_solve.
c     2.  iter_cg_c3d_solve.
c     3.  iter_gmr_c3d_solve.
c     4.  iter_c3d_init.
c     5.  iter_c3d_dealloc.
c     6.  iter_c3d_dot.
c     7.  iter_c3d_err.
c     8.  iter_c3d_ortho.
c     9.  iter_c3d_grot.
c     10. iter_c3d_min.
c     11. iter_c3d_hcheck.
c-----------------------------------------------------------------------
c     module containing the 3D cg routines.
c-----------------------------------------------------------------------
      MODULE iter_ky_c3d_mod
      USE vector_type_mod
      USE local
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      LOGICAL, PRIVATE :: symm,novec
      LOGICAL, PRIVATE, PARAMETER :: usefom=.false.,fixm=.false.,
     $         hcheck=.false.,leftprec=.false.,writecheck=.false.
      INTEGER(i4), PRIVATE, PARAMETER :: npoln=0_i4

      TYPE :: gm3_basis_type
        TYPE(cvector_type), POINTER, DIMENSION(:) :: base
      END TYPE gm3_basis_type

      INTEGER(i4), PRIVATE :: n_gm
      REAL(r8), DIMENSION(:,:), ALLOCATABLE, PRIVATE :: hess,givn
      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: errvec,minvec

      TYPE(gm3_basis_type), DIMENSION(:), POINTER, PRIVATE :: ep_d,pinv
      INTEGER(i4), PRIVATE :: nftotal

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iter_ky_c3d_solve
c     perform conjugate gradient iterations to solve a 3D matrix
c     equation using a matrix-free computation of the matrix-vector
c     product.  preconditioning is accomplished according to a user-
c     specificied subroutine.
c
c     the following are notes regarding parameters passed into the
c     iter_ky_c3d_solve routine:
 
c     1. dot_routine - the name of an external subroutine that finds
c                      the product of the matrix and the current
c                      direction vector.  its parameters are specified
c                      in the dot_routine interface block shown
c                      in the driver subroutines.

c     2. pre_routine - the name of an external subroutine that applies
c                      a preconditioning operation on a passed data
c                      structure.  its parameters are specied in the
c                      pre_routine interaface block.

c     3. sympass (in)- logical that indicates whether the system is
c                      Hermitian.

c     4. sln (inout) - cvector_type data structure that holds a user-
c                      supplied guess on input and the solution as
c                      output.

c     5. rhs (in)    - cvector_type that is the right side of the 
c                      A.sln=rhs algebraic system.

c     6. dir (inout) - cvector_type that that is used as the current
c                      direction vector during dot_routine calls.

c     7. nqty (in)   - integer point-vector size (1, 3, etc.)

c     8. poly_deg (in) - integer degree of continuous 2D polynomials
c                      used in the element expansion.

c     9. nbdsc (in)  - integer number of discontinuous bases per
c                      element.

c     10. nqdsc (in) - integer point-vector size for discontinuous
c                      bases.

c     11. nrbpass (in) - integer number of rblocks on this process.

c     12. ntbpass (in) - integer number of blocks on this process.

c     13. nfpass (in)  - integer number of Fourier components on this
c                        process.

c     14. etol (in)    - real error-tolerance specification.

c     15. maxits (in)  - integer maximum number of iterations.

c     16. nm_total (in) - integer number of Fourier components over all
c                         process-layers.

c     17. err (out)  - real norm of final residual, relative to the
c                      norm of the rhs.

c     18. seed (out) - character string showing which vector was used
c                      to start the iteration.

c     19. old_vec (in) - optional cvector_type that can be used for
c                        the reference when checking the tolerance.
 
c     previous notes on usage that are still relevant are:
c
c     a)  the direction vector (dir) must be the storage used in the
c     finite element computation of the matrix-vector product.
c
c     b)  the optional input old_vec is used to get a norm of the rhs.
c     this is used so that the solved equation can be of the form,
c     A.(x_new-x_old)=b, where the rhs for the x_new equation is
c     b+A.x_old.  the tolerance for convergence is then based on the
c     norm of the residual relative to the norm of rhs (not b).
c-----------------------------------------------------------------------
c     change log:

c     1. A restart capability has been added that checks if the error
c     increases after every nrst iterations.  If so, it goes back to the
c     old solution, throwing out the direction vector used at that step.
c     The code also compares the recursive and actual residuals.  If the
c     discrepancy is larger than resid_tol, the direction vector is
c     thrown out restarting the computation from that iterate as a
c     guess.
c       5/8/02, C. Sovinec

c     2. Calls to regular_ave and the new regular_zero_phi have been
c     added to handle n=1 regularity conditions on vectors.
c	6/11/02, C. Sovinec

c     3. A check of the actual residual has been added after the
c     iteration loop to see if the problem has converged when its
c     hits maxits.
c       6/8/03, C. Sovinec
 
c     4. This version "_ky_" uses either CG or QGMRES, depending
c     on symmetry of preconditioning matrix.
c       10/31/05 D. C. Barnes
  
c     5. The nonsymmetric solve has been modified to use the
c     GMRES algorithm from "Iterative Methods for Sparse Linear
c     Systems," 2nd ed., Y. Saad, SIAM (2003).
c     The old restart capability never worked well and has been removed
c     for better readability.
c       6/7/07, C. Sovinec
c
c     6. The iter_3d_cg_init routine now performs a seaming operation
c     to ensure that the provided guess has absolutely no discrepancies
c     among different images along block borders.  This prevents an
c     unusual error that otherwise appears in some computations.
c     11/26/07, C. Sovinec
 
c     7. Difficulties converging on a computation with time-step much
c     larger than tau_A prompted a check of the entire solver.  The
c     root problem appears to be that when the preconditioner matrix (P)
c     itself is ill-conditioned, P^-1 dotted into a linear combination
c     of the basis vectors is not equivalent to the same linear 
c     combination of P^-1 dotted into each basis vector due to finite-
c     precision math.  This is remedied by saving P^-1 dotted into
c     each basis in the new pinv data structure, as motivated by the
c     Barnes code that did not use the Hessenberg matrix but did not
c     exhibit the same problem.
 
c     In the process of sorting this out, other numerical options were
c     implemented as part of the testing.  The logical parameter usefom
c     may be coded to use Arnoldi's full orthogonalization method (FOM)
c     instead of the GMRES minimization when choosing the linear
c     combination of vectors for the solution.  Another, fixm, may be
c     set to skip the tolerance checks and run a full maxit steps
c     unless the orthogonalization runs out of new directions.  This
c     option also skips the Givens rotations and solves either the
c     orthogonalization equation (FOM) or the minimization problem
c     (GMRES) after all vectors are determined.  Finally, one may
c     choose left preconditioning by setting leftprec, but this
c     does not seem to offer the effectiveness of right preconditioning
c     when P is ill conditioned.
c     1/3/08, C. Sovinec

c     8. A polynomial approximation as been added to the preconditioning
c     routine.  As implemented, it amounts to Jacobi-like iteration
c     with the preconditioning matrix as the generalized diagonal block.
c     To apply the preconditioner more than once per GMRES iteration,
c     increase the npoln parameter.
c     1/17/08, C. Sovinec

c     9. Optional matrix structures can be supplied for off diagonal
c     in Fourier component contributions that are used in SOR and SSOR
c     preconditioning steps.
c     6/15/08, C. Sovinec

c     10. The GMRES algorithm in iter_3d_cg now includes coefficients
c     from auxiliary discontinuous fields.  The corresponding rhs
c     coefficients are passed in the arrtmp part of the vector structure
c     for the rhs.  There are also two new integers, nqdsc and nbdsc,
c     that are passed into iter_3d_ky_solve to specify, respectively,
c     the size of the quantity dimension and the number of bases for the
c     discontinuous fields.
c     8/19/13, C. Sovinec

c     11. Changes for discontinuous fields are being merged into the
c     main trunk.  At this time, the unused mat2d option is being
c     removed.
c     9/4/13, C. Sovinec

c     12. A higher-level layer of dot-product management routines is
c     now used for the get_rhs, regularity, and boundary-condition
c     calls.  This has been ported from the vblock branch, where it
c     has been in place since 7/2/14.  Also at this stage,
c     the iter_3d_ky solves no longer use input parameters through
c     the finite_element module.  With this and the external direction-
c     vector computation, the iter_3d_cg module does not depend on
c     the f_e module at all.
c     8/7/15, C. Sovinec

c     13. The polynomial approximation and SOR steps have been removed
c     from the preconditioning to facilitate 1D solves over the
c     Fourier components.

c     14. The GMRES and CG algorithms have been re-split into separate
c     subroutines.  The iter_ky_c3d_solve passes data into the correct
c     routine, according to the symmetry flag.  There are also changes,
c     throughout, for the new naming convention.
c     1/26/19

c     15. Preconditioning is now performed by an external subroutine.
c     The name of the subroutine is passed, and its parameters must
c     follow the "pre_routine" interface block that is shown below.
c     Other parameter passed into iter_ky_c3d_solve have also been
c     changed, because the matrix and factor structures are no longer
c     needed.
c-----------------------------------------------------------------------
      SUBROUTINE iter_ky_c3d_solve(dot_routine,pre_routine,sympass,
     $                             sln,rhs,dir,nqty,poly_deg,nbdsc,
     $                             nqdsc,nrbpass,ntbpass,nfpass,etol,
     $                             maxits,nm_total,err,its,seed,old_vec)
      USE matrix_type_mod
      USE factor_type_mod

      LOGICAL, INTENT(IN) :: sympass
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: sln,dir
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: rhs
      TYPE(cvector_type), DIMENSION(:), INTENT(IN), OPTIONAL :: old_vec
      INTEGER(i4), INTENT(IN) :: nqty,poly_deg,maxits,nm_total
      INTEGER(i4), INTENT(IN) :: nbdsc,nqdsc  ! nqty and # of bases for
      INTEGER(i4), INTENT(IN) :: nrbpass,ntbpass,nfpass
      REAL(r8), INTENT(IN) :: etol
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dot_routine(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE dot_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE pre_routine(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE pre_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     choose the appropriate Krylov solver.
c-----------------------------------------------------------------------
      IF (sympass) THEN
        CALL iter_cg_c3d_solve(dot_routine,pre_routine,
     $                         sln,rhs,dir,nqty,poly_deg,nbdsc,
     $                         nqdsc,nrbpass,ntbpass,nfpass,etol,
     $                         maxits,nm_total,err,its,seed,old_vec)
      ELSE
        CALL iter_gmr_c3d_solve(dot_routine,pre_routine,
     $                          sln,rhs,dir,nqty,poly_deg,nbdsc,
     $                          nqdsc,nrbpass,ntbpass,nfpass,etol,
     $                          maxits,nm_total,err,its,seed,old_vec)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_ky_c3d_solve
c-----------------------------------------------------------------------
c     subprogram 2. iter_cg_c3d_solve
c     This is the conjugate gradient solver for Hermitian complex
c     3D systems, solved with matrix-free dot products.
c-----------------------------------------------------------------------
      SUBROUTINE iter_cg_c3d_solve(dot_routine,pre_routine,
     $                             sln,rhs,dir,nqty,poly_deg,nbdsc,
     $                             nqdsc,nrbpass,ntbpass,nfpass,etol,
     $                             maxits,nm_total,err,its,seed,old_vec)
      USE matrix_mod
      USE factor_type_mod
      USE seam_storage_mod
      USE edge
      USE iter_cg
      USE time

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: sln,dir
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: rhs
      TYPE(cvector_type), DIMENSION(:), INTENT(IN), OPTIONAL :: old_vec
      INTEGER(i4), INTENT(IN) :: nqty,poly_deg,maxits,nm_total
      INTEGER(i4), INTENT(IN) :: nbdsc,nqdsc  ! nqty and # of bases for
                                              ! discontinuous fields
      INTEGER(i4), INTENT(IN) :: nrbpass,ntbpass,nfpass
      REAL(r8), INTENT(IN) :: etol
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(cvector_type), DIMENSION(:), POINTER :: ctmp,res

      INTEGER(i4) :: ibl,it3d,ierror,nfour,ntotb,nrb

      REAL(r8), PARAMETER :: resid_tol=1._r8,resid_margin=0.5_r8
      REAL(r8) :: ztr=0.,ztro,alp,er0,erg,rhsnorm
      REAL(r8) :: timestart_cg,timeend_cg

      LOGICAL, PARAMETER :: residcheck=.true.
      LOGICAL, PARAMETER :: everycheck=.false.
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dot_routine(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE dot_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE pre_routine(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE pre_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     start the timer, and set the symmetry flag and problem dimensions.
c-----------------------------------------------------------------------
      CALL timer(timestart_cg)
      symm=.true.
      nrb=nrbpass
      ntotb=ntbpass
      nfour=nfpass
c-----------------------------------------------------------------------
c     call the initialization routine.
c-----------------------------------------------------------------------
      nftotal=nm_total
      CALL iter_c3d_init(ctmp,nfour,nrb,ntotb,poly_deg,
     $                   nqty,nbdsc,nqdsc,sln,res)
      DO ibl=1,ntotb
        res(ibl)=rhs(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     find the norm for rhs.
c-----------------------------------------------------------------------
      IF (PRESENT(old_vec)) THEN
        DO ibl=1,ntotb
          dir(ibl)=old_vec(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.false.)
        DO ibl=1,ntotb
          CALL vector_add(ctmp(ibl),res(ibl))
        ENDDO
        CALL iter_c3d_err(rhsnorm,nrb,ntotb,poly_deg,ctmp)
      ELSE
        CALL iter_c3d_err(rhsnorm,nrb,ntotb,poly_deg,res)
      ENDIF
      rhsnorm=MAX(rhsnorm,SQRT(TINY(rhsnorm)))
      CALL iter_c3d_err(er0,nrb,ntotb,poly_deg,res)
c-----------------------------------------------------------------------
c     check if the guess provides a better starting point, i.e. minimize
c     error in the direction of the guess, a la DCB
c-----------------------------------------------------------------------
      CALL iter_c3d_dot(erg,sln,sln,seam,ntotb)
      IF (erg>0._r8) THEN
        DO ibl=1,ntotb
          dir(ibl)=sln(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.true.)

        CALL iter_c3d_dot(alp,ctmp,res,seam,ntotb)
        CALL iter_c3d_dot(erg,ctmp,ctmp,seam,ntotb)
        alp=alp/erg
        DO ibl=1,ntotb
          CALL vector_mult(sln(ibl),alp)
          CALL vector_add(ctmp(ibl),res(ibl),v1fac=-alp)
        ENDDO
      ELSE
        DO ibl=1,ntotb
          ctmp(ibl)=res(ibl)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     decide whether to use this improved guess.
c-----------------------------------------------------------------------
      CALL iter_c3d_err(erg,nrb,ntotb,poly_deg,ctmp)
      DO ibl=1,ntotb
        IF (erg<er0) THEN
          res(ibl)=ctmp(ibl)
c         sln(ibl)=dir(ibl)
          seed="guess"
          err=erg
        ELSE
          sln(ibl)=0._r8
          seed="0 vec"
          err=er0
        ENDIF
      ENDDO
      er0=rhsnorm
c-----------------------------------------------------------------------
c     catch trivial cases.
c-----------------------------------------------------------------------
      IF (err==0._r8) THEN
        CALL iter_c3d_dealloc(ctmp,res)
        its=0
        CALL timer(timeend_cg)
        time_iter = time_iter + timeend_cg-timestart_cg
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     write norms for diagnostics.
c-----------------------------------------------------------------------
      IF (residcheck.AND.writecheck.AND.node==0)
     $  WRITE(nim_wr,'(2(a,es15.8))')
     $  '   Initial resid=',err,' rhsnorm=',rhsnorm
c-----------------------------------------------------------------------
c     start iterations.
c-----------------------------------------------------------------------
      threed_loop: DO it3d=1,maxits
c-----------------------------------------------------------------------
c       precondition the 3D equation with the "user-supplied" routine.
c-----------------------------------------------------------------------
        CALL pre_routine(res,ctmp,it3d,.false.)
c-----------------------------------------------------------------------
c       CG-specific computations  find z^t.r
c-----------------------------------------------------------------------
        CALL iter_c3d_dot(ztr,ctmp,res,seam,ntotb)
c-----------------------------------------------------------------------
c       find the new direction vector.
c-----------------------------------------------------------------------
        IF (it3d==1) THEN
          DO ibl=1,ntotb
            dir(ibl)=ctmp(ibl)
          ENDDO
        ELSE
          DO ibl=1,ntotb
            CALL vector_add(dir(ibl),ctmp(ibl),ztr/ztro,1._r8)
          ENDDO
        ENDIF
        ztro=ztr
c-----------------------------------------------------------------------
c       find the dot product of the 3D matrix and the new direction
c       vector.
c-----------------------------------------------------------------------
        CALL dot_routine(dir,ctmp,.false.)
c-----------------------------------------------------------------------
c       CG-specific: find the new solution iterate and residual.
c-----------------------------------------------------------------------
        CALL iter_c3d_dot(alp,dir,ctmp,seam,ntotb)
        alp=ztr/alp
        DO ibl=1,ntotb
          CALL vector_add(sln(ibl),dir(ibl),1._r8,alp)
          CALL vector_add(res(ibl),ctmp(ibl),1._r8,-alp)
        ENDDO
c-----------------------------------------------------------------------
c       check error level.
c-----------------------------------------------------------------------
        CALL iter_c3d_err(err,nrb,ntotb,poly_deg,res)
        n_gm=it3d
        IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,i4,a,es15.8)')
     $     '   CG at step',it3d,' is ',err
        IF (err/er0<=resid_margin*etol) EXIT threed_loop
      ENDDO threed_loop 
      its=n_gm
c-----------------------------------------------------------------------
c     residual check.
c-----------------------------------------------------------------------
      IF (residcheck) THEN
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,es15.8)') '   Final iterate  3d resid=',err
        DO ibl=1,ntotb
          dir(ibl)=sln(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.false.)
        DO ibl=1,ntotb
          CALL vector_add(ctmp(ibl),rhs(ibl),v1fac=-1._r8)
        ENDDO
        CALL iter_c3d_err(err,nrb,ntotb,poly_deg,ctmp)
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,es15.8)') '   Final separate 3d resid=',err
      ENDIF
      err=err/rhsnorm
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      CALL iter_c3d_dealloc(ctmp,res)
      CALL timer(timeend_cg)
      time_iter = time_iter + timeend_cg-timestart_cg
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_cg_c3d_solve
c-----------------------------------------------------------------------
c     subprogram 3. iter_gmr_c3d_solve
c     The is the generalized minimum residual solver for non-Hermitian
c     complex 3D systems, solved with matrix-free dot products.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gmr_c3d_solve(dot_routine,pre_routine,sln,rhs,
     $                              dir,nqty,poly_deg,nbdsc,nqdsc,
     $                              nrbpass,ntbpass,nfpass,etol,maxits,
     $                              nm_total,err,its,seed,old_vec)
      USE matrix_mod
      USE factor_type_mod
      USE seam_storage_mod
      USE edge
      USE iter_cg
      USE time

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: sln,dir
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: rhs
      TYPE(cvector_type), DIMENSION(:), INTENT(IN), OPTIONAL :: old_vec
      INTEGER(i4), INTENT(IN) :: nqty,poly_deg,maxits,nm_total
      INTEGER(i4), INTENT(IN) :: nbdsc,nqdsc  ! nqty and # of bases for
                                              ! discontinuous fields
      INTEGER(i4), INTENT(IN) :: nrbpass,ntbpass,nfpass
      REAL(r8), INTENT(IN) :: etol
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(cvector_type), DIMENSION(:), POINTER :: ctmp,res

      INTEGER(i4) :: ibl,it3d,ierror,nfour,ntotb,nrb

      REAL(r8), PARAMETER :: resid_tol=1._r8,resid_margin=0.5_r8
      REAL(r8) :: alp,er0,erg,rhsnorm
      REAL(r8) :: timestart_cg,timeend_cg
      LOGICAL, SAVE :: first_call=.true.

      LOGICAL, PARAMETER :: residcheck=.true.
      LOGICAL, PARAMETER :: everycheck=.false.
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dot_routine(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE dot_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE pre_routine(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE pre_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     start the timer, set the symmetry flag, set problem dimensions,
c     and allocate space for GMRES algebra.
c-----------------------------------------------------------------------
      CALL timer(timestart_cg)
      symm=.false.
      nrb=nrbpass
      ntotb=ntbpass
      nfour=nfpass
      IF (first_call) THEN
        ALLOCATE(hess(maxits+1,maxits))
        ALLOCATE(errvec(maxits+1))
        ALLOCATE(minvec(maxits))
        ALLOCATE(givn(2,maxits+1))
        ALLOCATE(ep_d(maxits+1))
        IF (.NOT.leftprec) ALLOCATE(pinv(SIZE(ep_d)))
        first_call=.false.
      ENDIF
c-----------------------------------------------------------------------
c     call the initialization routine.
c-----------------------------------------------------------------------
      nftotal=nm_total
      CALL iter_c3d_init(ctmp,nfour,nrb,ntotb,poly_deg,
     $                   nqty,nbdsc,nqdsc,sln,res)
      DO ibl=1,ntotb
        res(ibl)=rhs(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     find the norm for rhs.
c-----------------------------------------------------------------------
      IF (PRESENT(old_vec)) THEN
        DO ibl=1,ntotb
          dir(ibl)=old_vec(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.false.)
        DO ibl=1,ntotb
          CALL vector_add(ctmp(ibl),res(ibl))
        ENDDO
        CALL iter_c3d_err(rhsnorm,nrb,ntotb,poly_deg,ctmp)
      ELSE
        CALL iter_c3d_err(rhsnorm,nrb,ntotb,poly_deg,res)
      ENDIF
      rhsnorm=MAX(rhsnorm,SQRT(TINY(rhsnorm)))
      CALL iter_c3d_err(er0,nrb,ntotb,poly_deg,res)
c-----------------------------------------------------------------------
c     check if the guess provides a better starting point, i.e. minimize
c     error in the direction of the guess, a la DCB
c-----------------------------------------------------------------------
      CALL iter_c3d_dot(erg,sln,sln,seam,ntotb)
      IF (erg>0._r8) THEN
        DO ibl=1,ntotb
          dir(ibl)=sln(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.true.)

        CALL iter_c3d_dot(alp,ctmp,res,seam,ntotb)
        CALL iter_c3d_dot(erg,ctmp,ctmp,seam,ntotb)
        alp=alp/erg
        DO ibl=1,ntotb
          CALL vector_mult(sln(ibl),alp)
          CALL vector_add(ctmp(ibl),res(ibl),v1fac=-alp)
        ENDDO
      ELSE
        DO ibl=1,ntotb
          ctmp(ibl)=res(ibl)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     decide whether to use this improved guess.
c-----------------------------------------------------------------------
      CALL iter_c3d_err(erg,nrb,ntotb,poly_deg,ctmp)
      DO ibl=1,ntotb
        IF (erg<er0) THEN
          res(ibl)=ctmp(ibl)
c         sln(ibl)=dir(ibl)
          seed="guess"
          err=erg
        ELSE
          sln(ibl)=0._r8
          seed="0 vec"
          err=er0
        ENDIF
      ENDDO
      er0=rhsnorm
c-----------------------------------------------------------------------
c     catch trivial cases.
c-----------------------------------------------------------------------
      IF (err==0._r8) THEN
        CALL iter_c3d_dealloc(ctmp,res)
        its=0
        CALL timer(timeend_cg)
        time_iter = time_iter + timeend_cg-timestart_cg
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     GMRES-specific initialization.  left preconditioning is applied
c     here before saving the initial residual and error.
c-----------------------------------------------------------------------
      IF (residcheck.AND.writecheck.AND.node==0)
     $  WRITE(nim_wr,'(2(a,es15.8))')
     $  '   Initial resid=',err,' rhsnorm=',rhsnorm
      IF (leftprec) THEN
        CALL pre_routine(res,ctmp,1_i4,.false.)
        DO ibl=1,ntotb
          res(ibl)=ctmp(ibl)
        ENDDO
        CALL iter_c3d_err(err,nrb,ntotb,poly_deg,res)
        err=MAX(err,SQRT(TINY(err)))
        IF (residcheck.AND.writecheck.AND.node==0)
     $    WRITE(nim_wr,'(2(a,es15.8))') '   Left resid=',err
        er0=err
      ENDIF

      errvec(1)=err
      DO ibl=1,ntotb
        CALL vector_mult(res(ibl),1._r8/err)
        ep_d(1)%base(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     start iterations.
c-----------------------------------------------------------------------
      threed_loop: DO it3d=1,maxits
c-----------------------------------------------------------------------
c       precondition the 3D equation with the "user-supplied" routine.
c       for GMRES, this is a right-preconditioning step.  flexible-
c       GMRES is possible with the pinv storage.
c-----------------------------------------------------------------------
        IF (.NOT.leftprec) THEN
          CALL pre_routine(res,ctmp,it3d,.true.)
          DO ibl=1,ntotb
            pinv(it3d)%base(ibl)=ctmp(ibl)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       find the new direction vector.
c-----------------------------------------------------------------------
        IF (leftprec) THEN
          DO ibl=1,ntotb
            dir(ibl)=res(ibl)
          ENDDO
        ELSE
          DO ibl=1,ntotb
            dir(ibl)=ctmp(ibl)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       find the dot product of the 3D matrix and the new direction
c       vector.
c-----------------------------------------------------------------------
        CALL dot_routine(dir,ctmp,.false.)
c-----------------------------------------------------------------------
c       GMRES: orthogonalize and add another column to the
c       Hessenberg matrix.  then, find the Givens rotation for the last
c       Hessenberg entry and determine the minimizing vector of
c       coefficients.
c
c       the left preconditioning step is applied here, after the matrix-
c       vector product.  here, ep_d is only used for temporary storage.
c-----------------------------------------------------------------------
        IF (leftprec) THEN
          DO ibl=1,ntotb
            ep_d(it3d+1)%base(ibl)=ctmp(ibl)
          ENDDO
          CALL pre_routine(ep_d(it3d+1)%base,ctmp,it3d,.false.)
        ENDIF
        CALL iter_c3d_ortho(it3d,ntotb,ctmp)
        IF (.NOT.fixm) CALL iter_c3d_grot(it3d)
c-----------------------------------------------------------------------
c       check error level.
c-----------------------------------------------------------------------
        err=ABS(errvec(n_gm+1))
        IF (everycheck.AND.node==0) WRITE(nim_wr,'(a,i4,a,es15.8)')
     $      '   GMR at step',it3d,' is ',err
        IF ((.NOT.fixm.OR.novec).AND.
     $      err/er0<=resid_margin*etol) THEN
          EXIT threed_loop
        ELSE
          DO ibl=1,ntotb
            res(ibl)=ctmp(ibl)
          ENDDO
        ENDIF
      ENDDO threed_loop 
c-----------------------------------------------------------------------
c     recompute the Hessenberg matrix for testing.
c-----------------------------------------------------------------------
      IF (hcheck.AND.node==0) 
     $  CALL iter_c3d_hcheck(dot_routine,pre_routine,ctmp,dir,
     $                       nqty,poly_deg,nbdsc,nqdsc,nrb,ntotb,nfour)
c-----------------------------------------------------------------------
c     for GMRES, determine the minimizing vector of coefficients, and
c     construct the solution from either the basis vectors for left
c     preconditioning or from P^-1 dotted into the basis vectors for
c     right preconditioning, where P^-1 is the preconditioning
c     operation.
c-----------------------------------------------------------------------
      IF (n_gm>=1) THEN
        CALL iter_c3d_min
        IF (leftprec) THEN
          DO it3d=1,n_gm
            DO ibl=1,ntotb
              CALL vector_add(sln(ibl),ep_d(it3d)%base(ibl),
     $                        v2fac=minvec(it3d))
            ENDDO
          ENDDO
        ELSE
          DO it3d=1,n_gm
            DO ibl=1,ntotb
              CALL vector_add(sln(ibl),pinv(it3d)%base(ibl),
     $                        v2fac=minvec(it3d))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      its=n_gm
c-----------------------------------------------------------------------
c     residual check.
c-----------------------------------------------------------------------
      IF (residcheck) THEN
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,es15.8)') '   Final iterate  3d resid=',err
        DO ibl=1,ntotb
          dir(ibl)=sln(ibl)
        ENDDO
        CALL dot_routine(dir,ctmp,.false.)
        DO ibl=1,ntotb
          CALL vector_add(ctmp(ibl),rhs(ibl),v1fac=-1._r8)
        ENDDO
        CALL iter_c3d_err(err,nrb,ntotb,poly_deg,ctmp)
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,es15.8)') '   Final separate 3d resid=',err
      ENDIF
      err=err/rhsnorm
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      CALL iter_c3d_dealloc(ctmp,res)
      CALL timer(timeend_cg)
      time_iter = time_iter + timeend_cg-timestart_cg
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gmr_c3d_solve
c-----------------------------------------------------------------------
c     subprogram 4. iter_c3d_init. 
c     allocate space for the linear algebra.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_init(ctmp,nfour,nrb,ntotb,
     $                         poly_deg,nqty,nbdsc,nqdsc,sln,res)
      USE matrix_type_mod
      USE factor_type_mod
      USE seam_storage_mod
      USE edge

      TYPE(cvector_type), DIMENSION(:), POINTER :: ctmp,res
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: sln
      INTEGER(i4), INTENT(IN) :: nfour,nrb,ntotb
      INTEGER(i4), INTENT(IN) :: poly_deg,nqty,nbdsc,nqdsc

      INTEGER(i4) :: ibl,mxb,myb,igm,iq,ix,iy,ibase,im,jm,iv,nbf,
     $               nsnd,isnd,ioff,nsize,ncomb
c-----------------------------------------------------------------------
c     find array dimensions.
c-----------------------------------------------------------------------
      ncomb=nqty*(poly_deg-1)**2+nqdsc*nbdsc
c-----------------------------------------------------------------------
c     ensure that the guess has exactly the same value in different
c     images of the same node along block borders.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL edge_load_carr(sln(ibl),nqty,1_i4,nfour,
     $                      poly_deg-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_cin(1:nqty*nfour)=
     $      seam(ibl)%vertex(iv)%seam_cin(1:nqty*nfour)*
     $      seam(ibl)%vertex(iv)%ave_factor
          IF (poly_deg>1) THEN
            seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1)*nfour)=
     $        seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1)*nfour)*
     $        seam(ibl)%segment(iv)%ave_factor
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(nqty,nfour,poly_deg-1_i4,.false.)
      DO ibl=1,ntotb
        CALL edge_unload_carr(sln(ibl),nqty,1_i4,nfour,
     $                        poly_deg-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     allocate vector structures.
c-----------------------------------------------------------------------
      ALLOCATE(ctmp(ntotb))
      ALLOCATE(res(ntotb))
      DO ibl=1,nrb
        mxb=SIZE(sln(ibl)%arr,2)-1
        myb=SIZE(sln(ibl)%arr,3)-1
        IF (nqdsc>0) THEN
          CALL vector_type_alloc(ctmp(ibl),poly_deg,mxb,myb,nqty,nfour,
     $                           nbdsc,nqdsc)
          CALL vector_type_alloc(res(ibl),poly_deg,mxb,myb,nqty,nfour,
     $                           nbdsc,nqdsc)
        ELSE
          CALL vector_type_alloc(ctmp(ibl),poly_deg,mxb,myb,nqty,nfour)
          CALL vector_type_alloc(res(ibl),poly_deg,mxb,myb,nqty,nfour)
        ENDIF
      ENDDO
      DO ibl=nrb+1,ntotb
        mxb=SIZE(sln(ibl)%arr,2)-1
        CALL vector_type_alloc(ctmp(ibl),1_i4,mxb,0_i4,nqty,nfour)
        CALL vector_type_alloc(res(ibl),1_i4,mxb,0_i4,nqty,nfour)
      ENDDO
c-----------------------------------------------------------------------
c     GMRES data structures.
c-----------------------------------------------------------------------
      IF (.NOT.symm) THEN
        DO igm=1,SIZE(ep_d)
          ALLOCATE(ep_d(igm)%base(ntotb))
          IF (.NOT.leftprec) ALLOCATE(pinv(igm)%base(ntotb))
          DO ibl=1,nrb
            mxb=SIZE(sln(ibl)%arr,2)-1
            myb=SIZE(sln(ibl)%arr,3)-1
            IF (nqdsc>0) THEN
              CALL vector_type_alloc(ep_d(igm)%base(ibl),
     $                               poly_deg,mxb,myb,nqty,nfour,
     $                               nbdsc,nqdsc)
              IF (.NOT.leftprec)
     $          CALL vector_type_alloc(pinv(igm)%base(ibl),
     $                                 poly_deg,mxb,myb,nqty,nfour,
     $                                 nbdsc,nqdsc)
            ELSE
              CALL vector_type_alloc(ep_d(igm)%base(ibl),
     $                               poly_deg,mxb,myb,nqty,nfour)
              IF (.NOT.leftprec)
     $          CALL vector_type_alloc(pinv(igm)%base(ibl),
     $                                 poly_deg,mxb,myb,nqty,nfour)
            ENDIF
          ENDDO
          DO ibl=nrb+1,ntotb
            mxb=SIZE(sln(ibl)%arr,2)-1
            CALL vector_type_alloc(ep_d(igm)%base(ibl),
     $                             1_i4,mxb,0_i4,nqty,nfour)
            IF (.NOT.leftprec)
     $        CALL vector_type_alloc(pinv(igm)%base(ibl),
     $                               1_i4,mxb,0_i4,nqty,nfour)
          ENDDO
        ENDDO
        hess=0._r8
        errvec=0._r8
        minvec=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     initialize data structures for 2D preconditioning operations for
c     each Fourier component.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d_init(nqty,poly_deg,nbdsc,nqdsc,nrb,ntotb,
     $                           sln)
c-----------------------------------------------------------------------
c     initialize the novec flag.
c-----------------------------------------------------------------------
      novec=.false.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_init
c-----------------------------------------------------------------------
c     subprogram 5. iter_c3d_dealloc. 
c     deallocate temporary space.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_dealloc(ctmp,res)

      TYPE(cvector_type), DIMENSION(:), POINTER :: ctmp,res

      INTEGER(i4) :: ibl,igm,im
c-----------------------------------------------------------------------
c     deallocate vector structures.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(ctmp)
        CALL vector_type_dealloc(ctmp(ibl))
        CALL vector_type_dealloc(res(ibl))
      ENDDO
      IF (.NOT.symm) THEN
        DO igm=1,SIZE(ep_d)
          DO ibl=1,SIZE(ctmp)
            CALL vector_type_dealloc(ep_d(igm)%base(ibl))
            IF (.NOT.leftprec)
     $        CALL vector_type_dealloc(pinv(igm)%base(ibl))
          ENDDO
          DEALLOCATE(ep_d(igm)%base)
          IF (.NOT.leftprec) DEALLOCATE(pinv(igm)%base)
        ENDDO
      ENDIF
      DEALLOCATE(ctmp,res)
c-----------------------------------------------------------------------
c     deallocate 2D data structures.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_dealloc
c-----------------------------------------------------------------------
c     subprogram 6. iter_c3d_dot. 
c     compute the inner product of two 3d vectors.
c
c     include the contribution of the conjugate terms.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_dot(dprod,v1,v2,seam,nbl)
      USE edge_type_mod

      REAL(r8), INTENT(OUT) :: dprod
      TYPE(cvector_type), DIMENSION(:) :: v1,v2
      TYPE(edge_type), DIMENSION(:), POINTER :: seam
      INTEGER(i4), INTENT(IN) :: nbl

      INTEGER(i4) :: ix,iy,iv,ni,ibl,ierror
      REAL(r8) :: tmp
c-----------------------------------------------------------------------
c     loop over blocks, and collect the on-processor contributions
c     from all Fourier components.
c-----------------------------------------------------------------------
      dprod=0._r8
      DO ibl=1,nbl
c-----------------------------------------------------------------------
c       find the dot product at each grid vertex.
c-----------------------------------------------------------------------
        IF (ilayer==0) THEN
          dprod=dprod+SUM(CONJG(v1(ibl)%arr(:,:,:,1))*
     $                    v2(ibl)%arr(:,:,:,1))+
     $          2._r8*SUM(CONJG(v1(ibl)%arr(:,:,:,2:))*
     $                    v2(ibl)%arr(:,:,:,2:))
        ELSE
          dprod=dprod+2._r8*SUM(CONJG(v1(ibl)%arr)*v2(ibl)%arr)
        ENDIF
c-----------------------------------------------------------------------
c       boundary points must be divided by the number of internal 
c       representations to get the correct sum over all blocks.
c       remove what has been contributed above, also.
c-----------------------------------------------------------------------
        IF (ilayer==0) THEN
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
              dprod=dprod+(seam(ibl)%vertex(iv)%ave_factor-1)*
     $               (      SUM(CONJG(v1(ibl)%arr(:,ix,iy,1))*
     $                          v2(ibl)%arr(:,ix,iy,1))+
     $                2._r8*SUM(CONJG(v1(ibl)%arr(:,ix,iy,2:))*
     $                          v2(ibl)%arr(:,ix,iy,2:)) )
          ENDDO
        ELSE
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
              dprod=dprod+2._r8*(seam(ibl)%vertex(iv)%ave_factor-1)*
     $                          SUM(CONJG(v1(ibl)%arr(:,ix,iy,:))*
     $                              v2(ibl)%arr(:,ix,iy,:))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       higher order contributions.
c-----------------------------------------------------------------------
        IF (ASSOCIATED(v1(ibl)%arrh)) THEN
          IF (ilayer==0) THEN
            dprod=dprod+SUM(CONJG(v1(ibl)%arrh(:,:,:,:,1))*
     $                            v2(ibl)%arrh(:,:,:,:,1))+
     $            2._r8*SUM(CONJG(v1(ibl)%arrh(:,:,:,:,2:))*
     $                            v2(ibl)%arrh(:,:,:,:,2:))+
     $                  SUM(CONJG(v1(ibl)%arrv(:,:,:,:,1))*
     $                            v2(ibl)%arrv(:,:,:,:,1))+
     $            2._r8*SUM(CONJG(v1(ibl)%arrv(:,:,:,:,2:))*
     $                            v2(ibl)%arrv(:,:,:,:,2:))+
     $                  SUM(CONJG(v1(ibl)%arri(:,:,:,:,1))*
     $                            v2(ibl)%arri(:,:,:,:,1))+
     $            2._r8*SUM(CONJG(v1(ibl)%arri(:,:,:,:,2:))*
     $                            v2(ibl)%arri(:,:,:,:,2:))
          ELSE
            dprod=dprod+2._r8*SUM(CONJG(v1(ibl)%arrh)*v2(ibl)%arrh)+
     $                  2._r8*SUM(CONJG(v1(ibl)%arrv)*v2(ibl)%arrv)+
     $                  2._r8*SUM(CONJG(v1(ibl)%arri)*v2(ibl)%arri)
          ENDIF
c-----------------------------------------------------------------------
c         coefficients for discontinuous bases.
c-----------------------------------------------------------------------
          IF (ASSOCIATED(v1(ibl)%arrtmp)) THEN
            IF (ilayer==0) THEN
              dprod=dprod+SUM(CONJG(v1(ibl)%arrtmp(:,:,:,:,1))*
     $                        v2(ibl)%arrtmp(:,:,:,:,1))+
     $              2._r8*SUM(CONJG(v1(ibl)%arrtmp(:,:,:,:,2:))*
     $                        v2(ibl)%arrtmp(:,:,:,:,2:))
            ELSE
              dprod=dprod+2._r8*SUM(CONJG(v1(ibl)%arrtmp)*
     $                              v2(ibl)%arrtmp)
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c         boundary segment-centered data must be divided by 2.
c-----------------------------------------------------------------------
          IF (ilayer==0) THEN
            DO iv=1,seam(ibl)%nvert
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              IF (seam(ibl)%segment(iv)%h_side) THEN
                dprod=dprod+(seam(ibl)%segment(iv)%ave_factor-1)*
     $              (      SUM(CONJG(v1(ibl)%arrh(:,:,ix,iy,1))*
     $                               v2(ibl)%arrh(:,:,ix,iy,1))+
     $               2._r8*SUM(CONJG(v1(ibl)%arrh(:,:,ix,iy,2:))*
     $                               v2(ibl)%arrh(:,:,ix,iy,2:)) )
              ELSE
                dprod=dprod+(seam(ibl)%segment(iv)%ave_factor-1)*
     $              (      SUM(CONJG(v1(ibl)%arrv(:,:,ix,iy,1))*
     $                               v2(ibl)%arrv(:,:,ix,iy,1))+
     $               2._r8*SUM(CONJG(v1(ibl)%arrv(:,:,ix,iy,2:))*
     $                               v2(ibl)%arrv(:,:,ix,iy,2:)) )
              ENDIF
            ENDDO
          ELSE
            DO iv=1,seam(ibl)%nvert
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              IF (seam(ibl)%segment(iv)%h_side) THEN
                dprod=dprod+2._r8*(seam(ibl)%segment(iv)%ave_factor-1)*
     $                SUM(CONJG(v1(ibl)%arrh(:,:,ix,iy,:))*
     $                          v2(ibl)%arrh(:,:,ix,iy,:))
              ELSE
                dprod=dprod+2._r8*(seam(ibl)%segment(iv)%ave_factor-1)*
     $                SUM(CONJG(v1(ibl)%arrv(:,:,ix,iy,:))*
     $                          v2(ibl)%arrv(:,:,ix,iy,:))
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     sum result over processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(dprod,tmp,1,mpi_nim_real,mpi_sum,
     $       mpi_comm_world,ierror)
        dprod=tmp
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_dot
c-----------------------------------------------------------------------
c     subprogram 7. iter_c3d_err. 
c     find the norm of a vector
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_err(er,nrb,ntotb,poly_deg,vec,norm_type)
      USE seam_storage_mod

      REAL(r8), INTENT(OUT) :: er
      INTEGER(i4), INTENT(IN) :: nrb,ntotb,poly_deg
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: vec
      CHARACTER(*), INTENT(IN), OPTIONAL :: norm_type

      REAL(r8) :: tmp
      INTEGER(i4) :: ibl,ierror
c-----------------------------------------------------------------------
c     the optional norm_type parameter can be used to set symm, which
c     facilitates external use of this routine.
c-----------------------------------------------------------------------
      IF (PRESENT(norm_type)) THEN
        IF (norm_type=='2-norm') THEN
          symm=.false.
        ELSE
          symm=.true.
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     determine the norm of the passed vector.
c     for consistency with older computations, CG still uses the
c     infinity norm, but GMRES uses the 2-norm.  the inf option may
c     be removed at some point.
c-----------------------------------------------------------------------
      IF (symm) THEN
        er=0
        DO ibl=1,ntotb
          er=MAX(MAXVAL(ABS(REAL(vec(ibl)%arr,r8))),
     $           MAXVAL(ABS(AIMAG(vec(ibl)%arr))),er)
          IF (poly_deg>1) THEN
            er=MAX(MAXVAL(ABS(REAL(vec(ibl)%arrh,r8))),
     $             MAXVAL(ABS(AIMAG(vec(ibl)%arrh))),er)
            er=MAX(MAXVAL(ABS(REAL(vec(ibl)%arrv,r8))),
     $             MAXVAL(ABS(AIMAG(vec(ibl)%arrv))),er)
            er=MAX(MAXVAL(ABS(REAL(vec(ibl)%arri,r8))),
     $             MAXVAL(ABS(AIMAG(vec(ibl)%arri))),er)
          ENDIF
          IF (ASSOCIATED(vec(ibl)%arrtmp)) THEN
            er=MAX(MAXVAL(ABS(REAL(vec(ibl)%arrtmp,r8))),
     $             MAXVAL(ABS(AIMAG(vec(ibl)%arrtmp))),er)
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       find the max over all processors.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          CALL mpi_allreduce(er,tmp,1,mpi_nim_real,mpi_max,
     $         mpi_comm_world,ierror)
          er=tmp
        ENDIF
c-----------------------------------------------------------------------
c     2-norm computation:
c-----------------------------------------------------------------------
      ELSE
        CALL iter_c3d_dot(tmp,vec,vec,seam,ntotb)
        er=SQRT(tmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_err
c-----------------------------------------------------------------------
c     subprogram 8. iter_c3d_ortho.
c     orthogonalize with modified Gram-Schmidt.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_ortho(jout,nbl,ep_n)
      USE seam_storage_mod

      INTEGER(i4), INTENT(IN) :: jout,nbl
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: ep_n

      INTEGER(i4) :: ibl,ii,ierror
      REAL(r8) :: rmag,rmag0,rtmp

      REAL(r8), PARAMETER :: small=1.e-14,orthtol=1.e-2
c-----------------------------------------------------------------------
c     compute the norm of A.v_jout prior to orthogonalization.
c-----------------------------------------------------------------------
      CALL iter_c3d_dot(rmag0,ep_n,ep_n,seam,nbl)
      rmag0=SQRT(rmag0)
c-----------------------------------------------------------------------
c     loop over existing basis vectors and save the inner products
c     in the hess array.
c-----------------------------------------------------------------------
      DO ii=1,jout
        CALL iter_c3d_dot(hess(ii,jout),ep_d(ii)%base,ep_n,seam,nbl)
        DO ibl=1,nbl
          CALL vector_add(ep_n(ibl),ep_d(ii)%base(ibl),
     $                    v1fac=1._r8,v2fac=-hess(ii,jout))
        ENDDO
      ENDDO
      CALL iter_c3d_dot(rmag,ep_n,ep_n,seam,nbl)
      rmag=SQRT(rmag)
c-----------------------------------------------------------------------
c     repeat orthogonalization of this vector if its magnitude is
c     reduced by more than a factor of orthtol.
c-----------------------------------------------------------------------
      IF (rmag/rmag0<orthtol) THEN
        DO ii=1,jout
          CALL iter_c3d_dot(rtmp,ep_d(ii)%base,ep_n,seam,nbl)
          hess(ii,jout)=hess(ii,jout)+rtmp
          DO ibl=1,nbl
            CALL vector_add(ep_n(ibl),ep_d(ii)%base(ibl),
     $                      v1fac=1._r8,v2fac=-rtmp)
          ENDDO
        ENDDO
        CALL iter_c3d_dot(rmag,ep_n,ep_n,seam,nbl)
        rmag=SQRT(rmag)
      ENDIF
c-----------------------------------------------------------------------
c     if there is nothing left of the new direction, the solution should
c     already be obtained.
c-----------------------------------------------------------------------
      IF (rmag<small*MAXVAL(ABS(hess(1:jout,jout)))) THEN
        n_gm=jout
        novec=.true.
        hess(jout+1,jout)=0._r8
        IF (writecheck.AND.node==0)
     $    WRITE(nim_wr,'(a,i4)') '   No new vector at step ',jout
      ELSE
        n_gm=jout
        hess(jout+1,jout)=rmag
        DO ibl=1,nbl
          CALL vector_mult(ep_n(ibl),1._r8/rmag)
          ep_d(jout+1)%base(ibl)=ep_n(ibl)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_ortho
c-----------------------------------------------------------------------
c     subprogram 9. iter_c3d_grot.
c     find and apply the Givens rotation to remove below-diagonal
c     entries from the Hessenberg matrix and update the error vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_grot(jout)

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
      hess(jout+1,jout)=gc*hess(jout+1,jout)-gs*hess(jout,jout)
      hess(jout,jout)=rtmp

      errvec(jout+1)=-gs*errvec(jout)
      errvec(jout)=gc*errvec(jout)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_grot
c-----------------------------------------------------------------------
c     subprogram 10. iter_c3d_min.
c     find the minimizing vector of coefficients for the orthogonal
c     basis vectors.
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_min
      USE math_tran

      INTEGER(i4) :: ii,jj
      LOGICAL :: singular

      REAL(r8) :: rtmp,gc,gs
c-----------------------------------------------------------------------
c     if the usefom flag is set, apply Arnoldi's full orthogonalization
c     method instead of minimizing the residual.  when the fixed
c     iteration count flag is not set, FOM is accomplished
c     by undoing the last Givens rotation and solving the truncated
c     system.
c-----------------------------------------------------------------------
      IF (usefom.AND..NOT.fixm) THEN
        gc=givn(1,n_gm)
        gs=givn(2,n_gm)
        rtmp=gc*hess(n_gm,n_gm)
        hess(n_gm+1,n_gm)=gs*hess(n_gm,n_gm)
        hess(n_gm,n_gm)=rtmp
        rtmp=gc*errvec(n_gm)-gs*errvec(n_gm+1)
        errvec(n_gm+1)=gs*errvec(n_gm)+gc*errvec(n_gm+1)
        errvec(n_gm)=rtmp
      ENDIF
c-----------------------------------------------------------------------
c     when the fixed iteration count flag is set, the Hessenberg matrix
c     is not reduced at this stage.  with FOM, we solve
c
c       H_m.minvec = errvec(1:n_gm)
c
c     and with GMRES we minimize including the last direction vector,
c
c       Hbar_m^T.Hbar_m.minvec = errvec(1:n_gm+1)
c
c     where rows 1:n_gm of Hbar_m are H_m, and the n_gm+1 row of Hbar_m
c     is zero except for the last entry.
c-----------------------------------------------------------------------
      IF (fixm) THEN
        IF (usefom) THEN
          DO ii=2,n_gm
            rtmp=hess(ii,ii-1)/hess(ii-1,ii-1)
            hess(ii,:)=hess(ii,:)-rtmp*hess(ii-1,:)
            errvec(ii)=errvec(ii)-rtmp*errvec(ii-1)
          ENDDO
        ELSE
          errvec(1:n_gm)=
     $      MATMUL(TRANSPOSE(hess(1:n_gm+1,1:n_gm)),errvec(1:n_gm+1))
          hess(1:n_gm,1:n_gm)=
     $      MATMUL(TRANSPOSE(hess(1:n_gm+1,1:n_gm)),
     $                       hess(1:n_gm+1,1:n_gm))
          CALL math_solve_q1_sym(n_gm,1_i4,hess(1:n_gm,1:n_gm),
     $                           minvec(1:n_gm),errvec(1:n_gm),'both',
     $                           singular)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     finish solving the system including the standard case where
c     usefom and fixm are both false.
c-----------------------------------------------------------------------
      IF (usefom.OR..NOT.fixm) THEN
        minvec(n_gm)=errvec(n_gm)/hess(n_gm,n_gm)
        DO ii=n_gm-1,1,-1
          minvec(ii)=errvec(ii)
          DO jj=ii+1,n_gm
            minvec(ii)=minvec(ii)-hess(ii,jj)*minvec(jj)
          ENDDO
          minvec(ii)=minvec(ii)/hess(ii,ii)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_min
c-----------------------------------------------------------------------
c     subprogram 11. iter_c3d_hcheck.
c     recompute the Vm^T.A.Vm matrix as a test after performing
c     the standard orthogonalization steps.  this should produce a
c     Hessenberg matrix if the orthogonalization is accurate.  
c-----------------------------------------------------------------------
      SUBROUTINE iter_c3d_hcheck(dot_routine,pre_routine,ctmp,
     $                           dir,nqty,poly_deg,nbdsc,nqdsc,nrb,
     $                           ntotb,nfour)
      USE matrix_mod
      USE factor_type_mod
      USE seam_storage_mod
      USE edge
      USE iter_cg

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: ctmp,dir
      INTEGER(i4), INTENT(IN) :: nqty,poly_deg,nrb,ntotb,nfour
      INTEGER(i4), INTENT(IN) :: nbdsc,nqdsc  !  for discontin. bases

      INTEGER(i4) :: ii,jj,ibl,mxb,myb
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: hnew
      TYPE(cvector_type), DIMENSION(:), ALLOCATABLE :: ctmp2
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dot_routine(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE dot_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE pre_routine(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE pre_routine
      END INTERFACE
c-----------------------------------------------------------------------
c     allocate storage for the inner products, then loop
c     over the exiting basis vectors and compute A.v_jj for each.
c-----------------------------------------------------------------------
      ALLOCATE(hnew(n_gm,n_gm))
 
      ALLOCATE(ctmp2(ntotb))
      DO ibl=1,nrb
        mxb=SIZE(dir(ibl)%arr,2)-1
        myb=SIZE(dir(ibl)%arr,3)-1
        IF (nqdsc>0) THEN
          CALL vector_type_alloc(ctmp2(ibl),poly_deg,mxb,myb,nqty,nfour,
     $                           nbdsc,nqdsc)
        ELSE
          CALL vector_type_alloc(ctmp2(ibl),poly_deg,mxb,myb,nqty,nfour)
        ENDIF
      ENDDO
      DO ibl=nrb+1,ntotb
        mxb=SIZE(dir(ibl)%arr,2)-1
        CALL vector_type_alloc(ctmp2(ibl),1_i4,mxb,0_i4,nqty,nfour)
      ENDDO
c-----------------------------------------------------------------------
c     with left preconditioning, A->P^-1.A, and with right precondition-
c     ing, A->A.P^-1, where P is the approximation of A.
c-----------------------------------------------------------------------

      DO jj=1,n_gm
        IF (leftprec) THEN
          DO ibl=1,ntotb
            dir(ibl)=ep_d(jj)%base(ibl)
          ENDDO
        ELSE
          DO ibl=1,ntotb
            ctmp(ibl)=ep_d(jj)%base(ibl)
          ENDDO
          CALL pre_routine(ctmp,ctmp2,jj,.true.)
          DO ibl=1,ntotb
            dir(ibl)=ctmp2(ibl)
          ENDDO
        ENDIF

        CALL dot_routine(dir,ctmp,.false.)
       
        IF (leftprec) THEN
          DO ibl=1,ntotb
            ctmp2(ibl)=ctmp(ibl)
          ENDDO
          CALL pre_routine(ctmp,ctmp2,1_i4,.false.)
        ENDIF
c-----------------------------------------------------------------------
c       find the inner product of each basis in Vm with A.v_jj.
c-----------------------------------------------------------------------
        DO ii=1,n_gm
          CALL iter_c3d_dot(hnew(ii,jj),ep_d(ii)%base,ctmp,seam,ntotb)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write the values of the inner products to a file.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(a,/)') " Hessenberg matrix Vm^T.A.Vm:"
      DO jj=1,n_gm
        DO ii=1,n_gm
          WRITE(nim_wr,'(2(a,i4),es15.8)')
     $      "   row ",ii,"  col ",jj,hnew(ii,jj)
        ENDDO
      ENDDO
      WRITE(nim_wr,'(/)')

      DO ibl=1,ntotb
        CALL vector_type_dealloc(ctmp2(ibl))
      ENDDO
      DEALLOCATE(ctmp2)
      DEALLOCATE(hnew)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_c3d_hcheck

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_ky_c3d_mod
