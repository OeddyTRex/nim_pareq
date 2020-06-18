c-----------------------------------------------------------------------
c     file f90_dum_lapack.f
c     stub versions of lapack routines that can be called by nimrod's
c     linear system solvers.  these may be used when the lapack
c     library is not available or one plans to use other linear
c     solvers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. zpbtrf.
c     2. cpbtrf.
c     3. dpbtrf.
c     4. spbtrf.
c     5. zpbtrs.
c     6. cpbtrs.
c     7. dpbtrs.
c     8. spbtrs.
c-----------------------------------------------------------------------
c     subprogram 1. zpbtrf.
c-----------------------------------------------------------------------
      SUBROUTINE zpbtrf(uplo,n,kd,a,lda,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda
      INTEGER, INTENT(OUT) :: ierr 
      COMPLEX :: a

      ierr=-999999

      RETURN
      END SUBROUTINE zpbtrf
c-----------------------------------------------------------------------
c     subprogram 2. cpbtrf.
c-----------------------------------------------------------------------
      SUBROUTINE cpbtrf(uplo,n,kd,a,lda,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda
      INTEGER, INTENT(OUT) :: ierr 
      COMPLEX :: a

      ierr=-999999

      RETURN
      END SUBROUTINE cpbtrf
c-----------------------------------------------------------------------
c     subprogram 3. dpbtrf.
c-----------------------------------------------------------------------
      SUBROUTINE dpbtrf(uplo,n,kd,a,lda,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda
      INTEGER, INTENT(OUT) :: ierr 
      REAL :: a

      ierr=-999999

      RETURN
      END SUBROUTINE dpbtrf
c-----------------------------------------------------------------------
c     subprogram 4. spbtrf.
c-----------------------------------------------------------------------
      SUBROUTINE spbtrf(uplo,n,kd,a,lda,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda
      INTEGER, INTENT(OUT) :: ierr 
      REAL :: a

      ierr=-999999

      RETURN
      END SUBROUTINE spbtrf
c-----------------------------------------------------------------------
c     subprogram 5. zpbtrs.
c-----------------------------------------------------------------------
      SUBROUTINE zpbtrs(uplo,n,kd,nrhs,a,lda,b,ldb,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda,nrhs,ldb
      INTEGER, INTENT(OUT) :: ierr 
      COMPLEX :: a,b

      ierr=-999999

      RETURN
      END SUBROUTINE zpbtrs
c-----------------------------------------------------------------------
c     subprogram 6. cpbtrs.
c-----------------------------------------------------------------------
      SUBROUTINE cpbtrs(uplo,n,kd,nrhs,a,lda,b,ldb,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda,nrhs,ldb
      INTEGER, INTENT(OUT) :: ierr 
      COMPLEX :: a,b

      ierr=-999999

      RETURN
      END SUBROUTINE cpbtrs
c-----------------------------------------------------------------------
c     subprogram 7. dpbtrs.
c-----------------------------------------------------------------------
      SUBROUTINE dpbtrs(uplo,n,kd,nrhs,a,lda,b,ldb,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda,nrhs,ldb
      INTEGER, INTENT(OUT) :: ierr 
      REAL :: a,b

      ierr=-999999

      RETURN
      END SUBROUTINE dpbtrs
c-----------------------------------------------------------------------
c     subprogram 8. spbtrs.
c-----------------------------------------------------------------------
      SUBROUTINE spbtrs(uplo,n,kd,nrhs,a,lda,b,ldb,ierr)
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: uplo
      INTEGER, INTENT(IN) :: n,kd,lda,nrhs,ldb
      INTEGER, INTENT(OUT) :: ierr 
      REAL :: a,b

      ierr=-999999

      RETURN
      END SUBROUTINE spbtrs
