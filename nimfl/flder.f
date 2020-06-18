      SUBROUTINE flder(neq,ind,y_lsode,dy)
      USE local

      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: ind
      REAL(r8), INTENT(INOUT), DIMENSION(neq+3) :: y_lsode
      REAL(r8), INTENT(OUT), DIMENSION(neq) :: dy

      REAL(r8), DIMENSION(3) :: b_xyz
      REAL(r8) :: bmag,bigr

      CALL get_bfield(neq,y_lsode,b_xyz,bigr)
      bmag=SQRT(SUM(b_xyz**2))+TINY(bmag)
      dy = b_xyz/bmag
      dy(3) = dy(3)/bigr

      RETURN
      END SUBROUTINE flder
