c-----------------------------------------------------------------------
c     file local_sgi.f
c     module containing defintions of real and integer kinds for sgi
c     computers.
c-----------------------------------------------------------------------
      MODULE local
      USE io
      IMPLICIT NONE

      INTEGER, PARAMETER ::
     $     ikind1=SELECTED_INT_KIND(2),
     $     ikind2=SELECTED_INT_KIND(4),
     $     i4=SELECTED_INT_KIND(9),
     $     i8=SELECTED_INT_KIND(18),
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi

      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.false.

      END MODULE local

c-----------------------------------------------------------------------
c SGI system-dependent timer - returns elapsed CPU seconds
c-----------------------------------------------------------------------
      SUBROUTINE timer(time)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(INOUT) :: time
      REAL(r4) :: ETIME
      REAL(r4), DIMENSION(2) :: tarray

      time = ETIME(tarray)
      
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     binary open--c90 needs to call assign, sgi doesn't.
c-----------------------------------------------------------------------
      SUBROUTINE open_bin(funit,fname,fstat,fpos,fbit)
      USE local

      CHARACTER(*), INTENT(IN) :: fname,fstat,fpos
      INTEGER, INTENT(IN) :: funit,fbit

      OPEN(UNIT=funit,FILE=fname,STATUS=fstat,POSITION=fpos,
     $     FORM='UNFORMATTED')

      RETURN
      END SUBROUTINE open_bin
c-----------------------------------------------------------------------
c     binary close--c90 needs to call assign, sgi doesn't.
c-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CLOSE(UNIT=funit)

      RETURN
      END SUBROUTINE close_bin
