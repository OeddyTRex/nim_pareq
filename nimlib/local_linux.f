c-----------------------------------------------------------------------
c     file local_linux.f
c     module containing defintions of real and integer kinds for linux
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
c     Linux system-dependent timer - returns elapsed CPU seconds
c     the 8-byte integer variables and count_tmp test are used to deal
c     with system-clock resettings and the limited maximum count. 
c-----------------------------------------------------------------------
      SUBROUTINE timer(time)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(INOUT) :: time

      INTEGER(i4), SAVE :: count_rate,count_max,count
      INTEGER(i8), SAVE :: last_count,count_tmp
      LOGICAL, SAVE :: first_call=.true.,warned=.false.

c     f90 intrinsic timer:

      IF (first_call) THEN
        CALL SYSTEM_CLOCK(count=count,count_rate=count_rate,
     $                    count_max=count_max)
        count_tmp=count
        first_call=.false.
      ELSE
        CALL SYSTEM_CLOCK(count=count)
        count_tmp=count
        DO WHILE (count_tmp<last_count-count_max/2)
          count_tmp=count_tmp+count_max
        ENDDO
      ENDIF
      time=REAL(count_tmp)/count_rate
      last_count=count_tmp
      
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     binary open--c90 needs to call assign, Linux doesn't.
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
c     binary close--c90 needs to call assign, Linux doesn't.
c-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CLOSE(UNIT=funit)

      RETURN
      END SUBROUTINE close_bin
c-----------------------------------------------------------------------
c     issue a shell command by calling a local c routine.
c-----------------------------------------------------------------------
      SUBROUTINE system(command)
      USE local

      CHARACTER(*), INTENT(IN) :: command

      CHARACTER(256) :: char_to_c 

      char_to_c=TRIM(command)

      CALL local_system(char_to_c)

      RETURN
      END SUBROUTINE system
