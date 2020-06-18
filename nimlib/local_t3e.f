c-----------------------------------------------------------------------
c     file local_t3e.f
c     module containing defintions of real and integer kinds.
c-----------------------------------------------------------------------
      MODULE local
      USE io
      IMPLICIT NONE

c SJP start - change 9 -> 18 for T3E
      INTEGER, PARAMETER ::
     $     ikind1=SELECTED_INT_KIND(2),
     $     ikind2=SELECTED_INT_KIND(4),
     $     i4=SELECTED_INT_KIND(18),
     $     i8=SELECTED_INT_KIND(18),
     $     r4=SELECTED_REAL_KIND(6,37),
     $     r8=SELECTED_REAL_KIND(13,307)
c SJP end
      REAL(r8), PARAMETER :: pi=3.1415926535897932385_r8,
     $     twopi=2*pi,pisq=pi*pi

      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.true.

      END MODULE local

c-----------------------------------------------------------------------
c Cray T3E system-dependent timer - returns elapsed CPU seconds
c-----------------------------------------------------------------------
      SUBROUTINE timer(time)
      USE local
      IMPLICIT NONE
      REAL(r8), INTENT(INOUT) :: time
      
      INCLUDE "mpif.h"

c     time = rtc() * 2.22E-9

c     time = mpi_wtime()
      
      CALL SECOND(time)

      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     binary open--t3e needs to call assign to generate standard data
c     formats:  plot files are written as 32 bit, but dump files are 
c     written as 64 bit.
c-----------------------------------------------------------------------
      SUBROUTINE open_bin(funit,fname,fstat,fpos,fbit)
      USE local

      CHARACTER(*), INTENT(IN) :: fname,fstat,fpos
      INTEGER, INTENT(IN) :: funit,fbit

      CHARACTER(256) :: assign_instr,msg
      CHARACTER(229) :: path
      CHARACTER(128) :: directory
      INTEGER(i4) :: ierror
      INTEGER :: ipxferr,ipxflen

c-----------------------------------------------------------------------
c     assign standard format to the full path name of the binary file.
c     do not call nim_stop if there is an error--this leads to an
c     infinite loop.
c-----------------------------------------------------------------------
      path=ADJUSTL(fname)
      IF (fbit==32.OR.fbit==64) THEN
        IF (path(1:1)/='/') THEN
          CALL PXFGETCWD(directory,ipxflen,ipxferr)
          IF (ipxferr/=0) THEN
            WRITE(msg,'(3a)') 'Open_bin unable to acquire path of ',
     $           TRIM(ADJUSTL(fname)),'.'
            OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $           POSITION='APPEND')
            WRITE(out_unit,'(2a)') 'OPEN_BIN => ', TRIM(msg)
            WRITE(nim_wr,'(2a)') 'OPEN_BIN => ', TRIM(msg)
            CLOSE(UNIT=out_unit)
            STOP
          ENDIF
c         CALL system('pwd > temporary_path_name')
c         OPEN(UNIT=temp_unit,FILE='temporary_path_name')
c         READ(temp_unit,'(a)') directory
c         CLOSE(temp_unit)
c         CALL system('rm temporary_path_name')
          path=TRIM(ADJUSTL(directory))//'/'//path
        ENDIF
        WRITE(assign_instr,'(2a)') 'assign -F f77 f:',TRIM(path)
        CALL assign(TRIM(assign_instr),ierror)
        IF (ierror/=0) THEN
          WRITE(msg,'(3a)') 'Open_bin unable to assign file ',
     $         TRIM(ADJUSTL(fname)),'.'
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          WRITE(out_unit,'(2a)') 'OPEN_BIN => ', TRIM(msg)
          WRITE(nim_wr,'(2a)') 'OPEN_BIN => ', TRIM(msg)
          CLOSE(UNIT=out_unit)
          STOP
        ENDIF
      ENDIF

      OPEN(UNIT=funit,FILE=TRIM(path),STATUS=fstat,POSITION=fpos,
     $     FORM='UNFORMATTED')

      RETURN
      END SUBROUTINE open_bin
c-----------------------------------------------------------------------
c     binary close--t3e needs to call assign.
c-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CHARACTER(256) :: assign_instr,msg
      CHARACTER(229) :: path
      CHARACTER(128) :: directory
      INTEGER(i4) :: ierror

      CLOSE(UNIT=funit)

      path=ADJUSTL(fname)
      IF (path(1:1)/='/') THEN
        CALL PXFGETCWD(directory,ipxflen,ipxferr)
        IF (ipxferr/=0) THEN
          WRITE(msg,'(3a)') 'Close_bin unable to acquire path of ',
     $         TRIM(ADJUSTL(fname)),'!!!'
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          WRITE(out_unit,'(2a)') 'CLOSE_BIN => ', TRIM(msg)
          WRITE(nim_wr,'(2a)') 'CLOSE_BIN => ', TRIM(msg)
          CLOSE(UNIT=out_unit)
        ENDIF
c       CALL system('pwd > temporary_path_name')
c       OPEN(UNIT=temp_unit,FILE='temporary_path_name')
c       READ(temp_unit,'(a)') directory
c       CLOSE(temp_unit)
c       CALL system('rm temporary_path_name')
        path=TRIM(ADJUSTL(directory))//'/'//path
      ENDIF
      WRITE(assign_instr,'(2a)') 'assign -R f:',TRIM(path)
      CALL assign(TRIM(assign_instr),ierror)

      IF (ierror/=0) THEN
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        WRITE(out_unit,'(3a)') 'Close_bin unable to reassign file ',
     $       TRIM(ADJUSTL(fname)),'!!!'
        CLOSE(UNIT=out_unit)
        WRITE(nim_wr,'(3a)') 'Close_bin unable to reassign file ',
     $       TRIM(ADJUSTL(fname)),'!!!'
      ENDIF

      RETURN
      END SUBROUTINE close_bin
c-----------------------------------------------------------------------
c     issue a shell command--t3e does not have the semi-standard
c     system command.
c-----------------------------------------------------------------------
      SUBROUTINE system(command)
      USE local

      CHARACTER(*), INTENT(IN) :: command

      INTEGER :: ierror
      CHARACTER(128) :: err_msg

      ierror=ISHELL(TRIM(command))

      IF (ierror<0) THEN
        WRITE(err_msg,'(a,i8,2a)') 'ISHELL error # ',ierror,' from ',
     $    TRIM(command)
        CALL nim_stop(TRIM(err_msg))
      ENDIF

      RETURN
      END SUBROUTINE system
