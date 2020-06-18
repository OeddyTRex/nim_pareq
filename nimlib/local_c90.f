c-----------------------------------------------------------------------
c     file local_c90.f
c     module containing defintions of real and integer kinds for the
c     cray c90.
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

      LOGICAL, PARAMETER :: rewind_namel=.false.,single_pr=.true.

      END MODULE local

c-----------------------------------------------------------------------
c Dummy double precision LaPack routines for the c90.
c-----------------------------------------------------------------------
      SUBROUTINE dpotrs(c,ii1,ii2,rr1,ii3,rr2,ii4,ii5)
      USE local
      INTEGER :: ii1,ii2,ii3,ii4,ii5
      CHARACTER(*) :: c
      REAL(r8), DIMENSION(:,:,:,:) :: rr1
      REAL(r8), DIMENSION(:,:,:) :: rr2
      RETURN
      END SUBROUTINE dpotrs

      SUBROUTINE dpotrf(c,ii1,rr1,ii2,ii3)
      USE local
      INTEGER :: ii1,ii2,ii3
      CHARACTER(*) :: c
      REAL(r8), DIMENSION(:,:,:,:) :: rr1
      RETURN
      END SUBROUTINE dpotrf

      SUBROUTINE dpbtrs(c,ii1,ii2,ii3,rr1,ii4,rr2,ii5,ii6)
      USE local
      INTEGER :: ii1,ii2,ii3,ii4,ii5,ii6
      CHARACTER(*) :: c
      REAL(r8), DIMENSION(:,:) :: rr1
      REAL(r8), DIMENSION(:,:,:,:,:) :: rr2
      RETURN
      END SUBROUTINE dpbtrs

      SUBROUTINE dpbtrf(c,ii1,ii2,rr1,ii3,ii4)
      USE local
      INTEGER :: ii1,ii2,ii3,ii4
      CHARACTER(*) :: c
      REAL(r8), DIMENSION(:,:) :: rr1
      RETURN
      END SUBROUTINE dpbtrf

      SUBROUTINE ddot(ii1,rr1,ii2,rr2,ii3)
      USE local
      INTEGER :: ii1,ii2,ii3
      REAL(r8), DIMENSION(:,:,:) :: rr1,rr2
      RETURN
      END SUBROUTINE ddot

c-----------------------------------------------------------------------
c Cray C90 system-dependent timer - returns elapsed CPU seconds
c-----------------------------------------------------------------------
      SUBROUTINE timer(time)
      USE local
      IMPLICIT NONE
      REAL(r8), INTENT(INOUT) :: time
      REAL :: SECOND
      
      time = SECOND()
      
      RETURN
      END SUBROUTINE timer
c-----------------------------------------------------------------------
c     binary open--c90 needs to call assign to generate standard data
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

c-----------------------------------------------------------------------
c     assign standard format to the full path name of the binary file.
c     do not call nim_stop if there is an error--this leads to an
c     infinite loop.
c-----------------------------------------------------------------------
      path=ADJUSTL(fname)
      IF (fbit==32.OR.fbit==64) THEN
        IF (path(1:1)/='/') THEN
          CALL system('pwd > temporary_path_name')
          OPEN(UNIT=temp_unit,FILE='temporary_path_name')
          READ(temp_unit,'(a)') directory
          CLOSE(temp_unit)
          CALL system('rm temporary_path_name')
          path=TRIM(ADJUSTL(directory))//'/'//path
        ENDIF
        WRITE(assign_instr,'(a,i2,2a)') 'assign -F f77 -N ieee_',fbit,
     $                                  ' f:',TRIM(path)
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
c     binary close--c90 needs to call assign.
c-----------------------------------------------------------------------
      SUBROUTINE close_bin(funit,fname)
      USE local

      CHARACTER(*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: funit

      CHARACTER(256) :: assign_instr,msg
      CHARACTER(244) :: path
      CHARACTER(128) :: directory
      INTEGER(i4) :: ierror

      CLOSE(UNIT=funit)

      path=ADJUSTL(fname)
      IF (path(1:1)/='/') THEN
        CALL system('pwd > temporary_path_name')
        OPEN(UNIT=temp_unit,FILE='temporary_path_name')
        READ(temp_unit,'(a)') directory
        CLOSE(temp_unit)
        CALL system('rm temporary_path_name')
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
c     issue a shell command--c90 does not have the semi-standard
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
