c--------------------------------------------------------------------
c     file time_stat.f
c     handles machine-dependent timing statistics.
c--------------------------------------------------------------------
c--------------------------------------------------------------------
c     delcarations.
c--------------------------------------------------------------------
      SUBROUTINE time_stat(mode,unit)
      USE local
      IMPLICIT NONE
      
      INTEGER(i4), INTENT(IN) :: mode,unit

      REAL(r8) :: time,time1
      REAL(r8), SAVE :: time0
c--------------------------------------------------------------------
c     get initial time.
c--------------------------------------------------------------------
      IF(mode == 0)THEN
         CALL timer(time0)
c--------------------------------------------------------------------
c     get final time.
c--------------------------------------------------------------------
      ELSE
         CALL timer(time1)
c--------------------------------------------------------------------
c     write elapsed time.
c--------------------------------------------------------------------
         time=time1-time0
         WRITE(unit,'(1X,A,1P,E11.3)')"total cpu time = ",time
         WRITE(nim_wr,'(1X,A,1P,E11.3,A)')"total cpu time = ",
     &            time,CHAR(7)
      ENDIF
c--------------------------------------------------------------------
c     terminate routine.
c--------------------------------------------------------------------
      RETURN
      END SUBROUTINE time_stat
