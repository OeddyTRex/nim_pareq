c-----------------------------------------------------------------------
c     file. nim_stop.f.
c     this is an abbreviated version of what appears in nimrod.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE nim_stop(message)
      USE io
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: message
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(2a)') 'NIM_STOP => ', TRIM(message)
      STOP
      END SUBROUTINE nim_stop
