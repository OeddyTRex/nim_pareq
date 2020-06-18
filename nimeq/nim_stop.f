c-----------------------------------------------------------------------
c     file nim_stop.f
c     closes output files, writes completion message and stops nimeq.
c     this is based on the parallel nim_stop from nimrod.
c-----------------------------------------------------------------------
      SUBROUTINE nim_stop(message)
      USE local
      USE input
      USE input_eq
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: message
      INTEGER(i4) :: ierror
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     for parallel runs, there are lots of cases where nim_stop is
c     called by a subset of the processors.  a less graceful exit is
c     required to stop all processes.
c-----------------------------------------------------------------------
      IF (nprocs>1.AND..NOT.(message=="Normal termination.".OR.
     $                       message=="CPU time limit reached.")) THEN
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        WRITE(out_unit,'(a,i3,2a)') 'NIM_STOP from node ',node,' => ',
     $    TRIM(message)
        CLOSE(UNIT=out_unit)
        WRITE(nim_wr,'(a,i3,2a)') 'NIM_STOP from node ',node,' => ',
     $    TRIM(message)
        CALL mpi_abort(mpi_comm_world,node,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     write completion message.
c-----------------------------------------------------------------------
      IF (node == 0) THEN
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        WRITE(out_unit,'(2a)') 'NIM_STOP => ', TRIM(message)
        WRITE(nim_wr,'(2a)') 'NIM_STOP => ', TRIM(message)
        CLOSE(UNIT=out_unit)
      ENDIF
c-----------------------------------------------------------------------
c     release the processor grid used by SuperLU.
c-----------------------------------------------------------------------
      IF (nimeq_solver(1:5)=='slu_d') THEN
        CALL c_fortran_slugrid(2_i4,comm_layer,node_layer,slu_nrowp,
     $                         slu_ncolp,slugrid_handle)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL mpi_finalize(ierror)
      STOP
      END SUBROUTINE nim_stop
