c module for interface to F77 MPI include file
c defines all MPI variables via parameter statements
c use this module to define machine-specific MPI datatypes
c   so that NIMROD source will not have to change for a new
c   machine or if system.f SELECTED_KINDS are changed

      MODULE mpi_nim

      USE local

      INCLUDE "mpif.h"

      INTEGER(i4), PARAMETER :: mpi_nim_int = mpi_integer
      INTEGER(i4), PARAMETER :: mpi_nim_real = mpi_double_precision
      INTEGER(i4), PARAMETER :: mpi_nim_comp = mpi_double_complex
      INTEGER(i4), PARAMETER :: mpi_nim_logical = mpi_logical
      INTEGER(i4), PARAMETER :: mpi_nim_char = mpi_character


      END MODULE mpi_nim

c     mpi wall-time should be more reliable than other timers when run
c     in parallel, unless the machine is brought down, or jobs are
c     swapped in and out.

c     locating the timer here will not load it preferentially over
c     the one in local.  if it is needed for a particular machine,
c     put this routine into local, but just for nimrod.

c     SUBROUTINE timer(time)
c     USE mpi_nim
c     IMPLICIT NONE

c     REAL(r8), INTENT(INOUT) :: time

c     time=mpi_wtime()
 
c     RETURN
c     END SUBROUTINE timer

