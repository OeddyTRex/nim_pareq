c a dummy timing module for nimplot.
c this helps us use nimrod modules without change.

      MODULE time
      USE local
      IMPLICIT NONE

c cummulative times and counters

      INTEGER(i4) :: seamcount = 0
      INTEGER(i4) :: segcount = 0
      INTEGER(i4) :: seamcount_loop,segcount_loop

      REAL(r8) :: time_io = 0.0
      REAL(r8) :: time_seam = 0.0
      REAL(r8) :: time_seg = 0.0
      REAL(r8) :: time_iter = 0.0
      REAL(r8) :: time_fac = 0.0
      REAL(r8) :: time_fft = 0.0
      REAL(r8) :: time_mat = 0.0
      REAL(r8) :: time_rhs = 0.0
      REAL(r8) :: time_line = 0.0
      REAL(r8) :: time_stcon = 0.0
      REAL(r8) :: time_io_loop,time_seam_loop,time_seg_loop
      REAL(r8) :: time_iter_loop,time_fft_loop,time_fe_loop
      REAL(r8) :: time_fac_loop

      REAL(r8) :: time_total_start,time_total_end
      REAL(r8) :: time_loop_start,time_loop_end
      REAL(r8) :: timestart,timeend

      END MODULE time
