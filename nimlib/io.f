c-----------------------------------------------------------------------
c     file io.f
c     module containing fortran unit numbers for io.
c-----------------------------------------------------------------------
      MODULE io
      IMPLICIT NONE

      INTEGER, PARAMETER :: in_unit = 1		! input file unit
      INTEGER, PARAMETER :: out_unit = 2	! output file unit
      INTEGER, PARAMETER :: xdr_unit = 3	! unit for xdraw output
      INTEGER, PARAMETER :: nim_rd = 5		! stdin unit for reads
      INTEGER, PARAMETER :: nim_wr = 6		! stdout unit for writes
      INTEGER, PARAMETER :: dx1_unit = 8	! nimrod.dx unit
      INTEGER, PARAMETER :: dx2_unit = 9	! nimrod.bin unit
      INTEGER, PARAMETER :: hst_unit = 11	! unit for nimhist.bin
      INTEGER, PARAMETER :: it_unit = 12	! unit for iter.out
      INTEGER, PARAMETER :: en_unit = 13	! unit for energy.bin
      INTEGER, PARAMETER :: dis_unit = 14	! unit for discharge.bin
      INTEGER, PARAMETER :: eq_unit = 16	! equilibrium file
      INTEGER, PARAMETER :: grd_unit = 17	! unit for grid check
      INTEGER, PARAMETER :: tec2d = 18          ! tecplot output unit
      INTEGER, PARAMETER :: tec1d=19
      INTEGER, PARAMETER :: dump_unit = 25	! dump file unit
      INTEGER, PARAMETER :: rstrt_unit = 26	! restart file unit
      INTEGER, PARAMETER :: pie_unit = 27	! pie file unit
      INTEGER, PARAMETER :: rim_unit = 28	! rim file unit
      INTEGER, PARAMETER :: xy_unit = 31 	! xy binary output
      INTEGER, PARAMETER :: xt_unit = 32 	! xt binary output
      INTEGER, PARAMETER :: yt_unit = 33 	! yt binary output
      INTEGER, PARAMETER :: grid_unit = 34 	! grid binary output
      INTEGER, PARAMETER :: con_unit = 35 	! contour plot output
      INTEGER, PARAMETER :: ascii_unit = 41	! ascii diagnostic file
      INTEGER, PARAMETER :: binary_unit = 42	! binary diagnostic file
      INTEGER, PARAMETER :: eq1_unit=43
      INTEGER, PARAMETER :: eq2_unit=44
      INTEGER, PARAMETER :: pack_ascii_unit=45
      INTEGER, PARAMETER :: pack_bin_unit=46
      INTEGER, PARAMETER :: pack_in_unit=47
      INTEGER, PARAMETER :: pack_out_unit=48
      INTEGER, PARAMETER :: pack_spline_unit=49
      INTEGER, PARAMETER :: temp_unit = 91 	! temporary files

      LOGICAL, SAVE :: out_opened  !  used by nimrod

      END MODULE io
