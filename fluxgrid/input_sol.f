c-----------------------------------------------------------------------
c     file input_sol.f
c     input for program sol.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. input_sol.
c     1. read_input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. input_sol.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE input_sol
      USE local
      IMPLICIT NONE

      INTEGER(i4) :: mr=65      ! number of radial grid zones
      INTEGER(i4) :: mz=65      ! number of axial grid zones
      INTEGER(i4) :: ma=64      ! number of flux grid zones
      REAL(r8) :: e=1           ! elongation
      REAL(r8) :: a=1           ! minor radius
      REAL(r8) :: r0=3          ! major radius
      REAL(r8) :: q0=1.26       ! safety factor at the o-point
      LOGICAL :: out=.TRUE.     ! produce ascii output?
      LOGICAL :: bin=.TRUE.     ! produce binary output?

      NAMELIST/all/mr,mz,ma,e,a,r0,q0,out,bin
      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. read_input.
c     reads input.
c-----------------------------------------------------------------------
      SUBROUTINE read_input
      OPEN(UNIT=eq1_unit,FILE='sol.in')
      READ(eq1_unit,NML=all)
      CLOSE(UNIT=eq1_unit)
      END SUBROUTINE read_input
      END MODULE input_sol
