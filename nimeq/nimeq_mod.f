c-----------------------------------------------------------------------
c     file nimeq_mod.f
c     contains a module for data needed in different parts of the
c     nimeq program.
c-----------------------------------------------------------------------
      MODULE nimeq_mod
      USE local
      USE spline
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     scalar data.
c-----------------------------------------------------------------------
      REAL(r8) :: psio_min  !  minimum value of open (surface) flux
      REAL(r8) :: psio_max  !  maximum value of open (surface) flux
      REAL(r8) :: psic_min  !  minimum value of closed flux
      REAL(r8) :: psic_max  !  maximum value of closed flux
      REAL(r8) :: cflux,oflux,sep_flux
      LOGICAL :: closed_flux,open_flux
c-----------------------------------------------------------------------
c     arrays for equilibrium data from the 1d.bin file (if used).
c-----------------------------------------------------------------------
      TYPE(spline_type) :: neqsq
c-----------------------------------------------------------------------
c     field-lines that integrate to the maximum length (relative length
c     is then 1) are considered closed, and this variable allows a
c     buffer.  with the tracing algorithm, 0.99 has been used.  with the
c     advection-like algorithm, a small value is used.
c-----------------------------------------------------------------------
      REAL(r8) :: rlencl=0.99  !  default for tracing
c-----------------------------------------------------------------------
c     the following are coefficients for Hermite cubic
c     expansions of P and F just inside the last closed surface to 
c     provide a smooth transition.  we assume that derivatives on open
c     flux are 0.
c-----------------------------------------------------------------------
      REAL(r8) :: pvalin,pdrvin,pvalout,pdpsi,
     $            fvalin,fdrvin,fvalout,fdpsi
      LOGICAL :: psmthset=.false.,fsmthset=.false.
c-----------------------------------------------------------------------
c     total number of external coils used in a computation, the
c     current multiplication factor for the feedback coils,
c     and allocatable arrays for all coils are defined.
c-----------------------------------------------------------------------
      INTEGER(i4) :: ncoil_tot=0_i4
      REAL(r8) :: rfbc_curr=0._r8
      REAL(r8) :: vfbc_curr=0._r8
      REAL(r8), DIMENSION(:), ALLOCATABLE :: allcoil_r,allcoil_z,
     $          allcoil_i
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimeq_mod
