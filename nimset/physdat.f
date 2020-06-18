c-----------------------------------------------------------------------
c     file physdat.f
c     module containing the physical parameters used in the 2-fluid
c     equations.
c-----------------------------------------------------------------------
      MODULE physdat
      USE local
      IMPLICIT NONE

      REAL(r8) :: elementary_q=1.60217733e-19_r8
      REAL(r8), DIMENSION(2) :: 
     $          qs=(/-1.60217733e-19_r8,1.60217733e-19_r8 /),
     $          ms=(/ 9.1093898e-31_r8, 3.3435860e-27_r8 /)
      REAL(r8) :: clight=2.99792458e8_r8
      REAL(r8) :: mu0=pi*4.e-7_r8 
      REAL(r8) :: gamma=5._r8/3._r8
      REAL(r8) :: kboltz=1.60217733e-19_r8
      REAL(r8) :: zeff,mtot,meomi,eps0,gamm1

      CONTAINS
c-----------------------------------------------------------------------
c     routine for setting parameters.  this must be called regardless
c     of whether constant are specified in the namelist.
c-----------------------------------------------------------------------
      SUBROUTINE physdat_set(chrg_inp,zeff_inp,mi_inp,me_inp,gam_inp,
     $                       kblz_inp,mu0_inp,c_inp)

      REAL(r8), INTENT(IN), OPTIONAL :: chrg_inp,zeff_inp,mi_inp,me_inp,
     $                                  gam_inp,kblz_inp,mu0_inp,c_inp 
c-----------------------------------------------------------------------
c     alter any parameters present in the call.
c-----------------------------------------------------------------------
      IF (PRESENT(chrg_inp)) THEN
        elementary_q=chrg_inp
        qs(1)=-elementary_q
        qs(2)= elementary_q
      ENDIF
      IF (PRESENT(zeff_inp)) qs(2)=zeff_inp*qs(2)
      IF (PRESENT(me_inp)) ms(1)=me_inp
      IF (PRESENT(mi_inp)) ms(2)=mi_inp
      IF (PRESENT(gam_inp)) gamma=gam_inp
      IF (PRESENT(kblz_inp)) kboltz=kblz_inp
      IF (PRESENT(mu0_inp)) mu0=mu0_inp
      IF (PRESENT(c_inp)) clight=c_inp
c-----------------------------------------------------------------------
c     compute secondary parameters.
c-----------------------------------------------------------------------
      zeff=-qs(2)/qs(1)
      mtot=ms(1)+ms(2)/zeff
      meomi=zeff*ms(1)/ms(2)
      eps0=1._r8/(clight**2*mu0)
      gamm1=gamma-1._r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE physdat_set
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE physdat
