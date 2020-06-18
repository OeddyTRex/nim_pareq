c-----------------------------------------------------------------------
c       this external function evaluates R*A_phi from a set of
c       coils that are modeled as circular loops of wire.  it uses
c       the complete elliptic integral function (cel) that has been
c       adapted from Numerical Recipes.
c-----------------------------------------------------------------------
        FUNCTION aph_eval(r,z,ncoil,coilr,coilz,coili,mu0) RESULT(aph)
        USE local
        IMPLICIT NONE
 
        REAL(r8) :: aph
        REAL(r8), INTENT(IN) :: r,z,mu0
        REAL(r8), DIMENSION(*), INTENT(IN) :: coilr,coilz,coili
        INTEGER(i4), INTENT(IN) :: ncoil
 
        REAL(r8) :: xc2,k2d,k2,efac,kc
        REAL(r8), EXTERNAL :: cel
        REAL(r8), PARAMETER :: aph_tol=1.e-12_r8
        INTEGER(i4), PARAMETER :: aph_maxn=10000
        INTEGER(i4) :: icoil
 
        aph=0._r8
        DO icoil=1,ncoil
          xc2=r**2+(z-coilz(icoil))**2              ! rcyl^2+zcyl^2
          k2d=coilr(icoil)**2+xc2+2*coilr(icoil)*r  ! rm^2
          k2=4*coilr(icoil)*r/k2d                   ! k^2
          IF (k2<=0) THEN
            efac=0.5_r8
          ELSE
            kc=SQRT(1._r8-k2)
            efac=(
     $        (2._r8-k2)*cel(kc,1._r8,1._r8,1._r8,aph_tol,aph_maxn)
     $         -2._r8*cel(kc,1._r8,1._r8,1._r8-k2,aph_tol,aph_maxn))/
     $        (pi*k2)
          ENDIF
          aph=aph+mu0*r*coili(icoil)*coilr(icoil)*efac/SQRT(k2d) ! rAphi
        ENDDO
 
        RETURN
        END FUNCTION aph_eval

c-----------------------------------------------------------------------
c       this external subroutine evaluates the Br and Bz components of
c       magnetic field from coils that are modeled as circular loops.
c-----------------------------------------------------------------------
        SUBROUTINE brz_eval(brz,r,z,ncoil,coilr,coilz,coili,mu0)
        USE local
        IMPLICIT NONE
 
        REAL(r8), DIMENSION(2), INTENT(OUT) :: brz
        REAL(r8), INTENT(IN) :: r,z,mu0
        REAL(r8), DIMENSION(*), INTENT(IN) :: coilr,coilz,coili
        INTEGER(i4), INTENT(IN) :: ncoil
 
        REAL(r8) :: xc2,k2d,k2,efac,kc,cele,celk,muotp,dz,cr2
        REAL(r8), EXTERNAL :: cel
        REAL(r8), PARAMETER :: aph_tol=1.e-12_r8
        INTEGER(i4), PARAMETER :: aph_maxn=10000
        INTEGER(i4) :: icoil

        brz(1:2)=0._r8
        muotp=mu0/twopi
        DO icoil=1,ncoil
          dz=z-coilz(icoil)
          cr2=coilr(icoil)**2
          xc2=r**2+dz**2                            ! rcyl^2+zcyl^2
          k2d=coilr(icoil)**2+xc2+2*coilr(icoil)*r  ! rm^2
          k2=4*coilr(icoil)*r/k2d                   ! k^2
          IF (k2<=0) THEN
            brz(1:2)=brz(1:2)+
     $               (/0._r8,-0.5_r8*mu0*coili(icoil)*coilr(icoil)**2/
     $                       (coilr(icoil)**2+dz**2)**1.5_r8/)
          ELSE
            kc=SQRT(1._r8-k2)
            celk=cel(kc,1._r8,1._r8,1._r8,aph_tol,aph_maxn)
            cele=cel(kc,1._r8,1._r8,1._r8-k2,aph_tol,aph_maxn)
            brz(1:2)=brz(1:2)-muotp*coili(icoil)/SQRT(k2d)*(/
     $        dz*((xc2+cr2)*cele/(xc2+cr2-2._r8*r*coilr(icoil))-celk)/r,
     $        (cr2-xc2)*cele/(xc2+cr2-2._r8*r*coilr(icoil))+celk /)
          ENDIF
        ENDDO
 
        RETURN
        END SUBROUTINE brz_eval
