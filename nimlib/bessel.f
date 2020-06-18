c-----------------------------------------------------------------------
c     subprogram bessel.
c     computes Bessel functions.
c     Abramowitz & Stegun, Eq. (9.1.10), p. 360.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. bessel_j.
c     2. bessel_y.
c     3. myfact.
c     4. mygamma.
c     5. annular_getco.
c     6. bessel_i
c     7. bessel_k
c     8. bessel_kn
c     9. bessel_in
c-----------------------------------------------------------------------
c     subprogram 1. bessel_j.
c     computes j bessel function of integer order n and real arg z.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bessel_j(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_j

      INTEGER(i4) :: n_abs,myfact,k
      REAL(r8) :: zfac,term
      REAL(r8), PARAMETER :: eps=1e-10_r8
c-----------------------------------------------------------------------
c     compute series.
c-----------------------------------------------------------------------
      n_abs=ABS(n)
      term=(z/2)**REAL(n_abs,r8)/myfact(n_abs)
      zfac=-z*z/4
      bessel_j=term
      k=0
      DO
         k=k+1
         term=term*zfac/(k*(n_abs+k))
         bessel_j=bessel_j+term
         IF(ABS(term).LE.eps*ABS(bessel_j))EXIT
      ENDDO
c-----------------------------------------------------------------------
c     change sign if necessary.
c-----------------------------------------------------------------------
      IF(n.LT.0.AND.MODULO(n,2_i4).EQ.1)bessel_j=-bessel_j
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bessel_j
c-----------------------------------------------------------------------
c     subprogram 2. bessel_y.
c     computes y bessel function of integer order n and real arg z.
c     Abramowitz & Stegun, Eq. (9.1.11), p. 360.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bessel_y(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_y

      INTEGER(i4) :: n_abs,myfact,k,k1
      REAL(r8) :: zfac,term,bessel_j,mypsi
      REAL(r8), PARAMETER :: eps=1e-10_r8
c-----------------------------------------------------------------------
c     compute first series.
c-----------------------------------------------------------------------
      n_abs=ABS(n)
      term=-(z/2)**REAL(n_abs,r8)/(pi*myfact(n_abs))
      zfac=-z*z/4
      bessel_y=term*(mypsi(1)+mypsi(n_abs+1))
      k=0
      DO
         k=k+1
         term=term*zfac/(k*(n_abs+k))
         bessel_y=bessel_y+term*(mypsi(k+1)+mypsi(n_abs+k+1))
         IF(ABS(term).LE.eps*ABS(bessel_y))EXIT
      ENDDO
c-----------------------------------------------------------------------
c     compute second series.
c-----------------------------------------------------------------------
      IF(n_abs.GT.0)THEN
         zfac=-zfac
         k=0
         k1=n_abs
         term=-(z/2)**REAL(-n_abs,r8)*myfact(n_abs-1)/pi
         DO
            bessel_y=bessel_y+term
            k1=k1-1
            IF(k1.EQ.0)EXIT
            k=k+1
            term=term*zfac/(k*k1)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     add J term and change sign if necessary.
c-----------------------------------------------------------------------
      bessel_y=bessel_y+bessel_j(n_abs,z)*LOG(z/2)*2/pi
      IF(n.LT.0.AND.MODULO(n,2_i4).EQ.1)bessel_y=-bessel_y
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bessel_y
c-----------------------------------------------------------------------
c     subprogram 3. myfact.
c     computes n!
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION myfact(n)
      USE local
      IMPLICIT NONE
      
      INTEGER(i4), INTENT(IN) :: n
      INTEGER(i4) :: myfact

      INTEGER(i4) :: i
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(n.GT.0)THEN
         myfact=n
         DO i=n-1,2,-1
            myfact=myfact*i
         ENDDO
      ELSEIF(n.EQ.0)THEN
         myfact=1
      ELSE
         myfact=HUGE(0)
         IF(MOD(n,2_i4).EQ.0)myfact=-myfact
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION myfact
c-----------------------------------------------------------------------
c     subprogram 4. mygamma.
c     computes little gamma function.
c     Abramowitz & Stegun, Eq. (6.3.2), p. 258.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION mypsi(n)
      USE local
      IMPLICIT NONE
      
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8) :: mypsi

      INTEGER(i4) :: k
      REAL(r8), PARAMETER :: gamma=0.57721566490153_r8
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      mypsi=-gamma
      DO k=1,n-1
         mypsi=mypsi+1.0_r8/k
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION mypsi
c-----------------------------------------------------------------------
c     subprogram 5. annular_getco.
c     finds wave number k and coefficients a, b for annular functions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE annular_getco(n,nr,zmin,zmax,k,a,b)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n,nr
      REAL(r8), INTENT(IN) :: zmin,zmax
      REAL(r8), DIMENSION(nr), INTENT(OUT) :: k,a,b

      INTEGER(i4) :: ir
      REAL(r8) :: bessel_j,bessel_y,norm,ki,dk,dk0,f,f1,
     $     j1max,j1min,jmax,jmin,y1max,y1min,ymax,ymin
      REAL(r8), PARAMETER :: eps=1e-9_r8
c-----------------------------------------------------------------------
c     find root wave numbers k.
c-----------------------------------------------------------------------
      ki=0
      dk0=pi/(zmax-zmin)
      DO ir=1,nr
         ki=ki+dk0
         IF(zmin.EQ.0)THEN
            DO
               f=bessel_j(n,ki*zmax)
               f1=f*n/ki-bessel_j(n+1_i4,ki*zmax)*zmax
               dk=-f/f1
               ki=ki+dk
               IF(ABS(dk).lt.eps*ABS(ki))EXIT
            ENDDO
         ELSE
            DO
               jmin=bessel_j(n,ki*zmin)
               jmax=bessel_j(n,ki*zmax)
               ymin=bessel_y(n,ki*zmin)
               ymax=bessel_y(n,ki*zmax)
               j1min=jmin*n/ki-bessel_j(n+1_i4,ki*zmin)*zmin
               j1max=jmax*n/ki-bessel_j(n+1_i4,ki*zmax)*zmax
               y1min=ymin*n/ki-bessel_y(n+1_i4,ki*zmin)*zmin
               y1max=ymax*n/ki-bessel_y(n+1_i4,ki*zmax)*zmax
               f=jmin*ymax-jmax*ymin
               f1=j1min*ymax+jmin*y1max-j1max*ymin-jmax*y1min
               dk=-f/f1
               ki=ki+dk
               IF(ABS(dk).lt.eps*ABS(ki))EXIT
            ENDDO
         ENDIF
         k(ir)=ki
      ENDDO
c-----------------------------------------------------------------------
c     compute coefficients a and b.
c-----------------------------------------------------------------------
      IF(zmin.EQ.0)THEN
         a=1
         b=0
      ELSE
         DO ir=1,nr
            a(ir)=bessel_y(n,k(ir)*zmin)
            b(ir)=-bessel_j(n,k(ir)*zmin)
            norm=sqrt(a(ir)**2+b(ir)**2)
            a(ir)=a(ir)/norm
            b(ir)=b(ir)/norm
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE annular_getco
c-----------------------------------------------------------------------
c     subprogram 6. bessel_i.
c     Abramowitz & Stegun, Eq. (9.6.10), p. 375.
c     Computes I modfied bessel function of integer order n and real arg z.
c-----------------------------------------------------------------------
      FUNCTION bessel_i(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_i

      INTEGER(i4) :: n_abs,myfact,k
      REAL(r8) :: zfac,term
      REAL(r8), PARAMETER :: eps=1e-40_r8 ! (around 8 digit precision)
c-----------------------------------------------------------------------
c     compute series.
c-----------------------------------------------------------------------
      n_abs=ABS(n)
      term=(z/2)**REAL(n_abs,r8)/myfact(n_abs)
      zfac=z*z/4
      bessel_i=term
      k=0
      DO
         k=k+1
         term=term*zfac/(k*(n_abs+k))
         bessel_i=bessel_i+term
         IF(ABS(term).LE.eps*ABS(bessel_i))EXIT
      ENDDO
      RETURN
      END FUNCTION bessel_i

c-----------------------------------------------------------------------
c     subprogram 7. bessel_k.
c     computes K modified bessel function of integer order n and real arg z.
c     Abramowitz & Stegun, Eq. (9.6.11), p. 360.
c     This is only really accurate for n < 10
c-----------------------------------------------------------------------
      FUNCTION bessel_k(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_k

      INTEGER(i4) :: n_abs,myfact,k,k1
      REAL(r8) :: zfac,term,bessel_i,mypsi
      REAL(r8), PARAMETER :: eps=1e-40_r8  ! (about 7 digits of precision)
c-----------------------------------------------------------------------
c     compute first series.
c-----------------------------------------------------------------------
      n_abs=ABS(n)
      term=(-z/2)**REAL(n_abs,r8)/(2*myfact(n_abs))
      zfac=z*z/4
      bessel_k=term*(mypsi(1)+mypsi(n_abs+1))
      k=0
      DO
         k=k+1
         term=term*zfac/(k*(n_abs+k))
         bessel_k=bessel_k+term*(mypsi(k+1)+mypsi(n_abs+k+1))
         IF(ABS(term).LE.eps*ABS(bessel_k))EXIT
      ENDDO
c-----------------------------------------------------------------------
c     compute second series.
c-----------------------------------------------------------------------
      IF(n_abs.GT.0)THEN
         zfac=-zfac
         k=0
         k1=n_abs
         term=(z/2)**REAL(-n_abs,r8)*myfact(n_abs-1)/2
         DO
            bessel_k=bessel_k+term
            k1=k1-1
            IF(k1.EQ.0)EXIT
            k=k+1
            term=term*zfac/(k*k1)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     add I term 
c-----------------------------------------------------------------------
      bessel_k=bessel_k+bessel_i(n_abs,z)*
     $                  LOG(z/2)*(-1)**REAL(n_abs+1,r8)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bessel_k

c-----------------------------------------------------------------------
c     subprogram 8. bessel_kn.
c     computes K modified bessel function of integer order n and real arg z.
c     by recursion from K(m=0) and K(m=1)
c     Abramowitz & Stegun, Eq. (9.6.26), p. 360.
c-----------------------------------------------------------------------
      FUNCTION bessel_kn(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_kn,bessel_km1,bessel_km0,bessel_k
      INTEGER(i4) :: i

      bessel_kn=bessel_k(0_i4,z)
      IF(n == 0)RETURN
      bessel_km1=bessel_kn
      bessel_kn=bessel_k(1_i4,z)
      IF(n == 1)RETURN
      bessel_km0=bessel_kn
      DO i=2,n
        bessel_kn=bessel_km0*2*REAL(i-1)/z+bessel_km1
        bessel_km1=bessel_km0
        bessel_km0=bessel_kn
      ENDDO

      RETURN
      END FUNCTION bessel_kn


c-----------------------------------------------------------------------
c     subprogram 9. bessel_in.
c     computes I modified bessel function of integer order n and real arg z.
c     by recursion from I(m=0) and I(m=1)
c     Abramowitz & Stegun, Eq. (9.6.26), p. 360.
c-----------------------------------------------------------------------
      FUNCTION bessel_in(n,z)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: z
      REAL(r8) :: bessel_in,bessel_im1,bessel_im0,bessel_i
      INTEGER(i4) :: i

      bessel_in=bessel_i(0_i4,z)
      IF(n == 0)RETURN
      bessel_im1=bessel_in
      bessel_in=bessel_i(1_i4,z)
      IF(n == 1)RETURN
      bessel_im0=bessel_in
      DO i=2,n
        bessel_in=-bessel_im0*2*REAL(i-1)/z+bessel_im1
        bessel_im1=bessel_im0
        bessel_im0=bessel_in
      ENDDO

      RETURN
      END FUNCTION bessel_in

