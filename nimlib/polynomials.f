c-----------------------------------------------------------------------
c     file polynomials.f
c     subprograms for the 1D polynomials used in the finite elements
c     and their numerical integration.  the poly_mod module only
c     contains variables and parameters, no subprograms.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. poly_mod.
c     2. poly_set.
c     3. poly_inquire.
c     4. poly_nodes.
c     5. lagr_1D.
c     6. gauleg.
c     7. lobleg.
c     8. radleg.
c     9. leg_poly1_sub.
c     10. legendre_poly.
c     11. legendre_polyd.
c     12. legendre_eval.
c     13. legendre_deriv.
c-----------------------------------------------------------------------
c     subprogram 1. poly_mod.
c     holds information that influences the location of basis functions
c     nodes.
c-----------------------------------------------------------------------
      MODULE poly_mod
      USE local
      IMPLICIT NONE

c     CHARACTER(7) :: node_dist='uniform'
      CHARACTER(7) :: node_dist='gll'

      INTEGER(i4), PARAMETER :: pd_max=24
      INTEGER(i4), SAVE :: n_last=-1
      CHARACTER(7), SAVE :: dist_last='none'
      REAL(r8), DIMENSION(0:pd_max), SAVE :: x_last
      REAL(r8), DIMENSION(0:pd_max,0:pd_max), SAVE :: cardcoefs

      END MODULE poly_mod

c-----------------------------------------------------------------------
c     subprogram 2. poly_set.
c     provides a mechanism to change subsequent basis node position
c     evaluations during the execution of a program.  
c     if this is used, it should be called before any allocation of
c     "lagr_quad" data structures; otherwise, subsequent evaluations
c     will be erroneous!
c-----------------------------------------------------------------------
      SUBROUTINE poly_set(new_dist)
      USE local
      USE poly_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: new_dist

      node_dist=new_dist

      RETURN
      END SUBROUTINE poly_set

c-----------------------------------------------------------------------
c     subprogram 3. poly_inquire.
c     returns the character variable that determines how nodes are
c     distributed.
c-----------------------------------------------------------------------
      SUBROUTINE poly_inquire(dist)
      USE local
      USE poly_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(OUT) :: dist

      dist=node_dist

      RETURN
      END SUBROUTINE poly_inquire

c-----------------------------------------------------------------------
c     subprogram 4. poly_nodes.
c     finds the n+1 location of nodes for 1D lagrange polynomials in the
c     domain 0<=x<=1.  
c
c     the distribution of nodes is set according to the module variable
c     node_dist.  The value 'uniform' gives uniform nodes, and 'gll'
c     gives a nonuniform distribution.  For the latter,
c     the location of nodes corresponds to the zeros of the 
c     Gauss-Lobatto polynomials, (1+x)(1-x) * d L_n(x)/dx, where
c     L_n is the n-th order Legendre polynomial.  [See "Spectral/hp
c     Element Methods for CFD," Karniadakis and Sherwin, for example.]
c-----------------------------------------------------------------------
      SUBROUTINE poly_nodes(n,x)
      USE local
      USE poly_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), DIMENSION(0:n), INTENT(OUT) :: x

      INTEGER(i4) :: i,j
      REAL(r8) :: wght
      REAL(r8), DIMENSION(0:n) :: wgll
      REAL(r8), EXTERNAL :: legendre_poly
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     don't repeat computation if the nodes haven't changed.
c-----------------------------------------------------------------------
      IF (n==n_last.AND.node_dist==dist_last) THEN
        x=x_last(0:n)
      ELSE
        IF (n>pd_max) THEN
          WRITE(msg,'(a,i3,a)') "POLY_NODES: requested ",n,
     $      " exceeds pd_max. Increase pd_max and recompile."
          CALL nim_stop(TRIM(msg))
        ENDIF
        n_last=n
        dist_last=node_dist
c-----------------------------------------------------------------------
c       uniform distribution.
c-----------------------------------------------------------------------
        IF (node_dist=='uniform') THEN
          DO i=0,n
            x(i)=REAL(i,r8)/REAL(n,r8)
          ENDDO
        ELSE
c-----------------------------------------------------------------------
c         for the Gauss-Lobatto-Legendre nodes, use the lobleg routine
c         and generate the coefficients of the Legendre polynomial
c         for each cardinal function.  note that cardcoefs is ordered
c         with the expansion index first and the cardinal-function
c         index second.
c-----------------------------------------------------------------------
          CALL lobleg(-1._r8,1._r8,x,wgll,n+1_i4)
          DO j=0,n
            wght=0.5_r8*REAL(2_i4*j+1_i4,r8)
            IF (j==n) wght=0.5_r8*REAL(j,r8)
            DO i=0,n
              cardcoefs(j,i)=legendre_poly(x(i),j)*wgll(i)*wght
            ENDDO
          ENDDO
          x=0.5_r8*x+0.5_r8
        ENDIF
        x_last(0:n)=x
      ENDIF

      END SUBROUTINE poly_nodes

c-----------------------------------------------------------------------
c     subprogram 5. lagr_1D.
c     evaluate the general 1D lagrange polynomials at a given point
c     x in the 0<=x<=1 interval.
c
c     note:  the assumed shape arrays al and dal require an interface
c     block in the calling routines.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
      USE local
      USE poly_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: pd,dmode
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal

      INTEGER(i4) :: i,j,k
      REAL(r8), DIMENSION(0:pd) :: c_norm,x_node
      REAL(r8) :: dxtmp

      REAL(r8), EXTERNAL :: legendre_eval,legendre_deriv
c-----------------------------------------------------------------------
c     get the locations of the nodes (zeros of the basis functions).
c-----------------------------------------------------------------------
      CALL poly_nodes(pd,x_node)
c-----------------------------------------------------------------------
c     if the node distribution is GLL, use the coefficients of the
c     Legendre-polynomial expansion of each cardinal function.
c-----------------------------------------------------------------------
      IF (node_dist=='gll') THEN
        DO i=0,pd
          al(i)=legendre_eval(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
        ENDDO
        IF (dmode>0) THEN
          DO i=0,pd
            dal(i)=legendre_deriv(pd,x,0._r8,1._r8,cardcoefs(0:pd,i))
          ENDDO
        ENDIF
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     for other distributions get the normalization constant for the
c     formal cardinal-function relation.
c-----------------------------------------------------------------------
      c_norm=1
      DO i=0,pd
        DO j=0,pd
          IF (j==i) CYCLE
          c_norm(i)=c_norm(i)/(x_node(i)-x_node(j))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute 1D basis values.
c-----------------------------------------------------------------------
      DO i=0,pd
        al(i)=c_norm(i)
        DO j=0,pd
          IF (j==i) CYCLE
          al(i)=al(i)*(x-x_node(j))
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute first derivatives.
c-----------------------------------------------------------------------
      IF (dmode<1) RETURN
      DO i=0,pd
        dal(i)=0
        DO k=0,pd
          IF (k==i) CYCLE
          dxtmp=c_norm(i)
          DO j=0,pd
            IF (j==i) CYCLE
            IF (j==k) CYCLE
            dxtmp=dxtmp*(x-x_node(j))
          ENDDO
          dal(i)=dal(i)+dxtmp
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate execution of lagr_1D.
c-----------------------------------------------------------------------
      RETURN

      CONTAINS

c-----------------------------------------------------------------------
c       factorial function.
c-----------------------------------------------------------------------
        FUNCTION lagr_1D_fac(n) RESULT(nfac)

        INTEGER(i4) :: nfac
        INTEGER(i4), INTENT(IN) :: n

        INTEGER(i4) :: jj

        nfac=1
        DO jj=2,n
          nfac=nfac*jj
        ENDDO

        END FUNCTION lagr_1D_fac

      END SUBROUTINE lagr_1D

c-----------------------------------------------------------------------
c     subprogram 6. gauleg.
c     abscissas and weights for Gauss-Legendre integration, adapted
c     from Numerical Recipies, 2nd ed., Cambridge Press.
c-----------------------------------------------------------------------
      SUBROUTINE gauleg(x1,x2,x,w,n)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w
    
      INTEGER(i4) :: i,j,m
      REAL(r8) :: p1,p2,p3,pp,xl,xm,z,z1
      REAL(r8), PARAMETER :: eps=1.e-14

      m=(n+1)/2
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      DO i=1,m
        z=COS(pi*(REAL(i,r8)-0.25_r8)/(REAL(n,r8)+0.5_r8))
        DO
          p1=1
          p2=0
          DO j=1,n
            p3=p2
            p2=p1
            p1=(REAL(2*j-1,r8)*z*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
          ENDDO
          pp=REAL(n,r8)*(z*p1-p2)/(z**2-1._r8)
          z1=z
          z=z1-p1/pp
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2._r8*xl/((1._r8-z**2)*pp**2)
        w(n+1-i)=w(i)
      ENDDO
    
      RETURN
      END SUBROUTINE gauleg

c-----------------------------------------------------------------------
c     subprogram 7. lobleg.
c     abscissas and weights for Lobatto-Legendre integration, based on
c     Abromowitz and Stegun.
c
c     this version is now independent of poly_nodes and is used by
c     poly_nodes.  its computation has been refined to get rid of the
c     numerical differentiation.
c-----------------------------------------------------------------------
      SUBROUTINE lobleg(x1,x2,x,w,n)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w
    
      INTEGER(i4) :: i,j,m,nleg
      REAL(r8) :: fac,xl,xm,z,z1,p1,p2,p3,gg,dg
      REAL(r8), PARAMETER :: eps=1.e-15
      REAL(r8), DIMENSION(n-1) :: xgau,wgau
c-----------------------------------------------------------------------
c     the zeros of the Legendre polynomial of degree n-1 are used to
c     start Newton's method for finding the zeros of the derivative of
c     the polynomial.
c-----------------------------------------------------------------------
      nleg=n-1_i4
      CALL gauleg(-1._r8,1._r8,xgau,wgau,nleg)
c-----------------------------------------------------------------------
c     the GLL node locations are the zeros of (x**2-1)*d(L_nleg)/dx for
c     the interval -1<=x<=1.  bisect the intervals between zeros of
c     L_nleg to start each Newton iteration and use recurrence to
c     evaluate
c         g(x)=(x**2-1)*d(L_nleg)/dx=nleg*(x*L_nleg-L_(nleg-1))
c     and
c         dg/dx=nleg*(nleg+1)*L_nleg
c-----------------------------------------------------------------------
      m=(n+1_i4)/2_i4
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      fac=2._r8*xl/REAL(n*(n-1_i4),r8)

      w(1)=fac
      x(1)=x1
      w(n)=fac
      x(n)=x2
      DO i=2,m
        z=0.5_r8*(xgau(i)+xgau(i-1))
        DO
          p1=1._r8
          p2=0._r8
          DO j=1,nleg
            p3=p2
            p2=p1
            p1=(REAL(2*j-1,r8)*z*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
          ENDDO
          gg=REAL(nleg,r8)*(z*p1-p2)
          dg=REAL(nleg*(nleg+1_i4),r8)*p1
          z1=z
          z=z1-gg/dg
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        w(i)=fac/p1**2
        w(n+1-i)=w(i)
        x(i)=xm+xl*z
        x(n+1-i)=xm-xl*z
      ENDDO
    
      RETURN
      END SUBROUTINE lobleg

c-----------------------------------------------------------------------
c     subprogram 8. radleg.
c     abscissas and weights for Radau-Legendre integration, based on
c     Abromowitz and Stegun.
c-----------------------------------------------------------------------
      SUBROUTINE radleg(x1,x2,x,w,n)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x1,x2
      REAL(r8), DIMENSION(n), INTENT(OUT) :: x,w
    
      INTEGER(i4) :: i,j,m,nleg
      REAL(r8) :: fac,xl,xm,z,z1,gg,dg,hh,dh,ff,df
      REAL(r8), PARAMETER :: eps=1.e-15
      REAL(r8), DIMENSION(n-1) :: xgau,wgau
      REAL(r8), DIMENSION(n) :: xgaun,wgaun
c-----------------------------------------------------------------------
c     the zeros of the Legendre polynomials of degree n-1 and n are used
c     to start Newton's method for finding the zeros of their sums.
c-----------------------------------------------------------------------
      nleg=n-1_i4
      CALL gauleg(-1._r8,1._r8,xgau,wgau,nleg)
      CALL gauleg(-1._r8,1._r8,xgaun,wgaun,n)
c-----------------------------------------------------------------------
c     the GRL node locations are the zeros of the polynomial
c
c      (L_nleg+L_(nleg+1))/(x+1) = L_nleg+(x-1)*[d(L_nleg)/dx]/(nleg+1)
c
c     the interval -1<=x<=1.  Use 1st-order Taylor expansion about zeros
c     of the two Legendre polynomials on the lhs for an initial guess.
c     to start each Newton iteration.  Use recurrence to
c     evaluate
c         g(x)=L_nleg+(x-1)*[d(L_nleg)/dx]/(nleg+1)
c     and
c         dg/dx={nleg*L_nleg+
c                [(1+x)*(nleg+1)+1-x]*[d(L_nleg)/dx]/(nleg+1)}/(1+x)
c-----------------------------------------------------------------------
      m=(n+1_i4)/2_i4
      xm=0.5_r8*(x2+x1)
      xl=0.5_r8*(x2-x1)
      fac=xl/REAL(n**2,r8)

      w(1)=2._r8*fac
      x(1)=x1
      DO i=2,n
        CALL leg_poly1_sub(xgau(i-1),nleg,gg,dg)
        CALL leg_poly1_sub(xgaun(i),n,ff,df)
        z=(xgau(i-1)*dg+xgaun(i)*df)/(dg+df)
        DO
          CALL leg_poly1_sub(z,nleg,ff,df)
          gg=ff+(z-1._r8)*df/REAL(n,r8)
          dg=(REAL(n,r8)*ff+(z+1._r8)*df-gg)/(z+1._r8)
          z1=z
          z=z1-gg/dg
          IF (ABS(z-z1)<=eps) EXIT
        ENDDO
        w(i)=fac*(1._r8-z)/ff**2
        x(i)=xm+xl*z
      ENDDO
    
      RETURN
      END SUBROUTINE radleg

c-----------------------------------------------------------------------
c     subprogram 9. leg_poly1_sub.
c
c     This subroutine evaluates a single Legendre polynomial and its
c     derivative, using recurrence.  The coding for the value is
c     extracted from gauleg.
c
c     The argument list is:
c
c     xx [real] {input} -- The value of the independent variable in
c        the standard -1 <= xx <= +1 range.
c     nn [integer] {input} -- The index of the polynomial.
c     val [real] {output} -- The value of L_n(xx).
c     drv [real] {output} -- The derivative of L_n(xx).
c-----------------------------------------------------------------------
      SUBROUTINE leg_poly1_sub(xx,nn,val,drv)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: xx
      INTEGER(i4), INTENT(IN) :: nn
      REAL(r8), INTENT(OUT) :: val,drv

      INTEGER(i4) :: j
      REAL(r8) :: p2,p3
      
      val=1._r8
      drv=0._r8
      p2=0._r8
      DO j=1,nn
        p3=p2
        p2=val
        val=(REAL(2*j-1,r8)*xx*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
        drv=xx*drv+REAL(j,r8)*p2
      ENDDO

      RETURN
      END SUBROUTINE leg_poly1_sub

c-----------------------------------------------------------------------
c     subprogram 10. legendre_poly.
c     compute the value of a legendre polynomial of non-negative
c     integer order at a specified position.
c
c     note that the abscissa is expected to be in the standard range
c     of -1 <= x <= +1.
c-----------------------------------------------------------------------
      FUNCTION legendre_poly(x,n) RESULT(ln)
      USE local
      IMPLICIT NONE

      REAL(r8) :: ln
      REAL(r8), INTENT(IN) :: x
      INTEGER(i4), INTENT(IN) :: n

      REAL(r8) :: lim1,lim2
      INTEGER(i4) :: i
c-----------------------------------------------------------------------
c     use standard recursion to generate the value of the desired
c     legendre polynomial.  [Schaum's outline, Mathematical Handbook,
c     by M. Spiegel, for example.]
c-----------------------------------------------------------------------
      IF (n==0) THEN
        ln=1._r8
        RETURN
      ELSE IF (n==1) THEN
        ln=x
        RETURN
      ELSE
        lim2=1._r8
        lim1=x
        DO i=2,n
          ln=(REAL(2*i-1,r8)*x*lim1-REAL(i-1,r8)*lim2)/REAL(i,r8)
          lim2=lim1
          lim1=ln
        ENDDO
      ENDIF
        
      END FUNCTION legendre_poly

c-----------------------------------------------------------------------
c     subprogram 11. legendre_polyd.
c
c     This is a function-version of leg_poly1_sub that returns just the
c     derivative of a Legendre polynomial of non-negative integer order.
c
c     The argument list is:
c
c     xx [real] {input} -- The value of the independent variable in
c        the standard -1 <= xx <= +1 range.
c     nn [integer] {input} -- The index of the polynomial.
c
c     The result of the function is:

c     drv [real]  -- The derivative of L_n(xx).
c-----------------------------------------------------------------------
      FUNCTION legendre_polyd(xx,nn) RESULT(drv)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: xx
      INTEGER(i4), INTENT(IN) :: nn
      REAL(r8) :: drv

      INTEGER(i4) :: j
      REAL(r8) :: p2,p3,val
      
      val=1._r8
      drv=0._r8
      p2=0._r8
      DO j=1,nn
        p3=p2
        p2=val
        val=(REAL(2*j-1,r8)*xx*p2-REAL(j-1,r8)*p3)/REAL(j,r8)
        drv=xx*drv+REAL(j,r8)*p2
      ENDDO

      RETURN
      END FUNCTION legendre_polyd

c-----------------------------------------------------------------------
c     subprogram 12. legendre_eval.
c
c     The function legendre_eval returns the value of a
c     Legendre series approximation with nmax+1 terms.  The routine
c     needs the nmax+1 expansion coefficients, and the independent
c     variable y in the domain ymin <= y <= ymax, which gets mapped
c     to the standard -1 <= x <= 1.
c
c     The argument list for the subroutine is:
c
c     nmax [integer] {input} - Maximum index for the expansion with
c	   (nmax+1) terms for indices 0 through nmax.
c     yy   [real] {input} - The value of the independent variable
c          for the evaluation.
c     ymin [real] {input} - The minimum value of the independent
c          variable for functions on an arbitrary domain.
c     ymax [real] {input} - The maximum value of the independent
c          variable for functions on an arbitrary domain.
c     coef [real(0:nmax)] {input} - Holds the coefficients of the 
c	   Legendre series.
c
c     The real legendre_eval function then evaluates to the sum
c
c	   f(x) = sum_n{ c_n*L_n(x), 0<=n<=nmax }
c
c     where x = x(y) is a linear function of yy.
c
c     This routine is adapted from cheb_eval from cyl_spec.
c-----------------------------------------------------------------------

      FUNCTION legendre_eval(nmax,yy,ymin,ymax,coef) RESULT(approx)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmax
      REAL(r8), INTENT(IN) :: yy,ymin,ymax
      REAL(r8), DIMENSION(0:nmax), INTENT(IN) :: coef
      REAL(r8) :: approx

      REAL(r8), EXTERNAL :: legendre_poly
c-----------------------------------------------------------------------
c     Local variables:
c-----------------------------------------------------------------------

      REAL(r8) :: xx
      INTEGER(i4) :: mm

c-----------------------------------------------------------------------
c     Evaluate the argument for the Legendre polynomials.
c-----------------------------------------------------------------------

      xx=(2._r8*yy-ymin-ymax)/(ymax-ymin)

c-----------------------------------------------------------------------
c     Use the legendre_poly function for each term in the series.
c-----------------------------------------------------------------------

      approx=coef(0)
      IF (nmax>0) approx=approx+coef(1)*xx

      DO mm=2,nmax
        approx=approx+coef(mm)*legendre_poly(xx,mm)
      ENDDO

      RETURN
      END FUNCTION legendre_eval

c-----------------------------------------------------------------------
c     subprogram 13. legendre_deriv.
c
c     The function legendre_deriv returns the derivative of a
c     Legendre series approximation with nmax+1 terms.  The routine
c     needs the nmax+1 expansion coefficients, and the independent
c     variable y in the domain ymin <= y <= ymax, which gets mapped
c     to the standard -1 <= x <= 1.
c
c     The argument list for the subroutine is:
c
c     nmax [integer] {input} - Maximum index for the expansion with
c	   (nmax+1) terms for indices 0 through nmax.
c     yy   [real] {input} - The value of the independent variable
c          for the evaluation.
c     ymin [real] {input} - The minimum value of the independent
c          variable for functions on an arbitrary domain.
c     ymax [real] {input} - The maximum value of the independent
c          variable for functions on an arbitrary domain.
c     coef [real(0:nmax)] {input} - Holds the coefficients of the 
c	   Legendre series.
c
c     The real legendre_eval function then evaluates to the sum
c
c	   f(x) = sum_n{ c_n*L_n(x), 0<=n<=nmax }
c
c     where x = x(y) is a linear function of yy.
c
c     This routine is adapted from cheb_eval from cyl_spec.
c-----------------------------------------------------------------------

      FUNCTION legendre_deriv(nmax,yy,ymin,ymax,coef) RESULT(approx)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmax
      REAL(r8), INTENT(IN) :: yy,ymin,ymax
      REAL(r8), DIMENSION(0:nmax), INTENT(IN) :: coef
      REAL(r8) :: approx

      REAL(r8), EXTERNAL :: legendre_polyd
c-----------------------------------------------------------------------
c     Local variables:
c-----------------------------------------------------------------------

      REAL(r8) :: xx
      INTEGER(i4) :: mm

c-----------------------------------------------------------------------
c     Evaluate the argument for the Legendre polynomials.
c-----------------------------------------------------------------------

      xx=(2._r8*yy-ymin-ymax)/(ymax-ymin)

c-----------------------------------------------------------------------
c     Use the legendre_polyd function for each term in the series.
c-----------------------------------------------------------------------

      approx=0._r8
      IF (nmax>0) approx=approx+coef(1)

      DO mm=2,nmax
        approx=approx+coef(mm)*legendre_polyd(xx,mm)
      ENDDO
      approx=approx*2._r8/(ymax-ymin)

      RETURN
      END FUNCTION legendre_deriv

