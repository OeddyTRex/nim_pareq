c-----------------------------------------------------------------------
c     file grid.f.
c     Routines that are useful for setting up the grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. grid.
c     1. psigrid.
c     2. qfind.
c     3. pack_run.
c     4. pack_pack_der.
c     5. pack_monitor.
c     6. drpack.
c     7. taugrid.
c     8. vacpack.
c-----------------------------------------------------------------------
c     subprogram 0. grid.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE grid
      USE local
      USE global
      USE input
      IMPLICIT NONE

      LOGICAL, PARAMETER, PRIVATE :: diagnose=.TRUE.
      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: xi,delta

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. psigrid.
c     defines the radial grid.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE psigrid(gt,sq,eqgrid,eqq)
      USE local
      USE spline
      USE input
      IMPLICIT NONE

      TYPE(global_type), INTENT(INOUT) :: gt
      TYPE(spline_type), INTENT(INOUT) :: sq
      REAL(r8), DIMENSION(0:), INTENT(IN) :: eqgrid, eqq

      INTEGER(i4) :: ipsi,m
      REAL(r8) :: dpsi,a,b,sinfac,cosfac

      REAL(r8), PARAMETER :: piby2=pi/2._r8
      LOGICAL, PARAMETER :: diagnose=.FALSE.
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     WRITE formats.
c-----------------------------------------------------------------------
 2010 FORMAT(3x,'mpsi',3x,'psilow',4x,'psihigh'//i6,1p,2e11.3)
 2020 FORMAT(/3x,'ipsi',4x,'psi',8x,'dpsi'/)
 2030 FORMAT(i6,1p,2e11.3)
c-----------------------------------------------------------------------
c     Select grid method
c-----------------------------------------------------------------------
      SELECT CASE(grid_method)
      CASE("original")
         sq%xs=eqgrid
c-----------------------------------------------------------------------
c     uniform grid.
c-----------------------------------------------------------------------
      CASE("uniform")
c         sq%xs=psilow+(psihigh-psilow)*(/(ipsi,ipsi=0,mpsi)/)/mpsi
         sq%xs=psihigh*(/(ipsi,ipsi=1,(mpsi+1))/)/(mpsi+1)
c-----------------------------------------------------------------------
c     stretched grid.
c-----------------------------------------------------------------------
      CASE("stretch")
         m=mpsi+2
         cosfac=COS(pi/(2*m))**2
         sinfac=SIN(pi/(2*m))**2
         a=(psilow*cosfac-psihigh*sinfac)/(cosfac-sinfac)
         b=(psihigh-psilow)/(cosfac-sinfac)
         sq%xs=a+b*SIN(piby2*(/(ipsi,ipsi=1,m-1)/)/m)**2
c-----------------------------------------------------------------------
c     sqrt grid.
c-----------------------------------------------------------------------
      CASE("sqrt")
         sq%xs=psihigh*(/(ipsi**2,ipsi=mxpie,mpsi+mxpie)/)
     $        /(mpsi+mxpie)**2
      CASE("use_gt")
        CALL nim_stop
     &     ("grid: Could not convert grid_type to grid_method")
c-----------------------------------------------------------------------
c     Stop if they haven't specified grid type
c-----------------------------------------------------------------------
      CASE default
         CALL nim_stop("grid: Unrecognized grid_method")
      END SELECT
c-----------------------------------------------------------------------
c     Calculate q stuff, and pack radial coordinate if asked for.
c-----------------------------------------------------------------------
      gt%nn=nn
      CALL qfind(gt,eqgrid,eqq)
      SELECT CASE(pack_method)
      CASE("none")
      CASE("gauss")
        CALL drpack(gt,sq)
      CASE("lorentz")
        CALL pack_run(gt,sq)
      CASE("no_gap")
        CALL pack_run(gt,sq)
      CASE default
        msg="Unknown pack_method "//TRIM(pack_method)
        CALL nim_stop(msg)
      END SELECT
c-----------------------------------------------------------------------
c     diagnose sq%xs array.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(UNIT=ascii_unit,FILE='psi.out',STATUS='UNKNOWN')
         WRITE(ascii_unit,2010)mpsi,psilow,psihigh
         WRITE(ascii_unit,2020)
         dpsi=0
         DO ipsi=0,mpsi
            WRITE(ascii_unit,2030)ipsi,sq%xs(ipsi),dpsi
            dpsi=sq%xs(ipsi+1)-sq%xs(ipsi)
         ENDDO
         WRITE(ascii_unit,2020)
         CLOSE(UNIT=ascii_unit)
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE psigrid
c-----------------------------------------------------------------------
c     subprogram 2. qfind.
c     finds positions of singular values of q.
c     eqgrid = psi_normal, eqq = q on eqgrid
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE qfind(gt,eqgrid,eqq)

      TYPE(global_type), INTENT(INOUT) :: gt
      REAL(r8), DIMENSION(0:), INTENT(IN) :: eqgrid, eqq
      TYPE(spline_type), TARGET :: qq

      INTEGER(i4) :: ipsi,ie,me,i,mnodes,ising,m,dm,it
      INTEGER(i4), PARAMETER :: itmax=200

      REAL(r8) :: b,c,d,dx,x0,x,xmax,
     $     singfac,psifac,dq,psifac0,psifac1
      REAL(r8), DIMENSION(0:30) :: qe,psie
      REAL(r8), DIMENSION(:), POINTER :: psifacs
      REAL(r8), DIMENSION(:,:), POINTER :: sqs
c-----------------------------------------------------------------------
c     Set up spline type
c-----------------------------------------------------------------------
      mnodes=SIZE(eqgrid)-1                             ! 1 b/c 0:mnodes
      CALL spline_alloc(qq,mnodes,1_i4)
      qq%xs = eqgrid
      qq%fs(:,1) = eqq
      CALL spline_fit(qq,"extrap")
c-----------------------------------------------------------------------
c     store left end point.  (e = extremum)
c-----------------------------------------------------------------------
      psifacs => qq%xs
      sqs => qq%fs
      ie=0
      psie(ie)=psifacs(0)
      qe(ie)=sqs(0,1)
c-----------------------------------------------------------------------
c     set up search for extrema of q.
c-----------------------------------------------------------------------
      DO ipsi=0,mnodes-1
         CALL spline_eval(qq,psifacs(ipsi),3_i4)
         xmax=psifacs(ipsi+1)-psifacs(ipsi)
         b=qq%f1(1)
         c=qq%f2(1)
         d=qq%f3(1)
         x0=-c/d
         dx=x0*x0-2*b/d
c-----------------------------------------------------------------------
c     store extremum.
c-----------------------------------------------------------------------
         IF(dx >= 0)THEN
            dx=SQRT(dx)
            DO i=1,2
               x=x0-dx
               IF(x >= 0 .AND. x < xmax)THEN
                  ie=ie+1
                  psie(ie)=psifacs(ipsi)+x
                  CALL spline_eval(qq,psie(ie),0_i4)
                  qe(ie)=qq%f(1)
               ENDIF
               dx=-dx
            ENDDO
         ENDIF
c-----------------------------------------------------------------------
c     complete search and store right end point.
c-----------------------------------------------------------------------
      ENDDO
      ie=ie+1
      psie(ie)=psifacs(mnodes)
      qe(ie)=sqs(mnodes,1)
      me=ie
c-----------------------------------------------------------------------
c     find special q values.
c-----------------------------------------------------------------------
      gt%q0=qq%fs(0,1)-qq%fs1(0,1)*qq%xs(0)
      gt%qmin=min(MINVAL(qe(0:me)),gt%q0)
      gt%qmax=qq%fs(mnodes,1)
      gt%qa=qq%fs(mnodes,1)+qq%fs1(mnodes,1)*(1-qq%xs(mnodes))
c-----------------------------------------------------------------------
c     start loop over extrema to find singular surfaces.
c-----------------------------------------------------------------------
      ising=0
      DO ie=1,me
         m=gt%nn*qe(ie-1)
         dq=qe(ie)-qe(ie-1)
         dm=sign(1._r8,dq)
         IF(dm > 0)m=m+1
c-----------------------------------------------------------------------
c     find singular surfaces by binary search.
c-----------------------------------------------------------------------
         DO
            IF((m-qe(ie-1))*(m-qe(ie)) > 0)EXIT
            it=0
            psifac0=psie(ie-1)
            psifac1=psie(ie)
            DO
               it=it+1
               psifac=(psifac0+psifac1)/2
               CALL spline_eval(qq,psifac,0_i4)
               singfac=(m-gt%nn*qq%f(1))*dm
               IF(singfac > 0)THEN
                  psifac0=psifac
                  psifac=(psifac+psifac1)/2
               ELSE
                  psifac1=psifac
                  psifac=(psifac+psifac0)/2
               ENDIF
               IF(ABS(singfac) <= 1e-12)EXIT
               IF(ABS(psifac-psifac1) <= 1e-12)EXIT
               IF(ABS(psifac-psifac0) <= 1e-12)EXIT
               IF(it > itmax)THEN
                  WRITE(nim_wr,*)"Searching between",psie(ie-1),psie(ie)
                  CALL nim_stop("qfind can't find root")
               ENDIF
            ENDDO
c-----------------------------------------------------------------------
c     store singular surfaces.
c-----------------------------------------------------------------------
            ising=ising+1
            CALL spline_eval(qq,psifac,1_i4)
            gt%msing(ising)=m
            gt%qsing(ising)=REAL(m,r8)/gt%nn
            gt%q1sing(ising)=qq%f1(1)
            gt%psising(ising)=psifac
            m=m+dm
         ENDDO
      ENDDO
      gt%nsing=ising
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE qfind
c-----------------------------------------------------------------------
c     subprogram 3. pack_run.
c     computes grid, function, and curvature.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pack_run(gt,sq)

      TYPE(global_type), INTENT(INOUT) :: gt
      TYPE(spline_type), INTENT(INOUT) :: sq

      INTEGER(i4), PARAMETER :: neq=1,liw=20,lrw=22+neq*16,nstep=10000
      INTEGER(i4) :: ix,iopt,istate,itask,itol,jac,mf,istep,
     $     ising,jsing,nsing
      INTEGER(i4), DIMENSION(liw) :: iwork=0

      REAL(r8), PARAMETER :: tol0=1e-8
      REAL(r8) :: x,xmax,atol,rtol,rnew
      REAL(r8), DIMENSION(neq) :: theta
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8), DIMENSION(0:1) :: f
      REAL(r8), DIMENSION(0:100000,0:1) :: fs
      REAL(r8), DIMENSION(0:mpsi) :: rold

      TYPE(spline_type) :: ff0,ff1, pack
c-----------------------------------------------------------------------
c     Temporary storage of the old grid for diagnostic purposes
c-----------------------------------------------------------------------
      rold(:) = sq%xs(:)
c-----------------------------------------------------------------------
c     format statements.
c-----------------------------------------------------------------------
 10   FORMAT(/5x,'i',6x,'xi',7x,'delta'/)
 20   FORMAT(i6,1p,2e11.3)
 30   FORMAT(/5x,'i',6x,'x',10x,'dx',7x,'theta',8x,'f'/)
 40   FORMAT(i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     count singular surfaces for packing.
c-----------------------------------------------------------------------
      nsing=0
      DO ising=1,gt%nsing
         IF(gt%qsing(ising) >= qsmin .AND. gt%qsing(ising) <= qsmax)
     $        nsing=nsing+1
      ENDDO
      IF(nsing == 0)THEN
       WRITE(nim_wr,*) "No psi packing: no singular surfaces in range"
       RETURN
      ENDIF
c-----------------------------------------------------------------------
c     define monitor function parameters.
c-----------------------------------------------------------------------
      ALLOCATE(xi(nsing))
      ALLOCATE(delta(nsing))
      jsing=0
      DO ising=1,gt%nsing
         IF(gt%qsing(ising) >= qsmin .AND. gt%qsing(ising) <= qsmax)THEN
            jsing=jsing+1
            xi(jsing)=SQRT(gt%psising(ising))
            delta(jsing)=delta0/(ABS(gt%q1sing(ising)))**.25
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     diagnose monitor function parameters.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(UNIT=pack_out_unit,FILE="pack.out",STATUS="UNKNOWN")
         WRITE(pack_out_unit,10)
         WRITE(pack_out_unit,20)(ix,xi(ix),delta(ix),ix=1,SIZE(xi))
         WRITE(pack_out_unit,10)
      ENDIF
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
      istep=0
      x=0
      xmax=SQRT(sq%xs(sq%nodes))
      theta=0
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istate=1
      itask=5
      iopt=1
      mf=10
      itol=1
      rtol=tol0
      atol=tol0
      rwork(1)=xmax
      rwork(11)=0
c-----------------------------------------------------------------------
c     write header.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         CALL open_bin(pack_bin_unit,"monitor.bin","UNKNOWN","REWIND",
     $                 32_i4)
         WRITE(pack_out_unit,30)
      ENDIF
c-----------------------------------------------------------------------
c     store and diagnose output.
c-----------------------------------------------------------------------
      DO
         fs(istep,:)=(/x,theta/)
         IF(diagnose)THEN
            f=pack_monitor(x)
            WRITE(pack_out_unit,40)istep,x,rwork(11),theta,f(0)
            WRITE(pack_bin_unit)REAL(x,4),REAL(theta,4),REAL(f(0),4)
         ENDIF
c-----------------------------------------------------------------------
c     advance.
c-----------------------------------------------------------------------
         IF(x >= xmax .OR. istep >= nstep .OR. istate < 0)EXIT
         istep=istep+1
         CALL dlsode(pack_der,neq,theta,x,xmax,itol,rtol,atol,itask,
     $               istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     write trailer.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         WRITE(pack_out_unit,30)
         WRITE(pack_bin_unit)
         CALL close_bin(pack_bin_unit,"monitor.bin")
      ENDIF
c-----------------------------------------------------------------------
c     fit ff0 to cubic splines.
c-----------------------------------------------------------------------
      CALL spline_alloc(ff0,istep,1_i4)
      ff0%xs=fs(0:istep,1)/fs(istep,1)
      ff0%fs(:,1)=fs(0:istep,0)
      CALL spline_fit(ff0,"extrap")
      ff0%title=(/"theta ","  x   "/)
c-----------------------------------------------------------------------
c     diagnose ff0.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
        CALL open_bin(pack_bin_unit,"ff0.bin","UNKNOWN","REWIND",32_i4)
        WRITE(pack_out_unit,*)"ff0:"
        CALL spline_write(ff0,.TRUE.,.TRUE.,pack_out_unit,pack_bin_unit)
        CALL close_bin(pack_bin_unit,"ff0.bin")
      ENDIF
c-----------------------------------------------------------------------
c     interpolate to vertices and fit ff1 to splines.
c-----------------------------------------------------------------------
      CALL spline_alloc(ff1,mx,2_i4)
      ff1%xs=(/(ix,ix=0,mx)/)/REAL(mx,r8)*ff0%xs(istep)
      DO ix=0,mx
         CALL spline_eval(ff0,ff1%xs(ix),0_i4)
         f=pack_monitor(ff0%f(1))
         ff1%fs(ix,:)=(/ff0%f,f(0)/)
      ENDDO
      CALL spline_fit(ff1,"extrap")
      ff1%title=(/"theta ","  x   ","  f   "/)
c-----------------------------------------------------------------------
c     diagnose ff1.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
        CALL open_bin(pack_bin_unit,"ff1.bin","UNKNOWN","REWIND",32_i4)
        WRITE(pack_out_unit,*)"ff1:"
        CALL spline_write(ff1,.TRUE.,.TRUE.,pack_out_unit,pack_bin_unit)
        CALL close_bin(pack_bin_unit,"ff1.bin")
      ENDIF
c-----------------------------------------------------------------------
c     Spline to new sq%xs which will be used in equilibrium mapping
c-----------------------------------------------------------------------
      CALL spline_alloc(sq,mx,4_i4)
      sq%title=gt%sq%title
      DO ix=0,mx
         CALL spline_eval(ff0,REAL(ix+1,r8)/(mx+1),0_i4)
         sq%xs(ix)=ff0%f(1)**2
      ENDDO
c-----------------------------------------------------------------------
c     Diagnostics - use drawpack.in file
c-----------------------------------------------------------------------
      CALL open_bin(pack_bin_unit,'pack.bin','UNKNOWN','REWIND',32_i4)
        WRITE(pack_bin_unit)REAL(0.,4),
     $                  REAL(rold(0),4),
     $                  REAL(sq%xs(0),4),
     $                  REAL((sq%xs(1) - sq%xs(0))*mx ,4)
      DO ix=1,mx
        WRITE(pack_bin_unit)REAL(REAL(ix)/REAL(mx),4),
     $                  REAL(rold(ix),4),
     $                  REAL(sq%xs(ix),4),
     $                  REAL((sq%xs(ix) - sq%xs(ix-1))*mx ,4)
      ENDDO
      CALL close_bin(pack_bin_unit,'pack.bin')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pack_run
c-----------------------------------------------------------------------
c     subprogram 4. pack_der.
c     computes differential equation.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE pack_der(n,x,theta,dtheta)

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(n), INTENT(in) :: theta
      REAL(r8), DIMENSION(n), INTENT(OUT) :: dtheta

      REAL(r8), DIMENSION(0:1) :: f
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      f=pack_monitor(x)
      dtheta(1)=SQRT(1+f(1)**2)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pack_der
c-----------------------------------------------------------------------
c     subprogram 5. pack_monitor.
c     computes monitor function and its derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION pack_monitor(x) RESULT(f)

      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(0:1) :: f

      INTEGER(i4) :: j
      REAL(r8) :: term0,term1,fac,xfac
c-----------------------------------------------------------------------
c     compute monitor function and derivatives, pack_method="lorentz".
c-----------------------------------------------------------------------
      f=0
      SELECT CASE(pack_method)
      CASE("lorentz")
         DO j=1,SIZE(xi)
            fac=delta(j)**2
            term0=1/((x-xi(j))**2+delta(j)**2)
            term1=-2*(x-xi(j))*term0**2
            f(0)=f(0)+term0*fac
            f(1)=f(1)+term1*fac
         ENDDO
c-----------------------------------------------------------------------
c     compute monitor function and derivatives, pack_method="no_gap".
c-----------------------------------------------------------------------
      CASE("no_gap")
         DO j=1,SIZE(xi)
            fac=delta(j)**2
            xfac=1/(ABS(x-xi(j))+delta(j))
            term0=xfac**2
            term1=-2*xfac*term0
            IF(x < xi(j))term1=-term1
            f(0)=f(0)+term0*fac
            f(1)=f(1)+term1*fac
         ENDDO
c-----------------------------------------------------------------------
c     compute monitor function and derivatives, pack_method="no_gap".
c-----------------------------------------------------------------------
      CASE DEFAULT
         CALL nim_stop
     $        ("Cannot recognize pack_method= "//TRIM(pack_method))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION pack_monitor
c-----------------------------------------------------------------------
c     subprogram 6. drpack.
c     Packs based on specifying the desired dr directly.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE drpack(gt,sq)

      TYPE(global_type), INTENT(IN) :: gt
      TYPE(spline_type), INTENT(INOUT) :: sq

      INTEGER(i4) :: ix, jx, iq, jmax, nsing
      INTEGER(i4), DIMENSION(gt%nsing+1) :: ipack
      REAL(r8) :: rexp, rext, dr0, factor, diff, diffold
      REAL(r8), DIMENSION(0:5*mpsi) :: rr
      REAL(r8), DIMENSION(0:mpsi) :: rgrid, rnew, rold
      TYPE(spline_type) :: grid, pack
c-----------------------------------------------------------------------
c     Temporary storage of the old grid for diagnostic purposes
c-----------------------------------------------------------------------
      rold(:) = sq%xs(:)
c-----------------------------------------------------------------------
c     count singular surfaces for packing.
c-----------------------------------------------------------------------
      nsing=0
      DO iq=1,gt%nsing
         IF(gt%qsing(iq) >= qsmin .AND. gt%qsing(iq) <= qsmax)
     $        nsing=nsing+1
      ENDDO
      IF (vac_pack) nsing=nsing+1
      IF(nsing == 0)THEN
       WRITE(nim_wr,*) "No psi packing: no singular surfaces in range"
       RETURN
      ENDIF
c-----------------------------------------------------------------------
c     Set up rgrid and find where rational surfaces are on grid
c-----------------------------------------------------------------------
      dr0 = 1._r8/REAL(mpsi,r8)                        ! equispaced dr
      diffold = 10.**6
      rgrid(0)=0.
      jx=0
      DO iq=1,gt%nsing
       IF(gt%qsing(iq) >= qsmin .AND. gt%qsing(iq) <= qsmax) THEN
          jx=jx+1
          DO ix=1,mpsi
             rgrid(ix) = REAL(ix,r8)/REAL(mpsi,r8)
             diff = ABS(sq%xs(ix) - gt%psising(iq) )
             IF (diff < diffold) ipack(jx) = ix
             diffold=diff
          ENDDO
       ENDIF
      ENDDO
      IF (vac_pack) ipack(nsing) = mpsi
c-----------------------------------------------------------------------
c     Pack.  To make dr vs. j look like Gaussians, multiply by a shaping
c       Gaussian shaping function.  This will cause the new r (rr) to be
c       on a larger grid (size of jmax).
c-----------------------------------------------------------------------
      rr(0)=0._r8
      factor = 1._r8
      jmax=0
      DO jx=1,5*mpsi
        rext = rr(jx-1) + dr0*factor                    ! Extrapolate r
        factor = 1._r8
        DO iq=1,nsing
          rexp = (rext - rgrid(ipack(iq)))
          factor = factor*(1.-(1. - 1./amp)*exp(-(2.*(rexp/wpack)**2))
     &                    )
        ENDDO
        rr(jx) = rr(jx-1) + dr0*factor
        IF (rr(jx) >= 1._r8) THEN
           jmax = jx
           EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     Rescale rr since rr(jmax) > 1. -- cleans up problems at boundary
c	Gaussian won't be exactly centered, but that isn't important
c-----------------------------------------------------------------------
      IF (jmax == 0) CALL nim_stop("ERRORS: Try decreasing packing")
      DO jx=1,jmax
          rr(jx) = rr(jx)/rr(jmax)
      ENDDO
c-----------------------------------------------------------------------
c     Need to remap to be on the same size grid as the original.
c-----------------------------------------------------------------------
      CALL spline_alloc(grid,jmax,1_i4)
      DO jx=0,jmax
          grid%xs(jx) = REAL(jx,r8)/REAL(jmax,r8)
          grid%fs(jx,1) = rr(jx)
      ENDDO
      CALL spline_fit(grid,"extrap")
      DO ix=0,mpsi
        CALL spline_eval(grid,rgrid(ix),1_i4)
        rnew(ix) = grid%f(1)
      ENDDO
      CALL spline_dealloc(grid)
c-----------------------------------------------------------------------
c     Spline to new sq%xs which will be used in equilibrium mapping
c-----------------------------------------------------------------------
      CALL spline_alloc(pack,mpsi,1_i4)
      pack%xs=rgrid
      pack%fs(:,1) = sq%xs(:)
      CALL spline_fit(pack,"extrap")

      DO ix=0,mpsi
         CALL spline_eval(pack,rnew(ix),1_i4)
         sq%xs(ix)=pack%f(1)
      ENDDO
      CALL spline_dealloc(pack)
c-----------------------------------------------------------------------
c     Diagnostics - use drawpack.in file
c-----------------------------------------------------------------------
      CALL open_bin(pack_bin_unit,'pack.bin','UNKNOWN','REWIND',32_i4)
        WRITE(pack_bin_unit)REAL(rgrid(0),4),
     $                  REAL(rold(0),4),
     $                  REAL(sq%xs(0),4),
     $                  REAL((sq%xs(1) - sq%xs(0))/dr0,4)
      DO ix=1,mpsi
        WRITE(pack_bin_unit)REAL(rgrid(ix),4),
     $                  REAL(rold(ix),4),
     $                  REAL(sq%xs(ix),4),
     $                  REAL((sq%xs(ix) - sq%xs(ix-1))/dr0,4)
      ENDDO
      CALL close_bin(pack_bin_unit,'pack.bin')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE drpack
c-----------------------------------------------------------------------
c     subprogram 7. taugrid.
c     defines the poloidal grid (taus).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE taugrid(thetas,eqarc,taus,mnode)

      INTEGER(i4), INTENT(IN) :: mnode
c      REAL(r8), DIMENSION(0:), INTENT(IN)  :: thetas,eqarc
      REAL(r8), DIMENSION(0:) :: thetas,eqarc
      REAL(r8), DIMENSION(0:), INTENT(OUT) :: taus
      
      INTEGER(i4) :: itau, ilo,k,kp,npts
      INTEGER(i4) :: po=3                              !Polynomial order
      REAL(r8) :: dtau, rtau, yint, dyint
      LOGICAL, PARAMETER :: out=.FALSE.
c-----------------------------------------------------------------------
c     Default is equally-spaced theta
c-----------------------------------------------------------------------
      taus=(/(itau,itau=0,mtheta)/)*twopi/mtheta
c-----------------------------------------------------------------------
c     Alter taus if we want to adapt
c-----------------------------------------------------------------------
      SELECT CASE(angle_pack)
      CASE("arclength")
         ilo=0
         npts=po+1
         DO itau=0,mtheta
            rtau=REAL(itau,r8)*twopi/mtheta
            CALL hunt(eqarc,mnode+1_i4,rtau,ilo)
            k=MIN(MAX(INT(ilo-(po-1)/2),0),mnode-po)
            kp=k+po
            CALL polint(eqarc(k:kp),thetas(k:kp),npts,rtau,yint,dyint)
            taus(itau)=yint
         ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     diagnose taus array.
c-----------------------------------------------------------------------
      IF(out)THEN
         OPEN(UNIT=ascii_unit,FILE='thetas.out',STATUS='UNKNOWN')
         DO itau=0,mnode
            WRITE(ascii_unit,*)itau,eqarc(itau),thetas(itau)
         ENDDO
         CLOSE(UNIT=ascii_unit)
         OPEN(UNIT=ascii_unit,FILE='taus.out',STATUS='UNKNOWN')
         WRITE(ascii_unit,2010)
         dtau=0
         DO itau=0,mtheta
            WRITE(ascii_unit,2020)itau,taus(itau),dtau
            dtau=taus(itau+1)-taus(itau)
         ENDDO
         WRITE(ascii_unit,2010)
         CLOSE(UNIT=ascii_unit)
      ENDIF
c-----------------------------------------------------------------------
c     WRITE formats.
c-----------------------------------------------------------------------
 2010 FORMAT(/3x,'itau',4x,'tau',8x,'dtau'/)
 2020 FORMAT(i6,1p,2e11.3)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE taugrid
c-----------------------------------------------------------------------
c     subprogram 8. vacpack.
c     Sets up packed grid for vacuum region.  Packs around i=0
c     Based on drpack since that is the best packing method ever ;)
c     NOTE (7/17/00): Actually, grids using this are kinda shitty, but 
c      leave it in for now since I don't know what we need.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE vacpack(newgrid,mnode)

      INTEGER(i4), INTENT(IN) :: mnode
      REAL(r8), DIMENSION(0:), INTENT(OUT) :: newgrid

      INTEGER(i4) :: ix, jx, jmax, ipack
      REAL(r8) :: rexp, rext, dr0, factor
      REAL(r8), DIMENSION(0:5*mnode) :: rr
      REAL(r8), DIMENSION(0:mnode) :: rgrid, rnew
      TYPE(spline_type) :: grid, pack
c-----------------------------------------------------------------------
c     Set up rgrid.
c-----------------------------------------------------------------------
      dr0 = 1._r8/REAL(mnode,r8)                        ! equispaced dr
      DO ix=0,mnode
         rgrid(ix) = REAL(ix,r8)*dr0
      ENDDO
      ipack=1
c-----------------------------------------------------------------------
c     Pack.  To make dr vs. j look like Gaussians, multiply by a shaping
c       Gaussian shaping function.  This will cause the new r (rr) to be
c       on a larger grid (size of jmax).
c-----------------------------------------------------------------------
      rr(0)=0._r8
      factor = 1._r8
      jmax=0
      DO jx=1,5*mnode
        rext = rr(jx-1) + dr0*factor                    ! Extrapolate r
        factor = 1._r8
        rexp = (rext - rgrid(ipack))
        factor = factor*(1.-(1. - 1./amp)*exp(-(2.*(rexp/wpack)**2)))
        rr(jx) = rr(jx-1) + dr0*factor
        IF (rr(jx) >= 1._r8) THEN
           jmax = jx
           EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     Rescale rr since rr(jmax) > 1. -- cleans up problems at boundary
c	Gaussian won't be exactly centered, but that isn't important
c-----------------------------------------------------------------------
      DO jx=1,jmax
          rr(jx) = rr(jx)/rr(jmax)
      ENDDO
c-----------------------------------------------------------------------
c     Need to remap to be on the same size grid as the original.
c-----------------------------------------------------------------------
      CALL spline_alloc(grid,jmax,1_i4)
      DO jx=0,jmax
          grid%xs(jx) = REAL(jx,r8)/REAL(jmax,r8)
          grid%fs(jx,1) = rr(jx)
      ENDDO
      CALL spline_fit(grid,"extrap")
      DO ix=0,mnode
        CALL spline_eval(grid,rgrid(ix),1_i4)
        rnew(ix) = grid%f(1)
      ENDDO
      CALL spline_dealloc(grid)
c-----------------------------------------------------------------------
c     Spline to new grid
c-----------------------------------------------------------------------
      CALL spline_alloc(pack,mnode,1_i4)
      pack%xs=rgrid
      pack%fs(:,1) = newgrid(:)
      CALL spline_fit(pack,"extrap")

      DO ix=0,mnode
         CALL spline_eval(pack,rnew(ix),1_i4)
         newgrid(ix)=pack%f(1)
      ENDDO
      CALL spline_dealloc(pack)
c-----------------------------------------------------------------------
c     Diagnostics - use drawpack.in file
c-----------------------------------------------------------------------
      CALL open_bin(pack_bin_unit,'vacpac.bin','UNKNOWN','REWIND',32_i4)
        WRITE(pack_bin_unit)REAL(rgrid(0),4),
     $                  REAL(rgrid(0),4),
     $                  REAL(newgrid(0),4),
     $                  REAL((newgrid(1) - newgrid(0))/dr0,4)
      DO ix=1,mnode
        WRITE(pack_bin_unit)REAL(rgrid(ix),4),
     $                  REAL(rgrid(ix),4),
     $                  REAL(newgrid(ix),4),
     $                  REAL((newgrid(ix) - newgrid(ix-1))/dr0,4)
      ENDDO
      CALL close_bin(pack_bin_unit,'vacpac.bin')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vacpack
c-----------------------------------------------------------------------
c     subprogram 10. polint
c     Polynomial interpolation based on Numerical Recipes.  Used because
c      of problems can occur with spline interpolation.
c-----------------------------------------------------------------------
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(:), INTENT(IN) :: xa,ya
      REAL(r8), INTENT(OUT) :: y, dy
      INTEGER(i4) :: i,m,ns
      REAL(r8) :: den,dif,dift,ho,hp,w,c(10),d(10)
c-----------------------------------------------------------------------
      ns=1
      dif=abs(x-xa(1))
      DO i=1,n
        dift=abs(x-xa(i))
        IF (dift.lt.dif) THEN
          ns=i
          dif=dift
        ENDIF
        c(i)=ya(i)
        d(i)=ya(i)
      ENDDO
      y=ya(ns)
      ns=ns-1
      DO m=1,n-1
        DO i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          IF(den.eq.0.) PAUSE 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        ENDDO   
        IF (2*ns.lt.n-m) THEN
          dy=c(ns+1)
        ELSE
          dy=d(ns)
          ns=ns-1
        ENDIF
        y=y+dy
      ENDDO   
      RETURN
      END SUBROUTINE polint
c-----------------------------------------------------------------------
c     subprogram 10. hunt  (from Numerical Recipes)
c     Given array xx(1:n) and given a value x returns jl0 such that x
c     is given between xx(jlo) and xx(jlo+1).  jlo is also used in initial
c     guess.  This is used so that lower order approximation can be used
c     for polint.
c-----------------------------------------------------------------------
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER(i4), INTENT(IN) :: n
      INTEGER(i4), INTENT(INOUT) :: jlo
      REAL(r8), INTENT(IN) :: x
      REAL(r8), DIMENSION(:), INTENT(IN) :: xx
      INTEGER inc,jhi,jm
      LOGICAL ascnd
c-----------------------------------------------------------------------
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1) RETURN
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END SUBROUTINE hunt
c-----------------------------------------------------------------------
      END MODULE grid
