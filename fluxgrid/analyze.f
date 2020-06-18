c     file analyze.f.
c     diagnosis equilibrium
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. analyze.
c     1. get_global.
c     2. get_stability.
c     3. fluxav.
c     4. qsearch.
c     5. qres.
c-----------------------------------------------------------------------
c     subprogram 0. analyze.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE analyze
      USE local
      USE global
      USE input
      USE physdat
      USE bicube
      IMPLICIT NONE

      REAL(r8) :: tauafac,taurfac,taua,taur,elecd,volume
      REAL(r8), DIMENSION(:),POINTER :: dideal, dres, hfactor
      REAL(r8), DIMENSION(:),POINTER :: dnc, alphas, fcirc
      REAL(r8), DIMENSION(:),POINTER :: lambdaprof
      REAL(r8), DIMENSION(:,:),POINTER :: lambda, delstr,gsrhs

      INTEGER(i4) :: num_root
      REAL(r8), DIMENSION(:), POINTER :: psi_root
      REAL(r8) :: q_search_value
      TYPE(spline_type) :: qprof

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. get_global.
c     computes global equilibrium parameters.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE get_global(gt)

      TYPE(global_type), INTENT(INOUT) :: gt

      INTEGER(i4) :: ipsi, iqty,itheta,iside
      TYPE(spline_type) :: gs,hs
      REAL(r8) :: dpsi,rfac,r,z,cosfac,eta,eta1,eta2,phase,
     $     phase1,phase2,r2,r21,r22,rfac1,rfac2,sinfac,z1,z2,bp0,
     $     zold,theta,dtheta
      REAL(r8) :: gss, gst, bsq, jdotb, ldelstr
      REAL(r8) :: p, pprime, f, fprime,psio,jacfac, mach, msprime
      REAL(r8), DIMENSION(2) :: rext,zext,thetaext
      REAL(r8), DIMENSION(0:mtheta) :: thetas
      REAL(r8), DIMENSION(0:mpsi) :: psifacs
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: g11,jac
      LOGICAL, PARAMETER :: diagnose=.FALSE.
      REAL(r8), PARAMETER :: mi=3.3435860e-27_r8
c-----------------------------------------------------------------------
c     WRITE formats.
c-----------------------------------------------------------------------
 2010 FORMAT(1x,'ipsi = ',i3,', psi = ',1p,e11.3)
 2015 FORMAT(1x,'ipsi = ',i3,', jpsi = ',i1,', psi = ',1p,e11.3)
 2020 FORMAT(/4x,'it',5x,'tau',9x,'r2',8x,'deta',7x,'eta',9x,'r',10x,
     $     'z'/)
 2030 FORMAT(i6,1p,6e11.3)

c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      thetas=gt%twod%xs
      psifacs=gt%twod%ys
      CALL spline_alloc(gs,mtheta,2_i4)
      jac(:,:)=gt%twod%fs(3,:,:)*gt%psio
      g11(:,:)= gt%twod%fs(4,:,:)/gt%psio**2
c-----------------------------------------------------------------------
c     compute poloidal integrands at plasma surface.
c-----------------------------------------------------------------------
      DO itheta=0,mtheta
         r=gt%twod%fs(1,itheta,mpsi)
         gs%fs(itheta,1)=jac(itheta,mpsi)*SQRT(g11(itheta,mpsi))/r
         gs%fs(itheta,2)=jac(itheta,mpsi)*g11(itheta,mpsi)/r**2
      ENDDO
      gs%xs=thetas
      CALL spline_fit(gs,"periodic")
      CALL spline_int(gs)
      gt%crnt=gt%psio*gs%fsi(mtheta,2)/(mu0*1.e6)
      bp0=gt%psio*gs%fsi(mtheta,2)/gs%fsi(mtheta,1)
c-----------------------------------------------------------------------
c     diagnose gs.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(UNIT=ascii_unit,FILE='gs.out')
         CALL open_bin(binary_unit,"gs.bin","UNKNOWN","REWIND",32_i4)
         WRITE(ascii_unit,2010)(iqty,iqty,iqty,iqty=1,2)
         WRITE(ascii_unit,2020)(itheta,thetas(itheta),
     $        (gs%fs(itheta,iqty),gs%fs1(itheta,iqty),
     $        gs%fsi(itheta,iqty),iqty=1,2),itheta=0,mtheta)
         WRITE(ascii_unit,2010)(iqty,iqty,iqty,iqty=1,2)
         gs%title=(/"theta ","  g1  ","  g2  "/)
         CALL spline_write(gs,.TRUE.,.TRUE.,ascii_unit,binary_unit)
         CLOSE(UNIT=ascii_unit)
         CALL close_bin(binary_unit,"gs.bin")
      ENDIF
c-----------------------------------------------------------------------
c     integrate poloidal magnetic field over each flux surface.
c-----------------------------------------------------------------------
      CALL spline_alloc(hs,mpsi,4_i4)
      hs%xs=psifacs
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            r=gt%twod%fs(1,itheta,ipsi)
            gs%fs(itheta,1)=jac(itheta,ipsi)
            gs%fs(itheta,2)=jac(itheta,ipsi)*g11(itheta,ipsi)/r**2
         ENDDO
         CALL spline_fit(gs,"periodic")
         CALL spline_int(gs)
         hs%fs(ipsi,3)=gs%fsi(mtheta,2)*gt%psio**2
         hs%fs(ipsi,4)=gs%fsi(mtheta,1)*twopi
      ENDDO
c-----------------------------------------------------------------------
c     integrate flux functions.
c-----------------------------------------------------------------------
      dpsi=1-psifacs(mpsi)
      hs%fs(:,1)=gt%sq%fs(:,2)*hs%fs(:,4)
      hs%fs(:,2)=hs%fs(:,4)
      CALL spline_fit(hs,"extrap")
      CALL spline_int(hs)
      hs%fsi(mpsi,:)=hs%fsi(mpsi,:)+(hs%fs(mpsi,:)
     $     +hs%fs1(mpsi,:)*dpsi/2)*dpsi
c-----------------------------------------------------------------------
c     diagnose hs.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(UNIT=ascii_unit,FILE='hs.out')
         CALL open_bin(binary_unit,"hs.bin","UNKNOWN","REWIND",32_i4)
         hs%title=(/" psi  ","  h1  ","  h2  ","  h3  ","  h4  "/)
         CALL spline_write(hs,.TRUE.,.TRUE.,ascii_unit,binary_unit)
         CLOSE(UNIT=ascii_unit)
         CALL close_bin(binary_unit,"hs.bin")
      ENDIF
c-----------------------------------------------------------------------
c     compute global quantities.
c-----------------------------------------------------------------------
      gt%rmean=(gt%rs1+gt%rs2)/2
      gt%amean=(gt%rs2-gt%rs1)/2
      gt%aratio=gt%rmean/gt%amean
      gt%bt0=(gt%sq%fs(mpsi,1)+gt%sq%fs1(mpsi,1)*dpsi/2)/gt%rmean
      gt%betat=2*(hs%fsi(mpsi,1)/hs%fsi(mpsi,2))/gt%bt0**2
      gt%betan=twopi*100*gt%amean*gt%bt0*gt%betat/gt%crnt/twopi
      gt%betap1=2*(hs%fsi(mpsi,1)/hs%fsi(mpsi,2))/bp0**2
      gt%betap2=4*hs%fsi(mpsi,1)/((1.e6*mu0*gt%crnt)**2*gt%ro)
      gt%betap3=4*hs%fsi(mpsi,1)/((1.e6*mu0*gt%crnt)**2*gt%rmean)
      gt%li1=hs%fsi(mpsi,3)/hs%fsi(mpsi,2)/bp0**2*twopi
      gt%li2=2*hs%fsi(mpsi,3)/((1.e6*mu0*gt%crnt)**2*gt%ro)*twopi
      gt%li3=2*hs%fsi(mpsi,3)/((1.e6*mu0*gt%crnt)**2*gt%rmean)*twopi
      gt%volume=hs%fsi(mpsi,2)
c-----------------------------------------------------------------------
c     compute scale factors.
c-----------------------------------------------------------------------
      volume=gt%volume
      tauafac=mu0*mi*(gt%ro**2*gt%q0/gt%sq%fs(0,1))**2
      taurfac=volume/(2*pi**2*gt%ro)
      taua=sqrt(tauafac*ndens)
      taur=sfac*taua
      elecd=taurfac/taur
c-----------------------------------------------------------------------
c     Calculate lambda=J.B/B^2 and Delstar(psi) 
c	J.B = delstar(psi) F/R^2 - F' gss/R^2
c	B^2 = (F^2 + gss)/R^2
c     We use 2 Delstar(psi's):
c	Delstar(psi)=delsquare(psi) - 2/R Grad(R) dot Grad(psi)
c	Delstar(psi)= - r**2 * pprime - f * fprime = RHS of GS Eq. 
c     We calculate both for checking and because nimrod uses one or the
c      other depending on the j_t flag
c-----------------------------------------------------------------------
      ALLOCATE(lambdaprof(0:mpsi))
      ALLOCATE(lambda(0:mtheta,0:mpsi))
      ALLOCATE(delstr(0:mtheta,0:mpsi))
      ALLOCATE(gsrhs(0:mtheta,0:mpsi))

      psio=gt%psio						! total flux
      DO ipsi=0,mpsi
              f = gt%sq%fs(ipsi,1)				! R B_tor
              fprime = gt%sq%fs1(ipsi,1)/psio			! dF/dpsi
              p = gt%sq%fs(ipsi,1)				! R B_tor
              pprime = gt%sq%fs1(ipsi,2)/psio			! mu0 pprime
              mach = gt%sq%fs(ipsi,4)				! Mach
              msprime = gt%sq%fs1(ipsi,4)/psio			! Mach'
          DO itheta=0,mtheta
              r = gt%twod%fs(1,itheta,ipsi)
              jacfac = gt%twod%fs(3,itheta,ipsi)
              gss = gt%twod%fs(4,itheta,ipsi)
              gst = gt%twod%fs(5,itheta,ipsi)
              delstr(itheta,ipsi) = 
     &       gt%twod%fsy(4,itheta,ipsi)/psio
     &     + gt%twod%fsx(5,itheta,ipsi)
     &     + gss/jacfac *  gt%twod%fsy(3,itheta,ipsi)/psio
     &     + gst/jacfac *  gt%twod%fsx(3,itheta,ipsi)
     &     - 2._r8/r * (gss *gt%twod%fsy(1,itheta,ipsi)/psio
     &               + gst *gt%twod%fsx(1,itheta,ipsi) )
              gsrhs(itheta,ipsi) = - f * fprime
     &            - r**2 * exp(5./6.*(mach*(r/gt%ro-1.))**2)
     &            * (pprime + p*5./3.*msprime*mach*(r/gt%ro)**2*mu0)

              ldelstr= delstr(itheta,ipsi)
              IF (j_t == "use_gs") ldelstr= gsrhs(itheta,ipsi)

              jdotb = (ldelstr*f - fprime*gss)/r**2
              bsq = (f**2 + gss)/r**2
              lambda(itheta,ipsi) = jdotb/bsq
          ENDDO
      ENDDO
      CALL fluxav(gt,lambda,lambdaprof)			! <J.B/B^2>
c-----------------------------------------------------------------------
c     find highest node.
c-----------------------------------------------------------------------
      itheta=0
      z=gt%zo
      DO
         itheta=itheta+1
         rfac=SQRT(gt%ob%fs(itheta,1))
         eta=thetas(itheta)+gt%ob%fs(itheta,2)
         zold=z
         z=gt%zo+rfac*SIN(eta)
         IF(z < zold)EXIT
      ENDDO
      thetaext(1)=thetas(itheta-1)
c-----------------------------------------------------------------------
c     find lowest node.
c-----------------------------------------------------------------------
      DO
         itheta=itheta+1
         rfac=SQRT(gt%ob%fs(itheta,1))
         eta=thetas(itheta)+gt%ob%fs(itheta,2)
         zold=z
         z=gt%zo+rfac*SIN(eta)
         IF(z > zold)EXIT
      ENDDO
      thetaext(2)=thetas(itheta-1)
c-----------------------------------------------------------------------
c     find top and bottom.
c-----------------------------------------------------------------------
      DO iside=1,2
         theta=thetaext(iside)
         DO
            CALL spline_eval(gt%ob,theta,2_i4)
            r2=gt%ob%f(1)
            r21=gt%ob%f1(1)
            r22=gt%ob%f2(1)
            eta=gt%ob%f(2)
            eta1=gt%ob%f1(2)
            eta2=gt%ob%f2(2)
            rfac=SQRT(r2)
            rfac1=r21/(2*rfac)
            rfac2=(r22-r21*rfac1/rfac)/(2*rfac)
            phase=theta+eta
            phase1=1+eta1
            phase2=eta2
            cosfac=COS(phase)
            sinfac=SIN(phase)
            r=gt%ro+rfac*cosfac
            z=gt%zo+rfac*sinfac
            z1=rfac*phase1*cosfac+rfac1*sinfac
            z2=(2*rfac1*phase1+rfac*phase2)*cosfac
     $           +(rfac2-rfac*phase1**2)*sinfac
            dtheta=-z1/z2
            theta=theta+dtheta
            IF(ABS(dtheta) <= 1.e-12*ABS(theta))EXIT
         ENDDO
         rext(iside)=r
         zext(iside)=z
      ENDDO
c-----------------------------------------------------------------------
c     compute kappa and delta.
c-----------------------------------------------------------------------
      gt%kappa=(zext(1)-zext(2))/(gt%rs2-gt%rs1)
      gt%delta1=(gt%rmean-rext(1))/gt%amean
      gt%delta2=(gt%rmean-rext(2))/gt%amean
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_global
c-----------------------------------------------------------------------
c     subprogram 2. get_stability.
c     Calculate useful stability parameters
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE get_stability(gt)

      TYPE(global_type), INTENT(IN) :: gt
      TYPE(spline_type) :: es

      INTEGER(i4) :: itheta, ipsi

      REAL(r8) :: r, psio, fcl, fcu, bmax
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: bmod, gss, jacinv, h
      REAL(r8), DIMENSION(0:mpsi) :: dpdpsi, dqdpsi, qsf, f
      REAL(r8), DIMENSION(0:mpsi) :: dvdpsi, d2vdpsi2, r2bar
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: quantity
      REAL(r8), DIMENSION(0:mpsi) :: aveh1, aveh2, aveh3, aveh4, aveh5
c-----------------------------------------------------------------------
c     Allocate stability arrays
c-----------------------------------------------------------------------
      ALLOCATE(dideal(0:mpsi))
      ALLOCATE(dres(0:mpsi))
      ALLOCATE(hfactor(0:mpsi))
      ALLOCATE(dnc(0:mpsi))
      ALLOCATE(alphas(0:mpsi))
      ALLOCATE(fcirc(0:mpsi))
c-----------------------------------------------------------------------
c     Load stuff up into arrays to make code more readable
c-----------------------------------------------------------------------
      psio = gt%psio
      f(:) = gt%sq%fs(:,1)                            	! R B_tor
      dpdpsi(:) = gt%sq%fs1(:,2)/psio			! mu0 dp/dpsi
      qsf(:) = gt%sq%fs(:,3)
      dqdpsi(:) = gt%sq%fs1(:,3)/psio

      DO ipsi = 0,mpsi
       DO itheta = 0,mtheta
        r = gt%twod%fs(1,itheta,ipsi)
        gss(itheta,ipsi) = gt%twod%fs(4,itheta,ipsi)
        bmod(itheta,ipsi) = (f(ipsi)**2 +gss(itheta,ipsi))**0.5/r
        jacinv(itheta,ipsi) = 1._r8/gt%twod%fs(3,itheta,ipsi)
       ENDDO
      ENDDO

c-----------------------------------------------------------------------
c     Calculate circulating particle fraction
c     References:
c       Lin-Liu & Miller Phys. Plas. 2 (1995) 1667
c		-- Provides remarkably good estimate of fcirc
c-----------------------------------------------------------------------
      DO ipsi = 0,mpsi
        bmax = bmod(INT(mtheta/2),ipsi)
        DO itheta = 0,mtheta
              bmax= MAX(bmax,bmod(itheta,ipsi))
        ENDDO
        DO itheta = 0,mtheta
              h(itheta,ipsi)= bmod(itheta,ipsi)/bmax
              h(itheta,ipsi)= MIN(h(itheta,ipsi),1.0_r8)
        ENDDO
      ENDDO

      CALL fluxav(gt,h,aveh1)
      quantity(:,:) = h(:,:)**2
      CALL fluxav(gt,quantity,aveh2)
      quantity(:,:) = h(:,:)**(-2)* (1._r8 - (1._r8 - h(:,:))**0.5
     &             * (0.5*h(:,:) + 1._r8)     )
      CALL fluxav(gt,quantity,aveh3)

      DO ipsi = 0,mpsi
           fcl = aveh2(ipsi)*aveh3(ipsi)			! Lower estimate
           fcu = aveh2(ipsi)/aveh1(ipsi)**2 *
     &       ( 1._r8 -(1._r8 -aveh1(ipsi))**0.5
     &          *   (1._r8 + aveh1(ipsi)*0.5) ) 		! Upper estimate
           fcirc(ipsi) = 0.75*fcu+0.25*fcl
      ENDDO

c-----------------------------------------------------------------------
c     Calculate D_I, D_R, H
c     References:
c	Glasser, Greene, Johnson, Phys. Fluids 18 (1975) 875
c		-- Original paper with D_R, H
c	Greene, Comments on Modern Physics Part E, 17 (1997) 
c		-- Nice form of D_I, D_R, and H
c-----------------------------------------------------------------------

      CALL fluxav(gt,jacinv,aveh1)
      dvdpsi = 4._r8 *pi**2/aveh1
      CALL spline_alloc(es,mpsi,1_i4)
      es%xs(:) = gt%sq%xs(:)
      es%fs(:,1) = dvdpsi(:)
      CALL spline_fit(es,"extrap")
      d2vdpsi2(:) = es%fs1(:,1)/psio
      CALL spline_dealloc(es)

      quantity = bmod**2/gss
      CALL fluxav(gt,quantity,aveh1)			!aveh1 = <B^2/gss>
      quantity = 1._r8/gss
      CALL fluxav(gt,quantity,aveh2)			!aveh2 = <1/gss>
      quantity = 1._r8/(bmod**2*gss)
      CALL fluxav(gt,quantity,aveh3)			!aveh3 = <1/(B^2 gss)>
      quantity = 1._r8/bmod**2
      CALL fluxav(gt,quantity,aveh4)			!aveh4 = <1/B^2>
      quantity = bmod**2
      CALL fluxav(gt,quantity,aveh5)			!aveh5 = <B^2>



      dideal = -1._r8/4._r8 + dpdpsi*dvdpsi/(4*pi**2*dqdpsi)**2 * (
     &	        - aveh1 * d2vdpsi2 + 4*pi**2*dqdpsi * f * aveh2
     &          + dpdpsi*dvdpsi * (aveh1*aveh4 
     &          +          f**2 * (aveh1*aveh3 - aveh2**2) )
     &    )

      hfactor = dpdpsi*dvdpsi/(4*pi**2*dqdpsi)*f*(aveh2 - aveh1/aveh5 )

      dres = dideal + (1._r8/2._r8 - hfactor)**2

c-----------------------------------------------------------------------
c     Calculate D_nc, alpha_s
c     References:
c       Kruger, Hegna, Callen, Phys. Plas. 5 (1998) 455
c		-- Paper discussing D_nc vs. D_R
c       Chris Hegna Phys. Plasmas 6 (1999) 3980
c		-- Correct theory of D_nc and D_R for nonlinear islands
c-----------------------------------------------------------------------

      r2bar(:)=f(:) *dvdpsi(:)/(4 * pi**2 * qsf(:))

      dnc(:) = -1.58*(1.-fcirc(:))/(1.58 - 0.58*fcirc(:)) *
     &   dpdpsi(:)/dqdpsi(:)  * qsf(:) * aveh1(:) * r2bar(:)/aveh5(:)
c     dnc(:) = 
c    &  -dpdpsi(:)/dqdpsi(:)  * qsf(:) * aveh1(:) * r2bar(:)/aveh5(:)

      DO ipsi = 0,mpsi
         IF (dideal(ipsi) .gt. 0 ) THEN
             alphas(ipsi) = 1._r8
         ELSE
             alphas(ipsi) = 1._r8/2._r8 + SQRT(-dideal(ipsi))	!small index
         ENDIF
      ENDDO

c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_stability



c-----------------------------------------------------------------------
c     subprogram 3. fluxav.
c     Takes flux-surface average of 2D quantities
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fluxav(gt,quantity,average)

      TYPE(global_type), INTENT(IN) :: gt
      REAL(r8), DIMENSION(0:mtheta,0:mpsi), INTENT(IN) :: quantity
      REAL(r8), DIMENSION(0:mpsi), INTENT(OUT) :: average
      TYPE(spline_type) :: intgd

      REAL(r8) :: jac
      INTEGER(i4) :: ipsi,itheta
c-----------------------------------------------------------------------
c     integrate over each flux surface.
c-----------------------------------------------------------------------
      CALL spline_alloc(intgd,mtheta,2_i4)
      intgd%xs= gt%twod%xs

      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            jac = gt%twod%fs(3,itheta,ipsi) 
            intgd%fs(itheta,1)= jac
            intgd%fs(itheta,2)=jac*quantity(itheta,ipsi)
         ENDDO
         CALL spline_fit(intgd,"periodic")
         CALL spline_int(intgd)
         average(ipsi)=intgd%fsi(mtheta,2)/intgd%fsi(mtheta,1)
      ENDDO

      CALL spline_dealloc(intgd)
c-----------------------------------------------------------------------
c     A routine just isn't a routine without the following comment:
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fluxav




c-----------------------------------------------------------------------
c     subprogram 4. qsearch.
c     Finds rational surfaces equal to qwant between psimin and psimax.
c-----------------------------------------------------------------------
      SUBROUTINE q_search(psimin,psimax)
      IMPLICIT NONE
      REAL(r8), INTENT(in) :: psimin,psimax
      REAL(r8) :: tstart,tend
      TYPE :: chain
        REAL(r8) :: psirs
        TYPE(chain), POINTER :: next    ! ptr to next node in queue
      END TYPE chain
      TYPE(chain), POINTER :: start,item
c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep
      INTEGER(i4), PARAMETER :: neq=1,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq+1) :: y
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8) :: tol0 = 1.0e-06
      INTEGER(i4) :: maxstep=5000
c-----------------------------------------------------------------------

c
c                                               set up integrator parameters.
      itask=5
      iopt=1
      mf=10
      itol=1
      rtol=tol0
      atol=tol0*tol0
      iwork=0
      rwork=0
      rwork(11)=1
      iwork(6)=2500

      istate=1
      rwork(1)=psimax
      tstart=psimin
      tend=psimax
      NULLIFY(start,item)
      num_root=0
      range: DO
        y(1)=0
        istep=0
        internal : DO
          CALL dlsode(qres,neq,y,tstart,tend,itol,rtol,atol,itask,
     $                istate,iopt,rwork,lrw,iwork,liw,jac,mf)
          istep=istep+1
          IF(istate < 0) STOP 100
          IF(istep > maxstep)THEN
              WRITE(6,*)"maxstep exceeded in integration"
              STOP 200
          ENDIF
          IF(ABS(y(2)) > HUGE(0)*0.01 )THEN
             IF (.NOT. ASSOCIATED (start) ) THEN
               ALLOCATE(start)
               item => start
             ELSE
               ALLOCATE(item%next)
               item => item%next
             ENDIF
             num_root=num_root+1
             item%psirs=tstart
             NULLIFY (item%next)
             EXIT internal
          ENDIF
          IF(tend-tstart < tol0)EXIT range
        ENDDO internal
        tstart=tstart+tol0
        istate=1
        IF(tend-tstart < tol0)EXIT range
      ENDDO range
c
c						Count the number of surfaces.

      IF(ASSOCIATED(psi_root))DEALLOCATE(psi_root)
      ALLOCATE(psi_root(num_root))
      num_root=0
      item => start
      DO WHILE(ASSOCIATED(item))
        num_root=num_root+1
        psi_root(num_root)=item%psirs
        item => item%next
      ENDDO
c
c						Deallocate step
      item => start
      DO WHILE(ASSOCIATED(item))
        start => item
        item => item%next
        DEALLOCATE(start)
      ENDDO

      RETURN
      END SUBROUTINE q_search
c-----------------------------------------------------------------------
c     subprogram 5. qres.
c     Integration function for lsode in qsearch
c-----------------------------------------------------------------------
      SUBROUTINE qres(neq,eta,y,dy)
      USE spline
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: eta
      REAL(r8), INTENT(INOUT) :: y(neq+1)
      REAL(r8), INTENT(OUT) :: dy(neq+1)

c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      CALL spline_eval(qprof,eta,0_i4)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1)=1.0/(q_search_value-qprof%f(1))
      y(2)=dy(1)
      RETURN
      END SUBROUTINE qres
c-----------------------------------------------------------------------
      END MODULE analyze
