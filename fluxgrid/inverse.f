c-----------------------------------------------------------------------
c     subprogram inverse.
c     processes inverse equilibria.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. inverse.
c     1. process_eq.
c-----------------------------------------------------------------------
c     subprogram 0. inverse.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE inverse
      USE local
      USE physdat
      USE spline
      USE bicube
      USE global
      USE input
      USE grid
      IMPLICIT none

      TYPE :: inverse_input_type
      LOGICAL :: symflag,p1flag
      INTEGER(i4) :: ma,mtau
      REAL(r8) :: ro,zo,psio
      REAL(r8), DIMENSION(:,:), POINTER :: rg,zg
      TYPE(spline_type) :: sq_in
      TYPE(bicube_type), POINTER :: eigvec
      END TYPE inverse_input_type

      TYPE(inverse_input_type) :: ii
      TYPE(spline_type) :: p1s

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. process_eq.
c     reads and writes equilibrium data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE process_eq(gt)
      
      TYPE(global_type), INTENT(OUT) :: gt

      INTEGER(i4) :: itau,itheta,ipsi,jtau,ia
      REAL(r8) :: r,eta,f0,f0fac,theta,theta1,ffac,psifac,bp,rfac,eta0
      REAL(r8) :: v11,v12,v21,v22,g11,g22,g12,tau,jacfac,r2, jac
      REAL(r8), DIMENSION(:), ALLOCATABLE :: taus,thetas,eqarc
      TYPE(spline_type) :: dtheta,ff
      TYPE(bicube_type) :: r2g0,twod0
c-----------------------------------------------------------------------
c     Set type of equilibria
c-----------------------------------------------------------------------
      gt%eqtype='inverse'
c-----------------------------------------------------------------------
c     copy axis quantities.
c-----------------------------------------------------------------------
      gt%ro=ii%ro
      gt%zo=ii%zo
      gt%psio=ii%psio
c-----------------------------------------------------------------------
c     integrate pressure profile and fit to cubic splines.
c-----------------------------------------------------------------------
      IF(ii%p1flag)THEN
         CALL spline_alloc(p1s,ii%ma,1_i4)
         p1s%xs=ii%sq_in%xs
         p1s%fs(:,1)=ii%sq_in%fs(:,2)
         CALL spline_fit(p1s,"extrap")
         CALL spline_int(p1s)
         ii%sq_in%fs(:,2)=(p1s%fsi(:,1)-p1s%fsi(ii%ma,1))*ii%psio
      ENDIF
      CALL spline_fit(ii%sq_in,"extrap")
c-----------------------------------------------------------------------
c     diagnose 1D input.
c-----------------------------------------------------------------------
      IF(out_1d.OR.bin_1d)THEN
        IF(out_1d)OPEN(UNIT=ascii_unit,FILE='1d.out',STATUS='unknown')
        IF(bin_1d)CALL open_bin(binary_unit,"1d.bin","UNKNOWN",
     $        "REWIND",32_i4)
        ii%sq_in%title=(/' psi  ','  f   ','  p   ','  q   ',' mach '/)
        CALL spline_write(ii%sq_in,out_1d,bin_1d,ascii_unit,binary_unit)
        IF(out_1d)CLOSE(ascii_unit)
        IF(bin_1d)CALL close_bin(binary_unit,"1d.bin")
      ENDIF
c-----------------------------------------------------------------------
c     symmetrize coordinates.
c-----------------------------------------------------------------------
      IF(ii%symflag)THEN
         DO ia=0,ii%ma
            jtau=2*ii%mtau
            DO itau=0,ii%mtau-1
               ii%rg(jtau,ia)=ii%rg(itau,ia)
               ii%zg(jtau,ia)=-ii%zg(itau,ia)
               jtau=jtau-1
            ENDDO
         ENDDO
         ii%mtau=2*ii%mtau
      ENDIF
c-----------------------------------------------------------------------
c     prepare poloidal grids.
c-----------------------------------------------------------------------
      ALLOCATE(taus(0:ii%mtau))
      taus=twopi*(/(itau,itau=0,ii%mtau)/)/DFLOAT(ii%mtau)
      ALLOCATE(thetas(0:ii%mtau))
      ALLOCATE(eqarc(0:ii%mtau))
      CALL spline_alloc(ff,ii%mtau,4_i4)
      CALL spline_alloc(dtheta,ii%mtau,3_i4)
      dtheta%xs=taus
c-----------------------------------------------------------------------
c     convert to radial coordinates.
c-----------------------------------------------------------------------
      CALL bicube_alloc(r2g0,ii%mtau,ii%ma,2_i4)
      DO ia=0,ii%ma
         eta0=-twopi
         DO itau=0,ii%mtau-1
            r2=(ii%rg(itau,ia)-ii%ro)**2+(ii%zg(itau,ia)-ii%zo)**2
            eta=ATAN2(ii%zg(itau,ia)-ii%zo,ii%rg(itau,ia)-ii%ro)
            IF(eta < eta0)eta = eta + twopi
            eta0=eta
            r2g0%fs(1,itau,ia)=r2
            r2g0%fs(2,itau,ia)=eta-taus(itau)
         ENDDO
         r2g0%fs(1,ii%mtau,ia)=r2g0%fs(1,0,ia)
         r2g0%fs(2,ii%mtau,ia)=r2g0%fs(2,0,ia)
      ENDDO
c-----------------------------------------------------------------------
c     fit r2 and eta to bicubic splines as functions of tau and a.
c-----------------------------------------------------------------------
      r2g0%xs=taus
      r2g0%ys=ii%sq_in%xs
      CALL bicube_fit(r2g0,"periodic","extrap")
c-----------------------------------------------------------------------
c     Fit R,Z in new spline type to make interpolation easier
c-----------------------------------------------------------------------
      CALL bicube_alloc(twod0,ii%mtau,ii%ma,2_i4)
      twod0%xs=taus
      twod0%ys=ii%sq_in%xs
      twod0%fs(1,:,:)=ii%rg(:,:)
      twod0%fs(2,:,:)=ii%zg(:,:)
      CALL bicube_fit(twod0,"periodic","extrap")
c-----------------------------------------------------------------------
c     compute surface quantities on new surface.
c     Doing something cheesy to get the jacobian - use rho, eta in center
c	and use R,Z on the outside.  Empirically works.
c     Need to  store r2g & twod info in ff since we don't know theta
c-----------------------------------------------------------------------
      IF(grid_method == "original")mpsi=ii%sq_in%nodes
      CALL spline_alloc(gt%sq,mpsi,4_i4)
      CALL psigrid(gt,gt%sq,ii%sq_in%xs,ii%sq_in%fs(:,3))
      mpsi = SIZE(gt%sq%xs) - 1
      CALL bicube_alloc(gt%r2g,mtheta,mpsi,2_i4)
      CALL bicube_alloc(gt%twod,mtheta,mpsi,6_i4)
      gt%r2g%ys=gt%sq%xs
      gt%twod%ys=gt%r2g%ys
      DO ipsi=mpsi,0,-1
         psifac=gt%sq%xs(ipsi)
         CALL spline_eval(ii%sq_in,psifac,1_i4)
         gt%sq%fs(ipsi,1)=ii%sq_in%f(1)			! F=R Btoroidal
         gt%sq%fs(ipsi,2)=ii%sq_in%f(2)			! mu0 Pressure
         gt%sq%fs(ipsi,4)=ii%sq_in%f(4)			! Mach

         DO itau=0,ii%mtau
            tau=taus(itau)
            CALL bicube_eval(r2g0,tau,psifac,2_i4)
            CALL bicube_eval(twod0,tau,psifac,2_i4)
            thetas(itau)=tau+r2g0%f(2)		! defines:angle_method='geom'
            rfac=SQRT(r2g0%f(1))
c            r=ii%ro+rfac*COS(tau+r2g0%f(2))
            r=twod0%f(1)
            IF (SQRT(gt%sq%xs(ipsi)) < rjac) THEN
              jacfac=r2g0%fy(1)*(1+r2g0%fx(2))-r2g0%fx(1)*r2g0%fy(2)
            ELSE
              jacfac = 2.*(
     &         twod0%fy(1)*twod0%fx(2)-twod0%fx(1)*twod0%fy(2))
            ENDIF
            v11=(1+r2g0%fx(2))*2*rfac/jacfac
            v12=-r2g0%fx(1)/(rfac*jacfac)
            bp=SQRT(v11*v11+v12*v12)/r
            dtheta%fs(itau,1)=r*jacfac*bp**ipb/r**ipr
            dtheta%fs(itau,2)=jacfac/r
            dtheta%fs(itau,3)=r*jacfac*bp

            ff%fs(itau,1)=r2g0%f(1)
            ff%fs(itau,2)=r2g0%f(2)
            ff%fs(itau,3)=r
            ff%fs(itau,4)=twod0%f(2)
         ENDDO        
c-----------------------------------------------------------------------
c        set up and spline to new poloidal grid
c-----------------------------------------------------------------------
         CALL spline_fit(dtheta,"periodic")
         CALL spline_int(dtheta)
         IF(TRIM(angle_method) == 'jac') THEN
           thetas=twopi*dtheta%fsi(:,1)/dtheta%fsi(ii%mtau,1)
         ENDIF
         IF (ipsi == mpsi) THEN
            eqarc=twopi*dtheta%fsi(:,3)/dtheta%fsi(ii%mtau,3)
            CALL taugrid(thetas,eqarc,gt%r2g%xs,ii%mtau)
         ENDIF

         ff%xs=thetas
         ff%fs(:,2)=ff%fs(:,2)+taus-thetas
         CALL spline_fit(ff,"periodic")
         DO itheta=0,mtheta
            CALL spline_eval(ff,gt%r2g%xs(itheta),2_i4)
            gt%r2g%fs(1:2,itheta,ipsi)=ff%f(1:2)
            gt%twod%fs(1:2,itheta,ipsi)=ff%f(3:4)
         ENDDO
c                                                            Recompute q
         gt%sq%fs(ipsi,3)=ii%sq_in%f(1)*dtheta%fsi(ii%mtau,2)
     $        /(4*pi*gt%psio)
      ENDDO
c-----------------------------------------------------------------------
c     fit r2 and eta to bicubic splines as functions of theta and psi.
c-----------------------------------------------------------------------
      CALL bicube_fit(gt%r2g,"periodic","extrap")
      CALL bicube_dealloc(r2g0)
      CALL bicube_dealloc(twod0)
c-----------------------------------------------------------------------
c     Now that we're on new grid, can finish gt%twod
c-----------------------------------------------------------------------
      gt%twod%xs=gt%r2g%xs
      CALL bicube_alloc(twod0,mtheta,mpsi,2_i4)
      twod0%xs=gt%twod%xs
      twod0%ys=gt%twod%ys
      twod0%fs(1:2,:,:)=gt%twod%fs(1:2,:,:)
      CALL bicube_fit(twod0,"periodic","extrap")
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            r=twod0%fs(1,itheta,ipsi)
            rfac=SQRT(gt%r2g%fs(1,itheta,ipsi))
            IF (SQRT(gt%sq%xs(ipsi)) < rjac) THEN
              jacfac=gt%r2g%fsy(1,itheta,ipsi)*
     &                 (1.+gt%r2g%fsx(2,itheta,ipsi))-
     &            gt%r2g%fsx(1,itheta,ipsi)*gt%r2g%fsy(2,itheta,ipsi)
            ELSE
              jacfac = 2.*(
     &         twod0%fsy(1,itheta,ipsi)*twod0%fsx(2,itheta,ipsi)-
     &          twod0%fsx(1,itheta,ipsi)*twod0%fsy(2,itheta,ipsi))
            ENDIF
            jac=r*jacfac/2./gt%psio
            v11=(1+gt%r2g%fsx(2,itheta,ipsi))*2*rfac/jacfac
            v12=-gt%r2g%fsx(1,itheta,ipsi)/(rfac*jacfac)
            v21= gt%r2g%fsy(2,itheta,ipsi) *2*rfac/jacfac
            v22=-gt%r2g%fsy(1,itheta,ipsi)/(rfac*jacfac)
            g11= (v11*v11+v12*v12)*gt%psio**2   ! g^(psi,psi)
            g12=-(v11*v21+v12*v22)*gt%psio      ! g^(psi,theta)
            g22= (v21*v21+v22*v22)              ! g^(theta,theta)

            gt%twod%fs(3,itheta,ipsi)=jac
            gt%twod%fs(4,itheta,ipsi)=g11
            gt%twod%fs(5,itheta,ipsi)=g12
            gt%twod%fs(6,itheta,ipsi)=g22
         ENDDO
      ENDDO
      CALL bicube_fit(gt%twod,"periodic","extrap")
c-----------------------------------------------------------------------
c     revise q profile.
c-----------------------------------------------------------------------
      gt%sq%title=(/' psi  ','  f   ','  p   ','  q   ',' mach '/)
      CALL spline_fit(gt%sq,"extrap")
      gt%q0=gt%sq%fs(0,3)-gt%sq%fs1(0,3)*gt%sq%xs(0)
      IF(newq0 /= 0)THEN
         f0=gt%sq%fs(0,1)-gt%sq%fs1(0,1)*gt%sq%xs(0)
         f0fac=f0**2*((newq0/gt%q0)**2-1)
         gt%q0=newq0
         DO ipsi=0,mpsi
            ffac=SQRT(1+f0fac/gt%sq%fs(ipsi,1)**2)
            gt%sq%fs(ipsi,1)=gt%sq%fs(ipsi,1)*ffac
            gt%sq%fs(ipsi,3)=gt%sq%fs(ipsi,3)*ffac
         ENDDO
         CALL spline_fit(gt%sq,"extrap")
      ENDIF
      gt%qa=gt%sq%fs(mpsi,3)+gt%sq%fs1(mpsi,3)*(1-gt%sq%xs(mpsi))
      CALL qfind(gt,gt%sq%xs,gt%sq%fs(:,3))
c-----------------------------------------------------------------------
c     fit outer boundary to cubic splines.
c-----------------------------------------------------------------------
      CALL spline_alloc(gt%ob,mtheta,4_i4)
      gt%ob%xs=gt%r2g%xs
      gt%ob%fs(:,1:2)=TRANSPOSE(gt%r2g%fs(:,:,mpsi))
      gt%ob%fs(:,3:4)=TRANSPOSE(gt%twod%fs(1:2,:,mpsi))
      CALL spline_fit(gt%ob,"periodic")
c-----------------------------------------------------------------------
c     For inverse codes, the separatrix position will be taken as
c      the last surface.  We may want to be more sophisticated
c      later, but for now let gt%rzsep = gt%ob
c-----------------------------------------------------------------------
      CALL spline_alloc(gt%rzsep,mtheta,2_i4)
      gt%rzsep%xs=gt%ob%xs
      gt%rzsep%fs(:,1:2)=gt%ob%fs(:,3:4)
      CALL spline_fit(gt%rzsep,"periodic")
c-----------------------------------------------------------------------
c     find separatrix positions: gt%rs1 and gt%rs2
c-----------------------------------------------------------------------
      theta=pi
      DO 
         CALL spline_eval(gt%ob,theta,1_i4) 
         theta1=(pi-gt%ob%f(2)-theta)/(1+gt%ob%f1(2))
         theta=theta+theta1
         IF(ABS(theta1) <= 1.e-10_r8) EXIT
      ENDDO
      CALL spline_eval(gt%ob,theta,0_i4) 
      gt%rs1=gt%ro-SQRT(gt%ob%f(1))
      gt%rs2=gt%ro+SQRT(gt%ob%fs(0,1))
c-----------------------------------------------------------------------
c     interpolate eigenvector to new grid.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(ii%eigvec))THEN
         ii%eigvec%xs=r2g0%xs
         ii%eigvec%ys=r2g0%ys
         CALL bicube_fit(ii%eigvec,"periodic","extrap")
         ALLOCATE(gt%eigvec(0:mpsi,0:mtheta,SIZE(ii%eigvec%fs,3)))
         DO ipsi=0,mpsi
            psifac=gt%r2g%ys(ipsi)
            DO itheta=0,mtheta
               tau=gt%r2g%xs(itheta)+gt%r2g%fs(2,itheta,ipsi)
               CALL bicube_eval(ii%eigvec,tau,psifac,0_i4)
               gt%eigvec(ipsi,itheta,:)=ii%eigvec%f
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE process_eq
      END MODULE inverse
