c-----------------------------------------------------------------------
c     file direct.f.
c     processes direct equilibria.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. direct.
c     1. process_eq.
c     2. get_bfield.
c     3. position.
c     4. flder.
c     5. refine.
c     6. write_rimgrid.
c     7. calc_vacuum.
c     8. write_piegrid.
c     9. separatrix.
c     10. insep.
c-----------------------------------------------------------------------
c     subprogram 0. direct.
c     common data for direct equilibrium.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE direct
      USE local
      USE spline
      USE bicube
      USE global
      USE input
      USE physdat
      USE grid
      IMPLICIT NONE

      TYPE :: direct_input_type
        REAL(r8) :: rmin,rmax,zmin,zmax,ro,zo,psio
      TYPE(spline_type) :: sq_in
      TYPE(bicube_type) :: psig
      END TYPE direct_input_type

      TYPE :: direct_bfield_type
        REAL(r8) :: psi,psir,psiz,psirz,psirr,psizz
        REAL(r8) :: br,bz,brr,brz,bzr,bzz
        REAL(r8) :: f,f1,p,p1,mach,mach1
      END TYPE direct_bfield_type

      TYPE(direct_input_type) :: di

      INTEGER(i4) :: nbsvert
      REAL(r8)  dr0,dz0, f0fac
      REAL(r8),DIMENSION(:),ALLOCATABLE :: rbs,zbs   ! outer boundary

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. process_eq.
c     gets equilibrium data and massages it.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE process_eq(gt)

      TYPE(global_type), INTENT(INOUT) :: gt
c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep
      INTEGER(i4), PARAMETER :: neq=4,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq) :: y
      REAL(r8), DIMENSION(lrw) :: rwork=0
c-----------------------------------------------------------------------
c     other local variables.
c-----------------------------------------------------------------------
      INTEGER(i4) :: itau,itheta,ipsi
      INTEGER(i4), PARAMETER :: nstep=16384
      REAL(r8) :: btor, jr,jz,jt,modpp
      REAL(r8) :: dr,rfac,r,z,psifac,eta,err,psi0,f0,q0
      REAL(r8) :: v11,v12,v21,v22,g11,g22,g12,jacob,jacfac
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ffac
      REAL(r8), DIMENSION(0:nstep,0:5) :: temp
      REAL(r8), PARAMETER :: eps=1e-12
      TYPE(direct_bfield_type) :: bf
      TYPE(spline_type) :: ff
      TYPE(bicube_type) :: rz
c-----------------------------------------------------------------------
c     Set type of equilibria
c-----------------------------------------------------------------------
      gt%eqtype='direct'
c-----------------------------------------------------------------------
c     Set up psi grid
c-----------------------------------------------------------------------
      CALL spline_alloc(gt%sq,mpsi,4_i4)
      gt%sq%name="gt%sq"
      CALL psigrid(gt,gt%sq,di%sq_in%xs,di%sq_in%fs(:,3))
      mpsi = SIZE(gt%sq%xs)-1
c-----------------------------------------------------------------------
c     Allocate global arrays
c-----------------------------------------------------------------------
      CALL bicube_alloc(gt%twod,mtheta,mpsi,6_i4)
      CALL bicube_alloc(rz,mtheta,mpsi,2_i4)
c-----------------------------------------------------------------------
c     find flux surface.
c-----------------------------------------------------------------------
      DO ipsi=mpsi,0,-1
         psifac=gt%sq%xs(ipsi)
         psi0=gt%psio*(1-psifac)
         r=gt%ro+SQRT(psifac)*(gt%rs2-gt%ro)
         z=gt%zo
         DO
            CALL get_bfield(r,z,bf,1_i4)
            dr=(psi0-bf%psi)/bf%psir
            r=r+dr
            IF(ABS(dr) <= eps*r)EXIT
         ENDDO
         psi0=bf%psi
c-----------------------------------------------------------------------
c     initialize variables.
c-----------------------------------------------------------------------
         istep=0
         eta=0.
         y(1)=0.
         y(2)=SQRT((r-gt%ro)**2+(z-gt%zo)**2)
         y(3)=0.
         y(4)=0.
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
         istate=1; itask=5; iopt=1; mf=10; itol=1; 
         rtol=tol0; atol=tol0*y(2); rwork(1)=twopi; rwork(11)=0
c-----------------------------------------------------------------------
c     advance differential equations  and store results for each step.
c			y(1) --> transformed angle
c			y(2) --> effective radial coordinate
c			y(3) --> geometric angle
c			y(4) --> equal-arc angle
c-----------------------------------------------------------------------
         DO
            rfac=y(2)
            CALL refine(rfac,eta,psi0)
            r=gt%ro+rfac*COS(eta)
            z=gt%zo+rfac*SIN(eta)
            CALL get_bfield(r,z,bf,2_i4)
            temp(istep,0)=y(1)			!   transformed angle
            temp(istep,1)=rfac**2		!   effective radius squared
            temp(istep,2)=eta			!   geometric angle.
            temp(istep,3)=y(4)			!   eq. arc angle
            temp(istep,4)=r			!   r
            temp(istep,5)=z			!   z
            err=(bf%psi-psi0)/bf%psi
            IF(eta >= twopi .OR. istep >= nstep  .OR.  istate < 0
     $            .OR. ABS(err) >= 1)EXIT
            istep=istep+1
            CALL dlsode(flder,neq,y,eta,twopi,itol,rtol,atol,itask,
     $                  istate,iopt,rwork,lrw,iwork,liw,jac,mf)
         ENDDO
c-----------------------------------------------------------------------
c     check for exceeding array bounds.
c-----------------------------------------------------------------------
         IF(eta < twopi)THEN
            WRITE(nim_wr,'(a,i4,a,1p,e10.3,a,i3)')
     $           ' In SUBROUTINE newsurf, istep = nstep =',nstep,
     $           ' at eta = ',eta,', ipsi =',ipsi
            STOP
         ENDIF
c-----------------------------------------------------------------------
c        Map temp(0:istep) -> r2g(0:mtheta,:,:), twod(0:mtheta,:,:)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c        Ensure periodicity
c-----------------------------------------------------------------------
         temp(0:istep,0)=temp(0:istep,0)*twopi/y(1)
         temp(0:istep,3)=temp(0:istep,3)*twopi/y(4)
c-----------------------------------------------------------------------
c        On first step, set xs (theta grid) and ys
c-----------------------------------------------------------------------
         IF(ipsi == mpsi)THEN
            IF(mtheta == 0)mtheta=istep
            CALL bicube_alloc(gt%r2g,mtheta,mpsi,2_i4)
            SELECT CASE(TRIM(angle_method))
            CASE('jac')
               CALL taugrid(temp(:,0),temp(:,3),gt%r2g%xs,istep)
            CASE('geom')
               CALL taugrid(temp(:,2),temp(:,3),gt%r2g%xs,istep)
            CASE default
              CALL nim_stop('Invalid angle_method')
            END SELECT
            gt%r2g%ys=gt%sq%xs
            gt%twod%ys = gt%sq%xs			! radial coord.
            gt%twod%xs = gt%r2g%xs			! theta
            rz%ys = gt%sq%xs				! radial coord.
            rz%xs = gt%r2g%xs				! theta
         ENDIF
c-----------------------------------------------------------------------
c        fit to cubic splines vs. theta and interpolate to tau grid.
c-----------------------------------------------------------------------
         CALL spline_alloc(ff,istep,4_i4)
         SELECT CASE(TRIM(angle_method))
         CASE('jac')
           ff%xs(0:istep)=temp(0:istep,0)               ! Calculated angle
         CASE('geom')
           ff%xs(0:istep)=temp(0:istep,2)               ! Geometric angle
         END SELECT
         ff%fs(0:istep,1)=temp(0:istep,1)
         ff%fs(0:istep,2)=temp(0:istep,2)-ff%xs(0:istep)
         ff%fs(0:istep,3)=temp(0:istep,4)
         ff%fs(0:istep,4)=temp(0:istep,5)
         CALL spline_fit(ff,"periodic")
         DO itau=0,mtheta
            CALL spline_eval(ff,gt%r2g%xs(itau),0_i4)
            gt%r2g%fs(1,itau,ipsi)=ff%f(1)
            gt%r2g%fs(2,itau,ipsi)=ff%f(2)
            gt%twod%fs(1,itau,ipsi)=ff%f(3)
            gt%twod%fs(2,itau,ipsi)=ff%f(4)
            rz%fs(1,itau,ipsi)=ff%f(3)
            rz%fs(2,itau,ipsi)=ff%f(4)
         ENDDO
         CALL spline_dealloc(ff)
c-----------------------------------------------------------------------
c     evaluate surface quantities.
c-----------------------------------------------------------------------
         gt%sq%fs(ipsi,1)=bf%f
         gt%sq%fs(ipsi,2)=bf%p
         gt%sq%fs(ipsi,3)=y(3)*bf%f/twopi
         gt%sq%fs(ipsi,4)=bf%mach
      ENDDO
c-----------------------------------------------------------------------
c     fit to splines.
c-----------------------------------------------------------------------
      CALL bicube_fit(gt%r2g,"periodic","extrap")
      CALL bicube_fit(rz,"periodic","extrap")
      CALL spline_fit(gt%sq,"extrap")
      gt%sq%title= (/' psi  ','  f   ','  p   ','  q   ',' mach '/)
c-----------------------------------------------------------------------
c     revise f and q profiles.
c-----------------------------------------------------------------------
      IF(newq0 /= 0)THEN
         gt%q0 = newq0
         q0=gt%sq%fs(0,3)-gt%sq%fs1(0,3)*gt%sq%xs(0)
         f0=gt%sq%fs(0,1)-gt%sq%fs1(0,1)*gt%sq%xs(0)
         f0fac=f0**2*((newq0/q0)**2-1)
         ALLOCATE(ffac(0:mpsi))
         ffac=SQRT(1+f0fac/gt%sq%fs(:,1)**2)
         gt%sq%fs(:,1)=gt%sq%fs(:,1)*ffac			! F=R Btor	
         gt%sq%fs(:,3)=gt%sq%fs(:,3)*ffac			! q
         DEALLOCATE(ffac)
         CALL spline_fit(gt%sq,"extrap")
      ENDIF
      CALL qfind(gt,gt%sq%xs,gt%sq%fs(:,3))
c-----------------------------------------------------------------------
c     compute metric quantities.
c     Doing something cheesy to get the jacobian - use rho, eta in center
c	and use R,Z on the outside.  Empirically works.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            rfac=SQRT(gt%r2g%fs(1,itheta,ipsi))
            eta=gt%r2g%fs(2,itheta,ipsi)+gt%r2g%xs(itheta)
c            r=gt%ro+rfac*COS(eta)
            r=gt%twod%fs(1,itheta,ipsi)
            IF (SQRT(gt%sq%xs(ipsi)) < rjac) THEN
              jacfac=gt%r2g%fsy(1,itheta,ipsi)
     $           *(1.+gt%r2g%fsx(2,itheta,ipsi))
     $           -gt%r2g%fsx(1,itheta,ipsi)
     $           *gt%r2g%fsy(2,itheta,ipsi)
            ELSE
              jacfac = 2.*(
     $         rz%fsy(1,itheta,ipsi)*rz%fsx(2,itheta,ipsi)-
     $         rz%fsx(1,itheta,ipsi)*rz%fsy(2,itheta,ipsi))
            ENDIF
            v11=(1.+gt%r2g%fsx(2,itheta,ipsi))*2.*rfac/jacfac
            v21=   gt%r2g%fsy(2,itheta,ipsi) *2.*rfac/jacfac
            v12=-gt%r2g%fsx(1,itheta,ipsi)/(rfac*jacfac)
            v22=-gt%r2g%fsy(1,itheta,ipsi)/(rfac*jacfac)
            g11=  v11*v11+v12*v12
            g12=-(v11*v21+v12*v22)
            g22=  v21*v21+v22*v22
            jacob=r*jacfac/2.
            gt%twod%fs(3,itheta,ipsi)= jacob/gt%psio      ! Jac. w/ psi
            gt%twod%fs(4,itheta,ipsi)= g11*gt%psio**2   ! g^(psi,psi)
            gt%twod%fs(5,itheta,ipsi)= g12*gt%psio      ! g^(psi,theta)
            gt%twod%fs(6,itheta,ipsi)= g22	        ! g^(theta,theta)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fit to splines.
c-----------------------------------------------------------------------
      CALL bicube_fit(gt%twod,"periodic","extrap")
c-----------------------------------------------------------------------
c     If requested, calculate the bfield directly from the di% data.
c     NOTE: I actually haven't debugged with flow coming through EQDSK
c-----------------------------------------------------------------------
      IF (calc_direct) THEN
       IF(newq0 /= 0) CALL nim_stop("calc_direct AND newq0 not allowed")
        CALL bicube_alloc(gt%dir,mtheta,mpsi,6_i4)
        gt%dir%ys = gt%twod%ys			! radial coord.
        gt%dir%xs = gt%twod%xs			! theta
        DO ipsi=0,mpsi
          DO itheta=0,mtheta
            r= gt%twod%fs(1,itheta,ipsi)
            z= gt%twod%fs(2,itheta,ipsi)
            CALL get_bfield(r,z,bf,1_i4)
            btor=bf%f/r
            jr = -bf%f1*bf%br/mu0/gt%psio
            jz = -bf%f1*bf%bz/mu0/gt%psio
            modpp=exp(5./6.*(bf%mach*(r/gt%ro-1.))**2)*(bf%p1 
     &            + bf%p*5./3.*bf%mach1*bf%mach*(r/gt%ro)**2*mu0)
            jt = -(bf%f*bf%f1 + r**2 * modpp)/gt%psio/(mu0*r**2)
            gt%dir%fs(1,itheta,ipsi)=bf%br
            gt%dir%fs(2,itheta,ipsi)=bf%bz
            gt%dir%fs(3,itheta,ipsi)=btor
            gt%dir%fs(4,itheta,ipsi)=jr
            gt%dir%fs(5,itheta,ipsi)=jz
            gt%dir%fs(6,itheta,ipsi)=jt
          ENDDO
        ENDDO
        CALL bicube_fit(gt%dir,"periodic","extrap")
      ENDIF
c-----------------------------------------------------------------------
c     fit outer boundary to cubic splines and diagnose 2D quantities.
c     Note: ob = outer boundary of rblock region, and not separatrix
c-----------------------------------------------------------------------
      CALL spline_alloc(gt%ob,mtheta,4_i4)
      gt%ob%xs=gt%r2g%xs
      gt%ob%fs(:,1:2)=TRANSPOSE(gt%r2g%fs(:,:,mpsi))
      gt%ob%fs(:,3:4)=TRANSPOSE(gt%twod%fs(1:2,:,mpsi))
      CALL spline_fit(gt%ob,"periodic")
c-----------------------------------------------------------------------
c     Find separatrix (gt%rzsep spline type)
c-----------------------------------------------------------------------
      CALL separatrix(gt)
c-----------------------------------------------------------------------
c     Generate rblocks for vacuum region if asked for (gt%vac spl type)
c-----------------------------------------------------------------------
      IF (extr > 1. .AND. mvac > 0) CALL calc_vacuum(gt)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE process_eq
c-----------------------------------------------------------------------
c     subprogram 2. get_bfield.
c     evaluates bicubic splines for field components and derivatives.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE get_bfield(r,z,bf,mode)

      INTEGER(i4), INTENT(IN) :: mode
      REAL(r8), INTENT(IN) :: r,z
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     compute spline interpolations.
c-----------------------------------------------------------------------
      CALL bicube_eval(di%psig,r,z,mode)
      bf%psi=di%psig%f(1)
      CALL spline_eval(di%sq_in,1-bf%psi/di%psio,1_i4)
      bf%f=di%sq_in%f(1)
      bf%f1=di%sq_in%f1(1)
      bf%p=di%sq_in%f(2)
      bf%p1=di%sq_in%f1(2)
      bf%mach=di%sq_in%f(4)
      bf%mach1=di%sq_in%f1(2)
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     evaluate magnetic fields.
c-----------------------------------------------------------------------
      bf%psir=di%psig%fx(1)
      bf%psiz=di%psig%fy(1)
      bf%br=bf%psiz/r
      bf%bz=-bf%psir/r
      IF(mode == 1)RETURN
c-----------------------------------------------------------------------
c     evaluate derivatives of magnetic fields.
c-----------------------------------------------------------------------
      bf%psirr=di%psig%fxx(1)
      bf%psirz=di%psig%fxy(1)
      bf%psizz=di%psig%fyy(1)
      bf%brr=(bf%psirz-bf%br)/r
      bf%brz=bf%psizz/r
      bf%bzr=-(bf%psirr+bf%bz)/r
      bf%bzz=-bf%psirz/r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_bfield
c-----------------------------------------------------------------------
c     subprogram 3. position.
c     finds radial positions of o-point and separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE position(gt)
      
      TYPE(global_type), INTENT(OUT) :: gt

      REAL(r8), PARAMETER :: eps=1e-12
      REAL(r8) :: ajac(2,2),det,dr,dz,fac,r,z
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     scan to find zero crossing of bz on midplane.
c-----------------------------------------------------------------------
      IF(di%ro == 0)THEN
         r=(di%rmax+di%rmin)/2
         z=(di%zmax+di%zmin)/2
         dr=(di%rmax-di%rmin)/20
         DO
            CALL get_bfield(r,z,bf,1_i4)
            IF(bf%bz >= 0)EXIT
            r=r+dr
         ENDDO
      ELSE
         r=di%ro
         z=di%zo
      ENDIF
c-----------------------------------------------------------------------
c     use newton iteration to find o-point.
c-----------------------------------------------------------------------
      DO
         CALL get_bfield(r,z,bf,2_i4)
         ajac(1,1)=bf%brr
         ajac(1,2)=bf%brz
         ajac(2,1)=bf%bzr
         ajac(2,2)=bf%bzz
         det=ajac(1,1)*ajac(2,2)-ajac(1,2)*ajac(2,1)
         dr=(ajac(1,2)*bf%bz-ajac(2,2)*bf%br)/det
         dz=(ajac(2,1)*bf%br-ajac(1,1)*bf%bz)/det
         r=r+dr
         z=z+dz
         IF(ABS(dr) <= eps*r .AND. ABS(dz) <= eps*r)EXIT
      ENDDO
      di%ro=r
      di%zo=z
      gt%ro=r
      gt%zo=z
c-----------------------------------------------------------------------
c     renormalize psi.
c-----------------------------------------------------------------------
      fac=di%psio/bf%psi
      di%psig%fs=di%psig%fs*fac
      di%psig%fsx=di%psig%fsx*fac
      di%psig%fsy=di%psig%fsy*fac
      di%psig%fsxy=di%psig%fsxy*fac
c-----------------------------------------------------------------------
c     use newton iteration to find inboard separatrix position.
c-----------------------------------------------------------------------
      r=(3*di%rmin+di%ro)/4
      z=di%zo
      DO
         CALL get_bfield(r,z,bf,1_i4)
         dr=-bf%psi/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r)EXIT
      ENDDO
      gt%rs1=r
c-----------------------------------------------------------------------
c     use newton iteration to find outboard separatrix position.
c-----------------------------------------------------------------------
      r=(di%ro+3*di%rmax)/4
      DO
         CALL get_bfield(r,z,bf,1_i4)
         dr=-bf%psi/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r)EXIT
      ENDDO
      gt%rs2=r
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE position
c-----------------------------------------------------------------------
c     subprogram 4. flder.
c     contains differential equations for field line averages.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE flder(neq,eta,y,dy)
      
      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: eta
      REAL(r8), INTENT(IN) :: y(neq)
      REAL(r8), INTENT(OUT) :: dy(neq)

      REAL(r8) :: cosfac,sinfac,bp,r,z
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      cosfac=COS(eta)
      sinfac=SIN(eta)
      r=di%ro+y(2)*cosfac
      z=di%zo+y(2)*sinfac
      CALL get_bfield(r,z,bf,1_i4)
      bp=SQRT(bf%br**2+bf%bz**2)
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      dy(1)=y(2)/(bf%bz*cosfac-bf%br*sinfac)
      dy(2)=dy(1)*(bf%br*cosfac+bf%bz*sinfac)
      dy(3)=dy(1)/(r*r)
      dy(4)=dy(1)*bp
      dy(1)=dy(1)*bp**ipb/r**ipr
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE flder
c-----------------------------------------------------------------------
c     subprogram 5. refine.
c     moves a point orthogonally to a specified flux surface.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE refine(rfac,eta,psi0)
      
      REAL(r8) :: rfac,eta,psi0
      
      REAL(r8) :: dpsi,cosfac,sinfac,drfac,r,z
      REAL(r8), PARAMETER :: eps=1.e-12
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     initialize iteration.
c-----------------------------------------------------------------------
      cosfac=COS(eta)
      sinfac=SIN(eta)
      r=di%ro+rfac*cosfac
      z=di%zo+rfac*sinfac
      CALL get_bfield(r,z,bf,1_i4)
      dpsi=bf%psi-psi0
c-----------------------------------------------------------------------
c     refine rfac by newton iteration.
c-----------------------------------------------------------------------
      DO
         drfac=-dpsi/(bf%psir*cosfac+bf%psiz*sinfac)
         rfac=rfac+drfac
         r=di%ro+rfac*cosfac
         z=di%zo+rfac*sinfac
         CALL get_bfield(r,z,bf,1_i4)
         dpsi=bf%psi-psi0
         IF(ABS(dpsi) <= eps*psi0 .OR. ABS(drfac) <= eps*rfac)EXIT
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE refine
c-----------------------------------------------------------------------
c     subprogram 6. write_rimgrid.
c     reads positions, writes equilibrium quantities.
c     This is used for triangular regions.  Needs work because it doesn't
c      work for vertex regions inside the separatrix and doesn't write
c      out the current density projections like rblock region.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_rimgrid(gt)

      TYPE(global_type), INTENT(IN) :: gt

      INTEGER(i4) :: iv,idum,mvert,k,i
      INTEGER :: ios
      REAL(r8) :: rval,zval,btor,fedge,ffac
      TYPE(direct_bfield_type) :: bf
      CHARACTER(64) :: nodefile,datfile
      CHARACTER(8) :: kstring
c-----------------------------------------------------------------------
c     Evaluate magnetic fields at each vertex location.
c-----------------------------------------------------------------------
      k=1
      fedge=di%sq_in%fs(di%psig%mx,1)
      DO
         WRITE(kstring,fmt='(i4)')k
         DO
            i= INDEX(TRIM(kstring)," ")
            IF(i == 0)EXIT
            kstring(i:)=kstring(i+1:)
         ENDDO
         nodefile="rim"//TRIM(kstring)//".1.node"
         datfile="rim"//TRIM(kstring)//".dat"
         OPEN(UNIT=eq1_unit,FILE=TRIM(nodefile), STATUS='old',
     $        iostat=ios)
         IF(ios.NE.0)EXIT
         WRITE(6,*)"Processing file ",TRIM(nodefile)
         WRITE(6,*)"Writing file ",TRIM(datfile)
         READ(eq1_unit,*)mvert,idum,idum,idum
         OPEN(UNIT=eq2_unit,FILE=TRIM(datfile), STATUS='replace')
         DO iv=1,mvert
            READ(eq1_unit,*)idum,rval,zval,idum
            CALL get_bfield(rval,zval,bf,1_i4)
            IF (insep(gt,rval,zval)) THEN 
                  ! Everything is OK here
c-TMP new q0?
              btor=bf%f/rval
            ELSE
               bf%p=0._r8
               IF(rval == 0.) THEN
                  btor=0
               ELSE
                  ffac=1.
                  IF(newq0 /= 0.) ffac=SQRT(1.+f0fac/fedge**2)
                  btor=fedge*ffac/rval
               ENDIF
            ENDIF
            write(eq2_unit,*)bf%br,bf%bz,btor,bf%p
         ENDDO
        CLOSE(UNIT=eq1_unit)
        CLOSE(UNIT=eq2_unit)
        k=k+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_rimgrid
c-----------------------------------------------------------------------
c     subprogram 7. calc_vacuum.
c     Generates the gt%vac array which contains the information needed
c	for the fluxgrid.dat file in the region between the outer 
c	boundary as defined by psihigh and the boundary defined by extr
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE calc_vacuum(gt)

      TYPE(global_type), INTENT(INOUT) :: gt

      INTEGER :: ivac, itheta
      REAL(r8) :: rcntr,zcntr, rob, zob, rnew,znew, rnext, znext
      REAL(r8) :: slope, dsep, ds, crossprod, small
      REAL(r8) :: conc, btor, jr,jz,jt, modpp, ffac,fedge
      REAL(r8) :: agrid,bgrid
      REAL(r8), DIMENSION(0:mvac) :: grid
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     Allocate vac array.  
c     Calculation ssumes no currents outside of separatrix
c-----------------------------------------------------------------------
      grid=(/(REAL(ivac,r8),ivac=0,mvac)/)/REAL(mvac,r8)
c      IF(vac_pack) CALL vacpack(grid,mvac)
      CALL bicube_alloc(gt%vac,mtheta,mvac,11_i4)
      gt%vac%xs=gt%ob%xs
      gt%vac%ys=grid
c-----------------------------------------------------------------------
c     agrid and bgrid provide linear variation in element sizes,
c     set by the firstv parameter (fraction of distance).
c-----------------------------------------------------------------------
      IF (firstv/=0.AND.firstv>=2._r8/mvac.OR.firstv<0) THEN
        CALL nim_stop("Inappropriate value for firstv.")
      ELSE IF (firstv==0) THEN
        agrid=0._r8
        bgrid=1._r8/mvac
      ELSE
        agrid=(1._r8/mvac-firstv)/(0.5_r8*(mvac-1_i4))
        bgrid=firstv-agrid
      ENDIF
c-----------------------------------------------------------------------
c     Define grid - Note that ivac=0 is going to give the same points
c       as the outer boundary.
c-----------------------------------------------------------------------

c     							CONFORMAL WALL
      IF (extr_type == 'conformal') THEN

        dsep=0.5*(1.-extr)*(gt%rs2-gt%rs1)
        small =1./(0.1*HUGE(0))

        DO ivac=0,mvac
c         ds=dsep*grid(ivac)
c         ds=dsep*(REAL(ivac,r8)/REAL(mvac,r8))**2

          IF (ivac==0) THEN
            ds=0._r8
          ELSE
            ds=(agrid*ivac+bgrid)*dsep+ds
          ENDIF

          DO itheta=0,mtheta
             rob= gt%ob%fs(itheta,3)
             zob= gt%ob%fs(itheta,4)
c                                                       R_ob,Z_ob ordered
             IF (itheta == mtheta) THEN
               rnext=gt%ob%fs(1,3) - rob
               znext=gt%ob%fs(1,4) - zob
             ELSE
               rnext=gt%ob%fs(itheta+1,3) - rob
               znext=gt%ob%fs(itheta+1,4) - zob
             ENDIF
c                                                       Calc slope-defines perp
             IF (ABS(gt%ob%fs1(itheta,3)) < small) THEN
                slope=0.1*HUGE(0)*gt%ob%fs1(itheta,3)/
     &                         ABS(gt%ob%fs1(itheta,3))
             ELSE
                slope=gt%ob%fs1(itheta,4)/gt%ob%fs1(itheta,3)  !dZ/dth/dR/dth
             ENDIF
c                                                        SQRT can be +or-
             znew=ds/SQRT(1.+slope**2)
             rnew=-slope*znew
c                                                        curl determines sign
             crossprod=rnew*znext-znew*rnext
             IF (crossprod < 0.) THEN
                 rnew=-rnew+rob; znew=-znew+zob
             ELSE
                 rnew=rnew+rob; znew=znew+zob
             ENDIF
c                                                        Must be in domain
             IF(rnew < di%rmin) rnew=di%rmin
             IF(rnew > di%rmax) rnew=di%rmax
             IF(znew < di%zmin) znew=di%zmin
             IF(znew > di%zmax) znew=di%zmax

             gt%vac%fs(1,itheta,ivac)=rnew
             gt%vac%fs(2,itheta,ivac)=znew
          ENDDO
        ENDDO


c     							SELF-SIMILAR WALL
      ELSEIF (extr_type == 'self_similar') THEN

        dsep=extr - 1.
        IF (extr_center == 'geom') THEN
           rcntr=0.5*(gt%rs2+gt%rs1)
           zcntr=0.
        ELSEIF (extr_center == 'mag_axis') THEN
           rcntr=gt%ro
           zcntr=gt%zo
        ELSE
           CALL nim_stop("extr_center choice not recognized")
        ENDIF

        DO ivac=0,mvac
c         ds=1. + dsep*REAL(ivac)/REAL(mvac)
c         ds=1. + dsep*(REAL(ivac,r8)/REAL(mvac,r8))**2

          IF (ivac==0) THEN
            ds=1._r8
          ELSE
            ds=(agrid*ivac+bgrid)*dsep+ds
          ENDIF

          DO itheta=0,mtheta
             rob= gt%ob%fs(itheta,3)
             zob= gt%ob%fs(itheta,4)
             rnew= rcntr + ds*(rob-rcntr)
             znew= zcntr + ds*(zob-zcntr)
c                                                        Must be in domain
             IF(rnew < di%rmin) rnew=di%rmin
             IF(rnew > di%rmax) rnew=di%rmax
             IF(znew < di%zmin) znew=di%zmin
             IF(znew > di%zmax) znew=di%zmax

             gt%vac%fs(1,itheta,ivac)=rnew
             gt%vac%fs(2,itheta,ivac)=znew
          ENDDO
        ENDDO


      ELSE 
            CALL nim_stop("extr_type choice not recognized")
      ENDIF
c-----------------------------------------------------------------------
c     Assign the equilibrium values to the grid
c-----------------------------------------------------------------------
        fedge=di%sq_in%fs(di%psig%mx,1)
        DO ivac=0,mvac
          DO itheta=0,mtheta
            rnew= gt%vac%fs(1,itheta,ivac)
            znew= gt%vac%fs(2,itheta,ivac)
            CALL get_bfield(rnew,znew,bf,1_i4)
            IF(newq0 /= 0) THEN
                ffac=SQRT(1.+f0fac/bf%f**2)
                bf%f =bf%f * ffac
                bf%f1=bf%f1/ ffac
            ENDIF
            IF (insep(gt,rnew,znew)) THEN 
               conc=0.
               btor=bf%f/rnew
               jr = -bf%f1*bf%br/mu0/gt%psio
               jz = -bf%f1*bf%bz/mu0/gt%psio
               modpp=exp(5./6.*(bf%mach*(rnew/gt%ro-1.))**2)*(bf%p1 
     &          +bf%p*5./3.*bf%mach1*bf%mach*(rnew/gt%ro)**2*mu0)
               jt = -(bf%f*bf%f1 + rnew**2*modpp)/gt%psio/(mu0*rnew**2)
            ELSE
               bf%p=0.
               jr = 0.
               jz = 0.
               jt = 0.
               conc=1.
                              ffac=1.
               IF(newq0 /= 0) ffac=SQRT(1.+f0fac/fedge**2)
               IF(rnew == 0) THEN
                  btor=0.
               ELSE
                  btor=fedge*ffac/rnew
               ENDIF
            ENDIF
            gt%vac%fs(3,itheta,ivac)=bf%br
            gt%vac%fs(4,itheta,ivac)=bf%bz
            gt%vac%fs(5,itheta,ivac)=btor
            gt%vac%fs(6,itheta,ivac)=jr
            gt%vac%fs(7,itheta,ivac)=jz
            gt%vac%fs(8,itheta,ivac)=jt
            gt%vac%fs(9,itheta,ivac)=bf%p
            gt%vac%fs(10,itheta,ivac)=conc
            gt%vac%fs(11,itheta,ivac)=bf%psi
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE calc_vacuum
c-----------------------------------------------------------------------
c     subprogram 8. write_piegrid.
c     reads positions, writes equilibrium quantities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_piegrid

      INTEGER(i4) :: iv,idum,mvert
      INTEGER :: ios
      REAL(r8) :: rval,zval
      TYPE(direct_bfield_type) :: bf
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      OPEN(UNIT=eq1_unit,FILE='pie.1.node', STATUS='old',iostat=ios)
      IF(ios.NE.0)RETURN
      READ(eq1_unit,*)mvert,idum,idum,idum
      OPEN(UNIT=eq2_unit,FILE='pie.dat', STATUS='replace')
      WRITE(6,*)"Processing file pie.1.node"
      WRITE(6,*)"Writing file pie.dat"
      DO iv=1,mvert
        READ(eq1_unit,*)idum,rval,zval,idum
        CALL get_bfield(rval,zval,bf,1_i4)
        write(eq2_unit,*)bf%br,bf%bz,bf%f,bf%p
      ENDDO
      CLOSE(UNIT=eq1_unit)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_piegrid
c-----------------------------------------------------------------------
c     subprogram 9. separatrix.
c     Finds RZ values of separatrix
c     Based on Alan's contour finding routine (See process_eq also)
c     To avoid X points and its problems, the separatrix is defined here
c      as the flux surface that is epsln away from actual separatrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE separatrix(gt)

      TYPE(global_type), INTENT(INOUT) :: gt
c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep
      INTEGER(i4), PARAMETER :: neq=4,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq) :: y
      REAL(r8), DIMENSION(lrw) :: rwork=0
c-----------------------------------------------------------------------
c     other local variables.
c-----------------------------------------------------------------------
      INTEGER(i4), PARAMETER :: nstep=4096
      REAL(r8), PARAMETER :: eps=1.e-12, epsln=1e-4
      INTEGER(i4) :: itau
      REAL(r8) :: dr,rfac,r,z,psifac,eta,err,psi0
      REAL(r8), DIMENSION(0:nstep,0:4) :: temp
      TYPE(direct_bfield_type) :: bf
      TYPE(spline_type) :: ff
c-----------------------------------------------------------------------
c     Initialize
c-----------------------------------------------------------------------
      psifac=1. - epsln                          ! Important
      psi0=gt%psio*(1-psifac)
      r=gt%ro+SQRT(psifac)*(gt%rs2-gt%ro)
      z=gt%zo
      DO
         CALL get_bfield(r,z,bf,1_i4)
         dr=(psi0-bf%psi)/bf%psir
         r=r+dr
         IF(ABS(dr) <= eps*r)EXIT
      ENDDO
      psi0=bf%psi
c-----------------------------------------------------------------------
c     initialize variables for contour integration.
c-----------------------------------------------------------------------
      istep=0
      eta=0.
      y(1)=0.
      y(2)=SQRT((r-gt%ro)**2+(z-gt%zo)**2)
      y(3)=0.
      y(4)=0.
c-----------------------------------------------------------------------
c     set up integrator parameters.
c-----------------------------------------------------------------------
      istate=1; itask=5; iopt=1; mf=10
      itol=1; rtol=tol0; atol=tol0*y(2); rwork(1)=twopi; rwork(11)=0
c-----------------------------------------------------------------------
c     store results for each step.
c-----------------------------------------------------------------------
      DO
         rfac=y(2)
         CALL refine(rfac,eta,psi0)
         r=gt%ro+rfac*COS(eta)
         z=gt%zo+rfac*SIN(eta)
         CALL get_bfield(r,z,bf,2_i4)
         temp(istep,0)=y(1)			!   transformed angle
         temp(istep,1)=rfac**2			!   eff. radius squared
         temp(istep,2)=eta			!   real geometric angle
         temp(istep,3)=r			!   r
         temp(istep,4)=z			!   z
         err=(bf%psi-psi0)/bf%psi
c-----------------------------------------------------------------------
         IF(eta >= twopi .OR. istep >= nstep  .OR.  istate < 0
     $            .OR. ABS(err) >= 1)EXIT
         istep=istep+1
         CALL dlsode(flder,neq,y,eta,twopi,itol,rtol,atol,itask,
     $               istate,iopt,rwork,lrw,iwork,liw,jac,mf)
      ENDDO
c-----------------------------------------------------------------------
c     check for exceeding array bounds.
c-----------------------------------------------------------------------
      IF(eta < twopi)THEN
         WRITE(nim_wr,'(a,i4,a,1p,e10.3,a,i3)')
     $           ' In SUBROUTINE separatrix, istep = nstep =',nstep,
     $           ' at eta = ',eta
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     fit to cubic splines vs. theta and interpolate to tau grid.
c-----------------------------------------------------------------------
      temp(0:istep,0)=temp(0:istep,0)*twopi/y(1)
      CALL spline_alloc(ff,istep,2_i4)
      SELECT CASE(TRIM(angle_method))
      CASE('jac')
        ff%xs(0:istep)=temp(0:istep,0)
      CASE('geom')
        ff%xs(0:istep)=temp(0:istep,2)
      END SELECT
      ff%fs(0:istep,1)=temp(0:istep,3)
      ff%fs(0:istep,2)=temp(0:istep,4)
      CALL spline_fit(ff,"periodic")
      CALL spline_alloc(gt%rzsep,mtheta,2_i4)
      gt%rzsep%xs=gt%r2g%xs
      DO itau=0,mtheta
         CALL spline_eval(ff,gt%rzsep%xs(itau),0_i4)
         gt%rzsep%fs(itau,1)=ff%f(1)                                 ! R_sep
         gt%rzsep%fs(itau,2)=ff%f(2)                                 ! Z_sep
      ENDDO
      CALL spline_fit(gt%rzsep,"periodic")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE separatrix
c-----------------------------------------------------------------------
c     subprogram 10. insep.
c     Tells whether a point is inside or outside the separatrix
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION insep(gt,r,z)
      
      TYPE(global_type), INTENT(IN) :: gt
      REAL(r8), INTENT(IN) :: r,z
      LOGICAL :: insep
      REAL(r8) :: rs,rt,zs,zt, distance,distmin, crossprod, zmin,zmax
      INTEGER(i4) :: itau, it
c-----------------------------------------------------------------------
      insep = .FALSE.
      IF (r < gt%rs1 .OR. r > gt%rs2) RETURN
      distmin = 10.e6
      zmin = 10.e6
      zmax = 0.
      DO itau=0,mtheta
         IF (gt%rzsep%fs(itau,2) < zmin) zmin = gt%rzsep%fs(itau,2)
         IF (gt%rzsep%fs(itau,2) > zmax) zmax = gt%rzsep%fs(itau,2)
         distance = (r - gt%rzsep%fs(itau,1))**2 
     &            + (z - gt%rzsep%fs(itau,2))**2
         IF (distance < distmin) THEN
            distmin = distance
            rs = gt%rzsep%fs(itau,1)
            zs = gt%rzsep%fs(itau,2)
            it = itau + 1
            IF (itau == mtheta) it = 1
            rt = gt%rzsep%fs(it,1)
            zt = gt%rzsep%fs(it,2)
         ENDIF
      ENDDO
      IF (z < zmin .OR. z > zmax) RETURN
      crossprod = (rs - r)*(zt - zs) - (rt - rs)*(zs - z)
      IF (crossprod >= 0.) insep = .TRUE.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION insep

c-----------------------------------------------------------------------
      END MODULE direct
