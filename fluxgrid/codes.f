c-----------------------------------------------------------------------
c     file codes.f.
c     Interface to other stability codes
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. write_nimrod.
c     2. write_pies.
c     3. write_pest.
c     4. write_far.
c     5. farfour.
c-----------------------------------------------------------------------




c-----------------------------------------------------------------------
c     subprogram 1. write_nimrod.
c     writes flux grid data to ascii file for use by nimrod
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_nimrod(gt,iub)

      USE analyze

      TYPE(global_type), INTENT(INOUT) :: gt
      INTEGER, INTENT(IN) :: iub

      INTEGER(i4) :: itheta,ipsi
      REAL(r8) :: jac,fprime, ldelstr

      REAL(r8) :: r,z,p,br,bz,bt,jr,jz,jt, conc

c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/'ipsi = ',i3,', psi = ',1p,e11.3)
 20   FORMAT(/3x,'ith',4x,'theta',8x,'r',10x,'z',
     $     10x,'br',9x,'bz',9x,'bt'/)
 30   FORMAT(i6,1p,6e11.3)
c-----------------------------------------------------------------------
c     compute and write coordinates and equilibrium magnetic field.
c     Note that conc=0 in plasma and =1 in vacuum.  
c       ipsi<mpsi is always plasma
c-----------------------------------------------------------------------
      OPEN(UNIT=iub,FILE="fluxgrid.dat",STATUS="UNKNOWN")
      WRITE(iub,*)mpsi+mvac
      WRITE(iub,*)mtheta
      WRITE(iub,*)mvac
      WRITE(iub,*)mxpie
      WRITE(iub,*)gt%sq%fs(0,:)-gt%sq%fs1(0,:)*gt%sq%xs(0)
      WRITE(iub,*)gt%ro,gt%zo,
     $     (gt%sq%fs(0,1)-gt%sq%fs1(0,1)*gt%sq%xs(0))/gt%ro
      DO ipsi=0,mpsi
           fprime = gt%sq%fs1(ipsi,1)/gt%psio			! dF/dpsi
           WRITE(iub,*)gt%sq%fs(ipsi,:)
         DO itheta=0,mtheta
           r = gt%twod%fs(1,itheta,ipsi)
           z = gt%twod%fs(2,itheta,ipsi)
           IF (.NOT. calc_direct) THEN
             jac = gt%twod%fs(3,itheta,ipsi)
             ldelstr = delstr(itheta,ipsi) 
             IF (j_t == "use_gs") ldelstr= gsrhs(itheta,ipsi)
             br = gt%twod%fsx(1,itheta,ipsi)/jac
             bz = gt%twod%fsx(2,itheta,ipsi)/jac
             bt = gt%sq%fs(ipsi,1)/gt%twod%fs(1,itheta,ipsi)
             jr = -fprime*br/mu0
             jz = -fprime*bz/mu0
             jt = ldelstr/(mu0*r**2)			!Contravariant
           ELSE
             br = gt%dir%fs(1,itheta,ipsi)
             bz = gt%dir%fs(2,itheta,ipsi)
             bt = gt%dir%fs(3,itheta,ipsi)
             jr = gt%dir%fs(4,itheta,ipsi)
             jz = gt%dir%fs(5,itheta,ipsi)
             jt = gt%dir%fs(6,itheta,ipsi)
           ENDIF

           WRITE(iub,*)r,z,br,bz,bt,jr,jz,jt

         ENDDO
      ENDDO

c                                                               Vacuum
      IF (extr > 1. .AND. mvac > 0) THEN
        DO ipsi=1,mvac                                          !Note the 1
         DO itheta=0,mtheta
           r  = gt%vac%fs(1,itheta,ipsi)
           z  = gt%vac%fs(2,itheta,ipsi)
           br = gt%vac%fs(3,itheta,ipsi)
           bz = gt%vac%fs(4,itheta,ipsi)
           bt = gt%vac%fs(5,itheta,ipsi)
           jr = gt%vac%fs(6,itheta,ipsi)
           jz = gt%vac%fs(7,itheta,ipsi)
           jt = gt%vac%fs(8,itheta,ipsi)
           p  = gt%vac%fs(9,itheta,ipsi)
           conc  = gt%vac%fs(10,itheta,ipsi)
           WRITE(iub,*)r,z,br,bz,bt,jr,jz,jt,p,conc
         ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write eigenfunction and close output file.
c-----------------------------------------------------------------------
      IF(ASSOCIATED(gt%eigvec))THEN
         WRITE(iub,*)gt%rquot,gt%omega0,gt%growth,gt%ntor
         WRITE(iub,*)gt%eigveco
         WRITE(iub,*)gt%eigvec
      ENDIF
      CLOSE(UNIT=iub)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_nimrod


c-----------------------------------------------------------------------
c     subprogram 2. write_pies.
c     writes input file for PIES code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_pies(gt,iua)

      USE analyze

      TYPE(global_type), INTENT(IN) :: gt
      INTEGER, INTENT(IN) :: iua

      INTEGER(i4) :: itheta,ipsi
      REAL(r8) :: r, gss, ldelstr, f, fprime, psio
      REAL(r8), DIMENSION(0:mpsi) :: jdotbprof, bdotgfprof
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: jdotb, bdotgf
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      OPEN(UNIT=iua,FILE="pies.dat",STATUS="UNKNOWN")
c-----------------------------------------------------------------------
c     Calculate <J.B>/<B.Grad phi>
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     	Calc. Delstar(psi) = delsquare(psi) - 2/R Grad(R) dot Grad(psi)
c-----------------------------------------------------------------------
      psio=gt%psio
      DO ipsi=0,mpsi
          f = gt%sq%fs(ipsi,1)					! R B_tor
          fprime = gt%sq%fs1(ipsi,1)/psio			! dF/dpsi
          DO itheta=0,mtheta
            r = gt%twod%fs(1,itheta,ipsi)
            gss = gt%twod%fs(4,itheta,ipsi)
            ldelstr = delstr(itheta,ipsi) 
            IF (j_t == "use_gs") ldelstr= gsrhs(itheta,ipsi)
            jdotb(itheta,ipsi) = (ldelstr*f - fprime*gss)/r**2
            bdotgf(itheta,ipsi) = f/r**2
          ENDDO
      ENDDO

      CALL fluxav(gt,jdotb,jdotbprof)
      CALL fluxav(gt,bdotgf,bdotgfprof)
c-----------------------------------------------------------------------
c     write flux surface quantities.
c     psi, mu0 Pressure, q, 1/mu0 <J dot B> /< B dot grad phi>
c-----------------------------------------------------------------------
      WRITE(iua,*)mpsi+1
      WRITE(iua,fmt='(5(2x,E9.4))')gt%sq%xs(:)
      WRITE(iua,fmt='(5(2x,E9.4))')gt%sq%fs(:,2)/gt%psio**2
      WRITE(iua,fmt='(5(2x,E9.4))')gt%sq%fs(:,3)
      WRITE(iua,*)mpsi+1
      WRITE(iua,fmt='(5(2x,E9.4))')jdotbprof(:)/bdotgfprof(:)/mu0
c-----------------------------------------------------------------------
c     Write boundary.
c-----------------------------------------------------------------------
      WRITE(iua,*)mtheta
      WRITE(iua,fmt='(5(5x,f8.5))')gt%ob%fs(0:mtheta-1,1)	! Rbnd
      WRITE(iua,fmt='(5(5x,f8.5))')gt%ob%fs(0:mtheta-1,2)	! Zbnd
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      CLOSE(UNIT=iua)
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_pies


c-----------------------------------------------------------------------
c     subprogram 3. write_pest.
c     writes input file for PEST code.
c     REFERENCES:
c	Grimm, Dewar and Manickam, JCP 49 (1983) 94  
c	 --Has defintions of variables that are being written out.  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_pest(gt,iub)
      USE local
      USE global
      USE input
      USE physdat
      USE analyze
      IMPLICIT NONE

      TYPE(global_type), INTENT(INOUT) :: gt
      TYPE(bicube_type) :: work, grid
      TYPE(spline_type) :: intgd
      INTEGER, INTENT(IN) :: iub

      INTEGER(i4) :: itau,ipsi, mth2, nths, nx, nz, mthpi
      REAL(r8) :: f, q, xmin,xmax, zmin,zmax,alx,alz
      REAL(r8) :: bunit, lunit, punit, funit, psiunit
      REAL(r8) :: psibigunit,psi2dunit

      REAL(r8), DIMENSION(:),ALLOCATABLE :: press,pressm,mesh_half
      REAL(r8), DIMENSION(:),ALLOCATABLE :: work1
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: delta,qdelps
      REAL(r8), DIMENSION(:),ALLOCATABLE :: chi,chim
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: jac,jacm,rg,zg
      INTEGER(i4), DIMENSION(:),ALLOCATABLE :: igrid
      INTEGER(i4) :: nadres,ith,icnt
      CHARACTER(80) :: pest_title=''
c-----------------------------------------------------------------------
c     Things required for the staggered meshes
c-----------------------------------------------------------------------
      ALLOCATE(press(0:mpsi),pressm(0:mpsi),mesh_half(0:mpsi))
      ALLOCATE(work1(0:mpsi))
      ALLOCATE(chi(0:mtheta+2),chim(0:mtheta+2))
      ALLOCATE(jac(0:mtheta+2,0:mpsi),jacm(0:mtheta+2,0:mpsi))
      ALLOCATE(rg(0:mtheta+2,0:mpsi),zg(0:mtheta+2,0:mpsi))
c-----------------------------------------------------------------------
c     Define psi and pressure half mesh
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi-1
        mesh_half(ipsi)=(gt%sq%xs(ipsi)+gt%sq%xs(ipsi+1))*0.5
        CALL spline_eval(gt%sq,mesh_half(ipsi),0_i4)
        pressm(ipsi)=gt%sq%f(2)
      ENDDO
      mesh_half(mpsi)=gt%sq%xs(mpsi)
      CALL spline_eval(gt%sq,mesh_half(mpsi),0_i4)
      pressm(mpsi)=gt%sq%f(2)
c-----------------------------------------------------------------------
c	PEST uses an abnormal grid: For x(i,j) where i is the
c	 theta coordinate, increasing i corresponds to going CLOCKWISE
c	 with the 3rd grid point corresponding to the outboard plane
c	 so that there is some overlap
c-----------------------------------------------------------------------
      mth2 = mtheta+2
      nths = mtheta+5
      CALL bicube_alloc(grid,nths,mpsi,6_i4)
      ALLOCATE(igrid(0:nths))
      igrid(0:3)=( (/3,2,1,0/) )
      DO itau=4,mtheta+3
          igrid(itau) = mtheta - itau + 3
      ENDDO
      igrid(mtheta+4) = mtheta-1
      igrid(mtheta+5) = mtheta-2

      grid%fs = 0
      grid%ys = gt%twod%ys                              ! radial coord.
      DO itau=0,nths					! theta coord.
        grid%xs(itau) = gt%twod%xs(igrid(itau))
      ENDDO
      grid%xs(0:3) = -grid%xs(0:3)
      grid%xs(4:nths-2) = twopi-grid%xs(4:nths-2)
      grid%xs(nths-1:nths) = twopi+ABS(twopi-grid%xs(nths-1:nths))

      DO ipsi=0,mpsi
        DO itau=0,nths
          grid%fs(:,itau,ipsi) = gt%twod%fs(:,igrid(itau),ipsi)
        ENDDO
      ENDDO
      CALL bicube_fit(grid,"extrap","extrap")
c-----------------------------------------------------------------------
c       grid node spacing for the x,z box  - PEST3 doesn't use but 
c	other versions of PEST do (for printing packages)
c-----------------------------------------------------------------------
      mthpi = AINT(REAL(mtheta)/2.) + 3
      nx = 100						! Somewhat 
      nz = 100						!  arbitrary
      xmin = MINVAL(grid%fs(:,:,1))
      xmax = MAXVAL(grid%fs(:,:,1))
      zmin = MINVAL(grid%fs(:,:,2))
      zmax = MAXVAL(grid%fs(:,:,2))
      alx = (xmax-xmin)/REAL(nx,r8)
      alz = (zmax-zmin)/REAL(nz,r8)
c-----------------------------------------------------------------------
c     Define units to normalize quantities
c-----------------------------------------------------------------------
      lunit=gt%rmean
      funit=gt%bt0*lunit
      bunit=funit/lunit
      psibigunit=funit*lunit/(twopi*gt%psio)
      psiunit=funit*lunit/gt%psio
      psi2dunit=funit*lunit		! Metric elements already have gt%psio
      punit=bunit**2
c-----------------------------------------------------------------------
c     Calculate needed stuff -> X^2 and q delta
c-----------------------------------------------------------------------
      CALL bicube_alloc(work,nths,mpsi,2_i4)
      CALL spline_alloc(intgd,nths,1_i4)
      work%ys = grid%ys					! radial coord.
      work%xs = grid%xs					! theta
      intgd%xs= grid%xs
      work%fs(1,:,:) = grid%fs(1,:,:)**2		! R^2
      DO ipsi=0,mpsi
        f = gt%sq%fs(ipsi,1)
        q = gt%sq%fs(ipsi,3) 
        intgd%fs(:,1) = f/q*grid%fs(3,:,ipsi)
     &       /work%fs(1,:,ipsi) - 1._r8			! F/q jac/X^2 - 1
        CALL spline_fit(intgd,"extrap")
        CALL spline_int(intgd)
        work%fs(2,:,ipsi) = intgd%fsi(:,1)*q		! delta*q
      ENDDO
      CALL bicube_fit(work,"extrap","extrap")
      ALLOCATE(qdelps(0:nths,0:mpsi))
      qdelps(0:nths,0:mpsi) = work%fsy(2,0:nths,0:mpsi)	! QDELPS

      DO ith=0,nths
        DO itau=1,3
          DO ipsi=1,mpsi-1
            work%fsy(2,ith,ipsi)= qdelps(ith,ipsi  )*0.50+
     &                            qdelps(ith,ipsi-1)*0.25+
     &                            qdelps(ith,ipsi+1)*0.25
          ENDDO
          DO ipsi=1,mpsi-1
            qdelps(ith,ipsi)= work%fsy(2,ith,ipsi  )*0.50+
     &                        work%fsy(2,ith,ipsi-1)*0.25+
     &                        work%fsy(2,ith,ipsi+1)*0.25
          ENDDO
        ENDDO
      ENDDO
 
 
      grid%fs(:,mtheta+5:mtheta+5,:) = 0		! Match Alex 
      grid%fsx(:,mtheta+5:mtheta+5,:) = 0		! Match Alex 
      grid%fsy(:,mtheta+5:mtheta+5,:) = 0		! Match Alex 
 
 							! Calculate DELTA
c       In pest coordinates delta is exactly zero, but computing it
c	introduces roundoff errors which the pest3 routines cannot fit.
      ALLOCATE(delta(1:nths,0:mpsi))
      IF((ipb == 0) .AND. (ipr == 2))THEN
        delta = 0
        qdelps = 0
        work%fsy(2,:,:)=0
      ELSE
        DO ipsi=0,mpsi
          q = gt%sq%fs(ipsi,3) 
          delta(1:nths,ipsi)=work%fs(2,1:nths,ipsi)/q
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c	Open and Write MPOUT1
c-----------------------------------------------------------------------
      OPEN(UNIT=iub,FILE='MPOUT1',ACCESS='direct',RECL=r8)
      nadres=1
      DO ipsi=1,20
         ith=(ipsi-1)*r4
         WRITE(iub,REC=ipsi) pest_title(ith:ith+r4)     ! PEST Variable Name
      ENDDO
      WRITE(iub,REC=21) 0._r8				! DAT
      WRITE(iub,REC=22) nx    				! NX 		nxx(1)
      WRITE(iub,REC=23) nz    				! NZ 		nxx(2)
      WRITE(iub,REC=24) mpsi+1    			! NOSURF 	nxx(3)
      WRITE(iub,REC=25) mtheta    			! MTH 		nxx(4)
      WRITE(iub,REC=26) 0    				! LPLESS 	nxx(5)
      WRITE(iub,REC=27) alx/lunit      			! ALX 		axx(1)
      WRITE(iub,REC=28) alz/lunit      			! ALZ 		axx(2)
      WRITE(iub,REC=29) 0._r8				! XZERO 	axx(3)
      WRITE(iub,REC=30) gt%ro/lunit			! XMA 		axx(4)
      WRITE(iub,REC=31) gt%rmean/lunit			! R  		axx(5)
      WRITE(iub,REC=32) gt%sq%fs(0,2)/punit		! P0 		axx(6)
      WRITE(iub,REC=33) 0._r8				! GP 		axx(7)
      WRITE(iub,REC=34)-1.0/psibigunit			! PSIMIN 	axx(8)
      WRITE(iub,REC=35) 0._r8				! PSILIM(by def.)axx(9)
      WRITE(iub,REC=36) 0._r8				! PSIPLS(by def.)axx(10)
      WRITE(iub,REC=37) 0._r8				! BETAP 	axx(11)
      WRITE(iub,REC=38) 0._r8				! BETAG 	axx(12)
      WRITE(iub,REC=39) 1._r8				! UPSILON 	axx(13)
      WRITE(iub,REC=40) ipb				! NJ2		nxy(1)
      WRITE(iub,REC=41) ipr				! MJ2		nxy(2)
      DO ipsi=0,7
         WRITE(iub,REC=42+ipsi)0 			!Records 42-49 nxy(3:10)
      ENDDO
c-----------------------------------------------------------------------
c     write flux surface quantities with normalizations
c-----------------------------------------------------------------------
      nadres=50
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1) pressm(ipsi)/punit		! p
      ENDDO
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1) 
     &                 gt%sq%fs1(ipsi,2)*psiunit/punit			! p'
      ENDDO

      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1) gt%sq%fs(ipsi,3)		! q
      ENDDO
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1)
     &                    gt%sq%fs1(ipsi,3)*psiunit			! q'
      ENDDO

      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1)
     &                     gt%sq%fs(ipsi,1)/funit			! F
      ENDDO
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1)
     &                  gt%sq%fs1(ipsi,1)/funit*psiunit			! F'
      ENDDO
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1)
     &          gt%sq%fs(ipsi,1)/gt%sq%fs(ipsi,3)/funit			! F/q
      ENDDO
      work1(:)=(gt%sq%fs1(:,1)/gt%sq%fs(:,3)	
     &     - gt%sq%fs(:,1)*gt%sq%fs1(:,3)/gt%sq%fs(:,3)**2 )
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1) 
     &                        work1(ipsi)/funit*psiunit		! (F/q)'
      ENDDO
c              (The stupid thing is the next quantity gets overwritten.)
      nadres = nadres + mpsi + 1
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1) dideal(ipsi)		! di
      ENDDO
c-----------------------------------------------------------------------
c	Write out 2D quantities - splines are defined 0:nths but PEST
c	is going from 1:nths
c-----------------------------------------------------------------------
      nadres = 50 + 8*(mpsi+1)
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                        grid%fs(1,ith,ipsi)/lunit		! R
      ENDDO
      ENDDO
      nadres=nadres+(mpsi+1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                              grid%fs(2,ith,ipsi)/lunit	! Z
      ENDDO
      ENDDO
      nadres=nadres+(mpsi+1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &                       grid%fsx(1,ith,ipsi)/lunit		! XDTH d/dth(R)
      ENDDO
      ENDDO
      nadres=nadres+(mpsi+1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                             grid%fsx(2,ith,ipsi)/lunit	! ZDTH d/dth(Z)
      ENDDO
      ENDDO
      nadres=nadres+(mpsi+1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                grid%fsy(1,ith,ipsi)*psiunit/lunit   ! XDPS  d/ds(R) 
      ENDDO
      ENDDO
      nadres=nadres+(mpsi+1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &                  grid%fsy(2,ith,ipsi)*psiunit/lunit  ! ZDPS d/ds(Z)
      ENDDO
      ENDDO
c	Close MPOUT1 and move on
c-----------------------------------------------------------------------
      CLOSE(iub)

c-----------------------------------------------------------------------
c	Open and Write MAPDSK
c-----------------------------------------------------------------------
      OPEN(UNIT=iub,FILE='MAPDSK',ACCESS='direct',RECL=r8)
c-----------------------------------------------------------------------
c	Write out dimension and other info 
c	Things no longer needed are set to 0
c-----------------------------------------------------------------------
      nadres=1
      DO ipsi=1,20
         ith=(ipsi-1)*r4
         WRITE(iub,REC=ipsi) pest_title(ith:ith+r4)     ! PEST Variable Name
      ENDDO
      WRITE(iub,REC=21) 0._r8				! DAT
      WRITE(iub,REC=22) nx				! NX 		nxx(1)
      WRITE(iub,REC=23) nz				! NZ 		nxx(2)
      WRITE(iub,REC=24) mpsi+1	  			! NOSURF 	nxx(3)
      WRITE(iub,REC=25) mtheta    			! MTH 		nxx(4)
      WRITE(iub,REC=26) 0				! LPLESS 	nxx(5)
      WRITE(iub,REC=27) alx/lunit			! ALX 		axx(1)
      WRITE(iub,REC=28) alz/lunit			! ALZ 		axx(2)
      WRITE(iub,REC=29) 0._r8				! XZERO 	axx(3)
      WRITE(iub,REC=30) gt%ro/lunit			! XMA 		axx(4)
      WRITE(iub,REC=31) gt%rmean/lunit			! R  		axx(5)
      WRITE(iub,REC=32) gt%sq%fs(0,2)/punit		! P0 		axx(6)
      WRITE(iub,REC=33) 0._r8				! GP 		axx(7)
      WRITE(iub,REC=34)-1.0/psibigunit			! PSIMIN 	axx(8)
      WRITE(iub,REC=35) 0._r8				! PSILIM(by def.)axx(9)
      WRITE(iub,REC=36) 0._r8				! PSIPLS(by def.)axx(10)
      WRITE(iub,REC=37) 0._r8				! BETAP 	axx(11)
      WRITE(iub,REC=38) 0._r8				! BETAG 	axx(12)
      WRITE(iub,REC=39) 1._r8				! UPSILON 	axx(13)
      WRITE(iub,REC=40) ipb				! NJ2		nxy(1)
      WRITE(iub,REC=41) ipr				! MJ2		nxy(2)
      DO ipsi=0,7
         WRITE(iub,REC=42+ipsi)0 			!Records 42-49 nxy(3:10)
      ENDDO
c-----------------------------------------------------------------------
c     Now write the rest of the stuff
c-----------------------------------------------------------------------
c  							outer boundary
      nadres = 50
      DO ith=1,nths-3
        WRITE(iub,REC=(nadres-1)*r8+ith) grid%fs(1,ith,mpsi)/lunit	! XINF 
      ENDDO
      nadres = nadres + nths-3
      DO ith=1,nths-3
        WRITE(iub,REC=(nadres-1)*r8+ith) grid%fs(2,ith,mpsi)/lunit	! ZINF
      ENDDO
      nadres = nadres + nths-3
      DO ipsi=0,mpsi
        WRITE(iub,REC=(nadres-1)*r8+ipsi+1)
     &                (gt%sq%xs(ipsi)-1.0)/psibigunit      	! PSIBIG
      ENDDO
      nadres = nadres + mpsi + 1
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &            grid%fs(4,ith,ipsi)*(lunit/psi2dunit)**2   	! GRPSSQ  (gss)
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                  grid%fs(1,ith,ipsi)**2/lunit**2 	! XSQ (R^2)
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &                         grid%fs(6,ith,ipsi)*lunit**2	! GRTHSQ  (gtt)
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &                -grid%fs(5,ith,ipsi)/psi2dunit*lunit**2	! GRPSTH  (gst)
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &       -grid%fsx(5,ith,ipsi)/psi2dunit*lunit**2	! GPTDTH  (d/dth(gst))
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &   work%fsy(1,ith,ipsi)/psiunit/lunit**2		! XSQDPS  (d/ds(R^2))
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &        grid%fsx(4,ith,ipsi)*(lunit/psi2dunit)**2	! GPSDTH  (d/dth(gss))
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &             work%fsx(1,ith,ipsi)/lunit**2 	! XSQDTH  (d/dt(R^2))
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &		grid%fs(3,ith,ipsi)/lunit**3*psi2dunit	! XJACOB  (jacobian)
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &   grid%fsy(3,ith,ipsi)/lunit**3*psiunit*psi2dunit ! XJPRYM  (d/ds(jac))
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt) 
     &                           delta(ith,ipsi)   	! DELTA   
      ENDDO
      ENDDO
      nadres = nadres + (mpsi + 1)*nths
      icnt=0
      DO ipsi=0,mpsi
      DO ith=1,nths
        icnt=icnt+1
        WRITE(iub,REC=(nadres-1)*r8+icnt)
     &    qdelps(ith,ipsi)*psiunit     	! QDELPS (d/ds(q delta))
!    &    work%fsy(2,ith,ipsi)*psiunit     	! QDELPS (d/ds(q delta))
      ENDDO
      ENDDO
c 					DELT is read by vacuum.f90
      DEALLOCATE(work1)
      ALLOCATE(work1(1:nths))
      nadres = nadres + (mpsi + 1)*nths
      work1=0.
      DO ith=1,nths
        WRITE(iub,REC=(nadres-1)*r8+ith)work1(ith)
      ENDDO
c-----------------------------------------------------------------------
c     Deallocate arrays and close output file
c-----------------------------------------------------------------------
      CALL spline_dealloc(intgd)
      CALL bicube_dealloc(work)
      DEALLOCATE(press,pressm,mesh_half,chi,chim,jac,jacm,rg,zg)
      DEALLOCATE(work1,delta,qdelps)
      CLOSE(UNIT=iub)
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_pest


c-----------------------------------------------------------------------
c     subprogram 4. write_far.
c     Subroutine to output FAR input file.  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_far(gt,iub)
      USE local
      USE global
      USE input
      USE physdat
      USE analyze
      IMPLICIT NONE

      INTEGER(i4), PARAMETER :: modes = 10 		! # of modes
      INTEGER(i4), INTENT(IN) :: iub
      TYPE(global_type), INTENT(IN) :: gt

      INTEGER(i4) :: i,j, idum, jdum

      REAL(r8) ::  anorm, eps, gee			! Conversion vars
      REAL(r8) ::  rstl,rstpsi,rstpr,rstf, rstb		! Conversion vars
      REAL(r8) ::  rstr,rstc1,rstc4,rstd1		! Conversion vars
      REAL(r8) :: fac1, fac2 				! Extrap. vars
      REAL(r8) :: arg,dum,drx

      REAL(r8), DIMENSION(0:mtheta) :: rint
      REAL(r8), DIMENSION(0:mpsi) :: rrst, dpsidr, qrst, psirst, mach
      REAL(r8), DIMENSION(0:mpsi) :: prst, pprst, frst, ffprst 

      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: quantity, grr, grt
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: xmap, bmod, brho, bthe
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: bigri, ccon1,ccon2,ccon3
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: ccon4,ccon5,ccon6
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: ddiv1,ddiv2,ddiv3,expfac
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: xrst, zrst, xrrst, zrrst
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: xtrst, ztrst
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: grrrst, grtrst, gttrst
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: rinv, exprst
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: con1,con2,con3,con4
      REAL(r8), DIMENSION(0:mpsi,1:modes) :: con5,con6,div1,div2,div3
c-----------------------------------------------------------------------
c     FAR is based on the PEST coordinate, so check for that
c-----------------------------------------------------------------------
      IF(.NOT.((ipb == 0).AND.(ipr == 2)) .AND. (angle_method /= 'jac'))
     &  CALL nim_stop("write_far:  Use PEST coordinates (ipb=0, ipr=2)")
c-----------------------------------------------------------------------
c     Define radial coordinate rrst. Calc. anorm, eps for norms.
c	rrst = Sqrt(Int(q/F,psimin,psi))
c-----------------------------------------------------------------------
      DO i=0,mpsi
        rrst(i)= 0.0
        rint(0)=gt%twod%fs(2,0,i)/gt%twod%fs(1,0,i)
        DO j=1, mtheta
          rint(j) = gt%twod%fs(2,j,i)/gt%twod%fs(1,j,i)
          drx = gt%twod%fs(1,j-1,i) - gt%twod%fs(1,j,i)
          rrst(i)=rrst(i) + 0.5*(rint(j-1)+rint(j) )*drx
        ENDDO
        rrst(i) = 2.0*pi*rrst(i)
      ENDDO
         gee = 4.*pi**2/rrst(mpsi)
      DO i=0,mpsi
           rrst(i) = (rrst(i)/rrst(mpsi))**0.5
      ENDDO
      anorm= (2. * gt%rmean/gee )**0.5 
      eps = anorm/gt%rmean

c								dpsidr=rF/q
      dpsidr(:)=rrst(:)*gt%sq%fs(:,1)/gt%sq%fs(:,3)/gt%psio
     &  *eps**2*gt%rmean
c-----------------------------------------------------------------------
c     Define some useful stuff 
c-----------------------------------------------------------------------
      fac1=(gt%sq%xs(2))**2/((gt%sq%xs(2))**2 - (gt%sq%xs(1))**2)
      fac2= -(gt%sq%xs(1))**2/((gt%sq%xs(2))**2 - (gt%sq%xs(1))**2)
      idum=0.0
      dum=0.0

c	Conversion Vars.				Quantities
c       ----------------				----------
      rstl = eps * gt%rmean				! length
      rstf = gt%sq%fs(mpsi,1)				! F (a.k.a. fwall)
      rstb = rstf/ gt%rmean				! B
      rstpsi = -rstb * rstl**2				! psi
      rstpr = gt%sq%fs(0,2) 				! pressure
      rstr= eps**2 * gt%rmean                       	! rho
      rstc1= 1.0/(gt%rmean**2 * rstb)			! con1, con2, con3
      rstc4= 1.0/(gt%rmean * rstb)			! con4, con5, con6
      rstd1= 1.0/(eps**2 * gt%rmean**2 * rstb)		! div1, div2, div3

      psirst(:) = (gt%sq%xs(:) - gt%sq%xs(mpsi) )*gt%psio/rstpsi
      frst(:) =   gt%sq%fs(:,1)/rstf
      prst(:) =   gt%sq%fs(:,2)/rstpr
      qrst(:) =   gt%sq%fs(:,3)
      mach(:) =   gt%sq%fs(:,4)
      pprst(:) =   gt%sq%fs1(:,2)/gt%psio/rstpr*rstpsi
      ffprst(:) =  frst(:) *gt%sq%fs1(:,1)/gt%psio/rstf*rstpsi

         DO i=0,mpsi
            grr(:,i) =gt%twod%fs(4,:,i)/gt%psio**2/dpsidr(i)**2
     &            *rstl**2			!grr=R^2 gss*(dr/ds)^2
     &            *(gt%twod%fs(1,:,i)/gt%rmean)**2
         ENDDO
         DO i=1,mpsi
          grt(:,i)=gt%twod%fs(5,:,i)/gt%psio*rrst(i)/dpsidr(i)
     &            *rstl**2			!grt=R^2 r gst*dr/ds
     &            *(gt%twod%fs(1,:,i)/gt%rmean)**2
         ENDDO
         grt(:,0) = fac1*grt(:,1) + fac2*grt(:,2)
c-----------------------------------------------------------------------
c	Calculate the Fourier coefficients for various quantities
c         Remember: 1=cos, 0 = sin
c-----------------------------------------------------------------------
       quantity(:,:) = (gt%twod%fs(1,:,:) - gt%rmean)/rstl
       CALL farfour(gt,quantity, xrst, modes, 1_i4)		!X

       quantity(:,:) = gt%twod%fs(2,:,:)/rstl
       CALL farfour(gt,quantity, zrst, modes, 0_i4)		!Z

       quantity(:,:)=gt%twod%fsy(1,:,:)
       CALL farfour(gt,quantity, xrrst, modes, 1_i4)		!dX/dr
       DO i=0,mpsi
           xrrst(i,:) = xrrst(i,:)*dpsidr(i)/rstl
       ENDDO

       quantity(:,:)=gt%twod%fsy(2,:,:)
       CALL farfour(gt,quantity, zrrst, modes, 0_i4)		!dZ/dr
       DO i=0,mpsi
           zrrst(i,:) = zrrst(i,:)*dpsidr(i)/rstl
       ENDDO

       quantity(:,:) = gt%twod%fsx(1,:,:)
       CALL farfour(gt,quantity, xtrst, modes, 0_i4)		!dX/dth
       DO i=0,mpsi
           xtrst(i,:) = xtrst(i,:)/rrst(i)/rstl
       ENDDO

       quantity(:,:) = gt%twod%fsx(2,:,:)
       CALL farfour(gt,quantity, ztrst, modes, 1_i4)		!dZ/dth
       DO i=0,mpsi
           ztrst(i,:) = ztrst(i,:)/rrst(i)/rstl
       ENDDO

       CALL farfour(gt,grr, grrrst, modes, 1_i4)
       CALL farfour(gt,grt, grtrst, modes, 0_i4)

       DO i=1,mpsi
          quantity(:,i) = gt%twod%fs(6,:,i)*rrst(i)**2
     &            *rstl**2
     &            *(gt%twod%fs(1,:,i)/gt%rmean)**2
       ENDDO
       quantity(:,0) = fac1*quantity(:,1) + fac2*quantity(:,2)
       CALL farfour(gt,quantity, gttrst, modes, 1_i4)		!gtt=gtt/(rR)^2

       xrrst(0,:) = fac1*xrrst(1,:) + fac2*xrrst(2,:)
       xrrst(0,:) = fac1*xrrst(1,:) + fac2*xrrst(2,:)
       ztrst(0,:) = fac1*ztrst(1,:) + fac2*ztrst(2,:)
       ztrst(0,:) = fac1*ztrst(1,:) + fac2*ztrst(2,:)
c-----------------------------------------------------------------------
c	Write out in format that FAR reads in
c 	dum, idum used for rsteq data the FAR doesn't really need
c-----------------------------------------------------------------------
       OPEN(UNIT=iub,FILE='fareq',STATUS='unknown',FORM='unformatted')
          WRITE(iub) idum,idum,idum,dum,dum,dum,dum,dum,dum,dum,idum
          WRITE(iub)dum, dum, dum, dum, idum
          WRITE(iub) mpsi+1, modes+1, gt%betat, eps
          WRITE(iub) 
     &               (rrst(i),        i=0,mpsi),
     &               (psirst(i),      i=0,mpsi),
     &               (frst(i),        i=0,mpsi),
     &               ((xrst(i,j),     i=0,mpsi),j=1,modes),
     &               ((zrst(i,j),     i=0,mpsi),j=1,modes)
          WRITE(iub) 
     &               (ffprst(i),      i=0,mpsi),
     &               ((xrrst(i,j),    i=0,mpsi),j=1,modes), 
     &               ((xtrst(i,j),    i=0,mpsi),j=1,modes),
     &               ((zrrst(i,j),    i=0,mpsi),j=1,modes),
     &               ((ztrst(i,j),    i=0,mpsi),j=1,modes),
     &               ((grrrst(i,j),   i=0,mpsi),j=1,modes),
     &               ((grtrst(i,j),   i=0,mpsi),j=1,modes),
     &               ((gttrst(i,j),   i=0,mpsi),j=1,modes)
          idum=9800
          jdum=6               		! Use these values.  It's easier
          WRITE(iub) idum, jdum
          WRITE(iub) (qrst(i),  i=0,mpsi)
          WRITE(iub) (prst(i),  i=0,mpsi)
          WRITE(iub) (pprst(i), i=0,mpsi)
       CLOSE(iub)
c-----------------------------------------------------------------------
c     Calculate the neoeq quantities 
c	Metric elements here have 1/R^2 in them
c     Note:            dF/drho(i) = -ffprst(i)*psic(i)/qrst(i)
c-----------------------------------------------------------------------
        xmap(:,:) = gt%twod%fs(1,:,:)

        DO i=0,mpsi
          DO j=0,mtheta
           arg=5._r8/6._r8*mach(i)**2*( (xmap(j,i)/gt%ro)**2-1.)
           expfac(j,i)= exp(arg)
            bmod(j,i)= SQRT(gt%sq%fs(i,1)**2 + gt%twod%fs(4,j,i))/
     &                       xmap(j,i)
            brho(j,i) = 0.5/bmod(j,i)/xmap(j,i)**2*(
     &           2*gt%sq%fs(i,1)*gt%sq%fs1(i,1) 
     &         + gt%twod%fsy(4,j,i)
     &         - 2.*(gt%sq%fs(i,1)**2+gt%twod%fs(4,j,i))
     &            *gt%twod%fsy(1,j,i)/xmap(j,i)   )*dpsidr(i)
            bthe(j,i) = 0.5/bmod(j,i)/xmap(j,i)**2*(
     &         + gt%twod%fsx(4,j,i)
     &         - 2.*(gt%sq%fs(i,1)**2+gt%twod%fs(4,j,i))
     &            *gt%twod%fsx(1,j,i)/xmap(j,i)   )
          ENDDO
        ENDDO

       DO i=0,mpsi
         DO j=0,mtheta
           bigri(j,i) = 1.0/xmap(j,i)
    
           ccon1(j,i) = 1.0/( xmap(j,i)**2 * bmod(j,i)   )
           ccon2(j,i) = bthe(j,i)/(xmap(j,i)**2*bmod(j,i)**2)
           ccon3(j,i) = brho(j,i)/(xmap(j,i)**2*bmod(j,i)**2)

           ccon4(j,i) = (frst(i)*rstf)/(xmap(j,i)**2*bmod(j,i)**2)
           ccon5(j,i) = (rrst(i)*frst(i)*rstf*grr(j,i)*rstl**2 )/
     &                     (qrst(i) * xmap(j,i)**4 * bmod(j,i)**2 )
           ccon6(j,i) = (rrst(i)*frst(i)*rstf*grt(j,i)*rstl**2 )/
     &                     (qrst(i) * xmap(j,i)**4 * bmod(j,i)**2 )

           ddiv1(j,i) = rstf/(rstr*xmap(j,i)**2) * (
     &              ffprst(i)*rrst(i)/qrst(i) / bmod(j,i)**2
     &                + frst(i)/bmod(j,i)**3 * 2.0 * brho(j,i) )

           ddiv2(j,i) = -1.0/bmod(j,i)**2*( 
     &          (pprst(i)*rstpr 
     &               + ffprst(i)/xmap(j,i)**2*rstf**2)/rstpsi
     &            + rrst(i)/qrst(i)*frst(i)*rstf*
     &                 2.0/bmod(j,i)/xmap(j,i)**4/rstr*
     &      rstl**2*(grr(j,i)*brho(j,i)+grt(j,i)/rrst(i)*bthe(j,i)))

           ddiv3(j,i) = (frst(i)*rstf)/(rrst(i)*rstr) *
     &            2.0/xmap(j,i)**2 * bthe(j,i)/bmod(j,i)**3

c							Normalize
           bigri(j,i) = bigri(j,i)*gt%rmean
           ccon1(j,i) = ccon1(j,i)/rstc1
           ccon2(j,i) = ccon2(j,i)/rstc1
           ccon3(j,i) = ccon3(j,i)/rstc1
           ccon4(j,i) = ccon4(j,i)/rstc4
           ccon5(j,i) = ccon5(j,i)/rstc4
           ccon6(j,i) = ccon6(j,i)/rstc4
           ddiv1(j,i) = ddiv1(j,i)/rstd1
           ddiv2(j,i) = ddiv2(j,i)/rstd1
           ddiv3(j,i) = ddiv3(j,i)/rstd1
         ENDDO
       ENDDO

c   Calculate Fourier coefficients
c       						1=cos, 0 = sin
         CALL farfour(gt,bigri, rinv, modes, 1_i4)
         CALL farfour(gt,ccon1, con1, modes, 1_i4)
         CALL farfour(gt,ccon2, con2, modes, 0_i4)
         CALL farfour(gt,ccon3, con3, modes, 1_i4)
         CALL farfour(gt,ccon4, con4, modes, 1_i4)
         CALL farfour(gt,ccon5, con5, modes, 1_i4)
         CALL farfour(gt,ccon6, con6, modes, 0_i4)
         CALL farfour(gt,ddiv1, div1, modes, 1_i4)
         CALL farfour(gt,ddiv2, div2, modes, 1_i4)
         CALL farfour(gt,ddiv3, div3, modes, 0_i4)
         CALL farfour(gt,expfac, exprst, modes, 1_i4)

       rinv(0,1)=fac1*rinv(1,1)+fac2*rinv(2,1)
       con1(0,1)=fac1*con1(1,1)+fac2*con1(2,1)
       con3(0,2)=fac1*con3(1,2)+fac2*con3(2,2)
       con4(0,1)=fac1*con4(1,1)+fac2*con4(2,1)
       div1(0,2)=fac1*div1(1,2)+fac2*div1(2,2)
       div2(0,1)=fac1*div2(1,1)+fac2*div2(2,1)
       div3(0,2)=fac1*div3(1,2)+fac2*div3(2,2)
       exprst(0,1)=fac1*exprst(1,1)+fac2*exprst(2,1)
c-----------------------------------------------------------------------
c     Write out neoeq
c-----------------------------------------------------------------------
      OPEN(UNIT=iub,FILE='neoeq',STATUS='unknown',FORM='unformatted')
      WRITE(iub) mpsi, modes, gt%betat, eps, gt%rmean, gt%bt0
      WRITE(iub) 
     &    ((rinv(i,j),i=0,mpsi),j=0,modes), 
     &    ((con1(i,j),i=0,mpsi),j=0,modes), 
     &    ((con2(i,j),i=0,mpsi),j=0,modes), 
     &    ((con3(i,j),i=0,mpsi),j=0,modes),
     &    ((con4(i,j),i=0,mpsi),j=0,modes),
     &    ((con5(i,j),i=0,mpsi),j=0,modes),
     &    ((con6(i,j),i=0,mpsi),j=0,modes),
     &    ((div1(i,j),i=0,mpsi),j=0,modes),
     &    ((div2(i,j),i=0,mpsi),j=0,modes),
     &    ((div3(i,j),i=0,mpsi),j=0,modes),
     &    ((exprst(i,j),i=0,mpsi),j=0,modes),
     &    (mach(i),i=0,mpsi)
       CLOSE(iub)
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_far




c-----------------------------------------------------------------------
c     subprogram 5. farfour.
c     Calculates the Fourier coefficients for FAR
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c     Flag says whether to use sin's (0) or cos's (1)
c-----------------------------------------------------------------------
      SUBROUTINE farfour(gt,quantity,fourcoef,modes,flag)
      USE global
      USE local
      USE input

      TYPE(global_type), INTENT(IN) :: gt
      INTEGER(i4), INTENT(IN) :: flag, modes
      REAL(r8), DIMENSION(0:mtheta,0:mpsi), INTENT(IN) :: quantity
      REAL(r8), DIMENSION(0:mpsi,1:modes), INTENT(OUT) :: fourcoef

      TYPE(spline_type) :: intgd
      INTEGER(i4) :: ipsi,md, m
c-----------------------------------------------------------------------
c     Integrate over theta for each flux surface
c-----------------------------------------------------------------------
      CALL spline_alloc(intgd,mtheta,modes)
      intgd%xs= gt%twod%xs
      DO ipsi=0,mpsi
        DO md = 1,modes
          m = real(md - 1)
          IF (flag.eq.1) THEN
            intgd%fs(:,md)= quantity(:,ipsi)*cos(m*intgd%xs(:)) 
          ELSE
            intgd%fs(:,md)= quantity(:,ipsi)*sin(m*intgd%xs(:)) 
          ENDIF
        ENDDO
        CALL spline_fit(intgd,"periodic")
        CALL spline_int(intgd)
        fourcoef(ipsi,:)=intgd%fsi(mtheta,:)/pi
      ENDDO
      CALL spline_dealloc(intgd)
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE farfour
