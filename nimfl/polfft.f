      MODULE polfft_mod
      USE local
      USE rfft_mod
      USE bicube
      IMPLICIT NONE
      REAL(r8), DIMENSION(:), ALLOCATABLE :: theta0
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rho0
      REAL(r8), DIMENSION(:), ALLOCATABLE :: vol0
      REAL(r8), DIMENSION(:), ALLOCATABLE :: q0
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: contour
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: transform
      TYPE(bicube_type) :: r2g
      END MODULE polfft_mod

      SUBROUTINE polfft
c					Calculation of poloidal fft's
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE start_positions
      USE polfft_mod
c     USE fft_mod
      IMPLICIT NONE

      LOGICAL :: local_debug=.FALSE.
      TYPE(cell_type), POINTER :: item,item_pre
      INTEGER(i4) :: ibl,iline,ipoloidal,mmcell,nqty,num_fourier
      INTEGER(i4) :: iqty,imode,im
      REAL(r8) :: x_start,y_start,jacobian,rval
      LOGICAL :: failure
      TYPE(location_type) :: p0
      COMPLEX(r8),DIMENSION(:,:,:),ALLOCATABLE :: bl_rho,bu_rho
      COMPLEX(r8),DIMENSION(:,:,:),ALLOCATABLE :: bl_the,bu_the
      COMPLEX(r8),DIMENSION(:,:,:),ALLOCATABLE :: bl_zed,bu_zed
      COMPLEX(r8),DIMENSION(:,:,:,:),ALLOCATABLE :: bu_nr,bl_nr
      COMPLEX(r8),DIMENSION(:,:,:,:),ALLOCATABLE :: bu_ni,bl_ni
      CHARACTER(72) :: zone_title
      CHARACTER(1) :: single_quote='"'
      REAL(r8) :: area_avg,norm_m,norm_n,vol_avg
      REAL(r8), DIMENSION(:),ALLOCATABLE  :: D,emsum

      item => start
      item_pre => start
c
c						Allocate storage 
      mmcell = 2*mm+1
      nqty=(SIZE(rbc(1)%be%f,1))
      num_fourier=(SIZE(rbc(1)%be%f,2))
      ALLOCATE(emsum(num_fourier))
      ALLOCATE(bl_rho(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bu_rho(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bl_the(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bu_the(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bl_zed(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bu_zed(num_fourier,0:mmcell-1,1:n_fieldlines))
      ALLOCATE(bl_nr(3,num_fourier,0:mm,1:n_fieldlines))
      ALLOCATE(bu_nr(3,num_fourier,0:mm,1:n_fieldlines))
      ALLOCATE(bl_ni(3,num_fourier,0:mm,1:n_fieldlines))
      ALLOCATE(bu_ni(3,num_fourier,0:mm,1:n_fieldlines))
c
c					Allocate work arrays for ffts
      nfft = mmcell
      ALLOCATE(WA(nfft))
      ALLOCATE(CH(nfft))
      ALLOCATE(C(nfft))
      WA = 0
      CH = 0
      C  = 0
      CALL rffti
      ALLOCATE(D(nfft))
      D  = 0
c
c					vol0=>volume/2*mu0

      IF (geom=='tor') THEN
        vol0=vol0/(3*4.e-7)			!  3 for div(r) in 3D
      ELSE
        vol0=per_length*vol0/(2*pi*8.e-7)	!  2 for div(r) in 2D
      ENDIF
c
c						FFT of the magnetic field

      DO iline=1,n_fieldlines
c							get config space data
        DO ipoloidal=0,mmcell-1				! do not repeat point 0
          p0%point = contour(1:2,ipoloidal,iline)
          CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
          CALL refine_cell(p0,item,failure,x_start,y_start,.TRUE.)
          ibl=item%ib
          CALL lagr_quad_eval(rbc(ibl)%be_eq,x_start,y_start,0_i4)
          CALL lagr_quad_eval(rbc(ibl)%be,x_start,y_start,0_i4)
c
c						add the equilibrium fields
c						(I should trap for linear)
          rbc(ibl)%be%f(1:2,1) = rbc(ibl)%be%f(1:2,1)
     &                          +rbc(ibl)%be_eq%f(1:2)
          IF (geom=='tor') THEN
            rval = p0%point(1)
            rbc(ibl)%be%f(3,:)=rbc(ibl)%be%f(3,:)*rval
          ELSE
            rval = 1._r8
          ENDIF
          rbc(ibl)%be%f(3,1)=rbc(ibl)%be%f(3,1)+rbc(ibl)%be_eq%f(3)

          jacobian = 
     &     transform(2,2,ipoloidal,iline)*transform(1,1,ipoloidal,iline)
     &    -transform(1,2,ipoloidal,iline)*transform(2,1,ipoloidal,iline)
     &        +TINY(jacobian)
          bu_rho(:,ipoloidal,iline)=
     &           (transform(2,2,ipoloidal,iline)*rbc(ibl)%be%f(1,:)
     &           -transform(1,2,ipoloidal,iline)*rbc(ibl)%be%f(2,:))
     &           /jacobian
          bu_the(:,ipoloidal,iline)=
     &           (transform(1,1,ipoloidal,iline)*rbc(ibl)%be%f(2,:)
     &           -transform(2,1,ipoloidal,iline)*rbc(ibl)%be%f(1,:))
     &           /jacobian
          bl_rho(:,ipoloidal,iline)=
     &            transform(1,1,ipoloidal,iline)*rbc(ibl)%be%f(1,:)
     &           +transform(2,1,ipoloidal,iline)*rbc(ibl)%be%f(2,:)
          bl_the(:,ipoloidal,iline)=
     &            transform(1,2,ipoloidal,iline)*rbc(ibl)%be%f(1,:)
     &           +transform(2,2,ipoloidal,iline)*rbc(ibl)%be%f(2,:)
          bl_zed(:,ipoloidal,iline)=rbc(ibl)%be%f(3,:)
          bu_zed(:,ipoloidal,iline)=rbc(ibl)%be%f(3,:)/rval**2
        ENDDO

c
c							convert to fft
        DO imode = 1,num_fourier
c
c							REAL and IMAG bu_rho
          c(1:mmcell) = REAL(bu_rho(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)
          c(1:mmcell) = AIMAG(bu_rho(imode,:,iline))
          CALL rfftf       
          bu_rho(imode,:,iline) = 
     &               (d(1:mmcell) + (0,1)*c(1:mmcell))/REAL(nfft,r8)
c
c							REAL and IMAG bu_the
          c(1:mmcell) = REAL(bu_the(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)/REAL(nfft,r8)
          c(1:mmcell) = AIMAG(bu_the(imode,:,iline))
          CALL rfftf       
          bu_the(imode,:,iline) = 
     &                d(1:mmcell) + (0,1)*c(1:mmcell)/REAL(nfft,r8)
c
c							REAL and IMAG bl_rho
          c(1:mmcell) = REAL(bl_rho(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)
          c(1:mmcell) = AIMAG(bl_rho(imode,:,iline))
          CALL rfftf       
          bl_rho(imode,:,iline) =
     &               (d(1:mmcell) + (0,1)*c(1:mmcell))/REAL(nfft,r8)
c
c							REAL and IMAG bl_the
          c(1:mmcell) = REAL(bl_the(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)/REAL(nfft,r8)
          c(1:mmcell) = AIMAG(bl_the(imode,:,iline))
          CALL rfftf       
          bl_the(imode,:,iline) =
     &                d(1:mmcell) + (0,1)*c(1:mmcell)/REAL(nfft,r8)
c
c							REAL and IMAG bl_zed
          c(1:mmcell) = REAL(bl_zed(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)/REAL(nfft,r8)
          c(1:mmcell) = AIMAG(bl_zed(imode,:,iline))
          CALL rfftf       
          bl_zed(imode,:,iline) =
     &                d(1:mmcell) + (0,1)*c(1:mmcell)/REAL(nfft,r8)
c
c							REAL and IMAG bu_zed
          c(1:mmcell) = REAL(bu_zed(imode,:,iline),r8)
          CALL rfftf       
          d(:) = c(:)/REAL(nfft,r8)
          c(1:mmcell) = AIMAG(bu_zed(imode,:,iline))
          CALL rfftf       
          bu_zed(imode,:,iline) =
     &                d(1:mmcell) + (0,1)*c(1:mmcell)/REAL(nfft,r8)
        ENDDO
      ENDDO

c
c						Collect m-coefficients
c						for real&imag parts of
c						expansion for each n.
      bu_nr(1,:,0,:)=REAL(bu_rho(:,0,:))
      bl_nr(1,:,0,:)=REAL(bl_rho(:,0,:))
      bu_nr(2,:,0,:)=REAL(bu_the(:,0,:))
      bl_nr(2,:,0,:)=REAL(bl_the(:,0,:))
      bu_nr(3,:,0,:)=REAL(bu_zed(:,0,:))
      bl_nr(3,:,0,:)=REAL(bl_zed(:,0,:))
      bu_ni(1,:,0,:)=AIMAG(bu_rho(:,0,:))
      bl_ni(1,:,0,:)=AIMAG(bl_rho(:,0,:))
      bu_ni(2,:,0,:)=AIMAG(bu_the(:,0,:))
      bl_ni(2,:,0,:)=AIMAG(bl_the(:,0,:))
      bu_ni(3,:,0,:)=AIMAG(bu_zed(:,0,:))
      bl_ni(3,:,0,:)=AIMAG(bl_zed(:,0,:))
      DO ipoloidal=1,mm
        bu_nr(1,:,ipoloidal,:)=REAL(bu_rho(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bu_rho(:,2*ipoloidal  ,:),r8)
        bl_nr(1,:,ipoloidal,:)=REAL(bl_rho(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bl_rho(:,2*ipoloidal  ,:),r8)
        bu_nr(2,:,ipoloidal,:)=REAL(bu_the(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bu_the(:,2*ipoloidal  ,:),r8)
        bl_nr(2,:,ipoloidal,:)=REAL(bl_the(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bl_the(:,2*ipoloidal  ,:),r8)
        bu_nr(3,:,ipoloidal,:)=REAL(bu_zed(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bu_zed(:,2*ipoloidal  ,:),r8)
        bl_nr(3,:,ipoloidal,:)=REAL(bl_zed(:,2*ipoloidal-1,:),r8)
     $                  +(0,1)*REAL(bl_zed(:,2*ipoloidal  ,:),r8)
        bu_ni(1,:,ipoloidal,:)=AIMAG(bu_rho(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bu_rho(:,2*ipoloidal  ,:))
        bl_ni(1,:,ipoloidal,:)=AIMAG(bl_rho(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bl_rho(:,2*ipoloidal  ,:))
        bu_ni(2,:,ipoloidal,:)=AIMAG(bu_the(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bu_the(:,2*ipoloidal  ,:))
        bl_ni(2,:,ipoloidal,:)=AIMAG(bl_the(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bl_the(:,2*ipoloidal  ,:))
        bu_ni(3,:,ipoloidal,:)=AIMAG(bu_zed(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bu_zed(:,2*ipoloidal  ,:))
        bl_ni(3,:,ipoloidal,:)=AIMAG(bl_zed(:,2*ipoloidal-1,:))
     $                  +(0,1)*AIMAG(bl_zed(:,2*ipoloidal  ,:))
      ENDDO
c
c						Volume average energies
c						decomposed by harmonic
        WRITE(nim_wr,*)"Writing output file 'b2fft.dat'"
        OPEN(UNIT=tec2d,FILE='b2fft.dat',STATUS='UNKNOWN')
        WRITE(tec2d,*) 'VARIABLES = "m", "n" "BP2" '

        WRITE(zone_title,fmt='(x,"ZONE i=",i4," j=",i4," F=POINT")')
     &      mmcell/2+1,
     &      2*num_fourier-1
        WRITE(tec2d,*)zone_title

c						Separate opposite
c 						helicities.  Plot by
c						-nmax<=n<=nmax m>=0.
        emsum=0

        DO im = -num_fourier,num_fourier
          IF (im==0.OR.im==-1) CYCLE
          imode=ABS(im)
c
c						m=0 special case
          area_avg = 0.0
          vol_avg = 0.0
          DO iline=1,n_fieldlines
            vol_avg = vol_avg
     &               + (vol0(iline)-vol0(iline-1))*
     &                SUM(REAL(bl_nr(:,imode,0,iline),r8)*
     &                    REAL(bu_nr(:,imode,0,iline),r8))
            IF (imode/=1)
     &        vol_avg = vol_avg
     &               + (vol0(iline)-vol0(iline-1))*
     &                SUM(REAL(bl_ni(:,imode,0,iline),r8)*
     &                    REAL(bu_ni(:,imode,0,iline),r8))
          ENDDO
c-TMP
c         WRITE(tec2d,*)0_i4,(im/imode)*(imode-1),ABS(vol_avg)
          WRITE(tec2d,*)0_i4,(im/imode)*(imode-1),vol_avg
          emsum(imode)=emsum(imode)+vol_avg
c
c						m\=0 cases
          DO ipoloidal=1,mm
            vol_avg = 0.0
            DO iline=1,n_fieldlines
              IF (im>0) THEN
                vol_avg = vol_avg
     &                 + (vol0(iline)-vol0(iline-1))*2*REAL(
     &                  SUM(      (bl_nr(:,imode,ipoloidal,iline)
     &                      +(0,1)*bl_ni(:,imode,ipoloidal,iline))*
     &                       CONJG(bu_nr(:,imode,ipoloidal,iline)
     &                      +(0,1)*bu_ni(:,imode,ipoloidal,iline))))
              ELSE
                vol_avg = vol_avg
     &                 + (vol0(iline)-vol0(iline-1))*2*REAL(
     &                  SUM( CONJG(bl_nr(:,imode,ipoloidal,iline)
     &                      -(0,1)*bl_ni(:,imode,ipoloidal,iline))*
     &                            (bu_nr(:,imode,ipoloidal,iline)
     &                      -(0,1)*bu_ni(:,imode,ipoloidal,iline))))
              ENDIF
            ENDDO
c-TMP
            WRITE(tec2d,*)ipoloidal,(im/imode)*(imode-1),
     &                    vol_avg
c    &                    ABS(vol_avg)
            emsum(imode)=emsum(imode)+vol_avg
          ENDDO
        ENDDO
        CLOSE(UNIT=tec2d)
c
        OPEN (UNIT=tec2d,FILE='emsum',STATUS='UNKNOWN')
        DO imode=1,num_fourier
          WRITE(tec2d,*) imode-1,emsum(imode)
        ENDDO
        WRITE(tec2d,*) 'total e ',SUM(emsum)
        CLOSE(UNIT=tec2d)
c
c						Deallocate
      DEALLOCATE(bl_rho,bu_rho,bl_the,bu_the,bu_zed,bl_zed)
      DEALLOCATE(bu_nr,bu_ni,bl_nr,bl_ni)
      DEALLOCATE(transform,contour,theta0,rho0,vol0,q0,emsum)
      CALL bicube_dealloc(r2g)
        RETURN
c     DEALLOCATE(D,WA,CH,C)			! Why is this the problem
 
      END SUBROUTINE polfft


      SUBROUTINE psigrid
c
c						Calculation of q value
c						at the given surface
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE start_positions
      USE node_type_mod
      USE polfft_mod
      IMPLICIT NONE

      LOGICAL :: local_debug=.FALSE.
      TYPE(location_type) :: p0,p1
      TYPE(node_type), POINTER :: contour_start,contour_cnt
      TYPE(node_type), POINTER :: contour_cnt0,contour_cnt1
      REAL(r8) :: qvalue,area,volume

      REAL(r8) :: r_min,z_min,r_max,z_max
      TYPE(cell_type), POINTER :: item,item_pre
      INTEGER(i4) :: ibl,maxstep=5000,iline,ipoloidal,mmcell
      REAL(r8) :: x_start,y_start
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-12
      LOGICAL :: failure
      REAL(r8) :: theta_test,theta_value,gcosv,gsinv,rval
       


c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep0
      INTEGER(i4), PARAMETER :: neq=5,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8) :: tstart,tend,tprior,norm_value
      EXTERNAL derivative1
c
c						Allocate storage space 
      mmcell = 2*mm+1
      ALLOCATE(theta0(0:mmcell))
      ALLOCATE(rho0(0:n_fieldlines))
      ALLOCATE(vol0(0:n_fieldlines))
      ALLOCATE(q0(0:n_fieldlines))
      ALLOCATE(contour(2,0:mmcell,0:n_fieldlines))
      theta0=(/(ipoloidal,ipoloidal=0,mmcell)/)*twopi/REAL(mmcell,r8)
      contour(1,0:mmcell,0)=rmaxis
      contour(2,0:mmcell,0)=zmaxis
      rho0 = 0
      vol0 = 0
      q0 = 0
    
c-----------------------------------------------------------------------
c     Loop through Starting Points.
c-----------------------------------------------------------------------
      DO iline = 1,n_fieldlines
c
c                       	                set the start positions
        ibl =  ibl_fieldlines(iline)
        x_start =  x_fieldlines(iline)
        y_start = y_fieldlines(iline)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        p0%point(1:2)=rbc(ibl)%rz%f             ! (R,Z) coordinate
        contour(1:2,0,iline)=p0%point(1:2)
        contour(1:2,mmcell,iline)=p0%point(1:2)
c
c                               	        Find qvalue
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        rho0(iline) = SQRT(area)/pi
        q0(iline) = qvalue
        WRITE(nim_wr,*)iline,area,qvalue
c
c                                               Compute (R,Z) positions
c                                               necessary for poloidal
c                                               decomposition of surface.

c
c					Initialize max and min on each contour
        r_min = HUGE(1.0)
        z_min = HUGE(1.0)
        r_max = -HUGE(1.0)
        z_max = -HUGE(1.0)
        p1%point = p0%point


c
c                                               set up integrator parameters.
        itask=5 			!  integration with intermediate output
        iopt=1 			!  optional inputs are used
        mf=10			!  non-stiff Adams method of integration
        itol=1			!  absolute tolerance is just a scalar
        rtol=ind_tol0		!  relative tolerance
        atol=ind_tol0*ind_tol0 	!  absolute tolerance
        iwork=0
        rwork=0
        rwork(11)=1
        rwork(6) = 0.01 	! set maximum lsode-internal step size
        rwork(7) = rtol	 	! set minimum lsode-internal step size
        iwork(6)=2500   	! set maximum lsode-internal steps
        iwork(7)=2		! set maximum lsode-internal error mesg printed
c-----------------------------------------------------------------------
c     					set the integrator parameters.
        istate=1			!  indicates the first lsode call
        istep0=0
ctom       ibl = item%ib
        ibl =  ibl_fieldlines(iline)
        y_lsode(1:2)=p0%point(1:2)		! (R,Z) coordinate
        y_lsode(3)=0.				! PEST angle
        y_lsode(4)=0.				! Area integration
        y_lsode(5)=0.				! Volume integration
        y_lsode(neq+1)=x_start 			! x coordinate in logical space
        y_lsode(neq+2)=y_start 			! y coordinate in logical space
        y_lsode(neq+3)=ibl 			! block label
        tstart=ATAN2(y_lsode(2)-zmaxis,y_lsode(1)-rmaxis)
        tend=tstart+twopi
        rwork(1)=tend            ! set maximum angle of integration


        NULLIFY(contour_start)
        inner : DO
             IF (.NOT. ASSOCIATED (contour_start) ) THEN
               ALLOCATE(contour_start)
               contour_cnt => contour_start
             ELSE
               ALLOCATE(contour_cnt%next)
               contour_cnt => contour_cnt%next
             ENDIF
c
c                                               Regular starting spots.
c            contour_cnt%area = y_lsode(3)/qvalue
             contour_cnt%area = y_lsode(3)
             contour_cnt%rz  = y_lsode(1:2)

           tprior=tstart
           CALL dlsode(derivative1,neq,y_lsode,tstart,tend,
     &                 itol,rtol,atol,itask,
     &                 istate,iopt,rwork,lrw,iwork,liw,jac,mf)
           istep0=istep0+1
           IF(y_lsode(1) < r_min)r_min = y_lsode(1)
           IF(y_lsode(1) > r_max)r_max = y_lsode(1)
           IF(y_lsode(2) < z_min)z_min = y_lsode(2)
           IF(y_lsode(2) > z_max)z_max = y_lsode(2)
           IF(istate < 0)THEN
             qvalue = HUGE(1.0)
             area = -HUGE(1.0)
             volume = -HUGE(1.0)
             STOP 100
           ENDIF
           IF(ABS(tstart-tend) < atol)THEN
              ALLOCATE(contour_cnt%next)
              contour_cnt => contour_cnt%next
              contour_cnt%area = y_lsode(3)
              norm_value = twopi/contour_cnt%area
              contour_cnt%rz  = y_lsode(1:2)
              NULLIFY(contour_cnt%next)
              EXIT inner
           ENDIF
           IF(istep0 > maxstep)THEN
              WRITE(6,*)"polfft:istep0 > maxstep "
              STOP 200
           ENDIF
         ENDDO inner
c
c						Fix angle normalization
         contour_cnt => contour_start
         DO
           IF(.NOT.ASSOCIATED(contour_cnt))EXIT
           contour_cnt%area=contour_cnt%area*norm_value
           contour_cnt => contour_cnt%next
         ENDDO
c
c						Write each contour on debug
         IF(local_debug)THEN
           WRITE(7,*)"ZONE"
           contour_cnt => contour_start
           DO
             IF(.NOT.ASSOCIATED(contour_cnt))EXIT
             WRITE(7,*)contour_cnt%area,contour_cnt%rz
             contour_cnt => contour_cnt%next
           ENDDO
         ENDIF
c
c						Search for the desired output
c						values for the poloidal mesh.
       DO ipoloidal=1,mmcell-1
         theta_value = theta0(ipoloidal)
         contour_cnt => contour_start
         DO
             IF(.NOT.ASSOCIATED(contour_cnt%next))EXIT
             contour_cnt0 => contour_cnt
             contour_cnt1 => contour_cnt%next
             theta_test = (contour_cnt0%area-theta_value)*
     &                    (contour_cnt1%area-theta_value)
             IF(theta_test < 0)THEN
                contour(1:2,ipoloidal,iline)= contour_cnt0%rz
     &                  +(theta_value      -contour_cnt0%area)
     &                  /(contour_cnt1%area-contour_cnt0%area)
     &                  *(contour_cnt1%rz-contour_cnt0%rz)
             ENDIF
             contour_cnt => contour_cnt%next
         ENDDO
       ENDDO
c
c						Deallocate the contour markers
         contour_cnt => contour_start
         DO
           IF(.NOT.ASSOCIATED(contour_cnt))EXIT
           contour_start => contour_cnt
           contour_cnt => contour_cnt%next
           DEALLOCATE(contour_start)
         ENDDO
c
c						Am I back at the same point?
c						Am I not at the boundary?
         IF((ABS(tstart-tend) < atol).AND.(y_lsode(neq+3) > 0))THEN
           qvalue = y_lsode(3)/twopi
           area = y_lsode(4)
           volume = ABS(y_lsode(5))
         ELSE
           qvalue = -HUGE(1.0)
           area = -HUGE(1.0)
           volume = -HUGE(1.0)
         ENDIF
         vol0(iline)=volume
      ENDDO
c----------------------------------------------------------------------
c
c				Compute the transformation matrix for 
c				conversion to the (rho,THETA) vector directions.
c
      CALL bicube_alloc(r2g,mmcell,n_fieldlines,2_i4)
      r2g%xs=theta0
      r2g%ys=rho0
      r2g%fs(1:2,:,0)=0.
      DO iline=1,n_fieldlines
      theta_test = -twopi
      DO ipoloidal=0,mmcell
c       r2g%fs(1,ipoloidal,iline)=
c    &           (contour(1,ipoloidal,iline)-rmaxis)**2
c    &          +(contour(2,ipoloidal,iline)-zmaxis)**2
c       theta_value = ATAN2( contour(2,ipoloidal,iline)-zmaxis
c    &                      ,contour(1,ipoloidal,iline)-rmaxis)
c       IF(theta_value < theta_test) theta_value = theta_value + twopi
c       theta_test = theta_value
c       r2g%fs(2,ipoloidal,iline)=theta_value-theta0(ipoloidal)
        r2g%fs(1:2,ipoloidal,iline)=contour(1:2,ipoloidal,iline)
      ENDDO
      ENDDO

      CALL bicube_fit(r2g,"periodic","extrap")

      ALLOCATE(transform(2,2,0:mmcell,0:n_fieldlines))
      DO iline=1,n_fieldlines
      DO ipoloidal=0,mmcell
        CALL bicube_eval(r2g,theta0(ipoloidal),rho0(iline),1_i4)
c       rval = SQRT(r2g%fs(1,ipoloidal,iline))
c       gcosv = COS(theta0(ipoloidal)+r2g%f(2))
c       gsinv = SIN(theta0(ipoloidal)+r2g%f(2))
c       transform(1,1,ipoloidal,iline) = 
c    &            0.5/rval*gcosv*r2g%fy(1) - rval*gsinv*r2g%fy(2)
c       transform(1,2,ipoloidal,iline) = 
c    &            0.5/rval*gcosv*r2g%fx(1) - rval*gsinv*(1.0+r2g%fx(2))
c       transform(2,1,ipoloidal,iline) = 
c    &            0.5/rval*gsinv*r2g%fy(1) + rval*gcosv*r2g%fy(2)
c       transform(2,2,ipoloidal,iline) = 
c    &            0.5/rval*gsinv*r2g%fx(1) + rval*gcosv*(1.0+r2g%fx(2))
        transform(1:2,1,ipoloidal,iline)=r2g%fy(1:2)
        transform(1:2,2,ipoloidal,iline)=r2g%fx(1:2)
      ENDDO
      ENDDO
c					 At magnetic axis use iline=1 transform
c					This may or may not be OK?
      transform(:,:,:,0) = transform(:,:,:,1)	

c----------------------------------------------------------------------
c
c						Xdraw output files
      WRITE(nim_wr,*)"Writing output file 'polgrid.bin'"
      CALL open_bin(xdr_unit,"polgrid.bin","UNKNOWN","REWIND",32_i4)
      DO iline=0,n_fieldlines
      DO ipoloidal=0,mmcell
        WRITE(xdr_unit)(/REAL(contour(1:2,ipoloidal,iline),4)
     &                 ,REAL(rho0(iline),4)
     &                 ,REAL(theta0(ipoloidal),4)/)
      ENDDO
      WRITE(xdr_unit)
      ENDDO
      CLOSE(UNIT=xdr_unit)
c
c						Tecplot output files
      IF(out_tecplot)THEN
        WRITE(nim_wr,*)"Writing output file 'polgrid.dat'"
        OPEN(UNIT=tec2d,FILE='polgrid.dat',STATUS='UNKNOWN')
        WRITE(tec2d,*) 'VARIABLES = "R", "Z", "rho", "tau", "r2", "eta"'
     &                ,'"a11","a21","a12","a22"'
        WRITE(tec2d,*) 'ZONE '
     &                ,',i=',mmcell+1
     &                ,' j=',n_fieldlines+1
     &                ,',F=POINT'
        DO iline=0,n_fieldlines
        DO ipoloidal=0,mmcell
          WRITE(tec2d,*)contour(1:2,ipoloidal,iline)
     &                 ,rho0(iline)
     &                 ,theta0(ipoloidal)
     &                 ,r2g%fs(1,ipoloidal,iline)
     &                 ,r2g%fs(2,ipoloidal,iline)
     &                 ,transform(1:2,1:2,ipoloidal,iline)
        ENDDO
        ENDDO
        CLOSE(UNIT=tec2d)
      ENDIF

      END SUBROUTINE psigrid

      SUBROUTINE derivative1(neq,ind,y_lsode,dy)
      USE local
      USE magnetic_axis

      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), INTENT(IN) :: ind
      REAL(r8), INTENT(INOUT), DIMENSION(neq+3) :: y_lsode
      REAL(r8), INTENT(OUT), DIMENSION(neq) :: dy

      REAL(r8), DIMENSION(3) :: b_xyz
      REAL(r8) :: bmag,bigr,btheta,rho

      CALL get_bfield0(neq,y_lsode,b_xyz,bigr)
      bmag=SQRT(SUM(b_xyz**2))+TINY(bmag)

      rho = SQRT((y_lsode(1)-rmaxis)**2+(y_lsode(2)-zmaxis)**2)
      btheta = COS(ind)*b_xyz(2) - SIN(ind)*b_xyz(1)
      dy(1:3) = b_xyz*rho/(btheta+TINY(btheta))
      dy(3) = rho/(bigr**2*(btheta+TINY(btheta)))
      dy(4) = dy(1)*y_lsode(2)+2._r8*dy(2)*y_lsode(1)
      dy(5) = bigr*(y_lsode(1)*dy(2)-y_lsode(2)*dy(1))

      RETURN
      END SUBROUTINE derivative1
