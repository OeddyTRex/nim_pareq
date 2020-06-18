      SUBROUTINE psigrid0
c
c				Integration of n=0 NIMROD fields for
c				poloidal flux and other diagnostics.
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
      IMPLICIT NONE

      LOGICAL :: local_debug=.FALSE.
      TYPE(location_type) :: p0
c     TYPE(cell_type), POINTER :: item
      INTEGER(i4) :: ibl,iline=0,maxline=1000,maxstep=5000
      REAL(r8) :: x_start,y_start,rmax
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-12
      REAL(r8) :: qvalue,area
      REAL(r8) :: r_min,z_min,r_max,z_max



      WRITE(nim_wr,*)"Compute poloidal decomposition grid."
c-----------------------------------------------------------------------
c     evaluate rmax.
c-----------------------------------------------------------------------
      rmax=0
      DO ibl=1,nrblc
        rmax=MAX(rmax,MAXVAL(rbc(ibl)%rz%fs(:,:,1)))
      ENDDO
c-----------------------------------------------------------------------
c     Loop through Starting Points.
c-----------------------------------------------------------------------
      DO iline = 1,n_fieldlines
c
c     					set the start positions
        ibl =  ibl_fieldlines(iline)
        x_start =  x_fieldlines(iline)
        y_start = y_fieldlines(iline)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        p0%point(1:2)=rbc(ibl)%rz%f		! (R,Z) coordinate
c
c     					Find qvalue
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        WRITE(nim_wr,*)iline,area,qvalue
c
c						Compute (R,Z) positions
c						necessary for poloidal
c						decomposition of surface.
      ENDDO
      END SUBROUTINE psigrid0

      SUBROUTINE get_bfield0(neq,y_lsode,b_xyz,bigr)
c
c						Intended to compute the
c						magnetic field strength
c						and presumes that x and y
c						are in the block with label ib.
c
      USE input0
      USE global
      USE dumpc
      USE cell_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), DIMENSION(neq+3), INTENT(INOUT) :: y_lsode 
      REAL(r8), DIMENSION(3), INTENT(OUT) :: b_xyz  
      REAL(r8), INTENT(OUT) :: bigr

      REAL(r8) :: jacobian,norm_factor
      REAL(r8) :: x,y,r,z,gcosv,gsinv,angle
      REAL(r8) :: btol = 1.e-12
      INTEGER(i4) :: ib,kfactor,imodes,isearch,jsearch
      TYPE(location_type):: p0
      TYPE(cell_type), POINTER :: item_pre
      LOGICAL :: failure


c
c						Translate y_lsode to local
c						representaion
      r=y_lsode(1)
      z=y_lsode(2)
      x=y_lsode(neq+1)
      y=y_lsode(neq+2)
      ib=NINT(y_lsode(neq+3))
c
c						Update ib,x,y for r,z
      p0%point(1)=r
      p0%point(2)=z
      isearch = NINT(search_idim*
     &              (r-search_rmin)/(search_rmax-search_rmin))
      jsearch = NINT(search_jdim*
     &              (z-search_zmin)/(search_zmax-search_zmin))
      lucky => search_map(isearch,jsearch)%next
      CALL rz_to_cell(p0,lucky,item_pre,failure,.FALSE.)
      CALL refine_cell(p0,lucky,failure,x,y,.FALSE.)
      ib=lucky%ib
      IF(failure)THEN
        WRITE(out_unit,*)" # FAILURE at ",p0%point(1:2)
        CALL nim_stop("ERROR in psigrid")
      ENDIF
      y_lsode(neq+1)=x
      y_lsode(neq+2)=y
      y_lsode(neq+3)=ib

      IF(ib  ==  0)THEN
c
c						no change, try to keep lsode
c						looping at some point.
        b_xyz=0
        IF (geom=='tor') THEN
          bigr=r
        ELSE
          bigr=1
        ENDIF
        RETURN
      ELSE IF(ib <= nrblc)THEN
c
c						Calculate the equilibrium
c						magnetic field from rblock
        call lagr_quad_eval(rbc(ib)%be_eq,x,y,0_i4)
        b_xyz=rbc(ib)%be_eq%f
c
c						Convert toroidal component
c						to cylindrical and set
c						angle.
!       angle=y_lsode(3)
        angle=0._r8
        IF (geom=='tor') THEN
          b_xyz(3)=b_xyz(3)/r
          bigr=r
        ELSE
          angle=twopi*angle/per_length
          bigr=1
        ENDIF
c
c
c						Sum the perturbed magnetic
c						fields over the number of
c						Fourier modes.
c	be%f(1,1)-------r-comp, 1st mode (n=0)
c	be%f(2,1)-------z-comp, 1st mode (n=0)
c	be%f(3,1)-------R*phi-comp, 1st mode
c	be%f(1,2)-------r-comp, 2nd mode (n=1)
c	be%f(2,2)-------z-comp, 2nd mode
c	be%f(3,2)-------R*phi-comp, 2nd mode
c	be%f(1,3)-------r-comp, 3rd mode (n=2)
c	
c
c       B_tilde(R,Z,phi)=sum_over_modes[B_tilde_n(R,Z)*exp( i*n*phi)
c                                +conjg(B_tilde_n(R,Z)*exp(-i*n*phi)]
c
        call lagr_quad_eval(rbc(ib)%be,x,y,0_i4)
c
        IF(nonlinear)THEN
c
c						n=0 Fourier component.
c
          b_xyz(:)=b_xyz(:)+rbc(ib)%be%f(:,1)
c
c						n/=0 Fourier components.
c
c         DO imodes=2,nmodes
c           gcosv=COS(angle*(imodes-1))
c           gsinv=SIN(angle*(imodes-1))
c           b_xyz(:)=b_xyz(:)+
c    &          2*(REAL(rbc(ib)%be%f(:,imodes),r8)*gcosv-
c    &             AIMAG(rbc(ib)%be%f(:,imodes)  )*gsinv)
c         ENDDO
        ELSE
c
c						Presumes lin_nmodes=1 
c         
c         gcosv=COS(angle*lin_nmax)
c         gsinv=SIN(angle*lin_nmax)
c         norm_factor=2._r8
c         IF(lin_nmax == 0)norm_factor=1._r8
c         b_xyz(:)=b_xyz(:)+norm_factor
c    &           *(REAL(rbc(ib)%be%f(:,1),r8)*gcosv-
c    &             AIMAG(rbc(ib)%be%f(:,1)  )*gsinv)
        ENDIF
      ELSE
c
c						Calculate the equilibrium
c						magnetic field from tblock
         CALL tri_linear_eval(tbc(lucky%ib)%be_eq,
     &                        tbc(lucky%ib)%tgeom,
     &                   p0%point(1),p0%point(2),lucky%icell,0_i4)
        b_xyz=tbc(ib)%be_eq%f
c
c						Convert toroidal component
c						to cylindrical and set angle.
        angle=y_lsode(3)
        IF (geom=='tor') THEN
          b_xyz(3)=b_xyz(3)/r
          bigr=r
        ELSE
          angle=twopi*angle/per_length
          bigr=1
        ENDIF
c
c						Calculate the equilibrium
c						magnetic field from tblock
        CALL tri_linear_eval(tbc(lucky%ib)%be,
     &                       tbc(lucky%ib)%tgeom,
     &                   p0%point(1),p0%point(2),lucky%icell,0_i4)
c
c
        IF(nonlinear)THEN
c
c						n=0 Fourier component.
c
          b_xyz(:)=b_xyz(:)+rbc(ib)%be%f(:,1)
c
c						n/=0 Fourier components.
c
c         DO imodes=2,nmodes
c           gcosv=COS(angle*(imodes-1))
c           gsinv=SIN(angle*(imodes-1))
c           b_xyz(:)=b_xyz(:)+
c    &          2*(REAL(tbc(ib)%be%f(:,imodes),r8)*gcosv-
c    &             AIMAG(tbc(ib)%be%f(:,imodes)  )*gsinv)
c         ENDDO
        ELSE
c
c						Presumes lin_nmodes=1 
c         
c         gcosv=COS(angle*lin_nmax)
c         gsinv=SIN(angle*lin_nmax)
c         norm_factor=2._r8
c         IF(lin_nmax == 0)norm_factor=1._r8
c         kfactor=0
c         b_xyz(:)=b_xyz(:)+norm_factor
c    &           *(REAL(tbc(ib)%be%f(:,1),r8)*gcosv-
c    &             AIMAG(tbc(ib)%be%f(:,1)  )*gsinv)
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE get_bfield0


      SUBROUTINE derivative(neq,ind,y_lsode,dy)
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
      dy(3) = dy(3)/bigr
      dy(4) = dy(1)*y_lsode(2)+2._r8*dy(2)*y_lsode(1)

      RETURN
      END SUBROUTINE derivative

      SUBROUTINE qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
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
      IMPLICIT NONE

      TYPE(location_type),INTENT(IN) :: p0
      REAL(r8), INTENT(OUT) :: qvalue,area
      REAL(r8),INTENT(OUT) :: r_min,z_min,r_max,z_max

      LOGICAL :: local_debug=.FALSE.
      TYPE(cell_type), POINTER :: item,item_pre
      INTEGER(i4) :: ibl,maxstep=50000
      REAL(r8) :: x_start,y_start
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-12
      LOGICAL :: failure
      TYPE(location_type) :: p1
       


c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep0
      INTEGER(i4), PARAMETER :: neq=4,liw=20,lrw=22+neq*16
c     INTEGER(i4), PARAMETER :: neq=4,liw=20+neq,lrw=22+neq*(9+neq)
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8) :: tstart,tend,tprior
      EXTERNAL derivative
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
c     mf=22			!  stiff integration
      itol=1			!  absolute tolerance is just a scalar
      rtol=ind_tol0		!  relative tolerance
      atol=ind_tol0*ind_tol0 	!  absolute tolerance
      iwork=0
      rwork=0
      rwork(11)=1
      rwork(6) = 0.01 		! set maximum lsode-internal step size
      rwork(7) = rtol	 	! set minimum lsode-internal step size
c-TMP
c     iwork(6)=2500   		! set maximum lsode-internal steps
      iwork(6)=maxstep   		! set maximum lsode-internal steps
      iwork(7)=2		! set maximum lsode-internal error mesg printed
c-----------------------------------------------------------------------
c     					set the integrator parameters.
        istate=1			!  indicates the first lsode call
        istep0=0
        item => start
        item_pre => start
        CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
        CALL refine_cell(p0,item,failure,x_start,y_start,.FALSE.)
        ibl = item%ib
        y_lsode(1:2)=p0%point(1:2)		! (R,Z) coordinate
        y_lsode(3)=0.				! PEST angle
        y_lsode(4)=0.				! Area integration
        y_lsode(neq+1)=x_start 			! x coordinate in logical space
        y_lsode(neq+2)=y_start 			! y coordinate in logical space
        y_lsode(neq+3)=ibl 			! block label
        tstart=ATAN2(y_lsode(2)-zmaxis,y_lsode(1)-rmaxis)
        tend=tstart+twopi
        rwork(1)=tend            ! set maximum angle of integration


        IF(local_debug)WRITE(7,*)"ZONE"
        inner : DO
           IF(local_debug)WRITE(7,*)tstart,y_lsode(1:2)
           tprior=tstart
           CALL dlsode(derivative,neq,y_lsode,tstart,tend,
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
             RETURN
c            STOP 100
           ENDIF
           IF(ABS(tstart-tend) < atol)EXIT inner
           IF(istep0 > maxstep)THEN
              WRITE(6,*)"istep0 > maxstep on LHS"
              STOP 200
           ENDIF
c          IF(ABS(rwork(11)-rtol) < atol )THEN
c            qvalue = 0.0
c            area = 0.0
c            RETURN
c          ENDIF
         ENDDO inner
c
c					Am I back at the same point?
c					Am I not at the boundary?
         IF((ABS(tstart-tend) < atol).AND.(y_lsode(neq+3) > 0))THEN
           qvalue = y_lsode(3)/twopi
           area = y_lsode(4)
         ELSE
           qvalue = -HUGE(1.0)
           area = -HUGE(1.0)
         ENDIF
      RETURN

      END SUBROUTINE qcompute


      SUBROUTINE q_to_start(found)
c
c					Find starting position based on q.
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE node_type_mod
      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: found

      LOGICAL :: local_debug=.FALSE.
      TYPE(cell_type), POINTER :: item
      INTEGER(i4) :: ibl,iline=0,maxline=1000,maxstep=5000,nline
      INTEGER(i4) :: mvalue,nvalue,i_bisect,max_bisect = 100,i_repeat
      REAL(r8) :: x_start,y_start,rmax
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-12
      REAL(r8) :: qvalue,area,qtest
      REAL(r8) :: r_min,z_min,r_max,z_max
      REAL(r8) :: x_left,x_right,q_left,q_right
      REAL(r8) :: qmin,qmax,area_last = 0
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: qprof
      TYPE(location_type) :: p0
      REAL(r8) :: q_search_value

      found = .FALSE.
      NULLIFY(surface_first)
      qmin=HUGE(1.0)
      qmax=-HUGE(1.0)

c
c					Number of surfaces.
      item => start
      DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO
      nline=iline
      ALLOCATE(qprof(nline,2))

c
c     					Set the coarse mesh
      iline=0
      item => start
      outer : DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        IF(local_debug)WRITE(7,*)"# Loop ",iline
        IF(local_debug)WRITE(8,*)"# Loop ",iline
c
c     					set the start positions
        ibl = item%ib
        x_start = item%p(1,1)
        y_start = item%p(2,1)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        p0%point(1:2)=rbc(ibl)%rz%f		! (R,Z) coordinate
c
c     					Find qvalue
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        WRITE(nim_wr,*)iline,area,qvalue
        IF(area < area_last)THEN		! How can this be?
          nline = iline - 1
          IF(nline < 1)THEN
            iline=iline-1
            item => item%face(1)%p
            CYCLE
          ENDIF
          EXIT outer
        ENDIF
        IF(area < 0) THEN		! Trap for hitting the wall.
          nline=iline-1
          IF(nline < 1)THEN
            iline=iline-1
            item => item%face(1)%p
            CYCLE
          ENDIF
          WRITE(nim_wr,*)"This n=0 fieldline intersects the wall."
          WRITE(nim_wr,*)"The calculation continues short of the wall."
          EXIT outer
        ENDIF
        area_last = area
        qprof(iline,1)=area
        qprof(iline,2)=qvalue
        IF(qvalue > qmax)qmax=qvalue
        IF(qvalue < qmin)qmin=qvalue
        

c
c						Advance to new start position
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO outer
      WRITE(nim_wr,*)'Range of q values',qmin,' to ',qmax
c-----------------------------------------------------------------
c						Loop through qvalues
      DO nvalue = 1,nn
       mloop: DO mvalue = FLOOR(qmin)*nvalue,CEILING(qmax)*nvalue
          q_search_value = REAL(mvalue,r8)/REAL(nvalue,r8)
          DO i_repeat=1,nvalue-1
            IF(MOD(q_search_value,REAL(i_repeat,r8)) == 0)CYCLE mloop
          ENDDO
c-----------------------------------------------------------------
c
c						Use the coarse grid to
c						bracket a local search.
      iline=0
      item => start
      DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        IF(iline+1 == nline) EXIT

        qtest =(qprof(iline,2)-q_search_value)
     &         *(qprof(iline+1,2)-q_search_value)
c
c						Apply bisection rule to refine
c						for the starting positions.
        IF(qtest < 0 ) THEN
          IF(local_debug)THEN
            WRITE(nim_wr,*)iline,item%id
            WRITE(nim_wr,*)qprof(iline,2),qprof(iline+1,2)
          ENDIF
          ibl = item%ib
          x_left  = item%p(1,1)
          x_right = item%p(1,1)+1.
          q_left  = qprof(iline,2)
          q_right = qprof(iline+1,2)
          x_start = item%p(1,1)+0.5
          y_start = item%p(2,1)
          i_bisect = 0
          DO
            i_bisect = i_bisect+1
            CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
            p0%point(1:2)=rbc(ibl)%rz%f
            CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
            IF(local_debug)WRITE(nim_wr,*)area,qvalue
            IF(qvalue > q_search_value)THEN
              q_right=qvalue 
              x_right=x_start
            ELSE
              q_left=qvalue 
              x_left=x_start
            ENDIF
c
c						Check convergence, add to list.
            IF(((qvalue - q_search_value)<ind_tol0).OR.
     &         (i_bisect > max_bisect))THEN
             found = .TRUE.
             IF (.NOT. ASSOCIATED (surface_first) ) THEN
               ALLOCATE(surface_first)
               surface_cnt => surface_first
             ELSE
               ALLOCATE(surface_cnt%next)
               surface_cnt => surface_cnt%next
             ENDIF
c
c						Regular starting spots.
             surface_cnt%area = area
             surface_cnt%ibl  = item%ib
             surface_cnt%x_start  = x_start
             surface_cnt%y_start  = y_start
             surface_cnt%rz  = p0%point(1:2)
c
c						Revised starting points.
             CALL bncompute(p0,area,qvalue,r_min,z_min,r_max,z_max)

             surface_cnt%area = node_bmin%area
             surface_cnt%ibl = node_bmin%ibl
             surface_cnt%x_start = node_bmin%x_start
             surface_cnt%y_start = node_bmin%y_start
             surface_cnt%rz = node_bmin%rz
             ALLOCATE(surface_cnt%next)
             surface_cnt => surface_cnt%next

             surface_cnt%area = node_bzero%area
             surface_cnt%ibl = node_bzero%ibl
             surface_cnt%x_start = node_bzero%x_start
             surface_cnt%y_start = node_bzero%y_start
             surface_cnt%rz = node_bzero%rz
             ALLOCATE(surface_cnt%next)
             surface_cnt => surface_cnt%next

             surface_cnt%area = node_bmax%area
             surface_cnt%ibl = node_bmax%ibl
             surface_cnt%x_start = node_bmax%x_start
             surface_cnt%y_start = node_bmax%y_start
             surface_cnt%rz = node_bmax%rz

             DEALLOCATE(node_bmin,node_bmax,node_bzero)

             NULLIFY (surface_cnt%next)
c            WRITE(nim_wr,*)area,qvalue
             EXIT
            ENDIF
            x_start=0.5*(x_left+x_right)
          ENDDO
        ENDIF
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO
      ENDDO mloop
      ENDDO
      RETURN

      END SUBROUTINE q_to_start

      SUBROUTINE bncompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
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
      USE node_type_mod
      IMPLICIT NONE

      TYPE(location_type),INTENT(IN) :: p0
      REAL(r8), INTENT(OUT) :: qvalue,area
      REAL(r8),INTENT(OUT) :: r_min,z_min,r_max,z_max

      LOGICAL :: local_debug=.FALSE.
      TYPE(cell_type), POINTER :: item,item_pre
      INTEGER(i4) :: ibl,maxstep=5000
      REAL(r8) :: x_start,y_start
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-12
      LOGICAL :: failure
      REAL(r8), DIMENSION(3) :: b0_xyz
      REAL(r8), DIMENSION(3) :: b_xyz
      REAL(r8), DIMENSION(2) :: gradpsi_hat
      REAL(r8) :: bigr,bmag,bnormal,angle
      REAL(r8) :: bnormal_min ,bnormal_max,bnormal_zero
       


c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep0
      INTEGER(i4), PARAMETER :: neq=4,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8) :: tstart,tend,tprior
      EXTERNAL derivative

      IF (.NOT. ASSOCIATED (node_bmin) ) ALLOCATE(node_bmin)
      IF (.NOT. ASSOCIATED (node_bmax) ) ALLOCATE(node_bmax)
      IF (.NOT. ASSOCIATED (node_bzero) ) ALLOCATE(node_bzero)
      bnormal_min=HUGE(1.0)
      bnormal_max=-HUGE(1.0)
      bnormal_zero=HUGE(1.0)
c
c					Initialize max and min on each contour
      r_min = HUGE(1.0)
      z_min = HUGE(1.0)
      r_max = -HUGE(1.0)
      z_max = -HUGE(1.0)


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
      rwork(6) = 0.01 		! set maximum lsode-internal step size
      rwork(7) = rtol	 	! set minimum lsode-internal step size
      iwork(6)=2500   		! set maximum lsode-internal steps
      iwork(7)=2		! set maximum lsode-internal error mesg printed
c-----------------------------------------------------------------------
c     					set the integrator parameters.
        istate=1			!  indicates the first lsode call
        istep0=0
        item => start
        item_pre => start
        CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
        CALL refine_cell(p0,item,failure,x_start,y_start,.FALSE.)
        ibl = item%ib
        y_lsode(1:2)=p0%point(1:2)		! (R,Z) coordinate
        y_lsode(3)=0.				! PEST angle
        y_lsode(4)=0.				! Area integration
        y_lsode(neq+1)=x_start 			! x coordinate in logical space
        y_lsode(neq+2)=y_start 			! y coordinate in logical space
        y_lsode(neq+3)=ibl 			! block label
        tstart=ATAN2(y_lsode(2)-zmaxis,y_lsode(1)-rmaxis)
        tend=tstart+twopi
        rwork(1)=tend            ! set maximum angle of integration


        inner : DO
           IF(local_debug)WRITE(7,*)y_lsode(1:2)
           IF(local_debug)WRITE(8,*)tstart,y_lsode(1:2)
           tprior=tstart
           CALL dlsode(derivative,neq,y_lsode,tstart,tend,
     &                 itol,rtol,atol,itask,
     &                 istate,iopt,rwork,lrw,iwork,liw,jac,mf)
c
c						Compute unit normal grad psi
           CALL get_bfield0(neq,y_lsode,b0_xyz,bigr)
           bmag = SQRT(b0_xyz(1)**2+b0_xyz(2)**2)+TINY(bmag)
           gradpsi_hat(1) =  b0_xyz(2)/bmag
           gradpsi_hat(2) = -b0_xyz(1)/bmag
c
c						Compute normal component of
c						total magnetic field.
c		Note that y_lsode(3) is the toroidal angle for get_bfield
           angle=y_lsode(3)
           y_lsode(3) = 0
           CALL get_bfield(neq,y_lsode,b_xyz,bigr)
           y_lsode(3) = angle
           bnormal = b_xyz(1)*gradpsi_hat(1)+b_xyz(2)*gradpsi_hat(2)
c
c						Trap for three start locactions
           IF(bnormal < bnormal_min)THEN
             node_bmin%area = 0.
             node_bmin%ibl  = NINT(y_lsode(neq+3))
             node_bmin%x_start  = y_lsode(neq+1)
             node_bmin%y_start  = y_lsode(neq+2)
             node_bmin%rz(1:2)  = y_lsode(1:2)
             bnormal_min = bnormal
           ENDIF
           IF(bnormal > bnormal_max)THEN
             node_bmax%area = 0.
             node_bmax%ibl  = NINT(y_lsode(neq+3))
             node_bmax%x_start  = y_lsode(neq+1)
             node_bmax%y_start  = y_lsode(neq+2)
             node_bmax%rz(1:2)  = y_lsode(1:2)
             bnormal_max = bnormal
           ENDIF
           IF(ABS(bnormal) < bnormal_zero)THEN
             node_bzero%area = 0.
             node_bzero%ibl  = NINT(y_lsode(neq+3))
             node_bzero%x_start  = y_lsode(neq+1)
             node_bzero%y_start  = y_lsode(neq+2)
             node_bzero%rz(1:2)  = y_lsode(1:2)
             bnormal_zero = ABS(bnormal)
           ENDIF
c          write(7,*)tstart,bnormal

c
c						Determine min and max radius.
           istep0=istep0+1
           IF(y_lsode(1) < r_min)r_min = y_lsode(1)
           IF(y_lsode(1) > r_max)r_max = y_lsode(1)
           IF(y_lsode(2) < z_min)z_min = y_lsode(2)
           IF(y_lsode(2) > z_max)z_max = y_lsode(2)
           IF(ABS(tstart-tend) < atol)EXIT inner
           IF(istate < 0)THEN
             qvalue = HUGE(1.0)
             area = -HUGE(1.0)
             EXIT INNER
c            STOP 100
           ENDIF
           IF(istep0 > maxstep)THEN
              WRITE(6,*)"istep0 > maxstep on LHS"
              STOP 200
           ENDIF
c          IF(ABS(rwork(11)-rtol) < atol )THEN
c            qvalue = 0.0
c            area = 0.0
c            RETURN
c          ENDIF
         ENDDO inner
         qvalue = y_lsode(3)/twopi
         area = y_lsode(4)
         node_bmin%area = area
         node_bmax%area = area
         node_bzero%area = area
      RETURN

      END SUBROUTINE bncompute
