      MODULE magnetic_axis

      USE local
      REAL(r8) :: rmaxis, zmaxis

      CONTAINS
      SUBROUTINE compute_magnetic_axis
      USE local
      USE io
      USE input
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      IMPLICIT NONE
 
      LOGICAL :: local_debug=.FALSE.
      TYPE(location_type) :: p0
      TYPE(cell_type), POINTER :: item,itemt,item_pre
      INTEGER(i4) :: ibl,iline=0,maxline=1000,maxstep=5000,nlines
      REAL(r8) :: rmax,x_start,y_start,bigr
      INTEGER(i4) :: i,i_plus_br,i_neg_br,i_plus_bz,i_neg_bz,i_outer
      REAL(r8), DIMENSION(3) :: b_xyz
      REAL(r8) :: rmaxis_new, zmaxis_new
      REAL(r8) :: rmaxis_new_min, zmaxis_new_min
      REAL(r8) :: rmaxis_new_max, zmaxis_new_max
      LOGICAL :: failure
c-----------------------------------------------------------------------
c     declarations pertaining to lsode.
c-----------------------------------------------------------------------
      REAL(r8) :: ind_tol0=1.0e-10
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,istep0
      INTEGER(i4), PARAMETER :: neq=4,liw=20,lrw=22+neq*16
      REAL(r8) :: atol,rtol
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(lrw) :: rwork=0
      REAL(r8) :: tstart,tend,tprior
      EXTERNAL derivative
c
c                                               set up integrator parameters.
      itask=5                   !  integration with intermediate output
      iopt=1                    !  optional inputs are used
      mf=10                     !  non-stiff Adams method of integration
      itol=1                    !  absolute tolerance is just a scalar
      rtol=ind_tol0             !  relative tolerance
      atol=ind_tol0*ind_tol0    !  absolute tolerance
      iwork=0
      rwork=0
      rwork(11)=1
      rwork(6) = 0.01           ! set maximum lsode-internal step size
      rwork(7) = rtol           ! set minimum lsode-internal step size
      iwork(6)=2500             ! set maximum lsode-internal steps
      iwork(7)=2                ! set maximum lsode-internal error mesg printed
c-----------------------------------------------------------------------
c     find rmax
c-----------------------------------------------------------------------
      WRITE(nim_wr,*)"Computing magnetic axis."
      rmax=0
      DO ibl=1,nrblc
        rmax=MAX(rmax,MAXVAL(rbc(ibl)%rz%fs(:,:,1)))
      ENDDO
c-----------------------------------------------------------------------
c     Search for the cell which likely contains the magnetic axis. (itemt)
c-----------------------------------------------------------------------
      itemt => start
      item => start
      nlines=0
      DO WHILE(ASSOCIATED(item%next))
        nlines=nlines+1
        ibl = item%ib
        i_plus_br = 0		! number of nodes with positive br
        i_neg_br = 0		! number of nodes with negative br
        i_plus_bz = 0		! number of nodes with positive bz
        i_neg_bz = 0		! number of nodes with negative bz
        DO i=1,item%nodes
          x_start = item%p(1,i)
          y_start = item%p(2,i)
          CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
          y_lsode(1:2)=rbc(ibl)%rz%f            ! (R,Z) coordinate
          y_lsode(3)=0.                         ! PEST angle
          y_lsode(4)=0.                         ! Area integration
          y_lsode(neq+1)=x_start                ! x coordinate in logical space
          y_lsode(neq+2)=y_start                ! y coordinate in logical space
          y_lsode(neq+3)=ibl                    ! block label
          CALL get_bfield0(neq,y_lsode,b_xyz,bigr)
          IF(b_xyz(1) > 0)THEN
             i_plus_br = i_plus_br+1
          ELSE
             i_neg_br = i_neg_br+1
          ENDIF
          IF(b_xyz(2) > 0)THEN
             i_plus_bz = i_plus_bz+1
          ELSE
             i_neg_bz = i_neg_bz+1
          ENDIF
        ENDDO
c       write(6,*)item%id,i_plus_br,i_neg_br,i_plus_bz,i_neg_bz
        IF((i_plus_br > 0).AND.(i_neg_br > 0).AND.
     &     (i_plus_bz > 0).AND.(i_neg_bz > 0))THEN
        itemt => item
        EXIT
        ENDIF
        item => item%next
        IF(item%id == 0 ) EXIT
      ENDDO
c
c					First axis guess
      item => itemt
      ibl = item%ib
      rmaxis=0
      zmaxis=0
      DO i=1,item%nodes
        x_start = item%p(1,i)
        y_start = item%p(2,i)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        rmaxis =  rmaxis+ rbc(ibl)%rz%f(1)
        zmaxis =  zmaxis+ rbc(ibl)%rz%f(2)
c-TMP       y_lsode(1:2)=rbc(ibl)%rz%f            ! (R,Z) coordinate
c       y_lsode(3)=0.                         ! PEST angle
c       y_lsode(4)=0.                         ! Area integration
c       y_lsode(neq+1)=x_start                ! x coordinate in logical space
c       y_lsode(neq+2)=y_start                ! y coordinate in logical space
c       y_lsode(neq+3)=ibl                    ! block label
c       CALL get_bfield0(neq,y_lsode,b_xyz,bigr)
c       write(6,*)rbc(ibl)%rz%f
c       write(6,*)b_xyz(1:2)
      ENDDO
      rmaxis =  rmaxis/REAL(item%nodes,r8)
      zmaxis =  zmaxis/REAL(item%nodes,r8)
      WRITE(nim_wr,*)"Initial axis guess ",rmaxis,zmaxis
c-----------------------------------------------------------------------
c     Begin integration 
c-----------------------------------------------------------------------
      item => itemt
      item_pre => itemt
      item => item%face(1)%p  	! Starting at the guess often fails.
      x_start = item%p(1,1)
      y_start = item%p(2,1)
      i_outer = 0
      outer : DO
      IF(local_debug)THEN
        WRITE(7,*)
        IF(local_debug)WRITE(7,*)"# Loop ", i_outer
        WRITE(7,*)rmaxis,zmaxis
        WRITE(7,*)
        WRITE(6,*)rmaxis,zmaxis
      ENDIF
        i_outer=i_outer+1
        IF(i_outer > maxstep)THEN
          WRITE(nim_wr,*)
     &         "Outer loop failure in magnetic axis integration"
          STOP 200
        ENDIF
        istate=1                        !  indicates the first lsode call
        istep0=0
        ibl = item%ib
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        y_lsode(1:2)=rbc(ibl)%rz%f            ! (R,Z) coordinate
        y_lsode(3)=0.                         ! PEST angle
        y_lsode(4)=0.                         ! Area integration
        y_lsode(neq+1)=x_start                ! x coordinate in logical space
        y_lsode(neq+2)=y_start                ! y coordinate in logical space
        y_lsode(neq+3)=ibl                    ! block label
        tstart=ATAN2(y_lsode(2)-zmaxis,y_lsode(1)-rmaxis)
        tend=tstart+twopi
        rwork(1)=tend            ! set maximum angle of integration

        rmaxis_new_min = HUGE(1.0)
        zmaxis_new_min = HUGE(1.0) 
        rmaxis_new_max = -HUGE(1.0)
        zmaxis_new_max = -HUGE(1.0) 
        DO
c          write(6,*)istep0
           CALL dlsode(derivative,neq,y_lsode,tstart,tend,
     &                 itol,rtol,atol,itask,
     &                 istate,iopt,rwork,lrw,iwork,liw,jac,mf)
           IF(local_debug)THEN
               WRITE(7,*)y_lsode(1:2)		! The contour
           ENDIF
           istep0=istep0+1
           IF(y_lsode(1) < rmaxis_new_min)rmaxis_new_min = y_lsode(1)
           IF(y_lsode(1) > rmaxis_new_max)rmaxis_new_max = y_lsode(1)
           IF(y_lsode(2) < zmaxis_new_min)zmaxis_new_min = y_lsode(2)
           IF(y_lsode(2) > zmaxis_new_max)zmaxis_new_max = y_lsode(2)
           IF(ABS(tstart-twopi) < atol)EXIT
c          write(7,*)istep0,tstart,y_lsode(1:2)
c          write(6,*)rwork(2),rwork(11)
           IF(ABS(rwork(11)-rtol) < atol )EXIT outer
           IF(istate < 0)THEN
              WRITE(nim_wr,*)"Error in lsode is ",istate
              STOP 100
           ENDIF
           IF(istep0 > maxstep)EXIT outer
       ENDDO
       rmaxis_new =  0.5*(rmaxis_new_min+rmaxis_new_max)
       zmaxis_new =  0.5*(zmaxis_new_min+zmaxis_new_max)
       p0%point(1) = rmaxis_new*0.1+rmaxis*0.9
       p0%point(2) = zmaxis_new*0.1+zmaxis*0.9
       rmaxis=rmaxis_new
       zmaxis=zmaxis_new
       CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
       CALL refine_cell(p0,item,failure,x_start,y_start,.FALSE.)

      ENDDO outer
      WRITE(nim_wr,*)"Magnetic axis found at",rmaxis,zmaxis
      IF(local_debug)THEN
        write(7,*)"#"
        write(7,*)rmaxis,zmaxis
      ENDIF
      END SUBROUTINE compute_magnetic_axis

      END MODULE magnetic_axis
