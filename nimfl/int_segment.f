c-----------------------------------------------------------------------
c     file int_segment.f
c     contains the subroutine that calls the integration package for
c     one segment of a field line.  recursive calls are allowed to 
c     converge on the crossings of the surfaces of section.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      RECURSIVE SUBROUTINE int_segment(ind_start,ind_end,nstep1,x,y,ib,
     $                                 dep_start,dep_end,plane_type1,
     $                                 levelp,rcp,zcp,percp,any_cross)
      USE local
      USE input0
      USE io
      USE dumpc
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nstep1,levelp
      INTEGER(i4), INTENT(INOUT) :: ib
      REAL(r8), INTENT(IN) :: ind_start,ind_end,rcp,zcp,percp
      REAL(r8), DIMENSION(3), INTENT(IN) :: dep_start
      REAL(r8), INTENT(INOUT) :: x,y
      REAL(r8), DIMENSION(3), INTENT(OUT) :: dep_end
      CHARACTER(*), INTENT(IN) :: plane_type1
      LOGICAL, INTENT(INOUT) :: any_cross
      REAL(r8) :: rnew,znew

c-TMP INTEGER(i4), PARAMETER :: neq=3,liw=20+neq,lrw=22+neq*(9+neq)
      INTEGER(i4), PARAMETER :: neq=3,liw=20,lrw=20+16*neq
      INTEGER(i4), PARAMETER :: nlevel_max=30
      INTEGER(i4), PARAMETER :: nstep_conv=2
      INTEGER(i4) :: iopt,istate,itask,itol,jac,mf,i_lsode
      INTEGER(i4), DIMENSION(liw) :: iwork=0
      REAL(r8) :: atol,rtol,ind,inde,d_ind,rc,zc,perc,perc_adj,xx,yy
      REAL(r8), DIMENSION(2) :: r,z,per
      REAL(r8), DIMENSION(3) :: lower_end
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-10
      REAL(r8), PARAMETER :: rz_tol = 1.e-10
      REAL(r8), PARAMETER :: over_shoot = 0.01
      REAL(r8), EXTERNAL :: flder
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(lrw) :: rwork=0
      LOGICAL :: cross,failure
      TYPE(location_type):: p0
      TYPE(cell_type), POINTER :: item_pre
c-----------------------------------------------------------------------
c     check that the maximum recursion level for converging plane
c     crossings has not been exceeded.
c-----------------------------------------------------------------------
      IF (levelp+1>nlevel_max) THEN
        WRITE(nim_wr,*) 'Maximum segment integration level exceeded.'
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     initialize variables.  note that the y_lsode array contains
c     the dependent variables of the function resulting from the
c     integration (r,z,periodic_coordinate), and additional information
c     regarding the logical position on the nimrod grid, which is
c     updated through the differential equation routine supplied to
c     lsode.

c	y_lsode(1) is the r coordinate
c	y_lsode(2) is the z coordinate
c	y_lsode(3) is the periodic coordinate (phi in toroidal geometry)
c	y_lsode(neq+1) is the x coordinate in logical space
c	y_lsode(neq+2) is the y coordinate in logical space
c	y_lsode(neq+3) is the grid block label,

c     where neq is the number of integrated equations.  the first
c     three are the field line positions.
c-----------------------------------------------------------------------
      y_lsode(1:3)=dep_start
      y_lsode(neq+1)=x
      y_lsode(neq+2)=y
      y_lsode(neq+3)=ib
c-----------------------------------------------------------------------
c     set the integrator parameters.

c	istate=1 :  indicates the first lsode call
c	itask=4 :  normal integration with limited over-shoot (set by
c	           rwork(1) in the i_lsode loop
c	iopt=1 :  optional inputs are used
c	rwork(6) :  set maximum lsode-internal step size
c	rwork(7) :  set minimum lsode-internal step size
c	iwork(6) :  set maximum lsode-internal steps
c	iwork(7) :  set maximum lsode-internal error messages printed
c	mf=10 :  non-stiff Adams method of integration
c	itol=1 :  indicates absolute tolerance is just a scalar
c	rtol :  relative tolerance
c	atol :  absolute tolerance
c-----------------------------------------------------------------------
      istate=1
      itask=4
      iopt=1
      rwork(6)=0.5
      rwork(7)=1.e-12
      iwork(6)=40000
      iwork(7)=2
c-TMP mf=22
      mf=10
      itol=1
      rtol=ind_tol0
      atol=ind_tol0
c-----------------------------------------------------------------------
c     advance the differential equation to ind_end.
c-----------------------------------------------------------------------
      ind=ind_start
      r(1)=dep_start(1)
      z(1)=dep_start(2)
      per(1)=dep_start(3)
      d_ind=(ind_end-ind_start)/nstep1
      DO i_lsode=1,nstep1
        IF (y_lsode(neq+3)==0)THEN 
          WRITE(nim_wr,'(a)') '  Field line exits computational domain.'
          EXIT
        ENDIF
        inde=ind_start+i_lsode*d_ind
        rwork(1)=inde+over_shoot*d_ind
        CALL dlsode(flder,neq,y_lsode,ind,inde,
     &              itol,rtol,atol,itask,istate,iopt,
     &              rwork,lrw,iwork,liw,jac,mf)
        IF (istate<0) THEN
          WRITE(nim_wr,*) 'Error in lsode:istate = ',istate
          EXIT
        ENDIF
        dep_end=y_lsode(1:3)
c-----------------------------------------------------------------------
c	check for surface crossing, and check for convergence if it
c       does.

c       revert to using input and normal correlation of x::r, y::z
c-----------------------------------------------------------------------
        IF (y_lsode(neq+3)/=0) THEN
          r(2)=dep_end(1)
          z(2)=dep_end(2)
          per(2)=dep_end(3)
          CALL plane_cross(r,z,per,cross,rc,zc,perc,perc_adj)
          IF(gridshape == 'rect')THEN
            SELECT CASE(periodicity)
            CASE("y-dir")
              IF(zc.GT.ymax)zc=ymin+MOD(zc-ymin,ymax-ymin)
              IF(zc.LT.ymin)zc=ymax+MOD(zc-ymin,ymax-ymin)
            CASE("both")
              IF(zc.GT.ymax)zc=ymin+MOD(zc-ymin,ymax-ymin)
              IF(zc.LT.ymin)zc=ymax+MOD(zc-ymin,ymax-ymin)
              IF(rc.GT.xmax)rc=xmin+MOD(rc-xmin,xmax-xmin)
              IF(rc.LT.xmin)rc=xmax+MOD(rc-xmin,xmax-xmin)
            END SELECT
          ENDIF
          IF (cross) THEN
            IF (ABS(rc-rcp)<=cross_tol.AND.
     $          ABS(zc-zcp)<=cross_tol.AND.
     $          ABS(perc-percp)<=cross_tol) THEN
c-----------------------------------------------------------------------
c             crossing tolerance met; write location.
c-----------------------------------------------------------------------
              SELECT CASE(plane_type1)
              CASE("r")
                WRITE(temp_unit) REAL(perc_adj,4),REAL(zc,4)
                IF(out_tecplot) 
     $            WRITE(tec2d,*) REAL(perc_adj,4),REAL(zc,4)
              CASE("z")
                IF (gridshape=='rect'.AND.geom=='tor') THEN
                  xx=rc*COS(perc_adj)
                  yy=rc*SIN(perc_adj)
                  WRITE(temp_unit) REAL(xx,4),REAL(yy,4)
                  IF(out_tecplot)
     $               WRITE(tec2d,*) REAL(xx,4),REAL(yy,4)
                ELSE
                  WRITE(temp_unit) REAL(rc,4),REAL(perc_adj,4)
                  IF(out_tecplot)
     $               WRITE(tec2d,*) REAL(rc,4),REAL(perc_adj,4)
                ENDIF
              CASE("periodic")
                WRITE(temp_unit) REAL(rc,4),REAL(zc,4)
                IF(out_tecplot)
     $            WRITE(tec2d,*) REAL(rc,4),REAL(zc,4)
              END SELECT
              any_cross=.TRUE.
            ELSE
c-----------------------------------------------------------------------
c             crossing tolerance is not met.  call int_segment
c             recursively to achieve tolerance.

c       revert to using input and normal correlation of x::r, y::z
c-----------------------------------------------------------------------
c
c-TMP                                           In rect geometry, r,z
c                                               may have periodicity
      znew=z(1)
      rnew=r(1)
      IF(gridshape == 'rect')THEN
        SELECT CASE(periodicity)
        CASE("y-dir")
          IF(z(1) > ymax)znew=ymin+MOD(z(1)-ymin,ymax-ymin)
          IF(z(1) < ymin)znew=ymax+MOD(z(1)-ymin,ymax-ymin)
        CASE("both")
          IF(z(1) > ymax)znew=ymin+MOD(z(1)-ymin,ymax-ymin)
          IF(z(1) < ymin)znew=ymax+MOD(z(1)-ymin,ymax-ymin)
          IF(r(1) > xmax)rnew=xmin+MOD(r(1)-xmin,xmax-xmin)
          IF(r(1) < xmin)rnew=xmax+MOD(r(1)-xmin,xmax-xmin)
        END SELECT
      ENDIF
      r(1)=rnew
      z(1)=znew

            p0%point(1) = r(1)
            p0%point(2) = z(1)
            CALL rz_to_cell(p0,lucky,item_pre,failure,.FALSE.)
            CALL refine_cell(p0,lucky,failure,x,y,.FALSE.)
            ib=lucky%ib
            IF(failure)THEN
               WRITE(out_unit,*)" # FAILURE at ",p0%point(1:2)
               CALL nim_stop("int_segment:ERROR")
            ENDIF
              lower_end=dep_end
              CALL int_segment(inde-d_ind,inde,nstep_conv,x,y,ib,
     $                         (/r(1),z(1),per(1)/),lower_end,
     $                        plane_type1,levelp+1_i4,rc,zc,perc,
     $                        any_cross)
            ENDIF
          ENDIF
        ENDIF
        r(1)=r(2)
        z(1)=z(2)
        per(1)=per(2)
        x=y_lsode(neq+1)
        y=y_lsode(neq+2)
        ib=NINT(y_lsode(neq+3))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE int_segment
