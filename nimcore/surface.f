c-----------------------------------------------------------------------
c     file surface.f
c     subprograms for handling finite element computations on the
c     surface of the domain.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. surface_set.
c     2. surface_dealloc.
c     3. surface_rbl_real_rhs.
c     4. surface_tbl_real_rhs.
c     5. surface_rbl_comp_rhs.
c     6. surface_tbl_comp_rhs.
c-----------------------------------------------------------------------
c     module declaration.
c-----------------------------------------------------------------------
      MODULE surface
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE vector_type_mod
      IMPLICIT NONE

      INTERFACE surface_real_rhs
        MODULE PROCEDURE surface_rbl_real_rhs,surface_tbl_real_rhs
      END INTERFACE

      INTERFACE surface_comp_rhs
        MODULE PROCEDURE surface_rbl_comp_rhs,surface_tbl_comp_rhs
      END INTERFACE

c-----------------------------------------------------------------------
c     private data for gaussian quadrature.
c-----------------------------------------------------------------------
      INTEGER(i4), PRIVATE :: ng
      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: xg,wg

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. surface_set.
c     set the locations and weights for quadratures.
c-----------------------------------------------------------------------
      SUBROUTINE surface_set(ng_input,poly_degree,int_formula)

      INTEGER(i4), INTENT(IN) :: ng_input,poly_degree
      CHARACTER(8), INTENT(IN) :: int_formula
c-----------------------------------------------------------------------
c     set number of quadrature points and weights according to input.
c     the number of points is now adjusted automatically with
c     poly_degree.
c-----------------------------------------------------------------------
      ng=ng_input+poly_degree-1
      ALLOCATE(xg(ng),wg(ng))
      IF (int_formula=='lobatto') THEN
        CALL lobleg(0._r8,1._r8,xg,wg,ng)
      ELSE
        CALL gauleg(0._r8,1._r8,xg,wg,ng)
      ENDIF
c-----------------------------------------------------------------------
c     terminate surface_set.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_set
c-----------------------------------------------------------------------
c     subprogram 2. surface_dealloc.
c     deallocate module arrays.
c-----------------------------------------------------------------------
      SUBROUTINE surface_dealloc

      DEALLOCATE(xg,wg)

      RETURN
      END SUBROUTINE surface_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. surface_rbl_real_rhs.
c     performs numerical integrations for surface integral contributions
c     to the rhs of an equation producing real data, where the
c     integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE surface_rbl_real_rhs(rb,rhs,get_integrand,nq,
     $                                intxys,intxyp,intxyn,h_side,
     $                                met_spl,poly_degree,geom)
      USE lagr_quad_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,poly_degree
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: intxys,intxyp,intxyn
      LOGICAL, INTENT(IN) :: h_side
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(vector_type), INTENT(INOUT) :: rhs
      CHARACTER(*), INTENT(IN) :: met_spl,geom

      INTEGER(i4) :: ix,iy,ixc,iyc,iq,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,n_grid,n_horz,n_vert
      REAL(r8) :: x,y,dx,dy,bigr,jac,dxdr,dydr,dxdz,dydz,ds
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(2) :: tang,norm,drzdx,drzdy,rz
      REAL(r8), DIMENSION(1:poly_degree+1) :: alpha,dalpha
      INTEGER(i4), DIMENSION(2) :: del_xy
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: integrand
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhsg
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: rhsh,rhsv
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,rb,tb,x,y,bigr,norm,
     $                           ijcell,alpha,dxdr,dydr,dxdz,dydz)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:), INTENT(OUT) :: integrand 
        TYPE(rblock_type), INTENT(INOUT) :: rb
        TYPE(tblock_type), INTENT(INOUT) :: tb
        REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
        REAL(r8), DIMENSION(2), INTENT(IN) :: norm
        REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
        INTEGER(i4), DIMENSION(2) :: ijcell
        
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for lagr_1D external routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
        USE local
        INTEGER(i4), INTENT(IN) :: pd,dmode
        REAL(r8), INTENT(IN) :: x
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
        END SUBROUTINE lagr_1D
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and determine the cell indices.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      del_xy=intxyn-intxyp
      IF (del_xy(1)==0) THEN  !  vertical side.
        ix=intxyn(1)
        iy=MAX(intxyn(2),intxyp(2))
        ixc=MAX(ix,1_i4)
        iyc=iy
      ELSE  !  horizontal side.
        ix=MAX(intxyn(1),intxyp(1))
        iy=intxyn(2)
        ixc=ix
        iyc=MAX(iy,1_i4)
      ENDIF
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=2
        rhsg=>rhs%arr
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid
      IF (ASSOCIATED(rhs%arrh).AND.h_side) THEN
        IF (del_xy(1)>0) THEN
          start_horz=2
        ELSE
          start_horz=poly_degree
        ENDIF
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
      ELSE
        n_horz=0
      ENDIF
      iv=iv+n_horz
      IF (ASSOCIATED(rhs%arrv).AND..NOT.h_side) THEN
        IF (del_xy(2)>0) THEN
          start_vert=2
        ELSE
          start_vert=poly_degree
        ENDIF
        n_vert=SIZE(rhs%arrv,2)
        rhsv=>rhs%arrv
      ELSE
        n_vert=0
      ENDIF
      iv=iv+n_vert
c-----------------------------------------------------------------------
c     flag the tblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integrand(nq,iv))
c-----------------------------------------------------------------------
c     start loop over Gaussian quadrature points along element surface.
c-----------------------------------------------------------------------
      DO ig=1,ng
        x=xg(ig)*del_xy(1)+intxyp(1)
        y=xg(ig)*del_xy(2)+intxyp(2)
        dx=x-ix+1
        dy=y-iy+1
c-----------------------------------------------------------------------
c       evaluate the necessary geometric quantities at the quadrature
c       point.
c
c-PRE   grid will be handled differently.
c-----------------------------------------------------------------------
        CALL lagr_1D(poly_degree,xg(ig),alpha,dalpha,0_i4)
        SELECT CASE(met_spl)
        CASE('pcnst')
          rz=0.25*(rb%rz%fs(:,ixc-1,iyc-1)
     $            +rb%rz%fs(:,ixc  ,iyc-1)
     $            +rb%rz%fs(:,ixc-1,iyc  )
     $            +rb%rz%fs(:,ixc  ,iyc  ))
          drzdx=0.5*(rb%rz%fs(:,ixc  ,iyc  )
     $              -rb%rz%fs(:,ixc-1,iyc  )
     $              +rb%rz%fs(:,ixc  ,iyc-1)
     $              -rb%rz%fs(:,ixc-1,iyc-1))
          drzdy=0.5*(rb%rz%fs(:,ixc  ,iyc  )
     $              -rb%rz%fs(:,ixc  ,iyc-1)
     $              +rb%rz%fs(:,ixc-1,iyc  )
     $              -rb%rz%fs(:,ixc-1,iyc-1))
        CASE('liner','linear','bilinear')
          rz=(1-dx)*(1-dy)*rb%rz%fs(:,ixc-1,iyc-1)
     $      +   dx *(1-dy)*rb%rz%fs(:,ixc  ,iyc-1)
     $      +(1-dx)*   dy *rb%rz%fs(:,ixc-1,iyc  )
     $      +   dx *   dy *rb%rz%fs(:,ixc  ,iyc  )
          drzdx=   dy *(rb%rz%fs(:,ixc  ,iyc  )
     $                 -rb%rz%fs(:,ixc-1,iyc  ))
     $         +(1-dy)*(rb%rz%fs(:,ixc  ,iyc-1)
     $                 -rb%rz%fs(:,ixc-1,iyc-1))
          drzdy=   dx *(rb%rz%fs(:,ixc  ,iyc  )
     $                 -rb%rz%fs(:,ixc  ,iyc-1))
     $         +(1-dx)*(rb%rz%fs(:,ixc-1,iyc  )
     $                 -rb%rz%fs(:,ixc-1,iyc-1))
        CASE('iso')
          CALL lagr_quad_eval(rb%rz,x,y,1_i4)
          rz=rb%rz%f
          drzdx=rb%rz%fx
          drzdy=rb%rz%fy
        CASE DEFAULT
          CALL nim_stop('Surface_rbl_real_rhs: '//TRIM(met_spl)//
     $                  ' is not a valid option for met_spl.')
        END SELECT
        jac=drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac
        dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac
        dydz= drzdx(1)/jac
        IF (geom=='tor') THEN
          bigr=rz(1)
        ELSE
          bigr=1
        ENDIF
        tang=drzdx*del_xy(1)+drzdy*del_xy(2)
        ds=SQRT(SUM(tang**2))
        norm=(/tang(2),-tang(1)/)/ds
        ds=ds*bigr
c-----------------------------------------------------------------------
c       evaluate the integrand at the quadrature point.
c-----------------------------------------------------------------------
        CALL get_integrand(integrand,rb,tdum,x,y,bigr,norm,
     $                     intxys,alpha,dxdr,dydr,dxdz,dydz)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from this
c       quadrature point into the correct arrays.
c       grid vertex-centered bases first.
c-----------------------------------------------------------------------
        IF (n_grid==2) THEN
          DO iq=1,nq
            rhsg(iq,intxyp(1),intxyp(2))=
     $        rhsg(iq,intxyp(1),intxyp(2))
     $        +wg(ig)*ds*integrand(iq,1)
            rhsg(iq,intxyn(1),intxyn(2))=
     $        rhsg(iq,intxyn(1),intxyn(2))
     $        +wg(ig)*ds*integrand(iq,poly_degree+1)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          iv=start_horz
          DO ib=1,n_horz
            DO iq=1,nq
              rhsh(iq,ib,intxys(1),intxys(2))=
     $          rhsh(iq,ib,intxys(1),intxys(2))
     $          +wg(ig)*ds*integrand(iq,iv)
            ENDDO
            iv=iv+del_xy(1)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          iv=start_vert
          DO ib=1,n_vert
            DO iq=1,nq
              rhsv(iq,ib,intxys(1),intxys(2))=
     $          rhsv(iq,ib,intxys(1),intxys(2))
     $          +wg(ig)*ds*integrand(iq,iv)
            ENDDO
            iv=iv+del_xy(2)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     finish loops over Gaussian quadrature points.
c-----------------------------------------------------------------------
      ENDDO
      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_rbl_real_rhs
c-----------------------------------------------------------------------
c     subprogram 4. surface_tbl_real_rhs.
c     performs numerical integrations for surface integral contributions
c     to the rhs of an equation producing real data, where the
c     integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE surface_tbl_real_rhs(tb,rhs,get_integrand,nq,
     $                                intxys,intxyp,intxyn,poly_degree,
     $                                geom)
      USE lagr_quad_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,poly_degree
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: intxys,intxyp,intxyn
      TYPE(tblock_type), INTENT(INOUT) :: tb
      TYPE(vector_type), INTENT(INOUT) :: rhs
      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: ix,iy,iq,jf,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,n_grid,n_horz,n_vert
      REAL(r8) :: x,y,dx,dy,bigr,jac,dxdr,dydr,dxdz,dydz,ds
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (rblock_type) :: rdum

      REAL(r8), DIMENSION(2) :: tang,norm,del_rz
      REAL(r8), DIMENSION(1:poly_degree+1) :: alpha,dalpha
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: integrand
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhsg
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: rhsh
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,rb,tb,x,y,bigr,norm,
     $                           ijcell,alpha,dxdr,dydr,dxdz,dydz)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:), INTENT(OUT) :: integrand 
        TYPE(rblock_type), INTENT(INOUT) :: rb
        TYPE(tblock_type), INTENT(INOUT) :: tb
        REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
        REAL(r8), DIMENSION(2), INTENT(IN) :: norm
        REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
        INTEGER(i4), DIMENSION(2) :: ijcell
        
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for lagr_1D external routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
        USE local
        INTEGER(i4), INTENT(IN) :: pd,dmode
        REAL(r8), INTENT(IN) :: x
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
        END SUBROUTINE lagr_1D
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=2
        rhsg=>rhs%arr
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid
      IF (ASSOCIATED(rhs%arrh)) THEN
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
      ELSE
        n_horz=0
      ENDIF
      iv=iv+n_horz
c-----------------------------------------------------------------------
c     flag the rblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      rdum%mx=-1
      ALLOCATE(integrand(nq,iv))
c-----------------------------------------------------------------------
c     start loop over Gaussian quadrature points along element surface.
c-----------------------------------------------------------------------
      DO ig=1,ng
        x=   xg(ig) *tb%tgeom%xs(intxyn(1))+
     $    (1-xg(ig))*tb%tgeom%xs(intxyp(1))
        y=   xg(ig) *tb%tgeom%ys(intxyn(1))+
     $    (1-xg(ig))*tb%tgeom%ys(intxyp(1))
c-----------------------------------------------------------------------
c       evaluate the necessary geometric quantities at the quadrature
c       point.
c
c-PRE   grid will be handled differently.
c-----------------------------------------------------------------------
        CALL lagr_1D(poly_degree,xg(ig),alpha,dalpha,0_i4)
        dxdr=1
        dydr=0
        dxdz=0
        dydz=1
        IF (geom=='tor') THEN
          bigr=x
        ELSE
          bigr=1
        ENDIF
        tang=(/tb%tgeom%xs(intxyn(1))-tb%tgeom%xs(intxyp(1)),
     $         tb%tgeom%ys(intxyn(1))-tb%tgeom%ys(intxyp(1))/)
        ds=SQRT(SUM(tang**2))
        norm=(/tang(2),-tang(1)/)/ds
        ds=ds*bigr
c-----------------------------------------------------------------------
c       evaluate the integrand at the quadrature point.
c-----------------------------------------------------------------------
        CALL get_integrand(integrand,rdum,tb,x,y,bigr,norm,
     $                     intxys,alpha,dxdr,dydr,dxdz,dydz)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from this
c       quadrature point into the correct arrays.
c       grid vertex-centered bases first.
c-----------------------------------------------------------------------
        IF (n_grid==2) THEN
          DO iq=1,nq
            rhsg(iq,intxyp(1),0)=rhsg(iq,intxyp(1),0)
     $        +wg(ig)*ds*integrand(iq,1)
            rhsg(iq,intxyn(1),0)=rhsg(iq,intxyn(1),0)
     $        +wg(ig)*ds*integrand(iq,poly_degree+1)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c-PRE   horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          iv=2
          DO ib=1,n_horz
            DO iq=1,nq
              rhsh(iq,ib,intxys(1),0)=rhsh(iq,ib,intxys(1),0)
     $          +wg(ig)*ds*integrand(iq,iv)
            ENDDO
            iv=iv+1
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     finish loops over Gaussian quadrature points.
c-----------------------------------------------------------------------
      ENDDO
      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_tbl_real_rhs
c-----------------------------------------------------------------------
c     subprogram 5. surface_rbl_comp_rhs.
c     performs numerical integrations for surface integral contributions
c     to the rhs of an equation producing complex data, where the
c     integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE surface_rbl_comp_rhs(rb,rhs,get_integrand,nq,nfour,
     $                                intxys,intxyp,intxyn,h_side,
     $                                met_spl,poly_degree,geom)
      USE lagr_quad_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,nfour,poly_degree
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: intxys,intxyp,intxyn
      LOGICAL, INTENT(IN) :: h_side
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(cvector_type), INTENT(INOUT) :: rhs
      CHARACTER(*), INTENT(IN) :: met_spl,geom

      INTEGER(i4) :: ix,iy,ixc,iyc,iq,jf,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,n_grid,n_horz,n_vert
      REAL(r8) :: x,y,dx,dy,bigr,jac,dxdr,dydr,dxdz,dydz,ds
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(2) :: tang,norm,drzdx,drzdy,rz
      REAL(r8), DIMENSION(1:poly_degree+1) :: alpha,dalpha
      INTEGER(i4), DIMENSION(2) :: del_xy
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: integrand
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhsh,rhsv
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,rb,tb,x,y,bigr,norm,
     $                           ijcell,alpha,dxdr,dydr,dxdz,dydz)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: integrand 
        TYPE(rblock_type), INTENT(INOUT) :: rb
        TYPE(tblock_type), INTENT(INOUT) :: tb
        REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
        REAL(r8), DIMENSION(2), INTENT(IN) :: norm
        REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
        INTEGER(i4), DIMENSION(2) :: ijcell
        
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for lagr_1D external routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
        USE local
        INTEGER(i4), INTENT(IN) :: pd,dmode
        REAL(r8), INTENT(IN) :: x
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
        END SUBROUTINE lagr_1D
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and determine the cell indices.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      del_xy=intxyn-intxyp
      IF (del_xy(1)==0) THEN  !  vertical side.
        ix=intxyn(1)
        iy=MAX(intxyn(2),intxyp(2))
        ixc=MAX(ix,1_i4)
        iyc=iy
      ELSE  !  horizontal side.
        ix=MAX(intxyn(1),intxyp(1))
        iy=intxyn(2)
        ixc=ix
        iyc=MAX(iy,1_i4)
      ENDIF
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=2
        rhsg=>rhs%arr
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid
      IF (ASSOCIATED(rhs%arrh).AND.h_side) THEN
        IF (del_xy(1)>0) THEN
          start_horz=2
        ELSE
          start_horz=poly_degree
        ENDIF
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
      ELSE
        n_horz=0
      ENDIF
      iv=iv+n_horz
      IF (ASSOCIATED(rhs%arrv).AND..NOT.h_side) THEN
        IF (del_xy(2)>0) THEN
          start_vert=2
        ELSE
          start_vert=poly_degree
        ENDIF
        n_vert=SIZE(rhs%arrv,2)
        rhsv=>rhs%arrv
      ELSE
        n_vert=0
      ENDIF
      iv=iv+n_vert
c-----------------------------------------------------------------------
c     flag the tblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integrand(nq,iv,nfour))
c-----------------------------------------------------------------------
c     start loop over Gaussian quadrature points along element surface.
c-----------------------------------------------------------------------
      DO ig=1,ng
        x=xg(ig)*del_xy(1)+intxyp(1)
        y=xg(ig)*del_xy(2)+intxyp(2)
        dx=x-ix+1
        dy=y-iy+1
c-----------------------------------------------------------------------
c       evaluate the necessary geometric quantities at the quadrature
c       point.
c
c-PRE   grid will be handled differently.
c-----------------------------------------------------------------------
        CALL lagr_1D(poly_degree,xg(ig),alpha,dalpha,0_i4)
        SELECT CASE(met_spl)
        CASE('pcnst')
          rz=0.25*(rb%rz%fs(:,ixc-1,iyc-1)
     $            +rb%rz%fs(:,ixc  ,iyc-1)
     $            +rb%rz%fs(:,ixc-1,iyc  )
     $            +rb%rz%fs(:,ixc  ,iyc  ))
          drzdx=0.5*(rb%rz%fs(:,ixc  ,iyc  )
     $              -rb%rz%fs(:,ixc-1,iyc  )
     $              +rb%rz%fs(:,ixc  ,iyc-1)
     $              -rb%rz%fs(:,ixc-1,iyc-1))
          drzdy=0.5*(rb%rz%fs(:,ixc  ,iyc  )
     $              -rb%rz%fs(:,ixc  ,iyc-1)
     $              +rb%rz%fs(:,ixc-1,iyc  )
     $              -rb%rz%fs(:,ixc-1,iyc-1))
        CASE('liner','linear','bilinear')
          rz=(1-dx)*(1-dy)*rb%rz%fs(:,ixc-1,iyc-1)
     $      +   dx *(1-dy)*rb%rz%fs(:,ixc  ,iyc-1)
     $      +(1-dx)*   dy *rb%rz%fs(:,ixc-1,iyc  )
     $      +   dx *   dy *rb%rz%fs(:,ixc  ,iyc  )
          drzdx=   dy *(rb%rz%fs(:,ixc  ,iyc  )
     $                 -rb%rz%fs(:,ixc-1,iyc  ))
     $         +(1-dy)*(rb%rz%fs(:,ixc  ,iyc-1)
     $                 -rb%rz%fs(:,ixc-1,iyc-1))
          drzdy=   dx *(rb%rz%fs(:,ixc  ,iyc  )
     $                 -rb%rz%fs(:,ixc  ,iyc-1))
     $         +(1-dx)*(rb%rz%fs(:,ixc-1,iyc  )
     $                 -rb%rz%fs(:,ixc-1,iyc-1))
        CASE('iso')
          CALL lagr_quad_eval(rb%rz,x,y,1_i4)
          rz=rb%rz%f
          drzdx=rb%rz%fx
          drzdy=rb%rz%fy
        CASE DEFAULT
          CALL nim_stop('Surface_rbl_comp_rhs: '//TRIM(met_spl)//
     $                  ' is not a valid option for met_spl.')
        END SELECT
        jac=drzdx(1)*drzdy(2)-drzdx(2)*drzdy(1)
        dxdr= drzdy(2)/jac
        dydr=-drzdx(2)/jac
        dxdz=-drzdy(1)/jac
        dydz= drzdx(1)/jac
        IF (geom=='tor') THEN
          bigr=rz(1)
        ELSE
          bigr=1
        ENDIF
        tang=drzdx*del_xy(1)+drzdy*del_xy(2)
        ds=SQRT(SUM(tang**2))
        norm=(/tang(2),-tang(1)/)/ds
        ds=ds*bigr
c-----------------------------------------------------------------------
c       evaluate the integrand at the quadrature point.
c-----------------------------------------------------------------------
        CALL get_integrand(integrand,rb,tdum,x,y,bigr,norm,
     $                     intxys,alpha,dxdr,dydr,dxdz,dydz)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from this
c       quadrature point into the correct arrays.
c       grid vertex-centered bases first.
c-----------------------------------------------------------------------
        IF (n_grid==2) THEN
          DO jf=1,nfour
            DO iq=1,nq
              rhsg(iq,intxyp(1),intxyp(2),jf)=
     $          rhsg(iq,intxyp(1),intxyp(2),jf)
     $          +wg(ig)*ds*integrand(iq,1,jf)
              rhsg(iq,intxyn(1),intxyn(2),jf)=
     $          rhsg(iq,intxyn(1),intxyn(2),jf)
     $          +wg(ig)*ds*integrand(iq,poly_degree+1,jf)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          DO jf=1,nfour
            iv=start_horz
            DO ib=1,n_horz
              DO iq=1,nq
                rhsh(iq,ib,intxys(1),intxys(2),jf)=
     $            rhsh(iq,ib,intxys(1),intxys(2),jf)
     $            +wg(ig)*ds*integrand(iq,iv,jf)
              ENDDO
              iv=iv+del_xy(1)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          DO jf=1,nfour
            iv=start_vert
            DO ib=1,n_vert
              DO iq=1,nq
                rhsv(iq,ib,intxys(1),intxys(2),jf)=
     $            rhsv(iq,ib,intxys(1),intxys(2),jf)
     $            +wg(ig)*ds*integrand(iq,iv,jf)
              ENDDO
              iv=iv+del_xy(2)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     finish loops over Gaussian quadrature points.
c-----------------------------------------------------------------------
      ENDDO
      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_rbl_comp_rhs
c-----------------------------------------------------------------------
c     subprogram 6. surface_tbl_comp_rhs.
c     performs numerical integrations for surface integral contributions
c     to the rhs of an equation producing complex data, where the
c     integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE surface_tbl_comp_rhs(tb,rhs,get_integrand,nq,nfour,
     $                                intxys,intxyp,intxyn,poly_degree,
     $                                geom)
      USE lagr_quad_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,nfour,poly_degree
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: intxys,intxyp,intxyn
      TYPE(tblock_type), INTENT(INOUT) :: tb
      TYPE(cvector_type), INTENT(INOUT) :: rhs
      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: ix,iy,iq,jf,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,n_grid,n_horz,n_vert
      REAL(r8) :: x,y,dx,dy,bigr,jac,dxdr,dydr,dxdz,dydz,ds
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (rblock_type) :: rdum

      REAL(r8), DIMENSION(2) :: tang,norm,del_rz
      REAL(r8), DIMENSION(1:poly_degree+1) :: alpha,dalpha
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: integrand
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhsh
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,rb,tb,x,y,bigr,norm,
     $                           ijcell,alpha,dxdr,dydr,dxdz,dydz)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: integrand 
        TYPE(rblock_type), INTENT(INOUT) :: rb
        TYPE(tblock_type), INTENT(INOUT) :: tb
        REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
        REAL(r8), DIMENSION(2), INTENT(IN) :: norm
        REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
        INTEGER(i4), DIMENSION(2) :: ijcell
        
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for lagr_1D external routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
        USE local
        INTEGER(i4), INTENT(IN) :: pd,dmode
        REAL(r8), INTENT(IN) :: x
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
        END SUBROUTINE lagr_1D
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=2
        rhsg=>rhs%arr
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid
      IF (ASSOCIATED(rhs%arrh)) THEN
        n_horz=SIZE(rhs%arrh,2)
        rhsh=>rhs%arrh
      ELSE
        n_horz=0
      ENDIF
      iv=iv+n_horz
c-----------------------------------------------------------------------
c     flag the rblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      rdum%mx=-1
      ALLOCATE(integrand(nq,iv,nfour))
c-----------------------------------------------------------------------
c     start loop over Gaussian quadrature points along element surface.
c-----------------------------------------------------------------------
      DO ig=1,ng
        x=   xg(ig) *tb%tgeom%xs(intxyn(1))+
     $    (1-xg(ig))*tb%tgeom%xs(intxyp(1))
        y=   xg(ig) *tb%tgeom%ys(intxyn(1))+
     $    (1-xg(ig))*tb%tgeom%ys(intxyp(1))
c-----------------------------------------------------------------------
c       evaluate the necessary geometric quantities at the quadrature
c       point.
c
c-PRE   grid will be handled differently.
c-----------------------------------------------------------------------
        CALL lagr_1D(poly_degree,xg(ig),alpha,dalpha,0_i4)
        dxdr=1
        dydr=0
        dxdz=0
        dydz=1
        IF (geom=='tor') THEN
          bigr=x
        ELSE
          bigr=1
        ENDIF
        tang=(/tb%tgeom%xs(intxyn(1))-tb%tgeom%xs(intxyp(1)),
     $         tb%tgeom%ys(intxyn(1))-tb%tgeom%ys(intxyp(1))/)
        ds=SQRT(SUM(tang**2))
        norm=(/tang(2),-tang(1)/)/ds
        ds=ds*bigr
c-----------------------------------------------------------------------
c       evaluate the integrand at the quadrature point.
c-----------------------------------------------------------------------
        CALL get_integrand(integrand,rdum,tb,x,y,bigr,norm,
     $                     intxys,alpha,dxdr,dydr,dxdz,dydz)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from this
c       quadrature point into the correct arrays.
c       grid vertex-centered bases first.
c-----------------------------------------------------------------------
        IF (n_grid==2) THEN
          DO jf=1,nfour
            DO iq=1,nq
              rhsg(iq,intxyp(1),0,jf)=rhsg(iq,intxyp(1),0,jf)
     $          +wg(ig)*ds*integrand(iq,1,jf)
              rhsg(iq,intxyn(1),0,jf)=rhsg(iq,intxyn(1),0,jf)
     $          +wg(ig)*ds*integrand(iq,poly_degree+1,jf)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c-PRE   horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          DO jf=1,nfour
            iv=2
            DO ib=1,n_horz
              DO iq=1,nq
                rhsh(iq,ib,intxys(1),0,jf)=rhsh(iq,ib,intxys(1),0,jf)
     $            +wg(ig)*ds*integrand(iq,iv,jf)
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     finish loops over Gaussian quadrature points.
c-----------------------------------------------------------------------
      ENDDO
      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_tbl_comp_rhs
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE surface
