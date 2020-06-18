c-----------------------------------------------------------------------
c     file generic_evals.f
c     module containing routines that find finite-element function
c     values and derivatives from appropriate splining and interpolation
c     routines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. generic_evals.
c     1. generic_eq_all_eval
c     2. generic_3D_all_eval.
c     3. generic_2D_all_eval.
c     4. generic_3D_modal_all_eval.
c     5. generic_2D_modal_all_eval.
c     6. generic_3D_ptr_set.
c     7. generic_2D_ptr_set.
c     8. generic_alpha_eval.
c     9. generic_3D_eval.
c     10. generic_2D_eval.
c     11. generic_3D_modal_eval.
c     12. generic_2D_modal_eval.
c-----------------------------------------------------------------------
c     0. module declarations for generic_evals.
c-----------------------------------------------------------------------
      MODULE generic_evals
      USE local
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     subprogram name interfaces
c-----------------------------------------------------------------------
      INTERFACE generic_all_eval
        MODULE PROCEDURE generic_2D_all_eval,generic_3D_all_eval,
     $                   generic_eq_all_eval,generic_2D_modal_all_eval,
     $                   generic_3D_modal_all_eval
      END INTERFACE

      INTERFACE generic_ptr_set
        MODULE PROCEDURE generic_2D_ptr_set,generic_3D_ptr_set
      END INTERFACE

      INTERFACE generic_eval
        MODULE PROCEDURE generic_2D_eval,generic_3D_eval,
     $                   generic_2D_modal_eval,generic_3D_modal_eval
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. generic_eq_all_eval.
c     get equilibrium field data at the requested quadrature points
c     in each cell in a block.
c
c     this routines evaluates data at the set of points specified by dx 
c     and dy for rblocks and at each of the quadrature-points (1:inode) 
c     for tblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_eq_all_eval(bc,tl,dxdr,dydr,dxdz,dydz,dx,dy,
     $                               tg,inode,data,data_r,data_z,
     $                               d_order)
      USE bicube
      USE tri_linear

      TYPE(bicube_type), INTENT(INOUT) :: bc
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(:), INTENT(IN) :: dx,dy
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: data,data_r,data_z

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dat,data_x,data_y
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: xg,yg
      INTEGER(i4) :: mx,my,iq,nq,ig,ng,ix,iy,ipol
c-----------------------------------------------------------------------
c     decide which block is the dummy.
c     rblock data:  splines are evaluated before the first time
c     step for efficiency, so just recover stored data.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        mx=bc%mx
        my=bc%my
        nq=bc%nqty
        nq=bc%nqty
        ng=SIZE(dx)
c-----------------------------------------------------------------------
c       allocate arrays for derivatives.
c-----------------------------------------------------------------------
        ALLOCATE(dat(nq,mx,my))
        IF (d_order>0) THEN
          ALLOCATE(data_x(nq,mx,my))
          ALLOCATE(data_y(nq,mx,my))
        ELSE
          ALLOCATE(data_x(1,1,1))
          ALLOCATE(data_y(1,1,1))
        ENDIF
c-----------------------------------------------------------------------
c       evaluate bicubic spline.
c-----------------------------------------------------------------------
        DO ig=1,ng
          xg=dx(ig)
          yg=dy(ig)
          CALL bicube_all_eval(bc,xg,yg,dat,data_x,data_y,
     $                         dv,dv,dv,d_order)
          ipol=1
          DO iy=1,my
            DO ix=1,mx
              data(:,ig,ipol)=dat(:,ix,iy)
              ipol=ipol+1
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         use the chain rule to get derivatives with respect to the 
c         cylindrical coordinates.
c-----------------------------------------------------------------------
          IF (d_order>0) THEN
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                data_r(:,ig,ipol)=dxdr(ig,ipol)*data_x(:,ix,iy)
     $                           +dydr(ig,ipol)*data_y(:,ix,iy)
                data_z(:,ig,ipol)=dxdz(ig,ipol)*data_x(:,ix,iy)
     $                           +dydz(ig,ipol)*data_y(:,ix,iy)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(bc%cmats)
c-----------------------------------------------------------------------
c     tblock data:  equilibrium quantities in tblocks are represented
c     by 2D linear elements.
c-----------------------------------------------------------------------
      ELSE block_type
        mx=tg%mcell
        nq=tl%nqty
        ng=inode
        ALLOCATE(dat(nq,mx,1))
        ALLOCATE(data_x(nq,mx,1))
        ALLOCATE(data_y(nq,mx,1))
        DO ig=1,ng
          CALL tri_linear_all_eval(tl,tg,ig,d_order,dat,data_x,data_y)
          data(:,ig,:)=dat(:,:,1)
          data_r(:,ig,:)=data_x(:,:,1)
          data_z(:,ig,:)=data_y(:,:,1)
        ENDDO
      ENDIF block_type
      DEALLOCATE(dat,data_x,data_y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_eq_all_eval
c-----------------------------------------------------------------------
c     subprogram 2. generic_3D_all_eval.
c     get complex dependent-variable field data at the requested
c     quadrature points in each cell in a block.
c
c     this routines evaluates data at the set of points specified by dx 
c     and dy for rblocks and at each of the quadrature-points (1:inode) 
c     for tblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_3D_all_eval(laq,tl,dxdr,dydr,dxdz,dydz,dx,dy,
     $                               tg,inode,data,data_r,data_z,
     $                               d_order)
      USE lagr_quad_mod
      USE tri_linear

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      TYPE(tri_linear_type), INTENT(IN) :: tl
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(:), INTENT(IN) :: dx,dy
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: data,data_r,data_z

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: dat,data_x,data_y
      REAL(r8) :: xg,yg
      INTEGER(i4) :: mx,my,iq,nq,ig,ng,ix,iy,ipol,jf,nf
c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        mx=laq%mx
        my=laq%my
        nq=laq%nqty
        nf=laq%nfour
        ng=SIZE(dx)
c-----------------------------------------------------------------------
c       allocate arrays for derivatives.
c-----------------------------------------------------------------------
        ALLOCATE(dat(nq,mx,my,nf))
        IF (d_order>0) THEN
          ALLOCATE(data_x(nq,mx,my,nf))
          ALLOCATE(data_y(nq,mx,my,nf))
        ELSE
          ALLOCATE(data_x(1,1,1,1))
          ALLOCATE(data_y(1,1,1,1))
        ENDIF
c-----------------------------------------------------------------------
c       rblock data:  call the lagr_quad interpolation routine at
c       the quadrature points.
c-----------------------------------------------------------------------
        DO ig=1,ng
          xg=dx(ig)
          yg=dy(ig)
          CALL lagr_quad_all_eval(laq,xg,yg,dat,data_x,data_y,d_order)
          DO jf=1,nf
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                data(:,ig,ipol,jf)=dat(:,ix,iy,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         use the chain rule to get derivatives with respect to the 
c         cylindrical coordinates.
c-----------------------------------------------------------------------
          IF (d_order>0) THEN
            DO jf=1,nf
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  data_r(:,ig,ipol,jf)=dxdr(ig,ipol)*data_x(:,ix,iy,jf)
     $                                +dydr(ig,ipol)*data_y(:,ix,iy,jf)
                  data_z(:,ig,ipol,jf)=dxdz(ig,ipol)*data_x(:,ix,iy,jf)
     $                                +dydz(ig,ipol)*data_y(:,ix,iy,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE block_type
c-----------------------------------------------------------------------
c       tblock data:  equilibrium quantities in tblocks are represented
c       by 2D linear elements.
c-----------------------------------------------------------------------
        mx=tg%mcell
        nq=tl%nqty
        ng=inode
        ALLOCATE(dat(nq,mx,1,nf))
        ALLOCATE(data_x(nq,mx,1,nf))
        ALLOCATE(data_y(nq,mx,1,nf))
        DO ig=1,ng
          CALL tri_linear_all_eval(tl,tg,ig,d_order,dat,data_x,data_y)
          data(:,ig,:,:)=dat(:,:,1,:)
          data_r(:,ig,:,:)=data_x(:,:,1,:)
          data_z(:,ig,:,:)=data_y(:,:,1,:)
        ENDDO
      ENDIF block_type
      DEALLOCATE(dat,data_x,data_y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_3D_all_eval
c-----------------------------------------------------------------------
c     subprogram 3. generic_2D_all_eval.
c     get real dependent-variable field data at the requested quadrature
c     points in each cell in a block.
c
c     this routines evaluates data at the set of points specified by dx 
c     and dy for rblocks and at each of the quadrature-points (1:inode) 
c     for tblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_2D_all_eval(laq,tl,dxdr,dydr,dxdz,dydz,dx,dy,
     $                               tg,inode,data,data_r,data_z,
     $                               d_order)
      USE lagr_quad_mod
      USE tri_linear

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(:), INTENT(IN) :: dx,dy
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: data,data_r,data_z

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dat,data_x,data_y
      REAL(r8) :: xg,yg
      INTEGER(i4) :: mx,my,iq,nq,ig,ng,ix,iy,ipol
c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        mx=laq%mx
        my=laq%my
        nq=laq%nqty
        ng=SIZE(dx)
c-----------------------------------------------------------------------
c       allocate arrays for derivatives.
c-----------------------------------------------------------------------
        ALLOCATE(dat(nq,mx,my))
        IF (d_order>0) THEN
          ALLOCATE(data_x(nq,mx,my))
          ALLOCATE(data_y(nq,mx,my))
        ELSE
          ALLOCATE(data_x(1,1,1))
          ALLOCATE(data_y(1,1,1))
        ENDIF
c-----------------------------------------------------------------------
c       rblock data:  call the lagr_quad interpolation routine at
c       the quadrature points.
c-----------------------------------------------------------------------
        DO ig=1,ng
          xg=dx(ig)
          yg=dy(ig)
          CALL lagr_quad_all_eval(laq,xg,yg,dat,data_x,data_y,d_order)
          ipol=1
          DO iy=1,my
            DO ix=1,mx
              data(:,ig,ipol)=dat(:,ix,iy)
              ipol=ipol+1
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         use the chain rule to get derivatives with respect to the 
c         cylindrical coordinates.
c-----------------------------------------------------------------------
          IF (d_order>0) THEN
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                data_r(:,ig,ipol)=dxdr(ig,ipol)*data_x(:,ix,iy)
     $                           +dydr(ig,ipol)*data_y(:,ix,iy)
                data_z(:,ig,ipol)=dxdz(ig,ipol)*data_x(:,ix,iy)
     $                           +dydz(ig,ipol)*data_y(:,ix,iy)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE block_type
c-----------------------------------------------------------------------
c       tblock data:  equilibrium quantities in tblocks are represented
c       by 2D linear elements.
c-----------------------------------------------------------------------
        mx=tg%mcell
        nq=tl%nqty
        ng=inode
        ALLOCATE(dat(nq,mx,1))
        ALLOCATE(data_x(nq,mx,1))
        ALLOCATE(data_y(nq,mx,1))
        DO ig=1,ng
          CALL tri_linear_all_eval(tl,tg,ig,d_order,dat,data_x,data_y)
          data(:,ig,:)=dat(:,:,1)
          data_r(:,ig,:)=data_x(:,:,1)
          data_z(:,ig,:)=data_y(:,:,1)
        ENDDO
      ENDIF block_type
      DEALLOCATE(dat,data_x,data_y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_2D_all_eval
c-----------------------------------------------------------------------
c     subprogram 4. generic_3D_modal_all_eval.
c     get complex dependent-variable field data at the requested
c     quadrature points in each cell in a block.
c
c     this routines evaluates data at the set of points specified by dx 
c     and dy for rblocks and at each of the quadrature-points (1:inode) 
c     for tblocks.
c
c     this modal version is invoked for discontinuous fields that
c     are represented by complete or incomplete modal expansions in
c     rblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_3D_modal_all_eval(modq,tl,dxdr,dydr,dxdz,dydz,
     $                                     dx,dy,tg,inode,data,data_r,
     $                                     data_z,d_order)
      USE modal_disc_mod
      USE tri_linear

      TYPE(modal_quad_type), INTENT(IN) :: modq
      TYPE(tri_linear_type), INTENT(IN) :: tl
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(:), INTENT(IN) :: dx,dy
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: data,data_r,data_z

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: dat,data_x,data_y
      REAL(r8) :: xg,yg
      INTEGER(i4) :: mx,my,iq,nq,ig,ng,ix,iy,ipol,jf,nf
c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        mx=modq%mx
        my=modq%my
        nq=modq%nqty
        nf=modq%nfour
        ng=SIZE(dx)
c-----------------------------------------------------------------------
c       allocate arrays for derivatives.
c-----------------------------------------------------------------------
        ALLOCATE(dat(nq,mx,my,nf))
        IF (d_order>0) THEN
          ALLOCATE(data_x(nq,mx,my,nf))
          ALLOCATE(data_y(nq,mx,my,nf))
        ELSE
          ALLOCATE(data_x(1,1,1,1))
          ALLOCATE(data_y(1,1,1,1))
        ENDIF
c-----------------------------------------------------------------------
c       rblock data:  call the modal interpolation routine at the
c       quadrature points.
c-----------------------------------------------------------------------
        DO ig=1,ng
          xg=dx(ig)
          yg=dy(ig)
          CALL modal_disc_all_eval(modq,xg,yg,dat,data_x,data_y,d_order)
          DO jf=1,nf
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                data(:,ig,ipol,jf)=dat(:,ix,iy,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         use the chain rule to get derivatives with respect to the 
c         cylindrical coordinates.
c-----------------------------------------------------------------------
          IF (d_order>0) THEN
            DO jf=1,nf
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  data_r(:,ig,ipol,jf)=dxdr(ig,ipol)*data_x(:,ix,iy,jf)
     $                                +dydr(ig,ipol)*data_y(:,ix,iy,jf)
                  data_z(:,ig,ipol,jf)=dxdz(ig,ipol)*data_x(:,ix,iy,jf)
     $                                +dydz(ig,ipol)*data_y(:,ix,iy,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE block_type
c-----------------------------------------------------------------------
c-PRE   This is just a placeholder at this point.
c-----------------------------------------------------------------------
      ENDIF block_type
      DEALLOCATE(dat,data_x,data_y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_3D_modal_all_eval
c-----------------------------------------------------------------------
c     subprogram 5. generic_2D_modal_all_eval.
c     get real dependent-variable field data at the requested
c     quadrature points in each cell in a block.
c
c     this routines evaluates data at the set of points specified by dx 
c     and dy for rblocks and at each of the quadrature-points (1:inode) 
c     for tblocks.
c
c     this modal version is invoked for discontinuous fields that
c     are represented by complete or incomplete modal expansions in
c     rblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_2D_modal_all_eval(modq,tl,dxdr,dydr,dxdz,dydz,
     $                                     dx,dy,tg,inode,data,data_r,
     $                                     data_z,d_order)
      USE modal_disc_mod
      USE tri_linear

      TYPE(modal_quad_2D_type), INTENT(IN) :: modq
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(:), INTENT(IN) :: dx,dy
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: data,data_r,data_z

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: dat,data_x,data_y
      REAL(r8) :: xg,yg
      INTEGER(i4) :: mx,my,iq,nq,ig,ng,ix,iy,ipol
c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        mx=modq%mx
        my=modq%my
        nq=modq%nqty
        ng=SIZE(dx)
c-----------------------------------------------------------------------
c       allocate arrays for derivatives.
c-PRE skip extra allocation in all of these.
c-----------------------------------------------------------------------
        ALLOCATE(dat(nq,mx,my))
        IF (d_order>0) THEN
          ALLOCATE(data_x(nq,mx,my))
          ALLOCATE(data_y(nq,mx,my))
        ELSE
          ALLOCATE(data_x(1,1,1))
          ALLOCATE(data_y(1,1,1))
        ENDIF
c-----------------------------------------------------------------------
c       rblock data:  call the modal interpolation routine at the
c       quadrature points.
c-----------------------------------------------------------------------
        DO ig=1,ng
          xg=dx(ig)
          yg=dy(ig)
          CALL modal_disc_all_eval(modq,xg,yg,dat,data_x,data_y,d_order)
          ipol=1
          DO iy=1,my
            DO ix=1,mx
              data(:,ig,ipol)=dat(:,ix,iy)
              ipol=ipol+1
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         use the chain rule to get derivatives with respect to the 
c         cylindrical coordinates.
c-----------------------------------------------------------------------
          IF (d_order>0) THEN
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                data_r(:,ig,ipol)=dxdr(ig,ipol)*data_x(:,ix,iy)
     $                           +dydr(ig,ipol)*data_y(:,ix,iy)
                data_z(:,ig,ipol)=dxdz(ig,ipol)*data_x(:,ix,iy)
     $                           +dydz(ig,ipol)*data_y(:,ix,iy)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE block_type
c-----------------------------------------------------------------------
c-PRE   This is just a placeholder at this point.
c-----------------------------------------------------------------------
      ENDIF block_type
      DEALLOCATE(dat,data_x,data_y)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_2D_modal_all_eval
c-----------------------------------------------------------------------
c     subprogram 6. generic_3D_ptr_set.
c     point to the 3D field data at the requested quadrature point in
c     each cell in a block.
c
c     this routine now just sets the pointer for data at all quadrature
c     points.
c-----------------------------------------------------------------------
      SUBROUTINE generic_3D_ptr_set(qr,qt,tg,inode,data,data_r,
     $                              data_z,d_order)
      USE rblock_type_mod
      USE tblock_type_mod

      TYPE(rb_comp_qp_type), INTENT(IN), TARGET :: qr
      TYPE(tb_comp_qp_type), INTENT(IN), TARGET :: qt
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: data,
     $             data_r,data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        data=>qr%qpf
        IF (d_order>0) THEN
          data_r=>qr%qpfr
          data_z=>qr%qpfz
        ENDIF
c-----------------------------------------------------------------------
c     tblock data:
c-----------------------------------------------------------------------
      ELSE block_type
        data=>qt%qpf
        IF (d_order>0) THEN
          data_r=>qt%qpfr
          data_z=>qt%qpfz
        ENDIF
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_3D_ptr_set
c-----------------------------------------------------------------------
c     subprogram 7. generic_2D_ptr_set.
c     point to the 2D field data at the requested quadrature point in
c     each cell in a block.
c
c     this routine now just sets the pointer for data at all quadrature
c     points.
c-----------------------------------------------------------------------
      SUBROUTINE generic_2D_ptr_set(qr,qt,tg,inode,data,data_r,
     $                              data_z,d_order)
      USE rblock_type_mod
      USE tblock_type_mod

      TYPE(rb_real_qp_type), INTENT(IN), TARGET :: qr
      TYPE(tb_real_qp_type), INTENT(IN), TARGET :: qt
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: data,data_r,
     $          data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
        data=>qr%qpf
        IF (d_order>0) THEN
          data_r=>qr%qpfr
          data_z=>qr%qpfz
        ENDIF
c-----------------------------------------------------------------------
c     tblock data:
c-----------------------------------------------------------------------
      ELSE block_type
        data=>qt%qpf
        IF (d_order>0) THEN
          data_r=>qt%qpfr
          data_z=>qt%qpfz
        ENDIF
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_2D_ptr_set
c-----------------------------------------------------------------------
c     subprogram 8. generic_alpha_eval.
c     get the test functions and their derivatives used for finite
c     element integrations.
c
c     this routine only locates data that has already been evaluated
c     evaluation at all quadrature points.
c-----------------------------------------------------------------------
      SUBROUTINE generic_alpha_eval(rb,tg,inode,alpha_type,alpha,
     $                              dalpdr,dalpdz,d_order,polyd,dalpdrc,
     $                              polydmin,polydmax)
      USE rblock_type_mod
      USE tri_linear

      TYPE(rblock_type), INTENT(IN), TARGET :: rb
      TYPE(tri_linear_geom_type), INTENT(IN), TARGET :: tg
      INTEGER(i4), INTENT(IN) :: inode,d_order,polyd
      INTEGER(i4), INTENT(IN), OPTIONAL :: polydmin,polydmax
      CHARACTER(*), INTENT(IN) :: alpha_type
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: alpha,dalpdr,
     $          dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS, OPTIONAL ::
     $          dalpdrc

      INTEGER(i4) :: iset,pmin,pmax
c-----------------------------------------------------------------------
c     determine if this is a tblock or an rblock based on whether
c     the tg structures has any vertices.
c
c     if the block is an rblock, identify the appropriate basis set
c     for the requested polynomial degree (polyd).
c
c     if alpha_type is 'rhs,' select the basis functions that are
c     already multiplied by the integration weight times the Jacobian.
c     if it is 'mat,' select the bases that are multiplied by the
c     square root of that same product.
c
c     if the prefix 'mod' appears in alpha_type, use the set of
c     modal bases. these bases may be incomplete, and the optional
c     polydmin and polydmax are needed to indicate polynomial range
c     of the limited expansion.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN   !   quadrilaterals
        IF (alpha_type(1:3)=='mod') THEN
          IF (PRESENT(polydmin)) THEN
            pmin=polydmin
          ELSE
            pmin=0
          ENDIF
          IF (PRESENT(polydmax)) THEN
            pmax=polydmax
          ELSE
            pmax=polyd
          ENDIF
          DO iset=1,SIZE(rb%base_modal)
            IF (rb%base_modal(iset)%poly_deg_basis==polyd.AND.
     $          rb%base_modal(iset)%poly_degmin_basis==pmin.AND.
     $          rb%base_modal(iset)%poly_degmax_basis==pmax) EXIT
          ENDDO
        ELSE
          DO iset=1,SIZE(rb%base_pd)
            IF (rb%base_pd(iset)%poly_deg_basis==polyd) EXIT
          ENDDO
        ENDIF

        SELECT CASE(alpha_type)
        CASE('rhs')
          alpha=>rb%base_pd(iset)%alpha
          IF (d_order>0) THEN
            dalpdr=>rb%base_pd(iset)%dalpdr
            dalpdz=>rb%base_pd(iset)%dalpdz
            IF (PRESENT(dalpdrc)) dalpdrc=>rb%base_pd(iset)%dalpdrc
          ENDIF
        CASE('mat')
          alpha=>rb%base_pd(iset)%alpham
          IF (d_order>0) THEN
            dalpdr=>rb%base_pd(iset)%dalpmdr
            dalpdz=>rb%base_pd(iset)%dalpmdz
            IF (PRESENT(dalpdrc)) dalpdrc=>rb%base_pd(iset)%dalpmdrc
          ENDIF
        CASE('modlrhs')
          alpha=>rb%base_modal(iset)%alpha
          IF (d_order>0) THEN
            dalpdr=>rb%base_modal(iset)%dalpdr
            dalpdz=>rb%base_modal(iset)%dalpdz
            IF (PRESENT(dalpdrc)) dalpdrc=>rb%base_modal(iset)%dalpdrc
          ENDIF
        CASE DEFAULT  !  'modlmat'
          alpha=>rb%base_modal(iset)%alpham
          IF (d_order>0) THEN
            dalpdr=>rb%base_modal(iset)%dalpmdr
            dalpdz=>rb%base_modal(iset)%dalpmdz
            IF (PRESENT(dalpdrc)) dalpdrc=>rb%base_modal(iset)%dalpmdrc
          ENDIF
        END SELECT
      ELSE block_type  !   triangles
        SELECT CASE(alpha_type)
        CASE('rhs')
          alpha=>tg%alpha_arr
          IF (d_order>0) THEN
            dalpdr=>tg%dalpdr
            dalpdz=>tg%dalpdz
            IF (PRESENT(dalpdrc)) dalpdrc=>tg%dalpdrc
          ENDIF
        CASE DEFAULT  !  'mat' or 'modlmat'
          alpha=>tg%alpham_arr
          IF (d_order>0) THEN
            dalpdr=>tg%dalpmdr
            dalpdz=>tg%dalpmdz
            IF (PRESENT(dalpdrc)) dalpdrc=>tg%dalpmdrc
          ENDIF
        END SELECT
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_alpha_eval
c-----------------------------------------------------------------------
c     subprogram 9. generic_3D_eval.
c     get complex dependent-variable field data at the requested
c     point in a block.
c
c     this routine performs a function evaluation at the specified
c     point.
c-----------------------------------------------------------------------
      SUBROUTINE generic_3D_eval(laq,tl,dxdr,dydr,dxdz,dydz,x,y,
     $                           tg,ijcell,data,data_r,data_z,d_order)
      USE lagr_quad_mod
      USE tri_linear

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), INTENT(IN) :: x,y
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: d_order
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: ijcell
      COMPLEX(r8), DIMENSION(:,:), INTENT(OUT) :: data,data_r,data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
c-----------------------------------------------------------------------
c       rblock data:  call the lagr_quad interpolation routine.
c-----------------------------------------------------------------------
        CALL lagr_quad_eval(laq,x,y,d_order)
        data=laq%f
c-----------------------------------------------------------------------
c       use the chain rule to get derivatives with respect to the 
c       cylindrical coordinates.
c-----------------------------------------------------------------------
        IF (d_order>0) THEN
          data_r=dxdr*laq%fx+dydr*laq%fy
          data_z=dxdz*laq%fx+dydz*laq%fy
        ENDIF
      ELSE block_type
c-----------------------------------------------------------------------
c       tblock data:  all quantities in tblocks are represented
c       by linear elements with Fourier series in the 3rd direction.
c       linear elements are interpolated directly in the physical
c       coordinates.
c-----------------------------------------------------------------------
        CALL tri_linear_eval(tl,tg,x,y,ijcell(1),d_order)
        data=tl%f
        IF (d_order>0) THEN
          data_r=tl%fx
          data_z=tl%fy
        ENDIF
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_3D_eval
c-----------------------------------------------------------------------
c     subprogram 10. generic_2D_eval.
c     get real-valued dependent-variable field data at the requested
c     point in a block.
c
c     this routine performs a function evaluation at the specified
c     point.
c-----------------------------------------------------------------------
      SUBROUTINE generic_2D_eval(laq,tl,dxdr,dydr,dxdz,dydz,x,y,
     $                           tg,ijcell,data,data_r,data_z,d_order)
      USE lagr_quad_mod
      USE tri_linear

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), INTENT(IN) :: x,y
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: d_order
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: ijcell
      REAL(r8), DIMENSION(:), INTENT(OUT) :: data,data_r,data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
c-----------------------------------------------------------------------
c       rblock data:  call the lagr_quad interpolation routine.
c-----------------------------------------------------------------------
        CALL lagr_quad_eval(laq,x,y,d_order)
        data=laq%f
c-----------------------------------------------------------------------
c       use the chain rule to get derivatives with respect to the 
c       cylindrical coordinates.
c-----------------------------------------------------------------------
        IF (d_order>0) THEN
          data_r=dxdr*laq%fx+dydr*laq%fy
          data_z=dxdz*laq%fx+dydz*laq%fy
        ENDIF
      ELSE block_type
c-----------------------------------------------------------------------
c       tblock data:  all quantities in tblocks are represented
c       by linear elements with Fourier series in the 3rd direction.
c       linear elements are interpolated directly in the physical
c       coordinates.
c-----------------------------------------------------------------------
        CALL tri_linear_eval(tl,tg,x,y,ijcell(1),d_order)
        data=tl%f
        IF (d_order>0) THEN
          data_r=tl%fx
          data_z=tl%fy
        ENDIF
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_2D_eval
c-----------------------------------------------------------------------
c     subprogram 11. generic_3D_modal_eval.
c     get complex dependent-variable field data at the requested
c     point in a block.
c
c     this modal version is invoked for discontinuous fields that
c     are represented by complete or incomplete modal expansions in
c     rblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_3D_modal_eval(modq,tl,dxdr,dydr,dxdz,dydz,x,y,
     $                                 tg,ijcell,data,data_r,data_z,
     $                                 d_order)
      USE modal_disc_mod
      USE tri_linear

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), INTENT(IN) :: x,y
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: d_order
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: ijcell
      COMPLEX(r8), DIMENSION(:,:), INTENT(OUT) :: data,data_r,data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
c-----------------------------------------------------------------------
c       rblock data:  call the modal interpolation routine.
c-----------------------------------------------------------------------
        CALL modal_disc_eval(modq,x,y,d_order)
        data=modq%f
c-----------------------------------------------------------------------
c       use the chain rule to get derivatives with respect to the 
c       cylindrical coordinates.
c-----------------------------------------------------------------------
        IF (d_order>0) THEN
          data_r=dxdr*modq%fx+dydr*modq%fy
          data_z=dxdz*modq%fx+dydz*modq%fy
        ENDIF
      ELSE block_type
c-----------------------------------------------------------------------
c-PRE   This is just a placeholder at this point.
c-----------------------------------------------------------------------
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_3D_modal_eval
c-----------------------------------------------------------------------
c     subprogram 12. generic_2D_modal_eval.
c     get real dependent-variable field data at the requested
c     point in a block.
c
c     this modal version is invoked for discontinuous fields that
c     are represented by complete or incomplete modal expansions in
c     rblocks.
c-----------------------------------------------------------------------
      SUBROUTINE generic_2D_modal_eval(modq,tl,dxdr,dydr,dxdz,dydz,x,y,
     $                                 tg,ijcell,data,data_r,data_z,
     $                                 d_order)
      USE modal_disc_mod
      USE tri_linear

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: dxdr,dydr,dxdz,dydz
      REAL(r8), INTENT(IN) :: x,y
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: d_order
      INTEGER(i4), DIMENSION(2), INTENT(IN) :: ijcell
      REAL(r8), DIMENSION(:), INTENT(OUT) :: data,data_r,data_z

c-----------------------------------------------------------------------
c     decide which block is the dummy.
c-----------------------------------------------------------------------
      block_type: IF (tg%mvert<0) THEN
c-----------------------------------------------------------------------
c       rblock data:  call the modal interpolation routine.
c-----------------------------------------------------------------------
        CALL modal_disc_eval(modq,x,y,d_order)
        data=modq%f
c-----------------------------------------------------------------------
c       use the chain rule to get derivatives with respect to the 
c       cylindrical coordinates.
c-----------------------------------------------------------------------
        IF (d_order>0) THEN
          data_r=dxdr*modq%fx+dydr*modq%fy
          data_z=dxdz*modq%fx+dydz*modq%fy
        ENDIF
      ELSE block_type
c-----------------------------------------------------------------------
c-PRE   This is just a placeholder at this point.
c-----------------------------------------------------------------------
      ENDIF block_type
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE generic_2D_modal_eval
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE generic_evals
