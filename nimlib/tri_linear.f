c-----------------------------------------------------------------------
c     file tri_linear.f
c     fits functions to tri_linear splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  tri_linear_type definition.
c     1.  tri_linear_geom_alloc.
c     2.  tri_linear_geom_dealloc.
c     3.  tri_linear_get_areas.
c     4.  tri_linear_3D_alloc.
c     5.  tri_linear_3D_dealloc.
c     6.  tri_linear_3D_all_eval.
c     7.  tri_linear_3D_eval.
c     8.  tri_linear_2D_alloc.
c     9.  tri_linear_2D_dealloc.
c     10. tri_linear_2D_all_eval.
c     11. tri_linear_2D_eval.
c     12. tri_linear_3D_assign_rsc.
c     13. tri_linear_3D_assign_csc.
c     14. tri_linear_3D_assign_tl3.
c     15. tri_linear_3D_assign_int.
c     16. tri_linear_2D_assign_rsc.
c     17. tri_linear_2D_assign_csc.
c     18. tri_linear_2D_assign_tl2.
c     19. tri_linear_2D_assign_int.
c-----------------------------------------------------------------------
c     subprogram 0. tri_linear_type definition.
c     defines tri_linear_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE tri_linear
      USE local
      IMPLICIT NONE
      
      TYPE :: neighbor_type
      INTEGER(i4), DIMENSION(:), POINTER :: vertex
      END TYPE neighbor_type
      
      TYPE :: tri_linear_geom_type
        INTEGER(i4) :: mvert,mcell
        INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: vertex
        INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: segment
        REAL(r8), DIMENSION(:), ALLOCATABLE :: xs,ys,area
        REAL(r8), DIMENSION(3,1:7) :: alpha,alphab,delta
        REAL(r8), DIMENSION(:,:), ALLOCATABLE :: wjac,bigr
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE ::
     $    alpha_x,alpha_y,alpha_arr,dalpdr,dalpdz,dalpdrc,
     $    dalpmdr,dalpmdz,alpham_arr,dalpmdrc
        TYPE(neighbor_type), DIMENSION(:), POINTER :: neighbor
      END TYPE tri_linear_geom_type
      
      TYPE :: tri_linear_type
        CHARACTER(6) :: name
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        INTEGER(i4) :: mvert,nqty,nfour
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fs
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: f,fx,fy
      END TYPE tri_linear_type

      TYPE :: tri_linear_2D_type
        CHARACTER(6) :: name
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        INTEGER(i4) :: mvert,nqty
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs
        REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy
      END TYPE tri_linear_2D_type

c-----------------------------------------------------------------------
c     subprogram name interfaces
c-----------------------------------------------------------------------
      INTERFACE tri_linear_alloc
        MODULE PROCEDURE tri_linear_2D_alloc,tri_linear_3D_alloc
      END INTERFACE

      INTERFACE tri_linear_dealloc
        MODULE PROCEDURE tri_linear_2D_dealloc,tri_linear_3D_dealloc
      END INTERFACE

      INTERFACE tri_linear_all_eval
        MODULE PROCEDURE tri_linear_2D_all_eval,tri_linear_3D_all_eval
      END INTERFACE

      INTERFACE tri_linear_eval
        MODULE PROCEDURE tri_linear_2D_eval,tri_linear_3D_eval
      END INTERFACE
c-----------------------------------------------------------------------
c     overloaded assignment.
c-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE tri_linear_3D_assign_csc,
     $    tri_linear_3D_assign_rsc,tri_linear_3D_assign_tl3,
     $    tri_linear_3D_assign_int,tri_linear_2D_assign_csc,
     $    tri_linear_2D_assign_rsc,tri_linear_2D_assign_tl2,
     $    tri_linear_2D_assign_int
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. tri_linear_geom_alloc.
c     allocates space for tri_linear geometry.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_geom_alloc(v,mvert,mcell)
      
      INTEGER(i4), INTENT(IN) :: mvert,mcell
      TYPE(tri_linear_geom_type), INTENT(OUT) :: v
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      v%mvert=mvert
      v%mcell=mcell
      ALLOCATE(v%xs(0:mvert))
      ALLOCATE(v%ys(0:mvert))
      ALLOCATE(v%vertex(mcell,3))
      ALLOCATE(v%alpha_x(mcell,1,3))
      ALLOCATE(v%alpha_y(mcell,1,3))
      ALLOCATE(v%neighbor(0:mvert))
      ALLOCATE(v%area(mcell))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_geom_alloc
c-----------------------------------------------------------------------
c     subprogram 2. tri_linear_geom_dealloc.
c     deallocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_geom_dealloc(v)
      
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: v
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(v%xs)
      DEALLOCATE(v%ys)
      DEALLOCATE(v%vertex)
      DEALLOCATE(v%alpha_x)
      DEALLOCATE(v%alpha_y)
      DEALLOCATE(v%neighbor)
      DEALLOCATE(v%area)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_geom_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. tri_linear_get_areas.
c     evaluates areas and basis functions of triangular cells.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_get_areas(v)
      
      TYPE(tri_linear_geom_type), INTENT(INOUT) :: v
      
      INTEGER(i4) :: icell,iqty
      REAL(r8), PARAMETER ::
     $     alpha1=0.05971587178976982_r8,
     $     alpha2=0.1012865073234563_r8,
     $     alpha3=0.3333333333333333_r8,
     $     alpha4=0.4701420641051151_r8,
     $     alpha5=0.7974269853530873_r8
c-----------------------------------------------------------------------
c     set weight functions.
c-----------------------------------------------------------------------
      v%alpha=TRANSPOSE(RESHAPE((/
     $     alpha3,alpha1,alpha4,alpha4,alpha5,alpha2,alpha2,
     $     alpha3,alpha4,alpha1,alpha4,alpha2,alpha5,alpha2,
     $     alpha3,alpha4,alpha4,alpha1,alpha2,alpha2,alpha5
     $     /),(/7,3/)))
      v%alphab=alpha3
      v%delta=TRANSPOSE(RESHAPE((/
     $     1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),(/7,3/)))
c-----------------------------------------------------------------------
c     compute areas.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         v%area(icell)=0.5_r8*
     $        (v%xs(v%vertex(icell,1))*v%ys(v%vertex(icell,2))
     $        -v%xs(v%vertex(icell,2))*v%ys(v%vertex(icell,1))
     $        +v%xs(v%vertex(icell,2))*v%ys(v%vertex(icell,3))
     $        -v%xs(v%vertex(icell,3))*v%ys(v%vertex(icell,2))
     $        +v%xs(v%vertex(icell,3))*v%ys(v%vertex(icell,1))
     $        -v%xs(v%vertex(icell,1))*v%ys(v%vertex(icell,3)))
      ENDDO
c-----------------------------------------------------------------------
c     compute factors for x derivatives.  the second array index is to
c     match the basis array structure in quadrilateral blocks.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         v%alpha_x(icell,1,1)=0.5_r8/v%area(icell)
     $        *(v%ys(v%vertex(icell,2))-v%ys(v%vertex(icell,3)))
         v%alpha_x(icell,1,2)=0.5_r8/v%area(icell)
     $        *(v%ys(v%vertex(icell,3))-v%ys(v%vertex(icell,1)))
         v%alpha_x(icell,1,3)=0.5_r8/v%area(icell)
     $        *(v%ys(v%vertex(icell,1))-v%ys(v%vertex(icell,2)))
      ENDDO
c-----------------------------------------------------------------------
c     compute factors for y derivatives.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         v%alpha_y(icell,1,1)=-0.5_r8/v%area(icell)
     $        *(v%xs(v%vertex(icell,2))-v%xs(v%vertex(icell,3)))
         v%alpha_y(icell,1,2)=-0.5_r8/v%area(icell)
     $        *(v%xs(v%vertex(icell,3))-v%xs(v%vertex(icell,1)))
         v%alpha_y(icell,1,3)=-0.5_r8/v%area(icell)
     $        *(v%xs(v%vertex(icell,1))-v%xs(v%vertex(icell,2)))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_get_areas
c-----------------------------------------------------------------------
c     subprogram 4. tri_linear_3D_alloc.
c     allocates space for tri_linear dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_alloc(u,mvert,nqty,nfour,name,title)
      
      INTEGER(i4), INTENT(IN) :: mvert,nqty,nfour
      TYPE(tri_linear_type), INTENT(OUT) :: u
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      u%mvert=mvert
      u%nqty=nqty
      u%nfour=nfour
      ALLOCATE(u%fs(nqty,0:mvert,0:0,nfour))
      ALLOCATE(u%title(nqty))
      ALLOCATE(u%f(nqty,nfour))
      ALLOCATE(u%fx(nqty,nfour))
      ALLOCATE(u%fy(nqty,nfour))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) u%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          u%title=title
        ELSE
          u%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_alloc
c-----------------------------------------------------------------------
c     subprogram 5. tri_linear_3D_dealloc.
c     deallocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_dealloc(u)

      TYPE(tri_linear_type), INTENT(INOUT) :: u
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(u%fs)
      DEALLOCATE(u%title)
      DEALLOCATE(u%f)
      DEALLOCATE(u%fx)
      DEALLOCATE(u%fy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_dealloc
c-----------------------------------------------------------------------
c     subprogram 6. tri_linear_3D_all_eval.
c     evaluates bicubic splines in all intervals for equal spacing. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_all_eval(u,v,node,mode,f,fx,fy)

      TYPE(tri_linear_type), INTENT(IN) :: u
      TYPE(tri_linear_geom_type), INTENT(IN) :: v
      INTEGER(i4), INTENT(IN) :: node,mode
      COMPLEX(r8), INTENT(OUT),DIMENSION(:,:,0:,:) :: f,fx,fy

      INTEGER(i4) :: icell,iqty
c-----------------------------------------------------------------------
c     verify consistency between u and v.
c-----------------------------------------------------------------------
      IF(u%mvert.NE.v%mvert)THEN
         WRITE(nim_wr,'(a)')
     $      "tri_linear_all_eval: geometry and dependent",
     $      " variables have inconsistent sizes."
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     interpolate functions.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         f(:,icell,0,:)
     $        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha(1,node)
     $        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha(2,node)
     $        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha(3,node)
      ENDDO
      IF(mode.EQ.0)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         fx(:,icell,0,:)
     $        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha_x(icell,1,1)
     $        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha_x(icell,1,2)
     $        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha_x(icell,1,3)
         fy(:,icell,0,:)
     $        =u%fs(:,v%vertex(icell,1),0,:)*v%alpha_y(icell,1,1)
     $        +u%fs(:,v%vertex(icell,2),0,:)*v%alpha_y(icell,1,2)
     $        +u%fs(:,v%vertex(icell,3),0,:)*v%alpha_y(icell,1,3)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_all_eval
c-----------------------------------------------------------------------
c     subprogram 7. tri_linear_3D_eval.
c     linear interpolation at (x,y) of u base on the
c     three points associated with cell icell on grid v.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_eval(u,v,x,y,icell,mode)

      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: icell,mode
      TYPE(tri_linear_type), INTENT(INOUT) :: u
      TYPE(tri_linear_geom_type), INTENT(IN) :: v

      INTEGER(i4) :: iqty,ifour
      TYPE :: position_type
        REAL(r8) :: r,z
      END TYPE position_type
      TYPE(position_type) :: p,p1,p2,p3
      REAL(r8) :: denom
      COMPLEX(r8) :: fx,fy

c-----------------------------------------------------------------------
c     verify consistency between u and v.
c-----------------------------------------------------------------------
      IF(u%mvert.NE.v%mvert)THEN
         WRITE(nim_wr,'(a)')
     $     "tri_linear_eval: geometry and dependent",
     $     " variables have inconsistent sizes."
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     verify icell is a valid cell
c-----------------------------------------------------------------------
      IF((icell < 0). OR. (icell > v%mcell))THEN
         CALL nim_stop("tri_linear_eval: invalid cell")
      ENDIF
c-----------------------------------------------------------------------
c     Define positions
c-----------------------------------------------------------------------
      p%r=x
      p%z=y
      p1%r=v%xs(v%vertex(icell,1))
      p1%z=v%ys(v%vertex(icell,1))
      p2%r=v%xs(v%vertex(icell,2))
      p2%z=v%ys(v%vertex(icell,2))
      p3%r=v%xs(v%vertex(icell,3))
      p3%z=v%ys(v%vertex(icell,3))
      denom=1.0_r8/((p2%r-p1%r)*(p3%z-p1%z)-(p3%r-p1%r)*(p2%z-p1%z))
c-----------------------------------------------------------------------
c     Evaluate derivatives and functions.
c-----------------------------------------------------------------------
      DO ifour=1,u%nfour
        DO iqty=1,u%nqty
          fx=denom*(
     &            (u%fs(iqty,v%vertex(icell,2),0,ifour)-
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p3%z-p1%z)-
     &            (u%fs(iqty,v%vertex(icell,3),0,ifour)-
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p2%z-p1%z))
          fy=denom*(
     &            (u%fs(iqty,v%vertex(icell,3),0,ifour)-
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p2%r-p1%r)-
     &            (u%fs(iqty,v%vertex(icell,2),0,ifour)-
     &             u%fs(iqty,v%vertex(icell,1),0,ifour))*(p3%r-p1%r))
          u%f(iqty,ifour)=
     &            u%fs(iqty,v%vertex(icell,1),0,ifour)+
     &            fy*(p%z-p1%z)+ 
     &            fx*(p%r-p1%r)
          IF(mode /= 0 )THEN
            u%fx(iqty,ifour) = fx
            u%fy(iqty,ifour) = fy
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE tri_linear_3D_eval
c-----------------------------------------------------------------------
c     subprogram 8. tri_linear_2D_alloc.
c     allocates space for tri_linear dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_alloc(u,mvert,nqty,name,title)
      
      INTEGER(i4), INTENT(IN) :: mvert,nqty
      TYPE(tri_linear_2D_type), INTENT(OUT) :: u
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      u%mvert=mvert
      u%nqty=nqty
      ALLOCATE(u%fs(nqty,0:mvert,0:0))
      ALLOCATE(u%title(nqty))
      ALLOCATE(u%f(nqty))
      ALLOCATE(u%fx(nqty))
      ALLOCATE(u%fy(nqty))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) u%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          u%title=title
        ELSE
          u%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_alloc
c-----------------------------------------------------------------------
c     subprogram 9. tri_linear_2D_dealloc.
c     deallocates space for spline_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_dealloc(u)

      TYPE(tri_linear_2D_type), INTENT(INOUT) :: u
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(u%fs)
      DEALLOCATE(u%title)
      DEALLOCATE(u%f)
      DEALLOCATE(u%fx)
      DEALLOCATE(u%fy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_dealloc
c-----------------------------------------------------------------------
c     subprogram 10. tri_linear_2D_all_eval.
c     evaluates bicubic splines in all intervals for equal spacing. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_all_eval(u,v,node,mode,f,fx,fy)

      TYPE(tri_linear_2D_type), INTENT(IN) :: u
      TYPE(tri_linear_geom_type), INTENT(IN) :: v
      INTEGER(i4), INTENT(IN) :: node,mode
      REAL(r8), INTENT(OUT),DIMENSION(:,:,0:) :: f,fx,fy

      INTEGER(i4) :: icell,iqty
c-----------------------------------------------------------------------
c     verify consistency between u and v.
c-----------------------------------------------------------------------
      IF(u%mvert.NE.v%mvert)THEN
         WRITE(nim_wr,'(a)')
     $     "tri_linear_all_eval: geometry and dependent",
     $     " variables have inconsistent sizes."
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     interpolate functions.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         f(:,icell,0)
     $        =u%fs(:,v%vertex(icell,1),0)*v%alpha(1,node)
     $        +u%fs(:,v%vertex(icell,2),0)*v%alpha(2,node)
     $        +u%fs(:,v%vertex(icell,3),0)*v%alpha(3,node)
      ENDDO
      IF(mode.EQ.0)RETURN
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      DO icell=1,v%mcell
         fx(:,icell,0)
     $        =u%fs(:,v%vertex(icell,1),0)*v%alpha_x(icell,1,1)
     $        +u%fs(:,v%vertex(icell,2),0)*v%alpha_x(icell,1,2)
     $        +u%fs(:,v%vertex(icell,3),0)*v%alpha_x(icell,1,3)
         fy(:,icell,0)
     $        =u%fs(:,v%vertex(icell,1),0)*v%alpha_y(icell,1,1)
     $        +u%fs(:,v%vertex(icell,2),0)*v%alpha_y(icell,1,2)
     $        +u%fs(:,v%vertex(icell,3),0)*v%alpha_y(icell,1,3)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_all_eval
c-----------------------------------------------------------------------
c     subprogram 11. tri_linear_2D_eval.
c     linear interpolation at (x,y) of u base on the
c     three points associated with cell icell on grid v.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_eval(u,v,x,y,icell,mode)

      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: icell,mode
      TYPE(tri_linear_2D_type), INTENT(INOUT) :: u
      TYPE(tri_linear_geom_type), INTENT(IN) :: v

      INTEGER(i4) :: iqty
      TYPE :: position_type
        REAL(r8) :: r,z
      END TYPE position_type
      TYPE(position_type) :: p,p1,p2,p3
      REAL(r8) :: denom,fx,fy

c-----------------------------------------------------------------------
c     verify consistency between u and v.
c-----------------------------------------------------------------------
      IF(u%mvert.NE.v%mvert)THEN
         WRITE(nim_wr,'(a)')
     $     "tri_linear_eval: geometry and dependent",
     $     " variables have inconsistent sizes."
         STOP
      ENDIF
c-----------------------------------------------------------------------
c     verify icell is a valid cell
c-----------------------------------------------------------------------
      IF((icell < 0). OR. (icell > v%mcell))THEN
         CALL nim_stop("tri_linear_eval: invalid cell")
      ENDIF
c-----------------------------------------------------------------------
c     Define positions
c-----------------------------------------------------------------------
      p%r=x
      p%z=y
      p1%r=v%xs(v%vertex(icell,1))
      p1%z=v%ys(v%vertex(icell,1))
      p2%r=v%xs(v%vertex(icell,2))
      p2%z=v%ys(v%vertex(icell,2))
      p3%r=v%xs(v%vertex(icell,3))
      p3%z=v%ys(v%vertex(icell,3))
      denom=1.0_r8/((p2%r-p1%r)*(p3%z-p1%z)-(p3%r-p1%r)*(p2%z-p1%z))
c-----------------------------------------------------------------------
c     Evaluate derivatives and functions.
c-----------------------------------------------------------------------
      DO iqty=1,u%nqty
        fx=denom*(
     &          (u%fs(iqty,v%vertex(icell,2),0)-
     &           u%fs(iqty,v%vertex(icell,1),0))*(p3%z-p1%z)-
     &          (u%fs(iqty,v%vertex(icell,3),0)-
     &           u%fs(iqty,v%vertex(icell,1),0))*(p2%z-p1%z))
        fy=denom*(
     &          (u%fs(iqty,v%vertex(icell,3),0)-
     &           u%fs(iqty,v%vertex(icell,1),0))*(p2%r-p1%r)-
     &          (u%fs(iqty,v%vertex(icell,2),0)-
     &           u%fs(iqty,v%vertex(icell,1),0))*(p3%r-p1%r))
        u%f(iqty)=
     &          u%fs(iqty,v%vertex(icell,1),0)+
     &          fy*(p%z-p1%z)+ 
     &          fx*(p%r-p1%r)
        IF(mode /= 0 )THEN
          u%fx(iqty) = fx
          u%fy(iqty) = fy
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE tri_linear_2D_eval
c-----------------------------------------------------------------------
c     subprogram 12. tri_linear_3D_assign_rsc.
c     assign a real scalar value to a 3D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_rsc(tl,rscalar)

      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=rscalar
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 13. tri_linear_3D_assign_csc.
c     assign a complex scalar value to a 3D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_csc(tl,cscalar)

      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=cscalar
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 14. tri_linear_3D_assign_tl3.
c     set one 3D tri_linear structure equal to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_tl3(tl1,tl2)

      TYPE(tri_linear_type), INTENT(INOUT) :: tl1
      TYPE(tri_linear_type), INTENT(IN) :: tl2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      tl1%fs=tl2%fs
c-PRE IF (ASSOCIATED(tl1%fss)) THEN
c       tl1%fss=tl2%fss
c       tl1%fsi=tl2%fsi
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_assign_tl3
c-----------------------------------------------------------------------
c     subprogram 15. tri_linear_3D_assign_int.
c     assign a integer value to a 3D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_3D_assign_int(tl,int)

      TYPE(tri_linear_type), INTENT(INOUT) :: tl
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=int
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_3D_assign_int
c-----------------------------------------------------------------------
c     subprogram 12. tri_linear_2D_assign_rsc.
c     assign a real scalar value to a 2D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_rsc(tl,rscalar)

      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=rscalar
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 12. tri_linear_2D_assign_csc.
c     assign a complex scalar value to a 2D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_csc(tl,cscalar)

      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=cscalar
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 14. tri_linear_2D_assign_tl2.
c     set one 2D tri_linear structure equal to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_tl2(tl1,tl2)

      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl1
      TYPE(tri_linear_2D_type), INTENT(IN) :: tl2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      tl1%fs=tl2%fs
c-PRE IF (ASSOCIATED(tl1%fss)) THEN
c       tl1%fss=tl2%fss
c       tl1%fsi=tl2%fsi
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_assign_tl2
c-----------------------------------------------------------------------
c     subprogram 15. tri_linear_2D_assign_int.
c     assign a integer value to a 2D tri_linear structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE tri_linear_2D_assign_int(tl,int)

      TYPE(tri_linear_2D_type), INTENT(INOUT) :: tl
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     grid vertex nodes completely represent the scalar.
c-----------------------------------------------------------------------
      tl%fs=int
c-PRE IF (ASSOCIATED(tl%fss)) THEN
c       tl%fss=0
c       tl%fsi=0
c     ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tri_linear_2D_assign_int
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE tri_linear
