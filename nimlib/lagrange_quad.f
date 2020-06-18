c-----------------------------------------------------------------------
c     file lagrange_quad.f
c     routines for evaluating Lagrange finite elements on blocks of
c     structured quadrilaterals.  this package is based on the 
c     previous bilinear-only version.
c
c     throughout this package, basis functions are catagorized into node
c     (grid vertex), side, and interior for each element.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0.  lagr_type_mod.
c      0.1 lagr_quad_mod.
c      1.  lagr_quad_bases.
c      2.  lagr_quad_3D_alloc.
c      3.  lagr_quad_3D_dealloc.
c      4.  lagr_quad_3D_eval.
c      5.  lagr_quad_3D_all_eval.
c      6.  lagr_quad_3D_assign_rsc.
c      7.  lagr_quad_3D_assign_csc.
c      8.  lagr_quad_3D_assign_laq.
c      9.  lagr_quad_3D_assign_int.
c      10. lagr_quad_3D_basis_assign_arr
c      11. lagr_quad_3D_basis_add_arr
c      12. lagr_quad_3D_basis_assign_loc
c      13. lagr_quad_3D_basis_add_loc
c      14. lagr_quad_2D_alloc.
c      15. lagr_quad_2D_dealloc.
c      16. lagr_quad_2D_eval.
c      17. lagr_quad_2D_all_eval.
c      18. lagr_quad_2D_assign_rsc.
c      19. lagr_quad_2D_assign_csc.
c      20. lagr_quad_2D_assign_laq.
c      21. lagr_quad_2D_assign_int.
c      22. lagr_quad_2D_basis_assign_arr
c      23. lagr_quad_2D_basis_add_arr
c      24. lagr_quad_2D_basis_assign_loc
c      25. lagr_quad_2D_basis_add_loc
c-----------------------------------------------------------------------
c     subprogram 0. lagr_type_mod definition.
c     defines the lagrangian quadrilateral types and may be used
c     separately from lagr_quad_mod.
c-----------------------------------------------------------------------
      MODULE lagr_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: lagr_quad_type
        INTEGER(i4) :: mx,my,nqty,nfour,n_side,n_int
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx,dy
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fs
        COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: fsh,fsv,fsi
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        CHARACTER(6) :: name
      END TYPE lagr_quad_type

      TYPE :: lagr_quad_2D_type
        INTEGER(i4) :: mx,my,nqty,n_side,n_int
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        REAL(r8), DIMENSION(:), ALLOCATABLE :: dx,dy
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs
        REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fsh,fsv,fsi
        REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        CHARACTER(6) :: name
      END TYPE lagr_quad_2D_type

      END MODULE lagr_type_mod
c-----------------------------------------------------------------------
c     subprogram 0.1 lagr_quad_mod definition.
c     contains the subprograms for manipulating lagrange_quad expansions
c     in structured blocks of quadrlaters, where the bases are
c     continuous across element borders.
c-----------------------------------------------------------------------
      MODULE lagr_quad_mod
      USE lagr_type_mod
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20
      REAL(r8), DIMENSION(0:npoly_max), PRIVATE :: alx,aly,dalx,daly
c-----------------------------------------------------------------------
c     subprogram name interfaces
c-----------------------------------------------------------------------
      INTERFACE lagr_quad_alloc
        MODULE PROCEDURE lagr_quad_2D_alloc,lagr_quad_3D_alloc
      END INTERFACE

      INTERFACE lagr_quad_dealloc
        MODULE PROCEDURE lagr_quad_2D_dealloc,lagr_quad_3D_dealloc
      END INTERFACE

      INTERFACE lagr_quad_all_eval
        MODULE PROCEDURE lagr_quad_2D_all_eval,lagr_quad_3D_all_eval
      END INTERFACE

      INTERFACE lagr_quad_eval
        MODULE PROCEDURE lagr_quad_2D_eval,lagr_quad_3D_eval
      END INTERFACE

      INTERFACE lagr_quad_basis_assign_arr
        MODULE PROCEDURE
     $    lagr_quad_2D_basis_assign_arr,lagr_quad_3D_basis_assign_arr
      END INTERFACE

      INTERFACE lagr_quad_basis_add_arr
        MODULE PROCEDURE
     $    lagr_quad_2D_basis_add_arr,lagr_quad_3D_basis_add_arr
      END INTERFACE

      INTERFACE lagr_quad_basis_assign_loc
        MODULE PROCEDURE
     $    lagr_quad_2D_basis_assign_loc,lagr_quad_3D_basis_assign_loc
      END INTERFACE

      INTERFACE lagr_quad_basis_add_loc
        MODULE PROCEDURE
     $    lagr_quad_2D_basis_add_loc,lagr_quad_3D_basis_add_loc
      END INTERFACE
c-----------------------------------------------------------------------
c     overloaded assignment.
c-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE
     $    lagr_quad_3D_assign_csc,lagr_quad_3D_assign_rsc,
     $    lagr_quad_3D_assign_laq,lagr_quad_3D_assign_int,
     $    lagr_quad_2D_assign_csc,lagr_quad_2D_assign_rsc,
     $    lagr_quad_2D_assign_laq,lagr_quad_2D_assign_int
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. lagr_quad_bases.
c     computes basis functions and their derivatives.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_bases(x,y,alpha,alphax,alphay,dmode)

      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax,alphay
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: pd,i,j,k
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
c     first determine the polynomial degree and set 1D basis functions.
c-----------------------------------------------------------------------
      pd=NINT(SQRT(REAL(SIZE(alpha))))-1
      CALL lagr_1D(pd,x,alx,dalx,dmode)
      CALL lagr_1D(pd,y,aly,daly,dmode)
c-----------------------------------------------------------------------
c     if derivatives are not needed, limit the computations.
c-----------------------------------------------------------------------
      IF (dmode==0) THEN
c-----------------------------------------------------------------------
c       compute basis functions with respect to the logical coordinates
c       of the unit square.  node (grid vertex-centered) bases first;
c       start in lower left corner, work left to right, bottom to top.
c-----------------------------------------------------------------------
        alpha(1)=alx(0) *aly(0)
        alpha(2)=alx(pd)*aly(0)
        alpha(3)=alx(0) *aly(pd)
        alpha(4)=alx(pd)*aly(pd)
c-----------------------------------------------------------------------
c       horizontal side centered bases if needed.  the ordering is
c       bottom to top, then left to right in pairs.
c-----------------------------------------------------------------------
        k=5
        DO i=1,pd-1
          alpha(k)= alx(i)* aly(0) 
          k=k+1
          alpha(k)= alx(i)* aly(pd) 
          k=k+1
        ENDDO
c-----------------------------------------------------------------------
c       vertical side centered bases if needed.  the ordering is
c       left to right, then bottom to top in pairs.
c-----------------------------------------------------------------------
        DO j=1,pd-1
          alpha(k)= alx(0)* aly(j) 
          k=k+1
          alpha(k)= alx(pd)* aly(j) 
          k=k+1
        ENDDO
c-----------------------------------------------------------------------
c       interior bases if needed.  the ordering is
c       left to right, bottom to top.
c-----------------------------------------------------------------------
        DO j=1,pd-1
          DO i=1,pd-1
            alpha(k)= alx(i)* aly(j) 
            k=k+1
          ENDDO
        ENDDO
      ELSE
c-----------------------------------------------------------------------
c       compute basis functions and derivatives with respect to the
c       logical coordinates of the unit square.  node (grid vertex-
c       centered) bases first; start in lower left corner, work left to
c       right, bottom to top.
c-----------------------------------------------------------------------
        alpha(1)=alx(0) *aly(0)
        alpha(2)=alx(pd)*aly(0)
        alpha(3)=alx(0) *aly(pd)
        alpha(4)=alx(pd)*aly(pd)
        alphax(1)=dalx(0) *aly(0)
        alphax(2)=dalx(pd)*aly(0)
        alphax(3)=dalx(0) *aly(pd)
        alphax(4)=dalx(pd)*aly(pd)
        alphay(1)=alx(0) *daly(0)
        alphay(2)=alx(pd)*daly(0)
        alphay(3)=alx(0) *daly(pd)
        alphay(4)=alx(pd)*daly(pd)
c-----------------------------------------------------------------------
c       horizontal side centered bases if needed.  the ordering is
c       bottom to top, then left to right in pairs.
c-----------------------------------------------------------------------
        k=5
        DO i=1,pd-1
          alpha (k)= alx(i)* aly(0) 
          alphax(k)=dalx(i)* aly(0) 
          alphay(k)= alx(i)*daly(0) 
          k=k+1
          alpha (k)= alx(i)* aly(pd) 
          alphax(k)=dalx(i)* aly(pd) 
          alphay(k)= alx(i)*daly(pd) 
          k=k+1
        ENDDO
c-----------------------------------------------------------------------
c       vertical side centered bases if needed.  the ordering is
c       left to right, then bottom to top in pairs.
c-----------------------------------------------------------------------
        DO j=1,pd-1
          alpha (k)= alx(0)* aly(j) 
          alphax(k)=dalx(0)* aly(j) 
          alphay(k)= alx(0)*daly(j) 
          k=k+1
          alpha (k)= alx(pd)* aly(j) 
          alphax(k)=dalx(pd)* aly(j) 
          alphay(k)= alx(pd)*daly(j) 
          k=k+1
        ENDDO
c-----------------------------------------------------------------------
c       interior bases if needed.  the ordering is
c       left to right, bottom to top.
c-----------------------------------------------------------------------
        DO j=1,pd-1
          DO i=1,pd-1
            alpha (k)= alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_bases
c-----------------------------------------------------------------------
c     subprogram 2. lagr_quad_3D_alloc.
c     allocates space for lagr_quad_3D_type.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_alloc(laq,mx,my,nqty,nfour,poly_degree,
     $                              name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,poly_degree
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(lagr_quad_type), INTENT(OUT) :: laq

      INTEGER(i4) :: ix,iy,ib,ip,pd1
      REAL(r8), DIMENSION(0:poly_degree) :: x_node
c-----------------------------------------------------------------------
c     store grid, vector, and fourier series dimensions, and set the
c     number of side and interior basis functions.
c-----------------------------------------------------------------------
      laq%mx=mx
      laq%my=my
      laq%nqty=nqty
      laq%nfour=nfour
      laq%n_side=poly_degree-1
      laq%n_int=(poly_degree-1)**2
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx,0:my,nfour))
      IF (poly_degree>1) THEN
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx,0:my,nfour))
        ALLOCATE(laq%fsv(nqty,laq%n_side,0:mx,1:my,nfour))
        ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my,nfour))
      ENDIF
      ALLOCATE(laq%title(nqty))
      ALLOCATE(laq%f(nqty,nfour))
      ALLOCATE(laq%fx(nqty,nfour))
      ALLOCATE(laq%fy(nqty,nfour))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) laq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          laq%title=title
        ELSE
          laq%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     compute basis starting indices and centerings in logical space.
c     note that the centerings are for standard lagrange polynomial
c     bases.
c
c     grid-vertex data has grid indices running (0:mx,0:my), so
c       ix0=0 and iy0=0.
c     horizontal side data has grid indices running (1:mx,0:my), so
c       ix0=1 and iy0=0.
c     vertical side data has grid indices running (0:mx,1:my), so
c       ix0=0 and iy0=1.
c     interior data has grid indices running (1:mx,1:my), so
c       ix0=1 and iy0=1.
c
c     dx and dy give the relative positions within a cell for the
c     different bases.  note that these quantities may be used to access
c     all data with loops like:
c
c     DO ibase=1,SIZE(laq%ix0)
c       DO iy=laq%iy0(ibase),laq%my
c         DO ix=laq%ix0(ibase),laq%mx
c           dx_data=ix-laq%ix0(ibase)+laq%dx(ibase)
c           dy_data=iy-laq%iy0(ibase)+laq%dy(ibase)
c           CALL lagr_quad_eval(laq,dx_data,dy_data,1_i4)
c           xxx=laq%f(yyy,zzz)
c         ENDDO
c       ENDDO
c     ENDDO
c-----------------------------------------------------------------------
      pd1=poly_degree-1
      ALLOCATE(laq%ix0(poly_degree**2))
      ALLOCATE(laq%iy0(poly_degree**2))
      ALLOCATE(laq%dx(poly_degree**2))
      ALLOCATE(laq%dy(poly_degree**2))
      CALL poly_nodes(poly_degree,x_node)
      DO ib=1,poly_degree**2
        IF (ib==1) THEN
          laq%ix0(ib)=0
          laq%iy0(ib)=0
          laq%dx(ib)=0
          laq%dy(ib)=0
        ELSE IF (ib<=poly_degree) THEN
          laq%ix0(ib)=1
          laq%iy0(ib)=0
          laq%dx(ib)=x_node(ib-1_i4)
          laq%dy(ib)=0
        ELSE IF (ib<2*poly_degree) THEN
          laq%ix0(ib)=0
          laq%iy0(ib)=1
          laq%dx(ib)=0
          laq%dy(ib)=x_node(ib-poly_degree)
        ELSE
          laq%ix0(ib)=1
          laq%iy0(ib)=1
          ip=ib-2*poly_degree
          ix=MODULO(ip,pd1)+1
          iy=ip/pd1+1
          laq%dx(ib)=x_node(ix)
          laq%dy(ib)=x_node(iy)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_alloc
c-----------------------------------------------------------------------
c     subprogram 3. lagr_quad_3D_dealloc.
c     deallocates space for lagr_quad_3D_type.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_dealloc(laq)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(laq%fs)
      IF (ALLOCATED(laq%fsh)) DEALLOCATE(laq%fsh,laq%fsv,laq%fsi)
      DEALLOCATE(laq%title)
      DEALLOCATE(laq%f)
      DEALLOCATE(laq%fx)
      DEALLOCATE(laq%fy)
      DEALLOCATE(laq%ix0)
      DEALLOCATE(laq%iy0)
      DEALLOCATE(laq%dx)
      DEALLOCATE(laq%dy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_dealloc
c-----------------------------------------------------------------------
c     subprogram 4. lagr_quad_3D_eval.
c     evaluates complex lagr_quad quantities at a single point within a
c     grid block.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_eval(laq,x,y,dmode)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,pd,i,j,k,im
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
c     find the interval, and compute 1D basis coefficients.
c-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),laq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),laq%my-1),0_i4)
      pd=laq%n_side+1
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode)
      CALL lagr_1D(pd,y-iy,aly,daly,dmode)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  bilinear
c     derivatives simplify.
c-----------------------------------------------------------------------
      IF (pd==1) THEN
        laq%f=(laq%fs(:,ix  ,iy  ,:)*alx(0)
     $        +laq%fs(:,ix+1,iy  ,:)*alx(1))*aly(0)+
     $        (laq%fs(:,ix  ,iy+1,:)*alx(0)
     $        +laq%fs(:,ix+1,iy+1,:)*alx(1))*aly(1)
        IF (dmode==0) RETURN
        laq%fx=(laq%fs(:,ix+1,iy  ,:)
     $         -laq%fs(:,ix  ,iy  ,:))*aly(0)+
     $         (laq%fs(:,ix+1,iy+1,:)
     $         -laq%fs(:,ix  ,iy+1,:))*aly(1)
        laq%fy=(laq%fs(:,ix  ,iy+1,:)
     $         -laq%fs(:,ix  ,iy  ,:))*alx(0)+
     $         (laq%fs(:,ix+1,iy+1,:)
     $         -laq%fs(:,ix+1,iy  ,:))*alx(1)
c-----------------------------------------------------------------------
c     other polynomials:
c-----------------------------------------------------------------------
      ELSE
        laq%f=(laq%fs(:,ix  ,iy  ,:)*alx(0)
     $        +laq%fs(:,ix+1,iy  ,:)*alx(pd))*aly(0)+
     $        (laq%fs(:,ix  ,iy+1,:)*alx(0)
     $        +laq%fs(:,ix+1,iy+1,:)*alx(pd))*aly(pd)
        DO i=1,pd-1
          laq%f=laq%f+
     $          (laq%fsh(:,i,ix+1,iy  ,:)*aly(0)
     $          +laq%fsh(:,i,ix+1,iy+1,:)*aly(pd))*alx(i)+
     $          (laq%fsv(:,i,ix  ,iy+1,:)*alx(0)
     $          +laq%fsv(:,i,ix+1,iy+1,:)*alx(pd))*aly(i)
        ENDDO
        k=1
        DO j=1,pd-1
          DO i=1,pd-1
            laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        IF (dmode==0) RETURN
        laq%fx=(laq%fs(:,ix  ,iy  ,:)*dalx(0)
     $         +laq%fs(:,ix+1,iy  ,:)*dalx(pd))*aly(0)+
     $         (laq%fs(:,ix  ,iy+1,:)*dalx(0)
     $         +laq%fs(:,ix+1,iy+1,:)*dalx(pd))*aly(pd)
        laq%fy=(laq%fs(:,ix  ,iy  ,:)*alx(0)
     $         +laq%fs(:,ix+1,iy  ,:)*alx(pd))*daly(0)+
     $         (laq%fs(:,ix  ,iy+1,:)*alx(0)
     $         +laq%fs(:,ix+1,iy+1,:)*alx(pd))*daly(pd)
        DO i=1,pd-1
          laq%fx=laq%fx+
     $           (laq%fsh(:,i,ix+1,iy  ,:)* aly(0)
     $           +laq%fsh(:,i,ix+1,iy+1,:)* aly(pd))*dalx(i)+
     $           (laq%fsv(:,i,ix  ,iy+1,:)*dalx(0)
     $           +laq%fsv(:,i,ix+1,iy+1,:)*dalx(pd))*aly(i)
          laq%fy=laq%fy+
     $           (laq%fsh(:,i,ix+1,iy  ,:)*daly(0)
     $           +laq%fsh(:,i,ix+1,iy+1,:)*daly(pd))*alx(i)+
     $           (laq%fsv(:,i,ix  ,iy+1,:)* alx(0)
     $           +laq%fsv(:,i,ix+1,iy+1,:)* alx(pd))*daly(i)
        ENDDO
        k=1
        DO j=1,pd-1
          DO i=1,pd-1
            laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_eval
c-----------------------------------------------------------------------
c     subprogram 5. lagr_quad_3D_all_eval.
c     evaluates complex lagr_quad quantities in all elements in a grid
c     block for equal spacing.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_all_eval(laq,x,y,f,fx,fy,dmode)

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      REAL(r8), INTENT(IN) :: x,y
      COMPLEX(r8), INTENT(OUT), DIMENSION(:,:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: mx,my,mx1,my1
      INTEGER(i4) :: pd,i,j,k,im,ix,iy
      REAL(r8), DIMENSION((laq%n_side+2)**2) :: alpha,dalpdx,dalpdy
c-----------------------------------------------------------------------
c     compute index limits.
c-----------------------------------------------------------------------
      mx=laq%mx
      mx1=mx-1
      my=laq%my
      my1=my-1
      pd=laq%n_side+1
      CALL lagr_quad_bases(x,y,alpha,dalpdx,dalpdy,dmode)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  bilinear
c     derivatives simplify.  [do-loops are for optimization, not
c     aesthetics.]
c-----------------------------------------------------------------------
      SELECT CASE(pd)
      CASE(1)
        IF (dmode==0) THEN
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                f(:,ix,iy,im)=laq%fs(:,ix-1,iy-1,im)*alpha(1)+
     $                        laq%fs(:,ix  ,iy-1,im)*alpha(2)+
     $                        laq%fs(:,ix-1,iy  ,im)*alpha(3)+
     $                        laq%fs(:,ix  ,iy  ,im)*alpha(4) 
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                f(:,ix,iy,im)=laq%fs(:,ix-1,iy-1,im)*alpha(1)+
     $                        laq%fs(:,ix  ,iy-1,im)*alpha(2)+
     $                        laq%fs(:,ix-1,iy  ,im)*alpha(3)+
     $                        laq%fs(:,ix  ,iy  ,im)*alpha(4) 
                fx(:,ix,iy,im)=(laq%fs(:,ix  ,iy-1,im)
     $                         -laq%fs(:,ix-1,iy-1,im))*aly(0)+
     $                         (laq%fs(:,ix  ,iy  ,im)
     $                         -laq%fs(:,ix-1,iy  ,im))*aly(1)
                fy(:,ix,iy,im)=(laq%fs(:,ix-1,iy  ,im)
     $                         -laq%fs(:,ix-1,iy-1,im))*alx(0)+
     $                         (laq%fs(:,ix  ,iy  ,im)
     $                         -laq%fs(:,ix  ,iy-1,im))*alx(1)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     biquadratics--no sums.
c-----------------------------------------------------------------------
      CASE(2)
        IF (dmode==0) THEN
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                f(:,ix,iy,im)=
     $               laq%fs(:,ix-1,iy-1,im)*alpha(1)+
     $               laq%fs(:,ix  ,iy-1,im)*alpha(2)+
     $               laq%fs(:,ix-1,iy  ,im)*alpha(3)+
     $               laq%fs(:,ix  ,iy  ,im)*alpha(4)+
     $               laq%fsh(:,1,ix,iy-1,im)*alpha(5)+
     $               laq%fsh(:,1,ix,iy  ,im)*alpha(6)+
     $               laq%fsv(:,1,ix-1,iy,im)*alpha(7)+
     $               laq%fsv(:,1,ix  ,iy,im)*alpha(8)+
     $               laq%fsi(:,1,ix,iy,im)*alpha(9)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                f(:,ix,iy,im)=
     $               laq%fs(:,ix-1,iy-1,im)*alpha(1)+
     $               laq%fs(:,ix  ,iy-1,im)*alpha(2)+
     $               laq%fs(:,ix-1,iy  ,im)*alpha(3)+
     $               laq%fs(:,ix  ,iy  ,im)*alpha(4)+
     $               laq%fsh(:,1,ix,iy-1,im)*alpha(5)+
     $               laq%fsh(:,1,ix,iy  ,im)*alpha(6)+
     $               laq%fsv(:,1,ix-1,iy,im)*alpha(7)+
     $               laq%fsv(:,1,ix  ,iy,im)*alpha(8)+
     $               laq%fsi(:,1,ix,iy,im)*alpha(9)
                fx(:,ix,iy,im)=
     $                laq%fs(:,ix-1,iy-1,im)*dalpdx(1)+
     $                laq%fs(:,ix  ,iy-1,im)*dalpdx(2)+
     $                laq%fs(:,ix-1,iy  ,im)*dalpdx(3)+
     $                laq%fs(:,ix  ,iy  ,im)*dalpdx(4)+ 
     $                laq%fsh(:,1,ix,iy-1,im)*dalpdx(5)+
     $                laq%fsh(:,1,ix,iy  ,im)*dalpdx(6)+
     $                laq%fsv(:,1,ix-1,iy,im)*dalpdx(7)+
     $                laq%fsv(:,1,ix  ,iy,im)*dalpdx(8)+
     $                laq%fsi(:,1,ix,iy,im)*dalpdx(9)
                fy(:,ix,iy,im)=
     $                laq%fs(:,ix-1,iy-1,im)*dalpdy(1)+
     $                laq%fs(:,ix  ,iy-1,im)*dalpdy(2)+
     $                laq%fs(:,ix-1,iy  ,im)*dalpdy(3)+
     $                laq%fs(:,ix  ,iy  ,im)*dalpdy(4)+ 
     $                laq%fsh(:,1,ix,iy-1,im)*dalpdy(5)+
     $                laq%fsh(:,1,ix,iy  ,im)*dalpdy(6)+
     $                laq%fsv(:,1,ix-1,iy,im)*dalpdy(7)+
     $                laq%fsv(:,1,ix  ,iy,im)*dalpdy(8)+
     $                laq%fsi(:,1,ix,iy,im)*dalpdy(9)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     other polynomials:
c-----------------------------------------------------------------------
      CASE DEFAULT
        j=2*pd+2
        k=4*pd
        IF (dmode==0) THEN
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                DO i=1,laq%nqty
                  f(i,ix,iy,im)=
     $                 laq%fs(i,ix-1,iy-1,im)*alpha(1)+
     $                 laq%fs(i,ix  ,iy-1,im)*alpha(2)+
     $                 laq%fs(i,ix-1,iy  ,im)*alpha(3)+
     $                 laq%fs(i,ix  ,iy  ,im)*alpha(4)+
     $             SUM(laq%fsh(i,:,ix,iy-1,im)*alpha(5:j-1:2)+
     $                 laq%fsh(i,:,ix,iy  ,im)*alpha(6:j  :2)+
     $                 laq%fsv(i,:,ix-1,iy,im)*alpha(j+1:k-1:2)+
     $                 laq%fsv(i,:,ix  ,iy,im)*alpha(j+2:k:2))+
     $             SUM(laq%fsi(i,:,ix,iy,im)*alpha(k+1:))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO im=1,laq%nfour
            DO iy=1,my
              DO ix=1,mx
                DO i=1,laq%nqty
                  f(i,ix,iy,im)=
     $                 laq%fs(i,ix-1,iy-1,im)*alpha(1)+
     $                 laq%fs(i,ix  ,iy-1,im)*alpha(2)+
     $                 laq%fs(i,ix-1,iy  ,im)*alpha(3)+
     $                 laq%fs(i,ix  ,iy  ,im)*alpha(4)+
     $             SUM(laq%fsh(i,:,ix,iy-1,im)*alpha(5:j-1:2)+
     $                 laq%fsh(i,:,ix,iy  ,im)*alpha(6:j  :2)+
     $                 laq%fsv(i,:,ix-1,iy,im)*alpha(j+1:k-1:2)+
     $                 laq%fsv(i,:,ix  ,iy,im)*alpha(j+2:k:2))+
     $             SUM(laq%fsi(i,:,ix,iy,im)*alpha(k+1:))
                  fx(i,ix,iy,im)=
     $                  laq%fs(i,ix-1,iy-1,im)*dalpdx(1)+
     $                  laq%fs(i,ix  ,iy-1,im)*dalpdx(2)+
     $                  laq%fs(i,ix-1,iy  ,im)*dalpdx(3)+
     $                  laq%fs(i,ix  ,iy  ,im)*dalpdx(4)+ 
     $              SUM(laq%fsh(i,:,ix,iy-1,im)*dalpdx(5:j-1:2)+
     $                  laq%fsh(i,:,ix,iy  ,im)*dalpdx(6:j  :2)+
     $                  laq%fsv(i,:,ix-1,iy,im)*dalpdx(j+1:k-1:2)+
     $                  laq%fsv(i,:,ix  ,iy,im)*dalpdx(j+2:k:2))+
     $              SUM(laq%fsi(i,:,ix,iy,im)*dalpdx(k+1:))
                  fy(i,ix,iy,im)=
     $                  laq%fs(i,ix-1,iy-1,im)*dalpdy(1)+
     $                  laq%fs(i,ix  ,iy-1,im)*dalpdy(2)+
     $                  laq%fs(i,ix-1,iy  ,im)*dalpdy(3)+
     $                  laq%fs(i,ix  ,iy  ,im)*dalpdy(4)+ 
     $              SUM(laq%fsh(i,:,ix,iy-1,im)*dalpdy(5:j-1:2)+
     $                  laq%fsh(i,:,ix,iy  ,im)*dalpdy(6:j  :2)+
     $                  laq%fsv(i,:,ix-1,iy,im)*dalpdy(j+1:k-1:2)+
     $                  laq%fsv(i,:,ix  ,iy,im)*dalpdy(j+2:k:2))+
     $              SUM(laq%fsi(i,:,ix,iy,im)*dalpdy(k+1:))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_all_eval
c-----------------------------------------------------------------------
c     subprogram 6. lagr_quad_3D_assign_rsc.
c     assign a real scalar value to a complex lagrange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_rsc(laq,rscalar)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: rscalar

      laq%fs=rscalar
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=rscalar
        laq%fsv=rscalar
        laq%fsi=rscalar
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 7. lagr_quad_3D_assign_csc.
c     assign a complex scalar value to a complex lagrange quad
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_csc(laq,cscalar)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), INTENT(IN) :: cscalar

      laq%fs=cscalar
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=cscalar
        laq%fsv=cscalar
        laq%fsi=cscalar
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 8. lagr_quad_3D_assign_laq.
c     set one complex lagrange quad structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_laq(laq1,laq2)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq1
      TYPE(lagr_quad_type), INTENT(IN) :: laq2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      laq1%fs=laq2%fs
      IF (ALLOCATED(laq1%fsh)) THEN
        laq1%fsh=laq2%fsh
        laq1%fsv=laq2%fsv
        laq1%fsi=laq2%fsi
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_assign_laq
c-----------------------------------------------------------------------
c     subprogram 9. lagr_quad_3D_assign_int.
c     assign a integer value to a complex lagrange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_assign_int(laq,int)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      INTEGER(i4), INTENT(IN) :: int

      laq%fs=int
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=int
        laq%fsv=int
        laq%fsi=int
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_assign_int
c-----------------------------------------------------------------------
c     subprogram 10. lagr_quad_3D_basis_assign_arr
c     assign data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_basis_assign_arr(laq,data,ibasis)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs=data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,:,:,:)=data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,:,:,:)=data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,:,:,:)=data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_basis_assign_arr
c-----------------------------------------------------------------------
c     subprogram 11. lagr_quad_3D_basis_add_arr
c     add data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_basis_add_arr(laq,data,ibasis)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs=laq%fs+data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,:,:,:)=laq%fsh(:,ibasis-1,:,:,:)+data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,:,:,:)=
     $    laq%fsv(:,ibasis-poly_degree,:,:,:)+data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,:,:,:)=
     $     laq%fsi(:,ibasis-2*poly_degree+1,:,:,:)+data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_basis_add_arr
c-----------------------------------------------------------------------
c     subprogram 12. lagr_quad_3D_basis_assign_loc
c     assign data into coefficient arrays for one basis function.
c
c     this is a local version of lagr_quad_basis_assign_arr, where the
c     the data is located at a given poloidal and fourier indices
c     triplet, only.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_basis_assign_loc(laq,data,ibasis,ix,iy,im)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs(:,ix,iy,im)=data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,ix,iy,im)=data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,ix,iy,im)=data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy,im)=data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_basis_assign_loc
c-----------------------------------------------------------------------
c     subprogram 13. lagr_quad_3D_basis_add_loc
c     add data into coefficient arrays for one basis function.
c
c     this is a local version of lagr_quad_basis_add_arr, where the
c     the data is located at a given poloidal and fourier indices
c     triplet, only.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_3D_basis_add_loc(laq,data,ibasis,ix,iy,im)

      TYPE(lagr_quad_type), INTENT(INOUT) :: laq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs(:,ix,iy,im)=laq%fs(:,ix,iy,im)+data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,ix,iy,im)=laq%fsh(:,ibasis-1,ix,iy,im)+data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,ix,iy,im)=
     $    laq%fsv(:,ibasis-poly_degree,ix,iy,im)+data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy,im)=
     $    laq%fsi(:,ibasis-2*poly_degree+1,ix,iy,im)+data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_3D_basis_add_loc
c-----------------------------------------------------------------------
c     subprogram 14. lagr_quad_2D_alloc.
c     allocates space for lagr_quad_2D_type.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_alloc(laq,mx,my,nqty,poly_degree,
     $                              name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq

      INTEGER(i4) :: ix,iy,ib,ip,pd1
      REAL(r8), DIMENSION(0:poly_degree) :: x_node
c-----------------------------------------------------------------------
c     store grid and vector dimensions, and set the
c     number of side and interior basis functions.
c-----------------------------------------------------------------------
      laq%mx=mx
      laq%my=my
      laq%nqty=nqty
      laq%n_side=poly_degree-1
      laq%n_int=(poly_degree-1)**2
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(laq%fs(nqty,0:mx,0:my))
      IF (poly_degree>1) THEN
        ALLOCATE(laq%fsh(nqty,laq%n_side,1:mx,0:my))
        ALLOCATE(laq%fsv(nqty,laq%n_side,0:mx,1:my))
        ALLOCATE(laq%fsi(nqty,laq%n_int,1:mx,1:my))
      ENDIF
      ALLOCATE(laq%title(nqty))
      ALLOCATE(laq%f(nqty))
      ALLOCATE(laq%fx(nqty))
      ALLOCATE(laq%fy(nqty))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) laq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          laq%title=title
        ELSE
          laq%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     compute basis starting indices and centerings in logical space.
c     note that the centerings are for standard lagrange polynomial
c     bases.
c
c     grid-vertex data has grid indices running (0:mx,0:my), so
c       ix0=0 and iy0=0.
c     horizontal side data has grid indices running (1:mx,0:my), so
c       ix0=1 and iy0=0.
c     vertical side data has grid indices running (0:mx,1:my), so
c       ix0=0 and iy0=1.
c     interior data has grid indices running (1:mx,1:my), so
c       ix0=1 and iy0=1.
c
c     dx and dy give the relative positions within a cell for the
c     different bases.  note that these quantities may be used to access
c     all data with loops like:
c
c     DO ibase=1,SIZE(laq%ix0)
c       DO iy=laq%iy0(ibase),laq%my
c         DO ix=laq%ix0(ibase),laq%mx
c           dx_data=ix-laq%ix0(ibase)+laq%dx(ibase)
c           dy_data=iy-laq%iy0(ibase)+laq%dy(ibase)
c           CALL lagr_quad_eval(laq,dx_data,dy_data,1_i4)
c           xxx=laq%f(yyy)
c         ENDDO
c       ENDDO
c     ENDDO
c-----------------------------------------------------------------------
      pd1=poly_degree-1
      ALLOCATE(laq%ix0(poly_degree**2))
      ALLOCATE(laq%iy0(poly_degree**2))
      ALLOCATE(laq%dx(poly_degree**2))
      ALLOCATE(laq%dy(poly_degree**2))
      CALL poly_nodes(poly_degree,x_node)
      DO ib=1,poly_degree**2
        IF (ib==1) THEN
          laq%ix0(ib)=0
          laq%iy0(ib)=0
          laq%dx(ib)=0
          laq%dy(ib)=0
        ELSE IF (ib<=poly_degree) THEN
          laq%ix0(ib)=1
          laq%iy0(ib)=0
          laq%dx(ib)=x_node(ib-1_i4)
          laq%dy(ib)=0
        ELSE IF (ib<2*poly_degree) THEN
          laq%ix0(ib)=0
          laq%iy0(ib)=1
          laq%dx(ib)=0
          laq%dy(ib)=x_node(ib-poly_degree)
        ELSE
          laq%ix0(ib)=1
          laq%iy0(ib)=1
          ip=ib-2*poly_degree
          ix=MODULO(ip,pd1)+1
          iy=ip/pd1+1
          laq%dx(ib)=x_node(ix)
          laq%dy(ib)=x_node(iy)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_alloc
c-----------------------------------------------------------------------
c     subprogram 15. lagr_quad_2D_dealloc.
c     deallocates space for lagr_quad_2D_type.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_dealloc(laq)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(laq%fs)
      IF (ALLOCATED(laq%fsh)) DEALLOCATE(laq%fsh,laq%fsv,laq%fsi)
      DEALLOCATE(laq%title)
      DEALLOCATE(laq%f)
      DEALLOCATE(laq%fx)
      DEALLOCATE(laq%fy)
      DEALLOCATE(laq%ix0)
      DEALLOCATE(laq%iy0)
      DEALLOCATE(laq%dx)
      DEALLOCATE(laq%dy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_dealloc
c-----------------------------------------------------------------------
c     subprogram 16. lagr_quad_2D_eval.
c     evaluates real lagr_quad quantities at a single point within a
c     grid block.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_eval(laq,x,y,dmode)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,pd,i,j,k,im
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
c     find the interval, and compute 1D basis coefficients.
c-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),laq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),laq%my-1),0_i4)
      pd=laq%n_side+1
      CALL lagr_1D(pd,x-ix,alx,dalx,dmode)
      CALL lagr_1D(pd,y-iy,aly,daly,dmode)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  bilinear
c     derivatives simplify.
c-----------------------------------------------------------------------
      IF (pd==1) THEN
        laq%f=(laq%fs(:,ix  ,iy  )*alx(0)
     $        +laq%fs(:,ix+1,iy  )*alx(1))*aly(0)+
     $        (laq%fs(:,ix  ,iy+1)*alx(0)
     $        +laq%fs(:,ix+1,iy+1)*alx(1))*aly(1)
        IF (dmode==0) RETURN
        laq%fx=(laq%fs(:,ix+1,iy  )
     $         -laq%fs(:,ix  ,iy  ))*aly(0)+
     $         (laq%fs(:,ix+1,iy+1)
     $         -laq%fs(:,ix  ,iy+1))*aly(1)
        laq%fy=(laq%fs(:,ix  ,iy+1)
     $         -laq%fs(:,ix  ,iy  ))*alx(0)+
     $         (laq%fs(:,ix+1,iy+1)
     $         -laq%fs(:,ix+1,iy  ))*alx(1)
c-----------------------------------------------------------------------
c     other polynomials:
c-----------------------------------------------------------------------
      ELSE
        laq%f=(laq%fs(:,ix  ,iy  )*alx(0)
     $        +laq%fs(:,ix+1,iy  )*alx(pd))*aly(0)+
     $        (laq%fs(:,ix  ,iy+1)*alx(0)
     $        +laq%fs(:,ix+1,iy+1)*alx(pd))*aly(pd)
        DO i=1,pd-1
          laq%f=laq%f+
     $          (laq%fsh(:,i,ix+1,iy  )*aly(0)
     $          +laq%fsh(:,i,ix+1,iy+1)*aly(pd))*alx(i)+
     $          (laq%fsv(:,i,ix  ,iy+1)*alx(0)
     $          +laq%fsv(:,i,ix+1,iy+1)*alx(pd))*aly(i)
        ENDDO
        k=1
        DO j=1,pd-1
          DO i=1,pd-1
            laq%f=laq%f+laq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        IF (dmode==0) RETURN
        laq%fx=(laq%fs(:,ix  ,iy  )*dalx(0)
     $         +laq%fs(:,ix+1,iy  )*dalx(pd))*aly(0)+
     $         (laq%fs(:,ix  ,iy+1)*dalx(0)
     $         +laq%fs(:,ix+1,iy+1)*dalx(pd))*aly(pd)
        laq%fy=(laq%fs(:,ix  ,iy  )*alx(0)
     $         +laq%fs(:,ix+1,iy  )*alx(pd))*daly(0)+
     $         (laq%fs(:,ix  ,iy+1)*alx(0)
     $         +laq%fs(:,ix+1,iy+1)*alx(pd))*daly(pd)
        DO i=1,pd-1
          laq%fx=laq%fx+
     $           (laq%fsh(:,i,ix+1,iy  )* aly(0)
     $           +laq%fsh(:,i,ix+1,iy+1)* aly(pd))*dalx(i)+
     $           (laq%fsv(:,i,ix  ,iy+1)*dalx(0)
     $           +laq%fsv(:,i,ix+1,iy+1)*dalx(pd))*aly(i)
          laq%fy=laq%fy+
     $           (laq%fsh(:,i,ix+1,iy  )*daly(0)
     $           +laq%fsh(:,i,ix+1,iy+1)*daly(pd))*alx(i)+
     $           (laq%fsv(:,i,ix  ,iy+1)* alx(0)
     $           +laq%fsv(:,i,ix+1,iy+1)* alx(pd))*daly(i)
        ENDDO
        k=1
        DO j=1,pd-1
          DO i=1,pd-1
            laq%fx=laq%fx+laq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            laq%fy=laq%fy+laq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_eval
c-----------------------------------------------------------------------
c     subprogram 17. lagr_quad_2D_all_eval.
c     evaluates real lagr_quad quantities in all elements in a grid
c     block for equal spacing.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_all_eval(laq,x,y,f,fx,fy,dmode)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: mx,my,mx1,my1
      INTEGER(i4) :: pd,i,j,k,ix,iy
      REAL(r8), DIMENSION((laq%n_side+2)**2) :: alpha,dalpdx,dalpdy
c-----------------------------------------------------------------------
c     compute index limits.
c-----------------------------------------------------------------------
      mx=laq%mx
      mx1=mx-1
      my=laq%my
      my1=my-1
      pd=laq%n_side+1
      CALL lagr_quad_bases(x,y,alpha,dalpdx,dalpdy,dmode)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  bilinear
c     derivatives simplify.  [do-loops are for optimization, not
c     aesthetics.]
c-----------------------------------------------------------------------
      SELECT CASE(pd)
      CASE(1)
        IF (dmode==0) THEN
          DO iy=1,my
            DO ix=1,mx
              f(:,ix,iy)=laq%fs(:,ix-1,iy-1)*alpha(1)+
     $                   laq%fs(:,ix  ,iy-1)*alpha(2)+
     $                   laq%fs(:,ix-1,iy  )*alpha(3)+
     $                   laq%fs(:,ix  ,iy  )*alpha(4) 
            ENDDO
          ENDDO
        ELSE
          DO iy=1,my
            DO ix=1,mx
              f(:,ix,iy)=laq%fs(:,ix-1,iy-1)*alpha(1)+
     $                   laq%fs(:,ix  ,iy-1)*alpha(2)+
     $                   laq%fs(:,ix-1,iy  )*alpha(3)+
     $                   laq%fs(:,ix  ,iy  )*alpha(4) 
              fx(:,ix,iy)=(laq%fs(:,ix  ,iy-1)
     $                    -laq%fs(:,ix-1,iy-1))*aly(0)+
     $                    (laq%fs(:,ix  ,iy  )
     $                    -laq%fs(:,ix-1,iy  ))*aly(1)
              fy(:,ix,iy)=(laq%fs(:,ix-1,iy  )
     $                    -laq%fs(:,ix-1,iy-1))*alx(0)+
     $                    (laq%fs(:,ix  ,iy  )
     $                    -laq%fs(:,ix  ,iy-1))*alx(1)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     biquadratics--no sums.
c-----------------------------------------------------------------------
      CASE(2)
        IF (dmode==0) THEN
          DO iy=1,my
            DO ix=1,mx
              f(:,ix,iy)=
     $             laq%fs(:,ix-1,iy-1)*alpha(1)+
     $             laq%fs(:,ix  ,iy-1)*alpha(2)+
     $             laq%fs(:,ix-1,iy  )*alpha(3)+
     $             laq%fs(:,ix  ,iy  )*alpha(4)+
     $             laq%fsh(:,1,ix,iy-1)*alpha(5)+
     $             laq%fsh(:,1,ix,iy  )*alpha(6)+
     $             laq%fsv(:,1,ix-1,iy)*alpha(7)+
     $             laq%fsv(:,1,ix  ,iy)*alpha(8)+
     $             laq%fsi(:,1,ix,iy)*alpha(9)
            ENDDO
          ENDDO
        ELSE
          DO iy=1,my
            DO ix=1,mx
              f(:,ix,iy)=
     $             laq%fs(:,ix-1,iy-1)*alpha(1)+
     $             laq%fs(:,ix  ,iy-1)*alpha(2)+
     $             laq%fs(:,ix-1,iy  )*alpha(3)+
     $             laq%fs(:,ix  ,iy  )*alpha(4)+
     $             laq%fsh(:,1,ix,iy-1)*alpha(5)+
     $             laq%fsh(:,1,ix,iy  )*alpha(6)+
     $             laq%fsv(:,1,ix-1,iy)*alpha(7)+
     $             laq%fsv(:,1,ix  ,iy)*alpha(8)+
     $             laq%fsi(:,1,ix,iy)*alpha(9)
              fx(:,ix,iy)=
     $              laq%fs(:,ix-1,iy-1)*dalpdx(1)+
     $              laq%fs(:,ix  ,iy-1)*dalpdx(2)+
     $              laq%fs(:,ix-1,iy  )*dalpdx(3)+
     $              laq%fs(:,ix  ,iy  )*dalpdx(4)+ 
     $              laq%fsh(:,1,ix,iy-1)*dalpdx(5)+
     $              laq%fsh(:,1,ix,iy  )*dalpdx(6)+
     $              laq%fsv(:,1,ix-1,iy)*dalpdx(7)+
     $              laq%fsv(:,1,ix  ,iy)*dalpdx(8)+
     $              laq%fsi(:,1,ix,iy)*dalpdx(9)
              fy(:,ix,iy)=
     $              laq%fs(:,ix-1,iy-1)*dalpdy(1)+
     $              laq%fs(:,ix  ,iy-1)*dalpdy(2)+
     $              laq%fs(:,ix-1,iy  )*dalpdy(3)+
     $              laq%fs(:,ix  ,iy  )*dalpdy(4)+ 
     $              laq%fsh(:,1,ix,iy-1)*dalpdy(5)+
     $              laq%fsh(:,1,ix,iy  )*dalpdy(6)+
     $              laq%fsv(:,1,ix-1,iy)*dalpdy(7)+
     $              laq%fsv(:,1,ix  ,iy)*dalpdy(8)+
     $              laq%fsi(:,1,ix,iy)*dalpdy(9)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     other polynomials:
c-----------------------------------------------------------------------
      CASE DEFAULT
        j=2*pd+2
        k=4*pd
        IF (dmode==0) THEN
          DO iy=1,my
            DO ix=1,mx
              DO i=1,laq%nqty
                f(i,ix,iy)=
     $               laq%fs(i,ix-1,iy-1)*alpha(1)+
     $               laq%fs(i,ix  ,iy-1)*alpha(2)+
     $               laq%fs(i,ix-1,iy  )*alpha(3)+
     $               laq%fs(i,ix  ,iy  )*alpha(4)+
     $           SUM(laq%fsh(i,:,ix,iy-1)*alpha(5:j-1:2)+
     $               laq%fsh(i,:,ix,iy  )*alpha(6:j  :2)+
     $               laq%fsv(i,:,ix-1,iy)*alpha(j+1:k-1:2)+
     $               laq%fsv(i,:,ix  ,iy)*alpha(j+2:k:2))+
     $           SUM(laq%fsi(i,:,ix,iy)*alpha(k+1:))
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO iy=1,my
            DO ix=1,mx
              DO i=1,laq%nqty
                f(i,ix,iy)=
     $               laq%fs(i,ix-1,iy-1)*alpha(1)+
     $               laq%fs(i,ix  ,iy-1)*alpha(2)+
     $               laq%fs(i,ix-1,iy  )*alpha(3)+
     $               laq%fs(i,ix  ,iy  )*alpha(4)+
     $           SUM(laq%fsh(i,:,ix,iy-1)*alpha(5:j-1:2)+
     $               laq%fsh(i,:,ix,iy  )*alpha(6:j  :2)+
     $               laq%fsv(i,:,ix-1,iy)*alpha(j+1:k-1:2)+
     $               laq%fsv(i,:,ix  ,iy)*alpha(j+2:k:2))+
     $           SUM(laq%fsi(i,:,ix,iy)*alpha(k+1:))
                fx(i,ix,iy)=
     $                laq%fs(i,ix-1,iy-1)*dalpdx(1)+
     $                laq%fs(i,ix  ,iy-1)*dalpdx(2)+
     $                laq%fs(i,ix-1,iy  )*dalpdx(3)+
     $                laq%fs(i,ix  ,iy  )*dalpdx(4)+ 
     $            SUM(laq%fsh(i,:,ix,iy-1)*dalpdx(5:j-1:2)+
     $                laq%fsh(i,:,ix,iy  )*dalpdx(6:j  :2)+
     $                laq%fsv(i,:,ix-1,iy)*dalpdx(j+1:k-1:2)+
     $                laq%fsv(i,:,ix  ,iy)*dalpdx(j+2:k:2))+
     $            SUM(laq%fsi(i,:,ix,iy)*dalpdx(k+1:))
                fy(i,ix,iy)=
     $                laq%fs(i,ix-1,iy-1)*dalpdy(1)+
     $                laq%fs(i,ix  ,iy-1)*dalpdy(2)+
     $                laq%fs(i,ix-1,iy  )*dalpdy(3)+
     $                laq%fs(i,ix  ,iy  )*dalpdy(4)+ 
     $            SUM(laq%fsh(i,:,ix,iy-1)*dalpdy(5:j-1:2)+
     $                laq%fsh(i,:,ix,iy  )*dalpdy(6:j  :2)+
     $                laq%fsv(i,:,ix-1,iy)*dalpdy(j+1:k-1:2)+
     $                laq%fsv(i,:,ix  ,iy)*dalpdy(j+2:k:2))+
     $            SUM(laq%fsi(i,:,ix,iy)*dalpdy(k+1:))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_all_eval
c-----------------------------------------------------------------------
c     subprogram 18. lagr_quad_2D_assign_rsc.
c     assign a real scalar value to a real lagrange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_rsc(laq,rscalar)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), INTENT(IN) :: rscalar

      laq%fs=rscalar
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=rscalar
        laq%fsv=rscalar
        laq%fsi=rscalar
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 19. lagr_quad_2D_assign_csc.
c     assign a real scalar value to a real lagrange quad
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_csc(laq,cscalar)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      COMPLEX(r8), INTENT(IN) :: cscalar

      laq%fs=cscalar
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=cscalar
        laq%fsv=cscalar
        laq%fsi=cscalar
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 20. lagr_quad_2D_assign_laq.
c     set one real lagrange quad structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_laq(laq1,laq2)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq1
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      laq1%fs=laq2%fs
      IF (ALLOCATED(laq1%fsh)) THEN
        laq1%fsh=laq2%fsh
        laq1%fsv=laq2%fsv
        laq1%fsi=laq2%fsi
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_assign_laq
c-----------------------------------------------------------------------
c     subprogram 21. lagr_quad_2D_assign_int.
c     assign a integer value to a real lagrange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_assign_int(laq,int)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      INTEGER(i4), INTENT(IN) :: int

      laq%fs=int
      IF (ALLOCATED(laq%fsh)) THEN
        laq%fsh=int
        laq%fsv=int
        laq%fsi=int
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_assign_int
c-----------------------------------------------------------------------
c     subprogram 22. lagr_quad_2D_basis_assign_arr
c     assign data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_basis_assign_arr(laq,data,ibasis)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs=data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,:,:)=data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,:,:)=data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,:,:)=data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_basis_assign_arr
c-----------------------------------------------------------------------
c     subprogram 23. lagr_quad_2D_basis_add_arr
c     add data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_basis_add_arr(laq,data,ibasis)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs=laq%fs+data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,:,:)=laq%fsh(:,ibasis-1,:,:)+data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,:,:)=
     $    laq%fsv(:,ibasis-poly_degree,:,:)+data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,:,:)=
     $     laq%fsi(:,ibasis-2*poly_degree+1,:,:)+data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_basis_add_arr
c-----------------------------------------------------------------------
c     subprogram 24. lagr_quad_2D_basis_assign_loc
c     assign data into coefficient arrays for one basis function.
c
c     this is a local version of lagr_quad_basis_assign_arr, where the
c     the data is located at given poloidal indices, only.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_basis_assign_loc(laq,data,ibasis,ix,iy)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs(:,ix,iy)=data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,ix,iy)=data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,ix,iy)=data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy)=data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_basis_assign_loc
c-----------------------------------------------------------------------
c     subprogram 25. lagr_quad_2D_basis_add_loc
c     add data into coefficient arrays for one basis function.
c
c     this is a local version of lagr_quad_basis_add_arr, where the
c     the data is located at given poloidal indices only.
c-----------------------------------------------------------------------
      SUBROUTINE lagr_quad_2D_basis_add_loc(laq,data,ibasis,ix,iy)

      TYPE(lagr_quad_2D_type), INTENT(INOUT) :: laq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy

      INTEGER(i4) :: poly_degree
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c     grid vertex-centered data is first.
c-----------------------------------------------------------------------
      poly_degree=laq%n_side+1
      IF (ibasis==1) THEN
        laq%fs(:,ix,iy)=laq%fs(:,ix,iy)+data
        RETURN
c-----------------------------------------------------------------------
c     horizontal sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<=poly_degree) THEN
        laq%fsh(:,ibasis-1,ix,iy)=laq%fsh(:,ibasis-1,ix,iy)+data
        RETURN
c-----------------------------------------------------------------------
c     vertical sides.
c-----------------------------------------------------------------------
      ELSE IF (ibasis<2*poly_degree) THEN
        laq%fsv(:,ibasis-poly_degree,ix,iy)=
     $    laq%fsv(:,ibasis-poly_degree,ix,iy)+data
        RETURN
c-----------------------------------------------------------------------
c     interior bases.
c-----------------------------------------------------------------------
      ELSE
        laq%fsi(:,ibasis-2*poly_degree+1,ix,iy)=
     $    laq%fsi(:,ibasis-2*poly_degree+1,ix,iy)+data
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE lagr_quad_2D_basis_add_loc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE lagr_quad_mod
