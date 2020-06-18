c-----------------------------------------------------------------------
c     file modal_disc_quad.f
c     routines for evaluating discontinuous modal bases in blocks of
c     structured quadrilateral elements.
c
c     unlike lagr_disc_quad, these modal bases have their own type
c     definition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0.  modal_type_mod.
c      0.1 modal_disc_mod.
c      1.  modal_disc_bases.
c      2.  modal_disc_3D_alloc.
c      3.  modal_disc_3D_dealloc.
c      4.  modal_disc_3D_eval.
c      5.  modal_disc_3D_all_eval.
c      6.  modal_disc_3D_assign_rsc.
c      7.  modal_disc_3D_assign_csc.
c      8.  modal_disc_3D_assign_modq.
c      9.  modal_disc_3D_assign_int.
c      10. modal_disc_3D_basis_assign_arr
c      11. modal_disc_3D_basis_add_arr
c      12. modal_disc_3D_basis_assign_loc
c      13. modal_disc_3D_basis_add_loc
c      14. modal_disc_2D_alloc.
c      15. modal_disc_2D_dealloc.
c      16. modal_disc_2D_eval.
c      17. modal_disc_2D_all_eval.
c      18. modal_disc_2D_assign_rsc.
c      19. modal_disc_2D_assign_csc.
c      20. modal_disc_2D_assign_modq.
c      21. modal_disc_2D_assign_int.
c      22. modal_disc_2D_basis_assign_arr
c      23. modal_disc_2D_basis_add_arr
c      24. modal_disc_2D_basis_assign_loc
c      25. modal_disc_2D_basis_add_loc
c-----------------------------------------------------------------------
c     subprogram 0. modal_type_mod definition.
c     contains the F90 user-type definition for modal bases in
c     quadrilateral elements.  the bases only have arrays for element-
c     interior coefficients, and the starting indices are always 1.
c
c     a distinction from standard bases is that the polynomials here
c     can be incomplete representations.  the expansion can be described
c     as the outer product of Legendre polynomials of degrees pdmin to
c     pdmax in the logical coordinate x and the complete basis from
c     0 to pd in y, 'unioned' with the same with x and y swapped.
c     it is assumed that pdmax<=pd.
c
c     for example, if pdmin=pdmax=3 and pd=4, bases include the outer
c     product of the cubic Legendre polynomial in x with the full
c     quartic expansion in y, plus the cubic Legendre polynomial in y
c     with the full quartic expansion in x, without duplication.
c-----------------------------------------------------------------------
      MODULE modal_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: modal_quad_type
        INTEGER(i4) :: mx,my,nqty,nfour,n_int,pd,pdmin,pdmax
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: fsi
        COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        CHARACTER(6) :: name
      END TYPE modal_quad_type

      TYPE :: modal_quad_2D_type
        INTEGER(i4) :: mx,my,nqty,n_int,pd,pdmin,pdmax
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ix0,iy0
        REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: fsi
        REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy
        CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
        CHARACTER(6) :: name
      END TYPE modal_quad_2D_type

      END MODULE modal_type_mod
c-----------------------------------------------------------------------
c     subprogram 0.1 modal_disc_mod definition.
c     contains the subprograms and interfaces.
c-----------------------------------------------------------------------
      MODULE modal_disc_mod
      USE modal_type_mod
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: npoly_max=20
      REAL(r8), DIMENSION(0:npoly_max), PRIVATE :: alx,aly,dalx,daly
c-----------------------------------------------------------------------
c     subprogram name interfaces
c-----------------------------------------------------------------------
      INTERFACE modal_disc_alloc
        MODULE PROCEDURE modal_disc_2D_alloc,modal_disc_3D_alloc
      END INTERFACE

      INTERFACE modal_disc_dealloc
        MODULE PROCEDURE modal_disc_2D_dealloc,modal_disc_3D_dealloc
      END INTERFACE

      INTERFACE modal_disc_all_eval
        MODULE PROCEDURE modal_disc_2D_all_eval,modal_disc_3D_all_eval
      END INTERFACE

      INTERFACE modal_disc_eval
        MODULE PROCEDURE modal_disc_2D_eval,modal_disc_3D_eval
      END INTERFACE

      INTERFACE modal_disc_basis_assign_arr
        MODULE PROCEDURE
     $    modal_disc_2D_basis_assign_arr,modal_disc_3D_basis_assign_arr
      END INTERFACE

      INTERFACE modal_disc_basis_add_arr
        MODULE PROCEDURE
     $    modal_disc_2D_basis_add_arr,modal_disc_3D_basis_add_arr
      END INTERFACE

      INTERFACE modal_disc_basis_assign_loc
        MODULE PROCEDURE
     $    modal_disc_2D_basis_assign_loc,modal_disc_3D_basis_assign_loc
      END INTERFACE

      INTERFACE modal_disc_basis_add_loc
        MODULE PROCEDURE
     $    modal_disc_2D_basis_add_loc,modal_disc_3D_basis_add_loc
      END INTERFACE

      INTERFACE modal_disc_assignment
        MODULE PROCEDURE
     $    modal_disc_3D_assign_csc,modal_disc_3D_assign_rsc,
     $    modal_disc_3D_assign_modq,modal_disc_3D_assign_int,
     $    modal_disc_2D_assign_csc,modal_disc_2D_assign_rsc,
     $    modal_disc_2D_assign_modq,modal_disc_2D_assign_int
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. modal_disc_bases.
c     computes modal basis functions and their derivatives at a given
c     logical position (x,y) within an element, based on nimrod's
c     0<=x,y<=1.  note that the minimum polynomial degree (pdmin)
c     and the maximum degree in each coordinate (pd) must be specified.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_bases(x,y,alpha,alphax,alphay,dmode,
     $                            pd,pdmin,pdmax)

      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), DIMENSION(:), INTENT(INOUT) :: alpha,alphax,alphay
      INTEGER(i4), INTENT(IN) :: dmode,pd,pdmin,pdmax

      INTEGER(i4) :: i,j,k
      REAL(r8) :: xstd,ystd
c-----------------------------------------------------------------------
c     first determine the logical coordinates in the standard,
c     -1<=xstd,ystd,<=+1 range and set 1D basis functions.
c-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_bases: pdmax>pd")
      xstd=2._r8*x-1._r8
      ystd=2._r8*y-1._r8
      DO i=0,pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
c-----------------------------------------------------------------------
c     if derivatives are not needed, limit the computations.
c-----------------------------------------------------------------------
      IF (dmode==0) THEN
c-----------------------------------------------------------------------
c       compute basis functions with respect to the logical coordinates
c       of the unit square.  the basis ordering starts from pdmin in y
c       with increasing degree in x from 0 to pd then completes from
c       pdmin in x and increasing degree in y.
c-----------------------------------------------------------------------
        k=1
        DO j=pdmin,pdmax
          DO i=0,pd
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
        ENDDO
        DO i=pdmin,pdmax
          DO j=0,pdmin-1
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
          DO j=pdmax+1,pd
            alpha(k)=alx(i)*aly(j) 
            k=k+1
          ENDDO
        ENDDO
      ELSE
c-----------------------------------------------------------------------
c       compute basis functions and derivatives.
c-----------------------------------------------------------------------
        k=1
        DO j=pdmin,pdmax
          DO i=0,pd
            alpha(k)=  alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1
          ENDDO
        ENDDO
        DO i=pdmin,pdmax
          DO j=0,pdmin-1
            alpha(k)=  alx(i)* aly(j) 
            alphax(k)=dalx(i)* aly(j) 
            alphay(k)= alx(i)*daly(j) 
            k=k+1
          ENDDO
          DO j=pdmax+1,pd
            alpha(k)=  alx(i)* aly(j) 
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
      END SUBROUTINE modal_disc_bases
c-----------------------------------------------------------------------
c     subprogram 2. modal_disc_3D_alloc.
c     allocates space for modal_disc_3D_type.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_alloc(modq,mx,my,nqty,nfour,pd,pdmin,
     $                               pdmax,name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,pd,pdmin,pdmax
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(modal_quad_type), INTENT(OUT) :: modq
c-----------------------------------------------------------------------
c     store grid, vector, and fourier series dimensions, and set the
c     number of interior basis functions.
c-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_alloc: pdmax>pd")
      modq%mx=mx
      modq%my=my
      modq%nqty=nqty
      modq%nfour=nfour
      modq%pd=pd
      modq%pdmin=pdmin
      modq%pdmax=pdmax
      modq%n_int=(pdmax-pdmin+1)*(2*pd+pdmin-pdmax+1)
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(modq%fsi(nqty,modq%n_int,1:mx,1:my,nfour))
      ALLOCATE(modq%title(nqty))
      ALLOCATE(modq%f(nqty,nfour))
      ALLOCATE(modq%fx(nqty,nfour))
      ALLOCATE(modq%fy(nqty,nfour))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) modq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          modq%title=title
        ELSE
          modq%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     all discontinuous bases are element-centered, so ix0=iy0=1.
c     for modal elements, their are no nodal positions (dx and dy in
c     nimrod's nodal representations).
c-----------------------------------------------------------------------
      ALLOCATE(modq%ix0(modq%n_int))
      ALLOCATE(modq%iy0(modq%n_int))
      modq%ix0=1
      modq%iy0=1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_alloc
c-----------------------------------------------------------------------
c     subprogram 3. modal_disc_3D_dealloc.
c     deallocates space for modal_disc_3D_type.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_dealloc(modq)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(modq%fsi)
      DEALLOCATE(modq%title)
      DEALLOCATE(modq%f)
      DEALLOCATE(modq%fx)
      DEALLOCATE(modq%fy)
      DEALLOCATE(modq%ix0)
      DEALLOCATE(modq%iy0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_dealloc
c-----------------------------------------------------------------------
c     subprogram 4. modal_disc_3D_eval.
c     evaluates complex modal_disc quantities at a single point within a
c     grid block.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_eval(modq,x,y,dmode)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,i,j,k,im
      REAL(r8) :: xstd,ystd
c-----------------------------------------------------------------------
c     find the interval, and compute 1D basis coefficients with logical
c     coordinates for a standard polynomial in -1<=xstd,ystd<=+1.
c-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),modq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),modq%my-1),0_i4)
      xstd=2._r8*(x-ix)-1._r8
      ystd=2._r8*(y-iy)-1._r8
      DO i=0,modq%pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  it is
c     important that the basis-function order is the same as defined
c     in modal_disc_bases.
c-----------------------------------------------------------------------
      k=1
      modq%f=0._r8
      modq%fx=0._r8
      modq%fy=0._r8

      IF (dmode==0) THEN
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1,:)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
      ELSE
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1,:)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1,:)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1,:)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_eval
c-----------------------------------------------------------------------
c     subprogram 5. modal_disc_3D_all_eval.
c     evaluates complex modal_disc quantities in all elements in a grid
c     block for equal spacing.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_all_eval(modq,x,y,f,fx,fy,dmode)

      TYPE(modal_quad_type), INTENT(IN) :: modq
      REAL(r8), INTENT(IN) :: x,y
      COMPLEX(r8), INTENT(OUT), DIMENSION(:,:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: i,j,k,im,ix,iy
      REAL(r8), DIMENSION(modq%n_int) :: alpha,dalpdx,dalpdy
c-----------------------------------------------------------------------
c     evaluate the bases at the specified offset (x,y) within each
c     element.
c-----------------------------------------------------------------------
      CALL modal_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode,
     $                      modq%pd,modq%pdmin,modq%pdmax)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.
c-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO im=1,modq%nfour
          DO iy=1,modq%my
            DO ix=1,modq%mx
              DO i=1,modq%nqty
                f(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*alpha(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO im=1,modq%nfour
          DO iy=1,modq%my
            DO ix=1,modq%mx
              DO i=1,modq%nqty
                f(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*alpha(:))
                fx(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*dalpdx(:))
                fy(i,ix,iy,im)=SUM(modq%fsi(i,:,ix,iy,im)*dalpdy(:))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_all_eval
c-----------------------------------------------------------------------
c     subprogram 6. modal_disc_3D_assign_rsc.
c     assign a real scalar value to a complex modalange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_rsc(modq,rscalar)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: rscalar

      modq%fsi=rscalar
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 7. modal_disc_3D_assign_csc.
c     assign a complex scalar value to a complex modalange quad
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_csc(modq,cscalar)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), INTENT(IN) :: cscalar

      modq%fsi=cscalar
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 8. modal_disc_3D_assign_modq.
c     set one complex modalange quad structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_modq(modq1,modq2)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq1
      TYPE(modal_quad_type), INTENT(IN) :: modq2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      modq1%fsi(:,:,:,:,:)=modq2%fsi(:,:,:,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_modq
c-----------------------------------------------------------------------
c     subprogram 9. modal_disc_3D_assign_int.
c     assign a integer value to a complex modalange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_assign_int(modq,int)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      INTEGER(i4), INTENT(IN) :: int

      modq%fsi=int
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_assign_int
c-----------------------------------------------------------------------
c     subprogram 10. modal_disc_3D_basis_assign_arr
c     assign data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_assign_arr(modq,data,ibasis)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:,:)=data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_assign_arr
c-----------------------------------------------------------------------
c     subprogram 11. modal_disc_3D_basis_add_arr
c     add data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_add_arr(modq,data,ibasis)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:,:)=modq%fsi(:,ibasis,:,:,:)+data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_add_arr
c-----------------------------------------------------------------------
c     subprogram 12. modal_disc_3D_basis_assign_loc
c     assign data into coefficient arrays for one basis function.
c
c     this is a local version of modal_disc_basis_assign_arr, where the
c     the data is located at a given poloidal and fourier indices
c     triplet, only.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_assign_loc(modq,data,ibasis,ix,iy,
     $                                          im)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy,im)=data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_assign_loc
c-----------------------------------------------------------------------
c     subprogram 13. modal_disc_3D_basis_add_loc
c     add data into coefficient arrays for one basis function.
c
c     this is a local version of modal_disc_basis_add_arr, where the
c     the data is located at a given poloidal and fourier indices
c     triplet, only.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_3D_basis_add_loc(modq,data,ibasis,ix,iy,im)

      TYPE(modal_quad_type), INTENT(INOUT) :: modq
      COMPLEX(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy,im
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy,im)=modq%fsi(:,ibasis,ix,iy,im)+data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_3D_basis_add_loc
c-----------------------------------------------------------------------
c     subprogram 14. modal_disc_2D_alloc.
c     allocates space for modal_disc_2D_type.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_alloc(modq,mx,my,nqty,pd,pdmin,pdmax,
     $                               name,title)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,pd,pdmin,pdmax
      CHARACTER(*), INTENT(IN), OPTIONAL :: name
      CHARACTER(*), DIMENSION(:), INTENT(IN), OPTIONAL :: title
      TYPE(modal_quad_2D_type), INTENT(OUT) :: modq
c-----------------------------------------------------------------------
c     store grid and vector dimensions, and set the
c     number of interior basis functions.
c-----------------------------------------------------------------------
      IF (pdmax>pd) CALL nim_stop("Modal_disc_alloc: pdmax>pd")
      modq%mx=mx
      modq%my=my
      modq%nqty=nqty
      modq%pd=pd
      modq%pdmin=pdmin
      modq%pdmax=pdmax
      modq%n_int=(pdmax-pdmin+1)*(2*pd+pdmin-pdmax+1)
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      ALLOCATE(modq%fsi(nqty,modq%n_int,1:mx,1:my))
      ALLOCATE(modq%title(nqty))
      ALLOCATE(modq%f(nqty))
      ALLOCATE(modq%fx(nqty))
      ALLOCATE(modq%fy(nqty))
c-----------------------------------------------------------------------
c     character descriptors, if present in input.
c-----------------------------------------------------------------------
      IF (PRESENT(name)) modq%name=name
      IF (PRESENT(title)) THEN
        IF (SIZE(title)==nqty) THEN
          modq%title=title
        ELSE
          modq%title=title(1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     all discontinuous bases are element-centered, so ix0=iy0=1.
c     for modal elements, their are no nodal positions (dx and dy in
c     nimrod's nodal representations).
c-----------------------------------------------------------------------
      ALLOCATE(modq%ix0(modq%n_int))
      ALLOCATE(modq%iy0(modq%n_int))
      modq%ix0=1
      modq%iy0=1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_alloc
c-----------------------------------------------------------------------
c     subprogram 15. modal_disc_2D_dealloc.
c     deallocates space for modal_disc_2D_type.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_dealloc(modq)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
c-----------------------------------------------------------------------
c     deallocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(modq%fsi)
      DEALLOCATE(modq%title)
      DEALLOCATE(modq%f)
      DEALLOCATE(modq%fx)
      DEALLOCATE(modq%fy)
      DEALLOCATE(modq%ix0)
      DEALLOCATE(modq%iy0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_dealloc
c-----------------------------------------------------------------------
c     subprogram 16. modal_disc_2D_eval.
c     evaluates complex modal_disc quantities at a single point within a
c     grid block.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_eval(modq,x,y,dmode)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: ix,iy,i,j,k,im
      REAL(r8) :: xstd,ystd
c-----------------------------------------------------------------------
c     find the interval, and compute 1D basis coefficients with logical
c     coordinates for a standard polynomial in -1<=xstd,ystd<=+1.
c-----------------------------------------------------------------------
      ix=MAX(MIN(INT(x),modq%mx-1),0_i4)
      iy=MAX(MIN(INT(y),modq%my-1),0_i4)
      xstd=2._r8*(x-ix)-1._r8
      ystd=2._r8*(y-iy)-1._r8
      DO i=0,modq%pd
        CALL leg_poly1_sub(xstd,i,alx(i),dalx(i))
        CALL leg_poly1_sub(ystd,i,aly(i),daly(i))
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.  it is
c     important that the basis-function order is the same as defined
c     in modal_disc_bases.
c-----------------------------------------------------------------------
      k=1
      modq%f=0._r8
      modq%fx=0._r8
      modq%fy=0._r8

      IF (dmode==0) THEN
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f=modq%f+modq%fsi(:,k,ix+1,iy+1)*alx(i)*aly(j)
            k=k+1
          ENDDO
        ENDDO
      ELSE
        DO j=modq%pdmin,modq%pdmax
          DO i=0,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
        DO i=modq%pdmin,modq%pdmax
          DO j=0,modq%pdmin-1
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
          DO j=modq%pdmax+1,modq%pd
            modq%f =modq%f +modq%fsi(:,k,ix+1,iy+1)* alx(i)* aly(j)
            modq%fx=modq%fx+modq%fsi(:,k,ix+1,iy+1)*dalx(i)* aly(j)
            modq%fy=modq%fy+modq%fsi(:,k,ix+1,iy+1)* alx(i)*daly(j)
            k=k+1
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_eval
c-----------------------------------------------------------------------
c     subprogram 17. modal_disc_2D_all_eval.
c     evaluates complex modal_disc quantities in all elements in a grid
c     block for equal spacing.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_all_eval(modq,x,y,f,fx,fy,dmode)

      TYPE(modal_quad_2D_type), INTENT(IN) :: modq
      REAL(r8), INTENT(IN) :: x,y
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) :: f,fx,fy
      INTEGER(i4), INTENT(IN) :: dmode

      INTEGER(i4) :: i,j,k,ix,iy
      REAL(r8), DIMENSION(modq%n_int) :: alpha,dalpdx,dalpdy
c-----------------------------------------------------------------------
c     evaluate the bases at the specified offset (x,y) within each
c     element.
c-----------------------------------------------------------------------
      CALL modal_disc_bases(x,y,alpha,dalpdx,dalpdy,dmode,
     $                      modq%pd,modq%pdmin,modq%pdmax)
c-----------------------------------------------------------------------
c     evaluate the function up to the requested derivative.
c-----------------------------------------------------------------------
      IF (dmode==0) THEN
        DO iy=1,modq%my
          DO ix=1,modq%mx
            DO i=1,modq%nqty
              f(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*alpha(:))
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO iy=1,modq%my
          DO ix=1,modq%mx
            DO i=1,modq%nqty
              f(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*alpha(:))
              fx(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*dalpdx(:))
              fy(i,ix,iy)=SUM(modq%fsi(i,:,ix,iy)*dalpdy(:))
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_all_eval
c-----------------------------------------------------------------------
c     subprogram 18. modal_disc_2D_assign_rsc.
c     assign a real scalar value to a real modalange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_rsc(modq,rscalar)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), INTENT(IN) :: rscalar

      modq%fsi=rscalar
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 19. modal_disc_2D_assign_csc.
c     assign a real scalar value to a real modalange quad
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_csc(modq,cscalar)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      COMPLEX(r8), INTENT(IN) :: cscalar

      modq%fsi=cscalar
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 20. modal_disc_2D_assign_modq.
c     set one real modalange quad structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_modq(modq1,modq2)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq1
      TYPE(modal_quad_2D_type), INTENT(IN) :: modq2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      modq1%fsi(:,:,:,:)=modq2%fsi(:,:,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_modq
c-----------------------------------------------------------------------
c     subprogram 21. modal_disc_2D_assign_int.
c     assign a integer value to a real modalange quad structure.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_assign_int(modq,int)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      INTEGER(i4), INTENT(IN) :: int

      modq%fsi=int
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_assign_int
c-----------------------------------------------------------------------
c     subprogram 22. modal_disc_2D_basis_assign_arr
c     assign data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_assign_arr(modq,data,ibasis)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:)=data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_assign_arr
c-----------------------------------------------------------------------
c     subprogram 23. modal_disc_2D_basis_add_arr
c     add data into coefficient arrays for one basis function.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_add_arr(modq,data,ibasis)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:,:,:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,:,:)=modq%fsi(:,ibasis,:,:)+data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_add_arr
c-----------------------------------------------------------------------
c     subprogram 24. modal_disc_2D_basis_assign_loc
c     assign data into coefficient arrays for one basis function.
c
c     this is a local version of modal_disc_basis_assign_arr, where the
c     the data is located at given poloidal indices, only.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_assign_loc(modq,data,ibasis,ix,iy)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy)=data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_assign_loc
c-----------------------------------------------------------------------
c     subprogram 25. modal_disc_2D_basis_add_loc
c     add data into coefficient arrays for one basis function.
c
c     this is a local version of modal_disc_basis_add_arr, where the
c     the data is located at given poloidal indices only.
c-----------------------------------------------------------------------
      SUBROUTINE modal_disc_2D_basis_add_loc(modq,data,ibasis,ix,iy)

      TYPE(modal_quad_2D_type), INTENT(INOUT) :: modq
      REAL(r8), DIMENSION(:), INTENT(IN) :: data
      INTEGER(i4), INTENT(IN) :: ibasis,ix,iy
c-----------------------------------------------------------------------
c     decide the proper storage location for this data according to
c     the value of ibasis.
c-----------------------------------------------------------------------
      modq%fsi(:,ibasis,ix,iy)=modq%fsi(:,ibasis,ix,iy)+data
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE modal_disc_2D_basis_add_loc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE modal_disc_mod
