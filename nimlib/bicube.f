c-----------------------------------------------------------------------
c     file bicube.f.
c     fits functions to bicubic splines.
c     Reference: H. Spaeth, "Spline Algorithms for Curves and Surfaces,"
c     Translated from the German by W. D. Hoskins and H. W. Sager.
c     Utilitas Mathematica Publishing Inc., Winnepeg, 1974.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c      0. bicube_type definition.
c      1. bicube_alloc.
c      2. bicube_dealloc.
c      3. bicube_fit.
c      4. bicube_eval.
c      5. bicube_getco.
c      6. bicube_all_eval.
c      7. bicube_all_getco.
c      8. bicube_write_xy.
c      9. bicube_write_yx.
c     10. bicube_write_arrays.
c     11. bicube_copy.
c     12. bicube_assign_rsc.
c     13. bicube_assign_bc.
c     14. bicube_assign_int.
c-----------------------------------------------------------------------
c     subprogram 0. bicube_type definition.
c     defines bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE bicube
      USE local
      USE spline
      IMPLICIT NONE

      TYPE :: bicube_type
      INTEGER(i4) :: mx,my,nqty,ix,iy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xs,ys
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: fs,fsx,fsy,fsxy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: f,fx,fy,fxx,fxy,fyy
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: cmats
      CHARACTER(6), DIMENSION(:), ALLOCATABLE :: title
      CHARACTER(6) :: name
      END TYPE bicube_type

c-----------------------------------------------------------------------
c     overload asignment and operators.
c-----------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE bicube_assign_rsc,bicube_assign_bc,
     $    bicube_assign_int
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. bicube_alloc.
c     allocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_alloc(bcs,mx,my,nqty)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty
      TYPE(bicube_type), INTENT(OUT) :: bcs
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      bcs%mx=mx
      bcs%my=my
      bcs%ix=0
      bcs%iy=0
      bcs%nqty=nqty
      ALLOCATE(bcs%xs(0:mx))
      ALLOCATE(bcs%ys(0:my))
      ALLOCATE(bcs%fs(nqty,0:mx,0:my))
      ALLOCATE(bcs%fsx(nqty,0:mx,0:my))
      ALLOCATE(bcs%fsy(nqty,0:mx,0:my))
      ALLOCATE(bcs%fsxy(nqty,0:mx,0:my))
      ALLOCATE(bcs%title(nqty))
      ALLOCATE(bcs%f(nqty))
      ALLOCATE(bcs%fx(nqty))
      ALLOCATE(bcs%fy(nqty))
      ALLOCATE(bcs%fxx(nqty))
      ALLOCATE(bcs%fxy(nqty))
      ALLOCATE(bcs%fyy(nqty))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_alloc
c-----------------------------------------------------------------------
c     subprogram 2. bicube_dealloc.
c     allocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_dealloc(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      DEALLOCATE(bcs%xs)
      DEALLOCATE(bcs%ys)
      DEALLOCATE(bcs%fs)
      DEALLOCATE(bcs%fsx)
      DEALLOCATE(bcs%fsy)
      DEALLOCATE(bcs%fsxy)
      DEALLOCATE(bcs%title)
      DEALLOCATE(bcs%f)
      DEALLOCATE(bcs%fx)
      DEALLOCATE(bcs%fy)
      DEALLOCATE(bcs%fxx)
      DEALLOCATE(bcs%fxy)
      DEALLOCATE(bcs%fyy)
      IF(ALLOCATED(bcs%cmats))DEALLOCATE(bcs%cmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. bicube_fit.
c     fits functions to bicubic splines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_fit(bcs,endmode1,endmode2)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      CHARACTER(*), INTENT(IN) :: endmode1,endmode2

      INTEGER(i4) :: iqty
      TYPE(spline_type) :: spl
c-----------------------------------------------------------------------
c     evaluate y derivatives.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%my,bcs%mx+1_i4)
      spl%xs=bcs%ys
      DO iqty=1,bcs%nqty
         spl%fs=TRANSPOSE(bcs%fs(iqty,:,:))
         CALL spline_fit(spl,endmode2)
         bcs%fsy(iqty,:,:)=TRANSPOSE(spl%fs1)
      ENDDO
      CALL spline_dealloc(spl)
c-----------------------------------------------------------------------
c     evaluate x derivatives.
c-----------------------------------------------------------------------
      CALL spline_alloc(spl,bcs%mx,bcs%my+1_i4)
      spl%xs=bcs%xs
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fs(iqty,:,:)
         CALL spline_fit(spl,endmode1)
         bcs%fsx(iqty,:,:)=spl%fs1
      ENDDO
c-----------------------------------------------------------------------
c     evaluate mixed derivatives.
c-----------------------------------------------------------------------
      DO iqty=1,bcs%nqty
         spl%fs=bcs%fsy(iqty,:,:)
         CALL spline_fit(spl,endmode1)
         bcs%fsxy(iqty,:,:)=spl%fs1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL spline_dealloc(spl)
      RETURN
      END SUBROUTINE bicube_fit
c-----------------------------------------------------------------------
c     subprogram 4. bicube_eval.
c     evaluates bicubic spline function.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_eval(bcs,x,y,mode)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      REAL(r8), INTENT(IN) :: x,y
      INTEGER(i4), INTENT(IN) :: mode

      INTEGER(i4) :: i
      REAL(r8) :: dx,dy
      REAL(r8), DIMENSION (4,4,bcs%nqty) :: c
c-----------------------------------------------------------------------
c     find x interval.
c-----------------------------------------------------------------------
      bcs%ix=max(bcs%ix,0_i4)
      bcs%ix=min(bcs%ix,bcs%mx-1_i4)
      DO
         IF(x >= bcs%xs(bcs%ix) .OR. bcs%ix <= 0)EXIT
         bcs%ix=bcs%ix-1
      ENDDO
      DO
         IF(x < bcs%xs(bcs%ix+1) .OR. bcs%ix >= bcs%mx-1)EXIT
         bcs%ix=bcs%ix+1
      ENDDO
c-----------------------------------------------------------------------
c     find y interval.
c-----------------------------------------------------------------------
      bcs%iy=max(bcs%iy,0_i4)
      bcs%iy=min(bcs%iy,bcs%my-1_i4)
      DO
         IF(y >= bcs%ys(bcs%iy) .OR. bcs%iy <= 0)EXIT
         bcs%iy=bcs%iy-1
      ENDDO
      DO
         IF(y < bcs%ys(bcs%iy+1) .OR. bcs%iy >= bcs%my-1)EXIT
         bcs%iy=bcs%iy+1
      ENDDO
c-----------------------------------------------------------------------
c     find offsets and compute local coefficients.
c-----------------------------------------------------------------------
      dx=x-bcs%xs(bcs%ix)
      dy=y-bcs%ys(bcs%iy)
      IF(ALLOCATED(bcs%cmats))THEN
         c=bcs%cmats(:,:,:,bcs%ix+1,bcs%iy+1)
      ELSE
         c=bicube_getco(bcs)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      bcs%f=0
      DO i=4,1,-1
         bcs%f=bcs%f*dx
     $        +((c(i,4,:)*dy
     $        +c(i,3,:))*dy
     $        +c(i,2,:))*dy
     $        +c(i,1,:)
      ENDDO
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     evaluate first derivatives of f
c-----------------------------------------------------------------------
      bcs%fx=0
      bcs%fy=0
      DO i=4,1,-1
         bcs%fy=bcs%fy*dx
     $        +(c(i,4,:)*3*dy
     $        +c(i,3,:)*2)*dy
     $        +c(i,2,:)
         bcs%fx=bcs%fx*dy
     $        +(c(4,i,:)*3*dx
     $        +c(3,i,:)*2)*dx
     $        +c(2,i,:)
      ENDDO
      IF(mode == 1)RETURN
c-----------------------------------------------------------------------
c     evaluate second derivatives of f
c-----------------------------------------------------------------------
      bcs%fxx=0
      bcs%fyy=0
      bcs%fxy=0
      DO i=4,1,-1
         bcs%fyy=bcs%fyy*dx
     $        +(c(i,4,:)*3*dy
     $        +c(i,3,:))*2
         bcs%fxx=bcs%fxx*dy
     $        +(c(4,i,:)*3*dx
     $        +c(3,i,:))*2
      ENDDO
      DO i=4,2,-1
         bcs%fxy=bcs%fxy*dx
     $        +((c(i,4,:)*3*dy
     $        +c(i,3,:)*2)*dy
     $        +c(i,2,:))*(i-1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_eval
c-----------------------------------------------------------------------
c     subprogram 5. bicube_getco.
c     computes coefficient matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      FUNCTION bicube_getco(bcs) RESULT(cmat)

      TYPE(bicube_type), INTENT(IN) :: bcs
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: cmat

      REAL(r8) :: hxfac,hxfac2,hxfac3
      REAL(r8) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4) :: gxmat,gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty) :: temp
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(bcs%ix+1)-bcs%xs(bcs%ix))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1)=-3*hxfac2
      gxmat(3,2)=-2*hxfac
      gxmat(3,3)=3*hxfac2
      gxmat(3,4)=-hxfac
      gxmat(4,1)=2*hxfac3
      gxmat(4,2)=hxfac2
      gxmat(4,3)=-2*hxfac3
      gxmat(4,4)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(bcs%iy+1)-bcs%ys(bcs%iy))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1)=-3*hyfac2
      gymat(3,2)=-2*hyfac
      gymat(3,3)=3*hyfac2
      gymat(3,4)=-hyfac
      gymat(4,1)=2*hyfac3
      gymat(4,2)=hyfac2
      gymat(4,3)=-2*hyfac3
      gymat(4,4)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      cmat(1,1,:)=bcs%fs(:,bcs%ix,bcs%iy)
      cmat(1,2,:)=bcs%fsy(:,bcs%ix,bcs%iy)
      cmat(1,3,:)=bcs%fs(:,bcs%ix,bcs%iy+1)
      cmat(1,4,:)=bcs%fsy(:,bcs%ix,bcs%iy+1)
      cmat(2,1,:)=bcs%fsx(:,bcs%ix,bcs%iy)
      cmat(2,2,:)=bcs%fsxy(:,bcs%ix,bcs%iy)
      cmat(2,3,:)=bcs%fsx(:,bcs%ix,bcs%iy+1)
      cmat(2,4,:)=bcs%fsxy(:,bcs%ix,bcs%iy+1)
      cmat(3,1,:)=bcs%fs(:,bcs%ix+1,bcs%iy)
      cmat(3,2,:)=bcs%fsy(:,bcs%ix+1,bcs%iy)
      cmat(3,3,:)=bcs%fs(:,bcs%ix+1,bcs%iy+1)
      cmat(3,4,:)=bcs%fsy(:,bcs%ix+1,bcs%iy+1)
      cmat(4,1,:)=bcs%fsx(:,bcs%ix+1,bcs%iy)
      cmat(4,2,:)=bcs%fsxy(:,bcs%ix+1,bcs%iy)
      cmat(4,3,:)=bcs%fsx(:,bcs%ix+1,bcs%iy+1)
      cmat(4,4,:)=bcs%fsxy(:,bcs%ix+1,bcs%iy+1)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:)=cmat(:,1:2,:)
      temp(:,3,:)
     $     =cmat(:,1,:)*gymat(3,1)
     $     +cmat(:,2,:)*gymat(3,2)
     $     +cmat(:,3,:)*gymat(3,3)
     $     +cmat(:,4,:)*gymat(3,4)
      temp(:,4,:)
     $     =cmat(:,1,:)*gymat(4,1)
     $     +cmat(:,2,:)*gymat(4,2)
     $     +cmat(:,3,:)*gymat(4,3)
     $     +cmat(:,4,:)*gymat(4,4)
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      cmat(1:2,:,:)=temp(1:2,:,:)
      cmat(3,:,:)
     $     =gxmat(3,1)*temp(1,:,:)
     $     +gxmat(3,2)*temp(2,:,:)
     $     +gxmat(3,3)*temp(3,:,:)
     $     +gxmat(3,4)*temp(4,:,:)
      cmat(4,:,:)
     $     =gxmat(4,1)*temp(1,:,:)
     $     +gxmat(4,2)*temp(2,:,:)
     $     +gxmat(4,3)*temp(3,:,:)
     $     +gxmat(4,4)*temp(4,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION bicube_getco
c-----------------------------------------------------------------------
c     subprogram 6. bicube_all_eval.
c     evaluates bicubic splines in all intervals for equal spacing. 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_all_eval(bcs,dx,dy,f,fx,fy,fxx,fyy,fxy,mode)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      REAL(r8), INTENT(IN) :: dx,dy
      REAL(r8), INTENT(OUT), DIMENSION(:,:,:) ::
     $     f,fx,fy,fxx,fyy,fxy
      INTEGER(i4), INTENT(IN) :: mode

      INTEGER(i4) :: i,ix,iy
      REAL(r8), DIMENSION(bcs%mx) :: dxv
      REAL(r8), DIMENSION(bcs%my) :: dyv
c-----------------------------------------------------------------------
c     compute local displacements and coefficients.
c-----------------------------------------------------------------------
      dxv=(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1))*dx
      dyv=(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1))*dy
      CALL bicube_all_getco(bcs)
c-----------------------------------------------------------------------
c     evaluate f.
c-----------------------------------------------------------------------
      f=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO ix=1,bcs%mx
               f(:,ix,:)=f(:,ix,:)*dxv(ix)
            ENDDO
         ENDIF
         DO iy=1,bcs%my
            f(:,:,iy)=f(:,:,iy)
     $           +((bcs%cmats(i,4,:,:,iy)*dyv(iy)
     $            +bcs%cmats(i,3,:,:,iy))*dyv(iy)
     $            +bcs%cmats(i,2,:,:,iy))*dy
     $            +bcs%cmats(i,1,:,:,iy)
         ENDDO
      ENDDO
      IF(mode == 0)RETURN
c-----------------------------------------------------------------------
c     evaluate fx.
c-----------------------------------------------------------------------
      fx=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO iy=1,bcs%my
               fx(:,:,iy)=fx(:,:,iy)*dyv(iy)
            ENDDO
         ENDIF
         DO ix=1,bcs%mx
            fx(:,ix,:)=fx(:,ix,:)
     $           +(bcs%cmats(4,i,:,ix,:)*3*dxv(ix)
     $           +bcs%cmats(3,i,:,ix,:)*2)*dxv(ix)
     $           +bcs%cmats(2,i,:,ix,:)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     evaluate fy.
c-----------------------------------------------------------------------
      fy=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO ix=1,bcs%mx
               fy(:,ix,:)=fy(:,ix,:)*dxv(ix)
            ENDDO
         ENDIF
         DO iy=1,bcs%my
            fy(:,:,iy)=fy(:,:,iy)
     $           +(bcs%cmats(i,4,:,:,iy)*3*dyv(iy)
     $           +bcs%cmats(i,3,:,:,iy)*2)*dyv(iy)
     $           +bcs%cmats(i,2,:,:,iy)
         ENDDO
      ENDDO
      IF(mode == 1)RETURN
c-----------------------------------------------------------------------
c     evaluate fxx.
c-----------------------------------------------------------------------
      fxx=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO iy=1,bcs%my
               fxx(:,:,iy)=fxx(:,:,iy)*dyv(iy)
            ENDDO
         ENDIF
         DO ix=1,bcs%mx
            fxx(:,ix,:)=fxx(:,ix,:)
     $           +(bcs%cmats(4,i,:,ix,:)*3*dxv(ix)
     $           +bcs%cmats(3,i,:,ix,:))*2
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     evaluate fyy and fxy
c-----------------------------------------------------------------------
      fyy=0
      DO i=4,1,-1
         IF(i /= 4)THEN
            DO ix=1,bcs%mx
               fyy(:,ix,:)=fyy(:,ix,:)*dxv(ix)
               fxy(:,ix,:)=fxy(:,ix,:)*dxv(ix)
            ENDDO
         ENDIF
         DO iy=1,bcs%my
            fyy(:,:,iy)=fyy(:,:,iy)
     $           +(bcs%cmats(i,4,:,:,iy)*3*dyv(iy)
     $           +bcs%cmats(i,3,:,:,iy))*2
            fxy(:,:,iy)=fxy(:,:,iy)
     $           +((bcs%cmats(i,4,:,:,iy)*3*dyv(iy)
     $           +bcs%cmats(i,3,:,:,iy)*2)*dyv(iy)
     $           +bcs%cmats(i,2,:,:,iy))*(i-1)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_all_eval
c-----------------------------------------------------------------------
c     subprogram 7. bicube_all_getco.
c     computes coefficient matrices.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_all_getco(bcs)

      TYPE(bicube_type), INTENT(INOUT) :: bcs

      INTEGER(i4) :: ix,iy
      REAL(r8), DIMENSION(bcs%mx) :: hxfac,hxfac2,hxfac3
      REAL(r8), DIMENSION(bcs%my) :: hyfac,hyfac2,hyfac3
      REAL(r8), DIMENSION(3:4,4,bcs%mx) :: gxmat
      REAL(r8), DIMENSION(3:4,4,bcs%my) :: gymat
      REAL(r8), DIMENSION(4,4,bcs%nqty,bcs%mx,bcs%my) :: temp
c-----------------------------------------------------------------------
c     allocate space.
c-----------------------------------------------------------------------
      IF(ALLOCATED(bcs%cmats))THEN
         RETURN
      ELSE
         ALLOCATE(bcs%cmats(4,4,bcs%nqty,bcs%mx,bcs%my))
      ENDIF
c-----------------------------------------------------------------------
c     compute gxmat.
c-----------------------------------------------------------------------
      hxfac=1/(bcs%xs(1:bcs%mx)-bcs%xs(0:bcs%mx-1))
      hxfac2=hxfac*hxfac
      hxfac3=hxfac2*hxfac
      gxmat(3,1,:)=-3*hxfac2
      gxmat(3,2,:)=-2*hxfac
      gxmat(3,3,:)=3*hxfac2
      gxmat(3,4,:)=-hxfac
      gxmat(4,1,:)=2*hxfac3
      gxmat(4,2,:)=hxfac2
      gxmat(4,3,:)=-2*hxfac3
      gxmat(4,4,:)=hxfac2
c-----------------------------------------------------------------------
c     compute gymat.
c-----------------------------------------------------------------------
      hyfac=1/(bcs%ys(1:bcs%my)-bcs%ys(0:bcs%my-1))
      hyfac2=hyfac*hyfac
      hyfac3=hyfac2*hyfac
      gymat(3,1,:)=-3*hyfac2
      gymat(3,2,:)=-2*hyfac
      gymat(3,3,:)=3*hyfac2
      gymat(3,4,:)=-hyfac
      gymat(4,1,:)=2*hyfac3
      gymat(4,2,:)=hyfac2
      gymat(4,3,:)=-2*hyfac3
      gymat(4,4,:)=hyfac2
c-----------------------------------------------------------------------
c     compute smat.
c-----------------------------------------------------------------------
      bcs%cmats(1,1,:,:,:)=bcs%fs(:,0:bcs%mx-1,0:bcs%my-1)
      bcs%cmats(1,2,:,:,:)=bcs%fsy(:,0:bcs%mx-1,0:bcs%my-1)
      bcs%cmats(1,3,:,:,:)=bcs%fs(:,0:bcs%mx-1,1:bcs%my)
      bcs%cmats(1,4,:,:,:)=bcs%fsy(:,0:bcs%mx-1,1:bcs%my)
      bcs%cmats(2,1,:,:,:)=bcs%fsx(:,0:bcs%mx-1,0:bcs%my-1)
      bcs%cmats(2,2,:,:,:)=bcs%fsxy(:,0:bcs%mx-1,0:bcs%my-1)
      bcs%cmats(2,3,:,:,:)=bcs%fsx(:,0:bcs%mx-1,1:bcs%my)
      bcs%cmats(2,4,:,:,:)=bcs%fsxy(:,0:bcs%mx-1,1:bcs%my)
      bcs%cmats(3,1,:,:,:)=bcs%fs(:,1:bcs%mx,0:bcs%my-1)
      bcs%cmats(3,2,:,:,:)=bcs%fsy(:,1:bcs%mx,0:bcs%my-1)
      bcs%cmats(3,3,:,:,:)=bcs%fs(:,1:bcs%mx,1:bcs%my)
      bcs%cmats(3,4,:,:,:)=bcs%fsy(:,1:bcs%mx,1:bcs%my)
      bcs%cmats(4,1,:,:,:)=bcs%fsx(:,1:bcs%mx,0:bcs%my-1)
      bcs%cmats(4,2,:,:,:)=bcs%fsxy(:,1:bcs%mx,0:bcs%my-1)
      bcs%cmats(4,3,:,:,:)=bcs%fsx(:,1:bcs%mx,1:bcs%my)
      bcs%cmats(4,4,:,:,:)=bcs%fsxy(:,1:bcs%mx,1:bcs%my)
c-----------------------------------------------------------------------
c     multiply by gymat^T.
c-----------------------------------------------------------------------
      temp(:,1:2,:,:,:)=bcs%cmats(:,1:2,:,:,:)
      DO iy=1,bcs%my
         temp(:,3,:,:,iy)
     $        =bcs%cmats(:,1,:,:,iy)*gymat(3,1,iy)
     $        +bcs%cmats(:,2,:,:,iy)*gymat(3,2,iy)
     $        +bcs%cmats(:,3,:,:,iy)*gymat(3,3,iy)
     $        +bcs%cmats(:,4,:,:,iy)*gymat(3,4,iy)
         temp(:,4,:,:,iy)
     $        =bcs%cmats(:,1,:,:,iy)*gymat(4,1,iy)
     $        +bcs%cmats(:,2,:,:,iy)*gymat(4,2,iy)
     $        +bcs%cmats(:,3,:,:,iy)*gymat(4,3,iy)
     $        +bcs%cmats(:,4,:,:,iy)*gymat(4,4,iy)
      ENDDO
c-----------------------------------------------------------------------
c     multiply by gxmat.
c-----------------------------------------------------------------------
      bcs%cmats(1:2,:,:,:,:)=temp(1:2,:,:,:,:)
      DO ix=1,bcs%mx
         bcs%cmats(3,:,:,ix,:)
     $        =gxmat(3,1,ix)*temp(1,:,:,ix,:)
     $        +gxmat(3,2,ix)*temp(2,:,:,ix,:)
     $        +gxmat(3,3,ix)*temp(3,:,:,ix,:)
     $        +gxmat(3,4,ix)*temp(4,:,:,ix,:)
         bcs%cmats(4,:,:,ix,:)
     $        =gxmat(4,1,ix)*temp(1,:,:,ix,:)
     $        +gxmat(4,2,ix)*temp(2,:,:,ix,:)
     $        +gxmat(4,3,ix)*temp(3,:,:,ix,:)
     $        +gxmat(4,4,ix)*temp(4,:,:,ix,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_all_getco
c-----------------------------------------------------------------------
c     subprogram 8. bicube_write_xy.
c     produces ascii and binarx output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_xy(bcs,out,bin,iua,iub)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub

      INTEGER(i4) :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      LOGICAL, PARAMETER :: interp=.TRUE.
      CHARACTER(80) :: format2,format1
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,'iy = ',i3,', y = ',1p,e11.3)
 20   FORMAT(1x,'iy = ',i3,', jy = ',i1,', y = ',1p,e11.3)
 30   FORMAT('(/4x,"ix",6x,"x",4x,',i3.3,'(4x,"f(",i3.3,")",1x)/)')
 40   FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,'(1x,a/)')"input data"
         WRITE(format1,30)bcs%nqty
         WRITE(format2,40)bcs%nqty
      ENDIF
      DO iy=0,bcs%my
         y=bcs%ys(iy)
         IF(out)then
            WRITE(iua,10)iy,bcs%ys(iy)
            WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
         ENDIF
         DO ix=0,bcs%mx
            x=bcs%xs(ix)
            bcs%f=bcs%fs(ix,iy,:)
            IF(out)WRITE(iua,format2)ix,x,bcs%f
            IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4)
         ENDDO
         IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
         IF(bin)WRITE(iub)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
      IF(interp)THEN
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data"
         DO iy=0,bcs%my-1
            dy=(bcs%ys(iy+1)-bcs%ys(iy))/4
            DO jy=0,4
               y=bcs%ys(iy)+dy*jy
               IF(out)then
                  WRITE(iua,20)iy,jy,y
                  WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
               DO ix=0,bcs%mx-1
                  dx=(bcs%xs(ix+1)-bcs%xs(ix))/4
                  DO jx=0,4
                     x=bcs%xs(ix)+dx*jx
                     CALL bicube_eval(bcs,y,x,0_i4)
                     IF(out)WRITE(iua,format2)ix,x,bcs%f
                     IF(bin)WRITE(iub)REAL(x,4),REAL(bcs%f,4)
                  ENDDO
                  IF(out)WRITE(iua,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over y.
c-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
               IF(bin)WRITE(iub)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_xy
c-----------------------------------------------------------------------
c     subprogram 9. bicube_write_yx.
c     produces ascii and binary output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_yx(bcs,out,bin,iua,iub)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out,bin
      INTEGER, INTENT(IN) :: iua,iub

      INTEGER(i4) :: ix,iy,jx,jy,iqty
      REAL(r8) :: x,y,dx,dy

      LOGICAL, PARAMETER :: interp=.TRUE.
      CHARACTER(80) :: format1,format2
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(1x,'ix = ',i3,', x = ',1p,e11.3)
 20   FORMAT(1x,'ix = ',i3,', jx = ',i1,', x = ',1p,e11.3)
 30   FORMAT('(/4x,"iy",6x,"y",4x,',i3.3,'(4x,"f(",i3.3,")",1x)/)')
 40   FORMAT('(i6,1p,e11.3,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. (out. OR. bin))RETURN
c-----------------------------------------------------------------------
c     write input data.
c-----------------------------------------------------------------------
      IF(out)THEN
         WRITE(iua,'(1x,a/)')"input data"
         WRITE(format1,30)bcs%nqty
         WRITE(format2,40)bcs%nqty
      ENDIF
      DO ix=0,bcs%mx
         x=bcs%xs(ix)
         IF(out)then
            WRITE(iua,10)ix,bcs%xs(ix)
            WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
         ENDIF
         DO iy=0,bcs%my
            y=bcs%ys(iy)
            bcs%f=bcs%fs(ix,iy,:)
            IF(out)WRITE(iua,format2)iy,y,bcs%f
            IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4)
         ENDDO
         IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
         IF(bin)WRITE(iub)
      ENDDO
c-----------------------------------------------------------------------
c     begin loops over x for interpolated data.
c-----------------------------------------------------------------------
      IF(interp)THEN
         IF(out)WRITE(iua,'(1x,a/)')"interpolated data"
         DO ix=0,bcs%mx-1
            dx=(bcs%xs(ix+1)-bcs%xs(ix))/4
            DO jx=0,4
               x=bcs%xs(ix)+dx*jx
               IF(out)then
                  WRITE(iua,20)ix,jx,x
                  WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
               ENDIF
c-----------------------------------------------------------------------
c     begin loops over y for interpolated data.
c-----------------------------------------------------------------------
               DO iy=0,bcs%my-1
                  dy=(bcs%ys(iy+1)-bcs%ys(iy))/4
                  DO jy=0,4
                     y=bcs%ys(iy)+dy*jy
                     CALL bicube_eval(bcs,x,y,0_i4)
                     IF(out)WRITE(iua,format2)iy,y,bcs%f
                     IF(bin)WRITE(iub)REAL(y,4),REAL(bcs%f,4)
                  ENDDO
                  IF(out)WRITE(iua,'()')
               ENDDO
c-----------------------------------------------------------------------
c     complete loops over x.
c-----------------------------------------------------------------------
               IF(out)WRITE(iua,format1)(iqty,iqty=1,bcs%nqty)
               IF(bin)WRITE(iub)
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_yx
c-----------------------------------------------------------------------
c     subprogram 10. bicube_write_arrays.
c     produces ascii and binary output for bicubic spline fits.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_write_arrays(bcs,out,iua,iqty)

      TYPE(bicube_type), INTENT(INOUT) :: bcs
      LOGICAL, INTENT(IN) :: out
      INTEGER, INTENT(IN) :: iua,iqty

      CHARACTER(80) :: format1,format2
      INTEGER(i4) :: ix,iy
c-----------------------------------------------------------------------
c     formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/2x,"ix/iy",',i3.3,'(3x,i3.3,5x)/)')
 20   FORMAT('(i5,1p,',i3.3,'e11.3)')
c-----------------------------------------------------------------------
c     abort.
c-----------------------------------------------------------------------
      IF(.NOT. out)RETURN
c-----------------------------------------------------------------------
c     write fs.
c-----------------------------------------------------------------------
      WRITE(format1,10)bcs%my+1
      WRITE(format2,20)bcs%my+1
      WRITE(iua,"(a)")"fs"
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,format2)(ix,(bcs%fs(iqty,ix,iy),iy=0,bcs%my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,"(a/)")"fs"
c-----------------------------------------------------------------------
c     write fsx.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsx"
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,format2)(ix,(bcs%fsx(iqty,ix,iy),iy=0,bcs%my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,"(a/)")"fsx"
c-----------------------------------------------------------------------
c     write fsy.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsy"
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,format2)(ix,(bcs%fsy(iqty,ix,iy),iy=0,bcs%my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,"(a/)")"fsy"
c-----------------------------------------------------------------------
c     write fsxy.
c-----------------------------------------------------------------------
      WRITE(iua,"(a)")"fsxy"
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,format2)(ix,(bcs%fsxy(iqty,ix,iy),iy=0,bcs%my),
     $     ix=0,bcs%mx)
      WRITE(iua,format1)(iy,iy=0,bcs%my)
      WRITE(iua,"(a/)")"fsxy"
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_write_arrays
c-----------------------------------------------------------------------
c     subprogram 11. bicube_copy.
c     allocates space for bicube_type.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_copy(bcs1,bcs2)

      TYPE(bicube_type), INTENT(IN) :: bcs1
      TYPE(bicube_type), INTENT(INOUT) :: bcs2
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      IF(ALLOCATED(bcs2%xs))CALL bicube_dealloc(bcs2)
      CALL bicube_alloc(bcs2,bcs1%mx,bcs1%my,bcs1%nqty)
      bcs2%xs=bcs1%xs
      bcs2%ys=bcs1%ys
      bcs2%fs=bcs1%fs
      bcs2%fsx=bcs1%fsx
      bcs2%fsy=bcs1%fsy
      bcs2%fsxy=bcs1%fsxy
      bcs2%name=bcs1%name
      bcs2%title=bcs1%title
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_copy
c-----------------------------------------------------------------------
c     subprogram 12. bicube_assign_rsc.
c     assign a real scalar value to a bicube structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_rsc(bc,rscalar)

      TYPE(bicube_type), INTENT(INOUT) :: bc
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     derivatives are 0.
c-----------------------------------------------------------------------
      bc%fs=rscalar
      bc%fsx=0
      bc%fsy=0
      bc%fsxy=0
      IF (ALLOCATED(bc%cmats)) DEALLOCATE(bc%cmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 13. bicube_assign_bc.
c     set one bicube structure array values equal to another.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_bc(bc1,bc2)

      TYPE(bicube_type), INTENT(INOUT) :: bc1
      TYPE(bicube_type), INTENT(IN) :: bc2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      bc1%fs=bc2%fs
      bc1%fsx=bc2%fsx
      bc1%fsy=bc2%fsy
      bc1%fsxy=bc2%fsxy
      IF (ALLOCATED(bc1%cmats)) DEALLOCATE(bc1%cmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_assign_bc
c-----------------------------------------------------------------------
c     subprogram 14. bicube_assign_int.
c     assign a integer value to a bicube structure.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE bicube_assign_int(bc,int)

      TYPE(bicube_type), INTENT(INOUT) :: bc
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     derivatives are 0.
c-----------------------------------------------------------------------
      bc%fs=int
      bc%fsx=0
      bc%fsy=0
      bc%fsxy=0
      IF (ALLOCATED(bc%cmats)) DEALLOCATE(bc%cmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bicube_assign_int
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE bicube
