c-----------------------------------------------------------------------
c     file tblock.f
c     contains the module for handling finite element computations
c     on unstructured triangular grid-blocks.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  tblock_set
c     2.  tblock_make_real_matrix.
c     3.  tblock_make_comp_matrix.
c     4.  tblock_get_real_rhs.
c     5.  tblock_get_comp_rhs.
c     6.  tblock_get_comp_rhs_q.
c     7.  tblock_basis_set.
c     8.  tblock_basis_dealloc.
c     9.  tblock_real_qp_update.
c     10. tblock_comp_qp_update.
c     11. tblock_real_qp_alloc.
c     12. tblock_comp_qp_alloc.
c     13. tblock_real_qp_dealloc.
c     14. tblock_comp_qp_dealloc.
c-----------------------------------------------------------------------
c     subprogram 0. tblock.
c     module definitions.
c-----------------------------------------------------------------------
      MODULE tblock
      USE local
      USE tblock_type_mod
      IMPLICIT NONE

      INTERFACE tblock_make_matrix
        MODULE PROCEDURE tblock_make_real_matrix,tblock_make_comp_matrix
      END INTERFACE

      INTERFACE tblock_get_rhs
        MODULE PROCEDURE tblock_get_real_rhs,tblock_get_comp_rhs
      END INTERFACE

      INTERFACE tblock_qp_update
        MODULE PROCEDURE tblock_real_qp_update,tblock_comp_qp_update
      END INTERFACE

      INTERFACE tblock_qp_alloc
        MODULE PROCEDURE tblock_real_qp_alloc,tblock_comp_qp_alloc
      END INTERFACE

      INTERFACE tblock_qp_dealloc
        MODULE PROCEDURE tblock_real_qp_dealloc,tblock_comp_qp_dealloc
      END INTERFACE

c-----------------------------------------------------------------------
c-PRE declarations for Gaussian quadrature weights.
c-----------------------------------------------------------------------
      INTEGER(i4), PARAMETER, PRIVATE :: ngauss=7
      REAL(r8), PARAMETER, PRIVATE :: wq0=0.225_r8,
     $     wq1=0.1259391805448271_r8,wq2=0.1323941527885062_r8
      REAL(r8), DIMENSION(ngauss), PARAMETER, PRIVATE :: 
     $     quad_weight=(/wq0,wq1,wq1,wq1,wq2,wq2,wq2/)

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. tblock_set.
c     set the weights for quadratures.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_set(tb)

      TYPE(tblock_type), INTENT(INOUT) :: tb
c-----------------------------------------------------------------------
c-PRE
c     set number of quadrature points and weights according to input.
c     the number of points is now adjusted automatically with
c     poly_degree.
c-----------------------------------------------------------------------
      tb%ng=ngauss     
      ALLOCATE(tb%wg(tb%ng))
      tb%wg=quad_weight
c-----------------------------------------------------------------------
c     terminate tblock_set.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_set
c-----------------------------------------------------------------------
c     subprogram 2. tblock_make_real_matrix.
c     computes a linear response matrix for a supplied integrand
c     subprogram that fits the interface block at the beginning
c     of this subprogram.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_make_real_matrix(tb,mat,get_integrand,nqty)
      USE rblock_type_mod
      USE matrix_type_mod
      USE time

      TYPE(tblock_type), INTENT(INOUT), TARGET :: tb
      TYPE(matrix_element_type3), DIMENSION(0:), INTENT(OUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      TYPE(rblock_type) :: rdum
      REAL(r8) :: dx=0,dy=0

      REAL(r8), DIMENSION(nqty,nqty,tb%mcell,3,3) :: integrand
      REAL(r8), DIMENSION(:,:), POINTER :: bigr

      INTEGER(i4) :: icell,iv,jv,iqty,jqty,inbr,nnbr,node,iv1,jv1,ierror
      REAL(r8) :: timestart_fe,timeend_fe
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,inode)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: inode
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     zero matrix.  quadrature-point loop is now within integrands.
c-----------------------------------------------------------------------
      DO iv=0,tb%mvert
        mat(iv)%element=0._r8
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the integrand at all quadrature point.
c-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr
      CALL get_integrand(integrand,bigr,rdum,dx,dy,tb,0_i4)
c-----------------------------------------------------------------------
c     quadrature-point looping is now done inside the integrand routine,
c     and integrand is summed for each basis function, element by
c     element.
c
c     factors of Jacobian and quadrature weight are already in the basis
c     function arrays.
c-----------------------------------------------------------------------
      DO icell=1,tb%mcell
        DO iv=1,3
          iv1=tb%tgeom%vertex(icell,iv)
          DO jv=1,3
            jv1=tb%tgeom%vertex(icell,jv)
            DO inbr=0,SIZE(tb%tgeom%neighbor(iv1)%vertex)-1
              IF (tb%tgeom%neighbor(iv1)%vertex(inbr)==jv1) THEN
                mat(iv1)%element(1:nqty,1:nqty,inbr)
     $            =mat(iv1)%element(1:nqty,1:nqty,inbr)
     $            +integrand(:,:,icell,jv,iv)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_mat = time_mat + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_make_real_matrix
c-----------------------------------------------------------------------
c     subprogram 3. tblock_make_comp_matrix.
c     computes a linear response matrix for a supplied integrand
c     subprogram that fits the interface block at the beginning
c     of this subprogram.
c
c     complex version
c-----------------------------------------------------------------------
      SUBROUTINE tblock_make_comp_matrix(tb,mat,get_integrand,nqty)
      USE rblock_type_mod
      USE matrix_type_mod
      USE time

      TYPE(tblock_type), INTENT(INOUT), TARGET :: tb
      TYPE(comp_matrix_element_type3), DIMENSION(0:), INTENT(OUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      TYPE(rblock_type) :: rdum
      REAL(r8) :: dx=0,dy=0

      COMPLEX(r8), DIMENSION(nqty,nqty,tb%mcell,3,3) :: integrand
      REAL(r8), DIMENSION(:,:), POINTER :: bigr

      INTEGER(i4) :: icell,iv,jv,iqty,jqty,inbr,nnbr,node,iv1,jv1,ierror
      REAL(r8) :: timestart_fe,timeend_fe
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,inode)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: inode
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     zero matrix.  quadrature-point loop is now within integrands.
c-----------------------------------------------------------------------
      DO iv=0,tb%mvert
        mat(iv)%element=0._r8
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the integrand at all quadrature points.
c-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr
      CALL get_integrand(integrand,bigr,rdum,dx,dy,tb,0_i4)
c-----------------------------------------------------------------------
c     quadrature-point looping is now done inside the integrand routine,
c     and integrand is summed for each basis function, element by
c     element.
c
c     factors of Jacobian and quadrature weight are already in the basis
c     function arrays.
c-----------------------------------------------------------------------
      DO icell=1,tb%mcell
        DO iv=1,3
          iv1=tb%tgeom%vertex(icell,iv)
          DO jv=1,3
            jv1=tb%tgeom%vertex(icell,jv)
            DO inbr=0,SIZE(tb%tgeom%neighbor(iv1)%vertex)-1
              IF (tb%tgeom%neighbor(iv1)%vertex(inbr)==jv1) THEN
                mat(iv1)%element(1:nqty,1:nqty,inbr)
     $            =mat(iv1)%element(1:nqty,1:nqty,inbr)
     $            +integrand(:,:,icell,jv,iv)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_mat = time_mat + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_make_comp_matrix
c-----------------------------------------------------------------------
c     subprogram 4. tblock_get_real_rhs.
c     performs finite-element integrations for the rhs of an equation
c     where the integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_get_real_rhs(tb,rhs,get_integrand,nq)
      USE rblock_type_mod
      USE vector_type_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq
      TYPE(tblock_type), INTENT(INOUT), TARGET :: tb
      TYPE(vector_type), INTENT(INOUT) :: rhs

      TYPE(rblock_type) :: rdum
      REAL(r8) :: dx=0,dy=0

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: integrand
      REAL(r8), DIMENSION(:,:), POINTER :: bigr
      REAL(r8), DIMENSION(:,:,:), POINTER :: rhsg
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: rhss,rhsi

      INTEGER(i4) :: node,ivert,icell,ierror,iseg
      INTEGER(i4) :: iq,iv,nv,ib,mcell
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int
      REAL(r8) :: timestart_fe,timeend_fe
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,inode)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: inode
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=3
        rhsg=>rhs%arr
        rhsg(1:nq,:,:)=0
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_seg=iv
        n_seg=SIZE(rhs%arrh,4)
        rhss=>rhs%arrh
        rhss(1:nq,:,:,:)=0
      ELSE
        n_seg=0
      ENDIF
      iv=iv+3*n_seg
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,4)
        rhsi=>rhs%arri
        rhsi(1:nq,:,:,:)=0
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
      ALLOCATE(integrand(nq,tb%mcell,iv))
c-----------------------------------------------------------------------
c     evaluate the integrand at all quadrature points.
c-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr
      CALL get_integrand(integrand,bigr,rdum,dx,dy,tb,node)
c-----------------------------------------------------------------------
c     accumulate element contributions.
c     grid vertex-centered bases first.
c
c     factors of Jacobian and quadrature weight are already in the
c     test function arrays.
c-----------------------------------------------------------------------
      IF (n_grid==3) THEN
        DO ib=1,3
          DO icell=1,tb%mcell
            ivert=tb%tgeom%vertex(icell,ib)
            rhsg(1:nq,ivert,0)=rhsg(1:nq,ivert,0)+integrand(:,icell,ib)
          ENDDO
        ENDDO
      ENDIF
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
      END SUBROUTINE tblock_get_real_rhs
c-----------------------------------------------------------------------
c     subprogram 5. tblock_get_comp_rhs.
c     performs finite-element integrations for the rhs of an equation
c     where the integrand is computed with a supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_get_comp_rhs(tb,rhs,get_integrand,nq,nfour)
      USE rblock_type_mod
      USE vector_type_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,nfour
      TYPE(tblock_type), INTENT(INOUT), TARGET :: tb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      TYPE(rblock_type) :: rdum
      REAL(r8) :: dx=0,dy=0

      REAL(r8), DIMENSION(:,:), POINTER :: bigr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integrand
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhss,rhsi

      INTEGER(i4) :: node,ivert,icell,ierror,iseg
      INTEGER(i4) :: iq,iv,nv,ib,mcell,jf
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int
      REAL(r8) :: timestart_fe,timeend_fe
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,inode)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: inode
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=3
        rhsg=>rhs%arr
        rhsg(1:nq,:,:,1:nfour)=0
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_seg=iv
        n_seg=SIZE(rhs%arrh,4)
        rhss=>rhs%arrh
        rhss(1:nq,:,:,:,1:nfour)=0
      ELSE
        n_seg=0
      ENDIF
      iv=iv+3*n_seg
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,4)
        rhsi=>rhs%arri
        rhsi(1:nq,:,:,:,1:nfour)=0
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
      ALLOCATE(integrand(nq,tb%mcell,iv,nfour))
c-----------------------------------------------------------------------
c     evaluate the integrand at the quadrature point.
c-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr
      CALL get_integrand(integrand,bigr,rdum,dx,dy,tb,0_i4)
c-----------------------------------------------------------------------
c     accumulate element contributions.
c     grid vertex-centered bases first.
c
c     factors of Jacobian and quadrature weight are already in the
c     test function arrays.
c-----------------------------------------------------------------------
      IF (n_grid==3) THEN
        DO jf=1,nfour
          DO ib=1,3
            DO icell=1,tb%mcell
              ivert=tb%tgeom%vertex(icell,ib)
              rhsg(1:nq,ivert,0,jf)=rhsg(1:nq,ivert,0,jf)
     $          +integrand(:,icell,ib,jf)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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
      END SUBROUTINE tblock_get_comp_rhs
c-----------------------------------------------------------------------
c     subprogram 6. tblock_get_comp_rhs_q.
c     performs finite-element integrations for the rhs of an equation
c     where the integrand is computed with a supplied subroutine.
c     same as tblock_get_comp_rhs but the quantity and Fourier component
c     indices are assumed to have the correct dimension for efficiency.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_get_comp_rhs_q(tb,rhs,get_integrand)
      USE rblock_type_mod
      USE vector_type_mod
      USE time

      TYPE(tblock_type), INTENT(INOUT), TARGET :: tb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      TYPE(rblock_type) :: rdum
      REAL(r8) :: dx=0,dy=0

      REAL(r8), DIMENSION(:,:), POINTER :: bigr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integrand
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: rhsg
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: rhss,rhsi

      INTEGER(i4) :: node,ivert,icell,ierror,iseg
      INTEGER(i4) :: iq,iv,nv,ib,mcell,jf
      INTEGER(i4) :: start_seg,start_int,n_grid,n_seg,n_int,nq,nfour
      REAL(r8) :: timestart_fe,timeend_fe
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,inode)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: inode
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=3
        rhsg=>rhs%arr
        nq=SIZE(rhsg,1)
        nfour=SIZE(rhsg,4)
        rhsg=0
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_seg=iv
        n_seg=SIZE(rhs%arrh,4)
        rhss=>rhs%arrh
        nq=SIZE(rhss,1)
        nfour=SIZE(rhss,5)
        rhss=0
      ELSE
        n_seg=0
      ENDIF
      iv=iv+3*n_seg
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,4)
        rhsi=>rhs%arri
        nq=SIZE(rhsi,1)
        nfour=SIZE(rhsi,5)
        rhsi=0
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
      ALLOCATE(integrand(nq,tb%mcell,iv,nfour))
c-----------------------------------------------------------------------
c     evaluate the integrand at all quadrature points.
c-----------------------------------------------------------------------
      bigr=>tb%tgeom%bigr
      CALL get_integrand(integrand,bigr,rdum,dx,dy,tb,node)
c-----------------------------------------------------------------------
c     accumulate element contributions.
c     grid vertex-centered bases first.
c
c     factors of Jacobian and quadrature weight are already in the
c     test function arrays.
c-----------------------------------------------------------------------
      IF (n_grid==3) THEN
        DO jf=1,nfour
          DO ib=1,3
            DO icell=1,tb%mcell
              ivert=tb%tgeom%vertex(icell,ib)
              rhsg(:,ivert,0,jf)=rhsg(:,ivert,0,jf)
     $          +integrand(:,icell,ib,jf)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
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
      END SUBROUTINE tblock_get_comp_rhs_q
c-----------------------------------------------------------------------
c     subprogram 7. tblock_basis_set.
c     set the locations and weights for quadratures.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_basis_set(tb,geom)

      TYPE(tblock_type), INTENT(INOUT) :: tb
      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: node,iv,icell
c-----------------------------------------------------------------------
c     copy the basis function values into arrays with dimensions to
c     match those in structured blocks (except wjac, which is the
c     quadrature weight times the Jacobian).
c
c     for toroidal geometry dalpdrc holds d(alpha)/dr + alpha/r, and
c     for linear geometry it's just d(alpha)/dx.
c-----------------------------------------------------------------------
      ALLOCATE(tb%tgeom%alpha_arr(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpdr(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpdz(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpdrc(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%alpham_arr(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpmdr(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpmdz(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%dalpmdrc(tb%ng,tb%tgeom%mcell,3))
      ALLOCATE(tb%tgeom%bigr(tb%ng,tb%tgeom%mcell))
      ALLOCATE(tb%tgeom%wjac(tb%ng,tb%tgeom%mcell))
      DO node=1,tb%ng
        DO iv=1,3
          tb%tgeom%alpha_arr(node,:,iv)=tb%tgeom%alpha(iv,node)
        ENDDO
        tb%tgeom%dalpdr(node,:,:)=tb%tgeom%alpha_x(:,1,:)
        tb%tgeom%dalpdz(node,:,:)=tb%tgeom%alpha_y(:,1,:)
        tb%tgeom%dalpdrc(node,:,:)=tb%tgeom%alpha_x(:,1,:)
        IF (geom=='tor') THEN
          DO icell=1,tb%tgeom%mcell
            tb%tgeom%bigr(node,icell)
     $        =tb%tgeom%xs(tb%tgeom%vertex(icell,1))*
     $         tb%tgeom%alpha(1,node)
     $        +tb%tgeom%xs(tb%tgeom%vertex(icell,2))*
     $         tb%tgeom%alpha(2,node)
     $        +tb%tgeom%xs(tb%tgeom%vertex(icell,3))*
     $         tb%tgeom%alpha(3,node)
          ENDDO
          DO iv=1,3
            tb%tgeom%dalpdrc(node,:,iv)=
     $        tb%tgeom%dalpdrc(node,:,iv)+
     $        tb%tgeom%alpha_arr(node,:,iv)/tb%tgeom%bigr(node,:)
          ENDDO
        ENDIF
        tb%tgeom%wjac(node,:)=tb%tgeom%area(:)*tb%wg(node)
      ENDDO
      IF (geom=='tor') THEN
        tb%tgeom%wjac=tb%tgeom%wjac*tb%tgeom%bigr
      ELSE
        tb%tgeom%bigr=1
      ENDIF
c-----------------------------------------------------------------------
c     generate copies of the basis/test functions for matrices, where
c     the square root of Jacobian times the quadrature weight is
c     already multiplied.
c
c     test functions used for rhs computations are multiplied by a full
c     factor of Jacobian time quadrature weight.
c-----------------------------------------------------------------------
      DO iv=1,3
        tb%tgeom%alpham_arr(:,:,iv)=tb%tgeom%alpha_arr(:,:,iv)*
     $                              SQRT(tb%tgeom%wjac)
        tb%tgeom%dalpmdr(:,:,iv)=tb%tgeom%dalpdr(:,:,iv)*
     $                           SQRT(tb%tgeom%wjac)
        tb%tgeom%dalpmdz(:,:,iv)=tb%tgeom%dalpdz(:,:,iv)*
     $                           SQRT(tb%tgeom%wjac)
        tb%tgeom%dalpmdrc(:,:,iv)=tb%tgeom%dalpdrc(:,:,iv)*
     $                            SQRT(tb%tgeom%wjac)
        tb%tgeom%alpha_arr(:,:,iv)=tb%tgeom%alpha_arr(:,:,iv)*
     $                             tb%tgeom%wjac
        tb%tgeom%dalpdr(:,:,iv)=tb%tgeom%dalpdr(:,:,iv)*tb%tgeom%wjac
        tb%tgeom%dalpdz(:,:,iv)=tb%tgeom%dalpdz(:,:,iv)*tb%tgeom%wjac
        tb%tgeom%dalpdrc(:,:,iv)=tb%tgeom%dalpdrc(:,:,iv)*tb%tgeom%wjac
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_basis_set
c-----------------------------------------------------------------------
c     subprogram 8. tblock_basis_dealloc.
c     deallocates space for tblock bases.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_basis_dealloc(tg)

      TYPE(tri_linear_geom_type), INTENT(INOUT) :: tg

      DEALLOCATE(tg%alpha_arr)
      DEALLOCATE(tg%dalpdr)
      DEALLOCATE(tg%dalpdz)
      DEALLOCATE(tg%dalpdrc)
      DEALLOCATE(tg%alpham_arr)
      DEALLOCATE(tg%dalpmdr)
      DEALLOCATE(tg%dalpmdz)
      DEALLOCATE(tg%dalpmdrc)
      DEALLOCATE(tg%bigr)
      DEALLOCATE(tg%wjac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_basis_dealloc
c-----------------------------------------------------------------------
c     subprogram 9. tblock_real_qp_update.
c     evaluate 2D tri_linear data at the gaussian quadrature points
c     for this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_update(tl,qtl,tb)

      TYPE(tri_linear_2D_type), INTENT(IN) :: tl
      TYPE(tb_real_qp_type), INTENT(INOUT) :: qtl
      TYPE(tblock_type), INTENT(IN) :: tb

      REAL(r8), DIMENSION(tl%nqty,tb%mcell,1) :: f,fr,fz
      
      INTEGER(i4) :: ig
c-----------------------------------------------------------------------
c     loop over quadrature points, evaluate, and save.
c-----------------------------------------------------------------------
      DO ig=1,tb%ng
        CALL tri_linear_all_eval(tl,tb%tgeom,ig,1_i4,f,fr,fz)
        qtl%qpf(:,ig,:)=f(:,:,1)
        qtl%qpfr(:,ig,:)=fr(:,:,1)
        qtl%qpfz(:,ig,:)=fz(:,:,1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_real_qp_update
c-----------------------------------------------------------------------
c     subprogram 10. tblock_comp_qp_update.
c     evaluate 3D tri_linear data at the gaussian quadrature points
c     for this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_update(tl,qtl,tb)

      TYPE(tri_linear_type), INTENT(IN) :: tl
      TYPE(tb_comp_qp_type), INTENT(INOUT) :: qtl
      TYPE(tblock_type), INTENT(IN) :: tb
      
      COMPLEX(r8), DIMENSION(tl%nqty,tb%mcell,1,tl%nfour) :: f,fr,fz

      INTEGER(i4) :: ig
c-----------------------------------------------------------------------
c     loop over quadrature points, evaluate, and save.
c-----------------------------------------------------------------------
      DO ig=1,tb%ng
        CALL tri_linear_all_eval(tl,tb%tgeom,ig,1_i4,f,fr,fz)
        qtl%qpf(:,ig,:,:)=f(:,:,1,:)
        qtl%qpfr(:,ig,:,:)=fr(:,:,1,:)
        qtl%qpfz(:,ig,:,:)=fz(:,:,1,:)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_comp_qp_update
c-----------------------------------------------------------------------
c     subprogram 11. tblock_real_qp_alloc.
c     allocate space for 2D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_alloc(qtl,tb,nqty)

      TYPE(tb_real_qp_type), INTENT(OUT) :: qtl
      TYPE(tblock_type), INTENT(IN) :: tb
      INTEGER(i4), INTENT(IN) :: nqty
      
      ALLOCATE(qtl%qpf(nqty,tb%ng,tb%mcell))
      ALLOCATE(qtl%qpfr(nqty,tb%ng,tb%mcell))
      ALLOCATE(qtl%qpfz(nqty,tb%ng,tb%mcell))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_real_qp_alloc
c-----------------------------------------------------------------------
c     subprogram 12. tblock_comp_qp_alloc.
c     allocate space for 3D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_alloc(qtl,tb,nqty,nfour)

      TYPE(tb_comp_qp_type), INTENT(OUT) :: qtl
      TYPE(tblock_type), INTENT(IN) :: tb
      INTEGER(i4), INTENT(IN) :: nqty,nfour
      
      ALLOCATE(qtl%qpf(nqty,tb%ng,tb%mcell,nfour))
      ALLOCATE(qtl%qpfr(nqty,tb%ng,tb%mcell,nfour))
      ALLOCATE(qtl%qpfz(nqty,tb%ng,tb%mcell,nfour))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_comp_qp_alloc
c-----------------------------------------------------------------------
c     subprogram 13. tblock_real_qp_dealloc.
c     deallocate space for 2D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_real_qp_dealloc(qtl)

      TYPE(tb_real_qp_type), INTENT(INOUT) :: qtl

      DEALLOCATE(qtl%qpf)
      IF (ALLOCATED(qtl%qpfr)) DEALLOCATE(qtl%qpfr)
      IF (ALLOCATED(qtl%qpfz)) DEALLOCATE(qtl%qpfz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_real_qp_dealloc
c-----------------------------------------------------------------------
c     subprogram 14. tblock_comp_qp_dealloc.
c     deallocate space for 3D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE tblock_comp_qp_dealloc(qtl)

      TYPE(tb_comp_qp_type), INTENT(INOUT) :: qtl
      
      DEALLOCATE(qtl%qpf)
      IF (ALLOCATED(qtl%qpfr)) DEALLOCATE(qtl%qpfr)
      IF (ALLOCATED(qtl%qpfz)) DEALLOCATE(qtl%qpfz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tblock_comp_qp_dealloc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE tblock
