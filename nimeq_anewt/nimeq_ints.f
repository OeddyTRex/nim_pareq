c-----------------------------------------------------------------------
c     file nimeq_ints.f
c     module that includes integrand routines used in the diagnostic
c     finite element integrations in nimplot.
c
c     note that nimplot does not use quadrature point storage, so the
c     interpolations are made on the fly with generic_all_eval instead
c     of setting pointers like in nimrod's integrands.f.  this also
c     means that the phi component of the interpolated equilibrium B is
c     R*B_eq_phi, while the phi component of equil J is J_eq_phi/R.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  get_mass.
c     2.  mag_op.
c     3.  delstar_op.
c     4.  delstlj_op.
c     5.  delstar_rhs.
c     6.  gs_rhs.
c     7.  err_gs.
c     8.  get_B.
c     9.  jphi_rhs.
c     10. delstarlj_dot.
c     11. delstar_dot.
c     12. err_ip.
c     13. scal_adv_op.
c     14. scal_adv_rhs.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE nimeq_ints

      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE lagr_quad_mod
      USE tri_linear
      USE generic_evals
      USE math_tran
      USE physdat
      USE input
      USE input_eq
      USE global
      USE fft_mod
      IMPLICIT NONE

      CONTAINS

c-----------------------------------------------------------------------
c     subprogram 1. get_mass.
c     compute the integrand used for the mass matrix.
c-----------------------------------------------------------------------
      SUBROUTINE get_mass(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,5)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
c-----------------------------------------------------------------------
c     find alpha(jv)*alpha(iv).
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM(alpha(:,:,jv)*alpha(:,:,iv),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_mass
c-----------------------------------------------------------------------
c     subprogram 2. mag_op.
c     this routine provides a mass matrix for a 3-vector, but the
c     first two components are multiplied by R.
c-----------------------------------------------------------------------
      SUBROUTINE mag_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,5)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
c-----------------------------------------------------------------------
c     find alpha(jv)*alpha(iv), multiplied by R for the first two
c     vector components.
c-----------------------------------------------------------------------
      int=0._r8
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM(alpha(:,:,jv)*alpha(:,:,iv)*bigr,1)
          int(2,2,:,jv,iv)=SUM(alpha(:,:,jv)*alpha(:,:,iv)*bigr,1)
          int(3,3,:,jv,iv)=SUM(alpha(:,:,jv)*alpha(:,:,iv),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mag_op
c-----------------------------------------------------------------------
c     subprogram 3. delstar_op.
c     the delstar operator acting on a field psi can be written as
c
c       div[ R**2 grad( psi/R**2 ) ]
c
c     the operator computed here is used to find psi/R**2 (not psi)
c     and in weak form after integration by parts, the integrand is
c
c       R**2*grad(test_function).grad(basis_function)
c
c     and the minus sign is moved to the right side of the equation.
c-----------------------------------------------------------------------
      SUBROUTINE delstar_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,5)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     find R^2* grad(alpha(jv)).grad(alpha(iv)).
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM((dalpdr(:,:,jv)*dalpdr(:,:,iv)+
     $           dalpdz(:,:,jv)*dalpdz(:,:,iv))*bigr(:,:)**2,1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstar_op
c-----------------------------------------------------------------------
c     subprogram 4. delstlj_op.
c     this version of the operator is used to solve a 2-vector system,
c     where the first component is lambda=psi/R**2 and the second is
c     mu0*J_phi/R.
c-----------------------------------------------------------------------
      SUBROUTINE delstlj_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,5)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     find R^2* grad(alpha(jv)).grad(alpha(iv)) for the first diagonal
c     entry and R^2*alpha(jv)*alpha(iv) for the second.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM((dalpdr(:,:,jv)*dalpdr(:,:,iv)+
     $           dalpdz(:,:,jv)*dalpdz(:,:,iv))*bigr(:,:)**2,1)
          int(2,1,:,jv,iv)=0._r8
          int(1,2,:,jv,iv)=0._r8
          int(2,2,:,jv,iv)=SUM(alpha(:,:,jv)*alpha(:,:,iv)*
     $                         bigr(:,:)**2,1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstlj_op
c-----------------------------------------------------------------------
c     subprogram 5. delstar_rhs.
c     compute the rhs integrand used for the del_star solve.
c     this is -alpha*mu0*Jphi*R, and the minus sign appears due to an
c     integration-by-parts on the lhs.
c-----------------------------------------------------------------------
      SUBROUTINE delstar_rhs(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ja,
     $     bez,ber,be
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: ja_eq
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: psi,
     $     psiz,psir
      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_all_eval(rb%ja_eq,tb%ja_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ja_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8)
      CALL generic_all_eval(rb%rwork1,tb%rwork1,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      psi,psir,psiz,1_i4)
c-----------------------------------------------------------------------
c     convert the phi components of equilibrium fields to cylindrical.
c-----------------------------------------------------------------------
      ja_eq(3,:,:)=mu0*ja_eq(3,:,:)*bigr
      DO im=1,nmodes
        IF (keff(im)==0) THEN
          ja(:,:,:,im)=REAL(ja(:,:,:,im))
          IF (nonlinear) ja_eq=ja_eq+ja(:,:,:,im)
          EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     find alpha*mu*Jphi*R + grad(alpha)*grad(psi_s/R^2)*R^2, where
c     psi_s is from the surface only.
c-----------------------------------------------------------------------
      DO iv=1,nv
        int(1,:,iv)=SUM(-(alpha(:,:,iv)*ja_eq(3,:,:)*bigr
     $                 +(dalpdr(:,:,iv)*psir(1,:,:)+
     $                   dalpdz(:,:,iv)*psiz(1,:,:))*bigr(:,:)**2),1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstar_rhs
c-----------------------------------------------------------------------
c     subprogram 6. gs_rhs.
c     compute the rhs integrand used for the gs solve,
c       -alpha ( F*F' + mu0*R**2*P' )
c-----------------------------------------------------------------------
      SUBROUTINE gs_rhs(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: f,press,dp
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: psib,
     $     psizb,psirb,psi
      REAL(r8), DIMENSION(1,1,1) :: dv
      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      IF (gs_type/='free')
     $  CALL generic_all_eval(rb%rwork1,tb%rwork1,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,psib,psirb,psizb,1_i4)
      CALL generic_ptr_set(rb%qrwork1,tb%qrwork1,tb%tgeom,
     $                     inode,f,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                     inode,press,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find alpha*(F*dF/dpsi+mu0*R**2*dP/dpsi) -
c     grad(alpha)*grad(psi_s/R^2)*R^2, where the second term is from
c     the surface flux.  if this is a free-boundary computation, there
c     is no separate inhomogeneous surface term.
c
c     note that press holds mu0*dP/dpsi at this point, not pressure, and
c     f has F*dF/dpsi.
c
c     the second quantity is without the surface-flux term for obtaining
c     a projection of mu0*J_phi/R, directly, during free-boundary
c     computations.
c-----------------------------------------------------------------------
      IF (SIZE(int,1)==1) THEN
        DO iv=1,nv
          int(1,:,iv)=SUM(
     $          alpha(:,:,iv)*(f(1,:,:)+bigr(:,:)**2*press(1,:,:))
     $       -(dalpdr(:,:,iv)*psirb(1,:,:)+
     $         dalpdz(:,:,iv)*psizb(1,:,:))*bigr(:,:)**2,1)
        ENDDO
      ELSE
        DO iv=1,nv
          int(1,:,iv)=SUM(
     $          alpha(:,:,iv)*(f(1,:,:)+bigr(:,:)**2*press(1,:,:)),1)
          int(2,:,iv)=-int(1,:,iv)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE gs_rhs
c-----------------------------------------------------------------------
c     subprogram 7. err_gs.
c     compute the integrand for rhs-delstar(psi).  this routine is
c     used for error estimates
c-----------------------------------------------------------------------
      SUBROUTINE err_gs(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz    
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: f,press,dp
      INTEGER(i4) :: iv,nv,im,i,j
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: psi,
     $     psiz,psir,psib,psizb,psirb
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: jaeq
      REAL(r8), DIMENSION(1,1,1) :: dv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_all_eval(rb%rwork2,tb%rwork2,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      psi,psir,psiz,1_i4)
      CALL generic_ptr_set(rb%qrwork1,tb%qrwork1,tb%tgeom,
     $                     inode,f,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                     inode,press,dp,dp,0_i4)
      IF (gs_type=='free') THEN
        CALL generic_all_eval(rb%ja_eq,tb%ja_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        jaeq,dv,dv,0_i4)
      ELSE
        CALL generic_all_eval(rb%rwork1,tb%rwork1,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,psib,psirb,psizb,1_i4)
c-----------------------------------------------------------------------
c       combine interior and surface contributions to psi/R**2
c-----------------------------------------------------------------------
        psi  = psi  + psib
        psir = psir + psirb
        psiz = psiz + psizb
      ENDIF
c-----------------------------------------------------------------------
c     find grad(alpha)*grad(psi/R^2)*R^2
c       -alpha*(F*dF/dpsi+mu0*R**2*dP/dpsi)
c
c     note that press holds mu0*dP/dpsi at this point, not pressure, and
c     f has F*dF/dpsi.
c
c     when there is a second vector component, mu0*J_phi/R, from the
c     linear relation with delstar(psi) for use in free-boundary solves.
c-----------------------------------------------------------------------
      IF (integrand_flag=="all terms") THEN
        IF (SIZE(int,1)==1) THEN
          DO iv=1,nv
            int(1,:,iv) = SUM(
     $        bigr(:,:)**2*(dalpdr(:,:,iv)*psir(1,:,:)+
     $                      dalpdz(:,:,iv)*psiz(1,:,:))-
     $        (f(1,:,:)+bigr(:,:)**2*press(1,:,:))*alpha(:,:,iv),1)
          ENDDO
        ELSE
          DO iv=1,nv
            int(1,:,iv) = SUM(
     $        bigr(:,:)**2*(dalpdr(:,:,iv)*psir(1,:,:)+
     $                      dalpdz(:,:,iv)*psiz(1,:,:))-
     $        (f(1,:,:)+bigr(:,:)**2*press(1,:,:))*alpha(:,:,iv),1)
            int(2,:,iv) = SUM(
     $        alpha(:,:,iv)*(bigr(:,:)**2*(jaeq(3,:,:)+press(1,:,:))+
     $                       f(1,:,:)),1)
          ENDDO
        ENDIF
      ELSE          !    just delstar part
        IF (SIZE(int,1)==1) THEN
          DO iv=1,nv
            int(1,:,iv) = SUM(
     $        bigr(:,:)**2*(dalpdr(:,:,iv)*psir(1,:,:)+
     $                      dalpdz(:,:,iv)*psiz(1,:,:)),1)
          ENDDO
        ELSE
          DO iv=1,nv
            int(1,:,iv) = SUM(
     $        bigr(:,:)**2*(dalpdr(:,:,iv)*psir(1,:,:)+
     $                      dalpdz(:,:,iv)*psiz(1,:,:)),1)
            int(2,:,iv)=0._r8
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE err_gs
c-----------------------------------------------------------------------
c     subprogram 8. get_B
c     calculate the RHS used to determine B from psi and F.
c       RBr = -dpsi/dz
c       RBz = dpsi/dr
c       RBp = F
c-----------------------------------------------------------------------
      SUBROUTINE get_B(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: psi,psiz,psir,
     $          rbp
      REAL(r8), DIMENSION(1,1,1) :: dv
      INTEGER(i4) :: iv,im,nv,ix,iy,nx,ny
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate the
c     incoming vector.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%rwork2,tb%rwork2,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      psi,psir,psiz,1_i4)
      CALL generic_all_eval(rb%rwork3,tb%rwork3,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      rbp,dv,dv,0_i4)
c-----------------------------------------------------------------------
c     B-poloidal is grad(phi)Xgrad(psi), but use the expansion for
c     psi/R**2 that includes surface contributions at this point
c     (local array psi):
c
c     R*Bpol=phi-hat X (R**2*grad(psi/R**2)+2*R*(psi/R**2)*R_hat)
c
c     the third component is R*B_phi = F.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)

      IF (geom=='tor') THEN
        DO iv=1,nv
          int(1,:,iv)=SUM(-alpha(:,:,iv)* psiz(1,:,:)*bigr**2,1)
          int(2,:,iv)=SUM( alpha(:,:,iv)*(psir(1,:,:)*bigr**2+
     $                         2._r8*bigr*psi (1,:,:)),1)
          int(3,:,iv)=SUM( alpha(:,:,iv)*rbp(1,:,:),1)
        ENDDO
      ELSE
        DO iv=1,nv
          int(1,:,iv)=SUM(-alpha(:,:,iv)*psiz(1,:,:),1)
          int(2,:,iv)=SUM( alpha(:,:,iv)*psir(1,:,:),1)
          int(3,:,iv)=SUM( alpha(:,:,iv)*rbp(1,:,:),1)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_B
c-----------------------------------------------------------------------
c     subprogram 9. jphi_rhs.
c     compute mu0*Jphi/R from the rhs of the GS equation.
c-----------------------------------------------------------------------
      SUBROUTINE jphi_rhs(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz    
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: psi,rlen
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: bigr2,psival,rlenval
      REAL(r8), EXTERNAL :: f_func,p_func
      INTEGER(i4) :: iv,nv,ix,iy,ncx,ncy
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.  rwork2 now has the
c     full psi/R**2 distribution across the interior and at the wall.
c
c     fllen holds the relative length of the field, as used to check
c     topology.
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_all_eval(rb%rwork2,tb%rwork2,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      psi,dv,dv,0_i4)
      CALL generic_all_eval(rb%fllen,tb%fllen,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      rlen,dv,dv,0_i4)
c-----------------------------------------------------------------------
c     find -alpha*(F*dF/dpsi/R**2+mu0*dP/dpsi).  here, we use the
c     external functions f_func and p_func which provide the user-
c     specified distributions according to psi.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO iy=1,ncy
          int(1,iy,iv)=0._r8
          DO ix=1,ncx
            bigr2=bigr(ix,iy)**2
            psival=psi(1,ix,iy)*bigr2
            rlenval=rlen(1,ix,iy)
            int(1,iy,iv)=int(1,iy,iv)-alpha(ix,iy,iv)*(
     $        f_func(psival,0_i4,rlenval,.false.)*
     $        f_func(psival,1_i4,rlenval,.false.)/bigr2+
     $        p_func(psival,1_i4,rlenval,.false.) )
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jphi_rhs
c-----------------------------------------------------------------------
c     subprogram 10. delstarlj_dot.
c     compute the effect of the delstar operator applied to one 
c     component of a direction vector from a matrix-free computation,
c     and apply a mass matrix to the second.
c-----------------------------------------------------------------------
      SUBROUTINE delstarlj_dot(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: delffp,delpp,dp
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: laj,
     $          lajr,lajz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2
      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     the direction vector from the solve is in the 2nd and 3rd 
c     components of be_eq.  the 2nd component is lambda=psi/R**2,
c     and the third is mu0*J_phi/R.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      laj,lajr,lajz,1_i4)
c-----------------------------------------------------------------------
c     if the linearized F*F' and P' terms are used, point to their
c     partials with respect to lambda, and form
c       grad(alpha)*grad(lambda)*R^2-
c         alpha*[(del(FF')/dellam)+mu0*R^2*(del(P')/dellam)]*lambda
c     and
c       alpha*mu0*J_phi*R+
c         alpha*[(del(FF')/dellam)+mu0*R^2*(del(P')/dellam)]*lambda
c-----------------------------------------------------------------------
      bigr2=bigr(:,:)**2
      IF (dellam>0._r8) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,delffp,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                       inode,delpp,dp,dp,0_i4)
        DO iv=1,nv
          int(1,:,iv)=SUM((dalpdr(:,:,iv)*lajr(2,:,:)+
     $                     dalpdz(:,:,iv)*lajz(2,:,:))*bigr2(:,:)-
     $                (delffp(1,:,:)+bigr2(:,:)*delpp(1,:,:))/dellam*
     $                alpha(:,:,iv)*laj(2,:,:),1)
          int(2,:,iv)=SUM(alpha(:,:,iv)*bigr2(:,:)*laj(3,:,:)+
     $                (delffp(1,:,:)+bigr2(:,:)*delpp(1,:,:))/dellam*
     $                alpha(:,:,iv)*laj(2,:,:),1)
        ENDDO
      ELSE
c-----------------------------------------------------------------------
c       just form grad(alpha)*grad(lambda)*R^2 and alpha*mu0*J_phi*R:
c-----------------------------------------------------------------------
        DO iv=1,nv
          int(1,:,iv)=SUM((dalpdr(:,:,iv)*lajr(2,:,:)+
     $                     dalpdz(:,:,iv)*lajz(2,:,:))*bigr2(:,:),1)
          int(2,:,iv)=SUM(alpha(:,:,iv)*bigr2(:,:)*laj(3,:,:),1)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstarlj_dot
c-----------------------------------------------------------------------
c     subprogram 11. delstar_dot.
c     compute the effect of the delstar operator applied to one 
c     component of a direction vector from a matrix-free computation.
c-----------------------------------------------------------------------
      SUBROUTINE delstar_dot(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: delffp,delpp,dp
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: lam,
     $          lamr,lamz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2
      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      bigr2=bigr(:,:)**2
c-----------------------------------------------------------------------
c     the direction vector from the solve is in rwork3.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%rwork3,tb%rwork3,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      lam,lamr,lamz,1_i4)
c-----------------------------------------------------------------------
c     if the linearized F*F' and P' terms are used, point to their
c     partials with respect to lambda, and form
c      grad(alpha)*grad(lambda)*R^2-
c      [(del(FF')/dellam)+mu0*R^2*(del(P')/dellam)]*lambda
c-----------------------------------------------------------------------
      IF (dellam>0._r8) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,delffp,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                       inode,delpp,dp,dp,0_i4)
        DO iv=1,nv
          int(1,:,iv)=SUM((dalpdr(:,:,iv)*lamr(1,:,:)+
     $                     dalpdz(:,:,iv)*lamz(1,:,:))*bigr2(:,:)-
     $                (delffp(1,:,:)+bigr(:,:)**2*delpp(1,:,:))/dellam*
     $                alpha(:,:,iv)*lam(1,:,:),1)
        ENDDO
      ELSE
c-----------------------------------------------------------------------
c       just grad(alpha)*grad(lambda)*R^2
c-----------------------------------------------------------------------
        DO iv=1,nv
          int(1,:,iv)=SUM((dalpdr(:,:,iv)*lamr(1,:,:)+
     $                     dalpdz(:,:,iv)*lamz(1,:,:))*bigr2(:,:),1)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE delstar_dot
c-----------------------------------------------------------------------
c     subprogram 12. err_ip.
c     when performing newton iteration with coefficients varied to
c     attain a target plasma current, the column of the Jacobian
c     matrix associated with the F-parameter (which is varied) and the
c     row associated with the error in current are eliminated via
c     matrix partitioning.  that row and column are created as vector
c     structures through this finite-element computation.
c-----------------------------------------------------------------------
      SUBROUTINE err_ip(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: dffdf,delffp,dp
      REAL(r8), DIMENSION(1,1,1) :: dv
      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_ptr_set(rb%qrwork4,tb%qrwork4,tb%tgeom,
     $                     inode,dffdf,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                     inode,delffp,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find -alpha*d(FF)/d(f1), where f1 is the parameter being used to
c     adjust current.  also, find -ip_weight*alpha*(d(FF)/d(lam))/R**2,
c     which is the influence of lambda on the plasma current.
c-----------------------------------------------------------------------
      DO iv=1,nv
        int(1,:,iv)=-SUM(alpha(:,:,iv)*dffdf(1,:,:),1)
        int(2,:,iv)=-ip_weight/dellam*
     $               SUM(alpha(:,:,iv)*delffp(1,:,:)/bigr(:,:)**2,1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE err_ip
c-----------------------------------------------------------------------
c     subprogram 13. scal_adv_op.
c     this routine produces the integrand for a passive scalar advection
c     along the poloidal magnetic-field direction vector, according to
c     the current distribution of lambda.
c-----------------------------------------------------------------------
      SUBROUTINE scal_adv_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: lam,lamr,lamz,
     $          blam,blamr,blamz
      REAL(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2)) :: bdir
      REAL(r8) :: bsign

      INTEGER(i4) :: nv,iv,jv
c-----------------------------------------------------------------------
c     convenience parameters.
c-----------------------------------------------------------------------
      nv=SIZE(int,5)
      IF (integrand_flag=="reverse") THEN
        bsign=-1._r8
      ELSE
        bsign= 1._r8
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     find the direction vector for the current value of lambda.
c     fixed solves also need the fixed border values that are stored
c     separately
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%rwork2,tb%rwork2,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      lam,lamr,lamz,1_i4)
      IF (gs_type=="fixed") THEN
        CALL generic_all_eval(rb%rwork1,tb%rwork1,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,blam,blamr,blamz,1_i4)
        lam=lam+blam
        lamr=lamr+blamr
        lamz=lamz+blamz
      ENDIF
      IF (geom=='tor') THEN
        bdir(1,:,:)=-lamz(1,:,:)*bigr**2
        bdir(2,:,:)= lamr(1,:,:)*bigr**2+2._r8*bigr*lam(1,:,:)
      ELSE
        bdir(1,:,:)=-lamz(1,:,:)
        bdir(2,:,:)= lamr(1,:,:)
      ENDIF
      bdir(1,:,:)=bsign*blen_check*bdir(1,:,:)/flux_n0
      bdir(2,:,:)=bsign*blen_check*bdir(2,:,:)/flux_n0
c-----------------------------------------------------------------------
c     find {alpha(i)+L*b.grad(alpha(i))}*{alpha(j)+L*b.grad(alpha(j))},
c     where L is a length scale over which to advect the scalar.
c     [At present L=blen_check.]
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
           int(1,1,:,jv,iv)=SUM(
     $       (alpha(:,:,iv)+bdir(1,:,:)*dalpdr(:,:,iv)+
     $                      bdir(2,:,:)*dalpdz(:,:,iv))*
     $       (alpha(:,:,jv)+bdir(1,:,:)*dalpdr(:,:,jv)+
     $                      bdir(2,:,:)*dalpdz(:,:,jv)),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE scal_adv_op
c-----------------------------------------------------------------------
c     subprogram 14. scal_adv_rhs.
c     this routine produces the rhs integrand for passive scalar
c     advection along the poloidal magnetic-field direction vector, 
c     according to the current distribution of lambda.
c-----------------------------------------------------------------------
      SUBROUTINE scal_adv_rhs(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: lam,lamr,lamz,
     $          blam,blamr,blamz
      REAL(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2)) :: bdir
      REAL(r8) :: bsign

      INTEGER(i4) :: nv,iv
c-----------------------------------------------------------------------
c     convenience parameters.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      IF (integrand_flag=="reverse") THEN
        bsign=-1._r8
      ELSE
        bsign= 1._r8
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     find the direction vector for the current value of lambda.
c     fixed solves also need the fixed border values that are stored
c     separately
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%rwork2,tb%rwork2,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      lam,lamr,lamz,1_i4)
      IF (gs_type=="fixed") THEN
        CALL generic_all_eval(rb%rwork1,tb%rwork1,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,blam,blamr,blamz,1_i4)
        lam=lam+blam
        lamr=lamr+blamr
        lamz=lamz+blamz
      ENDIF
      IF (geom=='tor') THEN
        bdir(1,:,:)=-lamz(1,:,:)*bigr**2
        bdir(2,:,:)= lamr(1,:,:)*bigr**2+2._r8*bigr*lam(1,:,:)
      ELSE
        bdir(1,:,:)=-lamz(1,:,:)
        bdir(2,:,:)= lamr(1,:,:)
      ENDIF
      bdir(1,:,:)=bsign*blen_check*bdir(1,:,:)/flux_n0
      bdir(2,:,:)=bsign*blen_check*bdir(2,:,:)/flux_n0
c-----------------------------------------------------------------------
c     find {alpha(i)+L*b.grad(alpha(i))}*old_scalar to make the rhs for
c     the backward Euler approximation using least-squares projection.
c     here, we set old_scalar=1 throughout the interior.
c-----------------------------------------------------------------------
      DO iv=1,nv
         int(1,:,iv)=SUM(
     $     (alpha(:,:,iv)+bdir(1,:,:)*dalpdr(:,:,iv)+
     $                    bdir(2,:,:)*dalpdz(:,:,iv)),1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE scal_adv_rhs

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimeq_ints
