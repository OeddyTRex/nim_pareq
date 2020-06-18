c-----------------------------------------------------------------------
c     file nimplot_ints.f
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
c     2.  curl.
c     3.  energy_density.
c     4.  divb_int.
c     5.  jpar_int.
c     6.  edotj_int.
c     7.  heat_flux_int.
c     8.  mach_int.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE nimplot_ints

      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE lagr_quad_mod
      USE tri_linear
      USE generic_evals
      USE math_tran
      USE physdat
      USE input
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
      nv=SIZE(int,4)
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
c     subprogram 2. curl
c     form the integrand for a curl operation.
c     the incoming vector is stored in the work1 bilinear_type structure
c     in Cartesian components or or cylindrical components.
c-----------------------------------------------------------------------
      SUBROUTINE curl(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,1,1) :: dv
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: vec,
     $             dvecr,dvecz
      REAL(r8) :: tor_fac
      INTEGER(i4) :: iv,im,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate the
c     incoming vector.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%work1,tb%work1,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      vec,dvecr,dvecz,1_i4)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set factor to collect extra terms according to geom.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        tor_fac=1
      ELSE
        tor_fac=0
      ENDIF
c-----------------------------------------------------------------------
c     construct alpha.curl(vector) in cylindrical components
c     (poloidal plane).
c     note:  vector storage is (x,y,z) or (r,z,phi) for each complex
c     Fourier component.
c-----------------------------------------------------------------------
      DO im=1,nmodes
        DO iv=1,nv
          int(1,:,iv,im)=SUM( alpha(:,:,iv)*( dvecz(3,:,:,im)
     $                          -(0,1)*keff(im)*vec(2,:,:,im) /bigr ),1)
          int(2,:,iv,im)=SUM(-alpha(:,:,iv)*( dvecr(3,:,:,im)
     $      +(tor_fac*vec(3,:,:,im)-
     $        (0,1)*keff(im)*vec(1,:,:,im))/bigr ),1)
          int(3,:,iv,im)=SUM( alpha(:,:,iv)*( dvecr(2,:,:,im)
     $                                       -dvecz(1,:,:,im) ),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curl
c-----------------------------------------------------------------------
c     subprogram 3. energy_density.
c     compute b**2/2*mu and ro*v**2/2 mode by mode.  they are placed in
c     the first two storage spaces in rhs for each mode.
c-----------------------------------------------------------------------
      SUBROUTINE energy_density(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be_eq,ve_eq
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: nd_eq,pres_eq,
     $          prese_eq,nd_n0,ti_eq,te_eq
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: nd_sym,presi_eq,
     $          wj
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: real_ve
      REAL(r8), DIMENSION(1,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: real_ke,
     $          real_nd
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kfac,pfac

      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ve,ve_tot
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             kin_e,tele,tion,nd
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: im,ncx,ncy,im0,mps
c-----------------------------------------------------------------------
c     convenience parameters
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      mps=ncx*ncy
c-----------------------------------------------------------------------
c     evaluate the perturbed and equilibrium fields.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,dc,dc,0_i4)
      CALL generic_all_eval(rb%ve,tb%ve,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,dc,dc,0_i4)
      IF (continuity=='none') THEN
        nd=0._r8
      ELSE
        CALL generic_all_eval(rb%nd,tb%nd,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        nd,dc,dc,0_i4)
      ENDIF
      ve_tot=ve
      IF (beta>0) THEN
        CALL generic_all_eval(rb%tele,tb%tele,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        tele,dc,dc,0_i4)
        CALL generic_all_eval(rb%tion,tb%tion,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        tion,dc,dc,0_i4)
        CALL generic_all_eval(rb%tele_eq,tb%tele_eq,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,te_eq,dv,dv,0_i4)
        CALL generic_all_eval(rb%tion_eq,tb%tion_eq,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,ti_eq,dv,dv,0_i4)
      ENDIF
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%pres_eq,tb%pres_eq,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                      tb%ng,pres_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%prese_eq,tb%prese_eq,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                      tb%ng,prese_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%nd_eq,tb%nd_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      nd_eq,dv,dv,0_i4)
      IF (eq_flow/='none') THEN
        CALL generic_all_eval(rb%ve_eq,tb%ve_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        ve_eq,dv,dv,0_i4)
        IF (nonlinear) THEN
          DO im=1,nmodes
            IF (keff(im)==0) THEN
              ve_tot(:,:,:,im)=ve_tot(:,:,:,im)+ve_eq
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     convert the phi component of the perturbed B to cylindrical.
c-----------------------------------------------------------------------
      be_eq(3,:,:)=be_eq(3,:,:)/bigr
c-----------------------------------------------------------------------
c     when number density is evolved and used in the advance, include
c     its n=0 part in the computation of kinetic energy for each
c     component.  remaining bits from non-symmetric number density
c     are added to the n=0 part when continuity=full.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO im=1,nmodes
          IF (keff(im)==0) nd_n0=nd(:,:,:,im)
        ENDDO
        IF (continuity=='full') THEN
          DO im=1,nmodes
            IF (keff(im)==0) nd(:,:,:,im)=nd(:,:,:,im)+nd_eq
          ENDDO
          CALL fft_nim('inverse',mps,mps,lphi,1_i4,nd,real_nd,dealiase)
          DO im=1,nmodes
            IF (keff(im)==0) nd(:,:,:,im)=nd(:,:,:,im)-nd_eq
          ENDDO
        ENDIF
      ENDIF
      IF (nonlinear.AND.
     $    (continuity=='full'.OR.continuity=='n=0 only')) THEN
        nd_sym=nd_eq(1,:,:)+nd_n0(1,:,:)
      ELSE
        nd_sym=nd_eq(1,:,:)
      ENDIF
    
      presi_eq=pres_eq(1,:,:)-prese_eq(1,:,:)
c-----------------------------------------------------------------------
c     save energy densities in the cell-index vertex.
c-----------------------------------------------------------------------
      mode_loop: DO im=1,nmodes
        IF (keff(im)==0) THEN
          IF (nonlinear) THEN
            int(1,:,1,im)=SUM(0.5*wj
     $                       *SUM((be(:,:,:,im)+be_eq)
     $                       *CONJG((be(:,:,:,im)+be_eq)),1)/mu0,1)
            IF (eq_flow/='none') THEN
              int(2,:,1,im)=SUM(0.5*nd_sym*mtot*wj
     $                         *SUM((ve(:,:,:,im)+ve_eq)
     $                         *CONJG((ve(:,:,:,im)+ve_eq)),1),1)
            ELSE
              int(2,:,1,im)=SUM(0.5*nd_sym*mtot*wj
     $                       *SUM(ve(:,:,:,im)*CONJG(ve(:,:,:,im)),1),1)
            ENDIF
            IF (beta>0) THEN
              int(3,:,1,im)=SUM(wj*(kboltz*
     $                          (tele(1,:,:,im)*nd_eq(1,:,:)+
     $                             nd(1,:,:,im)*te_eq(1,:,:)+
     $                           tele(1,:,:,im)*   nd(1,:,:,im))+
     $                         prese_eq(1,:,:))/gamm1,1)
              int(4,:,1,im)=SUM(wj*(kboltz/zeff*
     $                          (tion(1,:,:,im)*nd_eq(1,:,:)+
     $                             nd(1,:,:,im)*ti_eq(1,:,:)+
     $                           tion(1,:,:,im)*   nd(1,:,:,im))+
     $                         presi_eq)/gamm1,1)
            ELSE
              int(3:4,:,1,im)=0._r8
            ENDIF
c-----------------------------------------------------------------------
c         add 2*equilibrium*(n=0) for linear runs. 
c
c         here, nd_sym is just nd_eq.
c-----------------------------------------------------------------------
          ELSE
            int(1,:,1,im)=SUM(wj*
     $        (0.5*SUM(be(:,:,:,im)*CONJG(be(:,:,:,im)),1)
     $            +SUM(REAL(be(:,:,:,im))*be_eq,1))/mu0,1)
            IF (eq_flow/='none') THEN
              int(2,:,1,im)=SUM(wj*0.5*mtot*( (nd_sym+nd(1,:,:,im))*
     $                                  SUM((ve(:,:,:,im)+ve_eq)**2,1)-
     $                                  nd_sym*SUM(ve_eq**2,1) ),1)
            ELSE
              int(2,:,1,im)=SUM(wj*0.5*mtot*(nd_sym+nd(1,:,:,im))*
     $                        SUM(ve(:,:,:,im)*CONJG(ve(:,:,:,im)),1),1)
            ENDIF
            IF (beta>0) THEN
              int(3,:,1,im)=SUM(wj*kboltz*
     $                         (tele(1,:,:,im)*nd_eq(1,:,:)+
     $                            nd(1,:,:,im)*te_eq(1,:,:)+
     $                          tele(1,:,:,im)*   nd(1,:,:,im))/gamm1,1)
              int(4,:,1,im)=SUM(wj*kboltz/zeff*
     $                         (tion(1,:,:,im)*nd_eq(1,:,:)+
     $                            nd(1,:,:,im)*ti_eq(1,:,:)+
     $                          tion(1,:,:,im)*   nd(1,:,:,im))/gamm1,1)
            ELSE
              int(3:4,:,1,im)=0._r8
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c       n>0 Fourier components.  the factors of 2 take the complex
c       conjugate into account.
c-----------------------------------------------------------------------
        ELSE
          int(1,:,1,im)=
     $      SUM(wj*SUM(be(:,:,:,im)*CONJG(be(:,:,:,im)),1)/mu0,1)
          int(2,:,1,im)=
     $      SUM(wj*SUM(ve(:,:,:,im)*CONJG(ve(:,:,:,im)),1)*
     $      nd_sym*mtot,1)
          IF (beta>0._r8) THEN
            int(3,:,1,im)=SUM(wj*kboltz*
     $                      (tele(1,:,:,im)*CONJG(nd(1,:,:,im))+
     $                      CONJG(tele(1,:,:,im))*nd(1,:,:,im))/gamm1,1)
            int(4,:,1,im)=SUM(wj*kboltz/zeff*
     $                      (tion(1,:,:,im)*CONJG(nd(1,:,:,im))+
     $                      CONJG(tion(1,:,:,im))*nd(1,:,:,im))/gamm1,1)
          ELSE
            int(3:4,:,1,im)=0._r8
          ENDIF
        ENDIF
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     for nonlinear computations with continuity==full find the 
c     kinetic energy including number density perturbations.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.continuity=='full') THEN
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,ve_tot,real_ve,
     $               dealiase)
        real_ke(1,:,:)=0.5*mtot*
     $                 SUM(real_ve(:,:,:)**2,1)*real_nd(1,:,:)
        CALL fft_nim('forward',mps,mps,lphi,1_i4,kin_e,real_ke,dealiase)
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            int(2,:,1,im)=SUM(wj*kin_e(1,:,:,im),1)
            im0=im
            EXIT
          ENDIF
        ENDDO
        int(2,:,1,im0)=int(2,:,1,im0)-
     $                   SUM(int(2,:,1,:im0-1),2)-
     $                   SUM(int(2,:,1,im0+1:),2)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE energy_density
c-----------------------------------------------------------------------
c     subprogram 4. divb_int.
c     compute the integrand for div(b).
c-----------------------------------------------------------------------
      SUBROUTINE divb_int(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: divb
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ber,bez
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      INTEGER(i4) :: im,iv,nv
c-----------------------------------------------------------------------
c     set number of vertices per cell.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     evaluate the perturbed magnetic field.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,ber,bez,1_i4)
c-----------------------------------------------------------------------
c     find div(b) for this mode.
c-----------------------------------------------------------------------
      mode_loop: DO im=1,nmodes
        divb(:,:)=ber(1,:,:,im)+bez(2,:,:,im)
     $           +(0,1)*keff(im)*be(3,:,:,im)/bigr
        IF (geom=='tor') THEN
          divb(:,:)=divb(:,:)+be(1,:,:,im)/bigr
        ENDIF
c-----------------------------------------------------------------------
c       save alpha*div(b).
c-----------------------------------------------------------------------
        DO iv=1,nv
          int(1,:,iv,im)=SUM(alpha(:,:,iv)*divb(:,:),1)
        ENDDO
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE divb_int
c-----------------------------------------------------------------------
c     subprogram 5. jpar_int.
c     compute the integrand for mu0*J.B/B**2.
c-----------------------------------------------------------------------
      SUBROUTINE jpar_int(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: jpar
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ber,bez,ja
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be_eq,ja_eq
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: be_eq_mag2,
     $          ja_eq_dot_be_eq
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: real_be,
     $          real_ja
      REAL(r8), DIMENSION(1,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: real_jpar
      REAL(r8), DIMENSION(1,1,1) :: dv
      INTEGER(i4) :: im,iv,nv,iq,is,ivind,mps
c-----------------------------------------------------------------------
c     set number of vertices per cell.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
      mps=SIZE(bigr,1)*SIZE(bigr,2)
c-----------------------------------------------------------------------
c     evaluate the fields.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8)
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%ja_eq,tb%ja_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ja_eq,dv,dv,0_i4)
c-----------------------------------------------------------------------
c     convert the phi components of equilibrium fields to cylindrical.
c-----------------------------------------------------------------------
      be_eq(3,:,:)=be_eq(3,:,:)/bigr
      ja_eq(3,:,:)=ja_eq(3,:,:)*bigr
      ja_eq=mu0*ja_eq
c-----------------------------------------------------------------------
c     form the ratio directly for nonlinear computations:
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            be(:,:,:,im)=be(:,:,:,im)+be_eq
            ja(:,:,:,im)=ja(:,:,:,im)+ja_eq
          ENDIF
        ENDDO
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,be,real_be,dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,ja,real_ja,dealiase)
        real_jpar(1,:,:)=SUM(real_ja*real_be,1)/SUM(real_be**2,1)
        CALL fft_nim('forward',mps,mps,lphi,1_i4,jpar,real_jpar,
     $               dealiase)
c-----------------------------------------------------------------------
c     form the ratio from first-order terms for linear computations:
c-----------------------------------------------------------------------
      ELSE
        be_eq_mag2=SUM(be_eq**2,1)+TINY(mu0)
        ja_eq_dot_be_eq=SUM(ja_eq*be_eq,1)
        DO im=1,nmodes
          jpar(1,:,:,im)=(   SUM(ja_eq*be(:,:,:,im)
     $                          +be_eq*ja(:,:,:,im),1)
     $                    -2*SUM(be_eq*be(:,:,:,im),1)
     $                      *ja_eq_dot_be_eq/be_eq_mag2 ) / be_eq_mag2
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     save alpha*mu0*J.B/B**2
c-----------------------------------------------------------------------
      DO im=1,nmodes
        DO iv=1,nv
          int(1,:,iv,im)=SUM(alpha(:,:,iv)*jpar(1,:,:,im),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jpar_int
c-----------------------------------------------------------------------
c     subprogram 6. edotj_int.
c     compute the integrand for <-vXb>.<J> and eta*J**2
c-----------------------------------------------------------------------
      SUBROUTINE edotj_int(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ber,bez,ja,ve
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             prese,preser,presez,te
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: bcrossv,
     $             jcrossb
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $          eidotjp,eidotjt,ehdotjp,ehdotjt
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: ja_eq,be_eq,
     $          ve_eq,grad_prese_eq
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: diff,nd_eq,
     $          pe_eq,pe_eqr,pe_eqz,te_eq
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: etaj2,
     $          etaj2p,etaj2t
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kfac
      INTEGER(i4) :: im,iv,nv,iq,is,ivind
c-----------------------------------------------------------------------
c     set number of vertices per cell.
c-----------------------------------------------------------------------
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     evaluate the fields.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1/mu0)
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%ja_eq,tb%ja_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ja_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%ve,tb%ve,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,dc,dc,0_i4)
      IF (eq_flow/='none') THEN
        CALL generic_all_eval(rb%ve_eq,tb%ve_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        ve_eq,dv,dv,0_i4)
      ELSE
        ve_eq=0
      ENDIF
      IF (ohms/='mhd') THEN
        CALL generic_all_eval(rb%nd_eq,tb%nd_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        nd_eq,dv,dv,0_i4)
        IF (beta>0) THEN
          CALL generic_all_eval(rb%prese,tb%prese,rb%dxdr,rb%dydr,
     $                          rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                          tb%ng,prese,preser,presez,1_i4)
          CALL generic_all_eval(rb%prese_eq,tb%prese_eq,rb%dxdr,
     $                          rb%dydr,rb%dxdz,rb%dydz,rb%xg,rb%yg,
     $                          tb%tgeom,tb%ng,pe_eq,pe_eqr,pe_eqz,1_i4)
        ELSE
          preser=0
          presez=0
          pe_eqr=0
          pe_eqz=0
        ENDIF
      ELSE
        ehdotjp=0
        ehdotjt=0
      ENDIF
c-----------------------------------------------------------------------
c     electrical diffusivity may depend on Te, but the full 3D variation
c     is not considered.
c-----------------------------------------------------------------------
      IF (eta_model=='eta n=0 only'.OR.eta_model=='eta full') THEN
        CALL generic_all_eval(rb%tele_eq,tb%tele_eq,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,te_eq,dv,dv,0_i4)
        CALL generic_all_eval(rb%tele,tb%tele,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        te,dc,dc,0_i4)
        diff=MAX(0._r8,MIN(elecd_max,elecd*
     $       (eta_ref_t/MAX(smallnum,(te_eq+REAL(te(:,:,:,1)))))**1.5 ))
      ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
        CALL generic_all_eval(rb%diff_shape,tb%diff_shape,rb%dxdr,
     $                        rb%dydr,rb%dxdz,rb%dydz,rb%xg,rb%yg,
     $                        tb%tgeom,tb%ng,diff,dv,dv,0_i4)
        diff=elecd*diff
      ELSE
        diff=elecd
      ENDIF
c-----------------------------------------------------------------------
c     convert the phi comp of equilibrium fields to cylindrical. then,
c     add the n=0 comp of j to find <J>, and save this in both ja_eq
c     and the n=0 comp of j.
c-----------------------------------------------------------------------
      be_eq(3,:,:)=be_eq(3,:,:)/bigr
      ja_eq(3,:,:)=ja_eq(3,:,:)*bigr
      DO im=1,nmodes
        IF (keff(im)==0) THEN
          ja(:,:,:,im)=REAL(ja(:,:,:,im))
          be(:,:,:,im)=REAL(be(:,:,:,im))
          ve(:,:,:,im)=REAL(ve(:,:,:,im))
          IF (nonlinear) THEN
            ja_eq=ja_eq+ja(:,:,:,im)
            ja(:,:,:,im)=ja_eq
            be(:,:,:,im)=be_eq+be(:,:,:,im)
            ve(:,:,:,im)=ve_eq+ve(:,:,:,im)
          ENDIF
          EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     find the quasilinear product from each Fourier component directly.
c-PRE with T-dependent resistivity, eta is computed from T_n=0 only;
c     otherwise, the Fourier decomposition becomes complicated.
c-----------------------------------------------------------------------
      etaj2=0
      etaj2p=mu0*diff(1,:,:)*SUM(ja_eq(1:2,:,:)**2,1)
      etaj2t=mu0*diff(1,:,:)*ja_eq(3,:,:)**2
      DO im=1,nmodes
        IF (keff(im)==0) THEN
          kfac=1
        ELSE
          kfac=2
        ENDIF
c-----------------------------------------------------------------------
c       eta*j**2
c-----------------------------------------------------------------------
        etaj2=etaj2+kfac*mu0*diff(1,:,:)
     $                  *SUM(ja(:,:,:,im)*CONJG(ja(:,:,:,im)),1)
c-----------------------------------------------------------------------
c       -<vxb>.<J>
c-----------------------------------------------------------------------
        CALL math_cart_cross(bcrossv,be(:,:,:,im),
     $                       CONJG(ve(:,:,:,im)),1._r8)
        eidotjp(:,:,im)=kfac*SUM(bcrossv(1:2,:,:)*ja_eq(1:2,:,:),1)
        eidotjt(:,:,im)=kfac*bcrossv(3,:,:)*ja_eq(3,:,:)
c-----------------------------------------------------------------------
c       <jxb>.<J>/ne
c-----------------------------------------------------------------------
        IF (ohms/='mhd') THEN
          CALL math_cart_cross(jcrossb,ja(:,:,:,im),
     $                         CONJG(be(:,:,:,im)),1._r8)
          IF (keff(im)==0) THEN
            jcrossb(1,:,:)=jcrossb(1,:,:)-preser(1,:,:,im)
            jcrossb(2,:,:)=jcrossb(2,:,:)-presez(1,:,:,im)
            IF (nonlinear) THEN
              jcrossb(1,:,:)=jcrossb(1,:,:)-pe_eqr(1,:,:)
              jcrossb(2,:,:)=jcrossb(2,:,:)-pe_eqz(1,:,:)
            ENDIF
          ENDIF
          ehdotjp(:,:,im)=kfac*SUM(jcrossb(1:2,:,:)*ja_eq(1:2,:,:),1)
     $      *(1-meomi)/((1+meomi)*nd_eq(1,:,:)*elementary_q)
          ehdotjt(:,:,im)=kfac*jcrossb(3,:,:)*ja_eq(3,:,:)
     $      *(1-meomi)/((1+meomi)*nd_eq(1,:,:)*elementary_q)
        ENDIF

      ENDDO
c-----------------------------------------------------------------------
c     scatter the E.J contributions to the vertices.  note that Fourier
c     component storage in int is not standard.
c-----------------------------------------------------------------------
      DO im=1,nmodes
        DO iv=1,nv
          int(         im,:,iv)=SUM(
     $      alpha(:,:,iv)*eidotjp(:,:,im),1)
          int(  nmodes+im,:,iv)=SUM(
     $      alpha(:,:,iv)*eidotjt(:,:,im),1)
          int(2*nmodes+im,:,iv)=SUM(
     $      alpha(:,:,iv)*ehdotjp(:,:,im),1)
          int(3*nmodes+im,:,iv)=SUM(
     $      alpha(:,:,iv)*ehdotjt(:,:,im),1)
        ENDDO
      ENDDO
      DO iv=1,nv
        int(4*nmodes+1,:,iv)=SUM(alpha(:,:,iv)*etaj2p,1)
        int(4*nmodes+2,:,iv)=SUM(alpha(:,:,iv)*etaj2t,1)
        int(4*nmodes+3,:,iv)=SUM(alpha(:,:,iv)*etaj2,1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE edotj_int
c-----------------------------------------------------------------------
c     subprogram 7. heat_flux_int.
c     compute the integrand for the conductive part of toroidally 
c     averaged heat flux from the selected Fourier component of T.
c
c-PRE the present computation assumes separate_pe=F, i.e. a single-
c     temperature computation.
c-----------------------------------------------------------------------
      SUBROUTINE heat_flux_int(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: nd_eq,nd_eqr,
     $          nd_eqz,ti_eq,ti_eqr,ti_eqz,diff
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be_eq,ve_eq
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq,
     $          q_para_n,q_perp_n,chiprp
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: 
     $          real_work,real_be,bbgradt
      REAL(r8), DIMENSION(1,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: 
     $          real_scal,real_ndtot,chipar
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: nd,
     $             ti,tir,tiz,bdgr_ti
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             grad_ti,q_para,ve,convect
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,imode,iq,jq,ncx,ncy,mps,im0,ip
      REAL(r8) :: jfac
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      mps=ncx*ncy
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     set number density.  the n=0 component will have nd_eq added to 
c     it, and nd_eq will have nd_n=0 added if continuity=='n=0 only'.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%nd,tb%nd,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      nd,dc,dc,0_i4)
      CALL generic_all_eval(rb%nd_eq,tb%nd_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      nd_eq,nd_eqr,nd_eqz,1_i4)

      IF (nonlinear) THEN
        DO imode=1,nmodes
          IF (keff(imode)==0) THEN
            nd(:,:,:,imode)=nd(:,:,:,imode)+nd_eq
            IF (continuity=='n=0 only') THEN
              nd_eq=nd(:,:,:,imode)
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     set temperature and grad(T).  zero out all but jmode if jmode
c     represents the single chosen Fourier component.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%tion_eq,tb%tion_eq,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ti_eq,ti_eqr,ti_eqz,1_i4)
      CALL generic_all_eval(rb%tion,tb%tion,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ti,tir,tiz,1_i4)
      DO imode=1,nmodes
        IF (jmode>-1.AND.imode/=jmode) THEN
          grad_ti(:,:,:,imode)=0
        ELSE
          grad_ti(1,:,:,imode)=tir(1,:,:,imode)
          grad_ti(2,:,:,imode)=tiz(1,:,:,imode)
          grad_ti(3,:,:,imode)=(0,1)*keff(imode)*ti(1,:,:,imode)/bigr
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     for anisotropic thermal conduction, get the
c     magnetic field first.  note that we assume k_pll and k_perp have
c     the factor of (gamma-1) in them already.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,dc,dc,0_i4)
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be_eq,dv,dv,0_i4)
      be_eq(3,:,:)=be_eq(3,:,:)/bigr
c-----------------------------------------------------------------------
c     get flow velocity to compute the convective heat flux.
c-----------------------------------------------------------------------
      convect=0
      CALL generic_all_eval(rb%ve,tb%ve,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,dc,dc,0_i4)
      IF (eq_flow/='none')
     $  CALL generic_all_eval(rb%ve_eq,tb%ve_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        ve_eq,dv,dv,0_i4)
      DO imode=1,nmodes
        IF (jmode>-1.AND.imode/=jmode) ve(:,:,:,imode)=0
        IF (keff(imode)==0) THEN
          im0=imode
          q_para_n=ti(1,:,:,imode)
          ti(:,:,:,imode)=ti(:,:,:,imode)+ti_eq
        ENDIF
      ENDDO

      CALL fft_nim('inverse',mps,mps,lphi,1_i4,ti,real_scal,dealiase)
      CALL fft_nim('inverse',mps,mps,lphi,1_i4,nd,real_ndtot,dealiase)
      CALL fft_nim('inverse',mps,mps,lphi,3_i4,ve,real_work,dealiase)
      DO jq=1,3
        real_work(jq,:,:)=kboltz*(1+1/zeff)*real_work(jq,:,:)*
     $                    real_ndtot(1,:,:)*real_scal( 1,:,:)/gamm1
      ENDDO
      real_scal=kboltz*(1+1/zeff)*real_ndtot*real_scal/gamm1
      CALL fft_nim('forward',mps,mps,lphi,3_i4,convect,real_work,
     $             dealiase)
      CALL fft_nim('forward',mps,mps,lphi,1_i4,bdgr_ti,real_scal,
     $             dealiase)

      DO imode=1,nmodes
        IF (keff(imode)==0) ti(1,:,:,imode)=q_para_n
        IF (eq_flow/='none') THEN
          DO jq=1,3
            convect(jq,:,:,imode)=convect(jq,:,:,imode)+
     $              ve_eq(jq,:,:)*bdgr_ti( 1,:,:,imode)
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     determine diffusivity coefficients according to p_model.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%diff_shape,tb%diff_shape,rb%dxdr,
     $                      rb%dydr,rb%dxdz,rb%dydz,rb%xg,rb%yg,
     $                      tb%tgeom,tb%ng,diff,dv,dv,0_i4)
      IF (nonlinear.AND.p_model(1:6)=='aniso_') THEN

        DO imode=1,nmodes
          IF (keff(imode)==0) THEN
            im0=imode
            q_para_n=ti(1,:,:,imode)
            ti(:,:,:,imode)=ti(:,:,:,imode)+ti_eq
          ENDIF
        ENDDO

        CALL fft_nim('inverse',mps,mps,lphi,1_i4,ti,real_scal,dealiase)
        chipar=MAX(k_pll_min, MIN(k_pll_max, k_plli*
     $    (MAX(smallnum,real_scal)/k_pll_ref_t)**2.5 ))

        IF (p_model=='aniso_tdep') THEN
          chiprp=k_perpi*
     $           SQRT(MAX(1._r8+(kprp_mnrat*(diff(1,:,:)-1)/dvac)**2,
     $             k_pll_ref_t/REAL(ti(1,:,:,im0))))/
     $           SUM((be_eq+REAL(be(:,:,:,im0)))**2,1)
        ELSE
          chiprp=k_perpi
        ENDIF
        DO imode=1,nmodes
          IF (keff(imode)==0) ti(1,:,:,imode)=q_para_n
        ENDDO
      ELSE
        chipar=k_plle
        chiprp=k_perpi
      ENDIF
c-----------------------------------------------------------------------
c     For nonlinear calculations, the parallel tensor is 
c     B_tot B_tot / B^{2}. 
c     note that beq dot grad Teq is assumed to be zero, a slightly
c     stronger assumption than just a steady state relation.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO imode=1,nmodes
          bdgr_ti(1,:,:,imode)=be(1,:,:,imode)*ti_eqr(1,:,:)
     $                        +be(2,:,:,imode)*ti_eqz(1,:,:)
        ENDDO

        DO imode=1,nmodes
          IF (keff(imode)==0) THEN 
            be(:,:,:,imode)=be(:,:,:,imode)+be_eq
            bmagsq=SUM(be(:,:,:,imode)**2,1)
          ENDIF
        ENDDO

        CALL fft_nim('inverse',mps,mps,lphi,1_i4,bdgr_ti,real_scal,
     $               dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,be,real_be,dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,grad_ti,real_work,
     $               dealiase)
c-----------------------------------------------------------------------
c       Build the product B (B dot grad T).  we are adding 
c       (be_eq+be).grad(T) to be.grad(T_eq).
c-----------------------------------------------------------------------
        real_scal(1,:,:)=(real_scal(1,:,:)+SUM(real_be*real_work,1))
        bbgradt(1,:,:)=real_be(1,:,:)*real_scal(1,:,:)
        bbgradt(2,:,:)=real_be(2,:,:)*real_scal(1,:,:)
        bbgradt(3,:,:)=real_be(3,:,:)*real_scal(1,:,:)
c-----------------------------------------------------------------------
c	Convert the above vector to the Fourier representation
c	and store in q_para.  the denominator, B**2 is approximated
c       as (be_n0+be_eq)**2.
c-----------------------------------------------------------------------
        CALL fft_nim('forward',mps,mps,lphi,3_i4,q_para,bbgradt,
     $               dealiase)
        DO imode=1,nmodes
          q_para(1,:,:,imode)=q_para(1,:,:,imode)/bmagsq
          q_para(2,:,:,imode)=q_para(2,:,:,imode)/bmagsq
          q_para(3,:,:,imode)=q_para(3,:,:,imode)/bmagsq
        ENDDO
c-----------------------------------------------------------------------
c     For linear calculations, the parallel tensor is divided 
c     into linear contributions with B^{2} approximated as Beq^{2}.
c-----------------------------------------------------------------------
      ELSE
        bmagsq=SUM(be_eq*be_eq,1)
        DO imode=1,nmodes
          bdgr_ti(1,:,:,imode)=
     $      (be(1,:,:,imode)*ti_eqr(1,:,:)
     $      +be(2,:,:,imode)*ti_eqz(1,:,:)
     $      +SUM(grad_ti(:,:,:,imode)*be_eq,1))/bmagsq
          q_para(1,:,:,imode)=q_para(1,:,:,imode)*bdgr_ti(1,:,:,imode)
          q_para(2,:,:,imode)=q_para(2,:,:,imode)*bdgr_ti(1,:,:,imode)
          q_para(3,:,:,imode)=q_para(3,:,:,imode)*bdgr_ti(1,:,:,imode)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     separate parallel and perpendicular contributions for diagnostic
c     purposes.  the q_para array will hold the parallel part, while
c     grad_ti will hold the perpendicular part.  note that a factor
c     of gamma-1 is already in the chis, so it needs to be removed
c     to get heat flux.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        DO jq=1,3
          grad_ti(jq,:,:,imode)=-kboltz*(1+1/zeff)*chiprp*
     $        (grad_ti(jq,:,:,imode)-q_para(jq,:,:,imode))/gamm1
        ENDDO
      ENDDO
      q_para=-kboltz*(1+1/zeff)*q_para/gamm1
      IF (nonlinear.AND.p_model(1:6)=='aniso_') THEN
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,q_para,real_work,
     $               dealiase)
        DO jq=1,3
          real_work(jq,:,:)=real_work(jq,:,:)*chipar(1,:,:)
        ENDDO
        CALL fft_nim('forward',mps,mps,lphi,3_i4,q_para,real_work,
     $               dealiase)
      ELSE
        q_para=q_para*k_plle
      ENDIF
c-----------------------------------------------------------------------
c     sum thermal diffusion from isotropic and anisotropic terms and
c     multiply by number density if anisotropic conduction is used.
c-----------------------------------------------------------------------
      IF (.NOT.nonlinear.OR.continuity/='full') THEN
        DO imode=1,nmodes
          DO jq=1,3
            grad_ti(jq,:,:,imode)=grad_ti(jq,:,:,imode)*nd_eq(1,:,:)
            q_para(jq,:,:,imode)=q_para(jq,:,:,imode)*nd_eq(1,:,:)
          ENDDO
        ENDDO
      ELSE
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,grad_ti,real_work,
     $               dealiase)
        DO jq=1,3
          real_work(jq,:,:)=real_work(jq,:,:)*real_ndtot(1,:,:)
        ENDDO
        CALL fft_nim('forward',mps,mps,lphi,3_i4,grad_ti,real_work,
     $               dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,q_para,real_work,
     $               dealiase)
        DO jq=1,3
          real_work(jq,:,:)=real_work(jq,:,:)*real_ndtot(1,:,:)
        ENDDO
        CALL fft_nim('forward',mps,mps,lphi,3_i4,q_para,real_work,
     $               dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     find the component normal to the equilibrium + n=0 magnetic
c     field.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO imode=1,nmodes
          IF (keff(imode)/=0) CYCLE
          q_para_n=(q_para(1,:,:,imode)*be(2,:,:,imode)-
     $              q_para(2,:,:,imode)*be(1,:,:,imode))/SQRT(bmagsq)
          q_perp_n=(grad_ti(1,:,:,imode)*be(2,:,:,imode)-
     $              grad_ti(2,:,:,imode)*be(1,:,:,imode))/SQRT(bmagsq)
        ENDDO
      ELSE
        q_para_n=0
        q_perp_n=0
      ENDIF
c-----------------------------------------------------------------------
c     form the contributions to each basis function.  parallel
c     contributions are the first 3 components, and perpendicular
c     are the next 3.  only the toroidally averaged part of the
c     heat flux is plotted.
c-----------------------------------------------------------------------
      int=0
      DO imode=1,nmodes
        IF (keff(imode)/=0) CYCLE
        DO iv=1,nv
          int(1,:,iv)=SUM(alpha(:,:,iv)*q_para(1,:,:,imode),1)
          int(2,:,iv)=SUM(alpha(:,:,iv)*q_para(2,:,:,imode),1)
          int(3,:,iv)=SUM(alpha(:,:,iv)*q_para(3,:,:,imode),1)
          int(4,:,iv)=SUM(alpha(:,:,iv)*q_para_n,1)
          int(5,:,iv)=SUM(alpha(:,:,iv)*grad_ti(1,:,:,imode),1)
          int(6,:,iv)=SUM(alpha(:,:,iv)*grad_ti(2,:,:,imode),1)
          int(7,:,iv)=SUM(alpha(:,:,iv)*grad_ti(3,:,:,imode),1)
          int(8,:,iv)=SUM(alpha(:,:,iv)*q_perp_n,1)
          int(9,:,iv)=SUM(alpha(:,:,iv)*convect(1,:,:,imode),1)
          int(10,:,iv)=SUM(alpha(:,:,iv)*convect(2,:,:,imode),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE heat_flux_int
c-----------------------------------------------------------------------
c     subprogram 8. mach_int.
c     compute the integrand for the ion acoustic speed, mach number,
c     and parallel mach number for the selected Fourier component. 
c-----------------------------------------------------------------------
      SUBROUTINE mach_int(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be_eq,ve_eq
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: ti_eq,te_eq
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2),nphi) :: real_ve,
     $          real_be
      REAL(r8), DIMENSION(1,SIZE(bigr,1)*SIZE(bigr,2),nphi) ::
     $          real_tion,real_tele,real_machn,real_ve_tot,
     $          real_ion_aspd,real_parmach 
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ve
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             tele,tion,machn,ion_aspd,parmach
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: im,iv,nv,ncx,ncy,mps
      REAL(r8) :: cmin
c-----------------------------------------------------------------------
c     convenience parameters
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      mps=ncx*ncy
c-----------------------------------------------------------------------
c     evaluate the perturbed and equilibrium fields.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_all_eval(rb%tele,tb%tele,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      tele,dc,dc,0_i4)
      CALL generic_all_eval(rb%tion,tb%tion,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      tion,dc,dc,0_i4)
      CALL generic_all_eval(rb%tele_eq,tb%tele_eq,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                      tb%ng,te_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%tion_eq,tb%tion_eq,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                      tb%ng,ti_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%ve,tb%ve,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,dc,dc,0_i4)
      IF (eq_flow/='none')
     $  CALL generic_all_eval(rb%ve_eq,tb%ve_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                        rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                        ve_eq,dv,dv,0_i4)
      CALL generic_all_eval(rb%be,tb%be,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,dc,dc,0_i4)
      CALL generic_all_eval(rb%be_eq,tb%be_eq,rb%dxdr,rb%dydr,rb%dxdz,
     $                      rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be_eq,dv,dv,0_i4)
c-----------------------------------------------------------------------
c     convert the phi component of the perturbed B to cylindrical.
c-----------------------------------------------------------------------
      be_eq(3,:,:)=be_eq(3,:,:)/bigr
c-----------------------------------------------------------------------
c     add in equilibrium values
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            tele(:,:,:,im)=tele(:,:,:,im)+te_eq
            tion(:,:,:,im)=tion(:,:,:,im)+ti_eq
            IF (eq_flow/='none') ve(:,:,:,im)=ve(:,:,:,im)+ve_eq
            be(:,:,:,im)=be(:,:,:,im)+be_eq
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       fourier transform every mode, calculate mach number quantities 
c       and fourier transform back.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',mps,mps,lphi,1_i4,tele,real_tele,
     $               dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,1_i4,tion,real_tion,
     $               dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,ve,real_ve,
     $               dealiase)
        CALL fft_nim('inverse',mps,mps,lphi,3_i4,be,real_be,
     $               dealiase)
        real_parmach(1,:,:)=SUM(real_be*real_ve,1)/
     $                      MAX(SQRT(SUM(real_be**2,1)),smallnum)
        real_ve_tot(1,:,:)=SQRT(SUM(real_ve**2,1))
        real_ion_aspd=SQRT(MAX( zeff*kboltz*real_tele
     $                        +gamma*kboltz*real_tion,smallnum)/ms(2))
        cmin=1.e-8_r8*MAXVAL(real_ion_aspd)
        WHERE(real_ion_aspd<2._r8*cmin) real_ve_tot=0._r8
        WHERE(real_ion_aspd<2._r8*cmin) real_parmach=0._r8
        real_machn=real_ve_tot/real_ion_aspd
        real_parmach=real_parmach/real_ion_aspd
        CALL fft_nim('forward',mps,mps,lphi,1_i4,machn,real_machn,
     $               dealiase)
        CALL fft_nim('forward',mps,mps,lphi,1_i4,ion_aspd,real_ion_aspd,
     $               dealiase)
        CALL fft_nim('forward',mps,mps,lphi,1_i4,parmach,real_parmach,
     $               dealiase)
c-----------------------------------------------------------------------
c     for linear cases, just plot information from the equilibrium.
c-----------------------------------------------------------------------
      ELSE
        ion_aspd=0._r8 
        machn=0._r8 
        parmach=0._r8 
        ion_aspd(:,:,:,1)=SQRT(MAX(zeff*kboltz*te_eq
     $                           +gamma*kboltz*ti_eq,smallnum)/ms(2))
        cmin=1.e-8*MAXVAL(REAL(ion_aspd(:,:,:,1)))
        IF (eq_flow/='none') THEN
          machn(1,:,:,1)=SQRT(SUM(ve_eq**2,1))/
     $      MAX(cmin,REAL(ion_aspd(1,:,:,1)))
          parmach(1,:,:,1)=SUM(ve_eq*be_eq,1)/
     $      (MAX(SQRT(SUM(real_be**2,1)),smallnum)*
     $       MAX(cmin,REAL(ion_aspd(1,:,:,1))))
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     sum to get integral for each quantity
c-----------------------------------------------------------------------
      DO im=1,nmodes
        DO iv=1,nv
          int(1,:,iv,im)=SUM(alpha(:,:,iv)*ion_aspd(1,:,:,im),1)
          int(2,:,iv,im)=SUM(alpha(:,:,iv)*machn(1,:,:,im),1)
          int(3,:,iv,im)=SUM(alpha(:,:,iv)*parmach(1,:,:,im),1) 
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mach_int
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimplot_ints
