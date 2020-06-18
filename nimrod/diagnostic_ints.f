c-----------------------------------------------------------------------
c     file diagnostic_ints.f
c     module that includes integrand routines associated with physics
c     diagnostics.  all subroutines included here must use the same
c     interface block as those in integrands.f
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  div_b.
c     2.  energy_density.
c     3.  eq_i_phi.
c     4.  n0_i_phi.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE diagnostic_ints

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
c     subprogram 1. div_b.
c     compute the integrand for div(b)**2 and b**2.  this routine is
c     used for the div(b) diagnostic only and should not be confused
c     with the above routines for cleaning div(b) in integrands.f.
c-----------------------------------------------------------------------
      SUBROUTINE div_b(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq,wj
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kfac
      COMPLEX(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: divb
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: be,ber,bez
      INTEGER(i4) :: imode
c-----------------------------------------------------------------------
c     evaluate the perturbed magnetic field.
c
c     the sum of the test functions produces the quadrature weight times
c     the Jacobian.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
c-----------------------------------------------------------------------
c     begin mode loop.
c-----------------------------------------------------------------------
      int=0
      mode_loop: DO imode=1,nmodes
        IF (keff(imode)==0) THEN
          kfac=1
        ELSE
          kfac=2
        ENDIF
c-----------------------------------------------------------------------
c       find div(b) for this mode.
c-----------------------------------------------------------------------
        IF (geom=='tor') THEN
          divb=ber(1,:,:,imode)+bez(2,:,:,imode)
     $            +((0,1)*keff(imode)*be(3,:,:,imode)
     $                               +be(1,:,:,imode))/bigr
        ELSE
          divb=ber(1,:,:,imode)+bez(2,:,:,imode)
     $             +(0,1)*keff(imode)*be(3,:,:,imode)/bigr
        ENDIF
c-----------------------------------------------------------------------
c       find the real and imaginary contributions to b**2 for this mode.
c-----------------------------------------------------------------------
        bmagsq(:,:)=SUM(be(:,:,:,imode)*CONJG(be(:,:,:,imode)),1)
c-----------------------------------------------------------------------
c       save [div(b)]**2 and b**2 in the cell-index vertex.
c-----------------------------------------------------------------------
        int(1,:,1)=int(1,:,1)+SUM(wj*kfac*divb*CONJG(divb),1)
        int(2,:,1)=int(2,:,1)+SUM(wj*kfac*bmagsq,1)
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE div_b
c-----------------------------------------------------------------------
c     subprogram 2. energy_density.
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
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          be_eq,ve_eq,nd_eq,pres_eq,
     $          prese_eq,nd_n0,real_ndptr,ti_eq,te_eq,dp
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: nd_sym,presi_eq,
     $          wj
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_ke
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kfac,pfac

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,ve,tele,tion,nd,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             ve_tot
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             kin_e

      INTEGER(i4) :: im,ncx,ncy
c-----------------------------------------------------------------------
c     convenience parameters
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
c-----------------------------------------------------------------------
c     evaluate the perturbed and equilibrium fields.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,dcp,dcp,0_i4)
      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,dcp,dcp,0_i4)
      IF (continuity=='none') THEN
        ALLOCATE(nd(1,ncx,ncy,nmodes))
        nd=0._r8
      ELSE
        CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,nd,
     $                       dcp,dcp,0_i4)
      ENDIF
      ve_tot=ve
      IF (beta>0) THEN
        CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,inode,tele,dcp,
     $                       dcp,0_i4)
        CALL generic_ptr_set(rb%qtion,tb%qtion,tb%tgeom,inode,tion,dcp,
     $                       dcp,0_i4)
        CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                       inode,te_eq,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qtion_eq,tb%qtion_eq,tb%tgeom,
     $                       inode,ti_eq,dp,dp,0_i4)
      ENDIF
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                     inode,pres_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qprese_eq,tb%qprese_eq,tb%tgeom,
     $                     inode,prese_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      IF (eq_flow/='none') THEN
        CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
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
c     when number density is evolved and used in the advance, include
c     its n=0 part in the computation of kinetic energy for each
c     component.  remaining bits from non-symmetric number density
c     are added to the n=0 part when continuity=full.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
        IF (continuity=='full') THEN
          CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                         inode,real_ndptr,dp,dp,0_i4)
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
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve_tot,real_ve,
     $               dealiase)
        real_ke(1,:,:)=0.5*mtot*
     $                 SUM(real_ve(:,:,:)**2,1)*real_ndptr(1,:,:)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,kin_e,real_ke,
     $               dealiase)
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            int(2,:,1,im)=SUM(wj*kin_e(1,:,:,im),1)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF (continuity=='none') DEALLOCATE(nd)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE energy_density
c-----------------------------------------------------------------------
c     subprogram 3. eq_i_phi.
c     compute the plasma current and toroidal flux for _eq fields.
c-----------------------------------------------------------------------
      SUBROUTINE eq_i_phi(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: wj
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          be_eq,ja_eq,dp
c-----------------------------------------------------------------------
c     evaluate the equilibrium fields
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,ja_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     accumulate plasma current and toroidal flux in the non-conforming
c     vertex.  the jacobian in the rblock and tblock rhs routines will
c     eliminate the 1/R factor to give flux.
c-----------------------------------------------------------------------
      int(1,:,1)=SUM(wj*ja_eq(3,:,:)/bigr,1)
      int(2,:,1)=SUM(wj*be_eq(3,:,:)/bigr,1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE eq_i_phi
c-----------------------------------------------------------------------
c     subprogram 4. n0_i_phi.
c     compute the plasma current and toroidal flux for n=0 Fourier
c     component.  (only computed on layer 0)
c-----------------------------------------------------------------------
      SUBROUTINE n0_i_phi(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: wj
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: be,ber,bez
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ja
c-----------------------------------------------------------------------
c     evaluate the equilibrium fields
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1/mu0)
c-----------------------------------------------------------------------
c     place plasma current and toroidal flux in the non-conforming
c     vertex.  the jacobian in the rblock and tblock rhs routines will
c     eliminate the 1/R factor to give flux.
c-----------------------------------------------------------------------
      int(1,:,1)=SUM(wj*ja(3,:,:,1)/bigr,1)
      int(2,:,1)=SUM(wj*be(3,:,:,1)/bigr,1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE n0_i_phi
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE diagnostic_ints
