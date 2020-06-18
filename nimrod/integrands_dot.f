c-----------------------------------------------------------------------
c     file integrands.f
c     module that includes integrand routines used by matrix-free
c     solvers for computing the effect of multiplying a matrix
c     by a vector.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module integrands_dot
c     subprogram 1. b_3deta_dot.
c     subprogram 2. b_hmhd_dot.
c     subprogram 3. t_aniso_dot.
c     subprogram 4. v_aniso_dot.
c     subprogram 5. v_3dsi_dot.
c     subprogram 6. cont_dot.
c-----------------------------------------------------------------------
c     module containing integrands for the dot-product computations.
c-----------------------------------------------------------------------
      MODULE integrands_dot

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
c     subprogram 1. b_3deta_dot.
c     compute the action of the 3d resistive diffusion matrix (due to
c     3d resistivity) on a vector stored in work4.
c-----------------------------------------------------------------------
      SUBROUTINE b_3deta_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: elecd_phi,dp
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_e
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: be,
     $             ber,bez,dtcurlb
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2),nmodes) :: divb
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,nv,npol,ncx,ncy,b_derivs
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the magnetic
c     field iterate (from work4).  multiply the curl by feta*dt but not
c     1/mu0.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      CALL generic_all_eval(rb%work4,tb%work4,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,dtcurlb,feta*dt)
      CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                     inode,elecd_phi,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find div(b) from the old B if cleaning is not split.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb) THEN
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            divb(:,:,imode)=fdivb*divbd*dt*
     $        ( ber(1,:,:,imode)+bez(2,:,:,imode)
     $        +((0,1)*keff(imode)*be(3,:,:,imode)+be(1,:,:,imode))/bigr)
          ELSE
            divb(:,:,imode)=fdivb*divbd*dt*
     $        ( ber(1,:,:,imode)+bez(2,:,:,imode)
     $         +(0,1)*keff(imode)*be(3,:,:,imode)/bigr )
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     compute the resistive electric field.  dtcurlb becomes 
c     f*dt*eta*j.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,dtcurlb,real_e,
     $             dealiase)
      real_e(1,:,:)=elecd_phi(1,:,:)*real_e(1,:,:)
      real_e(2,:,:)=elecd_phi(1,:,:)*real_e(2,:,:)
      real_e(3,:,:)=elecd_phi(1,:,:)*real_e(3,:,:)
      CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,dtcurlb,real_e,
     $             dealiase)
c-----------------------------------------------------------------------
c     find the contributions that multiply the basis vector itself.
c     part of this is just the iterate (from the d/dt term).
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb) THEN
        DO imode=1,nmodes
          be(1,:,:,imode)=be(1,:,:,imode)-
     $                    (0,1)*keff(imode)/bigr*dtcurlb(2,:,:,imode)
          be(2,:,:,imode)=be(2,:,:,imode)+
     $                    (0,1)*keff(imode)/bigr*dtcurlb(1,:,:,imode)
          be(3,:,:,imode)=be(3,:,:,imode)-
     $                    (0,1)*keff(imode)*divb(:,:,imode)/bigr
        ENDDO
      ELSE
        DO imode=1,nmodes
          be(1,:,:,imode)=be(1,:,:,imode)-
     $                    (0,1)*keff(imode)/bigr*dtcurlb(2,:,:,imode)
          be(2,:,:,imode)=be(2,:,:,imode)+
     $                    (0,1)*keff(imode)/bigr*dtcurlb(1,:,:,imode)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     construct B+dt*curl(alpha^*).E [ +dt*div(alpha)*div(B) ].
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb) THEN
        DO imode=1,nmodes
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(-dalpdz (:,:,iv)*dtcurlb(3,:,:,imode)
     $                            +alpha  (:,:,iv)*be(1,:,:,imode)
     $                            +dalpdrc(:,:,iv)*divb(:,:,imode),1)
            int(2,:,iv,imode)=SUM( dalpdr(:,:,iv)*dtcurlb(3,:,:,imode)
     $                            +alpha (:,:,iv)*be(2,:,:,imode)
     $                            +dalpdz(:,:,iv)*divb(:,:,imode),1)
            int(3,:,iv,imode)=SUM( dalpdz (:,:,iv)*dtcurlb(1,:,:,imode)
     $                            -dalpdrc(:,:,iv)*dtcurlb(2,:,:,imode)
     $                            +alpha  (:,:,iv)*be(3,:,:,imode),1)
          ENDDO
        ENDDO
      ELSE
        DO imode=1,nmodes
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(-dalpdz (:,:,iv)*dtcurlb(3,:,:,imode)
     $                            +alpha  (:,:,iv)*be(1,:,:,imode),1)
            int(2,:,iv,imode)=SUM( dalpdr(:,:,iv)*dtcurlb(3,:,:,imode)
     $                            +alpha (:,:,iv)*be(2,:,:,imode),1)
            int(3,:,iv,imode)=SUM( dalpdz (:,:,iv)*dtcurlb(1,:,:,imode)
     $                          -dalpdrc(:,:,iv)*dtcurlb(2,:,:,imode)
     $                            +alpha(:,:,iv)*be(3,:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE b_3deta_dot
c-----------------------------------------------------------------------
c     subprogram 2. b_hmhd_dot.
c     compute the action of the implicit HMHD operator on the iterate
c     for the change in magnetic field (work4).
c
c     the terms that appear in this routine are only the ones that 
c     depend on changes in B over the time-split.
c
c     this has all such terms, including those that do not couple
c     Fourier components.
c
c       de=[dj X (Beq+b) + (Jeq+j) X db]/ne - (Veq+v) X db + eta*dj
c
c-----------------------------------------------------------------------
      SUBROUTINE b_hmhd_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          almod,dalmr,dalmz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          real_be,real_ja,
     $          nd_eq,nd_n0,real_ndptr,real_ve,ds2,elecd_n0,elecd_phi,dp
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: diff
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_dbe,real_de,real_tmp
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: hfac,hdt,fdb,fe,disdbfac

      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: aux
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             de,dja,dbe,dber,dbez,crlv2
      COMPLEX(r8), DIMENSION(6,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             d6v,d6vr,d6vz
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2),nmodes) :: divb
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: divbim,dv2
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,nv,ncx,ncy,in0,nvc,nvm
      REAL(r8) :: fhyp1,fhyp2,fhdb
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the 0-th
c     order magnetic field.  evaluate the magnetic field and current
c     density at the start of the time step, and evaluate the change
c     in each associated with the iterate of the matrix-free solve.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      IF (hyp_eta<=0._r8.AND.hyp_dbd<=0._r8.OR.split_hypeta) THEN
        CALL generic_all_eval(rb%work4,tb%work4,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,dbe,dber,dbez,1_i4)
      ELSE
        CALL generic_all_eval(rb%w6v2,tb%w6v2,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,d6v,d6vr,d6vz,1_i4)
        fhyp1=SQRT(hyp_eta*fhyp_eta*dt)
        fhyp2=fhyp1*mu0
        fhdb=SQRT(hyp_dbd*fhyp_dbd*dt)
        dbe =d6v (4:6,:,:,:)
        dber=d6vr(4:6,:,:,:)
        dbez=d6vz(4:6,:,:,:)
        CALL math_curl(nmodes,keff,geom,bigr,dbe,dber,dbez,
     $                 crlv2,fhyp1)
        dbe =d6v (1:3,:,:,:)
        dber=d6vr(1:3,:,:,:)
        dbez=d6vz(1:3,:,:,:)
      ENDIF
      CALL math_curl(nmodes,keff,geom,bigr,dbe,dber,dbez,dja,1._r8/mu0)

      CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,inode,
     $                     real_ve,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qja_tot,tb%qja_tot,tb%tgeom,inode,
     $                     real_ja,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,inode,
     $                     real_be,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     dissipation coefficients.
c-----------------------------------------------------------------------
      IF (elecd>0) THEN
        IF (threedeta) THEN
          CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                         inode,elecd_phi,dp,dp,0_i4)
          fe=dt*feta*mu0
        ELSE IF (eta_model=='eta n=0 only') THEN
          CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                         inode,elecd_n0,dp,dp,0_i4)
          diff=dt*feta*mu0*elecd_n0(1,:,:)
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          diff=dt*feta*mu0*elecd*MAX(ds2(1,:,:),0._r8)
        ELSE
          diff=dt*feta*mu0*elecd
        ENDIF
      ENDIF
      IF (split_divb.OR.divbd<0) THEN
        fdb=0._r8
      ELSE
        fdb=dt*fdivb*divbd
      ENDIF
c-----------------------------------------------------------------------
c     find div(b) for cleaning if not split.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb.OR.poly_divb>=0) THEN
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            divb(:,:,imode)=(dber(1,:,:,imode)+dbez(2,:,:,imode)
     $                     +((0,1)*keff(imode)* dbe(3,:,:,imode)
     $                                         +dbe(1,:,:,imode))/bigr)
          ELSE
            divb(:,:,imode)=(dber(1,:,:,imode)+dbez(2,:,:,imode)
     $                      +(0,1)*keff(imode)* dbe(3,:,:,imode)/bigr)
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the implicit contributions to the electric field that
c     couple Fourier components.
c
c     first find perturbed db and dj in real space, then compute
c     0.5*dt*(dj X B + J X db)/e, and transform back to Fourier
c     coefficients.  note that B and J include the equilibrium and
c     old-time-step perturbed fields.
c
c     also add the partial(J)/partial(t) term here.
c-----------------------------------------------------------------------
      hdt=0.5_r8*dt
      hfac=hdt*coefhll
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,dbe,real_dbe,
     $             dealiase)

      IF (ohms/='mhd') THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,dja,real_tmp,
     $               dealiase)
        CALL math_cart_cross(real_de,real_tmp,real_be,hfac)
        CALL math_cadd_cross(real_de,real_ja,real_dbe,hfac)
        IF (ohms=='2fl') THEN
          real_de=real_de+real_tmp*coefme1
        ENDIF
c-----------------------------------------------------------------------
c       all of the preceeding terms must be divided by number density.
c       the n that is used depends on the continuity input parameter.
c-----------------------------------------------------------------------
        IF (continuity/='full')
     $    CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,de,real_de,
     $                 dealiase)
        CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                       inode,nd_eq,dp,dp,0_i4)
        IF (continuity=='none'.OR.continuity=='fix profile') THEN
          DO imode=1,nmodes
            de(1,:,:,imode)=de(1,:,:,imode)/nd_eq(1,:,:)
            de(2,:,:,imode)=de(2,:,:,imode)/nd_eq(1,:,:)
            de(3,:,:,imode)=de(3,:,:,imode)/nd_eq(1,:,:)
          ENDDO
c-----------------------------------------------------------------------
c       use symmetric part of evolving density.
c-----------------------------------------------------------------------
        ELSE IF (continuity=='n=0 only') THEN
          CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                         nd_n0,dp,dp,0_i4)
          DO imode=1,nmodes
            de(1,:,:,imode)=de(1,:,:,imode)/(nd_eq(1,:,:)+nd_n0(1,:,:))
            de(2,:,:,imode)=de(2,:,:,imode)/(nd_eq(1,:,:)+nd_n0(1,:,:))
            de(3,:,:,imode)=de(3,:,:,imode)/(nd_eq(1,:,:)+nd_n0(1,:,:))
          ENDDO
c-----------------------------------------------------------------------
c       use 3D evolving density but do not forward transform until the
c       ideal term is added.  note that this no longer goes
c       back and forth to Fourier representation between the JxB
c       product and the division by n.
c-----------------------------------------------------------------------
        ELSE
          CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                         inode,real_ndptr,dp,dp,0_i4)
          real_de(1,:,:)=real_de(1,:,:)/real_ndptr(1,:,:)
          real_de(2,:,:)=real_de(2,:,:)/real_ndptr(1,:,:)
          real_de(3,:,:)=real_de(3,:,:)/real_ndptr(1,:,:)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     construct terms that do not depend on number density starting
c     with the ideal electric field, 0.5*dt * db X V.
c     [if continuity is not full, dber is used for temporary space.]
c     real_ve includes any equilibrium COM flow.
c
c     3D resistivity is also here to minimize the number of FFTs.
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN
        IF (threedeta) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,dja,real_tmp,
     $                 dealiase)
          real_tmp(1,:,:)=real_tmp(1,:,:)*elecd_phi(1,:,:)*fe
          real_tmp(2,:,:)=real_tmp(2,:,:)*elecd_phi(1,:,:)*fe
          real_tmp(3,:,:)=real_tmp(3,:,:)*elecd_phi(1,:,:)*fe
          CALL math_cadd_cross(real_tmp,real_dbe,real_ve,hdt)
        ELSE
          CALL math_cart_cross(real_tmp,real_dbe,real_ve,hdt)
        ENDIF
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,de,real_tmp,
     $               dealiase)

      ELSE IF (continuity=='full') THEN
        IF (threedeta) THEN
          real_de(1,:,:)=real_de(1,:,:)+
     $                   real_tmp(1,:,:)*elecd_phi(1,:,:)*fe
          real_de(2,:,:)=real_de(2,:,:)+
     $                   real_tmp(2,:,:)*elecd_phi(1,:,:)*fe
          real_de(3,:,:)=real_de(3,:,:)+
     $                   real_tmp(3,:,:)*elecd_phi(1,:,:)*fe
        ENDIF
        CALL math_cadd_cross(real_de,real_dbe,real_ve,hdt)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,de,real_de,
     $               dealiase)

      ELSE
        IF (threedeta) THEN
          real_tmp(1,:,:)=real_tmp(1,:,:)*elecd_phi(1,:,:)*fe
          real_tmp(2,:,:)=real_tmp(2,:,:)*elecd_phi(1,:,:)*fe
          real_tmp(3,:,:)=real_tmp(3,:,:)*elecd_phi(1,:,:)*fe
          CALL math_cadd_cross(real_tmp,real_dbe,real_ve,hdt)
        ELSE
          CALL math_cart_cross(real_tmp,real_dbe,real_ve,hdt)
        ENDIF
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,dber,real_tmp,
     $               dealiase)
        de=de+dber
      ENDIF
c-----------------------------------------------------------------------
c     resistive term for 2D eta_models.
c-----------------------------------------------------------------------
      IF (elecd>0.AND..NOT.threedeta) THEN
        DO imode=1,nmodes
          de(1,:,:,imode)=de(1,:,:,imode)+diff*dja(1,:,:,imode)
          de(2,:,:,imode)=de(2,:,:,imode)+diff*dja(2,:,:,imode)
          de(3,:,:,imode)=de(3,:,:,imode)+diff*dja(3,:,:,imode)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the contribution for the discontinuous scalar equation,
c     disdbfac*almod*div(delta-b), where almod is the value of the
c     discontinuous basis function and mwork2 holds div(delta-b).
c-----------------------------------------------------------------------
      IF (poly_divb>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',almod,dalmr,
     $                          dalmz,0_i4,poly_divb,
     $                          polydmin=poly_divb_min,
     $                          polydmax=poly_divb_max)
        nvm=SIZE(almod,3)
        CALL generic_all_eval(rb%mwork2,tb%mwork2,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,aux,dc,dc,0_i4)

        disdbfac=SQRT(disc_dbd*dt*fdivb)
      ELSE
        aux=0._r8
        disdbfac=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     construct +curl(A*).de_implicit plus div(dB) corrections
c-----------------------------------------------------------------------
      IF (hyp_eta<=0._r8.AND.hyp_dbd<=0._r8.OR.split_hypeta) THEN
        DO imode=1,nmodes
          divbim=fdb*divb(:,:,imode)+disdbfac*aux(1,:,:,imode)
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(   -dalpdz(:,:,iv)*de(3,:,:,imode)
     $         -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*de(2,:,:,imode)
     $             +dalpdrc(:,:,iv)*divbim
     $         +alpha(:,:,iv)*dbe(1,:,:,imode),1)
            int(2,:,iv,imode)=SUM(    dalpdr(:,:,iv)*de(3,:,:,imode)
     $         +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*de(1,:,:,imode)
     $             +dalpdz (:,:,iv)*divbim
     $         +alpha(:,:,iv)*dbe(2,:,:,imode),1)
            int(3,:,iv,imode)=SUM( dalpdz (:,:,iv)*de(1,:,:,imode)
     $                            -dalpdrc(:,:,iv)*de(2,:,:,imode)
     $                      -(0,1)*alpha(:,:,iv)*keff(imode)/bigr
     $                             *divbim
     $         +alpha(:,:,iv)*dbe(3,:,:,imode),1)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     duplicate version with contributions from the hyper-resistivity
c     auxiliary vector.  note that dja is multiplied by fhyp2 and
c     the hyper-resistive term is combined with de.  the divergence
c     of the auxiliary vector is computed one Fourier component at a
c     time.
c-----------------------------------------------------------------------
      ELSE
        DO imode=1,nmodes
          divbim=fdb*divb(:,:,imode)+disdbfac*aux(1,:,:,imode)
          dja(:,:,:,imode)=dja(:,:,:,imode)*fhyp2
          de(:,:,:,imode) =de (:,:,:,imode)-crlv2(:,:,:,imode)

          IF (geom=='tor') THEN
            divbim=fdb*divb(:,:,imode)+disdbfac*aux(1,:,:,imode)-
     $        fhdb*(d6vr(4,:,:,imode)+d6vz(5,:,:,imode)
     $             +((0,1)*keff(imode)*d6v(6,:,:,imode)
     $                                +d6v(4,:,:,imode))/bigr)
          ELSE
            divbim=fdb*divb(:,:,imode)+disdbfac*aux(1,:,:,imode)-
     $        fhdb*(d6vr(4,:,:,imode)+d6vz(5,:,:,imode)
     $             +(0,1)*keff(imode)* d6v(6,:,:,imode)/bigr)
          ENDIF

          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(   -dalpdz(:,:,iv)*de(3,:,:,imode)
     $         -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*de(2,:,:,imode)
     $             +dalpdrc(:,:,iv)*divbim
     $         +alpha(:,:,iv)*d6v(1,:,:,imode),1)
            int(2,:,iv,imode)=SUM(    dalpdr(:,:,iv)*de(3,:,:,imode)
     $         +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*de(1,:,:,imode)
     $             +dalpdz (:,:,iv)*divbim
     $         +alpha(:,:,iv)*d6v(2,:,:,imode),1)
            int(3,:,iv,imode)=SUM( dalpdz (:,:,iv)*de(1,:,:,imode)
     $                            -dalpdrc(:,:,iv)*de(2,:,:,imode)
     $         -(0,1)*alpha(:,:,iv)*keff(imode)*divbim/bigr
     $         +alpha(:,:,iv)*d6v(3,:,:,imode),1)

            int(4,:,iv,imode)=SUM(alpha(:,:,iv)*d6v(4,:,:,imode)-
     $         dalpdz(:,:,iv)*dja(3,:,:,imode)
     $        -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*dja(2,:,:,imode)
     $        +fhdb*dalpdrc(:,:,iv)*divb(:,:,imode),1)
            int(5,:,iv,imode)=SUM(alpha(:,:,iv)*d6v(5,:,:,imode)+
     $         dalpdr(:,:,iv)*dja(3,:,:,imode)
     $        +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*dja(1,:,:,imode)
     $        +fhdb*dalpdz(:,:,iv)*divb(:,:,imode),1)
            int(6,:,iv,imode)=SUM(alpha(:,:,iv)*d6v(6,:,:,imode)+
     $          dalpdz (:,:,iv)*dja(1,:,:,imode)
     $         -dalpdrc(:,:,iv)*dja(2,:,:,imode)
     $         -fhdb*(0,1)*keff(imode)*alpha(:,:,iv)*
     $               divb(:,:,imode)/bigr,1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete the integrand for the discontinuous scalar equation.
c-----------------------------------------------------------------------
      IF (poly_divb>=0) THEN
        divb=disdbfac*divb
        DO imode=1,nmodes
          DO iv=nvc+1,nvc+nvm
            int(1,:,iv,imode)=SUM(almod(:,:,iv-nvc)*
     $                            (aux(1,:,:,imode)-divb(:,:,imode)),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE b_hmhd_dot
c-----------------------------------------------------------------------
c     subprogram 3. t_aniso_dot.
c
c     compute the action of the anisotropic thermal conduction operator
c     (plus identity term) on a vector:
c
c     (ns*I-dt*f*D).dTs
c
c     where D is the conduction operator, f is a centering parameter,
c     dTs is an iterate, and n is number density.
c
c     when implicit advection is used, the operator is
c
c     (ns*I + 0.5*dt*ns*( V.grad(I)+(gamma-1)*div(V)*I ) - dt*f*D).dTs
c
c     this is now used for nonlinear computations only.
c-----------------------------------------------------------------------
      SUBROUTINE t_aniso_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          be_eq,b0,nd_eq,nd_n0,
     $          real_bptr,real_ndptr,kappl,kaprp,real_jptr,real_vptr,
     $          ndqr,ndqz,real_grdv,upwc,ndiff,dp
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: kprp_1d
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_vect,bbgradt,real_grn
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scl,real_t
      REAL(r8), DIMENSION(mpseudo,nphi) :: magBsq3D
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kdiff,zz,kpll,kprp,jfac,fdt,hdt,fupw,fupwi,gm1
      REAL(r8), DIMENSION(3) :: rvec
      REAL(r8) :: rscl2,rdot

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,nd,ndr,ndz,dc
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             grad_tie,grad_nd,be_tot
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: tie,
     $             tier,tiez,bdgr_tie

      INTEGER(i4) :: nv,iv,imode,ncx,ncy,ip,ix
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      fdt=fthc*dt
      hdt=0.5_r8*dt
      IF (integrand_flag(1:6)=='all ti') THEN
        zz=zeff
        IF (separate_pe) THEN
          kpll=k_plli
          jfac=hdt*coefjvi/zeff
        ELSE
          kpll=k_plle
          jfac=0._r8
        ENDIF
        kprp=k_perpi
      ELSE
        zz=1._r8
        kpll=k_plle
        kprp=k_perpe
        jfac=hdt*coefjve
      ENDIF
      IF (p_model=='isothermal') THEN
        gm1=0._r8
      ELSE
        gm1=gamm1
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     number densities as required.
c-PRE faster to form J.grad(ne) at quadrature points.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,ndqr,ndqz,1_i4)
      CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                     inode,real_ndptr,dp,dp,0_i4)
      IF (continuity=='n=0 only')
     $  CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
      IF (continuity/='none'.AND.impladv)
     $  CALL generic_ptr_set(rb%qndiff_phi,tb%qndiff_phi,tb%tgeom,
     $                       inode,ndiff,dp,dp,0_i4)
      IF (impladv.AND.separate_pe) THEN
        CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,
     $                       inode,nd,ndr,ndz,1_i4)
        DO imode=1,nmodes
          grad_nd(1,:,:,imode)=ndr(1,:,:,imode)
          grad_nd(2,:,:,imode)=ndz(1,:,:,imode)
          grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
          IF (keff(imode)==0) THEN
            grad_nd(1,:,:,imode)=grad_nd(1,:,:,imode)+ndqr(1,:,:)
            grad_nd(2,:,:,imode)=grad_nd(2,:,:,imode)+ndqz(1,:,:)
          ENDIF
        ENDDO
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $               real_grn,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     set dT_s and grad(dT_s) with the iterate stored in work2.  it
c     is used only once, so quadrature storage is not used--call
c     all_eval not ptr_set.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%work2,tb%work2,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      tie,tier,tiez,1_i4)
      DO imode=1,nmodes
        grad_tie(1,:,:,imode)=tier(1,:,:,imode)
        grad_tie(2,:,:,imode)=tiez(1,:,:,imode)
        grad_tie(3,:,:,imode)=(0,1)*keff(imode)*tie(1,:,:,imode)/bigr
      ENDDO
c-----------------------------------------------------------------------
c     terms associated with implicit advection and the n*dT/dt term
c     are here.  note that tie becomes
c
c     n_s*( dT_s + dt/2*(V_s.grad(dT_s) + (gamma-1)*n_s*div(V_s)*dT_s) )
c
c     time-centered implicit advection, n_s*V_s and n_s*div(V_s) for
c     each species s.  the Jeq.grad(neq)=0 approximation is used here.
c
c     for ions, we need
c
c       (dt/2)*n_i*V_i = (dt/2)*(n_e*V/z + J*meomi/(z*e*(1+meomi)) )
c
c     where meomi = z*me/mi.  jfac is the full coefficient for J.
c     we also need
c
c       (dt/2)*n_i*div(V_i) =
c         (dt/2)*(n_e*div(V)/z - J.grad(n_e)*meomi/(n_e*z*e*(1+meomi))
c
c     the electron terms have z->1 and the meomi in the numerator
c     replaced by -1.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,tie,real_t,
     $             dealiase)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_tie,
     $             real_vect,dealiase)

      IF (impladv) THEN
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,real_vptr,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                       inode,real_grdv,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qja_tot,tb%qja_tot,tb%tgeom,
     $                       inode,real_jptr,dp,dp,0_i4)

        IF (continuity=='fix profile') THEN
          real_scl(1,:,:)=
     $      hdt*(SUM(real_vptr*real_vect,1)+
     $           gm1*real_t(1,:,:)*(real_grdv(1,:,:)+real_grdv(5,:,:)+
     $                              real_grdv(9,:,:)) ) +
     $      fdt*ndiff(1,:,:)*real_t(1,:,:)/zz
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                 bdgr_tie,real_scl,dealiase)
          DO imode=1,nmodes
            tie(:,:,:,imode)=
     $        (tie(:,:,:,imode)+bdgr_tie(:,:,:,imode))*nd_eq/zz
          ENDDO
          IF (jfac/=0._r8) THEN
            real_scl(1,:,:)=jfac*(
     $        SUM(real_jptr*real_vect,1)-
     $        SUM(real_jptr*real_grn ,1)*
     $          gm1*real_t(1,:,:)/real_ndptr(1,:,:) )
            CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                   bdgr_tie,real_scl,dealiase)
            tie=tie+bdgr_tie
          ENDIF

        ELSE IF (continuity=='n=0 only') THEN
          real_scl(1,:,:)=
     $      hdt*(SUM(real_vptr*real_vect,1)+
     $           gm1*real_t(1,:,:)*(real_grdv(1,:,:)+real_grdv(5,:,:)+
     $                              real_grdv(9,:,:)) ) +
     $      fdt*ndiff(1,:,:)*real_t(1,:,:)/zz
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                 bdgr_tie,real_scl,dealiase)
          DO imode=1,nmodes
            tie(:,:,:,imode)=
     $        (tie(:,:,:,imode)+bdgr_tie(:,:,:,imode))*(nd_eq+nd_n0)/zz
          ENDDO
          IF (jfac/=0._r8) THEN
            real_scl(1,:,:)=jfac*(
     $        SUM(real_jptr*real_vect,1)-
     $        SUM(real_jptr*real_grn ,1)*
     $          gm1*real_t(1,:,:)/real_ndptr(1,:,:) )
            CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                   bdgr_tie,real_scl,dealiase)
            tie=tie+bdgr_tie
          ENDIF

        ELSE ! continuity = full
          real_scl(1,:,:)=real_ndptr(1,:,:)*( real_t(1,:,:)+
     $      hdt*(SUM(real_vptr*real_vect,1)+
     $           gm1*real_t(1,:,:)*(real_grdv(1,:,:)+real_grdv(5,:,:)+
     $                              real_grdv(9,:,:)) ) )/zz +
     $      fdt*ndiff(1,:,:)*real_t(1,:,:)/zz
          IF (jfac/=0._r8) THEN
            real_scl(1,:,:)=real_scl(1,:,:) + jfac*(
     $        SUM(real_jptr*real_vect,1)-
     $        SUM(real_jptr*real_grn ,1)*
     $          gm1*real_t(1,:,:)/real_ndptr(1,:,:) )
          ENDIF
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tie,real_scl,
     $                 dealiase)
        ENDIF

      ELSE
        IF (continuity=='fix profile') THEN
          DO imode=1,nmodes
            tie(:,:,:,imode)=tie(:,:,:,imode)*nd_eq/zz
          ENDDO
        ELSE IF (continuity=='n=0 only') THEN
          DO imode=1,nmodes
            tie(:,:,:,imode)=tie(:,:,:,imode)*(nd_eq+nd_n0)/zz
          ENDDO
        ELSE ! continuity = full
          real_t=real_t*real_ndptr/zz
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tie,real_t,
     $                 dealiase)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field.  note that we assume k_pll and k_perp have
c     the factor of (gamma-1) in them already.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,dc,dc,0_i4)
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      IF (p_model=="aniso_plltdep".OR.p_model=="aniso_tdep") THEN
        IF (integrand_flag(1:6)=='all ti') THEN
          CALL generic_ptr_set(rb%qkappli_phi,tb%qkappli_phi,tb%tgeom,
     $                         inode,kappl,dp,dp,0_i4)
          IF (.NOT.closure_n0_only) THEN
            CALL generic_ptr_set(rb%qkaprpi_phi,tb%qkaprpi_phi,tb%tgeom,
     $                           inode,kaprp,dp,dp,0_i4)
          ELSE
            CALL generic_ptr_set(rb%qkaprpi_n0,tb%qkaprpi_n0,tb%tgeom,
     $                           inode,kaprp,dp,dp,0_i4)
            kprp_1d=RESHAPE(kaprp,(/ncx*ncy/))
          ENDIF
        ELSE
          CALL generic_ptr_set(rb%qkapple_phi,tb%qkapple_phi,tb%tgeom,
     $                         inode,kappl,dp,dp,0_i4)
          IF (.NOT.closure_n0_only) THEN
            CALL generic_ptr_set(rb%qkaprpe_phi,tb%qkaprpe_phi,tb%tgeom,
     $                           inode,kaprp,dp,dp,0_i4)
          ELSE
            CALL generic_ptr_set(rb%qkaprpe_n0,tb%qkaprpe_n0,tb%tgeom,
     $                           inode,kaprp,dp,dp,0_i4)
            kprp_1d=RESHAPE(kaprp,(/ncx*ncy/))
          ENDIF
        ENDIF
      ELSE
        kdiff=kpll-kprp
      ENDIF
c-----------------------------------------------------------------------
c     For nonlinear calculations, the parallel tensor is 
c     B_tot B_tot / B^{2}. 
c     note that beq dot grad Teq is assumed to be zero.
c-----------------------------------------------------------------------
      IF (p_model(1:5)=='aniso') THEN
        CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                       inode,b0,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                       inode,real_bptr,dp,dp,0_i4)
c-----------------------------------------------------------------------
c       Build the product B (B dot grad T), (be_eq+be).grad(T) first.
c
c       if the thermal conductivity is temperature-dependent
c       multiply by its distribution in phi here.
c-----------------------------------------------------------------------
        IF (p_model=="aniso_plltdep".OR.p_model=="aniso_tdep") THEN
          IF (.NOT.closure_n0_only) THEN
            DO ip=1,nphi
              real_scl(1,:,ip)=(real_bptr(1,:,ip)*real_vect(1,:,ip)+
     $                          real_bptr(2,:,ip)*real_vect(2,:,ip)+
     $                          real_bptr(3,:,ip)*real_vect(3,:,ip))*
     $            MAX(0._r8,kappl(1,:,ip)-kaprp(1,:,ip))
            ENDDO
          ELSE
            DO ip=1,nphi
              real_scl(1,:,ip)=(real_bptr(1,:,ip)*real_vect(1,:,ip)+
     $                          real_bptr(2,:,ip)*real_vect(2,:,ip)+
     $                          real_bptr(3,:,ip)*real_vect(3,:,ip))*
     $            MAX(0._r8,kappl(1,:,ip)-kprp_1d(ipseust:ipseuen))
            ENDDO
          ENDIF
        ELSE
          real_scl(1,:,:)=kdiff*SUM(real_bptr*real_vect,1)
        ENDIF
        bbgradt(1,:,:)=real_bptr(1,:,:)*real_scl(1,:,:)
        bbgradt(2,:,:)=real_bptr(2,:,:)*real_scl(1,:,:)
        bbgradt(3,:,:)=real_bptr(3,:,:)*real_scl(1,:,:)
c-----------------------------------------------------------------------
c	Convert the above vector to the Fourier representation
c	and store in be_tot.  B**2 is approximated as (b_n0+be_eq)**2.
c
c       include fully 3D temperature_dependent perpendicular
c       thermal conduction if specified.
c-----------------------------------------------------------------------
        IF (.NOT.closure_n0_only) THEN
          magBsq3D(:,:)=SUM(real_bptr(:,:,:)**2,1)
          bbgradt(1,:,:)=bbgradt(1,:,:)/magBsq3D(:,:)
     $         +kaprp(1,:,:)*real_vect(1,:,:)
          bbgradt(2,:,:)=bbgradt(2,:,:)/magBsq3D(:,:)
     $         +kaprp(1,:,:)*real_vect(2,:,:)
          bbgradt(3,:,:)=bbgradt(3,:,:)/magBsq3D(:,:)
     $         +kaprp(1,:,:)*real_vect(3,:,:)
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $         bbgradt,dealiase)
        ELSE
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $                 bbgradt,dealiase)
          bmagsq=SUM((b0+be_eq)**2,1)
          DO imode=1,nmodes
            be_tot(1,:,:,imode)=be_tot(1,:,:,imode)/bmagsq
            be_tot(2,:,:,imode)=be_tot(2,:,:,imode)/bmagsq
            be_tot(3,:,:,imode)=be_tot(3,:,:,imode)/bmagsq
          ENDDO
        ENDIF
      ELSE
        be_tot=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     if the perpendicular thermal diffusivity is temperature-dependent,
c     multiply grad_t by its distribution.  also add parallel
c     contributions.
c
c     nonlinear temperature-dependent perpendicular thermal conduction
c     is now above.
c-----------------------------------------------------------------------
      IF (p_model=="aniso_tdep") THEN
        IF (.NOT.closure_n0_only) THEN
          grad_tie=dt*be_tot
        ELSE
          DO imode=1,nmodes
            grad_tie(1,:,:,imode)=fdt*
     $        (kaprp(1,:,:)*grad_tie(1,:,:,imode)+be_tot(1,:,:,imode))
            grad_tie(2,:,:,imode)=fdt*
     $        (kaprp(1,:,:)*grad_tie(2,:,:,imode)+be_tot(2,:,:,imode))
            grad_tie(3,:,:,imode)=fdt*
     $        (kaprp(1,:,:)*grad_tie(3,:,:,imode)+be_tot(3,:,:,imode))
          ENDDO
        ENDIF
      ELSE IF (p_model(1:5)=='aniso') THEN
        grad_tie=fdt*(kprp*grad_tie+be_tot)
      ELSE IF (p_model=='isotropic') THEN
        grad_tie=fdt*kprp*grad_tie
      ELSE
        grad_tie=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     find the upwinding-like artificial diffusion.  at this point
c     real_vect still holds grad_tie as a function of phi.  the
c     bbgradt array is used for temporary storage of species total
c     velocity and then the VV.grad(T) vector.
c
c-PRE reuse V.grad(T) computation earlier when continuity is unified.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        IF (integrand_flag(1:6)=='all ti') THEN
          CALL generic_ptr_set(rb%qupti_phi,tb%qupti_phi,tb%tgeom,
     $                         inode,upwc,dp,dp,0_i4)
        ELSE
          CALL generic_ptr_set(rb%qupte_phi,tb%qupte_phi,tb%tgeom,
     $                         inode,upwc,dp,dp,0_i4)
        ENDIF

        IF (jfac/=0._r8) THEN
          jfac=jfac*zz
          fupw =dt**2*t_dart_upw*upw_aniso
          fupwi=dt**2*t_dart_upw*(1._r8-upw_aniso)
          DO ip=1,nphi
            DO ix=1,mpseudo
              rvec(1)=hdt*real_vptr(1,ix,ip)+
     $               jfac*real_jptr(1,ix,ip)/real_ndptr(1,ix,ip)
              rvec(2)=hdt*real_vptr(2,ix,ip)+
     $               jfac*real_jptr(2,ix,ip)/real_ndptr(1,ix,ip)
              rvec(3)=hdt*real_vptr(3,ix,ip)+
     $               jfac*real_jptr(3,ix,ip)/real_ndptr(1,ix,ip)
              rdot=fupw*upwc(1,ix,ip)*
     $             SUM(rvec*real_vect(:,ix,ip)) /
     $           ( SUM(rvec*rvec) + smallnum )
              rscl2=fupwi*upwc(1,ix,ip)
              bbgradt(:,ix,ip)=rvec*rdot+real_vect(:,ix,ip)*rscl2
            ENDDO
          ENDDO
        ELSE
          fupw =dt**2*t_dart_upw*upw_aniso
          fupwi=dt**2*t_dart_upw*(1._r8-upw_aniso)
          DO ip=1,nphi
            DO ix=1,mpseudo
              rvec=real_vptr(:,ix,ip)
              rdot=fupw*upwc(1,ix,ip)*
     $             SUM(rvec*real_vect(:,ix,ip)) /
     $           ( SUM(rvec*rvec) + smallnum )
              rscl2=fupwi*upwc(1,ix,ip)
              bbgradt(:,ix,ip)=rvec*rdot+real_vect(:,ix,ip)*rscl2
            ENDDO
          ENDDO
        ENDIF

        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $               bbgradt,dealiase)
        grad_tie=grad_tie+be_tot
      ENDIF
c-----------------------------------------------------------------------
c     multiply by number density for all p_models.
c     the identity term also needs a factor of n in these cases.
c-----------------------------------------------------------------------
      IF (continuity=='fix profile') THEN
        DO imode=1,nmodes
          grad_tie(1,:,:,imode)=nd_eq(1,:,:)/zz*grad_tie(1,:,:,imode)
          grad_tie(2,:,:,imode)=nd_eq(1,:,:)/zz*grad_tie(2,:,:,imode)
          grad_tie(3,:,:,imode)=nd_eq(1,:,:)/zz*grad_tie(3,:,:,imode)
        ENDDO
      ELSE IF (continuity=='n=0 only') THEN
        DO imode=1,nmodes
          grad_tie(1,:,:,imode)=
     $      (nd_eq(1,:,:)+nd_n0(1,:,:))/zz*grad_tie(1,:,:,imode)
          grad_tie(2,:,:,imode)=
     $      (nd_eq(1,:,:)+nd_n0(1,:,:))/zz*grad_tie(2,:,:,imode)
          grad_tie(3,:,:,imode)=
     $      (nd_eq(1,:,:)+nd_n0(1,:,:))/zz*grad_tie(3,:,:,imode)
        ENDDO
      ELSE ! continuity = full
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_tie,
     $               real_vect,dealiase)
        real_vect(1,:,:)=real_vect(1,:,:)*real_ndptr(1,:,:)/zz
        real_vect(2,:,:)=real_vect(2,:,:)*real_ndptr(1,:,:)/zz
        real_vect(3,:,:)=real_vect(3,:,:)*real_ndptr(1,:,:)/zz
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,grad_tie,
     $               real_vect,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     sum thermal diffusion from isotropic and anisotropic terms,
c     then complete all contributions to each basis function.
c     the identity term and the gradt_phi term are combined before the
c     basis function loop for efficiency.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        grad_tie(3,:,:,imode)=tie(1,:,:,imode)
     $    -(0,1)*keff(imode)*grad_tie(3,:,:,imode)/bigr
        DO iv=1,nv
          int(1,:,iv,imode)=SUM( dalpdr(:,:,iv)*grad_tie(1,:,:,imode)
     $                          +dalpdz(:,:,iv)*grad_tie(2,:,:,imode)
     $                          + alpha(:,:,iv)*grad_tie(3,:,:,imode),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE t_aniso_dot
c-----------------------------------------------------------------------
c     subprogram 4. v_aniso_dot.
c     find the action of the (mass density*I + viscous) operator dotted
c     into a vector (work4) for conjugate gradient iterations.  this is
c     used when mass density has 3d variations and continuity is
c     set to 'full,' when parallel viscosity is used in a nonlinear
c     problem, or for nonlinear implicit advection.  In these cases, the
c     matrix couples Fourier components.
c
c     the semi-implicit terms are also computed here for the first
c     time since nimrod3_1.  the cost of finding matrices has exceeded
c     that of the extra dot computations.
c-----------------------------------------------------------------------
      SUBROUTINE v_aniso_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          almod,dalmr,dalmz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          real_ndptr,real_bptr,
     $          real_oldv,real_oldgv,beq,beqr,beqz,b0,b0r,b0z,nd_eq,
     $          nd_n0,j0,nl_pres,peq,peqr,peqz,p0,p0r,p0z,ti_eq,kappli,
     $          dp,real_titot,ndiff,ds2
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be,ber,bez,ja
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: pres,presr,
     $          presz
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: diff,n0,bigr2,
     $          iso,grd,bmagsq,sqptot
      REAL(r8), DIMENSION(3,3,SIZE(bigr,1),SIZE(bigr,2),ncontb) :: rcrl
      REAL(r8), DIMENSION(9,mpseudo,nphi) :: real_gvten
      REAL(r8), DIMENSION(6,mpseudo,nphi) :: real_pten
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_vec
      REAL(r8), DIMENSION(3,3) :: vtmpr,btmp,wdotr,bc_wdotr
      REAL(r8), DIMENSION(3) :: tmp1,tmp2

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: tion,dcp
      COMPLEX(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: aux
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ve,
     $             ver,vez,force
      COMPLEX(r8), DIMENSION(9,SIZE(bigr,1),SIZE(bigr,2),nmodes) ::
     $             gvten
      COMPLEX(r8), DIMENSION(6,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: pten
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: dbe,
     $             ja_cross_dbe,ve_cross_ja,forck
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: divv,
     $             ve_gpr,sqpaux,bvort
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: baux
      COMPLEX(r8), DIMENSION(3,3) :: vtmp,crl
      COMPLEX(r8) :: div

      INTEGER(i4) :: ncx,ncy,iv,imode,jt,ix,iy,ip,im,nvc,nvm
      REAL(r8) :: kin_coef,iso_coef,par_coef,gcoef,vcoef,dloc
      REAL(r8) :: rb2,divtmpr,rscal,hdt,g,ani,cj0,dt2,rr,dtdv,dtpv
      REAL(r8), PARAMETER :: third=1._r8/3._r8,twoth=2._r8/3._r8
      LOGICAL :: use_gv
c-----------------------------------------------------------------------
c     convenience
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      hdt=0.5_r8*dt
      bigr2=bigr**2
      dt2=dt*dt
      IF ((kin_visc>0.OR.iso_visc>0.OR.par_visc>0.OR.gyr_visc>0).AND.
     $ (integrand_flag(1:4)=='visc'.OR.integrand_flag(1:3)=='all')) THEN
        use_gv=.true.
      ELSE
        use_gv=.false.
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the
c     diffusivity shape function.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      IF ((ds_use=='kin_visc'.OR.ds_use=='both').AND.
     $    (kin_visc>0.OR.iso_visc>0)) THEN
        CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                       inode,ds2,dp,dp,0_i4)
        diff=MAX(ds2(1,:,:),0._r8)
      ELSE
        diff=1._r8
      ENDIF
      IF (par_visc>0.OR.gyr_visc>0) THEN
        CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                       inode,real_bptr,dp,dp,0_i4)
        par_coef=3._r8*dt*fvsc*mtot*par_visc
      ENDIF
      IF (gyr_visc>0) THEN
        CALL generic_ptr_set(rb%qtion,tb%qtion,tb%tgeom,
     $                       inode,tion,dcp,dcp,0_i4)
        CALL generic_ptr_set(rb%qtion_eq,tb%qtion_eq,tb%tgeom,
     $                       inode,ti_eq,dp,dp,0_i4)
      ENDIF
      kin_coef=dt*fvsc*mtot*kin_visc
      iso_coef=dt*fvsc*mtot*iso_visc
c-----------------------------------------------------------------------
c     evaluate the current direction vector for velocity change and
c     the number density.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%work4,tb%work4,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,ver,vez,1_i4)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $             dealiase)

      IF (continuity=='full') THEN
        CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)
      ELSE
        CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                       inode,nd_eq,dp,dp,0_i4)
      ENDIF
      IF (continuity=='n=0 only') THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
        n0=nd_eq(1,:,:)+nd_n0(1,:,:)
      ELSE IF (continuity/='full') THEN
        n0=nd_eq(1,:,:)
      ENDIF
      IF (impladv) THEN
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,real_oldv,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                       inode,real_oldgv,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     find the correction term for particle diffusion, if used.  it gets
c     added with the implicit advection terms later to avoid extra
c     ffts.
c-----------------------------------------------------------------------
      IF (continuity/='none'.AND.advect/='none'.AND.impladv.AND.
     $    nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        CALL generic_ptr_set(rb%qndiff_phi,tb%qndiff_phi,tb%tgeom,
     $                       inode,ndiff,dp,dp,0_i4)
        real_vec(1,:,:)=real_ve(1,:,:)*ndiff(1,:,:)*hdt*mtot
        real_vec(2,:,:)=real_ve(2,:,:)*ndiff(1,:,:)*hdt*mtot
        real_vec(3,:,:)=real_ve(3,:,:)*ndiff(1,:,:)*hdt*mtot
      ENDIF
c-----------------------------------------------------------------------
c     find rho*dve product for the d/dt term.
c-----------------------------------------------------------------------
      IF (continuity=='full') THEN
        real_ve(1,:,:)=real_ve(1,:,:)*real_ndptr(1,:,:)*mtot
        real_ve(2,:,:)=real_ve(2,:,:)*real_ndptr(1,:,:)*mtot
        real_ve(3,:,:)=real_ve(3,:,:)*real_ndptr(1,:,:)*mtot
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,force,real_ve,
     $               dealiase)
      ELSE
        DO imode=1,nmodes
          force(1,:,:,imode)=ve(1,:,:,imode)*n0*mtot
          force(2,:,:,imode)=ve(2,:,:,imode)*n0*mtot
          force(3,:,:,imode)=ve(3,:,:,imode)*n0*mtot
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluations for viscous dissipation:
c     first find and transform the nd*grad(V) tensor.
c-----------------------------------------------------------------------
      IF ((kin_visc>0.OR.iso_visc>0.OR.par_visc>0.OR.impladv).AND.
     $ (integrand_flag(1:4)=='visc'.OR.integrand_flag(1:3)=='all')) THEN

        CALL math_grad(nmodes,keff,3_i4,geom,ve,ver,vez,gvten,bigr)
        IF (continuity/='full') THEN
          DO imode=1,nmodes
            gvten(1,:,:,imode)=n0*gvten(1,:,:,imode)
            gvten(2,:,:,imode)=n0*gvten(2,:,:,imode)
            gvten(3,:,:,imode)=n0*gvten(3,:,:,imode)
            gvten(4,:,:,imode)=n0*gvten(4,:,:,imode)
            gvten(5,:,:,imode)=n0*gvten(5,:,:,imode)
            gvten(6,:,:,imode)=n0*gvten(6,:,:,imode)
            gvten(7,:,:,imode)=n0*gvten(7,:,:,imode)
            gvten(8,:,:,imode)=n0*gvten(8,:,:,imode)
            gvten(9,:,:,imode)=n0*gvten(9,:,:,imode)
          ENDDO
        ENDIF
        IF (continuity=='full'.OR.par_visc>0.OR.gyr_visc>0.OR.
     $      impladv) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,9_i4,
     $                 gvten,real_gvten,dealiase)
        ENDIF
        IF (continuity=='full') THEN
          DO jt=1,9
            real_gvten(jt,:,:)=real_gvten(jt,:,:)*real_ndptr(1,:,:)
          ENDDO
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,9_i4,
     $                 gvten,real_gvten,dealiase)
        ENDIF
c-----------------------------------------------------------------------
c       for nonlinear parallel viscosity, find the unique components of
c
c         [ B.grad(V).B - B**2*div(V) ] * (BB - I*B**2/3)
c
c       also divide by B**4 to make unit vectors.
c-----------------------------------------------------------------------
        IF (par_visc>0) THEN
          DO ip=1,nphi
            DO ix=1,mpseudo
              rb2=SUM(real_bptr(:,ix,ip)*real_bptr(:,ix,ip))+smallnum
              divtmpr=-third*rb2*(real_gvten(1,ix,ip)+
     $                real_gvten(5,ix,ip)+real_gvten(9,ix,ip))
              rscal=par_coef*(divtmpr+real_bptr(1,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(1,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(2,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(3,ix,ip))+
     $              real_bptr(2,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(4,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(5,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(6,ix,ip))+
     $              real_bptr(3,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(7,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(8,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(9,ix,ip)))/rb2**2
              rb2=-third*rb2
              real_pten(1,ix,ip)=rscal*
     $          (real_bptr(1,ix,ip)*real_bptr(1,ix,ip)+rb2)
              real_pten(2,ix,ip)=rscal*
     $           real_bptr(1,ix,ip)*real_bptr(2,ix,ip)
              real_pten(3,ix,ip)=rscal*
     $           real_bptr(1,ix,ip)*real_bptr(3,ix,ip)
              real_pten(4,ix,ip)=rscal*
     $          (real_bptr(2,ix,ip)*real_bptr(2,ix,ip)+rb2)
              real_pten(5,ix,ip)=rscal*
     $           real_bptr(2,ix,ip)*real_bptr(3,ix,ip)
              real_pten(6,ix,ip)=rscal*
     $          (real_bptr(3,ix,ip)*real_bptr(3,ix,ip)+rb2)
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         multiply by the temperature-dependent profile for Braginskii
c         parallel viscosity if desired.
c-----------------------------------------------------------------------
          IF (parvisc_model=='plltdep') THEN
            CALL generic_ptr_set(rb%qkappli_phi,tb%qkappli_phi,
     $                           tb%tgeom,inode,kappli,dp,dp,0_i4)
            DO ip=1,nphi
              DO ix=1,mpseudo
                rscal=kappli(1,ix,ip)/k_plli
                real_pten(1,ix,ip)=real_pten(1,ix,ip)*rscal
                real_pten(2,ix,ip)=real_pten(2,ix,ip)*rscal
                real_pten(3,ix,ip)=real_pten(3,ix,ip)*rscal
                real_pten(4,ix,ip)=real_pten(4,ix,ip)*rscal
                real_pten(5,ix,ip)=real_pten(5,ix,ip)*rscal
                real_pten(6,ix,ip)=real_pten(6,ix,ip)*rscal
              ENDDO
            ENDDO
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       for nonlinear gyroviscosity, find the unique components of
c
c         (n_e*m_i*T_i/Z**2*B**4)*{ [BxW.(B**2+3BB) ] + transpose[] }
c
c-----------------------------------------------------------------------
        IF (gyr_visc>0) THEN
          vcoef=
     $      -0.125_r8*gyr_visc*dt*ms(2)*kboltz/(zeff**2*elementary_q)
          CALL generic_ptr_set(rb%qti_tot,tb%qti_tot,tb%tgeom,
     $                         inode,real_titot,dp,dp,0_i4)
          IF (par_visc<=0) real_pten=0._r8
          DO ip=1,nphi
            DO ix=1,mpseudo
              rb2=SUM(real_bptr(:,ix,ip)**2)+smallnum
              gcoef=vcoef*real_titot(1,ix,ip)/rb2**2
              vtmpr(1:3,1)=real_gvten(1:3,ix,ip)+real_gvten(1:7:3,ix,ip)
              vtmpr(1:3,2)=real_gvten(4:6,ix,ip)+real_gvten(2:8:3,ix,ip)
              vtmpr(1:3,3)=real_gvten(7:9,ix,ip)+real_gvten(3:9:3,ix,ip)
              btmp(1:3,1)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(1  ,ix,ip)
              btmp(1:3,2)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(2  ,ix,ip)
              btmp(1:3,3)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(3  ,ix,ip)
              btmp(1,1)=rb2+btmp(1,1)
              btmp(2,2)=rb2+btmp(2,2)
              btmp(3,3)=rb2+btmp(3,3)
              wdotr(1,1)=SUM(vtmpr(1,:)*btmp(:,1))
              wdotr(2,1)=SUM(vtmpr(2,:)*btmp(:,1))
              wdotr(3,1)=SUM(vtmpr(3,:)*btmp(:,1))
              wdotr(1,2)=SUM(vtmpr(1,:)*btmp(:,2))
              wdotr(2,2)=SUM(vtmpr(2,:)*btmp(:,2))
              wdotr(3,2)=SUM(vtmpr(3,:)*btmp(:,2))
              wdotr(1,3)=SUM(vtmpr(1,:)*btmp(:,3))
              wdotr(2,3)=SUM(vtmpr(2,:)*btmp(:,3))
              wdotr(3,3)=SUM(vtmpr(3,:)*btmp(:,3))
              bc_wdotr(1,:)=real_bptr(2,ix,ip)*wdotr(3,:)-
     $                      real_bptr(3,ix,ip)*wdotr(2,:)
              bc_wdotr(2,:)=real_bptr(3,ix,ip)*wdotr(1,:)-
     $                      real_bptr(1,ix,ip)*wdotr(3,:)
              bc_wdotr(3,:)=real_bptr(1,ix,ip)*wdotr(2,:)-
     $                      real_bptr(2,ix,ip)*wdotr(1,:)
              real_pten(1:3,ix,ip)=real_pten(1:3,ix,ip)+
     $                           gcoef*(bc_wdotr(:,1)+bc_wdotr(1,:))
              real_pten(4  ,ix,ip)=real_pten(4  ,ix,ip)+
     $                           gcoef*(bc_wdotr(2,2)+bc_wdotr(2,2))
              real_pten(5:6,ix,ip)=real_pten(5:6,ix,ip)+
     $                           gcoef*(bc_wdotr(3,2:3)+bc_wdotr(2:3,3))
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       transform the anisotropic stresses.
c-----------------------------------------------------------------------
        IF (par_visc>0.OR.gyr_visc>0) THEN
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,6_i4,
     $                 pten,real_pten,dealiase)
        ENDIF
c-----------------------------------------------------------------------
c       complete -Pi(rho,V) and save it in gvten.
c-----------------------------------------------------------------------
        DO im=1,nmodes
          DO iy=1,ncy
            DO ix=1,ncx
              dloc=diff(ix,iy)
              vtmp(1:3,1)=gvten(1:3,ix,iy,im)
              vtmp(1:3,2)=gvten(4:6,ix,iy,im)
              vtmp(1:3,3)=gvten(7:9,ix,iy,im)
              div=-twoth*(vtmp(1,1)+vtmp(2,2)+vtmp(3,3))
              gvten(1,ix,iy,im)=dloc*((kin_coef+2*iso_coef)*vtmp(1,1)+
     $                          iso_coef*div)
              gvten(2,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(2,1)+
     $                          iso_coef*vtmp(1,2))
              gvten(3,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(3,1)+
     $                          iso_coef*vtmp(1,3))
              gvten(4,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(1,2)+
     $                          iso_coef*vtmp(2,1))
              gvten(5,ix,iy,im)=dloc*((kin_coef+2*iso_coef)*vtmp(2,2)+
     $                          iso_coef*div)
              gvten(6,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(3,2)+
     $                          iso_coef*vtmp(2,3))
              gvten(7,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(1,3)+
     $                          iso_coef*vtmp(3,1))
              gvten(8,ix,iy,im)=dloc*((kin_coef+  iso_coef)*vtmp(2,3)+
     $                          iso_coef*vtmp(3,2))
              gvten(9,ix,iy,im)=dloc*((kin_coef+2*iso_coef)*vtmp(3,3)+
     $                          iso_coef*div)
            ENDDO
          ENDDO
        ENDDO
        IF (par_visc>0.OR.gyr_visc>0) THEN
          gvten(1:3,:,:,:)=gvten(1:3,:,:,:)+pten(1:3,:,:,:)
          gvten(4,:,:,:)=gvten(4,:,:,:)+pten(2,:,:,:)
          gvten(5,:,:,:)=gvten(5,:,:,:)+pten(4,:,:,:)
          gvten(6,:,:,:)=gvten(6,:,:,:)+pten(5,:,:,:)
          gvten(7,:,:,:)=gvten(7,:,:,:)+pten(3,:,:,:)
          gvten(8,:,:,:)=gvten(8,:,:,:)+pten(5,:,:,:)
          gvten(9,:,:,:)=gvten(9,:,:,:)+pten(6,:,:,:)
        ENDIF
c-----------------------------------------------------------------------
c       terms for implicit advection:
c
c          0.5*dt*rho*( dV.grad(V_old) + V_old.grad(dV) )
c
c       where dV is the iterate from work4.
c
c       note that V_old already includes any equilibrium flow.  Also,
c       real_gvten still holds n*grad(dV) in real space, force still 
c       holds mtot*n*dV by Fourier component, and real_ve holds 
c       mtot*n*dV in real space if continuity is full.
c-----------------------------------------------------------------------
        IF (impladv.AND.(advect=='V only'.OR.advect=='all')) THEN
          IF (.NOT.continuity=='full') CALL fft_nim('inverse',ncx*ncy,
     $      mpseudo,lphi,3_i4,force,real_ve,dealiase)
          DO ip=1,nphi
            DO ix=1,mpseudo
               tmp1=hdt*real_ve(:,ix,ip)
               tmp2=hdt*mtot*real_oldv(:,ix,ip)
               real_ve(1,ix,ip)=real_ve(1,ix,ip)+
     $                   tmp1(1)*real_oldgv(1,ix,ip)+
     $                   tmp1(2)*real_oldgv(2,ix,ip)+
     $                   tmp1(3)*real_oldgv(3,ix,ip)+
     $                   tmp2(1)*real_gvten(1,ix,ip)+
     $                   tmp2(2)*real_gvten(2,ix,ip)+
     $                   tmp2(3)*real_gvten(3,ix,ip)
               real_ve(2,ix,ip)=real_ve(2,ix,ip)+
     $                   tmp1(1)*real_oldgv(4,ix,ip)+
     $                   tmp1(2)*real_oldgv(5,ix,ip)+
     $                   tmp1(3)*real_oldgv(6,ix,ip)+
     $                   tmp2(1)*real_gvten(4,ix,ip)+
     $                   tmp2(2)*real_gvten(5,ix,ip)+
     $                   tmp2(3)*real_gvten(6,ix,ip)
               real_ve(3,ix,ip)=real_ve(3,ix,ip)+
     $                   tmp1(1)*real_oldgv(7,ix,ip)+
     $                   tmp1(2)*real_oldgv(8,ix,ip)+
     $                   tmp1(3)*real_oldgv(9,ix,ip)+
     $                   tmp2(1)*real_gvten(7,ix,ip)+
     $                   tmp2(2)*real_gvten(8,ix,ip)+
     $                   tmp2(3)*real_gvten(9,ix,ip)
            ENDDO
          ENDDO
          IF (continuity/='none'.AND.nd_correrr.AND.
     $        (nd_diff>0.OR.nd_hypd>0)) real_ve=real_ve+real_vec
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,
     $                 force,real_ve,dealiase)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     evaluate fields for the semi-implicit operator.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,beq,beqr,beqz,1_i4)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,j0,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,inode,
     $                     b0,b0r,b0z,1_i4)
      ja(1,:,:)=j0(1,:,:)+ b0z(3,:,:)/mu0
      ja(2,:,:)=j0(2,:,:)- b0r(3,:,:)/mu0
      ja(3,:,:)=j0(3,:,:)+(b0r(2,:,:)-b0z(1,:,:))/mu0
      IF (geom=='tor') ja(2,:,:)=ja(2,:,:)-b0(3,:,:)/(mu0*bigr)
      be=b0+beq
      ber=b0r+beqr
      bez=b0z+beqz
      CALL generic_ptr_set(rb%qsi_nl_pres,tb%qsi_nl_pres,tb%tgeom,
     $                     inode,nl_pres,dp,dp,0_i4)
      bmagsq=SUM(be*be,1)
      IF (beta>0) THEN
        CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                       inode,peq,peqr,peqz,1_i4)
        CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                       inode,p0,p0r,p0z,1_i4)
        pres=p0+peq
        presr=p0r+peqr
        presz=p0z+peqz
      ELSE
        pres=0._r8
        presr=0._r8
        presz=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     find the spatial functions used as coefficients for the isotropic
c     and the anisotropic parts of the  operator.
c-----------------------------------------------------------------------
      ani=0.25_r8*si_fac_mhd* (1._r8-mhd_si_iso)*dt2/mu0
      cj0=0.125_r8*si_fac_j0*dt2
      IF (beta>0) THEN
        iso=0.25_r8*dt2*( mhd_si_iso*
     $    (si_fac_mhd*bmagsq/mu0+si_fac_pres*gamma*pres(1,:,:))+
     $      si_fac_nl*(nl_pres(1,:,:)/mu0+gamma*nl_pres(2,:,:)) )
        grd=0.25_r8*si_fac_pres*(1._r8-mhd_si_iso)*dt2*
     $              gamma*pres(1,:,:)
      ELSE
        iso=0.25_r8*dt2*( mhd_si_iso*si_fac_mhd*bmagsq+
     $                    si_fac_nl*nl_pres(1,:,:) )/mu0
        grd=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     the anisotropic operator has a term of the form
c     curl(alpha_iXB0).curl(alpha_jXB0) after the integration by parts,
c     and each curl is broken into B0.grad(alpha)-alpha.grad(B0)
c     -B0(div(alpha)).  B0 is in cylindrical components.
c     the first two indices of the temporary crl matrix are defined as
c     ( curl vector index, operand vector index ), where the vector
c     indices runs from 1 to 3 as (r,z,phi), and the last index is
c     the basis function index.
c
c     dbe is the ideal change in b, curl(veXB0).
c
c     this part of the operation is diagonal in Fourier component.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     this part of crl is k-independent.
c-PRE and could be saved for each step.
c-----------------------------------------------------------------------
      DO iv=1,nvc
        rcrl(1,1,:,:,iv)=
     $    dalpdz(:,:,iv)*  be(2,:,:)
     $    -alpha(:,:,iv)*(be(1,:,:)*g/bigr+ber(1,:,:))
        rcrl(1,2,:,:,iv)=
     $    -(dalpdz(:,:,iv)*be(1,:,:)+alpha(:,:,iv)*bez(1,:,:))
        rcrl(1,3,:,:,iv)=0._r8
        rcrl(2,1,:,:,iv)=
     $    -(dalpdr(:,:,iv)* be(2,:,:)
     $      +alpha(:,:,iv)*(be(2,:,:)*g/bigr+ber(2,:,:)))
        rcrl(2,2,:,:,iv)=
     $    dalpdr(:,:,iv)*be(1,:,:)-alpha(:,:,iv)*bez(2,:,:)
        rcrl(2,3,:,:,iv)=0._r8
        rcrl(3,1,:,:,iv)=
     $    -(dalpdr(:,:,iv)*be(3,:,:)+alpha(:,:,iv)*ber(3,:,:))
        rcrl(3,2,:,:,iv)=
     $    -(dalpdz(:,:,iv)*be(3,:,:)+alpha(:,:,iv)*bez(3,:,:))
        rcrl(3,3,:,:,iv)=
     $    dalpdr(:,:,iv)*be(1,:,:)+dalpdz(:,:,iv)*be(2,:,:)
     $    -alpha(:,:,iv)*be(1,:,:)*g/bigr
      ENDDO
c-----------------------------------------------------------------------
c     if auxiliary fields are used, find their contribution to the
c     velocity bases and the lhs of its own evolution equation.  the
c     equations for div(V) are expressed directly in weak form:
c
c       int[rho*A^*.delta(V)] = dt*int[A^*.forces]
c               - int[ sqrt(ddivv*nu*dt*f)*auxv*div(A^*) ]
c       int[(ups^*)*auxv]-int[sqrt(ddivv*nu*dt*f)*(ups^*)*div(delta(V))]
c         = int[ sqrt(ddivv*nu*dt/f)*(ups^*)*div(V_old)]
c
c     where f is the implicit centering, nu is a numerical viscosity
c     coefficient, nu=dt*cma**2, cma is sqrt(rho) times the
c     magneto-acoustic speed, and ups is the scalar test field for the
c     auxiliary field auxv.  note that the real pres array becomes
c     sqrt(gamma*P+B^2/mu0) and that mwork4 holds the div(delta-V) and
c     B.curl(delta-V) projections.
c
c     the auxiliary field to stabilize vorticity (auxw) is expressed
c     directly in weak form as
c
c       int[rho*A^*.delta(V)] = dt*int[A^*.forces]
c               - int[ dt*sqrt(dpvrt*f/mu0)*auxw*B.curl(A^*) ]
c       int[(lam^*)*auxw]
c        - dt*sqrt(dpvrt*f/mu0)*int[(lam^*)*B.curl(delta(V))]
c        = dt*sqrt(dpvrt/(f*mu0))*int[(lam^*)*B.curl(V_old)]
c
c     where A is the test vector for V, lam is the test function for
c     the discontinuous scalar field auxw, where f is the centering
c     coefficient.  here, the effective viscosity is dt*v_Alven**2.
c
c     divv is again used for temporary space, and force becomes
c     dt*cpvrt*curl(V)/sqrt(mu0).
c-----------------------------------------------------------------------
      IF (poly_divv>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',almod,dalmr,
     $                          dalmz,0_i4,poly_divv,
     $                          polydmin=poly_divv_min,
     $                          polydmax=poly_divv_max)
        nvm=SIZE(almod,3)
        CALL generic_all_eval(rb%mwork4,tb%mwork4,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,aux,dcp,dcp,0_i4)

        IF (beta>0) THEN
          sqptot=SQRT(MAX(0._r8,gamma*pres(1,:,:))+SUM(be**2,1)/mu0)
        ELSE
          sqptot=SQRT(SUM(be**2,1)/mu0)
        ENDIF
        dtdv=SQRT(ddivv*fdivv)*dt        !  for diffusive correction
        dtpv=dt*SQRT(dpvrt*fpvrt/(mu0))  !  for diffusive correction
      ELSE
        dtdv=0._r8
        dtpv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     start mode loop and find dbe.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        dbe(1,:,:)=
     $      vez(1,:,:,imode)*  be(2,:,:)
     $      -ve(1,:,:,imode)*((be(1,:,:)*g-(0,1)*keff(imode)*be(3,:,:))
     $                      /bigr+ber(1,:,:))
     $    -(vez(2,:,:,imode)*be(1,:,:)+ve(2,:,:,imode)*bez(1,:,:))
     $      -ve(3,:,:,imode)*(0,1)*keff(imode)*be(1,:,:)/bigr
        dbe(2,:,:)=
     $    -(ver(1,:,:,imode)* be(2,:,:)
     $      +ve(1,:,:,imode)*(be(2,:,:)*g/bigr+ber(2,:,:)))
     $    + ver(2,:,:,imode)*be(1,:,:)
     $    +  ve(2,:,:,imode)*((0,1)*keff(imode)*be(3,:,:)/bigr
     $                       -bez(2,:,:))
     $    -  ve(3,:,:,imode)*(0,1)*keff(imode)*be(2,:,:)/bigr
        dbe(3,:,:)=
     $    -(ver(1,:,:,imode)*be(3,:,:)+ ve(1,:,:,imode)*ber(3,:,:)
     $     +vez(2,:,:,imode)*be(3,:,:)+ ve(2,:,:,imode)*bez(3,:,:))
     $    + ver(3,:,:,imode)*be(1,:,:)+vez(3,:,:,imode)*be(2,:,:)
     $    -  ve(3,:,:,imode)*be(1,:,:)*g/bigr
        CALL math_cart_cross(ja_cross_dbe,ja,dbe,1._r8)
        CALL math_cart_cross(ve_cross_ja,ve(:,:,:,imode),ja,1._r8)
c-----------------------------------------------------------------------
c       isotropic terms that multiply alpha only:
c-----------------------------------------------------------------------
        IF (geom=='tor') THEN
          forck(1,:,:)=force(1,:,:,imode)+
     $      iso*((1+k2ef(imode))*ve(1,:,:,imode)+
     $        2*(0,1)*keff(imode)*ve(3,:,:,imode))/bigr2
          forck(2,:,:)=force(2,:,:,imode)+
     $      iso*k2ef(imode)*ve(2,:,:,imode)/bigr2
          forck(3,:,:)=force(3,:,:,imode)+
     $      iso*((1+k2ef(imode))*ve(3,:,:,imode)-
     $        2*(0,1)*keff(imode)*ve(1,:,:,imode))/bigr2
          divv=ver(1,:,:,imode)+ve(1,:,:,imode)/bigr+vez(2,:,:,imode)
     $        +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr
        ELSE
          forck(1,:,:)=force(1,:,:,imode)+
     $      iso*k2ef(imode)*ve(1,:,:,imode)/bigr2
          forck(2,:,:)=force(2,:,:,imode)+
     $      iso*k2ef(imode)*ve(2,:,:,imode)/bigr2
          forck(3,:,:)=force(3,:,:,imode)+
     $      iso*k2ef(imode)*ve(3,:,:,imode)/bigr2
          divv=ver(1,:,:,imode)+vez(2,:,:,imode)
     $        +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr
        ENDIF
        IF (use_gv) THEN
          IF (geom=='tor') THEN
            forck(1,:,:)=forck(1,:,:)
     $        +(                  gvten(9,:,:,imode)
     $         -(0,1)*keff(imode)*gvten(3,:,:,imode))/bigr
            forck(2,:,:)=forck(2,:,:)
     $         -(0,1)*keff(imode)*gvten(6,:,:,imode) /bigr
            forck(3,:,:)=forck(3,:,:)
     $        -(                  gvten(3,:,:,imode)
     $         +(0,1)*keff(imode)*gvten(9,:,:,imode))/bigr
          ELSE
            forck(1,:,:)=forck(1,:,:)
     $         -(0,1)*keff(imode)*gvten(3,:,:,imode)/bigr
            forck(2,:,:)=forck(2,:,:)
     $         -(0,1)*keff(imode)*gvten(6,:,:,imode)/bigr
            forck(3,:,:)=forck(3,:,:)
     $         -(0,1)*keff(imode)*gvten(9,:,:,imode)/bigr
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       compute isotropic, anisotropic, and pressure contributions
c       together for effeciency.
c
c       the "ani" contributions are the dot product of the curls
c       from the different vertices around each cell and the
c
c       symmetrized contributions from the J0Xcurl(vxB0) term are
c       the "cj0" terms without pressure gradient coefficients.
c
c       grad(P0) terms include the div(alpha)*div() part of the operator
c       and the symmetrized div(alpha)*alpha.grad(P0) part.
c
c       the mass density term is with the sink term.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          ve_gpr=ve(1,:,:,imode)*presr(1,:,:)
     $          +ve(2,:,:,imode)*presz(1,:,:)
        ELSE
          ve_gpr=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       find the contributions of P-total times div(dV) and B.curl(dV)
c       from the auxiliary discontinuous fields.
c-----------------------------------------------------------------------
        IF (poly_divv>0) THEN
          sqpaux=dtdv*sqptot*aux(1,:,:,imode)
          baux(1,:,:)=dtpv*be(1,:,:)*aux(2,:,:,imode)
          baux(2,:,:)=dtpv*be(2,:,:)*aux(2,:,:,imode)
          baux(3,:,:)=dtpv*be(3,:,:)*aux(2,:,:,imode)
          bvort=(be(1,:,:)*(vez(3,:,:,imode)
     $                     -(0,1)*keff(imode)*ve(2,:,:,imode)/bigr)+
     $           be(2,:,:)*(((0,1)*keff(imode)*ve(1,:,:,imode)-
     $                      g*ve(3,:,:,imode))/bigr-ver(3,:,:,imode))+
     $           be(3,:,:)*(ver(2,:,:,imode)-vez(1,:,:,imode)))*dtpv
        ELSE
          sqpaux=0._r8
          baux=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       accumulate all contributions, including viscous terms
c       -conjg(grad(alpha))^T : Pi(rho,V)
c
c       the k-dependent part of crl is completed here.  note that crl is
c       explicitly conjugated here since it is only used in conjugated
c       form.
c-----------------------------------------------------------------------
        IF (use_gv) THEN
          DO iv=1,nvc
            DO iy=1,ncy
              int(:,iy,iv,imode)=0._r8
              DO ix=1,ncx
                rr=bigr(ix,iy)
                crl(1,1)=rcrl(1,1,ix,iy,iv)
     $            -alpha(ix,iy,iv)*(0,1)*keff(imode)*be(3,ix,iy)/rr
                crl(2,1)=rcrl(2,1,ix,iy,iv)
                crl(3,1)=rcrl(3,1,ix,iy,iv)
                crl(1,2)=rcrl(1,2,ix,iy,iv)
                crl(2,2)=rcrl(2,2,ix,iy,iv)
     $            -alpha(ix,iy,iv)*(0,1)*keff(imode)*be(3,ix,iy)/rr
                crl(3,2)=rcrl(3,2,ix,iy,iv)
                crl(1,3)=rcrl(1,3,ix,iy,iv)
     $            +alpha(ix,iy,iv)*(0,1)*keff(imode)*be(1,ix,iy)/rr
                crl(2,3)=rcrl(2,3,ix,iy,iv)
     $            +alpha(ix,iy,iv)*(0,1)*keff(imode)*be(2,ix,iy)/rr
                crl(3,3)=rcrl(3,3,ix,iy,iv)
                div=grd(ix,iy)*divv(ix,iy)+sqpaux(ix,iy)

                int(1,iy,iv,imode)=int(1,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(1,ix,iy)+
     $                cj0*(divv(ix,iy)*presr(1,ix,iy)-
     $                          ja_cross_dbe(1,ix,iy))-
     $                baux(2,ix,iy)*(0,1)*keff(imode)/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,1))
     $            +dalpdrc(ix,iy,iv)*(div+cj0*ve_gpr(ix,iy))
     $            +ani*SUM(crl(:,1)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(1,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(1,ix,iy,imode))
     $            +dalpdr(ix,iy,iv)*gvten(1,ix,iy,imode)
     $            +dalpdz(ix,iy,iv)*(gvten(2,ix,iy,imode)-baux(3,ix,iy))
 
                int(2,iy,iv,imode)=int(2,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(2,ix,iy)+
     $                cj0*(divv(ix,iy)*presz(1,ix,iy)-
     $                          ja_cross_dbe(2,ix,iy))+
     $                baux(1,ix,iy)*(0,1)*keff(imode)/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,2))
     $            +dalpdz (ix,iy,iv)*(div+cj0*ve_gpr(ix,iy))
     $            +ani*SUM(crl(:,2)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(2,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(2,ix,iy,imode))
     $            +dalpdr(ix,iy,iv)*(gvten(4,ix,iy,imode)+baux(3,ix,iy))
     $            +dalpdz(ix,iy,iv)*gvten(5,ix,iy,imode)

                int(3,iy,iv,imode)=int(3,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(3,ix,iy)
     $               -cj0*ja_cross_dbe(3,ix,iy)
     $               -(0,1)*keff(imode)*(div+cj0*ve_gpr(ix,iy))/rr
     $               -baux(2,ix,iy)*g/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,3))
     $            +ani*SUM(crl(:,3)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(3,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(3,ix,iy,imode))
     $            +dalpdr(ix,iy,iv)*(gvten(7,ix,iy,imode)-baux(2,ix,iy))
     $            +dalpdz(ix,iy,iv)*(gvten(8,ix,iy,imode)+baux(1,ix,iy))
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c       same without viscous terms.
c-----------------------------------------------------------------------
        ELSE
          DO iv=1,nvc
            DO iy=1,ncy
              int(:,iy,iv,imode)=0._r8
              DO ix=1,ncx
                rr=bigr(ix,iy)
                crl(1,1)=rcrl(1,1,ix,iy,iv)
     $            -alpha(ix,iy,iv)*(0,1)*keff(imode)*be(3,ix,iy)/rr
                crl(2,1)=rcrl(2,1,ix,iy,iv)
                crl(3,1)=rcrl(3,1,ix,iy,iv)
                crl(1,2)=rcrl(1,2,ix,iy,iv)
                crl(2,2)=rcrl(2,2,ix,iy,iv)
     $            -alpha(ix,iy,iv)*(0,1)*keff(imode)*be(3,ix,iy)/rr
                crl(3,2)=rcrl(3,2,ix,iy,iv)
                crl(1,3)=rcrl(1,3,ix,iy,iv)
     $            +alpha(ix,iy,iv)*(0,1)*keff(imode)*be(1,ix,iy)/rr
                crl(2,3)=rcrl(2,3,ix,iy,iv)
     $            +alpha(ix,iy,iv)*(0,1)*keff(imode)*be(2,ix,iy)/rr
                crl(3,3)=rcrl(3,3,ix,iy,iv)
                div=grd(ix,iy)*divv(ix,iy)+sqpaux(ix,iy)

                int(1,iy,iv,imode)=int(1,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(1,ix,iy)+
     $                cj0*(divv(ix,iy)*presr(1,ix,iy)-
     $                          ja_cross_dbe(1,ix,iy))-
     $                baux(2,ix,iy)*(0,1)*keff(imode)/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,1))
     $            +dalpdrc(ix,iy,iv)*(div+cj0*ve_gpr(ix,iy))
     $            +ani*SUM(crl(:,1)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(1,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(1,ix,iy,imode))
     $            -dalpdz(ix,iy,iv)*baux(3,ix,iy)
 
                int(2,iy,iv,imode)=int(2,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(2,ix,iy)+
     $                cj0*(divv(ix,iy)*presz(1,ix,iy)-
     $                          ja_cross_dbe(2,ix,iy))+
     $                baux(1,ix,iy)*(0,1)*keff(imode)/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,2))
     $            +dalpdz (ix,iy,iv)*(div+cj0*ve_gpr(ix,iy))
     $            +ani*SUM(crl(:,2)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(2,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(2,ix,iy,imode))
     $            +dalpdr(ix,iy,iv)*baux(3,ix,iy)

                int(3,iy,iv,imode)=int(3,iy,iv,imode)+
     $            alpha(ix,iy,iv)*(forck(3,ix,iy)
     $               -cj0*ja_cross_dbe(3,ix,iy)
     $               -(0,1)*keff(imode)*(div+cj0*ve_gpr(ix,iy))/rr
     $               -baux(2,ix,iy)*g/rr)
     $            -cj0*SUM(ve_cross_ja(:,ix,iy)*crl(:,3))
     $            +ani*SUM(crl(:,3)*dbe(:,ix,iy))
     $            +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(3,ix,iy,imode)
     $                        +dalpdz(ix,iy,iv)*vez(3,ix,iy,imode))
     $            -dalpdr(ix,iy,iv)*baux(2,ix,iy)
     $            +dalpdz(ix,iy,iv)*baux(1,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       complete the integrand for the discontinuous scalar equations.
c-----------------------------------------------------------------------
        IF (poly_divv>0) THEN
          divv=divv*sqptot*dtdv
          DO iv=nvc+1,nvc+nvm  !  auxiliary equation
            int(1,:,iv,imode)=
     $        SUM(almod(:,:,iv-nvc)*(aux(1,:,:,imode)-divv),1)
            int(2,:,iv,imode)=
     $        SUM(almod(:,:,iv-nvc)*(aux(2,:,:,imode)-bvort),1)
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE v_aniso_dot
c-----------------------------------------------------------------------
c     subprogram 5. v_3dsi_dot.
c     find the action of the (mass density*I + viscous + semi-implicit)
c     operator dotted into a vector (work4) for CG or GMRES iterations.
c     this version is for computations with the 3D semi-implicit
c     operator, with or without implicit advection.  the matrix couples
c     Fourier components.
c
c     the semi-implicit terms are computed with 3D magnetic field,
c     current density, and pressure in this version.  note that after
c     integration by parts, the volumetric terms can be arranged into
c
c     A^*.[df-grad(B).db] + div(A^*)*dp + grad(A^*)^T:[B db], where
c
c     df=0.5*[div(dv)*grad(P)-JXcurl(dvXB)],
c     db=curl(dvXB)/mu0-0.5*dvXJ,
c     dp=gamma*P*div(dv)+0.5*dv.grad(P)-B.db, and
c     A^* is the complex conjugate of a test vector.
c
c     B, grad(B), J, P, and grad(P) have been stored as functions of
c     the periodic angle including steady-state parts.  the terms
c     that are multiplied by grad(A^*) are assembled in the stress
c     tensor before viscous effects are computed; this order of
c     operations differs from that in v_aniso_dot.
c-----------------------------------------------------------------------
      SUBROUTINE v_3dsi_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          almod,dalmr,dalmz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          real_ndptr,real_bptr,
     $          real_jptr,real_pptr,real_grdb,real_grdp,
     $          real_oldv,real_oldgv,beq,b0,nd_eq,nd_n0,nl_pres,peq,
     $          p0,kappli,dp,real_titot,ds2
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: pres
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: n0,bigr2,
     $          iso,bmagsq,sqptot
      REAL(r8), DIMENSION(SIZE(bigr,1)*SIZE(bigr,2)) :: n0_1d,df_1d
      REAL(r8), DIMENSION(9,mpseudo,nphi) :: real_gvten,real_pten
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_vtmp
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal
      REAL(r8), DIMENSION(3,3) :: vtmpr,btmp,wdotr,bc_wdotr
      REAL(r8), DIMENSION(3) :: tmp1,tmp2,dbeff,dfeff

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: dcp
      COMPLEX(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: aux
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ve,
     $             ver,vez,force
      COMPLEX(r8), DIMENSION(9,SIZE(bigr,1),SIZE(bigr,2),nmodes) ::
     $             gvten,pten
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: forck
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: divv,bvort
      COMPLEX(r8), DIMENSION(3) :: baux
      COMPLEX(r8) :: div,sqpaux

      INTEGER(i4) :: ncx,ncy,iv,imode,jt,ix,iy,ip,im,i1,nvc,nvm
      REAL(r8) :: kin_coef,iso_coef,par_coef,gcoef,vcoef,dloc
      REAL(r8) :: rb2,divtmpr,rscal,hdt,g,ani,cj0,dt2,cp0,rr,dtdv,dtpv
      REAL(r8) :: ten12,ten13,ten23,dpeff
      REAL(r8), PARAMETER :: third=1._r8/3._r8,twoth=2._r8/3._r8
c-----------------------------------------------------------------------
c     convenience
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      hdt=0.5_r8*dt
      bigr2=bigr*bigr
      dt2=dt*dt
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
c-----------------------------------------------------------------------
c     evaluate the current direction vector for velocity change and
c     the number density.
c
c     also multiply the diffusivity shape by n0.
c-----------------------------------------------------------------------
      CALL generic_all_eval(rb%work4,tb%work4,rb%dxdr,rb%dydr,
     $                      rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,tb%ng,
     $                      ve,ver,vez,1_i4)
      CALL math_grad(nmodes,keff,3_i4,geom,ve,ver,vez,gvten,bigr)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,9_i4,gvten,
     $             real_gvten,dealiase)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $             dealiase)

      IF (continuity=='full') THEN
        CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)
        real_scal=real_ndptr
      ELSE
        CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                       inode,nd_eq,dp,dp,0_i4)
        real_scal=1._r8
      ENDIF
      IF (continuity=='n=0 only') THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
        n0=nd_eq(1,:,:)+nd_n0(1,:,:)
        n0_1d=RESHAPE(n0,(/ncx*ncy/))
      ELSE IF (continuity=='full') THEN
        n0=1._r8
      ELSE
        n0=nd_eq(1,:,:)
        n0_1d=RESHAPE(n0,(/ncx*ncy/))
      ENDIF

      IF (impladv) THEN
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,real_oldv,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                       inode,real_oldgv,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     collect information for diffusivity shapes and coefficients.
c-----------------------------------------------------------------------
      IF ((ds_use=='kin_visc'.OR.ds_use=='both').AND.
     $    (kin_visc>0.OR.iso_visc>0)) THEN
        CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                       inode,ds2,dp,dp,0_i4)
        df_1d=RESHAPE(MAX(ds2(1,:,:),0._r8)*n0,(/ncx*ncy/))
      ELSE
        IF (continuity=='full') THEN
          df_1d=1._r8
        ELSE
          df_1d=n0_1d
        ENDIF
      ENDIF

      kin_coef=dt*fvsc*mtot*kin_visc
      iso_coef=dt*fvsc*mtot*iso_visc
      par_coef=3._r8*dt*fvsc*mtot*par_visc
c-----------------------------------------------------------------------
c     set pointers to fields used in the semi-implicit operator.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                     inode,real_bptr,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qgrdb,tb%qgrdb,tb%tgeom,
     $                     inode,real_grdb,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qja_tot,tb%qja_tot,tb%tgeom,
     $                     inode,real_jptr,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,beq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,inode,
     $                     b0,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qsi_nl_pres,tb%qsi_nl_pres,tb%tgeom,
     $                     inode,nl_pres,dp,dp,0_i4)
      be=b0+beq
      bmagsq=SUM(be**2,1)
      IF (beta>0) THEN
        CALL generic_ptr_set(rb%qpr_tot,tb%qpr_tot,tb%tgeom,
     $                       inode,real_pptr,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qgrdp,tb%qgrdp,tb%tgeom,
     $                       inode,real_grdp,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                       inode,peq,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                       inode,p0,dp,dp,0_i4)
        pres=p0+peq
      ELSE
c-PRE always allocate and zero-out in nimrod_init.
        ALLOCATE(real_pptr(1,mpseudo,nphi))
        ALLOCATE(real_grdp(3,mpseudo,nphi))
        real_pptr=0._r8
        real_grdp=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     if auxiliary fields are used, they have contributions to the
c     velocity bases and to the lhs of its own evolution equation.  the
c     equations for div(V) contributions are expressed directly in weak
c     form:
c
c       int[rho*A^*.delta(V)] = dt*int[A^*.forces]
c               - int[ sqrt(ddivv*nu*dt*f)*auxv*div(A^*) ]
c       int[(ups^*)*auxv]-int[sqrt(ddivv*nu*dt*f)*(ups^*)*div(delta(V))]
c         = int[ sqrt(ddivv*nu*dt/f)*(ups^*)*div(V_old)]
c
c     where f is the implicit centering, nu is a numerical viscosity
c     coefficient, nu=dt*cma**2, cma is sqrt(rho) times the
c     magneto-acoustic speed, and ups is the scalar test field for the
c     auxiliary field auxv.  note that the real pres array becomes
c     sqrt(gamma*P+B^2/mu0) and that mwork4 holds the div(delta-V) and
c     B.curl(delta-V) projections.
c
c     the auxiliary field to stabilize vorticity (auxw) is expressed
c     directly in weak form as
c
c       int[rho*A^*.delta(V)] = dt*int[A^*.forces]
c               - int[ dt*sqrt(dpvrt*f/mu0)*auxw*B.curl(A^*) ]
c       int[(lam^*)*auxw]
c        - dt*sqrt(dpvrt*f/mu0)*int[(lam^*)*B.curl(delta(V))]
c        = dt*sqrt(dpvrt/(f*mu0))*int[(lam^*)*B.curl(V_old)]
c
c     where A is the test vector for V, lam is the test function for
c     the discontinuous scalar field auxw, where f is the centering
c     coefficient.  here, the effective viscosity is dt*v_Alven**2.
c
c     at this point, just evaluate the auxiliary fields and coefficients
c     used in the equations.

cc     divv is again used for temporary space, and force becomes
cc     dt*cpvrt*curl(V)/sqrt(mu0).
c-----------------------------------------------------------------------
      IF (poly_divv>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',almod,dalmr,
     $                          dalmz,0_i4,poly_divv,
     $                          polydmin=poly_divv_min,
     $                          polydmax=poly_divv_max)
        nvm=SIZE(almod,3)
        CALL generic_all_eval(rb%mwork4,tb%mwork4,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,aux,dcp,dcp,0_i4)

        IF (beta>0) THEN
          sqptot=SQRT(gamma*pres(1,:,:)+SUM(be**2,1)/mu0)
        ELSE
          sqptot=SQRT(SUM(be**2,1)/mu0)
        ENDIF
        dtdv=SQRT(ddivv*fdivv)*dt        !  for diffusive correction
        dtpv=dt*SQRT(dpvrt*fpvrt/mu0)  !  for diffusive correction
      ENDIF
c-----------------------------------------------------------------------
c     find the spatial functions used as coefficients for the isotropic
c     and the anisotropic parts of the  operator.
c-----------------------------------------------------------------------
      ani=0.25_r8*si_fac_mhd* (1._r8-mhd_si_iso)*dt2/mu0
      cp0=0.25_r8*si_fac_pres*(1._r8-mhd_si_iso)*dt2*gamma
      cj0=0.125_r8*si_fac_j0*dt2
      IF (beta>0) THEN
        iso=0.25_r8*dt2*( mhd_si_iso*
     $    (si_fac_mhd*bmagsq/mu0+si_fac_pres*gamma*pres(1,:,:))+
     $      si_fac_nl*(nl_pres(1,:,:)/mu0+gamma*nl_pres(2,:,:)) )
      ELSE
        iso=0.25_r8*dt2*( mhd_si_iso*si_fac_mhd*bmagsq+
     $                    si_fac_nl*nl_pres(1,:,:) )/mu0
      ENDIF
c-----------------------------------------------------------------------
c     evaluate contributions to the 3D part of the semi-implicit
c     operator as a dyad and a vector to be multiplied by grad(A^*)^T
c     and A^*, respectively.  this is completed before the viscous
c     computations where grad(dv) gets multiplied by n.
c
c     start by finding the effective delta-f, delta-b, and delta-p that
c     are defined in the introductory comments:
c       db=curl(dvXB)/mu0-0.5*dvXJ
c       dp=gamma*P*div(dv)+0.5*dv.grad(P)-B.db
c       df=0.5*[div(dv)*grad(P)-JXcurl(dvXB)]
c     note that dbeff initially holds just curl(dvxB)
c
c     since isotropic viscous diffusion is used in nearly all nonlinear
c     computations, it is added here to avoid extra looping later.  the
c     factor of n is either in real_scal or in df_1d.
c-----------------------------------------------------------------------
      DO ip=1,nphi
        i1=ipseust-1
        DO ix=1,mpseudo
          i1=i1+1
          divtmpr=real_gvten(1,ix,ip)+real_gvten(5,ix,ip)+
     $            real_gvten(9,ix,ip)
          dbeff(1)=SUM(real_bptr(:,ix,ip)*real_gvten(1:3,ix,ip))-
     $             SUM(real_ve  (:,ix,ip)*real_grdb (1:3,ix,ip))-
     $             real_bptr(1,ix,ip)*divtmpr
          dbeff(2)=SUM(real_bptr(:,ix,ip)*real_gvten(4:6,ix,ip))-
     $             SUM(real_ve  (:,ix,ip)*real_grdb (4:6,ix,ip))-
     $             real_bptr(2,ix,ip)*divtmpr
          dbeff(3)=SUM(real_bptr(:,ix,ip)*real_gvten(7:9,ix,ip))-
     $             SUM(real_ve  (:,ix,ip)*real_grdb (7:9,ix,ip))-
     $             real_bptr(3,ix,ip)*divtmpr

          dfeff(1)=cj0*(divtmpr*real_grdp(1,ix,ip)+
     $                  dbeff(2)*real_jptr(3,ix,ip)-
     $                  dbeff(3)*real_jptr(2,ix,ip))
          dfeff(2)=cj0*(divtmpr*real_grdp(2,ix,ip)+
     $                  dbeff(3)*real_jptr(1,ix,ip)-
     $                  dbeff(1)*real_jptr(3,ix,ip))
          dfeff(3)=cj0*(divtmpr*real_grdp(3,ix,ip)+
     $                  dbeff(1)*real_jptr(2,ix,ip)-
     $                  dbeff(2)*real_jptr(1,ix,ip))

          dbeff(1)=ani*dbeff(1)+
     $             cj0*(real_ve(3,ix,ip)*real_jptr(2,ix,ip)-
     $                  real_ve(2,ix,ip)*real_jptr(3,ix,ip))
          dbeff(2)=ani*dbeff(2)+
     $             cj0*(real_ve(1,ix,ip)*real_jptr(3,ix,ip)-
     $                  real_ve(3,ix,ip)*real_jptr(1,ix,ip))
          dbeff(3)=ani*dbeff(3)+
     $             cj0*(real_ve(2,ix,ip)*real_jptr(1,ix,ip)-
     $                  real_ve(1,ix,ip)*real_jptr(2,ix,ip))

          dpeff=cp0*real_pptr(1,ix,ip)*divtmpr+
     $          cj0*SUM(real_ve(:,ix,ip)*real_grdp(:,ix,ip))-
     $          SUM(real_bptr(:,ix,ip)*dbeff)

          real_vtmp(1,ix,ip)=dfeff(1)-SUM(real_grdb(1:7:3,ix,ip)*dbeff)
          real_vtmp(2,ix,ip)=dfeff(2)-SUM(real_grdb(2:8:3,ix,ip)*dbeff)
          real_vtmp(3,ix,ip)=dfeff(3)-SUM(real_grdb(3:9:3,ix,ip)*dbeff)

          divtmpr=-twoth*divtmpr
          rscal=df_1d(i1)*real_scal(1,ix,ip)
          real_pten(1,ix,ip)=real_bptr(1,ix,ip)*dbeff(1)+dpeff+
     $      ((kin_coef+2._r8*iso_coef)*real_gvten(1,ix,ip)+
     $       iso_coef*divtmpr)*rscal
          real_pten(2,ix,ip)=real_bptr(2,ix,ip)*dbeff(1)+
     $      ((kin_coef+      iso_coef)*real_gvten(2,ix,ip)+
     $       iso_coef*real_gvten(4,ix,ip))*rscal
          real_pten(3,ix,ip)=real_bptr(3,ix,ip)*dbeff(1)+
     $      ((kin_coef+      iso_coef)*real_gvten(3,ix,ip)+
     $       iso_coef*real_gvten(7,ix,ip))*rscal
          real_pten(4,ix,ip)=real_bptr(1,ix,ip)*dbeff(2)+
     $      ((kin_coef+      iso_coef)*real_gvten(4,ix,ip)+
     $       iso_coef*real_gvten(2,ix,ip))*rscal
          real_pten(5,ix,ip)=real_bptr(2,ix,ip)*dbeff(2)+dpeff+
     $      ((kin_coef+2._r8*iso_coef)*real_gvten(5,ix,ip)+
     $       iso_coef*divtmpr)*rscal
          real_pten(6,ix,ip)=real_bptr(3,ix,ip)*dbeff(2)+
     $      ((kin_coef+      iso_coef)*real_gvten(6,ix,ip)+
     $       iso_coef*real_gvten(8,ix,ip))*rscal
          real_pten(7,ix,ip)=real_bptr(1,ix,ip)*dbeff(3)+
     $      ((kin_coef+      iso_coef)*real_gvten(7,ix,ip)+
     $       iso_coef*real_gvten(3,ix,ip))*rscal
          real_pten(8,ix,ip)=real_bptr(2,ix,ip)*dbeff(3)+
     $      ((kin_coef+      iso_coef)*real_gvten(8,ix,ip)+
     $       iso_coef*real_gvten(6,ix,ip))*rscal
          real_pten(9,ix,ip)=real_bptr(3,ix,ip)*dbeff(3)+dpeff+
     $      ((kin_coef+2._r8*iso_coef)*real_gvten(9,ix,ip)+
     $       iso_coef*divtmpr)*rscal
        ENDDO
      ENDDO
      IF (beta<=0) DEALLOCATE(real_pptr,real_grdp)
c-----------------------------------------------------------------------
c     find rho*dve product for the d/dt term, and add to real_vtmp.
c-----------------------------------------------------------------------
      IF (continuity=='full') THEN
        DO ip=1,nphi
          DO ix=1,mpseudo
            real_ve(:,ix,ip)=real_ve(:,ix,ip)*real_ndptr(1,ix,ip)*mtot
            real_vtmp(:,ix,ip)=real_vtmp(:,ix,ip)+real_ve(:,ix,ip)
          ENDDO
        ENDDO
      ELSE
        DO ip=1,nphi
          DO ix=1,mpseudo
            real_ve(:,ix,ip)=real_ve(:,ix,ip)*n0_1d(ix+ipseust-1)*mtot
            real_vtmp(:,ix,ip)=real_vtmp(:,ix,ip)+real_ve(:,ix,ip)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     evaluations for viscous dissipation:
c     first multiply grad(dv) by number density.
c-----------------------------------------------------------------------
      IF ((kin_visc>0.OR.iso_visc>0.OR.par_visc>0.OR.impladv).AND.
     $ (integrand_flag(1:4)=='visc'.OR.integrand_flag(1:3)=='all')) THEN

        IF (continuity=='full') THEN
          DO ip=1,nphi
            DO ix=1,mpseudo
              real_gvten(:,ix,ip)=
     $          real_gvten(:,ix,ip)*real_ndptr(1,ix,ip)
            ENDDO
          ENDDO
        ELSE
          DO ip=1,nphi
            DO ix=1,mpseudo
              real_gvten(:,ix,ip)=
     $          real_gvten(:,ix,ip)*n0_1d(ix+ipseust-1)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       for nonlinear parallel viscosity, find the unique components of
c
c         [ B.grad(V).B - B**2*div(V) ] * (BB - I*B**2/3)
c
c       also divide by B**4 to make unit vectors.
c
c-PRE this could be separated for efficiency.
c       multiply by the temperature-dependent profile for Braginskii
c       parallel viscosity if desired.
c-----------------------------------------------------------------------
        IF (par_visc>0) THEN
          IF (parvisc_model=='plltdep') THEN
            CALL generic_ptr_set(rb%qkappli_phi,tb%qkappli_phi,
     $                           tb%tgeom,inode,kappli,dp,dp,0_i4)
            real_scal=kappli/k_plli
          ELSE
            real_scal=1._r8
          ENDIF
          DO ip=1,nphi
            DO ix=1,mpseudo
              rb2=SUM(real_bptr(:,ix,ip)*real_bptr(:,ix,ip))+smallnum
              divtmpr=-third*rb2*(real_gvten(1,ix,ip)+
     $                real_gvten(5,ix,ip)+real_gvten(9,ix,ip))
              rscal=par_coef*real_scal(1,ix,ip)*(divtmpr+
     $              real_bptr(1,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(1,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(2,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(3,ix,ip))+
     $              real_bptr(2,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(4,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(5,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(6,ix,ip))+
     $              real_bptr(3,ix,ip)*(
     $                real_bptr(1,ix,ip)*real_gvten(7,ix,ip)+
     $                real_bptr(2,ix,ip)*real_gvten(8,ix,ip)+
     $                real_bptr(3,ix,ip)*real_gvten(9,ix,ip)))/rb2**2
              rb2=-third*rb2
              ten12=rscal*real_bptr(1,ix,ip)*real_bptr(2,ix,ip)
              ten13=rscal*real_bptr(1,ix,ip)*real_bptr(3,ix,ip)
              ten23=rscal*real_bptr(2,ix,ip)*real_bptr(3,ix,ip)
              real_pten(1,ix,ip)=real_pten(1,ix,ip)+rscal*
     $          (real_bptr(1,ix,ip)*real_bptr(1,ix,ip)+rb2)
              real_pten(2,ix,ip)=real_pten(2,ix,ip)+ten12
              real_pten(3,ix,ip)=real_pten(3,ix,ip)+ten13
              real_pten(4,ix,ip)=real_pten(4,ix,ip)+ten12
              real_pten(5,ix,ip)=real_pten(5,ix,ip)+rscal*
     $          (real_bptr(2,ix,ip)*real_bptr(2,ix,ip)+rb2)
              real_pten(6,ix,ip)=real_pten(6,ix,ip)+ten23
              real_pten(7,ix,ip)=real_pten(7,ix,ip)+ten13
              real_pten(8,ix,ip)=real_pten(8,ix,ip)+ten23
              real_pten(9,ix,ip)=real_pten(9,ix,ip)+rscal*
     $          (real_bptr(3,ix,ip)*real_bptr(3,ix,ip)+rb2)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       for nonlinear gyroviscosity, find the unique components of
c
c         (n_e*m_i*T_i/Z**2*B**4)*{ [BxW.(B**2+3BB) ] + transpose[] }
c
c-----------------------------------------------------------------------
        IF (gyr_visc>0) THEN
          vcoef=
     $      -0.125_r8*gyr_visc*dt*ms(2)*kboltz/(zeff**2*elementary_q)
          CALL generic_ptr_set(rb%qti_tot,tb%qti_tot,tb%tgeom,
     $                         inode,real_titot,dp,dp,0_i4)
          DO ip=1,nphi
            DO ix=1,mpseudo
              rb2=SUM(real_bptr(:,ix,ip)**2)+smallnum
              gcoef=vcoef*real_titot(1,ix,ip)/rb2**2
              vtmpr(1:3,1)=real_gvten(1:3,ix,ip)+real_gvten(1:7:3,ix,ip)
              vtmpr(1:3,2)=real_gvten(4:6,ix,ip)+real_gvten(2:8:3,ix,ip)
              vtmpr(1:3,3)=real_gvten(7:9,ix,ip)+real_gvten(3:9:3,ix,ip)
              btmp(1:3,1)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(1  ,ix,ip)
              btmp(1:3,2)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(2  ,ix,ip)
              btmp(1:3,3)=3._r8*real_bptr(1:3,ix,ip)*
     $                          real_bptr(3  ,ix,ip)
              btmp(1,1)=rb2+btmp(1,1)
              btmp(2,2)=rb2+btmp(2,2)
              btmp(3,3)=rb2+btmp(3,3)
              wdotr(1,1)=SUM(vtmpr(1,:)*btmp(:,1))
              wdotr(2,1)=SUM(vtmpr(2,:)*btmp(:,1))
              wdotr(3,1)=SUM(vtmpr(3,:)*btmp(:,1))
              wdotr(1,2)=SUM(vtmpr(1,:)*btmp(:,2))
              wdotr(2,2)=SUM(vtmpr(2,:)*btmp(:,2))
              wdotr(3,2)=SUM(vtmpr(3,:)*btmp(:,2))
              wdotr(1,3)=SUM(vtmpr(1,:)*btmp(:,3))
              wdotr(2,3)=SUM(vtmpr(2,:)*btmp(:,3))
              wdotr(3,3)=SUM(vtmpr(3,:)*btmp(:,3))
              bc_wdotr(1,:)=real_bptr(2,ix,ip)*wdotr(3,:)-
     $                      real_bptr(3,ix,ip)*wdotr(2,:)
              bc_wdotr(2,:)=real_bptr(3,ix,ip)*wdotr(1,:)-
     $                      real_bptr(1,ix,ip)*wdotr(3,:)
              bc_wdotr(3,:)=real_bptr(1,ix,ip)*wdotr(2,:)-
     $                      real_bptr(2,ix,ip)*wdotr(1,:)
              ten12=gcoef*(bc_wdotr(2,1)+bc_wdotr(1,2))
              ten13=gcoef*(bc_wdotr(3,1)+bc_wdotr(1,3))
              ten23=gcoef*(bc_wdotr(3,2)+bc_wdotr(2,3))
              real_pten(1,ix,ip)=real_pten(1,ix,ip)+
     $                           gcoef*(bc_wdotr(1,1)+bc_wdotr(1,1))
              real_pten(2,ix,ip)=real_pten(2,ix,ip)+ten12
              real_pten(3,ix,ip)=real_pten(3,ix,ip)+ten13
              real_pten(4,ix,ip)=real_pten(4,ix,ip)+ten12
              real_pten(5,ix,ip)=real_pten(5,ix,ip)+
     $                           gcoef*(bc_wdotr(2,2)+bc_wdotr(2,2))
              real_pten(6,ix,ip)=real_pten(6,ix,ip)+ten23
              real_pten(7,ix,ip)=real_pten(7,ix,ip)+ten13
              real_pten(8,ix,ip)=real_pten(8,ix,ip)+ten23
              real_pten(9,ix,ip)=real_pten(9,ix,ip)+
     $                           gcoef*(bc_wdotr(3,3)+bc_wdotr(3,3))
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       transform the stress/semi-implicit tensor.
c-----------------------------------------------------------------------
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,9_i4,
     $               pten,real_pten,dealiase)
c-----------------------------------------------------------------------
c       terms for implicit advection:
c
c          0.5*dt*rho*( dV.grad(V_old) + V_old.grad(dV) )
c
c       where dV is the iterate from work4.
c
c       note that V_old already includes any equilibrium flow.  Also,
c       real_gvten still holds n*grad(dV) in real space, and real_ve
c       holds  mtot*n*dV in real space.
c
c       add this advective force to the vector semi-implicit
c       contribution real_vtmp, which also has mtot*n*dV.
c-----------------------------------------------------------------------
        IF (impladv.AND.(advect=='V only'.OR.advect=='all')) THEN
          DO ip=1,nphi
            DO ix=1,mpseudo
               tmp1=hdt*real_ve(:,ix,ip)
               tmp2=hdt*mtot*real_oldv(:,ix,ip)
               real_vtmp(1,ix,ip)=real_vtmp(1,ix,ip)+
     $                   tmp1(1)*real_oldgv(1,ix,ip)+
     $                   tmp1(2)*real_oldgv(2,ix,ip)+
     $                   tmp1(3)*real_oldgv(3,ix,ip)+
     $                   tmp2(1)*real_gvten(1,ix,ip)+
     $                   tmp2(2)*real_gvten(2,ix,ip)+
     $                   tmp2(3)*real_gvten(3,ix,ip)
               real_vtmp(2,ix,ip)=real_vtmp(2,ix,ip)+
     $                   tmp1(1)*real_oldgv(4,ix,ip)+
     $                   tmp1(2)*real_oldgv(5,ix,ip)+
     $                   tmp1(3)*real_oldgv(6,ix,ip)+
     $                   tmp2(1)*real_gvten(4,ix,ip)+
     $                   tmp2(2)*real_gvten(5,ix,ip)+
     $                   tmp2(3)*real_gvten(6,ix,ip)
               real_vtmp(3,ix,ip)=real_vtmp(3,ix,ip)+
     $                   tmp1(1)*real_oldgv(7,ix,ip)+
     $                   tmp1(2)*real_oldgv(8,ix,ip)+
     $                   tmp1(3)*real_oldgv(9,ix,ip)+
     $                   tmp2(1)*real_gvten(7,ix,ip)+
     $                   tmp2(2)*real_gvten(8,ix,ip)+
     $                   tmp2(3)*real_gvten(9,ix,ip)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     transform the vector contribution from the semi-implicit operator,
c     the d/dt term, and implicit advection.
c-----------------------------------------------------------------------
      CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,
     $             force,real_vtmp,dealiase)
c-----------------------------------------------------------------------
c     the resultant is formed by finding
c
c       A^*.force + grad(A^*):pten plus isotropic contributions.
c
c     for each of the test vectors, one for each vector-direction /
c     spatial basis function.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
c-----------------------------------------------------------------------
c       find the contributions of P-total times div(dV) and B.curl(dV)
c       from the auxiliary discontinuous fields.  add the contributions
c       to the pten tensor so that grad(A^*)^T:pten includes
c       the div(A^*) and B.curl(A^*)=grad(A^*)^T:(BxI)^T from the
c       auxiliary forcing.
c
c-PRE   this could be done in the viscous-dissipation loop.
c-----------------------------------------------------------------------
        IF (poly_divv>0) THEN
          DO iy=1,ncy
            DO ix=1,ncx
              sqpaux=dtdv*sqptot(ix,iy)*aux(1,ix,iy,imode)
              baux(1)=dtpv*be(1,ix,iy)*aux(2,ix,iy,imode)
              baux(2)=dtpv*be(2,ix,iy)*aux(2,ix,iy,imode)
              baux(3)=dtpv*be(3,ix,iy)*aux(2,ix,iy,imode)
              pten(1,ix,iy,imode)=pten(1,ix,iy,imode)+sqpaux
              pten(2,ix,iy,imode)=pten(2,ix,iy,imode)-baux(3)
              pten(3,ix,iy,imode)=pten(3,ix,iy,imode)+baux(2)
              pten(4,ix,iy,imode)=pten(4,ix,iy,imode)+baux(3)
              pten(5,ix,iy,imode)=pten(5,ix,iy,imode)+sqpaux
              pten(6,ix,iy,imode)=pten(6,ix,iy,imode)-baux(1)
              pten(7,ix,iy,imode)=pten(7,ix,iy,imode)-baux(2)
              pten(8,ix,iy,imode)=pten(8,ix,iy,imode)+baux(1)
              pten(9,ix,iy,imode)=pten(9,ix,iy,imode)+sqpaux
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       complete the factors that are dotted into the test vectors, A.
c-----------------------------------------------------------------------
        IF (geom=='tor') THEN
          forck(1,:,:)=force(1,:,:,imode)+
     $      iso*((1+k2ef(imode))*ve(1,:,:,imode)+
     $        2*(0,1)*keff(imode)*ve(3,:,:,imode))/bigr2
     $      +(                  pten(9,:,:,imode)
     $       -(0,1)*keff(imode)*pten(3,:,:,imode))/bigr
          forck(2,:,:)=force(2,:,:,imode)+
     $      iso*k2ef(imode)*ve(2,:,:,imode)/bigr2
     $       -(0,1)*keff(imode)*pten(6,:,:,imode) /bigr
          forck(3,:,:)=force(3,:,:,imode)+
     $      iso*((1+k2ef(imode))*ve(3,:,:,imode)-
     $        2*(0,1)*keff(imode)*ve(1,:,:,imode))/bigr2
     $      -(                  pten(3,:,:,imode)
     $       +(0,1)*keff(imode)*pten(9,:,:,imode))/bigr
        ELSE
          forck(1,:,:)=force(1,:,:,imode)+
     $                 iso*k2ef(imode)*ve(1,:,:,imode)-
     $                 (0,1)*keff(imode)*pten(3,:,:,imode)
          forck(2,:,:)=force(2,:,:,imode)+
     $                 iso*k2ef(imode)*ve(2,:,:,imode)-
     $                 (0,1)*keff(imode)*pten(6,:,:,imode)
          forck(3,:,:)=force(3,:,:,imode)+
     $                 iso*k2ef(imode)*ve(3,:,:,imode)-
     $                 (0,1)*keff(imode)*pten(9,:,:,imode)
        ENDIF
c-----------------------------------------------------------------------
c       accumulate all contributions.
c-----------------------------------------------------------------------
        DO iv=1,nvc
          DO iy=1,ncy
            int(:,iy,iv,imode)=0._r8
            DO ix=1,ncx
              int(1,iy,iv,imode)=int(1,iy,iv,imode)+
     $          alpha(ix,iy,iv)*forck(1,ix,iy)
     $          +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(1,ix,iy,imode)
     $                      +dalpdz(ix,iy,iv)*vez(1,ix,iy,imode))
     $          +dalpdr(ix,iy,iv)*pten(1,ix,iy,imode)
     $          +dalpdz(ix,iy,iv)*pten(2,ix,iy,imode)

              int(2,iy,iv,imode)=int(2,iy,iv,imode)+
     $          alpha(ix,iy,iv)*forck(2,ix,iy)
     $          +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(2,ix,iy,imode)
     $                      +dalpdz(ix,iy,iv)*vez(2,ix,iy,imode))
     $          +dalpdr(ix,iy,iv)*pten(4,ix,iy,imode)
     $          +dalpdz(ix,iy,iv)*pten(5,ix,iy,imode)

              int(3,iy,iv,imode)=int(3,iy,iv,imode)+
     $          alpha(ix,iy,iv)*forck(3,ix,iy)
     $          +iso(ix,iy)*(dalpdr(ix,iy,iv)*ver(3,ix,iy,imode)
     $                      +dalpdz(ix,iy,iv)*vez(3,ix,iy,imode))
     $          +dalpdr(ix,iy,iv)*pten(7,ix,iy,imode)
     $          +dalpdz(ix,iy,iv)*pten(8,ix,iy,imode)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       complete the integrand for the discontinuous scalar equations.
c-----------------------------------------------------------------------
        IF (poly_divv>0) THEN
          divv=sqptot*dtdv*(ver(1,:,:,imode)+vez(2,:,:,imode)+
     $                      ((0,1)*keff(imode)*ve(3,:,:,imode)+
     $                                       g*ve(1,:,:,imode))/bigr)
          bvort=(be(1,:,:)*(vez(3,:,:,imode)
     $                     -(0,1)*keff(imode)*ve(2,:,:,imode)/bigr)+
     $           be(2,:,:)*(((0,1)*keff(imode)*ve(1,:,:,imode)-
     $                      g*ve(3,:,:,imode))/bigr-ver(3,:,:,imode))+
     $           be(3,:,:)*(ver(2,:,:,imode)-vez(1,:,:,imode)))*dtpv
          DO iv=nvc+1,nvc+nvm  !  auxiliary equation
            int(1,:,iv,imode)=
     $        SUM(almod(:,:,iv-nvc)*(aux(1,:,:,imode)-divv),1)
            int(2,:,iv,imode)=
     $        SUM(almod(:,:,iv-nvc)*(aux(2,:,:,imode)-bvort),1)
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE v_3dsi_dot
c-----------------------------------------------------------------------
c     subprogram 6. cont_dot.
c     compute the matrix-vector dot product for the continuity advance
c     with implicit advection.
c-----------------------------------------------------------------------
      SUBROUTINE cont_dot(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          nd_eq,real_ve,dart,upwc,dp
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: kap_art
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_vec,real_gn
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal

      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: nd,
     $             ndr,ndz
      COMPLEX(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: v2,
     $             v2r,v2z
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             nvel,grad_nd

      INTEGER(i4) :: nv,iv,imode,ncx,ncy,ix,ip
      REAL(r8) :: hdt,fdt,dfac,fupw,fupwi,rdot,rscl,fhyp
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      hdt=0.5_r8*dt
      dfac=fthc*dt*nd_diff
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate velocity.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,inode,
     $                     real_ve,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     density: the iterate for the change is put in nd,ndr,ndz.
c
c     if hyper-diffusivity is used, the iterate is a 2-vector.  number
c     density is the first quantity, and the auxiliary scalar, which is
c     proportional to grad**2(n) is the second quantity.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        fhyp=SQRT(nd_hypd*dt)
        CALL generic_all_eval(rb%work6,tb%work6,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,v2,v2r,v2z,2_i4)
        nd(1,:,:,:)=v2(1,:,:,:)
      ELSE
        CALL generic_all_eval(rb%work2,tb%work2,rb%dxdr,rb%dydr,
     $                        rb%dxdz,rb%dydz,rb%xg,rb%yg,tb%tgeom,
     $                        tb%ng,nd,ndr,ndz,1_i4)
      ENDIF
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      IF (nd_floor>0.AND.nd_diff>0) THEN
        CALL generic_ptr_set(rb%qdart,tb%qdart,tb%tgeom,
     $                       inode,dart,dp,dp,0_i4)
        kap_art=dart(1,:,:)+dfac
      ENDIF
      IF (nd_dart_upw>0) THEN
        CALL generic_ptr_set(rb%qupw_phi,tb%qupw_phi,tb%tgeom,
     $                       inode,upwc,dp,dp,0_i4)
      ENDIF
      IF (nd_hypd>0) THEN
        DO imode=1,nmodes
          grad_nd(1,:,:,imode)=v2r(1,:,:,imode)
          grad_nd(2,:,:,imode)=v2z(1,:,:,imode)
          grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
        ENDDO
      ELSE IF (nd_diff>0.OR.nd_dart_upw>0) THEN
        DO imode=1,nmodes
          grad_nd(1,:,:,imode)=ndr(1,:,:,imode)
          grad_nd(2,:,:,imode)=ndz(1,:,:,imode)
          grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the nonlinear contributions to n*V.
c     note that the real arrays have three dimensions, where the
c     first is the vector component, the second covers the poloidal
c     plane, and the third is the phi index,
c
c     real_ve already contains any equilibrium flow.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,nd,real_scal,
     $             dealiase)
c-----------------------------------------------------------------------
c     include the upwinding diffusive flux that is proportional to
c     VV.grad(n).
c-----------------------------------------------------------------------
      IF (nd_dart_upw>0) THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $               real_gn,dealiase)
        fupw =nd_dart_upw*dt**2*upw_aniso
        fupwi=nd_dart_upw*dt**2*(1._r8-upw_aniso)
        DO ip=1,nphi
          DO ix=1,mpseudo
            rdot=fupw*upwc(1,ix,ip)*
     $          SUM(real_ve(:,ix,ip)*real_gn(:,ix,ip)) /
     $        ( SUM(real_ve(:,ix,ip)*real_ve(:,ix,ip)) + smallnum )
            rscl=fupwi*upwc(1,ix,ip)
            real_vec(:,ix,ip)=hdt*real_ve(:,ix,ip)*real_scal(1,ix,ip)-
     $        ( rdot*real_ve(:,ix,ip)+rscl*real_gn(:,ix,ip) )
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     or just V*dn
c-----------------------------------------------------------------------
      ELSE
        real_vec(1,:,:)=hdt*real_ve(1,:,:)*real_scal(1,:,:)
        real_vec(2,:,:)=hdt*real_ve(2,:,:)*real_scal(1,:,:)
        real_vec(3,:,:)=hdt*real_ve(3,:,:)*real_scal(1,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     transform flux.
c-----------------------------------------------------------------------
      CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,nvel,real_vec,
     $             dealiase)
c-----------------------------------------------------------------------
c     begin mode loop.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       add the diffusive term -D*grad(n) and the hyper-diffusive
c       part of the flux +D_hyp*grad(f) if needed.
c-----------------------------------------------------------------------
        IF (nd_floor>0.AND.nd_diff>0) THEN
          nvel(1,:,:,imode)=nvel(1,:,:,imode)-
     $                kap_art*grad_nd(1,:,:,imode)
          nvel(2,:,:,imode)=nvel(2,:,:,imode)-
     $                kap_art*grad_nd(2,:,:,imode)
          nvel(3,:,:,imode)=nvel(3,:,:,imode)-
     $                kap_art*grad_nd(3,:,:,imode)
        ELSE IF (nd_diff>0) THEN
          nvel(1,:,:,imode)=nvel(1,:,:,imode)-dfac*grad_nd(1,:,:,imode)
          nvel(2,:,:,imode)=nvel(2,:,:,imode)-dfac*grad_nd(2,:,:,imode)
          nvel(3,:,:,imode)=nvel(3,:,:,imode)-dfac*grad_nd(3,:,:,imode)
        ENDIF
        IF (nd_hypd>0) THEN
          nvel(1,:,:,imode)=nvel(1,:,:,imode)+fhyp*v2r(2,:,:,imode)
          nvel(2,:,:,imode)=nvel(2,:,:,imode)+fhyp*v2z(2,:,:,imode)
          nvel(3,:,:,imode)=nvel(3,:,:,imode)+
     $                   fhyp*v2(2,:,:,imode)*(0,1)*keff(imode)/bigr
        ENDIF
c-----------------------------------------------------------------------
c       construct n-f*dt*grad(alpha*).( n*V - D*grad(n) ).
c
c       if hyper-diffusivity is used, also compute the off-diagonal
c       contributions for +fhyp*grad**2(aux-scalar) in the n rows and
c       -fhyp*grad**2(n) in the aux-scalar rows.
c-----------------------------------------------------------------------
        nd(1,:,:,imode)=  nd(1,:,:,imode)+
     $                  nvel(3,:,:,imode)*(0,1)*keff(imode)/bigr
        IF (nd_hypd>0) THEN
          v2(2,:,:,imode)=v2(2,:,:,imode)-
     $          fhyp*grad_nd(3,:,:,imode)*(0,1)*keff(imode)/bigr
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(alpha(:,:,iv)*  nd(1,:,:,imode)-
     $                           dalpdr(:,:,iv)*nvel(1,:,:,imode)-
     $                           dalpdz(:,:,iv)*nvel(2,:,:,imode),1)
            int(2,:,iv,imode)=SUM(alpha(:,:,iv)*     v2(2,:,:,imode)+
     $                     fhyp*(dalpdr(:,:,iv)*grad_nd(1,:,:,imode)+
     $                           dalpdz(:,:,iv)*grad_nd(2,:,:,imode)),1)
          ENDDO
        ELSE
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(alpha(:,:,iv)*  nd(1,:,:,imode)-
     $                           dalpdr(:,:,iv)*nvel(1,:,:,imode)-
     $                           dalpdz(:,:,iv)*nvel(2,:,:,imode),1)
          ENDDO
        ENDIF
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cont_dot
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE integrands_dot
