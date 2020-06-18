c-----------------------------------------------------------------------
c     file integrands_rhs.f
c     module that includes integrand routines used in the finite-
c     element integrations for creating the right-hand-side vectors
c     of the algebraic systems.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module integrands_rhs
c     subprogram 1. get_vol.
c     subprogram 2. curl
c     subprogram 3. brhs_mhd.
c     subprogram 4. brhs_hyp.
c     subprogram 5. brhs_hmhd.
c     subprogram 6. hall_cor.
c     subprogram 7. tirhs.
c     subprogram 8. terhs.
c     subprogram 9. vrhs.
c     subprogram 10. advect_cor.
c     subprogram 11. divb_rhs.
c     subprogram 12. ndrhs.
c     subprogram 13. scal_lapl_rhs.
c-----------------------------------------------------------------------
c     module containing integrands for the rhs of equations.
c-----------------------------------------------------------------------
      MODULE integrands_rhs

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
c     subprogram 1. get_vol.
c     compute the cell volume less a factor of 2*pi.
c-----------------------------------------------------------------------
      SUBROUTINE get_vol(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: wj
c-----------------------------------------------------------------------
c     put dVol and dArea in the vertex index used for cell-centered
c     integrals.  dArea/r is used for the computation of F.
c
c     the sum of the test functions produces the quadrature weight times
c     the Jacobian.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      wj=SUM(alpha,3)

      int(1,:,1)=SUM(wj,1)
      int(2,:,1)=SUM(wj/bigr,1)
      int(3,:,1)=SUM(wj/bigr**2,1)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_vol
c-----------------------------------------------------------------------
c     subprogram 2. curl
c     form the integrand for a curl operation.
c     the incoming vector and its derivatives have been evaluated at
c     the quadrature points and stored in the qwork1 structure.
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
      REAL(r8) :: tor_fac
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             vec,dvecr,dvecz
      INTEGER(i4) :: iv,im,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate the
c     incoming vector.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,0_i4,poly_degree)
      CALL generic_ptr_set(rb%qwork1,tb%qwork1,tb%tgeom,inode,
     $                     vec,dvecr,dvecz,1_i4)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set factor to collect extra terms according to geom.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        tor_fac=1._r8
      ELSE
        tor_fac=0._r8
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
     $                          -(0,1)*keff(im)*vec(2,:,:,im)/bigr ),1)
          int(2,:,iv,im)=SUM(-alpha(:,:,iv)*( dvecr(3,:,:,im)
     $      +(tor_fac*vec(3,:,:,im)-
     $        (0,1)*keff(im)*vec(1,:,:,im))/bigr ),1)
          int(3,:,iv,im)=SUM(alpha(:,:,iv)*( dvecr(2,:,:,im)
     $                                      -dvecz(1,:,:,im) ),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curl
c-----------------------------------------------------------------------
c     subprogram 3. brhs_mhd.
c     compute the integrand used in the rhs of the semi-implicit
c     magnetic field equation for the mhd time-split.
c-----------------------------------------------------------------------
      SUBROUTINE brhs_mhd(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          aldis,dalddr,dalddz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          be_eq,ve_eq,e_a,e_ar,
     $          e_az,elecd_phi,elecd_n0,elecd_eq,ja_eq,tele_eq,dp,ds2
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_be,real_eexp
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_j
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: diff
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: eleq_1d
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2)) :: jaeq_1d 
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,ber,bez,ve,tele,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             ja_old,eexp
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2),nmodes) :: divb
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: auxb
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,jq,nv,npol,ncx,ncy,ip,nvc,nvd
      REAL(r8) :: disdbfac,fhyp,fhdb
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and set pointers to
c     the velocity and magnetic field at the beginning of the time
c     split.  current density for resistive diffusion must be based on
c     this magnetic field.
c
c     if poly_divb is non-negative, divergence error is controlled with
c     an auxiliary field that is discontinuous at element borders, and
c     these bases are also needed.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      IF (poly_divb>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',aldis,
     $                          dalddr,dalddz,0_i4,poly_divb,
     $                          polydmin=poly_divb_min,
     $                          polydmax=poly_divb_max)
        nvd=SIZE(aldis,3)
      ENDIF

      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,dcp,dcp,0_i4)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
      IF ((elecd>0.OR.(hyp_eta>0.OR.hyp_dbd>0).AND.
     $    .NOT.split_hypeta).AND.
     $    integrand_flag(1:3)=='mhd') THEN
        CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,
     $                 ja_old,1._r8/mu0)
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_eq,tb%qelecd_eq,tb%tgeom,
     $                         inode,elecd_eq,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                         inode,ja_eq,dp,dp,0_i4)
        ENDIF
        IF (.NOT.nonlinear.AND.eta_model=='eta full') THEN
          CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                         inode,tele_eq,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,
     $                         inode,tele,dcp,dcp,0_i4)
        ELSE IF (threedeta) THEN
          CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                         inode,elecd_phi,dp,dp,0_i4)
          ALLOCATE(real_j(3,mpseudo,nphi))
          eleq_1d=RESHAPE(elecd_eq,(/ncx*ncy/))
          jaeq_1d=RESHAPE( ja_eq,(/3,ncx*ncy/))
        ELSE IF (eta_model=='eta n=0 only') THEN
          CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                         inode,elecd_n0,dp,dp,0_i4)
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          diff=elecd*MAX(ds2(1,:,:),0._r8)
        ELSE
          diff=elecd
        ENDIF
      ENDIF
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find div(b) from the old B if cleaning is not split.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb.OR.poly_divb>=0) THEN
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            divb(:,:,imode)=ber(1,:,:,imode)+bez(2,:,:,imode)
     $                    +((0,1)*keff(imode)*be(3,:,:,imode)
     $                                       +be(1,:,:,imode))/bigr
          ELSE
            divb(:,:,imode)=ber(1,:,:,imode)+bez(2,:,:,imode)
     $                     +(0,1)*keff(imode)*be(3,:,:,imode)/bigr
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     reset magnetic field to the fb_vxb-centered predicted level if it
c     is needed for an advection corrector step.
c-----------------------------------------------------------------------
      IF ((nonlinear.OR.eq_flow/='none').AND.
     $    integrand_flag/='mhd_predict'.AND..NOT.impladv) THEN
        CALL generic_ptr_set(rb%qwork1,tb%qwork1,tb%tgeom,inode,
     $                       be,ber,bez,1_i4)
      ENDIF
c-----------------------------------------------------------------------
c     find the explicit nonlinear contributions to the electric field.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
c-----------------------------------------------------------------------
c       first find b_perturbed and the perturbed velocity in real space.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,be,real_be,
     $               dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $               dealiase)
c-----------------------------------------------------------------------
c       find dt * b_perturbed X perturbed velocity in real space
c-----------------------------------------------------------------------
        CALL math_cart_cross(real_eexp,real_be,real_ve,dt)
c-----------------------------------------------------------------------
c       add the explicit contribution to resistive dissipation, if
c       3D resistivity is used.
c-----------------------------------------------------------------------
        IF (elecd>0.AND.threedeta) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,
     $                 ja_old,real_j,dealiase)
          DO ip=1,nphi
            DO jq=1,3
              real_eexp(jq,:,ip)=real_eexp(jq,:,ip)+
     $          dt*mu0*( elecd_phi(1,:,ip)*real_j(jq,:,ip)
     $                 +(elecd_phi(1,:,ip)-eleq_1d(ipseust:ipseuen))
     $                  *jaeq_1d(jq,ipseust:ipseuen))
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       transform to Fourier space.
c-----------------------------------------------------------------------
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,eexp,real_eexp,
     $               dealiase)
      ELSE
        eexp=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     begin the mode loop for linear contributions.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       find the explicit linear contributions to the electric field.
c-----------------------------------------------------------------------
        CALL math_cadd_cross(eexp(:,:,:,imode),be_eq,ve(:,:,:,imode),dt)
        IF (eq_flow/='none')
     $    CALL math_cadd_cross(eexp(:,:,:,imode),
     $                         be(:,:,:,imode),ve_eq,dt)
c-----------------------------------------------------------------------
c       add the explicit contribution to resistive dissipation, where
c       the electrical diffusivity is 2D.
c-----------------------------------------------------------------------
        IF (elecd>0.AND.integrand_flag(1:3)=='mhd') THEN
          IF (eta_model=='eta n=0 only') THEN
            DO jq=1,3
              eexp(jq,:,:,imode)=eexp(jq,:,:,imode)
     $          +dt*mu0*elecd_n0(1,:,:)*ja_old(jq,:,:,imode)
            ENDDO
            IF (nonlinear.AND.keff(imode)==0) THEN
              DO jq=1,3
                eexp(jq,:,:,imode)=eexp(jq,:,:,imode)+
     $            dt*mu0*(elecd_n0(1,:,:)-elecd_eq(1,:,:))*ja_eq(jq,:,:)
              ENDDO
            ENDIF
          ELSE IF (eta_model=='eta full'.OR.eta_model=='chodura') THEN
            IF (.NOT.nonlinear.AND.eta_model=='eta full') THEN
              DO jq=1,3
                eexp(jq,:,:,imode)=eexp(jq,:,:,imode)+
     $            dt*mu0*elecd_eq(1,:,:)*(ja_old(jq,:,:,imode)-
     $               1.5*ja_eq(jq,:,:)*tele(1,:,:,imode)/tele_eq(1,:,:))
              ENDDO
            ENDIF
          ELSE
            DO jq=1,3
              eexp(jq,:,:,imode)=eexp(jq,:,:,imode)
     $          +dt*mu0*diff*ja_old(jq,:,:,imode)
            ENDDO
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       construct -curl(alpha^*).E_explicit  [ -div(alpha)*div(B) ].
c-----------------------------------------------------------------------
        IF (divbd>0.AND..NOT.split_divb) THEN
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(   dalpdz(:,:,iv)*eexp(3,:,:,imode)
     $        +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(2,:,:,imode)
     $           -divbd*dt*dalpdrc(:,:,iv)*divb(:,:,imode),1)
            int(2,:,iv,imode)=SUM(  -dalpdr(:,:,iv)*eexp(3,:,:,imode)
     $        -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(1,:,:,imode)
     $        -divbd*dt*dalpdz (:,:,iv)*divb(:,:,imode),1)
            int(3,:,iv,imode)=SUM(  -dalpdz (:,:,iv)*eexp(1,:,:,imode)
     $                              +dalpdrc(:,:,iv)*eexp(2,:,:,imode)
     $                 +divbd*dt*(0,1)*alpha(:,:,iv)*keff(imode)
     $                          *divb(:,:,imode)/bigr,1)
          ENDDO
        ELSE
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(   dalpdz(:,:,iv)*eexp(3,:,:,imode)
     $        +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(2,:,:,imode),1)
            int(2,:,iv,imode)=SUM(  -dalpdr(:,:,iv)*eexp(3,:,:,imode)
     $        -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(1,:,:,imode),1)
            int(3,:,iv,imode)=SUM(  -dalpdz(:,:,iv)*eexp(1,:,:,imode)
     $                             +dalpdrc(:,:,iv)*eexp(2,:,:,imode),1)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       find the integrand for the discontinuous scalar equation,
c       disdbfac*aldis*div(b), where aldis is the value of the
c       discontinuous basis function.
c-----------------------------------------------------------------------
        IF (poly_divb>=0) THEN
          disdbfac=SQRT(disc_dbd*dt/fdivb)
          DO iv=nvc+1,nv
            int(1,:,iv,imode)=
     $        disdbfac*SUM(aldis(:,:,iv-nvc)*divb(:,:,imode),1)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       add -alpha*curl(E_applied) to the n=0 mode for cases without
c       an equilibrium field.  note:  the third component of e_applied
c       is 1/(2*pi) or 1/per_length at the surface.
c-----------------------------------------------------------------------
c-PRE   IF () THEN
c         CALL generic_ptr_set(rb%qe_applied,tb%qe_applied,tb%tgeom,
c    $                         inode,e_a,e_ar,e_az,1_i4)
c         DO iv=1,nv
c           int(1,:,iv,imode)=int(1,:,iv,imode)
c    $        -SUM(alpha(:,:,iv)*dt*volt*e_az(3,:,:)/bigr,1)
c           int(2,:,iv,imode)=int(2,:,iv,imode)
c    $        +SUM(alpha(:,:,iv)*dt*volt*e_ar(3,:,:)/bigr,1)
c           int(3,:,iv,imode)=int(3,:,iv,imode)
c    $        -SUM(alpha(:,:,iv)*dt*e_vert*e_ar(2,:,:),1)
c         ENDDO
c       ENDIF
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     if hyper-resistivity is used without a split step, an auxiliary
c     equation is solved simultaneously with dB.  the integrand for the
c     right side of this equation is
c
c       -SQRT(dt*hyp_eta/fhyp_eta)*curl(C^*).curl(B_old)
c       -SQRT(dt*hyp_dbd/fhyp_dbd)*div(C^*).div(B_old)
c
c     where C is the test vector for the auxiliary field.  note that
c     ja_old and divb are modified at this point.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0.OR.hyp_dbd>0).AND..NOT.split_hypeta) THEN
        fhyp=SQRT(dt*hyp_eta/fhyp_eta)*mu0
        fhdb=SQRT(dt*hyp_dbd/fhyp_dbd)
        DO imode=1,nmodes
          ja_old(:,:,:,imode)=ja_old(:,:,:,imode)*fhyp
          divb(:,:,imode)=divb(:,:,imode)*fhdb
          DO iv=1,nvc
            int(4,:,iv,imode)=
     $        SUM( dalpdz(:,:,iv)*ja_old(3,:,:,imode)+
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*
     $             ja_old(2,:,:,imode)-
     $             dalpdrc(:,:,iv)*divb(:,:,imode),1)
            int(5,:,iv,imode)=
     $        SUM(-dalpdr(:,:,iv)*ja_old(3,:,:,imode)-
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*
     $             ja_old(1,:,:,imode)-
     $             dalpdz(:,:,iv)*divb(:,:,imode),1)
            int(6,:,iv,imode)=
     $        SUM(-dalpdz (:,:,iv)*ja_old(1,:,:,imode)+
     $             dalpdrc(:,:,iv)*ja_old(2,:,:,imode)+
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*
     $             divb(:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     deallocate arrays as needed.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:3)=='mhd'.AND.threedeta) DEALLOCATE(real_j)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE brhs_mhd
c-----------------------------------------------------------------------
c     subprogram 4. brhs_hyp.
c     compute the integrand used in the rhs of the hyper-resistivity
c     time split.
c-----------------------------------------------------------------------
      SUBROUTINE brhs_hyp(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: be,ber,bez
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ja
      COMPLEX(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: divb
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,nvc
      REAL(r8) :: fhyp,fhdb
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and set pointers to
c     the magnetic field at the beginning of the time split.  
c     current density for resistive diffusion is based on
c     this magnetic field.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
c-----------------------------------------------------------------------
c     convenience parameters and curl(B):
c-----------------------------------------------------------------------
      nvc=SIZE(alpha,3)
      fhyp=SQRT(dt*hyp_eta/fhyp_eta)
      fhdb=SQRT(dt*hyp_dbd/fhyp_dbd)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,fhyp)
c-----------------------------------------------------------------------
c     start the mode loop with the divergence of this component.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (geom=='tor') THEN
          divb(:,:)=
     $      fhdb*( ber(1,:,:,imode)+bez(2,:,:,imode)
     $           +((0,1)*keff(imode)*be(3,:,:,imode)
     $                              +be(1,:,:,imode))/bigr)
        ELSE
          divb(:,:)=fhdb*(ber(1,:,:,imode)+bez(2,:,:,imode)
     $                   +(0,1)*keff(imode)*be(3,:,:,imode)/bigr)
        ENDIF
c-----------------------------------------------------------------------
c       if hyper-resistivity is used without a split step, an auxiliary
c       equation is solved simultaneously with dB.  the integrand for
c       the right side of this equation is
c
c         -SQRT(dt*hyp_eta/fhyp_eta)*curl(C^*).curl(B_old)
c         -SQRT(dt*hyp_dbd/fhyp_dbd)*div(C^*).div(B_old)
c
c       where C is the test vector for the auxiliary field.  note that
c-----------------------------------------------------------------------
        DO iv=1,nvc
          int(1:3,:,iv,imode)=0._r8
          int(4,:,iv,imode)=
     $      SUM( dalpdz(:,:,iv)*ja(3,:,:,imode)+
     $           (0,1)*alpha(:,:,iv)*keff(imode)/bigr*ja(2,:,:,imode)-
     $           dalpdrc(:,:,iv)*divb,1)
          int(5,:,iv,imode)=
     $      SUM(-dalpdr(:,:,iv)*ja(3,:,:,imode)-
     $           (0,1)*alpha(:,:,iv)*keff(imode)/bigr*ja(1,:,:,imode)-
     $           dalpdz(:,:,iv)*divb,1)
          int(6,:,iv,imode)=
     $      SUM(-dalpdz (:,:,iv)*ja(1,:,:,imode)+
     $           dalpdrc(:,:,iv)*ja(2,:,:,imode)+
     $           (0,1)*alpha(:,:,iv)*keff(imode)/bigr*divb,1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE brhs_hyp
c-----------------------------------------------------------------------
c     subprogram 5. brhs_hmhd.
c     compute the integrand used in the rhs of the semi-implicit
c     magnetic field equation for the combined Hall-MHD time-split.
c-----------------------------------------------------------------------
      SUBROUTINE brhs_hmhd(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          aldis,dalddr,dalddz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: be_eq,ja_eq,
     $          jaeq_r,jaeq_z,nd_eq,nd_n0,ve_eq,grdv_eq,
     $          real_ndptr,dp,ds2,elecd_n0,nd_eqr,nd_eqz,
     $          eq_force,elecd_eq,elecd_phi,tele_eq
      REAL(r8), DIMENSION(3,mpseudo,nphi), TARGET :: real_ja,real_be,
     $          real_eexp,real_ve
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_te
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: eleq_1d
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2)) :: jaeq_1d 
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_divv
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          real_gradv,real_gradj,divv_eq
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: diff
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: pfac,dbdt,dtm,disdbfac,fhyp,fhdb

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,ber,bez,ve,ver,vez,
     $             jafe,jar,jaz,nd,ndr,ndz,tele,teler,telez,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             jaore,eexp,ebxv,ja,grad_nd
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2),nmodes) :: divb
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: auxb
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: dbqd
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: jap,vep,divv
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,iq,jq,nv,ncx,ncy,in0,ip,nvc,nvd
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      pfac=coefgpe*kboltz
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the 0-th
c     order magnetic field.  evaluate the old current density
c     from the gradients of the old magnetic fields for resitive
c     dissipation.  for the mhd part, evaluate flow velocities.
c
c     if poly_divb is non-negative, divergence error is controlled with
c     an auxiliary field that is discontinuous at element borders, and
c     these bases are also needed.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      IF (poly_divb>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',aldis,
     $                          dalddr,dalddz,0_i4,poly_divb,
     $                          polydmin=poly_divb_min,
     $                          polydmax=poly_divb_max)
        nvd=SIZE(aldis,3)
      ENDIF
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,
     $                     be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8/mu0)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,ja_eq,jaeq_r,jaeq_z,1_i4)
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,nd_eqr,nd_eqz,1_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,
     $                     nd,ndr,ndz,1_i4)
      IF (beta>0) THEN
        CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,inode,
     $                       tele,teler,telez,1_i4)
        CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                       inode,tele_eq,dp,dp,0_i4)
      ENDIF
      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,ver,vez,1_i4)
      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     dissipation coefficients.
c-----------------------------------------------------------------------
      dbdt=dt*divbd
      dtm=dt*mu0
      IF (elecd>0) THEN
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_eq,tb%qelecd_eq,tb%tgeom,
     $                         inode,elecd_eq,dp,dp,0_i4)
        ENDIF
        IF (.NOT.nonlinear.AND.eta_model=='eta full') THEN
          CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,
     $                         inode,tele,dcp,dcp,0_i4)
        ELSE IF (threedeta) THEN
          CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                         inode,elecd_phi,dp,dp,0_i4)
          eleq_1d=RESHAPE(elecd_eq,(/ncx*ncy/))
          jaeq_1d=RESHAPE( ja_eq,(/3,ncx*ncy/))
        ELSE IF (eta_model=='eta n=0 only') THEN
          CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                         inode,elecd_n0,dp,dp,0_i4)
          diff=dtm*elecd_n0(1,:,:)
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          diff=dtm*elecd*MAX(ds2(1,:,:),0._r8)
        ELSE
          diff=dtm*elecd
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     find div(b) for cleaning if not split.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb.OR.poly_divb>=0) THEN
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            divb(:,:,imode)=(ber(1,:,:,imode)+bez(2,:,:,imode)
     $                     +((0,1)*keff(imode)*be(3,:,:,imode)
     $                                        +be(1,:,:,imode))/bigr)
          ELSE
            divb(:,:,imode)=(ber(1,:,:,imode)+bez(2,:,:,imode)
     $                      +(0,1)*keff(imode)*be(3,:,:,imode)/bigr)
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find grad(n) for the grad(pe) term.  this uses the time-averaged
c     number density.
c-----------------------------------------------------------------------
      IF (beta>0) THEN
        DO imode=1,nmodes
          grad_nd(1,:,:,imode)=pfac*ndr(1,:,:,imode)
          grad_nd(2,:,:,imode)=pfac*ndz(1,:,:,imode)
          grad_nd(3,:,:,imode)=(0,1)*pfac*keff(imode)*
     $                               nd(1,:,:,imode)/bigr
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the explicit nonlinear contributions to the electric field.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
c-----------------------------------------------------------------------
c       first find perturbed b, and j in real space.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,be,real_be,
     $               dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,real_ja,
     $               dealiase)
c-----------------------------------------------------------------------
c       find j_perturbed X b_perturbed / e in real space,
c-----------------------------------------------------------------------
        CALL math_cart_cross(real_eexp,real_ja,real_be,coefhll)
c-----------------------------------------------------------------------
c       include the -kboltz*Te*grad(n) term for finite pressure.  the
c       grad(Te) part does not contribute to dB/dt and is not included.
c       real_ve is used as temporary storage.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $                 real_ve,dealiase)
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,tele,
     $                 real_te,dealiase)
          real_eexp(1,:,:)=real_eexp(1,:,:)-
     $                       real_ve(1,:,:)*real_te(1,:,:)
          real_eexp(2,:,:)=real_eexp(2,:,:)-
     $                       real_ve(2,:,:)*real_te(1,:,:)
          real_eexp(3,:,:)=real_eexp(3,:,:)-
     $                       real_ve(3,:,:)*real_te(1,:,:)
        ENDIF
c-----------------------------------------------------------------------
c       transform to Fourier space.
c-----------------------------------------------------------------------
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,eexp,real_eexp,
     $               dealiase)
c-----------------------------------------------------------------------
c       create the nonlinear ideal mhd electric field.  these terms
c       are kept separate because they are not divided by n later.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $               dealiase)
        CALL math_cart_cross(real_eexp,real_be,real_ve,dt)
c-----------------------------------------------------------------------
c       add the explicit contribution to resistive dissipation if
c       3D resistivity is used.
c-----------------------------------------------------------------------
        IF (elecd>0.AND.threedeta) THEN
          DO ip=1,nphi
            DO jq=1,3
              real_eexp(jq,:,ip)=real_eexp(jq,:,ip)+
     $          dtm*( elecd_phi(1,:,ip)*real_ja(jq,:,ip)
     $              +(elecd_phi(1,:,ip)-eleq_1d(ipseust:ipseuen))
     $                *jaeq_1d(jq,ipseust:ipseuen))
            ENDDO
          ENDDO
        ENDIF
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,ebxv,real_eexp,
     $               dealiase)
      ELSE
        ebxv=0._r8
        eexp=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     begin the mode loop for linear contributions.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       find the explicit linear contributions to the electric field.
c-----------------------------------------------------------------------
        CALL math_cadd_cross(eexp(:,:,:,imode),
     $                       ja(:,:,:,imode),be_eq,coefhll)
        CALL math_cadd_cross(eexp(:,:,:,imode),
     $                       ja_eq,be(:,:,:,imode),coefhll)
c-----------------------------------------------------------------------
c       subtract the linear kbolz*grad(n)*Te terms.  the grad(Te)
c       should not contribute to dB/dt and is not included.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          eexp(1,:,:,imode)=eexp(1,:,:,imode)-
     $                   grad_nd(1,:,:,imode)*tele_eq(1,:,:)-
     $                 pfac*tele(1,:,:,imode)* nd_eqr(1,:,:)
          eexp(2,:,:,imode)=eexp(2,:,:,imode)-
     $                   grad_nd(2,:,:,imode)*tele_eq(1,:,:)-
     $                 pfac*tele(1,:,:,imode)* nd_eqz(1,:,:)
          eexp(3,:,:,imode)=eexp(3,:,:,imode)-
     $                   grad_nd(3,:,:,imode)*tele_eq(1,:,:)
        ENDIF
c-----------------------------------------------------------------------
c       find the linear contributions to the ideal mhd electric field.
c-----------------------------------------------------------------------
        CALL math_cadd_cross(ebxv(:,:,:,imode),be_eq,ve(:,:,:,imode),dt)
        IF (eq_flow/='none')
     $    CALL math_cadd_cross(ebxv(:,:,:,imode),
     $                         be(:,:,:,imode),ve_eq,dt)
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     add current advection terms ~ div(JV+VJ).  first
c     allocate and evaluate the needed fields, including the
c     finite element representation of J.
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all') THEN
        ALLOCATE(vep(3,ncx,ncy,nmodes),divv(1,ncx,ncy,nmodes),
     $           real_divv(1,mpseudo,nphi),jap(3,ncx,ncy,nmodes))
        CALL generic_ptr_set(rb%qja,tb%qja,tb%tgeom,inode,
     $                       jafe,jar,jaz,1_i4)
        CALL generic_ptr_set(rb%qdvv_eq,tb%qdvv_eq,tb%tgeom,inode,
     $                       divv_eq,dp,dp,0_i4)
        IF (eq_flow/='none') THEN
          CALL generic_ptr_set(rb%qgrdveq,tb%qgrdveq,tb%tgeom,
     $                         inode,grdv_eq,dp,dp,0_i4)
        ENDIF
c-----------------------------------------------------------------------
c       evaluate div(v) and the effective d(v)/dphi and d(j)/dphi.
c-----------------------------------------------------------------------
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            divv(1,:,:,imode)=ver(1,:,:,imode)+vez(2,:,:,imode)
     $                      +((0,1)*keff(imode)*ve(3,:,:,imode)
     $                                         +ve(1,:,:,imode))/bigr
          ELSE
            divv(1,:,:,imode)=ver(1,:,:,imode)+vez(2,:,:,imode)
     $                       +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr
          ENDIF
          vep(:,:,:,imode)=(0,1)*keff(imode)*ve(:,:,:,imode)
          jap(:,:,:,imode)=(0,1)*keff(imode)*jafe(:,:,:,imode)
          IF (geom=='tor') THEN
            vep(1,:,:,imode)=(vep(1,:,:,imode)-ve(3,:,:,imode))/bigr
            vep(2,:,:,imode)= vep(2,:,:,imode)/bigr
            vep(3,:,:,imode)=(vep(3,:,:,imode)+ve(1,:,:,imode))/bigr
            jap(1,:,:,imode)=(jap(1,:,:,imode)-jafe(3,:,:,imode))/bigr
            jap(2,:,:,imode)= jap(2,:,:,imode)/bigr
            jap(3,:,:,imode)=(jap(3,:,:,imode)+jafe(1,:,:,imode))/bigr
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       nonlinear products:
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,jafe,real_ja,
     $                 dealiase)
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,divv,
     $                 real_divv,dealiase)
c-----------------------------------------------------------------------
c         div(v)*j
c-----------------------------------------------------------------------
          DO iq=1,3
            real_eexp(iq,:,:)=real_divv(1,:,:)*real_ja(iq,:,:)
          ENDDO
c-----------------------------------------------------------------------
c         j.grad(v)--d(V)/dr, d(V)/dz and d(V)/dphi contributions.
c         pointers are used to save a little memory.
c-----------------------------------------------------------------------
          real_gradv=>real_be
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ver,
     $                 real_gradv,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ja(1,:,:)*real_gradv(iq,:,:)
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,vez,
     $                 real_gradv,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ja(2,:,:)*real_gradv(iq,:,:)
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,vep,
     $                 real_gradv,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ja(3,:,:)*real_gradv(iq,:,:)
          ENDDO
          NULLIFY(real_gradv)
c-----------------------------------------------------------------------
c         v.grad(j)--d(J)/dr, d(J)/dz and d(J)/dphi contributions.
c-----------------------------------------------------------------------
          real_gradj=>real_be
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,jar,
     $                 real_gradj,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ve(1,:,:)*real_gradj(iq,:,:)
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,jaz,
     $                 real_gradj,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ve(2,:,:)*real_gradj(iq,:,:)
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,jap,
     $                 real_gradj,dealiase)
          DO iq=1,3
            real_eexp(iq,:,:)=real_eexp(iq,:,:)
     $                +real_ve(3,:,:)*real_gradj(iq,:,:)
          ENDDO
          NULLIFY(real_gradj)
c-----------------------------------------------------------------------
c         transform the nonlinear contributions
c         (here, jaore is just a temporary array).
c-----------------------------------------------------------------------
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,jaore,
     $                 real_eexp,dealiase)
        ELSE
          jaore=0._r8
        ENDIF

        cur_adv: DO imode=1,nmodes
c-----------------------------------------------------------------------
c         div(v)*J_eq+div(V_eq)*j
c-----------------------------------------------------------------------
          DO jq=1,3
            jaore(jq,:,:,imode)=jaore(jq,:,:,imode)
     $                       +divv(1,:,:,imode)*ja_eq(jq,:,:)
     $                       +divv_eq(1,:,:)*jafe(jq,:,:,imode)
          ENDDO
c-----------------------------------------------------------------------
c         J_eq.grad(v)+j.grad(V_eq)
c-----------------------------------------------------------------------
          DO jq=1,3
            jaore(jq,:,:,imode)=jaore(jq,:,:,imode)
     $                         +ja_eq(1,:,:)*ver(jq,:,:,imode)
     $                         +ja_eq(2,:,:)*vez(jq,:,:,imode)
     $                         +ja_eq(3,:,:)*vep(jq,:,:,imode)
          ENDDO
          IF (eq_flow/='none') THEN
            jaore(1,:,:,imode)=jaore(1,:,:,imode)+
     $                          jafe(1,:,:,imode)*grdv_eq(1,:,:)+
     $                          jafe(2,:,:,imode)*grdv_eq(2,:,:)+
     $                          jafe(3,:,:,imode)*grdv_eq(3,:,:) 
            jaore(2,:,:,imode)=jaore(2,:,:,imode)+
     $                          jafe(1,:,:,imode)*grdv_eq(4,:,:)+
     $                          jafe(2,:,:,imode)*grdv_eq(5,:,:) 
            jaore(3,:,:,imode)=jaore(3,:,:,imode)+
     $                          jafe(1,:,:,imode)*grdv_eq(7,:,:)+
     $                          jafe(2,:,:,imode)*grdv_eq(8,:,:)+
     $                          jafe(3,:,:,imode)*grdv_eq(9,:,:) 
          ENDIF
c-----------------------------------------------------------------------
c         V_eq.grad(j)+v.grad(J_eq)
c-----------------------------------------------------------------------
          IF (eq_flow/='none') THEN
            DO jq=1,3
              jaore(jq,:,:,imode)=jaore(jq,:,:,imode)
     $                           +ve_eq(1,:,:)*jar(jq,:,:,imode)
     $                           +ve_eq(2,:,:)*jaz(jq,:,:,imode)
     $                           +ve_eq(3,:,:)*jap(jq,:,:,imode)
            ENDDO
          ENDIF
          DO jq=1,3
            jaore(jq,:,:,imode)=jaore(jq,:,:,imode)
     $                         +ve(1,:,:,imode)*jaeq_r(jq,:,:)
     $                         +ve(2,:,:,imode)*jaeq_z(jq,:,:)
          ENDDO
          IF (geom=='tor') THEN
            jaore(1,:,:,imode)=jaore(1,:,:,imode)
     $                        -ve(3,:,:,imode)*ja_eq(3,:,:)/bigr
            jaore(3,:,:,imode)=jaore(2,:,:,imode)
     $                        +ve(3,:,:,imode)*ja_eq(1,:,:)/bigr
          ENDIF
c-----------------------------------------------------------------------
c         multiply by me/e**2 and add to other terms.
c-----------------------------------------------------------------------
          DO jq=1,3
            eexp(jq,:,:,imode)=coefme2*jaore(jq,:,:,imode)
     $                        +eexp(jq,:,:,imode)
          ENDDO
        ENDDO cur_adv
      ENDIF
c-----------------------------------------------------------------------
c     all of the preceeding terms must be divided by number density.
c     the n that is used depends on the continuity input parameter.
c-----------------------------------------------------------------------
      IF (continuity=='none'.OR.
     $    nonlinear.AND.continuity=='fix profile') THEN
        DO imode=1,nmodes
          DO jq=1,3
            eexp(jq,:,:,imode)=dt*eexp(jq,:,:,imode)/nd_eq(1,:,:)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     when a perturbed n is used in the denominator, there is also
c     a term associated with equilibrium quantities in the numerator.
c     it is now computed once and saved.
c-----------------------------------------------------------------------
      ELSE
        CALL generic_ptr_set(rb%qeq_force,tb%qeq_force,tb%tgeom,
     $                       inode,eq_force,dp,dp,0_i4)
c-----------------------------------------------------------------------
c       linear density perturbation.
c-----------------------------------------------------------------------
        IF (.NOT.nonlinear) THEN
          DO imode=1,nmodes
            DO jq=1,3
              eexp(jq,:,:,imode)=dt/nd_eq(1,:,:)*
     $          (eexp(jq,:,:,imode)-
     $           eq_force(jq,:,:)*nd(1,:,:,imode)/nd_eq(1,:,:))
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c       use symmetric part of evolving density.
c-----------------------------------------------------------------------
        ELSE IF (continuity=='n=0 only') THEN
          CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                         nd_n0,dp,dp,0_i4)
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              DO jq=1,3
                eexp(jq,:,:,imode)=dt*(eexp(jq,:,:,imode)/
     $                             (nd_eq(1,:,:)+nd_n0(1,:,:))
     $              +eq_force(jq,:,:)*(1._r8/(nd_eq(1,:,:)+nd_n0(1,:,:))
     $                                -1._r8/ nd_eq(1,:,:)))
              ENDDO
            ELSE
              DO jq=1,3
                eexp(jq,:,:,imode)=dt*eexp(jq,:,:,imode)/
     $            (nd_eq(1,:,:)+nd_n0(1,:,:))
              ENDDO
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c       use 3D evolving density.
c-----------------------------------------------------------------------
        ELSE
          CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                         inode,real_ndptr,dp,dp,0_i4)
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              DO jq=1,3
                eexp(jq,:,:,imode)=eexp(jq,:,:,imode)+eq_force(jq,:,:)
              ENDDO
            ENDIF
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,
     $                 eexp,real_eexp,dealiase)
          DO jq=1,3
            real_eexp(jq,:,:)=dt*real_eexp(jq,:,:)/real_ndptr(1,:,:)
          ENDDO
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,
     $                 eexp,real_eexp,dealiase)
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              DO jq=1,3
                eexp(jq,:,:,imode)=eexp(jq,:,:,imode)-
     $                             dt*eq_force(jq,:,:)/nd_eq(1,:,:)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      IF (ohms=='2fl'.AND.advect=='all')
     $  DEALLOCATE(vep,divv,real_divv,jap)
c-----------------------------------------------------------------------
c     add the explicit contribution to resistive dissipation, where
c     the electrical diffusivity is 2D.
c-----------------------------------------------------------------------
      IF (elecd>0) THEN
        IF (.NOT.(eta_model=='eta full'.OR.eta_model=='chodura')) THEN
          DO imode=1,nmodes
            eexp(1,:,:,imode)=eexp(1,:,:,imode)+diff*ja(1,:,:,imode)
            eexp(2,:,:,imode)=eexp(2,:,:,imode)+diff*ja(2,:,:,imode)
            eexp(3,:,:,imode)=eexp(3,:,:,imode)+diff*ja(3,:,:,imode)
          ENDDO
          IF (eta_model=='eta n=0 only'.AND.nonlinear) THEN
            DO imode=1,nmodes
              IF (keff(imode)/=0) CYCLE
              DO jq=1,3
                eexp(jq,:,:,imode)=eexp(jq,:,:,imode)+
     $            dtm*(elecd_n0(1,:,:)-elecd_eq(1,:,:))*ja_eq(jq,:,:)
              ENDDO
            ENDDO
          ENDIF
        ELSE IF (.NOT.nonlinear.AND.eta_model=='eta full') THEN
          DO imode=1,nmodes
            DO jq=1,3
              eexp(jq,:,:,imode)=eexp(jq,:,:,imode)+
     $          dtm*elecd_eq(1,:,:)*(ja(jq,:,:,imode)-
     $             1.5*ja_eq(jq,:,:)*tele(1,:,:,imode)/tele_eq(1,:,:))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     construct -curl(alpha).E_explicit [ -div(alpha)*div(B) ].
c-----------------------------------------------------------------------
      eexp=eexp+ebxv
      IF (divbd>0.AND..NOT.split_divb) THEN
        DO imode=1,nmodes
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM( dalpdz(:,:,iv)*eexp(3,:,:,imode)
     $         +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(2,:,:,imode)
     $             -dbdt*dalpdrc(:,:,iv)*divb(:,:,imode),1)
            int(2,:,iv,imode)=SUM(-dalpdr(:,:,iv)*eexp(3,:,:,imode)
     $         -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(1,:,:,imode)
     $             -dbdt*dalpdz (:,:,iv)*divb(:,:,imode),1)
            int(3,:,iv,imode)=SUM(-dalpdz (:,:,iv)*eexp(1,:,:,imode)
     $                          +dalpdrc(:,:,iv)*eexp(2,:,:,imode)
     $                      +(0,1)*alpha(:,:,iv)*keff(imode)/bigr
     $                             *dbdt*divb(:,:,imode),1)
          ENDDO
        ENDDO
      ELSE
        DO imode=1,nmodes
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM( dalpdz(:,:,iv)*eexp(3,:,:,imode)
     $        +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(2,:,:,imode),1)
            int(2,:,iv,imode)=SUM(-dalpdr(:,:,iv)*eexp(3,:,:,imode)
     $        -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(1,:,:,imode),1)
            int(3,:,iv,imode)=SUM(-dalpdz (:,:,iv)*eexp(1,:,:,imode)
     $                            +dalpdrc(:,:,iv)*eexp(2,:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if hyper-resistivity is used without a split step, an auxiliary
c     equation is solved simultaneously with dB.  the integrand for the
c     right side of this equation is
c
c       -SQRT(dt*hyp_eta/fhyp_eta)*curl(C^*).curl(B_old)
c       -SQRT(dt*hyp_dbd/fhyp_dbd)*div(C^*).div(B_old)
c
c     where C is the test vector for the auxiliary field.  note that
c     ja is modified at this point.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0.OR.hyp_dbd>0).AND..NOT.split_hypeta) THEN
        fhyp=SQRT(dt*hyp_eta/fhyp_eta)*mu0
        fhdb=SQRT(dt*hyp_dbd/fhyp_dbd)
        DO imode=1,nmodes
          ja(:,:,:,imode)=ja(:,:,:,imode)*fhyp
          dbqd=fhdb*divb(:,:,imode)
          DO iv=1,nvc
            int(4,:,iv,imode)=
     $        SUM( dalpdz(:,:,iv)*ja(3,:,:,imode)+
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*ja(2,:,:,imode)-
     $             dalpdrc(:,:,iv)*dbqd,1)
            int(5,:,iv,imode)=
     $        SUM(-dalpdr(:,:,iv)*ja(3,:,:,imode)-
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*ja(1,:,:,imode)-
     $             dalpdz(:,:,iv)*dbqd,1)
            int(6,:,iv,imode)=
     $        SUM(-dalpdz (:,:,iv)*ja(1,:,:,imode)+
     $             dalpdrc(:,:,iv)*ja(2,:,:,imode)+
     $             (0,1)*alpha(:,:,iv)*keff(imode)/bigr*dbqd,1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the integrand for the discontinuous scalar equation,
c     disdbfac*aldis*div(b), where aldis is the value of the
c     discontinuous basis function.
c-----------------------------------------------------------------------
      IF (poly_divb>=0) THEN
        disdbfac=SQRT(disc_dbd*dt/fdivb)
        DO imode=1,nmodes
          DO iv=nvc+1,nv
            int(1,:,iv,imode)=
     $        disdbfac*SUM(aldis(:,:,iv-nvc)*divb(:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE brhs_hmhd
c-----------------------------------------------------------------------
c     subprogram 6. hall_cor.
c     find the correction term for the rhs of the bhmhd advance
c     when nonlinear iteration is used for dJxdB.
c
c     here, db is all of the change from the start of the time-step, and
c     because it's also used in the dot routine, the correction to the
c     rhs has the opposite sign of the explicit Hall term, i.e.
c
c       +0.25*dt*curl(test vector).(dJXdB)/ne
c-----------------------------------------------------------------------
      SUBROUTINE hall_cor(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          nd_eq,nd_n0,real_ndptr,dp
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ja,real_be,real_eexp
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: hfac

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,ber,bez,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             eexp,ja
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,iq,nvc,ncx,ncy
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      hfac=-0.25_r8*dt*coefhll
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the 0-th
c     order magnetic field.  evaluate the old current density
c     from the gradients of the old magnetic fields for resitive
c     dissipation.  for the mhd part, evaluate flow velocities.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      CALL generic_ptr_set(rb%qwork1,tb%qwork1,tb%tgeom,inode,
     $                     be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8/mu0)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find the explicit nonlinear contributions to the electric field.
c
c     first find perturbed b, and j in real space.
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,be,real_be,
     $             dealiase)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,real_ja,
     $             dealiase)
c-----------------------------------------------------------------------
c     find j_perturbed X b_perturbed / e in real space,
c-----------------------------------------------------------------------
      CALL math_cart_cross(real_eexp,real_ja,real_be,hfac)
c-----------------------------------------------------------------------
c     transform to Fourier space.
c-----------------------------------------------------------------------
      IF (continuity/='full') THEN
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,eexp,real_eexp,
     $               dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     all of the preceeding terms must be divided by number density.
c     the n that is used depends on the continuity input parameter.
c-----------------------------------------------------------------------
      IF (continuity=='none'.OR.continuity=='fix profile') THEN
        DO imode=1,nmodes
          eexp(1,:,:,imode)=eexp(1,:,:,imode)/nd_eq(1,:,:)
          eexp(2,:,:,imode)=eexp(2,:,:,imode)/nd_eq(1,:,:)
          eexp(3,:,:,imode)=eexp(3,:,:,imode)/nd_eq(1,:,:)
        ENDDO
c-----------------------------------------------------------------------
c       use symmetric part of evolving density.
c-----------------------------------------------------------------------
      ELSE IF (continuity=='n=0 only') THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                       nd_n0,dp,dp,0_i4)
        DO imode=1,nmodes
          eexp(1,:,:,imode)=eexp(1,:,:,imode)/
     $        (nd_eq(1,:,:)+nd_n0(1,:,:))
          eexp(2,:,:,imode)=eexp(2,:,:,imode)/
     $        (nd_eq(1,:,:)+nd_n0(1,:,:))
          eexp(3,:,:,imode)=eexp(3,:,:,imode)/
     $        (nd_eq(1,:,:)+nd_n0(1,:,:))
        ENDDO
c-----------------------------------------------------------------------
c     use 3D evolving density.
c-----------------------------------------------------------------------
      ELSE
        CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)
        real_eexp(1,:,:)=real_eexp(1,:,:)/real_ndptr(1,:,:)
        real_eexp(2,:,:)=real_eexp(2,:,:)/real_ndptr(1,:,:)
        real_eexp(3,:,:)=real_eexp(3,:,:)/real_ndptr(1,:,:)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,
     $               eexp,real_eexp,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     construct -curl(alpha).delta(E)
c
c     if hyper-resitivity is used without splitting, ensure that the
c     rows for the auxiliary field are 0.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        DO iv=1,nvc
          int(1,:,iv,imode)=SUM( dalpdz(:,:,iv)*eexp(3,:,:,imode)
     $       +(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(2,:,:,imode),1)
          int(2,:,iv,imode)=SUM(-dalpdr(:,:,iv)*eexp(3,:,:,imode)
     $       -(0,1)*alpha(:,:,iv)*keff(imode)/bigr*eexp(1,:,:,imode),1)
          int(3,:,iv,imode)=SUM(-dalpdz (:,:,iv)*eexp(1,:,:,imode)
     $                          +dalpdrc(:,:,iv)*eexp(2,:,:,imode),1)
        ENDDO
        IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
          DO iv=1,nvc
            int(4:6,:,iv,imode)=0._r8
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE hall_cor
c-----------------------------------------------------------------------
c     subprogram 7. tirhs.
c     compute the integrand used for the rhs of the ion temperature (or
c     temperature, if separate_pe=F) change equation with implicit
c     diffusion.
c
c     dt*(gamma-1)*(D(Ti_n)-ni*Ti*div(Vi))-dt*ni*Vi.grad(Ti)
c
c     where Ti_n is the temperature at the beginning of the time step,
c     Ti is either Ti_n or the predicted temperature, and D is the
c     diffusion operator.
c-----------------------------------------------------------------------
      SUBROUTINE tirhs(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          ti_eq,ti_eqr,ti_eqz,
     $          ve_eq,be_eq,ja_eq,b0,q_a,nd_eq,nd_eqr,
     $          nd_eqz,nd_n0,real_bptr,real_ndptr,elecd_t,elecd_eq,dp,
     $          ds2,kappli,kaprpi,ti_b2,bcrgti,divv_eq,real_vptr,
     $          real_jptr,upwc,real_ndiff
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_work,bbgradt
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal,real_divv
      REAL(r8), DIMENSION(mpseudo,nphi) :: magBsq3D
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: eleq_1d,elt_1d
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: tqr_1d,tqz_1d
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2)) :: jaeq_1d 
      REAL(r8), DIMENSION(3,mpseudo) :: rvec
      REAL(r8), DIMENSION(  mpseudo) :: jaeq2,rdot
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kdiff,pifrac
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             ve,dvr,dvz,be,ber,bez,ndiff,
     $             ti,tir,tiz,nd,ndr,ndz,vheat,hypheat,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             grad_ti,be_tot,grad_nd,vel,ja
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             vdgr_ti,divv,tdivv,bdgr_ti,ohm_en
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: b_cr_gt
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: cr_coef
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,imode,iq,jq,ncx,ncy,ip
      REAL(r8) :: fupw
      LOGICAL :: told_set
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      told_set=.false.
      IF (separate_pe) THEN
        pifrac=1._r8
      ELSE
        pifrac=1._r8-pe_frac
      ENDIF
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate velocity.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     temperature and grad(T) at the beginning of this time split for
c     predictor steps, or at the centered prediction for corrector
c     steps.  we assume ja_eq/=0, so if separate_pe=T, this integrand
c     is used with predictor/corrector steps.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qtion_eq,tb%qtion_eq,tb%tgeom,
     $                     inode,ti_eq,ti_eqr,ti_eqz,1_i4)
      IF (nonlinear.OR.eq_flow/='none'.OR.separate_pe) THEN
        IF (integrand_flag=='all ti terms pre'.OR.
     $      separate_pe.AND.k_cross>0.OR.impladv) THEN
          CALL generic_ptr_set(rb%qtion,tb%qtion,tb%tgeom,inode,
     $                         ti,tir,tiz,1_i4)
          told_set=.true.
        ELSE
          CALL generic_ptr_set(rb%qwork3,tb%qwork3,tb%tgeom,inode,
     $                         ti,tir,tiz,1_i4)
        ENDIF
        DO imode=1,nmodes
          grad_ti(1,:,:,imode)=tir(1,:,:,imode)
          grad_ti(2,:,:,imode)=tiz(1,:,:,imode)
          grad_ti(3,:,:,imode)=(0,1)*keff(imode)*ti(1,:,:,imode)/bigr
        ENDDO
      ENDIF
      IF (visc_heat) THEN
        CALL generic_ptr_set(rb%qvisc,tb%qvisc,tb%tgeom,inode,vheat,
     $                       dcp,dcp,0_i4)
      ENDIF
      IF (nonlinear.AND.ohm_heat.AND.(hyp_eta>0.OR.hyp_dbd>0)) THEN
        CALL generic_ptr_set(rb%qhyph,tb%qhyph,tb%tgeom,inode,hypheat,
     $                       dcp,dcp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate ion (if separate_pe=T) or COM velocity.
c
c     The ion velocity is computed from the COM velocity
c     and the total current density as
c
c          vel = v + j * nu / ( n e (nu+1) )
c
c     where "n" is the electron number density, "e" the electron
c     charge and "nu" is the electron-ion mass ratio times the atomic
c     number "Z":
c                  nu = Z me / mi
c
c     we also need div(v)=d(Vr)/dR + Vr/R + d(Vz)dZ + 1/R d(Vth)/dTh
c     here we assume that Jeq.grad(nd_eq)=0.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,dvr,dvz,1_i4)
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,real_vptr,dp,dp,0_i4)
      ENDIF
      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,nd,ndr,ndz,1_i4)
      IF (nonlinear.AND.continuity=='n=0 only')
     $  CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,nd_eqr,nd_eqz,1_i4)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,
     $                     be,ber,bez,1_i4)
      IF (separate_pe.OR.ohm_heat) THEN
        CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8/mu0)
        CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                       inode,ja_eq,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     extra pointers are needed for the upwinding-like artificial
c     diffusion.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        CALL generic_ptr_set(rb%qupti_phi,tb%qupti_phi,tb%tgeom,
     $                       inode,upwc,dp,dp,0_i4)
        IF (separate_pe)
     $    CALL generic_ptr_set(rb%qja_tot,tb%qja_tot,tb%tgeom,
     $                         inode,real_jptr,dp,dp,0_i4)
        fupw=t_dart_upw*dt**2
        tqr_1d=RESHAPE(ti_eqr,(/ncx*ncy/))
        tqz_1d=RESHAPE(ti_eqz,(/ncx*ncy/))
      ENDIF
c-----------------------------------------------------------------------
c     div(v-COM) and grad(n)
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (geom=='tor') THEN
          divv(1,:,:,imode)=dvr(1,:,:,imode)+dvz(2,:,:,imode)
     $                    +((0,1)*keff(imode)*ve(3,:,:,imode)
     $                                       +ve(1,:,:,imode))/bigr
        ELSE
          divv(1,:,:,imode)=dvr(1,:,:,imode)+dvz(2,:,:,imode)
     $                     +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr
        ENDIF
        grad_nd(1,:,:,imode)=ndr(1,:,:,imode)
        grad_nd(2,:,:,imode)=ndz(1,:,:,imode)
        grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
      ENDDO

      IF (nonlinear)
     $  CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)

      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qdvv_eq,tb%qdvv_eq,tb%tgeom,inode,
     $                       divv_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find the Ohmic heating energy before the ja array is modified.
c     only 1-pe_frac of it will be applied in the ion energy equation
c     when separate_pe=false; otherwise, none of it is applied.
c-----------------------------------------------------------------------
      ohm_en=0._r8
      IF (nonlinear.AND.(separate_pe.OR.ohm_heat))
     $  CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,real_ve,
     $               dealiase)
      IF (ohm_heat.AND..NOT.separate_pe) THEN
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_eq,tb%qelecd_eq,tb%tgeom,
     $                         inode,elecd_eq,dp,dp,0_i4)
          eleq_1d=RESHAPE(elecd_eq,(/ncx*ncy/))
          jaeq_1d=RESHAPE( ja_eq,(/3,ncx*ncy/))
          DO iq=1,mpseudo
            jaeq2(iq)=SUM(jaeq_1d(:,iq-1+ipseust)**2)
          ENDDO
        ENDIF
        DO imode=1,nmodes
          ohm_en(1,:,:,imode)=2*SUM(ja_eq*ja(:,:,:,imode),1)
        ENDDO
        IF (nonlinear) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,
     $                 ohm_en,real_scal,dealiase)
          real_scal(1,:,:)=real_scal(1,:,:)+SUM(real_ve**2,1)
          IF (threedeta) THEN
            CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            DO ip=1,nphi
              DO iq=1,mpseudo
                real_scal(1,iq,ip)=mu0*(1._r8-pe_frac)*
     $            (real_scal(1,iq,ip)*elecd_t(1,iq,ip)+
     $            (elecd_t(1,iq,ip)-eleq_1d(iq-1+ipseust))*jaeq2(iq))
              ENDDO
            ENDDO
          ELSE IF (eta_model=='eta n=0 only') THEN
            CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            elt_1d=RESHAPE(elecd_t,(/ncx*ncy/))
            DO iq=1,mpseudo
              jaeq2(iq)=jaeq2(iq)*
     $              (elt_1d(iq-1+ipseust)-eleq_1d(iq-1+ipseust))
            ENDDO
            DO ip=1,nphi
              DO iq=1,mpseudo
                real_scal(1,iq,ip)=mu0*(1._r8-pe_frac)*
     $            (real_scal(1,iq,ip)*elt_1d(iq-1+ipseust)+jaeq2(iq))
              ENDDO
            ENDDO
          ENDIF
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                 ohm_en,real_scal,dealiase)
        ENDIF
        IF (nonlinear.AND.threedeta) THEN  ! done
        ELSE IF (eta_model=='eta n=0 only') THEN  ! done if nonlinear
          IF (.NOT.nonlinear) THEN
            CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            DO imode=1,nmodes
              ohm_en(:,:,:,imode)=ohm_en(:,:,:,imode)*elecd_t*mu0*
     $          (1._r8-pe_frac)
            ENDDO
          ENDIF
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          DO imode=1,nmodes
            ohm_en(:,:,:,imode)=
     $        ohm_en(:,:,:,imode)*elecd*MAX(ds2,0._r8)*mu0*
     $          (1._r8-pe_frac)
          ENDDO
        ELSE
          ohm_en=ohm_en*elecd*mu0*(1._r8-pe_frac)
        ENDIF
        IF (nonlinear.AND.(hyp_eta>0.OR.hyp_dbd>0))
     $    ohm_en=ohm_en+hypheat*(1._r8-pe_frac)
      ENDIF
c-----------------------------------------------------------------------
c     if V_i differs from V, convert v and div(v) to ion velocity.  the
c     O(me/mi) terms are retained in case me/mi is modified.
c     j.grad(n) (saved in vdgr_ti) for nonlinear, full continuity
c     computations.  [keep div(j)=0 in mind here]
c-----------------------------------------------------------------------
      IF (coefjvi/=0._r8) THEN
        IF (nonlinear) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $                 real_work,dealiase)
          real_scal(1,:,:)=SUM(real_ve*real_work,1)
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $                 real_scal,dealiase)
        ELSE
          vdgr_ti=0._r8
        ENDIF

        IF (.NOT.nonlinear.OR.continuity=='fix profile') THEN
          DO imode=1,nmodes
            DO jq=1,3
              vel(jq,:,:,imode)=ve(jq,:,:,imode)
     $          +ja(jq,:,:,imode)*coefjvi/nd_eq(1,:,:)
            ENDDO
            divv(1,:,:,imode)=divv(1,:,:,imode)-
     $        coefjvi*(SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $                 +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $                 +ja(2,:,:,imode)*nd_eqz(1,:,:)
     $                 +vdgr_ti(1,:,:,imode))/(nd_eq(1,:,:)**2)
          ENDDO
        ELSE IF (continuity=='n=0 only') THEN
          bmagsq=1._r8/(nd_eq(1,:,:)+nd_n0(1,:,:))
          DO imode=1,nmodes
            DO jq=1,3
              vel(jq,:,:,imode)=ve(jq,:,:,imode)
     $          +ja(jq,:,:,imode)*coefjvi*bmagsq
            ENDDO
            divv(1,:,:,imode)=divv(1,:,:,imode)-
     $        coefjvi*(SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $                 +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $                 +ja(2,:,:,imode)*nd_eqz(1,:,:)
     $                 +vdgr_ti(1,:,:,imode))*bmagsq**2
          ENDDO
        ELSE  !  continuity is full  (ja will be left as ja/n).
          DO imode=1,nmodes
            vdgr_ti(1,:,:,imode)=-coefjvi*(vdgr_ti(1,:,:,imode)+
     $        SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $        +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $        +ja(2,:,:,imode)*nd_eqz(1,:,:))
          ENDDO
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $                 real_scal,dealiase)
          DO jq=1,3
            real_ve(jq,:,:)=real_ve(jq,:,:)/real_ndptr(1,:,:)
          ENDDO
          real_scal=real_scal/real_ndptr**2
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,ja,
     $                 real_ve,dealiase)
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $                 real_scal,dealiase)
          vel=ve+coefjvi*ja
          divv=divv+vdgr_ti
        ENDIF
      ELSE
        vel=ve
      ENDIF
c-----------------------------------------------------------------------
c     find the linear v_tilde.grad(Teq) for each mode--needed here
c     to complete the upwinding-like diffusion at the next step.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        vdgr_ti(1,:,:,imode)=vel(1,:,:,imode)*ti_eqr(1,:,:)+
     $                       vel(2,:,:,imode)*ti_eqz(1,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     find the energy density corrections for artificial particle
c     diffusivity and add to ohm_en.
c-----------------------------------------------------------------------
      IF (nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        CALL generic_ptr_set(rb%qndiffa,tb%qndiffa,tb%tgeom,inode,
     $                       ndiff,dcp,dcp,0_i4)
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qndiff_phi,tb%qndiff_phi,tb%tgeom,
     $                         inode,real_ndiff,dp,dp,0_i4)
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,ti,
     $                 real_scal,dealiase)
          IF (advect/='none') THEN
            real_divv(1,:,:)=real_ndiff(1,:,:)*
     $        (0.5_r8*mtot*SUM(real_vptr(:,:,:)**2,1)*pifrac-
     $         kboltz*real_scal(1,:,:)/(zeff*gamm1))
          ELSE
            real_divv(1,:,:)=-real_ndiff(1,:,:)*
     $        kboltz*real_scal(1,:,:)/(zeff*gamm1)
          ENDIF
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,bdgr_ti,
     $                 real_divv,dealiase)
          ohm_en=ohm_en+bdgr_ti
        ELSE IF (eq_flow/='none'.AND.advect/='none') THEN
          DO imode=1,nmodes
            ohm_en(1,:,:,imode)=ohm_en(1,:,:,imode)+
     $        ndiff(1,:,:,imode)*mtot*SUM(ve_eq(:,:,:)**2,1)*pifrac
          ENDDO
        ENDIF
        DO imode=1,nmodes
          ohm_en(1,:,:,imode)=ohm_en(1,:,:,imode)-
     $      ndiff(1,:,:,imode)*kboltz*ti_eq(1,:,:)/(zeff*gamm1)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the nonlinear contributions to v.grad(t) and t div(v)
c
c     transform v_perturbed and grad(t_perturbed) to real space,
c     find v.grad(t), then transform to Fourier space.  tdivv is
c     just used as temporary storage.
c
c     the upwinding-like diffusion term is also computed here.  the
c     predetermined coefficient is multiplied by VV.grad(T).
c     ** NOTE ** real_ve is now used to save this coefficient from
c     this point, so V_tilde(phi) array is lost.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,vel,real_ve,
     $               dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_ti,
     $               real_work,dealiase)
        real_scal(1,:,:)=SUM(real_ve*real_work,1)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tdivv,
     $               real_scal,dealiase)

        IF (impladv.AND.t_dart_upw>0) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $                 real_scal,dealiase)
          IF (coefjvi/=0._r8) THEN
            DO ip=1,nphi
              rvec(1,:)=real_vptr(1,:,ip)+
     $             coefjvi*real_jptr(1,:,ip)/real_ndptr(1,:,ip)
              rvec(2,:)=real_vptr(2,:,ip)+
     $             coefjvi*real_jptr(2,:,ip)/real_ndptr(1,:,ip)
              rvec(3,:)=real_vptr(3,:,ip)+
     $             coefjvi*real_jptr(3,:,ip)/real_ndptr(1,:,ip)
              rdot=fupw*upwc(1,:,ip)*upw_aniso*
     $          ( real_scal(1,:,ip) + SUM(rvec*real_work(:,:,ip),1) ) /
     $          ( SUM(rvec*rvec,1) + smallnum )
              jaeq2=fupw*upwc(1,:,ip)*(1._r8-upw_aniso)
              real_ve(1,:,ip)=rvec(1,:)*rdot+
     $          (real_work(1,:,ip)+tqr_1d(ipseust:ipseuen))*jaeq2
              real_ve(2,:,ip)=rvec(2,:)*rdot+
     $          (real_work(2,:,ip)+tqz_1d(ipseust:ipseuen))*jaeq2
              real_ve(3,:,ip)=rvec(3,:)*rdot+real_work(3,:,ip)*jaeq2
            ENDDO
          ELSE
            DO ip=1,nphi
              rdot=fupw*upwc(1,:,ip)*upw_aniso * ( real_scal(1,:,ip) +
     $            SUM(real_vptr(:,:,ip)*real_work(:,:,ip),1) ) /
     $          ( SUM(real_vptr(:,:,ip)*real_vptr(:,:,ip),1) + smallnum)
              jaeq2=fupw*upwc(1,:,ip)*(1._r8-upw_aniso)
              real_ve(1,:,ip)=real_vptr(1,:,ip)*rdot+
     $          (real_work(1,:,ip)+tqr_1d(ipseust:ipseuen))*jaeq2
              real_ve(2,:,ip)=real_vptr(2,:,ip)*rdot+
     $          (real_work(2,:,ip)+tqz_1d(ipseust:ipseuen))*jaeq2
              real_ve(3,:,ip)=real_vptr(3,:,ip)*rdot+
     $           real_work(3,:,ip)*jaeq2
            ENDDO
          ENDIF
        ENDIF
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     subtract the rest of the linear V.grad(T) for each mode.
c-----------------------------------------------------------------------
      IF (eq_flow/='none') THEN
        DO imode=1,nmodes
          vdgr_ti(1,:,:,imode)=-dt*( vdgr_ti(1,:,:,imode)+
     $      tdivv(1,:,:,imode)+SUM(ve_eq*grad_ti(:,:,:,imode),1) )
        ENDDO
      ELSE
        vdgr_ti=-dt*(vdgr_ti+tdivv)
      ENDIF
      IF (coefjvi/=0._r8) THEN
        DO imode=1,nmodes
          vdgr_ti(1,:,:,imode)=vdgr_ti(1,:,:,imode)-
     $       dt*coefjvi*SUM(ja_eq*grad_ti(:,:,:,imode),1)/nd_eq(1,:,:)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     multiply t_perturbed by div(v), then transform the
c     product to Fourier space.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.p_model/='isothermal') THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,ti,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,divv,
     $               real_divv,dealiase)
        real_scal=real_scal*real_divv
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tdivv,
     $               real_scal,dealiase)
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     subtract the linear (gamma-1)*T*div(V) for this mode.
c     also form the n_tilde*Teq*div(Veq) term and save in tdivv.
c-----------------------------------------------------------------------
      IF (p_model/='isothermal') THEN
        IF (eq_flow/='none') THEN
          DO imode=1,nmodes
            vdgr_ti(1,:,:,imode)=vdgr_ti(1,:,:,imode)
     $                         -dt*gamm1*(tdivv(1,:,:,imode)
     $                              +ti_eq(1,:,:)*divv(1,:,:,imode)
     $                              +ti(1,:,:,imode)*divv_eq(1,:,:))
            tdivv(1,:,:,imode)=-dt*gamm1*nd(1,:,:,imode)*
     $                               ti_eq(1,:,:)*divv_eq(1,:,:)/zeff
          ENDDO
        ELSE
          DO imode=1,nmodes
            vdgr_ti(1,:,:,imode)=vdgr_ti(1,:,:,imode)
     $                         -dt*gamm1*(tdivv(1,:,:,imode)
     $                              +ti_eq(1,:,:)*divv(1,:,:,imode))
            tdivv(1,:,:,imode)=0._r8
          ENDDO
        ENDIF
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     set temperature and grad(T) with T from the beginning of this
c     time split for thermal conduction if not already done.  
c-----------------------------------------------------------------------
      IF (p_model=='adiabat'.OR.p_model=='isothermal') THEN
        grad_ti=0._r8
      ELSE
        IF (.NOT.told_set) THEN
          CALL generic_ptr_set(rb%qtion,tb%qtion,tb%tgeom,inode,
     $                         ti,tir,tiz,1_i4)
          DO imode=1,nmodes
            grad_ti(1,:,:,imode)=tir(1,:,:,imode)
            grad_ti(2,:,:,imode)=tiz(1,:,:,imode)
            grad_ti(3,:,:,imode)=(0,1)*keff(imode)*ti(1,:,:,imode)/bigr
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c-PRE Heat soure term.
c     q_applied is set in q_applied_init.f and is normally zero
c-----------------------------------------------------------------------
c     DO imode=1,nmodes
c       IF (keff(imode)==0) THEN
c         CALL generic_ptr_set(rb%qq_applied,tb%qq_applied
c    $                         tb%tgeom,inode,q_a,dp,dp,0_i4)
c         vdgr_ti(1,:,:,imode)=vdgr_ti(1,:,:,imode)
c    $                         +dt*gamm1*q_a/(kboltz*bigr)
c       ENDIF
c     ENDDO
c-----------------------------------------------------------------------
c     For isotropic diffusion, multiply grad_ti by dt*diffusivity.
c     assume k_perp has gamma-1 factor in it.
c-----------------------------------------------------------------------
      pmodelif: IF (p_model=='isotropic') THEN
        grad_ti=grad_ti*dt*k_perpi
      ELSE IF (p_model(1:5)=='aniso') THEN
c-----------------------------------------------------------------------
c     for anisotropic thermal conduction, get the
c     magnetic field first.  note that we assume k_pll and k_perp have
c     the factor of (gamma-1) in them already.  for single-fluid
c     temperature, use the dominant k_plle.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      IF (nonlinear.AND.(p_model=="aniso_plltdep".OR.
     $                   p_model=="aniso_tdep")) THEN
        CALL generic_ptr_set(rb%qkappli_phi,tb%qkappli_phi,tb%tgeom,
     $                       inode,kappli,dp,dp,0_i4)

        IF (.NOT.closure_n0_only) THEN
          CALL generic_ptr_set(rb%qkaprpi_phi,tb%qkaprpi_phi,tb%tgeom,
     $                         inode,kaprpi,dp,dp,0_i4)
        ELSE
          CALL generic_ptr_set(rb%qkaprpi_n0,tb%qkaprpi_n0,tb%tgeom,
     $                         inode,kaprpi,dp,dp,0_i4)
          elt_1d=RESHAPE(kaprpi,(/ncx*ncy/))
        ENDIF
      ELSE
        IF (separate_pe) THEN
          kdiff=k_plli-k_perpi
        ELSE
          kdiff=k_plle-k_perpi
        ENDIF
      ENDIF
      IF (separate_pe.AND.k_cross>0) THEN
        CALL generic_ptr_set(rb%qti_b2,tb%qti_b2,tb%tgeom,
     $                       inode,ti_b2,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qbcrgti,tb%qbcrgti,tb%tgeom,
     $                       inode,bcrgti,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     For nonlinear calculations, the parallel tensor is 
c     B_tot B_tot / B^{2}. 
c     note that beq dot grad Teq is assumed to be zero, a slightly
c     stronger assumption than just a steady state relation.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                       inode,b0,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                       inode,real_bptr,dp,dp,0_i4)
        DO imode=1,nmodes
          bdgr_ti(1,:,:,imode)=be(1,:,:,imode)*ti_eqr(1,:,:)
     $                        +be(2,:,:,imode)*ti_eqz(1,:,:)
        ENDDO

        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,bdgr_ti,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_ti,
     $               real_work,dealiase)
c-----------------------------------------------------------------------
c       Build the product B (B dot grad T).  we are adding 
c       (be_eq+be).grad(T) to be.grad(T_eq).
c
c       if the thermal diffusivity is temperature-dependent,
c       multiply by its distribution in phi here.
c-----------------------------------------------------------------------
        IF (p_model=="aniso_plltdep".OR.p_model=="aniso_tdep") THEN
          IF (.NOT.closure_n0_only) THEN
            DO ip=1,nphi
              real_scal(1,:,ip)=(real_scal(1,:,ip)+
     $                           real_bptr(1,:,ip)*real_work(1,:,ip)+
     $                           real_bptr(2,:,ip)*real_work(2,:,ip)+
     $                           real_bptr(3,:,ip)*real_work(3,:,ip))*
     $            MAX(0._r8,kappli(1,:,ip)-kaprpi(1,:,ip))
            ENDDO
          ELSE
            DO ip=1,nphi
              real_scal(1,:,ip)=(real_scal(1,:,ip)+
     $                           real_bptr(1,:,ip)*real_work(1,:,ip)+
     $                           real_bptr(2,:,ip)*real_work(2,:,ip)+
     $                           real_bptr(3,:,ip)*real_work(3,:,ip))*
     $            MAX(0._r8,kappli(1,:,ip)-elt_1d(ipseust:ipseuen))
            ENDDO
          ENDIF
        ELSE
          real_scal(1,:,:)=
     $      kdiff*(real_scal(1,:,:)+SUM(real_bptr*real_work,1))
        ENDIF
        bbgradt(1,:,:)=real_bptr(1,:,:)*real_scal(1,:,:)
        bbgradt(2,:,:)=real_bptr(2,:,:)*real_scal(1,:,:)
        bbgradt(3,:,:)=real_bptr(3,:,:)*real_scal(1,:,:)
c-----------------------------------------------------------------------
c	Convert the above vector to the Fourier representation
c	and store in be_tot.  the denominator, B**2 is approximated
c       as (be_n0+be_eq)**2.
c
c       include fully 3D temperature_dependent perpendicular
c       thermal conduction if specified.
c-----------------------------------------------------------------------
        IF (.NOT.closure_n0_only) THEN
          magBsq3D(:,:)=SUM(real_bptr(:,:,:)**2,1)
          bbgradt(1,:,:)=bbgradt(1,:,:)/magBsq3D(:,:)
     $         +kaprpi(1,:,:)*real_work(1,:,:)
          bbgradt(2,:,:)=bbgradt(2,:,:)/magBsq3D(:,:)
     $         +kaprpi(1,:,:)*real_work(2,:,:)
          bbgradt(3,:,:)=bbgradt(3,:,:)/magBsq3D(:,:)
     $         +kaprpi(1,:,:)*real_work(3,:,:)
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
c-----------------------------------------------------------------------
c     For linear calculations, the parallel tensor is divided 
c     into linear contributions with B^{2} approximated as Beq^{2}.
c-----------------------------------------------------------------------
      ELSE
        bmagsq=SUM(be_eq*be_eq,1)
        DO imode=1,nmodes
          bdgr_ti(1,:,:,imode)=kdiff*
     $      (be(1,:,:,imode)*ti_eqr(1,:,:)
     $      +be(2,:,:,imode)*ti_eqz(1,:,:)
     $      +SUM(grad_ti(:,:,:,imode)*be_eq,1))/bmagsq
          be_tot(1,:,:,imode)=be_eq(1,:,:)*bdgr_ti(1,:,:,imode)
          be_tot(2,:,:,imode)=be_eq(2,:,:)*bdgr_ti(1,:,:,imode)
          be_tot(3,:,:,imode)=be_eq(3,:,:)*bdgr_ti(1,:,:,imode)
c-----------------------------------------------------------------------
c         add the Be_eqXgrad(Ti) and Be_eqXgrad(Ti_eq) terms here.
c-----------------------------------------------------------------------
          IF (separate_pe.AND.k_cross>0) THEN
            CALL math_cart_cross(b_cr_gt,be_eq,grad_ti(:,:,:,imode),
     $                           -1._r8)
            b_cr_gt(1,:,:)=b_cr_gt(1,:,:)+be(3,:,:,imode)*ti_eqz(1,:,:)
            b_cr_gt(2,:,:)=b_cr_gt(2,:,:)-be(3,:,:,imode)*ti_eqr(1,:,:)
            b_cr_gt(3,:,:)=b_cr_gt(3,:,:)-be(1,:,:,imode)*ti_eqz(1,:,:)+
     $                                    be(2,:,:,imode)*ti_eqr(1,:,:)
            cr_coef=ti_b2(1,:,:)*(2*SUM(be_eq*be(:,:,:,imode),1)/bmagsq-
     $                            nd(1,:,:,imode)/nd_eq(1,:,:)-
     $                            ti(1,:,:,imode)/ti_eq(1,:,:))
            be_tot(1,:,:,imode)=be_tot(1,:,:,imode)+
     $                          b_cr_gt(1,:,:)*ti_b2(1,:,:)+
     $                          bcrgti (1,:,:)*cr_coef
            be_tot(2,:,:,imode)=be_tot(2,:,:,imode)+
     $                          b_cr_gt(2,:,:)*ti_b2(1,:,:)+
     $                          bcrgti (2,:,:)*cr_coef
            be_tot(3,:,:,imode)=be_tot(3,:,:,imode)+
     $                          b_cr_gt(3,:,:)*ti_b2(1,:,:)+
     $                          bcrgti (3,:,:)*cr_coef
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if the perpendicular thermal diffusivity is temperature-dependent,
c     multiply grad_ti by its distribution.  also add parallel
c     contributions.
c
c     nonlinear temperature-dependent perpendicular thermal conduction
c     is coded with the parallel contribution above.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.p_model=="aniso_tdep") THEN
        IF (.NOT.closure_n0_only) THEN
          grad_ti=dt*be_tot
        ELSE
          DO imode=1,nmodes
            grad_ti(1,:,:,imode)=
     $        dt*(kaprpi(1,:,:)*grad_ti(1,:,:,imode)+
     $        be_tot(1,:,:,imode))
            grad_ti(2,:,:,imode)=
     $        dt*(kaprpi(1,:,:)*grad_ti(2,:,:,imode)+
     $            be_tot(2,:,:,imode))
            grad_ti(3,:,:,imode)=
     $        dt*(kaprpi(1,:,:)*grad_ti(3,:,:,imode)+
     $            be_tot(3,:,:,imode))
          ENDDO
        ENDIF
      ELSE
        grad_ti=dt*(k_perpi*grad_ti+be_tot)
      ENDIF

      ENDIF pmodelif
c-----------------------------------------------------------------------
c     also add the upwinding-like diffusion.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $               real_ve,dealiase)
        grad_ti=grad_ti+be_tot
      ENDIF
c-----------------------------------------------------------------------
c     multiply by number density.
c-----------------------------------------------------------------------
      IF (.NOT.nonlinear.OR.continuity=='fix profile') THEN
        DO imode=1,nmodes
          DO jq=1,3
            grad_ti(jq,:,:,imode)=
     $        nd_eq(1,:,:)/zeff*grad_ti(jq,:,:,imode)
          ENDDO
          vdgr_ti(1,:,:,imode)=nd_eq(1,:,:)/zeff*vdgr_ti(1,:,:,imode)
        ENDDO
      ELSE IF (continuity=='n=0 only') THEN
        DO imode=1,nmodes
          DO jq=1,3
            grad_ti(jq,:,:,imode)=
     $        (nd_eq(1,:,:)+nd_n0(1,:,:))/zeff*grad_ti(jq,:,:,imode)
          ENDDO
          vdgr_ti(1,:,:,imode)=(nd_eq(1,:,:)+nd_n0(1,:,:))/zeff*
     $                         vdgr_ti(1,:,:,imode)
        ENDDO
      ELSE ! continuity = full
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_ti,
     $               real_work,dealiase)
        real_scal=real_scal*real_ndptr/zeff
        DO jq=1,3
          real_work(jq,:,:)=real_work(jq,:,:)*real_ndptr(1,:,:)/zeff
        ENDDO
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_ti,
     $               real_scal,dealiase)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,grad_ti,
     $               real_work,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     add the divv_eq term that was multiplied by n-tilde earlier, and
c     form the contributions to each basis function.  ohm_en here has
c     the ndiff correction.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (visc_heat) THEN
          IF (separate_pe) THEN
            vdgr_ti(:,:,:,imode)=vdgr_ti(:,:,:,imode)+
     $        gamm1*vheat(:,:,:,imode)/kboltz
          ELSE
            vdgr_ti(:,:,:,imode)=vdgr_ti(:,:,:,imode)+
     $        (1._r8-pe_frac)*gamm1*vheat(:,:,:,imode)/kboltz
          ENDIF
        ENDIF
        vdgr_ti(1,:,:,imode)=vdgr_ti(1,:,:,imode)+
     $     (0,1)*keff(imode)*grad_ti(3,:,:,imode)/bigr+
     $     dt*gamm1*ohm_en(1,:,:,imode)/kboltz+tdivv(1,:,:,imode)
        DO iv=1,nv
          int(1,:,iv,imode)=SUM(  alpha(:,:,iv)*vdgr_ti(1,:,:,imode)
     $                          -dalpdr(:,:,iv)*grad_ti(1,:,:,imode)
     $                          -dalpdz(:,:,iv)*grad_ti(2,:,:,imode),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tirhs
c-----------------------------------------------------------------------
c     subprogram 8. terhs.
c     compute the integrand used for the rhs of the electron temperature
c     change equation with implicit diffusion.
c
c     dt*(gamma-1)*(D(Te_n)-ni*Te*div(Ve))-dt*ni*Ve.grad(Te)
c
c     where Te_n is the temperature at the beginning of the time step,
c     Te is either Te_n or the predicted temperature, and D is the
c     diffusion operator.
c-----------------------------------------------------------------------
      SUBROUTINE terhs(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          te_eq,te_eqr,te_eqz,
     $          ve_eq,be_eq,ja_eq,b0,q_a,nd_eq,nd_eqr,
     $          nd_eqz,nd_n0,real_bptr,real_ndptr,elecd_t,elecd_eq,dp,
     $          ds2,kapple,kaprpe,te_b2,bcrgte,divv_eq,real_vptr,
     $          real_jptr,upwc
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_work,bbgradt
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal,real_divv
      REAL(r8), DIMENSION(mpseudo,nphi) :: magBsq3D
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: eleq_1d,elt_1d
      REAL(r8), DIMENSION(  SIZE(bigr,1)*SIZE(bigr,2)) :: tqr_1d,tqz_1d
      REAL(r8), DIMENSION(3,SIZE(bigr,1)*SIZE(bigr,2)) :: jaeq_1d 
      REAL(r8), DIMENSION(3,mpseudo) :: rvec
      REAL(r8), DIMENSION(  mpseudo) :: jaeq2,rdot
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: kdiff
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             ve,dvr,dvz,be,ber,bez,
     $             te,ter,tez,nd,ndr,ndz,hypheat,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             grad_te,be_tot,grad_nd,vel,ja
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             vdgr_te,divv,tdivv,bdgr_te,ohm_en
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: b_cr_gt
      COMPLEX(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: cr_coef
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,imode,iq,jq,ncx,ncy,ip
      REAL(r8) :: fupw
      LOGICAL :: told_set
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      told_set=.false.
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     temperature and grad(T) at the beginning of this time split for
c     predictor steps, or at the centered prediction for corrector
c     steps.  we assume ja_eq/=0, so cases that call this integrand
c     always do predictor/corrector steps.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                     inode,te_eq,te_eqr,te_eqz,1_i4)
      IF (integrand_flag=='all te terms pre'.OR.
     $    k_cross>0.OR.impladv) THEN
        CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,inode,
     $                       te,ter,tez,1_i4)
        told_set=.true.
      ELSE
        CALL generic_ptr_set(rb%qwork3,tb%qwork3,tb%tgeom,inode,
     $                       te,ter,tez,1_i4)
      ENDIF
      DO imode=1,nmodes
        grad_te(1,:,:,imode)=ter(1,:,:,imode)
        grad_te(2,:,:,imode)=tez(1,:,:,imode)
        grad_te(3,:,:,imode)=(0,1)*keff(imode)*te(1,:,:,imode)/bigr
      ENDDO
      IF (nonlinear.AND.ohm_heat.AND.(hyp_eta>0.OR.hyp_dbd>0)) THEN
        CALL generic_ptr_set(rb%qhyph,tb%qhyph,tb%tgeom,inode,hypheat,
     $                       dcp,dcp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate electron velocity.
c
c     The electron velocity is computed from the COM velocity 
c     and the total current density as
c
c          vel = v - j / ( n e (nu+1) )
c
c     where "n" is the electron number density, "e" the electron
c     charge and "nu" is the electron-ion mass ratio times the atomic
c     number "Z":
c                  nu = Z me / mi
c
c     we also need div(v)=d(Vr)/dR + Vr/R + d(Vz)dZ + 1/R d(Vth)/dTh
c     here we assume that Jeq.grad(nd_eq)=0.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,dvr,dvz,1_i4)
      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,nd,ndr,ndz,1_i4)
      IF (nonlinear.AND.continuity=='n=0 only')
     $  CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,nd_eqr,nd_eqz,1_i4)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
      CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8/mu0)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,ja_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     extra pointers are needed for the upwinding-like artificial
c     diffusion.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,real_vptr,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qupte_phi,tb%qupte_phi,tb%tgeom,
     $                       inode,upwc,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qja_tot,tb%qja_tot,tb%tgeom,
     $                       inode,real_jptr,dp,dp,0_i4)
        fupw=t_dart_upw*dt**2
        tqr_1d=RESHAPE(te_eqr,(/ncx*ncy/))
        tqz_1d=RESHAPE(te_eqz,(/ncx*ncy/))
      ENDIF
c-----------------------------------------------------------------------
c     div(v-COM) and grad(n)
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (geom=='tor') THEN
          divv(1,:,:,imode)=dvr(1,:,:,imode)+dvz(2,:,:,imode)
     $                    +((0,1)*keff(imode)*ve(3,:,:,imode)
     $                                       +ve(1,:,:,imode))/bigr
        ELSE
          divv(1,:,:,imode)=dvr(1,:,:,imode)+dvz(2,:,:,imode)
     $                     +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr
        ENDIF
        grad_nd(1,:,:,imode)=ndr(1,:,:,imode)
        grad_nd(2,:,:,imode)=ndz(1,:,:,imode)
        grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
      ENDDO

      IF (nonlinear)
     $  CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)

      IF (eq_flow/='none')
     $  CALL generic_ptr_set(rb%qdvv_eq,tb%qdvv_eq,tb%tgeom,inode,
     $                       divv_eq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c     find the Ohmic heating energy before the ja array is modified.
c-----------------------------------------------------------------------
      IF (nonlinear)
     $  CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,real_ve,
     $               dealiase)
      IF (ohm_heat) THEN
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_eq,tb%qelecd_eq,tb%tgeom,
     $                         inode,elecd_eq,dp,dp,0_i4)
          eleq_1d=RESHAPE(elecd_eq,(/ncx*ncy/))
          jaeq_1d=RESHAPE( ja_eq,(/3,ncx*ncy/))
          DO iq=1,mpseudo
            jaeq2(iq)=SUM(jaeq_1d(:,iq-1+ipseust)**2)
          ENDDO
        ENDIF
        DO imode=1,nmodes
          ohm_en(1,:,:,imode)=2*SUM(ja_eq*ja(:,:,:,imode),1)
        ENDDO
        IF (nonlinear) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,
     $                 ohm_en,real_scal,dealiase)
          real_scal(1,:,:)=real_scal(1,:,:)+SUM(real_ve**2,1)
          IF (threedeta) THEN
            CALL generic_ptr_set(rb%qelecd_phi,tb%qelecd_phi,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            DO ip=1,nphi
              DO iq=1,mpseudo
                real_scal(1,iq,ip)=mu0*
     $            (real_scal(1,iq,ip)*elecd_t(1,iq,ip)+
     $            (elecd_t(1,iq,ip)-eleq_1d(iq-1+ipseust))*jaeq2(iq))
              ENDDO
            ENDDO
          ELSE IF (eta_model=='eta n=0 only') THEN
            CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            elt_1d=RESHAPE(elecd_t,(/ncx*ncy/))
            DO iq=1,mpseudo
              jaeq2(iq)=jaeq2(iq)*
     $              (elt_1d(iq-1+ipseust)-eleq_1d(iq-1+ipseust))
            ENDDO
            DO ip=1,nphi
              DO iq=1,mpseudo
                real_scal(1,iq,ip)=mu0*
     $            (real_scal(1,iq,ip)*elt_1d(iq-1+ipseust)+jaeq2(iq))
              ENDDO
            ENDDO
          ENDIF
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,
     $                 ohm_en,real_scal,dealiase)
        ENDIF
        IF (nonlinear.AND.threedeta) THEN  ! done
        ELSE IF (eta_model=='eta n=0 only') THEN  ! done if nonlinear
          IF (.NOT.nonlinear) THEN
            CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                           inode,elecd_t,dp,dp,0_i4)
            DO imode=1,nmodes
              ohm_en(:,:,:,imode)=ohm_en(:,:,:,imode)*elecd_t*mu0
            ENDDO
          ENDIF
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          DO imode=1,nmodes
            ohm_en(:,:,:,imode)=
     $        ohm_en(:,:,:,imode)*elecd*MAX(ds2,0._r8)*mu0
          ENDDO
        ELSE
          ohm_en=ohm_en*elecd*mu0
        ENDIF
        IF (nonlinear.AND.(hyp_eta>0.OR.hyp_dbd>0))
     $    ohm_en=ohm_en+hypheat
      ENDIF
c-----------------------------------------------------------------------
c     j.grad(n) (saved in vdgr_te) for nonlinear, full continuity
c     computations.  [keep div(j)=0 in mind here]
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $               real_work,dealiase)
        real_scal(1,:,:)=SUM(real_ve*real_work,1)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $               real_scal,dealiase)
      ELSE
        vdgr_te=0._r8
      ENDIF

      IF (coefjve==0._r8) THEN
        vel=ve
      ELSE IF (.NOT.nonlinear.OR.continuity=='fix profile') THEN
        DO imode=1,nmodes
          DO jq=1,3
            vel(jq,:,:,imode)=ve(jq,:,:,imode)
     $        +ja(jq,:,:,imode)*coefjve/nd_eq(1,:,:)
          ENDDO
          divv(1,:,:,imode)=divv(1,:,:,imode)-
     $      coefjve*(SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $               +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $               +ja(2,:,:,imode)*nd_eqz(1,:,:)
     $               +vdgr_te(1,:,:,imode))/(nd_eq(1,:,:)**2)
        ENDDO
      ELSE IF (continuity=='n=0 only') THEN
        bmagsq=1._r8/(nd_eq(1,:,:)+nd_n0(1,:,:))
        DO imode=1,nmodes
          DO jq=1,3
            vel(jq,:,:,imode)=ve(jq,:,:,imode)
     $        +ja(jq,:,:,imode)*coefjve*bmagsq
          ENDDO
          divv(1,:,:,imode)=divv(1,:,:,imode)-
     $      coefjve*(SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $               +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $               +ja(2,:,:,imode)*nd_eqz(1,:,:)
     $               +vdgr_te(1,:,:,imode))*bmagsq**2
        ENDDO
      ELSE  !  continuity is full  (ja will be left as ja/n).
        DO imode=1,nmodes
          vdgr_te(1,:,:,imode)=-coefjve*(vdgr_te(1,:,:,imode)+
     $      SUM(ja_eq*grad_nd(:,:,:,imode),1)
     $      +ja(1,:,:,imode)*nd_eqr(1,:,:)
     $      +ja(2,:,:,imode)*nd_eqz(1,:,:))
        ENDDO
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $               real_scal,dealiase)
        DO jq=1,3
          real_ve(jq,:,:)=real_ve(jq,:,:)/real_ndptr(1,:,:)
        ENDDO
        real_scal=real_scal/real_ndptr**2
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,ja,
     $               real_ve,dealiase)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $               real_scal,dealiase)
        vel=ve+coefjve*ja
        divv=divv+vdgr_te
      ENDIF
c-----------------------------------------------------------------------
c     find the linear v_tilde.grad(Teq) for each mode--needed here
c     to complete the upwinding-like diffusion at the next step.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        vdgr_te(1,:,:,imode)=vel(1,:,:,imode)*te_eqr(1,:,:)+
     $                       vel(2,:,:,imode)*te_eqz(1,:,:)
      ENDDO
c-----------------------------------------------------------------------
c     find the nonlinear contributions to v.grad(t) and t div(v)
c
c     transform v_perturbed and grad(t_perturbed) to real space,
c     find v.grad(t), then transform to Fourier space.  tdivv is
c     just used as temporary storage.
c
c     the upwinding-like diffusion term is also computed here.  the
c     predetermined coefficient is multiplied by VV.grad(T).
c     ** NOTE ** real_ve is now used to save this coefficient from
c     this point, so V_tilde(phi) array is lost.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,vel,real_ve,
     $               dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_te,
     $               real_work,dealiase)
        real_scal(1,:,:)=SUM(real_ve*real_work,1)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tdivv,
     $               real_scal,dealiase)

        IF (impladv.AND.t_dart_upw>0) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $                 real_scal,dealiase)
          DO ip=1,nphi
            rvec(1,:)=real_vptr(1,:,ip)+
     $           coefjve*real_jptr(1,:,ip)/real_ndptr(1,:,ip)
            rvec(2,:)=real_vptr(2,:,ip)+
     $           coefjve*real_jptr(2,:,ip)/real_ndptr(1,:,ip)
            rvec(3,:)=real_vptr(3,:,ip)+
     $           coefjve*real_jptr(3,:,ip)/real_ndptr(1,:,ip)
            rdot=fupw*upwc(1,:,ip)*upw_aniso*
     $        ( real_scal(1,:,ip) + SUM(rvec*real_work(:,:,ip),1) ) /
     $        ( SUM(rvec*rvec,1) + smallnum )
            jaeq2=fupw*upwc(1,:,ip)*(1._r8-upw_aniso)
            real_ve(1,:,ip)=rvec(1,:)*rdot+
     $        (real_work(1,:,ip)+tqr_1d(ipseust:ipseuen))*jaeq2
            real_ve(2,:,ip)=rvec(2,:)*rdot+
     $        (real_work(2,:,ip)+tqz_1d(ipseust:ipseuen))*jaeq2
            real_ve(3,:,ip)=rvec(3,:)*rdot+real_work(3,:,ip)*jaeq2
          ENDDO
        ENDIF
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     subtract the rest of the linear V.grad(T) for each mode.
c-----------------------------------------------------------------------
      IF (eq_flow/='none') THEN
        DO imode=1,nmodes
          vdgr_te(1,:,:,imode)=-dt*( vdgr_te(1,:,:,imode)+
     $      tdivv(1,:,:,imode)+SUM(ve_eq*grad_te(:,:,:,imode),1)+
     $      coefjve*SUM(ja_eq*grad_te(:,:,:,imode),1)/nd_eq(1,:,:) )
        ENDDO
      ELSE
        DO imode=1,nmodes
          vdgr_te(1,:,:,imode)=-dt*( vdgr_te(1,:,:,imode)+
     $                                 tdivv(1,:,:,imode)+
     $      coefjve*SUM(ja_eq*grad_te(:,:,:,imode),1)/nd_eq(1,:,:) )
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     multiply t_perturbed by div(v), then transform the
c     product to Fourier space.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.p_model/='isothermal') THEN
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,te,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,divv,
     $               real_divv,dealiase)
        real_scal=real_scal*real_divv
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,tdivv,
     $               real_scal,dealiase)
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     subtract the linear (gamma-1)*T*div(V) for this mode.
c     also form the n_tilde*Teq*div(Veq) term and save in tdivv.
c-----------------------------------------------------------------------
      IF (p_model/='isothermal') THEN
        IF (eq_flow/='none') THEN
          DO imode=1,nmodes
            vdgr_te(1,:,:,imode)=vdgr_te(1,:,:,imode)
     $                         -dt*gamm1*(tdivv(1,:,:,imode)
     $                              +te_eq(1,:,:)*divv(1,:,:,imode)
     $                              +te(1,:,:,imode)*divv_eq(1,:,:))
            tdivv(1,:,:,imode)=-dt*gamm1*nd(1,:,:,imode)*
     $                               te_eq(1,:,:)*divv_eq(1,:,:)
          ENDDO
        ELSE
          DO imode=1,nmodes
            vdgr_te(1,:,:,imode)=vdgr_te(1,:,:,imode)
     $                         -dt*gamm1*(tdivv(1,:,:,imode)
     $                              +te_eq(1,:,:)*divv(1,:,:,imode))
            tdivv(1,:,:,imode)=0._r8
          ENDDO
        ENDIF
      ELSE
        tdivv=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     set temperature and grad(T) with T from the beginning of this
c     time split for thermal conduction if not already done.  
c-----------------------------------------------------------------------
      IF (p_model=='adiabat'.OR.p_model=='isothermal') THEN
        grad_te=0._r8
      ELSE
        IF (.NOT.told_set) THEN
          CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,inode,
     $                         te,ter,tez,1_i4)
          DO imode=1,nmodes
            grad_te(1,:,:,imode)=ter(1,:,:,imode)
            grad_te(2,:,:,imode)=tez(1,:,:,imode)
            grad_te(3,:,:,imode)=(0,1)*keff(imode)*te(1,:,:,imode)/bigr
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c-PRE Heat soure term.
c     q_applied is set in q_applied_init.f and is normally zero
c-----------------------------------------------------------------------
c     DO imode=1,nmodes
c       IF (keff(imode)==0) THEN
c         CALL generic_ptr_set(rb%qq_applied,tb%qq_applied
c    $                         tb%tgeom,inode,q_a,dp,dp,0_i4)
c         vdgr_te(1,:,:,imode)=vdgr_te(1,:,:,imode)
c    $                         +dt*gamm1*q_a/(kboltz*bigr)
c       ENDIF
c     ENDDO
c-----------------------------------------------------------------------
c     For isotropic diffusion, add the grad k_perp grad T old with
c     the adiabatic terms.  assume k_perp has gamma-1 factor in it.
c-----------------------------------------------------------------------
      pmodelif: IF (p_model=='isotropic') THEN
        grad_te=grad_te*dt*k_perpe
      ELSE IF (p_model(1:5)=='aniso') THEN
c-----------------------------------------------------------------------
c     for anisotropic thermal conduction, get the
c     magnetic field first.  note that we assume k_pll and k_perp have
c     the factor of (gamma-1) in them already.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      IF (nonlinear.AND.(p_model=="aniso_plltdep".OR.
     $                   p_model=="aniso_tdep")) THEN
        CALL generic_ptr_set(rb%qkapple_phi,tb%qkapple_phi,tb%tgeom,
     $                       inode,kapple,dp,dp,0_i4)
        IF (.NOT.closure_n0_only) THEN
          CALL generic_ptr_set(rb%qkaprpe_phi,tb%qkaprpe_phi,tb%tgeom,
     $                         inode,kaprpe,dp,dp,0_i4)
        ELSE
          CALL generic_ptr_set(rb%qkaprpe_n0,tb%qkaprpe_n0,tb%tgeom,
     $                         inode,kaprpe,dp,dp,0_i4)
          elt_1d=RESHAPE(kaprpe,(/ncx*ncy/))
        ENDIF
      ELSE
        kdiff=k_plle-k_perpe
      ENDIF
      IF (k_cross>0) THEN
        CALL generic_ptr_set(rb%qte_b2,tb%qte_b2,tb%tgeom,
     $                       inode,te_b2,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qbcrgte,tb%qbcrgte,tb%tgeom,
     $                       inode,bcrgte,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     For nonlinear calculations, the parallel tensor is 
c     B_tot B_tot / B^{2}. 
c     note that beq dot grad Teq is assumed to be zero, a slightly
c     stronger assumption than just a steady state relation.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                       inode,b0,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                       inode,real_bptr,dp,dp,0_i4)
        DO imode=1,nmodes
          bdgr_te(1,:,:,imode)=be(1,:,:,imode)*te_eqr(1,:,:)
     $                        +be(2,:,:,imode)*te_eqz(1,:,:)
        ENDDO

        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,bdgr_te,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_te,
     $               real_work,dealiase)
c-----------------------------------------------------------------------
c       Build the product B (B dot grad T).  we are adding 
c       (be_eq+be).grad(T) to be.grad(T_eq).
c
c       if the thermal conductivity is temperature-dependent
c       multiply by its distribution in phi here.
c-----------------------------------------------------------------------
        IF (p_model=="aniso_plltdep".OR.p_model=="aniso_tdep") THEN
          IF (.NOT.closure_n0_only) THEN
            DO ip=1,nphi
              real_scal(1,:,ip)=(real_scal(1,:,ip)+
     $                           real_bptr(1,:,ip)*real_work(1,:,ip)+
     $                           real_bptr(2,:,ip)*real_work(2,:,ip)+
     $                           real_bptr(3,:,ip)*real_work(3,:,ip))*
     $            MAX(0._r8,kapple(1,:,ip)-kaprpe(1,:,ip))
            ENDDO
          ELSE
            DO ip=1,nphi
              real_scal(1,:,ip)=(real_scal(1,:,ip)+
     $                           real_bptr(1,:,ip)*real_work(1,:,ip)+
     $                           real_bptr(2,:,ip)*real_work(2,:,ip)+
     $                           real_bptr(3,:,ip)*real_work(3,:,ip))*
     $            MAX(0._r8,kapple(1,:,ip)-elt_1d(ipseust:ipseuen))
            ENDDO
          ENDIF
        ELSE
          real_scal(1,:,:)=
     $      kdiff*(real_scal(1,:,:)+SUM(real_bptr*real_work,1))
        ENDIF
        bbgradt(1,:,:)=real_bptr(1,:,:)*real_scal(1,:,:)
        bbgradt(2,:,:)=real_bptr(2,:,:)*real_scal(1,:,:)
        bbgradt(3,:,:)=real_bptr(3,:,:)*real_scal(1,:,:)
c-----------------------------------------------------------------------
c	Convert the above vector to the Fourier representation
c	and store in be_tot.  the denominator, B**2 is approximated
c       as (be_n0+be_eq)**2.
c
c       include fully 3D temperature_dependent perpendicular thermal
c       conduction if specified.
c-----------------------------------------------------------------------
        IF (.NOT.closure_n0_only) THEN
          magBsq3D(:,:)=SUM(real_bptr(:,:,:)**2,1)
          bbgradt(1,:,:)=bbgradt(1,:,:)/magBsq3D(:,:)
     $        +kaprpe(1,:,:)*real_work(1,:,:)
          bbgradt(2,:,:)=bbgradt(2,:,:)/magBsq3D(:,:)
     $        +kaprpe(1,:,:)*real_work(2,:,:)
          bbgradt(3,:,:)=bbgradt(3,:,:)/magBsq3D(:,:)
     $        +kaprpe(1,:,:)*real_work(3,:,:)
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $        bbgradt,dealiase)
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
c-----------------------------------------------------------------------
c     For linear calculations, the parallel tensor is divided 
c     into linear contributions with B^{2} approximated as Beq^{2}.
c-----------------------------------------------------------------------
      ELSE
        bmagsq=SUM(be_eq*be_eq,1)
        DO imode=1,nmodes
          bdgr_te(1,:,:,imode)=kdiff*
     $      (be(1,:,:,imode)*te_eqr(1,:,:)
     $      +be(2,:,:,imode)*te_eqz(1,:,:)
     $      +SUM(grad_te(:,:,:,imode)*be_eq,1))/bmagsq
          be_tot(1,:,:,imode)=be_eq(1,:,:)*bdgr_te(1,:,:,imode)
          be_tot(2,:,:,imode)=be_eq(2,:,:)*bdgr_te(1,:,:,imode)
          be_tot(3,:,:,imode)=be_eq(3,:,:)*bdgr_te(1,:,:,imode)
c-----------------------------------------------------------------------
c         add the Be_eqXgrad(Te) and Be_eqXgrad(Te_eq) terms here.
c-----------------------------------------------------------------------
          IF (k_cross>0) THEN
            CALL math_cart_cross(b_cr_gt,be_eq,grad_te(:,:,:,imode),
     $                           1._r8)
            b_cr_gt(1,:,:)=b_cr_gt(1,:,:)-be(3,:,:,imode)*te_eqz(1,:,:)
            b_cr_gt(2,:,:)=b_cr_gt(2,:,:)+be(3,:,:,imode)*te_eqr(1,:,:)
            b_cr_gt(3,:,:)=b_cr_gt(3,:,:)+be(1,:,:,imode)*te_eqz(1,:,:)-
     $                                    be(2,:,:,imode)*te_eqr(1,:,:)
            cr_coef=te_b2(1,:,:)*(nd(1,:,:,imode)/nd_eq(1,:,:)+
     $                            te(1,:,:,imode)/te_eq(1,:,:)-
     $                            2*SUM(be_eq*be(:,:,:,imode),1)/bmagsq)
            be_tot(1,:,:,imode)=be_tot(1,:,:,imode)+
     $                          b_cr_gt(1,:,:)*te_b2(1,:,:)+
     $                          bcrgte (1,:,:)*cr_coef
            be_tot(2,:,:,imode)=be_tot(2,:,:,imode)+
     $                          b_cr_gt(2,:,:)*te_b2(1,:,:)+
     $                          bcrgte (2,:,:)*cr_coef
            be_tot(3,:,:,imode)=be_tot(3,:,:,imode)+
     $                          b_cr_gt(3,:,:)*te_b2(1,:,:)+
     $                          bcrgte (3,:,:)*cr_coef
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if the perpendicular thermal diffusivity is temperature-dependent,
c     multiply grad_te by its distribution.  also add parallel
c     contributions.
c
c     nonlinear temperature-dependent perpendicular thermal conduction
c     is now above.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.p_model=="aniso_tdep") THEN
        IF (.NOT.closure_n0_only) THEN
          grad_te=dt*be_tot
        ELSE
          DO imode=1,nmodes
            grad_te(1,:,:,imode)=
     $         dt*(kaprpe(1,:,:)*grad_te(1,:,:,imode)+
     $             be_tot(1,:,:,imode))
            grad_te(2,:,:,imode)=
     $         dt*(kaprpe(1,:,:)*grad_te(2,:,:,imode)+
     $             be_tot(2,:,:,imode))
            grad_te(3,:,:,imode)=
     $         dt*(kaprpe(1,:,:)*grad_te(3,:,:,imode)+
     $             be_tot(3,:,:,imode))
          ENDDO
        ENDIF
      ELSE
        grad_te=dt*(k_perpe*grad_te+be_tot)
      ENDIF

      ENDIF pmodelif
c-----------------------------------------------------------------------
c     also add the upwinding-like diffusion.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,be_tot,
     $               real_ve,dealiase)
        grad_te=grad_te+be_tot
      ENDIF
c-----------------------------------------------------------------------
c     multiply by number density.
c-----------------------------------------------------------------------
      IF (.NOT.nonlinear.OR.continuity=='fix profile') THEN
        DO imode=1,nmodes
          DO jq=1,3
            grad_te(jq,:,:,imode)=nd_eq(1,:,:)*grad_te(jq,:,:,imode)
          ENDDO
          vdgr_te(1,:,:,imode)=nd_eq(1,:,:)*vdgr_te(1,:,:,imode)
        ENDDO
      ELSE IF (continuity=='n=0 only') THEN
        DO imode=1,nmodes
          DO jq=1,3
            grad_te(jq,:,:,imode)=
     $        (nd_eq(1,:,:)+nd_n0(1,:,:))*grad_te(jq,:,:,imode)
          ENDDO
          vdgr_te(1,:,:,imode)=(nd_eq(1,:,:)+nd_n0(1,:,:))*
     $                         vdgr_te(1,:,:,imode)
        ENDDO
      ELSE ! continuity = full
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $               real_scal,dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_te,
     $               real_work,dealiase)
        real_scal=real_scal*real_ndptr
        DO jq=1,3
          real_work(jq,:,:)=real_work(jq,:,:)*real_ndptr(1,:,:)
        ENDDO
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,vdgr_te,
     $               real_scal,dealiase)
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,grad_te,
     $               real_work,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     add the divv_eq term that was multiplied by n-tilde earlier, and
c     form the contributions to each basis function.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (ohm_heat) THEN
          vdgr_te(1,:,:,imode)=vdgr_te(1,:,:,imode)+
     $       (0,1)*keff(imode)*grad_te(3,:,:,imode)/bigr+
     $       dt*gamm1*ohm_en(1,:,:,imode)/kboltz+
     $       tdivv(1,:,:,imode)
        ELSE
          vdgr_te(1,:,:,imode)=vdgr_te(1,:,:,imode)+
     $       (0,1)*keff(imode)*grad_te(3,:,:,imode)/bigr+
     $       tdivv(1,:,:,imode)
        ENDIF
        DO iv=1,nv
          int(1,:,iv,imode)=SUM(  alpha(:,:,iv)*vdgr_te(1,:,:,imode)
     $                          -dalpdr(:,:,iv)*grad_te(1,:,:,imode)
     $                          -dalpdz(:,:,iv)*grad_te(2,:,:,imode),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE terhs
c-----------------------------------------------------------------------
c     subprogram 9. vrhs.
c     form the integrand for the rhs of the velocity equation
c     at a single gauss quadrature point for all cells in a block.
c-----------------------------------------------------------------------
      SUBROUTINE vrhs(int,bigr,rb,dx,dy,tb,inode)

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
     $          be_eq,ja_eq,nd_eq,
     $          ve_eq,veqten,nd_n0,real_ve,ti_eq,te_eq,
     $          real_gradv,real_ndptr,real_bptr,real_vten,
     $          pi_veq,pareq,gyreq,kappli,dp,beqr,beqz,be0,be0r,be0z,
     $          peq,p0,real_titot,real_ndiff,ds2
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal
      REAL(r8), DIMENSION(3,mpseudo,nphi), TARGET :: real_ja,real_be,
     $          real_force
      REAL(r8), DIMENSION(9,mpseudo,nphi), TARGET :: real_pten
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: diff,dtmp
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: teff_eq
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: eq_vadv,bes,
     $          besr,besz
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8), DIMENSION(3,3) :: vtmpr,ptmpr,btmp,wdotr,bc_wdotr
      REAL(r8) :: divtmpr,dtmt,rscal,rb2,vcoef,gcoef,g,dtpv,dtdv
      REAL(r8) :: third=1._r8/3._r8,twot=2._r8/3._r8

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             be,ber,bez,ve,ve_r,ve_z,
     $             nd,ndiff,vheat,tion,tele,prsptr,dcp
      COMPLEX(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             pres,teff
      COMPLEX(r8), DIMENSION(2,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: aux
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: ja,
     $             force
      COMPLEX(r8), DIMENSION(9,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             vten,piten,parten
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      COMPLEX(r8), DIMENSION(3,3) :: vtmp,ptmp,wdot,b_cr_wdot,btmpc
      COMPLEX(r8) :: divtmp,cscal

      INTEGER(i4) :: iv,imode,iq,jt,ivind,nv,ncx,ncy,ix,iy,ip,im
      INTEGER(i4) :: i1,i2,nvc,nvm,nvhyp
      LOGICAL :: ve_set,nd_set,aux_set
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
      ve_set=.false.
      nd_set=.false.
      aux_set=.false.
      dtmt=dt*mtot
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the perturbed
c     and 0-th order magnetic fields.  the current density is evaluated
c     from the gradients of the perturbed magnetic field
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,inode,
     $                     be_eq,beqr,beqz,1_i4)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
      IF (integrand_flag(1:3)=='all'.OR.integrand_flag(1:3)=='mhd') THEN
        CALL math_curl(nmodes,keff,geom,bigr,be,ber,bez,ja,1._r8/mu0)
        CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                       inode,ja_eq,dp,dp,0_i4)
        IF (beta>0) THEN
          CALL generic_ptr_set(rb%qtion,tb%qtion,tb%tgeom,inode,
     $                         tion,dcp,dcp,0_i4)
          CALL generic_ptr_set(rb%qtele,tb%qtele,tb%tgeom,inode,
     $                         tele,dcp,dcp,0_i4)
          CALL generic_ptr_set(rb%qtion_eq,tb%qtion_eq,tb%tgeom,
     $                         inode,ti_eq,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                         inode,te_eq,dp,dp,0_i4)
          teff=dt*kboltz*(tele+tion/zeff)
          teff_eq=dt*kboltz*(te_eq+ti_eq/zeff)
          IF (p_computation=='at nodes') THEN
            CALL generic_ptr_set(rb%qpres,tb%qpres,tb%tgeom,inode,
     $                           prsptr,dcp,dcp,0_i4)
            pres=dt*prsptr
          ENDIF
        ELSE
          pres=0._r8
        ENDIF
      ELSE
        pres=0._r8
      ENDIF
      IF (eq_flow/='none') THEN
        CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qgrdveq,tb%qgrdveq,tb%tgeom,
     $                       inode,veqten,dp,dp,0_i4)
      ENDIF
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,
     $                     nd,dcp,dcp,0_i4)
      IF (continuity=='n=0 only') THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                       nd_n0,dp,dp,0_i4)
      ENDIF
      IF (nonlinear.AND.continuity/='none') THEN
        CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     if poly_divv is non-negative, auxiliary fields are used to
c     stabilize the V-advance at the limit of element resolution.  the
c     modal bases for this field are then needed.
c-----------------------------------------------------------------------
      nvm=0
      IF (poly_divv>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlrhs',almod,dalmr,
     $                          dalmz,0_i4,poly_divv,
     $                          polydmin=poly_divv_min,
     $                          polydmax=poly_divv_max)
        nvm=SIZE(almod,3)
      ENDIF
c-----------------------------------------------------------------------
c     choose mhd forces, viscous & advective forces, or both according
c     to the time splitting.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:3)=='mhd'.OR.integrand_flag(1:3)=='all') THEN
c-----------------------------------------------------------------------
c       find the nonlinear contributions to the force density.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
c-----------------------------------------------------------------------
c         first find b_perturbed and j_perturbed in real space.
c-----------------------------------------------------------------------
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,be,real_be,
     $                 dealiase)
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,real_ja,
     $                 dealiase)
c-----------------------------------------------------------------------
c         find j_perturbed X b_perturbed in real space,
c         then transform to Fourier space.
c-----------------------------------------------------------------------
          CALL math_cart_cross(real_force,real_ja,real_be,dt)
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,force,
     $                 real_force,dealiase)
c-----------------------------------------------------------------------
c         find pressure from the product of n and Teff=Te+Ti/zeff.
c         the real_ndptr storage has nd_eq already added.  the
c         real_scal array is used for temporary storage.
c-----------------------------------------------------------------------
          IF (beta>0.AND.p_computation/='at nodes') THEN
            CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,teff,
     $                   real_scal,dealiase)
            real_scal=real_scal*real_ndptr
            CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,1_i4,pres,
     $                   real_scal,dealiase)
          ENDIF
        ELSE
          force=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       construct JXB and P separately
c-----------------------------------------------------------------------
        mode_loop: DO imode=1,nmodes
          CALL math_cadd_cross(force(:,:,:,imode),
     $                         ja(:,:,:,imode),be_eq,dt)
          CALL math_cadd_cross(force(:,:,:,imode),
     $                         ja_eq,be(:,:,:,imode),dt)
          IF (beta>0.AND.p_computation/='at nodes') THEN
            IF (nonlinear) THEN
              pres(:,:,:,imode)=pres(:,:,:,imode)+
     $                            nd(:,:,:,imode)*teff_eq
            ELSE
              pres(:,:,:,imode)=teff(:,:,:,imode)*nd_eq+
     $                            nd(:,:,:,imode)*teff_eq
            ENDIF
          ENDIF
        ENDDO mode_loop
      ELSE
        force=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     explicit gravitational acceleration.
c-----------------------------------------------------------------------
      IF (MAXVAL(ABS(gravity))>0._r8) THEN
        force(1,:,:,:)=force(1,:,:,:)+nd(1,:,:,:)*dtmt*gravity(1)
        force(2,:,:,:)=force(2,:,:,:)+nd(1,:,:,:)*dtmt*gravity(2)
        force(3,:,:,:)=force(3,:,:,:)+nd(1,:,:,:)*dtmt*gravity(3)
      ENDIF
c-----------------------------------------------------------------------
c     the vten grad(V) tensor is in Fourier representation and is left
c     without equilibrium contributions.  real_vten is in real space
c     with equilibrium contributions.
c-----------------------------------------------------------------------
      IF ((advect=='V only'.OR.advect=='all') .AND.
     $    (nonlinear.OR.eq_flow/='none') .AND.
     $    (integrand_flag(1:3)=='all'.OR.integrand_flag(9:12)=='vdgv')
     $   .OR.poly_divv>=0) THEN
        IF (integrand_flag(14:16)=='pre'.OR.
     $      integrand_flag(12:14)=='pre'.OR.impladv) THEN
          CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,
     $                         inode,ve,ve_r,ve_z,1_i4)
          CALL math_grad(nmodes,keff,3_i4,geom,ve,ve_r,ve_z,vten,bigr)
        ELSE     !     data for corrector step of p/c advection
          CALL generic_ptr_set(rb%qwork1,tb%qwork1,tb%tgeom,
     $                         inode,ve,ve_r,ve_z,1_i4)
          CALL math_grad(nmodes,keff,3_i4,geom,ve,ve_r,ve_z,vten,bigr)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     subtract alpha_i*rho*V.grad(V) when momentum density advection
c     is requested using the appropriate velocity (either old or
c     predicted, stored in work1). the copies are done here to reuse
c     memory from the viscous term, and pointers are used just to
c     rename the configuration-space variables.
c
c     real_vten is in real space with equilibrium contributions.
c-----------------------------------------------------------------------
      IF ((advect=='V only'.OR.advect=='all') .AND.
     $    (nonlinear.OR.eq_flow/='none') .AND.
     $    (integrand_flag(1:3)=='all'.OR.integrand_flag(9:12)=='vdgv')
     $   ) THEN
        IF (integrand_flag(14:16)=='pre'.OR.
     $      integrand_flag(12:14)=='pre'.OR.impladv) THEN
          IF (nonlinear) THEN
            CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                           inode,real_vten,dp,dp,0_i4)
          ENDIF
          ve_set=.true.
        ELSE     !     data for corrector step of p/c advection
          IF (nonlinear) THEN
            IF (eq_flow/='none') THEN
              DO imode=1,nmodes
                IF (keff(imode)==0)
     $            vten(:,:,:,imode)=vten(:,:,:,imode)+veqten
              ENDDO
            ENDIF
            real_vten=>real_pten
            CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,9_i4,vten,
     $                   real_vten,dealiase)
            IF (eq_flow/='none') THEN
              DO imode=1,nmodes
                IF (keff(imode)==0)
     $            vten(:,:,:,imode)=vten(:,:,:,imode)-veqten
              ENDDO
            ENDIF
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       evolving density requires a V0.grad(V0) term.
c-----------------------------------------------------------------------
        IF (eq_flow/='none'.AND.continuity/='none') THEN
          eq_vadv(1,:,:)=dtmt*SUM(ve_eq*veqten(1:3,:,:),1)
          eq_vadv(2,:,:)=dtmt*SUM(ve_eq*veqten(4:6,:,:),1)
          eq_vadv(3,:,:)=dtmt*SUM(ve_eq*veqten(7:9,:,:),1)
        ENDIF
c-----------------------------------------------------------------------
c       putting the linear advection loop first helps reduce the number
c       of FFTs.
c
c       compute linear terms V_eq.grad(V) and V.grad(V_eq) for linear
c       cases or just the former for nonlinear cases, since grad(V_eq)
c       is part of real_vten.
c
c       ja is used as temporary space.
c-----------------------------------------------------------------------
        IF (eq_flow/='none') THEN
          DO imode=1,nmodes
            ja(1,:,:,imode)=SUM(ve_eq*vten(1:3,:,:,imode),1)
            ja(2,:,:,imode)=SUM(ve_eq*vten(4:6,:,:,imode),1)
            ja(3,:,:,imode)=SUM(ve_eq*vten(7:9,:,:,imode),1)
          ENDDO
          IF (nonlinear) THEN
            CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ja,
     $                   real_force,dealiase)
          ELSE
            DO imode=1,nmodes
              ja(1,:,:,imode)=dtmt*(ja(1,:,:,imode)
     $                   +SUM(ve(:,:,:,imode)*veqten(1:3,:,:),1))
              ja(2,:,:,imode)=dtmt*(ja(2,:,:,imode)
     $                   +SUM(ve(:,:,:,imode)*veqten(4:6,:,:),1))
              ja(3,:,:,imode)=dtmt*(ja(3,:,:,imode)
     $                   +SUM(ve(:,:,:,imode)*veqten(7:9,:,:),1))
            ENDDO
          ENDIF
        ELSE
          IF (nonlinear) real_force=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       nonlinear terms:  determine V in configuration space, then
c       find V.grad(V+V_eq)
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          real_ve=>real_ja
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $                 dealiase)
          IF (continuity=='full') THEN
            DO ip=1,nphi
              DO ix=1,mpseudo
                rscal=dtmt*real_ndptr(1,ix,ip)
                real_force(1,ix,ip)=rscal*( real_force(1,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(1,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(2,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(3,ix,ip))
                real_force(2,ix,ip)=rscal*( real_force(2,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(4,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(5,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(6,ix,ip))
                real_force(3,ix,ip)=rscal*( real_force(3,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(7,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(8,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(9,ix,ip))
              ENDDO
            ENDDO
          ELSE
            DO ip=1,nphi
              DO ix=1,mpseudo
                real_force(1,ix,ip)=dtmt*( real_force(1,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(1,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(2,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(3,ix,ip))
                real_force(2,ix,ip)=dtmt*( real_force(2,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(4,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(5,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(6,ix,ip))
                real_force(3,ix,ip)=dtmt*( real_force(3,ix,ip)+
     $                              real_ve(1,ix,ip)*real_vten(7,ix,ip)+
     $                              real_ve(2,ix,ip)*real_vten(8,ix,ip)+
     $                              real_ve(3,ix,ip)*real_vten(9,ix,ip))
              ENDDO
            ENDDO
          ENDIF
          IF (.NOT.ve_set) NULLIFY(real_vten)
c-----------------------------------------------------------------------
c         transform to Fourier space (ja is used as a temporary array).
c-----------------------------------------------------------------------
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,ja,
     $                 real_force,dealiase)
        ENDIF
c-----------------------------------------------------------------------
c       multiply by mass density here for cases where it doesn't give
c       further Fourier-component coupling.
c-----------------------------------------------------------------------
        IF (.NOT.nonlinear.OR.continuity=='none'.OR.
     $      continuity=='fix profile') THEN
          IF (.NOT.nonlinear.AND.continuity/='none') THEN
            DO im=1,nmodes
              force(1,:,:,im)=force(1,:,:,im)-ja(1,:,:,im)*nd_eq(1,:,:)
     $                        -eq_vadv(1,:,:)*nd(1,:,:,im)
              force(2,:,:,im)=force(2,:,:,im)-ja(2,:,:,im)*nd_eq(1,:,:)
     $                        -eq_vadv(2,:,:)*nd(1,:,:,im)
              force(3,:,:,im)=force(3,:,:,im)-ja(3,:,:,im)*nd_eq(1,:,:)
     $                        -eq_vadv(3,:,:)*nd(1,:,:,im)
            ENDDO
          ELSE
            DO im=1,nmodes
              force(1,:,:,im)=force(1,:,:,im)-ja(1,:,:,im)*nd_eq(1,:,:)
              force(2,:,:,im)=force(2,:,:,im)-ja(2,:,:,im)*nd_eq(1,:,:)
              force(3,:,:,im)=force(3,:,:,im)-ja(3,:,:,im)*nd_eq(1,:,:)
            ENDDO
          ENDIF
        ELSE IF (continuity=='n=0 only') THEN
          DO im=1,nmodes
            force(1,:,:,im)=force(1,:,:,im)-
     $         ja(1,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
            force(2,:,:,im)=force(2,:,:,im)-
     $         ja(2,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
            force(3,:,:,im)=force(3,:,:,im)-
     $         ja(3,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
            IF (keff(im)==0.AND.eq_flow/='none') THEN
             force(1,:,:,im)=force(1,:,:,im)-eq_vadv(1,:,:)*nd_n0(1,:,:)
             force(2,:,:,im)=force(2,:,:,im)-eq_vadv(2,:,:)*nd_n0(1,:,:)
             force(3,:,:,im)=force(3,:,:,im)-eq_vadv(3,:,:)*nd_n0(1,:,:)
            ENDIF
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       complete the advective term with 3D mass density.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.continuity=='full') THEN
          IF (eq_flow/='none') THEN
            DO imode=1,nmodes
              force(1,:,:,imode)=force(1,:,:,imode)-ja(1,:,:,imode)
     $                              -eq_vadv(1,:,:)*nd(1,:,:,imode)
              force(2,:,:,imode)=force(2,:,:,imode)-ja(2,:,:,imode)
     $                              -eq_vadv(2,:,:)*nd(1,:,:,imode)
              force(3,:,:,imode)=force(3,:,:,imode)-ja(3,:,:,imode)
     $                              -eq_vadv(3,:,:)*nd(1,:,:,imode)
            ENDDO
          ELSE
            force=force-ja
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       include the momentum-density correction for artificial particle
c       diffusivity.  again ja is just used as temporary space.
c-----------------------------------------------------------------------
        IF (continuity/='none'.AND.nd_correrr.AND.
     $      (nd_diff>0.OR.nd_hypd>0)) THEN
          CALL generic_ptr_set(rb%qndiff,tb%qndiff,tb%tgeom,inode,
     $                         ndiff,dcp,dcp,0_i4)
          IF (nonlinear) THEN
            CALL generic_ptr_set(rb%qndiff_phi,tb%qndiff_phi,tb%tgeom,
     $                           inode,real_ndiff,dp,dp,0_i4)
            DO ip=1,nphi
              DO ix=1,mpseudo
                real_force(1,ix,ip)=
     $            dtmt*real_ve(1,ix,ip)*real_ndiff(1,ix,ip)
                real_force(2,ix,ip)=
     $            dtmt*real_ve(2,ix,ip)*real_ndiff(1,ix,ip)
                real_force(3,ix,ip)=
     $            dtmt*real_ve(3,ix,ip)*real_ndiff(1,ix,ip)
              ENDDO
            ENDDO
            CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,
     $                   ja,real_force,dealiase)
            force=force-ja
          ENDIF
          IF (eq_flow/='none') THEN
            DO imode=1,nmodes
              force(1,:,:,imode)=force(1,:,:,imode)-
     $          dtmt*ve_eq(1,:,:)*ndiff(1,:,:,imode)
              force(2,:,:,imode)=force(2,:,:,imode)-
     $          dtmt*ve_eq(2,:,:)*ndiff(1,:,:,imode)
              force(3,:,:,imode)=force(3,:,:,imode)-
     $          dtmt*ve_eq(3,:,:)*ndiff(1,:,:,imode)
            ENDDO
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     evaluations for viscous dissipation:
c-----------------------------------------------------------------------
      IF ((kin_visc>0.OR.iso_visc>0.OR.par_visc>0.OR.gyr_visc>0).AND.
     $ (integrand_flag(1:4)=='visc'.OR.integrand_flag(1:3)=='all')) THEN
        IF (.NOT.ve_set) THEN
          CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,
     $                         ve,ve_r,ve_z,1_i4)
          CALL math_grad(nmodes,keff,3_i4,geom,ve,ve_r,ve_z,vten,bigr)
        ENDIF
        IF (ds_use=='kin_visc'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          diff=dtmt*MAX(ds2(1,:,:),0._r8)
        ELSE
          diff=dtmt
        ENDIF
        IF (eq_flow/='none')
     $    CALL generic_ptr_set(rb%qpi_veq,tb%qpi_veq,tb%tgeom,
     $                         inode,pi_veq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c       the grad(V) tensor is in real space and includes any Veq terms.
c       find B+B_eq in real space.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.(par_visc>0.OR.gyr_visc>0)) THEN
          CALL generic_ptr_set(rb%qbe_tot,tb%qbe_tot,tb%tgeom,
     $                         inode,real_bptr,dp,dp,0_i4)
          IF (.NOT.ve_set)
     $      CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                           inode,real_vten,dp,dp,0_i4)
        ELSE IF (nonlinear.AND.visc_heat.AND..NOT.ve_set) THEN
          CALL generic_ptr_set(rb%qgrdv,tb%qgrdv,tb%tgeom,
     $                         inode,real_vten,dp,dp,0_i4)
        ENDIF
c-----------------------------------------------------------------------
c       the present strategy for the viscous dissipation terms in
c       nonlinear computations is to find the full stress, including
c       all equilibrium terms, then subtracting the purely equilibrium
c       part.  here, the steady-state is not considerably
c       larger than the perturbed contribution, so there should be
c       no loss of accuracy.
c  
c       kinematic and isotropic stresses are modulated by the shape
c       function diff.
c
c       add grad(V_eq) to the grad(V) tensor.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.eq_flow/='none') THEN
          DO imode=1,nmodes
            IF (keff(imode)==0)
     $        vten(:,:,:,imode)=vten(:,:,:,imode)+veqten
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       stress tensor for kinetematic and isotropic viscosity--now
c       in an external subroutine.
c-----------------------------------------------------------------------
        IF (iso_visc>0.OR.kin_visc>0) THEN
          CALL iso_stress(ncx,ncy,nmodes,kin_visc,iso_visc,
     $                    vten,diff,piten)
        ELSE
          piten=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       transform piten from kinematic and isotropic contributions.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.(par_visc>0.OR.gyr_visc>0.OR.
     $                     continuity=='full')) THEN
          IF (iso_visc>0.OR.kin_visc>0) THEN
            CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,9_i4,
     $                   piten,real_pten,dealiase)
          ELSE
            real_pten=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       stress tensor for parallel viscosity
c
c       the steady-state term is subtracted after forming the 
c       complete stress. 
c
c       the stored equilibrium stress does not include the factor of dt.
c-----------------------------------------------------------------------
        IF (par_visc>0) THEN
          vcoef=-3*dtmt*par_visc
          IF (eq_flow/='none')
     $      CALL generic_ptr_set(rb%qpi_pareq,tb%qpi_pareq,tb%tgeom,
     $                           inode,pareq,dp,dp,0_i4)
c-----------------------------------------------------------------------
c         nonlinear stress--now from the external subroutine
c         par_stress_nl.
c
c         find (B.grad(V).B - B**2*div(V))/B**4 and multiply the result
c         by (BB - I*B**2/3).
c
c         the Braginskii temperature dependence is the same as the
c         parallel ion thermal diffusivity, so use that data for
c         that coefficient if specified.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            IF (parvisc_model=='plltdep') THEN
              CALL generic_ptr_set(rb%qkappli_phi,tb%qkappli_phi,
     $                             tb%tgeom,inode,kappli,dp,dp,0_i4)
              real_scal=vcoef*kappli/k_plli
            ELSE
              real_scal=vcoef
            ENDIF
            CALL par_stress_nl(mpseudo,nphi,smallnum,real_scal,
     $                         real_vten,real_bptr,real_pten)
c-----------------------------------------------------------------------
c         linear version from the external subroutine par_stress_lin.
c-----------------------------------------------------------------------
          ELSE
            CALL par_stress_lin(ncx,ncy,nmodes,vcoef,smallnum,
     $                          vten,be_eq,parten)
c-----------------------------------------------------------------------
c           V_eq contributions to the linear parallel stress from
c           the external subroutine par_stress_veq.
c-PRE
c           Braginskii parallel does not have density.
c-----------------------------------------------------------------------
            IF (eq_flow/='none') THEN
              CALL par_stress_veq(ncx,ncy,nmodes,vcoef,dt,smallnum,vten,
     $                            be,be_eq,veqten,pareq,dtmp,parten)
            ENDIF
            piten=piten+parten
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       Braginskii gyroviscosity is 0.5*eta_3*{bxW.(I+3bb)-(I+3bb).Wxb}
c       or
c
c       0.25*[n_e*m_i*Ti/(Z**2*B**4)]*{ [(BxW.(B**2+3BB)]+transp[] }
c
c       leave the number density for later.  note that the div(V) terms
c       cancel analytically and are left out of the W used here.
c
c       reuse the parten array to save a little memory; this
c       section is strictly gyroviscosity, however.
c
c       the stored equilibrium part does not include the factor of dt.
c-----------------------------------------------------------------------
        IF (gyr_visc>0) THEN
          vcoef=0.25_r8*gyr_visc*dt*ms(2)*kboltz/(zeff**2*elementary_q)
          IF (eq_flow/='none') THEN
            CALL generic_ptr_set(rb%qpi_gyreq,tb%qpi_gyreq,tb%tgeom,
     $                           inode,gyreq,dp,dp,0_i4)
          ENDIF
c-----------------------------------------------------------------------
c         nonlinear gyroviscous stress including grad(V_eq)
c         contributions if there are any.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            CALL generic_ptr_set(rb%qti_tot,tb%qti_tot,tb%tgeom,
     $                           inode,real_titot,dp,dp,0_i4)
            DO ip=1,nphi
              DO ix=1,mpseudo
                rb2=SUM(real_bptr(:,ix,ip)**2)+smallnum
                gcoef=vcoef*real_titot(1,ix,ip)/rb2**2
                vtmpr(1:3,1)=real_vten(1:3,ix,ip)+real_vten(1:7:3,ix,ip)
                vtmpr(1:3,2)=real_vten(4:6,ix,ip)+real_vten(2:8:3,ix,ip)
                vtmpr(1:3,3)=real_vten(7:9,ix,ip)+real_vten(3:9:3,ix,ip)
                btmp(1:3,1)=3._r8*real_bptr(1:3,ix,ip)*
     $                            real_bptr(1  ,ix,ip)
                btmp(1:3,2)=3._r8*real_bptr(1:3,ix,ip)*
     $                            real_bptr(2  ,ix,ip)
                btmp(1:3,3)=3._r8*real_bptr(1:3,ix,ip)*
     $                            real_bptr(3  ,ix,ip)
                btmp(1,1)=rb2+btmp(1,1)
                btmp(2,2)=rb2+btmp(2,2)
                btmp(3,3)=rb2+btmp(3,3)
                DO i2=1,3
                  wdotr(1,i2)=SUM(vtmpr(1,:)*btmp(:,i2))
                  wdotr(2,i2)=SUM(vtmpr(2,:)*btmp(:,i2))
                  wdotr(3,i2)=SUM(vtmpr(3,:)*btmp(:,i2))
                ENDDO
                bc_wdotr(1,:)=real_bptr(2,ix,ip)*wdotr(3,:)-
     $                        real_bptr(3,ix,ip)*wdotr(2,:)
                bc_wdotr(2,:)=real_bptr(3,ix,ip)*wdotr(1,:)-
     $                        real_bptr(1,ix,ip)*wdotr(3,:)
                bc_wdotr(3,:)=real_bptr(1,ix,ip)*wdotr(2,:)-
     $                        real_bptr(2,ix,ip)*wdotr(1,:)
                real_pten(1:3,ix,ip)=real_pten(1:3,ix,ip)+
     $                               gcoef*(bc_wdotr(:,1)+bc_wdotr(1,:))
                real_pten(4:6,ix,ip)=real_pten(4:6,ix,ip)+
     $                               gcoef*(bc_wdotr(:,2)+bc_wdotr(2,:))
                real_pten(7:9,ix,ip)=real_pten(7:9,ix,ip)+
     $                               gcoef*(bc_wdotr(:,3)+bc_wdotr(3,:))
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c         stress tensor for linear gyroviscosity.
c-----------------------------------------------------------------------
          ELSE
            IF (eq_flow=='none') THEN
              DO im=1,nmodes
                DO iy=1,ncy
                  DO ix=1,ncx
                    rb2=SUM(be_eq(:,ix,iy)**2)+smallnum
                    gcoef=vcoef*ti_eq(1,ix,iy)/rb2**2
                    vtmp(1:3,1)=vten(1:3,ix,iy,im)+vten(1:7:3,ix,iy,im)
                    vtmp(1:3,2)=vten(4:6,ix,iy,im)+vten(2:8:3,ix,iy,im)
                    vtmp(1:3,3)=vten(7:9,ix,iy,im)+vten(3:9:3,ix,iy,im)
                    btmp(1:3,1)=3._8*be_eq(1:3,ix,iy)*be_eq(1,ix,iy)
                    btmp(1:3,2)=3._8*be_eq(1:3,ix,iy)*be_eq(2,ix,iy)
                    btmp(1:3,3)=3._8*be_eq(1:3,ix,iy)*be_eq(3,ix,iy)
                    btmp(1,1)=rb2+btmp(1,1)
                    btmp(2,2)=rb2+btmp(2,2)
                    btmp(3,3)=rb2+btmp(3,3)
                    DO i2=1,3
                      wdot(1,i2)=SUM(vtmp(1,:)*btmp(:,i2))
                      wdot(2,i2)=SUM(vtmp(2,:)*btmp(:,i2))
                      wdot(3,i2)=SUM(vtmp(3,:)*btmp(:,i2))
                    ENDDO
                    b_cr_wdot(1,:)=be_eq(2,ix,iy)*wdot(3,:)-
     $                             be_eq(3,ix,iy)*wdot(2,:)
                    b_cr_wdot(2,:)=be_eq(3,ix,iy)*wdot(1,:)-
     $                             be_eq(1,ix,iy)*wdot(3,:)
                    b_cr_wdot(3,:)=be_eq(1,ix,iy)*wdot(2,:)-
     $                             be_eq(2,ix,iy)*wdot(1,:)
                    piten(1:3,ix,iy,im)=piten(1:3,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,1)+b_cr_wdot(1,:))
                    piten(4:6,ix,iy,im)=piten(4:6,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,2)+b_cr_wdot(2,:))
                    piten(7:9,ix,iy,im)=piten(7:9,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,3)+b_cr_wdot(3,:))
                  ENDDO
                ENDDO
              ENDDO
c-----------------------------------------------------------------------
c           extra linear gyroviscous terms for equilibrium flow--not
c           for those with a weak stomach.  like the related parallel
c           stress, there are three terms here: 1) b-tilde in the cross
c           product with Weq, 2) perturbation of B**2+3BB, and 3)
c           perturbation of the Ti/B**4 coefficient.
c-----------------------------------------------------------------------
            ELSE
              DO im=1,nmodes
                DO iy=1,ncy
                  DO ix=1,ncx
                    rb2=SUM(be_eq(:,ix,iy)**2)+smallnum
                    gcoef=vcoef*ti_eq(1,ix,iy)/rb2**2
                    vtmp(1:3,1)=vten(1:3,ix,iy,im)+vten(1:7:3,ix,iy,im)
                    vtmp(1:3,2)=vten(4:6,ix,iy,im)+vten(2:8:3,ix,iy,im)
                    vtmp(1:3,3)=vten(7:9,ix,iy,im)+vten(3:9:3,ix,iy,im)
                    vtmpr(1:3,1)=veqten(1:3,ix,iy)+veqten(1:7:3,ix,iy)
                    vtmpr(1:3,2)=veqten(4:6,ix,iy)+veqten(2:8:3,ix,iy)
                    vtmpr(1:3,3)=veqten(7:9,ix,iy)+veqten(3:9:3,ix,iy)

                    btmp(1:3,1)=3._8*be_eq(1:3,ix,iy)*be_eq(1,ix,iy)
                    btmp(1:3,2)=3._8*be_eq(1:3,ix,iy)*be_eq(2,ix,iy)
                    btmp(1:3,3)=3._8*be_eq(1:3,ix,iy)*be_eq(3,ix,iy)
                    btmp(1,1)=rb2+btmp(1,1)
                    btmp(2,2)=rb2+btmp(2,2)
                    btmp(3,3)=rb2+btmp(3,3)
                    divtmp=2._r8*SUM(be_eq(:,ix,iy)*be(:,ix,iy,im))
                    btmpc(1:3,1)=3._8*(be_eq(1:3,ix,iy)*be(1,ix,iy,im)+
     $                                 be(1:3,ix,iy,im)*be_eq(1,ix,iy))
                    btmpc(1:3,2)=3._8*(be_eq(1:3,ix,iy)*be(2,ix,iy,im)+
     $                                 be(1:3,ix,iy,im)*be_eq(2,ix,iy))
                    btmpc(1:3,3)=3._8*(be_eq(1:3,ix,iy)*be(3,ix,iy,im)+
     $                                 be(1:3,ix,iy,im)*be_eq(3,ix,iy))
                    btmpc(1,1)=divtmp+btmpc(1,1)
                    btmpc(2,2)=divtmp+btmpc(2,2)
                    btmpc(3,3)=divtmp+btmpc(3,3)

                    DO i2=1,3
                      wdot(1,i2)=SUM(vtmp (1,:)*btmp (:,i2)+
     $                               vtmpr(1,:)*btmpc(:,i2))
                      wdot(2,i2)=SUM(vtmp (2,:)*btmp (:,i2)+
     $                               vtmpr(2,:)*btmpc(:,i2))
                      wdot(3,i2)=SUM(vtmp (3,:)*btmp (:,i2)+
     $                               vtmpr(3,:)*btmpc(:,i2))
                      wdotr(1,i2)=SUM(vtmpr(1,:)*btmp(:,i2))
                      wdotr(2,i2)=SUM(vtmpr(2,:)*btmp(:,i2))
                      wdotr(3,i2)=SUM(vtmpr(3,:)*btmp(:,i2))
                    ENDDO
                    b_cr_wdot(1,:)=be_eq(2,ix,iy)*wdot (3,:)-
     $                             be_eq(3,ix,iy)*wdot (2,:)+
     $                             be(2,ix,iy,im)*wdotr(3,:)-
     $                             be(3,ix,iy,im)*wdotr(2,:)
                    b_cr_wdot(2,:)=be_eq(3,ix,iy)*wdot (1,:)-
     $                             be_eq(1,ix,iy)*wdot (3,:)+
     $                             be(3,ix,iy,im)*wdotr(1,:)-
     $                             be(1,ix,iy,im)*wdotr(3,:)
                    b_cr_wdot(3,:)=be_eq(1,ix,iy)*wdot (2,:)-
     $                             be_eq(2,ix,iy)*wdot (1,:)+
     $                             be(1,ix,iy,im)*wdotr(2,:)-
     $                             be(2,ix,iy,im)*wdotr(1,:)

                    cscal=dt*(tion(1,ix,iy,im)/ti_eq(1,ix,iy)-
     $                        2._r8*divtmp/rb2)
                    piten(1:3,ix,iy,im)=piten(1:3,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,1)+b_cr_wdot(1,:))+
     $                cscal*gyreq(1:3,ix,iy)
                    piten(4:6,ix,iy,im)=piten(4:6,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,2)+b_cr_wdot(2,:))+
     $                cscal*gyreq(4:6,ix,iy)
                    piten(7:9,ix,iy,im)=piten(7:9,ix,iy,im)+
     $                gcoef*(b_cr_wdot(:,3)+b_cr_wdot(3,:))+
     $                cscal*gyreq(7:9,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       now multiply by number density.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.continuity=='full') THEN
          DO ip=1,nphi
            DO ix=1,mpseudo
              real_pten(:,ix,ip)=real_pten(:,ix,ip)*real_ndptr(1,ix,ip)
            ENDDO
          ENDDO
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,9_i4,
     $                 piten,real_pten,dealiase)
        ELSE IF (nonlinear.AND.(par_visc>0.OR.gyr_visc>0)) THEN
          CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,9_i4,
     $                 piten,real_pten,dealiase)
        ENDIF

        IF (nonlinear.AND.continuity=='n=0 only') THEN
          DO imode=1,nmodes
            DO iy=1,ncy
              DO ix=1,ncx
                piten(:,ix,iy,imode)=piten(:,ix,iy,imode)*
     $                               (nd_eq(1,ix,iy)+nd_n0(1,ix,iy))
              ENDDO
            ENDDO
          ENDDO
        ELSE IF (.NOT.nonlinear.OR.continuity/='full') THEN
          DO imode=1,nmodes
            DO iy=1,ncy
              DO ix=1,ncx
                piten(:,ix,iy,imode)=piten(:,ix,iy,imode)*nd_eq(1,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       now eliminate the n_eq*Pi(V_eq) term from the nonlinear stress
c       and add n_tilde*Pi(V0) to the linear stress if needed.  the
c       stored equilibrium stress now includes nd_eq.
c-----------------------------------------------------------------------
        IF (eq_flow/='none') THEN
          IF (nonlinear) THEN
            DO imode=1,nmodes
              IF (keff(imode)/=0) CYCLE
              DO iy=1,ncy
                DO ix=1,ncx
                  piten(:,ix,iy,imode)=piten(:,ix,iy,imode)-
     $              dt*pi_veq(:,ix,iy)
                ENDDO
              ENDDO
              IF (par_visc>0) THEN
                DO iy=1,ncy
                  DO ix=1,ncx
                    piten(:,ix,iy,imode)=piten(:,ix,iy,imode)-
     $                dt*pareq(:,ix,iy)
                  ENDDO
                ENDDO
              ENDIF
              IF (gyr_visc>0) THEN
                DO iy=1,ncy
                  DO ix=1,ncx
                    piten(:,ix,iy,imode)=piten(:,ix,iy,imode)-
     $                dt*gyreq(:,ix,iy)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE IF (continuity=='full') THEN
            DO imode=1,nmodes
              DO iy=1,ncy
                DO ix=1,ncx
                  piten(:,ix,iy,imode)=piten(:,ix,iy,imode)+
     $              dt*nd(1,ix,iy,imode)*pi_veq(:,ix,iy)
                ENDDO
              ENDDO
              IF (par_visc>0) THEN
                DO iy=1,ncy
                  DO ix=1,ncx
                    piten(:,ix,iy,imode)=piten(:,ix,iy,imode)+
     $                dt*nd(1,ix,iy,imode)*pareq(:,ix,iy)
                  ENDDO
                ENDDO
              ENDIF
              IF (gyr_visc>0) THEN
                DO iy=1,ncy
                  DO ix=1,ncx
                    piten(:,ix,iy,imode)=piten(:,ix,iy,imode)+
     $                dt*nd(1,ix,iy,imode)*gyreq(:,ix,iy)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       add total pressure to the diagonal of the stress tensor.
c-----------------------------------------------------------------------
        IF (beta>0.AND.(integrand_flag(1:3)=='mhd'.OR.
     $                  integrand_flag(1:3)=='all')) THEN
          DO imode=1,nmodes
            DO iy=1,ncy
              DO ix=1,ncx
                piten(1,ix,iy,imode)=piten(1,ix,iy,imode)+
     $                                pres(1,ix,iy,imode)
                piten(5,ix,iy,imode)=piten(5,ix,iy,imode)+
     $                                pres(1,ix,iy,imode)
                piten(9,ix,iy,imode)=piten(9,ix,iy,imode)+
     $                                pres(1,ix,iy,imode)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       add conjg(grad(alpha))^T : Pi(V)  transfer phi
c       derivatives to force before scattering contributions to the
c       nodes.
c-----------------------------------------------------------------------
        DO imode=1,nmodes
          IF (geom=='tor') THEN
            force(1,:,:,imode)=force(1,:,:,imode)
     $        +(                  piten(9,:,:,imode)
     $         -(0,1)*keff(imode)*piten(3,:,:,imode))/bigr
            force(2,:,:,imode)=force(2,:,:,imode)
     $         -(0,1)*keff(imode)*piten(6,:,:,imode) /bigr
            force(3,:,:,imode)=force(3,:,:,imode)
     $        -(                  piten(3,:,:,imode)
     $         +(0,1)*keff(imode)*piten(9,:,:,imode))/bigr
          ELSE
            force(1,:,:,imode)=force(1,:,:,imode)
     $         -(0,1)*keff(imode)*piten(3,:,:,imode)
            force(2,:,:,imode)=force(2,:,:,imode)
     $         -(0,1)*keff(imode)*piten(6,:,:,imode)
            force(3,:,:,imode)=force(3,:,:,imode)
     $         -(0,1)*keff(imode)*piten(9,:,:,imode)
          ENDIF
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(  alpha(:,:,iv)*force(1,:,:,imode)
     $                            +dalpdr(:,:,iv)*piten(1,:,:,imode)
     $                            +dalpdz(:,:,iv)*piten(2,:,:,imode),1)
            int(2,:,iv,imode)=SUM(  alpha(:,:,iv)*force(2,:,:,imode)
     $                            +dalpdr(:,:,iv)*piten(4,:,:,imode)
     $                            +dalpdz(:,:,iv)*piten(5,:,:,imode),1)
            int(3,:,iv,imode)=SUM(  alpha(:,:,iv)*force(3,:,:,imode)
     $                            +dalpdr(:,:,iv)*piten(7,:,:,imode)
     $                            +dalpdz(:,:,iv)*piten(8,:,:,imode),1)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     or just accumulate the force vertex contibutions.
c-----------------------------------------------------------------------
      ELSE
        DO imode=1,nmodes
          force(3,:,:,imode)=force(3,:,:,imode)
     $       -(0,1)*keff(imode)*pres(1,:,:,imode)/bigr
          DO iv=1,nvc
            int(1,:,iv,imode)=SUM(  alpha(:,:,iv)*force(1,:,:,imode)+
     $                            dalpdrc(:,:,iv)* pres(1,:,:,imode),1)
            int(2,:,iv,imode)=SUM(  alpha(:,:,iv)*force(2,:,:,imode)+
     $                             dalpdz(:,:,iv)* pres(1,:,:,imode),1)
            int(3,:,iv,imode)=SUM(alpha(:,:,iv)*force(3,:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if auxiliary fields are used, find their contribution to the
c     velocity bases and the rhs of its own evolution equation.  the
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
c     auxiliary field auxv.  pres and teff_eq are used as temporary
c     storage.
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
c     the discontinuous scalar field auxw, f is the centering
c     coefficient.  here, the effective viscosity is dt*v_Alven**2.
c
c     teff is used for temporary space, and ja becomes
c     dt*cpvrt*curl(V)/sqrt(mu0).
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
      IF (poly_divv>=0) THEN
        CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                       inode,peq,dp,dp,0_i4)
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                         inode,be0,be0r,be0z,1_i4)
          bes=be_eq+be0
          besr=beqr+be0r
          besz=beqz+be0z
          IF (beta>0) THEN
            CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                           inode,p0,dp,dp,0_i4)
            teff_eq(1,:,:)=SQRT(MAX(0._r8,gamma*(peq(1,:,:)+p0(1,:,:)))+
     $                          SUM(bes**2,1)/mu0)
          ELSE
            teff_eq(1,:,:)=SQRT(SUM(bes**2,1)/mu0)
          ENDIF
        ELSE
          bes=be_eq
          besr=beqr
          besz=beqz
          IF (beta>0) THEN
            teff_eq(1,:,:)=SQRT(MAX(0._r8,gamma*peq(1,:,:))+
     $                          SUM(bes**2,1)/mu0)
          ELSE
            teff_eq(1,:,:)=SQRT(SUM(bes**2,1)/mu0)
          ENDIF
        ENDIF

        dtdv=SQRT(ddivv/fdivv)*dt        !  for diffusive correction
        dtpv=dt*SQRT(dpvrt/(mu0*fpvrt))  !  for diffusive correction
        CALL math_curl(nmodes,keff,geom,bigr,ve,ve_r,ve_z,ja,dtpv)

        DO imode=1,nmodes
          IF (geom=='tor') THEN
            pres(1,:,:,imode)=(ve_r(1,:,:,imode)+ve_z(2,:,:,imode)
     $                         +((0,1)*keff(imode)*ve(3,:,:,imode)
     $                         +ve(1,:,:,imode))/bigr)*teff_eq(1,:,:)*
     $                          dtdv
          ELSE
            pres(1,:,:,imode)=(ve_r(1,:,:,imode)+ve_z(2,:,:,imode)
     $                         +(0,1)*keff(imode)*ve(3,:,:,imode)/bigr)
     $                        *teff_eq(1,:,:)*dtdv
          ENDIF
          teff(1,:,:,imode)=SUM(bes*ja(:,:,:,imode),1)

          DO iv=nvc+1,nv  !  auxiliary equation
            int(1,:,iv,imode)=SUM(almod(:,:,iv-nvc)*pres(1,:,:,imode),1)
            int(2,:,iv,imode)=
     $        SUM(almod(:,:,iv-nvc)*teff(1,:,:,imode),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vrhs
c-----------------------------------------------------------------------
c     subprogram 10. advect_cor.
c     this integrand is just the correction needed on the rhs
c     for iterating advection in a nonlinear implicit advance of V.
c-----------------------------------------------------------------------
      SUBROUTINE advect_cor(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          nd_eq,nd_n0,
     $          real_gradv,real_ndptr,dp
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_ve,real_force
      REAL(r8), DIMENSION(9,mpseudo,nphi) :: real_vten
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: dtmt,rscal

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             ve,ve_r,ve_z,nd,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) ::
     $             force
      COMPLEX(r8), DIMENSION(9,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: vten
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,imode,iq,nvc,ncx,ncy,ix,iy,ip,im
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      dtmt=0.25_r8*dt*mtot
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the perturbed
c     and 0-th order magnetic fields.  the current density is evaluated
c     from the gradients of the perturbed magnetic field
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nvc=SIZE(alpha,3)
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,
     $                     nd,dcp,dcp,0_i4)
      IF (continuity=='n=0 only') THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                       nd_n0,dp,dp,0_i4)
      ENDIF
      IF (continuity/='none') THEN
        CALL generic_ptr_set(rb%qnd_tot,tb%qnd_tot,tb%tgeom,
     $                       inode,real_ndptr,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     find 0.25*alpha_i*rho*dV.grad(dV) for nonlinear corrector steps
c     for implicit advection.  the dV vector is the difference between
c     the last nonlinear iteration and the beginning of the time-step,
c     and it is stored in work1.
c
c     the vten grad(V) tensor is in Fourier representation and is left
c     without equilibrium contributions.  real_vten is in real space
c     with equilibrium contributions.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qwork1,tb%qwork1,tb%tgeom,
     $                     inode,ve,ve_r,ve_z,1_i4)
      CALL math_grad(nmodes,keff,3_i4,geom,ve,ve_r,ve_z,vten,bigr)
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,9_i4,vten,
     $             real_vten,dealiase)
c-----------------------------------------------------------------------
c     determine dV in configuration space, then find dV.grad(dV).
c-----------------------------------------------------------------------
      CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_ve,
     $             dealiase)
      IF (continuity=='full') THEN
        DO ip=1,nphi
          DO ix=1,mpseudo
            rscal=dtmt*real_ndptr(1,ix,ip)
            real_force(1,ix,ip)=rscal*(
     $                          real_ve(1,ix,ip)*real_vten(1,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(2,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(3,ix,ip))
            real_force(2,ix,ip)=rscal*(
     $                          real_ve(1,ix,ip)*real_vten(4,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(5,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(6,ix,ip))
            real_force(3,ix,ip)=rscal*(
     $                          real_ve(1,ix,ip)*real_vten(7,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(8,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(9,ix,ip))
          ENDDO
        ENDDO
      ELSE
        DO ip=1,nphi
          DO ix=1,mpseudo
            real_force(1,ix,ip)=dtmt*(
     $                          real_ve(1,ix,ip)*real_vten(1,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(2,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(3,ix,ip))
            real_force(2,ix,ip)=dtmt*(
     $                          real_ve(1,ix,ip)*real_vten(4,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(5,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(6,ix,ip))
            real_force(3,ix,ip)=dtmt*(
     $                          real_ve(1,ix,ip)*real_vten(7,ix,ip)+
     $                          real_ve(2,ix,ip)*real_vten(8,ix,ip)+
     $                          real_ve(3,ix,ip)*real_vten(9,ix,ip))
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     transform to Fourier space.
c-----------------------------------------------------------------------
      CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,force,
     $             real_force,dealiase)
c-----------------------------------------------------------------------
c     multiply by mass density here for cases where it doesn't
c     involve Fourier-component coupling.
c-----------------------------------------------------------------------
      IF (continuity=='none'.OR.continuity=='fix profile') THEN
        DO im=1,nmodes
          force(1,:,:,im)=force(1,:,:,im)*nd_eq(1,:,:)
          force(2,:,:,im)=force(2,:,:,im)*nd_eq(1,:,:)
          force(3,:,:,im)=force(3,:,:,im)*nd_eq(1,:,:)
        ENDDO
      ELSE IF (continuity=='n=0 only') THEN
        DO im=1,nmodes
          force(1,:,:,im)=force(1,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
          force(2,:,:,im)=force(2,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
          force(3,:,:,im)=force(3,:,:,im)*(nd_eq(1,:,:)+nd_n0(1,:,:))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     accumulate the force vertex contibutions.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        DO iv=1,nvc
          int(1,:,iv,imode)=SUM(alpha(:,:,iv)*force(1,:,:,imode),1)
          int(2,:,iv,imode)=SUM(alpha(:,:,iv)*force(2,:,:,imode),1)
          int(3,:,iv,imode)=SUM(alpha(:,:,iv)*force(3,:,:,imode),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advect_cor
c-----------------------------------------------------------------------
c     subprogram 11. divb_rhs.
c     compute grad_div(b) for the rhs of the div(b) diffuser.
c-----------------------------------------------------------------------
      SUBROUTINE divb_rhs(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: be,ber,bez
      COMPLEX(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: divb

      INTEGER(i4) :: iv,imode,iq,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the perturbed
c     b for each mode.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      CALL generic_ptr_set(rb%qbe,tb%qbe,tb%tgeom,inode,be,ber,bez,1_i4)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     begin mode loop, and find the real and imaginary parts of div(b)
c     for this mode.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
        IF (geom=='tor') THEN
          divb=ber(1,:,:,imode)+bez(2,:,:,imode)
     $            +((0,1)*keff(imode)*be(3,:,:,imode)
     $                               +be(1,:,:,imode))/bigr
        ELSE
          divb=ber(1,:,:,imode)+bez(2,:,:,imode)
     $             +(0,1)*keff(imode)*be(3,:,:,imode)/bigr
        ENDIF
c-----------------------------------------------------------------------
c       evaluate -div(alpha)*div(B).
c-----------------------------------------------------------------------
        DO iv=1,nv
          int(1,:,iv,imode)=SUM(
     $           -divbd*dt*dalpdrc(:,:,iv)*divb,1)
          int(2,:,iv,imode)=SUM(
     $           -divbd*dt*dalpdz (:,:,iv)*divb,1)
          int(3,:,iv,imode)=SUM(
     $           +(0,1)*divbd*dt*alpha(:,:,iv)*keff(imode)/bigr*divb,1)
        ENDDO
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE divb_rhs
c-----------------------------------------------------------------------
c     subprogram 12. ndrhs.
c     compute the integrand used in the rhs of the continuity equation.
c-----------------------------------------------------------------------
      SUBROUTINE ndrhs(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          nd_eq,ve_eq,dart,upwc,
     $                                       vetot,ndeqr,ndeqz,dp
      REAL(r8), DIMENSION(3,mpseudo,nphi) :: real_vec,real_gn
      REAL(r8), DIMENSION(1,mpseudo,nphi) :: real_scal,real_n
      REAL(r8), DIMENSION(  mpseudo) :: rdot,rscl
      REAL(r8), DIMENSION(SIZE(bigr,1)*SIZE(bigr,2)) :: ndeqr_1d,
     $          ndeqz_1d

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             ve,nd,ndr,ndz,dcp
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),nmodes) :: 
     $             nvel,grad_nd

      INTEGER(i4) :: nv,iv,imode,iq,jq,ncx,ncy,ip
      REAL(r8) :: fupw,fhyp
c-----------------------------------------------------------------------
c     convenience parameters:
c-----------------------------------------------------------------------
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nv=SIZE(int,3)
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate velocity.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      CALL generic_ptr_set(rb%qve,tb%qve,tb%tgeom,inode,ve,dcp,dcp,0_i4)
      IF (eq_flow/='none') THEN
        CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,ve_eq,dp,dp,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     density
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,ndeqr,ndeqz,1_i4)
      CALL generic_ptr_set(rb%qnd,tb%qnd,tb%tgeom,inode,
     $                     nd,ndr,ndz,1_i4)
      IF (nd_diff>0.OR.nonlinear.AND.nd_dart_upw>0.AND.impladv.OR.
     $    nd_hypd>0.AND.impladv) THEN
        DO imode=1,nmodes
          grad_nd(1,:,:,imode)=ndr(1,:,:,imode)
          grad_nd(2,:,:,imode)=ndz(1,:,:,imode)
          grad_nd(3,:,:,imode)=(0,1)*keff(imode)*nd(1,:,:,imode)/bigr
        ENDDO
        fhyp=SQRT(nd_hypd*dt)
      ENDIF
      IF (nonlinear.AND.nd_diff>0.AND.nd_floor>0)
     $  CALL generic_ptr_set(rb%qdart,tb%qdart,tb%tgeom,
     $                       inode,dart,dp,dp,0_i4)
      IF (nonlinear.AND.nd_dart_upw>0.AND.impladv) THEN
        CALL generic_ptr_set(rb%qupw_phi,tb%qupw_phi,tb%tgeom,
     $                       inode,upwc,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qve_tot,tb%qve_tot,tb%tgeom,
     $                       inode,vetot,dp,dp,0_i4)
        ndeqr_1d=RESHAPE(ndeqr,(/ncx*ncy/))
        ndeqz_1d=RESHAPE(ndeqz,(/ncx*ncy/))
        fupw=nd_dart_upw*dt**2
      ENDIF
      IF ((nonlinear.OR.eq_flow/='none').AND.
     $    integrand_flag/='n predict'.AND..NOT.impladv) THEN
        CALL generic_ptr_set(rb%qwork3,tb%qwork3,tb%tgeom,inode,
     $                       nd,ndr,ndz,1_i4)
      ENDIF
c-----------------------------------------------------------------------
c     find the nonlinear contributions to n*V.
c     note that the real arrays have three dimensions, where the
c     first is the vector component, the second covers the poloidal
c     plane, and the third is the phi index,
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
c-----------------------------------------------------------------------
c       transform v_perturbed and n_perturbed to real space.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,ve,real_vec,
     $               dealiase)
        CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,1_i4,nd,real_n,
     $               dealiase)
c-----------------------------------------------------------------------
c       apply an upwinding-like artificial diffusion with the
c       pre-determined diffusivity multiplied by VV.grad(n) and add
c       this to n*V.
c-----------------------------------------------------------------------
        IF (nd_dart_upw>0.AND.impladv) THEN
          CALL fft_nim('inverse',ncx*ncy,mpseudo,lphi,3_i4,grad_nd,
     $                 real_gn,dealiase)
          DO ip=1,nphi
            rdot=fupw*upwc(1,:,ip)*upw_aniso * ( 
     $          SUM(vetot(:,:,ip)*real_gn(:,:,ip),1) +
     $              real_vec(1,:,ip)*ndeqr_1d(ipseust:ipseuen) +
     $              real_vec(2,:,ip)*ndeqz_1d(ipseust:ipseuen) ) /
     $        ( SUM(vetot(:,:,ip)*vetot(:,:,ip),1) + smallnum )
            rscl=fupw*upwc(1,:,ip)*(1._r8-upw_aniso)
            real_vec(1,:,ip)=dt*real_vec(1,:,ip)*real_n(1,:,ip)-
     $        ( rdot*vetot(1,:,ip)+
     $          rscl*(real_gn(1,:,ip)+ndeqr_1d(ipseust:ipseuen)) )
            real_vec(2,:,ip)=dt*real_vec(2,:,ip)*real_n(1,:,ip)-
     $        ( rdot*vetot(2,:,ip)+
     $          rscl*(real_gn(2,:,ip)+ndeqz_1d(ipseust:ipseuen)) )
            real_vec(3,:,ip)=dt*real_vec(3,:,ip)*real_n(1,:,ip)-
     $        ( rdot*vetot(3,:,ip)+rscl*real_gn(3,:,ip) )
          ENDDO
c-----------------------------------------------------------------------
c       or just find n*V
c-----------------------------------------------------------------------
        ELSE
          real_vec(1,:,:)=dt*real_vec(1,:,:)*real_n(1,:,:)
          real_vec(2,:,:)=dt*real_vec(2,:,:)*real_n(1,:,:)
          real_vec(3,:,:)=dt*real_vec(3,:,:)*real_n(1,:,:)
        ENDIF
c-----------------------------------------------------------------------
c       transform flux to Fourier space
c-----------------------------------------------------------------------
        CALL fft_nim('forward',ncx*ncy,mpseudo,lphi,3_i4,nvel,real_vec,
     $               dealiase)
      ELSE
        nvel=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     begin mode loop.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       find the linear n*V for this mode.
c-----------------------------------------------------------------------
        IF (eq_flow/='none') THEN
          DO iv=1,3
            nvel(iv,:,:,imode)=nvel(iv,:,:,imode)
     $        +dt*(ve(iv,:,:,imode)*nd_eq(1,:,:)
     $            +ve_eq(iv,:,:)*nd(1,:,:,imode))
          ENDDO
        ELSE
          DO iv=1,3
            nvel(iv,:,:,imode)=nvel(iv,:,:,imode)
     $        +dt*ve(iv,:,:,imode)*nd_eq(1,:,:)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       add diffusive term if needed.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.nd_floor>0.AND.nd_diff>0) THEN
          DO iv=1,3
            nvel(iv,:,:,imode)=nvel(iv,:,:,imode)-
     $        (dt*nd_diff+dart(1,:,:))*grad_nd(iv,:,:,imode)
          ENDDO
        ELSE IF (nd_diff>0) THEN
          nvel(:,:,:,imode)=nvel(:,:,:,imode)-
     $         dt*nd_diff*grad_nd(:,:,:,imode)
        ENDIF
c-----------------------------------------------------------------------
c       construct grad(alpha).( n*V )
c       with linear and nonlinear contributions, possibly including
c       the rhs for the hyper diffusivity auxiliary variable 
c       that is proportional to grad**2(n).
c-----------------------------------------------------------------------
        IF (nd_hypd>0.AND.impladv) THEN
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(
     $        +dalpdr(:,:,iv)*nvel(1,:,:,imode)
     $        +dalpdz(:,:,iv)*nvel(2,:,:,imode)
     $        -(0,1)*keff(imode)*alpha(:,:,iv)*nvel(3,:,:,imode)/bigr,1)
            int(2,:,iv,imode)=SUM(-fhyp*(
     $        dalpdr(:,:,iv)*grad_nd(1,:,:,imode)+
     $        dalpdz(:,:,iv)*grad_nd(2,:,:,imode)-
     $        (0,1)*keff(imode)*alpha(:,:,iv)*
     $                       grad_nd(3,:,:,imode)/bigr),1)
          ENDDO
        ELSE
          DO iv=1,nv
            int(1,:,iv,imode)=SUM(
     $        +dalpdr(:,:,iv)*nvel(1,:,:,imode)
     $        +dalpdz(:,:,iv)*nvel(2,:,:,imode)
     $        -(0,1)*keff(imode)*alpha(:,:,iv)*nvel(3,:,:,imode)/bigr,1)
          ENDDO
        ENDIF
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ndrhs

c-----------------------------------------------------------------------
c     subprogram 13. scal_lapl_rhs.
c     find the integrand for a projection of the Laplacian operator
c     acting on the scalar in qwork3.
c-----------------------------------------------------------------------
      SUBROUTINE scal_lapl_rhs(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             scal,scalr,scalz
      INTEGER(i4) :: iv,imode
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and evaluate velocity.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'rhs',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     scalar field:
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qwork3,tb%qwork3,tb%tgeom,inode,scal,
     $                     scalr,scalz,1_i4)
c-----------------------------------------------------------------------
c     begin mode loop and construct -grad(alpha^*).grad(scal).
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
        DO iv=1,SIZE(int,3)
          int(1,:,iv,imode)=-SUM(
     $       dalpdr(:,:,iv)*scalr(1,:,:,imode)
     $      +dalpdz(:,:,iv)*scalz(1,:,:,imode)
     $      +k2ef(imode)*alpha(:,:,iv)*scal(1,:,:,imode)/bigr**2,1)
        ENDDO
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE scal_lapl_rhs
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE integrands_rhs
