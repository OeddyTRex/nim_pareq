c-----------------------------------------------------------------------
c     file integrands.f
c     module that includes integrand routines used to create matrices
c     for linear algebraic operations over the finite-element mesh.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module integrands_mat
c     subprogram 1. get_mass.
c     subprogram 2. curl_de_iso.
c     subprogram 3. curl_de_ciso.
c     subprogram 4. b_hmhd_op.
c     subprogram 5. b_hyp_op.
c     subprogram 6. n_iso_op.
c     subprogram 7. t_aniso_op.
c     subprogram 8. vec_lap_op.
c     subprogram 9. v_aniso_op.
c     subprogram 10. grad_div.
c     subprogram 11. cont_op.
c     subprogram 12. vec_lap_op2.
c     subprogram 13. hyp_visc_op.
c-----------------------------------------------------------------------
c     module containing integrands for matrices on the lhs of equations.
c-----------------------------------------------------------------------
      MODULE integrands_mat

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
c
c     with separate polynomial degree specification for each field,
c     this is only used for the jfromb computation.
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
      REAL(r8), DIMENSION(1,1,1) :: dv
      INTEGER(i4) :: iv,jv,nv,iq,nq
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      nq=SIZE(int,1)
      nv=SIZE(int,5)
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
c-----------------------------------------------------------------------
c     find alpha(iv)*alpha(jv)
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(:,:,:,jv,iv)=0._r8
          DO iq=1,nq
            int(iq,iq,:,jv,iv)=SUM(alpha(:,:,iv)*alpha(:,:,jv),1)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE get_mass
c-----------------------------------------------------------------------
c     subprogram 2. curl_de_iso.
c     compute the integrand for the linearized matrix on the lhs of the 
c     magnetic field equation, where the effective impedance tensor
c     is isotropic.  This is now used for mhd computations only.
c-----------------------------------------------------------------------
      SUBROUTINE curl_de_iso(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: elecd_n0,dp,ds2
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: ziso,bigr2,diff
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc

      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      IF (elecd>0) THEN
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                         inode,elecd_n0,dp,dp,0_i4)
          diff=elecd_n0(1,:,:)
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          diff=elecd*MAX(ds2(1,:,:),0._r8)
        ELSE
          diff=elecd
        ENDIF
      ENDIF
      nv=SIZE(int,5)
c-----------------------------------------------------------------------
c     find R**2.
c-----------------------------------------------------------------------
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     evaluate coefficients for the impedance matrix--now just
c     the implicit contribution to resistive dissipation, where the
c     electrical diffusivity is modulated by the shape function diff.
c-----------------------------------------------------------------------
      IF (elecd>0) THEN
        ziso=dt*feta*diff
      ELSE
        ziso=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     compute the eta/mu0Xcurl(alpha_iv).curl(alpha_jv) operator.
c     Components are cylindrical.
c     note:  vector storage is (real1,real2,-imag3) or
c     (imag1,imag2,real3) and jmode is set in the calling routine.
c
c     if div(b) cleaning is not time split, add the grad(div) operator,
c     too.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb) THEN
        DO iv=1,nv
          DO jv=1,nv
            int(1,1,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                          ziso*(dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $              +  alpha(:,:,iv)*  alpha(:,:,jv)*k2ef(jmode)/bigr2)
     $              + dt*fdivb*divbd*dalpdrc(:,:,iv)*dalpdrc(:,:,jv),1)

            int(2,1,:,jv,iv)=SUM(
     $                         -ziso* dalpdz(:,:,iv)*dalpdr (:,:,jv)
     $              + dt*fdivb*divbd*dalpdrc(:,:,iv)*dalpdz (:,:,jv),1)

            int(1,2,:,iv,jv)=int(2,1,:,jv,iv)
            int(3,1,:,jv,iv)=SUM(
     $                          ziso*  alpha(:,:,iv)*keff(jmode)/bigr
     $                              *dalpdrc(:,:,jv)
     $              + dt*fdivb*divbd*dalpdrc(:,:,iv)
     $                              * alpha (:,:,jv)*keff(jmode)/bigr,1)

            int(1,3,:,iv,jv)=int(3,1,:,jv,iv)

            int(2,2,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                          ziso*(dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $              +  alpha(:,:,iv)*  alpha(:,:,jv)*k2ef(jmode)/bigr2)
     $              + dt*fdivb*divbd*dalpdz (:,:,iv)*dalpdz(:,:,jv),1)

            int(3,2,:,jv,iv)=SUM(
     $                          ziso*  alpha(:,:,iv)*keff(jmode)/bigr
     $                              * dalpdz(:,:,jv)
     $              + dt*fdivb*divbd*dalpdz (:,:,iv)
     $                              * alpha (:,:,jv)*keff(jmode)/bigr,1)

            int(2,3,:,iv,jv)=int(3,2,:,jv,iv)

            int(3,3,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                          ziso*(dalpdz(:,:,iv)* dalpdz(:,:,jv)
     $                              +dalpdrc(:,:,iv)*dalpdrc(:,:,jv))
     $              + dt*fdivb*divbd*alpha(:,:,iv)
     $                              *alpha(:,:,jv)*k2ef(jmode)/bigr2,1)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     same but no div(b) cleaning.
c-----------------------------------------------------------------------
      ELSE
        DO iv=1,nv
          DO jv=1,nv
            int(1,1,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                        ziso*(dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $              +  alpha(:,:,iv)*alpha(:,:,jv)*k2ef(jmode)/bigr2),1)
            int(2,1,:,jv,iv)=SUM(
     $                       -ziso* dalpdz(:,:,iv)*dalpdr(:,:,jv),1)
            int(1,2,:,iv,jv)=int(2,1,:,jv,iv)
            int(3,1,:,jv,iv)=SUM(
     $                        ziso*  alpha(:,:,iv)*keff(jmode)/bigr
     $                            *dalpdrc(:,:,jv),1)
            int(1,3,:,iv,jv)=int(3,1,:,jv,iv)

            int(2,2,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                        ziso*(dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $              +  alpha(:,:,iv)*alpha(:,:,jv)*k2ef(jmode)/bigr2),1)
            int(3,2,:,jv,iv)=SUM(
     $                        ziso*  alpha(:,:,iv)*keff(jmode)/bigr
     $                             *dalpdz(:,:,jv),1)
            int(2,3,:,iv,jv)=int(3,2,:,jv,iv)

            int(3,3,:,jv,iv)=SUM( alpha(:,:,iv)*alpha(:,:,jv)+
     $                        ziso*(dalpdz(:,:,iv)* dalpdz(:,:,jv)
     $                            +dalpdrc(:,:,iv)*dalpdrc(:,:,jv)),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curl_de_iso
c-----------------------------------------------------------------------
c     subprogram 3. curl_de_ciso.
c     compute the integrand for the symmetric part of the resistive
c     diffusion operator for magnetic field.  the matrix is used
c     for preconditioning matrix-free iterations for 3d resistivity.
c-----------------------------------------------------------------------
      SUBROUTINE curl_de_ciso(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: elecd_n0,dp
      REAL(r8), DIMENSION(1,1,1) :: dv
      REAL(r8) :: ziso,bigr2,massfac
      INTEGER(i4) :: iv,jv,nv,ix,iy,ncx,ncy
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and the symmetric
c     part of the electrical diffusivity.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                     inode,elecd_n0,dp,dp,0_i4)
      nv=SIZE(int,5)
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
c-----------------------------------------------------------------------
c     compute the curl(alpha_i)*Z*curl(alpha_j) operator.
c     Poloidal components are either Cartesian or cylindrical.
c     note:  jmode is set in the calling routine.
c
c     multiply electrical diffusivity by the centering parameter and dt.
c
c     if div(b) cleaning is not time split, add the grad(div) operator.
c-----------------------------------------------------------------------
      IF (divbd>0.AND..NOT.split_divb) THEN
       DO iv=1,nv
        DO jv=1,nv
         DO iy=1,ncy
          int(:  ,1,iy,jv,iv)=0._r8
          int(2:3,2,iy,jv,iv)=0._r8
          int(3  ,3,iy,jv,iv)=0._r8
          DO ix=1,ncx
            massfac=alpha(ix,iy,iv)*alpha(ix,iy,jv)
            ziso=dt*feta*elecd_n0(1,ix,iy)
            bigr2=bigr(ix,iy)*bigr(ix,iy)
            int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+massfac+
     $                  ziso*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $        +massfac*k2ef(jmode)/bigr2)
     $        +dt*fdivb*divbd*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
            int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-
     $                  ziso*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)
     $        +dt*fdivb*divbd*dalpdrc(ix,iy,iv)*dalpdz (ix,iy,jv)
            int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)+
     $                  (0,1)*keff(jmode)/bigr(ix,iy)*
     $        (          ziso*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv)
     $        +dt*fdivb*divbd*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv))

            int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+massfac+
     $                  ziso*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $          +massfac*k2ef(jmode)/bigr2)
     $          +dt*fdivb*divbd*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
            int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv)+
     $                  (0,1)*keff(jmode)/bigr(ix,iy)*
     $        (          ziso*dalpdz(ix,iy,jv)*alpha(ix,iy,iv)
     $        +dt*fdivb*divbd*dalpdz(ix,iy,iv)*alpha(ix,iy,jv))

            int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)+massfac+ziso*
     $        (dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
     $        +dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv))
     $        +dt*fdivb*divbd*massfac*
     $            k2ef(jmode)/bigr2
          ENDDO
          int(1,2,iy,iv,jv)=CONJG(int(2,1,iy,jv,iv))
          int(1,3,iy,iv,jv)=CONJG(int(3,1,iy,jv,iv))
          int(2,3,iy,iv,jv)=CONJG(int(3,2,iy,jv,iv))
         ENDDO
        ENDDO
       ENDDO
c-----------------------------------------------------------------------
c     same without div(b) cleaning.
c-----------------------------------------------------------------------
      ELSE
       DO iv=1,nv
        DO jv=1,nv
         DO iy=1,ncy
          int(:  ,1,iy,jv,iv)=0._r8
          int(2:3,2,iy,jv,iv)=0._r8
          int(3  ,3,iy,jv,iv)=0._r8
          DO ix=1,ncx
            massfac=alpha(ix,iy,iv)*alpha(ix,iy,jv)
            ziso=dt*feta*elecd_n0(1,ix,iy)
            bigr2=bigr(ix,iy)*bigr(ix,iy)
            int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+massfac+
     $                   ziso*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                  +massfac*k2ef(jmode)/bigr2)
            int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-
     $                   ziso*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)
            int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)+
     $                   (0,1)*keff(jmode)/bigr(ix,iy)*
     $                   ziso*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv)

            int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+massfac+
     $                   ziso*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                  +massfac*k2ef(jmode)/bigr2)
            int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv)+
     $                   (0,1)*keff(jmode)/bigr(ix,iy)*
     $                   ziso*dalpdz(ix,iy,jv)*alpha(ix,iy,iv)

            int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)+massfac+ziso*
     $        (dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
     $        +dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv))
          ENDDO
          int(1,2,iy,iv,jv)=CONJG(int(2,1,iy,jv,iv))
          int(1,3,iy,iv,jv)=CONJG(int(3,1,iy,jv,iv))
          int(2,3,iy,iv,jv)=CONJG(int(3,2,iy,jv,iv))
         ENDDO
        ENDDO
       ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE curl_de_ciso
c----------------------------------------------------------------------
c     subprogram 4. b_hmhd_op.
c
c     find the diagonal-in-Fourier-index part of the implicit 
c     operator used for the combined Hall-MHD advance of B.  this is
c     used as the implicit operator for linear HMHD computations and for
c     the preconditioner matrix with nonlinear HMHD.  In addition, part
c     of it can be used to reduce the number of matrix-free computations
c     in nonlinear advances.
c
c     the constructed operator is curl(A*).de-div(A*)div(db), where A 
c     is a test vector, de is
c
c       de=[dj X (Beq+b0) + (Jeq+j0) X db]/ne - (Veq+v0) X db + eta*dj
c
c     and the div(db) terms are for divergence cleaning.
c
c     this routine is also called for implicit advection for MHD; two-
c     fluid terms are switched-off according to integrand_flag.
c-----------------------------------------------------------------------
      SUBROUTINE b_hmhd_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
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
     $          b0,b0r,b0z,v0,nd_eq,n0,beq,
     $          jeq,veq,elecd_n0,ds2,nl_pres,dp
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be,ja
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: temp1,ziso

      COMPLEX(r8), DIMENSION(3,3,SIZE(bigr,1),SIZE(bigr,2)) :: curlxb0

      REAL(r8) :: fdb,hdt,hfac,bigr2,ziss,krow,massfac,disdbfac,fhyp,
     $            fhdb
      INTEGER(i4) :: iv,jv,nv,ix,iy,ncx,ncy,ioff,nvc,nvd,ivd,jvd
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions
c
c     if poly_divb is non-negative, divergence error is controlled with
c     an auxiliary field that is discontinuous at element borders, and
c     these bases are also needed.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      IF (poly_divb>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlmat',aldis,
     $                          dalddr,dalddz,0_i4,poly_divb,
     $                          polydmin=poly_divb_min,
     $                          polydmax=poly_divb_max)
      ENDIF
      nvc=ncontb
      nvd=ndiscb
      nv=nvc+nvd
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
c-----------------------------------------------------------------------
c     dissipation coefficients.
c-----------------------------------------------------------------------
      IF (elecd>0) THEN
        IF (eta_model/='fixed') THEN
          CALL generic_ptr_set(rb%qelecd_n0,tb%qelecd_n0,tb%tgeom,
     $                         inode,elecd_n0,dp,dp,0_i4)
          ziso=dt*feta*elecd_n0(1,:,:)
        ELSE IF (ds_use=='elecd'.OR.ds_use=='both') THEN
          CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
          ziso=dt*feta*elecd*MAX(ds2(1,:,:),0._r8)
        ELSE
          ziso=dt*feta*elecd
        ENDIF
      ELSE
        ziso=0._r8
      ENDIF
      IF (divbd>0.AND..NOT.split_divb) THEN
        fdb=dt*fdivb*divbd
      ELSE
        fdb=0._r8
      ENDIF
      disdbfac=SQRT(disc_dbd*dt*fdivb)  !  diffusive correction
c-----------------------------------------------------------------------
c     coefficients for the Hall term and the partial(J)/partial(t)
c     electron inertia term.
c-----------------------------------------------------------------------
      hdt=0.5_r8*dt
      IF (integrand_flag=='mhd') THEN
        hfac=0._r8
      ELSE
        hfac=hdt*coefhll
      ENDIF
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      IF (nonlinear.AND.
     $    (continuity=='full'.OR.continuity=='n=0 only')) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,inode,
     $                       n0,dp,dp,0_i4)
        temp1=hfac/(nd_eq(1,:,:)+n0(1,:,:))
        IF (ohms=='2fl') 
     $    ziso=ziso+si_fac_hall*coefme1/(mu0*(nd_eq(1,:,:)+n0(1,:,:)))
      ELSE
        temp1=hfac/nd_eq(1,:,:)
        IF (ohms=='2fl'.AND.nonlinear) THEN
          ziso=ziso+si_fac_hall*coefme1/(mu0*nd_eq(1,:,:))
        ELSE IF (ohms=='2fl') THEN
          ziso=ziso+coefme1/(mu0*nd_eq(1,:,:))
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     evaluate fields needed for the implicit operator.
c     note that be->0.5*dt*B/mu0*ne and ja->0.5*dt*(J/ne - V)
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,beq,dp,dp,0_i4)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,jeq,dp,dp,0_i4)
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,inode,
     $                       b0,b0r,b0z,1_i4)
        CALL generic_ptr_set(rb%qve_n0,tb%qve_n0,tb%tgeom,inode,
     $                       v0,dp,dp,0_i4)
        be(1,:,:)=temp1*(b0(1,:,:)+beq(1,:,:))/mu0
        be(2,:,:)=temp1*(b0(2,:,:)+beq(2,:,:))/mu0
        be(3,:,:)=temp1*(b0(3,:,:)+beq(3,:,:))/mu0
        ja(1,:,:)=temp1*(jeq(1,:,:)+b0z(3,:,:)/mu0)-hdt*v0(1,:,:)
        IF (geom=='tor') THEN
          ja(2,:,:)=temp1*(jeq(2,:,:)-(b0r(3,:,:)+b0(3,:,:)/bigr)/mu0)-
     $                    hdt*v0(2,:,:)
        ELSE
          ja(2,:,:)=temp1*(jeq(2,:,:)-b0r(3,:,:)/mu0)-hdt*v0(2,:,:)
        ENDIF
        ja(3,:,:)=temp1*(jeq(3,:,:)+(b0r(2,:,:)-b0z(1,:,:))/mu0)-
     $                  hdt*v0(3,:,:)
      ELSE
        be(1,:,:)=temp1*beq(1,:,:)/mu0
        be(2,:,:)=temp1*beq(2,:,:)/mu0
        be(3,:,:)=temp1*beq(3,:,:)/mu0
        ja(1,:,:)=temp1*jeq(1,:,:)
        ja(2,:,:)=temp1*jeq(2,:,:)
        ja(3,:,:)=temp1*jeq(3,:,:)
      ENDIF 
      IF (eq_flow/='none') THEN
        CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                       inode,veq,dp,dp,0_i4)
        ja=ja-hdt*veq
      ENDIF
c-----------------------------------------------------------------------
c     add a numerical term to the isotropic impedance that is
c     proportional to dt times a cyclotron frequency that is based
c     on the variation of B in the periodic coordinate.
c-----------------------------------------------------------------------
c     IF (nonlinear) THEN
c       CALL generic_ptr_set(rb%qsi_nl_pres,tb%qsi_nl_pres,tb%tgeom,
c    $                       inode,nl_pres,dp,dp,0_i4)
c       IF (continuity=='full'.OR.continuity=='n=0 only') THEN
c         ziso=ziso+si_fac_hall*dt*0.5_r8*nl_pres(1,:,:)/
c    $              (SQRT(SUM((b0+beq)**2,1))
c    $               *elementary_q*mu0*(nd_eq(1,:,:)+n0(1,:,:)))
c       ELSE
c         ziso=ziso+si_fac_hall*dt*0.5_r8*nl_pres(1,:,:)/
c    $              (SQRT(SUM((b0+beq)**2,1))
c    $               *elementary_q*mu0*nd_eq(1,:,:))
c       ENDIF
c     ENDIF
c-----------------------------------------------------------------------
c     to help construct the hall implicit operator, we create
c     a matrix that holds the curl of the three vector test functions
c     and a matrix that holds curl(Ar,Az,Aphi) X be + ja X (Ar,Az,Aphi)
c     the first two indices of the temporary matrices are defined as
c     ( result vector index, operand vector index ), where the vector
c     indices runs from 1 to 3 as (r,z,phi), and the last index is
c     the basis function index.  curl(Ar,Az,Aphi) goes first.
c
c     in this version, curl is not formed, and curlxb0 is within the
c     outer jv loop.
c-----------------------------------------------------------------------
      krow=keff(jmode)
      DO jv=1,nvc
        curlxb0(1,1,:,:)= be(3,:,:)*(0,1)*krow*alpha(:,:,jv)/bigr+
     $                    be(2,:,:)*dalpdz(:,:,jv)
        curlxb0(2,1,:,:)=-be(1,:,:)*dalpdz(:,:,jv)+
     $                    ja(3,:,:)*alpha(:,:,jv)
        curlxb0(3,1,:,:)=-be(1,:,:)*(0,1)*krow*alpha(:,:,jv)/bigr-
     $                    ja(2,:,:)*alpha(:,:,jv)
        curlxb0(1,2,:,:)=-be(2,:,:)*dalpdr(:,:,jv)-
     $                    ja(3,:,:)*alpha(:,:,jv)
        curlxb0(2,2,:,:)= be(1,:,:)*dalpdr(:,:,jv)+
     $                    be(3,:,:)*(0,1)*krow*alpha(:,:,jv)/bigr
        curlxb0(3,2,:,:)=-be(2,:,:)*(0,1)*krow*alpha(:,:,jv)/bigr+
     $                    ja(1,:,:)*alpha(:,:,jv)
        curlxb0(1,3,:,:)=-be(3,:,:)*dalpdrc(:,:,jv)+
     $                    ja(2,:,:)*alpha(:,:,jv)
        curlxb0(2,3,:,:)=-be(3,:,:)*dalpdz(:,:,jv)-
     $                    ja(1,:,:)*alpha(:,:,jv)
        curlxb0(3,3,:,:)= be(2,:,:)*dalpdz(:,:,jv)+
     $                    be(1,:,:)*dalpdrc(:,:,jv)
c-----------------------------------------------------------------------
c     compute the implicit hmhd operator.
c
c       CC(curl(A)).( (curl(dB))xB0/mu*ne + (J0/ne-V0)xdB + eta*dJ )
c     
c     where A is the vector test function, whose indices correspond
c     to row indices. the diagonal term is added elsewhere.
c-----------------------------------------------------------------------
       DO iv=1,nvc
        DO iy=1,ncy
         int(:,:,iy,jv,iv)=0._r8
         DO ix=1,ncx
          massfac=alpha(ix,iy,iv)*alpha(ix,iy,jv)
          bigr2=bigr(ix,iy)*bigr(ix,iy)
          ziss=ziso(ix,iy)

          int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+massfac-
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(2,1,ix,iy)-
     $        dalpdz(ix,iy,iv)*curlxb0(3,1,ix,iy)+
     $        fdb*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)+
     $        ziss*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $              +massfac*k2ef(jmode)/bigr2)

          int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(2,2,ix,iy)-
     $        dalpdz(ix,iy,iv)*curlxb0(3,2,ix,iy)+
     $        fdb*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)-
     $        ziss*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)             

          int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)-
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(2,3,ix,iy)-
     $        dalpdz(ix,iy,iv)*curlxb0(3,3,ix,iy)+
     $        (0,1)*krow/bigr(ix,iy)*
     $        (fdb *dalpdrc(ix,iy,iv)*alpha(ix,iy,jv)+
     $         ziss*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv))

          int(1,2,iy,jv,iv)=int(1,2,iy,jv,iv)+
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(1,1,ix,iy)+
     $        dalpdr(ix,iy,iv)*curlxb0(3,1,ix,iy)+
     $        fdb*dalpdrc(ix,iy,jv)*dalpdz (ix,iy,iv)-
     $        ziss*dalpdz(ix,iy,jv)*dalpdr(ix,iy,iv)             

          int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+massfac+
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(1,2,ix,iy)+
     $        dalpdr(ix,iy,iv)*curlxb0(3,2,ix,iy)+
     $        fdb*  dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $        ziss*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $              +massfac*k2ef(jmode)/bigr2)

          int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv)+
     $        (0,1)*krow*alpha(ix,iy,iv)/bigr(ix,iy)
     $                        *curlxb0(1,3,ix,iy)+
     $        dalpdr(ix,iy,iv)*curlxb0(3,3,ix,iy)+
     $        (0,1)*krow/bigr(ix,iy)*
     $        (fdb* dalpdz(ix,iy,iv)*alpha(ix,iy,jv)+
     $         ziss*dalpdz(ix,iy,jv)*alpha(ix,iy,iv))

          int(1,3,iy,jv,iv)=int(1,3,iy,jv,iv)+
     $        dalpdz (ix,iy,iv)*curlxb0(1,1,ix,iy)-
     $        dalpdrc(ix,iy,iv)*curlxb0(2,1,ix,iy)-
     $        (0,1)*krow/bigr(ix,iy)*
     $        (fdb *dalpdrc(ix,iy,jv)*alpha(ix,iy,iv)+
     $         ziss*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv))

          int(2,3,iy,jv,iv)=int(2,3,iy,jv,iv)+
     $        dalpdz (ix,iy,iv)*curlxb0(1,2,ix,iy)-
     $        dalpdrc(ix,iy,iv)*curlxb0(2,2,ix,iy)-
     $        (0,1)*krow/bigr(ix,iy)*
     $        (fdb* dalpdz(ix,iy,jv)*alpha(ix,iy,iv)+
     $         ziss*dalpdz(ix,iy,iv)*alpha(ix,iy,jv))

          int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)+massfac+
     $        dalpdz (ix,iy,iv)*curlxb0(1,3,ix,iy)-
     $        dalpdrc(ix,iy,iv)*curlxb0(2,3,ix,iy)+
     $        fdb*massfac*k2ef(jmode)/bigr2+
     $        ziss*(dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
     $             +dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv))
         ENDDO
        ENDDO
       ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     add contributions from the numerical hyper-resistivity if it
c     is used without a split step.  in the rows for dB, the following
c     adds
c
c       -fhyp*curl(A^*).curl(D)
c       -fhdb* div(A^*)* div(D)
c
c     and the rows for the auxiliary vector D are
c
c       (C^*).D + fhyp*curl(C^*).curl(dB)
c               + fhdb* div(C^*)* div(dB)
c
c     where A is the test vector for dB and C is the test vector for D.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0.OR.hyp_dbd>0).AND..NOT.split_hypeta) THEN
        fhyp=SQRT(hyp_eta*fhyp_eta*dt)
        fhdb=SQRT(hyp_dbd*fhyp_dbd*dt)
        DO iv=1,nvc
         DO jv=1,nvc
          DO iy=1,ncy
           DO ix=1,ncx
             massfac=alpha(ix,iy,iv)*alpha(ix,iy,jv)
             bigr2=bigr(ix,iy)*bigr(ix,iy)

             int(4,1,iy,jv,iv)=int(4,1,iy,jv,iv)-
     $           fhyp*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)-
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(5,1,iy,jv,iv)=int(5,1,iy,jv,iv)+
     $           fhyp*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)-
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(6,1,iy,jv,iv)=int(6,1,iy,jv,iv)-
     $           (fhyp*alpha(ix,iy,iv)*dalpdrc(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdrc(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)

             int(4,2,iy,jv,iv)=int(4,2,iy,jv,iv)+
     $           fhyp*dalpdr(ix,iy,iv)*dalpdz(ix,iy,jv)-
     $           fhdb*dalpdz(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(5,2,iy,jv,iv)=int(5,2,iy,jv,iv)-
     $           fhyp*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)-
     $           fhdb*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(6,2,iy,jv,iv)=int(6,2,iy,jv,iv)-
     $           (fhyp*alpha(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdz(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)

             int(4,3,iy,jv,iv)=int(4,3,iy,jv,iv)+
     $           (fhyp*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(5,3,iy,jv,iv)=int(5,3,iy,jv,iv)+
     $           (fhyp*dalpdz(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdz(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(6,3,iy,jv,iv)=int(6,3,iy,jv,iv)-
     $           fhyp*( dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv)
     $                 +dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv))-
     $           fhdb*massfac*k2ef(jmode)/bigr2

             int(1,4,iy,jv,iv)=int(1,4,iy,jv,iv)+
     $           fhyp*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)+
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(2,4,iy,jv,iv)=int(2,4,iy,jv,iv)-
     $           fhyp*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)+
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(3,4,iy,jv,iv)=int(3,4,iy,jv,iv)+
     $           (fhyp*alpha(ix,iy,iv)*dalpdrc(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdrc(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(4,4,iy,jv,iv)=int(4,4,iy,jv,iv)+massfac

             int(1,5,iy,jv,iv)=int(1,5,iy,jv,iv)-
     $           fhyp*dalpdr(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $           fhdb*dalpdz(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(2,5,iy,jv,iv)=int(2,5,iy,jv,iv)+
     $           fhyp*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)+
     $           fhdb*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(3,5,iy,jv,iv)=int(3,5,iy,jv,iv)+
     $           (fhyp*alpha(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdz(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(5,5,iy,jv,iv)=int(5,5,iy,jv,iv)+massfac

             int(1,6,iy,jv,iv)=int(1,6,iy,jv,iv)-
     $           (fhyp*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(2,6,iy,jv,iv)=int(2,6,iy,jv,iv)-
     $           (fhyp*dalpdz(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdz(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(3,6,iy,jv,iv)=int(3,6,iy,jv,iv)+
     $           fhyp*( dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv)
     $                 +dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv))+
     $           fhdb*massfac*k2ef(jmode)/bigr2
             int(6,6,iy,jv,iv)=int(6,6,iy,jv,iv)+massfac
           ENDDO
          ENDDO
         ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find the contributions associated with the rows and columns of the
c     discontinuous scalar for divb control.
c     
c     the equations are expressed directly in weak form:
c
c       int[A^*.delta(B)] = -dt*int[curl(A^*).E]
c               - int[ sqrt(disc_dbd*dt*f)*auxb*div(A^*) ]
c       int[(phi^*)*auxb]-int[sqrt(disc_dbd*dt*f)*(phi^*)*div(delta(B))]
c         = int[ sqrt(disc_dbd*nu*dt/f)*(phi^*)*div(B_old)]
c
c     where f is the implicit centering and phi is the scalar test
c     field for the auxiliary field auxb.
c-----------------------------------------------------------------------
      IF (poly_divb>=0) THEN
        DO iv=1,nvc
          DO jv=nvc+1,nv
            jvd=jv-nvc
            DO iy=1,ncy
              int(:,:,iy,jv,iv)=0._r8
              DO ix=1,ncx
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+
     $            disdbfac*dalpdrc(ix,iy,iv)*aldis(ix,iy,jvd)
                int(1,2,iy,jv,iv)=int(1,2,iy,jv,iv)+
     $            disdbfac*dalpdz(ix,iy,iv)*aldis(ix,iy,jvd)
                int(1,3,iy,jv,iv)=int(1,3,iy,jv,iv)-(0,1)*krow*
     $            disdbfac*alpha(ix,iy,iv)*aldis(ix,iy,jvd)/bigr(ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO iv=nvc+1,nv
          ivd=iv-nvc
          DO jv=1,nvc
            DO iy=1,ncy
              int(:,1,iy,jv,iv)=0._r8
              DO ix=1,ncx
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)-
     $            disdbfac*dalpdrc(ix,iy,jv)*aldis(ix,iy,ivd)
                int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-
     $            disdbfac*dalpdz(ix,iy,jv)*aldis(ix,iy,ivd)
                int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)-(0,1)*krow*
     $            disdbfac*alpha(ix,iy,jv)*aldis(ix,iy,ivd)/bigr(ix,iy)
              ENDDO
            ENDDO
          ENDDO
          DO jv=nvc+1,nv
            jvd=jv-nvc
            DO iy=1,ncy
              int(1,1,iy,jv,iv)=0._r8
              DO ix=1,ncx
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+
     $            aldis(ix,iy,jvd)*aldis(ix,iy,ivd)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE b_hmhd_op
c----------------------------------------------------------------------
c     subprogram 5. b_hyp_op.
c
c     find the implicit operator used for the time-split hyper-
c     diffusion advance of B, where the rows for delta-B are:
c
c       (A^*).dB - fhyp*curl(A^*).curl(D)
c                - fhdb* div(A^*)* div(D)
c
c     and the rows for the auxiliary vector D are
c
c       (C^*).D + fhyp*curl(C^*).curl(dB)
c               + fhdb* div(C^*)* div(dB)
c
c     where A is the test vector for dB and C is the test vector for D.
c-----------------------------------------------------------------------
      SUBROUTINE b_hyp_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8) :: fhyp,fhdb,krow,bigr2,massfac
      INTEGER(i4) :: iv,jv,ix,iy,ncx,ncy,nvc
c-----------------------------------------------------------------------
c     set index limits and real constants.
c-----------------------------------------------------------------------
      nvc=ncontb
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      krow=keff(jmode)
      fhyp=SQRT(hyp_eta*fhyp_eta*dt)
      fhdb=SQRT(hyp_dbd*fhyp_dbd*dt)
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
c-----------------------------------------------------------------------
c     loop over basis functions and add contributions for the
c     mass matrices and for the hyper-resistivity operator.
c-----------------------------------------------------------------------
      DO iv=1,nvc
        DO jv=1,nvc
          DO iy=1,ncy
           int(:,:,iy,jv,iv)=0._r8
           DO ix=1,ncx
             massfac=alpha(ix,iy,iv)*alpha(ix,iy,jv)
             bigr2=bigr(ix,iy)*bigr(ix,iy)

             int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+massfac
             int(4,1,iy,jv,iv)=int(4,1,iy,jv,iv)-
     $           fhyp*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)-
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(5,1,iy,jv,iv)=int(5,1,iy,jv,iv)+
     $           fhyp*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)-
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(6,1,iy,jv,iv)=int(6,1,iy,jv,iv)-
     $           (fhyp*alpha(ix,iy,iv)*dalpdrc(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdrc(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)

             int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+massfac
             int(4,2,iy,jv,iv)=int(4,2,iy,jv,iv)+
     $           fhyp*dalpdr(ix,iy,iv)*dalpdz(ix,iy,jv)-
     $           fhdb*dalpdz(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(5,2,iy,jv,iv)=int(5,2,iy,jv,iv)-
     $           fhyp*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)-
     $           fhdb*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(6,2,iy,jv,iv)=int(6,2,iy,jv,iv)-
     $           (fhyp*alpha(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdz(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)

             int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)+massfac
             int(4,3,iy,jv,iv)=int(4,3,iy,jv,iv)+
     $           (fhyp*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(5,3,iy,jv,iv)=int(5,3,iy,jv,iv)+
     $           (fhyp*dalpdz(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdz(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(6,3,iy,jv,iv)=int(6,3,iy,jv,iv)-
     $           fhyp*( dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv)
     $                 +dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv))-
     $           fhdb*massfac*k2ef(jmode)/bigr2

             int(1,4,iy,jv,iv)=int(1,4,iy,jv,iv)+
     $           fhyp*(dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)+
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(2,4,iy,jv,iv)=int(2,4,iy,jv,iv)-
     $           fhyp*dalpdz(ix,iy,iv)*dalpdr(ix,iy,jv)+
     $           fhdb*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(3,4,iy,jv,iv)=int(3,4,iy,jv,iv)+
     $           (fhyp*alpha(ix,iy,iv)*dalpdrc(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdrc(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(4,4,iy,jv,iv)=int(4,4,iy,jv,iv)+massfac

             int(1,5,iy,jv,iv)=int(1,5,iy,jv,iv)-
     $           fhyp*dalpdr(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $           fhdb*dalpdz(ix,iy,iv)*dalpdrc(ix,iy,jv)
             int(2,5,iy,jv,iv)=int(2,5,iy,jv,iv)+
     $           fhyp*(dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                 +massfac*k2ef(jmode)/bigr2)+
     $           fhdb*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
             int(3,5,iy,jv,iv)=int(3,5,iy,jv,iv)+
     $           (fhyp*alpha(ix,iy,iv)*dalpdz(ix,iy,jv)+
     $            fhdb*alpha(ix,iy,jv)*dalpdz(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(5,5,iy,jv,iv)=int(5,5,iy,jv,iv)+massfac

             int(1,6,iy,jv,iv)=int(1,6,iy,jv,iv)-
     $           (fhyp*dalpdrc(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdrc(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(2,6,iy,jv,iv)=int(2,6,iy,jv,iv)-
     $           (fhyp*dalpdz(ix,iy,iv)*alpha(ix,iy,jv)+
     $            fhdb*dalpdz(ix,iy,jv)*alpha(ix,iy,iv))*
     $           (0,1)*krow/bigr(ix,iy)
             int(3,6,iy,jv,iv)=int(3,6,iy,jv,iv)+
     $           fhyp*( dalpdz (ix,iy,iv)*dalpdz (ix,iy,jv)
     $                 +dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv))+
     $           fhdb*massfac*k2ef(jmode)/bigr2
             int(6,6,iy,jv,iv)=int(6,6,iy,jv,iv)+massfac
           ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE b_hyp_op
c-----------------------------------------------------------------------
c     subprogram 6. n_iso_op.
c     compute an isotropic operator for the number density advance.
c-----------------------------------------------------------------------
      SUBROUTINE n_iso_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: nd_eq,dart,dp
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2,kprp
      REAL(r8), DIMENSION(1,1,1) :: dv
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: nd,dcp
      INTEGER(i4) :: nv,iv,jv,iq,jq
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      nv=SIZE(int,5)
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     electrons and ions have different coefficients.  this routine
c     may also be used for a diffusive operator for the continuity eqn.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:2)=='n ') THEN
        IF (nonlinear.AND.nd_floor>0.AND.nd_diff>0) THEN
          CALL generic_ptr_set(rb%qdart,tb%qdart,tb%tgeom,
     $                         inode,dart,dp,dp,0_i4)
          kprp=fthc*dt*nd_diff+dart(1,:,:)
        ELSE
          kprp=fthc*dt*nd_diff
        ENDIF
      ELSE IF (integrand_flag(1:6)=='all ti') THEN
        kprp=fthc*dt*k_perpi
      ELSE
        kprp=fthc*dt*k_perpe
      ENDIF
c-----------------------------------------------------------------------
c     Compute grad(alpha(i))*grad(alpha(j)).  it is assumed that
c     k_perp has a factor of gamma-1, so that k_perp is the diffusivity
c     for temperature.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
           int(1,1,:,jv,iv)=SUM((dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $                          +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $                          +k2ef(jmode)*alpha(:,:,iv)*alpha(:,:,jv)
     $                          /bigr2)*kprp
     $                          +alpha(:,:,iv)*alpha(:,:,jv),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE n_iso_op
c-----------------------------------------------------------------------
c     subprogram 7. t_aniso_op.
c     compute the symmetric part of the operator for temperature advance
c
c     (n*I.-dt*f*D)
c
c     including an implicit isotropic diffusion operator and the
c     diagonal part of an anisotropic temperature diffusivity.
c     when implicit advection is used, the operator is
c
c     (n*I. + 0.5*dt*n*( V.grad(I.)+(gamma-1)*div(V)*I. ) - dt*f*D)
c
c-----------------------------------------------------------------------
      SUBROUTINE t_aniso_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          bb,nd_eq,nd_n0,kappl,ndiff,
     $          kaprp,be_eq,t_b2,ja_eq,ve_eq,divvq,v0,v0r,v0z,
     $          b0,b0r,b0z,nqr,nqz,nd0r,nd0z,vv,upave,bcrgt,t_eq,dp
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2,n0,
     $          prp_fac,pll_fac,ndivvs,fupw,fupwi
      REAL(r8), DIMENSION(6,SIZE(bigr,1),SIZE(bigr,2)) :: updyad
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: nvs,ja
      REAL(r8), DIMENSION(1,1,1) :: dv
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: ctmp
      COMPLEX(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: csclr
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,jv,iq,jq
      REAL(r8) :: zz,kpll,kprp,jfac,g,gm1
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions and the bb dyad.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      IF (p_model(1:5)=='aniso')
     $  CALL generic_ptr_set(rb%qbb,tb%qbb,tb%tgeom,inode,bb,dp,dp,0_i4)
      nv=SIZE(int,5)
      IF (integrand_flag(1:6)=='all ti') THEN
        zz=zeff
        IF (separate_pe) THEN
          kpll=k_plli
          jfac=0.5_r8*dt*coefjvi/zeff
        ELSE
          kpll=k_plle
          jfac=0._r8
        ENDIF
        kprp=k_perpi
      ELSE
        zz=1._r8
        kpll=k_plle
        kprp=k_perpe
        jfac=0.5_r8*dt*coefjve
      ENDIF
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
      IF (p_model=='isothermal') THEN
        gm1=0._r8
      ELSE
        gm1=gamm1
      ENDIF
      IF (p_model=='isotropic') kpll=kprp
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     number density, n, is needed for anisotropic cases, so the
c     equation is n*dT/dt= rather than dT/dt=.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,nqr,nqz,1_i4)
      IF (nonlinear.AND.
     $    (continuity=='n=0 only'.OR.continuity=='full') ) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,nd_n0,nd0r,nd0z,1_i4)
        n0=(nd_eq(1,:,:)+nd_n0(1,:,:))/zz
      ELSE
        n0=nd_eq(1,:,:)/zz
        NULLIFY(nd0r)
      ENDIF
c-----------------------------------------------------------------------
c     Compute the isotropic diffusion coefficient.  k_perp and k_pll
c     already include the factor of gamma-1.  multiply by the
c     appropriate z, too.  [this will not give the exact n=0 part of
c     the operator when continuity=full, but it should be sufficient
c     for the preconditioner in most cases.]
c
c-PRE also use dyad?
c-----------------------------------------------------------------------
      IF (nonlinear.AND.(p_model=="aniso_plltdep".OR.
     $                   p_model=="aniso_tdep")) THEN
        IF (integrand_flag(1:6)=='all ti') THEN
          CALL generic_ptr_set(rb%qkappli_n0,tb%qkappli_n0,tb%tgeom,
     $                         inode,kappl,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qkaprpi_n0,tb%qkaprpi_n0,tb%tgeom,
     $                         inode,kaprp,dp,dp,0_i4)
        ELSE
          CALL generic_ptr_set(rb%qkapple_n0,tb%qkapple_n0,tb%tgeom,
     $                         inode,kappl,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qkaprpe_n0,tb%qkaprpe_n0,tb%tgeom,
     $                         inode,kaprp,dp,dp,0_i4)
        ENDIF
        pll_fac=dt*fthc*MAX(0._r8,kappl(1,:,:)-kaprp(1,:,:))*n0
        prp_fac=dt*fthc*kaprp(1,:,:)*n0
      ELSE IF (.NOT.(p_model=='adiabat'.OR.p_model=='isothermal')) THEN
        pll_fac=dt*fthc*(kpll-kprp)*n0
        prp_fac=dt*fthc*kprp*n0
      ENDIF
c-----------------------------------------------------------------------
c     artificial thermal diffusivity for upwinding.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        IF (integrand_flag(1:6)=='all ti') THEN
          CALL generic_ptr_set(rb%qupti_n0,tb%qupti_n0,tb%tgeom,
     $                         inode,upave,dp,dp,0_i4)
        ELSE
          CALL generic_ptr_set(rb%qupte_n0,tb%qupte_n0,tb%tgeom,
     $                         inode,upave,dp,dp,0_i4)
        ENDIF
        CALL generic_ptr_set(rb%qvv,tb%qvv,tb%tgeom,
     $                       inode,vv,dp,dp,0_i4)
        fupw =t_dart_upw*dt**2*n0*upave(1,:,:)*upw_aniso/
     $        ( SUM(vv(1:3,:,:),1) + smallnum )
        fupwi=t_dart_upw*dt**2*n0*upave(1,:,:)*(1._r8-upw_aniso)
        updyad(1,:,:)= fupw*vv(1,:,:)+fupwi
        updyad(2,:,:)= fupw*vv(2,:,:)+fupwi
        updyad(3,:,:)=(fupw*vv(3,:,:)+fupwi)*k2ef(jmode)/bigr2
        updyad(4,:,:)= fupw*vv(4,:,:)
        updyad(5,:,:)= fupw*vv(5,:,:)*keff(jmode)/bigr
        updyad(6,:,:)= fupw*vv(6,:,:)*keff(jmode)/bigr
      ENDIF
c-----------------------------------------------------------------------
c     Compute
c       prp_fac*grad(alpha(i))*grad(alpha(j))+dt*(kpll-kprp)*
c         (grad(alpha(i)) dot B0) (B0 dot grad(alpha(j)) /B^{2}).
c-----------------------------------------------------------------------
      IF (p_model=='adiabat'.OR.p_model=='isothermal') THEN
        DO iv=1,nv
          DO jv=1,iv
            int(1,1,:,jv,iv)=SUM(alpha(:,:,iv)*alpha(:,:,jv)*n0,1)
          ENDDO
        ENDDO
      ELSE IF (p_model=='isotropic') THEN
        DO iv=1,nv
          DO jv=1,iv
            int(1,1,:,jv,iv)=SUM(
     $        (dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $        +dalpdz(:,:,iv)*dalpdz(:,:,jv))*prp_fac
     $        + alpha(:,:,iv)* alpha(:,:,jv)*(n0+
     $                   k2ef(jmode)*prp_fac/bigr2),1)
          ENDDO
        ENDDO
      ELSE
        DO iv=1,nv
          DO jv=1,iv
            int(1,1,:,jv,iv)=SUM(
     $         dalpdr(:,:,iv)*dalpdr(:,:,jv)*(prp_fac+pll_fac*bb(1,:,:))
     $        +dalpdz(:,:,iv)*dalpdz(:,:,jv)*(prp_fac+pll_fac*bb(2,:,:))
     $        + alpha(:,:,iv)* alpha(:,:,jv)*(n0+
     $                   k2ef(jmode)*(prp_fac+pll_fac*bb(3,:,:))/bigr2)
     $        +( (dalpdr(:,:,iv)*dalpdz(:,:,jv)
     $           +dalpdz(:,:,iv)*dalpdr(:,:,jv))*bb(4,:,:)
     $      +(0,1)*(
     $       (dalpdr(:,:,iv)*alpha(:,:,jv)-dalpdr(:,:,jv)*alpha(:,:,iv))
     $                      *keff(jmode)*bb(5,:,:)/bigr+
     $       (dalpdz(:,:,iv)*alpha(:,:,jv)-dalpdz(:,:,jv)*alpha(:,:,iv))
     $                      *keff(jmode)*bb(6,:,:)/bigr ) )*pll_fac,1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     add the anisotropic artificial diffusivity if needed.
c
c-PRE premultiply terms that do not change with jv?
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
        DO iv=1,nv
          DO jv=1,iv
            int(1,1,:,jv,iv)=int(1,1,:,jv,iv)+SUM(
     $          dalpdr(:,:,iv)*dalpdr(:,:,jv)*updyad(1,:,:)
     $         +dalpdz(:,:,iv)*dalpdz(:,:,jv)*updyad(2,:,:)
     $         + alpha(:,:,iv)* alpha(:,:,jv)*updyad(3,:,:)
     $        +(dalpdr(:,:,iv)*dalpdz(:,:,jv)
     $         +dalpdz(:,:,iv)*dalpdr(:,:,jv))*updyad(4,:,:)
     $         +(0,1)*(
     $       (dalpdr(:,:,iv)*alpha(:,:,jv)-dalpdr(:,:,jv)*alpha(:,:,iv))
     $                      *updyad(5,:,:)+
     $       (dalpdz(:,:,iv)*alpha(:,:,jv)-dalpdz(:,:,jv)*alpha(:,:,iv))
     $                      *updyad(6,:,:) ),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     correction factor for artificial particle diffusion.
c-----------------------------------------------------------------------
      IF (continuity/='none'.AND.nonlinear.AND.impladv.AND.
     $    nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        CALL generic_ptr_set(rb%qndiff_n0,tb%qndiff_n0,tb%tgeom,
     $                       inode,ndiff,dp,dp,0_i4)
        DO iv=1,nv
          DO jv=1,iv
            int(1,1,:,jv,iv)=int(1,1,:,jv,iv)+fthc*dt/zz*
     $        SUM(alpha(:,:,iv)*alpha(:,:,jv)*ndiff(1,:,:),1)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     use symmetry to complete the rest of the matrix elements.
c-----------------------------------------------------------------------
      DO iv=1,nv-1
        DO jv=iv+1,nv
          int(1,1,:,jv,iv)=CONJG(int(1,1,:,iv,jv))
        ENDDO
      ENDDO

      IF (.NOT.(separate_pe.AND.k_cross>0.OR.impladv)) RETURN
c-----------------------------------------------------------------------
c     add or subtract the grad(alpha(iv)).B_eqXgrad(alpha(jv)) terms
c     for drift effects.  complete the coefficient in the temporary 
c     prp_fac storage location.
c
c     note that n0->0.5*dt*n0
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,be_eq,dp,dp,0_i4)
      n0=0.5_r8*dt*n0

      IF (separate_pe.AND.k_cross>0.AND.p_model(1:5)=='aniso') THEN
        IF (integrand_flag(1:6)=='all ti') THEN
          CALL generic_ptr_set(rb%qti_b2,tb%qti_b2,tb%tgeom,
     $                         inode,t_b2,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qbcrgti,tb%qbcrgti,tb%tgeom,
     $                         inode,bcrgt,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qtion_eq,tb%qtion_eq,tb%tgeom,
     $                         inode,t_eq,dp,dp,0_i4)
          prp_fac=-t_b2(1,:,:)*n0
        ELSE
          CALL generic_ptr_set(rb%qte_b2,tb%qte_b2,tb%tgeom,
     $                         inode,t_b2,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qbcrgte,tb%qbcrgte,tb%tgeom,
     $                         inode,bcrgt,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qtele_eq,tb%qtele_eq,tb%tgeom,
     $                         inode,t_eq,dp,dp,0_i4)
          prp_fac= t_b2(1,:,:)*n0
        ENDIF
        ctmp(1,:,:)=bcrgt(1,:,:)/t_eq(1,:,:)+
     $              be_eq(2,:,:)*(0,1)*keff(jmode)/bigr
        ctmp(2,:,:)=bcrgt(2,:,:)/t_eq(1,:,:)-
     $              be_eq(1,:,:)*(0,1)*keff(jmode)/bigr
        ctmp(3,:,:)=bcrgt(3,:,:)/t_eq(1,:,:)
      ELSE
        prp_fac=0._r8
        ctmp=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     time-centered implicit advection, n_s*V_s and n_s*div(V_s) for
c     each species s.  the Jeq.grad(neq)=0 approximation is used here.
c
c     for ions, we need
c
c       (dt/2)*n_i*V_i = (dt/2)*(n_e*V/z + J*meomi/(z*e*(1+meomi)) )
c
c     where meomi = z*me/mi.  jfac is the full coefficient for J, and
c     n0 holds dt*n_e/2*z.
c
c     we also need
c
c       (dt/2)*n_i*div(V_i) =
c         (dt/2)*(n_e*div(V)/z - J.grad(n_e)*meomi/(n_e*z*e*(1+meomi))
c
c     the electron terms have z->1 and the meomi in the numerator
c     replaced by -1.
c-----------------------------------------------------------------------
      IF (impladv) THEN
        IF (eq_flow/='none') THEN
          CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                         inode,ve_eq,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qdvv_eq,tb%qdvv_eq,tb%tgeom,
     $                         inode,divvq,dp,dp,0_i4)
        ENDIF

        CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                       inode,ja_eq,dp,dp,0_i4)
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                         inode,b0,b0r,b0z,1_i4)
          IF (.NOT.ASSOCIATED(nd0r))
     $      CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                           inode,nd_n0,nd0r,nd0z,1_i4)
          ja(1,:,:)=  b0z(3,:,:)/mu0
          ja(2,:,:)=-(b0r(3,:,:)+g*b0(3,:,:)/bigr)/mu0
          ja(3,:,:)= (b0r(2,:,:)-b0z(1,:,:))/mu0
          nvs=jfac*(ja_eq+ja)
          jfac=0.5_r8*dt*jfac/zz
          ndivvs=-jfac*(ja_eq(1,:,:)*nd0r(1,:,:)+
     $                  ja_eq(2,:,:)*nd0z(1,:,:)+
     $                  ja(1,:,:)*(nd0r(1,:,:)+nqr(1,:,:))+
     $                  ja(2,:,:)*(nd0z(1,:,:)+nqz(1,:,:)))/n0
        ELSE
          nvs=jfac*ja_eq
          ndivvs=0._r8
        ENDIF

        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qve_n0,tb%qve_n0,tb%tgeom,
     $                         inode,v0,v0r,v0z,1_i4)
          IF (eq_flow/='none') THEN
            nvs(1,:,:)=nvs(1,:,:)+n0*(v0(1,:,:)+ve_eq(1,:,:))
            nvs(2,:,:)=nvs(2,:,:)+n0*(v0(2,:,:)+ve_eq(2,:,:))
            nvs(3,:,:)=nvs(3,:,:)+n0*(v0(3,:,:)+ve_eq(3,:,:))
            ndivvs=ndivvs+
     $        n0*(divvq(1,:,:)+v0r(1,:,:)+v0z(2,:,:)+g*v0(1,:,:)/bigr)
          ELSE
            nvs(1,:,:)=nvs(1,:,:)+n0*v0(1,:,:)
            nvs(2,:,:)=nvs(2,:,:)+n0*v0(2,:,:)
            nvs(3,:,:)=nvs(3,:,:)+n0*v0(3,:,:)
            ndivvs=ndivvs+n0*(v0r(1,:,:)+v0z(2,:,:)+g*v0(1,:,:)/bigr)
          ENDIF
        ELSE
          IF (eq_flow/='none') THEN
            nvs(1,:,:)=nvs(1,:,:)+n0*ve_eq(1,:,:)
            nvs(2,:,:)=nvs(2,:,:)+n0*ve_eq(2,:,:)
            nvs(3,:,:)=nvs(3,:,:)+n0*ve_eq(3,:,:)
            ndivvs=ndivvs+n0*divvq(1,:,:)
          ENDIF
        ENDIF
      ELSE
        nvs=0._r8
        ndivvs=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     add these contributions to the symmetric part of the operator.
c     use ctmp as storage for the coefficient of alpha_jv.
c-----------------------------------------------------------------------
      csclr(:,:)=gm1*ndivvs+(0,1)*keff(jmode)*nvs(3,:,:)/bigr
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=int(1,1,:,jv,iv)+SUM(
     $      prp_fac*(
     $      dalpdr(:,:,iv)*
     $        (ctmp(1,:,:)*alpha(:,:,jv)-be_eq(3,:,:)*dalpdz(:,:,jv))+
     $      dalpdz(:,:,iv)*
     $        (ctmp(2,:,:)*alpha(:,:,jv)+be_eq(3,:,:)*dalpdr(:,:,jv))-
     $      (0,1)*keff(jmode)*alpha(:,:,iv)/bigr*
     $        (be_eq(1,:,:)*dalpdz(:,:,jv)-be_eq(2,:,:)*dalpdr(:,:,jv)+
     $         ctmp(3,:,:)*alpha(:,:,jv)))
     $     +alpha(:,:,iv)*
     $         (nvs(1,:,:)*dalpdr(:,:,jv)+
     $          nvs(2,:,:)*dalpdz(:,:,jv)+
     $          csclr(:,:)* alpha(:,:,jv)),1)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE t_aniso_op
c-----------------------------------------------------------------------
c     subprogram 8. vec_lap_op.
c     find the Laplacian operator used in the velocity advances.  this
c     routine is used in different ways according the the number of
c     vector components at each vertex.  if the local block is 3x3,
c     then the vector order is (r,z,phi); this is used for si-ops in
c     inviscid calculations, where free-slip boundary conditions
c     couple the poloidal components.  if the local block is 2x2, the
c     vector order is (r,phi), where toroidal geometry couples the
c     two components, and z is decoupled with no-slip bcs.  if the
c     local block is 1x1, it's just the z vector component.
c
c     the equation is now rho*delta_V, so there is a mass density
c     scale to all terms, and the rho*mass_matrix term is added here.
c-----------------------------------------------------------------------
      SUBROUTINE vec_lap_op(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          beq,b0,peq,p0,nd_eq,n0,nl_pres,dp,ds2
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq,bigr2,
     $          func,md_sink,diff
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,jv,nq,np
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the
c     diffusivity shape function.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      IF ((ds_use=='kin_visc'.OR.ds_use=='both').AND.kin_visc>0) THEN
        CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                         inode,ds2,dp,dp,0_i4)
        diff=kin_visc*mtot*MAX(ds2(1,:,:),0._r8)
      ELSE
        diff=kin_visc*mtot
      ENDIF
c-----------------------------------------------------------------------
c     the mass density is the md_sink factor.  for nonlinear runs
c     with density evolution, the operator uses only the n=0 part of
c     the number density.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      IF (nonlinear.AND.
     $    (continuity=='n=0 only'.OR.continuity=='full')) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,n0,dp,dp,0_i4)
        md_sink=mtot*(n0(1,:,:)+nd_eq(1,:,:))
        diff=diff*(n0(1,:,:)+nd_eq(1,:,:))
      ELSE
        md_sink=mtot*nd_eq(1,:,:)
        diff=diff*nd_eq(1,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate fields for the semi-implicit operator if needed.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:3)=='mhd'.OR.integrand_flag(1:3)=='all') THEN
        CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                       inode,beq,dp,dp,0_i4)
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,
     $                         inode,b0,dp,dp,0_i4)
          bmagsq=SUM((beq+b0)**2,1)
          CALL generic_ptr_set(rb%qsi_nl_pres,tb%qsi_nl_pres,tb%tgeom,
     $                         inode,nl_pres,dp,dp,0_i4)
        ELSE
          bmagsq=SUM(beq**2,1)
        ENDIF
        IF (beta>0) THEN
          CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                         inode,peq,dp,dp,0_i4)
          IF (nonlinear) THEN
            CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,inode,
     $                           p0,dp,dp,0_i4)
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     find the spatial function that is used with the operator.
c-----------------------------------------------------------------------
      func=0._r8
      IF (kin_visc>0.AND.integrand_flag(1:3)/='mhd')
     $  func=dt*fvsc*diff
      IF (integrand_flag(1:4)/='visc') THEN
        IF (beta>0) THEN
          IF (nonlinear) THEN
            func=func+0.25_r8*dt**2*
     $             (si_fac_mhd*bmagsq/mu0
     $             +si_fac_pres*gamma*(peq(1,:,:)+p0(1,:,:))
     $             +si_fac_nl*(nl_pres(1,:,:)/mu0+gamma*nl_pres(2,:,:)))
          ELSE
            func=func+0.25_r8*dt**2*
     $             (si_fac_mhd*bmagsq/mu0+si_fac_pres*gamma*peq(1,:,:))
          ENDIF
        ELSE
          IF (nonlinear) THEN
            func=func+0.25_r8*dt**2*
     $             (si_fac_mhd*bmagsq/mu0+si_fac_nl*nl_pres(1,:,:)/mu0)
          ELSE
            func=func+0.25_r8*dt**2*si_fac_mhd*bmagsq/mu0
          ENDIF
        ENDIF
      ENDIF
      nv=SIZE(int,5)
      nq=SIZE(int,1)
      IF (nq==3) THEN
        np=3
      ELSE
        np=2
      ENDIF
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     compute func*grad(alpha(i))*grad(alpha(j)) and place this
c     in the 1-1 matrix elements.  [possibly r-r or z-z]
c
c     the mass density term uses the md_sink factor.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM( func
     $          *(dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $           +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $           +k2ef(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2)
     $         +md_sink*alpha(:,:,iv)*alpha(:,:,jv),1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     for an independent z-z matrix, that's it.
c-----------------------------------------------------------------------
      IF (nq==1) RETURN
c-----------------------------------------------------------------------
c     for a full 3x3 matrix, the z-z is the same as r-r prior to
c     toroidal effects.
c-----------------------------------------------------------------------
      IF (nq==3) THEN
        int(1,2,:,:,:)=0._r8
        int(2,2,:,:,:)=int(1,1,:,:,:)
        int(3,2,:,:,:)=0._r8
        int(2,1,:,:,:)=0._r8
        int(2,3,:,:,:)=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     for toroidal geometry, add the coupling between the real_r
c     and -imag_phi (or imag_r and real_phi) components, along with
c     the additional diagonal terms.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        DO iv=1,nv
          DO jv=1,nv
            int(1,1,:,jv,iv)=int(1,1,:,jv,iv)+SUM(
     $        +func*alpha(:,:,iv)*alpha(:,:,jv)/bigr2,1)
            int(np,1,:,jv,iv)=SUM(
     $         func*2*keff(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2,1)
            int(1,np,:,jv,iv)=int(np,1,:,jv,iv)
          ENDDO
        ENDDO
      ELSE
        int(np,1,:,:,:)=0._r8
        int(1,np,:,:,:)=0._r8
      ENDIF
      int(np,np,:,:,:)=int(1,1,:,:,:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vec_lap_op
c-----------------------------------------------------------------------
c     subprogram 9. v_aniso_op.
c     find the semi-implicit operator used in the velocity advances.
c     the integrand is now a complex array.
c
c     the equation is now rho*delta_V, so there is a mass density
c     scale to all terms, and the rho*mass_matrix term is added here.
c-----------------------------------------------------------------------
      SUBROUTINE v_aniso_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
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
     $          b0,b0r,b0z,beq,beqr,beqz,
     $          j0,nd_eq,n0,p0,p0r,p0z,peq,peqr,peqz,nl_pres,bb,dp,ds2,
     $          ti_sym,v0,v0r,v0z,veq,grdv_eq,kappli,ndiff
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: be,ber,bez,ja,
     $          ve,ver,vez
      REAL(r8), DIMENSION(1,SIZE(bigr,1),SIZE(bigr,2)) :: pres,presr,
     $          presz
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bmagsq,bigr2,
     $          iso,grd,md_sink,diff,alpr
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(3,3,SIZE(bigr,1),SIZE(bigr,2),ncontb)
     $             :: crl
      COMPLEX(r8), DIMENSION(9,3,SIZE(bigr,1),SIZE(bigr,2),ncontb)
     $             :: piten
      COMPLEX(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2),ncontb) ::
     $             gr_alpha
      COMPLEX(r8), DIMENSION(3,3,SIZE(bigr,1),SIZE(bigr,2)) :: ptmp
      COMPLEX(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: diva
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      COMPLEX(r8), DIMENSION(3,3) :: wdot,bc_wdot
      COMPLEX(r8), DIMENSION(3) :: ctmp

      INTEGER(i4) :: nv,iv,jv,ix,iy,ncx,ncy,iz,i1,i2,i3,nvc,nvm,ivm,jvm
      REAL(r8) :: g,ani,cj0,b2,gcoef,iso_diag,alprs,dtpv
      REAL(r8), PARAMETER :: third=1._r8/3._r8,twoth=2._r8/3._r8
      COMPLEX(r8) :: dvas,bdga
      LOGICAL :: bbset
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions, and evaluate the
c     diffusivity shape function.
c
c     if poly_divv is non-negative, auxiliary fields are used to
c     stabilize the flow-velocity advance.  the modal bases are
c     needed in this case.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)

      IF (poly_divv>=0) THEN
        CALL generic_alpha_eval(rb,tb%tgeom,inode,'modlmat',almod,dalmr,
     $                          dalmz,0_i4,poly_divv,
     $                          polydmin=poly_divv_min,
     $                          polydmax=poly_divv_max)
      ENDIF

      IF ((ds_use=='kin_visc'.OR.ds_use=='both').AND.
     $         (kin_visc>0.OR.iso_visc>0)) THEN
        CALL generic_ptr_set(rb%qdiff_shape,tb%qdiff_shape,tb%tgeom,
     $                       inode,ds2,dp,dp,0_i4)
        md_sink=mtot
        diff=mtot*MAX(ds2(1,:,:),0._r8)
      ELSE
        md_sink=mtot
        diff=mtot
      ENDIF
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      nvc=ncontb
      nvm=ndiscb
      nv=nvc+nvm
c-----------------------------------------------------------------------
c     the mass density factor is md_sink.  for nonlinear runs
c     with density evolution, the operator uses only the n=0 part of
c     the number density.  for full 3D density evolution, this operator
c     approximates the diagnonal-in-n part of the full operator and is
c     used to make the preconditioner.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qnd_eq,tb%qnd_eq,tb%tgeom,
     $                     inode,nd_eq,dp,dp,0_i4)
      IF (nonlinear.AND.
     $    (continuity=='n=0 only'.OR.continuity=='full')) THEN
        CALL generic_ptr_set(rb%qnd_n0,tb%qnd_n0,tb%tgeom,
     $                       inode,n0,dp,dp,0_i4)
        md_sink=md_sink*(n0(1,:,:)+nd_eq(1,:,:))
        diff=diff*(n0(1,:,:)+nd_eq(1,:,:))
      ELSE
        md_sink=md_sink*nd_eq(1,:,:)
        diff=diff*nd_eq(1,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     evaluate fields for the semi-implicit operator if needed.
c-----------------------------------------------------------------------
      CALL generic_ptr_set(rb%qbe_eq,tb%qbe_eq,tb%tgeom,
     $                     inode,beq,beqr,beqz,1_i4)
      CALL generic_ptr_set(rb%qja_eq,tb%qja_eq,tb%tgeom,
     $                     inode,j0,dp,dp,0_i4)
      IF (nonlinear) THEN
        CALL generic_ptr_set(rb%qbe_n0,tb%qbe_n0,tb%tgeom,inode,
     $                       b0,b0r,b0z,1_i4)
        ja(1,:,:)=j0(1,:,:)+ b0z(3,:,:)/mu0
        ja(2,:,:)=j0(2,:,:)- b0r(3,:,:)/mu0
        ja(3,:,:)=j0(3,:,:)+(b0r(2,:,:)-b0z(1,:,:))/mu0
        IF (geom=='tor') ja(2,:,:)=ja(2,:,:)-b0(3,:,:)/(mu0*bigr)
        be=b0+beq
        ber=b0r+beqr
        bez=b0z+beqz
        CALL generic_ptr_set(rb%qsi_nl_pres,tb%qsi_nl_pres,tb%tgeom,
     $                       inode,nl_pres,dp,dp,0_i4)
      ELSE
        ja=j0
        be=beq
        ber=beqr
        bez=beqz
      ENDIF
      bmagsq=SUM(be*be,1)
      IF (beta>0) THEN
        CALL generic_ptr_set(rb%qpres_eq,tb%qpres_eq,tb%tgeom,
     $                       inode,peq,peqr,peqz,1_i4)
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qpres_n0,tb%qpres_n0,tb%tgeom,
     $                         inode,p0,p0r,p0z,1_i4)
          pres=p0+peq
          presr=p0r+peqr
          presz=p0z+peqz
        ELSE
          pres=peq
          presr=peqr
          presz=peqz
        ENDIF
      ELSE
        pres=0._r8
        presr=0._r8
        presz=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     find the spatial functions used as coefficients for the isotropic
c     and the anisotropic parts of the  operator.
c-----------------------------------------------------------------------
      ani=0.25_r8*si_fac_mhd* (1-mhd_si_iso)*dt**2/mu0
      cj0=0.125_r8*si_fac_j0*dt**2
      IF (beta>0) THEN
        iso=mhd_si_iso*
     $    (si_fac_mhd*bmagsq/mu0+si_fac_pres*gamma*pres(1,:,:))
        grd=0.25_r8*si_fac_pres*(1-mhd_si_iso)*dt**2*gamma*pres(1,:,:)
      ELSE
        iso=mhd_si_iso*si_fac_mhd*bmagsq/mu0
        grd=0._r8
      ENDIF
      IF (nonlinear) THEN
        iso=iso+si_fac_nl*(nl_pres(1,:,:)/mu0+gamma*nl_pres(2,:,:))
      ENDIF
      iso=iso*(0.25*dt**2)
      IF (kin_visc>0.AND.integrand_flag(1:3)/='mhd')
     $  iso=iso+(dt*fvsc*kin_visc)*diff
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     fields for implicit advection.  note that flows and their
c     gradients are multiplied by the mass density and 0.5*dt here.
c-----------------------------------------------------------------------
      IF ((nonlinear.OR.eq_flow/='none').AND.
     $    (advect=='V only'.OR.advect=='all').AND.impladv) THEN
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qve_n0,tb%qve_n0,tb%tgeom,inode,
     $                         v0,v0r,v0z,1_i4)
          ve=v0
          ver=v0r
          vez=v0z
        ELSE
          ve=0._r8
          ver=0._r8
          vez=0._r8
        ENDIF
        IF (eq_flow/='none') THEN
          CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                         inode,veq,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qgrdveq,tb%qgrdveq,tb%tgeom,
     $                         inode,grdv_eq,dp,dp,0_i4)
          ve=ve+veq
          ver=ver+grdv_eq(1:7:3,:,:)
          vez=vez+grdv_eq(2:8:3,:,:)
        ENDIF
        DO i1=1,3
          ve(i1,:,:)=0.5_r8*dt*md_sink*ve(i1,:,:)
          ver(i1,:,:)=0.5_r8*dt*md_sink*ver(i1,:,:)
          vez(i1,:,:)=0.5_r8*dt*md_sink*vez(i1,:,:)
        ENDDO
      ELSE
        ve=0._r8
        ver=0._r8
        vez=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     the anisotropic operator has the form
c     curl(A_iXB0).curl(A_jXB0) after the integration by parts,
c     where A_i is the vector test function alpha*e_hat_i
c     and each curl is broken into B0.grad(A)-A.grad(B0)
c     -B0(div(A)).  B0 is in cylindrical components.
c     the first two indices of the temporary crl matrix are defined as
c     ( curl vector index, operand vector index ), where the vector
c     indices runs from 1 to 3 as (r,z,phi), and the last index is
c     the basis function index.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
      DO iv=1,nvc
        crl(1,1,:,:,iv)=
     $    dalpdz(:,:,iv)*  be(2,:,:)
     $    -alpha(:,:,iv)*((be(1,:,:)*g-(0,1)*keff(jmode)*be(3,:,:))/bigr
     $                   +ber(1,:,:))
        crl(1,2,:,:,iv)=
     $    -(dalpdz(:,:,iv)*be(1,:,:)+alpha(:,:,iv)*bez(1,:,:))
        crl(1,3,:,:,iv)=
     $    -alpha(:,:,iv)*(0,1)*keff(jmode)*be(1,:,:)/bigr
        crl(2,1,:,:,iv)=
     $    -(dalpdr(:,:,iv)* be(2,:,:)
     $      +alpha(:,:,iv)*(be(2,:,:)*g/bigr+ber(2,:,:)))
        crl(2,2,:,:,iv)=
     $    dalpdr(:,:,iv)*be(1,:,:)
     $    +alpha(:,:,iv)*((0,1)*keff(jmode)*be(3,:,:)/bigr-bez(2,:,:))
        crl(2,3,:,:,iv)=
     $    -alpha(:,:,iv)*(0,1)*keff(jmode)*be(2,:,:)/bigr
        crl(3,1,:,:,iv)=
     $    -(dalpdr(:,:,iv)*be(3,:,:)+alpha(:,:,iv)*ber(3,:,:))
        crl(3,2,:,:,iv)=
     $    -(dalpdz(:,:,iv)*be(3,:,:)+alpha(:,:,iv)*bez(3,:,:))
        crl(3,3,:,:,iv)=
     $    dalpdr(:,:,iv)*be(1,:,:)+dalpdz(:,:,iv)*be(2,:,:)
     $    -alpha(:,:,iv)*be(1,:,:)*g/bigr
      ENDDO
c-----------------------------------------------------------------------
c     create a stress tensor that allows for different forms of
c     viscous stress.  the isotropic contributions to the semi-implicit
c     operator now appear here, too.  we will also save the gradient
c     of the scalar function alpha.  the first index of piten runs
c     over rows first, then columns:
c
c       (1-9)=(r_r,z_r,phi_r,r_z,z_z,phi_z,r_phi,z_phi,phi_phi)
c
c-----------------------------------------------------------------------
      DO iv=1,nvc
        gr_alpha(1,:,:,iv)=dalpdr(:,:,iv)
        gr_alpha(2,:,:,iv)=dalpdz(:,:,iv)
        gr_alpha(3,:,:,iv)=(0,1)*keff(jmode)*alpha(:,:,iv)/bigr
        piten(1,1,:,:,iv)=-iso*gr_alpha(1,:,:,iv)
        piten(2,1,:,:,iv)=-iso*gr_alpha(2,:,:,iv)
        piten(3,1,:,:,iv)=-iso*gr_alpha(3,:,:,iv)
        piten(4:9,1,:,:,iv)=0._r8
        piten(1:3,2,:,:,iv)=0._r8
        piten(4:6,2,:,:,iv)=piten(1:3,1,:,:,iv)
        piten(7:9,2,:,:,iv)=0._r8
        piten(1:6,3,:,:,iv)=0._r8
        piten(7:9,3,:,:,iv)=piten(1:3,1,:,:,iv)
        IF (geom=='tor') THEN
          piten(9,1,:,:,iv)=-iso*alpha(:,:,iv)/bigr
          piten(3,3,:,:,iv)= iso*alpha(:,:,iv)/bigr
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     isotropic viscous stress, i.e. proportional to -W, where
c     W is grad(V)+grad(V)^T-2*I*div(V)/3
c-----------------------------------------------------------------------
        IF (iso_visc>0) THEN
          diff=dt*fvsc*iso_visc*diff
          DO iv=1,nvc
           DO iy=1,ncy
           DO ix=1,ncx
            alprs=g*alpha(ix,iy,iv)/bigr(ix,iy)
            dvas=twoth*(gr_alpha(1,ix,iy,iv)+alprs)
            piten(1,1,ix,iy,iv)=piten(1,1,ix,iy,iv)-diff(ix,iy)*
     $        (2._r8*gr_alpha(1,ix,iy,iv)-dvas)
            piten(2,1,ix,iy,iv)=piten(2,1,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(2,ix,iy,iv)
            piten(3,1,ix,iy,iv)=piten(3,1,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(3,ix,iy,iv)
            piten(4,1,ix,iy,iv)=piten(4,1,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(2,ix,iy,iv)
            piten(5,1,ix,iy,iv)=piten(5,1,ix,iy,iv)+diff(ix,iy)*dvas
            piten(7,1,ix,iy,iv)=piten(7,1,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(3,ix,iy,iv)
            piten(9,1,ix,iy,iv)=piten(9,1,ix,iy,iv)-diff(ix,iy)*
     $        (2._r8*alprs-dvas)

            dvas=twoth*gr_alpha(2,ix,iy,iv)
            piten(1,2,ix,iy,iv)=piten(1,2,ix,iy,iv)+diff(ix,iy)*dvas
            piten(2,2,ix,iy,iv)=piten(2,2,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(1,ix,iy,iv)
            piten(4,2,ix,iy,iv)=piten(4,2,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(1,ix,iy,iv)
            piten(5,2,ix,iy,iv)=piten(5,2,ix,iy,iv)-diff(ix,iy)*
     $        (2._r8*gr_alpha(2,ix,iy,iv)-dvas)
            piten(6,2,ix,iy,iv)=piten(6,2,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(3,ix,iy,iv)
            piten(8,2,ix,iy,iv)=piten(8,2,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(3,ix,iy,iv)
            piten(9,2,ix,iy,iv)=piten(9,2,ix,iy,iv)+diff(ix,iy)*dvas

            dvas=twoth*gr_alpha(3,ix,iy,iv)
            piten(1,3,ix,iy,iv)=piten(1,3,ix,iy,iv)+diff(ix,iy)*dvas
            piten(3,3,ix,iy,iv)=piten(3,3,ix,iy,iv)-diff(ix,iy)*
     $        (gr_alpha(1,ix,iy,iv)-alprs)
            piten(5,3,ix,iy,iv)=piten(5,3,ix,iy,iv)+diff(ix,iy)*dvas
            piten(6,3,ix,iy,iv)=piten(6,3,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(2,ix,iy,iv)
            piten(7,3,ix,iy,iv)=piten(7,3,ix,iy,iv)-diff(ix,iy)*
     $        (gr_alpha(1,ix,iy,iv)-alprs)
            piten(8,3,ix,iy,iv)=piten(8,3,ix,iy,iv)-diff(ix,iy)*
     $           gr_alpha(2,ix,iy,iv)
            piten(9,3,ix,iy,iv)=piten(9,3,ix,iy,iv)-diff(ix,iy)*
     $        (2._r8*gr_alpha(3,ix,iy,iv)-dvas)
           ENDDO
           ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       time-centered gyroviscosity operator.  (the div(V) terms
c       cancel analytically).
c-----------------------------------------------------------------------
        bbset=.false.
        IF (gyr_visc>0) THEN
          CALL generic_ptr_set(rb%qti_n0,tb%qti_n0,tb%tgeom,
     $                         inode,ti_sym,dp,dp,0_i4)
          CALL generic_ptr_set(rb%qbb,tb%qbb,tb%tgeom,inode,bb,
     $                         dp,dp,0_i4)
          bbset=.true.
          IF (nonlinear.AND.
     $        (continuity=='n=0 only'.OR.continuity=='full')) THEN
            diff=0.125_r8*gyr_visc*dt*ms(2)
     $          *kboltz*(n0(1,:,:)+nd_eq(1,:,:))
     $          *ti_sym(1,:,:)/(zeff**2*elementary_q*bmagsq)
          ELSE
            diff=0.125_r8*gyr_visc*dt*ms(2)*kboltz*nd_eq(1,:,:)*
     $           ti_sym(1,:,:)/(zeff**2*elementary_q*bmagsq)
          ENDIF
          DO iv=1,nvc
            IF (geom=='tor') THEN
              alpr=alpha(:,:,iv)/bigr
            ELSE
              alpr=0._r8
            ENDIF
            DO iz=1,3
              SELECT CASE(iz)
                CASE(1)
                  ptmp(1,1,:,:)=2._r8*gr_alpha(1,:,:,iv)
                  ptmp(2,1,:,:)=gr_alpha(2,:,:,iv)
                  ptmp(3,1,:,:)=gr_alpha(3,:,:,iv)
                  ptmp(1,2,:,:)=gr_alpha(2,:,:,iv)
                  ptmp(2:3,2,:,:)=0._r8
                  ptmp(1,3,:,:)=gr_alpha(3,:,:,iv)
                  ptmp(2,3,:,:)=0._r8
                  ptmp(3,3,:,:)=2._r8*alpr
                CASE(2)
                  ptmp(1,1,:,:)=0._r8
                  ptmp(2,1,:,:)=gr_alpha(1,:,:,iv)
                  ptmp(3,1,:,:)=0._r8
                  ptmp(1,2,:,:)=gr_alpha(1,:,:,iv)
                  ptmp(2,2,:,:)=2._r8*gr_alpha(2,:,:,iv)
                  ptmp(3,2,:,:)=gr_alpha(3,:,:,iv)
                  ptmp(1,3,:,:)=0._r8
                  ptmp(2,3,:,:)=gr_alpha(3,:,:,iv)
                  ptmp(3,3,:,:)=0._r8
                CASE(3)
                  ptmp(1:2,1,:,:)=0._r8
                  ptmp(3,1,:,:)=gr_alpha(1,:,:,iv)-alpr
                  ptmp(1:2,2,:,:)=0._r8
                  ptmp(3,2,:,:)=gr_alpha(2,:,:,iv)
                  ptmp(1,3,:,:)=gr_alpha(1,:,:,iv)-alpr
                  ptmp(2,3,:,:)=gr_alpha(2,:,:,iv)
                  ptmp(3,3,:,:)=2._r8*gr_alpha(3,:,:,iv)
              END SELECT
              DO iy=1,ncy
                DO ix=1,ncx
                  gcoef=diff(ix,iy)
                  wdot(:,1)=ptmp(:,1,ix,iy)+
     $                     3._r8*(ptmp(:,1,ix,iy)*bb(1,ix,iy)+
     $                            ptmp(:,2,ix,iy)*bb(4,ix,iy)+
     $                            ptmp(:,3,ix,iy)*bb(5,ix,iy))
                  wdot(:,2)=ptmp(:,2,ix,iy)+
     $                     3._r8*(ptmp(:,1,ix,iy)*bb(4,ix,iy)+
     $                            ptmp(:,2,ix,iy)*bb(2,ix,iy)+
     $                            ptmp(:,3,ix,iy)*bb(6,ix,iy))
                  wdot(:,3)=ptmp(:,3,ix,iy)+
     $                     3._r8*(ptmp(:,1,ix,iy)*bb(5,ix,iy)+
     $                            ptmp(:,2,ix,iy)*bb(6,ix,iy)+
     $                            ptmp(:,3,ix,iy)*bb(3,ix,iy))
                  bc_wdot(1,:)=be(2,ix,iy)*wdot(3,:)-
     $                         be(3,ix,iy)*wdot(2,:)
                  bc_wdot(2,:)=be(3,ix,iy)*wdot(1,:)-
     $                         be(1,ix,iy)*wdot(3,:)
                  bc_wdot(3,:)=be(1,ix,iy)*wdot(2,:)-
     $                         be(2,ix,iy)*wdot(1,:)
                  piten(1:3,iz,ix,iy,iv)=piten(1:3,iz,ix,iy,iv)+
     $                         gcoef*(bc_wdot(:,1)+bc_wdot(1,:))
                  piten(4:6,iz,ix,iy,iv)=piten(4:6,iz,ix,iy,iv)+
     $                         gcoef*(bc_wdot(:,2)+bc_wdot(2,:))
                  piten(7:9,iz,ix,iy,iv)=piten(7:9,iz,ix,iy,iv)+
     $                         gcoef*(bc_wdot(:,3)+bc_wdot(3,:))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       parallel viscous stress.
c-----------------------------------------------------------------------
        IF (par_visc>0) THEN
          IF (.NOT.bbset) CALL generic_ptr_set(rb%qbb,tb%qbb,tb%tgeom,
     $                                         inode,bb,dp,dp,0_i4)
          IF (nonlinear.AND.
     $        (continuity=='n=0 only'.OR.continuity=='full')) THEN
            diff=-3._r8*dt*fvsc*par_visc*mtot*(n0(1,:,:)+nd_eq(1,:,:))
          ELSE
            diff=-3._r8*dt*fvsc*par_visc*mtot*nd_eq(1,:,:)
          ENDIF
          IF (parvisc_model=='plltdep') THEN
            CALL generic_ptr_set(rb%qkappli_n0,tb%qkappli_n0,tb%tgeom,
     $                           inode,kappli,dp,dp,0_i4)
            diff=diff*kappli(1,:,:)/k_plli
          ENDIF
          DO iv=1,nvc
            IF (geom=='tor') THEN
              alpr=alpha(:,:,iv)/bigr
            ELSE
              alpr=0._r8
            ENDIF
            DO iy=1,ncy
              DO ix=1,ncx
                b2=-third*SUM(bb(1:3,ix,iy))
                ctmp(1)=(gr_alpha(1,ix,iy,iv)*(bb(1,ix,iy)+b2)+
     $                   gr_alpha(2,ix,iy,iv)* bb(4,ix,iy)+
     $                   gr_alpha(3,ix,iy,iv)* bb(5,ix,iy)+
     $                            alpr(ix,iy)*(bb(3,ix,iy)+b2))*
     $                  diff(ix,iy)
                ctmp(2)=(gr_alpha(1,ix,iy,iv)* bb(4,ix,iy)+
     $                   gr_alpha(2,ix,iy,iv)*(bb(2,ix,iy)+b2)+
     $                   gr_alpha(3,ix,iy,iv)* bb(6,ix,iy))*
     $                  diff(ix,iy)
                ctmp(3)=((gr_alpha(1,ix,iy,iv)-alpr(ix,iy))*bb(5,ix,iy)+
     $                    gr_alpha(2,ix,iy,iv)* bb(6,ix,iy)+
     $                    gr_alpha(3,ix,iy,iv)*(bb(3,ix,iy)+b2))*
     $                  diff(ix,iy)
                piten(1,:,ix,iy,iv)=piten(1,:,ix,iy,iv)+
     $                              ctmp(:)*(bb(1,ix,iy)+b2)
                piten(2,:,ix,iy,iv)=piten(2,:,ix,iy,iv)+
     $                              ctmp(:)* bb(4,ix,iy)
                piten(3,:,ix,iy,iv)=piten(3,:,ix,iy,iv)+
     $                              ctmp(:)* bb(5,ix,iy)
                piten(4,:,ix,iy,iv)=piten(4,:,ix,iy,iv)+
     $                              ctmp(:)* bb(4,ix,iy)
                piten(5,:,ix,iy,iv)=piten(5,:,ix,iy,iv)+
     $                              ctmp(:)*(bb(2,ix,iy)+b2)
                piten(6,:,ix,iy,iv)=piten(6,:,ix,iy,iv)+
     $                              ctmp(:)* bb(6,ix,iy)
                piten(7,:,ix,iy,iv)=piten(7,:,ix,iy,iv)+
     $                              ctmp(:)* bb(5,ix,iy)
                piten(8,:,ix,iy,iv)=piten(8,:,ix,iy,iv)+
     $                              ctmp(:)* bb(6,ix,iy)
                piten(9,:,ix,iy,iv)=piten(9,:,ix,iy,iv)+
     $                              ctmp(:)*(bb(3,ix,iy)+b2)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

c-----------------------------------------------------------------------
c     from here, gr_alpha will be used for the test function, so we
c     perform the conjugate operation.
c-----------------------------------------------------------------------
      gr_alpha=CONJG(gr_alpha)
c-----------------------------------------------------------------------
c     compute isotropic, anisotropic, and pressure contributions
c     together for efficiency.
c
c     the "ani" contributions are the dot product of the curls
c     from the different vertices around each cell and the
c
c     symmetrized contributions from the J0Xcurl(vxB0) term are
c     the "cj0" terms without pressure gradient coefficients.
c
c     grad(P0) terms include the div(A)*div() part of the operator
c     and the symmetrized div(A)*A.grad(P0) part.
c
c     the mass density term uses the md_sink factor and is modified
c     to include the correction factor for artificial particle
c     diffusivity.
c-----------------------------------------------------------------------
      IF (continuity/='none'.AND.nonlinear.AND.impladv.AND.
     $    nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0).AND.
     $    advect/='none') THEN
        CALL generic_ptr_set(rb%qndiff_n0,tb%qndiff_n0,tb%tgeom,
     $                       inode,ndiff,dp,dp,0_i4)
        md_sink=md_sink+0.5_r8*dt*mtot*ndiff(1,:,:)
      ENDIF
c-----------------------------------------------------------------------
c     symmetric contributions
c-----------------------------------------------------------------------
      DO iv=1,nvc
        DO jv=1,nvc
         DO iy=1,ncy
         int(:  ,1,iy,jv,iv)=0._r8
         int(2:3,2,iy,jv,iv)=0._r8
         int(3  ,3,iy,jv,iv)=0._r8
         DO ix=1,ncx
          iso_diag=md_sink(ix,iy)*alpha(ix,iy,iv)*alpha(ix,iy,jv)

          int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+iso_diag+
     $        ani*(CONJG(crl(1,1,ix,iy,iv))*crl(1,1,ix,iy,jv)+
     $             CONJG(crl(2,1,ix,iy,iv))*crl(2,1,ix,iy,jv)+
     $             CONJG(crl(3,1,ix,iy,iv))*crl(3,1,ix,iy,jv))
     $      -cj0*(alpha(ix,iy,iv)*(-ja(3,ix,iy)*crl(2,1,ix,iy,jv)
     $                             +ja(2,ix,iy)*crl(3,1,ix,iy,jv)
     $                           -dalpdrc(ix,iy,jv)*presr(1,ix,iy))
     $           +alpha(ix,iy,jv)*(-ja(3,ix,iy)*CONJG(crl(2,1,ix,iy,iv))
     $                             +ja(2,ix,iy)*CONJG(crl(3,1,ix,iy,iv))
     $                             -dalpdrc(ix,iy,iv)*presr(1,ix,iy)))
     $       +grd(ix,iy)*dalpdrc(ix,iy,iv)*dalpdrc(ix,iy,jv)

          int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)+
     $        ani*(CONJG(crl(1,1,ix,iy,iv))*crl(1,2,ix,iy,jv)+
     $             CONJG(crl(2,1,ix,iy,iv))*crl(2,2,ix,iy,jv)+
     $             CONJG(crl(3,1,ix,iy,iv))*crl(3,2,ix,iy,jv))
     $      -cj0*(alpha(ix,iy,iv)*(-ja(3,ix,iy)*crl(2,2,ix,iy,jv)
     $                             +ja(2,ix,iy)*crl(3,2,ix,iy,jv)
     $                           -dalpdz(ix,iy,jv)*presr(1,ix,iy))
     $           +alpha(ix,iy,jv)*( ja(3,ix,iy)*CONJG(crl(1,1,ix,iy,iv))
     $                             -ja(1,ix,iy)*CONJG(crl(3,1,ix,iy,iv))
     $                            -dalpdrc(ix,iy,iv)*presz(1,ix,iy)))
     $       +grd(ix,iy)*dalpdrc(ix,iy,iv)*dalpdz(ix,iy,jv)


          int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)+
     $        ani*(CONJG(crl(1,1,ix,iy,iv))*crl(1,3,ix,iy,jv)+
     $             CONJG(crl(2,1,ix,iy,iv))*crl(2,3,ix,iy,jv)+
     $             CONJG(crl(3,1,ix,iy,iv))*crl(3,3,ix,iy,jv))
     $      -cj0*(alpha(ix,iy,iv)*(-ja(3,ix,iy)*crl(2,3,ix,iy,jv)
     $                             +ja(2,ix,iy)*crl(3,3,ix,iy,jv)
     $                           +gr_alpha(3,ix,iy,jv)*presr(1,ix,iy))
     $           +alpha(ix,iy,jv)*(-ja(2,ix,iy)*CONJG(crl(1,1,ix,iy,iv))
     $                           +ja(1,ix,iy)*CONJG(crl(2,1,ix,iy,iv))))
     $       -grd(ix,iy)*dalpdrc(ix,iy,iv)*gr_alpha(3,ix,iy,jv)


          int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+iso_diag+
     $        ani*(CONJG(crl(1,2,ix,iy,iv))*crl(1,2,ix,iy,jv)+
     $             CONJG(crl(2,2,ix,iy,iv))*crl(2,2,ix,iy,jv)+
     $             CONJG(crl(3,2,ix,iy,iv))*crl(3,2,ix,iy,jv))
     $       -cj0*(alpha(ix,iy,iv)*( ja(3,ix,iy)*crl(1,2,ix,iy,jv)
     $                              -ja(1,ix,iy)*crl(3,2,ix,iy,jv)
     $                            -dalpdz(ix,iy,jv)*presz(1,ix,iy))
     $           +alpha(ix,iy,jv)*( ja(3,ix,iy)*CONJG(crl(1,2,ix,iy,iv))
     $                             -ja(1,ix,iy)*CONJG(crl(3,2,ix,iy,iv))
     $                            -dalpdz(ix,iy,iv)*presz(1,ix,iy)))
     $       +grd(ix,iy)*dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)

          int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv)+
     $        ani*(CONJG(crl(1,2,ix,iy,iv))*crl(1,3,ix,iy,jv)+
     $             CONJG(crl(2,2,ix,iy,iv))*crl(2,3,ix,iy,jv)+
     $             CONJG(crl(3,2,ix,iy,iv))*crl(3,3,ix,iy,jv))
     $       -cj0*(alpha(ix,iy,iv)*( ja(3,ix,iy)*crl(1,3,ix,iy,jv)
     $                              -ja(1,ix,iy)*crl(3,3,ix,iy,jv)
     $                            +gr_alpha(3,ix,iy,jv)*presz(1,ix,iy))
     $           +alpha(ix,iy,jv)*(-ja(2,ix,iy)*CONJG(crl(1,2,ix,iy,iv))
     $                           +ja(1,ix,iy)*CONJG(crl(2,2,ix,iy,iv))))
     $       -grd(ix,iy)*dalpdz(ix,iy,iv)*gr_alpha(3,ix,iy,jv)


          int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)+iso_diag+
     $        ani*(CONJG(crl(1,3,ix,iy,iv))*crl(1,3,ix,iy,jv)+
     $             CONJG(crl(2,3,ix,iy,iv))*crl(2,3,ix,iy,jv)+
     $             CONJG(crl(3,3,ix,iy,iv))*crl(3,3,ix,iy,jv))
     $       -cj0*(alpha(ix,iy,iv)*(-ja(2,ix,iy)*crl(1,3,ix,iy,jv)
     $                              +ja(1,ix,iy)*crl(2,3,ix,iy,jv))
     $           +alpha(ix,iy,jv)*(-ja(2,ix,iy)*CONJG(crl(1,3,ix,iy,iv))
     $                           +ja(1,ix,iy)*CONJG(crl(2,3,ix,iy,iv))))
     $       +grd(ix,iy)*alpha(ix,iy,iv)*alpha(ix,iy,jv)*
     $                   k2ef(jmode)/bigr2(ix,iy)

         ENDDO
         int(1,2,iy,iv,jv)=CONJG(int(2,1,iy,jv,iv))
         int(1,3,iy,iv,jv)=CONJG(int(3,1,iy,jv,iv))
         int(2,3,iy,iv,jv)=CONJG(int(3,2,iy,jv,iv))
         ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     possibly asymmetric stresses and implicit advection.  dvas and
c     iso_diag are used as temporary storage.
c
c     implicit advection includes Ai*.V.grad(Aj)+Ai*.Aj.grad(V), where
c     V is the n=0 part of the solution plus equilibrium center-of-mass
c     flow velocity.  [V is also multiplied by the symmetric part of
c     rho above, so it's really momentum density.]  Ai* and Aj are the
c     vector test- and basis-functions, where i and j represent the
c     finite element expansion-function row and column indices.  thus,
c     Aj are the basis functions for the delta-V solution vector.
c
c     with ei being the unit direction vector and alpha_i being
c     the scalar expansion function for the i-th basis, the
c     Ai*.V.grad(Aj) term includes:
c
c       (alpha_i*)V.grad(alpha_j)(ei.ej)  --  the dvas terms below
c
c       -(Vphi)(alpha_i*)(alpha_j)(ei.er)(ej.ephi)/r  --  in int(3,1)
c
c       +(Vphi)(alpha_i*)(alpha_j)(ei.ephi)(ej.er)/r  --  in int(1,3)
c
c     Ai*.Aj.grad(V) term includes:
c
c       (alpha_i*)(alpha_j)(ej.grad(V.ei))  --  ver & vez terms
c
c       -(Vphi)(alpha_i*)(alpha_j)(ei.er)(ej.ephi)/r  --  in int(3,1)
c
c       +(Vr)(alpha_i*)(alpha_j)(ei.ephi)(ej.ephi)/r  --  in int(3,3)
c-----------------------------------------------------------------------
      DO iv=1,nvc
        DO jv=1,nvc
         DO iy=1,ncy
         DO ix=1,ncx
          alprs=g*alpha(ix,iy,iv)/bigr(ix,iy)
          iso_diag=alpha(ix,iy,iv)*alpha(ix,iy,jv)
          dvas=ve(1,ix,iy)*dalpdr(ix,iy,jv)+
     $         ve(2,ix,iy)*dalpdz(ix,iy,jv)+
     $         ve(3,ix,iy)*(0,1)*keff(jmode)*alpha(ix,iy,jv)/bigr(ix,iy)

          int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(1,1,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(2,1,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(3,1,ix,iy,jv)+
     $         alprs*piten(9,1,ix,iy,jv))
     $       +alpha(ix,iy,iv)*(dvas+alpha(ix,iy,jv)*ver(1,ix,iy))

          int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(1,2,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(2,2,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(3,2,ix,iy,jv)+
     $         alprs*piten(9,2,ix,iy,jv))
     $       +iso_diag*vez(1,ix,iy)

          int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(1,3,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(2,3,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(3,3,ix,iy,jv)+
     $         alprs*(piten(9,3,ix,iy,jv)+
     $                2.*ve(3,ix,iy)*alpha(ix,iy,jv)))

          int(1,2,iy,jv,iv)=int(1,2,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(4,1,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(5,1,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(6,1,ix,iy,jv))
     $       +iso_diag*ver(2,ix,iy)

          int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(4,2,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(5,2,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(6,2,ix,iy,jv))
     $       +alpha(ix,iy,iv)*(dvas+alpha(ix,iy,jv)*vez(2,ix,iy))

          int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(4,3,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(5,3,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(6,3,ix,iy,jv))

          int(1,3,iy,jv,iv)=int(1,3,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(7,1,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(8,1,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(9,1,ix,iy,jv))
     $       +alprs*(piten(3,1,ix,iy,jv)+ve(3,ix,iy)*alpha(ix,iy,jv))
     $       +iso_diag*ver(3,ix,iy)

          int(2,3,iy,jv,iv)=int(2,3,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(7,2,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(8,2,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(9,2,ix,iy,jv))
     $       +alprs*piten(3,2,ix,iy,jv)
     $       +iso_diag*vez(3,ix,iy)

          int(3,3,iy,jv,iv)=int(3,3,iy,jv,iv)
     $       -(gr_alpha(1,ix,iy,iv)*piten(7,3,ix,iy,jv)+
     $         gr_alpha(2,ix,iy,iv)*piten(8,3,ix,iy,jv)+
     $         gr_alpha(3,ix,iy,iv)*piten(9,3,ix,iy,jv))
     $       +alprs*(piten(3,3,ix,iy,jv)+ve(1,ix,iy)*alpha(ix,iy,jv))
     $       +alpha(ix,iy,iv)*dvas
         ENDDO
         ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find the contributions associated with the rows and columns of the
c     discontinuous scalar for stabilizing the compression with an
c     auxiliary field.  the equations are expressed directly in weak
c     form:
c
c       int[rho*A^*.delta(V)] 
c         + int[ sqrt(ddivv*nu*dt*f)*auxv*div(A^*) ]=dt*int[A^*.forces]
c       int[(ups^*)*auxv]
c         - int[ sqrt(ddivv*nu*dt*f)*(ups^*)*div(delta(V))]
c         = int[ sqrt(ddivv*nu*dt/f)*(ups^*)*div(V_old)]
c
c     where cma is sqrt(rho) times the magneto-acoustic speed, ups
c     is the scalar test field for the auxiliary field auxv, f is the
c     implicit centering and nu is a numerical viscosity coefficient,
c     nu=dt*cma**2.
c
c     md_sink is redefined as sqrt(rho) times the magneto-
c     acoustic speed, multiplied by numerical factors of dt and the
c     divv coefficient.
c
c     note that gr_alpha here is grad(alpha^*).
c
c     find the contributions associated with the rows and columns of the
c     discontinuous scalar for stabilizing parallel vorticity using an
c     auxiliary field.
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
c     here we create the time-centered operator for the implicit part
c     of the discretized system.
c
c     the implicit centering parameter is fpvrt.
c-----------------------------------------------------------------------
      IF (poly_divv>=0) THEN
        IF (beta>0) THEN
          md_sink=dt*SQRT(fdivv*ddivv)*
     $            SQRT(MAX(0._r8,gamma*pres(1,:,:))+SUM(be**2,1)/mu0)
        ELSE
          md_sink=dt*SQRT(fdivv*ddivv)*SQRT(SUM(be**2,1)/mu0)
        ENDIF
        dtpv=dt*SQRT(fpvrt*dpvrt/mu0)

        DO iv=1,nvc
          DO jv=nvc+1,nv
            jvm=jv-nvc
            DO iy=1,ncy
              int(:,:,iy,jv,iv)=0._r8
              DO ix=1,ncx
                gcoef=bigr(ix,iy)
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+
     $            md_sink(ix,iy)*dalpdrc(ix,iy,iv)*almod(ix,iy,jvm)
                int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv) - dtpv*
     $            almod(ix,iy,jvm)*(
     $              be(2,ix,iy)*(0,1)*alpha(ix,iy,iv)*keff(jmode)/gcoef+
     $              be(3,ix,iy)*dalpdz(ix,iy,iv) )
                int(1,2,iy,jv,iv)=int(1,2,iy,jv,iv)+
     $            md_sink(ix,iy)*dalpdz(ix,iy,iv)*almod(ix,iy,jvm)
                int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv) + dtpv*
     $            almod(ix,iy,jvm)*(
     $              be(1,ix,iy)*(0,1)*alpha(ix,iy,iv)*keff(jmode)/gcoef+
     $              be(3,ix,iy)*dalpdr(ix,iy,iv) )
                int(1,3,iy,jv,iv)=int(1,3,iy,jv,iv)-(0,1)*keff(jmode)*
     $            md_sink(ix,iy)*alpha(ix,iy,iv)*almod(ix,iy,jvm)/gcoef
                int(2,3,iy,jv,iv)=int(2,3,iy,jv,iv) - dtpv*
     $            almod(ix,iy,jvm)*(
     $              be(2,ix,iy)*dalpdrc(ix,iy,iv)-
     $              be(1,ix,iy)*dalpdz (ix,iy,iv) )
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DO iv=nvc+1,nv
          ivm=iv-nvc
          DO jv=1,nvc
            DO iy=1,ncy
              int(:,:,iy,jv,iv)=0._r8
              DO ix=1,ncx
                gcoef=bigr(ix,iy)
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)-
     $            md_sink(ix,iy)*dalpdrc(ix,iy,jv)*almod(ix,iy,ivm)
                int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-
     $            md_sink(ix,iy)*dalpdz(ix,iy,jv)*almod(ix,iy,ivm)
                int(3,1,iy,jv,iv)=int(3,1,iy,jv,iv)-(0,1)*keff(jmode)*
     $            md_sink(ix,iy)*alpha(ix,iy,jv)*almod(ix,iy,ivm)/gcoef
                int(1,2,iy,jv,iv)=int(1,2,iy,jv,iv) - dtpv*
     $            almod(ix,iy,ivm)*(
     $              be(2,ix,iy)*(0,1)*alpha(ix,iy,jv)*keff(jmode)/gcoef-
     $              be(3,ix,iy)*dalpdz(ix,iy,jv) )
                int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv) + dtpv*
     $            almod(ix,iy,ivm)*(
     $              be(1,ix,iy)*(0,1)*alpha(ix,iy,jv)*keff(jmode)/gcoef-
     $              be(3,ix,iy)*dalpdr(ix,iy,jv) )
                int(3,2,iy,jv,iv)=int(3,2,iy,jv,iv) + dtpv*
     $            almod(ix,iy,ivm)*(
     $              be(2,ix,iy)*dalpdrc(ix,iy,jv)-
     $              be(1,ix,iy)*dalpdz (ix,iy,jv) )
              ENDDO
            ENDDO
          ENDDO
          DO jv=nvc+1,nv
            jvm=jv-nvc
            DO iy=1,ncy
              int(:,:,iy,jv,iv)=0._r8
              DO ix=1,ncx
                int(1,1,iy,jv,iv)=int(1,1,iy,jv,iv)+
     $            almod(ix,iy,jvm)*almod(ix,iy,ivm)
                int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+
     $            almod(ix,iy,jvm)*almod(ix,iy,ivm)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE v_aniso_op
c-----------------------------------------------------------------------
c     subprogram 10. grad_div.
c     compute the integrand for the lhs of the time split div(b)
c     diffuser for a single mode.
c-----------------------------------------------------------------------
      SUBROUTINE grad_div(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2
      INTEGER(i4) :: iv,jv,nv
c-----------------------------------------------------------------------
c     set the rblock or tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nv=SIZE(int,5)
c-----------------------------------------------------------------------
c     evaluate div(alpha)*div() for divergence cleaning.
c     note that the three components of B are either 
c     (real_r,real_z,-imag_phi) or (imag_r,imag_z,real_phi).
c-----------------------------------------------------------------------
      bigr2=bigr**2
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM(alpha(:,:,iv)*alpha(:,:,jv)+
     $       dt*fdivb*divbd*dalpdrc(:,:,iv)*dalpdrc(:,:,jv),1)
          int(2,1,:,jv,iv)=SUM(
     $       dt*fdivb*divbd*dalpdrc(:,:,iv)*dalpdz(:,:,jv),1)
          int(1,2,:,iv,jv)=int(2,1,:,jv,iv)
          int(3,1,:,jv,iv)=SUM(
     $       dt*fdivb*divbd*dalpdrc(:,:,iv)
     $                     *  alpha(:,:,jv)*keff(jmode)/bigr,1)
          int(1,3,:,iv,jv)=int(3,1,:,jv,iv)

          int(2,2,:,jv,iv)=SUM(alpha(:,:,iv)*alpha(:,:,jv)+
     $       dt*fdivb*divbd*dalpdz(:,:,iv)*dalpdz(:,:,jv),1)
          int(3,2,:,jv,iv)=SUM(
     $       dt*fdivb*divbd*dalpdz(:,:,iv)
     $                     * alpha(:,:,jv)*keff(jmode)/bigr,1)
          int(2,3,:,iv,jv)=int(3,2,:,jv,iv)

          int(3,3,:,jv,iv)=SUM(alpha(:,:,iv)*alpha(:,:,jv)+
     $       dt*fdivb*divbd*alpha(:,:,iv)*alpha(:,:,jv)*k2ef(jmode)
     $                     /bigr2,1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE grad_div
c-----------------------------------------------------------------------
c     subprogram 11. cont_op.
c     compute a nonsymmetric advective/diffusion operator for number 
c     density.  integration by parts is applied to both terms.
c-----------------------------------------------------------------------
      SUBROUTINE cont_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          ve_eq,nd_eq,v0,dart,upave,vv,dp
      REAL(r8), DIMENSION(3,SIZE(bigr,1),SIZE(bigr,2)) :: ve
      REAL(r8), DIMENSION(  SIZE(bigr,1),SIZE(bigr,2)) :: bigr2,dfac,pll

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: nd,dcp

      INTEGER(i4) :: nv,iv,jv,iq,jq,ix,iy,ncx,ncy
      REAL(r8) :: fhyp,alpij
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      IF (impladv) THEN
        IF (nonlinear) THEN
          CALL generic_ptr_set(rb%qve_n0,tb%qve_n0,tb%tgeom,
     $                         inode,v0,dp,dp,0_i4)
          ve=0.5_r8*dt*v0
        ELSE
          ve=0._r8
        ENDIF
        IF (eq_flow/='none') THEN
          CALL generic_ptr_set(rb%qve_eq,tb%qve_eq,tb%tgeom,
     $                         inode,ve_eq,dp,dp,0_i4)
          ve=ve+0.5_r8*dt*ve_eq
        ENDIF
        fhyp=SQRT(nd_hypd*dt)
      ELSE
        ve=0._r8
      ENDIF
      nv=SIZE(int,5)
      ncx=SIZE(bigr,1)
      ncy=SIZE(bigr,2)
      bigr2=bigr**2
c-----------------------------------------------------------------------
c     artificial particle diffusivities.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.nd_floor>0.AND.nd_diff>0) THEN
        CALL generic_ptr_set(rb%qdart,tb%qdart,tb%tgeom,
     $                       inode,dart,dp,dp,0_i4)
        dfac=fthc*dt*nd_diff+dart(1,:,:)
      ELSE
        dfac=fthc*dt*nd_diff
      ENDIF
      IF (nonlinear.AND.nd_dart_upw>0) THEN
        CALL generic_ptr_set(rb%qupw_n0,tb%qupw_n0,tb%tgeom,
     $                       inode,upave,dp,dp,0_i4)
        CALL generic_ptr_set(rb%qvv,tb%qvv,tb%tgeom,
     $                       inode,vv,dp,dp,0_i4)
        pll=nd_dart_upw*dt**2*upave(1,:,:)*upw_aniso/
     $                  ( SUM(vv(1:3,:,:),1) + smallnum )
        dfac=dfac+nd_dart_upw*dt**2*upave(1,:,:)*(1._r8-upw_aniso)
      ENDIF
c-----------------------------------------------------------------------
c     Find -grad(alpha(i)).V*alpha(j)
c          -nd_diff*grad(alpha(i)).grad(alpha(j)).
c     fthc is the centering parameter for diffusion, but implicit
c     advection is time-centered.
c
c     the second part is for anisotropic upwinding-like diffusion.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     the third part applies a fourth-order hyper-diffusivity using
c     an auxiliary field.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
           int(1,1,:,jv,iv)=SUM(-(ve(1,:,:)*dalpdr(:,:,iv)
     $                           +ve(2,:,:)*dalpdz(:,:,iv)
     $                           -ve(3,:,:)*alpha(:,:,iv)*
     $                                    (0,1)*keff(jmode)/bigr)*
     $                         alpha(:,:,jv)
     $                        +(dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $                         +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $                         +k2ef(jmode)*alpha(:,:,iv)*alpha(:,:,jv)
     $                         /bigr2)*dfac
     $                        + alpha(:,:,iv)*alpha(:,:,jv),1)

          IF (.NOT.(nonlinear.AND.nd_dart_upw>0)) CYCLE

            int(1,1,:,jv,iv)=int(1,1,:,jv,iv)+SUM(
     $       ( dalpdr(:,:,iv)*dalpdr(:,:,jv)*vv(1,:,:)
     $        +dalpdz(:,:,iv)*dalpdz(:,:,jv)*vv(2,:,:)
     $        + alpha(:,:,iv)* alpha(:,:,jv)*k2ef(jmode)*vv(3,:,:)/bigr2
     $        + (dalpdr(:,:,iv)*dalpdz(:,:,jv)
     $          +dalpdz(:,:,iv)*dalpdr(:,:,jv))*vv(4,:,:)
     $      +(0,1)*(
     $       (dalpdr(:,:,iv)*alpha(:,:,jv)-dalpdr(:,:,jv)*alpha(:,:,iv))
     $                              *keff(jmode)*vv(5,:,:)/bigr+
     $       (dalpdz(:,:,iv)*alpha(:,:,jv)-dalpdz(:,:,jv)*alpha(:,:,iv))
     $                              *keff(jmode)*vv(6,:,:)/bigr))*pll,1)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     add contributions for a fourth-order hyper-diffusivity using
c     an auxiliary field.
c-----------------------------------------------------------------------
      IF (nd_hypd>0.AND.impladv) THEN
        DO iv=1,nv
          DO jv=1,nv
            DO iy=1,ncy
              int(2,1,iy,jv,iv)=0._r8
              int(2,2,iy,jv,iv)=0._r8
              DO ix=1,ncx
                alpij=alpha(ix,iy,iv)*alpha(ix,iy,jv)
                int(2,1,iy,jv,iv)=int(2,1,iy,jv,iv)-fhyp*
     $                         (dalpdr(ix,iy,iv)*dalpdr(ix,iy,jv)
     $                         +dalpdz(ix,iy,iv)*dalpdz(ix,iy,jv)
     $                         +k2ef(jmode)*alpij/bigr2(ix,iy))
                int(2,2,iy,jv,iv)=int(2,2,iy,jv,iv)+alpij
              ENDDO
              int(1,2,iy,jv,iv)=-int(2,1,iy,jv,iv)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cont_op
c-----------------------------------------------------------------------
c     subprogram 12. vec_lap_op2.
c     construct the integrand for a basic vector Laplacian operator.
c-----------------------------------------------------------------------
      SUBROUTINE vec_lap_op2(int,bigr,rb,dx,dy,tb,inode)

      REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz,dalpdrc
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2

      INTEGER(i4) :: iv,jv,nv
      REAL(r8) :: g
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree,dalpdrc)
      nv=SIZE(int,5)
      bigr2=bigr**2
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     construct the integrand, grad(Ai^*).grad(Aj), where the
c     vector As are the basis/test functions.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM(
     $           dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $          +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $          +(k2ef(jmode)+g)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2 ,1)
          int(2,1,:,jv,iv)=0._r8
          int(3,1,:,jv,iv)=SUM(g*
     $         2._r8*keff(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2,1)
          int(1,2,:,jv,iv)=0._r8
          int(2,2,:,jv,iv)=SUM(
     $           dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $          +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $          +k2ef(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2 ,1)
          int(3,2,:,jv,iv)=0._r8
          int(1,3,:,jv,iv)=int(3,1,:,jv,iv)
          int(2,3,:,jv,iv)=0._r8
          int(3,3,:,jv,iv)=int(1,1,:,jv,iv)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vec_lap_op2
c-----------------------------------------------------------------------
c     subprogram 13. hyp_visc_op.
c     the matrix created through this routine applies hyper-viscous
c     diffusion as a composite of two Laplacians.  each Laplacian
c     only couples (real-r,real-z,-imag-phi) or (imag-r,imag-z,real-phi)
c     components, and the equation for delta(V) and the auxiliary field
c     are packed into one 3x3 complex system that can be applied to
c     either of those real delta(V)-vector systems.  the resulting
c     matrix is symmetric and complex, not Hermitian.
c-----------------------------------------------------------------------
      SUBROUTINE hyp_visc_op(int,bigr,rb,dx,dy,tb,inode)

      COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: int
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
      TYPE(rblock_type), INTENT(INOUT) :: rb
      REAL(r8), INTENT(IN) :: dx,dy
      TYPE(tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: inode

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          alpha,dalpdr,dalpdz
      REAL(r8), DIMENSION(SIZE(bigr,1),SIZE(bigr,2)) :: bigr2
      REAL(r8), DIMENSION(1,1,1) :: dv

      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: nv,iv,jv,nq,np
      COMPLEX(r8) :: cfac
      REAL(r8) :: g
c-----------------------------------------------------------------------
c     set rblock and tblock weight functions.
c-----------------------------------------------------------------------
      CALL generic_alpha_eval(rb,tb%tgeom,inode,'mat',alpha,dalpdr,
     $                        dalpdz,1_i4,poly_degree)
      nv=SIZE(int,5)
      bigr2=bigr**2
      IF (geom=='tor') THEN
        g=1._r8
      ELSE
        g=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     coefficient for the differential operation.  the imaginary part
c     couples the delta-V and Laplacian(V) components.
c-----------------------------------------------------------------------
      cfac=(0,1)*SQRT(hyp_visc*fhyp_visc*dt)
c-----------------------------------------------------------------------
c     toroidal and straight geometries include
c     grad(alph_i*).grad(alpha_j) for each of the three vector
c     components.  toroidal geometry has additional r/phi terms.
c-----------------------------------------------------------------------
      DO iv=1,nv
        DO jv=1,nv
          int(1,1,:,jv,iv)=SUM( cfac*
     $           ( dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $            +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $            +(k2ef(jmode)+g)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2)
     $           +alpha(:,:,iv)*alpha(:,:,jv),1)
          int(2,1,:,jv,iv)=0._r8
          int(3,1,:,jv,iv)=SUM(g*cfac*(0,1)*
     $           2._r8*keff(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2,1)
          int(1,2,:,jv,iv)=0._r8
          int(2,2,:,jv,iv)=SUM( cfac*
     $           ( dalpdr(:,:,iv)*dalpdr(:,:,jv)
     $            +dalpdz(:,:,iv)*dalpdz(:,:,jv)
     $            +k2ef(jmode)*alpha(:,:,iv)*alpha(:,:,jv)/bigr2)
     $           +alpha(:,:,iv)*alpha(:,:,jv),1)
          int(3,2,:,jv,iv)=0._r8
          int(1,3,:,jv,iv)=-int(3,1,:,jv,iv)
          int(2,3,:,jv,iv)=0._r8
          int(3,3,:,jv,iv)=int(1,1,:,jv,iv)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE hyp_visc_op
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE integrands_mat
