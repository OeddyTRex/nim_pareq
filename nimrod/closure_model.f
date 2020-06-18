c-----------------------------------------------------------------------
c     file closure_model.f:  contains subprograms that perform
c     closure-related computations on physical-field data located at
c     the quadrature points for the finite elements.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  closure_model_mod
c     1.  closure_model_set
c     2.  kappa_standard_n0.
c     3.  kappa_standard_phi.
c     4.  kappa_braginskii_phi.
c     5.  kappa_k2_phi.
c     6.  calc_tdep_coul_log.
c     7.  tdep_thermal_equil.
c-----------------------------------------------------------------------
c     0. module declaration for closure_model_mod
c-----------------------------------------------------------------------
      MODULE closure_model_mod
      USE local
      IMPLICIT NONE

      REAL(r8) :: tequil_coeff

      REAL(r8), DIMENSION(2) :: chi_coeff,mag_coeff
      LOGICAL :: set_closure_coeff=.false.

      REAL(r8) :: khat_plle
      REAL(r8), DIMENSION(5) :: kprpe_coeff

      REAL(r8), PARAMETER, PRIVATE :: smalltemp=1.e-10_r8,
     $          coul_min=0.1_r8

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. closure_model_set
c-----------------------------------------------------------------------
      SUBROUTINE closure_model_set
      USE local
      USE physdat
      USE input
      IMPLICIT NONE

      set_closure_coeff=.true.

      tequil_coeff=gamm1*0.25*SQRT(2.)*SQRT(ms(1))/ms(2)
     $    /(pi**1.5*eps0**2*kboltz**1.5*zeff)
     $    *ABS(qs(1)**2*qs(2)**2)

      IF (closure_model=="k2") THEN
        mag_coeff(1)=6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**1.5
     $      /(ABS(qs(1))**3*SQRT(ms(1)))
        mag_coeff(2)=6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**1.5
     $      *zeff/(ABS(qs(2))**3*SQRT(ms(2)))

        chi_coeff(1)=gamm1*6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**2.5
     $      /(ABS(qs(1))**4*SQRT(ms(1)))
        chi_coeff(2)=gamm1*6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**2.5
     $      *zeff/(ABS(qs(2))**4*SQRT(ms(2)))

      ELSE ! braginskii
        mag_coeff(1)=6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**1.5
     $      *zeff/(ABS(qs(1)*qs(2)**2)*SQRT(ms(1)))
        mag_coeff(2)=12.*pi**1.5*eps0**2*kboltz**1.5
     $      *zeff/(ABS(qs(2))**3*SQRT(ms(2)))

        chi_coeff(1)=gamm1*6.*SQRT(2.)*pi**1.5*eps0**2*kboltz**2.5
     $      *zeff/(ABS(qs(1)**2*qs(2)**2)*SQRT(ms(1)))
        chi_coeff(2)=gamm1*12.*pi**1.5*eps0**2*kboltz**2.5
     $      *zeff/(ABS(qs(2))**4*SQRT(ms(2)))
      ENDIF

      IF (closure_model=="k2") THEN
        khat_plle = (13.5*zeff**2+54.4*zeff+25.2)/
     $      (zeff**3+8.35*zeff**2+15.2*zeff+4.51)
        kprpe_coeff(1) = (9.91*zeff**3+75.3*zeff**2
     $      +518.*zeff+333.)/1000.
        kprpe_coeff(2) = (0.211*zeff**3+12.7*zeff**2
     $      +48.4*zeff+6.45)/(zeff+57.1)
        kprpe_coeff(3) = (0.932*zeff**(7._r8/3._r8)
     $      +0.135*zeff**2+12.3*zeff+8.77)/(zeff+4.84)
        kprpe_coeff(4) = (0.246*zeff**3+2.65*zeff**2
     $      -92.8*zeff-1.96)/(zeff**2+19.9*zeff+35.3)
        kprpe_coeff(5) = (2.76*zeff**(5._r8/3._r8)
     $      -0.836*zeff**(2._r8/3._r8)-0.0611)/(zeff-0.214)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE closure_model_set
c-----------------------------------------------------------------------
c     subprogram 2. kappa_standard_n0
c     compute temperature-dependent thermal diffusivities at the
c     quadrature points as a function of toroidal grid position
c     using the standard high-magnetization limit formulation of
c     Braginskii.
c-----------------------------------------------------------------------
      SUBROUTINE kappa_standard_n0
      USE local
      USE fft_mod
      USE fields
      USE global
      USE input
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ig,lx,ly,im,ierror
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: tn0
c-----------------------------------------------------------------------
c     temperature is saved by Fourier component only.  add
c     equilibrium temperature temporarily.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        lx=SIZE(rb(ibl)%qtele%qpf,2)
        ly=SIZE(rb(ibl)%qtele%qpf,3)
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            ALLOCATE(tn0(1,lx,ly))
            tn0=rb(ibl)%qtion%qpf(:,:,:,im)
            rb(ibl)%qtion%qpf(:,:,:,im)=
     $        rb(ibl)%qtion%qpf(:,:,:,im)+rb(ibl)%qtion_eq%qpf
            IF (p_model=='aniso_tdep') THEN
              rb(ibl)%qkaprpi_n0%qpf(1,:,:)=k_perpi*
     $          SQRT(MAX(1._r8+(kprp_mnrat*
     $            (rb(ibl)%qdiff_shape%qpf(1,:,:)-1)/dvac)**2,
     $             k_pll_ref_t/REAL(rb(ibl)%qtion%qpf(1,:,:,im))))/
     $          SUM((     rb(ibl)%qbe_eq%qpf+
     $               REAL(rb(ibl)%qbe%qpf(:,:,:,im)))**2,1)
            ENDIF
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               rb(ibl)%qtion%qpf,rb(ibl)%qkappli_phi%qpf,dealiase)
        rb(ibl)%qkappli_phi%qpf=MAX(k_pll_min, MIN(k_pll_max, k_plli*
     $    (MAX(smallnum,rb(ibl)%qkappli_phi%qpf)/k_pll_ref_t)**2.5 ))
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            rb(ibl)%qtion%qpf(:,:,:,im)=tn0
            IF (.NOT.separate_pe) DEALLOCATE(tn0)
            EXIT
          ENDIF
        ENDDO
        IF (separate_pe) THEN
          DO im=1,nmodes
            IF (keff(im)==0) THEN
              tn0=rb(ibl)%qtele%qpf(:,:,:,im)
              rb(ibl)%qtele%qpf(:,:,:,im)=
     $          rb(ibl)%qtele%qpf(:,:,:,im)+rb(ibl)%qtele_eq%qpf
              IF (p_model=='aniso_tdep') THEN
                rb(ibl)%qkaprpe_n0%qpf(1,:,:)=k_perpe*
     $            SQRT(MAX(1._r8+(kprp_mnrat*
     $              (rb(ibl)%qdiff_shape%qpf(1,:,:)-1)/dvac)**2,
     $               k_pll_ref_t/REAL(rb(ibl)%qtele%qpf(1,:,:,im))))/
     $            SUM((     rb(ibl)%qbe_eq%qpf+
     $                 REAL(rb(ibl)%qbe%qpf(:,:,:,im)))**2,1)
              ENDIF
              EXIT
            ENDIF
          ENDDO
          CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $                 rb(ibl)%qtele%qpf,
     $                 rb(ibl)%qkapple_phi%qpf,dealiase)
          rb(ibl)%qkapple_phi%qpf=MAX(k_pll_min, MIN(k_pll_max, k_plle*
     $      (MAX(smallnum,rb(ibl)%qkapple_phi%qpf)/k_pll_ref_t)**2.5 ))
          DO im=1,nmodes
            IF (keff(im)==0) THEN
              rb(ibl)%qtele%qpf(:,:,:,im)=tn0
              DEALLOCATE(tn0)
              EXIT
            ENDIF
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       broadcast k_perp if needed
c-----------------------------------------------------------------------
        IF (nlayers>1.AND.p_model=='aniso_tdep') THEN
          CALL mpi_bcast(rb(ibl)%qkaprpi_n0%qpf,lx*ly,
     $                   mpi_nim_real,0,comm_mode,ierror)
          IF (separate_pe) THEN
            CALL mpi_bcast(rb(ibl)%qkaprpe_n0%qpf,lx*ly,
     $                     mpi_nim_real,0,comm_mode,ierror)
          ENDIF
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        lx=SIZE(tb(ibl)%qtele%qpf,2)
        ly=SIZE(tb(ibl)%qtele%qpf,3)
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            ALLOCATE(tn0(1,lx,ly))
            tn0=tb(ibl)%qtion%qpf(:,:,:,im)
            tb(ibl)%qtion%qpf(:,:,:,im)=
     $        tb(ibl)%qtion%qpf(:,:,:,im)+tb(ibl)%qtion_eq%qpf
            IF (p_model=='aniso_tdep') THEN
              tb(ibl)%qkaprpi_n0%qpf(1,:,:)=k_perpi*
     $          SQRT(MAX(1._r8+(kprp_mnrat*
     $            (tb(ibl)%qdiff_shape%qpf(1,:,:)-1)/dvac)**2,
     $             k_pll_ref_t/REAL(tb(ibl)%qtion%qpf(1,:,:,im))))/
     $          SUM((     tb(ibl)%qbe_eq%qpf+
     $               REAL(tb(ibl)%qbe%qpf(:,:,:,im)))**2,1)
            ENDIF
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               tb(ibl)%qtion%qpf,tb(ibl)%qkappli_phi%qpf,dealiase)
        tb(ibl)%qkappli_phi%qpf=MAX(k_pll_min, MIN(k_pll_max, k_plli*
     $    (MAX(smallnum,tb(ibl)%qkappli_phi%qpf)/k_pll_ref_t)**2.5 ))
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            tb(ibl)%qtion%qpf(:,:,:,im)=tn0
            IF (.NOT.separate_pe) DEALLOCATE(tn0)
            EXIT
          ENDIF
        ENDDO
        IF (separate_pe) THEN
          DO im=1,nmodes
            IF (keff(im)==0) THEN
              tn0=tb(ibl)%qtele%qpf(:,:,:,im)
              tb(ibl)%qtele%qpf(:,:,:,im)=
     $          tb(ibl)%qtele%qpf(:,:,:,im)+tb(ibl)%qtele_eq%qpf
              IF (p_model=='aniso_tdep') THEN
                tb(ibl)%qkaprpe_n0%qpf(1,:,:)=k_perpe*
     $            SQRT(MAX(1._r8+(kprp_mnrat*
     $              (tb(ibl)%qdiff_shape%qpf(1,:,:)-1)/dvac)**2,
     $               k_pll_ref_t/REAL(tb(ibl)%qtele%qpf(1,:,:,im))))/
     $            SUM((     tb(ibl)%qbe_eq%qpf+
     $                 REAL(tb(ibl)%qbe%qpf(:,:,:,im)))**2,1)
              ENDIF
              EXIT
            ENDIF
          ENDDO
          CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $                 tb(ibl)%qtele%qpf,
     $                 tb(ibl)%qkapple_phi%qpf,dealiase)
          tb(ibl)%qkapple_phi%qpf=MAX(k_pll_min, MIN(k_pll_max, k_plle*
     $      (MAX(smallnum,tb(ibl)%qkapple_phi%qpf)/k_pll_ref_t)**2.5 ))
          DO im=1,nmodes
            IF (keff(im)==0) THEN
              tb(ibl)%qtele%qpf(:,:,:,im)=tn0
              DEALLOCATE(tn0)
              EXIT
            ENDIF
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       broadcast k_perp if needed
c-----------------------------------------------------------------------
        IF (nlayers>1.AND.p_model=='aniso_tdep') THEN
          CALL mpi_bcast(tb(ibl)%qkaprpi_n0%qpf,lx*ly,
     $                   mpi_nim_real,0,comm_mode,ierror)
          IF (separate_pe) THEN
            CALL mpi_bcast(tb(ibl)%qkaprpe_n0%qpf,lx*ly,
     $                     mpi_nim_real,0,comm_mode,ierror)
          ENDIF
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE kappa_standard_n0
c-----------------------------------------------------------------------
c     subprogram 3. kappa_standard_phi
c     compute temperature-dependent thermal diffusivities at the
c     quadrature points as a function of toroidal grid position
c     using the standard high-magnetization limit formulation of
c     Braginskii.
c
c     this routine is only called when p_model is "aniso_tdep" and
c     closure_n0_only is then automatically false.
c-----------------------------------------------------------------------
      SUBROUTINE kappa_standard_phi
      USE local
      USE fft_mod
      USE fields
      USE global
      USE input
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      INTEGER(i4) :: ibl,ig,im,ierror,iph,ipst,ipen,npol
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: tie
      REAL(r8), DIMENSION(:), ALLOCATABLE :: qdiff_1d
c-----------------------------------------------------------------------
c     temperature is saved by Fourier component only.  add
c     equilibrium temperature temporarily.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ipst=ipqst_block(ibl)
        ipen=ipqen_block(ibl)
        npol=rb(ibl)%ng*rb(ibl)%mx*rb(ibl)%my

        ALLOCATE(tie(mpsq_block(ibl),nphi))
        ALLOCATE(qdiff_1d(npol))
        qdiff_1d=RESHAPE(rb(ibl)%qdiff_shape%qpf,(/npol/))

        CALL qp_fft_save(rb(ibl)%qtion%qpf,tie,
     $       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $       1_i4,rb(ibl)%ng,rb(ibl)%qtion_eq%qpf)
        tie=tie/k_pll_ref_t
        DO iph=1,nphi
          rb(ibl)%qkaprpi_phi%qpf(1,:,iph)=k_perpi*
     $      SQRT(MAX(1._r8+(kprp_mnrat*
     $      (qdiff_1d(ipst:ipen)-1._r8)/dvac)**2,1._r8/tie(:,iph)))/
     $      SUM((rb(ibl)%qbe_tot%qpf(:,:,iph))**2,1)
        ENDDO
        rb(ibl)%qkappli_phi%qpf(1,:,:)=MAX(k_pll_min,
     $      MIN(k_pll_max,k_plli*(MAX(smallnum,tie)**2.5)))

        rb_sep_pe: IF (separate_pe) THEN
          CALL qp_fft_save(rb(ibl)%qtele%qpf,tie,
     $         rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $         1_i4,rb(ibl)%ng,rb(ibl)%qtele_eq%qpf)
          tie=tie/k_pll_ref_t
          DO iph=1,nphi
            rb(ibl)%qkaprpe_phi%qpf(1,:,iph)=k_perpe*
     $        SQRT(MAX(1._r8+(kprp_mnrat*
     $        (qdiff_1d(ipst:ipen)-1._r8)/dvac)**2,1._r8/tie(:,iph)))/
     $        SUM((rb(ibl)%qbe_tot%qpf(:,:,iph))**2,1)
          ENDDO
          rb(ibl)%qkapple_phi%qpf(1,:,:)=MAX(k_pll_min,
     $        MIN(k_pll_max,k_plle*(MAX(smallnum,tie)**2.5)))
        ENDIF rb_sep_pe
        DEALLOCATE(tie,qdiff_1d)
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        ipst=ipqst_block(ibl)
        ipen=ipqen_block(ibl)
        npol=tb(ibl)%ng*tb(ibl)%mcell

        ALLOCATE(tie(mpsq_block(ibl),nphi))
        ALLOCATE(qdiff_1d(npol))
        qdiff_1d=RESHAPE(tb(ibl)%qdiff_shape%qpf,(/npol/))

        CALL qp_fft_save(tb(ibl)%qtion%qpf,tie,
     $       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $       1_i4,tb(ibl)%ng,tb(ibl)%qtion_eq%qpf)
        tie=tie/k_pll_ref_t
        DO iph=1,nphi
          tb(ibl)%qkaprpi_phi%qpf(1,:,iph)=k_perpi*
     $      SQRT(MAX(1._r8+(kprp_mnrat*
     $      (qdiff_1d(ipst:ipen)-1._r8)/dvac)**2,1._r8/tie(:,iph)))/
     $      SUM((tb(ibl)%qbe_tot%qpf(:,:,iph))**2,1)
        ENDDO
        tb(ibl)%qkappli_phi%qpf(1,:,:)=MAX(k_pll_min,
     $      MIN(k_pll_max,k_plli*(MAX(smallnum,tie)**2.5)))

        tb_sep_pe: IF (separate_pe) THEN
          CALL qp_fft_save(tb(ibl)%qtele%qpf,tie,
     $         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $         1_i4,tb(ibl)%ng,tb(ibl)%qtele_eq%qpf)
          tie=tie/k_pll_ref_t
          DO iph=1,nphi
            tb(ibl)%qkaprpe_phi%qpf(1,:,iph)=k_perpe*
     $        SQRT(MAX(1._r8+(kprp_mnrat*
     $        (qdiff_1d(ipst:ipen)-1._r8)/dvac)**2,1._r8/tie(:,iph)))/
     $        SUM((tb(ibl)%qbe_tot%qpf(:,:,iph))**2,1)
          ENDDO
          tb(ibl)%qkapple_phi%qpf(1,:,:)=MAX(k_pll_min,
     $        MIN(k_pll_max,k_plle*(MAX(smallnum,tie)**2.5)))
        ENDIF tb_sep_pe
        DEALLOCATE(tie,qdiff_1d)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE kappa_standard_phi
c-----------------------------------------------------------------------
c     subprogram 4. kappa_braginskii_phi
c     compute temperature-dependent thermal diffusivities at the
c     quadrature points as a function of toroidal grid position
c     using the full Braginskii formulation.
c-----------------------------------------------------------------------
      SUBROUTINE kappa_braginskii_phi
      USE local
      USE fft_mod
      USE fields
      USE global
      USE input
      USE mpi_nim
      USE pardata
      USE physdat
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ig,im,ierror
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: te,ti,xs,bmag
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: coul_vec
c-----------------------------------------------------------------------
c     set closure coefficients if not already set
c-----------------------------------------------------------------------
      IF (.NOT.set_closure_coeff) CALL closure_model_set
c-----------------------------------------------------------------------
c     temperature is saved by Fourier component only.  add
c     equilibrium temperature temporarily.
c
c     the Coulomb-log is evaluated in coul_vec, which then becomes
c     the factor with number density used in the magnetization and
c     conductivity computations.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ALLOCATE(coul_vec(2,mpsq_block(ibl),nphi),
     $           te(mpsq_block(ibl),nphi),
     $           ti(mpsq_block(ibl),nphi),
     $           bmag(mpsq_block(ibl),nphi),
     $           xs(mpsq_block(ibl),nphi))

        CALL qp_fft_save(rb(ibl)%qtion%qpf,ti,
     $       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $       1_i4,rb(ibl)%ng,rb(ibl)%qtion_eq%qpf)
        ti=MAX(ti/k_pll_ref_t,smallnum)
        IF (separate_pe) THEN
          CALL qp_fft_save(rb(ibl)%qtele%qpf,te,
     $         rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $         1_i4,rb(ibl)%ng,rb(ibl)%qtele_eq%qpf)
          te=MAX(te/k_pll_ref_t,smallnum)
        ELSE
          te=ti/zeff*pe_frac/MAX(1.-pe_frac,smallnum)
        ENDIF

        IF (tdep_coul_log) THEN
          CALL calc_tdep_coul_log(coul_vec,
     $        rb(ibl)%qnd_tot%qpf(1,:,:),ti,te)
        ELSE
          coul_vec(:,:,:)=coulomb_logarithm
        ENDIF
        coul_vec(1,:,:)=1._r8/
     $    MAX(rb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(1,:,:),smallnum)
        coul_vec(2,:,:)=1._r8/
     $    MAX(rb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(2,:,:),smallnum)

        bmag=SQRT(SUM(rb(ibl)%qbe_tot%qpf**2,1))
        xs=magfac_ion*mag_coeff(2)*ti**1.5*bmag*coul_vec(2,:,:)
        rb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $      k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $      *(gamma1_ion*xs**2+gamma0_ion)
     $      /(xs**4+delta1_ion*xs**2+delta0_ion)

        IF (.NOT.separate_pe.AND.pe_frac/=0._r8) THEN
          rb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *gamma0_ele/delta0_ele))
        ELSE
          rb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *gamma0_ion/delta0_ion))
        ENDIF

        rb_sep_pe: IF (separate_pe) THEN
          xs=magfac_ele*mag_coeff(1)*te**1.5*bmag*coul_vec(1,:,:)
          rb(ibl)%qkaprpe_phi%qpf(1,:,:)=
     $        k_perpe*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *(gamma1_ele*xs**2+gamma0_ele)
     $        /(xs**4+delta1_ele*xs**2+delta0_ele)
          rb(ibl)%qkapple_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *gamma0_ele/delta0_ele))
        ENDIF rb_sep_pe
        DEALLOCATE(coul_vec,ti,te,xs,bmag)
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        ALLOCATE(coul_vec(2,mpsq_block(ibl),nphi),
     $           te(mpsq_block(ibl),nphi),
     $           ti(mpsq_block(ibl),nphi),
     $           bmag(mpsq_block(ibl),nphi),
     $           xs(mpsq_block(ibl),nphi))

        CALL qp_fft_save(tb(ibl)%qtion%qpf,ti,
     $       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $       1_i4,tb(ibl)%ng,tb(ibl)%qtion_eq%qpf)
        ti=MAX(ti/k_pll_ref_t,smallnum)
        IF (separate_pe) THEN
          CALL qp_fft_save(tb(ibl)%qtele%qpf,te,
     $         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $         1_i4,tb(ibl)%ng,tb(ibl)%qtele_eq%qpf)
          te=MAX(te/k_pll_ref_t,smallnum)
        ELSE
          te=ti/zeff*pe_frac/MAX(1.-pe_frac,smallnum)
        ENDIF

        IF (tdep_coul_log) THEN
          CALL calc_tdep_coul_log(coul_vec,
     $        tb(ibl)%qnd_tot%qpf(1,:,:),ti,te)
        ELSE
          coul_vec(:,:,:)=coulomb_logarithm
        ENDIF
        coul_vec(1,:,:)=1._r8/
     $    MAX(tb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(1,:,:),smallnum)
        coul_vec(2,:,:)=1._r8/
     $    MAX(tb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(2,:,:),smallnum)

        bmag=SQRT(SUM(tb(ibl)%qbe_tot%qpf**2,1))
        xs=magfac_ion*chi_coeff(2)*ti**1.5*bmag*coul_vec(2,:,:)
        tb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $      k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $      *(gamma1_ion*xs**2+gamma0_ion)
     $      /(xs**4+delta1_ion*xs**2+delta0_ion)

        IF (.NOT.separate_pe.AND.pe_frac/=0._r8) THEN
          tb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *gamma0_ele/delta0_ele))
        ELSE
          tb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *gamma0_ion/delta0_ion))
        ENDIF

        tb_sep_pe: IF (separate_pe) THEN
          CALL qp_fft_save(tb(ibl)%qtele%qpf,te,
     $         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $         1_i4,tb(ibl)%ng,tb(ibl)%qtele_eq%qpf)
          te=te/k_pll_ref_t

          xs=magfac_ele*mag_coeff(1)*te**1.5*bmag*coul_vec(1,:,:)
          tb(ibl)%qkaprpe_phi%qpf(1,:,:)=
     $        k_perpe*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *(gamma1_ele*xs**2+gamma0_ele)
     $        /(xs**4+delta1_ele*xs**2+delta0_ele)
          tb(ibl)%qkapple_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *gamma0_ele/delta0_ele))
        ENDIF tb_sep_pe
        DEALLOCATE(coul_vec,ti,te,xs,bmag)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE kappa_braginskii_phi
c-----------------------------------------------------------------------
c     subprogram 5. kappa_k2_phi
c     compute temperature-dependent thermal diffusivities at the
c     quadrature points as a function of toroidal grid position
c     using the k2 closure scheme.
c-----------------------------------------------------------------------
      SUBROUTINE kappa_k2_phi
      USE local
      USE fft_mod
      USE fields
      USE global
      USE input
      USE mpi_nim
      USE pardata
      USE physdat
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ig,im,ierror

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: 
     $    delta_plli,delta_prpi,khat_plli,te,ti,xs,zeta
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: coul_vec
c-----------------------------------------------------------------------
c     set closure coefficients if not already set
c-----------------------------------------------------------------------
      IF (.NOT.set_closure_coeff) CALL closure_model_set
c-----------------------------------------------------------------------
c     temperature is saved by Fourier component only.  add
c     equilibrium temperature temporarily.
c
c     the Coulomb-log is evaluated in coul_vec, which then becomes
c     the factor with number density used in the magnetization and
c     conductivity computations.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ALLOCATE(coul_vec(2,mpsq_block(ibl),nphi),
     $           delta_plli(mpsq_block(ibl),nphi),
     $           delta_prpi(mpsq_block(ibl),nphi),
     $           khat_plli(mpsq_block(ibl),nphi),
     $           te(mpsq_block(ibl),nphi),
     $           ti(mpsq_block(ibl),nphi),
     $           xs(mpsq_block(ibl),nphi),
     $           zeta(mpsq_block(ibl),nphi))

        CALL qp_fft_save(rb(ibl)%qtion%qpf,ti,
     $      rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $      1_i4,rb(ibl)%ng,rb(ibl)%qtion_eq%qpf)
        ti=MAX(ti/k_pll_ref_t,smallnum)

        IF (separate_pe) THEN
          CALL qp_fft_save(rb(ibl)%qtele%qpf,te,
     $        rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $        1_i4,rb(ibl)%ng,rb(ibl)%qtele_eq%qpf)
          te=MAX(te/k_pll_ref_t,smallnum)
          zeta = SQRT(ms(1)/ms(2)*ti/te)/zeff
        ELSE
          zeta = SQRT(ms(1)/ms(2)*zeff*(1.-pe_frac)/
     $        MAX(pe_frac,smallnum))/zeff
          te = ti/zeff*pe_frac/MAX(1.-pe_frac,smallnum)
        ENDIF

        IF (tdep_coul_log) THEN
          CALL calc_tdep_coul_log(coul_vec,
     $        rb(ibl)%qnd_tot%qpf(1,:,:),ti,te)
        ELSE
          coul_vec(:,:,:)=coulomb_logarithm
        ENDIF
        coul_vec(1,:,:)=1._r8/
     $    MAX(rb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(1,:,:),smallnum)
        coul_vec(2,:,:)=1._r8/
     $    MAX(rb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(2,:,:),smallnum)

        xs=magfac_ion*mag_coeff(2)*ti**1.5
     $      *SQRT(SUM(rb(ibl)%qbe_tot%qpf**2,1))*coul_vec(2,:,:)

        delta_plli = 1.+13.50*zeta+36.46*zeta**2
        delta_prpi = xs**4+(1.352+12.49*zeta+34.*zeta**2)*xs**2
     $      +0.1693*delta_plli**2
        khat_plli = (5.524+30.38*zeta)/delta_plli

        IF (separate_pe) THEN
          rb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $        k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *((SQRT(2.)+7.5*zeta)*xs**2+
     $            0.1693*khat_plli*delta_plli**2)/delta_prpi
          rb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)*khat_plli))

          xs=magfac_ele*mag_coeff(1)*te**1.5
     $        *SQRT(SUM(rb(ibl)%qbe_tot%qpf**2,1))*coul_vec(1,:,:)
          rb(ibl)%qkaprpe_phi%qpf(1,:,:)=
     $        k_perpe*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *((14./3.*zeta+SQRT(2.))*xs+kprpe_coeff(1)*khat_plle)
     $        /(xs**3+kprpe_coeff(5)*xs**(7./3.)+kprpe_coeff(4)*xs**2
     $            +kprpe_coeff(3)*xs**(5./3.)
     $            +kprpe_coeff(2)*xs+kprpe_coeff(1))
          rb(ibl)%qkapple_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)*khat_plle))
        ELSE
          rb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $        k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *((SQRT(2.)+7.5*zeta)*xs**2+
     $            0.1693*khat_plli*delta_plli**2)/delta_prpi
          IF (pe_frac/=0._r8) THEN
            rb(ibl)%qkappli_phi%qpf(1,:,:)=
     $          MAX(k_pll_min,MIN(k_pll_max,
     $          k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)*khat_plle))
          ELSE
            rb(ibl)%qkappli_phi%qpf(1,:,:)=
     $          MAX(k_pll_min,MIN(k_pll_max,
     $          k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)*khat_plli))
          ENDIF
        ENDIF

        DEALLOCATE(coul_vec,delta_plli,delta_prpi,khat_plli)
        DEALLOCATE(te,ti,xs,zeta)
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        ALLOCATE(coul_vec(2,mpsq_block(ibl),nphi),
     $           delta_plli(mpsq_block(ibl),nphi),
     $           delta_prpi(mpsq_block(ibl),nphi),
     $           khat_plli(mpsq_block(ibl),nphi),
     $           te(mpsq_block(ibl),nphi),
     $           ti(mpsq_block(ibl),nphi),
     $           xs(mpsq_block(ibl),nphi),
     $           zeta(mpsq_block(ibl),nphi))

        CALL qp_fft_save(tb(ibl)%qtion%qpf,ti,
     $      tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $      1_i4,tb(ibl)%ng,tb(ibl)%qtion_eq%qpf)
        ti=MAX(ti/k_pll_ref_t,smallnum)

        IF (separate_pe) THEN
          CALL qp_fft_save(tb(ibl)%qtele%qpf,te,
     $        tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $        1_i4,tb(ibl)%ng,tb(ibl)%qtele_eq%qpf)
          te=MAX(te/k_pll_ref_t,smallnum)
          zeta = SQRT(ms(1)/ms(2)*ti/te)/zeff
        ELSE
          zeta = SQRT(ms(1)/ms(2)*zeff*(1.-pe_frac)/
     $        MAX(pe_frac,smallnum))/zeff
          te = ti/zeff*pe_frac/MAX(1.-pe_frac,smallnum)
        ENDIF

        IF (tdep_coul_log) THEN
          CALL calc_tdep_coul_log(coul_vec,
     $        tb(ibl)%qnd_tot%qpf(1,:,:),ti,te)
        ELSE
          coul_vec(:,:,:)=coulomb_logarithm
        ENDIF
        coul_vec(1,:,:)=1._r8/
     $    MAX(tb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(1,:,:),smallnum)
        coul_vec(2,:,:)=1._r8/
     $    MAX(tb(ibl)%qnd_tot%qpf(1,:,:)*coul_vec(2,:,:),smallnum)

        xs=magfac_ion*mag_coeff(2)*ti**1.5
     $      *SQRT(SUM(tb(ibl)%qbe_tot%qpf**2,1))*coul_vec(2,:,:)

        delta_plli = 1.+13.50*zeta+36.46*zeta**2
        delta_prpi = xs**4+(1.352+12.49*zeta+34.*zeta**2)*xs**2
     $      +0.1693*delta_plli**2
        khat_plli = (5.524+30.38*zeta)/delta_plli

        IF (separate_pe) THEN
          tb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $        k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *((SQRT(2.)+7.5*zeta)*xs**2+
     $            0.1693*khat_plli*delta_plli**2)/delta_prpi
          tb(ibl)%qkappli_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)*khat_plli))

          xs=magfac_ele*mag_coeff(1)*te**1.5
     $      *SQRT(SUM(tb(ibl)%qbe_tot%qpf**2,1))*coul_vec(1,:,:)
          tb(ibl)%qkaprpe_phi%qpf(1,:,:)=
     $        k_perpe*chi_coeff(1)*te**2.5*coul_vec(1,:,:)
     $        *((14./3.*zeta+SQRT(2.))*xs+kprpe_coeff(1)*khat_plle)
     $        /(xs**3+kprpe_coeff(5)*xs**(7./3.)+kprpe_coeff(4)*xs**2
     $            +kprpe_coeff(3)*xs**(5./3.)
     $            +kprpe_coeff(2)*xs+kprpe_coeff(1))
          tb(ibl)%qkapple_phi%qpf(1,:,:)=
     $        MAX(k_pll_min,MIN(k_pll_max,
     $        k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)*khat_plle))
        ELSE
          tb(ibl)%qkaprpi_phi%qpf(1,:,:)=
     $        k_perpi*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)
     $        *((SQRT(2.)+7.5*zeta)*xs**2+
     $            0.1693*khat_plli*delta_plli**2)/delta_prpi
          IF (pe_frac/=0._r8) THEN
            tb(ibl)%qkappli_phi%qpf(1,:,:)=
     $          MAX(k_pll_min,MIN(k_pll_max,
     $          k_plle*chi_coeff(1)*te**2.5*coul_vec(1,:,:)*khat_plle))
          ELSE
            tb(ibl)%qkappli_phi%qpf(1,:,:)=
     $          MAX(k_pll_min,MIN(k_pll_max,
     $          k_plli*chi_coeff(2)*ti**2.5*coul_vec(2,:,:)*khat_plli))
          ENDIF
        ENDIF

        DEALLOCATE(coul_vec,delta_plli,delta_prpi,khat_plli)
        DEALLOCATE(te,ti,xs,zeta)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE kappa_k2_phi
c-----------------------------------------------------------------------
c     subprogram 6. calc_tdep_coul_log
c     compute temperature-dependent coulomb logarithm coefficients
c     used for computing thermal diffusivities
c-----------------------------------------------------------------------
      SUBROUTINE calc_tdep_coul_log(coul_vec,nd,te,ti)
      USE local
      USE input
      USE physdat
      IMPLICIT NONE

      INTEGER(i4) :: lx,ly
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: nd,te,ti
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: coul_vec

      REAL(r8), DIMENSION(SIZE(nd,1),SIZE(nd,2)) :: debye_len
      REAL(r8) :: smalln
c-----------------------------------------------------------------------
c     computes debye length and corresponding coulomb logarithms
c-----------------------------------------------------------------------
      smalln=ndens/1.e20_r8
      debye_len=1./SQRT(MAX(nd,smalln)/(eps0*kboltz)
     $    *(qs(1)**2/te+qs(2)**2/(zeff*ti)))

      IF (closure_model=="k2") THEN
c       electron-electron coulomb logarithm
        coul_vec(1,:,:)=LOG(24.*pi*eps0/qs(1)**2*kboltz*te*debye_len)
c       ion-ion coulomb logarithm
        coul_vec(2,:,:)=LOG(24.*pi*eps0/qs(2)**2*kboltz*ti*debye_len)

      ELSE  ! braginskii
c       electron-ion coulomb logarithm
        coul_vec(1,:,:)=LOG(12.*pi*eps0/ABS(qs(1)*qs(2))*kboltz
     $      *(ms(1)*ti+ms(2)*te)/(ms(1)+ms(2))*debye_len)
c       ion-ion coulomb logarithm
        coul_vec(2,:,:)=LOG(24.*pi*eps0/qs(2)**2*kboltz*ti*debye_len)
      ENDIF
c-----------------------------------------------------------------------
c     limits on temps alone are not sufficient for robust operation.
c-----------------------------------------------------------------------
      coul_vec=MAX(coul_vec,coul_min)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE calc_tdep_coul_log
c-----------------------------------------------------------------------
c     subprogram 7. tdep_thermal_equil
c     Thermal equilibration between ions and electrons,
c     based on temperature dependent equilibration rate.
c
c     this version sacrifices some readability for parallel scaling
c     with respect to processor layers.  different fields are packed
c     into one Fourier-component array and one real-space array.
c
c     note that unlike the non-T-dependent equilibration, if the
c     equilibrium ion and electron temperatures differ, this routine
c     will tend to make (perturbed+equilibrium) temperatures the same.
c-----------------------------------------------------------------------
      SUBROUTINE tdep_thermal_equil
      USE local
      USE fields
      USE global
      USE input
      USE physdat
      USE fft_mod
      USE pardata
      IMPLICIT NONE

      REAL(r8) :: tfrac,tscale,smalln
      INTEGER(i4) :: ibl,iphi,lx,ly,ix,iy,npolf,npolr,ib,nq,ist,ien,
     $               im,rem,ipol
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rarr
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: coul_vec
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: carr
c-----------------------------------------------------------------------
c     set closure coefficients if not already set
c-----------------------------------------------------------------------
      IF (.NOT.set_closure_coeff) CALL closure_model_set
      smalln=ndens/1.e20_r8
c-----------------------------------------------------------------------
c     set the total number of different fields that need to be
c     transformed.  also, establish a temperature value for scaling.
c-----------------------------------------------------------------------
      nq=3
      tscale=beta*be0**2/(mu0*ndens*kboltz)
c-----------------------------------------------------------------------
c     determine the number of data nodes for this block for both the
c     Fourier and real-space arrays, taking layer decomposition into
c     account for the latter.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          lx=rb(ibl)%mx
          ly=rb(ibl)%my
        ELSE
          lx=tb(ibl)%mvert
          ly=1
        ENDIF

        npolf=(lx+1)*(ly+1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          npolf=npolf+(poly_degree-1)*(lx*(ly+1)+ly*(lx+1))+
     $                (poly_degree-1)**2*lx*ly
        ENDIF
        rem=MODULO(npolf,nlayers)
        npolr=npolf/nlayers
        ist=ilayer*npolr+1+MIN(rem,ilayer)
        ien=(ilayer+1)*npolr+MIN(rem,ilayer+1)
        IF (ilayer<rem) npolr=npolr+1
c-----------------------------------------------------------------------
c       in this array allocation, total density will be the first
c       quantity, total ion temperature is the second, and total
c       electron temperature is the third.
c-----------------------------------------------------------------------
        ALLOCATE(carr(nq,npolf,nmodes))
        ALLOCATE(rarr(nq,npolr,nphi))
        ALLOCATE(coul_vec(npolr,nphi))
c-----------------------------------------------------------------------
c       gather number density and the two temperatures.
c-----------------------------------------------------------------------
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            ipol=1
            DO iy=0,ly
              DO ix=0,lx
                carr(1,ipol,im)=nd(ibl)%arr(1,ix,iy,im)
     $                         +nd_eq(ibl)%arr(1,ix,iy)
                carr(2,ipol,im)=tion(ibl)%arr(1,ix,iy,im)
     $                         +tion_eq(ibl)%arr(1,ix,iy)
                carr(3,ipol,im)=tele(ibl)%arr(1,ix,iy,im)
     $                         +tele_eq(ibl)%arr(1,ix,iy)
                ipol=ipol+1
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=0,ly
                DO ix=1,lx
                  carr(1,ipol,im)=nd(ibl)%arrh(1,ib,ix,iy,im)
     $                           +nd_eq(ibl)%arrh(1,ib,ix,iy)
                  carr(2,ipol,im)=tion(ibl)%arrh(1,ib,ix,iy,im)
     $                           +tion_eq(ibl)%arrh(1,ib,ix,iy)
                  carr(3,ipol,im)=tele(ibl)%arrh(1,ib,ix,iy,im)
     $                           +tele_eq(ibl)%arrh(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=1,ly
                DO ix=0,lx
                  carr(1,ipol,im)=nd(ibl)%arrv(1,ib,ix,iy,im)
     $                           +nd_eq(ibl)%arrv(1,ib,ix,iy)
                  carr(2,ipol,im)=tion(ibl)%arrv(1,ib,ix,iy,im)
     $                           +tion_eq(ibl)%arrv(1,ib,ix,iy)
                  carr(3,ipol,im)=tele(ibl)%arrv(1,ib,ix,iy,im)
     $                           +tele_eq(ibl)%arrv(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,(poly_degree-1)**2
              DO iy=1,ly
                DO ix=1,lx
                  carr(1,ipol,im)=nd(ibl)%arri(1,ib,ix,iy,im)
     $                           +nd_eq(ibl)%arri(1,ib,ix,iy)
                  carr(2,ipol,im)=tion(ibl)%arri(1,ib,ix,iy,im)
     $                           +tion_eq(ibl)%arri(1,ib,ix,iy)
                  carr(3,ipol,im)=tele(ibl)%arri(1,ib,ix,iy,im)
     $                           +tele_eq(ibl)%arri(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ELSE
            ipol=1
            DO iy=0,ly
              DO ix=0,lx
                carr(1,ipol,im)=nd(ibl)%arr(1,ix,iy,im)
                carr(2,ipol,im)=tion(ibl)%arr(1,ix,iy,im)
                carr(3,ipol,im)=tele(ibl)%arr(1,ix,iy,im)
                ipol=ipol+1
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=0,ly
                DO ix=1,lx
                  carr(1,ipol,im)=nd(ibl)%arrh(1,ib,ix,iy,im)
                  carr(2,ipol,im)=tion(ibl)%arrh(1,ib,ix,iy,im)
                  carr(3,ipol,im)=tele(ibl)%arrh(1,ib,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=1,ly
                DO ix=0,lx
                  carr(1,ipol,im)=nd(ibl)%arrv(1,ib,ix,iy,im)
                  carr(2,ipol,im)=tion(ibl)%arrv(1,ib,ix,iy,im)
                  carr(3,ipol,im)=tele(ibl)%arrv(1,ib,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,(poly_degree-1)**2
              DO iy=1,ly
                DO ix=1,lx
                  carr(1,ipol,im)=nd(ibl)%arri(1,ib,ix,iy,im)
                  carr(2,ipol,im)=tion(ibl)%arri(1,ib,ix,iy,im)
                  carr(3,ipol,im)=tele(ibl)%arri(1,ib,ix,iy,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       packing in the fft routine will mix number density and
c       temperature, so scaling prior to transformation is important.
c-----------------------------------------------------------------------
        carr(1,:,:)=carr(1,:,:)/ndens
        carr(2:3,:,:)=carr(2:3,:,:)/tscale
c-----------------------------------------------------------------------
c       transform all fields simultaneously and return to physical
c       units.
c-----------------------------------------------------------------------
        CALL fft_nim('inverse',npolf,npolr,lphi,nq,carr,rarr,dealiase)

        rarr(1,:,:)=ndens*rarr(1,:,:)
        rarr(2:3,:,:)=tscale*rarr(2:3,:,:)
c-----------------------------------------------------------------------
c       calculate the Coulomb logarithm with a reasonable lower limit.
c       the first part of the temperature-dependent computation is the
c       Debye length, and the second part finds the Coulomb log.
c
c       l_Debye=1/SQRT(nd/(eps0*kboltz)*(qe**2/Te+qi**2/(Z*Ti)))
c-----------------------------------------------------------------------
        IF (tdep_coul_log) THEN
          coul_vec=1./SQRT(MAX(rarr(1,:,:),smalln)/(eps0*kboltz)
     $        *(qs(1)**2/MAX(rarr(3,:,:),smalltemp)+
     $          qs(2)**2/(zeff*MAX(rarr(2,:,:),smalltemp))))

          coul_vec=MAX(LOG(12.*pi*eps0/ABS(qs(1)*qs(2))*kboltz
     $        *(ms(1)*MAX(rarr(2,:,:),smalltemp)+
     $          ms(2)*MAX(rarr(3,:,:),smalltemp))
     $        /(ms(1)+ms(2))*coul_vec),coul_min)
        ELSE
          coul_vec=coulomb_logarithm
        ENDIF
c-----------------------------------------------------------------------
c       apply thermal equilibration rate factors.  the Coulomb log
c       data is over-written by the spatially dependent equilibration
c       rates.
c-----------------------------------------------------------------------
        tfrac=dt*tequil_rate*tequil_coeff

        coul_vec=tfrac*rarr(1,:,:)*coul_vec
     $      /MAX(rarr(3,:,:),smalltemp)**1.5
c-----------------------------------------------------------------------
c       calculates updated ion and electron temperatures
c-----------------------------------------------------------------------
        rarr(2,:,:)=(rarr(2,:,:)*(1._r8+coul_vec)
     $    +rarr(3,:,:)*zeff*coul_vec)/(1._r8+coul_vec*(zeff+1._r8))

        rarr(3,:,:)=(rarr(3,:,:)+rarr(2,:,:)*coul_vec)/(1._r8+coul_vec)
c-----------------------------------------------------------------------
c       transform back to the Fourier representation.  the number-
c       density entries are not used, so just set them to 0 to avoid
c       rescaling.
c-----------------------------------------------------------------------
        rarr(1,:,:)=0
        CALL fft_nim('forward',npolf,npolr,lphi,nq,carr,rarr,dealiase)
c-----------------------------------------------------------------------
c       put the updated temperatures back into their nodal arrays.
c-----------------------------------------------------------------------
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            ipol=1
            DO iy=0,ly
              DO ix=0,lx
                tion(ibl)%arr(1,ix,iy,im)=carr(2,ipol,im)
     $                         -tion_eq(ibl)%arr(1,ix,iy)
                tele(ibl)%arr(1,ix,iy,im)=carr(3,ipol,im)
     $                         -tele_eq(ibl)%arr(1,ix,iy)
                ipol=ipol+1
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=0,ly
                DO ix=1,lx
                  tion(ibl)%arrh(1,ib,ix,iy,im)=carr(2,ipol,im)
     $                           -tion_eq(ibl)%arrh(1,ib,ix,iy)
                  tele(ibl)%arrh(1,ib,ix,iy,im)=carr(3,ipol,im)
     $                           -tele_eq(ibl)%arrh(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=1,ly
                DO ix=0,lx
                  tion(ibl)%arrv(1,ib,ix,iy,im)=carr(2,ipol,im)
     $                           -tion_eq(ibl)%arrv(1,ib,ix,iy)
                  tele(ibl)%arrv(1,ib,ix,iy,im)=carr(3,ipol,im)
     $                           -tele_eq(ibl)%arrv(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,(poly_degree-1)**2
              DO iy=1,ly
                DO ix=1,lx
                  tion(ibl)%arri(1,ib,ix,iy,im)=carr(2,ipol,im)
     $                           -tion_eq(ibl)%arri(1,ib,ix,iy)
                  tele(ibl)%arri(1,ib,ix,iy,im)=carr(3,ipol,im)
     $                           -tele_eq(ibl)%arri(1,ib,ix,iy)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ELSE
            ipol=1
            DO iy=0,ly
              DO ix=0,lx
                tion(ibl)%arr(1,ix,iy,im)=carr(2,ipol,im)
                tele(ibl)%arr(1,ix,iy,im)=carr(3,ipol,im)
                ipol=ipol+1
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=0,ly
                DO ix=1,lx
                  tion(ibl)%arrh(1,ib,ix,iy,im)=carr(2,ipol,im)
                  tele(ibl)%arrh(1,ib,ix,iy,im)=carr(3,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,poly_degree-1
              DO iy=1,ly
                DO ix=0,lx
                  tion(ibl)%arrv(1,ib,ix,iy,im)=carr(2,ipol,im)
                  tele(ibl)%arrv(1,ib,ix,iy,im)=carr(3,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
            DO ib=1,(poly_degree-1)**2
              DO iy=1,ly
                DO ix=1,lx
                  tion(ibl)%arri(1,ib,ix,iy,im)=carr(2,ipol,im)
                  tele(ibl)%arri(1,ib,ix,iy,im)=carr(3,ipol,im)
                  ipol=ipol+1
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       deallocates quantities
c-----------------------------------------------------------------------
        DEALLOCATE(carr,rarr,coul_vec)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE tdep_thermal_equil
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE closure_model_mod
