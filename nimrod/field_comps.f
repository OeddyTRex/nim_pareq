c-----------------------------------------------------------------------
c     file field_comps.f:  contains subprograms that perform
c     computations directly on physical-field data.  these routines
c     were formerly part of utilities.f.
c
c     note that the reordering of quadrature-point array indices to be
c     before a combined element index effects many computations in this
c     file.  in some instances, array syntax makes the changes
c     transparent, and in most others, it simplifies the coding.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  surface_exb.
c     2.  p_from_nt.
c     3.  ave_field_check.
c     4.  ave_n_check.
c     5.  ave_v_check.
c     6.  ave_eta_check.
c     7.  find_eta_t.
c     8.  find_bb.
c     9.  ave_kappa_check.
c     10. find_kappa_t.
c     11. n_store.
c     12. temp_store.
c     13. vcom_store.
c     14. find_vv.
c     15. b_store.
c     17. pieq_comp.
c     18. iso_stress.
c     19. par_stress_lin.
c     20. par_stress_veq.
c     21. par_stress_nl.
c     22. viscous_heating.
c     23. thermal_equil.
c     24. si_store.
c     25. hyper_bheat.
c     26. ndiff_store.
c-----------------------------------------------------------------------
c     subprogram 1. surface_exb
c     set the n=0 flow at the surface according to the applied e and
c     local b.
c-----------------------------------------------------------------------
      SUBROUTINE surface_exb(vdum,vcent)
      USE local
      USE fields
      USE seam_storage_mod
      USE global
      USE input
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vdum
      REAL(r8), INTENT(IN) :: vcent

      REAL(r8), DIMENSION(2) :: normal,tangent,exb
      REAL(r8), DIMENSION(3) :: bes
      REAL(r8) :: ncomp,tcomp,bigr=1,dx,dy,gfac,lvolt,evert
      INTEGER(i4) :: ibe,ibl,imode,iv,ix,iy,ivp,ibase
c-----------------------------------------------------------------------
c     find the n=0 mode if there is one.
c-----------------------------------------------------------------------
      imode=1
      DO
        IF (imode>nmodes) RETURN
        IF (keff(imode)==0) EXIT
        imode=imode+1
      ENDDO
      IF (geom=='tor') THEN
        gfac=1/twopi
      ELSE
        gfac=1/per_length
      ENDIF
c-----------------------------------------------------------------------
c     center voltages according to vcent.
c-----------------------------------------------------------------------
      lvolt=vcent*volt+(1-vcent)*volt_old
      evert=vcent*e_vert+(1-vcent)*e_vert_old
c-----------------------------------------------------------------------
c     loop over all blocks touching the surface of the domain.
c-----------------------------------------------------------------------
      block: DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
c-----------------------------------------------------------------------
c       loop over the block boundary, and find exb for points on the
c       external surface (the third component of be_eq includes a 
c       factor of R).
c-----------------------------------------------------------------------
        vert: DO iv=1,seam(ibe)%nvert
          IF (.NOT.seam(ibe)%expoint(iv).OR.
     $             seam(ibe)%r0point(iv).OR.
     $             seam(ibe)%excorner(iv)) CYCLE
          ix=seam(ibe)%vertex(iv)%intxy(1)
          iy=seam(ibe)%vertex(iv)%intxy(2)
          normal =seam(ibe)%vertex(iv)%norm
          tangent=seam(ibe)%vertex(iv)%tang
          bes=be_eq(ibe)%arr(:,ix,iy)
          IF (geom=='tor') THEN
            IF (ibe<=nrbl) THEN
              bigr=rb(ibe)%rz%fs(1,ix,iy)
              bes(3)=bes(3)/bigr
            ELSE
              bigr=tb(ibe)%tgeom%xs(ix)
              bes(3)=bes(3)/bigr
            ENDIF
          ENDIF
          bes=bes+be(ibe)%arr(:,ix,iy,imode)
          exb(1)=evert*ABS(tangent(2))*bes(3)
     $            -lvolt*gfac*bes(2)/bigr
          exb(2)=lvolt*gfac*bes(1)/bigr
          ncomp=SUM(exb*normal)/SUM(bes**2)
          tcomp=SUM(vdum(ibe)%arr(1:2,ix,iy,imode)*tangent)
          vdum(ibe)%arr(1:2,ix,iy,imode)=ncomp*normal+tcomp*tangent
        ENDDO vert
c-----------------------------------------------------------------------
c       side-centered bases:
c-----------------------------------------------------------------------
        IF (poly_degree>1) THEN
          seg: DO iv=1,seam(ibe)%nvert
            ivp=iv-1
            IF (ivp==0) ivp=seam(ibe)%nvert
            IF (.NOT.(seam(ibe)%expoint(iv).AND.seam(ibe)%expoint(ivp))
     $          .OR.(seam(ibe)%r0point(iv).AND.seam(ibe)%r0point(ivp)))
     $        CYCLE
            ix=seam(ibe)%segment(iv)%intxys(1)
            iy=seam(ibe)%segment(iv)%intxys(2)
            DO ibase=1,poly_degree-1
              normal =seam(ibe)%segment(iv)%norm(:,ibase)
              tangent=seam(ibe)%segment(iv)%tang(:,ibase)
c-PRE triangles
              IF (seam(ibe)%segment(iv)%h_side) THEN
                dx=rb(ibe)%be%dx(ibase+1)+
     $             MIN(seam(ibe)%segment(iv)%intxyn(1),
     $                 seam(ibe)%segment(iv)%intxyp(1))
                dy=iy
              ELSE
                dx=ix
                dy=rb(ibe)%be%dy(ibase+poly_degree)+
     $             MIN(seam(ibe)%segment(iv)%intxyn(2),
     $                 seam(ibe)%segment(iv)%intxyp(2))
              ENDIF
              CALL lagr_quad_eval(rb(ibe)%be_eq,dx,dy,0_i4)
              bes=rb(ibe)%be_eq%f
              IF (geom=='tor') THEN
                CALL lagr_quad_eval(rb(ibe)%rz,dx,dy,0_i4)
                bigr=rb(ibe)%rz%f(1)
                bes(3)=bes(3)/bigr
              ENDIF
              IF (seam(ibe)%segment(iv)%h_side) THEN
                bes=bes+be(ibe)%arrh(:,ibase,ix,iy,imode)
                exb(1)=evert*ABS(tangent(2))*bes(3)
     $                  -lvolt*gfac*bes(2)/bigr
                exb(2)=lvolt*gfac*bes(1)/bigr
                ncomp=SUM(exb*normal)/SUM(bes**2)
                tcomp=SUM(vdum(ibe)%arrh(1:2,ibase,ix,iy,imode)*tangent)
                vdum(ibe)%arrh(1:2,ibase,ix,iy,imode)=ncomp*normal
     $                                               +tcomp*tangent
              ELSE
                bes=bes+be(ibe)%arrv(:,ibase,ix,iy,imode)
                exb(1)=evert*ABS(tangent(2))*bes(3)
     $                  -lvolt*gfac*bes(2)/bigr
                exb(2)=lvolt*gfac*bes(1)/bigr
                ncomp=SUM(exb*normal)/SUM(bes**2)
                tcomp=SUM(vdum(ibe)%arrv(1:2,ibase,ix,iy,imode)*tangent)
                vdum(ibe)%arrv(1:2,ibase,ix,iy,imode)=ncomp*normal
     $                                               +tcomp*tangent
              ENDIF
            ENDDO
          ENDDO seg
        ENDIF
      ENDDO block
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_exb
c-----------------------------------------------------------------------
c     subprogram 2. p_from_nt
c     compute the n*k*T product at quadrature points for the symmetric
c     storage or at node locations for all basis functions.
c-----------------------------------------------------------------------
      SUBROUTINE p_from_nt(flag)
      USE local
      USE pardata
      USE mpi_nim
      USE fields
      USE global
      USE input
      USE fft_mod
      USE physdat
      USE rblock
      USE tblock
      USE time
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: ibl,imode,ierror,iq,m1,m2,npol,rep,ipolst,ipolen,ns
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_n,real_p
      REAL(r8) :: timest_qp,timend_qp
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: comp_p
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: comp_p2
c-----------------------------------------------------------------------
c     multiply symmetric parts of n and T to get an approximate n=0 
c     pressure at the quadrature points.  this is only used to compute
c     the 'nonlinear pressure' coefficient for the semi-implicit
c     operator.  note that ti_n0 and te_n0 already have the equilibrium
c     part added.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      IF (flag=='qp sym') THEN
        DO ibl=1,nrbl
          rb(ibl)%qpres_n0%qpf=kboltz*(
     $      rb(ibl)%qnd_n0%qpf*(rb(ibl)%qte_n0%qpf+
     $                          rb(ibl)%qti_n0%qpf/zeff)+
     $      rb(ibl)%qnd_eq%qpf*
     $        ( rb(ibl)%qte_n0%qpf-rb(ibl)%qtele_eq%qpf+
     $         (rb(ibl)%qti_n0%qpf-rb(ibl)%qtion_eq%qpf)/zeff) )
          rb(ibl)%qpres_n0%qpfr=kboltz*(
     $      rb(ibl)%qnd_n0%qpfr*(rb(ibl)%qte_n0%qpf +
     $                           rb(ibl)%qti_n0%qpf /zeff)+
     $      rb(ibl)%qnd_n0%qpf *(rb(ibl)%qte_n0%qpfr+
     $                           rb(ibl)%qti_n0%qpfr/zeff)+
     $      rb(ibl)%qnd_eq%qpfr*
     $        ( rb(ibl)%qte_n0%qpf -rb(ibl)%qtele_eq%qpf +
     $         (rb(ibl)%qti_n0%qpf -rb(ibl)%qtion_eq%qpf )/zeff)+
     $      rb(ibl)%qnd_eq%qpf *
     $        ( rb(ibl)%qte_n0%qpfr-rb(ibl)%qtele_eq%qpfr+
     $         (rb(ibl)%qti_n0%qpfr-rb(ibl)%qtion_eq%qpfr)/zeff) )
          rb(ibl)%qpres_n0%qpfz=kboltz*(
     $      rb(ibl)%qnd_n0%qpfz*(rb(ibl)%qte_n0%qpf +
     $                           rb(ibl)%qti_n0%qpf /zeff)+
     $      rb(ibl)%qnd_n0%qpf *(rb(ibl)%qte_n0%qpfz+
     $                           rb(ibl)%qti_n0%qpfz/zeff)+
     $      rb(ibl)%qnd_eq%qpfz*
     $        ( rb(ibl)%qte_n0%qpf -rb(ibl)%qtele_eq%qpf +
     $         (rb(ibl)%qti_n0%qpf -rb(ibl)%qtion_eq%qpf )/zeff)+
     $      rb(ibl)%qnd_eq%qpf *
     $        ( rb(ibl)%qte_n0%qpfz-rb(ibl)%qtele_eq%qpfz+
     $         (rb(ibl)%qti_n0%qpfz-rb(ibl)%qtion_eq%qpfz)/zeff) )
        ENDDO
        DO ibl=nrbl+1,nbl
          tb(ibl)%qpres_n0%qpf=kboltz*(
     $      tb(ibl)%qnd_n0%qpf*(tb(ibl)%qte_n0%qpf+
     $                          tb(ibl)%qti_n0%qpf/zeff)+
     $      tb(ibl)%qnd_eq%qpf*
     $        ( tb(ibl)%qte_n0%qpf-tb(ibl)%qtele_eq%qpf+
     $         (tb(ibl)%qti_n0%qpf-tb(ibl)%qtion_eq%qpf)/zeff) )
          tb(ibl)%qpres_n0%qpfr=kboltz*(
     $      tb(ibl)%qnd_n0%qpfr*(tb(ibl)%qte_n0%qpf +
     $                           tb(ibl)%qti_n0%qpf /zeff)+
     $      tb(ibl)%qnd_n0%qpf *(tb(ibl)%qte_n0%qpfr+
     $                           tb(ibl)%qti_n0%qpfr/zeff)+
     $      tb(ibl)%qnd_eq%qpfr*
     $        ( tb(ibl)%qte_n0%qpf -tb(ibl)%qtele_eq%qpf +
     $         (tb(ibl)%qti_n0%qpf -tb(ibl)%qtion_eq%qpf )/zeff)+
     $      tb(ibl)%qnd_eq%qpf *
     $        ( tb(ibl)%qte_n0%qpfr-tb(ibl)%qtele_eq%qpfr+
     $         (tb(ibl)%qti_n0%qpfr-tb(ibl)%qtion_eq%qpfr)/zeff) )
          tb(ibl)%qpres_n0%qpfz=kboltz*(
     $      tb(ibl)%qnd_n0%qpfz*(tb(ibl)%qte_n0%qpf +
     $                           tb(ibl)%qti_n0%qpf /zeff)+
     $      tb(ibl)%qnd_n0%qpf *(tb(ibl)%qte_n0%qpfz+
     $                           tb(ibl)%qti_n0%qpfz/zeff)+
     $      tb(ibl)%qnd_eq%qpfz*
     $        ( tb(ibl)%qte_n0%qpf -tb(ibl)%qtele_eq%qpf +
     $         (tb(ibl)%qti_n0%qpf -tb(ibl)%qtion_eq%qpf )/zeff)+
     $      tb(ibl)%qnd_eq%qpf *
     $        ( tb(ibl)%qte_n0%qpfz-tb(ibl)%qtele_eq%qpfz+
     $         (tb(ibl)%qti_n0%qpfz-tb(ibl)%qtion_eq%qpfz)/zeff) )
        ENDDO
        CALL timer(timend_qp)
        time_qpfld=time_qpfld+timend_qp-timest_qp
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     multiply n and T to get electron and total pressures with
c     nonlinear and linear parts separated.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (nonlinear) THEN
          m1=SIZE(tele(ibl)%arr,2)
          m2=SIZE(tele(ibl)%arr,3)
          npol=(m1*m2)/nlayers
          rep=MODULO(m1*m2,nlayers)
          ipolst=ilayer*npol+1+MIN(rep,ilayer)
          ipolen=(ilayer+1)*npol+MIN(rep,ilayer+1)
          npol=ipolen-ipolst+1
          ALLOCATE(real_n(1,npol,nphi),real_p(2,npol,nphi))
          ALLOCATE(comp_p(2,m1,m2,nmodes))
          CALL fft_nim('inverse',m1*m2,npol,lphi,1_i4,nd(ibl)%arr,
     $                 real_n,dealiase)
          comp_p(1,:,:,:)=tele(ibl)%arr(1,:,:,:)
          comp_p(2,:,:,:)=tion(ibl)%arr(1,:,:,:)
          CALL fft_nim('inverse',m1*m2,npol,lphi,2_i4,comp_p,
     $                 real_p,dealiase)
          real_p(1,:,:)=kboltz*real_n(1,:,:)*real_p(1,:,:)
          real_p(2,:,:)=kboltz*real_n(1,:,:)*real_p(2,:,:)/zeff
          CALL fft_nim('forward',m1*m2,npol,lphi,2_i4,comp_p,
     $                 real_p,dealiase)
          prese(ibl)%arr(1,:,:,:)=comp_p(1,:,:,:)
          pres (ibl)%arr(1,:,:,:)=comp_p(2,:,:,:)
          DEALLOCATE(real_n,real_p,comp_p)

          IF (poly_degree>1.AND.ibl<=nrbl) THEN
            ns=SIZE(tele(ibl)%arrh,2)
            m1=m1-1
            npol=(ns*m1*m2)/nlayers
            rep=MODULO(ns*m1*m2,nlayers)
            ipolst=ilayer*npol+1+MIN(rep,ilayer)
            ipolen=(ilayer+1)*npol+MIN(rep,ilayer+1)
            npol=ipolen-ipolst+1
            ALLOCATE(real_n(1,npol,nphi),real_p(2,npol,nphi))
            ALLOCATE(comp_p2(2,ns,m1,m2,nmodes))
            comp_p2(1,:,:,:,:)=tele(ibl)%arrh(1,:,:,:,:)
            comp_p2(2,:,:,:,:)=tion(ibl)%arrh(1,:,:,:,:)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,1_i4,
     $                   nd(ibl)%arrh,real_n,dealiase)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            real_p(1,:,:)=kboltz*real_n(1,:,:)*real_p(1,:,:)
            real_p(2,:,:)=kboltz*real_n(1,:,:)*real_p(2,:,:)/zeff
            CALL fft_nim('forward',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            prese(ibl)%arrh(1,:,:,:,:)=comp_p2(1,:,:,:,:)
            pres (ibl)%arrh(1,:,:,:,:)=comp_p2(2,:,:,:,:)
            DEALLOCATE(real_n,real_p,comp_p2)

            ns=SIZE(tele(ibl)%arrv,2)
            m1=m1+1
            m2=m2-1
            npol=(ns*m1*m2)/nlayers
            rep=MODULO(ns*m1*m2,nlayers)
            ipolst=ilayer*npol+1+MIN(rep,ilayer)
            ipolen=(ilayer+1)*npol+MIN(rep,ilayer+1)
            npol=ipolen-ipolst+1
            ALLOCATE(real_n(1,npol,nphi),real_p(2,npol,nphi))
            ALLOCATE(comp_p2(2,ns,m1,m2,nmodes))
            comp_p2(1,:,:,:,:)=tele(ibl)%arrv(1,:,:,:,:)
            comp_p2(2,:,:,:,:)=tion(ibl)%arrv(1,:,:,:,:)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,1_i4,
     $                   nd(ibl)%arrv,real_n,dealiase)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            real_p(1,:,:)=kboltz*real_n(1,:,:)*real_p(1,:,:)
            real_p(2,:,:)=kboltz*real_n(1,:,:)*real_p(2,:,:)/zeff
            CALL fft_nim('forward',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            prese(ibl)%arrv(1,:,:,:,:)=comp_p2(1,:,:,:,:)
            pres (ibl)%arrv(1,:,:,:,:)=comp_p2(2,:,:,:,:)
            DEALLOCATE(real_n,real_p,comp_p2)

            ns=SIZE(tele(ibl)%arri,2)
            m1=m1-1
            npol=(ns*m1*m2)/nlayers
            rep=MODULO(ns*m1*m2,nlayers)
            ipolst=ilayer*npol+1+MIN(rep,ilayer)
            ipolen=(ilayer+1)*npol+MIN(rep,ilayer+1)
            npol=ipolen-ipolst+1
            ALLOCATE(real_n(1,npol,nphi),real_p(2,npol,nphi))
            ALLOCATE(comp_p2(2,ns,m1,m2,nmodes))
            comp_p2(1,:,:,:,:)=tele(ibl)%arri(1,:,:,:,:)
            comp_p2(2,:,:,:,:)=tion(ibl)%arri(1,:,:,:,:)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,1_i4,
     $                   nd(ibl)%arri,real_n,dealiase)
            CALL fft_nim('inverse',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            real_p(1,:,:)=kboltz*real_n(1,:,:)*real_p(1,:,:)
            real_p(2,:,:)=kboltz*real_n(1,:,:)*real_p(2,:,:)/zeff
            CALL fft_nim('forward',ns*m1*m2,npol,lphi,2_i4,
     $                   comp_p2,real_p,dealiase)
            prese(ibl)%arri(1,:,:,:,:)=comp_p2(1,:,:,:,:)
            pres (ibl)%arri(1,:,:,:,:)=comp_p2(2,:,:,:,:)
            DEALLOCATE(real_n,real_p,comp_p2)
          ENDIF
        ELSE
          prese(ibl)=0
          pres(ibl)=0
        ENDIF
        DO imode=1,nmodes
          prese(ibl)%arr(:,:,:,imode)=prese(ibl)%arr(:,:,:,imode)+
     $      kboltz*(nd_eq(ibl)%arr*tele(ibl)%arr(:,:,:,imode)+
     $               tele_eq(ibl)%arr*nd(ibl)%arr(:,:,:,imode))
          pres(ibl)%arr(:,:,:,imode)=
     $      pres(ibl)%arr(:,:,:,imode)+prese(ibl)%arr(:,:,:,imode)+
     $      kboltz*(nd_eq(ibl)%arr*tion(ibl)%arr(:,:,:,imode)+
     $               tion_eq(ibl)%arr*nd(ibl)%arr(:,:,:,imode))/zeff
          IF (poly_degree>1.AND.ibl<=nrbl) THEN
            prese(ibl)%arrh(:,:,:,:,imode)=
     $        prese(ibl)%arrh(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arrh*tele(ibl)%arrh(:,:,:,:,imode)+
     $                 tele_eq(ibl)%arrh*nd(ibl)%arrh(:,:,:,:,imode))
            pres(ibl)%arrh(:,:,:,:,imode)=
     $        pres(ibl)%arrh(:,:,:,:,imode)+
     $        prese(ibl)%arrh(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arrh*tion(ibl)%arrh(:,:,:,:,imode)+
     $               tion_eq(ibl)%arrh*nd(ibl)%arrh(:,:,:,:,imode))/zeff
            prese(ibl)%arrv(:,:,:,:,imode)=
     $        prese(ibl)%arrv(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arrv*tele(ibl)%arrv(:,:,:,:,imode)+
     $                 tele_eq(ibl)%arrv*nd(ibl)%arrv(:,:,:,:,imode))
            pres(ibl)%arrv(:,:,:,:,imode)=
     $        pres(ibl)%arrv(:,:,:,:,imode)+
     $        prese(ibl)%arrv(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arrv*tion(ibl)%arrv(:,:,:,:,imode)+
     $               tion_eq(ibl)%arrv*nd(ibl)%arrv(:,:,:,:,imode))/zeff
            prese(ibl)%arri(:,:,:,:,imode)=
     $        prese(ibl)%arri(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arri*tele(ibl)%arri(:,:,:,:,imode)+
     $                 tele_eq(ibl)%arri*nd(ibl)%arri(:,:,:,:,imode))
            pres(ibl)%arri(:,:,:,:,imode)=
     $        pres(ibl)%arri(:,:,:,:,imode)+
     $        prese(ibl)%arri(:,:,:,:,imode)+
     $        kboltz*(nd_eq(ibl)%arri*tion(ibl)%arri(:,:,:,:,imode)+
     $               tion_eq(ibl)%arri*nd(ibl)%arri(:,:,:,:,imode))/zeff
          ENDIF
        ENDDO
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE p_from_nt
c-----------------------------------------------------------------------
c     subprogram 3. ave_field_check
c     find the fractional change in average B and P to determine if
c     the linear part of the semi-implicit operators need updating in
c     nonlinear runs.  also check the change in maximum pressures for
c     updating the nonlinear part of the operator.
c-----------------------------------------------------------------------
      SUBROUTINE ave_field_check
      USE local
      USE fields
      USE input
      USE global
      USE physdat
      USE fft_mod
      USE mpi_nim
      USE pardata
      USE rblock
      USE tblock
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,m1,m2,ierror,imode,ip,ig,ng,mxb,myb,ii,ipg
      INTEGER(i4) :: nl,npol,ipolst,ipolen,rep,ix,iy,ipol,ipolm,il
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls
      REAL(r8) :: tmp,tm2,ave_b_change,ave_p_change,nl_change
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: aveb,teff
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: nltmp2
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          btot,ntot,tiq,teq,
     $          ben0,beq,ten0,tin0,ndn0,prq,prn0,nl_pr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: cteff
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: cti,cte
      LOGICAL, SAVE :: first_call=.true.
c-----------------------------------------------------------------------
c     check the n=0 fields, which are in layer 0 if this is a parallel
c     run with multiple Fourier layers.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      ave_b_change=0
      ave_p_change=0
      nl_change=0
c-----------------------------------------------------------------------
c     find the relative change for the n=0 + equilibrium fields.
c     just use grid vertex-centered information for this check.
c-----------------------------------------------------------------------
      IF (ilayer==0.AND..NOT.first_call) THEN
        DO ibl=1,nbl
          IF (geom=='tor') THEN
            m1=SIZE(be(ibl)%arr,2)
            m2=SIZE(be(ibl)%arr,3)
            ALLOCATE(aveb(3,m1,m2))
            aveb=be_eq(ibl)%arr
            IF (ibl<=nrbl) THEN
              aveb(3,:,:)=aveb(3,:,:)/(rb(ibl)%rz%fs(1,:,:)+smallnum)
            ELSE
              aveb(3,:,1)=aveb(3,:,1)/(tb(ibl)%tgeom%xs+smallnum)
            ENDIF
            aveb=aveb+be_n0(ibl)%arr
            ave_b_change=MAX(ave_b_change,MAXVAL(SUM((be_n0(ibl)%arr
     $        -REAL(be(ibl)%arr(:,:,:,1),r8))**2,1)/
     $               MAX(SUM(aveb**2,1),smallnum) ) )
            DEALLOCATE(aveb)
          ELSE
            ave_b_change=MAX(ave_b_change,MAXVAL(SUM((be_n0(ibl)%arr
     $         -REAL(be(ibl)%arr(:,:,:,1),r8))**2,1)
     $         /MAX(SUM((be_n0(ibl)%arr+be_eq(ibl)%arr)**2,1),
     $              smallnum)))
          ENDIF
c-----------------------------------------------------------------------
c         the pressure coefficients are not used in the time-stepping, 
c         so compute an approximate of n=0 p from the n=0 n and T.
c-----------------------------------------------------------------------
          IF (beta>0) THEN
            work3(ibl)%arr(:,:,:,1)=kboltz*(
     $        (nd_eq(ibl)%arr+nd(ibl)%arr(:,:,:,1))*
     $          (tele(ibl)%arr(:,:,:,1)+tion(ibl)%arr(:,:,:,1)/zeff)+
     $        nd(ibl)%arr(:,:,:,1)*
     $          (tele_eq(ibl)%arr+tion_eq(ibl)%arr/zeff) )
            ave_p_change=MAX(ave_p_change,MAXVAL(ABS(pres_n0(ibl)%arr
     $        -REAL(work3(ibl)%arr(:,:,:,1),r8))
     $        /ABS(pres_n0(ibl)%arr+pres_eq(ibl)%arr+smallnum)))
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate max change to all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(ave_b_change,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        ave_b_change=tmp
        IF (beta>0) THEN
          CALL mpi_allreduce(ave_p_change,tmp,1,mpi_nim_real,mpi_max,
     $         mpi_comm_world,ierror)
          ave_p_change=tmp
        ENDIF
      ENDIF
      ave_b_change=SQRT(ave_b_change)
c-----------------------------------------------------------------------
c     if changes are > specified limit, matrices will be updated during
c     the time step.  the updated n=0 grid vertex data is saved here
c     just for layer 0, because this data is only used for testing.
c     the matrices use quadrature-points data, which is now updated
c     during every step.
c
c     also save the difference between the total pressures and
c     the average pressures for the semi-implicit operators that are
c     used to stabilize nonlinear activity.
c-----------------------------------------------------------------------
      b0_changed=.false.
      p0_changed=.false.
      IF (ave_b_change>ave_change_limit.OR.
     $    ave_p_change>ave_change_limit.AND.beta>0.OR.first_call) THEN
        DO ibl=1,nbl
          IF (ilayer==0) be_n0(ibl)%arr=be(ibl)%arr(:,:,:,1)
          b0_changed=.true.
          IF (beta>0) THEN
            IF (ilayer==0) pres_n0(ibl)%arr=work3(ibl)%arr(:,:,:,1)
            p0_changed=.true.
          ENDIF
        ENDDO
      ENDIF
      first_call=.false.
c-----------------------------------------------------------------------
c     update si_nl_pres, which holds the difference between the maximum
c     pressure over the toroidal direction and the average pressure.
c     do this for both magnetic and internal pressures.
c
c     this is now performed at quadrature points, and only the first
c     quad point of old data per cell is saved for the matrix check.
c     this uses as much existing quadrature-point data as possible.
c-----------------------------------------------------------------------
      IF (si_fac_nl<=0) THEN
        CALL timer(timend_qp)
        time_qpfld=time_qpfld+timend_qp-timest_qp
        RETURN
      ENDIF
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          ng=rb(ibl)%ng
          btot=>rb(ibl)%qbe_tot%qpf
          ntot=>rb(ibl)%qnd_tot%qpf
          tiq=>rb(ibl)%qtion_eq%qpf
          teq=>rb(ibl)%qtele_eq%qpf
          cti=>rb(ibl)%qtion%qpf
          cte=>rb(ibl)%qtele%qpf
          beq=>rb(ibl)%qbe_eq%qpf
          ben0=>rb(ibl)%qbe_n0%qpf
          ndn0=>rb(ibl)%qnd_n0%qpf
          tin0=>rb(ibl)%qti_n0%qpf
          ten0=>rb(ibl)%qte_n0%qpf
          prq=>rb(ibl)%qpres_eq%qpf
          prn0=>rb(ibl)%qpres_n0%qpf
          nl_pr=>rb(ibl)%qsi_nl_pres%qpf
        ELSE
          ng=tb(ibl)%ng
          btot=>tb(ibl)%qbe_tot%qpf
          ntot=>tb(ibl)%qnd_tot%qpf
          tiq=>tb(ibl)%qtion_eq%qpf
          teq=>tb(ibl)%qtele_eq%qpf
          cti=>tb(ibl)%qtion%qpf
          cte=>tb(ibl)%qtele%qpf
          beq=>tb(ibl)%qbe_eq%qpf
          ben0=>tb(ibl)%qbe_n0%qpf
          ndn0=>tb(ibl)%qnd_n0%qpf
          tin0=>tb(ibl)%qti_n0%qpf
          ten0=>tb(ibl)%qte_n0%qpf
          prq=>tb(ibl)%qpres_eq%qpf
          prn0=>tb(ibl)%qpres_n0%qpf
          nl_pr=>tb(ibl)%qsi_nl_pres%qpf
        ENDIF
        m1=SIZE(tiq,2)
        m2=SIZE(tiq,3)
        mxb=SIZE(si_nl_pres(ibl)%arri,3)
        myb=SIZE(si_nl_pres(ibl)%arri,4)
        npol=mpsq_block(ibl)
        ipolst=ipqst_block(ibl)
        ipolen=ipqen_block(ibl)
c-----------------------------------------------------------------------
c       create and transform the effective temperature.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          ALLOCATE(cteff(1,m1,m2,nmodes))
          DO imode=1,nmodes
            cteff(1,:,:,imode)=cte(1,:,:,imode)+
     $                         cti(1,:,:,imode)/zeff
            IF (keff(imode)==0) cteff(1,:,:,imode)=
     $        cteff(1,:,:,imode)+teq(1,:,:)+tiq(1,:,:)/zeff
          ENDDO
          ALLOCATE(teff(1,npol,nphi))
          CALL fft_nim('inverse',m1*m2,npol,lphi,1_i4,cteff,teff,
     $                 dealiase)
        ENDIF
c-----------------------------------------------------------------------
c       determine the nonlinear pressures at every poloidal location.
c       also check the change for flagging matrix updates.
c
c       the second index in quadrature-point arrays is for quadrature
c       points, but it is folded into the poloidal index with ffts.
c-----------------------------------------------------------------------
        ALLOCATE(nltmp2(2,npol))
        nltmp2=0._r8
        ipol=0
        yloop: DO iy=1,myb
          DO ix=1,mxb
            DO ig=1,ng
              ipol=ipol+1
              IF (ipol<ipolst) CYCLE
              IF (ipol>ipolen) EXIT yloop
              ii=(iy-1)*mxb+ix
              ipg=ipol-ipolst+1
              IF (beta>0) THEN
                tmp=SUM((beq(:,ig,ii)+ben0(:,ig,ii))**2)+smallnum
                tm2=prq(1,ig,ii)+prn0(1,ig,ii)+smallnum
                DO ip=1,nphi
                  nltmp2(1,ipg)=MAX(nltmp2(1,ipg),
     $              ABS(SUM(btot(:,ipg,ip)**2)-tmp) )
                  nltmp2(2,ipg)=MAX(nltmp2(2,ipg),
     $              ABS(kboltz*ntot(1,ipg,ip)*teff(1,ipg,ip)-tm2))
                ENDDO
              ELSE
                tmp=SUM((beq(:,ig,ii)+ben0(:,ig,ii))**2)+smallnum
                DO ip=1,nphi
                  nltmp2(1,ipg)=MAX(nltmp2(1,ipg),
     $              ABS(SUM(btot(:,ipg,ip)**2)-tmp) )
                ENDDO
              ENDIF
              IF (ig==1.AND.beta>0) THEN
                nl_change=MAX( nl_change,
     $            ABS(si_nl_pres(ibl)%arri(1,1,ix,iy)
     $                -nltmp2(1,ipg))/tmp,
     $            ABS(si_nl_pres(ibl)%arri(2,1,ix,iy)
     $               -nltmp2(2,ipg))/tm2 )
              ELSE IF (ig==1) THEN
                nl_change=MAX( nl_change,
     $            ABS(si_nl_pres(ibl)%arri(1,1,ix,iy)
     $                -nltmp2(1,ipg))/tmp )
              ENDIF
            ENDDO
          ENDDO
        ENDDO yloop
c-----------------------------------------------------------------------
c       communicate the nonlinear pressures to all layers.
c-----------------------------------------------------------------------
        IF (nlayers==1) THEN
          nl_pr=RESHAPE(nltmp2,(/2,m1,m2/))
        ELSE
          npol=(m1*m2)/nlayers
          rep=MODULO(m1*m2,nlayers)
          DO il=0,nlayers-1
            ipolst=il*npol+1+MIN(rep,il)
            ipolen=(il+1)*npol+MIN(rep,il+1)
            rdispls(il)=2*(ipolst-1)
            rcounts(il)=2*(ipolen-ipolst+1)
          ENDDO
          CALL mpi_allgatherv(nltmp2(1,1),
     $         rcounts(ilayer),mpi_nim_real,nl_pr(1,1,1),
     $         rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
        ENDIF
        DEALLOCATE(nltmp2)
        IF (beta>0) DEALLOCATE(cteff,teff)
      ENDDO
c-----------------------------------------------------------------------
c     find the maximum nl_change over all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(nl_change,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        nl_change=tmp
      ENDIF
      IF (nl_change>ave_change_limit) THEN
        nl_changed=.true.
        DO ibl=1,nrbl
          mxb=SIZE(si_nl_pres(ibl)%arri,3)
          myb=SIZE(si_nl_pres(ibl)%arri,4)
          si_nl_pres(ibl)%arri(:,1,:,:)=
     $      RESHAPE(rb(ibl)%qsi_nl_pres%qpf(:,1,:),(/2,mxb,myb/))
        ENDDO
        DO ibl=nrbl+1,nbl
          mxb=SIZE(si_nl_pres(ibl)%arri,3)
          myb=SIZE(si_nl_pres(ibl)%arri,4)
          si_nl_pres(ibl)%arri(:,1,:,:)=
     $      RESHAPE(tb(ibl)%qsi_nl_pres%qpf(:,1,:),(/2,mxb,myb/))
        ENDDO
      ELSE
        nl_changed=.false.
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ave_field_check
c-----------------------------------------------------------------------
c     subprogram 4. ave_n_check
c     find the fractional change in average number density to determine
c     if matrices need updating in nonlinear runs.  the operations
c     performed here are similar to those in ave_field_check.
c-----------------------------------------------------------------------
      SUBROUTINE ave_n_check
      USE local
      USE fields
      USE input
      USE global
      USE mpi_nim
      USE pardata
      USE rblock
      USE tblock
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ierror,n,nh,nv,ni
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls
      REAL(r8) :: tmp,ave_n_change
      REAL(r8) :: timest_qp,timend_qp
      LOGICAL, SAVE :: first_call=.true.
c-----------------------------------------------------------------------
c     check the n=0 fields, which are in layer 0 if this is a parallel
c     run with multiple Fourier layers.
c-----------------------------------------------------------------------
      IF (.NOT.(continuity=='n=0 only'.OR.continuity=='full')) RETURN
      CALL timer(timest_qp)
      ave_n_change=0
c-----------------------------------------------------------------------
c     find the relative change for the n=0 + equilibrium fields.
c     just use grid vertex-centered information for this check.
c-----------------------------------------------------------------------
      IF (ilayer==0.AND..NOT.first_call) THEN
        DO ibl=1,nbl
          ave_n_change=MAX(ave_n_change,MAXVAL(ABS(nd_n0(ibl)%arr
     $      -REAL(nd(ibl)%arr(:,:,:,1),r8))
     $      /ABS(nd_n0(ibl)%arr+nd_eq(ibl)%arr+smallnum)))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate max change to all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(ave_n_change,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        ave_n_change=tmp
      ENDIF
c-----------------------------------------------------------------------
c     if changes are > specified limit, matrices will be updated during
c     the time step, and the updated n=0 fields are saved here for 
c     layer 0 only
c-----------------------------------------------------------------------
      n0_changed=.false.
      IF (ave_n_change>ave_change_limit.OR.first_call) THEN
        DO ibl=1,nbl
          IF (ilayer==0) nd_n0(ibl)%arr=nd(ibl)%arr(:,:,:,1)
          n0_changed=.true.
        ENDDO
      ENDIF
      first_call=.false.
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ave_n_check
c-----------------------------------------------------------------------
c     subprogram 5. ave_v_check
c     find the fractional change in average flow velocity to determine
c     if matrices need updating in nonlinear runs.  the operations
c     performed here are similar to those in ave_field_check.
c-----------------------------------------------------------------------
      SUBROUTINE ave_v_check
      USE local
      USE fields
      USE input
      USE global
      USE mpi_nim
      USE pardata
      USE rblock
      USE tblock
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ierror,n,nh,nv,ni
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls
      REAL(r8) :: tmp,ave_v_change
      REAL(r8) :: timest_qp,timend_qp
      LOGICAL, SAVE :: first_call=.true.
c-----------------------------------------------------------------------
c     check the n=0 fields, which are in layer 0 if this is a parallel
c     run with multiple Fourier layers.
c-----------------------------------------------------------------------
      IF (.NOT.impladv) RETURN
      CALL timer(timest_qp)
      ave_v_change=0
c-----------------------------------------------------------------------
c     find the relative change for the v=0 + equilibrium fields.
c     just use grid vertex-centered information for this check.
c-----------------------------------------------------------------------
      IF (ilayer==0.AND..NOT.first_call) THEN
        IF (eq_flow/='none') THEN
          DO ibl=1,nbl
            ave_v_change=MAX(ave_v_change,MAXVAL(SUM((ve_n0(ibl)%arr
     $         -REAL(ve(ibl)%arr(:,:,:,1),r8))**2,1)
     $         /MAX(SUM((ve_n0(ibl)%arr+ve_eq(ibl)%arr)**2,1),
     $              smallnum)))
          ENDDO
        ELSE
          DO ibl=1,nbl
            ave_v_change=MAX(ave_v_change,MAXVAL(SUM((ve_n0(ibl)%arr
     $         -REAL(ve(ibl)%arr(:,:,:,1),r8))**2,1)
     $         /MAX(SUM(ve_n0(ibl)%arr**2,1),smallnum)))
          ENDDO
        ENDIF
      ENDIF
      ave_v_change=SQRT(ave_v_change)
c-----------------------------------------------------------------------
c     communicate max change to all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(ave_v_change,tmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        ave_v_change=tmp
      ENDIF
c-----------------------------------------------------------------------
c     if changes are > specified limit, matrices will be updated during
c     the time step, and the updated n=0 fields are saved here for
c     layer 0.
c-----------------------------------------------------------------------
      v0_changed=.false.
      IF (ave_v_change>ave_change_limit.OR.first_call) THEN
        DO ibl=1,nbl
          IF (ilayer==0) ve_n0(ibl)%arr=ve(ibl)%arr(:,:,:,1)
          v0_changed=.true.
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     recomputing matrices for preconditioning due to every change in V
c     hasn't proven necessary.
c-----------------------------------------------------------------------
      v0_changed=.false.
      first_call=.false.
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ave_v_check
c-----------------------------------------------------------------------
c     subprogram 6. ave_eta_check
c     find the fractional change in the average electrical diffusivity
c     if matrices need updating in nonlinear runs.  the computation
c     performed here depends on eta_model.  for "eta full" or "chodura"
c     we compare the average of the toroidally dependent elecd with the 
c     stored value.  for "eta n=0 only" we compute eta based on the n=0
c     part of Te.
c-----------------------------------------------------------------------
      SUBROUTINE ave_eta_check
      USE local
      USE fields
      USE input
      USE global
      USE mpi_nim
      USE pardata
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,m1,m2,ierror,ig,ng
      INTEGER(i4) :: npol,ipolst,ipolen,rep,ix,iy,ipol,il
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls
      REAL(r8) :: tmp,ave_eta_change
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: eta0
      LOGICAL, SAVE :: first_call=.true.
      TYPE(vector_type) :: eta_ave(nbl)
c-----------------------------------------------------------------------
c     Always update the _n0 fields at the quadrature points.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      IF (eta_model=='eta n=0 only') THEN
        DO ibl=1,nrbl
          IF (ilayer==0) THEN
            rb(ibl)%qelecd_n0%qpf=
     $        MAX(elecd_min,MIN(elecd_max,elecd*
     $           (eta_ref_t/MAX(smallnum,(rb(ibl)%qtele_eq%qpf+
     $                  REAL(rb(ibl)%qtele%qpf(:,:,:,1)))))**1.5))
          ENDIF
          IF (nlayers>1) THEN
            CALL mpi_bcast(rb(ibl)%qelecd_n0%qpf,
     $                     SIZE(rb(ibl)%qelecd_n0%qpf),
     $                     mpi_nim_real,0,comm_mode,ierror)
          ENDIF
          ALLOCATE(eta_ave(ibl)%arri(SIZE(elecd_n0(ibl)%arri),1,1,1))
        ENDDO
        DO ibl=nrbl+1,nbl
          IF (ilayer==0) THEN
            tb(ibl)%qelecd_n0%qpf=
     $        MAX(elecd_min,MIN(elecd_max,elecd*
     $           (eta_ref_t/MAX(smallnum,(tb(ibl)%qtele_eq%qpf+
     $                  REAL(tb(ibl)%qtele%qpf(:,:,:,1)))))**1.5))
          ENDIF
          IF (nlayers>1) THEN
            CALL mpi_bcast(tb(ibl)%qelecd_n0%qpf,
     $                     SIZE(tb(ibl)%qelecd_n0%qpf),
     $                     mpi_nim_real,0,comm_mode,ierror)
          ENDIF
          ALLOCATE(eta_ave(ibl)%arri(SIZE(elecd_n0(ibl)%arri),1,1,1))
        ENDDO
c-----------------------------------------------------------------------
c     for full we need to gather computations from the different
c     subdivisions of the block.
c-----------------------------------------------------------------------
      ELSE IF (threedeta) THEN
        DO ibl=1,nrbl
          ALLOCATE(eta_ave(ibl)%arr(mpsq_block(ibl),1,1))
          ALLOCATE(eta_ave(ibl)%arri(SIZE(elecd_n0(ibl)%arri),1,1,1))
          eta_ave(ibl)%arr(:,1,1)=
     $      SUM(rb(ibl)%qelecd_phi%qpf(1,:,:),2)/nphi
        ENDDO
        DO ibl=nrbl+1,nbl
          ALLOCATE(eta_ave(ibl)%arr(mpsq_block(ibl),1,1))
          ALLOCATE(eta_ave(ibl)%arri(SIZE(elecd_n0(ibl)%arri),1,1,1))
          eta_ave(ibl)%arr(:,1,1)=
     $      SUM(tb(ibl)%qelecd_phi%qpf(1,:,:),2)/nphi
        ENDDO
c-----------------------------------------------------------------------
c       communicate the n=0 resistivity to all layers.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            eta0=>rb(ibl)%qelecd_n0%qpf
          ELSE
            eta0=>tb(ibl)%qelecd_n0%qpf
          ENDIF
          m1=SIZE(eta0,2)
          m2=SIZE(eta0,3)
          IF (nlayers==1) THEN
            eta0=RESHAPE(eta_ave(ibl)%arr,(/1,m1,m2/))
          ELSE
            npol=(m1*m2)/nlayers
            rep=MODULO(m1*m2,nlayers)
            DO il=0,nlayers-1
              ipolst=il*npol+1+MIN(rep,il)
              ipolen=(il+1)*npol+MIN(rep,il+1)
              rdispls(il)=ipolst-1
              rcounts(il)=ipolen-ipolst+1
            ENDDO
            CALL mpi_allgatherv(eta_ave(ibl)%arr(1,1,1),
     $           rcounts(ilayer),mpi_nim_real,eta0(1,1,1),
     $           rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
          ENDIF
          DEALLOCATE(eta_ave(ibl)%arr)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     Check for changes in the n=0 fields, which are in all layers.
c     this computation is performed at the first quadrature point.
c-----------------------------------------------------------------------
      ave_eta_change=0
      IF (.NOT.first_call) THEN
        DO ibl=1,nrbl
          eta_ave(ibl)%arri=RESHAPE
     $      (elecd_n0(ibl)%arri,(/SIZE(eta_ave(ibl)%arri),1,1,1/))
          ave_eta_change=MAX(ave_eta_change,
     $      MAXVAL(ABS(rb(ibl)%qelecd_n0%qpf(1,1,:)
     $                 -eta_ave(ibl)%arri(:,1,1,1))
     $        /(rb(ibl)%qelecd_n0%qpf(1,1,:)+smallnum)))
        ENDDO
        DO ibl=nrbl+1,nbl
          eta_ave(ibl)%arri=RESHAPE
     $      (elecd_n0(ibl)%arri,(/SIZE(eta_ave(ibl)%arri),1,1,1/))
          ave_eta_change=MAX(ave_eta_change,
     $      MAXVAL(ABS(tb(ibl)%qelecd_n0%qpf(1,1,:)
     $                 -eta_ave(ibl)%arri(:,1,1,1))
     $        /(tb(ibl)%qelecd_n0%qpf(1,1,:)+smallnum)))
        ENDDO
c-----------------------------------------------------------------------
c       communicate max change to all processors.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          CALL mpi_allreduce(ave_eta_change,tmp,1,mpi_nim_real,mpi_max,
     $         mpi_comm_world,ierror)
          ave_eta_change=tmp
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     if changes are > specified limit, matrices will be updated during
c     the time step.
c-----------------------------------------------------------------------
      eta_changed=.false.
      IF (ave_eta_change>ave_change_limit.OR.first_call) THEN
        eta_changed=.true.
        DO ibl=1,nrbl
          m1=SIZE(elecd_n0(ibl)%arri,3)
          m2=SIZE(elecd_n0(ibl)%arri,4)
          elecd_n0(ibl)%arri=RESHAPE
     $      (rb(ibl)%qelecd_n0%qpf(1,1,:),(/1,1,m1,m2/))
        ENDDO
        DO ibl=nrbl+1,nbl
          m1=SIZE(elecd_n0(ibl)%arri,3)
          m2=SIZE(elecd_n0(ibl)%arri,4)
          elecd_n0(ibl)%arri=RESHAPE
     $      (tb(ibl)%qelecd_n0%qpf(1,1,:),(/1,1,m1,m2/))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     deallocate temporary storage.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        DEALLOCATE(eta_ave(ibl)%arri)
      ENDDO
      first_call=.false.
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ave_eta_check
c-----------------------------------------------------------------------
c     subprogram 7. find_eta_t
c     compute electron temperature-dependent electrical diffusivity
c     at the quadrature points as a function of toroidal grid
c     position.
c-----------------------------------------------------------------------
      SUBROUTINE find_eta_t
      USE local
      USE fft_mod
      USE fields
      USE global
      USE input
      USE physdat
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ig,lx,ly,im
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: tn0
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: teff
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     electron temperature is saved by Fourier component only.  add
c     equilibrium temperature temporarily.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO ibl=1,nrbl
        lx=SIZE(rb(ibl)%qtele%qpf,2)
        ly=SIZE(rb(ibl)%qtele%qpf,3)
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            ALLOCATE(tn0(1,lx,ly))
            tn0=rb(ibl)%qtele%qpf(:,:,:,im)
            rb(ibl)%qtele%qpf(:,:,:,im)=
     $        rb(ibl)%qtele%qpf(:,:,:,im)+rb(ibl)%qtele_eq%qpf
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               rb(ibl)%qtele%qpf,rb(ibl)%qelecd_phi%qpf,dealiase)
        rb(ibl)%qelecd_phi%qpf=MAX( elecd_min, MIN( elecd_max,
     $    elecd*(eta_ref_t/MAX(smallnum,rb(ibl)%qelecd_phi%qpf))**1.5) )
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            rb(ibl)%qtele%qpf(:,:,:,im)=tn0
            DEALLOCATE(tn0)
            EXIT
          ENDIF
        ENDDO
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
            tn0=tb(ibl)%qtele%qpf(:,:,:,im)
            tb(ibl)%qtele%qpf(:,:,:,im)=
     $        tb(ibl)%qtele%qpf(:,:,:,im)+tb(ibl)%qtele_eq%qpf
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               tb(ibl)%qtele%qpf,tb(ibl)%qelecd_phi%qpf,dealiase)
        tb(ibl)%qelecd_phi%qpf=MAX( elecd_min, MIN( elecd_max,
     $    elecd*(eta_ref_t/MAX(smallnum,tb(ibl)%qelecd_phi%qpf))**1.5) )
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            tb(ibl)%qtele%qpf(:,:,:,im)=tn0
            DEALLOCATE(tn0)
            EXIT
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     the chodura model adds a drift-speed-dependent model to
c     the Spitzer resistivity.  the electrical diffusivity from this
c     model, as implemented here is
c
c     elecd_chodura*SQRT(ndens/n)*
c       {1-EXP[-f_chodura*SQRT(j**2*mtot/gamma*p*n)]}
c-----------------------------------------------------------------------
      IF (eta_model/="chodura") THEN
        CALL timer(timend_qp)
        time_qpfld=time_qpfld+timend_qp-timest_qp
        RETURN
      ENDIF
      DO ibl=1,nrbl
        lx=SIZE(rb(ibl)%qtele%qpf,2)
        ly=SIZE(rb(ibl)%qtele%qpf,3)
        ALLOCATE(teff(1,mpsq_block(ibl),nphi))
        rb(ibl)%qwork3%qpf=rb(ibl)%qtele%qpf+rb(ibl)%qtion%qpf/zeff
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            rb(ibl)%qwork3%qpf(:,:,:,im)=
     $        rb(ibl)%qwork3%qpf(:,:,:,im)+
     $        rb(ibl)%qtele_eq%qpf+rb(ibl)%qtion_eq%qpf/zeff
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               rb(ibl)%qwork3%qpf,teff,dealiase)
        teff=kboltz*MAX(teff,smallnum)
        rb(ibl)%qelecd_phi%qpf(1,:,:)=rb(ibl)%qelecd_phi%qpf(1,:,:)+
     $    elecd_chodura*SQRT(ndens/rb(ibl)%qnd_tot%qpf(1,:,:))*
     $    (1._r8-EXP(-f_chodura*
     $               SQRT(SUM(rb(ibl)%qja_tot%qpf**2,1)*mtot/
     $               (gamma*teff(1,:,:)))/
     $               (rb(ibl)%qnd_tot%qpf(1,:,:)*elementary_q)))
        DEALLOCATE(teff)
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        lx=SIZE(tb(ibl)%qtele%qpf,2)
        ly=SIZE(tb(ibl)%qtele%qpf,3)
        ALLOCATE(teff(1,mpsq_block(ibl),nphi))
        tb(ibl)%qwork3%qpf=tb(ibl)%qtele%qpf+tb(ibl)%qtion%qpf/zeff
        DO im=1,nmodes
          IF (keff(im)==0) THEN
            tb(ibl)%qwork3%qpf(:,:,:,im)=
     $        tb(ibl)%qwork3%qpf(:,:,:,im)+
     $        tb(ibl)%qtele_eq%qpf+tb(ibl)%qtion_eq%qpf/zeff
            EXIT
          ENDIF
        ENDDO
        CALL fft_nim('inverse',lx*ly,mpsq_block(ibl),lphi,1_i4,
     $               tb(ibl)%qwork3%qpf,teff,dealiase)
        teff=kboltz*MAX(teff,smallnum)
        tb(ibl)%qelecd_phi%qpf(1,:,:)=tb(ibl)%qelecd_phi%qpf(1,:,:)+
     $    elecd_chodura*SQRT(ndens/tb(ibl)%qnd_tot%qpf(1,:,:))*
     $    (1._r8-EXP(-f_chodura*
     $               SQRT(SUM(tb(ibl)%qja_tot%qpf**2,1)*mtot/
     $               (gamma*teff(1,:,:)))/
     $               (tb(ibl)%qnd_tot%qpf(1,:,:)*elementary_q)))
        DEALLOCATE(teff)
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE find_eta_t
c-----------------------------------------------------------------------
c     subprogram 8. find_bb
c     compute the (b_hat)(b_hat) dyad at the quadrature points.
c-----------------------------------------------------------------------
      SUBROUTINE find_bb
      USE local
      USE pardata
      USE mpi_nim
      USE fields
      USE global
      USE input
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,imode,ierror,iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmp
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     linear computations use equilibrium fields only.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      IF (.NOT.nonlinear) THEN
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            rb(ibl)%qbb%qpf(1:3,:,:)=rb(ibl)%qbe_eq%qpf**2
            rb(ibl)%qbb%qpf(4  ,:,:)=rb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               rb(ibl)%qbe_eq%qpf(2,:,:)
            rb(ibl)%qbb%qpf(5  ,:,:)=rb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               rb(ibl)%qbe_eq%qpf(3,:,:)
            rb(ibl)%qbb%qpf(6  ,:,:)=rb(ibl)%qbe_eq%qpf(2,:,:)*
     $                               rb(ibl)%qbe_eq%qpf(3,:,:)
            DO iq=1,6
              rb(ibl)%qbb%qpf(iq,:,:)=rb(ibl)%qbb%qpf(iq,:,:)/
     $          SUM(rb(ibl)%qbe_eq%qpf**2,1)
            ENDDO
          ELSE
            tb(ibl)%qbb%qpf(1:3,:,:)=tb(ibl)%qbe_eq%qpf**2
            tb(ibl)%qbb%qpf(4  ,:,:)=tb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               tb(ibl)%qbe_eq%qpf(2,:,:)
            tb(ibl)%qbb%qpf(5  ,:,:)=tb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               tb(ibl)%qbe_eq%qpf(3,:,:)
            tb(ibl)%qbb%qpf(6  ,:,:)=tb(ibl)%qbe_eq%qpf(2,:,:)*
     $                               tb(ibl)%qbe_eq%qpf(3,:,:)
            DO iq=1,6
              tb(ibl)%qbb%qpf(iq,:,:)=tb(ibl)%qbb%qpf(iq,:,:)/
     $          SUM(tb(ibl)%qbe_eq%qpf**2,1)
            ENDDO
          ENDIF
        ENDDO
        CALL timer(timend_qp)
        time_qpfld=time_qpfld+timend_qp-timest_qp
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     work directly with the data stored at quadrature points.
c     n=0 contribution has equilibrium field added, and n/=0
c     contribution results from +-|n|.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          rb(ibl)%qbb%qpf=0
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              rb(ibl)%qbb%qpf(1:3,:,:)=rb(ibl)%qbb%qpf(1:3,:,:)+
     $          (rb(ibl)%qbe%qpf(1:3,:,:,imode)+rb(ibl)%qbe_eq%qpf)**2
              rb(ibl)%qbb%qpf(4,:,:)=rb(ibl)%qbb%qpf(4,:,:)+
     $          (rb(ibl)%qbe%qpf(1,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(1,:,:))*
     $          (rb(ibl)%qbe%qpf(2,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(2,:,:))
              rb(ibl)%qbb%qpf(5,:,:)=rb(ibl)%qbb%qpf(5,:,:)+
     $          (rb(ibl)%qbe%qpf(1,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(1,:,:))*
     $          (rb(ibl)%qbe%qpf(3,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(3,:,:))
              rb(ibl)%qbb%qpf(6,:,:)=rb(ibl)%qbb%qpf(6,:,:)+
     $          (rb(ibl)%qbe%qpf(2,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(2,:,:))*
     $          (rb(ibl)%qbe%qpf(3,:,:,imode)+
     $           rb(ibl)%qbe_eq%qpf(3,:,:))
            ELSE
              rb(ibl)%qbb%qpf(1:3,:,:)=rb(ibl)%qbb%qpf(1:3,:,:)+
     $              2*rb(ibl)%qbe%qpf(1:3,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(1:3,:,:,imode))
              rb(ibl)%qbb%qpf(4,:,:)=rb(ibl)%qbb%qpf(4,:,:)+
     $                rb(ibl)%qbe%qpf(1,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(2,:,:,imode))+
     $                rb(ibl)%qbe%qpf(2,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(1,:,:,imode))
              rb(ibl)%qbb%qpf(5,:,:)=rb(ibl)%qbb%qpf(5,:,:)+
     $                rb(ibl)%qbe%qpf(1,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(3,:,:,imode))+
     $                rb(ibl)%qbe%qpf(3,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(1,:,:,imode))
              rb(ibl)%qbb%qpf(6,:,:)=rb(ibl)%qbb%qpf(6,:,:)+
     $                rb(ibl)%qbe%qpf(2,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(3,:,:,imode))+
     $                rb(ibl)%qbe%qpf(3,:,:,imode)*
     $          CONJG(rb(ibl)%qbe%qpf(2,:,:,imode))
            ENDIF
          ENDDO
          ALLOCATE(tmp(6,rb(ibl)%ng,rb(ibl)%mx*rb(ibl)%my))
          IF (nlayers>1) THEN
            CALL mpi_allreduce(rb(ibl)%qbb%qpf(1,1,1),tmp,SIZE(tmp),
     $                         mpi_nim_real,mpi_sum,comm_mode,ierror)
            rb(ibl)%qbb%qpf=tmp
          ENDIF
          DO iq=1,6
            rb(ibl)%qbb%qpf(iq,:,:)=rb(ibl)%qbb%qpf(iq,:,:)/
     $        SUM((rb(ibl)%qbe_eq%qpf+rb(ibl)%qbe_n0%qpf)**2,1)
          ENDDO
          DEALLOCATE(tmp)
c-----------------------------------------------------------------------
c       same for tblocks.
c-----------------------------------------------------------------------
        ELSE
          tb(ibl)%qbb%qpf=0
          DO imode=1,nmodes
            IF (keff(imode)==0) THEN
              tb(ibl)%qbb%qpf(1:3,:,:)=tb(ibl)%qbb%qpf(1:3,:,:)+
     $          (tb(ibl)%qbe%qpf(1:3,:,:,imode)+tb(ibl)%qbe_eq%qpf)**2
              tb(ibl)%qbb%qpf(4,:,:)=tb(ibl)%qbb%qpf(4,:,:)+
     $          (tb(ibl)%qbe%qpf(1,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(1,:,:))*
     $          (tb(ibl)%qbe%qpf(2,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(2,:,:))
              tb(ibl)%qbb%qpf(5,:,:)=tb(ibl)%qbb%qpf(5,:,:)+
     $          (tb(ibl)%qbe%qpf(1,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(1,:,:))*
     $          (tb(ibl)%qbe%qpf(3,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(3,:,:))
              tb(ibl)%qbb%qpf(6,:,:)=tb(ibl)%qbb%qpf(6,:,:)+
     $          (tb(ibl)%qbe%qpf(2,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(2,:,:))*
     $          (tb(ibl)%qbe%qpf(3,:,:,imode)+
     $           tb(ibl)%qbe_eq%qpf(3,:,:))
            ELSE
              tb(ibl)%qbb%qpf(1:3,:,:)=tb(ibl)%qbb%qpf(1:3,:,:)+
     $              2*tb(ibl)%qbe%qpf(1:3,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(1:3,:,:,imode))
              tb(ibl)%qbb%qpf(4,:,:)=tb(ibl)%qbb%qpf(4,:,:)+
     $                tb(ibl)%qbe%qpf(1,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(2,:,:,imode))+
     $                tb(ibl)%qbe%qpf(2,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(1,:,:,imode))
              tb(ibl)%qbb%qpf(5,:,:)=tb(ibl)%qbb%qpf(5,:,:)+
     $                tb(ibl)%qbe%qpf(1,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(3,:,:,imode))+
     $                tb(ibl)%qbe%qpf(3,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(1,:,:,imode))
              tb(ibl)%qbb%qpf(6,:,:)=tb(ibl)%qbb%qpf(6,:,:)+
     $                tb(ibl)%qbe%qpf(2,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(3,:,:,imode))+
     $                tb(ibl)%qbe%qpf(3,:,:,imode)*
     $          CONJG(tb(ibl)%qbe%qpf(2,:,:,imode))
            ENDIF
          ENDDO
          ALLOCATE(tmp(6,tb(ibl)%ng,tb(ibl)%mcell))
          IF (nlayers>1) THEN
            CALL mpi_allreduce(tb(ibl)%qbb%qpf(1,1,1),tmp,SIZE(tmp),
     $                         mpi_nim_real,mpi_sum,comm_mode,ierror)
            tb(ibl)%qbb%qpf=tmp
          ENDIF
          DO iq=1,6
            tb(ibl)%qbb%qpf(iq,:,:)=tb(ibl)%qbb%qpf(iq,:,:)/
     $        SUM((tb(ibl)%qbe_eq%qpf+tb(ibl)%qbe_n0%qpf)**2,1)
          ENDDO
          DEALLOCATE(tmp)
        ENDIF
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE find_bb
c-----------------------------------------------------------------------
c     subprogram 9. ave_kappa_check
c     find the fractional change in the average parallel thermal
c     conductivity.  this is only used for p_model="aniso_plltdep"
c     or "aniso_tdep" in nonlinear computations.
c-----------------------------------------------------------------------
      SUBROUTINE ave_kappa_check
      USE local
      USE fields
      USE input
      USE global
      USE mpi_nim
      USE pardata
      USE fft_mod
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,m1,m2,ierror,ig,ng,nkp,ikp
      INTEGER(i4) :: npol,ipolst,ipolen,rep,ix,iy,ipol,il
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls
      REAL(r8) :: tmp,ave_kap_change
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          plli,prpi,plle,prpe
      LOGICAL, SAVE :: first_call=.true.

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: kap_ave1
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: kap_ave2
c-----------------------------------------------------------------------
c     first determine the number of conductivities that need to be
c     averaged and communicated.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      nkp=1
      IF (separate_pe) nkp=nkp+1
      IF (.NOT.closure_n0_only) nkp=2*nkp
c-----------------------------------------------------------------------
c     for 3D conductivity, find and save the average to avoid
c     recomputing later.
c
c     the arr part of the kap_ave structure is used for the average in
c     the non-layered decomposition, while the arri part is used for the
c     layered decomposition.
c
c     perpendicular thermal conductivity terms have been added for the
c     Braginskii magnetization model.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          ng=rb(ibl)%ng
          plli=>rb(ibl)%qkappli_phi%qpf
          prpi=>rb(ibl)%qkaprpi_phi%qpf
          plle=>rb(ibl)%qkapple_phi%qpf
          prpe=>rb(ibl)%qkaprpe_phi%qpf
        ELSE
          ng=tb(ibl)%ng
          plli=>tb(ibl)%qkappli_phi%qpf
          prpi=>tb(ibl)%qkaprpi_phi%qpf
          plle=>tb(ibl)%qkapple_phi%qpf
          prpe=>tb(ibl)%qkaprpe_phi%qpf
        ENDIF
        ALLOCATE(kap_ave1(nkp,mpsq_block(ibl)))
        kap_ave1(1,:)=SUM(plli(1,:,:),2)/nphi
        ikp=2
        IF (.NOT.closure_n0_only) THEN
          kap_ave1(ikp,:)=SUM(prpi(1,:,:),2)/nphi
          ikp=ikp+1
        ENDIF
        IF (separate_pe) THEN
          kap_ave1(ikp,:)=SUM(plle(1,:,:),2)/nphi
          ikp=ikp+1
          IF (.NOT.closure_n0_only) 
     $      kap_ave1(ikp,:)=SUM(prpe(1,:,:),2)/nphi
        ENDIF
c-----------------------------------------------------------------------
c       communicate the n=0 conductivies to all layers.  this now
c       includes perpendicular thermal conductivity terms.
c-----------------------------------------------------------------------
        m1=SIZE(kappli_n0(ibl)%arri,3)
        m2=SIZE(kappli_n0(ibl)%arri,4)
        ALLOCATE(kap_ave2(nkp,ng,m1*m2))
        IF (nlayers==1) THEN
          kap_ave2=RESHAPE(kap_ave1,(/nkp,ng,m1*m2/))
        ELSE
          npol=(ng*m1*m2)/nlayers
          rep=MODULO(ng*m1*m2,nlayers)
          DO il=0,nlayers-1
            ipolst=il*npol+1+MIN(rep,il)
            ipolen=(il+1)*npol+MIN(rep,il+1)
            rdispls(il)=nkp*(ipolst-1)
            rcounts(il)=nkp*(ipolen-ipolst+1)
          ENDDO
          CALL mpi_allgatherv(kap_ave1(1,1),
     $         rcounts(ilayer),mpi_nim_real,kap_ave2(1,1,1),
     $         rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
        ENDIF
        IF (ibl<=nrbl) THEN
          plli=>rb(ibl)%qkappli_n0%qpf
          prpi=>rb(ibl)%qkaprpi_n0%qpf
          plle=>rb(ibl)%qkapple_n0%qpf
          prpe=>rb(ibl)%qkaprpe_n0%qpf
        ELSE
          plli=>tb(ibl)%qkappli_n0%qpf
          prpi=>tb(ibl)%qkaprpi_n0%qpf
          plle=>tb(ibl)%qkapple_n0%qpf
          prpe=>tb(ibl)%qkaprpe_n0%qpf
        ENDIF
        plli(1,:,:)=kap_ave2(1,:,:)
        ikp=2
        IF (.NOT.closure_n0_only) THEN
          prpi(1,:,:)=kap_ave2(ikp,:,:)
          ikp=ikp+1
        ENDIF
        IF (separate_pe) THEN
          plle(1,:,:)=kap_ave2(ikp,:,:)
          ikp=ikp+1
          IF (.NOT.closure_n0_only) prpe(1,:,:)=kap_ave2(ikp,:,:)
        ENDIF
        DEALLOCATE(kap_ave1,kap_ave2)
      ENDDO
c-----------------------------------------------------------------------
c     Check for changes in the n=0 fields, which are in all layers.
c     this computation is performed at the first quadrature point.
c-----------------------------------------------------------------------
      ave_kap_change=0
      IF (.NOT.first_call) THEN
        DO ibl=1,nrbl
          m1=rb(ibl)%mx*rb(ibl)%my
          ave_kap_change=MAX(ave_kap_change,
     $      MAXVAL(ABS(rb(ibl)%qkappli_n0%qpf(1,1,:)
     $                 -RESHAPE(kappli_n0(ibl)%arri,(/m1/)))
     $        /(rb(ibl)%qkappli_n0%qpf(1,1,:)+smallnum)))
          IF (separate_pe) THEN
            ave_kap_change=MAX(ave_kap_change,
     $        MAXVAL(ABS(rb(ibl)%qkapple_n0%qpf(1,1,:)
     $                   -RESHAPE(kapple_n0(ibl)%arri,(/m1/)))
     $          /(rb(ibl)%qkapple_n0%qpf(1,1,:)+smallnum)))
          ENDIF
        ENDDO
        DO ibl=nrbl+1,nbl
          m1=tb(ibl)%mcell
          ave_kap_change=MAX(ave_kap_change,
     $      MAXVAL(ABS(tb(ibl)%qkappli_n0%qpf(1,1,:)
     $                 -RESHAPE(kappli_n0(ibl)%arri,(/m1/)))
     $        /(tb(ibl)%qkappli_n0%qpf(1,1,:)+smallnum)))
          IF (separate_pe) THEN
            ave_kap_change=MAX(ave_kap_change,
     $        MAXVAL(ABS(tb(ibl)%qkapple_n0%qpf(1,1,:)
     $                   -RESHAPE(kapple_n0(ibl)%arri,(/m1/)))
     $          /(tb(ibl)%qkapple_n0%qpf(1,1,:)+smallnum)))
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       communicate max change to all processors.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          CALL mpi_allreduce(ave_kap_change,tmp,1,mpi_nim_real,mpi_max,
     $         mpi_comm_world,ierror)
          ave_kap_change=tmp
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     if changes are > specified limit, matrices will be updated during
c     the time step, and the updated n=0 fields are saved here for all
c     layers.
c-----------------------------------------------------------------------
      kpll_changed=.false.
      IF (ave_kap_change>ave_change_limit.OR.first_call) THEN
        kpll_changed=.true.
        DO ibl=1,nrbl
          m1=rb(ibl)%mx
          m2=rb(ibl)%my
          kappli_n0(ibl)%arri(1,1,:,:)=
     $      RESHAPE(rb(ibl)%qkappli_n0%qpf(1,1,:),(/m1,m2/))
          IF (separate_pe)
     $      kapple_n0(ibl)%arri(1,1,:,:)=
     $        RESHAPE(rb(ibl)%qkapple_n0%qpf(1,1,:),(/m1,m2/))
        ENDDO
        DO ibl=nrbl+1,nbl
          kappli_n0(ibl)%arri(1,1,:,1)=tb(ibl)%qkappli_n0%qpf(1,1,:)
          IF (separate_pe)
     $      kapple_n0(ibl)%arri(1,1,:,1)=tb(ibl)%qkapple_n0%qpf(1,1,:)
        ENDDO
      ENDIF
      first_call=.false.
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ave_kappa_check
c-----------------------------------------------------------------------
c     subprogram 10. find_kappa_t
c     compute temperature-dependent thermal diffusivities at the
c     quadrature points as a function of toroidal grid position.
c
c     this subroutine now acts as a shell to call the appropriate
c     closure models located in closure_model.f.
c-----------------------------------------------------------------------
      SUBROUTINE find_kappa_t
      USE local
      USE input
      USE closure_model_mod
      USE time
      IMPLICIT NONE

      REAL(r8) :: timest_qp,timend_qp

      CALL timer(timest_qp)
      IF (closure_model=='std kprp n0') THEN
        CALL kappa_standard_n0
      ELSE IF (closure_model=='standard') THEN
        CALL kappa_standard_phi
      ELSE IF (closure_model=='braginskii') THEN
        CALL kappa_braginskii_phi
      ELSE IF (closure_model=='k2') THEN
        CALL kappa_k2_phi
      ELSE
        CALL nim_stop('Invalid closure_model option')
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE find_kappa_t
c-----------------------------------------------------------------------
c     subprogram 11. n_store
c     save number density at the quadrature points.  the saved data
c     will either be the density at the end of the step or the average
c     of the beginning and ending values, as determined by the passed
c     parameter.
c-----------------------------------------------------------------------
      SUBROUTINE n_store(nchoice,newdart)
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE pardata
      USE mpi_nim
      USE fft_mod
      USE time
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: nchoice
      LOGICAL, INTENT(OUT) :: newdart

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ntmp
      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_g1
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          dart,upwc,ndtot,vetot,upave,ndeqr,ndeqz
      REAL(r8) :: ndfl,ndexp,dmax,dtmp,umax
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), SAVE :: dmax_old=-1._r8,umax_old=-1._r8

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_nd
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             ndn,ndr,ndz,vept

      INTEGER(i4) :: ibl,ncx,ncy,ng,ip,ig,ix,iy,ipolst,ipolen,npol,rep,
     $               il,ierror,nc,im,iq
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls

      TYPE(vector_type), DIMENSION(:), POINTER :: vtmp
c-----------------------------------------------------------------------
c     interface block for qp0_bcast.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE qp0_bcast(threedata,twodata,nl)
        USE local
        COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
        REAL(r8), DIMENSION(:,:,:) :: twodata
        INTEGER(i4), INTENT(IN) :: nl
        END SUBROUTINE qp0_bcast
      END INTERFACE
c-----------------------------------------------------------------------
c     part 1: updating quadrature-point storage.
c
c     determine what needs to be saved.  use the work3 temporary
c     space to find the average.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      SELECT CASE(nchoice)
      CASE('average')
        DO ibl=1,nbl
          work3(ibl)=nd_old(ibl)
          CALL vector_add(work3(ibl),nd(ibl),v1fac=0.5_r8,v2fac=0.5_r8)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qnd,rb(ibl))
            IF (nonlinear) THEN
              CALL qp_fft_save(rb(ibl)%qnd%qpf,rb(ibl)%qnd_tot%qpf,
     $                         rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                         1_i4,rb(ibl)%ng,rb(ibl)%qnd_eq%qpf)
            ENDIF
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qnd,tb(ibl))
            IF (nonlinear) THEN
              CALL qp_fft_save(tb(ibl)%qnd%qpf,tb(ibl)%qnd_tot%qpf,
     $                         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                         1_i4,tb(ibl)%ng,tb(ibl)%qnd_eq%qpf)
            ENDIF
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c     this is for nonlinear cases where the n-advance may be
c     recomputed based on a new artificial diffusivity, so nd quad-pt
c     storage is preserved as n at the start of the time-step.  use
c     the qwork3 space for temporary quad-point storage to get the
c     new nd_tot for the diffusivity computation itself.  also use the
c     vtmp vector to save the old qnd_tot data.
c-----------------------------------------------------------------------
      CASE('check ave')
        ALLOCATE(vtmp(nbl))
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            ALLOCATE(vtmp(ibl)%arr(1,mpsq_block(ibl),nphi))
            vtmp(ibl)%arr=rb(ibl)%qnd_tot%qpf
          ELSE
            ALLOCATE(vtmp(ibl)%arr(1,mpsq_block(ibl),nphi))
            vtmp(ibl)%arr=tb(ibl)%qnd_tot%qpf
          ENDIF
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%nd,rb(ibl)%qwork3,rb(ibl))
            CALL qp_fft_save(rb(ibl)%qwork3%qpf,rb(ibl)%qnd_tot%qpf,
     $                       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                       1_i4,rb(ibl)%ng,rb(ibl)%qnd_eq%qpf)
          ELSE
            CALL tblock_qp_update(tb(ibl)%nd,tb(ibl)%qwork3,tb(ibl))
            CALL qp_fft_save(tb(ibl)%qwork3%qpf,tb(ibl)%qnd_tot%qpf,
     $                       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                       1_i4,tb(ibl)%ng,tb(ibl)%qnd_eq%qpf)
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c     store the field at the end of the step.
c-----------------------------------------------------------------------
      CASE('end')
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%nd,rb(ibl)%qnd,rb(ibl))
            IF (nonlinear) THEN
              CALL qp_fft_save(rb(ibl)%qnd%qpf,rb(ibl)%qnd_tot%qpf,
     $                         rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                         1_i4,rb(ibl)%ng,rb(ibl)%qnd_eq%qpf)
            ENDIF
          ELSE
            CALL tblock_qp_update(tb(ibl)%nd,tb(ibl)%qnd,tb(ibl))
            IF (nonlinear) THEN
              CALL qp_fft_save(tb(ibl)%qnd%qpf,tb(ibl)%qnd_tot%qpf,
     $                         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                         1_i4,tb(ibl)%ng,tb(ibl)%qnd_eq%qpf)
            ENDIF
          ENDIF
        ENDDO
      END SELECT

c-----------------------------------------------------------------------
c     part 2:  local artificial particle diffusivity based on n.
c
c     check if the minimum number density (over toroidal angle) falls
c     below nd_floor, and if so, create an enhancement to the artificial
c     diffusivity to fill-in the hole.
c
c     first find the minimum--collect over all layers if necessary.
c-----------------------------------------------------------------------
      dmax=-1._r8
      IF (nonlinear.AND.nd_floor>0.AND.nd_diff>0.AND.
     $    nchoice/='average') THEN
        ndfl=nd_floor*ndens
        ndexp=nd_exp*ndens
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            ncx=rb(ibl)%mx
            ncy=rb(ibl)%my
            ng=rb(ibl)%ng
            dart=>rb(ibl)%qdart%qpf
            ndtot=>rb(ibl)%qnd_tot%qpf
          ELSE
            ncx=tb(ibl)%mcell
            ncy=1
            ng=tb(ibl)%ng
            dart=>tb(ibl)%qdart%qpf
            ndtot=>tb(ibl)%qnd_tot%qpf
          ENDIF
          nc=ncx*ncy
          IF (nlayers==1) THEN
            ALLOCATE(ntmp(mpsq_block(ibl),1))
            ntmp(:,1)=ndtot(1,:,1)
            DO ip=2,nphi
              ntmp(:,1)=MIN(ntmp(:,1),ndtot(1,:,ip))
            ENDDO
            dart=RESHAPE(ntmp,(/1,ng,ncx*ncy/))
          ELSE
            ALLOCATE(ntmp(mpsq_block(ibl),1))
            ntmp(:,1)=ndtot(1,:,1)
            DO ip=2,nphi
              ntmp(:,1)=MIN(ntmp(:,1),ndtot(1,:,ip))
            ENDDO
            npol=(ng*nc)/nlayers
            rep=MODULO(ng*nc,nlayers)
            DO il=0,nlayers-1
              ipolst=il*npol+1+MIN(rep,il)
              ipolen=(il+1)*npol+MIN(rep,il+1)
              rdispls(il)=ipolst-1
              rcounts(il)=ipolen-ipolst+1
            ENDDO
            CALL mpi_allgatherv(ntmp(1,1),
     $           rcounts(ilayer),mpi_nim_real,dart(1,1,1),
     $           rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
          ENDIF
          DEALLOCATE(ntmp)
c-----------------------------------------------------------------------
c         now determine the 2D arificial diffusivity enhancement (times
c         dt) based on the toroidal-minimum number density.
c-----------------------------------------------------------------------
          dart=nd_dart_fac*cross_section/(mx*my)*
     $      MAX(0._r8,TANH( (ndfl-dart)/ndexp ) )
          dmax=MAX(dmax,MAXVAL(dart))
        ENDDO
        IF (nprocs_layer>1) THEN
          CALL mpi_allreduce(dmax,dtmp,1,mpi_nim_real,mpi_max,
     $                       comm_layer,ierror)
          dmax=dtmp
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     set the newdart flag if n-total is dropping below nd_floor*ndens
c     at any point.
c-----------------------------------------------------------------------
      IF (dmax>0._r8.OR.dmax_old>0._r8) THEN
        newdart=.true.
      ELSE
        newdart=.false.
      ENDIF
      IF (nchoice/='average') dmax_old=dmax

c-----------------------------------------------------------------------
c     part 3: local particle diffusivity based on v.grad(n).
c-----------------------------------------------------------------------
      umax=-1._r8
      IF (nonlinear.AND.impladv.AND.nd_dart_upw>0.AND.
     $    nchoice/='average') THEN
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            ncx=rb(ibl)%mx
            ncy=rb(ibl)%my
            ng=rb(ibl)%ng
            upwc=>rb(ibl)%qupw_phi%qpf
            upave=>rb(ibl)%qupw_n0%qpf
            ndtot=>rb(ibl)%qnd_tot%qpf
            ndeqr=>rb(ibl)%qnd_eq%qpfr
            ndeqz=>rb(ibl)%qnd_eq%qpfz
            vetot=>rb(ibl)%qve_tot%qpf
            vept=>rb(ibl)%qve%qpf
            bigr=>rb(ibl)%bigr
            IF (nchoice=='check ave') THEN
              ndn=>rb(ibl)%qwork3%qpf
              ndr=>rb(ibl)%qwork3%qpfr
              ndz=>rb(ibl)%qwork3%qpfz
            ELSE
              ndn=>rb(ibl)%qnd%qpf
              ndr=>rb(ibl)%qnd%qpfr
              ndz=>rb(ibl)%qnd%qpfz
            ENDIF
          ELSE
            ncx=tb(ibl)%mcell
            ncy=1
            ng=tb(ibl)%ng
            upwc=>tb(ibl)%qupw_phi%qpf
            upave=>tb(ibl)%qupw_n0%qpf
            ndtot=>tb(ibl)%qnd_tot%qpf
            ndeqr=>tb(ibl)%qnd_eq%qpfr
            ndeqz=>tb(ibl)%qnd_eq%qpfz
            vetot=>tb(ibl)%qve_tot%qpf
            vept=>tb(ibl)%qve%qpf
            bigr=>tb(ibl)%tgeom%bigr
            IF (nchoice=='check ave') THEN
              ndn=>tb(ibl)%qwork3%qpf
              ndr=>tb(ibl)%qwork3%qpfr
              ndz=>tb(ibl)%qwork3%qpfz
            ELSE
              ndn=>tb(ibl)%qnd%qpf
              ndr=>tb(ibl)%qnd%qpfr
              ndz=>tb(ibl)%qnd%qpfz
            ENDIF
          ENDIF
          nc=ncx*ncy
          ALLOCATE(real_g1(4,mpsq_block(ibl),nphi))
          ALLOCATE(grad_nd(4,ng,nc,nmodes))
          ALLOCATE(ntmp(mpsq_block(ibl),nphi))
c-----------------------------------------------------------------------
c         collect and transform grad(nd) and v-tilde.grad(nd_eq).
c         to minimize the number of mpi calls.
c-----------------------------------------------------------------------
          DO im=1,nmodes
            DO ix=1,nc
              DO ig=1,ng
                grad_nd(1,ig,ix,im)=ndr(1,ig,ix,im)
                grad_nd(2,ig,ix,im)=ndz(1,ig,ix,im)
                grad_nd(3,ig,ix,im)=
     $            (0,1)*keff(im)*ndn(1,ig,ix,im)/bigr(ig,ix)
                grad_nd(4,ig,ix,im)=
     $            vept(1,ig,ix,im)*ndeqr(1,ig,ix)+
     $            vept(2,ig,ix,im)*ndeqz(1,ig,ix)
              ENDDO
            ENDDO
          ENDDO
          CALL fft_nim('inverse',ng*nc,mpsq_block(ibl),lphi,4_i4,
     $                 grad_nd,real_g1,dealiase)
c-----------------------------------------------------------------------
c         now determine the coefficient at all quadrature points.
c         it is
c
c         jac*(V.grad(n))**2/( ntot**2 ),
c
c         where V.grad(n) includes all but Veq.grad(neq).  the 2d
c         Jacobian provides an area per element.  the coefficient
c         is also multiplied by nd_dart_upw*dt in the integrand routine.
c-----------------------------------------------------------------------
          npol=mpsq_block(ibl)
          ipolst=ipqst_block(ibl)
          ipolen=ipqen_block(ibl)
          dtmp=dt**2+smallnum
          ntmp=(SUM(vetot*real_g1(1:3,:,:),1)+real_g1(4,:,:))**2/
     $            ( ndtot(1,:,:)**2 + smallnum )
          ntmp=MIN(ntmp,upw_limit/dtmp)+
     $         MAX(0._r8,TANH((nd_floor_upw-ndtot(1,:,:))/
     $                         nd_width_upw))/dtmp
          umax=MAX(umax,MAXVAL(ntmp))
          IF (ibl<=nrbl) THEN
            DO ip=1,npol
              ix=(ip+ipolst-2)/ng+1
              ig=ip+ipolst-1-(ix-1)*ng
              ntmp(ip,:)=ntmp(ip,:)*rb(ibl)%jac2d(ig,ix)
            ENDDO
          ELSE
            DO ip=1,npol
              ix=(ip+ipolst-2)/ng+1
              ntmp(ip,:)=ntmp(ip,:)*tb(ibl)%tgeom%area(ix)
            ENDDO
          ENDIF
          IF (nchoice=='check ave') THEN
c-TMP
c           upwc(1,:,:)=0.5_r8*( upwc(1,:,:) + ntmp )
            upwc(1,:,:)=MAX( upwc(1,:,:) , ntmp )
          ELSE
            upwc(1,:,:)=ntmp
          ENDIF
          DEALLOCATE(real_g1,grad_nd,ntmp)
c-----------------------------------------------------------------------
c         find the toroidal average coefficient for preconditioning.
c-----------------------------------------------------------------------
          IF (nlayers==1) THEN
            ALLOCATE(ntmp(mpsq_block(ibl),1))
            ntmp(:,1)=SUM(upwc(1,:,:),2)/nphi
            upave=RESHAPE(ntmp,(/1,ng,nc/))
          ELSE
            ALLOCATE(ntmp(mpsq_block(ibl),1))
            ntmp(:,1)=SUM(upwc(1,:,:),2)/nphi
            npol=(ng*nc)/nlayers
            rep=MODULO(ng*nc,nlayers)
            DO il=0,nlayers-1
              ipolst=il*npol+1+MIN(rep,il)
              ipolen=(il+1)*npol+MIN(rep,il+1)
              rdispls(il)=ipolst-1
              rcounts(il)=ipolen-ipolst+1
            ENDDO
            CALL mpi_allgatherv(ntmp(1,1),
     $           rcounts(ilayer),mpi_nim_real,upave(1,1,1),
     $           rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
          ENDIF
          DEALLOCATE(ntmp)
        ENDDO
c-----------------------------------------------------------------------
c       find max (dt*V.grad(n)/n)**2 over all blocks.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          CALL mpi_allreduce(umax,dtmp,1,mpi_nim_real,mpi_max,
     $                       mpi_comm_world,ierror)
          umax=dtmp
        ENDIF
        umax=dt**2*umax
      ENDIF
c-----------------------------------------------------------------------
c     set the newdart flag if the maximum (dt*V.grad(n)/ntot)**2 exceeds
c     ave_change_limit.
c-----------------------------------------------------------------------
      IF (umax>ave_change_limit.OR.umax_old>ave_change_limit) 
     $  newdart=.true.
      IF (nchoice/='average') umax_old=umax

c-----------------------------------------------------------------------
c     part 4: reset fields for the 'check ave' step.
c
c     if the first n-advance for a nonlinear implicit leapfrog leads
c     to a new diffusivity coefficient, leave the nodal values as the
c     new n for predicting a solution to the recompute, and leave the
c     quadrature point storage of nd as the old values for the rhs
c     of the recompute.  nd_tot is not used for advancing n, so it's
c     storage doesn't matter.
c
c     if a recompute isn't needed, save the average in qnd and qnd_tot.
c-----------------------------------------------------------------------
      IF (nchoice=='check ave') THEN
        IF (newdart) THEN
          DO ibl=1,nbl
            DEALLOCATE(vtmp(ibl)%arr)
          ENDDO
          DEALLOCATE(vtmp)
        ELSE
          DO ibl=1,nbl
            work3(ibl)=nd_old(ibl)
            CALL vector_add(work3(ibl),nd(ibl),v1fac=0.5_r8,
     $                      v2fac=0.5_r8)
            IF (ibl<=nrbl) THEN
              CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qnd,rb(ibl))
              rb(ibl)%qnd_tot%qpf=0.5_r8*(rb(ibl)%qnd_tot%qpf+
     $                                    vtmp(ibl)%arr)
            ELSE
              CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qnd,tb(ibl))
              tb(ibl)%qnd_tot%qpf=0.5_r8*(tb(ibl)%qnd_tot%qpf+
     $                                    vtmp(ibl)%arr)
            ENDIF
            DEALLOCATE(vtmp(ibl)%arr)
          ENDDO
          DEALLOCATE(vtmp)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     part 5: update storage of symmetric component.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        DO ibl=1,nrbl
          CALL qp0_bcast(rb(ibl)%qnd%qpf,rb(ibl)%qnd_n0%qpf,nlayers)
          CALL qp0_bcast(rb(ibl)%qnd%qpfr,rb(ibl)%qnd_n0%qpfr,nlayers)
          CALL qp0_bcast(rb(ibl)%qnd%qpfz,rb(ibl)%qnd_n0%qpfz,nlayers)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL qp0_bcast(tb(ibl)%qnd%qpf,tb(ibl)%qnd_n0%qpf,nlayers)
          CALL qp0_bcast(tb(ibl)%qnd%qpfr,tb(ibl)%qnd_n0%qpfr,nlayers)
          CALL qp0_bcast(tb(ibl)%qnd%qpfz,tb(ibl)%qnd_n0%qpfz,nlayers)
        ENDDO
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE n_store
c-----------------------------------------------------------------------
c     subprogram 12. temp_store
c     save temperatures at the quadrature points.  the saved data
c     will either be from the end of the step or the average
c     of the beginning and ending values, as determined by the passed
c     parameter.
c-----------------------------------------------------------------------
      SUBROUTINE temp_store(tchoice,newkart)
      USE local
      USE fields
      USE global
      USE input
      USE physdat
      USE rblock
      USE tblock
      USE pardata
      USE mpi_nim
      USE fft_mod
      USE math_tran
      USE time
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: tchoice
      LOGICAL, INTENT(OUT) :: newkart

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ttmp
      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_g1
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          dart,upwc,vetot,vet2,
     $          upave,tieq,teeq,teqr,teqz,ndtot,ndq,nd0,jatot
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), TARGET :: umaxi,umaxe
      REAL(r8), TARGET, SAVE :: umaxi_old=-1._r8,umaxe_old=-1._r8
      REAL(r8), POINTER :: umax,umax_old

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: grad_t
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             temp,tempr,tempz,vept,ve2,bee,ber,bez

      INTEGER(i4) :: ibl,ncx,ncy,ng,ip,ig,ix,iy,ipolst,ipolen,npol,rep,
     $               il,ierror,nc,im,iq
      INTEGER(i4), DIMENSION(0:nlayers-1) :: rcounts,rdispls

      REAL(r8) :: deni,dene,tfrac,jfac,utmp,dt2
      LOGICAL :: do_n0
c-----------------------------------------------------------------------
c     interface block for qp0_bcast.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE qp0_bcast(threedata,twodata,nl)
        USE local
        COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
        REAL(r8), DIMENSION(:,:,:) :: twodata
        INTEGER(i4), INTENT(IN) :: nl
        END SUBROUTINE qp0_bcast
      END INTERFACE
c-----------------------------------------------------------------------
c     part 1: updating quadrature-point storage.
c
c     determine what needs to be saved.  use the work3 temporary
c     space to find the average.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      SELECT CASE(tchoice(1:7))
      CASE('ion ave')
        DO ibl=1,nbl
          work3(ibl)=ti_old(ibl)
          CALL vector_add(work3(ibl),tion(ibl),v1fac=0.5_r8,
     $                    v2fac=0.5_r8)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qtion,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qtion,tb(ibl))
          ENDIF
        ENDDO
      CASE('ele ave')
        DO ibl=1,nbl
          work3(ibl)=te_old(ibl)
          CALL vector_add(work3(ibl),tele(ibl),v1fac=0.5_r8,
     $                    v2fac=0.5_r8)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qtele,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qtele,tb(ibl))
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c     store the field at the end of the step.
c-----------------------------------------------------------------------
      CASE('ion end','ion chk')
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%tion,rb(ibl)%qtion,rb(ibl))
            IF (nonlinear.AND.gyr_visc>0.AND.tchoice(5:7)=='end')
     $        CALL qp_fft_save(rb(ibl)%qtion%qpf,rb(ibl)%qti_tot%qpf,
     $                         rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                         1_i4,rb(ibl)%ng,rb(ibl)%qtion_eq%qpf)
          ELSE
            CALL tblock_qp_update(tb(ibl)%tion,tb(ibl)%qtion,tb(ibl))
            IF (nonlinear.AND.gyr_visc>0.AND.tchoice(5:7)=='end')
     $        CALL qp_fft_save(tb(ibl)%qtion%qpf,tb(ibl)%qti_tot%qpf,
     $                         tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                         1_i4,tb(ibl)%ng,tb(ibl)%qtion_eq%qpf)
          ENDIF
        ENDDO
      CASE('ele end','ele chk')
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%tele,rb(ibl)%qtele,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%tele,tb(ibl)%qtele,tb(ibl))
          ENDIF
        ENDDO
      END SELECT

c-----------------------------------------------------------------------
c     part 2: local thermal diffusivity based on v.grad(T).
c     if electron temperature is computed separately, use the
c     appropriate species velocity.
c-----------------------------------------------------------------------
      IF (tchoice(1:3)=='ion') THEN
        umax=>umaxi
        umax_old=>umaxi_old
        jfac=coefjvi
      ELSE
        umax=>umaxe
        umax_old=>umaxe_old
        jfac=coefjve
      ENDIF

      umax=-1._r8
      IF (nonlinear.AND.impladv.AND.t_dart_upw>0.AND.
     $     tchoice(5:7)/='ave'.AND.
     $    (tchoice(1:3)=='ion'.OR.separate_pe)) THEN
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            ncx=rb(ibl)%mx
            ncy=rb(ibl)%my
            ng=rb(ibl)%ng
            bigr=>rb(ibl)%bigr
            ndtot=>rb(ibl)%qnd_tot%qpf
            ndq=>rb(ibl)%qnd_eq%qpf
            nd0=>rb(ibl)%qnd_n0%qpf
            bee=>rb(ibl)%qbe%qpf
            ber=>rb(ibl)%qbe%qpfr
            bez=>rb(ibl)%qbe%qpfz
            vetot=>rb(ibl)%qve_tot%qpf
            vept=>rb(ibl)%qve%qpf
            jatot=>rb(ibl)%qja_tot%qpf
            IF (tchoice(1:3)=='ion') THEN
              temp=>rb(ibl)%qtion%qpf
              tempr=>rb(ibl)%qtion%qpfr
              tempz=>rb(ibl)%qtion%qpfz
              tieq=>rb(ibl)%qtion_eq%qpf
              teqr=>rb(ibl)%qtion_eq%qpfr
              teqz=>rb(ibl)%qtion_eq%qpfz
              upwc=>rb(ibl)%qupti_phi%qpf
              upave=>rb(ibl)%qupti_n0%qpf
            ELSE
              temp=>rb(ibl)%qtele%qpf
              tempr=>rb(ibl)%qtele%qpfr
              tempz=>rb(ibl)%qtele%qpfz
              tieq=>rb(ibl)%qtele_eq%qpf
              teqr=>rb(ibl)%qtele_eq%qpfr
              teqz=>rb(ibl)%qtele_eq%qpfz
              upwc=>rb(ibl)%qupte_phi%qpf
              upave=>rb(ibl)%qupte_n0%qpf
            ENDIF
          ELSE
            ncx=tb(ibl)%mcell
            ncy=1
            ng=tb(ibl)%ng
            bigr=>tb(ibl)%tgeom%bigr
            ndtot=>tb(ibl)%qnd_tot%qpf
            ndq=>tb(ibl)%qnd_eq%qpf
            nd0=>tb(ibl)%qnd_n0%qpf
            bee=>tb(ibl)%qbe%qpf
            ber=>tb(ibl)%qbe%qpfr
            bez=>tb(ibl)%qbe%qpfz
            vetot=>tb(ibl)%qve_tot%qpf
            vept=>tb(ibl)%qve%qpf
            jatot=>tb(ibl)%qja_tot%qpf
            IF (tchoice(1:3)=='ion') THEN
              temp=>tb(ibl)%qtion%qpf
              tempr=>tb(ibl)%qtion%qpfr
              tempz=>tb(ibl)%qtion%qpfz
              tieq=>tb(ibl)%qtion_eq%qpf
              teqr=>tb(ibl)%qtion_eq%qpfr
              teqz=>tb(ibl)%qtion_eq%qpfz
              upwc=>tb(ibl)%qupti_phi%qpf
              upave=>tb(ibl)%qupti_n0%qpf
            ELSE
              temp=>tb(ibl)%qtele%qpf
              tempr=>tb(ibl)%qtele%qpfr
              tempz=>tb(ibl)%qtele%qpfz
              tieq=>tb(ibl)%qtele_eq%qpf
              teqr=>tb(ibl)%qtele_eq%qpfr
              teqz=>tb(ibl)%qtele_eq%qpfz
              upwc=>tb(ibl)%qupte_phi%qpf
              upave=>tb(ibl)%qupte_n0%qpf
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c         compute the species velocity (or an approximate) if using
c         separate electron temperature.
c-----------------------------------------------------------------------
          IF (jfac/=0._r8) THEN
            ALLOCATE(ve2(3,ng,ncx*ncy,nmodes))
            ALLOCATE(vet2(3,mpsq_block(ibl),nphi))
            CALL math_curl(nmodes,keff,geom,bigr,
     $                     bee,ber,bez,ve2,jfac/mu0)
            DO im=1,nmodes
              DO iq=1,3
                ve2(iq,:,:,im)=ve2(iq,:,:,im)/
     $            (ndq(1,:,:)+nd0(1,:,:))+vept(iq,:,:,im)
              ENDDO
            ENDDO
            DO iq=1,3
              vet2(iq,:,:)=vetot(iq,:,:)+jfac*jatot(iq,:,:)/ndtot(1,:,:)
            ENDDO
          ELSE
            ve2=>vept
            vet2=>vetot
          ENDIF
c-----------------------------------------------------------------------
c         allocate temporary space for this block
c-----------------------------------------------------------------------
          nc=ncx*ncy
          ALLOCATE(real_g1(5,mpsq_block(ibl),nphi))
          ALLOCATE(grad_t(5,ng,ncx*ncy,nmodes))
          ALLOCATE(ttmp(mpsq_block(ibl),nphi))
c-----------------------------------------------------------------------
c         collect and transform grad(T), v-tilde.grad(T_eq), and T.
c         to minimize the number of mpi calls.
c-----------------------------------------------------------------------
          DO im=1,nmodes
            DO ix=1,nc
              DO ig=1,ng
                grad_t(1,ig,ix,im)=tempr(1,ig,ix,im)
                grad_t(2,ig,ix,im)=tempz(1,ig,ix,im)
                grad_t(3,ig,ix,im)=
     $            (0,1)*keff(im)*temp(1,ig,ix,im)/bigr(ig,ix)
                grad_t(4,ig,ix,im)=
     $            ve2(1,ig,ix,im)*teqr(1,ig,ix)+
     $            ve2(2,ig,ix,im)*teqz(1,ig,ix)
                grad_t(5,ig,ix,im)=temp(1,ig,ix,im)
              ENDDO
            ENDDO
          ENDDO
          CALL fft_nim('inverse',ng*ncx*ncy,mpsq_block(ibl),lphi,5_i4,
     $                 grad_t,real_g1,dealiase)
c-----------------------------------------------------------------------
c         now determine the coefficient at all quadrature points.
c         it is
c
c         jac*(V.grad(T))**2/( T**2 ),
c
c         where V.grad(T) includes all but Veq.grad(Teq).  the 2d
c         Jacobian provides an area per element.  the coefficient is
c         also multiplied by t_dart_upw*dt in the integrand routine.
c-----------------------------------------------------------------------
          npol=mpsq_block(ibl)
          ipolst=ipqst_block(ibl)
          ipolen=ipqen_block(ibl)
          dt2=dt**2+smallnum
          DO ip=1,npol
            ix=(ip+ipolst-2)/ng+1
            ig=ip+ipolst-1-(ix-1)*ng
            real_g1(5,ip,:)=real_g1(5,ip,:)+tieq(1,ig,ix)
          ENDDO
          ttmp=(SUM(vet2*real_g1(1:3,:,:),1)+real_g1(4,:,:))**2/
     $            ( real_g1(5,:,:)**2 + smallnum )
          ttmp=MIN(ttmp,upw_limit/dt2)+
     $         MAX(0._r8,TANH((t_floor_upw-real_g1(5,:,:))/
     $                         t_width_upw))/dt2
          umax=MAX(umax,MAXVAL(ttmp))
          IF (ibl<=nrbl) THEN
            DO ip=1,npol
              ix=(ip+ipolst-2)/ng+1
              ig=ip+ipolst-1-(ix-1)*ng
              ttmp(ip,:)=ttmp(ip,:)*rb(ibl)%jac2d(ig,ix)
            ENDDO
          ELSE
            DO ip=1,npol
              ix=(ip+ipolst-2)/ng+1
              ttmp(ip,:)=ttmp(ip,:)*tb(ibl)%tgeom%area(ix)
            ENDDO
          ENDIF
          IF (tchoice(5:7)=='chk') THEN
c-TMP
c           upwc(1,:,:)=0.5_r8*( upwc(1,:,:) + ttmp )
            upwc(1,:,:)=MAX( upwc(1,:,:) , ttmp )
          ELSE
            upwc(1,:,:)=ttmp
          ENDIF

          DEALLOCATE(real_g1,grad_t,ttmp)
          IF (jfac/=0._r8) DEALLOCATE(ve2,vet2)
c-----------------------------------------------------------------------
c         find the toroidal average coefficient for preconditioning.
c-----------------------------------------------------------------------
          IF (nlayers==1) THEN
            ALLOCATE(ttmp(mpsq_block(ibl),1))
            ttmp(:,1)=SUM(upwc(1,:,:),2)/nphi
            upave=RESHAPE(ttmp,(/1,ng,nc/))
          ELSE
            ALLOCATE(ttmp(mpsq_block(ibl),1))
            ttmp(:,1)=SUM(upwc(1,:,:),2)/nphi
            npol=(ng*nc)/nlayers
            rep=MODULO(ng*nc,nlayers)
            DO il=0,nlayers-1
              ipolst=il*npol+1+MIN(rep,il)
              ipolen=(il+1)*npol+MIN(rep,il+1)
              rdispls(il)=ipolst-1
              rcounts(il)=ipolen-ipolst+1
            ENDDO
            CALL mpi_allgatherv(ttmp(1,1),
     $           rcounts(ilayer),mpi_nim_real,upave(1,1,1),
     $           rcounts,rdispls,mpi_nim_real,comm_mode,ierror)
          ENDIF
          DEALLOCATE(ttmp)
        ENDDO
c-----------------------------------------------------------------------
c       find max (dt*V.grad(T)/T)**2 over all blocks.
c-----------------------------------------------------------------------
        IF (nprocs>1) THEN
          CALL mpi_allreduce(umax,utmp,1,mpi_nim_real,mpi_max,
     $                       mpi_comm_world,ierror)
          umax=utmp
        ENDIF
        umax=dt**2*umax
      ENDIF
c-----------------------------------------------------------------------
c     set the newkart flag if the maximum (dt*V.grad(T)/T)**2 exceeds
c     ave_change_limit.
c-----------------------------------------------------------------------
      IF (umax>ave_change_limit.OR.umax_old>ave_change_limit) THEN
        newkart=.true.
      ELSE
        newkart=.false.
      ENDIF
      IF (tchoice(5:7)/='ave') umax_old=umax

c-----------------------------------------------------------------------
c     part 3: reset fields for the check step.
c
c     if the first T-advance for a nonlinear implicit leapfrog leads
c     to a new diffusivity coefficient, reset the quadrature point
c     storage to the old values.  leave the nodal values as the result
c     of the first pass to help make a guess for the solution of the
c     recompute.
c
c     a recompute of T will be done if a check step is called.
c-----------------------------------------------------------------------
      IF (tchoice(1:7)=='ion chk') THEN
        DO ibl=1,nbl
          work3(ibl)=ti_old(ibl)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qtion,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qtion,tb(ibl))
          ENDIF
        ENDDO
      ELSE IF (tchoice(1:7)=='ele chk') THEN
        DO ibl=1,nbl
          work3(ibl)=te_old(ibl)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qtele,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qtele,tb(ibl))
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     part 4: update storage of symmetric component, provided this call
c     is not for diagnostic purposes only.
c-----------------------------------------------------------------------
      do_n0=.true.
      IF (LEN_TRIM(tchoice)>=13) THEN
        IF (tchoice(9:13)=='diagn') do_n0=.false.
      ENDIF
      IF (nonlinear.AND.do_n0) THEN
        IF (tchoice(1:3)=='ion') THEN
          DO ibl=1,nrbl
            CALL qp0_bcast(rb(ibl)%qtion%qpf,rb(ibl)%qti_n0%qpf,nlayers)
            CALL qp0_bcast(rb(ibl)%qtion%qpfr,rb(ibl)%qti_n0%qpfr,
     $                     nlayers)
            CALL qp0_bcast(rb(ibl)%qtion%qpfz,rb(ibl)%qti_n0%qpfz,
     $                     nlayers)
            rb(ibl)%qti_n0%qpf=rb(ibl)%qti_n0%qpf+rb(ibl)%qtion_eq%qpf
            rb(ibl)%qti_n0%qpfr=rb(ibl)%qti_n0%qpfr+
     $                          rb(ibl)%qtion_eq%qpfr
            rb(ibl)%qti_n0%qpfz=rb(ibl)%qti_n0%qpfz+
     $                          rb(ibl)%qtion_eq%qpfz
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL qp0_bcast(tb(ibl)%qtion%qpf,tb(ibl)%qti_n0%qpf,nlayers)
            CALL qp0_bcast(tb(ibl)%qtion%qpfr,tb(ibl)%qti_n0%qpfr,
     $                     nlayers)
            CALL qp0_bcast(tb(ibl)%qtion%qpfz,tb(ibl)%qti_n0%qpfz,
     $                     nlayers)
            tb(ibl)%qti_n0%qpf=tb(ibl)%qti_n0%qpf+tb(ibl)%qtion_eq%qpf
            tb(ibl)%qti_n0%qpfr=tb(ibl)%qti_n0%qpfr+
     $                          tb(ibl)%qtion_eq%qpfr
            tb(ibl)%qti_n0%qpfz=tb(ibl)%qti_n0%qpfz+
     $                          tb(ibl)%qtion_eq%qpfz
          ENDDO
        ELSE
          DO ibl=1,nrbl
            CALL qp0_bcast(rb(ibl)%qtele%qpf,rb(ibl)%qte_n0%qpf,nlayers)
            CALL qp0_bcast(rb(ibl)%qtele%qpfr,rb(ibl)%qte_n0%qpfr,
     $                     nlayers)
            CALL qp0_bcast(rb(ibl)%qtele%qpfz,rb(ibl)%qte_n0%qpfz,
     $                     nlayers)
            rb(ibl)%qte_n0%qpf=rb(ibl)%qte_n0%qpf+rb(ibl)%qtele_eq%qpf
            rb(ibl)%qte_n0%qpfr=rb(ibl)%qte_n0%qpfr+
     $                          rb(ibl)%qtele_eq%qpfr
            rb(ibl)%qte_n0%qpfz=rb(ibl)%qte_n0%qpfz+
     $                          rb(ibl)%qtele_eq%qpfz
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL qp0_bcast(tb(ibl)%qtele%qpf,tb(ibl)%qte_n0%qpf,nlayers)
            CALL qp0_bcast(tb(ibl)%qtele%qpfr,tb(ibl)%qte_n0%qpfr,
     $                     nlayers)
            CALL qp0_bcast(tb(ibl)%qtele%qpfz,tb(ibl)%qte_n0%qpfz,
     $                     nlayers)
            tb(ibl)%qte_n0%qpf=tb(ibl)%qte_n0%qpf+tb(ibl)%qtele_eq%qpf
            tb(ibl)%qte_n0%qpfr=tb(ibl)%qte_n0%qpfr+
     $                          tb(ibl)%qtele_eq%qpfr
            tb(ibl)%qte_n0%qpfz=tb(ibl)%qte_n0%qpfz+
     $                          tb(ibl)%qtele_eq%qpfz
          ENDDO
        ENDIF
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE temp_store
c-----------------------------------------------------------------------
c     subprogram 13. vcom_store
c     save center-of-mass velocity at the quadrature points.  also,
c     create and save its gradient if necessary for the stress
c     computations.
c-----------------------------------------------------------------------
      SUBROUTINE vcom_store(v_choice)
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE pardata
      USE fft_mod
      USE math_tran
      USE mpi_nim
      USE time
      IMPLICIT NONE

      CHARACTER(*) :: v_choice

      INTEGER(i4) :: ibl,mxb,myb,ig,ng,nv,ierror,im
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     interface block for qp0_bcast.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE qp0_bcast(threedata,twodata,nl)
        USE local
        COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
        REAL(r8), DIMENSION(:,:,:) :: twodata
        INTEGER(i4), INTENT(IN) :: nl
        END SUBROUTINE qp0_bcast
      END INTERFACE
c-----------------------------------------------------------------------
c     all cases need to save Fourier components.  total velocity as
c     a function of phi is needed for nonlinear implicit advection, and
c     total grad(V) as a function of phi is needed for nonlinear
c     implicit advection or nonlinear anisotropic viscosity.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO ibl=1,nrbl
        IF (v_choice=='standard')
     $    CALL rblock_qp_update(rb(ibl)%ve,rb(ibl)%qve,rb(ibl))
        IF (nonlinear) THEN
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          ng=rb(ibl)%ng
          ALLOCATE(ctmp(9,ng,mxb*myb,nmodes))
          IF (v_choice=='advect it') THEN
            CALL math_grad(nmodes,keff,3_i4,geom,
     $                     rb(ibl)%qve%qpf +0.5_r8*rb(ibl)%qwork1%qpf,
     $                     rb(ibl)%qve%qpfr+0.5_r8*rb(ibl)%qwork1%qpfr,
     $                     rb(ibl)%qve%qpfz+0.5_r8*rb(ibl)%qwork1%qpfz,
     $                     ctmp,bigr=rb(ibl)%bigr)
          ELSE
            CALL math_grad(nmodes,keff,3_i4,geom,rb(ibl)%qve%qpf,
     $                   rb(ibl)%qve%qpfr,rb(ibl)%qve%qpfz,
     $                   ctmp,bigr=rb(ibl)%bigr)
          ENDIF
          IF (eq_flow/='none') THEN
            IF (v_choice=='advect it') THEN
              CALL qp_fft_save(rb(ibl)%qve%qpf+
     $                         0.5_r8*rb(ibl)%qwork1%qpf,
     $                         rb(ibl)%qve_tot%qpf,
     $                         mxb,myb,mpsq_block(ibl),3_i4,ng,
     $                         rb(ibl)%qve_eq%qpf)
            ELSE
              CALL qp_fft_save(rb(ibl)%qve%qpf,rb(ibl)%qve_tot%qpf,
     $                         mxb,myb,mpsq_block(ibl),3_i4,ng,
     $                         rb(ibl)%qve_eq%qpf)
            ENDIF
            DO im=1,nmodes
              IF (keff(im)==0) THEN
                ctmp(:,:,:,im)=ctmp(:,:,:,im)+rb(ibl)%qgrdveq%qpf
                EXIT
              ENDIF
            ENDDO
          ELSE
            IF (v_choice=='advect it') THEN
              CALL qp_fft_noeq_save(rb(ibl)%qve%qpf+
     $                              0.5_r8*rb(ibl)%qwork1%qpf,
     $                              rb(ibl)%qve_tot%qpf,
     $                              mxb,myb,mpsq_block(ibl),3_i4,ng)
            ELSE
              CALL qp_fft_noeq_save(rb(ibl)%qve%qpf,rb(ibl)%qve_tot%qpf,
     $                              mxb,myb,mpsq_block(ibl),3_i4,ng)
            ENDIF
          ENDIF
          CALL qp_fft_noeq_save(ctmp,rb(ibl)%qgrdv%qpf,mxb,myb,
     $                          mpsq_block(ibl),9_i4,ng)
          DEALLOCATE(ctmp)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        IF (v_choice=='standard')
     $    CALL tblock_qp_update(tb(ibl)%ve,tb(ibl)%qve,tb(ibl))
        IF (nonlinear) THEN
          mxb=tb(ibl)%mcell
          myb=1_i4
          ng=tb(ibl)%ng
          ALLOCATE(ctmp(9,ng,mxb*myb,nmodes))
          IF (v_choice=='advect it') THEN
            CALL math_grad(nmodes,keff,3_i4,geom,
     $                     tb(ibl)%qve%qpf +0.5_r8*tb(ibl)%qwork1%qpf,
     $                     tb(ibl)%qve%qpfr+0.5_r8*tb(ibl)%qwork1%qpfr,
     $                     tb(ibl)%qve%qpfz+0.5_r8*tb(ibl)%qwork1%qpfz,
     $                     ctmp,bigr=tb(ibl)%tgeom%bigr)
          ELSE
            CALL math_grad(nmodes,keff,3_i4,geom,tb(ibl)%qve%qpf,
     $                     tb(ibl)%qve%qpfr,tb(ibl)%qve%qpfz,
     $                     ctmp,bigr=tb(ibl)%tgeom%bigr)
          ENDIF
          IF (eq_flow/='none') THEN
            IF (v_choice=='advect it') THEN
              CALL qp_fft_save(tb(ibl)%qve%qpf+
     $                         0.5_r8*tb(ibl)%qwork1%qpf,
     $                         tb(ibl)%qve_tot%qpf,
     $                         mxb,myb,mpsq_block(ibl),3_i4,ng,
     $                         tb(ibl)%qve_eq%qpf)
            ELSE
              CALL qp_fft_save(tb(ibl)%qve%qpf,tb(ibl)%qve_tot%qpf,
     $                         mxb,myb,mpsq_block(ibl),3_i4,ng,
     $                         tb(ibl)%qve_eq%qpf)
            ENDIF
            DO im=1,nmodes
              IF (keff(im)==0) THEN
                ctmp(:,:,:,im)=ctmp(:,:,:,im)+tb(ibl)%qgrdveq%qpf
                EXIT
              ENDIF
            ENDDO
          ELSE
            IF (v_choice=='advect it') THEN
              CALL qp_fft_noeq_save(tb(ibl)%qve%qpf+
     $                              0.5_r8*tb(ibl)%qwork1%qpf,
     $                              tb(ibl)%qve_tot%qpf,
     $                              mxb,myb,mpsq_block(ibl),3_i4,ng)
            ELSE
              CALL qp_fft_noeq_save(tb(ibl)%qve%qpf,tb(ibl)%qve_tot%qpf,
     $                              mxb,myb,mpsq_block(ibl),3_i4,ng)
            ENDIF
          ENDIF
          CALL qp_fft_noeq_save(ctmp,tb(ibl)%qgrdv%qpf,mxb,myb,
     $                          mpsq_block(ibl),9_i4,ng)
          DEALLOCATE(ctmp)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     broadcast the symmetric component here--derivatives are needed.
c     also broadcast asymmetric components that are used for coupling
c     in the preconditioner.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.impladv.AND.v_choice=='standard') THEN
        DO ibl=1,nrbl
          CALL qp0_bcast(rb(ibl)%qve%qpf,rb(ibl)%qve_n0%qpf,nlayers)
          CALL qp0_bcast(rb(ibl)%qve%qpfr,rb(ibl)%qve_n0%qpfr,nlayers)
          CALL qp0_bcast(rb(ibl)%qve%qpfz,rb(ibl)%qve_n0%qpfz,nlayers)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL qp0_bcast(tb(ibl)%qve%qpf,tb(ibl)%qve_n0%qpf,nlayers)
          CALL qp0_bcast(tb(ibl)%qve%qpfr,tb(ibl)%qve_n0%qpfr,nlayers)
          CALL qp0_bcast(tb(ibl)%qve%qpfz,tb(ibl)%qve_n0%qpfz,nlayers)
        ENDDO
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vcom_store
c-----------------------------------------------------------------------
c     subprogram 14. find_vv
c     compute the VV dyad at the quadrature points.  this is almost
c     a duplicate of find_bb, but there is no need in linear
c     computations, and the dyad is not normalized.
c
c     this routine is only called for implicit advection
c     when the upwinding-like diffusivity is on.
c-----------------------------------------------------------------------
      SUBROUTINE find_vv
      USE local
      USE pardata
      USE mpi_nim
      USE fields
      USE global
      USE input
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,imode,ierror,iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmp
      REAL(r8) :: kf
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     work directly with the data stored at quadrature points.
c     n=0 contribution has equilibrium field added, and n/=0
c     contribution results from +-|n|.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          rb(ibl)%qvv%qpf=0
          DO imode=1,nmodes
            kf=1._r8
            IF (keff(imode)==0.AND.eq_flow/='none') THEN
              rb(ibl)%qvv%qpf(1:3,:,:)=rb(ibl)%qvv%qpf(1:3,:,:)+
     $          (rb(ibl)%qve%qpf(1:3,:,:,imode)+rb(ibl)%qve_eq%qpf)**2
              rb(ibl)%qvv%qpf(4,:,:)=rb(ibl)%qvv%qpf(4,:,:)+
     $          (rb(ibl)%qve%qpf(1,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(1,:,:))*
     $          (rb(ibl)%qve%qpf(2,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(2,:,:))
              rb(ibl)%qvv%qpf(5,:,:)=rb(ibl)%qvv%qpf(5,:,:)+
     $          (rb(ibl)%qve%qpf(1,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(1,:,:))*
     $          (rb(ibl)%qve%qpf(3,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(3,:,:))
              rb(ibl)%qvv%qpf(6,:,:)=rb(ibl)%qvv%qpf(6,:,:)+
     $          (rb(ibl)%qve%qpf(2,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(2,:,:))*
     $          (rb(ibl)%qve%qpf(3,:,:,imode)+
     $           rb(ibl)%qve_eq%qpf(3,:,:))
            ELSE
              IF (keff(imode)==0) kf=0.5_r8
              rb(ibl)%qvv%qpf(1:3,:,:)=rb(ibl)%qvv%qpf(1:3,:,:)+
     $           2*kf*rb(ibl)%qve%qpf(1:3,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(1:3,:,:,imode))
              rb(ibl)%qvv%qpf(4,:,:)=rb(ibl)%qvv%qpf(4,:,:)+
     $            kf*(rb(ibl)%qve%qpf(1,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(2,:,:,imode))+
     $                rb(ibl)%qve%qpf(2,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(1,:,:,imode)))
              rb(ibl)%qvv%qpf(5,:,:)=rb(ibl)%qvv%qpf(5,:,:)+
     $            kf*(rb(ibl)%qve%qpf(1,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(3,:,:,imode))+
     $                rb(ibl)%qve%qpf(3,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(1,:,:,imode)))
              rb(ibl)%qvv%qpf(6,:,:)=rb(ibl)%qvv%qpf(6,:,:)+
     $            kf*(rb(ibl)%qve%qpf(2,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(3,:,:,imode))+
     $                rb(ibl)%qve%qpf(3,:,:,imode)*
     $          CONJG(rb(ibl)%qve%qpf(2,:,:,imode)))
            ENDIF
          ENDDO
          IF (nlayers>1) THEN
            ALLOCATE(tmp(6,rb(ibl)%ng,rb(ibl)%mx*rb(ibl)%my))
            CALL mpi_allreduce(rb(ibl)%qvv%qpf(1,1,1),tmp,SIZE(tmp),
     $                         mpi_nim_real,mpi_sum,comm_mode,ierror)
            rb(ibl)%qvv%qpf=tmp
            DEALLOCATE(tmp)
          ENDIF
c-----------------------------------------------------------------------
c       same for tblocks.
c-----------------------------------------------------------------------
        ELSE
          tb(ibl)%qvv%qpf=0
          DO imode=1,nmodes
            kf=1._r8
            IF (keff(imode)==0.AND.eq_flow/='none') THEN
              tb(ibl)%qvv%qpf(1:3,:,:)=tb(ibl)%qvv%qpf(1:3,:,:)+
     $          (tb(ibl)%qve%qpf(1:3,:,:,imode)+tb(ibl)%qve_eq%qpf)**2
              tb(ibl)%qvv%qpf(4,:,:)=tb(ibl)%qvv%qpf(4,:,:)+
     $          (tb(ibl)%qve%qpf(1,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(1,:,:))*
     $          (tb(ibl)%qve%qpf(2,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(2,:,:))
              tb(ibl)%qvv%qpf(5,:,:)=tb(ibl)%qvv%qpf(5,:,:)+
     $          (tb(ibl)%qve%qpf(1,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(1,:,:))*
     $          (tb(ibl)%qve%qpf(3,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(3,:,:))
              tb(ibl)%qvv%qpf(6,:,:)=tb(ibl)%qvv%qpf(6,:,:)+
     $          (tb(ibl)%qve%qpf(2,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(2,:,:))*
     $          (tb(ibl)%qve%qpf(3,:,:,imode)+
     $           tb(ibl)%qve_eq%qpf(3,:,:))
            ELSE
              IF (keff(imode)==0) kf=0.5_r8
              tb(ibl)%qvv%qpf(1:3,:,:)=tb(ibl)%qvv%qpf(1:3,:,:)+
     $           2*kf*tb(ibl)%qve%qpf(1:3,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(1:3,:,:,imode))
              tb(ibl)%qvv%qpf(4,:,:)=tb(ibl)%qvv%qpf(4,:,:)+
     $            kf*(tb(ibl)%qve%qpf(1,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(2,:,:,imode))+
     $                tb(ibl)%qve%qpf(2,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(1,:,:,imode)))
              tb(ibl)%qvv%qpf(5,:,:)=tb(ibl)%qvv%qpf(5,:,:)+
     $            kf*(tb(ibl)%qve%qpf(1,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(3,:,:,imode))+
     $                tb(ibl)%qve%qpf(3,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(1,:,:,imode)))
              tb(ibl)%qvv%qpf(6,:,:)=tb(ibl)%qvv%qpf(6,:,:)+
     $            kf*(tb(ibl)%qve%qpf(2,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(3,:,:,imode))+
     $                tb(ibl)%qve%qpf(3,:,:,imode)*
     $          CONJG(tb(ibl)%qve%qpf(2,:,:,imode)))
            ENDIF
          ENDDO
          IF (nlayers>1) THEN
            ALLOCATE(tmp(6,tb(ibl)%ng,tb(ibl)%mcell))
            CALL mpi_allreduce(tb(ibl)%qvv%qpf(1,1,1),tmp,SIZE(tmp),
     $                         mpi_nim_real,mpi_sum,comm_mode,ierror)
            tb(ibl)%qvv%qpf=tmp
            DEALLOCATE(tmp)
          ENDIF
        ENDIF
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE find_vv
c-----------------------------------------------------------------------
c     subprogram 15. b_store
c     save magnetic field at the quadrature points.  the saved data
c     will either be from the end of the step or the average
c     of the beginning and ending values, as determined by the passed
c     parameter.
c-----------------------------------------------------------------------
      SUBROUTINE b_store(bchoice)
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE math_tran
      USE physdat
      USE time
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: bchoice

      INTEGER(i4) :: ibl,mxb,myb,ng,ig
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     interface block for qp0_bcast.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE qp0_bcast(threedata,twodata,nl)
        USE local
        COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
        REAL(r8), DIMENSION(:,:,:) :: twodata
        INTEGER(i4), INTENT(IN) :: nl
        END SUBROUTINE qp0_bcast
      END INTERFACE
c-----------------------------------------------------------------------
c     update quadrature-point storage.
c
c     use the work1 temporary space to find the average.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      SELECT CASE(bchoice)
      CASE('ave')
        DO ibl=1,nbl
          work1(ibl)=be_old(ibl)
          CALL vector_add(work1(ibl),be(ibl),v1fac=0.5_r8,v2fac=0.5_r8)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qbe,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qbe,tb(ibl))
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c     store the field at the end of the step, including the separate
c     storage for the symmetric part.
c-----------------------------------------------------------------------
      CASE('end')
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%be,rb(ibl)%qbe,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%be,tb(ibl)%qbe,tb(ibl))
        ENDDO
      END SELECT
c-----------------------------------------------------------------------
c     qp operations that do not depend on bchoice.
c     [however, the current density broadcast is here to avoid
c     recomputing it below.]
c
c     for the nonlinear Hall iteration, update qbe_tot and qja_tot
c     with qbe + 0.5*qwork1, where qwork1 has the most recent
c     iterate of delta-B.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        IF (nonlinear) THEN
          IF (bchoice=='hall it') THEN
            CALL qp_fft_save(rb(ibl)%qbe%qpf+0.5_r8*rb(ibl)%qwork1%qpf,
     $                       rb(ibl)%qbe_tot%qpf,
     $                       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                       3_i4,rb(ibl)%ng,rb(ibl)%qbe_eq%qpf)
          ELSE
            CALL qp_fft_save(rb(ibl)%qbe%qpf,rb(ibl)%qbe_tot%qpf,
     $                       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                       3_i4,rb(ibl)%ng,rb(ibl)%qbe_eq%qpf)
          ENDIF
        ENDIF
        IF (nonlinear.AND.(impladv.OR.eta_model=="chodura".OR.
     $                     siop_type=='3D')) THEN
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          ng=rb(ibl)%ng
          ALLOCATE(ctmp(3,ng,mxb*myb,nmodes))
          IF (bchoice=='hall it') THEN
            CALL math_curl(nmodes,keff,geom,rb(ibl)%bigr,
     $                     rb(ibl)%qbe%qpf +0.5_r8*rb(ibl)%qwork1%qpf,
     $                     rb(ibl)%qbe%qpfr+0.5_r8*rb(ibl)%qwork1%qpfr,
     $                     rb(ibl)%qbe%qpfz+0.5_r8*rb(ibl)%qwork1%qpfz,
     $                     ctmp,1._r8/mu0)
          ELSE
            CALL math_curl(nmodes,keff,geom,rb(ibl)%bigr,
     $                     rb(ibl)%qbe%qpf,rb(ibl)%qbe%qpfr,
     $                     rb(ibl)%qbe%qpfz,ctmp,1._r8/mu0)
          ENDIF
          CALL qp_fft_save(ctmp,rb(ibl)%qja_tot%qpf,mxb,myb,
     $                     mpsq_block(ibl),3_i4,ng,rb(ibl)%qja_eq%qpf)
          DEALLOCATE(ctmp)
        ENDIF
      ENDDO
      DO ibl=nrbl+1,nbl
        IF (nonlinear) THEN
          IF (bchoice=='hall it') THEN
            CALL qp_fft_save(tb(ibl)%qbe%qpf+0.5_r8*tb(ibl)%qwork1%qpf,
     $                       tb(ibl)%qbe_tot%qpf,
     $                       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                       3_i4,tb(ibl)%ng,tb(ibl)%qbe_eq%qpf)
          ELSE
            CALL qp_fft_save(tb(ibl)%qbe%qpf,tb(ibl)%qbe_tot%qpf,
     $                       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                       3_i4,tb(ibl)%ng,tb(ibl)%qbe_eq%qpf)
          ENDIF
        ENDIF
        IF (nonlinear.AND.(impladv.OR.eta_model=="chodura".OR.
     $                     siop_type=='3D')) THEN
          mxb=tb(ibl)%mcell
          myb=1_i4
          ng=tb(ibl)%ng
          ALLOCATE(ctmp(3,ng,mxb*myb,nmodes))
          IF (bchoice=='hall it') THEN
            CALL math_curl(nmodes,keff,geom,
     $                     tb(ibl)%tgeom%bigr,
     $                     tb(ibl)%qbe%qpf +0.5_r8*tb(ibl)%qwork1%qpf,
     $                     tb(ibl)%qbe%qpfr+0.5_r8*tb(ibl)%qwork1%qpfr,
     $                     tb(ibl)%qbe%qpfz+0.5_r8*tb(ibl)%qwork1%qpfz,
     $                     ctmp,1._r8/mu0)
          ELSE
            CALL math_curl(nmodes,keff,geom,tb(ibl)%tgeom%bigr,
     $                     tb(ibl)%qbe%qpf,tb(ibl)%qbe%qpfr,
     $                     tb(ibl)%qbe%qpfz,ctmp,1._r8/mu0)
          ENDIF
          CALL qp_fft_save(ctmp,tb(ibl)%qja_tot%qpf,mxb,myb,
     $                     mpsq_block(ibl),3_i4,ng,tb(ibl)%qja_eq%qpf)
          DEALLOCATE(ctmp)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     broadcast the symmetric component here--derivatives are needed.
c     also broadcast asymmetric components that are used for coupling
c     in the preconditioner.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.bchoice/='hall it') THEN
        DO ibl=1,nrbl
          CALL qp0_bcast(rb(ibl)%qbe%qpf,rb(ibl)%qbe_n0%qpf,nlayers)
          CALL qp0_bcast(rb(ibl)%qbe%qpfr,rb(ibl)%qbe_n0%qpfr,nlayers)
          CALL qp0_bcast(rb(ibl)%qbe%qpfz,rb(ibl)%qbe_n0%qpfz,nlayers)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL qp0_bcast(tb(ibl)%qbe%qpf,tb(ibl)%qbe_n0%qpf,nlayers)
          CALL qp0_bcast(tb(ibl)%qbe%qpfr,tb(ibl)%qbe_n0%qpfr,nlayers)
          CALL qp0_bcast(tb(ibl)%qbe%qpfz,tb(ibl)%qbe_n0%qpfz,nlayers)
        ENDDO
      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE b_store
c-----------------------------------------------------------------------
c     subprogram 17. pieq_comp
c     find the equilibrium contributions to the center-of-mass stress
c     tensor.  the results is left in its separate physical
c     contributions.
c
c     nonlinear computations include factors of equilibrium number
c     density, but linear computations save the coefficient of the
c     perturbed linear n that depends on V_eq and other equilibrium
c     fields.
c
c     this routine is only called during start-up, so readability can be
c     favored over optimization.
c-----------------------------------------------------------------------
      SUBROUTINE pieq_comp
      USE local
      USE fields
      USE input
      USE physdat
      USE global
      USE time
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          grdveq,beq,ndq,diff,pi_veq,pareq,gyreq,tiq
      REAL(r8), DIMENSION(3,3) :: vtmpr,ptmpr,btmp,wdotr,bc_wdotr
      REAL(r8) :: rscal,divtmpr,rb2,vcoef,gcoef,kpmin,kpmax
      REAL(r8) :: timest_qp,timend_qp
      REAL(r8), PARAMETER :: third=1._r8/3._r8,twot=2._r8/3._r8
      INTEGER(i4) :: ibl,ncx,ncy,ix,iy,ng,ig,i2
c-----------------------------------------------------------------------
c     return if there is no equilibrium flow.
c-----------------------------------------------------------------------
      IF (eq_flow=='none') RETURN
      CALL timer(timest_qp)
c-----------------------------------------------------------------------
c     loop over all blocks and start by assigning pointers to the
c     appropriate data.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          ncx=rb(ibl)%mx
          ncy=rb(ibl)%my
          ng=rb(ibl)%ng
          ndq=>rb(ibl)%qnd_eq%qpf
          beq=>rb(ibl)%qbe_eq%qpf
          tiq=>rb(ibl)%qtion_eq%qpf
          grdveq=>rb(ibl)%qgrdveq%qpf
          pi_veq=>rb(ibl)%qpi_veq%qpf
          IF (par_visc>0) pareq=>rb(ibl)%qpi_pareq%qpf
          IF (gyr_visc>0) gyreq=>rb(ibl)%qpi_gyreq%qpf
          IF (ds_use=='both'.OR.ds_use=='kin_visc')
     $      diff=>rb(ibl)%qdiff_shape%qpf
        ELSE
          ncx=tb(ibl)%mcell
          ncy=1
          ng=tb(ibl)%ng
          ndq=>tb(ibl)%qnd_eq%qpf
          beq=>tb(ibl)%qbe_eq%qpf
          tiq=>tb(ibl)%qtion_eq%qpf
          grdveq=>tb(ibl)%qgrdveq%qpf
          pi_veq=>tb(ibl)%qpi_veq%qpf
          IF (par_visc>0) pareq=>tb(ibl)%qpi_pareq%qpf
          IF (gyr_visc>0) gyreq=>tb(ibl)%qpi_gyreq%qpf
          IF (ds_use=='both'.OR.ds_use=='kin_visc')
     $     diff=>tb(ibl)%qdiff_shape%qpf
        ENDIF
        IF (.NOT.(ds_use=='both'.OR.ds_use=='kin_visc')) THEN
          ALLOCATE(diff(1,ng,ncx*ncy))
          diff=1._r8
        ENDIF
c-----------------------------------------------------------------------
c       kinematic part:
c-----------------------------------------------------------------------
        IF (kin_visc>0) THEN
          DO ix=1,ncx*ncy
            DO ig=1,ng
              rscal=-mtot*kin_visc*diff(1,ig,ix)
              pi_veq(:,ig,ix)=rscal*grdveq(:,ig,ix)
            ENDDO
          ENDDO
        ELSE
          pi_veq=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       isotropic part:
c-----------------------------------------------------------------------
        IF (iso_visc>0) THEN
          DO ix=1,ncx*ncy
            DO ig=1,ng
              rscal=mtot*iso_visc*diff(1,ig,ix)
              vtmpr(1:3,1)=grdveq(1:3,ig,ix)
              vtmpr(1:3,2)=grdveq(4:6,ig,ix)
              vtmpr(1:3,3)=grdveq(7:9,ig,ix)
              ptmpr=rscal*(vtmpr+TRANSPOSE(vtmpr))
              divtmpr=-twot*rscal*(vtmpr(1,1)+vtmpr(2,2)+vtmpr(3,3))
              ptmpr(1,1)=ptmpr(1,1)+divtmpr
              ptmpr(2,2)=ptmpr(2,2)+divtmpr
              ptmpr(3,3)=ptmpr(3,3)+divtmpr
              pi_veq(1:3,ig,ix)=pi_veq(1:3,ig,ix)-ptmpr(1:3,1)
              pi_veq(4:6,ig,ix)=pi_veq(4:6,ig,ix)-ptmpr(1:3,2)
              pi_veq(7:9,ig,ix)=pi_veq(7:9,ig,ix)-ptmpr(1:3,3)
            ENDDO
          ENDDO
        ENDIF
        IF (.NOT.(ds_use=='both'.OR.ds_use=='kin_visc'))
     $    DEALLOCATE(diff)
c-----------------------------------------------------------------------
c       multiply by nd_eq for nonlinear computations.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          DO ix=1,ncx*ncy
            DO ig=1,ng
              pi_veq(:,ig,ix)=pi_veq(:,ig,ix)*ndq(1,ig,ix)
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       parallel part:
c-----------------------------------------------------------------------
        IF (par_visc>0) THEN
          vcoef=-3._r8*mtot*par_visc
          DO ix=1,ncx*ncy
            DO ig=1,ng
              rb2=SUM(beq(:,ig,ix)**2)+smallnum
              divtmpr=-third*rb2*(grdveq(1,ig,ix)+
     $                            grdveq(5,ig,ix)+
     $                            grdveq(9,ig,ix))
              rscal=vcoef*(divtmpr+beq(1,ig,ix)*(
     $                beq(1,ig,ix)*grdveq(1,ig,ix)+
     $                beq(2,ig,ix)*grdveq(2,ig,ix)+
     $                beq(3,ig,ix)*grdveq(3,ig,ix))+
     $              beq(2,ig,ix)*(
     $                beq(1,ig,ix)*grdveq(4,ig,ix)+
     $                beq(2,ig,ix)*grdveq(5,ig,ix)+
     $                beq(3,ig,ix)*grdveq(6,ig,ix))+
     $              beq(3,ig,ix)*(
     $                beq(1,ig,ix)*grdveq(7,ig,ix)+
     $                beq(2,ig,ix)*grdveq(8,ig,ix)+
     $                beq(3,ig,ix)*grdveq(9,ig,ix)))/rb2**2
              rb2=-third*rb2
              vtmpr(1,1)=rscal*(beq(1,ig,ix)*beq(1,ig,ix)+rb2)
              vtmpr(2,1)=rscal* beq(1,ig,ix)*beq(2,ig,ix)
              vtmpr(3,1)=rscal* beq(1,ig,ix)*beq(3,ig,ix)
              vtmpr(1,2)=vtmpr(2,1)
              vtmpr(2,2)=rscal*(beq(2,ig,ix)*beq(2,ig,ix)+rb2)
              vtmpr(3,2)=rscal* beq(2,ig,ix)*beq(3,ig,ix)
              vtmpr(1,3)=vtmpr(3,1)
              vtmpr(2,3)=vtmpr(3,2)
              vtmpr(3,3)=rscal*(beq(3,ig,ix)*beq(3,ig,ix)+rb2)
              pareq(:,ig,ix)=RESHAPE(vtmpr,(/9/))
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         multiply by nd_eq for nonlinear computations.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            DO ix=1,ncx*ncy
              DO ig=1,ng
                pareq(:,ig,ix)=pareq(:,ig,ix)*ndq(1,ig,ix)
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           multiply by the temperature-dependence for parallel
c           transport coefficients, if specified.
c-----------------------------------------------------------------------
            IF (parvisc_model=='plltdep') THEN
              kpmin=k_pll_min/k_plli
              kpmax=k_pll_max/k_plli
              DO ix=1,ncx*ncy
                DO ig=1,ng
                  pareq(:,ig,ix)=pareq(:,ig,ix)*
     $               MAX(kpmin, MIN(kpmax, 
     $                  (MAX(smallnum,tiq(1,ig,ix))/k_pll_ref_t)**2.5))
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       gyroviscous part:
c-----------------------------------------------------------------------
        IF (gyr_visc>0) THEN
          vcoef=0.25*gyr_visc*ms(2)*kboltz/(zeff**2*elementary_q)
          DO ix=1,ncx*ncy
            DO ig=1,ng
              rb2=SUM(beq(:,ig,ix)**2)+smallnum
              gcoef=vcoef*tiq(1,ig,ix)/rb2**2
              vtmpr(1:3,1)=grdveq(1:3,ig,ix)+grdveq(1:7:3,ig,ix)
              vtmpr(1:3,2)=grdveq(4:6,ig,ix)+grdveq(2:8:3,ig,ix)
              vtmpr(1:3,3)=grdveq(7:9,ig,ix)+grdveq(3:9:3,ig,ix)
              btmp(1:3,1)=3._r8*beq(1:3,ig,ix)*beq(1,ig,ix)
              btmp(1:3,2)=3._r8*beq(1:3,ig,ix)*beq(2,ig,ix)
              btmp(1:3,3)=3._r8*beq(1:3,ig,ix)*beq(3,ig,ix)
              btmp(1,1)=rb2+btmp(1,1)
              btmp(2,2)=rb2+btmp(2,2)
              btmp(3,3)=rb2+btmp(3,3)
              DO i2=1,3
                wdotr(1,i2)=SUM(vtmpr(1,:)*btmp(:,i2))
                wdotr(2,i2)=SUM(vtmpr(2,:)*btmp(:,i2))
                wdotr(3,i2)=SUM(vtmpr(3,:)*btmp(:,i2))
              ENDDO
              bc_wdotr(1,:)=beq(2,ig,ix)*wdotr(3,:)-
     $                      beq(3,ig,ix)*wdotr(2,:)
              bc_wdotr(2,:)=beq(3,ig,ix)*wdotr(1,:)-
     $                      beq(1,ig,ix)*wdotr(3,:)
              bc_wdotr(3,:)=beq(1,ig,ix)*wdotr(2,:)-
     $                      beq(2,ig,ix)*wdotr(1,:)
              gyreq(1:3,ig,ix)=gcoef*(bc_wdotr(:,1)+bc_wdotr(1,:))
              gyreq(4:6,ig,ix)=gcoef*(bc_wdotr(:,2)+bc_wdotr(2,:))
              gyreq(7:9,ig,ix)=gcoef*(bc_wdotr(:,3)+bc_wdotr(3,:))
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         multiply by nd_eq for nonlinear computations.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            DO ix=1,ncx*ncy
              DO ig=1,ng
                gyreq(:,ig,ix)=gyreq(:,ig,ix)*ndq(1,ig,ix)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pieq_comp
c-----------------------------------------------------------------------
c     subprogram 18. iso_stress
c     assemble the sum of the kinematic and isotropic stresses, except
c     for any number density dependence.  the kinematic part is
c     proportional to grad(V), and the isotropic part is proportional to
c     the rate of strain tensor.  both are multiplied by a scalar 2D
c     field, and the computations are for for one grid-block of data.
c
c     all information is passed in the argument list for flexibility.
c     arrays are passed f77-style to avoid interface blocks.
c-----------------------------------------------------------------------
      SUBROUTINE iso_stress(ncx,ncy,nm,kvsc,ivsc,grdv,scal,piten)
      USE local
      USE time
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: ncx,ncy,nm
      REAL(r8), INTENT(IN) :: kvsc,ivsc
      REAL(r8), DIMENSION(ncx,ncy), INTENT(IN) :: scal
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(IN) :: grdv
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(OUT) :: piten

      INTEGER(i4) :: imode,ix,iy
      COMPLEX(r8), DIMENSION(3,3) :: vtmp,ptmp
      COMPLEX(r8) :: divtmp
      REAL(r8), PARAMETER :: twot=2._r8/3._r8
      REAL(r8) :: dvc,sumc,ic
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     loop over Fourier components and elements, and find
c       scal*( kvsc*grad(V)+ivsc*(grad(V)+grad(V)^T - 2/3*div(V)) )
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      ic=MAX(ivsc,0._r8)
      dvc=twot*ic
      sumc=ic+MAX(kvsc,0._r8)
      DO imode=1,nm
        DO iy=1,ncy
          DO ix=1,ncx
            vtmp(1:3,1)=grdv(1:3,ix,iy,imode)
            vtmp(1:3,2)=grdv(4:6,ix,iy,imode)
            vtmp(1:3,3)=grdv(7:9,ix,iy,imode)
            ptmp=scal(ix,iy)*(sumc*vtmp+ic*TRANSPOSE(vtmp))
            divtmp=-dvc*scal(ix,iy)*(vtmp(1,1)+vtmp(2,2)+vtmp(3,3))
            ptmp(1,1)=ptmp(1,1)+divtmp
            ptmp(2,2)=ptmp(2,2)+divtmp
            ptmp(3,3)=ptmp(3,3)+divtmp
            piten(1:3,ix,iy,imode)=-ptmp(1:3,1)
            piten(4:6,ix,iy,imode)=-ptmp(1:3,2)
            piten(7:9,ix,iy,imode)=-ptmp(1:3,3)
          ENDDO
        ENDDO
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iso_stress
c-----------------------------------------------------------------------
c     subprogram 19. par_stress_lin
c     find the linear parallel stress except for V_eq contributions and
c     the number density factor.  the passed scalar coefficient vcoef
c     is usually -3*dt*mtot*par_visc.
c
c     the computations are for for one grid-block of data.
c
c     all information is passed in the argument list for flexibility.
c     arrays are passed f77-style to avoid interface blocks.
c-----------------------------------------------------------------------
      SUBROUTINE par_stress_lin(ncx,ncy,nm,vcoef,sml,grdv,be_eq,parten)
      USE local
      USE time
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: ncx,ncy,nm
      REAL(r8), INTENT(IN) :: vcoef,sml
      REAL(r8), DIMENSION(3,ncx,ncy), INTENT(IN) :: be_eq
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(IN) :: grdv
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(OUT) :: parten

      INTEGER(i4) :: imode,ix,iy
      COMPLEX(r8) :: divtmp,cscal
      REAL(r8), PARAMETER :: third=1._r8/3._r8
      REAL(r8) :: rb2
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     find (Beq.grad(V).Beq - Beq**2*div(V))/Beq**4 and multiply
c     the result by (BeqBeq - I*Beq**2/3), where V is the perturbed
c     velocity.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO imode=1,nm
        DO iy=1,ncy
          DO ix=1,ncx
            rb2=SUM(be_eq(:,ix,iy)**2)+sml
            divtmp=-third*rb2*(grdv(1,ix,iy,imode)+
     $              grdv(5,ix,iy,imode)+grdv(9,ix,iy,imode))
            cscal=vcoef*(divtmp+be_eq(1,ix,iy)*(
     $              be_eq(1,ix,iy)*grdv(1,ix,iy,imode)+
     $              be_eq(2,ix,iy)*grdv(2,ix,iy,imode)+
     $              be_eq(3,ix,iy)*grdv(3,ix,iy,imode))+
     $            be_eq(2,ix,iy)*(
     $              be_eq(1,ix,iy)*grdv(4,ix,iy,imode)+
     $              be_eq(2,ix,iy)*grdv(5,ix,iy,imode)+
     $              be_eq(3,ix,iy)*grdv(6,ix,iy,imode))+
     $            be_eq(3,ix,iy)*(
     $              be_eq(1,ix,iy)*grdv(7,ix,iy,imode)+
     $              be_eq(2,ix,iy)*grdv(8,ix,iy,imode)+
     $              be_eq(3,ix,iy)*grdv(9,ix,iy,imode)))/rb2**2
            rb2=-third*rb2
            parten(1,ix,iy,imode)=
     $        cscal*(be_eq(1,ix,iy)*be_eq(1,ix,iy)+rb2)
            parten(2,ix,iy,imode)=
     $        cscal*be_eq(1,ix,iy)*be_eq(2,ix,iy)
            parten(3,ix,iy,imode)=
     $        cscal*be_eq(1,ix,iy)*be_eq(3,ix,iy)
            parten(4,ix,iy,imode)=parten(2,ix,iy,imode)
            parten(5,ix,iy,imode)=
     $        cscal*(be_eq(2,ix,iy)*be_eq(2,ix,iy)+rb2)
            parten(6,ix,iy,imode)=
     $        cscal*be_eq(2,ix,iy)*be_eq(3,ix,iy)
            parten(7,ix,iy,imode)=parten(3,ix,iy,imode)
            parten(8,ix,iy,imode)=parten(6,ix,iy,imode)
            parten(9,ix,iy,imode)=
     $        cscal*(be_eq(3,ix,iy)*be_eq(3,ix,iy)+rb2)
          ENDDO
        ENDDO
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE par_stress_lin
c-----------------------------------------------------------------------
c     subprogram 20. par_stress_veq
c     find the Veq contributions to the linear parallel stress, except
c     the number density factor.  the passed scalar coefficient vcoef
c     is usually -3*dt*mtot*par_visc.
c
c     the computations are for for one grid-block of data, and results
c     are added to an existing stress array.
c
c     all information is passed in the argument list for flexibility.
c     arrays are passed f77-style to avoid interface blocks.  scal
c     is space for a 2D work array.
c-----------------------------------------------------------------------
      SUBROUTINE par_stress_veq(ncx,ncy,nm,vcoef,dt,sml,grdv,be,be_eq,
     $                          grdveq,pareq,scal,parten)
      USE local
      USE time
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: ncx,ncy,nm
      REAL(r8), INTENT(IN) :: vcoef,dt,sml
      REAL(r8), DIMENSION(ncx,ncy) :: scal
      REAL(r8), DIMENSION(3,ncx,ncy), INTENT(IN) :: be_eq
      REAL(r8), DIMENSION(9,ncx,ncy), INTENT(IN) :: grdveq,pareq
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(IN) :: grdv
      COMPLEX(r8), DIMENSION(3,ncx,ncy,nm), INTENT(IN) :: be
      COMPLEX(r8), DIMENSION(9,ncx,ncy,nm), INTENT(OUT) :: parten

      INTEGER(i4) :: imode,ix,iy
      COMPLEX(r8) :: divtmp,cscal,eqfac
      REAL(r8), DIMENSION(3,3) :: vtmpr
      REAL(r8), PARAMETER :: third=1._r8/3._r8,twot=2._r8/3._r8
      REAL(r8) :: rb2,divtmpr,rscal
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     first find the coefficient that is used for all Fourier
c     components.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO iy=1,ncy
        DO ix=1,ncx
          rb2=SUM(be_eq(:,ix,iy)**2)+sml
          divtmpr=-third*rb2*(grdveq(1,ix,iy)+grdveq(5,ix,iy)+
     $                        grdveq(9,ix,iy))
          scal(ix,iy)=vcoef*(divtmpr+be_eq(1,ix,iy)*(
     $            be_eq(1,ix,iy)*grdveq(1,ix,iy)+
     $            be_eq(2,ix,iy)*grdveq(2,ix,iy)+
     $            be_eq(3,ix,iy)*grdveq(3,ix,iy))+
     $          be_eq(2,ix,iy)*(
     $            be_eq(1,ix,iy)*grdveq(4,ix,iy)+
     $            be_eq(2,ix,iy)*grdveq(5,ix,iy)+
     $            be_eq(3,ix,iy)*grdveq(6,ix,iy))+
     $          be_eq(3,ix,iy)*(
     $            be_eq(1,ix,iy)*grdveq(7,ix,iy)+
     $            be_eq(2,ix,iy)*grdveq(8,ix,iy)+
     $            be_eq(3,ix,iy)*grdveq(9,ix,iy)))/rb2**2
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     V_eq contributions to the linear parallel stress in three
c     pieces: 1) scalar product terms using the symmetry of Weq,
c     2) dyad terms, and 3) magnitude terms.  density
c     perturbations are added separately if continuity is full.
c-----------------------------------------------------------------------
      DO imode=1,nm
        DO iy=1,ncy
          DO ix=1,ncx
            rb2=SUM(be_eq(:,ix,iy)**2)+sml
            divtmpr=-twot*
     $        (grdveq(1,ix,iy)+grdveq(5,ix,iy)+grdveq(9,ix,iy))
            vtmpr(1:3,1)=grdveq(1:3,ix,iy)+grdveq(1:7:3,ix,iy)
            vtmpr(1:3,2)=grdveq(4:6,ix,iy)+grdveq(2:8:3,ix,iy)
            vtmpr(1:3,3)=grdveq(7:9,ix,iy)+grdveq(3:9:3,ix,iy)
            vtmpr(1,1)=vtmpr(1,1)+divtmpr
            vtmpr(2,2)=vtmpr(2,2)+divtmpr
            vtmpr(3,3)=vtmpr(3,3)+divtmpr
            cscal=vcoef*(
     $        be_eq(1,ix,iy)*SUM(be(:,ix,iy,imode)*vtmpr(:,1))+
     $        be_eq(2,ix,iy)*SUM(be(:,ix,iy,imode)*vtmpr(:,2))+
     $        be_eq(3,ix,iy)*SUM(be(:,ix,iy,imode)*vtmpr(:,3)))/rb2**2
            rscal=scal(ix,iy)
            divtmp=SUM(be_eq(:,ix,iy)*be(:,ix,iy,imode))
            eqfac=4._r8*divtmp*dt/rb2

            parten(1,ix,iy,imode)=parten(1,ix,iy,imode)+
     $        cscal*(be_eq(1,ix,iy)*be_eq(1,ix,iy)-third*rb2)+
     $        rscal*2*(be_eq(1,ix,iy)*be(1,ix,iy,imode)-
     $                 third*divtmp)-
     $        eqfac*pareq(1,ix,iy)
            parten(2,ix,iy,imode)=parten(2,ix,iy,imode)+
     $        cscal*be_eq(1,ix,iy)*be_eq(2,ix,iy)+
     $        rscal*(be_eq(2,ix,iy)*be(1,ix,iy,imode)+
     $               be_eq(1,ix,iy)*be(2,ix,iy,imode))-
     $        eqfac*pareq(2,ix,iy)
            parten(3,ix,iy,imode)=parten(3,ix,iy,imode)+
     $        cscal*be_eq(1,ix,iy)*be_eq(3,ix,iy)+
     $        rscal*(be_eq(3,ix,iy)*be(1,ix,iy,imode)+
     $               be_eq(1,ix,iy)*be(3,ix,iy,imode))-
     $        eqfac*pareq(3,ix,iy)
            parten(4,ix,iy,imode)=parten(2,ix,iy,imode)
            parten(5,ix,iy,imode)=parten(5,ix,iy,imode)+
     $        cscal*(be_eq(2,ix,iy)*be_eq(2,ix,iy)-third*rb2)+
     $        rscal*2*(be_eq(2,ix,iy)*be(2,ix,iy,imode)-
     $                 third*divtmp)-
     $        eqfac*pareq(5,ix,iy)
            parten(6,ix,iy,imode)=parten(6,ix,iy,imode)+
     $        cscal*be_eq(2,ix,iy)*be_eq(3,ix,iy)+
     $        rscal*(be_eq(3,ix,iy)*be(2,ix,iy,imode)+
     $               be_eq(2,ix,iy)*be(3,ix,iy,imode))-
     $        eqfac*pareq(6,ix,iy)
            parten(7,ix,iy,imode)=parten(3,ix,iy,imode)
            parten(8,ix,iy,imode)=parten(6,ix,iy,imode)
            parten(9,ix,iy,imode)=parten(9,ix,iy,imode)+
     $        cscal*(be_eq(3,ix,iy)*be_eq(3,ix,iy)-third*rb2)+
     $        rscal*2*(be_eq(3,ix,iy)*be(3,ix,iy,imode)-
     $                 third*divtmp)-
     $        eqfac*pareq(9,ix,iy)
          ENDDO
        ENDDO
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE par_stress_veq
c-----------------------------------------------------------------------
c     subprogram 21. par_stress_nl
c     find the nonlinear parallel stress, except the number density
c     factor, as a function of toroidal angle.  the passed coefficient
c     array vcoef is -3*dt*mtot*par_visc multiplied by a spatial
c     profile.
c
c     the computations are for for one grid-block of data, and results
c     are added to an existing stress array.
c
c     all information is passed in the argument list for flexibility.
c     arrays are passed f77-style to avoid interface blocks.
c-----------------------------------------------------------------------
      SUBROUTINE par_stress_nl(npol,nphi,sml,vcoef,real_grdv,real_be,
     $                         real_pten)
      USE local
      USE time
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: npol,nphi
      REAL(r8), INTENT(IN) :: sml
      REAL(r8), DIMENSION(npol,*), INTENT(IN) :: vcoef
      REAL(r8), DIMENSION(3,npol,nphi), INTENT(IN) :: real_be
      REAL(r8), DIMENSION(9,npol,nphi), INTENT(IN) :: real_grdv
      REAL(r8), DIMENSION(9,npol,nphi), INTENT(INOUT) :: real_pten

      INTEGER(i4) :: ip,ix
      REAL(r8), DIMENSION(3,3) :: btmp
      REAL(r8), PARAMETER :: third=1._r8/3._r8
      REAL(r8) :: rb2,divtmpr,rscal
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     find (B.grad(V).B - B**2*div(V))/B**4 and multiply the result
c     by (BB - I*B**2/3).
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO ip=1,nphi
        DO ix=1,npol
          rb2=SUM(real_be(:,ix,ip)**2)+sml
          divtmpr=-third*rb2*(real_grdv(1,ix,ip)+
     $             real_grdv(5,ix,ip)+real_grdv(9,ix,ip))
          rscal=vcoef(ix,ip)*(divtmpr+real_be(1,ix,ip)*(
     $            real_be(1,ix,ip)*real_grdv(1,ix,ip)+
     $            real_be(2,ix,ip)*real_grdv(2,ix,ip)+
     $            real_be(3,ix,ip)*real_grdv(3,ix,ip))+
     $          real_be(2,ix,ip)*(
     $            real_be(1,ix,ip)*real_grdv(4,ix,ip)+
     $            real_be(2,ix,ip)*real_grdv(5,ix,ip)+
     $            real_be(3,ix,ip)*real_grdv(6,ix,ip))+
     $          real_be(3,ix,ip)*(
     $            real_be(1,ix,ip)*real_grdv(7,ix,ip)+
     $            real_be(2,ix,ip)*real_grdv(8,ix,ip)+
     $            real_be(3,ix,ip)*real_grdv(9,ix,ip)))/rb2**2
          rb2=-third*rb2

          btmp(1,1)=rscal*real_be(1,ix,ip)*real_be(2,ix,ip)
          btmp(2,1)=rscal*real_be(1,ix,ip)*real_be(3,ix,ip)
          btmp(3,1)=rscal*real_be(2,ix,ip)*real_be(3,ix,ip)

          real_pten(1,ix,ip)=real_pten(1,ix,ip)+rscal*
     $      (real_be(1,ix,ip)*real_be(1,ix,ip)+rb2)
          real_pten(2,ix,ip)=real_pten(2,ix,ip)+btmp(1,1)
          real_pten(3,ix,ip)=real_pten(3,ix,ip)+btmp(2,1)
          real_pten(4,ix,ip)=real_pten(4,ix,ip)+btmp(1,1)
          real_pten(5,ix,ip)=real_pten(5,ix,ip)+rscal*
     $      (real_be(2,ix,ip)*real_be(2,ix,ip)+rb2)
          real_pten(6,ix,ip)=real_pten(6,ix,ip)+btmp(3,1)
          real_pten(7,ix,ip)=real_pten(7,ix,ip)+btmp(2,1)
          real_pten(8,ix,ip)=real_pten(8,ix,ip)+btmp(3,1)
          real_pten(9,ix,ip)=real_pten(9,ix,ip)+rscal*
     $      (real_be(3,ix,ip)*real_be(3,ix,ip)+rb2)
        ENDDO
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE par_stress_nl
c-----------------------------------------------------------------------
c     subprogram 22. viscous_heating
c     compute the viscous heating density using the present state
c     stored in the quadrature-point data.  the computations are
c     performed for all quadrature points in all grid blocks,
c     and the only quad-point data that should be modified is the
c     storage for viscous heating.
c
c     the stress computation necessarily mimics vrhs in integrands,
c     and the external stress computations help minimize duplicate
c     coding.
c
c     note that gyroviscosity does not contribute to heating, so the
c     gyroviscous stress is not computed here.
c-----------------------------------------------------------------------
      SUBROUTINE viscous_heating
      USE local
      USE physdat
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE math_tran
      USE pardata
      USE mpi_nim
      USE fft_mod
      IMPLICIT NONE

      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $             vpt,ve_r,ve_z,bpt,vheat,ndpt
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: grdv,piten,
     $             parten
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          beq,grdveq,pareq,pi_veq,diff,ndq,nd0
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $          real_bptr,real_grdv,real_ndptr,kappli
      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_pten,real_scal
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: dtmp

      INTEGER(i4) :: ibl,im,ip,ig,ncx,ncy,ng,mps
      REAL(r8) :: vcoef,dtmt
c-----------------------------------------------------------------------
c     find the data from the quadrature-point storage.
c-----------------------------------------------------------------------
      dtmt=dt*mtot
      vcoef=-3._r8*dtmt*par_visc
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          ncx=rb(ibl)%mx
          ncy=rb(ibl)%my
          ng=rb(ibl)%ng
        ELSE
          ncx=tb(ibl)%mcell
          ncy=1
          ng=tb(ibl)%ng
        ENDIF
        mps=mpsq_block(ibl)
c-----------------------------------------------------------------------
c       allocate space for this block.
c-----------------------------------------------------------------------
        ALLOCATE(grdv(9,ng,ncx*ncy,nmodes))
        ALLOCATE(piten(9,ng,ncx*ncy,nmodes))
        ALLOCATE(dtmp(ng,ncx*ncy))
        IF (nonlinear) THEN
          ALLOCATE(real_pten(9,mps,nphi))
          ALLOCATE(real_scal(1,mps,nphi))
        ENDIF
        IF (par_visc>0) ALLOCATE(parten(9,ng,ncx*ncy,nmodes))
c-----------------------------------------------------------------------
c       find the data from the quadrature-point storage.
c-----------------------------------------------------------------------
          IF (ibl<=nrbl) THEN
            bigr=>rb(ibl)%bigr
            vpt=>rb(ibl)%qve%qpf
            ve_r=>rb(ibl)%qve%qpfr
            ve_z=>rb(ibl)%qve%qpfz
            real_grdv=>rb(ibl)%qgrdv%qpf
            grdveq=>rb(ibl)%qgrdveq%qpf
            pi_veq=>rb(ibl)%qpi_veq%qpf
            pareq=>rb(ibl)%qpi_pareq%qpf
            bpt=>rb(ibl)%qbe%qpf
            real_bptr=>rb(ibl)%qbe_tot%qpf
            beq=>rb(ibl)%qbe_eq%qpf
            diff=>rb(ibl)%qdiff_shape%qpf
            ndpt=>rb(ibl)%qnd%qpf
            ndq=>rb(ibl)%qnd_eq%qpf
            nd0=>rb(ibl)%qnd_n0%qpf
            real_ndptr=>rb(ibl)%qnd_tot%qpf
            vheat=>rb(ibl)%qvisc%qpf
            IF (par_visc>0.AND.parvisc_model=='plltdep')
     $        kappli=>rb(ibl)%qkappli_phi%qpf
          ELSE
            bigr=>tb(ibl)%tgeom%bigr
            vpt=>tb(ibl)%qve%qpf
            ve_r=>tb(ibl)%qve%qpfr
            ve_z=>tb(ibl)%qve%qpfz
            real_grdv=>tb(ibl)%qgrdv%qpf
            grdveq=>tb(ibl)%qgrdveq%qpf
            pi_veq=>tb(ibl)%qpi_veq%qpf
            pareq=>tb(ibl)%qpi_pareq%qpf
            bpt=>tb(ibl)%qbe%qpf
            real_bptr=>tb(ibl)%qbe_tot%qpf
            beq=>tb(ibl)%qbe_eq%qpf
            diff=>tb(ibl)%qdiff_shape%qpf
            ndpt=>tb(ibl)%qnd%qpf
            ndq=>tb(ibl)%qnd_eq%qpf
            nd0=>tb(ibl)%qnd_n0%qpf
            real_ndptr=>tb(ibl)%qnd_tot%qpf
            vheat=>tb(ibl)%qvisc%qpf
            IF (par_visc>0.AND.parvisc_model=='plltdep')
     $        kappli=>tb(ibl)%qkappli_phi%qpf
          ENDIF
c-----------------------------------------------------------------------
c         start with the grad(V) tensor and add grad(V_eq).
c-----------------------------------------------------------------------
          CALL math_grad(nmodes,keff,3_i4,geom,vpt,ve_r,ve_z,grdv,bigr)
          IF (nonlinear.AND.eq_flow/='none') THEN
            DO im=1,nmodes
              IF (keff(im)==0) grdv(:,:,:,im)=grdv(:,:,:,im)+grdveq
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         stress tensor for kinetematic and isotropic viscosity--now
c         in an external subroutine.
c-----------------------------------------------------------------------
          IF (iso_visc>0.OR.kin_visc>0) THEN
            IF (ds_use=='kin_visc'.OR.ds_use=='both') THEN
              dtmp=dtmt*MAX(diff(1,:,:),0._r8)
            ELSE
              dtmp=dtmt
            ENDIF
            CALL iso_stress(ng,ncx*ncy,nmodes,kin_visc,iso_visc,
     $                      grdv,dtmp,piten)
          ELSE
            piten=0._r8
          ENDIF
c-----------------------------------------------------------------------
c         transform piten from kinematic and isotropic contributions.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            IF (iso_visc>0.OR.kin_visc>0) THEN
              CALL fft_nim('inverse',ng*ncx*ncy,mps,lphi,9_i4,
     $                     piten,real_pten,dealiase)
            ELSE
              real_pten=0._r8
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c         stress tensor for parallel viscosity
c
c         the steady-state term is subtracted after forming the 
c         complete stress. 
c
c         stored equilibrium stress does not include the factor of dt.
c-----------------------------------------------------------------------
          IF (par_visc>0) THEN
c-----------------------------------------------------------------------
c           nonlinear stress from the external subroutine par_stress_nl.
c
c           the Braginskii temperature dependence is the same as the
c           parallel ion thermal diffusivity, so use that data for
c           the coefficient if specified.
c-----------------------------------------------------------------------
            IF (nonlinear) THEN
              IF (parvisc_model=='plltdep') THEN
                real_scal=vcoef*kappli/k_plli
              ELSE
                real_scal=vcoef
              ENDIF
              CALL par_stress_nl(mps,nphi,smallnum,real_scal,real_grdv,
     $                           real_bptr,real_pten)
c-----------------------------------------------------------------------
c           linear version from the external subroutine par_stress_lin.
c-----------------------------------------------------------------------
            ELSE
              CALL par_stress_lin(ng,ncx*ncy,nmodes,vcoef,smallnum,
     $                            grdv,beq,parten)
c-----------------------------------------------------------------------
c             V_eq contributions to the linear parallel stress from
c             the external subroutine par_stress_veq.
c-PRE
c             Braginskii parallel does not have density.
c-----------------------------------------------------------------------
              IF (eq_flow/='none') THEN
                CALL par_stress_veq(ng,ncx*ncy,nmodes,vcoef,dt,smallnum,
     $                              grdv,bpt,beq,grdveq,pareq,dtmp,
     $                              parten)
              ENDIF
              piten=piten+parten
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c         now find the viscous heating density,
c
c         -dt*(grad(V))^T:Pi(V)
c
c         is formed here and saved as quadrature point data.
c         Pi is the stress tensor function of V which satisfies
c
c         (grad(A))^T:Pi(V) = (grad(V))^T:Pi(A)
c
c         where boundary conditions eliminate the surfaces terms
c         when integrating over the volume.
c
c         note that V and rho have to be divided into equilibrium and
c         perturbed parts.  for linear computations, writing Pi(v)
c         as everything but the density factor, we form
c
c         -r_eq*2*(grad(V_eq))^T:Pi(v) - r*(grad(V_eq))^T:Pi(V_eq)
c
c         for nonlinear computations, real_grdv and real_pten include
c         V_eq contributions, so the heating we want looks like
c
c         -(r_eq+r)*((grad(V_eq+v))^T:Pi(V_eq+v))
c                +r_eq*(grad(V_eq))^T:Pi(V_eq)
c
c         what perturbed density, r, is used depends on the continuity
c         input parameter.
c-----------------------------------------------------------------------
          IF (nonlinear) THEN
            real_scal(1,:,:)=-SUM(real_grdv*real_pten,1)
            IF (continuity=='full') real_scal=real_scal*real_ndptr
            CALL fft_nim('forward',ng*ncx*ncy,mps,lphi,1_i4,
     $                   vheat,real_scal,dealiase)
          ELSE
            IF (eq_flow/='none') THEN
              DO im=1,nmodes
                vheat(1,:,:,im)=-2._r8*SUM(grdveq*piten(:,:,:,im),1)
              ENDDO
            ELSE
              vheat=0
            ENDIF
          ENDIF
          IF (.NOT.nonlinear.OR.continuity=='fix profile') THEN
            DO im=1,nmodes
              vheat(:,:,:,im)=vheat(:,:,:,im)*ndq
            ENDDO
          ELSE IF (nonlinear.AND.continuity=='n=0 only') THEN
            DO im=1,nmodes
              vheat(:,:,:,im)=vheat(:,:,:,im)*(ndq+nd0)
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         now for the pesky -(grad(V_eq))^T:Pi(V_eq) term.  for 
c         nonlinear computations only, Pi(V_eq) now includes the
c         nd_eq factor.
c-----------------------------------------------------------------------
          IF (eq_flow/='none') THEN
            IF (nonlinear) THEN
              DO im=1,nmodes
                IF (keff(im)/=0) CYCLE
                IF (par_visc>0) THEN
                  vheat(1,:,:,im)=vheat(1,:,:,im)+
     $              dt*SUM(grdveq*(pi_veq+pareq),1)
                ELSE
                  vheat(1,:,:,im)=vheat(1,:,:,im)+
     $              dt*SUM(grdveq*pi_veq,1)
                ENDIF
              ENDDO
            ELSE
              IF (par_visc>0) THEN
                DO im=1,nmodes
                  vheat(1,:,:,im)=vheat(1,:,:,im)-
     $              dt*ndpt(1,:,:,im)*SUM(grdveq*(pi_veq+pareq),1)
                ENDDO
              ELSE
                DO im=1,nmodes
                  vheat(1,:,:,im)=vheat(1,:,:,im)-
     $              dt*ndpt(1,:,:,im)*SUM(grdveq*pi_veq,1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
c-----------------------------------------------------------------------
c       deallocate space for this block.
c-----------------------------------------------------------------------
        DEALLOCATE(grdv,piten,dtmp)
        IF (nonlinear) DEALLOCATE(real_pten,real_scal)
        IF (par_visc>0) DEALLOCATE(parten)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE viscous_heating
c-----------------------------------------------------------------------
c     subprogram 23. thermal_equil
c     Thermal equilibration between ions and electrons
c-----------------------------------------------------------------------
      SUBROUTINE thermal_equil
      USE local
      USE fields
      USE global
      USE input
      USE physdat
      USE closure_model_mod
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl
      REAL(r8) :: deni,dene,tfrac
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     adjust temperatures to model equilibration between ions and
c     electrons.  doing it here is simple, energy conserving, and
c     monotonic, even if it represents another form of time splitting.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      IF (tdep_tequil.AND.nonlinear) THEN
        CALL tdep_thermal_equil
      ELSE
        tfrac=dt*tequil_rate
        deni=1+tfrac*(zeff+1)
        dene=1+tfrac
        DO ibl=1,nbl
          CALL vector_add(tion(ibl),tele(ibl),
     $                    v1fac=dene/deni,v2fac=zeff*tfrac/deni)
          CALL vector_add(tele(ibl),tion(ibl),
     $                    v1fac=1/dene,v2fac=tfrac/dene)
        ENDDO

      ENDIF
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE thermal_equil
c-----------------------------------------------------------------------
c     subprogram 24. si_store
c     save fields and gradients that are used exclusively for the 3D
c     semi-implicit operator in the advance of V.
c-----------------------------------------------------------------------
      SUBROUTINE si_store
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE pardata
      USE fft_mod
      USE math_tran
      USE mpi_nim
      USE time
      IMPLICIT NONE

      INTEGER(i4) :: ibl,mxb,myb,ig,ng,nv,ierror,im
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     find expansion coefficients for pressure if this has not been
c     done.
c-----------------------------------------------------------------------
      IF (beta>0.AND.p_computation/='at nodes') CALL p_from_nt('all')
c-----------------------------------------------------------------------
c     store the gradient of magnetic field, the gradient of total
c     pressure, and total pressure as a function of toroidal-angle
c     index.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        ng=rb(ibl)%ng
        ALLOCATE(ctmp(9,ng,mxb*myb,nmodes))
        CALL math_grad(nmodes,keff,3_i4,geom,rb(ibl)%qbe%qpf,
     $                 rb(ibl)%qbe%qpfr,rb(ibl)%qbe%qpfz,
     $                 ctmp,bigr=rb(ibl)%bigr,vec0=rb(ibl)%qbe_eq%qpf,
     $                 dvec0r=rb(ibl)%qbe_eq%qpfr,
     $                 dvec0z=rb(ibl)%qbe_eq%qpfz)
        CALL qp_fft_noeq_save(ctmp,rb(ibl)%qgrdb%qpf,mxb,myb,
     $                        mpsq_block(ibl),9_i4,ng)
        DEALLOCATE(ctmp)
        IF (beta>0) THEN
          ALLOCATE(ctmp(3,ng,mxb*myb,nmodes))
          CALL math_grad(nmodes,keff,1_i4,geom,rb(ibl)%qpres%qpf,
     $                   rb(ibl)%qpres%qpfr,rb(ibl)%qpres%qpfz,
     $                   ctmp,bigr=rb(ibl)%bigr,
     $                   vec0  =rb(ibl)%qpres_eq%qpf,
     $                   dvec0r=rb(ibl)%qpres_eq%qpfr,
     $                   dvec0z=rb(ibl)%qpres_eq%qpfz)
          CALL qp_fft_noeq_save(ctmp,rb(ibl)%qgrdp%qpf,mxb,myb,
     $                          mpsq_block(ibl),3_i4,ng)
          DEALLOCATE(ctmp)
          CALL qp_fft_save(rb(ibl)%qpres%qpf,rb(ibl)%qpr_tot%qpf,
     $                     mxb,myb,mpsq_block(ibl),1_i4,
     $                     ng,rb(ibl)%qpres_eq%qpf)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        mxb=tb(ibl)%mcell
        myb=1_i4
        ng=tb(ibl)%ng
        ALLOCATE(ctmp(9,ng,mxb*myb,nmodes))
        CALL math_grad(nmodes,keff,3_i4,geom,tb(ibl)%qbe%qpf,
     $                 tb(ibl)%qbe%qpfr,tb(ibl)%qbe%qpfz,
     $                 ctmp,bigr=tb(ibl)%tgeom%bigr,
     $                 vec0  =tb(ibl)%qbe_eq%qpf,
     $                 dvec0r=tb(ibl)%qbe_eq%qpfr,
     $                 dvec0z=tb(ibl)%qbe_eq%qpfz)
        CALL qp_fft_noeq_save(ctmp,tb(ibl)%qgrdb%qpf,mxb,myb,
     $                        mpsq_block(ibl),9_i4,ng)
        DEALLOCATE(ctmp)
        IF (beta>0) THEN
          ALLOCATE(ctmp(3,ng,mxb*myb,nmodes))
          CALL math_grad(nmodes,keff,1_i4,geom,tb(ibl)%qpres%qpf,
     $                   tb(ibl)%qpres%qpfr,tb(ibl)%qpres%qpfz,
     $                   ctmp,bigr=tb(ibl)%tgeom%bigr,
     $                   vec0  =tb(ibl)%qpres_eq%qpf,
     $                   dvec0r=tb(ibl)%qpres_eq%qpfr,
     $                   dvec0z=tb(ibl)%qpres_eq%qpfz)
          CALL qp_fft_noeq_save(ctmp,tb(ibl)%qgrdp%qpf,mxb,myb,
     $                          mpsq_block(ibl),3_i4,ng)
          CALL qp_fft_save(tb(ibl)%qpres%qpf,tb(ibl)%qpr_tot%qpf,
     $                     mxb,myb,mpsq_block(ibl),1_i4,
     $                     ng,tb(ibl)%qpres_eq%qpf)
          DEALLOCATE(ctmp)
        ENDIF
      ENDDO
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE si_store
c-----------------------------------------------------------------------
c     subprogram 25. hyper_bheat
c     compute and save the energy loss from hyper-resistive and
c     hyper-divergence dissipation.  it will be added to the ohmic
c     heating.
c-----------------------------------------------------------------------
      SUBROUTINE hyper_bheat
      USE local
      USE physdat
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE math_tran
      USE pardata
      USE mpi_nim
      USE fft_mod
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ncx,mps
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: real_vec,real_scl
c-----------------------------------------------------------------------
c     the auxiliary vector from magnetic advances with hyper terms
c     is temporarily saved in work1.  interpolate it to the quad-point
c     storage and save its squared magnitude in qhyph.
c
c     complete rblocks first.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        ncx=rb(ibl)%mx*rb(ibl)%my*rb(ibl)%ng
        mps=mpsq_block(ibl)
        ALLOCATE(real_vec(3,mps,nphi))
        ALLOCATE(real_scl(1,mps,nphi))
        CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        CALL fft_nim('inverse',ncx,mps,lphi,3_i4,
     $               rb(ibl)%qwork1%qpf,real_vec,dealiase)
        real_scl(1,:,:)=SUM(real_vec(:,:,:)**2,1)*fhyp_eta/(mu0*dt)
        CALL fft_nim('forward',ncx,mps,lphi,1_i4,
     $               rb(ibl)%qhyph%qpf,real_scl,dealiase)
        DEALLOCATE(real_vec,real_scl)
      ENDDO
c-----------------------------------------------------------------------
c     tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        ncx=tb(ibl)%mcell*tb(ibl)%ng
        ALLOCATE(real_vec(3,mps,nphi))
        ALLOCATE(real_scl(1,mps,nphi))
        CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        CALL fft_nim('inverse',ncx,mps,lphi,3_i4,
     $               tb(ibl)%qwork1%qpf,real_vec,dealiase)
        real_scl(1,:,:)=SUM(real_vec(:,:,:)**2,1)*fhyp_eta/(mu0*dt)
        CALL fft_nim('forward',ncx,mps,lphi,1_i4,
     $               tb(ibl)%qhyph%qpf,real_scl,dealiase)
        DEALLOCATE(real_vec,real_scl)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE hyper_bheat
c-----------------------------------------------------------------------
c     subprogram 26. ndiff_store
c     when using the correction factor for artificial particle
c     diffusivity, quad-point storage over phi and the all-layer shared
c     symmetric (n0) part are saved twice during the implicit leapfrog.
c     those saves are collected here.
c-----------------------------------------------------------------------
      SUBROUTINE ndiff_store(ndchoice)
      USE local
      USE fields
      USE global
      USE input
      USE rblock
      USE tblock
      USE math_tran
      USE physdat
      USE fft_mod
      USE time
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: ndchoice

      INTEGER(i4) :: ibl
      REAL(r8) :: timest_qp,timend_qp
c-----------------------------------------------------------------------
c     interface block for qp0_bcast.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE qp0_bcast(threedata,twodata,nl)
        USE local
        COMPLEX(r8), DIMENSION(:,:,:,:) :: threedata
        REAL(r8), DIMENSION(:,:,:) :: twodata
        INTEGER(i4), INTENT(IN) :: nl
        END SUBROUTINE qp0_bcast
      END INTERFACE
c-----------------------------------------------------------------------
c     average works with qndiffa.
c-----------------------------------------------------------------------
      CALL timer(timest_qp)
      SELECT CASE(ndchoice)
      CASE('ave')
        DO ibl=1,nrbl
          CALL fft_nim('inverse',rb(ibl)%ng*rb(ibl)%mx*rb(ibl)%my,
     $                 mpsq_block(ibl),lphi,1_i4,
     $                 rb(ibl)%qndiffa%qpf,
     $                 rb(ibl)%qndiff_phi%qpf,dealiase)
          CALL qp0_bcast(rb(ibl)%qndiffa%qpf,rb(ibl)%qndiff_n0%qpf,
     $                   nlayers)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL fft_nim('inverse',tb(ibl)%ng*tb(ibl)%mcell,
     $                 mpsq_block(ibl),lphi,1_i4,
     $                 tb(ibl)%qndiffa%qpf,
     $                 tb(ibl)%qndiff_phi%qpf,dealiase)
          CALL qp0_bcast(tb(ibl)%qndiffa%qpf,tb(ibl)%qndiff_n0%qpf,
     $                   nlayers)
        ENDDO
c-----------------------------------------------------------------------
c     end works with qndiff.
c-----------------------------------------------------------------------
      CASE('end')
        DO ibl=1,nrbl
          CALL fft_nim('inverse',rb(ibl)%ng*rb(ibl)%mx*rb(ibl)%my,
     $                 mpsq_block(ibl),lphi,1_i4,
     $                 rb(ibl)%qndiff%qpf,
     $                 rb(ibl)%qndiff_phi%qpf,dealiase)
          CALL qp0_bcast(rb(ibl)%qndiff%qpf,rb(ibl)%qndiff_n0%qpf,
     $                   nlayers)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL fft_nim('inverse',tb(ibl)%ng*tb(ibl)%mcell,
     $                 mpsq_block(ibl),lphi,1_i4,
     $                 tb(ibl)%qndiff%qpf,
     $                 tb(ibl)%qndiff_phi%qpf,dealiase)
          CALL qp0_bcast(tb(ibl)%qndiff%qpf,tb(ibl)%qndiff_n0%qpf,
     $                   nlayers)
        ENDDO
      END SELECT
      CALL timer(timend_qp)
      time_qpfld=time_qpfld+timend_qp-timest_qp
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ndiff_store
