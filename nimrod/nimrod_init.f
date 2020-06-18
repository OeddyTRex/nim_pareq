c-----------------------------------------------------------------------
c     file nimrod_init.f:  contains the initialization routines for
c     nimrod (routines needed after reading a restart dump).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  cell_init.
c     2.  boundary_init.
c     3.  variable_alloc.
c     4.  mass_mat_init.
c     5.  quadrature_save.
c     6.  e_applied_init.
c     7.  boundary_vals_init.
c     8.  matrix_init.
c     9.  pointer_init.
c     10. q_applied_init.
c     11. fourier_init.
c     12. proj_mat_init.
c     13. hv_mat_init.
c-----------------------------------------------------------------------
c     subprogram 1. cell_init.
c     computes cell-related geometric quantities.
c-----------------------------------------------------------------------
      SUBROUTINE cell_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE input
      USE global
      USE integrands
      USE computation_pointers
      USE mpi_nim
      USE pardata
      IMPLICIT NONE
      
      INTEGER(i4) :: ibl,mxb,myb,ix,iy,mc,ic,iv0,iv1,ierror,rem,ng,ig,
     $               ipol
      REAL(r8) :: kmax,tmp,dnode,minr
      REAL(r8), DIMENSION(2) :: drzdx,drzdy
c-----------------------------------------------------------------------
c     find the centroid radius and cell volume.  also save the total
c     volume and cross section area.  the integral of dA/r is used
c     for the computation of the reversal parameter.
c-----------------------------------------------------------------------
      total_volume=0._r8
      cross_section=0._r8
      cross_s_overr=0._r8
      DO ibl=1,nrbl
        CALL rblock_get_rhs(rb(ibl),cell_rhs(ibl),get_vol,3_i4)
        rb(ibl)%cell_vol(:,:)=cell_rhs(ibl)%arri(1,1,:,:)
        rb(ibl)%r_cent(:,:)=cell_rhs(ibl)%arri(1,1,:,:)
     $                     /cell_rhs(ibl)%arri(2,1,:,:)
        total_volume=total_volume+SUM(cell_rhs(ibl)%arri(1,1,:,:))
        cross_section=cross_section+SUM(cell_rhs(ibl)%arri(2,1,:,:))
        cross_s_overr=cross_s_overr+SUM(cell_rhs(ibl)%arri(3,1,:,:))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_get_rhs(tb(ibl),cell_rhs(ibl),get_vol,3_i4)
        tb(ibl)%cell_vol(:)=cell_rhs(ibl)%arri(1,1,:,1)
        tb(ibl)%r_cent(:)=cell_rhs(ibl)%arri(1,1,:,1)
     $                   /cell_rhs(ibl)%arri(2,1,:,1)
        total_volume=total_volume+SUM(cell_rhs(ibl)%arri(1,1,:,1))
        cross_section=cross_section+SUM(cell_rhs(ibl)%arri(2,1,:,1))
        cross_s_overr=cross_s_overr+SUM(cell_rhs(ibl)%arri(3,1,:,1))
      ENDDO
c-----------------------------------------------------------------------
c     sum over processors and multiply by a geometric factor.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(total_volume,tmp,1,mpi_nim_real,
     $       mpi_sum,comm_layer,ierror)
        total_volume=tmp
        CALL mpi_allreduce(cross_section,tmp,1,mpi_nim_real,
     $       mpi_sum,comm_layer,ierror)
        cross_section=tmp
        CALL mpi_allreduce(cross_s_overr,tmp,1,mpi_nim_real,
     $       mpi_sum,comm_layer,ierror)
        cross_s_overr=tmp
      ENDIF
      IF (geom=='lin') THEN
        total_volume=total_volume*per_length
      ELSE
        total_volume=total_volume*twopi
      ENDIF
c-----------------------------------------------------------------------
c     find the equivalent square of the minimum grid-dimension for
c     CFL calculations.  the computations include a factor of pi**2, so
c     delta2 is now more like max(k**2).
c
c     this is used for the explicit hyper-viscosity.  do not include
c     elements that touch a degenerate point, where the hyper-viscosity
c     coefficient is set to 0.
c-----------------------------------------------------------------------
      delta2=HUGE(delta2)
      DO ibl=1,nrbl
        IF (poly_degree==1) THEN
          dnode=1._r8
        ELSE
          dnode=MINVAL(rb(ibl)%rz%dx(2:poly_degree)-
     $                 rb(ibl)%rz%dx(1:poly_degree-1))
        ENDIF
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        IF (rb(ibl)%degenerate) THEN
          iv0=2
        ELSE
          iv0=1
        ENDIF
        ipol=1
        DO iy=1,myb
          DO ix=iv0,mxb
            drzdx=0.5_r8*dnode*
     $        (rb(ibl)%rz%fs(:,ix  ,iy-1)+rb(ibl)%rz%fs(:,ix  ,iy)
     $        -rb(ibl)%rz%fs(:,ix-1,iy-1)-rb(ibl)%rz%fs(:,ix-1,iy))
            drzdy=0.5_r8*dnode*
     $        (rb(ibl)%rz%fs(:,ix-1,iy  )+rb(ibl)%rz%fs(:,ix,iy  )
     $        -rb(ibl)%rz%fs(:,ix-1,iy-1)-rb(ibl)%rz%fs(:,ix,iy-1))
            minr=HUGE(minr)
            DO ig=1,SIZE(rb(ibl)%bigr,1)
              IF (rb(ibl)%bigr(ig,ipol)>0._r8)
     $          minr=MIN(minr,rb(ibl)%bigr(ig,ipol))
            ENDDO
            kmax=MAXVAL(ABS(keff))/minr
            delta2=MIN(delta2,1._r8/
     $             (pi**2/SUM(drzdx**2)+pi**2/SUM(drzdy**2)+kmax**2))
            ipol=ipol+1
          ENDDO
        ENDDO
      ENDDO
      DO ibl=nrbl+1,nbl
        mc=tb(ibl)%mcell
        DO ic=1,mc
          minr=HUGE(minr)
          DO ig=1,SIZE(tb(ibl)%tgeom%bigr,1)
            IF (tb(ibl)%tgeom%bigr(ig,ic)>0._r8)
     $        minr=MIN(minr,tb(ibl)%tgeom%bigr(ig,ic))
          ENDDO
          kmax=MAXVAL(ABS(keff))/minr
          DO iv1=1,3
            iv0=iv1-1
            IF (iv0==0) iv0=3
            drzdx=(/tb(ibl)%tgeom%xs(tb(ibl)%tgeom%vertex(ic,iv1))
     $             -tb(ibl)%tgeom%xs(tb(ibl)%tgeom%vertex(ic,iv0)),
     $              tb(ibl)%tgeom%ys(tb(ibl)%tgeom%vertex(ic,iv1))
     $             -tb(ibl)%tgeom%ys(tb(ibl)%tgeom%vertex(ic,iv0))/)
            delta2=MIN(delta2,
     $             1._r8/(4._r8/SUM(drzdx**2)+kmax**2))
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     minimize over processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(delta2,tmp,1,mpi_nim_real,
     $       mpi_min,comm_layer,ierror)
        delta2=tmp
      ENDIF
c-----------------------------------------------------------------------
c     finally, determine the poloidal dimension size for pseudo-spectral
c     computations.  the ipolst_block and ipolen_block are the start and
c     end index for each layer out of the 1D list of cells (1:mc).
c
c     the second set of arrays (with 'q' in the names) are similar but
c     fold the quadrature-point indices into the poloidal index.  they
c     are now used in integrands and field_comps to reduce the number
c     of communications for the ffts.
c-----------------------------------------------------------------------
      ALLOCATE(mps_block(nbl))
      ALLOCATE(ipolst_block(nbl))
      ALLOCATE(ipolen_block(nbl))
      ALLOCATE(mpsq_block(nbl))
      ALLOCATE(ipqst_block(nbl))
      ALLOCATE(ipqen_block(nbl))
      DO ibl=1,nbl
        IF (ibl<=nrbl) THEN
          mc=rb(ibl)%mx*rb(ibl)%my
          ng=rb(ibl)%ng
        ELSE
          mc=tb(ibl)%mcell
          ng=tb(ibl)%ng
        ENDIF
        rem=MODULO(mc,nlayers)
        mps_block(ibl)=mc/nlayers
        ipolst_block(ibl)=ilayer*mps_block(ibl)+1+MIN(rem,ilayer)
        ipolen_block(ibl)=(ilayer+1)*mps_block(ibl)+MIN(rem,ilayer+1)
        IF (ilayer<rem) mps_block(ibl)=mps_block(ibl)+1
        rem=MODULO(mc*ng,nlayers)
        mpsq_block(ibl)=(mc*ng)/nlayers
        ipqst_block(ibl)=ilayer*mpsq_block(ibl)+1+MIN(rem,ilayer)
        ipqen_block(ibl)=(ilayer+1)*mpsq_block(ibl)+MIN(rem,ilayer+1)
        IF (ilayer<rem) mpsq_block(ibl)=mpsq_block(ibl)+1
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cell_init
c-----------------------------------------------------------------------
c     subprogram 2. boundary_init.
c     initializes arrays used for setting boundary and regularity
c     conditions.
c-----------------------------------------------------------------------
      SUBROUTINE boundary_init(geom)
      USE local
      USE fields
      USE seam_storage_mod
      USE pardata
      USE mpi_nim
      USE computation_pointers
      USE regularity
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: iv,ip,np,ibl,ibv,mx,my,ix,iy,nexbl,nv,ibe,
     $               iblocal,nr0bl,in,ierror
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ivp
      LOGICAL, DIMENSION(nbl) :: blmask,blr0mask
      REAL(r8), DIMENSION(2) :: rpandn
c-----------------------------------------------------------------------
c     create an external block mask and allocate communication arrays.  
c     set blmask to .true. ONLY if I own the block on the exterior 
c     boundary
c-----------------------------------------------------------------------
      blmask=.false.
      DO iv=1,seam0%nvert
        np=SIZE(seam0%vertex(iv)%ptr,2)
        DO ip=1,np
          ibl=seam0%vertex(iv)%ptr(1,ip)
          IF (block2proc(ibl) == node) blmask(global2local(ibl))=.true.
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     save the major radius for toroidal geometry--it is part of the
c     basis vector in certain equations and it gets used for setting up
c     regularity conditions.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          IF (geom=='tor') THEN
            IF(ibl>nrbl)THEN
              seam(ibl)%vertex(iv)%rgeom=tb(ibl)%tgeom%xs(ix)
            ELSE
              seam(ibl)%vertex(iv)%rgeom=rb(ibl)%rz%fs(1,ix,iy)
            ENDIF
          ELSE
            seam(ibl)%vertex(iv)%rgeom=1
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     create an array of block numbers which touch the exterior
c     boundary. 
c     NOTE - my local block numbers are stored in exblock_list
c            NOT global block numbers
c-----------------------------------------------------------------------
      nexbl=0
      DO ibl=1,nbl
        IF (blmask(ibl)) nexbl=nexbl+1
      ENDDO
      ALLOCATE(exblock_list(nexbl))
      nexbl=0
      DO ibl=1,nbl
        IF (blmask(ibl)) THEN
          nexbl=nexbl+1
          exblock_list(nexbl)=ibl
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     create external point and corner masks in seams that touch the
c     external seam.  work on a seam's vertex ONLY if I own the block 
c     to which seam0 is connected.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        nv=seam(ibl)%nvert
        ALLOCATE(seam(ibl)%expoint(nv))
        seam(ibl)%expoint=.false.
        ALLOCATE(seam(ibl)%excorner(nv))
        seam(ibl)%excorner=.false.
      ENDDO
      DO iv=1,seam0%nvert
        np=SIZE(seam0%vertex(iv)%ptr,2)
        DO ip=1,np
          ibl=seam0%vertex(iv)%ptr(1,ip)
          IF (block2proc(ibl) == node) THEN
            iblocal = global2local(ibl)
            ibv=seam0%vertex(iv)%ptr(2,ip)
            seam(iblocal)%expoint(ibv)=.true.
            seam(iblocal)%excorner(ibv)=seam0%excorner(iv)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find r=0 points, set the r0point flags, and reset expoint.
c     the r=0 points are assumed to lie on the edge of the domain;
c     hence the blocks are in exblock_list at this point.
c     this needs to be before block_create_tang.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        ALLOCATE(seam(ibl)%r0point(seam(ibl)%nvert))
        seam(ibl)%r0point=.false.
        IF (ibl<=nrbl) rb(ibl)%r0block=.false.
      ENDDO
      ALLOCATE(r0block_list(nbl))
      r0block_list=0
      blr0mask=.false.
      IF (geom=='tor') THEN
        DO ibe=1,SIZE(exblock_list)
          ibl=exblock_list(ibe)
          nv=seam(ibl)%nvert
          DO iv=1,nv
            ip=iv-1
            in=iv+1
            IF (in>nv) in=1
            IF (ip<1 ) ip=nv
            IF (seam(ibl)%vertex(iv)%rgeom<=1.e-12) THEN
              seam(ibl)%r0point(iv)=.true.
              seam(ibl)%excorner(iv)=.false.
              r0block_list(ibl)=1
              IF (ibl<=nrbl) rb(ibl)%r0block=.true.
              blr0mask(ibl)=.true.
              rpandn=(/seam(ibl)%vertex(ip)%rgeom,
     $                 seam(ibl)%vertex(in)%rgeom/)
              IF (SUM(rpandn)<=1.e-12.OR.
     $          (rpandn(1)<=1.e-12.AND..NOT.seam(ibl)%expoint(in)).OR.
     $          (rpandn(2)<=1.e-12.AND..NOT.seam(ibl)%expoint(ip))) THEN
                seam(ibl)%expoint(iv)=.false.
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     set the r0block list, and reset the external block list.
c-----------------------------------------------------------------------
      nr0bl=SUM(r0block_list)
      DEALLOCATE(r0block_list)
      ALLOCATE(r0block_list(nr0bl))
      nexbl=SIZE(exblock_list)
      IF (geom=='tor') THEN
        ibv=1
        r0bl_loop: DO ibl=1,nbl
          IF (blr0mask(ibl)) THEN
            r0block_list(ibv)=ibl
            ibv=ibv+1
            DO iv=1,seam(ibl)%nvert
              IF (seam(ibl)%expoint(iv)) CYCLE r0bl_loop
            ENDDO
            nexbl=nexbl-1
            blmask(ibl)=.false.
          ENDIF
        ENDDO r0bl_loop
        IF (nexbl<SIZE(exblock_list)) THEN
          DEALLOCATE(exblock_list)
          ALLOCATE(exblock_list(nexbl))
          nexbl=0
          DO ibl=1,nbl
            IF (blmask(ibl)) THEN
              nexbl=nexbl+1
              exblock_list(nexbl)=ibl
            ENDIF
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     determine if any processors have r0 blocks.
c-----------------------------------------------------------------------
      CALL mpi_allreduce(nr0bl,ibl,1,mpi_nim_int,
     $     mpi_sum,mpi_comm_world,ierror)
      IF (ibl>0) THEN
        any_r0blocks=.true.
      ELSE
        any_r0blocks=.false.
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE boundary_init
c-----------------------------------------------------------------------
c     subprogram 3. variable_alloc.
c     allocates space for the dependent variables.
c-----------------------------------------------------------------------
      SUBROUTINE variable_alloc
      USE local
      USE fields
      USE global
      USE input
      USE computation_pointers
      USE rblock
      USE tblock
      USE physdat
      IMPLICIT NONE

      INTEGER(i4) :: ibl,mxb,myb,im,mv,mc,iv,iq,jq,nn,ix,iy,ibasis,
     $               ix0,iy0
      REAL(r8), DIMENSION(1) :: te,ti
      LOGICAL :: q3alloc=.false.
c-----------------------------------------------------------------------
c     allocate generic blocks for finite element computations.
c-----------------------------------------------------------------------
      ALLOCATE(rhs(nbl))
      ALLOCATE(crhs(nbl))
      ALLOCATE(cvecn(nbl))
      ALLOCATE(cell_rhs(nbl))
      ALLOCATE(cell_crhs(nbl))
      ALLOCATE(sln(nbl))
      ALLOCATE(csln(nbl))
      ALLOCATE(vectr(nbl))
      ALLOCATE(cvectr(nbl))
      ALLOCATE(si_nl_pres(nbl))
c-----------------------------------------------------------------------
c     nullify structures used for algebraic operations during the
c     finite element computations.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        NULLIFY(crhs(ibl)%arr ,crhs(ibl)%arrh,
     $          crhs(ibl)%arrv,crhs(ibl)%arri,crhs(ibl)%arrtmp)
        NULLIFY(cvecn(ibl)%arr ,cvecn(ibl)%arrh,
     $          cvecn(ibl)%arrv,cvecn(ibl)%arri,cvecn(ibl)%arrtmp)
        NULLIFY(sln(ibl)%arr ,sln(ibl)%arrh,
     $          sln(ibl)%arrv,sln(ibl)%arri,sln(ibl)%arrtmp)
        NULLIFY(csln(ibl)%arr ,csln(ibl)%arrh,
     $          csln(ibl)%arrv,csln(ibl)%arri,csln(ibl)%arrtmp)
        NULLIFY(vectr(ibl)%arr ,vectr(ibl)%arrh,
     $          vectr(ibl)%arrv,vectr(ibl)%arri,vectr(ibl)%arrtmp)
        NULLIFY(cvectr(ibl)%arr ,cvectr(ibl)%arrh,
     $          cvectr(ibl)%arrv,cvectr(ibl)%arri,cvectr(ibl)%arrtmp)
      ENDDO
c-----------------------------------------------------------------------
c     allocate storage for the real and imaginary parts of non-
c     fundamental dependent variables.  start with rectangular blocks.
c-----------------------------------------------------------------------
      rblocks: DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL lagr_quad_alloc(rb(ibl)%ja,mxb,myb,3_i4,nmodes,
     $                       poly_degree,'ja',(/'  ja  '/))
        rb(ibl)%ja=0
c-----------------------------------------------------------------------
c       allocate storage for dependent variables at the gaussian
c       quadrature points.
c-----------------------------------------------------------------------
        CALL rblock_qp_alloc(rb(ibl)%qbe,rb(ibl),3_i4,nmodes)
        IF (ohms=='2fl'.AND.advect=='all'.OR.
     $      separate_pe.AND.nonlinear) THEN
          CALL rblock_qp_alloc(rb(ibl)%qja,rb(ibl),3_i4,nmodes)
        ENDIF
        CALL rblock_qp_alloc(rb(ibl)%qve,rb(ibl),3_i4,nmodes)
        IF (beta>0) THEN
          CALL rblock_qp_alloc(rb(ibl)%qpres,rb(ibl),1_i4,nmodes)
          CALL rblock_qp_alloc(rb(ibl)%qprese,rb(ibl),1_i4,nmodes)
          CALL rblock_qp_alloc(rb(ibl)%qtion,rb(ibl),1_i4,nmodes)
          CALL rblock_qp_alloc(rb(ibl)%qtele,rb(ibl),1_i4,nmodes)
          CALL rblock_qp_alloc(rb(ibl)%qwork3,rb(ibl),1_i4,nmodes)
          q3alloc=.true.
        ENDIF
        IF (beta>0.AND.p_model(1:5)=='aniso'.OR.
     $      par_visc>0.OR.gyr_visc>0)
     $    CALL rblock_qp_alloc(rb(ibl)%qbb,rb(ibl),6_i4)
        IF (nonlinear.OR.eq_flow/='none'.OR.ohms/='mhd')
     $    CALL rblock_qp_alloc(rb(ibl)%qwork1,rb(ibl),3_i4,nmodes)
        CALL rblock_qp_alloc(rb(ibl)%qnd,rb(ibl),1_i4,nmodes)
        IF (continuity/='none') THEN
          IF (.NOT.q3alloc)
     $      CALL rblock_qp_alloc(rb(ibl)%qwork3,rb(ibl),1_i4,nmodes)
        ENDIF
c-----------------------------------------------------------------------
c       allocate space for the fields used for creating semi-implicit
c       operators.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          CALL lagr_quad_alloc(rb(ibl)%be_n0,mxb,myb,3_i4,
     $                         1_i4,'b0',(/' be_n0'/))
          CALL rblock_qp_alloc(rb(ibl)%qbe_n0,rb(ibl),3_i4)
          IF (beta>0) THEN
            CALL lagr_quad_alloc(rb(ibl)%pres_n0,mxb,myb,1_i4,
     $                           1_i4,'p0',(/'prs_n0'/))
            CALL rblock_qp_alloc(rb(ibl)%qpres_n0,rb(ibl),1_i4)
            CALL rblock_qp_alloc(rb(ibl)%qti_n0,rb(ibl),1_i4)
            CALL rblock_qp_alloc(rb(ibl)%qte_n0,rb(ibl),1_i4)
          ENDIF
          CALL vector_type_alloc(si_nl_pres(ibl),0_i4,mxb,myb,2_i4)
          CALL rblock_qp_alloc(rb(ibl)%qsi_nl_pres,rb(ibl),2_i4)
          si_nl_pres(ibl)%arri=0._r8
          rb(ibl)%qsi_nl_pres%qpf=0._r8
          IF (continuity/='none') THEN
            CALL lagr_quad_alloc(rb(ibl)%nd_n0,mxb,myb,1_i4,
     $                           1_i4,'n0',(/' nd_n0'/))
            CALL rblock_qp_alloc(rb(ibl)%qnd_n0,rb(ibl),1_i4)
          ENDIF
          IF (impladv) THEN
            CALL lagr_quad_alloc(rb(ibl)%ve_n0,mxb,myb,3_i4,
     $                           1_i4,'v0',(/' ve_n0'/))
            CALL rblock_qp_alloc(rb(ibl)%qve_n0,rb(ibl),3_i4)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       allocate storage for temporary computations.
c-----------------------------------------------------------------------
        CALL lagr_quad_alloc(rb(ibl)%work1,mxb,myb,3_i4,nmodes,
     $                       poly_degree,'w1',(/'work1 '/))
        CALL lagr_quad_alloc(rb(ibl)%work2,mxb,myb,1_i4,nmodes,
     $                       poly_degree,'w2',(/'work2 '/))
        CALL lagr_quad_alloc(rb(ibl)%work3,mxb,myb,1_i4,nmodes,
     $                       poly_degree,'w3',(/'work3 '/))
        CALL lagr_quad_alloc(rb(ibl)%work4,mxb,myb,3_i4,nmodes,
     $                       poly_degree,'w4',(/'work4 '/))
        CALL lagr_quad_alloc(rb(ibl)%rwork1,mxb,myb,1_i4,
     $                       poly_degree,'r1',(/'rwork1'/))
        rb(ibl)%work1=0
        rb(ibl)%work2=0
        rb(ibl)%work3=0
        rb(ibl)%work4=0
        rb(ibl)%rwork1=0
        IF (nonlinear.AND.impladv.AND.nd_hypd>0) THEN
          CALL lagr_quad_alloc(rb(ibl)%work5,mxb,myb,2_i4,nmodes,
     $                         poly_degree,'w5',(/'work5 '/))
          CALL lagr_quad_alloc(rb(ibl)%work6,mxb,myb,2_i4,nmodes,
     $                         poly_degree,'w6',(/'work6 '/))
          rb(ibl)%work5=0
          rb(ibl)%work6=0
        ENDIF
        IF (nonlinear.AND.(hyp_eta>0.OR.hyp_dbd>0).AND.
     $    .NOT.split_hypeta) THEN
          CALL lagr_quad_alloc(rb(ibl)%w6v1,mxb,myb,6_i4,nmodes,
     $                         poly_degree,'v1',(/' w6v1 '/))
          CALL lagr_quad_alloc(rb(ibl)%w6v2,mxb,myb,6_i4,nmodes,
     $                         poly_degree,'v2',(/' w6v2 '/))
          rb(ibl)%w6v1=0
          rb(ibl)%w6v2=0
        ENDIF
c-----------------------------------------------------------------------
c       discontinuous modal structures and work structures for 3d solves
c       are needed for auxiliary fields.
c-----------------------------------------------------------------------
        IF (poly_divb>=0) THEN
          CALL modal_disc_alloc(rb(ibl)%auxb,mxb,myb,1_i4,nmodes,
     $                          poly_divb,poly_divb_min,
     $                          poly_divb_max,name='ab',
     $                          title=(/' auxb '/))
          CALL modal_disc_assignment(rb(ibl)%auxb,0._r8)
        ELSE
          rb(ibl)%auxb%nqty=0
          rb(ibl)%auxb%name='ab'
        ENDIF
        IF (poly_divv>=0) THEN
          CALL modal_disc_alloc(rb(ibl)%auxv,mxb,myb,2_i4,nmodes,
     $                          poly_divv,poly_divv_min,
     $                          poly_divv_max,name='av',
     $                          title=(/' auxv '/))
          CALL modal_disc_assignment(rb(ibl)%auxv,0._r8)
        ELSE
          rb(ibl)%auxv%nqty=0
          rb(ibl)%auxv%name='av'
        ENDIF

        IF (poly_divb>0.AND.nonlinear) THEN
          CALL modal_disc_alloc(rb(ibl)%mwork1,rb(ibl)%mx,rb(ibl)%my,
     $                          1_i4,nmodes,poly_divb,poly_divb_min,
     $                          poly_divb_max)
          CALL modal_disc_alloc(rb(ibl)%mwork2,rb(ibl)%mx,rb(ibl)%my,
     $                          1_i4,nmodes,poly_divb,poly_divb_min,
     $                          poly_divb_max)
        ENDIF
        IF (poly_divv>0.AND.nonlinear) THEN
          CALL modal_disc_alloc(rb(ibl)%mwork3,rb(ibl)%mx,rb(ibl)%my,
     $                          2_i4,nmodes,poly_divv,poly_divv_min,
     $                          poly_divv_max)
          CALL modal_disc_alloc(rb(ibl)%mwork4,rb(ibl)%mx,rb(ibl)%my,
     $                          2_i4,nmodes,poly_divv,poly_divv_min,
     $                          poly_divv_max)
        ENDIF
c-----------------------------------------------------------------------
c       create lagrange_quad structures for equilibrium
c       electron and ion temperatures.
c-----------------------------------------------------------------------
        CALL lagr_quad_alloc(rb(ibl)%tele_eq,mxb,myb,1_i4,
     $                       poly_degree,'teleeq',(/'teleeq'/))
        CALL lagr_quad_alloc(rb(ibl)%tion_eq,mxb,myb,1_i4,
     $                       poly_degree,'tioneq',(/'tioneq'/))
        DO ibasis=1,SIZE(rb(ibl)%tele_eq%dx)
          ix0=rb(ibl)%tele_eq%ix0(ibasis)
          iy0=rb(ibl)%tele_eq%iy0(ibasis)
          DO iy=iy0,myb
            DO ix=ix0,mxb
              CALL lagr_quad_eval(rb(ibl)%pres_eq,
     $                         ix-ix0+rb(ibl)%tele_eq%dx(ibasis),
     $                         iy-iy0+rb(ibl)%tele_eq%dy(ibasis),0_i4)
              CALL lagr_quad_eval(rb(ibl)%prese_eq,
     $                         ix-ix0+rb(ibl)%tele_eq%dx(ibasis),
     $                         iy-iy0+rb(ibl)%tele_eq%dy(ibasis),0_i4)
              CALL lagr_quad_eval(rb(ibl)%nd_eq,
     $                         ix-ix0+rb(ibl)%tele_eq%dx(ibasis),
     $                         iy-iy0+rb(ibl)%tele_eq%dy(ibasis),0_i4)
              te=MAX(smallnum,
     $               rb(ibl)%prese_eq%f(1)/(kboltz*rb(ibl)%nd_eq%f(1)))
              ti=MAX(smallnum,
     $               (rb(ibl)%pres_eq%f(1)-rb(ibl)%prese_eq%f(1))*
     $               zeff/(kboltz*rb(ibl)%nd_eq%f(1)))
              CALL lagr_quad_basis_assign_loc
     $          (rb(ibl)%tele_eq,te,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc
     $          (rb(ibl)%tion_eq,ti,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       allocate storage for rhs and solution of equations.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(rhs(ibl),poly_degree,mxb,myb,3_i4)
        CALL vector_type_alloc(cell_rhs(ibl),0_i4,mxb,myb,3_i4)
        CALL vector_type_alloc(cell_crhs(ibl),0_i4,mxb,myb,4_i4,nmodes)
c-----------------------------------------------------------------------
c       allocate storage for cell volume and centroid.
c-----------------------------------------------------------------------
        ALLOCATE(rb(ibl)%cell_vol(mxb,myb))
        ALLOCATE(rb(ibl)%r_cent(mxb,myb))
c-----------------------------------------------------------------------
c       storage for viscous heating density and hyper-resistive heating.
c-----------------------------------------------------------------------
        IF (visc_heat)
     $    CALL rblock_qp_alloc(rb(ibl)%qvisc,rb(ibl),1_i4,nmodes)
        IF (ohm_heat.AND.(hyp_eta>0._r8.OR.hyp_dbd>0._r8))
     $    CALL rblock_qp_alloc(rb(ibl)%qhyph,rb(ibl),1_i4,nmodes)
      ENDDO rblocks
c-----------------------------------------------------------------------
c     blocks of unstructured triangles, non-fundamentals first:
c-----------------------------------------------------------------------
      tblocks: DO ibl=nrbl+1,nbl
        mv=tb(ibl)%mvert
        mc=tb(ibl)%mcell
        CALL tri_linear_alloc(tb(ibl)%ja,mv,3_i4,nmodes,'ja',
     $                        (/'  ja  '/))
        tb(ibl)%ja=0
c-----------------------------------------------------------------------
c       allocate storage for dependent variables at the gaussian
c       quadrature points.
c-----------------------------------------------------------------------
        CALL tblock_qp_alloc(tb(ibl)%qbe,tb(ibl),3_i4,nmodes)
        IF (ohms=='2fl'.AND.advect=='all'.OR.
     $      separate_pe.AND.nonlinear) THEN
          CALL tblock_qp_alloc(tb(ibl)%qja,tb(ibl),3_i4,nmodes)
        ENDIF
        CALL tblock_qp_alloc(tb(ibl)%qve,tb(ibl),3_i4,nmodes)
        IF (beta>0) THEN
          CALL tblock_qp_alloc(tb(ibl)%qpres,tb(ibl),1_i4,nmodes)
          CALL tblock_qp_alloc(tb(ibl)%qprese,tb(ibl),1_i4,nmodes)
          CALL tblock_qp_alloc(tb(ibl)%qtion,tb(ibl),1_i4,nmodes)
          CALL tblock_qp_alloc(tb(ibl)%qtele,tb(ibl),1_i4,nmodes)
          CALL tblock_qp_alloc(tb(ibl)%qwork3,tb(ibl),1_i4,nmodes)
          q3alloc=.true.
        ENDIF
        IF (beta>0.AND.p_model(1:5)=='aniso'.OR.
     $      par_visc>0.OR.gyr_visc>0)
     $    CALL tblock_qp_alloc(tb(ibl)%qbb,tb(ibl),6_i4)
        IF (nonlinear.OR.eq_flow/='none'.OR.ohms/='mhd')
     $    CALL tblock_qp_alloc(tb(ibl)%qwork1,tb(ibl),3_i4,nmodes)
        CALL tblock_qp_alloc(tb(ibl)%qnd,tb(ibl),1_i4,nmodes)
        IF (continuity/='none') THEN
          IF (.NOT.q3alloc)
     $      CALL tblock_qp_alloc(tb(ibl)%qwork3,tb(ibl),1_i4,nmodes)
        ENDIF
c-----------------------------------------------------------------------
c       allocate space for the fields used for creating semi-implicit
c       operators.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          CALL tri_linear_alloc(tb(ibl)%be_n0,mv,3_i4,'b0',
     $                          (/' be_n0'/))
          CALL tblock_qp_alloc(tb(ibl)%qbe_n0,tb(ibl),3_i4)
          IF (beta>0) THEN
            CALL tri_linear_alloc(tb(ibl)%pres_n0,mv,1_i4,'p0',
     $                            (/'prs_n0'/))
            CALL tblock_qp_alloc(tb(ibl)%qpres_n0,tb(ibl),1_i4)
            CALL tblock_qp_alloc(tb(ibl)%qti_n0,tb(ibl),1_i4)
            CALL tblock_qp_alloc(tb(ibl)%qte_n0,tb(ibl),1_i4)
          ENDIF
          CALL vector_type_alloc(si_nl_pres(ibl),0_i4,mc,1_i4,2_i4)
          CALL tblock_qp_alloc(tb(ibl)%qsi_nl_pres,tb(ibl),2_i4)
          si_nl_pres(ibl)%arri=0._r8
          tb(ibl)%qsi_nl_pres%qpf=0._r8
          IF (continuity/='none') THEN
            CALL tri_linear_alloc(tb(ibl)%nd_n0,mv,1_i4,'n0',
     $                            (/' nd_n0'/))
            CALL tblock_qp_alloc(tb(ibl)%qnd_n0,tb(ibl),1_i4)
          ENDIF
          IF (impladv) THEN
            CALL tri_linear_alloc(tb(ibl)%ve_n0,mv,3_i4,'v0',
     $                            (/' ve_n0'/))
            CALL tblock_qp_alloc(tb(ibl)%qve_n0,tb(ibl),3_i4)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       allocate storage for temporary computations.
c-----------------------------------------------------------------------
        CALL tri_linear_alloc(tb(ibl)%work1,mv,3_i4,nmodes,'w1',
     $                        (/'work1 '/))
        CALL tri_linear_alloc(tb(ibl)%work2,mv,1_i4,nmodes,'w2',
     $                        (/'work2 '/))
        CALL tri_linear_alloc(tb(ibl)%work3,mv,1_i4,nmodes,'w3',
     $                        (/'work3 '/))
        CALL tri_linear_alloc(tb(ibl)%work4,mv,3_i4,nmodes,'w4',
     $                        (/'work4 '/))
        CALL tri_linear_alloc(tb(ibl)%rwork1,mv,1_i4,'r1',(/'rwork1'/))
        tb(ibl)%work1=0
        tb(ibl)%work2=0
        tb(ibl)%work3=0
        tb(ibl)%work4=0
        tb(ibl)%rwork1=0
        IF (nonlinear.AND.impladv.AND.nd_hypd>0) THEN
          CALL tri_linear_alloc(tb(ibl)%work5,mv,2_i4,nmodes,'w5',
     $                          (/'work5 '/))
          CALL tri_linear_alloc(tb(ibl)%work6,mv,2_i4,nmodes,'w6',
     $                          (/'work6 '/))
          tb(ibl)%work5=0
          tb(ibl)%work6=0
        ENDIF
        IF (nonlinear.AND.(hyp_eta>0.OR.hyp_dbd>0).AND.
     $    .NOT.split_hypeta) THEN
          CALL tri_linear_alloc(tb(ibl)%w6v1,mv,6_i4,nmodes,
     $                         'v1',(/' w6v1 '/))
          CALL tri_linear_alloc(tb(ibl)%w6v2,mv,6_i4,nmodes,
     $                         'v2',(/' w6v2 '/))
          tb(ibl)%w6v1=0
          tb(ibl)%w6v2=0
        ENDIF
c-----------------------------------------------------------------------
c       create tri_linear structures for equilibrium
c       electron and ion temperatures.
c-----------------------------------------------------------------------
        CALL tri_linear_alloc(tb(ibl)%tele_eq,mv,1_i4,'teleeq',
     $                        (/'teleeq'/))
        CALL tri_linear_alloc(tb(ibl)%tion_eq,mv,1_i4,'tioneq',
     $                        (/'tioneq'/))
        tb(ibl)%tele_eq%fs(1,:,:)=tb(ibl)%prese_eq%fs(1,:,:)/
     $    (kboltz*tb(ibl)%nd_eq%fs(1,:,:))
        tb(ibl)%tion_eq%fs(1,:,:)=
     $    (tb(ibl)%pres_eq%fs(1,:,:)-tb(ibl)%prese_eq%fs(1,:,:))*zeff/
     $    (kboltz*tb(ibl)%nd_eq%fs(1,:,:))
c-----------------------------------------------------------------------
c       allocate storage for rhs and solution of equations.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(rhs(ibl),1_i4,mv,0_i4,3_i4)
        CALL vector_type_alloc(cell_rhs(ibl),0_i4,mc,1_i4,3_i4)
        CALL vector_type_alloc(cell_crhs(ibl),0_i4,mc,1_i4,4_i4,nmodes)
c-----------------------------------------------------------------------
c       allocate storage for cell volume and centroid.
c-----------------------------------------------------------------------
        ALLOCATE(tb(ibl)%cell_vol(mc))
        ALLOCATE(tb(ibl)%r_cent(mc))
c-----------------------------------------------------------------------
c       storage for viscous heating density and hyper-resistive heating.
c-----------------------------------------------------------------------
        IF (visc_heat)
     $    CALL tblock_qp_alloc(tb(ibl)%qvisc,tb(ibl),1_i4,nmodes)
        IF (ohm_heat.AND.(hyp_eta>0._r8.OR.hyp_dbd>0._r8))
     $    CALL tblock_qp_alloc(tb(ibl)%qhyph,tb(ibl),1_i4,nmodes)
      ENDDO tblocks
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE variable_alloc
c-----------------------------------------------------------------------
c     subprogram 4. mass_mat_init.
c     generate the mass matrix that is used for the jfromb matrix.
c-----------------------------------------------------------------------
      SUBROUTINE mass_mat_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE matrix_mod
      USE matrix_storage_mod
      USE integrands
      USE regularity
      USE input
      IMPLICIT NONE

      INTEGER(i4) :: ibl,iv,ivp,ibe,ix,iy,iq,is,mxb,myb,mv
c-----------------------------------------------------------------------
c     for toroidal simulations that include the geometric axis within
c     the computational domain, we also need a mass matrix with 
c     Dirichlet boundaries applied to the R=0 border.
c-----------------------------------------------------------------------
      ALLOCATE(mass_mat(1))
      mass_mat(1)%fcomp=0
      mass_mat(1)%foff=0
      mass_mat(1)%eliminated=.false.
      mass_mat(1)%symmetric=.true.
c-----------------------------------------------------------------------
c     create the mass matrix.
c-----------------------------------------------------------------------
      DO iq=1,SIZE(mass_mat)
        mass_mat(iq)%nqty=1
        mass_mat(iq)%nqdis=0
        mass_mat(iq)%diag_scale=1._r8
        mass_mat(iq)%essential_cond="none"
        ALLOCATE(mass_mat(iq)%vcomp(1))
        mass_mat(iq)%vcomp(1)='s'
        ALLOCATE(mass_mat(iq)%rbl_mat(nrbl))
        DO ibl=1,nrbl
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          CALL matrix_rbl_alloc(mass_mat(iq)%rbl_mat(ibl),mxb,myb,
     $                          1_i4,poly_degree,0_i4,0_i4)
          CALL rblock_make_matrix(rb(ibl),mass_mat(iq)%rbl_mat(ibl),
     $                            get_mass,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       tblocks:
c-----------------------------------------------------------------------
        ALLOCATE(mass_mat(iq)%tbl_mat(nrbl+1:nbl))
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_alloc(mass_mat(iq)%tbl_mat(ibl),tb(ibl)%tgeom,
     $                          1_i4)
          CALL tblock_make_matrix(tb(ibl),mass_mat(iq)%
     $                            tbl_mat(ibl)%lmat,get_mass,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       apply the dirichlet r=0 conditions if needed.  (not used)
c-----------------------------------------------------------------------
        IF (iq==2)  CALL regular_op(mass_mat(2))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mass_mat_init
c-----------------------------------------------------------------------
c     subprogram 5. quadrature_save.
c     evaluate equilibrium quantities and needed derivatives at the
c     Gaussian quadrature points for efficiency during finite element
c     computations.  also evaluate the initial conditions.
c-----------------------------------------------------------------------
      SUBROUTINE quadrature_save
      USE local
      USE input
      USE fields
      USE rblock
      USE tblock
      USE global
      USE input
      USE physdat
      USE math_tran
      IMPLICIT NONE
      
      INTEGER(i4) :: ibl,ncx,ncy,ng,ig
      REAL(r8) :: hfac,kcrfac
      REAL(r8), DIMENSION(:,:,:), POINTER :: nd_eqr,nd_eqz,te_eq,
     $          b_eq,j_eq,jaeq_r,jaeq_z,v_eq,veq_r,veq_z,dvveq,eq_force
      LOGICAL :: booltmp,parvtdep
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
c     these vector structures are used to determine whether matrices
c     need updating.
c-----------------------------------------------------------------------
      ALLOCATE(elecd_n0(nbl))
      ALLOCATE(kappli_n0(nbl))
      ALLOCATE(kapple_n0(nbl))
c-----------------------------------------------------------------------
c     loop over blocks and evaluate quantities needed in the finite
c     element integrations.  this will need to be modified if new
c     computations require more derivatives.
c
c     change equilibrium be_phi and ja_phi to cylindrical.
c
c     for equilibrium flow, it is more useful to save grad(V_eq) than
c     the separate R- and Z- derivatives.
c
c     compute equilibrium temperatures at the quadrature points.
c-----------------------------------------------------------------------
      hfac=(1._r8-meomi)/(elementary_q*(1._r8+meomi))
      parvtdep=(par_visc>0.AND.parvisc_model=='plltdep'.AND.beta>0)
      DO ibl=1,nrbl
        ncx=rb(ibl)%mx
        ncy=rb(ibl)%my
        ng=rb(ibl)%ng

        CALL rblock_qp_alloc(rb(ibl)%qbe_eq,rb(ibl),3_i4)
        CALL rblock_qp_update(rb(ibl)%be_eq,rb(ibl)%qbe_eq,rb(ibl))
        IF (geom=='tor') THEN
          rb(ibl)%qbe_eq%qpf(3,:,:)=rb(ibl)%qbe_eq%qpf(3,:,:)
     $                             /rb(ibl)%bigr
          rb(ibl)%qbe_eq%qpfr(3,:,:)=
     $      (rb(ibl)%qbe_eq%qpfr(3,:,:)-rb(ibl)%qbe_eq%qpf(3,:,:))
     $                                 /rb(ibl)%bigr
          rb(ibl)%qbe_eq%qpfz(3,:,:)=rb(ibl)%qbe_eq%qpfz(3,:,:)
     $                              /rb(ibl)%bigr
        ENDIF

        CALL rblock_qp_alloc(rb(ibl)%qja_eq,rb(ibl),3_i4)
        CALL rblock_qp_update(rb(ibl)%ja_eq,rb(ibl)%qja_eq,rb(ibl))
        IF (geom=='tor') THEN
          rb(ibl)%qja_eq%qpf(3,:,:)=rb(ibl)%qja_eq%qpf(3,:,:)
     $                             *rb(ibl)%bigr
          rb(ibl)%qja_eq%qpfr(3,:,:)=
     $       rb(ibl)%qja_eq%qpfr(3,:,:)*rb(ibl)%bigr
     $      +rb(ibl)%qja_eq%qpf (3,:,:)/rb(ibl)%bigr
          rb(ibl)%qja_eq%qpfz(3,:,:)=rb(ibl)%qja_eq%qpfz(3,:,:)
     $                              *rb(ibl)%bigr
        ENDIF
c-----------------------------------------------------------------------
c       the equilibrium current density may be recomputed from the
c       curl of the equilibrium magnetic field.
c-----------------------------------------------------------------------
        IF (tor_eqja_fe) THEN
          rb(ibl)%qja_eq%qpf(1,:,:)=rb(ibl)%qbe_eq%qpfz(3,:,:)/mu0
          IF (geom=='tor') THEN
            rb(ibl)%qja_eq%qpf(2,:,:)=
     $        -(rb(ibl)%qbe_eq%qpfr(3,:,:)+
     $          rb(ibl)%qbe_eq%qpf (3,:,:)/rb(ibl)%bigr)/mu0
          ELSE
            rb(ibl)%qja_eq%qpf(2,:,:)=
     $        -rb(ibl)%qbe_eq%qpfr(3,:,:)/mu0
          ENDIF
          rb(ibl)%qja_eq%qpf(3,:,:)=
     $       (rb(ibl)%qbe_eq%qpfr(2,:,:)-
     $        rb(ibl)%qbe_eq%qpfz(1,:,:))/mu0
        ENDIF

        ALLOCATE(rb(ibl)%qdvv_eq%qpf(1,ng,ncx*ncy))
        IF (eq_flow/='none') THEN
          CALL rblock_qp_alloc(rb(ibl)%qve_eq,rb(ibl),3_i4)
          CALL rblock_qp_update(rb(ibl)%ve_eq,rb(ibl)%qve_eq,rb(ibl))
          ALLOCATE(rb(ibl)%qgrdveq%qpf(9,ng,ncx*ncy))
          ALLOCATE(rb(ibl)%qpi_veq%qpf(9,ng,ncx*ncy))
          IF (par_visc>0) ALLOCATE(rb(ibl)%qpi_pareq%qpf(9,ng,ncx*ncy))
          IF (gyr_visc>0) ALLOCATE(rb(ibl)%qpi_gyreq%qpf(9,ng,ncx*ncy))
          rb(ibl)%qgrdveq%qpf(1,:,:)=rb(ibl)%qve_eq%qpfr(1,:,:)
          rb(ibl)%qgrdveq%qpf(2,:,:)=rb(ibl)%qve_eq%qpfz(1,:,:)
          rb(ibl)%qgrdveq%qpf(4,:,:)=rb(ibl)%qve_eq%qpfr(2,:,:)
          rb(ibl)%qgrdveq%qpf(5,:,:)=rb(ibl)%qve_eq%qpfz(2,:,:)
          rb(ibl)%qgrdveq%qpf(6,:,:)=0._r8
          rb(ibl)%qgrdveq%qpf(7,:,:)=rb(ibl)%qve_eq%qpfr(3,:,:)
          rb(ibl)%qgrdveq%qpf(8,:,:)=rb(ibl)%qve_eq%qpfz(3,:,:)
          IF (geom=='tor') THEN
            rb(ibl)%qgrdveq%qpf(3,:,:)=
     $        -rb(ibl)%qve_eq%qpf(3,:,:)/rb(ibl)%bigr
            rb(ibl)%qgrdveq%qpf(9,:,:)=
     $         rb(ibl)%qve_eq%qpf(1,:,:)/rb(ibl)%bigr
          ELSE
            rb(ibl)%qgrdveq%qpf(3,:,:)=0._r8
            rb(ibl)%qgrdveq%qpf(9,:,:)=0._r8
          ENDIF
          rb(ibl)%qdvv_eq%qpf(1,:,:)=
     $      SUM(rb(ibl)%qgrdveq%qpf(1:9:4,:,:),1)
          DEALLOCATE(rb(ibl)%qve_eq%qpfr,rb(ibl)%qve_eq%qpfz)
        ELSE
          rb(ibl)%qdvv_eq%qpf=0._r8
        ENDIF

        CALL rblock_qp_alloc(rb(ibl)%qpres_eq,rb(ibl),1_i4)
        CALL rblock_qp_update(rb(ibl)%pres_eq,rb(ibl)%qpres_eq,rb(ibl))
        CALL rblock_qp_alloc(rb(ibl)%qprese_eq,rb(ibl),1_i4)
        CALL rblock_qp_alloc(rb(ibl)%qdiff_shape,rb(ibl),1_i4)
c-----------------------------------------------------------------------
c       when Jeq is recomputed, have the perpendicular component satisfy
c       static force-balance at the quadrature points exactly.
c       the qdiff_shape structure is used for temporary storage of
c       Beq**2, and the qprese_eq structure is used for temporary
c       storage of Jeq.Beq/Beq**2.
c-----------------------------------------------------------------------
        IF (tor_eqja_fe) THEN
          rb(ibl)%qdiff_shape%qpf(1,:,:)=
     $      MAX(SUM(rb(ibl)%qbe_eq%qpf**2,1),smallnum)
          rb(ibl)%qprese_eq%qpf(1,:,:)=
     $      SUM(rb(ibl)%qbe_eq%qpf*rb(ibl)%qja_eq%qpf,1)/
     $      rb(ibl)%qdiff_shape%qpf(1,:,:)
          rb(ibl)%qja_eq%qpf(1,:,:)=-rb(ibl)%qbe_eq%qpf(3,:,:)*
     $                               rb(ibl)%qpres_eq%qpfz(1,:,:)/
     $                               rb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      rb(ibl)%qprese_eq%qpf(1,:,:)*rb(ibl)%qbe_eq%qpf(1,:,:)
          rb(ibl)%qja_eq%qpf(2,:,:)= rb(ibl)%qbe_eq%qpf(3,:,:)*
     $                               rb(ibl)%qpres_eq%qpfr(1,:,:)/
     $                               rb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      rb(ibl)%qprese_eq%qpf(1,:,:)*rb(ibl)%qbe_eq%qpf(2,:,:)
          rb(ibl)%qja_eq%qpf(3,:,:)=(rb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               rb(ibl)%qpres_eq%qpfz(1,:,:)-
     $                               rb(ibl)%qbe_eq%qpf(2,:,:)*
     $                               rb(ibl)%qpres_eq%qpfr(1,:,:))/
     $                               rb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      rb(ibl)%qprese_eq%qpf(1,:,:)*rb(ibl)%qbe_eq%qpf(3,:,:)
        ENDIF

        CALL rblock_qp_update(rb(ibl)%prese_eq,rb(ibl)%qprese_eq,
     $                        rb(ibl))
        CALL rblock_qp_alloc(rb(ibl)%qnd_eq,rb(ibl),1_i4)
        CALL rblock_qp_update(rb(ibl)%nd_eq,rb(ibl)%qnd_eq,rb(ibl))
        CALL rblock_qp_update(rb(ibl)%diff_shape,rb(ibl)%qdiff_shape,
     $                        rb(ibl))
c-----------------------------------------------------------------------
c       compute equilibrium temperatures at the quadrature points.
c-----------------------------------------------------------------------
        CALL rblock_qp_alloc(rb(ibl)%qtele_eq,rb(ibl),1_i4)
        CALL rblock_qp_alloc(rb(ibl)%qtion_eq,rb(ibl),1_i4)
        CALL rblock_qp_alloc(rb(ibl)%qte_n0,rb(ibl),1_i4)
        CALL rblock_qp_alloc(rb(ibl)%qti_n0,rb(ibl),1_i4)
        rb(ibl)%qtion_eq%qpf=MAX(smallnum,zeff*
     $    (rb(ibl)%qpres_eq%qpf-rb(ibl)%qprese_eq%qpf)/
     $    (kboltz*rb(ibl)%qnd_eq%qpf))
        rb(ibl)%qtion_eq%qpfr=(zeff*
     $    (rb(ibl)%qpres_eq%qpfr-rb(ibl)%qprese_eq%qpfr)/kboltz-
     $     rb(ibl)%qtion_eq%qpf*rb(ibl)%qnd_eq%qpfr)/rb(ibl)%qnd_eq%qpf
        rb(ibl)%qtion_eq%qpfz=(zeff*
     $    (rb(ibl)%qpres_eq%qpfz-rb(ibl)%qprese_eq%qpfz)/kboltz-
     $     rb(ibl)%qtion_eq%qpf*rb(ibl)%qnd_eq%qpfz)/rb(ibl)%qnd_eq%qpf
        rb(ibl)%qtele_eq%qpf=MAX(smallnum,
     $    rb(ibl)%qprese_eq%qpf/(kboltz*rb(ibl)%qnd_eq%qpf))
        rb(ibl)%qtele_eq%qpfr=
     $    (rb(ibl)%qprese_eq%qpfr/kboltz-
     $     rb(ibl)%qtele_eq%qpf*rb(ibl)%qnd_eq%qpfr)/rb(ibl)%qnd_eq%qpf
        rb(ibl)%qtele_eq%qpfz=
     $    (rb(ibl)%qprese_eq%qpfz/kboltz-
     $     rb(ibl)%qtele_eq%qpf*rb(ibl)%qnd_eq%qpfz)/rb(ibl)%qnd_eq%qpf
c-----------------------------------------------------------------------
c       quadrature point storage and initial values for perturbed
c       fields.  start with magnetic field and current density.
c-----------------------------------------------------------------------
        IF (ohms=='2fl'.AND.advect=='all'.OR.
     $      separate_pe.AND.nonlinear)
     $    CALL rblock_qp_update(rb(ibl)%ja,rb(ibl)%qja,rb(ibl))
        IF (nonlinear)
     $    ALLOCATE(rb(ibl)%qbe_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.(impladv.OR.eta_model=="chodura".OR.
     $                     siop_type=="3D"))
     $    ALLOCATE(rb(ibl)%qja_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.gyr_visc>0)
     $    ALLOCATE(rb(ibl)%qti_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
c-----------------------------------------------------------------------
c       space for flow velocity and grad(V) for implicit advection.
c-----------------------------------------------------------------------
        CALL rblock_qp_update(rb(ibl)%ve,rb(ibl)%qve,rb(ibl))
        IF (nonlinear) THEN
          ALLOCATE(rb(ibl)%qve_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
          ALLOCATE(rb(ibl)%qgrdv%qpf(9_i4,mpsq_block(ibl),2**lphi))
        ENDIF
c-----------------------------------------------------------------------
c       pressures and number density.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          CALL rblock_qp_update(rb(ibl)%pres,rb(ibl)%qpres,rb(ibl))
          CALL rblock_qp_update(rb(ibl)%prese,rb(ibl)%qprese,rb(ibl))
          CALL rblock_qp_update(rb(ibl)%tion,rb(ibl)%qtion,rb(ibl))
          CALL rblock_qp_update(rb(ibl)%tele,rb(ibl)%qtele,rb(ibl))
          IF (nonlinear) THEN
            CALL qp0_bcast(rb(ibl)%qtion%qpf,rb(ibl)%qti_n0%qpf,nlayers)
            CALL qp0_bcast(rb(ibl)%qtele%qpf,rb(ibl)%qte_n0%qpf,nlayers)
            rb(ibl)%qti_n0%qpf=rb(ibl)%qti_n0%qpf+rb(ibl)%qtion_eq%qpf
            rb(ibl)%qte_n0%qpf=rb(ibl)%qte_n0%qpf+rb(ibl)%qtele_eq%qpf
          ELSE
            rb(ibl)%qti_n0%qpf=rb(ibl)%qtion_eq%qpf
            rb(ibl)%qte_n0%qpf=rb(ibl)%qtele_eq%qpf
          ENDIF
        ENDIF
        IF (continuity/='none') THEN
          CALL rblock_qp_update(rb(ibl)%nd,rb(ibl)%qnd,rb(ibl))
          IF (nonlinear) THEN
            ALLOCATE(rb(ibl)%qnd_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
            CALL qp_fft_save(rb(ibl)%qnd%qpf,rb(ibl)%qnd_tot%qpf,
     $                       rb(ibl)%mx,rb(ibl)%my,mpsq_block(ibl),
     $                       1_i4,rb(ibl)%ng,rb(ibl)%qnd_eq%qpf)
            IF (nd_diff>0.AND.nd_floor>0) THEN
              CALL rblock_qp_alloc(rb(ibl)%qdart,rb(ibl),1_i4)
              rb(ibl)%qdart%qpf=0._r8
            ENDIF
            IF (impladv.AND.nd_dart_upw>0) THEN
              ALLOCATE(rb(ibl)%qupw_phi%qpf(1_i4,mpsq_block(ibl),
     $                                      2**lphi))
              CALL rblock_qp_alloc(rb(ibl)%qupw_n0,rb(ibl),1_i4)
              CALL rblock_qp_alloc(rb(ibl)%qvv,rb(ibl),6_i4)
            ENDIF
          ENDIF
          IF ((nd_diff>0.OR.nd_hypd>0).AND.nd_correrr) THEN
            CALL rblock_qp_alloc(rb(ibl)%qndiff,rb(ibl),1_i4,nmodes)
            CALL rblock_qp_alloc(rb(ibl)%qndiffa,rb(ibl),1_i4,nmodes)
          ENDIF
          IF (nonlinear) THEN
            ALLOCATE(rb(ibl)%qndiff_phi%qpf
     $        (1_i4,mpsq_block(ibl),2**lphi))
            rb(ibl)%qndiff_phi%qpf=0._r8
            CALL rblock_qp_alloc(rb(ibl)%qndiff_n0,rb(ibl),1_i4)
            rb(ibl)%qndiff_n0%qpf=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       storage for factors in 3D semi-implicit operator.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.(siop_type=='3D')) THEN
          ALLOCATE(rb(ibl)%qgrdb%qpf(9_i4,mpsq_block(ibl),2**lphi))
          IF (beta>0) THEN
            ALLOCATE(rb(ibl)%qgrdp%qpf(3_i4,mpsq_block(ibl),2**lphi))
            ALLOCATE(rb(ibl)%qpr_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       coefficients for thermal conductivities and resistivity.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          IF ((p_model(1:5)=='aniso'.OR.parvtdep).AND.nonlinear) THEN
            IF (p_model=='aniso_plltdep'.OR.p_model=='aniso_tdep'.OR.
     $          parvtdep) THEN
              CALL vector_type_alloc(kappli_n0(ibl),0_i4,ncx,ncy,1_i4)
              CALL rblock_qp_alloc(rb(ibl)%qkappli_n0,rb(ibl),1_i4)
              CALL rblock_qp_alloc(rb(ibl)%qkaprpi_n0,rb(ibl),1_i4)
              rb(ibl)%qkaprpi_n0%qpf=k_perpi
              ALLOCATE(rb(ibl)%qkappli_phi%qpf(1_i4,mpsq_block(ibl),
     $                                         2**lphi))
              IF (.NOT.closure_n0_only)
     $          ALLOCATE(rb(ibl)%qkaprpi_phi%qpf(1_i4,mpsq_block(ibl),
     $                                           2**lphi))
              IF (separate_pe) THEN
                CALL vector_type_alloc(kapple_n0(ibl),0_i4,ncx,ncy,1_i4)
                CALL rblock_qp_alloc(rb(ibl)%qkapple_n0,rb(ibl),1_i4)
                CALL rblock_qp_alloc(rb(ibl)%qkaprpe_n0,rb(ibl),1_i4)
                rb(ibl)%qkaprpe_n0%qpf=k_perpe
                ALLOCATE(rb(ibl)%qkapple_phi%qpf(1_i4,mpsq_block(ibl),
     $                                           2**lphi))
                IF (.NOT.closure_n0_only)
     $            ALLOCATE(rb(ibl)%qkaprpe_phi%qpf(1_i4,mpsq_block(ibl),
     $                                             2**lphi))
              ENDIF
            ENDIF
          ENDIF
          IF (p_model(1:5)=='aniso'.AND.separate_pe.AND.k_cross>0) THEN
            kcrfac=2.5_r8*k_cross*gamm1*kboltz/elementary_q
            ALLOCATE(rb(ibl)%qte_b2%qpf(1_i4,ng,ncx*ncy))
            ALLOCATE(rb(ibl)%qti_b2%qpf(1_i4,ng,ncx*ncy))
            ALLOCATE(rb(ibl)%qbcrgte%qpf(3_i4,ng,ncx*ncy))
            ALLOCATE(rb(ibl)%qbcrgti%qpf(3_i4,ng,ncx*ncy))
            rb(ibl)%qte_b2%qpf(1,:,:)=kcrfac*
     $        rb(ibl)%qtele_eq%qpf(1,:,:)/SUM(rb(ibl)%qbe_eq%qpf**2,1)
            rb(ibl)%qti_b2%qpf(1,:,:)=kcrfac/zeff*
     $        rb(ibl)%qtion_eq%qpf(1,:,:)/SUM(rb(ibl)%qbe_eq%qpf**2,1)
            rb(ibl)%qbcrgte%qpf(1,:,:)=
     $       -rb(ibl)%qbe_eq%qpf(3,:,:)*rb(ibl)%qtele_eq%qpfz(1,:,:)
            rb(ibl)%qbcrgte%qpf(2,:,:)=
     $        rb(ibl)%qbe_eq%qpf(3,:,:)*rb(ibl)%qtele_eq%qpfr(1,:,:)
            rb(ibl)%qbcrgte%qpf(3,:,:)=
     $        (rb(ibl)%qbe_eq%qpf(1,:,:)*
     $         rb(ibl)%qtele_eq%qpfz(1,:,:)
     $        -rb(ibl)%qbe_eq%qpf(2,:,:)*
     $         rb(ibl)%qtele_eq%qpfr(1,:,:))
            rb(ibl)%qbcrgti%qpf(1,:,:)=
     $       -rb(ibl)%qbe_eq%qpf(3,:,:)*rb(ibl)%qtion_eq%qpfz(1,:,:)
            rb(ibl)%qbcrgti%qpf(2,:,:)=
     $        rb(ibl)%qbe_eq%qpf(3,:,:)*rb(ibl)%qtion_eq%qpfr(1,:,:)
            rb(ibl)%qbcrgti%qpf(3,:,:)=
     $        (rb(ibl)%qbe_eq%qpf(1,:,:)*
     $         rb(ibl)%qtion_eq%qpfz(1,:,:)
     $        -rb(ibl)%qbe_eq%qpf(2,:,:)*
     $         rb(ibl)%qtion_eq%qpfr(1,:,:))
          ENDIF
        ENDIF
        IF (eta_model=='eta n=0 only'.OR.eta_model=='eta full'.OR.
     $      eta_model=='chodura') THEN
          CALL vector_type_alloc(elecd_n0(ibl),0_i4,ncx,ncy,1_i4)
          CALL rblock_qp_alloc(rb(ibl)%qelecd_n0,rb(ibl),1_i4)
          CALL rblock_qp_alloc(rb(ibl)%qelecd_eq,rb(ibl),1_i4)
          rb(ibl)%qelecd_eq%qpf=MAX( elecd_min, MIN( elecd_max,
     $      elecd*(eta_ref_t/MAX(smallnum,rb(ibl)%qtele_eq%qpf))**1.5))
          IF (eta_model=='chodura') THEN  !  phenomenological Chodura
            rb(ibl)%qelecd_eq%qpf(1,:,:)=
     $        rb(ibl)%qelecd_eq%qpf(1,:,:)+
     $        elecd_chodura*SQRT(ndens/rb(ibl)%qnd_eq%qpf(1,:,:))*
     $        (1._r8-EXP(-f_chodura*
     $          SQRT(SUM(rb(ibl)%qja_eq%qpf**2,1)*mtot/
     $            (gamma*MAX(rb(ibl)%qpres_eq%qpf(1,:,:),smallnum)*
     $                   rb(ibl)%qnd_eq%qpf(1,:,:)))/elementary_q))
          ENDIF
          IF (.NOT.nonlinear)
     $      rb(ibl)%qelecd_n0%qpf=rb(ibl)%qelecd_eq%qpf
        ENDIF
        IF (threedeta)
     $    ALLOCATE(rb(ibl)%qelecd_phi%qpf(1_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
          ALLOCATE(rb(ibl)%qupti_phi%qpf(1_i4,mpsq_block(ibl),2**lphi))
          CALL rblock_qp_alloc(rb(ibl)%qupti_n0,rb(ibl),1_i4)
          IF (separate_pe) THEN
            ALLOCATE(rb(ibl)%qupte_phi%qpf(1_i4,mpsq_block(ibl),
     $                                     2**lphi))
            CALL rblock_qp_alloc(rb(ibl)%qupte_n0,rb(ibl),1_i4)
          ENDIF
          IF (continuity=='none'.OR.nd_dart_upw<=0)
     $      CALL rblock_qp_alloc(rb(ibl)%qvv,rb(ibl),6_i4)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     same for tblocks.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        ncx=tb(ibl)%mcell
        ncy=1_i4
        ng=tb(ibl)%ng

        CALL tblock_qp_alloc(tb(ibl)%qbe_eq,tb(ibl),3_i4)
        CALL tblock_qp_update(tb(ibl)%be_eq,tb(ibl)%qbe_eq,tb(ibl))
        IF (geom=='tor') THEN
          tb(ibl)%qbe_eq%qpf(3,:,:)=tb(ibl)%qbe_eq%qpf(3,:,:)
     $                             /tb(ibl)%tgeom%bigr
          tb(ibl)%qbe_eq%qpfr(3,:,:)=
     $      (tb(ibl)%qbe_eq%qpfr(3,:,:)-tb(ibl)%qbe_eq%qpf(3,:,:))
     $                                 /tb(ibl)%tgeom%bigr
          tb(ibl)%qbe_eq%qpfz(3,:,:)=tb(ibl)%qbe_eq%qpfz(3,:,:)
     $                                /tb(ibl)%tgeom%bigr
        ENDIF
        CALL tblock_qp_alloc(tb(ibl)%qja_eq,tb(ibl),3_i4)
        CALL tblock_qp_update(tb(ibl)%ja_eq,tb(ibl)%qja_eq,tb(ibl))
        IF (geom=='tor') THEN
          tb(ibl)%qja_eq%qpf(3,:,:)=tb(ibl)%qja_eq%qpf(3,:,:)
     $                             *tb(ibl)%tgeom%bigr
          tb(ibl)%qja_eq%qpfr(3,:,:)=
     $       tb(ibl)%qja_eq%qpfr(3,:,:)*tb(ibl)%tgeom%bigr
     $      +tb(ibl)%qja_eq%qpf (3,:,:)/tb(ibl)%tgeom%bigr
          tb(ibl)%qja_eq%qpfz(3,:,:)=tb(ibl)%qja_eq%qpfz(3,:,:)
     $                              *tb(ibl)%tgeom%bigr
        ENDIF
c-----------------------------------------------------------------------
c       the equilibrium current density may be recomputed from the
c       curl of the equilibrium magnetic field.
c-----------------------------------------------------------------------
        IF (tor_eqja_fe) THEN
          tb(ibl)%qja_eq%qpf(1,:,:)=tb(ibl)%qbe_eq%qpfz(3,:,:)/mu0
          IF (geom=='tor') THEN
            tb(ibl)%qja_eq%qpf(2,:,:)=
     $        -(tb(ibl)%qbe_eq%qpfr(3,:,:)+
     $          tb(ibl)%qbe_eq%qpf (3,:,:)/tb(ibl)%tgeom%bigr)/mu0
          ELSE
            tb(ibl)%qja_eq%qpf(2,:,:)=
     $        -tb(ibl)%qbe_eq%qpfr(3,:,:)/mu0
          ENDIF
          tb(ibl)%qja_eq%qpf(3,:,:)=
     $       (tb(ibl)%qbe_eq%qpfr(2,:,:)-
     $        tb(ibl)%qbe_eq%qpfz(1,:,:))/mu0
        ENDIF


        ALLOCATE(tb(ibl)%qdvv_eq%qpf(1,ng,ncx*ncy))
        IF (eq_flow/='none') THEN
          CALL tblock_qp_alloc(tb(ibl)%qve_eq,tb(ibl),3_i4)
          CALL tblock_qp_update(tb(ibl)%ve_eq,tb(ibl)%qve_eq,tb(ibl))
          ALLOCATE(tb(ibl)%qgrdveq%qpf(9,ng,ncx*ncy))
          ALLOCATE(tb(ibl)%qpi_veq%qpf(9,ng,ncx*ncy))
          IF (par_visc>0) ALLOCATE(tb(ibl)%qpi_pareq%qpf(9,ng,ncx*ncy))
          IF (gyr_visc>0) ALLOCATE(tb(ibl)%qpi_gyreq%qpf(9,ng,ncx*ncy))
          tb(ibl)%qgrdveq%qpf(1,:,:)=tb(ibl)%qve_eq%qpfr(1,:,:)
          tb(ibl)%qgrdveq%qpf(2,:,:)=tb(ibl)%qve_eq%qpfz(1,:,:)
          tb(ibl)%qgrdveq%qpf(4,:,:)=tb(ibl)%qve_eq%qpfr(2,:,:)
          tb(ibl)%qgrdveq%qpf(5,:,:)=tb(ibl)%qve_eq%qpfz(2,:,:)
          tb(ibl)%qgrdveq%qpf(6,:,:)=0._r8
          tb(ibl)%qgrdveq%qpf(7,:,:)=tb(ibl)%qve_eq%qpfr(3,:,:)
          tb(ibl)%qgrdveq%qpf(8,:,:)=tb(ibl)%qve_eq%qpfz(3,:,:)
          IF (geom=='tor') THEN
            tb(ibl)%qgrdveq%qpf(3,:,:)=
     $        -tb(ibl)%qve_eq%qpf(3,:,:)/tb(ibl)%tgeom%bigr
            tb(ibl)%qgrdveq%qpf(9,:,:)=
     $         tb(ibl)%qve_eq%qpf(1,:,:)/tb(ibl)%tgeom%bigr
          ELSE
            tb(ibl)%qgrdveq%qpf(3,:,:)=0._r8
            tb(ibl)%qgrdveq%qpf(9,:,:)=0._r8
          ENDIF
          tb(ibl)%qdvv_eq%qpf(1,:,:)=
     $      SUM(tb(ibl)%qgrdveq%qpf(1:9:4,:,:),1)
          DEALLOCATE(tb(ibl)%qve_eq%qpfr,tb(ibl)%qve_eq%qpfz)
        ELSE
          tb(ibl)%qdvv_eq%qpf=0._r8
        ENDIF

        CALL tblock_qp_alloc(tb(ibl)%qpres_eq,tb(ibl),1_i4)
        CALL tblock_qp_update(tb(ibl)%pres_eq,tb(ibl)%qpres_eq,
     $                        tb(ibl))
        CALL tblock_qp_alloc(tb(ibl)%qprese_eq,tb(ibl),1_i4)
        CALL tblock_qp_alloc(tb(ibl)%qdiff_shape,tb(ibl),1_i4)
c-----------------------------------------------------------------------
c       when Jeq is recomputed, have the perpendicular component satisfy
c       static force-balance at the quadrature points exactly.
c       the qdiff_shape structure is used for temporary storage of
c       Beq**2, and the qprese_eq structure is used for temporary
c       storage of Jeq.Beq/Beq**2.
c-----------------------------------------------------------------------
        IF (tor_eqja_fe) THEN
          tb(ibl)%qdiff_shape%qpf(1,:,:)=
     $      SUM(tb(ibl)%qbe_eq%qpf**2,1)
          tb(ibl)%qprese_eq%qpf(1,:,:)=
     $      SUM(tb(ibl)%qbe_eq%qpf*tb(ibl)%qja_eq%qpf,1)/
     $      tb(ibl)%qdiff_shape%qpf(1,:,:)
          tb(ibl)%qja_eq%qpf(1,:,:)=-tb(ibl)%qbe_eq%qpf(3,:,:)*
     $                               tb(ibl)%qpres_eq%qpfz(1,:,:)/
     $                               tb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      tb(ibl)%qprese_eq%qpf(1,:,:)*tb(ibl)%qbe_eq%qpf(1,:,:)
          tb(ibl)%qja_eq%qpf(2,:,:)= tb(ibl)%qbe_eq%qpf(3,:,:)*
     $                               tb(ibl)%qpres_eq%qpfr(1,:,:)/
     $                               tb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      tb(ibl)%qprese_eq%qpf(1,:,:)*tb(ibl)%qbe_eq%qpf(2,:,:)
          tb(ibl)%qja_eq%qpf(3,:,:)=(tb(ibl)%qbe_eq%qpf(1,:,:)*
     $                               tb(ibl)%qpres_eq%qpfz(1,:,:)-
     $                               tb(ibl)%qbe_eq%qpf(2,:,:)*
     $                               tb(ibl)%qpres_eq%qpfr(1,:,:))/
     $                               tb(ibl)%qdiff_shape%qpf(1,:,:)+
     $      tb(ibl)%qprese_eq%qpf(1,:,:)*tb(ibl)%qbe_eq%qpf(3,:,:)
        ENDIF

        CALL tblock_qp_update(tb(ibl)%prese_eq,tb(ibl)%qprese_eq,
     $                        tb(ibl))
        CALL tblock_qp_alloc(tb(ibl)%qnd_eq,tb(ibl),1_i4)
        CALL tblock_qp_update(tb(ibl)%nd_eq,tb(ibl)%qnd_eq,tb(ibl))
        CALL tblock_qp_update(tb(ibl)%diff_shape,tb(ibl)%qdiff_shape,
     $                        tb(ibl))
c-----------------------------------------------------------------------
c       compute equilibrium temperatures at the quadrature points.
c-----------------------------------------------------------------------
        CALL tblock_qp_alloc(tb(ibl)%qtele_eq,tb(ibl),1_i4)
        CALL tblock_qp_alloc(tb(ibl)%qtion_eq,tb(ibl),1_i4)
        CALL tblock_qp_alloc(tb(ibl)%qte_n0,tb(ibl),1_i4)
        CALL tblock_qp_alloc(tb(ibl)%qti_n0,tb(ibl),1_i4)
        tb(ibl)%qtion_eq%qpf=MAX(smallnum,zeff*
     $    (tb(ibl)%qpres_eq%qpf-tb(ibl)%qprese_eq%qpf)/
     $    (kboltz*tb(ibl)%qnd_eq%qpf))
        tb(ibl)%qtion_eq%qpfr=(zeff*
     $    (tb(ibl)%qpres_eq%qpfr-tb(ibl)%qprese_eq%qpfr)/kboltz-
     $     tb(ibl)%qtion_eq%qpf*tb(ibl)%qnd_eq%qpfr)/tb(ibl)%qnd_eq%qpf
        tb(ibl)%qtion_eq%qpfz=(zeff*
     $    (tb(ibl)%qpres_eq%qpfz-tb(ibl)%qprese_eq%qpfz)/kboltz-
     $     tb(ibl)%qtion_eq%qpf*tb(ibl)%qnd_eq%qpfz)/tb(ibl)%qnd_eq%qpf
        tb(ibl)%qtele_eq%qpf=MAX(smallnum,
     $    tb(ibl)%qprese_eq%qpf/(kboltz*tb(ibl)%qnd_eq%qpf))
        tb(ibl)%qtele_eq%qpfr=
     $    (tb(ibl)%qprese_eq%qpfr/kboltz-
     $     tb(ibl)%qtele_eq%qpf*tb(ibl)%qnd_eq%qpfr)/tb(ibl)%qnd_eq%qpf
        tb(ibl)%qtele_eq%qpfz=
     $    (tb(ibl)%qprese_eq%qpfz/kboltz-
     $     tb(ibl)%qtele_eq%qpf*tb(ibl)%qnd_eq%qpfz)/tb(ibl)%qnd_eq%qpf
c-----------------------------------------------------------------------
c       storage at quadrature points for perturbed fields.
c       start with magnetic field and current density.
c-----------------------------------------------------------------------
        IF (ohms=='2fl'.AND.advect=='all'.OR.
     $      separate_pe.AND.nonlinear)
     $    CALL tblock_qp_update(tb(ibl)%ja,tb(ibl)%qja,tb(ibl))
        IF (nonlinear)
     $    ALLOCATE(tb(ibl)%qbe_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.(impladv.OR.eta_model=="chodura".OR.
     $                     siop_type=="3D"))
     $    ALLOCATE(tb(ibl)%qja_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.gyr_visc>0)
     $    ALLOCATE(tb(ibl)%qti_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
c-----------------------------------------------------------------------
c       space for flow velocity and grad(V) for implicit advection.
c-----------------------------------------------------------------------
        CALL tblock_qp_update(tb(ibl)%ve,tb(ibl)%qve,tb(ibl))
        IF (nonlinear) THEN
          ALLOCATE(tb(ibl)%qve_tot%qpf(3_i4,mpsq_block(ibl),2**lphi))
          ALLOCATE(tb(ibl)%qgrdv%qpf(9_i4,mpsq_block(ibl),2**lphi))
        ENDIF
c-----------------------------------------------------------------------
c       pressures and number density.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          CALL tblock_qp_update(tb(ibl)%pres,tb(ibl)%qpres,tb(ibl))
          CALL tblock_qp_update(tb(ibl)%prese,tb(ibl)%qprese,tb(ibl))
          CALL tblock_qp_update(tb(ibl)%tion,tb(ibl)%qtion,tb(ibl))
          CALL tblock_qp_update(tb(ibl)%tele,tb(ibl)%qtele,tb(ibl))
          IF (nonlinear) THEN
            CALL qp0_bcast(tb(ibl)%qtion%qpf,tb(ibl)%qti_n0%qpf,nlayers)
            CALL qp0_bcast(tb(ibl)%qtele%qpf,tb(ibl)%qte_n0%qpf,nlayers)
            tb(ibl)%qti_n0%qpf=tb(ibl)%qti_n0%qpf+tb(ibl)%qtion_eq%qpf
            tb(ibl)%qte_n0%qpf=tb(ibl)%qte_n0%qpf+tb(ibl)%qtele_eq%qpf
          ELSE
            tb(ibl)%qti_n0%qpf=tb(ibl)%qtion_eq%qpf
            tb(ibl)%qte_n0%qpf=tb(ibl)%qtele_eq%qpf
          ENDIF
        ENDIF
        IF (continuity/='none') THEN
          CALL tblock_qp_update(tb(ibl)%nd,tb(ibl)%qnd,tb(ibl))
          IF (nonlinear) THEN
            ALLOCATE(tb(ibl)%qnd_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
            CALL qp_fft_save(tb(ibl)%qnd%qpf,tb(ibl)%qnd_tot%qpf,
     $                       tb(ibl)%mcell,1_i4,mpsq_block(ibl),
     $                       1_i4,tb(ibl)%ng,tb(ibl)%qnd_eq%qpf)
            IF (nd_diff>0.AND.nd_floor>0) THEN
              CALL tblock_qp_alloc(tb(ibl)%qdart,tb(ibl),1_i4)
              tb(ibl)%qdart%qpf=0._r8
            ENDIF
            IF (impladv.AND.nd_dart_upw>0) THEN
              ALLOCATE(tb(ibl)%qupw_phi%qpf(1_i4,mpsq_block(ibl),
     $                                      2**lphi))
              CALL tblock_qp_alloc(tb(ibl)%qupw_n0,tb(ibl),1_i4)
              CALL tblock_qp_alloc(tb(ibl)%qvv,tb(ibl),6_i4)
            ENDIF
          ENDIF
          IF ((nd_diff>0.OR.nd_hypd>0).AND.nd_correrr) THEN
            CALL tblock_qp_alloc(tb(ibl)%qndiff,tb(ibl),1_i4,nmodes)
            CALL tblock_qp_alloc(tb(ibl)%qndiffa,tb(ibl),1_i4,nmodes)
          ENDIF
          IF (nonlinear) THEN
            ALLOCATE(tb(ibl)%qndiff_phi%qpf
     $        (1_i4,mpsq_block(ibl),2**lphi))
            tb(ibl)%qndiff_phi%qpf=0._r8
            CALL tblock_qp_alloc(tb(ibl)%qndiff_n0,tb(ibl),1_i4)
            tb(ibl)%qndiff_n0%qpf=0._r8
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       storage for factors in 3D semi-implicit operator.
c-----------------------------------------------------------------------
        IF (nonlinear.AND.(siop_type=='3D')) THEN
          ALLOCATE(tb(ibl)%qgrdb%qpf(9_i4,mpsq_block(ibl),2**lphi))
          IF (beta>0) THEN
            ALLOCATE(tb(ibl)%qgrdp%qpf(3_i4,mpsq_block(ibl),2**lphi))
            ALLOCATE(tb(ibl)%qpr_tot%qpf(1_i4,mpsq_block(ibl),2**lphi))
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c       coefficients for thermal conductivities and resistivity.
c-----------------------------------------------------------------------
        IF (beta>0) THEN
          IF ((p_model(1:5)=='aniso'.OR.parvtdep).AND.nonlinear) THEN
            IF (p_model=='aniso_plltdep'.OR.p_model=='aniso_tdep'.OR.
     $          parvtdep) THEN
              CALL vector_type_alloc(kappli_n0(ibl),0_i4,ncx,ncy,1_i4)
              CALL tblock_qp_alloc(tb(ibl)%qkappli_n0,tb(ibl),1_i4)
              CALL tblock_qp_alloc(tb(ibl)%qkaprpi_n0,tb(ibl),1_i4)
              tb(ibl)%qkaprpi_n0%qpf=k_perpi
              ALLOCATE(tb(ibl)%qkappli_phi%qpf(1_i4,mpsq_block(ibl),
     $                                         2**lphi))
              IF (.NOT.closure_n0_only)
     $          ALLOCATE(tb(ibl)%qkaprpi_phi%qpf(1_i4,mpsq_block(ibl),
     $                                           2**lphi))
              IF (separate_pe) THEN
                CALL vector_type_alloc(kapple_n0(ibl),0_i4,ncx,ncy,1_i4)
                CALL tblock_qp_alloc(tb(ibl)%qkapple_n0,tb(ibl),1_i4)
                CALL tblock_qp_alloc(tb(ibl)%qkaprpe_n0,tb(ibl),1_i4)
                tb(ibl)%qkaprpe_n0%qpf=k_perpe
                ALLOCATE(tb(ibl)%qkapple_phi%qpf(1_i4,mpsq_block(ibl),
     $                                           2**lphi))
                IF (.NOT.closure_n0_only)
     $            ALLOCATE(tb(ibl)%qkaprpe_phi%qpf(1_i4,mpsq_block(ibl),
     $                                             2**lphi))
              ENDIF
            ENDIF
          ENDIF
          IF (p_model(1:5)=='aniso'.AND.separate_pe.AND.k_cross>0) THEN
            ALLOCATE(tb(ibl)%qte_b2%qpf(1_i4,ng,ncx*ncy))
            ALLOCATE(tb(ibl)%qti_b2%qpf(1_i4,ng,ncx*ncy))
            ALLOCATE(tb(ibl)%qbcrgte%qpf(3_i4,ng,ncx*ncy))
            ALLOCATE(tb(ibl)%qbcrgti%qpf(3_i4,ng,ncx*ncy))
            tb(ibl)%qte_b2%qpf(1,:,:)=2.5_r8*k_cross*gamm1*
     $        tb(ibl)%qtele_eq%qpf(1,:,:)/SUM(tb(ibl)%qbe_eq%qpf**2,1)
            tb(ibl)%qti_b2%qpf(1,:,:)=2.5_r8*k_cross*gamm1/zeff*
     $        tb(ibl)%qtion_eq%qpf(1,:,:)/SUM(tb(ibl)%qbe_eq%qpf**2,1)
            tb(ibl)%qbcrgte%qpf(1,:,:)=
     $       -tb(ibl)%qbe_eq%qpf(3,:,:)*tb(ibl)%qtele_eq%qpfz(1,:,:)
            tb(ibl)%qbcrgte%qpf(2,:,:)=
     $        tb(ibl)%qbe_eq%qpf(3,:,:)*tb(ibl)%qtele_eq%qpfr(1,:,:)
            tb(ibl)%qbcrgte%qpf(3,:,:)=
     $        (tb(ibl)%qbe_eq%qpf(1,:,:)*
     $         tb(ibl)%qtele_eq%qpfz(1,:,:)
     $        -tb(ibl)%qbe_eq%qpf(2,:,:)*
     $         tb(ibl)%qtele_eq%qpfr(1,:,:))
            tb(ibl)%qbcrgti%qpf(1,:,:)=
     $       -tb(ibl)%qbe_eq%qpf(3,:,:)*tb(ibl)%qtion_eq%qpfz(1,:,:)
            tb(ibl)%qbcrgti%qpf(2,:,:)=
     $        tb(ibl)%qbe_eq%qpf(3,:,:)*tb(ibl)%qtion_eq%qpfr(1,:,:)
            tb(ibl)%qbcrgti%qpf(3,:,:)=
     $        (tb(ibl)%qbe_eq%qpf(1,:,:)*
     $         tb(ibl)%qtion_eq%qpfz(1,:,:)
     $        -tb(ibl)%qbe_eq%qpf(2,:,:)*
     $         tb(ibl)%qtion_eq%qpfr(1,:,:))
          ENDIF
        ENDIF
        IF (eta_model=='eta n=0 only'.OR.eta_model=='eta full'.OR.
     $      eta_model=='chodura') THEN
          CALL vector_type_alloc(elecd_n0(ibl),0_i4,ncx,ncy,1_i4)
          CALL tblock_qp_alloc(tb(ibl)%qelecd_n0,tb(ibl),1_i4)
          CALL tblock_qp_alloc(tb(ibl)%qelecd_eq,tb(ibl),1_i4)
          tb(ibl)%qelecd_eq%qpf=MAX( elecd_min, MIN( elecd_max,
     $      elecd*(eta_ref_t/MAX(smallnum,tb(ibl)%qtele_eq%qpf))**1.5))
          IF (eta_model=='chodura') THEN  !  phenomenological Chodura
            tb(ibl)%qelecd_eq%qpf(1,:,:)=
     $        tb(ibl)%qelecd_eq%qpf(1,:,:)+
     $        elecd_chodura*SQRT(ndens/tb(ibl)%qnd_eq%qpf(1,:,:))*
     $        (1._r8-EXP(-f_chodura*
     $          SQRT(SUM(tb(ibl)%qja_eq%qpf**2,1)*mtot/
     $            (gamma*MAX(tb(ibl)%qpres_eq%qpf(1,:,:),smallnum)*
     $                   tb(ibl)%qnd_eq%qpf(1,:,:)))/elementary_q))
          ENDIF
          IF (.NOT.nonlinear)
     $      tb(ibl)%qelecd_n0%qpf=tb(ibl)%qelecd_eq%qpf
        ENDIF
        IF (threedeta)
     $    ALLOCATE(tb(ibl)%qelecd_phi%qpf(1_i4,mpsq_block(ibl),2**lphi))
        IF (nonlinear.AND.impladv.AND.t_dart_upw>0) THEN
          ALLOCATE(tb(ibl)%qupti_phi%qpf(1_i4,mpsq_block(ibl),2**lphi))
          CALL tblock_qp_alloc(tb(ibl)%qupti_n0,tb(ibl),1_i4)
          IF (separate_pe) THEN
            ALLOCATE(tb(ibl)%qupte_phi%qpf(1_i4,mpsq_block(ibl),
     $                                     2**lphi))
            CALL tblock_qp_alloc(tb(ibl)%qupte_n0,tb(ibl),1_i4)
          ENDIF
          IF (continuity=='none'.OR.nd_dart_upw<=0)
     $      CALL tblock_qp_alloc(tb(ibl)%qvv,tb(ibl),6_i4)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     call storage routines to set initial quad-point data for magnetic
c     field and flow velocity and grad(V) for implicit advection.
c     do the same for number density and temperatures but use the 
c     average flag to avoid nonlinear diffusivity computations; the
c     *_old field are set by this point.
c-----------------------------------------------------------------------
      CALL b_store('end')
      CALL vcom_store('standard')
      IF (continuity/='none') CALL n_store('average',booltmp)
      IF (beta>0) THEN
        CALL temp_store('ion ave',booltmp)
        CALL temp_store('ele ave',booltmp)
      ENDIF
c-----------------------------------------------------------------------
c     compute and save the equilibrium stress tensors.
c-----------------------------------------------------------------------
      CALL pieq_comp
c-----------------------------------------------------------------------
c     save the equilibrium force acting on electrons for the Hall term. 
c     note that the -kboltz*grad(Te_eq)/e part is not needed.
c-----------------------------------------------------------------------
      IF (ohms/='mhd') THEN
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            ncx=rb(ibl)%mx
            ncy=rb(ibl)%my
            ng=rb(ibl)%ng
            ALLOCATE(rb(ibl)%qeq_force%qpf(3,ng,ncx*ncy))
            eq_force=>rb(ibl)%qeq_force%qpf
            nd_eqr=>rb(ibl)%qnd_eq%qpfr
            nd_eqz=>rb(ibl)%qnd_eq%qpfz
            te_eq=>rb(ibl)%qtele_eq%qpf
            b_eq=>rb(ibl)%qbe_eq%qpf
            j_eq=>rb(ibl)%qja_eq%qpf
            jaeq_r=>rb(ibl)%qja_eq%qpfr
            jaeq_z=>rb(ibl)%qja_eq%qpfz
            v_eq=>rb(ibl)%qve_eq%qpf
            veq_r=>rb(ibl)%qve_eq%qpfr
            veq_z=>rb(ibl)%qve_eq%qpfz
            dvveq=>rb(ibl)%qdvv_eq%qpf
          ELSE
            ncx=tb(ibl)%mcell
            ncy=1_i4
            ng=tb(ibl)%ng
            ALLOCATE(tb(ibl)%qeq_force%qpf(3,ng,ncx*ncy))
            eq_force=>tb(ibl)%qeq_force%qpf
            nd_eqr=>rb(ibl)%qnd_eq%qpfr
            nd_eqz=>rb(ibl)%qnd_eq%qpfz
            te_eq=>rb(ibl)%qtele_eq%qpf
            b_eq=>tb(ibl)%qbe_eq%qpf
            j_eq=>tb(ibl)%qja_eq%qpf
            jaeq_r=>tb(ibl)%qja_eq%qpfr
            jaeq_z=>tb(ibl)%qja_eq%qpfz
            v_eq=>tb(ibl)%qve_eq%qpf
            veq_r=>tb(ibl)%qve_eq%qpfr
            veq_z=>tb(ibl)%qve_eq%qpfz
            dvveq=>tb(ibl)%qdvv_eq%qpf
          ENDIF
          CALL math_cart_cross(eq_force,j_eq,b_eq,hfac)
          IF (beta>0) THEN
            eq_force(1,:,:)=eq_force(1,:,:)-
     $        nd_eqr(1,:,:)*te_eq(1,:,:)*kboltz*hfac
            eq_force(2,:,:)=eq_force(2,:,:)-
     $        nd_eqz(1,:,:)*te_eq(1,:,:)*kboltz*hfac
          ENDIF
          IF (ohms=='2fl'.AND.advect=='all'.AND.eq_flow/='none') THEN
            eq_force(1,:,:)=eq_force(1,:,:)+
     $         (v_eq(1,:,:)*jaeq_r(1,:,:)+
     $          v_eq(2,:,:)*jaeq_z(1,:,:)+
     $          j_eq(1,:,:)*veq_r(1,:,:)+
     $          j_eq(2,:,:)*veq_z(1,:,:)-
     $          2._r8*v_eq(3,:,:)*j_eq(3,:,:)+
     $          dvveq(1,:,:)*j_eq(1,:,:))/elementary_q**2
            eq_force(2,:,:)=eq_force(2,:,:)+
     $         (v_eq(1,:,:)*jaeq_r(2,:,:)+
     $          v_eq(2,:,:)*jaeq_z(2,:,:)+
     $          j_eq(1,:,:)*veq_r(2,:,:)+
     $          j_eq(2,:,:)*veq_z(2,:,:)+
     $          dvveq(1,:,:)*j_eq(2,:,:))/elementary_q**2
            eq_force(3,:,:)=eq_force(3,:,:)+
     $         (v_eq(1,:,:)*jaeq_r(3,:,:)+
     $          v_eq(2,:,:)*jaeq_z(3,:,:)+
     $          j_eq(1,:,:)*veq_r(3,:,:)+
     $          j_eq(2,:,:)*veq_z(3,:,:)+
     $          v_eq(3,:,:)*j_eq(1,:,:)+
     $          v_eq(1,:,:)*j_eq(3,:,:)+
     $          dvveq(1,:,:)*j_eq(3,:,:))/elementary_q**2
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE quadrature_save
c-----------------------------------------------------------------------
c     subprogram 6. e_applied_init.
c     create a spatial function for applying an electric field to the
c     n=0 mode.  this is usually only used when an equilibrium field is
c     not available.
c-----------------------------------------------------------------------
      SUBROUTINE e_applied_init
      USE local
      USE fields
      USE seam_storage_mod
      USE input
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ibe,iv,ix,iy,ivp
c-----------------------------------------------------------------------
c     allocate arrays for an applied electric field.
c     this is treated as a bilinear quantity to prevent non-monotonic
c     profiles.
c-----------------------------------------------------------------------
      ALLOCATE(e_applied(nbl))
      DO ibl=1,nrbl
        CALL lagr_quad_alloc(rb(ibl)%e_applied,rb(ibl)%mx,rb(ibl)%my,
     $                       3_i4,1_i4,'ea',(/'e_appl'/))
        rb(ibl)%e_applied=0
        CALL vector_ptassign_laq2(e_applied(ibl),rb(ibl)%e_applied)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tri_linear_alloc(tb(ibl)%e_applied,tb(ibl)%mvert,
     $                        3_i4,'ea',(/'e_appl'/))
        tb(ibl)%e_applied=0
        CALL vector_ptassign_tl2(e_applied(ibl),tb(ibl)%e_applied)
      ENDDO
c-----------------------------------------------------------------------
c-PRE a block loop should specify the shape of e_applied here.
c     if used, it will need to be added to electric fields in
c     integrands, too.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE e_applied_init
c-----------------------------------------------------------------------
c     subprogram 7. boundary_vals_init.
c     ensure that the initial conditions satisfy appropriate boundary
c     conditions at the start of a run.
c-----------------------------------------------------------------------
      SUBROUTINE boundary_vals_init
      USE local
      USE fields
      USE seam_storage_mod
      USE input
      USE global
      USE boundary
      USE regularity
      IMPLICIT NONE
      
      INTEGER(i4) :: ibl,ibe
c-----------------------------------------------------------------------
c     interface block for surface_exb.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE surface_exb(vdum,vcent)
        USE local
        USE fields
        USE seam_storage_mod
        USE global
        USE input
        IMPLICIT NONE

        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vdum
        REAL(r8), INTENT(IN) :: vcent
        END SUBROUTINE surface_exb
      END INTERFACE
c-----------------------------------------------------------------------
c     loop over blocks, and set boundary conditions on B and V:
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        IF (zero_bnorm)
     $    CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
        IF (flow_bc=='no-slip') THEN
          CALL dirichlet_rhs(ve(ibe),seam(ibe),'all',3_i4)
        ELSE
          CALL dirichlet_rhs(ve(ibe),seam(ibe),'3vn',3_i4)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     loop over blocks, and set regularity conditions at R=0 points 
c     on B, V, and Ps:
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(r0block_list)
        ibe=r0block_list(ibl)
        CALL regular_vec(be(ibe),seam(ibe),'init_cyl_vec',3_i4,
     $                   nmodes,nindex)
        CALL regular_vec(ve(ibe),seam(ibe),'init_cyl_vec',3_i4,
     $                   nmodes,nindex)
        CALL regular_vec(pres(ibe),seam(ibe),'scalar',1_i4,
     $                   nmodes,nindex)
        CALL regular_vec(prese(ibe),seam(ibe),'scalar',1_i4,
     $                   nmodes,nindex)
        CALL regular_vec(tion(ibe),seam(ibe),'scalar',1_i4,
     $                   nmodes,nindex)
        CALL regular_vec(tele(ibe),seam(ibe),'scalar',1_i4,
     $                   nmodes,nindex)
      ENDDO
c-----------------------------------------------------------------------
c     also ensure that n=1 V,B_phi=i*V,B_r at R=0 points.
c-----------------------------------------------------------------------
      CALL regular_ave(be,3_i4,nmodes,nindex)
      CALL regular_ave(ve,3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     reset surface exb velocity if there is an applied voltage.
c-----------------------------------------------------------------------
      IF ((loop_volt/=0.OR.i_desired/=0.OR.e_vertical/=0).AND.
     $    norm_flow(1:3)=='exb') CALL surface_exb(ve,0._r8)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE boundary_vals_init
c-----------------------------------------------------------------------
c     subprogram 8. matrix_init.
c     allocates arrays for saving the matrix of each equation that
c     requires more than a mass matrix on its left side.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_init
      USE local
      USE matrix_storage_mod
      USE fields
      USE boundary
      USE input
      USE global
      USE iter_cg
      USE seam_storage_mod
      USE regularity
      IMPLICIT NONE

      INTEGER(i4) :: ibl,imode,iq,iv,nn,nregmin=1,nregmax=1,imat,jmat,
     $               jq,iq2,ndis,nqd,nqc
      CHARACTER(1), DIMENSION(3) :: vstand=(/'r','z','p'/)
c-----------------------------------------------------------------------
c     allocate global matrix structures:
c     if the effective impedance is 3D or the implicit operator is non-
c     Hermitian, we need a complex 2D matrix structure, otherwise
c     use real 3x3 vector blocks.
c
c     if poly_divb is non-negative, there is an auxiliary discontinuous
c     scalar field that is solved simultaneously with the advance of B.
c
c     if hyper-resistivity is used in without time-splitting, it is
c     solved simultaneously with the rest of the B-advance.  a larger
c     system is needed for the auxiliary vector.
c-----------------------------------------------------------------------
      IF (threedeta.OR.ohms/='mhd'.OR.impladv) THEN
        ALLOCATE(bhmhd_cmat(nmodes))
        ALLOCATE(bhmhd_cfac(nmodes))
        IF (poly_divb>=0) THEN
          nqd=1
          ndis=rb(1)%auxb%n_int
        ELSE
          nqd=0
          ndis=0
        ENDIF
        IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
          nqc=6
        ELSE
          nqc=3
        ENDIF
        DO imode=1,nmodes
          bhmhd_cmat(imode)%fcomp=nindex(imode)
          bhmhd_cmat(imode)%foff=0
          ALLOCATE(bhmhd_cmat(imode)%vcomp(nqc))
          bhmhd_cmat(imode)%vcomp(1:3)=vstand(1:3)
          IF (nqc==6) bhmhd_cmat(imode)%vcomp(4:6)=vstand(1:3)
          CALL comp_matrix_init_alloc(bhmhd_cmat(imode),nqc,
     $                                poly_degree,nqd,ndis)
          IF (ohms/='mhd'.OR.impladv)
     $        bhmhd_cmat(imode)%hermitian=.false.
          CALL iter_fac_alloc(bhmhd_cmat(imode),bhmhd_cfac(imode),
     $                        nqc,bmhd_solver,.true.)
        ENDDO
      ELSE
        ALLOCATE(bmhd_mat(nmodes))
        ALLOCATE(bmhd_fac(nmodes))
        DO imode=1,nmodes
          bmhd_mat(imode)%fcomp=nindex(imode)
          bmhd_mat(imode)%foff=0
          ALLOCATE(bmhd_mat(imode)%vcomp(3))
          bmhd_mat(imode)%vcomp=vstand(1:3)
          CALL matrix_init_alloc(bmhd_mat(imode),3_i4,poly_degree,
     $                           0_i4,0_i4)
          CALL iter_fac_alloc(bmhd_mat(imode),bmhd_fac(imode),
     $                        3_i4,bmhd_solver,.true.)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     if hyper-resistivity is used through a time-split advance,
c     it needs a separate matrix structure.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND.split_hypeta) THEN
        ALLOCATE(bhyp_cmat(nmodes))
        ALLOCATE(bhyp_cfac(nmodes))
        DO imode=1,nmodes
          bhyp_cmat(imode)%fcomp=nindex(imode)
          bhyp_cmat(imode)%foff=0
          ALLOCATE(bhyp_cmat(imode)%vcomp(6))
          bhyp_cmat(imode)%vcomp(1:3)=vstand(1:3)
          bhyp_cmat(imode)%vcomp(4:6)=vstand(1:3)
          CALL comp_matrix_init_alloc(bhyp_cmat(imode),6_i4,
     $                                poly_degree,0_i4,0_i4)
          bhyp_cmat(imode)%hermitian=.false.
          CALL iter_fac_alloc(bhyp_cmat(imode),bhyp_cfac(imode),
     $                        6_i4,bmhd_solver,.true.)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     divergence of b cleaner; factors are reduced to 2x2 for n=0:
c-----------------------------------------------------------------------
      IF (divbd>0.AND.split_divb) THEN
        ALLOCATE(divb_mat(nmodes))
        ALLOCATE(divb_fac(nmodes))
        DO imode=1,nmodes
          divb_mat(imode)%fcomp=nindex(imode)
          divb_mat(imode)%foff=0
          CALL matrix_init_alloc(divb_mat(imode),3_i4,poly_degree,
     $                           0_i4,0_i4)
          CALL iter_fac_alloc(divb_mat(imode),
     $                        divb_fac(imode),3_i4,solver,.true.)
          ALLOCATE(divb_mat(imode)%vcomp(3))
          divb_mat(imode)%vcomp=vstand(1:3)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     nodal current density (essentially a vector mass matrix with
c     possible coupling from regularity conditions):
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all'.OR.
     $    separate_pe.AND.nonlinear) THEN
        IF (any_r0blocks) THEN
          nregmin=3
          DO imode=1,nmodes
            SELECT CASE(nindex(imode))
            CASE(0)
              nregmin=1
            CASE(1)
              nregmin=MIN(nregmin,2)
              nregmax=MAX(nregmax,2)
            CASE DEFAULT
              nregmax=3
            END SELECT
          ENDDO
        ENDIF
        ALLOCATE(j_mat(nregmin:nregmax))
        ALLOCATE(j_fac(nregmin:nregmax))
        DO imode=nregmin,nregmax
          ALLOCATE(j_mat(imode)%vcomp(3))
          j_mat(imode)%vcomp=vstand(1:3)
          CALL matrix_init_alloc(j_mat(imode),3_i4,poly_degree,
     $                           0_i4,0_i4)
          j_mat(imode)%fcomp=imode-1
          j_mat(imode)%foff=0
c-----------------------------------------------------------------------
c         form:
c-----------------------------------------------------------------------
          DO ibl=1,nrbl
            DO imat=1,j_mat(imode)%rbl_mat(ibl)%nbtype
              DO jmat=1,j_mat(imode)%rbl_mat(ibl)%nbtype
                j_mat(imode)%rbl_mat(ibl)%mat(jmat,imat)%arr=0
              ENDDO
            ENDDO
          ENDDO
          DO ibl=nrbl+1,nbl
            DO iv=0,tb(ibl)%mvert
              j_mat(imode)%tbl_mat(ibl)%lmat(iv)%element=0
            ENDDO
          ENDDO
          CALL add_mass(j_mat(imode),3_i4)
          CALL regular_op(j_mat(imode))
        ENDDO
        DO ibl=1,nrbl
          CALL matrix_rbl_dealloc(mass_mat(1)%rbl_mat(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_dealloc(mass_mat(1)%tbl_mat(ibl))
        ENDDO
        DEALLOCATE(mass_mat)
c-----------------------------------------------------------------------
c       factor:
c-----------------------------------------------------------------------
        DO imode=nregmin,nregmax
          CALL iter_fac_alloc(j_mat(imode),j_fac(imode),3_i4,solver,
     $                        .false.)
          CALL iter_factor(j_mat(imode),j_fac(imode),3_i4,solver,
     $                     off_diag_fac)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     velocity advances:
c-----------------------------------------------------------------------
      IF (ohms/='hall') THEN
c-----------------------------------------------------------------------
c       time split viscous diffusion.
c-----------------------------------------------------------------------
        IF (split_visc) THEN
          ALLOCATE(visc_mat(nmodes))
          ALLOCATE(visc_fac(nmodes))
          DO imode=1,nmodes
            visc_mat(imode)%fcomp=nindex(imode)
            visc_mat(imode)%foff=0
            ALLOCATE(visc_mat(imode)%vcomp(3))
            visc_mat(imode)%vcomp=vstand(1:3)
            CALL matrix_init_alloc(visc_mat(imode),3_i4,poly_degree,
     $                             0_i4,0_i4)
            CALL iter_fac_alloc(visc_mat(imode),
     $                          visc_fac(imode),3_i4,solver,.true.)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       the isotropic semi-implicit operator for velocity has the same
c       form as the viscosity operator.  the anisotropic operator
c       couples all six (real and imag) vector components.
c
c       full continuity evolution requires the complex operator
c       regardless of mhd_si_iso and a separate operator with the semi-
c       implicit terms only to help speed matrix-vector product
c       computations.
c
c       if poly_divv is non-negative, there are auxiliary discontinuous
c       scalar fields that are solved simultaneously with the advance
c       of V.
c-----------------------------------------------------------------------
        IF (mhd_si_iso<1.OR.nonlinear.AND.continuity=='full'.OR.
     $      par_visc>0.OR.iso_visc>0.OR.impladv) THEN
          ALLOCATE(vmhd_cmat(nmodes))
          ALLOCATE(vmhd_cfac(nmodes))
          IF (poly_divv>=0.AND.nrbl>0) THEN
            nqd=2
            ndis=rb(1)%auxv%n_int
          ELSE
            nqd=0
            ndis=0
          ENDIF
          DO imode=1,nmodes
            vmhd_cmat(imode)%fcomp=nindex(imode)
            vmhd_cmat(imode)%foff=0
            ALLOCATE(vmhd_cmat(imode)%vcomp(3))
            vmhd_cmat(imode)%vcomp=vstand(1:3)
            CALL comp_matrix_init_alloc(vmhd_cmat(imode),3_i4,
     $                                  poly_degree,nqd,ndis)
            IF (impladv) vmhd_cmat(imode)%hermitian=.false.
            CALL iter_fac_alloc(vmhd_cmat(imode),vmhd_cfac(imode),
     $                          3_i4,vmhd_solver,.true.)
          ENDDO
        ELSE
          ALLOCATE(vmhd_mat(nmodes))
          ALLOCATE(vmhd_fac(nmodes))
          DO imode=1,nmodes
            vmhd_mat(imode)%fcomp=nindex(imode)
            vmhd_mat(imode)%foff=0
            ALLOCATE(vmhd_mat(imode)%vcomp(3))
            vmhd_mat(imode)%vcomp=vstand(1:3)
            CALL matrix_init_alloc(vmhd_mat(imode),3_i4,poly_degree,
     $                             0_i4,0_i4)
            CALL iter_fac_alloc(vmhd_mat(imode),vmhd_fac(imode),
     $                          3_i4,vmhd_solver,.true.)
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     temperature advance:
c-----------------------------------------------------------------------
      IF (beta>0) THEN
        ALLOCATE(ti_cmat(nmodes))
        ALLOCATE(ti_cfac(nmodes))
        DO imode=1,nmodes
          ti_cmat(imode)%fcomp=nindex(imode)
          ti_cmat(imode)%foff=0
          ALLOCATE(ti_cmat(imode)%vcomp(1))
          ti_cmat(imode)%vcomp(1)='s'
          CALL comp_matrix_init_alloc(ti_cmat(imode),1_i4,poly_degree,
     $                                0_i4,0_i4)
          IF (separate_pe.AND.k_cross>0.OR.impladv)
     $      ti_cmat(imode)%hermitian=.false.
          CALL iter_fac_alloc(ti_cmat(imode),ti_cfac(imode),
     $                        1_i4,temp_solver,.true.)
        ENDDO
        IF (separate_pe) THEN
          ALLOCATE(te_cmat(nmodes))
          ALLOCATE(te_cfac(nmodes))
          DO imode=1,nmodes
            te_cmat(imode)%fcomp=nindex(imode)
            te_cmat(imode)%foff=0
            ALLOCATE(te_cmat(imode)%vcomp(1))
            te_cmat(imode)%vcomp(1)='s'
            CALL comp_matrix_init_alloc(te_cmat(imode),1_i4,poly_degree,
     $                                  0_i4,0_i4)
            IF (k_cross>0.OR.impladv) te_cmat(imode)%hermitian=.false.
            CALL iter_fac_alloc(te_cmat(imode),te_cfac(imode),
     $                          1_i4,temp_solver,.true.)
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     particle density diffusion.
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        IF (impladv) THEN
          ALLOCATE(nd_cmat(nmodes))
          ALLOCATE(nd_cfac(nmodes))
          IF (nd_hypd>0) THEN
            jq=2
          ELSE
            jq=1
          ENDIF
          DO imode=1,nmodes
            nd_cmat(imode)%fcomp=nindex(imode)
            nd_cmat(imode)%foff=0
            ALLOCATE(nd_cmat(imode)%vcomp(jq))
            nd_cmat(imode)%vcomp(:)='s'
            CALL comp_matrix_init_alloc(nd_cmat(imode),jq,poly_degree,
     $                                  0_i4,0_i4)
            nd_cmat(imode)%hermitian=.false.
            CALL iter_fac_alloc(nd_cmat(imode),
     $                          nd_cfac(imode),jq,solver,.true.)
          ENDDO
        ELSE
          ALLOCATE(ndiso_mat(nmodes))
          ALLOCATE(ndiso_fac(nmodes))
          DO imode=1,nmodes
            ndiso_mat(imode)%fcomp=nindex(imode)
            ndiso_mat(imode)%foff=0
            ALLOCATE(ndiso_mat(imode)%vcomp(1))
            ndiso_mat(imode)%vcomp(1)='s'
            CALL matrix_init_alloc(ndiso_mat(imode),1_i4,poly_degree,
     $                             0_i4,0_i4)
            CALL iter_fac_alloc(ndiso_mat(imode),
     $                          ndiso_fac(imode),1_i4,solver,.true.)
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     if hyper-viscosity is used, the matrix initialization is done
c     in a separate routine.
c-----------------------------------------------------------------------
      IF (hyp_visc>0._r8) CALL hv_mat_init
c-----------------------------------------------------------------------
c     finally, create the weights along block borders that are used by
c     the iterative solver.
c-----------------------------------------------------------------------
      CALL border_weights(solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
c-----------------------------------------------------------------------
c     subprograms internal to matrix_init for convenience.
c-----------------------------------------------------------------------
        CONTAINS

c-----------------------------------------------------------------------
c       perform allocations needed for a matrix structure.
c-----------------------------------------------------------------------
        SUBROUTINE matrix_init_alloc(mat_str,mq,polyd,mqd,nbd)

        TYPE(global_matrix_type), INTENT(OUT) :: mat_str
        INTEGER(i4), INTENT(IN) :: mq,polyd,mqd,nbd

        ALLOCATE(mat_str%rbl_mat(nrbl))
        DO ibl=1,nrbl
          CALL matrix_rbl_alloc(mat_str%rbl_mat(ibl),rb(ibl)%mx,
     $                          rb(ibl)%my,mq,polyd,mqd,nbd)
        ENDDO
        ALLOCATE(mat_str%tbl_mat(nrbl+1:nbl))
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_alloc(mat_str%tbl_mat(ibl),tb(ibl)%tgeom,mq)
        ENDDO

        mat_str%nqty=mq
        mat_str%nqdis=mqd
        mat_str%eliminated=.false.
        mat_str%symmetric=.true.
        mat_str%diag_scale=1._r8
        mat_str%essential_cond="none"

        RETURN
        END SUBROUTINE matrix_init_alloc

c-----------------------------------------------------------------------
c       perform allocations needed for a complex matrix structure.
c-----------------------------------------------------------------------
        SUBROUTINE comp_matrix_init_alloc(mat_str,mq,polyd,mqd,nbd)

        TYPE(complex_matrix_type), INTENT(OUT) :: mat_str
        INTEGER(i4), INTENT(IN) :: mq,polyd,mqd,nbd

        ALLOCATE(mat_str%rbl_mat(nrbl))
        DO ibl=1,nrbl
          CALL matrix_rbl_alloc(mat_str%rbl_mat(ibl),rb(ibl)%mx,
     $                          rb(ibl)%my,mq,polyd,mqd,nbd)
        ENDDO
        ALLOCATE(mat_str%tbl_mat(nrbl+1:nbl))
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_alloc(mat_str%tbl_mat(ibl),tb(ibl)%tgeom,mq)
        ENDDO

        mat_str%nqty=mq
        mat_str%nqdis=mqd
        mat_str%eliminated=.false.
        mat_str%hermitian=.true.
        mat_str%diag_scale=1._r8
        mat_str%essential_cond="none"

        RETURN
        END SUBROUTINE comp_matrix_init_alloc

      END SUBROUTINE matrix_init
c-----------------------------------------------------------------------
c     subprogram 9. pointer_init.
c     allocate and assigns pointers to nodal data, which are independent
c     of block type.
c-----------------------------------------------------------------------
      SUBROUTINE pointer_init
      USE local
      USE fields
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ix,iy,ibase,m1,m2,ib1
      REAL(r8) :: dx,dy
c-----------------------------------------------------------------------
c     allocate block-wise arrays for all spatially dependent fields.
c-----------------------------------------------------------------------
      ALLOCATE(be(nbl),be_n0(nbl),ja(nbl),ve(nbl),pres(nbl),
     $  pres_n0(nbl),prese(nbl),nd(nbl),nd_n0(nbl),conc(nbl),
     $  be_eq(nbl),ja_eq(nbl),ve_eq(nbl),pres_eq(nbl),prese_eq(nbl),
     $  nd_eq(nbl),diff_shape(nbl),rwork1(nbl),work1(nbl),
     $  work2(nbl),work3(nbl),work4(nbl),work5(nbl),work6(nbl),
     $  tele(nbl),tion(nbl),tele_eq(nbl),tion_eq(nbl),ve_n0(nbl),
     $  auxb(nbl),auxv(nbl),w6v1(nbl),w6v2(nbl))
c-----------------------------------------------------------------------
c     make the pointer assignments to the rblock and tblock data.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL cvector_ptassign_laq(be(ibl),rb(ibl)%be)
        CALL cvector_ptassign_laq(ja(ibl),rb(ibl)%ja)
        CALL cvector_ptassign_laq(ve(ibl),rb(ibl)%ve)
        CALL cvector_ptassign_laq(pres(ibl),rb(ibl)%pres)
        CALL cvector_ptassign_laq(prese(ibl),rb(ibl)%prese)
        CALL cvector_ptassign_laq(nd(ibl),rb(ibl)%nd)
        CALL cvector_ptassign_laq(conc(ibl),rb(ibl)%conc)
        CALL cvector_ptassign_laq(tele(ibl),rb(ibl)%tele)
        CALL cvector_ptassign_laq(tion(ibl),rb(ibl)%tion)
        CALL cvector_ptassign_laq(work1(ibl),rb(ibl)%work1)
        CALL cvector_ptassign_laq(work2(ibl),rb(ibl)%work2)
        CALL cvector_ptassign_laq(work3(ibl),rb(ibl)%work3)
        CALL cvector_ptassign_laq(work4(ibl),rb(ibl)%work4)
        CALL cvector_ptassign_laq(work5(ibl),rb(ibl)%work5)
        CALL cvector_ptassign_laq(work6(ibl),rb(ibl)%work6)
        CALL cvector_ptassign_laq(w6v1(ibl),rb(ibl)%w6v1)
        CALL cvector_ptassign_laq(w6v2(ibl),rb(ibl)%w6v2)
        CALL cvector_ptassign_modq(auxb(ibl),rb(ibl)%auxb)
        CALL cvector_ptassign_modq(auxv(ibl),rb(ibl)%auxv)
        CALL vector_ptassign_laq2(be_n0(ibl),rb(ibl)%be_n0)
        CALL vector_ptassign_laq2(pres_n0(ibl),rb(ibl)%pres_n0)
        CALL vector_ptassign_laq2(nd_n0(ibl),rb(ibl)%nd_n0)
        CALL vector_ptassign_laq2(ve_n0(ibl),rb(ibl)%ve_n0)
        CALL vector_ptassign_laq2(rwork1(ibl),rb(ibl)%rwork1)
        CALL vector_ptassign_laq2(be_eq(ibl),rb(ibl)%be_eq)
        CALL vector_ptassign_laq2(ja_eq(ibl),rb(ibl)%ja_eq)
        CALL vector_ptassign_laq2(ve_eq(ibl),rb(ibl)%ve_eq)
        CALL vector_ptassign_laq2(pres_eq(ibl),rb(ibl)%pres_eq)
        CALL vector_ptassign_laq2(prese_eq(ibl),rb(ibl)%prese_eq)
        CALL vector_ptassign_laq2(tele_eq(ibl),rb(ibl)%tele_eq)
        CALL vector_ptassign_laq2(tion_eq(ibl),rb(ibl)%tion_eq)
        CALL vector_ptassign_laq2(nd_eq(ibl),rb(ibl)%nd_eq)
        CALL vector_ptassign_laq2(diff_shape(ibl),rb(ibl)%diff_shape)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL cvector_ptassign_tl(be(ibl),tb(ibl)%be)
        CALL cvector_ptassign_tl(ja(ibl),tb(ibl)%ja)
        CALL cvector_ptassign_tl(ve(ibl),tb(ibl)%ve)
        CALL cvector_ptassign_tl(pres(ibl),tb(ibl)%pres)
        CALL cvector_ptassign_tl(prese(ibl),tb(ibl)%prese)
        CALL cvector_ptassign_tl(nd(ibl),tb(ibl)%nd)
        CALL cvector_ptassign_tl(conc(ibl),tb(ibl)%conc)
        CALL cvector_ptassign_tl(tele(ibl),tb(ibl)%tele)
        CALL cvector_ptassign_tl(tion(ibl),tb(ibl)%tion)
        CALL cvector_ptassign_tl(work1(ibl),tb(ibl)%work1)
        CALL cvector_ptassign_tl(work2(ibl),tb(ibl)%work2)
        CALL cvector_ptassign_tl(work3(ibl),tb(ibl)%work3)
        CALL cvector_ptassign_tl(work4(ibl),tb(ibl)%work4)
        CALL cvector_ptassign_tl(work5(ibl),tb(ibl)%work5)
        CALL cvector_ptassign_tl(work6(ibl),tb(ibl)%work6)
        CALL cvector_ptassign_tl(w6v1(ibl),tb(ibl)%w6v1)
        CALL cvector_ptassign_tl(w6v2(ibl),tb(ibl)%w6v2)
        CALL vector_ptassign_tl2(be_n0(ibl),tb(ibl)%be_n0)
        CALL vector_ptassign_tl2(pres_n0(ibl),tb(ibl)%pres_n0)
        CALL vector_ptassign_tl2(nd_n0(ibl),tb(ibl)%nd_n0)
        CALL vector_ptassign_tl2(ve_n0(ibl),tb(ibl)%ve_n0)
        CALL vector_ptassign_tl2(rwork1(ibl),tb(ibl)%rwork1)
        CALL vector_ptassign_tl2(be_eq(ibl),tb(ibl)%be_eq)
        CALL vector_ptassign_tl2(ja_eq(ibl),tb(ibl)%ja_eq)
        CALL vector_ptassign_tl2(ve_eq(ibl),tb(ibl)%ve_eq)
        CALL vector_ptassign_tl2(pres_eq(ibl),tb(ibl)%pres_eq)
        CALL vector_ptassign_tl2(prese_eq(ibl),tb(ibl)%prese_eq)
        CALL vector_ptassign_tl2(tele_eq(ibl),tb(ibl)%tele_eq)
        CALL vector_ptassign_tl2(tion_eq(ibl),tb(ibl)%tion_eq)
        CALL vector_ptassign_tl2(nd_eq(ibl),tb(ibl)%nd_eq)
        CALL vector_ptassign_tl2(diff_shape(ibl),tb(ibl)%diff_shape)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pointer_init
c-----------------------------------------------------------------------
c     subprogram 10. q_applied_init.
c     create a spatial function for applying a heat flux to the
c     n=0 mode.  this is usually only used when an equilibrium field is
c     not available.

c-PRE uncomment here and in prhs in integrands.f.
c-----------------------------------------------------------------------
      SUBROUTINE q_applied_init
      USE local
      USE fields
      USE seam_storage_mod
      USE input
      IMPLICIT NONE

      INTEGER(i4) :: ibl,ibe,iv,ix,iy
c-----------------------------------------------------------------------
c     allocate arrays for an applied heat flux
c-----------------------------------------------------------------------
c     ALLOCATE(q_applied(nbl))
c     DO ibl=1,nrbl
c       CALL lagr_quad_alloc(rb(ibl)%q_applied,rb(ibl)%mx,rb(ibl)%my,
c    $                       1_i4,poly_degree,'qa',(/'q_appl'/))
c       rb(ibl)%q_applied=0
c       CALL vector_ptassign_laq2(q_applied(ibl),rb(ibl)%q_applied)
c     ENDDO
c     DO ibl=nrbl+1,nbl
c       CALL tri_linear_alloc(tb(ibl)%q_applied,tb(ibl)%mvert,
c    $                        1_i4,'qa',(/'q_appl'/))
c       CALL vector_ptassign_tl2(q_applied(ibl),tb(ibl)%q_applied)
c     ENDDO
c-----------------------------------------------------------------------
c     loop over blocks touching the edge of the domain for applying
c     heat flux.
c-----------------------------------------------------------------------
c     block: DO ibe=1,SIZE(exblock_list)
c       ibl=exblock_list(ibe)
c-----------------------------------------------------------------------
c       loop over the block boundary.
c-----------------------------------------------------------------------
c       vert: DO iv=1,seam(ibl)%nvert
c         IF (.NOT.seam(ibl)%expoint(iv)) CYCLE
c         ix=seam(ibl)%vertex(iv)%intxy(1)
c         iy=seam(ibl)%vertex(iv)%intxy(2)
c         IF(ix == 0)q_applied(ibl)%arr(1,ix,iy)=-1/per_length
c       ENDDO vert
c     ENDDO block
c     DO ibl=1,nrbl
c       DO ibasis=1,SIZE(rb(ibl)%diff_shape%dx)
c         ix0=rb(ibl)%diff_shape%ix0(ibasis)
c         iy0=rb(ibl)%diff_shape%iy0(ibasis)
c         DO iy=iy0,rb(ibl)%my
c           DO ix=ix0,rb(ibl)%mx
c             CALL lagr_quad_eval(rb(ibl)%rz,
c    &                            ix-ix0+rb(ibl)%be%dx(ibasis),
c    &                            iy-iy0+rb(ibl)%be%dy(ibasis),0_i4)
c             IF((ABS(rb(ibl)%rz%f(2)) < 0.125).AND.
c    &           (ABS(rb(ibl)%rz%f(1)-0.5) < 0.125))THEN
c               tmp=-1
c             ELSE
c               tmp=0
c             ENDIF
c             CALL lagr_quad_basis_assign_loc(rb(ibl)%q_applied,
c    &                                        tmp,ibasis,ix,iy)
c           ENDDO
c         ENDDO
c       ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE q_applied_init
c-----------------------------------------------------------------------
c     subprogram 11. fourier_init.
c     when dealiasing is not used in a nonlinear computation, the last
c     Fourier component must be real.
c-----------------------------------------------------------------------
      SUBROUTINE fourier_init
      USE local
      USE fields
      USE input
      USE global
      IMPLICIT NONE
      
      INTEGER(i4) :: ibl,imode
c-----------------------------------------------------------------------
c     use the nindex array to check for the last Fourier component in
c     parallel computations.
c-----------------------------------------------------------------------
      IF (nindex(nmodes)==nphi/2) THEN
        DO ibl=1,nbl
          CALL cvector_real_comp(ve(ibl),nmodes)
          CALL cvector_real_comp(be(ibl),nmodes)
          CALL cvector_real_comp(nd(ibl),nmodes)
          CALL cvector_real_comp(tele(ibl),nmodes)
          CALL cvector_real_comp(tion(ibl),nmodes)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fourier_init
c-----------------------------------------------------------------------
c     subprogram 12. proj_mat_init.
c     this routine generates a mass matrix that is condensed and
c     factored for use when projecting scalars.
c-----------------------------------------------------------------------
      SUBROUTINE proj_mat_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE matrix_mod
      USE matrix_storage_mod
      USE integrands
      USE regularity
      USE boundary
      USE input
      USE iter_cg
      IMPLICIT NONE

      INTEGER(i4) :: ibl,iv,ivp,ibe,ix,iy,iq,is,mxb,myb,mv
      REAL(r8), DIMENSION(2) :: df
c-----------------------------------------------------------------------
c     nmodes copies are created in order to apply static condensation
c     with standard routines.
c-----------------------------------------------------------------------
      ALLOCATE(sproj_mat(nmodes))
      ALLOCATE(sproj_fac(nmodes))
c-----------------------------------------------------------------------
c     create the sproj matrix.
c-----------------------------------------------------------------------
      DO iq=1,nmodes
        sproj_mat(iq)%fcomp=nindex(iq)
        sproj_mat(iq)%foff=0
        sproj_mat(iq)%eliminated=.false.
        sproj_mat(iq)%symmetric=.true.
        sproj_mat(iq)%nqty=1
        sproj_mat(iq)%nqdis=0
        ALLOCATE(sproj_mat(iq)%vcomp(1))
        sproj_mat(iq)%vcomp(1)='s'
        ALLOCATE(sproj_mat(iq)%rbl_mat(nrbl))
        DO ibl=1,nrbl
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          CALL matrix_rbl_alloc(sproj_mat(iq)%rbl_mat(ibl),mxb,myb,
     $                          1_i4,poly_degree,0_i4,0_i4)
          CALL rblock_make_matrix(rb(ibl),sproj_mat(iq)%rbl_mat(ibl),
     $                            get_mass,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       tblocks:
c-----------------------------------------------------------------------
        ALLOCATE(sproj_mat(iq)%tbl_mat(nrbl+1:nbl))
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_alloc(sproj_mat(iq)%tbl_mat(ibl),
     $                          tb(ibl)%tgeom,1_i4)
          CALL tblock_make_matrix(tb(ibl),sproj_mat(iq)%
     $                            tbl_mat(ibl)%lmat,get_mass,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       apply the dirichlet r=0 conditions if needed.
c-----------------------------------------------------------------------
        CALL regular_op(sproj_mat(iq))
c-----------------------------------------------------------------------
c       apply essential conditions according to what is used for
c       particle density.
c-----------------------------------------------------------------------
        IF (nd_bc=='dirichlet')
     $    CALL dirichlet_op(sproj_mat(iq),'sd',1._r8)
c-----------------------------------------------------------------------
c       eliminate connections to cell-centered coefficients.
c-----------------------------------------------------------------------
        CALL matelim_real_inv_int(sproj_mat(iq),1_i4)
c-----------------------------------------------------------------------
c       collect matrix elements at degenerate points.
c-----------------------------------------------------------------------
        deg_bl: DO ibl=1,nrbl
          IF (rb(ibl)%degenerate) THEN
            IF (solver(1:8)=='bl_drect') THEN
              myb=rb(ibl)%my
              DO ix=1,rb(ibl)%mx
                df=ABS(rb(ibl)%rz%fs(:,ix,0)-rb(ibl)%rz%fs(:,ix,myb))
                IF (df(1)<1.e-12.AND.df(2)<1.e-12) CYCLE deg_bl
              ENDDO
            ENDIF
            CALL matrix_degen_collect_real
     $        (sproj_mat(iq)%rbl_mat(ibl),1_i4,.true.)
          ENDIF
        ENDDO deg_bl
c-----------------------------------------------------------------------
c       factor the mass matrix.
c-----------------------------------------------------------------------
        CALL iter_fac_alloc(sproj_mat(iq),sproj_fac(iq),
     $                      1_i4,solver,.true.)
        CALL iter_factor(sproj_mat(iq),sproj_fac(iq),1_i4,solver,
     $                   off_diag_fac)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE proj_mat_init
c-----------------------------------------------------------------------
c     subprogram 13. hv_mat_init.
c     generate the grad**2 and implicit matrices that get used for
c     time-split hyper-viscous dissipation.
c-----------------------------------------------------------------------
      SUBROUTINE hv_mat_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE matrix_mod
      USE matrix_storage_mod
      USE integrands
      USE regularity
      USE boundary
      USE input
      USE iter_cg
      USE global
      IMPLICIT NONE

      INTEGER(i4) :: ibl,mxb,myb
c-----------------------------------------------------------------------
c     nmodes copies are created in order to apply static condensation
c     with standard routines.
c-----------------------------------------------------------------------
      ALLOCATE(hypv_expmat(nmodes))
      ALLOCATE(hypv_cmat(nmodes))
      ALLOCATE(hypv_cfac(nmodes))
c-----------------------------------------------------------------------
c     allocate the hypv_expmat (grad**2) and implicit matrices, and
c     create the former if it is used.
c-----------------------------------------------------------------------
      DO jmode=1,nmodes
        hypv_expmat(jmode)%fcomp=nindex(jmode)
        hypv_expmat(jmode)%foff=0
        hypv_expmat(jmode)%eliminated=.false.
        hypv_expmat(jmode)%symmetric=.true.
        hypv_expmat(jmode)%nqty=3
        hypv_expmat(jmode)%nqdis=0
        ALLOCATE(hypv_expmat(jmode)%vcomp(3))
        hypv_expmat(jmode)%vcomp=(/'r','z','p'/)
        ALLOCATE(hypv_expmat(jmode)%rbl_mat(nrbl))
        ALLOCATE(hypv_expmat(jmode)%tbl_mat(nrbl+1:nbl))

        hypv_cmat(jmode)%fcomp=nindex(jmode)
        hypv_cmat(jmode)%foff=0
        hypv_cmat(jmode)%eliminated=.false.
        hypv_cmat(jmode)%hermitian=.false.
        hypv_cmat(jmode)%nqty=3
        hypv_cmat(jmode)%nqdis=0
        ALLOCATE(hypv_cmat(jmode)%vcomp(3))
        hypv_cmat(jmode)%vcomp=(/'r','z','p'/)
        ALLOCATE(hypv_cmat(jmode)%rbl_mat(nrbl))
        ALLOCATE(hypv_cmat(jmode)%tbl_mat(nrbl+1:nbl))

        DO ibl=1,nrbl
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          CALL matrix_rbl_alloc(hypv_cmat(jmode)%rbl_mat(ibl),mxb,myb,
     $                          3_i4,poly_degree,0_i4,0_i4)
          CALL matrix_rbl_alloc(hypv_expmat(jmode)%rbl_mat(ibl),mxb,
     $                          myb,3_i4,poly_degree,0_i4,0_i4)
          CALL rblock_make_matrix(rb(ibl),hypv_expmat(jmode)%
     $                            rbl_mat(ibl),vec_lap_op2,3_i4)
        ENDDO
c-----------------------------------------------------------------------
c       tblocks:
c-----------------------------------------------------------------------
        DO ibl=nrbl+1,nbl
          CALL matrix_tbl_alloc(hypv_cmat(jmode)%tbl_mat(ibl),
     $                          tb(ibl)%tgeom,3_i4)
          CALL matrix_tbl_alloc(hypv_expmat(jmode)%tbl_mat(ibl),
     $                          tb(ibl)%tgeom,3_i4)
          CALL tblock_make_matrix(tb(ibl),hypv_expmat(jmode)%
     $                            tbl_mat(ibl)%lmat,vec_lap_op2,3_i4)
        ENDDO
c-----------------------------------------------------------------------
c       allocate the factor structure.
c-----------------------------------------------------------------------
        CALL iter_fac_alloc(hypv_cmat(jmode),hypv_cfac(jmode),
     $                      3_i4,solver,.true.)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE hv_mat_init
