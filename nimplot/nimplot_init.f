c-----------------------------------------------------------------------
c     file nimplot_init.f:  contains the initialization routines for
c     the finite element computations performed in nimplot.  the
c     routines in this file are, for the most part, abbreviated, serial
c     versions of what is found in nimrod_init.f.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  variable_alloc.
c     2.  mass_mat_init.
c     3.  pointer_init.
c     4.  boundary_init.
c-----------------------------------------------------------------------
c     subprogram 1. variable_alloc.
c     allocates space for the dependent variables.
c-----------------------------------------------------------------------
      SUBROUTINE variable_alloc(nmodes)
      USE local
      USE fields
      USE input
      USE computation_pointers
      USE physdat
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes

      INTEGER(i4) :: ibl,mxb,myb,im,mv,mc,iv,iq,jq,nn,ibasis,ix,iy,
     $               ix0,iy0
      COMPLEX(r8), DIMENSION(1) :: ctmp
      REAL(r8), DIMENSION(1) :: te,ti
c-----------------------------------------------------------------------
c     allocate generic blocks for finite element computations.
c-----------------------------------------------------------------------
      ALLOCATE(rhs(nbl))
      ALLOCATE(crhs(nbl))
      ALLOCATE(cell_crhs(nbl))
      ALLOCATE(cell_rhs(nbl))
      ALLOCATE(sln(nbl))
      ALLOCATE(csln(nbl))
      ALLOCATE(vectr(nbl))
      ALLOCATE(cvectr(nbl))
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
              te=rb(ibl)%prese_eq%f(1)/(kboltz*rb(ibl)%nd_eq%f(1))
              ti=(rb(ibl)%pres_eq%f(1)-rb(ibl)%prese_eq%f(1))*
     $           zeff/(kboltz*rb(ibl)%nd_eq%f(1))
              CALL lagr_quad_basis_assign_loc
     $          (rb(ibl)%tele_eq,te,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc
     $          (rb(ibl)%tion_eq,ti,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       allocate storage for temporary computations.
c-----------------------------------------------------------------------
        CALL lagr_quad_alloc(rb(ibl)%work1,mxb,myb,3_i4,nmodes,
     $                       poly_degree,'w1',(/'work1 '/))
        CALL lagr_quad_alloc(rb(ibl)%work2,mxb,myb,5_i4,1_i4,
     $                       poly_degree,'w2',(/'work2 '/))
        CALL lagr_quad_alloc(rb(ibl)%work3,mxb,myb,1_i4,nmodes,
     $                       poly_degree,'w3',(/'work3 '/))
        rb(ibl)%work1=0
c-----------------------------------------------------------------------
c       allocate storage for the finite element vector computations.
c       rhs structures are now allocated and deallocated on the fly to
c       save memory.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(sln(ibl),poly_degree,mxb,myb,1_i4)
        CALL vector_type_alloc(csln(ibl),poly_degree,mxb,myb,1_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,mxb,myb,1_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,mxb,myb,1_i4)
      ENDDO rblocks
c-----------------------------------------------------------------------
c     blocks of unstructured triangles, non-fundamentals first:
c-----------------------------------------------------------------------
      tblocks: DO ibl=nrbl+1,nbl
        mv=tb(ibl)%mvert
        mc=tb(ibl)%mcell
        CALL tri_linear_alloc(tb(ibl)%ja,mv,3_i4,nmodes,name='ja',
     $                        title=(/'  ja  '/))
        tb(ibl)%ja=0
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
c       allocate storage for temporary computations.
c-----------------------------------------------------------------------
        CALL tri_linear_alloc(tb(ibl)%work1,mv,3_i4,nmodes,name='w1',
     $                        title=(/'work1 '/))
        CALL tri_linear_alloc(tb(ibl)%work2,mv,5_i4,1_i4,name='w2',
     $                        title=(/'work2 '/))
        CALL tri_linear_alloc(tb(ibl)%work3,mv,1_i4,nmodes,name='w3',
     $                        title=(/'work3 '/))
        tb(ibl)%work1=0
c-----------------------------------------------------------------------
c       allocate storage for the finite element vector computations.
c       rhs structures are now allocated and deallocated on the fly to
c       save memory.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(sln(ibl),1_i4,mv,0_i4,1_i4)
        CALL vector_type_alloc(csln(ibl),1_i4,mv,0_i4,1_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,mv,0_i4,1_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,mv,0_i4,1_i4)
      ENDDO tblocks
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE variable_alloc
c-----------------------------------------------------------------------
c     subprogram 2. mass_mat_init.
c     generate the mass matrices.
c-----------------------------------------------------------------------
      SUBROUTINE mass_mat_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE nimplot_ints
      USE seam_storage_mod
      USE input
      USE matrix_storage_mod
      USE matrix_mod
      USE iter_cg
      IMPLICIT NONE

      INTEGER(i4) :: ibl,iv,nn,ibe,ix,iy,mxb,myb,imat,jmat
c-----------------------------------------------------------------------
c     create the mass matrix in real and complex array storage.
c     rblocks:
c-----------------------------------------------------------------------
      mass_mat%nqty=1
      mass_mat%eliminated=.false.
      mass_mat%symmetric=.true.
      cmass_mat%nqty=1
      cmass_mat%eliminated=.false.
      cmass_mat%hermitian=.true.
      ALLOCATE(mass_mat%rbl_mat(nrbl))
      ALLOCATE(cmass_mat%rbl_mat(nrbl))
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL matrix_rbl_alloc(mass_mat%rbl_mat(ibl),mxb,myb,
     $                        1_i4,poly_degree,0_i4,0_i4)
        CALL matrix_rbl_alloc(cmass_mat%rbl_mat(ibl),mxb,myb,
     $                        1_i4,poly_degree,0_i4,0_i4)
        CALL rblock_make_matrix(rb(ibl),mass_mat%rbl_mat(ibl),
     $                          get_mass,1_i4)
        IF (lump_all) CALL matrix_lump(mass_mat%rbl_mat(ibl),1_i4)
        IF (rb(ibl)%degenerate)
     $    CALL matrix_degen_collect_real(mass_mat%rbl_mat(ibl),1_i4,
     $                                   mass_mat%symmetric)
        DO imat=1,MIN(poly_degree**2,4_i4)
          DO jmat=1,MIN(poly_degree**2,4_i4)
            cmass_mat%rbl_mat(ibl)%mat(imat,jmat)%arr=
     $        mass_mat%rbl_mat(ibl)%mat(imat,jmat)%arr
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     tblocks:
c-----------------------------------------------------------------------
      ALLOCATE(mass_mat%tbl_mat(nrbl+1:nbl))
      ALLOCATE(cmass_mat%tbl_mat(nrbl+1:nbl))
      DO ibl=nrbl+1,nbl
        CALL matrix_tbl_alloc(mass_mat%tbl_mat(ibl),tb(ibl)%tgeom,1_i4)
        CALL matrix_tbl_alloc(cmass_mat%tbl_mat(ibl),tb(ibl)%tgeom,1_i4)
        CALL tblock_make_matrix(tb(ibl),mass_mat%tbl_mat(ibl)%lmat,
     $                          get_mass,1_i4)
        IF (lump_all) CALL matrix_lump(mass_mat%tbl_mat(ibl)%lmat,1_i4)
        DO iv=0,tb(ibl)%mvert
          cmass_mat%tbl_mat(ibl)%lmat(iv)%element=
     $      mass_mat%tbl_mat(ibl)%lmat(iv)%element
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set-up the preconditioner structure. 
c-----------------------------------------------------------------------
      CALL iter_fac_alloc(mass_mat,mass_fac,1_i4,solver,.false.)
      CALL iter_fac_alloc(cmass_mat,cmass_fac,1_i4,solver,.false.)
      DO ibl=1,nrbl
        mass_fac%bl_fac(ibl)%degenerate=rb(ibl)%degenerate
        cmass_fac%bl_fac(ibl)%degenerate=rb(ibl)%degenerate
      ENDDO
      DO ibl=nrbl+1,nbl
        mass_fac%bl_fac(ibl)%degenerate=.false.
        cmass_fac%bl_fac(ibl)%degenerate=.false.
      ENDDO
      CALL iter_factor(mass_mat,mass_fac,1_i4,solver,off_diag_fac)
      CALL iter_factor(cmass_mat,cmass_fac,1_i4,solver,off_diag_fac)
c-----------------------------------------------------------------------
c     compute averaging factors used in the iterative solve.
c-----------------------------------------------------------------------
      CALL border_weights(solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mass_mat_init
c-----------------------------------------------------------------------
c     subprogram 3. pointer_init.
c     allocate and assigns pointers to nodal data, which are independent
c     of block type.
c-----------------------------------------------------------------------
      SUBROUTINE pointer_init
      USE local
      USE fields

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     allocate block-wise arrays for all spatially dependent fields.
c-----------------------------------------------------------------------
      ALLOCATE(be(nbl),be_n0(nbl),ja(nbl),ve(nbl),pres(nbl),
     $  pres_n0(nbl),prese(nbl),nd(nbl),conc(nbl),
     $  work1(nbl),work2(nbl),work3(nbl),
     $  be_eq(nbl),ja_eq(nbl),ve_eq(nbl),pres_eq(nbl),prese_eq(nbl),
     $  tele_eq(nbl),tion_eq(nbl),
     $  nd_eq(nbl),diff_shape(nbl),si_nl_pres(nbl),
     $  tele(nbl),tion(nbl))
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
        CALL vector_ptassign_laq2(be_n0(ibl),rb(ibl)%be_n0)
        CALL vector_ptassign_laq2(pres_n0(ibl),rb(ibl)%pres_n0)
        CALL vector_ptassign_laq2(si_nl_pres(ibl),rb(ibl)%si_nl_pres)
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
        CALL vector_ptassign_tl2(be_n0(ibl),tb(ibl)%be_n0)
        CALL vector_ptassign_tl2(pres_n0(ibl),tb(ibl)%pres_n0)
        CALL vector_ptassign_tl2(si_nl_pres(ibl),tb(ibl)%si_nl_pres)
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
c     subprogram 4. boundary_init.
c     initializes arrays used for indicating boundary and r=0
c     vertices.
c-----------------------------------------------------------------------
      SUBROUTINE boundary_init(geom)
      USE local
      USE fields
      USE seam_storage_mod
      USE edge
      USE computation_pointers
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: iv,ip,np,ibl,ibv,mx,my,ix,iy,nexbl,nv,ibe,
     $               in,ierror
      LOGICAL, DIMENSION(nbl) :: blmask
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
          blmask(ibl)=.true.
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
c     create an external point mask.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        nv=seam(ibl)%nvert
        ALLOCATE(seam(ibl)%expoint(nv))
        seam(ibl)%expoint=.false.
      ENDDO
      DO iv=1,seam0%nvert
        np=SIZE(seam0%vertex(iv)%ptr,2)
        DO ip=1,np
          ibl=seam0%vertex(iv)%ptr(1,ip)
          ibv=seam0%vertex(iv)%ptr(2,ip)
          seam(ibl)%expoint(ibv)=.true.
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find r=0 points, set the r0point flags, and reset expoint.
c     the r=0 points are assumed to lie on the edge of the domain;
c     hence the blocks are in exblock_list at this point.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        ALLOCATE(seam(ibl)%r0point(seam(ibl)%nvert))
        seam(ibl)%r0point=.false.
      ENDDO
      IF (geom=='tor') THEN
        DO ibe=1,SIZE(exblock_list)
          ibl=exblock_list(ibe)
          nv=seam(ibl)%nvert
          DO iv=1,nv
            ip=iv-1
            in=iv+1
            IF (in>nv) in=1
            IF (ip<1 ) ip=nv
            IF (seam(ibl)%vertex(iv)%rgeom==0) THEN
              seam(ibl)%r0point(iv)=.true.
              rpandn=(/seam(ibl)%vertex(ip)%rgeom,
     $                 seam(ibl)%vertex(in)%rgeom/)
              IF (SUM(rpandn)==0.OR.
     $            (rpandn(1)==0.AND..NOT.seam(ibl)%expoint(in)).OR.
     $            (rpandn(2)==0.AND..NOT.seam(ibl)%expoint(ip))) THEN
                seam(ibl)%expoint(iv)=.false.
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE boundary_init
