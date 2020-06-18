c-----------------------------------------------------------------------
c     file nimeq_init.f:  contains the initialization routines for
c     the finite element computations performed in nimeq.  the
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
      USE input_eq
      USE computation_pointers
      USE physdat
      USE rblock
      USE tblock
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
        CALL rblock_qp_alloc(rb(ibl)%qpres_eq,rb(ibl),1_i4)
        CALL rblock_qp_alloc(rb(ibl)%qja_eq,rb(ibl),3_i4)
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
        CALL lagr_quad_alloc(rb(ibl)%rwork1,mxb,myb,1_i4,
     $                       poly_degree,'rw1',(/'rwork1 '/))
        CALL lagr_quad_alloc(rb(ibl)%rwork2,mxb,myb,1_i4,
     $                       poly_degree,'rw2',(/'rwork2 '/))
        CALL lagr_quad_alloc(rb(ibl)%rwork3,mxb,myb,1_i4,
     $                       poly_degree,'rw3',(/'rwork3 '/))
        CALL lagr_quad_alloc(rb(ibl)%pflux,mxb,myb,1_i4,
     $                       poly_degree,'pflux',(/'pflux '/))
        CALL lagr_quad_alloc(rb(ibl)%fllen,mxb,myb,1_i4,
     $                       poly_degree,'fllen',(/'fllen '/))
        CALL lagr_quad_alloc(rb(ibl)%be_n0,mxb,myb,4_i4,
     $                       poly_degree,'be_n0',(/'be_n0 '/))
        CALL rblock_qp_alloc(rb(ibl)%qrwork1,rb(ibl),1_i4)
        rb(ibl)%work1=0
        rb(ibl)%work2=0
        rb(ibl)%pflux=0
c-----------------------------------------------------------------------
c       allocate storage for the finite element vector computations.
c       rhs structures are now allocated and deallocated on the fly to
c       save memory.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(sln(ibl),poly_degree,mxb,myb,1_i4)
        CALL vector_type_alloc(csln(ibl),poly_degree,mxb,myb,1_i4)
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
        CALL tblock_qp_alloc(tb(ibl)%qpres_eq,tb(ibl),1_i4)
        CALL tblock_qp_alloc(tb(ibl)%qja_eq,tb(ibl),3_i4)
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
        CALL tri_linear_alloc(tb(ibl)%rwork1,mv,1_i4,
     $                        name='rw1',title=(/'rwork1 '/))
        CALL tri_linear_alloc(tb(ibl)%rwork2,mv,1_i4,
     $                        name='rw2',title=(/'rwork2 '/))
        CALL tri_linear_alloc(tb(ibl)%rwork3,mv,1_i4,
     $                        name='rw3',title=(/'rwork3 '/))
        CALL tri_linear_alloc(tb(ibl)%pflux,mv,1_i4,
     $                        name='pflux',title=(/'pflux '/))
        CALL tri_linear_alloc(tb(ibl)%fllen,mv,1_i4,
     $                        name='fllen',title=(/'fllen '/))
        CALL tri_linear_alloc(tb(ibl)%be_n0,mv,4_i4,
     $                        name='be_n0',title=(/'be_n0 '/))
        CALL tblock_qp_alloc(tb(ibl)%qrwork1,tb(ibl),1_i4)
        tb(ibl)%work1=0
        tb(ibl)%work2=0
        tb(ibl)%pflux=0
c-----------------------------------------------------------------------
c       allocate storage for the finite element vector computations.
c       rhs structures are now allocated and deallocated on the fly to
c       save memory.
c-----------------------------------------------------------------------
        CALL vector_type_alloc(sln(ibl),1_i4,mv,0_i4,1_i4)
        CALL vector_type_alloc(csln(ibl),1_i4,mv,0_i4,1_i4)
      ENDDO tblocks
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE variable_alloc
c-----------------------------------------------------------------------
c     subprogram 2. mass_mat_init.
c     generate the mass and other matrices used by nimeq.
c-----------------------------------------------------------------------
      SUBROUTINE mass_mat_init
      USE local
      USE rblock
      USE tblock
      USE fields
      USE nimeq_ints
      USE seam_storage_mod
      USE input
      USE input_eq
      USE matrix_storage_mod
      USE matrix_mod
      USE iter_cg
      USE boundary
      USE regularity
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ibl,iv,nn,ibe,ix,iy,mxb,myb,imat,jmat,iqt,jqt,iq,
     $               iqm,jqm,ierror
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat,mass
      REAL(r8) :: dscale,dtmp
c-----------------------------------------------------------------------
c     compute averaging factors used in the iterative solve.
c-----------------------------------------------------------------------
      CALL border_weights(nimeq_solver)
c-----------------------------------------------------------------------
c     create the mass matrix in real and complex array storage.
c     rblocks:
c-----------------------------------------------------------------------
      mass_mat%nqty=1
      mass_mat%nqdis=0
      mass_mat%fcomp=0
      mass_mat%foff=0
      ALLOCATE(mass_mat%vcomp(1))
      mass_mat%vcomp(1)='s'
      mass_mat%eliminated=.false.
      mass_mat%symmetric=.true.
      mass_mat%diag_scale=1._r8
      mass_mat%essential_cond="none"
      ALLOCATE(mass_mat%rbl_mat(nrbl))

      delstar_mat%nqty=1
      delstar_mat%nqdis=0
      delstar_mat%fcomp=0
      delstar_mat%foff=0
      ALLOCATE(delstar_mat%vcomp(1))
      delstar_mat%vcomp(1)='s'
      delstar_mat%eliminated=.false.
      delstar_mat%symmetric=.true.
      delstar_mat%diag_scale=1._r8
      delstar_mat%essential_cond="none"

      bfield_mat%nqty=3
      bfield_mat%nqdis=0
      bfield_mat%fcomp=0
      bfield_mat%foff=0
      ALLOCATE(bfield_mat%vcomp(3))
      bfield_mat%vcomp=(/'r','z','p'/)
      bfield_mat%eliminated=.false.
      bfield_mat%symmetric=.true.
      bfield_mat%diag_scale=1._r8
      bfield_mat%essential_cond="none"

      IF (gs_type=='free') THEN
        delstlj_mat%nqty=2
        delstlj_mat%nqdis=0
        delstlj_mat%fcomp=0
        delstlj_mat%foff=0
        ALLOCATE(delstlj_mat%vcomp(2))
        delstlj_mat%vcomp=(/'s','s'/)
        delstlj_mat%eliminated=.false.
        delstlj_mat%symmetric=.true.
        delstlj_mat%diag_scale=1._r8
        delstlj_mat%essential_cond="none"
      ENDIF

      IF (btop_check=='passive adv') THEN
        fladv_mat%nqty=1
        fladv_mat%nqdis=0
        fladv_mat%fcomp=0
        fladv_mat%foff=0
        ALLOCATE(fladv_mat%vcomp(1))
        fladv_mat%vcomp(1)='s'
        fladv_mat%eliminated=.false.
        fladv_mat%symmetric=.true.
        fladv_mat%diag_scale=1._r8
        fladv_mat%essential_cond="none"
        ALLOCATE(fladv_mat%rbl_mat(nrbl))
      ENDIF

      ALLOCATE(bfield_mat%rbl_mat(nrbl))
      ALLOCATE(delstar_mat%rbl_mat(nrbl))
      IF (gs_type=='free') ALLOCATE(delstlj_mat%rbl_mat(nrbl))

      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        CALL matrix_rbl_alloc(mass_mat%rbl_mat(ibl),mxb,myb,
     $                        1_i4,poly_degree,0_i4,0_i4)
        CALL matrix_rbl_alloc(delstar_mat%rbl_mat(ibl),mxb,myb,
     $                        1_i4,poly_degree,0_i4,0_i4)
        CALL matrix_rbl_alloc(bfield_mat%rbl_mat(ibl),mxb,myb,
     $                        3_i4,poly_degree,0_i4,0_i4)
        IF (gs_type=='free')
     $    CALL matrix_rbl_alloc(delstlj_mat%rbl_mat(ibl),mxb,myb,
     $                          2_i4,poly_degree,0_i4,0_i4)

        IF (btop_check=='passive adv')
     $    CALL matrix_rbl_alloc(fladv_mat%rbl_mat(ibl),mxb,myb,
     $                          1_i4,poly_degree,0_i4,0_i4)

        CALL rblock_make_matrix(rb(ibl),mass_mat%rbl_mat(ibl),
     $                          get_mass,1_i4)
        CALL rblock_make_matrix(rb(ibl),delstar_mat%rbl_mat(ibl),
     $                          delstar_op,1_i4)
        CALL rblock_make_matrix(rb(ibl),bfield_mat%rbl_mat(ibl),
     $                          mag_op,3_i4)
        IF (gs_type=='free')
     $    CALL rblock_make_matrix(rb(ibl),delstlj_mat%rbl_mat(ibl),
     $                            delstlj_op,2_i4)
      ENDDO
c-----------------------------------------------------------------------
c     tblocks:
c-----------------------------------------------------------------------
      ALLOCATE(mass_mat%tbl_mat(nrbl+1:nbl))
      ALLOCATE(delstar_mat%tbl_mat(nrbl+1:nbl))
      ALLOCATE(bfield_mat%tbl_mat(nrbl+1:nbl))
      IF (gs_type=='free') ALLOCATE(delstlj_mat%tbl_mat(nrbl+1:nbl))
      IF (btop_check=='passive adv')
     $  ALLOCATE(fladv_mat%tbl_mat(nrbl+1:nbl))
      DO ibl=nrbl+1,nbl
        CALL matrix_tbl_alloc(mass_mat%tbl_mat(ibl),tb(ibl)%tgeom,1_i4)

        CALL matrix_tbl_alloc(delstar_mat%tbl_mat(ibl),tb(ibl)%tgeom,
     $                        1_i4)
        CALL matrix_tbl_alloc(bfield_mat%tbl_mat(ibl),tb(ibl)%tgeom,
     $                        3_i4)
        IF (gs_type=='free')
     $    CALL matrix_tbl_alloc(delstlj_mat%tbl_mat(ibl),tb(ibl)%tgeom,
     $                          2_i4)

        IF (btop_check=='passive adv')
     $    CALL matrix_tbl_alloc(fladv_mat%tbl_mat(ibl),tb(ibl)%tgeom,
     $                          1_i4)

        CALL tblock_make_matrix(tb(ibl),mass_mat%tbl_mat(ibl)%lmat,
     $                          get_mass,1_i4)
        CALL tblock_make_matrix(tb(ibl),delstar_mat%tbl_mat(ibl)%lmat,
     $                          delstar_op,1_i4)
        CALL tblock_make_matrix(tb(ibl),bfield_mat%tbl_mat(ibl)%lmat,
     $                          mag_op,3_i4)
        IF (gs_type=='free')
     $    CALL tblock_make_matrix(tb(ibl),delstlj_mat%tbl_mat(ibl)%
     $                            lmat,delstlj_op,2_i4)
      ENDDO
c----------------------------------------------------------------------
c     find the scaling for the diagonal of the del-star operator.
c----------------------------------------------------------------------
      dscale=0._r8
      DO ibl=1,nrbl
        dscale=MAX(dscale,MAXVAL(ABS(delstar_mat%rbl_mat(ibl)%mat(1,1)%
     $                               arr(1,0,0,1,:,:))))
      ENDDO
      DO ibl=nrbl+1,nbl
        DO iv=0,SIZE(delstar_mat%tbl_mat(ibl)%lmat)-1
          dscale=MAX(dscale,ABS(delstar_mat%tbl_mat(ibl)%lmat(iv)%
     $                          element(1,1,0)))
        ENDDO
      ENDDO
      IF (nprocs>1) THEN
        CALL mpi_allreduce(dscale,dtmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierror)
        dscale=dtmp
      ENDIF
c----------------------------------------------------------------------
c     apply regularity conditions to both operators and essential
c     conditions to the delstar operator.
c----------------------------------------------------------------------
      CALL regular_op(delstar_mat)
      CALL regular_op(bfield_mat)
      IF (gs_type=='free') CALL regular_op(delstlj_mat)
c----------------------------------------------------------------------
c     the symm_region input indicates whether to apply a symmetry
c     condition on the top or bottom of the domain.
c----------------------------------------------------------------------
      SELECT CASE(symm_region)
      CASE("top")
        CALL dirichlet_op(delstar_mat,'all',dscale,"b")
      CASE("bottom")
        CALL dirichlet_op(delstar_mat,'all',dscale,"t")
      CASE DEFAULT
        CALL dirichlet_op(delstar_mat,'all',dscale)
      END SELECT
      IF (gs_type=='free') CALL dirichlet_op(delstlj_mat,'all',dscale)
c-----------------------------------------------------------------------
c     eliminate connections to cell centered data.
c-----------------------------------------------------------------------
      CALL matelim_real_inv_int(delstar_mat,1_i4)
      CALL matelim_real_inv_int(bfield_mat,3_i4)
      IF (gs_type=='free') CALL matelim_real_inv_int(delstlj_mat,2_i4)
c-----------------------------------------------------------------------
c     set-up the preconditioner structure.
c-----------------------------------------------------------------------
      CALL iter_fac_alloc(mass_mat,mass_fac,1_i4,nimeq_solver,.false.)
      CALL iter_fac_alloc(delstar_mat,delstar_fac,1_i4,nimeq_solver,
     $                    .true.)
      CALL iter_fac_alloc(bfield_mat,bfield_fac,3_i4,nimeq_solver,
     $                    .true.)
      IF (gs_type=='free')
     $  CALL iter_fac_alloc(delstlj_mat,delstlj_fac,2_i4,nimeq_solver,
     $                     .true.)

      IF (btop_check=='passive adv')
     $  CALL iter_fac_alloc(fladv_mat,fladv_fac,1_i4,nimeq_solver,
     $                      .true.)

      DO ibl=1,nrbl
        IF (rb(ibl)%degenerate) THEN
          CALL matrix_degen_collect_real(mass_mat%rbl_mat(ibl),1_i4,
     $                                   mass_mat%symmetric)
          CALL matrix_degen_collect_real(delstar_mat%rbl_mat(ibl),1_i4,
     $                                   delstar_mat%symmetric)
          CALL matrix_degen_collect_real(bfield_mat%rbl_mat(ibl),3_i4,
     $                                   bfield_mat%symmetric)
          IF (gs_type=='free')
     $      CALL matrix_degen_collect_real(delstlj_mat%rbl_mat(ibl),
     $                                     2_i4,delstlj_mat%symmetric)
        ENDIF
      ENDDO

      CALL iter_factor(mass_mat,mass_fac,1_i4,nimeq_solver,off_diag_fac)
      CALL iter_factor(delstar_mat,delstar_fac,1_i4,nimeq_solver,
     $                 off_diag_fac)
      CALL iter_factor(bfield_mat,bfield_fac,3_i4,nimeq_solver,
     $                 off_diag_fac)
      IF (gs_type=='free')
     $  CALL iter_factor(delstlj_mat,delstlj_fac,2_i4,nimeq_solver,
     $                   off_diag_fac)
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
     $  nd_eq(nbl),diff_shape(nbl),si_nl_pres(nbl),rwork3(nbl),
     $  rwork1(nbl),rwork2(nbl),tele(nbl),tion(nbl),pflux(nbl),
     $  fllen(nbl))
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
        CALL vector_ptassign_laq2(rwork1(ibl),rb(ibl)%rwork1)
        CALL vector_ptassign_laq2(rwork2(ibl),rb(ibl)%rwork2)
        CALL vector_ptassign_laq2(rwork3(ibl),rb(ibl)%rwork3)
        CALL vector_ptassign_laq2(pflux(ibl),rb(ibl)%pflux)
        CALL vector_ptassign_laq2(fllen(ibl),rb(ibl)%fllen)
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
        CALL vector_ptassign_tl2(rwork1(ibl),tb(ibl)%rwork1)
        CALL vector_ptassign_tl2(rwork2(ibl),tb(ibl)%rwork2)
        CALL vector_ptassign_tl2(rwork3(ibl),tb(ibl)%rwork3)
        CALL vector_ptassign_tl2(pflux(ibl),tb(ibl)%pflux)
        CALL vector_ptassign_tl2(fllen(ibl),tb(ibl)%fllen)
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
      USE regularity
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: iv,ip,np,ibl,ibv,mx,my,ix,iy,nexbl,nv,ibe,
     $               in,ierror,nr0bl,iblocal
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ive
      LOGICAL, DIMENSION(nbl) :: blmask, blr0mask
      REAL(r8), DIMENSION(2) :: rpandn

      LOGICAL, EXTERNAL :: per_block
c-----------------------------------------------------------------------
c     create an external block mask and allocate communication arrays.  
c     set blmask to .true. ONLY if I own the block on the exterior 
c     boundary
c -ech blmask is true for any block that has a vertex on edge
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
        ALLOCATE(seam(ibl)%excorner(nv))
        seam(ibl)%excorner=.false.
      ENDDO
      DO iv=1,seam0%nvert
        np=SIZE(seam0%vertex(iv)%ptr,2)
        DO ip=1,np
          ibl=seam0%vertex(iv)%ptr(1,ip)
          IF (block2proc(ibl)==node) THEN
            iblocal=global2local(ibl)
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
              IF (SUM(rpandn)<=2.e-12.OR.
     $          (rpandn(1)<=1.e-12.AND..NOT.seam(ibl)%expoint(in)).OR.
     $          (rpandn(2)<=1.e-12.AND..NOT.seam(ibl)%expoint(ip)))
     $          seam(ibl)%expoint(iv)=.false.
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c----------------------------------------------------------------------
c     set the r0block list, and reset the external block list.
c----------------------------------------------------------------------
      nr0bl=SUM(r0block_list)
      DEALLOCATE(r0block_list)
      ALLOCATE(r0block_list(nr0bl))
      nexbl=SIZE(exblock_list)
      IF (geom=='tor') THEN
         ibv=1
         r0bl_loop: DO ibl=1, nbl
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
             DEALLOCATE (exblock_list)
             ALLOCATE(exblock_list(nexbl))
             nexbl=0
             DO ibl=1,nbl
                IF (blmask(ibl)) THEN
                   nexbl=nexbl +1
                   exblock_list(nexbl)=ibl
                ENDIF
             ENDDO
          ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     save logical variables that indicate whether each rblock is
c     periodic in the y-direction.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        rb(ibl)%periodicy=per_block(ibl,rb(ibl)%mx)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE boundary_init
