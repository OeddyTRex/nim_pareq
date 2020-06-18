c-----------------------------------------------------------------------
c     file iter_utils.f
c
c     the iter_utils module contains utility routines for performing
c     common math functions during the iterative (and direct solves).
c     module procedure statments are used to give a common name to real
c     and complex versions.
c
c     the iter_dir_fac module contains allocation routines for the
c     sparisty pattern used by the direct solvers.  J. King moved
c     these operations from coding that had been duplicated in
c     iter_cg_f90.f and iter_cg_comp.f, and the information is now
c     saved in the spp data structures.  He also created a new version
c     of the sparsity communication, alloc_sparsity_pattern_dist, that
c     does not require large all_reduce operations, saving memory.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_utils
c     1. iter_resid_real
c     2. iter_resid_comp
c     3. iter_inf_norm_real
c     4. iter_inf_norm_comp
c     5. iter_2_norm_real
c     6. iter_2_norm_comp
c     7. iter_dot_real.
c     8. iter_dot_comp.
c     module iter_dir_fac
c     9. iter_dirfac_alloc.
c     10. assign_node_indices.
c     11. alloc_sparsity_pattern.
c     12. alloc_sparsity_pattern_dist.
c     13. pardir_setup.
c     14. iter_dirfac_dealloc.
c     15. set_precon_opts.
c-----------------------------------------------------------------------
c     module containing mathematical utility routines for the Krylov-
c     space solvers.
c-----------------------------------------------------------------------
      MODULE iter_utils
      USE local
      USE vector_type_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     module-procedure statements to make common routine names.
c-----------------------------------------------------------------------
      INTERFACE iter_resid
        MODULE PROCEDURE iter_resid_real,iter_resid_comp
      END INTERFACE

      INTERFACE iter_inf_norm
        MODULE PROCEDURE iter_inf_norm_real,iter_inf_norm_comp
      END INTERFACE

      INTERFACE iter_2_norm
        MODULE PROCEDURE iter_2_norm_real,iter_2_norm_comp
      END INTERFACE

      INTERFACE iter_dot
        MODULE PROCEDURE iter_dot_real,iter_dot_comp
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iter_resid_real
c     compute the residual (rhs - A.x) for a real system.
c-----------------------------------------------------------------------
      SUBROUTINE iter_resid_real(mat_str,xx,re,rh,nqty,ntotb,nrb,
     $                           poly_deg)
      USE matrix_mod

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty,ntotb,nrb,poly_deg
      TYPE(vector_type), DIMENSION(:), INTENT(IN) :: xx,rh
      TYPE(vector_type), DIMENSION(:), INTENT(OUT) :: re

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     call the matrix-vector multiplication routine then subtract from
c     rhs.
c-----------------------------------------------------------------------
      CALL matvec(mat_str,xx,re,nqty)
      DO ibl=1,ntotb
        re(ibl)%arr=rh(ibl)%arr-re(ibl)%arr
        IF (poly_deg>1.AND.ibl<=nrb) THEN
          re(ibl)%arrh=rh(ibl)%arrh-re(ibl)%arrh
          re(ibl)%arrv=rh(ibl)%arrv-re(ibl)%arrv
          IF (.NOT.mat_str%eliminated)
     $      re(ibl)%arri=rh(ibl)%arri-re(ibl)%arri
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_resid_real
c-----------------------------------------------------------------------
c     subprogram 2. iter_resid_comp
c     compute the residual (rhs - A.x) for a complex system.
c-----------------------------------------------------------------------
      SUBROUTINE iter_resid_comp(mat_str,xx,re,rh,nqty,ntotb,nrb,
     $                           poly_deg)
      USE matrix_mod

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty,ntotb,nrb,poly_deg
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: xx,rh
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(OUT) :: re

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     call the matrix-vector multiplication routine then subtract from
c     rhs.
c-----------------------------------------------------------------------
      CALL matvec(mat_str,xx,re,nqty)
      DO ibl=1,ntotb
        re(ibl)%arr=rh(ibl)%arr-re(ibl)%arr
        IF (poly_deg>1.AND.ibl<=nrb) THEN
          re(ibl)%arrh=rh(ibl)%arrh-re(ibl)%arrh
          re(ibl)%arrv=rh(ibl)%arrv-re(ibl)%arrv
          IF (.NOT.mat_str%eliminated)
     $      re(ibl)%arri=rh(ibl)%arri-re(ibl)%arri
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_resid_comp
c-----------------------------------------------------------------------
c     subprogram 3. iter_inf_norm_real
c     compute the infinity norm of a vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_inf_norm_real(vec,err,er0,ntotb,nrb,poly_deg,
     $                              use_int)
      USE mpi_nim
      USE pardata

      TYPE(vector_type), DIMENSION(:), INTENT(IN) :: vec
      REAL(r8), INTENT(OUT) :: err
      REAL(r8), INTENT(IN) :: er0
      INTEGER(i4), INTENT(IN) :: ntotb,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      REAL(r8) :: tmp
      INTEGER(i4) :: ibl,ierror
c-----------------------------------------------------------------------
c     check each basis type in the calculation.
c-----------------------------------------------------------------------
      err=0
      DO ibl=1,ntotb
        err=MAX(err,MAXVAL(ABS(vec(ibl)%arr)))
        IF (poly_deg>1.AND.ibl<=nrb) THEN
          err=MAX(err,MAXVAL(ABS(vec(ibl)%arrh)))
          err=MAX(err,MAXVAL(ABS(vec(ibl)%arrv)))
          IF (use_int) THEN
            err=MAX(err,MAXVAL(ABS(vec(ibl)%arri)))
          ENDIF
        ENDIF
      ENDDO
      IF (nprocs_layer > 1) THEN
        CALL mpi_allreduce(err,tmp,1,mpi_nim_real,mpi_max,
     $       comm_layer,ierror)
        err = tmp
      ENDIF
      err=err/er0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_inf_norm_real
c-----------------------------------------------------------------------
c     subprogram 4. iter_inf_norm_comp
c     compute the infinity norm of a vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_inf_norm_comp(vec,err,er0,ntotb,nrb,poly_deg,
     $                              use_int)
      USE mpi_nim
      USE pardata

      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: vec
      REAL(r8), INTENT(OUT) :: err
      REAL(r8), INTENT(IN) :: er0
      INTEGER(i4), INTENT(IN) :: ntotb,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      REAL(r8) :: tmp
      INTEGER(i4) :: ibl,ierror
c-----------------------------------------------------------------------
c     check each basis type in the calculation.
c-----------------------------------------------------------------------
      err=0
      DO ibl=1,ntotb
        err=MAX(err,MAXVAL(ABS(REAL(vec(ibl)%arr,r8))))
        err=MAX(err,MAXVAL(ABS(AIMAG(vec(ibl)%arr))))
        IF (poly_deg>1.AND.ibl<=nrb) THEN
          err=MAX(err,MAXVAL(ABS(REAL(vec(ibl)%arrh,r8))))
          err=MAX(err,MAXVAL(ABS(AIMAG(vec(ibl)%arrh))))
          err=MAX(err,MAXVAL(ABS(REAL(vec(ibl)%arrv,r8))))
          err=MAX(err,MAXVAL(ABS(AIMAG(vec(ibl)%arrv))))
          IF (use_int) THEN
            err=MAX(err,MAXVAL(ABS(REAL(vec(ibl)%arri,r8))))
            err=MAX(err,MAXVAL(ABS(AIMAG(vec(ibl)%arri))))
          ENDIF
        ENDIF
      ENDDO
      IF (nprocs_layer > 1) THEN
        CALL mpi_allreduce(err,tmp,1,mpi_nim_real,mpi_max,
     $       comm_layer,ierror)
        err = tmp
      ENDIF
      err=err/er0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_inf_norm_comp
c-----------------------------------------------------------------------
c     subprogram 5. iter_2_norm_real
c     compute the 2-norm of a real vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_2_norm_real(vec,err,er0,ntotb,nrb,poly_deg,
     $                            use_int)
      USE mpi_nim
      USE pardata
      USE seam_storage_mod

      TYPE(vector_type), DIMENSION(:), INTENT(IN) :: vec
      REAL(r8), INTENT(OUT) :: err
      REAL(r8), INTENT(IN) :: er0
      INTEGER(i4), INTENT(IN) :: ntotb,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      REAL(r8) :: tmp
      INTEGER(i4) :: ibl,ierror
c-----------------------------------------------------------------------
c     check each basis type in the calculation.
c-----------------------------------------------------------------------
      err=0
      DO ibl=1,ntotb
        CALL iter_dot(err,vec(ibl),vec(ibl),seam(ibl),ibl,nrb,poly_deg,
     $                use_int)
      ENDDO
      IF (nprocs_layer > 1) THEN
        CALL mpi_allreduce(err,tmp,1,mpi_nim_real,mpi_sum,
     $                     comm_layer,ierror)
        err = tmp
      ENDIF
      err=SQRT(err)/er0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_2_norm_real
c-----------------------------------------------------------------------
c     subprogram 6. iter_2_norm_comp
c     compute the 2-norm of a complex vector.
c-----------------------------------------------------------------------
      SUBROUTINE iter_2_norm_comp(vec,err,er0,ntotb,nrb,poly_deg,
     $                            use_int)
      USE mpi_nim
      USE pardata
      USE seam_storage_mod

      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: vec
      REAL(r8), INTENT(OUT) :: err
      REAL(r8), INTENT(IN) :: er0
      INTEGER(i4), INTENT(IN) :: ntotb,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      REAL(r8) :: tmp
      COMPLEX(r8) :: prod
      INTEGER(i4) :: ibl,ierror
c-----------------------------------------------------------------------
c     check each basis type in the calculation.
c-----------------------------------------------------------------------
      prod=0
      DO ibl=1,ntotb
        CALL iter_dot(prod,vec(ibl),vec(ibl),seam(ibl),ibl,nrb,
     $                poly_deg,use_int)
      ENDDO
      err=prod
      IF (nprocs_layer > 1) THEN
        CALL mpi_allreduce(err,tmp,1,mpi_nim_real,mpi_sum,
     $                     comm_layer,ierror)
        err = tmp
      ENDIF
      err=SQRT(err)/er0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_2_norm_comp
c-----------------------------------------------------------------------
c     subprogram 7. iter_dot_real.
c     compute the inner product of two vectors over a block.
c-----------------------------------------------------------------------
      SUBROUTINE iter_dot_real(dprod,v1,v2,bseam,ibl,nrb,poly_deg,
     $                         use_int)
      USE edge_type_mod

      REAL(r8), INTENT(INOUT) :: dprod
      TYPE(vector_type), INTENT(IN) :: v1,v2
      TYPE(edge_type), INTENT(IN) :: bseam
      INTEGER(i4), INTENT(IN) :: ibl,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      INTEGER(i4) :: ix,iy,iv,ni
c-----------------------------------------------------------------------
c     find the dot product at each grid vertex.
c-----------------------------------------------------------------------
      dprod=dprod+SUM(v1%arr*v2%arr)
c-----------------------------------------------------------------------
c     boundary points must be divided by the number of internal
c     representations to get the correct sum over all blocks.
c     remove what has been contributed above, also.
c-----------------------------------------------------------------------
      DO iv=1,bseam%nvert
        ix=bseam%vertex(iv)%intxy(1)
        iy=bseam%vertex(iv)%intxy(2)
        dprod=dprod+(bseam%vertex(iv)%ave_factor-1._r8)*
     $               SUM(v1%arr(:,ix,iy)*v2%arr(:,ix,iy))
      ENDDO
c-----------------------------------------------------------------------
c-PRE higher order contributions.
c-----------------------------------------------------------------------
      IF (poly_deg>1.AND.ibl<=nrb) THEN
        dprod=dprod+SUM(v1%arrh*v2%arrh)+SUM(v1%arrv*v2%arrv)
        IF (use_int) THEN
          dprod=dprod+SUM(v1%arri*v2%arri)
        ENDIF
c-----------------------------------------------------------------------
c-PRE   boundary segment-centered data must be divided by 2;
c       however, the average factor precludes contributions from
c       tblocks.
c-----------------------------------------------------------------------
        DO iv=1,bseam%nvert
          ix=bseam%segment(iv)%intxys(1)
          iy=bseam%segment(iv)%intxys(2)
          IF (bseam%segment(iv)%h_side) THEN
            dprod=dprod+(bseam%segment(iv)%ave_factor-1._r8)*
     $                   SUM(v1%arrh(:,:,ix,iy)*v2%arrh(:,:,ix,iy))
          ELSE
            dprod=dprod+(bseam%segment(iv)%ave_factor-1._r8)*
     $                   SUM(v1%arrv(:,:,ix,iy)*v2%arrv(:,:,ix,iy))
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_dot_real
c-----------------------------------------------------------------------
c     subprogram 8. iter_dot_comp.
c     compute the inner product of two complex vectors over a block.
c-----------------------------------------------------------------------
      SUBROUTINE iter_dot_comp(dprod,v1,v2,bseam,ibl,nrb,poly_deg,
     $                         use_int)
      USE edge_type_mod

      COMPLEX(r8), INTENT(INOUT) :: dprod
      TYPE(cvector_2D_type), INTENT(IN) :: v1,v2
      TYPE(edge_type), INTENT(IN) :: bseam
      INTEGER(i4), INTENT(IN) :: ibl,nrb,poly_deg
      LOGICAL, INTENT(IN) :: use_int

      INTEGER(i4) :: ix,iy,iv,ni
c-----------------------------------------------------------------------
c     find the dot product at each grid vertex.
c-----------------------------------------------------------------------
      dprod=dprod+SUM(CONJG(v1%arr)*v2%arr)
c-----------------------------------------------------------------------
c     boundary points must be divided by the number of internal
c     representations to get the correct sum over all blocks.
c     remove what has been contributed above, also.
c-----------------------------------------------------------------------
      DO iv=1,bseam%nvert
        ix=bseam%vertex(iv)%intxy(1)
        iy=bseam%vertex(iv)%intxy(2)
        dprod=dprod+(bseam%vertex(iv)%ave_factor-1._r8)*
     $               SUM(CONJG(v1%arr(:,ix,iy))*v2%arr(:,ix,iy))
      ENDDO
c-----------------------------------------------------------------------
c-PRE higher order contributions.
c-----------------------------------------------------------------------
      IF (poly_deg>1.AND.ibl<=nrb) THEN
        dprod=dprod+SUM(CONJG(v1%arrh)*v2%arrh)+
     $              SUM(CONJG(v1%arrv)*v2%arrv)
        IF (use_int) THEN
          dprod=dprod+SUM(CONJG(v1%arri)*v2%arri)
        ENDIF
c-----------------------------------------------------------------------
c-PRE   boundary segment-centered data must be divided by 2;
c       however, the average factor precludes contributions from
c       tblocks.
c-----------------------------------------------------------------------
        DO iv=1,bseam%nvert
          ix=bseam%segment(iv)%intxys(1)
          iy=bseam%segment(iv)%intxys(2)
          IF (bseam%segment(iv)%h_side) THEN
            dprod=dprod+(bseam%segment(iv)%ave_factor-1._r8)*
     $                 SUM(CONJG(v1%arrh(:,:,ix,iy))*v2%arrh(:,:,ix,iy))
          ELSE
            dprod=dprod+(bseam%segment(iv)%ave_factor-1._r8)*
     $                 SUM(CONJG(v1%arrv(:,:,ix,iy))*v2%arrv(:,:,ix,iy))
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_dot_comp
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_utils

c-----------------------------------------------------------------------
c     module containing the allocation routines for the sparsity
c     pattern used by the direct solvers.
c-PRE
c     triangles could be included!
c-----------------------------------------------------------------------
      MODULE iter_dir_fac
      USE local
      IMPLICIT NONE
      PUBLIC iter_dirfac_alloc
      PUBLIC iter_dirfac_dealloc
      PRIVATE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 9. iter_dirfac_alloc.
c
c     control routine for the allocation of the sparsity pattern factor
c     structure for use by the direct solvers.
c-----------------------------------------------------------------------
      SUBROUTINE iter_dirfac_alloc(mat_info,spp,bl_spp,int_elim,
     $                             nrb,nqty,precon)
      USE factor_type_mod

      TYPE(rbl_mat_info),     INTENT(IN)    :: mat_info(:)
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      TYPE(bl_sparsity_type), INTENT(INOUT) :: bl_spp(:)
      LOGICAL, INTENT(IN) :: int_elim
      INTEGER(i4), INTENT(IN) :: nrb,nqty
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: nbmax
      LOGICAL, EXTERNAL :: direct_check

c-----------------------------------------------------------------------
c     nullify pointers that are tested in the dealloc routine.
c-----------------------------------------------------------------------
      NULLIFY(spp%j_acc_save)
      NULLIFY(spp%start_acc_save)
      NULLIFY(spp%sendrowst)
      NULLIFY(spp%sendstacc)
      spp%solve_nopts=0
c-----------------------------------------------------------------------
c     determine if the matrix and sparsity pattern are distributed.
c-----------------------------------------------------------------------
      spp%matrix_distributed=.FALSE.
      IF (precon(1:7)=="slu_dst")
     $  spp%matrix_distributed=.TRUE.

      spp%sparsity_distributed=.FALSE.
      IF (precon=="slu_dsta")
     $  spp%sparsity_distributed=.TRUE.

      spp%direct_solver=.FALSE.
      IF (direct_check(precon))  spp%direct_solver=.TRUE.

      spp%acc_lustored=.FALSE.

      IF (.NOT.spp%direct_solver) RETURN
c-----------------------------------------------------------------------
c     determine the number of basis functions.
c-----------------------------------------------------------------------
      IF (int_elim) THEN
        nbmax=MIN(mat_info(1)%nbtype,3_i4)
      ELSE
        nbmax=MIN(mat_info(1)%nbtype,4_i4)
      ENDIF
c-----------------------------------------------------------------------
c     assign each node a global index.
c-----------------------------------------------------------------------
      CALL assign_node_indices(mat_info,spp,bl_spp,nrb,nqty,nbmax)
c-----------------------------------------------------------------------
c     set the preconditioning options for external solvers.
c-----------------------------------------------------------------------
      CALL set_precon_opts(precon,spp%solve_nopts,spp%iopts,spp%dopts)
c-----------------------------------------------------------------------
c     allocate and fill the sparsity pattern, and setup the arrays
c     needed for the distributed communication.
c-----------------------------------------------------------------------
      IF (spp%sparsity_distributed) THEN
        CALL alloc_sparsity_pattern_dist(mat_info,spp,bl_spp,nrb,nqty,
     $                                   nbmax)
      ELSE
        CALL alloc_sparsity_pattern(mat_info,spp,bl_spp,nrb,nqty,nbmax)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_dirfac_alloc
c-----------------------------------------------------------------------
c     subprogram 10. assign_node_indices.
c     creates a generic row index for nimrod rblock data structures.
c
c     the following is done for the compressed row storage:
c
c       1. define a unique index (0:spp%nrow-1) for each quantity
c          component at each unique node of the finite element mesh.
c          nodes duplicated by the block structure are claimed by the
c          block with the lowest global ordering that is determined
c          through nproc_layer.
c       2. accumulate a running total of unique nodes, block-by-block:
c            spp%irowst_block(0:nbl_total)
c       3. create a nimrod-matrix-format data structure (i.e. in terms
c          of blocks, basis types, iqty, ix, and iy) that holds the
c          corresponding index from 1 above.  this structure is
c            fac%bl_fac(ibl)%row_ind(ibasis)%rarr(iq,ix,iy)
c
c-----------------------------------------------------------------------
      SUBROUTINE assign_node_indices(mat_info,spp,bl_spp,nrb,nqty,nbmax)
      USE pardata
      USE mpi_nim
      USE factor_type_mod
      USE seam_storage_mod, ONLY: seam
      USE edge, ONLY: edge_network,edge_load_limits

      TYPE(rbl_mat_info),     INTENT(IN)    :: mat_info(:)
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      TYPE(bl_sparsity_type), INTENT(INOUT) :: bl_spp(:)
      INTEGER(i4), INTENT(IN) :: nrb,nqty,nbmax

      INTEGER(i4) :: ix,iy,ix0,iy0,nq,ib,nrank,yst,ibl,mx,my,iv,ip,
     $               nb_total,ierr,nside,id,iq,ibst,iben,iqoff
      INTEGER(i4), ALLOCATABLE :: tmp1(:)
      INTEGER(i4), DIMENSION(1) :: iloc
      CHARACTER(128) :: msg
      LOGICAL, DIMENSION(SIZE(global2local)) :: bmask

c-----------------------------------------------------------------------
c     first allocate the row index table.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        ALLOCATE(bl_spp(ibl)%row_ind(nbmax))
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        DO ib=1,nbmax
          ix0=mat_info(ibl)%ix0(ib)
          iy0=mat_info(ibl)%iy0(ib)
          nq=mat_info(ibl)%nq_type(ib)
          ALLOCATE(bl_spp(ibl)%row_ind(ib)%rarr(nq,ix0:mx,iy0:my))
          bl_spp(ibl)%row_ind(ib)%rarr=1
        ENDDO
      ENDDO
      nb_total=SIZE(global2local)
      ALLOCATE(spp%irowst_block(0:nb_total))
      spp%irowst_block=0
c-----------------------------------------------------------------------
c     generate global block ordering that keeps all rows (hence blocks)
c     on a processor together.
c-----------------------------------------------------------------------
      ALLOCATE(spp%gblock_order(nb_total))
      bmask=.true.
      DO id=1,nb_total
        iloc=MINLOC(block2proc,MASK=bmask)
        spp%gblock_order(iloc(1))=id
        bmask(iloc(1))=.false.
      ENDDO
c-----------------------------------------------------------------------
c     determine which nodes are considered to be owned by each block,
c     based on the global block number.
c
c     spp%irowst_block(0:nb_total) will hold a running total of unique
c     vector-rows from each block.  it is used to test whether a block
c     is the owner of a given vector-row, since the first row index of
c     global-block-# ib is spp%irowst_block(spp%gblock_order(ib)-1).
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          IF (bl_spp(ibl)%perblock.AND.
     $        (iv<=mx.OR.iv==seam(ibl)%nvert).OR.
     $        bl_spp(ibl)%degenerate.AND.iv>2*mx+my) THEN
            bl_spp(ibl)%row_ind(1)%rarr(:,ix,iy)=0
            CYCLE
          ENDIF
          DO ip=1,seam(ibl)%vertex(iv)%nimage
            IF (spp%gblock_order(seam(ibl)%vertex(iv)%ptr2(1,ip))<id)
     $        bl_spp(ibl)%row_ind(1)%rarr(:,ix,iy)=0
          ENDDO
          IF (bl_spp(ibl)%degenerate.AND.iv==2*mx+my) THEN
            DO ip=1,seam(ibl)%vertex(seam(ibl)%nvert)%nimage
              IF (spp%gblock_order(
     $            seam(ibl)%vertex(seam(ibl)%nvert)%ptr2(1,ip))<id)
     $          bl_spp(ibl)%row_ind(1)%rarr(:,ix,iy)=0
            ENDDO
          ENDIF
        ENDDO
        IF (nbmax>1) THEN
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            IF (bl_spp(ibl)%perblock.AND.iv<=mx) THEN
              bl_spp(ibl)%row_ind(2)%rarr(:,ix,iy)=0
              CYCLE
            ENDIF
            IF (bl_spp(ibl)%degenerate.AND.iv>2*mx+my) THEN
              bl_spp(ibl)%row_ind(3)%rarr(:,ix,iy)=0
              CYCLE
            ENDIF
            IF (seam(ibl)%segment(iv)%ptr(1)==0) THEN !  keep boundary
            ELSE IF (spp%gblock_order(seam(ibl)%segment(iv)%ptr(1))<id)
     $        THEN
              IF (seam(ibl)%segment(iv)%h_side) THEN
                bl_spp(ibl)%row_ind(2)%rarr(:,ix,iy)=0
              ELSE
                bl_spp(ibl)%row_ind(3)%rarr(:,ix,iy)=0
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        DO ib=1,nbmax
          spp%irowst_block(id)=spp%irowst_block(id)+
     $      SUM(bl_spp(ibl)%row_ind(ib)%rarr)
        ENDDO
      ENDDO
      IF (nprocs_layer>1) THEN
        ALLOCATE(tmp1(nb_total))
        CALL mpi_allreduce(spp%irowst_block(1),tmp1(1),nb_total,
     $                     mpi_nim_int,mpi_sum,comm_layer,ierr)
        spp%irowst_block(1:)=tmp1
        DEALLOCATE(tmp1)
      ENDIF
      DO id=1,nb_total
        spp%irowst_block(id)=spp%irowst_block(id-1)+
     $                       spp%irowst_block(id)
      ENDDO
      spp%nrow=spp%irowst_block(nb_total)
c-----------------------------------------------------------------------
c     loop over all finite element nodes and assign indices.
c
c     the desired basis type order is numerically decreasing so that
c     the interior, vertical side, horizontal side, and grid vertex
c     bases are picked up in that order for each element.
c
c     all representations of the degenerate point are mapped to
c     a unique set of vector component indices that are ordered
c     after the ix=1 nodes for less fill-in.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        nrank=spp%irowst_block(id-1)
        IF (bl_spp(ibl)%perblock) THEN
          yst=1
        ELSE
          yst=0
        ENDIF
        DO ix=0,mx
          IF (ix==0.AND.bl_spp(ibl)%degenerate) CYCLE
          DO iy=yst,my
            DO ib=nbmax,1,-1
              ix0=mat_info(ibl)%ix0(ib)
              iy0=mat_info(ibl)%iy0(ib)
              nq =mat_info(ibl)%nq_type(ib)
              IF (ix0>ix.OR.iy0>iy) CYCLE
              DO iq=1,nq
                IF (nrank==0.OR.
     $              bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)>0) THEN
                  bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)=nrank
                  nrank=nrank+1
                ENDIF
              ENDDO
            ENDDO
          ENDDO
          IF (ix==1.AND.bl_spp(ibl)%degenerate) THEN
            DO iy=yst,my
              DO iq=1,nqty
                IF (bl_spp(ibl)%row_ind(1)%rarr(iq,0,iy)>0) THEN
                  bl_spp(ibl)%row_ind(1)%rarr(iq,0,iy)=nrank
                  nrank=nrank+1
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF (nrank/=spp%irowst_block(id)) THEN
          WRITE(msg,'(a,i7,a,i7,a,i3)') "Assign_node_indices: "//
     $      "block-rank ",nrank," doesn't match irowst ",
     $      spp%irowst_block(id)," in block ordered ",id
          CALL nim_stop(msg)
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     communicate global row indices.
c-----------------------------------------------------------------------
      IF (nbmax>1) THEN
        nside=mat_info(1)%nq_type(2)/nqty
      ELSE
        nside=0
      ENDIF
      DO ibl=1,nrb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          seam(ibl)%vertex(iv)%seam_in(1:nqty)=
     $      bl_spp(ibl)%row_ind(1)%rarr(:,ix,iy)
          IF (nbmax>1) THEN
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                            nside,ibst,iben)
            id=2
            IF (.NOT.seam(ibl)%segment(iv)%h_side) id=3
            nq=1
            DO ib=ibst,iben,seam(ibl)%segment(iv)%load_dir
              iqoff=nqty*(ib-1)+1
              seam(ibl)%segment(iv)%seam_in(nq:nq+nqty-1)=
     $          bl_spp(ibl)%row_ind(id)%rarr(iqoff:iqoff+nqty-1,ix,iy)
              nq=nq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(nqty,0_i4,nside,.false.)
      DO ibl=1,nrb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          bl_spp(ibl)%row_ind(1)%rarr(:,ix,iy)=
     $      NINT(seam(ibl)%vertex(iv)%seam_out(1:nqty))
          IF (nbmax>1) THEN
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                            nside,ibst,iben)
            id=2
            IF (.NOT.seam(ibl)%segment(iv)%h_side) id=3
            nq=1
            DO ib=ibst,iben,seam(ibl)%segment(iv)%load_dir
              iqoff=nqty*(ib-1)+1
              bl_spp(ibl)%row_ind(id)%rarr(iqoff:iqoff+nqty-1,ix,iy)=
     $          NINT(seam(ibl)%segment(iv)%seam_out(nq:nq+nqty-1))
              nq=nq+nqty
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE assign_node_indices
c-----------------------------------------------------------------------
c     subprogram 11. alloc_sparsity_pattern.
c     allocation routine for the sparsity pattern. this
c     routine picks up after step 3 of assign_node_indices and
c     continues to initialize and allocate the sparsity_pattern spp
c     structure.
c
c     the following additional steps are done for the compressed row
c     storage:
c
c       4. determine the number of nonzero entries in each row of the
c          matrix.  a sum that is decomposed by block is placed in
c            ientst_block(1:nbl_total,0:spp%nrow-1)
c          and the running total is placed in
c            ientst_block(0,0:spp%nrow-1)
c       5. the running total from 4 is equivalent to the index where
c          each row starts in the compressed-row storage for the matrix.
c          it is copied to spp%start_acc(0:spp%nrow-1).
c       6. the column indices for each nonzero matrix entry are
c          determined row-by-row.  they are placed in the array
c            spp%j_acc(0:nnz-1)
c       7. rows owned by this processor are defined.  spp%mloc is the
c          number of these rows and spp%start_loc(0:spp%mloc-1) is the
c          starting location in compressed-row storage for these rows
c          only.  [if the full-storage interface is used,
c          spp%mloc=nrow, and start_loc is not allocated.]
c       8. lists of compressed row indices and rows for all matrix
c          elements with contributions from blocks on this processor
c          are also created and saved in
c            spp%jentry(:)
c            irw(:)
c          to assist communication for distributed-memory solves.
c
c-----------------------------------------------------------------------
      SUBROUTINE alloc_sparsity_pattern(mat_info,spp,bl_spp,
     $                                  nrb,nqty,nbmax)
      USE pardata
      USE mpi_nim
      USE factor_type_mod

      TYPE(rbl_mat_info),     INTENT(IN)    :: mat_info(:)
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      TYPE(bl_sparsity_type), INTENT(INOUT) :: bl_spp(:)
      INTEGER(i4), INTENT(IN) :: nrb,nqty,nbmax

      INTEGER(i4) :: iq,jq,lx,ly,ix,iy,ix0,iy0,jx0,jy0,nq,ib,jb,
     $               nnz,inq,jnq,jx,jy,iind,jind,ientry,nbw,off,
     $               yst,ibl,mx,my,iv,ip,nb_total,ierr,nside,id,
     $               jdeg,nnz2
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ientst_block,tmp2
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: tmp1,nentry,ndispl,
     $             sbd_jacc,sbd_ind,glb_jacc,glb_ind,irw
      CHARACTER(160) :: msg
      LOGICAL, EXTERNAL :: concave_corner

c-----------------------------------------------------------------------
c     determine the total number of nonzero entries that are possible.
c
c     ientst_block will have the starting entry index of the
c     compressed row format for the nonzero entries in each unique
c     nimrod-matrix row, organized by the global block ordering that
c     is defined by spp%gblock_order.
c
c     upon completion, global-block-# ib has unique entries for row
c     irow (numbered from 0:spp%nrow-1) starting at cc index
c     ientst_block(spp%gblock_order(ib)-1,irow).
c-----------------------------------------------------------------------
      nb_total=SIZE(global2local)
      nbw=0
      nnz=0
      nnz2=0
      ALLOCATE(ientst_block(0:nb_total,0:spp%nrow-1))
      ientst_block=0
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        IF (bl_spp(ibl)%perblock) THEN
          yst=1
        ELSE
          yst=0
        ENDIF
        DO ib=1,nbmax
          ix0=mat_info(ibl)%ix0(ib)
          iy0=mat_info(ibl)%iy0(ib)
          inq=mat_info(ibl)%nq_type(ib)
          DO jb=1,nbmax
            jx0=mat_info(ibl)%ix0(jb)
            jy0=mat_info(ibl)%iy0(jb)
            jnq=mat_info(ibl)%nq_type(jb)
            DO ix=ix0,mx
              DO iy=MAX(yst,iy0),my
                DO lx=jx0-1,1-ix0
                  jx=ix+lx
                  IF (jx<jx0.OR.jx>mx) CYCLE
                  DO ly=jy0-1,1-iy0
                    jy=iy+ly
                    IF (bl_spp(ibl)%perblock.AND.jy==my+1) jy=1
                    IF (jy<jy0.OR.jy>my) CYCLE
                    nnz2=nnz2+inq*jnq
                    IF (bl_spp(ibl)%row_ind(ib)%rarr(1,ix,iy)>=
     $                  spp%irowst_block(id-1).OR.
     $                  bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)>=
     $                  spp%irowst_block(id-1).OR.
     $                  concave_corner(lx,ly,ix0,iy0,jx0,jy0)) THEN
                      off=bl_spp(ibl)%row_ind(ib)%rarr(inq,ix,iy)-
     $                    bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)
                      nbw=MAX(nbw,ABS(off))
                      off=bl_spp(ibl)%row_ind(ib)%rarr(1,ix,iy)-
     $                    bl_spp(ibl)%row_ind(jb)%rarr(jnq,jx,jy)
                      nbw=MAX(nbw,ABS(off))
                      IF (bl_spp(ibl)%degenerate.AND.ix*jx==0) CYCLE
                      nnz=nnz+inq*jnq
                      DO iq=1,inq
                        iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)
                        ientst_block(id,iind)=ientst_block(id,iind)+jnq
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       adjust the count for degenerate blocks: add the diagonal
c       and twice the number of off-diagonals in the degenerate point's
c       rows (storage doesn't economize for symmetry).
c-----------------------------------------------------------------------
        IF (bl_spp(ibl)%degenerate) THEN
          nnz2=nnz2+nqty**2
          IF (bl_spp(ibl)%row_ind(1)%rarr(1,0,my)>=
     $        spp%irowst_block(id-1)) THEN
            nnz=nnz+nqty**2
            DO iq=1,nqty
              iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
              ientst_block(id,iind)=ientst_block(id,iind)+nqty
            ENDDO
          ENDIF
          DO ib=1,nbmax
            iy0=mat_info(ibl)%iy0(ib)
            inq=mat_info(ibl)%nq_type(ib)
            DO iy=MAX(yst,iy0),my
              nnz2=nnz2+2*nqty*inq
              IF (bl_spp(ibl)%row_ind(ib)%rarr(1,1,iy)>=
     $            spp%irowst_block(id-1)) THEN
                nnz=nnz+2*nqty*inq
                DO iq=1,nqty
                  iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
                  ientst_block(id,iind)=ientst_block(id,iind)+inq
                ENDDO
                DO iq=1,inq
                  iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,1,iy)
                  ientst_block(id,iind)=ientst_block(id,iind)+nqty
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      IF (nprocs_layer>1) THEN
        CALL mpi_allreduce(nnz,iind,1_i4,mpi_nim_int,
     $                     mpi_sum,comm_layer,ierr)
        nnz=iind
        CALL mpi_allreduce(nbw,iind,1_i4,mpi_nim_int,
     $                     mpi_max,comm_layer,ierr)
        nbw=iind
      ENDIF
c-----------------------------------------------------------------------
c     determine where each block contributes to the compressed row
c     format by accumulating the number of unique entries.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        ALLOCATE(tmp2(0:nb_total,0:spp%nrow-1))
        CALL mpi_allreduce(ientst_block(0,0),tmp2(0,0),SIZE(tmp2),
     $                     mpi_nim_int,mpi_sum,comm_layer,ierr)
        ientst_block=tmp2
        DEALLOCATE(tmp2)
      ENDIF
      nq=0
      DO iind=0,spp%nrow-1
        DO id=0,nb_total
          nq=nq+ientst_block(id,iind)
          ientst_block(id,iind)=nq
        ENDDO
      ENDDO
      IF (ientst_block(nb_total,spp%nrow-1)/=nnz)
     $  CALL nim_stop("Alloc_sparsity_pattern: miscount on # of"//
     $                " entries.")
c-----------------------------------------------------------------------
c     allocate space for the nimrod column index and the starting index
c     for each row.  space for the array element values will be created
c     and destroyed on the fly to conserve memory.
c-----------------------------------------------------------------------
      spp%nnz=nnz
      spp%nbw=nbw
      ALLOCATE(spp%j_acc(0:nnz-1))
      ALLOCATE(spp%start_acc(0:spp%nrow))
      spp%start_acc(0:spp%nrow-1)=ientst_block(0,:)
      spp%start_acc(spp%nrow)=nnz
c-----------------------------------------------------------------------
c     temporarily allocate and fill row and column lists for all matrix
c     elements with contributions from this block.
c-----------------------------------------------------------------------
      ALLOCATE(irw(nnz2))
      ALLOCATE(spp%jentry(nnz2))
c-----------------------------------------------------------------------
c     the sparsity pattern is determined by the mesh and the number of
c     vector components, so it can be created here and saved.
c
c     here, we use nimrod's organization by considering
c     start_acc(irow) as the starting entry-index for row irow
c     and j_acc(ientry) as the column index for nonzero entry
c     ientry.
c-----------------------------------------------------------------------
      ALLOCATE(tmp1(-nbw:nbw))
      spp%j_acc=-1
      nnz2=0
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        IF (bl_spp(ibl)%perblock) THEN
          yst=1
        ELSE
          yst=0
        ENDIF
        DO ix=0,mx
          IF (bl_spp(ibl)%degenerate) THEN
            IF (ix==0) CYCLE
            jdeg=bl_spp(ibl)%row_ind(1)%rarr(1,0,my)
          ELSE
            jdeg=-1
          ENDIF
          DO iy=yst,my
            DO ib=nbmax,1,-1
              ix0=mat_info(ibl)%ix0(ib)
              iy0=mat_info(ibl)%iy0(ib)
              IF (ix0>ix.OR.iy0>iy) CYCLE
              DO iq=1,mat_info(ibl)%nq_type(ib)
                iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)
                ientry=ientst_block(id-1,iind)
                tmp1=-1
                DO jb=1,nbmax
                  jx0=mat_info(ibl)%ix0(jb)
                  jy0=mat_info(ibl)%iy0(jb)
                  jnq=mat_info(ibl)%nq_type(jb)
                  DO ly=jy0-1,1-iy0
                    jy=iy+ly
                    IF (bl_spp(ibl)%perblock.AND.jy==my+1) jy=1
                    IF (jy<jy0.OR.jy>my) CYCLE
                    DO lx=jx0-1,1-ix0
                      jx=ix+lx
                      IF (jx<jx0.OR.jx>mx) CYCLE
                      DO jq=1,jnq
                        jind=bl_spp(ibl)%row_ind(jb)%rarr(jq,jx,jy)
                        nnz2=nnz2+1
                        irw(nnz2)=iind
                        spp%jentry(nnz2)=jind
                        IF (iind>=spp%irowst_block(id-1).OR.
     $                      jind>=spp%irowst_block(id-1).OR.
     $                      concave_corner(lx,ly,ix0,iy0,jx0,jy0).AND.
     $                      bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)/=
     $                      jdeg) THEN
                          off=jind-iind
                          tmp1(off)=jind
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                DO jq=-nbw,nbw
                  IF (tmp1(jq)>=0) THEN
                    IF (ientry>=ientst_block(id,iind)) THEN
                      WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                  "Alloc_sparsity_pattern: miscounted ientry ",
     $                  ientry," nim row ",iind," in block ordered ",id
                      CALL nim_stop(msg)
                    ENDIF
                    IF (spp%j_acc(ientry)>=0) THEN
                      WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                  "Alloc_sparsity_pattern: re-counted ientry ",
     $                  ientry," nim col ",tmp1(jq),
     $                  " in block ordered",id
                      CALL nim_stop(msg)
                    ENDIF
                    spp%j_acc(ientry)=tmp1(jq)
                    ientry=ientry+1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         the degenerate point has only nqty distinct rows, but the loop
c         over iy is needed to get all columns.
c-----------------------------------------------------------------------
          IF (ix==1.AND.bl_spp(ibl)%degenerate) THEN
            DO iq=1,nqty
              tmp1=-1
              iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
              ientry=ientst_block(id-1,iind)
              DO jy=yst,my
                DO jb=nbmax,1,-1
                  jy0=mat_info(ibl)%iy0(jb)
                  jnq=mat_info(ibl)%nq_type(jb)
                  IF (jy<jy0) CYCLE
                  DO jq=1,jnq
                    jind=bl_spp(ibl)%row_ind(jb)%rarr(jq,1,jy)
                    nnz2=nnz2+1
                    irw(nnz2)=iind
                    spp%jentry(nnz2)=jind
                    IF (iind>=spp%irowst_block(id-1).OR.
     $                  jind>=spp%irowst_block(id-1)) THEN
                      off=jind-iind
                      tmp1(off)=jind
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              DO jq=1,nqty
                jind=bl_spp(ibl)%row_ind(1)%rarr(jq,0,my)
                nnz2=nnz2+1
                irw(nnz2)=iind
                spp%jentry(nnz2)=jind
                IF (iind>=spp%irowst_block(id-1).OR.
     $              jind>=spp%irowst_block(id-1)) THEN
                  off=jind-iind
                  tmp1(off)=jind
                ENDIF
              ENDDO
              DO jq=-nbw,nbw
                IF (tmp1(jq)>=0) THEN
                  IF (ientry>=ientst_block(id,iind)) THEN
                    WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                "Alloc_sparsity_pattern: deg miscounted ientry ",
     $                ientry," nim row ",iind," in block ordered ",id
                    CALL nim_stop(msg)
                  ENDIF
                  IF (spp%j_acc(ientry)>=0) THEN
                    WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                "Alloc_sparsity_pattern: deg re-counted ientry ",
     $                ientry," nim col ",tmp1(jq),
     $                " in block ordered ",id
                    CALL nim_stop(msg)
                  ENDIF
                  spp%j_acc(ientry)=tmp1(jq)
                  ientry=ientry+1
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(tmp1)
      DEALLOCATE(ientst_block)
c-----------------------------------------------------------------------
c     communicate j_acc.  to avoid a large allreduce, communicate the
c     number of nonzero contributions, their indices, and their values
c     via gather operations.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        ALLOCATE(nentry(0:nprocs_layer-1))
        ALLOCATE(ndispl(0:nprocs_layer))
        nq=SUM(MIN(spp%j_acc+1_i4,1_i4))
        ALLOCATE(sbd_jacc(nq))
        ALLOCATE(sbd_ind(nq))
        ientry=1
        DO iq=0,nnz-1
          IF (spp%j_acc(iq)>-1) THEN
            sbd_jacc(ientry)=spp%j_acc(iq)
            sbd_ind(ientry)=iq
            ientry=ientry+1
          ENDIF
        ENDDO
        CALL mpi_allgather(nq,1,mpi_nim_int,nentry(0),1,mpi_nim_int,
     $                     comm_layer,ierr)

        ndispl(0)=0
        DO id=1,nprocs_layer
          ndispl(id)=ndispl(id-1)+nentry(id-1)
        ENDDO
        jnq=SUM(nentry)
        ALLOCATE(glb_jacc(jnq))
        ALLOCATE(glb_ind(jnq))
        CALL mpi_allgatherv(sbd_jacc(1),nq,mpi_nim_int,glb_jacc(1),
     $       nentry(0),ndispl(0),mpi_nim_int,comm_layer,ierr)
        CALL mpi_allgatherv(sbd_ind(1),nq,mpi_nim_int,glb_ind(1),
     $       nentry(0),ndispl(0),mpi_nim_int,comm_layer,ierr)
        DO iq=1,jnq
          spp%j_acc(glb_ind(iq))=glb_jacc(iq)
        ENDDO
        DEALLOCATE(nentry,ndispl,sbd_jacc,sbd_ind,glb_jacc,glb_ind)
      ENDIF

      IF (MINVAL(spp%j_acc)<0.AND.nprocs_layer==1)
     $  CALL nim_stop("Alloc_sparsity_pattern: unfilled j_acc")
c-----------------------------------------------------------------------
c     nimrod needs to keep copies of j_acc and start_acc that are not
c     disturbed by SuperLU.
c
c     also use this space to trim duplicate entries in the irw and
c     jentry arrays before the final copy.
c-----------------------------------------------------------------------
      ALLOCATE(spp%j_acc_save(0:nnz-1))
      ALLOCATE(spp%start_acc_save(0:spp%nrow))

      spp%j_acc_save=0
      DO iq=1,nnz2
        iind=irw(iq)
        jind=spp%jentry(iq)
        DO jq=spp%start_acc(iind),spp%start_acc(iind+1)
          IF (jq>=spp%start_acc(iind+1)) THEN
            WRITE(msg,'(a,i8,a,i7,a,i8)')
     $      "Alloc_sparsity_pattern: jentry miscounted jq ",
     $      jq," nim row ",iind," cc ind ",jind
            CALL nim_stop(msg)
          ELSE IF (spp%j_acc(jq)==jind) THEN
            spp%j_acc_save(jq)=1
            EXIT
          ENDIF
        ENDDO
      ENDDO
      nnz2=SUM(spp%j_acc_save)
      DEALLOCATE(irw,spp%jentry)
      ALLOCATE(spp%jentry(nnz2),irw(nnz2))

      iq=0
      ix=0
      DO jq=0,nnz-1
        IF (jq==spp%start_acc(ix+1)) ix=ix+1
        IF (spp%j_acc_save(jq)==1) THEN
          iq=iq+1
          spp%jentry(iq)=jq
          irw(iq)=ix
        ENDIF
      ENDDO

      spp%j_acc_save(:)=spp%j_acc(:)
      spp%start_acc_save(:)=spp%start_acc(:)
c-----------------------------------------------------------------------
c     save information for distributed-memory storage of rows owned
c     by this processor. otherwise, have the local limits duplicate
c     the global limits.
c-----------------------------------------------------------------------
      IF (spp%matrix_distributed) THEN
        spp%fstrow=spp%irowst_block(spp%gblock_order(loc2glob(1))-1)
        spp%lstrow=spp%irowst_block(spp%gblock_order(loc2glob(nrb)))-1
        spp%mloc=0
        DO ibl=1,nrb
          id=spp%gblock_order(loc2glob(ibl))
          spp%mloc=spp%mloc+spp%irowst_block(id)-spp%irowst_block(id-1)
        ENDDO
        spp%indst=spp%start_acc_save(spp%fstrow)
        spp%indend=spp%start_acc_save(spp%fstrow+spp%mloc)-1
        spp%nnzloc=spp%indend-spp%indst+1
        ALLOCATE(spp%start_loc(0:spp%mloc))
        spp%start_loc=spp%start_acc(spp%fstrow:spp%fstrow+spp%mloc)-
     $                spp%start_acc(spp%fstrow)
      ELSE
        spp%fstrow=0
        spp%lstrow=spp%nrow-1
        spp%mloc=spp%nrow
        spp%indst=0
        spp%indend=spp%nnz-1
        spp%nnzloc=spp%nnz
      ENDIF
c-----------------------------------------------------------------------
c     acquire communication information for distributed-memory storage.
c-----------------------------------------------------------------------
      CALL pardir_setup(spp,nrb,nnz2,irw,spp%jentry)
      DEALLOCATE(irw)
      IF (spp%nsnd>0) DEALLOCATE(spp%sendrowst)
c-----------------------------------------------------------------------
c     if using distributed storage, we are finished.
c-----------------------------------------------------------------------
      IF (spp%matrix_distributed) RETURN
c-----------------------------------------------------------------------
c     generate pointer information for full-storage allgathers.
c-----------------------------------------------------------------------
      ALLOCATE(spp%algcount(0:nprocs_layer-1),
     $         spp%algdispl(0:nprocs_layer),
     $         spp%algcntr(0:nprocs_layer-1),
     $         spp%algdsplr(0:nprocs_layer))
      spp%algdispl=0
      spp%algcount=0
      spp%algdsplr=0
      spp%algcntr=0
      DO ibl=1,nb_total
        ip=block2proc(ibl)-ilayer*nprocs_layer
        id=spp%gblock_order(ibl)
        spp%algcount(ip)=spp%algcount(ip)+
     $                   spp%start_acc(spp%irowst_block(id))-
     $                   spp%start_acc(spp%irowst_block(id-1))
        spp%algcntr(ip)=spp%algcntr(ip)+
     $                  spp%irowst_block(id)-
     $                  spp%irowst_block(id-1)
      ENDDO
      DO ip=1,nprocs_layer
        spp%algdispl(ip)=spp%algdispl(ip-1)+spp%algcount(ip-1)
        spp%algdsplr(ip)=spp%algdsplr(ip-1)+spp%algcntr(ip-1)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE alloc_sparsity_pattern
c-----------------------------------------------------------------------
c     subprogram 12. alloc_sparsity_pattern_dist.
c     distributed allocation routine for the sparsity pattern. this
c     routine picks up after step 3 of assign_node_indices and
c     continues to initialize and allocate the sparsity_pattern spp
c     structure.
c
c     the following additional steps are done for the compressed row
c     storage:
c
c       4. determine the number of nonzero entries in each row of the
c          matrix.  a sum that is decomposed by block is placed in
c            ientst_block(0:nbl,0:spp%nrow-1)
c          and the running total is placed in
c            ientst_block(0,0:spp%nrow-1) and
c            spp%start_acc(0:spp%nrow-1)
c          the running total is equivalent to the index where
c          each row starts in the compressed-row storage for the matrix.
c          this step allocates arrays of size spp%nrow, and thus may
c          become a future memory scaling barrier.
c       5. rows owned by this processor are defined.  spp%mloc is the
c          number of these rows and spp%start_loc(0:spp%mloc-1) is the
c          starting location in compressed-row storage for these rows
c          only. start_acc is restricted to
c            spp%start_acc(spp%spp%fstrow:spp%lstrow+1)
c       6. the column indices for each nonzero matrix entry are
c          determined row-by-row.  they are placed in the array
c            spp%j_acc(spp%indst:spp%indend)
c          which is resricted to the node local values.
c       7. lists of compressed row indices and rows for all matrix
c          elements with contributions from blocks on this processor
c          are also created and saved in
c            spp%jentry(:)
c            irw(:)
c          to assist communication for distributed-memory solves.
c
c-----------------------------------------------------------------------
      SUBROUTINE alloc_sparsity_pattern_dist(mat_info,spp,bl_spp,
     $                                       nrb,nqty,nbmax)
      USE pardata
      USE mpi_nim
      USE factor_type_mod
      USE io, ONLY: nim_wr

      TYPE(rbl_mat_info),     INTENT(IN)    :: mat_info(:)
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      TYPE(bl_sparsity_type), INTENT(INOUT) :: bl_spp(:)
      INTEGER(i4), INTENT(IN) :: nrb,nqty,nbmax

      INTEGER(i4) :: iq,jq,lx,ly,ix,iy,ix0,iy0,jx0,jy0,nq,ib,jb,n,
     $               nnz,inq,jnq,jx,jy,iind,jind,ientry,nbw,off,yst,
     $               ibl,mx,my,ierr,id,nnz2,iindhi,iindlo,isnd,nsnd,
     $               recvhi,jdeg,icomm,nb_total,nnzsnd,ircv,iv,jv,ns,
     $               sz,pd
      INTEGER(i4), PARAMETER :: chunk_size=2000
      INTEGER(i4), ALLOCATABLE :: recv_ientst(:),send_ientst(:),tmp1(:),
     $                            j_acc_comm_val(:),j_acc_comm_rind(:),
     $                            j_acc_comm_ientry(:),irw(:),tmpn(:,:),
     $                            recvrowst(:),recvirw(:,:),
     $                            ientst_block(:,:)
      TYPE :: comm_data
        INTEGER(i4), ALLOCATABLE :: data(:)
      END TYPE comm_data
      TYPE(comm_data), ALLOCATABLE :: recvstr(:),sendstr(:)
      INTEGER(i4), DIMENSION(mpi_status_size) :: mpi_status
      LOGICAL :: sdone
      INTEGER(i4) :: mpi_request
      CHARACTER(128) :: msg
      LOGICAL, EXTERNAL :: concave_corner

c-----------------------------------------------------------------------
c     write out a warning if ientst_block will be larger than 100Mb.
c     this will help identify if this is a problem for future scaling.
c     the mpi barrier prevents an OOM event before the write.
c-----------------------------------------------------------------------
      IF (spp%nrow*nrb*4_i4 > 100000000_i4) THEN
        IF (node==0) THEN
          WRITE(nim_wr,'(a,i4,a,i8,a,i6,a)')
     $                    'Alloc_sparsity_pattern_dist: Warning: the '
     $                  //'temporary array ientst_block of size '
     $                  //'(nbl+1)*nrows (nrb=',nrb,'nrows=',spp%nrow,
     $                    ') will use ',spp%nrow*(nrb+1)*4_i4/
     $                    1000000_i4,'Mb of space'
        ENDIF
        CALL mpi_barrier(comm_layer,ierr)
      ENDIF

c-----------------------------------------------------------------------
c     determine the total number of nonzero entries that are possible.
c
c     ientst_block will have the starting entry index of the
c     compressed row format for the nonzero entries in each unique
c     nimrod-matrix row, organized by block number.
c-----------------------------------------------------------------------
      nb_total=SIZE(global2local)
      nbw=0
      nnz=0
      nnz2=0
      ALLOCATE(ientst_block(0:nrb,0:spp%nrow-1))
      ientst_block=0
      ALLOCATE(spp%start_acc(0:spp%nrow))
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        IF (bl_spp(ibl)%perblock) THEN
          yst=1
        ELSE
          yst=0
        ENDIF
        DO ib=1,nbmax
          ix0=mat_info(ibl)%ix0(ib)
          iy0=mat_info(ibl)%iy0(ib)
          inq=mat_info(ibl)%nq_type(ib)
          DO jb=1,nbmax
            jx0=mat_info(ibl)%ix0(jb)
            jy0=mat_info(ibl)%iy0(jb)
            jnq=mat_info(ibl)%nq_type(jb)
            DO ix=ix0,mx
              DO iy=MAX(yst,iy0),my
                DO lx=jx0-1,1-ix0
                  jx=ix+lx
                  IF (jx<jx0.OR.jx>mx) CYCLE
                  DO ly=jy0-1,1-iy0
                    jy=iy+ly
                    IF (bl_spp(ibl)%perblock.AND.jy==my+1) jy=1
                    IF (jy<jy0.OR.jy>my) CYCLE
                    nnz2=nnz2+inq*jnq
                    IF (bl_spp(ibl)%row_ind(ib)%rarr(1,ix,iy)>=
     $                  spp%irowst_block(id-1).OR.
     $                  bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)>=
     $                  spp%irowst_block(id-1).OR.
     $                  concave_corner(lx,ly,ix0,iy0,jx0,jy0)) THEN
                      off=bl_spp(ibl)%row_ind(ib)%rarr(inq,ix,iy)-
     $                    bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)
                      nbw=MAX(nbw,ABS(off))
                      off=bl_spp(ibl)%row_ind(ib)%rarr(1,ix,iy)-
     $                    bl_spp(ibl)%row_ind(jb)%rarr(jnq,jx,jy)
                      nbw=MAX(nbw,ABS(off))
                      IF (bl_spp(ibl)%degenerate.AND.ix*jx==0) CYCLE
                      nnz=nnz+inq*jnq
                      DO iq=1,inq
                        iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)
                        ientst_block(ibl,iind)=
     $                    ientst_block(ibl,iind)+jnq
                      ENDDO
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       adjust the count for degenerate blocks: add the diagonal
c       and twice the number of off-diagonals in the degenerate point's
c       rows (storage doesn't economize for symmetry).
c-----------------------------------------------------------------------
        IF (bl_spp(ibl)%degenerate) THEN
          nnz2=nnz2+nqty**2
          IF (bl_spp(ibl)%row_ind(1)%rarr(1,0,my)>=
     $        spp%irowst_block(id-1)) THEN
            nnz=nnz+nqty**2
            DO iq=1,nqty
              iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
              ientst_block(ibl,iind)=ientst_block(ibl,iind)+nqty
            ENDDO
          ENDIF
          DO ib=1,nbmax
            iy0=mat_info(ibl)%iy0(ib)
            inq=mat_info(ibl)%nq_type(ib)
            DO iy=MAX(yst,iy0),my
              nnz2=nnz2+2*nqty*inq
              IF (bl_spp(ibl)%row_ind(ib)%rarr(1,1,iy)>=
     $            spp%irowst_block(id-1)) THEN
                nnz=nnz+2*nqty*inq
                DO iq=1,nqty
                  iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
                  ientst_block(ibl,iind)=ientst_block(ibl,iind)+inq
                ENDDO
                DO iq=1,inq
                  iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,1,iy)
                  ientst_block(ibl,iind)=ientst_block(ibl,iind)+nqty
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      IF (nprocs_layer>1) THEN
        CALL mpi_allreduce(nnz,iind,1_i4,mpi_nim_int,
     $                     mpi_sum,comm_layer,ierr)
        nnz=iind
        CALL mpi_allreduce(nbw,iind,1_i4,mpi_nim_int,
     $                     mpi_max,comm_layer,ierr)
        nbw=iind
      ENDIF
c-----------------------------------------------------------------------
c     determine where each block contributes to the compressed row
c     format by accumulating the number of unique entries. the
c     mpi sends are done in chunks of chunk_size.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        nsnd=spp%nrow/chunk_size+1
        IF (MOD(spp%nrow,chunk_size)==0) nsnd=nsnd-1
        ALLOCATE(recv_ientst(0:chunk_size-1))
        recvhi=chunk_size-1
        DO isnd=1,nsnd
          recv_ientst=0
          IF (node_layer/=0) THEN
            CALL mpi_recv(recv_ientst(0),chunk_size,mpi_nim_int,
     $                    node_layer-1,0_i4,comm_layer,mpi_status,ierr)
          ENDIF
          iindlo=(isnd-1)*chunk_size
          iindhi=isnd*chunk_size-1
          IF(isnd==nsnd) THEN
            recvhi=MOD(spp%nrow,chunk_size)-1
            iindhi=iindlo+MOD(spp%nrow,chunk_size)-1
          ENDIF
          DO ibl=1,nrb
            recv_ientst(0:recvhi)=recv_ientst(0:recvhi)
     $                           +ientst_block(ibl,iindlo:iindhi)
            ientst_block(ibl,iindlo:iindhi)=recv_ientst(0:recvhi)
          ENDDO
          IF (node_layer/=nprocs_layer-1) THEN
            CALL mpi_send(recv_ientst(0),chunk_size,mpi_nim_int,
     $                    node_layer+1,0_i4,comm_layer,ierr)
          ELSE
            spp%start_acc(iindlo+1:iindhi+1)=
     $                                 ientst_block(nrb,iindlo:iindhi)
          ENDIF
        ENDDO
        DEALLOCATE(recv_ientst)
        spp%start_acc(0)=0
        CALL mpi_bcast(spp%start_acc(0),spp%nrow,mpi_nim_int,
     $                 nprocs_layer-1,comm_layer,ierr)
        DO iind=1,spp%nrow-1
          spp%start_acc(iind)=spp%start_acc(iind)+spp%start_acc(iind-1)
          ientst_block(1:nrb,iind)=ientst_block(1:nrb,iind)
     $                            +spp%start_acc(iind)
        ENDDO
c-----------------------------------------------------------------------
c       communicate the previous block row contributions to the
c       next process.
c-----------------------------------------------------------------------
        ALLOCATE(recv_ientst(0:spp%nrow-1),send_ientst(0:spp%nrow-1))
        mpi_request=mpi_request_null
        IF (node_layer/=0) THEN
          CALL mpi_irecv(recv_ientst(0),spp%nrow,mpi_nim_int,
     $                   node_layer-1,0_i4,comm_layer,mpi_request,ierr)
        ENDIF
        IF (node_layer/=nprocs_layer-1) THEN
          send_ientst=ientst_block(nrb,:)
          CALL mpi_send(send_ientst(0),spp%nrow,mpi_nim_int,
     $                  node_layer+1,0_i4,comm_layer,ierr)
        ENDIF
        IF (node_layer/=0) THEN
          CALL mpi_wait(mpi_request,mpi_status,ierr)
          ientst_block(0,:)=recv_ientst(0:spp%nrow-1)
        ELSE
          ientst_block(0,:)=spp%start_acc(0:spp%nrow-1)
        ENDIF
        DEALLOCATE(recv_ientst,send_ientst)
      ELSE
        nq=0
        DO iind=0,spp%nrow-1
          DO ibl=0,nb_total
            nq=nq+ientst_block(ibl,iind)
            ientst_block(ibl,iind)=nq
          ENDDO
        ENDDO
        spp%start_acc(0:spp%nrow-1)=ientst_block(0,:)
      ENDIF
      spp%start_acc(spp%nrow)=nnz
      IF (node_layer==nprocs_layer-1) THEN
        IF (ientst_block(nrb,spp%nrow-1)/=nnz)
     $    CALL nim_stop("Alloc_sparsity_pattern: miscount on # of "//
     $                  "entries.")
      ENDIF
c-----------------------------------------------------------------------
c     restrict the row pointer array (start_acc) to the local rows
c     owned by this processor (start_loc). set the row bounds.
c-----------------------------------------------------------------------
      id=spp%gblock_order(loc2glob(1))
      spp%fstrow=spp%irowst_block(id-1)
      id=spp%gblock_order(loc2glob(nrb))
      spp%lstrow=spp%irowst_block(id)-1
      spp%mloc=spp%lstrow-spp%fstrow+1
      spp%indst=spp%start_acc(spp%fstrow)
      spp%indend=spp%start_acc(spp%lstrow+1)-1
      ALLOCATE(spp%start_loc(0:spp%mloc))
      spp%start_loc=spp%start_acc(spp%fstrow:spp%lstrow+1)
      DEALLOCATE(spp%start_acc)
      ALLOCATE(spp%start_acc(spp%fstrow:spp%lstrow+1))
      spp%start_acc=spp%start_loc(0:spp%mloc)
      spp%start_loc=spp%start_loc-spp%start_acc(spp%fstrow)
c-----------------------------------------------------------------------
c     allocate space for the nimrod column index (j_acc), only store
c     the local values. space for the array element values will be
c     created and destroyed on the fly to conserve memory.
c-----------------------------------------------------------------------
      spp%nnz=nnz
      spp%nnzloc=spp%indend-spp%indst+1
      spp%nbw=nbw
      ALLOCATE(spp%j_acc(spp%indst:spp%indend))
      spp%j_acc=-1
c-----------------------------------------------------------------------
c     temporarily allocate and fill row and column lists for all matrix
c     elements with contributions from this processor.
c-----------------------------------------------------------------------
      ALLOCATE(irw(nnz2),spp%jentry(nnz2))
      irw=-1
      ALLOCATE(j_acc_comm_ientry(nnz2),j_acc_comm_val(nnz2))
      ALLOCATE(j_acc_comm_rind(nnz2))
      j_acc_comm_rind=-1
c TMP
c      DO jq=0,nprocs_layer-1
c        call mpi_barrier(comm_layer,ierr)
c        if (jq==node_layer) then
c          write(*,*) node_layer,nnz,nnz2,nbw
c          write(*,*) node_layer,spp%nrow,spp%fstrow,spp%lstrow
c          write(*,*) node_layer,spp%nnzloc,spp%indst,spp%indend,
c     $               nnz2-spp%nnzloc
c          write(*,*) '--------'
c        endif
c      enddo
c-----------------------------------------------------------------------
c     the sparsity pattern is determined by the mesh and the number of
c     vector components, so it can be created here and saved.
c
c     here, we use nimrod's organization by considering
c     start_acc(irow) as the starting entry-index for row irow
c     and j_acc(ientry) as the column index for nonzero entry ientry.
c-----------------------------------------------------------------------
      ALLOCATE(tmp1(-nbw:nbw))
      nnz2=0
      icomm=0
      DO ibl=1,nrb
        mx=mat_info(ibl)%mx
        my=mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        IF (bl_spp(ibl)%perblock) THEN
          yst=1
        ELSE
          yst=0
        ENDIF
        DO ix=0,mx
          IF (bl_spp(ibl)%degenerate) THEN
            IF (ix==0) CYCLE
            jdeg=bl_spp(ibl)%row_ind(1)%rarr(1,0,my)
          ELSE
            jdeg=-1
          ENDIF
          DO iy=yst,my
            DO ib=nbmax,1,-1
              ix0=mat_info(ibl)%ix0(ib)
              iy0=mat_info(ibl)%iy0(ib)
              IF (ix0>ix.OR.iy0>iy) CYCLE
              DO iq=1,mat_info(ibl)%nq_type(ib)
                iind=bl_spp(ibl)%row_ind(ib)%rarr(iq,ix,iy)
                ientry=ientst_block(ibl-1,iind)
                tmp1=-1
                DO jb=1,nbmax
                  jx0=mat_info(ibl)%ix0(jb)
                  jy0=mat_info(ibl)%iy0(jb)
                  jnq=mat_info(ibl)%nq_type(jb)
                  DO ly=jy0-1,1-iy0
                    jy=iy+ly
                    IF (bl_spp(ibl)%perblock.AND.jy==my+1) jy=1
                    IF (jy<jy0.OR.jy>my) CYCLE
                    DO lx=jx0-1,1-ix0
                      jx=ix+lx
                      IF (jx<jx0.OR.jx>mx) CYCLE
                      DO jq=1,jnq
                        jind=bl_spp(ibl)%row_ind(jb)%rarr(jq,jx,jy)
                        IF (iind < spp%fstrow) THEN
                          nnz2=nnz2+1
                          irw(nnz2)=iind
                          spp%jentry(nnz2)=jind
                        ENDIF
                        IF (iind>=spp%irowst_block(id-1).OR.
     $                      jind>=spp%irowst_block(id-1).OR.
     $                      concave_corner(lx,ly,ix0,iy0,jx0,jy0).AND.
     $                      bl_spp(ibl)%row_ind(jb)%rarr(1,jx,jy)/=
     $                      jdeg) THEN
                          off=jind-iind
                          tmp1(off)=jind
                        ENDIF
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
                DO jq=-nbw,nbw
                  IF (tmp1(jq)>=0) THEN
                    IF (ientry>=ientst_block(ibl,iind)) THEN
                      WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                  "Alloc_sparsity_pattern: miscounted ientry ",
     $                  ientry," nim row ",iind," in block ordered ",id
                      CALL nim_stop(msg)
                    ENDIF
                    IF (ientry>=spp%indst.AND.ientry<=spp%indend) THEN
                      IF (spp%j_acc(ientry)>=0) THEN
                        WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                    "Alloc_sparsity_pattern: re-counted ientry ",
     $                    ientry," nim col ",tmp1(jq),
     $                    " in block ordered ",id
                        CALL nim_stop(msg)
                      ENDIF
                      spp%j_acc(ientry)=tmp1(jq)
                    ELSE
                      icomm=icomm+1
                      j_acc_comm_ientry(icomm)=ientry
                      j_acc_comm_val(icomm)=tmp1(jq)
                      j_acc_comm_rind(icomm)=iind
                    ENDIF
                    ientry=ientry+1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         the degenerate point has only nqty distinct rows, but the loop
c         over iy is needed to get all columns.
c-----------------------------------------------------------------------
          IF (ix==1.AND.bl_spp(ibl)%degenerate) THEN
            DO iq=1,nqty
              tmp1=-1
              iind=bl_spp(ibl)%row_ind(1)%rarr(iq,0,my)
              ientry=ientst_block(ibl-1,iind)
              DO jy=yst,my
                DO jb=nbmax,1,-1
                  jy0=mat_info(ibl)%iy0(jb)
                  jnq=mat_info(ibl)%nq_type(jb)
                  IF (jy<jy0) CYCLE
                  DO jq=1,jnq
                    jind=bl_spp(ibl)%row_ind(jb)%rarr(jq,1,jy)
                    IF (iind < spp%fstrow) THEN
                      nnz2=nnz2+1
                      irw(nnz2)=iind
                      spp%jentry(nnz2)=jind
                    ENDIF
                    IF (iind>=spp%irowst_block(id-1).OR.
     $                  jind>=spp%irowst_block(id-1)) THEN
                      off=jind-iind
                      tmp1(off)=jind
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
              DO jq=1,nqty
                jind=bl_spp(ibl)%row_ind(1)%rarr(jq,0,my)
                IF (iind < spp%fstrow) THEN
                  nnz2=nnz2+1
                  irw(nnz2)=iind
                  spp%jentry(nnz2)=jind
                ENDIF
                IF (iind>=spp%irowst_block(id-1).OR.
     $              jind>=spp%irowst_block(id-1)) THEN
                  off=jind-iind
                  tmp1(off)=jind
                ENDIF
              ENDDO
              DO jq=-nbw,nbw
                IF (tmp1(jq)>=0) THEN
                  IF (ientry>=ientst_block(ibl,iind)) THEN
                    WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                "Alloc_sparsity_pattern: deg miscounted ientry ",
     $                ientry," nim row ",iind," in block ordered ",id
                    CALL nim_stop(msg)
                  ENDIF
                  IF (ientry>=spp%indst.AND.ientry<=spp%indend) THEN
                    IF (spp%j_acc(ientry)>=0) THEN
                      WRITE(msg,'(a,i7,a,i7,a,i3)')
     $                  "Alloc_sparsity_pattern: deg re-counted ientry",
     $                  ientry," nim col ",tmp1(jq),
     $                  " in block ordered ",id
                      CALL nim_stop(msg)
                    ENDIF
                    spp%j_acc(ientry)=tmp1(jq)
                  ELSE
                    icomm=icomm+1
                    j_acc_comm_ientry(icomm)=ientry
                    j_acc_comm_val(icomm)=tmp1(jq)
                    j_acc_comm_rind(icomm)=iind
                  ENDIF
                  ientry=ientry+1
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(tmp1)
      DEALLOCATE(ientst_block)
c-----------------------------------------------------------------------
c     sort the j_acc_comm_* arrays by row.
c-----------------------------------------------------------------------
      nnzsnd=icomm
      IF (nprocs_layer>1) THEN
        CALL qsort(j_acc_comm_rind(1:nnzsnd),j_acc_comm_val(1:nnzsnd),
     $             j_acc_comm_ientry(1:nnzsnd))
c-----------------------------------------------------------------------
c       communicate j_acc entries owned by other processors.
c-----------------------------------------------------------------------
        CALL pardir_setup(spp,nrb,nnzsnd,j_acc_comm_rind(1:nnzsnd),
     $                         j_acc_comm_ientry(1:nnzsnd))

        IF (spp%nrcv>0) ALLOCATE(recvstr(spp%nrcv))
        IF (spp%nsnd>0) ALLOCATE(sendstr(spp%nsnd))
c-----------------------------------------------------------------------
c       post receives.
c-----------------------------------------------------------------------
        DO ircv=1,spp%nrcv
          ALLOCATE(recvstr(ircv)%data(spp%recvtot(ircv)))
          CALL mpi_irecv(recvstr(ircv)%data(1),spp%recvtot(ircv),
     $                   mpi_nim_int,spp%recvlist(ircv),0_i4,
     $                   comm_layer,spp%recv_req(ircv),ierr)
        ENDDO
c-----------------------------------------------------------------------
c       send j_acc values.
c-----------------------------------------------------------------------
        DO isnd=1,spp%nsnd
          ALLOCATE(sendstr(isnd)%data(spp%sendtot(isnd)))
          iv=spp%sendst(isnd)
          DO jv=1,spp%sendtot(isnd)
            sendstr(isnd)%data(jv)=j_acc_comm_val(iv)
            iv=iv+1
          ENDDO
          CALL mpi_isend(sendstr(isnd)%data(1),spp%sendtot(isnd),
     $                   mpi_nim_int,spp%sendlist(isnd),0_i4,
     $                   comm_layer,spp%send_req(isnd),ierr)
        ENDDO
c-----------------------------------------------------------------------
c       collect j_acc values from other processors as they become
c       available. then check on sends.
c-----------------------------------------------------------------------
        ns=0
        DO iv=1,spp%nrcv
          CALL mpi_waitany(spp%nrcv,spp%recv_req,ircv,mpi_status,ierr)
          DO jv=1,spp%recvtot(ircv)
            spp%j_acc(spp%recvjentry(jv,ircv))=recvstr(ircv)%data(jv)
          ENDDO
          DEALLOCATE(recvstr(ircv)%data)

          IF (spp%nsnd>0) THEN
            CALL mpi_testany(spp%nsnd,spp%send_req,isnd,sdone,
     $                       mpi_status,ierr)
            IF (sdone.AND.ns<spp%nsnd) THEN
              ns=ns+1
              DEALLOCATE(sendstr(isnd)%data)
            ENDIF
          ENDIF
        ENDDO

        DO iv=ns+1,spp%nsnd
          CALL mpi_waitany(spp%nsnd,spp%send_req,isnd,mpi_status,ierr)
          DEALLOCATE(sendstr(isnd)%data)
        ENDDO

        IF (spp%nrcv>0) DEALLOCATE(recvstr)
        IF (spp%nsnd>0) DEALLOCATE(sendstr)
c-----------------------------------------------------------------------
c       deallocate the send/recv storage space used to communicate
c       values of j_acc to other processors.
c-----------------------------------------------------------------------
        IF (spp%nrcv>0) THEN
          DEALLOCATE(spp%recvlist,spp%recvtot,
     $               spp%recvjentry,spp%recv_req)
        ENDIF
        IF (spp%nsnd>0) THEN
          DEALLOCATE(spp%sendlist,spp%sendtot,spp%sendst,
     $               spp%sendrowst,spp%send_req)
        ENDIF
      ENDIF

      DEALLOCATE(j_acc_comm_rind,j_acc_comm_ientry,j_acc_comm_val)

      IF (MINVAL(spp%j_acc)<0.AND.nprocs_layer==1)
     $  CALL nim_stop("Alloc_sparsity_pattern: unfilled j_acc")
c-----------------------------------------------------------------------
c     trim duplicate values then order irw and jentry by row index.
c-----------------------------------------------------------------------
      DO iq=1,nnz2
        IF (irw(iq)==-1) CYCLE
        DO jq=iq+1,nnz2
          IF (irw(iq)==irw(jq).AND.
     $        spp%jentry(iq)==spp%jentry(jq)) irw(jq)=-1
        ENDDO
      ENDDO

      jq=0
      DO iq=1,nnz2
        DO WHILE (iq+jq<=nnz2.AND.irw(iq+jq)==-1)
          jq=jq+1
        ENDDO
        IF (iq+jq>nnz2) EXIT
        IF (jq>0) THEN
          irw(iq)=irw(iq+jq)
          spp%jentry(iq)=spp%jentry(iq+jq)
        ENDIF
      ENDDO
      nnz2=nnz2-jq

      CALL qsort(irw(1:nnz2),spp%jentry(1:nnz2))
c-----------------------------------------------------------------------
c     acquire communication information for distributed-memory storage,
c     then assemble the sendstacc structure for loading the distributed
c     sends.
c-----------------------------------------------------------------------
      CALL pardir_setup(spp,nrb,nnz2,irw,spp%jentry)

      IF (spp%nsnd>0) THEN
        jq=1
        ALLOCATE(spp%sendstacc(spp%nsnd))
        DO isnd=1,spp%nsnd
          sz=1
          iind=-1
          DO iq=1,spp%sendtot(isnd)
            IF (irw(jq)/=iind) THEN
              iind=irw(jq)
              sz=sz+1
            ENDIF
            jq=jq+1
          ENDDO
          ALLOCATE(spp%sendstacc(isnd)%data(sz,2))
        ENDDO
        jq=1
        DO isnd=1,spp%nsnd
          sz=1
          iind=-1
          DO iq=1,spp%sendtot(isnd)
            IF (irw(jq)/=iind) THEN
              iind=irw(jq)
              spp%sendstacc(isnd)%data(sz,1)=iind
              spp%sendstacc(isnd)%data(sz,2)=jq
              sz=sz+1
            ENDIF
            jq=jq+1
          ENDDO
          spp%sendstacc(isnd)%data(sz,1)=-1
          spp%sendstacc(isnd)%data(sz,2)=jq
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     in this version spp%jentry is the column index, not a global
c     matrix index. convert spp%recvjentry to a global matrix index
c     for the distributed matrix loads of received data.
c-----------------------------------------------------------------------
      IF (spp%nrcv>0) THEN
        ALLOCATE(recvrowst(spp%nrcv))
        DO ircv=1,spp%nrcv
          CALL mpi_irecv(recvrowst(ircv),1_i4,mpi_nim_int,
     $                   spp%recvlist(ircv),0_i4,comm_layer,
     $                   spp%recv_req(ircv),ierr)
        ENDDO
      ENDIF

      DO isnd=1,spp%nsnd
        CALL mpi_send(spp%sendrowst(isnd),1_i4,mpi_nim_int,
     $                spp%sendlist(isnd),0_i4,comm_layer,ierr)
      ENDDO

      IF (spp%nrcv>0) THEN
        ALLOCATE(tmpn(mpi_status_size,spp%nrcv))
        CALL mpi_waitall(spp%nrcv,spp%recv_req,tmpn,ierr)

        ALLOCATE(recvirw(SIZE(spp%recvjentry,DIM=1),spp%nrcv))
        DO ircv=1,spp%nrcv
          CALL mpi_irecv(recvirw(1,ircv),spp%recvtot(ircv),
     $                   mpi_nim_int,spp%recvlist(ircv),0_i4,
     $                   comm_layer,spp%recv_req(ircv),ierr)
        ENDDO
      ENDIF

      DO isnd=1,spp%nsnd
        CALL mpi_send(irw(spp%sendst(isnd)),spp%sendtot(isnd),
     $                mpi_nim_int,spp%sendlist(isnd),0_i4,
     $                comm_layer,ierr)
      ENDDO

      IF (spp%nrcv>0) THEN
        CALL mpi_waitall(spp%nrcv,spp%recv_req,tmpn,ierr)
        DEALLOCATE(tmpn)

        DO ircv=1,spp%nrcv
          DO jv=1,spp%recvtot(ircv)
            iind=recvirw(jv,ircv)
            DO n=spp%start_acc(iind),spp%start_acc(iind+1)-1
              jind=spp%recvjentry(jv,ircv)
              IF (jind==spp%j_acc(n)) THEN
                spp%recvjentry(jv,ircv)=n
                EXIT
              ENDIF
              IF (n==spp%start_acc(iind+1)-1) THEN
                WRITE(msg,'(a,i7,a,i6,a,i6,a,i6)')
     $              'Alloc_sparsity_pattern_dist: spp%recvjentry not '//
     $              'found, jv=',jv,' ircv=',ircv,
     $              'iind=',iind,'jind=',jind
                CALL nim_stop(msg)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(recvrowst,recvirw)
      ENDIF
c-----------------------------------------------------------------------
c     resize jentry before deallocating irw.
c-----------------------------------------------------------------------
      irw(1:nnz2)=spp%jentry(1:nnz2)
      DEALLOCATE(spp%jentry)
      ALLOCATE(spp%jentry(nnz2))
      spp%jentry=irw(1:nnz2)
      DEALLOCATE(irw)
c-----------------------------------------------------------------------
c     create j_acc_save and start_acc_save and store saved arrays.
c-----------------------------------------------------------------------
      ALLOCATE(spp%j_acc_save(spp%indst:spp%indend))
      ALLOCATE(spp%start_acc_save(0:spp%mloc))
      spp%j_acc_save(:)=spp%j_acc(:)
      spp%start_acc_save(:)=spp%start_loc(:)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      CONTAINS
c-----------------------------------------------------------------------
c     recursive quick-sort algorithm.
c-----------------------------------------------------------------------
      PURE RECURSIVE SUBROUTINE qsort(arrs,arr1,arr2)
      INTEGER(i4), INTENT(INOUT) :: arrs(:)
      INTEGER(i4), INTENT(INOUT) :: arr1(:)
      INTEGER(i4), OPTIONAL, INTENT(INOUT) :: arr2(:)

      INTEGER(i4) :: iq

      IF (SIZE(arrs)>1) THEN
        IF (PRESENT(arr2)) THEN
          CALL partition(arrs,arr1,iq,arr2)
          CALL qsort(arrs(:iq-1),arr1(:iq-1),arr2(:iq-1))
          CALL qsort(arrs(iq:),  arr1(iq:),  arr2(iq:)  )
        ELSE
          CALL partition(arrs,arr1,iq)
          CALL qsort(arrs(:iq-1),arr1(:iq-1))
          CALL qsort(arrs(iq:),  arr1(iq:)  )
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE qsort

c-----------------------------------------------------------------------
c     partition into two sorted sub-arrays greater and less than the
c     pivot.
c-----------------------------------------------------------------------
      PURE SUBROUTINE partition(arrs,arr1,iq,arr2)
      INTEGER(i4), INTENT(INOUT) :: arrs(:)
      INTEGER(i4), INTENT(INOUT) :: arr1(:)
      INTEGER(i4), INTENT(OUT) :: iq    ! split index
      INTEGER(i4), OPTIONAL, INTENT(INOUT) :: arr2(:)

      INTEGER(i4) :: ii,jj,tmp,pivot

      pivot=arrs(SIZE(arrs)/2)
      ii=1
      jj=SIZE(arrs)
      DO WHILE (ii <= jj)
        DO WHILE (arrs(jj) > pivot)
          jj=jj-1
        ENDDO
        DO WHILE (arrs(ii) < pivot)
          ii=ii+1
        ENDDO
        IF (ii < jj) THEN ! exchange ii and jj
          tmp=arrs(ii);  arrs(ii)=arrs(jj);  arrs(jj)=tmp
          tmp=arr1(ii);  arr1(ii)=arr1(jj);  arr1(jj)=tmp
          IF (PRESENT(arr2)) THEN
            tmp=arr2(ii);  arr2(ii)=arr2(jj);  arr2(jj)=tmp
          ENDIF
          jj=jj-1
          ii=ii+1
        ELSEIF (ii==jj) THEN
          jj=jj-1
          ii=ii+1
        ENDIF
      ENDDO
      iq=ii
      END SUBROUTINE partition

      END SUBROUTINE alloc_sparsity_pattern_dist
c-----------------------------------------------------------------------
c     subprogram 13. pardir_setup.
c     collect information used for completing rows of the compressed
c     row storage when blocks are allocated to different processors.
c
c     the spp%gblock_order array keeps ranges of row indices together
c     on a processor, and on each processor, compressed row indices
c     increase with increasing local/global block number.
c
c     upon exit:
c       spp%nrcv is the number of other processors that contribute
c         to the rows owned by this processor.
c       spp%recvlist(1:spp%nrcv) is the list of those processor ranks.
c       spp%recvtot(1:spp%nrcv) is the number of potentially nonzero
c         matrix elements from each contributing processor.
c       spp%recvjentry(:,1:spp%nrcv) is the compressed-row index of each
c         datum from contributing processors.
c
c       spp%nsnd is the number of other processors that receive data
c         from this processor.
c       spp%sendlist(1:spp%nsnd) is the list of those processor ranks.
c       spp%sendtot(1:spp%nsnd) is the number of potentially nonzero
c         matrix elements going to each receving processor.
c       spp%sendst(1:spp%nsnd) is the first index in the spp%jentry list
c         that is sent to each receiving processor.
c       spp%sendst(1:spp%nsnd) is the first row global index that is
c         sent to each receiving processor.
c-----------------------------------------------------------------------
      SUBROUTINE pardir_setup(spp,nb,nnz2,irw,jentry)
      USE mpi_nim
      USE pardata
      USE factor_type_mod, ONLY:sparsity_pattern

      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      INTEGER(i4), INTENT(IN) :: nb,nnz2
      INTEGER(i4), DIMENSION(nnz2), INTENT(IN) :: irw,jentry

      INTEGER(i4) :: ib,jb,ig,jg,ierror,ipr,jpr,ncol,irow,ii,jj,maxrcv
      INTEGER(i4) :: ircv,isnd,nb_total,id
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: tmpncol,tmpst,tmprowst
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: tmpn,tmp2

      IF (nprocs_layer<=1) THEN
        spp%nsnd=0; spp%nrcv=0
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     use the list of indices in the compressed row storage for
c     contributions from this block to determine where to send data.
c
c     loop over other processor's blocks and look for matching rows.
c-----------------------------------------------------------------------
      ALLOCATE(tmpncol(0:nprocs_layer-1),tmpst(0:nprocs_layer-1),
     $         tmprowst(0:nprocs_layer-1))

      tmpncol=0
      tmpst=nnz2+1
      tmprowst=spp%nrow

      nb_total=SIZE(global2local)
      DO ib=1,nb_total
        ipr=block2proc(ib)-ilayer*nprocs_layer
        id=spp%gblock_order(ib)
        DO jj=1,nnz2
          irow=irw(jj)
          IF (irow>=spp%irowst_block(id-1)
     $        .AND.irow<spp%irowst_block(id)) THEN
            tmpncol(ipr)=tmpncol(ipr)+1
            tmpst(ipr)=MIN(tmpst(ipr),jj)
            tmprowst(ipr)=MIN(tmprowst(ipr),irow)
          ENDIF
        ENDDO
      ENDDO

      spp%nsnd=0
      DO ipr=0,nprocs_layer-1
        IF (tmpncol(ipr)>0.AND.ipr/=node_layer) THEN
          spp%nsnd=spp%nsnd+1
        ENDIF
      ENDDO

      IF (spp%nsnd>0) THEN
        ALLOCATE(spp%sendlist(spp%nsnd))
        ALLOCATE(spp%sendtot(spp%nsnd))
        ALLOCATE(spp%sendst(spp%nsnd))
        ALLOCATE(spp%sendrowst(spp%nsnd))
        ALLOCATE(spp%send_req(spp%nsnd))
        isnd=1
        DO ipr=0,nprocs_layer-1
          IF (tmpncol(ipr)>0.AND.ipr/=node_layer) THEN
            spp%sendlist(isnd)=ipr
            spp%sendtot(isnd)=tmpncol(ipr)
            spp%sendst(isnd)=tmpst(ipr)
            spp%sendrowst(isnd)=tmprowst(ipr)
            isnd=isnd+1
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     inform other processors to expect communication.
c-----------------------------------------------------------------------
      spp%nrcv=0
      ALLOCATE(tmpn(0:nprocs_layer-1,0:nprocs_layer-1))
      ALLOCATE(tmp2(0:nprocs_layer-1,0:nprocs_layer-1))
      tmpn=0
      DO ii=1,spp%nsnd
        tmpn(node_layer,spp%sendlist(ii))=spp%sendtot(ii)
      ENDDO
      CALL mpi_allreduce(tmpn(0,0),tmp2(0,0),nprocs_layer**2,
     $                   mpi_nim_int,mpi_sum,comm_layer,ierror)
      tmpn=tmp2
      DO ii=0,nprocs_layer-1
        IF (tmpn(ii,node_layer)>0) spp%nrcv=spp%nrcv+1
      ENDDO

      IF (spp%nrcv>0) THEN
        ALLOCATE(spp%recvlist(spp%nrcv))
        ALLOCATE(spp%recvtot(spp%nrcv))
        ALLOCATE(spp%recv_req(spp%nrcv))
        ircv=0
        DO ii=0,nprocs_layer-1
          IF (tmpn(ii,node_layer)>0) THEN
            ircv=ircv+1
            spp%recvlist(ircv)=ii
            spp%recvtot(ircv)=tmpn(ii,node_layer)
          ENDIF
        ENDDO
        ALLOCATE(spp%recvjentry(MAXVAL(spp%recvtot),spp%nrcv))
      ENDIF

      DEALLOCATE(tmpn,tmp2)
c-----------------------------------------------------------------------
c     post receives for the list of entries coming from different
c     processors.
c-----------------------------------------------------------------------
      DO ircv=1,spp%nrcv
        CALL mpi_irecv(spp%recvjentry(1,ircv),spp%recvtot(ircv),
     $                 mpi_nim_int,spp%recvlist(ircv),0_i4,
     $                 comm_layer,spp%recv_req(ircv),ierror)
      ENDDO
c-----------------------------------------------------------------------
c     send the list of compressed-row indices that will be contributed
c     to other processors.
c-----------------------------------------------------------------------
      DO isnd=1,spp%nsnd
        CALL mpi_send(jentry(spp%sendst(isnd)),spp%sendtot(isnd),
     $                mpi_nim_int,spp%sendlist(isnd),0_i4,
     $                comm_layer,ierror)
      ENDDO
      IF (spp%nrcv>0) THEN
        ALLOCATE(tmpn(mpi_status_size,spp%nrcv))
        CALL mpi_waitall(spp%nrcv,spp%recv_req,tmpn,ierror)
        DEALLOCATE(tmpn)
      ENDIF

      DEALLOCATE(tmpncol,tmpst)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE pardir_setup

c-----------------------------------------------------------------------
c     subprogram 14. iter_dirfac_dealloc.
c     deallocates storage for preconditioner matrix factors.
c-----------------------------------------------------------------------
      SUBROUTINE iter_dirfac_dealloc(spp)
      USE pardata
      USE factor_type_mod, ONLY:sparsity_pattern
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      INTEGER(i4) :: isnd

      DEALLOCATE(spp%irowst_block)
      DEALLOCATE(spp%gblock_order)
      DEALLOCATE(spp%j_acc)
      DEALLOCATE(spp%start_acc)
      IF (ASSOCIATED(spp%j_acc_save)) DEALLOCATE(spp%j_acc_save)
      IF (ASSOCIATED(spp%start_acc_save)) DEALLOCATE(spp%start_acc_save)
      DEALLOCATE(spp%jentry)
      IF (spp%matrix_distributed) THEN
        DEALLOCATE(spp%start_loc)
      ELSE
        DEALLOCATE(spp%algcount,spp%algdispl)
        DEALLOCATE(spp%algcntr,spp%algdsplr)
      ENDIF
      IF (spp%nrcv>0) THEN
        DEALLOCATE(spp%recvlist,spp%recvtot,
     $             spp%recvjentry,spp%recv_req)
      ENDIF
      IF (spp%nsnd>0) THEN
        DEALLOCATE(spp%sendlist,spp%sendtot,spp%sendst,spp%send_req)
        IF (ASSOCIATED(spp%sendrowst)) DEALLOCATE(spp%sendrowst)
        IF (ASSOCIATED(spp%sendstacc)) THEN
          DO isnd=1,spp%nsnd
            DEALLOCATE(spp%sendstacc(isnd)%data)
          ENDDO
          DEALLOCATE(spp%sendstacc)
        ENDIF
      ENDIF
      IF (spp%solve_nopts>0) DEALLOCATE(spp%iopts,spp%dopts)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_dirfac_dealloc

c-----------------------------------------------------------------------
c     subprogram 15. set_precon_opts
c     set preconditioner options.  in this version, they are hardcoded
c     here.
c-----------------------------------------------------------------------
      SUBROUTINE set_precon_opts(solve_flag,solve_nopts,iopts,dopts)

      CHARACTER(*), INTENT(IN) :: solve_flag
      INTEGER(i4), INTENT(OUT) :: solve_nopts
      INTEGER(i4), DIMENSION(:), POINTER :: iopts
      REAL(r8), DIMENSION(:), POINTER :: dopts

      INTEGER(i4) :: ii

      SELECT CASE(solve_flag)
      CASE ("slu_dist", "slu_dstm", "slu_dsta")
        solve_nopts=5
        ALLOCATE(iopts(solve_nopts),dopts(solve_nopts))
        iopts(1)=0  !  no reporting
        iopts(2)=0  !  equilibrate option is false
        iopts(3)=1  !  diagonal scaling, see c_fortran_pzloc.c
        iopts(4)=2  !  column permutation, see c_fortran_pzloc.c
        iopts(5)=0  !  print stats is false
        dopts(1:solve_nopts)=0     !  not used
      CASE DEFAULT
        solve_nopts=0
        NULLIFY(iopts,dopts)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE set_precon_opts
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_dir_fac
