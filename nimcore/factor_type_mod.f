c-----------------------------------------------------------------------
c     file factor_type_mod.f
c     contains a module that defines structures for saving matrix
c     factors.  (previously part of matrix_type_mod).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. factor_ptr_copy_rmat_info
c     2. factor_dealloc_rmat_info
c     3. factor_ptr_copy_cmat_info
c     4. factor_dealloc_cmat_info
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE factor_type_mod
      USE local
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     type for holding the distributed col/row pointer for sent data.
c-----------------------------------------------------------------------
      TYPE :: send_start_acc
        INTEGER(i4), POINTER :: data(:,:)
      END TYPE send_start_acc
c-----------------------------------------------------------------------
c     type for holding the sparsity pattern for direct solves.
c-----------------------------------------------------------------------
      TYPE :: sparsity_pattern
        !> ordering of global blocks in the sparsity pattern
        INTEGER(i4), POINTER :: gblock_order(:)
        !> running total of unique nodes by block (1:nbl_total)
        !> where nrow=irowst_block(nbl_total) and the starting row for
        !> any block of global number id is given by
        !> irowst_block(spp%gblock_order(id)-1)
        INTEGER(i4), POINTER :: irowst_block(:)
        !> global (local if sparsity distributed) row/col index
        INTEGER(i4), POINTER :: j_acc(:)
        !> global (local if sparsity distributed) matrix col/row ptr
        INTEGER(i4), POINTER :: start_acc(:)
        !> distributed (CSC/CSR) matrix: 0-index col/row ptr
        INTEGER(i4), POINTER :: start_loc(:)
        !> save location for j_acc if overwritten by the external solver
        INTEGER(i4), POINTER :: j_acc_save(:)
        !> save location if overwritten by the external solver
        INTEGER(i4), POINTER :: start_acc_save(:)
        INTEGER(i4) :: nrow   !< number of rows (global mat)
        INTEGER(i4) :: nnz    !< number of non-zeros (global mat)
        INTEGER(i4) :: nbw    !< maximum width between irow and icol
        INTEGER(i4) :: mloc   !< number of local rows
        INTEGER(i4) :: fstrow !< first local row
        INTEGER(i4) :: lstrow !< last local row
        INTEGER(i4) :: nnzloc !< local number of non-zeros
        INTEGER(i4) :: indst      !< local starting index of val/j_acc
        INTEGER(i4) :: indend     !< local ending index of val/j_acc
        INTEGER(i4) :: nsnd       !< node number of sends
        INTEGER(i4) :: nrcv       !< node number of receives
        !> list of nodes that send to this node (size nrcv)
        INTEGER(i4), POINTER :: recvlist(:)
        !> list of size of receives (size nrcv)
        INTEGER(i4), POINTER :: recvtot(:)
        !> location in acc of received data (size max(recvtot),nrcv)
        INTEGER(i4), POINTER :: recvjentry(:,:)
        !> list of nodes that this node is sends to (size nsnd)
        INTEGER(i4), POINTER :: sendlist(:)
        !> list of size of sends (size nsnd)
        INTEGER(i4), POINTER :: sendtot(:)
        !> location in jentry() where send starts (size nsnd)
        INTEGER(i4), POINTER :: sendst(:)
        !> row index where send starts (size nsnd)
        INTEGER(i4), POINTER :: sendrowst(:)
        !> col/row ptr for sends
        TYPE(send_start_acc), POINTER :: sendstacc(:)
        !> location in acc of sent data (col ind if dist sparsity)
        INTEGER(i4), POINTER :: jentry(:)
        !> asynchronous receive request handle
        INTEGER(i4), POINTER :: recv_req(:)
        !> asynchronous send request handle
        INTEGER(i4), POINTER :: send_req(:)
        !> for global matrix gather communication
        INTEGER(i4), DIMENSION(:), POINTER :: algdispl,algcount
        INTEGER(i4), DIMENSION(:), POINTER :: algdsplr,algcntr
        !> true if matrix has been factored previously
        LOGICAL :: acc_lustored
        !> true if the external solver distributed interface is used
        LOGICAL :: matrix_distributed
        !> true if the sparsity pattern is distributed
        LOGICAL :: sparsity_distributed
        !> true if a direct solver is being used
        LOGICAL :: direct_solver
        !> SuperLU handle for bridge routine storage
        INTEGER(i8) :: acc_handle
        !> integer options
        INTEGER(i4) :: solve_nopts
        INTEGER(i4), DIMENSION(:), POINTER :: iopts
        !> double precision options
        REAL(r8), DIMENSION(:), POINTER :: dopts
      END TYPE sparsity_pattern
c-----------------------------------------------------------------------
c     copies of matrix specifications that are packaged with the
c     factor structures.
c-----------------------------------------------------------------------
      TYPE :: rbl_mat_info
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
      END TYPE rbl_mat_info
c-----------------------------------------------------------------------
c     types used to define arrays for matrix factors and working
c     pointers.
c-----------------------------------------------------------------------
      TYPE :: matrix_factor_type
        TYPE(bl_fac_type),      DIMENSION(:), POINTER :: bl_fac
        TYPE(bl_sparsity_type), DIMENSION(:), POINTER :: bl_spp
        TYPE(rbl_mat_info),     DIMENSION(:), POINTER :: mat_info
        TYPE(sparsity_pattern) :: spp
        INTEGER(i4) :: nqty,ilu_fill
        !> factors for precon=='lapack'
        REAL(r8), DIMENSION(:,:), POINTER :: a11
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE matrix_factor_type
      TYPE :: complex_factor_type
        TYPE(bl_comp_fac_type), DIMENSION(:), POINTER :: bl_fac
        TYPE(bl_sparsity_type), DIMENSION(:), POINTER :: bl_spp
        TYPE(rbl_mat_info),     DIMENSION(:), POINTER :: mat_info
        TYPE(sparsity_pattern) :: spp
        INTEGER(i4) :: nqty,ilu_fill
        !> factors for precon=='lapack'
        COMPLEX(r8), DIMENSION(:,:), POINTER :: a11
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE complex_factor_type
c-----------------------------------------------------------------------
c     block by block matrix factor arrays.  no one preconditioner uses
c     all of these arrays, and only those needed are allocated.
c-----------------------------------------------------------------------
      TYPE :: bl_fac_type
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: lt_dg0,ut_dg0
        REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: lt_dg1,ut_dg1,
     $            xelim,yelim,xlc,xcl,ylc,ycl
        REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: lt_cen,lt_off,
     $            lt_per,ut_cen,ut_off,ut_per,adix,adiy
        TYPE(arr_4d_type), DIMENSION(:), POINTER :: pmat
        INTEGER(i4) :: nv,mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(ikind2), DIMENSION(:), ALLOCATABLE :: ixm,iym
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_fac_type
      TYPE :: bl_comp_fac_type
        COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: lt_dg0,ut_dg0
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: lt_dg1,ut_dg1,
     $               xelim,yelim,xlc,xcl,ylc,ycl
        COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: lt_cen,lt_off,
     $            lt_per,ut_cen,ut_off,ut_per,adix,adiy
        TYPE(comp_arr_4d_type), DIMENSION(:), POINTER :: pmat
        INTEGER(i4) :: nv,mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(ikind2), DIMENSION(:), ALLOCATABLE :: ixm,iym
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_comp_fac_type
c-----------------------------------------------------------------------
c     types used for holding inverses of diagonal point-block.
c-----------------------------------------------------------------------
      TYPE :: arr_4d_type
        REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: arr
      END TYPE arr_4d_type
      TYPE :: comp_arr_4d_type
        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: arr
      END TYPE comp_arr_4d_type
c-----------------------------------------------------------------------
c     type for block specific information for direct solves.
c-----------------------------------------------------------------------
      TYPE :: bl_sparsity_type
        TYPE(row_type), DIMENSION(:), POINTER :: row_ind
        LOGICAL :: perblock,degenerate
        !> memory id for memory logging.
        INTEGER(i4), POINTER :: mem_id
      END TYPE bl_sparsity_type
c-----------------------------------------------------------------------
c     type for assigning row numbers and organizing direct solves.
c-----------------------------------------------------------------------
      TYPE :: row_type
        INTEGER(i4), DIMENSION(:,:,:), POINTER :: rarr
      END TYPE row_type

c-----------------------------------------------------------------------
c     interfaces for the allocation and deallocation routines.
c-----------------------------------------------------------------------
      INTERFACE factor_ptr_copy_mat_info
        MODULE PROCEDURE factor_ptr_copy_rmat_info,
     $                   factor_ptr_copy_cmat_info
      END INTERFACE

      INTERFACE factor_dealloc_mat_info
        MODULE PROCEDURE factor_dealloc_rmat_info,
     $                   factor_dealloc_cmat_info
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. factor_ptr_copy_rmat_info
c     log the memory of the matrix factor structure.
c-----------------------------------------------------------------------
      SUBROUTINE factor_ptr_copy_rmat_info(fac,rmat)
      USE matrix_type_mod

      TYPE(matrix_factor_type), INTENT(OUT) :: fac
      TYPE(rbl_mat_type), INTENT(IN) :: rmat(:)
      INTEGER(i4) :: ibl,nbl

      nbl=SIZE(rmat)
      ALLOCATE(fac%mat_info(nbl))
      DO ibl=1,nbl
        fac%mat_info(ibl)%nbasis_el=rmat(ibl)%nbasis_el
        fac%mat_info(ibl)%nbtype=rmat(ibl)%nbtype
        fac%mat_info(ibl)%mx=rmat(ibl)%mx
        fac%mat_info(ibl)%my=rmat(ibl)%my
        fac%mat_info(ibl)%ix0=>rmat(ibl)%ix0
        fac%mat_info(ibl)%iy0=>rmat(ibl)%iy0
        fac%mat_info(ibl)%nb_type=>rmat(ibl)%nb_type
        fac%mat_info(ibl)%nq_type=>rmat(ibl)%nq_type
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_ptr_copy_rmat_info
c-----------------------------------------------------------------------
c     subprogram 2. factor_dealloc_rmat_info
c     log the memory of the matrix factor structure.
c-----------------------------------------------------------------------
      SUBROUTINE factor_dealloc_rmat_info(fac)
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: ibl

      DEALLOCATE(fac%mat_info)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_dealloc_rmat_info
c-----------------------------------------------------------------------
c     subprogram 3. factor_ptr_copy_cmat_info
c     log the memory of the matrix factor structure.
c-----------------------------------------------------------------------
      SUBROUTINE factor_ptr_copy_cmat_info(fac,cmat)
      USE matrix_type_mod

      TYPE(complex_factor_type), INTENT(OUT) :: fac
      TYPE(rbl_comp_mat_type), INTENT(IN) :: cmat(:)
      INTEGER(i4) :: ibl,nbl

      nbl=SIZE(cmat)
      ALLOCATE(fac%mat_info(nbl))
      DO ibl=1,nbl
        fac%mat_info(ibl)%nbasis_el=cmat(ibl)%nbasis_el
        fac%mat_info(ibl)%nbtype=cmat(ibl)%nbtype
        fac%mat_info(ibl)%mx=cmat(ibl)%mx
        fac%mat_info(ibl)%my=cmat(ibl)%my
        fac%mat_info(ibl)%ix0=>cmat(ibl)%ix0
        fac%mat_info(ibl)%iy0=>cmat(ibl)%iy0
        fac%mat_info(ibl)%nb_type=>cmat(ibl)%nb_type
        fac%mat_info(ibl)%nq_type=>cmat(ibl)%nq_type
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_ptr_copy_cmat_info
c-----------------------------------------------------------------------
c     subprogram 4. factor_dealloc_cmat_info
c     log the memory of the matrix factor structure.
c-----------------------------------------------------------------------
      SUBROUTINE factor_dealloc_cmat_info(fac)
      TYPE(complex_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4) :: ibl

      DEALLOCATE(fac%mat_info)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE factor_dealloc_cmat_info

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE factor_type_mod
