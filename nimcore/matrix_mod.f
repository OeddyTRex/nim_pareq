c-----------------------------------------------------------------------
c     file matrix_mod.f
c     module containing routines that performs standard operations on
c     global matrices.
c
c     this module has been broken into a number of separate modules
c     to facilitate compilation on some computers.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module matvec_real_mod
c     1. matvec_real.
c     2. matvecgg_real_rbl.
c     3. matvecgh_real_rbl.
c     4. matvecgv_real_rbl.
c     5. matvecgi_real_rbl.
c     6. matvechg_real_rbl.
c     7. matvechh_real_rbl.
c     8. matvechv_real_rbl.
c     9. matvechi_real_rbl.
c     10. matvecvg_real_rbl.
c     11. matvecvh_real_rbl.
c     12. matvecvv_real_rbl.
c     13. matvecvi_real_rbl.
c     14. matvecig_real_rbl.
c     15. matvecih_real_rbl.
c     16. matveciv_real_rbl.
c     17. matvecii_real_rbl.
c     18. matvec_real_tbl.
c     module matvec_comp_mod
c     19. matvec_comp.
c     20. matvec_2D_comp.
c     21. matvecgg_comp_rbl.
c     22. matvecgh_comp_rbl.
c     23. matvecgv_comp_rbl.
c     24. matvecgi_comp_rbl.
c     25. matvechg_comp_rbl.
c     26. matvechh_comp_rbl.
c     27. matvechv_comp_rbl.
c     28. matvechi_comp_rbl.
c     29. matvecvg_comp_rbl.
c     30. matvecvh_comp_rbl.
c     31. matvecvv_comp_rbl.
c     32. matvecvi_comp_rbl.
c     33. matvecig_comp_rbl.
c     34. matvecih_comp_rbl.
c     35. matveciv_comp_rbl.
c     36. matvecii_comp_rbl.
c     37. matvec_comp_tbl.
c     module matelim_mod
c     38. matelim_real_inv_int.
c     39. matelim_real_presolve.
c     40. matelim_real_postsolve.
c     41. matelim_comp_inv_int.
c     42. matelim_comp_presolve.
c     43. matelim_comp_postsolve.
c     module matrix_mod
c     44. matrix_lump_real_rmat
c     45. matrix_lump_real_tmat
c     46. matrix_degen_collect_real
c     47. matrix_lump_comp_rmat
c     48. matrix_lump_comp_tmat
c     49. matrix_degen_collect_comp
c-----------------------------------------------------------------------
c     matvec_real_mod contains routines for matrix multiplication of
c     real data.
c-----------------------------------------------------------------------
      MODULE matvec_real_mod
      USE local
      USE matrix_type_mod
      USE vector_type_mod
      USE edge
      USE seam_storage_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. matvec_real.
c     multiply the real mat matrix by the real operand vector,
c     and return the product.  each matrix array contains all
c     connections between two sets of basis types, so operand and
c     product arrays have basis index lumped with quantity index.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvec_real(mat,operand,product,nq,do_net)
      USE time

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(vector_type), DIMENSION(:), INTENT(IN) :: operand
      TYPE(vector_type), DIMENSION(:), INTENT(OUT) :: product
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN), OPTIONAL :: do_net

      INTEGER(i4) :: ibl,iv,mx,my,mv,in,nrb=0,ntb=0,ix,iy,pd,
     $               ityp,jtyp,ibasis,nbtype,typ_max,nqp,nqo
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: op_ptr,pr_ptr
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS :: mat_ptr
      REAL(r8) :: timematv_st,timematv_en
      LOGICAL :: new_product,netw
c-----------------------------------------------------------------------
c     if the optional do_net input is provided and set to false, skip
c     the networking operation.
c-----------------------------------------------------------------------
      netw=.true.
      IF (PRESENT(do_net)) THEN
        IF (.NOT.do_net) netw=.false.
      ENDIF
c-----------------------------------------------------------------------
c     compute the number of blocks.
c-----------------------------------------------------------------------
      CALL timer(timematv_st)
      nrb=SIZE(mat%rbl_mat)
      IF (nrb>0) THEN
        nbtype=mat%rbl_mat(1)%nbtype
        IF (nbtype>1) THEN
          pd=SUM(mat%rbl_mat(1)%nb_type(1:2))
        ELSE
          pd=1
        ENDIF
      ELSE
        nbtype=1
        pd=1
      ENDIF
      ntb=SIZE(mat%tbl_mat)
c-----------------------------------------------------------------------
c     if interior nodes are eliminated, basis loops only cover grid
c     and element side centered nodes.
c-----------------------------------------------------------------------
      IF (mat%eliminated) THEN
        typ_max=MIN(nbtype,3_i4)
      ELSE
        typ_max=nbtype
      ENDIF
c-----------------------------------------------------------------------
c     rblocks: call the structured-data matrix operation for each
c     basis type to basis type pair, and load seam array.
c     every basis-type pair now uses a separate multiplication routine
c     for optimization.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
        my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
        DO ityp=1,typ_max
          new_product=.true.
          nqp=mat%rbl_mat(ibl)%nq_type(ityp)
          SELECT CASE(ityp)
          CASE(1)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgg_real_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechg_real_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvg_real_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecig_real_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(2)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgh_real_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechh_real_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvh_real_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecih_real_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(3)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgv_real_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechv_real_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvv_real_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matveciv_real_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(4)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgi_real_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechi_real_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvi_real_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecii_real_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          END SELECT
        ENDDO
        IF (netw) CALL edge_load_arr(product(ibl),nq,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     tblocks: call the unstructured-data matrix operation and load
c     seam array.
c-----------------------------------------------------------------------
      DO ibl=nrb+1,nrb+ntb
        mv=SIZE(mat%tbl_mat(ibl)%lmat)-1
        pr_ptr=>product(ibl)%arr
        pr_ptr=0
        op_ptr=>operand(ibl)%arr
        CALL matvec_real_tbl(pr_ptr,mat%tbl_mat(ibl)%lmat,op_ptr,mv,nq)
        IF (netw) CALL edge_load_arr(product(ibl),nq,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     communicate across block boundaries and unload seam.
c-----------------------------------------------------------------------
      IF (netw) THEN
        CALL edge_network(nq,0_i4,pd-1_i4,.false.)
        DO ibl=1,nrb+ntb
          CALL edge_unload_arr(product(ibl),nq,pd-1_i4,seam(ibl))
        ENDDO
      ENDIF
      CALL timer(timematv_en)
      time_matvec=time_matvec+timematv_en-timematv_st
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvec_real
c-----------------------------------------------------------------------
c     subprogram 2. matvecgg_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgg_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,-1:,:,0:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,0:,0:,iq,0,0)*vector(:,:1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,:0,0:,iq,mx,0)*vector(:,mx-1:,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,0:,:0,iq,0,my)*vector(:,:1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,:0,:0,iq,mx,my)*vector(:,mx-1:,my-1:my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,0:,0:,iq,0,0)*vector(:,:1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,:0,0:,iq,mx,0)*vector(:,mx-1:,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,0:,:0,iq,0,my)*vector(:,:1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,:0,:0,iq,mx,my)*vector(:,mx-1:,my-1:my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgg_real_rbl
c-----------------------------------------------------------------------
c     subprogram 3. matvecgh_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid
c     vertex to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgh_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,-1:,:,1:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix,my-1:my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix,my-1:my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgh_real_rbl
c-----------------------------------------------------------------------
c     subprogram 4. matvecgv_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgv_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,-1:,:,0:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgv_real_rbl
c-----------------------------------------------------------------------
c     subprogram 5. matvecgi_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgi_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,-1:,:,1:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgi_real_rbl
c-----------------------------------------------------------------------
c     subprogram 6. matvechg_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechg_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,-1:,:,0:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,1,0:,iq,0,0)*vector(:,1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,0,0:,iq,mx,0)*vector(:,mx,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,1,:0,iq,0,my)*vector(:,1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,0,:0,iq,mx,my)*vector(:,mx,my-1:my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,1,0:,iq,0,0)*vector(:,1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,0,0:,iq,mx,0)*vector(:,mx,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,1,:0,iq,0,my)*vector(:,1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,0,:0,iq,mx,my)*vector(:,mx,my-1:my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechg_real_rbl
c-----------------------------------------------------------------------
c     subprogram 7. matvechh_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechh_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,-1:,:,1:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix,my-1:my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix,my-1:my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechh_real_rbl
c-----------------------------------------------------------------------
c     subprogram 8. matvechv_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechv_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,-1:,:,0:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechv_real_rbl
c-----------------------------------------------------------------------
c     subprogram 9. matvechi_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechi_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,-1:,:,1:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechi_real_rbl
c-----------------------------------------------------------------------
c     subprogram 10. matvecvg_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvg_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,0:,:,0:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,0:,1,iq,0,0)*vector(:,:1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,:0,1,iq,mx,0)*vector(:,mx-1:,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,0:,0,iq,0,my)*vector(:,:1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,:0,0,iq,mx,my)*vector(:,mx-1:,my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,0:,1,iq,0,0)*vector(:,:1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,:0,1,iq,mx,0)*vector(:,mx-1:,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,0:,0,iq,0,my)*vector(:,:1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,:0,0,iq,mx,my)*vector(:,mx-1:,my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvg_real_rbl
c-----------------------------------------------------------------------
c     subprogram 11. matvecvh_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvh_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,0:,:,1:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix,my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix,my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvh_real_rbl
c-----------------------------------------------------------------------
c     subprogram 12. matvecvv_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvv_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,0:,:,0:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,0,iq,0,iy)*vector(:,:1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,0,iq,mx,iy)*vector(:,mx-1:,iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,0,iq,0,iy)*vector(:,:1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,0,iq,mx,iy)*vector(:,mx-1:,iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvv_real_rbl
c-----------------------------------------------------------------------
c     subprogram 13. matvecvi_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvi_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,-1:,0:,:,1:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix,iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix,iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvi_real_rbl
c-----------------------------------------------------------------------
c     subprogram 14. matvecig_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecig_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,0:,:,0:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,1,1,iq,0,0)*vector(:,1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,0,1,iq,mx,0)*vector(:,mx,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,1,0,iq,0,my)*vector(:,1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,0,0,iq,mx,my)*vector(:,mx,my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,1,1,iq,0,0)*vector(:,1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,0,1,iq,mx,0)*vector(:,mx,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,1,0,iq,0,my)*vector(:,1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,0,0,iq,mx,my)*vector(:,mx,my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecig_real_rbl
c-----------------------------------------------------------------------
c     subprogram 15. matvecih_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecih_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,0:,:,1:,0:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,0,1,iq,ix,0)*vector(:,ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,0,0,iq,ix,my)*vector(:,ix,my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,0,1,iq,ix,0)*vector(:,ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,0,0,iq,ix,my)*vector(:,ix,my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecih_real_rbl
c-----------------------------------------------------------------------
c     subprogram 16. matveciv_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matveciv_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,0:,:,0:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,0,iq,0,iy)*vector(:,1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,0,iq,mx,iy)*vector(:,mx,iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,0,iq,0,iy)*vector(:,1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,0,iq,mx,iy)*vector(:,mx,iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matveciv_real_rbl
c-----------------------------------------------------------------------
c     subprogram 17. matvecii_real_rbl.
c     perform a real matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to interior only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecii_real_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      REAL(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      REAL(r8), DIMENSION(:,0:,0:,:,1:,1:), INTENT(IN) ::
     $          matrix
      REAL(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,0,0,iq,ix,iy)*vector(:,ix,iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,0,0,iq,ix,iy)*vector(:,ix,iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecii_real_rbl
c-----------------------------------------------------------------------
c     subprogram 18. matvec_real_tbl.
c     perform a matrix/vector multiplication for a tblock and return
c     the product in the array result.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvec_real_tbl(result,matrix,vector,mv,nqty)

      INTEGER(i4), INTENT(IN) :: mv,nqty
      TYPE(matrix_element_type3), DIMENSION(:), POINTER :: matrix
      REAL(r8), DIMENSION(:,0:,0:), INTENT(IN) :: vector
      REAL(r8), DIMENSION(:,0:,0:), INTENT(OUT) :: result

      INTEGER(i4) :: iqty,jqty,ivert,inbr,mvert
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: matel
c-----------------------------------------------------------------------
c     find the matrix/vector product.
c-----------------------------------------------------------------------
      DO ivert=0,mv
        matel=>matrix(ivert)%element
        DO inbr=0,SIZE(matrix(ivert)%from_vert)-1
          DO iqty=1,nqty
            result(iqty,ivert,0)=result(iqty,ivert,0)
     $        +SUM(matel(:,iqty,inbr)
     $            *vector(1:nqty,matrix(ivert)%from_vert(inbr),0))
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvec_real_tbl
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matvec_real_mod


c-----------------------------------------------------------------------
c     matvec_comp_mod contains routines for matrix multiplication of
c     complex data.
c-----------------------------------------------------------------------
      MODULE matvec_comp_mod
      USE local
      USE matrix_type_mod
      USE vector_type_mod
      USE edge
      USE seam_storage_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 19. matvec_comp.
c     multiply the complex mat matrix by the complex operand vector,
c     and return the product.  this version uses cvectors_types, and
c     the Fourier component index must be specified.
c
c     each matrix array contains all connections between two sets of
c     basis types, so operand and product arrays have basis index
c     lumped with quantity index.
c
c     it is now assumed that the operand
c     and product have vector component dimensioned nq.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvec_comp(mat,operand,product,nq,ifo,ifp,do_net)
      USE time

      TYPE(complex_matrix_type), INTENT(IN) :: mat
      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: operand
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: product
      INTEGER(i4), INTENT(IN) :: nq,ifo,ifp
      LOGICAL, INTENT(IN), OPTIONAL :: do_net

      INTEGER(i4) :: ibl,iv,mx,my,mv,in,nrb=0,ntb=0,ix,iy,pd,
     $               ityp,jtyp,ibasis,nbtype,typ_max,nqp,nqo
      COMPLEX(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $             op_ptr,pr_ptr
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS ::
     $             mat_ptr
      REAL(r8) :: timematv_st,timematv_en
      LOGICAL :: new_product,netw
c-----------------------------------------------------------------------
c     if the optional do_net input is provided and set to false, skip
c     the networking operation.
c-----------------------------------------------------------------------
      netw=.true.
      IF (PRESENT(do_net)) THEN
        IF (.NOT.do_net) netw=.false.
      ENDIF
c-----------------------------------------------------------------------
c     compute the number of blocks.
c-----------------------------------------------------------------------
      CALL timer(timematv_st)
      nrb=SIZE(mat%rbl_mat)
      IF (nrb>0) THEN
        nbtype=mat%rbl_mat(1)%nbtype
        IF (nbtype>1) THEN
          pd=SUM(mat%rbl_mat(1)%nb_type(1:2))
        ELSE
          pd=1
        ENDIF
      ELSE
        nbtype=1
        pd=1
      ENDIF
      ntb=SIZE(mat%tbl_mat)
c-----------------------------------------------------------------------
c     if interior nodes are eliminated, basis loops only cover grid
c     and element side centered nodes.
c-----------------------------------------------------------------------
      IF (mat%eliminated) THEN
        typ_max=MIN(nbtype,3_i4)
      ELSE
        typ_max=nbtype
      ENDIF
c-----------------------------------------------------------------------
c     rblocks: call the structured-data matrix operation for each
c     basis type to basis type pair, and load seam array.
c     every basis-type pair now uses a separate multiplication routine
c     for optimization.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
        my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
        DO ityp=1,typ_max
          new_product=.true.
          nqp=mat%rbl_mat(ibl)%nq_type(ityp)
          SELECT CASE(ityp)
          CASE(1)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgg_comp_rbl(product(ibl)%arr(:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arr(:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechg_comp_rbl(product(ibl)%arr(:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrh(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvg_comp_rbl(product(ibl)%arr(:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrv(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecig_comp_rbl(product(ibl)%arr(:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arri(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(2)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgh_comp_rbl(product(ibl)%arrh(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arr(:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechh_comp_rbl(product(ibl)%arrh(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrh(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvh_comp_rbl(product(ibl)%arrh(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrv(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecih_comp_rbl(product(ibl)%arrh(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arri(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(3)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgv_comp_rbl(product(ibl)%arrv(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arr(:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechv_comp_rbl(product(ibl)%arrv(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrh(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvv_comp_rbl(product(ibl)%arrv(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrv(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matveciv_comp_rbl(product(ibl)%arrv(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arri(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(4)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgi_comp_rbl(product(ibl)%arri(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arr(:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechi_comp_rbl(product(ibl)%arri(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrh(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvi_comp_rbl(product(ibl)%arri(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arrv(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecii_comp_rbl(product(ibl)%arri(:,:,:,:,ifp),
     $                                 mat_ptr,
     $                                 operand(ibl)%arri(:,:,:,:,ifo),
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          END SELECT
        ENDDO
        IF (netw)
     $    CALL edge_load_carr(product(ibl),nq,ifp,ifp,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     tblocks: call the unstructured-data matrix operation and load
c     seam array.
c-----------------------------------------------------------------------
      DO ibl=nrb+1,nrb+ntb
        mv=SIZE(mat%tbl_mat(ibl)%lmat)-1
        pr_ptr=>product(ibl)%arr(:,:,:,ifp)
        pr_ptr=0
        op_ptr=>operand(ibl)%arr(:,:,:,ifo)
        CALL matvec_comp_tbl(pr_ptr,mat%tbl_mat(ibl)%lmat,op_ptr,mv,nq)
        IF (netw)
     $    CALL edge_load_carr(product(ibl),nq,ifp,ifp,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     communicate across block boundaries and unload seam.
c-----------------------------------------------------------------------
      IF (netw) THEN
        CALL edge_network(nq,1_i4,pd-1_i4,.false.)
        DO ibl=1,nrb+ntb
          CALL edge_unload_carr(product(ibl),nq,ifp,ifp,pd-1_i4,
     $                          seam(ibl))
        ENDDO
      ENDIF
      CALL timer(timematv_en)
      time_matvec=time_matvec+timematv_en-timematv_st
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvec_comp
c-----------------------------------------------------------------------
c     subprogram 20. matvec_2D_comp.
c     multiply the complex mat matrix by the complex operand vector,
c     and return the product.  this version uses cvectors_2D_types, 
c     which do not have a Fourier component index.
c
c     each matrix array contains all connections between two sets of
c     basis types, so operand and product arrays have basis index
c     lumped with quantity index.

c     it is now assumed that the operand
c     and product have vector component dimensioned nq.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvec_2D_comp(mat,operand,product,nq,do_net)
      USE time

      TYPE(complex_matrix_type), INTENT(IN) :: mat
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(IN) :: operand
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(OUT) :: product
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN), OPTIONAL :: do_net

      INTEGER(i4) :: ibl,iv,mx,my,mv,in,nrb=0,ntb=0,ix,iy,pd,
     $               ityp,jtyp,ibasis,nbtype,typ_max,nqp,nqo
      COMPLEX(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS ::
     $             op_ptr,pr_ptr
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS ::
     $             mat_ptr
      LOGICAL :: new_product,netw
      REAL(r8) :: timematv_st,timematv_en
c-----------------------------------------------------------------------
c     if the optional do_net input is provided and set to false, skip
c     the networking operation.
c-----------------------------------------------------------------------
      netw=.true.
      IF (PRESENT(do_net)) THEN
        IF (.NOT.do_net) netw=.false.
      ENDIF
c-----------------------------------------------------------------------
c     compute the number of blocks.
c-----------------------------------------------------------------------
      CALL timer(timematv_st)
      nrb=SIZE(mat%rbl_mat)
      IF (nrb>0) THEN
        nbtype=mat%rbl_mat(1)%nbtype
        IF (nbtype>1) THEN
          pd=SUM(mat%rbl_mat(1)%nb_type(1:2))
        ELSE
          pd=1
        ENDIF
      ELSE
        nbtype=1
        pd=1
      ENDIF
      ntb=SIZE(mat%tbl_mat)
c-----------------------------------------------------------------------
c     if interior nodes are eliminated, basis loops only cover grid
c     and element side centered nodes.
c-----------------------------------------------------------------------
      IF (mat%eliminated) THEN
        typ_max=MIN(nbtype,3_i4)
      ELSE
        typ_max=nbtype
      ENDIF
c-----------------------------------------------------------------------
c     rblocks: call the structured-data matrix operation for each
c     basis type to basis type pair, and load seam array.
c     every basis-type pair now uses a separate multiplication routine
c     for optimization.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
        my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
        DO ityp=1,typ_max
          new_product=.true.
          nqp=mat%rbl_mat(ibl)%nq_type(ityp)
          SELECT CASE(ityp)
          CASE(1)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgg_comp_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechg_comp_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvg_comp_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecig_comp_rbl(product(ibl)%arr,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(2)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgh_comp_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechh_comp_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvh_comp_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecih_comp_rbl(product(ibl)%arrh,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(3)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgv_comp_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechv_comp_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvv_comp_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matveciv_comp_rbl(product(ibl)%arrv,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          CASE(4)
            DO jtyp=1,typ_max
              nqo=mat%rbl_mat(ibl)%nq_type(jtyp)
              mat_ptr=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
              SELECT CASE(jtyp)
              CASE(1)
                CALL matvecgi_comp_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arr,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(2)
                CALL matvechi_comp_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arrh,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(3)
                CALL matvecvi_comp_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arrv,
     $                                 mx,my,nqp,nqo,new_product)
              CASE(4)
                CALL matvecii_comp_rbl(product(ibl)%arri,mat_ptr,
     $                                 operand(ibl)%arri,
     $                                 mx,my,nqp,nqo,new_product)
              END SELECT
              new_product=.false.
            ENDDO
          END SELECT
        ENDDO
        IF (netw)
     $    CALL edge_load_2D_carr(product(ibl),nq,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     tblocks: call the unstructured-data matrix operation and load
c     seam array.
c-----------------------------------------------------------------------
      DO ibl=nrb+1,nrb+ntb
        mv=SIZE(mat%tbl_mat(ibl)%lmat)-1
        pr_ptr=>product(ibl)%arr(:,:,:)
        pr_ptr=0
        op_ptr=>operand(ibl)%arr(:,:,:)
        CALL matvec_comp_tbl(pr_ptr,mat%tbl_mat(ibl)%lmat,op_ptr,mv,nq)
        IF (netw)
     $    CALL edge_load_2D_carr(product(ibl),nq,pd-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     communicate across block boundaries and unload seam.
c-----------------------------------------------------------------------
      IF (netw) THEN
        CALL edge_network(nq,1_i4,pd-1_i4,.false.)
        DO ibl=1,nrb+ntb
          CALL edge_unload_2D_carr(product(ibl),nq,pd-1_i4,seam(ibl))
        ENDDO
      ENDIF
      CALL timer(timematv_en)
      time_matvec=time_matvec+timematv_en-timematv_st
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvec_2D_comp
c-----------------------------------------------------------------------
c     subprogram 21. matvecgg_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgg_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,-1:,:,0:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,0:,0:,iq,0,0)*vector(:,:1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,:0,0:,iq,mx,0)*vector(:,mx-1:,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,0:,:0,iq,0,my)*vector(:,:1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,:0,:0,iq,mx,my)*vector(:,mx-1:,my-1:my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,0:,0:,iq,0,0)*vector(:,:1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,:0,0:,iq,mx,0)*vector(:,mx-1:,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,0:,:0,iq,0,my)*vector(:,:1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,:0,:0,iq,mx,my)*vector(:,mx-1:,my-1:my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgg_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 22. matvecgh_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid
c     vertex to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgh_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,-1:,:,1:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix,my-1:my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix-1:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix-1:ix,my-1:my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgh_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 23. matvecgv_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgv_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,-1:,:,0:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy-1:iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgv_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 24. matvecgi_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles grid vertex
c     to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecgi_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,-1:,:,1:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecgi_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 25. matvechg_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechg_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,-1:,:,0:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,1,0:,iq,0,0)*vector(:,1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,0,0:,iq,mx,0)*vector(:,mx,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,1,:0,iq,0,my)*vector(:,1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,0,:0,iq,mx,my)*vector(:,mx,my-1:my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,1,0:,iq,0,0)*vector(:,1,:1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix+1,:1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,0,0:,iq,mx,0)*vector(:,mx,:1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,1,:0,iq,0,my)*vector(:,1,my-1:my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix+1,my-1:my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,0,:0,iq,mx,my)*vector(:,mx,my-1:my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechg_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 26. matvechh_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechh_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,-1:,:,1:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix,my-1:my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,0:,iq,ix,0)*vector(:,ix:ix,:1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,:0,iq,ix,my)*vector(:,ix:ix,my-1:my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechh_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 27. matvechv_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechv_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,-1:,:,0:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy-1:iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy-1:iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy-1:iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechv_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 28. matvechi_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles horizontal
c     side to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvechi_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,-1:,:,1:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,0:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy-1:iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvechi_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 29. matvecvg_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvg_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,0:,:,0:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,0:,1,iq,0,0)*vector(:,:1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,:0,1,iq,mx,0)*vector(:,mx-1:,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,0:,0,iq,0,my)*vector(:,:1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,:0,0,iq,mx,my)*vector(:,mx-1:,my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,0:,1,iq,0,0)*vector(:,:1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,:0,1,iq,mx,0)*vector(:,mx-1:,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,:,iq,0,iy)*vector(:,:1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,:,iq,mx,iy)*vector(:,mx-1:,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,0:,0,iq,0,my)*vector(:,:1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,:0,0,iq,mx,my)*vector(:,mx-1:,my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvg_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 30. matvecvh_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvh_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,0:,:,1:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix,my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix-1:ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix-1:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix-1:ix,my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvh_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 31. matvecvv_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvv_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,0:,:,0:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,0:,0,iq,0,iy)*vector(:,:1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,:0,0,iq,mx,iy)*vector(:,mx-1:,iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,0:,0,iq,0,iy)*vector(:,:1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,:0,0,iq,mx,iy)*vector(:,mx-1:,iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvv_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 32. matvecvi_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles vertical
c     side to interior operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecvi_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,-1:,0:,:,1:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,0:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix,iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix-1:ix,iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecvi_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 33. matvecig_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to grid vertex operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecig_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,0:,:,0:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create or sum result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=
     $      SUM(matrix(:,1,1,iq,0,0)*vector(:,1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=
     $      SUM(matrix(:,0,1,iq,mx,0)*vector(:,mx,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=
     $      SUM(matrix(:,1,0,iq,0,my)*vector(:,1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=
     $      SUM(matrix(:,0,0,iq,mx,my)*vector(:,mx,my))
        ENDDO

      ELSE res_if
c-----------------------------------------------------------------------
c       bottom
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,0)=result(iq,0,0)+
     $      SUM(matrix(:,1,1,iq,0,0)*vector(:,1,1))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,:,1,iq,ix,0)*vector(:,ix:ix+1,1))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,0)=result(iq,mx,0)+
     $      SUM(matrix(:,0,1,iq,mx,0)*vector(:,mx,1))
        ENDDO
c-----------------------------------------------------------------------
c       interior
c-----------------------------------------------------------------------
        DO iy=1,my-1
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,:,iq,0,iy)*vector(:,1,iy:iy+1))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $         SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix+1,iy:iy+1))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,:,iq,mx,iy)*vector(:,mx,iy:iy+1))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       top
c-----------------------------------------------------------------------
        DO iq=1,nqr
          result(iq,0,my)=result(iq,0,my)+
     $      SUM(matrix(:,1,0,iq,0,my)*vector(:,1,my))
        ENDDO
        DO ix=1,mx-1
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,:,0,iq,ix,my)*vector(:,ix:ix+1,my))
          ENDDO
        ENDDO
        DO iq=1,nqr
          result(iq,mx,my)=result(iq,mx,my)+
     $      SUM(matrix(:,0,0,iq,mx,my)*vector(:,mx,my))
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecig_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 34. matvecih_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to horizontal side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecih_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,0:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,0:,:,1:,0:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=
     $        SUM(matrix(:,0,1,iq,ix,0)*vector(:,ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=
     $        SUM(matrix(:,0,0,iq,ix,my)*vector(:,ix,my))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,0)=result(iq,ix,0)+
     $        SUM(matrix(:,0,1,iq,ix,0)*vector(:,ix,1))
          ENDDO
        ENDDO
        DO iy=1,my-1
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,:,iq,ix,iy)*vector(:,ix:ix,iy:iy+1))
            ENDDO
          ENDDO
        ENDDO
        DO ix=1,mx
          DO iq=1,nqr
            result(iq,ix,my)=result(iq,ix,my)+
     $        SUM(matrix(:,0,0,iq,ix,my)*vector(:,ix,my))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecih_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 35. matveciv_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to vertical side operations only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matveciv_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,0:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,0:,:,0:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=
     $        SUM(matrix(:,1,0,iq,0,iy)*vector(:,1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=
     $        SUM(matrix(:,0,0,iq,mx,iy)*vector(:,mx,iy))
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO iq=1,nqr
            result(iq,0,iy)=result(iq,0,iy)+
     $        SUM(matrix(:,1,0,iq,0,iy)*vector(:,1,iy))
          ENDDO
          DO ix=1,mx-1
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,:,0,iq,ix,iy)*vector(:,ix:ix+1,iy))
            ENDDO
          ENDDO
          DO iq=1,nqr
            result(iq,mx,iy)=result(iq,mx,iy)+
     $        SUM(matrix(:,0,0,iq,mx,iy)*vector(:,mx,iy))
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matveciv_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 36. matvecii_comp_rbl.
c     perform a comp matrix/vector multiplication for an rblock and
c     return the product in the array result.  this handles interior
c     to interior only.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvecii_comp_rbl(result,matrix,vector,
     $                             mx,my,nqr,nqv,new_result)

      INTEGER(i4), INTENT(IN) :: mx,my,nqr,nqv
      COMPLEX(r8), DIMENSION(nqr,1:mx,1:*), INTENT(INOUT) :: result
      COMPLEX(r8), DIMENSION(:,0:,0:,:,1:,1:), INTENT(IN) ::
     $             matrix
      COMPLEX(r8), DIMENSION(nqv,1:mx,1:*), INTENT(IN) :: vector
      LOGICAL, INTENT(IN) :: new_result

      INTEGER(i4) :: ix,iy,iq
c-----------------------------------------------------------------------
c     create result vector.
c-----------------------------------------------------------------------
      res_if: IF (new_result) THEN
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=
     $         SUM(matrix(:,0,0,iq,ix,iy)*vector(:,ix,iy))
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     sum result vector.
c-----------------------------------------------------------------------
      ELSE res_if
        DO iy=1,my
          DO ix=1,mx
            DO iq=1,nqr
              result(iq,ix,iy)=result(iq,ix,iy)+
     $          SUM(matrix(:,0,0,iq,ix,iy)*vector(:,ix,iy))
            ENDDO
          ENDDO
        ENDDO
      ENDIF res_if
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvecii_comp_rbl
c-----------------------------------------------------------------------
c     subprogram 37. matvec_comp_tbl.
c     perform a matrix/vector multiplication for a tblock and return
c     the product in the array result.
c     complex algebra version.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matvec_comp_tbl(result,matrix,vector,mv,nqty)

      INTEGER(i4), INTENT(IN) :: mv,nqty
      TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: matrix
      COMPLEX(r8), DIMENSION(:,0:,0:), INTENT(IN) :: vector
      COMPLEX(r8), DIMENSION(:,0:,0:), INTENT(INOUT) :: result

      INTEGER(i4) :: iqty,jqty,ivert,inbr,mvert
      COMPLEX(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: matel
c-----------------------------------------------------------------------
c     find the matrix/vector product.
c-----------------------------------------------------------------------
      DO ivert=0,mv
        matel=>matrix(ivert)%element
        DO inbr=0,SIZE(matrix(ivert)%from_vert)-1
          DO iqty=1,nqty
            result(iqty,ivert,0)=result(iqty,ivert,0)
     $        +SUM(matel(:,iqty,inbr)
     $            *vector(1:nqty,matrix(ivert)%from_vert(inbr),0))
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matvec_comp_tbl
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matvec_comp_mod


c-----------------------------------------------------------------------
c     matelim_mod contains routines for reducing matrices by eliminating
c     connections to cell-interior data.
c-----------------------------------------------------------------------
      MODULE matelim_mod
      USE local
      USE matrix_type_mod
      USE vector_type_mod
      USE matvec_real_mod
      USE matvec_comp_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 38. matelim_real_inv_int.
c     invert the connections within cell interiors for basis functions
c     of degree 2 or more.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_real_inv_int(mat_str,nq)
      USE math_tran

      TYPE(global_matrix_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq

      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr,int_ptr
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: dense
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: xx,id
      INTEGER(i4) :: ibl,nrb,pd,pd2,iof,jof,nqint,
     $               mxb,myb,ix,iy,iq,jq,ixo,iyo,jxo,jyo,kmat
      LOGICAL :: sing
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
      nqint=mat_str%rbl_mat(1)%nq_type(4)
c-----------------------------------------------------------------------
c-PRE for each rblock, factor the interior to interior matrix.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        ALLOCATE(dense(nqint,nqint,mxb,myb))
        ALLOCATE(xx(nqint,mxb,myb))
        ALLOCATE(id(nqint,mxb,myb))
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        dense=mat_ptr(:,0,0,:,:,:)
        IF (mat_str%symmetric) THEN
          CALL math_solve_q1_sym(nqint,mxb*myb,dense,xx,id,'factor',
     $                           sing)
        ELSE
          CALL math_solve_q1_nxn(nqint,mxb*myb,dense,xx,id,'factor',
     $                           sing)
        ENDIF
        IF (sing) CALL nim_stop
     $    ('Matelim_real_inv_int: dense interior does not factor.')
c-----------------------------------------------------------------------
c       find and save each column of the inverse.
c-----------------------------------------------------------------------
        DO jq=1,nqint
          id=0
          id(jq,:,:)=1
          IF (mat_str%symmetric) THEN
            CALL math_solve_q1_sym(nqint,mxb*myb,dense,xx,id,'solve',
     $                             sing)
          ELSE
            CALL math_solve_q1_nxn(nqint,mxb*myb,dense,xx,id,'solve',
     $                             sing)
          ENDIF
          mat_ptr(jq,0,0,:,:,:)=xx
        ENDDO
        DEALLOCATE(dense)
c-----------------------------------------------------------------------
c       if nq is zero, there are no bases that are continuous across
c       element borders.       
c-----------------------------------------------------------------------
        IF (nq==0) THEN
          DEALLOCATE(id,xx)
          CYCLE
        ENDIF
c-----------------------------------------------------------------------
c       now create the Schur complement, A_oo-A_oi.A_ii**-1.A_io,
c       where i refers to interior data, and o refers to other data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c       start with a computation of A_ii**-1.A_io taking one grid
c       offset and column vector-component-index (but all ix and iy) at
c       a time to use matrix/vector algebra. 
c-----------------------------------------------------------------------
        DO jyo=-1,0
          DO jxo=-1,0
            DO jq=1,nq
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,4)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
              xx=mat_ptr(jq,jxo,jyo,:,:,:)
              CALL matvecii_real_rbl(id,int_ptr,xx,mxb,myb,
     $                               nqint,nqint,.true.)
c-----------------------------------------------------------------------
c             complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c             one column (but all ix and iy) of A_io.
c
c             modify the rows connecting to grid nodes (i.e. A_oo=A_gg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,1)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,nq
                    DO iyo=0,1
                      jof=jyo+iyo
                      DO ixo=0,1
                        iof=jxo+ixo
                        mat_ptr(jq,iof,jof,iq,ix-ixo,iy-iyo)=
     $                    mat_ptr(jq,iof,jof,iq,ix-ixo,iy-iyo)-
     $                      SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                          id(:,ix,iy))
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
c-----------------------------------------------------------------------
c             modify the rows connecting to horizontal side nodes
c             (i.e. A_oo=A_hg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,2)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                    DO iyo=0,1
                      jof=jyo+iyo
                      mat_ptr(jq,jxo,jof,iq,ix,iy-iyo)=
     $                  mat_ptr(jq,jxo,jof,iq,ix,iy-iyo)-
     $                    SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
c-----------------------------------------------------------------------
c             modify the rows connecting to vertical side nodes
c             (i.e. A_oo=A_vg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,3)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                    DO ixo=0,1
                      iof=jxo+ixo
                      mat_ptr(jq,iof,jyo,iq,ix-ixo,iy)=
     $                  mat_ptr(jq,iof,jyo,iq,ix-ixo,iy)-
     $                    SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       now take one horizontal side node and vector-component-index
c       as a column element.
c-----------------------------------------------------------------------
        DO jyo=-1,0
          DO jq=1,mat_str%rbl_mat(ibl)%nq_type(2)
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,4)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
            xx=mat_ptr(jq,0,jyo,:,:,:)
            CALL matvecii_real_rbl(id,int_ptr,xx,mxb,myb,
     $                             nqint,nqint,.true.)
c-----------------------------------------------------------------------
c           complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c           one column (but all ix and iy) of A_io.
c
c           modify the rows connecting to grid nodes (i.e. A_oo=A_gh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,1)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,nq
                  DO iyo=0,1
                    jof=jyo+iyo
                    DO ixo=0,1
                      mat_ptr(jq,ixo,jof,iq,ix-ixo,iy-iyo)=
     $                  mat_ptr(jq,ixo,jof,iq,ix-ixo,iy-iyo)-
     $                    SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                        id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to horizontal side nodes
c           (i.e. A_oo=A_hh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,2)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                  DO iyo=0,1
                    jof=jyo+iyo
                    mat_ptr(jq,0,jof,iq,ix,iy-iyo)=
     $                mat_ptr(jq,0,jof,iq,ix,iy-iyo)-
     $                  SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to vertical side nodes
c           (i.e. A_oo=A_vh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,3)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                  DO ixo=0,1
                    mat_ptr(jq,ixo,jyo,iq,ix-ixo,iy)=
     $                mat_ptr(jq,ixo,jyo,iq,ix-ixo,iy)-
     $                  SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       finally take one vertical side node and vector-component-index
c       as a column element.
c-----------------------------------------------------------------------
        DO jxo=-1,0
          DO jq=1,mat_str%rbl_mat(ibl)%nq_type(3)
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,4)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
            xx=mat_ptr(jq,jxo,0,:,:,:)
            CALL matvecii_real_rbl(id,int_ptr,xx,mxb,myb,
     $                             nqint,nqint,.true.)
c-----------------------------------------------------------------------
c           complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c           one column (but all ix and iy) of A_io.
c
c           modify the rows connecting to grid nodes (i.e. A_oo=A_gv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,1)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,nq
                  DO iyo=0,1
                    DO ixo=0,1
                      iof=jxo+ixo
                      mat_ptr(jq,iof,iyo,iq,ix-ixo,iy-iyo)=
     $                  mat_ptr(jq,iof,iyo,iq,ix-ixo,iy-iyo)-
     $                    SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                        id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to horizontal side nodes
c           (i.e. A_oo=A_hv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,2)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                  DO iyo=0,1
                    mat_ptr(jq,jxo,iyo,iq,ix,iy-iyo)=
     $                mat_ptr(jq,jxo,iyo,iq,ix,iy-iyo)-
     $                  SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to vertical side nodes
c           (i.e. A_oo=A_vv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,3)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                  DO ixo=0,1
                    iof=jxo+ixo
                    mat_ptr(jq,iof,0,iq,ix-ixo,iy)=
     $                mat_ptr(jq,iof,0,iq,ix-ixo,iy)-
     $                  SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(id,xx)
      ENDDO
c-----------------------------------------------------------------------
c     set flag.
c-----------------------------------------------------------------------
      mat_str%eliminated=.true.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_real_inv_int
c-----------------------------------------------------------------------
c     subprogram 39. matelim_real_presolve.
c     for matrices partitioned into cell-interior / other, find
c     A_ii**-1.b_i and b_o - A_oi.A_ii**-1.b_i (where i means interior
c     and o is other).  assume that matelim_real_inv_int has been
c     called, so A_ii**-1 is available in the A_ii storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_real_presolve(mat_str,operand,product,nq)

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: operand,product
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: ibl,nrb,mxb,myb,nqi,nqh,nqv
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr
      REAL(r8), DIMENSION(:,:,:), POINTER :: op_ptr,pr_ptr
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: op_ptr4,pr_ptr4
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
c-----------------------------------------------------------------------
c     check that A_ii has been inverted.
c-----------------------------------------------------------------------
      IF (.NOT.mat_str%eliminated) CALL nim_stop
     $  ('Matelim_real_presolve: dense interior is not factored.')
c-----------------------------------------------------------------------
c-PRE loop over rblocks, and compute and save A_ii**-1.b_i first.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        nqh=mat_str%rbl_mat(ibl)%nq_type(2)
        nqv=mat_str%rbl_mat(ibl)%nq_type(3)
        nqi=mat_str%rbl_mat(ibl)%nq_type(4)
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        CALL matvecii_real_rbl(product(ibl)%arri,mat_ptr,
     $                         operand(ibl)%arri,
     $                         mxb,myb,nqi,nqi,.true.)
c-----------------------------------------------------------------------
c       if nqh is zero, there are no bases that are continuous across
c       element borders.       
c-----------------------------------------------------------------------
        IF (nqh==0) CYCLE
c-----------------------------------------------------------------------
c       now find b_o-A_oi.A_ii**-1.b_i for each of the "other" types
c       of bases.
c-----------------------------------------------------------------------
        op_ptr4=>product(ibl)%arri
        pr_ptr=>product(ibl)%arr
        pr_ptr=-operand(ibl)%arr
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
        CALL matvecig_real_rbl(pr_ptr,mat_ptr,op_ptr4,
     $                         mxb,myb,nq,nqi,.false.)
        pr_ptr=-pr_ptr

        pr_ptr4=>product(ibl)%arrh
        pr_ptr4=-operand(ibl)%arrh
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
        CALL matvecih_real_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqh,nqi,.false.)
        pr_ptr4=-pr_ptr4

        pr_ptr4=>product(ibl)%arrv
        pr_ptr4=-operand(ibl)%arrv
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
        CALL matveciv_real_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqv,nqi,.false.)
        pr_ptr4=-pr_ptr4
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_real_presolve
c-----------------------------------------------------------------------
c     subprogram 40. matelim_real_postsolve.
c     for matrices partitioned into cell-interior / other,
c     subtract A_ii**-1.A_io.x_o from A_ii**-1.b_i to get x_i.
c     the product contains A_ii**-1.b_i in the interior storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_real_postsolve(mat_str,operand,product,nq)

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: operand,product
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: ibl,nrb,mxb,myb,nqi,nqh,nqv
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr
      REAL(r8), DIMENSION(:,:,:), POINTER :: op_ptr,pr_ptr
      REAL(r8), DIMENSION(:,:,:,:), POINTER :: op_ptr4,pr_ptr4
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
c-----------------------------------------------------------------------
c-PRE loop over rblocks, and compute -A_io.x_o first.  use operand
c     interior arrays for temporary storage.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        nqh=mat_str%rbl_mat(ibl)%nq_type(2)
        nqv=mat_str%rbl_mat(ibl)%nq_type(3)
        nqi=mat_str%rbl_mat(ibl)%nq_type(4)

        pr_ptr4=>operand(ibl)%arri
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,4)%arr
        op_ptr=>operand(ibl)%arr
        CALL matvecgi_real_rbl(pr_ptr4,mat_ptr,op_ptr,
     $                         mxb,myb,nqi,nq,.true.)

        mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,4)%arr
        op_ptr4=>operand(ibl)%arrh
        CALL matvechi_real_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqh,.false.)

        mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,4)%arr
        op_ptr4=>operand(ibl)%arrv
        CALL matvecvi_real_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqv,.false.)

        pr_ptr4=-pr_ptr4
c-----------------------------------------------------------------------
c       now find A_ii**-1.b_i-A_ii**-1.A_io.x_o.
c-----------------------------------------------------------------------
        pr_ptr4=>product(ibl)%arri
        op_ptr4=>operand(ibl)%arri
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        CALL matvecii_real_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqi,.false.)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_real_postsolve
c-----------------------------------------------------------------------
c     subprogram 41. matelim_comp_inv_int.
c     invert the connections within cell interiors for basis functions
c     of degree 2 or more, and form the Schur complement for the
c     rest of the matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_comp_inv_int(mat_str,nq)
      USE math_tran

      TYPE(complex_matrix_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq

      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr,int_ptr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: dense
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: xx,id
      INTEGER(i4) :: ibl,nrb,pd,pd2,iof,jof,nqint,
     $               mxb,myb,ix,iy,iq,jq,ixo,iyo,jxo,jyo,kmat
      LOGICAL :: sing
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
      nqint=mat_str%rbl_mat(1)%nq_type(4)
c-----------------------------------------------------------------------
c-PRE for each rblock, factor the interior to interior matrix.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        ALLOCATE(dense(nqint,nqint,mxb,myb))
        ALLOCATE(xx(nqint,mxb,myb))
        ALLOCATE(id(nqint,mxb,myb))
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        dense=mat_ptr(:,0,0,:,:,:)
        IF (mat_str%hermitian) THEN
          CALL math_solve_q1_herm(nqint,mxb*myb,dense,xx,id,'factor',
     $                            sing)
        ELSE
          CALL math_solve_q1_cnxn(nqint,mxb*myb,dense,xx,id,'factor',
     $                            sing)
        ENDIF
        IF (sing) CALL nim_stop
     $    ('Matelim_comp_inv_int: dense interior does not factor.')
c-----------------------------------------------------------------------
c       find and save each column of the inverse.
c-----------------------------------------------------------------------
        DO jq=1,nqint
          id=0
          id(jq,:,:)=1
          IF (mat_str%hermitian) THEN
            CALL math_solve_q1_herm(nqint,mxb*myb,dense,xx,id,'solve',
     $                              sing)
          ELSE
            CALL math_solve_q1_cnxn(nqint,mxb*myb,dense,xx,id,'solve',
     $                              sing)
          ENDIF
          mat_ptr(jq,0,0,:,:,:)=xx
        ENDDO
        DEALLOCATE(dense)
c-----------------------------------------------------------------------
c       if nq is zero, there are no bases that are continuous across
c       element borders.       
c-----------------------------------------------------------------------
        IF (nq==0) THEN
          DEALLOCATE(id,xx)
          CYCLE
        ENDIF
c-----------------------------------------------------------------------
c       now create the Schur complement, A_oo-A_oi.A_ii**-1.A_io,
c       where i refers to interior data, and o refers to other data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c       start with a computation of A_ii**-1.A_io taking one grid
c       offset and column vector-component-index (but all ix and iy) at
c       a time to use matrix/vector algebra. 
c-----------------------------------------------------------------------
        DO jyo=-1,0
          DO jxo=-1,0
            DO jq=1,nq
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,4)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
              xx=mat_ptr(jq,jxo,jyo,:,:,:)
              CALL matvecii_comp_rbl(id,int_ptr,xx,mxb,myb,
     $                               nqint,nqint,.true.)
c-----------------------------------------------------------------------
c             complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c             one column (but all ix and iy) of A_io.
c
c             modify the rows connecting to grid nodes (i.e. A_oo=A_gg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,1)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,nq
                    DO iyo=0,1
                      jof=jyo+iyo
                      DO ixo=0,1
                        iof=jxo+ixo
                        mat_ptr(jq,iof,jof,iq,ix-ixo,iy-iyo)=
     $                    mat_ptr(jq,iof,jof,iq,ix-ixo,iy-iyo)-
     $                      SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                          id(:,ix,iy))
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
c-----------------------------------------------------------------------
c             modify the rows connecting to horizontal side nodes
c             (i.e. A_oo=A_hg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,2)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                    DO iyo=0,1
                      jof=jyo+iyo
                      mat_ptr(jq,jxo,jof,iq,ix,iy-iyo)=
     $                  mat_ptr(jq,jxo,jof,iq,ix,iy-iyo)-
     $                    SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
c-----------------------------------------------------------------------
c             modify the rows connecting to vertical side nodes
c             (i.e. A_oo=A_vg):
c-----------------------------------------------------------------------
              mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,3)%arr
              int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                    DO ixo=0,1
                      iof=jxo+ixo
                      mat_ptr(jq,iof,jyo,iq,ix-ixo,iy)=
     $                  mat_ptr(jq,iof,jyo,iq,ix-ixo,iy)-
     $                    SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       now take one horizontal side node and vector-component-index
c       as a column element.
c-----------------------------------------------------------------------
        DO jyo=-1,0
          DO jq=1,mat_str%rbl_mat(ibl)%nq_type(2)
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,4)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
            xx=mat_ptr(jq,0,jyo,:,:,:)
            CALL matvecii_comp_rbl(id,int_ptr,xx,mxb,myb,
     $                             nqint,nqint,.true.)
c-----------------------------------------------------------------------
c           complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c           one column (but all ix and iy) of A_io.
c
c           modify the rows connecting to grid nodes (i.e. A_oo=A_gh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,1)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,nq
                  DO iyo=0,1
                    jof=jyo+iyo
                    DO ixo=0,1
                      mat_ptr(jq,ixo,jof,iq,ix-ixo,iy-iyo)=
     $                  mat_ptr(jq,ixo,jof,iq,ix-ixo,iy-iyo)-
     $                    SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                        id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to horizontal side nodes
c           (i.e. A_oo=A_hh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,2)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                  DO iyo=0,1
                    jof=jyo+iyo
                    mat_ptr(jq,0,jof,iq,ix,iy-iyo)=
     $                mat_ptr(jq,0,jof,iq,ix,iy-iyo)-
     $                  SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to vertical side nodes
c           (i.e. A_oo=A_vh):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,3)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                  DO ixo=0,1
                    mat_ptr(jq,ixo,jyo,iq,ix-ixo,iy)=
     $                mat_ptr(jq,ixo,jyo,iq,ix-ixo,iy)-
     $                  SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       finally take one vertical side node and vector-component-index
c       as a column element.
c-----------------------------------------------------------------------
        DO jxo=-1,0
          DO jq=1,mat_str%rbl_mat(ibl)%nq_type(3)
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,4)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
            xx=mat_ptr(jq,jxo,0,:,:,:)
            CALL matvecii_comp_rbl(id,int_ptr,xx,mxb,myb,
     $                             nqint,nqint,.true.)
c-----------------------------------------------------------------------
c           complete the computation of A_oo-A_oi.A_ii**-1.A_io for
c           one column (but all ix and iy) of A_io.
c
c           modify the rows connecting to grid nodes (i.e. A_oo=A_gv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,1)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,nq
                  DO iyo=0,1
                    DO ixo=0,1
                      iof=jxo+ixo
                      mat_ptr(jq,iof,iyo,iq,ix-ixo,iy-iyo)=
     $                  mat_ptr(jq,iof,iyo,iq,ix-ixo,iy-iyo)-
     $                    SUM(int_ptr(:,ixo,iyo,iq,ix-ixo,iy-iyo)*
     $                        id(:,ix,iy))
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to horizontal side nodes
c           (i.e. A_oo=A_hv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,2)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(2)
                  DO iyo=0,1
                    mat_ptr(jq,jxo,iyo,iq,ix,iy-iyo)=
     $                mat_ptr(jq,jxo,iyo,iq,ix,iy-iyo)-
     $                  SUM(int_ptr(:,0,iyo,iq,ix,iy-iyo)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           modify the rows connecting to vertical side nodes
c           (i.e. A_oo=A_vv):
c-----------------------------------------------------------------------
            mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,3)%arr
            int_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
            DO iy=1,myb
              DO ix=1,mxb
                DO iq=1,mat_str%rbl_mat(ibl)%nq_type(3)
                  DO ixo=0,1
                    iof=jxo+ixo
                    mat_ptr(jq,iof,0,iq,ix-ixo,iy)=
     $                mat_ptr(jq,iof,0,iq,ix-ixo,iy)-
     $                  SUM(int_ptr(:,ixo,0,iq,ix-ixo,iy)*id(:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(id,xx)
      ENDDO
c-----------------------------------------------------------------------
c     set flag.
c-----------------------------------------------------------------------
      mat_str%eliminated=.true.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_comp_inv_int
c-----------------------------------------------------------------------
c     subprogram 42. matelim_comp_presolve.
c     for matrices partitioned into cell-interior / other, find
c     A_ii**-1.b_i and b_o - A_oi.A_ii**-1.b_i (where i means interior
c     and o is other).  assume that matelim_comp_inv_int has been
c     called, so A_ii**-1 is available in the A_ii storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_comp_presolve(mat_str,operand,product,nq)

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(INOUT) :: operand,
     $                       product
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: ibl,nrb,mxb,myb,nqi,nqh,nqv
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: op_ptr,pr_ptr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: op_ptr4,pr_ptr4
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
c-----------------------------------------------------------------------
c     check that A_ii has been inverted.
c-----------------------------------------------------------------------
      IF (.NOT.mat_str%eliminated) CALL nim_stop
     $  ('Matelim_comp_presolve: dense interior is not factored.')
c-----------------------------------------------------------------------
c-PRE loop over rblocks, and compute and save A_ii**-1.b_i first.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        nqh=mat_str%rbl_mat(ibl)%nq_type(2)
        nqv=mat_str%rbl_mat(ibl)%nq_type(3)
        nqi=mat_str%rbl_mat(ibl)%nq_type(4)
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        CALL matvecii_comp_rbl(product(ibl)%arri,mat_ptr,
     $                         operand(ibl)%arri,
     $                         mxb,myb,nqi,nqi,.true.)
c-----------------------------------------------------------------------
c       if nqh is zero, there are no bases that are continuous across
c       element borders.       
c-----------------------------------------------------------------------
        IF (nqh==0) CYCLE
c-----------------------------------------------------------------------
c       now find b_o-A_oi.A_ii**-1.b_i for each of the "other" types
c       of bases.
c-----------------------------------------------------------------------
        op_ptr4=>product(ibl)%arri
        pr_ptr=>product(ibl)%arr
        pr_ptr=-operand(ibl)%arr
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,1)%arr
        CALL matvecig_comp_rbl(pr_ptr,mat_ptr,op_ptr4,
     $                         mxb,myb,nq,nqi,.false.)
        pr_ptr=-pr_ptr

        pr_ptr4=>product(ibl)%arrh
        pr_ptr4=-operand(ibl)%arrh
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,2)%arr
        CALL matvecih_comp_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqh,nqi,.false.)
        pr_ptr4=-pr_ptr4

        pr_ptr4=>product(ibl)%arrv
        pr_ptr4=-operand(ibl)%arrv
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,3)%arr
        CALL matveciv_comp_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqv,nqi,.false.)
        pr_ptr4=-pr_ptr4
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_comp_presolve
c-----------------------------------------------------------------------
c     subprogram 43. matelim_comp_postsolve.
c     for matrices partitioned into cell-interior / other,
c     subtract A_ii**-1.A_io.x_o from A_ii**-1.b_i to get x_i.
c     the product contains A_ii**-1.b_i in the interior storage.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matelim_comp_postsolve(mat_str,operand,product,nq)

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      TYPE(cvector_2D_type), DIMENSION(:), INTENT(INOUT) :: operand,
     $                       product
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: ibl,nrb,mxb,myb,nqi,nqh,nqv
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: mat_ptr
      COMPLEX(r8), DIMENSION(:,:,:), POINTER :: op_ptr,pr_ptr
      COMPLEX(r8), DIMENSION(:,:,:,:), POINTER :: op_ptr4,pr_ptr4
c-----------------------------------------------------------------------
c     compute the number of blocks and polynomial degree.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      IF (nrb<1) RETURN
      IF (mat_str%rbl_mat(1)%nbtype==1) RETURN
c-----------------------------------------------------------------------
c-PRE loop over rblocks, and compute -A_io.x_o first.  use operand
c     interior arrays for temporary storage.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,5)
        myb=SIZE(mat_str%rbl_mat(ibl)%mat(4,4)%arr,6)
        nqh=mat_str%rbl_mat(ibl)%nq_type(2)
        nqv=mat_str%rbl_mat(ibl)%nq_type(3)
        nqi=mat_str%rbl_mat(ibl)%nq_type(4)

        pr_ptr4=>operand(ibl)%arri
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(1,4)%arr
        op_ptr=>operand(ibl)%arr
        CALL matvecgi_comp_rbl(pr_ptr4,mat_ptr,op_ptr,
     $                         mxb,myb,nqi,nq,.true.)

        mat_ptr=>mat_str%rbl_mat(ibl)%mat(2,4)%arr
        op_ptr4=>operand(ibl)%arrh
        CALL matvechi_comp_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqh,.false.)

        mat_ptr=>mat_str%rbl_mat(ibl)%mat(3,4)%arr
        op_ptr4=>operand(ibl)%arrv
        CALL matvecvi_comp_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqv,.false.)

        pr_ptr4=-pr_ptr4
c-----------------------------------------------------------------------
c       now find A_ii**-1.b_i-A_ii**-1.A_io.x_o.
c-----------------------------------------------------------------------
        pr_ptr4=>product(ibl)%arri
        op_ptr4=>operand(ibl)%arri
        mat_ptr=>mat_str%rbl_mat(ibl)%mat(4,4)%arr
        CALL matvecii_comp_rbl(pr_ptr4,mat_ptr,op_ptr4,
     $                         mxb,myb,nqi,nqi,.false.)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matelim_comp_postsolve
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matelim_mod


c-----------------------------------------------------------------------
c     generic module selects the appropriate routine for real or 
c     complex algebra and contains some utility routines for matrices.
c-----------------------------------------------------------------------
      MODULE matrix_mod
      USE local
      USE matrix_type_mod
      USE matvec_real_mod
      USE matvec_comp_mod
      USE matelim_mod
      IMPLICIT NONE

      INTERFACE matvec
        MODULE PROCEDURE matvec_real,matvec_comp,matvec_2D_comp
      END INTERFACE

      INTERFACE matrix_lump
        MODULE PROCEDURE matrix_lump_real_rmat,matrix_lump_real_tmat,
     $                   matrix_lump_comp_rmat,matrix_lump_comp_tmat
      END INTERFACE

      INTERFACE matelim_presolve
        MODULE PROCEDURE matelim_real_presolve,matelim_comp_presolve
      END INTERFACE

      INTERFACE matelim_postsolve
        MODULE PROCEDURE matelim_real_postsolve,matelim_comp_postsolve
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 44. matrix_lump_real_rmat
c     lump the rblock matrix stored in mat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_lump_real_rmat(mat_str,nq,lumpleft)

      TYPE(rbl_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN), OPTIONAL :: lumpleft

      INTEGER(i4) :: ix,iy,ityp,jtyp,imat,jmat,ix0,iy0,jx0,jy0,iq,jq,
     $               iqo,jqoo,jqdo,ixs
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: omat,dmat
      REAL(r8) :: tmp
c-----------------------------------------------------------------------
c     loop over matrices for different basis function types and bases
c     of each type.
c-----------------------------------------------------------------------
      DO ityp=1,SIZE(mat_str%mat,2)
        ix0=mat_str%ix0(ityp)
        ixs=ix0
        IF (PRESENT(lumpleft)) THEN
          IF (.NOT.lumpleft) ixs=1
        ENDIF
        iy0=mat_str%iy0(ityp)
        dmat=>mat_str%mat(ityp,ityp)%arr
        DO jtyp=1,SIZE(mat_str%mat,1)
          jx0=mat_str%ix0(jtyp)
          jy0=mat_str%iy0(jtyp)
          omat=>mat_str%mat(jtyp,ityp)%arr
          DO imat=1,mat_str%nb_type(ityp)
            iqo=nq*(imat-1)
            DO jmat=1,mat_str%nb_type(jtyp)
              jqoo=nq*(jmat-1)
              jqdo=nq*(imat-1)
c-----------------------------------------------------------------------
c             loop over grid points, and lump the operator.
c-----------------------------------------------------------------------
              IF (imat==jmat.AND.ityp==jtyp) THEN
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO ix=ixs,SIZE(dmat,5)-(1-ix0)
                    DO iq=iqo+1,iqo+nq
                      DO jq=iqo+1,iqo+nq
                        tmp=SUM(dmat(jq,:,:,iq,ix,iy))
                        dmat(jq,:,:,iq,ix,iy)=0
                        dmat(jq,0,0,iq,ix,iy)=tmp
                      ENDDO
                    ENDDO
                  ENDDO
                  IF (ixs>ix0) THEN
                    DO iq=iqo+1,iqo+nq
                      DO jq=iqo+1,iqo+nq
                        tmp=SUM(dmat(jq,0,:,iq,ix0,iy))
                        dmat(jq,:,:,iq,ix0,iy)=0
                        dmat(jq,0,0,iq,ix0,iy)=tmp
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ELSE
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO ix=ixs,SIZE(dmat,5)-(1-ix0)
                    DO iq=iqo+1,iqo+nq
                      DO jq=1,nq
                        tmp=SUM(omat(jqoo+jq,:,:,iq,ix,iy))
                        omat(jqoo+jq,:,:,iq,ix,iy)=0
                        dmat(jqdo+jq,0,0,iq,ix,iy)=
     $                    dmat(jqdo+jq,0,0,iq,ix,iy)+tmp
                      ENDDO
                    ENDDO
                  ENDDO
                  IF (ixs>ix0) THEN
                    DO iq=iqo+1,iqo+nq
                      DO jq=1,nq
                        tmp=SUM(omat(jqoo+jq,0,:,iq,ix0,iy))
                        omat(jqoo+jq,:,:,iq,ix0,iy)=0
                        dmat(jqdo+jq,0,0,iq,ix0,iy)=
     $                    dmat(jqdo+jq,0,0,iq,ix0,iy)+tmp
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_lump_real_rmat
c-----------------------------------------------------------------------
c     subprogram 45. matrix_lump_real_tmat
c     lump the tblock matrix stored in lmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_lump_real_tmat(mat,nqty)
      USE tri_linear

      TYPE(matrix_element_type3), DIMENSION(0:), INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: iqty,jqty,ivert
c-----------------------------------------------------------------------
c     lump elements into the diagonal and zero out off-diagonals.
c-----------------------------------------------------------------------
      DO ivert=0,SIZE(mat)-1
         mat(ivert)%element(1:nqty,1:nqty,0)=
     $      SUM(mat(ivert)%element(1:nqty,1:nqty,:),3)
         mat(ivert)%element(1:nqty,1:nqty,1:)=0
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_lump_real_tmat
c-----------------------------------------------------------------------
c     subprogram 46. matrix_degen_collect_real
c     collects off-diagonal elements along a degenerate boundary of an
c     rblock and places them in the diagonal elements.  this is done to
c     assist the partial factorizations and does not affect the actual
c     matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_degen_collect_real(mat_str,nq,symmetric)

      TYPE(rbl_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN) :: symmetric

      INTEGER(i4) :: iy,iym,imat,jmat,ityp,jtyp,iy0,jy0,iq,jq,iqo,jqo
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: omat,dmat
      REAL(r8), DIMENSION(nq,nq) :: off
      REAL(r8) :: tmp
c-----------------------------------------------------------------------
c     loop over the different basis types and bases (for each type) with
c     nodes along the degenerate side.
c-----------------------------------------------------------------------
      DO ityp=1,MIN(SIZE(mat_str%mat,2),3_i4),2
        iy0=mat_str%iy0(ityp)
        dmat=>mat_str%mat(ityp,ityp)%arr
        DO jtyp=1,MIN(SIZE(mat_str%mat,1),3_i4),2
          jy0=mat_str%iy0(jtyp)
          omat=>mat_str%mat(jtyp,ityp)%arr
          DO imat=1,mat_str%nb_type(ityp)
            iqo=nq*(imat-1)
            DO jmat=1,mat_str%nb_type(jtyp)
              jqo=nq*(jmat-1)
c-----------------------------------------------------------------------
c             loop over nodes along the degenerate side, and lump matrix
c             elements connecting other nodes along the degenerate side.
c-----------------------------------------------------------------------
              IF (imat==jmat.AND.ityp==jtyp) THEN
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO iq=iqo+1,iqo+nq
                    DO jq=iqo+1,iqo+nq
                      tmp=SUM(dmat(jq,0,:,iq,0,iy))
                      dmat(jq,0,:,iq,0,iy)=0
                      dmat(jq,0,0,iq,0,iy)=tmp
                    ENDDO
                  ENDDO
                  IF (ityp/=1.OR.iy==0) CYCLE
                  iym=iy-1
                  IF (symmetric) THEN
                    off=0.5*(dmat(:,1,-1,:,0,iy)
     $                      +TRANSPOSE(dmat(:,-1,1,:,1,iym)))
                    dmat(:, 1, 0,:,0,iym)=dmat(:, 1,0,:,0,iym)+off
                    dmat(:,-1, 0,:,1,iym)=dmat(:,-1,0,:,1,iym)+
     $                                    TRANSPOSE(off)
                    dmat(:, 1,-1,:,0,iy )=0
                    dmat(:,-1, 1,:,1,iym)=0
                    off=0.5*(dmat(:,1,1,:,0,iym)
     $                      +TRANSPOSE(dmat(:,-1,-1,:,1,iy)))
                    dmat(:, 1, 0,:,0,iy )=dmat(:, 1,0,:,0,iy )+off
                    dmat(:,-1, 0,:,1,iy )=dmat(:,-1,0,:,1,iy )+
     $                                    TRANSPOSE(off)
                    dmat(:, 1, 1,:,0,iym)=0
                    dmat(:,-1,-1,:,1,iy )=0
                  ELSE
                    dmat(:, 1, 0,:,0,iym)=dmat(:,1, 0,:,0,iym)+
     $                                    dmat(:,1,-1,:,0,iy )
                    dmat(:, 1,-1,:,0,iy )=0
                    dmat(:,-1, 0,:,1,iym)=dmat(:,-1,0,:,1,iym)+
     $                                    dmat(:,-1,1,:,1,iym)
                    dmat(:,-1, 1,:,1,iym)=0
                    dmat(:, 1, 0,:,0,iy )=dmat(:,1,0,:,0,iy )+
     $                                    dmat(:,1,1,:,0,iym)
                    dmat(:, 1, 1,:,0,iym)=0
                    dmat(:,-1, 0,:,1,iy )=dmat(:,-1, 0,:,1,iy)+
     $                                    dmat(:,-1,-1,:,1,iy)
                    dmat(:,-1,-1,:,1,iy )=0
                  ENDIF
                ENDDO
              ELSE
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO iq=iqo+1,iqo+nq
                    DO jq=jqo+1,jqo+nq
                      tmp=SUM(omat(jq,0,:,iq,0,iy))
                      omat(jq,0,:,iq,0,iy)=0
                      dmat(iqo-jqo+jq,0,0,iq,0,iy)=
     $                  dmat(iqo-jqo+jq,0,0,iq,0,iy)+tmp
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_degen_collect_real
c-----------------------------------------------------------------------
c     subprogram 47. matrix_lump_comp_rmat
c     lump the complex rblock matrix stored in mat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_lump_comp_rmat(mat_str,nq,lumpleft)

      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN), OPTIONAL :: lumpleft

      INTEGER(i4) :: ix,iy,ityp,jtyp,imat,jmat,ix0,iy0,jx0,jy0,iq,jq,
     $               iqo,jqoo,jqdo,ixs
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: omat,dmat
      COMPLEX(r8) :: tmp
c-----------------------------------------------------------------------
c     loop over matrices for different basis function types and bases
c     of each type.
c-----------------------------------------------------------------------
      DO ityp=1,SIZE(mat_str%mat,2)
        ix0=mat_str%ix0(ityp)
        ixs=ix0
        IF (PRESENT(lumpleft)) THEN
          IF (.NOT.lumpleft) ixs=1
        ENDIF
        iy0=mat_str%iy0(ityp)
        dmat=>mat_str%mat(ityp,ityp)%arr
        DO jtyp=1,SIZE(mat_str%mat,1)
          jx0=mat_str%ix0(jtyp)
          jy0=mat_str%iy0(jtyp)
          omat=>mat_str%mat(jtyp,ityp)%arr
          DO imat=1,mat_str%nb_type(ityp)
            iqo=nq*(imat-1)
            DO jmat=1,mat_str%nb_type(jtyp)
              jqoo=nq*(jmat-1)
              jqdo=nq*(imat-1)
c-----------------------------------------------------------------------
c             loop over grid points, and lump the operator.
c-----------------------------------------------------------------------
              IF (imat==jmat.AND.ityp==jtyp) THEN
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO ix=ixs,SIZE(dmat,5)-(1-ix0)
                    DO iq=iqo+1,iqo+nq
                      DO jq=iqo+1,iqo+nq
                        tmp=SUM(dmat(jq,:,:,iq,ix,iy))
                        dmat(jq,:,:,iq,ix,iy)=0
                        dmat(jq,0,0,iq,ix,iy)=tmp
                      ENDDO
                    ENDDO
                  ENDDO
                  IF (ixs>ix0) THEN
                    DO iq=iqo+1,iqo+nq
                      DO jq=iqo+1,iqo+nq
                        tmp=SUM(dmat(jq,0,:,iq,ix0,iy))
                        dmat(jq,:,:,iq,ix0,iy)=0
                        dmat(jq,0,0,iq,ix0,iy)=tmp
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ELSE
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO ix=ixs,SIZE(dmat,5)-(1-ix0)
                    DO iq=iqo+1,iqo+nq
                      DO jq=1,nq
                        tmp=SUM(omat(jqoo+jq,:,:,iq,ix,iy))
                        omat(jqoo+jq,:,:,iq,ix,iy)=0
                        dmat(jqdo+jq,0,0,iq,ix,iy)=
     $                    dmat(jqdo+jq,0,0,iq,ix,iy)+tmp
                      ENDDO
                    ENDDO
                  ENDDO
                  IF (ixs>ix0) THEN
                    DO iq=iqo+1,iqo+nq
                      DO jq=1,nq
                        tmp=SUM(omat(jqoo+jq,0,:,iq,ix0,iy))
                        omat(jqoo+jq,:,:,iq,ix0,iy)=0
                        dmat(jqdo+jq,0,0,iq,ix0,iy)=
     $                    dmat(jqdo+jq,0,0,iq,ix0,iy)+tmp
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_lump_comp_rmat
c-----------------------------------------------------------------------
c     subprogram 48. matrix_lump_comp_tmat
c     lump the complex tblock matrix stored in lmat.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_lump_comp_tmat(mat,nqty)
      USE tri_linear
      USE matrix_type_mod

      TYPE(comp_matrix_element_type3), DIMENSION(0:),
     $                                 INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: iqty,jqty,ivert
c-----------------------------------------------------------------------
c     lump elements into the diagonal and zero out off-diagonals.
c-----------------------------------------------------------------------
      DO ivert=0,SIZE(mat)-1
         mat(ivert)%element(1:nqty,1:nqty,0)=
     $      SUM(mat(ivert)%element(1:nqty,1:nqty,:),3)
         mat(ivert)%element(1:nqty,1:nqty,1:)=0
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_lump_comp_tmat
c-----------------------------------------------------------------------
c     subprogram 49. matrix_degen_collect_comp
c     collects off-diagonal elements along a degenerate boundary of an
c     rblock and places them in the diagonal elements.  this is done to
c     assist the partial factorizations and does not affect the actual
c     matrix.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_degen_collect_comp(mat_str,nq,hermitian)

      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nq
      LOGICAL, INTENT(IN) :: hermitian

      INTEGER(i4) :: iy,iym,imat,jmat,ityp,jtyp,iy0,jy0,iq,jq,iqo,jqo
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: omat,dmat
      COMPLEX(r8), DIMENSION(nq,nq) :: off
      COMPLEX(r8) :: tmp
c-----------------------------------------------------------------------
c     loop over the different basis types and bases (for each type) with
c     nodes along the degenerate side.
c-----------------------------------------------------------------------
      DO ityp=1,MIN(SIZE(mat_str%mat,2),3_i4),2
        iy0=mat_str%iy0(ityp)
        dmat=>mat_str%mat(ityp,ityp)%arr
        DO jtyp=1,MIN(SIZE(mat_str%mat,1),3_i4),2
          jy0=mat_str%iy0(jtyp)
          omat=>mat_str%mat(jtyp,ityp)%arr
          DO imat=1,mat_str%nb_type(ityp)
            iqo=nq*(imat-1)
            DO jmat=1,mat_str%nb_type(jtyp)
              jqo=nq*(jmat-1)
c-----------------------------------------------------------------------
c             loop over nodes along the degenerate side, and lump matrix
c             elements connecting other nodes along the degenerate side.
c-----------------------------------------------------------------------
              IF (imat==jmat.AND.ityp==jtyp) THEN
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO iq=iqo+1,iqo+nq
                    DO jq=iqo+1,iqo+nq
                      tmp=SUM(dmat(jq,0,:,iq,0,iy))
                      dmat(jq,0,:,iq,0,iy)=0
                      dmat(jq,0,0,iq,0,iy)=tmp
                    ENDDO
                  ENDDO
                  IF (ityp/=1.OR.iy==0) CYCLE
                  iym=iy-1
                  IF (hermitian) THEN
                    off=0.5*(dmat(:,1,-1,:,0,iy)
     $                 +TRANSPOSE(CONJG(dmat(:,-1,1,:,1,iym))))
                    dmat(:, 1, 0,:,0,iym)=dmat(:, 1,0,:,0,iym)+off
                    dmat(:,-1, 0,:,1,iym)=dmat(:,-1,0,:,1,iym)
     $                                   +TRANSPOSE(CONJG(off))
                    dmat(:, 1,-1,:,0,iy )=0
                    dmat(:,-1, 1,:,1,iym)=0
                    off=0.5*(dmat(:,1,1,:,0,iym)
     $                 +TRANSPOSE(CONJG(dmat(:,-1,-1,:,1,iy))))
                    dmat(:, 1, 0,:,0,iy )=dmat(:, 1,0,:,0,iy )+off
                    dmat(:,-1, 0,:,1,iy )=dmat(:,-1,0,:,1,iy )
     $                                   +TRANSPOSE(CONJG(off))
                    dmat(:, 1, 1,:,0,iym)=0
                    dmat(:,-1,-1,:,1,iy )=0
                  ELSE
                    dmat(:, 1, 0,:,0,iym)=dmat(:,1, 0,:,0,iym)+
     $                                    dmat(:,1,-1,:,0,iy )
                    dmat(:, 1,-1,:,0,iy )=0
                    dmat(:,-1, 0,:,1,iym)=dmat(:,-1,0,:,1,iym)+
     $                                    dmat(:,-1,1,:,1,iym)
                    dmat(:,-1, 1,:,1,iym)=0
                    dmat(:, 1, 0,:,0,iy )=dmat(:,1,0,:,0,iy )+
     $                                    dmat(:,1,1,:,0,iym)
                    dmat(:, 1, 1,:,0,iym)=0
                    dmat(:,-1, 0,:,1,iy )=dmat(:,-1, 0,:,1,iy)+
     $                                    dmat(:,-1,-1,:,1,iy)
                    dmat(:,-1,-1,:,1,iy )=0
                  ENDIF
                ENDDO
              ELSE
                DO iy=iy0,SIZE(dmat,6)-(1-iy0)
                  DO iq=iqo+1,iqo+nq
                    DO jq=jqo+1,jqo+nq
                      tmp=SUM(omat(jq,0,:,iq,0,iy))
                      omat(jq,0,:,iq,0,iy)=0
                      dmat(iqo-jqo+jq,0,0,iq,0,iy)=
     $                  dmat(iqo-jqo+jq,0,0,iq,0,iy)+tmp
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        IF (hermitian) THEN
          DO iy=iy0,SIZE(dmat,6)-(1-iy0)
            DO iq=1,SIZE(dmat,1)
              dmat(iq,0,0,iq,0,iy)=REAL(dmat(iq,0,0,iq,0,iy),r8)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_degen_collect_comp
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matrix_mod
