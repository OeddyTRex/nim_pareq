c-----------------------------------------------------------------------
c     file data_dealloc.f:  contains subroutines for deallocating
c     nimrod data structures used by nimplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. deallocate_data.
c     2. bl_rb_dealloc.
c     3. bl_rb_dealloc2.
c     4. bl_tb_dealloc.
c     5. bl_tb_dealloc2.
c     6. bl_seam_dealloc.
c-----------------------------------------------------------------------
c     subprogram 1. deallocate_data
c     deallocates the seam, rb and tb data structures.
c-----------------------------------------------------------------------
      SUBROUTINE deallocate_data
      USE local
      USE fields
      USE computation_pointers
      USE matrix_storage_mod
      USE seam_storage_mod
      USE rblock
      USE tblock
      USE surface
      USE input
      USE input_eq
      USE iter_cg
      USE regularity
      USE nim_locate
      USE nimeq_all
      USE nimeq_mod
      IMPLICIT NONE

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     deallocate data pointers.
c-----------------------------------------------------------------------
      DEALLOCATE(be,be_n0,ja,ve,pres,pres_n0,nd_n0,prese,nd,conc,
     $  work1,work2,work3,be_eq,ja_eq,ve_eq,pres_eq,prese_eq,
     $  tele_eq,tion_eq,nd_eq,diff_shape,si_nl_pres,rwork1,rwork2,
     $  rwork3,rwork4,tele,tion,pflux,fllen)
c-----------------------------------------------------------------------
c     deallocate seam data.
c-----------------------------------------------------------------------
      CALL bl_seam_dealloc(seam0)
      DO ibl=1,nbl
        CALL bl_seam_dealloc(seam(ibl))
      ENDDO
      DEALLOCATE(seam)
      DEALLOCATE(exblock_list)
      DEALLOCATE(r0block_list)
c-----------------------------------------------------------------------
c     deallocate rblock data.
c-----------------------------------------------------------------------
      IF (nrbl>0) THEN
        DO ibl=1,nrbl
          CALL bl_rb_dealloc(rb(ibl))
          CALL bl_rb_dealloc2(rb(ibl))
          CALL rblock_basis_dealloc(rb(ibl))
          CALL vector_type_dealloc(sln(ibl))
          CALL vector_type_dealloc(csln(ibl))
          CALL vector_type_dealloc(rhs(ibl))
          CALL vector_type_dealloc(crhs(ibl))
          CALL vector_type_dealloc(diperr(ibl))
          CALL matrix_rbl_dealloc(mass_mat%rbl_mat(ibl))
          CALL matrix_rbl_dealloc(delstar_mat%rbl_mat(ibl))
          IF (gs_type=='free')
     $      CALL matrix_rbl_dealloc(delstlj_mat%rbl_mat(ibl))
          IF (btop_check/='none')
     $      CALL matrix_rbl_dealloc(fladv_mat%rbl_mat(ibl))
          CALL matrix_rbl_dealloc(bfield_mat%rbl_mat(ibl))
          DEALLOCATE(rb(ibl)%xg)
          DEALLOCATE(rb(ibl)%yg)
          DEALLOCATE(rb(ibl)%wg)
        ENDDO
        DEALLOCATE(rb)
        CALL nim_rb_loc_dealloc
      ENDIF
c-----------------------------------------------------------------------
c     deallocate tblock data.
c-----------------------------------------------------------------------
      IF (nbl>nrbl) THEN
        DO ibl=nrbl+1,nbl
          CALL bl_tb_dealloc(tb(ibl))
          CALL bl_tb_dealloc2(tb(ibl))
          CALL tblock_basis_dealloc(tb(ibl)%tgeom)
          CALL vector_type_dealloc(sln(ibl))
          CALL vector_type_dealloc(csln(ibl))
          CALL vector_type_dealloc(rhs(ibl))
          CALL vector_type_dealloc(crhs(ibl))
          CALL vector_type_dealloc(diperr(ibl))
          CALL matrix_tbl_dealloc(mass_mat%tbl_mat(ibl))
          CALL matrix_tbl_dealloc(delstar_mat%tbl_mat(ibl))
          IF (gs_type=='free')
     $      CALL matrix_tbl_dealloc(delstlj_mat%tbl_mat(ibl))
          IF (btop_check/='none')
     $      CALL matrix_tbl_dealloc(fladv_mat%tbl_mat(ibl))
          CALL matrix_tbl_dealloc(bfield_mat%tbl_mat(ibl))
          DEALLOCATE(tb(ibl)%wg)
        ENDDO
        DEALLOCATE(tb)
      ENDIF
c-----------------------------------------------------------------------
c     deallocate surface data.
c-----------------------------------------------------------------------
      CALL surface_dealloc
c-----------------------------------------------------------------------
c     deallocate computation pointers.
c-----------------------------------------------------------------------
      DEALLOCATE(rhs,crhs,sln,csln,diperr)
c-----------------------------------------------------------------------
c     deallocate preconditioning factor structures.
c-----------------------------------------------------------------------
      CALL iter_fac_dealloc(mass_fac,nrbl,nimeq_solver)
      CALL iter_fac_dealloc(delstar_fac,nrbl,nimeq_solver)
      IF (gs_type=='free')
     $  CALL iter_fac_dealloc(delstlj_fac,nrbl,nimeq_solver)
      IF (btop_check/='none')
     $  CALL iter_fac_dealloc(fladv_fac,nrbl,nimeq_solver)
      CALL iter_fac_dealloc(bfield_fac,nrbl,nimeq_solver)
c-----------------------------------------------------------------------
c     deallocate global data structures and data structures for
c     parallel communication.
c-----------------------------------------------------------------------
      IF (nprocs>1) CALL nimeq_all_dealloc
      CALL parallel_block_dealloc
      CALL parallel_seam_dealloc
      CALL parallel_seg_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE deallocate_data
c-----------------------------------------------------------------------
c     subprogram 2. bl_rb_dealloc.
c     deallocates space for a single rblock.
c-----------------------------------------------------------------------
      SUBROUTINE bl_rb_dealloc(rb)
      USE local
      USE rblock_type_mod
      IMPLICIT NONE

      TYPE(rblock_type), INTENT(INOUT) :: rb

      CALL lagr_quad_dealloc(rb%rz)
      CALL lagr_quad_dealloc(rb%be_eq)
      CALL lagr_quad_dealloc(rb%be)
      CALL lagr_quad_dealloc(rb%ja_eq)
      CALL lagr_quad_dealloc(rb%ve_eq)
      CALL lagr_quad_dealloc(rb%ve)
      CALL lagr_quad_dealloc(rb%pres_eq)
      CALL lagr_quad_dealloc(rb%pres)
      CALL lagr_quad_dealloc(rb%prese_eq)
      CALL lagr_quad_dealloc(rb%prese)
      CALL lagr_quad_dealloc(rb%nd_eq)
      CALL lagr_quad_dealloc(rb%nd)
      CALL lagr_quad_dealloc(rb%diff_shape)
      CALL lagr_quad_dealloc(rb%conc)
      CALL lagr_quad_dealloc(rb%tele_eq)
      CALL lagr_quad_dealloc(rb%tele)
      CALL lagr_quad_dealloc(rb%tion_eq)
      CALL lagr_quad_dealloc(rb%tion)

      RETURN
      END SUBROUTINE bl_rb_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. bl_rb_dealloc2.
c     deallocates non-fundamental arrays for an rblock.
c-----------------------------------------------------------------------
      SUBROUTINE bl_rb_dealloc2(rb)
      USE local
      USE rblock
      IMPLICIT NONE

      TYPE(rblock_type), INTENT(INOUT) :: rb

      INTEGER(i4) :: imat

      CALL lagr_quad_dealloc(rb%ja)
      CALL lagr_quad_dealloc(rb%work1)
      CALL lagr_quad_dealloc(rb%work2)
      CALL lagr_quad_dealloc(rb%work3)
      CALL lagr_quad_dealloc(rb%rwork1)
      CALL lagr_quad_dealloc(rb%rwork2)
      CALL lagr_quad_dealloc(rb%rwork3)
      CALL lagr_quad_dealloc(rb%rwork4)
      CALL lagr_quad_dealloc(rb%pflux)
      CALL lagr_quad_dealloc(rb%fllen)
      CALL lagr_quad_dealloc(rb%be_n0)
      CALL lagr_quad_dealloc(rb%pres_n0)
      CALL lagr_quad_dealloc(rb%nd_n0)
      CALL rblock_qp_dealloc(rb%qrwork1)
      CALL rblock_qp_dealloc(rb%qrwork4)
      CALL rblock_qp_dealloc(rb%qpres_eq)
      CALL rblock_qp_dealloc(rb%qja_eq)
      CALL rblock_qp_dealloc(rb%qpres_n0)
      CALL rblock_qp_dealloc(rb%qnd_n0)
     
      RETURN
      END SUBROUTINE bl_rb_dealloc2
c-----------------------------------------------------------------------
c     subprogram 4. bl_tb_dealloc.
c     deallocates space for a single tblock.
c-----------------------------------------------------------------------
      SUBROUTINE bl_tb_dealloc(tb)
      USE local
      USE tblock_type_mod
      IMPLICIT NONE

      TYPE(tblock_type), INTENT(INOUT) :: tb

      CALL tri_linear_geom_dealloc(tb%tgeom)
      CALL tri_linear_dealloc(tb%be_eq)
      CALL tri_linear_dealloc(tb%be)
      CALL tri_linear_dealloc(tb%ja_eq)
      CALL tri_linear_dealloc(tb%ve_eq)
      CALL tri_linear_dealloc(tb%ve)
      CALL tri_linear_dealloc(tb%pres_eq)
      CALL tri_linear_dealloc(tb%pres)
      CALL tri_linear_dealloc(tb%prese_eq)
      CALL tri_linear_dealloc(tb%prese)
      CALL tri_linear_dealloc(tb%nd_eq)
      CALL tri_linear_dealloc(tb%nd)
      CALL tri_linear_dealloc(tb%diff_shape)
      CALL tri_linear_dealloc(tb%conc)
      CALL tri_linear_dealloc(tb%tele_eq)
      CALL tri_linear_dealloc(tb%tele)
      CALL tri_linear_dealloc(tb%tion_eq)
      CALL tri_linear_dealloc(tb%tion)

      RETURN
      END SUBROUTINE bl_tb_dealloc
c-----------------------------------------------------------------------
c     subprogram 5. bl_tb_dealloc2.
c     deallocates non-fundamental arrays for a tblock.
c-----------------------------------------------------------------------
      SUBROUTINE bl_tb_dealloc2(tb)
      USE local
      USE tblock
      IMPLICIT NONE

      TYPE(tblock_type), INTENT(INOUT) :: tb

      CALL tri_linear_dealloc(tb%ja)
      CALL tri_linear_dealloc(tb%work1)
      CALL tri_linear_dealloc(tb%work2)
      CALL tri_linear_dealloc(tb%work3)
      CALL tri_linear_dealloc(tb%rwork1)
      CALL tri_linear_dealloc(tb%rwork2)
      CALL tri_linear_dealloc(tb%rwork3)
      CALL tri_linear_dealloc(tb%rwork4)
      CALL tri_linear_dealloc(tb%pflux)
      CALL tri_linear_dealloc(tb%fllen)
      CALL tri_linear_dealloc(tb%be_n0)
      CALL tri_linear_dealloc(tb%pres_n0)
      CALL tri_linear_dealloc(tb%nd_n0)
      CALL tblock_qp_dealloc(tb%qrwork1)
      CALL tblock_qp_dealloc(tb%qrwork4)
      CALL tblock_qp_dealloc(tb%qpres_eq)
      CALL tblock_qp_dealloc(tb%qja_eq)
      CALL tblock_qp_dealloc(tb%qpres_n0)
      CALL tblock_qp_dealloc(tb%qnd_n0)

      RETURN
      END SUBROUTINE bl_tb_dealloc2
c-----------------------------------------------------------------------
c     subprogram 6. bl_seam_dealloc.
c     deallocates space for a single block seam.
c-----------------------------------------------------------------------
      SUBROUTINE bl_seam_dealloc(seam)
      USE local
      USE edge_type_mod
      IMPLICIT NONE

      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ivert,ierror

      ierror=0_i4
      DO ivert=1,seam%nvert
        DEALLOCATE(seam%vertex(ivert)%ptr)
        IF (seam%id==0) CYCLE
        DEALLOCATE(seam%vertex(ivert)%ptr2)
        DEALLOCATE(seam%vertex(ivert)%order)
        DEALLOCATE(seam%vertex(ivert)%seam_in)
        DEALLOCATE(seam%vertex(ivert)%seam_out)
        DEALLOCATE(seam%vertex(ivert)%seam_hold)
        DEALLOCATE(seam%vertex(ivert)%seam_save)
        DEALLOCATE(seam%vertex(ivert)%seam_cin)
        DEALLOCATE(seam%vertex(ivert)%seam_cout)
        DEALLOCATE(seam%vertex(ivert)%seam_chold)
        DEALLOCATE(seam%vertex(ivert)%seam_csave)
        IF (ASSOCIATED(seam%vertex(ivert)%tang))
     $    DEALLOCATE(seam%vertex(ivert)%tang)
        IF (ASSOCIATED(seam%vertex(ivert)%norm))
     $    DEALLOCATE(seam%vertex(ivert)%norm)
        DEALLOCATE(seam%segment(ivert)%seam_in)
        DEALLOCATE(seam%segment(ivert)%seam_out)
        DEALLOCATE(seam%segment(ivert)%seam_save)
        DEALLOCATE(seam%segment(ivert)%seam_cin)
        DEALLOCATE(seam%segment(ivert)%seam_cout)
        DEALLOCATE(seam%segment(ivert)%seam_csave)
        DEALLOCATE(seam%segment(ivert)%seam_mat_in)
        DEALLOCATE(seam%segment(ivert)%seam_mat_out)
        DEALLOCATE(seam%segment(ivert)%seam_mat_save)
        DEALLOCATE(seam%segment(ivert)%seam_mat_cin)
        DEALLOCATE(seam%segment(ivert)%seam_mat_cout)
        DEALLOCATE(seam%segment(ivert)%seam_mat_csave)
        IF (ASSOCIATED(seam%segment(ivert)%tang))
     $    DEALLOCATE(seam%segment(ivert)%tang)
        IF (ASSOCIATED(seam%segment(ivert)%norm))
     $    DEALLOCATE(seam%segment(ivert)%norm)
      ENDDO
      DEALLOCATE(seam%vertex)
      IF (seam%id==0) THEN
        IF (seam%nvert>0) DEALLOCATE(seam%excorner)
        RETURN
      ENDIF
      DEALLOCATE(seam%segment)
      DEALLOCATE(seam%expoint)
      DEALLOCATE(seam%r0point)

      RETURN
      END SUBROUTINE bl_seam_dealloc
