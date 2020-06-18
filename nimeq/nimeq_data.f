c-----------------------------------------------------------------------
c     file nimeq_data.f
c     contains dump file and data management for nimeq. it is based on
c     plot_data for nimplot, and that module name is retained for
c     other shared utilities.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  plot_data.
c     1.  acquire_data.
c-----------------------------------------------------------------------
c     subprogram 0. plot_data.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE plot_data
      USE local
      USE fields
      USE edge
      USE global
      IMPLICIT NONE

      REAL(r8) :: phi,dphi
      INTEGER(i4) :: mode_plot_start,mode_plot_end,iphi
      CHARACTER(1) :: flag,add_eq,do_flux='n'
      CHARACTER(16) :: selection_type
      CHARACTER(128) :: last_dump='initial_no_dump'
      LOGICAL :: read_file,inquire,inquire_flux

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. acquire_data.
c     acquires data from a dump file or from an existing data structure.
c-----------------------------------------------------------------------
      SUBROUTINE acquire_data(dump_data,dump_step,dump_time)
      USE input
      USE input_eq
      USE pardata
      USE rblock
      USE tblock
      USE surface
      USE dump
      USE nim_locate
      USE nimeq_all

      CHARACTER(*), INTENT(IN) :: dump_data
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time

      INTEGER(i4) :: ibl,imode,max_nqty,offm
      INTEGER(i4), SAVE :: im=0
      LOGICAL, SAVE :: paralloc=.false.,line_init=.false.,
     $                 slud_init=.false.
      LOGICAL :: new_dump
c-----------------------------------------------------------------------
c     read dump file if needed.  deallocate previous data if it exists.
c     option c corrupts the data.
c-----------------------------------------------------------------------
      dump_file=dump_data
      new_dump=.false.
      IF (read_file.AND.dump_file/=last_dump) THEN
        new_dump=.true.
        IF (last_dump/='initial_no_dump') CALL deallocate_data
        CALL dump_read(nmodes,nmodes_total,keff,keff_total,
     $                 dump_time,dump_step,dump_file)
        last_dump=dump_file
      ENDIF
c-----------------------------------------------------------------------
c     determine the total poloidal dimension for pseudo-spectral
c     coding and ffts.
c-----------------------------------------------------------------------
      IF (.NOT.paralloc) THEN
        ALLOCATE(mps_block(nbl))
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            mps_block(ibl)=rb(ibl)%mx*rb(ibl)%my
          ELSE
            mps_block(ibl)=tb(ibl)%mcell
          ENDIF
        ENDDO
        paralloc=.true.
      ENDIF
c-----------------------------------------------------------------------
c     initialize communication and block data, and
c     compute current density as a vertex quantity.
c-----------------------------------------------------------------------
      IF (new_dump) THEN
        CALL surface_set(ngr,poly_degree,integration_formula)
        DO ibl=1,nrbl
          CALL rblock_set(ngr,poly_degree,integration_formula,rb(ibl))
          CALL rblock_basis_set(rb(ibl),(/poly_degree/),
     $                          (/-1_i4/),(/-1_i4/),(/-1_i4/),
     $                          met_spl,geom)
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_set(tb(ibl))
          CALL tri_linear_get_areas(tb(ibl)%tgeom)
          CALL tblock_basis_set(tb(ibl),geom)
        ENDDO
        CALL variable_alloc(nmodes)
        CALL pointer_init
        max_nqty=MAX(9_i4,6_i4*nmodes)
        CALL edge_init(max_nqty,9_i4)
        IF (nprocs>1) CALL parallel_seam_init(max_nqty)
        CALL boundary_init(geom)
        max_nqty=MAX(2_i4,MAX(9_i4,6_i4*nmodes)*(poly_degree-1))
        offm=poly_degree**2+poly_degree
        CALL edge_segment_init(max_nqty,3_i4,offm)
        IF (nprocs>1) CALL parallel_seg_init(max_nqty,3_i4,offm)
        IF (solver=='gl_diaga'.AND..NOT.line_init) THEN
          CALL parallel_line_init(rb,nrbl,nbl,poly_degree)
          line_init=.true.
        ENDIF
        CALL block_create_tang(poly_degree)
c-----------------------------------------------------------------------
c       initialize the processor grid used by SuperLU_DIST.
c-----------------------------------------------------------------------
        IF (nimeq_solver(1:5)=='slu_d'.AND..NOT.slud_init) THEN
          slu_nrowp=SQRT(REAL(nprocs_layer))
          DO
            IF (MODULO(nprocs_layer,slu_nrowp)==0) EXIT
            slu_nrowp=slu_nrowp-1
          ENDDO
          slu_ncolp=nprocs_layer/slu_nrowp
          CALL c_fortran_slugrid(1_i4,comm_layer,node_layer,slu_nrowp,
     $                           slu_ncolp,slugrid_handle)
          slud_init=.true.
        ENDIF
        CALL mass_mat_init
      ENDIF
c-----------------------------------------------------------------------
c     initialize global data structures for searching.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL nimeq_all_init(poly_degree)
        CALL nim_rb_locate_init(nrbl_total,rball)
      ELSE
        CALL nim_rb_locate_init(nrbl,rb)
      ENDIF
c-----------------------------------------------------------------------
c     formats:
c-----------------------------------------------------------------------
  10  format(/'>>> ',a)
  20  format(2(/'>>> ',a))
  11  format('>>>? ',$)
 100  format(2(/'>>> ',a),/'>>> ',10a)
 110  format(/'>>> ',a,/'>>> ',10a,2(/'>>> ',a),/'>>>? ',$)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE acquire_data
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE plot_data
