c-----------------------------------------------------------------------
c     file dump.f
c     module containing all routines reading and writing data dumps.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module dump_write
c     1.  dump_write_mod.
c     2.  dump_write_rblock.
c     3.  dump_write_tblock.
c     4.  dump_write_seam.
c     5.  dump_write_bicube.
c     6.  dump_write_lagr_quad.
c     7.  dump_write_lagr_quad_2D.
c     8.  dump_write_tri_linear.
c     9.  dump_write_tri_linear_2D.
c     module dump_read
c     10.  dump_read_mod.
c     11. dump_read_rblock.
c     12. dump_read_tblock.
c     13. dump_read_seam.
c     14. dump_read_bicube.
c     15. dump_read_lagr_quad.
c     16. dump_read_lagr_quad_2D.
c     17. dump_read_tri_linear.
c     18. dump_read_tri_linear_2D.
c     module dump
c     19. dump_rblock_dealloc.
c     20. dump_tblock_dealloc.
c     21. dump_seam_dealloc.
c-----------------------------------------------------------------------
c     module dump_write_mod.
c     contains routines for writing dump files.
c-----------------------------------------------------------------------
      MODULE dump_write_mod
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      USE mpi_nim
      USE pardata
      USE time
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. dump_write.
c     writes dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write(nmodes,nmodes_total,keff_total,dump_time,
     $                      dump_step,dump_nm)
      
      INTEGER(i4), INTENT(IN) :: nmodes,nmodes_total,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff_total
      CHARACTER(*), INTENT(IN) :: dump_nm

      INTEGER(i4) :: ib,i,ierror
      INTEGER(i4) :: dump_code
      LOGICAL :: dump_test,unit_test
      TYPE(edge_type) :: seamtmp
      TYPE(rblock_type) :: rbtmp,rbtmp2
      TYPE(tblock_type) :: tbtmp,tbtmp2
      CHARACTER(64) :: step_name
      CHARACTER(128) :: dump_frd
c-----------------------------------------------------------------------
c     test for existing dump file.
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)
      IF (node == 0) THEN
        dump_code = 0
        IF (dump_step>999999) THEN
          WRITE(step_name,fmt='(i7.7)') dump_step
        ELSE IF (dump_step>99999) THEN
          WRITE(step_name,fmt='(i6.6)') dump_step
        ELSE
          WRITE(step_name,fmt='(i5.5)') dump_step
        ENDIF
        IF(dump_dir.EQ.".")THEN
          dump_frd=TRIM(dump_nm)//"."//TRIM(step_name)
        ELSE
          dump_frd=TRIM(dump_dir)//"/"//TRIM(dump_nm)//"."
     $         //TRIM(step_name)
        ENDIF
        INQUIRE(file=TRIM(dump_frd),exist=dump_test)
c-----------------------------------------------------------------------
c       open dump file in appropriate mode.
c-----------------------------------------------------------------------
        IF(dump_test)THEN
          SELECT CASE (dump_over)
          CASE(0)
            CALL open_bin
     $           (dump_unit,TRIM(dump_frd),'REPLACE','ASIS',64_i4)
          CASE(1)
            CALL open_bin
     $           (dump_unit,TRIM(dump_frd),'OLD','APPEND',64_i4)
          CASE(2)
            dump_code = 1
            WRITE(nim_wr,*)'Unable to open dump file, '//TRIM(dump_frd)
            RETURN
          END SELECT
        ELSE
          CALL open_bin
     $         (dump_unit,TRIM(dump_frd),'NEW','ASIS',64_i4)
        ENDIF
        WRITE(nim_wr,'(/,2a)') 'Writing to restart file ',
     $                            TRIM(dump_frd)
        IF (.NOT.out_opened) THEN
          OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $         POSITION='APPEND')
          out_opened=.true.
        ENDIF
        WRITE(out_unit,'(/,2a)') 'Writing to restart file ',
     $                              TRIM(dump_frd)
      ENDIF
c-----------------------------------------------------------------------
c     check dump_code flag
c-----------------------------------------------------------------------
      IF (nprocs > 1) CALL mpi_bcast(dump_code,1,mpi_nim_int,0,
     $     mpi_comm_world,ierror)
      IF (dump_code == 1) RETURN
c-----------------------------------------------------------------------
c     write global data, converting integers to 64 bit reals to have
c     the same # of bits as the reals and to avoid 64 bit integers
c     for linux machines.
c-----------------------------------------------------------------------
      IF (node==0) THEN
        WRITE(dump_unit) dump_time
        WRITE(dump_unit) REAL(dump_step,r8)
        WRITE(dump_unit) REAL(nbl_total,r8)
        WRITE(dump_unit) REAL(nrbl_total,r8)
        WRITE(dump_unit) REAL(poly_degree,r8)
        WRITE(dump_unit) REAL(nmodes_total,r8)
        WRITE(dump_unit) keff_total
      ENDIF
c-----------------------------------------------------------------------
c     write seam data.  all written by node 0
c-----------------------------------------------------------------------
      IF (node == 0) CALL dump_write_seam(seam0)
      IF (nprocs==1) THEN
        DO ib=1,nbl_total
          CALL dump_write_seam(seam(ib))
        ENDDO
      ELSE IF (ilayer==0) THEN
        DO ib=1,nbl_total
          IF (block2proc(ib)==0) THEN
            IF (node==0) CALL dump_write_seam(seam(global2local(ib)))
          ELSE
            IF (node==0.OR.node==block2proc(ib))
     $        CALL send_seam(block2proc(ib),seam(global2local(ib)),
     $                       0,seamtmp)
            IF (node==0) THEN
              CALL dump_write_seam(seamtmp)
              CALL dump_seam_dealloc(seamtmp)
            ENDIF
          ENDIF
          CALL mpi_barrier(comm_layer,ierror)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write rblock data.  all written by node 0, and all Fourier layers
c     are temporarily packed into layer 0.
c-----------------------------------------------------------------------
      IF (nprocs == 1) THEN
        DO ib=1,nrbl_total
          CALL dump_write_rblock(rb(ib))
        ENDDO
      ELSE
        DO ib=1,nrbl_total
          IF (block2proc(ib)==node)
     $      CALL gather_rblock(rbtmp2,rb(global2local(ib)),
     $                         nmodes,nmodes_total)
          IF (block2proc(ib)==0) THEN
            IF (node==0) CALL dump_write_rblock(rbtmp2)
          ELSE IF (ilayer==0) THEN
            IF (node==0.OR.node==block2proc(ib))
     $        CALL send_rblock(block2proc(ib),rbtmp2,0,rbtmp)
            IF (node==0) THEN
              CALL dump_write_rblock(rbtmp)
              CALL dump_rblock_dealloc(rbtmp)
            ENDIF
          ENDIF
          IF (block2proc(ib)==node.AND.ilayer==0)
     $      CALL dump_rblock_dealloc(rbtmp2)
        ENDDO
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     write tblock data.  all written by node 0
c-----------------------------------------------------------------------
      IF (nprocs == 1) THEN
        DO ib=nrbl_total+1,nbl_total
          CALL dump_write_tblock(tb(ib))
        ENDDO
      ELSE
        DO ib=nrbl_total+1,nbl_total
          IF (block2proc(ib)==node)
     $      CALL gather_tblock(tbtmp2,tb(global2local(ib)),
     $                         nmodes,nmodes_total)
          IF (block2proc(ib)==0) THEN
            IF (node==0) CALL dump_write_tblock(tbtmp2)
          ELSE IF (ilayer==0) THEN
            IF (node==0.OR.node==block2proc(ib))
     $        CALL send_tblock(block2proc(ib),tbtmp2,0,tbtmp)
            IF (node==0) THEN
              CALL dump_write_tblock(tbtmp)
              CALL dump_tblock_dealloc(tbtmp)
            ENDIF
          ENDIF
          IF (block2proc(ib)==node.AND.ilayer==0)
     $      CALL dump_tblock_dealloc(tbtmp2)
        ENDDO
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      IF (node == 0) CALL close_bin(dump_unit,TRIM(dump_frd))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
      RETURN
      END  SUBROUTINE dump_write
c-----------------------------------------------------------------------
c     subprogram 2. dump_write_rblock.
c     writes rblock data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_rblock(rb)

      TYPE (rblock_type), INTENT(IN) :: rb         
      INTEGER(i4) :: ix,iy,iqty
c-----------------------------------------------------------------------
c     write block descriptors.
c-----------------------------------------------------------------------
      WRITE(dump_unit) REAL(rb%id,r8)
      WRITE(dump_unit) REAL(rb%mx,r8)
      WRITE(dump_unit) REAL(rb%my,r8)
      IF (rb%degenerate) THEN
         WRITE(dump_unit) 1._r8
      ELSE
         WRITE(dump_unit)-1._r8
      ENDIF
c-----------------------------------------------------------------------
c     write coordinates.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%rz)
c-----------------------------------------------------------------------
c     write magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%be_eq)
      CALL dump_write_lagr_quad(rb%be)
c-----------------------------------------------------------------------
c     write equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%ja_eq)
c-----------------------------------------------------------------------
c     write fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%ve_eq)
      CALL dump_write_lagr_quad(rb%ve)
c-----------------------------------------------------------------------
c     write pressures.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%pres_eq)
      CALL dump_write_lagr_quad(rb%pres)
      CALL dump_write_lagr_quad_2D(rb%prese_eq)
      CALL dump_write_lagr_quad(rb%prese)
c-----------------------------------------------------------------------
c     write number densities.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%nd_eq)
      CALL dump_write_lagr_quad(rb%nd)
c-----------------------------------------------------------------------
c     write diffusivity shape function and material concentration.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(rb%diff_shape)
      CALL dump_write_lagr_quad(rb%conc)
c-----------------------------------------------------------------------
c     write temperatures.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad(rb%tele)
      CALL dump_write_lagr_quad(rb%tion)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_write_rblock
c-----------------------------------------------------------------------
c     subprogram 3. dump_write_tblock.
c     writes tblock data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_tblock(tb)

      TYPE (tblock_type), INTENT(IN) :: tb         
      INTEGER(i4) :: ivert,icell,iqty,inbr,nnbr
c-----------------------------------------------------------------------
c     write geometry.
c-----------------------------------------------------------------------
      WRITE(dump_unit) REAL(tb%id,r8)
      WRITE(dump_unit) REAL(tb%tgeom%mvert,r8)
      WRITE(dump_unit) REAL(tb%tgeom%mcell,r8)
      WRITE(dump_unit)tb%tgeom%xs
      WRITE(dump_unit)tb%tgeom%ys
      WRITE(dump_unit) REAL(tb%tgeom%vertex,r8)
      DO ivert=0,tb%tgeom%mvert
         nnbr=SIZE(tb%tgeom%neighbor(ivert)%vertex)-1
         WRITE(dump_unit) REAL(nnbr,r8)
         WRITE(dump_unit) REAL(tb%tgeom%neighbor(ivert)%vertex,r8)
      ENDDO
c-----------------------------------------------------------------------
c     write magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%be_eq)
      CALL dump_write_tri_linear(tb%be)
c-----------------------------------------------------------------------
c     write equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%ja_eq)
c-----------------------------------------------------------------------
c     write fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%ve_eq)
      CALL dump_write_tri_linear(tb%ve)
c-----------------------------------------------------------------------
c     write pressures.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%pres_eq)
      CALL dump_write_tri_linear(tb%pres)
      CALL dump_write_tri_linear_2D(tb%prese_eq)
      CALL dump_write_tri_linear(tb%prese)
c-----------------------------------------------------------------------
c     write number densities.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%nd_eq)
      CALL dump_write_tri_linear(tb%nd)
c-----------------------------------------------------------------------
c     write diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear_2D(tb%diff_shape)
      CALL dump_write_tri_linear(tb%conc)
c-----------------------------------------------------------------------
c     write temperatures.
c-----------------------------------------------------------------------
      CALL dump_write_tri_linear(tb%tele)
      CALL dump_write_tri_linear(tb%tion)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_write_tblock
c-----------------------------------------------------------------------
c     subprogram 4. dump_write_seam.
c     writes seam data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_seam(seam)

      TYPE (edge_type), INTENT(IN) :: seam

      INTEGER(i4) :: i,nv,np,ip,iq
c-----------------------------------------------------------------------
c     write block descriptors.
c-----------------------------------------------------------------------
      WRITE(dump_unit) REAL(seam%id,r8)
      WRITE(dump_unit) REAL(seam%nvert,r8)
c-----------------------------------------------------------------------
c     write seam data.
c-----------------------------------------------------------------------
      DO i=1,seam%nvert
         np=SIZE(seam%vertex(i)%ptr,2)
         WRITE(dump_unit) REAL(np,r8)
         WRITE(dump_unit) REAL(seam%vertex(i)%ptr,r8)
         IF (seam%id>0) WRITE(dump_unit) REAL(seam%vertex(i)%intxy,r8)
      ENDDO
c-----------------------------------------------------------------------
c     write external corner data.
c-----------------------------------------------------------------------
      IF(seam%id==0.AND.seam%nvert>0)THEN
         DO i=1,seam%nvert
            IF (seam%excorner(i)) THEN
               WRITE(dump_unit) 1._r8
            ELSE
               WRITE(dump_unit) 0._r8
            ENDIF
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE dump_write_seam
c-----------------------------------------------------------------------
c     subprogram 5. dump_write_bicube.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_bicube(rb_bc)

      TYPE(bicube_type), INTENT(IN) :: rb_bc

      WRITE(dump_unit) REAL(rb_bc%nqty,r8)
      WRITE(dump_unit) rb_bc%xs
      WRITE(dump_unit) rb_bc%ys
      WRITE(dump_unit) rb_bc%fs
      WRITE(dump_unit) rb_bc%fsx
      WRITE(dump_unit) rb_bc%fsy
      WRITE(dump_unit) rb_bc%fsxy
     
      RETURN
      END SUBROUTINE dump_write_bicube
c-----------------------------------------------------------------------
c     subprogram 6.  dump_write_lagr_quad.
c     real and imaginary parts are witten separately for possible
c     data conversion.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_lagr_quad(laq)

      TYPE(lagr_quad_type), INTENT(IN) :: laq

      WRITE(dump_unit) REAL(laq%nqty,r8)
      WRITE(dump_unit) REAL(laq%nfour,r8)
      WRITE(dump_unit) REAL(laq%n_side,r8)
      WRITE(dump_unit) REAL(laq%n_int,r8)
      WRITE(dump_unit) REAL(laq%fs,r8)
      WRITE(dump_unit) REAL(-(0,1)*laq%fs,r8)
      IF (ALLOCATED(laq%fsh)) THEN
        WRITE(dump_unit) REAL(laq%fsh,r8)
        WRITE(dump_unit) REAL(-(0,1)*laq%fsh,r8)
        WRITE(dump_unit) REAL(laq%fsv,r8)
        WRITE(dump_unit) REAL(-(0,1)*laq%fsv,r8)
        WRITE(dump_unit) REAL(laq%fsi,r8)
        WRITE(dump_unit) REAL(-(0,1)*laq%fsi,r8)
      ENDIF
     
      RETURN
      END SUBROUTINE dump_write_lagr_quad
c-----------------------------------------------------------------------
c     subprogram 7.  dump_write_lagr_quad_2D.
c     2D version.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_lagr_quad_2D(laq)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq

      WRITE(dump_unit) REAL(laq%nqty,r8)
      WRITE(dump_unit) REAL(laq%n_side,r8)
      WRITE(dump_unit) REAL(laq%n_int,r8)
      WRITE(dump_unit) REAL(laq%fs,r8)
      IF (ALLOCATED(laq%fsh)) THEN
        WRITE(dump_unit) REAL(laq%fsh,r8)
        WRITE(dump_unit) REAL(laq%fsv,r8)
        WRITE(dump_unit) REAL(laq%fsi,r8)
      ENDIF

      RETURN
      END SUBROUTINE dump_write_lagr_quad_2D
c-----------------------------------------------------------------------
c     subprogram 8. dump_write_tri_linear.
c     real and imaginary parts are witten separately for possible
c     data conversion.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear(tb_l)

      TYPE(tri_linear_type), INTENT(IN) :: tb_l

      WRITE(dump_unit) REAL(tb_l%nqty,r8)
      WRITE(dump_unit) REAL(tb_l%nfour,r8)
      WRITE(dump_unit) REAL(tb_l%fs,r8)
      WRITE(dump_unit) REAL(-(0,1)*tb_l%fs,r8)
     
      RETURN
      END SUBROUTINE dump_write_tri_linear
c-----------------------------------------------------------------------
c     subprogram 9. dump_write_tri_linear_2D.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear_2D(tb_l)

      TYPE(tri_linear_2D_type), INTENT(IN) :: tb_l

      WRITE(dump_unit) REAL(tb_l%nqty,r8)
      WRITE(dump_unit) tb_l%fs
     
      RETURN
      END SUBROUTINE dump_write_tri_linear_2D
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE dump_write_mod
c-----------------------------------------------------------------------
c     module dump_read_mod.
c     contains routines for writing dump files.
c-----------------------------------------------------------------------
      MODULE dump_read_mod
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      USE mpi_nim
      USE pardata
      USE time
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 10. dump_read.
c     reads dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read(nmodes,nmodes_total,keff,keff_total,
     $                     dump_time,dump_step,dump_frd,rname)

      INTEGER(i4), INTENT(OUT) :: dump_step,nmodes
      INTEGER(i4), INTENT(IN) :: nmodes_total
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff,keff_total
      CHARACTER(*), INTENT(IN) :: dump_frd
      CHARACTER(*), INTENT(IN), OPTIONAL :: rname

      INTEGER(i4) :: ib,nb,dump_code=0_i4,i,ierror,nmodes_dump,prev,
     $               poly_deg_dump
      TYPE(edge_type) :: seamtmp
      TYPE(rblock_type) :: rbtmp
      TYPE(tblock_type) :: tbtmp
      REAL(r8) :: iread
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timestart)
      IF (node == 0) THEN
        INQUIRE(FILE=TRIM(dump_frd),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop
     $    ('Dump file '//TRIM(dump_frd)//' does not exist.')
        CALL open_bin(rstrt_unit,TRIM(dump_frd),'OLD','REWIND',64_i4)
        WRITE(nim_wr,'(/,2a)') 'Reading from restart file ',
     $                          TRIM(dump_frd)
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        WRITE(out_unit,'(/,2a)') 'Reading from restart file ',
     $                            TRIM(dump_frd)
        CLOSE(UNIT=out_unit)
        out_opened=.false.
        READ(rstrt_unit) dump_time
        READ(rstrt_unit) iread
        dump_step=NINT(iread)
        READ(rstrt_unit) iread
        nbl_total=NINT(iread)
        READ(rstrt_unit) iread
        nrbl_total=NINT(iread)
        READ(rstrt_unit) iread
        poly_deg_dump=NINT(iread)
        READ(rstrt_unit) iread
        nmodes_dump=NINT(iread)
      ENDIF
c-----------------------------------------------------------------------
c     broadcast global information.
c-----------------------------------------------------------------------
      IF (nprocs > 1) THEN
        CALL mpi_bcast(dump_time,1,mpi_nim_real,0,
     $       mpi_comm_world,ierror)
        CALL mpi_bcast(dump_step,1,mpi_nim_int,0,mpi_comm_world,ierror)
        CALL mpi_bcast(nbl_total,1,mpi_nim_int,0,mpi_comm_world,ierror)
        CALL mpi_bcast(nrbl_total,1,mpi_nim_int,0,mpi_comm_world,ierror)
        CALL mpi_bcast(poly_deg_dump,1,mpi_nim_int,0,
     $       mpi_comm_world,ierror)
        CALL mpi_bcast(nmodes_dump,1,mpi_nim_int,0,
     $       mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     check for basis consistency.
c-----------------------------------------------------------------------
      IF (poly_deg_dump/=poly_degree) THEN
         CALL nim_stop
     $     ('NIMROD input is not consistent with dumped poly_degree.')
      ENDIF
c-----------------------------------------------------------------------
c     check for fourier rep. consistency, then read wavenumber array.
c-----------------------------------------------------------------------
      IF (nmodes_dump/=nmodes_total) CALL nim_stop
     $  ('NIMROD input is not consistent with dumped nmodes.')
c-----------------------------------------------------------------------
c     assign block(s) and Fourier components to each processor--set
c     local nbl, nrbl, and nmodes for each proc.
c-----------------------------------------------------------------------
      IF (PRESENT(rname)) THEN
        CALL parallel_block_init(nmodes,nmodes_total,nlayers,decompflag,
     $                           rname)
      ELSE
        CALL parallel_block_init(nmodes,nmodes_total,nlayers,decompflag,
     $                           'plasma')
      ENDIF
c-----------------------------------------------------------------------
c     now read and broadcast keff values.
c-----------------------------------------------------------------------
      ALLOCATE(keff(nmodes),keff_total(nmodes_total))
      IF (node==0) READ(rstrt_unit) keff_total
      IF (nprocs>1) CALL mpi_bcast
     $  (keff_total,nmodes_total,mpi_nim_real,0,mpi_comm_world,ierror)
      keff=keff_total(mode_lo:mode_hi)
c-----------------------------------------------------------------------
c     read seam data.  all read by node 0, then distributed across
c     layer 0.
c-----------------------------------------------------------------------
      ALLOCATE(seam0)
      ALLOCATE(seam(nbl))
      IF (nprocs==1) THEN
        CALL dump_read_seam(seam0)
      ELSE
        IF (node==0) CALL dump_read_seam(seam0)
        DO i=1,nprocs-1
          CALL send_seam(0,seam0,i,seam0)
        ENDDO
      ENDIF
      IF (nprocs==1) THEN
        DO ib=1,nbl_total
          CALL dump_read_seam(seam(ib))
        ENDDO
      ELSE IF (ilayer==0) THEN
        DO ib=1,nbl_total
          IF (block2proc(ib)==0) THEN
            IF (node==0) CALL dump_read_seam(seam(global2local(ib)))
          ELSE
            IF (node==0) CALL dump_read_seam(seamtmp)
            IF (node==0.OR.node==block2proc(ib)) THEN
              CALL send_seam(0,seamtmp,
     $             block2proc(ib),seam(global2local(ib)))
            ENDIF
            IF (node==0) CALL dump_seam_dealloc(seamtmp)
          ENDIF
        ENDDO
        CALL mpi_barrier(comm_layer,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     distribute seams to all layers.
c-----------------------------------------------------------------------
      IF (nlayers>1) THEN
        DO ib=1,nbl
          DO i=1,nlayers-1
            CALL send_seam(layer2proc(0),seam(ib),layer2proc(i),
     $                     seam(ib))
          ENDDO
        ENDDO
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     read rblock data.  all read by node 0
c-----------------------------------------------------------------------
      ALLOCATE(rb(nrbl))
      IF (nprocs==1) THEN
        DO ib=1,nrbl_total
          CALL dump_read_rblock(rb(ib))
        ENDDO
      ELSE IF (ilayer==0) THEN
        DO ib=1,nrbl_total
          IF (block2proc(ib)==0) THEN
            IF (node == 0) CALL dump_read_rblock(rb(global2local(ib)))
          ELSE
            IF (node == 0) CALL dump_read_rblock(rbtmp)
            IF (node==0.OR.node==block2proc(ib)) THEN
              CALL send_rblock(0,rbtmp,
     $             block2proc(ib),rb(global2local(ib)))
            ENDIF
            IF (node == 0) CALL dump_rblock_dealloc(rbtmp)
          ENDIF
        ENDDO
        CALL mpi_barrier(comm_layer,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     distribute rblocks to all layers, then reduce to the local modes.
c-----------------------------------------------------------------------
      IF (nlayers>1) THEN
        DO ib=1,nrbl
          DO i=1,nlayers-1
            CALL send_rblock(layer2proc(0),rb(ib),layer2proc(i),rb(ib))
          ENDDO
          CALL trim_rblock(rb(ib),nmodes,nmodes_total)
        ENDDO
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     read tblock data.  all read by node 0
c-----------------------------------------------------------------------
      ALLOCATE(tb(nrbl+1:nbl))
      IF (nprocs==1) THEN
        DO ib=nrbl_total+1,nbl_total
          CALL dump_read_tblock(tb(ib))
        ENDDO
      ELSE IF (ilayer==0) THEN
        DO ib=nrbl_total+1,nbl_total
          IF (block2proc(ib)==0) THEN
            IF (node == 0) CALL dump_read_tblock(tb(global2local(ib)))
          ELSE
            IF (node == 0) CALL dump_read_tblock(tbtmp)
            IF (node==0.OR.node==block2proc(ib)) THEN
              CALL send_tblock(0,tbtmp,
     $             block2proc(ib),tb(global2local(ib)))
            ENDIF
            IF (node == 0) CALL dump_tblock_dealloc(tbtmp)
          ENDIF
        ENDDO
        CALL mpi_barrier(comm_layer,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     distribute tblocks to all layers, then reduce to the local modes.
c-----------------------------------------------------------------------
      IF (nlayers>1) THEN
        DO ib=nrbl+1,nbl
          DO i=1,nlayers-1
            CALL send_tblock(layer2proc(0),tb(ib),layer2proc(i),tb(ib))
          ENDDO
          CALL trim_tblock(tb(ib),nmodes,nmodes_total)
        ENDDO
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDIF
c-----------------------------------------------------------------------
c     now that block sizes are know, assemble block_sizes vector that
c     stores global sizes of all blocks - copy of it on every proc
c-----------------------------------------------------------------------
      CALL parallel_block_init2
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      IF (node==0) CALL close_bin(rstrt_unit,TRIM(dump_frd))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(timeend)
      time_io = time_io + timeend-timestart
      RETURN
      END SUBROUTINE dump_read
c-----------------------------------------------------------------------
c     subprogram 11. dump_read_rblock.
c     reads rblock data from dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_rblock(rb)

      TYPE (rblock_type), INTENT(OUT) :: rb         
      INTEGER(i4) :: ix,iy,iqty
      REAL(r8) :: iread
      REAL(r8) :: rread
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      rb%name='rblock'
      READ(rstrt_unit) iread
      rb%id=NINT(iread)
      READ(rstrt_unit) iread
      rb%mx=NINT(iread)
      READ(rstrt_unit) iread
      rb%my=NINT(iread)
      READ(rstrt_unit) rread
      IF (rread>0) THEN
        rb%degenerate=.true.
      ELSE
        rb%degenerate=.false.
      ENDIF
c-----------------------------------------------------------------------
c     read coordinates.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%rz,rb%mx,rb%my,'rz',
     $                            (/'   r  ','   z  '/))
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%be_eq,rb%mx,rb%my,'bq',
     $                            (/' be_eq'/))
      CALL dump_read_lagr_quad(rb%be,rb%mx,rb%my,'be',(/'  be  '/))
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ja_eq,rb%mx,rb%my,'jq',
     $                            (/' ja_eq'/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ve_eq,rb%mx,rb%my,'vq',
     $                            (/' ve_eq'/))
      CALL dump_read_lagr_quad(rb%ve,rb%mx,rb%my,'ve',(/'  ve  '/))
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%pres_eq,rb%mx,rb%my,'pq',
     $                            (/' pr_eq'/))
      CALL dump_read_lagr_quad(rb%pres,rb%mx,rb%my,'pr',(/' pres '/))
      CALL dump_read_lagr_quad_2D(rb%prese_eq,rb%mx,rb%my,'pq',
     $                            (/'pre_eq'/))
      CALL dump_read_lagr_quad(rb%prese,rb%mx,rb%my,'pe',(/' prese'/))
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%nd_eq,rb%mx,rb%my,'nq',
     $                            (/' nd_eq'/))
      CALL dump_read_lagr_quad(rb%nd,rb%mx,rb%my,'nd',(/'  nd  '/))
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%diff_shape,rb%mx,rb%my,'ds',
     $                            (/'dif_sh'/))
      CALL dump_read_lagr_quad(rb%conc,rb%mx,rb%my,'co',(/' conc '/))
c-----------------------------------------------------------------------
c     read temperatures.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad(rb%tele,rb%mx,rb%my,'te',(/' tele '/))
      CALL dump_read_lagr_quad(rb%tion,rb%mx,rb%my,'ti',(/' tion '/))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_rblock
c-----------------------------------------------------------------------
c     subprogram 12. dump_read_tblock.
c     reads tblock data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tblock(tb)

      TYPE (tblock_type), INTENT(OUT) :: tb
      INTEGER(i4) :: ivert,icell,iqty,inbr,nnbr
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ia1read
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ia2read
c-----------------------------------------------------------------------
c     read geometry.
c-----------------------------------------------------------------------
      tb%name='tblock'
      READ(rstrt_unit) iread
      tb%id=NINT(iread)
      READ(rstrt_unit) iread
      tb%tgeom%mvert=NINT(iread)
      READ(rstrt_unit) iread
      tb%tgeom%mcell=NINT(iread)
      CALL tri_linear_geom_alloc(tb%tgeom,tb%tgeom%mvert,tb%tgeom%mcell)
      READ(rstrt_unit) tb%tgeom%xs
      READ(rstrt_unit) tb%tgeom%ys
      ALLOCATE(ia2read(tb%tgeom%mcell,3))
      READ(rstrt_unit) ia2read
      tb%tgeom%vertex=NINT(ia2read)
      DEALLOCATE(ia2read)
      DO ivert=0,tb%tgeom%mvert
         READ(rstrt_unit) iread
         nnbr=NINT(iread)
         ALLOCATE(tb%tgeom%neighbor(ivert)%vertex(0:nnbr))
         ALLOCATE(ia1read(0:nnbr))
         READ(rstrt_unit) ia1read
         tb%tgeom%neighbor(ivert)%vertex=NINT(ia1read)
         DEALLOCATE(ia1read)
      ENDDO
      tb%mvert=tb%tgeom%mvert
      tb%mcell=tb%tgeom%mcell
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%be_eq,tb%mvert,'bq',(/' be_eq'/))
      CALL dump_read_tri_linear(tb%be,tb%mvert,'be',(/'  be  '/))
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ja_eq,tb%mvert,'jq',(/' ja_eq'/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ve_eq,tb%mvert,'vq',(/' ve_eq'/))
      CALL dump_read_tri_linear(tb%ve,tb%mvert,'ve',(/'  ve  '/))
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%pres_eq,tb%mvert,'pq',
     $                             (/' pr_eq'/))
      CALL dump_read_tri_linear(tb%pres,tb%mvert,'pr',(/' pres '/))
      CALL dump_read_tri_linear_2D(tb%prese_eq,tb%mvert,'pq',
     $                             (/'pre_eq'/))
      CALL dump_read_tri_linear(tb%prese,tb%mvert,'pe',(/' prese'/))
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%nd_eq,tb%mvert,'nq',(/' nd_eq'/))
      CALL dump_read_tri_linear(tb%nd,tb%mvert,'nd',(/'  nd  '/))
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%diff_shape,tb%mvert,'ds',
     $                          (/'dif_sh'/))
      CALL dump_read_tri_linear(tb%conc,tb%mvert,'co',(/' conc '/))
c-----------------------------------------------------------------------
c     read temperatures.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear(tb%tele,tb%mvert,'te',(/' tele '/))
      CALL dump_read_tri_linear(tb%tion,tb%mvert,'ti',(/' tion '/))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_tblock
c-----------------------------------------------------------------------
c     subprogram 13. dump_read_seam.
c     reads seam data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_seam(seam)

      TYPE (edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: i,nv,np,ip,iq
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ia1read
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ia2read
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      seam%name='seam'
      READ(rstrt_unit) iread
      seam%id=NINT(iread)
      READ(rstrt_unit) iread
      seam%nvert=NINT(iread)
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      ALLOCATE(seam%vertex(seam%nvert))
      ALLOCATE(ia1read(2))
      DO i=1,seam%nvert
         READ(rstrt_unit) iread
         np=NINT(iread)
         ALLOCATE(seam%vertex(i)%ptr(2,np))
         ALLOCATE(ia2read(2,np))
         READ(rstrt_unit) ia2read
         seam%vertex(i)%ptr=NINT(ia2read)
         DEALLOCATE(ia2read)
         IF (seam%id>0) THEN
           READ(rstrt_unit) ia1read
           seam%vertex(i)%intxy=NINT(ia1read)
         ENDIF
      ENDDO
      DEALLOCATE(ia1read)
c-----------------------------------------------------------------------
c     read external corner data.
c-----------------------------------------------------------------------
      IF(seam%id==0.AND.seam%nvert>0)THEN
         ALLOCATE(seam%excorner(seam%nvert))
         DO i=1,seam%nvert
            READ(rstrt_unit) iread
            seam%excorner(i)=(NINT(iread)>0)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE dump_read_seam
c-----------------------------------------------------------------------
c     subprogram 14. dump_read_bicube.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_bicube(rb_bc,lx,ly,name,title)

      INTEGER(i4), INTENT(IN) :: lx,ly
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(bicube_type), INTENT(OUT) :: rb_bc

      INTEGER(i4) :: nqloc
      REAL(r8) :: iread

      READ(rstrt_unit) iread
      nqloc=NINT(iread)
      CALL bicube_alloc(rb_bc,lx,ly,nqloc)
      rb_bc%name=name
      IF (SIZE(title)<SIZE(rb_bc%title)) THEN
        rb_bc%title=title(1)
      ELSE
        rb_bc%title=title
      ENDIF
      READ(rstrt_unit) rb_bc%xs
      READ(rstrt_unit) rb_bc%ys
      READ(rstrt_unit) rb_bc%fs
      READ(rstrt_unit) rb_bc%fsx
      READ(rstrt_unit) rb_bc%fsy
      READ(rstrt_unit) rb_bc%fsxy

      RETURN
      END SUBROUTINE dump_read_bicube
c-----------------------------------------------------------------------
c     subprogram 15. dump_read_lagr_quad.
c     real and imaginary parts are read separately for possible
c     data conversion.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_quad(laq,lx,ly,name,title)

      INTEGER(i4), INTENT(IN) :: lx,ly
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(lagr_quad_type), INTENT(OUT) :: laq

      INTEGER(i4) :: nqloc,nsloc,niloc,nfloc
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rread2

      READ(rstrt_unit) iread
      nqloc=NINT(iread)
      READ(rstrt_unit) iread
      nfloc=NINT(iread)
      READ(rstrt_unit) iread
      nsloc=NINT(iread)
      READ(rstrt_unit) iread
      niloc=NINT(iread)
      CALL lagr_quad_alloc(laq,lx,ly,nqloc,nfloc,nsloc+1_i4)
      laq%name=name
      IF (SIZE(title)<SIZE(laq%title)) THEN
        laq%title=title(1)
      ELSE
        laq%title=title
      ENDIF

      ALLOCATE(rread(laq%nqty,0:lx,0:ly,laq%nfour))
      READ(rstrt_unit) rread
      laq%fs=rread
      READ(rstrt_unit) rread
      laq%fs=laq%fs+(0,1)*rread
      DEALLOCATE(rread)

      IF (laq%n_side>0) THEN
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx,0:ly,laq%nfour))
        READ(rstrt_unit) rread2
        laq%fsh=rread2
        READ(rstrt_unit) rread2
        laq%fsh=laq%fsh+(0,1)*rread2
        DEALLOCATE(rread2)
        ALLOCATE(rread2(laq%nqty,laq%n_side,0:lx,1:ly,laq%nfour))
        READ(rstrt_unit) rread2
        laq%fsv=rread2
        READ(rstrt_unit) rread2
        laq%fsv=laq%fsv+(0,1)*rread2
        DEALLOCATE(rread2)
        ALLOCATE(rread2(laq%nqty,laq%n_int,1:lx,1:ly,laq%nfour))
        READ(rstrt_unit) rread2
        laq%fsi=rread2
        READ(rstrt_unit) rread2
        laq%fsi=laq%fsi+(0,1)*rread2
        DEALLOCATE(rread2)
      ENDIF

      RETURN
      END SUBROUTINE dump_read_lagr_quad
c-----------------------------------------------------------------------
c     subprogram 16. dump_read_lagr_quad_2D.
c     2D version.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_quad_2D(laq,lx,ly,name,title)

      INTEGER(i4), INTENT(IN) :: lx,ly
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(lagr_quad_2D_type), INTENT(OUT) :: laq

      INTEGER(i4) :: nqloc,nsloc,niloc
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread2

      READ(rstrt_unit) iread
      nqloc=NINT(iread)
      READ(rstrt_unit) iread
      nsloc=NINT(iread)
      READ(rstrt_unit) iread
      niloc=NINT(iread)

      CALL lagr_quad_alloc(laq,lx,ly,nqloc,nsloc+1_i4)
      laq%name=name
      IF (SIZE(title)<SIZE(laq%title)) THEN
        laq%title=title(1)
      ELSE
        laq%title=title
      ENDIF

      ALLOCATE(rread(laq%nqty,0:lx,0:ly))
      READ(rstrt_unit) rread
      laq%fs=rread
      DEALLOCATE(rread)

      IF (laq%n_side>0) THEN
        ALLOCATE(rread2(laq%nqty,laq%n_side,1:lx,0:ly))
        READ(rstrt_unit) rread2
        laq%fsh=rread2
        DEALLOCATE(rread2)
        ALLOCATE(rread2(laq%nqty,laq%n_side,0:lx,1:ly))
        READ(rstrt_unit) rread2
        laq%fsv=rread2
        DEALLOCATE(rread2)
        ALLOCATE(rread2(laq%nqty,laq%n_int,1:lx,1:ly))
        READ(rstrt_unit) rread2
        laq%fsi=rread2
        DEALLOCATE(rread2)
      ENDIF

      RETURN
      END SUBROUTINE dump_read_lagr_quad_2D
c-----------------------------------------------------------------------
c     subprogram 17. dump_read_tri_linear.
c     real and imaginary parts are read separately for possible
c     data conversion.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear(tb_l,lv,name,title)

      INTEGER(i4), INTENT(IN) :: lv
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_type), INTENT(OUT) :: tb_l

      INTEGER(i4) :: nqloc,nfloc
      REAL(r8) :: iread
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread

      READ(rstrt_unit) iread
      nqloc=NINT(iread)
      READ(rstrt_unit) iread
      nfloc=NINT(iread)
      CALL tri_linear_alloc(tb_l,lv,nqloc,nfloc)
      tb_l%name=name
      IF (SIZE(title)<SIZE(tb_l%title)) THEN
        tb_l%title=title(1)
      ELSE
        tb_l%title=title
      ENDIF

      ALLOCATE(rread(tb_l%nqty,0:lv,0:0,tb_l%nfour))
      READ(rstrt_unit) rread
      tb_l%fs=rread
      READ(rstrt_unit) rread
      tb_l%fs=tb_l%fs+(0,1)*rread
      DEALLOCATE(rread)
     
      RETURN
      END SUBROUTINE dump_read_tri_linear
c-----------------------------------------------------------------------
c     subprogram 18. dump_read_tri_linear_2D.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear_2D(tb_l,lv,name,title)

      INTEGER(i4), INTENT(IN) :: lv
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_2D_type), INTENT(OUT) :: tb_l

      INTEGER(i4) :: nqloc
      REAL(r8) :: iread

      READ(rstrt_unit) iread
      nqloc=NINT(iread)
      CALL tri_linear_alloc(tb_l,lv,nqloc)
      tb_l%name=name
      IF (SIZE(title)<SIZE(tb_l%title)) THEN
        tb_l%title=title(1)
      ELSE
        tb_l%title=title
      ENDIF
      READ(rstrt_unit) tb_l%fs
     
      RETURN
      END SUBROUTINE dump_read_tri_linear_2D
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE dump_read_mod
c-----------------------------------------------------------------------
c     module dump.
c     access to read and write routines.
c-----------------------------------------------------------------------
      MODULE dump
      USE dump_write_mod
      USE dump_read_mod

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE dump

c-----------------------------------------------------------------------
c     subprogram 19. dump_rblock_dealloc.
c     deallocates space for temporary block used in rddump,wrdump
c-----------------------------------------------------------------------
      SUBROUTINE dump_rblock_dealloc(rb)
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
      CALL lagr_quad_dealloc(rb%tele)
      CALL lagr_quad_dealloc(rb%tion)

      RETURN
      END SUBROUTINE dump_rblock_dealloc
c-----------------------------------------------------------------------
c     subprogram 20. dump_tblock_dealloc.
c     deallocates space for temporary block used in rddump,wrdump
c-----------------------------------------------------------------------
      SUBROUTINE dump_tblock_dealloc(tb)
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
      CALL tri_linear_dealloc(tb%tele)
      CALL tri_linear_dealloc(tb%tion)

      RETURN
      END SUBROUTINE dump_tblock_dealloc
c-----------------------------------------------------------------------
c     subprogram 21. dump_seam_dealloc.
c     deallocates space for temporary seam used in rddump,wrdump
c-----------------------------------------------------------------------
      SUBROUTINE dump_seam_dealloc(seam)
      USE local
      USE edge_type_mod
      IMPLICIT NONE

      TYPE(edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: ivert

      DO ivert=1,seam%nvert
        DEALLOCATE(seam%vertex(ivert)%ptr)
      ENDDO
      DEALLOCATE(seam%vertex)
      IF (seam%id == 0) DEALLOCATE(seam%excorner)

      RETURN
      END SUBROUTINE dump_seam_dealloc
