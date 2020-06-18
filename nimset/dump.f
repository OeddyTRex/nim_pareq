c-----------------------------------------------------------------------
c     file dump.f
c     module containing all routines reading and writing data dumps.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  dump
c     1.  dump_write.
c     2.  dump_read.
c     3.  dump_write_rblock.
c     4.  dump_read_rblock.
c     5.  dump_write_tblock.
c     6.  dump_read_tblock.
c     7.  dump_write_seam.
c     8.  dump_read_seam.
c     9.  dump_write_bicube.
c     10. dump_write_lagr_quad.
c     11. dump_write_lagr_quad_2D.
c     12. dump_write_tri_linear.
c     13. dump_write_tri_linear_2D.
c     14. dump_read_bicube.
c     15. dump_read_lagr_quad.
c     16. dump_read_lagr_quad_2D.
c     17. dump_read_tri_linear.
c     18. dump_read_tri_linear_2D.
c     19. dump_read_reset.
c     20. dump_read_rblock_reset.
c     21. dump_read_tblock_reset.
c     22. dump_read_seam_reset.
c-----------------------------------------------------------------------
c     subprogram 0. dump.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE dump
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      IMPLICIT NONE

      LOGICAL :: reset_equil=.false.   !    flag for loading equilibrium
                                       !    fields from reset files.

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. dump_write.
c     writes dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write(nmodes,keff,dump_time,dump_step)
      
      INTEGER(i4), INTENT(IN) :: nmodes,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff

      INTEGER(i4) :: ib,i
      INTEGER(i4) :: dump_code = 0_i4
      LOGICAL :: dump_test,unit_test
      CHARACTER(64) :: step_name
c-----------------------------------------------------------------------
c     test for existing dump file.
c-----------------------------------------------------------------------
      IF (dump_step>999999) THEN
        WRITE(step_name,fmt='(i7.7)') dump_step
      ELSE IF (dump_step>99999) THEN
        WRITE(step_name,fmt='(i6.6)') dump_step
      ELSE
        WRITE(step_name,fmt='(i5.5)') dump_step
      ENDIF
      IF(dump_dir.EQ.".")THEN
         dump_file=TRIM(dump_name)//"."//TRIM(step_name)
      ELSE
         dump_file=TRIM(dump_dir)//"/"//TRIM(dump_name)//"."
     $        //TRIM(step_name)
      ENDIF
      INQUIRE(file=TRIM(dump_file),exist=dump_test)
c-----------------------------------------------------------------------
c     open dump file in appropriate mode.
c-----------------------------------------------------------------------
      IF(dump_test)THEN
         SELECT CASE (dump_over)
      CASE(0)
         CALL open_bin(dump_unit,TRIM(dump_file),'REPLACE','ASIS',64_i4)
      CASE(1)
         CALL open_bin(dump_unit,TRIM(dump_file),'OLD','APPEND',64_i4)
      CASE(2)
         dump_code = 1
         WRITE(nim_wr,*)'Unable to open dump file'
         RETURN
      END SELECT
      ELSE
         CALL open_bin(dump_unit,TRIM(dump_file),'NEW','ASIS',64_i4)
      ENDIF
c-----------------------------------------------------------------------
c     write global data, converting integers to 64 bit reals to have
c     the same # of bits as the reals and to avoid 64 bit integers
c     for linux machines.
c-----------------------------------------------------------------------
      WRITE(dump_unit) dump_time
      WRITE(dump_unit) REAL(dump_step,r8)
      WRITE(dump_unit) REAL(nbl,r8)
      WRITE(dump_unit) REAL(nrbl,r8)
      WRITE(dump_unit) REAL(poly_degree,r8)
      WRITE(dump_unit) REAL(nmodes,r8)
      WRITE(dump_unit) keff
c-----------------------------------------------------------------------
c     write seam data.
c-----------------------------------------------------------------------
      CALL dump_write_seam(seam0)
      DO ib=1,nbl
         CALL dump_write_seam(seam(ib))
      ENDDO
c-----------------------------------------------------------------------
c     write rblock data.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         CALL dump_write_rblock(rb(ib))
      ENDDO
c-----------------------------------------------------------------------
c     write tblock data.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         CALL dump_write_tblock(tb(ib))
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(dump_unit,TRIM(dump_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END  SUBROUTINE dump_write
c-----------------------------------------------------------------------
c     subprogram 2. dump_read.
c     reads dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read(nmodes,keff,dump_time,dump_step,readk,rdfile)

      INTEGER(i4), INTENT(INOUT) :: nmodes
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER, INTENT(INOUT) :: keff
      LOGICAL, INTENT(IN), OPTIONAL :: readk
      CHARACTER(*), INTENT(IN), OPTIONAL :: rdfile

      INTEGER(i4) :: ib,nb,dump_code = 0_i4,nmodes_dump,i
      REAL(r8) :: iread
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      IF (PRESENT(rdfile)) THEN
        INQUIRE(FILE=TRIM(rdfile),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop
     $    ('Dump file '//TRIM(rdfile)//' does not exist.')
        CALL open_bin(rstrt_unit,TRIM(rdfile),'OLD','REWIND',64_i4)
      ELSE
        INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat)
        IF (.NOT.file_stat) CALL nim_stop
     $    ('Dump file '//TRIM(dump_file)//' does not exist.')
        CALL open_bin(rstrt_unit,TRIM(dump_file),'OLD','REWIND',64_i4)
      ENDIF
      READ(rstrt_unit) dump_time
      READ(rstrt_unit) iread
      dump_step=NINT(iread)
      READ(rstrt_unit) iread
      nbl=NINT(iread)
      READ(rstrt_unit) iread
      nrbl=NINT(iread)
c-----------------------------------------------------------------------
c     check for basis consistency.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      IF (NINT(iread)/=poly_degree) THEN
         CALL nim_stop
     $     ('NIMROD input is not consistent with dumped poly_degree.')
      ENDIF
c-----------------------------------------------------------------------
c     check for fourier rep. consistency, then read wavenumber array.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      IF (PRESENT(readk).AND.readk) THEN
        nmodes=NINT(iread)
        ALLOCATE(keff(nmodes))
      ELSE
        nmodes_dump=NINT(iread)
        IF (nmodes_dump/=nmodes) THEN
           CALL nim_stop
     $       ('NIMROD input is not consistent with dumped nmodes.')
        ENDIF
      ENDIF
      READ(rstrt_unit) keff
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      ALLOCATE(seam0)
      ALLOCATE(seam(nbl))
      CALL dump_read_seam(seam0)
      DO ib=1,nbl
         CALL dump_read_seam(seam(ib))
      ENDDO
c-----------------------------------------------------------------------
c     read rblock data.
c-----------------------------------------------------------------------
      ALLOCATE(rb(nrbl))
      DO ib=1,nrbl
         CALL dump_read_rblock(rb(ib))
      ENDDO
c-----------------------------------------------------------------------
c     read tblock data.
c-----------------------------------------------------------------------
      ALLOCATE(tb(nrbl+1:nbl))
      DO ib=nrbl+1,nbl
         CALL dump_read_tblock(tb(ib))
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      IF (PRESENT(rdfile)) THEN
        CALL close_bin(rstrt_unit,TRIM(rdfile))
      ELSE
        CALL close_bin(rstrt_unit,TRIM(dump_file))
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read
c-----------------------------------------------------------------------
c     subprogram 3. dump_write_rblock.
c     writes rblock data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_rblock(rb)

      TYPE (rblock_type), INTENT(IN) :: rb         
      INTEGER(i4) :: ix,iy,iqty,igx,igy
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
c     subprogram 4. dump_read_rblock.
c     reads rblock data from dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_rblock(rb)

      TYPE (rblock_type), INTENT(OUT) :: rb         
      INTEGER(i4) :: ix,iy,iqty,igx,igy
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
      CALL dump_read_lagr_quad_2D(rb%rz,rb%mx,rb%my,'lrz',
     $                            (/' lrz  '/))
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%be_eq,rb%mx,rb%my,'lbq',
     $                            (/' lbq  '/))
      CALL dump_read_lagr_quad(rb%be,rb%mx,rb%my,'be',(/'  be  '/))
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ja_eq,rb%mx,rb%my,'ljq',
     $                            (/' ljq  '/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%ve_eq,rb%mx,rb%my,'lvq',
     $                            (/' lvq  '/))
      CALL dump_read_lagr_quad(rb%ve,rb%mx,rb%my,'ve',(/'  ve  '/))
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%pres_eq,rb%mx,rb%my,'lpr',
     $                            (/' lpr  '/))
      CALL dump_read_lagr_quad(rb%pres,rb%mx,rb%my,'pr',(/' pres '/))
      CALL dump_read_lagr_quad_2D(rb%prese_eq,rb%mx,rb%my,'lpe',
     $                            (/' lpe  '/))
      CALL dump_read_lagr_quad(rb%prese,rb%mx,rb%my,'pe',(/' prese'/))
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%nd_eq,rb%mx,rb%my,'lnq',
     $                            (/' lndq '/))
      CALL dump_read_lagr_quad(rb%nd,rb%mx,rb%my,'nd',(/'  nd  '/))
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_2D(rb%diff_shape,rb%mx,rb%my,'lds',
     $                            (/'ldf_sh'/))
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
c     subprogram 5. dump_write_tblock.
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
c     subprogram 6. dump_read_tblock.
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
     $                             (/'dif_sh'/))
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
c     subprogram 7. dump_write_seam.
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
c     subprogram 8. dump_read_seam.
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
c     subprogram 9. dump_write_bicube.
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
c     subprogram 10.  dump_write_lagr_quad.
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
c     subprogram 11.  dump_write_lagr_quad_2D.
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
c     subprogram 12. dump_write_tri_linear.
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
c     subprogram 13. dump_write_tri_linear_2D.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear_2D(tb_l)

      TYPE(tri_linear_2D_type), INTENT(IN) :: tb_l

      WRITE(dump_unit) REAL(tb_l%nqty,r8)
      WRITE(dump_unit) tb_l%fs
     
      RETURN
      END SUBROUTINE dump_write_tri_linear_2D
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
c     subprogram 19. dump_read_reset.
c     reads the perturbed fields from a dump file with an arbitrary
c     number of Fourier components, and over-writes as much of the new
c     fields as possible.  the geometry represented in the reset_file
c     must match what is created in nimset, however.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_reset(nmodes)

      INTEGER(i4), INTENT(IN) :: nmodes

      INTEGER(i4) :: ib,nb,dump_code = 0_i4,nmodes_dump,i,nm_copy
      REAL(r8) :: iread
      REAL(r8) :: rread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ra1read
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     a 64 bit variable then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(reset_file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Reset dump file '//TRIM(reset_file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(reset_file),'OLD','REWIND',64_i4)
      READ(rstrt_unit) rread
      READ(rstrt_unit) iread
      READ(rstrt_unit) iread
      IF (NINT(iread)/=nbl) CALL nim_stop
     $  (TRIM(dump_file)//' nbl does not match nbl.')
      READ(rstrt_unit) iread
      IF (NINT(iread)/=nrbl) CALL nim_stop
     $  (TRIM(dump_file)//' nrbl does not match nrbl.')
c-----------------------------------------------------------------------
c     check for basis consistency.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      IF (NINT(iread)/=poly_degree) THEN
         CALL nim_stop
     $     ('NIMROD input is not consistent with dumped poly_degree.')
      ENDIF
c-----------------------------------------------------------------------
c     skip wavenumber array.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      nmodes_dump=NINT(iread)
      nm_copy=MIN(nmodes_dump,nmodes)
      ALLOCATE(ra1read(nmodes_dump))
      READ(rstrt_unit) ra1read
      DEALLOCATE(ra1read)
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      CALL dump_read_seam_reset(seam0)
      DO ib=1,nbl
         CALL dump_read_seam_reset(seam(ib))
      ENDDO
c-----------------------------------------------------------------------
c     read rblock data.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         CALL dump_read_rblock_reset(rb(ib),nm_copy)
      ENDDO
c-----------------------------------------------------------------------
c     read tblock data.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         CALL dump_read_tblock_reset(tb(ib),nm_copy)
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(rstrt_unit,TRIM(dump_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_reset
c-----------------------------------------------------------------------
c     subprogram 20. dump_read_rblock_reset.
c     reads perturbed rblock quantities from a dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_rblock_reset(rb,nm_copy)

      TYPE (rblock_type), INTENT(INOUT) :: rb         
      INTEGER(i4), INTENT(IN) :: nm_copy

      INTEGER(i4) :: ix,iy,iqty,igx,igy
      REAL(r8) :: iread
      REAL(r8) :: rread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ra1read
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: ra3read
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ra4read
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ra5read
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      rb%id=NINT(iread)
      READ(rstrt_unit) iread
      IF (NINT(iread)/=rb%mx) THEN
        WRITE(msg,'(2a,i3,a)') TRIM(dump_file),' mx for block ',rb%id,
     $       ' does not match.'
        CALL nim_stop(msg)
      ENDIF
      READ(rstrt_unit)iread
      IF (NINT(iread)/=rb%my) THEN
        WRITE(msg,'(2a,i3,a)') TRIM(dump_file),' my for block ',rb%id,
     $       ' does not match.'
        CALL nim_stop(msg)
      ENDIF
      READ(rstrt_unit) rread
c-----------------------------------------------------------------------
c     coordinates and equilibrium fields may be skipped or read
c     according to the reset_equil flag.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%rz,2_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%rz,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     equilibrium magnetic field.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%be_eq,3_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%be_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed magnetic field.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%be,3_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium current density and fluid velocity.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%ja_eq,3_i4)
        CALL dump_read_lagr_quad_2D_reset(rb%ve_eq,3_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%ja_eq,0_i4)
        CALL dump_read_lagr_quad_2D_reset(rb%ve_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed fluid velocity.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%ve,3_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium pressure.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%pres_eq,1_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%pres_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed pressure.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%pres,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium electron pressure.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%prese_eq,1_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%prese_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed electron pressure.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%prese,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium number density.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%nd_eq,1_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%nd_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed number density.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%nd,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     diffusivity shape factor.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_lagr_quad_2D_reset(rb%diff_shape,1_i4)
      ELSE
        CALL dump_read_lagr_quad_2D_reset(rb%diff_shape,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%conc,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     reset the perturbed temperatures.
c-----------------------------------------------------------------------
      CALL dump_read_lagr_quad_reset(rb%tele,1_i4,nm_copy)
      CALL dump_read_lagr_quad_reset(rb%tion,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     end of executable statements.
c-----------------------------------------------------------------------
      RETURN
c-----------------------------------------------------------------------
c     internal subroutines used in dump_read_rblock_reset:
c-----------------------------------------------------------------------
      CONTAINS
c-----------------------------------------------------------------------
c       modified lagrange 2D structure read:
c-----------------------------------------------------------------------
        SUBROUTINE dump_read_lagr_quad_2D_reset(rb_laq,nq)

        INTEGER(i4), INTENT(IN) :: nq
        TYPE(lagr_quad_2d_type), INTENT(INOUT) :: rb_laq

        READ(rstrt_unit) iread
        rb_laq%nqty=NINT(iread)
        READ(rstrt_unit) iread
        rb_laq%n_side=NINT(iread)
        READ(rstrt_unit) iread
        rb_laq%n_int=NINT(iread)
 
        ALLOCATE(ra3read(rb_laq%nqty,0:rb%mx,0:rb%my))
        READ(rstrt_unit) ra3read
        IF (nq>0) rb_laq%fs=ra3read
        DEALLOCATE(ra3read)
 
        IF (rb_laq%n_side>0) THEN
          ALLOCATE(ra4read(rb_laq%nqty,rb_laq%n_side,1:rb%mx,0:rb%my))
          READ(rstrt_unit) ra4read
          IF (nq>0) rb_laq%fsh=ra4read
          DEALLOCATE(ra4read)
          ALLOCATE(ra4read(rb_laq%nqty,rb_laq%n_side,0:rb%mx,1:rb%my))
          READ(rstrt_unit) ra4read
          IF (nq>0) rb_laq%fsv=ra4read
          DEALLOCATE(ra4read)
          ALLOCATE(ra4read(rb_laq%nqty,rb_laq%n_int,1:rb%mx,1:rb%my))
          READ(rstrt_unit) ra4read
          IF (nq>0) rb_laq%fsi=ra4read
          DEALLOCATE(ra4read)
        ENDIF

        END SUBROUTINE dump_read_lagr_quad_2D_reset
c-----------------------------------------------------------------------
c       modified lagrange quadrilateral structure read:
c-----------------------------------------------------------------------
        SUBROUTINE dump_read_lagr_quad_reset(rb_laq,nq,nm)

        INTEGER(i4), INTENT(IN) :: nq,nm
        TYPE(lagr_quad_type), INTENT(INOUT) :: rb_laq

        INTEGER(i4) :: laq_nq,laq_nf

        READ(rstrt_unit) iread
        laq_nq=NINT(iread)
        READ(rstrt_unit) iread
        laq_nf=NINT(iread)
        READ(rstrt_unit) iread
        READ(rstrt_unit) iread

        ALLOCATE(ra4read(laq_nq,0:rb%mx,0:rb%my,laq_nf))
        READ(rstrt_unit) ra4read
        rb_laq%fs(1:nq,:,:,1:nm)=ra4read(1:nq,:,:,1:nm)
        READ(rstrt_unit) ra4read
        rb_laq%fs(1:nq,:,:,1:nm)=rb_laq%fs(1:nq,:,:,1:nm)+
     $                           (0,1)*ra4read(1:nq,:,:,1:nm)
        DEALLOCATE(ra4read)

        IF (rb_laq%n_side>0) THEN
          ALLOCATE(ra5read(laq_nq,rb_laq%n_side,1:rb%mx,0:rb%my,laq_nf))
          READ(rstrt_unit) ra5read
          rb_laq%fsh(1:nq,:,:,:,1:nm)=ra5read(1:nq,:,:,:,1:nm)
          READ(rstrt_unit) ra5read
          rb_laq%fsh(1:nq,:,:,:,1:nm)=rb_laq%fsh(1:nq,:,:,:,1:nm)+
     $                                (0,1)*ra5read(1:nq,:,:,:,1:nm)
          DEALLOCATE(ra5read)
          ALLOCATE(ra5read(laq_nq,rb_laq%n_side,0:rb%mx,1:rb%my,laq_nf))
          READ(rstrt_unit) ra5read
          rb_laq%fsv(1:nq,:,:,:,1:nm)=ra5read(1:nq,:,:,:,1:nm)
          READ(rstrt_unit) ra5read
          rb_laq%fsv(1:nq,:,:,:,1:nm)=rb_laq%fsv(1:nq,:,:,:,1:nm)+
     $                                (0,1)*ra5read(1:nq,:,:,:,1:nm)
          DEALLOCATE(ra5read)
          ALLOCATE(ra5read(laq_nq,rb_laq%n_int,1:rb%mx,1:rb%my,laq_nf))
          READ(rstrt_unit) ra5read
          rb_laq%fsi(1:nq,:,:,:,1:nm)=ra5read(1:nq,:,:,:,1:nm)
          READ(rstrt_unit) ra5read
          rb_laq%fsi(1:nq,:,:,:,1:nm)=rb_laq%fsi(1:nq,:,:,:,1:nm)+
     $                                (0,1)*ra5read(1:nq,:,:,:,1:nm)
          DEALLOCATE(ra5read)
        ENDIF

        END SUBROUTINE dump_read_lagr_quad_reset
c-----------------------------------------------------------------------
c     end of included subroutines.
c-----------------------------------------------------------------------
      END SUBROUTINE dump_read_rblock_reset
c-----------------------------------------------------------------------
c     subprogram 21. dump_read_tblock_reset.
c     reads perturbed tblock quantities from a dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tblock_reset(tb,nm_copy)

      TYPE (tblock_type), INTENT(INOUT) :: tb
      INTEGER(i4), INTENT(IN) :: nm_copy

      INTEGER(i4) :: ivert,icell,iqty,inbr,nnbr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ia1read
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ia2read
      REAl(r8) :: iread
      REAl(r8) :: rread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ra1read
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: ra3read
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ra4read
      CHARACTER(64) :: msg
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      READ(rstrt_unit) iread
      IF (NINT(iread)/=tb%mvert) THEN
        WRITE(msg,'(2a,i3,a)')TRIM(dump_file),' mvert for block ',tb%id,
     $       ' does not match.'
        CALL nim_stop(msg)
      ENDIF
      READ(rstrt_unit) iread
      IF (NINT(iread)/=tb%mcell) THEN
        WRITE(msg,'(2a,i3,a)')TRIM(dump_file),' mcell for block ',tb%id,
     $       ' does not match.'
        CALL nim_stop(msg)
      ENDIF
c-----------------------------------------------------------------------
c     skip geometry.
c-----------------------------------------------------------------------
      ALLOCATE(ra1read(0:tb%tgeom%mvert))
      READ(rstrt_unit) ra1read
      READ(rstrt_unit) ra1read
      DEALLOCATE(ra1read)
      ALLOCATE(ia2read(tb%tgeom%mcell,3))
      READ(rstrt_unit) ia2read
      DEALLOCATE(ia2read)
      DO ivert=0,tb%tgeom%mvert
         READ(rstrt_unit) iread
         ALLOCATE(ia1read(0:NINT(iread)))
         READ(rstrt_unit) ia1read
         DEALLOCATE(ia1read)
      ENDDO
c-----------------------------------------------------------------------
c     equilibrium magnetic field.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%be_eq,3_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%be_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed magnetic field.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%be,3_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium current density and fluid velocity.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%ja_eq,3_i4)
        CALL dump_read_tri_linear_reset_2D(tb%ve_eq,3_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%ja_eq,0_i4)
        CALL dump_read_tri_linear_reset_2D(tb%ve_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed fluid velocity.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%ve,3_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium pressure.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%pres_eq,1_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%pres_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed pressure.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%pres,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium electron pressure.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%prese_eq,1_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%prese_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed electron pressure.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%prese,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     equilibrium number density.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%nd_eq,1_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%nd_eq,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed number density.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%nd,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     skip the diffusivity shape factor.
c-----------------------------------------------------------------------
      IF (reset_equil) THEN
        CALL dump_read_tri_linear_reset_2D(tb%diff_shape,1_i4)
      ELSE
        CALL dump_read_tri_linear_reset_2D(tb%diff_shape,0_i4)
      ENDIF
c-----------------------------------------------------------------------
c     reset the perturbed material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%conc,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     reset the perturbed temperatures.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_reset(tb%tele,1_i4,nm_copy)
      CALL dump_read_tri_linear_reset(tb%tion,1_i4,nm_copy)
c-----------------------------------------------------------------------
c     end of executable statements.
c-----------------------------------------------------------------------
      RETURN
c-----------------------------------------------------------------------
c     internal subroutines used in dump_read_tblock_reset:
c-----------------------------------------------------------------------
      CONTAINS
c-----------------------------------------------------------------------
c       modified tri_linear structure read:
c-----------------------------------------------------------------------
        SUBROUTINE dump_read_tri_linear_reset(tb_l,nq,nm)

        INTEGER(i4), INTENT(IN) :: nq,nm
        TYPE(tri_linear_type), INTENT(INOUT) :: tb_l

        INTEGER(i4) :: tb_nq,tb_nf

        READ(rstrt_unit) iread
        tb_nq=NINT(iread)
        READ(rstrt_unit) iread
        tb_nf=NINT(iread)

        ALLOCATE(ra4read(tb_nq,0:tb%mvert,0,tb_nf))
        READ(rstrt_unit) ra4read
        tb_l%fs(1:nq,:,:,1:nm)=ra4read(1:nq,:,:,1:nm)
        READ(rstrt_unit) ra4read
        tb_l%fs(1:nq,:,:,1:nm)=tb_l%fs(1:nq,:,:,1:nm)+
     $                         (0,1)*ra4read(1:nq,:,:,1:nm)
        DEALLOCATE(ra4read)

        END SUBROUTINE dump_read_tri_linear_reset
c-----------------------------------------------------------------------
c       modified tri_linear_2D structure read:
c-----------------------------------------------------------------------
        SUBROUTINE dump_read_tri_linear_reset_2D(tb_l,nq)

        INTEGER(i4), INTENT(IN) :: nq
        TYPE(tri_linear_2D_type), INTENT(INOUT) :: tb_l

        INTEGER(i4) :: tb_nq

        READ(rstrt_unit) iread
        tb_nq=NINT(iread)

        ALLOCATE(ra3read(tb_nq,0:tb%mvert,0))
        READ(rstrt_unit) ra3read
        IF (nq>0) tb_l%fs(1:nq,:,:)=ra3read(1:nq,:,:)
        DEALLOCATE(ra3read)

        END SUBROUTINE dump_read_tri_linear_reset_2D
c-----------------------------------------------------------------------
c     end of included subroutines.
c-----------------------------------------------------------------------
      END SUBROUTINE dump_read_tblock_reset
c-----------------------------------------------------------------------
c     subprogram 22. dump_read_seam_reset.
c     reads seam data from dump file for cases that over-write n=0.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_seam_reset(seam)

      TYPE (edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: i,nv,np,ip,iq
      CHARACTER(64) :: msg
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
      IF (NINT(iread)/=seam%nvert) THEN
        WRITE(msg,'(2a,i3,a)')TRIM(dump_file),' seam%nvert for block ',
     $       seam%id,' does not match.'
        CALL nim_stop(msg)
      ENDIF
c-----------------------------------------------------------------------
c     skip seam data.
c-----------------------------------------------------------------------
      ALLOCATE(ia1read(2))
      DO i=1,seam%nvert
         READ(rstrt_unit) iread
         ALLOCATE(ia2read(2,NINT(iread)))
         READ(rstrt_unit) ia2read
         DEALLOCATE(ia2read)
         IF (seam%id>0) READ(rstrt_unit) ia1read
      ENDDO
      DEALLOCATE(ia1read)
c-----------------------------------------------------------------------
c     skip external corner data.
c-----------------------------------------------------------------------
      IF (seam%id == 0_i4 .AND. seam%nvert>0) THEN
         DO i=1,seam%nvert
            READ(rstrt_unit) iread
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN 
      END SUBROUTINE dump_read_seam_reset
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE dump
