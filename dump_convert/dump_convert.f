c-----------------------------------------------------------------------
c     file dump_covert.f
c     modified routines for reading nimrod2_* dump
c     files and for writing nimrod3_1 dump files.
c
c     main program is at the end of this file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  dump_con_mod
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
c     11. dump_write_tri_linear.
c     12. dump_write_tri_linear_2D.
c     13. dump_read_bicube.
c     14. dump_read_lagr_quad.
c     15. dump_read_tri_linear.
c     16. dump_read_tri_linear_2D.
c     17. dump_convert.
c-----------------------------------------------------------------------
c     subprogram 0. dump_con_mod.
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE dump_con_mod
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      IMPLICIT NONE

      INTEGER(i4) :: poly_degree

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. dump_write.
c     writes dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write(nmodes,keff,dump_time,dump_step,new_file)
      
      INTEGER(i4), INTENT(IN) :: nmodes,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff
      CHARACTER(*), INTENT(IN) :: new_file

      INTEGER(i4) :: ib,i,im
      INTEGER(i4) :: dump_code = 0_i4
      LOGICAL :: dump_test,unit_test
      CHARACTER(64) :: step_name
c-----------------------------------------------------------------------
c     open dump file in appropriate mode.
c-----------------------------------------------------------------------
      CALL open_bin(dump_unit,TRIM(new_file),'UNKNOWN','ASIS',64_i4)
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
c     allocate and fill fields that have been added since version 2.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
        CALL lagr_quad_alloc(rb(ib)%tele,rb(ib)%mx,rb(ib)%my,1_i4,
     $                       nmodes,poly_degree,name='te',
     $                       title=(/' tele '/))
        CALL lagr_quad_alloc(rb(ib)%tion,rb(ib)%mx,rb(ib)%my,1_i4,
     $                       nmodes,poly_degree,name='ti',
     $                       title=(/' tion '/))
        DO im=1,nmodes
          rb(ib)%tele%fs(:,:,:,im)=
     $           rb(ib)%prese%fs(:,:,:,im)/rb(ib)%nd_eq%fs
          rb(ib)%tion%fs(:,:,:,im)=
     $          (rb(ib)%pres %fs(:,:,:,im)
     $          -rb(ib)%prese%fs(:,:,:,im))/rb(ib)%nd_eq%fs
        ENDDO
      ENDDO
      DO ib=nrbl+1,nbl
        CALL tri_linear_alloc(tb(ib)%tele,tb(ib)%tgeom%mvert,1_i4,
     $                        nmodes,name='te',title=(/' tele '/))
        CALL tri_linear_alloc(tb(ib)%tion,tb(ib)%tgeom%mvert,1_i4,
     $                        nmodes,name='ti',title=(/' tion '/))
        DO im=1,nmodes
          tb(ib)%tele%fs(:,:,:,im)=
     $           tb(ib)%prese%fs(:,:,:,im)/tb(ib)%nd_eq%fs
          tb(ib)%tion%fs(:,:,:,im)=
     $          (tb(ib)%pres %fs(:,:,:,im)
     $          -tb(ib)%prese%fs(:,:,:,im))/tb(ib)%nd_eq%fs
        ENDDO
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
      CALL close_bin(dump_unit,TRIM(new_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END  SUBROUTINE dump_write
c-----------------------------------------------------------------------
c     subprogram 2. dump_read.
c     reads nimrod2_* dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read(nmodes,keff,dump_time,dump_step,old_file)

      INTEGER(i4), INTENT(OUT) :: nmodes
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff
      CHARACTER(*), INTENT(IN) :: old_file

      INTEGER(i4) :: ib,nb,dump_code = 0_i4,i
      INTEGER(i8) :: iread
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(old_file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Dump file '//TRIM(old_file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(old_file),'OLD','REWIND',64_i4)
      READ(rstrt_unit) dump_time
      READ(rstrt_unit) iread
      dump_step=iread
      READ(rstrt_unit) iread
      nbl=iread
      READ(rstrt_unit) iread
      nrbl=iread
c-----------------------------------------------------------------------
c     nimrod2_* has linear elements only.
c-----------------------------------------------------------------------
      poly_degree=1
      READ(rstrt_unit) iread
      nmodes=iread
      ALLOCATE(keff(nmodes))
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
         CALL dump_read_rblock(rb(ib),nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     read tblock data.
c-----------------------------------------------------------------------
      ALLOCATE(tb(nrbl+1:nbl))
      DO ib=nrbl+1,nbl
         CALL dump_read_tblock(tb(ib),nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(rstrt_unit,TRIM(old_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read
c-----------------------------------------------------------------------
c     subprogram 3. dump_write_rblock.
c     writes rblock data to dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
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
      CALL dump_write_bicube(rb%rz)
c-----------------------------------------------------------------------
c     write magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%be_eq)
      CALL dump_write_lagr_quad(rb%be)
c-----------------------------------------------------------------------
c     write equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%ja_eq)
c-----------------------------------------------------------------------
c     write fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%ve_eq)
      CALL dump_write_lagr_quad(rb%ve)
c-----------------------------------------------------------------------
c     write pressures.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%pres_eq)
      CALL dump_write_lagr_quad(rb%pres)
      CALL dump_write_bicube(rb%prese_eq)
      CALL dump_write_lagr_quad(rb%prese)
c-----------------------------------------------------------------------
c     write number densities.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%nd_eq)
      CALL dump_write_lagr_quad(rb%nd)
c-----------------------------------------------------------------------
c     write diffusivity shape function and material concentration.
c-----------------------------------------------------------------------
      CALL dump_write_bicube(rb%diff_shape)
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
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_rblock(rb,nmodes)

      TYPE (rblock_type), INTENT(OUT) :: rb         
      INTEGER(i4), INTENT(IN) :: nmodes

      INTEGER(i4) :: ix,iy,iqty,igx,igy
      INTEGER(i8) :: iread
      REAL(r8) :: rread
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      rb%name='rblock'
      READ(rstrt_unit) iread
      rb%id=iread
      READ(rstrt_unit) iread
      rb%mx=iread
      READ(rstrt_unit) iread
      rb%my=iread
      READ(rstrt_unit) rread
      IF (rread>0) THEN
        rb%degenerate=.true.
      ELSE
        rb%degenerate=.false.
      ENDIF
c-----------------------------------------------------------------------
c     read coordinates.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%rz,rb%mx,rb%my,'rz',
     $                      (/'   r  ','   z  '/))
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%be_eq,rb%mx,rb%my,'bq',(/' be_eq'/))
      CALL dump_read_lagr_quad(rb%be,rb%mx,rb%my,'be',(/'  be  '/),
     $                         nmodes)
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%ja_eq,rb%mx,rb%my,'jq',(/' ja_eq'/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%ve_eq,rb%mx,rb%my,'vq',(/' ve_eq'/))
      CALL dump_read_lagr_quad(rb%ve,rb%mx,rb%my,'ve',(/'  ve  '/),
     $                         nmodes)
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%pres_eq,rb%mx,rb%my,'pq',(/' pr_eq'/))
      CALL dump_read_lagr_quad(rb%pres,rb%mx,rb%my,'pr',(/' pres '/),
     $                         nmodes)
      CALL dump_read_bicube(rb%prese_eq,rb%mx,rb%my,'pq',(/'pre_eq'/))
      CALL dump_read_lagr_quad(rb%prese,rb%mx,rb%my,'pe',(/' prese'/),
     $                         nmodes)
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%nd_eq,rb%mx,rb%my,'nq',(/' nd_eq'/))
      CALL dump_read_lagr_quad(rb%nd,rb%mx,rb%my,'nd',(/'  nd  '/),
     $                         nmodes)
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%diff_shape,rb%mx,rb%my,'ds',(/'dif_sh'/))
      CALL dump_read_lagr_quad(rb%conc,rb%mx,rb%my,'co',(/' conc '/),
     $                         nmodes)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_rblock
c-----------------------------------------------------------------------
c     subprogram 5. dump_write_tblock.
c     writes tblock data to dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
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
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tblock(tb,nmodes)

      TYPE (tblock_type), INTENT(OUT) :: tb
      INTEGER(i4), INTENT(IN) :: nmodes

      INTEGER(i4) :: ivert,icell,iqty,inbr,nnbr
      INTEGER(i8) :: iread
      INTEGER(i8), DIMENSION(:), ALLOCATABLE :: ia1read
      INTEGER(i8), DIMENSION(:,:), ALLOCATABLE :: ia2read
c-----------------------------------------------------------------------
c     read geometry.
c-----------------------------------------------------------------------
      tb%name='tblock'
      READ(rstrt_unit) iread
      tb%id=iread
      READ(rstrt_unit) iread
      tb%tgeom%mvert=iread
      READ(rstrt_unit) iread
      tb%tgeom%mcell=iread
      CALL tri_linear_geom_alloc(tb%tgeom,tb%tgeom%mvert,tb%tgeom%mcell)
      READ(rstrt_unit) tb%tgeom%xs
      READ(rstrt_unit) tb%tgeom%ys
      ALLOCATE(ia2read(tb%tgeom%mcell,3))
      READ(rstrt_unit) ia2read
      tb%tgeom%vertex=ia2read
      DEALLOCATE(ia2read)
      DO ivert=0,tb%tgeom%mvert
         READ(rstrt_unit) iread
         nnbr=iread
         ALLOCATE(tb%tgeom%neighbor(ivert)%vertex(0:nnbr))
         ALLOCATE(ia1read(0:nnbr))
         READ(rstrt_unit) ia1read
         tb%tgeom%neighbor(ivert)%vertex=ia1read
         DEALLOCATE(ia1read)
      ENDDO
      tb%mvert=tb%tgeom%mvert
      tb%mcell=tb%tgeom%mcell
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%be_eq,tb%mvert,'bq',(/' be_eq'/))
      CALL dump_read_tri_linear(tb%be,tb%mvert,'be',(/'  be  '/),
     $                          nmodes)
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ja_eq,tb%mvert,'jq',(/' ja_eq'/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%ve_eq,tb%mvert,'vq',(/' ve_eq'/))
      CALL dump_read_tri_linear(tb%ve,tb%mvert,'ve',(/'  ve  '/),
     $                          nmodes)
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%pres_eq,tb%mvert,'pq',
     $                             (/' pr_eq'/))
      CALL dump_read_tri_linear(tb%pres,tb%mvert,'pr',(/' pres '/),
     $                          nmodes)
      CALL dump_read_tri_linear_2D(tb%prese_eq,tb%mvert,'pq',
     $                             (/'pre_eq'/))
      CALL dump_read_tri_linear(tb%prese,tb%mvert,'pe',(/' prese'/),
     $                          nmodes)
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%nd_eq,tb%mvert,'nq',(/' nd_eq'/))
      CALL dump_read_tri_linear(tb%nd,tb%mvert,'nd',(/'  nd  '/),
     $                          nmodes)
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_tri_linear_2D(tb%diff_shape,tb%mvert,'ds',
     $                             (/'dif_sh'/))
      CALL dump_read_tri_linear(tb%conc,tb%mvert,'co',(/' conc '/),
     $                          nmodes)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_tblock
c-----------------------------------------------------------------------
c     subprogram 7. dump_write_seam.
c     writes seam data to dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
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
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_seam(seam)

      TYPE (edge_type), INTENT(INOUT) :: seam

      INTEGER(i4) :: i,nv,np,ip,iq
      INTEGER(i8) :: iread
      INTEGER(i8), DIMENSION(:), ALLOCATABLE :: ia1read
      INTEGER(i8), DIMENSION(:,:), ALLOCATABLE :: ia2read
c-----------------------------------------------------------------------
c     read block descriptors.
c-----------------------------------------------------------------------
      seam%name='seam'
      READ(rstrt_unit) iread
      seam%id=iread
      READ(rstrt_unit) iread
      seam%nvert=iread
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      ALLOCATE(seam%vertex(seam%nvert))
      ALLOCATE(ia1read(2))
      DO i=1,seam%nvert
         READ(rstrt_unit) iread
         np=iread
         ALLOCATE(seam%vertex(i)%ptr(2,np))
         ALLOCATE(ia2read(2,np))
         READ(rstrt_unit) ia2read
         seam%vertex(i)%ptr=ia2read
         DEALLOCATE(ia2read)
         IF (seam%id>0) THEN
           READ(rstrt_unit) ia1read
           seam%vertex(i)%intxy=ia1read
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
            seam%excorner(i)=(iread>0)
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
      IF (ASSOCIATED(laq%fsh)) THEN
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
c     subprogram 11. dump_write_tri_linear.
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
c     subprogram 12. dump_write_tri_linear_2D.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_tri_linear_2D(tb_l)

      TYPE(tri_linear_2D_type), INTENT(IN) :: tb_l

      WRITE(dump_unit) REAL(tb_l%nqty,r8)
      WRITE(dump_unit) tb_l%fs
     
      RETURN
      END SUBROUTINE dump_write_tri_linear_2D
c-----------------------------------------------------------------------
c     subprogram 13. dump_read_bicube.
c     old bicubes have quantity index last.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_bicube(rb_bc,lx,ly,name,title)

      INTEGER(i4), INTENT(IN) :: lx,ly
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(bicube_type), INTENT(OUT) :: rb_bc

      INTEGER(i8) :: iread
      INTEGER(i4) :: iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread

      READ(rstrt_unit) iread
      rb_bc%nqty=iread

      CALL bicube_alloc(rb_bc,lx,ly,rb_bc%nqty)
      ALLOCATE(rread(0:lx,0:ly,rb_bc%nqty))
      rb_bc%name=name
      IF (SIZE(title)<SIZE(rb_bc%title)) THEN
        rb_bc%title=title(1)
      ELSE
        rb_bc%title=title
      ENDIF
      READ(rstrt_unit) rb_bc%xs
      READ(rstrt_unit) rb_bc%ys
      READ(rstrt_unit) rread
      DO iq=1,rb_bc%nqty
        rb_bc%fs(iq,:,:)=rread(:,:,iq)
      ENDDO
      READ(rstrt_unit) rread
      DO iq=1,rb_bc%nqty
        rb_bc%fsx(iq,:,:)=rread(:,:,iq)
      ENDDO
      READ(rstrt_unit) rread
      DO iq=1,rb_bc%nqty
        rb_bc%fsy(iq,:,:)=rread(:,:,iq)
      ENDDO
      READ(rstrt_unit) rread
      DO iq=1,rb_bc%nqty
        rb_bc%fsxy(iq,:,:)=rread(:,:,iq)
      ENDDO
      DEALLOCATE(rread)

      RETURN
      END SUBROUTINE dump_read_bicube
c-----------------------------------------------------------------------
c     subprogram 14. dump_read_lagr_quad.
c     actually reads bilinear structures and convertes.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_lagr_quad(laq,lx,ly,name,title,nm)

      INTEGER(i4), INTENT(IN) :: lx,ly,nm
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(lagr_quad_type), INTENT(OUT) :: laq

      INTEGER(i8) :: iread
      INTEGER(i4) :: nq_m,im,iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread
      REAL(r8), DIMENSION(:), ALLOCATABLE :: sread

      READ(rstrt_unit) iread
      nq_m=iread/(2*nm)
      laq%nqty=nq_m
      laq%nfour=nm
      laq%n_side=0
      laq%n_int=0
      CALL lagr_quad_alloc(laq,lx,ly,laq%nqty,laq%nfour,1_i4)
      laq%name=name
      IF (SIZE(title)<SIZE(laq%title)) THEN
        laq%title=title(1)
      ELSE
        laq%title=title
      ENDIF

      ALLOCATE(sread(0:lx))
      READ(rstrt_unit) sread
      DEALLOCATE(sread)
      ALLOCATE(sread(0:ly))
      READ(rstrt_unit) sread
      DEALLOCATE(sread)

      ALLOCATE(rread(0:lx,0:ly,iread))
      READ(rstrt_unit) rread
      DO im=1,nm
        DO iq=1,nq_m
          laq%fs(iq,:,:,im)=rread(:,:,2*nq_m*(im-1)+iq)
     $               +(0,1)*rread(:,:,2*nq_m*(im-1)+iq+nq_m)
        ENDDO
      ENDDO
      DEALLOCATE(rread)

      RETURN
      END SUBROUTINE dump_read_lagr_quad
c-----------------------------------------------------------------------
c     subprogram 15. dump_read_tri_linear.
c     old tri_linear data has quantity index last and mixed with Fourier
c     index.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear(tb_l,lv,name,title,nm)

      INTEGER(i4), INTENT(IN) :: lv,nm
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_type), INTENT(OUT) :: tb_l

      INTEGER(i8) :: iread
      INTEGER(i4) :: nq_m,im,iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread

      READ(rstrt_unit) iread
      nq_m=iread/(2*nm)
      tb_l%nqty=nq_m
      tb_l%nfour=nm
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty,tb_l%nfour)
      tb_l%name=name
      IF (SIZE(title)<SIZE(tb_l%title)) THEN
        tb_l%title=title(1)
      ELSE
        tb_l%title=title
      ENDIF

      ALLOCATE(rread(0:lv,0:0,iread))
      READ(rstrt_unit) rread
      DO im=1,nm
        DO iq=1,nq_m
          tb_l%fs(iq,:,:,im)=rread(:,:,2*nq_m*(im-1)+iq)
     $                +(0,1)*rread(:,:,2*nq_m*(im-1)+iq+nq_m)
        ENDDO
      ENDDO
      DEALLOCATE(rread)
     
      RETURN
      END SUBROUTINE dump_read_tri_linear
c-----------------------------------------------------------------------
c     subprogram 16. dump_read_tri_linear_2D.
c     this reads an old tri_linear with quantity index last.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear_2D(tb_l,lv,name,title)

      INTEGER(i4), INTENT(IN) :: lv
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_2D_type), INTENT(OUT) :: tb_l

      INTEGER(i8) :: iread
      INTEGER(i4) :: iq
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rread

      READ(rstrt_unit) iread
      tb_l%nqty=iread
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty)
      tb_l%name=name
      IF (SIZE(title)<SIZE(tb_l%title)) THEN
        tb_l%title=title(1)
      ELSE
        tb_l%title=title
      ENDIF
     
      ALLOCATE(rread(0:lv,0:0,iread))
      READ(rstrt_unit) rread
      DO iq=1,iread
        tb_l%fs(iq,:,:)=rread(:,:,iq)
      ENDDO
      DEALLOCATE(rread)
     
      RETURN
      END SUBROUTINE dump_read_tri_linear_2D
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE dump_con_mod


c-----------------------------------------------------------------------
c     17.  main program, dump_convert
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM dump_convert
      USE local
      USE dump_con_mod
      USE fields
      IMPLICIT NONE

      CHARACTER(128) :: old_file,new_file
      INTEGER(i4) :: nmodes,dump_step,ib,ix,iy,iv,jv,icell
      REAL(r8) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff
c-----------------------------------------------------------------------
c     get files names and read old dump file.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) "Enter the names of the old (nimrod2_*) dump file"
      WRITE(nim_wr,*) "and the new (nimrod3_*) dump file."
      READ(nim_rd,*) old_file,new_file

      CALL dump_read(nmodes,keff,dump_time,dump_step,old_file)
c-----------------------------------------------------------------------
c     draw grid--lifted from diagnose.f
c-----------------------------------------------------------------------
      CALL open_bin(grid_unit,"grid.bin","UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     draw rblocks.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         DO iy=0,rb(ib)%my
            DO ix=0,rb(ib)%mx
               WRITE(grid_unit)REAL(rb(ib)%rz%fs(:,ix,iy),4)
            ENDDO
            WRITE(grid_unit)
         ENDDO
         DO ix=0,rb(ib)%mx
            DO iy=0,rb(ib)%my
               WRITE(grid_unit)REAL(rb(ib)%rz%fs(:,ix,iy),4)
            ENDDO
            WRITE(grid_unit)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     draw tblocks.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         DO icell=1,tb(ib)%mcell
            DO iv=1,3
               jv=tb(ib)%tgeom%vertex(icell,iv)
               WRITE(grid_unit)REAL(tb(ib)%tgeom%xs(jv),4),
     $              REAL(tb(ib)%tgeom%ys(jv),4)
            ENDDO
            jv=tb(ib)%tgeom%vertex(icell,1)
            WRITE(grid_unit)REAL(tb(ib)%tgeom%xs(jv),4),
     $           REAL(tb(ib)%tgeom%ys(jv),4)
            WRITE(grid_unit)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CALL close_bin(grid_unit,"grid.bin")
c-----------------------------------------------------------------------
c     write new dump file.
c-----------------------------------------------------------------------
      CALL dump_write(nmodes,keff,dump_time,dump_step,new_file)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')
      END PROGRAM dump_convert

