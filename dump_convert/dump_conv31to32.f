c-----------------------------------------------------------------------
c     file dump_conv30to31.f
c     converts 3_0* style dump files to 3_1 style files by adding the
c     temperature fields.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  dump_con_mod.
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
c     16. dump_read_tri_linear.
c     17. dump_read_tri_linear_2D.
c     18. dump_conv30to31.
c-----------------------------------------------------------------------
c     subprogram 0. dump_conv.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE dump_con_mod
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      IMPLICIT NONE

      TYPE(lagr_quad_2D_type), DIMENSION(:), POINTER :: lrz,lbq,ljq,lvq,
     $                                       lpq,lpe,lnd,lds

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. dump_write.
c     writes dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write(nmodes,keff,dump_time,dump_step,new_file)
      USE physdat
      
      INTEGER(i4), INTENT(IN) :: nmodes,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(:), INTENT(IN) :: keff
      CHARACTER(*), INTENT(IN) :: new_file

      INTEGER(i4) :: ib,i,ibasis,im,ix,iy,ix0,iy0,lx,ly
      INTEGER(i4) :: dump_code = 0_i4
      REAL(r8) :: dx,dy
      LOGICAL :: dump_test,unit_test
      CHARACTER(64) :: step_name
c-----------------------------------------------------------------------
c     open dump file.
c-----------------------------------------------------------------------
      CALL open_bin(dump_unit,TRIM(new_file),'UNKNOWN','REWIND',64_i4)
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
c     convert bicube mapping and equilibrium fields to lagrange
c     polynomials for version 3_2.
c
c     note that simple evaluations of the bicubic splines are used,
c     so overshoots in the splines may be preserved or worsened.
c-----------------------------------------------------------------------
      ALLOCATE(lrz(nrbl),lbq(nrbl),ljq(nrbl),lvq(nrbl),lpq(nrbl),
     $         lpe(nrbl),lnd(nrbl),lds(nrbl))
      DO ib=1,nrbl
        lx=rb(ib)%mx
        ly=rb(ib)%my
        CALL lagr_quad_alloc(lrz(ib),lx,ly,2_i4,poly_degree,
     $        name='lagrrz',title=(/'  r   ','  z   '/))
        CALL lagr_quad_alloc(lbq(ib),lx,ly,3_i4,poly_degree,
     $        name='lbe_eq',title=(/' be_r ',' be_z ',' be_p '/))
        CALL lagr_quad_alloc(ljq(ib),lx,ly,3_i4,poly_degree,
     $        name='lja_eq',title=(/' ja_r ',' ja_z ',' ja_p '/))
        CALL lagr_quad_alloc(lvq(ib),lx,ly,3_i4,poly_degree,
     $        name='lve_eq',title=(/' ve_r ',' ve_z ',' ve_p '/))
        CALL lagr_quad_alloc(lpq(ib),lx,ly,1_i4,poly_degree,
     $        name='lpres',title=(/'lpres '/))
        CALL lagr_quad_alloc(lpe(ib),lx,ly,1_i4,poly_degree,
     $        name='lprese',title=(/'lprese'/))
        CALL lagr_quad_alloc(lnd(ib),lx,ly,1_i4,poly_degree,
     $        name='lnd_eq',title=(/' lndeq'/))
        CALL lagr_quad_alloc(lds(ib),lx,ly,1_i4,
     $        poly_degree,name='ldiff',title=(/' ldiff'/))
        DO ibasis=1,poly_degree**2
          ix0=rb(ib)%ve%ix0(ibasis)
          iy0=rb(ib)%ve%iy0(ibasis)
          dx=rb(ib)%ve%dx(ibasis)
          dy=rb(ib)%ve%dy(ibasis)
          DO iy=iy0,ly
            DO ix=ix0,lx
              CALL bicube_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lrz(ib),
     $             rb(ib)%rz%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%be_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lbq(ib),
     $             rb(ib)%be_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%ja_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(ljq(ib),
     $             rb(ib)%ja_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%ve_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lvq(ib),
     $             rb(ib)%ve_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%pres_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lpq(ib),
     $             rb(ib)%pres_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%prese_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lpe(ib),
     $             rb(ib)%prese_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%nd_eq,ix-ix0+dx,iy-iy0+dy,0_i4)
              CALL lagr_quad_basis_assign_loc(lnd(ib),
     $             rb(ib)%nd_eq%f,ibasis,ix,iy)
              CALL bicube_eval(rb(ib)%diff_shape,ix-ix0+dx,iy-iy0+dy,
     $                         0_i4)
              CALL lagr_quad_basis_assign_loc(lds(ib),
     $             rb(ib)%diff_shape%f,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write rblock data.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         CALL dump_write_rblock(rb(ib),lrz(ib),lbq(ib),ljq(ib),
     $                          lvq(ib),lpq(ib),lpe(ib),lnd(ib),lds(ib))
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
      SUBROUTINE dump_read(nmodes,keff,dump_time,dump_step,old_file)

      INTEGER(i4), INTENT(OUT) :: nmodes
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff
      CHARACTER(*), INTENT(IN) :: old_file

      INTEGER(i4) :: ib,nb,dump_code = 0_i4,nmodes_dump,i
      REAL(r8) :: iread
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
      dump_step=NINT(iread)
      READ(rstrt_unit) iread
      nbl=NINT(iread)
      READ(rstrt_unit) iread
      nrbl=NINT(iread)
c-----------------------------------------------------------------------
c     spatial representation parameters.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      poly_degree=iread
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
      CALL close_bin(rstrt_unit,TRIM(dump_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read
c-----------------------------------------------------------------------
c     subprogram 3. dump_write_rblock.
c     writes rblock data to dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_rblock(rb,lrz,lbq,ljq,lvq,lpq,lpe,lnd,lds)

      TYPE (rblock_type), INTENT(IN) :: rb         
      TYPE(lagr_quad_2D_type) :: lrz,lbq,ljq,lvq,lpq,lpe,lnd,lds

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
      CALL dump_write_lagr_quad_2D(lrz)
c-----------------------------------------------------------------------
c     write magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(lbq)
      CALL dump_write_lagr_quad(rb%be)
c-----------------------------------------------------------------------
c     write equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(ljq)
c-----------------------------------------------------------------------
c     write fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(lvq)
      CALL dump_write_lagr_quad(rb%ve)
c-----------------------------------------------------------------------
c     write pressures.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(lpq)
      CALL dump_write_lagr_quad(rb%pres)
      CALL dump_write_lagr_quad_2D(lpe)
      CALL dump_write_lagr_quad(rb%prese)
c-----------------------------------------------------------------------
c     write number densities.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(lnd)
      CALL dump_write_lagr_quad(rb%nd)
c-----------------------------------------------------------------------
c     write diffusivity shape function and material concentration.
c-----------------------------------------------------------------------
      CALL dump_write_lagr_quad_2D(lds)
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
      CALL dump_read_bicube(rb%rz,rb%mx,rb%my,'rz',
     $                      (/'   r  ','   z  '/))
c-----------------------------------------------------------------------
c     read magnetic fields.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%be_eq,rb%mx,rb%my,'bq',(/' be_eq'/))
      CALL dump_read_lagr_quad(rb%be,rb%mx,rb%my,'be',(/'  be  '/))
c-----------------------------------------------------------------------
c     read equilibrium current density.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%ja_eq,rb%mx,rb%my,'jq',(/' ja_eq'/))
c-----------------------------------------------------------------------
c     read fluid velocities.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%ve_eq,rb%mx,rb%my,'vq',(/' ve_eq'/))
      CALL dump_read_lagr_quad(rb%ve,rb%mx,rb%my,'ve',(/'  ve  '/))
c-----------------------------------------------------------------------
c     read pressures.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%pres_eq,rb%mx,rb%my,'pq',(/' pr_eq'/))
      CALL dump_read_lagr_quad(rb%pres,rb%mx,rb%my,'pr',(/' pres '/))
      CALL dump_read_bicube(rb%prese_eq,rb%mx,rb%my,'pq',(/'pre_eq'/))
      CALL dump_read_lagr_quad(rb%prese,rb%mx,rb%my,'pe',(/' prese'/))
c-----------------------------------------------------------------------
c     read number densities.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%nd_eq,rb%mx,rb%my,'nq',(/' nd_eq'/))
      CALL dump_read_lagr_quad(rb%nd,rb%mx,rb%my,'nd',(/'  nd  '/))
c-----------------------------------------------------------------------
c     read diffusivity shape factor and material concentration.
c-----------------------------------------------------------------------
      CALL dump_read_bicube(rb%diff_shape,rb%mx,rb%my,'ds',(/'dif_sh'/))
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
c     subprogram 11.  dump_write_lagr_quad_2D.
c     2D version.
c-----------------------------------------------------------------------
      SUBROUTINE dump_write_lagr_quad_2D(laq)
 
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
 
      WRITE(dump_unit) REAL(laq%nqty,r8)
      WRITE(dump_unit) REAL(laq%n_side,r8)
      WRITE(dump_unit) REAL(laq%n_int,r8)
      WRITE(dump_unit) REAL(laq%fs,r8)
      IF (ASSOCIATED(laq%fsh)) THEN
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

      REAL(r8) :: iread

      READ(rstrt_unit) iread
      rb_bc%nqty=NINT(iread)
      CALL bicube_alloc(rb_bc,lx,ly,rb_bc%nqty)
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

      REAL(r8) :: iread
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread
      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: rread2

      READ(rstrt_unit) iread
      laq%nqty=NINT(iread)
      READ(rstrt_unit) iread
      laq%nfour=NINT(iread)
      READ(rstrt_unit) iread
      laq%n_side=NINT(iread)
      READ(rstrt_unit) iread
      laq%n_int=NINT(iread)
      CALL lagr_quad_alloc(laq,lx,ly,laq%nqty,laq%nfour,
     $                     laq%n_side+1_i4)
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
c     subprogram 16. dump_read_tri_linear.
c     real and imaginary parts are read separately for possible
c     data conversion.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear(tb_l,lv,name,title)

      INTEGER(i4), INTENT(IN) :: lv
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_type), INTENT(OUT) :: tb_l

      REAL(r8) :: iread
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: rread

      READ(rstrt_unit) iread
      tb_l%nqty=NINT(iread)
      READ(rstrt_unit) iread
      tb_l%nfour=NINT(iread)
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty,tb_l%nfour)
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
c     subprogram 17. dump_read_tri_linear_2D.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_tri_linear_2D(tb_l,lv,name,title)

      INTEGER(i4), INTENT(IN) :: lv
      CHARACTER(*), INTENT(IN) :: name
      CHARACTER(6), DIMENSION(:), INTENT(IN) :: title
      TYPE(tri_linear_2D_type), INTENT(OUT) :: tb_l

      REAL(r8) :: iread

      READ(rstrt_unit) iread
      tb_l%nqty=NINT(iread)
      CALL tri_linear_alloc(tb_l,lv,tb_l%nqty)
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
c     close module
c-----------------------------------------------------------------------
      END MODULE dump_con_mod


c-----------------------------------------------------------------------
c     18.  main program, dump_conv31to32
c-----------------------------------------------------------------------
      PROGRAM dump_conv31to32
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
      WRITE(nim_wr,*) "Enter the names of the old (nimrod3_1) dump file"
      WRITE(nim_wr,*) "and the new (nimrod3_2) dump file."
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
      END PROGRAM dump_conv31to32

