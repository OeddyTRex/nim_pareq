c-----------------------------------------------------------------------
c     file dump_reset.f
c     module contains routines for reading dump files for regridding
c     purposes.  The routines read the dump file information into
c     the variables rbc, tbc, seam0c, seamc.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. dumpc
c     1. dumpc_read
c     2. rz_to_cell
c     3. make_cell
c     4. map_cell_position
c     5. print_cell_position
c     6. seam_to_cell
c     7. refine_cell
c     9. make_search_map
c-----------------------------------------------------------------------
c     subprogram 0. dumpc.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE dumpc
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE cell_type_mod
      USE edge_type_mod
      USE vector_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      USE dump

      IMPLICIT NONE

      TYPE :: cheat
         TYPE(cell_type),POINTER :: p
      END TYPE cheat
      TYPE :: remap
        TYPE(cheat), DIMENSION(:,:), POINTER :: p
      END TYPE remap
      INTEGER(i4) :: nblc,nrblc,maxcell,poly_degreec
      TYPE(rblock_type), DIMENSION(:), POINTER :: rbc
      TYPE(tblock_type), DIMENSION(:), POINTER :: tbc
      TYPE(edge_type) :: seam0c
      TYPE(edge_type), DIMENSION(:), POINTER :: seamc

      TYPE(cell_type),POINTER :: start,extbnd,lucky
      INTEGER(i4),DIMENSION(3,2) :: ic

      TYPE :: fail_list
         TYPE(fail_list),POINTER :: next
         INTEGER(i4) :: ib,ibasis,ix,iy
      END TYPE fail_list
 
      TYPE(fail_list),POINTER :: fail_start,fail_cnt

c-----------------------------------------------------------------------
c     the bl2cl structure had been called temp and was deallocated at
c     the end of the make_cell routine.  however, it is a fast way
c     to find a cell pointer if the block and element indices are known,
c     so it is now retained.
c-----------------------------------------------------------------------
      TYPE(remap), DIMENSION(:), ALLOCATABLE :: bl2cl

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. dumpc_read.
c     reads dump file.
c-----------------------------------------------------------------------
      SUBROUTINE dumpc_read(nmodes,keff)

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), INTENT(IN), DIMENSION(nmodes):: keff
      REAL(r8), DIMENSION(:), ALLOCATABLE :: keffc
      LOGICAL :: extrap_uninit=.TRUE.	! If FALSE, stops data extrapolation.

      INTEGER(i4) :: ib,ix,iy,nmodesc,ibold,ibx,iby
      INTEGER(i4) :: imode,jmode,iq,jq,iqs,jqs,ivold,icell
      INTEGER(i4) :: iv,iv1,iblast,iclast
      REAL(r8) :: iread
      LOGICAL :: file_stat,uninit_flag=.FALSE.
      LOGICAL :: local_debug=.FALSE.
      REAL(r8) :: dump_timec,x,y,dx,dy
      REAL(r8) :: tol=1.0e-07_r8
      TYPE(cell_type), POINTER :: item,itemt,item_pre
      TYPE(location_type):: p0
      INTEGER(i4) :: ip,ix0,ix1,iy0,iy1,ibasis,max_basisr,max_basist
      LOGICAL :: failure,data_uninit=.FALSE.
      COMPLEX(r8), DIMENSION(:,:), ALLOCATABLE :: lbe,lve,lpr,lpe,
     $             lnd,lco,lte,lti
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     a 64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(reset_file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Dump file '//TRIM(reset_file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(reset_file),'OLD','REWIND',64_i4)
      READ(rstrt_unit) dump_timec
      READ(rstrt_unit) iread
      READ(rstrt_unit) iread
      nblc=NINT(iread)
      READ(rstrt_unit) iread
      nrblc=NINT(iread)
      READ(rstrt_unit) iread
      poly_degreec=NINT(iread)
      READ(rstrt_unit) iread
      nmodesc=NINT(iread)
      ALLOCATE(keffc(nmodesc))
      READ(rstrt_unit) keffc
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      ALLOCATE(seamc(nblc))
      CALL dump_read_seam(seam0c)
      DO ib=1,nblc
         CALL dump_read_seam(seamc(ib))
      ENDDO
c-----------------------------------------------------------------------
c     read rblock data.
c-----------------------------------------------------------------------
      ALLOCATE(rbc(nrblc))
      DO ib=1,nrblc
         CALL dump_read_rblock(rbc(ib))
      ENDDO
c-----------------------------------------------------------------------
c     read tblock data.
c-----------------------------------------------------------------------
      ALLOCATE(tbc(nrblc+1:nblc))
      DO ib=nrblc+1,nblc
         CALL dump_read_tblock(tbc(ib))
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(rstrt_unit,TRIM(reset_file))
c-----------------------------------------------------------------------
c     Build cell information
c-----------------------------------------------------------------------
      CALL make_cell
      CALL make_search_map
      item => start
      DO WHILE(ASSOCIATED(item))
        CALL map_cell_position(item)
        item => item%next
      ENDDO
c-----------------------------------------------------------------------
c     DEBUG code:  This is used to plot each cell and its connecting cells.
c-----------------------------------------------------------------------
c     CALL open_bin(xdr_unit,'temp.bin','REPLACE','ASIS',64_i4)
c     item => start
c     DO WHILE(ASSOCIATED(item))
c       WRITE(6,*)" # CELL",item%id-1
c       CALL print_cell_position(item)
c       DO ip=1,item%nodes
c         itemt => item%face(ip)%p
c         WRITE(6,*)" # face",ip,itemt%id
c         CALL print_cell_position(itemt)
c       ENDDO
c       WRITE(xdr_unit)
c       item => item%next
c     ENDDO
c     CLOSE(UNIT=xdr_unit)
c-----------------------------------------------------------------------
c     DEBUG code:  This is used to plot the search path when the call
c     to rz_to_cell is true.  The first rz_to_cell puts the code in
c     the same start position as where the code should be for the problem
c     cell associated with the second cell.   The output is path.bin.
c-----------------------------------------------------------------------
      IF(local_debug)THEN
        item => start
        ib=1
        ibasis = 2
        ix=5
        iy=5
        ix0=rb(ib)%be%ix0(ibasis)
        iy0=rb(ib)%be%iy0(ibasis)
        dx=rb(ib)%be%dx(ibasis)
        dy=rb(ib)%be%dy(ibasis)
        CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
        p0%point(1:2) = rb(ib)%rz%f
        CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
        CALL refine_cell(p0,item,failure,x,y,.FALSE.)
        WRITE(6,*)
        WRITE(6,*)"Logical Coordinates at start:",x,y
        WRITE(6,*)"Position at start",p0%point
c       item => start
c       item => item%next
        ib=1
        ibasis = 1
        ix=3
        iy=3
        ix0=rb(ib)%be%ix0(ibasis)
        iy0=rb(ib)%be%iy0(ibasis)
        dx=rb(ib)%be%dx(ibasis)
        dy=rb(ib)%be%dy(ibasis)
        CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
        p0%point(1:2) = rb(ib)%rz%f
        CALL rz_to_cell(p0,item,item_pre,failure,.TRUE.) !failure=.TRUE. is bad
        WRITE(6,*)"After rz_to_cell",failure
        CALL refine_cell(p0,item,failure,x,y,.FALSE.)
        WRITE(6,*)"After refine_cell",failure
        WRITE(6,*)"Logical Coordinates:",x,y
        WRITE(6,*)"Position",p0%point
cTMP
c       CALL open_bin(xdr_unit,'temp.bin','REPLACE','ASIS',64_i4)
c       CALL print_cell_position(item)
c       DO ip=1,item%nodes
c         itemt => item%face(ip)%p
c         WRITE(6,*)" # face",ip,itemt%id
c         CALL print_cell_position(itemt)
c       ENDDO
c       WRITE(xdr_unit)
c       CLOSE(UNIT=xdr_unit)

        CALL nim_stop("dump_reset:  local_debug is true.")
      ENDIF
c-----------------------------------------------------------------------
c     Search
c-----------------------------------------------------------------------
      item => start
c-----------------------------------------------------------------------
c     set limits on the basis type index for different block types,
c     and allocate space to copy point data from the reset file.
c-PRE triangles are linear for now.
c-----------------------------------------------------------------------
      max_basisr=poly_degree**2
      max_basist=1
      ALLOCATE(lbe(3,nmodesc),lve(3,nmodesc),lpr(1,nmodesc),
     $         lpe(1,nmodesc),lnd(1,nmodesc),lco(1,nmodesc),
     $         lte(1,nmodesc),lti(1,nmodesc))
c-----------------------------------------------------------------------
c     loop over basis types and blocks.
c-PRE linear and quad in rblocks only for now.
c-----------------------------------------------------------------------
      ib=1
      ibasis=1
      bb_loop: DO
        IF (ib<=nrbl) THEN
          IF (ibasis>max_basisr) THEN
            ibasis=1
            ib=ib+1
          ENDIF
        ELSE 
          IF (ibasis>max_basist) THEN
            ibasis=1
            ib=ib+1
          ENDIF
        ENDIF
        IF (ib>nbl) EXIT
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       quadrilateral elements.
c-----------------------------------------------------------------------
        IF (ib<=nrbl) THEN
          ix1=rb(ib)%mx
          iy1=rb(ib)%my
          ix0=rb(ib)%be%ix0(ibasis)
          iy0=rb(ib)%be%iy0(ibasis)
          dx=rb(ib)%be%dx(ibasis)
          dy=rb(ib)%be%dy(ibasis)
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       triangular elements.
c-PRE   only linear implemented so far.
c-----------------------------------------------------------------------
        ELSE
          ix0=0
          ix1=tb(ib)%mvert
          iy0=0
          iy1=0
        ENDIF
c-----------------------------------------------------------------------
c       loop over elements for a given block and basis.
c-----------------------------------------------------------------------
        ix_loop: DO ix=ix0,ix1
          iy_loop: DO iy=iy0,iy1
            IF (ib<=nrbl) THEN
              CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
              p0%point(1:2) = rb(ib)%rz%f
            ELSE
              p0%point(1) = tb(ib)%tgeom%xs(ix)
              p0%point(2) = tb(ib)%tgeom%ys(ix)
            ENDIF
            IF(item%id == 0)item => start
            CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
c-TMP           WRITE(18,*) "#"
            CALL refine_cell(p0,item,failure,x,y,.FALSE.)
c-TMP           call bicube_eval(rbc(item%ib)%rz,x,y,0_i4)
c-TMP           WRITE(18,*)p0%point(1),rbc(item%ib)%rz%f(1)
c-TMP           WRITE(18,*)p0%point(2),rbc(item%ib)%rz%f(2)
            IF(failure)THEN
              IF (.NOT. ASSOCIATED (fail_start) ) THEN
                ALLOCATE(fail_start)
                fail_cnt => fail_start
              ELSE
                ALLOCATE(fail_cnt%next)
                fail_cnt => fail_cnt%next
              ENDIF
              fail_cnt%ib=ib
              fail_cnt%ibasis=ibasis
              fail_cnt%ix=ix
              fail_cnt%iy=iy
              NULLIFY(fail_cnt%next)
              WRITE(out_unit,*)" # EXTRAPOLATION at ib,ibasis,ix,iy",
     &          ib,ibasis,ix,iy
              WRITE(out_unit,*)" # EXTRAPOLATION at ",p0%point(1:2)
              data_uninit=.TRUE.
              CYCLE
            ELSE
              IF (item%nodes==0) THEN
                IF (.NOT. ASSOCIATED (fail_start) ) THEN
                  ALLOCATE(fail_start)
                  fail_cnt => fail_start
                ELSE
                  ALLOCATE(fail_cnt%next)
                  fail_cnt => fail_cnt%next
                ENDIF
                fail_cnt%ib=ib
                fail_cnt%ibasis=ibasis
                fail_cnt%ix=ix
                fail_cnt%iy=iy
                NULLIFY(fail_cnt%next)
                WRITE(out_unit,*)" # Exterior point at ",p0%point(1:2)
                data_uninit=.TRUE.
              ELSE
                CALL dumpc_read_data(item)
              ENDIF
            ENDIF
          ENDDO iy_loop
        ENDDO ix_loop
        ibasis=ibasis+1
      ENDDO bb_loop 
c-----------------------------------------------------------------------
c     Extrapolate failed points
c-----------------------------------------------------------------------
      IF(data_uninit.AND.extrap_uninit)THEN
        CALL open_bin(xdr_unit,'failure.bin','REPLACE','ASIS',64_i4)
        fail_cnt => fail_start
        DO WHILE(ASSOCIATED(fail_cnt))
          item => start
          ib=fail_cnt%ib
          ibasis=fail_cnt%ibasis
          ix=fail_cnt%ix
          iy=fail_cnt%iy
          IF (ib<=nrbl) THEN
            ix0=rb(ib)%be%ix0(ibasis)
            iy0=rb(ib)%be%iy0(ibasis)
            dx=rb(ib)%be%dx(ibasis)
            dy=rb(ib)%be%dy(ibasis)
            CALL lagr_quad_eval(rb(ib)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
            p0%point(1:2) = rb(ib)%rz%f
          ELSE
            p0%point(1) = tb(ib)%tgeom%xs(ix)
            p0%point(2) = tb(ib)%tgeom%ys(ix)
          ENDIF
          WRITE(xdr_unit)REAL(p0%point,4)
          WRITE(xdr_unit)
          CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
          CALL print_cell_position(item_pre)
          WRITE(xdr_unit)
          CALL refine_cell(p0,item_pre,failure,x,y,.TRUE.)
c-TMP
c           WRITE(6,*)"Position"
c           WRITE(6,*)"Logical Coordinates:",x,y
c           CALL bicube_eval(rbc(item_pre%ib)%rz,x,y,0_i4)
c           WRITE(6,*)p0%point(1),rbc(item_pre%ib)%rz%f(1)
c           WRITE(6,*)p0%point(2),rbc(item_pre%ib)%rz%f(2)
c-TMP
          CALL dumpc_read_data(item_pre)
          fail_cnt => fail_cnt%next
        ENDDO
        CALL close_bin(xdr_unit,'failure.bin')
      ENDIF
c-----------------------------------------------------------------------
c     Send Warning Message if data uninitialized.
c-----------------------------------------------------------------------
      IF(data_uninit)THEN
         WRITE(nim_wr,*) ' '
         WRITE(nim_wr,*)
     &    "WARNING: The reset_file extrapolated boundary points."
         WRITE(nim_wr,*)
     &    "Check nimset.out for locations, or plot data markers in failu
     &re.bin."
         WRITE(nim_wr,*) ' '
      ENDIF
      RETURN

c-----------------------------------------------------------------------
c     internal subroutine for collecting data.
c-----------------------------------------------------------------------
      CONTAINS

        SUBROUTINE dumpc_read_data(item)

        TYPE(cell_type), POINTER :: item
        REAL(r8) :: rescale = 1.000

        SELECT CASE(item%nodes)
        CASE(3)
          CALL tri_linear_eval(tbc(item%ib)%be,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%ve,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%pres,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%prese,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%nd,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%conc,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%tele,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          CALL tri_linear_eval(tbc(item%ib)%tion,
     &                         tbc(item%ib)%tgeom,
     &              p0%point(1),p0%point(2),item%icell,0_i4)
          lbe=tbc(item%ib)%be%f*rescale
          lve=tbc(item%ib)%ve%f*rescale
          lpr=tbc(item%ib)%pres%f*rescale
          lpe=tbc(item%ib)%prese%f*rescale
          lnd=tbc(item%ib)%nd%f*rescale
          lco=tbc(item%ib)%conc%f*rescale
          lte=tbc(item%ib)%tele%f*rescale
          lti=tbc(item%ib)%tion%f*rescale
        CASE(4)
c                                  Evaluate fields in quad cells.
          CALL lagr_quad_eval(rbc(item%ib)%be   ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%ve   ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%pres ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%prese,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%nd   ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%conc ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%tele ,x,y,0_i4)
          CALL lagr_quad_eval(rbc(item%ib)%tion ,x,y,0_i4)
          lbe=rbc(item%ib)%be%f*rescale
          lve=rbc(item%ib)%ve%f*rescale
ctom         lve=0.
          lpr=rbc(item%ib)%pres%f*rescale
          lpe=rbc(item%ib)%prese%f*rescale
          lnd=rbc(item%ib)%nd%f*rescale
          lco=rbc(item%ib)%conc%f*rescale
          lte=rbc(item%ib)%tele%f*rescale
          lti=rbc(item%ib)%tion%f*rescale
        END SELECT
c-----------------------------------------------------------------------
c       transfer the data to the appropriate basis and Fourier
c       component index.
c-----------------------------------------------------------------------
        DO imode=1,nmodes
          DO jmode=1,nmodesc
            IF ( ABS(keff(imode)-keffc(jmode))/
     $          (ABS(keff(imode))+100*TINY(x))<1.e-11) THEN
              IF (ib<=nrbl) THEN
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%be,lbe(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%ve,lve(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%pres,lpr(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%prese,lpe(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%nd,lnd(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%conc,lco(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%tele,lte(:,jmode),ibasis,ix,iy,imode)
                CALL lagr_quad_basis_assign_loc
     $            (rb(ib)%tion,lti(:,jmode),ibasis,ix,iy,imode)
              ELSE
                tb(ib)%be%fs(:,ix,iy,imode)=lbe(:,jmode)
                tb(ib)%ve%fs(:,ix,iy,imode)=lve(:,jmode)
                tb(ib)%pres%fs(:,ix,iy,imode)=lpr(:,jmode)
                tb(ib)%prese%fs(:,ix,iy,imode)=lpe(:,jmode)
                tb(ib)%nd%fs(:,ix,iy,imode)=lnd(:,jmode)
                tb(ib)%conc%fs(:,ix,iy,imode)=lco(:,jmode)
                tb(ib)%tele%fs(:,ix,iy,imode)=lte(:,jmode)
                tb(ib)%tion%fs(:,ix,iy,imode)=lti(:,jmode)
              ENDIF
            ENDIF
          ENDDO
        ENDDO

        RETURN
        END SUBROUTINE dumpc_read_data


      END SUBROUTINE dumpc_read


      SUBROUTINE rz_to_cell(p0,item,item_pre,failure,debugpath)
c-----------------------------------------------------------------------
c     Routine returns the cell associated with the point p0 based
c     only on the corner points of each cells.  If the routine is
c     unable to find the cell, the failure becomes true.
c-----------------------------------------------------------------------
      USE local
      IMPLICIT NONE
      LOGICAL, INTENT(out) :: failure
      TYPE(location_type), INTENT(IN) :: p0
      LOGICAL, INTENT(in) :: debugpath
      TYPE(cell_type), POINTER :: item,itemt,item_pre
      INTEGER(i4) :: ip,ip1,icnt,face
      INTEGER(i4) :: maxcnt=1000
      REAL(r8) :: etest, emax,reval,zeval,dist,length
      REAL(r8) :: tol_corner = 1.0e-07_r8

      failure=.FALSE.
      maxcnt=maxcell
      
      IF(debugpath)
     &    CALL open_bin(xdr_unit,'path.bin','REPLACE','ASIS',64_i4)
      icnt=0
      itemt => item
      IF(debugpath)WRITE(xdr_unit)REAL(p0%point,4)
      IF(debugpath)WRITE(xdr_unit)
      search_loop : DO WHILE(ASSOCIATED(itemt))
         icnt=icnt+1
         IF(icnt > maxcnt)THEN
            failure=.TRUE.
            item => extbnd
            EXIT search_loop
         ENDIF
         IF(debugpath)CALL print_cell_position(itemt)
         IF(debugpath)WRITE(xdr_unit)
         emax=-HUGE(1)
         DO ip=1,itemt%nodes
            ip1=ip+1
            IF(ip1 > itemt%nodes)ip1=1
            reval=0.5*(itemt%pp(1,ip )+itemt%pp(1,ip1))	! center pt on boundary
            zeval=0.5*(itemt%pp(2,ip )+itemt%pp(2,ip1))
            dist=(p0%point(1)-reval)**2+(p0%point(2)-zeval)**2
            IF(dist < tol_corner) THEN
               item => itemt
               EXIT search_loop
            ENDIF
            length=SQRT((itemt%pp(1,ip1)-itemt%pp(1,ip ))**2
     &                 +(itemt%pp(2,ip1)-itemt%pp(2,ip ))**2)
            IF(length == 0) CYCLE
            etest=(itemt%pp(2,ip1)-itemt%pp(2,ip ))*(reval-p0%point(1))
     &           -(itemt%pp(1,ip1)-itemt%pp(1,ip ))*(zeval-p0%point(2))
            etest=etest*SIGN(1.0_r8,itemt%area)/(dist*length)
            IF(etest > emax)THEN
               face=ip
               emax=etest
            ENDIF
         ENDDO
c
c						Trap for host cell
         IF(emax <= tol_corner)THEN
            item => itemt
            EXIT search_loop
         ENDIF
         item_pre => itemt
         itemt => itemt%face(face)%p
      ENDDO search_loop
      IF(debugpath)CALL close_bin(xdr_unit,'path.bin')
      RETURN
      END SUBROUTINE rz_to_cell
     

      SUBROUTINE make_cell
      IMPLICIT NONE
      INTEGER(i4) :: icell,ib,ix,iy,iglobal,iedge,iedget
      INTEGER(i4) :: iva,ivb,iv,iv1,iv2,ip1,ib1,ix1,iy1,itemp
      INTEGER(i4) :: iba,ibb,iva1,iva2,ip,itempa,ivb2
      TYPE(cell_type),POINTER :: item,itemt
      ic(1,1:2) =( (/1,2/))
      ic(2,1:2) =( (/2,3/))
      ic(3,1:2) =( (/3,1/))
c---------------------------------------------------------------------
c						External boundary cell
c---------------------------------------------------------------------
      NULLIFY(extbnd)
      ALLOCATE(extbnd)
      extbnd%nodes=0
      extbnd%id=0
      extbnd%ib=0
c---------------------------------------------------------------------
c     Build rblock cells
c---------------------------------------------------------------------
c                                                ^  .2    .1=(ib,ix,iy)  
c                                                |    (cell=ix,iy)  
c                                                y  .3    .4
c                                                  x --->
      ALLOCATE(bl2cl(nblc))
      NULLIFY(start,item)
      iglobal=0
      DO ib=1,nrblc
        ALLOCATE(bl2cl(ib)%p(rbc(ib)%mx,rbc(ib)%my))
        DO ix=1,rbc(ib)%mx
          DO iy=1,rbc(ib)%my
            iglobal=iglobal+1
            IF (.NOT. ASSOCIATED (start) ) THEN
               ALLOCATE(start)
               item => start
            ELSE
               ALLOCATE(item%next)
               item => item%next
            ENDIF
            ALLOCATE(bl2cl(ib)%p(ix,iy)%p)
            bl2cl(ib)%p(ix,iy)%p => item
            item%nodes=4
            item%id=iglobal
            item%icell=iy+(ix-1)*rbc(ib)%my
            item%representation='cubic'
            item%degenerate=.FALSE.
            ALLOCATE(item%face(4))
            ALLOCATE(item%p(2,4))
            item%ib=ib
            item%p(1,1)=ix
            item%p(2,1)=iy
            item%p(1,2)=ix
            item%p(2,2)=iy-1
            item%p(1,3)=ix-1
            item%p(2,3)=iy-1
            item%p(1,4)=ix-1
            item%p(2,4)=iy
            NULLIFY (item%next)
          ENDDO
        ENDDO
      ENDDO
c---------------------------------------------------------------------
c     Build tblock cells
c---------------------------------------------------------------------
      DO ib=nrblc+1,nblc
        ALLOCATE(bl2cl(ib)%p(tbc(ib)%mcell,1))
        DO icell=1,tbc(ib)%mcell
            iglobal=iglobal+1
            IF (.NOT. ASSOCIATED (start) ) THEN
               ALLOCATE(start)
               item => start
            ELSE
               ALLOCATE(item%next)
               item => item%next
            ENDIF
            ALLOCATE(bl2cl(ib)%p(icell,1)%p)
            bl2cl(ib)%p(icell,1)%p => item
            item%nodes=3
            item%ib=ib
            item%id=iglobal
            item%icell=icell
            item%representation='liner'
            item%degenerate=.FALSE.
            ALLOCATE(item%p(1,3))
            ALLOCATE(item%face(3))
            item%p(1,1:3)=tbc(ib)%tgeom%vertex(icell,1:3)
            item%face(1)%p => extbnd
            item%face(2)%p => extbnd
            item%face(3)%p => extbnd
        ENDDO
      ENDDO
      maxcell=iglobal
c---------------------------------------------------------------------
c     Internal faces of rblocks
c---------------------------------------------------------------------
c                                                       f4
c                                                ^  .4    .1=(ib,ix,iy)  
c                                           f3   |    (cell=ix,iy)   f1
c                                                y  .3    .2
c                                                  x --->
c                                                       f2
      DO ib=1,nrblc
        DO ix=2,rbc(ib)%mx-1
        DO iy=2,rbc(ib)%my-1
           item => bl2cl(ib)%p(ix,iy)%p
           item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
           item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
           item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
           item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 
        ENDDO
        ENDDO
c					Internal part of cells with seam
        ix=1
        iy=1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
          item%face(2)%p => extbnd
          item%face(3)%p => extbnd
          item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 
        DO iy=2,rbc(ib)%my-1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
          item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
          item%face(3)%p => extbnd
          item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 
        ENDDO
        iy=rbc(ib)%my
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
          item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
          item%face(3)%p => extbnd
          item%face(4)%p => extbnd

        DO ix=2,rbc(ib)%mx-1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
          item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
          item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
          item%face(4)%p => extbnd
        ENDDO

        ix=rbc(ib)%mx
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => extbnd
          item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
          item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
          item%face(4)%p => extbnd
        DO iy=2,rbc(ib)%my-1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => extbnd
          item%face(2)%p => bl2cl(ib)%p(ix  ,iy-1)%p 
          item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
          item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 
        ENDDO
        iy=1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => extbnd
          item%face(2)%p => extbnd
          item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
          item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 

        DO ix=2,rbc(ib)%mx-1
          item => bl2cl(ib)%p(ix,iy)%p
          item%face(1)%p => bl2cl(ib)%p(ix+1,iy  )%p 
          item%face(2)%p => extbnd
          item%face(3)%p => bl2cl(ib)%p(ix-1,iy  )%p 
          item%face(4)%p => bl2cl(ib)%p(ix  ,iy+1)%p 
        ENDDO
      ENDDO
c---------------------------------------------------------------------
c     tblock internal faces
c---------------------------------------------------------------------
      DO ib=nrblc+1,nblc
        DO icell=1,tbc(ib)%mcell
          item => bl2cl(ib)%p(icell,1)%p
c
c                                       	Loop through vertex pairs
          vert: DO iv=1,3
             iva=tbc(ib)%tgeom%vertex(icell,ic(iv,1))
             ivb=tbc(ib)%tgeom%vertex(icell,ic(iv,2))
c
c                                       Search for matching vertices.
             DO ix=1,tbc(ib)%mcell
               IF (ix == icell) CYCLE
               DO iv1=1,3
                 iva2 = tbc(ib)%tgeom%vertex(ix,ic(iv1,1))
                 ivb2 = tbc(ib)%tgeom%vertex(ix,ic(iv1,2))
                 IF(((iva == iva2) .AND. (ivb == ivb2)) .OR.
     &              ((iva == ivb2) .AND. (ivb == iva2))) THEN
                   item%face(iv)%p => bl2cl(ib)%p(ix,1)%p
                   CYCLE vert
                 ENDIF
               ENDDO
             ENDDO
          ENDDO vert
        ENDDO
      ENDDO
c
c						seam based faces
      DO ib=1,nblc
        DO iv1=1,seamc(ib)%nvert
          iv2=iv1+1
          IF(iv2 > seamc(ib)%nvert)iv2=1
c						convert seam pair to a cell
          CALL seam_to_cell(ib,iv1,iv2,itemp,item)
c         WRITE(6,fmt='(x,"local",8(x,i6))') ib,iv1,iv2,item%id,itemp
c						convert seam pair to seam pair
c						across the seam with same block
          seamtoseam : DO ip=1,SIZE(seamc(ib)%vertex(iv1)%ptr,2)
c            WRITE(6,*)"ip ",ip,
c    &                 seamc(ib)%vertex(iv1)%ptr(1,ip),
c    &                 seamc(ib)%vertex(iv1)%ptr(2,ip)
             DO ip1=1,SIZE(seamc(ib)%vertex(iv2)%ptr,2)
c            WRITE(6,*)"ip1",ip1,
c    &                 seamc(ib)%vertex(iv2)%ptr(1,ip1),
c    &                 seamc(ib)%vertex(iv2)%ptr(2,ip1)
c
c						compare block labels first.
                iba=seamc(ib)%vertex(iv1)%ptr(1,ip )
                ibb=seamc(ib)%vertex(iv2)%ptr(1,ip1)
                IF(iba == ibb) THEN
                   iva1=seamc(ib)%vertex(iv1)%ptr(2,ip )
                   iva2=seamc(ib)%vertex(iv2)%ptr(2,ip1)
c
c						check the seam values.
c                  write(6,*)iva1,iva2
                   IF (iba == 0) EXIT seamtoseam
                   IF(ABS(iva1-iva2) == 1) EXIT seamtoseam
                   IF(((iva1== 1) .AND. (iva2 == seamc(iba)%nvert))
     &               .OR.
     &                ((iva2== 1) .AND. (iva1 == seamc(iba)%nvert)))
     &             EXIT seamtoseam
                ENDIF
             ENDDO
          ENDDO seamtoseam
c
c						convert seam pair to a cell
          CALL seam_to_cell(iba,iva1,iva2,itempa,itemt)
c         WRITE(6,fmt='(8(x,i8))') iba,iva1,iva2,itemt%id,itempa
c         WRITE(6,*)
          item%face(itemp)%p => itemt
c
c						Special Fixups
          IF((iba == 0) .AND. (iva1 == iva2))THEN
c           WRITE(6,fmt='(x,"local",8(x,i6))') ib,iv1,iv2,item%id,itemp
c           WRITE(6,fmt='(8(x,i8))') iba,iva1,iva2,itemt%id,itempa
c
c							Degenerate
            IF(itemp == 3) THEN
              item%degenerate=.TRUE.
              IF(itemp == item%nodes)THEN
                item%face(itemp)%p => item%face(1)%p
              ELSE
                item%face(itemp)%p => item%face(itemp+1)%p
              ENDIF
            ELSE
c
c							Boundary
              item%face(itemp)%p => extbnd
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE make_cell


      SUBROUTINE map_cell_position(ptr)
c-----------------------------------------------------------------------
c     Routine adds the actual physical points of the corners to each cell
c     and determines if the cell order is up or down.
c-----------------------------------------------------------------------
      TYPE(cell_type), POINTER :: ptr
      INTEGER(i4) :: ip,ib,ix,iy,iv
      REAL(r8) :: g
c
c						assign values to nodes.
      ib=ptr%ib
      IF(ib == 0 )THEN
         RETURN
      ELSE IF (ib <= nrblc) THEN
        ALLOCATE(ptr%pp(2,ptr%nodes))
        DO ip=1,ptr%nodes
          ix=ptr%p(1,ip)
          iy=ptr%p(2,ip)
          ptr%pp(1,ip) = rbc(ib)%rz%fs(1,ix,iy)
          ptr%pp(2,ip) = rbc(ib)%rz%fs(2,ix,iy)
        ENDDO
      ELSE 
        ALLOCATE(ptr%pp(2,ptr%nodes))
        DO ip=1,ptr%nodes
          iv=ptr%p(1,ip)
          ptr%pp(1,ip) = tbc(ib)%tgeom%xs(iv)
          ptr%pp(2,ip) = tbc(ib)%tgeom%ys(iv)
        ENDDO
      ENDIF
c
c						compute area.
      g=0
      DO ip=1,ptr%nodes-1
        g=g+ptr%pp(1,ip+1)*ptr%pp(2,ip)
     &     -ptr%pp(1,ip  )*ptr%pp(2,ip)
     &  +0.5*(ptr%pp(1,ip+1)-ptr%pp(1,ip))*(ptr%pp(2,ip+1)-ptr%pp(2,ip))
      ENDDO
      ip=ptr%nodes
      g=g+ptr%pp(1,   1)*ptr%pp(2,ip)
     &   -ptr%pp(1,ip  )*ptr%pp(2,ip)
     &  +0.5*(ptr%pp(1,   1)-ptr%pp(1,ip))*(ptr%pp(2,   1)-ptr%pp(2,ip))
      ptr%area = g
      RETURN
      END SUBROUTINE map_cell_position


      SUBROUTINE print_cell_position(ptr)
c-----------------------------------------------------------------------
c     Create an xdraw entry that contains the cell and neighbors.
c-----------------------------------------------------------------------
      TYPE(cell_type), POINTER :: ptr,item
      INTEGER(i4) :: ip,ib,ix,iy,iv
      ib=ptr%ib
      IF(ib == 0 )RETURN
      DO ip=1,ptr%nodes
        WRITE(xdr_unit)REAL(ptr%pp(1,ip),4),REAL(ptr%pp(2,ip),4)
      ENDDO
      ip=1
      WRITE(xdr_unit)REAL(ptr%pp(1,ip),4),REAL(ptr%pp(2,ip),4)
      RETURN
      END SUBROUTINE print_cell_position


      SUBROUTINE seam_to_cell(ib,iv1,iv2,face,ptr)
c-----------------------------------------------------------------------
c     Return the cell ptr for the two seam points (ib,iv1) to (ib,iv2).
c-----------------------------------------------------------------------
      TYPE(cell_type), POINTER :: ptr,item1,item2,itema,itemb
      INTEGER(i4), INTENT(in) :: ib,iv1,iv2
      INTEGER(i4), INTENT(out) :: face
      INTEGER(i4) :: mx,my,ix,iy
      INTEGER(i4) :: iv,icell,iva2,ivb2,iva,ivb

c
c				Exterior boundary.
      IF(ib == 0)THEN
        face=0
        ptr => extbnd
c
c				Determine common cell of iv1 and iv2 in tblocks.
      ELSE IF (ib > nrblc) THEN
         iva=seamc(ib)%vertex(iv1)%intxy(1)
         ivb=seamc(ib)%vertex(iv2)%intxy(1)
c        WRITE(6,*)" #",iva,ivb
         DO icell=1,tbc(ib)%mcell
               DO iv=1,3
                 iva2 = tbc(ib)%tgeom%vertex(icell,ic(iv,1))
                 ivb2 = tbc(ib)%tgeom%vertex(icell,ic(iv,2))
                 IF(((iva == iva2) .AND. (ivb == ivb2)) .OR.
     &              ((iva == ivb2) .AND. (ivb == iva2))) THEN
                   ptr => bl2cl(ib)%p(icell,1)%p
                   face=iv
                   RETURN
                 ENDIF
               ENDDO
         ENDDO
c        ptr => extbnd
         CALL nim_stop("seam_to_cell:  Could not match tblock")
      ELSE
c
c				Determine possible cells of iv1 in rblocks.
        mx=rbc(ib)%mx
        my=rbc(ib)%my
        IF(iv1 < mx)THEN
         iy=0
         ix=iv1
         item1 => bl2cl(ib)%p(ix  ,iy+1)%p
         item2 => bl2cl(ib)%p(ix+1,iy+1)%p
         face=2
        ELSEIF(iv1 == mx)THEN
         iy=0
         ix=mx
         item1 => bl2cl(ib)%p(ix,iy+1)%p
         item2 => item1
         IF(iv2 < iv1)THEN
           face=2
         ELSE
           face=1
         ENDIF
        ELSEIF(iv1 < my+mx)THEN
         ix=mx
         iy=iv1-mx
         item1 => bl2cl(ib)%p(ix,iy  )%p
         item2 => bl2cl(ib)%p(ix,iy+1)%p
         face=1
        ELSEIF(iv1 == my+mx)THEN
         ix=mx
         iy=my
         item1 => bl2cl(ib)%p(ix  ,iy)%p
         item2 => item1
         IF(iv2 < iv1)THEN
           face=1
         ELSE
           face=4
         ENDIF
        ELSEIF(iv1 < 2*mx+my)THEN
         ix=2*mx+my-iv1
         iy=my
         item1 => bl2cl(ib)%p(ix  ,iy)%p
         item2 => bl2cl(ib)%p(ix+1,iy)%p
         face=4
        ELSEIF(iv1 == 2*mx+my)THEN
         ix=0
         iy=my
         item1 => bl2cl(ib)%p(ix+1,iy)%p
         item2 => item1
         IF(iv2 < iv1)THEN
           face=4
         ELSE
           face=3
         ENDIF
        ELSEIF(iv1 < (2*(my+mx)))THEN
         ix=0
         iy=2*(mx+my)-iv1
         item1 => bl2cl(ib)%p(ix+1,iy  )%p
         item2 => bl2cl(ib)%p(ix+1,iy+1)%p
         face=3
        ELSEIF(iv1 == 2*(my+mx))THEN
         ix=0
         iy=0
         item1 => bl2cl(ib)%p(ix+1,iy+1)%p
         item2 => item1
         IF(ABS(iv2-iv1) > 1)THEN
           face=2
         ELSE
           face=3
         ENDIF
        ELSE
         CALL nim_stop("seam_to_cell: invalid iv1")
        ENDIF
c
c					Determine possible cells of iv2
        IF(iv2 < mx)THEN
         iy=0
         ix=iv2
         itema => bl2cl(ib)%p(ix  ,iy+1)%p
         itemb => bl2cl(ib)%p(ix+1,iy+1)%p
        ELSEIF(iv2 == mx)THEN
         iy=0
         ix=mx
         itema => bl2cl(ib)%p(ix,iy+1)%p
         itemb => itema
        ELSEIF(iv2 < my+mx)THEN
         ix=mx
         iy=iv2-mx
         itema => bl2cl(ib)%p(ix,iy  )%p
         itemb => bl2cl(ib)%p(ix,iy+1)%p
        ELSEIF(iv2 == my+mx)THEN
         ix=mx
         iy=my
         itema => bl2cl(ib)%p(ix  ,iy)%p
         itemb => itema
        ELSEIF(iv2 < 2*mx+my)THEN
         ix=2*mx+my-iv2
         iy=my
         itema => bl2cl(ib)%p(ix  ,iy)%p
         itemb => bl2cl(ib)%p(ix+1,iy)%p
        ELSEIF(iv2 == 2*mx+my)THEN
         ix=0
         iy=my
         itema => bl2cl(ib)%p(ix+1,iy)%p
         itemb => itema
        ELSEIF(iv2 < (2*(my+mx)))THEN
         ix=0
         iy=2*(mx+my)-iv2
         itema => bl2cl(ib)%p(ix+1,iy  )%p
         itemb => bl2cl(ib)%p(ix+1,iy+1)%p
        ELSEIF(iv2 == 2*(my+mx))THEN
         ix=0
         iy=0
         itema => bl2cl(ib)%p(ix+1,iy+1)%p
         itemb => item1
        ELSE
         CALL nim_stop("seam_to_cell: invalid iv2")
        ENDIF
c
c						Return the common cell.
        IF((item1%id == itema%id).OR.(item1%id == itemb%id))THEN
         ptr => item1
        ELSEIF((item2%id == itema%id).OR.(item2%id == itemb%id))THEN
         ptr => item2
        ELSE
         face=0
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE seam_to_cell


      SUBROUTINE refine_cell(p0,item,failure,x,y,remain_local)

      IMPLICIT NONE

      INTEGER(i4) :: ib,icell
      TYPE(location_type), INTENT(in) :: p0
      TYPE(cell_type), POINTER :: item,itemt
      LOGICAL, INTENT(inout) :: failure
      LOGICAL, INTENT(in) :: remain_local
      REAL(r8),INTENT(out) :: x,y
      REAL(r8) :: tol=1.0e-12
      REAL(r8) :: xcell,ycell,denom
      REAL(r8) :: dr,dz,dx,dy,err,jacobian,xold,yold,slow_down
      INTEGER(i4) :: max_loop = 1000
      INTEGER(i4) :: iter_loop
      REAL(r8) :: slow_down_frac=0.2
 

      SELECT CASE(item%nodes)
      CASE(0)					! Exterior boundary
        RETURN
      CASE(3)					! Triangles
        SELECT CASE(item%representation)
        CASE('liner')
           RETURN
        CASE('pcnst')
           CALL nim_stop
     &      ("refine_cell: pcnst is unimplemented for triangles")
        CASE('cubic')
           CALL nim_stop
     &      ("refine_cell: cubic is unimplemented for triangles")
        CASE default
           CALL nim_stop
     &      ("refine_cell: unrecognized 'represenation' for triangles")
        END SELECT

      CASE(4)					! quadrilateral cells
        SELECT CASE(item%representation)
        CASE('liner')
           CALL nim_stop
     &      ("refine_cell: liner is unimplemented for quads")
        CASE('pcnst')
           CALL nim_stop
     &      ("refine_cell: pcnst is unimplemented for quads")
        CASE('cubic')
c         denom=(item%pp(1,1)-item%pp(1,2))*(item%pp(2,1)-item%pp(2,4))
c    &         +(item%pp(1,1)-item%pp(1,4))*(item%pp(2,2)-item%pp(2,1))
c         xcell=((p0%point(1) -item%pp(1,2))*(item%pp(2,1)-item%pp(2,4))
c    &         +(item%pp(1,1)-item%pp(1,4))*(item%pp(2,1)-p0%point(2)))
c    &         /denom
c         ycell=((p0%point(1) -item%pp(1,2))*(item%pp(2,2)-item%pp(2,1))
c    &         +(item%pp(1,2)-item%pp(1,1))*(item%pp(2,1)-p0%point(2)))
c    &         /denom
          iter_loop = 0_i4
          xcell=-0.5
          ycell=-0.5
          refine: do
            x=REAL(item%p(1,1))+xcell
            y=REAL(item%p(2,1))+ycell
            call lagr_quad_eval(rbc(item%ib)%rz,x,y,1_i4)
c
c					 	Convergence test of dr and dz
            dr=p0%point(1)-rbc(item%ib)%rz%f(1)
            dz=p0%point(2)-rbc(item%ib)%rz%f(2)
            err=ABS(dr)+ABS(dz)
c           write(6,*)iter_loop,err
c-TMP           IF(failure)THEN
c-TMP             WRITE(18,fmt='(" LOOP: ",i4,x,i4,x,e12.5)')
c-TMP    &              iter_loop,item%id,err
c-TMP             WRITE(18,*)p0%point(1),rbc(item%ib)%rz%f(1)
c-TMP             WRITE(18,*)p0%point(2),rbc(item%ib)%rz%f(2)
c-TMP           ENDIF
            IF(err < tol)EXIT refine
c
c						limit the number of iterations.
            iter_loop=iter_loop+1
            IF(iter_loop > max_loop)THEN
              failure=.TRUE.
c-TMP             CALL print_cell_position(item)
c-TMP             WRITE(xdr_unit)
              RETURN
            ENDIF
c-TMP           IF(iter_loop > max_loop+10) THEN
c-TMP             failure=.TRUE.
c-TMP             RETURN
c-TMP           ENDIF
c
c						Update x and y
            jacobian=rbc(item%ib)%rz%fx(1)*rbc(item%ib)%rz%fy(2)
     &              -rbc(item%ib)%rz%fy(1)*rbc(item%ib)%rz%fx(2)
            IF(jacobian /= 0)THEN
              jacobian=1.0/jacobian
              IF (iter_loop<=4) jacobian=jacobian/2**(5-iter_loop)

              dx = (rbc(item%ib)%rz%fy(2) * dr 
     &            - rbc(item%ib)%rz%fy(1) * dz)
     &             *jacobian
              dy = (rbc(item%ib)%rz%fx(1) * dz 
     &            - rbc(item%ib)%rz%fx(2) * dr)
     &             *jacobian
            ELSE
c						send the solution off in
c						a random direction.
              dx = (ABS(dr)+ABS(dz))
              dy = dx
            ENDIF
            IF(ABS(dx) > 1.0)dx=SIGN(0.5_r8,dx)
            IF(ABS(dy) > 1.0)dy=SIGN(0.5_r8,dy)
            itemt => item
            xold=xcell
            yold=ycell
            slow_down=ABS(REAL(max_loop-iter_loop))/REAL(max_loop)
            xcell=xcell+dx*slow_down*slow_down_frac
            ycell=ycell+dy*slow_down*slow_down_frac

           IF(.NOT.remain_local)THEN
            IF(xcell > 0) THEN
              item => item%face(1)%p
              IF(item%ib > nrblc)THEN
                item => itemt
                xcell = 0
              ELSE
                xcell=xcell-1.0
              ENDIF
            ELSEIF(xcell < -1.)THEN
              item => item%face(3)%p
              IF(item%ib > nrblc)THEN
                item => itemt
                xcell = -1.
              ELSE
                xcell=xcell+1.0
              ENDIF
            ENDIF
            IF(item%id == 0 ) THEN
               failure = .TRUE.
               RETURN
            ENDIF
            IF(ycell > 0) THEN
              item => item%face(4)%p
              IF(item%ib > nrblc)THEN
                item => itemt
                ycell = 0
              ELSE
                ycell=ycell-1.0
              ENDIF
            ELSEIF(ycell < -1.)THEN
              item => item%face(2)%p
              IF(item%ib > nrblc)THEN
                item => itemt
                ycell = -1.0
              ELSE
                ycell=ycell+1.0
              ENDIF
            ENDIF
            IF(item%id == 0 ) THEN
               failure = .TRUE.
               RETURN
            ENDIF
           ENDIF
c
c       Note: Usually a point should not refine out of the rblock.
c       (This may change with other than bilinear triangles!)
          ENDDO refine
        CASE default
           CALL nim_stop
     &      ("refine_cell: unrecognized 'represenation' for quads")
        END SELECT
      CASE default
         WRITE(nim_wr,*)" CELL",item%id," in block ",item%ib
         WRITE(nim_wr,*)" NODES",item%nodes
         CALL nim_stop("refine_cell: Invalid number of nodes.")
      END SELECT
      RETURN
      END SUBROUTINE refine_cell

      SUBROUTINE make_search_map
      IMPLICIT NONE
      TYPE(cell_type), POINTER :: item,item_pre
      TYPE(location_type):: p0
      INTEGER(i4) :: ip,num_cell,ix,iy,ix_test,iy_test
      LOGICAL :: failure
      REAL(r8) :: x,y
      INTEGER(i4),DIMENSION(8) :: ix_offset,iy_offset

c
c						Biggest and smallest R and Z.
      search_rmin = HUGE(1.0)
      search_rmax = -HUGE(1.0)
      search_zmin = HUGE(1.0)
      search_zmax = -HUGE(1.0)
      num_cell = 0
      item => start
      DO WHILE(ASSOCIATED(item))
        CALL map_cell_position(item)
        DO ip=1,item%nodes
           IF(item%pp(1,ip) > search_rmax)search_rmax=item%pp(1,ip)
           IF(item%pp(1,ip) < search_rmin)search_rmin=item%pp(1,ip)
           IF(item%pp(2,ip) > search_zmax)search_zmax=item%pp(2,ip)
           IF(item%pp(2,ip) < search_zmin)search_zmin=item%pp(2,ip)
        ENDDO
        num_cell=num_cell + 1
        item => item%next
      ENDDO
c
c						Determine search map dimension

      search_idim  = NINT(SQRT(REAL(num_cell,r8)))
      search_jdim  = search_idim
      ALLOCATE(search_map(0:search_idim,0:search_jdim))

c
c						Determine map to cell info.
      DO ix=0,search_idim
      DO iy=0,search_jdim
      item => start
      item_pre => start
        p0%point(1) = search_rmin 
     &    + REAL(ix,r8)*(search_rmax-search_rmin)/REAL(search_idim,r8)
        p0%point(2) = search_zmin
     &    + REAL(iy,r8)*(search_zmax-search_zmin)/REAL(search_jdim,r8)
        CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
        CALL refine_cell(p0,item,failure,x,y,.FALSE.)
        search_map(ix,iy)%next => item
      ENDDO
      ENDDO
c
c						Reassign cells on the 
c						computational boundary to an 
c						interior cell.
      ix_offset=(/-1, 0, 1,1,1,0,-1,0/)
      iy_offset=(/-1,-1,-1,0,1,1, 1,0/)
      DO ix=0,search_idim
      DO iy=0,search_jdim
        item => search_map(ix,iy)%next
        IF(item%ib == 0 ) THEN
          DO ip=1,8
            ix_test = ix+ix_offset(ip)
            IF(ix_test < 0) ix_test = 0
            IF(ix_test > search_idim) ix_test = search_idim
            iy_test = iy+iy_offset(ip)
            IF(iy_test < 0) iy_test = 0
            IF(iy_test > search_jdim) iy_test = search_jdim
            item_pre => search_map(ix_test,iy_test)%next
            IF(item_pre%ib .NE. 0) THEN
               search_map(ix,iy)%next => item_pre
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     DEBUG code:  This is used to plot each cell.
c-----------------------------------------------------------------------
c     CALL open_bin(xdr_unit,'temp.bin','REPLACE','ASIS',64_i4)
c     DO ix=0,search_idim
c     DO iy=0,search_jdim
c       p0%point(1) = search_rmin 
c    &    + REAL(ix,r8)*(search_rmax-search_rmin)/REAL(search_idim,r8)
c       p0%point(2) = search_zmin
c    &    + REAL(iy,r8)*(search_zmax-search_zmin)/REAL(search_jdim,r8)
c       WRITE(xdr_unit)REAL(p0%point,4)
c       WRITE(xdr_unit)
c       item => search_map(ix,iy)%next
c       IF(item%id == 0) cycle
c       CALL print_cell_position(item)
c       WRITE(xdr_unit)
c     ENDDO
c     ENDDO
c     CLOSE(UNIT=xdr_unit)
c-----------------------------------------------------------------------
      END SUBROUTINE make_search_map

c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE dumpc
