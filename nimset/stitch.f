c-----------------------------------------------------------------------
c     program stitch.
c     nimset-like initialization program for multiple regions.  each
c     region is defined with standard nimset operations and then
c     the regions are stitched together.
c
c     there are two important limitations to the present version of
c     stitch.  first, it only handles rectangular gridshapes.  having
c     stitch add external blocks to an existing polar region should be
c     a straightforward modification.  the second limitation is that
c     stitch does not work with (and cannot create) periodic meshes
c     where the external border has completely separate pieces.
c
c     CHANGE LOG
c
c     assembly of the composite seam0 structure has been changed to
c     remove any duplications in the ptr (block-seam pointer) arrays.
c     also, dangling block-seam connections are removed as a new final
c     step to avoid problems after a zipper-like stitch to mend a
c     partial cut left after other stitching operations.
c       CRS, 9/22/09
c
c     the stitch_seams_polar routine and modifications to the driver
c     now allow stitching of complete annular regions.  it can be used
c     to assemble plasma and internal-vacuum regions with distinct
c     properties.  it also works on rectangular regions with y-
c     direction periodicity.  another change is that dump resets are
c     called region-by-region and read equilibrium information from
c     each reset file (specified by reset_file in the nimrod.in for
c     that region).
c       KJB and CRS, 7/18/14
c
c     stitch has been fixed for performing three or more stitchings
c     operations that meet at a vertex.  previously, the seam0 pointer
c     information was not updated in all of the linked list entries, so
c     all connections did not appear in the final seam0.  similarly, the
c     "next" pointer was not connected in all linked list entries.
c     finally, the trimming of seam0 connections from blocks is now
c     using the correct indexing for the new seam0 after all stitching.
c       CRS and KJB, 4/22/15
c     
c     revisions for creating annular regions (as a result of stitching
c     aperiodic regions) are also in this new version.  a
c     check for active region-seam0 vertices that are not in the
c     composite seam0 list identifies stranded sections that are then
c     spliced into the composite.
c       KJB and CRS, 5/1/15
c
c     a new version of stitch_seams_polar has been installed.  it is
c     only used for assembling annular regions, or a disk and annular
c     regions, so it can use stitching operations that differ
c     significantly from those in stitch_seams.  this version now
c     works when nybl>1 in the separate regions.  it also appears to do
c     the right thing when nybl differs among regions; nimrod runs but
c     is not able to create slice plots.
c       CRS, 12/28/16

c     additional pathological cases when creating annular regions
c     prompted further changes.  this type of operation is now handled
c     as a special case that is flagged when a user sets the number
c     of regions to 1.  the "in linked listed" check is no longer used
c     and has been removed from the end of the stitch_seams routine.
c       CRS, 5/20/18

c     extra connections at the beginning and end of a stitch have been
c     added to better handle multiple operations that stitch over
c     previously stitched vertices.
c       CRS, 5/21/18
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. stitch_mod
c     0.1. stitch_ll_init
c     0.2. stitch_bl_adj
c     1. stitch
c     2. write_header
c     3. draw_seam0
c     4. stitch_seams
c     5. stitch_seams_polar
c     6. reset_input
c-----------------------------------------------------------------------
c     0.  stitch_mod module helps avoid duplication between stitch_seams
c     and stitch_seams_polar.
c-----------------------------------------------------------------------
      MODULE stitch_mod
      USE local
      IMPLICIT NONE

      TYPE :: list_type
        TYPE(list_type), POINTER :: prev
        TYPE(list_type), POINTER :: next
        TYPE(list_type), POINTER :: prevorig
        TYPE(list_type), POINTER :: nextorig
        INTEGER(i4) :: iv,ir
        INTEGER(i4), DIMENSION(:,:), POINTER :: ptr
        LOGICAL :: active,excorner,stflag
      END TYPE list_type

      CONTAINS

c-----------------------------------------------------------------------
c     0.1. module subprogram stitch_ll_init initializes the linked
c     lists that are used for combining seam0 structures.
c-----------------------------------------------------------------------
      SUBROUTINE stitch_ll_init(reg,s0starts)
      USE region_type_mod

      TYPE(region_type), DIMENSION(:), INTENT(INOUT) :: reg
      TYPE(list_type), DIMENSION(:), INTENT(INOUT), POINTER :: s0starts

      TYPE(list_type), POINTER :: current,next
      INTEGER(i4) :: ireg,nreg,iv,nv,np
c-----------------------------------------------------------------------
c     duplicate the seam0 information for each region into
c     linked lists that can be assembled relatively easily.
c-----------------------------------------------------------------------
      nreg=SIZE(reg)
      DO ireg=1,nreg
        nv=reg(ireg)%seam0%nvert
        current=>s0starts(ireg)
        DO iv=1,nv-1
          current%iv=iv
          current%ir=ireg
          current%active=.true.
          current%excorner=reg(ireg)%seam0%excorner(iv)
          current%stflag=.false.
          np=SIZE(reg(ireg)%seam0%vertex(iv)%ptr,2)
          ALLOCATE(current%ptr(2,np))
          current%ptr(:,:)=reg(ireg)%seam0%vertex(iv)%ptr(:,:)
          ALLOCATE(current%next)
          next=>current%next
          current%nextorig=>next
          next%prev=>current
          next%prevorig=>current
          current=>next
        ENDDO
        current%iv=nv
        current%ir=ireg
        current%active=.true.
        current%excorner=reg(ireg)%seam0%excorner(nv)
        current%stflag=.false.
        np=SIZE(reg(ireg)%seam0%vertex(nv)%ptr,2)
        ALLOCATE(current%ptr(2,np))
        current%ptr(:,:)=reg(ireg)%seam0%vertex(nv)%ptr(:,:)
        current%next=>s0starts(ireg)
        current%nextorig=>s0starts(ireg)
        s0starts(ireg)%prev=>current
        s0starts(ireg)%prevorig=>current
      ENDDO

      RETURN
      END SUBROUTINE stitch_ll_init
c-----------------------------------------------------------------------
c     0.2. module subprogram stitch_bl_adj combines pointer information
c     in block-seams at a seam0 vertex being stitched together.    
c-----------------------------------------------------------------------
      SUBROUTINE stitch_bl_adj(reg,ir0,iv0,ir1,iv1)
      USE region_type_mod

      TYPE(region_type), DIMENSION(:), INTENT(INOUT) :: reg
      INTEGER(i4), INTENT(IN) :: ir0,iv0,ir1,iv1

      INTEGER(i4) :: np0,np1,ip0,ip1,np,ip,ib,iv,jb
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ptmp
c-----------------------------------------------------------------------
c     identify the pointer structures from each region and their sizes.
c-PRE
c             IF (ib<=nrblt) THEN
c               jb=ib-reg(irout)%irst+1
c             ELSE
c               jb=ib-reg(irout)%itst+1
c             ENDIF
c-----------------------------------------------------------------------
      np0=SIZE(reg(ir0)%seam0%vertex(iv0)%ptr,2)
      np1=SIZE(reg(ir1)%seam0%vertex(iv1)%ptr,2)

      DO ip0=1,np0
        ib=reg(ir0)%seam0%vertex(iv0)%ptr(1,ip0)
        iv=reg(ir0)%seam0%vertex(iv0)%ptr(2,ip0)
        jb=ib-reg(ir0)%irst+1
        np=0
        DO ip=1,SIZE(reg(ir0)%seam(jb)%vertex(iv)%ptr,2)
          IF (reg(ir0)%seam(jb)%vertex(iv)%ptr(1,ip)>0) np=np+1
        ENDDO
        ALLOCATE(ptmp(2,np))
        np=0
        DO ip=1,SIZE(reg(ir0)%seam(jb)%vertex(iv)%ptr,2)
          IF (reg(ir0)%seam(jb)%vertex(iv)%ptr(1,ip)>0) THEN
            np=np+1
            ptmp(:,np)=reg(ir0)%seam(jb)%vertex(iv)%ptr(:,ip)
          ENDIF
        ENDDO
        DEALLOCATE(reg(ir0)%seam(jb)%vertex(iv)%ptr)
        ALLOCATE(reg(ir0)%seam(jb)%vertex(iv)%ptr(2,np+np1))
        reg(ir0)%seam(jb)%vertex(iv)%ptr(:,1:np)=ptmp
        reg(ir0)%seam(jb)%vertex(iv)%ptr(:,np+1:)=
     $    reg(ir1)%seam0%vertex(iv1)%ptr(:,:)
        DEALLOCATE(ptmp)
      ENDDO

      DO ip1=1,np1
        ib=reg(ir1)%seam0%vertex(iv1)%ptr(1,ip1)
        iv=reg(ir1)%seam0%vertex(iv1)%ptr(2,ip1)
        jb=ib-reg(ir1)%irst+1
        np=0
        DO ip=1,SIZE(reg(ir1)%seam(jb)%vertex(iv)%ptr,2)
          IF (reg(ir1)%seam(jb)%vertex(iv)%ptr(1,ip)>0) np=np+1
        ENDDO
        ALLOCATE(ptmp(2,np))
        np=0
        DO ip=1,SIZE(reg(ir1)%seam(jb)%vertex(iv)%ptr,2)
          IF (reg(ir1)%seam(jb)%vertex(iv)%ptr(1,ip)>0) THEN
            np=np+1
            ptmp(:,np)=reg(ir1)%seam(jb)%vertex(iv)%ptr(:,ip)
          ENDIF
        ENDDO
        DEALLOCATE(reg(ir1)%seam(jb)%vertex(iv)%ptr)
        ALLOCATE(reg(ir1)%seam(jb)%vertex(iv)%ptr(2,np+np0))
        reg(ir1)%seam(jb)%vertex(iv)%ptr(:,1:np)=ptmp
        reg(ir1)%seam(jb)%vertex(iv)%ptr(:,np+1:)=
     $    reg(ir0)%seam0%vertex(iv0)%ptr(:,:)
        DEALLOCATE(ptmp)
      ENDDO

      reg(ir0)%seam0%excorner(iv0)=.false.
      reg(ir1)%seam0%excorner(iv1)=.false.

      RETURN
      END SUBROUTINE stitch_bl_adj

      END MODULE stitch_mod

c-----------------------------------------------------------------------
c     1.  main program, stitch
c-----------------------------------------------------------------------
      PROGRAM stitch
      USE local
      USE input
      USE fields
      USE dump
      USE seam_storage_mod
      USE nimset_init
      USE polar_init
      USE diagnose
      USE region_type_mod
      IMPLICIT NONE

      INTEGER(i4) :: istep=0,nmodes,ib,in,nv,ip,np,iv,ibg,ibgt,
     $               nstitch,jv,jb,jp,np1,iv0
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ptmp,count
      LOGICAL :: file_stat,allper=.true.
      REAL(r8):: t=0
      REAL(r8), DIMENSION(:), POINTER :: keff
      CHARACTER(32), DIMENSION(:), ALLOCATABLE :: infile
      CHARACTER(128) :: file_save

      TYPE(region_type), DIMENSION(:), POINTER :: region
      INTEGER(i4) :: ireg,nreg
c-----------------------------------------------------------------------
c     interface block for draw_seam0.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE draw_seam0(nreg,nrblt,nblt,reg)
        USE local
        USE input
        USE region_type_mod
        IMPLICIT NONE
        INTEGER(i4), INTENT(IN) :: nreg,nrblt,nblt
        TYPE(region_type), DIMENSION(:), POINTER :: reg
        END SUBROUTINE draw_seam0
      END INTERFACE
c-----------------------------------------------------------------------
c     nullify pointers.
c-----------------------------------------------------------------------
      NULLIFY(rb)
      NULLIFY(tb)
      NULLIFY(seam0)
      NULLIFY(seam)
      NULLIFY(region)
      OPEN(UNIT=out_unit,FILE='nimset.out',STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     begin dialog to acquire the separate regions.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(4(/,a))')
     $  "Note: to create a topologically annular region from",
     $  "simply connected regions, first use stitch to create",
     $  "a single simply connected region that has all blocks. Then,",
     $  "run stitch separately to perform the last operation."
      WRITE(nim_wr,'(/a)')
     $  'Enter the number of seperate regions for the mesh.'
      READ(nim_rd,*) nreg
      ALLOCATE(infile(nreg))
      ALLOCATE(region(nreg))
      IF (nreg<=0) THEN
        WRITE(nim_wr,'(/a)') 'Need at least 1 region.'
      ENDIF
      IF (nreg<=0) STOP

      WRITE(nim_wr,'(/a)')
     $  'Enter the file name of the nimrod.in file for each region.'
      DO ireg=1,nreg
        WRITE(nim_wr,'(a,i4)',ADVANCE='YES')
     $    'Enter the nimrod.in file name for region ',ireg-1
        WRITE(nim_wr,'(a)',ADVANCE='NO')
     $    "> "
        READ(nim_rd,*) infile(ireg)
      ENDDO
c-----------------------------------------------------------------------
c     loop over each region, performing normal nimset operations
c     for each region, and then setting the region structure pointers to
c     retain the data.
c-----------------------------------------------------------------------
      nbl_total=0
      nrbl_total=0
      region_loop: DO ireg=1,nreg
c-----------------------------------------------------------------------
c       check for the namelist input file.
c-----------------------------------------------------------------------
        INQUIRE(FILE=infile(ireg),EXIST=file_stat)
        IF (.NOT.file_stat) THEN
          WRITE(nim_wr,*) 'The input file, ',TRIM(infile(ireg)),
     $                    ', does not exist.'
          STOP
        ENDIF
c-----------------------------------------------------------------------
c       reset essential defaults, then read namelist input and set
c       constants.
c-----------------------------------------------------------------------
        CALL reset_input
        CALL read_namelist(infile(ireg),.true.)
        IF (set_phys_constants) THEN
          CALL physdat_set(chrg_input,zeff_input,mi_input,me_input,
     $      gam_input,kblz_input,mu0_input,c_input)
        ELSE
          CALL physdat_set()
        ENDIF
c-----------------------------------------------------------------------
c       catch disabled and undefined options.
c-----------------------------------------------------------------------
        IF (.NOT.conform)
     $    CALL nim_stop('The conform=F option is not available.')
        IF (nbl_rim<2)
     $    CALL nim_stop('nbl_rim must be >= 2 for now.')
        IF (pieflag/='rblock'.AND.pieflag/='tblock0'.AND.
     $      pieflag/='tblock1') THEN
          CALL nim_stop('pieflag input '//TRIM(pieflag)//
     $                  ' is undefined.')
        ENDIF
        IF (rimflag/='rblock'.AND.rimflag/='tblock'.AND.
     $      rimflag/='none') THEN
          CALL nim_stop('rimflag input '//TRIM(rimflag)//
     $                  ' is undefined.')
        ENDIF
        IF ((pieflag=='tblock1'.OR.rimflag=='tblock').AND.
     $      gridshape/='flux') THEN
          CALL nim_stop
     $      ('gridshape must be flux when triangulation is used.')
        ENDIF
c-----------------------------------------------------------------------
c       check if the current region has y-direction periodicity.
c-----------------------------------------------------------------------
        IF (periodicity/='y-dir') allper=.false.
c-----------------------------------------------------------------------
c       set the number of Fourier components.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          IF (dealiase) THEN
            nmodes=2**lphi/3+1
          ELSE  !  only n up to nphi/2 - 1 are retained.
            nmodes=2**lphi/2
          ENDIF
        ELSE
          nmodes=lin_nmodes
        ENDIF
c-----------------------------------------------------------------------
c       if a reset file is specified, read this region from the reset
c       file instead of going through the nimset initialization.  this
c       differs from nimset's use of reset.
c-----------------------------------------------------------------------
        read_reg: IF (TRIM(reset_file)/='none') THEN
          file_save=dump_file
          dump_file=reset_file
          ALLOCATE(keff(nmodes))
          CALL dump_read(nmodes,keff,t,istep)
          dump_file=file_save
        ELSE read_reg
c-----------------------------------------------------------------------
c       initialize grid, rectangular cross-sections:
c-----------------------------------------------------------------------
        IF (gridshape/='rect'.AND.gridshape/='circ'.AND.
     $      gridshape/='rect_cir'.AND.gridshape/='flux'.AND.
     $      gridshape/='rectquad')
     $    CALL nim_stop("Nimset: gridshape "//TRIM(gridshape)//
     $                  " is not used by stitch.")
        IF (gridshape(1:4)=='rect') THEN
          SELECT CASE (gridshape)
          CASE ('rect')
            CALL rect_init
          CASE ('rect_cir')
            CALL rect_shaped_init
          CASE ('rectquad')
            CALL rect_quad_init
          END SELECT
          CALL block_init(nrbl)
          nbl =nrbl
          CALL seam0_init
          CALL seam_init
c-----------------------------------------------------------------------
c     grids based on polar meshes:
c-----------------------------------------------------------------------
        ELSE
           IF(gridshape.EQ.'circ')THEN
              CALL polar_circlegrid_init
           ELSEIF(gridshape.EQ.'flux')THEN
              CALL polar_fluxgrid_init
           ENDIF
           CALL polar_rblock_init(rb,nrbl)
c-----------------------------------------------------------------------
c        count the number of blocks excluding any rim blocks.
c-----------------------------------------------------------------------
           IF (TRIM(pieflag)=='rblock') THEN
             nbl=nrbl
           ELSE
             nbl=nrbl+1
           ENDIF
c-----------------------------------------------------------------------
c        add rim blocks, then allocate tb and seam.
c-----------------------------------------------------------------------
           SELECT CASE(TRIM(rimflag))
           CASE('rblock')
             nbl=nbl+nbl_rim
             CALL nim_stop('Unimplimented rimflag option.')
           CASE('tblock')
             nbl=nbl+nbl_rim
           END SELECT
           ALLOCATE(tb(nrbl+1:nbl))
           ALLOCATE(seam(nbl))
c-----------------------------------------------------------------------
c        perform core initialization.
c-----------------------------------------------------------------------
           SELECT CASE(TRIM(pieflag))
           CASE('rblock')
             CALL polar_seam_init(rb,tb,seam,nrbl,nbl)
           CASE('tblock0')
             CALL polar_tblock_init(tb,nrbl)
             CALL polar_seam_init(rb,tb,seam,nrbl,nbl)
           END SELECT
c-----------------------------------------------------------------------
c        rim initialization.
c-----------------------------------------------------------------------
           SELECT CASE(TRIM(rimflag))
           CASE('tblock')
             CALL nim_stop('Rimflag options have been disabled.')
           CASE('none')
             CALL seam0_init
           END SELECT
        ENDIF
c-----------------------------------------------------------------------
c       name the blocks and seams
c-----------------------------------------------------------------------
        DO ib=1,nrbl
           rb(ib)%id=ib
           rb(ib)%name='rblock'
        ENDDO
        DO ib=nrbl+1,nbl
           tb(ib)%id=ib
           tb(ib)%name='tblock'
        ENDDO
        seam0%id=0
        seam0%name='seam'
        DO ib=1,nbl
           seam(ib)%id=ib
           seam(ib)%name='seam'
        ENDDO
c-----------------------------------------------------------------------
c       initialize wavenumber array.
c-----------------------------------------------------------------------
        ALLOCATE(keff(nmodes))
        keff=(/(in,in=0,nmodes-1)/)
        IF (.NOT.nonlinear.AND.lin_nmax/=0) keff=keff+lin_nmax-nmodes+1
        IF (geom=='lin') keff=keff*twopi/per_length
        keff=keff*zperiod
c-----------------------------------------------------------------------
c       allocate dependent variables and nodal quantity initialization.
c-----------------------------------------------------------------------
        CALL var0_alloc(nmodes)
        IF (gridshape(1:4)=='rect') THEN
          CALL phys_init(nmodes,keff)
        ELSE IF (ASSOCIATED(eigvec)) THEN
          CALL eig_phys_init(nmodes,keff)
        ELSE
          CALL polar_phys_init(nmodes,keff)
        ENDIF
c-----------------------------------------------------------------------
c       transfer equilibrium fields to n=0 if requested.
c-----------------------------------------------------------------------
        IF (transfer_eq.AND.nonlinear)
     $    CALL eq_swap(nmodes,keff,geom,eq_flow,nedge,0._r8)
c-----------------------------------------------------------------------
c       add vacuum field-error perturbations.
c-----------------------------------------------------------------------
        CALL field_err_init(nmodes,keff)
c-----------------------------------------------------------------------
c       diagnose this region.  it has to come before a possible
c       reset--if nimplot contour/vector plots afterwards to verify.
c-----------------------------------------------------------------------
        CALL xy_slice(ireg,t)
c-----------------------------------------------------------------------
c       end of reset option.
c-----------------------------------------------------------------------
        ENDIF read_reg
        nbl_total=nbl_total+nbl
        nrbl_total=nrbl_total+nbl
c-----------------------------------------------------------------------
c       set the block and seam pointers for this region.
c-----------------------------------------------------------------------
        region(ireg)%nbl=nbl
        region(ireg)%nrbl=nrbl
        region(ireg)%rb=>rb
        region(ireg)%tb=>tb
        region(ireg)%seam=>seam
        region(ireg)%seam0=>seam0
        IF (ireg<nreg) DEALLOCATE(keff)
      ENDDO region_loop
c-----------------------------------------------------------------------
c     reset the ids within each region.
c-----------------------------------------------------------------------
      ibg=0
      ibgt=nrbl_total
      DO ireg=1,nreg
        rb=>region(ireg)%rb
        tb=>region(ireg)%tb
        seam=>region(ireg)%seam
        nbl=region(ireg)%nbl
        nrbl=region(ireg)%nrbl
        region(ireg)%irst=ibg+1
        region(ireg)%itst=ibgt+1
        DO ib=1,nrbl
          rb(ib)%id=rb(ib)%id+ibg
          seam(ib)%id=seam(ib)%id+ibg
          DO iv=1,seam(ib)%nvert
            np=SIZE(seam(ib)%vertex(iv)%ptr,2)
            DO ip=1,np
              IF (seam(ib)%vertex(iv)%ptr(1,ip)<=nrbl.AND.
     $            seam(ib)%vertex(iv)%ptr(1,ip)>0) THEN
                seam(ib)%vertex(iv)%ptr(1,ip)=
     $            seam(ib)%vertex(iv)%ptr(1,ip)+ibg
              ELSE IF (seam(ib)%vertex(iv)%ptr(1,ip)>nrbl) THEN
                seam(ib)%vertex(iv)%ptr(1,ip)=
     $            seam(ib)%vertex(iv)%ptr(1,ip)+ibgt
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DO ib=nrbl+1,nbl
          tb(ib)%id=tb(ib)%id+ibg
          seam(ib)%id=seam(ib)%id+ibgt
          DO iv=1,seam(ib)%nvert
            np=SIZE(seam(ib)%vertex(iv)%ptr,2)
            DO ip=1,np
              IF (seam(ib)%vertex(iv)%ptr(1,ip)<=nrbl.AND.
     $            seam(ib)%vertex(iv)%ptr(1,ip)>0) THEN
                seam(ib)%vertex(iv)%ptr(1,ip)=
     $            seam(ib)%vertex(iv)%ptr(1,ip)+ibg
              ELSE IF (seam(ib)%vertex(iv)%ptr(1,ip)>nrbl) THEN
                seam(ib)%vertex(iv)%ptr(1,ip)=
     $            seam(ib)%vertex(iv)%ptr(1,ip)+ibgt
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        DO iv=1,region(ireg)%seam0%nvert
          np=SIZE(region(ireg)%seam0%vertex(iv)%ptr,2)
          DO ip=1,np
            IF (region(ireg)%seam0%vertex(iv)%ptr(1,ip)<=nrbl.AND.
     $          region(ireg)%seam0%vertex(iv)%ptr(1,ip)>0) THEN
              region(ireg)%seam0%vertex(iv)%ptr(1,ip)=
     $          region(ireg)%seam0%vertex(iv)%ptr(1,ip)+ibg
            ELSE IF (region(ireg)%seam0%vertex(iv)%ptr(1,ip)>nrbl) THEN
              region(ireg)%seam0%vertex(iv)%ptr(1,ip)=
     $          region(ireg)%seam0%vertex(iv)%ptr(1,ip)+ibgt
            ENDIF
          ENDDO
        ENDDO
        ibg=ibg+nrbl
        ibgt=ibgt+nbl-nrbl
      ENDDO
c-----------------------------------------------------------------------
c     plot the seam0 structures to help the user define stitching
c     operations.
c-----------------------------------------------------------------------
      CALL draw_seam0(nreg,nrbl_total,nbl_total,region)
c-----------------------------------------------------------------------
c     perform the stitching operations on the seams.
c-----------------------------------------------------------------------
      IF (gridshape=='circ'.OR.gridshape=='flux'.OR.
     $    gridshape=='rect'.AND.allper) THEN
        CALL stitch_seams_polar(nreg,nrbl_total,nbl_total,nstitch,
     $                          region,seam0)
      ELSE
        CALL stitch_seams(nreg,nrbl_total,nbl_total,nstitch,
     $                    region,seam0)
      ENDIF
c-----------------------------------------------------------------------
c     write the regions into a dump file as a means of collecting
c     the rblock, tblock, and seam structures.  after that
c     read the dump file into global rbl, tbl, and seam structures.
c     all rblocks and rblock seams must be listed first.
c-----------------------------------------------------------------------
      nbl=nbl_total
      nrbl=nrbl_total
      CALL write_header(nmodes,keff,t,istep)
      CALL dump_write_seam(seam0)
      CALL seam_dealloc(seam0)
      DO ireg=1,nreg
        seam=>region(ireg)%seam
        DO ib=1,region(ireg)%nrbl
          CALL dump_write_seam(seam(ib))
        ENDDO
      ENDDO
      DO ireg=1,nreg
        seam=>region(ireg)%seam
        DO ib=region(ireg)%nrbl+1,region(ireg)%nbl
          CALL dump_write_seam(seam(ib))
        ENDDO
      ENDDO
      DO ireg=1,nreg
        rb=>region(ireg)%rb
        DO ib=1,region(ireg)%nrbl
          CALL dump_write_rblock(rb(ib))
        ENDDO
      ENDDO
      DO ireg=1,nreg
        tb=>region(ireg)%tb
        DO ib=region(ireg)%nrbl+1,region(ireg)%nbl
          CALL dump_write_tblock(tb(ib))
        ENDDO
      ENDDO

      CALL close_bin(dump_unit,TRIM(dump_file))
      NULLIFY(rb)
      NULLIFY(tb)
      NULLIFY(seam)
      CALL dump_read(nmodes,keff,t,istep)
c-----------------------------------------------------------------------
c     enter the new seam0 vertex numbers into the block seam structures.
c-----------------------------------------------------------------------
      DO iv=1,seam0%nvert
        DO ip=1,SIZE(seam0%vertex(iv)%ptr,2)
          jb=seam0%vertex(iv)%ptr(1,ip)
          jv=seam0%vertex(iv)%ptr(2,ip)
          DO jp=1,SIZE(seam(jb)%vertex(jv)%ptr,2)
            IF (seam(jb)%vertex(jv)%ptr(1,jp)==0)
     $          seam(jb)%vertex(jv)%ptr(2,jp)=iv
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     finally, run through all seams and eliminate any stranded
c     connections to seam0, which may occur from two or more stitches
c     across the same physical position.
c-----------------------------------------------------------------------
      IF (gridshape=='rect'.OR.gridshape=='rect-cir') THEN
        DO ib=1,nbl
          DO iv=1,seam(ib)%nvert
            np=SIZE(seam(ib)%vertex(iv)%ptr(1,:))
            ALLOCATE(count(1,np),ptmp(2,np))
            count(1,:)=1
            DO ip=1,np
              IF (seam(ib)%vertex(iv)%ptr(1,ip)==0) THEN
                count(1,ip)=0
                jb=seam(ib)%id
                ipcheck: DO iv0=1,seam0%nvert
                  DO jp=1,SIZE(seam0%vertex(iv0)%ptr(1,:))
                    IF (seam0%vertex(iv0)%ptr(1,jp)==jb.AND.
     $                  seam0%vertex(iv0)%ptr(2,jp)==iv) THEN
                      count(1,ip)=1
                      EXIT ipcheck
                    ENDIF
                  ENDDO
                ENDDO ipcheck
              ENDIF
            ENDDO
            IF (SUM(count(1,:))<np) THEN
              ptmp(:,:)=seam(ib)%vertex(iv)%ptr(:,:)
              DEALLOCATE(seam(ib)%vertex(iv)%ptr)
              np1=SUM(count(1,:))
              ALLOCATE(seam(ib)%vertex(iv)%ptr(2,np1))
              ip=1
              DO jp=1,np
                IF (count(1,jp)==1) THEN
                  seam(ib)%vertex(iv)%ptr(:,ip)=ptmp(:,jp)
                  ip=ip+1
                ENDIF
              ENDDO
            ENDIF
            DEALLOCATE(count,ptmp)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write new seam information if detflag=.true.
c-----------------------------------------------------------------------
      detflag=.true.
      IF (detflag) THEN
        WRITE(out_unit,'(/,/a/)')
     $    '*** seaming information after stitching ***'
        WRITE(out_unit,*) 'ibl ', 0, ' nvert ',seam0%nvert
        DO iv=1,seam0%nvert
         ip=SIZE(seam0%vertex(iv)%ptr,2)
         WRITE(out_unit,*) iv, ' block ptrs ',
     $                     seam0%vertex(iv)%ptr(1,1:ip)
         WRITE(out_unit,*) iv, ' vertx ptrs ',
     $                     seam0%vertex(iv)%ptr(2,1:ip)
        ENDDO
        DO ib=1,nbl
          WRITE(out_unit,*) 'ibl ', ib, ' nvert ',seam(ib)%nvert
          DO iv=1,seam(ib)%nvert
           ip=SIZE(seam(ib)%vertex(iv)%ptr,2)
           WRITE(out_unit,*) iv, ' block ptrs ',
     $                       seam(ib)%vertex(iv)%ptr(1,1:ip)
           WRITE(out_unit,*) iv, ' vertx ptrs ',
     $                       seam(ib)%vertex(iv)%ptr(2,1:ip)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write text output for the complete data structure.
c-----------------------------------------------------------------------
      CALL detail(istep,t)
c-----------------------------------------------------------------------
c     plot the grid and write the dump file.
c-----------------------------------------------------------------------
      CALL drawgrid
      CALL dump_write(nmodes,keff,t,istep)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')

c-----------------------------------------------------------------------
c     convenience subroutines for stitch only.
c-----------------------------------------------------------------------
      CONTAINS

c-----------------------------------------------------------------------
c       transfer seam0 from one structure to another.
c-----------------------------------------------------------------------
        SUBROUTINE seam0_tx(seamin,seamout)

        TYPE(edge_type), INTENT(IN) :: seamin
        TYPE(edge_type), INTENT(OUT) :: seamout

        INTEGER(i4) :: iv,nv,np

        nv=seamin%nvert
        seamout%nvert=nv
        ALLOCATE(seamout%vertex(nv))
        ALLOCATE(seamout%excorner(nv))
        seamout%excorner=seamin%excorner
        DO iv=1,nv
          np=SIZE(seamin%vertex(iv)%ptr,2)
          ALLOCATE(seamout%vertex(iv)%ptr(2,np))
          seamout%vertex(iv)%ptr(1:2,1:np)=
     $      seamin%vertex(iv)%ptr(1:2,1:np)
        ENDDO

        END SUBROUTINE seam0_tx

c-----------------------------------------------------------------------
c       deallocate seam arrays.
c-----------------------------------------------------------------------
        SUBROUTINE seam_dealloc(seam)

        TYPE(edge_type) :: seam

        INTEGER(i4) :: iv,nv,np

        nv=seam%nvert
        DO iv=1,nv
          DEALLOCATE(seam%vertex(iv)%ptr)
        ENDDO
        DEALLOCATE(seam%vertex)
        IF (seam%id==0) DEALLOCATE(seam%excorner)

        END SUBROUTINE seam_dealloc

      END PROGRAM stitch
c-----------------------------------------------------------------------
c     subprogram 2. write_header.
c     write the header information for a  dump file.
c-----------------------------------------------------------------------
      SUBROUTINE write_header(nmodes,keff,dump_time,dump_step)
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      USE dump
      IMPLICIT NONE
      
      INTEGER(i4), INTENT(IN) :: nmodes,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

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
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END  SUBROUTINE write_header
c-----------------------------------------------------------------------
c     subprogram 3. draw_seam0.
c     make xdraw plots of the seams around each region to help the
c     user define his stitching operations.
c-----------------------------------------------------------------------
      SUBROUTINE draw_seam0(nreg,nrblt,nblt,reg)
      USE local
      USE input
      USE region_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nreg,nrblt,nblt
      TYPE(region_type), DIMENSION(:), POINTER :: reg

      INTEGER(i4) :: ib,ix,iy,icell,iv,jv,kv,nv,ireg,ipass
      REAL(r8) :: smx,smy,small=1.e-4
      CHARACTER(16) :: fname
c-----------------------------------------------------------------------
c     open the first binary file.
c-----------------------------------------------------------------------
      CALL open_bin(grid_unit,"seam0s.bin","UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     draw the seam0 vertices for each region.
c-----------------------------------------------------------------------
      DO ireg=1,nreg
        nv=reg(ireg)%seam0%nvert
        DO iv=1,nv
          ib=reg(ireg)%seam0%vertex(iv)%ptr(1,1)
          jv=reg(ireg)%seam0%vertex(iv)%ptr(2,1)
          
          IF (ib<=nrblt) THEN
            ib=ib-reg(ireg)%irst+1
          ELSE
            ib=ib-reg(ireg)%itst+1
          ENDIF
          ix=reg(ireg)%seam(ib)%vertex(jv)%intxy(1)
          iy=reg(ireg)%seam(ib)%vertex(jv)%intxy(2)
          IF (ib<=reg(ireg)%nrbl) THEN
            WRITE(grid_unit) REAL(reg(ireg)%rb(ib)%rz%fs(:,ix,iy),4)
          ELSE
            WRITE(grid_unit) REAL(reg(ireg)%tb(ib)%tgeom%xs(ix),4),
     $                       REAL(reg(ireg)%tb(ib)%tgeom%ys(ix),4)
          ENDIF
        ENDDO
        WRITE(grid_unit)
      ENDDO
c-----------------------------------------------------------------------
c     close the file that plots all regions.
c-----------------------------------------------------------------------
      CALL close_bin(grid_unit,"seam0s.bin")
c-----------------------------------------------------------------------
c     now write separate graphics files for each region, so that the
c     user will be able to see vertex numbers.  the shift in region
c     index makex the xdraw labels in seam0s.bin plots correspond to
c     the file names for individual regions.
c-----------------------------------------------------------------------
      DO ireg=1,nreg
        IF (ireg<11) THEN
          WRITE(fname,'(a,i1,a)') "seam_reg0",ireg-1,".bin"
        ELSE
          WRITE(fname,'(a,i2,a)') "seam_reg",ireg-1,".bin"
        ENDIF
        CALL open_bin(grid_unit,TRIM(fname),"UNKNOWN","REWIND",32_i4)
        nv=reg(ireg)%seam0%nvert
        DO kv=0,nv
          IF (kv==0) THEN
            iv=nv
          ELSE
            iv=kv
          ENDIF
          ib=reg(ireg)%seam0%vertex(iv)%ptr(1,1)
          jv=reg(ireg)%seam0%vertex(iv)%ptr(2,1)
          
          IF (ib<=nrblt) THEN
            ib=ib-reg(ireg)%irst+1
          ELSE
            ib=ib-reg(ireg)%itst+1
          ENDIF
          ix=reg(ireg)%seam(ib)%vertex(jv)%intxy(1)
          iy=reg(ireg)%seam(ib)%vertex(jv)%intxy(2)
          DO ipass=1,4
            SELECT CASE(ipass)
            CASE(1,4)
              smx=0
              smy=0
            CASE(2)
              smx= small
              smy=-small
            CASE(3)
              smx=-small
              smy=-small
            END SELECT
            IF (ib<=reg(ireg)%nrbl) THEN
              WRITE(grid_unit)
     $          REAL(reg(ireg)%rb(ib)%rz%fs(1,ix,iy)+smx,4),
     $          REAL(reg(ireg)%rb(ib)%rz%fs(2,ix,iy)+smy,4)
            ELSE
              WRITE(grid_unit)
     $          REAL(reg(ireg)%tb(ib)%tgeom%xs(ix)+smx,4),
     $          REAL(reg(ireg)%tb(ib)%tgeom%ys(ix)+smy,4)
            ENDIF
          ENDDO
          WRITE(grid_unit)
        ENDDO
        CALL close_bin(grid_unit,TRIM(fname))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE draw_seam0
c-----------------------------------------------------------------------
c     subprogram 4. stitch_seams.
c     perform a set of stitching operations that connect two regions
c     at a time.
c-----------------------------------------------------------------------
      SUBROUTINE stitch_seams(nreg,nrblt,nblt,nst,reg,seam0)
      USE local
      USE input
      USE region_type_mod
      USE edge_type_mod
      USE stitch_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nreg,nrblt,nblt
      INTEGER(i4), INTENT(OUT) :: nst
      TYPE(region_type), DIMENSION(nreg), INTENT(INOUT) :: reg
      TYPE(edge_type) :: seam0

      INTEGER(i4) :: ib,iv,jb,jv,nv,ireg,ipass,ip,jp,np,ip1,np1
      INTEGER(i4) :: ircw,js0cw,ircc,js0cc,nvst,ivst,ipcw,npcw,ipcc,npcc
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ptmp,iunq
      TYPE(list_type), DIMENSION(:), POINTER :: s0starts
      TYPE(list_type), POINTER :: current,next,prev,cwsm0,ccsm0,
     $                 fstcc,fstcw,lstcc,lstcw,regpt
c-----------------------------------------------------------------------
c     reused output format
c-----------------------------------------------------------------------
  100 FORMAT(/,'  Connection of deactivated seam0 vertex at stitch',/,
     $       ' index ',i4,' is not allowed.')
  101 FORMAT(/,'  Starting vertex of region ',i4,/,
     $       ' has already been stitched out.')
c-----------------------------------------------------------------------
c     start by duplicating the seam0 information for each region into
c     linked lists that can be assembled relatively easily.
c-----------------------------------------------------------------------
      ALLOCATE(s0starts(nreg))
      CALL stitch_ll_init(reg,s0starts)
c-----------------------------------------------------------------------
c     instructions for the user.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(8(/,2x,a))')
     $  "At this point, you will begin to enter information for a set",
     $  "of stitching operations that connect two regions at a time.",
     $  "Use xdraw with drawsm0.in to plot the external seam0 for each",
     $  "region and hit l to see region labels numbered from 0.",
     $  "Use xdraw with drawsmreg.in to plot one external seam0 at a",
     $  "time (comment-out all but one binary in drawsmreg.in) and",
     $  "hit l to see seam vertex labels.  Here, 0 and the",
     $  "largest index are equivalent."
      WRITE(nim_wr,'(6(/,2x,a))')
     $  "For each operation, you will be prompted to enter the number",
     $  "of vertices in the stitch, the two region labels, and the",
     $  "starting seam0 vertex for each region.  Enter the region",
     $  "where the stitch proceeds with increasing vertex number (ccw",
     $  "around the region) first.  Note that a stitch can be defined",
     $  "to cross a seam origin."
c-----------------------------------------------------------------------
c     do loop over stitching operations.
c-----------------------------------------------------------------------
      nst=0
      DO
        nst=nst+1
        WRITE(nim_wr,'(/,2x,a,i3,a,/,2x,a)',ADVANCE='YES')
     $    "For stitch operation ",nst,", enter the number of seam",
     $    "vertices in the stitch, or enter -1 to finish."
        WRITE(nim_wr,'(2x,a)',ADVANCE='NO') "> "
        READ(nim_rd,*) nvst
        IF (nvst==-1) EXIT

        WRITE(nim_wr,'(/,2x,a,i3,a,3(/,2x,a))',ADVANCE='YES')
     $     "For stitch operation ",nst," enter the label of the",
     $    "region where the stitch runs counter-clockwise around the",
     $     "region and enter its seam0 vertex label where the stitch",
     $     "begins."
         WRITE(nim_wr,'(2x,a)',ADVANCE='NO') "> "
         READ(nim_rd,*) ircc,js0cc
         WRITE(nim_wr,'(/,2x,a,i3,a,2(/,2x,a))',ADVANCE='YES')
     $     "For stitch operation ",nst,", enter the region where the",
     $     "stitch runs clockwise around the region and enter its",
     $     "seam0 vertex label where the stitch begins."
         WRITE(nim_wr,'(2x,a)',ADVANCE='NO') "> "
        READ(nim_rd,*) ircw,js0cw

        ircc=ircc+1
        ircw=ircw+1
        IF (nvst>0) THEN
          IF (js0cc==0) js0cc=reg(ircc)%seam0%nvert
          IF (js0cw==0) js0cw=reg(ircw)%seam0%nvert
c-----------------------------------------------------------------------
c       first modify the seam0 linked lists to have the new topology.
c       start by finding the starting vertex in each list.
c-----------------------------------------------------------------------
          ccsm0=>s0starts(ircc)
          DO iv=1,js0cc-1
            ccsm0=>ccsm0%nextorig
          ENDDO
          IF (ircc/=ccsm0%ir) THEN
            WRITE(nim_wr,101) ircc
            nst=nst-1
            CYCLE
          ENDIF
          fstcc=>ccsm0
          cwsm0=>s0starts(ircw)
          DO iv=1,js0cw-1
            cwsm0=>cwsm0%nextorig
          ENDDO
          IF (ircw/=cwsm0%ir) THEN
            WRITE(nim_wr,101) ircw
            nst=nst-1
            CYCLE
          ENDIF
          fstcw=>cwsm0
c-----------------------------------------------------------------------
c       check if requested seam0 vertices are still active, and then
c       deactivate the vertices that are not at the ends of the stitch.
c-----------------------------------------------------------------------
          DO iv=1,nvst
            IF (.NOT.(ccsm0%active.AND.cwsm0%active)) THEN
              WRITE(nim_wr,100) iv
              nst=nst-1
              CYCLE
            ENDIF
            IF (iv/=1.AND.iv/=nvst) THEN
              ccsm0%active=.false.
              cwsm0%active=.false.
            ENDIF
            IF (iv/=nvst) THEN
              ccsm0=>ccsm0%nextorig
              cwsm0=>cwsm0%prevorig
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c       combine distinct pointer information in the last vertex.
c-----------------------------------------------------------------------
          lstcc=>ccsm0
          lstcw=>cwsm0
          CALL pointer_combine
          DEALLOCATE(ccsm0%ptr,cwsm0%ptr)
          ALLOCATE(ccsm0%ptr(2,np),cwsm0%ptr(2,np))
          ccsm0%ptr(:,:)=ptmp(:,:)
          cwsm0%ptr(:,:)=ptmp(:,:)
          ccsm0%excorner=.false.
          cwsm0%excorner=.false.
          ccsm0%prev=>cwsm0%prev
          cwsm0%next=>ccsm0%next
c-----------------------------------------------------------------------
c         if this was the end of a previous stitch, the previously
c         stitched region also needs to connect to the new next.
c-----------------------------------------------------------------------
          prev=>cwsm0%prev
          next=>prev%next
          next%next=>ccsm0%next
          DEALLOCATE(next%ptr)
          ALLOCATE(next%ptr(2,np))
          next%ptr(:,:)=ptmp(:,:)
          DEALLOCATE(ptmp)
          next=>ccsm0%next
          prev=>next%prev
          prev%prev=>cwsm0%prev

          IF (.NOT.(ccsm0%prev%active.OR.ccsm0%next%active))
     $      ccsm0%active=.false.
          IF (.NOT.(cwsm0%prev%active.OR.cwsm0%next%active))
     $      cwsm0%active=.false.
c-----------------------------------------------------------------------
c       combine distinct pointer information in the first vertex.
c-----------------------------------------------------------------------
          ccsm0=>fstcc
          cwsm0=>fstcw
          CALL pointer_combine
          DEALLOCATE(ccsm0%ptr,cwsm0%ptr)
          ALLOCATE(ccsm0%ptr(2,np),cwsm0%ptr(2,np))
          ccsm0%ptr(:,:)=ptmp(:,:)
          cwsm0%ptr(:,:)=ptmp(:,:)
          ccsm0%excorner=.false.
          cwsm0%excorner=.false.
          cwsm0%prev=>ccsm0%prev
          ccsm0%next=>cwsm0%next
c-----------------------------------------------------------------------
c         if this was the start of a previous stitch, the previously
c         stitched region also needs to connect to the new next.
c-----------------------------------------------------------------------
          prev=>ccsm0%prev
          next=>prev%next
          next%next=>cwsm0%next
          DEALLOCATE(next%ptr)
          ALLOCATE(next%ptr(2,np))
          next%ptr(:,:)=ptmp(:,:)
          DEALLOCATE(ptmp)
          next=>cwsm0%next
          prev=>next%prev
          prev%prev=>ccsm0%prev

          IF (.NOT.(ccsm0%prev%active.OR.ccsm0%next%active))
     $      ccsm0%active=.false.
          IF (.NOT.(cwsm0%prev%active.OR.cwsm0%next%active))
     $      cwsm0%active=.false.
c-----------------------------------------------------------------------
c       write the status of the seam0 linked list after each seaming.
c-----------------------------------------------------------------------
          IF (detflag) THEN
            CALL find_start
            prev=>current%prev
            next=>prev%next
            next%stflag=.true.
            WRITE(out_unit,*)
     $        "  **seam0 status after stitching step ",nst," **"
            WRITE(out_unit,*)
     $        "  current iv ",current%iv," linked iv ",next%iv
            DO 
              WRITE(out_unit,*) "  iv=",current%iv
              np=SIZE(current%ptr,2)
              WRITE(out_unit,*) "    jb=",(current%ptr(1,ip),ip=1,np)
              WRITE(out_unit,*) "    jv=",(current%ptr(2,ip),ip=1,np)
              current=>current%next
              IF (current%stflag) THEN
                current%stflag=.false.
                EXIT
              ENDIF
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c       use the pointer information at the beginning and end
c       of the linked vertices of the stitch to update the block-pointer
c       information for each block seam.
c-----------------------------------------------------------------------
          WRITE(nim_wr,'(/,2x,a,i3,a)')
     $      "Performing stitch operation ",nst,"."

          npcc=SIZE(fstcc%ptr,2)
          DO ipcc=1,npcc
            ib=fstcc%ptr(1,ipcc)
            iv=fstcc%ptr(2,ipcc)
            DO ireg=1,nreg
              IF (ib<reg(ireg)%irst+reg(ireg)%nbl) THEN
                jb=ib-reg(ireg)%irst+1
                EXIT
              ENDIF
            ENDDO
            DEALLOCATE(reg(ireg)%seam(jb)%vertex(iv)%ptr)
            ALLOCATE(reg(ireg)%seam(jb)%vertex(iv)%ptr(2,npcc))
            DO ip=1,npcc
              IF (ip==ipcc) THEN
                reg(ireg)%seam(jb)%vertex(iv)%ptr(:,ip)=(/0_i4,0_i4/)
              ELSE
                reg(ireg)%seam(jb)%vertex(iv)%ptr(:,ip)=fstcc%ptr(:,ip)
              ENDIF
            ENDDO
          ENDDO

          npcw=SIZE(lstcw%ptr,2)
          DO ipcw=1,npcw
            ib=lstcw%ptr(1,ipcw)
            iv=lstcw%ptr(2,ipcw)
            DO ireg=1,nreg
              IF (ib<reg(ireg)%irst+reg(ireg)%nbl) THEN
                jb=ib-reg(ireg)%irst+1
                EXIT
              ENDIF
            ENDDO
            DEALLOCATE(reg(ireg)%seam(jb)%vertex(iv)%ptr)
            ALLOCATE(reg(ireg)%seam(jb)%vertex(iv)%ptr(2,npcw))
            DO ip=1,npcw
              IF (ip==ipcw) THEN
                reg(ireg)%seam(jb)%vertex(iv)%ptr(:,ip)=(/0_i4,0_i4/)
              ELSE
                reg(ireg)%seam(jb)%vertex(iv)%ptr(:,ip)=lstcw%ptr(:,ip)
              ENDIF
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c       loop over the connected vertices that are not at the ends of the
c       stitch and modify the block seam pointers to establish all
c       connections.  a global seam0 will be defined later.
c
c       use ptmp for temporary space.
c-----------------------------------------------------------------------
          reg(ircc)%seam0%excorner(js0cc)=.false.
          reg(ircw)%seam0%excorner(js0cw)=.false.
          js0cc=js0cc+1
          js0cw=js0cw-1
          IF (js0cc>reg(ircc)%seam0%nvert) js0cc=1
          IF (js0cw==0) js0cw=reg(ircw)%seam0%nvert

          DO ivst=2,nvst-1
            CALL stitch_bl_adj(reg,ircc,js0cc,ircw,js0cw)
            js0cc=js0cc+1
            js0cw=js0cw-1
            IF (js0cc>reg(ircc)%seam0%nvert) js0cc=1
            IF (js0cw==0) js0cw=reg(ircw)%seam0%nvert
          ENDDO
          reg(ircc)%seam0%excorner(js0cc)=.false.
          reg(ircw)%seam0%excorner(js0cw)=.false.
        ENDIF
c-----------------------------------------------------------------------
c     For no connection cases, "connect" blocks by hand.
c-----------------------------------------------------------------------
        IF (nvst==0.AND.nreg>1) THEN
          next=>s0starts(ircc)%next
          prev=>s0starts(ircw)%next
          prev=>prev%next

          current=>next%next
          next%next=>prev
          prev%prev=>next

          prev=>s0starts(ircw)%next
          prev%next=>current
          current%prev=>prev
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     a one-region stitch is another special case.  it is used to create
c     a topologically annular region from a topologically circular one.
c     now that the block pointers have been changed for the stitch,
c     reset the linked list to connect the first cc vertex in the stitch
c     to the vertex after the last cc vertex in the stitch, and
c     connect the last cw vertex to the vertex after the first cw
c     vertex.
c-----------------------------------------------------------------------
      IF (nreg==1) THEN
        fstcc%next=>lstcc%next
        lstcc%next%prev=>fstcc
        lstcc%active=.false.
        lstcw%next=>fstcw%next
        fstcw%next%prev=>lstcw
        fstcw%active=.false.
      ENDIF
c-----------------------------------------------------------------------
c     now construct the composite seam0 structure based on the linked
c     list.  also eliminate any duplication.
c-----------------------------------------------------------------------
      CALL find_start
      prev=>current%prev
      next=>prev%next
      next%stflag=.true.
      nv=1
      DO 
        current=>current%next
        IF (current%stflag) THEN
          current%stflag=.false.
          EXIT
        ELSE
          nv=nv+1
          current%iv=nv
        ENDIF
      ENDDO
      seam0%id=0
      seam0%nvert=nv
      ALLOCATE(seam0%vertex(nv))
      ALLOCATE(seam0%excorner(nv))
      DO iv=1,nv
        np=SIZE(current%ptr,2)
        ALLOCATE(ptmp(1,np))
        ptmp(1,:)=1
        DO ip=2,np
          DO jp=1,ip-1
            IF (current%ptr(1,jp)==current%ptr(1,ip).AND.
     $          current%ptr(2,jp)==current%ptr(2,ip)) ptmp(1,ip)=0
          ENDDO
        ENDDO
        ALLOCATE(seam0%vertex(iv)%ptr(2,SUM(ptmp(1,:))))
        ip=1
        DO jp=1,np
          IF (ptmp(1,jp)==1) THEN
            seam0%vertex(iv)%ptr(:,ip)=current%ptr(:,jp)
            ip=ip+1
          ENDIF
        ENDDO
        DEALLOCATE(ptmp)
        seam0%excorner(iv)=current%excorner
        current=>current%next
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

c-----------------------------------------------------------------------
c     internal subroutine find_start finds an active starting point
c     in the linked list.  its result leaves the 'current' pointer
c     set at this location.
c
c     pointer_combine combines the distinct pointer information that
c     is in ccsm0 and cwsm0.  it leaves the ptmp array allocated.
c-----------------------------------------------------------------------
      CONTAINS

        SUBROUTINE find_start

        searchdo: DO ireg=1,nreg
          current=>s0starts(ireg)
          DO
            IF (current%active) EXIT searchdo
            current=>current%next
            IF (current%iv==1) EXIT
          ENDDO
        ENDDO searchdo

        RETURN
        END SUBROUTINE find_start

        SUBROUTINE pointer_combine

        np1=SIZE(ccsm0%ptr,2)
        np =SIZE(cwsm0%ptr,2)
        ALLOCATE(iunq(np,1))
        iunq=1
        DO ip1=1,np1
          DO ip=1,np
            IF (ccsm0%ptr(1,ip1)==cwsm0%ptr(1,ip).AND.
     $          ccsm0%ptr(2,ip1)==cwsm0%ptr(2,ip)) iunq(ip,1)=0
          ENDDO
        ENDDO
        np=np1+SUM(iunq)
        ALLOCATE(ptmp(2,np))
        ptmp(:,1:np1)=ccsm0%ptr(:,:)
        ip1=np1
        DO ip=1,SIZE(iunq,1)
          IF (iunq(ip,1)>0) THEN
            ip1=ip1+1
            ptmp(:,ip1)=cwsm0%ptr(:,ip)
          ENDIF
        ENDDO
        DEALLOCATE(iunq)

        RETURN
        END SUBROUTINE pointer_combine

      END SUBROUTINE stitch_seams
c-----------------------------------------------------------------------
c     subprogram 5. stitch_seams_polar.
c     perform a set of stitching operations that connect two regions
c     at a time in polar meshes.
c-----------------------------------------------------------------------
      SUBROUTINE stitch_seams_polar(nreg,nrblt,nblt,nst,reg,seam0)
      USE local
      USE input
      USE region_type_mod
      USE edge_type_mod
      USE stitch_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nreg,nrblt,nblt
      INTEGER(i4), INTENT(OUT) :: nst
      TYPE(region_type), DIMENSION(nreg), INTENT(INOUT) :: reg
      TYPE(edge_type) :: seam0

      INTEGER(i4) :: ib,iv,jb,jv,nv,ireg,ipass,ip,jp,np,np1
      INTEGER(i4) :: irout,js0out,irin,js0in,numch,readst,ivi,ivo
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: ptmp
      INTEGER(i4), DIMENSION(2) :: iread
      TYPE(list_type), DIMENSION(:), POINTER :: s0starts
      TYPE(list_type), POINTER :: current,next,prev,outsm0,insm0

      CHARACTER(72) :: cinput
      LOGICAL :: match
c-----------------------------------------------------------------------
c     reused output format
c-----------------------------------------------------------------------
  100 FORMAT(/,'  Connection of deactivated seam0 vertex at stitch',/,
     $       ' index ',i4,' is not allowed.')
  101 FORMAT(/,'  Starting vertex of region ',i4,/,
     $       ' has already been stitched out.')
c-----------------------------------------------------------------------
c     start by duplicating the seam0 information for each region into
c     linked lists that can be assembled relatively easily.
c-----------------------------------------------------------------------
      ALLOCATE(s0starts(nreg))
      CALL stitch_ll_init(reg,s0starts)
c-----------------------------------------------------------------------
c     instructions for the user.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(8(/,2x,a))')
     $  "At this point, you will begin to enter information for a set",
     $  "of stitching operations that connect two regions at a time.",
     $  "Use xdraw with drawsm0.in to plot the external seam0 for each",
     $  "region and hit l to see region labels numbered from 0.",
     $  "Use xdraw with drawsmreg.in to plot one external seam0 at a",
     $  "time (comment-out all but one binary in drawsmreg.in) and",
     $  "hit l to see seam vertex labels.  Here, 0 and the",
     $  "largest index are equivalent."
      WRITE(nim_wr,'(3(/,2x,a))')
     $  "For each operation, you will be prompted to enter the",
     $  "two region labels and one seam0 vertex from each to",
     $  "start the stitching operation."
c-----------------------------------------------------------------------
c     do loop over stitching operations.
c-----------------------------------------------------------------------
      nst=0
      perop: DO
        irin=-1
        nst=nst+1
        IF (nreg>1) THEN
          WRITE(nim_wr,'(/,2x,a,i3,a,2(/,2x,a))',ADVANCE='YES')
     $      "For stitch operation ",nst," enter the label of the",
     $      "inner region and any seam0 vertex within the stitch,",
     $      "or hit return if stitching is complete."
          WRITE(nim_wr,'(2x,a)',ADVANCE='NO')
     $       "> "
          READ(nim_rd,FMT='(a)',IOSTAT=readst,SIZE=numch,
     $         ADVANCE='NO') cinput
          IF (numch<=0) EXIT
          READ(cinput,*) iread
          irin=iread(1)
          js0in=iread(2)

          WRITE(nim_wr,'(/,2x,a,i3,a,2(/,2x,a))',ADVANCE='YES')
     $      "For stitch operation ",nst," enter the label of the",
     $      "outer region and its seam0 vertex that matches the,",
     $      "specified inner-region vertex."
          WRITE(nim_wr,'(2x,a)',ADVANCE='NO')
     $       "> "
          READ(nim_rd,*) irout,js0out

          irin=irin+1
          irout=irout+1
        ENDIF
c-----------------------------------------------------------------------
c       for these periodic regions, the seaming strategy is to refer to
c       block-seam connectivity to trace a stitched region.  connections
c       for the remaining sections of the inner and outer seam0s are
c       joined afterwards.
c-----------------------------------------------------------------------
        IF (irin>0) THEN
          IF (js0in==0) js0in=reg(irin)%seam0%nvert
          IF (js0out==0) js0out=reg(irout)%seam0%nvert
c-----------------------------------------------------------------------
c         find the specified vertex within the separate seam0 lists.
c-----------------------------------------------------------------------
          insm0=>s0starts(irin)
          DO iv=1,js0in-1
            insm0=>insm0%nextorig
          ENDDO
          IF (irin/=insm0%ir.OR..NOT.insm0%active) THEN
            WRITE(nim_wr,101) irin
            nst=nst-1
            CYCLE
          ENDIF
          outsm0=>s0starts(irout)
          DO iv=1,js0out-1
            outsm0=>outsm0%nextorig
          ENDDO
          IF (irout/=outsm0%ir.OR..NOT.outsm0%active) THEN
            WRITE(nim_wr,101) irout
            nst=nst-1
            CYCLE
          ENDIF
c-----------------------------------------------------------------------
c         loop through both regions, combining block-seam
c         information and deactivating the seam0 list members.
c-----------------------------------------------------------------------
          perdeact: DO
            ivi=insm0%iv
            ivo=outsm0%iv
            insm0%active=.false.
            outsm0%active=.false.
            CALL stitch_bl_adj(reg,irin,ivi,irout,ivo)

            nextin: DO
              insm0=>insm0%next
              match=.false.
              DO ip=1,SIZE(reg(irin)%seam0%vertex(ivi)%ptr,2)
                ib=reg(irin)%seam0%vertex(ivi)%ptr(1,ip)
                iv=reg(irin)%seam0%vertex(ivi)%ptr(2,ip)
                DO jp=1,SIZE(reg(irin)%seam0%vertex(insm0%iv)%ptr,2)
                  jb=reg(irin)%seam0%vertex(insm0%iv)%ptr(1,jp)
                  jv=reg(irin)%seam0%vertex(insm0%iv)%ptr(2,jp)-1
                  IF (jv<1) jv=reg(irin)%seam(jb-reg(irin)%irst+1)%nvert
                  IF (ib==jb.AND.iv==jv) match=.true.
                ENDDO
              ENDDO
              IF (match) EXIT nextin
            ENDDO nextin
            IF (.NOT.insm0%active) EXIT perdeact

            prevout: DO
              outsm0=>outsm0%prev
              match=.false.
              DO ip=1,SIZE(reg(irout)%seam0%vertex(ivo)%ptr,2)
                ib=reg(irout)%seam0%vertex(ivo)%ptr(1,ip)
                iv=reg(irout)%seam0%vertex(ivo)%ptr(2,ip)
                DO jp=1,SIZE(reg(irout)%seam0%vertex(outsm0%iv)%ptr,2)
                  jb=reg(irout)%seam0%vertex(outsm0%iv)%ptr(1,jp)
                  jv=reg(irout)%seam0%vertex(outsm0%iv)%ptr(2,jp)+1
                  IF (jv>reg(irout)%seam(jb-reg(irout)%irst+1)%nvert)
     $              jv=1
                  IF (ib==jb.AND.iv==jv) match=.true.
                ENDDO
              ENDDO
              IF (match) EXIT prevout
            ENDDO prevout
            IF (.NOT.outsm0%active) EXIT perdeact
          ENDDO perdeact
        ENDIF

      ENDDO perop
c-----------------------------------------------------------------------
c     create the new seam0 structure.  for periodic cases, entire
c     sections are left intact, so this is just a matter of collecting
c     the active vertices.
c-----------------------------------------------------------------------
      nv=0
      DO ireg=1,nreg
        current=>s0starts(ireg)
        DO
          IF (current%stflag) EXIT
          IF (current%active) nv=nv+1
          current%stflag=.true.
          current=>current%next
        ENDDO
      ENDDO

      seam0%id=0
      seam0%nvert=nv
      ALLOCATE(seam0%vertex(nv))
      ALLOCATE(seam0%excorner(nv))

      iv=0
      DO ireg=1,nreg
        current=>s0starts(ireg)
        DO
          IF (.NOT.current%stflag) EXIT
          IF (current%active) THEN
            iv=iv+1
            np=SIZE(current%ptr,2)
            ALLOCATE(ptmp(1,np))
            ptmp(1,:)=1
            DO ip=2,np
              DO jp=1,ip-1
                IF (current%ptr(1,jp)==current%ptr(1,ip).AND.
     $              current%ptr(2,jp)==current%ptr(2,ip)) ptmp(1,ip)=0
              ENDDO
            ENDDO
            ALLOCATE(seam0%vertex(iv)%ptr(2,SUM(ptmp(1,:))))
            ip=1
            DO jp=1,np
              IF (ptmp(1,jp)==1) THEN
                seam0%vertex(iv)%ptr(:,ip)=current%ptr(:,jp)
                ip=ip+1
              ENDIF
            ENDDO
            DEALLOCATE(ptmp)
            seam0%excorner(iv)=current%excorner
          ENDIF
          current%stflag=.false.
          current=>current%next
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE stitch_seams_polar
c-----------------------------------------------------------------------
c     subprogram 6. reset_input.
c     stitch reads a nimrod.in file for each region.  for the second
c     and subsequent reads, values not specificied in a nimrod.in file
c     are left at their values from a previous read, which can have
c     unintended consequences.  this operation resets the input
c     parameters that are most likely to cause problems to their
c     default values.
c-----------------------------------------------------------------------
      SUBROUTINE reset_input
      USE local
      USE input
      IMPLICIT NONE

      quadarc=.false.
      amp=(/0._r8,0._r8,0._r8,0._r8,0._r8,0._r8,0._r8,0._r8,0._r8,0._r8,
     $      0._r8,0._r8,0._r8,0._r8,0._r8/)
c-----------------------------------------------------------------------
      dvac=1._r8
      xvac=0._r8
      firstx=0._r8
      firsty=0._r8
      pres_offset=0._r8
      gamma_nimset=0._r8
      lam0=0._r8
      eq_flow='none'
      ncoil=0
      reset_file="none"

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE reset_input
