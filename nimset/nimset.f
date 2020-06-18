c-----------------------------------------------------------------------
c     program nimset.
c     performs initialization and writes a dump file for nimrod.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     1.  main program, nimset
c-----------------------------------------------------------------------
      PROGRAM nimset
      USE local
      USE input
      USE fields
      USE dump
      USE dumpc
      USE seam_storage_mod
      USE nimset_init
      USE polar_init
      USE diagnose

      IMPLICIT NONE

      INTEGER(i4) :: istep=0,nmodes,ib,in
      LOGICAL :: file_stat
      REAL(r8):: t=0
      REAL(r8), DIMENSION(:), POINTER :: keff
c-----------------------------------------------------------------------
c     check for the namelist input file.
c-----------------------------------------------------------------------
      INQUIRE(FILE='nimrod.in',EXIST=file_stat)
      IF (.NOT.file_stat) THEN
        WRITE(nim_wr,*) 'The input file, nimrod.in, does not exist.'
        STOP
      ENDIF
c-----------------------------------------------------------------------
c     open files.
c-----------------------------------------------------------------------
      OPEN(UNIT=out_unit,FILE='nimset.out',STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     read namelist input and set constants.
c-----------------------------------------------------------------------
      CALL read_namelist('nimrod.in',.true.)
      IF (set_phys_constants) THEN
        CALL physdat_set(chrg_input,zeff_input,mi_input,me_input,
     $    gam_input,kblz_input,mu0_input,c_input)
      ELSE
        CALL physdat_set()
      ENDIF
c-----------------------------------------------------------------------
c     catch disabled and undefined options.
c-----------------------------------------------------------------------
      IF (.NOT.conform)
     $  CALL nim_stop('The conform=F option is not available.')
      IF (nbl_rim<2)
     $  CALL nim_stop('nbl_rim must be >= 2 for now.')
      IF (pieflag/='rblock'.AND.pieflag/='tblock0'.AND.
     $    pieflag/='tblock1') THEN
        CALL nim_stop('pieflag input '//TRIM(pieflag)//
     $                ' is undefined.')
      ENDIF
      IF (rimflag/='rblock'.AND.rimflag/='tblock'.AND.
     $    rimflag/='none') THEN
        CALL nim_stop('rimflag input '//TRIM(rimflag)//
     $                ' is undefined.')
      ENDIF
      IF ((pieflag=='tblock1'.OR.rimflag=='tblock').AND.
     $    gridshape/='flux') THEN
        CALL nim_stop
     $    ('gridshape must be flux when triangulation is used.')
      ENDIF
c-----------------------------------------------------------------------
c     if this is a reset that uses the mesh and equilibrium fields
c     from the reset file, start with a normal dump read, then
c     deallocate the solution fields before proceeding.
c-----------------------------------------------------------------------
      IF (reset_file/='none'.AND.reset_eq_mesh) THEN
        CALL dump_read(nmodes,keff,t,istep,readk=.true.,
     $                 rdfile=reset_file)
        CALL dealloc_soln
        DEALLOCATE(keff)
c-----------------------------------------------------------------------
c     initialize grid, rectangular cross-sections:
c-----------------------------------------------------------------------
      ELSE IF (gridshape/='rect'.AND.gridshape/='circ'.AND.
     $         gridshape/='rect_cir'.AND.gridshape/='flux'.AND.
     $         gridshape/='rectquad') THEN
        CALL nim_stop("Nimset: gridshape "//TRIM(gridshape)//
     $                " is not recognized.")
      ELSE IF (gridshape(1:4)=='rect') THEN
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
c     name the blocks and seams
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
c     initialize wavenumber array.
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
      ALLOCATE(keff(nmodes))
      keff=(/(in,in=0,nmodes-1)/)
      IF (.NOT.nonlinear.AND.lin_nmax/=0) keff=keff+lin_nmax-nmodes+1
      IF (geom=='lin') keff=keff*twopi/per_length
      keff=keff*zperiod
c-----------------------------------------------------------------------
c     allocate dependent variables and nodal quantity initialization.
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
c     read the perturbed components from another dump file (possibly
c     with a different number of Fourier components.
c-----------------------------------------------------------------------
      IF (TRIM(reset_file)/='none') THEN
        CALL dump_read_reset(nmodes)
        t=reset_time
        istep=reset_step
      ENDIF
c-----------------------------------------------------------------------
c     transfer equilibrium fields to n=0 if requested.
c-----------------------------------------------------------------------
      IF (transfer_eq.AND.nonlinear)
     $  CALL eq_swap(nmodes,keff,geom,eq_flow,nedge,0._r8)
c-----------------------------------------------------------------------
c     add vacuum field-error perturbations.
c-----------------------------------------------------------------------
      CALL field_err_init(nmodes,keff)
c-----------------------------------------------------------------------
c     write the dump file.
c-----------------------------------------------------------------------
      CALL dump_write(nmodes,keff,t,istep)
c-----------------------------------------------------------------------
c     diagnose.
c-----------------------------------------------------------------------
      CALL detail(istep,t)
      CALL drawgrid
      IF (reset_file=="none".OR..NOT.reset_eq_mesh)
     $  CALL xy_slice(istep,t)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')
      END PROGRAM nimset

c-----------------------------------------------------------------------
c     subprogram 2. dealloc_soln.
c     This subroutine deallocates the data structures for solution
c     solution fields but not those for the mesh or equilibrium.
c-----------------------------------------------------------------------
      SUBROUTINE dealloc_soln
      USE fields

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     loop over blocks and deallocate the arrays that are created by
c     var0_init.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL lagr_quad_dealloc(rb(ibl)%be)
        CALL lagr_quad_dealloc(rb(ibl)%ve)
        CALL lagr_quad_dealloc(rb(ibl)%pres)
        CALL lagr_quad_dealloc(rb(ibl)%prese)
        CALL lagr_quad_dealloc(rb(ibl)%tele)
        CALL lagr_quad_dealloc(rb(ibl)%tion)
        CALL lagr_quad_dealloc(rb(ibl)%nd)
        CALL lagr_quad_dealloc(rb(ibl)%conc)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tri_linear_dealloc(tb(ibl)%be)
        CALL tri_linear_dealloc(tb(ibl)%ve)
        CALL tri_linear_dealloc(tb(ibl)%pres)
        CALL tri_linear_dealloc(tb(ibl)%prese)
        CALL tri_linear_dealloc(tb(ibl)%tele)
        CALL tri_linear_dealloc(tb(ibl)%tion)
        CALL tri_linear_dealloc(tb(ibl)%nd)
        CALL tri_linear_dealloc(tb(ibl)%conc)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dealloc_soln

