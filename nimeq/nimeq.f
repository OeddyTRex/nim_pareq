c-----------------------------------------------------------------------
c     file nimeq.f
c     contains a program that solves the Grad-Shafranov equation or
c     just inverts the del-star operator for NIMROD data.  the
c     program structure is based on nimplot.
c-----------------------------------------------------------------------
c     main program nimeq.
c-----------------------------------------------------------------------
      PROGRAM nimeq
      USE local
      USE fields
      USE input
      USE input_eq
      USE global
      USE time
      USE rblock
      USE tblock
      USE computation_pointers
      USE mpi_nim
      USE pardata
      USE contour_mod
      USE matrix_storage_mod
      USE plot_data
      USE regularity
      USE physdat
      USE nimeq_mod
      USE dump
      USE cell_type_mod
      USE dumpc
      USE spline
      USE nimeq_all
      IMPLICIT NONE

      INTEGER(i4) :: isel,ifile,nfiles,nfiles_old,nq,ndcon,isq,msq
      INTEGER :: read_stat,num_chars,ierror
      LOGICAL :: file_stat,plotit=.true.,first_list=.true.,
     $           first_gs=.true.,first_pass,exit_loop,fldat,
     $           gslast=.false.,first_cell=.true.
      REAL(r8), DIMENSION(:), POINTER :: setlabel
      REAL(r4), DIMENSION(5) :: sq
      CHARACTER(128), DIMENSION(1000) :: file_list
      CHARACTER(128) :: input_file='nimrod.in',drawgrad='drawgrad.in',
     $   drawpflux='drawpflux.in',xdrawex,dump_data,ctmp
      CHARACTER(1) :: repeat_sequence,mgt_flag
      CHARACTER(8), DIMENSION(2:4) ::
     $  choice=(/'xt_slice','yt_slice','contour '/)
      CHARACTER(16) :: rep_type
      TYPE(cell_type), POINTER :: item
c-----------------------------------------------------------------------
c     interface block.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE xdraw_plot(xdrawex,draw_file,dphi,flag,setlabel)
        USE local
        USE tecplot_mod
        USE vtk_mod
        IMPLICIT NONE

        CHARACTER(*), INTENT(IN) :: xdrawex,flag
        CHARACTER(*), INTENT(INOUT) :: draw_file
        REAL(r8), INTENT(IN) :: dphi
        REAL(r8), DIMENSION(:), POINTER :: setlabel
        END SUBROUTINE xdraw_plot
      END INTERFACE
c-----------------------------------------------------------------------
c     parallel initialization
c-----------------------------------------------------------------------
      CALL mpi_init(ierror)
      CALL mpi_comm_rank(mpi_comm_world,node,ierror)
      CALL mpi_comm_size(mpi_comm_world,nprocs,ierror)
c-----------------------------------------------------------------------
c     nullify block pointers to comply with the F90 standard.
c-----------------------------------------------------------------------
      NULLIFY(rb)
      NULLIFY(tb)
      NULLIFY(setlabel)
c-----------------------------------------------------------------------
c     dialogue descriptor.
c-----------------------------------------------------------------------
      open_dialogue: IF (node==0) THEN

      WRITE(nim_wr,'(7(/a))')
     $  'Welcome to NIMEQ, the Grad-Shafranov solver built within the',
     $  'NIMROD finite-element framework.  As with NIMPLOT, cursors',
     $  'serve as reminder of the program level:',
     $  '     = pre-looping level',
     $  '>    = computation option',
     $  '>>   = dump file selection',
     $  '>>>> = plot control'
c-----------------------------------------------------------------------
c     get xdraw executable path name.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(/5a,/a/a)')
     $  'Enter ',"'",'xdraw',"'",' (or the full path name for xdraw if',
     $  'not locally available) to plot after each selection, or',
     $  'hit return to just create the xdraw binary files.'
      WRITE(nim_wr,'(a)',ADVANCE='NO') '? '
      READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $           SIZE=num_chars) xdrawex
      IF (num_chars<1) THEN
        plotit=.false.
      ELSE
        INQUIRE(FILE=TRIM('draw.in'),EXIST=file_stat)
        IF (file_stat) WRITE(nim_wr,'(/a)')
     $    'CAUTION: draw.in will be overwitten when calling xdraw.'
      ENDIF
c-----------------------------------------------------------------------
c     read namelist input for file names and directories.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(input_file),EXIST=file_stat)
      DO
        IF (file_stat) THEN
          WRITE(nim_wr,'(/2a,/a,/a)')
     $      TRIM(input_file),' exists.',
     $      'Hit return if this input file is appropriate for the dump',
     $      'files that will be read, or enter another input file name.'
        ELSE
          WRITE(nim_wr,'(/2a,/a,/a)')
     $      TRIM(input_file),' does not exist.',
     $      'Enter the name of an input file appropriate for the dump',
     $      'files that will be read.'
        ENDIF
        WRITE(nim_wr,'(a)',ADVANCE='NO') '? '
        READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $             SIZE=num_chars) ctmp
        IF (num_chars>0) input_file=ADJUSTL(ctmp)
        INQUIRE(FILE=TRIM(input_file),EXIST=file_stat)
        IF (file_stat) EXIT
      ENDDO
      CALL read_namelist(TRIM(input_file),.false.)

      ENDIF open_dialogue

      nlayers=1  !  nimeq does not use layer decomposition.
      IF (nprocs > 1) CALL broadcast_input

      IF (set_phys_constants) THEN
        CALL physdat_set(chrg_input,zeff_input,mi_input,me_input,
     $    gam_input,kblz_input,mu0_input,c_input)
      ELSE
        CALL physdat_set()
      ENDIF
      ndcon=poly_degree
      inquire=.false.
      inquire_flux=.false.
c-----------------------------------------------------------------------
c     set nmodes and allocate keff.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        IF (dealiase) THEN
          nmodes=2**lphi/3+1
        ELSE  !  only n up to nphi/2 - 1 are retained.
          nmodes=2**lphi/2
        ENDIF
        nphi=2**lphi
      ELSE
        nmodes=lin_nmodes
        nphi=1
      ENDIF
      nmodes_total=nmodes
      ALLOCATE(keff(nmodes))
c----------------------------------------------------------------------
c     set nindex for psi diagnostics and gs solve.
c----------------------------------------------------------------------
      ALLOCATE(nindex(1))
      nindex=0
c-----------------------------------------------------------------------
c     interactive loop.
c-----------------------------------------------------------------------
      inter: DO
c-----------------------------------------------------------------------
c       output selection.
c-----------------------------------------------------------------------
        choice_dialogue: IF (node==0) THEN
        DO
          IF (.NOT.gslast) THEN
            WRITE(nim_wr,'(4(/a))',ADVANCE='YES')
     $         "> Plot Options:",
     $         ">   0: exit",
     $         ">   1: solve GS equation",
     $         ">   2: poloidal flux contours (for existing R*J_phi)"
          ELSE
            WRITE(nim_wr,'(5(/a))',ADVANCE='YES')
     $         "> Plot Options:",
     $         ">   0: exit",
     $         ">   1: solve GS equation",
     $         ">   2: poloidal flux contours (for existing R*J_phi)",
     $         ">   3: write dump file from previous GS solution"
          ENDIF
          WRITE(nim_wr,'(a)',ADVANCE='NO')
     $       ">? "
          READ(nim_rd,*) ctmp
          SELECT CASE(ctmp)
          CASE('0','1','2')
            READ(ctmp,*) isel
            EXIT
          CASE('3')
            IF (.NOT.gslast) THEN
              WRITE(nim_wr,*) '>  Option ', "'",TRIM(ctmp),"'",
     $          ' is used after solving the GS eqn.  Try again.'
            ELSE
              READ(ctmp,*) isel
              EXIT
            ENDIF
          CASE DEFAULT
            WRITE(nim_wr,*) '>  Option ', "'",TRIM(ctmp),"'",
     $        ' is not recognized.  Try again.'
          END SELECT
        ENDDO

        ENDIF choice_dialogue
        IF (nprocs>1)
     $    CALL mpi_bcast(isel,1,mpi_nim_int,0,mpi_comm_world,ierror)
c-----------------------------------------------------------------------
c       always read a file if not the same as previous.
c-----------------------------------------------------------------------
        read_file=.true.
c-----------------------------------------------------------------------
c       exit loop.
c-----------------------------------------------------------------------
        SELECT CASE(isel)
        CASE(0)
          EXIT inter
c-----------------------------------------------------------------------
c       solve the Grad-Shafranov equation.
c-----------------------------------------------------------------------
        CASE(1) 
          IF (node==0) CALL read_namelist_eq("nimeq.in",.false.)
          IF (nprocs>1) CALL broadcast_eq_input
          selection_type='computation'
          IF (node==0) CALL acquire_dump_name(dump_data,first_gs)
          CALL bcast_str(dump_data,128)
          CALL mpi_bcast(first_gs,1,mpi_nim_logical,0,
     $                   mpi_comm_world,ierror)
          flag='a'
c-----------------------------------------------------------------------
c         acquire the data for the next plot.
c-----------------------------------------------------------------------
          CALL acquire_data(dump_data,istep,t)
c-----------------------------------------------------------------------
c         read the data from fluxgrid's 1d.bin if desired.
c-----------------------------------------------------------------------
          IF ((f_model=="read1d".OR.pres_model=="read1d").AND.mfile<=0)
     $      CALL nim_stop("Nimeq: reading F and P profiles "//
     $                    "from 1d.bin requires mfile>0.") 
          IF (f_model=="read1d".OR.pres_model=="read1d") THEN
            INQUIRE(FILE=TRIM("1d.bin"),EXIST=file_stat)
            IF (.NOT.file_stat) CALL nim_stop
     $        ("Nimeq: the 1d.bin file from fluxgrid is missing.")
            IF (f0file/=0._r8.OR.p0file/=0._r8) THEN
              msq=mfile
              CALL spline_alloc(neqsq,msq/muse,3_i4)
              neqsq%xs(0)=0._r8
              neqsq%fs(0,1)=0._r8
              neqsq%fs(0,2)=f0file
              neqsq%fs(0,3)=p0file
            ELSE
              msq=mfile-1
              CALL spline_alloc(neqsq,msq/muse,3_i4)
            ENDIF
            CALL open_bin(xdr_unit,'1d.bin','OLD','REWIND',32_i4)
            DO isq=msq-mfile+1,msq
              READ(xdr_unit,IOSTAT=read_stat) sq
              IF (MODULO(isq,muse)==0) THEN
                neqsq%xs(isq/muse)=sq(1)
                neqsq%fs(isq/muse,1)=sq(1)
                neqsq%fs(isq/muse,2)=sq(2)
                neqsq%fs(isq/muse,3)=sq(3)
              ENDIF
            ENDDO
            CALL spline_fit(neqsq,"extrap")
            CALL close_bin(xdr_unit,'1d.bin')
          ENDIF
c-----------------------------------------------------------------------
c         build linked lists of cells if needed.
c-----------------------------------------------------------------------
c-TMP
c         IF (first_cell.AND.btop_check=="beq-trace") THEN
c           nblc=nbl  !  dump_reset routines use "c" data structures
c           nrblc=nrbl
c           poly_degreec=poly_degree
c           rbc=>rb
c           tbc=>tb
c           seamc=>seam
c           CALL make_cell
c           CALL make_search_map
c           item => start
c           DO WHILE(ASSOCIATED(item))
c             CALL map_cell_position(item)
c             item => item%next
c           ENDDO
c           WRITE(nim_wr,*) "Global Cell Link List Built"
c           first_cell=.false.
c         ENDIF
c-----------------------------------------------------------------------
c         compute and write delstar psi
c-----------------------------------------------------------------------
          IF (gs_type=='free') THEN
            CALL gsfree(.TRUE.,flag,phi,ndcon,.false.)
          ELSE
            CALL gssolve(.TRUE.,flag,phi,ndcon,.false.)
          ENDIF
          IF (plotit.AND.node==0)
     $      CALL xdraw_plot(xdrawex,drawgrad,dphi,flag,setlabel)
          last_dump='xxx_bogus_xxx'
          IF ((pres_model=='read1d'.OR.f_model=='read1d').AND.mfile>0)
     $      CALL spline_dealloc(neqsq)
          gslast=.true.
c-----------------------------------------------------------------------
c       contours of poloidal flux uses modified delstar operator.
c-----------------------------------------------------------------------
        CASE(2)
          IF (node==0) CALL read_namelist_eq("nimeq.in",.false.)
          IF (nprocs>1) CALL broadcast_eq_input
c-----------------------------------------------------------------------
c         calcuation of poloidal flux from current does not use 
c         matrix centering
c-----------------------------------------------------------------------
          selection_type='computation'
          IF (node==0)
     $      CALL acquire_dump_list(file_list,first_list,nfiles,
     $                             nfiles_old)
          CALL mpi_bcast(nfiles,1,mpi_nim_int,0,mpi_comm_world,ierror)
          CALL mpi_bcast(nfiles_old,1,mpi_nim_int,0,mpi_comm_world,
     $                   ierror)
          CALL mpi_bcast(first_list,1,mpi_nim_logical,0,
     $                   mpi_comm_world,ierror)
          DO ifile=1,nfiles
            CALL bcast_str(file_list(ifile),128)
          ENDDO
          flag='a'
          IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c         loop over file set, then plot
c-----------------------------------------------------------------------
          DO ifile=1,nfiles
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(file_list(ifile),istep,t)
            IF (nfiles>1) setlabel(ifile)=t
c-----------------------------------------------------------------------
c           compute and write delstar psi
c-----------------------------------------------------------------------
            CALL delstar((ifile==1),flag,phi,ndcon)
          ENDDO
          IF (plotit.AND.node==0)
     $      CALL xdraw_plot(xdrawex,drawpflux,dphi,flag,setlabel)
          IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          gslast=.false.
c-----------------------------------------------------------------------
c       write a dump file for the data that is save in the grid-block
c       data structures.
c-----------------------------------------------------------------------
        CASE(3)
          IF (node==0) CALL read_namelist_eq("nimeq.in",.false.)
          IF (nprocs>1) CALL broadcast_eq_input
          dump_name=dumpgs_name
          CALL dump_write(nmodes,nmodes_total,keff_total,t,istep,
     $                    dump_name)
          IF (node==0)
     $      WRITE(nim_wr,'(/,3a)')'Dump file ',TRIM(dump_name),
     $        ' has been written.'
c-----------------------------------------------------------------------
c       error message for user.
c-----------------------------------------------------------------------
        CASE DEFAULT
          IF (node==0) WRITE(nim_wr,'(/,a,i2,a)')
     $      'Option ',isel,' is not recognized.'
        END SELECT
      ENDDO inter
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')

      END PROGRAM nimeq
