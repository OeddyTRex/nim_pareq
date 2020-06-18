c-----------------------------------------------------------------------
c     file nimplot.f
c     contains a post-processing program and modified nimrod routines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  nimplot.
c-----------------------------------------------------------------------
c     main program:  nimplot
c     reads a nimrod dump file and generates graphics files from it.
c-----------------------------------------------------------------------
      PROGRAM nimplot
      USE local
      USE diagnose
      USE input
      USE contour_mod
      USE plot_data
      USE global
      USE pardata
      USE physdat
      IMPLICIT NONE

      INTEGER(i4) :: isel,ifile,nfiles,nfiles_old,nq,ndcon
      INTEGER :: read_stat,num_chars
      LOGICAL :: file_stat,plotit=.true.,first_list=.true.,
     $           first_xy=.true.,first_pass,exit_loop,fldat
      REAL(r8) :: tmp
      REAL(r8), DIMENSION(:), POINTER :: setlabel
      CHARACTER(128), DIMENSION(1000) :: file_list
      CHARACTER(128) :: input_file='nimrod.in',drawxy='drawxy.in',
     $   drawxt='drawxt.in',drawyt='drawyt.in',drawc='drawcon.in',
     $   xdrawex,dump_data,ctmp,drawfl='drawfl.in',
     $   drawpolfl='drawpolfl.in',drawen='drawen.in',
     $   drawdivb='drawdivb.in',drawvec='drawvec.in',
     $   drawjp='drawjpar.in',drawedj='drawedj.in',
     $   drawxypf='drawxypf.in',drawheat='drawheat.in',
     $   drawmach='drawmach.in'
      CHARACTER(1) :: repeat_sequence,mgt_flag
      CHARACTER(8), DIMENSION(2:4) ::
     $  choice=(/'xt_slice','yt_slice','contour '/)
      CHARACTER(16) :: rep_type
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
c     nullify block pointers to comply with the F90 standard.
c-----------------------------------------------------------------------
      NULLIFY(rb)
      NULLIFY(tb)
      NULLIFY(setlabel)
c-----------------------------------------------------------------------
c     nimplot is a serial program.
c-----------------------------------------------------------------------
      nprocs=1
      smallnum=SQRT(TINY(smallnum))
c-----------------------------------------------------------------------
c     dialogue descriptor.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(9(/a))')
     $  'Welcome to NIMPLOT, the interactive dump file processor for',
     $  'NIMROD.  There are five levels of input associated with',
     $  'making plots from the dumped data.  The cursors will serve',
     $  'as reminder of the level:',
     $  '     = pre-looping level',
     $  '>    = plot option',
     $  '>>   = dump file selection',
     $  '>>>  = data representation',
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
      IF (set_phys_constants) THEN
        CALL physdat_set(chrg_input,zeff_input,mi_input,me_input,
     $    gam_input,kblz_input,mu0_input,c_input)
      ELSE
        CALL physdat_set()
      ENDIF
      ndcon=poly_degree
c-----------------------------------------------------------------------
c     reset the linear system solver to an option that does not require
c     an external library.
c-----------------------------------------------------------------------
      solver='bl_diaga'
c-----------------------------------------------------------------------
c     set nmodes and allocate keff.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        nphi=2**lphi
        IF (dealiase) THEN
          nmodes=nphi/3+1
        ELSE  !  only n up to nphi/2 - 1 are retained.
          nmodes=nphi/2
        ENDIF
      ELSE
        nmodes=lin_nmodes
        nphi=1
      ENDIF
      ALLOCATE(keff(nmodes))
c-----------------------------------------------------------------------
c     interactive loop.
c-----------------------------------------------------------------------
      inter: DO
c-----------------------------------------------------------------------
c       output selection.
c-----------------------------------------------------------------------
        DO
          WRITE(nim_wr,'(19(/a))',ADVANCE='YES')
     $       "> Plot Options:",
     $       ">   0: exit",
     $       ">   1: xy_slices of solution fields",
     $       ">   2: xt_slices of solution fields",
     $       ">   3: yt_slices of solution fields",
     $       ">   4: contour plots of solution fields",
     $       ">   5: vector plots of solution fields",
     $       ">   6: toroidal ave. safety factor and parallel current",
     $       ">   7: toroidal ave. poloidal flux contours",
     $       ">   8: xy_slices vs. sqrt of poloidal flux",
     $       ">   9: energy spectra",
     $       ">  10: node-averaged div(b) contours",
     $       ">  11: mu0*J.B/B**2 contours",
     $       ">  12: <E>.<J> contours",
     $       ">  13: conductive heat flux vectors",
     $       ">  14: div(b) from elements (data output only)",
     $       ">  15: J from elements (data output only)",
     $       ">  16: dynamo from elements (data output only)",
     $       ">  17: ion acoustic speed, Mach number, par. Mach number"
          WRITE(nim_wr,'((a))',ADVANCE='YES')
     $       ">  90: Change data/element in contour plots"
          IF (do_jfromb=='y') THEN
            WRITE(nim_wr,'(a)',ADVANCE='YES')
     $         ">  91: Toggle jfromb call OFF for data plots"
          ELSE
            WRITE(nim_wr,'(a)',ADVANCE='YES')
     $         ">  91: Toggle jfromb call ON for data plots"
          ENDIF
          WRITE(nim_wr,'(a)',ADVANCE='NO')
     $       ">? "
          READ(nim_rd,*) ctmp
          SELECT CASE(ctmp)
          CASE('0','1','2','3','4','5','6','7','8','9','10','11','12',
     $         '13','14','15','16','17','90','91')
            READ(ctmp,*) isel
            EXIT
          CASE DEFAULT
            WRITE(nim_wr,*) '>  Option ', "'",TRIM(ctmp),"'",
     $        ' is not recognized.  Try again.'
          END SELECT
        ENDDO
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
c       xy_slices.
c-----------------------------------------------------------------------
        CASE(1,8)
          selection_type='data'
          inquire=.true.
          inquire_flux=.true.
c-----------------------------------------------------------------------
c         acquire the data for the next plot.
c-----------------------------------------------------------------------
          CALL acquire_dump_name(dump_data,first_xy)
          first_pass=.true.
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          xy_plots: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT xy_plots
            CALL acquire_data(dump_data,istep,t)
c-----------------------------------------------------------------------
c           compute poloidal flux if requested, then write slice files.
c-----------------------------------------------------------------------
            IF (isel==8) THEN
              fldat=.true.
            ELSE
              fldat=.false.
            ENDIF
            IF (plotit) THEN
              CALL xy_slice(0_i4,t,mode_plot_start,mode_plot_end,fldat)
              IF (isel==8) THEN
                CALL xdraw_plot(xdrawex,drawxypf,dphi,flag,setlabel)
              ELSE
                CALL xdraw_plot(xdrawex,drawxy,dphi,flag,setlabel)
              ENDIF
            ELSE
              CALL xy_slice(istep,t,mode_plot_start,mode_plot_end,fldat)
            ENDIF
            first_pass=.false.
          ENDDO xy_plots
c-----------------------------------------------------------------------
c       xt time slices:
c-----------------------------------------------------------------------
        CASE(2)
          selection_type='data'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         find the relative position.
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          xt_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT xt_rep
            WRITE(nim_wr,'(3(/a),f6.3,/a)',ADVANCE='NO')
     $        '>>> Enter the fractional y-position of the',
     $        '>>> x-direction slice, (0 <= y0fac <= 1).',
     $        '>>> Hit return to keep y0fac = ',y0fac,
     $        '>>>? '
            READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $           SIZE=num_chars) ctmp
            IF (num_chars>0) READ(ctmp,*) y0fac
            CALL time_slice_init
            inquire=.true.
            inquire_flux=.true.
c-----------------------------------------------------------------------
c           loop over file set.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              inquire=.false.
              inquire_flux=.false.
c-----------------------------------------------------------------------
c             diagnose.
c-----------------------------------------------------------------------
              CALL xt_slice(y0fac,istep,t,
     $                      mode_plot_start,mode_plot_end)
            ENDDO
            CALL time_slice_close
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawxt,dphi,flag,setlabel)
            first_pass=.false.
          ENDDO xt_rep
c-----------------------------------------------------------------------
c       yt time slices:
c-----------------------------------------------------------------------
        CASE(3)
          selection_type='data'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         find the relative position.
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          yt_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT yt_rep
            WRITE(nim_wr,'(3(/a),f6.3,/a)',ADVANCE='NO')
     $        '>>> Enter the fractional x-position of the',
     $        '>>> y-direction slice, (0 <= x0fac <= 1).',
     $        '>>> Hit return to keep x0fac = ',x0fac,
     $        '>>>? '
            READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $           SIZE=num_chars) ctmp
            IF (num_chars>0) READ(ctmp,*) x0fac
            CALL time_slice_init
            inquire=.true.
            inquire_flux=.true.
c-----------------------------------------------------------------------
c           loop over file set.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              inquire=.false.
              inquire_flux=.false.
c-----------------------------------------------------------------------
c             diagnose.
c-----------------------------------------------------------------------
              CALL yt_slice(x0fac,istep,t,
     $                      mode_plot_start,mode_plot_end)
            ENDDO
            CALL time_slice_close
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawyt,dphi,flag,setlabel)
            first_pass=.false.
          ENDDO yt_rep
c-----------------------------------------------------------------------
c       contour plots:
c-----------------------------------------------------------------------
        CASE(4)
          selection_type='data'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         loop for different data representations:
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          con_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,nfiles,geom)
            IF (exit_loop) EXIT con_rep
            inquire=.true.
            inquire_flux=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file or periodic coordinate for different data
c           sets.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
              IF (nfiles==1.AND.(flag=='p'.OR.flag=='t')) THEN
                iphi=0
                DO
                  CALL acquire_data(file_list(ifile),istep,t)
                  inquire=.false.
                  inquire_flux=.false.
                  nq=43
                  IF (iphi==0) CALL contour_init(nq,ndcon,"contour.bin")
                  CALL contour_step(1_i4,1_i4,ndcon)
                  IF (iphi*dphi>=1) EXIT
                  iphi=iphi+1
                ENDDO
              ELSE
                CALL acquire_data(file_list(ifile),istep,t)
                IF (nfiles>1) setlabel(ifile)=t
                inquire=.false.
                inquire_flux=.false.
                nq=13+30*(mode_plot_end-mode_plot_start+1)
                IF (ifile==1) CALL contour_init(nq,ndcon,"contour.bin")
                CALL contour_step(mode_plot_start,mode_plot_end,ndcon)
              ENDIF
            ENDDO
            CALL close_bin(con_unit,"contour.bin")
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawc,dphi,flag,setlabel)
            first_pass=.false.
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          ENDDO con_rep
c-----------------------------------------------------------------------
c       vector plots--use the same data file as contour plots, but
c       don't allow conversion to flux-surface components:
c-----------------------------------------------------------------------
        CASE(5)
          selection_type='data'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
          do_flux='n'
          inquire_flux=.false.
c-----------------------------------------------------------------------
c         loop for different data representations:
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          vec_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,nfiles,geom)
            IF (exit_loop) EXIT vec_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file or periodic coordinate for different data
c           sets.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
              IF (nfiles==1.AND.(flag=='p'.OR.flag=='t')) THEN
                iphi=0
                DO
                  CALL acquire_data(file_list(ifile),istep,t)
                  inquire=.false.
                  nq=43
                  IF (iphi==0) CALL contour_init(nq,ndcon,"contour.bin")
                  CALL contour_step(1_i4,1_i4,ndcon)
                  IF (iphi*dphi>=1) EXIT
                  iphi=iphi+1
                ENDDO
              ELSE
                CALL acquire_data(file_list(ifile),istep,t)
                IF (nfiles>1) setlabel(ifile)=t
                inquire=.false.
                nq=13+30*(mode_plot_end-mode_plot_start+1)
                IF (ifile==1) CALL contour_init(nq,ndcon,"contour.bin")
                CALL contour_step(mode_plot_start,mode_plot_end,ndcon)
              ENDIF
            ENDDO
            CALL close_bin(con_unit,"contour.bin")
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawvec,dphi,flag,setlabel)
            first_pass=.false.
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          ENDDO vec_rep
c-----------------------------------------------------------------------
c       safety factor and parallel current.
c-----------------------------------------------------------------------
        CASE(6)
          selection_type='computation'
          flag='a'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
c-----------------------------------------------------------------------
c         loop over file set.
c-----------------------------------------------------------------------
          DO ifile=1,nfiles
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(file_list(ifile),istep,t)
c-----------------------------------------------------------------------
c           compute q and lambda as functions of poloidal flux.
c-----------------------------------------------------------------------
            CALL fl_surface("q",(ifile==1),t)
          ENDDO
          IF (plotit) CALL xdraw_plot(xdrawex,drawfl,dphi,flag,setlabel)
c-----------------------------------------------------------------------
c       flux surface plots.
c-PRE   at some point, this should be generalized to include tblocks.
c-----------------------------------------------------------------------
        CASE(7)
          selection_type='computation'
          flag='a'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c         loop over file set.
c-----------------------------------------------------------------------
          DO ifile=1,nfiles
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(file_list(ifile),istep,t)
            IF (nfiles>1) setlabel(ifile)=t
c-----------------------------------------------------------------------
c           find the poloidal flux function and write it to the end
c           of the contour plot file.
c-----------------------------------------------------------------------
            CALL fl_surface("polflux",(ifile==1),t)
          ENDDO
          IF (plotit)
     $      CALL xdraw_plot(xdrawex,drawpolfl,dphi,flag,setlabel)
          IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
c-----------------------------------------------------------------------
c       energy spectra.
c-----------------------------------------------------------------------
        CASE(9)
          selection_type='computation'
          flag='a'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
c-----------------------------------------------------------------------
c         loop over file set, then plot the spectra.
c-----------------------------------------------------------------------
          DO ifile=1,nfiles
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(file_list(ifile),istep,t)
c-----------------------------------------------------------------------
c           compute and write energy spectra.
c-----------------------------------------------------------------------
            CALL energies((ifile==1))
          ENDDO
          IF (plotit) CALL xdraw_plot(xdrawex,drawen,dphi,flag,setlabel)
c-----------------------------------------------------------------------
c       node-averaged div(b) contours.
c-----------------------------------------------------------------------
        CASE(10)
          selection_type='computation'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop.
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          divb_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,nfiles,geom)
            IF (exit_loop) EXIT divb_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file set, then plot the contours.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              IF (nfiles>1) setlabel(ifile)=t
              inquire=.false.
c-----------------------------------------------------------------------
c             compute and write div(b).
c-----------------------------------------------------------------------
              IF (nfiles==1.AND.(flag=='p'.OR.flag=='t')) phi=dphi
              CALL divb((ifile==1),flag,mode_plot_start,
     $                  mode_plot_end,phi,ndcon)
            ENDDO
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawdivb,dphi,flag,setlabel)
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
            first_pass=.false.
          ENDDO divb_rep
c-----------------------------------------------------------------------
c       mu0*J.B/B**2 contours.
c-----------------------------------------------------------------------
        CASE(11)
          selection_type='computation'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop.
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          jpar_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,nfiles,geom)
            IF (exit_loop) EXIT jpar_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file set, then plot the contours.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              IF (nfiles>1) setlabel(ifile)=t
              inquire=.false.
c-----------------------------------------------------------------------
c             compute and write parallel j.
c-----------------------------------------------------------------------
              IF (nfiles==1.AND.(flag=='p'.OR.flag=='t')) phi=dphi
              CALL parallel_current((ifile==1),flag,mode_plot_start,
     $                              mode_plot_end,phi,ndcon)
            ENDDO
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawjp,dphi,flag,setlabel)
            first_pass=.false.
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          ENDDO jpar_rep
c-----------------------------------------------------------------------
c       <E>.<J> contours, where mode selection determines what part
c       of <E> is displayed.
c-----------------------------------------------------------------------
        CASE(12)
          selection_type='computation'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop.
c-----------------------------------------------------------------------
          edj_rep: DO
            WRITE(nim_wr,'(2(/a))')
     $       ">>> The Data Options for this plot select controls which",
     $       ">>> Fourier components appear in the <E> computation."
            CALL data_rep(first_pass,exit_loop,'correlation',0_i4,geom)
            IF (exit_loop) EXIT edj_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file set, then plot the contours.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              IF (nfiles>1) setlabel(ifile)=t
              inquire=.false.
c-----------------------------------------------------------------------
c             compute and write the power densities.
c-----------------------------------------------------------------------
              CALL e_dot_j((ifile==1),flag,mode_plot_start,ndcon)
            ENDDO
            IF (plotit) 
     $        CALL xdraw_plot(xdrawex,drawedj,dphi,flag,setlabel)
            first_pass=.false.
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          ENDDO edj_rep
c-----------------------------------------------------------------------
c       heat flux vectors, where mode selection determines what part
c       of q is displayed.
c-----------------------------------------------------------------------
        CASE(13)
          selection_type='computation'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop.
c-----------------------------------------------------------------------
          heat_rep: DO
            WRITE(nim_wr,'(2(/a))')
     $       ">>> The Data Options for this plot select controls which",
     $       ">>> Fourier components of T appear in the q computation."
            CALL data_rep(first_pass,exit_loop,'correlation',0_i4,geom)
            IF (exit_loop) EXIT heat_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file set, then plot the contours.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              IF (nfiles>1) setlabel(ifile)=t
              inquire=.false.
c-----------------------------------------------------------------------
c             compute and write heat flux vectors.
c-----------------------------------------------------------------------
              CALL heat_flux((ifile==1),flag,mode_plot_start,ndcon)
            ENDDO
            IF (plotit) 
     $        CALL xdraw_plot(xdrawex,drawheat,dphi,flag,setlabel)
            first_pass=.false.
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
          ENDDO heat_rep
c-----------------------------------------------------------------------
c       discontinuous div(b) contours from elements.
c-----------------------------------------------------------------------
        CASE(14)
          selection_type='computation'
          CALL acquire_dump_name(dump_data,first_xy)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop, but only use the first file.
c-----------------------------------------------------------------------
          rep_type='linear'
          divb_el_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT divb_el_rep
            inquire=.true.
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(dump_data,istep,t)
            inquire=.false.
c-----------------------------------------------------------------------
c           compute and write div(b) directly to a data file.
c-----------------------------------------------------------------------
            CALL divb_element(flag,mode_plot_start,mode_plot_end)
            first_pass=.false.
          ENDDO divb_el_rep
c-----------------------------------------------------------------------
c       discontinuous J contours from elements.
c-----------------------------------------------------------------------
        CASE(15)
          selection_type='computation'
          CALL acquire_dump_name(dump_data,first_xy)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop, but only use the first file.
c-----------------------------------------------------------------------
          rep_type='linear'
          j_el_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT j_el_rep
            inquire=.true.
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(dump_data,istep,t)
            inquire=.false.
c-----------------------------------------------------------------------
c           compute and write J directly to a data file.
c-----------------------------------------------------------------------
            CALL j_element(flag,mode_plot_start,mode_plot_end)
            first_pass=.false.
          ENDDO j_el_rep
c-----------------------------------------------------------------------
c       discontinuous dynamo contours from elements.
c-----------------------------------------------------------------------
        CASE(16)
          selection_type='computation'
          CALL acquire_dump_name(dump_data,first_xy)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop, but only use the first file.
c-----------------------------------------------------------------------
          rep_type='linear'
          dyn_el_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,0_i4,geom)
            IF (exit_loop) EXIT dyn_el_rep
            inquire=.true.
c-----------------------------------------------------------------------
c           acquire the data for the next plot.
c-----------------------------------------------------------------------
            CALL acquire_data(dump_data,istep,t)
            inquire=.false.
c-----------------------------------------------------------------------
c           compute and write J directly to a data file.
c-----------------------------------------------------------------------
            CALL dyn_element(flag,mode_plot_start,mode_plot_end)
            first_pass=.false.
          ENDDO dyn_el_rep
c-----------------------------------------------------------------------
c       Mach number calculated along with ion acoustic speed
c-----------------------------------------------------------------------
        CASE(17)
          selection_type='computation'
          CALL acquire_dump_list(file_list,first_list,nfiles,nfiles_old)
          first_pass=.true.
c-----------------------------------------------------------------------
c         start a data-representation loop.
c-----------------------------------------------------------------------
          rep_type='linear'
          IF (nonlinear) rep_type='nonlinear'
          mach_rep: DO
            CALL data_rep(first_pass,exit_loop,rep_type,nfiles,geom)
            IF (exit_loop) EXIT mach_rep
            inquire=.true.
            IF (nfiles>1) ALLOCATE(setlabel(nfiles))
c-----------------------------------------------------------------------
c           loop over file set, then plot the contours.
c-----------------------------------------------------------------------
            DO ifile=1,nfiles
c-----------------------------------------------------------------------
c             acquire the data for the next plot.
c-----------------------------------------------------------------------
              CALL acquire_data(file_list(ifile),istep,t)
              IF (nfiles>1) setlabel(ifile)=t
              inquire=.false.
c-----------------------------------------------------------------------
c             compute and write mach number.
c-----------------------------------------------------------------------
              IF (nfiles==1.AND.(flag=='p'.OR.flag=='t')) phi=dphi
              CALL mach((ifile==1),flag,mode_plot_start,
     $                  mode_plot_end,phi,ndcon)
            ENDDO
            IF (plotit)
     $        CALL xdraw_plot(xdrawex,drawmach,dphi,flag,setlabel)
            IF (ASSOCIATED(setlabel)) DEALLOCATE(setlabel)
            first_pass=.false.
          ENDDO mach_rep
c-----------------------------------------------------------------------
c       change the amount of data plotted per element in contour and
c       vector plots.
c-----------------------------------------------------------------------
        CASE(90)
          WRITE(nim_wr,'(3(/a),2(i2,a),/a)',ADVANCE='NO')
     $      '>> Please enter the number of data points per element ',
     $      '>> for contour and vector plots.  Poly_degree (default',
     $      '>> value) is ',poly_degree,' and the present value is ',
     $          ndcon,'.',
     $      '>>? '
          READ(nim_rd,*) ndcon
c-----------------------------------------------------------------------
c       turn the nodal J computation on or off for field-data plots such
c       as slice and vector plots.
c-----------------------------------------------------------------------
        CASE(91)
          IF (do_jfromb=='y') THEN
            do_jfromb='n'
          ELSE
            do_jfromb='y'
          ENDIF
c-----------------------------------------------------------------------
c       catch invalid plot options.
c-----------------------------------------------------------------------
        CASE DEFAULT
          WRITE(nim_wr,'(/,a,i2,a)')'Option ',isel,' is not recognized.'
        END SELECT
      ENDDO inter
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')
      END PROGRAM nimplot
