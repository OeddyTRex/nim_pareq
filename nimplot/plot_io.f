c-----------------------------------------------------------------------
c     file plot_io.f
c     contains routines that query the user for plotting information.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  xdraw_plot.
c     2.  data_rep.
c     3.  acquire_dump_name.
c     4.  acquire_dump_list.
c-----------------------------------------------------------------------
c     subprogram 1. xdraw_plot.
c     call the xdraw executable for the current selection.
c-----------------------------------------------------------------------
      SUBROUTINE xdraw_plot(xdrawex,draw_file,dphi,flag,setlabel)
      USE local
      USE tecplot_mod
      USE vtk_mod
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: xdrawex,flag
      CHARACTER(*), INTENT(INOUT) :: draw_file
      REAL(r8), INTENT(IN) :: dphi
      REAL(r8), DIMENSION(:), POINTER :: setlabel

      INTEGER :: read_stat,num_chars
      INTEGER(i4) :: iset,nsets
      LOGICAL :: file_stat
      CHARACTER(128) :: new_draw,path,plt_file,binfile
      CHARACTER(64) :: setname
      CHARACTER(1) :: pl_type
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
c-----------------------------------------------------------------------
c     suggested draw*.in file check.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(draw_file),EXIST=file_stat)
      IF (.NOT.file_stat) THEN
        new_draw='../draw/'//TRIM(draw_file)
        INQUIRE(FILE=TRIM(new_draw),EXIST=file_stat)
        IF (file_stat) draw_file=new_draw
      ENDIF
c-----------------------------------------------------------------------
c     query and test for the draw*.in file.
c-----------------------------------------------------------------------
      DO
        IF (file_stat) THEN
          WRITE(nim_wr,'(/a,/3a,/a)',ADVANCE='NO')
     $      '>>>> Enter the name of the draw*.in file, or hit return',
     $      '>>>> to use ',TRIM(draw_file),'.',
     $      '>>>>? '
        ELSE
          WRITE(nim_wr,'(/a,/2a,/a)',ADVANCE='NO')
     $      '>>>> Enter the name of the draw*.in file.',
     $      '>>>>? '
        ENDIF
        READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $       SIZE=num_chars) new_draw
        IF (num_chars>0) draw_file=TRIM(ADJUSTL(new_draw))
        INQUIRE(FILE=TRIM(draw_file),EXIST=file_stat)
        IF (file_stat) EXIT
        WRITE(nim_wr,'(/3a)') '>>>> File ',TRIM(draw_file),
     $    ' does not exist.'
      ENDDO
c-----------------------------------------------------------------------
c     drawing loop.
c-----------------------------------------------------------------------
      DO
        IF (file_stat) THEN
          CALL system('cp '//TRIM(draw_file)//' draw.in')
          CALL system(TRIM(xdrawex))
c-----------------------------------------------------------------------
c         convert xdraw binary to ascii tecplot data.
c-----------------------------------------------------------------------
          WRITE(nim_wr,'(5(/a))',ADVANCE='NO')
     $      '>>>> Enter a file name for an ascii version of the xdraw',
     $      '>>>> binary just displayed, or hit return to avoid this',
     $      '>>>> data conversion.  If the file name ends in vtk,',
     $      '>>>> then the file is in VTK format; otherwise, tecplot.',
     $      '>>>>? '
          READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $         SIZE=num_chars) plt_file
          IF (num_chars>0) THEN
            CALL read_draw('draw.in',pl_type,varnames,setname,binfile)
            IF (setname==' ') setname='x set'
            IF (pl_type=='g') THEN
              IF (plt_file(num_chars-2:num_chars)=='vtk') THEN
                CALL vtk_slice_conv(TRIM(binfile),TRIM(plt_file),.true.,
     $                              varnames)
              ELSE
                CALL slice_conv(TRIM(binfile),TRIM(plt_file),.true.,
     $                          varnames)
              ENDIF
            ELSE
              IF (ASSOCIATED(setlabel)) THEN
                nsets=SIZE(setlabel)
              ELSE IF (flag=='p'.OR.flag=='t') THEN
                IF (flag=='p') setname='periodic coordinate'
                nsets=0
                DO 
                  nsets=nsets+1
                  IF ((nsets-1)*dphi>=1) EXIT
                ENDDO
                ALLOCATE(setlabel(nsets))
                DO iset=1,nsets
                  setlabel(iset)=MIN((iset-1)*dphi,1._r8)
                ENDDO
              ELSE
                nsets=1
              ENDIF
              IF (flag=='t') THEN
                IF (plt_file(num_chars-2:num_chars)=='vtk') THEN
                  CALL vtk_contour_conv(TRIM(binfile),TRIM(plt_file),
     $              .true.,varnames,nsets,setname,setlabel,'y')
                ELSE
                  CALL contour_conv(TRIM(binfile),TRIM(plt_file),.true.,
     $                              varnames,nsets,setname,setlabel,'y')
                ENDIF
              ELSE
                IF (plt_file(num_chars-2:num_chars)=='vtk') THEN
                  CALL vtk_contour_conv(TRIM(binfile),TRIM(plt_file),
     $              .true.,varnames,nsets,setname,setlabel,'n')
                ELSE
                  CALL contour_conv(TRIM(binfile),TRIM(plt_file),.true.,
     $                              varnames,nsets,setname,setlabel,'n')
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          CALL system('rm draw.in')
        ENDIF
c-----------------------------------------------------------------------
c       back to dialogue.
c-----------------------------------------------------------------------
        WRITE(nim_wr,'(4(/a))',ADVANCE='NO')
     $    '>>>> Hit return to return to the data selection level, or',
     $    '>>>> enter the name of a new draw*.in file for viewing',
     $    '>>>> different quantities in the same binary file.',
     $    '>>>>? '
        READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $       SIZE=num_chars) new_draw
        IF (num_chars<1) EXIT
        draw_file=ADJUSTL(new_draw)
        INQUIRE(FILE=TRIM(draw_file),EXIST=file_stat)
        IF (.NOT.file_stat) WRITE(nim_wr,'(/3a)')
     $    '>>>> File ',TRIM(draw_file),' does not exist.'
      ENDDO

      RETURN
      END SUBROUTINE xdraw_plot
c-----------------------------------------------------------------------
c     subprogram 2. data_rep.
c     make a request for the data representation.
c-----------------------------------------------------------------------
      SUBROUTINE data_rep(first_pass,exit_loop,pl_type,nfiles,geom)
      USE local
      USE plot_data
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: pl_type,geom
      LOGICAL, INTENT(IN) :: first_pass
      LOGICAL, INTENT(OUT) :: exit_loop
      INTEGER(i4), INTENT(IN) :: nfiles

      INTEGER :: read_stat,num_chars
c-----------------------------------------------------------------------
c     set the data representation flag or get an exit signal.
c-----------------------------------------------------------------------
      exit_loop=.false.
      datar: DO 
        IF (first_pass) THEN
          SELECT CASE(pl_type)
          CASE ('nonlinear')
            IF (nfiles==1.AND.geom=='tor') THEN
              WRITE(nim_wr,'(12(/a))',ADVANCE='NO')
     $             ">>> Data Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>      (Requires writing special draw*.in files ",
     $             ">>>       if dump file has >1 Fourier comp.)",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>   p: plot sets of evenly spaced data in",
     $             ">>>      the periodic coordinate.",
     $             ">>>   t: same as p, but transform vector",
     $             ">>>      components to cartesian for toroidal",
     $             ">>>      3D plots (not xdraw).",
     $             ">>>?  "
              READ(nim_rd,*) flag
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c'.OR.
     $            flag=='p'.OR.flag=='t') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ELSE IF (nfiles==1) THEN
              WRITE(nim_wr,'(9(/a))',ADVANCE='NO')
     $             ">>> Data Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>      (Requires writing special draw*.in files ",
     $             ">>>       if dump file has >1 Fourier comp.)",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>   p: plot sets of evenly spaced data in",
     $             ">>>      the periodic coordinate.",
     $             ">>>?  "
              READ(nim_rd,*) flag
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c'.OR.
     $            flag=='p') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ELSE
              WRITE(nim_wr,'(7(/a))',ADVANCE='NO')
     $             ">>> Data Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>      (Requires writing special draw*.in files ",
     $             ">>>       if dump file has >1 Fourier comp.)",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>?  "
              READ(nim_rd,*) flag
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ENDIF
          CASE('linear')
            WRITE(nim_wr,'(6(/a))',ADVANCE='NO')
     $           ">>> Data Options:",
     $           ">>>   a: plot all Fourier components",
     $           ">>>      (Requires writing special draw*.in files ",
     $           ">>>       if dump file has >1 Fourier comp.)",
     $           ">>>   o: plot one Fourier component",
     $           ">>>?  "
            READ(nim_rd,*) flag
            IF (flag=='a'.OR.flag=='o') THEN
              EXIT
            ELSE IF (flag=='c') THEN
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $          ' may only be used for nonlinear cases .  Try again.'
            ELSE
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $         ' is not recognized.  Try again.'
            ENDIF
          CASE('correlation')
            WRITE(nim_wr,'(4(/a))',ADVANCE='NO')
     $           ">>> Data Options:",
     $           ">>>   a: sum over all Fourier components",
     $           ">>>   o: correlation from one Fourier component",
     $           ">>>?  "
            READ(nim_rd,*) flag
            IF (flag=='a'.OR.flag=='o') THEN
              EXIT
            ELSE
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $         ' is not recognized.  Try again.'
            ENDIF
          END SELECT
        ELSE
          SELECT CASE(pl_type)
          CASE ('nonlinear')
            IF (nfiles==1.AND.geom=='tor') THEN
              WRITE(nim_wr,'(11(/a))',ADVANCE='NO')
     $             ">>> Data Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>   p: plot sets of evenly spaced data in",
     $             ">>>      the periodic coordinate.",
     $             ">>>   t: same as p, but transform vector",
     $             ">>>      components to cartesian for toroidal",
     $             ">>>      3D plots (not xdraw).",
     $             ">>>   <return>: return to plot option level",
     $             ">>>?  "
              READ(nim_rd,'(a)') flag
              num_chars=LEN_TRIM(flag)
              IF (num_chars<1) THEN
                exit_loop=.true.
                EXIT datar
              ENDIF
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c'.OR.
     $            flag=='p'.OR.flag=='t') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ELSE IF (nfiles==1) THEN
              WRITE(nim_wr,'(8(/a))',ADVANCE='NO')
     $             ">>> Data Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>   p: plot sets of evenly spaced data in",
     $             ">>>      the periodic coordinate.",
     $             ">>>   <return>: return to plot option level",
     $             ">>>?  "
              READ(nim_rd,'(a)') flag
              num_chars=LEN_TRIM(flag)
              IF (num_chars<1) THEN
                exit_loop=.true.
                EXIT datar
              ENDIF
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c'.OR.
     $            flag=='p') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ELSE
              WRITE(nim_wr,'(6(/a))',ADVANCE='NO')
     $             ">>> Options:",
     $             ">>>   a: plot all Fourier components",
     $             ">>>   o: plot one Fourier component",
     $             ">>>   c: plot point in configuration space",
     $             ">>>   <return>: return to plot option level",
     $             ">>>?  "
              READ(nim_rd,'(a)') flag
              num_chars=LEN_TRIM(flag)
              IF (num_chars<1) THEN
                exit_loop=.true.
                EXIT datar
              ENDIF
              IF (flag=='a'.OR.flag=='o'.OR.flag=='c') THEN
                EXIT
              ELSE
                WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $            ' is not recognized.  Try again.'
              ENDIF
            ENDIF
          CASE('linear')
            WRITE(nim_wr,'(5(/a))',ADVANCE='NO')
     $           ">>> Options:",
     $           ">>>   a: plot all Fourier components",
     $           ">>>   o: plot one Fourier component",
     $           ">>>   <return>: return to plot option level",
     $           ">>>?  "
            READ(nim_rd,'(a)') flag
            num_chars=LEN_TRIM(flag)
            IF (num_chars<1) THEN
              exit_loop=.true.
              EXIT datar
            ENDIF
            IF (flag=='a'.OR.flag=='o') THEN
              EXIT
            ELSE IF (flag=='c') THEN
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $          ' may only be used for nonlinear cases .  Try again.'
            ELSE
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $          ' is not recognized.  Try again.'
            ENDIF
          CASE('correlation')
            WRITE(nim_wr,'(5(/a))',ADVANCE='NO')
     $           ">>> Data Options:",
     $           ">>>   a: sum over all Fourier components",
     $           ">>>   o: correlation from one Fourier component",
     $           ">>>   <return>: return to plot option level",
     $           ">>>?  "
            READ(nim_rd,'(a)') flag
            num_chars=LEN_TRIM(flag)
            IF (num_chars<1) THEN
              exit_loop=.true.
              EXIT datar
            ENDIF
            IF (flag=='a'.OR.flag=='o') THEN
              EXIT
            ELSE
              WRITE(nim_wr,*) '>>>  ', "'",flag,"'",
     $          ' is not recognized.  Try again.'
            ENDIF
          END SELECT
        ENDIF
      ENDDO datar

      RETURN
      END SUBROUTINE data_rep
c-----------------------------------------------------------------------
c     subprogram 3. acquire_dump_name
c     get a single dump file name.
c-----------------------------------------------------------------------
      SUBROUTINE acquire_dump_name(dump_data,first)
      USE local
      IMPLICIT NONE

      CHARACTER(128), INTENT(INOUT) :: dump_data
      LOGICAL, INTENT(INOUT) :: first

      LOGICAL :: file_stat
      INTEGER :: read_stat
      INTEGER(i4) :: num_chars
      CHARACTER(128) :: ctmp
c-----------------------------------------------------------------------
c     query and check dump file status.
c-----------------------------------------------------------------------
      DO
        IF (first) THEN
          WRITE(nim_wr,'(/a,/a)',ADVANCE='NO')
     $      '>> Please enter the name of the dump file to be read.',
     $      '>>? '
          READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $         SIZE=num_chars) ctmp
          IF (num_chars>0) THEN
            dump_data=ADJUSTL(ctmp)
            first=.false.
          ELSE
            CYCLE
          ENDIF
        ELSE
          WRITE(nim_wr,'(/a,/2a,/a)',ADVANCE='NO')
     $      '>> Please enter the name of the dump file to be read,',
     $      '>> or hit return to read from ',TRIM(dump_data),
     $      '>>? '
          READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $         SIZE=num_chars) ctmp
          IF (num_chars>0) dump_data=ADJUSTL(ctmp)
        ENDIF
        INQUIRE(FILE=TRIM(dump_data),EXIST=file_stat)
        IF (file_stat) EXIT
        WRITE(nim_wr,'(/,3a)')
     $    '>>  File ',TRIM(dump_data), ' does not exist.'
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE acquire_dump_name
c-----------------------------------------------------------------------
c     subprogram 4. acquire_dump_list
c     get an entire dump file list.
c-----------------------------------------------------------------------
      SUBROUTINE acquire_dump_list(file_list,first,nfiles,nfiles_old)
      USE local
      IMPLICIT NONE

      CHARACTER(128), DIMENSION(1000), INTENT(INOUT) :: file_list
      LOGICAL, INTENT(INOUT) :: first
      INTEGER(i4), INTENT(INOUT) :: nfiles_old
      INTEGER(i4), INTENT(OUT) :: nfiles

      LOGICAL :: file_stat
      INTEGER :: read_stat
      INTEGER(i4) :: num_chars
      CHARACTER(128) :: ctmp
c-----------------------------------------------------------------------
c     query and check each dump file status.
c-----------------------------------------------------------------------
      IF (first) THEN
        WRITE(nim_wr,'(4(/a),/4a,2(/a))')
     $    '>> Please enter a list of dump file names, each',
     $    '>> followed by a carriage return. Hit an additional',
     $    '>> carriage return when the list is complete.',
     $    '>>',
     $    '>> For a long list, suspend nimplot, type ',"'",
     $    'ls -1 dump.*',"'",
     $    '>> mark the list, put nimplot back in foreground,',
     $    '>> then paste the list, and hit return.'
      ELSE
        WRITE(nim_wr,'(5(/a))')
     $    '>> Please enter a list of dump file names, each',
     $    '>> followed by a carriage return. Hit an additional',
     $    '>> carriage return when the list is complete.',
     $    '>> If no files are entered, the previous file list',
     $    '>> is used.'
      ENDIF
      nfiles=0
      file_list_read: DO
        WRITE(nim_wr,'(a)',ADVANCE='NO') '>>? '
        READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $       SIZE=num_chars) ctmp
        IF (num_chars<1) THEN
          IF (nfiles>0.OR..NOT.first) THEN
            EXIT file_list_read
          ELSE
            CYCLE file_list_read
          ENDIF
        ELSE
          file_list(nfiles+1)=ADJUSTL(ctmp)
        ENDIF
        INQUIRE(FILE=TRIM(file_list(nfiles+1)),EXIST=file_stat)
        IF (.NOT.file_stat) THEN
          WRITE(nim_wr,*)
          WRITE(nim_wr,*) '  File ',TRIM(file_list(nfiles+1)),
     $      ' does not exist.  Please continue.'
        ELSE
          nfiles=nfiles+1
        ENDIF
      ENDDO file_list_read
      IF (.NOT.first.AND.nfiles==0) nfiles=nfiles_old
      first=.false.
      nfiles_old=nfiles
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE acquire_dump_list

