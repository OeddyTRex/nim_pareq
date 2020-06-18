c-----------------------------------------------------------------------
c     file tecplot_conv
c     contains a module with subprograms for converting xdraw graph
c     and xdraw contour plots to ascii tecplot format.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  tecplot_mod
c     1.  slice_conv
c     2.  contour_conv
c     3.  read_draw
c     4.  tec_eldata_write
c-----------------------------------------------------------------------
c     module tecplot_mod
c-----------------------------------------------------------------------
      MODULE tecplot_mod
      USE local
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1.  slice_conv
c     converts graph-type xdraw binary information to IJ-Ordered
c     POINT data for tecplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE slice_conv(binfile,pltfile,use_draw,varnames)

      CHARACTER(*), INTENT(IN) :: binfile,pltfile
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
      LOGICAL, INTENT(IN) :: use_draw

      INTEGER :: read_stat
      INTEGER(i4) :: ix,nx,nxlim=200,g1=0,g2=0,
     $               irec,nrec,nreclim=50000,ibl,ndbl
      REAL(r4), DIMENSION(:), ALLOCATABLE :: data
      REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: data_bl
      LOGICAL, DIMENSION(:), ALLOCATABLE :: isuitable,jsuitable
c-----------------------------------------------------------------------
c     determine record size.
c-----------------------------------------------------------------------
      DO nx=1,nxlim+1
        ALLOCATE(data(nx))
        CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
        READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
        CALL close_bin(xdr_unit,TRIM(binfile))
        DEALLOCATE(data)
        IF (read_stat>0) THEN
          EXIT     !     one beyond record length
        ELSE IF (read_stat<0) THEN
          WRITE(nim_wr,*) 'Prematurely hit the end of file.'
          STOP
        ENDIF
      ENDDO

      IF (nx>nxlim) THEN
        WRITE(nim_wr,*) 'Record size too large; increase nxlim.'
        STOP
      ELSE IF (nx==1) THEN
        WRITE(nim_wr,*) 'File '//TRIM(binfile)//' is empty???  (nx=0)'
        STOP
      ENDIF
      nx=nx-1

      IF (use_draw) THEN
        IF (nx/=SIZE(varnames)) THEN
          WRITE(nim_wr,*)
     $      'Variable list does not match with record size.'
          STOP
        ENDIF
      ENDIF

      WRITE(nim_wr,*) 'Record length: ',nx
      ALLOCATE(data(nx))
c-----------------------------------------------------------------------
c     determine number of records per data block.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
      DO nrec=1,nreclim+1
        READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
        IF (read_stat/=0) EXIT
      ENDDO
      CALL close_bin(xdr_unit,TRIM(binfile))
      IF (nrec>nreclim) THEN
        WRITE(nim_wr,*) 'Too many records; increase nreclim.'
        STOP
      ELSE IF (nrec==1) THEN
        WRITE(nim_wr,*) 'File '//TRIM(binfile)//' is empty???  (nrec=0)'
        STOP
      ENDIF
      nrec=nrec-1

      WRITE(nim_wr,*) 'Number of records: ',nrec
c-----------------------------------------------------------------------
c     determine number of data blocks.
c     also find an appropriate variables for I and J data.
c-----------------------------------------------------------------------
      ALLOCATE(data_bl(nx,nrec,3))
      ALLOCATE(isuitable(nx),jsuitable(nx))
      isuitable=.true.
      jsuitable=.true.
      CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
      ndbl=0
      ibl=1
      DO
        READ(xdr_unit,IOSTAT=read_stat) data_bl(1:nx,1,ibl)
        IF (read_stat==0) THEN
          ndbl=ndbl+1
        ELSE
          EXIT
        ENDIF
        DO irec=2,nrec
          READ(xdr_unit,IOSTAT=read_stat) data_bl(1:nx,irec,ibl)
          IF (read_stat/=0) THEN
            WRITE(nim_wr,*) 'Record block ',ndbl,' incomplete.'
            STOP
          ENDIF
          IF (irec>2) THEN
            DO ix=1,nx
              IF ((data_bl(ix,irec  ,ibl)-data_bl(ix,irec-1,ibl))*
     $            (data_bl(ix,irec-1,ibl)-data_bl(ix,irec-2,ibl))<=0)
     $          isuitable(ix)=.false.
            ENDDO
          ENDIF
        ENDDO
        IF (ibl==3) THEN
          DO ix=1,nx
            DO irec=1,nrec
              IF ((data_bl(ix,irec,3)-data_bl(ix,irec,2))*
     $            (data_bl(ix,irec,2)-data_bl(ix,irec,1))<=0)
     $          jsuitable(ix)=.false.
            ENDDO
          ENDDO
          data_bl=EOSHIFT(data_bl,1,DIM=3)
        ENDIF
        ibl=MIN(ibl+1,3_i4)
        READ(xdr_unit,IOSTAT=read_stat)
        IF (read_stat/=0) EXIT
      ENDDO

      DEALLOCATE(data_bl)
      CALL close_bin(xdr_unit,TRIM(binfile))

      WRITE(nim_wr,*) 'Number of data blocks: ',ndbl
c-----------------------------------------------------------------------
c     acquire variable names.
c-----------------------------------------------------------------------
      IF (.NOT.use_draw) THEN
        ALLOCATE(varnames(nx))
        WRITE(nim_wr,*) 'Enter the names of the ',nx,' variables ',
     $    '(from the appropriate draw*.in file).'
        READ(nim_rd,*) varnames
      ENDIF
c-----------------------------------------------------------------------
c     select which variables are appropriate for the grid.
c-----------------------------------------------------------------------
      IF (nrec/=1) THEN
        DO g1=1,nx
          IF (isuitable(g1)) EXIT
        ENDDO
      ENDIF
      IF (ndbl/=1) THEN
        DO g2=1,nx
          IF (jsuitable(g2).AND.g2/=g1) EXIT
        ENDDO
      ENDIF

      IF (ndbl==1) THEN
        WRITE(nim_wr,*) TRIM(varnames(g1)),
     $    ' will be used as default I data.'
      ELSE IF (nrec==1) THEN
        WRITE(nim_wr,*) TRIM(varnames(g2)),
     $    ' will be used as default I data.'
        g1=g2
        g2=0
      ELSE
        WRITE(nim_wr,*) TRIM(varnames(g1)),
     $    ' will be used as default I data.'
        WRITE(nim_wr,*) TRIM(varnames(g2)),
     $    ' will be used as default J data.'
      ENDIF

      DEALLOCATE(isuitable,jsuitable)
c-----------------------------------------------------------------------
c     open files.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
      OPEN(UNIT=out_unit,FILE=TRIM(pltfile),STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     plt-file header information.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(a)',ADVANCE="NO") "VARIABLES="

      WRITE(out_unit,'(x,3a)',ADVANCE="NO") '"',TRIM(varnames(g1)),'"'
      IF (g2/=0) 
     $  WRITE(out_unit,'(x,3a)',ADVANCE="NO") '"',TRIM(varnames(g2)),'"'
      DO ix=1,nx
        IF (ix==g1.OR.ix==g2) CYCLE
        WRITE(out_unit,'(x,3a)',ADVANCE="NO") '"',TRIM(varnames(ix)),'"'
      ENDDO
      WRITE(out_unit,'(x)')

      IF (ndbl==1) THEN
        WRITE(out_unit,'(a,i4,a)') "ZONE I=",nrec," F=POINT"
      ELSE IF (nrec==1) THEN
        WRITE(out_unit,'(a,i4,a)') "ZONE I=",ndbl," F=POINT"
      ELSE
        WRITE(out_unit,'(2(a,i4),a)')
     $    "ZONE I=",nrec," J=",ndbl," F=POINT"
      ENDIF
c-----------------------------------------------------------------------
c     read and write data.
c-----------------------------------------------------------------------
      DO ibl=1,ndbl
        DO irec=1,nrec
          READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
          IF (g2==0) THEN
            WRITE(out_unit,*) data(g1),data(:g1-1),data(g1+1:)
          ELSE
            WRITE(out_unit,*) data(g1),data(g2),data(:MIN(g1,g2)-1),
     $        data(MIN(g1,g2)+1:MAX(g1,g2)-1),data(MAX(g1,g2)+1:)
          ENDIF
        ENDDO
        READ(xdr_unit,IOSTAT=read_stat)
      ENDDO
c-----------------------------------------------------------------------
c     deallocate arrays, close files, and terminate.
c-----------------------------------------------------------------------
      DEALLOCATE(data)
      DEALLOCATE(varnames)
      CALL close_bin(xdr_unit,TRIM(binfile))
      CLOSE(out_unit)

      RETURN
      END SUBROUTINE slice_conv

c-----------------------------------------------------------------------
c     subprogram 2. contour_conv
c     converts contour-type xdraw binary information to FEPOINT
c     data for tecplot.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE contour_conv(binfile,pltfile,use_draw,varnames,
     $                        nsets,setname,setlabel,cart_conv)

      CHARACTER(*), INTENT(IN) :: binfile,pltfile,setname,cart_conv
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
      LOGICAL, INTENT(IN) :: use_draw
      INTEGER(i4), INTENT(IN) :: nsets
      REAL(r8), DIMENSION(:), POINTER :: setlabel

      CHARACTER(64) :: format,tmpdir
      CHARACTER(64), DIMENSION(:), ALLOCATABLE :: zonefile
      INTEGER :: read_stat
      INTEGER(i4) :: ibl,nrbl,ntbl,nbl,nx
      INTEGER(i4) :: icell,ix,ivert
      INTEGER(i4) :: iset
      LOGICAL :: cartesian

      TYPE :: bl_type
        CHARACTER(16) :: el_type
        INTEGER(i4) :: mx
        INTEGER(i4) :: my
        INTEGER(i4) :: mvert
        INTEGER(i4) :: mcell
        INTEGER(i4), DIMENSION(:,:), POINTER :: vertex
        REAL(r4), DIMENSION(:,:), POINTER :: rz
        REAL(r4), DIMENSION(:,:), POINTER :: data
      END TYPE bl_type

      TYPE(bl_type), DIMENSION(:), POINTER :: bl
c-----------------------------------------------------------------------
c     convenience parameters.
c-----------------------------------------------------------------------
      IF (cart_conv=='y') THEN
        cartesian=.true.
      ELSE
        cartesian=.false.
      ENDIF
c-----------------------------------------------------------------------
c     read block sizes and number of quantities.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
      READ(xdr_unit) nrbl,ntbl,nx
      IF (use_draw) THEN
        IF (nx+2/=SIZE(varnames)) THEN
          WRITE(nim_wr,*)
     $     'Variable list does not match with record size.'
          STOP
        ENDIF
      ENDIF
      nbl=nrbl+ntbl
      ALLOCATE(bl(nbl))
c-----------------------------------------------------------------------
c     read rblock grid.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        READ(xdr_unit) bl(ibl)%mx,bl(ibl)%my
        bl(ibl)%mvert=(bl(ibl)%mx+1)*(bl(ibl)%my+1)
        bl(ibl)%mcell=bl(ibl)%mx*bl(ibl)%my
c-TMP
c       WRITE(nim_wr,*) 'block ',ibl
c       WRITE(nim_wr,*) 'mx ',bl(ibl)%mx
c       WRITE(nim_wr,*) 'my ',bl(ibl)%my
c       WRITE(nim_wr,*) 'mvert ',bl(ibl)%mvert
c       WRITE(nim_wr,*) 'mcell ',bl(ibl)%mcell
        ALLOCATE(bl(ibl)%rz(bl(ibl)%mvert,2))
        READ(xdr_unit) bl(ibl)%rz
        IF (nsets==1) THEN
          bl(ibl)%el_type="QUADRILATERAL"
        ELSE
          bl(ibl)%el_type="BRICK"
        ENDIF
        ALLOCATE(bl(ibl)%vertex(bl(ibl)%mcell,4))
        DO icell=1,bl(ibl)%mcell
          bl(ibl)%vertex(icell,1)=((icell-1)/bl(ibl)%mx)+icell
          bl(ibl)%vertex(icell,2)=bl(ibl)%vertex(icell,1)+1
          bl(ibl)%vertex(icell,3)=bl(ibl)%vertex(icell,1)+bl(ibl)%mx+2
          bl(ibl)%vertex(icell,4)=bl(ibl)%vertex(icell,1)+bl(ibl)%mx+1
        ENDDO 
      ENDDO
c-----------------------------------------------------------------------
c     read tblock grid.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        READ(xdr_unit) bl(ibl)%mvert,bl(ibl)%mcell
        bl(ibl)%mvert=bl(ibl)%mvert+1
c-TMP
c       WRITE(nim_wr,*) 'block ',ibl
c       WRITE(nim_wr,*) 'mvert ',bl(ibl)%mvert
c       WRITE(nim_wr,*) 'mcell ',bl(ibl)%mcell
        ALLOCATE(bl(ibl)%rz(bl(ibl)%mvert,2))
        READ(xdr_unit) bl(ibl)%rz
        IF (nsets==1) THEN
          bl(ibl)%el_type="TRIANGLE"
        ELSE
          bl(ibl)%el_type="BRICK"
        ENDIF
        ALLOCATE(bl(ibl)%vertex(bl(ibl)%mcell,3))
        READ(xdr_unit) bl(ibl)%vertex
        bl(ibl)%vertex=bl(ibl)%vertex+1
      ENDDO
c-----------------------------------------------------------------------
c     create a temporary directory.
c-----------------------------------------------------------------------
      tmpdir="bin2plt_"//pltfile//"_tmpdir"
      CALL system("mkdir "//TRIM(tmpdir))
      ALLOCATE(zonefile(nbl))
      DO ibl=1,nbl
        IF (ibl<=9) THEN
          WRITE(zonefile(ibl),'(a,i1)')TRIM(tmpdir)//"/"//"ZONE00",ibl
        ELSE IF (ibl<=99) THEN
          WRITE(zonefile(ibl),'(a,i2)')TRIM(tmpdir)//"/"//"ZONE0",ibl
        ELSE
          WRITE(zonefile(ibl),'(a,i3)')TRIM(tmpdir)//"/"//"ZONE",ibl
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     acquire variable names.
c-----------------------------------------------------------------------
      IF (.NOT.use_draw) THEN
        ALLOCATE(varnames(nx+2))
        WRITE(nim_wr,*) ' '
        WRITE(nim_wr,*) 'Enter the names of the ',nx+2,' variables ',
     $    '(from the appropriate draw*.in file).'
        READ(nim_rd,*) varnames
      ENDIF
c-----------------------------------------------------------------------
c     open files.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        OPEN(UNIT=con_unit+ibl,FILE=TRIM(zonefile(ibl)),
     $       STATUS='UNKNOWN')
      ENDDO
c-----------------------------------------------------------------------
c     plt-file header information.
c-----------------------------------------------------------------------
      WRITE(con_unit+1,'(a)',ADVANCE="NO") "VARIABLES="
      IF (cartesian) THEN
        WRITE(con_unit+1,'(3(x,3a))',ADVANCE="NO")
     $      '"','x','"','"','y','"','"','z','"'
      ELSE
        DO ix=1,2
          WRITE(con_unit+1,'(x,3a)',ADVANCE="NO")
     $      '"',TRIM(varnames(ix)),'"'
        ENDDO
        IF (nsets>1) WRITE(con_unit+1,'(x,3a)',ADVANCE="NO")
     $      '"',TRIM(setname),'"'
      ENDIF
      DO ix=3,nx+2
        WRITE(con_unit+1,'(x,3a)',ADVANCE="NO")
     $    '"',TRIM(varnames(ix)),'"'
      ENDDO
      WRITE(con_unit+1,'(x)')
c-----------------------------------------------------------------------
c     read and write data.
c-----------------------------------------------------------------------
      IF (nsets==1) THEN
        DO ibl=1,nbl
          WRITE(con_unit+ibl,'(x)')
          WRITE(con_unit+ibl,'(a,i7,a,i7,2a)') "ZONE N=",bl(ibl)%mvert,
     $      " E=",bl(ibl)%mcell," F=FEPOINT ET=",TRIM(bl(ibl)%el_type)
c-----------------------------------------------------------------------
c         data
c-----------------------------------------------------------------------
          ALLOCATE(bl(ibl)%data(bl(ibl)%mvert,nx))
          DO ix=1,nx
            READ(xdr_unit,IOSTAT=read_stat) bl(ibl)%data(:,ix)
            IF (read_stat/=0) THEN
              WRITE(nim_wr,*) 'Error reading quantity ',
     $          TRIM(varnames(ix+2)),' in block ',ibl
            ENDIF
          ENDDO
          DO ivert=1,bl(ibl)%mvert
            WRITE(con_unit+ibl,*) bl(ibl)%rz(ivert,:),
     $                            bl(ibl)%data(ivert,:)
          ENDDO
          DEALLOCATE(bl(ibl)%data)
c-----------------------------------------------------------------------
c         element connections
c-----------------------------------------------------------------------
          WRITE(con_unit+ibl,'(x)')
          DO icell=1,bl(ibl)%mcell
            WRITE(con_unit+ibl,*) bl(ibl)%vertex(icell,:)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       multi-set data written as 3D brick finite elements.
c-----------------------------------------------------------------------
      ELSE
        DO ibl=1,nbl
          WRITE(con_unit+ibl,'(x)')
          WRITE(con_unit+ibl,'(a,i7,a,i7,2a)') 
     $      "ZONE N=",bl(ibl)%mvert*nsets," E=",bl(ibl)%mcell*(nsets-1),
     $      " F=FEPOINT ET=",TRIM(bl(ibl)%el_type)
        ENDDO
c-----------------------------------------------------------------------
c       data
c-----------------------------------------------------------------------
        DO iset=1,nsets
          DO ibl=1,nbl
            ALLOCATE(bl(ibl)%data(bl(ibl)%mvert,nx))
            DO ix=1,nx
              READ(xdr_unit,IOSTAT=read_stat) bl(ibl)%data(:,ix)
              IF (read_stat/=0) THEN
                WRITE(nim_wr,*) 'Error reading quantity ',
     $            TRIM(varnames(ix+2)),' in block ',ibl
              ENDIF
            ENDDO
            IF (cartesian) THEN
              DO ivert=1,bl(ibl)%mvert
                WRITE(con_unit+ibl,*)
     $             bl(ibl)%rz(ivert,1)*COS(twopi*setlabel(iset)),
     $            -bl(ibl)%rz(ivert,1)*SIN(twopi*setlabel(iset)),
     $             bl(ibl)%rz(ivert,2),
     $             bl(ibl)%data(ivert,:)
              ENDDO
            ELSE
              DO ivert=1,bl(ibl)%mvert
                WRITE(con_unit+ibl,*) bl(ibl)%rz(ivert,:),
     $                                setlabel(iset),
     $                                bl(ibl)%data(ivert,:)
              ENDDO
            ENDIF
            DEALLOCATE(bl(ibl)%data)
          ENDDO
        ENDDO
        DEALLOCATE(setlabel)
c-----------------------------------------------------------------------
c       element connections
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          WRITE(con_unit+ibl,'(x)')
          DO iset=1,nsets-1
            IF (SIZE(bl(ibl)%vertex,2)==3) THEN
              DO icell=1,bl(ibl)%mcell
                WRITE(con_unit+ibl,*)
     $            bl(ibl)%vertex(icell,1)+bl(ibl)%mvert*(iset-1),
     $            bl(ibl)%vertex(icell,:)+bl(ibl)%mvert*(iset-1),
     $            bl(ibl)%vertex(icell,1)+bl(ibl)%mvert*(iset),
     $            bl(ibl)%vertex(icell,:)+bl(ibl)%mvert*(iset)
              ENDDO
            ELSE
              DO icell=1,bl(ibl)%mcell
                WRITE(con_unit+ibl,*)
     $            bl(ibl)%vertex(icell,:)+bl(ibl)%mvert*(iset-1),
     $            bl(ibl)%vertex(icell,:)+bl(ibl)%mvert*(iset)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close files
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CLOSE(con_unit+ibl)
      ENDDO
      CALL close_bin(xdr_unit,TRIM(binfile))
c-----------------------------------------------------------------------
c     concatenate separate zones into one plt file.
c-----------------------------------------------------------------------
      CALL system("cat "//TRIM(tmpdir)//"/ZONE* > "//TRIM(pltfile))
      CALL system("rm -r "//TRIM(tmpdir))
c-----------------------------------------------------------------------
c     deallocate arrays and terminate.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        DEALLOCATE(bl(ibl)%vertex)
        DEALLOCATE(bl(ibl)%rz)
        CLOSE(con_unit+ibl)
      ENDDO
      DEALLOCATE(bl)
      DEALLOCATE(varnames,zonefile)

      RETURN
      END SUBROUTINE contour_conv

c-----------------------------------------------------------------------
c     subprogram 3.  read_draw.
c     read an xdraw draw*.in file to acquire variable names.
c-----------------------------------------------------------------------
      SUBROUTINE read_draw(dfile,pl_type,varnames,tvar,binfile)

      CHARACTER(*), INTENT(IN) :: dfile
      CHARACTER(*), INTENT(INOUT) :: tvar
      CHARACTER(1), INTENT(OUT) :: pl_type
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
      CHARACTER(*), INTENT(OUT) :: binfile

      CHARACTER(128) :: fileline
      CHARACTER(64) :: tmpvar,xvar='x',yvar='y'
      INTEGER(i4) :: nlen,ipass,ivar,nvar,nivar,jvar
      INTEGER :: read_stat,numchars
c-----------------------------------------------------------------------
c     open draw*.in file, and determine plot type and the binary
c     file name in the draw*.in file.
c-----------------------------------------------------------------------
      OPEN(temp_unit,FILE=TRIM(dfile),STATUS="OLD",ACTION="READ")
      READ(temp_unit,'(a)') fileline
      nlen=LEN_TRIM(fileline)
      IF (fileline(nlen:nlen)=='G'.OR.fileline(nlen:nlen)=='g') THEN
        pl_type='g'
      ELSE
        pl_type='c'
      ENDIF
      DO
        READ(temp_unit,'(a)',IOSTAT=read_stat) fileline
        IF (read_stat/=0) THEN
          WRITE(nim_wr,*) 'Error reading file ',TRIM(dfile),'.'
          STOP
        ENDIF
        fileline=ADJUSTL(fileline)
        IF (fileline(1:8)=='filename'.OR.
     $      fileline(1:8)=='Filename'.OR.
     $      fileline(1:8)=='FILENAME') EXIT
      ENDDO
      DO
        READ(temp_unit,'(a)',IOSTAT=read_stat) fileline
        IF (read_stat/=0) THEN
          WRITE(nim_wr,*) 'Error reading file ',TRIM(dfile),'.'
          STOP
        ENDIF
        binfile=ADJUSTL(fileline)
        IF (binfile(1:1)/=';') EXIT
      ENDDO
      CLOSE (temp_unit)
c-----------------------------------------------------------------------
c     for graph-type plots, skip to the variable list.
c-----------------------------------------------------------------------
      IF (pl_type=='g') THEN
        gpasses: DO ipass=1,2
          OPEN(temp_unit,FILE=TRIM(dfile),STATUS="OLD",ACTION="READ")
          IF (ipass==2) ALLOCATE(varnames(nvar))
          DO 
            READ(temp_unit,'(a)',IOSTAT=read_stat) fileline
            IF (read_stat/=0) THEN
              WRITE(nim_wr,*) 'Error reading file ',TRIM(dfile),'.'
              STOP
            ENDIF
            fileline=ADJUSTL(fileline)
            IF (fileline(1:14)=='variable names'.OR.
     $          fileline(1:14)=='Variable Names'.OR.
     $          fileline(1:14)=='Variable names'.OR.
     $          fileline(1:14)=='VARIABLE NAMES') EXIT
          ENDDO
          nvar=0
          DO
            READ(temp_unit,'(a)',IOSTAT=read_stat,ADVANCE="NO",
     $           SIZE=numchars) fileline
            IF (numchars<1) THEN
              CYCLE
            ELSE IF (fileline(1:2)=='ix') THEN
              EXIT
            ELSE IF (ipass==2) THEN
              nvar=nvar+1
              READ(fileline,*) ivar,varnames(nvar)
            ELSE
              nvar=nvar+1
            ENDIF
          ENDDO
          CLOSE (temp_unit)
        ENDDO gpasses
c-----------------------------------------------------------------------
c     for contour-type plots, skip to the variable list.
c-----------------------------------------------------------------------
      ELSE
        cpasses: DO ipass=1,2
          OPEN(temp_unit,FILE=TRIM(dfile),STATUS="OLD",ACTION="READ")
          IF (ipass==2) THEN
            ALLOCATE(varnames(nvar+2))
            varnames(1:2)=(/xvar,yvar/)
          ENDIF
          DO 
            READ(temp_unit,'(a)',IOSTAT=read_stat) fileline
            IF (read_stat/=0) THEN
              WRITE(nim_wr,*) 'Error reading file ',TRIM(dfile),'.'
              STOP
            ENDIF
            fileline=ADJUSTL(fileline)
            IF (fileline(1:26)=='independent variable names'.OR.
     $          fileline(1:26)=='Independent Variable Names'.OR.
     $          fileline(1:26)=='Independent variable names'.OR.
     $          fileline(1:26)=='INDEPENDENT VARIABLE NAMES') EXIT
          ENDDO
          DO 
            READ(temp_unit,'(a)',IOSTAT=read_stat) fileline
            IF (read_stat/=0) THEN
              WRITE(nim_wr,*) 'Error reading file ',TRIM(dfile),'.'
              STOP
            ENDIF
            fileline=ADJUSTL(fileline)
            IF (fileline(1:1)=='t') THEN
              tvar=ADJUSTL(fileline(2:))
            ELSE IF (fileline(1:1)=='x') THEN
              xvar=ADJUSTL(fileline(2:))
            ELSE IF (fileline(1:1)=='y') THEN
              yvar=ADJUSTL(fileline(2:))
            ELSE IF (fileline(1:24)=='dependent variable names'.OR.
     $               fileline(1:24)=='Dependent Variable Names'.OR.
     $               fileline(1:24)=='Dependent variable names'.OR.
     $               fileline(1:24)=='DEPENDENT VARIABLE NAMES') THEN
              EXIT
            ENDIF
          ENDDO
          nvar=0
          DO
            READ(temp_unit,'(a)',IOSTAT=read_stat,ADVANCE="NO",
     $           SIZE=numchars) fileline
            IF (numchars<1) THEN
              CYCLE
            ELSE IF (fileline(1:4)=='iqty') THEN
              EXIT
            ELSE IF (ipass==2) THEN
              nvar=nvar+1
              READ(fileline,*) ivar,varnames(nvar+2)
            ELSE
              nvar=nvar+1
            ENDIF
          ENDDO
          CLOSE (temp_unit)
        ENDDO cpasses
      ENDIF
c-----------------------------------------------------------------------
c     check for repeated variable names.
c-----------------------------------------------------------------------
      DO ipass=1,SIZE(varnames)
        jvar=1
        DO ivar=ipass+1,SIZE(varnames)
          IF (varnames(ivar)==varnames(ipass)) THEN
            jvar=jvar+1
            IF (jvar<=9) THEN
              WRITE(varnames(ivar),'(2a,i1)') 
     $          TRIM(varnames(ivar)),' # ',jvar
            ELSE IF (jvar<=99) THEN
              WRITE(varnames(ivar),'(2a,i2)') 
     $          TRIM(varnames(ivar)),' # ',jvar
            ELSE
              WRITE(varnames(ivar),'(2a,i3)') 
     $          TRIM(varnames(ivar)),' # ',jvar
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_draw

c-----------------------------------------------------------------------
c     subprogram 4. tec_eldata_write
c     this subroutine writes fields that are discontinuous across
c     element borders in Tecplot's element format.  all of the required
c     information is passed in the r4eldata_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE tec_eldata_write(eldata,filename,unitno)
      USE eldata_type_mod

      TYPE(r4eldata_type), DIMENSION(:), INTENT(IN) :: eldata
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER(i4), INTENT(IN) :: unitno

      INTEGER(i4) :: npt,ipt,nbl,ibl,iel,nel,isub,ivar,nvar
c-----------------------------------------------------------------------
c     open the tecplot file.
c-----------------------------------------------------------------------
      OPEN(unitno,FILE=TRIM(filename),STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     each element is written as a separate tecplot zone to plot
c     discontinuous data.
c-----------------------------------------------------------------------
      nbl=SIZE(eldata)
      nvar=eldata(1)%nfields
      DO ibl=1,nbl
        DO iel=1,eldata(ibl)%nele
          WRITE(unitno,'(/5a)',ADVANCE="NO") ' VARIABLES= "',
     $      eldata(ibl)%coordlabels(1:1),'" "',
     $      eldata(ibl)%coordlabels(2:2),'"'
          DO ivar=1,nvar-1
            WRITE(unitno,'(3a)',ADVANCE="NO") ' "',
     $        TRIM(eldata(ibl)%datalabels(ivar)),'"'
          ENDDO
          WRITE(unitno,'(3a)') ' "',
     $      TRIM(eldata(ibl)%datalabels(nvar)),'"'
          WRITE(unitno,*) 'ZONE N=',eldata(ibl)%nptperel,
     $      ' E=',SIZE(eldata(ibl)%connect,2),
     $      ' F=FEPOINT ET=QUADRILATERAL'
          DO ipt=1,eldata(ibl)%nptperel
            WRITE(unitno,'(2es16.8)',ADVANCE="NO") 
     $        eldata(ibl)%coords(:,ipt,iel)
            DO ivar=1,nvar-1
              WRITE(unitno,'(es16.8)',ADVANCE="NO") 
     $          eldata(ibl)%data(ivar,ipt,iel)
            ENDDO
            WRITE(unitno,'(es16.8)') eldata(ibl)%data(nvar,ipt,iel)
          ENDDO
          DO isub=1,SIZE(eldata(ibl)%connect,2)
            WRITE(unitno,*) eldata(ibl)%connect(:,isub,iel)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CLOSE(unitno)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE tec_eldata_write

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE tecplot_mod
