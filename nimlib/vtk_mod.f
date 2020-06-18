c-----------------------------------------------------------------------
c     file vtk_mod
c     contains a module with subprograms for converting xdraw graph
c     and xdraw contour plots to the VTK legacy format, which can be
c     read by VisIt and VTK.  For more information on the format see
c     http://www.vtk.org/img/file-formats.pdf .
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  vtk_mod
c     1.  vtk_slice_conv
c     2.  vtk_contour_conv
c     3.  vtk_vec_check
c     4.  vtk_label_suffix
c     5.  vtk_eldata_write
c-----------------------------------------------------------------------
c     module vtk_mod
c-----------------------------------------------------------------------
      MODULE vtk_mod
      USE local
      IMPLICIT NONE

      CHARACTER(6), PARAMETER :: vtk_format="BINARY"  !  or "ASCII"
      CHARACTER(1), PARAMETER :: linef=char(10)  !  line feed

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1.  vtk_slice_conv
c     converts graph-type xdraw binary information to IJ-Ordered
c     POINT data for vtk.
c-----------------------------------------------------------------------
      SUBROUTINE vtk_slice_conv(binfile,pltfile,use_draw,varnames)

      CHARACTER(*), INTENT(IN) :: binfile,pltfile
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
      LOGICAL, INTENT(IN) :: use_draw

      INTEGER :: read_stat
      INTEGER(i4), PARAMETER :: nxlim=200,nreclim=50000
      INTEGER(i4) :: ix,nx,g1,g2,irec,nrec,ibl,ndbl,i1,i2,npt,nlab,
     $               ii,jj
      INTEGER(i4), DIMENSION(nxlim) :: iind
      REAL(r4), DIMENSION(:), ALLOCATABLE :: data
      REAL(r4), DIMENSION(:,:,:), ALLOCATABLE :: data_bl
      LOGICAL, DIMENSION(:), ALLOCATABLE :: isuitable,jsuitable
      CHARACTER(128) :: field_label
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

      WRITE(nim_wr,'(/,a,i5)') 'Record length: ',nx
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

      WRITE(nim_wr,'(a,i5)') 'Number of records: ',nrec
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

      WRITE(nim_wr,'(a,i5)') 'Number of data blocks: ',ndbl
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
c     for VisIt, have the user select one independent coordinate.
c-----------------------------------------------------------------------
      WRITE(nim_wr,'(/,a)')
     $  'Select the desired independent variable from the following:'
      i2=0
      iind=0
      DO i1=1,nx
        IF (isuitable(i1).AND.nrec/=1) THEN
          i2=i2+1
          iind(i2)=i1
          WRITE(nim_wr,'(i3,2a)') i2,') ',varnames(i1)
        ELSE IF (jsuitable(i1).AND.ndbl/=1) THEN
          i2=i2+1
          iind(i2)=-i1
          WRITE(nim_wr,'(i3,2a)') i2,') ',varnames(i1)
        ENDIF
      ENDDO
      WRITE(nim_wr,'(a)',ADVANCE='NO') 'Enter index >? '
      READ(nim_rd,*) g1
      i1=ABS(iind(g1))
      DEALLOCATE(isuitable,jsuitable)
c-----------------------------------------------------------------------
c     select another variable for the labeled index.
c-----------------------------------------------------------------------
      IF (iind(g1)>0) THEN 
        DO g2=1,nx
          IF (iind(g2)==0.OR.g2/=g1.AND.iind(g2)<0) EXIT
        ENDDO
      ELSE
        DO g2=1,nx
          IF (iind(g2)==0.OR.g2/=g1.AND.iind(g2)>0) EXIT
        ENDDO
      ENDIF
      i2=0
      IF (iind(g2)/=0) i2=ABS(iind(g2))
c-----------------------------------------------------------------------
c     open and read the xdraw file.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(binfile),'OLD','REWIND',32_i4)
      ALLOCATE(data_bl(nx,nrec,ndbl))
      DO ibl=1,ndbl
        DO irec=1,nrec
          READ(xdr_unit,IOSTAT=read_stat) data_bl(:,irec,ibl)
        ENDDO
        READ(xdr_unit,IOSTAT=read_stat)
      ENDDO
      CALL close_bin(xdr_unit,TRIM(binfile))
c-----------------------------------------------------------------------
c     vtk output may be ascii or binary.  write the vtk-file header
c     information in the appropriate format.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        OPEN(UNIT=out_unit,FILE=TRIM(pltfile),STATUS='UNKNOWN',
     $       ACCESS="STREAM")
        WRITE(out_unit) "# vtk DataFile Version 2.0"
        WRITE(out_unit) linef//"Plots over "//TRIM(varnames(i1))//
     $    " indexed by "//TRIM(varnames(i2))
        WRITE(out_unit) linef//"BINARY"
        WRITE(out_unit) linef//"DATASET RECTILINEAR_GRID"
      ELSE
        OPEN(UNIT=out_unit,FILE=TRIM(pltfile),STATUS='UNKNOWN')
        WRITE(out_unit,'(a)') "# vtk DataFile Version 2.0"
        WRITE(out_unit,'(4a)') "Plots over ",TRIM(varnames(i1)),
     $    " indexed by ",TRIM(varnames(i2))
        WRITE(out_unit,'(a,/,/,a)') "ASCII","DATASET RECTILINEAR_GRID"
      ENDIF
c-----------------------------------------------------------------------
c     the x-y plots in VisIt work best with one coordinate.  loop
c     and write over records in a block or over data blocks according
c     to the user's choice, and turn the other index into multiple
c     variable names.
c
c     give VisIt dummy variables for the other coordinates.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        IF (iind(g1)>=0) THEN  !  coordinate is over records
          WRITE(field_label,'(a,3i8)') "DIMENSIONS ",nrec,1,1
          WRITE(out_unit) linef//TRIM(field_label)
          WRITE(field_label,'(a,i8,a)') "X_COORDINATES ",nrec," float"
          WRITE(out_unit) linef//TRIM(field_label)//linef
          WRITE(out_unit) data_bl(i1,1:nrec,1)
          npt=nrec
          nlab=ndbl
        ELSE                  !  coordinate is over data blocks
          WRITE(field_label,'(a,3i8)') "DIMENSIONS ",ndbl,1,1
          WRITE(out_unit) linef//TRIM(field_label)
          WRITE(field_label,'(a,i8,a)') "X_COORDINATES ",ndbl," float"
          WRITE(out_unit) linef//TRIM(field_label)//linef
          WRITE(out_unit) data_bl(i1,1,1:ndbl)
          npt=nrec
          npt=ndbl
          nlab=nrec
        ENDIF
        WRITE(out_unit) linef//"Y_COORDINATES 1 float"//linef
        WRITE(out_unit) 0._r4
        WRITE(out_unit) linef//"Z_COORDINATES 1 float"//linef
        WRITE(out_unit) 0._r4
      ELSE
        IF (iind(g1)>=0) THEN  !  coordinate is over records
          WRITE(out_unit,'(a,3i8)') "DIMENSIONS ",nrec,1,1
          WRITE(out_unit,'(a,i8,a)') "X_COORDINATES ",nrec," double"
          DO irec=1,nrec
            WRITE(out_unit,*) data_bl(i1,irec,1)
          ENDDO
          npt=nrec
          nlab=ndbl
        ELSE                  !  coordinate is over data blocks
          WRITE(out_unit,'(a,3i8)') "DIMENSIONS ",ndbl,1,1
          WRITE(out_unit,'(a,i8,a)') "X_COORDINATES ",ndbl," double"
          DO ibl=1,ndbl
            WRITE(out_unit,*) data_bl(i1,1,ibl)
          ENDDO
          npt=ndbl
          nlab=nrec
        ENDIF
        WRITE(out_unit,'(2(/,a,es9.2))') "Y_COORDINATES 1 double",0._r4,
     $                                   "Z_COORDINATES 1 double",0._r4
      ENDIF
c-----------------------------------------------------------------------
c     loop over the outer label index, then loop over the coordinate
c     index.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(field_label,'(a,i9)') "POINT_DATA ",npt
        WRITE(out_unit) linef//TRIM(field_label)//linef
      ELSE
        WRITE(out_unit,'(/,a,i9)') "POINT_DATA ",npt
      ENDIF
      DO jj=1,nlab
        DO ii=1,nx
          IF (ii==i1.OR.ii==i2) CYCLE
          CALL vtk_label_suffix(varnames(ii),.false.,field_label)
          SELECT CASE (jj)
          CASE(1:9)
            WRITE(field_label,'(4a,i1)') TRIM(field_label),'_',
     $        TRIM(varnames(i2)),'00',jj
          CASE(10:99)
            WRITE(field_label,'(4a,i2)') TRIM(field_label),'_',
     $        TRIM(varnames(i2)),'0',jj
          CASE DEFAULT
            WRITE(field_label,'(3a,i3)') TRIM(field_label),'_',
     $        TRIM(varnames(i2)),jj
          END SELECT
          CALL vtk_data_write(ii,iind(g1)>=0,jj)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     deallocate arrays, close vtk file, and terminate.
c-----------------------------------------------------------------------
      CLOSE(out_unit)
      DEALLOCATE(data,data_bl)
      DEALLOCATE(varnames)
      RETURN

      CONTAINS
c-----------------------------------------------------------------------
c       internal subroutine for writing one field.
c-----------------------------------------------------------------------
        SUBROUTINE vtk_data_write(ifld,recwrite,jfld)

        INTEGER(i4), INTENT(IN) :: ifld,jfld
        LOGICAL, INTENT(IN) :: recwrite

        IF (vtk_format=="BINARY") THEN
          WRITE(out_unit) linef//"SCALARS "//TRIM(field_label)//
     $      " float 1"
          WRITE(out_unit) linef//"LOOKUP_TABLE default"//linef
          IF (recwrite) THEN ! writing index is over records in a block
            WRITE(out_unit) data_bl(ifld,1:nrec,jfld)
          ELSE  !  writing index is over blocks for a given record
            WRITE(out_unit) data_bl(ifld,jfld,1:ndbl)
          ENDIF
        ELSE
          WRITE(out_unit,'(/,a)')"SCALARS "//TRIM(field_label)//
     $      " double 1"
          WRITE(out_unit,'(a)') "LOOKUP_TABLE default"
          IF (recwrite) THEN ! writing index is over records in a block
            DO irec=1,nrec
              WRITE(out_unit,*) data_bl(ifld,irec,jfld)
            ENDDO
          ELSE  !  writing index is over blocks for a given record
            DO ibl=1,ndbl
              WRITE(out_unit,*) data_bl(ifld,jfld,ibl)
            ENDDO
          ENDIF
        ENDIF

        RETURN
        END SUBROUTINE vtk_data_write

      END SUBROUTINE vtk_slice_conv

c-----------------------------------------------------------------------
c     subprogram 2. vtk_contour_conv
c     converts contour-type xdraw binary information to VTK's
c     UNSTRUCTURED_GRID format.  It is similar to Tecplot's FEPOINT,
c     except that all points and elements are defined before data is
c     written, i.e. there are no separate zones.
c-----------------------------------------------------------------------
      SUBROUTINE vtk_contour_conv(binfile,pltfile,use_draw,varnames,
     $                            nsets,setname,setlabel,cart_conv)

      CHARACTER(*), INTENT(IN) :: binfile,pltfile,setname,cart_conv
      CHARACTER(64), DIMENSION(:), POINTER :: varnames
      LOGICAL, INTENT(IN) :: use_draw
      INTEGER(i4), INTENT(IN) :: nsets
      REAL(r8), DIMENSION(:), POINTER :: setlabel

      CHARACTER(1) :: c1,c2,c3
      CHARACTER(8) :: suffix
      CHARACTER(72) :: field_label
      CHARACTER(128), DIMENSION(SIZE(varnames)) :: header
      INTEGER :: read_stat
      INTEGER(i4) :: ibl,nrbl,ntbl,nbl,nx,npt,ncell,ncdata
      INTEGER(i4) :: icell,ix,ivert,ioff,nfield,inam
      INTEGER(i4) :: iset,ifield,narr
      INTEGER(i4), PARAMETER :: nline=1_i4
      LOGICAL :: cartesian,vector
      LOGICAL, DIMENSION(SIZE(varnames)) :: veclist
      REAL(r4), DIMENSION(nsets) :: setl4
      REAL(r4) :: twopi4

      TYPE :: bl_type
        INTEGER(i4) :: el_type  !  integer unlike Tecplot
        INTEGER(i4) :: mx
        INTEGER(i4) :: my
        INTEGER(i4) :: mvert
        INTEGER(i4) :: mcell
        INTEGER(i4), DIMENSION(:,:), POINTER :: vertex
        REAL(r4), DIMENSION(:,:), POINTER :: rz
        REAL(r4), DIMENSION(:,:), POINTER :: data
      END TYPE bl_type

      TYPE(bl_type), DIMENSION(:,:), POINTER :: bl
c-----------------------------------------------------------------------
c     convenience parameters.
c-----------------------------------------------------------------------
      IF (cart_conv=='y') THEN
        cartesian=.true.
      ELSE
        cartesian=.false.
      ENDIF
      IF (ASSOCIATED(setlabel)) THEN
        setl4=setlabel(1:nsets)
      ELSE
        setl4=0._r4
      ENDIF
      twopi4=twopi
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
      ALLOCATE(bl(nbl,nsets))
c-----------------------------------------------------------------------
c     read rblock grid into the first-set index.  the data is already
c     organized as 4-node quads, so the polynomial bases are disguised.
c     VTK indexes vertices from 0.
c-----------------------------------------------------------------------
      npt=0
      ncell=0
      ncdata=0
      DO ibl=1,nrbl
        READ(xdr_unit) bl(ibl,1)%mx,bl(ibl,1)%my
        bl(ibl,1)%mvert=(bl(ibl,1)%mx+1)*(bl(ibl,1)%my+1)
        bl(ibl,1)%mcell=bl(ibl,1)%mx*bl(ibl,1)%my
        ALLOCATE(bl(ibl,1)%rz(bl(ibl,1)%mvert,2))
        READ(xdr_unit) bl(ibl,1)%rz
        IF (nsets==1) THEN
          bl(ibl,1)%el_type=9_i4  !  VTK_QUAD
        ELSE
          bl(ibl,1)%el_type=12_i4 !  VTK_HEXAHEDRON
        ENDIF
        ALLOCATE(bl(ibl,1)%vertex(bl(ibl,1)%mcell,4))
        DO icell=1,bl(ibl,1)%mcell
          bl(ibl,1)%vertex(icell,1)=
     $      ((icell-1)/bl(ibl,1)%mx)+icell+npt-1_i4
          bl(ibl,1)%vertex(icell,2)=
     $      bl(ibl,1)%vertex(icell,1)+1
          bl(ibl,1)%vertex(icell,3)=
     $      bl(ibl,1)%vertex(icell,1)+bl(ibl,1)%mx+2
          bl(ibl,1)%vertex(icell,4)=
     $      bl(ibl,1)%vertex(icell,1)+bl(ibl,1)%mx+1
          IF (nsets==1) THEN
            ncdata=ncdata+5
          ELSE
            ncdata=ncdata+9*(nsets-1)
          ENDIF
        ENDDO 
        npt=npt+bl(ibl,1)%mvert
        ncell=ncell+bl(ibl,1)%mcell
      ENDDO
c-----------------------------------------------------------------------
c     read tblock grid into the first-set index.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        READ(xdr_unit) bl(ibl,1)%mvert,bl(ibl,1)%mcell
        bl(ibl,1)%mvert=bl(ibl,1)%mvert+1
        ALLOCATE(bl(ibl,1)%rz(bl(ibl,1)%mvert,2))
        READ(xdr_unit) bl(ibl,1)%rz
        IF (nsets==1) THEN
          bl(ibl,1)%el_type=5_i4  !  VTK_TRIANGLE
          ncdata=ncdata+bl(ibl,1)%mcell*4
        ELSE
          bl(ibl,1)%el_type=13_i4 !  VTK_WEDGE
          ncdata=ncdata+bl(ibl,1)%mcell*7*(nsets-1)
        ENDIF
        ALLOCATE(bl(ibl,1)%vertex(bl(ibl,1)%mcell,3))
        READ(xdr_unit) bl(ibl,1)%vertex
        bl(ibl,1)%vertex=bl(ibl,1)%vertex+npt
        npt=npt+bl(ibl,1)%mvert
        ncell=ncell+bl(ibl,1)%mcell
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
c     vtk output may be ascii or binary.  write the vtk-file header
c     information in the appropriate format.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        OPEN(UNIT=con_unit,FILE=TRIM(pltfile),STATUS='UNKNOWN',
     $       ACCESS="STREAM")
        WRITE(con_unit) "# vtk DataFile Version 2.0"
        IF (nsets==1) THEN
          WRITE(con_unit) linef//"Unstructured 2D Grid and Data"
        ELSE
          WRITE(con_unit) linef//"Unstructured 3D Grid and Data"
        ENDIF
        WRITE(con_unit) linef//"BINARY"
        WRITE(con_unit) linef//"DATASET UNSTRUCTURED_GRID"
        WRITE(field_label,'(a,i9,a)') "POINTS ",npt*nsets," float"
        WRITE(con_unit) linef//TRIM(field_label)//linef
      ELSE
        OPEN(UNIT=con_unit,FILE=TRIM(pltfile),STATUS='UNKNOWN')
        WRITE(con_unit,'(a)') "# vtk DataFile Version 2.0"
        IF (nsets==1) THEN
          WRITE(con_unit,'(a)') "Unstructured 2D Grid and Data"
        ELSE
          WRITE(con_unit,'(a)') "Unstructured 3D Grid and Data"
        ENDIF
        WRITE(con_unit,'(a,/,/,a)') "ASCII","DATASET UNSTRUCTURED_GRID"
        WRITE(con_unit,'(/,a,i9,a)') "POINTS ",npt*nsets," double"
      ENDIF
c-----------------------------------------------------------------------
c     write the coordinates of each data point.
c-----------------------------------------------------------------------
      DO iset=1,nsets
        DO ibl=1,nbl
          IF (cartesian) THEN
            DO ivert=1,bl(ibl,1)%mvert
              IF (vtk_format=="BINARY") THEN
                WRITE(con_unit)
     $             bl(ibl,1)%rz(ivert,1)*COS(twopi4*setl4(iset)),
     $            -bl(ibl,1)%rz(ivert,1)*SIN(twopi4*setl4(iset)),
     $             bl(ibl,1)%rz(ivert,2)
              ELSE
                WRITE(con_unit,*)
     $             bl(ibl,1)%rz(ivert,1)*COS(twopi*setlabel(iset)),
     $            -bl(ibl,1)%rz(ivert,1)*SIN(twopi*setlabel(iset)),
     $             bl(ibl,1)%rz(ivert,2)
              ENDIF
            ENDDO
          ELSE
            IF (nsets>1) THEN
              DO ivert=1,bl(ibl,1)%mvert
                IF (vtk_format=="BINARY") THEN
                  WRITE(con_unit) bl(ibl,1)%rz(ivert,:),setl4(iset)
                ELSE
                  WRITE(con_unit,*) bl(ibl,1)%rz(ivert,:),setlabel(iset)
                ENDIF
              ENDDO
            ELSE
              DO ivert=1,bl(ibl,1)%mvert
                IF (vtk_format=="BINARY") THEN
                  WRITE(con_unit) bl(ibl,1)%rz(ivert,:),0._r4
                ELSE
                  WRITE(con_unit,*) bl(ibl,1)%rz(ivert,:),0._r4
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write the element connections.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(field_label,'(a,2i9)') "CELLS ",
     $     ncell*MAX(1_i4,nsets-1_i4),ncdata
        WRITE(con_unit) linef//TRIM(field_label)//linef
      ELSE
        WRITE(con_unit,'(/,a,2i9)') "CELLS ",ncell*MAX(1_i4,nsets-1_i4),
     $     ncdata
      ENDIF
      IF (nsets==1) THEN
        DO ibl=1,nbl
          DO icell=1,bl(ibl,1)%mcell
            IF (vtk_format=="BINARY") THEN
              WRITE(con_unit) SIZE(bl(ibl,1)%vertex,2),
     $          bl(ibl,1)%vertex(icell,:)
            ELSE
              WRITE(con_unit,*) SIZE(bl(ibl,1)%vertex,2),
     $          bl(ibl,1)%vertex(icell,:)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO ibl=1,nbl
          DO iset=1,nsets-1
            DO icell=1,bl(ibl,1)%mcell
              IF (vtk_format=="BINARY") THEN
                WRITE(con_unit) 2_i4*SIZE(bl(ibl,1)%vertex,2),
     $            bl(ibl,1)%vertex(icell,:)+npt*(iset-1),
     $            bl(ibl,1)%vertex(icell,:)+npt*(iset)
              ELSE
                WRITE(con_unit,*) 2_i4*SIZE(bl(ibl,1)%vertex,2),
     $            bl(ibl,1)%vertex(icell,:)+npt*(iset-1),
     $            bl(ibl,1)%vertex(icell,:)+npt*(iset)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     write element types.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(field_label,'(a,i9)') "CELL_TYPES ",
     $    ncell*MAX(1_i4,nsets-1_i4)
        WRITE(con_unit) linef//TRIM(field_label)//linef
      ELSE
        WRITE(con_unit,'(/,a,i9)') "CELL_TYPES ",
     $    ncell*MAX(1_i4,nsets-1_i4)
      ENDIF
      DO ibl=1,nbl
        DO iset=1,MAX(1_i4,nsets-1_i4)
          DO icell=1,bl(ibl,1)%mcell
            IF (vtk_format=="BINARY") THEN
              WRITE(con_unit) bl(ibl,1)%el_type
            ELSE
              WRITE(con_unit,*) bl(ibl,1)%el_type
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     read data and save in memory.
c-----------------------------------------------------------------------
      DO iset=1,nsets
        DO ibl=1,nbl
          ALLOCATE(bl(ibl,iset)%data(bl(ibl,1)%mvert,nx))
          DO ix=1,nx
            READ(xdr_unit,IOSTAT=read_stat) bl(ibl,iset)%data(:,ix)
            IF (read_stat/=0) THEN
              WRITE(nim_wr,*) 'Error reading quantity ',
     $          TRIM(varnames(ix+2)),' in block ',ibl,' set ',iset
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     organize scalars and vectors separately.  check for vector fields
c     by whether the first character of consecutive names matches.
c-----------------------------------------------------------------------
      ix=1
      nfield=0
      DO
        inam=ix+2
        nfield=nfield+1
        veclist(nfield)=.false.
        IF (ix+2<=nx) THEN
          veclist(nfield)=vtk_vec_check(varnames(inam),varnames(inam+1),
     $                                  varnames(inam+2),field_label)
        ENDIF
        CALL vtk_label_suffix(varnames(inam),veclist(nfield),
     $                        field_label)

        IF (vtk_format=="BINARY") THEN
          field_label=TRIM(field_label)//" float"
        ELSE
          field_label=TRIM(field_label)//" double"
        ENDIF

        IF (veclist(nfield)) THEN
          header(nfield)="VECTORS "//TRIM(field_label)
          narr=2
        ELSE
          header(nfield)="SCALARS "//TRIM(field_label)//" 1"
          narr=0
        ENDIF
        ix=ix+narr+1
        IF (ix>nx) EXIT
      ENDDO
c-----------------------------------------------------------------------
c     write the data in point format.  make two passes writing the
c     vectors first and then the scalars.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(field_label,'(a,i9)') "POINT_DATA ",npt*nsets
        WRITE(con_unit) linef//TRIM(field_label)
      ELSE
        WRITE(con_unit,'(/,a,i9)') "POINT_DATA ",npt*nsets
      ENDIF

      ix=1
      DO ifield=1,nfield
        IF (veclist(ifield)) THEN
          IF (vtk_format=="BINARY") THEN
            WRITE(con_unit) linef//TRIM(header(ifield))//linef
          ELSE
            WRITE(con_unit,'(/,a)') TRIM(header(ifield))
          ENDIF
          narr=2
          DO iset=1,nsets
            DO ibl=1,nbl
              DO ivert=1,bl(ibl,1)%mvert
                IF (vtk_format=="BINARY") THEN
                  WRITE(con_unit) bl(ibl,iset)%data(ivert,ix:ix+narr)
                ELSE
                  WRITE(con_unit,*) bl(ibl,iset)%data(ivert,ix:ix+narr)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          narr=0
        ENDIF
        ix=ix+narr+1
      ENDDO

      ix=1
      DO ifield=1,nfield
        IF (.NOT.veclist(ifield)) THEN
          IF (vtk_format=="BINARY") THEN
            WRITE(con_unit) linef//TRIM(header(ifield))
            WRITE(con_unit) linef//"LOOKUP_TABLE default"//linef
          ELSE
            WRITE(con_unit,'(/,a)') TRIM(header(ifield))
            WRITE(con_unit,'(a)') "LOOKUP_TABLE default"
          ENDIF
          narr=0
          DO iset=1,nsets
            DO ibl=1,nbl
              DO ivert=1,bl(ibl,1)%mvert
                IF (vtk_format=="BINARY") THEN
                  WRITE(con_unit) bl(ibl,iset)%data(ivert,ix:ix+narr)
                ELSE
                  WRITE(con_unit,*) bl(ibl,iset)%data(ivert,ix:ix+narr)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          narr=2
        ENDIF
        ix=ix+narr+1
      ENDDO
c-----------------------------------------------------------------------
c     close files
c-----------------------------------------------------------------------
      CLOSE(con_unit)
      CALL close_bin(xdr_unit,TRIM(binfile))
c-----------------------------------------------------------------------
c     deallocate arrays and terminate.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        DEALLOCATE(bl(ibl,1)%vertex)
        DEALLOCATE(bl(ibl,1)%rz)
      ENDDO
      DO iset=1,nsets
        DO ibl=1,nbl
          DEALLOCATE(bl(ibl,iset)%data)
        ENDDO
      ENDDO
      DEALLOCATE(bl)
      DEALLOCATE(varnames)

      RETURN
      END SUBROUTINE vtk_contour_conv
c-----------------------------------------------------------------------
c     subprogram 3. vtk_vec_check
c     this logical function parses the variable names to determine
c     if a set of three forms a vector field.
c-----------------------------------------------------------------------
      LOGICAL FUNCTION vtk_vec_check(var1,var2,var3,vecroot)

      CHARACTER(*), INTENT(IN) :: var1,var2,var3
      CHARACTER(*), INTENT(OUT) :: vecroot

      INTEGER(i4) :: ichar,jchar,kchar,numch,icoord
      INTEGER(i4), PARAMETER :: ncoord=4_i4
      CHARACTER(1), DIMENSION(3,ncoord), PARAMETER ::
     $  ccoord=RESHAPE((/'R','Z','P','r','z','p','X','Y','Z',
     $                   'x','y','z'/),(/3_i4,ncoord/))
c-----------------------------------------------------------------------
c     sift through the characters of var1-3 to check for consecutive
c     coordinate labels.
c-----------------------------------------------------------------------
      vtk_vec_check=.false.
      numch=LEN_TRIM(var1)
      ich: DO ichar=1,numch
        DO icoord=1,ncoord
          IF (var1(ichar:ichar)==ccoord(1,icoord)) EXIT ich
        ENDDO
      ENDDO ich
      IF (ichar>numch) RETURN

      DO jchar=1,numch
        IF (var2(jchar:jchar)==ccoord(2,icoord)) EXIT
      ENDDO
      IF (jchar>numch) RETURN

      DO kchar=1,numch
        IF (var3(kchar:kchar)==ccoord(3,icoord)) EXIT
      ENDDO
      IF (kchar>numch) RETURN

      vtk_vec_check=.true.
      vecroot=TRIM(var1(1:ichar-1))

      RETURN
      END FUNCTION vtk_vec_check
c-----------------------------------------------------------------------
c     subprogram 4. vtk_label_suffix
c     this subroutine adds a numerical suffix to the field label if
c     the variable name indicates a repeated label.  vectors and
c     scalars are handled separately, because field_label already
c     contains the root name for vectors.
c-----------------------------------------------------------------------
      SUBROUTINE vtk_label_suffix(varname,vec,field_label)

      LOGICAL, INTENT(IN) :: vec
      CHARACTER(*), INTENT(IN) :: varname
      CHARACTER(*), INTENT(INOUT) :: field_label

      INTEGER(i4) :: ichar,numch
      LOGICAL :: addnum
c-----------------------------------------------------------------------
c     sift through the characters of varname to check for the 
c     occurrence of "#."
c-----------------------------------------------------------------------
      numch=LEN_TRIM(varname)
      addnum=.false.
      DO ichar=1,numch
        IF (varname(ichar:ichar)=='#') THEN
          addnum=.true.
          EXIT
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     collect the numerical suffix and modify or create the field_label.
c-----------------------------------------------------------------------
      IF (addnum) THEN
        IF (vec) THEN
          field_label=TRIM(field_label)//"_"//
     $                TRIM(ADJUSTL(varname(ichar+1:)))
        ELSE
          field_label=TRIM(varname(1:ichar-1))//"_"//
     $                TRIM(ADJUSTL(varname(ichar+1:)))
        ENDIF
      ELSE IF (.NOT.vec) THEN
        field_label=TRIM(varname)
      ENDIF

      RETURN
      END SUBROUTINE vtk_label_suffix
c-----------------------------------------------------------------------
c     subprogram 4. vtk_eldata_write
c     this subroutine writes fields that are discontinuous across
c     element borders in vtk's unstructed format.  all of the required
c     information is passed in the r4eldata_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vtk_eldata_write(eldata,filename,unitno)
      USE eldata_type_mod

      TYPE(r4eldata_type), DIMENSION(:), INTENT(IN) :: eldata
      CHARACTER(*), INTENT(IN) :: filename
      INTEGER(i4), INTENT(IN) :: unitno

      INTEGER(i4) :: npt,ipt,nbl,ibl,iel,nel,isub,ivar,ncdata,numch,ich
      CHARACTER(128) :: varname,chstring
c-----------------------------------------------------------------------
c     open the vtk file, and write the header information.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        OPEN(UNIT=unitno,FILE=TRIM(filename),STATUS='UNKNOWN',
     $       ACCESS="STREAM")
        WRITE(unitno) "# vtk DataFile Version 2.0"
        WRITE(unitno) linef//"Unstructured "//
     $    TRIM(eldata(1)%coordlabels)//" Grid and"
        DO ivar=1,eldata(1)%nfields
          WRITE(unitno) " "//TRIM(eldata(1)%datalabels(ivar))
        ENDDO
        WRITE(unitno) " Data"
        WRITE(unitno) linef//"BINARY"
        WRITE(unitno) linef//"DATASET UNSTRUCTURED_GRID"
      ELSE
        OPEN(unitno,FILE=TRIM(filename),STATUS='UNKNOWN')
        WRITE(unitno,'(a)') "# vtk DataFile Version 2.0"
        WRITE(unitno,'(3a)',ADVANCE='NO') "Unstructured ",
     $    TRIM(eldata(1)%coordlabels)," Grid and"
        DO ivar=1,eldata(1)%nfields
          WRITE(unitno,'(2a)',ADVANCE='NO') " ",
     $      TRIM(eldata(1)%datalabels(ivar))
        ENDDO
        WRITE(unitno,'(a)') " Data"
        WRITE(unitno,'(a,/,/,a)') "ASCII","DATASET UNSTRUCTURED_GRID"
      ENDIF
c-----------------------------------------------------------------------
c     write the coordinates of each data point.  also determine the
c     amount of vtk cell data.
c-----------------------------------------------------------------------
      nbl=SIZE(eldata)
      npt=0
      nel=0
      ncdata=0
      DO ibl=1,nbl
        npt=npt+SIZE(eldata(ibl)%coords(1,:,:))
      ENDDO
      IF (vtk_format=="BINARY") THEN
        WRITE(chstring,'(a,i9,a)') "POINTS ",npt," float"
        WRITE(unitno) linef//TRIM(chstring)//linef
      ELSE
        WRITE(unitno,'(/,a,i9,a)') "POINTS ",npt," double"
      ENDIF
      DO ibl=1,nbl
        DO iel=1,eldata(ibl)%nele
          DO ipt=1,eldata(ibl)%nptperel
            IF (vtk_format=="BINARY") THEN
              WRITE(unitno) eldata(ibl)%coords(:,ipt,iel),0._r4
            ELSE
              WRITE(unitno,*) eldata(ibl)%coords(:,ipt,iel),0._r4
            ENDIF
          ENDDO
          nel=nel+SIZE(eldata(ibl)%connect,2)
          SELECT CASE(eldata(ibl)%eltype)
          CASE ("quad")  !  sub-element representation for vtk
            ncdata=ncdata+5*SIZE(eldata(ibl)%connect,2)
          CASE ("tri")
            ncdata=ncdata+4*SIZE(eldata(ibl)%connect,2)
          END SELECT
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write element connections.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(chstring,'(a,2i9)') "CELLS ",nel,ncdata
        WRITE(unitno) linef//TRIM(chstring)//linef
      ELSE
        WRITE(unitno,'(/,a,2i9)') "CELLS ",nel,ncdata
      ENDIF
      DO ibl=1,nbl
        DO iel=1,eldata(ibl)%nele
          DO isub=1,SIZE(eldata(ibl)%connect,2)
            IF (vtk_format=="BINARY") THEN
              WRITE(unitno) SIZE(eldata(ibl)%connect,1),
     $          eldata(ibl)%connect(:,isub,iel)
            ELSE
              WRITE(unitno,*) SIZE(eldata(ibl)%connect,1),
     $          eldata(ibl)%connect(:,isub,iel)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write element types.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(chstring,'(a,2i9)') "CELL_TYPES ",nel
        WRITE(unitno) linef//TRIM(chstring)//linef
      ELSE
        WRITE(unitno,'(/,a,2i9)') "CELL_TYPES ",nel
      ENDIF
      DO ibl=1,nbl
        SELECT CASE(eldata(ibl)%eltype)
        CASE ("quad")
          DO iel=1,eldata(ibl)%nele
            DO isub=1,SIZE(eldata(ibl)%connect,2)
              IF (vtk_format=="BINARY") THEN
                WRITE(unitno) 9_i4
              ELSE
                WRITE(unitno,*) 9_i4
              ENDIF
            ENDDO
          ENDDO
        CASE ("tri")
          DO iel=1,eldata(ibl)%nele
            DO isub=1,SIZE(eldata(ibl)%connect,2)
              IF (vtk_format=="BINARY") THEN
                WRITE(unitno) 5_i4
              ELSE
                WRITE(unitno,*) 5_i4
              ENDIF
            ENDDO
          ENDDO
        END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     write data as scalars, as contours are best for showing
c     discontinuity.
c-----------------------------------------------------------------------
      IF (vtk_format=="BINARY") THEN
        WRITE(chstring,'(a,i9)') "POINT_DATA ",npt
        WRITE(unitno) linef//TRIM(chstring)
      ELSE
        WRITE(unitno,'(/,a,i9)') "POINT_DATA ",npt
      ENDIF

      DO ivar=1,eldata(1)%nfields
        varname=eldata(1)%datalabels(ivar)
        numch=LEN_TRIM(varname)
        IF (vtk_format=="BINARY") THEN
          WRITE(unitno) linef//"SCALARS "
          DO ich=1,numch
            IF (varname(ich:ich)==" ") THEN
              WRITE(unitno) "_"
            ELSE
              WRITE(unitno) varname(ich:ich)
            ENDIF
          ENDDO
          WRITE(unitno) " float 1"
          WRITE(unitno) linef//"LOOKUP_TABLE default"//linef
          DO ibl=1,nbl
            DO iel=1,eldata(ibl)%nele
              DO ipt=1,eldata(ibl)%nptperel
                WRITE(unitno) eldata(ibl)%data(ivar,ipt,iel)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          WRITE(unitno,'(/,a)',ADVANCE='NO') "SCALARS "
          DO ich=1,numch
            IF (varname(ich:ich)==" ") THEN
              WRITE(unitno,'(a)',ADVANCE='NO') "_"
            ELSE
              WRITE(unitno,'(a)',ADVANCE='NO') varname(ich:ich)
            ENDIF
          ENDDO
          WRITE(unitno,'(a)') " double 1"
          WRITE(unitno,'(a)') "LOOKUP_TABLE default"
          DO ibl=1,nbl
            DO iel=1,eldata(ibl)%nele
              DO ipt=1,eldata(ibl)%nptperel
                WRITE(unitno,*) eldata(ibl)%data(ivar,ipt,iel)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CLOSE(unitno)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vtk_eldata_write
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE vtk_mod
