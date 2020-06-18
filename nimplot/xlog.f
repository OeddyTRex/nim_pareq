c-----------------------------------------------------------------------
c     file xlog.f:  reads output written for xdraw and takes the log
c     of the data.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM xlog
      USE local
      IMPLICIT NONE

      CHARACTER(64) :: rfile,wfile,format
      INTEGER(i4) :: ix,it=-1,iq,ni,nx,nxlim=2000,
     $               irec,nrec,nreclim=2000,num_chars,nwrite,nwr=0,ndbl
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: comp_write
      REAL(r4), DIMENSION(:), ALLOCATABLE :: data,out_data
      REAL(r4), DIMENSION(:,:), ALLOCATABLE :: data_old
      INTEGER :: read_stat,old_stat=0
      CHARACTER(1) :: anss='n',ans10='n'
      CHARACTER(32) :: ctmp
      REAL(r8) :: ran
      LOGICAL :: write_blank=.true.,write_data=.true.
c-----------------------------------------------------------------------
c     interactive dialogue.
c-----------------------------------------------------------------------
      WRITE(*,*) 'Enter the file name of the existing xdraw binary'
      WRITE(*,*) 'and the name of the new output file.'
      READ(*,*) rfile,wfile
      IF (rfile(1:7)=='nimhist'.OR.rfile(1:6)=='energy') THEN
        ni=4
        it=2
      ELSE IF (rfile(1:9)=='discharge') THEN
        ni=2
        it=2
      ELSE
        WRITE(*,*) 'Enter the number of independent variables.'
        READ(*,*) ni
        WRITE(*,*) 'Enter the index for time.'
        READ(*,*) it
      ENDIF
      WRITE(*,*)
     $    'Do you want slopes instead of values. [Type y or n (ret=n).]'
      READ(*,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $     SIZE=num_chars) ctmp
      IF (num_chars>0) anss=ctmp
      WRITE(*,*) 'Do you want log base 10? [Type y or n (ret=n).]'
      READ(*,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $     SIZE=num_chars) ctmp
      IF (num_chars>0) ans10=ctmp
      IF (rfile(1:7)=='nimhist'.OR.rfile(1:6)=='energy') THEN
        WRITE(*,*)
     $    'Enter the # of Fourier components to be written from the '
     $    //'data. [Hit ret for all comps.]'
        READ(*,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $       SIZE=num_chars) ctmp
        IF (num_chars>0) THEN
          READ(ctmp,*) nwr
          ALLOCATE(comp_write(nwr))
          WRITE(*,*) 'Enter the list of Fourier components to write.'
          READ(*,*) comp_write
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     determine record size.
c-----------------------------------------------------------------------
      DO nx=1,nxlim+1
        ALLOCATE(data(nx))
        CALL open_bin(xdr_unit,TRIM(rfile),'OLD','REWIND',32_i4)
        READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
        CALL close_bin(xdr_unit,TRIM(rfile))
        DEALLOCATE(data)
        IF (read_stat>0) THEN
          EXIT     !     one beyond record length
        ELSE IF (read_stat<0) THEN
          WRITE(*,*) 'Prematurely hit the end of file.'
          STOP
        ENDIF
      ENDDO

      IF (nx>nxlim) THEN
        WRITE(*,*) 'Record size too large.'
        STOP
      ELSE IF (nx==1) THEN
        WRITE(*,*) 'File is empty???'
        STOP
      ENDIF
      nx=nx-1

      WRITE(*,*) 'Record length: ',nx
      WRITE(format,'(a,i2,a)') '(',nx,'es10.3)'
      ALLOCATE(data(nx),out_data(nx))
c-----------------------------------------------------------------------
c     determine number of records per data block.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(rfile),'OLD','REWIND',32_i4)
      DO nrec=1,nreclim+1
        READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
        IF (read_stat/=0) EXIT
      ENDDO
      CALL close_bin(xdr_unit,TRIM(rfile))
      IF (nrec>nreclim) THEN
        WRITE(*,*) 'Too many records.'
        STOP
      ELSE IF (nrec==1) THEN
        WRITE(*,*) 'File is empty???'
        STOP
      ENDIF
      nrec=nrec-1
      IF (anss=='y') THEN
        ALLOCATE(data_old(nx,nrec))
        data_old=-1.e20
        write_blank=.false.
      ENDIF

      WRITE(*,*) 'Number of records: ',nrec
c-----------------------------------------------------------------------
c     open files.
c-----------------------------------------------------------------------
      CALL open_bin(xdr_unit,TRIM(rfile),'OLD','REWIND',32_i4)
      CALL open_bin(out_unit,TRIM(wfile),'UNKNOWN','REWIND',32_i4)
c     OPEN(UNIT=out_unit,FILE=TRIM(wfile),STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     read each record, write log of dependent variables, and skip 
c     blank records.
c-----------------------------------------------------------------------
      ndbl=0
      DO 
        READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
        IF (read_stat==0) THEN
          ndbl=ndbl+1
          irec=1
          CALL process_data
          DO irec=2,nrec
            READ(xdr_unit,IOSTAT=read_stat) data(1:nx)
            IF (read_stat==0) THEN
              CALL process_data
            ELSE
              WRITE(*,*) 'Record block ',ndbl,' incomplete.'
              STOP
            ENDIF
          ENDDO
        ELSE
          EXIT
        ENDIF
        READ(xdr_unit,IOSTAT=read_stat)
        IF (read_stat/=0) EXIT
        IF (write_blank) WRITE(out_unit)
      ENDDO

      CALL close_bin(xdr_unit,TRIM(rfile))
      CALL close_bin(out_unit,TRIM(wfile))
c     CLOSE(out_unit)
      
      STOP

c     Internal function for data processing.

      CONTAINS

        SUBROUTINE process_data

          DO ix=ni+1,nx
            IF (ABS(data(ix))<1.e-20) THEN
              CALL RANDOM_NUMBER(ran)
              data(ix)=1.e-20+1.e-23*ran
            ENDIF
          ENDDO
          IF (anss=='y') THEN
            IF (irec>nrec) THEN
              WRITE(*,*) 'Record index error.'
              STOP
            ENDIF
            IF (data_old(1,irec)>-1.e20) THEN
              write_data=.true.
              write_blank=.true.
              IF (ans10=='y') THEN
                out_data=(/data(1:ni),
     $            (LOG10(ABS(data(ni+1:nx)))
     $            -LOG10(ABS(data_old(ni+1:nx,irec))))
     $            /(data(it)-data_old(it,irec))/)
              ELSE
                out_data=(/data(1:ni),
     $            (LOG(ABS(data(ni+1:nx)))
     $            -LOG(ABS(data_old(ni+1:nx,irec))))
     $            /(data(it)-data_old(it,irec))/)
              ENDIF
            ELSE
              write_data=.false.
            ENDIF
            data_old(:,irec)=data
          ELSE
            IF (ans10=='y') THEN
              out_data=(/data(1:ni),LOG10(ABS(data(ni+1:nx)))/)
            ELSE
              out_data=(/data(1:ni),LOG(ABS(data(ni+1:nx)))/)
            ENDIF
          ENDIF
          IF (write_data) THEN
            IF (nwr==0) THEN
              WRITE(out_unit) REAL(out_data,4)
            ELSE
              IF (MINVAL(ABS(comp_write-irec+1))==0)
     $          WRITE(out_unit) REAL(out_data,4)
            ENDIF
          ENDIF

        END SUBROUTINE process_data

      END PROGRAM xlog
