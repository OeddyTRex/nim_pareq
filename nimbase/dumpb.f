c-----------------------------------------------------------------------
c     file dumpb.f
c     this file contains a dump reader and a dump writer that are
c     external programs.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  dumpb_read.
c     2.  dumpb_write.
c-----------------------------------------------------------------------
c     subprogram 1. dumpb_read.
c     reads a dump file without the consistency checks.
c-----------------------------------------------------------------------
      SUBROUTINE dumpb_read(nmodes,keff,dump_time,dump_step)
      USE local
      USE fields
      USE dump
      IMPLICIT NONE

      INTEGER(i4), INTENT(OUT) :: nmodes
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff

      INTEGER(i4) :: ib,nb,i
      REAL(r8) :: iread
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Dump file '//TRIM(dump_file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(dump_file),'OLD','REWIND',64_i4)
      READ(rstrt_unit) dump_time
      READ(rstrt_unit) iread
      dump_step=NINT(iread)
      READ(rstrt_unit) iread
      nbl=NINT(iread)
      READ(rstrt_unit) iread
      nrbl=NINT(iread)
      READ(rstrt_unit) iread
      poly_degree=NINT(iread)
      READ(rstrt_unit) iread
      nmodes=NINT(iread)
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
      END SUBROUTINE dumpb_read
c-----------------------------------------------------------------------
c     subprogram 2. dumpb_write.
c     writes a dump file from the data structures in the argument list.
c-----------------------------------------------------------------------
      SUBROUTINE dumpb_write(nmodes,keff,dump_time,dump_step,rbl,tbl)
      USE local
      USE input
      USE fields
      USE dump
      IMPLICIT NONE
      
      INTEGER(i4), INTENT(IN) :: nmodes,dump_step
      REAL(r8), INTENT(IN) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff
      TYPE(rblock_type), DIMENSION(:), POINTER :: rbl
      TYPE(tblock_type), DIMENSION(:), POINTER :: tbl

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
      dump_file=TRIM(dump_name)//"."//TRIM(step_name)
c-----------------------------------------------------------------------
c     open dump file.
c-----------------------------------------------------------------------
      CALL open_bin(dump_unit,TRIM(dump_file),'REPLACE','ASIS',64_i4)
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
c     write rblock data.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         CALL dump_write_rblock(rbl(ib))
      ENDDO
c-----------------------------------------------------------------------
c     write tblock data.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         CALL dump_write_tblock(tbl(ib))
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(dump_unit,TRIM(dump_file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END  SUBROUTINE dumpb_write
