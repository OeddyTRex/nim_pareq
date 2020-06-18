c-----------------------------------------------------------------------
c     Randomly change the phase of the Fourier components in a
c     dump file and rewrite the file.
c-----------------------------------------------------------------------
      PROGRAM mode_scramble
      USE local
      USE dump
      USE input
      USE fields
      IMPLICIT NONE

      CHARACTER(128) :: old_file,new_file,op,ctmp
      INTEGER(i4) :: nmodes,dump_step,ib,im,numch,read_stat
      INTEGER :: nseed
      INTEGER, DIMENSION(:), ALLOCATABLE :: newseed
      REAL(r8) :: dump_time,rphase
      REAL(r8), DIMENSION(:), POINTER :: keff
      COMPLEX(r8) :: cfactor

c-----------------------------------------------------------------------
c     interface block for dump_read_ms.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dump_read_ms(nmodes,keff,dump_time,dump_step,file)
        USE local
        USE dump
        USE fields
        USE input
        IMPLICIT NONE
 
        INTEGER(i4), INTENT(OUT) :: nmodes
        INTEGER(i4), INTENT(OUT) :: dump_step
        REAL(r8), INTENT(OUT) :: dump_time
        REAL(r8), DIMENSION(:), POINTER :: keff
        CHARACTER(*), INTENT(IN) :: file
        END SUBROUTINE dump_read_ms
      END INTERFACE
c-----------------------------------------------------------------------
c     get files names and read old dump file.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) "Enter the names of the original dump file."
      READ(nim_rd,*) old_file
      new_file="mix_"//TRIM(dump_name)

      CALL dump_read_ms(nmodes,keff,dump_time,dump_step,old_file)
c-----------------------------------------------------------------------
c     select whether to set a random seed, skip randomization, or 
c     set a random seed.
c-----------------------------------------------------------------------
      op='skip'
      WRITE(*,*)
      WRITE(*,*) "Enter ","'","get","'"," to get a random seed,"
      WRITE(*,*) "'","set","'"," to set a random seed, or hit"
      WRITE(*,*) "return to avoid random mixing for checking."
      READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $     SIZE=numch) ctmp
      IF (numch>0) op=ctmp
      IF (op/='get'.AND.op/='set') op='skip'
c-----------------------------------------------------------------------
c     set a random phase seed.
c-----------------------------------------------------------------------
      IF (op=='get') THEN
        CALL RANDOM_SEED
        CALL RANDOM_SEED(size=nseed)
        ALLOCATE(newseed(nseed))
        CALL RANDOM_SEED(get=newseed)
        WRITE(*,*) "Random seed is ",newseed
      ELSE IF (op=='set') THEN
        CALL RANDOM_SEED(size=nseed)
        WRITE(*,*) "Enter the new ",nseed," integer seed number(s)."
        ALLOCATE(newseed(nseed))
        READ(*,*) newseed
        CALL RANDOM_SEED(put=newseed)
      ENDIF
c-----------------------------------------------------------------------
c     generate a new random phase for each Fourier component (phir), and
c     multiply the complex coefficients by exp(-i*k*phir).  This moves
c     each sine wave of the original data along phir radians in toroidal
c     angle.
c-----------------------------------------------------------------------
      rphase=0._r8
      mode_loop: DO im=1,nmodes
        IF (op/='skip') CALL RANDOM_NUMBER(rphase)
        rphase=rphase*twopi
        cfactor=COS(keff(im)*rphase)+(0,1)*SIN(keff(im)*rphase)
c-TMP
        WRITE(*,*) im,rphase,cfactor

        DO ib=1,nrbl
          rb(ib)%be%fs(:,:,:,im)=rb(ib)%be%fs(:,:,:,im)*cfactor
          rb(ib)%ve%fs(:,:,:,im)=rb(ib)%ve%fs(:,:,:,im)*cfactor
          rb(ib)%pres%fs(:,:,:,im)=rb(ib)%pres%fs(:,:,:,im)*cfactor
          rb(ib)%prese%fs(:,:,:,im)=rb(ib)%prese%fs(:,:,:,im)*cfactor
          rb(ib)%nd%fs(:,:,:,im)=rb(ib)%nd%fs(:,:,:,im)*cfactor
          rb(ib)%conc%fs(:,:,:,im)=rb(ib)%conc%fs(:,:,:,im)*cfactor
          rb(ib)%tele%fs(:,:,:,im)=rb(ib)%tele%fs(:,:,:,im)*cfactor
          rb(ib)%tion%fs(:,:,:,im)=rb(ib)%tion%fs(:,:,:,im)*cfactor
          IF (poly_degree>1) THEN
            rb(ib)%be%fsh(:,:,:,:,im)=rb(ib)%be%fsh(:,:,:,:,im)*cfactor
            rb(ib)%ve%fsh(:,:,:,:,im)=rb(ib)%ve%fsh(:,:,:,:,im)*cfactor
            rb(ib)%pres%fsh(:,:,:,:,im)=
     $        rb(ib)%pres%fsh(:,:,:,:,im)*cfactor
            rb(ib)%prese%fsh(:,:,:,:,im)=
     $        rb(ib)%prese%fsh(:,:,:,:,im)*cfactor
            rb(ib)%nd%fsh(:,:,:,:,im)=
     $        rb(ib)%nd%fsh(:,:,:,:,im)*cfactor
            rb(ib)%conc%fsh(:,:,:,:,im)=
     $        rb(ib)%conc%fsh(:,:,:,:,im)*cfactor
            rb(ib)%tele%fsh(:,:,:,:,im)=
     $        rb(ib)%tele%fsh(:,:,:,:,im)*cfactor
            rb(ib)%tion%fsh(:,:,:,:,im)=
     $        rb(ib)%tion%fsh(:,:,:,:,im)*cfactor

            rb(ib)%be%fsv(:,:,:,:,im)=rb(ib)%be%fsv(:,:,:,:,im)*cfactor
            rb(ib)%ve%fsv(:,:,:,:,im)=rb(ib)%ve%fsv(:,:,:,:,im)*cfactor
            rb(ib)%pres%fsv(:,:,:,:,im)=
     $        rb(ib)%pres%fsv(:,:,:,:,im)*cfactor
            rb(ib)%prese%fsv(:,:,:,:,im)=
     $        rb(ib)%prese%fsv(:,:,:,:,im)*cfactor
            rb(ib)%nd%fsv(:,:,:,:,im)=
     $        rb(ib)%nd%fsv(:,:,:,:,im)*cfactor
            rb(ib)%conc%fsv(:,:,:,:,im)=
     $        rb(ib)%conc%fsv(:,:,:,:,im)*cfactor
            rb(ib)%tele%fsv(:,:,:,:,im)=
     $        rb(ib)%tele%fsv(:,:,:,:,im)*cfactor
            rb(ib)%tion%fsv(:,:,:,:,im)=
     $        rb(ib)%tion%fsv(:,:,:,:,im)*cfactor

            rb(ib)%be%fsi(:,:,:,:,im)=rb(ib)%be%fsi(:,:,:,:,im)*cfactor
            rb(ib)%ve%fsi(:,:,:,:,im)=rb(ib)%ve%fsi(:,:,:,:,im)*cfactor
            rb(ib)%pres%fsi(:,:,:,:,im)=
     $        rb(ib)%pres%fsi(:,:,:,:,im)*cfactor
            rb(ib)%prese%fsi(:,:,:,:,im)=
     $        rb(ib)%prese%fsi(:,:,:,:,im)*cfactor
            rb(ib)%nd%fsi(:,:,:,:,im)=
     $        rb(ib)%nd%fsi(:,:,:,:,im)*cfactor
            rb(ib)%conc%fsi(:,:,:,:,im)=
     $        rb(ib)%conc%fsi(:,:,:,:,im)*cfactor
            rb(ib)%tele%fsi(:,:,:,:,im)=
     $        rb(ib)%tele%fsi(:,:,:,:,im)*cfactor
            rb(ib)%tion%fsi(:,:,:,:,im)=
     $        rb(ib)%tion%fsi(:,:,:,:,im)*cfactor
          ENDIF
        ENDDO

        DO ib=nrbl+1,nbl
          tb(ib)%be%fs(:,:,:,im)=tb(ib)%be%fs(:,:,:,im)*cfactor
          tb(ib)%ve%fs(:,:,:,im)=tb(ib)%ve%fs(:,:,:,im)*cfactor
          tb(ib)%pres%fs(:,:,:,im)=tb(ib)%pres%fs(:,:,:,im)*cfactor
          tb(ib)%prese%fs(:,:,:,im)=tb(ib)%prese%fs(:,:,:,im)*cfactor
          tb(ib)%nd%fs(:,:,:,im)=tb(ib)%nd%fs(:,:,:,im)*cfactor
          tb(ib)%conc%fs(:,:,:,im)=tb(ib)%conc%fs(:,:,:,im)*cfactor
          tb(ib)%tele%fs(:,:,:,im)=tb(ib)%tele%fs(:,:,:,im)*cfactor
          tb(ib)%tion%fs(:,:,:,im)=tb(ib)%tion%fs(:,:,:,im)*cfactor
        ENDDO
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     write new dump file.
c-----------------------------------------------------------------------
      dump_name=new_file
      CALL dump_write(nmodes,keff,dump_time,dump_step)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')
      END PROGRAM mode_scramble

c-----------------------------------------------------------------------
c     subprogram 1. dump_read_ms.
c     reads dump file without normal checks.
c-----------------------------------------------------------------------
      SUBROUTINE dump_read_ms(nmodes,keff,dump_time,dump_step,file)
      USE local
      USE dump
      USE fields
      USE input
      IMPLICIT NONE

      INTEGER(i4), INTENT(OUT) :: nmodes
      INTEGER(i4), INTENT(OUT) :: dump_step
      REAL(r8), INTENT(OUT) :: dump_time
      REAL(r8), DIMENSION(:), POINTER :: keff
      CHARACTER(*), INTENT(IN) :: file

      INTEGER(i4) :: ib,nb,dump_code = 0_i4,nmodes_dump,i
      REAL(r8) :: iread
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     64 bit reals then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Dump file '//TRIM(file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(file),'OLD','REWIND',64_i4)
      READ(rstrt_unit) dump_time
      READ(rstrt_unit) iread
      dump_step=NINT(iread)
      READ(rstrt_unit) iread
      nbl=NINT(iread)
      READ(rstrt_unit) iread
      nrbl=NINT(iread)
c-----------------------------------------------------------------------
c     check for basis consistency.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      poly_degree=NINT(iread)
c-----------------------------------------------------------------------
c     check for fourier rep. consistency, then read wavenumber array.
c-----------------------------------------------------------------------
      READ(rstrt_unit) iread
      nmodes_dump=NINT(iread)
      nmodes=nmodes_dump
      ALLOCATE(keff(nmodes))
      READ(rstrt_unit) keff
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
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
      CALL close_bin(rstrt_unit,TRIM(file))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dump_read_ms
