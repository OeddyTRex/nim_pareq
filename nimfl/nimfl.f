      PROGRAM nimfl
      USE local
      USE io
      USE input
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE start_positions
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ierror
      INTEGER(i4) :: nmodesc
      REAL(r8) :: iread
      INTEGER(i4), DIMENSION(:,:,:), ALLOCATABLE :: ijstart
      REAL(r8), DIMENSION(:), ALLOCATABLE :: keffc
      REAL(r8) :: dump_timec
      LOGICAL :: file_stat
      INTEGER(i4) :: ibl,iv
      LOGICAL:: found,failure,any_cross
      TYPE(location_type) :: p0
      REAL(r8) :: x_start,y_start,d_ind,rmax
      TYPE(cell_type), POINTER :: item,item_pre

c-----------------------------------------------------------------------
c     parallel initialization
c-----------------------------------------------------------------------
      CALL mpi_init(ierror)
      CALL mpi_comm_rank(mpi_comm_world,node,ierror)
      CALL mpi_comm_size(mpi_comm_world,nprocs,ierror)
c-----------------------------------------------------------------------
c     read input file.
c-----------------------------------------------------------------------
      IF (node==0) CALL read_namelist("nimrod.in",.FALSE.)
      IF (nprocs>1) CALL broadcast_input
      CALL read_namelist0("nimfl.in",.FALSE.)
c-----------------------------------------------------------------------
c     read restart dump.
c-----------------------------------------------------------------------
c     nphi=2**lphi
c     nmodes_total=nphi/3+1
c     nmodes=nmodes_total
c     ALLOCATE(keff(nmodes))
c     CALL dump_read(nmodes,keff,t,istep)
c-----------------------------------------------------------------------
c     open dump file and read global data.  integers are read into
c     a 64 bit variable then converted upon copying.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(dump_file),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('Dump file '//TRIM(dump_file)//' does not exist.')
      CALL open_bin(rstrt_unit,TRIM(dump_file),'OLD','REWIND',64_i4)
      READ(rstrt_unit)dump_timec
      READ(rstrt_unit)iread
      READ(rstrt_unit)iread
      nblc=NINT(iread)
      READ(rstrt_unit)iread
      nrblc=NINT(iread)
      READ(rstrt_unit) iread
      poly_degreec=NINT(iread)
      READ(rstrt_unit)iread
      nmodesc=NINT(iread)
      nmodes=nmodesc				! Needed for get_bfield
      ALLOCATE(keffc(nmodesc))
      READ(rstrt_unit)keffc
c-----------------------------------------------------------------------
c     read seam data.
c-----------------------------------------------------------------------
      ALLOCATE(seamc(nblc))
      CALL dump_read_seam(seam0c)
      DO ibl=1,nblc
         CALL dump_read_seam(seamc(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     read rblock data.
c-----------------------------------------------------------------------
      ALLOCATE(rbc(nrblc))
      DO ibl=1,nrblc
         CALL dump_read_rblock(rbc(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     read tblock data.
c-----------------------------------------------------------------------
      ALLOCATE(tbc(nrblc+1:nblc))
      DO ibl=nrblc+1,nblc
         CALL dump_read_tblock(tbc(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     close dump file.
c-----------------------------------------------------------------------
      CALL close_bin(rstrt_unit,TRIM(dump_file))
      WRITE(nim_wr,*)"Dump file read ",TRIM(dump_file)
c-----------------------------------------------------------------------
c     Preliminary bounds and periodicity determination for rectangular grids.
c
c     revert to using input and normal correlation of x::r, y::z
c-----------------------------------------------------------------------
c     IF(gridshape == 'rect')THEN
c       ymin=rbc(1)%rz%fs(1,0,0)
c       xmin=rbc(1)%rz%fs(2,0,0)
c       ymax=rbc(nrblc)%rz%fs(1,rbc(nrblc)%mx,rbc(nrblc)%my)
c       xmax=rbc(nrblc)%rz%fs(2,rbc(nrblc)%my,rbc(nrblc)%my)
c       IF(seamc(1)%vertex(1)%ptr(1,1) > 0 )THEN
c         iv=2*(rbc(1)%mx+rbc(1)%my)-1
c         IF(seamc(1)%vertex(iv)%ptr(1,1) > 0 )THEN
c            periodicity='both'
c         ELSE
c            periodicity='y-dir'
c         ENDIF
c       ENDIF
c     ENDIF
c-----------------------------------------------------------------------
c     Build cell information
c-----------------------------------------------------------------------
      CALL make_cell
      CALL make_search_map
      item => start
      DO WHILE(ASSOCIATED(item))
        CALL map_cell_position(item)
        item => item%next
      ENDDO
      WRITE(nim_wr,*)"Global Cell Link List Built"
c-----------------------------------------------------------------------
c
c						Determine the magnetic axis.
      IF(gridshape /= 'rect')THEN
        SELECT CASE(find_axis)
        CASE('read')
           rmaxis=rmguess
           zmaxis=zmguess
        CASE default
           CALL compute_magnetic_axis
        END SELECT
      ENDIF
c
c						Select the computational task.
      SELECT CASE(task)
      CASE('poincare')
        SELECT CASE(poincare_positions)
        CASE('read')
          CALL read_start_positions
        CASE('qval')
           CALL q_start_positions
        CASE('cell')
          CALL cell_start_positions
        CASE('axis')
          CALL axis_start_positions
        CASE default
          CALL cell_start_positions
        END SELECT
        CALL poincare
      CASE('ntm')
        CALL threshold
      CASE('poloidal')
        CALL axis_start_positions
        CALL psigrid
        CALL polfft
      CASE default
        WRITE(nim_wr,*)"Unknown task ",task
        WRITE(nim_wr,*)"Valid tasks are :"
        WRITE(nim_wr,*)"'poincare','ntm','poloidal'"
      END SELECT


      CALL nim_stop('Normal termination.')

      END PROGRAM nimfl


c-----------------------------------------------------------------------
c     nim_stop utility routine.
c-----------------------------------------------------------------------
      SUBROUTINE nim_stop(message)
      USE local
      USE input
      USE global
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: message
      INTEGER(i4) :: ierror
c
c
c     					write completion message.
      IF (node==0)
     $  WRITE(nim_wr,'(2a)') 'NIM_STOP => ', TRIM(message)
c
      IF (nprocs>1) CALL mpi_finalize(ierror)
      STOP
      END SUBROUTINE nim_stop

