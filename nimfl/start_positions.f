      MODULE node_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: node_type
c       INTEGER(i4) :: m,n
        REAL(r8) :: area                        ! In lieu of psi
        INTEGER(i4) :: ibl
        REAL(r8) :: x_start,y_start
        REAL(r8),DIMENSION(2) :: rz
        TYPE(node_type), POINTER :: next    ! ptr to next node in queue
      END TYPE node_type
      TYPE(node_type), POINTER :: surface_first,surface_cnt
      TYPE(node_type), POINTER :: node_bmin,node_bmax,node_bzero

      END MODULE node_type_mod
c-----------------------------------------------------------------------
c     This module contains routines used to initialize start positions
c     for the field line integration.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  start_positions
c     1.  cell_start_positions
c     2.  read_start_positions
c     3.  q_start_positions
c     4.  axis_start_positions
c-----------------------------------------------------------------------

      MODULE start_positions
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: n_fieldlines		! actual number of field lines
      INTEGER(i4),DIMENSION(:), ALLOCATABLE  :: ibl_fieldlines	! block number
      REAL(r8),DIMENSION(:), ALLOCATABLE  :: x_fieldlines	! logical x
      REAL(r8),DIMENSION(:), ALLOCATABLE  :: y_fieldlines	! logical y

      CONTAINS


      SUBROUTINE cell_start_positions
c					Simple cell based start locations
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      IMPLICIT NONE

      INTEGER(i4) :: iline
      TYPE(cell_type), POINTER :: item


c
c     					Count number of lines
      n_fieldlines = 0
      item => start
      DO WHILE(ASSOCIATED(item%face(1)%p))
        n_fieldlines=n_fieldlines+1
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO
   
      ALLOCATE(ibl_fieldlines(n_fieldlines))
      ALLOCATE(x_fieldlines(n_fieldlines))
      ALLOCATE(y_fieldlines(n_fieldlines))

      item => start
      iline=0
      DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        ibl_fieldlines(iline) = item%ib
        x_fieldlines(iline) = item%p(1,1)
        y_fieldlines(iline) = item%p(2,1)
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO

      END SUBROUTINE cell_start_positions

      SUBROUTINE read_start_positions
c					read (r,z) positions 
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      IMPLICIT NONE
      CHARACTER(19) :: infile='start_positions.dat'
      INTEGER :: open_stat
      INTEGER(i4) :: iline,ibl
      REAL(r8) :: rtemp,ztemp
      TYPE(location_type) :: p0
      TYPE(cell_type), POINTER :: item,item_pre
      REAL(r8) :: x_start,y_start
      LOGICAL :: failure

      OPEN(UNIT=in_unit,FILE=infile,STATUS='OLD',POSITION='REWIND',
     &      IOSTAT=open_stat)

      IF (open_stat/=0)THEN 
         CALL cell_start_positions
         OPEN(UNIT=in_unit,FILE=infile,STATUS='UNKNOWN',
     &      POSITION='REWIND',IOSTAT=open_stat)
         WRITE(in_unit,*)n_fieldlines
         DO iline=1,n_fieldlines
            ibl = ibl_fieldlines(iline)
            x_start = x_fieldlines(iline)
            y_start = y_fieldlines(iline)
            CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
            WRITE(in_unit,*)rbc(ibl)%rz%f(1:2)
         ENDDO
      ELSE
         READ(in_unit,*)n_fieldlines
         ALLOCATE(ibl_fieldlines(n_fieldlines))
         ALLOCATE(x_fieldlines(n_fieldlines))
         ALLOCATE(y_fieldlines(n_fieldlines))
         item => start
         DO iline=1,n_fieldlines
           READ(in_unit,*)rtemp,ztemp
           p0%point(1)=rtemp
           p0%point(2)=ztemp
           CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
           CALL refine_cell(p0,item,failure,x_start,y_start,.TRUE.)
           ibl_fieldlines(iline)=item%ib
           x_fieldlines(iline) = x_start
           y_fieldlines(iline) = y_start
         ENDDO
      ENDIF
      CLOSE(UNIT=in_unit)

      END SUBROUTINE read_start_positions

      SUBROUTINE q_start_positions
c					Locate q surfaces and pick points
c					on that surface.
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE node_type_mod
      IMPLICIT NONE
      TYPE(location_type) :: p0
      TYPE(cell_type), POINTER :: item,item_pre
      REAL(r8) :: x_start,y_start
      LOGICAL :: failure,found
      INTEGER(i4) :: iline

c
c						Build link list of surfaces
      CALL q_to_start(found)
c
c						Transfer to an array
      IF(found)THEN
        surface_cnt => surface_first
        iline = 0
        DO 
          IF(.NOT.ASSOCIATED(surface_cnt))EXIT
          iline=iline+1
          surface_cnt => surface_cnt%next
        ENDDO
         
        n_fieldlines = iline
        ALLOCATE(ibl_fieldlines(n_fieldlines))
        ALLOCATE(x_fieldlines(n_fieldlines))
        ALLOCATE(y_fieldlines(n_fieldlines))

        iline = 0
        surface_cnt => surface_first
        DO 
          IF(.NOT.ASSOCIATED(surface_cnt))EXIT
          iline = iline+1
          ibl_fieldlines(iline)=surface_cnt%ibl
          x_fieldlines(iline) = surface_cnt%x_start
          y_fieldlines(iline) = surface_cnt%y_start
          surface_cnt => surface_cnt%next
        ENDDO

        surface_cnt => surface_first
        DO
          IF(.NOT.ASSOCIATED(surface_cnt))EXIT
          surface_first => surface_cnt
          surface_cnt => surface_cnt%next
          DEALLOCATE(surface_first)
        ENDDO
      ENDIF

      END SUBROUTINE q_start_positions


      SUBROUTINE axis_start_positions
c					Starting positions from the magnetic
c					axis to the last "good" surface 
c					with z=zmaxis for all points.
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE node_type_mod
      IMPLICIT NONE
      TYPE(location_type) :: p0,start0,end
      TYPE(cell_type), POINTER :: item,item_pre
      REAL(r8) :: x_start,y_start,rmax,deltar,r_right,r_left,r_left_old
      REAL(r8) :: area,qvalue,r_max,z_max,r_min,z_min
      REAL(r8) :: r_start,r_last,area_last
      LOGICAL :: failure,found
      INTEGER(i4) :: iline,i,ibl,i_bisect,max_bisect=1000
      REAL(r8) :: ind_tol0 = 1.0e-08
      LOGICAL :: local_debug = .FALSE.
c
c						Find the largest radius 
c						and use it to set a deltar
c						search.
      rmax = -HUGE(1.0)
      item => start
      iline=0
      DO WHILE(ASSOCIATED(item))
        iline=iline+1
        ibl= item%ib
        DO i = 1,item%nodes
          x_start = item%p(1,i)
          y_start = item%p(2,i)
          ibl = item%ib
          CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
          IF(rbc(ibl)%rz%f(1) > rmax)rmax = rbc(ibl)%rz%f(1)
        ENDDO
        item => item%next
        IF(item%id == maxcell ) EXIT		! This shouldn't be needed
        IF(item%id == 0 ) EXIT
      ENDDO
c
c					Debug a specific position.
c     
c     p0%point = (/2.742598880161820,-2.940385823193813E-03/)
c     CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
c     write(6,*)area,qvalue,p0%point(1)
c     stop 111
c
c					Determine first valid surfaces
      deltar = ABS(rmax-rmaxis)/REAL(num_surface*2,r8)
      start0%point(1)=(rmaxis+rmax)*0.5
      start0%point(2)=zmaxis
      p0%point(2)=zmaxis
      i_bisect = 0
      area_last = -HUGE(1.0)
      DO
        p0%point(1)=start0%point(1) - REAL(i_bisect,r8)*deltar
        IF(p0%point(1) < rmaxis)EXIT
        i_bisect = i_bisect+1
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        IF(area < ind_tol0)EXIT
        area_last = area
      ENDDO
      start0%point(1) = p0%point(1)+deltar
c
c					Refine the search
      deltar = ABS(rmax-rmaxis)/REAL(num_surface*num_surface**0.5,r8)
      i_bisect = 0
      DO
        p0%point(1)=start0%point(1) - REAL(i_bisect,r8)*deltar
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        IF((area < 0) .OR. (area > area_last))EXIT
        i_bisect = i_bisect+1
        area_last = area
      ENDDO
      start0%point(1) = p0%point(1)+deltar
      CALL qcompute(start0,area,qvalue,r_min,z_min,r_max,z_max)
      WRITE(nim_wr,*)"First Integrable Surface at"
      WRITE(nim_wr,*)"area,qvalue,r,z"
      WRITE(nim_wr,fmt='(4(x,E12.5))')area,qvalue,start0%point
c
c					Determine last valid surface
      deltar = ABS(rmax-rmaxis)/REAL(num_surface,r8)
      i_bisect = 1
      area_last = -HUGE(1.)
      DO
        p0%point(1)=start0%point(1) + REAL(i_bisect,r8)*deltar
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
c       write(6,*)i_bisect,area,qvalue,p0%point(1)
        IF(area < area_last)EXIT
        area_last = area
        i_bisect = i_bisect+1
      ENDDO
      end%point(1)=p0%point(1) - deltar
c
c					Refine the search for the last surface
      deltar = ABS(rmax-rmaxis)/REAL(num_surface*num_surface**0.5,r8)
      i_bisect = 0
      DO
        p0%point(1)=end%point(1) + REAL(i_bisect,r8)*deltar
        i_bisect = i_bisect+1
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        IF(area < area_last)EXIT
        area_last = area
      ENDDO
      end%point(1)=p0%point(1) - deltar
      end%point(2)=zmaxis
      CALL qcompute(end,area,qvalue,r_min,z_min,r_max,z_max)
      WRITE(nim_wr,*)"Last Integrable Surface at"
      WRITE(nim_wr,fmt='(4(x,E12.5))')area,qvalue,end%point
c
c						Diagnostic output.
      IF(local_debug)THEN
        WRITE(nim_wr,*)start0%point
        WRITE(nim_wr,*)end%point
        deltar = ABS(end%point(1)-start0%point(1))/REAL(num_surface,r8)
        DO i_bisect = 0,num_surface
          p0%point(2)=zmaxis
          p0%point(1)=start0%point(1) + REAL(i_bisect,r8)*deltar
          CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
          WRITE(nim_wr,*)p0%point(1),area,qvalue
        ENDDO
      ENDIF
c
c						Store start positions
      n_fieldlines = num_surface+1
      ALLOCATE(ibl_fieldlines(n_fieldlines))
      ALLOCATE(x_fieldlines(n_fieldlines))
      ALLOCATE(y_fieldlines(n_fieldlines))
      deltar = ABS(end%point(1)-start0%point(1))/REAL(num_surface,r8)
      item => start
      item_pre => start
      DO i_bisect = 0,num_surface
        p0%point(2)=zmaxis
        p0%point(1)=start0%point(1) + REAL(i_bisect,r8)*deltar
        CALL rz_to_cell(p0,item,item_pre,failure,.FALSE.)
        CALL refine_cell(p0,item,failure,x_start,y_start,.TRUE.)
        ibl_fieldlines(i_bisect+1)=item%ib
        x_fieldlines(i_bisect+1) = x_start
        y_fieldlines(i_bisect+1) = y_start
      ENDDO
      END SUBROUTINE axis_start_positions

      END MODULE start_positions
 
