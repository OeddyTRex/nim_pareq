c-----------------------------------------------------------------------
c     file nimbase.f
c     this program allows a user to read a nimrod dump file, change the
c     degree or node distribution of the polynomial basis functions, and
c     writes the result in a new dump file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  nimbase.
c     2.  laq_comp.
c     3.  laq2_comp.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     1.  main program, nimbase
c-----------------------------------------------------------------------
      PROGRAM nimbase
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE dump
      USE fields
      USE input
      IMPLICIT NONE

      CHARACTER(128) :: old_file,new_file
      CHARACTER(8) :: old_dist,new_dist,step_name
      CHARACTER(1) :: ans
      INTEGER(i4) :: nmodes,dump_step,ib,ix,iy,new_pd,mxb,myb,ibasis,im
      INTEGER :: num_chars,read_stat
      REAL(r8) :: dump_time,dx,dy
      REAL(r8), DIMENSION(:), POINTER :: keff
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x_node

      TYPE(rblock_type), DIMENSION(:), POINTER :: rbn
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rsarr,rzarr,rvarr
      REAL(r8), DIMENSION(1,1,1) :: rdarr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: csarr,cvarr
      COMPLEX(r8), DIMENSION(1,1,1,1) :: cdarr
c-----------------------------------------------------------------------
c     interface block for dumpb_read.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dumpb_read(nmodes,keff,dump_time,dump_step)
        USE local
        INTEGER(i4), INTENT(OUT) :: nmodes
        INTEGER(i4), INTENT(OUT) :: dump_step
        REAL(r8), INTENT(OUT) :: dump_time
        REAL(r8), DIMENSION(:), POINTER :: keff
        END SUBROUTINE dumpb_read
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for dumpb_write.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE dumpb_write(nmodes,keff,dump_time,dump_step,rbl,tbl)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        INTEGER(i4), INTENT(IN) :: nmodes,dump_step
        REAL(r8), INTENT(IN) :: dump_time
        REAL(r8), DIMENSION(:), POINTER :: keff
        TYPE(rblock_type), DIMENSION(:), POINTER :: rbl
        TYPE(tblock_type), DIMENSION(:), POINTER :: tbl
        END SUBROUTINE dumpb_write
      END INTERFACE
c-----------------------------------------------------------------------
c     get files names and read old dump file.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) " "
      WRITE(nim_wr,*) "Enter the complete name of the old dump file,"
      WRITE(nim_wr,*) "and the root name for the new dump file."
      READ(nim_rd,*) old_file,new_file

      dump_file=old_file
c-----------------------------------------------------------------------
c     check the default basis distribution for this build.
c-----------------------------------------------------------------------
      CALL poly_inquire(old_dist)
      WRITE(nim_wr,*) " "
      WRITE(nim_wr,*)
     $  "The default node distribution for this code-suite build is ",
     $  "'",TRIM(old_dist),"'"
      WRITE(nim_wr,*) "Has the data in ",TRIM(old_file)
      WRITE(nim_wr,*) "been generated with this distribution? (","'",
     $  "y","'"," or ","'","n","'",")"
      READ(nim_rd,*) ans
      WRITE(nim_wr,*) " "

      IF (ans=='n') THEN
        IF (old_dist=='uniform') THEN
          old_dist='gll'
        ELSE
          old_dist='uniform'
        ENDIF
        CALL poly_set(old_dist)
      ENDIF

      WRITE(nim_wr,*) "'",TRIM(old_dist),"'",
     $  " nodes will be used to interpret ",TRIM(old_file)
      WRITE(nim_wr,*) " "

c-----------------------------------------------------------------------
c     read the old dump file.
c-----------------------------------------------------------------------
      CALL dumpb_read(nmodes,keff,dump_time,dump_step)
      IF (nrbl<=0) CALL nim_stop("There are no rblocks.")
c-----------------------------------------------------------------------
c     report poly_degree and request the new value.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) TRIM(old_file)," has poly_degree=",poly_degree
      WRITE(nim_wr,*) " "
      WRITE(nim_wr,*) "Enter the desired value of poly_degree and the"
      WRITE(nim_wr,*) "new node distribution (","'","uniform","'",
     $  " or ","'","gll","'",")."
      READ(nim_rd,*) new_pd,new_dist
      WRITE(nim_wr,*) " "

      IF (dump_step>999999) THEN
        WRITE(step_name,fmt='(i7.7)') dump_step
      ELSE IF (dump_step>99999) THEN
        WRITE(step_name,fmt='(i6.6)') dump_step
      ELSE
        WRITE(step_name,fmt='(i5.5)') dump_step
      ENDIF

      WRITE(nim_wr,*) TRIM(new_file)//"."//TRIM(step_name),
     $  " will be written for poly_degree=",new_pd," polynomials"
      WRITE(nim_wr,*) "with ","'",TRIM(new_dist),"'",
     $  " node distribution."
      WRITE(nim_wr,*) " "
c-----------------------------------------------------------------------
c     allocate the new rblocks with node locations set according to
c     new_dist.
c-----------------------------------------------------------------------
      CALL poly_set(new_dist)
      ALLOCATE(rbn(nrbl))
      DO ib=1,nrbl
        mxb=rb(ib)%mx
        myb=rb(ib)%my
        rbn(ib)%id=rb(ib)%id
        rbn(ib)%mx=rb(ib)%mx
        rbn(ib)%my=rb(ib)%my
        rbn(ib)%degenerate=rb(ib)%degenerate

        CALL lagr_quad_alloc(rbn(ib)%rz,mxb,myb,2_i4,new_pd,
     $    name=rb(ib)%rz%name,title=rb(ib)%rz%title)
        CALL lagr_quad_alloc(rbn(ib)%be_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%be_eq%name,title=rb(ib)%be_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%ja_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%ja_eq%name,title=rb(ib)%ja_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%ve_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%ve_eq%name,title=rb(ib)%ve_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%pres_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%pres_eq%name,title=rb(ib)%pres_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%prese_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%prese_eq%name,title=rb(ib)%prese_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%nd_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%nd_eq%name,title=rb(ib)%nd_eq%title)
        CALL lagr_quad_alloc(rbn(ib)%diff_shape,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%diff_shape%name,title=rb(ib)%diff_shape%title)

        CALL lagr_quad_alloc(rbn(ib)%be,mxb,myb,3_i4,nmodes,new_pd,
     $    name=rb(ib)%be%name,title=rb(ib)%be%title)
        CALL lagr_quad_alloc(rbn(ib)%ve,mxb,myb,3_i4,nmodes,new_pd,
     $    name=rb(ib)%ve%name,title=rb(ib)%ve%title)
        CALL lagr_quad_alloc(rbn(ib)%pres,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%pres%name,title=rb(ib)%pres%title)
        CALL lagr_quad_alloc(rbn(ib)%prese,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%prese%name,title=rb(ib)%prese%title)
        CALL lagr_quad_alloc(rbn(ib)%tele,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%tele%name,title=rb(ib)%tele%title)
        CALL lagr_quad_alloc(rbn(ib)%tion,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%tion%name,title=rb(ib)%tion%title)
        CALL lagr_quad_alloc(rbn(ib)%nd,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%nd%name,title=rb(ib)%nd%title)
        CALL lagr_quad_alloc(rbn(ib)%conc,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%conc%name,title=rb(ib)%conc%title)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the old data at the locations of the new nodes.  the
c     node_dist variable is set back to what's appropriate for the old
c     dump file, but the node locations in the new data structures
c     determine the evaluation points.
c-----------------------------------------------------------------------
      CALL poly_set(old_dist)
      DO ib=1,nrbl
        mxb=rbn(ib)%mx
        myb=rbn(ib)%my

        DO ibasis=1,SIZE(rbn(ib)%rz%ix0)
          DO iy=rbn(ib)%rz%iy0(ibasis),myb
            DO ix=rbn(ib)%rz%ix0(ibasis),mxb
              dx=ix-rbn(ib)%rz%ix0(ibasis)+rbn(ib)%rz%dx(ibasis)
              dy=iy-rbn(ib)%rz%iy0(ibasis)+rbn(ib)%rz%dy(ibasis)

              CALL lagr_quad_eval(rb(ib)%rz,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%be_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%ja_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%ve_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%pres_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%prese_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%nd_eq,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%diff_shape,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%be,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%ve,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%pres,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%prese,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%tele,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%tion,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%nd,dx,dy,0_i4)
              CALL lagr_quad_eval(rb(ib)%conc,dx,dy,0_i4)

              CALL lagr_quad_basis_assign_loc(rbn(ib)%rz,
     $          rb(ib)%rz%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%be_eq,
     $          rb(ib)%be_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%ja_eq,
     $          rb(ib)%ja_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%ve_eq,
     $          rb(ib)%ve_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%pres_eq,
     $          rb(ib)%pres_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%prese_eq,
     $          rb(ib)%prese_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%nd_eq,
     $          rb(ib)%nd_eq%f,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rbn(ib)%diff_shape,
     $          rb(ib)%diff_shape%f,ibasis,ix,iy)

              DO im=1,nmodes
                CALL lagr_quad_basis_assign_loc(rbn(ib)%be,
     $            rb(ib)%be%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%ve,
     $            rb(ib)%ve%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%pres,
     $            rb(ib)%pres%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%prese,
     $            rb(ib)%prese%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%tele,
     $            rb(ib)%tele%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%tion,
     $            rb(ib)%tion%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%nd,
     $            rb(ib)%nd%f(:,im),ibasis,ix,iy,im)
                CALL lagr_quad_basis_assign_loc(rbn(ib)%conc,
     $            rb(ib)%conc%f(:,im),ibasis,ix,iy,im)
              ENDDO

            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     data comparison option.  re-read old dump file afterwards.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*)
     $  "As a test, you may enter the name of another dump"
      WRITE(nim_wr,*)
     $  "file to compare against the newly interpolated data."
      WRITE(nim_wr,*)
     $  "Hit return to skip this."
      WRITE(nim_wr,'(a)',ADVANCE='NO') 'Test dump file? '
      READ(nim_rd,IOSTAT=read_stat,FMT='(a)',ADVANCE='NO',
     $           SIZE=num_chars) dump_file

      IF (num_chars>0) THEN
        WRITE(nim_wr,*) " "
        WRITE(nim_wr,*) "Enter the appropriate node distribution for"
        WRITE(nim_wr,*) TRIM(dump_file)," (","'","uniform","'",
     $    " or ","'","gll","'",")."
        READ(nim_rd,*) new_dist
        CALL poly_set(new_dist)
        CALL dumpb_read(nmodes,keff,dump_time,dump_step)
        OPEN(UNIT=out_unit,FILE='nimbase.out',STATUS='UNKNOWN')
        DO ib=1,nrbl
          CALL laq2_comp(rb(ib)%rz,rbn(ib)%rz,ib)
          CALL laq2_comp(rb(ib)%be_eq,rbn(ib)%be_eq,ib)
          CALL laq2_comp(rb(ib)%ja_eq,rbn(ib)%ja_eq,ib)
          CALL laq2_comp(rb(ib)%ve_eq,rbn(ib)%ve_eq,ib)
          CALL laq2_comp(rb(ib)%pres_eq,rbn(ib)%pres_eq,ib)
          CALL laq2_comp(rb(ib)%prese_eq,rbn(ib)%prese_eq,ib)
          CALL laq2_comp(rb(ib)%nd_eq,rbn(ib)%nd_eq,ib)
          CALL laq2_comp(rb(ib)%diff_shape,rbn(ib)%diff_shape,ib)
          CALL laq_comp(rb(ib)%be,rbn(ib)%be,ib)
          CALL laq_comp(rb(ib)%ve,rbn(ib)%ve,ib)
          CALL laq_comp(rb(ib)%pres,rbn(ib)%pres,ib)
          CALL laq_comp(rb(ib)%prese,rbn(ib)%prese,ib)
          CALL laq_comp(rb(ib)%tele,rbn(ib)%tele,ib)
          CALL laq_comp(rb(ib)%tion,rbn(ib)%tion,ib)
          CALL laq_comp(rb(ib)%nd,rbn(ib)%nd,ib)
          CALL laq_comp(rb(ib)%conc,rbn(ib)%conc,ib)
        ENDDO
        CLOSE(out_unit)
        dump_file=old_file
        CALL dumpb_read(nmodes,keff,dump_time,dump_step)
      ENDIF
c-----------------------------------------------------------------------
c     write new dump file.
c-----------------------------------------------------------------------
      dump_name=new_file
      poly_degree=new_pd
      CALL dumpb_write(nmodes,keff,dump_time,dump_step,rbn,tb)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')

      END PROGRAM nimbase

c-----------------------------------------------------------------------
c     subprogram 2. laq_comp.
c     compare two 3D lagrange_quad data structures.
c-----------------------------------------------------------------------
      SUBROUTINE laq_comp(laq1,laq2,ibl)
      USE local
      USE lagr_quad_mod
      IMPLICIT NONE
      
      TYPE(lagr_quad_type), INTENT(IN) :: laq1,laq2
      INTEGER(i4), INTENT(IN) :: ibl
    
      INTEGER(i4) :: ii
c-----------------------------------------------------------------------
c     header info for 1.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/,a,i4,/)') "BLOCK ",ibl
      WRITE(out_unit,'(2a)') "laq1 ",laq1%name
      WRITE(out_unit,*) laq1%title
      WRITE(out_unit,*) laq1%mx,laq1%my,laq1%nqty,laq1%nfour,
     $                  laq1%n_side,laq1%n_int
      DO ii=1,SIZE(laq1%ix0)
        WRITE(out_unit,*) ii,laq1%ix0(ii),laq1%iy0(ii),laq1%dx(ii),
     $                    laq1%dy(ii)
      ENDDO
c-----------------------------------------------------------------------
c     header info for 2.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(2a)') "laq2 ",laq2%name
      WRITE(out_unit,*) laq2%title
      WRITE(out_unit,*) laq2%mx,laq2%my,laq2%nqty,laq2%nfour,
     $                  laq2%n_side,laq2%n_int
      DO ii=1,SIZE(laq2%ix0)
        WRITE(out_unit,*) ii,laq2%ix0(ii),laq2%iy0(ii),laq2%dx(ii),
     $                    laq2%dy(ii)
      ENDDO
c-----------------------------------------------------------------------
c     compare data.
c-----------------------------------------------------------------------
      WRITE(out_unit,*) " "
      WRITE(out_unit,*) "Grid vertices ",MAXVAL(ABS(laq1%fs-laq2%fs)),
     $                  MAXLOC(ABS(laq1%fs-laq2%fs))
      IF (laq1%n_side>0) THEN
        WRITE(out_unit,*) "Horiz. sides  ",
     $                    MAXVAL(ABS(laq1%fsh-laq2%fsh)),
     $                    MAXLOC(ABS(laq1%fsh-laq2%fsh))
        WRITE(out_unit,*) "Vert.  sides  ",
     $                    MAXVAL(ABS(laq1%fsv-laq2%fsv)),
     $                    MAXLOC(ABS(laq1%fsv-laq2%fsv))
      ENDIF
      IF (laq1%n_int>0) THEN
        WRITE(out_unit,*) "Inter.  nodes ",
     $                    MAXVAL(ABS(laq1%fsi-laq2%fsi)),
     $                    MAXLOC(ABS(laq1%fsi-laq2%fsi))
      ENDIF
      WRITE(out_unit,*) " "
c-----------------------------------------------------------------------
c     finish.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE laq_comp

c-----------------------------------------------------------------------
c     subprogram 3. laq2_comp.
c     compare two 2D lagrange_quad data structures.
c-----------------------------------------------------------------------
      SUBROUTINE laq2_comp(laq1,laq2,ibl)
      USE local
      USE lagr_quad_mod
      IMPLICIT NONE
      
      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq1,laq2
      INTEGER(i4), INTENT(IN) :: ibl
    
      INTEGER(i4) :: ii
c-----------------------------------------------------------------------
c     header info for 1.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(/,a,i4,/)') "BLOCK ",ibl
      WRITE(out_unit,'(2a)') "laq1 ",laq1%name
      WRITE(out_unit,*) laq1%title
      WRITE(out_unit,*) laq1%mx,laq1%my,laq1%nqty,
     $                  laq1%n_side,laq1%n_int
      DO ii=1,SIZE(laq1%ix0)
        WRITE(out_unit,*) ii,laq1%ix0(ii),laq1%iy0(ii),laq1%dx(ii),
     $                    laq1%dy(ii)
      ENDDO
c-----------------------------------------------------------------------
c     header info for 2.
c-----------------------------------------------------------------------
      WRITE(out_unit,'(2a)') "laq2 ",laq2%name
      WRITE(out_unit,*) laq2%title
      WRITE(out_unit,*) laq2%mx,laq2%my,laq2%nqty,
     $                  laq2%n_side,laq2%n_int
      DO ii=1,SIZE(laq2%ix0)
        WRITE(out_unit,*) ii,laq2%ix0(ii),laq2%iy0(ii),laq2%dx(ii),
     $                    laq2%dy(ii)
      ENDDO
c-----------------------------------------------------------------------
c     compare data.
c-----------------------------------------------------------------------
      WRITE(out_unit,*) " "
      WRITE(out_unit,*) "Grid vertices ",MAXVAL(ABS(laq1%fs-laq2%fs)),
     $                  MAXLOC(ABS(laq1%fs-laq2%fs))
      IF (laq1%n_side>0) THEN
        WRITE(out_unit,*) "Horiz. sides  ",
     $                    MAXVAL(ABS(laq1%fsh-laq2%fsh)),
     $                    MAXLOC(ABS(laq1%fsh-laq2%fsh))
        WRITE(out_unit,*) "Vert.  sides  ",
     $                    MAXVAL(ABS(laq1%fsv-laq2%fsv)),
     $                    MAXLOC(ABS(laq1%fsv-laq2%fsv))
      ENDIF
      IF (laq1%n_int>0) THEN
        WRITE(out_unit,*) "Inter.  nodes ",
     $                    MAXVAL(ABS(laq1%fsi-laq2%fsi)),
     $                    MAXLOC(ABS(laq1%fsi-laq2%fsi))
      ENDIF
      WRITE(out_unit,*) " "
c-----------------------------------------------------------------------
c     finish.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE laq2_comp
