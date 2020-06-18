c-----------------------------------------------------------------------
c     file nimcomb.f
c     this program allows a user to read a nimrod dump file, and then
c     add the perturbed fields from another dump file.  that the
c     two dump files have matching data structures (blocks, elements,
c     Fourier components, polynomial degree, etc.) is assumed.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  nimcomb.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c     1.  main program, nimcomb
c-----------------------------------------------------------------------
      PROGRAM nimcomb
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
      REAL(r8) :: dump_time,dx,dy,pscale
      REAL(r8), DIMENSION(:), POINTER :: keff
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x_node

      TYPE(rblock_type), DIMENSION(:), POINTER :: rbn
      TYPE(cvector_type) :: cvec1,cvec2
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
      WRITE(nim_wr,*) "Enter the complete name of the complete dump,"
      WRITE(nim_wr,*) "file and the root name for the new output dump."
      READ(nim_rd,*) old_file,new_file
      dump_file=old_file
c-----------------------------------------------------------------------
c     read the old dump file.
c-----------------------------------------------------------------------
      CALL dumpb_read(nmodes,keff,dump_time,dump_step)
      IF (nrbl<=0) CALL nim_stop("There are no rblocks.")
c-----------------------------------------------------------------------
c     report poly_degree and number of Fourier components.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) TRIM(old_file)," has poly_degree=",poly_degree
      WRITE(nim_wr,*) "and ",nmodes," Fourier components."

      IF (dump_step>999999) THEN
        WRITE(step_name,fmt='(i7.7)') dump_step
      ELSE IF (dump_step>99999) THEN
        WRITE(step_name,fmt='(i6.6)') dump_step
      ELSE
        WRITE(step_name,fmt='(i5.5)') dump_step
      ENDIF
c-----------------------------------------------------------------------
c     allocate the new rblocks and copy all equilibrium and perturbed
c     data.
c-----------------------------------------------------------------------
      ALLOCATE(rbn(nrbl))
      new_pd=poly_degree
      DO ib=1,nrbl
        mxb=rb(ib)%mx
        myb=rb(ib)%my
        rbn(ib)%id=rb(ib)%id
        rbn(ib)%mx=rb(ib)%mx
        rbn(ib)%my=rb(ib)%my
        rbn(ib)%degenerate=rb(ib)%degenerate

        CALL lagr_quad_alloc(rbn(ib)%rz,mxb,myb,2_i4,new_pd,
     $    name=rb(ib)%rz%name,title=rb(ib)%rz%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%rz,rb(ib)%rz)
        CALL lagr_quad_alloc(rbn(ib)%be_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%be_eq%name,title=rb(ib)%be_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%be_eq,rb(ib)%be_eq)
        CALL lagr_quad_alloc(rbn(ib)%ja_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%ja_eq%name,title=rb(ib)%ja_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%ja_eq,rb(ib)%ja_eq)
        CALL lagr_quad_alloc(rbn(ib)%ve_eq,mxb,myb,3_i4,new_pd,
     $    name=rb(ib)%ve_eq%name,title=rb(ib)%ve_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%ve_eq,rb(ib)%ve_eq)
        CALL lagr_quad_alloc(rbn(ib)%pres_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%pres_eq%name,title=rb(ib)%pres_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%pres_eq,rb(ib)%pres_eq)
        CALL lagr_quad_alloc(rbn(ib)%prese_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%prese_eq%name,title=rb(ib)%prese_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%prese_eq,rb(ib)%prese_eq)
        CALL lagr_quad_alloc(rbn(ib)%nd_eq,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%nd_eq%name,title=rb(ib)%nd_eq%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%nd_eq,rb(ib)%nd_eq)
        CALL lagr_quad_alloc(rbn(ib)%diff_shape,mxb,myb,1_i4,new_pd,
     $    name=rb(ib)%diff_shape%name,title=rb(ib)%diff_shape%title)
          CALL lagr_quad_2D_assign_laq(rbn(ib)%diff_shape,
     $                                 rb(ib)%diff_shape)

        CALL lagr_quad_alloc(rbn(ib)%be,mxb,myb,3_i4,nmodes,new_pd,
     $    name=rb(ib)%be%name,title=rb(ib)%be%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%be,rb(ib)%be)
        CALL lagr_quad_alloc(rbn(ib)%ve,mxb,myb,3_i4,nmodes,new_pd,
     $    name=rb(ib)%ve%name,title=rb(ib)%ve%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%ve,rb(ib)%ve)
        CALL lagr_quad_alloc(rbn(ib)%pres,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%pres%name,title=rb(ib)%pres%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%pres,rb(ib)%pres)
        CALL lagr_quad_alloc(rbn(ib)%prese,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%prese%name,title=rb(ib)%prese%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%prese,rb(ib)%prese)
        CALL lagr_quad_alloc(rbn(ib)%tele,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%tele%name,title=rb(ib)%tele%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%tele,rb(ib)%tele)
        CALL lagr_quad_alloc(rbn(ib)%tion,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%tion%name,title=rb(ib)%tion%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%tion,rb(ib)%tion)
        CALL lagr_quad_alloc(rbn(ib)%nd,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%nd%name,title=rb(ib)%nd%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%nd,rb(ib)%nd)
        CALL lagr_quad_alloc(rbn(ib)%conc,mxb,myb,1_i4,nmodes,new_pd,
     $    name=rb(ib)%conc%name,title=rb(ib)%conc%title)
          CALL lagr_quad_3D_assign_laq(rbn(ib)%conc,rb(ib)%conc)
      ENDDO
c-----------------------------------------------------------------------
c     collect the perturbation data and a scaling factor for it.
c-----------------------------------------------------------------------
      WRITE(nim_wr,*) "Enter the name of the dump file that holds"
      WRITE(nim_wr,*) "solution fields to add to the first dump."
      READ(nim_rd,IOSTAT=read_stat,FMT='(a)') dump_file
      CALL dumpb_read(nmodes,keff,dump_time,dump_step)
      WRITE(nim_wr,*) " "
      WRITE(nim_wr,*) "Enter the scaling factor for the new data."
      READ(nim_rd,*) pscale
      WRITE(nim_wr,*) " "
c-----------------------------------------------------------------------
c     loop over blocks and combine the solution fields in the rbn
c     blocks.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%be)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%be)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%ve)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%ve)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%pres)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%pres)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%prese)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%prese)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%tele)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%tele)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%tion)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%tion)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%nd)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%nd)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
        CALL cvector_ptassign_laq(cvec1,rbn(ib)%conc)
        CALL cvector_ptassign_laq(cvec2,rb (ib)%conc)
        CALL vector_add(cvec1,cvec2,v2fac=pscale)
      ENDDO
c-----------------------------------------------------------------------
c     write new dump file.
c-----------------------------------------------------------------------
      dump_name=new_file
      CALL dumpb_write(nmodes,keff,dump_time,dump_step,rbn,tb)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL nim_stop('Normal termination.')

      END PROGRAM nimcomb
