c-----------------------------------------------------------------------
c     subprogram contour.
c     writes data for generalized contour plots.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. contour_mod.
c     1. contour_init.
c     2. contour_step.
c     3. contour_laq_write.
c     4. contour_laq2_write.
c     5. contour_bicube_write.
c     6. contour_tl_write.
c     7. contour_tl2_write.
c-----------------------------------------------------------------------
c     subprogram 0. contour_mod.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE contour_mod
      USE local
      USE fields
      IMPLICIT NONE

      INTEGER(i4), PRIVATE :: nrb,ntb,nqty

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. contour_init.
c     writes coordinates for generalized contour plots.
c-----------------------------------------------------------------------
      SUBROUTINE contour_init(nqty,ndcon,filename)
      
      INTEGER(i4), INTENT(IN) :: nqty,ndcon
      CHARACTER(*), INTENT(IN) :: filename

      INTEGER(i4) :: ib,pd
c-----------------------------------------------------------------------
c     open file, and write the number of blocks and quantities.
c-----------------------------------------------------------------------
      CALL open_bin(con_unit,filename,"UNKNOWN","REWIND",32_i4)
      WRITE(con_unit) INT(nrbl,4),INT(nbl-nrbl,4),
     $                INT(nqty,4)
c-----------------------------------------------------------------------
c     write sizes and coordinates for rblocks with the specified
c     amount of data per element.
c-----------------------------------------------------------------------
      pd=ndcon
      DO ib=1,nrbl
        WRITE(con_unit) INT(pd*rb(ib)%mx,4),INT(pd*rb(ib)%my,4)
        CALL contour_laq2_write(rb(ib)%rz,pd,one_record=.true.)
      ENDDO
c-----------------------------------------------------------------------
c     write sizes and coordinates for tblocks.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         WRITE(con_unit) INT(tb(ib)%mvert,4),INT(tb(ib)%mcell,4)
         WRITE(con_unit) REAL(tb(ib)%tgeom%xs,4),REAL(tb(ib)%tgeom%ys,4)
         WRITE(con_unit) INT(tb(ib)%tgeom%vertex,4)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_init
c-----------------------------------------------------------------------
c     subprogram 2. contour_step.
c     writes dependent variables for generalized contour plots
c     for each time step.
c-----------------------------------------------------------------------
      SUBROUTINE contour_step(mode_st,mode_en,ndcon)
      
      INTEGER(i4), INTENT(IN) :: mode_st,mode_en,ndcon

      INTEGER(i4) :: ib,pd
c-----------------------------------------------------------------------
c     write components for rblocks with the specified amount of data
c     per element.
c-----------------------------------------------------------------------
      pd=ndcon
      DO ib=1,nrbl
        CALL contour_laq2_write(rb(ib)%be_eq,pd)
        CALL contour_laq2_write(rb(ib)%ja_eq,pd)
        CALL contour_laq2_write(rb(ib)%ve_eq,pd)
        CALL contour_laq2_write(rb(ib)%pres_eq,pd)
        CALL contour_laq2_write(rb(ib)%prese_eq,pd)
        CALL contour_laq2_write(rb(ib)%nd_eq,pd)
        CALL contour_laq2_write(rb(ib)%diff_shape,pd)
        CALL contour_laq_write(rb(ib)%be,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%ja,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%ve,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%pres,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%prese,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%nd,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%conc,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%tele,mode_st,mode_en,pd)
        CALL contour_laq_write(rb(ib)%tion,mode_st,mode_en,pd)
      ENDDO
c-----------------------------------------------------------------------
c     write components for tblocks.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
        CALL contour_tl2_write(tb(ib)%be_eq)
        CALL contour_tl2_write(tb(ib)%ja_eq)
        CALL contour_tl2_write(tb(ib)%ve_eq)
        CALL contour_tl2_write(tb(ib)%pres_eq)
        CALL contour_tl2_write(tb(ib)%prese_eq)
        CALL contour_tl2_write(tb(ib)%nd_eq)
        CALL contour_tl2_write(tb(ib)%diff_shape)
        CALL contour_tl_write(tb(ib)%be,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%ja,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%ve,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%pres,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%prese,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%nd,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%conc,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%tele,mode_st,mode_en)
        CALL contour_tl_write(tb(ib)%tion,mode_st,mode_en)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_step
c-----------------------------------------------------------------------
c     subprogram 3. contour_laq_write.
c     writes one block of lagrange quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE contour_laq_write(laq,im0,im1,pd,nqlim)

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      INTEGER(i4), INTENT(IN) :: im0,im1,pd
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqlim

      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp,carr
      COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
      INTEGER(i4) :: iv,nv,nvwr,im,nm,ix,iy,mx,my
      REAL(r8) :: dx,dy,dl
c-----------------------------------------------------------------------
c     copy to the old array index format.
c-----------------------------------------------------------------------
      dl=1._r8/pd
      nv=laq%nqty
      IF (PRESENT(nqlim)) THEN
        nvwr=nqlim
      ELSE
        nvwr=nv
      ENDIF
      nm=laq%nfour
      mx=laq%mx
      my=laq%my
      ALLOCATE(ctmp(nv,mx,my,nm),carr(nv,0:pd*mx,0:pd*my,nm))
      DO iy=0,pd
        dy=dl*iy
        DO ix=0,pd
          dx=dl*ix
          CALL lagr_quad_all_eval(laq,dx,dy,ctmp,dc,dc,0_i4)
          carr(:,ix:pd*(mx-1)+ix:pd,iy:pd*(my-1)+iy:pd,:)=ctmp
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write contour data.
c-----------------------------------------------------------------------
      DO im=im0,im1
        DO iv=1,nvwr
          WRITE(con_unit) REAL(carr(iv,:,:,im),4)
        ENDDO
        DO iv=1,nvwr
          WRITE(con_unit) REAL(-(0,1)*carr(iv,:,:,im),4)
        ENDDO
      ENDDO
      DEALLOCATE(ctmp,carr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_laq_write
c-----------------------------------------------------------------------
c     subprogram 4. contour_laq2_write.
c     writes one block of 2D lagrange quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE contour_laq2_write(laq,pd,nqlim,one_record)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      INTEGER(i4), INTENT(IN) :: pd
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqlim
      LOGICAL, OPTIONAL :: one_record

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rtmp,rarr
      REAL(r8), DIMENSION(1,1,1) :: dr
      INTEGER(i4) :: iv,nv,nvwr,im,nm,ix,iy,mx,my
      REAL(r8) :: dx,dy,dl
c-----------------------------------------------------------------------
c     copy to the old array index format.
c-----------------------------------------------------------------------
      dl=1._r8/pd
      nv=laq%nqty
      IF (PRESENT(nqlim)) THEN
        nvwr=nqlim
      ELSE
        nvwr=nv
      ENDIF
      mx=laq%mx
      my=laq%my
      ALLOCATE(rtmp(nv,mx,my),rarr(nv,0:pd*mx,0:pd*my))
      DO iy=0,pd
        dy=dl*iy
        DO ix=0,pd
          dx=dl*ix
          CALL lagr_quad_all_eval(laq,dx,dy,rtmp,dr,dr,0_i4)
          rarr(:,ix:pd*(mx-1)+ix:pd,iy:pd*(my-1)+iy:pd)=rtmp
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write contour data.
c-----------------------------------------------------------------------
      IF (PRESENT(one_record)) THEN
        WRITE(con_unit) (REAL(rarr(iv,:,:),4),iv=1,nv)
      ELSE
        DO iv=1,nvwr
          WRITE(con_unit) REAL(rarr(iv,:,:),4)
        ENDDO
      ENDIF
      DEALLOCATE(rtmp,rarr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_laq2_write
c-----------------------------------------------------------------------
c     subprogram 5. contour_bicube_write.
c     writes one block of bicubic data.
c-----------------------------------------------------------------------
      SUBROUTINE contour_bicube_write(bc,pd,one_record)

      TYPE(bicube_type), INTENT(INOUT) :: bc
      INTEGER(i4), INTENT(IN) :: pd
      LOGICAL, OPTIONAL :: one_record

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rtmp,rarr
      REAL(r8), DIMENSION(1,1,1) :: db
      INTEGER(i4) :: iv,nv,ix,iy,mx,my
      REAL(r8) :: dx,dy,dl
c-----------------------------------------------------------------------
c     write contour data.
c-----------------------------------------------------------------------
      dl=1._r8/pd
      nv=bc%nqty
      mx=bc%mx
      my=bc%my
      ALLOCATE(rtmp(nv,mx,my),rarr(nv,0:pd*mx,0:pd*my))
      DO iy=0,pd
        dy=dl*iy
        DO ix=0,pd
          dx=dl*ix
          CALL bicube_all_eval(bc,dx,dy,rtmp,db,db,db,db,db,0_i4)
          rarr(:,ix:pd*(mx-1)+ix:pd,iy:pd*(my-1)+iy:pd)=rtmp
        ENDDO
      ENDDO
      DEALLOCATE(bc%cmats)
      IF (PRESENT(one_record)) THEN
        WRITE(con_unit) (REAL(rarr(iv,:,:),4),iv=1,nv)
      ELSE
        DO iv=1,nv
          WRITE(con_unit) REAL(rarr(iv,:,:),4)
        ENDDO
      ENDIF
      DEALLOCATE(rtmp,rarr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_bicube_write
c-----------------------------------------------------------------------
c     subprogram 6. contour_tl_write.
c     writes one block of 3D tri_linear data.
c-----------------------------------------------------------------------
      SUBROUTINE contour_tl_write(tl,im0,im1,nqlim)

      TYPE(tri_linear_type), INTENT(IN) :: tl
      INTEGER(i4), INTENT(IN) :: im0,im1
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqlim

      INTEGER(i4) :: iv,nv,im
c-----------------------------------------------------------------------
c     write contour data.
c-----------------------------------------------------------------------
      IF (PRESENT(nqlim)) THEN
        nv=nqlim
      ELSE
        nv=tl%nqty
      ENDIF
      DO im=im0,im1
        DO iv=1,nv
          WRITE(con_unit) REAL(tl%fs(iv,:,:,im),4)
        ENDDO
        DO iv=1,nv
          WRITE(con_unit) REAL(-(0,1)*tl%fs(iv,:,:,im),4)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_tl_write
c-----------------------------------------------------------------------
c     subprogram 7. contour_tl2_write.
c     writes one block of 2D tri_linear data.
c-----------------------------------------------------------------------
      SUBROUTINE contour_tl2_write(tl)

      TYPE(tri_linear_2D_type), INTENT(IN) :: tl

      INTEGER(i4) :: iv,nv
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      nv=tl%nqty
      DO iv=1,nv
        WRITE(con_unit) REAL(tl%fs(iv,:,:),4)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE contour_tl2_write
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE contour_mod
