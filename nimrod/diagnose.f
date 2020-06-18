c-----------------------------------------------------------------------
c     file diagnose.f
c     contains diagnostic routines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  global declarations.
c     1.  drawgrid.
c     2.  detail.
c     3.  xy_slice.
c     4.  xt_slice.
c     5.  yt_slice.
c     6.  time_slice_init.
c     7.  time_slice_close.
c     8.  struct_set.
c     9.  struct_dealloc.
c-----------------------------------------------------------------------
c     subprogram 0. global declarations.
c     declares everything needed for the whole module.
c-----------------------------------------------------------------------
      MODULE diagnose
      USE local
      USE input
      USE fields
      USE time
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: sbe,sja,sve,spr,spe,
     $          snd,sco,ste,sti,sbq,sjq,svq,spq,speq,sndq,srz,sdfs
      CHARACTER(64) :: xt_file,yt_file

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. drawgrid.
c     draws the grid.
c-----------------------------------------------------------------------
      SUBROUTINE drawgrid

      INTEGER(i4) :: ib,ix,iy,icell,iv,jv
c-----------------------------------------------------------------------
c     open file.
c-----------------------------------------------------------------------
      CALL timer(timestart)

      CALL open_bin(grid_unit,"grid.bin","UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     draw rblocks.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
         DO iy=0,rb(ib)%my
            DO ix=0,rb(ib)%mx
               WRITE(grid_unit)REAL(rb(ib)%rz%fs(:,ix,iy),4)
            ENDDO
            WRITE(grid_unit)
         ENDDO
         DO ix=0,rb(ib)%mx
            DO iy=0,rb(ib)%my
               WRITE(grid_unit)REAL(rb(ib)%rz%fs(:,ix,iy),4)
            ENDDO
            WRITE(grid_unit)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     draw tblocks.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         DO icell=1,tb(ib)%mcell
            DO iv=1,3
               jv=tb(ib)%tgeom%vertex(icell,iv)
               WRITE(grid_unit)REAL(tb(ib)%tgeom%xs(jv),4),
     $              REAL(tb(ib)%tgeom%ys(jv),4)
            ENDDO
            jv=tb(ib)%tgeom%vertex(icell,1)
            WRITE(grid_unit)REAL(tb(ib)%tgeom%xs(jv),4),
     $           REAL(tb(ib)%tgeom%ys(jv),4)
            WRITE(grid_unit)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CALL close_bin(grid_unit,"grid.bin")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE drawgrid
c-----------------------------------------------------------------------
c     subprogram 2. detail.
c     produces detailed ascii diagnosis of solution.
c-----------------------------------------------------------------------
      SUBROUTINE detail(istep,t)

      INTEGER(i4), INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t

      CHARACTER(128) :: format1,format2
      INTEGER(i4) :: ib,icell,ip,iqty,iv,ix,iy,iybl,nb
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT('(/2x,"iy",',i2.2,'(4x,"ix=",i2.2,2x)/)')
 20   FORMAT('(i4,1p,',i2.2,'e11.3)')
 30   FORMAT('tblock(',i1,')')
 40   FORMAT(/4x,'iv',4x,'ip',3x,'vert'/)
 50   FORMAT(3i6)
 60   FORMAT(/4x,'ic',4x,'iv1',3x,'iv2',3x,'iv3'/)
 70   FORMAT(4i6,2x)
 80   FORMAT(/4x,'iv',6x,'x',10x,'y'/)
 90   FORMAT(i6,1p,2e11.3)
c-----------------------------------------------------------------------
c     set up formats and write header.
c-----------------------------------------------------------------------
      CALL timer(timestart)

      IF(.NOT.detflag)RETURN
      IF (.NOT.out_opened) THEN
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        out_opened=.true.
      ENDIF
      nb=nxbl*nybl
      ix=0
      DO ib=1,nb,nybl
         ix=ix+rb(ib)%mx+1
      ENDDO
      WRITE(format1,10)ix
      WRITE(format2,20)ix
      WRITE(out_unit,'(/a,i3,a,1p,e9.3/)')"istep = ",istep,", t = ",t
c-----------------------------------------------------------------------
c     write rblock grid positions.
c-----------------------------------------------------------------------
      DO iqty=1,2
         WRITE(out_unit,'(/a,i1,a/)')"rz(",iqty,"):"
         WRITE(out_unit,format1)((ix,ix=0,rb(ib)%mx),ib=1,nb,nybl)
         WRITE(out_unit,format2)((iy,((rb(ib)%rz%fs(iqty,ix,iy),
     $        ix=0,rb(ib)%mx),ib=iybl,nb,nybl),
     $        iy=rb(iybl)%my,0,-1),iybl=nybl,1,-1)
         WRITE(out_unit,format1)((ix,ix=0,rb(ib)%mx),ib=1,nb,nybl)
      ENDDO
c-----------------------------------------------------------------------
c     write tblock connectivity arrays.
c-----------------------------------------------------------------------
      DO ib=nrbl+1,nbl
         WRITE(out_unit,30)ib
         WRITE(out_unit,40)
         WRITE(out_unit,50)((iv,ip,tb(ib)%tgeom%neighbor(iv)%vertex(ip),
     $        ip=0,SIZE(tb(ib)%tgeom%neighbor(iv)%vertex)-1),
     $        iv=0,tb(ib)%mvert)
         WRITE(out_unit,40)
         WRITE(out_unit,60)
         WRITE(out_unit,70)(icell,(tb(ib)%tgeom%vertex(icell,iv),
     $        iv=1,3),icell=1,tb(ib)%mcell)
         WRITE(out_unit,60)
c-----------------------------------------------------------------------
c     write tblock grid positions.
c-----------------------------------------------------------------------
         WRITE(out_unit,80)
         WRITE(out_unit,90)(iv,tb(ib)%tgeom%xs(iv),tb(ib)%tgeom%ys(iv),
     $        iv=0,tb(ib)%mvert)
         WRITE(out_unit,80)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE detail
c-----------------------------------------------------------------------
c     subprogram 3. xy_slice.
c     produces 2D binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE xy_slice(istep,t)

      INTEGER(i4), INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t

      CHARACTER(5) :: stepname
      CHARACTER(64) :: filename
      INTEGER(i4) :: ix,iy,mpx,mpy
c-----------------------------------------------------------------------
c     open file and fill pointers.
c-----------------------------------------------------------------------
      CALL timer(timestart)

      WRITE(stepname,fmt='(i5.5)')istep
      filename=TRIM(xdraw_dir)//"/xy_slice."//TRIM(stepname)//".bin"
      CALL open_bin(xy_unit,TRIM(filename),"UNKNOWN","REWIND",32_i4)
      CALL struct_set
c-----------------------------------------------------------------------
c     write binary data.
c-----------------------------------------------------------------------
      mpx=SIZE(srz,1)-1
      mpy=SIZE(srz,2)-1
      DO iy=0,mpy
        DO ix=0,mpx
          WRITE(xy_unit) (/
     $      REAL(ix,4)/mpx,REAL(iy,4)/mpy,
     $      REAL(srz(ix,iy,:),4),REAL(sbq(ix,iy,:),4),
     $      REAL(sjq(ix,iy,:),4),REAL(svq(ix,iy,:),4),
     $      REAL(spq(ix,iy,:),4),REAL(speq(ix,iy,:),4),
     $      REAL(sndq(ix,iy,:),4),REAL(sdfs(ix,iy,:),4),
     $      REAL(sbe(ix,iy,:),4),REAL(sja(ix,iy,:),4),
     $      REAL(sve(ix,iy,:),4),REAL(spr(ix,iy,:),4),
     $      REAL(spe(ix,iy,:),4),REAL(snd(ix,iy,:),4),
     $      REAL(sco(ix,iy,:),4),REAL(ste(ix,iy,:),4),
     $      REAL(sti(ix,iy,:),4) /)
        ENDDO
        WRITE(xy_unit)
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CALL close_bin(xy_unit,TRIM(filename))
      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE xy_slice
c-----------------------------------------------------------------------
c     subprogram 4. xt_slice.
c     produces time-series binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE xt_slice(y0fac,istep,t)

      INTEGER(i4), INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t,y0fac

      INTEGER(i4) :: ix,mpx,iy,mpy
c-----------------------------------------------------------------------
c     write binary data.
c-TMP array constructors are needed to avoid a write error on the c90
c     and j90.
c-----------------------------------------------------------------------
      CALL timer(timestart)
      CALL struct_set

      mpx=SIZE(srz,1)-1
      mpy=SIZE(srz,2)-1
      iy=NINT(y0fac*mpy)
      DO ix=0,mpx
        WRITE(xt_unit) (/
     $    REAL(ix,4)/mpx,REAL(istep,4),REAL(t,4),
     $    REAL(srz(ix,iy,:),4),REAL(sbq(ix,iy,:),4),
     $    REAL(sjq(ix,iy,:),4),REAL(svq(ix,iy,:),4),
     $    REAL(spq(ix,iy,:),4),REAL(speq(ix,iy,:),4),
     $    REAL(sndq(ix,iy,:),4),REAL(sdfs(ix,iy,:),4),
     $    REAL(sbe(ix,iy,:),4),REAL(sja(ix,iy,:),4),
     $    REAL(sve(ix,iy,:),4),REAL(spr(ix,iy,:),4),
     $    REAL(spe(ix,iy,:),4),REAL(snd(ix,iy,:),4),
     $    REAL(sco(ix,iy,:),4),REAL(ste(ix,iy,:),4),
     $    REAL(sti(ix,iy,:),4) /)
      ENDDO
      WRITE(xt_unit)

      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE xt_slice
c-----------------------------------------------------------------------
c     subprogram 5. yt_slice.
c     produces time-series binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE yt_slice(x0fac,istep,t)

      INTEGER(i4), INTENT(IN) :: istep
      REAL(r8), INTENT(IN) :: t,x0fac

      INTEGER(i4) :: ix,mpx,iy,mpy
c-----------------------------------------------------------------------
c     write binary data.
c-TMP array constructors are needed to avoid a write error on the c90
c     and j90.
c-----------------------------------------------------------------------
      CALL timer(timestart)
      CALL struct_set

      mpx=SIZE(srz,1)-1
      mpy=SIZE(srz,2)-1
      ix=NINT(x0fac*mpx)
      DO iy=0,mpy
        WRITE(yt_unit) (/
     $    REAL(iy,4)/mpy,REAL(istep,4),REAL(t,4),
     $    REAL(srz(ix,iy,:),4),REAL(sbq(ix,iy,:),4),
     $    REAL(sjq(ix,iy,:),4),REAL(svq(ix,iy,:),4),
     $    REAL(spq(ix,iy,:),4),REAL(speq(ix,iy,:),4),
     $    REAL(sndq(ix,iy,:),4),REAL(sdfs(ix,iy,:),4),
     $    REAL(sbe(ix,iy,:),4),REAL(sja(ix,iy,:),4),
     $    REAL(sve(ix,iy,:),4),REAL(spr(ix,iy,:),4),
     $    REAL(spe(ix,iy,:),4),REAL(snd(ix,iy,:),4),
     $    REAL(sco(ix,iy,:),4),REAL(ste(ix,iy,:),4),
     $    REAL(sti(ix,iy,:),4) /)
      ENDDO
      WRITE(yt_unit)

      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE yt_slice
c-----------------------------------------------------------------------
c     subprogram 6. time_slice_init
c     initialize the time_slice files.
c-----------------------------------------------------------------------
      SUBROUTINE time_slice_init

c-----------------------------------------------------------------------
c     create file names, and open files.
c-----------------------------------------------------------------------
      CALL timer(timestart)

      xt_file=TRIM(xdraw_dir)//"/xt_slice.bin"
      yt_file=TRIM(xdraw_dir)//"/yt_slice.bin"
      CALL open_bin(xt_unit,TRIM(xt_file),"UNKNOWN","REWIND",32_i4)
      CALL open_bin(yt_unit,TRIM(yt_file),"UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE time_slice_init
c-----------------------------------------------------------------------
c     subprogram 7. time_slice_close
c     initialize the time_slice files.
c-----------------------------------------------------------------------
      SUBROUTINE time_slice_close


      CALL timer(timestart)

      CALL close_bin(xt_unit,TRIM(xt_file))
      CALL close_bin(yt_unit,TRIM(yt_file))
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      CALL timer(timeend)
      time_io = time_io + timeend-timestart

      RETURN
      END SUBROUTINE time_slice_close
c-----------------------------------------------------------------------
c     subprogram 8. struct_set.
c     create arrays containing data for the structured rblock region.
c-----------------------------------------------------------------------
      SUBROUTINE struct_set

      INTEGER(i4) :: ix,iy,mpx,mpy,nm,mx1,mx2,my1,my2,ixbl,iybl,mxb,myb,
     $               pd,ixs,ixe,iys,iye,ntb,toff,iv,im,ibl
      REAL(r8) :: dx,dy,dl
c-----------------------------------------------------------------------
c     determine the resolution of data used for the slice plots,
c     and allocate arrays.
c-----------------------------------------------------------------------
      IF (nrbl<1) RETURN
      pd=poly_degree
      mpx=pd*mx
      mpy=pd*my
      dl=1._r8/pd
      nm=rb(1)%be%nfour
c-----------------------------------------------------------------------
c     offset for a tblock0 pie of tblocks at the magnetic axis.
c-PRE this can be improved later for higher order triangles.
c-----------------------------------------------------------------------
      IF (pieflag=='tblock0'.AND.xmin==0.AND.
     $    (gridshape=='circ'.OR.gridshape=='flux')) THEN
        toff=1
      ELSE
        toff=0
      ENDIF
      ALLOCATE(sbe(0:mpx+toff,0:mpy,6*nm),sja(0:mpx+toff,0:mpy,6*nm),
     $         sve(0:mpx+toff,0:mpy,6*nm),spr(0:mpx+toff,0:mpy,2*nm),
     $         spe(0:mpx+toff,0:mpy,2*nm),snd(0:mpx+toff,0:mpy,2*nm),
     $         sco(0:mpx+toff,0:mpy,2*nm),ste(0:mpx+toff,0:mpy,2*nm),
     $         sti(0:mpx+toff,0:mpy,2*nm))
      ALLOCATE(sbq(0:mpx+toff,0:mpy,3),sjq(0:mpx+toff,0:mpy,3),
     $         svq(0:mpx+toff,0:mpy,3),spq(0:mpx+toff,0:mpy,1),
     $         speq(0:mpx+toff,0:mpy,1),sndq(0:mpx+toff,0:mpy,1),
     $         srz(0:mpx+toff,0:mpy,2),sdfs(0:mpx+toff,0:mpy,1))
c-----------------------------------------------------------------------
c     fill the global arrays with data from each block.  there is
c     redundancy in the grid-vertex assignments.
c-----------------------------------------------------------------------
      ibl=0
      mx2=0
      DO ixbl=1,nxbl
        mx1=mx2
        mx2=mx*ixbl/nxbl
        my2=0
        DO iybl=1,nybl
          my1=my2
          my2=my*iybl/nybl
          ibl=ibl+1
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          DO iy=0,pd
            dy=dl*iy
            iys=pd*my1+iy
            iye=pd*(my2-1)+iy
            DO ix=0,pd
              dx=dl*ix
              ixs=toff+pd*mx1+ix
              ixe=toff+pd*(mx2-1)+ix
              CALL struct_set_laq2(rb(ibl)%rz,srz)
              CALL struct_set_laq2(rb(ibl)%be_eq,sbq)
              CALL struct_set_laq2(rb(ibl)%ja_eq,sjq)
              CALL struct_set_laq2(rb(ibl)%ve_eq,svq)
              CALL struct_set_laq2(rb(ibl)%pres_eq,spq)
              CALL struct_set_laq2(rb(ibl)%prese_eq,speq)
              CALL struct_set_laq2(rb(ibl)%nd_eq,sndq)
              CALL struct_set_laq2(rb(ibl)%diff_shape,sdfs)
              CALL struct_set_laq(rb(ibl)%be,sbe)
              CALL struct_set_laq(rb(ibl)%ja,sja)
              CALL struct_set_laq(rb(ibl)%ve,sve)
              CALL struct_set_laq(rb(ibl)%pres,spr)
              CALL struct_set_laq(rb(ibl)%prese,spe)
              CALL struct_set_laq(rb(ibl)%nd,snd)
              CALL struct_set_laq(rb(ibl)%conc,sco)
              CALL struct_set_laq(rb(ibl)%tele,ste)
              CALL struct_set_laq(rb(ibl)%tion,sti)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     fill the space reserved for a tblock0
c-----------------------------------------------------------------------
      IF (toff>0) THEN
        ntb=nrbl+1
        srz(0,:,1)=tb(ntb)%tgeom%xs(0)
        srz(0,:,2)=tb(ntb)%tgeom%ys(0)
        spq(0,:,1)=tb(ntb)%pres_eq%fs(1,0,0)
        speq(0,:,1)=tb(ntb)%prese_eq%fs(1,0,0)
        sndq(0,:,1)=tb(ntb)%nd_eq%fs(1,0,0)
        sdfs(0,:,1)=tb(ntb)%diff_shape%fs(1,0,0)
        DO iv=1,3
          sbq(0,:,iv)=tb(ntb)%be_eq%fs(iv,0,0)
          sjq(0,:,iv)=tb(ntb)%ja_eq%fs(iv,0,0)
          svq(0,:,iv)=tb(ntb)%ve_eq%fs(iv,0,0)
        ENDDO
        DO im=1,nm
          DO iv=1,3
            sbe(0,:,iv+6*(im-1))=REAL(tb(ntb)%be%fs(iv,0,0,im))
            sbe(0,:,iv+6*(im-1)+3)=AIMAG(tb(ntb)%be%fs(iv,0,0,im))
            sja(0,:,iv+6*(im-1))=REAL(tb(ntb)%ja%fs(iv,0,0,im))
            sja(0,:,iv+6*(im-1)+3)=AIMAG(tb(ntb)%ja%fs(iv,0,0,im))
            sve(0,:,iv+6*(im-1))=REAL(tb(ntb)%ve%fs(iv,0,0,im))
            sve(0,:,iv+6*(im-1)+3)=AIMAG(tb(ntb)%ve%fs(iv,0,0,im))
          ENDDO
          spr(0,:,1+2*(im-1))=REAL(tb(ntb)%pres%fs(1,0,0,im))
          spr(0,:,2+2*(im-1))=AIMAG(tb(ntb)%pres%fs(1,0,0,im))
          spe(0,:,1+2*(im-1))=REAL(tb(ntb)%prese%fs(1,0,0,im))
          spe(0,:,2+2*(im-1))=AIMAG(tb(ntb)%prese%fs(1,0,0,im))
          snd(0,:,1+2*(im-1))=REAL(tb(ntb)%nd%fs(1,0,0,im))
          snd(0,:,2+2*(im-1))=AIMAG(tb(ntb)%nd%fs(1,0,0,im))
          sco(0,:,1+2*(im-1))=REAL(tb(ntb)%conc%fs(1,0,0,im))
          sco(0,:,2+2*(im-1))=AIMAG(tb(ntb)%conc%fs(1,0,0,im))
          ste(0,:,1+2*(im-1))=REAL(tb(ntb)%tele%fs(1,0,0,im))
          ste(0,:,2+2*(im-1))=AIMAG(tb(ntb)%tele%fs(1,0,0,im))
          sti(0,:,1+2*(im-1))=REAL(tb(ntb)%tion%fs(1,0,0,im))
          sti(0,:,2+2*(im-1))=AIMAG(tb(ntb)%tion%fs(1,0,0,im))
        ENDDO
      ENDIF

      RETURN
c-----------------------------------------------------------------------
c     evaluate a complex lagrange quadrilateral function, and transfer
c     it into real strorage format.
c-----------------------------------------------------------------------
      CONTAINS
        SUBROUTINE struct_set_laq(laq,rarr)

        TYPE(lagr_quad_type), INTENT(IN) :: laq
        REAL(r8), DIMENSION(0:,0:,:), INTENT(OUT) :: rarr

        COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: ctmp
        COMPLEX(r8), DIMENSION(1,1,1,1) :: dc
        INTEGER(i4) :: iv,nv,im,iq

        nv=laq%nqty
        ALLOCATE(ctmp(nv,mxb,myb,nm))
        CALL lagr_quad_all_eval(laq,dx,dy,ctmp,dc,dc,0_i4)
        DO im=1,nm
          DO iv=1,nv
            iq=2*nv*(im-1)+iv
            rarr(ixs:ixe:pd,iys:iye:pd,iq)=REAL(ctmp(iv,:,:,im))
            iq=iq+nv
            rarr(ixs:ixe:pd,iys:iye:pd,iq)=AIMAG(ctmp(iv,:,:,im))
          ENDDO
        ENDDO
        DEALLOCATE(ctmp)

        RETURN
        END SUBROUTINE struct_set_laq


        SUBROUTINE struct_set_laq2(laq,rarr)

        TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
        REAL(r8), DIMENSION(0:,0:,:), INTENT(OUT) :: rarr

        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: tmp
        REAL(r8), DIMENSION(1,1,1) :: dr
        INTEGER(i4) :: iv,nv,im,iq

        nv=laq%nqty
        ALLOCATE(tmp(nv,mxb,myb))
        CALL lagr_quad_all_eval(laq,dx,dy,tmp,dr,dr,0_i4)
        DO iv=1,nv
          rarr(ixs:ixe:pd,iys:iye:pd,iv)=REAL(tmp(iv,:,:))
        ENDDO
        DEALLOCATE(tmp)

        RETURN
        END SUBROUTINE struct_set_laq2

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE struct_set
c-----------------------------------------------------------------------
c     subprogram 9. struct_dealloc.
c     create a cell-centered grid and copies appropriate j and v.
c-----------------------------------------------------------------------
      SUBROUTINE struct_dealloc

      IF (nrbl<1) RETURN
      DEALLOCATE(sbe,sja,sve,spr,spe,snd,sco,ste,sti)
      DEALLOCATE(sbq,sjq,svq,spq,speq,sndq,srz,sdfs)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE struct_dealloc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE DIAGNOSE
