c-----------------------------------------------------------------------
c     file diagnose.f
c     contains xdraw output routines for nimplot, converted from the
c     comparable nimset file.
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
c     8.  struct_set
c     9.  struct_dealloc
c     10. fl_surface.
c-----------------------------------------------------------------------
c     subprogram 0. global declarations.
c     declares everything needed for the whole module.
c-----------------------------------------------------------------------
      MODULE diagnose
      USE local
      USE input
      USE fields
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: sbe,sja,sve,spr,spe,
     $          snd,sco,ste,sti,sbq,sjq,svq,spq,speq,sndq,srz,sw1,sdfs
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: pol_flux
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
 10   FORMAT('(/2x,"iy",',i3.3,'(4x,"ix=",i3.3,1x)/)')
 20   FORMAT('(i4,1p,',i3.3,'e11.3)')
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
      IF(.NOT.detflag)RETURN
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
      RETURN
      END SUBROUTINE detail
c-----------------------------------------------------------------------
c     subprogram 3. xy_slice.
c     produces 2D binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE xy_slice(istep,t,mode_st,mode_en,fldata)

      INTEGER(i4), INTENT(IN) :: istep,mode_st,mode_en
      REAL(r8), INTENT(IN) :: t
      LOGICAL, INTENT(IN) :: fldata

      CHARACTER(5) :: stepname
      CHARACTER(64) :: filename
      INTEGER(i4) :: ix,iy,iq,eq,is,es,mpx,mpy
c-----------------------------------------------------------------------
c     open file and fill pointers.
c-----------------------------------------------------------------------
      filename=TRIM(xdraw_dir)//"/xy_slice"
      IF (fldata) filename=TRIM(filename)//"pf"
      IF (istep>0) THEN
        WRITE(stepname,fmt='(i5.5)')istep
        filename=TRIM(filename)//"."//TRIM(stepname)//".bin"
      ELSE
        filename=TRIM(filename)//".bin"
      ENDIF
      CALL open_bin(xy_unit,TRIM(filename),"UNKNOWN","REWIND",32_i4)
      CALL struct_set
c-----------------------------------------------------------------------
c     set quantity limits.  this is unique to nimplot.
c-----------------------------------------------------------------------
      iq=6*(mode_st-1)+1
      eq=6*(mode_en)
      is=2*(mode_st-1)+1
      es=2*(mode_en)
c-----------------------------------------------------------------------
c     write binary data.  if poloidal flux is written, compute and write
c     it as part of each record, too, then deallocate pol_flux.  this is
c     a nimplot-only option.
c-----------------------------------------------------------------------
      mpx=SIZE(srz,1)-1
      mpy=SIZE(srz,2)-1
      IF (fldata) THEN
        CALL fl_surface("no plot",.true.,t)
        DO iy=0,mpy
          DO ix=0,mpx
            WRITE(xy_unit) (/
     $        REAL(ix,4)/mpx,REAL(iy,4)/mpy,
     $        REAL(srz(ix,iy,:),4),REAL(sbq(ix,iy,:),4),
     $        REAL(sjq(ix,iy,:),4),REAL(svq(ix,iy,:),4),
     $        REAL(spq(ix,iy,:),4),REAL(speq(ix,iy,:),4),
     $        REAL(sndq(ix,iy,:),4),REAL(sdfs(ix,iy,:),4),
     $        REAL(sbe(ix,iy,iq:eq),4),REAL(sja(ix,iy,iq:eq),4),
     $        REAL(sve(ix,iy,iq:eq),4),REAL(spr(ix,iy,is:es),4),
     $        REAL(spe(ix,iy,is:es),4),REAL(snd(ix,iy,is:es),4),
     $        REAL(sco(ix,iy,is:es),4),REAL(ste(ix,iy,is:es),4),
     $        REAL(sti(ix,iy,is:es),4),REAL(pol_flux(ix,iy),4),
     $        REAL(SQRT(pol_flux(ix,iy)),4) /)
     $        
          ENDDO
          WRITE(xy_unit)
        ENDDO
        DEALLOCATE(pol_flux)
      ELSE
        DO iy=0,mpy
          DO ix=0,mpx
            WRITE(xy_unit) (/
     $        REAL(ix,4)/mpx,REAL(iy,4)/mpy,
     $        REAL(srz(ix,iy,:),4),REAL(sbq(ix,iy,:),4),
     $        REAL(sjq(ix,iy,:),4),REAL(svq(ix,iy,:),4),
     $        REAL(spq(ix,iy,:),4),REAL(speq(ix,iy,:),4),
     $        REAL(sndq(ix,iy,:),4),REAL(sdfs(ix,iy,:),4),
     $        REAL(sbe(ix,iy,iq:eq),4),REAL(sja(ix,iy,iq:eq),4),
     $        REAL(sve(ix,iy,iq:eq),4),REAL(spr(ix,iy,is:es),4),
     $        REAL(spe(ix,iy,is:es),4),REAL(snd(ix,iy,is:es),4),
     $        REAL(sco(ix,iy,is:es),4),REAL(ste(ix,iy,is:es),4),
     $        REAL(sti(ix,iy,is:es),4) /)
     $        
          ENDDO
          WRITE(xy_unit)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CALL close_bin(xy_unit,TRIM(filename))
      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE xy_slice
c-----------------------------------------------------------------------
c     subprogram 4. xt_slice.
c     produces time-series binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE xt_slice(y0fac,istep,t,mode_st,mode_en)

      INTEGER(i4), INTENT(IN) :: istep,mode_st,mode_en
      REAL(r8), INTENT(IN) :: t,y0fac

      INTEGER(i4) :: ix,iy,iq,eq,is,es,mpx,mpy
c-----------------------------------------------------------------------
c     set quantity limits.  this is unique to nimplot.
c-----------------------------------------------------------------------
      iq=6*(mode_st-1)+1
      eq=6*(mode_en)
      is=2*(mode_st-1)+1
      es=2*(mode_en)
c-----------------------------------------------------------------------
c     write binary data.
c-TMP array constructors are needed to avoid a write error on the c90
c     and j90.
c-----------------------------------------------------------------------
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
     $    REAL(sbe(ix,iy,iq:eq),4),REAL(sja(ix,iy,iq:eq),4),
     $    REAL(sve(ix,iy,iq:eq),4),REAL(spr(ix,iy,is:es),4),
     $    REAL(spe(ix,iy,is:es),4),REAL(snd(ix,iy,is:es),4),
     $    REAL(sco(ix,iy,is:es),4),REAL(ste(ix,iy,is:es),4),
     $    REAL(sti(ix,iy,is:es),4) /)
      ENDDO
      WRITE(xt_unit)

      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE xt_slice
c-----------------------------------------------------------------------
c     subprogram 5. yt_slice.
c     produces time-series binary output for xdraw.
c-----------------------------------------------------------------------
      SUBROUTINE yt_slice(x0fac,istep,t,mode_st,mode_en)

      INTEGER(i4), INTENT(IN) :: istep,mode_st,mode_en
      REAL(r8), INTENT(IN) :: t,x0fac

      INTEGER(i4) :: ix,iy,iq,eq,is,es,mpx,mpy
c-----------------------------------------------------------------------
c     set quantity limits.  this is unique to nimplot.
c-----------------------------------------------------------------------
      iq=6*(mode_st-1)+1
      eq=6*(mode_en)
      is=2*(mode_st-1)+1
      es=2*(mode_en)
c-----------------------------------------------------------------------
c     write binary data.
c-TMP array constructors are needed to avoid a write error on the c90
c     and j90.
c-----------------------------------------------------------------------
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
     $    REAL(sbe(ix,iy,iq:eq),4),REAL(sja(ix,iy,iq:eq),4),
     $    REAL(sve(ix,iy,iq:eq),4),REAL(spr(ix,iy,is:es),4),
     $    REAL(spe(ix,iy,is:es),4),REAL(snd(ix,iy,is:es),4),
     $    REAL(sco(ix,iy,is:es),4),REAL(ste(ix,iy,is:es),4),
     $    REAL(sti(ix,iy,is:es),4) /)
      ENDDO
      WRITE(yt_unit)

      CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
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
      xt_file=TRIM(xdraw_dir)//"/xt_slice.bin"
      yt_file=TRIM(xdraw_dir)//"/yt_slice.bin"
      CALL open_bin(xt_unit,TRIM(xt_file),"UNKNOWN","REWIND",32_i4)
      CALL open_bin(yt_unit,TRIM(yt_file),"UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE time_slice_init
c-----------------------------------------------------------------------
c     subprogram 7. time_slice_close
c     initialize the time_slice files.
c-----------------------------------------------------------------------
      SUBROUTINE time_slice_close


      CALL close_bin(xt_unit,TRIM(xt_file))
      CALL close_bin(yt_unit,TRIM(yt_file))
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE time_slice_close
c-----------------------------------------------------------------------
c     subprogram 8. struct_set.
c     create arrays containing data for the structured rblock region.
c-----------------------------------------------------------------------
      SUBROUTINE struct_set

      INTEGER(i4) :: ix,iy,mpx,mpy,nm,mx1,mx2,my1,my2,ixbl,iybl,mxb,myb,
     $               pd,ixs,ixe,iys,iye,ntb,toff,iv,im,ibl
      REAL(r8), DIMENSION(1,1,1) :: db
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
     $         sco(0:mpx+toff,0:mpy,2*nm),sw1(0:mpx+toff,0:mpy,6*nm),
     $         ste(0:mpx+toff,0:mpy,2*nm),sti(0:mpx+toff,0:mpy,2*nm))
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
              CALL struct_set_laq(rb(ibl)%be,sbe)
              CALL struct_set_laq(rb(ibl)%ja,sja)
              CALL struct_set_laq(rb(ibl)%ve,sve)
              CALL struct_set_laq(rb(ibl)%pres,spr)
              CALL struct_set_laq(rb(ibl)%prese,spe)
              CALL struct_set_laq(rb(ibl)%nd,snd)
              CALL struct_set_laq(rb(ibl)%conc,sco)
              CALL struct_set_laq(rb(ibl)%work1,sw1)
              CALL struct_set_laq(rb(ibl)%tele,ste)
              CALL struct_set_laq(rb(ibl)%tion,sti)
              CALL struct_set_laq2(rb(ibl)%rz,srz)
              CALL struct_set_laq2(rb(ibl)%be_eq,sbq)
              CALL struct_set_laq2(rb(ibl)%ja_eq,sjq)
              CALL struct_set_laq2(rb(ibl)%ve_eq,svq)
              CALL struct_set_laq2(rb(ibl)%pres_eq,spq)
              CALL struct_set_laq2(rb(ibl)%prese_eq,speq)
              CALL struct_set_laq2(rb(ibl)%nd_eq,sndq)
              CALL struct_set_laq2(rb(ibl)%diff_shape,sdfs)
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
            sw1(0,:,iv+6*(im-1))=REAL(tb(ntb)%work1%fs(iv,0,0,im))
            sw1(0,:,iv+6*(im-1)+3)=AIMAG(tb(ntb)%work1%fs(iv,0,0,im))
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
      DEALLOCATE(sbe,sja,sve,spr,spe,snd,sco,sw1,ste,sti)
      DEALLOCATE(sbq,sjq,svq,spq,speq,sndq,srz,sdfs)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE struct_dealloc
c-----------------------------------------------------------------------
c     subprogram 11. fl_surface.
c     produce limited flux-surface information based on equilibrium
c     and n=0 fields.  this is only for the polar grid region.
c-----------------------------------------------------------------------
      SUBROUTINE fl_surface(plot_type,first,t)

      CHARACTER(*), INTENT(IN) :: plot_type
      LOGICAL, INTENT(IN) :: first
      REAL(r8), INTENT(IN) :: t

      INTEGER(i4) :: ix,iy,iyc,ixc,ntb,isst,isen,ixst,mx,my,iyst
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: be_arr,ja_arr
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: torarea,cellvol
      REAL(r8), DIMENSION(:), ALLOCATABLE :: tor_surf,pol_surf,q_surf,
     $          par_surf,dvol_surf
      REAL(r8), DIMENSION(3) :: be,xext,ext_prod,dbdx,dbdy
      REAL(r8), DIMENSION(2) :: polarea
      REAL(r8), DIMENSION(4,2) :: dl
      REAL(r8) :: dflux,max_flux,min_flux,tor_flux,
     $            polp,polm,tor_den,lo,hi,drdx,dzdx,drdy,dzdy,jac,
     $            dxdr,dxdz,dydr,dydz,dvol,mu_den,bigr
      LOGICAL :: str_needed
c-----------------------------------------------------------------------
c     fill global arrays if needed.
c-----------------------------------------------------------------------
      IF (ALLOCATED(sbe)) THEN
        str_needed=.false.
      ELSE
        str_needed=.true.
        CALL struct_set
      ENDIF
c-----------------------------------------------------------------------
c     create global current and magnetic field arrays--now copied from
c     the structured data arrays to take advantage of higher-order
c     elements when they are used. the phi component of the be_arr array
c     is covariant.
c
c     if plot_type is "no plot", the n=0 magnetic field may
c     be converted to flux components.  in that case, work1 has the
c     r, z components of the solution.
c-----------------------------------------------------------------------
      ixst=0
      mx=SIZE(srz,1)-1
      my=SIZE(srz,2)-1
      IF (gridshape=='circ'.OR.gridshape=='flux') THEN
        iyst=1
      ELSE
        iyst=0
      ENDIF
      ALLOCATE(be_arr(ixst:mx,0:my,3))
      ALLOCATE(torarea(ixst+1:mx,1:my),cellvol(ixst+1:mx,1:my))
      be_arr=0
      IF (nonlinear) THEN
        IF (plot_type=="no plot") THEN
          be_arr=sw1(:,:,1:3)
        ELSE
          be_arr=sbe(:,:,1:3)
        ENDIF
      ENDIF
      IF (geom=='tor') be_arr(:,:,3)=be_arr(:,:,3)*srz(:,:,1)
      be_arr=be_arr+sbq
      DO iyc=1,my
        iy=iyc-1
        DO ixc=ixst+1,mx
          ix=ixc-1
          dl(1,:)=srz(ix+1,iy  ,:)-srz(ix  ,iy  ,:)
          dl(2,:)=srz(ix+1,iy+1,:)-srz(ix+1,iy  ,:)
          dl(3,:)=srz(ix  ,iy+1,:)-srz(ix+1,iy+1,:)
          dl(4,:)=srz(ix  ,iy  ,:)-srz(ix  ,iy+1,:)
          torarea(ixc,iyc)=
     $      0.5*(ABS(dl(2,1)*dl(1,2)-dl(2,2)*dl(1,1))
     $          +ABS(dl(4,1)*dl(3,2)-dl(4,2)*dl(3,1)))
          IF (geom=='tor') THEN
            cellvol(ixc,iyc)=4*twopi*torarea(ixc,iyc)
     $        *(srz(ixc,iyc  ,1)+srz(ixc-1,iyc  ,1)
     $         +srz(ixc,iyc-1,1)+srz(ixc-1,iyc-1,1))
          ELSE
            cellvol(ixc,iyc)=per_length*torarea(ixc,iyc)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     mu0 * current density for parallel current.
c-----------------------------------------------------------------------
      ALLOCATE(ja_arr(ixst+1:mx,my,3))
      DO iyc=1,my
        DO ixc=ixst+1,mx
          drdx=(srz(ixc,iyc  ,1)-srz(ixc-1,iyc  ,1)
     $         +srz(ixc,iyc-1,1)-srz(ixc-1,iyc-1,1))/2
          dzdx=(srz(ixc,iyc  ,2)-srz(ixc-1,iyc  ,2)
     $         +srz(ixc,iyc-1,2)-srz(ixc-1,iyc-1,2))/2
          drdy=(srz(ixc  ,iyc,1)-srz(ixc  ,iyc-1,1)
     $         +srz(ixc-1,iyc,1)-srz(ixc-1,iyc-1,1))/2
          dzdy=(srz(ixc  ,iyc,2)-srz(ixc  ,iyc-1,2)
     $         +srz(ixc-1,iyc,2)-srz(ixc-1,iyc-1,2))/2
          jac=drdx*dzdy-drdy*dzdx
          dxdr= dzdy/jac
          dxdz=-drdy/jac
          dydr=-dzdx/jac
          dydz= drdx/jac
          dbdx=(be_arr(ixc,iyc  ,:)-be_arr(ixc-1,iyc  ,:)
     $         +be_arr(ixc,iyc-1,:)-be_arr(ixc-1,iyc-1,:))/2
          dbdy=(be_arr(ixc  ,iyc,:)-be_arr(ixc  ,iyc-1,:)
     $         +be_arr(ixc-1,iyc,:)-be_arr(ixc-1,iyc-1,:))/2
          bigr=(srz(ixc,iyc  ,1)+srz(ixc-1,iyc  ,1)
     $         +srz(ixc,iyc-1,1)+srz(ixc-1,iyc-1,1))/4
          ja_arr(ixc,iyc,1)= dbdx(3)*dxdz+dbdy(3)*dydz
          ja_arr(ixc,iyc,2)=-dbdx(3)*dxdr+dbdy(3)*dydr
          ja_arr(ixc,iyc,3)= dbdx(2)*dxdr+dbdy(2)*dydr
     $                      -dbdx(1)*dxdz-dbdy(1)*dydz
          IF (geom=='tor') ja_arr(ixc,iyc,1:2)=ja_arr(ixc,iyc,1:2)/bigr
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     allocate surface arrays.
c-----------------------------------------------------------------------
      ALLOCATE(pol_flux(ixst:mx,0:my),pol_surf(ixst:mx),
     $         q_surf(ixst:mx),tor_surf(ixst+1:mx),dvol_surf(ixst+1:mx),
     $         par_surf(ixst+1:mx))
c-----------------------------------------------------------------------
c     compute poloidal flux from global vertex (mx,1).
c-----------------------------------------------------------------------
      pol_flux=0
      IF (plot_type=="no plot") THEN
c lower left, working up
        DO ix=ixst,mx
          IF (ix>ixst) THEN
            polarea=(/srz(ix  ,iyst,2)-srz(ix-1,iyst,2),
     $                srz(ix-1,iyst,1)-srz(ix  ,iyst,1)/)
            be=(be_arr(ix,iyst,:)+be_arr(ix-1,iyst,:))/2
            IF (geom=='tor') THEN
              polarea=pi*polarea*(srz(ix,iyst,1)+srz(ix-1,iyst,1))
            ELSE
              polarea=per_length*polarea
            ENDIF
            pol_flux(ix,iyst)=pol_flux(ix-1,iyst)+SUM(be(1:2)*polarea)
          ENDIF
          DO iy=iyst+1,my
            polarea=(/srz(ix,iy,2)-srz(ix,iy-1,2),
     $                srz(ix,iy-1,1)-srz(ix,iy,1)/)
            be=(be_arr(ix,iy,:)+be_arr(ix,iy-1,:))/2
            IF (geom=='tor') THEN
              polarea=pi*polarea*(srz(ix,iy-1,1)+srz(ix,iy,1))
            ELSE
              polarea=per_length*polarea
            ENDIF
            pol_flux(ix,iy)=pol_flux(ix,iy-1)+SUM(be(1:2)*polarea)
          ENDDO
        ENDDO
      ELSE
c from lower right
        DO iy=iyst,my
          IF (iy>iyst) THEN
            polarea=(/srz(mx,iy  ,2)-srz(mx,iy-1,2),
     $                srz(mx,iy-1,1)-srz(mx,iy  ,1)/)
            be=(be_arr(mx,iy,:)+be_arr(mx,iy-1,:))/2
            IF (geom=='tor') THEN
              polarea=pi*polarea*(srz(mx,iy,1)+srz(mx,iy-1,1))
            ELSE
              polarea=per_length*polarea
            ENDIF
            pol_flux(mx,iy)=pol_flux(mx,iy-1)+SUM(be(1:2)*polarea)
          ENDIF
          DO ix=mx-1,ixst,-1
            polarea=(/srz(ix  ,iy,2)-srz(ix+1,iy,2),
     $                srz(ix+1,iy,1)-srz(ix  ,iy,1)/)
            be=(be_arr(ix,iy,:)+be_arr(ix+1,iy,:))/2
            IF (geom=='tor') THEN
              polarea=pi*polarea*(srz(ix+1,iy,1)+srz(ix,iy,1))
            ELSE
              polarea=per_length*polarea
            ENDIF
            pol_flux(ix,iy)=pol_flux(ix+1,iy)+SUM(be(1:2)*polarea)
          ENDDO
        ENDDO
      ENDIF
c from top left
c     DO iy=my,iyst,-1
c       IF (iy<my) THEN
c         polarea=(/srz(ixst,iy  ,2)-srz(ixst,iy+1,2),
c    $              srz(ixst,iy+1,1)-srz(ixst,iy  ,1)/)
c         be=(be_arr(ixst,iy,:)+be_arr(ixst,iy+1,:))/2
c         IF (geom=='tor') THEN
c           polarea=pi*polarea*(srz(ixst,iy,1)+srz(ixst,iy+1,1))
c         ELSE
c           polarea=per_length*polarea
c         ENDIF
c         pol_flux(ixst,iy)=pol_flux(ixst,iy+1)+SUM(be(1:2)*polarea)
c       ENDIF
c       DO ix=ixst+1,mx
c         polarea=(/srz(ix  ,iy,2)-srz(ix-1,iy,2),
c    $              srz(ix-1,iy,1)-srz(ix  ,iy,1)/)
c         be=(be_arr(ix,iy,:)+be_arr(ix-1,iy,:))/2
c         IF (geom=='tor') THEN
c           polarea=pi*polarea*(srz(ix-1,iy,1)+srz(ix,iy,1))
c         ELSE
c           polarea=per_length*polarea
c         ENDIF
c         pol_flux(ix,iy)=pol_flux(ix-1,iy)+SUM(be(1:2)*polarea)
c       ENDDO
c     ENDDO
      IF (iyst==1) THEN
        pol_flux(:,0)=pol_flux(:,my)
        pol_flux(ixst,:)=0.5*(MAXVAL(pol_flux(ixst,:))
     $                       +MINVAL(pol_flux(ixst,:)))
      ENDIf
c-----------------------------------------------------------------------
c     choose the type of computation/plot that's appropriate.
c-----------------------------------------------------------------------
      SELECT CASE (plot_type)
c-----------------------------------------------------------------------
c     just compute poloidal flux for use elsewhere; the pol_flux array
c     must then be deallocated elsewhere.  reset values to increase from
c     the grid axis and re-normalize.
c-----------------------------------------------------------------------
      CASE("no plot")
        IF (SUM(pol_flux(ixst,:))>SUM(pol_flux(mx,:)))pol_flux=-pol_flux
        pol_flux=(pol_flux-MINVAL(pol_flux))/
     $           (MAXVAL(pol_flux)-MINVAL(pol_flux))
c-----------------------------------------------------------------------
c     write contour plots of poloidal flux and r*bphi.
c-----------------------------------------------------------------------
      CASE ("polflux")
        IF (first) THEN
          CALL open_bin(con_unit,"polflux.bin","UNKNOWN","REWIND",32_i4)
          WRITE(con_unit) 1_4,0_4,2_4
          WRITE(con_unit) INT(mx-ixst,4),INT(my,4)
          WRITE(con_unit) REAL(srz,4)
        ELSE
          CALL open_bin(con_unit,"polflux.bin","OLD","APPEND",32_i4)
        ENDIF
        WRITE(con_unit) REAL(pol_flux,4)
        WRITE(con_unit) REAL(be_arr(:,:,3),4)
        CALL close_bin(con_unit,"polflux.bin")
        DEALLOCATE(pol_flux)
c-----------------------------------------------------------------------
c     make bins of flux as a function of poloidal flux.
c-----------------------------------------------------------------------
      CASE ("q")
        max_flux=MAXVAL(pol_flux)
        min_flux=MINVAL(pol_flux)
        IF (pol_flux(ixst,0)>0) THEN
          dflux=(min_flux-max_flux)
          pol_surf(ixst)=max_flux
        ELSE
          dflux=(max_flux-min_flux)
          pol_surf(ixst)=min_flux
        ENDIF
        DO ix=ixst+1,mx
          pol_surf(ix)=pol_surf(ixst)+dflux*REAL(ix-ixst)/REAL(mx-ixst)
        ENDDO
c-----------------------------------------------------------------------
c       accumulate toroidal flux within each bin.  loop over cells, find
c       the toroidal flux in the cell, and associate it with the
c       poloidal flux increment between the max and min values
c       surrounding the cell.
c-----------------------------------------------------------------------
        tor_surf=0
        par_surf=0
        dvol_surf=0
        DO iy=1,my
          DO ix=ixst+1,mx
            be=(be_arr(ix  ,iy  ,:)+be_arr(ix-1,iy  ,:)
     $         +be_arr(ix  ,iy-1,:)+be_arr(ix-1,iy-1,:))/4
            IF (geom=='tor') THEN
              be(3)=4*be(3)/(srz(ix,iy  ,1)+srz(ix-1,iy  ,1)
     $                      +srz(ix,iy-1,1)+srz(ix-1,iy-1,1))
            ENDIF
            tor_flux=be(3)*torarea(ix,iy)
            polm=MINVAL(pol_flux(ix-1:ix,iy-1:iy))
            polp=MAXVAL(pol_flux(ix-1:ix,iy-1:iy))
            tor_den=tor_flux/(polp-polm)
            dvol=cellvol(ix,iy)/(polp-polm)
            mu_den=dvol*SUM(ja_arr(ix,iy,:)*be)/SUM(be**2)
            IF (pol_surf(ixst)>0.AND.polm>=pol_surf(ixst).OR.
     $          pol_surf(ixst)<0.AND.polm<=pol_surf(ixst)) isst=ixst+1
            IF (pol_surf(ixst)>0.AND.polp>=pol_surf(ixst).OR.
     $          pol_surf(ixst)<0.AND.polp<=pol_surf(ixst)) isen=ixst+1
            DO ixc=ixst+1,mx
              IF ((pol_surf(ixc-1)-polm)*(pol_surf(ixc)-polm)<0)isst=ixc
              IF ((pol_surf(ixc-1)-polp)*(pol_surf(ixc)-polp)<0)isen=ixc
            ENDDO
            IF (pol_surf(ixst)>0.AND.polm<=pol_surf(mx).OR.
     $          pol_surf(ixst)<0.AND.polm>=pol_surf(mx)) isst=mx
            IF (pol_surf(ixst)>0.AND.polp<=pol_surf(mx).OR.
     $          pol_surf(ixst)<0.AND.polp>=pol_surf(mx)) isen=mx
            DO ixc=isst,isen,SIGN(1_i4,isen-isst)
              IF (pol_surf(ixc)<pol_surf(ixc-1)) THEN
                lo=MAX(polm,pol_surf(ixc  ))
                hi=MIN(polp,pol_surf(ixc-1))
              ELSE
                lo=MAX(polm,pol_surf(ixc-1))
                hi=MIN(polp,pol_surf(ixc  ))
              ENDIF
              tor_surf(ixc)=tor_surf(ixc)+(hi-lo)*tor_den
              dvol_surf(ixc)=dvol_surf(ixc)+(hi-lo)*dvol
              par_surf(ixc)=par_surf(ixc)+(hi-lo)*mu_den
            ENDDO
          ENDDO
        ENDDO
        par_surf=par_surf/dvol_surf
c-----------------------------------------------------------------------
c       reset pol_surf to increase from the magnetic axis, then
c       accumulate net toroidal flux from magnetic axis, then find
c       q=d(tor_flux)/d(pol_flux).
c-----------------------------------------------------------------------
        IF (pol_surf(ixst)>0) pol_surf=-pol_surf
        pol_surf=pol_surf-pol_surf(ixst)
        DO ixc=ixst+2,mx
          tor_surf(ixc)=tor_surf(ixc)+tor_surf(ixc-1)
        ENDDO
        dflux=max_flux-min_flux
        DO ix=ixst+1,mx-1
          q_surf(ix)=2*(tor_surf(ix+1)-tor_surf(ix))
     $                /(pol_surf(ix+1)-pol_surf(ix-1))
        ENDDO
        xext=pol_surf(ixst+1:ixst+3)
        ext_prod(1)=xext(2)*xext(3)/
     $    ((xext(1)-xext(2))*(xext(1)-xext(3)))
        ext_prod(2)=xext(1)*xext(3)/
     $    ((xext(2)-xext(1))*(xext(2)-xext(3)))
        ext_prod(3)=xext(1)*xext(2)/
     $    ((xext(3)-xext(1))*(xext(3)-xext(2)))
        q_surf(ixst)=SUM(q_surf(ixst+1:ixst+3)*ext_prod)
        IF (first) THEN
          CALL open_bin(xy_unit,"flsurf.bin","UNKNOWN","REWIND",32_i4)
        ELSE
          CALL open_bin(xy_unit,"flsurf.bin","OLD","APPEND",32_i4)
        ENDIF
        DO ix=ixst,mx-1
          polm=pol_surf(ix)/dflux
          polp=(polm+pol_surf(ix+1)/dflux)/2
          WRITE(xy_unit) REAL(t,4),REAL(ix,4),REAL(polm,4),
     $      REAL(SIGN(1._r8,polm)*SQRT(ABS(polm)),4),REAL(polp,4),
     $      REAL(SIGN(1._r8,polp)*SQRT(ABS(polp)),4),
     $      REAL(q_surf(ix),4),REAL(par_surf(ix+1),4)
        ENDDO
        WRITE(xy_unit)
        CALL close_bin(xy_unit,"flsurf.bin")
        DEALLOCATE(pol_flux)
      CASE DEFAULT
        CALL nim_stop("Unrecognized plot type in fl_surface.")
      END SELECT
c-----------------------------------------------------------------------
c     deallocations.
c-----------------------------------------------------------------------
      DEALLOCATE(be_arr,tor_surf,pol_surf,q_surf,torarea,
     $           ja_arr,par_surf,cellvol,dvol_surf)
      IF (str_needed) CALL struct_dealloc
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fl_surface
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE DIAGNOSE
