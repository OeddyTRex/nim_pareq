c-----------------------------------------------------------------------
c     file physics_init.f:  contains routines for setting the initial
c     conditions on fundmental dependent variables.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. polar_phys_init.
c     2. phys_init.
c     3. eig_phys_init.
c     4. field_err_init.
c-----------------------------------------------------------------------
c     subprogram 1. polar_phys_init.
c     initializes the physical fields on polar-like grids.
c-----------------------------------------------------------------------
      SUBROUTINE polar_phys_init(nmodes,keff)
      USE local
      USE input
      USE physdat
      USE fields
      USE lagr_quad_mod
      USE tri_linear
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

      INTEGER(i4) :: ix,iy,ibl,iv,ibv,nv,jv,nlim,ivstep,iqty,
     $               ix0,ix1,iy0,iy1,mm,i,j,nbl_init,iny,nyend,ibasis,
     $               max_basisr,max_basist
      REAL(r8) :: xfac,xfac1,yfac,yfac2,yfac1,rfac,bex,bey,bez,omega,
     $            rv,zv,tv,rad,rad_out,dth,dth_edge,vamp,bigr,dx,dy
      REAL(r8), DIMENSION(nmodes) :: kk
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ro_edge,th_edge
      REAL(r8), DIMENSION(nx) :: kr,acoe,bcoe
      REAL(r8) :: bessel_j,bessel_y
      COMPLEX(r8), DIMENSION(nmodes) :: vr,vt,vz
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: temp
      CHARACTER(80) :: msg
c-----------------------------------------------------------------------
c     define format for error message.
c-----------------------------------------------------------------------
  100 FORMAT ('polar_phys_init index too ',a,' ibl=',i3,' ix=',i3,
     $        ' iy=',i3,' iv=',i3)
c-----------------------------------------------------------------------
c     check requested wavenumber and initialization type.
c-----------------------------------------------------------------------
      IF (MODULO(ny,2_i4)/=0)
     $  CALL nim_stop('ny must be even for polar grids.')
      IF (init_type(1:8)=='linear b' .AND. geom/='lin') CALL nim_stop
     $  ('The init_type linear b is only suitable for linear geometry.')
c-----------------------------------------------------------------------
c     create lists of edge radius and angle vs. boundary vertex
c     from the outer-most rblocks in the polar grid grid.
c-----------------------------------------------------------------------
      nv=0
      DO i=1,nybl
        ibl=(nxbl-1)*nybl+i
        nv=nv+rb(ibl)%my
      ENDDO
      ALLOCATE(ro_edge(nv+1))
      ALLOCATE(th_edge(nv+1))
      ibl=(nxbl-1)*nybl+1
      ix=rb(ibl)%mx
      iy=1
      rv=rb(ibl)%rz%fs(1,ix,iy)
      zv=rb(ibl)%rz%fs(2,ix,iy)
      xfac=rv-xo
      yfac=zv-yo
      omega=pi+ATAN2(yfac,xfac)
      iv=0
      DO i=1,nybl
        ibl=(nxbl-1)*nybl+i
        DO j=1,rb(ibl)%my
          iv=iv+1
          rv=rb(ibl)%rz%fs(1,rb(ibl)%mx,j)
          zv=rb(ibl)%rz%fs(2,rb(ibl)%mx,j)
          xfac=rv-xo
          yfac=zv-yo
          ro_edge(iv)=SQRT(xfac**2+yfac**2)
          rv= xfac*COS(omega)+yfac*SIN(omega)
          zv=-xfac*SIN(omega)+yfac*COS(omega)
          th_edge(iv)=ATAN2(zv,rv)
        ENDDO
      ENDDO
      ro_edge(nv+1)=ro_edge(1)
      th_edge(nv+1)=th_edge(1)
      dth=th_edge(3)-th_edge(2)
      IF (dth>0..OR.dth<-pi) THEN      !       ccw sequencing
        th_edge(1)=-pi
        th_edge(nv+1)=pi
      ELSE 				     !        cw sequencing
        th_edge(1)=pi
        th_edge(nv+1)=-pi
      ENDIF
c-----------------------------------------------------------------------
c     initialize perturbations
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        rb(ibl)%be=0
        rb(ibl)%ve=0
        rb(ibl)%pres=0
        rb(ibl)%prese=0
        rb(ibl)%tele=0
        rb(ibl)%tion=0
        rb(ibl)%nd=0
        rb(ibl)%conc=0
      ENDDO
      DO ibl=nrbl+1,nbl
        tb(ibl)%be=0
        tb(ibl)%ve=0
        tb(ibl)%pres=0
        tb(ibl)%prese=0
        tb(ibl)%tele=0
        tb(ibl)%tion=0
        tb(ibl)%nd=0
        tb(ibl)%conc=0
      ENDDO
c-----------------------------------------------------------------------
c     initialize waves in the non-rim blocks only.
c-----------------------------------------------------------------------
      vamp=bamp/SQRT(mu0*mtot*ndens)
      IF (TRIM(pieflag)=='rblock') THEN
        nbl_init=nrbl
      ELSE
        nbl_init=nrbl+1
      ENDIF
c-----------------------------------------------------------------------
c     set limits on the basis type index for different block types.
c-----------------------------------------------------------------------
      max_basisr=poly_degree**2
c-PRE
      max_basist=1
c-----------------------------------------------------------------------
c     loop over basis types and blocks.
c-----------------------------------------------------------------------
      ibl=1
      ibasis=1
      bb_loop: DO
        IF (ibl<=nrbl) THEN
          IF (ibasis>max_basisr) THEN
            ibasis=1
            ibl=ibl+1
          ENDIF
        ELSE 
          IF (ibasis>max_basist) THEN
            ibasis=1
            ibl=ibl+1
          ENDIF
        ENDIF
        IF (ibl>nbl_init) EXIT
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       quadrilateral elements.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          ix1=rb(ibl)%mx
          iy1=rb(ibl)%my
          ix0=rb(ibl)%be%ix0(ibasis)
          iy0=rb(ibl)%be%iy0(ibasis)
          dx=rb(ibl)%be%dx(ibasis)
          dy=rb(ibl)%be%dy(ibasis)
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       triangular elements.
c-PRE   only linear implemented so far.
c-----------------------------------------------------------------------
        ELSE
          ix0=0
          ix1=tb(ibl)%mvert
          iy0=0
          iy1=0
        ENDIF
c-----------------------------------------------------------------------
c       loop over elements for a given block and basis.
c-----------------------------------------------------------------------
        ALLOCATE(temp(3,ix0:ix1,iy0:iy1,nmodes))
        temp=0
        DO ix=ix0,ix1
          DO iy=iy0,iy1
            IF (ibl<=nrbl) THEN
              CALL lagr_quad_eval(rb(ibl)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
              rv=rb(ibl)%rz%f(1)
              zv=rb(ibl)%rz%f(2)
            ELSE
              rv=tb(ibl)%tgeom%xs(ix)
              zv=tb(ibl)%tgeom%ys(ix)
            ENDIF
            IF (geom=='tor') THEN
              bigr=rv
            ELSE
              bigr=1
            ENDIF
            xfac=rv-xo
            yfac=zv-yo
            rv= xfac*COS(omega)+yfac*SIN(omega)
            zv=-xfac*SIN(omega)+yfac*COS(omega)
            IF (rv==0.AND.zv==0) rv=1
            tv=ATAN2(zv,rv)
            dth_edge=th_edge(2)-th_edge(1)
            IF (dth_edge>0.) THEN
              IF (tv>=pi) tv=-pi
              jv=nv*((tv+pi)/twopi)+1
            ELSE
              IF (tv<=-pi) tv=pi
              jv=nv*((pi-tv)/twopi)+1
            ENDIF
            jv=MAX(1_i4,jv)
            jv=MIN(nv,jv)
            dth=tv-th_edge(jv)
            IF (dth*dth_edge>0.) THEN
              nlim=nv
              ivstep=1
            ELSE
              nlim=1
              ivstep=-1
            ENDIF
            DO iv=jv,nlim,ivstep
              dth_edge=th_edge(iv+1)-th_edge(iv)
              dth=tv-th_edge(iv)
              IF (ABS(dth)<=ABS(dth_edge).AND.dth*dth_edge>=0.) EXIT
            ENDDO
            IF (iv<1) THEN
              WRITE(msg,100) 'small',ibl,ix,iy,iv
              CALL nim_stop(msg)
            ELSEIF (iv>nv) THEN
              WRITE(msg,100) 'large',ibl,ix,iy,iv
              CALL nim_stop(msg)
            ENDIF
            rad_out=dth/dth_edge*ro_edge(iv+1)
     $             +(1.-dth/dth_edge)*ro_edge(iv)
c-----------------------------------------------------------------------
c           use a loop over m to perturb different components if
c           desired.  then the input, ny is just a limit.
c-----------------------------------------------------------------------
            iny=LEN_TRIM(init_type)
            IF (init_type(iny-3:iny)=='mult') THEN
              nyend=-ny
            ELSE
              nyend=ny
            ENDIF
            m_loop: DO iny=MIN(ny,nyend),MAX(ny,nyend),2
              xfac1=xfac
              vr=0
              vt=0
              vz=0
              IF (init_type(1:9)=='shear alf' .OR.
     $            init_type(1:8)=='linear b') THEN
c-----------------------------------------------------------------------
c               assume k>=0 with exp[i*k*z-i*m*theta]
c-----------------------------------------------------------------------
                mm=iny/2
                kk=rad_out*keff/bigr
                rad=SQRT(xfac**2+yfac**2)/rad_out
                IF (rad==0) THEN
                  tv=0
                ELSE
                  tv=ATAN2(yfac,xfac)
                ENDIF
                yfac1=COS(tv*mm)
                yfac2=SIN(tv*mm)
                IF (mm==0) THEN
                  vr=-kk*rad*(1-rad**2)
                  vt= 2*rad
                  vz= 2-4*rad**2
                ELSE
                  IF (ABS(mm)/=1.OR.rad>0) THEN
                    vr=-rad**(ABS(mm)-1)*(kk+mm)*(1-rad**2)
                    vt= 2*rad**(ABS(mm)+1)
     $                -SIGN(1_i4,mm)*rad**(ABS(mm)-1)*(kk+mm)*(1-rad**2)
                    vz=-2*rad**ABS(mm)
                  ELSE
                    vr=-(kk+mm)*(1-rad**2)
                    vt= 2*rad**(ABS(mm)+1)
     $                 -SIGN(1_i4,mm)*(kk+mm)*(1-rad**2)
                    vz= 0
                  ENDIF
                ENDIF
                vr=vr*(yfac2+(0,1)*yfac1)
                vt=vt*(yfac1-(0,1)*yfac2)
                vz=vz*(yfac1-(0,1)*yfac2)
                DO iqty=1,nmodes
                  IF (keff(iqty)==0) THEN
                    vr(iqty)=REAL(vr(iqty),r8)
                    vt(iqty)=REAL(vt(iqty),r8)
                    vz(iqty)=REAL(vz(iqty),r8)
                  ENDIF
                ENDDO
              ELSE
                yfac1=yfac
                IF (xfac1==0.AND.yfac1==0) THEN
                  tv=0
                ELSE
                  tv=ATAN2(yfac1,xfac1)
                  rad=rad_out
                  IF (gridshape=='circ') THEN
                    CALL annular_getco(1_i4,nx,xmin,rad,kr,acoe,bcoe)
                    rad=SQRT(xfac1**2+yfac1**2)*kr(nx)
                    xfac1=acoe(nx)*bessel_j(1_i4,rad)
                    IF (xmin/=0) xfac1=xfac1+bcoe(nx)*bessel_y(1_i4,rad)
                  ELSE
                    CALL annular_getco(1_i4,nx,0._r8,rad,kr,acoe,bcoe)
                    rad=SQRT(xfac1**2+yfac1**2)*kr(nx)
                    xfac1=bessel_j(1_i4,rad)
                  ENDIF
                ENDIF
                yfac1=COS(tv*iny/2._r8)
                yfac2=SIN(tv*iny/2._r8)
                vr=xfac1*yfac1
                vt=-xfac1*yfac2
                DO iqty=1,nmodes
                  IF (keff(iqty)/=0) THEN
                    vr(iqty)=vr(iqty)-(0,1)*vt(iqty)
                    vt(iqty)=vt(iqty)+(0,1)*REAL(vr(iqty),r8)
                  ELSE
                    vr(iqty)=REAL(vr(iqty),r8)
                    vt(iqty)=REAL(vt(iqty),r8)
                  ENDIF
                  IF (nz>=0.AND.NINT(keff(iqty))/=nz) THEN
                    vr(iqty)=0._r8
                    vt(iqty)=0._r8
                  ENDIF
                ENDDO
                vz=0
              ENDIF
              temp(1,ix,iy,:)=temp(1,ix,iy,:)+vr*COS(tv)-vt*SIN(tv)
              temp(2,ix,iy,:)=temp(2,ix,iy,:)+vr*SIN(tv)+vt*COS(tv)
              temp(3,ix,iy,:)=temp(3,ix,iy,:)+vz
            ENDDO m_loop
          ENDDO
        ENDDO
        IF (init_type(1:8)=='linear b') THEN
          temp=bamp*temp
        ELSE
          temp=vamp*temp
        ENDIF
c-----------------------------------------------------------------------
c       transfer the data to the appropriate basis.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          IF (init_type(1:8)=='linear b') THEN
            CALL lagr_quad_basis_assign_arr(rb(ibl)%be,temp,ibasis)
          ELSE
            CALL lagr_quad_basis_assign_arr(rb(ibl)%ve,temp,ibasis)
          ENDIF
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       triangular elements.
c-PRE   only linear implemented so far.
c-----------------------------------------------------------------------
        ELSE
          IF (init_type(1:8)=='linear b') THEN
            tb(ibl)%be%fs=temp
          ELSE
            tb(ibl)%ve%fs=temp
          ENDIF
        ENDIF
        DEALLOCATE(temp)
        ibasis=ibasis+1
      ENDDO bb_loop 
      DEALLOCATE(ro_edge,th_edge)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_phys_init
c-----------------------------------------------------------------------
c     subprogram 2. phys_init.
c     initializes the physical fields on rectangular domains.
c-----------------------------------------------------------------------
      SUBROUTINE phys_init(nmodes,keff)
      USE local
      USE input
      USE physdat
      USE fields
      USE lagr_quad_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

      INTEGER(i4) :: ix,iy,ibl,mxb,myb,iqty,nyy,nyend,nyinc,ibasis,
     $               ix0,iy0
      REAL(r8) :: kx,ky,xfac,yfac,rfac1,rfac2,bey,bez,
     $            bmag,yfac2,kkx,kky,dx,dy,bex,rc,zc,valf
      REAL(r8) :: rpd,ipd,kyd,rod,iod
      REAL(r8), DIMENSION(3) :: rbd=0._r8,ibd=0._r8,rvd,ivd
      REAL(r8), DIMENSION(nmodes) :: kkz,kkmag,kpermag,kdb,bxw,byw,bzw,
     $            vxw,vyw,vzw,vamp,kperx,kpery,kperz,omega,pelw
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: coskr,sinkr,x,y
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: rz,dumv
      REAL(r8), DIMENSION(nx) :: kr,acoe,bcoe
      REAL(r8) :: bessel_j,bessel_y,rand=0
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: bw,vw,tew,tiw,nw
      CHARACTER(16) :: temp_type
      LOGICAL :: file_stat
c-----------------------------------------------------------------------
c     check init_type.
c-----------------------------------------------------------------------
      IF (init_type(1:8)=='linear b') CALL nim_stop
     $  ('The init_type linear b is not suitable for gridshape=rect.')
c-----------------------------------------------------------------------
c     initialize the perturbations.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        rb(ibl)%be=0
        rb(ibl)%ve=0
        rb(ibl)%pres=0
        rb(ibl)%prese=0
        rb(ibl)%tele=0
        rb(ibl)%tion=0
        rb(ibl)%nd=0
        rb(ibl)%conc=0
      ENDDO
c-----------------------------------------------------------------------
c     allow a loop over ny values to initialize multiple ky's.
c-----------------------------------------------------------------------
      temp_type=ADJUSTR(init_type)
      IF (temp_type(13:16)=='mult') THEN
        nyend=-ny
      ELSE
        nyend=ny
      ENDIF
      IF (periodicity/='none') THEN
        nyinc=2
      ELSE
        nyinc=1
      ENDIF
      nyloop: DO nyy=ny,nyend,SIGN(nyinc,nyend-ny)
c-----------------------------------------------------------------------
c       determine wavenumbers.
c-----------------------------------------------------------------------
        IF (nyy/=ny) CALL RANDOM_NUMBER(rand)
        IF (periodicity/='none') THEN
          IF (MODULO(nyy,2_i4)/=0) THEN
            CALL nim_stop('ny must be even for periodic geometry.')
          ENDIF
        ENDIF
        IF (periodicity=='both') THEN
          IF (MODULO(nx,2_i4)/=0) THEN
            CALL nim_stop('nx must be even for periodic geometry.')
          ENDIF
        ENDIF
        kx=pi/(xmax-xmin)
        ky=pi/(ymax-ymin)
        kkx=nx*kx
        kky=nyy*ky
        kkz=keff
        bex=be0*SIN(pi*thetab)*COS(pi*phib)
        bey=be0*SIN(pi*thetab)*SIN(pi*phib)
        bez=be0*COS(pi*thetab)
        kkmag=SQRT(kkx**2+kky**2+kkz**2)
        kdb=kkx*bex+kky*bey+kkz*bez
c-----------------------------------------------------------------------
c       set Bessel function coefficients for toroidal geometries.
c-----------------------------------------------------------------------
        IF (nx>0.and.geom=='tor') THEN
          CALL annular_getco(1_i4,nx,xmin,xmax,kr,acoe,bcoe)
        ENDIF
c-----------------------------------------------------------------------
c       initialize Alfven and whistler waves.
c-----------------------------------------------------------------------
        bmag=be0
c-----------------------------------------------------------------------
c       shear waves (traveling):
c-----------------------------------------------------------------------
        IF (init_type(1:9)=='shear alf') THEN
          IF (geom=='tor')  CALL nim_stop(
     $      "shear alf requires linear geometry for now.")
          IF (TRIM(periodicity)/='both') CALL nim_stop(
     $      "shear alf presently requires a doubly periodic domain.")
          omega=SQRT(kdb**2/(mu0*mtot*ndens))
          bxw=bamp/(bmag*kkmag)*(bey*kkz-bez*kky)
          byw=bamp/(bmag*kkmag)*(bez*kkx-bex*kkz)
          bzw=bamp/(bmag*kkmag)*(bex*kky-bey*kkx)
          IF (MINVAL(ABS(kdb))<100*TINY(kdb))
     $      CALL nim_stop("Phys_init: Wave vector"
     $                    //" is perpendicular to B0; no shear wave.")
          vamp=omega*bamp/kdb
          vxw=vamp*(bey*kkz-bez*kky)/(kkmag*bmag)
          vyw=vamp*(bez*kkx-bex*kkz)/(kkmag*bmag)
          vzw=vamp*(bex*kky-bey*kkx)/(kkmag*bmag)
c-----------------------------------------------------------------------
c       compressional waves (standing):
c-----------------------------------------------------------------------
        ELSE IF (init_type(1:9)=='compr alf') THEN
          omega=bmag*kkmag/SQRT(mu0*mtot*ndens)
          bxw=0
          byw=0
          bzw=0
          vamp=omega*bamp/(bmag*kkmag)
          DO iqty=1,nmodes
            IF (kdb(iqty)==kkmag(iqty)*bmag) THEN
              vxw(iqty)=vamp(iqty)
              vyw(iqty)=0
              vzw(iqty)=0
            ELSE
              kperx(iqty)=kkx-bex*kdb(iqty)/bmag**2
              kpery(iqty)=kky-bey*kdb(iqty)/bmag**2
              kperz(iqty)=kkz(iqty)-bez*kdb(iqty)/bmag**2
              kpermag(iqty)=
     $           SQRT(kperx(iqty)**2+kpery(iqty)**2+kperz(iqty)**2)
              vxw(iqty)=vamp(iqty)*kperx(iqty)/kpermag(iqty)
              vyw(iqty)=vamp(iqty)*kpery(iqty)/kpermag(iqty)
              vzw(iqty)=vamp(iqty)*kperz(iqty)/kpermag(iqty)
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c       whistler waves (standing):
c-----------------------------------------------------------------------
        ELSE IF (init_type(1:8)=='whistler') THEN
          bxw=bamp*SIN(pi*(thetab+0.5))*COS(pi*phib)
          byw=bamp*SIN(pi*(thetab+0.5))*SIN(pi*phib)
          bzw=bamp*COS(pi*(thetab+0.5))
          vxw=0
          vyw=0
          vzw=0
c-----------------------------------------------------------------------
c       testing g-mode 
c-----------------------------------------------------------------------
        ELSE IF (init_type(1:5)=='gmode') THEN
          omega=bmag*kkmag/SQRT(mu0*mtot*ndens)
          vamp=omega*bamp/(bmag*kkmag)
          bxw=0
          byw=0
          bzw=0
          kperx=0 
          kpery=kky 
          kperz=0
          kpermag=kky
          vxw=vamp*kpery/kpermag
          vyw=0
          vzw=0
c-----------------------------------------------------------------------
c       electron pressure perturbation
c-----------------------------------------------------------------------
        ELSE IF (init_type=='el pres') THEN
          bxw=0
          byw=0
          bzw=0
          vxw=0
          vyw=0
          vzw=0
c-----------------------------------------------------------------------
c       read wave mode from processed output of dispersion code.
c       undo the normalizations used by dispersion after reading mode.
c-----------------------------------------------------------------------
        ELSE IF (init_type=='disp mode') THEN
          INQUIRE(FILE='setdisp.txt',EXIST=file_stat)
          IF (.NOT.file_stat) THEN
            WRITE(nim_wr,*) 'The file setdisp.txt does not exist.'
            STOP
          ENDIF
          OPEN(UNIT=temp_unit,FILE='setdisp.txt',STATUS='OLD')
          READ(temp_unit,*) kyd
          READ(temp_unit,*) rod,iod
          READ(temp_unit,*) rvd(1),ivd(1)
          READ(temp_unit,*) rvd(2),ivd(2)
          READ(temp_unit,*) rvd(3),ivd(3)
          READ(temp_unit,*) rbd(1),ibd(1)
          READ(temp_unit,*) rbd(3),ibd(3)
          READ(temp_unit,*) rpd,ipd
          IF ( ABS((kky-kyd)/kyd)> 1.e-6 )
     $      CALL nim_stop('Phys_init: error matching wavenumber.')
          valf=SQRT(be0**2/(mu0*mtot*ndens))
          rvd=rvd*valf
          ivd=ivd*valf
          rbd=rbd*be0
          ibd=ibd*be0
          rpd=rpd*0.5_r8*beta*be0**2/mu0
          ipd=ipd*0.5_r8*beta*be0**2/mu0
c-----------------------------------------------------------------------
c       vertical push, used for exciting vdes.
c-----------------------------------------------------------------------
        ELSE IF (init_type(1:8)=='vertical') THEN
          omega=bmag*kkmag/SQRT(mu0*mtot*ndens)
          bxw=0
          byw=0
          bzw=0
          vamp=omega*bamp/(bmag*kkmag)
          DO iqty=1,nmodes
            vxw(iqty)=0
            vyw(iqty)=vamp(iqty)
            vzw(iqty)=0
          ENDDO
        ELSE
          CALL nim_stop('Physics_init:  unrecognized init_type.')
        ENDIF
c-----------------------------------------------------------------------
c       if nz is set, only perturb that Fourier component.
c-----------------------------------------------------------------------
        DO iqty=1,nmodes
          IF (nz>=0.AND.NINT(keff(iqty))/=nz) THEN
            bxw(iqty)=0._r8
            byw(iqty)=0._r8
            bzw(iqty)=0._r8
            vxw(iqty)=0._r8
            vyw(iqty)=0._r8
            vzw(iqty)=0._r8
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       apply to all basis functions in all blocks.
c-----------------------------------------------------------------------
        ibl=1
        ibasis=1
        bb_loop: DO
          IF (ibasis>poly_degree**2) THEN
            ibasis=1
            ibl=ibl+1
          ENDIF
          IF (ibl>nbl) EXIT
c-----------------------------------------------------------------------
c         select the offsets for the different basis functions in 
c         quadrilateral elements.
c-----------------------------------------------------------------------
          mxb=rb(ibl)%mx
          myb=rb(ibl)%my
          ix0=rb(ibl)%be%ix0(ibasis)
          iy0=rb(ibl)%be%iy0(ibasis)
          dx=rb(ibl)%be%dx(ibasis)
          dy=rb(ibl)%be%dy(ibasis)
          ALLOCATE(coskr(ix0:mxb,iy0:myb))
          ALLOCATE(sinkr(ix0:mxb,iy0:myb))
          ALLOCATE(x(ix0:mxb,iy0:myb))
          ALLOCATE(y(ix0:mxb,iy0:myb))
          ALLOCATE(bw(3,ix0:mxb,iy0:myb,nmodes))
          ALLOCATE(vw(3,ix0:mxb,iy0:myb,nmodes))
          ALLOCATE(tew(1,ix0:mxb,iy0:myb,nmodes))
          ALLOCATE(tiw(1,ix0:mxb,iy0:myb,nmodes))
          ALLOCATE(nw(1,ix0:mxb,iy0:myb,nmodes))
          bw=0._r8
          vw=0._r8
          tew=0._r8
          tiw=0._r8
          nw=0._r8
c-----------------------------------------------------------------------
c         set the wave shape at the appropriate locations.
c-----------------------------------------------------------------------
          DO ix=ix0,mxb
            DO iy=iy0,myb
              CALL lagr_quad_eval(rb(ibl)%rz,ix-ix0+dx,iy-iy0+dy,0_i4)
              x(ix,iy)=rb(ibl)%rz%f(1)
              y(ix,iy)=rb(ibl)%rz%f(2)
            ENDDO
          ENDDO
          coskr=COS(kky*(y-ymin)-twopi*rand
     $             +kkx*(x-xmin))
          IF (init_type(1:8)=='vertical') THEN
            sinkr=COS(0.5_r8*kky*y-twopi*rand
     $               +kkx*(x-xmin))
          ELSE
            sinkr=SIN(kky*(y-ymin)-twopi*rand
     $               +kkx*(x-xmin))
          ENDIF
c-----------------------------------------------------------------------
c         toroidal initialization (intended for compressional with
c         n=0 only).
c-----------------------------------------------------------------------
          IF (geom=='tor'.AND.nx/=0) THEN
            DO ix=ix0,mxb
              DO iy=iy0,myb
                sinkr(ix,iy)=
     $            acoe(nx)*bessel_j(1_i4,kr(nx)*x(ix,iy))
                IF (xmin/=0) sinkr(ix,iy)=sinkr(ix,iy)
     $           +bcoe(nx)*bessel_y(1_i4,kr(nx)*x(ix,iy))
                IF (.NOT.nonlinear.AND.lin_nmax==0.AND.lin_nmodes==1
     $              .AND.kky/=0.)
     $            sinkr(ix,iy)=sinkr(ix,iy)*SIN(kky*(y(ix,iy)-ymin))
              ENDDO
            ENDDO
            coskr=sinkr  !  avoids cos(r) when used for n>0
          ENDIF
c-----------------------------------------------------------------------
c         compute the b and v perturbations for this block, wave, and
c         wavenumber.
c
c         the initialization is different when reading dispersion
c         output.  here, real and imaginary phases are used in the
c         y-direction, and kz=0 only.
c-----------------------------------------------------------------------
          IF (init_type=='disp mode') THEN
            vw(1,:,:,1)=rvd(1)*coskr-ivd(1)*sinkr
            vw(2,:,:,1)=rvd(2)*coskr-ivd(2)*sinkr
            vw(3,:,:,1)=rvd(3)*coskr-ivd(3)*sinkr
            bw(1,:,:,1)=rbd(1)*coskr-ibd(1)*sinkr
            bw(3,:,:,1)=rbd(3)*coskr-ibd(3)*sinkr
            nw(1,:,:,1)=(rpd*coskr-ipd*sinkr)*
     $                   ndens/(gamma*beta*0.5_r8*be0**2/mu0)
            tew(1,:,:,1)=(rpd*coskr-ipd*sinkr)*gamm1*pe_frac/
     $                   (gamma*kboltz*ndens)
            tiw(1,:,:,1)=(rpd*coskr-ipd*sinkr)*gamm1*(1._r8-pe_frac)/
     $                   (gamma*kboltz*ndens)
          ELSE
            DO iqty=1,nmodes
              bw(1,:,:,iqty)=bxw(iqty)*coskr+(0,1)*bxw(iqty)*sinkr
              bw(2,:,:,iqty)=byw(iqty)*coskr+(0,1)*byw(iqty)*sinkr
              bw(3,:,:,iqty)=bzw(iqty)*coskr+(0,1)*bzw(iqty)*sinkr
              vw(1,:,:,iqty)=vxw(iqty)*sinkr-(0,1)*vxw(iqty)*coskr
              vw(2,:,:,iqty)=vyw(iqty)*sinkr-(0,1)*vyw(iqty)*coskr
              vw(3,:,:,iqty)=vzw(iqty)*sinkr-(0,1)*vzw(iqty)*coskr
              IF (keff(iqty)==0) THEN
                bw(:,:,:,iqty)=REAL(bw(:,:,:,iqty),r8)
                vw(:,:,:,iqty)=REAL(vw(:,:,:,iqty),r8)
              ENDIF
              IF (init_type(1:5)=='gmode') THEN
                vw(1,:,:,iqty)=vxw(iqty)*
     $              SIN((x-xmin)/(xmax-xmin)*pi)*COS(kky*(y-ymin))
                vw(2,:,:,iqty)=vxw(iqty)*pi/(kky*(xmax-xmin))*
     $              COS((x-xmin)/(xmax-xmin)*pi)*SIN(kky*(y-ymin))
              ENDIF
              IF (init_type=='el pres')
     $          tew(1,:,:,iqty)=sinkr*0.5*bamp*beta*be0**2/
     $                          (kboltz*ndens*mu0)
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         transfer to the storage for this block and basis.
c-----------------------------------------------------------------------
          CALL lagr_quad_basis_add_arr(rb(ibl)%be,bw,ibasis)
          CALL lagr_quad_basis_add_arr(rb(ibl)%ve,vw,ibasis)
          CALL lagr_quad_basis_add_arr(rb(ibl)%tele,tew,ibasis)
          CALL lagr_quad_basis_add_arr(rb(ibl)%tion,tiw,ibasis)
          CALL lagr_quad_basis_add_arr(rb(ibl)%nd,nw,ibasis)
          DEALLOCATE(coskr,sinkr,x,y,bw,vw,tew,tiw,nw)
          ibasis=ibasis+1
        ENDDO bb_loop
      ENDDO nyloop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE phys_init
c-----------------------------------------------------------------------
c     subprogram 3. eig_phys_init.
c     initializes the physical fields from input eigenfunction.
c-----------------------------------------------------------------------
      SUBROUTINE eig_phys_init(nmodes,keff)
      USE local
      USE input
      USE physdat
      USE fields
      USE polar_init
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

      INTEGER (i4) :: ibl,ixbl,iybl,ix0,ix1,iy0,iy1,mvert
      REAL(r8) :: gamma_r,gamma_i
      COMPLEX(r8) :: gamma_c
c-----------------------------------------------------------------------
c     initialize the perturbations.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        rb(ibl)%be=0
        rb(ibl)%ve=0
        rb(ibl)%pres=0
        rb(ibl)%prese=0
        rb(ibl)%tele=0
        rb(ibl)%tion=0
        rb(ibl)%nd=0
        rb(ibl)%conc=0
      ENDDO
      DO ibl=nrbl+1,nbl
        tb(ibl)%be=0
        tb(ibl)%ve=0
        tb(ibl)%pres=0
        tb(ibl)%prese=0
        tb(ibl)%tele=0
        tb(ibl)%tion=0
        tb(ibl)%nd=0
        tb(ibl)%conc=0
      ENDDO
c-----------------------------------------------------------------------
c     compute complex growth rate.
c-----------------------------------------------------------------------
      gamma_c=SQRT(CMPLX(-growth/(mtot*ndens),0._r8))
      gamma_r=REAL(gamma_c)
      gamma_i=REAL(gamma_c*CMPLX(0._r8,-1._r8))
c-----------------------------------------------------------------------
c     initialize perturbed quantities for rblocks.
c     The gato and nimrod toroidal mode number have opposite signs.
c     factoring in a pi phase shift for comparison with other nimrod
c     initializations, leads to a real-imaginary swap.
c-----------------------------------------------------------------------
      ibl=1
      ix1=0
      DO ixbl=1,nxbl
         ix0=ix1
         ix1=ix1+rb(ibl)%mx
         iy1=0
         DO iybl=1,nybl
            iy0=iy1
            iy1=iy1+rb(ibl)%my
            rb(ibl)%pres%fs(1,:,:,1)=eigvec(2,ix0:ix1,iy0:iy1)
     $                        +(0,1)*eigvec(1,ix0:ix1,iy0:iy1)
            rb(ibl)%prese%fs(1,:,:,1)=pe_frac*rb(ibl)%pres%fs(1,:,:,1)
            rb(ibl)%tele%fs(:,:,:,1)=rb(ibl)%prese%fs(:,:,:,1)/
     $                               (ndens*kboltz)
            rb(ibl)%tion%fs(:,:,:,1)=(1-pe_frac)*zeff*
     $        rb(ibl)%pres%fs(:,:,:,1)/(ndens*kboltz)
            rb(ibl)%ve%fs(1:3,:,:,1)=
     $                eigvec(6:8,ix0:ix1,iy0:iy1)*gamma_r
     $               +eigvec(3:5,ix0:ix1,iy0:iy1)*gamma_i
     $        +(0,1)*(eigvec(3:5,ix0:ix1,iy0:iy1)*gamma_r
     $               -eigvec(6:8,ix0:ix1,iy0:iy1)*gamma_i)
            rb(ibl)%be%fs(1:3,:,:,1)=eigvec(12:14,ix0:ix1,iy0:iy1)
     $                        +(0,1)*eigvec( 9:11,ix0:ix1,iy0:iy1)
            ibl=ibl+1
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE eig_phys_init
c-----------------------------------------------------------------------
c     subprogram 4. field_err_init.
c     adds vacuum perturbations to the magnetic-field data.
c-----------------------------------------------------------------------
      SUBROUTINE field_err_init(nmodes,keff)
      USE local
      USE input
      USE physdat
      USE fields
      USE lagr_quad_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff

      INTEGER(i4) :: im,nn,ibl,ix,iy,ib,iferr
      REAL(r8) :: xx,dx,dy,det
      COMPLEX(r8) :: bcoef,b1,b2
      COMPLEX(r8), DIMENSION(3) :: bvec
c-----------------------------------------------------------------------
c     loop over field-error input.
c-----------------------------------------------------------------------
      DO iferr=1,nfield_err
        IF (ferr_amp(iferr)==0._r8) CYCLE
c-----------------------------------------------------------------------
c       loop over Fourier components and find the coefficients for the
c       vacuum distribution.
c-----------------------------------------------------------------------
        DO im=1,nmodes
          IF (geom=='tor') THEN
            nn=NINT(keff(im))
            IF (nn==0_i4) THEN
              bcoef=ferr_amp(iferr)
            ELSE
              bcoef=  0.5_r8*ferr_amp(iferr)*COS(pi*ferr_phase(iferr))+
     $          (0,1)*0.5_r8*ferr_amp(iferr)*SIN(pi*ferr_phase(iferr))
            ENDIF
          ELSE
            nn=NINT(per_length*keff(im)/twopi)
            det=EXP(keff(im)*(xmax-xmin))-EXP(keff(im)*(xmin-xmax))
            b1=0.5_r8*ferr_amp(iferr)/det*
     $       (EXP(-keff(im)*xmin)-EXP(-keff(im)*xmax))*
     $       (COS(pi*ferr_phase(iferr))+(0,1)*SIN(pi*ferr_phase(iferr)))
            b2=0.5_r8*ferr_amp(iferr)/det*
     $       (EXP(keff(im)*xmin)-EXP(keff(im)*xmax))*
     $       (COS(pi*ferr_phase(iferr))+(0,1)*SIN(pi*ferr_phase(iferr)))
          ENDIF
          IF (nn/=ferr_n(iferr)) CYCLE
c-----------------------------------------------------------------------
c         loop over node points, determine the local perturbation field,
c         and add to the data storage.
c-----------------------------------------------------------------------
          DO ibl=1,nrbl
            DO ib=1,SIZE(rb(ibl)%be%ix0)
              DO iy=rb(ibl)%be%iy0(ib),rb(ibl)%be%my
                DO ix=rb(ibl)%be%ix0(ib),rb(ibl)%be%mx
                  dx=ix-rb(ibl)%be%ix0(ib)+rb(ibl)%be%dx(ib)
                  dy=iy-rb(ibl)%be%iy0(ib)+rb(ibl)%be%dy(ib)
                  CALL lagr_quad_eval(rb(ibl)%rz,dx,dy,0_i4)
                  IF (geom=='tor') THEN
                    IF (nn==0_i4) THEN
                      bvec=bcoef*(/1._r8,0._r8,0._r8/)*
     $                           (xo/rb(ibl)%rz%f(1))
                    ELSE
                      bvec=bcoef*(/1._r8,0._r8,0._r8/)*
     $                           (rb(ibl)%rz%f(1)/xo)**(nn-1_i4)
                      bvec(3)=(0,1)*bvec(1)
                    ENDIF
                  ELSE
                    IF (nn==0_i4) THEN
                      bvec=(/ferr_amp(iferr),0._r8,0._r8/)
                    ELSE
                      bvec(1)=b1*EXP( keff(im)*rb(ibl)%rz%f(1))
     $                       -b2*EXP(-keff(im)*rb(ibl)%rz%f(1))
                      bvec(2)=0._r8
                      bvec(3)=(0,1)*(b1*EXP( keff(im)*rb(ibl)%rz%f(1))
     $                              +b2*EXP(-keff(im)*rb(ibl)%rz%f(1)))
                    ENDIF
                  ENDIF
                  CALL lagr_quad_basis_add_loc(rb(ibl)%be,bvec,ib,
     $                                         ix,iy,im)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c-PRE     triangles
c-----------------------------------------------------------------------
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE field_err_init
