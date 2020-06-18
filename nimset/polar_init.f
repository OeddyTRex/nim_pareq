c-----------------------------------------------------------------------
c     file polar_init.f
c     performs initialization for nimrod polar-like grids.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. polar_init
c     1. polar_circlegrid_init.
c     2. polar_circle_b0.
c     3. polar_circle_lam.
c     4. polar_circle_pp.
c     5. polar_fluxgrid_init.
c     6. polar_rblock_init.
c     7. polar_tblock_init.
c     8. polar_seam_init.
c     9. polar_circlegrid_pack.
c     10. polar_assign_conc
c-----------------------------------------------------------------------
c     subprogram 0. polar_init.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE polar_init
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE input
      USE physdat
      USE lagr_quad_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(4), PRIVATE :: sqo
      REAL(r8), PRIVATE :: bo,jo,po,vo
      REAL(r8), DIMENSION(:,:), POINTER, PRIVATE :: sq
      TYPE(lagr_quad_2d_type), PRIVATE :: conceq,lxy,lbq,ljq,lpq,lvq,
     $                                    lnd,ldiff

      REAL(r8) :: rquot,omega0,growth
      REAL(r8), DIMENSION(14) :: eigveco
      REAL(r8), DIMENSION(:,:,:), POINTER :: eigvec

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. polar_circlegrid_init.
c     generates grid for nimrod code.
c-----------------------------------------------------------------------
      SUBROUTINE polar_circlegrid_init

      INTEGER(i4) :: ix,iy,iir,pmx,pmy,jx,jy,ib,ix0,iy0,ip,ixm,ng
      REAL(r8), DIMENSION(:), ALLOCATABLE :: r,bth_eq,bz_eq,p_eq,
     $          jth_eq,jz_eq,theta,q_eq
      REAL(r8), DIMENSION(0:mx) :: rvert,rin
      REAL(r8), DIMENSION(2,0:poly_degree*mx,0:poly_degree*my) :: tmprz
      REAL(r8) :: rayl,sina,cosa,ag,bg,rlast,dr1,xm,flag,tmpr,tmprm
      REAL(r8) :: jt0,jz0,bt,bz,jt,jz,thangle,dx,dy,gc1,gc22,jda,jdb,
     $            jdc,jdk,rmaj,bti
      REAL(r8), DIMENSION(0:poly_degree) :: x_node
      REAL(r8), DIMENSION(npack) :: asave,qpsave,wpsave
      REAL(r8), DIMENSION(ngr+poly_degree) :: wg,rg2
      REAL(r8), DIMENSION(0:ngr+poly_degree+1) :: rg1,bzg

      REAL(r8), DIMENSION(2) :: xxyy
      REAL(r8), DIMENSION(3) :: bvec
      INTEGER(i4) :: icoil
c-----------------------------------------------------------------------
c     interface block for lagr_1D external routine.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE lagr_1D(pd,x,al,dal,dmode)
        USE local
        INTEGER(i4), INTENT(IN) :: pd,dmode
        REAL(r8), INTENT(IN) :: x
        REAL(r8), DIMENSION(0:), INTENT(OUT) :: al,dal
        END SUBROUTINE lagr_1D
      END INTERFACE
c-----------------------------------------------------------------------
c     check domain limits.
c-----------------------------------------------------------------------
      IF (xmax<=xmin) CALL nim_stop
     $  ("Polar_circlegrid_init: xmax must be > xmin.")
c-----------------------------------------------------------------------
c     for grids with triangles at the magnetic axis, reduce mx so the
c     input is still the number of radial cells.
c-----------------------------------------------------------------------
      IF (pieflag=='tblock0'.AND.xmin==0) mx=mx-1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lxy,mx,my,2_i4,poly_degree)
      CALL lagr_quad_alloc(lbq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(ljq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(lpq,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(lnd,mx,my,1_i4,poly_degree)
      lnd=ndens
      x_node=(/lxy%dx(1:poly_degree),1._r8/)
c-----------------------------------------------------------------------
c     check if packing is needed.
c-----------------------------------------------------------------------
      If (MAXVAL(amp)>0._r8) THEN
        flag=1
      ELSE
        flag=0
      ENDIF
c-----------------------------------------------------------------------
c     start geometric if-then block.
c-----------------------------------------------------------------------
      geom_if: IF (geom=='lin') THEN
c-----------------------------------------------------------------------
c       compute grid and equilibrium fields, twice if needed for grid
c       packing.
c-----------------------------------------------------------------------
        eq_do: DO iir=1,2
          IF (iir==1) THEN
            pmx=mx
            pmy=my
          ELSE
            pmx=poly_degree*mx
            pmy=poly_degree*my
          ENDIF
          ALLOCATE(r(0:pmx),bth_eq(0:pmx),bz_eq(0:pmx),p_eq(0:pmx))
          ALLOCATE(jth_eq(0:pmx),jz_eq(0:pmx))
          ALLOCATE(theta(0:pmy))
          jth_eq=0
          jz_eq=0
c-----------------------------------------------------------------------
c         grid computations.
c-----------------------------------------------------------------------
          IF (iir==1) THEN
            theta=(/(iy,iy=0,pmy)/)*twopi/pmy
            IF (pieflag/='rblock') THEN
              IF (xmin/=0) CALL nim_stop
     $          ('xmin must be 0 if pieflag/=rblock.')
              xm=mxpie*xmax/(mxpie+mx)
            ELSE
              xm=xmin
            ENDIF
            rayl=xmax-xm
            IF (firsty/=0) THEN
              dr1=firsty*rayl/mx
            ELSE
              dr1=rayl/mx
            ENDIF
            r(0)=xm
            IF (dr1>2*rayl/mx) CALL nim_stop
     $        ('The linear dx distribution cannot handle this firsty.')
            ag=2*(rayl/mx-dr1)/(mx-1)
            bg=dr1-ag
            DO ix=1,mx-1
              r(ix)=r(ix-1)+(ag*ix+bg)
            ENDDO
            r(mx)=xmax
            rvert=r
            IF (flag==0) THEN
              DEALLOCATE(r,theta,bth_eq,bz_eq,jth_eq,jz_eq,p_eq)
              CYCLE eq_do
            ENDIF
          ELSE
            r(0)=rvert(0)
            DO iy=1,mx
              DO ix=1,poly_degree
                r(poly_degree*(iy-1)+ix)=
     $                          (rvert(iy-1)*(1-x_node(ix))+
     $                           rvert(iy  )*x_node(ix))
              ENDDO
            ENDDO
            theta(0)=0
            DO iy=1,my
              DO ix=1,poly_degree
                theta(poly_degree*(iy-1)+ix)=twopi*(iy-1+x_node(ix))/my
              ENDDO
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         start analytic-model-dependent computations.
c-----------------------------------------------------------------------
          IF (lamprof=='pitprs' .OR . lamprof=='qspec') THEN
            IF (pieflag/='rblock'.AND.xmin==0) THEN
              CALL polar_circle_pp((/0._r8,(r(ix),ix=0,pmx)/),bth_eq,
     $                             bz_eq,p_eq,jth_eq,jz_eq,jt0,jz0,
     $                             xmin,xmax,geom)
              jo=jz0
              po=beta*be0**2/(2*mu0)
            ELSE
              CALL polar_circle_pp(r,bth_eq(1:pmx),bz_eq(1:pmx),
     $                   p_eq(1:pmx),jth_eq(1:pmx),jz_eq(1:pmx),
     $                   jth_eq(0),jz_eq(0),xmin,xmax,geom)
              bth_eq(0)=r(0)*be0/(pit_0*xmax)
              bz_eq(0)=be0
              p_eq(0)=beta*be0**2/(2*mu0)
              jo=jz_eq(0)
              po=p_eq(0)
            ENDIF
c-----------------------------------------------------------------------
c         the 'grbrap' profile is specified by simple rational
c         functions.  the two constants are specified with pit_0 (the
c         normalized pitch on axis) and pit_2.
c-----------------------------------------------------------------------
          ELSE IF (lamprof=='grbrap') THEN
            gc1=1._r8/((xmax-xmin)*pit_0)
            gc22=pit_2**2
            bth_eq=be0*gc1*r/(1._r8+gc22*r**2)
            bz_eq=be0
            jth_eq=0._r8
            jz_eq=2._r8*be0*gc1/(1._r8+gc22*r**2)**2/mu0
            p_eq=0.5_r8*be0**2*gc1**2/(mu0*gc22)*
     $           (1._r8/(1._r8+gc22*r**2)**2-1._r8/(1._r8+gc22)**2)
            po=p_eq(0)
            bo=be0
c-----------------------------------------------------------------------
c         the 'jardel' profile needs integrals over radius.  use
c         Gaussian quadrature between nimrod-node locations to integrate
c         pressure.  between those quadrature points, use Gaussian
c         integration to obtain the b-theta integral using r/R for the
c         variable of integration.  note that this profile is only
c         for 0<=r<=1.
c-----------------------------------------------------------------------
          ELSE IF (lamprof=='jardel') THEN
            ng=ngr+poly_degree
            rmaj=per_length/twopi
            jda=(2._r8-0.5_r8*alpha)*pit_0**2*rmaj**2/(1._r8-pit_0**2)
            jdb=(pit_0**2/rmaj**2)/(1._r8-2._r8*pit_0**2)
            jdc=(pit_0**2*rmaj**2)/(1._r8-2._r8*pit_0**2)
            jdk=(1._r8-pit_0**2)/(1._r8-2._r8*pit_0**2)
            bth_eq(0)=0._r8
            bz_eq(0)=be0
            jz_eq(0)=2._r8*be0/(mu0*pit_0)
            p_eq(0)=0._r8
            bti=0._r8
            rg1(0)=0._r8
            DO ix=1,pmx    !   first gauleg is to find grid for 2nd int
              CALL gauleg(r(ix-1)/rmaj,r(ix)/rmaj,rg1(1:ng),wg,ng)
              rg1(ng+1)=r(ix)/rmaj
              DO iy=1,ng+1
                CALL gauleg(rg1(iy-1),rg1(iy),rg2,wg,ng)
                bti=bti+jdk*SUM(wg*(2._r8*rg2+jda*rg2**3)/
     $                            (jdb+rg2**2+jdc*rg2**4))
                bzg(iy)=be0*(1._r8-(rmaj*rg1(iy))**2)*exp(-bti)
              ENDDO
              bth_eq(ix)=be0*r(ix)*exp(-bti)/pit_0
              bz_eq(ix)=bzg(ng+1)
              tmpr=r(ix)/rmaj
              tmpr=jdk*((2._r8*tmpr+jda*tmpr**3)/
     $                 (jdb+tmpr**2+jdc*tmpr**4))/rmaj
              jth_eq(ix)=pit_0*bth_eq(ix)/(mu0*r(ix))*
     $                   ((1._r8-r(ix)**2)*tmpr+2._r8*r(ix))
              jz_eq(ix)=(bth_eq(ix)/mu0)*(2._r8/r(ix)-tmpr)
              CALL gauleg(r(ix-1)/rmaj,r(ix)/rmaj,rg1(1:ng),wg,ng)
              p_eq(ix)=p_eq(ix-1)-0.5_r8*alpha*rmaj**4*
     $          SUM(wg*rg1(1:ng)*(rg1(1:ng)*bzg(1:ng)/
     $                           (1._r8-(rmaj*rg1(1:ng))**2))**2)
              rg1(0)=rg1(ng+1)
            ENDDO
            p_eq=(p_eq-p_eq(pmx))/mu0
            bo=be0
            po=p_eq(0)
c-----------------------------------------------------------------------
c         cases with specified parallel current profiles
c-----------------------------------------------------------------------
          ELSE
            IF (pieflag/='rblock'.AND.xmin==0) THEN
              CALL polar_circle_b0((/0._r8,(r(ix),ix=0,pmx)/),bth_eq,
     $                             bz_eq,p_eq,jth_eq,jz_eq,xmin,xmax)
            ELSE
              CALL polar_circle_b0(r,bth_eq(1:pmx),bz_eq(1:pmx),
     $                             p_eq(1:pmx),jth_eq(1:pmx),
     $                             jz_eq(1:pmx),xmin,xmax)
              IF (r(0)>0) THEN
                bth_eq(0)=be0*SIN(pi*thetab)
                bz_eq(0)=be0*COS(pi*thetab)
                jth_eq(0)=polar_circle_lam(r(0),bth_eq(0),bz_eq(0),
     $                                     xmin,xmax)*bth_eq(0)/mu0
                jz_eq(0) =polar_circle_lam(r(0),bth_eq(0),bz_eq(0),
     $                                     xmin,xmax)*bz_eq(0)/mu0
              ELSE
                bth_eq(0)=0
                bz_eq(0)=be0
                jth_eq(0)=0
                jz_eq(0)=polar_circle_lam(r(0),bth_eq(0),bz_eq(0),
     $                                    xmin,xmax)*bz_eq(0)/mu0
              ENDIF
              p_eq(0)=beta*be0**2/(2*mu0)
            ENDIF 
            po=beta*be0**2/(2*mu0)
          ENDIF
          bo=be0
c-----------------------------------------------------------------------
c         compute and write the Suydam parameter and R*q.  this version
c         is general.
c-----------------------------------------------------------------------
          IF (iir==2) THEN
            CALL open_bin(xy_unit,"suydam.bin","UNKNOWN","REWIND",
     $                    32_i4)
            DO ix=1,pmx
              tmpr=r(ix)*(bz_eq(ix)*2/r(ix)-mu0*jth_eq(ix)
     $                    -mu0*bz_eq(ix)*jz_eq(ix)/bth_eq(ix) )**2
              IF (ABS(tmpr)<1.e-14) tmpr=1.e-14
              WRITE(xy_unit) (/REAL(ix,4),REAL(r(ix),4),REAL(
     $          2*mu0*(jz_eq(ix)*bth_eq(ix)-jth_eq(ix)*bz_eq(ix))/
     $                 tmpr,4),
     $          REAL(r(ix)*bz_eq(ix)/bth_eq(ix),4)/)
            ENDDO
            CALL close_bin(xy_unit,"suydam.bin")
          ENDIF
c-----------------------------------------------------------------------
c         Pack around the rational surfaces if required.
c         After packing, recalculate 2d fields by looping again.
c-----------------------------------------------------------------------
          IF (flag/=0._r8.AND.iir<2) THEN
c-----------------------------------------------------------------------
c           evaluate q at the cell-centers of a uniform grid.
c-----------------------------------------------------------------------
            ALLOCATE(q_eq(mx))
            DO ix=1,mx
              q_eq(ix)=0.5*(r(ix)+r(ix-1))*(bz_eq(ix)+bz_eq(ix-1))/
     $                    ((bth_eq(ix)+bth_eq(ix-1))*per_length/(2.*pi))
            ENDDO
            CALL polar_circlegrid_pack(r,q_eq)
            rvert=r
            DEALLOCATE(r,theta,bth_eq,bz_eq,jth_eq,jz_eq,p_eq,q_eq)
          ENDIF
        ENDDO eq_do
c-----------------------------------------------------------------------
c       generate 2d lagrange-type grid.
c-----------------------------------------------------------------------
        DO ix=0,mx
          DO iy=0,my
            lxy%fs(1,ix,iy)=
     $        r(poly_degree*ix)*COS(theta(poly_degree*iy))+xo
            lxy%fs(2,ix,iy)=
     $        r(poly_degree*ix)*SIN(theta(poly_degree*iy))+yo
            DO ib=1,poly_degree-1
              IF (ix>0) THEN
                lxy%fsh(1,ib,ix,iy)=
     $            r(poly_degree*(ix-1)+ib)*COS(theta(poly_degree*iy))+xo
                lxy%fsh(2,ib,ix,iy)=
     $            r(poly_degree*(ix-1)+ib)*SIN(theta(poly_degree*iy))+yo
              ENDIF
              IF (iy>0) THEN
                lxy%fsv(1,ib,ix,iy)=
     $            r(poly_degree*ix)*COS(theta(poly_degree*(iy-1)+ib))+xo
                lxy%fsv(2,ib,ix,iy)=
     $            r(poly_degree*ix)*SIN(theta(poly_degree*(iy-1)+ib))+yo
              ENDIF
            ENDDO
            IF (ix>0.AND.iy>0) THEN
              DO ib=1,(poly_degree-1)**2
                jx=MODULO(ib-1_i4,poly_degree-1_i4)+1
                jy=(ib-1)/(poly_degree-1)+1
                lxy%fsi(1,ib,ix,iy)=r(poly_degree*(ix-1)+jx)*
     $            COS(theta(poly_degree*(iy-1)+jy))+xo
                lxy%fsi(2,ib,ix,iy)=r(poly_degree*(ix-1)+jx)*
     $            SIN(theta(poly_degree*(iy-1)+jy))+yo
              ENDDO
            ENDIF
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       set 2d lagrange-type arrays for the equilibrium fields.
c-----------------------------------------------------------------------
        DO iy=0,my
          jy=poly_degree*iy
          lbq%fs(1,:,iy)=-bth_eq(0:pmx:poly_degree)*SIN(theta(jy))
          lbq%fs(2,:,iy)= bth_eq(0:pmx:poly_degree)*COS(theta(jy))
          lbq%fs(3,:,iy)= bz_eq(0:pmx:poly_degree)
          ljq%fs(1,:,iy)=-jth_eq(0:pmx:poly_degree)*SIN(theta(jy))
          ljq%fs(2,:,iy)= jth_eq(0:pmx:poly_degree)*COS(theta(jy))
          ljq%fs(3,:,iy)= jz_eq(0:pmx:poly_degree)
          lpq%fs(1,:,iy)= p_eq(0:pmx:poly_degree)
        ENDDO
        DO iy=0,my
          DO ix=0,mx
            DO ib=1,poly_degree-1
              IF (ix>0) THEN
                lbq%fsh(1,ib,ix,iy)=-bth_eq(poly_degree*(ix-1)+ib)*
     $            SIN(theta(poly_degree*iy))
                lbq%fsh(2,ib,ix,iy)= bth_eq(poly_degree*(ix-1)+ib)*
     $            COS(theta(poly_degree*iy))
                lbq%fsh(3,ib,ix,iy)= bz_eq (poly_degree*(ix-1)+ib)
                ljq%fsh(1,ib,ix,iy)=-jth_eq(poly_degree*(ix-1)+ib)*
     $            SIN(theta(poly_degree*iy))
                ljq%fsh(2,ib,ix,iy)= jth_eq(poly_degree*(ix-1)+ib)*
     $            COS(theta(poly_degree*iy))
                ljq%fsh(3,ib,ix,iy)= jz_eq (poly_degree*(ix-1)+ib)
                lpq%fsh(1,ib,ix,iy)= p_eq (poly_degree*(ix-1)+ib)
              ENDIF
              IF (iy>0) THEN
                lbq%fsv(1,ib,ix,iy)=-bth_eq(poly_degree*ix)*
     $            SIN(theta(poly_degree*(iy-1)+ib))
                lbq%fsv(2,ib,ix,iy)= bth_eq(poly_degree*ix)*
     $            COS(theta(poly_degree*(iy-1)+ib))
                lbq%fsv(3,ib,ix,iy)= bz_eq (poly_degree*ix)
                ljq%fsv(1,ib,ix,iy)=-jth_eq(poly_degree*ix)*
     $            SIN(theta(poly_degree*(iy-1)+ib))
                ljq%fsv(2,ib,ix,iy)= jth_eq(poly_degree*ix)*
     $            COS(theta(poly_degree*(iy-1)+ib))
                ljq%fsv(3,ib,ix,iy)= jz_eq (poly_degree*ix)
                lpq%fsv(1,ib,ix,iy)= p_eq (poly_degree*ix)
              ENDIF
            ENDDO
            IF (ix>0.AND.iy>0) THEN
              DO ib=1,(poly_degree-1)**2
                jx=MODULO(ib-1_i4,poly_degree-1_i4)+1
                jy=(ib-1)/(poly_degree-1)+1
                lbq%fsi(1,ib,ix,iy)=-bth_eq(poly_degree*(ix-1)+jx)*
     $            SIN(theta(poly_degree*(iy-1)+jy))
                lbq%fsi(2,ib,ix,iy)= bth_eq(poly_degree*(ix-1)+jx)*
     $            COS(theta(poly_degree*(iy-1)+jy))
                lbq%fsi(3,ib,ix,iy)= bz_eq (poly_degree*(ix-1)+jx)
                ljq%fsi(1,ib,ix,iy)=-jth_eq(poly_degree*(ix-1)+jx)*
     $            SIN(theta(poly_degree*(iy-1)+jy))
                ljq%fsi(2,ib,ix,iy)= jth_eq(poly_degree*(ix-1)+jx)*
     $            COS(theta(poly_degree*(iy-1)+jy))
                ljq%fsi(3,ib,ix,iy)= jz_eq (poly_degree*(ix-1)+jx)
                lpq%fsi(1,ib,ix,iy)= p_eq (poly_degree*(ix-1)+jx)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     toroidal cases--firstx allows the grid to shift from the center
c     of the circle.
c-----------------------------------------------------------------------
      ELSE geom_if
        pmx=poly_degree*mx
        pmy=poly_degree*my
        ALLOCATE(theta(0:pmy))
        theta(0)=0
        DO iy=1,my
          DO ix=1,poly_degree
            theta(poly_degree*(iy-1)+ix)=twopi*(iy-1+x_node(ix))/my
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       allow radial packing by normalized minor-radius coordinate.
c       if any wpack values are < 0, separate packing is used on the
c       inboard side.  packing on the inboard side is specified
c	by packing indices i where wpack(i)<0, and each width is
c	|wpack(i)|.
c-----------------------------------------------------------------------
        IF (flag>0) THEN
          ALLOCATE(r(0:mx))
          asave=amp
          qpsave=qpack
          wpsave=wpack
          DO ip=1,npack
            IF (wpack(ip)<0) THEN
              amp(ip)=0._r8
              qpack(ip)=0._r8
              wpack(ip)=0._r8
            ENDIF
          ENDDO
          DO ix=0,mx
            rvert(ix)=REAL(ix,r8)/REAL(mx,r8)
            r(ix)=rvert(ix)
            rin(ix)=rvert(ix)
          ENDDO
          CALL polar_circlegrid_pack(r,rvert)
          IF (MINVAL(wpsave)==0._r8) THEN
            rin=r
          ELSE
            amp=0._r8
            qpack=0._r8
            wpack=0._r8
            ix=1
            DO ip=1,npack
              IF (wpsave(ip)<0) THEN
                amp(ix)=asave(ip)
                qpack(ix)=qpsave(ip)
                wpack(ix)=ABS(wpsave(ip))
                ix=ix+1
              ENDIF
            ENDDO
            CALL polar_circlegrid_pack(rin,rvert)
            amp=asave
            qpack=qpsave
            wpack=wpsave
          ENDIF
        ENDIF
        DO iy=0,pmy
          rayl=SQRT(xmax**2+firstx**2-2*xmax*firstx*COS(theta(iy)))
          sina= xmax*SIN(theta(iy))/rayl
          cosa=(xmax*COS(theta(iy))-firstx)/rayl
          IF (iy==0) THEN
            IF (pieflag=='rblock') THEN
              xm=xmin
            ELSE
              xm=mxpie*rayl/(mx+mxpie)
            ENDIF
            IF (firsty/=0) THEN
              dr1=firsty*(rayl-xm)/mx
            ELSE
              dr1=(rayl-xm)/mx
            ENDIF
          ENDIF
          rlast=xm
          IF (dr1>2*(rayl-xm)/mx) CALL nim_stop
     $      ('The linear dx distribution cannot handle this firstx.')
          ag=2*((rayl-xm)/mx-dr1)/(mx-1)
          bg=dr1-ag
          tmprm=xm
          tmprz(1,0,iy)=cosa*tmprm+xo+firstx
          tmprz(2,0,iy)=sina*tmprm+yo
          DO ix=1,mx
            IF (flag>0) THEN
              tmpr=xm+(rayl-xm)*(  r(ix)*COS(0.5_r8*theta(iy))**2+
     $                           rin(ix)*SIN(0.5_r8*theta(iy))**2)
            ELSE
              tmpr=tmprm+(ag*ix+bg)
            ENDIF
            DO ib=1,poly_degree
              tmprz(1,poly_degree*(ix-1)+ib,iy)=xo+firstx+
     $          cosa*(tmprm*(1-x_node(ib))+tmpr*x_node(ib))
              tmprz(2,poly_degree*(ix-1)+ib,iy)=yo+
     $          sina*(tmprm*(1-x_node(ib))+tmpr*x_node(ib))
            ENDDO
            tmprm=tmpr
          ENDDO
        ENDDO
        IF (MINVAL(tmprz(1,:,:))<=0) CALL nim_stop
     $    ('Toroidal geometry circular cases need to have R>0.')
c-----------------------------------------------------------------------
c       load mesh into a lagrange-quad type structure.
c-----------------------------------------------------------------------
        lxy%fs=tmprz(:,0:pmx:poly_degree,0:pmy:poly_degree)
        IF (poly_degree>1) THEN
          DO ix=1,mx
            DO ib=1,poly_degree-1
              lxy%fsh(:,ib,ix,:)=
     $          tmprz(:,poly_degree*(ix-1)+ib,0:pmy:poly_degree)
            ENDDO
          ENDDO
          DO iy=1,my
            DO ib=1,poly_degree-1
              lxy%fsv(:,ib,:,iy)=
     $          tmprz(:,0:pmx:poly_degree,poly_degree*(iy-1)+ib)
            ENDDO
          ENDDO
          DO iy=1,my
            DO ix=1,mx
              DO ib=1,(poly_degree-1)**2
                jx=MODULO(ib-1_i4,poly_degree-1_i4)+1
                jy=(ib-1)/(poly_degree-1)+1
                lxy%fsi(:,ib,ix,iy)=
     $            tmprz(:,poly_degree*(ix-1)+jx,poly_degree*(iy-1)+jy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vacuum fields for toroidal cases without a Grad-Shafranov soln.
c-----------------------------------------------------------------------
        bo=xo*be0
        po=beta*be0**2/(2*mu0)
        DO ib=1,poly_degree**2
          ix0=lxy%ix0(ib)
          iy0=lxy%iy0(ib)
          DO iy=iy0,my
            DO ix=ix0,mx
              CALL lagr_quad_basis_assign_loc(lbq,
     $          (/0._r8,be0*SIN(pi*thetab),
     $            xo*be0*COS(pi*thetab)/),ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        lpq=po
        ljq=0
        jo=0
        IF (pieflag/='rblock') THEN
          xo=xo+firstx
          yo=0
        ENDIF
      ENDIF geom_if
c-----------------------------------------------------------------------
c     add pressure offset.
c-----------------------------------------------------------------------
      IF (pres_offset/=0) THEN
        lpq%fs=lpq%fs+pres_offset*po
        IF (poly_degree>1) THEN
          lpq%fsh=lpq%fsh+pres_offset*po
          lpq%fsv=lpq%fsv+pres_offset*po
          lpq%fsi=lpq%fsi+pres_offset*po
        ENDIF
        po=po*(1._r8+pres_offset)
      ENDIF
c-----------------------------------------------------------------------
c     reset the number density profile to satisfy 
c     p/p0=(n/ndens)**gamma_nimset if gamma_nimset>0.
c-----------------------------------------------------------------------
      IF (gamma_nimset>0) THEN
        DO ib=1,SIZE(lnd%ix0)
          DO iy=lnd%iy0(ib),lnd%my
            DO ix=lnd%ix0(ib),lnd%mx
              dx=ix-lnd%ix0(ib)+lnd%dx(ib)
              dy=iy-lnd%iy0(ib)+lnd%dy(ib)
              CALL lagr_quad_eval(lpq,dx,dy,0_i4)
              tmpr=ndens*(lpq%f(1)/po)**(1._r8/gamma_nimset)
              CALL lagr_quad_basis_assign_loc(lnd,(/tmpr/),ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find vacuum magnetic field components due to solenoids.
c-----------------------------------------------------------------------
      IF (ncoil>0.AND.geom=='tor') THEN
        bvec(3)=0._r8
        DO ib=1,SIZE(lxy%ix0)
          DO iy=lxy%iy0(ib),lxy%my
            DO ix=lxy%ix0(ib),lxy%mx
              dx=ix-lxy%ix0(ib)+lxy%dx(ib)
              dy=iy-lxy%iy0(ib)+lxy%dy(ib)
              CALL lagr_quad_eval(lxy,dx,dy,1_i4)
              CALL brz_eval(bvec(1:2),lxy%f(1),lxy%f(2),ncoil,
     $                      coil_r,coil_z,coil_current,mu0)
              CALL lagr_quad_basis_add_loc(lbq,bvec,ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_circlegrid_init
c-----------------------------------------------------------------------
c     subprogram 2. polar_circle_b0.
c     generates equilibrium fields in linear geometry with circular
c     cross section for various parallel current profiles.
c-----------------------------------------------------------------------
      SUBROUTINE polar_circle_b0(r,bth,bz,pr,jt,jz,r0,r1)

      REAL(r8), DIMENSION(1:), INTENT(OUT) :: bth,bz,pr,jt,jz
      REAL(r8), DIMENSION(0:), INTENT(IN) :: r
      REAL(r8), INTENT(IN) :: r0,r1

      INTEGER(i4) :: nr,i,isub,nsub
      INTEGER(i4), PARAMETER :: nsub_max=1024
      REAL(r8) :: lam,rbthp,bzp,bthm,bzm,rc,dr,rp,rm,bthl,bzl,bmag,xdim,
     $            x1,dpr,bmag2
      REAL(r8), PARAMETER :: btol=1.e-7
      REAL(r8), DIMENSION(4) :: rjz,jth 
c-----------------------------------------------------------------------
c     check for vacuum region:
c-----------------------------------------------------------------------
      IF (xvac>0) THEN
        xdim=xvac-r0
        x1=xvac
      ELSE
        xdim=r1-r0
        x1=r1
      ENDIF
c-----------------------------------------------------------------------
c     do a Runge-Kutta integration to find the
c     equilibrium field from the specified current profile.
c     subdivide the intervals as needed to acheive convergence.
c-----------------------------------------------------------------------
      nr=SIZE(bth)
      IF (r(0)>0) THEN
        bthl=be0*SIN(pi*thetab)
        bzl=be0*COS(pi*thetab)
      ELSE
        bthl=0
        bzl=be0
      ENDIF
      DO i=1,nr
        nsub=1
        sub_size: DO
          IF (i>1) THEN
            bthm=bth(i-1)
            bzm=bz(i-1)
          ELSE
            IF (r(0)>0) THEN
              bthm=be0*SIN(pi*thetab)
              bzm=be0*COS(pi*thetab)
            ELSE
              bthm=0
              bzm=be0
            ENDIF
          ENDIF
          dr=(r(i)-r(i-1))/nsub
          rm=r(i-1)
          sub_step: DO isub=1,nsub
            rp=rm+dr
            rc=rm+0.5*dr
c-----------------------------------------------------------------------
c           Step 1:
c-----------------------------------------------------------------------
            lam=polar_circle_lam(rm,bthm,bzm,r0,r1)
            dpr=dp_fun(rm)
            bmag2=bthm**2+bzm**2
            jth(1)=lam*bthm+bzm*dpr/bmag2
            rjz(1)=rm*(lam*bzm-bthm*dpr/bmag2)
            rbthp=rm*bthm+0.5*dr*rjz(1)
            bzp=bzm-0.5*dr*jth(1)
c-----------------------------------------------------------------------
c           Step 2:
c-----------------------------------------------------------------------
            lam=polar_circle_lam(rc,rbthp/rc,bzp,r0,r1)
            dpr=dp_fun(rc)
            bmag2=(rbthp/rc)**2+bzp**2
            jth(2)=lam*rbthp/rc+bzp*dpr/bmag2
            rjz(2)=rc*(lam*bzp-rbthp/rc*dpr/bmag2)
            rbthp=rm*bthm+0.5*dr*rjz(2)
            bzp=bzm-0.5*dr*jth(2)
c-----------------------------------------------------------------------
c           Step 3:
c-----------------------------------------------------------------------
            lam=polar_circle_lam(rc,rbthp/rc,bzp,r0,r1)
            bmag2=(rbthp/rc)**2+bzp**2
            jth(3)=lam*rbthp/rc+bzp*dpr/bmag2
            rjz(3)=rc*(lam*bzp-rbthp/rc*dpr/bmag2)
            rbthp=rm*bthm+dr*rjz(3)
            bzp=bzm-dr*jth(3)
c-----------------------------------------------------------------------
c           Step 4:
c-----------------------------------------------------------------------
            lam=polar_circle_lam(rp,rbthp/rp,bzp,r0,r1)
            dpr=dp_fun(rp)
            bmag2=(rbthp/rp)**2+bzp**2
            jth(4)=lam*rbthp/rp+bzp*dpr/bmag2
            rjz(4)=rp*(lam*bzp-rbthp/rp*dpr/bmag2)
            rbthp=rm*bthm+dr/6.*(SUM(rjz(1:4:3))+2.*SUM(rjz(2:3)))
            bzp=bzm-dr/6.*(SUM(jth(1:4:3))+2.*SUM(jth(2:3)))
            rm=rp
            bthm=rbthp/rp
            bzm=bzp
          ENDDO sub_step
c-----------------------------------------------------------------------
c         check convergence.
c-----------------------------------------------------------------------
          bmag=SQRT(bthl**2+bzl**2)
          IF ( ABS(bthm-bthl)/bmag<=btol .AND.
     $         ABS(bzm-bzl)/bmag<=btol ) THEN
            bth(i)=bthm
            bz(i)=bzm
            EXIT sub_size
          ENDIF
          bthl=bthm
          bzl=bzm
          nsub=2*nsub
          IF (nsub>nsub_max)
     $      CALL nim_stop('Polar_circle_b0 is not converging.')
        ENDDO sub_size
      ENDDO

c-----------------------------------------------------------------------
c     pressure profile.
c-----------------------------------------------------------------------
      WHERE (r(1:)>x1)
        pr=0
      ELSE WHERE
        pr=beta*be0**2*(1+pres_2*((r(1:)-r0)/xdim)**2
     $                   +pres_4*((r(1:)-r0)/xdim)**4)/(2*mu0)
      END WHERE
c-----------------------------------------------------------------------
c     current density.
c-----------------------------------------------------------------------
      DO i=1,nr
        bmag2=bth(i)**2+bz(i)**2
        dpr=dp_fun(r(i))
        lam=polar_circle_lam(r(i),bth(i),bz(i),r0,r1)
        jt(i)=(lam*bth(i)+bz(i)*dpr/bmag2)/mu0
        jz(i)=(lam*bz(i)-bth(i)*dpr/bmag2)/mu0
      ENDDO
      IF (r(0)>0) THEN
        bthl=be0*SIN(pi*thetab)
        bzl=be0*COS(pi*thetab)
        jo=polar_circle_lam(r(0),bthl,bzl,r0,r1)*be0/
     $     (mu0*COS(pi*thetab))
      ELSE
        jo=polar_circle_lam(0._r8,0._r8,be0,r0,r1)*be0/mu0
      ENDIF
c-----------------------------------------------------------------------
c     compute and write the Suydam parameter and R*q.
c-----------------------------------------------------------------------
      CALL open_bin(xy_unit,"suydam.bin","UNKNOWN","REWIND",32_i4)
      DO i=1,nr
        WRITE(xy_unit) (/REAL(i,4),REAL(r(i),4),REAL(
     $    -2*r(i)*dp_fun(r(i))/xdim/
     $    (2*bz(i)-r(i)*mu0*(jt(i)+jz(i)*bz(i)/bth(i)))**2,4),
     $    REAL(r(i)*bz(i)/bth(i),4)/)
      ENDDO
      CALL close_bin(xy_unit,"suydam.bin")
c-----------------------------------------------------------------------
c     end of executable statements.
c-----------------------------------------------------------------------
      RETURN

      CONTAINS
c-----------------------------------------------------------------------
c       compute mu0*dp/dr.
c
c       if there is rigid rotation, incorporate its radial force
c       density  =>  mu0*(dp/dr - ndens*mtot*r*v0**2/rmax**2)
c-----------------------------------------------------------------------
        FUNCTION dp_fun(rad)

        REAL(r8), INTENT(IN) :: rad
        REAL(r8) :: dp_fun

        IF (rad>xvac.AND.xvac>0) THEN
          dp_fun=0
          RETURN
        ENDIF
        dp_fun=beta*be0**2/(2*xdim)
     $             *(2*pres_2*((rad-r0)/xdim)
     $              +4*pres_4*((rad-r0)/xdim)**3)

        IF (eq_flow=='rigid') THEN
          dp_fun=dp_fun-mu0*ndens*mtot*rad*ve0**2/xmax**2
        ENDIF

        END FUNCTION dp_fun
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE polar_circle_b0
c-----------------------------------------------------------------------
c     subprogram 3. polar_circle_lam
c     finds the parallel current density for a given radius
c     based on a specified profile.  note that the returned value
c     is in m**(-1).
c-----------------------------------------------------------------------
      FUNCTION polar_circle_lam(r,bth,bz,r0,r1)

      REAL(r8), INTENT(IN) :: r,bth,bz,r0,r1
      REAL(r8) :: polar_circle_lam,xdim,x1
c-----------------------------------------------------------------------
c     check for vacuum region:
c-----------------------------------------------------------------------
      IF (xvac>0) THEN
        xdim=xvac-r0
        x1=xvac
      ELSE
        xdim=r1-r0
        x1=r1
      ENDIF
c-----------------------------------------------------------------------
c     modified bessel function model:
c-----------------------------------------------------------------------
      SELECT CASE(TRIM(lamprof))
      CASE('mbfm')
        IF ((r-r0)/xdim<=rbreak) THEN
          polar_circle_lam=lam0/xdim
        ELSE
          polar_circle_lam=lam0*( (x1-r)/(1-rbreak) )/xdim**2
        ENDIF
c-----------------------------------------------------------------------
c     alpha model:
c-----------------------------------------------------------------------
      CASE('alpha')
        polar_circle_lam=lam0*(1-((r-r0)/xdim)**alpha)/xdim
c-----------------------------------------------------------------------
c     paramagnetic pinch:
c-----------------------------------------------------------------------
      CASE('para')
        polar_circle_lam=lam0*be0*bz/(xdim*(bz**2+bth**2))
        IF (ds_use=='equil'.OR.xvac>0)
     $    polar_circle_lam=polar_circle_lam/
     $      MIN(dvac,(1+(SQRT(dvac)-1)*((r-r0)/xdim)**dexp)**2)
      CASE DEFAULT
        CALL nim_stop('Unrecognized lambda profile.')
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION polar_circle_lam
c-----------------------------------------------------------------------
c     subprogram 4. polar_circle_pp.
c     generates equilibrium fields in linear geometry with circular
c     cross section for pitch (q) & pressure profiles.  
c       For lamprof="pitprs":
c       Pitch is defined as (r*Bz/a*Bth), i.e. q*R/a, and it's
c       specified through the input variables pit_0, pit_2, and pit_4:
c
c       pitch(r)=pit_0*(1+pit_2*x**2+pit_4*x**4) ,
c
c       where x=(r-r0)/(r1-r0).
c       The pressure is specified through beta, pres_2, pres_4 as
c
c       p(r)=beta*be0**2/(2*mu0)*(1+pres_2*x**2+pres_4*x**4) .
c
c       For lamprof="qspec":
c               See Holmes et.al., Phys. Fluid 26 (1983) 2569
c       The q profile (really the pitch when coded) is specified through
c	the input variables pit_0, pit_2, and pit_4:
c
c            q(r)=pit_0*(1+(x/pit_2)**(2*pit_4))**(1/pit_4) .
c
c	See the paper for details.
c       Note that the qspec computation assumes that r0=0.
c
c       jt0 and jz0 return the current densities at r(0).
c       For geom='tor', jt0 is J_theta/r.
c-----------------------------------------------------------------------
      SUBROUTINE polar_circle_pp(r,bth,bz,pr,jt,jz,jt0,jz0,r0,r1,geom)

      REAL(r8), DIMENSION(1:), INTENT(OUT) :: bth,bz,pr,jt,jz
      REAL(r8), DIMENSION(0:), INTENT(IN) :: r
      REAL(r8), INTENT(OUT) :: jt0,jz0
      REAL(r8), INTENT(IN) :: r0,r1
      CHARACTER(*), INTENT(IN) :: geom

      INTEGER(i4) :: nr,i,isub,nsub
      INTEGER(i4), PARAMETER :: nsub_max=1024
      REAL(r8) :: xdim,ydim,bthp,bthm,rc,dr,prm,drp,rcp,jpar,jper,bmag,
     $            kpr,pit0,pit2,rm,rp,prl,prp,bthl
      REAL(r8), PARAMETER :: tol=1.e-10
      REAL(r8), DIMENSION(4) :: dbt, dpr
c-----------------------------------------------------------------------
c     preliminaries.
c-----------------------------------------------------------------------
      xdim=r1-r0
      kpr=1._r8
      pit0=pit_0
      pit2=pit_2
      IF (gridshape=='circ'.OR.gridshape=='rect_cir') THEN
        ydim=per_length
      ELSE
        ydim=ymax-ymin
      ENDIF
      nr=SIZE(bth)
c-----------------------------------------------------------------------
c     For qspec, user specifies q, but we really want to work with pitch
c     The pit_2 variable is given as r_o/a in FAR notation.  
c     Translate to match these definitions.
c-----------------------------------------------------------------------
      IF (lamprof=='qspec') THEN
        pit0=(ydim/(2.*pi)/xdim)*pit_0
        pit2=r1*pit_2
c-----------------------------------------------------------------------
c       For qspec (only), integrate the pressure profile to determine
c       the normalization constant, kpr, which will be used in dp_fun
c       when integrating (a*btheta/r) below.
c-----------------------------------------------------------------------
        prl=mu0
        DO i=1,nr
          nsub=1
          sub_sizep: DO
            IF (i>1) THEN
              prm=pr(i-1)
            ELSE
              prm=mu0
            ENDIF
            dr=(r(i)-r(i-1))/nsub
            rm=r(i-1)
            sub_stepp: DO isub=1,nsub
              rc=rm+0.5*dr
              rp=rm+dr
              dpr(1)=dp_fun(rm)			!  Step 1
              dpr(2:3)=dp_fun(rc)		!  Steps 2 & 3
              dpr(4)=dp_fun(rp)			!  Step 4
              rm=rp
              prm=prm+dr/(6.*r1)*(SUM(dpr(1:4:3))+2.*SUM(dpr(2:3)))
            ENDDO sub_stepp
c-----------------------------------------------------------------------
c           check convergence.
c-----------------------------------------------------------------------
            IF (ABS(prm-prl)/(0.5*(ABS(prm)+ABS(prl)))<=tol) THEN
              pr(i)=prm
              EXIT sub_sizep
            ENDIF
            prl=prm
            nsub=2*nsub
            IF (nsub>nsub_max)
     $        CALL nim_stop('Polar_circle_pp is not converging on p.')
          ENDDO sub_sizep
        ENDDO
        kpr=(be0**2*beta)/(2._r8*(mu0-pr(nr)))
        WRITE(nim_wr,*) "Polar_circle_pp: qspec kpr= ",kpr
        pr=(pr-pr(nr))*(kpr/mu0)
c-----------------------------------------------------------------------
c     just evaluate pressure for pitprs.
c-----------------------------------------------------------------------
      ELSE
        pr=beta*be0**2*(1+pres_2*((r(1:nr)-r0)/xdim)**2
     $                   +pres_4*((r(1:nr)-r0)/xdim)**4)/(2*mu0)
      ENDIF
c-----------------------------------------------------------------------
c     do a Runge-Kutta integration of (a*b_theta/r)**2 using the pitch
c     profile to eliminate b_z and using the pressure profile to
c     determine the perpendicular current.  note that during the
c     integration, bth is (a*b_theta/r)**2
c-----------------------------------------------------------------------
      DO i=1,nr
        nsub=1
        sub_sizeb: DO
          dr=(r(i)-r(i-1))/nsub
          IF (i>1) THEN
            bthm=bth(i-1)
          ELSE
            bthm=(be0/pit0)**2
          ENDIF
          rm=r(i-1)
          sub_stepb: DO isub=1,nsub
            rc=rm+0.5*dr
            rp=rm+dr
c-----------------------------------------------------------------------
c           Step 1:
c-----------------------------------------------------------------------
            dbt(1)=-( bthm*(4*rm/r1+dpit2_fun(rm)) + 2*dp_fun(rm) )
     $             /( (rm/r1)**2 + pit_fun(rm)**2 )
c-----------------------------------------------------------------------
c           Step 2:
c-----------------------------------------------------------------------
            bthp=bthm+0.5*dr*dbt(1)/r1
            dbt(2)=-( bthp*(4*rc/r1+dpit2_fun(rc)) + 2*dp_fun(rc) )
     $             /( (rc/r1)**2 + pit_fun(rc)**2 )
c-----------------------------------------------------------------------
c           Step 3:
c-----------------------------------------------------------------------
            bthp=bthm+0.5*dr*dbt(2)/r1
            dbt(3)=-( bthp*(4*rc/r1+dpit2_fun(rc)) + 2*dp_fun(rc) )
     $             /( (rc/r1)**2 + pit_fun(rc)**2 )
c-----------------------------------------------------------------------
c           Step 4:
c-----------------------------------------------------------------------
            bthp=bthm+dr*dbt(3)/r1
            dbt(4)=-( bthp*(4*rp/r1+dpit2_fun(rp)) + 2*dp_fun(rp) )
     $             /( (rp/r1)**2 + pit_fun(rp)**2 )
            bthm=bthm+dr/(6.*r1)*(SUM(dbt(1:4:3))+2*SUM(dbt(2:3)))
            rm=rp
          ENDDO sub_stepb
c-----------------------------------------------------------------------
c         check convergence.
c-----------------------------------------------------------------------
          IF (ABS(bthm-bthl)/(0.5*(ABS(bthm)+ABS(bthl)))<=tol) THEN
            bth(i)=bthm
            EXIT sub_sizeb
          ENDIF
          bthl=bthm
          nsub=2*nsub
          IF (nsub>nsub_max)
     $      CALL nim_stop('Polar_circle_pp is not converging on b.')
        ENDDO sub_sizeb
      ENDDO
c-----------------------------------------------------------------------
c     set the equilibrium fields.
c-----------------------------------------------------------------------
      bth=SQRT(bth)
      DO i=1,nr
        bz(i)=bth(i)*pit_fun(r(i))
        bth(i)=bth(i)*r(i)/r1
      ENDDO
c-----------------------------------------------------------------------
c     current density.
c-----------------------------------------------------------------------
      DO i=1,nr
        bmag=SQRT(bth(i)**2+bz(i)**2)
        jpar=bth(i)**2*(2*r1*pit_fun(r(i))/r(i)
     $                  -0.5*dpit2_fun(r(i))/pit_fun(r(i)))/
     $                 (mu0*r(i)*bmag)
        jt(i)=jpar*bth(i)/bmag+bz(i)*dp_fun(r(i))/(mu0*r1*bmag**2)
        jz(i)=jpar*bz(i)/bmag-bth(i)*dp_fun(r(i))/(mu0*r1*bmag**2)
      ENDDO
      IF (r(0)==0) THEN
        jz0=2*be0/(mu0*r1*pit0)
        IF (geom=='tor') THEN		!   return jth/r
          jt0=(0.5*mu0*jz0**2+
     $         dp_fun(0.001_r8*r(1))/(mu0*xmax*0.001*r(1)))/be0
        ELSE
          jt0=0
        ENDIF
      ELSE
        bmag=be0*SQRT(1+(r(0)/(pit0*r1))**2)
        jpar=(bmag**2-be0**2)*(2*r1*pit_fun(r(0))/r(0)
     $                         -0.5*dpit2_fun(r(0))/pit_fun(r(0)))/
     $                        (mu0*r(0)*bmag)
        jt0=jpar*SQRT(1-(be0/bmag)**2)+be0*dp_fun(r(0))/
     $      (mu0*r1*bmag**2)
        IF (geom=='tor') jt0=jt0/r(0)
        jz0=jpar*be0/bmag-SQRT(1-(be0/bmag)**2)*dp_fun(r(0))/
     $      (mu0*r1*bmag)
      ENDIF
c-----------------------------------------------------------------------
c     compute and write the Suydam parameter and R*q.
c-----------------------------------------------------------------------
      CALL open_bin(xy_unit,"suydam.bin","UNKNOWN","REWIND",32_i4)
      DO i=1,nr
        WRITE(xy_unit) (/REAL(i,4),REAL(r(i),4),REAL(
     $    -2*r(i)*dp_fun(r(i))/xdim/
     $    (2*bz(i)-r(i)*mu0*(jt(i)+jz(i)*bz(i)/bth(i)))**2,4),
     $    REAL(r(i)*bz(i)/bth(i),4)/)
      ENDDO
      CALL close_bin(xy_unit,"suydam.bin")
c-----------------------------------------------------------------------
c     end of executable statements.
c-----------------------------------------------------------------------
      RETURN
c-----------------------------------------------------------------------
c     internal functions used in polar_circle_pp:
c-----------------------------------------------------------------------
      CONTAINS
c-----------------------------------------------------------------------
c       compute pitch at an arbitrary radius.
c-----------------------------------------------------------------------
        FUNCTION pit_fun(rad)

        REAL(r8), INTENT(IN) :: rad
        REAL(r8) :: pit_fun

        IF (lamprof=='qspec') THEN
               pit_fun= pit0*
     $           (  1+  ((rad-r0)/pit2)**(2*pit_4)  )**(1./pit_4)
        ELSE
               pit_fun=pit0*(1+pit2*((rad-r0)/xdim)**2
     $                  +pit_4*((rad-r0)/xdim)**4)
        ENDIF

        END FUNCTION pit_fun
c-----------------------------------------------------------------------
c       compute a*d(pitch**2)/dr.
c-----------------------------------------------------------------------
        FUNCTION dpit2_fun(rad)

        REAL(r8), INTENT(IN) :: rad
        REAL(r8) :: dpit2_fun

        IF (lamprof=='qspec') THEN
          dpit2_fun=4.*pit_fun(rad)* pit0 * r1/pit2
     $       *((rad-r0)/pit2)**(2.*pit_4-1.) 
     $       * (1+ ((rad-r0)/pit2)**(2*pit_4) )**((1.-pit_4)/pit_4)
        ELSE
          dpit2_fun=2*pit_fun(rad)*pit0*r1/xdim
     $             *(2*pit2*((rad-r0)/xdim)
     $              +4*pit_4*((rad-r0)/xdim)**3)
        ENDIF


        END FUNCTION dpit2_fun
c-----------------------------------------------------------------------
c       compute a*mu0*dp/dr.
c-----------------------------------------------------------------------
        FUNCTION dp_fun(rad)

        REAL(r8), INTENT(IN) :: rad
        REAL(r8) :: dp_fun

        IF (lamprof=='qspec') THEN
          dp_fun=be0**2/(2.*pit_fun(rad)**4) * ((rad-r0)/xdim)
     $     * ( ((rad-r0)/xdim)*dpit2_fun(rad) - 4.*pit_fun(rad)**2 )
     $     * kpr
        ELSE
          dp_fun=beta*be0**2*r1/(2*xdim)
     $             *(2*pres_2*((rad-r0)/xdim)
     $              +4*pres_4*((rad-r0)/xdim)**3)
        ENDIF

        END FUNCTION dp_fun
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE polar_circle_pp
c-----------------------------------------------------------------------
c     subprogram 5. polar_fluxgrid_init.
c     reads ascii grid file.
c-----------------------------------------------------------------------
      SUBROUTINE polar_fluxgrid_init
      USE spline
      USE bicube

      INTEGER(i4) :: ix,iy,mx_fldat,my_fldat,ntor,mvac,mpsi,ib
      INTEGER(i4) :: it
      INTEGER :: ios
      INTEGER(i4) :: mxf_rat,myf_rat,mpsi_pol
      LOGICAL :: file_stat
      TYPE(bicube_type) :: bc_xy,bc_be,bc_ja,bc_pr
      TYPE(spline_type) :: q_spline
      REAL(r8) :: x,y,dx,dy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rnode,q_eq

      REAL(r8), DIMENSION(2) :: xxyy
      REAL(r8), DIMENSION(3) :: bvec
      INTEGER(i4) :: icoil
c-----------------------------------------------------------------------
c     read scalar data.
c-----------------------------------------------------------------------
      INQUIRE(FILE=TRIM(eqfile),EXIST=file_stat)
      IF (.NOT.file_stat) CALL nim_stop
     $  ('FLUXGRID output file '//TRIM(eqfile)//' does not exist.')
      OPEN(UNIT=eq_unit,FILE=TRIM(eqfile),STATUS="OLD")
      READ(eq_unit,*)mx_fldat
      READ(eq_unit,*)my_fldat
      READ(eq_unit,*)mvac
      READ(eq_unit,*)mxpie
      READ(eq_unit,*)sqo
      READ(eq_unit,*)xo,yo,bo
      mpsi=mx_fldat-mvac
c-----------------------------------------------------------------------
c     write warning message if input grid parameters are different
c     from what is in fluxgrid.dat.
c-----------------------------------------------------------------------
      IF (mx/=mx_fldat+1) THEN
        WRITE(nim_wr,*) ' '
        WRITE(nim_wr,*) "WARNING: nimrod.in input mx = ",mx," /= grid ",
     $    "mx = ",mx_fldat+1," from fluxgrid.dat."
        WRITE(nim_wr,*) ' '
        WRITE(out_unit,*) ' '
        WRITE(out_unit,*)
     $    "WARNING: nimrod.in input mx = ",mx," /= grid ",
     $    "mx = ",mx_fldat+1," from fluxgrid.dat."
        WRITE(out_unit,*) ' '
      ENDIF
      IF (my/=my_fldat) THEN
        WRITE(nim_wr,*) ' '
        WRITE(nim_wr,*) "WARNING: nimrod.in input my = ",my," /= grid ",
     $    "my = ",my_fldat," from fluxgrid.dat."
        WRITE(nim_wr,*) ' '
        WRITE(out_unit,*) ' '
        WRITE(out_unit,*)
     $    "WARNING: nimrod.in input my = ",my," /= grid ",
     $    "my = ",my_fldat," from fluxgrid.dat."
        WRITE(out_unit,*) ' '
      ENDIF
c-----------------------------------------------------------------------
c     if mx_fldat/=mx, it must be a multiple of mx.
c-----------------------------------------------------------------------
      IF (MODULO(mx_fldat+1_i4,mx)/=0) THEN
        CALL nim_stop
     $    ("Polar_fluxgrid_init: mx_fldat is not a multiple of mx.")
      ELSE
        mxf_rat=(mx_fldat+1)/mx
      ENDIF
c-----------------------------------------------------------------------
c     if my_fldat/=my, it must be a multiple of my.
c-----------------------------------------------------------------------
      IF (MODULO(my_fldat,my)/=0) THEN
        CALL nim_stop
     $    ("Polar_fluxgrid_init: my_fldat is not a multiple of my.")
      ELSE
        myf_rat=my_fldat/my
      ENDIF
c-----------------------------------------------------------------------
c     read equilibrium arrays 
c     In plasma: ix=0:mpsi, conc=0
c     In vacuum: ix=1:mvac, sq() not defined since grid not flux-aligned
c     Note that the "vacuum" here, may actually contain plasma
c     and that conceq is a temporary variable of poly_degree = 1.
c
c     Separate splines for plasma & vacuum regions are possible.
c-----------------------------------------------------------------------
c     read equilibrium arrays with tblock.
c-----------------------------------------------------------------------
      IF(TRIM(pieflag)/='rblock')THEN
         CALL bicube_alloc(bc_xy,mx_fldat,my_fldat,2_i4)
         CALL bicube_alloc(bc_be,mx_fldat,my_fldat,3_i4)
         CALL bicube_alloc(bc_ja,mx_fldat,my_fldat,3_i4)
         CALL bicube_alloc(bc_pr,mx_fldat,my_fldat,1_i4)
         CALL lagr_quad_alloc(conceq,mx_fldat,my_fldat,1_i4,1_i4,
     $                       name='co',title=(/'conceq'/))
         ALLOCATE(sq(0:mpsi,4))
c                                                               Plasma
         conceq%fs=0.
         DO ix=0,mpsi
            READ(eq_unit,*)sq(ix,:)
            DO iy=0,my_fldat
               READ(eq_unit,*)bc_xy%fs(1,ix,iy),bc_xy%fs(2,ix,iy),
     $            bc_be%fs(:,ix,iy),bc_ja%fs(:,ix,iy)
            ENDDO
         ENDDO
c                                                               Vacuum
         DO ix=mpsi+1,mx_fldat
            DO iy=0,my_fldat
               READ(eq_unit,*)bc_xy%fs(1,ix,iy),bc_xy%fs(2,ix,iy),
     $            bc_be%fs(:,ix,iy),bc_ja%fs(:,ix,iy),
     $            bc_pr%fs(1,ix,iy),conceq%fs(1,ix,iy)
            ENDDO
         ENDDO
c-----------------------------------------------------------------------
c     read equilibrium arrays without tblock.
c-----------------------------------------------------------------------
      ELSE
         mx_fldat=mx_fldat+1
         mpsi=mpsi+1
         CALL bicube_alloc(bc_xy,mx_fldat,my_fldat,2_i4)
         CALL bicube_alloc(bc_be,mx_fldat,my_fldat,3_i4)
         CALL bicube_alloc(bc_ja,mx_fldat,my_fldat,3_i4)
         CALL bicube_alloc(bc_pr,mx_fldat,my_fldat,1_i4)
         CALL lagr_quad_alloc(conceq,mx_fldat,my_fldat,1_i4,1_i4,
     $                        name='co',title=(/' conc '/))
         ALLOCATE(sq(0:mpsi,4))
c                                                               Plasma
         conceq%fs=0.
         DO ix=1,mpsi
            READ(eq_unit,*)sq(ix,:)
            DO iy=0,my_fldat
               READ(eq_unit,*)bc_xy%fs(1,ix,iy),bc_xy%fs(2,ix,iy),
     $            bc_be%fs(:,ix,iy),bc_ja%fs(:,ix,iy)
            ENDDO
         ENDDO
c                                                               Vacuum
         DO ix=mpsi+1,mx_fldat
            DO iy=0,my_fldat
               READ(eq_unit,*)bc_xy%fs(1,ix,iy),bc_xy%fs(2,ix,iy),
     $            bc_be%fs(:,ix,iy),bc_ja%fs(:,ix,iy),
     $            bc_pr%fs(1,ix,iy),conceq%fs(1,ix,iy)
            ENDDO
         ENDDO
c                                               Degenerate pt at center
         bc_xy%fs(1,0,:)=xo
         bc_xy%fs(2,0,:)=yo
         bc_be%fs(:,0,:)=0
         bc_be%fs(3,0,:)=bo
         bc_ja%fs(1:2,0,0:my_fldat)=0
         bc_ja%fs(3,0,0:my_fldat)=1/REAL(2*my_fldat)
     $     *SUM(2*bc_ja%fs(3,1,1:my_fldat)
     $     -mxpie**2*(2*bc_ja%fs(3,2,1:my_fldat)
     $                -bc_ja%fs(3,1,1:my_fldat)
     $                -bc_ja%fs(3,3,1:my_fldat))
     $     -mxpie*(4*bc_ja%fs(3,2,1:my_fldat)
     $             -3*bc_ja%fs(3,1,1:my_fldat)
     $               -bc_ja%fs(3,3,1:my_fldat)))
         sq(0,:)=sqo
      ENDIF
c-----------------------------------------------------------------------
c     phi component is covariant.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
         bc_be%fs(3,:,:)=bc_be%fs(3,:,:)*bc_xy%fs(1,:,:)
         bo=bo*xo
      ENDIF
c-----------------------------------------------------------------------
c     set MKS pressure from the normalized flux-surface quantities.
c       p=p_(no_flow)*exp(gamma/2 Ms^2 ((R/Raxis)^2 -1) )
c-----------------------------------------------------------------------
      DO iy=0,my_fldat
        bc_pr%fs(1,0:mpsi,iy)=sq(:,2)/mu0
     $   *EXP(5._r8/6._r8*sq(:,4)**2*((bc_xy%fs(1,0:mpsi,iy)/xo)**2-1))
      ENDDO
      po=sqo(2)/mu0
c-----------------------------------------------------------------------
c     ensure exact periodicity.
c-----------------------------------------------------------------------
      bc_xy%fs(:,:,my_fldat)=bc_xy%fs(:,:,0)
      bc_be%fs(:,:,my_fldat)=bc_be%fs(:,:,0)
      bc_ja%fs(:,:,my_fldat)=bc_ja%fs(:,:,0)
      bc_pr%fs(:,:,my_fldat)=bc_pr%fs(:,:,0)
c-----------------------------------------------------------------------
c     fit the bicubics using R along the radially outward grid-line as
c     the independent variable.
c-----------------------------------------------------------------------
      bc_xy%xs=bc_xy%fs(1,:,0)
      bc_xy%ys=(/(iy,iy=0_i4,my_fldat)/)
      bc_be%xs=bc_xy%xs
      bc_be%ys=bc_xy%ys
      bc_ja%xs=bc_xy%xs
      bc_ja%ys=bc_xy%ys
      bc_pr%xs=bc_xy%xs
      bc_pr%ys=bc_xy%ys
      CALL bicube_fit(bc_xy,"extrap","periodic")
      CALL bicube_fit(bc_be,"extrap","periodic")
      CALL bicube_fit(bc_ja,"extrap","periodic")
      CALL bicube_fit(bc_pr,"extrap","periodic")
c-----------------------------------------------------------------------
c     also fit a spline of the surface quantities with R as the
c     independent variable.
c-----------------------------------------------------------------------
      CALL spline_alloc(q_spline,mpsi,1_i4)
      q_spline%xs=bc_xy%xs(0:mpsi)
      q_spline%fs(:,1)=sq(:,3)
      CALL spline_fit(q_spline,"extrap")
c-----------------------------------------------------------------------
c     apply this module's packing based on q(R_outboard) or
c     R(R_outboard), as specified by packbigr.
c-----------------------------------------------------------------------
      ALLOCATE(rnode(0:mx))
      rnode=bc_xy%xs(0:mx_fldat:mxf_rat)
      IF (MAXVAL(wpack*amp*qpack)>0._r8) THEN
        mpsi_pol=mpsi/mxf_rat
        ALLOCATE(q_eq(mpsi_pol))
        DO ix=1,mpsi_pol
          IF (packbigr) THEN
            q_eq(ix)=0.5_r8*(rnode(ix)+rnode(ix-1))
          ELSE
            CALL spline_eval(q_spline,0.5*(rnode(ix)+rnode(ix-1)),0_i4)
            q_eq(ix)=q_spline%f(1)
          ENDIF
        ENDDO
        CALL polar_circlegrid_pack(rnode(0:mpsi_pol),q_eq)
      ENDIF
c-----------------------------------------------------------------------
c     now transfer the smoothly fitted data to the lagrange_quad
c     types for the equilibrium fields.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lxy,mx,my,2_i4,poly_degree)
      CALL lagr_quad_alloc(lbq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(ljq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(lpq,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(lnd,mx,my,1_i4,poly_degree)
      DO iy=0,my
        DO ix=0,mx
          ib=1
          x=rnode(ix)
          y=iy*myf_rat
          CALL bicube_eval(bc_xy,x,y,0_i4)
          CALL bicube_eval(bc_be,x,y,0_i4)
          CALL bicube_eval(bc_ja,x,y,0_i4)
          CALL bicube_eval(bc_pr,x,y,0_i4)
          CALL lagr_quad_basis_assign_loc(lxy,bc_xy%f,ib,ix,iy)
          CALL lagr_quad_basis_assign_loc(lbq,bc_be%f,ib,ix,iy)
          CALL lagr_quad_basis_assign_loc(ljq,bc_ja%f,ib,ix,iy)
          CALL lagr_quad_basis_assign_loc(lpq,bc_pr%f+pres_offset,
     $                                      ib,ix,iy)
        ENDDO
      ENDDO
      lnd=ndens
c-----------------------------------------------------------------------
c     vertical sides follow uniform spacing in the element coordinate
c     for the poloidal angle.
c-----------------------------------------------------------------------
      IF (poly_degree>1) THEN
        DO iy=1,my
          DO ix=0,mx
            DO ib=poly_degree+1,2*poly_degree-1
              x=rnode(ix)
              y=(iy-1+lxy%dy(ib))*myf_rat
              CALL bicube_eval(bc_xy,x,y,0_i4)
              CALL bicube_eval(bc_be,x,y,0_i4)
              CALL bicube_eval(bc_ja,x,y,0_i4)
              CALL bicube_eval(bc_pr,x,y,0_i4)
              CALL lagr_quad_basis_assign_loc(lxy,bc_xy%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lbq,bc_be%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,bc_ja%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,bc_pr%f+pres_offset,
     $                                          ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       for horizontal sides, divide the independent variable for the
c       bicube fits into poly_degree equal sements, evaluate position
c       and other fields here, and set the lagrange fields to these
c       values.
c-----------------------------------------------------------------------
        DO iy=0,my
          DO ix=1,mx
            DO ib=2,poly_degree
              x=(1-lxy%dx(ib))*rnode(ix-1)+lxy%dx(ib)*rnode(ix)
              y=iy*myf_rat
              CALL bicube_eval(bc_xy,x,y,0_i4)
              CALL bicube_eval(bc_be,x,y,0_i4)
              CALL bicube_eval(bc_ja,x,y,0_i4)
              CALL bicube_eval(bc_pr,x,y,0_i4)
              CALL lagr_quad_basis_assign_loc(lxy,bc_xy%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lbq,bc_be%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,bc_ja%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,bc_pr%f+pres_offset,
     $                                          ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       interiors are treated similarly.
c-----------------------------------------------------------------------
        DO iy=1,my
          DO ix=1,mx
            DO ib=2*poly_degree,poly_degree**2
              x=(1-lxy%dx(ib))*rnode(ix-1)+lxy%dx(ib)*rnode(ix)
              y=(iy-1+lxy%dy(ib))*myf_rat
              CALL bicube_eval(bc_xy,x,y,0_i4)
              CALL bicube_eval(bc_be,x,y,0_i4)
              CALL bicube_eval(bc_ja,x,y,0_i4)
              CALL bicube_eval(bc_pr,x,y,0_i4)
              CALL lagr_quad_basis_assign_loc(lxy,bc_xy%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lbq,bc_be%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,bc_ja%f,ib,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,bc_pr%f+pres_offset,
     $                                          ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     reset the number density profile to satisfy 
c     p/p0=(n/ndens)**gamma_nimset if gamma_nimset>0.
c-----------------------------------------------------------------------
      IF (gamma_nimset>0) THEN
        DO ib=1,SIZE(lnd%ix0)
          DO iy=lnd%iy0(ib),lnd%my
            DO ix=lnd%ix0(ib),lnd%mx
              dx=ix-lnd%ix0(ib)+lnd%dx(ib)
              dy=iy-lnd%iy0(ib)+lnd%dy(ib)
              CALL lagr_quad_eval(lpq,dx,dy,0_i4)
c              tmpr=ndens*(lpq%f(1)/po)**(1._r8/gamma_nimset)
              CALL lagr_quad_basis_assign_loc(lnd,
     $          (/ndens*(lpq%f(1)/po)**(1._r8/gamma_nimset)/),ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     superpose a vacuum magnetic fields due to solenoids.
c-----------------------------------------------------------------------
      IF (ncoil>0) THEN
        bvec(3)=0._r8
        DO ib=1,SIZE(lxy%ix0)
          DO iy=lxy%iy0(ib),lxy%my
            DO ix=lxy%ix0(ib),lxy%mx
              dx=ix-lxy%ix0(ib)+lxy%dx(ib)
              dy=iy-lxy%iy0(ib)+lxy%dy(ib)
              CALL lagr_quad_eval(lxy,dx,dy,1_i4)
              CALL brz_eval(bvec(1:2),lxy%f(1),lxy%f(2),ncoil,
     $                      coil_r,coil_z,coil_current,mu0)
              CALL lagr_quad_basis_add_loc(lbq,bvec,ib,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_fluxgrid_init
c-----------------------------------------------------------------------
c     subprogram 6. polar_rblock_init.
c     initializes rectangular grid blocks.
c-----------------------------------------------------------------------
      SUBROUTINE polar_rblock_init(rb,nrbl)

      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      INTEGER(i4), INTENT(OUT) :: nrbl

      INTEGER(i4) :: ix,iy,ib,ixbl,iybl,mx1,my1,mx2,my2,lx,ly,x0,ibasis
      REAL(r8), DIMENSION(3) :: vec
      REAL(r8) :: eta,ephi,r,xv,dx,dy,rmx,b2,lam,dfs
      INTEGER(i4) :: mx_temp,ios
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: beg,rzg,pg
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: beh,rzh,ph
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: bev,rzv,pv
      REAL(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: bei,rzi,pi
      INTEGER(i4) :: polar_shift
c-----------------------------------------------------------------------
c     find the diffusivity shaping factor and concentration var.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(ldiff,mx,my,1_i4,poly_degree)
      IF (gridshape=='circ') THEN
        CALL lagr_quad_alloc(conceq,mx,my,1_i4,1_i4,
     $                       name='co',title=(/' conc '/))
        conceq%fs=0.
        IF (xvac>0) THEN
          xv=xvac
        ELSE
          xv=xmax
        ENDIF
        DO ibasis=1,SIZE(ldiff%ix0)
          DO iy=ldiff%iy0(ibasis),ldiff%my
            DO ix=ldiff%ix0(ibasis),ldiff%mx
              dx=ix-ldiff%ix0(ibasis)+ldiff%dx(ibasis)
              dy=iy-ldiff%iy0(ibasis)+ldiff%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=SQRT((lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2)
              IF (lamprof=='qspec') THEN		! Holmes const E
                CALL lagr_quad_eval(ljq,dx,dy,0_i4)
                CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/jo/ljq%f(3)/),ibasis,ix,iy)
              ELSE					! std profile
                IF (dvac>=1) THEN
                  dfs=MIN(dvac,(1+(SQRT(dvac)-1)*
     $                ((r-xmin)/(xv-xmin))**dexp)**2)
                ELSE
                  dfs=(1+(SQRT(dvac)-1)*
     $                ((r-xmin)/(xv-xmin))**dexp)**2
                ENDIF
                CALL lagr_quad_basis_assign_loc(ldiff,
     $              (/dfs/),ibasis,ix,iy)
              ENDIF
              IF (ibasis==1.AND.r>xv) conceq%fs(1,ix,iy)=1.
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     For fluxgrid gridshapes, use the radial logical coordinate for
c     computing the diffusivity shape.
c-----------------------------------------------------------------------
      ELSE
        DO ibasis=1,SIZE(ldiff%ix0)
          DO iy=ldiff%iy0(ibasis),ldiff%my
            DO ix=ldiff%ix0(ibasis),ldiff%mx
              dx=ix-ldiff%ix0(ibasis)+ldiff%dx(ibasis)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1+(SQRT(dvac)-1)*(dx/mx)**dexp)**2/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     create the equilibrium flow profile.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lvq,mx,my,3_i4,poly_degree)
      SELECT CASE(eq_flow)
c-----------------------------------------------------------------------
c     no flow--set lvq to 0.
c-----------------------------------------------------------------------
      CASE('none')
        lvq=0
c-----------------------------------------------------------------------
c     uniform axial flow or rigid rotation for toroidal geometry.
c-----------------------------------------------------------------------
      CASE('uniform')
        vo=ve0
        IF (geom=='tor') THEN
          lvq%fs(1:2,:,:)=0
          lvq%fs(3,:,:)=ve0*lxy%fs(1,:,:)/xo
          IF (poly_degree>1) THEN
            lvq%fsh(1:2,:,:,:)=0
            lvq%fsh(3,:,:,:)=ve0*lxy%fsh(1,:,:,:)/xo
            lvq%fsv(1:2,:,:,:)=0
            lvq%fsv(3,:,:,:)=ve0*lxy%fsv(1,:,:,:)/xo
            lvq%fsi(1:2,:,:,:)=0
            lvq%fsi(3,:,:,:)=ve0*lxy%fsi(1,:,:,:)/xo
          ENDIF
        ELSE
          lvq%fs(1:2,:,:)=0
          lvq%fs(3,:,:)=ve0
          IF (poly_degree>1) THEN
            lvq%fsh(1:2,:,:,:)=0
            lvq%fsh(3,:,:,:)=ve0
            lvq%fsv(1:2,:,:,:)=0
            lvq%fsv(3,:,:,:)=ve0
            lvq%fsi(1:2,:,:,:)=0
            lvq%fsi(3,:,:,:)=ve0
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     rigid poloidal rotation for a cylinder.
c-----------------------------------------------------------------------
      CASE('rigid')
        vo=ve0
        vec(3)=0
        DO ibasis=1,SIZE(lvq%ix0)
          DO iy=lvq%iy0(ibasis),lvq%my
            DO ix=lvq%ix0(ibasis),lvq%mx
              dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
              dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
              eta=twopi*dy/REAL(lvq%my)
              CALL lagr_quad_eval(lxy,REAL(mx,r8),dy,0_i4)
              rmx=(lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=SQRT((lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2)
              IF (r>0) THEN
                vec(1)=-ve0*r/rmx*SIN(eta)
                vec(2)= ve0*r/rmx*COS(eta)
              ELSE
                vec(1:2)=0
              ENDIF
              CALL lagr_quad_basis_assign_loc(lvq,vec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     something like Kissick's sheared poloidal rotation profile for ET.
c     thetav is used as rho1, and phiv is used as rho2.
c-----------------------------------------------------------------------
      CASE('sheared')
        vo=ve0
        vec(3)=0
        DO ibasis=1,SIZE(lvq%ix0)
          DO iy=lvq%iy0(ibasis),lvq%my
            DO ix=lvq%ix0(ibasis),lvq%mx
              dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
              dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
              eta=twopi*dy/REAL(lvq%my)
              CALL lagr_quad_eval(lxy,REAL(mx,r8),dy,0_i4)
              rmx=(lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=SQRT((lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2)
              IF (r>0) THEN
                vec(1)=-ve0*SIN(eta)*
     $                 (r/rmx)/thetav*
     $                 EXP(-(r/rmx)/thetav)/EXP(-1._r8)
                vec(2)= ve0*COS(eta)*
     $                 (r/rmx)/thetav*
     $                 EXP(-(r/rmx)/thetav)/EXP(-1._r8)
              ELSE
                vec(1:2)=0
              ENDIF
              CALL lagr_quad_basis_assign_loc(lvq,vec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     uniform case modulated by an approximate gaussian profile along
c     the minor axis rays.
c-----------------------------------------------------------------------
      CASE('gauss')
        vo=ve0
        vec(1:2)=0
        DO ibasis=1,SIZE(lvq%ix0)
          DO iy=lvq%iy0(ibasis),lvq%my
            DO ix=lvq%ix0(ibasis),lvq%mx
              dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
              dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
              CALL lagr_quad_eval(lxy,REAL(mx,r8),dy,0_i4)
              rmx=(lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=(lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2
              vec(3)=ve0*EXP(-r/(eqflow_width**2*rmx))
              IF (geom=='tor') vec(3)=vec(3)*lxy%f(1)/xo 
              CALL lagr_quad_basis_assign_loc(lvq,vec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     mhd pinch flow based on inferred E at the magnetic axis.
c-----------------------------------------------------------------------
      CASE('pinch')
        vo=0
        ephi=mu0*elecd*jo
        vec(3)=0
        IF (geom=='tor') ephi=ephi*xo
        DO ibasis=1,SIZE(lvq%ix0)
          DO iy=lvq%iy0(ibasis),lvq%my
            DO ix=lvq%ix0(ibasis),lvq%mx
              dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
              dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
              eta=mu0*elecd
              IF (ds_use=='equil') THEN
                CALL lagr_quad_eval(ldiff,dx,dy,0_i4)
                eta=eta*ldiff%f(1)
              ENDIF
              CALL lagr_quad_eval(lbq,dx,dy,0_i4)
              CALL lagr_quad_eval(ljq,dx,dy,0_i4)
              IF (geom=='tor') THEN
                CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                lbq%f(3)=lbq%f(3)/lxy%f(1)
                ljq%f(3)=ljq%f(3)*lxy%f(1)
              ENDIF
              vec(1)= (     -eta*ljq%f(2) *lbq%f(3)
     $                -(ephi-eta*ljq%f(3))*lbq%f(2))/SUM(lbq%f**2)
              vec(2)= (      eta*ljq%f(1) *lbq%f(3)
     $                +(ephi-eta*ljq%f(3))*lbq%f(1))/SUM(lbq%f**2)
              CALL lagr_quad_basis_assign_loc(lvq,vec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     ion diamagnetic flow as (1/(1+meomi)-pe_frac)*J_perp/(n_e*e)
c-----------------------------------------------------------------------
      CASE('diamagnetic')
        vo=0
        DO ibasis=1,SIZE(lvq%ix0)
          DO iy=lvq%iy0(ibasis),lvq%my
            DO ix=lvq%ix0(ibasis),lvq%mx
              dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
              dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
              CALL lagr_quad_eval(lbq,dx,dy,0_i4)
              CALL lagr_quad_eval(ljq,dx,dy,0_i4)
              CALL lagr_quad_eval(lnd,dx,dy,0_i4)
              IF (geom=='tor') THEN
                CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                lbq%f(3)=lbq%f(3)/lxy%f(1)
                ljq%f(3)=ljq%f(3)*lxy%f(1)
              ENDIF
              b2=SUM(lbq%f**2)
              lam=SUM(lbq%f*ljq%f)/b2
              vec=(1._r8/(1._r8+meomi)-pe_frac)*(ljq%f-lam*lbq%f)/
     $            (lnd%f(1)*elementary_q)
              CALL lagr_quad_basis_assign_loc(lvq,vec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     Toroidal flow profile read from equilibrium file.
c-----------------------------------------------------------------------
c     Toroidal flow profile read from equilibrium file.
c	Veq(phi) =(R/Raxis)  Ms Vsound =(R/Raxis) Ms (gamma p/rho)^(1/2)
c	In TOQ, I have ndens=1.0 to make things easier.  
c	Perhaps these should be consistent, but then you have to make
c	 sure that your TOQ density matches the NIMROD density and that
c	 seems like a pain since what we really want is a simple 
c	 consistency in normalizations.  There is no easy way of doing 
c	 that since equilibrium codes aren't normally geared for 
c	 giving density information.
c     Note that this should probably be set up so that one can both read
c	in information from files for toroidal flows, and specify a 
c	poloidal flow profile.  Still needs work.
c-----------------------------------------------------------------------
      CASE('eq_file')
        vo=0
c-TMP
c       DO ix=0,mx
c         DO iy=0,my
c            ve_eq%fs(3,ix,iy)=xy%fs(1,ix,iy)/xo * sq(ix,4)
c    $      *(5._r8/6._r8* sq(ix,2)/mu0 /(mtot))**0.5
c         ENDDO
c       ENDDO
      CASE DEFAULT
        CALL nim_stop
     $  ('Equilibrium flow profile '//TRIM(eq_flow)//' not recognized.')
      END SELECT
c-----------------------------------------------------------------------
c     initialize parameters and allocate rblocks.
c-----------------------------------------------------------------------
      nrbl=nxbl*nybl
      ALLOCATE(rb(nrbl))
c-----------------------------------------------------------------------
c     store degeneracy flag for the iterative solver.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
        rb(ib)%degenerate=.false.
      ENDDO
      IF (xmin==0.AND.TRIM(pieflag)=='rblock') THEN
        DO ib=1,nybl
          rb(ib)%degenerate=.true.
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     shift data to mis-align the periodic connections for the mesh
c     and the data.
c-----------------------------------------------------------------------
      ALLOCATE(beg(3,0:mx,0:my),rzg(2,0:mx,0:my),pg(1,0:mx,0:my))
      IF (poly_degree>1) THEN
        ALLOCATE(beh(3,poly_degree-1,1:mx,0:my),
     $           rzh(2,poly_degree-1,1:mx,0:my),
     $            ph(1,poly_degree-1,1:mx,0:my))
        ALLOCATE(bev(3,poly_degree-1,0:mx,1:my),
     $           rzv(2,poly_degree-1,0:mx,1:my),
     $            pv(1,poly_degree-1,0:mx,1:my))
        ALLOCATE(bei(3,(poly_degree-1)**2,1:mx,1:my),
     $           rzi(2,(poly_degree-1)**2,1:mx,1:my),
     $            pi(1,(poly_degree-1)**2,1:mx,1:my))
      ENDIF
      polar_shift=1

      rzg(:,:,polar_shift:my)=lxy%fs(:,:,0:my-polar_shift)
      rzg(:,:,0:polar_shift-1)=lxy%fs(:,:,my-polar_shift:my-1)
      lxy%fs=rzg
      IF (poly_degree>1) THEN
       rzh(:,:,:,polar_shift:my)=lxy%fsh(:,:,:,0:my-polar_shift)
       rzh(:,:,:,0:polar_shift-1)=lxy%fsh(:,:,:,my-polar_shift:my-1)
       rzv(:,:,:,polar_shift+1:my)=lxy%fsv(:,:,:,1:my-polar_shift)
       rzv(:,:,:,1:polar_shift)=lxy%fsv(:,:,:,my-polar_shift+1:my)
       rzi(:,:,:,polar_shift+1:my)=lxy%fsi(:,:,:,1:my-polar_shift)
       rzi(:,:,:,1:polar_shift)=lxy%fsi(:,:,:,my-polar_shift+1:my)
       lxy%fsh=rzh
       lxy%fsv=rzv
       lxy%fsi=rzi
      ENDIF
      beg(:,:,polar_shift:my)=lbq%fs(:,:,0:my-polar_shift)
      beg(:,:,0:polar_shift-1)=lbq%fs(:,:,my-polar_shift:my-1)
      lbq%fs=beg
      IF (poly_degree>1) THEN
        beh(:,:,:,polar_shift:my)=lbq%fsh(:,:,:,0:my-polar_shift)
        beh(:,:,:,0:polar_shift-1)=lbq%fsh(:,:,:,my-polar_shift:my-1)
        bev(:,:,:,polar_shift+1:my)=lbq%fsv(:,:,:,1:my-polar_shift)
        bev(:,:,:,1:polar_shift)=lbq%fsv(:,:,:,my-polar_shift+1:my)
        bei(:,:,:,polar_shift+1:my)=lbq%fsi(:,:,:,1:my-polar_shift)
        bei(:,:,:,1:polar_shift)=lbq%fsi(:,:,:,my-polar_shift+1:my)
        lbq%fsh=beh
        lbq%fsv=bev
        lbq%fsi=bei
      ENDIF
      beg(:,:,polar_shift:my)=ljq%fs(:,:,0:my-polar_shift)
      beg(:,:,0:polar_shift-1)=ljq%fs(:,:,my-polar_shift:my-1)
      ljq%fs=beg
      IF (poly_degree>1) THEN
        beh(:,:,:,polar_shift:my)=ljq%fsh(:,:,:,0:my-polar_shift)
        beh(:,:,:,0:polar_shift-1)=ljq%fsh(:,:,:,my-polar_shift:my-1)
        bev(:,:,:,polar_shift+1:my)=ljq%fsv(:,:,:,1:my-polar_shift)
        bev(:,:,:,1:polar_shift)=ljq%fsv(:,:,:,my-polar_shift+1:my)
        bei(:,:,:,polar_shift+1:my)=ljq%fsi(:,:,:,1:my-polar_shift)
        bei(:,:,:,1:polar_shift)=ljq%fsi(:,:,:,my-polar_shift+1:my)
        ljq%fsh=beh
        ljq%fsv=bev
        ljq%fsi=bei
      ENDIF
      beg(:,:,polar_shift:my)=lvq%fs(:,:,0:my-polar_shift)
      beg(:,:,0:polar_shift-1)=lvq%fs(:,:,my-polar_shift:my-1)
      lvq%fs=beg
      IF (poly_degree>1) THEN
        beh(:,:,:,polar_shift:my)=lvq%fsh(:,:,:,0:my-polar_shift)
        beh(:,:,:,0:polar_shift-1)=lvq%fsh(:,:,:,my-polar_shift:my-1)
        bev(:,:,:,polar_shift+1:my)=lvq%fsv(:,:,:,1:my-polar_shift)
        bev(:,:,:,1:polar_shift)=lvq%fsv(:,:,:,my-polar_shift+1:my)
        bei(:,:,:,polar_shift+1:my)=lvq%fsi(:,:,:,1:my-polar_shift)
        bei(:,:,:,1:polar_shift)=lvq%fsi(:,:,:,my-polar_shift+1:my)
        lvq%fsh=beh
        lvq%fsv=bev
        lvq%fsi=bei
      ENDIF
      pg(:,:,polar_shift:my)=lpq%fs(:,:,0:my-polar_shift)
      pg(:,:,0:polar_shift-1)=lpq%fs(:,:,my-polar_shift:my-1)
      lpq%fs=pg
      IF (poly_degree>1) THEN
        ph(:,:,:,polar_shift:my)=lpq%fsh(:,:,:,0:my-polar_shift)
        ph(:,:,:,0:polar_shift-1)=lpq%fsh(:,:,:,my-polar_shift:my-1)
        pv(:,:,:,polar_shift+1:my)=lpq%fsv(:,:,:,1:my-polar_shift)
        pv(:,:,:,1:polar_shift)=lpq%fsv(:,:,:,my-polar_shift+1:my)
        pi(:,:,:,polar_shift+1:my)=lpq%fsi(:,:,:,1:my-polar_shift)
        pi(:,:,:,1:polar_shift)=lpq%fsi(:,:,:,my-polar_shift+1:my)
        lpq%fsh=ph
        lpq%fsv=pv
        lpq%fsi=pi
      ENDIF
      pg(:,:,polar_shift:my)=lnd%fs(:,:,0:my-polar_shift)
      pg(:,:,0:polar_shift-1)=lnd%fs(:,:,my-polar_shift:my-1)
      lnd%fs=pg
      IF (poly_degree>1) THEN
        ph(:,:,:,polar_shift:my)=lnd%fsh(:,:,:,0:my-polar_shift)
        ph(:,:,:,0:polar_shift-1)=lnd%fsh(:,:,:,my-polar_shift:my-1)
        pv(:,:,:,polar_shift+1:my)=lnd%fsv(:,:,:,1:my-polar_shift)
        pv(:,:,:,1:polar_shift)=lnd%fsv(:,:,:,my-polar_shift+1:my)
        pi(:,:,:,polar_shift+1:my)=lnd%fsi(:,:,:,1:my-polar_shift)
        pi(:,:,:,1:polar_shift)=lnd%fsi(:,:,:,my-polar_shift+1:my)
        lnd%fsh=ph
        lnd%fsv=pv
        lnd%fsi=pi
      ENDIF
      pg(:,:,polar_shift:my)=ldiff%fs(:,:,0:my-polar_shift)
      pg(:,:,0:polar_shift-1)=ldiff%fs(:,:,my-polar_shift:my-1)
      ldiff%fs=pg
      IF (poly_degree>1) THEN
        ph(:,:,:,polar_shift:my)=ldiff%fsh(:,:,:,0:my-polar_shift)
        ph(:,:,:,0:polar_shift-1)=ldiff%fsh(:,:,:,my-polar_shift:my-1)
        pv(:,:,:,polar_shift+1:my)=ldiff%fsv(:,:,:,1:my-polar_shift)
        pv(:,:,:,1:polar_shift)=ldiff%fsv(:,:,:,my-polar_shift+1:my)
        pi(:,:,:,polar_shift+1:my)=ldiff%fsi(:,:,:,1:my-polar_shift)
        pi(:,:,:,1:polar_shift)=ldiff%fsi(:,:,:,my-polar_shift+1:my)
        ldiff%fsh=ph
        ldiff%fsv=pv
        ldiff%fsi=pi
      ENDIF
c-----------------------------------------------------------------------
c     start loop over blocks.
c-----------------------------------------------------------------------
      ib=0
      mx2=0
      DO ixbl=1,nxbl
         mx1=mx2
         mx2=mx*ixbl/nxbl
         lx=mx2-mx1
         my2=0
         DO iybl=1,nybl
            my1=my2
            my2=my*iybl/nybl
            ly=my2-my1
            ib=ib+1
            rb(ib)%mx=lx
            rb(ib)%my=ly
c-----------------------------------------------------------------------
c     allocate and fill coordinates.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%rz,lx,ly,2_i4,poly_degree,
     $        name='lagrrz',title=(/'  r   ','  z   '/))
            rb(ib)%rz%fs=lxy%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%rz%fsh=lxy%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%rz%fsv=lxy%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%rz%fsi=lxy%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill be_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%be_eq,lx,ly,3_i4,poly_degree,
     $        name='lbe_eq',title=(/' be_r ',' be_z ',' be_p '/))
            rb(ib)%be_eq%fs=lbq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%be_eq%fsh=lbq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%be_eq%fsv=lbq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%be_eq%fsi=lbq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill ja_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%ja_eq,lx,ly,3_i4,poly_degree,
     $        name='lja_eq',title=(/' ja_r ',' ja_z ',' ja_p '/))
            rb(ib)%ja_eq%fs=ljq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%ja_eq%fsh=ljq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%ja_eq%fsv=ljq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%ja_eq%fsi=ljq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill ve_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%ve_eq,lx,ly,3_i4,poly_degree,
     $        name='lve_eq',title=(/' ve_r ',' ve_z ',' ve_p '/))
            rb(ib)%ve_eq%fs=lvq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%ve_eq%fsh=lvq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%ve_eq%fsv=lvq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%ve_eq%fsi=lvq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill pres_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%pres_eq,lx,ly,1_i4,poly_degree,
     $        name='lpres',title=(/'lpres '/))
            rb(ib)%pres_eq%fs=lpq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%pres_eq%fsh=lpq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%pres_eq%fsv=lpq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%pres_eq%fsi=lpq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill prese_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%prese_eq,lx,ly,1_i4,
     $                           poly_degree,name='lprese',
     $                           title=(/'lprese'/))
            rb(ib)%prese_eq%fs=pe_frac*lpq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%prese_eq%fsh=pe_frac*
     $          lpq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%prese_eq%fsv=pe_frac*
     $          lpq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%prese_eq%fsi=pe_frac*
     $          lpq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill nd_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%nd_eq,lx,ly,1_i4,poly_degree,
     $        name='lnd_eq',title=(/' lndeq'/))
            rb(ib)%nd_eq%fs=lnd%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%nd_eq%fsh=lnd%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%nd_eq%fsv=lnd%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%nd_eq%fsi=lnd%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill the diffusivity shaping factor.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%diff_shape,lx,ly,1_i4,
     $        poly_degree,name='ldiff',title=(/' ldiff'/))
            rb(ib)%diff_shape%fs=ldiff%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%diff_shape%fsh=ldiff%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%diff_shape%fsv=ldiff%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%diff_shape%fsi=ldiff%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     end loop over blocks.
c-----------------------------------------------------------------------
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_rblock_init
c-----------------------------------------------------------------------
c     subprogram 7. polar_tblock_init.
c     initializes triangular grid blocks.
c-----------------------------------------------------------------------
      SUBROUTINE polar_tblock_init(tb,nrbl)

      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
      INTEGER(i4), INTENT(IN) :: nrbl

      INTEGER(i4) :: ix,iy,ivert,mvert,icell,mcell,kb
c-----------------------------------------------------------------------
c     abort for annular grid.
c-----------------------------------------------------------------------
      IF (gridshape=='circ'.AND.xmin>0) THEN
        IF (pieflag/='rblock') THEN
          CALL nim_stop('xmin must be 0 if pieflag/=rblock.')
        ELSE
          RETURN
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     compute sizes and allocate tblocks.
c-----------------------------------------------------------------------
      mvert=my
      mcell=my
      kb=nrbl+1
      CALL tri_linear_geom_alloc(tb(kb)%tgeom,mvert,mcell)
      tb(kb)%mvert=mvert
      tb(kb)%mcell=mcell
      CALL tri_linear_alloc(tb(kb)%be_eq,mvert,3_i4)
      tb(kb)%be_eq%name='bq'
      tb(kb)%be_eq%title=' be_eq'
      CALL tri_linear_alloc(tb(kb)%ja_eq,mvert,3_i4)
      tb(kb)%ja_eq%name='jq'
      tb(kb)%ja_eq%title=' ja_eq'
      CALL tri_linear_alloc(tb(kb)%ve_eq,mvert,3_i4)
      tb(kb)%ve_eq%name='vq'
      tb(kb)%ve_eq%title=' ve_eq'
      CALL tri_linear_alloc(tb(kb)%pres_eq,mvert,1_i4)
      tb(kb)%pres_eq%name='pq'
      tb(kb)%pres_eq%title=' pr_eq'
      CALL tri_linear_alloc(tb(kb)%prese_eq,mvert,1_i4)
      tb(kb)%prese_eq%name='pq'
      tb(kb)%prese_eq%title='pre_eq'
      CALL tri_linear_alloc(tb(kb)%nd_eq,mvert,1_i4)
      tb(kb)%nd_eq%name='nq'
      tb(kb)%nd_eq%title=' nd_eq'
      CALL tri_linear_alloc(tb(kb)%diff_shape,mvert,1_i4)
      tb(kb)%diff_shape%name='ds'
      tb(kb)%diff_shape%title='dif_sh'
c-----------------------------------------------------------------------
c     define pointers from vertices to neighbors.
c-----------------------------------------------------------------------
      ALLOCATE(tb(kb)%tgeom%neighbor(0)%vertex(0:mvert))
      tb(kb)%tgeom%neighbor(0)%vertex(0)=0
      DO ivert=1,mvert
         tb(kb)%tgeom%neighbor(0)%vertex(ivert)=ivert
      ENDDO
      DO ivert=1,mvert
         ALLOCATE(tb(kb)%tgeom%neighbor(ivert)%vertex(0:3))
         tb(kb)%tgeom%neighbor(ivert)%vertex(0)=ivert
         tb(kb)%tgeom%neighbor(ivert)%vertex(1)=ivert+1
         tb(kb)%tgeom%neighbor(ivert)%vertex(2)=0
         tb(kb)%tgeom%neighbor(ivert)%vertex(3)=ivert-1
      ENDDO
      tb(kb)%tgeom%neighbor(mvert)%vertex(1)=1
      tb(kb)%tgeom%neighbor(1)%vertex(3)=mvert
c-----------------------------------------------------------------------
c     fill coordinates.
c-----------------------------------------------------------------------
      tb(kb)%tgeom%xs(0)=xo
      tb(kb)%tgeom%ys(0)=yo
      ix=0
      iy=0
      DO ivert=1,mvert
         iy=iy+1
         tb(kb)%tgeom%xs(ivert)=lxy%fs(1,ix,iy)
         tb(kb)%tgeom%ys(ivert)=lxy%fs(2,ix,iy)
      ENDDO
c-----------------------------------------------------------------------
c     define pointers from cells to vertices.
c-----------------------------------------------------------------------
      DO icell=1,mcell
         tb(kb)%tgeom%vertex(icell,1)=0
         tb(kb)%tgeom%vertex(icell,2)=icell
         tb(kb)%tgeom%vertex(icell,3)=icell+1
      ENDDO
      tb(kb)%tgeom%vertex(mcell,3)=1
c-----------------------------------------------------------------------
c     fill the equilibrium magnetic field, pressure, density, etc.
c-----------------------------------------------------------------------
      ix=0
      iy=0
      tb(kb)%be_eq%fs=0
      tb(kb)%be_eq%fs(3,0,0)=bo
      tb(kb)%ja_eq%fs=0
      tb(kb)%ja_eq%fs(3,0,0)=jo
      tb(kb)%ve_eq%fs=0
      tb(kb)%ve_eq%fs(3,0,0)=vo
      tb(kb)%pres_eq%fs(1,0,0)=po
      tb(kb)%prese_eq%fs(1,0,0)=pe_frac*po
      DO ivert=1,mvert
         iy=iy+1
         tb(kb)%be_eq%fs(:,ivert,0)=lbq%fs(:,ix,iy)
         tb(kb)%ja_eq%fs(:,ivert,0)=ljq%fs(:,ix,iy)
         tb(kb)%ve_eq%fs(:,ivert,0)=lvq%fs(:,ix,iy)
         tb(kb)%pres_eq%fs(1,ivert,0)=lpq%fs(1,ix,iy)
         tb(kb)%prese_eq%fs(1,ivert,0)=pe_frac*lpq%fs(1,ix,iy)
      ENDDO
      tb(kb)%nd_eq=ndens
      tb(kb)%diff_shape=1
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_tblock_init
c-----------------------------------------------------------------------
c     subprogram 8. polar_seam_init.
c     initializes the grid seams.
c-----------------------------------------------------------------------
      SUBROUTINE polar_seam_init(rb,tb,seam,nrbl,nbl)

      TYPE(rblock_type), DIMENSION(:), POINTER :: rb
      TYPE(tblock_type), DIMENSION(:), POINTER :: tb
      TYPE(edge_type), DIMENSION(:), POINTER :: seam
      INTEGER(i4), INTENT(IN) :: nrbl,nbl

      INTEGER(i4) :: ixbl,iybl,ib,lx,ly,iv,ix,iy,ip,nv,kv,kb,
     $     jvold,jbold,mvert,lb
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: np
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: jb,jv
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ixv,iyv
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 2010 FORMAT('rblock(',i2.2,'), lx = ',i3,', ly = ',i3,', mseam = ',i3)
 2020 FORMAT(/4x,'iv',4x,'ip',4x,'ix',4x,'iy',3x,'outb',2x,'outv'/)
 2030 FORMAT(6i6)
 2040 FORMAT('tblock(',i2.2,'), mvert = ',i3,', mseam = ',i3)
c-----------------------------------------------------------------------
c     start loops over rblocks.
c-----------------------------------------------------------------------
      ib=0
      DO ixbl=1,nxbl
         DO iybl=1,nybl
            ib=ib+1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
            lx=rb(ib)%mx
            ly=rb(ib)%my
            nv=2*lx+2*ly
            seam(ib)%nvert=nv
            ALLOCATE(seam(ib)%vertex(nv))
            ALLOCATE(jb(nv,3))
            ALLOCATE(jv(nv,3))
            ALLOCATE(np(nv))
            ALLOCATE(ixv(2*(lx+ly)))
            ALLOCATE(iyv(2*(lx+ly)))
            np=(/(1,ix=1,lx-1),3,(1,iy=1,ly-1),3,
     $           (1,ix=1,lx-1),3,(1,iy=1,ly-1),3/)
            ixv=(/(ix,ix=1,lx),(lx,iy=1,ly),
     $           (ix,ix=lx-1,0,-1),(0_i4,iy=ly-1,0,-1)/)
            iyv=(/(0_i4,ix=1,lx),(iy,iy=1,ly),
     $           (ly,ix=lx-1,0,-1),(iy,iy=ly-1,0,-1)/)
c-----------------------------------------------------------------------
c     fill block pointers.
c-----------------------------------------------------------------------
            jb=0
            jb(:,1)=(/(ib-1_i4,ix=1,lx),(ib+nybl,iy=1,ly),
     $           (ib+1_i4,ix=1,lx),(ib-nybl,iy=1,ly)/)
            jb(lx,2)=ib-1+nybl
            jb(lx,3)=ib+nybl
            jb(lx+ly,2)=ib+nybl+1
            jb(lx+ly,3)=ib+1
            jb(2*lx+ly,2)=ib+1-nybl
            jb(2*lx+ly,3)=ib-nybl
            jb(2*lx+2*ly,2)=ib-nybl-1
            jb(2*lx+2*ly,3)=ib-1
c-----------------------------------------------------------------------
c     trim left pointers.
c-----------------------------------------------------------------------
            IF(ixbl.EQ.1)THEN
               jb(2*lx+ly+1:2*(lx+ly),1)=0
               jb(2*lx+ly,2:3)=0
               jb(2*(lx+ly),2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim right pointers.
c-----------------------------------------------------------------------
            IF(ixbl.EQ.nxbl)THEN
               jb(lx+1:lx+ly,1)=0
               jb(lx,2:3)=0
               jb(lx+ly,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim bottom pointers.
c-----------------------------------------------------------------------
            IF(iybl.EQ.1)THEN
               jb(1:lx,1)=0
               jb(2*(lx+ly),2:3)=0
               jb(lx,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     trim top pointers.
c-----------------------------------------------------------------------
            IF(iybl.EQ.nybl)THEN
               jb(1+lx+ly:2*lx+ly,1)=0
               jb(lx+ly,2:3)=0
               jb(2*lx+ly,2)=0
            ENDIF
c-----------------------------------------------------------------------
c     fill vertex pointers, non-corners.
c-----------------------------------------------------------------------
            jv=0
            DO iv=1,lx
               IF(jb(iv,1).NE.0)
     $              jv(iv,1)=2*lx+rb(jb(iv,1))%my-iv
               IF(jb(iv+lx+ly,1).NE.0)
     $              jv(iv+lx+ly,1)=lx-iv
            ENDDO
            DO iv=1,ly
               IF(jb(iv+lx,1).NE.0)
     $              jv(iv+lx,1)=2*rb(jb(iv+lx,1))%mx+2*ly-iv
               IF(jb(iv+2*lx+ly,1).NE.0)
     $              jv(iv+2*lx+ly,1)=rb(jb(iv+2*lx+ly,1))%mx+ly-iv
            ENDDO
c-----------------------------------------------------------------------
c     fill vertex pointers, corners.
c-----------------------------------------------------------------------
            IF(jb(lx,2).NE.0)
     $           jv(lx,2)=2*rb(jb(lx,2))%mx+rb(jb(lx,2))%my
            IF(jb(lx,3).NE.0)
     $           jv(lx,3)=2*rb(jb(lx,3))%mx+2*rb(jb(lx,3))%my
            IF(jb(lx+ly,2).NE.0)
     $           jv(lx+ly,2)=2*rb(jb(lx+ly,2))%mx+2*rb(jb(lx+ly,2))%my
            IF(jb(lx+ly,3).NE.0)
     $           jv(lx+ly,3)=rb(jb(lx+ly,3))%mx
            IF(jb(2*lx+ly,1).NE.0)
     $           jv(2*lx+ly,1)=2*(rb(jb(2*lx+ly,1))%mx
     $           +rb(jb(2*lx+ly,1))%my)
            IF(jb(2*lx+ly,2).NE.0)
     $           jv(2*lx+ly,2)=rb(jb(2*lx+ly,2))%mx
            IF(jb(2*lx+ly,3).NE.0)
     $           jv(2*lx+ly,3)=rb(jb(2*lx+ly,3))%mx+rb(jb(2*lx+ly,3))%my
            IF(jb(2*(lx+ly),2).NE.0)
     $           jv(2*(lx+ly),2)=rb(jb(2*(lx+ly),2))%mx
     $           +rb(jb(2*(lx+ly),2))%my
            IF(jb(2*(lx+ly),3).NE.0)
     $           jv(2*(lx+ly),3)=2*rb(jb(2*(lx+ly),3))%mx
     $           +rb(jb(2*(lx+ly),3))%my
c-----------------------------------------------------------------------
c     fill seam pointers.
c-----------------------------------------------------------------------
            DO iv=1,nv
               ALLOCATE(seam(ib)%vertex(iv)%ptr(2,np(iv)))
               seam(ib)%vertex(iv)%intxy=(/ixv(iv),iyv(iv)/)
               DO ip=1,np(iv)
                  seam(ib)%vertex(iv)%ptr(1,ip)=jb(iv,ip)
                  seam(ib)%vertex(iv)%ptr(2,ip)=jv(iv,ip)
               ENDDO
            ENDDO
c-----------------------------------------------------------------------
c     finish loop over blocks.
c-----------------------------------------------------------------------
            DEALLOCATE(np)
            DEALLOCATE(jb)
            DEALLOCATE(jv)
            DEALLOCATE(ixv)
            DEALLOCATE(iyv)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     periodic boundary conditions in theta.
c-----------------------------------------------------------------------
      ib=1
      kb=nybl
      kv=2*rb(kb)%mx+rb(kb)%my
      DO ixbl=1,nxbl
         iv=SIZE(seam(ib)%vertex)
         seam(ib)%vertex(iv)%ptr(1,3)=kb
         seam(ib)%vertex(iv)%ptr(2,3)=kv
         seam(kb)%vertex(kv)%ptr(1,1)=ib
         seam(kb)%vertex(kv)%ptr(2,1)=iv
         DO iv=1,rb(ib)%mx-1
            kv=kv-1
            seam(ib)%vertex(iv)%ptr(1,1)=kb
            seam(ib)%vertex(iv)%ptr(2,1)=kv
            seam(kb)%vertex(kv)%ptr(1,1)=ib
            seam(kb)%vertex(kv)%ptr(2,1)=iv
         ENDDO
         iv=rb(ib)%mx
         kv=kv-1
         seam(ib)%vertex(iv)%ptr(1,1)=kb
         seam(ib)%vertex(iv)%ptr(2,1)=kv
         seam(kb)%vertex(kv)%ptr(1,3)=ib
         seam(kb)%vertex(kv)%ptr(2,3)=iv
         IF(ixbl.lt.nxbl)THEN
            seam(ib+nybl)%vertex(SIZE(seam(ib+nybl)%vertex))%ptr(1,2)=kb
            seam(ib+nybl)%vertex(SIZE(seam(ib+nybl)%vertex))%ptr(2,2)=kv
            seam(kb)%vertex(kv)%ptr(1,2)=ib+nybl
            seam(kb)%vertex(kv)%ptr(2,2)=SIZE(seam(ib+nybl)%vertex)
            kb=kb+nybl
            kv=2*rb(kb)%mx+rb(kb)%my
            seam(ib)%vertex(iv)%ptr(1,2)=kb
            seam(ib)%vertex(iv)%ptr(2,2)=kv
            seam(kb)%vertex(kv)%ptr(1,2)=ib
            seam(kb)%vertex(kv)%ptr(2,2)=iv
            ib=ib+nybl
         ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     connect all last block vertices at r=0 for degenerate
c     rectangles; the nimrod network routine now contains block-internal
c     communication for these blocks.  also leave the the top left
c     corners connected for the edge segment initialization.
c-----------------------------------------------------------------------
      IF(xmin == 0 .AND. TRIM(pieflag)=='rblock' .AND. nybl>1)THEN
         DO ib=1,nybl
            nv=SIZE(seam(ib)%vertex)
            DO iv=2*rb(ib)%mx+rb(ib)%my+1,nv-1
              seam(ib)%vertex(iv)%ptr=0
            ENDDO
            DEALLOCATE(seam(ib)%vertex(nv)%ptr)
            ALLOCATE(seam(ib)%vertex(nv)%ptr(2,nybl))
            DO kb=1,nybl
               kv=SIZE(seam(kb)%vertex)
               IF(kb /= ib)THEN
                  seam(ib)%vertex(nv)%ptr(1,kb)=kb
                  seam(ib)%vertex(nv)%ptr(2,kb)=kv
               ELSE
                  lb=ib-1
                  IF (lb<1) lb=nybl
                  seam(ib)%vertex(nv)%ptr(1,kb)=lb
                  seam(ib)%vertex(nv)%ptr(2,kb)=2*rb(lb)%mx+rb(lb)%my
               ENDIF
            ENDDO
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     diagnose rblock seams before pie.
c-----------------------------------------------------------------------
      IF(detflag)THEN
         IF(nbl>nrbl)WRITE(out_unit,'(a)')"Before pie:"
         DO ib=1,nrbl
            WRITE(out_unit,2010)ib,rb(ib)%mx,
     $               rb(ib)%my,SIZE(seam(ib)%vertex)
            WRITE(out_unit,2020)
            DO iv=1,SIZE(seam(ib)%vertex)
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     $                 seam(ib)%vertex(iv)%intxy,
     $                 seam(ib)%vertex(iv)%ptr(1,ip),
     $                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
         ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     No pie initialization should take place for degenerate rblocks.
c-----------------------------------------------------------------------
      IF(TRIM(pieflag).EQ.'rblock')THEN
        RETURN
      ELSE
c-----------------------------------------------------------------------
c     set up pointers between rblocks and pie, non-corners.
c-----------------------------------------------------------------------
         kb=nrbl+1
         mvert=tb(kb)%mvert
         ALLOCATE(seam(kb)%vertex(mvert))
         seam(kb)%nvert=mvert
         ib=1
         kv=2*(rb(1)%mx+rb(1)%my)
         DO iv=1,mvert
            kv=kv-1
            seam(kb)%vertex(iv)%intxy=(/iv,0_i4/)
            IF(kv.gt.2*rb(ib)%mx+rb(ib)%my)THEN
               ALLOCATE(seam(kb)%vertex(iv)%ptr(2,1))
               seam(kb)%vertex(iv)%ptr(1,1)=ib
               seam(kb)%vertex(iv)%ptr(2,1)=kv
               seam(ib)%vertex(kv)%ptr(1,1)=kb
               seam(ib)%vertex(kv)%ptr(2,1)=iv
c-----------------------------------------------------------------------
c     set up pointers between rblocks and pie, corners, previous rblock.
c-----------------------------------------------------------------------
            ELSE
               ALLOCATE(seam(kb)%vertex(iv)%ptr(2,2))
               seam(kb)%vertex(iv)%ptr(1,1)=ib
               seam(kb)%vertex(iv)%ptr(2,1)=kv
               jbold=seam(ib)%vertex(kv)%ptr(1,1)
               jvold=seam(ib)%vertex(kv)%ptr(2,1)
               DEALLOCATE(seam(ib)%vertex(kv)%ptr)
               ALLOCATE(seam(ib)%vertex(kv)%ptr(2,2))
               seam(ib)%vertex(kv)%ptr(1,1)=jbold
               seam(ib)%vertex(kv)%ptr(2,1)=jvold
               seam(ib)%vertex(kv)%ptr(1,2)=kb
               seam(ib)%vertex(kv)%ptr(2,2)=iv
c-----------------------------------------------------------------------
c     set up pointers between rblocks and pie, corners, next rblock.
c-----------------------------------------------------------------------
               ib=ib+1
               IF(ib.gt.nybl)ib=1
               kv=2*(rb(ib)%mx+rb(ib)%my)
               seam(kb)%vertex(iv)%ptr(1,2)=ib
               seam(kb)%vertex(iv)%ptr(2,2)=kv
               jbold=seam(ib)%vertex(kv)%ptr(1,3)
               jvold=seam(ib)%vertex(kv)%ptr(2,3)
               DEALLOCATE(seam(ib)%vertex(kv)%ptr)
               ALLOCATE(seam(ib)%vertex(kv)%ptr(1,2))
               ALLOCATE(seam(ib)%vertex(kv)%ptr(2,2))
               seam(ib)%vertex(kv)%ptr(1,1)=kb
               seam(ib)%vertex(kv)%ptr(2,1)=iv
               seam(ib)%vertex(kv)%ptr(1,2)=jbold
               seam(ib)%vertex(kv)%ptr(2,2)=jvold
            ENDIF
         ENDDO
c-----------------------------------------------------------------------
c     diagnose seams after pie.
c-----------------------------------------------------------------------
        IF(detflag)THEN
           WRITE(out_unit,'(a)')"After pie:"
           DO ib=1,nrbl
              WRITE(out_unit,2010)ib,rb(ib)%mx,
     $            rb(ib)%my,SIZE(seam(ib)%vertex)
              WRITE(out_unit,2020)
              DO iv=1,SIZE(seam(ib)%vertex)
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     $                 seam(ib)%vertex(iv)%intxy,
     $                 seam(ib)%vertex(iv)%ptr(1,ip),
     $                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
           ENDDO
c-----------------------------------------------------------------------
c          tblock seams after pie
c-----------------------------------------------------------------------
           DO ib=nrbl+1,nrbl+1
            WRITE(out_unit,2040)ib,tb(ib)%mvert,seam(ib)%nvert
            WRITE(out_unit,2020)
            DO iv=1,seam(ib)%nvert
               DO ip=1,SIZE(seam(ib)%vertex(iv)%ptr,2)
                  WRITE(out_unit,2030)iv,ip,
     $                 seam(ib)%vertex(iv)%intxy,
     $                 seam(ib)%vertex(iv)%ptr(1,ip),
     $                 seam(ib)%vertex(iv)%ptr(2,ip)
               ENDDO
            ENDDO
            WRITE(out_unit,2020)
           ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_seam_init
c-----------------------------------------------------------------------
c     subprogram 9. polar_circlegrid_pack.
c     Packs grid around rational surfaces for circlular cross-section,
c	linear geometries.  
c     Note that q must be provided midway between each r-value.
c-----------------------------------------------------------------------
      SUBROUTINE polar_circlegrid_pack(r,q_eq)

      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: r
      REAL(r8), DIMENSION(1:), INTENT(IN) :: q_eq

      INTEGER(i4), DIMENSION(0:16) :: ipack
      INTEGER(i4), DIMENSION(0:mx) :: iclaim
      INTEGER(i4) :: ix,iq,jq,nrat,jx,jmax,istart,iend
      REAL(r8), DIMENSION(0:mx) :: rnew
      REAL(r8) :: rexp, rext, dr0, a, qspread, dr0hr
      REAL(r8) :: factor, err, weight, wnorm
      TYPE(spline_type) :: grid, pack
      CHARACTER(128) :: msg

      INTEGER(i4), PARAMETER :: packres=99999
      INTEGER(i4) :: res 
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rhr,xsi,q_eqhr
c-----------------------------------------------------------------------
c     Define some basic parameters
c     Find rhr (high res r) and q_eqhr by interpolation
c-----------------------------------------------------------------------
      a = r(mx) - r(0)
      dr0 = a/REAL(mx)					! equispaced dr
      nrat=0						! # of q_rat's
      DO iq = 1,npack
        IF (amp(iq)/=0._r8) nrat=nrat+1
      ENDDO
      qspread=MAXVAL(q_eq)-MINVAL(q_eq)+SQRT(TINY(q_eq))
      dr0hr = dr0/res
      IF(MOD(packres,2) == 0) THEN !res must be odd
        res=packres+1
      ELSE 
        res=packres
      ENDIF
      dr0hr = dr0/res

      ALLOCATE(rhr(0:res*mx))
      ALLOCATE(q_eqhr(1:res*mx))
      DO ix = 0,mx-1
        DO jx = 0,res-1
          rhr(ix*res+jx)=r(ix)+(r(ix+1)-r(ix))*jx/res
        ENDDO
      ENDDO
      rhr(res*mx)=r(mx)
      DO ix =1,mx-1
        DO jx =0,res-1
          q_eqhr(ix*res+jx-(res-1)/2)=
     $      q_eq(ix)+(q_eq(ix+1)-q_eq(ix))*jx/res
        ENDDO   
      ENDDO
      DO jx=0,(res-1)/2
        q_eqhr(mx*res+jx-(res-1)/2)=
     $    q_eq(mx)+(q_eq(mx)-q_eq(mx-1))*jx/res
        q_eqhr((res-1)/2-jx+1)=q_eq(1)-(q_eq(2)-q_eq(1))*jx/res
      ENDDO
c-----------------------------------------------------------------------
c     find d(xsi)/dr for every cell center of a uniform mesh, where
c     xsi will be the new logical coordinate as a function of the old
c     logical coordinate, ix.
c-----------------------------------------------------------------------
      ALLOCATE(xsi(0:res*mx))
      xsi(0)=0
      DO ix=1,mx*res
        weight=1
        DO iq=1,nrat
          IF (qpack(iq)>=999) THEN
            weight=weight+amp(iq)*
     $          EXP(-((0.5*(rhr(ix)+rhr(ix-1))-xvac)/(a*wpack(iq)))**2)
          ELSE
            weight=weight+amp(iq)*
     $          EXP(-((q_eqhr(ix)-qpack(iq))/(qspread*wpack(iq)))**2)
          ENDIF
        ENDDO
        xsi(ix)=xsi(ix-1)+weight
      ENDDO
      wnorm=xsi(mx*res)
      xsi=xsi*mx/wnorm

c-----------------------------------------------------------------------
c     map uniform spacing in xsi back to logical coordinate and hence
c     rold.
c-----------------------------------------------------------------------
      rnew(0)=r(0)
      rnew(mx)=r(mx)
      newi: DO ix=1,mx-1
        DO jx=1,mx*res
          IF ((ix-xsi(jx-1))*(ix-xsi(jx))<=0) THEN
            rnew(ix)=rhr(jx-1)+(ix-xsi(jx-1))*(rhr(jx)-rhr(jx-1))/
     $                (xsi(jx)-xsi(jx-1))
            CYCLE newi
          ENDIF
        ENDDO
        WRITE(msg,'(a,i3)') "Polar_circlegrid_pack: not finding new"
     $    //" logical coordinate ",ix
        CALL nim_stop(msg)
      ENDDO newi
c-----------------------------------------------------------------------
c     Diagnostics - use drawpack.in file
c-----------------------------------------------------------------------
      CALL open_bin(temp_unit,'pack.bin','UNKNOWN','REWIND',32_i4)
        WRITE(temp_unit)REAL(0._r4,4),
     $                  REAL(r(0),4),
     $                  REAL(rnew(0),4),
     $                  REAL((rnew(1) - rnew(0))/dr0,4),
     $                  REAL(0.5*(r(0)+r(1)),4),
     $                  REAL(q_eq(1),4)
      DO ix=1,mx
        WRITE(temp_unit)REAL(REAL(ix)/REAL(mx),4),
     $                  REAL(r(ix),4),
     $                  REAL(rnew(ix),4),
     $                  REAL((rnew(ix) - rnew(ix-1))/dr0,4),
     $                  REAL(0.5*(r(ix-1)+r(ix)),4),
     $                  REAL(q_eq(ix),4)
      ENDDO
      CALL close_bin(temp_unit,'pack.bin')
c-----------------------------------------------------------------------
c     change the mesh to rnew.
c-----------------------------------------------------------------------
      r=rnew
      DEALLOCATE(rhr,xsi,q_eqhr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_circlegrid_pack
c-----------------------------------------------------------------------
c     subprogram 10. polar_assign_conc
c     Because the value of the conc variable depends on the equilibrium,
c     it was set up in the temporary variable (conceq) in this module.  
c     Because it is evolved, it is declared in physics_init.f and set=0.
c     After that, we can assign its value here.
c-----------------------------------------------------------------------
      SUBROUTINE polar_assign_conc(nmodes)
      USE fields
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      INTEGER(i4) :: ix,iy,ibl,ix0,ix1,iy0,iy1,ibasis,
     $               max_basisr,max_basist
      REAL(r8) :: dx,dy,dx_data,dy_data
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: temp
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
        IF (ibl>nbl) EXIT
c-----------------------------------------------------------------------
c       select the offsets for the different basis functions in 
c       quadrilateral elements.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          ix1=rb(ibl)%mx
          iy1=rb(ibl)%my
          ix0=rb(ibl)%conc%ix0(ibasis)
          iy0=rb(ibl)%conc%iy0(ibasis)
          dx=rb(ibl)%conc%dx(ibasis)
          dy=rb(ibl)%conc%dy(ibasis)
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
        ALLOCATE(temp(1,ix0:ix1,iy0:iy1,nmodes))
        temp=0.
        DO ix=ix0,ix1
          DO iy=iy0,iy1
            dx_data=ix-ix0+dx
            dy_data=iy-iy0+dy
            CALL lagr_quad_eval(conceq,dx_data,dy_data,0_i4)
            temp(1,ix,iy,1)=conceq%f(1)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       transfer the data to the appropriate basis.
c-----------------------------------------------------------------------
        IF (ibl<=nrbl) THEN
          CALL lagr_quad_basis_assign_arr(rb(ibl)%conc,temp,ibasis)
c-----------------------------------------------------------------------
c-PRE   Assume triangles (except for pie) are in vacuum for now
c-----------------------------------------------------------------------
        ELSE
           IF (TRIM(pieflag)/='rblock' .AND. ibl == nrbl+1) THEN
              tb(ibl)%conc%fs=0.
           ELSE
              tb(ibl)%conc%fs=1.
           ENDIF
        ENDIF
        DEALLOCATE(temp)
        ibasis=ibasis+1
      ENDDO bb_loop 
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE polar_assign_conc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE polar_init
