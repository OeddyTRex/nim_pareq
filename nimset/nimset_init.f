c-----------------------------------------------------------------------
c     file nimset_init.f:  contains non-physics initialization 
c     routines for nimset.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  nimset_init.
c     1.  rect_init.
c     2.  rect_shaped_init.
c     3.  rect_sh_slab.
c     4.  rect_sh_lam.
c     5.  block_init.
c     6.  seam_init.
c     7.  seam0_init.
c     8.  var0_alloc.
c     9.  ptr_alloc.
c     10. edge_match.
c     11. seam0_poly_init
c     12. rect_quad_init.
c-----------------------------------------------------------------------
c     subprogram 0. nimset_init.
c     module declarations.
c-----------------------------------------------------------------------
      MODULE nimset_init
      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE edge_type_mod
      USE fields
      USE seam_storage_mod
      USE input
      USE physdat
      USE polar_init
      IMPLICIT NONE

      TYPE(lagr_quad_2d_type) :: lxy,lbq,ljq,lpq,lvq,lnd,ldiff

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. rect_init.
c     generates grid for testing nimrod code.
c-----------------------------------------------------------------------
      SUBROUTINE rect_init

      INTEGER(i4) :: ix,iy,jx,jy
      REAL(r8) :: agrid,bgrid
      REAL(r8), DIMENSION(0:mx) :: xx
      REAL(r8), DIMENSION(0:my) :: yy
      CHARACTER(64) :: message

      REAL(r8) :: bex,bey,bez,xh,ctb,stb,xcen,yh,jy0,jz0,dx,dy,xv,ndq,
     $            gc1,gc22,p0,dfs
      REAL(r8), DIMENSION(2) :: xxyy
      REAL(r8), DIMENSION(3) :: bvec,jvec
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,by_eq,bz_eq,jy_eq,jz_eq,
     $          pr_eq,ndeq,b_0,tanh_pr,dbdx
      REAL(r8) :: jda,jdb,jdc,jdk,rmaj,bti,tmpr
      REAL(r8), DIMENSION(ngr+poly_degree) :: wg,rg2
      REAL(r8), DIMENSION(0:ngr+poly_degree+1) :: rg1,bzg

      INTEGER(i4) :: mxh,myh,ixm,iym,ibasis,pmx,ng

      INTEGER(i4) :: icoil
c-----------------------------------------------------------------------
c     check domain limits.
c-----------------------------------------------------------------------
      IF (xmax<=xmin) CALL nim_stop("Rect_init: xmax must be > xmin.")
      IF (ymax<=ymin) CALL nim_stop("Rect_init: ymax must be > ymin.")
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lxy,mx,my,2_i4,poly_degree)
c-----------------------------------------------------------------------
c     preliminary computations.
c-----------------------------------------------------------------------
      IF (firstx==0) THEN 
        xx=(/(ix,ix=0,mx)/)*(xmax-xmin)/mx+xmin
        IF (MAXVAL(ABS(amp))/=0._r8) THEN
          CALL polar_circlegrid_pack(xx,0.5_r8*(xx(0:mx-1)+xx(1:mx)))
        ENDIF
      ELSE IF (firstx>=2._r8/mx.OR.firstx<0) THEN
        WRITE(message,*) 'Inappropriate value for firstx:  ',firstx
        CALL nim_stop(message)
      ELSE IF (xmin==-xmax) THEN
        IF (MODULO(mx,2_i4)/=0) CALL nim_stop('mx should be even.')
        mxh=mx/2
        agrid=(1._r8/mxh-firstx)/(0.5_r8*(mxh-1_i4))
        bgrid=firstx-agrid
        xx(mxh)=0._r8
        DO ix=1,mxh
          xx(mxh+ix)=(agrid*ix+bgrid)*(xmax)+xx(mxh+ix-1)
          xx(mxh-ix)=-xx(mxh+ix)
        ENDDO
      ELSE
        agrid=(1._r8/mx-firstx)/(0.5_r8*(mx-1_i4))
        bgrid=firstx-agrid
        xx(0)=xmin
        DO ix=1,mx
          xx(ix)=(agrid*ix+bgrid)*(xmax-xmin)+xx(ix-1)
        ENDDO
      ENDIF
      IF (firsty==0) THEN 
        yy=(/(iy,iy=0,my)/)*(ymax-ymin)/my+ymin
      ELSE
        IF (lam0/=0.OR.lamprof=='pitprs'.OR.ymin==-ymax) THEN
          myh=my/2
        ELSE
          myh=my
        ENDIF
        IF (firsty>=2._r8/myh.OR.firsty<0) THEN
          WRITE(message,*) 'Inappropriate value for firsty:  ',firsty
          CALL nim_stop(message)
        ELSE
          yh=myh*(ymax-ymin)/my
          agrid=(1._r8/myh-firsty)/(0.5_r8*(myh-1_i4))
          bgrid=firsty-agrid
          yy(0)=ymin
          DO iy=1,myh
            yy(iy)=(agrid*iy+bgrid)*yh+yy(iy-1)
          ENDDO
          IF (myh<my) THEN
            yh=ymax-ymin-yh
            myh=my-myh
            agrid=(1._r8/myh-firsty)/(0.5_r8*(myh-1_i4))
            bgrid=firsty-agrid
            yy(my)=ymax
            DO iy=1,myh
              yy(my-iy)=yy(my-iy+1)-(agrid*iy+bgrid)*yh
            ENDDO
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     generate grid.
c-----------------------------------------------------------------------
      DO jx=0,mx
        DO jy=0,my
          lxy%fs(1,jx,jy)=xx(jx)
          lxy%fs(2,jx,jy)=yy(jy)+skew*(lxy%fs(1,jx,jy)-xmin)
        ENDDO
      ENDDO
      IF (poly_degree>1) THEN
        DO ibasis=2,SIZE(lxy%ix0)
          dx=lxy%dx(ibasis)
          dy=lxy%dy(ibasis)
          DO iy=lxy%iy0(ibasis),lxy%my
            iym=MAX(iy-1,0_i4)
            IF (dy==0) iym=iy
            DO ix=lxy%ix0(ibasis),lxy%mx
              ixm=MAX(ix-1,0_i4)
              IF (dx==0) ixm=ix
              xxyy=(1-dy)*((1-dx)*lxy%fs(:,ixm,iym)+dx*lxy%fs(:,ix,iym))
     $               +dy *((1-dx)*lxy%fs(:,ixm,iy )+dx*lxy%fs(:,ix,iy))
              CALL lagr_quad_basis_assign_loc(lxy,xxyy,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     set equilibrium fields.
c
c     the first block is used for very simple equilibria, now including
c     uniformly rotating cylinders with inertia balanced by JthXBz.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lbq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(ljq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(lpq,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(lnd,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(ldiff,mx,my,1_i4,poly_degree)
      ljq=0
      lnd=ndens
      p0=0
      IF (lam0==0.AND.lamprof/='pitprs'.AND.lamprof/='qspec'.AND.
     $    lamprof(1:5)/='gmode'.AND.lamprof/='tanh'.AND.
     $    lamprof/='grbrap'.AND.lamprof/='jardel'.AND.
     $    lamprof/='cos') THEN
        bex=be0*SIN(pi*thetab)*COS(pi*phib)
        bey=be0*SIN(pi*thetab)*SIN(pi*phib)
        bez=be0*COS(pi*thetab)
        IF (xvac>0) THEN
          xv=xvac
        ELSE
          xv=xmax
        ENDIF
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              IF (geom=='tor') THEN
                dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
                dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
                CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                IF (xmin>0) THEN
                  CALL lagr_quad_basis_assign_loc(lbq,
     $              (/bex*xmin/lxy%f(1),bey,bez*xmin/),ibasis,ix,iy)
                ELSE
                  bvec=(/0._r8,bey,0._r8/)
                  jvec=0._r8
                  IF (eq_flow=="uniform".AND.thetav==0.) THEN
                    bvec(2)=
     $                SQRT(bey**2+mu0*ndens*mtot*(ve0*lxy%f(1)/xmax)**2)
                    jvec(3)=
     $                ndens*mtot*(ve0/xmax)**2/bvec(2)   !   Jph/R
                  ENDIF
                  CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
                  CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
                ENDIF
              ELSE
                CALL lagr_quad_basis_assign_loc(lbq,
     $            (/bex,bey,bez/),ibasis,ix,iy)
              ENDIF
              IF (dvac>=1) THEN
                dfs=MIN(dvac,(1+(SQRT(dvac)-1)*
     $              ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2)
              ELSE
                dfs=(1+(SQRT(dvac)-1)*
     $              ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2
              ENDIF
              CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/dfs/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        lpq=beta*be0**2/(2._r8*mu0)
        p0=beta*be0**2/(2._r8*mu0)
c-----------------------------------------------------------------------
c     gravitational instability equilibrium
c-----------------------------------------------------------------------
      ELSE IF (lamprof=='gmode0') THEN
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              ndq=ndens*EXP(-lxy%f(1)/glength)
              CALL lagr_quad_basis_assign_loc(lnd,(/ndq/),ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1.+(SQRT(dvac)-1)*(-lxy%f(1)/xmax)**(2*dexp))**2/),
     $          ibasis,ix,iy)
              ndq=ndq-ndens*EXP(xmax/glength)
              CALL lagr_quad_basis_assign_loc(lpq,
     $          (/beta*be0**2/(2._r8*mu0)-mtot*gravity(1)*glength*ndq/),
     $          ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lbq,
     $          (/0._r8,0._r8,be0/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     hyperbolic tangent equilibrium with asymmetric pressure for
c     sheared slabs.
c-----------------------------------------------------------------------
      ELSE IF (geom=='lin'.AND.lamprof=='tanh') THEN
        pmx=poly_degree*mx
        ALLOCATE(x(0:pmx),by_eq(0:pmx),bz_eq(0:pmx),
     $           jy_eq(0:pmx),jz_eq(0:pmx),pr_eq(0:pmx))
        x(0:pmx:poly_degree)=lxy%fs(1,0:mx,1)
        DO ibasis=1,poly_degree-1
          x(ibasis:pmx-poly_degree+ibasis:poly_degree)=
     $      lxy%fsh(1,ibasis,1:mx,1)
        ENDDO
        DO ix=0,pmx
          pr_eq(ix)=0.5_r8*beta*be0**2/mu0*(1._r8+tanh_pfrac
     $                 *(EXP(x(ix)/tanh_prl)-EXP(-x(ix)/tanh_prl))
     $                 /(EXP(x(ix)/tanh_prl)+EXP(-x(ix)/tanh_prl)))
          by_eq(ix)=be0*SIN(phib)
     $                 *(EXP(x(ix)/tanh_byl)-EXP(-x(ix)/tanh_byl))
     $                 /(EXP(x(ix)/tanh_byl)+EXP(-x(ix)/tanh_byl))
          bz_eq(ix)=
     $      SQRT(be0**2*(1.+beta)-by_eq(ix)**2-2._r8*mu0*pr_eq(ix))
          jz_eq(ix)=4._r8*be0*SIN(phib)/
     $      (mu0*tanh_byl*(EXP(x(ix)/tanh_byl)+EXP(-x(ix)/tanh_byl))**2)
          jy_eq(ix)=2._r8*beta*be0**2*tanh_pfrac/
     $      (mu0*tanh_prl*(EXP(x(ix)/tanh_prl)+EXP(-x(ix)/tanh_prl))**2)
          jy_eq(ix)=(by_eq(ix)*jz_eq(ix)+jy_eq(ix))/bz_eq(ix)
        ENDDO
        ctb=COS(pi*thetab)
        stb=SIN(pi*thetab)
        xh=(xmax-xmin)/2
        xcen=(xmax+xmin)/2
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              IF (ibasis==1) THEN
                ixm=poly_degree*ix
              ELSE IF (ibasis<=poly_degree) THEN
                ixm=poly_degree*(ix-1)+ibasis-1
              ELSE IF (ibasis<2*poly_degree) THEN
                ixm=poly_degree*ix
              ELSE
                ixm=poly_degree*(ix-1)+1+
     $              MODULO(ibasis-2_i4*poly_degree,poly_degree-1_i4)
              ENDIF
              bvec=(/0._r8, by_eq(ixm)*ctb+bz_eq(ixm)*stb,
     $                      bz_eq(ixm)*ctb-by_eq(ixm)*stb/)
              jvec=(/0._r8, jy_eq(ixm)*ctb+jz_eq(ixm)*stb,
     $                      jz_eq(ixm)*ctb-jy_eq(ixm)*stb/)
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/pr_eq(ixm)/),
     $                                        ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1+(SQRT(dvac)-1)*ABS((lxy%f(1)-xcen)/xh)**dexp)**2/),
     $          ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        p0=0.5_r8*beta*be0**2/mu0
        DEALLOCATE(x,by_eq,bz_eq,jy_eq,jz_eq,pr_eq)
c-----------------------------------------------------------------------
c       B_eq   = B_0(x)*[              phib*sin(x/tanh_byl) s^
c                         +sqrt(1-(phib*sin(x/tanh_byl))^2) g^]
c       g^ = cos(pi*thetab) phi^ + sin(pi*thetab) Z^
c       s^ =-sin(pi*thetab) phi^ + cos(pi*thetab) Z^
c The guide field and sheared field amplitudes are controlled by
c the phib parameter and tanh_byl sets the shear-length scale
c If variation is allowed only in Z (i.e. lin_nmodes=1,lin_nmax=0)
c the k_parallel fraction can be set using thetab to rotate the magetic
c field such that a small portion is in the page
c Alternatively, we can keep thetab=0 and control k_par with lin_nmax=1
c and changing per_length
c       B_0(x) = B_0(0)*[1 - beta*tanh_pfrac*TANH(x/tanh_prl) ]^(1/2)
c       p_0(x) = beta*B_0**2/(2*mu0)*[1 + tanh_pfrac*TANH(x/tanh_prl) ]
c       n_0(x) = ndens*[1 + tanh_nfrac*TANH(x/tanh_ndl) ]
c PREVIOUSLY HAD
c       If beta>0 then n(x)=ndens*p_0(x)/p_0(0), else n(x)=ndens
c-----------------------------------------------------------------------
      ELSE IF (geom=='lin'.AND.lamprof=='cos') THEN
        pmx=poly_degree*mx
        ALLOCATE(x(0:pmx),by_eq(0:pmx),bz_eq(0:pmx),
     $           jy_eq(0:pmx),jz_eq(0:pmx),pr_eq(0:pmx),ndeq(0:pmx))
        ALLOCATE(b_0(0:pmx),tanh_pr(0:pmx),dbdx(0:pmx))
        x(0:pmx:poly_degree)=lxy%fs(1,0:mx,1)
        DO ibasis=1,poly_degree-1
          x(ibasis:pmx-poly_degree+ibasis:poly_degree)=
     $      lxy%fsh(1,ibasis,1:mx,1)
        ENDDO
        DO ix=0,pmx
          tanh_pr(ix)=tanh_pfrac*TANH(x(ix)/tanh_prl)
          pr_eq(ix)=0.5_r8*beta*be0**2/mu0*(1._r8+tanh_pr(ix))
          b_0(ix)=be0*SQRT(1-beta*tanh_pr(ix))
c         note here dbdx is actually 1/B_0 dB_0/dx
          dbdx(ix)=-beta*tanh_pfrac/(2*tanh_prl)
     $             *(1-TANH(x(ix)/tanh_prl)**2)/(1-beta*tanh_pr(ix))
          by_eq(ix)=b_0(ix)*phib*SIN(x(ix)/tanh_byl)
          bz_eq(ix)=b_0(ix)*SQRT(1-(phib*SIN(x(ix)/tanh_byl))**2)
          jz_eq(ix)=(by_eq(ix)*dbdx(ix)
     $               +b_0(ix)*phib/tanh_byl*COS(x(ix)/tanh_byl) ) / mu0 
          jy_eq(ix)=bz_eq(ix)/mu0*(
     $        (phib**2*SIN(x(ix)/tanh_byl)*COS(x(ix)/tanh_byl)/tanh_byl)
     $         /(1-(phib*SIN(x(ix)/tanh_byl))**2)
     $        -dbdx(ix) ) 
c         density profile
          ndeq(ix)=ndens*(1+tanh_nfrac*TANH(x(ix)/tanh_ndl))
        ENDDO
        ctb=COS(pi*thetab)
        stb=SIN(pi*thetab)
        xh=(xmax-xmin)/2
        xcen=(xmax+xmin)/2
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              IF (ibasis==1) THEN
                ixm=poly_degree*ix
              ELSE IF (ibasis<=poly_degree) THEN
                ixm=poly_degree*(ix-1)+ibasis-1
              ELSE IF (ibasis<2*poly_degree) THEN
                ixm=poly_degree*ix
              ELSE
                ixm=poly_degree*(ix-1)+1+
     $              MODULO(ibasis-2_i4*poly_degree,poly_degree-1_i4)
              ENDIF
              bvec=(/0._r8, by_eq(ixm)*ctb+bz_eq(ixm)*stb,
     $                      bz_eq(ixm)*ctb-by_eq(ixm)*stb/)
              jvec=(/0._r8, jy_eq(ixm)*ctb+jz_eq(ixm)*stb,
     $                      jz_eq(ixm)*ctb-jy_eq(ixm)*stb/)
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/pr_eq(ixm)/),
     $                                        ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lnd,(/ndeq(ixm)/),
     $                                        ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1+(SQRT(dvac)-1)*ABS((lxy%f(1)-xcen)/xh)**dexp)**2/),
     $          ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        DEALLOCATE(x,by_eq,bz_eq,jy_eq,jz_eq,pr_eq,ndeq)
        DEALLOCATE(b_0,tanh_pr,dbdx)
c-----------------------------------------------------------------------
c     sheared slab model based on specified current profiles.
c-----------------------------------------------------------------------
      ELSE IF (geom=='lin'.AND.lamprof/='pitprs'.AND.lamprof/='qspec')
     $  THEN
        IF (MODULO(mx,2_i4)/=0) CALL nim_stop
     $    ('Sheared slab model needs mx to be even.')
        pmx=poly_degree*mx
        mxh=pmx/2
        xh=(xmax-xmin)/2
        xcen=(xmax+xmin)/2
        ALLOCATE(x(0:mxh),by_eq(0:mxh),bz_eq(0:mxh),
     $           jy_eq(0:mxh),jz_eq(0:mxh),pr_eq(0:mxh))
        x(0:mxh:poly_degree)=lxy%fs(1,mx/2:mx,1)-xcen
        DO ibasis=1,poly_degree-1
          x(ibasis:mxh-poly_degree+ibasis:poly_degree)=
     $      lxy%fsh(1,ibasis,mx/2+1:mx,1)-xcen
        ENDDO
        bz_eq(0)=be0
        by_eq(0)=0
        jz_eq(0)=lam0*be0/(mu0*xh)
        jy_eq(0)=0
        CALL rect_sh_slab(x,by_eq,bz_eq,jy_eq,jz_eq,pr_eq)
        ctb=COS(pi*thetab)
        stb=SIN(pi*thetab)
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              IF (ibasis==1) THEN
                ixm=ABS(ix-mx/2_i4)
                ixm=poly_degree*ixm
              ELSE IF (ibasis<=poly_degree) THEN
                ixm=ix-mx/2_i4
                IF (ixm<=0) THEN
                  ixm=-ixm+1_i4
                  ixm=poly_degree*ixm-ibasis+1_i4
                ELSE
                  ixm=poly_degree*(ixm-1_i4)+ibasis-1_i4
                ENDIF
              ELSE IF (ibasis<2_i4*poly_degree) THEN
                ixm=ABS(ix-mx/2_i4)
                ixm=poly_degree*ixm
              ELSE
                ixm=ix-mx/2_i4
                IF (ixm<=0_i4) THEN
                  ixm=-ixm+1_i4
                  ixm=poly_degree*ixm-1_i4-
     $                MODULO(ibasis-2_i4*poly_degree,poly_degree-1_i4)
                ELSE
                  ixm=poly_degree*(ixm-1_i4)+1_i4+
     $                MODULO(ibasis-2_i4*poly_degree,poly_degree-1_i4)
                ENDIF
              ENDIF
              IF (lxy%f(1)-xcen<0) THEN
                bvec=(/0._r8,-by_eq(ixm)*ctb+bz_eq(ixm)*stb,
     $                        bz_eq(ixm)*ctb+by_eq(ixm)*stb/)
                jvec=(/0._r8,-jy_eq(ixm)*ctb+jz_eq(ixm)*stb,
     $                        jz_eq(ixm)*ctb+jy_eq(ixm)*stb/)
              ELSE
                bvec=(/0._r8, by_eq(ixm)*ctb+bz_eq(ixm)*stb,
     $                        bz_eq(ixm)*ctb-by_eq(ixm)*stb/)
                jvec=(/0._r8, jy_eq(ixm)*ctb+jz_eq(ixm)*stb,
     $                        jz_eq(ixm)*ctb-jy_eq(ixm)*stb/)
              ENDIF
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/pr_eq(ixm)/),
     $                                        ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1+(SQRT(dvac)-1)*ABS((lxy%f(1)-xcen)/xh)**dexp)**2/),
     $          ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        p0=pr_eq(0)
        DEALLOCATE(x,by_eq,bz_eq,jy_eq,jz_eq,pr_eq)
c-----------------------------------------------------------------------
c     cylinder or annulus centered on the coordinate-system z-axis.
c     for xmin>0, there is no diamagnetic current at xmin due to the
c     form of p(r-xmin).
c
c     the 'grbrap' profile is specified by simple rational functions.
c     the two constants are specified with pit_0 (the normalized
c     pitch on axis) and pit_2.
c-----------------------------------------------------------------------
      ELSE IF (geom=='tor'.AND.lamprof=='grbrap') THEN
        gc1=1._r8/((xmax-xmin)*pit_0)
        gc22=pit_2**2
        bvec=0._r8
        jvec=0._r8
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lbq%ix0(ibasis),lbq%mx
              dx=ix-lbq%ix0(ibasis)+lbq%dx(ibasis)
              dy=iy-lbq%iy0(ibasis)+lbq%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              bvec(2)=be0
              bvec(3)=be0*gc1*lxy%f(1)**2/(1._r8+gc22*lxy%f(1)**2)
              jvec(2)=-2._r8*be0*gc1/(1._r8+gc22*lxy%f(1)**2)**2/mu0
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/
     $          0.5_r8*be0**2*gc1**2/(mu0*gc22)*
     $          (1._r8/(1._r8+gc22*lxy%f(1)**2)**2-
     $           1._r8/(1._r8+gc22)**2)/),ibasis,ix,iy)
              IF (dvac>=1) THEN
                dfs=MIN(dvac,(1+(SQRT(dvac)-1)*
     $              ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2)
              ELSE
                dfs=(1+(SQRT(dvac)-1)*
     $              ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2
              ENDIF
              CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/dfs/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
        p0=lpq%fs(1,0,0)
c-----------------------------------------------------------------------
c     equilibria that result from ODE solutions.  the mesh is assumed to
c     be uniform in the periodic direction (z).
c-----------------------------------------------------------------------
      ELSE  !  geom='tor'
        pmx=poly_degree*mx
        ALLOCATE(x(0:pmx),by_eq(0:pmx),bz_eq(0:pmx),
     $           jy_eq(0:pmx),jz_eq(0:pmx),pr_eq(0:pmx))
        x(0:pmx:poly_degree)=lxy%fs(1,:,1)
        DO ibasis=1,poly_degree-1
          x(ibasis:pmx-poly_degree+ibasis:poly_degree)=
     $      lxy%fsh(1,ibasis,:,1)
        ENDDO
        IF (geom/='tor') CALL nim_stop
     $    ('The specified equilibrium needs toroidal geometry for a '//
     $     'rectangular gridshape.')
        IF (lamprof=='pitprs'.OR.lamprof=='qspec') THEN
          CALL polar_circle_pp(x,bz_eq(1:pmx),by_eq(1:pmx),
     $                         pr_eq(1:pmx),jz_eq(1:pmx),
     $                         jy_eq(1:pmx),jz_eq(0),jy_eq(0),
     $                         xmin,xmax,geom)
          by_eq(0)=be0
          bz_eq(0)=x(0)*be0/(pit_0*xmax)
          pr_eq(0)=beta*be0**2/(2._r8*mu0)
        ELSE IF (lamprof=='jardel') THEN
          ng=ngr+poly_degree
          rmaj=(ymax-ymin)/twopi
          jda=(2._r8-0.5_r8*alpha)*pit_0**2*rmaj**2/(1._r8-pit_0**2)
          jdb=(pit_0**2/rmaj**2)/(1._r8-2._r8*pit_0**2)
          jdc=(pit_0**2*rmaj**2)/(1._r8-2._r8*pit_0**2)
          jdk=(1._r8-pit_0**2)/(1._r8-2._r8*pit_0**2)
          by_eq(0)=be0
          bz_eq(0)=0._r8
          jy_eq(0)=2._r8*be0/(mu0*pit_0)
          jz_eq(0)=2._r8*be0/(mu0*pit_0**2)
          pr_eq(0)=0._r8
          bti=0._r8
          rg1(0)=0._r8
          DO ix=1,pmx    !   first gauleg is to find grid for 2nd int
            CALL gauleg(x(ix-1)/rmaj,x(ix)/rmaj,rg1(1:ng),wg,ng)
            rg1(ng+1)=x(ix)/rmaj
            DO iy=1,ng+1
              CALL gauleg(rg1(iy-1),rg1(iy),rg2,wg,ng)
              bti=bti+jdk*SUM(wg*(2._r8*rg2+jda*rg2**3)/
     $                          (jdb+rg2**2+jdc*rg2**4))
              bzg(iy)=be0*(1._r8-(rmaj*rg1(iy))**2)*exp(-bti)
            ENDDO
            bz_eq(ix)=be0*x(ix)*exp(-bti)/pit_0
            by_eq(ix)=bzg(ng+1)
            tmpr=x(ix)/rmaj
            tmpr=jdk*((2._r8*tmpr+jda*tmpr**3)/
     $               (jdb+tmpr**2+jdc*tmpr**4))/rmaj
            jz_eq(ix)=pit_0*bz_eq(ix)/(mu0*x(ix))*
     $                ((1._r8-x(ix)**2)*tmpr+2._r8*x(ix))
            jy_eq(ix)=(bz_eq(ix)/mu0)*(2._r8/x(ix)-tmpr)
            CALL gauleg(x(ix-1)/rmaj,x(ix)/rmaj,rg1(1:ng),wg,ng)
            pr_eq(ix)=pr_eq(ix-1)-0.5_r8*alpha*rmaj**4*
     $        SUM(wg*rg1(1:ng)*(rg1(1:ng)*bzg(1:ng)/
     $                         (1._r8-(rmaj*rg1(1:ng))**2))**2)
            rg1(0)=rg1(ng+1)
          ENDDO
          pr_eq=(pr_eq+0.5_r8*beta*be0**2)/mu0
        ELSE
          CALL polar_circle_b0(x,bz_eq(1:pmx),by_eq(1:pmx),
     $                         pr_eq(1:pmx),jz_eq(1:pmx),
     $                         jy_eq(1:pmx),xmin,xmax)
          IF (x(0)>0) THEN
            by_eq(0)=be0*COS(pi*thetab)
            bz_eq(0)=be0*SIN(pi*thetab)
            jy_eq(0)=polar_circle_lam(x(0),bz_eq(0),by_eq(0),xmin,xmax)
     $                 *by_eq(0)/mu0
            jz_eq(0)=polar_circle_lam(x(0),bz_eq(0),by_eq(0),xmin,xmax)
     $                 *bz_eq(0)/(x(0)*mu0)
          ELSE
            by_eq(0)=be0
            bz_eq(0)=0
            jy_eq(0)=polar_circle_lam(0._r8,0._r8,be0,xmin,xmax)*be0/mu0
            jz_eq(0)=0.5*polar_circle_lam(0._r8,0._r8,be0,xmin,xmax)**2*
     $                   be0/mu0+pres_2*beta*be0/(xmax**2*mu0)
          ENDIF
          pr_eq(0)=beta*be0**2/(2._r8*mu0)
        ENDIF
        bz_eq( :)=bz_eq( :)*x
        jz_eq(1:)=jz_eq(1:)/x(1:)
        IF (xvac>0) THEN
          xv=xvac
        ELSE
          xv=xmax
        ENDIF
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              IF (ibasis==1) THEN
                ixm=poly_degree*ix
              ELSE IF (ibasis<=poly_degree) THEN
                ixm=poly_degree*(ix-1)+ibasis-1
              ELSE IF (ibasis<2*poly_degree) THEN
                ixm=poly_degree*ix
              ELSE
                ixm=poly_degree*(ix-1)+1+
     $              MODULO(ibasis-2_i4*poly_degree,poly_degree-1_i4)
              ENDIF
              bvec=(/0._r8,by_eq(ixm),-bz_eq(ixm)/)
              jvec=(/0._r8,jy_eq(ixm),-jz_eq(ixm)/)
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/pr_eq(ixm)/),
     $                                        ibasis,ix,iy)
              IF (lamprof=='qspec') THEN	!  Holmes uniform E	
                CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/jy_eq(0)/jy_eq(ixm)/),ibasis,ix,iy)
              ELSE				!  standard profile
                IF (dvac>=1) THEN
                  dfs=MIN(dvac,(1+(SQRT(dvac)-1)*
     $                ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2)
                ELSE
                  dfs=(1+(SQRT(dvac)-1)*
     $                ((lxy%f(1)-xmin)/(xv-xmin))**dexp)**2
                ENDIF
                CALL lagr_quad_basis_assign_loc(ldiff,
     $              (/dfs/),ibasis,ix,iy)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        p0=pr_eq(0)
      ENDIF
c-----------------------------------------------------------------------
c     add pressure offset.
c-----------------------------------------------------------------------
      IF (pres_offset/=0.AND.p0>0) THEN
        lpq%fs=lpq%fs+pres_offset*p0
        IF (poly_degree>1) THEN
          lpq%fsh=lpq%fsh+pres_offset*p0
          lpq%fsv=lpq%fsv+pres_offset*p0
          lpq%fsi=lpq%fsi+pres_offset*p0
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     reset the number density profile to satisfy 
c     p/p0=(n/ndens)**gamma_nimset if gamma_nimset>0.
c-----------------------------------------------------------------------
      IF (gamma_nimset>0.AND.p0>0) THEN
        DO ibasis=1,SIZE(lnd%ix0)
          DO iy=lnd%iy0(ibasis),lnd%my
            DO ix=lnd%ix0(ibasis),lnd%mx
              dx=ix-lnd%ix0(ibasis)+lnd%dx(ibasis)
              dy=iy-lnd%iy0(ibasis)+lnd%dy(ibasis)
              CALL lagr_quad_eval(lpq,dx,dy,0_i4)
              ndq=ndens*(lpq%f(1)/p0)**(1._r8/gamma_nimset)
              CALL lagr_quad_basis_assign_loc(lnd,(/ndq/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find vacuum magnetic field components due to solenoids.
c-----------------------------------------------------------------------
      IF (ncoil>0.AND.geom=='tor') THEN
        bvec(3)=0._r8
        DO ibasis=1,SIZE(lxy%ix0)
          DO iy=lxy%iy0(ibasis),lxy%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,1_i4)
              CALL brz_eval(bvec(1:2),lxy%f(1),lxy%f(2),ncoil,
     $                      coil_r,coil_z,coil_current,mu0)
              CALL lagr_quad_basis_add_loc(lbq,bvec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rect_init
c-----------------------------------------------------------------------
c     subprogram 2. rect_shaped_init.
c     generates a logically rectangular but physically elliptical mesh.
c-----------------------------------------------------------------------
      SUBROUTINE rect_shaped_init

      INTEGER(i4) :: ix,iy,jx,jy,pmx,pmy
      REAL(r8) :: xdim,ydim,xc,yc,r,theta,dr
      REAL(r8), DIMENSION(2,0:poly_degree*mx,0:poly_degree*my,2) :: xy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,bth_eq,bz_eq,pr_eq,
     $          jth_eq,jz_eq
      REAL(r8), DIMENSION(2) :: xxyy
      REAL(r8), DIMENSION(3) :: bvec,jvec
      REAL(r8), DIMENSION(0:poly_degree) :: al,dal

      REAL(r8) :: bex,bey,bez,jy0,jz0,pr0,dx,dy,xv,gc1,gc22
      INTEGER(i4) :: ibasis,ixm,iym,irl
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
      IF (xmax<=xmin) CALL nim_stop("Rect_init: xmax must be > xmin.")
      IF (ymax<=ymin) CALL nim_stop("Rect_init: ymax must be > ymin.")
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lxy,mx,my,2_i4,poly_degree)
c-----------------------------------------------------------------------
c     generate initial grid.
c-----------------------------------------------------------------------
      xc=0.5*(xmin+xmax)
      yc=0.5*(ymin+ymax)
      xdim=0.5*(xmax-xmin)
      ydim=0.5*(ymax-ymin)
      DO jx=0,mx
        DO jy=0,my
          xy(:,jx,jy,1)=(/xc+xdim*SIN(pi*(jx-0.5*mx)/(mx+my)),
     $                    yc+ydim*SIN(pi*(jy-0.5*my)/(mx+my))/)
        ENDDO
      ENDDO
      DO jx=1,mx-1
        xy(2,jx, 0,1)=yc-ydim*COS(pi*(jx-0.5*mx)/(mx+my))
        xy(2,jx,my,1)=yc+ydim*COS(pi*(jx-0.5*mx)/(mx+my))
      ENDDO
      DO jy=1,my-1
        xy(1, 0,jy,1)=xc-xdim*COS(pi*(jy-0.5*my)/(mx+my))
        xy(1,mx,jy,1)=xc+xdim*COS(pi*(jy-0.5*my)/(mx+my))
      ENDDO
c-----------------------------------------------------------------------
c     relax mesh
c-----------------------------------------------------------------------
      DO irl=1,npc
        xy(:,:,:,2)=xy(:,:,:,1)
        DO jy=1,my-1
          DO jx=1,mx-1
            xy(:,jx,jy,1)=SUM(SUM(xy(:,jx-1:jx+1,jy-1:jy+1,1),2),2)/9
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     interpolate locations of side and interior nodes.
c-----------------------------------------------------------------------
      DO ibasis=1,SIZE(lxy%ix0)
        dx=lxy%dx(ibasis)
        dy=lxy%dy(ibasis)
        DO iy=lxy%iy0(ibasis),lxy%my
          iym=MAX(iy-1,0_i4)
          IF (dy==0) iym=iy
          DO ix=lxy%ix0(ibasis),lxy%mx
            ixm=MAX(ix-1,0_i4)
            IF (dx==0) ixm=ix
            xxyy=(1-dy)*((1-dx)*xy(:,ixm,iym,1)+dx*xy(:,ix,iym,1))
     $             +dy *((1-dx)*xy(:,ixm,iy ,1)+dx*xy(:,ix,iy,1))
            CALL lagr_quad_basis_assign_loc(lxy,xxyy,ibasis,ix,iy)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set equilibrium fields.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lbq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(ljq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(lpq,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(lnd,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(ldiff,mx,my,1_i4,poly_degree)
      ljq=0
      lnd=ndens
      IF (lam0==0.AND.lamprof/='pitprs'.AND.lamprof/='qspec'.AND.
     $    lamprof/='grbrap'.OR.geom=='tor') THEN
        IF (thetab/=0) CALL nim_stop('For gridshape=rect_cir'//
     $                               ' thetab must be 0.')
        bex=0
        bey=0
        bez=be0
        DO ibasis=1,SIZE(lbq%ix0)
          DO iy=lbq%iy0(ibasis),lbq%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              IF (geom=='tor') THEN
                IF (xmin>0) THEN
                  dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
                  dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
                  CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                  CALL lagr_quad_basis_assign_loc(lbq,
     $              (/bex*xc/lxy%f(1),bey,bez*xc/),ibasis,ix,iy)
                ELSE
                  CALL lagr_quad_basis_assign_loc(lbq,
     $              (/0._r8,bey,0._r8/),ibasis,ix,iy)
                ENDIF
              ELSE
                CALL lagr_quad_basis_assign_loc(lbq,
     $            (/bex,bey,bez/),ibasis,ix,iy)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        lpq=beta*be0**2/(2._r8*mu0)
c-----------------------------------------------------------------------
c       diffusivity profile
c-----------------------------------------------------------------------
        DO ibasis=1,SIZE(ldiff%ix0)
          DO iy=ldiff%iy0(ibasis),ldiff%my
            DO ix=ldiff%ix0(ibasis),ldiff%mx
              dx=ix-ldiff%ix0(ibasis)+ldiff%dx(ibasis)
              dy=iy-ldiff%iy0(ibasis)+ldiff%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=SQRT((lxy%f(1)-xc)**2+(lxy%f(2)-yc)**2)
              IF (r>SQRT(TINY(r))) THEN
                theta=ATAN2(ABS(lxy%f(2)-yc),ABS(lxy%f(1)-xc))
              ELSE
                theta=0
              ENDIF
              r=r*SQRT((xdim*SIN(theta))**2+(ydim*COS(theta))**2)/
     $          (xdim*ydim)
              CALL lagr_quad_basis_assign_loc(ldiff,
     $          (/(1+(SQRT(dvac)-1)*r**dexp)**2/),ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     the following allows non-trivial equilibria as long as the cross-
c     section is circular.  
c-----------------------------------------------------------------------
      ELSE IF (ABS(xdim-ydim)<1.e-10*xdim.AND.geom=='lin') THEN
        pmx=poly_degree*mx
        ALLOCATE(x(0:pmx),bth_eq(0:pmx),bz_eq(0:pmx),
     $           jth_eq(0:pmx),jz_eq(0:pmx),pr_eq(0:pmx))
        x(0)=0._r8
        DO ix=poly_degree,pmx,poly_degree
          x(ix)=ix*(xdim/pmx)
          DO jx=2,poly_degree
            dx=lxy%dx(jx)
            x(ix-poly_degree-1+jx)=(1._r8-dx)*x(ix-poly_degree)+dx*x(ix)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       pitch and pressure profiles:
c-----------------------------------------------------------------------
        IF (lamprof=='pitprs'.OR.lamprof=='qspec') THEN
          CALL polar_circle_pp(x,bth_eq(1:pmx),bz_eq(1:pmx),
     $                         pr_eq(1:pmx),jth_eq(1:pmx),
     $                         jz_eq(1:pmx),jth_eq(0),jz_eq(0),
     $                         0._r8,xdim,geom)
          bth_eq(0)=0._r8
          bz_eq(0)=be0
          jth_eq(0)=0._r8
          pr_eq(0)=beta*be0**2/(2._r8*mu0)
c-----------------------------------------------------------------------
c       the 'grbrap' profile is specified by simple rational
c       functions.  the two constants are specified with pit_0 (the
c       normalized pitch on axis) and pit_2.
c-----------------------------------------------------------------------
        ELSE IF (lamprof=='grbrap') THEN
          gc1=1._r8/(xdim*pit_0)
          gc22=pit_2**2
          bth_eq=be0*gc1*x/(1._r8+gc22*x**2)
          bz_eq=be0
          jth_eq=0._r8
          jz_eq=2._r8*be0*gc1/(1._r8+gc22*x**2)**2/mu0
          pr_eq=0.5_r8*be0**2*gc1**2/(mu0*gc22)*
     $          (1._r8/(1._r8+gc22*x**2)**2-1._r8/(1._r8+gc22)**2)
c-----------------------------------------------------------------------
c       cases with specified parallel current profiles
c-----------------------------------------------------------------------
        ELSE
          CALL polar_circle_b0(x,bth_eq(1:pmx),bz_eq(1:pmx),
     $                         pr_eq(1:pmx),jth_eq(1:pmx),
     $                         jz_eq(1:pmx),0._r8,xdim)
          bth_eq(0)=0._r8
          bz_eq(0)=be0
          jth_eq(0)=0._r8
          jz_eq(0)=polar_circle_lam(0._r8,0._r8,be0,0._r8,xdim)*be0/mu0
          pr_eq(0)=beta*be0**2/(2._r8*mu0)
        ENDIF
        IF (xvac>0) THEN
          xv=xvac
        ELSE
          xv=xdim
        ENDIF
c-----------------------------------------------------------------------
c       compute and write the Suydam parameter and R*q.
c-----------------------------------------------------------------------
        CALL open_bin(xy_unit,"suydam.bin","UNKNOWN","REWIND",32_i4)
        DO ix=1,pmx
          WRITE(xy_unit) (/REAL(ix,4),REAL(x(ix),4),REAL(
     $      +2*x(ix)*jz_eq(ix)*bth_eq(ix)/xdim/
     $      (2*bz_eq(ix)-x(ix)*mu0*
     $       (jth_eq(ix)+jz_eq(ix)*bz_eq(ix)/bth_eq(ix)))**2,4),
     $      REAL(x(ix)*bz_eq(ix)/bth_eq(ix),4)/)
        ENDDO
        CALL close_bin(xy_unit,"suydam.bin")
c-----------------------------------------------------------------------
c       now run through nodes on the finite element mesh and use
c       lagr_1D to evaluate interpolation coefficients based on the
c       radius.
c-----------------------------------------------------------------------
        DO ibasis=1,SIZE(lxy%ix0)
          DO iy=lxy%iy0(ibasis),lxy%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,0_i4)
              r=SQRT((lxy%f(1)-xc)**2+(lxy%f(2)-yc)**2)/xdim
              IF (r>SQRT(TINY(r))) THEN
                theta=ATAN2(lxy%f(2)-yc,lxy%f(1)-xc)
              ELSE
                theta=0
              ENDIF
              IF (r<=SQRT(TINY(r))) THEN
                bey=bth_eq(0)
                bez=bz_eq(0)
                jy0=jth_eq(0)
                jz0=jz_eq(0)
                pr0=pr_eq(0)
              ELSE IF (r>=0.99999999_r8) THEN
                bey=bth_eq(pmx)
                bez=bz_eq(pmx)
                jy0=jth_eq(pmx)
                jz0=jz_eq(pmx)
                pr0=pr_eq(pmx)
              ELSE
                ixm=r*mx
                dr=(r*mx-ixm)
                ixm=ixm*poly_degree
                CALL lagr_1D(poly_degree,dr,al,dal,0_i4)
                bey=SUM(al*bth_eq(ixm:ixm+poly_degree))
                bez=SUM(al*bz_eq(ixm:ixm+poly_degree))
                jy0=SUM(al*jth_eq(ixm:ixm+poly_degree))
                jz0=SUM(al*jz_eq(ixm:ixm+poly_degree))
                pr0=SUM(al*pr_eq(ixm:ixm+poly_degree))
              ENDIF

              bvec=(/-SIN(theta)*bey,COS(theta)*bey,bez/)
              jvec=(/-SIN(theta)*jy0,COS(theta)*jy0,jz0/)
              CALL lagr_quad_basis_assign_loc(lbq,bvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(ljq,jvec,ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(lpq,(/pr0/),ibasis,ix,iy)
              IF (lamprof=='qspec') THEN		! Holmes E
                CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/jz_eq(0)/jz0/),ibasis,ix,iy)
              ELSE					! std profile
                CALL lagr_quad_basis_assign_loc(ldiff,
     $            (/MIN(dvac,(1+(SQRT(dvac)-1)*(r*xdim/xv)**dexp)**2)/),
     $            ibasis,ix,iy)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        CALL nim_stop('No suitable rect_cir initialization.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rect_shaped_init
c-----------------------------------------------------------------------
c     subprogram 3. rect_sh_slab.
c     generates equilibrium fields in linear geometry with rectangular
c     cross section for various parallel current profiles.
c-----------------------------------------------------------------------
      SUBROUTINE rect_sh_slab(x,by,bz,jy,jz,pr)

      REAL(r8), DIMENSION(0:), INTENT(INOUT) :: x,by,bz,jy,jz,pr

      INTEGER(i4) :: nr,i
      REAL(r8) :: lam,byp,bzp,bym,bzm,rc,dr,xdim,dpr,bmag2
      REAL(r8), DIMENSION(4) :: rjz,jth 
c-----------------------------------------------------------------------
c     do a Runge-Kutta integration to find the
c     equilibrium field from the specified current profile.
c-----------------------------------------------------------------------
      nr=SIZE(by)-1
      xdim=x(nr)-x(0)
      DO i=1,nr
        dr=x(i)-x(i-1)
        rc=0.5*(x(i)+x(i-1))
        bym=by(i-1)
        bzm=bz(i-1)
c-----------------------------------------------------------------------
c       Step 1:
c-----------------------------------------------------------------------
        lam=rect_sh_lam(x(i-1),bym,bzm)
        dpr=dp_fun(x(i-1))
        bmag2=bym**2+bzm**2
        jth(1)=lam*bym+bzm*dpr/bmag2
        rjz(1)=lam*bzm-bym*dpr/bmag2
        byp=bym+0.5*dr*rjz(1)
        bzp=bzm-0.5*dr*jth(1)
c-----------------------------------------------------------------------
c       Step 2:
c-----------------------------------------------------------------------
        lam=rect_sh_lam(rc,byp,bzp)
        dpr=dp_fun(rc)
        bmag2=byp**2+bzp**2
        jth(2)=lam*byp+bzp*dpr/bmag2
        rjz(2)=lam*bzp-byp*dpr/bmag2
        byp=bym+0.5*dr*rjz(2)
        bzp=bzm-0.5*dr*jth(2)
c-----------------------------------------------------------------------
c       Step 3:
c-----------------------------------------------------------------------
        lam=rect_sh_lam(rc,byp,bzp)
        bmag2=byp**2+bzp**2
        jth(3)=lam*byp+bzp*dpr/bmag2
        rjz(3)=lam*bzp-byp*dpr/bmag2
        byp=bym+dr*rjz(3)
        bzp=bzm-dr*jth(3)
c-----------------------------------------------------------------------
c       Step 4:
c-----------------------------------------------------------------------
        lam=rect_sh_lam(x(i),byp,bzp)
        dpr=dp_fun(x(i))
        bmag2=byp**2+bzp**2
        jth(4)=lam*byp+bzp*dpr/bmag2
        rjz(4)=lam*bzp-byp*dpr/bmag2
        by(i)=bym+dr/6.*(SUM(rjz(1:4:3))+2.*SUM(rjz(2:3)))
        bz(i)=bzm-dr/6.*(SUM(jth(1:4:3))+2.*SUM(jth(2:3)))
        jy(i)=lam*by(i)/mu0
        jz(i)=lam*bz(i)/mu0
      ENDDO

c-----------------------------------------------------------------------
c     pressure profile.
c-----------------------------------------------------------------------
      pr=beta*be0**2*(1+pres_2*((x-x(0))/xdim)**2
     $                 +pres_4*((x-x(0))/xdim)**4)/(2._r8*mu0)
c-----------------------------------------------------------------------
c     current density.
c-----------------------------------------------------------------------
      DO i=0,nr
        bmag2=by(i)**2+bz(i)**2
        dpr=dp_fun(x(i))
        lam=rect_sh_lam(x(i),by(i),bz(i))
        jy(i)=(lam*by(i)+bz(i)*dpr/bmag2)/mu0
        jz(i)=(lam*bz(i)-by(i)*dpr/bmag2)/mu0
      ENDDO
      RETURN

      CONTAINS
c-----------------------------------------------------------------------
c       compute mu0*dp/dr.
c-----------------------------------------------------------------------
        FUNCTION dp_fun(rad)

        REAL(r8), INTENT(IN) :: rad
        REAL(r8) :: dp_fun

        dp_fun=beta*be0**2/(2*xdim)
     $             *(2*pres_2*((rad-x(0))/xdim)
     $              +4*pres_4*((rad-x(0))/xdim)**3)

        END FUNCTION dp_fun

      END SUBROUTINE rect_sh_slab
c-----------------------------------------------------------------------
c     subprogram 4. rect_sh_lam
c     finds the parallel current density for a given position
c     based on a specified profile.  note that the returned value
c     is in m**(-1).
c-----------------------------------------------------------------------
      FUNCTION rect_sh_lam(r,bth,bz)

      REAL(r8), INTENT(IN) :: r,bth,bz
      REAL(r8) :: rect_sh_lam,xdim
c-----------------------------------------------------------------------
c     modified bessel function model:
c-----------------------------------------------------------------------
      xdim=(xmax-xmin)/2
      SELECT CASE(TRIM(lamprof))
      CASE('mbfm')
        IF ((r-xmin)/xdim<=rbreak) THEN
          rect_sh_lam=lam0/xdim
        ELSE
          rect_sh_lam=lam0*( (xdim-r)/(1-rbreak) )/xdim**2
        ENDIF
c-----------------------------------------------------------------------
c     alpha model:
c-----------------------------------------------------------------------
      CASE('alpha')
        rect_sh_lam=lam0*(1-(r/xdim)**alpha)/xdim
c-----------------------------------------------------------------------
c     paramagnetic pinch:
c-----------------------------------------------------------------------
      CASE('para')
        rect_sh_lam=lam0*be0*bz/(xdim*(bz**2+bth**2))
        IF (ds_use=='equil') rect_sh_lam=rect_sh_lam/
     $      (1+(SQRT(dvac)-1)*(r/xdim)**dexp)**2
c----------------------------------------------------------------------
c    sheet pinch (lam0=mu0*xdmim*J/B sets the parallel current
c      density on axis, and SIN(phib)*be0 sets By at the walls):
c-----------------------------------------------------------------------
      CASE('sheet')
        rect_sh_lam=(lam0/xdim)*
     $    ( 2/( EXP( lam0*r/(xdim*SIN(phib)+TINY(lam0)))
     $         +EXP(-lam0*r/(xdim*SIN(phib)+TINY(lam0))) ) )**2
      CASE DEFAULT
        CALL nim_stop('Unrecognized lambda profile.')
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION rect_sh_lam
c-----------------------------------------------------------------------
c     subprogram 5. block_init.
c     initializes the grid blocks with geometry and equilibrium fields.
c-----------------------------------------------------------------------
      SUBROUTINE block_init(nrbl)

      INTEGER(i4), INTENT(OUT) :: nrbl

      INTEGER(i4) :: ix,iy,ib,ixbl,iybl,mx1,my1,mx2,my2,mxb,myb,iq,imin,
     $               ibasis
      REAL(r8), DIMENSION(3) :: veq,xext,ext_prod
      REAL(r8) :: eta,bmag2,lam
      REAL(r8) :: ey,ez,x_grmx,gg,xdim,xcen,dx,dy
c-----------------------------------------------------------------------
c     create the equilibrium flow profile.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lvq,mx,my,3_i4,poly_degree)
      veq(1)=ve0*SIN(pi*thetav)*COS(pi*phiv)
      veq(2)=ve0*SIN(pi*thetav)*SIN(pi*phiv)
      veq(3)=ve0*COS(pi*thetav)
      IF (poly_degree>1) THEN
        x_grmx=MAX(MAXVAL(lxy%fs(1,:,:)),MAXVAL(lxy%fsv(1,:,:,:)))
      ELSE
        x_grmx=MAXVAL(lxy%fs(1,:,:))
      ENDIF
      SELECT CASE(eq_flow)
c-----------------------------------------------------------------------
c     no flow--set ve_eq to 0.
c-----------------------------------------------------------------------
      CASE('none')
       lvq=0
c-----------------------------------------------------------------------
c     uniform axial flow or rigid rotation for toroidal geometry.
c-----------------------------------------------------------------------
      CASE('uniform','gauss')
        IF (geom=='tor') THEN
          lvq%fs(1,:,:)=veq(1)
          lvq%fs(2,:,:)=veq(2)
          lvq%fs(3,:,:)=veq(3)*lxy%fs(1,:,:)/x_grmx
          IF (poly_degree>1) THEN
            lvq%fsh(1,:,:,:)=veq(1)
            lvq%fsh(2,:,:,:)=veq(2)
            lvq%fsh(3,:,:,:)=veq(3)*lxy%fsh(1,:,:,:)/x_grmx
            lvq%fsv(1,:,:,:)=veq(1)
            lvq%fsv(2,:,:,:)=veq(2)
            lvq%fsv(3,:,:,:)=veq(3)*lxy%fsv(1,:,:,:)/x_grmx
            lvq%fsi(1,:,:,:)=veq(1)
            lvq%fsi(2,:,:,:)=veq(2)
            lvq%fsi(3,:,:,:)=veq(3)*lxy%fsi(1,:,:,:)/x_grmx
          ENDIF
        ELSE
          lvq%fs(1,:,:)=veq(1)
          lvq%fs(2,:,:)=veq(2)
          lvq%fs(3,:,:)=veq(3)
          IF (poly_degree>1) THEN
            lvq%fsh(1,:,:,:)=veq(1)
            lvq%fsh(2,:,:,:)=veq(2)
            lvq%fsh(3,:,:,:)=veq(3)
            lvq%fsv(1,:,:,:)=veq(1)
            lvq%fsv(2,:,:,:)=veq(2)
            lvq%fsv(3,:,:,:)=veq(3)
            lvq%fsi(1,:,:,:)=veq(1)
            lvq%fsi(2,:,:,:)=veq(2)
            lvq%fsi(3,:,:,:)=veq(3)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c     gauss is the uniform case modulated by a gaussian profile in the
c     horizontal direction.
c-----------------------------------------------------------------------
        IF (eq_flow=='gauss') THEN
          xdim=0.5*(xmax-xmin)
          xcen=0.5*(xmax+xmin)
          DO ibasis=1,SIZE(lvq%ix0)
            DO iy=lvq%iy0(ibasis),lvq%my
              DO ix=lvq%ix0(ibasis),lvq%mx
                dx=ix-lvq%ix0(ibasis)+lvq%dx(ibasis)
                dy=iy-lvq%iy0(ibasis)+lvq%dy(ibasis)
                CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                CALL lagr_quad_eval(lvq,dx,dy,0_i4)
                gg=EXP(-(lxy%f(1)-xcen)**2/(eqflow_width**2*xdim**2))
                veq=gg*lvq%f
                CALL lagr_quad_basis_assign_loc(lvq,veq,ibasis,ix,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c     mhd pinch option for 'rect' grids based on an inferred E.
c-----------------------------------------------------------------------
      CASE('pinch')
        eta=mu0*elecd
        IF (lam0==0.AND.lamprof/='pitprs') THEN
          ey=0
          ez=0
        ELSE IF (geom=='lin'.AND.lamprof/='pitprs') THEN
          ey=ljq%fs(2,mx/2,0)*eta
          ez=ljq%fs(3,mx/2,0)*eta
        ELSE  !  geom='tor'
          ey=ljq%fs(2,0,0)*eta
          ez=0
        ENDIF
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
                IF (lxy%f(1)>0) lbq%f(3)=lbq%f(3)/lxy%f(1)
                ljq%f(3)=ljq%f(3)*lxy%f(1)
              ENDIF
              bmag2=SUM(lbq%f**2)
              veq(1)= ( (ey-eta*ljq%f(2))*lbq%f(3)
     $                 -(ez-eta*ljq%f(3))*lbq%f(2))/bmag2
              veq(2)= (     eta*ljq%f(1) *lbq%f(3)
     $                 +(ez-eta*ljq%f(3))*lbq%f(1))/bmag2
              veq(3)= (    -eta*ljq%f(1) *lbq%f(2)
     $                 -(ey-eta*ljq%f(2))*lbq%f(1))/bmag2
              CALL lagr_quad_basis_assign_loc(lvq,veq,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     ion diamagnetic flow as (1/(1+meomi)-pe_frac)*J_perp/(n_e*e)
c-----------------------------------------------------------------------
      CASE('diamagnetic')
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
                IF (lxy%f(1)>1.e-6_r8*xmax/(poly_degree*mx))
     $            lbq%f(3)=lbq%f(3)/lxy%f(1)
                ljq%f(3)=ljq%f(3)*lxy%f(1)
              ENDIF
              bmag2=SUM(lbq%f**2)
              lam=SUM(lbq%f*ljq%f)/bmag2
              veq=(1._r8/(1._r8+meomi)-pe_frac)*(ljq%f-lam*lbq%f)/
     $            (lnd%f(1)*elementary_q)
              CALL lagr_quad_basis_assign_loc(lvq,veq,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     equilibrium files are only for flux gridshapes.
c-----------------------------------------------------------------------
      CASE('eq_file')
      CASE DEFAULT
        CALL nim_stop
     $  ('Equilibrium flow profile '//TRIM(eq_flow)//' not recognized.')
      END SELECT
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      nrbl=nxbl*nybl
      ALLOCATE(rb(nrbl))
c-----------------------------------------------------------------------
c     rectangular grids contain no degenerate cells.
c-----------------------------------------------------------------------
      DO ib=1,nrbl
        rb(ib)%degenerate=.false.
      ENDDO
c-----------------------------------------------------------------------
c     start loop over blocks.
c-----------------------------------------------------------------------
      ib=0
      mx2=0
      DO ixbl=1,nxbl
         mx1=mx2
         mx2=mx*ixbl/nxbl
         mxb=mx2-mx1
         my2=0
         DO iybl=1,nybl
            my1=my2
            my2=my*iybl/nybl
            myb=my2-my1
            ib=ib+1
            rb(ib)%mx=mxb
            rb(ib)%my=myb
c-----------------------------------------------------------------------
c     fill coordinates.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%rz,mxb,myb,2_i4,poly_degree,
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
            CALL lagr_quad_alloc(rb(ib)%be_eq,mxb,myb,3_i4,poly_degree,
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
            CALL lagr_quad_alloc(rb(ib)%ja_eq,mxb,myb,3_i4,poly_degree,
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
            CALL lagr_quad_alloc(rb(ib)%ve_eq,mxb,myb,3_i4,poly_degree,
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
            CALL lagr_quad_alloc(rb(ib)%pres_eq,mxb,myb,1_i4,
     $        poly_degree,name='lpres',title=(/'lpres '/))
            rb(ib)%pres_eq%fs=lpq%fs(:,mx1:mx2,my1:my2)
            IF (poly_degree>1) THEN
              rb(ib)%pres_eq%fsh=lpq%fsh(:,:,1+mx1:mx2,my1:my2)
              rb(ib)%pres_eq%fsv=lpq%fsv(:,:,mx1:mx2,1+my1:my2)
              rb(ib)%pres_eq%fsi=lpq%fsi(:,:,1+mx1:mx2,1+my1:my2)
            ENDIF
c-----------------------------------------------------------------------
c     allocate and fill prese_eq.
c-----------------------------------------------------------------------
            CALL lagr_quad_alloc(rb(ib)%prese_eq,mxb,myb,1_i4,
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
            CALL lagr_quad_alloc(rb(ib)%nd_eq,mxb,myb,1_i4,poly_degree,
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
            CALL lagr_quad_alloc(rb(ib)%diff_shape,mxb,myb,1_i4,
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
      END SUBROUTINE block_init
c-----------------------------------------------------------------------
c     subprogram 6. seam_init.
c     initializes the grid seams for rectangular blocks.
c-----------------------------------------------------------------------
      SUBROUTINE seam_init

      INTEGER(i4) :: ixbl,iybl,ib,mxb,myb,iv,ix,iy,ip,nv,ivlim
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: np,jd
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: jb,jv
      INTEGER(i4), DIMENSION(2,-nybl:nxbl*nybl+nybl+1) :: kb
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ixv,iyv
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(seam(nxbl*nybl))
c-----------------------------------------------------------------------
c     create a temporary array containing block dimensions with
c     a rim of ghost blocks to prevent operand range errors later.
c-----------------------------------------------------------------------
      kb=0
      ib=0
      DO ixbl=1,nxbl
        DO iybl=1,nybl
          ib=ib+1
          kb(1,ib)=rb(ib)%mx
          kb(2,ib)=rb(ib)%my
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     start loop over interior blocks.
c-----------------------------------------------------------------------
      ib=0
      DO ixbl=1,nxbl
         DO iybl=1,nybl
            ib=ib+1
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
            mxb=rb(ib)%mx
            myb=rb(ib)%my
            nv=2*mxb+2*myb
            seam(ib)%nvert=nv
            ALLOCATE(seam(ib)%vertex(nv))
            ALLOCATE(np(nv))
            ALLOCATE(jb(nv,3))
            ALLOCATE(jv(nv,3))
            ALLOCATE(jd(nv))
            np=(/(1,ix=1,mxb-1),3,(1,iy=1,myb-1),3,
     $           (1,ix=1,mxb-1),3,(1,iy=1,myb-1),3/)
c-----------------------------------------------------------------------
c     fill block pointers.
c-----------------------------------------------------------------------
            jb=0
            jb(:,1)=(/(ib-1_i4,ix=1,mxb),(ib+nybl,iy=1,myb),
     $                (ib+1_i4,ix=1,mxb),(ib-nybl,iy=1,myb)/)
            IF (TRIM(periodicity)/='none') THEN
              IF (nybl==1) THEN
                jb(1:mxb,1)=ib
                jb(mxb+myb+1:2*mxb+myb,1)=ib
              ELSE IF (iybl==1) THEN
                jb(1:mxb,1)=ib+nybl-1
              ELSE IF (iybl==nybl) THEN
                jb(mxb+myb+1:2*mxb+myb,1)=ib-nybl+1
              ENDIF
              IF (TRIM(periodicity)=='both') THEN
                IF (nxbl==1) THEN
                  jb(mxb+1:mxb+myb,1)=ib
                  jb(2*mxb+myb+1:2*(mxb+myb),1)=ib
                ELSE IF (ixbl==1) THEN
                  jb(2*mxb+myb+1:2*(mxb+myb),1)=ib+nybl*(nxbl-1)
                ELSE IF (ixbl==nxbl) THEN
                  jb(mxb+1:mxb+myb,1)=ib-nybl*(nxbl-1)
                ENDIF
              ENDIF
            ENDIF
            jb(mxb,2)=jb(mxb,1)+nybl
            IF (jb(mxb,2)>nbl) jb(mxb,2)=jb(mxb,2)-nbl
            jb(mxb,3)=ib+nybl
            IF (jb(mxb,3)>nbl) jb(mxb,3)=jb(mxb,3)-nbl
            jb(mxb+myb,2)=jb(mxb+myb+1,1)+nybl
            IF (jb(mxb+myb,2)>nbl) jb(mxb+myb,2)=jb(mxb+myb,2)-nbl
            jb(mxb+myb,3)=jb(mxb+myb+1,1)
            jb(2*mxb+myb,2)=jb(2*mxb+myb,1)-nybl
            IF (jb(2*mxb+myb,2)<1) jb(2*mxb+myb,2)=jb(2*mxb+myb,2)+nbl
            jb(2*mxb+myb,3)=ib-nybl
            IF (jb(2*mxb+myb,3)<1) jb(2*mxb+myb,3)=jb(2*mxb+myb,3)+nbl
            jb(2*mxb+2*myb,2)=jb(1,1)-nybl
            IF (jb(2*mxb+2*myb,2)<1) jb(2*mxb+2*myb,2)=
     $                               jb(2*mxb+2*myb,2)+nbl
            jb(2*mxb+2*myb,3)=jb(1,1)
c-----------------------------------------------------------------------
c     fill vertex pointers.
c-----------------------------------------------------------------------
            jv=0
            jv(:,1)=(/(2*mxb+kb(2,jb(iv,1))-iv,iv=1,mxb),
     $           (2*kb(1,jb(iv,1))+2*myb+mxb-iv,iv=mxb+1,mxb+myb),
     $           (2*mxb+myb-iv,iv=1+mxb+myb,2*mxb+myb),
     $           (kb(1,jb(iv,1))+2*myb+2*mxb-iv,
     $                                     iv=2*mxb+myb+1,2*(mxb+myb))/)
            jv(mxb,2)=2*kb(1,jb(mxb,2))+kb(2,jb(mxb,2))
            jv(mxb,3)=2*kb(1,jb(mxb,3))+2*kb(2,jb(mxb,3))
            jv(mxb+myb,2)=2*kb(1,jb(mxb+myb,2))+2*kb(2,jb(mxb+myb,2))
            jv(mxb+myb,3)=kb(1,jb(mxb+myb,3))
            jv(2*mxb+myb,1)=2*(kb(1,jb(2*mxb+myb,1))
     $                        +kb(2,jb(2*mxb+myb,1)))
            jv(2*mxb+myb,2)=kb(1,jb(2*mxb+myb,2))
            jv(2*mxb+myb,3)=kb(1,jb(2*mxb+myb,3))+kb(2,jb(2*mxb+myb,3))
            jv(2*(mxb+myb),2)=kb(1,jb(2*(mxb+myb),2))
     $           +kb(2,jb(2*(mxb+myb),2))
            jv(2*(mxb+myb),3)=2*kb(1,jb(2*(mxb+myb),3))
     $           +kb(2,jb(2*(mxb+myb),3))
c-----------------------------------------------------------------------
c     initialize pointer delete indicator for external boundaries.
c-----------------------------------------------------------------------
            jd=0
c-----------------------------------------------------------------------
c     trim left pointers.
c-----------------------------------------------------------------------
            IF(ixbl==1.AND.TRIM(periodicity)/='both')THEN
               iv=edge_match(seam0%vertex,seam0%nvert,
     $                       ib,2_i4*mxb+myb+1_i4)
               jb(2*mxb+myb+1:2*(mxb+myb),1)=0
               jv(2*mxb+myb+1:2*(mxb+myb),1)=
     $                      (/(iv+iy,iy=0_i4,myb-1_i4)/)
               jd(2*(mxb+myb))=1
               IF (iybl==1.AND.TRIM(periodicity)/='y-dir') THEN
                 jd(2*(mxb+myb))=2
               ENDIF
               IF (iybl==nybl) THEN
                 IF (TRIM(periodicity)=='y-dir') THEN
                   jb(2*mxb+myb,3)=0
                   jv(2*mxb+myb,3)=edge_match(seam0%vertex,
     $                           seam0%nvert,ib,2_i4*mxb+myb)
                   jd(2*mxb+myb)=1
                 ENDIF
               ELSE
                 jb(2*mxb+myb,3)=0
                 jv(2*mxb+myb,3)=iv-1
                 jd(2*mxb+myb)=1
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c     trim right pointers.
c-----------------------------------------------------------------------
            IF(ixbl==nxbl.AND.TRIM(periodicity)/='both')THEN
               iv=edge_match(seam0%vertex,seam0%nvert,
     $                       ib,mxb+1_i4)
               jb(mxb+1:mxb+myb,1)=0
               jv(mxb+1:mxb+myb,1)=(/(iv+iy,iy=0_i4,myb-1_i4)/)
               jd(mxb+myb)=1
               IF (iybl==nybl.AND.TRIM(periodicity)/='y-dir') THEN
                 jd(mxb+myb)=2
               ENDIF
               IF (iybl==1) THEN
                 IF (TRIM(periodicity)=='y-dir') THEN
                   jb(mxb,3)=0
                   jv(mxb,3)=edge_match(seam0%vertex,
     $                      seam0%nvert,ib,mxb)
                   jd(mxb)=1
                 ENDIF
               ELSE
                 jb(mxb,3)=0
                 jv(mxb,3)=iv-1
                 jd(mxb)=1
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c     trim bottom pointers.
c-----------------------------------------------------------------------
            IF(iybl==1.AND.TRIM(periodicity)=='none') THEN
               iv=edge_match(seam0%vertex,seam0%nvert,ib,1_i4)
               jb(1:mxb,1)=0
               jv(1:mxb,1)=(/(iv+ix,ix=0_i4,mxb-1_i4)/)
               jd(mxb)=1
               IF(ixbl==nxbl)THEN
                 jd(mxb)=2
               ENDIF
               IF(ixbl/=1)THEN
                 jb(2*(mxb+myb),3)=0
                 jv(2*(mxb+myb),3)=iv-1
                 jd(2*(mxb+myb))=1
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c     trim top pointers.
c-----------------------------------------------------------------------
            IF(iybl==nybl.AND.TRIM(periodicity)=='none') THEN
               iv=edge_match(seam0%vertex,seam0%nvert,
     $                       ib,mxb+myb+1_i4)
               jb(1+mxb+myb:2*mxb+myb,1)=0
               jv(1+mxb+myb:2*mxb+myb,1)=(/(iv+ix,ix=0_i4,mxb-1_i4)/)
               jd(2*mxb+myb)=1
               IF(ixbl==1)THEN
                 jd(2*mxb+myb)=2
               ENDIF
               IF(ixbl/=nxbl)THEN
                 jb(mxb+myb,3)=0
                 jv(mxb+myb,3)=iv-1
                 jd(mxb+myb)=1
               ENDIF
            ENDIF
c-----------------------------------------------------------------------
c     fill pointers.
c-----------------------------------------------------------------------
            DO iv=1,nv
               IF (np(iv)-jd(iv)<=0) THEN
                 CALL nim_stop('Error assigning seam pointers.')
               ENDIF
               CALL ptr_alloc(seam(ib)%vertex(iv),2_i4,np(iv)-jd(iv))
               seam(ib)%vertex(iv)%ptr(1,1)=jb(iv,1)
               seam(ib)%vertex(iv)%ptr(2,1)=jv(iv,1)
               IF (jd(iv)==0) THEN
                 DO ip=2,np(iv)
                   seam(ib)%vertex(iv)%ptr(1,ip)=jb(iv,ip)
                   seam(ib)%vertex(iv)%ptr(2,ip)=jv(iv,ip)
                 ENDDO
               ELSE IF (jd(iv)==1) THEN
                 seam(ib)%vertex(iv)%ptr(1,2)=jb(iv,3)
                 seam(ib)%vertex(iv)%ptr(2,2)=jv(iv,3)
               ENDIF
            ENDDO
c-----------------------------------------------------------------------
c     finish loop over blocks.
c-----------------------------------------------------------------------
            DEALLOCATE(np)
            DEALLOCATE(jb)
            DEALLOCATE(jv)
            DEALLOCATE(jd)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     write seam information if detflag=.true.
c-----------------------------------------------------------------------
      IF (detflag) THEN
        WRITE(out_unit,*) 'ibl ', 0, ' nvert ',seam0%nvert
        DO iv=1,seam0%nvert
         ip=SIZE(seam0%vertex(iv)%ptr,2)
         WRITE(out_unit,*) iv, ' block ptrs ',
     $                     seam0%vertex(iv)%ptr(1,1:ip)
         WRITE(out_unit,*) iv, ' vertx ptrs ',
     $                     seam0%vertex(iv)%ptr(2,1:ip)
        ENDDO
        DO ib=1,nbl
          WRITE(out_unit,*) 'ibl ', ib, ' nvert ',seam(ib)%nvert
          DO iv=1,seam(ib)%nvert
           ip=SIZE(seam(ib)%vertex(iv)%ptr,2)
           WRITE(out_unit,*) iv, ' block ptrs ',
     $                       seam(ib)%vertex(iv)%ptr(1,1:ip)
           WRITE(out_unit,*) iv, ' vertx ptrs ',
     $                       seam(ib)%vertex(iv)%ptr(2,1:ip)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     find and save block-internal vertex coordinates.
c-----------------------------------------------------------------------
      DO ib=1,nbl
        mxb=rb(ib)%mx
        myb=rb(ib)%my
        ALLOCATE(ixv(2*(mxb+myb)))
        ALLOCATE(iyv(2*(mxb+myb)))
        ixv=(/(ix,ix=1,mxb),(mxb,iy=1,myb),
     $       (ix,ix=mxb-1,0,-1),(0_i4,iy=myb-1,0,-1)/)
        iyv=(/(0_i4,ix=1,mxb),(iy,iy=1,myb),
     $       (myb,ix=mxb-1,0,-1),(iy,iy=myb-1,0,-1)/)
        DO iv=1,seam(ib)%nvert
          seam(ib)%vertex(iv)%intxy(1)=ixv(iv)
          seam(ib)%vertex(iv)%intxy(2)=iyv(iv)
        ENDDO
        DEALLOCATE(ixv)
        DEALLOCATE(iyv)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE seam_init
c-----------------------------------------------------------------------
c     subprogram 7. seam0_init.
c     create the seam representing the external boundary.
c-----------------------------------------------------------------------
      SUBROUTINE seam0_init

      INTEGER(i4) :: ixbl,iybl,ib,iv,ix,iy,ip,nv,ivlim
c-----------------------------------------------------------------------
c     evaluate pointers for the external boundary.  this is a separate 
c     structure which is not used after nimrod is initialized.
c-----------------------------------------------------------------------
      ALLOCATE(seam0)
      IF (gridshape(1:4)=='rect') THEN
        IF (TRIM(periodicity)=='y-dir') THEN
          nv=2*my
        ELSE IF (TRIM(periodicity)=='both') THEN
          nv=0
        ELSE
          nv=2*mx+2*my
        ENDIF
      ELSE IF (gridshape=='circ'.AND.xmin/=0) THEN
        nv=2*my
      ELSE
        nv=my
      ENDIF
      seam0%nvert=nv

      IF (TRIM(periodicity)/='both'.OR.gridshape(1:4)/='rect') THEN
c-----------------------------------------------------------------------
c>>>> start of periodicity/='both' IF-THEN block.
c-----------------------------------------------------------------------
      ALLOCATE(seam0%vertex(nv))
      ALLOCATE(seam0%excorner(nv))
      seam0%excorner=.false.
      ib=1
      iv=1
c-----------------------------------------------------------------------
c     traverse bottom of blocks.
c-----------------------------------------------------------------------
      yper1: IF (TRIM(periodicity)/='y-dir'.AND.
     $           gridshape(1:4)=='rect') THEN
        xblock1: DO ixbl=1,nxbl-1
          xvert1: DO ix=1,rb(ib)%mx-1
            CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
            seam0%vertex(iv)%ptr(1,1)=ib
            seam0%vertex(iv)%ptr(2,1)=ix
            iv=iv+1
          ENDDO xvert1
          CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx
          ib=ib+nybl
          seam0%vertex(iv)%ptr(1,2)=ib
          seam0%vertex(iv)%ptr(2,2)=2*(rb(ib)%mx+rb(ib)%my)
          iv=iv+1
        ENDDO xblock1
c-----------------------------------------------------------------------
c     corner block.
c-----------------------------------------------------------------------
        xvert2: DO ix=1,rb(ib)%mx
          CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=ix
          iv=iv+1
        ENDDO xvert2
c-----------------------------------------------------------------------
c     no bottom if periodic in y-direction.
c-----------------------------------------------------------------------
      ELSE yper1
        ib=ib+nybl*(nxbl-1)
      ENDIF yper1
c-----------------------------------------------------------------------
c     up right side of blocks.
c-----------------------------------------------------------------------
      IF (TRIM(periodicity)/='y-dir'.AND.gridshape=='rect') 
     $       seam0%excorner(iv-1)=.true.
      yblock1: DO iybl=1,nybl-1
        yvert1: DO iy=1,rb(ib)%my-1
          CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+iy
          iv=iv+1
        ENDDO yvert1
        CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
        seam0%vertex(iv)%ptr(1,1)=ib
        seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+rb(ib)%my
        ib=ib+1
        seam0%vertex(iv)%ptr(1,2)=ib
        seam0%vertex(iv)%ptr(2,2)=rb(ib)%mx
        iv=iv+1
      ENDDO yblock1
c-----------------------------------------------------------------------
c     corner block.
c-----------------------------------------------------------------------
      IF (TRIM(periodicity)=='y-dir'.OR.gridshape(1:4)/='rect') THEN
        ivlim=rb(ib)%my-1
      ELSE
        ivlim=rb(ib)%my
      ENDIF
      yvert2: DO iy=1,ivlim
        CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
        seam0%vertex(iv)%ptr(1,1)=ib
        seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+iy
        iv=iv+1
      ENDDO yvert2
      IF (TRIM(periodicity)=='y-dir'.OR.gridshape(1:4)/='rect') THEN
        CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
        seam0%vertex(iv)%ptr(1,1)=ib
        seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+rb(ib)%my
        seam0%vertex(iv)%ptr(1,2)=ib+1-nybl
        seam0%vertex(iv)%ptr(2,2)=rb(ib+1-nybl)%mx
        iv=iv+1
      ELSE IF (gridshape=='rect') THEN
        seam0%excorner(iv-1)=.true.
      ENDIF
c-----------------------------------------------------------------------
c     traverse top of blocks.
c-----------------------------------------------------------------------
      yper2: IF (TRIM(periodicity)/='y-dir'.AND.
     $           gridshape(1:4)=='rect') THEN
        xblock2: DO ixbl=nxbl,2,-1
          xvert3: DO ix=1,rb(ib)%mx-1
            CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
            seam0%vertex(iv)%ptr(1,1)=ib
            seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+rb(ib)%my+ix
            iv=iv+1
          ENDDO xvert3
          CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=2*rb(ib)%mx+rb(ib)%my
          ib=ib-nybl
          seam0%vertex(iv)%ptr(1,2)=ib
          seam0%vertex(iv)%ptr(2,2)=rb(ib)%mx+rb(ib)%my
          iv=iv+1
        ENDDO xblock2
c-----------------------------------------------------------------------
c     corner block.
c-----------------------------------------------------------------------
        xvert4: DO ix=1,rb(ib)%mx
          CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=rb(ib)%mx+rb(ib)%my+ix
          iv=iv+1
        ENDDO xvert4
c-----------------------------------------------------------------------
c     no top if periodic in y-direction.
c-----------------------------------------------------------------------
      ELSE yper2
        ib=ib-nybl*(nxbl-1)
      ENDIF yper2
c-----------------------------------------------------------------------
c     down left side of blocks.
c-----------------------------------------------------------------------
      IF (TRIM(periodicity)/='y-dir'.AND.gridshape=='rect') 
     $      seam0%excorner(iv-1)=.true.
      left: IF (gridshape(1:4)=='rect'.OR.
     $          gridshape=='circ'.AND.xmin/=0) THEN
        yblock2: DO iybl=nybl,2,-1
          yvert3: DO iy=1,rb(ib)%my-1
            CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
            seam0%vertex(iv)%ptr(1,1)=ib
            seam0%vertex(iv)%ptr(2,1)=2*rb(ib)%mx+rb(ib)%my+iy
            iv=iv+1
          ENDDO yvert3
          CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=2*(rb(ib)%mx+rb(ib)%my)
          ib=ib-1
          seam0%vertex(iv)%ptr(1,2)=ib
          seam0%vertex(iv)%ptr(2,2)=2*rb(ib)%mx+rb(ib)%my
          iv=iv+1
        ENDDO yblock2
c-----------------------------------------------------------------------
c       corner block.
c-----------------------------------------------------------------------
        IF (TRIM(periodicity)=='y-dir'.AND.gridshape=='rect'
     $      .OR.gridshape=='circ'.AND.xmin/=0) THEN
          ivlim=rb(ib)%my-1
        ELSE
          ivlim=rb(ib)%my
        ENDIF
        yvert4: DO iy=1,ivlim
          CALL ptr_alloc(seam0%vertex(iv),2_i4,1_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=2*rb(ib)%mx+rb(ib)%my+iy
          iv=iv+1
        ENDDO yvert4
        IF (TRIM(periodicity)=='y-dir'.AND.gridshape=='rect'
     $      .OR.gridshape=='circ'.AND.xmin/=0) THEN
          CALL ptr_alloc(seam0%vertex(iv),2_i4,2_i4)
          seam0%vertex(iv)%ptr(1,1)=ib
          seam0%vertex(iv)%ptr(2,1)=2*(rb(ib)%mx+rb(ib)%my)
          seam0%vertex(iv)%ptr(1,2)=ib+nybl-1
          seam0%vertex(iv)%ptr(2,2)=2*rb(ib+nybl-1)%mx
     $                               +rb(ib+nybl-1)%my
          iv=iv+1
        ELSE IF (gridshape=='rect') THEN
          seam0%excorner(iv-1)=.true.
        ENDIF
      ENDIF left
      IF (iv.ne.nv+1) THEN
        CALL nim_stop('Error in seam0 initialization.')
      ENDIF
c-----------------------------------------------------------------------
c>>>> end of periodicity/='both' IF-THEN block.
c-----------------------------------------------------------------------
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE seam0_init
c-----------------------------------------------------------------------
c     subprogram 8. var0_alloc.
c     allocates space for the fundamental dependent varibles and 
c     matrices, i.e., those written to restart files.
c-----------------------------------------------------------------------
      SUBROUTINE var0_alloc(nmodes)

      INTEGER(i4), INTENT(IN) :: nmodes

      INTEGER(i4) :: ibl,mxb,myb,mv
c-----------------------------------------------------------------------
c     allocate storage for the real and imaginary parts of the
c     dependent varibles.  rblocks first.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        mxb=rb(ibl)%mx
        myb=rb(ibl)%my
        IF (geom=='tor') THEN
          CALL lagr_quad_alloc(rb(ibl)%be,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='be',title=(/' be_r ',' be_z ',' be_p '/))
          CALL lagr_quad_alloc(rb(ibl)%ja,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='ja',title=(/' ja_r ',' ja_z ',' ja_p '/))
          CALL lagr_quad_alloc(rb(ibl)%ve,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='ve',title=(/' ve_r ',' ve_z ',' ve_p '/))
        ELSE
          CALL lagr_quad_alloc(rb(ibl)%be,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='be',title=(/' be_x ',' be_y ',' be_z '/))
          CALL lagr_quad_alloc(rb(ibl)%ja,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='ja',title=(/' ja_x ',' ja_y ',' ja_z '/))
          CALL lagr_quad_alloc(rb(ibl)%ve,mxb,myb,3_i4,nmodes,
     $      poly_degree,name='ve',title=(/' ve_x ',' ve_y ',' ve_z '/))
        ENDIF
        rb(ibl)%ja=0
        CALL lagr_quad_alloc(rb(ibl)%pres,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='pr',title=(/' pres '/))
        CALL lagr_quad_alloc(rb(ibl)%prese,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='pe',title=(/' prese'/))
        CALL lagr_quad_alloc(rb(ibl)%tele,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='te',title=(/' tele '/))
        CALL lagr_quad_alloc(rb(ibl)%tion,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='ti',title=(/' tion '/))
        CALL lagr_quad_alloc(rb(ibl)%nd,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='nd',title=(/' nden '/))
        rb(ibl)%nd=0
        CALL lagr_quad_alloc(rb(ibl)%conc,mxb,myb,1_i4,nmodes,
     $                       poly_degree,name='co',title=(/' conc '/))
        rb(ibl)%conc=0
      ENDDO
c-----------------------------------------------------------------------
c     tblocks
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        mv=tb(ibl)%mvert
        IF (geom=='tor') THEN
          CALL tri_linear_alloc(tb(ibl)%be,mv,3_i4,nmodes,name='be',
     $                          title=(/' be_r ',' be_z ',' be_p '/))
          CALL tri_linear_alloc(tb(ibl)%ja,mv,3_i4,nmodes,name='ja',
     $                          title=(/' ja_r ',' ja_z ',' ja_p '/))
          CALL tri_linear_alloc(tb(ibl)%ve,mv,3_i4,nmodes,name='ve',
     $                          title=(/' ve_r ',' ve_z ',' ve_p '/))
        ELSE
          CALL tri_linear_alloc(tb(ibl)%be,mv,3_i4,nmodes,name='be',
     $                          title=(/' be_x ',' be_y ',' be_z '/))
          CALL tri_linear_alloc(tb(ibl)%ja,mv,3_i4,nmodes,name='ja',
     $                          title=(/' ja_x ',' ja_y ',' ja_z '/))
          CALL tri_linear_alloc(tb(ibl)%ve,mv,3_i4,nmodes,name='ve',
     $                          title=(/' ve_x ',' ve_y ',' ve_z '/))
        ENDIF
        tb(ibl)%ja=0
        CALL tri_linear_alloc(tb(ibl)%pres,mv,1_i4,nmodes,name='pr',
     $                        title=(/' pres '/))
        CALL tri_linear_alloc(tb(ibl)%prese,mv,1_i4,nmodes,name='pe',
     $                        title=(/' prese'/))
        CALL tri_linear_alloc(tb(ibl)%tele,mv,1_i4,nmodes,name='te',
     $                        title=(/' tele '/))
        CALL tri_linear_alloc(tb(ibl)%tion,mv,1_i4,nmodes,name='ti',
     $                        title=(/' tion '/))
        CALL tri_linear_alloc(tb(ibl)%nd,mv,1_i4,nmodes,name='nd',
     $                        title=(/' nden '/))
        tb(ibl)%nd=0
        CALL tri_linear_alloc(tb(ibl)%conc,mv,1_i4,nmodes,name='co',
     $                        title=(/' conc '/))
        tb(ibl)%conc=0
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE var0_alloc
c-----------------------------------------------------------------------
c     subprogram 9. ptr_alloc.
c     allocate ptr under a seam vertex.  this seems to be necessary due
c     to a bug in the ibm compiler.
c-----------------------------------------------------------------------
      SUBROUTINE ptr_alloc(vert,idim1,idim2)

      TYPE(vertex_type), INTENT(INOUT) :: vert
      INTEGER(i4), INTENT(IN) :: idim1,idim2
c-----------------------------------------------------------------------
c     allocate ptr.
c-----------------------------------------------------------------------
      ALLOCATE(vert%ptr(idim1,idim2))
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ptr_alloc
c-----------------------------------------------------------------------
c     subprogram 10. edge_match
c     matches the external boundary pointer with internal block and
c     vertex indices.  returns the external seam vertex.
c-----------------------------------------------------------------------
      FUNCTION edge_match(vert0,nvert,intbl,intvrt) RESULT(iv0)

      INTEGER(i4), INTENT(IN) :: nvert,intbl,intvrt
      TYPE(vertex_type), DIMENSION(nvert), INTENT(IN) :: vert0

      INTEGER(i4) :: iv0,ip
      LOGICAL :: match
      CHARACTER(64) :: message
c-----------------------------------------------------------------------
c     loop over external vertices
c-----------------------------------------------------------------------
      iv0=1
      match=.false.
      ext: DO 
        ptr: DO ip=1,SIZE(vert0(iv0)%ptr,2)
          IF (vert0(iv0)%ptr(1,ip)==intbl) THEN
            IF (vert0(iv0)%ptr(2,ip)==intvrt) THEN
              match=.true.
            ENDIF
          ENDIF
        ENDDO ptr
        IF (match) EXIT ext
        iv0=iv0+1
        IF (iv0.gt.nvert) THEN
          WRITE(message,'(2a)')'Error matching external seam',
     $                         ' in edge_match.'
          CALL nim_stop(message)
        ENDIF
      ENDDO ext
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION edge_match
c-----------------------------------------------------------------------
c     subprogram 11. seam0_poly_init
c     initializes seam0 for the wall boundary
c     when that boundary is one or more tblocks.
c-----------------------------------------------------------------------
      SUBROUTINE seam0_poly_init(nblstart)

      INTEGER(i4), INTENT(IN) :: nblstart

      INTEGER(i4) :: ibl,iv,jb,jv,ip,kp,jp,jvold,tv
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: ivm,ivp
      REAL(r8) :: r1,z1,r2,z2,r3,z3,theta,andotc
      REAL(r8) :: theta0=1.91986                ! maximum angle for sharp corner
c
c                                               Determine the number of 
c						seam0 vertici by looping
c						through the rim seams.
      jv=0_i4
      DO ibl=nblstart+1,nbl_rim+nblstart
        DO iv=1,seam(ibl)%nvert
          DO ip=1,SIZE(seam(ibl)%vertex(iv)%ptr,2)
          IF(seam(ibl)%vertex(iv)%ptr(1,ip).EQ.0_i4)THEN
            jv=jv+1_i4
          ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF(nbl_rim.NE.1)jv=jv-nbl_rim
      seam0%nvert=jv
      ALLOCATE(seam0%vertex(seam0%nvert))
      ALLOCATE(seam0%excorner(seam0%nvert))
      DO jv=1,seam0%nvert
        ALLOCATE(seam0%vertex(jv)%ptr(2,1))
      ENDDO
c
c						Loop through all rim tblock
c						seam vertici.
      DO ibl=nblstart+1,nblstart+nbl_rim
        DO iv=1,seam(ibl)%nvert
          jp=SIZE(seam(ibl)%vertex(iv)%ptr,2)
          DO ip=1,jp
c						Search for points which
c						point to seam0
            IF(seam(ibl)%vertex(iv)%ptr(1,ip).EQ.0_i4)THEN
              jvold=seam(ibl)%vertex(iv)%ptr(2,ip)
              IF(jp.EQ.1)THEN
c							noncorner seam0
                 seam0%vertex(jvold)%ptr(1,1)=ibl
                 seam0%vertex(jvold)%ptr(2,1)=iv
              ELSE
c							corner seam0
                DEALLOCATE(seam0%vertex(jvold)%ptr)
                ALLOCATE(seam0%vertex(jvold)%ptr(2,jp))
                DO kp=1,jp
                  seam0%vertex(jvold)%ptr(1,kp)=
     &                seam(ibl)%vertex(iv)%ptr(1,kp)
                  seam0%vertex(jvold)%ptr(2,kp)=
     &                seam(ibl)%vertex(iv)%ptr(2,kp)
                ENDDO
                seam0%vertex(jvold)%ptr(1,ip)=ibl
                seam0%vertex(jvold)%ptr(2,ip)=iv
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
c						Check for sharp angles on
c						the boundary.
c						A : vector 1 to 2 
c						B : vector 2 to 3 
c						C : vector perp bisector
c						 between 1 and 3 from 2
      seam0%excorner=.false.
      ALLOCATE(ivm(seam0%nvert))
      ALLOCATE(ivp(seam0%nvert))
      ivm=(/seam0%nvert,(iv,iv=1,seam0%nvert-1)/)
      ivp=(/(iv,iv=2,seam0%nvert),1_i4/)
c
c						Loop over seam0
c						points.

      DO iv=1,seam0%nvert
        jb=seam0%vertex(ivm(iv))%ptr(1,1)
        jv=seam0%vertex(ivm(iv))%ptr(2,1)
        tv=seam(jb)%vertex(jv)%intxy(1)
        r1=tb(jb)%tgeom%xs(tv)
        z1=tb(jb)%tgeom%ys(tv)

        jb=seam0%vertex(iv)%ptr(1,1)
        jv=seam0%vertex(iv)%ptr(2,1)
        tv=seam(jb)%vertex(jv)%intxy(1)
        r2=tb(jb)%tgeom%xs(tv)
        z2=tb(jb)%tgeom%ys(tv)

        jb=seam0%vertex(ivp(iv))%ptr(1,1)
        jv=seam0%vertex(ivp(iv))%ptr(2,1)
        tv=seam(jb)%vertex(jv)%intxy(1)
        r3=tb(jb)%tgeom%xs(tv)
        z3=tb(jb)%tgeom%ys(tv)
c
c						Dot product of C with normal
c						vector of A.
        andotc=-(z2-z1)*(0.5_r8*(r3+r1)-r2)
     &         +(r2-r1)*(0.5_r8*(z3+z1)-z2)
 
        theta=((r1-r2)*(r3-r2)+(z1-z2)*(z3-z2))
     &        /(((r1-r2)**2+(z1-z2)**2)
     &         *((r3-r2)**2+(z3-z2)**2))**0.5_r8
        IF(theta.GT.1._r8)THEN
          theta=0._r8
        ELSEIF(theta.LT.-1._r8)THEN
          theta=pi
        ELSE
          theta=ACOS(theta)
        ENDIF
        IF(andotc.LT.0._r8)THEN
          theta=twopi-theta
        ENDIF
        IF(theta.LT.theta0)seam0%excorner(iv)=.true.

      ENDDO
c
c						Diagnose seam0
      IF(detflag)THEN
        WRITE(out_unit,*)"seam0, mseam =",seam0%nvert
        WRITE(out_unit,fmt='(6x,"iv",4x,"ip",3x,"outb",2x,"outv")')
        DO iv=1,SIZE(seam0%vertex)
          DO ip=1,SIZE(seam0%vertex(iv)%ptr,2)
            WRITE(out_unit,fmt='(x,4i6)')iv,ip,
     &                 seam0%vertex(iv)%ptr(1,ip),
     &                 seam0%vertex(iv)%ptr(2,ip)
          ENDDO
        ENDDO
        WRITE(out_unit,*)"seam0 vertici with angles less than ",theta0
        DO iv=1,seam0%nvert
          WRITE(out_unit,*)iv,seam0%excorner(iv)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE seam0_poly_init
c-----------------------------------------------------------------------
c     subprogram 12. rect_quad_init.
c     generates a logically rectangular mesh with arbitrarily located
c     corners.  one side may also be bowed to match a circular arc.
c-----------------------------------------------------------------------
      SUBROUTINE rect_quad_init

      INTEGER(i4) :: ix,iy,jx,jy,iint
      REAL(r8) :: xdim,ydim,xc,yc,r,xcen,ycen,lone,ltwo,lbot,
     $            cross,sinth,xm,ym,xcp,ycp
      REAL(r8), DIMENSION(2) :: lonev,ltwov
      REAL(r8) :: agridx,agridy,bgridx,bgridy
      REAL(r8), DIMENSION(2,0:poly_degree*mx,0:poly_degree*my) :: xy
      REAL(r8), DIMENSION(:), ALLOCATABLE :: x,bth_eq,bz_eq,pr_eq,
     $          jth_eq,jz_eq
      REAL(r8), DIMENSION(2) :: xxyy,dir,dip
      REAL(r8), DIMENSION(3) :: bvec,jvec
      REAL(r8), DIMENSION(0:poly_degree) :: al,dal,x_node

      REAL(r8) :: bex,bey,bez,jy0,jz0,pr0,dx,dy
      REAL(r8), PARAMETER :: smallnum=1.e-5_r8
      INTEGER(i4) :: ibasis,ixm,iym
      CHARACTER(64) :: message

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
      IF (xmax<=xmin) CALL nim_stop("Rect_init: xmax must be > xmin.")
      IF (ymax<=ymin) CALL nim_stop("Rect_init: ymax must be > ymin.")
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lxy,mx,my,2_i4,poly_degree)
c-----------------------------------------------------------------------
c     find the parameters for linearly varying mesh spacing.
c-----------------------------------------------------------------------
      IF (firstx==0) THEN 
        agridx=0._r8
        bgridx=1._r8/mx
      ELSE IF (firstx>=2._r8/mx.OR.firstx<0) THEN
        WRITE(message,*) 'Inappropriate value for firstx:  ',firstx
        CALL nim_stop(message)
      ELSE
        agridx=(1._r8/mx-firstx)/(0.5_r8*(mx-1_i4))
        bgridx=firstx-agridx
      ENDIF

      IF (firsty==0) THEN 
        agridy=0._r8
        bgridy=1._r8/my
      ELSE IF (firsty>=2._r8/my.OR.firsty<0) THEN
        WRITE(message,*) 'Inappropriate value for firsty:  ',firsty
        CALL nim_stop(message)
      ELSE
        agridy=(1._r8/my-firsty)/(0.5_r8*(my-1_i4))
        bgridy=firsty-agridy
      ENDIF
c-----------------------------------------------------------------------
c     if one side of the region is curved, find the center of curvature
c     and place mesh points using trigonometry.
c-----------------------------------------------------------------------
      IF (quadarc) THEN
        lonev=(/quadxx(4)-quadxx(1),quadyy(4)-quadyy(1)/)
        ltwov=(/quadxx(3)-quadxx(2),quadyy(3)-quadyy(2)/)
        lone=SQRT(SUM(lonev**2))
        ltwo=SQRT(SUM(ltwov**2))
        lonev=lonev/lone
        ltwov=ltwov/ltwo
        cross=(ltwov(1)*lonev(2)-lonev(1)*ltwov(2))
        IF (ABS(cross)<=smallnum)
     $    CALL nim_stop('Quadarc option needs non-parallel sides.')

        ltwo=(lonev(2)*(quadxx(1)-quadxx(2))-
     $        lonev(1)*(quadyy(1)-quadyy(2)))/cross
        lone=(ltwov(2)*(quadxx(1)-quadxx(2))-
     $        ltwov(1)*(quadyy(1)-quadyy(2)))/cross
        xcen=quadxx(1)+lone*lonev(1)
        ycen=quadyy(1)+lone*lonev(2)

        r=SQRT((quadxx(4)-xcen)**2+(quadyy(4)-ycen)**2)
        IF (cross>0._r8) r=-r
        x_node=(/lxy%dx(1:poly_degree),1._r8/)

        xc=-bgridx  !  straight segment between corners 1 and 2
        DO jx=0,mx
          xc=xc+agridx*jx+bgridx
          lxy%fs(:,jx,0)=(/(xc*quadxx(2)+(1._r8-xc)*quadxx(1)),
     $                     (xc*quadyy(2)+(1._r8-xc)*quadyy(1))/)
          IF (jx==0) CYCLE
          DO ix=1,poly_degree-1
            lxy%fsh(:,ix,jx,0)=x_node(ix) *lxy%fs(:,jx,0)+
     $                  (1._r8-x_node(ix))*lxy%fs(:,jx-1,0)
          ENDDO
        ENDDO

        xc=-bgridx  !  curved segment between corners 4 and 3
        DO jx=0,mx
          xc=xc+agridx*jx+bgridx
          dir=(1._r8-xc)*lonev+xc*ltwov
          dir=dir/SQRT(SUM(dir**2))  !  direction vector to node
          lxy%fs(:,jx,my)=(/xcen,ycen/)-r*dir
          IF (jx>0) THEN
            DO ix=1,poly_degree-1
              xxyy=x_node(ix)*dir+(1._r8-x_node(ix))*dip
              xxyy=xxyy/SQRT(SUM(xxyy**2))
              lxy%fsh(:,ix,jx,my)=(/xcen,ycen/)-r*xxyy
            ENDDO
          ENDIF
          xcp=xc
          dip=dir
        ENDDO

        DO jx=0,mx  !  other nodes are interpolated between the
          yc=0      !  1-2 and 4-3 segments.
          DO jy=1,my
            yc=yc+agridy*jy+bgridy
            IF (jy<my) THEN
              lxy%fs(:,jx,jy)=
     $          (1._r8-yc)*lxy%fs(:,jx,0)+yc*lxy%fs(:,jx,my)
              IF (jx>0) THEN
                DO ix=1,poly_degree-1
                  lxy%fsh(:,ix,jx,jy)=
     $              (1._r8-yc)*lxy%fsh(:,ix,jx,0)+yc*lxy%fsh(:,ix,jx,my)
                ENDDO
              ENDIF
            ENDIF
            iint=0
            DO iy=1,poly_degree-1
              lxy%fsv(:,iy,jx,jy)=x_node(iy) *lxy%fs(:,jx,jy)+
     $                     (1._r8-x_node(iy))*lxy%fs(:,jx,jy-1)
              IF (jx>0) THEN
                DO ix=1,poly_degree-1
                  iint=iint+1
                  lxy%fsi(:,iint,jx,jy)=x_node(iy)*lxy%fsh(:,ix,jx,jy)+
     $                              (1-x_node(iy))*lxy%fsh(:,ix,jx,jy-1)
                ENDDO
              ENDIF
            ENDDO
            ycp=yc
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c     generate the grid of vertices using bilinear interpolation of the
c     corner positions.
c-----------------------------------------------------------------------
      ELSE
        xc=-bgridx
        DO jx=0,mx
          xc=xc+agridx*jx+bgridx
          yc=-bgridy
          DO jy=0,my
            yc=yc+agridy*jy+bgridy
            xy(:,jx,jy)=
     $        (/(1._r8-yc)*(xc*quadxx(2)+(1._r8-xc)*quadxx(1))+
     $                 yc *(xc*quadxx(3)+(1._r8-xc)*quadxx(4)),
     $          (1._r8-yc)*(xc*quadyy(2)+(1._r8-xc)*quadyy(1))+
     $                 yc *(xc*quadyy(3)+(1._r8-xc)*quadyy(4))/)
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       interpolate locations of side and interior nodes.
c-----------------------------------------------------------------------
        DO ibasis=1,SIZE(lxy%ix0)
          dx=lxy%dx(ibasis)
          dy=lxy%dy(ibasis)
          DO iy=lxy%iy0(ibasis),lxy%my
            iym=MAX(iy-1,0_i4)
            IF (dy==0) iym=iy
            DO ix=lxy%ix0(ibasis),lxy%mx
              ixm=MAX(ix-1,0_i4)
              IF (dx==0) ixm=ix
              xxyy=(1-dy)*((1-dx)*xy(:,ixm,iym)+dx*xy(:,ix,iym))
     $               +dy *((1-dx)*xy(:,ixm,iy )+dx*xy(:,ix,iy))
              CALL lagr_quad_basis_assign_loc(lxy,xxyy,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     set equilibrium fields to uniform values.
c-----------------------------------------------------------------------
      CALL lagr_quad_alloc(lbq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(ljq,mx,my,3_i4,poly_degree)
      CALL lagr_quad_alloc(lpq,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(lnd,mx,my,1_i4,poly_degree)
      CALL lagr_quad_alloc(ldiff,mx,my,1_i4,poly_degree)
      ljq=0
      lnd=ndens
      bex=be0*SIN(pi*thetab)*COS(pi*phib)
      bey=be0*SIN(pi*thetab)*SIN(pi*phib)
      bez=be0*COS(pi*thetab)
      xc=MINVAL(quadxx)
      DO ibasis=1,SIZE(lbq%ix0)
        DO iy=lbq%iy0(ibasis),lbq%my
          DO ix=lxy%ix0(ibasis),lxy%mx
            IF (geom=='tor') THEN
              IF (xc>0) THEN
                dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
                dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
                CALL lagr_quad_eval(lxy,dx,dy,0_i4)
                CALL lagr_quad_basis_assign_loc(lbq,
     $            (/bex*xc/lxy%f(1),bey,bez*xc/),ibasis,ix,iy)
              ELSE
                CALL lagr_quad_basis_assign_loc(lbq,
     $            (/0._r8,bey,0._r8/),ibasis,ix,iy)
              ENDIF
            ELSE
              CALL lagr_quad_basis_assign_loc(lbq,
     $          (/bex,bey,bez/),ibasis,ix,iy)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      lpq=beta*be0**2/(2._r8*mu0)
c-----------------------------------------------------------------------
c     diffusivity profile -- base on distance from (xo,yo) with 
c     scales relative to xvac.  these are not standard uses of these
c     input variables, but they are not used, otherwise, for 
c     rectquad meshes.
c-----------------------------------------------------------------------
      DO ibasis=1,SIZE(ldiff%ix0)
        DO iy=ldiff%iy0(ibasis),ldiff%my
          DO ix=ldiff%ix0(ibasis),ldiff%mx
            dx=ix-ldiff%ix0(ibasis)+ldiff%dx(ibasis)
            dy=iy-ldiff%iy0(ibasis)+ldiff%dy(ibasis)
            CALL lagr_quad_eval(lxy,dx,dy,0_i4)
            r=SQRT((lxy%f(1)-xo)**2+(lxy%f(2)-yo)**2)/(xvac+smallnum**2)
            CALL lagr_quad_basis_assign_loc(ldiff,
     $        (/(1+(SQRT(dvac)-1)*r**dexp)**2/),ibasis,ix,iy)
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find vacuum magnetic field components due to solenoids.
c-----------------------------------------------------------------------
      IF (ncoil>0.AND.geom=='tor') THEN
        bvec(3)=0._r8
        DO ibasis=1,SIZE(lxy%ix0)
          DO iy=lxy%iy0(ibasis),lxy%my
            DO ix=lxy%ix0(ibasis),lxy%mx
              dx=ix-lxy%ix0(ibasis)+lxy%dx(ibasis)
              dy=iy-lxy%iy0(ibasis)+lxy%dy(ibasis)
              CALL lagr_quad_eval(lxy,dx,dy,1_i4)
              CALL brz_eval(bvec(1:2),lxy%f(1),lxy%f(2),ncoil,
     $                      coil_r,coil_z,coil_current,mu0)
              CALL lagr_quad_basis_add_loc(lbq,bvec,ibasis,ix,iy)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rect_quad_init
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimset_init
