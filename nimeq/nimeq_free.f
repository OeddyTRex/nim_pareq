c-----------------------------------------------------------------------
c     file nimeq_free.f
c     contains a module for updating poloidal flux along the boundary
c     of the Grad-Shafranov domain when performing free-boundary
c     equilibrium computations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  nimeq_free
c     1.  nimeq_free_init
c     2.  nimeq_free_eval
c     3.  nimeq_free_iwrite
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     0. module declaration for nimeq_free
c-----------------------------------------------------------------------
      MODULE nimeq_free
      USE local
      USE nimeq_mod
      USE nim_locate
      USE vector_type_mod
      IMPLICIT NONE

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: psi_from_coil,psi_from_j
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: rzieff
      REAL(r8), DIMENSION(:), ALLOCATABLE :: ieff,wieff
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: surfpr,symieff

      INTEGER(i4), PRIVATE :: ib_rfb1,ib_rfb2,ib_vfb1,ib_vfb2
      INTEGER(i4), PRIVATE :: node1,node2,node3,node4
      REAL(r8), PRIVATE :: x_rfb1,x_rfb2,x_vfb1,x_vfb2 
      REAL(r8), PRIVATE :: y_rfb1,y_rfb2,y_vfb1,y_vfb2 
      REAL(r8), PRIVATE :: riamp,ridamp,riiamp
      REAL(r8), PRIVATE :: viamp,vidamp,viiamp
      REAL(r8), PRIVATE :: rfb_err_prev,rfb_err_int
      REAL(r8), PRIVATE :: vfb_err_prev,vfb_err_int

      REAL(r8), DIMENSION(:), ALLOCATABLE, PRIVATE :: invr2
      TYPE(vector_type), DIMENSION(:), ALLOCATABLE, PRIVATE :: icvect

      INTEGER(i4) :: fst_rfbc=-1 ! index of fst_rfbc
      INTEGER(i4) :: fst_vfbc=-1 ! index of fst_vfbc

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. nimeq_free_init creates the matrices that are used
c     to compute surface-lambda=psi/R**2 from the external coils and
c     from the toroidal current density in the interior of the domain.
c
c-PRE linear geometry
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_free_init
      USE input
      USE input_eq
      USE fields
      USE seam_storage_mod
      USE physdat
      USE mpi_nim
      USE pardata
      USE edge

      INTEGER(i4) :: ib,ibp,ix,ixp,iy,iyp,iv0,iv0p,iv,ivp,ip,ipp,
     $               dix,diy,is,icoil,icount,ncount,ijphi,njphi,ierr,
     $               ibgl,itmp,ierror
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: it1d
      INTEGER(i4) :: status(mpi_status_size)
      REAL(r8) :: x,y,dx,dy,denom,rtmp
      REAL(r8), DIMENSION(2) :: rzsurf
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xnode
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: r2d,r2d2,rzsarr
      REAL(r8), PARAMETER :: xy_buff=5.e-2_r8
      REAL(r8), EXTERNAL :: aph_eval
c-----------------------------------------------------------------------
c     use seam0 to determine how many surface values of lambda are
c     needed.
c-----------------------------------------------------------------------
      ALLOCATE(xnode(0:poly_degree))
      CALL poly_nodes(poly_degree,xnode)
      ncount=0
      iv0p=seam0%nvert
      DO iv0=1,seam0%nvert
        ip1: DO ip=1,SIZE(seam0%vertex(iv0)%ptr,2)
          ibgl=seam0%vertex(iv0)%ptr(1,ip)
          iv  =seam0%vertex(iv0)%ptr(2,ip)
          ip2: DO ipp=1,SIZE(seam0%vertex(iv0p)%ptr,2)
            ibp=seam0%vertex(iv0p)%ptr(1,ipp)
            ivp=seam0%vertex(iv0p)%ptr(2,ipp)
            IF (ibp==ibgl) EXIT ip1
          ENDDO ip2
        ENDDO ip1
        IF (ibp/=ibgl) CALL nim_stop("Nimeq_free_init: prev not found.")
        IF (block2proc(ibgl)==node) THEN
          ib=global2local(ibgl)
          IF (seam(ib)%expoint(iv)) ncount=ncount+1
          IF (seam(ib)%expoint(iv ).AND.
     $        seam(ib)%expoint(ivp).AND.ib<=nrbl)
     $      ncount=ncount+poly_degree-1
        ENDIF
        iv0p=iv0
      ENDDO
      IF (nprocs>1) THEN
        CALL mpi_allreduce(ncount,itmp,1,mpi_nim_int,mpi_sum,
     $       mpi_comm_world,ierr)
        ncount=itmp
      ENDIF
c-----------------------------------------------------------------------
c     create the matrix for the external coil contributions.  the
c     surface index is first.  the icvect vector type is used to save
c     the surface-node row index for mapping back to nimrod data.
c-----------------------------------------------------------------------
      ALLOCATE(psi_from_coil(ncount,ncoil_tot),r2d(ncount,ncoil_tot))
      ALLOCATE(rzsarr(2,ncount),r2d2(2,ncount))
      ALLOCATE(surfpr(ncount),it1d(ncount))
      ALLOCATE(icvect(nbl))
      psi_from_coil=0._r8
      rzsarr=0._r8
      surfpr=0
      DO ib=1,nbl
        IF (ib<=nrbl) THEN
          CALL vector_type_alloc(icvect(ib),poly_degree,rb(ib)%mx,
     $                           rb(ib)%my,1_i4)
        ELSE
          CALL vector_type_alloc(icvect(ib),1_i4,tb(ib)%mvert,0_i4,1_i4)
        ENDIF
        icvect(ib)=0
      ENDDO

      iv0p=seam0%nvert
      icount=0
      DO iv0=1,seam0%nvert
        ip3: DO ip=1,SIZE(seam0%vertex(iv0)%ptr,2)
          ibgl=seam0%vertex(iv0)%ptr(1,ip)
          iv  =seam0%vertex(iv0)%ptr(2,ip)
          ip4: DO ipp=1,SIZE(seam0%vertex(iv0p)%ptr,2)
            ibp=seam0%vertex(iv0p)%ptr(1,ipp)
            ivp=seam0%vertex(iv0p)%ptr(2,ipp)
            IF (ibp==ibgl) EXIT ip3
          ENDDO ip4
        ENDDO ip3
        IF (block2proc(ibgl)==node) THEN
          ib=global2local(ibgl)
          IF (seam(ib)%expoint(iv).AND.seam(ib)%expoint(ivp).AND.
     $        ib<=nrbl.AND.poly_degree>1) THEN
            ix=seam(ib)%segment(iv)%intxys(1)
            iy=seam(ib)%segment(iv)%intxys(2)
            dix=ABS(seam(ib)%segment(iv)%intxyn(1)-
     $              seam(ib)%segment(iv)%intxyp(1))
            diy=ABS(seam(ib)%segment(iv)%intxyn(2)-
     $              seam(ib)%segment(iv)%intxyp(2))
            DO is=1,poly_degree-1
              icount=icount+1
              x=ix+dix*(xnode(is)-1._r8)
              y=iy+diy*(xnode(is)-1._r8)
              CALL lagr_quad_eval(rb(ib)%rz,x,y,0_i4)
              rzsarr(:,icount)=rb(ib)%rz%f
              surfpr(icount)=node
              DO icoil=1,ncoil_tot
                psi_from_coil(icount,icoil)=
     $            -aph_eval(rb(ib)%rz%f(1),rb(ib)%rz%f(2),1_i4,
     $                      allcoil_r(icoil),allcoil_z(icoil),1._r8,mu0)
              ENDDO
              IF (seam(ib)%segment(iv)%h_side) THEN
                icvect(ib)%arrh(1,is,ix,iy)=icount
              ELSE
                icvect(ib)%arrv(1,is,ix,iy)=icount
              ENDIF
            ENDDO
          ENDIF
          IF (seam(ib)%expoint(iv)) THEN
            ix=seam(ib)%vertex(iv)%intxy(1)
            iy=seam(ib)%vertex(iv)%intxy(2)
            icount=icount+1
            surfpr(icount)=node
            icvect(ib)%arr(1,ix,iy)=icount
            IF (ib<=nrbl) THEN
              CALL lagr_quad_eval(rb(ib)%rz,REAL(ix,r8),REAL(iy,r8),
     $                            0_i4)
              rzsarr(:,icount)=rb(ib)%rz%f
              DO icoil=1,ncoil_tot
                psi_from_coil(icount,icoil)=
     $            -aph_eval(rb(ib)%rz%f(1),rb(ib)%rz%f(2),1_i4,
     $                      allcoil_r(icoil),allcoil_z(icoil),1._r8,mu0)
              ENDDO
            ELSE
              rzsarr(:,icount)=(/tb(ib)%tgeom%xs(ix),
     $                           tb(ib)%tgeom%ys(ix)/)
              DO icoil=1,ncoil_tot
                psi_from_coil(icount,icoil)=
     $            -aph_eval(tb(ib)%tgeom%xs(ix),
     $                      tb(ib)%tgeom%ys(ix),1_i4,
     $                      allcoil_r(icoil),allcoil_z(icoil),1._r8,mu0)
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        IF (nprocs>1) THEN
          CALL mpi_allreduce(icount,itmp,1,mpi_nim_int,mpi_max,
     $         mpi_comm_world,ierr)
          icount=itmp
        ENDIF
        iv0p=iv0
      ENDDO
c-----------------------------------------------------------------------
c     only one process should complete a matrix entry in the above
c     looping, so summing assembles the correct array on all processors.
c-----------------------------------------------------------------------
      IF (nprocs>1) THEN
        CALL mpi_allreduce(psi_from_coil,r2d,ncount*ncoil_tot,
     $       mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
        psi_from_coil=r2d
        CALL mpi_allreduce(rzsarr,r2d2,2*ncount,
     $       mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
        rzsarr=r2d2
        CALL mpi_allreduce(surfpr,it1d,ncount,
     $       mpi_nim_int,mpi_sum,mpi_comm_world,ierr)
        surfpr=it1d
      ENDIF
      DEALLOCATE(r2d,r2d2,it1d)
c-----------------------------------------------------------------------
c     in the above looping, only the first matched block receives the
c     edge-node row index, so network this information.
c-----------------------------------------------------------------------
      DO ib=1,nbl
        CALL edge_load_arr(icvect(ib),1_i4,poly_degree-1_i4,seam(ib))
      ENDDO
      CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.false.)
      DO ib=1,nbl
        CALL edge_unload_arr(icvect(ib),1_i4,poly_degree-1_i4,
     $                       seam(ib))
      ENDDO
c-----------------------------------------------------------------------
c     find the normalized amplification and logical location of the
c     poloidal flux sensors for radial feedback.
c-----------------------------------------------------------------------
      node1=-1
      node2=-1
      node3=-1
      node4=-1
      IF (fst_rfbc>0) THEN
        CALL sensor_find(rfb_r1,rfb_z1,'first radial',node1,ib_rfb1,
     $                   x_rfb1,y_rfb1)
        CALL sensor_find(rfb_r2,rfb_z2,'second radial',node2,ib_rfb2,
     $                   x_rfb2,y_rfb2)
c-----------------------------------------------------------------------
c       calculate the normalized radial feedback amplification factors.
c       ramp = -rfb_amp /(G(x1, x_coils) - G(x2, x_coils)) ;
c       aph_eval returns negative G .
c-----------------------------------------------------------------------
        denom=aph_eval(rfb_r1,rfb_z1,nrfbcoil,rfbcoil_r,
     $                    rfbcoil_z,rfbcoil_i,mu0)-
     $        aph_eval(rfb_r2,rfb_z2,nrfbcoil,rfbcoil_r,
     $                    rfbcoil_z,rfbcoil_i,mu0)
c-DEBUG
c       IF (node==0) WRITE(*,*) "rfb denom",denom
        IF (denom==0) denom=1._r8
        riamp=rfb_amp/denom
        ridamp=rfb_damp/denom
        riiamp=rfb_iamp/denom
        rfb_err_prev=0._r8
        rfb_err_int=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     find the normalized amplification and logical location of the
c     poloidal flux sensors for vertical feedback.
c-----------------------------------------------------------------------
      IF (fst_vfbc>0) THEN
        CALL sensor_find(vfb_r1,vfb_z1,'first vertical',node3,ib_vfb1,
     $                   x_vfb1,y_vfb1)
        CALL sensor_find(vfb_r2,vfb_z2,'second vertical',node4,ib_vfb2,
     $                   x_vfb2,y_vfb2)
c-----------------------------------------------------------------------
c       calculate the normalized vertical feedback amplification
c       factors.
c-----------------------------------------------------------------------
        denom = (aph_eval(vfb_r1,vfb_z1,nvfbcoil,vfbcoil_r,
     $                    vfbcoil_z,vfbcoil_i,mu0)-
     $           aph_eval(vfb_r2,vfb_z2,nvfbcoil,vfbcoil_r,
     $                    vfbcoil_z,vfbcoil_i,mu0))
c-DEBUG
c       IF (node==0) WRITE(*,*) "vfb denom",denom
        IF (denom==0) denom=1._r8
        viamp=vfb_amp/denom
        vidamp=vfb_damp/denom
        viiamp=vfb_iamp/denom
c-DEBUG
c        WRITE(*,*) "node,viamp,vidamp,viiamp"
c        WRITE(*,*) node,viamp,vidamp,viiamp
        vfb_err_prev=0._r8
        vfb_err_int=0._r8
      ENDIF
c-----------------------------------------------------------------------
c     the next operations create the matrix to generate the surface
c     flux values from the internal toroidal current density.  the
c     internal current density is effectively treated as current loops
c     at the quadrature points of the finite-element integration
c     that is used for the finite-element projections.
c-----------------------------------------------------------------------
      njphi=0
      DO ib=1,nrbl
        njphi=njphi+SIZE(rb(ib)%wjac)
      ENDDO
      DO ib=nrbl+1,nbl
        njphi=njphi+tb(ib)%mcell*tb(ib)%ng
      ENDDO
      ALLOCATE(psi_from_j(ncount,njphi))
      ALLOCATE(invr2(ncount))
c-----------------------------------------------------------------------
c     the position of each quadrature point and the final effective
c     currents will be written for initializing external regions.
c     create space for this information and initialize it.
c-----------------------------------------------------------------------
      CALL ieff_init
c-----------------------------------------------------------------------
c     loop over all nodes along the surface and find the contribution
c     to poloidal flux from the normalized current at each quadrature
c     point.  here, normalized current means per unit of mu0*J_phi/R.
c-----------------------------------------------------------------------
      DO icount=1,ncount
        invr2(icount)=1._r8/rzsarr(1,icount)**2
        CALL psi_row(icount,rzsarr(:,icount))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

      CONTAINS

c-----------------------------------------------------------------------
c       internal subroutine to complete a row of the psi_from_j
c       matrix.  matrix elements will be multiplied by +mu0*J_phi/R
c       in nimeq_free_eval.
c-----------------------------------------------------------------------
        SUBROUTINE psi_row(irow,rzsf)

        INTEGER(i4), INTENT(IN) :: irow
        REAL(r8), DIMENSION(2), INTENT(IN) :: rzsf

        INTEGER(i4) :: jb,jg,jx,jy,jcol,jelm
        REAL(r8) :: d2
        REAL(r8) :: dmin=1.e-8_r8

        jcol=0
        DO jb=1,nrbl
          jelm=1
          DO jy=0,rb(jb)%my-1
            DO jx=0,rb(jb)%mx-1
              DO jg=1,rb(jb)%ng
                jcol=jcol+1
                CALL lagr_quad_eval(rb(jb)%rz,jx+rb(jb)%xg(jg),
     $                              jy+rb(jb)%yg(jg),0_i4)
                d2=SUM((rzsf-rb(jb)%rz%f)**2)
                IF (d2<dmin*rb(jb)%jac2d(jg,jelm)) THEN
                  psi_from_j(irow,jcol)=0._r8
                  CYCLE
                ENDIF
                IF (symm_region/="neither".OR.free_symmj) THEN
                  d2=(rzsf(1)-rb(jb)%rz%f(1))**2+
     $               (rzsf(2)+rb(jb)%rz%f(2))**2
                  IF (d2<dmin*rb(jb)%jac2d(jg,jelm)) THEN
                    psi_from_j(irow,jcol)=0._r8
                    CYCLE
                  ENDIF
                ENDIF
                IF (rb(jb)%rz%f(2)>0._r8) THEN
                  IF (symm_region=="top".OR.free_symmj) THEN
                    psi_from_j(irow,jcol)=
     $                -(aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                           rb(jb)%rz%f(2),1._r8,1._r8)+
     $                  aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                           -rb(jb)%rz%f(2),1._r8,1._r8))*
     $                 rb(jb)%wjac(jg,jelm)
                  ELSE IF (symm_region=="bottom") THEN
                    psi_from_j(irow,jcol)=0._r8
                  ELSE
                    psi_from_j(irow,jcol)=
     $                -aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                          rb(jb)%rz%f(2),1._r8,1._r8)*
     $                 rb(jb)%wjac(jg,jelm)
                  ENDIF
                ELSE IF (rb(jb)%rz%f(2)==0._r8) THEN
                  psi_from_j(irow,jcol)=
     $              -aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                        rb(jb)%rz%f(2),1._r8,1._r8)*
     $               rb(jb)%wjac(jg,jelm)
                ELSE
                  IF (symm_region=="top".OR.free_symmj) THEN
                    psi_from_j(irow,jcol)=0._r8
                  ELSE IF (symm_region=="bottom") THEN
                    psi_from_j(irow,jcol)=
     $                -(aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                           rb(jb)%rz%f(2),1._r8,1._r8)+
     $                  aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                           -rb(jb)%rz%f(2),1._r8,1._r8))*
     $                 rb(jb)%wjac(jg,jelm)
                  ELSE
                    psi_from_j(irow,jcol)=
     $                -aph_eval(rzsf(1),rzsf(2),1_i4,rb(jb)%rz%f(1),
     $                          rb(jb)%rz%f(2),1._r8,1._r8)*
     $                 rb(jb)%wjac(jg,jelm)
                  ENDIF
                ENDIF
              ENDDO
              jelm=jelm+1
            ENDDO
          ENDDO
        ENDDO

c-PRE
        DO jb=nrbl+1,nbl
        ENDDO

        RETURN
        END SUBROUTINE psi_row

c-----------------------------------------------------------------------
c       this internal subroutine loops over quad points in the same
c       order as psi_row and saves information that is needed for
c       writing quad-point currents as individual coils.
c         rzieff= poloidal position
c         wieff = weighting factor to get I from mu0*Jphi/R at quad pts
c         symieff= 1 --  active location, no symmetry
c         symieff=-1 --  active location, along with is up-down image
c         symieff= 0 --  inactive location, no contribution
c-----------------------------------------------------------------------
        SUBROUTINE ieff_init

        INTEGER(i4) :: jb,jg,jx,jy,jcol,jelm

        ALLOCATE(rzieff(2,njphi),ieff(njphi),symieff(njphi),
     $           wieff(njphi))

        jcol=0
        DO jb=1,nrbl
          jelm=1
          DO jy=0,rb(jb)%my-1
            DO jx=0,rb(jb)%mx-1
              DO jg=1,rb(jb)%ng
                jcol=jcol+1
                CALL lagr_quad_eval(rb(jb)%rz,jx+rb(jb)%xg(jg),
     $                              jy+rb(jb)%yg(jg),0_i4)
                rzieff(:,jcol)=rb(jb)%rz%f  !  save for output
                IF (rb(jb)%rz%f(2)>0._r8) THEN
                  IF (symm_region=="top".OR.free_symmj) THEN
                    symieff(jcol)=-1_i4
                  ELSE IF (symm_region=="bottom") THEN
                    symieff(jcol)=0_i4
                  ELSE
                    symieff(jcol)=1_i4
                  ENDIF
                ELSE IF (rb(jb)%rz%f(2)==0._r8) THEN
                  symieff(jcol)=1_i4
                ELSE
                  IF (symm_region=="top".OR.free_symmj) THEN
                    symieff(jcol)=0_i4
                  ELSE IF (symm_region=="bottom") THEN
                    symieff(jcol)=-1_i4
                  ELSE
                    symieff(jcol)=1_i4
                  ENDIF
                ENDIF
                wieff(jcol)=ABS(symieff(jcol))*rb(jb)%wjac(jg,jelm)/mu0
              ENDDO
              jelm=jelm+1
            ENDDO
          ENDDO
        ENDDO

c-PRE
        DO jb=nrbl+1,nbl
        ENDDO

        RETURN
        END SUBROUTINE ieff_init

c-----------------------------------------------------------------------
c       this internal subroutine uses nim_locate routines to find
c       the block and logical coordinates of the feedback sensor points.
c-----------------------------------------------------------------------
        SUBROUTINE sensor_find(rpos,zpos,slabel,snode,ibpos,xpos,ypos)

        REAL(r8), INTENT(IN) :: rpos,zpos
        CHARACTER(*), INTENT(IN) :: slabel
        REAL(r8), INTENT(OUT) :: xpos,ypos
        INTEGER(i4), INTENT(OUT) :: snode,ibpos

        INTEGER(i4) :: itmp,ierror
c-----------------------------------------------------------------------
c       determine and communicate the node that has the feedback point
c       in its rblock.
c-----------------------------------------------------------------------
        snode=0
        DO ibgl=1,nrbl_total
          IF (block2proc(ibgl)==node) THEN
            ib=global2local(ibgl) 
            CALL nim_rb_contrz(rpos,zpos,rb(ib),seam(ib),ibpos,iv)
            IF (ibpos==0) THEN
              CALL nim_stop(slabel//'feedback point outside domain.')
            ELSE IF (ibpos==ibgl) THEN
              ibpos=global2local(ibpos)
              snode=node
              EXIT 
            ENDIF
          ENDIF
        ENDDO 
        IF (nprocs>1) THEN
          CALL mpi_allreduce(snode,itmp,1,mpi_nim_int,mpi_sum,
     $                       mpi_comm_world,ierror)
          snode=itmp
        ENDIF
c-----------------------------------------------------------------------
c       find the logical coordinates within the block.
c-----------------------------------------------------------------------
        IF (node==snode) THEN
          CALL nim_rb_locate(rpos,zpos,rpos,rb(ibpos),xpos,ypos,ierr)
c-DEBUG
c         WRITE(*,'(a,3i3,4es11.4)') "node ib ierr rpos zpos xpos ypos",
c    $      node,ibpos,ierr,rpos,zpos,xpos,ypos
          IF (ierr/=0.OR.xpos>rb(ibpos)%mx+xy_buff.OR.
     $                   ypos>rb(ibpos)%my+xy_buff.OR.
     $                   xpos<-xy_buff.OR.ypos<-xy_buff)
     $      CALL nim_stop('Error finding '//slabel//' feedback point.')
        ENDIF

        RETURN
        END SUBROUTINE sensor_find

      END SUBROUTINE nimeq_free_init
c-----------------------------------------------------------------------
c     subprogram 2. nimeq_free_eval evaluates boundary contributions to
c     the poloidal flux function that result from external coils and
c     from internal current density, which is also treated as discrete
c     wires that are located at the quadrature points within elements.
c     this is used within the nonlinear iteration, so the lambda=psi/R^2
c     field is computed.
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_free_eval(lamvec,ipsi,poly_degree,frcn,pmax,pmin,
     $                           dcoil_fb)
      USE input_eq
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata

      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: lamvec
      INTEGER(i4), INTENT(IN) :: ipsi,poly_degree
      REAL(r8), INTENT(IN) :: frcn
      REAL(r8), INTENT(OUT) :: pmax,pmin
      LOGICAL, INTENT(IN), OPTIONAL :: dcoil_fb

      INTEGER(i4) :: ib,ibp,ix,iy,iv0,iv0p,iv,ivp,ip,ipp,np,
     $               icoil,icount,ns,jb,jg,jx,jy,jcol,jelm,ncount,ierr
      REAL(r8), DIMENSION(SIZE(psi_from_j,1)) :: spsi,stmp
      REAL(r8) :: omgsc,bz,jac,dxdr,dydr,rtmp
      REAL(r8) :: psi1,psi2,this_err
      LOGICAL :: coil_fb
c-----------------------------------------------------------------------
c     the mu0*J_phi/R values are in the quadrature-point storage
c     qja_eq(3,:,:).  the rows of the psi_from_j matrix are organized
c     according to this storage.
c
c     the distinct quadrature points owned by each process contribute to
c     all surface nodes, so sum after this operation.
c-----------------------------------------------------------------------
      IF (PRESENT(dcoil_fb)) THEN
        coil_fb=dcoil_fb
      ELSE
        coil_fb=.false.
      ENDIF
      ncount=SIZE(psi_from_j,1)
      spsi=0._r8
      jcol=0
      DO jb=1,nrbl
        jelm=1
        DO jy=0,rb(jb)%my-1
          DO jx=0,rb(jb)%mx-1
            DO jg=1,rb(jb)%ng
              jcol=jcol+1
              spsi=spsi+psi_from_j(:,jcol)*rb(jb)%qja_eq%qpf(3,jg,jelm)
              ieff(jcol)=wieff(jcol)*rb(jb)%qja_eq%qpf(3,jg,jelm)
            ENDDO
            jelm=jelm+1
          ENDDO
        ENDDO
      ENDDO
c-PRE
      DO jb=nrbl+1,nbl
      ENDDO

      IF (nprocs>1) THEN
        CALL mpi_allreduce(spsi,stmp,ncount,
     $       mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
        spsi=stmp
      ENDIF
c-----------------------------------------------------------------------
c     if radial feedback is used, calculate the change in the
c     radial feedback coil currents
c-----------------------------------------------------------------------
      IF (coil_fb) THEN
        IF (fst_rfbc>0) THEN
          IF (node==node1) THEN
            CALL lagr_quad_eval(rb(ib_rfb1)%rwork2,x_rfb1,y_rfb1,0_i4)
            psi1=rb(ib_rfb1)%rwork2%f(1)*rfb_r1**2
          ELSE 
            psi1=0._r8 
          ENDIF
          IF (node==node2) THEN
            CALL lagr_quad_eval(rb(ib_rfb2)%rwork2,x_rfb2,y_rfb2,0_i4)
            psi2=rb(ib_rfb2)%rwork2%f(1)*rfb_r2**2
          ELSE 
            psi2=0._r8
          ENDIF
          IF (nprocs>1) THEN
            CALL mpi_allreduce(psi1,rtmp,1,
     $           mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
            psi1=rtmp
            CALL mpi_allreduce(psi2,rtmp,1,
     $           mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
            psi2=rtmp
          ENDIF
          
          this_err=psi1-psi2
          rfb_err_int=rfb_err_int+this_err

          rfbc_curr=rfbc_curr+riamp*this_err+
     $        ridamp*(this_err-rfb_err_prev) + riiamp*rfb_err_int
          rfb_err_prev=this_err
c-DEBUG
c          WRITE(*,*) "New rfb curr ",rfbc_curr
c          WRITE(*,*) "psi1,psi2 ",psi1,psi2
          allcoil_i(fst_rfbc:fst_rfbc+nrfbcoil-1)=
     $      rfbcoil_i(1:nrfbcoil)*rfbc_curr
        ENDIF
c-----------------------------------------------------------------------
c       if vertical feedback is used, calculate the change in the
c       vertical feedback coil currents
c-----------------------------------------------------------------------
        IF (fst_vfbc>0) THEN
          IF (node==node3) THEN
            CALL lagr_quad_eval(rb(ib_vfb1)%rwork2,x_vfb1,y_vfb1,0_i4)
            psi1=rb(ib_vfb1)%rwork2%f(1)*vfb_r1**2
          ELSE
            psi1=0._r8
          ENDIF
          IF (node==node4) THEN
            CALL lagr_quad_eval(rb(ib_vfb2)%rwork2,x_vfb2,y_vfb2,0_i4)
            psi2=rb(ib_vfb2)%rwork2%f(1)*vfb_r2**2
          ELSE
            psi2=0._r8
          ENDIF
          IF (nprocs>1) THEN
            CALL mpi_allreduce(psi1,rtmp,1,
     $           mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
            psi1=rtmp
            CALL mpi_allreduce(psi2,rtmp,1,
     $           mpi_nim_real,mpi_sum,mpi_comm_world,ierr)
            psi2=rtmp
          ENDIF
          
          this_err=psi1-psi2
          vfb_err_int=vfb_err_int+this_err

          vfbc_curr=vfbc_curr+viamp*this_err+
     $        vidamp*(this_err-vfb_err_prev)+viiamp*vfb_err_int
          vfb_err_prev=this_err
c-DEBUG
c          WRITE(*,*) "New vfb curr ",vfbc_curr
c          WRITE(*,*) "psi1,psi2 ",psi1,psi2
          allcoil_i(fst_vfbc:fst_vfbc+nvfbcoil-1)=
     $      vfbcoil_i(1:nvfbcoil)*vfbc_curr
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     compute the external coil contributions from the previously
c     created matrix.
c-----------------------------------------------------------------------
      DO icoil=1,ncoil_tot
        spsi=spsi+psi_from_coil(:,icoil)*allcoil_i(icoil)
      ENDDO
c-----------------------------------------------------------------------
c     apply the relaxation factor.
c-----------------------------------------------------------------------
      spsi=spsi*frcn
      omgsc=1._r8-frcn
c-----------------------------------------------------------------------
c     add these contributions to the respective locations in the passed
c     data structure.  also determine the new extrema for the open-
c     flux range.
c-----------------------------------------------------------------------
      iv0p=seam0%nvert
      icount=1
      ns=poly_degree-1
      pmin= HUGE(pmin)
      pmax=-HUGE(pmax)
      DO ib=1,nbl
        ivp=seam(ib)%nvert
        DO iv=1,seam(ib)%nvert
          IF (seam(ib)%expoint(iv).AND.seam(ib)%expoint(ivp).AND.
     $        ib<=nrbl.AND.poly_degree>1) THEN
            ix=seam(ib)%segment(iv)%intxys(1)
            iy=seam(ib)%segment(iv)%intxys(2)
            IF (seam(ib)%segment(iv)%h_side) THEN
              icount=NINT(icvect(ib)%arrh(1,1,ix,iy))
              lamvec(ib)%arrh(ipsi,1:ns,ix,iy)=
     $          omgsc*lamvec(ib)%arrh(ipsi,1:ns,ix,iy)+
     $          spsi(icount:icount+ns-1)*invr2(icount:icount+ns-1)
              pmax=MAX(pmax,MAXVAL(lamvec(ib)%arrh(ipsi,1:ns,ix,iy)/
     $                             invr2(icount:icount+ns-1)))
              pmin=MIN(pmin,MINVAL(lamvec(ib)%arrh(ipsi,1:ns,ix,iy)/
     $                             invr2(icount:icount+ns-1)))
            ELSE
              icount=NINT(icvect(ib)%arrv(1,1,ix,iy))
              lamvec(ib)%arrv(ipsi,1:ns,ix,iy)=
     $          omgsc*lamvec(ib)%arrv(ipsi,1:ns,ix,iy)+
     $          spsi(icount:icount+ns-1)*invr2(icount:icount+ns-1)
              pmax=MAX(pmax,MAXVAL(lamvec(ib)%arrv(ipsi,1:ns,ix,iy)/
     $                             invr2(icount:icount+ns-1)))
              pmin=MIN(pmin,MINVAL(lamvec(ib)%arrv(ipsi,1:ns,ix,iy)/
     $                             invr2(icount:icount+ns-1)))
            ENDIF
          ENDIF
          IF (seam(ib)%expoint(iv)) THEN
            ix=seam(ib)%vertex(iv)%intxy(1)
            iy=seam(ib)%vertex(iv)%intxy(2)
            icount=NINT(icvect(ib)%arr(1,ix,iy))
            lamvec(ib)%arr(ipsi,ix,iy)=omgsc*lamvec(ib)%arr(ipsi,ix,iy)+
     $                                 spsi(icount)*invr2(icount)
            pmax=MAX(pmax,lamvec(ib)%arr(ipsi,ix,iy)/invr2(icount))
            pmin=MIN(pmin,lamvec(ib)%arr(ipsi,ix,iy)/invr2(icount))
          ENDIF
          ivp=iv
        ENDDO
      ENDDO

      IF (nprocs>1) THEN
        CALL mpi_allreduce(pmax,rtmp,1,mpi_nim_real,mpi_max,
     $       mpi_comm_world,ierr)
        pmax=rtmp
        CALL mpi_allreduce(pmin,rtmp,1,mpi_nim_real,mpi_min,
     $       mpi_comm_world,ierr)
        pmin=rtmp
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nimeq_free_eval
c-----------------------------------------------------------------------
c     subprogram 3. nimeq_free_iwrite outputs the quad-point
c     current-density values, weighted to make current filaments, 
c     and their positions to a binary file.
c-----------------------------------------------------------------------
      SUBROUTINE nimeq_free_iwrite
      USE fields
      USE mpi_nim
      USE pardata

      INTEGER(i4) ii,ipr,ierror
c-----------------------------------------------------------------------
c     have one processor write at a time.  each line of the data
c     is self-contained, so order is not important.
c-----------------------------------------------------------------------
      DO ipr=0,nprocs-1
        IF (node==ipr) THEN
          IF (ipr==0) THEN
            CALL open_bin(binary_unit,"nimeq_iint.bin",
     $                    "UNKNOWN","REWIND",64_i4)
          ELSE
            CALL open_bin(binary_unit,"nimeq_iint.bin","OLD",
     $                    "APPEND",64_i4)
          ENDIF
          DO ii=1,SIZE(ieff)
            SELECT CASE(symieff(ii))
            CASE(1_i4)
              WRITE(binary_unit) rzieff(:,ii),ieff(ii)
            CASE(-1_i4)
              WRITE(binary_unit) rzieff(:,ii),ieff(ii)
              WRITE(binary_unit) rzieff(1,ii),-rzieff(2,ii),ieff(ii)
            END SELECT
          ENDDO
          CALL close_bin(binary_unit,"nimeq_iint.bin")
        ENDIF
        CALL mpi_barrier(mpi_comm_world,ierror)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE nimeq_free_iwrite
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimeq_free
