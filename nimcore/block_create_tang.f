c-----------------------------------------------------------------------
c     file block_create_tang.f
c     create Cartesian components of surface tangents for all blocks.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  block_create_tang.
c     2.  block_dealloc_tang.
c-----------------------------------------------------------------------
c     subprogram 1. block_create_tang
c     create the norm and tang arrays for blocks along the domain
c     boundary.
c-----------------------------------------------------------------------
      SUBROUTINE block_create_tang(poly_degree)
      USE local
      USE seam_storage_mod
      USE fields
      USE edge
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: poly_degree

      REAL(r8) :: tempm,tempp,dx,dy
      REAL(r8), DIMENSION(2) :: tempv
      REAL(r8), PARAMETER :: offset=1.e-12_r8
      INTEGER(i4) :: ib,ibe,iv,mx,my,ix,iy,minblock,i,ixp,ixm,is
      INTEGER(i4), DIMENSION(:), ALLOCATABLE  :: ivp, ivm
c-----------------------------------------------------------------------
c     zero out the seam input array of all blocks.
c-----------------------------------------------------------------------
      DO ib=1,nbl
        DO iv=1,seam(ib)%nvert
          seam(ib)%vertex(iv)%seam_in=0._r8
          NULLIFY(seam(ib)%vertex(iv)%norm)
          NULLIFY(seam(ib)%vertex(iv)%tang)
          NULLIFY(seam(ib)%segment(iv)%norm)
          NULLIFY(seam(ib)%segment(iv)%tang)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find external boundaries and compute surface tangents.
c     note:  overlap between blocks will be averaged and
c     renormalized in edge_init.
c-----------------------------------------------------------------------
      block: DO ibe=1,SIZE(exblock_list)
        ib=exblock_list(ibe)
        ALLOCATE(ivp(seam(ib)%nvert))
        ALLOCATE(ivm(seam(ib)%nvert))
        ivp=(/(i,i=2,seam(ib)%nvert),1_i4/)
        ivm=(/seam(ib)%nvert,(i,i=1,seam(ib)%nvert-1_i4)/)
c-----------------------------------------------------------------------
c       tblocks:
c-PRE   side centered vectors will be needed here, too.
c-----------------------------------------------------------------------
        IF(ib>nrbl)THEN
          verttri: DO iv=1,seam(ib)%nvert
            IF (seam(ib)%expoint(iv)) THEN
              ix=seam(ib)%vertex(iv)%intxy(1)
              iy=seam(ib)%vertex(iv)%intxy(2)
              ixp=seam(ib)%vertex(ivp(iv))%intxy(1)
              ixm=seam(ib)%vertex(ivm(iv))%intxy(1)
              tempm=SQRT((tb(ib)%tgeom%xs(ix)-tb(ib)%tgeom%xs(ixm))**2
     &              +(tb(ib)%tgeom%ys(ix)-tb(ib)%tgeom%ys(ixm))**2)
              tempp=SQRT((tb(ib)%tgeom%xs(ix)-tb(ib)%tgeom%xs(ixp))**2
     &              +(tb(ib)%tgeom%ys(ix)-tb(ib)%tgeom%ys(ixp))**2)
              tempv(1)=(tb(ib)%tgeom%xs(ix)-tb(ib)%tgeom%xs(ixm))/tempm+
     &                 (tb(ib)%tgeom%xs(ixp)-tb(ib)%tgeom%xs(ix))/tempp
              tempv(2)=(tb(ib)%tgeom%ys(ix)-tb(ib)%tgeom%ys(ixm))/tempm+
     &                 (tb(ib)%tgeom%ys(ixp)-tb(ib)%tgeom%ys(ix))/tempp
              CALL block_crnormalize(tempv)
              seam(ib)%vertex(iv)%seam_in(1:2)=tempv
            ENDIF
          ENDDO verttri
c-----------------------------------------------------------------------
c       rblocks:
c-----------------------------------------------------------------------
        ELSE
          mx=rb(ib)%mx
          my=rb(ib)%my
c-----------------------------------------------------------------------
c         bottom:
c-----------------------------------------------------------------------
          vertb: DO iv=1,mx
            IF (seam(ib)%expoint(iv)) THEN
              ix=seam(ib)%vertex(iv)%intxy(1)
              iy=seam(ib)%vertex(iv)%intxy(2)
              tempv=0
              IF (seam(ib)%expoint(ivm(iv))) THEN
                dx=ix-offset
                dy=iy
                CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                tempv=rb(ib)%rz%fx
                CALL block_crnormalize(tempv)
              ENDIF
              IF (seam(ib)%expoint(ivp(iv))) THEN
                IF (iv<mx) THEN
                  dx=ix+offset
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fx)
                  tempv=tempv+rb(ib)%rz%fx
                ELSE
                  dx=ix
                  dy=iy+offset
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fy)
                  tempv=tempv+rb(ib)%rz%fy
                ENDIF
              ENDIF
              CALL block_crnormalize(tempv)
              seam(ib)%vertex(iv)%seam_in(1:2)=tempv
              IF (seam(ib)%expoint(ivm(iv)).AND.poly_degree>1) THEN
                ALLOCATE(seam(ib)%segment(iv)%tang(2,poly_degree-1))
                ALLOCATE(seam(ib)%segment(iv)%norm(2,poly_degree-1))
                DO is=1,rb(ib)%be%n_side
                  ix=seam(ib)%segment(iv)%intxys(1)
                  iy=seam(ib)%segment(iv)%intxys(2)
                  dx=ix-1+rb(ib)%be%dx(1+is)
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fx)
                  tempv=rb(ib)%rz%fx
                  seam(ib)%segment(iv)%tang(:,is)=tempv
                  seam(ib)%segment(iv)%norm(:,is)=(/tempv(2),-tempv(1)/)
                ENDDO
              ENDIF
            ENDIF
          ENDDO vertb
c-----------------------------------------------------------------------
c         right side:
c-----------------------------------------------------------------------
          vertr: DO iv=mx+1,mx+my
            IF (seam(ib)%expoint(iv)) THEN
              ix=seam(ib)%vertex(iv)%intxy(1)
              iy=seam(ib)%vertex(iv)%intxy(2)
              tempv=0
              IF (seam(ib)%expoint(ivm(iv))) THEN
                dx=ix
                dy=iy-offset
                CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                tempv=rb(ib)%rz%fy
                CALL block_crnormalize(tempv)
              ENDIF
              IF (seam(ib)%expoint(ivp(iv))) THEN
                IF (iv<mx+my) THEN
                  dx=ix
                  dy=iy+offset
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fy)
                  tempv=tempv+rb(ib)%rz%fy
                ELSE
                  dx=ix-offset
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fx)
                  tempv=tempv-rb(ib)%rz%fx
                ENDIF
              ENDIF
              CALL block_crnormalize(tempv)
              seam(ib)%vertex(iv)%seam_in(1:2)=tempv
              IF (seam(ib)%expoint(ivm(iv)).AND.poly_degree>1) THEN
                ALLOCATE(seam(ib)%segment(iv)%tang(2,poly_degree-1))
                ALLOCATE(seam(ib)%segment(iv)%norm(2,poly_degree-1))
                DO is=1,rb(ib)%be%n_side
                  ix=seam(ib)%segment(iv)%intxys(1)
                  iy=seam(ib)%segment(iv)%intxys(2)
                  dx=ix
                  dy=iy-1+rb(ib)%be%dy(poly_degree+is)
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  tempv=rb(ib)%rz%fy
                  CALL block_crnormalize(tempv)
                  seam(ib)%segment(iv)%tang(:,is)=tempv
                  seam(ib)%segment(iv)%norm(:,is)=(/tempv(2),-tempv(1)/)
                ENDDO
              ENDIF
            ENDIF
          ENDDO vertr
c-----------------------------------------------------------------------
c         top:
c-----------------------------------------------------------------------
          vertt: DO iv=mx+my+1,2*mx+my
            IF (seam(ib)%expoint(iv)) THEN
              ix=seam(ib)%vertex(iv)%intxy(1)
              iy=seam(ib)%vertex(iv)%intxy(2)
              tempv=0
              IF (seam(ib)%expoint(ivm(iv))) THEN
                dx=ix+offset
                dy=iy
                CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                tempv=-rb(ib)%rz%fx
                CALL block_crnormalize(tempv)
              ENDIF
              IF (seam(ib)%expoint(ivp(iv))) THEN
                IF (iv<2*mx+my) THEN
                  dx=ix-offset
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fx)
                  tempv=tempv-rb(ib)%rz%fx
                ELSE
                  dx=ix
                  dy=iy-offset
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fy)
                  tempv=tempv-rb(ib)%rz%fy
                ENDIF
              ENDIF
              CALL block_crnormalize(tempv)
              seam(ib)%vertex(iv)%seam_in(1:2)=tempv
              IF (seam(ib)%expoint(ivm(iv)).AND.poly_degree>1) THEN
                ALLOCATE(seam(ib)%segment(iv)%tang(2,poly_degree-1))
                ALLOCATE(seam(ib)%segment(iv)%norm(2,poly_degree-1))
                DO is=1,rb(ib)%be%n_side
                  ix=seam(ib)%segment(iv)%intxys(1)
                  iy=seam(ib)%segment(iv)%intxys(2)
                  dx=ix-1+rb(ib)%be%dx(1+is)
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  tempv=-rb(ib)%rz%fx
                  CALL block_crnormalize(tempv)
                  seam(ib)%segment(iv)%tang(:,is)=tempv
                  seam(ib)%segment(iv)%norm(:,is)=(/tempv(2),-tempv(1)/)
                ENDDO
              ENDIF
            ENDIF
          ENDDO vertt
c-----------------------------------------------------------------------
c         left side (exceptional case is a degenerate point):
c-----------------------------------------------------------------------
          vertl: DO iv=2*mx+my+1,2*mx+2*my
            IF (seam(ib)%expoint(iv)) THEN
              ix=seam(ib)%vertex(iv)%intxy(1)
              iy=seam(ib)%vertex(iv)%intxy(2)
              tempv=0
              IF (seam(ib)%expoint(ivm(iv))) THEN
                dx=ix
                dy=iy+offset
                CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                tempv=-rb(ib)%rz%fy
                CALL block_crnormalize(tempv)
              ENDIF
              IF (seam(ib)%expoint(ivp(iv))) THEN
                IF (iv<2*mx+2*my) THEN
                  dx=ix
                  dy=iy-offset
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fy)
                  tempv=tempv-rb(ib)%rz%fy
                ELSE
                  dx=ix+offset
                  dy=iy
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  CALL block_crnormalize(rb(ib)%rz%fx)
                  tempv=tempv+rb(ib)%rz%fx
                ENDIF
              ENDIF
              tempm=SQRT(tempv(1)**2+tempv(2)**2)
              IF (tempm>0._r8) THEN
                tempv=tempv/tempm
              ELSE IF (seam(ib)%expoint(ivp(iv)).OR.
     $                 seam(ib)%expoint(ivm(iv))) THEN
                tempv(1)=0
                tempv(2)=-1
                tempm=1
              ENDIF
              seam(ib)%vertex(iv)%seam_in(1:2)=tempv
              IF (seam(ib)%expoint(ivm(iv)).AND.poly_degree>1) THEN
                ALLOCATE(seam(ib)%segment(iv)%tang(2,poly_degree-1))
                ALLOCATE(seam(ib)%segment(iv)%norm(2,poly_degree-1))
                DO is=1,rb(ib)%be%n_side
                  ix=seam(ib)%segment(iv)%intxys(1)
                  iy=seam(ib)%segment(iv)%intxys(2)
                  dx=ix
                  dy=iy-1+rb(ib)%be%dy(poly_degree+is)
                  CALL lagr_quad_eval(rb(ib)%rz,dx,dy,1_i4)
                  tempv=-rb(ib)%rz%fy
                  CALL block_crnormalize(tempv)
                  seam(ib)%segment(iv)%tang(:,is)=tempv
                  seam(ib)%segment(iv)%norm(:,is)=(/tempv(2),-tempv(1)/)
                ENDDO
              ENDIF
            ENDIF
          ENDDO vertl
        ENDIF
        DEALLOCATE(ivp,ivm)
      ENDDO block
c-----------------------------------------------------------------------
c     seam the vertex tangent vectors to get an average where there
c     may be ambiguity.
c-----------------------------------------------------------------------
      CALL edge_network(2_i4,0_i4,0_i4,.false.)
c-----------------------------------------------------------------------
c     copy the tangent vector from networking to seam, renormalize
c     for block averaging, and find the unit normal.  this is
c     only applied those pts that are external vertices.
c-----------------------------------------------------------------------
      DO ibe=1,SIZE(exblock_list)
        ib=exblock_list(ibe)
        DO iv=1,seam(ib)%nvert
          ix=seam(ib)%vertex(iv)%intxy(1)
          iy=seam(ib)%vertex(iv)%intxy(2)
          IF (seam(ib)%expoint(iv)) THEN
            ALLOCATE(seam(ib)%vertex(iv)%tang(2))
            ALLOCATE(seam(ib)%vertex(iv)%norm(2))
            CALL block_crnormalize(seam(ib)%vertex(iv)%seam_out(1:2))
            seam(ib)%vertex(iv)%tang(:)=
     $        seam(ib)%vertex(iv)%seam_out(1:2)
            seam(ib)%vertex(iv)%norm(1)= seam(ib)%vertex(iv)%tang(2)
            seam(ib)%vertex(iv)%norm(2)=-seam(ib)%vertex(iv)%tang(1)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN

      CONTAINS

c-----------------------------------------------------------------------
c       this internal subroutine normalizes two-vectors.
c-----------------------------------------------------------------------

        SUBROUTINE block_crnormalize(vec)

        REAL(r8), DIMENSION(2), INTENT(INOUT) :: vec
        REAL(r8) :: vmag

        vmag=SQRT(vec(1)**2+vec(2)**2)
        IF (vmag>0._r8) THEN
          vec=vec/vmag
        ENDIF

        END SUBROUTINE block_crnormalize

      END SUBROUTINE block_create_tang
c-----------------------------------------------------------------------
c     subprogram 2. block_dealloc_tang
c     deallocate the norm and tang arrays for each block.
c-----------------------------------------------------------------------
      SUBROUTINE block_dealloc_tang
      USE local
      USE seam_storage_mod
      USE fields
      IMPLICIT NONE

      INTEGER(i4) :: ib,iv
c-----------------------------------------------------------------------
c     deallocate the norm and tang arrays in places where the storage is
c     associated.
c-----------------------------------------------------------------------
      DO ib=1,nbl
        DO iv=1,seam(ib)%nvert
          IF (ASSOCIATED(seam(ib)%vertex(iv)%norm))
     $      DEALLOCATE(seam(ib)%vertex(iv)%norm)
          IF (ASSOCIATED(seam(ib)%vertex(iv)%tang))
     $      DEALLOCATE(seam(ib)%vertex(iv)%tang)
          IF (ASSOCIATED(seam(ib)%segment(iv)%norm))
     $      DEALLOCATE(seam(ib)%segment(iv)%norm)
          IF (ASSOCIATED(seam(ib)%segment(iv)%tang))
     $      DEALLOCATE(seam(ib)%segment(iv)%tang)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE block_dealloc_tang
