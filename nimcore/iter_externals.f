c-----------------------------------------------------------------------
c     file iter_externals.f
c     this collection of subroutines allows access to block-specific
c     information without using the fields module in the calling
c     routine.  these subprograms had been in iter_cg.f and have been
c     separated for modularity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. per_block.
c     2. border_weights.
c     3. deg_block.
c     4. concave_corner.
c     5. direct_check.
c-----------------------------------------------------------------------
c     subprogram 1. per_block.
c function to determine whether a single block is self-periodic
c return .TRUE. if ALL edge vertices 0:mx have a seam connection to self
c return .FALSE. if ANY do not.  this is used in both real and complex
c versions of the cg solver and must be outside the iter_cg module.
c-----------------------------------------------------------------------
      LOGICAL FUNCTION per_block(ibl,mx)
      USE local
      USE pardata
      USE seam_storage_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: ibl,mx

      INTEGER(i4) :: ivert,image,nv,ib_image
c-----------------------------------------------------------------------

c for each vertex, set per_block to FALSE
c reset to TRUE if any seam image points to this block
c when one FALSE is found, return

      nv=seam(ibl)%nvert
      per_block = .FALSE.
      do image = 1,seam(ibl)%vertex(nv)%nimage
         ib_image=seam(ibl)%vertex(nv)%ptr2(1,image)
         if (block2proc(ib_image) == node .AND. 
     $       global2local(ib_image) == ibl) then
            per_block = .TRUE.
            exit
         endif
      enddo
      if (.NOT. per_block) return

      do ivert = 1,mx
         per_block = .FALSE.
         do image = 1,seam(ibl)%vertex(ivert)%nimage
           ib_image=seam(ibl)%vertex(ivert)%ptr2(1,image)
           if (block2proc(ib_image) == node .AND. 
     $         global2local(ib_image) == ibl) then
               per_block = .TRUE.
               exit
            endif
         enddo
         if (.NOT. per_block) return
      enddo

c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION per_block

c-----------------------------------------------------------------------
c     subprogram 2. border_weights.
c     initializes the weight factors along block borders that are
c     used during the iterative solves.
c-----------------------------------------------------------------------
      SUBROUTINE border_weights(precon)
      USE local
      USE fields
      USE seam_storage_mod
      USE edge
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: iv,ibl,nv
      LOGICAL, EXTERNAL :: per_block,direct_check
c-----------------------------------------------------------------------
c     also create block-boundary node weights for use by the iterative
c     solver.  bottoms of periodic blocks are zeroed, so tops
c     get full weight, and for degenerate blocks, only the iy=1
c     position is weighted.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        nv=seam(ibl)%nvert
        DO iv=1,nv
          seam(ibl)%vertex(iv)%seam_in(1:2)=1._r8
          IF (per_block(ibl,rb(ibl)%mx).AND.(iv==nv.OR.iv<=rb(ibl)%mx))
     $      seam(ibl)%vertex(iv)%seam_in(1:2)=0._r8
          IF (rb(ibl)%degenerate.AND.iv>2*rb(ibl)%mx+rb(ibl)%my)
     $      seam(ibl)%vertex(iv)%seam_in(1:2)=0._r8
          seam(ibl)%segment(iv)%seam_in(1:2)=1._r8
          IF (per_block(ibl,rb(ibl)%mx).AND.(iv<=rb(ibl)%mx))
     $      seam(ibl)%segment(iv)%seam_in(1:2)=0._r8
          IF (rb(ibl)%degenerate.AND.iv>2*rb(ibl)%mx+rb(ibl)%my)
     $      seam(ibl)%segment(iv)%seam_in(1:2)=0._r8
        ENDDO
      ENDDO
      DO ibl=nrbl+1,nbl
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_in(1)=1._r8
          seam(ibl)%vertex(iv)%seam_in(2)=0._r8
          seam(ibl)%segment(iv)%seam_in(1:2)=0._r8
        ENDDO
      ENDDO
      CALL edge_network(2_i4,0_i4,1_i4,.false.)
c-----------------------------------------------------------------------
c     the averaging factors used in preconditioning are placed in
c     a separate variable, ave_factor_pre.
c     if global preconditioning is used in rblocks, tblock
c     contributions along borders with rblocks are set to 0, and
c     only tblock-tblock factors use sqrt.
c
c     the segment average gives 1/2 between rblocks, 1 on the rblock
c     side of block borders with tblocks and 0 on the tblock side.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        nv=seam(ibl)%nvert
        DO iv=1,nv
          seam(ibl)%vertex(iv)%ave_factor=
     $      1._r8/seam(ibl)%vertex(iv)%seam_out(1)
          seam(ibl)%vertex(iv)%ave_factor_pre=
     $      SQRT(seam(ibl)%vertex(iv)%ave_factor)
          IF (ibl<=nrbl) THEN
            IF (precon=='gl_diaga'.OR.direct_check(precon))
     $        seam(ibl)%vertex(iv)%ave_factor_pre=
     $          1._r8/seam(ibl)%vertex(iv)%seam_out(2)
            IF (per_block(ibl,rb(ibl)%mx).AND.
     $          (iv==nv.OR.iv<=rb(ibl)%mx)) THEN
              seam(ibl)%vertex(iv)%ave_factor=0._r8
              seam(ibl)%vertex(iv)%ave_factor_pre=0._r8
            ENDIF
            IF (rb(ibl)%degenerate.AND.
     $          iv>=2*rb(ibl)%mx+rb(ibl)%my.AND.iv/=nv-1) THEN
              seam(ibl)%vertex(iv)%ave_factor=0._r8
              seam(ibl)%vertex(iv)%ave_factor_pre=0._r8
            ENDIF
          ELSE IF ((precon=='gl_diaga'.OR.direct_check(precon)).AND.
     $             seam(ibl)%vertex(iv)%seam_out(2)>0) THEN
            seam(ibl)%vertex(iv)%ave_factor_pre=0._r8
          ENDIF

          IF (ibl<=nrbl) THEN
            seam(ibl)%segment(iv)%ave_factor=
     $        1/seam(ibl)%segment(iv)%seam_out(1)
            IF (rb(ibl)%degenerate.AND.iv>2*rb(ibl)%mx+rb(ibl)%my)
     $        seam(ibl)%segment(iv)%ave_factor=0._r8
            IF (per_block(ibl,rb(ibl)%mx).AND.iv<=rb(ibl)%mx)
     $        seam(ibl)%segment(iv)%ave_factor=0._r8
          ELSE
            seam(ibl)%segment(iv)%ave_factor=0._r8
          ENDIF
          IF (precon=='gl_diaga'.OR.direct_check(precon)) THEN
            seam(ibl)%segment(iv)%ave_factor_pre=
     $        seam(ibl)%segment(iv)%ave_factor
          ELSE
            seam(ibl)%segment(iv)%ave_factor_pre=
     $        SQRT(seam(ibl)%segment(iv)%ave_factor)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE border_weights

c-----------------------------------------------------------------------
c     subprogram 3. deg_block.
c     determine whether a block has a degenerate side.  return T if
c     degenerate, F if not.
c-----------------------------------------------------------------------
      LOGICAL FUNCTION deg_block(ibl)
      USE local
      USE fields
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: ibl

      deg_block=.false.
      IF (ibl<=nrbl) THEN
        IF (rb(ibl)%degenerate) deg_block=.true.
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION deg_block
c-----------------------------------------------------------------------
c     subprogram 4. concave_corner.
c     test if the offset indices passed into this function define
c     a matrix element across a block's corner.
c
c     lx = x-direction column offset (from row index)
c     ly = y-direction column offset (from row index)
c     ix0, iy0 = block starting indices for the row basis type.
c     jx0, jy0 = block starting indices for the column basis type.
c-----------------------------------------------------------------------
      LOGICAL FUNCTION concave_corner(lx,ly,ix0,iy0,jx0,jy0)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: lx,ly,ix0,iy0,jx0,jy0

      concave_corner=(lx/=0).AND.(ly/=0).OR.
     $               (lx/=0).AND.(iy0/=jy0).OR.
     $               (ly/=0).AND.(ix0/=jx0).OR.
     $               (lx==0).AND.(ly==0).AND.(ix0/=jx0).AND.
     $               (iy0/=jy0)
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION concave_corner
c-----------------------------------------------------------------------
c     subprogram 5. direct_check.
c     test if the preconditioner option is a direct solve over the
c     fe mesh.
c-----------------------------------------------------------------------
      LOGICAL FUNCTION direct_check(precon)
      USE local
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: precon

      IF (precon=='lapack'.OR.precon=='seq_slu'.OR.
     $    precon(1:5)=='slu_d') THEN
        direct_check=.true.
      ELSE
        direct_check=.false.
      ENDIF
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION direct_check
