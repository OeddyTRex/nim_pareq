c-----------------------------------------------------------------------
c     file regularity.f:  a module containing routines that apply
c     regularity conditions to an R=0 domain border.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declaration.
c-----------------------------------------------------------------------
      MODULE regularity
      USE local
      USE matrix_type_mod
      USE vector_type_mod
      IMPLICIT NONE

      INTEGER(i4), DIMENSION(:), POINTER :: r0block_list
      INTEGER(i4) :: ndr_pass
      LOGICAL :: any_r0blocks

      CONTAINS
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. regular_vec.
c     2. regular_ave.
c     3. regular_zero_phi.
c     4. regular_pre_feop.
c     5. regular_op.
c     6. regular_comp_op.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 1. regular_vec.
c     apply regularity conditions to a vector. 
c-----------------------------------------------------------------------
      SUBROUTINE regular_vec(vec,blseam,flag,nq,nmodes,nindex)
      USE edge_type_mod

      INTEGER(i4), INTENT(IN) :: nq,nmodes
      INTEGER(i4), DIMENSION(nmodes), INTENT(IN) :: nindex
      TYPE(cvector_type), INTENT(INOUT) :: vec
      TYPE(edge_type), INTENT(IN) :: blseam
      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: iv,ivp,ix,iy,iqty,imode,imn1,iside,ivec,nvec
      REAL(r8), DIMENSION(nq,nmodes) :: mult
c-----------------------------------------------------------------------
c     apply regularity conditions to the different Fourier components
c     of a vector.  for a scalar quantity, n>0 are set to 0.  for a
c     vector, thd following components are set to 0
c       n=0:  r and phi
c       n=1:  z
c       n>1:  r, z, and phi
c
c     create an array of 1s and 0s to zero out the appropriate
c     components.  the geometry is always toroidal at this point.
c-----------------------------------------------------------------------
      mult=0
      imn1=0
      SELECT CASE(nq)
      CASE(1,2)              !   scalars
        DO imode=1,nmodes
          IF (nindex(imode)==0) THEN
            mult(:,imode)=1
            EXIT
          ENDIF
        ENDDO
      CASE(3,6,9,12,15,18)   !   3-vectors
        nvec=nq/3
        DO ivec=0,nvec-1
          DO imode=1,nmodes
            IF (nindex(imode)==0) THEN
              mult(3*ivec+2,imode)=1
              CYCLE
            ENDIF
            IF (nindex(imode)==1) THEN
              mult(3*ivec+1,imode)=1
              IF (flag=='cyl_vec') imn1=imode
              IF (flag=='init_cyl_vec') mult(3*ivec+3,imode)=1
            ENDIF
          ENDDO
        ENDDO
      CASE DEFAULT
        CALL nim_stop("Regular_vec: inconsistent # of components.")
      END SELECT
c-----------------------------------------------------------------------
c     loop over block border elements and apply mult to the vertices
c     at R=0.  also combine rhs for n=1 r and phi comps.
c-----------------------------------------------------------------------
      vert: DO iv=1,blseam%nvert
        IF (.NOT.blseam%r0point(iv)) CYCLE
        ix=blseam%vertex(iv)%intxy(1)
        iy=blseam%vertex(iv)%intxy(2)
        IF (imn1/=0) THEN
          DO ivec=0,nvec-1
            vec%arr(3*ivec+1,ix,iy,imn1)=vec%arr(3*ivec+1,ix,iy,imn1)
     $                            -(0,1)*vec%arr(3*ivec+3,ix,iy,imn1)
          ENDDO
        ENDIF
        DO imode=1,nmodes
          vec%arr(1:nq,ix,iy,imode)=
     $      vec%arr(1:nq,ix,iy,imode)*mult(:,imode)
        ENDDO
        ivp=iv-1
        IF (ivp==0) ivp=blseam%nvert
        IF (.NOT.blseam%r0point(ivp)) CYCLE
        ix=blseam%segment(iv)%intxys(1)
        iy=blseam%segment(iv)%intxys(2)
        IF (blseam%segment(iv)%h_side.AND.ASSOCIATED(vec%arrh)) THEN
          IF (imn1/=0) THEN
            DO ivec=0,nvec-1
              vec%arrh(3*ivec+1,:,ix,iy,imn1)=
     $                 vec%arrh(3*ivec+1,:,ix,iy,imn1)
     $          -(0,1)*vec%arrh(3*ivec+3,:,ix,iy,imn1)
            ENDDO
          ENDIF
          DO imode=1,nmodes
            DO iside=1,SIZE(vec%arrh,2)
              vec%arrh(1:nq,iside,ix,iy,imode)=
     $          vec%arrh(1:nq,iside,ix,iy,imode)*mult(:,imode)
            ENDDO
          ENDDO
        ENDIF
        IF (.NOT.blseam%segment(iv)%h_side.AND.
     $      ASSOCIATED(vec%arrv)) THEN
          IF (imn1/=0) THEN
            DO ivec=0,nvec-1
              vec%arrv(3*ivec+1,:,ix,iy,imn1)=
     $                 vec%arrv(3*ivec+1,:,ix,iy,imn1)
     $          -(0,1)*vec%arrv(3*ivec+3,:,ix,iy,imn1)
            ENDDO
          ENDIF
          DO imode=1,nmodes
            DO iside=1,SIZE(vec%arrv,2)
              vec%arrv(1:nq,iside,ix,iy,imode)=
     $          vec%arrv(1:nq,iside,ix,iy,imode)*mult(:,imode)
            ENDDO
          ENDDO
        ENDIF
      ENDDO vert
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_vec
c-----------------------------------------------------------------------
c     subprogram 2. regular_ave.
c     set the real phi component of n=1 to minus the imaginary radial
c     component, and set the imaginary phi component to the real
c     radial component.  this should only be called for vectors with
c     cylindrical components (not covariant or contravariant).
c     it is assumed that the matrix equation has already solved for
c     the appropriate averages, placing the result in the radial
c     index locations.
c-----------------------------------------------------------------------
      SUBROUTINE regular_ave(vec,nqty,nmodes,nindex)
      USE seam_storage_mod
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: nqty,nmodes
      INTEGER(i4), DIMENSION(nmodes), INTENT(IN) :: nindex
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vec

      INTEGER(i4) :: ibl,ib0,iv,ivp,ix,iy,imode,nvec,ivec
c-----------------------------------------------------------------------
c     return if no blocks touch R=0.
c-----------------------------------------------------------------------
      IF (SIZE(r0block_list)<=0) RETURN
c-----------------------------------------------------------------------
c     check the number of quantities per Fourier component.
c-----------------------------------------------------------------------
      IF (MODULO(nqty,3_i4)/=0) CALL nim_stop
     $  ("Regular_ave: passed array is not a vector.")
      nvec=nqty/3
c-----------------------------------------------------------------------
c     find the index of the n=1 Fourier component.
c-----------------------------------------------------------------------
      imode=1
      DO
        IF (imode>nmodes) RETURN
        IF (nindex(imode)==1) EXIT
        imode=imode+1
      ENDDO
c-----------------------------------------------------------------------
c     for the n=1 Fourier component, set the phi component
c     to i times the radial component at all R=0 points.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(r0block_list)
        ib0=r0block_list(ibl)
        DO iv=1,seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(iv)) CYCLE
          ix=seam(ib0)%vertex(iv)%intxy(1)
          iy=seam(ib0)%vertex(iv)%intxy(2)
          DO ivec=0,nvec-1
            vec(ib0)%arr(3*ivec+3,ix,iy,imode)=
     $        (0,1)*vec(ib0)%arr(3*ivec+1,ix,iy,imode)
          ENDDO
          ivp=iv-1
          IF (ivp==0) ivp=seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(ivp)) CYCLE
          ix=seam(ib0)%segment(iv)%intxys(1)
          iy=seam(ib0)%segment(iv)%intxys(2)
          IF (seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrh)) THEN
            DO ivec=0,nvec-1
              vec(ib0)%arrh(3*ivec+3,:,ix,iy,imode)=
     $          (0,1)*vec(ib0)%arrh(3*ivec+1,:,ix,iy,imode)
            ENDDO
          ENDIF
          IF (.NOT.seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrv)) THEN
            DO ivec=0,nvec-1
              vec(ib0)%arrv(3*ivec+3,:,ix,iy,imode)=
     $          (0,1)*vec(ib0)%arrv(3*ivec+1,:,ix,iy,imode)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_ave
c-----------------------------------------------------------------------
c     subprogram 3. regular_zero_phi.
c     zero-out the n=1 phi component of vectors at nodes along the R=0
c     edge of the computational mesh.  this is used after matrix-vector
c     multiplication in 3d cg solves, where that storage temporarily
c     represents the half difference between A_r and -i*A_phi.
c-----------------------------------------------------------------------
      SUBROUTINE regular_zero_phi(vec,nqty,nmodes,nindex)
      USE seam_storage_mod
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: nqty,nmodes
      INTEGER(i4), DIMENSION(nmodes), INTENT(IN) :: nindex
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vec

      INTEGER(i4) :: ibl,ib0,iv,ivp,ix,iy,imode
c-----------------------------------------------------------------------
c     return if no blocks touch R=0.
c-----------------------------------------------------------------------
      IF (SIZE(r0block_list)<=0) RETURN
c-----------------------------------------------------------------------
c     check the number of quantities per Fourier component.
c-----------------------------------------------------------------------
      IF (MODULO(nqty,3_i4)/=0) CALL nim_stop
     $  ("Regular_zero_phi: passed array is not a vector.")
c-----------------------------------------------------------------------
c     find the index of the n=1 Fourier component.
c-----------------------------------------------------------------------
      imode=1
      DO
        IF (imode>nmodes) RETURN
        IF (nindex(imode)==1) EXIT
        imode=imode+1
      ENDDO
c-----------------------------------------------------------------------
c     for the n=1 Fourier component, set the phi component
c     to 0 at all R=0 points.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(r0block_list)
        ib0=r0block_list(ibl)
        DO iv=1,seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(iv)) CYCLE
          ix=seam(ib0)%vertex(iv)%intxy(1)
          iy=seam(ib0)%vertex(iv)%intxy(2)
          vec(ib0)%arr(3:nqty:3,ix,iy,imode)=0
          ivp=iv-1
          IF (ivp==0) ivp=seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(ivp)) CYCLE
          ix=seam(ib0)%segment(iv)%intxys(1)
          iy=seam(ib0)%segment(iv)%intxys(2)
          IF (seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrh)) THEN
            vec(ib0)%arrh(3:nqty:3,:,ix,iy,imode)=0
          ENDIF
          IF (.NOT.seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrv)) THEN
            vec(ib0)%arrv(3:nqty:3,:,ix,iy,imode)=0
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_zero_phi
c-----------------------------------------------------------------------
c     subprogram 4. regular_pre_feop.
c     this routine is a variant of regular_vec, as needed before
c     using finite-element operations for matrix-free Kylov-space
c     products.  unlike regular_vec values are saved in seam_csave when
c     setting degrees of freedom to zero in the vector.  also, for
c     n=1 vectors, the phi-component of the vector is set to i times
c     the r-component, as in regular_ave.
c-----------------------------------------------------------------------
      SUBROUTINE regular_pre_feop(vec,flag,nq,nmodes,nindex)
      USE seam_storage_mod
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: nq,nmodes
      INTEGER(i4), DIMENSION(nmodes), INTENT(IN) :: nindex
      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vec
      CHARACTER(*), INTENT(IN) :: flag

      INTEGER(i4) :: iv,ivp,ix,iy,iqty,imode,imn1,iside,ivec,nvec,
     $               ibl,ib0,iq
      REAL(r8), DIMENSION(nq,nmodes) :: mult,esave
c-----------------------------------------------------------------------
c     apply regularity conditions to the different Fourier components
c     of a vector.  for a scalar quantity, n>0 are set to 0.  for a
c     vector, thd following components are set to 0
c       n=0:  r and phi
c       n=1:  z
c       n>1:  r, z, and phi
c
c     create an array of 1s and 0s to zero out the appropriate
c     components.  the geometry is always toroidal at this point.
c-----------------------------------------------------------------------
      mult=0
      imn1=0
      SELECT CASE(nq)
      CASE(1,2)              !   scalars
        DO imode=1,nmodes
          IF (nindex(imode)==0) THEN
            mult(:,imode)=1
            EXIT
          ENDIF
        ENDDO
      CASE(3,6,9,12,15,18)   !   3-vectors
        nvec=nq/3
        DO ivec=0,nvec-1
          DO imode=1,nmodes
            IF (nindex(imode)==0) THEN
              mult(3*ivec+2,imode)=1
              CYCLE
            ENDIF
            IF (nindex(imode)==1) THEN
              mult(3*ivec+1,imode)=1
              IF (flag=='cyl_vec') imn1=imode
              IF (flag=='init_cyl_vec') mult(3*ivec+3,imode)=1
            ENDIF
          ENDDO
        ENDDO
      CASE DEFAULT
        CALL nim_stop("Regular_vec: inconsistent # of components.")
      END SELECT
      esave=1_i4-mult
c-----------------------------------------------------------------------
c     loop over block border elements, save the values that mult
c     eliminates, then apply mult to the vertices
c     at R=0.  also set the phi-component to i times the r-component.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(r0block_list)
        ib0=r0block_list(ibl)
        vert: DO iv=1,seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(iv)) CYCLE
          ix=seam(ib0)%vertex(iv)%intxy(1)
          iy=seam(ib0)%vertex(iv)%intxy(2)
          iq=1
          DO imode=1,nmodes
            seam(ib0)%vertex(iv)%seam_csave(iq:iq+nq-1)=
     $        seam(ib0)%vertex(iv)%seam_csave(iq:iq+nq-1)+
     $        vec(ib0)%arr(1:nq,ix,iy,imode)*esave(:,imode)
            vec(ib0)%arr(1:nq,ix,iy,imode)=
     $        vec(ib0)%arr(1:nq,ix,iy,imode)*mult(:,imode)
            iq=iq+nq
          ENDDO
          IF (imn1/=0) THEN
            DO ivec=0,nvec-1
              vec(ib0)%arr(3*ivec+3,ix,iy,imn1)=
     $          (0,1)*vec(ib0)%arr(3*ivec+1,ix,iy,imn1)
            ENDDO
          ENDIF
          ivp=iv-1
          IF (ivp==0) ivp=seam(ib0)%nvert
          IF (.NOT.seam(ib0)%r0point(ivp)) CYCLE
          ix=seam(ib0)%segment(iv)%intxys(1)
          iy=seam(ib0)%segment(iv)%intxys(2)
          IF (seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrh)) THEN
            iq=1
            DO imode=1,nmodes
              DO iside=1,SIZE(vec(ib0)%arrh,2)
                seam(ib0)%segment(iv)%seam_csave(iq:iq+nq-1)=
     $            seam(ib0)%segment(iv)%seam_csave(iq:iq+nq-1)+
     $            vec(ib0)%arrh(1:nq,iside,ix,iy,imode)*esave(:,imode)
                vec(ib0)%arrh(1:nq,iside,ix,iy,imode)=
     $            vec(ib0)%arrh(1:nq,iside,ix,iy,imode)*mult(:,imode)
                iq=iq+nq
              ENDDO
            ENDDO
            IF (imn1/=0) THEN
              DO ivec=0,nvec-1
                vec(ib0)%arrh(3*ivec+3,:,ix,iy,imn1)=
     $            (0,1)*vec(ib0)%arrh(3*ivec+1,:,ix,iy,imn1)
              ENDDO
            ENDIF
          ENDIF
          IF (.NOT.seam(ib0)%segment(iv)%h_side.AND.
     $        ASSOCIATED(vec(ib0)%arrv)) THEN
            iq=1
            DO imode=1,nmodes
              DO iside=1,SIZE(vec(ib0)%arrv,2)
                seam(ib0)%segment(iv)%seam_csave(iq:iq+nq-1)=
     $            seam(ib0)%segment(iv)%seam_csave(iq:iq+nq-1)+
     $            vec(ib0)%arrv(1:nq,iside,ix,iy,imode)*esave(:,imode)
                vec(ib0)%arrv(1:nq,iside,ix,iy,imode)=
     $            vec(ib0)%arrv(1:nq,iside,ix,iy,imode)*mult(:,imode)
                iq=iq+nq
              ENDDO
            ENDDO
            IF (imn1/=0) THEN
              DO ivec=0,nvec-1
                vec(ib0)%arrv(3*ivec+3,:,ix,iy,imn1)=
     $            (0,1)*vec(ib0)%arrv(3*ivec+1,:,ix,iy,imn1)
              ENDDO
            ENDIF
          ENDIF
        ENDDO vert
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_pre_feop
c-----------------------------------------------------------------------
c     subprogram 5. regular_op.
c     apply regularity conditions to R=0 points of an operator.
c     the action on the matrix is (I-uu).M.(I-uu)+uu, where u is the 
c     unit vector for the components that need elimination.  the
c     vcomp list of component descriptions in the matrix structure
c     and the fcomp Fourier component index are used to decide
c     where to apply this action.
c
c     if dscale is provided, it is used to scale the added uu entry, 
c     instead of determining the scaling on the fly.
c-----------------------------------------------------------------------
      SUBROUTINE regular_op(mat,dscale)
      USE seam_storage_mod

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      REAL(r8), INTENT(IN), OPTIONAL :: dscale

      REAL(r8), DIMENSION(mat%nqty,mat%nqty) :: mult,multj
      INTEGER(i4) :: iv,nv,ip,ipn,ibl,ibe,ibv,jxmin,jxmax,jymin,jymax,
     $               jx,jy,ix,iy,ijx,ijy,mx,my,iq,jq,clim,ixn,
     $               ix0,iy0,jx0,jy0,imat,jmat,ivp,nvec,ivec,iv0,iv1,
     $               iqt,jqt,iqp,jqp,itype,jtype,iq0,iq1,jq0,jq1,jcomp
      REAL(r8), DIMENSION(mat%nqty) :: mult_v,mult_j
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      REAL(r8) :: diagp
      TYPE(matrix_element_type3), DIMENSION(:), POINTER :: tmat
c-----------------------------------------------------------------------
c     apply regularity conditions to the different vector components
c     of a matrix.  the elements at R=0 are uncoupled for n>0, if the
c     matrix acts on a scalar field.  the elements at R=0 of a matrix 
c     for a vector field are uncoupled for:
c       n=0:  r and phi
c       n=1:  z
c       n>1:  r, z, and phi
c     a covariant phi component (a_phi*r) is uncoupled for all n.  these
c     components are indicated with 'c'.
c
c-PRE
c     this version works for multiple 3-vectors in a system, but that
c     generalization is not extended to tblocks.
c-----------------------------------------------------------------------
      clim=mat%nqty
      nvec=clim/3
c-----------------------------------------------------------------------
c     create an array of 1s and 0s to zero out the appropriate
c     components.  the geometry is always toroidal at this point.
c
c     the second array, mult_j, is used for the columns.
c     [preconditioning may require matrix elements with different
c     Fourier indices for rows and columns.]
c-----------------------------------------------------------------------
      mult_v=0
      SELECT CASE(mat%fcomp)
      CASE(0)
        WHERE(mat%vcomp=='s') mult_v=1
        WHERE(mat%vcomp=='z') mult_v=1
      CASE(1)
        WHERE(mat%vcomp=='r') mult_v=1
      END SELECT
      mult_j=0
      jcomp=mat%fcomp+mat%foff
      SELECT CASE(ABS(jcomp))
      CASE(0)
        WHERE(mat%vcomp=='s') mult_j=1
        WHERE(mat%vcomp=='z') mult_j=1
      CASE(1)
        WHERE(mat%vcomp=='r') mult_j=1
      END SELECT
      DO jq=1,clim
        mult(jq,:)=mult_v
        multj(:,jq)=mult_j
      ENDDO
c-----------------------------------------------------------------------
c     loop over the borders of blocks touching R=0, and decouple the
c     appropriate matrix elements.
c-----------------------------------------------------------------------
      block: DO ibe=1,SIZE(r0block_list)
        ibl=r0block_list(ibe)
        nv=seam(ibl)%nvert
c-----------------------------------------------------------------------
c       loop over the block boundary.
c-----------------------------------------------------------------------
        vert: DO iv=1,nv
          ivp=iv-1
          IF (ivp==0) ivp=nv
c-----------------------------------------------------------------------
c         begin rblock-specific coding.
c-----------------------------------------------------------------------
          block_type: IF (ibl<=SIZE(mat%rbl_mat)) THEN
            mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
            my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
c-----------------------------------------------------------------------
c           loop over couplings from a given node at R=0 for n=1
c           vectors only.
c-----------------------------------------------------------------------
            IF (ABS(jcomp)==1.AND.MODULO(clim,3_i4)==0) THEN
              DO itype=1,mat%rbl_mat(ibl)%nbtype
                ix0=mat%rbl_mat(ibl)%ix0(itype)
                iy0=mat%rbl_mat(ibl)%iy0(itype)
                DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                  rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                  jx0=mat%rbl_mat(ibl)%ix0(jtype)
                  jy0=mat%rbl_mat(ibl)%iy0(jtype)
                  DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                    IF (jtype==1) THEN
                      IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                      ix=seam(ibl)%vertex(iv)%intxy(1)
                      iy=seam(ibl)%vertex(iv)%intxy(2)
                    ELSE
                      IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                          seam(ibl)%r0point(ivp))) CYCLE
                      IF (seam(ibl)%segment(iv)%h_side) THEN
                        IF (jtype>2) CYCLE
                      ELSE
                        IF (jtype<3) CYCLE
                      ENDIF
                      ix=seam(ibl)%segment(iv)%intxys(1)
                      iy=seam(ibl)%segment(iv)%intxys(2)
                    ENDIF
                    jxmin=MAX(ix0-1,ix0-ix)
                    jxmax=MIN(1-jx0,mx -ix)
                    jymin=MAX(iy0-1,iy0-iy)
                    jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                   combine r and phi equations to find averages in the
c                   matrix that enforce Vec_phi=i*Vec_r for n=1.  for
c                   preconditioning, we may also need Vec_phi=-i*Vec_r
c                   for n=-1.
c
c                   mathematically, the steps are 1) change variables 
c                   to X1=(Vec_r+i*Vec_phi)/2 and X2=(Vec_r-i*Vec_phi)/2
c                   for each vertex at R=0, 2) set X1=0 (remove
c                   couplings to and from X1) giving an overdetermined
c                   system of equations, 3) add -i*(the X2-equation) to
c                   the X1 equation.
c
c                   for the n=-1 preconditioning matrix case, we will
c                   set X2 to zero and keep X1 in the r-comp location.
c-----------------------------------------------------------------------
                    DO ivec=0,nvec-1
                      jq0=(jmat-1)*clim+3*ivec+1
                      jq1=jq0+2
                      DO jy=jymin,jymax
                        ijy=iy+jy
                        DO jx=jxmin,jxmax
                          ijx=ix+jx
                          rmat(jq0,-jx,-jy,:,ijx,ijy)=
     $                       rmat(jq0,-jx,-jy,:,ijx,ijy)
     $                      -rmat(jq1,-jx,-jy,:,ijx,ijy)*jcomp
                        ENDDO
                      ENDDO 
                    ENDDO 
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
c-----------------------------------------------------------------------
c           loop over couplings to a given node at R=0.
c-----------------------------------------------------------------------
            DO itype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,mat%rbl_mat(ibl)%nbtype
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO imat=1,mat%rbl_mat(ibl)%nb_type(itype)
                  IF (itype==1) THEN
                    IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                  ELSE
                    IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                        seam(ibl)%r0point(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (itype>2) CYCLE
                    ELSE
                      IF (itype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                  ENDIF
                  jxmin=MAX(jx0-1,jx0-ix)
                  jxmax=MIN(1-ix0,mx -ix)
                  jymin=MAX(jy0-1,jy0-iy)
                  jymax=MIN(1-iy0,my -iy)
c-----------------------------------------------------------------------
c                 combine r and phi equations to find averages in the
c                 matrix that enforce Vec_phi=i*Vec_r for n=1.  [note
c                 rows always represent n>=0.]
c
c                 then, eliminate couplings to this node.
c-----------------------------------------------------------------------
                  iq0=(imat-1)*clim+1
                  iq1=imat*clim
                  DO jy=jymin,jymax
                    DO jx=jxmin,jxmax
                      IF (mat%fcomp==1.AND.MODULO(clim,3_i4)==0) THEN
                        DO ivec=0,nvec-1
                          iv0=(imat-1)*clim+3*ivec+1
                          iv1=iv0+2
                          rmat(:,jx,jy,iv0,ix,iy)=
     $                      rmat(:,jx,jy,iv0,ix,iy)
     $                     -rmat(:,jx,jy,iv1,ix,iy)
                        ENDDO
                      ENDIF
                      DO jq=1,mat%rbl_mat(ibl)%nq_type(jtype)
                        rmat(jq,jx,jy,iq0:iq1,ix,iy)=
     $                    rmat(jq,jx,jy,iq0:iq1,ix,iy)*mult_v
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           again loop over couplings from a given node at R=0.
c-----------------------------------------------------------------------
            DO itype=1,mat%rbl_mat(ibl)%nbtype
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                  IF (jtype==1) THEN
                    IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                  ELSE
                    IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                        seam(ibl)%r0point(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (jtype>2) CYCLE
                    ELSE
                      IF (jtype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                  ENDIF
                  jxmin=MAX(ix0-1,ix0-ix)
                  jxmax=MIN(1-jx0,mx -ix)
                  jymin=MAX(iy0-1,iy0-iy)
                  jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                 eliminate couplings from the node on R=0.
c-----------------------------------------------------------------------
                  jq0=(jmat-1)*clim+1
                  jq1=jmat*clim
                  DO jy=jymin,jymax
                    ijy=iy+jy
                    DO jx=jxmin,jxmax
                      ijx=ix+jx
                      DO iq=1,mat%rbl_mat(ibl)%nq_type(itype)
                        rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)=
     $                    rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)*mult_j
                      ENDDO
                    ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c                 take care of uu and zz.
c-----------------------------------------------------------------------
                  IF (jtype==itype.AND.mat%foff==0) THEN
                    IF (PRESENT(dscale)) THEN
                      diagp=dscale
                    ELSE
                      diagp=0 
                      DO jq=0,clim-1
                        diagp=MAX(diagp,rmat(jq+jq0,0,0,jq+jq0,ix,iy))
                      ENDDO
                      IF (diagp==0) diagp=1
                    ENDIF
                    DO jq=0,clim-1
                      rmat(jq+jq0,0,0,jq+jq0,ix,iy)=
     $                  rmat(jq+jq0,0,0,jq+jq0,ix,iy)
     $                 +(1-mult_v(jq+1))*diagp
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c         begin tblock-specific coding.
c-----------------------------------------------------------------------
          ELSE block_type
            tmat=>mat%tbl_mat(ibl)%lmat
c-----------------------------------------------------------------------
c           combine r and phi equations to find averages in the matrix
c           that enforce Vec_phi=i*Vec_r for n=1.  if clim is 3,
c           components are either (r_r,z_r,-phi_i) or (r_i,z_i,phi_r).
c-----------------------------------------------------------------------
            IF (mat%fcomp==1.AND.clim==3) THEN
              DO ip=0,SIZE(tmat(ix)%from_vert)-1
                tmat(ix)%element(:,1,ip)=tmat(ix)%element(:,1,ip)
     $                                  -tmat(ix)%element(:,3,ip)
              ENDDO
            ENDIF
            IF (ABS(jcomp)==1.AND.clim==3) THEN
              DO ip=0,SIZE(tmat(ix)%from_vert)-1
                ixn=tmat(ix)%from_vert(ip)
                DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                  IF (tmat(ixn)%from_vert(ipn)==ix) THEN
                    tmat(ixn)%element(1,:,ipn)=
     $                 tmat(ixn)%element(1,:,ipn)
     $                -tmat(ixn)%element(3,:,ipn)*jcomp
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
c-----------------------------------------------------------------------
c           eliminate couplings.
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              tmat(ix)%element(:,:,ip)=
     $          tmat(ix)%element(:,:,ip)*mult
              ixn=tmat(ix)%from_vert(ip)
              DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                IF (tmat(ixn)%from_vert(ipn)==ix)
     $            tmat(ixn)%element(:,:,ipn)=
     $              tmat(ixn)%element(:,:,ipn)*multj
              ENDDO
            ENDDO
            IF (mat%foff==0) THEN
              IF (PRESENT(dscale)) THEN
                diagp=dscale
              ELSE
                diagp=0 
                DO jq=0,clim-1
                  diagp=MAX(diagp,tmat(ix)%element(jq,jq,0))
                ENDDO
                IF (diagp==0) diagp=1
              ENDIF
              DO jq=1,clim
                tmat(ix)%element(jq,jq,0)=
     $            tmat(ix)%element(jq,jq,0)+(1-mult(jq,jq))*diagp
              ENDDO
            ENDIF
          ENDIF block_type
c-----------------------------------------------------------------------
c       close outer loops.
c-----------------------------------------------------------------------
        ENDDO vert
      ENDDO block
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_op
c-----------------------------------------------------------------------
c     subprogram 6. regular_comp_op.
c     apply regularity conditions to R=0 points of an operator.
c     the action on the matrix is (I-uu).M.(I-uu)+uu, where u is the 
c     unit vector for the components that need elimination.  the
c     vcomp list of component descriptions in the matrix structure
c     and the fcomp Fourier component index are used to decide
c     where to apply this action.
c
c     this version handles complex 3-vectors directly.
c
c-PRE
c     this version works for multiple 3-vectors in a system, but that
c     generalization is not extended to tblocks.
c
c     if dscale is provided, it is used to scale the added uu entry, 
c     instead of determining the scaling on the fly.
c-----------------------------------------------------------------------
      SUBROUTINE regular_comp_op(mat,dscale)
      USE seam_storage_mod

      TYPE(complex_matrix_type), INTENT(INOUT) :: mat
      REAL(r8), INTENT(IN), OPTIONAL :: dscale

      INTEGER(i4) :: iv,nv,ip,ipn,ibl,ibe,ibv,jxmin,jxmax,jymin,jymax,
     $               jx,jy,ix,iy,ijx,ijy,mx,my,iq,jq,clim,ixn,
     $               ix0,iy0,jx0,jy0,imat,jmat,ivp,nvec,ivec,iv0,iv1,
     $               iqt,jqt,iqp,jqp,itype,jtype,iq0,iq1,jq0,jq1,jcomp
      COMPLEX(r8), DIMENSION(mat%nqty,mat%nqty) :: mult,multj
      REAL(r8), DIMENSION(mat%nqty) :: mult_v,mult_j
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      REAL(r8) :: diagp
      TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: tmat
c-----------------------------------------------------------------------
c     apply regularity conditions to the different vector components
c     of a matrix.  the elements at R=0 are uncoupled for n>0, if the
c     matrix acts on a scalar field.  the elements at R=0 of a matrix 
c     for a vector field are uncoupled for:
c       n=0:  r and phi
c       n=1:  z
c       n>1:  r, z, and phi
c     a covariant phi component (a_phi*r) is uncoupled for all n.  these
c     components are indicated with 'c'.
c-----------------------------------------------------------------------
      clim=mat%nqty
      nvec=clim/3
c-----------------------------------------------------------------------
c     create an array of 1s and 0s to zero out the appropriate
c     components.  the geometry is always toroidal at this point.
c
c     the second array, mult_j, is used for the columns.
c     [preconditioning may require matrix elements with different
c     Fourier indices for rows and columns.]
c-----------------------------------------------------------------------
      mult_v=0
      SELECT CASE(mat%fcomp)
      CASE(0)
        WHERE(mat%vcomp=='s') mult_v=1
        WHERE(mat%vcomp=='z') mult_v=1
      CASE(1)
        WHERE(mat%vcomp=='r') mult_v=1
      END SELECT
      mult_j=0
      jcomp=mat%fcomp+mat%foff
      SELECT CASE(ABS(jcomp))
      CASE(0)
        WHERE(mat%vcomp=='s') mult_j=1
        WHERE(mat%vcomp=='z') mult_j=1
      CASE(1)
        WHERE(mat%vcomp=='r') mult_j=1
      END SELECT
      DO jq=1,clim
        mult(jq,:)=mult_v
        multj(:,jq)=mult_j
      ENDDO
c-----------------------------------------------------------------------
c     loop over the borders of blocks touching R=0, and decouple the
c     appropriate matrix elements.
c-----------------------------------------------------------------------
      block: DO ibe=1,SIZE(r0block_list)
        ibl=r0block_list(ibe)
        nv=seam(ibl)%nvert
c-----------------------------------------------------------------------
c       loop over the block boundary.
c-----------------------------------------------------------------------
        vert: DO iv=1,nv
          ivp=iv-1
          IF (ivp==0) ivp=nv
c-----------------------------------------------------------------------
c         begin rblock-specific coding.
c-----------------------------------------------------------------------
          block_type: IF (ibl<=SIZE(mat%rbl_mat)) THEN
            mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
            my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
c-----------------------------------------------------------------------
c           loop over couplings from a given node at R=0 for n=1
c           vectors only.
c-----------------------------------------------------------------------
            IF (ABS(jcomp)==1.AND.MODULO(clim,3_i4)==0) THEN
              DO itype=1,mat%rbl_mat(ibl)%nbtype
                ix0=mat%rbl_mat(ibl)%ix0(itype)
                iy0=mat%rbl_mat(ibl)%iy0(itype)
                DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                  rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                  jx0=mat%rbl_mat(ibl)%ix0(jtype)
                  jy0=mat%rbl_mat(ibl)%iy0(jtype)
                  DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                    IF (jtype==1) THEN
                      IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                      ix=seam(ibl)%vertex(iv)%intxy(1)
                      iy=seam(ibl)%vertex(iv)%intxy(2)
                    ELSE
                      IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                          seam(ibl)%r0point(ivp))) CYCLE
                      IF (seam(ibl)%segment(iv)%h_side) THEN
                        IF (jtype>2) CYCLE
                      ELSE
                        IF (jtype<3) CYCLE
                      ENDIF
                      ix=seam(ibl)%segment(iv)%intxys(1)
                      iy=seam(ibl)%segment(iv)%intxys(2)
                    ENDIF
                    jxmin=MAX(ix0-1,ix0-ix)
                    jxmax=MIN(1-jx0,mx -ix)
                    jymin=MAX(iy0-1,iy0-iy)
                    jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                   combine r and phi equations to find averages in the
c                   matrix that enforce Vec_phi=i*Vec_r for n=1.  for
c                   preconditioning, we may also need Vec_phi=-i*Vec_r
c                   for n=-1.
c
c                   mathematically, the steps are 1) change variables 
c                   to X1=(Vec_r+i*Vec_phi)/2 and X2=(Vec_r-i*Vec_phi)/2
c                   for each vertex at R=0, 2) set X1=0 (remove
c                   couplings to and from X1) giving an overdetermined
c                   system of equations, 3) add -i*(the X2-equation) to
c                   the X1 equation.
c
c                   for the n=-1 preconditioning matrix case, we will
c                   set X2 to zero and keep X1 in the r-comp location.
c-----------------------------------------------------------------------
                    DO ivec=0,nvec-1
                      jq0=(jmat-1)*clim+3*ivec+1
                      jq1=jq0+2
                      DO jy=jymin,jymax
                        ijy=iy+jy
                        DO jx=jxmin,jxmax
                          ijx=ix+jx
                          rmat(jq0,-jx,-jy,:,ijx,ijy)=
     $                             rmat(jq0,-jx,-jy,:,ijx,ijy)
     $                      +(0,1)*rmat(jq1,-jx,-jy,:,ijx,ijy)*jcomp
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
c-----------------------------------------------------------------------
c           loop over couplings to a given node at R=0.
c-----------------------------------------------------------------------
            DO itype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,mat%rbl_mat(ibl)%nbtype
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO imat=1,mat%rbl_mat(ibl)%nb_type(itype)
                  IF (itype==1) THEN
                    IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                  ELSE
                    IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                        seam(ibl)%r0point(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (itype>2) CYCLE
                    ELSE
                      IF (itype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                  ENDIF
                  jxmin=MAX(jx0-1,jx0-ix)
                  jxmax=MIN(1-ix0,mx -ix)
                  jymin=MAX(jy0-1,jy0-iy)
                  jymax=MIN(1-iy0,my -iy)
c-----------------------------------------------------------------------
c                 combine r and phi equations to find averages in the
c                 matrix that enforce Vec_phi=i*Vec_r for n=1.  [note
c                 rows always represent n>=0.]
c
c                 then, eliminate couplings to this node.
c-----------------------------------------------------------------------
                  iq0=(imat-1)*clim+1
                  iq1=imat*clim
                  DO jy=jymin,jymax
                    DO jx=jxmin,jxmax
                      IF (mat%fcomp==1.AND.MODULO(clim,3_i4)==0) THEN
                        DO ivec=0,nvec-1
                          iv0=(imat-1)*clim+3*ivec+1
                          iv1=iv0+2
                          rmat(:,jx,jy,iv0,ix,iy)=
     $                            rmat(:,jx,jy,iv0,ix,iy)
     $                     -(0,1)*rmat(:,jx,jy,iv1,ix,iy)
                        ENDDO
                      ENDIF
                      DO jq=1,mat%rbl_mat(ibl)%nq_type(jtype)
                        rmat(jq,jx,jy,iq0:iq1,ix,iy)=
     $                    rmat(jq,jx,jy,iq0:iq1,ix,iy)*mult_v
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           again loop over couplings from a given node at R=0.
c-----------------------------------------------------------------------
            DO itype=1,mat%rbl_mat(ibl)%nbtype
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                  IF (jtype==1) THEN
                    IF (.NOT.seam(ibl)%r0point(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                  ELSE
                    IF (.NOT.(seam(ibl)%r0point(iv).AND.
     $                        seam(ibl)%r0point(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (jtype>2) CYCLE
                    ELSE
                      IF (jtype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                  ENDIF
                  jxmin=MAX(ix0-1,ix0-ix)
                  jxmax=MIN(1-jx0,mx -ix)
                  jymin=MAX(iy0-1,iy0-iy)
                  jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                 eliminate couplings from the node on R=0.
c-----------------------------------------------------------------------
                  jq0=(jmat-1)*clim+1
                  jq1=jmat*clim
                  DO jy=jymin,jymax
                    ijy=iy+jy
                    DO jx=jxmin,jxmax
                      ijx=ix+jx
                      DO iq=1,mat%rbl_mat(ibl)%nq_type(itype)
                        rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)=
     $                    rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)*mult_j
                      ENDDO
                    ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c                 take care of uu and zz.
c-----------------------------------------------------------------------
                  IF (jtype==itype.AND.mat%foff==0) THEN
                    IF (PRESENT(dscale)) THEN
                      diagp=dscale
                    ELSE
                      diagp=0 
                      DO jq=0,clim-1
                        diagp=MAX(diagp,
     $                            REAL(rmat(jq+jq0,0,0,jq+jq0,ix,iy)))
                      ENDDO
                      IF (diagp==0) diagp=1
                    ENDIF
                    DO jq=0,clim-1
                      rmat(jq+jq0,0,0,jq+jq0,ix,iy)=
     $                  rmat(jq+jq0,0,0,jq+jq0,ix,iy)
     $                 +(1-mult_v(jq+1))*diagp
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c         begin tblock-specific coding.
c-----------------------------------------------------------------------
          ELSE block_type
            tmat=>mat%tbl_mat(ibl)%lmat
c-----------------------------------------------------------------------
c           combine r and phi equations to find averages in the matrix
c           that enforce Vec_phi=i*Vec_r for n=1.  see the comment
c           above.
c-----------------------------------------------------------------------
            IF (mat%fcomp==1.AND.clim==3) THEN
              DO ip=0,SIZE(tmat(ix)%from_vert)-1
                  tmat(ix)%element(:,1,ip)=tmat(ix)%element(:,1,ip)
     $                              -(0,1)*tmat(ix)%element(:,3,ip)
              ENDDO
            ENDIF
            IF (ABS(jcomp)==1.AND.clim==3) THEN
              DO ip=0,SIZE(tmat(ix)%from_vert)-1
                ixn=tmat(ix)%from_vert(ip)
                DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                  IF (tmat(ixn)%from_vert(ipn)==ix) THEN
                      tmat(ixn)%element(1,:,ipn)=
     $                         tmat(ixn)%element(1,:,ipn)
     $                  +(0,1)*tmat(ixn)%element(3,:,ipn)*jcomp
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
c-----------------------------------------------------------------------
c           eliminate couplings.
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              tmat(ix)%element(:,:,ip)=
     $          tmat(ix)%element(:,:,ip)*mult
              ixn=tmat(ix)%from_vert(ip)
              DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                IF (tmat(ixn)%from_vert(ipn)==ix)
     $            tmat(ixn)%element(:,:,ipn)=
     $              tmat(ixn)%element(:,:,ipn)*multj
              ENDDO
            ENDDO
            IF (mat%foff==0) THEN
              IF (PRESENT(dscale)) THEN
                diagp=dscale
              ELSE
                diagp=0 
                DO jq=0,clim-1
                  diagp=MAX(diagp,REAL(tmat(ix)%element(jq,jq,0)))
                ENDDO
                IF (diagp==0) diagp=1
              ENDIF
              DO jq=1,clim
                tmat(ix)%element(jq,jq,0)=
     $            tmat(ix)%element(jq,jq,0)+(1-mult(jq,jq))*diagp
              ENDDO
            ENDIF
          ENDIF block_type
c-----------------------------------------------------------------------
c       close outer loops.
c-----------------------------------------------------------------------
        ENDDO vert
      ENDDO block
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE regular_comp_op
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE regularity
