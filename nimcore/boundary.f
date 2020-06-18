c-----------------------------------------------------------------------
c     file boundary.f:  contains all subprograms related to the external
c     boundary conditions for nimrod.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declaration.
c-----------------------------------------------------------------------
      MODULE boundary
      USE local
      IMPLICIT NONE

      INTEGER(i4), PARAMETER, PRIVATE :: max3v=32_i4
      INTEGER(i4), DIMENSION(max3v), PRIVATE :: start3v

      CONTAINS
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. dirichlet_rhs.
c     2. no_rhs_bc.
c     3. dirichlet_op.
c     4. no_mat_bc.
c     5. dirichlet_comp_op.
c     6. no_comp_mat_bc.
c     7. bcflag_parse
c     8. bcdir_set
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 1. dirichlet_rhs.
c     apply external boundary conditions to blocks in exblock_list.
c
c     when the optional svess parameter is provided, the data that is
c     subtracted from rhs is saved in the seam_csave data structure.
c     note that while svess is a character, only its presence is tested.
c-----------------------------------------------------------------------
      SUBROUTINE dirichlet_rhs(rhs,blseam,component,nqty,symm,svess)
      USE edge_type_mod
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: nqty
      TYPE(cvector_type), INTENT(INOUT) :: rhs
      TYPE(edge_type), INTENT(IN) :: blseam
      CHARACTER(*), INTENT(IN) :: component
      CHARACTER(*), INTENT(IN), OPTIONAL :: symm,svess

      COMPLEX(r8), DIMENSION(nqty) :: cproj
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: arr5d

      REAL(r8), DIMENSION(nqty,nqty) :: bcpmat
      INTEGER(i4), DIMENSION(nqty,2) :: bciarr
      INTEGER(i4) :: iv,ix,iy,imode,ivp,nside,is,isymm,iq,nm
c-----------------------------------------------------------------------
c     set the symmetry flag.  if symm starts with "t", the top boundary
c     is a symmetry condition.  if symm starts with "b", the bottom
c     boundary is a symmetry condition.
c-----------------------------------------------------------------------
      isymm=0
      IF (PRESENT(symm)) THEN
        SELECT CASE(symm(1:1))
        CASE('t','T')
          isymm=1
        CASE('b','B')
          isymm=-1
        END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     parse the component flag to create an integer array, which
c     indicates which scalar and 3-vector components have essential
c     conditions.
c-----------------------------------------------------------------------
      CALL bcflag_parse(component,nqty,bciarr)
c-----------------------------------------------------------------------
c     the bcdir_set routine combines the bciarr information with local
c     surface-normal and tangential unit directions.  mesh vertices
c     are treated first.
c-----------------------------------------------------------------------
      nm=SIZE(rhs%arr,4)
      vert: DO iv=1,blseam%nvert
        IF (.NOT.blseam%expoint(iv)) CYCLE
        ix=blseam%vertex(iv)%intxy(1)
        iy=blseam%vertex(iv)%intxy(2)
        CALL bcdir_set(nqty,bciarr,blseam%excorner(iv),isymm,
     $                 blseam%vertex(iv)%norm,blseam%vertex(iv)%tang,
     $                 bcpmat)
        IF (PRESENT(svess)) THEN
          iq=1
          DO imode=1,nm
            cproj=MATMUL(bcpmat,rhs%arr(1:nqty,ix,iy,imode))
            rhs%arr(1:nqty,ix,iy,imode)=
     $        rhs%arr(1:nqty,ix,iy,imode)-cproj
            blseam%vertex(iv)%seam_csave(iq:iq+nqty-1)=
     $        blseam%vertex(iv)%seam_csave(iq:iq+nqty-1)+cproj
            iq=iq+nqty
          ENDDO
        ELSE
          DO imode=1,nm
            cproj=MATMUL(bcpmat,rhs%arr(1:nqty,ix,iy,imode))
            rhs%arr(1:nqty,ix,iy,imode)=
     $        rhs%arr(1:nqty,ix,iy,imode)-cproj
          ENDDO
        ENDIF
      ENDDO vert
c-----------------------------------------------------------------------
c     side-centered nodes.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arrh)) THEN
        nside=SIZE(rhs%arrh,2)
        seg: DO iv=1,blseam%nvert
          ivp=iv-1
          IF (ivp==0) ivp=blseam%nvert
          IF (.NOT.(blseam%expoint(iv).AND.blseam%expoint(ivp)))
     $      CYCLE seg
          ix=blseam%segment(iv)%intxys(1)
          iy=blseam%segment(iv)%intxys(2)
          IF (blseam%segment(iv)%h_side) THEN
            arr5d=>rhs%arrh
          ELSE
            arr5d=>rhs%arrv
          ENDIF

          IF (PRESENT(svess)) THEN
            iq=1
            DO imode=1,nm
              DO is=1,nside
                CALL bcdir_set(nqty,bciarr,.false.,isymm,
     $                   blseam%segment(iv)%norm(:,is),
     $                   blseam%segment(iv)%tang(:,is),bcpmat)
                cproj=MATMUL(bcpmat,arr5d(1:nqty,is,ix,iy,imode))
                arr5d(1:nqty,is,ix,iy,imode)=
     $            arr5d(1:nqty,is,ix,iy,imode)-cproj
                blseam%segment(iv)%seam_csave(iq:iq+nqty-1)=
     $            blseam%segment(iv)%seam_csave(iq:iq+nqty-1)+cproj
                iq=iq+nqty
              ENDDO
            ENDDO
          ELSE
            DO imode=1,nm
              DO is=1,nside
                CALL bcdir_set(nqty,bciarr,.false.,isymm,
     $                   blseam%segment(iv)%norm(:,is),
     $                   blseam%segment(iv)%tang(:,is),bcpmat)
                cproj=MATMUL(bcpmat,arr5d(1:nqty,is,ix,iy,imode))
                arr5d(1:nqty,is,ix,iy,imode)=
     $            arr5d(1:nqty,is,ix,iy,imode)-cproj
              ENDDO
            ENDDO
          ENDIF
        ENDDO seg
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dirichlet_rhs
c-----------------------------------------------------------------------
c     subprogram 2. no_rhs_bc.
c     a dummy boundary condition routine with the same interface
c     block as dirichlet_rhs (and future boundary condition routines).
c-----------------------------------------------------------------------
      SUBROUTINE no_rhs_bc(rhs,blseam,component,nqty,symm,svess)
      USE edge_type_mod
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: nqty
      TYPE(cvector_type), INTENT(INOUT) :: rhs
      TYPE(edge_type), INTENT(IN) :: blseam
      CHARACTER(*), INTENT(IN) :: component
      CHARACTER(*), INTENT(IN), OPTIONAL :: symm,svess

      RETURN
      END SUBROUTINE no_rhs_bc
c-----------------------------------------------------------------------
c     subprogram 3. dirichlet_op.
c     apply Dirichlet boundary conditions to a Cartesian operator.
c     If the specified component is tangent, the resultant matrix is 
c     (I-tt-zz).M.(I-tt-zz)+tt+zz, where t is the surface tangent in 
c     the computational plane, z is the unit vector normal to the plane,
c     and I is the identity matrix.  If the specified component is 
c     normal, the resultant matrix is (I-nn).M.(I-nn)+nn, where n is 
c     the surface normal.  If the input is all, then the result is
c     (I-tt-zz-nn).M.(I-tt-zz-nn)+tt+zz+nn.
c
c     when the end of the component parameter is 'offdiag,' the passed
c     matrix structure contains only off-diagonal matrix entries, and
c     there is no diagonal entry to enter after couplings are
c     eliminated.
c-----------------------------------------------------------------------
      SUBROUTINE dirichlet_op(mat,component,dscale,symm)
      USE seam_storage_mod
      USE matrix_type_mod

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      CHARACTER(*), INTENT(IN) :: component
      REAL(r8), INTENT(IN) :: dscale
      CHARACTER(*), INTENT(IN), OPTIONAL :: symm

      INTEGER(i4) :: iv,nv,ip,ipn,ibl,ibe,ibv,jxmin,jxmax,jymin,jymax,
     $               jx,jy,ix,iy,ijx,ijy,mx,my,iq,jq,clim,ixn,
     $               imat,jmat,ivp,ix0,iy0,jx0,jy0,ibase_st,nqtyp,
     $               itype,jtype,iq0,iq1,jq0,jq1,is,js,isymm
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      TYPE(matrix_element_type3), DIMENSION(:), POINTER :: tmat
      LOGICAL :: diag

      REAL(r8), DIMENSION(:), ALLOCATABLE :: proj
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: bcpmat
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: bciarr
c-----------------------------------------------------------------------
c     check vector component range.
c-----------------------------------------------------------------------
      IF (SIZE(mat%rbl_mat)>0) THEN
        clim=SIZE(mat%rbl_mat(1)%mat(1,1)%arr,1)
      ELSE
        clim=SIZE(mat%tbl_mat(1)%lmat(0)%element,1)
      ENDIF
      ALLOCATE(proj(clim),bcpmat(clim,clim),bciarr(clim,2))
c-----------------------------------------------------------------------
c     set the symmetry flag.  if symm starts with "t", the top boundary
c     is a symmetry condition.  if symm starts with "b", the bottom
c     boundary is a symmetry condition.
c-----------------------------------------------------------------------
      isymm=0
      IF (PRESENT(symm)) THEN
        SELECT CASE(symm(1:1))
        CASE('t','T')
          isymm=1
        CASE('b','B')
          isymm=-1
        END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     parse the component flag to create an integer array, which
c     indicates which scalar and 3-vector components have essential
c     conditions.
c-----------------------------------------------------------------------
      CALL bcflag_parse(component,clim,bciarr,diag)
c-----------------------------------------------------------------------
c     save the diagonal scaling and essential-condition flag.
c-----------------------------------------------------------------------
      mat%diag_scale=dscale
      mat%essential_cond=component
c-----------------------------------------------------------------------
c     loop over all external boundary points and zero the couplings
c     to the specified components along the boundary.
c-----------------------------------------------------------------------
      block: DO ibe=1,SIZE(exblock_list)
        ibl=exblock_list(ibe)
        nv=seam(ibl)%nvert
c-----------------------------------------------------------------------
c       begin rblock-specific coding.  rblock basis types indices are
c       1 == grid vertex centered
c       2 == horizontal side centered
c       3 == vertical side centered
c       4 == interior centered
c       1 : 3 are affected by boundary conditions.
c       note that itype is the to basis type and jtype is the from
c       basis type.
c-----------------------------------------------------------------------
        block_type: IF (ibl<=SIZE(mat%rbl_mat)) THEN
c-----------------------------------------------------------------------
c         connections to the rblock boundary vertices.
c-----------------------------------------------------------------------
          rb_vert: DO iv=1,nv
            ivp=iv-1
            IF (ivp==0) ivp=seam(ibl)%nvert
            mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
            my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
c-----------------------------------------------------------------------
c           eliminate couplings from the specified components only.
c-----------------------------------------------------------------------
            DO itype=1,mat%rbl_mat(ibl)%nbtype
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              nqtyp=mat%rbl_mat(ibl)%nq_type(itype)
              DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                  IF (jtype==1) THEN
                    IF (.NOT.seam(ibl)%expoint(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                    CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),
     $                             isymm,seam(ibl)%vertex(iv)%norm,
     $                             seam(ibl)%vertex(iv)%tang,bcpmat)
                  ELSE
                    IF (.NOT.(seam(ibl)%expoint(iv).AND.
     $                        seam(ibl)%expoint(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (jtype>2) CYCLE
                    ELSE
                      IF (jtype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                    CALL bcdir_set(clim,bciarr,.false.,isymm,
     $                        seam(ibl)%segment(iv)%norm(:,jmat),
     $                        seam(ibl)%segment(iv)%tang(:,jmat),bcpmat)
                  ENDIF
                  jxmin=MAX(ix0-1,ix0-ix)
                  jxmax=MIN(1-jx0,mx -ix)
                  jymin=MAX(iy0-1,iy0-iy)
                  jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                 find the vector elements of M.uu then subtract 
c                 M.uu from M
c-----------------------------------------------------------------------
                  jq0=(jmat-1)*clim+1
                  jq1=jmat*clim
                  DO jy=jymin,jymax
                    ijy=iy+jy
                    DO jx=jxmin,jxmax
                      ijx=ix+jx
                      DO iq=1,nqtyp
                        proj=
     $                   MATMUL(bcpmat,rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy))
                        rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)=
     $                    rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)-proj
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           eliminate couplings to the specified component only.
c-----------------------------------------------------------------------
            DO itype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,mat%rbl_mat(ibl)%nbtype
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                nqtyp=mat%rbl_mat(ibl)%nq_type(jtype)
                DO imat=1,mat%rbl_mat(ibl)%nb_type(itype)
                  IF (itype==1) THEN
                    IF (.NOT.seam(ibl)%expoint(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                    CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),
     $                             isymm,seam(ibl)%vertex(iv)%norm,
     $                             seam(ibl)%vertex(iv)%tang,bcpmat)
                  ELSE
                    IF (.NOT.(seam(ibl)%expoint(iv).AND.
     $                        seam(ibl)%expoint(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (itype>2) CYCLE
                    ELSE
                      IF (itype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                    CALL bcdir_set(clim,bciarr,.false.,isymm,
     $                        seam(ibl)%segment(iv)%norm(:,imat),
     $                        seam(ibl)%segment(iv)%tang(:,imat),bcpmat)
                  ENDIF
                  jxmin=MAX(jx0-1,jx0-ix)
                  jxmax=MIN(1-ix0,mx -ix)
                  jymin=MAX(jy0-1,jy0-iy)
                  jymax=MIN(1-iy0,my -iy)
c-----------------------------------------------------------------------
c                 find the vector elements of uu.M then subtract 
c                 uu.M from M
c-----------------------------------------------------------------------
                  iq0=(imat-1)*clim+1
                  iq1=imat*clim
                  DO jy=jymin,jymax
                    DO jx=jxmin,jxmax
                      DO jq=1,nqtyp
                        proj=
     $                   MATMUL(bcpmat,rmat(jq,jx,jy,iq0:iq1,ix,iy))
                        rmat(jq,jx,jy,iq0:iq1,ix,iy)=
     $                    rmat(jq,jx,jy,iq0:iq1,ix,iy)-proj
                      ENDDO
                    ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c                 take care of uu and zz.
c-----------------------------------------------------------------------
                  IF (jtype==itype.AND.diag) THEN
                    rmat(iq0:iq1,0,0,iq0:iq1,ix,iy)=
     $                rmat(iq0:iq1,0,0,iq0:iq1,ix,iy)+bcpmat*dscale
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO rb_vert
c-----------------------------------------------------------------------
c       begin tblock-specific coding.
c-----------------------------------------------------------------------
        ELSE block_type
          tb_vert: DO iv=1,nv
            IF (.NOT.seam(ibl)%expoint(iv)) CYCLE tb_vert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            tmat=>mat%tbl_mat(ibl)%lmat
            CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),isymm,
     $                     seam(ibl)%vertex(iv)%norm,
     $                     seam(ibl)%vertex(iv)%tang,bcpmat)
c-----------------------------------------------------------------------
c           find the vector elements of M.uu then subtract 
c           M.uu from M
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              ixn=tmat(ix)%from_vert(ip)
              DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                IF (tmat(ixn)%from_vert(ipn)==ix) THEN
                  DO iq=1,clim
                    proj=MATMUL(bcpmat,tmat(ixn)%element(:,iq,ipn))
                    tmat(ixn)%element(:,iq,ipn)=
     $                tmat(ixn)%element(:,iq,ipn)-proj
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           find the vector elements of u.M then subtract 
c           uu.M and (and zz.M for specified tangent) from M
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              DO jq=1,clim
                proj=MATMUL(bcpmat,tmat(ixn)%element(jq,:,ip))
                tmat(ixn)%element(jq,:,ip)=
     $            tmat(ixn)%element(jq,:,ip)-proj
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           take care of uu and zz.
c-----------------------------------------------------------------------
            IF (diag) THEN
              tmat(ix)%element(:,:,0)=
     $          tmat(ix)%element(:,:,0)+bcpmat*dscale
            ENDIF
          ENDDO tb_vert
        ENDIF block_type
      ENDDO block
      DEALLOCATE(proj,bcpmat,bciarr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dirichlet_op
c-----------------------------------------------------------------------
c     subprogram 4. no_mat_bc.
c     a dummy boundary condition routine with the same interface
c     block as dirichlet_op (and future boundary condition routines).
c-----------------------------------------------------------------------
      SUBROUTINE no_mat_bc(mat,component,dscale,symm)
      USE seam_storage_mod
      USE matrix_type_mod

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      CHARACTER(*), INTENT(IN) :: component
      REAL(r8), INTENT(IN) :: dscale
      CHARACTER(*), INTENT(IN), OPTIONAL :: symm

      RETURN
      END SUBROUTINE no_mat_bc
c-----------------------------------------------------------------------
c     subprogram 5. dirichlet_comp_op.
c     apply Dirichlet boundary conditions to a Cartesian operator.
c     If the specified component is tangent, the resultant matrix is 
c     (I-tt-zz).M.(I-tt-zz)+tt+zz, where t is the surface tangent in 
c     the computational plane, z is the unit vector normal to the plane,
c     and I is the identity matrix.  If the specified component is 
c     normal, the resultant matrix is (I-nn).M.(I-nn)+nn, where n is 
c     the surface normal.  If the input is all, then the result is
c     (I-tt-zz-nn).M.(I-tt-zz-nn)+tt+zz+nn.
c
c     when the end of the component parameter is 'offdiag,' the passed
c     matrix structure contains only off-diagonal matrix entries, and
c     there is no diagonal entry to enter after couplings are
c     eliminated.
c
c     this generalized version handles combinations of 3-vectors and
c     scalars.
c-----------------------------------------------------------------------
      SUBROUTINE dirichlet_comp_op(mat,component,dscale)
      USE seam_storage_mod
      USE matrix_type_mod

      TYPE(complex_matrix_type), INTENT(INOUT) :: mat
      CHARACTER(*), INTENT(IN) :: component
      REAL(r8), INTENT(IN) :: dscale

      INTEGER(i4) :: iv,nv,ip,ipn,ibl,ibe,ibv,jxmin,jxmax,jymin,jymax,
     $               jx,jy,ix,iy,ijx,ijy,mx,my,iq,jq,clim,ixn,
     $               imat,jmat,ivp,ix0,iy0,jx0,jy0,
     $               itype,jtype,iq0,iq1,jq0,jq1,is,js,nqtyp
      REAL(r8), DIMENSION(2) :: unit
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: tmat
      LOGICAL :: diag

      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: cproj
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: bcpmat
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: bciarr
c-----------------------------------------------------------------------
c     check vector component range.
c-----------------------------------------------------------------------
      IF (SIZE(mat%rbl_mat)>0) THEN
        clim=SIZE(mat%rbl_mat(1)%mat(1,1)%arr,1)
      ELSE
        clim=SIZE(mat%tbl_mat(1)%lmat(0)%element,1)
      ENDIF
      ALLOCATE(cproj(clim),bcpmat(clim,clim),bciarr(clim,2))
c-----------------------------------------------------------------------
c     parse the component flag to create an integer array, which
c     indicates which scalar and 3-vector components have essential
c     conditions.
c-----------------------------------------------------------------------
      CALL bcflag_parse(component,clim,bciarr,diag)
c-----------------------------------------------------------------------
c     save the diagonal scaling and essential-condition flag.
c-----------------------------------------------------------------------
      mat%diag_scale=dscale
      mat%essential_cond=component
c-----------------------------------------------------------------------
c     loop over all external boundary points and zero the couplings
c     to the specified components along the boundary.
c-----------------------------------------------------------------------
      block: DO ibe=1,SIZE(exblock_list)
        ibl=exblock_list(ibe)
        nv=seam(ibl)%nvert
c-----------------------------------------------------------------------
c       begin rblock-specific coding.  rblock basis types indices are
c       1 == grid vertex centered
c       2 == horizontal side centered
c       3 == vertical side centered
c       4 == interior centered
c       1 : 3 are affected by boundary conditions.
c       note that itype is the to basis type and jtype is the from
c       basis type.
c-----------------------------------------------------------------------
        block_type: IF (ibl<=SIZE(mat%rbl_mat)) THEN
c-----------------------------------------------------------------------
c         connections to the rblock boundary vertices.
c-----------------------------------------------------------------------
          rb_vert: DO iv=1,nv
            ivp=iv-1
            IF (ivp==0) ivp=seam(ibl)%nvert
            mx=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,5)-1
            my=SIZE(mat%rbl_mat(ibl)%mat(1,1)%arr,6)-1
c-----------------------------------------------------------------------
c           eliminate couplings from the specified components only.
c-----------------------------------------------------------------------
            DO itype=1,mat%rbl_mat(ibl)%nbtype
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              nqtyp=mat%rbl_mat(ibl)%nq_type(itype)
              DO jtype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                DO jmat=1,mat%rbl_mat(ibl)%nb_type(jtype)
                  IF (jtype==1) THEN
                    IF (.NOT.seam(ibl)%expoint(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                    CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),
     $                             0_i4,seam(ibl)%vertex(iv)%norm,
     $                             seam(ibl)%vertex(iv)%tang,bcpmat)
                  ELSE
                    IF (.NOT.(seam(ibl)%expoint(iv).AND.
     $                        seam(ibl)%expoint(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (jtype>2) CYCLE
                    ELSE
                      IF (jtype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                    CALL bcdir_set(clim,bciarr,.false.,0_i4,
     $                        seam(ibl)%segment(iv)%norm(:,jmat),
     $                        seam(ibl)%segment(iv)%tang(:,jmat),bcpmat)
                  ENDIF
                  jxmin=MAX(ix0-1,ix0-ix)
                  jxmax=MIN(1-jx0,mx -ix)
                  jymin=MAX(iy0-1,iy0-iy)
                  jymax=MIN(1-jy0,my -iy)
c-----------------------------------------------------------------------
c                 find the vector elements of M.uu then subtract 
c                 M.uu from M
c-----------------------------------------------------------------------
                  jq0=(jmat-1)*clim+1
                  jq1=jmat*clim
                  DO jy=jymin,jymax
                    ijy=iy+jy
                    DO jx=jxmin,jxmax
                      ijx=ix+jx
                      DO iq=1,nqtyp
                        cproj=
     $                   MATMUL(bcpmat,rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy))
                        rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)=
     $                    rmat(jq0:jq1,-jx,-jy,iq,ijx,ijy)-cproj
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           eliminate couplings to the specified component only.
c-----------------------------------------------------------------------
            DO itype=1,MIN(mat%rbl_mat(ibl)%nbtype,3_i4)
              ix0=mat%rbl_mat(ibl)%ix0(itype)
              iy0=mat%rbl_mat(ibl)%iy0(itype)
              DO jtype=1,mat%rbl_mat(ibl)%nbtype
                rmat=>mat%rbl_mat(ibl)%mat(jtype,itype)%arr
                jx0=mat%rbl_mat(ibl)%ix0(jtype)
                jy0=mat%rbl_mat(ibl)%iy0(jtype)
                nqtyp=mat%rbl_mat(ibl)%nq_type(jtype)
                DO imat=1,mat%rbl_mat(ibl)%nb_type(itype)
                  IF (itype==1) THEN
                    IF (.NOT.seam(ibl)%expoint(iv)) CYCLE
                    ix=seam(ibl)%vertex(iv)%intxy(1)
                    iy=seam(ibl)%vertex(iv)%intxy(2)
                    CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),
     $                             0_i4,seam(ibl)%vertex(iv)%norm,
     $                             seam(ibl)%vertex(iv)%tang,bcpmat)
                  ELSE
                    IF (.NOT.(seam(ibl)%expoint(iv).AND.
     $                        seam(ibl)%expoint(ivp))) CYCLE
                    IF (seam(ibl)%segment(iv)%h_side) THEN
                      IF (itype>2) CYCLE
                    ELSE
                      IF (itype<3) CYCLE
                    ENDIF
                    ix=seam(ibl)%segment(iv)%intxys(1)
                    iy=seam(ibl)%segment(iv)%intxys(2)
                    CALL bcdir_set(clim,bciarr,.false.,0_i4,
     $                        seam(ibl)%segment(iv)%norm(:,imat),
     $                        seam(ibl)%segment(iv)%tang(:,imat),bcpmat)
                  ENDIF
                  jxmin=MAX(jx0-1,jx0-ix)
                  jxmax=MIN(1-ix0,mx -ix)
                  jymin=MAX(jy0-1,jy0-iy)
                  jymax=MIN(1-iy0,my -iy)
c-----------------------------------------------------------------------
c                 find the vector elements of uu.M then subtract 
c                 uu.M from M
c-----------------------------------------------------------------------
                  iq0=(imat-1)*clim+1
                  iq1=imat*clim
                  DO jy=jymin,jymax
                    DO jx=jxmin,jxmax
                      DO jq=1,nqtyp
                        cproj=
     $                   MATMUL(bcpmat,rmat(jq,jx,jy,iq0:iq1,ix,iy))
                        rmat(jq,jx,jy,iq0:iq1,ix,iy)=
     $                    rmat(jq,jx,jy,iq0:iq1,ix,iy)-cproj
                      ENDDO
                    ENDDO
                  ENDDO
c-----------------------------------------------------------------------
c                 take care of uu and zz.
c-----------------------------------------------------------------------
                  IF (jtype==itype.AND.diag) THEN
                    rmat(iq0:iq1,0,0,iq0:iq1,ix,iy)=
     $                rmat(iq0:iq1,0,0,iq0:iq1,ix,iy)+bcpmat*dscale
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO rb_vert
c-----------------------------------------------------------------------
c       begin tblock-specific coding.
c-----------------------------------------------------------------------
        ELSE block_type
          tb_vert: DO iv=1,nv
            IF (.NOT.seam(ibl)%expoint(iv)) CYCLE tb_vert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            tmat=>mat%tbl_mat(ibl)%lmat
            CALL bcdir_set(clim,bciarr,seam(ibl)%excorner(iv),0_i4,
     $                     seam(ibl)%vertex(iv)%norm,
     $                     seam(ibl)%vertex(iv)%tang,bcpmat)
c-----------------------------------------------------------------------
c           find the vector elements of M.uu then subtract 
c           M.uu from M
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              ixn=tmat(ix)%from_vert(ip)
              DO ipn=0,SIZE(tmat(ixn)%from_vert)-1
                IF (tmat(ixn)%from_vert(ipn)==ix) THEN
                  DO iq=1,clim
                    cproj=MATMUL(bcpmat,tmat(ixn)%element(:,iq,ipn))
                    tmat(ixn)%element(:,iq,ipn)=
     $                tmat(ixn)%element(:,iq,ipn)-cproj
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           find the vector elements of u.M then subtract 
c           uu.M and (and zz.M for specified tangent) from M
c-----------------------------------------------------------------------
            DO ip=0,SIZE(tmat(ix)%from_vert)-1
              DO jq=1,clim
                cproj=MATMUL(bcpmat,tmat(ixn)%element(jq,:,ip))
                tmat(ixn)%element(jq,:,ip)=
     $            tmat(ixn)%element(jq,:,ip)-cproj
              ENDDO
            ENDDO
c-----------------------------------------------------------------------
c           take care of uu and zz.
c-----------------------------------------------------------------------
            IF (diag) THEN
              tmat(ix)%element(:,:,0)=
     $          tmat(ix)%element(:,:,0)+bcpmat*dscale
            ENDIF
          ENDDO tb_vert
        ENDIF block_type
      ENDDO block
      DEALLOCATE(cproj,bcpmat,bciarr)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE dirichlet_comp_op
c-----------------------------------------------------------------------
c     subprogram 6. no_comp_mat_bc.
c     a dummy boundary condition routine with the same interface
c     block as dirichlet_op (and future boundary condition routines).
c-----------------------------------------------------------------------
      SUBROUTINE no_comp_mat_bc(mat,component,dscale)
      USE seam_storage_mod
      USE matrix_type_mod

      TYPE(complex_matrix_type), INTENT(INOUT) :: mat
      CHARACTER(*), INTENT(IN) :: component
      REAL(r8), INTENT(IN) :: dscale

      RETURN
      END SUBROUTINE no_comp_mat_bc
c-----------------------------------------------------------------------
c     subprogram 7. bcflag_parse.
c     this routine parses the component flag that is sent into one of
c     the dirichlet_ routines.  input character strings are expected
c     to have a set of substrings separated by spaces, where each
c     substring indicates which components of a 3-vector or whether
c     a scalar needs to be set to zero.  valid substrings and their
c     meanings are:
c
c       "3vn" - zero-out the normal component of a 3-vector
c       "3vt" - zero-out the tangential components of a 3-vector
c	"3vf" - 3-vector without a Dichlet condition is applied (free)
c       "sd"  - zero-out a scalar
c       "sf"  - scalar without a Dirichlet condition
c       
c     the output is the two-dimensional integer array bciarr that codes
c     which components have Dirichlet conditions.  bciarr(:,1) is the
c     coding for scalars, and bciarr(:,2) is the coding for vectors.
c     bciarr(:,1) is set as:
c
c       0 = do not alter (scalar or vector quantity index)
c       1 = set this scalar to zero
c
c     bciarr(:,2) has zero for each scalar quantity, and for each
c     3-vector, it can have the following sequences:
c
c       000 = do not alter
c       100 = set the normal component to zero
c       011 = set all tangential components to zero
c       010 = set the in-plane tangential components to zero (not used)
c       001 = set the the out-of plane (third) component to zero (ditto)
c
c     this routine also sets the module array start3v, which is a list
c     of starting indices for all 3-vectors within bciarr(:,2).
c-----------------------------------------------------------------------
      SUBROUTINE bcflag_parse(bcfl,nqty,bciarr,diag)
      USE seam_storage_mod
      USE matrix_type_mod

      CHARACTER(*), INTENT(IN) :: bcfl
      INTEGER(i4), INTENT(IN) :: nqty
      INTEGER(i4), DIMENSION(nqty,2), INTENT(OUT) :: bciarr
      LOGICAL, INTENT(OUT), OPTIONAL :: diag

      INTEGER(i4) :: ivec,ichar,lenstr,iind
c-----------------------------------------------------------------------
c     intialize local and output variables.
c-----------------------------------------------------------------------
      lenstr=LEN_TRIM(bcfl)
      ichar=0
      iind=0
      ivec=0
      bciarr=0
      start3v=0
c-----------------------------------------------------------------------
c     set the diag logical for calls from routines for matrices.  the
c     choice is based on the last characters of bcfl.
c-----------------------------------------------------------------------
      IF (PRESENT(diag)) THEN
        diag=.true.
        IF (lenstr>6) THEN
          IF (bcfl(lenstr-6:lenstr)=="offdiag") THEN
            diag=.false.
            lenstr=lenstr-7
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     if bcflag indicates that all components have Dirichlet conditions,
c     treat all system components as scalars.
c-----------------------------------------------------------------------
      IF (lenstr>=3) THEN
        IF (bcfl(1:3)=="all") THEN
          bciarr(:,1)=1_i4
          RETURN
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     determine the number of characters in bcfl and start a loop for
c     parsing individual characters.
c-----------------------------------------------------------------------
      char_loop: DO
        ichar=ichar+1
        SELECT CASE (bcfl(ichar:ichar))
        CASE (" ")
          CYCLE char_loop
        CASE ("s")
          ichar=ichar+1
          iind=iind+1
          IF (bcfl(ichar:ichar)=="d") bciarr(iind,1)=1_i4
        CASE ("3")
          ivec=ivec+1
          start3v(ivec)=iind+1
          ichar=ichar+2
          iind=iind+3
          IF (bcfl(ichar:ichar)=="n") THEN
            bciarr(iind-2:iind,2)=(/1_i4,0_i4,0_i4/)
          ELSE IF (bcfl(ichar:ichar)=="t") THEN
            bciarr(iind-2:iind,2)=(/0_i4,1_i4,1_i4/)
          ENDIF
        CASE DEFAULT
          CALL nim_stop("Bcflag_parse: "//bcfl(ichar:ichar)//
     $                  " not recognized.")
        END SELECT
        IF (iind==nqty) EXIT char_loop
      ENDDO char_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bcflag_parse
c-----------------------------------------------------------------------
c     subprogram 8. bcdir_set.
c     this routine sets the algebraic projection matrix
c     for projecting the components of the system that have essential,
c     i.e. Dirichlet conditions. 
c-----------------------------------------------------------------------
      SUBROUTINE bcdir_set(nqty,bciarr,excorner,isymm,norm,tang,bcpmat)
      USE seam_storage_mod
      USE matrix_type_mod

      INTEGER(i4), INTENT(IN) :: nqty,isymm
      INTEGER(i4), DIMENSION(nqty,2), INTENT(IN) :: bciarr
      LOGICAL, INTENT(IN) :: excorner
      REAL(r8), DIMENSION(2), INTENT(IN) :: norm,tang

      REAL(r8), DIMENSION(nqty,nqty), INTENT(OUT) :: bcpmat
      INTEGER(i4) :: ivec,ivst,iq
c-----------------------------------------------------------------------
c     intialize output.
c-----------------------------------------------------------------------
      bcpmat=0._r8
c-----------------------------------------------------------------------
c     catch the symmetry conditions for scalar systems.
c-----------------------------------------------------------------------
      IF (isymm/=0.AND.norm(2)==REAL(isymm,r8).AND.
     $    MINVAL(bciarr(:,1))==1_i4) RETURN
c-----------------------------------------------------------------------
c     set the directions for scalars and cases with "all".
c-----------------------------------------------------------------------
      DO iq=1,nqty
        bcpmat(iq,iq)=bciarr(iq,1)
      ENDDO
c-----------------------------------------------------------------------
c     set directions for 3-vectors.
c-----------------------------------------------------------------------
      DO ivec=1,SUM(MIN(start3v,1_i4))
        ivst=start3v(ivec)
        bcpmat(ivst+2,ivst+2)=bciarr(ivst+2,2)
        IF (excorner) THEN
          bcpmat(ivst,ivst)=1._r8
          bcpmat(ivst+1,ivst+1)=1._r8
          CYCLE
        ENDIF
        bcpmat(ivst:ivst+1,ivst)=
     $    norm*norm(1)*bciarr(ivst,2)+tang*tang(1)*bciarr(ivst+1,2)
        bcpmat(ivst:ivst+1,ivst+1)=
     $    norm*norm(2)*bciarr(ivst,2)+tang*tang(2)*bciarr(ivst+1,2)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE bcdir_set
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE boundary
