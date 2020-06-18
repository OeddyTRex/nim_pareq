c-----------------------------------------------------------------------
c     file matrix_type_mod.f
c     contains a module that defines structures for saving matrices.
c     (structures for matrix factors are now in factor_type_mod.)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  matrix_type_mod.
c     1.  matrix_rbl_real_alloc
c     2.  matrix_rbl_real_dealloc
c     3.  matrix_rbl_comp_alloc
c     4.  matrix_rbl_comp_dealloc
c     5.  matrix_rbl_make_real
c     6.  matrix_rbl_dof_init
c     7.  matrix_tbl_real_alloc
c     8.  matrix_tbl_real_dealloc
c     9.  matrix_tbl_comp_alloc
c     10. matrix_tbl_comp_dealloc
c     11. matrix_tbl_make_real
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE matrix_type_mod
      USE local
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     types used for defining matrix structures for a single Fourier
c     component.
c-----------------------------------------------------------------------
      TYPE :: global_matrix_type
        TYPE(rbl_mat_type), DIMENSION(:), POINTER :: rbl_mat
        TYPE(tbl_mat_type), DIMENSION(:), POINTER :: tbl_mat
        INTEGER(i4) :: fcomp,foff
        INTEGER(i4) :: nqty,nqdis
        CHARACTER(1), DIMENSION(:), POINTER :: vcomp
        LOGICAL :: eliminated
        LOGICAL :: symmetric
        REAL(r8) :: diag_scale
        CHARACTER(16) :: essential_cond
      END TYPE global_matrix_type
      TYPE :: complex_matrix_type
        TYPE(rbl_comp_mat_type), DIMENSION(:), POINTER :: rbl_mat
        TYPE(tbl_comp_mat_type), DIMENSION(:), POINTER :: tbl_mat
        INTEGER(i4) :: fcomp,foff
        INTEGER(i4) :: nqty,nqdis
        CHARACTER(1), DIMENSION(:), POINTER :: vcomp
        LOGICAL :: eliminated
        LOGICAL :: hermitian
        REAL(r8) :: diag_scale
        CHARACTER(16) :: essential_cond
      END TYPE complex_matrix_type
c-----------------------------------------------------------------------
c     types used for holding matrix elements in rblocks.
c     nbasis_el is the number of nonzero basis functions per element.
c     nbasis_cont and nbasis_disc are the numbers of continuous and
c     discontinuous bases, respectively.
c
c     the mat arrays are dimensioned nbtypeXnbtype to cover
c     the matrix coupling among the different basis types.
c
c     ndof_el is the number of degrees of freedom per element, and
c     the pointer arrays dof_ib, dof_iv, dof_ix, dof_iy, and dof_iq hold
c     information that helps transfer integrals from element-based
c     storage to structured rblock storage.  for each degree of freedom,
c
c       dof_ib = element basis-function index
c       dof_iv = element (physical-field) vector-component index
c       dof_ix = block ix offset relative to element (=-1 or 0)
c       dof_iy = block iy offset relative to element (=-1 or 0)
c       dof_iq = block quantity index
c
c     also, for each basis type, den_type hold the ending DOF index.
c-----------------------------------------------------------------------
      TYPE :: rbl_mat_type
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbasis_cont
        INTEGER(i4) :: nbasis_disc
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: ndof_el
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(i4), DIMENSION(:), POINTER :: dof_ib,dof_iv,dof_ix,
     $               dof_iy,dof_iq
        INTEGER(i4), DIMENSION(:), POINTER :: den_type
        TYPE(arr_6d_type), DIMENSION(:,:), POINTER :: mat
        INTEGER(i4), POINTER :: mem_id
      END TYPE rbl_mat_type
      TYPE :: rbl_comp_mat_type
        INTEGER(i4) :: nbasis_el
        INTEGER(i4) :: nbasis_cont
        INTEGER(i4) :: nbasis_disc
        INTEGER(i4) :: nbtype
        INTEGER(i4) :: ndof_el
        INTEGER(i4) :: mx,my
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        INTEGER(i4), DIMENSION(:), POINTER :: nb_type,nq_type
        INTEGER(i4), DIMENSION(:), POINTER :: dof_ib,dof_iv,dof_ix,
     $               dof_iy,dof_iq
        INTEGER(i4), DIMENSION(:), POINTER :: den_type
        TYPE(comp_arr_6d_type), DIMENSION(:,:), POINTER :: mat
        INTEGER(i4), POINTER :: mem_id
      END TYPE rbl_comp_mat_type
c-----------------------------------------------------------------------
c     types used for holding rblock matrix elements.
c-----------------------------------------------------------------------
      TYPE :: arr_6d_type
        REAL(r8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: arr
      END TYPE arr_6d_type
      TYPE :: comp_arr_6d_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: arr
      END TYPE comp_arr_6d_type
c-----------------------------------------------------------------------
c     types used for holding matrix elements in tblocks.
c-----------------------------------------------------------------------
      TYPE :: tbl_mat_type
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        TYPE(matrix_element_type3), DIMENSION(:), POINTER :: lmat
      END TYPE tbl_mat_type
      TYPE :: tbl_comp_mat_type
        INTEGER(i4), DIMENSION(:), POINTER :: ix0,iy0
        TYPE(comp_matrix_element_type3), DIMENSION(:), POINTER :: lmat
      END TYPE tbl_comp_mat_type
c-----------------------------------------------------------------------
c     types used for holding tblock matrix elements.  it was formerly
c     part of tblock_type_mod.
c-----------------------------------------------------------------------
      TYPE :: matrix_element_type3
        REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: element
        INTEGER(i4), DIMENSION(:), POINTER :: from_vert
      END TYPE matrix_element_type3
      TYPE :: comp_matrix_element_type3
        COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: element
        INTEGER(i4), DIMENSION(:), POINTER :: from_vert
      END TYPE comp_matrix_element_type3

c-----------------------------------------------------------------------
c     interfaces for the allocation and deallocation routines.
c-----------------------------------------------------------------------
      INTERFACE matrix_rbl_alloc
        MODULE PROCEDURE matrix_rbl_real_alloc,matrix_rbl_comp_alloc
      END INTERFACE

      INTERFACE matrix_rbl_dealloc
        MODULE PROCEDURE matrix_rbl_real_dealloc,matrix_rbl_comp_dealloc
      END INTERFACE

      INTERFACE matrix_tbl_alloc
        MODULE PROCEDURE matrix_tbl_real_alloc,matrix_tbl_comp_alloc
      END INTERFACE

      INTERFACE matrix_tbl_dealloc
        MODULE PROCEDURE matrix_tbl_real_dealloc,matrix_tbl_comp_dealloc
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. matrix_rbl_real_alloc
c     allocate arrays needed for a real rblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_real_alloc(rmat,mx,my,nqc,pdc,nqd,nbd)

      TYPE(rbl_mat_type), INTENT(OUT) :: rmat
      INTEGER(i4), INTENT(IN) :: mx,my,nqc,pdc,nqd,nbd

      INTEGER(i4) :: id,jd,ixst,iyst,jxst,jyst,x0off,x1off,y0off,y1off,
     $               sz
c-----------------------------------------------------------------------
c     the mat array allows multiple basis function types (grid vertex,
c     horizontal element side, vertical element side, and interior-
c     centered).
c
c     the structures accommodate equations with nqc continuous fields
c     of polynomial degree pdc and nqd discontinuous fields with
c     nbd basis functions.  the storage for 'element interiors'
c     is used for both the interior coefficients of continuous fields
c     and for all coefficients of discontinuous fields.
c
c     if poly_degree=1, there is only one basis type.  otherwise, there
c     are 4.
c-----------------------------------------------------------------------
      rmat%nbasis_el=MIN(nqc,1_i4)*(pdc+1)**2+MIN(nqd,1_i4)*nbd
      rmat%nbasis_cont=MIN(nqc,1_i4)*(pdc+1)**2
      rmat%nbasis_disc=MIN(nqd,1_i4)*nbd
      rmat%ndof_el=nqc*(pdc+1)**2+nqd*nbd
      IF (pdc==1.AND.nqd==0) THEN
        rmat%nbtype=1
      ELSE
        rmat%nbtype=4
      ENDIF
      rmat%mx=mx
      rmat%my=my
      ALLOCATE(rmat%mat(rmat%nbtype,rmat%nbtype))
      ALLOCATE(rmat%ix0(rmat%nbtype))
      ALLOCATE(rmat%iy0(rmat%nbtype))
      ALLOCATE(rmat%nb_type(rmat%nbtype))
      ALLOCATE(rmat%nq_type(rmat%nbtype))
      ALLOCATE(rmat%den_type(rmat%nbtype))
      ALLOCATE(rmat%dof_ib(rmat%ndof_el))
      ALLOCATE(rmat%dof_iv(rmat%ndof_el))
      ALLOCATE(rmat%dof_ix(rmat%ndof_el))
      ALLOCATE(rmat%dof_iy(rmat%ndof_el))
      ALLOCATE(rmat%dof_iq(rmat%ndof_el))
c-----------------------------------------------------------------------
c     logical indices for each of the basis types.
c       index 1 is grid vertex-centered.
c       index 2 is horizontal side-centered.
c       index 3 is vertical side-centered.
c       index 4 is interior-centered and all discontinuous bases.
c
c     this version has the 6D matrix array indices defined
c     (col_comp,col_x_off,col_y_off,row_comp,row_x_index,row_y_index),
c     where comp is vector component and basis for types with multiple
c     bases (vector component varying faster) and off is the offset
c     from the row index.
c-----------------------------------------------------------------------
      DO id=1,rmat%nbtype
        SELECT CASE(id)
        CASE(1)
          rmat%ix0(id)=0
          rmat%iy0(id)=0
          rmat%nb_type(id)=MIN(nqc,1_i4)
          rmat%nq_type(id)=nqc
          rmat%den_type(id)=4*nqc
        CASE(2)
          rmat%ix0(id)=1
          rmat%iy0(id)=0
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          rmat%nq_type(id)=nqc*(pdc-1)
          rmat%den_type(id)=rmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(3)
          rmat%ix0(id)=0
          rmat%iy0(id)=1
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          rmat%nq_type(id)=nqc*(pdc-1)
          rmat%den_type(id)=rmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(4)
          rmat%ix0(id)=1
          rmat%iy0(id)=1
          rmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)**2+
     $                     MIN(nqd,1_i4)*nbd
          rmat%nq_type(id)=nqc*(pdc-1)**2+nqd*nbd
          rmat%den_type(id)=rmat%den_type(id-1)+rmat%nq_type(id)
        END SELECT
      ENDDO
      DO id=1,rmat%nbtype
        ixst=rmat%ix0(id)
        iyst=rmat%iy0(id)
        DO jd=1,rmat%nbtype
          jxst=rmat%ix0(jd)
          jyst=rmat%iy0(jd)
          x0off=jxst-1
          x1off=1-ixst
          y0off=jyst-1
          y1off=1-iyst
          ALLOCATE(rmat%mat(jd,id)%
     $      arr(rmat%nq_type(jd),x0off:x1off,y0off:y1off,
     $          rmat%nq_type(id),ixst:mx,iyst:my))
        ENDDO
      ENDDO
      CALL matrix_rbl_dof_init(rmat%dof_ib,rmat%dof_iv,rmat%dof_ix,
     $                         rmat%dof_iy,rmat%dof_iq,nqc,pdc,nqd,nbd)
c-----------------------------------------------------------------------
c     register this object.
c-----------------------------------------------------------------------
      NULLIFY(rmat%mem_id)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_real_alloc
c-----------------------------------------------------------------------
c     subprogram 2. matrix_rbl_real_dealloc
c     deallocate arrays needed for a real rblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_real_dealloc(rmat)

      TYPE(rbl_mat_type), INTENT(INOUT) :: rmat

      INTEGER(i4) :: id,jd

c-----------------------------------------------------------------------
c     loop over different basis combinations and deallocte.
c-----------------------------------------------------------------------
      DO jd=1,SIZE(rmat%mat,2)
        DO id=1,SIZE(rmat%mat,1)
          DEALLOCATE(rmat%mat(id,jd)%arr)
        ENDDO
      ENDDO
      DEALLOCATE(rmat%mat,rmat%ix0,rmat%iy0)
      DEALLOCATE(rmat%nb_type,rmat%nq_type)
      DEALLOCATE(rmat%den_type,rmat%dof_ib,rmat%dof_iv,
     $           rmat%dof_ix,rmat%dof_iy,rmat%dof_iq)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_real_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. matrix_rbl_comp_alloc
c     allocate arrays needed for a comp rblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_comp_alloc(cmat,mx,my,nqc,pdc,nqd,nbd)

      TYPE(rbl_comp_mat_type), INTENT(OUT) :: cmat
      INTEGER(i4), INTENT(IN) :: mx,my,nqc,pdc,nqd,nbd

      INTEGER(i4) :: id,jd,ixst,iyst,jxst,jyst,x0off,x1off,y0off,y1off,
     $               sz
c-----------------------------------------------------------------------
c     the mat array allows multiple basis function types (grid vertex,
c     horizontal element side, vertical element side, and interior-
c     centered).
c
c     the structures accommodate equations with nqc continuous fields
c     of polynomial degree pdc and nqd discontinuous fields with
c     nbd basis functions.  the storage for 'element interiors'
c     is used for both the interior coefficients of continuous fields
c     and for all coefficients of discontinuous fields.
c
c     if poly_degree=1, there is only one basis type.  otherwise, there
c     are 4.
c-----------------------------------------------------------------------
      cmat%nbasis_el=MIN(nqc,1_i4)*(pdc+1)**2+MIN(nqd,1_i4)*nbd
      cmat%nbasis_cont=MIN(nqc,1_i4)*(pdc+1)**2
      cmat%nbasis_disc=MIN(nqd,1_i4)*nbd
      cmat%ndof_el=nqc*(pdc+1)**2+nqd*nbd
      cmat%mx=mx
      cmat%my=my
      IF (pdc==1.AND.nqd==0) THEN
        cmat%nbtype=1
      ELSE
        cmat%nbtype=4
      ENDIF
      ALLOCATE(cmat%mat(cmat%nbtype,cmat%nbtype))
      ALLOCATE(cmat%ix0(cmat%nbtype))
      ALLOCATE(cmat%iy0(cmat%nbtype))
      ALLOCATE(cmat%nb_type(cmat%nbtype))
      ALLOCATE(cmat%nq_type(cmat%nbtype))
      ALLOCATE(cmat%den_type(cmat%nbtype))
      ALLOCATE(cmat%dof_ib(cmat%ndof_el))
      ALLOCATE(cmat%dof_iv(cmat%ndof_el))
      ALLOCATE(cmat%dof_ix(cmat%ndof_el))
      ALLOCATE(cmat%dof_iy(cmat%ndof_el))
      ALLOCATE(cmat%dof_iq(cmat%ndof_el))
c-----------------------------------------------------------------------
c     logical indices for each of the basis types.
c       index 1 is grid vertex-centered.
c       index 2 is horizontal side-centered.
c       index 3 is vertical side-centered.
c       index 4 is interior-centered and all discontinuous bases.
c
c     this version has the 6D matrix array indices defined
c     (col_comp,col_x_off,col_y_off,row_comp,row_x_index,row_y_index),
c     where comp is vector component and basis for types with multiple
c     bases (vector component varying faster) and off is the offset
c     from the row index.
c-----------------------------------------------------------------------
      DO id=1,cmat%nbtype
        SELECT CASE(id)
        CASE(1)
          cmat%ix0(id)=0
          cmat%iy0(id)=0
          cmat%nb_type(id)=MIN(nqc,1_i4)
          cmat%nq_type(id)=nqc
          cmat%den_type(id)=4*nqc
        CASE(2)
          cmat%ix0(id)=1
          cmat%iy0(id)=0
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          cmat%nq_type(id)=nqc*(pdc-1)
          cmat%den_type(id)=cmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(3)
          cmat%ix0(id)=0
          cmat%iy0(id)=1
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)
          cmat%nq_type(id)=nqc*(pdc-1)
          cmat%den_type(id)=cmat%den_type(id-1)+2*nqc*(pdc-1)
        CASE(4)
          cmat%ix0(id)=1
          cmat%iy0(id)=1
          cmat%nb_type(id)=MIN(nqc,1_i4)*(pdc-1)**2+
     $                     MIN(nqd,1_i4)*nbd
          cmat%nq_type(id)=nqc*(pdc-1)**2+nqd*nbd
          cmat%den_type(id)=cmat%den_type(id-1)+cmat%nq_type(id)
        END SELECT
      ENDDO
      DO id=1,cmat%nbtype
        ixst=cmat%ix0(id)
        iyst=cmat%iy0(id)
        DO jd=1,cmat%nbtype
          jxst=cmat%ix0(jd)
          jyst=cmat%iy0(jd)
          x0off=jxst-1
          x1off=1-ixst
          y0off=jyst-1
          y1off=1-iyst
          ALLOCATE(cmat%mat(jd,id)%
     $      arr(cmat%nq_type(jd),x0off:x1off,y0off:y1off,
     $          cmat%nq_type(id),ixst:mx,iyst:my))
        ENDDO
      ENDDO
      CALL matrix_rbl_dof_init(cmat%dof_ib,cmat%dof_iv,cmat%dof_ix,
     $                         cmat%dof_iy,cmat%dof_iq,nqc,pdc,nqd,nbd)
c-----------------------------------------------------------------------
c     register this object.
c-----------------------------------------------------------------------
      NULLIFY(cmat%mem_id)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_comp_alloc
c-----------------------------------------------------------------------
c     subprogram 4. matrix_rbl_comp_dealloc
c     deallocate arrays needed for a comp rblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_comp_dealloc(cmat)

      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: cmat

      INTEGER(i4) :: id,jd

c-----------------------------------------------------------------------
c     loop over different basis combinations and deallocte.
c-----------------------------------------------------------------------
      DO jd=1,SIZE(cmat%mat,2)
        DO id=1,SIZE(cmat%mat,1)
          DEALLOCATE(cmat%mat(id,jd)%arr)
        ENDDO
      ENDDO
      DEALLOCATE(cmat%mat,cmat%ix0,cmat%iy0)
      DEALLOCATE(cmat%nb_type,cmat%nq_type)
      DEALLOCATE(cmat%den_type,cmat%dof_ib,cmat%dof_iv,
     $           cmat%dof_ix,cmat%dof_iy,cmat%dof_iq)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_comp_dealloc
c-----------------------------------------------------------------------
c     subprogram 5. matrix_rbl_make_real
c     set matrix elements of a comp rblock matrix to their real parts.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_make_real(cmat)

      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: cmat

      INTEGER(i4) :: id,jd

c-----------------------------------------------------------------------
c     loop over different basis combinations.
c-----------------------------------------------------------------------
      DO jd=1,SIZE(cmat%mat,2)
        DO id=1,SIZE(cmat%mat,1)
          cmat%mat(id,jd)%arr=REAL(cmat%mat(id,jd)%arr,r8)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_rbl_make_real
c-----------------------------------------------------------------------
c     subprogram 6. matrix_rbl_dof_init
c     initialize arrays that help translate element-based data.
c
c     see matrix_rbl_real_alloc for definitions of nqc, pdc, nqd,
c     and nbd.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_rbl_dof_init(dib,div,dix,diy,diq,
     $                               nqc,pdc,nqd,nbd)

      INTEGER(i4), DIMENSION(:), INTENT(OUT) :: dib,div,dix,diy,diq
      INTEGER(i4), INTENT(IN) :: nqc,pdc,nqd,nbd

      INTEGER(i4) :: ix,iy,iq,ii,id
c-----------------------------------------------------------------------
c     fill the dof arrays that are used to generate indices
c     when transferring information into rmat arrays.
c-----------------------------------------------------------------------
      ii=1
      DO iy=0,1        !     4 continuous grid-centered bases
        DO ix=0,1
          DO iq=1,nqc
            dib(ii)=ix+2*iy+1
            div(ii)=iq
            dix(ii)=ix-1
            diy(ii)=iy-1
            diq(ii)=iq
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,pdc-2    !     continuous horizontal side bases
        DO iy=0,1
          DO iq=1,nqc
            dib(ii)=iy+2*id+5
            div(ii)=iq
            dix(ii)=0
            diy(ii)=iy-1
            diq(ii)=iq+nqc*id
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,pdc-2    !     continuous vertical side bases
        DO ix=0,1
          DO iq=1,nqc
            dib(ii)=ix+2*id+5+2*(pdc-1)
            div(ii)=iq
            dix(ii)=ix-1
            diy(ii)=0
            diq(ii)=iq+nqc*id
            ii=ii+1
          ENDDO
        ENDDO
      ENDDO
      DO id=0,(pdc-1)**2-1  !  interior bases for continuous expansions
        DO iq=1,nqc
          dib(ii)=id+1+4*pdc
          div(ii)=iq
          dix(ii)=0
          diy(ii)=0
          diq(ii)=iq+nqc*id
          ii=ii+1
        ENDDO
      ENDDO
      DO id=0,nbd-1    !  bases for discontinuous expansions
        DO iq=1,nqd
          dib(ii)=id+1+MIN(nqc,1_i4)*(pdc+1)**2  ! continue interior #s
          div(ii)=iq
          dix(ii)=0
          diy(ii)=0
          diq(ii)=iq+nqd*id+nqc*(pdc-1)**2
          ii=ii+1
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE matrix_rbl_dof_init
c-----------------------------------------------------------------------
c     subprogram 7. matrix_tbl_real_alloc
c     allocate arrays needed for a real tblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_real_alloc(tmat,tg,nq)
      USE tri_linear

      TYPE(tbl_mat_type), INTENT(OUT) :: tmat
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: iv,nn
c-----------------------------------------------------------------------
c-PRE allocate matrix elements for linear bases only.
c-----------------------------------------------------------------------
      ALLOCATE(tmat%lmat(0:tg%mvert))
      DO iv=0,tg%mvert
        nn=SIZE(tg%neighbor(iv)%vertex)-1
        ALLOCATE(tmat%lmat(iv)%element(nq,nq,0:nn))
        ALLOCATE(tmat%lmat(iv)%from_vert(0:nn))
        tmat%lmat(iv)%from_vert=tg%neighbor(iv)%vertex
      ENDDO
      ALLOCATE(tmat%ix0(1))
      ALLOCATE(tmat%iy0(1))
      tmat%ix0=0
      tmat%iy0=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_real_alloc
c-----------------------------------------------------------------------
c     subprogram 8. matrix_tbl_real_dealloc
c     deallocate arrays needed for a real tblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_real_dealloc(tmat)

      TYPE(tbl_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
c-----------------------------------------------------------------------
c-PRE deallocate linear bases only.
c-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        DEALLOCATE(tmat%lmat(iv)%element)
        DEALLOCATE(tmat%lmat(iv)%from_vert)
      ENDDO
      DEALLOCATE(tmat%lmat,tmat%ix0,tmat%iy0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_real_dealloc
c-----------------------------------------------------------------------
c     subprogram 9. matrix_tbl_comp_alloc
c     allocate arrays needed for a comp tblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_comp_alloc(tmat,tg,nq)
      USE tri_linear

      TYPE(tbl_comp_mat_type), INTENT(OUT) :: tmat
      TYPE(tri_linear_geom_type), INTENT(IN) :: tg
      INTEGER(i4), INTENT(IN) :: nq

      INTEGER(i4) :: iv,nn
c-----------------------------------------------------------------------
c-PRE allocate matrix elements for linear bases only.
c-----------------------------------------------------------------------
      ALLOCATE(tmat%lmat(0:tg%mvert))
      DO iv=0,tg%mvert
        nn=SIZE(tg%neighbor(iv)%vertex)-1
        ALLOCATE(tmat%lmat(iv)%element(nq,nq,0:nn))
        ALLOCATE(tmat%lmat(iv)%from_vert(0:nn))
        tmat%lmat(iv)%from_vert=tg%neighbor(iv)%vertex
      ENDDO
      ALLOCATE(tmat%ix0(1))
      ALLOCATE(tmat%iy0(1))
      tmat%ix0=0
      tmat%iy0=0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_comp_alloc
c-----------------------------------------------------------------------
c     subprogram 10. matrix_tbl_comp_dealloc
c     deallocate arrays needed for a comp tblock matrix.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_comp_dealloc(tmat)

      TYPE(tbl_comp_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
c-----------------------------------------------------------------------
c-PRE deallocate linear bases only.
c-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        DEALLOCATE(tmat%lmat(iv)%element)
        DEALLOCATE(tmat%lmat(iv)%from_vert)
      ENDDO
      DEALLOCATE(tmat%lmat,tmat%ix0,tmat%iy0)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_comp_dealloc
c-----------------------------------------------------------------------
c     subprogram 11. matrix_tbl_make_real
c     set matrix elements of a comp tblock matrix to their real parts.
c-----------------------------------------------------------------------
      SUBROUTINE matrix_tbl_make_real(tmat)

      TYPE(tbl_comp_mat_type), INTENT(INOUT) :: tmat

      INTEGER(i4) :: iv
c-----------------------------------------------------------------------
c-PRE using linear bases only.
c-----------------------------------------------------------------------
      DO iv=0,SIZE(tmat%lmat)-1
        tmat%lmat(iv)%element=REAL(tmat%lmat(iv)%element,r8)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE matrix_tbl_make_real

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matrix_type_mod
