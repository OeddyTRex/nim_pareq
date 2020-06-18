c-----------------------------------------------------------------------
c     file rblock.f
c     subprograms for handling finite element computations on logically
c     rectangular grid-blocks.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  rblock_set.
c     2.  rblock_make_real_matrix.
c     3.  rblock_make_comp_matrix.
c     4.  rblock_get_real_rhs.
c     5.  rblock_get_comp_rhs.
c     6.  rblock_get_comp_rhs_q.
c     7.  rblock_basis_set.
c     8.  rblock_basis_dealloc.
c     9.  rblock_bicube_set.
c     10. rblock_real_qp_update.
c     11. rblock_comp_qp_update.
c     12. rblock_real_qp_alloc.
c     13. rblock_comp_qp_alloc.
c     14. rblock_real_qp_dealloc.
c     15. rblock_comp_qp_dealloc.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE rblock
      USE local
      USE rblock_type_mod
      IMPLICIT NONE

      INTERFACE rblock_make_matrix
        MODULE PROCEDURE rblock_make_real_matrix,rblock_make_comp_matrix
      END INTERFACE

      INTERFACE rblock_get_rhs
        MODULE PROCEDURE rblock_get_real_rhs,rblock_get_comp_rhs
      END INTERFACE

      INTERFACE rblock_qp_update
        MODULE PROCEDURE rblock_real_qp_update,rblock_comp_qp_update
      END INTERFACE

      INTERFACE rblock_qp_alloc
        MODULE PROCEDURE rblock_real_qp_alloc,rblock_comp_qp_alloc
      END INTERFACE

      INTERFACE rblock_qp_dealloc
        MODULE PROCEDURE rblock_real_qp_dealloc,rblock_comp_qp_dealloc
      END INTERFACE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. rblock_set.
c     set the locations and weights for quadratures.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_set(ngr,poly_degree,int_formula,rb)

      INTEGER(i4), INTENT(IN) :: ngr,poly_degree
      CHARACTER(8), INTENT(IN) :: int_formula
      TYPE(rblock_type), INTENT(INOUT) :: rb

      REAL(r8), DIMENSION(ngr+poly_degree-1) :: xg1d,wg1d
      INTEGER(i4) :: ix,iy,ig
c-----------------------------------------------------------------------
c     set number of quadrature points and weights according to input.
c     the number of points is now adjusted automatically with
c     poly_degree.
c-----------------------------------------------------------------------
      IF (int_formula=='gaussian') THEN
        CALL gauleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ELSE
        CALL lobleg(0._r8,1._r8,xg1d,wg1d,ngr+poly_degree-1_i4)
      ENDIF
      rb%ng=(ngr+poly_degree-1)**2
      ALLOCATE(rb%xg(rb%ng),rb%yg(rb%ng),rb%wg(rb%ng))
      ig=0
      DO iy=1,ngr+poly_degree-1
        DO ix=1,ngr+poly_degree-1
          ig=ig+1
          rb%xg(ig)=xg1d(ix)
          rb%yg(ig)=xg1d(iy)
          rb%wg(ig)=wg1d(ix)*wg1d(iy)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate rblock_set.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_set
c-----------------------------------------------------------------------
c     subprogram 2. rblock_make_real_matrix.
c     computes a linear response matrix for a supplied integrand
c     subprogram that fits the interface block at the beginning
c     of this subprogram.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_make_real_matrix(rb,mat_str,get_integrand,nqty)
      USE tblock_type_mod
      USE matrix_type_mod
      USE time

      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(rbl_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ityp,jtyp,idof,jdof,idofst,jdofst,
     $               jxoff,jyoff,ipol,ix,iy,ix0,iy0,mx,my
      REAL(r8) :: dx,dy
      REAL(r8) :: timestart_fe,timeend_fe,teststart,testend
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS :: mat
      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: integrand
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and preliminary computations.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      mx=rb%mx
      my=rb%my
c-----------------------------------------------------------------------
c     initialize matrix arrays.
c-----------------------------------------------------------------------
      DO jtyp=1,mat_str%nbtype
        DO ityp=1,mat_str%nbtype
          mat_str%mat(ityp,jtyp)%arr=0._r8
        ENDDO 
      ENDDO 
c-----------------------------------------------------------------------
c     flag the tblock as a dummy.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
c-----------------------------------------------------------------------
c     quadrature-point looping is now done inside the integrand routine,
c     and integrand is summed for each basis function, element by
c     element.
c
c     factors of Jacobian and quadrature weight are already in the basis
c     function arrays.
c-----------------------------------------------------------------------
      ALLOCATE(integrand(nqty,nqty,rb%mx*rb%my,mat_str%nbasis_el,
     $                   mat_str%nbasis_el))
      bigr=>rb%bigr
      dx=0
      dy=0
      CALL get_integrand(integrand,bigr,rb,dx,dy,tdum,0_i4)
c-----------------------------------------------------------------------
c     assemble and accumulate the contributions from each element.
c     the degree-of-freedom arrays are now used to replace
c     complicated looping with pre-computed information.
c-----------------------------------------------------------------------
      idofst=1
      DO ityp=1,mat_str%nbtype
        IF (mat_str%nq_type(ityp)==0) CYCLE
        jdofst=1
        DO jtyp=1,mat_str%nbtype
          IF (mat_str%nq_type(jtyp)==0) CYCLE
          mat=>mat_str%mat(jtyp,ityp)%arr
          DO ipol=1,mx*my
            iy0=(ipol-1)/mx+1
            ix0=ipol-mx*(iy0-1)
            DO idof=idofst,mat_str%den_type(ityp)
              iy=iy0+mat_str%dof_iy(idof)
              ix=ix0+mat_str%dof_ix(idof)
              DO jdof=jdofst,mat_str%den_type(jtyp)
                jxoff=mat_str%dof_ix(jdof)-mat_str%dof_ix(idof)
                jyoff=mat_str%dof_iy(jdof)-mat_str%dof_iy(idof)
                mat(mat_str%dof_iq(jdof),jxoff,jyoff,
     $              mat_str%dof_iq(idof),ix,iy)=
     $            mat(mat_str%dof_iq(jdof),jxoff,jyoff,
     $                mat_str%dof_iq(idof),ix,iy)+
     $            integrand(mat_str%dof_iv(jdof),mat_str%dof_iv(idof),
     $                 ipol,mat_str%dof_ib(jdof),mat_str%dof_ib(idof))
              ENDDO
            ENDDO
          ENDDO
          jdofst=mat_str%den_type(jtyp)+1
        ENDDO
        idofst=mat_str%den_type(ityp)+1
      ENDDO

      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_mat = time_mat + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_make_real_matrix
c-----------------------------------------------------------------------
c     subprogram 3. rblock_make_comp_matrix.
c     computes a linear response matrix for a supplied integrand
c     subprogram that fits the interface block at the beginning
c     of this subprogram.
c
c     complex version.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_make_comp_matrix(rb,mat_str,get_integrand,nqty)
      USE tblock_type_mod
      USE matrix_type_mod
      USE time

      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(rbl_comp_mat_type), INTENT(INOUT) :: mat_str
      INTEGER(i4), INTENT(IN) :: nqty

      INTEGER(i4) :: ityp,jtyp,idof,jdof,idofst,jdofst,
     $               jxoff,jyoff,ipol,ix,iy,ix0,iy0,mx,my
      REAL(r8) :: dx,dy
      REAL(r8) :: timestart_fe,timeend_fe,teststart,testend
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS :: mat
      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr

      TYPE (tblock_type) :: tdum

      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: integrand
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and preliminary computations.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      mx=rb%mx
      my=rb%my
c-----------------------------------------------------------------------
c     initialize matrix arrays.
c-----------------------------------------------------------------------
      DO jtyp=1,mat_str%nbtype
        DO ityp=1,mat_str%nbtype
          mat_str%mat(ityp,jtyp)%arr=0._r8
        ENDDO 
      ENDDO 
c-----------------------------------------------------------------------
c     flag the tblock as a dummy.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
c-----------------------------------------------------------------------
c     quadrature-point looping is now done inside the integrand routine,
c     and integrand is summed for each basis function, element by
c     element.
c
c     factors of Jacobian and quadrature weight are already in the basis
c     function arrays.
c-----------------------------------------------------------------------
      ALLOCATE(integrand(nqty,nqty,rb%mx*rb%my,mat_str%nbasis_el,
     $                   mat_str%nbasis_el))
      bigr=>rb%bigr
      dx=0
      dy=0
      CALL get_integrand(integrand,bigr,rb,dx,dy,tdum,0_i4)
c-----------------------------------------------------------------------
c     assemble and accumulate the contributions from each element.
c     the degree-of-freedom arrays are now used to replace
c     complicated looping with pre-computed information.
c-----------------------------------------------------------------------
      idofst=1
      DO ityp=1,mat_str%nbtype
        IF (mat_str%nq_type(ityp)==0) CYCLE
        jdofst=1
        DO jtyp=1,mat_str%nbtype
          IF (mat_str%nq_type(jtyp)==0) CYCLE
          mat=>mat_str%mat(jtyp,ityp)%arr
          DO ipol=1,mx*my
            iy0=(ipol-1)/mx+1
            ix0=ipol-mx*(iy0-1)
            DO idof=idofst,mat_str%den_type(ityp)
              iy=iy0+mat_str%dof_iy(idof)
              ix=ix0+mat_str%dof_ix(idof)
              DO jdof=jdofst,mat_str%den_type(jtyp)
                jxoff=mat_str%dof_ix(jdof)-mat_str%dof_ix(idof)
                jyoff=mat_str%dof_iy(jdof)-mat_str%dof_iy(idof)
                mat(mat_str%dof_iq(jdof),jxoff,jyoff,
     $              mat_str%dof_iq(idof),ix,iy)=
     $            mat(mat_str%dof_iq(jdof),jxoff,jyoff,
     $                mat_str%dof_iq(idof),ix,iy)+
     $            integrand(mat_str%dof_iv(jdof),mat_str%dof_iv(idof),
     $                 ipol,mat_str%dof_ib(jdof),mat_str%dof_ib(idof))
              ENDDO
            ENDDO
          ENDDO
          jdofst=mat_str%den_type(jtyp)+1
        ENDDO
        idofst=mat_str%den_type(ityp)+1
      ENDDO

      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_mat = time_mat + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_make_comp_matrix
c-----------------------------------------------------------------------
c     subprogram 4. rblock_get_real_rhs.
c     performs finite-element integrations for a rhs of an equation
c     producing real data, where the integrand is computed with a
c     supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_get_real_rhs(rb,rhs,get_integrand,nq)
      USE tblock_type_mod
      USE vector_type_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq
      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(vector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,ipol,mx,my,mxm1,mym1,iq,ig,iv,nv,ib
      INTEGER(i4) :: start_horz,start_vert,start_int,
     $               n_grid,n_horz,n_vert,n_int
      REAL(r8) :: dx,dy
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: integrand
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        REAL(r8), DIMENSION(:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and preliminary computations.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        rhs%arr(1:nq,:,:)=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        rhs%arrh(1:nq,:,:,:)=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        rhs%arrv(1:nq,:,:,:)=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        rhs%arri(1:nq,:,:,:)=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
c-----------------------------------------------------------------------
c     flag the tblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integrand(nq,mx*my,iv))
c-----------------------------------------------------------------------
c     looping over quadrature points is now within integrand routines.
c-----------------------------------------------------------------------
        dx=0
        dy=0
c-----------------------------------------------------------------------
c       evaluate the integrand all quadrature points.
c-----------------------------------------------------------------------
        bigr=>rb%bigr
        CALL get_integrand(integrand,bigr,rb,dx,dy,tdum,ig)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from each element
c       into the correct arrays.
c       grid vertex-centered bases first.
c
c       factors of Jacobian and quadrature weight are already in the
c       test function arrays.
c-----------------------------------------------------------------------
        IF (n_grid==4) THEN
          ipol=1
          DO iy=0,mym1
            DO ix=0,mxm1
              rhs%arr(1:nq,ix,iy)=rhs%arr(1:nq,ix,iy)
     $          +integrand(:,ipol,1)
              rhs%arr(1:nq,ix+1,iy)=rhs%arr(1:nq,ix+1,iy)
     $          +integrand(:,ipol,2)
              rhs%arr(1:nq,ix,iy+1)=rhs%arr(1:nq,ix,iy+1)
     $          +integrand(:,ipol,3)
              rhs%arr(1:nq,ix+1,iy+1)=rhs%arr(1:nq,ix+1,iy+1)
     $          +integrand(:,ipol,4)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          iv=start_horz
          DO ib=1,n_horz
            ipol=1
            DO iy=0,mym1
              DO ix=1,mx
                rhs%arrh(1:nq,ib,ix,iy)=rhs%arrh(1:nq,ib,ix,iy)
     $            +integrand(:,ipol,iv)
                rhs%arrh(1:nq,ib,ix,iy+1)=rhs%arrh(1:nq,ib,ix,iy+1)
     $            +integrand(:,ipol,iv+1)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          iv=start_vert
          DO ib=1,n_vert
            ipol=1
            DO iy=1,my
              DO ix=0,mxm1
                rhs%arrv(1:nq,ib,ix,iy)=rhs%arrv(1:nq,ib,ix,iy)
     $            +integrand(:,ipol,iv)
                rhs%arrv(1:nq,ib,ix+1,iy)=rhs%arrv(1:nq,ib,ix+1,iy)
     $            +integrand(:,ipol,iv+1)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+2
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       element interior-centered bases.
c-----------------------------------------------------------------------
        IF (n_int>0) THEN
          iv=start_int
          DO ib=1,n_int
            ipol=1
            DO iy=1,my
              DO ix=1,mx
                rhs%arri(1:nq,ib,ix,iy)=rhs%arri(1:nq,ib,ix,iy)
     $            +integrand(:,ipol,iv)
                ipol=ipol+1
              ENDDO
            ENDDO
            iv=iv+1
          ENDDO
        ENDIF

      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_real_rhs
c-----------------------------------------------------------------------
c     subprogram 5. rblock_get_comp_rhs.
c     performs finite-element integrations for a rhs of an equation
c     producing complex data, where the integrand is computed with a
c     supplied subroutine.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_get_comp_rhs(rb,rhs,get_integrand,nq,nfour)
      USE tblock_type_mod
      USE vector_type_mod
      USE time

      INTEGER(i4), INTENT(IN) :: nq,nfour
      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,mx,my,mxm1,mym1,iq,jf,ig,iv,nv,ib,ipol
      INTEGER(i4) :: start_horz,start_vert,start_int,
     $               n_grid,n_horz,n_vert,n_int
      REAL(r8) :: dx,dy
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integrand
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and preliminary computations.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        rhs%arr(1:nq,:,:,1:nfour)=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        rhs%arrh(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        rhs%arrv(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        rhs%arri(1:nq,:,:,:,1:nfour)=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int-1
c-----------------------------------------------------------------------
c     flag the tblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integrand(nq,mx*my,iv,nfour))
c-----------------------------------------------------------------------
c     looping over quadrature points is now within integrand routines.
c-----------------------------------------------------------------------
        dx=0
        dy=0
c-----------------------------------------------------------------------
c       evaluate the integrand all quadrature points.
c-----------------------------------------------------------------------
        bigr=>rb%bigr
        CALL get_integrand(integrand,bigr,rb,dx,dy,tdum,ig)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from each element
c       into the correct arrays.
c       grid vertex-centered bases first.
c
c       factors of Jacobian and quadrature weight are already in the
c       test function arrays.
c-----------------------------------------------------------------------
        IF (n_grid==4) THEN
          DO jf=1,nfour
            ipol=1
            DO iy=0,mym1
              DO ix=0,mxm1
                rhs%arr(1:nq,ix,iy,jf)=rhs%arr(1:nq,ix,iy,jf)
     $            +integrand(:,ipol,1,jf)
                rhs%arr(1:nq,ix+1,iy,jf)=rhs%arr(1:nq,ix+1,iy,jf)
     $            +integrand(:,ipol,2,jf)
                rhs%arr(1:nq,ix,iy+1,jf)=rhs%arr(1:nq,ix,iy+1,jf)
     $            +integrand(:,ipol,3,jf)
                rhs%arr(1:nq,ix+1,iy+1,jf)=rhs%arr(1:nq,ix+1,iy+1,jf)
     $            +integrand(:,ipol,4,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          DO jf=1,nfour
            iv=start_horz
            DO ib=1,n_horz
              ipol=1
              DO iy=0,mym1
                DO ix=1,mx
                  rhs%arrh(1:nq,ib,ix,iy,jf)=
     $              rhs%arrh(1:nq,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  rhs%arrh(1:nq,ib,ix,iy+1,jf)=
     $              rhs%arrh(1:nq,ib,ix,iy+1,jf)
     $              +integrand(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          DO jf=1,nfour
            iv=start_vert
            DO ib=1,n_vert
              ipol=1
              DO iy=1,my
                DO ix=0,mxm1
                  rhs%arrv(1:nq,ib,ix,iy,jf)=
     $              rhs%arrv(1:nq,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  rhs%arrv(1:nq,ib,ix+1,iy,jf)=
     $              rhs%arrv(1:nq,ib,ix+1,iy,jf)
     $              +integrand(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       element interior-centered bases.
c-----------------------------------------------------------------------
        IF (n_int>0) THEN
          DO jf=1,nfour
            iv=start_int
            DO ib=1,n_int
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  rhs%arri(1:nq,ib,ix,iy,jf)=rhs%arri(1:nq,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF

      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_comp_rhs
c-----------------------------------------------------------------------
c     subprogram 6. rblock_get_comp_rhs_q.
c     performs finite-element integrations for a rhs of an equation
c     producing complex data, where the integrand is computed with a
c     supplied subroutine.  this is the same as rblock_get_comp_rhs,
c     except that the rhs array is assumed to have the quantity and
c     Fourier component indices dimensioned correctly for this equation
c     for efficiency.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_get_comp_rhs_q(rb,rhs,get_integrand)
      USE tblock_type_mod
      USE vector_type_mod
      USE time

      TYPE(rblock_type), INTENT(INOUT), TARGET :: rb
      TYPE(cvector_type), INTENT(INOUT) :: rhs

      INTEGER(i4) :: ix,iy,mx,my,mxm1,mym1,iq,jf,ig,iv,nv,ib,ipol
      INTEGER(i4) :: start_horz,start_vert,start_int,start_disc,
     $               n_grid,n_horz,n_vert,n_int,nq,nfour,nqd,n_disc
      REAL(r8) :: dx,dy
      REAL(r8) :: timestart_fe,timeend_fe

      TYPE (tblock_type) :: tdum

      REAL(r8), DIMENSION(:,:), POINTER, CONTIGUOUS :: bigr
      COMPLEX(r8), DIMENSION(:,:,:,:), ALLOCATABLE :: integrand
c-----------------------------------------------------------------------
c     interface block for the integrand computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE get_integrand(integrand,bigr,rb,dx,dy,tb,ig)
        USE local
        USE rblock_type_mod
        USE tblock_type_mod
        COMPLEX(r8), DIMENSION(:,:,:,:), INTENT(OUT) :: integrand 
        REAL(r8), DIMENSION(:,:), INTENT(IN) :: bigr
        TYPE(rblock_type), INTENT(INOUT) :: rb
        REAL(r8), INTENT(IN) :: dx,dy
        TYPE(tblock_type), INTENT(INOUT) :: tb
        INTEGER(i4), INTENT(IN) :: ig
        END SUBROUTINE get_integrand
      END INTERFACE
c-----------------------------------------------------------------------
c     timer start and preliminary computations.
c-----------------------------------------------------------------------
      CALL timer(timestart_fe)
      mx=rb%mx
      my=rb%my
      mxm1=mx-1
      mym1=my-1
c-----------------------------------------------------------------------
c     examine the vector to determine what bases are used.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arr)) THEN
        n_grid=4
        nq=SIZE(rhs%arr,1)
        nfour=SIZE(rhs%arr,4)
        rhs%arr=0._r8
      ELSE
        n_grid=0
      ENDIF
      iv=n_grid+1
      IF (ASSOCIATED(rhs%arrh)) THEN
        start_horz=iv
        n_horz=SIZE(rhs%arrh,2)
        nq=SIZE(rhs%arrh,1)
        nfour=SIZE(rhs%arrh,5)
        rhs%arrh=0._r8
      ELSE
        n_horz=0
      ENDIF
      iv=iv+2*n_horz
      IF (ASSOCIATED(rhs%arrv)) THEN
        start_vert=iv
        n_vert=SIZE(rhs%arrv,2)
        nq=SIZE(rhs%arrv,1)
        nfour=SIZE(rhs%arrv,5)
        rhs%arrv=0._r8
      ELSE
        n_vert=0
      ENDIF
      iv=iv+2*n_vert
      IF (ASSOCIATED(rhs%arri)) THEN
        start_int=iv
        n_int=SIZE(rhs%arri,2)
        nq=SIZE(rhs%arri,1)
        nfour=SIZE(rhs%arri,5)
        rhs%arri=0._r8
      ELSE
        n_int=0
      ENDIF
      iv=iv+n_int
c-----------------------------------------------------------------------
c     a separate set of element-centered computations (for discontinuous
c     auxiliary fields) may be used.  if so, it is assumed that the
c     number of components at each node is no larger than nq.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(rhs%arrtmp)) THEN
        start_disc=iv
        n_disc=SIZE(rhs%arrtmp,2)
        nqd=SIZE(rhs%arrtmp,1)
        rhs%arrtmp=0._r8
      ELSE
        n_disc=0
        nqd=0
      ENDIF
      iv=iv+n_disc-1
c-----------------------------------------------------------------------
c     flag the tblock as a dummy and allocate int.
c-----------------------------------------------------------------------
      tdum%tgeom%mvert=-1
      ALLOCATE(integrand(nq,mx*my,iv,nfour))
c-----------------------------------------------------------------------
c     looping over quadrature points is now within integrand routines.
c-----------------------------------------------------------------------
        dx=0
        dy=0
c-----------------------------------------------------------------------
c       evaluate the integrand all quadrature points.
c-----------------------------------------------------------------------
        bigr=>rb%bigr
        CALL get_integrand(integrand,bigr,rb,dx,dy,tdum,ig)
c-----------------------------------------------------------------------
c       assemble and accumulate the contributions from each element
c       into the correct arrays.
c       grid vertex-centered bases first.
c
c       factors of Jacobian and quadrature weight are already in the
c       test function arrays.
c-----------------------------------------------------------------------
        IF (n_grid==4) THEN
          DO jf=1,nfour
            ipol=1
            DO iy=0,mym1
              DO ix=0,mxm1
                rhs%arr(:,ix,iy,jf)=rhs%arr(:,ix,iy,jf)
     $            +integrand(:,ipol,1,jf)
                rhs%arr(:,ix+1,iy,jf)=rhs%arr(:,ix+1,iy,jf)
     $            +integrand(:,ipol,2,jf)
                rhs%arr(:,ix,iy+1,jf)=rhs%arr(:,ix,iy+1,jf)
     $            +integrand(:,ipol,3,jf)
                rhs%arr(:,ix+1,iy+1,jf)=rhs%arr(:,ix+1,iy+1,jf)
     $            +integrand(:,ipol,4,jf)
                ipol=ipol+1
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       horizontal side-centered bases.
c-----------------------------------------------------------------------
        IF (n_horz>0) THEN
          DO jf=1,nfour
            iv=start_horz
            DO ib=1,n_horz
              ipol=1
              DO iy=0,mym1
                DO ix=1,mx
                  rhs%arrh(:,ib,ix,iy,jf)=rhs%arrh(:,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  rhs%arrh(:,ib,ix,iy+1,jf)=rhs%arrh(:,ib,ix,iy+1,jf)
     $              +integrand(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       vertical side-centered bases.
c-----------------------------------------------------------------------
        IF (n_vert>0) THEN
          DO jf=1,nfour
            iv=start_vert
            DO ib=1,n_vert
              ipol=1
              DO iy=1,my
                DO ix=0,mxm1
                  rhs%arrv(:,ib,ix,iy,jf)=rhs%arrv(:,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  rhs%arrv(:,ib,ix+1,iy,jf)=rhs%arrv(:,ib,ix+1,iy,jf)
     $              +integrand(:,ipol,iv+1,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+2
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       element interior-centered bases.
c-----------------------------------------------------------------------
        IF (n_int>0) THEN
          DO jf=1,nfour
            iv=start_int
            DO ib=1,n_int
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  rhs%arri(:,ib,ix,iy,jf)=rhs%arri(:,ib,ix,iy,jf)
     $              +integrand(:,ipol,iv,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       auxiliary computations for discontinuous auxiliary fields.
c-----------------------------------------------------------------------
        IF (n_disc>0) THEN
          DO jf=1,nfour
            iv=start_disc
            DO ib=1,n_disc
              ipol=1
              DO iy=1,my
                DO ix=1,mx
                  rhs%arrtmp(:,ib,ix,iy,jf)=rhs%arrtmp(:,ib,ix,iy,jf)
     $              +integrand(1:nqd,ipol,iv,jf)
                  ipol=ipol+1
                ENDDO
              ENDDO
              iv=iv+1
            ENDDO
          ENDDO
        ENDIF

      DEALLOCATE(integrand)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fe)
      time_rhs = time_rhs + timeend_fe - timestart_fe
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_get_comp_rhs_q
c-----------------------------------------------------------------------
c     subprogram 7. rblock_basis_set.
c     evaluate basis values and derivatives at quadrature points.
c
c     the input array, pd_arr, is a list of basis polynomial degree
c     values needed for continuous (nodal) expansions.
c
c     the input arrays, pdm_arr and pdmmin_arr, list of basis polynomial
c     degree values needed for discontinuous modal expansions.  these
c     expansions may be incomplete, so the minimum degree is provided
c     in the pdmmin_arr list. 
c-----------------------------------------------------------------------
      SUBROUTINE rblock_basis_set(rb,pd_arr,pdm_arr,pdmmin_arr,
     $                            pdmmax_arr,met_spl,geom)
      USE math_tran

      TYPE(rblock_type), INTENT(INOUT) :: rb
      INTEGER(i4), DIMENSION(:), INTENT(IN) :: pd_arr,pdm_arr,
     $             pdmmin_arr,pdmmax_arr
      CHARACTER(*), INTENT(IN) :: met_spl,geom

      REAL(r8), DIMENSION(:), ALLOCATABLE :: alpha,alphax,alphay
      REAL(r8), DIMENSION(2,rb%mx,rb%my) :: rz,drzdx,drzdy
      REAL(r8), DIMENSION(2,rb%ng,rb%mx*rb%my) :: drzdxq,drzdyq
      REAL(r8), DIMENSION(1,1,1) :: dc
      REAL(r8) :: dx,dy
      INTEGER(i4), DIMENSION(SIZE(pd_arr)) :: use_pd
      INTEGER(i4), DIMENSION(SIZE(pdm_arr)) :: use_pdm
      INTEGER(i4) :: iv,nv,ig,mx,my,ib,ix,iy,ipol,num_bases,iset,
     $               num_basesm
c-----------------------------------------------------------------------
c     determine the number of unique basis sets are needed according
c     to the pd_arr array, which lists the polynomial degree of
c     different fields.
c-----------------------------------------------------------------------
      use_pd=1_i4
      DO ib=2,SIZE(pd_arr)
        IF (MINVAL(ABS(pd_arr(1:ib-1)-pd_arr(ib)))==0) use_pd(ib)=0
      ENDDO
      num_bases=SUM(use_pd)
      ALLOCATE(rb%base_pd(num_bases))
      NULLIFY(rb%base_disc)  !  not used
c-----------------------------------------------------------------------
c     do the same for modal discontinuous bases.  negative values
c     for the polynomial degree indicated a field that is not used.
c     a unique basis is determined by the pd value and by the minimum
c     and maximum degrees for the limited representation (all three
c     values have to match for two expansions to have the same basis).
c-----------------------------------------------------------------------
      use_pdm=1_i4
      IF (pdm_arr(1)<0) use_pdm(1)=0
      DO ib=2,SIZE(pdm_arr)
        IF (MINVAL(ABS(pdm_arr(1:ib-1)-pdm_arr(ib)))==0.AND.
     $      MINVAL(ABS(pdmmin_arr(1:ib-1)-pdmmin_arr(ib)))==0.AND.
     $      MINVAL(ABS(pdmmax_arr(1:ib-1)-pdmmax_arr(ib)))==0)
     $    use_pdm(ib)=0
        IF (pdm_arr(ib)<0) use_pdm(ib)=0
      ENDDO
      num_basesm=SUM(use_pdm)
      ALLOCATE(rb%base_modal(num_basesm))
c-----------------------------------------------------------------------
c     allocate arrays used for saving the basis function values, basis
c     function gradient values, and grid derivatives at each of the
c     Gaussian quadrature points.  the loop provides storage for all
c     unique basis sets.
c
c     the quadrature points are now indexed first, and (ix,iy) indices
c     are combined.
c-----------------------------------------------------------------------
      mx=rb%mx
      my=rb%my
      iset=0
      DO ib=1,SIZE(pd_arr)
        IF (use_pd(ib)==0) CYCLE
        iset=iset+1
        rb%base_pd(iset)%poly_deg_basis=pd_arr(ib)
        nv=(pd_arr(ib)+1)**2
        ALLOCATE(rb%base_pd(iset)%alpha(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpdrc(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%alpham(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_pd(iset)%dalpmdrc(rb%ng,mx*my,nv))
      ENDDO
      iset=0
      DO ib=1,SIZE(pdm_arr)
        IF (use_pdm(ib)==0) CYCLE
        iset=iset+1
        rb%base_modal(iset)%poly_deg_basis=pdm_arr(ib)
        rb%base_modal(iset)%poly_degmin_basis=pdmmin_arr(ib)
        rb%base_modal(iset)%poly_degmax_basis=pdmmax_arr(ib)
        nv=(pdmmax_arr(ib)-pdmmin_arr(ib)+1)*
     $     (2*pdm_arr(ib)+pdmmin_arr(ib)-pdmmax_arr(ib)+1)
        ALLOCATE(rb%base_modal(iset)%alpha(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpdrc(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%alpham(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdr(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdz(rb%ng,mx*my,nv))
        ALLOCATE(rb%base_modal(iset)%dalpmdrc(rb%ng,mx*my,nv))
      ENDDO
      ALLOCATE(rb%dxdr(rb%ng,mx*my))
      ALLOCATE(rb%dxdz(rb%ng,mx*my))
      ALLOCATE(rb%dydr(rb%ng,mx*my))
      ALLOCATE(rb%dydz(rb%ng,mx*my))
      ALLOCATE(rb%bigr(rb%ng,mx*my))
      ALLOCATE(rb%wjac(rb%ng,mx*my))
      ALLOCATE(rb%jac2d(rb%ng,mx*my))
c-----------------------------------------------------------------------
c     evaluate and save the basis function weights at the quadrature
c     points.  evaluate the grid derivatives to construct the
c     gradients of the basis functions.  finally, save the quadrature
c     weight times the coordinate-mapping Jacobian for efficiency during
c     finite element computations.
c
c     for toroidal geometry dalpdrc holds d(alpha)/dr + alpha/r, and
c     for linear geometry it's just d(alpha)/dx.
c-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        SELECT CASE(met_spl)
        CASE('pcnst')
          rz=0.25*(rb%rz%fs(:,0:mx-1,0:my-1)
     $            +rb%rz%fs(:,1:mx  ,0:my-1)
     $            +rb%rz%fs(:,0:mx-1,1:my  )
     $            +rb%rz%fs(:,1:mx  ,1:my  ))
          drzdx=0.5*(rb%rz%fs(:,1:mx  ,1:my  )
     $              -rb%rz%fs(:,0:mx-1,1:my  )
     $              +rb%rz%fs(:,1:mx  ,0:my-1)
     $              -rb%rz%fs(:,0:mx-1,0:my-1))
          drzdy=0.5*(rb%rz%fs(:,1:mx  ,1:my  )
     $              -rb%rz%fs(:,1:mx  ,0:my-1)
     $              +rb%rz%fs(:,0:mx-1,1:my  )
     $              -rb%rz%fs(:,0:mx-1,0:my-1))
        CASE('liner','linear','bilinear')
          rz=(1-dx)*(1-dy)*rb%rz%fs(:,0:mx-1,0:my-1)
     $      +   dx *(1-dy)*rb%rz%fs(:,1:mx  ,0:my-1)
     $      +(1-dx)*   dy *rb%rz%fs(:,0:mx-1,1:my  )
     $      +   dx *   dy *rb%rz%fs(:,1:mx  ,1:my  )
          drzdx=   dy *(rb%rz%fs(:,1:mx  ,1:my  )
     $                 -rb%rz%fs(:,0:mx-1,1:my  ))
     $         +(1-dy)*(rb%rz%fs(:,1:mx  ,0:my-1)
     $                 -rb%rz%fs(:,0:mx-1,0:my-1))
          drzdy=   dx *(rb%rz%fs(:,1:mx  ,1:my  )
     $                 -rb%rz%fs(:,1:mx  ,0:my-1))
     $         +(1-dx)*(rb%rz%fs(:,0:mx-1,1:my  )
     $                 -rb%rz%fs(:,0:mx-1,0:my-1))
        CASE('iso','lagrz')
          CALL lagr_quad_all_eval(rb%rz,dx,dy,rz,drzdx,drzdy,1_i4)
        CASE DEFAULT
          CALL nim_stop('Rblock_basis_set: '//TRIM(met_spl)//
     $                  ' is not a valid option for met_spl.')
        END SELECT
        
        ipol=1
        DO iy=1,rb%my
          DO ix=1,rb%mx
            drzdxq(:,ig,ipol)=drzdx(:,ix,iy)
            drzdyq(:,ig,ipol)=drzdy(:,ix,iy)
            rb%bigr(ig,ipol)=rz(1,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     invert the Jacobian matrix and save partial derivatives.
c-----------------------------------------------------------------------
      CALL math_grid('all',drzdxq,drzdyq,rb%jac2d(:,:),
     $               rb%dxdr(:,:),rb%dxdz(:,:),
     $               rb%dydr(:,:),rb%dydz(:,:))
      rb%jac2d=MAX(rb%jac2d,0._r8)
      rb%wjac(:,:)=rb%jac2d(:,:)
      IF (geom=='tor') THEN
        rb%wjac=rb%wjac*rb%bigr
      ELSE
        rb%bigr=1._r8
      ENDIF
c-----------------------------------------------------------------------
c     multiply the jacobian by the quadrature weight.
c-----------------------------------------------------------------------
      DO ig=1,rb%ng
        rb%wjac(ig,:)=rb%wg(ig)*rb%wjac(ig,:)
      ENDDO
c-----------------------------------------------------------------------
c     load the alpha arrays for all quadrature points with the index
c     reordering.  this is needed for each unique basis set.
c
c     basis/test functions used in the matrix integrand routines are
c     multipied by the square root of (Jacobian*quadrature-weight).
c
c     test functions used in rhs integrand routines are multiplied by
c     Jacobian*quadrature-weight (no square root).
c-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        DO iset=1,num_bases
          nv=(rb%base_pd(iset)%poly_deg_basis+1)**2
          ALLOCATE(alpha(nv),alphax(nv),alphay(nv))
          CALL lagr_quad_bases(dx,dy,alpha,alphax,alphay,1_i4)

          DO iv=1,nv
            rb%base_pd(iset)%alpha (ig,:,iv)=alpha(iv)
            rb%base_pd(iset)%dalpdr(ig,:,iv)=rb%dxdr(ig,:)*alphax(iv)
     $                                      +rb%dydr(ig,:)*alphay(iv)
            rb%base_pd(iset)%dalpdz(ig,:,iv)=rb%dxdz(ig,:)*alphax(iv)
     $                                      +rb%dydz(ig,:)*alphay(iv)
            IF (geom=='tor') THEN
              rb%base_pd(iset)%dalpdrc(ig,:,iv)=
     $          rb%base_pd(iset)%dalpdr(ig,:,iv)+alpha(iv)/rb%bigr(ig,:)
            ELSE
              rb%base_pd(iset)%dalpdrc(ig,:,iv)=
     $          rb%base_pd(iset)%dalpdr(ig,:,iv)
            ENDIF

            rb%base_pd(iset)%alpham(ig,:,iv)=
     $        rb%base_pd(iset)%alpha(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdr(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdr(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdz(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdz(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%dalpmdrc(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdrc(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_pd(iset)%alpha(ig,:,iv)=
     $        rb%base_pd(iset)%alpha(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdr(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdr(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdz(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdz(ig,:,iv)*rb%wjac(ig,:)
            rb%base_pd(iset)%dalpdrc(ig,:,iv)=
     $        rb%base_pd(iset)%dalpdrc(ig,:,iv)*rb%wjac(ig,:)
          ENDDO
          DEALLOCATE(alpha,alphax,alphay)
        ENDDO

        DO iset=1,num_basesm
          nv=(rb%base_modal(iset)%poly_degmax_basis
     $       -rb%base_modal(iset)%poly_degmin_basis+1)*
     $       (2*rb%base_modal(iset)%poly_deg_basis+1
     $       -rb%base_modal(iset)%poly_degmax_basis
     $       +rb%base_modal(iset)%poly_degmin_basis)
          ALLOCATE(alpha(nv),alphax(nv),alphay(nv))
          CALL modal_disc_bases(dx,dy,alpha,alphax,alphay,1_i4,
     $                          rb%base_modal(iset)%poly_deg_basis,
     $                          rb%base_modal(iset)%poly_degmin_basis,
     $                          rb%base_modal(iset)%poly_degmax_basis)

          DO iv=1,nv
            rb%base_modal(iset)%alpha (ig,:,iv)=alpha(iv)
            rb%base_modal(iset)%dalpdr(ig,:,iv)=rb%dxdr(ig,:)*alphax(iv)
     $                                         +rb%dydr(ig,:)*alphay(iv)
            rb%base_modal(iset)%dalpdz(ig,:,iv)=rb%dxdz(ig,:)*alphax(iv)
     $                                         +rb%dydz(ig,:)*alphay(iv)
            IF (geom=='tor') THEN
              rb%base_modal(iset)%dalpdrc(ig,:,iv)=
     $          rb%base_modal(iset)%dalpdr(ig,:,iv)+
     $          alpha(iv)/rb%bigr(ig,:)
            ELSE
              rb%base_modal(iset)%dalpdrc(ig,:,iv)=
     $          rb%base_modal(iset)%dalpdr(ig,:,iv)
            ENDIF

            rb%base_modal(iset)%alpham(ig,:,iv)=
     $        rb%base_modal(iset)%alpha(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdr(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdr(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdz(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdz(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%dalpmdrc(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdrc(ig,:,iv)*SQRT(rb%wjac(ig,:))
            rb%base_modal(iset)%alpha(ig,:,iv)=
     $        rb%base_modal(iset)%alpha(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdr(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdr(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdz(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdz(ig,:,iv)*rb%wjac(ig,:)
            rb%base_modal(iset)%dalpdrc(ig,:,iv)=
     $        rb%base_modal(iset)%dalpdrc(ig,:,iv)*rb%wjac(ig,:)
          ENDDO
          DEALLOCATE(alpha,alphax,alphay)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_basis_set
c-----------------------------------------------------------------------
c     subprogram 8. rblock_basis_dealloc.
c     deallocates arrays for rblock bases.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_basis_dealloc(rb)

      TYPE(rblock_type), INTENT(INOUT) :: rb
      INTEGER(i4) :: iset

      DO iset=1,SIZE(rb%base_pd)
        DEALLOCATE(rb%base_pd(iset)%alpha)
        DEALLOCATE(rb%base_pd(iset)%dalpdr)
        DEALLOCATE(rb%base_pd(iset)%dalpdz)
        DEALLOCATE(rb%base_pd(iset)%dalpdrc)
        DEALLOCATE(rb%base_pd(iset)%alpham)
        DEALLOCATE(rb%base_pd(iset)%dalpmdr)
        DEALLOCATE(rb%base_pd(iset)%dalpmdz)
        DEALLOCATE(rb%base_pd(iset)%dalpmdrc)
      ENDDO
      DEALLOCATE(rb%base_pd)
      DO iset=1,SIZE(rb%base_modal)
        DEALLOCATE(rb%base_modal(iset)%alpha)
        DEALLOCATE(rb%base_modal(iset)%dalpdr)
        DEALLOCATE(rb%base_modal(iset)%dalpdz)
        DEALLOCATE(rb%base_modal(iset)%dalpdrc)
        DEALLOCATE(rb%base_modal(iset)%alpham)
        DEALLOCATE(rb%base_modal(iset)%dalpmdr)
        DEALLOCATE(rb%base_modal(iset)%dalpmdz)
        DEALLOCATE(rb%base_modal(iset)%dalpmdrc)
      ENDDO
      DEALLOCATE(rb%base_modal)
      DEALLOCATE(rb%dxdr)
      DEALLOCATE(rb%dxdz)
      DEALLOCATE(rb%dydr)
      DEALLOCATE(rb%dydz)
      DEALLOCATE(rb%bigr)
      DEALLOCATE(rb%wjac)
      DEALLOCATE(rb%jac2d)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_basis_dealloc
c-----------------------------------------------------------------------
c     subprogram 9. rblock_bicube_set.
c     evaluate bicubic spline data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_bicube_set(bc,qbc,rb,d_order)

      TYPE(bicube_type), INTENT(INOUT) :: bc
      TYPE(rb_real_qp_type), INTENT(OUT) :: qbc
      TYPE(rblock_type), INTENT(INOUT) :: rb
      CHARACTER(*), INTENT(IN) :: d_order
      
      INTEGER(i4) :: ig,bcmode,mxb,myb,iq,ix,iy,ipol
      REAL(r8), DIMENSION(bc%nqty,rb%mx,rb%my) :: f,fx,fy
      REAL(r8), DIMENSION(1,1,1) :: db
      REAL(r8) :: dx,dy
c-----------------------------------------------------------------------
c     allocate quadrature storage space as needed.
c-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      ALLOCATE(qbc%qpf(bc%nqty,rb%ng,mxb*myb))
      IF (d_order/='values') THEN
        ALLOCATE(qbc%qpfr(bc%nqty,rb%ng,mxb*myb))
        ALLOCATE(qbc%qpfz(bc%nqty,rb%ng,mxb*myb))
      ELSE
        ALLOCATE(qbc%qpfr(1,1,1))
        ALLOCATE(qbc%qpfz(1,1,1))
      ENDIF
c-----------------------------------------------------------------------
c     set evaluation mode.
c-----------------------------------------------------------------------
      SELECT CASE(d_order)
      CASE('values')
        bcmode=0
      CASE('1st derivs')
        bcmode=1
      CASE DEFAULT
        CALL nim_stop
     $    ('Unrecognized derivative order in rblock_bicube_set.')
      END SELECT
c-----------------------------------------------------------------------
c     loop over quadrature points, evaluate, and save.
c-----------------------------------------------------------------------
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL bicube_all_eval(bc,dx,dy,f,fx,fy,db,db,db,bcmode)
        ipol=1
        DO iy=1,myb
          DO ix=1,mxb
            qbc%qpf(:,ig,ipol)=f(:,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
        IF (bcmode>0) THEN
          ipol=1
          DO iy=1,myb
            DO ix=1,mxb
              qbc%qpfr(:,ig,ipol)=rb%dxdr(ig,ipol)*fx(:,ix,iy)
     $                           +rb%dydr(ig,ipol)*fy(:,ix,iy)
              qbc%qpfz(:,ig,ipol)=rb%dxdz(ig,ipol)*fx(:,ix,iy)
     $                           +rb%dydz(ig,ipol)*fy(:,ix,iy)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     deallocate bicube evaluation matrix.
c-----------------------------------------------------------------------
      DEALLOCATE(bc%cmats)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_bicube_set
c-----------------------------------------------------------------------
c     subprogram 10. rblock_real_qp_update.
c     evaluate 2D lagrange_quad data at the gaussian quadrature points
c     for this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_update(laq,qlaq,rb)

      TYPE(lagr_quad_2D_type), INTENT(IN) :: laq
      TYPE(rb_real_qp_type), INTENT(INOUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      
      INTEGER(i4) :: ig,mxb,myb,iq,ix,iy,ipol
      REAL(r8), DIMENSION(laq%nqty,rb%mx,rb%my) :: f,fx,fy
      REAL(r8) :: dx,dy
c-----------------------------------------------------------------------
c     loop over quadrature points, evaluate, and save.
c-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL lagr_quad_all_eval(laq,dx,dy,f,fx,fy,1_i4)
        ipol=1
        DO iy=1,myb
          DO ix=1,mxb
            qlaq%qpf (:,ig,ipol)=f(:,ix,iy)
            qlaq%qpfr(:,ig,ipol)=rb%dxdr(ig,ipol)*fx(:,ix,iy)
     $                          +rb%dydr(ig,ipol)*fy(:,ix,iy)
            qlaq%qpfz(:,ig,ipol)=rb%dxdz(ig,ipol)*fx(:,ix,iy)
     $                          +rb%dydz(ig,ipol)*fy(:,ix,iy)
            ipol=ipol+1
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_update
c-----------------------------------------------------------------------
c     subprogram 11. rblock_comp_qp_update.
c     evaluate 3D lagrange_quad data at the gaussian quadrature points
c     for this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_update(laq,qlaq,rb)

      TYPE(lagr_quad_type), INTENT(IN) :: laq
      TYPE(rb_comp_qp_type), INTENT(INOUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      
      INTEGER(i4) :: ig,mxb,myb,iq,im,ix,iy,ipol
      COMPLEX(r8), DIMENSION(laq%nqty,rb%mx,rb%my,laq%nfour) :: f,fx,fy
      REAL(r8) :: dx,dy
c-----------------------------------------------------------------------
c     loop over quadrature points, evaluate, and save.
c-----------------------------------------------------------------------
      mxb=rb%mx
      myb=rb%my
      DO ig=1,rb%ng
        dx=rb%xg(ig)
        dy=rb%yg(ig)
        CALL lagr_quad_all_eval(laq,dx,dy,f,fx,fy,1_i4)
        DO im=1,laq%nfour
          ipol=1
          DO iy=1,myb
            DO ix=1,mxb
              qlaq%qpf (:,ig,ipol,im)=f(:,ix,iy,im)
              qlaq%qpfr(:,ig,ipol,im)=rb%dxdr(ig,ipol)*fx(:,ix,iy,im)
     $                               +rb%dydr(ig,ipol)*fy(:,ix,iy,im)
              qlaq%qpfz(:,ig,ipol,im)=rb%dxdz(ig,ipol)*fx(:,ix,iy,im)
     $                               +rb%dydz(ig,ipol)*fy(:,ix,iy,im)
              ipol=ipol+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_update
c-----------------------------------------------------------------------
c     subprogram 12. rblock_real_qp_alloc.
c     allocate space for 2D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_alloc(qlaq,rb,nqty)

      TYPE(rb_real_qp_type), INTENT(OUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty
      
      ALLOCATE(qlaq%qpf (nqty,rb%ng,rb%mx*rb%my)); qlaq%qpf=0._r8
      ALLOCATE(qlaq%qpfr(nqty,rb%ng,rb%mx*rb%my)); qlaq%qpfr=0._r8
      ALLOCATE(qlaq%qpfz(nqty,rb%ng,rb%mx*rb%my)); qlaq%qpfz=0._r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_alloc
c-----------------------------------------------------------------------
c     subprogram 13. rblock_comp_qp_alloc.
c     allocate space for 3D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_alloc(qlaq,rb,nqty,nfour)

      TYPE(rb_comp_qp_type), INTENT(OUT) :: qlaq
      TYPE(rblock_type), INTENT(IN) :: rb
      INTEGER(i4), INTENT(IN) :: nqty,nfour
      
      ALLOCATE(qlaq%qpf(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpf=0._r8
      ALLOCATE(qlaq%qpfr(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpfr=0._r8
      ALLOCATE(qlaq%qpfz(nqty,rb%ng,rb%mx*rb%my,nfour)); qlaq%qpfz=0._r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_alloc
c-----------------------------------------------------------------------
c     subprogram 14. rblock_real_qp_dealloc.
c     deallocate space for 2D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_real_qp_dealloc(qlaq)

      TYPE(rb_real_qp_type), INTENT(INOUT) :: qlaq

      DEALLOCATE(qlaq%qpf)
      IF (ALLOCATED(qlaq%qpfr)) DEALLOCATE(qlaq%qpfr)
      IF (ALLOCATED(qlaq%qpfz)) DEALLOCATE(qlaq%qpfz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_real_qp_dealloc
c-----------------------------------------------------------------------
c     subprogram 15. rblock_comp_qp_dealloc.
c     deallocate space for 3D data at the gaussian quadrature points for
c     this block.
c-----------------------------------------------------------------------
      SUBROUTINE rblock_comp_qp_dealloc(qlaq)

      TYPE(rb_comp_qp_type), INTENT(INOUT) :: qlaq
      
      DEALLOCATE(qlaq%qpf)
      IF (ALLOCATED(qlaq%qpfr)) DEALLOCATE(qlaq%qpfr)
      IF (ALLOCATED(qlaq%qpfz)) DEALLOCATE(qlaq%qpfz)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE rblock_comp_qp_dealloc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE rblock
