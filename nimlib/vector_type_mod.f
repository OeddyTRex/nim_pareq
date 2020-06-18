c-----------------------------------------------------------------------
c     file vector_type_mod.f
c     contains a module that defines structures for the linear algebra
c     vectors.
c     
c     this has been broken into separate modules for different data
c     types to facilitate compilation on some machines.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module vector_defn_type_mod
c     module rvector_type_mod
c     1. vector_rtype_alloc.
c     2. vector_rtype_dealloc.
c     3. vector_assign_rsc.
c     4. vector_assign_csc.
c     5. vector_assign_int.
c     6. vector_assign_vec.
c     6.1. vector_assignq_vec.
c     7. vector_assign_cvec.
c     7.1. vector_assign_cvec3.
c     7.2. vector_assign_cvec2.
c     8. vector_ptassign_bc.
c     9. vector_ptassign_laq2.
c     9.1. vector_ptassign_modq2.
c     10. vector_ptassign_tl2.
c     11. vector_add_vec.
c     12. vector_mult_rsc.
c     13. vector_mult_int.
c     module cvector_type_mod
c     14. vector_ctype_alloc.
c     15. vector_ctype_dealloc.
c     16. cvector_assign_rsc.
c     17. cvector_assign_csc.
c     18. cvector_assign_int.
c     19. cvector_assign_cvec.
c     19.1. cvector_assignq_cvec.
c     20. cvector_assign_vec.
c     21. cvector_assign_cvec2.
c     22. cvector_ptassign_laq.
c     22.1. cvector_ptassign_modq.
c     23. cvector_ptassign_tl.
c     24. cvector_add_cvec.
c     25. cvector_addc_cvec.
c     26. cvector_mult_rsc.
c     27. cvector_mult_csc.
c     28. cvector_mult_int.
c     28.1. cvector_real_comp.
c     module cvector_2D_type_mod
c     29. vector_2D_ctype_alloc.
c     30. vector_2D_ctype_dealloc.
c     31. cvector_2D_assign_rsc.
c     32. cvector_2D_assign_csc.
c     33. cvector_2D_assign_int.
c     34. cvector_2D_assign_cvec2.
c     34.1. cvector_2D_assign_vec.
c     35. cvector_2D_add_cvec2.
c     36. cvector_2D_addc_cvec2.
c     37. cvector_2D_assign_cvec.
c     38. cvector_2D_mult_rsc.
c     39. cvector_2D_mult_csc.
c     40. cvector_2D_mult_int.
c     40.1. cvector_2D_pack_cvec.
c     40.2. cvector_2D_pack_cvec2.
c     40.3. cvector_2D_unpack_cvec.
c     40.4. cvector_2D_unpack_add_cvec.
c     40.5. cvector_2D_unpack_cvec2.
c     module rvector_3D_type_mod
c     41. vector_3D_rtype_alloc.
c     42. vector_3D_rtype_dealloc.
c     43. vector_3D_assign_rsc.
c     44. vector_3D_assign_csc.
c     45. vector_3D_assign_int.
c     46. vector_3D_assign_vec.
c     47. vector_3D_assign_vec3.
c     48. vector_3D_add_vec.
c     49. vector_3D_add_vec3.
c     50. vector_3D_mult_rsc.
c     51. vector_3D_mult_int.
c     module vector_type_mod
c-----------------------------------------------------------------------
c     module for type definitions.
c-----------------------------------------------------------------------
      MODULE vector_defn_type_mod
      USE local
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     the vector_type is set-up for 2D arrays of real vector quantities
c     with element side and interior centerings, as well as grid
c     vertices.
c
c     the arrtmp arrays are extra space for temporary fields
c     and are not affected by assignment and algebraic operations.
c-----------------------------------------------------------------------
      TYPE :: vector_type
        REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: arr
        REAL(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $            arrh,arrv,arri,arrtmp
      END TYPE vector_type
c-----------------------------------------------------------------------
c     the cvector_type is set-up for 3D arrays of complex vector
c     quantities with element side and interior centerings, as well
c     as grid vertices.
c-----------------------------------------------------------------------
      TYPE :: cvector_type
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: arr
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER, CONTIGUOUS ::
     $               arrh,arrv,arri,arrtmp
      END TYPE cvector_type
c-----------------------------------------------------------------------
c     the cvector_2D_type is 2D complex vector array used for interim
c     computations.
c-----------------------------------------------------------------------
      TYPE :: cvector_2D_type
        COMPLEX(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: arr
        COMPLEX(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS ::
     $               arrh,arrv,arri,arrtmp
      END TYPE cvector_2D_type
c-----------------------------------------------------------------------
c     the vector_3D_type is set-up for 3D arrays of real vector
c     quantities with the last index running over indices for the
c     periodic coordinate.
c-----------------------------------------------------------------------
      TYPE :: vector_3D_type
        REAL(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: arr
        REAL(r8), DIMENSION(:,:,:,:,:), POINTER, CONTIGUOUS ::
     $            arrh,arrv,arri,arrd,arrtmp
      END TYPE vector_3D_type
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE vector_defn_type_mod


c-----------------------------------------------------------------------
c     module for 2D-in-space arrays of real vectors quantities.
c-----------------------------------------------------------------------
      MODULE rvector_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. vector_rtype_alloc.
c     allocates space for a real vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_rtype_alloc(rvt,poly_degree,mx,my,nqty,nbt,nqt,
     $                              alloc_int)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(vector_type), INTENT(OUT) :: rvt
      LOGICAL, INTENT(IN), OPTIONAL :: alloc_int
c-----------------------------------------------------------------------
c     allocate space according to the basis functions needed.
c     if poly_degree is non-positive, storage for discontinuous-field
c     coefficients is allocated.
c
c     interior-basis arrays are allocated by default unless prevented
c     by the optional alloc_int being set to false.
c
c-PRE triangles will need something here, too.
c-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,mx,my))
        NULLIFY(rvt%arr,rvt%arrh,rvt%arrv)
      CASE(1)  !  linear elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my))
        NULLIFY(rvt%arri,rvt%arrh,rvt%arrv)
      CASE(2:)  !  higher-order elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my))
        ALLOCATE(rvt%arrh(nqty,poly_degree-1,1:mx,0:my))
        ALLOCATE(rvt%arrv(nqty,poly_degree-1,0:mx,1:my))
        IF (.NOT.PRESENT(alloc_int)) THEN
          ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
        ELSE IF (alloc_int) THEN
          ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
        ELSE
          NULLIFY(rvt%arri)
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     if the optional input for temporary arrays is present, allocate
c     the arrtmp arrays; otherwise, nullify them.  note that nbt is
c     the number of temporary bases, and nqt is the number of quantities
c     at each basis.
c-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(rvt%arrtmp(nqt,nbt,1:mx,1:my))
      ELSE
        NULLIFY(rvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_rtype_alloc
c-----------------------------------------------------------------------
c     subprogram 2. vector_rtype_dealloc.
c     deallocates space for a real vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_rtype_dealloc(rvt)

      TYPE(vector_type), INTENT(INOUT) :: rvt

      IF (ASSOCIATED(rvt%arr)) THEN
        DEALLOCATE(rvt%arr)
        NULLIFY(rvt%arr)
      ENDIF
      IF (ASSOCIATED(rvt%arrh)) THEN
        DEALLOCATE(rvt%arrh)
        NULLIFY(rvt%arrh)
      ENDIF
      IF (ASSOCIATED(rvt%arrv)) THEN
        DEALLOCATE(rvt%arrv)
        NULLIFY(rvt%arrv)
      ENDIF
      IF (ASSOCIATED(rvt%arri)) THEN
        DEALLOCATE(rvt%arri)
        NULLIFY(rvt%arri)
      ENDIF
      IF (ASSOCIATED(rvt%arrtmp)) THEN
        DEALLOCATE(rvt%arrtmp)
        NULLIFY(rvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_rtype_dealloc
c-----------------------------------------------------------------------
c     subprogram 3. vector_assign_rsc.
c     assign a real scalar value to a vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_rsc(vec,rscalar)

      TYPE(vector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop
     $    ('Vector_assign_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 4. vector_assign_csc.
c     assign a complex scalar value to a vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_csc(vec,cscalar)

      TYPE(vector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop
     $    ('Vector_assign_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_csc
c-----------------------------------------------------------------------
c     subprogram 5. vector_assign_int.
c     assign a integer value to a vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_int(vec,int)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop
     $    ('Vector_assign_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_int
c-----------------------------------------------------------------------
c     subprogram 6. vector_assign_vec.
c     set one vector structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_vec(vec1,vec2)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Vector_assign_vec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_vec
c-----------------------------------------------------------------------
c     subprogram 6.1 vector_assignq_vec.
c     set a quantity range of one real vector structure equal to
c     that of another.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assignq_vec(vec1,vec2,nqty,nstart1,nstart2)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty
      INTEGER(i4), INTENT(IN), OPTIONAL :: nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (PRESENT(nstart1)) THEN
        nst1=nstart1
      ELSE
        nst1=1
      ENDIF
      IF (PRESENT(nstart2)) THEN
        nst2=nstart2
      ELSE
        nst2=1
      ENDIF
      nend1=nst1+nqty-1
      nend2=nst2+nqty-1
      IF (ASSOCIATED(vec1%arr)) THEN
        vec1%arr(nst1:nend1,:,:)=vec2%arr(nst2:nend2,:,:)
        IF (ASSOCIATED(vec1%arrh))
     $    vec1%arrh(nst1:nend1,:,:,:)=vec2%arrh(nst2:nend2,:,:,:)
        IF (ASSOCIATED(vec1%arrv))
     $    vec1%arrv(nst1:nend1,:,:,:)=vec2%arrv(nst2:nend2,:,:,:)
        IF (ASSOCIATED(vec1%arri))
     $    vec1%arri(nst1:nend1,:,:,:)=vec2%arri(nst2:nend2,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(nst1:nend1,:,:,:)=
     $    vec2%arrtmp(nst2:nend2,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri)) THEN
        vec1%arri(nst1:nend1,:,:,:)=vec2%arri(nst2:nend2,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Vector_assignq_vec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assignq_vec
c-----------------------------------------------------------------------
c     subprogram 7. vector_assign_cvec.
c     set one real vector structure equal to part of a complex vector
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_cvec(vec1,vec2,r_i,fcomp,nqty,nstart1,
     $                              nstart2)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.  nstart1 and nstart2 can be used to set starting
c     indices.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(nst1:nend1,:,:)=vec2%arr(nst2:nend2,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $        vec2%arrh(nst2:nend2,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $        vec2%arrv(nst2:nend2,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $        vec2%arri(nst2:nend2,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $        vec2%arrtmp(nst2:nend2,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arr(nst1:nend1,:,:)=
     $        AIMAG(vec2%arr(nst2:nend2,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrh(nst2:nend2,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrv(nst2:nend2,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arri(nst2:nend2,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrtmp(nst2:nend2,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(nst1:nend1,:,:,:)=
     $        vec2%arri(nst2:nend2,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arri(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arri(nst2:nend2,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Vector_assign_cvec: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c     the r12mi3 and i12r3 flags are used in several of the
c     nimrod management routines.  r12mi3 means transfer the real
c     1 & 2 vector components and minus the third imaginary  comp.
c     i12r3 means transfer the imaginary 1 & 2 and the real third.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr=vec2%arr(:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh=vec2%arrh(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv=vec2%arrv(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri=vec2%arri(:,:,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp=vec2%arrtmp(:,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arr=AIMAG(vec2%arr(:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh=AIMAG(vec2%arrh(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv=AIMAG(vec2%arrv(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri=AIMAG(vec2%arri(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp=AIMAG(vec2%arrtmp(:,:,:,:,fcomp))
          CASE ('r12mi3')
            vec1%arr(1:2,:,:)=vec2%arr(1:2,:,:,fcomp)
            vec1%arr(3,:,:)=-AIMAG(vec2%arr(3,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=vec2%arrh(1:2,:,:,:,fcomp)
              vec1%arrh(3,:,:,:)=-AIMAG(vec2%arrh(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=vec2%arrv(1:2,:,:,:,fcomp)
              vec1%arrv(3,:,:,:)=-AIMAG(vec2%arrv(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:,fcomp)
              vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=vec2%arrtmp(1:2,:,:,:,fcomp)
              vec1%arrtmp(3,:,:,:)=-AIMAG(vec2%arrtmp(3,:,:,:,fcomp))
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:)=AIMAG(vec2%arr(1:2,:,:,fcomp))
            vec1%arr(3,:,:)=vec2%arr(3,:,:,fcomp)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=AIMAG(vec2%arrh(1:2,:,:,:,fcomp))
              vec1%arrh(3,:,:,:)=vec2%arrh(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=AIMAG(vec2%arrv(1:2,:,:,:,fcomp))
              vec1%arrv(3,:,:,:)=vec2%arrv(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:,fcomp))
              vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:,fcomp)
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=AIMAG(vec2%arrtmp(1:2,:,:,:,fcomp))
              vec1%arrtmp(3,:,:,:)=vec2%arrtmp(3,:,:,:,fcomp)
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri=vec2%arri(:,:,:,:,fcomp)
          CASE ('imag','IMAG')
            vec1%arri=AIMAG(vec2%arri(:,:,:,:,fcomp))
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:,fcomp)
            vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:,fcomp))
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:,fcomp))
            vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:,fcomp)
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Vector_assign_cvec: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_cvec
c-----------------------------------------------------------------------
c     subprogram 7.1. vector_assign_vec3.
c     set a vector structure equal to a 3D vector structure at the
c     specified index of the periodic coordinate.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_vec3(vec1,vec2,ip)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: ip
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c--------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:,ip)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) 
     $    vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:,ip)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:,ip)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:,ip)
      ELSE
        CALL nim_stop
     $    ('Vector_assign_vec3: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_vec3
c-----------------------------------------------------------------------
c     subprogram 7.2. vector_assign_cvec2.
c     set one real vector structure equal to part of a complex 2D vector
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_assign_cvec2(vec1,vec2,r_i,nqty,nstart1,nstart2)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.  nstart1 and nstart2 can be used to set starting
c     indices.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(nst1:nend1,:,:)=vec2%arr(nst2:nend2,:,:)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $        vec2%arrh(nst2:nend2,:,:,:)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $        vec2%arrv(nst2:nend2,:,:,:)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $        vec2%arri(nst2:nend2,:,:,:)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $        vec2%arrtmp(nst2:nend2,:,:,:)
          CASE ('imag','IMAG')
            vec1%arr(nst1:nend1,:,:)=
     $        AIMAG(vec2%arr(nst2:nend2,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrh(nst2:nend2,:,:,:))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrv(nst2:nend2,:,:,:))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arri(nst2:nend2,:,:,:))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arrtmp(nst2:nend2,:,:,:))
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec2: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(nst1:nend1,:,:,:)=
     $        vec2%arri(nst2:nend2,:,:,:)
          CASE ('imag','IMAG')
            vec1%arri(nst1:nend1,:,:,:)=
     $        AIMAG(vec2%arri(nst2:nend2,:,:,:))
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec2: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Vector_assign_cvec2: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c     the r12mi3 and i12r3 flags are used in several of the
c     nimrod management routines.  r12mi3 means transfer the real
c     1 & 2 vector components and minus the third imaginary  comp.
c     i12r3 means transfer the imaginary 1 & 2 and the real third.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr=vec2%arr(:,:,:)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh=vec2%arrh(:,:,:,:)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv=vec2%arrv(:,:,:,:)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri=vec2%arri(:,:,:,:)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp=vec2%arrtmp(:,:,:,:)
          CASE ('imag','IMAG')
            vec1%arr=AIMAG(vec2%arr(:,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh=AIMAG(vec2%arrh(:,:,:,:))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv=AIMAG(vec2%arrv(:,:,:,:))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri=AIMAG(vec2%arri(:,:,:,:))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp=AIMAG(vec2%arrtmp(:,:,:,:))
          CASE ('r12mi3')
            vec1%arr(1:2,:,:)=vec2%arr(1:2,:,:)
            vec1%arr(3,:,:)=-AIMAG(vec2%arr(3,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=vec2%arrh(1:2,:,:,:)
              vec1%arrh(3,:,:,:)=-AIMAG(vec2%arrh(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=vec2%arrv(1:2,:,:,:)
              vec1%arrv(3,:,:,:)=-AIMAG(vec2%arrv(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:)
              vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=vec2%arrtmp(1:2,:,:,:)
              vec1%arrtmp(3,:,:,:)=-AIMAG(vec2%arrtmp(3,:,:,:))
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:)=AIMAG(vec2%arr(1:2,:,:))
            vec1%arr(3,:,:)=vec2%arr(3,:,:)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=AIMAG(vec2%arrh(1:2,:,:,:))
              vec1%arrh(3,:,:,:)=vec2%arrh(3,:,:,:)
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=AIMAG(vec2%arrv(1:2,:,:,:))
              vec1%arrv(3,:,:,:)=vec2%arrv(3,:,:,:)
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:))
              vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:)
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=AIMAG(vec2%arrtmp(1:2,:,:,:))
              vec1%arrtmp(3,:,:,:)=vec2%arrtmp(3,:,:,:)
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec2: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri=vec2%arri(:,:,:,:)
          CASE ('imag','IMAG')
            vec1%arri=AIMAG(vec2%arri(:,:,:,:))
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:)
            vec1%arri(3,:,:,:)=-AIMAG(vec2%arri(3,:,:,:))
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:)=AIMAG(vec2%arri(1:2,:,:,:))
            vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:)
          CASE DEFAULT
            CALL nim_stop
     $        ('Vector_assign_cvec2: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Vector_assign_cvec2: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_assign_cvec2
c-----------------------------------------------------------------------
c     subprogram 8. vector_ptassign_bc.
c     make a pointer assignment of real vector data to bicube data.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_bc(vec,bc)
      USE bicube

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(bicube_type), INTENT(IN), TARGET :: bc

      IF (ALLOCATED(bc%fs)) THEN
        vec%arr=>bc%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri,vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_bc
c-----------------------------------------------------------------------
c     subprogram 9. vector_ptassign_laq2.
c     make a pointer assignment of real vector data to 2D lagrange
c     quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_laq2(vec,laq)
      USE lagr_quad_mod

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(lagr_quad_2D_type), INTENT(IN), TARGET :: laq

      IF (ALLOCATED(laq%fs)) THEN
        vec%arr=>laq%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF

      IF (ALLOCATED(laq%fsh)) THEN
        vec%arrh=>laq%fsh
      ELSE
        NULLIFY(vec%arrh)
      ENDIF

      IF (ALLOCATED(laq%fsv)) THEN
        vec%arrv=>laq%fsv
      ELSE
        NULLIFY(vec%arrv)
      ENDIF

      IF (ALLOCATED(laq%fsi)) THEN
        vec%arri=>laq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_laq2
c-----------------------------------------------------------------------
c     subprogram 9.1. vector_ptassign_modq2.
c     make a pointer assignment of real vector data to 2D modal
c     quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_modq2(vec,modq)
      USE modal_type_mod

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(modal_quad_2D_type), INTENT(IN), TARGET :: modq

      IF (ALLOCATED(modq%fsi)) THEN
        vec%arri=>modq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arr,vec%arrh,vec%arrv,vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_modq2
c-----------------------------------------------------------------------
c     subprogram 10. vector_ptassign_tl2.
c     make a pointer assignment of real vector data to 2D tri_linear
c     data.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ptassign_tl2(vec,tl2)
      USE tri_linear

      TYPE(vector_type), INTENT(OUT) :: vec
      TYPE(tri_linear_2D_type), INTENT(IN), TARGET :: tl2

      IF (ALLOCATED(tl2%fs)) THEN
        vec%arr=>tl2%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri,vec%arrtmp)

      RETURN
      END SUBROUTINE vector_ptassign_tl2
c-----------------------------------------------------------------------
c     subprogram 11. vector_add_vec.
c     add one vector structure to another.
c-----------------------------------------------------------------------
      SUBROUTINE vector_add_vec(vec1,vec2,v1fac,v2fac)

      TYPE(vector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
c-----------------------------------------------------------------------
c     set coefficients to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Vector_add_vec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_add_vec
c-----------------------------------------------------------------------
c     subprogram 12. vector_mult_rsc.
c     multiply a vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE vector_mult_rsc(vec,rsc)

      TYPE(vector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Vector_mult_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_mult_rsc
c-----------------------------------------------------------------------
c     subprogram 13. vector_mult_int.
c     multiply a vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE vector_mult_int(vec,int)

      TYPE(vector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop
     $    ('Vector_mult_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_mult_int
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE rvector_type_mod


c-----------------------------------------------------------------------
c     module for 3D-in-space arrays of vectors quantities.
c-----------------------------------------------------------------------
      MODULE cvector_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 14. vector_ctype_alloc.
c     allocates space for a complex vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ctype_alloc(cvt,poly_degree,mx,my,nqty,nfour,
     $                              nbt,nqt,alloc_int)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,nfour,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(cvector_type), INTENT(OUT) :: cvt
      LOGICAL, INTENT(IN), OPTIONAL :: alloc_int
c-----------------------------------------------------------------------
c     allocate space according to the basis functions needed.
c     if poly_degree is non-positive, storage for discontinuous-field
c     coefficients is allocated.
c
c     interior-basis arrays are allocated by default unless prevented
c     by the optional alloc_int being set to false.
c-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,mx,my,nfour))
        NULLIFY(cvt%arr,cvt%arrh,cvt%arrv)
      CASE(1)  !  linear elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my,nfour))
        NULLIFY(cvt%arri,cvt%arrh,cvt%arrv)
      CASE(2:)  !  higher-order elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my,nfour))
        ALLOCATE(cvt%arrh(nqty,poly_degree-1,1:mx,0:my,nfour))
        ALLOCATE(cvt%arrv(nqty,poly_degree-1,0:mx,1:my,nfour))
        IF (.NOT.PRESENT(alloc_int)) THEN
          ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nfour))
        ELSE IF (alloc_int) THEN
          ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nfour))
        ELSE
          NULLIFY(cvt%arri)
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     if the optional input for temporary arrays is present, allocate
c     the arrtmp arrays; otherwise, nullify them.  note that nbt is
c     the number of temporary bases, and nqt is the number of quantities
c     at each basis.
c-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(cvt%arrtmp(nqt,nbt,1:mx,1:my,nfour))
      ELSE
        NULLIFY(cvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_ctype_alloc
c-----------------------------------------------------------------------
c     subprogram 15. vector_ctype_dealloc.
c     deallocates space for a complex vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_ctype_dealloc(cvt)

      TYPE(cvector_type), INTENT(INOUT) :: cvt

      IF (ASSOCIATED(cvt%arr)) THEN
        DEALLOCATE(cvt%arr)
        NULLIFY(cvt%arr)
      ENDIF
      IF (ASSOCIATED(cvt%arrh)) THEN
        DEALLOCATE(cvt%arrh)
        NULLIFY(cvt%arrh)
      ENDIF
      IF (ASSOCIATED(cvt%arrv)) THEN
        DEALLOCATE(cvt%arrv)
        NULLIFY(cvt%arrv)
      ENDIF
      IF (ASSOCIATED(cvt%arri)) THEN
        DEALLOCATE(cvt%arri)
        NULLIFY(cvt%arri)
      ENDIF
      IF (ASSOCIATED(cvt%arrtmp)) THEN
        DEALLOCATE(cvt%arrtmp)
        NULLIFY(cvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_ctype_dealloc
c-----------------------------------------------------------------------
c     subprogram 16. cvector_assign_rsc.
c     assign a real scalar value to a complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_rsc(vec,rscalar)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop
     $    ('Cvector_assign_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 17. cvector_assign_csc.
c     assign a complex scalar value to a complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_csc(vec,cscalar)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop
     $    ('Cvector_assign_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_csc
c-----------------------------------------------------------------------
c     subprogram 18. cvector_assign_int.
c     assign a integer value to a complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_int(vec,int)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop
     $    ('Cvector_assign_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_int
c-----------------------------------------------------------------------
c     subprogram 19. cvector_assign_cvec.
c     set one complex vector structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_cvec(vec1,vec2)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:,:)=vec2%arr(:,:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:,:)=vec2%arrh(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv(:,:,:,:,:)=vec2%arrv(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:,:)=vec2%arrtmp(:,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Cvector_assign_cvec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_cvec
c-----------------------------------------------------------------------
c     subprogram 19.1 cvector_assignq_cvec.
c     set a quantity range of one complex vector structure equal to
c     that of another.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assignq_cvec(vec1,vec2,nqty,nstart1,nstart2)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: nqty
      INTEGER(i4), INTENT(IN), OPTIONAL :: nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (PRESENT(nstart1)) THEN
        nst1=nstart1
      ELSE
        nst1=1
      ENDIF
      IF (PRESENT(nstart2)) THEN
        nst2=nstart2
      ELSE
        nst2=1
      ENDIF
      nend1=nst1+nqty-1
      nend2=nst2+nqty-1
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(nst1:nend1,:,:,:)=vec2%arr(nst2:nend2,:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(nst1:nend1,:,:,:,:)=vec2%arrh(nst2:nend2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv(nst1:nend1,:,:,:,:)=vec2%arrv(nst2:nend2,:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(nst1:nend1,:,:,:,:)=vec2%arri(nst2:nend2,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(nst1:nend1,:,:,:,:)=
     $    vec2%arrtmp(nst2:nend2,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(nst1:nend1,:,:,:,:)=vec2%arri(nst2:nend2,:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Cvector_assignq_cvec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assignq_cvec
c-----------------------------------------------------------------------
c     subprogram 20. cvector_assign_vec.
c     set one component of a complex vector structure equal to a real
c     vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_vec(vec1,vec2,r_i,fcomp,nqty,nstart1,
     $                              nstart2)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(nst1:nend1,:,:,fcomp)=vec2%arr(nst2:nend2,:,:)
     $        +(0,1)*AIMAG(vec1%arr(nst1:nend1,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:,fcomp)=
     $          vec2%arrh(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrh(nst1:nend1,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:,fcomp)=
     $          vec2%arrv(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrv(nst1:nend1,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:,fcomp)=
     $          vec2%arri(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arri(nst1:nend1,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:,fcomp)=
     $          vec2%arrtmp(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrtmp(nst1:nend1,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arr(nst1:nend1,:,:,fcomp)=
     $        (0,1)*vec2%arr(nst2:nend2,:,:)
     $        +REAL(vec1%arr(nst1:nend1,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:,fcomp)=
     $          (0,1)*vec2%arrh(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrh(nst1:nend1,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:,fcomp)=
     $          (0,1)*vec2%arrv(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrv(nst1:nend1,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:,fcomp)=
     $          (0,1)*vec2%arri(nst2:nend2,:,:,:)
     $          +REAL(vec1%arri(nst1:nend1,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:,fcomp)=
     $          (0,1)*vec2%arrtmp(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrtmp(nst1:nend1,:,:,:,fcomp),r8)
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(nst1:nend1,:,:,:,fcomp)=
     $        vec2%arri(nst2:nend2,:,:,:)
     $        +(0,1)*AIMAG(vec1%arri(nst1:nend1,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arri(nst1:nend1,:,:,:,fcomp)=
     $        (0,1)*vec2%arri(nst2:nend2,:,:,:)
     $        +REAL(vec1%arri(nst1:nend1,:,:,:,fcomp),r8)
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Cvector_assign_vec: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c     the r12mi3 and i12r3 flags are used in several of the
c     nimrod management routines.  r12mi3 means transfer the real
c     1 & 2 vector components and minus the third imaginary  comp.
c     i12r3 means transfer the imaginary 1 & 2 and the real third.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(:,:,:,fcomp)=vec2%arr
     $               +(0,1)*AIMAG(vec1%arr(:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(:,:,:,:,fcomp)=vec2%arrh
     $                    +(0,1)*AIMAG(vec1%arrh(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(:,:,:,:,fcomp)=vec2%arrv
     $                    +(0,1)*AIMAG(vec1%arrv(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(:,:,:,:,fcomp)=vec2%arri
     $                    +(0,1)*AIMAG(vec1%arri(:,:,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(:,:,:,:,fcomp)=vec2%arrtmp
     $                    +(0,1)*AIMAG(vec1%arrtmp(:,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arr(:,:,:,fcomp)=(0,1)*vec2%arr
     $                            +REAL(vec1%arr(:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(:,:,:,:,fcomp)=(0,1)*vec2%arrh
     $                                +REAL(vec1%arrh(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(:,:,:,:,fcomp)=(0,1)*vec2%arrv
     $                                +REAL(vec1%arrv(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(:,:,:,:,fcomp)=(0,1)*vec2%arri
     $                                +REAL(vec1%arri(:,:,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(:,:,:,:,fcomp)=(0,1)*vec2%arrtmp
     $                              +REAL(vec1%arrtmp(:,:,:,:,fcomp),r8)
          CASE ('r12mi3')
            vec1%arr(1:2,:,:,fcomp)=vec2%arr(1:2,:,:)
     $               +(0,1)*AIMAG(vec1%arr(1:2,:,:,fcomp))
            vec1%arr(3,:,:,fcomp)=-(0,1)*vec2%arr(3,:,:)
     $               +REAL(vec1%arr(3,:,:,fcomp),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:,fcomp)=vec2%arrh(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrh(1:2,:,:,:,fcomp))
              vec1%arrh(3,:,:,:,fcomp)=-(0,1)*vec2%arrh(3,:,:,:)
     $                 +REAL(vec1%arrh(3,:,:,:,fcomp),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:,fcomp)=vec2%arrv(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrv(1:2,:,:,:,fcomp))
              vec1%arrv(3,:,:,:,fcomp)=-(0,1)*vec2%arrv(3,:,:,:)
     $                 +REAL(vec1%arrv(3,:,:,:,fcomp),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:,fcomp)=vec2%arri(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arri(1:2,:,:,:,fcomp))
              vec1%arri(3,:,:,:,fcomp)=-(0,1)*vec2%arri(3,:,:,:)
     $                 +REAL(vec1%arri(3,:,:,:,fcomp),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:,fcomp)=vec2%arrtmp(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrtmp(1:2,:,:,:,fcomp))
              vec1%arrtmp(3,:,:,:,fcomp)=-(0,1)*vec2%arrtmp(3,:,:,:)
     $                 +REAL(vec1%arrtmp(3,:,:,:,fcomp),r8)
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:,fcomp)=(0,1)*vec2%arr(1:2,:,:)
     $               +REAL(vec1%arr(1:2,:,:,fcomp),r8)
            vec1%arr(3,:,:,fcomp)=vec2%arr(3,:,:)
     $               +(0,1)*AIMAG(vec1%arr(3,:,:,fcomp))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:,fcomp)=(0,1)*vec2%arrh(1:2,:,:,:)
     $                  +REAL(vec1%arrh(1:2,:,:,:,fcomp),r8)
              vec1%arrh(3,:,:,:,fcomp)=vec2%arrh(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrh(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:,fcomp)=(0,1)*vec2%arrv(1:2,:,:,:)
     $                  +REAL(vec1%arrv(1:2,:,:,:,fcomp),r8)
              vec1%arrv(3,:,:,:,fcomp)=vec2%arrv(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrv(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:,fcomp)=(0,1)*vec2%arri(1:2,:,:,:)
     $                  +REAL(vec1%arri(1:2,:,:,:,fcomp),r8)
              vec1%arri(3,:,:,:,fcomp)=vec2%arri(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arri(3,:,:,:,fcomp))
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:,fcomp)=(0,1)*vec2%arrtmp(1:2,:,:,:)
     $                  +REAL(vec1%arrtmp(1:2,:,:,:,fcomp),r8)
              vec1%arrtmp(3,:,:,:,fcomp)=vec2%arrtmp(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrtmp(3,:,:,:,fcomp))
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(:,:,:,:,fcomp)=vec2%arri
     $                  +(0,1)*AIMAG(vec1%arri(:,:,:,:,fcomp))
          CASE ('imag','IMAG')
            vec1%arri(:,:,:,:,fcomp)=(0,1)*vec2%arri
     $                               +REAL(vec1%arri(:,:,:,:,fcomp),r8)
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:,fcomp)=vec2%arri(1:2,:,:,:)
     $               +(0,1)*AIMAG(vec1%arri(1:2,:,:,:,fcomp))
            vec1%arri(3,:,:,:,fcomp)=-(0,1)*vec2%arri(3,:,:,:)
     $               +REAL(vec1%arri(3,:,:,:,fcomp),r8)
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:,fcomp)=(0,1)*vec2%arri(1:2,:,:,:)
     $                +REAL(vec1%arri(1:2,:,:,:,fcomp),r8)
            vec1%arri(3,:,:,:,fcomp)=vec2%arri(3,:,:,:)
     $                +(0,1)*AIMAG(vec1%arri(3,:,:,:,fcomp))
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Cvector_assign_vec: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_vec
c-----------------------------------------------------------------------
c     subprogram 21. cvector_assign_cvec2.
c     set one component of a complex vector structure equal to a complex
c     2D vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_assign_cvec2(vec1,vec2,fcomp,nqty,nstart1,
     $                                nstart2)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          vec1%arr(nst1:nend1,:,:,fcomp)=vec2%arr(nst2:nend2,:,:)
          IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $      vec1%arrh(nst1:nend1,:,:,:,fcomp)=
     $      vec2%arrh(nst2:nend2,:,:,:)
          IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $      vec1%arrv(nst1:nend1,:,:,:,fcomp)=
     $      vec2%arrv(nst2:nend2,:,:,:)
          IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $      vec1%arri(nst1:nend1,:,:,:,fcomp)=
     $      vec2%arri(nst2:nend2,:,:,:)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $      vec1%arrtmp(nst1:nend1,:,:,:,fcomp)=
     $      vec2%arrtmp(nst2:nend2,:,:,:)
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          vec1%arri(nst1:nend1,:,:,:,fcomp)=vec2%arri(nst2:nend2,:,:,:)
        ELSE
          CALL nim_stop
     $      ('Cvector_assign_cvec2: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          vec1%arr(:,:,:,fcomp)=vec2%arr
          IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $      vec1%arrh(:,:,:,:,fcomp)=vec2%arrh
          IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $      vec1%arrv(:,:,:,:,fcomp)=vec2%arrv
          IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $      vec1%arri(:,:,:,:,fcomp)=vec2%arri
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $      vec1%arrtmp(:,:,:,:,fcomp)=vec2%arrtmp
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          vec1%arri(:,:,:,:,fcomp)=vec2%arri
        ELSE
          CALL nim_stop
     $      ('Cvector_assign_cvec2: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_assign_cvec2
c-----------------------------------------------------------------------
c     subprogram 22. cvector_ptassign_laq.
c     make a pointer assignment of complex vector data to lagrange
c     quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_laq(vec,laq)
      USE lagr_quad_mod

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(lagr_quad_type), INTENT(IN), TARGET :: laq

      IF (ALLOCATED(laq%fs)) THEN
        vec%arr=>laq%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF

      IF (ALLOCATED(laq%fsh)) THEN
        vec%arrh=>laq%fsh
      ELSE
        NULLIFY(vec%arrh)
      ENDIF

      IF (ALLOCATED(laq%fsv)) THEN
        vec%arrv=>laq%fsv
      ELSE
        NULLIFY(vec%arrv)
      ENDIF

      IF (ALLOCATED(laq%fsi)) THEN
        vec%arri=>laq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_laq
c-----------------------------------------------------------------------
c     subprogram 22.1. cvector_ptassign_modq.
c     make a pointer assignment of complex vector data to modal
c     quadrilateral data.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_modq(vec,modq)
      USE modal_type_mod

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(modal_quad_type), INTENT(IN), TARGET :: modq

      IF (ALLOCATED(modq%fsi)) THEN
        vec%arri=>modq%fsi
      ELSE
        NULLIFY(vec%arri)
      ENDIF

      NULLIFY(vec%arr,vec%arrh,vec%arrv,vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_modq
c-----------------------------------------------------------------------
c     subprogram 23. cvector_ptassign_tl.
c     make a pointer assignment of complex vector data to 3D tri_linear
c     data.
c-PRE this will eventually need the additional arrays.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_ptassign_tl(vec,tl)
      USE tri_linear

      TYPE(cvector_type), INTENT(OUT) :: vec
      TYPE(tri_linear_type), INTENT(IN), TARGET :: tl

      IF (ALLOCATED(tl%fs)) THEN
        vec%arr=>tl%fs
      ELSE
        NULLIFY(vec%arr)
      ENDIF
      NULLIFY(vec%arrh,vec%arrv,vec%arri,vec%arrtmp)

      RETURN
      END SUBROUTINE cvector_ptassign_tl
c-----------------------------------------------------------------------
c     subprogram 24. cvector_add_cvec.
c     add one complex vector structure to another.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_add_cvec(vec1,vec2,v1fac,v2fac)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
c-----------------------------------------------------------------------
c     set coefficients to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_add_cvec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_add_cvec
c-----------------------------------------------------------------------
c     subprogram 25. cvector_addc_cvec.
c     add one complex vector structure to another a complex scalar
c     with complex coefficients.  here, only v2fac is optional, so that
c     there is no confusion with cvector_add_cvec in the module
c     procedure definition.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_addc_cvec(vec1,vec2,v1f,v2fac)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      COMPLEX(r8), INTENT(IN) :: v1f
      COMPLEX(r8), INTENT(IN), OPTIONAL :: v2fac

      COMPLEX(r8) :: v2f
c-----------------------------------------------------------------------
c     set coefficient to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1._r8
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_addc_cvec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_addc_cvec
c-----------------------------------------------------------------------
c     subprogram 26. cvector_mult_rsc.
c     multiply a complex vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_rsc(vec,rsc)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_mult_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_rsc
c-----------------------------------------------------------------------
c     subprogram 27. cvector_mult_csc.
c     multiply a complex vector structure by a complex scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_csc(vec,csc)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: csc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=csc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=csc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=csc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=csc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=csc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=csc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_mult_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_csc
c-----------------------------------------------------------------------
c     subprogram 28. cvector_mult_int.
c     multiply a complex vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_mult_int(vec,int)

      TYPE(cvector_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_mult_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_mult_int
c-----------------------------------------------------------------------
c     subprogram 28.1. cvector_real_comp.
c     set one component of a complex vector structure its real part.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_real_comp(vec1,fcomp)

      TYPE(cvector_type), INTENT(INOUT) :: vec1
      INTEGER(i4), INTENT(IN) :: fcomp

c-----------------------------------------------------------------------
c     use r8 in the intrinsic calls to maintain accuracy.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr))
     $  vec1%arr(:,:,:,fcomp)=REAL(vec1%arr(:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrh))
     $  vec1%arrh(:,:,:,:,fcomp)=REAL(vec1%arrh(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrv))
     $  vec1%arrv(:,:,:,:,fcomp)=REAL(vec1%arrv(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arri))
     $  vec1%arri(:,:,:,:,fcomp)=REAL(vec1%arri(:,:,:,:,fcomp),r8)
      IF (ASSOCIATED(vec1%arrtmp))
     $  vec1%arrtmp(:,:,:,:,fcomp)=REAL(vec1%arrtmp(:,:,:,:,fcomp),r8)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_real_comp
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE cvector_type_mod


c-----------------------------------------------------------------------
c     module for 2D-in-space arrays of complex vectors quantities.
c-----------------------------------------------------------------------
      MODULE cvector_2D_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 29. vector_2D_ctype_alloc.
c     allocates space for a 2D complex vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_2D_ctype_alloc(cvt,poly_degree,mx,my,nqty,
     $                                 nbt,nqt,alloc_int)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(cvector_2D_type), INTENT(OUT) :: cvt
      LOGICAL, INTENT(IN), OPTIONAL :: alloc_int
c-----------------------------------------------------------------------
c     allocate space according to the basis functions needed.
c     if poly_degree is non-positive, storage for discontinuous-field
c     coefficients is allocated.
c
c     interior-basis arrays are allocated by default unless prevented
c     by the optional alloc_int being set to false.
c-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,mx,my))
        NULLIFY(cvt%arr,cvt%arrh,cvt%arrv)
      CASE(1)  !  linear elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my))
        NULLIFY(cvt%arri,cvt%arrh,cvt%arrv)
      CASE(2:)  !  higher-order elements
        ALLOCATE(cvt%arr(nqty,0:mx,0:my))
        ALLOCATE(cvt%arrh(nqty,poly_degree-1,1:mx,0:my))
        ALLOCATE(cvt%arrv(nqty,poly_degree-1,0:mx,1:my))
        IF (.NOT.PRESENT(alloc_int)) THEN
          ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
        ELSE IF (alloc_int) THEN
          ALLOCATE(cvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my))
        ELSE
          NULLIFY(cvt%arri)
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     if the optional input for temporary arrays is present, allocate
c     the arrtmp arrays; otherwise, nullify them.  note that nbt is
c     the number of temporary bases, and nqt is the number of quantities
c     at each basis.
c-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(cvt%arrtmp(nqt,nbt,1:mx,1:my))
      ELSE
        NULLIFY(cvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_2D_ctype_alloc
c-----------------------------------------------------------------------
c     subprogram 30. vector_2D_ctype_dealloc.
c     deallocates space for a 2D complex vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_2D_ctype_dealloc(cvt)

      TYPE(cvector_2D_type), INTENT(INOUT) :: cvt

      IF (ASSOCIATED(cvt%arr)) THEN
        DEALLOCATE(cvt%arr)
        NULLIFY(cvt%arr)
      ENDIF
      IF (ASSOCIATED(cvt%arrh)) THEN
        DEALLOCATE(cvt%arrh)
        NULLIFY(cvt%arrh)
      ENDIF
      IF (ASSOCIATED(cvt%arrv)) THEN
        DEALLOCATE(cvt%arrv)
        NULLIFY(cvt%arrv)
      ENDIF
      IF (ASSOCIATED(cvt%arri)) THEN
        DEALLOCATE(cvt%arri)
        NULLIFY(cvt%arri)
      ENDIF
      IF (ASSOCIATED(cvt%arrtmp)) THEN
        DEALLOCATE(cvt%arrtmp)
        NULLIFY(cvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_2D_ctype_dealloc
c-----------------------------------------------------------------------
c     subprogram 31. cvector_2D_assign_rsc.
c     assign a real scalar value to a 2D complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_rsc(vec,rscalar)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_assign_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 32. cvector_2D_assign_csc.
c     assign a 2D complex scalar value to a 2D complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_csc(vec,cscalar)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_assign_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 33. cvector_2D_assign_int.
c     assign a integer value to a 2D complex vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_int(vec,int)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_assign_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_int
c-----------------------------------------------------------------------
c     subprogram 34. cvector_2D_assign_cvec2.
c     set one 2D complex vector structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_cvec2(vec1,vec2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv(:,:,:,:)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_assign_cvec2: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_cvec2
c-----------------------------------------------------------------------
c     subprogram 34.1. cvector_2D_assign_vec.
c     set a 2D complex vector structure equal to a real vector with the
c     the real/imaginary placement options of cvector_assign_vec.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_vec(vec1,vec2,r_i,nqty,nstart1,
     $                                 nstart2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      CHARACTER(*), INTENT(IN) :: r_i
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(nst1:nend1,:,:)=vec2%arr(nst2:nend2,:,:)
     $        +(0,1)*AIMAG(vec1%arr(nst1:nend1,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $          vec2%arrh(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrh(nst1:nend1,:,:,:))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $          vec2%arrv(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrv(nst1:nend1,:,:,:))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $          vec2%arri(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arri(nst1:nend1,:,:,:))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $          vec2%arrtmp(nst2:nend2,:,:,:)
     $          +(0,1)*AIMAG(vec1%arrtmp(nst1:nend1,:,:,:))
          CASE ('imag','IMAG')
            vec1%arr(nst1:nend1,:,:)=
     $        (0,1)*vec2%arr(nst2:nend2,:,:)
     $        +REAL(vec1%arr(nst1:nend1,:,:),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(nst1:nend1,:,:,:)=
     $          (0,1)*vec2%arrh(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrh(nst1:nend1,:,:,:),r8)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(nst1:nend1,:,:,:)=
     $          (0,1)*vec2%arrv(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrv(nst1:nend1,:,:,:),r8)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(nst1:nend1,:,:,:)=
     $          (0,1)*vec2%arri(nst2:nend2,:,:,:)
     $          +REAL(vec1%arri(nst1:nend1,:,:,:),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(nst1:nend1,:,:,:)=
     $          (0,1)*vec2%arrtmp(nst2:nend2,:,:,:)
     $          +REAL(vec1%arrtmp(nst1:nend1,:,:,:),r8)
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_2D_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(nst1:nend1,:,:,:)=
     $        vec2%arri(nst2:nend2,:,:,:)
     $        +(0,1)*AIMAG(vec1%arri(nst1:nend1,:,:,:))
          CASE ('imag','IMAG')
            vec1%arri(nst1:nend1,:,:,:)=
     $        (0,1)*vec2%arri(nst2:nend2,:,:,:)
     $        +REAL(vec1%arri(nst1:nend1,:,:,:),r8)
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_2D_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Cvector_2D_assign_vec: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c     the r12mi3 and i12r3 flags are used in several of the
c     nimrod management routines.  r12mi3 means transfer the real
c     1 & 2 vector components and minus the third imaginary  comp.
c     i12r3 means transfer the imaginary 1 & 2 and the real third.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arr(:,:,:)=vec2%arr
     $               +(0,1)*AIMAG(vec1%arr(:,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(:,:,:,:)=vec2%arrh
     $                    +(0,1)*AIMAG(vec1%arrh(:,:,:,:))
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(:,:,:,:)=vec2%arrv
     $                    +(0,1)*AIMAG(vec1%arrv(:,:,:,:))
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(:,:,:,:)=vec2%arri
     $                    +(0,1)*AIMAG(vec1%arri(:,:,:,:))
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(:,:,:,:)=vec2%arrtmp
     $                    +(0,1)*AIMAG(vec1%arrtmp(:,:,:,:))
          CASE ('imag','IMAG')
            vec1%arr(:,:,:)=(0,1)*vec2%arr
     $                            +REAL(vec1%arr(:,:,:),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $        vec1%arrh(:,:,:,:)=(0,1)*vec2%arrh
     $                                +REAL(vec1%arrh(:,:,:,:),r8)
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $        vec1%arrv(:,:,:,:)=(0,1)*vec2%arrv
     $                                +REAL(vec1%arrv(:,:,:,:),r8)
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $        vec1%arri(:,:,:,:)=(0,1)*vec2%arri
     $                                +REAL(vec1%arri(:,:,:,:),r8)
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        vec1%arrtmp(:,:,:,:)=(0,1)*vec2%arrtmp
     $                              +REAL(vec1%arrtmp(:,:,:,:),r8)
          CASE ('r12mi3')
            vec1%arr(1:2,:,:)=vec2%arr(1:2,:,:)
     $               +(0,1)*AIMAG(vec1%arr(1:2,:,:))
            vec1%arr(3,:,:)=-(0,1)*vec2%arr(3,:,:)
     $               +REAL(vec1%arr(3,:,:),r8)
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=vec2%arrh(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrh(1:2,:,:,:))
              vec1%arrh(3,:,:,:)=-(0,1)*vec2%arrh(3,:,:,:)
     $                 +REAL(vec1%arrh(3,:,:,:),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=vec2%arrv(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrv(1:2,:,:,:))
              vec1%arrv(3,:,:,:)=-(0,1)*vec2%arrv(3,:,:,:)
     $                 +REAL(vec1%arrv(3,:,:,:),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arri(1:2,:,:,:))
              vec1%arri(3,:,:,:)=-(0,1)*vec2%arri(3,:,:,:)
     $                 +REAL(vec1%arri(3,:,:,:),r8)
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=vec2%arrtmp(1:2,:,:,:)
     $                 +(0,1)*AIMAG(vec1%arrtmp(1:2,:,:,:))
              vec1%arrtmp(3,:,:,:)=-(0,1)*vec2%arrtmp(3,:,:,:)
     $                 +REAL(vec1%arrtmp(3,:,:,:),r8)
            ENDIF
          CASE ('i12r3')
            vec1%arr(1:2,:,:)=(0,1)*vec2%arr(1:2,:,:)
     $               +REAL(vec1%arr(1:2,:,:),r8)
            vec1%arr(3,:,:)=vec2%arr(3,:,:)
     $               +(0,1)*AIMAG(vec1%arr(3,:,:))
            IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
              vec1%arrh(1:2,:,:,:)=(0,1)*vec2%arrh(1:2,:,:,:)
     $                  +REAL(vec1%arrh(1:2,:,:,:),r8)
              vec1%arrh(3,:,:,:)=vec2%arrh(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrh(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
              vec1%arrv(1:2,:,:,:)=(0,1)*vec2%arrv(1:2,:,:,:)
     $                  +REAL(vec1%arrv(1:2,:,:,:),r8)
              vec1%arrv(3,:,:,:)=vec2%arrv(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrv(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
              vec1%arri(1:2,:,:,:)=(0,1)*vec2%arri(1:2,:,:,:)
     $                  +REAL(vec1%arri(1:2,:,:,:),r8)
              vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arri(3,:,:,:))
            ENDIF
            IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $        THEN
              vec1%arrtmp(1:2,:,:,:)=(0,1)*vec2%arrtmp(1:2,:,:,:)
     $                  +REAL(vec1%arrtmp(1:2,:,:,:),r8)
              vec1%arrtmp(3,:,:,:)=vec2%arrtmp(3,:,:,:)
     $                  +(0,1)*AIMAG(vec1%arrtmp(3,:,:,:))
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_2D_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
c-----------------------------------------------------------------------
c       cell-centered data only.
c-----------------------------------------------------------------------
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          SELECT CASE(r_i)
          CASE ('real','REAL')
            vec1%arri(:,:,:,:)=vec2%arri
     $                  +(0,1)*AIMAG(vec1%arri(:,:,:,:))
          CASE ('imag','IMAG')
            vec1%arri(:,:,:,:)=(0,1)*vec2%arri
     $                               +REAL(vec1%arri(:,:,:,:),r8)
          CASE ('r12mi3')
            vec1%arri(1:2,:,:,:)=vec2%arri(1:2,:,:,:)
     $               +(0,1)*AIMAG(vec1%arri(1:2,:,:,:))
            vec1%arri(3,:,:,:)=-(0,1)*vec2%arri(3,:,:,:)
     $               +REAL(vec1%arri(3,:,:,:),r8)
          CASE ('i12r3')
            vec1%arri(1:2,:,:,:)=(0,1)*vec2%arri(1:2,:,:,:)
     $                +REAL(vec1%arri(1:2,:,:,:),r8)
            vec1%arri(3,:,:,:)=vec2%arri(3,:,:,:)
     $                +(0,1)*AIMAG(vec1%arri(3,:,:,:))
          CASE DEFAULT
            CALL nim_stop
     $        ('Cvector_2D_assign_vec: '//r_i//' flag not recognized.')
          END SELECT
        ELSE
          CALL nim_stop
     $      ('Cvector_2D_assign_vec: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_vec
c-----------------------------------------------------------------------
c     subprogram 35. cvector_2D_add_cvec2.
c     add one 2D complex vector structure to another.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_add_cvec2(vec1,vec2,v1fac,v2fac)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
c-----------------------------------------------------------------------
c     set coefficients to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_add_cvec2: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_add_cvec2
c-----------------------------------------------------------------------
c     subprogram 36. cvector_2D_addc_cvec2.
c     add one 2D complex vector structure to another with a complex
c     coefficients.  here, only v2fac is optional, so that
c     there is no confusion with cvector_add_cvec in the module
c     procedure definition.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_addc_cvec2(vec1,vec2,v1f,v2fac)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_2D_type), INTENT(IN) :: vec2
      COMPLEX(r8), INTENT(IN) :: v1f
      COMPLEX(r8), INTENT(IN), OPTIONAL :: v2fac

      COMPLEX(r8) :: v2f
c-----------------------------------------------------------------------
c     set coefficient to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1._r8
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_addc_cvec2: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_addc_cvec2
c-----------------------------------------------------------------------
c     subprogram 37. cvector_2D_assign_cvec.
c     set one complex 2D vector structure equal to part of a complex
c     vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_assign_cvec(vec1,vec2,fcomp,nqty,nstart1,
     $                                  nstart2)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec1
      TYPE(cvector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: fcomp
      INTEGER(i4), INTENT(IN), OPTIONAL :: nqty,nstart1,nstart2

      INTEGER(i4) :: nend1,nend2,nst1,nst2
c-----------------------------------------------------------------------
c     if the number of vector components is specified with nqty, limit
c     transfer.
c-----------------------------------------------------------------------
      IF (PRESENT(nqty)) THEN
        IF (PRESENT(nstart1)) THEN
          nst1=nstart1
        ELSE
          nst1=1
        ENDIF
        IF (PRESENT(nstart2)) THEN
          nst2=nstart2
        ELSE
          nst2=1
        ENDIF
        nend1=nst1+nqty-1
        nend2=nst2+nqty-1
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          vec1%arr(nst1:nend1,:,:)=vec2%arr(nst2:nend2,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $      vec1%arrh(nst1:nend1,:,:,:)=
     $      vec2%arrh(nst2:nend2,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $      vec1%arrv(nst1:nend1,:,:,:)=
     $      vec2%arrv(nst2:nend2,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $      vec1%arri(nst1:nend1,:,:,:)=
     $      vec2%arri(nst2:nend2,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $      vec1%arrtmp(nst1:nend1,:,:,:)=
     $      vec2%arrtmp(nst2:nend2,:,:,:,fcomp)
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          vec1%arri(nst1:nend1,:,:,:)=vec2%arri(nst2:nend2,:,:,:,fcomp)
        ELSE
          CALL nim_stop
     $      ('Cvector_2D_assign_cvec: vector arrays not associated.')
        ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      ELSE
        IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
          vec1%arr=vec2%arr(:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $      vec1%arrh=vec2%arrh(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $      vec1%arrv=vec2%arrv(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $      vec1%arri=vec2%arri(:,:,:,:,fcomp)
          IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $      vec1%arrtmp=vec2%arrtmp(:,:,:,:,fcomp)
        ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          vec1%arri=vec2%arri(:,:,:,:,fcomp)
        ELSE
          CALL nim_stop
     $      ('Cvector_2D_assign_cvec: vector arrays not associated.')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_assign_cvec
c-----------------------------------------------------------------------
c     subprogram 38. cvector_2D_mult_rsc.
c     multiply a complex 2D vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_rsc(vec,rsc)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_mult_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_rsc
c-----------------------------------------------------------------------
c     subprogram 39. cvector_2D_mult_csc.
c     multiply a complex 2D vector structure by a complex scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_csc(vec,csc)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: csc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=csc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=csc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=csc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=csc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=csc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=csc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_mult_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_csc
c-----------------------------------------------------------------------
c     subprogram 40. cvector_2D_mult_int.
c     multiply a complex 2D vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_mult_int(vec,int)

      TYPE(cvector_2D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop
     $    ('Cvector_2D_mult_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_mult_int
c-----------------------------------------------------------------------
c     subprogram 40.1. cvector_2D_pack_cvec.
c     packs the interior and discontinuous coefficients of a cvector
c     into the interior coefficients of a complex 2D vector.
c
c     the arri array of the 2D vector should be allocated to hold all
c     of the quantity and basis indices as a flat array in its first
c     array index.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_pack_cvec(cvt,cv2p,nqi,nbi,nqd,nbd,fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_type), INTENT(IN) :: cvt
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2p

      INTEGER(i4) :: ix,iy,iv,ib,iq
c-----------------------------------------------------------------------
c     at each element, loop over the interior and discontinuous bases
c     separately.
c-----------------------------------------------------------------------
      DO iy=1,SIZE(cvt%arri,4)
        DO ix=1,SIZE(cvt%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2p%arri(iv,1,ix,iy)=cvt%arri(iq,ib,ix,iy,fcomp)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cv2p%arri(iv,1,ix,iy)=cvt%arrtmp(iq,ib,ix,iy,fcomp)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_pack_cvec
c-----------------------------------------------------------------------
c     subprogram 40.2. cvector_2D_pack_cvec2.
c     packs the interior and discontinuous coefficients of a 2D complex
c     vector into the interior coefficients of another complex 2D
c     vector.
c
c     the arri array of the 2D vector should be allocated to hold all
c     of the quantity and basis indices as a flat array in its first
c     array index.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_pack_cvec2(cv2in,cv2p,nqi,nbi,nqd,nbd)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd
      TYPE(cvector_2D_type), INTENT(IN) :: cv2in
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2p

      INTEGER(i4) :: ix,iy,iv,ib,iq
c-----------------------------------------------------------------------
c     at each element, loop over the interior and discontinuous bases
c     separately.
c-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2in%arri,4)
        DO ix=1,SIZE(cv2in%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2p%arri(iv,1,ix,iy)=cv2in%arri(iq,ib,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cv2p%arri(iv,1,ix,iy)=cv2in%arrtmp(iq,ib,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_pack_cvec2
c-----------------------------------------------------------------------
c     subprogram 40.3. cvector_2D_unpack_cvec.
c     unpacks the interior coefficients of a complex 2D vector into
c     another complex 2D vector and a cvector, respectively.
c
c     the arri array of the packed 2D vector is allocated to hold all
c     of the quantity and basis indices as a flat array in its first
c     array index.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_cvec(cv2p,cv2i,cvt,nqi,nbi,nqd,nbd,
     $                                  fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i
      TYPE(cvector_type), INTENT(OUT) :: cvt

      INTEGER(i4) :: ix,iy,iv,ib,iq
c-----------------------------------------------------------------------
c     at each element, loop over the interior and discontinuous bases
c     separately.
c-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cvt%arrtmp(iq,ib,ix,iy,fcomp)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_cvec
c-----------------------------------------------------------------------
c     subprogram 40.4. cvector_2D_unpack_add_cvec.
c     unpacks the interior coefficients of a complex 2D vector into
c     another complex 2D vector and a cvector, respectively.
c     here the result is added to the interior of the cvector.
c
c     the arri array of the packed 2D vector is allocated to hold all
c     of the quantity and basis indices as a flat array in its first
c     array index.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_add_cvec(cv2p,cv2i,cvt,nqi,nbi,nqd,
     $                                      nbd,fcomp)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd,fcomp
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i
      TYPE(cvector_type), INTENT(INOUT) :: cvt

      INTEGER(i4) :: ix,iy,iv,ib,iq
c-----------------------------------------------------------------------
c     at each element, loop over the interior and discontinuous bases
c     separately.
c-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          DO ib=1,nbd
            DO iq=1,nqd
              cvt%arri(iq,ib,ix,iy,fcomp)=
     $          cvt%arri(iq,ib,ix,iy,fcomp)+cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_add_cvec
c-----------------------------------------------------------------------
c     subprogram 40.5. cvector_2D_unpack_cvec2.
c     unpacks the interior coefficients of a complex 2D vector into
c     another complex 2D vector, skipping the indices from a
c     discontinuous field if unpd is false.
c
c     the arri array of the packed 2D vector is allocated to hold all
c     of the quantity and basis indices as a flat array in its first
c     array index.
c-----------------------------------------------------------------------
      SUBROUTINE cvector_2D_unpack_cvec2(cv2p,cv2i,nqi,nbi,nqd,nbd,unpd)

      INTEGER(i4), INTENT(IN) :: nqi,nbi,nqd,nbd
      LOGICAL, INTENT(IN) :: unpd
      TYPE(cvector_2D_type), INTENT(IN) :: cv2p
      TYPE(cvector_2D_type), INTENT(OUT) :: cv2i

      INTEGER(i4) :: ix,iy,iv,ib,iq
c-----------------------------------------------------------------------
c     at each element, loop over the interior bases.
c-----------------------------------------------------------------------
      DO iy=1,SIZE(cv2p%arri,4)
        DO ix=1,SIZE(cv2p%arri,3)
          iv=1
          DO ib=1,nbi
            DO iq=1,nqi
              cv2i%arri(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
              iv=iv+1
            ENDDO
          ENDDO
          IF (unpd) THEN
            DO ib=1,nbd
              DO iq=1,nqd
                cv2i%arrtmp(iq,ib,ix,iy)=cv2p%arri(iv,1,ix,iy)
                iv=iv+1
              ENDDO
            ENDDO
          ELSE
            iv=iv+nqd*nbd
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cvector_2D_unpack_cvec2
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE cvector_2D_type_mod


c-----------------------------------------------------------------------
c     module for 3D-in-space arrays of real vectors quantities.
c-----------------------------------------------------------------------
      MODULE rvector_3D_type_mod
      USE local
      USE vector_defn_type_mod
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 41. vector_3D_rtype_alloc.
c     allocates space for a real 3D vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_rtype_alloc(rvt,poly_degree,mx,my,nqty,nph,
     $                                 nbt,nqt,alloc_int)

      INTEGER(i4), INTENT(IN) :: mx,my,nqty,poly_degree,nph
      INTEGER(i4), INTENT(IN), OPTIONAL :: nbt,nqt
      TYPE(vector_3D_type), INTENT(OUT) :: rvt
      LOGICAL, INTENT(IN), OPTIONAL :: alloc_int
c-----------------------------------------------------------------------
c     allocate space according to the basis functions needed.
c     if poly_degree is non-positive, storage for discontinuous-field
c     coefficients is allocated.
c
c     interior-basis arrays are allocated by default unless prevented
c     by the optional alloc_int being set to false.
c
c-PRE triangles will need something here, too.
c-----------------------------------------------------------------------
      SELECT CASE(poly_degree)
      CASE(:0)  !  piecewise continuous fields
        ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,mx,my,nph))
        NULLIFY(rvt%arr,rvt%arrh,rvt%arrv)
      CASE(1)  !  linear elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my,nph))
        NULLIFY(rvt%arri,rvt%arrh,rvt%arrv)
      CASE(2:)  !  higher-order elements
        ALLOCATE(rvt%arr(nqty,0:mx,0:my,nph))
        ALLOCATE(rvt%arrh(nqty,poly_degree-1,1:mx,0:my,nph))
        ALLOCATE(rvt%arrv(nqty,poly_degree-1,0:mx,1:my,nph))
        IF (.NOT.PRESENT(alloc_int)) THEN
          ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nph))
        ELSE IF (alloc_int) THEN
          ALLOCATE(rvt%arri(nqty,(poly_degree-1)**2,1:mx,1:my,nph))
        ELSE
          NULLIFY(rvt%arri)
        ENDIF
      END SELECT
c-----------------------------------------------------------------------
c     if the optional input for temporary arrays is present, allocate
c     the arrtmp arrays; otherwise, nullify them.  note that nbt is
c     the number of temporary bases, and nqt is the number of quantities
c     at each basis.
c-----------------------------------------------------------------------
      IF (PRESENT(nbt)) THEN
        ALLOCATE(rvt%arrtmp(nqt,nbt,1:mx,1:my,nph))
      ELSE
        NULLIFY(rvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_rtype_alloc
c-----------------------------------------------------------------------
c     subprogram 42. vector_3D_rtype_dealloc.
c     deallocates space for a real 3D vector_type structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_rtype_dealloc(rvt)

      TYPE(vector_3D_type), INTENT(INOUT) :: rvt

      IF (ASSOCIATED(rvt%arr)) THEN
        DEALLOCATE(rvt%arr)
        NULLIFY(rvt%arr)
      ENDIF
      IF (ASSOCIATED(rvt%arrh)) THEN
        DEALLOCATE(rvt%arrh)
        NULLIFY(rvt%arrh)
      ENDIF
      IF (ASSOCIATED(rvt%arrv)) THEN
        DEALLOCATE(rvt%arrv)
        NULLIFY(rvt%arrv)
      ENDIF
      IF (ASSOCIATED(rvt%arri)) THEN
        DEALLOCATE(rvt%arri)
        NULLIFY(rvt%arri)
      ENDIF
      IF (ASSOCIATED(rvt%arrtmp)) THEN
        DEALLOCATE(rvt%arrtmp)
        NULLIFY(rvt%arrtmp)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_rtype_dealloc
c-----------------------------------------------------------------------
c     subprogram 43. vector_3D_assign_rsc.
c     assign a real scalar value to a 3D vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_rsc(vec,rscalar)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=rscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rscalar
      ELSE
        CALL nim_stop
     $    ('Vector_3D_assign_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_rsc
c-----------------------------------------------------------------------
c     subprogram 44. vector_3D_assign_csc.
c     assign the real part of a complex scalar value to a 3D vector
c     structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_csc(vec,cscalar)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      COMPLEX(r8), INTENT(IN) :: cscalar

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=cscalar
        IF (ASSOCIATED(vec%arrh)) vec%arrh=cscalar
        IF (ASSOCIATED(vec%arrv)) vec%arrv=cscalar
        IF (ASSOCIATED(vec%arri)) vec%arri=cscalar
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=cscalar
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=cscalar
      ELSE
        CALL nim_stop
     $    ('Vector_3D_assign_csc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_csc
c-----------------------------------------------------------------------
c     subprogram 45. vector_3D_assign_int.
c     assign an integer value to a 3D vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_int(vec,int)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     if the grid vertex-centered array is allocated, treat as a 
c     standard element.  If, not the structure represents piecewise
c     constant.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int
        IF (ASSOCIATED(vec%arri)) vec%arri=int
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int
      ELSE
        CALL nim_stop
     $    ('Vector_3D_assign_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_int
c-----------------------------------------------------------------------
c     subprogram 46. vector_3D_assign_vec.
c     set a 3D vector structure equal to a vector structure at the
c     specified index of the periodic coordinate.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_vec(vec1,vec2,ip)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      INTEGER(i4), INTENT(IN) :: ip
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:,ip)=vec2%arr(:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:,ip)=vec2%arrh(:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) 
     $    vec1%arrv(:,:,:,:,ip)=vec2%arrv(:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:,ip)=vec2%arri(:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:,ip)=vec2%arrtmp(:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:,ip)=vec2%arri(:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Vector_3D_assign_vec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_vec
c-----------------------------------------------------------------------
c     subprogram 47. vector_3D_assign_vec3.
c     set a 3D vector structure equal to another.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_assign_vec3(vec1,vec2)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr(:,:,:,:)=vec2%arr(:,:,:,:)
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh(:,:,:,:,:)=vec2%arrh(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv(:,:,:,:,:)=vec2%arrv(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp(:,:,:,:,:)=vec2%arrtmp(:,:,:,:,:)
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri(:,:,:,:,:)=vec2%arri(:,:,:,:,:)
      ELSE
        CALL nim_stop
     $    ('Vector_3D_assign_vec3: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_assign_vec3
c-----------------------------------------------------------------------
c     subprogram 48. vector_3D_add_vec.
c     add a 2D vector structure to every index of the periodic
c     coordinate of a 3D vector structure.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_add_vec(vec1,vec2,v1fac,v2fac)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      INTEGER(i4) :: ip,np
      REAL(r8) :: v1f,v2f
c-----------------------------------------------------------------------
c     set coefficients to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        np=SIZE(vec1%arr,4)
        DO ip=1,np
          vec1%arr(:,:,:,ip)=v1f*vec1%arr(:,:,:,ip)+v2f*vec2%arr
        ENDDO
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh)) THEN
          DO ip=1,np
            vec1%arrh(:,:,:,:,ip)=v1f*vec1%arrh(:,:,:,:,ip)+
     $                            v2f*vec2%arrh
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv)) THEN
          DO ip=1,np
            vec1%arrv(:,:,:,:,ip)=v1f*vec1%arrv(:,:,:,:,ip)+
     $                            v2f*vec2%arrv
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
          DO ip=1,np
            vec1%arri(:,:,:,:,ip)=v1f*vec1%arri(:,:,:,:,ip)+
     $                            v2f*vec2%arri
          ENDDO
        ENDIF
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp)) THEN
          DO ip=1,np
            vec1%arrtmp(:,:,:,:,ip)=v1f*vec1%arrtmp(:,:,:,:,ip)+
     $                              v2f*vec2%arrtmp
          ENDDO
        ENDIF
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        np=SIZE(vec1%arri,5)
        DO ip=1,np
          vec1%arri(:,:,:,:,ip)=v1f*vec1%arri(:,:,:,:,ip)+v2f*vec2%arri
        ENDDO
      ELSE
        CALL nim_stop
     $    ('Vector_3D_add_vec: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_add_vec
c-----------------------------------------------------------------------
c     subprogram 49. vector_3D_add_vec3.
c     add one 3D vector structure to another.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_add_vec3(vec1,vec2,v1fac,v2fac)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec1
      TYPE(vector_3D_type), INTENT(IN) :: vec2
      REAL(r8), INTENT(IN), OPTIONAL :: v1fac,v2fac

      REAL(r8) :: v1f,v2f
c-----------------------------------------------------------------------
c     set coefficients to input if used.
c-----------------------------------------------------------------------
      IF (PRESENT(v1fac)) THEN
        v1f=v1fac
      ELSE
        v1f=1
      ENDIF
      IF (PRESENT(v2fac)) THEN
        v2f=v2fac
      ELSE
        v2f=1
      ENDIF
c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec1%arr).AND.ASSOCIATED(vec2%arr)) THEN
        vec1%arr=v1f*vec1%arr+v2f*vec2%arr
        IF (ASSOCIATED(vec1%arrh).AND.ASSOCIATED(vec2%arrh))
     $    vec1%arrh=v1f*vec1%arrh+v2f*vec2%arrh
        IF (ASSOCIATED(vec1%arrv).AND.ASSOCIATED(vec2%arrv))
     $    vec1%arrv=v1f*vec1%arrv+v2f*vec2%arrv
        IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri))
     $    vec1%arri=v1f*vec1%arri+v2f*vec2%arri
        IF (ASSOCIATED(vec1%arrtmp).AND.ASSOCIATED(vec2%arrtmp))
     $    vec1%arrtmp=v1f*vec1%arrtmp+v2f*vec2%arrtmp
      ELSE IF (ASSOCIATED(vec1%arri).AND.ASSOCIATED(vec2%arri)) THEN
        vec1%arri=v1f*vec1%arri+v2f*vec2%arri
      ELSE
        CALL nim_stop
     $    ('Vector_3D_add_vec3: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_add_vec3
c-----------------------------------------------------------------------
c     subprogram 50. vector_3D_mult_rsc.
c     multiply a 3D vector structure by a real scalar.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_mult_rsc(vec,rsc)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      REAL(r8), INTENT(IN) :: rsc

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=rsc*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=rsc*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=rsc*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=rsc*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=rsc*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=rsc*vec%arri
      ELSE
        CALL nim_stop
     $    ('Vector_3D_mult_rsc: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_mult_rsc
c-----------------------------------------------------------------------
c     subprogram 51. vector_3D_mult_int.
c     multiply a 3D vector structure by an integer.
c-----------------------------------------------------------------------
      SUBROUTINE vector_3D_mult_int(vec,int)

      TYPE(vector_3D_type), INTENT(INOUT) :: vec
      INTEGER(i4), INTENT(IN) :: int

c-----------------------------------------------------------------------
c     warning:  there are no checks on compatibility.
c-----------------------------------------------------------------------
      IF (ASSOCIATED(vec%arr)) THEN
        vec%arr=int*vec%arr
        IF (ASSOCIATED(vec%arrh)) vec%arrh=int*vec%arrh
        IF (ASSOCIATED(vec%arrv)) vec%arrv=int*vec%arrv
        IF (ASSOCIATED(vec%arri)) vec%arri=int*vec%arri
        IF (ASSOCIATED(vec%arrtmp)) vec%arrtmp=int*vec%arrtmp
      ELSE IF (ASSOCIATED(vec%arri)) THEN
        vec%arri=int*vec%arri
      ELSE
        CALL nim_stop
     $    ('Vector_3D_mult_int: vector arrays not associated.')
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE vector_3D_mult_int
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE rvector_3D_type_mod


c-----------------------------------------------------------------------
c     the generic module contains overloaded assignment and operator
c     and interfaces.
c-----------------------------------------------------------------------
      MODULE vector_type_mod
      USE rvector_type_mod
      USE cvector_type_mod
      USE cvector_2D_type_mod
      USE rvector_3D_type_mod
      IMPLICIT NONE

      INTERFACE ASSIGNMENT (=)
        MODULE PROCEDURE vector_assign_rsc,vector_assign_csc,
     $    vector_assign_int,vector_assign_vec,cvector_assign_rsc,
     $    cvector_assign_csc,cvector_assign_int,cvector_assign_cvec,
     $    cvector_2D_assign_rsc,cvector_2D_assign_csc,
     $    cvector_2D_assign_int,cvector_2D_assign_cvec2,
     $    vector_3D_assign_rsc,vector_3D_assign_csc,
     $    vector_3D_assign_int,vector_3D_assign_vec3
      END INTERFACE

      INTERFACE vector_add
        MODULE PROCEDURE vector_add_vec,cvector_add_cvec,
     $    cvector_2D_add_cvec2,cvector_addc_cvec,
     $    cvector_2D_addc_cvec2,vector_3D_add_vec,vector_3D_add_vec3
      END INTERFACE

      INTERFACE vector_mult
        MODULE PROCEDURE vector_mult_rsc,vector_mult_int,
     $    cvector_mult_rsc,cvector_mult_int,
     $    cvector_2D_mult_rsc,cvector_2D_mult_int,
     $    cvector_mult_csc,cvector_2D_mult_csc,vector_3D_mult_rsc,
     $    vector_3D_mult_int
      END INTERFACE

      INTERFACE vector_type_alloc
        MODULE PROCEDURE vector_rtype_alloc,vector_ctype_alloc,
     $    vector_2D_ctype_alloc,vector_3D_rtype_alloc
      END INTERFACE

      INTERFACE vector_type_dealloc
        MODULE PROCEDURE vector_rtype_dealloc,vector_ctype_dealloc,
     $    vector_2D_ctype_dealloc,vector_3D_rtype_dealloc
      END INTERFACE
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE vector_type_mod
