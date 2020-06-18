c-----------------------------------------------------------------------
c     file cell_type_mod.f
c     module containing cell_type definition.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE cell_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: location_type
        REAL(r8), DIMENSION(2) :: point 
      END TYPE location_type

      TYPE :: face_type
        TYPE(cell_type), POINTER :: p
      END TYPE face_type
      TYPE :: cell_type
        CHARACTER(8) ::  representation		! see met_spl: 'liner', 'pcnst'
						! or 'cubic'
        INTEGER(i4) :: 	id			! Unique global identifier
        INTEGER(i4) :: nodes			! 3 is a triangle, 4 is a quad
        INTEGER(i4) :: ib,icell			! mapping to original cell
        INTEGER(i4), DIMENSION(:,:), POINTER :: p ! local block cell corners
c              ?%p(1,1)  is ix of point 1 and ?%p(2,1)  is iy of point 1
        REAL(r8), DIMENSION(:,:), POINTER :: pp   ! physical positions
        TYPE(cell_type), POINTER :: next
        TYPE(face_type), DIMENSION(:), POINTER :: face	
					! cell associated with that face
        LOGICAL :: degenerate			! degenerate quad?
        REAL(r8) :: area	! Cell area, plus sign indicates sequence.
      END TYPE cell_type
c-----------------------------------------------------------------------
c    				 Useful info for speeding searches
      INTEGER(i4) :: search_idim,search_jdim
      REAL(r8) :: search_rmin		! minimum R of all cells
      REAL(r8) :: search_rmax		! maximum R of all cells
      REAL(r8) :: search_zmin		! minimum Z of all cells
      REAL(r8) :: search_zmax		! maximum Z of all cells
      TYPE(cell_type),DIMENSION(:,:), POINTER :: search_map
c     isearch = NINT(search_idim*
c    &              (r-search_rmin)/(search_rmax-search_rmin))
c     jsearch = NINT(search_jdim*
c    &              (z-search_zmin)/(search_zmax-search_zmin))
c       Usage item => search_map(isearch,jsearch)%next
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE cell_type_mod 
