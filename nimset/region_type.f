c-----------------------------------------------------------------------
c     file region_type.f
c     module that defines a region data structure for holding separate
c     sets of blocks and seams.
c-----------------------------------------------------------------------
      MODULE region_type_mod
      USE local
      USE fields
      USE edge_type_mod
      USE vector_type_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     type definition for each region
c-----------------------------------------------------------------------
      TYPE :: region_type
        INTEGER(i4) :: nbl_total,nrbl_total
        INTEGER(i4) :: nbl,nrbl,irst,itst
        TYPE(rblock_type), DIMENSION(:), POINTER :: rb
        TYPE(tblock_type), DIMENSION(:), POINTER :: tb
        TYPE(edge_type), POINTER :: seam0
        TYPE(edge_type), DIMENSION(:), POINTER :: seam
      END TYPE region_type
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE region_type_mod
