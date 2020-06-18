c-----------------------------------------------------------------------
c     file seam_storage.f
c     module containing the seam structures used for communications
c     among grid blocks.
c-----------------------------------------------------------------------
      MODULE seam_storage_mod
      USE local
      USE edge_type_mod
      IMPLICIT NONE

      TYPE(edge_type), POINTER :: seam0
      TYPE(edge_type), DIMENSION(:), POINTER :: seam
      INTEGER(i4), DIMENSION(:), POINTER :: exblock_list
      INTEGER(i4) :: max_imags
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE seam_storage_mod
