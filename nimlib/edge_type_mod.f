c-----------------------------------------------------------------------
c     file edge_type_mod.f
c     module containing edge_type and vertex_type definitions.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE edge_type_mod
      USE local
      IMPLICIT NONE

      TYPE :: vertex_type
        INTEGER(i4), DIMENSION(2) :: intxy
        INTEGER(i4) :: nimage
        INTEGER(i4), DIMENSION(:), POINTER :: order
        INTEGER(i4), DIMENSION(:,:), POINTER :: ptr,ptr2
        REAL(r8), DIMENSION(:), POINTER :: seam_in,seam_out,tang,norm,
     $            seam_save
        REAL(r8), DIMENSION(:,:), POINTER :: seam_hold
        COMPLEX(r8), DIMENSION(:), POINTER :: seam_cin,seam_cout,
     $            seam_csave
        COMPLEX(r8), DIMENSION(:,:), POINTER :: seam_chold
        REAL(r8) :: ave_factor,ave_factor_pre,rgeom
      END TYPE vertex_type

      TYPE :: segment_type
        INTEGER(i4), DIMENSION(2) :: intxyn,intxyp,intxys,ptr
        INTEGER(i4) :: load_dir
        LOGICAL :: h_side
        REAL(r8), DIMENSION(:), POINTER :: seam_in,seam_out,
     $            seam_save
        REAL(r8), DIMENSION(:,:), POINTER :: tang,norm
        REAL(r8), DIMENSION(:,:,:), POINTER :: seam_mat_in,seam_mat_out,
     $            seam_mat_save
        COMPLEX(r8), DIMENSION(:), POINTER :: seam_cin,
     $            seam_cout,seam_csave
        COMPLEX(r8), DIMENSION(:,:,:), POINTER :: seam_mat_cin,
     $            seam_mat_cout,seam_mat_csave
        REAL(r8) :: ave_factor,ave_factor_pre
      END TYPE segment_type

      TYPE :: edge_type
        CHARACTER(64) :: name
        INTEGER(i4) :: id
        INTEGER(i4) :: nvert
        TYPE(vertex_type), DIMENSION(:), POINTER :: vertex
        TYPE(segment_type), DIMENSION(:), POINTER :: segment
        LOGICAL, DIMENSION(:), POINTER :: expoint
        LOGICAL, DIMENSION(:), POINTER :: excorner
        LOGICAL, DIMENSION(:), POINTER :: r0point
      END TYPE edge_type
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE edge_type_mod
