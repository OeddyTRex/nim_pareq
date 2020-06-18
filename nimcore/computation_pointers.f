c-----------------------------------------------------------------------
c     file computation_pointer.f
c     module containing pointers for use during finite element
c     computations.  the rhs, sln, and cell_rhs pointers are
c     associated with (otherwise) unnamed memory locations.
c-----------------------------------------------------------------------
      MODULE computation_pointers
      USE local
      USE vector_type_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     working array structures.
c-----------------------------------------------------------------------
      TYPE(vector_type), DIMENSION(:), POINTER :: rhs,cell_rhs,sln,
     $                   vectr,lump_mass,lump_summed
      TYPE(cvector_type), DIMENSION(:), POINTER :: crhs,cell_crhs,cvecn
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: csln,cvectr
c-----------------------------------------------------------------------
c     close module
c-----------------------------------------------------------------------
      END MODULE computation_pointers
