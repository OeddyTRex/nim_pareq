c-----------------------------------------------------------------------
c     file matrix_storage.f
c     contains a module for storage matrices and their approximate
c     factors which are used for preconditioning.
c     this is a stripped version of the kernel's matrix storage module.
c-----------------------------------------------------------------------

      MODULE matrix_storage_mod
      USE local
      USE matrix_type_mod
      USE factor_type_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     storage for matrices used for each equation.
c-----------------------------------------------------------------------
      TYPE(global_matrix_type) :: mass_mat,delstar_mat,bfield_mat,
     $  delstlj_mat,fladv_mat
c-----------------------------------------------------------------------
c     storage for matrix factors used for each equation.
c-----------------------------------------------------------------------
      TYPE(matrix_factor_type) :: mass_fac,delstar_fac,bfield_fac,
     $  delstlj_fac,fladv_fac
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matrix_storage_mod
