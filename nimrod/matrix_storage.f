c-----------------------------------------------------------------------
c     file matrix_storage.f
c     contains a module for storage matrices and their approximate
c     factors which are used for preconditioning.
c-----------------------------------------------------------------------

      MODULE matrix_storage_mod
      USE local
      USE matrix_type_mod
      USE factor_type_mod
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     storage for matrices used for each equation.
c-----------------------------------------------------------------------
      TYPE(global_matrix_type), DIMENSION(:), POINTER ::
     $  bmhd_mat,bhall_mat,bhalp_mat,divb_mat,visc_mat,
     $  vmhd_mat,mass_mat,j_mat,teiso_mat,tiiso_mat,ndiso_mat,sproj_mat,
     $  hypv_expmat
      TYPE(complex_matrix_type), DIMENSION(:), POINTER ::
     $  bhmhd_cmat,vmhd_cmat,vmhdsi_cmat,te_cmat,ti_cmat,nd_cmat,
     $  bhyp_cmat,hypv_cmat
      TYPE(complex_matrix_type), DIMENSION(:,:), POINTER ::
     $  bhmhd_cplus,bhmhd_cminus,ti_cminus,ti_cplus,te_cminus,te_cplus
c-----------------------------------------------------------------------
c     storage for matrix factors used for each equation.
c-----------------------------------------------------------------------
      TYPE(matrix_factor_type), DIMENSION(:), POINTER ::
     $  bmhd_fac,bhall_fac,bhalp_fac,divb_fac,visc_fac,
     $  vmhd_fac,mass_fac,j_fac,teiso_fac,tiiso_fac,ndiso_fac,sproj_fac
      TYPE(complex_factor_type), DIMENSION(:), POINTER ::
     $  bhmhd_cfac,vmhd_cfac,te_cfac,ti_cfac,nd_cfac,extra_cfac,
     $  bhyp_cfac,hypv_cfac
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE matrix_storage_mod
