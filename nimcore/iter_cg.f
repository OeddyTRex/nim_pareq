c-----------------------------------------------------------------------
c     file iter_cg.f
c     module that serves as a switchyard for different versions of
c     the 2D conjugate gradient routines that are for real or imaginary
c     linear systems.
c
c     the external subroutines that had been in this file have been
c     moved to iter_externals.f for modularity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     0. module declarations for iter_cg:  the entry routines into the
c     cg suites are determined by the type of structures (real or 
c     complex) that are passed.
c-----------------------------------------------------------------------
      MODULE iter_cg
      USE iter_cg_real
      USE iter_cg_comp

      INTERFACE iter_cg_2d_solve
        MODULE PROCEDURE iter_cg_r2d_solve,iter_cg_c2d_solve
      END INTERFACE

      INTERFACE iter_fac_alloc
        MODULE PROCEDURE iter_fac_alloc_real,iter_fac_alloc_comp
      END INTERFACE

      INTERFACE iter_fac_dealloc
        MODULE PROCEDURE iter_fac_dealloc_real,iter_fac_dealloc_comp
      END INTERFACE

      INTERFACE iter_factor
        MODULE PROCEDURE iter_factor_real,iter_factor_comp
      END INTERFACE
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_cg
