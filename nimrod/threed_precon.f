c-----------------------------------------------------------------------
c     file threed_precon.f:  contains external subprograms that apply
c     preconditioning operations for 'matrix-free' iteration in the
c     iter_ky_c3d_solve.
c
c     writing the calls here enhances flexbility in the operations that
c     can be used for preconditioning, though it adds an extra layer
c     when calling the standard Fourier-diagonal operations through
c     iter_pre_c3dto2d.
c
c     the subroutine interfaces all have the common form of:
c
c     INTERFACE
c       SUBROUTINE pre_routine(resd,zeed,iiter,flex)
c       USE vector_type_mod
c       USE local
c       TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
c       TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
c       INTEGER(i4), INTENT(IN) :: iiter
c       LOGICAL, INTENT(IN) :: flex
c       END SUBROUTINE pre_routine
c     END INTERFACE
c-----------------------------------------------------------------------
c     code organization (intentionally mirrors threed_dot_mgt).
c-----------------------------------------------------------------------
c     1. threed_preb_3deta.
c     2. threed_preb_3dnsym.
c     3. threed_prev_3dn.
c     4. threed_pret_aniso.
c     5. threed_pren_3dnsym.
c-----------------------------------------------------------------------
c     subprogram 1. threed_preb_3deta
c     apply the preconditioner for the magnetic advance with 3D
c     resistivity, without implicit advection.
c-----------------------------------------------------------------------
      SUBROUTINE threed_preb_3deta(res,zee,itcount,flex)
      USE local
      USE vector_type_mod
      USE input
      USE fields
      USE global
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: res
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zee
      INTEGER(i4), INTENT(IN) :: itcount
      LOGICAL, INTENT(IN) :: flex

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     copy res here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        zee(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     call the standard diagonal-in-Fourier operations with the
c     appropriate matrix and factor structures.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d(bhmhd_cmat,bhmhd_cfac,zee,poly_degree,
     $                      3_i4,0_i4,0_i4,nrbl,nbl,nmodes,bmhd_solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_preb_3deta
c-----------------------------------------------------------------------
c     subprogram 2. threed_preb_3dnsym
c     apply the preconditioner for the 3D magnetic advance with
c     asymmetric operators from implicit advection or two-fluid effects.
c-----------------------------------------------------------------------
      SUBROUTINE threed_preb_3dnsym(res,zee,itcount,flex)
      USE local
      USE vector_type_mod
      USE input
      USE fields
      USE global
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: res
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zee
      INTEGER(i4), INTENT(IN) :: itcount
      LOGICAL, INTENT(IN) :: flex

      INTEGER(i4) :: nq,nqdsc,nbdsc,ibl
c-----------------------------------------------------------------------
c     copy res here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        zee(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     set the quantity and basis dimensions.  this must match the
c     adv_b_3dnsym management routine.
c-----------------------------------------------------------------------
      IF (poly_divb>=0.AND.nrbl>0) THEN
        nqdsc=rb(1)%auxb%nqty
        nbdsc=rb(1)%auxb%n_int
      ELSE
        nqdsc=0_i4
        nbdsc=0_i4
      ENDIF
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
        nq=6
      ELSE
        nq=3
      ENDIF
c-----------------------------------------------------------------------
c     call the standard diagonal-in-Fourier operations with the
c     appropriate matrix and factor structures.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d(bhmhd_cmat,bhmhd_cfac,zee,poly_degree,
     $                      nq,nbdsc,nqdsc,nrbl,nbl,nmodes,bmhd_solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_preb_3dnsym
c-----------------------------------------------------------------------
c     subprogram 3. threed_prev_3dn
c     apply the preconditioner for the flow-velocity advance with
c     spatially varying density and asymmetric contributions.
c-----------------------------------------------------------------------
      SUBROUTINE threed_prev_3dn(res,zee,itcount,flex)
      USE local
      USE vector_type_mod
      USE input
      USE fields
      USE global
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: res
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zee
      INTEGER(i4), INTENT(IN) :: itcount
      LOGICAL, INTENT(IN) :: flex

      INTEGER(i4) :: nqdsc,nbdsc,ibl
c-----------------------------------------------------------------------
c     copy res here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        zee(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     set the quantity and basis dimensions.  this must match the
c     adv_v_3dn management routine.
c-----------------------------------------------------------------------
      IF (poly_divv>=0.AND.nrbl>0) THEN
        nqdsc=rb(1)%auxv%nqty
        nbdsc=rb(1)%auxv%n_int
      ELSE
        nqdsc=0_i4
        nbdsc=0_i4
      ENDIF
c-----------------------------------------------------------------------
c     call the standard diagonal-in-Fourier operations with the
c     appropriate matrix and factor structures.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d(vmhd_cmat,vmhd_cfac,zee,poly_degree,3_i4,
     $                      nbdsc,nqdsc,nrbl,nbl,nmodes,vmhd_solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_prev_3dn
c-----------------------------------------------------------------------
c     subprogram 4. threed_pret_aniso
c     apply the preconditioner for the temperature advance with
c     spatially varying density, conductivity tensor, or implicit
c     convection.
c-----------------------------------------------------------------------
      SUBROUTINE threed_pret_aniso(res,zee,itcount,flex)
      USE local
      USE vector_type_mod
      USE input
      USE fields
      USE global
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: res
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zee
      INTEGER(i4), INTENT(IN) :: itcount
      LOGICAL, INTENT(IN) :: flex

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     copy res here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        zee(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     call the standard diagonal-in-Fourier operations with the
c     appropriate matrix and factor structures.  this depends on
c     whether the system is for electron temperature or ion temperature.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:6)=="all ti") THEN
        CALL iter_pre_c3dto2d(ti_cmat,ti_cfac,zee,poly_degree,1_i4,
     $                        0_i4,0_i4,nrbl,nbl,nmodes,temp_solver)
      ELSE
        CALL iter_pre_c3dto2d(te_cmat,te_cfac,zee,poly_degree,1_i4,
     $                        0_i4,0_i4,nrbl,nbl,nmodes,temp_solver)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_pret_aniso
c-----------------------------------------------------------------------
c     subprogram 5. threed_pren_3dnsym
c     apply the preconditioner for the density advance, possibly
c     with the auxiliary scalar for hyper-diffusivity.
c-----------------------------------------------------------------------
      SUBROUTINE threed_pren_3dnsym(res,zee,itcount,flex)
      USE local
      USE vector_type_mod
      USE input
      USE fields
      USE global
      USE matrix_storage_mod
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: res
      TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zee
      INTEGER(i4), INTENT(IN) :: itcount
      LOGICAL, INTENT(IN) :: flex

      INTEGER(i4) :: nq,ibl
c-----------------------------------------------------------------------
c     copy res here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        zee(ibl)=res(ibl)
      ENDDO
c-----------------------------------------------------------------------
c     set the quantity dimension.  this must match the
c     adv_nd_3dnsym management routine.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        nq=2
      ELSE
        nq=1
      ENDIF
c-----------------------------------------------------------------------
c     call the standard diagonal-in-Fourier operations with the
c     appropriate matrix and factor structures.
c-----------------------------------------------------------------------
      CALL iter_pre_c3dto2d(nd_cmat,nd_cfac,zee,poly_degree,
     $                      nq,0_i4,0_i4,nrbl,nbl,nmodes,solver)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_pren_3dnsym
