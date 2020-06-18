c-----------------------------------------------------------------------
c     file threed_dot_mgt.f:  contains external subprograms that compute
c     the matrix-vector dot products needed for 'matrix-free' iteration
c     in the iter_ky_c3d_solve.
c
c     writing the calls here simplifies the coding in iter_krylov_c3d.f.
c     the subroutine interfaces all have the common form of:

c     INTERFACE
c       SUBROUTINE dot_routine(oper,prod,bc_oper)
c       USE vector_type_mod
c       USE local
c       TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
c       LOGICAL, INTENT(IN) :: bc_oper
c       END SUBROUTINE dot_routine
c     END INTERFACE
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. threed_b_3deta.
c     2. threed_b_3dnsym.
c     3. threed_v_3dn.
c     4. threed_t_aniso.
c     5. threed_n_3dnsym.
c-----------------------------------------------------------------------
c     subprogram 1. threed_b_3deta
c     find the matrix-vector product for the magnetic advance with 3D
c     resistivity, without implicit advection.
c-----------------------------------------------------------------------
      SUBROUTINE threed_b_3deta(oper,prod,bc_oper)
      USE local
      USE input
      USE integrands
      USE boundary
      USE surface_ints
      USE finite_element_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      USE matrix_storage_mod
      USE edge
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,im,is,iqv,iqs,iv
c-----------------------------------------------------------------------
c     project-out the degrees of freedom that would violate the
c     essential conditions.  save this information in the seam_csave
c     parts of the seam structure.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL edge_zero_csave(seam(ibl))
      ENDDO
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(oper(ibe),seam(ibe),'3vn',3_i4,
     $                     svess='yes')
      ENDDO
      CALL regular_pre_feop(oper,'cyl_vec',3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     find the matrix-vector product via FE computations.
c-----------------------------------------------------------------------
      CALL get_rhs(b_3deta_dot,prod,dirichlet_rhs,'3vn','cyl_vec',
     $             .false.,no_surf_int)
c-----------------------------------------------------------------------
c     add the scaled degrees of freedom into the prod structure, and
c     restore the original factors in the oper structure.
c-----------------------------------------------------------------------
      CALL regular_zero_phi(oper,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL edge_add_csave(oper(ibl),3_i4,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          iqv=1
          iqs=1
          DO im=1,nmodes
            seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+2)=
     $        seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+2)*
     $        bhmhd_cmat(im)%diag_scale
            iqv=iqv+3
            DO is=1,poly_degree-1
              seam(ibl)%segment(iv)%seam_csave(iqs:iqs+2)=
     $          seam(ibl)%segment(iv)%seam_csave(iqs:iqs+2)*
     $          bhmhd_cmat(im)%diag_scale
              iqs=iqs+3
            ENDDO
          ENDDO
        ENDDO
        CALL edge_add_csave(prod(ibl),3_i4,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_b_3deta
c-----------------------------------------------------------------------
c     subprogram 2. threed_b_3dnsym
c     find the matrix-vector product for the 3D magnetic advance with
c     asymmetric operators from implicit advection or two-fluid effects.
c-----------------------------------------------------------------------
      SUBROUTINE threed_b_3dnsym(oper,prod,bc_oper)
      USE local
      USE input
      USE integrands
      USE boundary
      USE surface_ints
      USE finite_element_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      USE matrix_storage_mod
      USE edge
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,nqvec,ibe,im,is,iqv,iqs,iv,nqm
c-----------------------------------------------------------------------
c     choose the number of vector components in the system.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
        nqvec=6
      ELSE
        nqvec=3
      ENDIF
      nqm=nqvec-1
c-----------------------------------------------------------------------
c     project-out the degrees of freedom that would violate the
c     essential conditions.  save this information in the seam_csave
c     parts of the seam structure.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL edge_zero_csave(seam(ibl))
      ENDDO
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(oper(ibe),seam(ibe),'3vn sd sd sd',nqvec,
     $                     svess='yes')
      ENDDO
      CALL regular_pre_feop(oper,'cyl_vec',nqvec,nmodes,nindex)
c-----------------------------------------------------------------------
c     find the matrix-vector product via FE computations.
c-----------------------------------------------------------------------
      CALL get_rhs(b_hmhd_dot,prod,dirichlet_rhs,'3vn sd sd sd',
     $             'cyl_vec',.false.,no_surf_int)
c-----------------------------------------------------------------------
c     add the scaled degrees of freedom into the prod structure, and
c     restore the original factors in the oper structure.
c-----------------------------------------------------------------------
      CALL regular_zero_phi(oper,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL edge_add_csave(oper(ibl),nqvec,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          iqv=1
          iqs=1
          DO im=1,nmodes
            seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+nqm)=
     $        seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+nqm)*
     $        bhmhd_cmat(im)%diag_scale
            iqv=iqv+nqvec
            DO is=1,poly_degree-1
              seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nqm)=
     $          seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nqm)*
     $          bhmhd_cmat(im)%diag_scale
              iqs=iqs+nqvec
            ENDDO
          ENDDO
        ENDDO
        CALL edge_add_csave(prod(ibl),nqvec,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_b_3dnsym
c-----------------------------------------------------------------------
c     subprogram 3. threed_v_3dn
c     find the matrix-vector product for the flow-velocity advance with
c     spatially varying density and asymmetric contributions.
c-----------------------------------------------------------------------
      SUBROUTINE threed_v_3dn(oper,prod,bc_oper)
      USE local
      USE input
      USE integrands
      USE boundary
      USE surface_ints
      USE finite_element_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      USE matrix_storage_mod
      USE edge
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,im,is,iqv,iqs,iv
      CHARACTER(8) :: bc_flag
c-----------------------------------------------------------------------
c     set the boundary option.
c-----------------------------------------------------------------------
      IF (flow_bc=='free-slip') THEN
        bc_flag='3vn'
      ELSE
        bc_flag='all'
      ENDIF
c-----------------------------------------------------------------------
c     project-out the degrees of freedom that would violate the
c     essential conditions.  save this information in the seam_csave
c     parts of the seam structure.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL edge_zero_csave(seam(ibl))
      ENDDO
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(oper(ibe),seam(ibe),bc_flag,3_i4,svess='yes')
      ENDDO
      CALL regular_pre_feop(oper,'cyl_vec',3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     use the 3D semi-implicit operator if it is specified.
c-----------------------------------------------------------------------
      IF (siop_type=='3D') THEN
        CALL get_rhs(v_3dsi_dot,prod,dirichlet_rhs,bc_flag,'cyl_vec',
     $               .false.,no_surf_int)
      ELSE
        CALL get_rhs(v_aniso_dot,prod,dirichlet_rhs,bc_flag,'cyl_vec',
     $               .false.,no_surf_int)
      ENDIF
c-----------------------------------------------------------------------
c     add the scaled degrees of freedom into the prod structure, and
c     restore the original factors in the oper structure.
c-----------------------------------------------------------------------
      CALL regular_zero_phi(oper,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL edge_add_csave(oper(ibl),3_i4,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          iqv=1
          iqs=1
          DO im=1,nmodes
            seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+2)=
     $        seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+2)*
     $        vmhd_cmat(im)%diag_scale
            iqv=iqv+3
            DO is=1,poly_degree-1
              seam(ibl)%segment(iv)%seam_csave(iqs:iqs+2)=
     $          seam(ibl)%segment(iv)%seam_csave(iqs:iqs+2)*
     $          vmhd_cmat(im)%diag_scale
              iqs=iqs+3
            ENDDO
          ENDDO
        ENDDO
        CALL edge_add_csave(prod(ibl),3_i4,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_v_3dn
c-----------------------------------------------------------------------
c     subprogram 4. threed_t_aniso
c     find the matrix-vector product for the temperature advance with
c     spatially varying density, conductivity tensor, or implicit
c     convection.
c-----------------------------------------------------------------------
      SUBROUTINE threed_t_aniso(oper,prod,bc_oper)
      USE local
      USE input
      USE integrands
      USE boundary
      USE surface_ints
      USE finite_element_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      USE matrix_storage_mod
      USE edge
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,im,is,iqs,iv,nsm
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: t_mat
c-----------------------------------------------------------------------
c     the boundary condition for temperature depends on p_model and
c     on insulate, consistent with adv_t_aniso.
c-----------------------------------------------------------------------
      IF (p_model=='adiabat'.OR.p_model=='isothermal'.OR.insulate) THEN
        CALL get_rhs(t_aniso_dot,prod,no_rhs_bc,' ','scalar',
     $               .false.,no_surf_int)
      ELSE
c-----------------------------------------------------------------------
c       project-out the degrees of freedom that would violate the
c       essential conditions.  save this information in the seam_csave
c       parts of the seam structure.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL edge_zero_csave(seam(ibl))
        ENDDO
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(oper(ibe),seam(ibe),'sd',1_i4,svess='yes')
        ENDDO
        CALL regular_pre_feop(oper,'scalar',1_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c       find the matrix-vector product via FE computations.
c-----------------------------------------------------------------------
        CALL get_rhs(t_aniso_dot,prod,dirichlet_rhs,'sd','scalar',
     $               .false.,no_surf_int)
c-----------------------------------------------------------------------
c       add the scaled degrees of freedom into the prod structure, and
c       restore the original factors in the oper structure.
c-----------------------------------------------------------------------
        IF (integrand_flag(2:2)=='e'.OR.integrand_flag(6:6)=='e') THEN
          t_mat=>te_cmat
        ELSE
          t_mat=>ti_cmat
        ENDIF
        nsm=poly_degree-2
        DO ibl=1,nbl
          CALL edge_add_csave(oper(ibl),1_i4,1_i4,nmodes,
     $                        poly_degree-1_i4,seam(ibl))
          DO iv=1,seam(ibl)%nvert
            iqs=1
            DO im=1,nmodes
              seam(ibl)%vertex(iv)%seam_csave(im)=
     $          seam(ibl)%vertex(iv)%seam_csave(im)*
     $          t_mat(im)%diag_scale
              seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nsm)=
     $          seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nsm)*
     $          t_mat(im)%diag_scale
              iqs=iqs+poly_degree-1
            ENDDO
          ENDDO
          CALL edge_add_csave(prod(ibl),1_i4,1_i4,nmodes,
     $                        poly_degree-1_i4,seam(ibl))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_t_aniso
c-----------------------------------------------------------------------
c     subprogram 5. threed_n_3dnsym
c     find the matrix-vector product for the density advance, possibly
c     with the auxiliary scalar for hyper-diffusivity.
c-----------------------------------------------------------------------
      SUBROUTINE threed_n_3dnsym(oper,prod,bc_oper)
      USE local
      USE input
      USE integrands
      USE boundary
      USE surface_ints
      USE finite_element_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE mpi_nim
      USE pardata
      USE matrix_storage_mod
      USE edge
      IMPLICIT NONE

      TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,im,is,iqv,iqs,iv,nq,nqm
c-----------------------------------------------------------------------
c     the boundary condition for n is indicated by nd_bc, and
c     use of the auxiliary field is determined by nd_hypd.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        nq=2
      ELSE
        nq=1
      ENDIF
      nqm=nq-1
c-----------------------------------------------------------------------
c     project-out the degrees of freedom that would violate the
c     essential conditions.  save this information in the seam_csave
c     parts of the seam structure.
c-----------------------------------------------------------------------
      IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        DO ibl=1,nbl
          CALL edge_zero_csave(seam(ibl))
        ENDDO
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(oper(ibe),seam(ibe),'sd sd',nq,svess='yes')
        ENDDO
        CALL regular_pre_feop(oper,'scalar',nq,nmodes,nindex)
c-----------------------------------------------------------------------
c       find the matrix-vector product via FE computations.
c-----------------------------------------------------------------------
        CALL get_rhs(cont_dot,prod,dirichlet_rhs,'sd sd','scalar',
     $               .false.,no_surf_int)
c-----------------------------------------------------------------------
c       add the scaled degrees of freedom into the prod structure, and
c       restore the original factors in the oper structure.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL edge_add_csave(oper(ibl),nq,1_i4,nmodes,
     $                        poly_degree-1_i4,seam(ibl))
          DO iv=1,seam(ibl)%nvert
            iqv=1
            iqs=1
            DO im=1,nmodes
              seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+nqm)=
     $          seam(ibl)%vertex(iv)%seam_csave(iqv:iqv+nqm)*
     $          nd_cmat(im)%diag_scale
              iqv=iqv+nq
              DO is=1,poly_degree-1
                seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nqm)=
     $            seam(ibl)%segment(iv)%seam_csave(iqs:iqs+nqm)*
     $            nd_cmat(im)%diag_scale
                iqs=iqs+nq
              ENDDO
            ENDDO
          ENDDO
          CALL edge_add_csave(prod(ibl),nq,1_i4,nmodes,
     $                        poly_degree-1_i4,seam(ibl))
        ENDDO
c-----------------------------------------------------------------------
c     if there are no essential conditions on density, just call the
c     FE routines.
c-----------------------------------------------------------------------
      ELSE
        CALL get_rhs(cont_dot,prod,no_rhs_bc,' ','scalar',
     $               .false.,no_surf_int)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE threed_n_3dnsym
