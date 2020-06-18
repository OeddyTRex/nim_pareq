c-----------------------------------------------------------------------
c     file nimeq_dot_mgt.f:  contains external subprograms that compute
c     the matrix-vector dot products needed for 'matrix-free' iteration.
c
c     the subroutine interfaces all have the common form of:

c     INTERFACE
c       SUBROUTINE dot_routine(oper,prod,bc_oper)
c       USE vector_type_mod
c       USE local
c       TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
c       LOGICAL, INTENT(IN) :: bc_oper
c       END SUBROUTINE dot_routine
c     END INTERFACE
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. eqdot_delstlj.
c-----------------------------------------------------------------------
c     subprogram 1. eqdot_delstlj
c     this dot-product management routine computes the result of
c     applying the del-star operator on a scalar (lambda=psi/R**2) and
c     applying the mass matrix on a second scalar (mu0*J_phi/R), as
c     needed for the free-boundary Grad-Shafranov solves.
c-----------------------------------------------------------------------
      SUBROUTINE eqdot_delstlj(oper,prod,bc_oper)
      USE local
      USE input
      USE nimeq_ints
      USE boundary
      USE regularity
      USE global
      USE fields
      USE computation_pointers
      USE rblock
      USE tblock
      USE nimeq_free
      USE seam_storage_mod
      USE edge
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,nmodes_save,ncoil_save,iv,ivp,ix,iy,is
      REAL(r8) :: sflmin,sflmax
c-----------------------------------------------------------------------
c     the integrand computation uses ja_eq(2:3) as the operand, so
c     transfer the passed data to ja_eq.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_assignq_vec(ja_eq(ibl),oper(ibl),nqty=2_i4,
     $                          nstart1=2_i4,nstart2=1_i4)
        CALL rblock_qp_update(rb(ibl)%ja_eq,rb(ibl)%qja_eq,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_assignq_vec(ja_eq(ibl),oper(ibl),nqty=2_i4,
     $                          nstart1=2_i4,nstart2=1_i4)
        CALL tblock_qp_update(tb(ibl)%ja_eq,tb(ibl)%qja_eq,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     the following evaluates the surface flux from the existing
c     internal currents represented in the direction vector and
c     without external current.  here, sln is just temporary storage.
c-----------------------------------------------------------------------
      ncoil_save=ncoil_tot
      ncoil_tot=0
      DO ibl=1,nbl
        sln(ibl)=0._r8
      ENDDO
      CALL nimeq_free_eval(sln,1_i4,poly_degree,1._r8,
     $                     sflmax,sflmin)
      ncoil_tot=ncoil_save
c-----------------------------------------------------------------------
c     the boundary condition is applied to the operand only during
c     initialization of the Krylov solve, as indicated by bc_oper.
c
c     use copies to the complex crhs data structure to use the standard
c     boundary and regularity routines.
c-----------------------------------------------------------------------
      IF (bc_oper) THEN
        nmodes_save=nmodes
        nmodes=1
        DO ibl=1,SIZE(r0block_list)
          ibe=r0block_list(ibl)
          CALL cvector_assign_vec(crhs(ibe),ja_eq(ibe),'real',1_i4,
     $                            nqty=2_i4,nstart1=1_i4,nstart2=2_i4)
          CALL regular_vec(crhs(ibe),seam(ibe),'all',2_i4,
     $                     nmodes,nindex)
          CALL vector_assign_cvec(ja_eq(ibe),crhs(ibe),'real',1_i4,
     $                            nqty=2_i4,nstart1=2_i4,nstart2=1_i4)
        ENDDO
        nmodes=nmodes_save
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the matrix vector product directly through the
c     finite-element computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL rblock_get_rhs(rb(ibl),prod(ibl),delstarlj_dot,2_i4)
        CALL edge_load_arr(prod(ibl),2_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_get_rhs(tb(ibl),prod(ibl),delstarlj_dot,2_i4)
        CALL edge_load_arr(prod(ibl),2_i4,0_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     sum contributions over block boundaries.
c-----------------------------------------------------------------------
      CALL edge_network(2_i4,0_i4,poly_degree-1_i4,.false.)
      DO ibl=1,nrbl
        CALL edge_unload_arr(prod(ibl),2_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL edge_unload_arr(prod(ibl),2_i4,0_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     apply regularity conditions to the product.
c-----------------------------------------------------------------------
      nmodes_save=nmodes
      nmodes=1
      DO ibl=1,SIZE(r0block_list)
        ibe=r0block_list(ibl)
        CALL cvector_assign_vec(crhs(ibe),prod(ibe),'real',1_i4)
        CALL regular_vec(crhs(ibe),seam(ibe),'all',2_i4,
     $                   nmodes,nindex)
        CALL vector_assign_cvec(prod(ibe),crhs(ibe),'real',1_i4)
      ENDDO
      nmodes=nmodes_save
c-----------------------------------------------------------------------
c     replace surface-lambda values with the result based on internal
c     currents.  each surface-lamda row effectively amounts to
c     lambda_i-SUM(mat_ij*I_j,j-interior)=SUM(mat_ij*I_j,j-exterior) 
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        ivp=seam(ibe)%nvert
        DO iv=1,seam(ibe)%nvert
          IF (seam(ibe)%expoint(iv)) THEN
            ix=seam(ibe)%vertex(iv)%intxy(1)
            iy=seam(ibe)%vertex(iv)%intxy(2)
            prod(ibe)%arr(1,ix,iy)=oper(ibe)%arr(1,ix,iy)-
     $                             sln(ibe)%arr(1,ix,iy)
            IF (seam(ibe)%expoint(ivp).AND.ibe<=nrbl.AND.
     $          poly_degree>1) THEN
              ix=seam(ibe)%segment(iv)%intxys(1)
              iy=seam(ibe)%segment(iv)%intxys(2)
              IF (seam(ibe)%segment(iv)%h_side) THEN
                prod(ibe)%arrh(1,:,ix,iy)=oper(ibe)%arrh(1,:,ix,iy)-
     $                                    sln(ibe)%arrh(1,:,ix,iy)
              ELSE
                prod(ibe)%arrv(1,:,ix,iy)=oper(ibe)%arrv(1,:,ix,iy)-
     $                                    sln(ibe)%arrv(1,:,ix,iy)
              ENDIF
            ENDIF
          ENDIF
          ivp=iv
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE eqdot_delstlj
