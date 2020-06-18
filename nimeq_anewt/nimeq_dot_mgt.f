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
c     2. eqdot_delst.
c-----------------------------------------------------------------------
c     subprogram 1. eqdot_delstlj
c     this dot-product management routine computes the result of
c     applying the del-star operator on a scalar (lambda=psi/R**2) and
c     applying the mass matrix on a second scalar (mu0*J_phi/R), as
c     needed for the free-boundary Grad-Shafranov solves.
c
c     regularity conditions are not needed, because lambda and J_phi/R
c     have the same limiting behavior as n=0 scalars, which are not
c     modified by our regularity routines.
c-----------------------------------------------------------------------
      SUBROUTINE eqdot_delstlj(oper,prod,bc_oper)
      USE local
      USE input
      USE nimeq_ints
      USE boundary
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
      USE iter_utils
      USE nimeq_mod
      IMPLICIT NONE

      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,ncoil_save,iv,ivp,ix,iy,is,ierror
      REAL(r8) :: sflmin,sflmax,afldl,tmp
c-----------------------------------------------------------------------
c     the integrand computation interpolates be_eq(2:3) as the operand,
c     so transfer the passed data to be_eq, but set the ja quad-point
c     storage for the surface-lambda computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_assignq_vec(be_eq(ibl),oper(ibl),nqty=2_i4,
     $                          nstart1=2_i4,nstart2=1_i4)
        CALL rblock_qp_update(rb(ibl)%be_eq,rb(ibl)%qja_eq,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_assignq_vec(be_eq(ibl),oper(ibl),nqty=2_i4,
     $                          nstart1=2_i4,nstart2=1_i4)
        CALL tblock_qp_update(tb(ibl)%be_eq,tb(ibl)%qja_eq,tb(ibl))
      ENDDO
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
c     when Newton iteration is active and a target Ip is used, complete
c     the matrix-partitioning operations.  the off-diagonal elements
c     are stored in diperr.  first find A_f_lam.dlam, then divide by
c     -A_f_f, then multiply A_lam_f and add the result from the rest of
c     the product vector.
c-----------------------------------------------------------------------
      IF (ip_equil/=0._r8.AND.relerr<1._r8.AND.ip_weight>0._r8) THEN
        afldl=0._r8
        DO ibl=1,nbl
          CALL vector_assignq_vec(sln(ibl),oper(ibl),1_i4,
     $                            nstart1=1_i4,nstart2=1_i4)
          CALL vector_assignq_vec(rwork3(ibl),diperr(ibl),1_i4,
     $                            nstart1=1_i4,nstart2=2_i4)
          CALL iter_dot(afldl,rwork3(ibl),sln(ibl),seam(ibl),ibl,nrbl,
     $                  poly_degree,.true.)
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(afldl,tmp,1,mpi_nim_real,mpi_sum,
     $         comm_layer,ierror)
          afldl=tmp
        ENDIF
        afldl=-afldl/(dipdf*ip_weight)
        DO ibl=1,nbl
          prod(ibl)%arr(1,:,:)=prod(ibl)%arr(1,:,:)+
     $                         afldl*diperr(ibl)%arr(1,:,:)
          IF (ASSOCIATED(prod(ibl)%arrh))
     $      prod(ibl)%arrh(1,:,:,:)=prod(ibl)%arrh(1,:,:,:)+
     $                              afldl*diperr(ibl)%arrh(1,:,:,:)
          IF (ASSOCIATED(prod(ibl)%arrv))
     $      prod(ibl)%arrv(1,:,:,:)=prod(ibl)%arrv(1,:,:,:)+
     $                              afldl*diperr(ibl)%arrv(1,:,:,:)
          IF (ASSOCIATED(prod(ibl)%arri))
     $      prod(ibl)%arri(1,:,:,:)=prod(ibl)%arri(1,:,:,:)+
     $                              afldl*diperr(ibl)%arri(1,:,:,:)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     the following evaluates the surface flux from the change in
c     internal currents represented in the direction vector,
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
c     replace surface-lambda values with the result based on internal
c     currents.  each surface-lamda row effectively amounts to
c     lambda_i-SUM(mat_ij*I_j,j-interior)=SUM(mat_ij*I_j,j-exterior),
c     and the entire relation is multiplied by a weight factor.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        ivp=seam(ibe)%nvert
        DO iv=1,seam(ibe)%nvert
          IF (seam(ibe)%expoint(iv)) THEN
            ix=seam(ibe)%vertex(iv)%intxy(1)
            iy=seam(ibe)%vertex(iv)%intxy(2)
            prod(ibe)%arr(1,ix,iy)=
     $        bweight*(oper(ibe)%arr(1,ix,iy)-sln(ibe)%arr(1,ix,iy))
            IF (seam(ibe)%expoint(ivp).AND.ibe<=nrbl.AND.
     $          poly_degree>1) THEN
              ix=seam(ibe)%segment(iv)%intxys(1)
              iy=seam(ibe)%segment(iv)%intxys(2)
              IF (seam(ibe)%segment(iv)%h_side) THEN
                prod(ibe)%arrh(1,:,ix,iy)=
     $            bweight*(oper(ibe)%arrh(1,:,ix,iy)-
     $                       sln(ibe)%arrh(1,:,ix,iy))
              ELSE
                prod(ibe)%arrv(1,:,ix,iy)=
     $            bweight*(oper(ibe)%arrv(1,:,ix,iy)-
     $                      sln(ibe)%arrv(1,:,ix,iy))
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
c-----------------------------------------------------------------------
c     subprogram 2. eqdot_delst
c     this dot-product management routine computes the result of
c     applying the del-star operator on a scalar (lambda=psi/R**2) for
c     fixed-boundary Grad-Shafranov solves.
c
c     regularity conditions are not needed, because lambda has
c     the same limiting behavior as n=0 scalars, which are not
c     modified by our regularity routines.
c-----------------------------------------------------------------------
      SUBROUTINE eqdot_delst(oper,prod,bc_oper)
      USE local
      USE input
      USE nimeq_ints
      USE boundary
      USE global
      USE fields
      USE computation_pointers
      USE rblock
      USE tblock
      USE seam_storage_mod
      USE edge
      USE mpi_nim
      USE pardata
      USE iter_utils
      USE nimeq_mod
      IMPLICIT NONE

      TYPE(vector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
      LOGICAL, INTENT(IN) :: bc_oper

      INTEGER(i4) :: ibl,ibe,ix,iy,iv,ivp,ierror
      REAL(r8) :: afldl,tmp
c-----------------------------------------------------------------------
c     use the rwork3 data for work space with block interpolation.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_assign_vec(rwork3(ibl),oper(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     use copies to the complex crhs data structure to use the standard
c     boundary routines.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL cvector_assign_vec(crhs(ibe),rwork3(ibe),'real',1_i4)
        CALL dirichlet_rhs(crhs(ibe),seam(ibe),'all',1_i4)
        CALL vector_assign_cvec(rwork3(ibe),crhs(ibe),'real',1_i4)
      ENDDO
c-----------------------------------------------------------------------
c     evaluate the matrix vector product directly through the
c     finite-element computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL rblock_get_rhs(rb(ibl),prod(ibl),delstar_dot,1_i4)
        CALL edge_load_arr(prod(ibl),1_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_get_rhs(tb(ibl),prod(ibl),delstar_dot,1_i4)
        CALL edge_load_arr(prod(ibl),1_i4,0_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     sum contributions over block boundaries.
c-----------------------------------------------------------------------
      CALL edge_network(1_i4,0_i4,poly_degree-1_i4,.false.)
      DO ibl=1,nrbl
        CALL edge_unload_arr(prod(ibl),1_i4,poly_degree-1_i4,seam(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL edge_unload_arr(prod(ibl),1_i4,0_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     when Newton iteration is active and a target Ip is used, complete
c     the matrix-partitioning operations.  the off-diagonal elements
c     are stored in diperr.  first find A_f_lam.dlam, then divide by
c     -A_f_f, then multiply A_lam_f and add the result from the rest of
c     the product vector.
c-----------------------------------------------------------------------
      IF (ip_equil/=0._r8.AND.relerr<1._r8.AND.ip_weight>0._r8) THEN
        afldl=0._r8
        DO ibl=1,nbl
          CALL vector_assignq_vec(rwork4(ibl),diperr(ibl),1_i4,
     $                            nstart1=1_i4,nstart2=2_i4)
          CALL iter_dot(afldl,rwork3(ibl),rwork4(ibl),seam(ibl),ibl,
     $                  nrbl,poly_degree,.true.)
        ENDDO
        IF (nprocs>1) THEN
          CALL mpi_allreduce(afldl,tmp,1,mpi_nim_real,mpi_sum,
     $         comm_layer,ierror)
          afldl=tmp
        ENDIF
        afldl=-afldl/(dipdf*ip_weight)
        DO ibl=1,nbl
c         CALL vector_assignq_vec(rwork3(ibl),diperr(ibl),1_i4,
c    $                            nstart1=1_i4,nstart2=1_i4)
c         CALL vector_add(prod(ibl),rwork3(ibl),v2fac=afldl)
          prod(ibl)%arr(1,:,:)=prod(ibl)%arr(1,:,:)+
     $                         afldl*diperr(ibl)%arr(1,:,:)
          IF (ASSOCIATED(prod(ibl)%arrh))
     $      prod(ibl)%arrh(1,:,:,:)=prod(ibl)%arrh(1,:,:,:)+
     $                              afldl*diperr(ibl)%arrh(1,:,:,:)
          IF (ASSOCIATED(prod(ibl)%arrv))
     $      prod(ibl)%arrv(1,:,:,:)=prod(ibl)%arrv(1,:,:,:)+
     $                              afldl*diperr(ibl)%arrv(1,:,:,:)
          IF (ASSOCIATED(prod(ibl)%arri))
     $      prod(ibl)%arri(1,:,:,:)=prod(ibl)%arri(1,:,:,:)+
     $                              afldl*diperr(ibl)%arri(1,:,:,:)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     replace surface product values with the respective operand
c     values to reproduce the effect of having 1 in the diagonal element
c     of the matrix at surface nodes with the rest of the row being 0.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        ivp=seam(ibe)%nvert
        DO iv=1,seam(ibe)%nvert
          IF (seam(ibe)%expoint(iv)) THEN
            ix=seam(ibe)%vertex(iv)%intxy(1)
            iy=seam(ibe)%vertex(iv)%intxy(2)
            prod(ibe)%arr(1,ix,iy)=oper(ibe)%arr(1,ix,iy)
            IF (seam(ibe)%expoint(ivp).AND.ibe<=nrbl.AND.
     $          poly_degree>1) THEN
              ix=seam(ibe)%segment(iv)%intxys(1)
              iy=seam(ibe)%segment(iv)%intxys(2)
              IF (seam(ibe)%segment(iv)%h_side) THEN
                prod(ibe)%arrh(1,:,ix,iy)=oper(ibe)%arrh(1,:,ix,iy)
              ELSE
                prod(ibe)%arrv(1,:,ix,iy)=oper(ibe)%arrv(1,:,ix,iy)
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
      END SUBROUTINE eqdot_delst
