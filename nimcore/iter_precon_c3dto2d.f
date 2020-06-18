c-----------------------------------------------------------------------
c     file iter_precon_c3dto2d.f
c     this file has the standard preconditioner for 3D solves, which
c     does static condensation and calls 2D preconditioning routines
c     for each Fourier component.  the module just holds data structures
c     that need to be retained during the course of one 3D solve.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0.  iter_pre_c3dto2d_mod
c     1.  iter_pre_c3dto2d.
c     2.  iter_pre_c3dto2d_init.
c     3.  iter_pre_c3dto2d_dealloc.
c-----------------------------------------------------------------------
c     module for the 2D data structures that can be allocated at the
c     start of a 3D solve and deallocated after completion.  naming
c     convention is r=rhs, z=preconditioned residual, t=temporary, and
c     p=pack.
c-----------------------------------------------------------------------
      MODULE iter_pre_c3dto2d_mod
      USE local
      USE vector_type_mod
      IMPLICIT NONE

      TYPE(cvector_2D_type), DIMENSION(:), POINTER, SAVE ::
     $                       c2r,c2z,c2t,c2p
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_pre_c3dto2d_mod

c-----------------------------------------------------------------------
c     subprogram 1. iter_pre_c3dto2d. 
c     this routine collect the steps for calling the 2D preconditioner
c     for each Fourier components.  each component also goes through
c     static condensation to make the 2D preconditioners efficient.
c-----------------------------------------------------------------------
      SUBROUTINE iter_pre_c3dto2d(mat_str,fac_str,ctmp,poly_deg,nqty,
     $                            nbdsc,nqdsc,nrb,ntotb,nfour,precon)

      USE matrix_mod
      USE factor_type_mod
      USE seam_storage_mod
      USE edge
      USE iter_pre_c3dto2d_mod
      USE time

      TYPE(complex_matrix_type), DIMENSION(nfour), INTENT(IN) :: mat_str
      TYPE(complex_factor_type), DIMENSION(nfour), INTENT(IN) :: fac_str
      TYPE(cvector_type), DIMENSION(ntotb), INTENT(INOUT) :: ctmp
      INTEGER(i4), INTENT(IN) :: poly_deg,nfour,nqty,nrb,ntotb
      INTEGER(i4), INTENT(IN) :: nbdsc,nqdsc  !  for discontin. bases
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: ibl,ix,iy,iv,im
c-----------------------------------------------------------------------
c     interface block for 2D preconditioner driver.
c-----------------------------------------------------------------------
      INCLUDE "iter_precon_intf.inc"
c-----------------------------------------------------------------------
c     return if there is no preconditioner.
c-----------------------------------------------------------------------
      IF (precon=='no prec') RETURN
c-----------------------------------------------------------------------
c     apply the preconditioning operation over the poloidal plane
c     for each Fourier component.
c-----------------------------------------------------------------------
        mode_loop: DO im=1,nfour
c-----------------------------------------------------------------------
c         transfer the ctmp to a cvector_2D_type, and eliminate 
c         interior elements.  first multiply block-border coefficients
c         by their weights to mimic finite_elements operations.
c-----------------------------------------------------------------------
          IF (poly_deg>1) THEN
            DO ibl=1,ntotb
              CALL cvector_2D_assign_cvec(c2z(ibl),ctmp(ibl),im)
              DO iv=1,seam(ibl)%nvert
                ix=seam(ibl)%vertex(iv)%intxy(1)
                iy=seam(ibl)%vertex(iv)%intxy(2)
                c2z(ibl)%arr(:,ix,iy)=c2z(ibl)%arr(:,ix,iy)
     $                                 *seam(ibl)%vertex(iv)%ave_factor
                ix=seam(ibl)%segment(iv)%intxys(1)
                iy=seam(ibl)%segment(iv)%intxys(2)
                IF (seam(ibl)%segment(iv)%h_side) THEN
                  c2z(ibl)%arrh(:,:,ix,iy)=c2z(ibl)%arrh(:,:,ix,iy)
     $                                 *seam(ibl)%segment(iv)%ave_factor
                ELSE
                  c2z(ibl)%arrv(:,:,ix,iy)=c2z(ibl)%arrv(:,:,ix,iy)
     $                                 *seam(ibl)%segment(iv)%ave_factor
                ENDIF
              ENDDO
c-----------------------------------------------------------------------
c             if there are discontinuous fields, pack their coefficients
c             together with the interior coefficients of the continuous
c             fields.
c-----------------------------------------------------------------------
              IF (nqdsc>0.AND.ibl<=nrb) THEN
                CALL cvector_2D_pack_cvec2(c2z(ibl),c2p(ibl),nqty,
     $                                     (poly_deg-1_i4)**2,nqdsc,
     $                                     nbdsc)
                c2p(ibl)%arrh=>c2z(ibl)%arri ! holds original arri mem
                c2p(ibl)%arrv=>c2r(ibl)%arri
                c2z(ibl)%arri=>c2p(ibl)%arri ! flat arrays for packing
                c2r(ibl)%arri=>c2p(ibl)%arrtmp
              ENDIF
            ENDDO
c-----------------------------------------------------------------------
c           eliminate interior (and discontinuous) coefficients.
c-----------------------------------------------------------------------
            CALL timer(timestart)
            CALL matelim_presolve(mat_str(im),c2z,c2r,nqty)
            CALL timer(timeend)
            time_stcon=time_stcon+timeend-timestart
c-----------------------------------------------------------------------
c           combine contributions across block borders.
c-----------------------------------------------------------------------
            DO ibl=1,ntotb
c-----------------------------------------------------------------------
c             unpack the discontinuous fields into the c2r%arrtmp
c             storage.  the first copy overwrites the old packed data.
c-----------------------------------------------------------------------
              IF (nqdsc>0.AND.ibl<=nrb) THEN
                c2p(ibl)%arri(:,:,:,:)=c2p(ibl)%arrtmp(:,:,:,:)
                c2r(ibl)%arri=>c2p(ibl)%arrv   ! original arri assigmnts
                c2z(ibl)%arri=>c2p(ibl)%arrh
                CALL cvector_2D_unpack_cvec2(c2p(ibl),c2r(ibl),nqty,
     $                                       (poly_deg-1_i4)**2,nqdsc,
     $                                       nbdsc,.true.)
                NULLIFY(c2p(ibl)%arrh,c2p(ibl)%arrv)
              ENDIF
              CALL edge_load_2D_carr(c2r(ibl),nqty,poly_deg-1_i4,
     $                               seam(ibl))
            ENDDO
            CALL edge_network(nqty,1_i4,poly_deg-1_i4,.false.)
            DO ibl=1,ntotb
              CALL edge_unload_2D_carr(c2r(ibl),nqty,
     $                                 poly_deg-1_i4,seam(ibl))
            ENDDO
          ELSE
            DO ibl=1,ntotb
              CALL cvector_2D_assign_cvec(c2r(ibl),ctmp(ibl),im)
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         call 2D preconditioner.
c-----------------------------------------------------------------------
          CALL iter_pre_comp(mat_str(im),fac_str(im),c2r,c2z,
     $                       c2t,nqty,nrb,ntotb,precon)
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.
c-----------------------------------------------------------------------
          IF (poly_deg>1) THEN
            DO ibl=1,ntotb
c-----------------------------------------------------------------------
c             if there are discontinuous fields, information from their
c             elimination is in c2r%arrtmp.  pack this together with
c             information from eliminating the interior coefficients
c             (in c2r%arri) before the postsolve step.
c-----------------------------------------------------------------------
              IF (nqdsc>0.AND.ibl<=nrb) THEN
                CALL cvector_2D_pack_cvec2(c2r(ibl),c2p(ibl),nqty,
     $                                     (poly_deg-1_i4)**2,nqdsc,
     $                                     nbdsc)
                c2p(ibl)%arrh=>c2z(ibl)%arri ! hold original arri mem
                c2p(ibl)%arrv=>c2r(ibl)%arri
                c2r(ibl)%arri=>c2p(ibl)%arri ! flat packed array for c2r
                c2z(ibl)%arri=>c2p(ibl)%arrtmp ! flat array for c2z
              ENDIF
            ENDDO
c-----------------------------------------------------------------------
c           reassemble the coefficients that had been eliminated.
c-----------------------------------------------------------------------
            CALL timer(timestart)
            CALL matelim_postsolve(mat_str(im),c2z,c2r,nqty)
            CALL timer(timeend)
            time_stcon=time_stcon+timeend-timestart
            DO ibl=1,ntotb
c-----------------------------------------------------------------------
c             reset pointers to their original memory assignments, then
c             unpack the discontinuous fields into the c2z%arrtmp
c             storage.
c-----------------------------------------------------------------------
              IF (nqdsc>0.AND.ibl<=nrb) THEN
                c2z(ibl)%arri=>c2p(ibl)%arrh   ! original arri assigmnts
                c2r(ibl)%arri=>c2p(ibl)%arrv
                CALL cvector_2D_unpack_cvec2(c2p(ibl),c2z(ibl),nqty,
     $                                       (poly_deg-1_i4)**2,nqdsc,
     $                                       nbdsc,.true.)
                NULLIFY(c2p(ibl)%arrh,c2p(ibl)%arrv)
              ELSE
                c2z(ibl)%arri(:,:,:,:)=c2r(ibl)%arri(:,:,:,:)
              ENDIF
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         make sure that the n=0 part is real before moving to the next
c         component.
c-----------------------------------------------------------------------
          IF (mat_str(im)%fcomp==0) THEN
            DO ibl=1,ntotb
              c2z(ibl)%arr=REAL(c2z(ibl)%arr,r8)
              IF (ASSOCIATED(c2z(ibl)%arrh)) THEN
                c2z(ibl)%arrh=REAL(c2z(ibl)%arrh,r8)
                c2z(ibl)%arrv=REAL(c2z(ibl)%arrv,r8)
                c2z(ibl)%arri=REAL(c2z(ibl)%arri,r8)
              ENDIF
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         combine the result of this preconditioning step as needed for
c         Gauss-Seidel or Jacobi iteration.
c-----------------------------------------------------------------------
          DO ibl=1,ntotb
            CALL cvector_assign_cvec2(ctmp(ibl),c2z(ibl),im)
          ENDDO
        ENDDO mode_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_pre_c3dto2d
c-----------------------------------------------------------------------
c     subprogram 2. iter_pre_c3dto2d_init. 
c     allocate space for the linear algebra.
c-----------------------------------------------------------------------
      SUBROUTINE iter_pre_c3dto2d_init(nqty,poly_deg,nbdsc,nqdsc,nrb,
     $                                 nbl,anycv)
      USE iter_pre_c3dto2d_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nqty,poly_deg,nbdsc,nqdsc,nrb,nbl
      TYPE(cvector_type), DIMENSION(nbl), INTENT(IN) :: anycv

      INTEGER(i4) :: ibl,mxb,myb,ncomb
c-----------------------------------------------------------------------
c     find array dimensions.
c-----------------------------------------------------------------------
      ncomb=nqty*(poly_deg-1)**2+nqdsc*nbdsc
c-----------------------------------------------------------------------
c     allocate vector structures.
c-----------------------------------------------------------------------
      ALLOCATE(c2z(nbl))
      ALLOCATE(c2r(nbl))
      ALLOCATE(c2t(nbl))
      ALLOCATE(c2p(nbl))
      DO ibl=1,nrb
        mxb=SIZE(anycv(ibl)%arr,2)-1
        myb=SIZE(anycv(ibl)%arr,3)-1
        IF (nqdsc>0) THEN
          CALL vector_type_alloc(c2z(ibl),poly_deg,mxb,myb,nqty,
     $                           nbdsc,nqdsc)
          CALL vector_type_alloc(c2r(ibl),poly_deg,mxb,myb,nqty,
     $                           nbdsc,nqdsc)
          CALL vector_type_alloc(c2t(ibl),poly_deg,mxb,myb,nqty,
     $                           nbdsc,nqdsc)
          ALLOCATE(c2p(ibl)%arri(ncomb,1,mxb,myb))
          ALLOCATE(c2p(ibl)%arrtmp(ncomb,1,mxb,myb))
        ELSE
          CALL vector_type_alloc(c2z(ibl),poly_deg,mxb,myb,nqty)
          CALL vector_type_alloc(c2r(ibl),poly_deg,mxb,myb,nqty)
          CALL vector_type_alloc(c2t(ibl),poly_deg,mxb,myb,nqty)
          NULLIFY(c2p(ibl)%arri,c2p(ibl)%arrtmp)
        ENDIF
        NULLIFY(c2p(ibl)%arr,c2p(ibl)%arrh,c2p(ibl)%arrv)
      ENDDO
      DO ibl=nrb+1,nbl
        mxb=SIZE(anycv(ibl)%arr,2)-1
        CALL vector_type_alloc(c2z(ibl),1_i4,mxb,0_i4,nqty)
        CALL vector_type_alloc(c2r(ibl),1_i4,mxb,0_i4,nqty)
        CALL vector_type_alloc(c2t(ibl),1_i4,mxb,0_i4,nqty)
        NULLIFY(c2p(ibl)%arr,c2p(ibl)%arrh,c2p(ibl)%arrv)
        NULLIFY(c2p(ibl)%arri,c2p(ibl)%arrtmp)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_pre_c3dto2d_init
c-----------------------------------------------------------------------
c     subprogram 3. iter_pre_c3dto2d_dealloc. 
c     deallocate temporary space.
c-----------------------------------------------------------------------
      SUBROUTINE iter_pre_c3dto2d_dealloc
      USE iter_pre_c3dto2d_mod
      IMPLICIT NONE

      INTEGER(i4) :: ibl
c-----------------------------------------------------------------------
c     deallocate vector structures.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(c2z)
        CALL vector_type_dealloc(c2z(ibl))
        CALL vector_type_dealloc(c2r(ibl))
        CALL vector_type_dealloc(c2t(ibl))
        IF (ASSOCIATED(c2p(ibl)%arri)) THEN
          DEALLOCATE(c2p(ibl)%arri)
          DEALLOCATE(c2p(ibl)%arrtmp)
        ENDIF
      ENDDO
      DEALLOCATE(c2z,c2r,c2t,c2p)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_pre_c3dto2d_dealloc
