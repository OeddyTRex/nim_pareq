c-----------------------------------------------------------------------
c     file iter_precon_real.f
c     This set of modules holds nimrod's traditional preconditioning
c     routines for real 2D algebraic systems over the meshed plane.
c     these routines had been the bulk of iter_cg_f90.f.  The CG solver
c     module for real 2D systems, itself, now appears separately in
c     the iter_cg_real.f file.  See that file for documentation of
c     older changes to the preconditioners that are now here.
c     1/9/19, C. Sovinec
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_real_direct
c     1. fac_dir.
c     2. solve_dir.
c     module iter_real_line
c     3. iter_line_fac
c     4. iter_line_solve
c     5. iter_ld_fac
c     6. iter_ld_solve
c     7. iter_solve_sym.
c     8. iter_gl_ld_fac
c     9. iter_gl_ld_solve
c     10. iter_ld_elim
c     11. iter_inner_pr
c     module iter_real_fac
c     12. iter_fac_alloc_real.
c     13. iter_fac_dealloc_real.
c     14. iter_factor_real.
c     15. iter_mat_com.
c     16. iter_mat_rest.
c     17. iter_pre_real (external).
c-----------------------------------------------------------------------
c     module for the block-direct type preconditioner.
c-----------------------------------------------------------------------
      MODULE iter_real_direct
      USE local
      USE matrix_type_mod
      USE factor_type_mod
      USE pardata
      USE mpi_nim
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fac_dir.
c     form and factor matrices for a global direct solve.
c-----------------------------------------------------------------------
      SUBROUTINE fac_dir(mat,fac,spp,nqty,nrb,precon)
      USE io
      USE time

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      TYPE(sparsity_pattern), INTENT(INOUT) :: spp
      INTEGER(i4), INTENT(IN) :: nqty,nrb
      CHARACTER(*), INTENT(IN) :: precon

      TYPE :: it1d_rtype
        REAL(r8), DIMENSION(:), ALLOCATABLE :: data
      END TYPE it1d_rtype

      INTEGER(i4) :: iopt,ierr,ist,ien
      INTEGER(i4), PARAMETER :: nrhs=1
      REAL(r8), DIMENSION(:), ALLOCATABLE :: acc,atmp
      REAL(r8) :: bdum
      REAL(r8) :: timefac_st,timefac_en
      INTEGER(i4), PARAMETER :: mpi_nim_iter = mpi_nim_real

      CHARACTER(128) :: msg

      TYPE(it1d_rtype), DIMENSION(:), POINTER :: recvstr,sendstr
      INTEGER(i4), DIMENSION(mpi_status_size) :: status
      LOGICAL :: sdone

c-----------------------------------------------------------------------
c     reset the index arrays as needed.
c-----------------------------------------------------------------------
      IF (spp%acc_lustored) THEN !  subsequent factorization
        iopt=2
        IF (spp%sparsity_distributed) THEN
          spp%j_acc(:)=spp%j_acc_save(:)
          spp%start_loc(:)=spp%start_acc_save(:)
        ELSE
          spp%j_acc(:)=spp%j_acc_save(:)
          spp%start_acc(:)=spp%start_acc_save(:)
        ENDIF
      ELSE !  first factorization
        iopt=1
        spp%acc_lustored=.true.
      ENDIF
c-----------------------------------------------------------------------
c     create the compressed-row data array for the nonzero matrix
c     elements where this processor's grid-blocks contribute. then
c     load the matrix values into acc. these functions are defined in
c     iter_cg_shared.F
c-----------------------------------------------------------------------
      IF (spp%sparsity_distributed) THEN
        ALLOCATE(acc(spp%indst:spp%indend))
        CALL load_matrix_dist
      ELSE
        ist=MIN(spp%indst,MINVAL(spp%jentry))
        ien=MAX(spp%indend,MAXVAL(spp%jentry))
        ALLOCATE(acc(ist:ien))
        CALL load_matrix
      ENDIF
c-----------------------------------------------------------------------
c     Factor main block with the appropriate library.
c-----------------------------------------------------------------------
      CALL timer(timefac_st)
      IF (precon=='lapack') THEN
        IF (.NOT.single_pr)  THEN
          CALL dpbtrf('L',spp%nrow,spp%nbw,fac%a11(1,1),
     $                spp%nbw+1_i4,ierr)
        ELSE
          CALL spbtrf('L',spp%nrow,spp%nbw,fac%a11(1,1),
     $                spp%nbw+1_i4,ierr)
        ENDIF
        IF (ierr/=0) THEN
          WRITE(msg,'(a,i9)') 'Fac_dir real: **dpbtrf** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ELSE IF (precon=='seq_slu') THEN
        CALL c_fortran_dgssv(
     $         iopt,spp%nrow,spp%nnz,nrhs,acc,spp%j_acc,
     $         spp%start_acc,bdum,spp%nrow,spp%acc_handle,ierr)
        IF (ierr==-999999) THEN
          CALL nim_stop
     $      ('Fac_dir real: the sequential SLU library is not linked.')
        ELSE IF (ierr/=0) THEN
          WRITE(msg,'(a,i9)') 'Fac_dir real: **dgssv** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ELSE IF (precon(1:5)=='slu_d') THEN
        IF (precon(1:7)=='slu_dst') THEN
          CALL c_fortran_pdloc(
     $           iopt,spp%nrow,spp%mloc,spp%nnzloc,spp%fstrow,nrhs,
     $           acc(spp%indst),spp%j_acc(spp%indst),spp%start_loc(0),
     $           bdum,INT(spp%mloc,i4),spp%acc_handle,slugrid_handle,
     $           spp%iopts,spp%dopts,ierr)
        ELSE
          CALL c_fortran_pdglob(
     $           iopt,spp%nrow,spp%nnz,nrhs,acc,spp%j_acc,
     $           spp%start_acc,bdum,INT(spp%nrow,i4),spp%acc_handle,
     $           slugrid_handle,ierr)
        ENDIF
        IF (ierr==-999999) THEN
          CALL nim_stop
     $      ('Fac_dir real: the SLU_DIST library is not linked.')
        ELSE IF (ierr/=0) THEN
          WRITE(msg,'(a,i9)') 'Fac_dir real: **pdloc fac** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ENDIF
      CALL timer(timefac_en)
      time_exfact=time_exfact+timefac_en-timefac_st
      DEALLOCATE(acc)
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1.1. load_matrix.
c     load the matrix with a global sparsity pattern.
c-----------------------------------------------------------------------
      SUBROUTINE load_matrix

      INTEGER(i4) :: ibl,nv,mx,my,n,iq,jq,iv,jv,lx,ly,ix,iy,it,jt,
     $               ix0,iy0,jx0,jy0,iind,jind,jx,jy,inq,jnq,ierr,
     $               isnd,ircv,ns,nbmax
      INTEGER(i4), DIMENSION(:,:,:), POINTER :: iarr,jarr
      INTEGER(i4), DIMENSION(:), POINTER :: jcol,irst

      INTEGER(i4), SAVE :: ilabel=0
      INTEGER(i4), PARAMETER :: iwrcc=0
      CHARACTER(24) :: cdat,crow,ccpt
      CHARACTER(9) :: cnnz

c-----------------------------------------------------------------------
c     post receives for contributions coming from other processors to
c     rows owned by this processor.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        IF (spp%nrcv>0) ALLOCATE(recvstr(spp%nrcv))
        IF (spp%nsnd>0) ALLOCATE(sendstr(spp%nsnd))
        DO ircv=1,spp%nrcv
          ALLOCATE(recvstr(ircv)%data(spp%recvtot(ircv)))
          CALL mpi_irecv(recvstr(ircv)%data(1),spp%recvtot(ircv),
     $                   mpi_nim_iter,spp%recvlist(ircv),0_i4,
     $                   comm_layer,spp%recv_req(ircv),ierr)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     load the nonzero matrix elements into the SuperLU compressed row/
c     column storage.  full-storage SLU_DIST requires CC.
c-----------------------------------------------------------------------
      acc=0
      jcol=>spp%j_acc
      irst=>spp%start_acc
      DO ibl=1,nrb
        mx=fac%mat_info(ibl)%mx
        my=fac%mat_info(ibl)%my
        nbmax=SIZE(fac%bl_spp(ibl)%row_ind)
        DO it=1,nbmax
          ix0=mat%rbl_mat(ibl)%ix0(it)
          iy0=mat%rbl_mat(ibl)%iy0(it)
          inq=mat%rbl_mat(ibl)%nq_type(it)
          iarr=>fac%bl_spp(ibl)%row_ind(it)%rarr
          DO iy=iy0,my
            DO ix=ix0,mx
              DO iq=1,inq
                iind=iarr(iq,ix,iy)
                DO jt=1,nbmax
                  jx0=mat%rbl_mat(ibl)%ix0(jt)
                  jy0=mat%rbl_mat(ibl)%iy0(jt)
                  jnq=mat%rbl_mat(ibl)%nq_type(jt)
                  jarr=>fac%bl_spp(ibl)%row_ind(jt)%rarr
                  DO ly=jy0-1,1-iy0
                    jy=iy+ly
                    IF (jy<jy0.OR.jy>my) CYCLE
                    DO lx=jx0-1,1-ix0
                      jx=ix+lx
                      IF (jx<jx0.OR.jx>mx) CYCLE
                      IF (precon=='slu_dist' .or.
     $                    precon=='seq_slu') THEN  !  compressed col.
                        DO jq=1,jnq
                          jind=jarr(jq,jx,jy)
                          DO n=irst(iind),irst(iind+1)-1
                            IF (jind==jcol(n)) THEN
                              acc(n)=acc(n)+
     $                               mat%rbl_mat(ibl)%mat(it,jt)%
     $                                   arr(iq,-lx,-ly,jq,jx,jy)
                              EXIT
                            ENDIF
                          ENDDO
                        ENDDO
                      ELSE !  compressed row
                        DO jq=1,jnq
                          jind=jarr(jq,jx,jy)
                          DO n=irst(iind),irst(iind+1)-1
                            IF (jind==jcol(n)) THEN
                              acc(n)=acc(n)+
     $                               mat%rbl_mat(ibl)%mat(jt,it)%
     $                                   arr(jq,lx,ly,iq,ix,iy)
                              EXIT
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     communicate if there is more than one processor per layer.
c     first, send the contributions to rows that belong to other
c     processors.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        DO isnd=1,spp%nsnd
          ALLOCATE(sendstr(isnd)%data(spp%sendtot(isnd)))
          iv=spp%sendst(isnd)
          DO jv=1,spp%sendtot(isnd)
            sendstr(isnd)%data(jv)=acc(spp%jentry(iv))
            iv=iv+1
          ENDDO
          CALL mpi_isend(sendstr(isnd)%data(1),spp%sendtot(isnd),
     $                   mpi_nim_iter,spp%sendlist(isnd),0_i4,
     $                   comm_layer,spp%send_req(isnd),ierr)
        ENDDO
c-----------------------------------------------------------------------
c       now, collect information from other processors for the rows
c       owned by this processor, as it becomes available.  also check
c       on sends.
c-----------------------------------------------------------------------
        ns=0
        DO iv=1,spp%nrcv
          CALL mpi_waitany(spp%nrcv,spp%recv_req,ircv,status,ierr)
          DO jv=1,spp%recvtot(ircv)
            acc(spp%recvjentry(jv,ircv))=acc(spp%recvjentry(jv,ircv))+
     $                                   recvstr(ircv)%data(jv)
          ENDDO
          DEALLOCATE(recvstr(ircv)%data)

          IF (spp%nsnd>0) THEN
            CALL mpi_testany(spp%nsnd,spp%send_req,isnd,sdone,status,
     $                       ierr)
            IF (sdone.AND.ns<spp%nsnd) THEN
              ns=ns+1
              DEALLOCATE(sendstr(isnd)%data)
            ENDIF
          ENDIF
        ENDDO

        DO iv=ns+1,spp%nsnd
          CALL mpi_waitany(spp%nsnd,spp%send_req,isnd,status,ierr)
          DEALLOCATE(sendstr(isnd)%data)
        ENDDO
        IF (spp%nrcv>0) DEALLOCATE(recvstr)
        IF (spp%nsnd>0) DEALLOCATE(sendstr)
      ENDIF
c-----------------------------------------------------------------------
c     full-storage communication uses the following to acquire rows
c     owned by other processors.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1 .AND. .NOT.spp%matrix_distributed) THEN
        ALLOCATE(atmp(0:spp%nnz-1))
        CALL mpi_allgatherv(acc(spp%algdispl(node_layer)),
     $       spp%algcount(node_layer),mpi_nim_iter,atmp(0),
     $       spp%algcount(0),spp%algdispl(0),mpi_nim_iter,
     $       comm_layer,ierr)
        acc=atmp
        DEALLOCATE(atmp)
      ENDIF
c-----------------------------------------------------------------------
c     The compressed row/column storage is written to text files when
c     the iwrcc parameter matches the numbered matrix label.
c-----------------------------------------------------------------------
      ilabel=ilabel+1
      IF (ilabel==iwrcc) THEN
        WRITE(cnnz,'(i9)') spp%nnz
        WRITE(cdat,'(a,2i1,a)') 'cc_data',ilabel,node,'_'//ADJUSTL(cnnz)
        WRITE(crow,'(a,2i1,a)') 'cc_rind',ilabel,node,'_'//ADJUSTL(cnnz)
        WRITE(ccpt,'(a,2i1,a)') 'cc_cptr',ilabel,node,'_'//ADJUSTL(cnnz)
        OPEN(UNIT=xy_unit,FILE=TRIM(cdat),STATUS='UNKNOWN')
        OPEN(UNIT=xt_unit,FILE=TRIM(crow),STATUS='UNKNOWN')
        OPEN(UNIT=yt_unit,FILE=TRIM(ccpt),STATUS='UNKNOWN')

        DO lx=spp%indst,spp%indst+spp%nnzloc-1
          WRITE(xy_unit,*) acc(lx)
          WRITE(xt_unit,*) spp%j_acc(lx)
        ENDDO
        DO lx=0,spp%nrow
          WRITE(yt_unit,*) spp%start_acc(lx)
        ENDDO
        CLOSE(xy_unit)
        CLOSE(xt_unit)
        CLOSE(yt_unit)
      ENDIF
c-----------------------------------------------------------------------
c     Load the banded symmetric storage array if needed.  Note that acc
c     row and column indices start from 0, whereas a11 starts from 1.
c-----------------------------------------------------------------------
      IF (precon=='lapack') THEN
        fac%a11=0
        iq=0
        DO lx=0,spp%nnz-1
          jq=spp%j_acc(lx)
          IF (lx==spp%start_acc(iq+1)) iq=iq+1
          IF (jq<iq) CYCLE
          fac%a11(jq-iq+1,iq+1)=acc(lx)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE load_matrix
c-----------------------------------------------------------------------
c     subprogram 1.2. load_matrix_dist.
c     load the matrix with a distributed sparsity pattern.
c-----------------------------------------------------------------------
      SUBROUTINE load_matrix_dist

      INTEGER(i4) :: ibl,nv,mx,my,n,iq,jq,iv,jv,lx,ly,ix,iy,it,jt,
     $               ix0,iy0,jx0,jy0,iind,jind,jx,jy,inq,jnq,ierr,
     $               isnd,ircv,ns,sndst,isndind,nbmax
      INTEGER(i4), DIMENSION(:,:,:), POINTER :: iarr,jarr
      INTEGER(i4), DIMENSION(:), POINTER :: jcol,irst

      INTEGER(i4), SAVE :: ilabel=0
      INTEGER(i4), PARAMETER :: iwrcc=0
      CHARACTER(24) :: cdat,crow,ccpt
      CHARACTER(9) :: cnnz

c-----------------------------------------------------------------------
c     post receives for contributions coming from other processors to
c     rows owned by this processor.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        IF (spp%nrcv>0) ALLOCATE(recvstr(spp%nrcv))
        DO ircv=1,spp%nrcv
          ALLOCATE(recvstr(ircv)%data(spp%recvtot(ircv)))
          CALL mpi_irecv(recvstr(ircv)%data(1),spp%recvtot(ircv),
     $                   mpi_nim_iter,spp%recvlist(ircv),0_i4,
     $                   comm_layer,spp%recv_req(ircv),ierr)
        ENDDO
        IF (spp%nsnd>0) ALLOCATE(sendstr(spp%nsnd))
        DO isnd=1,spp%nsnd
          sndst=spp%sendst(isnd)
          ALLOCATE(sendstr(isnd)%data(sndst:sndst+spp%sendtot(isnd)))
          sendstr(isnd)%data=0
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     load the nonzero matrix elements into the SuperLU compressed row/
c     column storage.
c-----------------------------------------------------------------------
      acc=0
      jcol=>spp%j_acc
      irst=>spp%start_acc
      DO ibl=1,nrb
        mx=fac%mat_info(ibl)%mx
        my=fac%mat_info(ibl)%my
        nbmax=SIZE(fac%bl_spp(ibl)%row_ind)
        DO it=1,nbmax
          ix0=mat%rbl_mat(ibl)%ix0(it)
          iy0=mat%rbl_mat(ibl)%iy0(it)
          inq=mat%rbl_mat(ibl)%nq_type(it)
          iarr=>fac%bl_spp(ibl)%row_ind(it)%rarr
          DO iy=iy0,my
            DO ix=ix0,mx
              DO iq=1,inq
                iind=iarr(iq,ix,iy)
                IF (iind >= spp%fstrow) THEN ! this row is local
                  DO jt=1,nbmax
                    jx0=mat%rbl_mat(ibl)%ix0(jt)
                    jy0=mat%rbl_mat(ibl)%iy0(jt)
                    jnq=mat%rbl_mat(ibl)%nq_type(jt)
                    jarr=>fac%bl_spp(ibl)%row_ind(jt)%rarr
                    DO ly=jy0-1,1-iy0
                      jy=iy+ly
                      IF (jy<jy0.OR.jy>my) CYCLE
                      DO lx=jx0-1,1-ix0
                        jx=ix+lx
                        IF (jx<jx0.OR.jx>mx) CYCLE
                        DO jq=1,jnq !  compressed row
                          jind=jarr(jq,jx,jy)
                          DO n=irst(iind),irst(iind+1)-1
                            IF (jind==jcol(n)) THEN
                              acc(n)=acc(n)+
     $                               mat%rbl_mat(ibl)%mat(jt,it)%
     $                                   arr(jq,lx,ly,iq,ix,iy)
                              EXIT
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ELSE ! this row is sent
                  IF (spp%sendrowst(spp%nsnd)<=iind) THEN
                    isnd=spp%nsnd
                  ELSE
                    DO isnd=1,spp%nsnd-1
                      IF (spp%sendrowst(isnd)<=iind   .AND.
     $                    spp%sendrowst(isnd+1)>iind) EXIT
                    ENDDO
                  ENDIF
                  DO isndind=1,SIZE(spp%sendstacc(isnd)%data,DIM=1)
                    IF (spp%sendstacc(isnd)%data(isndind,1)==iind) EXIT
                  ENDDO
                  DO jt=1,nbmax
                    jx0=mat%rbl_mat(ibl)%ix0(jt)
                    jy0=mat%rbl_mat(ibl)%iy0(jt)
                    jnq=mat%rbl_mat(ibl)%nq_type(jt)
                    jarr=>fac%bl_spp(ibl)%row_ind(jt)%rarr
                    DO ly=jy0-1,1-iy0
                      jy=iy+ly
                      IF (jy<jy0.OR.jy>my) CYCLE
                      DO lx=jx0-1,1-ix0
                        jx=ix+lx
                        IF (jx<jx0.OR.jx>mx) CYCLE
                        DO jq=1,jnq !  compressed row
                          jind=jarr(jq,jx,jy)
                          DO n=spp%sendstacc(isnd)%data(isndind,2),
     $                         spp%sendstacc(isnd)%data(isndind+1,2)-1
                            IF (jind==spp%jentry(n)) THEN
                               sendstr(isnd)%data(n)=
     $                               sendstr(isnd)%data(n)+
     $                               mat%rbl_mat(ibl)%mat(jt,it)%
     $                                   arr(jq,lx,ly,iq,ix,iy)
                              EXIT
                            ENDIF
                          ENDDO
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     communicate if there is more than one processor per layer.
c     first, send the contributions to rows that belong to other
c     processors.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        DO isnd=1,spp%nsnd
          CALL mpi_isend(sendstr(isnd)%data(spp%sendst(isnd)),
     $                   spp%sendtot(isnd),mpi_nim_iter,
     $                   spp%sendlist(isnd),0_i4,comm_layer,
     $                   spp%send_req(isnd),ierr)
        ENDDO
c-----------------------------------------------------------------------
c       now, collect information from other processors for the rows
c       owned by this processor, as it becomes available.  also check
c       on sends.
c-----------------------------------------------------------------------
        ns=0
        DO iv=1,spp%nrcv
          CALL mpi_waitany(spp%nrcv,spp%recv_req,ircv,status,ierr)
          DO jv=1,spp%recvtot(ircv)
            acc(spp%recvjentry(jv,ircv))=acc(spp%recvjentry(jv,ircv))+
     $                                   recvstr(ircv)%data(jv)
          ENDDO
          DEALLOCATE(recvstr(ircv)%data)

          IF (spp%nsnd>0) THEN
            CALL mpi_testany(spp%nsnd,spp%send_req,isnd,sdone,status,
     $                       ierr)
            IF (sdone.AND.ns<spp%nsnd) THEN
              ns=ns+1
              DEALLOCATE(sendstr(isnd)%data)
            ENDIF
          ENDIF
        ENDDO

        DO iv=ns+1,spp%nsnd
          CALL mpi_waitany(spp%nsnd,spp%send_req,isnd,status,ierr)
          DEALLOCATE(sendstr(isnd)%data)
        ENDDO
        IF (spp%nrcv>0) DEALLOCATE(recvstr)
        IF (spp%nsnd>0) DEALLOCATE(sendstr)
      ENDIF

      RETURN
      END SUBROUTINE load_matrix_dist

      END SUBROUTINE fac_dir
c-----------------------------------------------------------------------
c     subprogram 2. solve_dir.
c     perform the direct solve.
c-----------------------------------------------------------------------
      SUBROUTINE solve_dir(fac,spp,zee,res,nqty,pd,nbl,precon)
      USE vector_type_mod
      USE edge
      USE seam_storage_mod
      USE time

      INTEGER(i4), INTENT(IN):: nqty,pd,nbl
      CHARACTER(*), INTENT(IN) :: precon
      TYPE(matrix_factor_type), INTENT(IN) :: fac
      TYPE(sparsity_pattern), INTENT(IN) :: spp
      TYPE(vector_type), DIMENSION(nbl), INTENT(INOUT) :: res,zee

      INTEGER(i4) :: ix,iy,mx,my,ib,iq,jq,pdm1,pdm1s,jo,ierr,nbmax,
     $               iq0,iq1,ir0,ir1,ibl,id,nqm1
      INTEGER(i4), PARAMETER :: iopt=3,nrhs=1
      REAL(r8), DIMENSION(:), ALLOCATABLE :: btmp,bb
      REAL(r8) :: accdum
      CHARACTER(64) :: msg
      REAL(r8) :: timesolv_st,timesolv_en
c----------------------------------------------------------------------
c     transfer the rhs to 1-vector storage.
c-----------------------------------------------------------------------
      ALLOCATE(bb(0:spp%mloc-1))
      nqm1=nqty-1
      pdm1=pd-1
      pdm1s=pdm1**2
      DO ibl=1,nbl
        mx=fac%mat_info(ibl)%mx
        my=fac%mat_info(ibl)%my
        id=spp%gblock_order(loc2glob(ibl))
        nbmax=SIZE(fac%bl_spp(ibl)%row_ind)
        DO iy=0,my
          DO ix=0,mx
            ir0=fac%bl_spp(ibl)%row_ind(1)%rarr(1,ix,iy)
            IF (ir0<spp%irowst_block(id-1)) CYCLE
            ir0=ir0-spp%fstrow
            bb(ir0:ir0+nqm1)=res(ibl)%arr(:,ix,iy)
          ENDDO
        ENDDO
        IF (nbmax>1) THEN
          DO iy=0,my
            DO ix=1,mx
              DO ib=1,pdm1
                iq0=(ib-1)*nqty+1
                ir0=fac%bl_spp(ibl)%row_ind(2)%rarr(iq0,ix,iy)
                IF (ir0<spp%irowst_block(id-1)) CYCLE
                ir0=ir0-spp%fstrow
                bb(ir0:ir0+nqm1)=res(ibl)%arrh(:,ib,ix,iy)
              ENDDO
            ENDDO
          ENDDO
          DO iy=1,my
            DO ix=0,mx
              DO ib=1,pdm1
                iq0=(ib-1)*nqty+1
                ir0=fac%bl_spp(ibl)%row_ind(3)%rarr(iq0,ix,iy)
                IF (ir0<spp%irowst_block(id-1)) CYCLE
                ir0=ir0-spp%fstrow
                bb(ir0:ir0+nqm1)=res(ibl)%arrv(:,ib,ix,iy)
              ENDDO
            ENDDO
          ENDDO
          IF (nbmax>3) THEN
            DO iy=1,my
              DO ix=1,mx
                DO ib=1,pdm1s
                  iq0=(ib-1)*nqty+1
                  ir0=fac%bl_spp(ibl)%row_ind(4)%rarr(iq0,ix,iy)-
     $                spp%fstrow
                  bb(ir0:ir0+nqm1)=res(ibl)%arri(:,ib,ix,iy)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     communicate if there is more than one processor per layer.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1 .AND. .NOT.spp%matrix_distributed) THEN
        ALLOCATE(btmp(0:spp%nrow-1))
        CALL mpi_allgatherv(bb(spp%algdsplr(node_layer)),
     $       spp%algcntr(node_layer),mpi_nim_real,btmp(0),
     $       spp%algcntr(0),spp%algdsplr(0),mpi_nim_real,
     $       comm_layer,ierr)
        bb=btmp
        DEALLOCATE(btmp)
      ENDIF
c-----------------------------------------------------------------------
c     solve the linear system with the appropriate library.
c-----------------------------------------------------------------------
      CALL timer(timesolv_st)
      IF (precon=='lapack') THEN
        IF (.NOT.single_pr)  THEN
          CALL dpbtrs('L',spp%nrow,spp%nbw,1_i4,fac%a11(1,1),
     $                spp%nbw+1_i4,bb(0),spp%nrow,ierr)
        ELSE
          CALL spbtrs('L',spp%nrow,spp%nbw,1_i4,fac%a11(1,1),
     $                spp%nbw+1_i4,bb(0),spp%nrow,ierr)
        ENDIF
        IF (ierr/=0) THEN
          WRITE (msg,'(a,i9)') 'Solve_dir real: **dpbtrs** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ELSE IF (precon=='seq_slu') THEN
        CALL c_fortran_dgssv(iopt,spp%nrow,spp%nnz,nrhs,accdum,
     $                       spp%j_acc,spp%start_acc,bb,spp%nrow,
     $                       spp%acc_handle,ierr)
        IF (ierr/=0) THEN
          WRITE (msg,'(a,i9)') 'Solve_dir real: **dgssv** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ELSE IF (precon(1:5)=='slu_d') THEN
        IF (precon(1:7)=='slu_dst') THEN
          CALL c_fortran_pdloc(
     $           iopt,spp%nrow,spp%mloc,spp%nnzloc,spp%fstrow,nrhs,
     $           accdum,spp%j_acc(spp%indst),spp%start_loc(0),
     $           bb(0),INT(spp%mloc,i4),spp%acc_handle,
     $           slugrid_handle,spp%iopts,spp%dopts,ierr)
        ELSE
          CALL c_fortran_pdglob(iopt,spp%nrow,spp%nnz,nrhs,accdum,
     $                      spp%j_acc,spp%start_acc,bb,INT(spp%nrow,i4),
     $                      spp%acc_handle,slugrid_handle,ierr)
        ENDIF
        IF (ierr/=0) THEN
          WRITE (msg,'(a,i9)') 'Solve_dir real: **pdloc** ierr=',ierr
          CALL nim_stop(msg)
        ENDIF
      ENDIF
      CALL timer(timesolv_en)
      time_exsolv=time_exsolv+timesolv_en-timesolv_st
c-----------------------------------------------------------------------
c     transfer the solution from 1-vector storage.  distributed
c     memory version first.
c-----------------------------------------------------------------------
      IF (spp%matrix_distributed) THEN
        DO ibl=1,nbl
          mx=fac%mat_info(ibl)%mx
          my=fac%mat_info(ibl)%my
          id=spp%gblock_order(loc2glob(ibl))
          nbmax=SIZE(fac%bl_spp(ibl)%row_ind)
          DO iy=0,my
            DO ix=0,mx
              ir0=fac%bl_spp(ibl)%row_ind(1)%rarr(1,ix,iy)
              IF (ir0< spp%irowst_block(id-1).OR.
     $            ir0>=spp%irowst_block(id).OR.
     $            iy==0.AND.fac%bl_spp(ibl)%perblock) THEN
                zee(ibl)%arr(:,ix,iy)=0._r8
              ELSE
                ir0=ir0-spp%fstrow
                zee(ibl)%arr(:,ix,iy)=bb(ir0:ir0+nqm1)
              ENDIF
            ENDDO
          ENDDO
          IF (nbmax>1) THEN
            DO iy=0,my
              DO ix=1,mx
                DO ib=1,pdm1
                  iq0=(ib-1)*nqty+1
                  ir0=fac%bl_spp(ibl)%row_ind(2)%rarr(iq0,ix,iy)
                  IF (ir0< spp%irowst_block(id-1).OR.
     $                ir0>=spp%irowst_block(id).OR.
     $                iy==0.AND.fac%bl_spp(ibl)%perblock) THEN
                    zee(ibl)%arrh(:,ib,ix,iy)=0._r8
                  ELSE
                    ir0=ir0-spp%fstrow
                    zee(ibl)%arrh(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            DO iy=1,my
              DO ix=0,mx
                DO ib=1,pdm1
                  iq0=(ib-1)*nqty+1
                  ir0=fac%bl_spp(ibl)%row_ind(3)%rarr(iq0,ix,iy)
                  IF (ir0< spp%irowst_block(id-1).OR.
     $                ir0>=spp%irowst_block(id)) THEN
                    zee(ibl)%arrv(:,ib,ix,iy)=0._r8
                  ELSE
                    ir0=ir0-spp%fstrow
                    zee(ibl)%arrv(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            IF (nbmax>3) THEN
              DO iy=1,my
                DO ix=1,mx
                  DO ib=1,pdm1s
                    iq0=(ib-1)*nqty+1
                    ir0=fac%bl_spp(ibl)%row_ind(4)%rarr(iq0,ix,iy)-
     $                  spp%fstrow
                    zee(ibl)%arri(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
          IF (fac%bl_spp(ibl)%degenerate) THEN
            zee(ibl)%arr(:,0,0:my-1)=0._r8
            IF (nbmax>1) zee(ibl)%arrv(:,:,0,:)=0._r8
          ENDIF
          CALL edge_load_arr(zee(ibl),nqty,pdm1,seam(ibl))
        ENDDO
        CALL edge_network(nqty,0_i4,pdm1,.false.)
        DO ibl=1,nbl
          CALL edge_unload_arr(zee(ibl),nqty,pdm1,seam(ibl))
        ENDDO
c-----------------------------------------------------------------------
c     global storage version.
c-----------------------------------------------------------------------
      ELSE
        DO ibl=1,nbl
          mx=fac%mat_info(ibl)%mx
          my=fac%mat_info(ibl)%my
          nbmax=SIZE(fac%bl_spp(ibl)%row_ind)
          DO iy=0,my
            DO ix=0,mx
              ir0=fac%bl_spp(ibl)%row_ind(1)%rarr(1,ix,iy)
              zee(ibl)%arr(:,ix,iy)=bb(ir0:ir0+nqm1)
            ENDDO
          ENDDO
          IF (nbmax>1) THEN
            DO iy=0,my
              DO ix=1,mx
                DO ib=1,pdm1
                  iq0=(ib-1)*nqty+1
                  ir0=fac%bl_spp(ibl)%row_ind(2)%rarr(iq0,ix,iy)
                  zee(ibl)%arrh(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                ENDDO
              ENDDO
            ENDDO
            DO iy=1,my
              DO ix=0,mx
                DO ib=1,pdm1
                  iq0=(ib-1)*nqty+1
                  ir0=fac%bl_spp(ibl)%row_ind(3)%rarr(iq0,ix,iy)
                  zee(ibl)%arrv(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                ENDDO
              ENDDO
            ENDDO
            IF (nbmax>3) THEN
              DO iy=1,my
                DO ix=1,mx
                  DO ib=1,pdm1s
                    iq0=(ib-1)*nqty+1
                    ir0=fac%bl_spp(ibl)%row_ind(4)%rarr(iq0,ix,iy)
                    zee(ibl)%arri(:,ib,ix,iy)=bb(ir0:ir0+nqm1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      DEALLOCATE(bb)
c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE solve_dir
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_real_direct


c-----------------------------------------------------------------------
c     module for line-Jacobi type preconditioners.
c-----------------------------------------------------------------------
      MODULE iter_real_line
      USE local
      USE math_tran
      USE edge_type_mod
      USE matrix_type_mod
      USE factor_type_mod
      USE vector_type_mod
      USE edge
      USE seam_storage_mod
      IMPLICIT NONE

      INTEGER(i4), PRIVATE :: nrb_allpr

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 3. iter_line_fac
c     perform a unidirectional factorization of a matrix at each index
c     perpendicular to the factorization.  the off diagonal elements
c     are expected in the +-1 offset indices of mat (4th index).  if
c     the matrix is periodic, the periodic factors are transfered to the
c     +-2 offset locations.
c
c     note: iter_line_fac and iter_line_solve and the preconditioning
c     schemes that use them have been designed for vector machines.
c     they should not be optimized on scalar machines.
c-----------------------------------------------------------------------
      SUBROUTINE iter_line_fac(nqty,nmx,offlim,off_fac,mat)

      INTEGER(i4), INTENT(IN) :: nqty,nmx,offlim
      REAL(r8), INTENT(IN) :: off_fac
      REAL(r8), DIMENSION(:,:,:,-offlim:,:), INTENT(INOUT) :: mat

      INTEGER(i4) :: iq,jq,ix,nend,ioff
      INTEGER(i4), PARAMETER :: noff=5
      REAL(r8) :: offinc
      REAL(r8), DIMENSION(SIZE(mat,1)) :: off
      REAL(r8), DIMENSION(SIZE(mat,1),1:1,nqty,nqty) :: lu,bb
      REAL(r8), DIMENSION(SIZE(mat,1),1:1,nqty) :: xx
      REAL(r8), DIMENSION(SIZE(mat,1),SIZE(mat,2),SIZE(mat,3),
     $                    -offlim:offlim,SIZE(mat,5)) ::matsav
      LOGICAL :: singular
c-----------------------------------------------------------------------
c     set the loop limit, and transfer periodic elements.
c-----------------------------------------------------------------------
      matsav=mat
      IF (offlim>1) THEN
        nend=nmx-1
        matsav(:,:,:, 2,1    )=matsav(:,:,:,-1,1)
        matsav(:,:,:,-1,1    )=0
        matsav(:,:,:,-2,1    )=matsav(:,:,:, 1,nmx)
        matsav(:,:,:, 1,nmx  )=0
        matsav(:,:,:,-2,nmx-1)=matsav(:,:,:,-1,nmx)
        matsav(:,:,:,-1,nmx  )=0
        matsav(:,:,:, 2,nmx-1)=matsav(:,:,:, 1,nmx-1)
        matsav(:,:,:, 1,nmx-1)=0
      ELSE
        nend=nmx
      ENDIF
      bb=0
      DO iq=1,nqty
        bb(:,:,iq,iq)=1
      ENDDO
c-----------------------------------------------------------------------
c     factor the strided direction.  the upper triangle has nqtyXnqty
c     identity sub-matrices on the diagonal, and the diagonal of the mat
c     array stores the inverse of the lower triangle.
c
c     a loop for testing off-diagonal factors (to find the largest
c     permissible) has been added.
c-----------------------------------------------------------------------
      lu=1
      off=1
      offinc=(1-ABS(off_fac))/noff

      off_loop: DO
        DO ix=1,SIZE(mat,1)
          IF (REAL(lu(ix,1,1,1),r8)<=0) off=off-offinc
        ENDDO
        IF (MINVAL(off)<ABS(off_fac)) EXIT off_loop
        DO ix=1,SIZE(mat,1)
          mat(ix,:,:,:,:)=off(ix)*matsav(ix,:,:,:,:)
        ENDDO
        IF (off_fac>0) THEN
          DO iq=1,nqty
            mat(:,iq,iq,0,:)=matsav(:,iq,iq,0,:)
          ENDDO
        ELSE
          mat(:,:,:,0,:)=matsav(:,:,:,0,:)
        ENDIF
        singular=.false.
        lu(:,1,:,:)=mat(:,:,:,0,1)
        CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,1),'factor',singular)
        IF (singular) CALL nim_stop
     $    ('Iter_line_fac unable to factor first point-block.')
        DO iq=1,nqty
          CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,iq),'solve',singular)
          mat(:,:,iq,0,1)=xx(:,1,:)
        ENDDO
        lu(:,1,:,:)=mat(:,:,:,1,1)
        DO jq=1,nqty
          DO iq=1,nqty
            mat(:,iq,jq,1,1)=SUM(mat(:,iq,:,0,1)*lu(:,1,:,jq),2)
          ENDDO
        ENDDO
        main_loop: DO ix=2,nend
          DO jq=1,nqty
            DO iq=1,nqty
              mat(:,iq,jq,0,ix)=mat(:,iq,jq,0,ix)
     $           -SUM(mat(:,iq,:,-1,ix)*mat(:,:,jq, 1,ix-1),2)
            ENDDO
          ENDDO
          lu(:,1,:,:)=mat(:,:,:,0,ix)
          CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,1),'factor',singular)
          IF (singular) CYCLE off_loop
          DO iq=1,nqty
            CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,iq),
     $                          'solve',singular)
            mat(:,:,iq,0,ix)=xx(:,1,:)
          ENDDO
          lu(:,1,:,:)=mat(:,:,:,1,ix)
          DO jq=1,nqty
            DO iq=1,nqty
              mat(:,iq,jq,1,ix)=SUM(mat(:,iq,:,0,ix)*lu(:,1,:,jq),2)
            ENDDO
          ENDDO
        ENDDO main_loop
c-----------------------------------------------------------------------
c       if the matrix is periodic, the last row of the L factor is
c       placed in the -2 offset, and the last column of the U factor is
c       placed in the +2 offset.
c-----------------------------------------------------------------------
        IF (offlim>1) THEN
          lu(:,1,:,:)=mat(:,:,:,2,1)
          DO jq=1,nqty
            DO iq=1,nqty
              mat(:,iq,jq,2,1)=SUM(mat(:,iq,:,0,1)*lu(:,1,:,jq),2)
            ENDDO
          ENDDO
          DO jq=1,nqty
            DO iq=1,nqty
              mat(:,iq,jq,0,nmx)=mat(:,iq,jq,0,nmx)
     $                -SUM(mat(:,iq,:,-2,1)*mat(:,:,jq, 2,1),2)
            ENDDO
          ENDDO
          per_loop: DO ix=2,nend
            DO jq=1,nqty
              DO iq=1,nqty
                mat(:,iq,jq,-2,ix)=mat(:,iq,jq,-2,ix)
     $                 -SUM(mat(:,iq,:,-2,ix-1)*mat(:,:,jq, 1,ix-1),2)
                mat(:,iq,jq, 2,ix)=mat(:,iq,jq, 2,ix)
     $                 -SUM(mat(:,iq,:,-1,ix)*mat(:,:,jq, 2,ix-1),2)
              ENDDO
            ENDDO
            lu(:,1,:,:)=mat(:,:,:,2,ix)
            DO jq=1,nqty
              DO iq=1,nqty
                mat(:,iq,jq,2,ix)=SUM(mat(:,iq,:,0,ix)*lu(:,1,:,jq),2)
              ENDDO
            ENDDO
            DO jq=1,nqty
              DO iq=1,nqty
                mat(:,iq,jq,0,nmx)=mat(:,iq,jq,0,nmx)
     $                 -SUM(mat(:,iq,:,-2,ix)*mat(:,:,jq, 2,ix),2)
              ENDDO
            ENDDO
          ENDDO per_loop
          lu(:,1,:,:)=mat(:,:,:,0,nmx)
          CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,1),'factor',singular)
          IF (singular) CYCLE off_loop
          DO iq=1,nqty
            CALL math_solve_sym(nqty,lu,xx,bb(:,:,:,iq),
     $                          'solve',singular)
            mat(:,:,iq,0,nmx)=xx(:,1,:)
          ENDDO
        ENDIF
        IF (.NOT.singular) EXIT off_loop
      ENDDO off_loop

      IF (singular) CALL nim_stop
     $  ('Iter_line_fac unable to factor.  Try a smaller off_diag_fac.')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_line_fac
c-----------------------------------------------------------------------
c     subprogram 4. iter_line_solve
c     uni-directional solver.
c-----------------------------------------------------------------------
      SUBROUTINE iter_line_solve(nqty,nmx,offlim,mat,zz)

      INTEGER(i4), INTENT(IN) :: nqty,nmx,offlim
      REAL(r8), DIMENSION(:,:,:,-offlim:,:), INTENT(IN) :: mat
      REAL(r8), DIMENSION(:,:,:), INTENT(INOUT) :: zz

      INTEGER(i4) :: ix,iq,jq,nend
      REAL(r8), DIMENSION(SIZE(mat,1),nqty) :: tmp
c-----------------------------------------------------------------------
c     set the loop limit.
c-----------------------------------------------------------------------
      IF (offlim>1) THEN
        nend=nmx-1
      ELSE
        nend=nmx
      ENDIF
c-----------------------------------------------------------------------
c     Solve Lyy=rhs.  (using zz for yy storage)
c-----------------------------------------------------------------------
      tmp(:,:)=zz(:,:,1)
      DO iq=1,nqty
        zz(:,iq,1)=0
        DO jq=1,nqty
          zz(:,iq,1)=zz(:,iq,1)+mat(:,iq,jq,0,1)*tmp(:,jq)
        ENDDO
      ENDDO
      DO ix=2,nend
        DO iq=1,nqty
          tmp(:,iq)=zz(:,iq,ix)-mat(:,iq,1,-1,ix)*zz(:,1,ix-1)
          DO jq=2,nqty
            tmp(:,iq)=tmp(:,iq)-mat(:,iq,jq,-1,ix)*zz(:,jq,ix-1)
          ENDDO
        ENDDO
        DO iq=1,nqty
          zz(:,iq,ix)=0
          DO jq=1,nqty
            zz(:,iq,ix)=zz(:,iq,ix)+mat(:,iq,jq,0,ix)*tmp(:,jq)
          ENDDO
        ENDDO
      ENDDO
      IF (offlim>1) THEN
        DO iq=1,nqty
          tmp(:,iq)=zz(:,iq,nmx)
          DO ix=1,nend
            DO jq=1,nqty
              tmp(:,iq)=tmp(:,iq)-mat(:,iq,jq,-2,ix)*zz(:,jq,ix)
            ENDDO
          ENDDO
        ENDDO
        DO iq=1,nqty
          zz(:,iq,nmx)=0
          DO jq=1,nqty
            zz(:,iq,nmx)=zz(:,iq,nmx)+mat(:,iq,jq,0,nmx)*tmp(:,jq)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     Solve Uzz=yy.
c-----------------------------------------------------------------------
      IF (offlim>1) THEN
        DO ix=nmx-1,1,-1
          DO iq=1,nqty
            DO jq=1,nqty
              zz(:,iq,ix)=zz(:,iq,ix)-mat(:,iq,jq,1,ix)*zz(:,jq,ix+1)
     $                               -mat(:,iq,jq,2,ix)*zz(:,jq,nmx )
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ix=nmx-1,1,-1
          DO iq=1,nqty
            DO jq=1,nqty
              zz(:,iq,ix)=zz(:,iq,ix)-mat(:,iq,jq,1,ix)*zz(:,jq,ix+1)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_line_solve
c-----------------------------------------------------------------------
c     subprogram 5. iter_ld_fac
c     factor directions of the two-dimensional matrix separately to
c     produce a line-diagonal preconditioner.  if the direction is 'a',
c     factor both directions for an average.
c     the first dimension of the arrays is the perpendicular direction
c     to permit vectorization on appropriate machines.
c
c     'elements' of the factors are nqtyXnqty sub-matrices.
c-----------------------------------------------------------------------
      SUBROUTINE iter_ld_fac(mat,fac,nq,precon,off_diag_fac,poly_deg,
     $                       poly_d2,nrb)

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4), INTENT(IN) :: nq,poly_deg,poly_d2,nrb
      CHARACTER(*), INTENT(IN) :: precon
      REAL(r8), INTENT(IN) :: off_diag_fac

      INTEGER(i4) :: ibl,mxb,myb,iq,jq,ix,iy,iyst,jx,jy,itype,
     $               jst,iof,jof,xad,yad,xst,yst,npdm,ibase,elen
      REAL(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: lc,cl
      REAL(r8), DIMENSION(:,:,:,:,:), POINTER, CONTIGUOUS :: adix,adiy
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS :: lmat
      LOGICAL :: sing
c-----------------------------------------------------------------------
c     check specified direction ordering.
c-----------------------------------------------------------------------
      IF (precon(8:8)/='x'.AND.precon(8:8)/='y'.AND.precon(8:8)/='a')
     $  CALL nim_stop
     $  ('Iter_ld_fac does not recognize direction '//precon(8:8)//'.')
c-----------------------------------------------------------------------
c     get the factors needed for the eliminations.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        CALL iter_ld_elim(mat,fac,nq,precon,poly_deg,poly_d2,nrb)
        npdm=nq*(poly_deg-1)
      ENDIF
c-----------------------------------------------------------------------
c     loop over blocks and set vertex limits.
c-----------------------------------------------------------------------
      block_loop: DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        IF (fac%bl_fac(ibl)%perblock) THEN
          iyst=1
        ELSE
          iyst=0
        ENDIF
c-----------------------------------------------------------------------
c       for the grid x-direction, copy the line to line connections.
c-----------------------------------------------------------------------
        IF (precon(8:8)=='x'.OR.precon(8:8)=='a') THEN
          adix=>fac%bl_fac(ibl)%adix
          DO ibase=1,poly_deg
            IF (ibase==1) THEN
              itype=1
            ELSE
              itype=3
            ENDIF
            yst=mat%rbl_mat(ibl)%iy0(itype)
            yad=(ibase-1)*myb
            lmat=>mat%rbl_mat(ibl)%mat(itype,itype)%arr
            iof=nq*MAX(0_i4,ibase-2_i4)
            DO ix=0,mxb
              DO jx=-1,1
                DO iy=MAX(iyst,yst),myb
                  DO iq=1,nq
                    DO jq=1,nq
                      adix(iy+yad,iq,jq,jx,ix)=
     $                  lmat(iof+jq,jx,0,iof+iq,ix,iy)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         eliminate the line-to-cell and cell-to-line connections.
c-----------------------------------------------------------------------
          IF (poly_deg>1) THEN
            cl=>fac%bl_fac(ibl)%xcl
            lc=>fac%bl_fac(ibl)%xlc
            IF (mat%eliminated) THEN
              elen=myb
            ELSE
              elen=poly_deg*myb
            ENDIF
            DO ix=1,mxb
              DO iq=1,nq
                DO jq=1,nq
                  adix(:elen,iq,jq, 0,ix-1)=adix(:elen,iq,jq, 0,ix-1)-
     $              SUM(cl(:,iq   ,ix,:)*lc(jq   ,:,ix,:),1)
                  adix(:elen,iq,jq, 1,ix-1)=adix(:elen,iq,jq, 1,ix-1)-
     $              SUM(cl(:,iq   ,ix,:)*lc(jq+nq,:,ix,:),1)
                  adix(:elen,iq,jq,-1,ix  )=adix(:elen,iq,jq,-1,ix  )-
     $              SUM(cl(:,iq+nq,ix,:)*lc(jq   ,:,ix,:),1)
                  adix(:elen,iq,jq, 0,ix  )=adix(:elen,iq,jq, 0,ix  )-
     $              SUM(cl(:,iq+nq,ix,:)*lc(jq+nq,:,ix,:),1)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         leaving the degenerate point uncoupled seems to help.  reset
c         the diagonal elements.
c-----------------------------------------------------------------------
          IF (fac%bl_fac(ibl)%degenerate) THEN
            adix(:,:,:, 1,0)=0
            adix(:,:,:,-1,1)=0
          ENDIF
          CALL iter_line_fac(nq,mxb+1_i4,1_i4,off_diag_fac,adix)
        ENDIF
c-----------------------------------------------------------------------
c       for the grid y-direction, copy the line to line connections.
c-----------------------------------------------------------------------
        IF (precon(8:8)=='y'.OR.precon(8:8)=='a') THEN
          adiy=>fac%bl_fac(ibl)%adiy
          DO ibase=1,poly_deg
            IF (ibase==1) THEN
              itype=1
            ELSE
              itype=2
            ENDIF
            xst=mat%rbl_mat(ibl)%ix0(itype)
            xad=(ibase-1)*mxb
            lmat=>mat%rbl_mat(ibl)%mat(itype,itype)%arr
            iof=nq*MAX(0_i4,ibase-2_i4)
            DO iy=iyst,myb
              DO jy=-1,1
                IF (fac%bl_fac(ibl)%perblock.AND.iy==myb.AND.jy==1) THEN
                  DO ix=xst,mxb
                    DO iq=1,nq
                      DO jq=1,nq
                        adiy(ix+xad,iq,jq,jy,iy)=
     $                    lmat(iof+jq,0,jy,iof+iq,ix,0)
                      ENDDO
                    ENDDO
                  ENDDO
                ELSE
                  DO ix=xst,mxb
                    DO iq=1,nq
                      DO jq=1,nq
                        adiy(ix+xad,iq,jq,jy,iy)=
     $                    lmat(iof+jq,0,jy,iof+iq,ix,iy)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         eliminate the line-to-cell and cell-to-line connections.
c-----------------------------------------------------------------------
          IF (poly_deg>1) THEN
            cl=>fac%bl_fac(ibl)%ycl
            lc=>fac%bl_fac(ibl)%ylc
            IF (mat%eliminated) THEN
              elen=mxb
            ELSE
              elen=poly_deg*mxb
            ENDIF
            DO iy=1,myb
              jof=iy-1
              IF (fac%bl_fac(ibl)%perblock.AND.jof==0) jof=myb
              DO iq=1,nq
                DO jq=1,nq
                  adiy(:elen,iq,jq, 0,jof)=adiy(:elen,iq,jq, 0,jof)-
     $              SUM(cl(:,iq   ,:,iy)*lc(jq   ,:,:,iy),1)
                  adiy(:elen,iq,jq, 1,jof)=adiy(:elen,iq,jq, 1,jof)-
     $              SUM(cl(:,iq   ,:,iy)*lc(jq+nq,:,:,iy),1)
                  adiy(:elen,iq,jq,-1,iy)=adiy(:elen,iq,jq,-1,iy)-
     $              SUM(cl(:,iq+nq,:,iy)*lc(jq   ,:,:,iy),1)
                  adiy(:elen,iq,jq, 0,iy)=adiy(:elen,iq,jq, 0,iy)-
     $              SUM(cl(:,iq+nq,:,iy)*lc(jq+nq,:,:,iy),1)
                ENDDO
              ENDDO
            ENDDO
          ENDIF

          IF (fac%bl_fac(ibl)%degenerate) adiy(0,:,:,-1:1:2,:)=0
          CALL iter_line_fac(nq,myb+1_i4-iyst,1_i4+iyst,
     $                       off_diag_fac,adiy)
        ENDIF
      ENDDO block_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_ld_fac
c-----------------------------------------------------------------------
c     subprogram 6. iter_ld_solve
c     perform line-Jacobi solves over the specified direction in each
c     rblock.
c-----------------------------------------------------------------------
      SUBROUTINE iter_ld_solve(bl_fac,zeebl,resbl,nq,precon,poly_deg,
     $                         poly_d2,use_int)

      TYPE(bl_fac_type), INTENT(IN) :: bl_fac
      TYPE(vector_type), INTENT(INOUT) :: zeebl,resbl
      INTEGER(i4), INTENT(IN) :: nq,poly_deg,poly_d2
      CHARACTER(*), INTENT(IN) :: precon
      LOGICAL, INTENT(IN) :: use_int

      INTEGER(i4) :: mxb,myb,iyst,ix,iy,iq,jq,ibase,ibarr,jbase,jbarr,
     $               xad,yad,npdm,iym1,jbar2,perpb_max,jx,jy,ib0,ib1,iv
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: zz,zh
c-----------------------------------------------------------------------
c     block-wise operations.
c-----------------------------------------------------------------------
      mxb=bl_fac%mx
      myb=bl_fac%my
      IF (bl_fac%perblock) THEN
        iyst=1
      ELSE
        iyst=0
      ENDIF
      npdm=nq*(poly_deg-1)
c-----------------------------------------------------------------------
c     if cell-interior data has been eliminated from the matrix equation
c     before calling the iterative solver, it is not used here.
c-----------------------------------------------------------------------
      IF (use_int) THEN
        perpb_max=poly_deg
      ELSE
        perpb_max=1
      ENDIF
c-----------------------------------------------------------------------
c     perform a direct solve over the x-direction of all
c     rblocks, using elimination for connections interior to vertices.
c
c     (only one point along degenerate border has nonzero res, so
c     no averaging is needed.)
c-----------------------------------------------------------------------
      IF (precon(8:8)=='x'.OR.precon(8:8)=='a') THEN
        ALLOCATE(zz(iyst:myb*poly_deg,nq,0:mxb))
        ALLOCATE(zh(npdm,1:mxb,iyst:myb*perpb_max))
        DO ix=0,mxb
          DO iy=iyst,myb
            zz(iy,:,ix)=resbl%arr(:,ix,iy)
          ENDDO
        ENDDO
        IF (poly_deg>1) THEN
          DO ix=0,mxb
            DO ibase=1,poly_deg-1
              yad=ibase*myb
              DO iy=1,myb
                zz(iy+yad,:,ix)=resbl%arrv(:,ibase,ix,iy)
              ENDDO
            ENDDO
          ENDDO
          DO iy=iyst,myb
            DO ix=1,mxb
              DO iq=1,npdm
                zh(iq,ix,iy)=0
                iv=0
                DO ibase=1,poly_deg-1
                  DO jq=1,nq
                    iv=iv+1
                    zh(iq,ix,iy)=zh(iq,ix,iy)+
     $                bl_fac%xelim(iv,iq,ix,iy)*
     $                resbl%arrh(jq,ibase,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (use_int) THEN
            DO jy=1,poly_deg-1
              yad=jy*myb
              ib0=(jy-1)*(poly_deg-1)+1
              ib1=jy*(poly_deg-1)
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,npdm
                    zh(iq,ix,iy+yad)=0
                    iv=0
                    DO ibase=ib0,ib1
                      DO jq=1,nq
                        iv=iv+1
                        zh(iq,ix,iy+yad)=zh(iq,ix,iy+yad)+
     $                      bl_fac%xelim(iv,iq,ix,iy+yad)*
     $                      resbl%arri(jq,ibase,ix,iy)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DO iy=iyst,myb*perpb_max
            DO ix=1,mxb
              DO iq=1,nq
                zz(iy,iq,ix-1)=zz(iy,iq,ix-1)-
     $            SUM(bl_fac%xcl(:,iq   ,ix,iy)*zh(:,ix,iy))
                zz(iy,iq,ix  )=zz(iy,iq,ix  )-
     $            SUM(bl_fac%xcl(:,iq+nq,ix,iy)*zh(:,ix,iy))
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        CALL iter_line_solve(nq,mxb+1_i4,1_i4,bl_fac%adix,zz)

        IF (poly_deg>1) THEN
          DO iy=iyst,myb*perpb_max
            DO ix=1,mxb
              DO jq=1,nq
                zh(:,ix,iy)=zh(:,ix,iy)
     $            -bl_fac%xlc(jq   ,:,ix,iy)*zz(iy,jq,ix-1)
     $            -bl_fac%xlc(jq+nq,:,ix,iy)*zz(iy,jq,ix)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        DO ix=0,mxb
          DO iy=iyst,myb
            zeebl%arr(:,ix,iy)=zz(iy,:,ix)
          ENDDO
        ENDDO
        IF (poly_deg>1) THEN
          DO ix=0,mxb
            DO ibase=1,poly_deg-1
              yad=ibase*myb
              DO iy=1,myb
                zeebl%arrv(:,ibase,ix,iy)=zz(iy+yad,:,ix)
              ENDDO
            ENDDO
          ENDDO
          zeebl%arrh(:,:,:,iyst:myb)=
     $      RESHAPE(zh(:,:,iyst:myb),
     $          (/nq,poly_deg-1_i4,mxb,myb-iyst+1_i4/))
          IF (use_int) THEN
            DO iy=1,myb
              DO ix=1,mxb
                DO ibase=2*poly_deg,poly_d2
                  ibarr=ibase-2*poly_deg+1
                  ib0=nq*(bl_fac%ixm(ibase)-2)+1
                  ib1=nq*(bl_fac%ixm(ibase)-1)
                  yad=(bl_fac%iym(ibase)-1)*myb
                  zeebl%arri(:,ibarr,ix,iy)=zh(ib0:ib1,ix,yad+iy)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        DEALLOCATE(zz,zh)
      ENDIF
c-----------------------------------------------------------------------
c     y-direction:
c-----------------------------------------------------------------------
      IF (precon(8:8)=='y'.OR.precon(8:8)=='a') THEN
        ALLOCATE(zz(0:mxb*poly_deg,nq,iyst:myb))
        ALLOCATE(zh(npdm,0:mxb*perpb_max,1:myb))

        DO iy=iyst,myb
          DO ix=0,mxb
            zz(ix,:,iy)=resbl%arr(:,ix,iy)
          ENDDO
        ENDDO
        IF (poly_deg>1) THEN
          DO iy=iyst,myb
            DO ibase=1,poly_deg-1
              xad=ibase*mxb
              DO ix=1,mxb
                zz(ix+xad,:,iy)=resbl%arrh(:,ibase,ix,iy)
              ENDDO
            ENDDO
          ENDDO
          DO iy=1,myb
            DO ix=0,mxb
              DO iq=1,npdm
                zh(iq,ix,iy)=0
                iv=0
                DO ibase=1,poly_deg-1
                  DO jq=1,nq
                    iv=iv+1
                    zh(iq,ix,iy)=zh(iq,ix,iy)+
     $                bl_fac%yelim(iv,iq,ix,iy)*
     $                resbl%arrv(jq,ibase,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (use_int) THEN
            DO iy=1,myb
              DO jbase=2*poly_deg,poly_d2
                jbarr=jbase-2*poly_deg+1
                jx=bl_fac%ixm(jbase)-1
                ib0=nq*(bl_fac%iym(jbase)-2)+1
                ib1=nq*(bl_fac%iym(jbase)-1)
                xad=jx*mxb
                DO ix=1,mxb
                  DO iq=1,npdm
                    zh(iq,ix+xad,iy)=
     $                SUM(bl_fac%yelim(ib0:ib1,iq,ix+xad,iy)
     $                    *resbl%arri(:,jbarr,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DO iy=1,myb
            iym1=iy-1
            IF (bl_fac%perblock.AND.iym1==0) iym1=myb
            DO ix=0,mxb*perpb_max
              DO iq=1,nq
                zz(ix,iq,iym1)=zz(ix,iq,iym1)-
     $            SUM(bl_fac%ycl(:,iq   ,ix,iy)*zh(:,ix,iy))
                zz(ix,iq,iy  )=zz(ix,iq,iy  )-
     $            SUM(bl_fac%ycl(:,iq+nq,ix,iy)*zh(:,ix,iy))
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        CALL iter_line_solve(nq,myb+1_i4-iyst,iyst+1_i4,bl_fac%adiy,zz)

        IF (poly_deg>1) THEN
          DO iy=1,myb
            iym1=iy-1
            IF (bl_fac%perblock.AND.iym1==0) iym1=myb
            DO ix=0,mxb*perpb_max
              DO jq=1,nq
                zh(:,ix,iy)=zh(:,ix,iy)
     $            -bl_fac%ylc(jq   ,:,ix,iy)*zz(ix,jq,iym1)
     $            -bl_fac%ylc(jq+nq,:,ix,iy)*zz(ix,jq,iy)
              ENDDO
            ENDDO
          ENDDO
        ENDIF

        IF (precon(8:8)=='a') THEN
          DO iy=iyst,myb
            DO ix=0,mxb
              zeebl%arr(:,ix,iy)=0.5*(zeebl%arr(:,ix,iy)+zz(ix,:,iy))
            ENDDO
          ENDDO
          IF (poly_deg>1) THEN
            DO iy=iyst,myb
              DO ibase=1,poly_deg-1
                xad=ibase*mxb
                DO ix=1,mxb
                  zeebl%arrh(:,ibase,ix,iy)=
     $              0.5*(zeebl%arrh(:,ibase,ix,iy)+zz(ix+xad,:,iy))
                ENDDO
              ENDDO
            ENDDO
            zeebl%arrv=0.5*(zeebl%arrv+
     $        RESHAPE(zh(:,0:mxb,:),(/nq,poly_deg-1_i4,mxb+1_i4,myb/)))
            IF (use_int) THEN
              DO iy=1,myb
                DO ix=1,mxb
                  DO ibase=2*poly_deg,poly_d2
                    ibarr=ibase-2*poly_deg+1
                    ib0=nq*(bl_fac%iym(ibase)-2)+1
                    ib1=nq*(bl_fac%iym(ibase)-1)
                    xad=(bl_fac%ixm(ibase)-1)*mxb
                    zeebl%arri(:,ibarr,ix,iy)=0.5*
     $                (zeebl%arri(:,ibarr,ix,iy)+zh(ib0:ib1,xad+ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ELSE
          DO iy=iyst,myb
            DO ix=0,mxb
              zeebl%arr(:,ix,iy)=zz(ix,:,iy)
            ENDDO
          ENDDO
          IF (poly_deg>1) THEN
            DO iy=iyst,myb
              DO ibase=1,poly_deg-1
                xad=ibase*mxb
                DO ix=1,mxb
                  zeebl%arrh(:,ibase,ix,iy)=zz(ix+xad,:,iy)
                ENDDO
              ENDDO
            ENDDO
            zeebl%arrv=
     $        RESHAPE(zh(:,0:mxb,:),(/nq,poly_deg-1_i4,mxb+1_i4,myb/))
            IF (use_int) THEN
              DO iy=1,myb
                DO ix=1,mxb
                  DO ibase=2*poly_deg,poly_d2
                    ibarr=ibase-2*poly_deg+1
                    ib0=nq*(bl_fac%iym(ibase)-2)+1
                    ib1=nq*(bl_fac%iym(ibase)-1)
                    xad=(bl_fac%ixm(ibase)-1)*mxb
                    zeebl%arri(:,ibarr,ix,iy)=zh(ib0:ib1,xad+ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        DEALLOCATE(zz,zh)
      ENDIF

      IF (bl_fac%perblock) THEN
        zeebl%arr(:,:,0)=zeebl%arr(:,:,myb)
        IF (poly_deg>1) zeebl%arrh(:,:,:,0)=zeebl%arrh(:,:,:,myb)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_ld_solve
c-----------------------------------------------------------------------
c     subprogram 7. iter_solve_sym.
c     use Cholesky factorization to solve a symmetric nxn system.
c     if the instruction is 'factor', the lower triangular part of the
c     decomposition is returned in matrix.  if the instruction is
c     'solve', the incoming matrix is expected to be the decomposition.
c     this is similar to the math_tran version, but the array dimensions
c     are different.  there is only one nxn matrix (single grid point),
c     but it is used on multiple rhs's.
c-----------------------------------------------------------------------
      SUBROUTINE iter_solve_sym(n1,n2,lu,xx,bb,instruction,singular)

      INTEGER(i4), INTENT(IN) :: n1,n2
      REAL(r8), DIMENSION(:,:), INTENT(INOUT) :: lu
      REAL(r8), DIMENSION(:,:), INTENT(IN) :: bb
      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: xx
      CHARACTER(*), INTENT(IN) :: instruction
      LOGICAL, INTENT(OUT) :: singular

      REAL(r8), DIMENSION(n1,n2) :: y
      INTEGER(i4) :: jcol,i,k
c-----------------------------------------------------------------------
c     factor matrix into L*transpose(L) (first column is treated first).
c-----------------------------------------------------------------------
      singular=.false.
      IF (instruction/='solve') THEN
        IF (lu(1,1)<=0) THEN
          singular=.true.
          RETURN
        ENDIF
        lu(1,1)=SQRT(lu(1,1))
        DO i=2,n1
          lu(i,1)=lu(i,1)/lu(1,1)
        ENDDO
c-----------------------------------------------------------------------
c       column loop for L.
c-----------------------------------------------------------------------
        DO jcol=2,n1
          lu(jcol,jcol)=lu(jcol,jcol)
     $          -SUM(lu(jcol,1:jcol-1)*lu(jcol,1:jcol-1))
          IF (lu(jcol,jcol)<=0) THEN
            singular=.true.
            RETURN
          ENDIF
          lu(jcol,jcol)=SQRT(lu(jcol,jcol))
          DO i=jcol+1,n1
            lu(i,jcol)=(lu(i,jcol)
     $          -SUM(lu(i,1:jcol-1)*lu(jcol,1:jcol-1)))/lu(jcol,jcol)
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     return if solution is not required.
c-----------------------------------------------------------------------
      IF (instruction=='factor') RETURN
c-----------------------------------------------------------------------
c     solve for xx by getting y first.  {Ly=bb then transpose(L)xx=y}
c-----------------------------------------------------------------------
      y(1,:)=bb(1,:)/lu(1,1)
      DO i=2,n1
        DO k=1,n2
          y(i,k)=(bb(i,k)-SUM(lu(i,1:i-1)*y(1:i-1,k)))/lu(i,i)
        ENDDO
      ENDDO
      xx(n1,:)=y(n1,:)/lu(n1,n1)
      DO i=n1-1,1,-1
        DO k=1,n2
          xx(i,k)=(y(i,k)-SUM(lu(i+1:n1,i)*xx(i+1:n1,k)))/lu(i,i)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_solve_sym
c-----------------------------------------------------------------------
c     subprogram 8. iter_gl_ld_fac
c     factor directions of the two-dimensional matrix separately to
c     produce a line-diagonal preconditioner.  this version uses
c     global lines that extend across the rblocks as much as possible.
c     both directions are always factored with this option.
c
c     'elements' of the factors are nqXnq sub-matrices.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gl_ld_fac(mat,fac,nq,nrb,off_diag_fac,poly_deg,
     $                          poly_d2,ntotb)
      USE pardata
      USE mpi_nim
      USE time

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4), INTENT(IN) :: nq,nrb,poly_deg,poly_d2,ntotb
      REAL(r8), INTENT(IN) :: off_diag_fac

      INTEGER(i4) :: ibl,mxb,myb,iq,jq,ix,iy,iv,jx,jy,npe,npa,past,paen,
     $               ierror,nq2,ibase,imat,yst,xst,yad,xad,perpb_max,
     $               itype,jtype,iof,ibst,iben
      REAL(r8), DIMENSION(:,:,:), POINTER, CONTIGUOUS :: btmp
      REAL(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: lc,cl
      REAL(r8), DIMENSION(:,:,:,:,:), POINTER, CONTIGUOUS :: adi
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER, CONTIGUOUS :: lmat
      REAL(r8) :: timeline_st,timeline_en
      TYPE(vector_type), DIMENSION(nrb) :: blck_tmp,line_tmp
c-----------------------------------------------------------------------
c     first determine the total number of rblocks.
c-----------------------------------------------------------------------
      IF (nprocs_layer>1) THEN
        CALL mpi_allreduce(nrb,nrb_allpr,1,mpi_nim_int,mpi_sum,
     $       comm_layer,ierror)
      ELSE
        nrb_allpr=nrb
      ENDIF
      nq2=nq**2
c-----------------------------------------------------------------------
c     load the x-direction matrix elements into a 3-dimensional array
c     suitable for the line communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        npe=linex(ibl)%perpen-linex(ibl)%perpst+1
        paen=linex(ibl)%paraen
        ALLOCATE(blck_tmp(ibl)%arr(3*nq2,0:mxb,0:myb*poly_deg))
        ALLOCATE(line_tmp(ibl)%arr(npe,3*nq2,0:paen))
        DO ibase=1,poly_deg
          IF (ibase==1) THEN
            itype=1
          ELSE
            itype=3
          ENDIF
          yst=mat%rbl_mat(ibl)%iy0(itype)
          yad=(ibase-1)*myb
          lmat=>mat%rbl_mat(ibl)%mat(itype,itype)%arr
          iof=nq*MAX(0_i4,ibase-2_i4)
          jy=0
          DO jx=-1,1
            DO jq=1,nq
              DO iq=1,nq
                jy=jy+1
                blck_tmp(ibl)%arr(jy,:,yad+yst:yad+myb)=
     $            lmat(iof+jq,jx,0,iof+iq,:,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     get the factors needed for the eliminations, and communicate
c     border contributions.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        CALL iter_ld_elim(mat,fac,nq,"gl_diaga",poly_deg,poly_d2,nrb)
        DO ibl=1,ntotb
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:nq2)=0
            seam(ibl)%segment(iv)%seam_in(1:(poly_deg-1)*nq2)=0
          ENDDO
        ENDDO
        DO ibl=1,nrb
          myb=fac%bl_fac(ibl)%my
          cl=>fac%bl_fac(ibl)%xcl
          lc=>fac%bl_fac(ibl)%xlc
          DO iy=1,myb
            iv=seam(ibl)%nvert-iy
            jy=0
            DO jq=1,nq
              DO iq=1,nq
                jy=jy+1
                seam(ibl)%vertex(iv)%seam_in(jy)=
     $            SUM(cl(:,iq,1,iy)*lc(jq,:,1,iy))
              ENDDO
            ENDDO
            iv=iv+1
            IF (.NOT.mat%eliminated) THEN
              jy=0
              CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                              poly_deg-1_i4,ibst,iben)
              DO ibase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                yad=ibase*myb
                DO jq=1,nq
                  DO iq=1,nq
                    jy=jy+1
                    seam(ibl)%segment(iv)%seam_in(jy)=
     $                SUM(cl(:,iq,1,iy+yad)*lc(jq,:,1,iy+yad))
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          IF (seam(ibl)%segment(1)%ptr(1)<=0.OR.
     $        seam(ibl)%segment(1)%ptr(1)>nrb_allpr) THEN
            iv=seam(ibl)%nvert
            jy=0
            DO jq=1,nq
              DO iq=1,nq
                jy=jy+1
                seam(ibl)%vertex(iv)%seam_in(jy)=
     $            SUM(cl(:,iq,1,0)*lc(jq,:,1,0))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF (mat%eliminated) THEN
          CALL edge_network(nq2,0_i4,0_i4,.false.)
          perpb_max=poly_deg
        ELSE
          CALL edge_network(nq2,0_i4,poly_deg-1_i4,.false.)
          perpb_max=2*poly_deg-1
        ENDIF
c-----------------------------------------------------------------------
c       now modify the matrix with the eliminations including
c       communicated contributions for the right side.
c-----------------------------------------------------------------------
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          cl=>fac%bl_fac(ibl)%xcl
          lc=>fac%bl_fac(ibl)%xlc
          btmp=>blck_tmp(ibl)%arr
          DO ibase=poly_deg,perpb_max
            imat=ibase
            IF (imat==poly_deg) imat=1
            IF (imat==1) THEN
              yst=0
            ELSE
              yst=1
            ENDIF
            yad=(fac%bl_fac(ibl)%iym(imat)-1)*myb
            DO iy=yad+yst,yad+myb
              DO ix=1,mxb
                jy=0
                DO jq=1,nq
                  DO iq=1,nq
                    jy=jy+1
                    btmp(jy+  nq2,ix-1,iy)=btmp(jy+  nq2,ix-1,iy)-
     $                SUM(cl(:,iq   ,ix,iy)*lc(jq   ,:,ix,iy))
                    btmp(jy+2*nq2,ix-1,iy)=btmp(jy+2*nq2,ix-1,iy)-
     $                SUM(cl(:,iq   ,ix,iy)*lc(jq+nq,:,ix,iy))
                    btmp(jy      ,ix  ,iy)=btmp(jy      ,ix  ,iy)-
     $                SUM(cl(:,iq+nq,ix,iy)*lc(jq   ,:,ix,iy))
                    btmp(jy+  nq2,ix  ,iy)=btmp(jy+  nq2,ix  ,iy)-
     $                SUM(cl(:,iq+nq,ix,iy)*lc(jq+nq,:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO iv=mxb,mxb+myb
            iy=iv-mxb
            btmp(nq2+1:2*nq2,mxb,iy)=btmp(nq2+1:2*nq2,mxb,iy)-
     $        seam(ibl)%vertex(iv)%seam_out(1:nq2)
          ENDDO
          IF (.NOT.mat%eliminated) THEN
            CALL edge_load_limits(seam(ibl)%segment(mxb+1)%load_dir,
     $                            poly_deg-1_i4,ibst,iben)
            jq=0
            DO ibase=ibst,iben,seam(ibl)%segment(mxb+1)%load_dir
              yad=ibase*myb
              DO iv=mxb+1,mxb+myb
                iy=iv-mxb
                btmp(nq2+1:2*nq2,mxb,iy+yad)=
     $            btmp(nq2+1:2*nq2,mxb,iy+yad)-
     $              seam(ibl)%segment(iv)%seam_out(jq+1:jq+nq2)
              ENDDO
              jq=jq+nq2
            ENDDO
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate off-diagonals extending one vertex into the start of
c     the block to the end of the previous block.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_in(1:nq2)=0
          seam(ibl)%segment(iv)%seam_in(1:(poly_deg-1)*nq2)=0
        ENDDO
      ENDDO
      DO ibl=1,nrb
        btmp=>blck_tmp(ibl)%arr
        myb=fac%bl_fac(ibl)%my
        DO iy=1,myb
          iv=seam(ibl)%nvert-iy
          seam(ibl)%vertex(iv)%seam_in(1:nq2)=btmp(2*nq2+1:3*nq2,0,iy)
        ENDDO
        IF (seam(ibl)%segment(1)%ptr(1)<=0.OR.
     $      seam(ibl)%segment(1)%ptr(1)>nrb_allpr) THEN
          iv=seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_in(1:nq2)=btmp(2*nq2+1:3*nq2,0,0)
        ENDIF
        CALL edge_load_limits(seam(ibl)%segment(seam(ibl)%nvert)%
     $                        load_dir,poly_deg-1_i4,ibst,iben)
        jq=0
        DO ibase=ibst,iben,seam(ibl)%segment(seam(ibl)%nvert)%load_dir
          DO iy=ibase*myb+1,(ibase+1)*myb
            iv=seam(ibl)%nvert-iy+ibase*myb+1
            seam(ibl)%segment(iv)%seam_in(jq+1:jq+nq2)=
     $        btmp(2*nq2+1:3*nq2,0,iy)
          ENDDO
          jq=jq+nq2
        ENDDO
      ENDDO
      CALL edge_network(nq2,0_i4,poly_deg-1_i4,.false.)
c-----------------------------------------------------------------------
c     transfer the off-diagonal connection to the right side.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        btmp=>blck_tmp(ibl)%arr
        DO iy=0,myb
          iv=mxb+iy
          btmp(2*nq2+1:3*nq2,mxb,iy)=
     $       seam(ibl)%vertex(iv)%seam_out(1:nq2)
        ENDDO
        CALL edge_load_limits(seam(ibl)%segment(mxb+1)%load_dir,
     $                        poly_deg-1_i4,ibst,iben)
        jq=0
        DO ibase=ibst,iben,seam(ibl)%segment(mxb+1)%load_dir
          yad=ibase*myb
          DO iv=mxb+1,mxb+myb
            iy=iv-mxb
            btmp(2*nq2+1:3*nq2,mxb,iy+yad)=
     $        seam(ibl)%segment(iv)%seam_out(jq+1:jq+nq2)
          ENDDO
          jq=jq+nq2
        ENDDO
        IF (fac%bl_fac(ibl)%degenerate) THEN
          btmp(2*nq2+1:3*nq2,0,:)=0
          btmp(      1:  nq2,1,:)=0
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     transfer to line format then to the 5D line-mat array.  finally,
c     factor the lines.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      CALL parallel_line_comm_bl2line(blck_tmp,line_tmp,nrb,ntotb,
     $                                3_i4*nq2,linex)
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
      DO ibl=1,nrb
        past=linex(ibl)%parast
        npa=linex(ibl)%paraen-past+1
        npe=linex(ibl)%perpen-linex(ibl)%perpst+1
        adi=>fac%bl_fac(ibl)%adix
        adi(:,:,:,-1:1,:)=
     $    RESHAPE(line_tmp(ibl)%arr(:,:,past:),(/npe,nq,nq,3_i4,npa/))
        DEALLOCATE(blck_tmp(ibl)%arr,line_tmp(ibl)%arr)
        CALL iter_line_fac(nq,npa,past+1_i4,off_diag_fac,adi)
      ENDDO

c-----------------------------------------------------------------------
c     load the y-direction matrix elements into a 3-dimensional array
c     suitable for the line communication.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        npe=liney(ibl)%perpen-liney(ibl)%perpst+1
        paen=liney(ibl)%paraen
        ALLOCATE(blck_tmp(ibl)%arr(3*nq2,0:mxb*poly_deg,0:myb))
        ALLOCATE(line_tmp(ibl)%arr(npe,3*nq2,0:paen))
        DO ibase=1,poly_deg
          IF (ibase==1) THEN
            itype=1
          ELSE
            itype=2
          ENDIF
          xst=mat%rbl_mat(ibl)%ix0(itype)
          xad=(ibase-1)*mxb
          lmat=>mat%rbl_mat(ibl)%mat(itype,itype)%arr
          iof=nq*MAX(0_i4,ibase-2_i4)
          jx=0
          DO jy=-1,1
            DO jq=1,nq
              DO iq=1,nq
                jx=jx+1
                blck_tmp(ibl)%arr(jx,xad+xst:xad+mxb,:)=
     $                    lmat(iof+jq,0,jy,iof+iq,:,:)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     get the factors needed for the eliminations, and communicate
c     border contributions.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        DO ibl=1,ntotb
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:nq2)=0
            seam(ibl)%segment(iv)%seam_in(1:(poly_deg-1)*nq2)=0
          ENDDO
        ENDDO
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          cl=>fac%bl_fac(ibl)%ycl
          lc=>fac%bl_fac(ibl)%ylc
          DO ix=1,mxb
            jx=0
            DO jq=1,nq
              DO iq=1,nq
                jx=jx+1
                seam(ibl)%vertex(ix)%seam_in(jx)=
     $            SUM(cl(:,iq,ix,1)*lc(jq,:,ix,1))
              ENDDO
            ENDDO
            IF (.NOT.mat%eliminated) THEN
              jx=0
              CALL edge_load_limits(seam(ibl)%segment(ix)%load_dir,
     $                              poly_deg-1_i4,ibst,iben)
              DO ibase=ibst,iben,seam(ibl)%segment(ix)%load_dir
                xad=ibase*mxb
                DO jq=1,nq
                  DO iq=1,nq
                    jx=jx+1
                    seam(ibl)%segment(ix)%seam_in(jx)=
     $                SUM(cl(:,iq,ix+xad,1)*lc(jq,:,ix+xad,1))
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          iv=seam(ibl)%nvert
          IF (seam(ibl)%segment(iv)%ptr(1)<=0.OR.
     $        seam(ibl)%segment(iv)%ptr(1)>nrb_allpr) THEN
            jx=0
            DO jq=1,nq
              DO iq=1,nq
                jx=jx+1
                seam(ibl)%vertex(iv)%seam_in(jx)=
     $            SUM(cl(:,iq,0,1)*lc(jq,:,0,1))
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF (mat%eliminated) THEN
          CALL edge_network(nq2,0_i4,0_i4,.false.)
          perpb_max=1
        ELSE
          CALL edge_network(nq2,0_i4,poly_deg-1_i4,.false.)
          perpb_max=poly_deg
        ENDIF
c-----------------------------------------------------------------------
c       now modify the matrix with the eliminations including
c       communicated contributions for the top.
c-----------------------------------------------------------------------
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          cl=>fac%bl_fac(ibl)%ycl
          lc=>fac%bl_fac(ibl)%ylc
          btmp=>blck_tmp(ibl)%arr
          DO imat=1,perpb_max
            IF (imat==1) THEN
              xst=0
            ELSE
              xst=1
            ENDIF
            xad=(imat-1)*mxb
            DO ix=xad+xst,xad+mxb
              DO iy=1,myb
                jx=0
                DO jq=1,nq
                  DO iq=1,nq
                    jx=jx+1
                    btmp(jx+  nq2,ix,iy-1)=btmp(jx+  nq2,ix,iy-1)-
     $                SUM(cl(:,iq   ,ix,iy)*lc(jq   ,:,ix,iy))
                    btmp(jx+2*nq2,ix,iy-1)=btmp(jx+2*nq2,ix,iy-1)-
     $                SUM(cl(:,iq   ,ix,iy)*lc(jq+nq,:,ix,iy))
                    btmp(jx      ,ix,iy  )=btmp(jx      ,ix,iy  )-
     $                SUM(cl(:,iq+nq,ix,iy)*lc(jq   ,:,ix,iy))
                    btmp(jx+  nq2,ix,iy  )=btmp(jx+  nq2,ix,iy  )-
     $                SUM(cl(:,iq+nq,ix,iy)*lc(jq+nq,:,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO iv=mxb+myb,2*mxb+myb
            ix=2*mxb+myb-iv
            btmp(nq2+1:2*nq2,ix,myb)=btmp(nq2+1:2*nq2,ix,myb)-
     $        seam(ibl)%vertex(iv)%seam_out(1:nq2)
          ENDDO
          IF (.NOT.mat%eliminated) THEN
            CALL edge_load_limits(seam(ibl)%segment(2*mxb+myb)%load_dir,
     $                            poly_deg-1_i4,ibst,iben)
            jq=0
            DO ibase=ibst,iben,seam(ibl)%segment(2*mxb+myb)%load_dir
              xad=ibase*mxb
              DO iv=mxb+myb+1,2*mxb+myb
                ix=2*mxb+myb+1-iv
                btmp(nq2+1:2*nq2,ix+xad,myb)=
     $            btmp(nq2+1:2*nq2,ix+xad,myb)-
     $              seam(ibl)%segment(iv)%seam_out(jq+1:jq+nq2)
              ENDDO
              jq=jq+nq2
            ENDDO
          ENDIF
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate off-diagonals extending one vertex into the start of
c     the block to the end of the previous block.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_in(1:nq2)=0
          seam(ibl)%segment(iv)%seam_in(1:(poly_deg-1)*nq2)=0
        ENDDO
      ENDDO
      DO ibl=1,nrb
        btmp=>blck_tmp(ibl)%arr
        mxb=fac%bl_fac(ibl)%mx
        DO ix=1,mxb
          seam(ibl)%vertex(ix)%seam_in(1:nq2)=btmp(2*nq2+1:3*nq2,ix,0)
        ENDDO
        iv=seam(ibl)%nvert
        IF (seam(ibl)%segment(iv)%ptr(1)<=0.OR.
     $      seam(ibl)%segment(iv)%ptr(1)>nrb_allpr) THEN
          seam(ibl)%vertex(iv)%seam_in(1:nq2)=btmp(2*nq2+1:3*nq2,0,0)
        ENDIF
        CALL edge_load_limits(seam(ibl)%segment(1)%load_dir,
     $                        poly_deg-1_i4,ibst,iben)
        jq=0
        DO ibase=ibst,iben,seam(ibl)%segment(1)%load_dir
          DO iv=1,mxb
            ix=ibase*mxb+iv
            seam(ibl)%segment(iv)%seam_in(jq+1:jq+nq2)=
     $        btmp(2*nq2+1:3*nq2,ix,0)
          ENDDO
          jq=jq+nq2
        ENDDO
      ENDDO
      CALL edge_network(nq2,0_i4,poly_deg-1_i4,.false.)
c-----------------------------------------------------------------------
c     transfer the off-diagonal connection to the top.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        btmp=>blck_tmp(ibl)%arr
        DO iv=mxb+myb,2*mxb+myb
          ix=2*mxb+myb-iv
          btmp(2*nq2+1:3*nq2,ix,myb)=
     $       seam(ibl)%vertex(iv)%seam_out(1:nq2)
        ENDDO
        CALL edge_load_limits(seam(ibl)%segment(2*mxb+myb)%load_dir,
     $                        poly_deg-1_i4,ibst,iben)
        jq=0
        DO ibase=ibst,iben,seam(ibl)%segment(2*mxb+myb)%load_dir
          xad=ibase*mxb
          DO iv=mxb+myb+1,2*mxb+myb
            ix=2*mxb+myb+1-iv
            btmp(2*nq2+1:3*nq2,ix+xad,myb)=
     $        seam(ibl)%segment(iv)%seam_out(jq+1:jq+nq2)
          ENDDO
          jq=jq+nq2
        ENDDO
        IF (fac%bl_fac(ibl)%degenerate) THEN
          btmp(      1:  nq2,0,:)=0
          btmp(2*nq2+1:3*nq2,0,:)=0
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     transfer to line format then to the 5D line-mat array.  finally,
c     factor the lines.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      CALL parallel_line_comm_bl2line(blck_tmp,line_tmp,nrb,ntotb,
     $                                3_i4*nq2,liney)
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
      DO ibl=1,nrb
        past=liney(ibl)%parast
        npa=liney(ibl)%paraen-past+1
        npe=liney(ibl)%perpen-liney(ibl)%perpst+1
        adi=>fac%bl_fac(ibl)%adiy
        adi(:,:,:,-1:1,:)=
     $    RESHAPE(line_tmp(ibl)%arr(:,:,past:),(/npe,nq,nq,3_i4,npa/))
        DEALLOCATE(blck_tmp(ibl)%arr,line_tmp(ibl)%arr)
        CALL iter_line_fac(nq,npa,past+1_i4,off_diag_fac,adi)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gl_ld_fac
c-----------------------------------------------------------------------
c     subprogram 9. iter_gl_ld_solve
c     perform rblock-global line Jacobi solves over both directions.
c-----------------------------------------------------------------------
      SUBROUTINE iter_gl_ld_solve(fac,nq,nrb,zee,res,adr,poly_deg,
     $                            poly_d2,ntotb,use_int)
      USE pardata
      USE time

      TYPE(matrix_factor_type), INTENT(IN) :: fac
      TYPE(vector_type), DIMENSION(ntotb), INTENT(INOUT) :: res,zee,adr
      INTEGER(i4), INTENT(IN) :: nq,nrb,poly_deg,poly_d2,ntotb
      LOGICAL, INTENT(IN) :: use_int

      INTEGER(i4) :: ibl,mxb,myb,iyst,ix,iy,iq,npe,past,paen,ibase,
     $               jbase,ibarr,xad,yad,jq,jof,iof,iv,nv,perpb_max,
     $               itype,jtype,ib0,ib1,jx,jy,jbarr,npdm,ibst,iben
      TYPE(vector_type), DIMENSION(nrb) :: zz,rx,rxh,ry,ryv
      REAL(r8) :: timeline_st,timeline_en
c-----------------------------------------------------------------------
c     if cell-interior data has been eliminated from the matrix equation
c     before calling the iterative solver, it is not used here.
c-----------------------------------------------------------------------
      IF (use_int) THEN
        perpb_max=poly_deg
      ELSE
        perpb_max=1
      ENDIF
      npdm=nq*(poly_deg-1)
c-----------------------------------------------------------------------
c     set-up block data for x-direction solves.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        npe=linex(ibl)%perpen-linex(ibl)%perpst+1
        paen=linex(ibl)%paraen
        ALLOCATE(zz(ibl)%arr(npe,nq,0:paen))
c-----------------------------------------------------------------------
c       eliminate data in the 1D cell interiors for the x-direction.
c-----------------------------------------------------------------------
        IF (poly_deg>1) THEN
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          ALLOCATE(rx(ibl)%arr(nq,0:mxb,0:myb*poly_deg))
          ALLOCATE(rxh(ibl)%arr(npdm,1:mxb,0:myb*perpb_max))
          rx(ibl)%arr(:,:,0:myb)=res(ibl)%arr
          DO ibase=1,poly_deg-1
            yad=ibase*myb
            rx(ibl)%arr(:,:,yad+1:yad+myb)=res(ibl)%arrv(:,ibase,:,:)
          ENDDO
          DO iy=0,myb
            DO ix=1,mxb
              DO iq=1,npdm
                rxh(ibl)%arr(iq,ix,iy)=0
                iv=0
                DO ibase=1,poly_deg-1
                  DO jq=1,nq
                    iv=iv+1
                    rxh(ibl)%arr(iq,ix,iy)=rxh(ibl)%arr(iq,ix,iy)+
     $                fac%bl_fac(ibl)%xelim(iv,iq,ix,iy)*
     $                res(ibl)%arrh(jq,ibase,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (use_int) THEN
            DO jy=1,poly_deg-1
              yad=jy*myb
              ib0=(jy-1)*(poly_deg-1)+1
              ib1=jy*(poly_deg-1)
              DO iy=1,myb
                DO ix=1,mxb
                  DO iq=1,npdm
                    rxh(ibl)%arr(iq,ix,iy+yad)=0
                    iv=0
                    DO ibase=ib0,ib1
                      DO jq=1,nq
                        iv=iv+1
                        rxh(ibl)%arr(iq,ix,iy+yad)=
     $                    rxh(ibl)%arr(iq,ix,iy+yad)+
     $                      fac%bl_fac(ibl)%xelim(iv,iq,ix,iy+yad)*
     $                      res(ibl)%arri(jq,ibase,ix,iy)
                      ENDDO
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DO iy=0,myb*perpb_max
            DO ix=1,mxb-1
              DO jq=1,nq
                rx(ibl)%arr(jq,ix,iy)=rx(ibl)%arr(jq,ix,iy)-
     $            SUM(fac%bl_fac(ibl)%xcl(:,jq+nq,ix  ,iy)*
     $                rxh(ibl)%arr(:,ix  ,iy)
     $               +fac%bl_fac(ibl)%xcl(:,jq   ,ix+1,iy)*
     $                rxh(ibl)%arr(:,ix+1,iy))
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         load eliminations along the left and right borders for the
c         x-direction solves into the seams.
c-----------------------------------------------------------------------
          xad=1
          IF (seam(ibl)%segment(1)%ptr(1)<=0.OR.
     $        seam(ibl)%segment(1)%ptr(1)>nrb_allpr) xad=0
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            seam(ibl)%vertex(iv)%seam_in(1:2*nq)=0
            IF (use_int)
     $        seam(ibl)%segment(iv)%seam_in(1:2*nq*(poly_deg-1))=0
            IF (ix==0.AND.iy>=xad) THEN
              DO jq=1,nq
                seam(ibl)%vertex(iv)%seam_in(jq)=
     $            seam(ibl)%vertex(iv)%seam_in(jq)+
     $            SUM(fac%bl_fac(ibl)%xcl(:,jq,1,iy)*
     $                rxh(ibl)%arr(:,1,iy))
              ENDDO
            ELSE IF (ix==mxb.AND.iy>=xad) THEN
              DO jq=1,nq
                seam(ibl)%vertex(iv)%seam_in(jq)=
     $            seam(ibl)%vertex(iv)%seam_in(jq)+
     $            SUM(fac%bl_fac(ibl)%xcl(:,jq+nq,mxb,iy)*
     $                rxh(ibl)%arr(:,mxb,iy))
              ENDDO
            ENDIF
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                            poly_deg-1_i4,ibst,iben)
            IF (use_int.AND..NOT.seam(ibl)%segment(iv)%h_side) THEN
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              IF (ix==0) THEN
                jof=0
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  yad=myb*jbase
                  DO jq=1,nq
                    seam(ibl)%segment(iv)%seam_in(jof+jq)=
     $                seam(ibl)%segment(iv)%seam_in(jof+jq)+
     $                SUM(fac%bl_fac(ibl)%xcl(:,jq,1,yad+iy)*
     $                    rxh(ibl)%arr(:,1,yad+iy))
                  ENDDO
                  jof=jof+2*nq
                ENDDO
              ELSE
                jof=0
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  yad=myb*jbase
                  DO jq=1,nq
                    seam(ibl)%segment(iv)%seam_in(jof+jq)=
     $                seam(ibl)%segment(iv)%seam_in(jof+jq)+
     $                SUM(fac%bl_fac(ibl)%xcl(:,jq+nq,mxb,yad+iy)*
     $                    rxh(ibl)%arr(:,mxb,yad+iy))
                  ENDDO
                  jof=jof+2*nq
                ENDDO
              ENDIF
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c         eliminate data in the 1D cell interiors for the y-direction
c         solves.  starting the y-direction reduces the seam
c         communication by one call per cg iteration.
c-----------------------------------------------------------------------
          ALLOCATE(ry(ibl)%arr(nq,0:mxb*poly_deg,0:myb))
          ALLOCATE(ryv(ibl)%arr(npdm,0:mxb*perpb_max,1:myb))
          ry(ibl)%arr(:,0:mxb,:)=res(ibl)%arr
          DO ibase=1,poly_deg-1
            xad=ibase*mxb
            ry(ibl)%arr(:,xad+1:xad+mxb,:)=res(ibl)%arrh(:,ibase,:,:)
          ENDDO
          DO iy=1,myb
            DO ix=0,mxb
              DO iq=1,npdm
                ryv(ibl)%arr(iq,ix,iy)=0
                iv=0
                DO ibase=1,poly_deg-1
                  DO jq=1,nq
                    iv=iv+1
                    ryv(ibl)%arr(iq,ix,iy)=ryv(ibl)%arr(iq,ix,iy)+
     $                fac%bl_fac(ibl)%yelim(iv,iq,ix,iy)*
     $                res(ibl)%arrv(jq,ibase,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          IF (use_int) THEN
            DO iy=1,myb
              DO jbase=2*poly_deg,poly_d2
                jbarr=jbase-2*poly_deg+1
                jx=fac%bl_fac(ibl)%ixm(jbase)-1
                ib0=nq*(fac%bl_fac(ibl)%iym(jbase)-2)+1
                ib1=nq*(fac%bl_fac(ibl)%iym(jbase)-1)
                xad=jx*mxb
                DO ix=1,mxb
                  DO iq=1,npdm
                    ryv(ibl)%arr(iq,ix+xad,iy)=
     $                SUM(fac%bl_fac(ibl)%yelim(ib0:ib1,iq,ix+xad,iy)
     $                    *res(ibl)%arri(:,jbarr,ix,iy))
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DO iy=1,myb-1
            DO ix=0,mxb*perpb_max
              DO jq=1,nq
                ry(ibl)%arr(jq,ix,iy)=ry(ibl)%arr(jq,ix,iy)-
     $            SUM(fac%bl_fac(ibl)%ycl(:,jq+nq,ix,iy  )*
     $                ryv(ibl)%arr(:,ix,iy  )
     $               +fac%bl_fac(ibl)%ycl(:,jq   ,ix,iy+1)*
     $                ryv(ibl)%arr(:,ix,iy+1))
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         load eliminations along the bottom and top borders into the
c         seams for the y-direction solves.
c-----------------------------------------------------------------------
          yad=1
          nv=seam(ibl)%nvert
          IF (seam(ibl)%segment(nv)%ptr(1)<=0.OR.
     $        seam(ibl)%segment(nv)%ptr(1)>nrb_allpr) yad=0
          DO iv=1,nv
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            IF (iy==0.AND.ix>=yad) THEN
              DO jq=1,nq
                seam(ibl)%vertex(iv)%seam_in(jq+nq)=
     $            seam(ibl)%vertex(iv)%seam_in(jq+nq)+
     $            SUM(fac%bl_fac(ibl)%ycl(:,jq,ix,1)*
     $                ryv(ibl)%arr(:,ix,1))
              ENDDO
            ELSE IF (iy==myb.AND.ix>=yad) THEN
              DO jq=1,nq
                seam(ibl)%vertex(iv)%seam_in(jq+nq)=
     $            seam(ibl)%vertex(iv)%seam_in(jq+nq)+
     $            SUM(fac%bl_fac(ibl)%ycl(:,jq+nq,ix,myb)*
     $                ryv(ibl)%arr(:,ix,myb))
              ENDDO
            ENDIF
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                            poly_deg-1_i4,ibst,iben)
            IF (use_int.AND.seam(ibl)%segment(iv)%h_side) THEN
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              IF (iy==0) THEN
                jof=nq
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  xad=mxb*jbase
                  DO jq=1,nq
                    seam(ibl)%segment(iv)%seam_in(jof+jq)=
     $                seam(ibl)%segment(iv)%seam_in(jof+jq)+
     $                SUM(fac%bl_fac(ibl)%ycl(:,jq,ix+xad,1)*
     $                    ryv(ibl)%arr(:,xad+ix,1))
                  ENDDO
                  jof=jof+2*nq
                ENDDO
              ELSE
                jof=nq
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  xad=mxb*jbase
                  DO jq=1,nq
                    seam(ibl)%segment(iv)%seam_in(jof+jq)=
     $                seam(ibl)%segment(iv)%seam_in(jof+jq)+
     $                SUM(fac%bl_fac(ibl)%ycl(:,jq+nq,ix+xad,myb)*
     $                    ryv(ibl)%arr(:,xad+ix,myb))
                  ENDDO
                  jof=jof+2*nq
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     sum border elimination contributions, and distribute.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        DO ibl=nrb+1,ntotb
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:2*nq)=0
            IF (use_int)
     $        seam(ibl)%segment(iv)%seam_in(1:2*nq*(poly_deg-1))=0
          ENDDO
        ENDDO
        IF (use_int) THEN
          CALL edge_network(2_i4*nq,0_i4,poly_deg-1_i4,.false.)
        ELSE
          CALL edge_network(2_i4*nq,0_i4,0_i4,.false.)
        ENDIF
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            IF (ix==0.OR.ix==mxb) THEN
              rx(ibl)%arr(:,ix,iy)=rx(ibl)%arr(:,ix,iy)-
     $          seam(ibl)%vertex(iv)%seam_out(1:nq)
            ENDIF
            IF (iy==0.OR.iy==myb) THEN
              ry(ibl)%arr(:,ix,iy)=ry(ibl)%arr(:,ix,iy)-
     $          seam(ibl)%vertex(iv)%seam_out(nq+1:2*nq)
            ENDIF
            IF (use_int) THEN
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,
     $                              poly_deg-1_i4,ibst,iben)
              IF (seam(ibl)%segment(iv)%h_side) THEN
                jof=nq
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  xad=mxb*jbase
                  ry(ibl)%arr(:,xad+ix,iy)=ry(ibl)%arr(:,xad+ix,iy)-
     $               seam(ibl)%segment(iv)%seam_out(jof+1:jof+nq)
                  jof=jof+2*nq
                ENDDO
              ELSE
                jof=0
                DO jbase=ibst,iben,seam(ibl)%segment(iv)%load_dir
                  yad=myb*jbase
                  rx(ibl)%arr(:,ix,yad+iy)=rx(ibl)%arr(:,ix,yad+iy)-
     $               seam(ibl)%segment(iv)%seam_out(jof+1:jof+nq)
                  jof=jof+2*nq
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     communicate the residual to global lines for the x-direction
c     solves.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      IF (poly_deg==1) THEN
        CALL parallel_line_comm_bl2line(res,zz,nrb,ntotb,nq,linex)
      ELSE
        CALL parallel_line_comm_bl2line(rx,zz,nrb,ntotb,nq,linex)
      ENDIF
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
c-----------------------------------------------------------------------
c     perform direct solves over the x-direction of all rblocks.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        past=linex(ibl)%parast
        paen=linex(ibl)%paraen
        CALL iter_line_solve(nq,paen-past+1_i4,past+1_i4,
     $    fac%bl_fac(ibl)%adix,zz(ibl)%arr(:,:,past:))
        IF (past>0) zz(ibl)%arr(:,:,0)=zz(ibl)%arr(:,:,paen)
      ENDDO
c-----------------------------------------------------------------------
c     communicate back to block organization.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      IF (poly_deg==1) THEN
        CALL parallel_line_comm_line2bl(zee,zz,nrb,ntotb,nq,linex)
      ELSE
        CALL parallel_line_comm_line2bl(rx,zz,nrb,ntotb,nq,linex)
      ENDIF
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
c-----------------------------------------------------------------------
c     compute 1D cell interior data.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          DO iy=0,myb*perpb_max
            DO ix=1,mxb
              DO jq=1,npdm
                rxh(ibl)%arr(jq,ix,iy)=rxh(ibl)%arr(jq,ix,iy)-
     $              SUM(fac%bl_fac(ibl)%xlc(  :nq,jq,ix,iy)*
     $                  rx(ibl)%arr(:,ix-1,iy)
     $                 +fac%bl_fac(ibl)%xlc(nq+1:,jq,ix,iy)*
     $                  rx(ibl)%arr(:,ix  ,iy))
              ENDDO
            ENDDO
          ENDDO
          zee(ibl)%arr=rx(ibl)%arr(:,:,0:myb)
          DO ibase=1,poly_deg-1
            yad=ibase*myb
            zee(ibl)%arrv(:,ibase,:,:)=rx(ibl)%arr(:,:,yad+1:yad+myb)
          ENDDO
          zee(ibl)%arrh=RESHAPE(rxh(ibl)%arr(:,:,0:myb),
     $                          (/nq,poly_deg-1_i4,mxb,myb+1_i4/))
          IF (use_int) THEN
            DO iy=1,myb
              DO ix=1,mxb
                DO ibase=2*poly_deg,poly_d2
                  ibarr=ibase-2*poly_deg+1
                  ib0=nq*(fac%bl_fac(ibl)%ixm(ibase)-2)+1
                  ib1=nq*(fac%bl_fac(ibl)%ixm(ibase)-1)
                  yad=(fac%bl_fac(ibl)%iym(ibase)-1)*myb
                  zee(ibl)%arri(:,ibarr,ix,iy)=
     $              rxh(ibl)%arr(ib0:ib1,ix,yad+iy)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DEALLOCATE(rx(ibl)%arr,rxh(ibl)%arr)
        ENDDO
      ENDIF

c-----------------------------------------------------------------------
c     set-up block data for y-direction solves.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        npe=liney(ibl)%perpen-liney(ibl)%perpst+1
        paen=liney(ibl)%paraen
        DEALLOCATE(zz(ibl)%arr)
        ALLOCATE(zz(ibl)%arr(npe,nq,0:paen))
      ENDDO
c-----------------------------------------------------------------------
c     communicate the residual to global lines for the y-direction
c     solves.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      IF (poly_deg==1) THEN
        CALL parallel_line_comm_bl2line(res,zz,nrb,ntotb,nq,liney)
      ELSE
        CALL parallel_line_comm_bl2line(ry,zz,nrb,ntotb,nq,liney)
      ENDIF
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
c-----------------------------------------------------------------------
c     perform direct solves over the y-direction of all rblocks.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        past=liney(ibl)%parast
        paen=liney(ibl)%paraen
        CALL iter_line_solve(nq,paen-past+1_i4,past+1_i4,
     $    fac%bl_fac(ibl)%adiy,zz(ibl)%arr(:,:,past:))
        IF (past>0) zz(ibl)%arr(:,:,0)=zz(ibl)%arr(:,:,paen)
      ENDDO
c-----------------------------------------------------------------------
c     communicate back to block organization, using adr for temporary
c     storage.
c-----------------------------------------------------------------------
      CALL timer(timeline_st)
      IF (poly_deg==1) THEN
        CALL parallel_line_comm_line2bl(adr,zz,nrb,ntotb,nq,liney)
      ELSE
        CALL parallel_line_comm_line2bl(ry,zz,nrb,ntotb,nq,liney)
      ENDIF
      CALL timer(timeline_en)
      time_line = time_line + timeline_en-timeline_st
c-----------------------------------------------------------------------
c     compute 1D cell interior data.
c-----------------------------------------------------------------------
      IF (poly_deg>1) THEN
        DO ibl=1,nrb
          mxb=fac%bl_fac(ibl)%mx
          myb=fac%bl_fac(ibl)%my
          DO iy=1,myb
            DO ix=0,mxb*perpb_max
              DO jq=1,npdm
                ryv(ibl)%arr(jq,ix,iy)=ryv(ibl)%arr(jq,ix,iy)-
     $              SUM(fac%bl_fac(ibl)%ylc(  :nq,jq,ix,iy)*
     $                  ry(ibl)%arr(:,ix,iy-1)
     $                 +fac%bl_fac(ibl)%ylc(nq+1:,jq,ix,iy)*
     $                  ry(ibl)%arr(:,ix,iy  ))
              ENDDO
            ENDDO
          ENDDO
          adr(ibl)%arr=ry(ibl)%arr(:,0:mxb,:)
          DO ibase=1,poly_deg-1
            xad=ibase*mxb
            adr(ibl)%arrh(:,ibase,:,:)=ry(ibl)%arr(:,xad+1:xad+mxb,:)
          ENDDO
          adr(ibl)%arrv=RESHAPE(ryv(ibl)%arr(:,0:mxb,:),
     $                          (/nq,poly_deg-1_i4,mxb+1_i4,myb/))
          IF (use_int) THEN
            DO iy=1,myb
              DO ix=1,mxb
                DO ibase=2*poly_deg,poly_d2
                  ibarr=ibase-2*poly_deg+1
                  ib0=nq*(fac%bl_fac(ibl)%iym(ibase)-2)+1
                  ib1=nq*(fac%bl_fac(ibl)%iym(ibase)-1)
                  xad=(fac%bl_fac(ibl)%ixm(ibase)-1)*mxb
                  adr(ibl)%arri(:,ibarr,ix,iy)=
     $              ryv(ibl)%arr(ib0:ib1,xad+ix,iy)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          DEALLOCATE(ry(ibl)%arr,ryv(ibl)%arr)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     deallocate line storage and sum the results of the two directions.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        DEALLOCATE(zz(ibl)%arr)
        zee(ibl)%arr=0.5*(zee(ibl)%arr+adr(ibl)%arr)
        IF (poly_deg>1) THEN
          zee(ibl)%arrh=0.5*(zee(ibl)%arrh+adr(ibl)%arrh)
          zee(ibl)%arrv=0.5*(zee(ibl)%arrv+adr(ibl)%arrv)
          IF (use_int) zee(ibl)%arri=0.5*(zee(ibl)%arri+adr(ibl)%arri)
        ENDIF
        myb=fac%bl_fac(ibl)%my*poly_deg
        IF (fac%bl_fac(ibl)%degenerate) THEN
          DO iq=1,nq
            IF (poly_deg==1) THEN
              zee(ibl)%arr(iq,0,:)=SUM(zee(ibl)%arr(iq,0,:))/(myb+1)
            ELSE
              zee(ibl)%arr(iq,0,:)=(SUM(zee(ibl)%arr(iq,0,:))
     $                       +SUM(zee(ibl)%arrv(iq,:,0,:)))/(myb+1)
              zee(ibl)%arrv(iq,:,0,:)=zee(ibl)%arr(iq,0,0)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_gl_ld_solve
c-----------------------------------------------------------------------
c     subprogram 10. iter_ld_elim
c     find the factors needed to eliminate connections begtween grid
c     vertices and cell centered data (1D sense) for the line Jacobi
c     preconditioners.  this should only be called if poly_deg>1.
c-----------------------------------------------------------------------
      SUBROUTINE iter_ld_elim(mat,fac,nq,precon,poly_deg,poly_d2,nrb)

      TYPE(global_matrix_type), INTENT(IN) :: mat
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac
      INTEGER(i4), INTENT(IN) :: nq,poly_deg,poly_d2,nrb
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: ibl,mxb,myb,iq,jq,ix,iy,iyst,jx,jy,
     $               iof,jof,xad,yad,xst,yst,npdm,
     $               line_elmax,perpb_max,itype,jtype,nqo,iqo,jqo
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: xx
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: cc,id
      LOGICAL :: sing
c-----------------------------------------------------------------------
c     if cell-interior elements have already been eliminated from the
c     matrix, compute 1D eliminations along grid lines only.
c-----------------------------------------------------------------------
      IF (mat%eliminated) THEN
        line_elmax=1
        perpb_max=2*poly_deg-1
      ELSE
        line_elmax=poly_deg
        perpb_max=poly_d2
      ENDIF
c-----------------------------------------------------------------------
c     loop over blocks and set vertex limits.
c-----------------------------------------------------------------------
      npdm=nq*(poly_deg-1)
      block_loop: DO ibl=1,nrb
        mxb=fac%bl_fac(ibl)%mx
        myb=fac%bl_fac(ibl)%my
        IF (fac%bl_fac(ibl)%perblock.AND.precon(1:1)=='b') THEN
          iyst=1
        ELSE
          iyst=0
        ENDIF
c-----------------------------------------------------------------------
c       for the x-direction, invert matrix connections within grid
c       lines.
c-----------------------------------------------------------------------
        IF (precon(8:8)=='x'.OR.precon(8:8)=='a') THEN
          ALLOCATE(cc(npdm,npdm,iyst:myb*line_elmax))
          ALLOCATE(id(npdm,iyst:myb*line_elmax,npdm))
          ALLOCATE(xx(npdm,iyst:myb*line_elmax))
          id=0
          DO iq=1,npdm
            id(iq,:,iq)=1
          ENDDO
          DO ix=1,mxb
            DO jy=1,line_elmax
              itype=MIN(2_i4*jy,4_i4)
              yst=mat%rbl_mat(ibl)%iy0(itype)
              yad=(jy-1)*myb
              iof=MAX(0_i4,npdm*(jy-2))
              DO iy=MAX(iyst,yst),myb
                cc(:,:,iy+yad)=
     $            mat%rbl_mat(ibl)%mat(itype,itype)%
     $                arr(iof+1:iof+npdm,0,0,iof+1:iof+npdm,ix,iy)
              ENDDO
            ENDDO
            CALL math_solve_q1_sym(npdm,SIZE(cc,3),cc,xx,id(:,:,1),
     $                             'factor',sing)
            IF (sing) CALL nim_stop
     $        ('Iter_ld_elim: high-order diag does not factor.')
            DO iq=1,npdm
              CALL math_solve_q1_sym(npdm,SIZE(cc,3),cc,xx,id(:,:,iq),
     $                               'solve',sing)
              fac%bl_fac(ibl)%xelim(iq,:,ix,:)=xx
            ENDDO
          ENDDO
          DEALLOCATE(cc,id,xx)
c-----------------------------------------------------------------------
c         collect the line-to-cell and cell-to-line connections.
c         the former are stored after applying the cell-to-cell
c         inverse.
c-----------------------------------------------------------------------
          ALLOCATE(cc(2*nq,npdm,iyst:myb*line_elmax))
          DO ix=1,mxb
            DO jy=1,line_elmax
              itype=MIN(2_i4*jy,4_i4)
              yst=mat%rbl_mat(ibl)%iy0(itype)
              yad=(jy-1)*myb
              iof=MAX(0_i4,npdm*(jy-2))
              jof=MAX(0_i4,nq*(jy-2))
              jtype=itype-1
              IF (fac%bl_fac(ibl)%degenerate.AND.ix==1) THEN
                DO iy=MAX(iyst,yst),myb
                  cc(1:nq,:,   iy+yad)=0
                  fac%bl_fac(ibl)%xcl(:,1:nq,ix,iy+yad)=0
                ENDDO
              ELSE
                DO iy=MAX(iyst,yst),myb
                  cc(1:nq,:,   iy+yad)=
     $              mat%rbl_mat(ibl)%mat(jtype,itype)%
     $                arr(jof+1:jof+nq,-1,0,iof+1:iof+npdm,ix,iy)
                  fac%bl_fac(ibl)%xcl(:,1:nq,ix,iy+yad)=
     $              mat%rbl_mat(ibl)%mat(itype,jtype)%
     $                arr(iof+1:iof+npdm,1,0,jof+1:jof+nq,ix-1,iy)
                ENDDO
              ENDIF
              DO iy=MAX(iyst,yst),myb
                cc(nq+1:2*nq,:,   iy+yad)=
     $            mat%rbl_mat(ibl)%mat(jtype,itype)%
     $              arr(jof+1:jof+nq,0,0,iof+1:iof+npdm,ix,iy)
                fac%bl_fac(ibl)%xcl(:,nq+1:2*nq,ix,iy+yad)=
     $            mat%rbl_mat(ibl)%mat(itype,jtype)%
     $              arr(iof+1:iof+npdm,0,0,jof+1:jof+nq,ix,iy)
              ENDDO
            ENDDO
            DO iy=iyst,myb*line_elmax
              DO iq=1,npdm
                DO jq=1,2*nq
                  fac%bl_fac(ibl)%xlc(jq,iq,ix,iy)=
     $              SUM(fac%bl_fac(ibl)%xelim(:,iq,ix,iy)*cc(jq,:,iy))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(cc)
        ENDIF
c-----------------------------------------------------------------------
c       for the y-direction, invert matrix connections within grid
c       lines.
c-----------------------------------------------------------------------
        IF (precon(8:8)=='y'.OR.precon(8:8)=='a') THEN
          ALLOCATE(cc(npdm,npdm,0:mxb*line_elmax))
          ALLOCATE(id(npdm,0:mxb*line_elmax,npdm))
          ALLOCATE(xx(npdm,0:mxb*line_elmax))
          id=0
          DO iq=1,npdm
            id(iq,:,iq)=1
          ENDDO
          DO iy=1,myb
            DO jx=1,line_elmax
              itype=MIN(3_i4*jx,4_i4)
              xst=mat%rbl_mat(ibl)%ix0(itype)
              xad=(jx-1)*mxb
              IF (jx==1) THEN
                nqo=nq
              ELSE
                nqo=npdm
              ENDIF
              DO ix=xst,mxb
                DO iof=0,poly_deg-2
                  iqo=iof*nqo+MAX(0_i4,jx-2_i4)*nq
                  DO jof=0,poly_deg-2
                    jqo=jof*nqo+MAX(0_i4,jx-2_i4)*nq
                    cc(jof*nq+1:(jof+1)*nq,iof*nq+1:(iof+1)*nq,ix+xad)=
     $                mat%rbl_mat(ibl)%mat(itype,itype)%
     $                  arr(jqo+1:jqo+nq,0,0,iqo+1:iqo+nq,ix,iy)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
            CALL math_solve_q1_sym(npdm,SIZE(cc,3),cc,xx,id(:,:,1),
     $                             'factor',sing)
            IF (sing) CALL nim_stop
     $        ('Iter_ld_elim: high-order diag does not factor.')
            DO iq=1,npdm
              CALL math_solve_q1_sym(npdm,SIZE(cc,3),cc,xx,id(:,:,iq),
     $                               'solve',sing)
              fac%bl_fac(ibl)%yelim(iq,:,:,iy)=xx
            ENDDO
          ENDDO
          DEALLOCATE(cc,id,xx)
c-----------------------------------------------------------------------
c         collect the line-to-cell and cell-to-line connections.
c         the former are stored after applying the cell-to-cell
c         inverse.
c-----------------------------------------------------------------------
          ALLOCATE(cc(2*nq,npdm,0:mxb*line_elmax))
          DO iy=1,myb
            DO jx=1,line_elmax
              itype=MIN(3_i4*jx,4_i4)
              xst=mat%rbl_mat(ibl)%ix0(itype)
              xad=(jx-1)*mxb
              jtype=itype-2
              IF (jx==1) THEN
                nqo=nq
                jqo=0
              ELSE
                nqo=npdm
                jqo=nq*(jx-2)
              ENDIF
              DO ix=xst,mxb
                DO iof=0,poly_deg-2
                  iqo=iof*nqo+MAX(0_i4,jx-2_i4)*nq
                  cc(1:nq,iof*nq+1:(iof+1)*nq,ix+xad)=
     $              mat%rbl_mat(ibl)%mat(jtype,itype)%
     $                arr(jqo+1:jqo+nq,0,-1,iqo+1:iqo+nq,ix,iy)
                  fac%bl_fac(ibl)%ycl(iof*nq+1:(iof+1)*nq,1:nq,
     $                                ix+xad,iy)=
     $              mat%rbl_mat(ibl)%mat(itype,jtype)%
     $                arr(iqo+1:iqo+nq,0, 1,jqo+1:jqo+nq,ix,iy-1)
                  cc(nq+1:2*nq,iof*nq+1:(iof+1)*nq,ix+xad)=
     $              mat%rbl_mat(ibl)%mat(jtype,itype)%
     $                arr(jqo+1:jqo+nq,0,0,iqo+1:iqo+nq,ix,iy)
                  fac%bl_fac(ibl)%ycl(iof*nq+1:(iof+1)*nq,nq+1:2*nq,
     $                                ix+xad,iy)=
     $              mat%rbl_mat(ibl)%mat(itype,jtype)%
     $                arr(iqo+1:iqo+nq,0,0,jqo+1:jqo+nq,ix,iy)
                ENDDO
              ENDDO
            ENDDO
            DO ix=0,mxb*line_elmax
              DO iq=1,npdm
                DO jq=1,2*nq
                  fac%bl_fac(ibl)%ylc(jq,iq,ix,iy)=
     $              SUM(fac%bl_fac(ibl)%yelim(:,iq,ix,iy)*cc(jq,:,ix))
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(cc)
        ENDIF
      ENDDO block_loop
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_ld_elim
c-----------------------------------------------------------------------
c     subprogram 11. iter_inner_pr
c     find the inner product of two arrays.  this is used for two
c     nonconforming arrays (different dimensions but same number of
c     components) to avoid slow reshape calls.
c-----------------------------------------------------------------------
      FUNCTION iter_inner_pr(arr1,arr2,n) RESULT(pr)

      REAL(r8) :: pr

      INTEGER(i4), INTENT(IN) :: n
      REAL(r8), DIMENSION(*), INTENT(IN) :: arr1,arr2

      pr=SUM(arr1(1:n)*arr2(1:n))

      RETURN
      END FUNCTION iter_inner_pr
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_real_line


c-----------------------------------------------------------------------
c     module containing routines for setting-up preconditioners.
c-----------------------------------------------------------------------
      MODULE iter_real_fac
      USE local
      USE edge_type_mod
      USE matrix_type_mod
      USE factor_type_mod
      USE vector_type_mod
      USE matrix_mod
      USE edge
      USE seam_storage_mod
      USE iter_real_direct
      USE iter_real_line
      IMPLICIT NONE

      INTEGER(i4), DIMENSION(:,:), PRIVATE, ALLOCATABLE :: offx,offy
      REAL(r8), DIMENSION(:,:,:,:), PRIVATE, ALLOCATABLE ::
     $             dg11t,dg11f,dg12,dg21,dg13,dg31,dg23,dg32
      LOGICAL, DIMENSION(:), PRIVATE, ALLOCATABLE :: degflags

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 12. iter_fac_alloc_real.
c     allocates storage for preconditioner matrix factors.  these
c     operations were previously part of iter_init.  note that the
c     elimination flag is passed in, since the flag in mat_str is
c     not set at this point.
c-----------------------------------------------------------------------
      SUBROUTINE iter_fac_alloc_real(mat_str,fac_str,nqty,precon,
     $                               int_elim)
      USE pardata
      USE iter_dir_fac, ONLY: iter_dirfac_alloc

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(matrix_factor_type), INTENT(OUT) :: fac_str
      INTEGER(i4), INTENT(IN) :: nqty
      CHARACTER(*), INTENT(IN) :: precon
      LOGICAL, INTENT(IN) :: int_elim

      INTEGER(i4) :: ibl,mx,my,n,nd,mv,il,iyst,past,paen,npe,
     $  itype_dst,itype_den,ibase,ix0,iy0,pdm1,line_elmax,itype,nqd,
     $  nrb,ntotb,poly_deg,poly_d2,nbtype,ientry
      LOGICAL, EXTERNAL :: per_block,deg_block,direct_check
c-----------------------------------------------------------------------
c     verify that preconitioner choice is available.
c-----------------------------------------------------------------------
      IF (.NOT.(precon=='lapack'      .OR. precon=='diagonal'
     $     .OR. precon(1:7)=='bl_diag'
     $     .OR. precon=='gl_diaga'    .OR. precon=='seq_slu'
     $     .OR. precon(1:5)=='slu_d'  .OR. precon=='no prec'))
     $  CALL nim_stop
     $    ('iter_fac_alloc_real: preconditioner choice is not
     $      recognized.')
c-----------------------------------------------------------------------
c     determine block information.
c-----------------------------------------------------------------------
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      fac_str%nqty=nqty
c-----------------------------------------------------------------------
c     nullify pointers.
c-----------------------------------------------------------------------
      fac_str%spp%j_acc=>NULL()
      fac_str%spp%start_acc=>NULL()
      fac_str%spp%start_loc=>NULL()
      fac_str%bl_fac=>NULL()
      fac_str%bl_spp=>NULL()
      fac_str%mem_id=>NULL()
      IF (precon=='no prec') RETURN
c-----------------------------------------------------------------------
c     if using a direct solver, call the dir_alloc routine and return.
c-----------------------------------------------------------------------
      IF (direct_check(precon)) THEN
        CALL factor_ptr_copy_mat_info(fac_str,mat_str%rbl_mat)
        ALLOCATE(fac_str%bl_spp(ntotb))
        DO ibl=1,ntotb
          fac_str%bl_spp(ibl)%row_ind=>NULL()
          fac_str%bl_spp(ibl)%mem_id=>NULL()
          fac_str%bl_spp(ibl)%degenerate=deg_block(ibl)
          fac_str%bl_spp(ibl)%perblock=
     $      per_block(ibl,fac_str%mat_info%mx)
        ENDDO
        CALL iter_dirfac_alloc(fac_str%mat_info,fac_str%spp,
     $                         fac_str%bl_spp,int_elim,nrb,nqty,precon)
c-----------------------------------------------------------------------
c       when using lapack, allocate the space for the factors.
c-----------------------------------------------------------------------
        IF (precon=="lapack")
     $    ALLOCATE(fac_str%a11(fac_str%spp%nbw+1,fac_str%spp%nrow))
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     determine which rblocks are periodic.
c-----------------------------------------------------------------------
      ALLOCATE(fac_str%bl_fac(ntotb))
      DO ibl=1,nrb
        mx=SIZE(mat_str%rbl_mat(ibl)%mat(1,1)%arr,5)-1
        fac_str%bl_fac(ibl)%perblock=per_block(ibl,mx)
      ENDDO
      DO ibl=1,ntotb
        fac_str%bl_fac(ibl)%degenerate=deg_block(ibl)
c-----------------------------------------------------------------------
c       nullify bl_fac pointers.
c-----------------------------------------------------------------------
        fac_str%bl_fac(ibl)%pmat=>NULL()
        fac_str%bl_fac(ibl)%mem_id=>NULL()
      ENDDO
c-----------------------------------------------------------------------
c     determine the number of different basis types.
c-----------------------------------------------------------------------
      IF (nrb>0) THEN
        nbtype=mat_str%rbl_mat(1)%nbtype
        IF (nbtype>1) THEN
          poly_deg=SUM(mat_str%rbl_mat(1)%nb_type(1:2))
        ELSE
          poly_deg=1
        ENDIF
        poly_d2=poly_deg**2
      ELSE
c-PRE
        nbtype=1
        poly_deg=1
        poly_d2=1
      ENDIF
      pdm1=poly_deg-1
c-----------------------------------------------------------------------
c     allocate the arrays needed for diagonal preconditioning.
c-PRE this will need modification for triangles, too.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        IF (precon=='diagonal'.OR.ibl>nrb) THEN
          itype_dst=1
        ELSE IF (precon(2:7)=='l_diag'.OR.direct_check(precon)) THEN
          CYCLE
        ELSE
          itype_dst=2
        ENDIF
        IF (int_elim) THEN
          itype_den=MIN(nbtype,3_i4)
        ELSE
          itype_den=MIN(nbtype,4_i4)
        ENDIF
        ALLOCATE(fac_str%bl_fac(ibl)%pmat(itype_dst:nbtype))
        ALLOCATE(fac_str%bl_fac(ibl)%ix0(nbtype))
        ALLOCATE(fac_str%bl_fac(ibl)%iy0(nbtype))
        fac_str%bl_fac(ibl)%ix0=0
        fac_str%bl_fac(ibl)%iy0=0
        IF (ibl<=nrb) THEN
          mx=SIZE(mat_str%rbl_mat(ibl)%mat(1,1)%arr,5)-1
          my=SIZE(mat_str%rbl_mat(ibl)%mat(1,1)%arr,6)-1
          DO itype=itype_dst,nbtype
            ix0=mat_str%rbl_mat(ibl)%ix0(itype)
            iy0=mat_str%rbl_mat(ibl)%iy0(itype)
            fac_str%bl_fac(ibl)%ix0(itype)=ix0
            fac_str%bl_fac(ibl)%iy0(itype)=iy0
            nqd=mat_str%rbl_mat(ibl)%nq_type(itype)
            IF (itype<=itype_den) THEN
              ALLOCATE(fac_str%bl_fac(ibl)%pmat(itype)%
     $                 arr(nqd,nqd,ix0:mx,iy0:my))
            ENDIF
          ENDDO
        ELSE
          mx=SIZE(mat_str%tbl_mat(ibl)%lmat)-1
          IF (itype_dst==1) THEN
            ALLOCATE(fac_str%bl_fac(ibl)%pmat(1)%
     $               arr(nqty,nqty,0:mx,0:0))
          ENDIF
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     allocate storage for other preconditioning options.
c-----------------------------------------------------------------------
      IF (int_elim) THEN
        line_elmax=1
      ELSE
        line_elmax=poly_deg
      ENDIF
      DO ibl=1,ntotb
        IF (ibl<=nrb) THEN
          mx=SIZE(mat_str%rbl_mat(ibl)%mat(1,1)%arr,5)-1
          my=SIZE(mat_str%rbl_mat(ibl)%mat(1,1)%arr,6)-1
        ELSE
          mx=SIZE(mat_str%tbl_mat(ibl)%lmat)-1
          my=0
        ENDIF
        fac_str%bl_fac(ibl)%mx=mx
        fac_str%bl_fac(ibl)%my=my
        IF (precon=='bl_diaga'.AND.ibl<=nrb) THEN
          iyst=0
          IF (fac_str%bl_fac(ibl)%perblock) iyst=1
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adix(iyst:my*poly_deg,nqty,nqty,-1:1,0:mx))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xelim(nqty*pdm1,nqty*pdm1,mx,iyst:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xlc(2*nqty,nqty*pdm1,mx,iyst:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xcl(nqty*pdm1,2*nqty,mx,iyst:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adiy(0:mx*poly_deg,nqty,nqty,-(1+iyst):(1+iyst),
     $                  iyst:my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             yelim(nqty*pdm1,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ylc(2*nqty,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ycl(nqty*pdm1,2*nqty,0:mx*line_elmax,my))
        ELSE IF (precon=='bl_diagx'.AND.ibl<=nrb) THEN
          iyst=0
          IF (fac_str%bl_fac(ibl)%perblock) iyst=1
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adix(iyst:my*poly_deg,nqty,nqty,-1:1,0:mx))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xelim(nqty*pdm1,nqty*pdm1,mx,iyst:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xlc(2*nqty,nqty*pdm1,mx,iyst:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xcl(nqty*pdm1,2*nqty,mx,iyst:my*line_elmax))
        ELSE IF (precon=='bl_diagy'.AND.ibl<=nrb) THEN
          iyst=0
          IF (fac_str%bl_fac(ibl)%perblock) iyst=1
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adiy(0:mx*poly_deg,nqty,nqty,-(1+iyst):(1+iyst),
     $                  iyst:my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             yelim(nqty*pdm1,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ylc(2*nqty,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ycl(nqty*pdm1,2*nqty,0:mx*line_elmax,my))
        ELSE IF (precon=='gl_diaga'.AND.ibl<=nrb) THEN
          npe=linex(ibl)%perpen-linex(ibl)%perpst+1
          past=linex(ibl)%parast
          paen=linex(ibl)%paraen
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adix(npe,nqty,nqty,-(1+past):(1+past),past:paen))
          npe=liney(ibl)%perpen-liney(ibl)%perpst+1
          past=liney(ibl)%parast
          paen=liney(ibl)%paraen
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             adiy(npe,nqty,nqty,-(1+past):(1+past),past:paen))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xelim(nqty*pdm1,nqty*pdm1,mx,0:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xlc(2*nqty,nqty*pdm1,mx,0:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             xcl(nqty*pdm1,2*nqty,mx,0:my*line_elmax))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             yelim(nqty*pdm1,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ylc(2*nqty,nqty*pdm1,0:mx*line_elmax,my))
          ALLOCATE(fac_str%bl_fac(ibl)%
     $             ycl(nqty*pdm1,2*nqty,0:mx*line_elmax,my))
        ENDIF
c-----------------------------------------------------------------------
c       also create the ixm and iym index arrays that help off-diagonal
c       references in the line-Jacobi routines.
c-----------------------------------------------------------------------
        IF (precon(4:7)=='diag'.AND.ibl<=nrb) THEN
          ALLOCATE(fac_str%bl_fac(ibl)%ixm(poly_d2))
          ALLOCATE(fac_str%bl_fac(ibl)%iym(poly_d2))
          DO il=1,poly_deg
            fac_str%bl_fac(ibl)%ixm(il)=il
            fac_str%bl_fac(ibl)%ixm(il-1+poly_deg)=1
            fac_str%bl_fac(ibl)%iym(il)=1
            fac_str%bl_fac(ibl)%iym(il-1+poly_deg)=il
          ENDDO
          DO il=2*poly_deg,poly_d2
            fac_str%bl_fac(ibl)%ixm(il)=
     $        MODULO((il-2_i4*poly_deg),poly_deg-1_i4)+2
            fac_str%bl_fac(ibl)%iym(il)=(il-2*poly_deg)/(poly_deg-1)+2
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_fac_alloc_real
c-----------------------------------------------------------------------
c     subprogram 13. iter_fac_dealloc_real.
c     deallocates storage for preconditioner matrix factors.
c-----------------------------------------------------------------------
      SUBROUTINE iter_fac_dealloc_real(fac_str,nrb,precon)
      USE pardata
      USE iter_dir_fac, ONLY: iter_dirfac_dealloc

      TYPE(matrix_factor_type), INTENT(INOUT) :: fac_str
      CHARACTER(*), INTENT(IN) :: precon
      INTEGER(i4), INTENT(IN) :: nrb

      INTEGER(i4) :: ibl,ibasis_dst,ibase,ntotb,ierr
      INTEGER(i4), PARAMETER :: iopt=4,nrhs=1
      REAL(r8) :: accdum,bdum
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     get block sizes.
c-----------------------------------------------------------------------
      ntotb=SIZE(fac_str%bl_fac)
c-----------------------------------------------------------------------
c     deallocate storage for the preconditioners over the linear and
c     bilinear basis types.
c-----------------------------------------------------------------------
      IF (direct_check(precon)) THEN
        IF (precon=='lapack') THEN
          DEALLOCATE(fac_str%a11)
        ELSE IF (precon=='seq_slu') THEN
          CALL c_fortran_dgssv(
     $      iopt,fac_str%spp%nrow,fac_str%spp%nnz,nrhs,accdum,
     $      fac_str%spp%j_acc,fac_str%spp%start_acc,bdum,
     $      fac_str%spp%nrow,fac_str%spp%acc_handle,ierr)
        ELSE IF (precon(1:7)=='slu_dst') THEN
          CALL c_fortran_pdloc(
     $      iopt,fac_str%spp%nrow,fac_str%spp%mloc,fac_str%spp%nnzloc,
     $      fac_str%spp%fstrow,nrhs,accdum,fac_str%spp%j_acc,
     $      fac_str%spp%start_loc,bdum,INT(fac_str%spp%mloc,i4),
     $      fac_str%spp%acc_handle,slugrid_handle,
     $      fac_str%spp%iopts,fac_str%spp%dopts,ierr)
        ELSE IF (precon=='slu_dist') THEN
          CALL c_fortran_pdglob(iopt,fac_str%spp%nrow,fac_str%spp%nnz,
     $      nrhs,accdum,fac_str%spp%j_acc,fac_str%spp%start_acc,bdum,
     $      INT(fac_str%spp%nrow,i4),fac_str%spp%acc_handle,
     $      slugrid_handle,ierr)
        ENDIF
        CALL iter_dirfac_dealloc(fac_str%spp)
        CALL factor_dealloc_mat_info(fac_str)
        DO ibl=1,nrb
          DO ibase=1,SIZE(fac_str%bl_spp(ibl)%row_ind)
            DEALLOCATE(fac_str%bl_spp(ibl)%row_ind(ibase)%rarr)
          ENDDO
          DEALLOCATE(fac_str%bl_spp(ibl)%row_ind)
        ENDDO
        DEALLOCATE(fac_str%bl_spp)
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     deallocate the arrays needed for diagonal preconditioning.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        IF (precon=='diagonal'.OR.ibl>nrb) THEN
          ibasis_dst=1
        ELSE IF (precon(2:7)=='l_diag'.OR.direct_check(precon)) THEN
          CYCLE
        ELSE
          ibasis_dst=2
        ENDIF
        DEALLOCATE(fac_str%bl_fac(ibl)%ix0)
        DEALLOCATE(fac_str%bl_fac(ibl)%iy0)
        DO ibase=ibasis_dst,SIZE(fac_str%bl_fac(ibl)%pmat)+ibasis_dst-1
          IF (ALLOCATED(fac_str%bl_fac(ibl)%pmat(ibase)%arr)) THEN
            DEALLOCATE(fac_str%bl_fac(ibl)%pmat(ibase)%arr)
          ENDIF
        ENDDO
        DEALLOCATE(fac_str%bl_fac(ibl)%pmat)
      ENDDO
      DO ibl=1,ntotb
        IF (precon=='bl_diaga'.AND.ibl<=nrb) THEN
          DEALLOCATE(fac_str%bl_fac(ibl)%adix)
          DEALLOCATE(fac_str%bl_fac(ibl)%xelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%xlc)
          DEALLOCATE(fac_str%bl_fac(ibl)%xcl)
          DEALLOCATE(fac_str%bl_fac(ibl)%adiy)
          DEALLOCATE(fac_str%bl_fac(ibl)%yelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%ylc)
          DEALLOCATE(fac_str%bl_fac(ibl)%ycl)
        ELSE IF (precon=='bl_diagx'.AND.ibl<=nrb) THEN
          DEALLOCATE(fac_str%bl_fac(ibl)%adix)
          DEALLOCATE(fac_str%bl_fac(ibl)%xelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%xlc)
          DEALLOCATE(fac_str%bl_fac(ibl)%xcl)
        ELSE IF (precon=='bl_diagy'.AND.ibl<=nrb) THEN
          DEALLOCATE(fac_str%bl_fac(ibl)%adiy)
          DEALLOCATE(fac_str%bl_fac(ibl)%yelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%ylc)
          DEALLOCATE(fac_str%bl_fac(ibl)%ycl)
        ELSE IF (precon=='gl_diaga'.AND.ibl<=nrb) THEN
          DEALLOCATE(fac_str%bl_fac(ibl)%adix)
          DEALLOCATE(fac_str%bl_fac(ibl)%xelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%xlc)
          DEALLOCATE(fac_str%bl_fac(ibl)%xcl)
          DEALLOCATE(fac_str%bl_fac(ibl)%adiy)
          DEALLOCATE(fac_str%bl_fac(ibl)%yelim)
          DEALLOCATE(fac_str%bl_fac(ibl)%ylc)
          DEALLOCATE(fac_str%bl_fac(ibl)%ycl)
        ENDIF
        IF (precon(4:7)=='diag'.AND.ibl<=nrb) THEN
          DEALLOCATE(fac_str%bl_fac(ibl)%ixm)
          DEALLOCATE(fac_str%bl_fac(ibl)%iym)
        ENDIF
      ENDDO
      DEALLOCATE(fac_str%bl_fac)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_fac_dealloc_real
c-----------------------------------------------------------------------
c     subprogram 14. iter_factor_real.
c     compute the approximate matrix factors within each block for
c     preconditioning.  this was formerly part of iter_init, but it
c     has been separated so that these factors can be saved and
c     reused.
c-----------------------------------------------------------------------
      SUBROUTINE iter_factor_real(mat_str,fac_str,nqty,precon,
     $                            off_diag_fac)
      USE pardata
      USE time

      TYPE(global_matrix_type), INTENT(INOUT) :: mat_str
      TYPE(matrix_factor_type), INTENT(INOUT) :: fac_str
      INTEGER(i4), INTENT(IN) :: nqty
      CHARACTER(*), INTENT(IN) :: precon
      REAL(r8), INTENT(IN) :: off_diag_fac

      INTEGER(i4) :: ibl,mx,my,iv,iq,mv,ix,iy,imat,jmat,jx,jy,
     $               itype_dst,itype_den,itype,ibase,ix0,iy0,nqd,
     $               poly_deg,poly_d2,nbtype,nrb,ntotb
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: shif
      INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: offt,offt2
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: lu,iden,tmp
      REAL(r8) :: timestart_fac,timeend_fac
      LOGICAL :: singular=.false.
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     return if there is no preconditioning.
c-----------------------------------------------------------------------
      IF (precon=='no prec') RETURN
c-----------------------------------------------------------------------
c     determine the number of basis types.
c-----------------------------------------------------------------------
      CALL timer(timestart_fac)
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      IF (nrb>0) THEN
        nbtype=mat_str%rbl_mat(1)%nbtype
        IF (nbtype>1) THEN
          poly_deg=SUM(mat_str%rbl_mat(1)%nb_type(1:2))
        ELSE
          poly_deg=1
        ENDIF
        poly_d2=poly_deg**2
      ELSE
c-PRE
        nbtype=1
        poly_deg=1
        poly_d2=1
      ENDIF
c-----------------------------------------------------------------------
c     construct arrays giving the basis numbers of neighbors.
c-----------------------------------------------------------------------
      ALLOCATE(shif(poly_deg))
      ALLOCATE(offt(poly_deg,poly_deg))
      ALLOCATE(offt2(poly_deg,poly_deg))
      ALLOCATE(offx(-poly_deg:poly_deg,poly_d2))
      ALLOCATE(offy(-poly_deg:poly_deg,poly_d2))
      DO jmat=1,poly_deg
        offt(1,jmat)=jmat-1+poly_deg
        offt(jmat,1)=jmat
      ENDDO
      DO jx=2*poly_deg,poly_d2
        imat=MODULO((jx-2_i4*poly_deg),poly_deg-1_i4)+2
        jmat=(jx-2*poly_deg)/(poly_deg-1)+2
        offt(imat,jmat)=jx
      ENDDO
      DO jx=-poly_deg,poly_deg
        shif=jx
        offt2=CSHIFT(offt,shif,DIM=1)
        offx(jx,1:poly_deg)=offt2(:,1)
        IF (poly_deg>1) THEN
          offx(jx,poly_deg+1:2*poly_deg-1)=offt2(1,2:)
          offx(jx,2*poly_deg:)=RESHAPE(offt2(2:,2:),(/(poly_deg-1)**2/))
        ENDIF
      ENDDO
      DO jy=-poly_deg,poly_deg
        shif=jy
        offt2=CSHIFT(offt,shif,DIM=2)
        offy(jy,1:poly_deg)=offt2(:,1)
        IF (poly_deg>1) THEN
          offy(jy,poly_deg+1:2*poly_deg-1)=offt2(1,2:)
          offy(jy,2*poly_deg:)=RESHAPE(offt2(2:,2:),(/(poly_deg-1)**2/))
        ENDIF
      ENDDO
      DEALLOCATE(shif,offt,offt2)
c-----------------------------------------------------------------------
c     communicate diagonal elements of the matrix.
c-----------------------------------------------------------------------
      IF (.NOT.direct_check(precon))
     $  CALL iter_mat_com(mat_str,nqty,precon,nrb,ntotb,poly_deg,
     $                    poly_d2,nbtype)
c-----------------------------------------------------------------------
c     zero out the preconditioner arrays as needed.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        IF (direct_check(precon).AND.ibl<=nrb) THEN
        ELSE IF (precon(2:8)=='l_diaga'.AND.ibl<=nrb) THEN
          fac_str%bl_fac(ibl)%adix=0
          fac_str%bl_fac(ibl)%adiy=0
        ELSE IF (precon=='bl_diagx'.AND.ibl<=nrb) THEN
          fac_str%bl_fac(ibl)%adix=0
        ELSE IF (precon=='bl_diagy'.AND.ibl<=nrb) THEN
          fac_str%bl_fac(ibl)%adiy=0
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     initialize preconditioner.  for domain decomposition, this means
c     factoring the matrix within each rblock.  for tblocks and
c     diagonal preconditioning, just factor the pointwise blocks.
c-----------------------------------------------------------------------
      IF (direct_check(precon)) THEN
        CALL fac_dir(mat_str,fac_str,fac_str%spp,nqty,nrb,precon)
        itype_dst=5
      ELSE IF (precon(1:7)=='bl_diag') THEN
        CALL iter_ld_fac(mat_str,fac_str,nqty,precon,off_diag_fac,
     $                   poly_deg,poly_d2,nrb)
        itype_dst=5
      ELSE IF (precon=='gl_diaga') THEN
        CALL iter_gl_ld_fac(mat_str,fac_str,nqty,nrb,off_diag_fac,
     $                      poly_deg,poly_d2,ntotb)
        itype_dst=5
      ELSE
        itype_dst=1
      ENDIF
      IF (mat_str%eliminated) THEN
        itype_den=MIN(nbtype,3_i4)
      ELSE
        itype_den=MIN(nbtype,4_i4)
      ENDIF
c-----------------------------------------------------------------------
c     copy all diagonal elements for diagonal preconditioning.
c-----------------------------------------------------------------------
      DO ibl=1,nrb
        DO itype=itype_dst,itype_den
          fac_str%bl_fac(ibl)%pmat(itype)%arr=
     $      mat_str%rbl_mat(ibl)%mat(itype,itype)%arr(:,0,0,:,:,:)
        ENDDO
      ENDDO
c-PRE
      DO ibl=nrb+1,ntotb
        DO iv=0,fac_str%bl_fac(ibl)%mx
          fac_str%bl_fac(ibl)%pmat(1)%arr(:,:,iv,0)
     $       =mat_str%tbl_mat(ibl)%lmat(iv)%element(1:nqty,1:nqty,0)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     restore diagonal elements of the matrix.
c-----------------------------------------------------------------------
      IF (.NOT.direct_check(precon))
     $  CALL iter_mat_rest(mat_str,nqty,precon,nrb,ntotb,
     $                     poly_deg,nbtype)
c-----------------------------------------------------------------------
c     invert diagonal block for diagonal preconditioning
c     (zee and res are used as temporary storage).
c-----------------------------------------------------------------------
        DO ibl=1,ntotb
          IF (precon=='diagonal'.OR.ibl>nrb) THEN
            itype_dst=1
          ELSE IF (precon(2:7)=='l_diag'.OR.direct_check(precon)) THEN
            CYCLE
          ELSE
            itype_dst=2
          ENDIF
          mx=fac_str%bl_fac(ibl)%mx
          my=fac_str%bl_fac(ibl)%my
          DO itype=itype_dst,itype_den
c-PRE
            IF (ibl>nrb.AND.itype>1) CYCLE
            IF (ibl<=nrb) THEN
              nqd=mat_str%rbl_mat(ibl)%nq_type(itype)
            ELSE
              nqd=nqty
            ENDIF
            ALLOCATE(iden(nqd,nqd),lu(nqd,nqd),tmp(nqd,nqd))
            iden=0
            DO iq=1,nqd
              iden(iq,iq)=1
            ENDDO
            ix0=fac_str%bl_fac(ibl)%ix0(itype)
            iy0=fac_str%bl_fac(ibl)%iy0(itype)
            DO iy=iy0,my
              DO ix=ix0,mx
                lu=fac_str%bl_fac(ibl)%pmat(itype)%arr(:,:,ix,iy)
                CALL iter_solve_sym(nqd,nqd,lu,tmp,
     $                              iden,'both',singular)
                IF (singular) CALL nim_stop
     $            ('Iter_factor unable to factor point-block.')
                fac_str%bl_fac(ibl)%pmat(itype)%arr(:,:,ix,iy)=
     $            TRANSPOSE(tmp)
              ENDDO
            ENDDO
            DEALLOCATE(iden,lu,tmp)
          ENDDO
        ENDDO

      DEALLOCATE(offx,offy)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      CALL timer(timeend_fac)
      time_fac = time_fac + timeend_fac-timestart_fac
      RETURN
      END SUBROUTINE iter_factor_real
c-----------------------------------------------------------------------
c     subprogram 15. iter_mat_com.
c     communicate the border elements of the matrix along block
c     boundaries.  the initial elements are saved in seam_save for
c     later use.
c
c     do not assume symmetry so that the preconditioners can be used
c     for non-symmetric solves.
c-----------------------------------------------------------------------
      SUBROUTINE iter_mat_com(mat,nqty,precon,nrb,ntotb,poly_deg,
     $                        poly_d2,nbtype)

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty,nrb,ntotb,poly_deg,poly_d2,nbtype
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: ibl,n,nm,iqty,jqty,iv,ix,iy,ixp,iyp,ixs,iys,mx,my,
     $               ibase,jbase,imat,jmat,ioff,joff,ityp,jtyp,iqo,jqo,
     $               idg,ndg,nqpm1,nqp1p,nqpd,ibst,iben,nb,nq2
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
      LOGICAL, EXTERNAL :: per_block,deg_block
c-----------------------------------------------------------------------
c     for degenerate blocks, sum any on-block off-diagonal matrix
c     elements that are communicated to neighboring degenerate
c     blocks.  they occur in the corner elements along the degenerate
c     side.  use dg* to save the original entries.  a 1 in the
c     third index indicates data for the lower corner, and a 2 indicates
c     data for the upper corner.
c-----------------------------------------------------------------------
      ndg=0
      ALLOCATE(degflags(nrb))
      DO ibl=1,nrb
        degflags(ibl)=deg_block(ibl)
        IF (degflags(ibl)) ndg=ndg+1
      ENDDO
      nqpm1=nqty*(poly_deg-1)
      nqp1p=nqpm1+1
      nqpd=nqty*poly_deg
      nq2=nqty*nqty
      IF (ndg>0) THEN
        ALLOCATE(dg11t(nqty,nqty,2,ndg))
        ALLOCATE(dg11f(nqty,nqty,2,ndg))
        IF (poly_deg>1) THEN
          ALLOCATE(dg12(nqty,nqpm1,2,ndg))
          ALLOCATE(dg21(nqpm1,nqty,2,ndg))
          ALLOCATE(dg13(nqty,nqpm1,2,ndg))
          ALLOCATE(dg31(nqpm1,nqty,2,ndg))
          ALLOCATE(dg23(nqpm1,nqpm1,2,ndg))
          ALLOCATE(dg32(nqpm1,nqpm1,2,ndg))
        ENDIF
      ENDIF

      idg=0
      DO ibl=1,nrb
        IF (.NOT.degflags(ibl)) CYCLE
        idg=idg+1
c-----------------------------------------------------------------------
c       copy the off-diagonal connections in corner elements and zero-
c       out that storage in rbl_mat
c-----------------------------------------------------------------------
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        my=SIZE(rmat,6)-1
        dg11t(:,:,1,idg)=rmat(1:nqty,1,-1,1:nqty,0,1)
        rmat(1:nqty,1,-1,1:nqty,0,1)=0._r8
        dg11t(:,:,2,idg)=rmat(1:nqty,1, 1,1:nqty,0,my-1)
        rmat(1:nqty,1, 1,1:nqty,0,my-1)=0._r8

        dg11f(:,:,1,idg)=rmat(1:nqty,-1, 1,1:nqty,1,0)
        rmat(1:nqty,-1, 1,1:nqty,1,0)=0._r8
        dg11f(:,:,2,idg)=rmat(1:nqty,-1,-1,1:nqty,1,my)
        rmat(1:nqty,-1,-1,1:nqty,1,my)=0._r8

        IF (poly_deg>1) THEN
          rmat=>mat%rbl_mat(ibl)%mat(2,1)%arr
          dg21(:,:,1,idg)=rmat(1:nqpm1,1,-1,1:nqty,0,1)
          rmat(1:nqpm1,1,-1,1:nqty,0,1)=0._r8
          dg21(:,:,2,idg)=rmat(1:nqpm1,1, 1,1:nqty,0,my-1)
          rmat(1:nqpm1,1, 1,1:nqty,0,my-1)=0._r8

          rmat=>mat%rbl_mat(ibl)%mat(3,1)%arr
          dg31(:,:,1,idg)=rmat(1:nqpm1,-1,1,1:nqty,1,0)
          rmat(1:nqpm1,-1,1,1:nqty,1,0)=0._r8
          dg31(:,:,2,idg)=rmat(1:nqpm1,-1,0,1:nqty,1,my)
          rmat(1:nqpm1,-1,0,1:nqty,1,my)=0._r8

          rmat=>mat%rbl_mat(ibl)%mat(1,2)%arr
          dg12(:,:,1,idg)=rmat(1:nqty,-1, 1,1:nqpm1,1,0)
          rmat(1:nqty,-1, 1,1:nqpm1,1,0)=0._r8
          dg12(:,:,2,idg)=rmat(1:nqty,-1,-1,1:nqpm1,1,my)
          rmat(1:nqty,-1,-1,1:nqpm1,1,my)=0._r8

          rmat=>mat%rbl_mat(ibl)%mat(3,2)%arr
          dg32(:,:,1,idg)=rmat(1:nqpm1,-1,1,1:nqpm1,1,0)
          rmat(1:nqpm1,-1,1,1:nqpm1,1,0)=0._r8
          dg32(:,:,2,idg)=rmat(1:nqpm1,-1,0,1:nqpm1,1,my)
          rmat(1:nqpm1,-1,0,1:nqpm1,1,my)=0._r8

          rmat=>mat%rbl_mat(ibl)%mat(1,3)%arr
          dg13(:,:,1,idg)=rmat(1:nqty,1,-1,1:nqpm1,0,1)
          rmat(1:nqty,1,-1,1:nqpm1,0,1)=0._r8
          dg13(:,:,2,idg)=rmat(1:nqty,1, 0,1:nqpm1,0,my)
          rmat(1:nqty,1, 0,1:nqpm1,0,my)=0._r8

          rmat=>mat%rbl_mat(ibl)%mat(2,3)%arr
          dg23(:,:,1,idg)=rmat(1:nqpm1,1,-1,1:nqpm1,0,1)
          rmat(1:nqpm1,1,-1,1:nqpm1,0,1)=0._r8
          dg23(:,:,2,idg)=rmat(1:nqpm1,1, 0,1:nqpm1,0,my)
          rmat(1:nqpm1,1, 0,1:nqpm1,0,my)=0._r8
        ENDIF
c-----------------------------------------------------------------------
c       add the off-diagonal corner connections to off-diagonal
c       entries along the block border.
c-----------------------------------------------------------------------
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        rmat(1:nqty,1,0,1:nqty,0,0)=
     $    rmat(1:nqty,1,0,1:nqty,0,0)+dg11t(:,:,1,idg)
        rmat(1:nqty,1,0,1:nqty,0,my)=
     $    rmat(1:nqty,1,0,1:nqty,0,my)+dg11t(:,:,2,idg)
        rmat(1:nqty,-1,0,1:nqty,1,0)=
     $    rmat(1:nqty,-1,0,1:nqty,1,0)+dg11f(:,:,1,idg)
        rmat(1:nqty,-1,0,1:nqty,1,my)=
     $    rmat(1:nqty,-1,0,1:nqty,1,my)+dg11f(:,:,2,idg)

        IF (poly_deg>1) THEN
          rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqty,1,0,1:nqty,0,0)=rmat(1:nqty,1,0,1:nqty,0,0)+
     $        dg13(:,joff+1:joff+nqty,1,idg)
            rmat(1:nqty,1,0,1:nqty,0,my)=rmat(1:nqty,1,0,1:nqty,0,my)+
     $        dg13(:,joff+1:joff+nqty,2,idg)

            rmat(1:nqty,-1,0,1:nqty,1,0)=rmat(1:nqty,-1,0,1:nqty,1,0)+
     $        dg31(joff+1:joff+nqty,:,1,idg)
            rmat(1:nqty,-1,0,1:nqty,1,my)=rmat(1:nqty,-1,0,1:nqty,1,my)+
     $        dg31(joff+1:joff+nqty,:,2,idg)
          ENDDO

          rmat=>mat%rbl_mat(ibl)%mat(2,1)%arr
          rmat(1:nqpm1,1,0,1:nqty,0,0)=
     $      rmat(1:nqpm1,1,0,1:nqty,0,0)+dg21(:,:,1,idg)
          rmat(1:nqpm1,1,0,1:nqty,0,my)=
     $      rmat(1:nqpm1,1,0,1:nqty,0,my)+dg21(:,:,2,idg)
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqpm1,1,0,1:nqty,0,0)=
     $        rmat(1:nqpm1,1,0,1:nqty,0,0)+
     $        dg23(:,joff+1:joff+nqty,1,idg)
            rmat(1:nqpm1,1,0,1:nqty,0,my)=
     $        rmat(1:nqpm1,1,0,1:nqty,0,my)+
     $        dg23(:,joff+1:joff+nqty,2,idg)
          ENDDO

          rmat=>mat%rbl_mat(ibl)%mat(1,2)%arr
          rmat(1:nqty,-1,0,1:nqpm1,1,0)=
     $      rmat(1:nqty,-1,0,1:nqpm1,1,0)+dg12(:,:,1,idg)
          rmat(1:nqty,-1,0,1:nqpm1,1,my)=
     $      rmat(1:nqty,-1,0,1:nqpm1,1,my)+dg12(:,:,2,idg)
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqty,-1,0,1:nqpm1,1,0)=
     $        rmat(1:nqty,-1,0,1:nqpm1,1,0)+
     $        dg32(joff+1:joff+nqty,:,1,idg)
            rmat(1:nqty,-1,0,1:nqpm1,1,my)=
     $        rmat(1:nqty,-1,0,1:nqpm1,1,my)+
     $        dg32(joff+1:joff+nqty,:,2,idg)
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     load seam_in from perimeter of matrix arrays, rblocks first.
c-----------------------------------------------------------------------
      rblock_in: DO ibl=1,nrb
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              seam(ibl)%vertex(iv)%seam_in(n)=rmat(jqty,0,0,iqty,ix,iy)
              seam(ibl)%vertex(iv)%seam_save(n)=
     $          seam(ibl)%vertex(iv)%seam_in(n)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       when basis functions of degree>1 are used, also load diagonal
c       elements for these bases.  seam_in ordering for these bases
c       is determined by load_dir (but not seam_save).
c-----------------------------------------------------------------------
        IF (nbtype>1) THEN
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            IF (seam(ibl)%segment(iv)%h_side) THEN
              rmat=>mat%rbl_mat(ibl)%mat(2,2)%arr
              nb=mat%rbl_mat(ibl)%nb_type(2)
            ELSE
              rmat=>mat%rbl_mat(ibl)%mat(3,3)%arr
              nb=mat%rbl_mat(ibl)%nb_type(3)
            ENDIF
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,nb,
     $                            ibst,iben)
            DO ibase=1,nb
              iqo=nqty*(ibase-1)
              nm=n*(ibase-1)
              DO iqty=1,nqty
                DO jqty=1,nqty
                  nm=nm+1
                  seam(ibl)%segment(iv)%seam_save(nm)=
     $              rmat(iqo+jqty,0,0,iqo+iqty,ix,iy)
                ENDDO
              ENDDO
            ENDDO
            nm=1
            DO ibase=ibst,iben,seam(ibl)%segment(iv)%load_dir
              seam(ibl)%segment(iv)%seam_in(nq2*(ibase-1)+1:nq2*ibase)=
     $          seam(ibl)%segment(iv)%seam_save(nm:nm+nq2-1)
              nm=nm+nq2
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       off-diagonal connections.
c-----------------------------------------------------------------------
        DO iv=1,seam(ibl)%nvert
          ix =seam(ibl)%segment(iv)%intxyn(1)
          iy =seam(ibl)%segment(iv)%intxyn(2)
          ixp=seam(ibl)%segment(iv)%intxyp(1)
          iyp=seam(ibl)%segment(iv)%intxyp(2)
          ixs=seam(ibl)%segment(iv)%intxys(1)
          iys=seam(ibl)%segment(iv)%intxys(2)
          IF (seam(ibl)%segment(iv)%h_side) THEN
            ioff=1+ixs-ix
            DO ibase=poly_deg+1,2,-1
              imat=ibase
              IF (imat>poly_deg) imat=1
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=-1,-(ibase-1),-1
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%ix0(jtyp)-1
                seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,ioff)
     $            =rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,ixs,iys)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+ixs-ixp
            DO imat=1,poly_deg
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=1,poly_deg+1-imat
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%ix0(ityp)
                seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,ioff)
     $            =rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,
     $                  ixs-joff,iys)
                ioff=ioff+2
              ENDDO
            ENDDO
          ELSE
            ioff=1+iys-iy
            DO ibase=2*poly_deg,poly_deg+1,-1
              imat=ibase
              IF (imat==2*poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=-1,-(ibase-poly_deg),-1
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%iy0(jtyp)-1
                seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,ioff)
     $            =rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,ixs,iys)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+iys-iyp
            DO ibase=poly_deg,2*poly_deg-1
              imat=ibase
              IF (imat==poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=1,2*poly_deg-ibase
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%iy0(ityp)
                seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,ioff)
     $            =rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,
     $                  ixs,iys-joff)
                ioff=ioff+2
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO rblock_in
c-----------------------------------------------------------------------
c     tblocks
c-----------------------------------------------------------------------
      tblock_in: DO ibl=nrb+1,ntotb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              seam(ibl)%vertex(iv)%seam_in(n)=
     $          mat%tbl_mat(ibl)%lmat(ix)%element(jqty,iqty,0)
            ENDDO
          ENDDO
        ENDDO
        DO iv=1,seam(ibl)%nvert
          ix =seam(ibl)%segment(iv)%intxyn(1)
          iy =seam(ibl)%segment(iv)%intxyn(2)
          ixp=seam(ibl)%segment(iv)%intxyp(1)
          iyp=seam(ibl)%segment(iv)%intxyp(2)
          seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,1)=
     $      mat%tbl_mat(ibl)%lmat(ix)%element(1:nqty,1:nqty,iy)
          seam(ibl)%segment(iv)%seam_mat_in(1:nqty,1:nqty,2)=
     $      mat%tbl_mat(ibl)%lmat(ixp)%element(1:nqty,1:nqty,iyp)
        ENDDO
      ENDDO tblock_in
c-----------------------------------------------------------------------
c     copy seam_mat_in information into the seam_save arrays for
c     restoration later.  this is necessary since at least one of the
c     global precondition options requires seaming during factorization.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%segment(iv)%seam_mat_save
     $      (1:nqty,1:nqty,1:poly_d2+poly_deg)=
     $        seam(ibl)%segment(iv)%seam_mat_in
     $          (1:nqty,1:nqty,1:poly_d2+poly_deg)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     seam diagonal and off-diagonal block-border matrix elements.
c-----------------------------------------------------------------------
      CALL edge_network(n,0_i4,poly_deg-1_i4,.false.)
      CALL edge_seg_network(nqty,0_i4,poly_d2+poly_deg,.false.)
c-----------------------------------------------------------------------
c     unload seam_out into matrix arrays.
c-----------------------------------------------------------------------
      rblock_out: DO ibl=1,nrb
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              rmat(jqty,0,0,iqty,ix,iy)=
     $             seam(ibl)%vertex(iv)%seam_out(n)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       unload rblock higher-order diagonals.
c-----------------------------------------------------------------------
        IF (nbtype>1) THEN
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            IF (seam(ibl)%segment(iv)%h_side) THEN
              rmat=>mat%rbl_mat(ibl)%mat(2,2)%arr
              nb=mat%rbl_mat(ibl)%nb_type(2)
            ELSE
              rmat=>mat%rbl_mat(ibl)%mat(3,3)%arr
              nb=mat%rbl_mat(ibl)%nb_type(3)
            ENDIF
            CALL edge_load_limits(seam(ibl)%segment(iv)%load_dir,nb,
     $                            ibst,iben)
            nm=0
            DO ibase=ibst,iben,seam(ibl)%segment(iv)%load_dir
              iqo=nqty*(ibase-1)
              DO iqty=1,nqty
                DO jqty=1,nqty
                  nm=nm+1
                  rmat(iqo+jqty,0,0,iqo+iqty,ix,iy)=
     $              seam(ibl)%segment(iv)%seam_out(nm)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       unload off-diagonal elements.
c-----------------------------------------------------------------------
        DO iv=1,seam(ibl)%nvert
          ix =seam(ibl)%segment(iv)%intxyn(1)
          iy =seam(ibl)%segment(iv)%intxyn(2)
          ixp=seam(ibl)%segment(iv)%intxyp(1)
          iyp=seam(ibl)%segment(iv)%intxyp(2)
          ixs=seam(ibl)%segment(iv)%intxys(1)
          iys=seam(ibl)%segment(iv)%intxys(2)
          IF (seam(ibl)%segment(iv)%h_side) THEN
            ioff=1+ixs-ix
            DO ibase=poly_deg+1,2,-1
              imat=ibase
              IF (imat>poly_deg) imat=1
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=-1,-(ibase-1),-1
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%ix0(jtyp)-1
                rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,ixs,iys)=
     $           seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+ixs-ixp
            DO imat=1,poly_deg
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=1,poly_deg+1-imat
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%ix0(ityp)
                rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,ixs-joff,iys)=
     $           seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
          ELSE
            ioff=1+iys-iy
            DO ibase=2*poly_deg,poly_deg+1,-1
              imat=ibase
              IF (imat==2*poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=-1,-(ibase-poly_deg),-1
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%iy0(jtyp)-1
                rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,ixs,iys)=
     $           seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+iys-iyp
            DO ibase=poly_deg,2*poly_deg-1
              imat=ibase
              IF (imat==poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=1,2*poly_deg-ibase
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%iy0(ityp)
                rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,ixs,iys-joff)=
     $           seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO rblock_out
c-----------------------------------------------------------------------
c     tblocks
c-----------------------------------------------------------------------
      tblock_out: DO ibl=nrb+1,ntotb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              mat%tbl_mat(ibl)%lmat(ix)%element(jqty,iqty,0)=
     $             seam(ibl)%vertex(iv)%seam_out(n)
            ENDDO
          ENDDO
        ENDDO
        IF (precon/='diagonal') THEN
          DO iv=1,seam(ibl)%nvert
            ix =seam(ibl)%segment(iv)%intxyn(1)
            iy =seam(ibl)%segment(iv)%intxyn(2)
            ixp=seam(ibl)%segment(iv)%intxyp(1)
            iyp=seam(ibl)%segment(iv)%intxyp(2)
            mat%tbl_mat(ibl)%lmat(ix)%element(1:nqty,1:nqty,iy)=
     $         seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,1)
            mat%tbl_mat(ibl)%lmat(ixp)%element(1:nqty,1:nqty,iyp)=
     $         seam(ibl)%segment(iv)%seam_mat_out(1:nqty,1:nqty,2)
          ENDDO
        ENDIF
      ENDDO tblock_out
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_mat_com
c-----------------------------------------------------------------------
c     subprogram 16. iter_mat_rest.
c     restore the border elements of the matrix from seam_save.
c-----------------------------------------------------------------------
      SUBROUTINE iter_mat_rest(mat,nqty,precon,nrb,ntotb,poly_deg,
     $                         nbtype)

      TYPE(global_matrix_type), INTENT(INOUT) :: mat
      INTEGER(i4), INTENT(IN) :: nqty,nrb,ntotb,poly_deg,nbtype
      CHARACTER(*), INTENT(IN) :: precon

      INTEGER(i4) :: ibl,n,nm,iqty,jqty,iv,ix,iy,ixp,iyp,ixs,iys,
     $               ioff,joff,ibase,jbase,imat,jmat,ityp,jtyp,iqo,jqo,
     $               idg,ndg,nqpm1,nqp1p,nqpd,my
      REAL(r8), DIMENSION(:,:,:,:,:,:), POINTER :: rmat
c-----------------------------------------------------------------------
c     unload seam_save into matrix arrays.
c-----------------------------------------------------------------------
      rblock_out: DO ibl=1,nrb
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              rmat(jqty,0,0,iqty,ix,iy)=
     $          seam(ibl)%vertex(iv)%seam_save(n)
            ENDDO
          ENDDO
        ENDDO
c-----------------------------------------------------------------------
c       restore rblock higher-order diagonals.
c-----------------------------------------------------------------------
        IF (nbtype>1) THEN
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            IF (seam(ibl)%segment(iv)%h_side) THEN
              rmat=>mat%rbl_mat(ibl)%mat(2,2)%arr
              DO ibase=1,mat%rbl_mat(ibl)%nb_type(2)
                iqo=nqty*(ibase-1)
                nm=n*(ibase-1)
                DO iqty=1,nqty
                  DO jqty=1,nqty
                    nm=nm+1
                    rmat(iqo+jqty,0,0,iqo+iqty,ix,iy)=
     $                seam(ibl)%segment(iv)%seam_save(nm)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              rmat=>mat%rbl_mat(ibl)%mat(3,3)%arr
              DO ibase=1,mat%rbl_mat(ibl)%nb_type(3)
                iqo=nqty*(ibase-1)
                nm=n*(ibase-1)
                DO iqty=1,nqty
                  DO jqty=1,nqty
                    nm=nm+1
                    rmat(iqo+jqty,0,0,iqo+iqty,ix,iy)=
     $                seam(ibl)%segment(iv)%seam_save(nm)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       restore off-diagonal elements.
c-----------------------------------------------------------------------
        DO iv=1,seam(ibl)%nvert
          ix =seam(ibl)%segment(iv)%intxyn(1)
          iy =seam(ibl)%segment(iv)%intxyn(2)
          ixp=seam(ibl)%segment(iv)%intxyp(1)
          iyp=seam(ibl)%segment(iv)%intxyp(2)
          ixs=seam(ibl)%segment(iv)%intxys(1)
          iys=seam(ibl)%segment(iv)%intxys(2)
          IF (seam(ibl)%segment(iv)%h_side) THEN
            ioff=1+ixs-ix
            DO ibase=poly_deg+1,2,-1
              imat=ibase
              IF (imat>poly_deg) imat=1
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=-1,-(ibase-1),-1
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%ix0(jtyp)-1
                rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,ixs,iys)=
     $            seam(ibl)%segment(iv)%seam_mat_save
     $              (1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+ixs-ixp
            DO imat=1,poly_deg
              ityp=MIN(imat,2_i4)
              iqo=nqty*(ityp/2)*(imat-2)
              DO jbase=1,poly_deg+1-imat
                jmat=offx(jbase,imat)
                jtyp=MIN(jmat,2_i4)
                jqo=nqty*(jtyp/2)*(jmat-2)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%ix0(ityp)
                rmat(jqo+1:jqo+nqty,joff,0,iqo+1:iqo+nqty,ixs-joff,iys)=
     $            seam(ibl)%segment(iv)%seam_mat_save
     $              (1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
          ELSE
            ioff=1+iys-iy
            DO ibase=2*poly_deg,poly_deg+1,-1
              imat=ibase
              IF (imat==2*poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=-1,-(ibase-poly_deg),-1
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=mat%rbl_mat(ibl)%iy0(jtyp)-1
                rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,ixs,iys)=
     $            seam(ibl)%segment(iv)%seam_mat_save
     $              (1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
            ioff=1+iys-iyp
            DO ibase=poly_deg,2*poly_deg-1
              imat=ibase
              IF (imat==poly_deg) imat=1
              ityp=MIN(imat,3_i4)
              iqo=nqty*(ityp/3)*(imat-poly_deg-1)
              DO jbase=1,2*poly_deg-ibase
                jmat=offy(jbase,imat)
                jtyp=MIN(jmat,3_i4)
                jqo=nqty*(jtyp/3)*(jmat-poly_deg-1)
                rmat=>mat%rbl_mat(ibl)%mat(jtyp,ityp)%arr
                joff=1-mat%rbl_mat(ibl)%iy0(ityp)
                rmat(jqo+1:jqo+nqty,0,joff,iqo+1:iqo+nqty,ixs,iys-joff)=
     $            seam(ibl)%segment(iv)%seam_mat_save
     $              (1:nqty,1:nqty,ioff)
                ioff=ioff+2
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO rblock_out
c-----------------------------------------------------------------------
c     tblocks
c-----------------------------------------------------------------------
      tblock_out: DO ibl=nrb+1,ntotb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          n=0
          DO iqty=1,nqty
            DO jqty=1,nqty
              n=n+1
              mat%tbl_mat(ibl)%lmat(ix)%element(jqty,iqty,0)=
     $          seam(ibl)%vertex(iv)%seam_save(n)
            ENDDO
          ENDDO
        ENDDO
        IF (precon/='diagonal') THEN
          DO iv=1,seam(ibl)%nvert
            ix =seam(ibl)%segment(iv)%intxyn(1)
            iy =seam(ibl)%segment(iv)%intxyn(2)
            ixp=seam(ibl)%segment(iv)%intxyp(1)
            iyp=seam(ibl)%segment(iv)%intxyp(2)
            mat%tbl_mat(ibl)%lmat(ix)%element(1:nqty,1:nqty,iy)=
     $        seam(ibl)%segment(iv)%seam_mat_save(1:nqty,1:nqty,1)
            mat%tbl_mat(ibl)%lmat(ixp)%element(1:nqty,1:nqty,iyp)=
     $        seam(ibl)%segment(iv)%seam_mat_save(1:nqty,1:nqty,2)
          ENDDO
        ENDIF
      ENDDO tblock_out
c-----------------------------------------------------------------------
c     for degenerate blocks, undo the summation of on-block off-diagonal
c     matrix elements that have been communicated to neighboring
c     degenerate blocks.  original entries are in the dg* arrays.
c     a 1 in the third index indicates data for the lower corner, and a
c     2 indicates data for the upper corner.
c-----------------------------------------------------------------------
      nqpm1=nqty*(poly_deg-1)
      nqp1p=nqpm1+1
      nqpd=nqty*poly_deg

      idg=0
      DO ibl=1,nrb
        IF (.NOT.degflags(ibl)) CYCLE
        idg=idg+1
c-----------------------------------------------------------------------
c       subtract the off-diagonal corner connections to off-diagonal
c       entries along the block border.
c-----------------------------------------------------------------------
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        my=SIZE(rmat,6)-1
        rmat(1:nqty,1,0,1:nqty,0,0)=
     $    rmat(1:nqty,1,0,1:nqty,0,0)-dg11t(:,:,1,idg)
        rmat(1:nqty,1,0,1:nqty,0,my)=
     $    rmat(1:nqty,1,0,1:nqty,0,my)-dg11t(:,:,2,idg)
        rmat(1:nqty,-1,0,1:nqty,1,0)=
     $    rmat(1:nqty,-1,0,1:nqty,1,0)-dg11f(:,:,1,idg)
        rmat(1:nqty,-1,0,1:nqty,1,my)=
     $    rmat(1:nqty,-1,0,1:nqty,1,my)-dg11f(:,:,2,idg)

        IF (poly_deg>1) THEN
          rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqty,1,0,1:nqty,0,0)=rmat(1:nqty,1,0,1:nqty,0,0)-
     $        dg13(:,joff+1:joff+nqty,1,idg)
            rmat(1:nqty,1,0,1:nqty,0,my)=rmat(1:nqty,1,0,1:nqty,0,my)-
     $        dg13(:,joff+1:joff+nqty,2,idg)

            rmat(1:nqty,-1,0,1:nqty,1,0)=rmat(1:nqty,-1,0,1:nqty,1,0)-
     $        dg31(joff+1:joff+nqty,:,1,idg)
            rmat(1:nqty,-1,0,1:nqty,1,my)=rmat(1:nqty,-1,0,1:nqty,1,my)-
     $        dg31(joff+1:joff+nqty,:,2,idg)
          ENDDO

          rmat=>mat%rbl_mat(ibl)%mat(2,1)%arr
          rmat(1:nqpm1,1,0,1:nqty,0,0)=
     $      rmat(1:nqpm1,1,0,1:nqty,0,0)-dg21(:,:,1,idg)
          rmat(1:nqpm1,1,0,1:nqty,0,my)=
     $      rmat(1:nqpm1,1,0,1:nqty,0,my)-dg21(:,:,2,idg)
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqpm1,1,0,1:nqty,0,0)=
     $        rmat(1:nqpm1,1,0,1:nqty,0,0)-
     $        dg23(:,joff+1:joff+nqty,1,idg)
            rmat(1:nqpm1,1,0,1:nqty,0,my)=
     $        rmat(1:nqpm1,1,0,1:nqty,0,my)-
     $        dg23(:,joff+1:joff+nqty,2,idg)
          ENDDO

          rmat=>mat%rbl_mat(ibl)%mat(1,2)%arr
          rmat(1:nqty,-1,0,1:nqpm1,1,0)=
     $      rmat(1:nqty,-1,0,1:nqpm1,1,0)-dg12(:,:,1,idg)
          rmat(1:nqty,-1,0,1:nqpm1,1,my)=
     $      rmat(1:nqty,-1,0,1:nqpm1,1,my)-dg12(:,:,2,idg)
          DO jbase=1,poly_deg-1
            joff=(jbase-1)*nqty
            rmat(1:nqty,-1,0,1:nqpm1,1,0)=
     $        rmat(1:nqty,-1,0,1:nqpm1,1,0)-
     $        dg32(joff+1:joff+nqty,:,1,idg)
            rmat(1:nqty,-1,0,1:nqpm1,1,my)=
     $        rmat(1:nqty,-1,0,1:nqpm1,1,my)-
     $        dg32(joff+1:joff+nqty,:,2,idg)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       replace the off-diagonal connections in corner elements.
c-----------------------------------------------------------------------
        rmat=>mat%rbl_mat(ibl)%mat(1,1)%arr
        rmat(1:nqty,1,-1,1:nqty,0,1)=dg11t(:,:,1,idg)
        rmat(1:nqty,1, 1,1:nqty,0,my-1)=dg11t(:,:,2,idg)
        rmat(1:nqty,-1, 1,1:nqty,1,0)=dg11f(:,:,1,idg)
        rmat(1:nqty,-1,-1,1:nqty,1,my)=dg11f(:,:,2,idg)
        IF (poly_deg>1) THEN
          rmat=>mat%rbl_mat(ibl)%mat(2,1)%arr
          rmat(1:nqpm1,1,-1,1:nqty,0,1)=dg21(:,:,1,idg)
          rmat(1:nqpm1,1, 1,1:nqty,0,my-1)=dg21(:,:,2,idg)
          rmat=>mat%rbl_mat(ibl)%mat(3,1)%arr
          rmat(1:nqpm1,-1,1,1:nqty,1,0)=dg31(:,:,1,idg)
          rmat(1:nqpm1,-1,0,1:nqty,1,my)=dg31(:,:,2,idg)
          rmat=>mat%rbl_mat(ibl)%mat(1,2)%arr
          rmat(1:nqty,-1, 1,1:nqpm1,1,0)=dg12(:,:,1,idg)
          rmat(1:nqty,-1,-1,1:nqpm1,1,my)=dg12(:,:,2,idg)
          rmat=>mat%rbl_mat(ibl)%mat(3,2)%arr
          rmat(1:nqpm1,-1,1,1:nqpm1,1,0)=dg32(:,:,1,idg)
          rmat(1:nqpm1,-1,0,1:nqpm1,1,my)=dg32(:,:,2,idg)
          rmat=>mat%rbl_mat(ibl)%mat(1,3)%arr
          rmat(1:nqty,1,-1,1:nqpm1,0,1)=dg13(:,:,1,idg)
          rmat(1:nqty,1, 0,1:nqpm1,0,my)=dg13(:,:,2,idg)
          rmat=>mat%rbl_mat(ibl)%mat(2,3)%arr
          rmat(1:nqpm1,1,-1,1:nqpm1,0,1)=dg23(:,:,1,idg)
          rmat(1:nqpm1,1, 0,1:nqpm1,0,my)=dg23(:,:,2,idg)
        ENDIF
      ENDDO

      DEALLOCATE(degflags)
      IF (idg>0) THEN
        DEALLOCATE(dg11t,dg11f)
        IF (poly_deg>1) DEALLOCATE(dg12,dg21,dg13,dg31,dg23,dg32)
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_mat_rest
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_real_fac

c-----------------------------------------------------------------------
c     subprogram 17. iter_pre_real
c     do the preconditioning operations on the residual.  a symmetric
c     averaging is used where sqrt(1/# of blocks owning vertex) is
c     multiplied with the residual before the block-based operation
c     and again after the operation.
c-----------------------------------------------------------------------
      SUBROUTINE iter_pre_real(mat_str,fac_str,resp,zeep,adrp,nqty,nrb,
     $                         ntotb,precon)
      USE local
      USE matrix_type_mod
      USE factor_type_mod
      USE vector_type_mod
      USE edge
      USE seam_storage_mod
      USE iter_real_direct
      USE iter_real_line
      USE iter_real_fac
      USE iter_utils
      IMPLICIT NONE

      TYPE(global_matrix_type), INTENT(IN) :: mat_str
      TYPE(matrix_factor_type), INTENT(IN) :: fac_str
      INTEGER(i4), INTENT(IN) :: nqty,nrb,ntotb
      CHARACTER(*), INTENT(IN) :: precon
      TYPE(vector_type), DIMENSION(ntotb), INTENT(INOUT) :: resp,zeep,
     $                   adrp

      INTEGER(i4) :: ibl,iv,ix,iy,iq,rsstart,itype_dst,ibase,nv,itype,
     $               nqd,mxb,myb,poly_deg,poly_d2
      CHARACTER(64) :: msg
      LOGICAL :: xym_set,use_int
      REAL(r8), DIMENSION(:,:,:,:), POINTER, CONTIGUOUS :: pmatptr
      REAL(r8), DIMENSION(:,:,:), ALLOCATABLE :: zeeptr,resptr
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     return a copy of resp in zeep if there is no preconditioning.
c-----------------------------------------------------------------------
      IF (precon=='no prec') THEN
        DO ibl=1,ntotb
          zeep(ibl)=resp(ibl)
        ENDDO
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     set a module flag that indicates whether the matrix includes
c     data at cell interiors.  duplication with iter_pre_real allows
c     direct calls to iter_pre_real.
c-----------------------------------------------------------------------
      use_int=.NOT.mat_str%eliminated
c-----------------------------------------------------------------------
c     determine the number of basis functions.
c-----------------------------------------------------------------------
      IF (nrb>0) THEN
        ibase=mat_str%rbl_mat(1)%nbtype
        IF (ibase>1) THEN
          poly_deg=SUM(mat_str%rbl_mat(1)%nb_type(1:2))
        ELSE
          poly_deg=1
        ENDIF
        poly_d2=poly_deg**2
      ELSE
c-PRE
        poly_deg=1
        poly_d2=1
      ENDIF
c-----------------------------------------------------------------------
c     save the residual at the block boundaries, then multiply by
c     the averaging factor.  this is only appropriate for block-based
c     preconditioning options.  thus, when global preconditioning is
c     used in the rblocks, only tblock-tblock borders get averaged
c     here.  this is set up in the ave_factor_pre variable.
c-----------------------------------------------------------------------
      IF (precon=='gl_diaga'.OR.direct_check(precon)) THEN
        rsstart=nrb+1
      ELSE
        rsstart=1
      ENDIF
      DO ibl=rsstart,ntotb
        DO iv=1,seam(ibl)%nvert
          ix=seam(ibl)%vertex(iv)%intxy(1)
          iy=seam(ibl)%vertex(iv)%intxy(2)
          seam(ibl)%vertex(iv)%seam_save(1:nqty)=resp(ibl)%arr(:,ix,iy)
          resp(ibl)%arr(:,ix,iy)=resp(ibl)%arr(:,ix,iy)*
     $      seam(ibl)%vertex(iv)%ave_factor_pre
          IF (poly_deg>1) THEN
            ix=seam(ibl)%segment(iv)%intxys(1)
            iy=seam(ibl)%segment(iv)%intxys(2)
            iq=1
            IF (seam(ibl)%segment(iv)%h_side) THEN
              DO ibase=1,poly_deg-1
                seam(ibl)%segment(iv)%seam_save(iq:iq+nqty-1)=
     $            resp(ibl)%arrh(:,ibase,ix,iy)
                resp(ibl)%arrh(:,ibase,ix,iy)=
     $            resp(ibl)%arrh(:,ibase,ix,iy)*
     $              seam(ibl)%segment(iv)%ave_factor_pre
                iq=iq+nqty
              ENDDO
            ELSE
              DO ibase=1,poly_deg-1
                seam(ibl)%segment(iv)%seam_save(iq:iq+nqty-1)=
     $            resp(ibl)%arrv(:,ibase,ix,iy)
                resp(ibl)%arrv(:,ibase,ix,iy)=
     $            resp(ibl)%arrv(:,ibase,ix,iy)*
     $              seam(ibl)%segment(iv)%ave_factor_pre
                iq=iq+nqty
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     find 'preconditioned' residual.
c-----------------------------------------------------------------------
      IF (direct_check(precon)) THEN
        CALL solve_dir(fac_str,fac_str%spp,zeep,resp,nqty,poly_deg,
     $                 ntotb,precon)
      ELSE IF (precon(1:7)=='bl_diag') THEN
        DO ibl=1,nrb
          CALL iter_ld_solve(fac_str%bl_fac(ibl),zeep(ibl),resp(ibl),
     $                       nqty,precon,poly_deg,poly_d2,use_int)
        ENDDO
      ELSE IF (precon=='gl_diaga') THEN
        CALL iter_gl_ld_solve(fac_str,nqty,nrb,zeep,resp,adrp,poly_deg,
     $                        poly_d2,ntotb,use_int)
      ENDIF
c-----------------------------------------------------------------------
c     apply diagonal preconditioning where needed.
c
c     [don't delete (:,:,:,:) in pmatptr assigment.  this needs to be
c     interpretted as an array segment.]
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        IF (precon=='diagonal'.OR.ibl>nrb) THEN
          itype_dst=1
        ELSE IF (precon(2:7)=='l_diag'.OR.direct_check(precon)) THEN
          CYCLE
        ELSE
          itype_dst=2
        ENDIF
        mxb=fac_str%bl_fac(ibl)%mx
        myb=fac_str%bl_fac(ibl)%my
        DO itype=itype_dst,SIZE(fac_str%bl_fac(ibl)%pmat)+itype_dst-1
          pmatptr=>fac_str%bl_fac(ibl)%pmat(itype)%arr(:,:,:,:)
          IF (.NOT.ASSOCIATED(pmatptr)) CYCLE
          IF (ibl<=nrb) THEN
            nqd=mat_str%rbl_mat(ibl)%nq_type(itype)
          ELSE
            nqd=nqty
          ENDIF
          SELECT CASE(itype)
          CASE(1)
            ALLOCATE(zeeptr(nqd,1:mxb+1,1:myb+1))
            ALLOCATE(resptr(nqd,1:mxb+1,1:myb+1))
            resptr=resp(ibl)%arr
          CASE(2)
            ALLOCATE(zeeptr(nqd,1:mxb,1:myb+1))
            ALLOCATE(resptr(nqd,1:mxb,1:myb+1))
            resptr=RESHAPE(resp(ibl)%arrh,(/nqd,mxb,myb+1_i4/))
          CASE(3)
            ALLOCATE(zeeptr(nqd,1:mxb+1,1:myb))
            ALLOCATE(resptr(nqd,1:mxb+1,1:myb))
            resptr=RESHAPE(resp(ibl)%arrv,(/nqd,mxb+1_i4,myb/))
          CASE(4)
            ALLOCATE(zeeptr(nqd,1:mxb,1:myb))
            ALLOCATE(resptr(nqd,1:mxb,1:myb))
            resptr=RESHAPE(resp(ibl)%arri,(/nqd,mxb,myb/))
          END SELECT
          nv=SIZE(zeeptr,2)
          DO iy=1,SIZE(zeeptr,3)
            DO ix=1,nv
              DO iq=1,nqd
                zeeptr(iq,ix,iy)=
     $            SUM(pmatptr(:,iq,ix,iy)*resptr(:,ix,iy))
              ENDDO
            ENDDO
          ENDDO
          SELECT CASE(itype)
          CASE(1)
            zeep(ibl)%arr=zeeptr
          CASE(2)
            zeep(ibl)%arrh=
     $        RESHAPE(zeeptr,(/nqty,nqd/nqty,mxb,myb+1_i4/))
          CASE(3)
            zeep(ibl)%arrv=
     $        RESHAPE(zeeptr,(/nqty,nqd/nqty,mxb+1_i4,myb/))
          CASE(4)
            zeep(ibl)%arri=RESHAPE(zeeptr,(/nqty,nqd/nqty,mxb,myb/))
          END SELECT
          DEALLOCATE(zeeptr,resptr)
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     for block-based preconditioning, restore the residual, and
c     average at block boundaries.
c
c     we now assume that direct solves are global, so averaging is not
c     needed.
c-----------------------------------------------------------------------
      IF (.NOT.direct_check(precon)) THEN
        DO ibl=1,rsstart-1
          CALL edge_load_arr(zeep(ibl),nqty,poly_deg-1_i4,seam(ibl))
          DO iv=1,seam(ibl)%nvert
            seam(ibl)%vertex(iv)%seam_in(1:nqty)=
     $        seam(ibl)%vertex(iv)%seam_in(1:nqty)*
     $        seam(ibl)%vertex(iv)%ave_factor_pre
            IF (poly_deg>1) THEN
              seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))=
     $          seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))*
     $          seam(ibl)%segment(iv)%ave_factor_pre
            ENDIF
          ENDDO
        ENDDO
        DO ibl=rsstart,ntotb
          CALL edge_load_arr(zeep(ibl),nqty,poly_deg-1_i4,seam(ibl))
          DO iv=1,seam(ibl)%nvert
            ix=seam(ibl)%vertex(iv)%intxy(1)
            iy=seam(ibl)%vertex(iv)%intxy(2)
            resp(ibl)%arr(:,ix,iy)=
     $        seam(ibl)%vertex(iv)%seam_save(1:nqty)
            seam(ibl)%vertex(iv)%seam_in(1:nqty)=
     $        seam(ibl)%vertex(iv)%seam_in(1:nqty)*
     $        seam(ibl)%vertex(iv)%ave_factor_pre
            IF (poly_deg>1) THEN
              ix=seam(ibl)%segment(iv)%intxys(1)
              iy=seam(ibl)%segment(iv)%intxys(2)
              iq=1
              IF (seam(ibl)%segment(iv)%h_side) THEN
                DO ibase=1,poly_deg-1
                  resp(ibl)%arrh(:,ibase,ix,iy)=
     $              seam(ibl)%segment(iv)%seam_save(iq:iq+nqty-1)
                  iq=iq+nqty
                ENDDO
              ELSE
                DO ibase=1,poly_deg-1
                  resp(ibl)%arrv(:,ibase,ix,iy)=
     $              seam(ibl)%segment(iv)%seam_save(iq:iq+nqty-1)
                  iq=iq+nqty
                ENDDO
              ENDIF
              seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))=
     $          seam(ibl)%segment(iv)%seam_in(1:nqty*(poly_deg-1))*
     $          seam(ibl)%segment(iv)%ave_factor_pre
            ENDIF
          ENDDO
        ENDDO
        CALL edge_network(nqty,0_i4,poly_deg-1_i4,.false.)
        DO ibl=1,ntotb
          CALL edge_unload_arr(zeep(ibl),nqty,poly_deg-1_i4,seam(ibl))
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_pre_real
