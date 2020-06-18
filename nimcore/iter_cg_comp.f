c-----------------------------------------------------------------------
c     file iter_cg_comp.f
c     module containing a conjugate gradient solver with either a simple
c     nxn-block diagonal preconditioner or domain-decomposition
c     preconditioner for complex matrices and vectors.  this module
c     mirrors the real arithmetic version found in iter_cg_f90.f
c     12/11/98  C. Sovinec
c
c     the off-diagonal factor used to improve diagonal dominance in the
c     preconditioners is now being applied within the diagonal point-
c     block (vector block) for the line solves (not ILU).  this allows
c     better convergence on ill-conditioned matrices.
c     5/11/99  C. Sovinec
c
c     a bug in a line diagonal routine, setting up the line
c     matrices when periodic, has been fixed.  this makes a big
c     difference in performance, and 'bl_diaga' now does much better
c     than 'bl_ilu_1' in many cases with periodic blocks.
c     5/24/99  C. Sovinec
c
c     a global line Jacobi preconditioning option has been added.  the
c     option is 'gl_diaga', and the iter_gl_ld_* routines are used
c     to factor and solve the lines, which extend across as many
c     adjacent rblocks as possible.  they use the new parallel_line*
c     routines for both serial and parallel communication of line data.
c     this communication represents decomposition swaps similar to
c     the processor layers used for computing Fourier components and
c     FFTs/pseudo-spectral operations in parallel at different stages.
c     the residual and preconditioned residual are not averaged along
c     block borders, since the operation is redundant on
c     different processors.  in addition, the communication
c     and restoration of matrix elements along block borders has been
c     modified, since seaming is done within this preconditioner
c     factorization, and seam_in gets over-written.  the new seam_save
c     array has been added to the seams for this purpose.
c     6/25/99  C. Sovinec
c
c     the array orders for vectors and matrices have been switched to
c     the quantity-first arrangement.  the seaming arrays and calls
c     have also been changed in prepration for higher-order elements.
c     9/30/99  C. Sovinec
c
c     the solution and rhs vectors are now assumed to have the same
c     number of quantities as the matrix for optimization.
c     2/29/00  C. Sovinec
c
c     the block and global line Jacobi preconditioning schemes have a
c     new elimination step.  Along each line, data from the cell
c     interiors and from the segments parallel to the line being solved
c     is eliminated to reduce the 1D matrices to block-tridiagonal,
c     regardless of basis function degree.
c     4/10/00  C. Sovinec
c
c     provisions have been made for matrix computations with cell-
c     interior data eliminated prior to calling the iter_solve routine.
c     whether or not this data is used in the iterative solve dependends
c     on the matrix structure flag, eliminated.
c     4/14/00 C. Sovinec
c
c     the basis index for side and interior centered data now follows
c     the quantity index, i.e. it's the 2nd index.  rblock matrices
c     have also been changed so that all horizontal side basis fall
c     under a single basis type, and the quantity index is expanded
c     to have all bases within the type.  the same has been done for
c     vertical sides and interiors.
c     6/9/00 C. Sovinec
c
c     routines have been broken into separate modules to facilitate
c     compilation on some machines.  the conjugate gradient routines
c     now appear in the last module in this file; preconditioner
c     modules are first.
c     12/13/00 C. Sovinec
c
c     information for performing a preconditioning pass is now computed
c     or passed into iter_pre_comp, so that it may be called from
c     routines other than iter_solve.
c     1/22/01 C. Sovinec
c
c     a restart capability and a check to compare the actual and
c     recursive residuals have been added.  there is also a little
c     reorganization to make the residual computation and residual norm
c     computation more modular.
c     5/8/02, C. Sovinec
c
c     based on discussions with Mark Fahey at ORNL and Nathan Wichmann
c     of Cray Inc., the number of function calls at the inner level
c     of loops has been reduced to help optimization.
c     5/8/03, C. Sovinec
c
c     the block direct solve (via lapack) has been upgraded to handle
c     high-order elements and the degenerate point.
c     5/23/03, C. Sovinec
c
c     the lapack solve has been modified to cover all rblocks, and
c     options to substitute the sparse solvers, Sequential SuperLU
c     and SuperLU_DIST, have been added.  SSLU is a serial library,
c     whereas SLU_D performs parallel factorizations and solves.
c     [The SuperLU packages were written by Xiaoye Li (LBL),
c     James Demmel (UC-Berkeley), and John Gilbert (UC-SB).  To find
c     documentation, see the ACTS Collection web page at NERSC,
c     http://acts.nersc.gov/superlu/]
c     6/18/03, C. Sovinec
c
c     the norm, residual, and dot product routines have been moved
c     into a separate module and file, iter_utils, so that they are
c     generally available.
c     9/10/03, C. Sovinec
c
c     the matrix loads into compressed-column format for SLU have been
c     changed to read the columns of the nimrod matrices directly,
c     instead of taking the conjugate of the rows.  this allows us
c     to use the same routine for non-Hermitian matrices (with something
c     other than cg); though, there may be a slight computational
c     penalty for skipping around in memory.
c     9/26/03, C. Sovinec
c
c     the data for SLU is getting switched to compressed-row format so
c     that we can use the less memory intensive interface for SLU_DIST.
c     this is also a better ordering for nimrod data.
c     8/30/06, C. Sovinec
c
c     diagonal preconditioning for non-Hermitian matrices and a no-
c     preconditioning option have been added for GMRES tests.
c     also, the nqty=1 special case has been removed for diagonal
c     preconditioning, because it interfered with poly_degree>2.
c     6/7/07,  C. Sovinec
c
c     communication for the distributed-memory SLU_DIST interface has
c     finally been completed.  it uses point-to-point mpi routines in
c     fac_dir, and the new iter_pardir_setup routine is called from
c     iter_dir_alloc to initialize communication information.
c     this interface uses less memory, but it also leads to slower
c     factorizations and solves, so it is now a separate preconditioner
c     (solver parameter) option of 'slu_dstm.'  the full-storage
c     interface is still called when 'slu_dist' is specified.  it
c     necessarily uses compressed-column format, requiring the
c     transposing copy in fac_dir.  however, communication for this
c     option has been changed to pt-to-pt for completing columns
c     from different processors and allgatherv instead of allreduce
c     for assembling full-storage.  finally, the sequential option now
c     uses compressed-row format with some minor changes in the
c     c_fortran bridge routines.
c     11/5/07, C. Sovinec
c
c     the iter_init routine now performs a seaming operation
c     to ensure that the provided guess has absolutely no discrepancies
c     among different images along block borders.  This prevents an
c     unusual error that otherwise appears in some computations.
c     11/26/07, C. Sovinec
c
c     the upgrade to sequential SuperLU3 requires compressed column
c     format.
c     7/3/09, J. King
c
c     The sparsity pattern of the *_factor_type structures has been
c     moved into a new sparsity_pattern type. This type contains only
c     the (integer) sparsity pattern common to both factor types. This
c     allowed the iter_dirfac_alloc routines to be moved out of the
c     iter_cg_*.f files and combined as a single routine in the
c     iter_util.f file. The iter_dirfac_alloc has become a control
c     routine where its old contents have been split between an
c     initial routine assign_node_indices (common to all allocations)
c     which assigns a global index to each node and
c     alloc_sparsity_pattern which creates the j_acc, start_acc and
c     start_loc arrays and sets up the communication pattern. The new
c     distributed sparsity pattern uses the routine
c     alloc_sparsity_pattern_dist. Both of the the
c     alloc_sparsity_pattern routines call the pardir_setup routine
c     (previously present in iter_utils.f) which has had some minor
c     modifications.
c
c     With respect to the sparsity_pattern type there are two new
c     additions: sendrowst which holds the starting row for each send,
c     and sendstacc which is the equivalent of start_acc for sent data.
c     The irw array is only needed during allocation and has been
c     removed from the sparsity_pattern type and is now local to the
c     alloc_sparsity_pattern routines.
c
c     The deallocation routine, iter_dirfac_dealloc, has also been moved
c     from the iter_cg_*.f files and into the iter_util.f file. Together
c     these routines in the iter_util.f file have become a new module
c     iter_dir_fac of which only the _alloc and _dealloc routines are
c     public.
c     8/13/13, J. King
c
c     The ilu routines were never upgraded for high-order elements and
c     have not been used for years.  They have been removed, but the
c     the iter_solve_herm routine is retained with the 'ilu' module.
c     8/26/13, C. Sovinec
c
c     The 2D preconditioning routine has been doing block averaging
c     after the preconditioning operation, even when the operation is
c     a global direct solve.  The averaging is now removed when
c     direct_check is true to eliminate unnecessary seaming, since the
c     direct solves in use today are all global.
c     6/10/16, C. Sovinec
c
c     The traditional preconditioning routines and the driver,
c     iter_pre_real have been moved into a separate file,
c     iter_precon_real.f.  The driver is now an external subroutine,
c     and it is now named iter_cg_c2d_solve.
c     1/13/19, C. Sovinec
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     module iter_cg_comp
c     1. iter_cg_c2d_solve.
c     2. iter_init.
c     3. iter_dealloc.
c-----------------------------------------------------------------------
c     module containing the conjugate gradient solver.
c-----------------------------------------------------------------------
      MODULE iter_cg_comp
      USE local
      USE matrix_type_mod
      USE factor_type_mod
      USE vector_type_mod
      USE matrix_mod
      USE edge
      USE seam_storage_mod
      USE iter_comp_direct
      USE iter_comp_line
      USE iter_comp_fac
      USE iter_utils
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     module-specific definition.
c-----------------------------------------------------------------------
      INTEGER(i4), PRIVATE :: nrb,ntotb,poly_deg,poly_d2
      LOGICAL, PRIVATE :: use_int

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. iter_cg_c2d_solve.
c     controls the conjugate gradient iterations.
c-----------------------------------------------------------------------
      SUBROUTINE iter_cg_c2d_solve(mat_str,fac_str,sln,rhs,nqty,tol,
     $                             maxit,precon,err,its,seed)
      USE mpi_nim
      USE pardata
      USE time

      TYPE(complex_matrix_type), INTENT(IN) :: mat_str
      TYPE(complex_factor_type), INTENT(IN) :: fac_str
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: sln,rhs
      INTEGER(i4), INTENT(IN) :: nqty,maxit
      REAL(r8), INTENT(IN) :: tol
      CHARACTER(*), INTENT(IN) :: precon
      INTEGER(i4), INTENT(OUT) :: its
      REAL(r8), INTENT(OUT) :: err
      CHARACTER(8), INTENT(OUT) :: seed

      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res,zee,dir,adr,
     $  slnl
      INTEGER(i4) :: ibl,it,irst
      INTEGER(i4) :: ix,iy,iv,mv,iq
      INTEGER(i4), PARAMETER :: nrst=10000
      COMPLEX(r8) :: alpha,beta,ztr,ztro,ctmp
      REAL(r8) :: er0,erl,ercheck
      REAL(r8), PARAMETER :: resid_tol=1.,small=TINY(1._r4)
      CHARACTER(64) :: msg
      REAL(r8) :: tmp,timestart_cg,timeend_cg
      INTEGER(i4) :: ierror
c-----------------------------------------------------------------------
c     allocate storage; determine best guess and initial error.
c-----------------------------------------------------------------------
      CALL timer(timestart_cg)
      nrb=SIZE(mat_str%rbl_mat)
      ntotb=nrb+SIZE(mat_str%tbl_mat)
      CALL iter_init(mat_str,fac_str,sln,rhs,nqty,precon,err,er0,
     $               tol,seed,res,zee,dir,adr,slnl)
c-----------------------------------------------------------------------
c     begin cg iterations.
c-----------------------------------------------------------------------
      erl=10*err
      ztro=1
      beta=0
      irst=0
      itloop: DO it=1,maxit
c-----------------------------------------------------------------------
c       check actual residual and restart condition.
c       if both err and ercheck < tol, there is no need to make
c       corrections.
c-----------------------------------------------------------------------
        IF (MODULO(irst,nrst)==0.OR.err<=tol) THEN
          CALL iter_resid(mat_str,sln,adr,rhs,nqty,ntotb,nrb,poly_deg)
          CALL iter_inf_norm(adr,ercheck,er0,ntotb,nrb,poly_deg,use_int)
        ENDIF
        IF (MODULO(irst,nrst)==0.OR.(err<=tol.AND.ercheck>tol)) THEN
          IF (ercheck>erl) THEN
            DO ibl=1,ntotb
              sln(ibl)=slnl(ibl)
              dir(ibl)=0
              ztro=1
              beta=0
            ENDDO
            CALL iter_resid(mat_str,sln,res,rhs,nqty,ntotb,nrb,poly_deg)
            CALL iter_inf_norm(res,err,er0,ntotb,nrb,poly_deg,use_int)
c-TMP
c           IF (node==0) WRITE(nim_wr,'(/,a,i5,a,es11.4)')
c    $        "Iter_cg_comp: restarting at it ",it," err=",err
          ELSE
            IF (ABS(ercheck-err)/MAX(MIN(err,ercheck),small)>resid_tol
     $          .OR.ercheck<err/(1+resid_tol)) THEN
c-TMP
c             IF (node==0) WRITE(nim_wr,'(/,a,i5,2(a,es9.2))')
c    $          "Iter_cg_comp: correcting residual it ",
c    $          it," recursive=",err," actual=",ercheck
              DO ibl=1,ntotb
                dir(ibl)=0
                res(ibl)=adr(ibl)
              ENDDO
              ztro=1
              beta=0
              err=ercheck
            ENDIF
            DO ibl=1,ntotb
              slnl(ibl)=sln(ibl)
            ENDDO
          ENDIF
          erl=ercheck
          irst=0
        ENDIF
        irst=irst+1
c-----------------------------------------------------------------------
c       check tolerance.
c-----------------------------------------------------------------------
        IF (err<=tol) THEN
          its=it-1
          CALL iter_dealloc(res,zee,dir,adr,slnl)
c         CALL mpi_barrier(comm_layer,ierror)
          CALL timer(timeend_cg)
          time_iter = time_iter + timeend_cg-timestart_cg
          RETURN
        ENDIF
c-----------------------------------------------------------------------
c       find 'preconditioned' residual.
c-----------------------------------------------------------------------
        CALL iter_pre_comp(mat_str,fac_str,res,zee,adr,nqty,nrb,ntotb,
     $                     precon)
c-----------------------------------------------------------------------
c       find the new orthogonality coefficient.
c-----------------------------------------------------------------------
        ztr=0
        DO ibl=1,ntotb
          CALL iter_dot(ztr,zee(ibl),res(ibl),seam(ibl),ibl,nrb,
     $                  poly_deg,use_int)
        ENDDO
        IF (nprocs_layer > 1) THEN
          CALL mpi_allreduce(ztr,ctmp,1,mpi_nim_comp,mpi_sum,
     $         comm_layer,ierror)
          ztr=ctmp
        ENDIF
        beta=ztr/ztro
        ztro=ztr
c-----------------------------------------------------------------------
c       find the new orthogonal direction.
c-----------------------------------------------------------------------
        DO ibl=1,ntotb
          dir(ibl)%arr=zee(ibl)%arr+beta*dir(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            dir(ibl)%arrh=zee(ibl)%arrh+beta*dir(ibl)%arrh
            dir(ibl)%arrv=zee(ibl)%arrv+beta*dir(ibl)%arrv
            IF (use_int) THEN
              dir(ibl)%arri=zee(ibl)%arri+beta*dir(ibl)%arri
            ENDIF
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       find the product of the full matrix and the new direction.
c-----------------------------------------------------------------------
        CALL matvec(mat_str,dir,adr,nqty)
c-----------------------------------------------------------------------
c       find the step size for the new direction.
c-----------------------------------------------------------------------
        alpha=0
        DO ibl=1,ntotb
          CALL iter_dot(alpha,dir(ibl),adr(ibl),seam(ibl),ibl,nrb,
     $                  poly_deg,use_int)
        ENDDO
        IF (nprocs_layer > 1) THEN
          CALL mpi_allreduce(alpha,ctmp,1,mpi_nim_comp,mpi_sum,
     $         comm_layer,ierror)
          alpha = ctmp
        ENDIF
        alpha=ztr/alpha
c-----------------------------------------------------------------------
c       update the residual, error, and solution.
c-----------------------------------------------------------------------
        DO ibl=1,ntotb
          sln(ibl)%arr=sln(ibl)%arr+alpha*dir(ibl)%arr
          res(ibl)%arr=res(ibl)%arr-alpha*adr(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            sln(ibl)%arrh=sln(ibl)%arrh+alpha*dir(ibl)%arrh
            res(ibl)%arrh=res(ibl)%arrh-alpha*adr(ibl)%arrh
            sln(ibl)%arrv=sln(ibl)%arrv+alpha*dir(ibl)%arrv
            res(ibl)%arrv=res(ibl)%arrv-alpha*adr(ibl)%arrv
            IF (use_int) THEN
              sln(ibl)%arri=sln(ibl)%arri+alpha*dir(ibl)%arri
              res(ibl)%arri=res(ibl)%arri-alpha*adr(ibl)%arri
            ENDIF
          ENDIF
        ENDDO
        CALL iter_inf_norm(res,err,er0,ntotb,nrb,poly_deg,use_int)
c-----------------------------------------------------------------------
c     complete iteration loop.
c-----------------------------------------------------------------------
      ENDDO itloop
      its=maxit
c-----------------------------------------------------------------------
c     if err>tol return the residual in the solution vector for
c     diagnosis.
c-----------------------------------------------------------------------
      CALL iter_resid(mat_str,sln,adr,rhs,nqty,ntotb,nrb,poly_deg)
      CALL iter_inf_norm(adr,err,er0,ntotb,nrb,poly_deg,use_int)
      IF (err>tol) THEN
        DO ibl=1,ntotb
          sln(ibl)%arr=res(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            sln(ibl)%arrh=res(ibl)%arrh
            sln(ibl)%arrv=res(ibl)%arrv
            IF (use_int) sln(ibl)%arri=res(ibl)%arri
          ENDIF
        ENDDO
      ENDIF
      CALL iter_dealloc(res,zee,dir,adr,slnl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
c     CALL mpi_barrier(comm_layer,ierror)
      CALL timer(timeend_cg)
      time_iter = time_iter + timeend_cg-timestart_cg

      RETURN
      END SUBROUTINE iter_cg_c2d_solve
c-----------------------------------------------------------------------
c     subprogram 2. iter_init.
c     allocates storage used during conjugate gradient solve.
c     also initialize preconditioner from full matrix.
c-----------------------------------------------------------------------
      SUBROUTINE iter_init(mat,fac,sln,rhs,nqty,precon,err,er0,tol,seed,
     $                     res,zee,dir,adr,slnl)
      USE mpi_nim
      USE pardata

      TYPE(complex_matrix_type), INTENT(IN) :: mat
      TYPE(complex_factor_type), INTENT(IN) :: fac
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: sln,rhs
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res,zee,dir,adr,
     $   slnl
      INTEGER(i4), INTENT(IN) :: nqty
      CHARACTER(*), INTENT(IN) :: precon
      REAL(r8), INTENT(OUT) :: err,er0
      REAL(r8), INTENT(IN) :: tol
      CHARACTER(8), INTENT(OUT) :: seed

      INTEGER(i4) :: ibl,iv,mx,my,iq,jq,ix,iy,n,nd,mv,ierror,ibase
      REAL(r8) :: er_sl,er_pc,er_0v,maxrhs,maxres,maxrpc
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     determine the number of basis functions.
c-----------------------------------------------------------------------
      IF (nrb>0) THEN
        ibase=mat%rbl_mat(1)%nbtype
        IF (ibase>1) THEN
          poly_deg=SUM(mat%rbl_mat(1)%nb_type(1:2))
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
c     set a module flag that indicates whether the matrix includes
c     data at cell interiors.
c-----------------------------------------------------------------------
      use_int=.NOT.mat%eliminated
c-----------------------------------------------------------------------
c     ensure that the guess has exactly the same value in different
c     images of the same node along block borders.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL edge_load_2D_carr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
        DO iv=1,seam(ibl)%nvert
          seam(ibl)%vertex(iv)%seam_cin(1:nqty)=
     $      seam(ibl)%vertex(iv)%seam_cin(1:nqty)*
     $      seam(ibl)%vertex(iv)%ave_factor
          IF (poly_deg>1) THEN
            seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1))=
     $        seam(ibl)%segment(iv)%seam_cin(1:nqty*(poly_deg-1))*
     $        seam(ibl)%segment(iv)%ave_factor
          ENDIF
        ENDDO
      ENDDO
      CALL edge_network(nqty,1_i4,poly_deg-1_i4,.false.)
      DO ibl=1,ntotb
        CALL edge_unload_2D_carr(sln(ibl),nqty,poly_deg-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     allocate storage for the cg operations.
c-----------------------------------------------------------------------
      ALLOCATE(res(ntotb),zee(ntotb),dir(ntotb),adr(ntotb),slnl(ntotb))
      DO ibl=1,ntotb
        mx=SIZE(rhs(ibl)%arr,2)-1
        my=SIZE(rhs(ibl)%arr,3)-1
        IF (SIZE(rhs(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_init: rhs dimension inconsistent with nqty')
        IF (SIZE(sln(ibl)%arr,1)/=nqty) CALL nim_stop
     $    ('Iter_init: sln dimension inconsistent with nqty')
        CALL vector_type_alloc(res(ibl),poly_deg,mx,my,nqty)
        CALL vector_type_alloc(zee(ibl),poly_deg,mx,my,nqty)
        CALL vector_type_alloc(dir(ibl),poly_deg,mx,my,nqty)
        CALL vector_type_alloc(adr(ibl),poly_deg,mx,my,nqty)
        CALL vector_type_alloc(slnl(ibl),poly_deg,mx,my,nqty)
      ENDDO
c-----------------------------------------------------------------------
c     determine the best initial guess from 1) the zero vector,
c     2) what is supplied, and 3) applying the preconditioner to the
c     rhs.
c
c     for 3) use the preconditioner to find a first guess,
c     then find the respective residual (the guess is in zee and the
c     respective residual is in adr).
c-----------------------------------------------------------------------
      IF (.NOT.direct_check(precon)) THEN
        DO ibl=1,ntotb
          res(ibl)%arr=rhs(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            res(ibl)%arrh=rhs(ibl)%arrh
            res(ibl)%arrv=rhs(ibl)%arrv
            IF (use_int) res(ibl)%arri=rhs(ibl)%arri
          ENDIF
        ENDDO
        CALL iter_pre_comp(mat,fac,res,zee,adr,nqty,nrb,ntotb,precon)
      ENDIF
c-----------------------------------------------------------------------
c     for 2) compute the residual with the supplied solution.
c-----------------------------------------------------------------------
      CALL iter_resid(mat,sln,res,rhs,nqty,ntotb,nrb,poly_deg)
c-----------------------------------------------------------------------
c     if the preconditioner is a direct solve, use the guess before
c     preconditioning for option 3).
c-----------------------------------------------------------------------
      IF (direct_check(precon)) THEN
        CALL iter_pre_comp(mat,fac,res,zee,adr,nqty,nrb,ntotb,precon)
        DO ibl=1,ntotb
          CALL vector_add(zee(ibl),sln(ibl))
        ENDDO
      ENDIF
      CALL iter_resid(mat,zee,adr,rhs,nqty,ntotb,nrb,poly_deg)
      DO ibl=1,ntotb
        dir(ibl)=0
      ENDDO
c-----------------------------------------------------------------------
c     check errors and find reference.
c-----------------------------------------------------------------------
      er_sl=0._r8; er_pc=0._r8; er_0v=0._r8
      CALL iter_inf_norm(res,er_sl,1._r8,ntotb,nrb,poly_deg,use_int)
      CALL iter_inf_norm(rhs,er_0v,1._r8,ntotb,nrb,poly_deg,use_int)
      CALL iter_inf_norm(adr,er_pc,1._r8,ntotb,nrb,poly_deg,use_int)
      er0=MAX(er_sl,er_0v,SQRT(TINY(er0)))
c-----------------------------------------------------------------------
c     select the best guess and copy respective residual.
c-----------------------------------------------------------------------
      err=MIN(er_sl,er_pc,er_0v)
      IF (err==er_sl) THEN
c       do nothing but keep this solution and residual.
        seed='guess'
      ELSE IF (err==er_pc) THEN
        seed='pinv.rhs'
        DO ibl=1,ntotb
          sln(ibl)%arr=zee(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            sln(ibl)%arrh=zee(ibl)%arrh
            sln(ibl)%arrv=zee(ibl)%arrv
            IF (use_int) sln(ibl)%arri=zee(ibl)%arri
          ENDIF
          res(ibl)=adr(ibl)
        ENDDO
      ELSE IF (err==er_0v) THEN
        seed='0 vector'
        DO ibl=1,ntotb
          sln(ibl)%arr=0
          res(ibl)%arr=rhs(ibl)%arr
          IF (poly_deg>1.AND.ibl<=nrb) THEN
            sln(ibl)%arrh=0
            sln(ibl)%arrv=0
            res(ibl)%arrh=rhs(ibl)%arrh
            res(ibl)%arrv=rhs(ibl)%arrv
            IF (use_int) THEN
              sln(ibl)%arri=0
              res(ibl)%arri=rhs(ibl)%arri
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      err=err/er0
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_init
c-----------------------------------------------------------------------
c     subprogram 3. iter_dealloc.
c     deallocates storage used during conjugate gradient solve.
c-----------------------------------------------------------------------
      SUBROUTINE iter_dealloc(res,zee,dir,adr,slnl)

      INTEGER(i4) :: ibl
      TYPE(cvector_2D_type), DIMENSION(:), POINTER :: res,zee,dir,adr,
     $  slnl
c-----------------------------------------------------------------------
c     deallocate storage.
c-----------------------------------------------------------------------
      DO ibl=1,ntotb
        CALL vector_type_dealloc(res(ibl))
        CALL vector_type_dealloc(zee(ibl))
        CALL vector_type_dealloc(dir(ibl))
        CALL vector_type_dealloc(adr(ibl))
        CALL vector_type_dealloc(slnl(ibl))
      ENDDO
      DEALLOCATE(res,zee,dir,adr,slnl)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE iter_dealloc
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE iter_cg_comp
