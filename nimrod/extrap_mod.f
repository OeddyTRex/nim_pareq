c-----------------------------------------------------------------------
c     file extrap_mod.f
c     module containing definitions and subprograms for extrapolating
c     solutions.  they are used to produce seeds for the iterative
c     solver.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. extrap_mod.
c     1. extrap_init.
c     2. extrap_init_array.
c     3. extrap_sln.
c     4. extrap_eval_cvec2D.
c     5. extrap_eval_vec.
c     6. extrap_eval_cvec.
c     7. extrap_update.
c     8. extrap_index.
c     9. extrap_correct.
c-----------------------------------------------------------------------
c     0. module declarations.
c-----------------------------------------------------------------------
      MODULE extrap_mod
      USE local
      IMPLICIT NONE

c-----------------------------------------------------------------------
c     types used for storing solutions at multiple previous times. 
c-----------------------------------------------------------------------
      TYPE :: extrap_type
        COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: field
        COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: fieldh,fieldv,
     $               fieldi
      END TYPE extrap_type
      TYPE :: extrap_eqn_type
        TYPE(extrap_type), DIMENSION(:), POINTER :: exeqn
      END TYPE extrap_eqn_type
c-----------------------------------------------------------------------
c     variables local to extrap_mod.
c
c     extrap_t keeps the times at which data is stored.
c     extrap_q is a flag identifying each field and time split.
c     extrap_nq is the vector component index dimension for a field.
c     extrap_int is a logical that determines whether cell-centered
c                data is stored and extrapolated.  (it's not if this
c                data is eliminated before a matrix solve.)
c-----------------------------------------------------------------------
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: extrap_t
      CHARACTER(9), DIMENSION(15) :: extrap_q
      INTEGER(i4), DIMENSION(15) :: extrap_nq
      LOGICAL, DIMENSION(15) :: extrap_int
      INTEGER(i4) :: neqn
      TYPE(extrap_eqn_type), DIMENSION(:), POINTER :: extrap_bl

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. extrap_init.
c     allocates arrays and initializes field data for the extrapolation
c     routine.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_init
      USE global
      USE fields
      USE input

      INTEGER(i4) :: exdim,it,ibl,lx,ly,iq,iy,ns,ni,nq
      CHARACTER(8) :: pass
c-----------------------------------------------------------------------
c     set equation-pointer arrays.
c-----------------------------------------------------------------------
      neqn=0
c-----------------------------------------------------------------------
c     velocity
c-----------------------------------------------------------------------
      IF (nonlinear.AND.
     $    (advect=='V only'.OR.advect=='all').AND..NOT.impladv) THEN
        neqn=neqn+1
        extrap_q(neqn)='v pre'
        extrap_nq(neqn)=3
        IF (nonlinear.AND.(continuity=='full'.OR.par_visc>0)) THEN
          extrap_int(neqn)=.true.
        ELSE
          extrap_int(neqn)=.false.
        ENDIF
      ENDIF
      neqn=neqn+1
      extrap_q(neqn)='v cor'
      extrap_nq(neqn)=3
      IF (nonlinear.AND.
     $    (continuity=='full'.OR.par_visc>0.OR.impladv)) THEN
        extrap_int(neqn)=.true.
      ELSE
        extrap_int(neqn)=.false.
      ENDIF
      IF (split_visc) THEN
        neqn=neqn+1
        extrap_q(neqn)='visc'
        extrap_nq(neqn)=3
        extrap_int(neqn)=.false.
      ENDIF
c-----------------------------------------------------------------------
c     particle continuity
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        IF (nonlinear.AND..NOT.impladv) THEN
          neqn=neqn+1
          extrap_q(neqn)='n pre'
          extrap_nq(neqn)=1
          extrap_int(neqn)=.false.
        ENDIF
        neqn=neqn+1
        extrap_q(neqn)='n cor'
        extrap_nq(neqn)=1
        IF (nonlinear.AND.impladv) THEN
          extrap_int(neqn)=.true.
        ELSE
          extrap_int(neqn)=.false.
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field (mhd advance)
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN
        IF (nonlinear.AND..NOT.impladv) THEN
          neqn=neqn+1
          extrap_q(neqn)='bmhd pre'
          extrap_nq(neqn)=3
          IF (threedeta) THEN
            extrap_int(neqn)=.true.
          ELSE
            extrap_int(neqn)=.false.
          ENDIF
        ENDIF
        neqn=neqn+1
        extrap_q(neqn)='bmhd cor'
        extrap_nq(neqn)=3
        IF (threedeta.OR.nonlinear.AND.impladv) THEN
          extrap_int(neqn)=.true.
        ELSE
          extrap_int(neqn)=.false.
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field (implicit HMHD)
c-----------------------------------------------------------------------
      IF (ohms=='hall'.OR.ohms=='mhd&hall'.OR.ohms=='2fl') THEN
        neqn=neqn+1
        extrap_q(neqn)='bhmhd cor'
        extrap_nq(neqn)=3
        IF (nonlinear) THEN
          extrap_int(neqn)=.true.
        ELSE
          extrap_int(neqn)=.false.
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field divergence cleaner
c-----------------------------------------------------------------------
      IF (divbd>0.AND.split_divb) THEN
        neqn=neqn+1
        extrap_q(neqn)='divb diff'
        extrap_nq(neqn)=3
        extrap_int(neqn)=.false.
      ENDIF
c-----------------------------------------------------------------------
c     matrix solve for current density
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all'.OR.
     $    separate_pe.AND.nonlinear) THEN
        neqn=neqn+1
        extrap_q(neqn)='jfromb'
        extrap_nq(neqn)=3
        extrap_int(neqn)=.true.
      ENDIF
c-----------------------------------------------------------------------
c     temperatures (separate_pe implies Vi is computed with equilibrium
c     J contributions.)
c-----------------------------------------------------------------------
      IF (beta>0) THEN
        IF (nonlinear.AND.
     $      .NOT.(impladv.OR.separate_pe.AND.k_cross>0)) THEN
          neqn=neqn+1
          extrap_q(neqn)='ti pre'
          extrap_nq(neqn)=1
          IF (nonlinear) THEN
            extrap_int(neqn)=.true.
          ELSE
            extrap_int(neqn)=.false.
          ENDIF
        ENDIF
        neqn=neqn+1
        extrap_q(neqn)='ti cor'
        extrap_nq(neqn)=1
        IF (nonlinear) THEN
          extrap_int(neqn)=.true.
        ELSE
          extrap_int(neqn)=.false.
        ENDIF

        IF (separate_pe) THEN
          IF (nonlinear.AND..NOT.(impladv.OR.k_cross>0)) THEN
            neqn=neqn+1
            extrap_q(neqn)='te pre'
            extrap_nq(neqn)=1
            IF (nonlinear) THEN
              extrap_int(neqn)=.true.
            ELSE
              extrap_int(neqn)=.false.
            ENDIF
          ENDIF
          neqn=neqn+1
          extrap_q(neqn)='te cor'
          extrap_nq(neqn)=1
          IF (nonlinear) THEN
            extrap_int(neqn)=.true.
          ELSE
            extrap_int(neqn)=.false.
          ENDIF
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     allocate solution array.
c-----------------------------------------------------------------------
      ALLOCATE(extrap_bl(nbl))
      exdim=extrap_order+1
      DO ibl=1,nbl
        ALLOCATE(extrap_bl(ibl)%exeqn(neqn))
        IF (ibl<=nrbl) THEN
          lx=rb(ibl)%mx
          iy=0
          ly=rb(ibl)%my
        ELSE
          lx=tb(ibl)%mvert
          iy=0
          ly=0
        ENDIF
        DO iq=1,neqn
          nq=extrap_nq(iq)
          ALLOCATE(extrap_bl(ibl)%exeqn(iq)%
     $             field(nq,0:lx,iy:ly,nmodes,exdim))
          IF (poly_degree>1.AND.ibl<=nrbl) THEN
            ns=poly_degree-1
            ni=ns**2
            ALLOCATE(extrap_bl(ibl)%exeqn(iq)%
     $               fieldh(nq,ns,1:lx,0:ly,nmodes,exdim))
            ALLOCATE(extrap_bl(ibl)%exeqn(iq)%
     $               fieldv(nq,ns,0:lx,1:ly,nmodes,exdim))
            IF (extrap_int(iq)) THEN
              ALLOCATE(extrap_bl(ibl)%exeqn(iq)%
     $                 fieldi(nq,ni,1:lx,1:ly,nmodes,exdim))
            ELSE
              NULLIFY(extrap_bl(ibl)%exeqn(iq)%fieldi)
            ENDIF
          ELSE
            NULLIFY(extrap_bl(ibl)%exeqn(iq)%fieldh)
            NULLIFY(extrap_bl(ibl)%exeqn(iq)%fieldv)
            NULLIFY(extrap_bl(ibl)%exeqn(iq)%fieldi)
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     allocate time array.
c-----------------------------------------------------------------------
      ALLOCATE(extrap_t(extrap_order+1,neqn))
c-----------------------------------------------------------------------
c     call extrap_init_array to set initial values.
c-----------------------------------------------------------------------
      CALL extrap_init_array
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_init
c-----------------------------------------------------------------------
c     subprogram 2. extrap_init_array.
c     help the first couple of time steps by loading the initial
c     values for the predictor steps of the predictor-corrector
c     advance and for the implicit steps of the implicit leapfrog, which
c     use the "cor" label.
c
c     in 3d solves, the old solution is subtracted from the extrapolated
c     field to create the appropriate guess.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_init_array
      USE global
      USE fields
      USE input

      INTEGER(i4) :: exdim,it,ibl,lx,ly,iq,iy,ns,ni,nq
      CHARACTER(16) :: pass
c-----------------------------------------------------------------------
c     zero out the data array.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        DO iq=1,neqn
          extrap_bl(ibl)%exeqn(iq)%field=0
          IF (poly_degree>1.AND.ibl<=nrbl) THEN
            extrap_bl(ibl)%exeqn(iq)%fieldh=0
            extrap_bl(ibl)%exeqn(iq)%fieldv=0
            IF (extrap_int(iq)) THEN
              extrap_bl(ibl)%exeqn(iq)%fieldi=0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     set the times.
c-----------------------------------------------------------------------
      DO iq=1,neqn
        extrap_t(:,iq)=(/(t-it*dtm/1000,it=0,extrap_order)/)
      ENDDO
c-----------------------------------------------------------------------
c     velocity
c-----------------------------------------------------------------------
      IF (nonlinear.AND.
     $    (advect=='V only'.OR.advect=='all').AND..NOT.impladv) THEN
        pass='v pre'
      ELSE
        pass='v cor'
      ENDIF
      DO it=1,extrap_order+1
        DO ibl=1,nbl
          CALL extrap_update(ve(ibl),ibl,TRIM(pass))
        ENDDO 
      ENDDO 
      IF (split_visc) THEN
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(ve(ibl),ibl,'visc')
          ENDDO 
        ENDDO 
      ENDIF
c-----------------------------------------------------------------------
c     particle continuity
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        IF (nonlinear.AND..NOT.impladv) THEN
          pass='n pre'
        ELSE
          pass='n cor'
        ENDIF
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(nd(ibl),ibl,TRIM(pass))
          ENDDO 
        ENDDO 
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field (mhd advance)
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN
        IF (nonlinear.AND..NOT.impladv) THEN
          pass='bmhd pre'
        ELSE
          pass='bmhd cor'
        ENDIF
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(be(ibl),ibl,TRIM(pass))
          ENDDO 
        ENDDO 
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field (HMHD advance)
c-----------------------------------------------------------------------
      IF (ohms=='hall'.OR.ohms=='mhd&hall'.OR.ohms=='2fl') THEN
        pass='bhmhd cor'
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(be(ibl),ibl,TRIM(pass))
          ENDDO 
        ENDDO 
      ENDIF
c-----------------------------------------------------------------------
c     magnetic field divergence cleaner--always a corrector, leave 0s.
c
c     matrix solve for current density
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all'.OR.
     $    separate_pe.AND.nonlinear) THEN
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(ja(ibl),ibl,'jfromb')
          ENDDO 
        ENDDO 
      ENDIF
c-----------------------------------------------------------------------
c     temperatures (separate_pe implies Vi is computed with equilibrium
c     J contributions.)
c-----------------------------------------------------------------------
      IF (beta>0) THEN
        IF (nonlinear.AND.
     $      .NOT.(impladv.OR.separate_pe.AND.k_cross>0)) THEN
          pass='ti pre'
        ELSE
          pass='ti cor'
        ENDIF
        DO it=1,extrap_order+1
          DO ibl=1,nbl
            CALL extrap_update(tion(ibl),ibl,TRIM(pass))
          ENDDO 
        ENDDO 
        IF (separate_pe) THEN
          IF (nonlinear.AND..NOT.(impladv.OR.k_cross>0)) THEN
            pass='te pre'
          ELSE
            pass='te cor'
          ENDIF
          DO it=1,extrap_order+1
            DO ibl=1,nbl
              CALL extrap_update(tele(ibl),ibl,TRIM(pass))
            ENDDO 
          ENDDO 
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_init_array
c-----------------------------------------------------------------------
c     subprogram 3. extrap_sln.
c     extrapolates a solution with the specified polynomial order.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_sln(equation,t)
      USE fields
      USE input

      CHARACTER(*), INTENT(IN) :: equation
      REAL(r8), INTENT(IN) :: t

      REAL(r8), DIMENSION(extrap_order+1) :: dtquotient
      COMPLEX(r8), DIMENSION(:,:,:,:,:), POINTER :: fptr
      COMPLEX(r8), DIMENSION(:,:,:,:,:,:), POINTER :: fptr2
      INTEGER(i4) :: i,il,ig,iq,ibl
c-----------------------------------------------------------------------
c     find quantity index.
c-----------------------------------------------------------------------
      iq=extrap_index(equation)
c-----------------------------------------------------------------------
c     loop over previous saves to find delta-t products.
c-----------------------------------------------------------------------
      IF (extrap_order>0) THEN
        DO i=1,extrap_order+1
          dtquotient(i)=
     $      PRODUCT(t             -(/(extrap_t(il,iq),il=1,i-1),
     $         (extrap_t(ig,iq),ig=i+1,extrap_order+1)/)) /
     $      PRODUCT(extrap_t(i,iq)-(/(extrap_t(il,iq),il=1,i-1),
     $         (extrap_t(ig,iq),ig=i+1,extrap_order+1)/))
        ENDDO
      ELSE
        dtquotient=1
      ENDIF
c-----------------------------------------------------------------------
c     extrapolate the appropriate data, placing the result in the last
c     save index, where data is discarded once the new solution is
c     computed.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        fptr=>extrap_bl(ibl)%exeqn(iq)%field
        fptr(:,:,:,:,extrap_order+1)=
     $    dtquotient(extrap_order+1)*fptr(:,:,:,:,extrap_order+1)
        DO i=1,extrap_order 
          fptr(:,:,:,:,extrap_order+1)=fptr(:,:,:,:,extrap_order+1)
     $      +dtquotient(i)*fptr(:,:,:,:,i)
        ENDDO
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          fptr2=>extrap_bl(ibl)%exeqn(iq)%fieldh
          fptr2(:,:,:,:,:,extrap_order+1)=
     $      dtquotient(extrap_order+1)*fptr2(:,:,:,:,:,extrap_order+1)
          DO i=1,extrap_order 
            fptr2(:,:,:,:,:,extrap_order+1)=
     $        fptr2(:,:,:,:,:,extrap_order+1)
     $        +dtquotient(i)*fptr2(:,:,:,:,:,i)
          ENDDO
          fptr2=>extrap_bl(ibl)%exeqn(iq)%fieldv
          fptr2(:,:,:,:,:,extrap_order+1)=
     $      dtquotient(extrap_order+1)*fptr2(:,:,:,:,:,extrap_order+1)
          DO i=1,extrap_order 
            fptr2(:,:,:,:,:,extrap_order+1)=
     $        fptr2(:,:,:,:,:,extrap_order+1)
     $        +dtquotient(i)*fptr2(:,:,:,:,:,i)
          ENDDO
          IF (.NOT.extrap_int(iq)) CYCLE
          fptr2=>extrap_bl(ibl)%exeqn(iq)%fieldi
          fptr2(:,:,:,:,:,extrap_order+1)=
     $      dtquotient(extrap_order+1)*fptr2(:,:,:,:,:,extrap_order+1)
          DO i=1,extrap_order 
            fptr2(:,:,:,:,:,extrap_order+1)=
     $        fptr2(:,:,:,:,:,extrap_order+1)
     $        +dtquotient(i)*fptr2(:,:,:,:,:,i)
          ENDDO
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     shift and update the time array here, since the field update is
c     called block by block.
c-----------------------------------------------------------------------
      extrap_t(:,iq)=EOSHIFT(extrap_t(:,iq),-1)
      extrap_t(1,iq)=t
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_sln
c-----------------------------------------------------------------------
c     subprogram 4. extrap_eval_cvec2D.
c     evaluate the guess for the linear solver for a particular
c     equation in one block, and place it in a cvector_2D_type.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_eval_cvec2D(guess,ibl,im,equation)
      USE input
      USE fields
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: ibl,im
      CHARACTER(*), INTENT(IN) :: equation
      TYPE(cvector_2D_type), INTENT(INOUT) :: guess

      INTEGER(i4) :: iq,eq_len,ic,ip
c-----------------------------------------------------------------------
c     find quantity indices.
c-----------------------------------------------------------------------
      iq=extrap_index(equation)
      ip=iq-1
c-----------------------------------------------------------------------
c     if this is a corrector step for advection, the guess is the sum
c     of the predictor solution and the extrapolated difference between
c     predictor and corrector steps.  otherwise, the guess is just
c     what has been extrapolated.
c-----------------------------------------------------------------------
      IF (extrap_correct(equation)) THEN
        guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,im,extrap_order+1)
     $           +extrap_bl(ibl)%exeqn(ip)%field(:,:,:,im,1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,im,extrap_order+1)
     $              +extrap_bl(ibl)%exeqn(ip)%fieldh(:,:,:,:,im,1)
          guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,im,extrap_order+1)
     $              +extrap_bl(ibl)%exeqn(ip)%fieldv(:,:,:,:,im,1)
          IF (extrap_int(iq))
     $      guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,im,extrap_order+1)
     $                +extrap_bl(ibl)%exeqn(ip)%fieldi(:,:,:,:,im,1)
        ENDIF
      ELSE
        guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,im,extrap_order+1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,im,extrap_order+1)
          guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,im,extrap_order+1)
          IF (extrap_int(iq))
     $      guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,im,extrap_order+1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_eval_cvec2D
c-----------------------------------------------------------------------
c     subprogram 5. extrap_eval_vec.
c     evaluate the guess for the linear solver for a particular
c     equation in one block, and place the requested data in
c     a vector_type.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_eval_vec(guess,ibl,im,equation,r_i)
      USE input
      USE fields
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: ibl,im
      CHARACTER(*), INTENT(IN) :: equation,r_i
      TYPE(vector_type), INTENT(INOUT) :: guess

      INTEGER(i4) :: iq,eq_len,ic,ip
c-----------------------------------------------------------------------
c     find quantity indices.
c-----------------------------------------------------------------------
      iq=extrap_index(equation)
      ip=iq-1
c-----------------------------------------------------------------------
c     if this is a corrector step for advection, the guess is the sum
c     of the predictor solution and the extrapolated difference between
c     predictor and corrector steps.  otherwise, the guess is just
c     what has been extrapolated.
c-----------------------------------------------------------------------
      IF (extrap_correct(equation)) THEN
        SELECT CASE(r_i)
        CASE('real','REAL')
          guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,im,extrap_order+1)
     $             +extrap_bl(ibl)%exeqn(ip)%field(:,:,:,im,1)
        CASE('imag','IMAG')
          guess%arr=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                    field(:,:,:,im,extrap_order+1)
     $                   +extrap_bl(ibl)%exeqn(ip)%field(:,:,:,im,1))
        CASE('r12mi3')
          guess%arr(1:2,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                       field(1:2,:,:,im,extrap_order+1)
     $                     +extrap_bl(ibl)%exeqn(ip)%field(1:2,:,:,im,1)
          guess%arr(3,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                            field(3,:,:,im,extrap_order+1)
     $                     +extrap_bl(ibl)%exeqn(ip)%field(3,:,:,im,1))
        CASE('i12r3')
          guess%arr(1:2,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                             field(1:2,:,:,im,extrap_order+1)
     $                    +extrap_bl(ibl)%exeqn(ip)%field(1:2,:,:,im,1))
          guess%arr(3,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                     field(3,:,:,im,extrap_order+1)
     $                    +extrap_bl(ibl)%exeqn(ip)%field(3,:,:,im,1)
        CASE DEFAULT
          CALL nim_stop
     $      ('Extrap_eval_vec: '//r_i//' flag not recognized.')
        END SELECT
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          SELECT CASE(r_i)
          CASE('real','REAL')
            guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,im,extrap_order+1)
     $                +extrap_bl(ibl)%exeqn(ip)%fieldh(:,:,:,:,im,1)
            guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,im,extrap_order+1)
     $                +extrap_bl(ibl)%exeqn(ip)%fieldv(:,:,:,:,im,1)
            IF (extrap_int(iq))
     $        guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,im,extrap_order+1)
     $                  +extrap_bl(ibl)%exeqn(ip)%fieldi(:,:,:,:,im,1)
          CASE('imag','IMAG')
            guess%arrh=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldh(:,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldh(:,:,:,:,im,1))
            guess%arrv=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldv(:,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldv(:,:,:,:,im,1))
            IF (extrap_int(iq))
     $        guess%arri=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $           fieldi(:,:,:,:,im,extrap_order+1)
     $          +extrap_bl(ibl)%exeqn(ip)%fieldi(:,:,:,:,im,1))
          CASE('r12mi3')
            guess%arrh(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $         fieldh(1:2,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldh(1:2,:,:,:,im,1)
            guess%arrh(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldh(3,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldh(3,:,:,:,im,1))
            guess%arrv(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $         fieldv(1:2,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldv(1:2,:,:,:,im,1)
            guess%arrv(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldv(3,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldv(3,:,:,:,im,1))
            IF (extrap_int(iq)) THEN
              guess%arri(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $           fieldi(1:2,:,:,:,im,extrap_order+1)
     $          +extrap_bl(ibl)%exeqn(ip)%fieldi(1:2,:,:,:,im,1)
              guess%arri(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $           fieldi(3,:,:,:,im,extrap_order+1)
     $          +extrap_bl(ibl)%exeqn(ip)%fieldi(3,:,:,:,im,1))
            ENDIF
          CASE('i12r3')
            guess%arrh(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldh(1:2,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldh(1:2,:,:,:,im,1))
            guess%arrh(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $         fieldh(3,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldh(3,:,:,:,im,1)
            guess%arrv(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $         fieldv(1:2,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldv(1:2,:,:,:,im,1))
            guess%arrv(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $         fieldv(3,:,:,:,im,extrap_order+1)
     $        +extrap_bl(ibl)%exeqn(ip)%fieldv(3,:,:,:,im,1)
            IF (extrap_int(iq)) THEN
              guess%arri(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $           fieldi(1:2,:,:,:,im,extrap_order+1)
     $          +extrap_bl(ibl)%exeqn(ip)%fieldi(1:2,:,:,:,im,1))
              guess%arri(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $           fieldi(3,:,:,:,im,extrap_order+1)
     $          +extrap_bl(ibl)%exeqn(ip)%fieldi(3,:,:,:,im,1)
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Extrap_eval_vec: '//r_i//' flag not recognized.')
          END SELECT
        ENDIF
c-----------------------------------------------------------------------
c     select the appropriate data for a stand-alone guess.
c-----------------------------------------------------------------------
      ELSE
        SELECT CASE(r_i)
        CASE('real','REAL')
          guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,im,extrap_order+1)
        CASE('imag','IMAG')
          guess%arr=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                    field(:,:,:,im,extrap_order+1))
        CASE('r12mi3')
          guess%arr(1:2,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                       field(1:2,:,:,im,extrap_order+1)
          guess%arr(3,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                            field(3,:,:,im,extrap_order+1))
        CASE('i12r3')
          guess%arr(1:2,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                             field(1:2,:,:,im,extrap_order+1))
          guess%arr(3,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                     field(3,:,:,im,extrap_order+1)
        CASE DEFAULT
          CALL nim_stop
     $      ('Extrap_eval_vec: '//r_i//' flag not recognized.')
        END SELECT
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          SELECT CASE(r_i)
          CASE('real','REAL')
            guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,im,extrap_order+1)
            guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,im,extrap_order+1)
            IF (extrap_int(iq))
     $        guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,im,extrap_order+1)
          CASE('imag','IMAG')
            guess%arrh=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                       fieldh(:,:,:,:,im,extrap_order+1))
            guess%arrv=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                       fieldv(:,:,:,:,im,extrap_order+1))
            IF (extrap_int(iq))
     $        guess%arri=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                         fieldi(:,:,:,:,im,extrap_order+1))
          CASE('r12mi3')
            guess%arrh(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                            fieldh(1:2,:,:,:,im,extrap_order+1)
            guess%arrh(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                                fieldh(3,:,:,:,im,extrap_order+1))
            guess%arrv(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                            fieldv(1:2,:,:,:,im,extrap_order+1)
            guess%arrv(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                                fieldv(3,:,:,:,im,extrap_order+1))
            IF (extrap_int(iq)) THEN
              guess%arri(1:2,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                              fieldi(1:2,:,:,:,im,extrap_order+1)
              guess%arri(3,:,:,:)=-AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                              fieldi(3,:,:,:,im,extrap_order+1))
            ENDIF
          CASE('i12r3')
            guess%arrh(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                              fieldh(1:2,:,:,:,im,extrap_order+1))
            guess%arrh(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                          fieldh(3,:,:,:,im,extrap_order+1)
            guess%arrv(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                              fieldv(1:2,:,:,:,im,extrap_order+1))
            guess%arrv(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                          fieldv(3,:,:,:,im,extrap_order+1)
            IF (extrap_int(iq)) THEN
              guess%arri(1:2,:,:,:)=AIMAG(extrap_bl(ibl)%exeqn(iq)%
     $                              fieldi(1:2,:,:,:,im,extrap_order+1))
              guess%arri(3,:,:,:)=extrap_bl(ibl)%exeqn(iq)%
     $                            fieldi(3,:,:,:,im,extrap_order+1)
            ENDIF
          CASE DEFAULT
            CALL nim_stop
     $        ('Extrap_eval_vec: '//r_i//' flag not recognized.')
          END SELECT
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_eval_vec
c-----------------------------------------------------------------------
c     subprogram 6. extrap_eval_cvec.
c     evaluate the guess for the linear solver for a particular
c     equation in one block, and place it in a cvector_type.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_eval_cvec(guess,ibl,equation)
      USE input
      USE fields
      USE vector_type_mod

      INTEGER(i4), INTENT(IN) :: ibl
      CHARACTER(*), INTENT(IN) :: equation
      TYPE(cvector_type), INTENT(INOUT) :: guess

      INTEGER(i4) :: iq,eq_len,ic,ip
c-----------------------------------------------------------------------
c     find quantity indices.
c-----------------------------------------------------------------------
      iq=extrap_index(equation)
      ip=iq-1
c-----------------------------------------------------------------------
c     if this is a corrector step for advection, the guess is the sum
c     of the predictor solution and the extrapolated difference between
c     predictor and corrector steps.  otherwise, the guess is just
c     what has been extrapolated.
c-----------------------------------------------------------------------
      IF (extrap_correct(equation)) THEN
        guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,:,extrap_order+1)
     $           +extrap_bl(ibl)%exeqn(ip)%field(:,:,:,:,1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,:,extrap_order+1)
     $              +extrap_bl(ibl)%exeqn(ip)%fieldh(:,:,:,:,:,1)
          guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,:,extrap_order+1)
     $              +extrap_bl(ibl)%exeqn(ip)%fieldv(:,:,:,:,:,1)
          IF (extrap_int(iq))
     $      guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,:,extrap_order+1)
     $                +extrap_bl(ibl)%exeqn(ip)%fieldi(:,:,:,:,:,1)
        ENDIF
      ELSE
        guess%arr=extrap_bl(ibl)%exeqn(iq)%
     $              field(:,:,:,:,extrap_order+1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          guess%arrh=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldh(:,:,:,:,:,extrap_order+1)
          guess%arrv=extrap_bl(ibl)%exeqn(iq)%
     $                 fieldv(:,:,:,:,:,extrap_order+1)
          IF (extrap_int(iq))
     $      guess%arri=extrap_bl(ibl)%exeqn(iq)%
     $                   fieldi(:,:,:,:,:,extrap_order+1)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_eval_cvec
c-----------------------------------------------------------------------
c     subprogram 7. extrap_update.
c     shift old data block by block, and save the new solution in the
c     field array.
c-----------------------------------------------------------------------
      SUBROUTINE extrap_update(new,ibl,equation)
      USE input
      USE fields
      USE vector_type_mod

      TYPE(cvector_type), INTENT(IN) :: new
      INTEGER(i4), INTENT(IN) :: ibl
      CHARACTER(*), INTENT(IN) :: equation

      INTEGER(i4) :: iq,eq_len,ic,ip
c-----------------------------------------------------------------------
c     find quantity index and array index limits.
c-----------------------------------------------------------------------
      iq=extrap_index(equation)
      ip=iq-1
c-----------------------------------------------------------------------
c     shift old data.
c-----------------------------------------------------------------------
      IF (extrap_order>0) THEN
        extrap_bl(ibl)%exeqn(iq)%field(:,:,:,:,:)=
     $    EOSHIFT(extrap_bl(ibl)%exeqn(iq)%field(:,:,:,:,:),-1,dim=5)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
         extrap_bl(ibl)%exeqn(iq)%fieldh(:,:,:,:,:,:)=
     $    EOSHIFT(extrap_bl(ibl)%exeqn(iq)%fieldh(:,:,:,:,:,:),-1,dim=6)
         extrap_bl(ibl)%exeqn(iq)%fieldv(:,:,:,:,:,:)=
     $    EOSHIFT(extrap_bl(ibl)%exeqn(iq)%fieldv(:,:,:,:,:,:),-1,dim=6)
         IF (extrap_int(iq))
     $    extrap_bl(ibl)%exeqn(iq)%fieldi(:,:,:,:,:,:)=
     $    EOSHIFT(extrap_bl(ibl)%exeqn(iq)%fieldi(:,:,:,:,:,:),-1,dim=6)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     if this is a corrector step for advection, save the difference
c     between the new field and the predictor step, not just the new
c     field.
c-----------------------------------------------------------------------
      IF (extrap_correct(equation)) THEN
        extrap_bl(ibl)%exeqn(iq)%field(:,:,:,:,1)=
     $    new%arr-extrap_bl(ibl)%exeqn(ip)%field(:,:,:,:,1)
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          extrap_bl(ibl)%exeqn(iq)%fieldh(:,:,:,:,:,1)=
     $      new%arrh-extrap_bl(ibl)%exeqn(ip)%fieldh(:,:,:,:,:,1)
          extrap_bl(ibl)%exeqn(iq)%fieldv(:,:,:,:,:,1)=
     $      new%arrv-extrap_bl(ibl)%exeqn(ip)%fieldv(:,:,:,:,:,1)
          IF (extrap_int(iq))
     $      extrap_bl(ibl)%exeqn(iq)%fieldi(:,:,:,:,:,1)=
     $      new%arri-extrap_bl(ibl)%exeqn(ip)%fieldi(:,:,:,:,:,1)
        ENDIF
      ELSE
        extrap_bl(ibl)%exeqn(iq)%field(:,:,:,:,1)=new%arr
        IF (ibl<=nrbl.AND.poly_degree>1) THEN
          extrap_bl(ibl)%exeqn(iq)%fieldh(:,:,:,:,:,1)=new%arrh
          extrap_bl(ibl)%exeqn(iq)%fieldv(:,:,:,:,:,1)=new%arrv
          IF (extrap_int(iq))
     $      extrap_bl(ibl)%exeqn(iq)%fieldi(:,:,:,:,:,1)=new%arri
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE extrap_update
c-----------------------------------------------------------------------
c     subprogram 8. extrap_index.
c     return the equation-pointer for the extrapolation fields.
c-----------------------------------------------------------------------
      FUNCTION extrap_index(eqn) RESULT(ind)

      CHARACTER(*), INTENT(IN) :: eqn
      INTEGER(i4) :: ind

      INTEGER(i4) :: i
c-----------------------------------------------------------------------
c     check equation list.
c-----------------------------------------------------------------------
      DO i=1,neqn
        IF (eqn==extrap_q(i)) THEN
          ind=i
          RETURN
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     error message.
c-----------------------------------------------------------------------
      CALL nim_stop('Extrap_index cannot find '//TRIM(eqn)//'.')
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION extrap_index
c-----------------------------------------------------------------------
c     subprogram 9. extrap_correct.
c     check the equation flag and input to determine if the step is a
c     corrector step.
c-----------------------------------------------------------------------
      FUNCTION extrap_correct(eqn) RESULT(cor_step)
      USE input
      USE global

      CHARACTER(*), INTENT(IN) :: eqn
      LOGICAL :: cor_step
c-----------------------------------------------------------------------
c     the logic here must correspond to that in subroutine
c     advance_pc, as cor_step=T refers to the corrector step of
c     predictor/corrector advances.
c-----------------------------------------------------------------------
      cor_step=.false.
      IF (eqn=='divb diff') cor_step=.true.
      IF (impladv.OR.eqn=='divb diff') RETURN

      SELECT CASE(eqn)
      CASE('v cor')
        IF (nonlinear.AND.(advect=='V only'.OR.advect=='all'))
     $    cor_step=.true.
      CASE('n cor')
        IF (nonlinear) cor_step=.true.
      CASE('ti cor')
        IF (nonlinear.AND..NOT.(separate_pe.AND.k_cross>0))
     $    cor_step=.true.
      CASE('te cor')
        IF (nonlinear.AND..NOT.k_cross>0) cor_step=.true.
      CASE('bmhd cor')
        IF (nonlinear) cor_step=.true.
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END FUNCTION extrap_correct
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE extrap_mod
