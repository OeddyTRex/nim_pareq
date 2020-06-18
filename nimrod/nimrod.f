c-----------------------------------------------------------------------
c     file nimrod.f:  contains the main nimrod program and management
c     routines for the finite-element solves.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  nimrod.
c     2.  advance_il.
c     3.  advance_lin.
c     4.  advance_pc.
c     5.  adv_b_iso.
c     6.  adv_b_nsym.
c     6.1 adv_b_hyp.
c     7.  adv_b_3deta.
c     8.  adv_b_3dnsym.
c     9.  clean_divb.
c     10. jfromb.
c     11. adv_v_clap.
c     12. adv_v_aniso.
c     13. adv_v_3dn.
c     14. adv_t_aniso.
c     15. adv_t_nsym.
c     16. adv_nd.
c     17. adv_nd_nsym.
c     18. adv_nd_3dnsym.
c     19. project_ndiff.
c     20. adv_v_hypv_matv.
c-----------------------------------------------------------------------
c     subprogram 1. nimrod.
c     main program, controls I/O, integrates equations.
c-----------------------------------------------------------------------
      PROGRAM nimrod
      USE local
      USE input
      USE physdat
      USE global
      USE fields
      USE rblock
      USE tblock
      USE surface
      USE dump
      USE edge
      USE extrap_mod
      USE mpi_nim
      USE pardata
      USE time
      USE hist_mod
      USE diagnose
      USE seam_storage_mod
      USE regularity
      IMPLICIT NONE

      INTEGER(i4) :: ierror,ibl,i,cpu_tstop=0,max_nqty,nqsav,offm,istep0
      LOGICAL :: file_stat,blk_pre,glb_pre,converged=.true.
      REAL(r8) :: cpu_tcheck
      CHARACTER(64) :: stop_msg
      LOGICAL, EXTERNAL :: direct_check
c-----------------------------------------------------------------------
c     parallel initialization
c-----------------------------------------------------------------------
      CALL mpi_init(ierror)
      CALL mpi_comm_rank(mpi_comm_world,node,ierror)
      CALL mpi_comm_size(mpi_comm_world,nprocs,ierror)
c-----------------------------------------------------------------------
c     start timer
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(time_total_start)
c-----------------------------------------------------------------------
c     initialize and write the source code version number to the output
c     file.
c-----------------------------------------------------------------------
      CALL nim_version
c-----------------------------------------------------------------------
c     read the namelist file and broadcast to all processors.
c-----------------------------------------------------------------------
      INQUIRE(FILE='nimrod.in',EXIST=file_stat)
      IF (.NOT.file_stat) THEN
        IF (node == 0)
     $    WRITE(nim_wr,*) 'The input file, nimrod.in, does not exist.'
        STOP
      ENDIF
      IF (node == 0) THEN
        OPEN(UNIT=out_unit,FILE='nimrod.out',STATUS='UNKNOWN',
     $       POSITION='APPEND')
        out_opened=.true.
        CALL read_namelist('nimrod.in',.true.)
      ENDIF
      IF (nprocs > 1) CALL broadcast_input
c-----------------------------------------------------------------------
c     set physical constants based on input if necessary.
c-----------------------------------------------------------------------
      IF (set_phys_constants) THEN
        CALL physdat_set(chrg_input,zeff_input,mi_input,me_input,
     $    gam_input,kblz_input,mu0_input,c_input)
      ELSE
        CALL physdat_set()
      ENDIF
      CALL set_2fl_coefs(ohms,advect,separate_pe,
     $                   ms(1),meomi,elementary_q)
      IF (node == 0 .AND. itflag)
     $  OPEN(UNIT=it_unit,FILE='iter.out',STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     read restart dump and complete wavenumber initialization.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        nphi=2**lphi
        IF (dealiase) THEN
          nmodes_total=nphi/3+1
        ELSE  !  only n up to nphi/2 - 1 are retained.
          nmodes_total=nphi/2
        ENDIF
      ELSE
        nphi=0
        nmodes_total=lin_nmodes
      ENDIF
      CALL dump_read(nmodes,nmodes_total,keff,keff_total,t,istep,
     $               dump_file)
      ALLOCATE(k2ef(nmodes),nindex(nmodes),nindex_total(nmodes_total))
      k2ef=keff**2
      IF (geom=='tor') THEN
        nindex=NINT(keff)
        nindex_total=NINT(keff_total)
      ELSE
        nindex=NINT(per_length*keff/twopi)
        nindex_total=NINT(per_length*keff_total/twopi)
      ENDIF
      nstop=nstep+istep
      smallnum=SQRT(TINY(smallnum))
c-----------------------------------------------------------------------
c     compute mx and my in a parallel-safe manner in case they are
c     not consistent with the fluxgrid input.
c-----------------------------------------------------------------------
      mx=0
      my=0
      DO ibl=1,nrbl
        mx=mx+rb(ibl)%mx
        my=my+rb(ibl)%my
      ENDDO
      CALL mpi_allreduce(mx,i,1,mpi_nim_int,mpi_sum,comm_layer,
     $                   ierror)
      mx=i/nybl
      CALL mpi_allreduce(my,i,1,mpi_nim_int,mpi_sum,comm_layer,
     $                   ierror)
      my=i/nxbl
c-----------------------------------------------------------------------
c     set the desired number of quadrature points, and
c     determine basis values and gradients at those points.
c-----------------------------------------------------------------------
      CALL surface_set(ngr,poly_degree,integration_formula)
      DO ibl=1,nrbl
        CALL rblock_set(ngr,poly_degree,integration_formula,rb(ibl))
        CALL rblock_basis_set(rb(ibl),(/poly_degree/),
     $                        (/poly_divv,poly_divb/),
     $                        (/poly_divv_min,poly_divb_min/),
     $                        (/poly_divv_max,poly_divb_max/),
     $                        met_spl,geom)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_set(tb(ibl))
        CALL tri_linear_get_areas(tb(ibl)%tgeom)
        CALL tblock_basis_set(tb(ibl),geom)
      ENDDO
c-----------------------------------------------------------------------
c     set the internal flag for implicit advection here, since it
c     affects several initializations.  also set internal flags
c     for 3D resistive diffusion and preconditioning.
c-----------------------------------------------------------------------
      IF (ohms/='mhd'.OR.mhdadv_alg=='centered'.OR.gyr_visc>0.OR.
     $    (eq_flow/='none'.OR.separate_pe).AND..NOT.nonlinear.OR.
     $     nd_hypd>0.OR.poly_divb>=0.OR.poly_divv>=0.OR.hyp_eta>0.OR.
     $     hyp_dbd>0)
     $  impladv=.true.
      IF (impladv.AND.split_visc) THEN
        CALL nim_stop('Split viscosity is not compatible with implicit'
     $                //' advection.')
      ENDIF
      IF (elecd>0.AND.nonlinear.AND.
     $    (eta_model=='eta full'.OR.eta_model=='chodura'))
     $  threedeta=.true.
      IF (p_model=='aniso_tdep'.AND.closure_model/='std kprp n0')
     $  closure_n0_only=.false.
c-----------------------------------------------------------------------
c     allocate dependent variables in the rb and tb structures.
c-----------------------------------------------------------------------
      CALL variable_alloc
      CALL pointer_init
c-----------------------------------------------------------------------
c     initialize data structures for parallel and serial vertex data
c     communication.
c-----------------------------------------------------------------------
      IF (hyp_eta>0._r8.OR.hyp_dbd>0._r8) THEN
        max_nqty=MAX(36_i4,6_i4*nmodes)
        nqsav=max_nqty
      ELSE
        max_nqty=MAX(9_i4,3_i4*nmodes)
        nqsav=max_nqty
      ENDIF
      CALL edge_init(max_nqty,nqsav)
      IF (nprocs>1) CALL parallel_seam_init(max_nqty)
      CALL boundary_init(geom)
      IF (hyp_eta>0._r8.OR.hyp_dbd>0._r8) THEN
        max_nqty=MAX(2_i4,MAX(36_i4,6_i4*nmodes)*(poly_degree-1_i4))
        nqsav=6_i4
      ELSE
        max_nqty=MAX(2_i4,MAX(9_i4,3_i4*nmodes)*(poly_degree-1_i4))
        nqsav=3_i4
      ENDIF
      offm=poly_degree**2+poly_degree
      CALL edge_segment_init(max_nqty,nqsav,offm)
      IF (nprocs>1)
     $  CALL parallel_seg_init(max_nqty,nqsav,offm)
c-----------------------------------------------------------------------
c     initializate communication for global line preconditioning.
c-----------------------------------------------------------------------
      IF (solver(1:2)=='bl'.OR.vmhd_solver(1:2)=='bl'.OR.
     $    bmhd_solver(1:2)=='bl'.OR.temp_solver(1:2)=='bl') THEN
        blk_pre=.true.
      ELSE
        blk_pre=.false.
      ENDIF
      IF (solver(1:2)=='gl'.OR.vmhd_solver(1:2)=='gl'.OR.
     $    bmhd_solver(1:2)=='gl'.OR.temp_solver(1:2)=='gl'.OR.
     $    direct_check(solver).OR.direct_check(vmhd_solver).OR.
     $    direct_check(bmhd_solver).OR.direct_check(temp_solver)) THEN
        glb_pre=.true.
      ELSE
        glb_pre=.false.
      ENDIF
      IF (blk_pre.AND.glb_pre) CAll nim_stop
     $  ('Solver choices must be either all global or all block-based.')
      IF (solver=='gl_diaga'.OR.vmhd_solver=='gl_diaga'.OR.
     $    bmhd_solver=='gl_diaga'.OR.temp_solver=='gl_diaga')
     $  CALL parallel_line_init(rb,nrbl,nbl,poly_degree)
c-----------------------------------------------------------------------
c     compute external block surface tangents.
c-----------------------------------------------------------------------
      CALL block_create_tang(poly_degree)
c-----------------------------------------------------------------------
c     verify that input does not conflict with an R=0 domain.
c-----------------------------------------------------------------------
      IF (any_r0blocks) THEN
        IF (nonlinear.AND..NOT.dealiase.AND.lphi<=2) CALL nim_stop
     $    ("Input lphi must be>2 with dealiase=F when domain has R=0.")
c       IF (neoe_flag/='none') CALL nim_stop
c    $    ("At least one block contacts R=0, set neo_flag=none.")
c       IF (neoi_flag/='none') CALL nim_stop
c    $    ("At least one block contacts R=0, set neo_flag=none.")
      ENDIF
c-----------------------------------------------------------------------
c     find cell-based geometric quantities.
c-----------------------------------------------------------------------
      CALL cell_init
c-----------------------------------------------------------------------
c     initialize the processor grid used by SuperLU.
c-----------------------------------------------------------------------
      IF (direct_check(solver).OR.direct_check(vmhd_solver).OR.
     $    direct_check(bmhd_solver).OR.direct_check(temp_solver)) THEN
        slu_nrowp=SQRT(REAL(nprocs_layer))
        DO
          IF (MODULO(nprocs_layer,slu_nrowp)==0) EXIT
          slu_nrowp=slu_nrowp-1
        ENDDO
        slu_ncolp=nprocs_layer/slu_nrowp
        CALL c_fortran_slugrid(1_i4,comm_layer,node_layer,slu_nrowp,
     $                         slu_ncolp,slugrid_handle)
      ENDIF
c-----------------------------------------------------------------------
c     create the mass matrix, if needed, and storage locations for
c     other matrices.
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all'.OR.
     $    separate_pe.AND.nonlinear) CALL mass_mat_init
      IF (continuity/='none'.AND.nd_correrr.AND.
     $    (nd_diff>0.OR.nd_hypd>0)) CALL proj_mat_init
      CALL matrix_init
c-----------------------------------------------------------------------
c     set the shape function for an applied heat source
c-----------------------------------------------------------------------
      CALL q_applied_init
c-----------------------------------------------------------------------
c     ensure that the initial conditions satisfy the boundary
c     conditions.
c-----------------------------------------------------------------------
      IF (loop_volt/=0.OR.i_desired/=0) CALL loop_voltage
      IF (e_vertical/=0) CALL vertical_efield
      CALL boundary_vals_init
c-----------------------------------------------------------------------
c     initialize the *_old data structures that hold pre-advance data.
c     store the equilibrium quantities and necessary derivatives at the
c     quadrature points, and interpolate the initial conditions.
c-----------------------------------------------------------------------
      CALL soln_save(converged)
      CALL quadrature_save
c-----------------------------------------------------------------------
c     find grid point for history output and open file if nhist > 0.
c     also initialize other graphics if necessary.
c-----------------------------------------------------------------------
      IF (nhist>0) CALL probe_hist_init
      IF (nprocs>1) THEN
        IF (xt_stride>0) CALL nim_stop
     $    ('Set xt_stride=0 when running in parallel and use nimplot.')
        IF (yt_stride>0) CALL nim_stop
     $    ('Set yt_stride=0 when running in parallel and use nimplot.')
        IF (xy_stride>0) CALL nim_stop
     $    ('Set xy_stride=0 when running in parallel and use nimplot.')
      ELSE
        IF (xt_stride>0 .OR. yt_stride>0) CALL time_slice_init
      ENDIF
c-----------------------------------------------------------------------
c     warnings about nonconforming elements.
c-----------------------------------------------------------------------
      IF ( (nonlinear.OR.eq_flow/='none') .AND.
     $     (advect=='V only'.OR.advect =='all') .AND. .NOT. conform)
     $  CALL nim_stop('Conforming elements are needed for V advection.')
      IF ((kin_visc>0.OR.iso_visc>0) .AND. .NOT.conform)
     $  CALL nim_stop('Conforming elements are needed for viscosity.')
      IF (lump_all.OR.lump_b)
     $  CALL nim_stop('Lumping options have been disabled.')
c-----------------------------------------------------------------------
c     initialize the perturbed current density for electron inertia.
c-----------------------------------------------------------------------
      IF (ohms=='2fl'.AND.advect=='all'.OR.
     $    separate_pe.AND.nonlinear) CALL jfromb('j pre')
c-----------------------------------------------------------------------
c     initialize bb dyad.
c-----------------------------------------------------------------------
      IF (beta>0.AND.nonlinear) THEN
        IF (p_computation=='at nodes') THEN
          CALL p_from_nt('all')
        ELSE
          CALL p_from_nt('qp sym')
        ENDIF
      ENDIF
      IF (nonlinear) CALL ave_field_check
      IF (beta>0.AND.p_model(1:5)=='aniso'.OR.
     $    par_visc>0.OR.gyr_visc>0)  CALL find_bb
c-----------------------------------------------------------------------
c     compute the cfls based on initial conditions for history
c     diagnostics.  also set initial resistivity and transport
c     coefficients that are temperature-dependent.
c-----------------------------------------------------------------------
      IF (nonlinear) THEN
        CALL ave_n_check
        CALL ave_v_check
        IF (threedeta) CALL find_eta_t
        IF (eta_model/='fixed'.AND.elecd>0._r8) CALL ave_eta_check
        IF (p_model=='aniso_plltdep'.OR.p_model=='aniso_tdep'.OR.
     $      par_visc>0.AND.parvisc_model=='plltdep') THEN
          CALL find_kappa_t
          CALL ave_kappa_check
        ENDIF
        CALL vcom_store('standard')
      ENDIF
      CALL new_dt(converged)
c-----------------------------------------------------------------------
c     initialize the particle-diffusivity error-correction factor.
c-----------------------------------------------------------------------
      IF (continuity/='none'.AND.nd_correrr.AND.
     $    (nd_diff>0.OR.nd_hypd>0))
     $  CALL project_ndiff(converged,.false.)
c-----------------------------------------------------------------------
c     initialize arrays for extrapolating guesses for the linear
c     solver. 
c-----------------------------------------------------------------------
      CALL extrap_init
c-----------------------------------------------------------------------
c     time-step loop.
c-----------------------------------------------------------------------
      CALL timer_init
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(time_loop_start)
      istep0=istep
      time_step_loop: DO
        IF (t>=tmax.OR.istep>=nstop.OR.cpu_tstop>0) EXIT
c-----------------------------------------------------------------------
c       check the divergence(B).
c-----------------------------------------------------------------------
        CALL divb_check
c-----------------------------------------------------------------------
c       diagnose and write results.
c-----------------------------------------------------------------------
        IF (converged) THEN
          IF (istep>istep0) CALL nim_output(.true.)
          istep=istep+1
        ENDIF
c-----------------------------------------------------------------------
c       advance the solution with the set of management routines
c       for the selected type of calculation.
c-----------------------------------------------------------------------
        IF (nonlinear) THEN
          IF (impladv) THEN
            CALL advance_il(converged)  !  nonlinear implicit leapfrog
          ELSE
            CALL advance_pc(converged)  !  semi-imp mhd with p/c advect
          ENDIF
        ELSE
          CALL advance_lin(converged)   !  linear computations
        ENDIF
c-----------------------------------------------------------------------
c       check if all of the iterative solves have converged.  if not,
c       soln_save will restore the fields at the beginning of the
c       time-step and decrease dt.
c-----------------------------------------------------------------------
        IF (converged) t=t+dt
        CALL soln_save(converged)
        IF (.NOT.converged) CALL vcom_store('standard')
c-----------------------------------------------------------------------
c       check the amount of cpu time used.
c-----------------------------------------------------------------------
        CALL timer(cpu_tcheck)
        IF (cpu_tcheck-time_total_start>cpu_tmax) cpu_tstop=1
        IF (nprocs>1) THEN
          CALL mpi_allreduce(cpu_tstop,i,1,mpi_nim_int,mpi_max,
     $         mpi_comm_world,ierror)
          cpu_tstop=i
        ENDIF
      ENDDO time_step_loop
c-----------------------------------------------------------------------
c     finish timing.
c-----------------------------------------------------------------------
      CALL mpi_barrier(mpi_comm_world,ierror)
      CALL timer(time_loop_end)
      CALL timer_close
c-----------------------------------------------------------------------
c     terminal diagnosis and cathartic dump.
c-----------------------------------------------------------------------
      CALL nim_output(.false.)
c-----------------------------------------------------------------------
c     expiration.
c-----------------------------------------------------------------------
      IF (cpu_tstop>0) THEN
        stop_msg='CPU time limit reached.'
      ELSE
        stop_msg='Normal termination.'
      ENDIF
      CALL nim_stop(stop_msg)
      END PROGRAM nimrod
c-----------------------------------------------------------------------
c     subprogram 2. advance_il.
c     main control routine of the nonlinear implicit leapfrog advance.
c
c     the management routines used in this time-stepping are
c
c     velocity:  adv_v_3dn   - for non-Hermitian 3D operators
c
c     number density:  adv_nd_3dnsym - non-Hermitian 3D advect operator
c
c     temperature: adv_t_aniso - handles 3D and anisotropic operators
c
c     magnetic field: adv_b_3dnsym - non-Hermitian 3D operators from
c                                    implicit advection or HMHD.
c                     adv_b_hyp    - split hyper-resistivity step.
c-----------------------------------------------------------------------
      SUBROUTINE advance_il(converged)
      USE local
      USE physdat
      USE global
      USE fields
      USE input
      USE rblock
      USE tblock
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: converged

      INTEGER(i4) :: ibl,ierror
      LOGICAL :: force_recomp,nmat_comp,conv_tmp,do_ti,do_te
      LOGICAL, SAVE :: newdart=.false.,newkarti=.false.,newkarte=.false.
c-----------------------------------------------------------------------
c     interface for set_nodal_min:
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE set_nodal_min(nb,minv,compv,realv)
        USE local
        USE vector_type_mod
        IMPLICIT NONE
        INTEGER(i4), INTENT(IN) :: nb
        REAL(r8), INTENT(IN) :: minv
        TYPE(cvector_type), DIMENSION(:), POINTER :: compv
        TYPE(vector_type), DIMENSION(:), POINTER :: realv
        END SUBROUTINE
      END INTERFACE
c-----------------------------------------------------------------------
c     replace the quadrature-point storage of time-average number
c     density, temperature, and magnetic field with the fields at the
c     end of the  previous advance.  also, multiply n and T to get
c     pressure for the si coefficient.
c-----------------------------------------------------------------------
      CALL b_store('end')
      IF (continuity/='none') CALL n_store('end',newdart)
      IF (beta>0) THEN
        CALL temp_store('ion end',newkarti)
        CALL temp_store('ele end',newkarte)
        CALL p_from_nt('qp sym')
        IF (p_computation=='at nodes'.OR.siop_type=='3D') THEN
          CALL p_from_nt('all')
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%pres,rb(ibl)%qpres,rb(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%pres,tb(ibl)%qpres,tb(ibl))
          ENDDO
        ENDIF
        IF (p_model(1:5)=='aniso'.OR.par_visc>0.OR.
     $      gyr_visc>0) CALL find_bb
        IF (p_model=='aniso_plltdep'.OR.p_model=='aniso_tdep'.OR.
     $      par_visc>0.AND.parvisc_model=='plltdep') THEN
          CALL find_kappa_t
          CALL ave_kappa_check
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     check the change in average B and P for the semi-implicit
c     operators.
c-----------------------------------------------------------------------
      CALL ave_field_check
c-----------------------------------------------------------------------
c     determine new time-step.
c-----------------------------------------------------------------------
      CALL new_dt(converged)
      IF (converged) THEN
        force_recomp=.false.
      ELSE
        force_recomp=.true.
      ENDIF

      IF (MODULO(istep,n_mat_update)==0) force_recomp=.true.
      converged=.true.
c-----------------------------------------------------------------------
c     find the applied voltages if needed.
c-----------------------------------------------------------------------
      IF (loop_volt/=0.OR.i_desired/=0) CALL loop_voltage
      IF (e_vertical/=0) CALL vertical_efield
c-----------------------------------------------------------------------
c     advance flow velocity with implicit advection.  if it is used,
c     first update the fields for the 3D semi-implicit operator.
c-----------------------------------------------------------------------
      IF (ohms/='hall') THEN
        integrand_flag='all forces cor'
        IF (siop_type=='3D') CALL si_store
c-----------------------------------------------------------------------
c       set diffusive-correction factor.
c-----------------------------------------------------------------------
        IF (continuity/='none'.AND.nd_correrr.AND.
     $      (nd_diff>0.OR.nd_hypd>0)) CALL ndiff_store('end')
        CALL adv_v_3dn('v cor    ',converged,force_recomp)
        IF (.NOT.converged) RETURN
c-----------------------------------------------------------------------
c       apply the time-split implicit hyper-viscous diffusion.
c       also, apply the Fourier damping.
c-----------------------------------------------------------------------
        IF (hyp_visc>0._r8) THEN
          CALL adv_v_hypv_matv(converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN
        ENDIF
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,ve)
c-----------------------------------------------------------------------
c       update the quadrature-point data and ve_n0 data for advection
c       coefficients.
c-----------------------------------------------------------------------
        CALL vcom_store('standard')
        CALL ave_v_check
c-----------------------------------------------------------------------
c       update the <VV> tensor if needed for artificial diffusivity.
c-----------------------------------------------------------------------
        IF (continuity/='none'.AND.nd_dart_upw>0.OR.
     $                  beta>0.AND. t_dart_upw>0) CALL find_vv
      ENDIF
c-----------------------------------------------------------------------
c     advance number density, then save the average of the old and new
c     number densities at the quadrature points.  it is used for the
c     temperature and HMHD magnetic field advances.
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        integrand_flag='n correct'
        nmat_comp=force_recomp.OR.newdart
        CALL adv_nd_3dnsym('n cor    ',converged,nmat_comp)
        IF (.NOT.converged) RETURN
        CALL n_store('check ave',newdart)

        IF (newdart) THEN
          CALL adv_nd_3dnsym('n cor dart',converged,.false.)
          IF (.NOT.converged) RETURN
        ENDIF
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,nd)
        IF (nd_nodal_floor>=0._r8)
     $    CALL set_nodal_min(nbl,nd_nodal_floor,nd,nd_eq)
        IF (newdart.OR.fmx_drate>0._r8.OR.nd_nodal_floor>=0._r8)
     $    CALL n_store('average',newdart)
        IF (nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0))
     $    CALL project_ndiff(converged,.true.)
      ENDIF
c-----------------------------------------------------------------------
c     reset the nd_n0 data for coefficients and n=0 only continuity as
c     needed.
c-----------------------------------------------------------------------
      CALL ave_n_check
c-----------------------------------------------------------------------
c     advance temperatures and pressures.
c-----------------------------------------------------------------------
      IF (beta>0) THEN
c-----------------------------------------------------------------------
c       ion temperature advance.  separate_pe implies
c       equilibrium J is used in the computation of Vi.
c
c       also save the average of the old and new temperatures at the
c       quadrature points.  it is used for electron pressure in the
c       HMHD advance of B.
c-----------------------------------------------------------------------
        IF (visc_heat) CALL viscous_heating

        integrand_flag='all ti terms'
        CALL adv_t_aniso('ti cor     ',converged,force_recomp)
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tion)
        IF (ti_nodal_floor>=0._r8)
     $    CALL set_nodal_min(nbl,ti_nodal_floor,tion,tion_eq)
        CALL temp_store('ion ave',newkarti)
c-----------------------------------------------------------------------
c       electron temperature as either a separate computation or as a
c       copy of ion temperature.
c-----------------------------------------------------------------------
        IF (separate_pe) THEN
          integrand_flag='all te terms'
          CALL adv_t_aniso('te cor     ',converged,force_recomp)
          IF (.NOT.converged) RETURN
          IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tele)
          IF (te_nodal_floor>=0._r8)
     $      CALL set_nodal_min(nbl,te_nodal_floor,tele,tele_eq)
        ELSE
          DO ibl=1,nbl
            tele(ibl)=tion(ibl)
            CALL vector_mult(tele(ibl),pe_frac/(zeff*(1-pe_frac)))
          ENDDO
        ENDIF
        CALL temp_store('ele ave',newkarte)
      ENDIF
c-----------------------------------------------------------------------
c     compute Te-dependent resistivity at the quadrature points.
c-----------------------------------------------------------------------
      IF (threedeta) CALL find_eta_t
      IF (eta_model/='fixed'.AND.elecd>0._r8) CALL ave_eta_check
c-----------------------------------------------------------------------
c     update the magnetic field and current density with the mhd
c     Ohm's law.
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN
        integrand_flag='mhd'
        CALL adv_b_3dnsym('bmhd cor',converged,force_recomp)
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,be)
      ENDIF
c-----------------------------------------------------------------------
c     update the magnetic field with a Hall-MHD Ohm's law.  if
c     current advection is used, solve for J as a nodal quantity.  this
c     is also needed for Pe advection to find the CFL condition.
c-----------------------------------------------------------------------
      IF (ohms=='hall'.OR.ohms=='mhd&hall'.OR.ohms=='2fl') THEN
        integrand_flag='hmhd'
        IF (ohms=='2fl'.AND.advect=='all') CALL jfromb('j pre')
        CALL adv_b_3dnsym('bhmhd cor',converged,force_recomp)
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,be)
      ELSE IF (separate_pe.AND.nonlinear) THEN
        CALL jfromb('j pre')
      ENDIF
c-----------------------------------------------------------------------
c     apply the split numerical hyper-resistivity step if needed.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND.split_hypeta) THEN
        integrand_flag='hyp eta'
        CALL adv_b_hyp('bhyp',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     diffuse div(B) and find the bb dyad if necessary.
c-----------------------------------------------------------------------
      IF (divbd>0.AND.MODULO(istep,ndivb)==0.AND.split_divb) THEN
        CALL b_store('end')
        CALL clean_divb(converged)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
      CALL b_store('ave')
      IF (beta>0.AND.p_model(1:5)=='aniso')  CALL find_bb
c-----------------------------------------------------------------------
c     recompute the temperature fields if needed for nonlinear centering
c     of coefficients.  leave the average temperature in the quad-point
c     storage for the energy diagnostics.
c-----------------------------------------------------------------------
      do_ti=.false.
      do_te=.false.
      IF (beta>0.AND.
     $    (t_dart_upw>0.OR.p_model(1:5)=='aniso'.OR.ohm_heat.OR.
     $     visc_heat.AND.par_visc>0)) THEN

        IF (visc_heat.AND.par_visc>0) CALL viscous_heating
        IF (nonlinear.AND.(p_model=='aniso_plltdep'.OR.
     $                     p_model=='aniso_tdep')) THEN
          CALL find_kappa_t
          CALL ave_kappa_check
        ENDIF
        CALL temp_store('ion chk',newkarti)
        CALL temp_store('ele chk',newkarte)

        integrand_flag='all ti terms'
        CALL adv_t_aniso('ti cor kart',converged,.false.)
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tion)
        IF (ti_nodal_floor>=0._r8)
     $    CALL set_nodal_min(nbl,ti_nodal_floor,tion,tion_eq)

        IF (separate_pe) THEN
          integrand_flag='all te terms'
          CALL adv_t_aniso('te cor kart',converged,.false.)
          IF (.NOT.converged) RETURN
          IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tele)
        ELSE
          DO ibl=1,nbl
            tele(ibl)=tion(ibl)
            CALL vector_mult(tele(ibl),pe_frac/(zeff*(1-pe_frac)))
          ENDDO
        ENDIF
        do_ti=.true.
        do_te=.true.
      ELSE IF (beta>0.AND.separate_pe) THEN  !  time-centered J
        CALL temp_store('ele chk',newkarte)
        integrand_flag='all te terms'
        CALL adv_t_aniso('te cor kart',converged,.false.)
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tele)
        do_te=.true.
      ENDIF
c-----------------------------------------------------------------------
c     Perform the time-split thermal equilibration and update time-
c     averaged temperatures at quad points for diagnostics.
c-----------------------------------------------------------------------
      IF (beta>0) THEN
        IF (separate_pe.AND.tequil_rate>0) THEN
          CALL thermal_equil
          do_ti=.true.
          do_te=.true.
        ENDIF
        IF (do_ti) CALL temp_store('ion ave diagn',newkarti)
        IF (do_te) THEN
          IF (te_nodal_floor>=0._r8)
     $      CALL set_nodal_min(nbl,te_nodal_floor,tele,tele_eq)
          CALL temp_store('ele ave diagn',newkarte)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advance_il
c-----------------------------------------------------------------------
c     subprogram 3. advance_lin.
c     main control routine of the linear implicit leapfrog advance.
c     this is identical to the old semi-implicit mhd advance for linear
c     computations without equilibrium flow.
c
c     the management routines used in this time-stepping are
c
c     velocity:  adv_v_aniso - anisotropic 2D operators 
c                adv_v_clap  - real isotropic 2D operators 
c
c     number density:  adv_nd_nsym - non-Hermitian or aniso 2D operator
c                      adv_nd      - real isotropic 2D operators
c
c     temperature: adv_t_nsym - handles non-Hermitian aniso 2D operators
c
c     magnetic field: adv_b_nsym - handles non-Hermitian 2D operators
c                     adv_b_iso  - real 2D diffusion operator
c                     adv_b_hyp  - 6-vector solve for hyper-resistivity
c-----------------------------------------------------------------------
      SUBROUTINE advance_lin(converged)
      USE local
      USE physdat
      USE global
      USE fields
      USE input
      USE rblock
      USE tblock
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: converged

      INTEGER(i4) :: ibl,ierror
      LOGICAL :: force_recomp,newdart,newkart,conv_tmp
c-----------------------------------------------------------------------
c     determine new time-step.
c-----------------------------------------------------------------------
      CALL new_dt(converged)
      IF (converged) THEN
        force_recomp=.false.
      ELSE
        force_recomp=.true.
      ENDIF
      converged=.true.
c-----------------------------------------------------------------------
c     replace the quadrature-point storage of number density,
c     temperature, and magnetic field with the fields at the end of the
c     previous advance.
c-----------------------------------------------------------------------
      CALL b_store('end')
      IF (continuity/='none') CALL n_store('end',newdart)
      IF (beta>0) THEN
        CALL temp_store('ion end',newkart)
        CALL temp_store('ele end',newkart)
        IF (p_computation=='at nodes') THEN
          CALL p_from_nt('all')
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%pres,rb(ibl)%qpres,rb(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%pres,tb(ibl)%qpres,tb(ibl))
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     advance flow velocity.
c-----------------------------------------------------------------------
      IF (ohms/='hall') THEN
        IF (split_visc.AND.par_visc<=0.AND..NOT.impladv) THEN
          integrand_flag='mhd and vdgv cor'
        ELSE
          integrand_flag='all forces cor'
        ENDIF
        IF (mhd_si_iso<1.OR.iso_visc>0.OR.par_visc>0.OR.impladv) THEN
          CALL adv_v_aniso('v cor    ',converged,force_recomp)
        ELSE
          CALL adv_v_clap('v cor    ',converged,force_recomp)
        ENDIF
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
c-----------------------------------------------------------------------
c       apply the time-split implicit hyper-viscous diffusion.
c-----------------------------------------------------------------------
        IF (hyp_visc>0._r8) THEN
          CALL adv_v_hypv_matv(converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN
        ENDIF
c-----------------------------------------------------------------------
c       update the quadrature-point data.
c-----------------------------------------------------------------------
        CALL vcom_store('standard')
c-----------------------------------------------------------------------
c       separate viscosity time split.
c-----------------------------------------------------------------------
        IF (split_visc.AND.par_visc<=0.AND..NOT.impladv) THEN
          integrand_flag='visc'
          CALL adv_v_clap('visc     ',converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN
          CALL vcom_store('standard')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     advance number density, then save the average of the old and new
c     number densities at the quadrature points.  it is used for the
c     temperature and HMHD magnetic field advances.
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        integrand_flag='n correct'
        IF (impladv) THEN
          CALL adv_nd_nsym('n cor    ',converged,force_recomp)
        ELSE
          CALL adv_nd('n cor    ',converged,force_recomp)
        ENDIF
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
        CALL n_store('average',newdart)
        IF (nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0))
     $    CALL project_ndiff(converged,.true.)
      ENDIF
c-----------------------------------------------------------------------
c     advance temperatures and pressures.
c-----------------------------------------------------------------------
      IF (beta>0) THEN
c-----------------------------------------------------------------------
c       ion temperature advance.  separate_pe implies
c       equilibrium J is used in the computation of Vi.
c-----------------------------------------------------------------------
        IF (visc_heat) CALL viscous_heating
        integrand_flag='all ti terms'
        CALL adv_t_nsym('ti cor    ',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
c-----------------------------------------------------------------------
c       electron temperature as either a separate computation or as a
c       copy of ion temperature.
c-----------------------------------------------------------------------
        IF (separate_pe) THEN
          integrand_flag='all te terms'
          CALL adv_t_nsym('te cor    ',converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN
        ELSE
          DO ibl=1,nbl
            tele(ibl)=tion(ibl)
            CALL vector_mult(tele(ibl),pe_frac/(zeff*(1-pe_frac)))
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       save the average of the old and new temperatures at the
c       quadrature points.  it is used for electron pressure in the
c       HMHD advance of B.
c-----------------------------------------------------------------------
        CALL temp_store('ion ave',newkart)
        CALL temp_store('ele ave',newkart)
      ENDIF
c-----------------------------------------------------------------------
c     update the magnetic field and current density with the mhd
c     Ohm's law.
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN
        integrand_flag='mhd'
        IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
          CALL adv_b_hyp('bmhd cor ',converged,force_recomp)
        ELSE IF (impladv) THEN
          CALL adv_b_nsym('bmhd cor ',converged,force_recomp)
        ELSE
          CALL adv_b_iso('bmhd cor ',converged,force_recomp)
        ENDIF
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     update the magnetic field with a Hall-MHD Ohm's law.  if
c     current advection is used, solve for J as a nodal quantity.  this
c     is also needed for Pe advection to find the CFL condition.
c-----------------------------------------------------------------------
      IF (ohms=='hall'.OR.ohms=='mhd&hall'.OR.ohms=='2fl') THEN
        integrand_flag='hmhd'
        IF (ohms=='2fl'.AND.advect=='all') CALL jfromb('j pre')
        IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
          CALL adv_b_hyp('bhmhd cor',converged,force_recomp)
        ELSE
          CALL adv_b_nsym('bhmhd cor',converged,force_recomp)
        ENDIF
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     apply the split numerical hyper-resistivity step if needed.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND.split_hypeta) THEN
        integrand_flag='hyp eta'
        CALL adv_b_hyp('bhyp',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     diffuse div(B)
c-----------------------------------------------------------------------
      IF (divbd>0.AND.MODULO(istep,ndivb)==0.AND.split_divb) THEN
        CALL b_store('end')
        CALL clean_divb(converged)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
      CALL b_store('ave')
c-----------------------------------------------------------------------
c     recompute electron temperature so that the electron-V used in 
c     the Te advance has time-centered J-contributions, consistent with
c     the implicit leapfrog algorithm.
c-----------------------------------------------------------------------
      IF (beta>0.AND.separate_pe) THEN
        integrand_flag='all te terms'
        CALL temp_store('ele chk',newkart)
        CALL adv_t_nsym('te cor rec',converged,.false.)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     Perform the time-split thermal equilibration and update time-
c     averaged temperatures at quad points for diagnostics.
c-----------------------------------------------------------------------
      IF (beta>0.AND.separate_pe) THEN
        IF (tequil_rate>0) THEN
          CALL thermal_equil
          CALL temp_store('ion ave diagn',newkart)
        ENDIF
        CALL temp_store('ele ave diagn',newkart)
      ENDIF
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advance_lin
c-----------------------------------------------------------------------
c     subprogram 4. advance_pc.
c     main control routine of nonlinear semi-implicit mhd advance
c     with predictor/corrector advection.
c
c     the management routines used in this time-stepping are
c
c     velocity:  adv_v_3dn   - for 3D operators
c                adv_v_aniso - anisotropic 2D operators 
c                adv_v_clap  - real isotropic 2D operators 
c
c     number density:  adv_nd - real isotropic 2D operators
c
c     temperature: adv_t_aniso - handles 3D and anisotropic operators
c
c     magnetic field: adv_b_3deta - for 3D resistive diffusion operator
c                     adv_b_iso   - real 2D diffusion operator
c-----------------------------------------------------------------------
      SUBROUTINE advance_pc(converged)
      USE local
      USE physdat
      USE global
      USE fields
      USE input
      USE rblock
      USE tblock
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: converged

      INTEGER(i4) :: ibl,ierror
      LOGICAL :: force_recomp,nmat_comp,conv_tmp
      LOGICAL, SAVE :: newdart=.false.,newkarti=.false.,newkarte=.false.
c-----------------------------------------------------------------------
c     interface for set_nodal_min:
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE set_nodal_min(nb,minv,compv,realv)
        USE local
        USE vector_type_mod
        IMPLICIT NONE
        INTEGER(i4), INTENT(IN) :: nb
        REAL(r8), INTENT(IN) :: minv
        TYPE(cvector_type), DIMENSION(:), POINTER :: compv
        TYPE(vector_type), DIMENSION(:), POINTER :: realv
        END SUBROUTINE
      END INTERFACE
c-----------------------------------------------------------------------
c     replace the quadrature-point storage of number density,
c     temperature, and magnetic field with the fields at the end of the
c     previous advance.  also, multiply n and T to get pressure for
c     the si coefficient.
c-----------------------------------------------------------------------
      CALL b_store('end')
      IF (continuity/='none') CALL n_store('end',newdart)
      IF (beta>0) THEN
        CALL temp_store('ion end',newkarti)
        CALL temp_store('ele end',newkarte)
        CALL p_from_nt('qp sym')
        IF (p_computation=='at nodes'.OR.siop_type=='3D') THEN
          CALL p_from_nt('all')
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%pres,rb(ibl)%qpres,rb(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%pres,tb(ibl)%qpres,tb(ibl))
          ENDDO
        ENDIF
        IF (p_model(1:5)=='aniso'.OR.par_visc>0)  CALL find_bb
        IF (p_model=='aniso_plltdep'.OR.p_model=='aniso_tdep'.OR.
     $      par_visc>0.AND.parvisc_model=='plltdep') THEN
          CALL find_kappa_t
          CALL ave_kappa_check
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     check the change in average B and P for the semi-implicit
c     operators.
c-----------------------------------------------------------------------
      CALL ave_field_check
c-----------------------------------------------------------------------
c     determine new time-step.
c-----------------------------------------------------------------------
      CALL new_dt(converged)
      IF (converged) THEN
        force_recomp=.false.
      ELSE
        force_recomp=.true.
      ENDIF

      IF (MODULO(istep,n_mat_update)==0) force_recomp=.true.
      converged=.true.
c-----------------------------------------------------------------------
c     find the applied voltages if needed.
c-----------------------------------------------------------------------
      IF (loop_volt/=0.OR.i_desired/=0) CALL loop_voltage
      IF (e_vertical/=0) CALL vertical_efield
c-----------------------------------------------------------------------
c     Compute pressure anisotopy.
c-----------------------------------------------------------------------
c     IF ((neoe_flag/='none').OR.(neoi_flag/='none'))
c    $            CALL advance_neo('explicit')
c-----------------------------------------------------------------------
c     set diffusive-correction factor.
c-----------------------------------------------------------------------
      IF (continuity/='none'.AND.nd_correrr.AND.
     $    (nd_diff>0.OR.nd_hypd>0)) CALL ndiff_store('end')
c-----------------------------------------------------------------------
c     velocity predictor step.
c-----------------------------------------------------------------------
      IF (advect=='V only'.OR.advect=='all') THEN
        IF (split_visc.AND.par_visc<=0) THEN
          integrand_flag='mhd and vdgv pre'
        ELSE
          integrand_flag='all forces pre'
        ENDIF
        IF (continuity=='full'.OR.par_visc>0.OR.siop_type=='3D') THEN
          IF (siop_type=='3D') CALL si_store
          CALL adv_v_3dn('v pre    ',converged,force_recomp)
        ELSE IF (mhd_si_iso<1.OR.iso_visc>0.OR.par_visc>0) THEN
          CALL adv_v_aniso('v pre    ',converged,force_recomp)
        ELSE
          CALL adv_v_clap('v pre    ',converged,force_recomp)
        ENDIF
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     velocity corrector step.
c-----------------------------------------------------------------------
      IF (split_visc.AND.par_visc<=0) THEN
        integrand_flag='mhd and vdgv cor'
      ELSE
        integrand_flag='all forces cor'
      ENDIF
      IF (continuity=='full'.OR.par_visc>0.OR.siop_type=='3D') THEN
        CALL adv_v_3dn('v cor    ',converged,force_recomp)
      ELSE IF (mhd_si_iso<1.OR.iso_visc>0.OR.par_visc>0) THEN
        CALL adv_v_aniso('v cor    ',converged,force_recomp)
      ELSE
        CALL adv_v_clap('v cor    ',converged,force_recomp)
      ENDIF
      IF (nlayers>1) THEN
        CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                     mpi_land,comm_mode,ierror)
        converged=conv_tmp
      ENDIF
      IF (.NOT.converged) RETURN
      IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,ve)
c-----------------------------------------------------------------------
c     apply the time-split implicit hyper-viscous diffusion.
c-----------------------------------------------------------------------
      IF (hyp_visc>0._r8) THEN
        CALL adv_v_hypv_matv(converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
c-----------------------------------------------------------------------
c     update the quadrature-point data.
c-----------------------------------------------------------------------
      CALL vcom_store('standard')
c-----------------------------------------------------------------------
c     separate viscosity time split.
c-----------------------------------------------------------------------
      IF (split_visc.AND.par_visc<=0) THEN
        integrand_flag='visc'
        CALL adv_v_clap('visc     ',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
        CALL vcom_store('standard')
      ENDIF
c-----------------------------------------------------------------------
c     advance density, then save it at the quadrature points.
c-----------------------------------------------------------------------
      IF (continuity/='none') THEN
        integrand_flag='n predict'
        nmat_comp=force_recomp.OR.newdart
        CALL adv_nd('n pre    ',converged,nmat_comp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN

        integrand_flag='n correct'
        CALL adv_nd('n cor    ',converged,nmat_comp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,nd)
        IF (nd_nodal_floor>=0._r8)
     $    CALL set_nodal_min(nbl,nd_nodal_floor,nd,nd_eq)
        CALL n_store('average',newdart)
        IF (nd_correrr.AND.(nd_diff>0.OR.nd_hypd>0))
     $    CALL project_ndiff(converged,.true.)
      ENDIF
c-----------------------------------------------------------------------
c     reset the nd_n0 data for coefficients and n=0 only continuity as
c     needed.
c-----------------------------------------------------------------------
      CALL ave_n_check
c-----------------------------------------------------------------------
c     advance temperatures and pressures.
c-----------------------------------------------------------------------
      IF (beta>0) THEN
c-----------------------------------------------------------------------
c       ion temperature prediction & correction.  separate_pe implies
c       equilibrium J is used in the computation of Vi.
c-----------------------------------------------------------------------
        IF (visc_heat) CALL viscous_heating

        integrand_flag='all ti terms pre'
        CALL adv_t_aniso('ti pre     ',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN

        integrand_flag='all ti terms'
        CALL adv_t_aniso('ti cor     ',converged,force_recomp)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
        IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tion)
        IF (ti_nodal_floor>=0._r8)
     $    CALL set_nodal_min(nbl,ti_nodal_floor,tion,tion_eq)
c-----------------------------------------------------------------------
c       electron temperature as either a separate computation or as a
c       copy of ion temperature.
c-----------------------------------------------------------------------
        IF (separate_pe) THEN
          integrand_flag='all te terms pre'
          CALL adv_t_aniso('te pre     ',converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN

          integrand_flag='all te terms'
          CALL adv_t_aniso('te cor     ',converged,force_recomp)
          IF (nlayers>1) THEN
            CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                         mpi_land,comm_mode,ierror)
            converged=conv_tmp
          ENDIF
          IF (.NOT.converged) RETURN
          IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,tele)
          IF (te_nodal_floor>=0._r8)
     $      CALL set_nodal_min(nbl,te_nodal_floor,tele,tele_eq)
          IF (tequil_rate>0) CALL thermal_equil
        ELSE
          DO ibl=1,nbl
            tele(ibl)=tion(ibl)
            CALL vector_mult(tele(ibl),pe_frac/(zeff*(1-pe_frac)))
          ENDDO
        ENDIF
        CALL temp_store('ion ave',newkarti)
        CALL temp_store('ele ave',newkarte)
      ENDIF
c-----------------------------------------------------------------------
c     compute Te-dependent resistivity at the quadrature points.
c-----------------------------------------------------------------------
      IF (threedeta) CALL find_eta_t
      IF (eta_model/='fixed'.AND.elecd>0._r8) CALL ave_eta_check
c-----------------------------------------------------------------------
c     update the magnetic field and current density with the mhd
c     Ohm's law.  the predictor step uses all terms.
c
c     B predictor step.
c-----------------------------------------------------------------------
      integrand_flag='mhd_predict'
      IF (threedeta) THEN
        CALL adv_b_3deta('bmhd pre ',converged,force_recomp)
      ELSE
        CALL adv_b_iso('bmhd pre ',converged,force_recomp)
      ENDIF
      IF (nlayers>1) THEN
        CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                     mpi_land,comm_mode,ierror)
        converged=conv_tmp
      ENDIF
      IF (.NOT.converged) RETURN
c-----------------------------------------------------------------------
c     B corrector step.
c-----------------------------------------------------------------------
      integrand_flag='mhd'
      IF (threedeta) THEN
        CALL adv_b_3deta('bmhd cor ',converged,force_recomp)
      ELSE
        CALL adv_b_iso('bmhd cor ',converged,force_recomp)
      ENDIF
      IF (nlayers>1) THEN
        CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                     mpi_land,comm_mode,ierror)
        converged=conv_tmp
      ENDIF
      IF (.NOT.converged) RETURN
      IF (fmx_drate>0._r8) CALL fourier_damp(nbl,dt,be)
c-----------------------------------------------------------------------
c     diffuse div(B) if necessary.
c-----------------------------------------------------------------------
      IF (divbd>0.AND.MODULO(istep,ndivb)==0.AND.split_divb) THEN
        CALL b_store('end')
        CALL clean_divb(converged)
        IF (nlayers>1) THEN
          CALL mpi_allreduce(converged,conv_tmp,1_i4,mpi_nim_logical,
     $                       mpi_land,comm_mode,ierror)
          converged=conv_tmp
        ENDIF
        IF (.NOT.converged) RETURN
      ENDIF
      CALL b_store('ave')
c-----------------------------------------------------------------------
c     termination.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE advance_pc
c-----------------------------------------------------------------------
c     subprogram 5. adv_b_iso
c     control-routine for magnetic field advances where the effective
c     impedance tensor is isotropic.  for a predictor
c     step, the result is averaged with the old b and saved in work1.
c     for a corrector step, the result is placed in be.
c-----------------------------------------------------------------------
      SUBROUTINE adv_b_iso(b_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE vector_type_mod
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: b_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,i_ri,ibe,iv,bits,imode,nc
      INTEGER(i4), SAVE :: istep_mhd=-1
      REAL(r8) :: err,center
      CHARACTER(8) :: seed,flag
      TYPE(global_matrix_type), DIMENSION(:), POINTER :: b_mat
      TYPE(matrix_factor_type), DIMENSION(:), POINTER :: b_fac
      LOGICAL :: new_mat,do_surf
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(b_pass,t+dt)
c-----------------------------------------------------------------------
c     for predictor steps, save the fb_vxb centered result for mhd.
c-----------------------------------------------------------------------
      IF (b_pass=='bmhd pre') THEN
        center=fb_vxb
      ENDIF
c-----------------------------------------------------------------------
c     create the matrix for each Fourier component.
c-----------------------------------------------------------------------
      new_mat=.false.
      b_mat=>bmhd_mat
      b_fac=>bmhd_fac
      nc=4
      IF ((dt/=dt_old.OR.eta_model=='eta n=0 only'.AND.eta_changed)
     $    .AND.istep>istep_mhd.OR.force_recomp) THEN
        new_mat=.true.
        istep_mhd=istep
      ENDIF
      IF (new_mat) THEN
        CALL matrix_create(b_mat,b_fac,curl_de_iso,dirichlet_op,
     $                    '3vn',b_pass(1:nc),.true.,bmhd_solver,
     $                    mass_type='none')
      ENDIF
c-----------------------------------------------------------------------
c     determine if n=0 tangential electric field should be applied.
c-----------------------------------------------------------------------
      IF (volt/=0.OR.e_vert/=0) THEN
        do_surf=.true.
      ELSE
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for the selected time split, and
c     impose homogeneous dirichlet boundary conditions on the normal
c     component of the change in B.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      CALL get_rhs(brhs_mhd ,crhs,dirichlet_rhs,'3vn','cyl_vec',
     $             do_surf,e_tangential,rmat_elim=b_mat)
c-----------------------------------------------------------------------
c     loop over the modes, which are independent in this version.
c-----------------------------------------------------------------------
      bits=0
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       separate the solves for (real_r,real_z,-imag_phi) and
c       (imag_r,imag_z,real_phi), which use the same lhs.  an exception
c       is made for n=0, where the phi component is uncoupled and all
c       real components are solved at the same time.
c-----------------------------------------------------------------------
        real_imag: DO i_ri=0,1 
          IF (i_ri==0) THEN
            flag='r12mi3'
          ELSE
            flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loop
            flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         evaluate the product of the matrix and the old magnetic field.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),be(ibl),flag,imode)
          ENDDO
          CALL matvec(b_mat(imode),vectr,sln,3_i4)
c-----------------------------------------------------------------------
c         add the product of the matrix and the old vector to the rhs,
c         then set the initial guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),crhs(ibl),flag,imode)
            IF (ASSOCIATED(sln(ibl)%arri)) sln(ibl)%arri=0
            CALL vector_add(vectr(ibl),sln(ibl))
            CALL extrap_eval_vec(sln(ibl),ibl,imode,b_pass,flag)
          ENDDO
c-----------------------------------------------------------------------
c         invert the matrix.
c-----------------------------------------------------------------------
          CALL iter_cg_2d_solve(b_mat(imode),b_fac(imode),sln,vectr,
     $                          3_i4,tol,maxit,bmhd_solver,err,its,seed)
          IF (node == 0 .AND. itflag) CALL iter_out(b_pass,seed,its,err)
          bits=its+bits
c-----------------------------------------------------------------------
c         save the solution.
c-----------------------------------------------------------------------
          convergence: IF (err>tol) THEN
            converged=.false.
            RETURN
          ELSE convergence
c-----------------------------------------------------------------------
c           determine the interior data if eliminated for the matrix
c           solve.
c-----------------------------------------------------------------------
            IF (poly_degree>1) CALL fe_postsolve(b_mat(imode),vectr,sln,
     $                                          be,crhs,3_i4,imode,flag)
            update: DO ibl=1,nbl
              SELECT CASE(b_pass)
c-----------------------------------------------------------------------
c             put the appropriately centered prediction in work1.
c-----------------------------------------------------------------------
              CASE ('bmhd pre')
                CALL vector_assign_cvec(vectr(ibl),be(ibl),flag,imode)
                CALL vector_add(sln(ibl),vectr(ibl),v1fac=center,
     $                          v2fac=1-center)
                CALL cvector_assign_vec(work1(ibl),sln(ibl),flag,imode)
c-----------------------------------------------------------------------
c             update be for this time split.
c-----------------------------------------------------------------------
              CASE ('bmhd cor')
                CALL cvector_assign_vec(be(ibl),sln(ibl),flag,imode)
c-----------------------------------------------------------------------
c           trap unexpected pass flags.
c-----------------------------------------------------------------------
              CASE DEFAULT
                CALL nim_stop
     $          (TRIM(b_pass)//' is inappropriate pass for adv_b_iso.')
              END SELECT
            ENDDO update
          ENDIF convergence
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     save iteration counts.
c-----------------------------------------------------------------------
      mhdits=mhdits+bits
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c     we're done with crhs and it has the same nq as work1, so it is
c     used for temporary storage here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        SELECT CASE(b_pass)
        CASE('bmhd pre')
          crhs(ibl)=work1(ibl)
          CALL vector_add(crhs(ibl),be(ibl),1/center,1-1/center)
          CALL extrap_update(crhs(ibl),ibl,b_pass)
        CASE('bmhd cor')
          CALL extrap_update(be(ibl),ibl,b_pass)
        END SELECT
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          SELECT CASE(b_pass)
          CASE('bmhd pre')
            CALL dirichlet_rhs(work1(ibe),seam(ibe),'3vn',3_i4)
          CASE('bmhd cor')
            CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
          END SELECT
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.
c-----------------------------------------------------------------------
      SELECT CASE(b_pass)
      CASE('bmhd pre')
        CALL regular_ave(work1,3_i4,nmodes,nindex)
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      CASE('bmhd cor')
        CALL regular_ave(be,3_i4,nmodes,nindex)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_b_iso
c-----------------------------------------------------------------------
c     subprogram 6. adv_b_nsym
c     control-routine for magnetic field advances where the effective
c     impedance tensor is non-Hermitian.
c-----------------------------------------------------------------------
      SUBROUTINE adv_b_nsym(b_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE vector_type_mod
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE iter_gmres_c2d
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: b_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,ibe,iv,imode,nbdis
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,flag
      LOGICAL :: do_surf
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(csln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        IF (poly_divb>=0) THEN
          nbdis=rb(ibl)%auxb%n_int
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes,nbdis,1_i4)
        ELSE
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes)
        ENDIF
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(b_pass,t+dt)
c-----------------------------------------------------------------------
c     create the matrix for each Fourier component.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.b0_changed.OR.n0_changed.OR.v0_changed).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        CALL matrix_create(bhmhd_cmat,bhmhd_cfac,b_hmhd_op,
     $                    dirichlet_comp_op,'3vn','bhmhd',.true.,
     $                    bmhd_solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     determine if n=0 tangential electric field should be applied.
c-----------------------------------------------------------------------
      IF (volt/=0.OR.e_vert/=0.OR.ohms/='mhd'.AND.beta>0) THEN
        do_surf=.true.
      ELSE
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side and impose homogeneous dirichlet
c     boundary conditions on the normal component of the change in B.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN   !   for implicit advection
        CALL get_rhs(brhs_mhd,crhs,dirichlet_rhs,'3vn','cyl_vec',
     $               do_surf,e_tangential,cmat_elim=bhmhd_cmat)
      ELSE
        CALL get_rhs(brhs_hmhd,crhs,dirichlet_rhs,'3vn','cyl_vec',
     $               do_surf,e_tangential,cmat_elim=bhmhd_cmat)
      ENDIF
c-----------------------------------------------------------------------
c     loop over the modes, which are independent in this version.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       evaluate the product of the matrix and the old magnetic field.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),be(ibl),imode)
        ENDDO
        CALL matvec(bhmhd_cmat(imode),cvectr,csln,3_i4)
c-----------------------------------------------------------------------
c       add the product of the matrix and the old vector to the rhs,
c       then set the initial guess.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
          IF (ASSOCIATED(csln(ibl)%arri)) csln(ibl)%arri=0
          CALL vector_add(cvectr(ibl),csln(ibl))
          CALL extrap_eval_cvec2D(csln(ibl),ibl,imode,b_pass)
        ENDDO
c-----------------------------------------------------------------------
c       invert the matrix.
c-----------------------------------------------------------------------
        CALL iter_gmr_c2d_solve(bhmhd_cmat(imode),bhmhd_cfac(imode),
     $                          csln,cvectr,3_i4,tol,maxit,bmhd_solver,
     $                          err,its,seed)
        IF (node == 0 .AND. itflag) CALL iter_out(b_pass,seed,its,err)
        hallits=its+hallits
c-----------------------------------------------------------------------
c       save the solution.
c-----------------------------------------------------------------------
        convergence: IF (err>tol) THEN
          converged=.false.
          RETURN
        ELSE convergence
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.
c-----------------------------------------------------------------------
          IF (poly_degree>1) THEN
            IF (poly_divb<0) THEN
              CALL fe_postsolve(bhmhd_cmat(imode),cvectr,csln,
     $                          be,crhs,3_i4,imode)
            ELSE
              DO ibl=1,nrbl            !  diffusive, not hyperbolic
                auxb(ibl)%arri=0._r8   !  correction, so no old value
              ENDDO
              CALL fe_postsolve(bhmhd_cmat(imode),cvectr,csln,
     $                          be,crhs,3_i4,imode,auxb)
            ENDIF
          ENDIF
          SELECT CASE(b_pass)
c-----------------------------------------------------------------------
c         update be for this time split.
c-----------------------------------------------------------------------
          CASE ('bhmhd cor','bmhd cor')
            DO ibl=1,nbl
              CALL cvector_assign_cvec2(be(ibl),csln(ibl),imode)
            ENDDO
c-----------------------------------------------------------------------
c         trap unexpected pass flags.
c-----------------------------------------------------------------------
          CASE DEFAULT
            CALL nim_stop
     $        (TRIM(b_pass)//' pass is inappropriate for adv_b_nsym.')
          END SELECT
        ENDIF convergence
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation and we're done with crhs.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL extrap_update(be(ibl),ibl,b_pass)
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.
c-----------------------------------------------------------------------
      CALL regular_ave(be,3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_b_nsym
c-----------------------------------------------------------------------
c     subprogram 6.1 adv_b_hyp
c     control-routine for magnetic field advances with an auxiliary
c     vector for numerical hyper-resistivity.
c-----------------------------------------------------------------------
      SUBROUTINE adv_b_hyp(b_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE vector_type_mod
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE iter_gmres_c2d
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: b_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,ibe,iv,imode,nbdis
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,flag
      LOGICAL :: do_surf,first_call=.true.
      TYPE(cvector_type), DIMENSION(:), POINTER, SAVE :: cvec6
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: b_mat
      TYPE(complex_factor_type), DIMENSION(:), POINTER :: b_fac
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      IF (first_call) ALLOCATE(cvec6(nbl))
      DO ibl=1,nrbl
        CALL vector_type_alloc(csln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,6_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,6_i4)
        IF (poly_divb>=0.AND.b_pass/='bhyp') THEN
          nbdis=rb(ibl)%auxb%n_int
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,6_i4,nmodes,nbdis,1_i4)
        ELSE
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,6_i4,nmodes)
        ENDIF
        IF (first_call)
     $    CALL vector_type_alloc(cvec6(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,6_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,6_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,6_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,6_i4,
     $                         nmodes)
      ENDDO
      first_call=.false.
c-----------------------------------------------------------------------
c     for time-split hyper-diffusion in nonlinear computations, update
c     magnetic field at the quadrature points from the previous split.
c     otherwise, update the polynomial fit for the linear solve.
c-----------------------------------------------------------------------
      IF (b_pass=='bhyp') THEN
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%be,rb(ibl)%qbe,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%be,tb(ibl)%qbe,tb(ibl))
        ENDDO
        b_mat=>bhyp_cmat
        b_fac=>bhyp_cfac
      ELSE
        CALL extrap_sln(b_pass,t+dt)
        b_mat=>bhmhd_cmat
        b_fac=>bhmhd_cfac
      ENDIF
c-----------------------------------------------------------------------
c     create the matrix for each Fourier component.  for linear
c     computations, use b_hmhd_op for the unsplit advance.  for
c     nonlinear computations, use b_hyp_op for the split advance.
c
c     set all components of the auxiliary field to zero at the boundary
c     to control energy flux.
c-----------------------------------------------------------------------
      IF (dt/=dt_old.AND.istep>istep_mat.OR.force_recomp) THEN
        IF (b_pass/='bhyp') THEN
          CALL matrix_create(b_mat,b_fac,b_hmhd_op,
     $                       dirichlet_comp_op,'3vn sd sd sd','bhmhd',
     $                       .true.,bmhd_solver,mass_type='none')
        ELSE
          CALL matrix_create(b_mat,b_fac,b_hyp_op,
     $                       dirichlet_comp_op,'3vn sd sd sd','bhyp',
     $                       .true.,bmhd_solver,mass_type='none')
        ENDIF
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     determine if n=0 tangential electric field should be applied.
c-----------------------------------------------------------------------
      IF ((volt/=0.OR.e_vert/=0.OR.ohms/='mhd'.AND.beta>0).AND.
     $    b_pass/='bhyp') THEN
        do_surf=.true.
      ELSE
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side and impose homogeneous dirichlet
c     boundary conditions on the normal component of the change in B.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      IF (b_pass=='bhyp') THEN     !   split hyper-resistive step
        CALL get_rhs(brhs_hyp,crhs,dirichlet_rhs,'3vn sd sd sd',
     $               'cyl_vec',do_surf,e_tangential,cmat_elim=b_mat)
      ELSE IF (ohms=='mhd') THEN   !   for implicit advection
        CALL get_rhs(brhs_mhd,crhs,dirichlet_rhs,'3vn sd sd sd',
     $               'cyl_vec',do_surf,e_tangential,cmat_elim=b_mat)
      ELSE
        CALL get_rhs(brhs_hmhd,crhs,dirichlet_rhs,'3vn sd sd sd',
     $               'cyl_vec',do_surf,e_tangential,cmat_elim=b_mat)
      ENDIF
c-----------------------------------------------------------------------
c     load the cvec6 structure with B at the start of the time-split,
c     padded with 0s in the rows of the auxiliary field.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        cvec6(ibl)=0._r8
        CALL cvector_assignq_cvec(cvec6(ibl),be(ibl),3_i4)
      ENDDO
c-----------------------------------------------------------------------
c     loop over the modes, which are independent in this version.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       evaluate the product of the matrix and the old magnetic field.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),cvec6(ibl),imode)
        ENDDO
        CALL matvec(b_mat(imode),cvectr,csln,6_i4)
c-----------------------------------------------------------------------
c       add the product of the matrix and the old vector to the rhs,
c       then set the initial guess.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
          IF (ASSOCIATED(csln(ibl)%arri)) csln(ibl)%arri=0
          CALL vector_add(cvectr(ibl),csln(ibl))
          csln(ibl)=0._r8
        ENDDO
c-----------------------------------------------------------------------
c       invert the matrix.
c-----------------------------------------------------------------------
        CALL iter_gmr_c2d_solve(b_mat(imode),b_fac(imode),csln,cvectr,
     $                          6_i4,tol,maxit,bmhd_solver,err,its,seed)
        IF (node == 0 .AND. itflag) CALL iter_out(b_pass,seed,its,err)
        hallits=its+hallits
c-----------------------------------------------------------------------
c       save the solution.
c-----------------------------------------------------------------------
        convergence: IF (err>tol) THEN
          converged=.false.
          RETURN
        ELSE convergence
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.
c-----------------------------------------------------------------------
          IF (poly_degree>1) THEN
            IF (poly_divb<0.OR.b_pass=='bhyp') THEN
              CALL fe_postsolve(b_mat(imode),cvectr,csln,
     $                          cvec6,crhs,6_i4,imode)
            ELSE
              DO ibl=1,nrbl            !  diffusive, not hyperbolic
                auxb(ibl)%arri=0._r8   !  correction, so no old value
              ENDDO
              CALL fe_postsolve(b_mat(imode),cvectr,csln,
     $                          cvec6,crhs,6_i4,imode,auxb)
            ENDIF
          ENDIF
          SELECT CASE(b_pass)
c-----------------------------------------------------------------------
c         update be for this time split.
c-----------------------------------------------------------------------
          CASE ('bhmhd cor','bmhd cor','bhyp')
            DO ibl=1,nbl
              CALL cvector_assign_cvec2(cvec6(ibl),csln(ibl),imode)
            ENDDO
c-----------------------------------------------------------------------
c         trap unexpected pass flags.
c-----------------------------------------------------------------------
          CASE DEFAULT
            CALL nim_stop
     $        (TRIM(b_pass)//' pass is inappropriate for adv_b_hyp.')
          END SELECT
        ENDIF convergence
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL cvector_assignq_cvec(be(ibl),cvec6(ibl),3_i4)
        IF (nonlinear.AND.ohm_heat)
     $    CALL cvector_assignq_cvec(work1(ibl),cvec6(ibl),3_i4,
     $                              nstart1=1_i4,nstart2=4_i4)
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation and we're done with crhs.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        IF (b_pass/='bhyp') CALL extrap_update(be(ibl),ibl,b_pass)
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.
c-----------------------------------------------------------------------
      CALL regular_ave(be,3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     find heating from the hyper-dissipation.
c-----------------------------------------------------------------------
      IF (nonlinear.AND.ohm_heat) CALL hyper_bheat
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_b_hyp
c-----------------------------------------------------------------------
c     subprogram 7. adv_b_3deta.
c     advance magnetic field for cases with 3D variations in the
c     electrical diffusivity (due to temperature dependence).
c     these cases require a 3d matrix solution.
c
c     the 3D matrix is never formed; its action on a vector is found by
c     a finite element rhs computation, with the integrand specified in
c     the iter_ky_c3d_solve input.
c-----------------------------------------------------------------------
      SUBROUTINE adv_b_3deta(b_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE iter_ky_c3d_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: b_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,ibe,its,iv,ierror,imode,kmode
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed
      LOGICAL :: do_surf
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_b_3deta(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE threed_b_3deta
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_preb_3deta(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE threed_preb_3deta
      END INTERFACE
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.eta_changed).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        CALL matrix_create(bhmhd_cmat,bhmhd_cfac,curl_de_ciso,
     $                     dirichlet_comp_op,'3vn',TRIM(b_pass),
     $                     .true.,bmhd_solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     determine if n=0 tangential electric field should be applied.
c-----------------------------------------------------------------------
      IF (volt/=0.OR.e_vert/=0) THEN
        do_surf=.true.
      ELSE
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes.
c
c     cell-interior data is not eliminated for 3d solves.
c-----------------------------------------------------------------------
      CALL get_rhs(brhs_mhd,crhs,dirichlet_rhs,'3vn','cyl_vec',
     $             do_surf,e_tangential)
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve,
c     and transfer it to work1 as the initial guess.
c-----------------------------------------------------------------------
      CALL extrap_sln(b_pass,t+dt)
      DO ibl=1,nbl
        CALL extrap_eval_cvec(work1(ibl),ibl,b_pass)
        CALL vector_add(work1(ibl),be(ibl),v2fac=-1._r8)
        CALL dirichlet_rhs(work1(ibl),seam(ibl),'3vn',3_i4)
        CALL regular_vec(work1(ibl),seam(ibl),'cyl_vec',3_i4,nmodes,
     $                   nindex)
      ENDDO
c-----------------------------------------------------------------------
c     Solve the 3D equation.
c-----------------------------------------------------------------------
      CALL iter_ky_c3d_solve(threed_b_3deta,threed_preb_3deta,.true.,
     $                       work1,crhs,work4,3_i4,poly_degree,0_i4,
     $                       0_i4,nrbl,nbl,nmodes,tol,maxit,
     $                       nmodes_total,err,its,seed,old_vec=be)
      IF (node==0.AND.itflag) CALL iter_out(b_pass,seed,its,err)
      IF (node==0) mhdits=mhdits+its
c-----------------------------------------------------------------------
c     test convergence of 3D equation.
c-----------------------------------------------------------------------
      IF (err>tol) THEN
        converged=.false.
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(work1(ibe),seam(ibe),'3vn',3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components, then combine delta-B and B, and update data for
c     polynomial extrapolation.
c-----------------------------------------------------------------------
      CALL regular_ave(work1,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
        CALL vector_add(work1(ibl),be(ibl))
        CALL extrap_update(work1(ibl),ibl,b_pass)
      ENDDO
c-----------------------------------------------------------------------
c     predictor steps save an average of the new and old fields.
c     only predictor steps need quad-point updates here.
c-----------------------------------------------------------------------
      SELECT CASE(b_pass)
      CASE('bmhd pre')
        DO ibl=1,nrbl
          CALL vector_add(work1(ibl),be(ibl),fb_vxb,1._r8-fb_vxb)
          CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_add(work1(ibl),be(ibl),fb_vxb,1._r8-fb_vxb)
          CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      CASE('bmhd cor')
        DO ibl=1,nbl
          be(ibl)=work1(ibl)
        ENDDO
      CASE DEFAULT
        CALL nim_stop
     $    (TRIM(b_pass)//' pass is inappropriate for adv_b_3deta.')
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_b_3deta
c-----------------------------------------------------------------------
c     subprogram 8. adv_b_3dnsym.
c     advance magnetic field for cases with 3D variations in the
c     effective impedance matrix that is not a Hermitian operator
c     (due to implicit Hall and/or advection).
c
c     the 3D matrix is never formed; its action on a vector is found by
c     a finite element rhs computation, with the integrand specified in
c     the 3D iterative solver input.
c-----------------------------------------------------------------------
      SUBROUTINE adv_b_3dnsym(b_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE pardata
      USE iter_ky_c3d_mod
      USE mpi_nim
      USE pardata
      USE math_tran
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: b_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,ibe,its,iv,ierror,imode,kmode,ig,ng,
     $               off,hall_it,hall_mx,nqdsc,nbdsc,nqc
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err,rhsnrm,nlnrm
      CHARACTER(8) :: seed
      CHARACTER(32) :: save_flag
      LOGICAL :: do_surf
      LOGICAL, SAVE :: first_call=.true.
      COMPLEX(r8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ctmp
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: boff_cmat
      TYPE(cvector_type), DIMENSION(:), POINTER :: slptr,drptr
      TYPE(cvector_type), DIMENSION(:), POINTER, SAVE :: oldb
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_b_3dnsym(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE threed_b_3dnsym
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_preb_3dnsym(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE threed_preb_3dnsym
      END INTERFACE
c-----------------------------------------------------------------------
c     set integer values for the number of discontinuous bases and
c     the number of fields for each.
c-----------------------------------------------------------------------
      IF (poly_divb>=0.AND.nrbl>0) THEN
        nqdsc=rb(1)%auxb%nqty
        nbdsc=rb(1)%auxb%n_int
      ELSE
        nqdsc=0_i4
        nbdsc=0_i4
      ENDIF
c-----------------------------------------------------------------------
c     the structures used for solution and direction vectors differ
c     when unsplit hyper-resistivity is used.  set pointers to
c     facilitate the different possibilities.
c-----------------------------------------------------------------------
      IF ((hyp_eta>0._r8.OR.hyp_dbd>0._r8).AND..NOT.split_hypeta) THEN
        nqc=6
        slptr=>w6v1
        drptr=>w6v2
      ELSE
        nqc=3
        slptr=>work1
        drptr=>work4
      ENDIF
c-----------------------------------------------------------------------
c     passing be for the old vector of the 3d linear solve is not
c     compatible with 6-vector solves, so an extra data structure is
c     needed, despite being inefficient memory-wise.
c-----------------------------------------------------------------------
      IF (first_call) THEN
        ALLOCATE(oldb(nbl))
        DO ibl=1,nbl
          IF (ibl<=nrbl) THEN
            CALL vector_type_alloc(oldb(ibl),poly_degree,rb(ibl)%mx,
     $                             rb(ibl)%my,nqc,nmodes)
          ELSE
            CALL vector_type_alloc(oldb(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                             nqc,nmodes)
          ENDIF
        ENDDO 
        first_call=.false.
      ENDIF
    
      DO ibl=1,nbl
        oldb(ibl)=0._r8
        CALL cvector_assignq_cvec(oldb(ibl),be(ibl),3_i4)
      ENDDO 
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      IF (ohms/='mhd'.AND.maxit_nl>1) THEN
        hall_mx=maxit_nl
      ELSE
        hall_mx=1
      ENDIF
      DO ibl=1,nrbl
        IF (nqdsc>0) THEN
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,nqc,nmodes,nbdsc,nqdsc)
          slptr(ibl)%arrtmp=>rb(ibl)%mwork1%fsi  !  mwork 1 & 2 are
          drptr(ibl)%arrtmp=>rb(ibl)%mwork2%fsi  !  allocated for div(B)
          be(ibl)%arrtmp=>auxb(ibl)%arri
          work3(ibl)%arrtmp=>crhs(ibl)%arrtmp    !  see comments below
        ELSE
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,nqc,nmodes)
        ENDIF
        IF (hall_mx>1)
     $    CALL vector_type_alloc(cvecn(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,nqc,nmodes,nbdsc,nqdsc)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,nqc,
     $                         nmodes)
        IF (hall_mx>1)
     $    CALL vector_type_alloc(cvecn(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                           nqc,nmodes)
      ENDDO
      DO ibl=1,nbl
        slptr(ibl)=0._r8
        drptr(ibl)=0._r8
      ENDDO
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c
c     set all components of the auxiliary field to zero at the boundary
c     to control energy flux.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.v0_changed.OR.b0_changed.OR.n0_changed).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        CALL matrix_create(bhmhd_cmat,bhmhd_cfac,b_hmhd_op,
     $                     dirichlet_comp_op,'3vn sd sd sd','bhmhd',
     $                     .true.,bmhd_solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     determine if n=0 tangential electric field should be applied.
c-----------------------------------------------------------------------
      IF (volt/=0.OR.e_vert/=0.OR.ohms/='mhd'.AND.beta>0) THEN
        do_surf=.true.
      ELSE
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes.
c
c     cell-interior data is not eliminated for 3d solves.
c-----------------------------------------------------------------------
      IF (ohms=='mhd') THEN  !  for implicit advection
        CALL get_rhs(brhs_mhd,crhs,dirichlet_rhs,'3vn sd sd sd',
     $               'cyl_vec',do_surf,e_tangential)
      ELSE
        CALL get_rhs(brhs_hmhd,crhs,dirichlet_rhs,'3vn sd sd sd',
     $               'cyl_vec',do_surf,e_tangential)
      ENDIF
c-----------------------------------------------------------------------
c     if nonlinear iteration is used for the Hall term, copy the
c     original rhs and find the norm for the dB part (vector comps 1:3),
c     and not the div(b) part that uses the discontinuous bases.
c-----------------------------------------------------------------------
      IF (hall_mx>1) THEN
        DO ibl=1,nbl
          cvecn(ibl)=crhs(ibl)
          IF (nqdsc>0.AND.ibl<=nrbl) NULLIFY(crhs(ibl)%arrtmp)
        ENDDO
        CALL iter_c3d_err(rhsnrm,nrbl,nbl,poly_degree,crhs,'2-norm')
        IF (nqdsc>0) THEN
          DO ibl=1,nrbl
            crhs(ibl)%arrtmp=>work3(ibl)%arrtmp
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve,
c     and transfer it to work1 as the initial guess.
c
c     use the old solution as a guess for the discontinuous fields.
c-----------------------------------------------------------------------
      CALL extrap_sln(b_pass,t+dt)
      DO ibl=1,nbl
        CALL extrap_eval_cvec(work1(ibl),ibl,b_pass)
        CALL vector_add(work1(ibl),be(ibl),v2fac=-1._r8)
        CALL dirichlet_rhs(work1(ibl),seam(ibl),'3vn',3_i4)
        CALL regular_vec(work1(ibl),seam(ibl),'cyl_vec',3_i4,nmodes,
     $                   nindex)
        IF (nqdsc>0.AND.ibl<=nrbl) THEN
          slptr(ibl)%arrtmp(:,:,:,:,:)=auxb(ibl)%arri(:,:,:,:,:)
          auxb(ibl)%arri=0._r8 ! diffusive correction, so 0 old value
        ENDIF
      ENDDO
      IF (nqc>3) THEN
        DO ibl=1,nbl
          slptr(ibl)=0
          CALL cvector_assignq_cvec(slptr(ibl),work1(ibl),3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     Nonlinear iteration for the Hall term may be done here.
c-----------------------------------------------------------------------
      hall_loop: DO hall_it=1,hall_mx
c-----------------------------------------------------------------------
c       Solve the 3D equation.
c-----------------------------------------------------------------------
        CALL iter_ky_c3d_solve(threed_b_3dnsym,threed_preb_3dnsym,
     $                         .false.,slptr,crhs,drptr,nqc,poly_degree,
     $                         nbdsc,nqdsc,nrbl,nbl,nmodes,tol,maxit,
     $                         nmodes_total,err,its,seed,old_vec=oldb)
        IF (node==0.AND.itflag) CALL iter_out(b_pass,seed,its,err)
        IF (node==0) hallits=hallits+its
c-----------------------------------------------------------------------
c       test convergence of 3D equation.
c-----------------------------------------------------------------------
        IF (err>tol) THEN
          converged=.false.
          DO ibl=1,nbl
            NULLIFY(be(ibl)%arrtmp,slptr(ibl)%arrtmp,drptr(ibl)%arrtmp)
            NULLIFY(work3(ibl)%arrtmp)
          ENDDO
          RETURN
        ENDIF
c-----------------------------------------------------------------------
c       copy the dB part of the solution into work1 if this is a
c       6-vector solve.
c-----------------------------------------------------------------------
        IF (nqc>3) THEN
          DO ibl=1,nbl
            CALL cvector_assignq_cvec(work1(ibl),slptr(ibl),3_i4)
          ENDDO
        ENDIF
c-----------------------------------------------------------------------
c       if nonlinear iterations are used for the implicit Hall term,
c       compute the correction to the original rhs.  also assess the
c       norm of the change to the rhs.
c
c       this alters coefficients of the continuous bases of crhs, only,
c       so use work3 to locate the memory for the crhs discontinuous
c       bases, while its pointers are temporarily nullified.  also
c       temporarily nullify the work4 arrtmp pointer, so that nlnrm
c       is the norm for the continuous bases only.  note that the
c       vector_add routine only adds arrtmp arrays when both pointers
c       are associated.
c-----------------------------------------------------------------------
        IF (hall_it<hall_mx) THEN
          CALL regular_ave(work1,3_i4,nmodes,nindex)
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
            drptr(ibl)=crhs(ibl)
            IF (nqdsc>0) NULLIFY(crhs(ibl)%arrtmp,work4(ibl)%arrtmp)
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
            drptr(ibl)=crhs(ibl)
          ENDDO
          CALL get_rhs(hall_cor,crhs,dirichlet_rhs,'3vn 3vn','cyl_vec',
     $                 .false.,no_surf_int)
          DO ibl=1,nbl
            CALL vector_add(crhs(ibl),cvecn(ibl))
            CALL vector_add(drptr(ibl),crhs(ibl),v2fac=-1._r8)
            IF (nqc>3)
     $        CALL cvector_assignq_cvec(work4(ibl),drptr(ibl),3_i4)
          ENDDO
          CALL iter_c3d_err(nlnrm,nrbl,nbl,poly_degree,work4,'2-norm')
          IF (node==0.AND.itflag)
     $      CALL iter_out("Nlin Hall"," ",hall_it,nlnrm/rhsnrm)
          IF (nqdsc>0) THEN
            DO ibl=1,nrbl
              crhs(ibl)%arrtmp=>work3(ibl)%arrtmp
              IF (nqc==3) work4(ibl)%arrtmp=>rb(ibl)%mwork2%fsi
            ENDDO
          ENDIF
          IF (nlnrm/rhsnrm<=tol_nl) EXIT hall_loop
          CALL b_store('hall it') 
        ENDIF
      ENDDO hall_loop
c-----------------------------------------------------------------------
c     save the solution.
c-----------------------------------------------------------------------
      SELECT CASE(b_pass)
c-----------------------------------------------------------------------
c     'cor' is historical; update storage to end of step.
c-----------------------------------------------------------------------
      CASE('bhmhd cor','bmhd cor')
        DO ibl=1,nbl
          CALL vector_add(be(ibl),work1(ibl))
          CALL vector_type_dealloc(crhs(ibl))
          IF (hall_mx>1) CALL vector_type_dealloc(cvecn(ibl))
        ENDDO
c-----------------------------------------------------------------------
c     trap unexpected pass flags.
c-----------------------------------------------------------------------
      CASE DEFAULT
        CALL nim_stop
     $    (TRIM(b_pass)//' pass is inappropriate for adv_b_3dnsym.')
      END SELECT
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components, and update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      CALL regular_ave(be,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL extrap_update(be(ibl),ibl,b_pass)
        IF (nqc==6.AND.ohm_heat)
     $    CALL cvector_assignq_cvec(work1(ibl),slptr(ibl),3_i4,
     $                              nstart1=1_i4,nstart2=4_i4)
      ENDDO
c-----------------------------------------------------------------------
c     find heating from the hyper-dissipation.
c-----------------------------------------------------------------------
      IF (nqc==6.AND.ohm_heat) CALL hyper_bheat
c-----------------------------------------------------------------------
c     nullify the arrtmp pointers.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        NULLIFY(be(ibl)%arrtmp,slptr(ibl)%arrtmp,drptr(ibl)%arrtmp)
        NULLIFY(work3(ibl)%arrtmp)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_b_3dnsym
c-----------------------------------------------------------------------
c     subprogram 9. clean_divb.
c     control-routine for the divergence(B) diffusion.
c-----------------------------------------------------------------------
      SUBROUTINE clean_divb(converged)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: converged

      INTEGER(i4) :: ibl,its,i_ri,ibe,imode
      REAL(r8) :: err
      CHARACTER(8) :: seed,flag
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln('divb diff',t+dt)
c-----------------------------------------------------------------------
c     create the matrix for each Fourier component.
c-----------------------------------------------------------------------
      IF (dt/=dt_old) THEN
        CALL matrix_create(divb_mat,divb_fac,grad_div,dirichlet_op,
     $                     '3vn','divb_diff',.true.,solver,
     $                     mass_type='none')
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes, and impose
c     homogeneous dirichlet boundary conditions on the normal component
c     of the change in B.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      CALL get_rhs(divb_rhs,crhs,dirichlet_rhs,'3vn','cyl_vec',
     $             .false.,no_surf_int,rmat_elim=divb_mat)
c-----------------------------------------------------------------------
c     mode loop for the linear equation.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       for each mode treat the (real_r,real_z,-imag_phi) and
c       the (imag_r,imag_z,real_phi) components separately.
c-----------------------------------------------------------------------
        real_imag: DO i_ri=0,1 
          IF (i_ri==0) THEN
            flag='r12mi3'
          ELSE
            flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loop
            flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         evaluate the product of the matrix and the magnetic field at
c         the start of the time split.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),be(ibl),flag,imode)
          ENDDO
          CALL matvec(divb_mat(imode),vectr,sln,3_i4)
c-----------------------------------------------------------------------
c         add the product of the matrix and the old vector to the rhs,
c         then set the initial guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),crhs(ibl),flag,imode)
            IF (ASSOCIATED(sln(ibl)%arri)) sln(ibl)%arri=0
            CALL vector_add(vectr(ibl),sln(ibl))
            CALL extrap_eval_vec(sln(ibl),ibl,imode,'divb diff',flag)
          ENDDO
c-----------------------------------------------------------------------
c         invert the matrix.
c-----------------------------------------------------------------------
          CALL iter_cg_2d_solve(divb_mat(imode),divb_fac(imode),sln,
     $                          vectr,3_i4,tol,maxit,solver,err,its,
     $                          seed)
          IF (node == 0 .AND. itflag)
     $      CALL iter_out('divb diff',seed,its,err)
          dbits=its+dbits
c-----------------------------------------------------------------------
c         save solution.
c-----------------------------------------------------------------------
          IF (err>tol) THEN
            converged=.false.
            RETURN
          ELSE
c-----------------------------------------------------------------------
c           determine the interior data if eliminated for the matrix
c           solve.
c-----------------------------------------------------------------------
            IF (poly_degree>1)
     $        CALL fe_postsolve(divb_mat(imode),vectr,sln,
     $                          be,crhs,3_i4,imode,flag)
            DO ibl=1,nbl
              CALL cvector_assign_vec(be(ibl),sln(ibl),flag,imode)
            ENDDO
          ENDIF
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL extrap_update(be(ibl),ibl,'divb diff')
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      IF (zero_bnorm) THEN
        DO ibl=1,SIZE(exblock_list)
          ibe=exblock_list(ibl)
          CALL dirichlet_rhs(be(ibe),seam(ibe),'3vn',3_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.  also update storage at quadrature points.
c-----------------------------------------------------------------------
      CALL regular_ave(be,3_i4,nmodes,nindex)
      DO ibl=1,nrbl
        CALL rblock_qp_update(rb(ibl)%be,rb(ibl)%qbe,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_qp_update(tb(ibl)%be,tb(ibl)%qbe,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE clean_divb
c-----------------------------------------------------------------------
c     subprogram 10. jfromb.
c     find total current density as a bilinear quantity from magnetic
c     field.  in this version of nimrod, the bilinear j is only used
c     for the electron inertia in the 2-fluid ohm's law.  in all other
c     places it is computed locally from the curl of b.
c-----------------------------------------------------------------------
      SUBROUTINE jfromb(j_pass)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: j_pass

      INTEGER(i4) :: ibl,its,i_ri,imode,imat
      REAL(r8) :: err
      CHARACTER(8) :: seed,flag
      TYPE(global_matrix_type), POINTER :: matptr
      TYPE(matrix_factor_type), POINTER :: facptr
      LOGICAL, SAVE :: first_solve=.true.
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      IF (.NOT.first_solve.AND.j_pass=='j pre')
     $  CALL extrap_sln('jfromb',t+dt)
c-----------------------------------------------------------------------
c     evaluate J=curl(B)/mu0 in weak form.  if this is a corrector
c     step, the predicted B and its derivatives are already evaluated
c     in qwork1.
c-----------------------------------------------------------------------
      IF (j_pass=='j pre') THEN
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%be,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%be,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      ENDIF
      CALL get_rhs(curl,crhs,no_rhs_bc,' ','cyl_vec',.false.,
     $             no_surf_int)
c-----------------------------------------------------------------------
c     loop over modes.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       separate the solves for (real_r,real_z,-imag_phi) and
c       (imag_r,imag_z,real_phi), which use the same lhs.  [regularity
c       conditions lead to coupling.]  an exception is made for n=0,
c       where all real components are solved at the same time.
c-----------------------------------------------------------------------
        real_imag: DO i_ri=0,1 
          IF (i_ri==0) THEN
            flag='r12mi3'
          ELSE
            flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loop
            flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         select the mass matrix with appropriate regularity conditions
c         if necessary.  zero out r and phi comps for all but n=1;
c         zero out z for n>0.
c-----------------------------------------------------------------------
          IF (any_r0blocks) THEN
            imat=MIN(nindex(imode)+1,3)
          ELSE
            imat=1
          ENDIF
          matptr=>j_mat(imat)
          facptr=>j_fac(imat)
c-----------------------------------------------------------------------
c         copy the appropriate part of rhs for this solve, then set the
c         initial guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),crhs(ibl),flag,imode)
            CALL vector_mult(vectr(ibl),1/mu0)
            IF (.NOT.first_solve.AND.j_pass=='j pre') THEN
              CALL extrap_eval_vec(sln(ibl),ibl,imode,'jfromb',flag)
            ELSE
              CALL vector_assign_cvec(sln(ibl),ja(ibl),flag,imode)
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c         call solver.
c-----------------------------------------------------------------------
          CALL iter_cg_2d_solve(matptr,facptr,sln,vectr,3_i4,tol,
     $                          maxit,solver,err,its,seed)
          IF (node == 0 .AND. itflag)
     $      CALL iter_out('jfromb   ',seed,its,err)
          jaits=its+jaits
          IF (err>tol) THEN
            CALL nim_output(.false.)
            CALL nim_stop
     $         ('Jfromb: no convergence for mass matrix solve.')
          ELSE
            DO ibl=1,nbl
              CALL cvector_assign_vec(ja(ibl),sln(ibl),flag,imode)
            ENDDO
          ENDIF
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      IF (.NOT.first_solve.AND.j_pass=='j pre') THEN
        DO ibl=1,nbl
          CALL extrap_update(ja(ibl),ibl,'jfromb')
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.  also update storage at quadrature points.
c-----------------------------------------------------------------------
      CALL regular_ave(ja,3_i4,nmodes,nindex)
      DO ibl=1,nrbl
        CALL rblock_qp_update(rb(ibl)%ja,rb(ibl)%qja,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL tblock_qp_update(tb(ibl)%ja,tb(ibl)%qja,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
      first_solve=.false.
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE jfromb
c-----------------------------------------------------------------------
c     subprogram 11. adv_v_clap.
c     advance single-fluid velocity (momentum density/mass density)
c     with a Laplacian operator, where free-slip boundary conditions
c     are applied.  the free-slip bcs couple the poloidal components
c     and prevent splitting the r-phi components from the z component.
c     at present, it's assumed that this is only used for a semi-
c     implicit operator in v--not for viscosity.
c-----------------------------------------------------------------------
      SUBROUTINE adv_v_clap(v_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: v_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,imode,ibe,i_ri,nq
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,op_flag,bc_flag,comp_flag,precon
      LOGICAL :: new_mat,do_surf
      TYPE(global_matrix_type), DIMENSION(:), POINTER :: mat
      TYPE(matrix_factor_type), DIMENSION(:), POINTER :: fac
c-----------------------------------------------------------------------
c     interface block for surface_exb.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE surface_exb(vdum,vcent)
        USE local
        USE fields
        USE seam_storage_mod
        USE global
        USE input
        IMPLICIT NONE

        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vdum
        REAL(r8), INTENT(IN) :: vcent
        END SUBROUTINE surface_exb
      END INTERFACE
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(v_pass,t+dt)
c-----------------------------------------------------------------------
c     set the boundary condition flag.
c-----------------------------------------------------------------------
      IF (flow_bc=='free-slip') THEN
        bc_flag='3vn'
        do_surf=.true.
      ELSE
        bc_flag='all'
        do_surf=.false.
      ENDIF
c-----------------------------------------------------------------------
c     select the appropriate matrix and factor,
c     and determine if they need to be recomputed.
c-----------------------------------------------------------------------
      IF (integrand_flag(1:4)=='visc') THEN
        mat=>visc_mat
        fac=>visc_fac
        op_flag='visc'
        precon=solver
      ELSE
        mat=>vmhd_mat
        fac=>vmhd_fac
        op_flag='vmhd'
        precon=vmhd_solver
      ENDIF
      new_mat=.false.
      IF (dt/=dt_old.AND.(istep>istep_mat.OR.v_pass=='visc').OR.
     $    force_recomp)
     $  new_mat=.true.
      IF (integrand_flag(1:4)/='visc'.AND.istep>istep_mat.AND.
     $    (b0_changed.OR.p0_changed.OR.n0_changed.OR.nl_changed))
     $   new_mat=.true.
c-----------------------------------------------------------------------
c     create the separate matrices for each Fourier component with
c     storage order (r,z,+-phi), where +- indicates (real_pol,-imag_phi)
c     or (imag_pol,+real_phi).
c
c     the mass matrix is not added, since the separate term
c     proportional to mass density is added in vec_lap_op.
c-----------------------------------------------------------------------
      IF (new_mat) THEN
        CALL matrix_create(mat,fac,vec_lap_op,dirichlet_op,
     $                     bc_flag,TRIM(op_flag),.true.,
     $                     precon,mass_type='none')
        IF (v_pass/='vmhd') istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     integrate the force density over the finite elements, and impose
c     free-slip boundary conditions if needed.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      CALL get_rhs(vrhs,crhs,dirichlet_rhs,bc_flag,'cyl_vec',
     $             do_surf,mag_tension,rmat_elim=mat)
c-----------------------------------------------------------------------
c     begin the loop over modes.  at most, the (real_r,real_z,-imag_phi)
c     and (imag_r,imag_z,real_phi) components couple and each set uses
c     the same operator.  for n=0 phi is uncoupled and all real
c     components are done simultaneously.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
        real_imag: DO i_ri=0,1 
          IF (i_ri==0) THEN
            comp_flag='r12mi3'
          ELSE
            comp_flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loop
            comp_flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         evaluate the product of the matrix and the old velocity field.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),ve(ibl),comp_flag,imode)
          ENDDO
          CALL matvec(mat(imode),vectr,sln,3_i4)
c-----------------------------------------------------------------------
c         add the product of the matrix and the old vector to the rhs,
c         then set the initial guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),crhs(ibl),comp_flag,
     $                              imode)
            IF (ASSOCIATED(sln(ibl)%arri)) sln(ibl)%arri=0
            CALL vector_add(vectr(ibl),sln(ibl))
            CALL extrap_eval_vec(sln(ibl),ibl,imode,v_pass,comp_flag)
          ENDDO
c-----------------------------------------------------------------------
c         call solver.
c-----------------------------------------------------------------------
          CALL iter_cg_2d_solve(mat(imode),fac(imode),sln,vectr,3_i4,
     $                          tol,maxit,precon,err,its,seed)
          IF (node == 0 .AND. itflag)
     $      CALL iter_out(TRIM(v_pass),seed,its,err)
          IF (op_flag=='vmhd') THEN
            vmhdits=its+vmhdits
          ELSE
            viscits=its+viscits
          ENDIF
c-----------------------------------------------------------------------
c         save solution.
c-----------------------------------------------------------------------
          IF (err>tol) THEN
            converged=.false.
            RETURN
          ELSE
c-----------------------------------------------------------------------
c           determine the interior data if eliminated for the matrix
c           solve.
c-----------------------------------------------------------------------
            IF (poly_degree>1)
     $        CALL fe_postsolve(mat(imode),vectr,sln,
     $                          ve,crhs,3_i4,imode,comp_flag)
            udpate_loop: DO ibl=1,nbl
              SELECT CASE(v_pass)
c-----------------------------------------------------------------------
c             put the fv_vdgv-centered prediction in work1.
c-----------------------------------------------------------------------
              CASE('v pre')
                CALL vector_assign_cvec(vectr(ibl),ve(ibl),comp_flag,
     $                                  imode)
                CALL vector_add(sln(ibl),vectr(ibl),v1fac=fv_vdgv,
     $                          v2fac=1-fv_vdgv)
                CALL cvector_assign_vec(work1(ibl),sln(ibl),comp_flag,
     $                                  imode)
c-----------------------------------------------------------------------
c             update ve.
c-----------------------------------------------------------------------
              CASE('v cor','vmhd','visc')
                CALL cvector_assign_vec(ve(ibl),sln(ibl),comp_flag,
     $                                  imode)
c-----------------------------------------------------------------------
c             trap unexpected pass flags.
c-----------------------------------------------------------------------
              CASE DEFAULT
                CALL nim_stop
     $            (TRIM(v_pass)//' pass is inappropriate for'//
     $             ' adv_v_clap.')
              END SELECT
            ENDDO udpate_loop
          ENDIF
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c     we're done with crhs and it has the same nq as work1, so it is
c     used for temporary storage here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        SELECT CASE(v_pass)
        CASE('v pre')
          crhs(ibl)=work1(ibl)
          CALL vector_add(crhs(ibl),ve(ibl),1/fv_vdgv,1-1/fv_vdgv)
          CALL extrap_update(crhs(ibl),ibl,v_pass)
        CASE('v cor','vmhd','visc')
          CALL extrap_update(ve(ibl),ibl,v_pass)
        END SELECT
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to
c     avoid accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        SELECT CASE(v_pass)
        CASE('v pre')
          CALL dirichlet_rhs(work1(ibe),seam(ibe),bc_flag,3_i4)
        CASE('v cor','vmhd','visc')
          CALL dirichlet_rhs(ve(ibe),seam(ibe),bc_flag,3_i4)
        END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     if a surface electric field is applied, set the normal component
c     of flow to exb/b**2 for the n=0 mode.
c-----------------------------------------------------------------------
      IF ((loop_volt/=0.OR.i_desired/=0.OR.e_vertical/=0).AND.
     $    norm_flow(1:3)=='exb') THEN
        SELECT CASE(v_pass)
        CASE('v pre')
          CALL surface_exb(work1,fv_vdgv)
        CASE('v cor','vmhd','visc')
          CALL surface_exb(ve,1._r8)
        END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.  also, update temporary storage at quad points.
c-----------------------------------------------------------------------
      SELECT CASE(v_pass)
      CASE('v pre')
        CALL regular_ave(work1,3_i4,nmodes,nindex)
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      CASE('v cor','vmhd','visc')
        CALL regular_ave(ve,3_i4,nmodes,nindex)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_v_clap
c-----------------------------------------------------------------------
c     subprogram 12. adv_v_aniso.
c     control-routine for a velocity advance with an anisotropic
c     semi-implicit operator.  for a predictor step, the result is
c     averaged with the old v and saved in work1.  for a corrector
c     step, the result is placed in ve.
c-----------------------------------------------------------------------
      SUBROUTINE adv_v_aniso(v_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE iter_gmres_c2d
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: v_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,iv,its,ibe,mv,imode,nvdis
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,bc_flag
      LOGICAL :: new_mat,do_surf
c-----------------------------------------------------------------------
c     interface block for surface_exb.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE surface_exb(vdum,vcent)
        USE local
        USE fields
        USE seam_storage_mod
        USE global
        USE input
        IMPLICIT NONE

        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vdum
        REAL(r8), INTENT(IN) :: vcent
        END SUBROUTINE surface_exb
      END INTERFACE
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(csln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        IF (poly_divv>=0) THEN
          nvdis=rb(ibl)%auxv%n_int
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes,nvdis,2_i4)
        ELSE
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes)
        ENDIF
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(v_pass,t+dt)
c-----------------------------------------------------------------------
c     create the matrix for each Fourier component including
c     homogeneous dirichlet boundary conditions on the normal component
c     for inviscid conditions or on all components for viscous cases.
c
c     the mass matrix is not added, since the separate term
c     proportional to mass density is added in v_aniso_op.
c-----------------------------------------------------------------------
      IF (flow_bc=='free-slip') THEN
        bc_flag='3vn'
        do_surf=.true.
      ELSE
        bc_flag='all'
        do_surf=.false.
      ENDIF
      new_mat=.false.
      IF ( force_recomp .OR. istep>istep_mat .AND. (dt/=dt_old.OR.
     $     b0_changed.OR.p0_changed.OR.n0_changed.OR.nl_changed.OR.
     $     impladv.AND.v0_changed) ) THEN
        CALL matrix_create(vmhd_cmat,vmhd_cfac,v_aniso_op,
     $                     dirichlet_comp_op,bc_flag,'vmhd',
     $                     .true.,vmhd_solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes.
c
c     cell-interior data is also eliminated for poly_degree>1 to reduce
c     the matrix equations solved below.
c-----------------------------------------------------------------------
      CALL get_rhs(vrhs,crhs,dirichlet_rhs,bc_flag,'cyl_vec',
     $             do_surf,mag_tension,cmat_elim=vmhd_cmat)
c-----------------------------------------------------------------------
c     loop over the modes, which are independent in this version.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       evaluate the product of the matrix and the old velocity.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),ve(ibl),imode)
        ENDDO
        CALL matvec(vmhd_cmat(imode),cvectr,csln,3_i4)
c-----------------------------------------------------------------------
c       add the product of the matrix and the old vector to the rhs,
c       then set the initial guess.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
          IF (ASSOCIATED(csln(ibl)%arri)) csln(ibl)%arri=0
          CALL vector_add(cvectr(ibl),csln(ibl))
          CALL extrap_eval_cvec2D(csln(ibl),ibl,imode,v_pass)
        ENDDO
c-----------------------------------------------------------------------
c       invert the matrix.
c-----------------------------------------------------------------------
        IF (impladv) THEN
          CALL iter_gmr_c2d_solve(vmhd_cmat(imode),vmhd_cfac(imode),
     $                            csln,cvectr,3_i4,tol,maxit,
     $                            vmhd_solver,err,its,seed)
        ELSE
          CALL iter_cg_2d_solve(vmhd_cmat(imode),vmhd_cfac(imode),csln,
     $                          cvectr,3_i4,tol,maxit,vmhd_solver,err,
     $                          its,seed)
        ENDIF
        IF (node == 0 .AND. itflag) CALL iter_out(v_pass,seed,its,err)
        vmhdits=its+vmhdits
c-----------------------------------------------------------------------
c       save the solution.
c-----------------------------------------------------------------------
        convergence: IF (err>tol) THEN
          converged=.false.
          RETURN
        ELSE convergence
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.  this will include the auxiliary fields for stabilizing
c         perpendicular divergence and parallel vorticity if they are
c         used.
c-----------------------------------------------------------------------
          IF (poly_degree>1) THEN
            IF (poly_divv<0) THEN
              CALL fe_postsolve(vmhd_cmat(imode),cvectr,csln,
     $                          ve,crhs,3_i4,imode)
            ELSE
              DO ibl=1,nrbl  !  diffusive correction, so no old value
                auxv(ibl)%arri(:,:,:,:,:)=0._r8
              ENDDO
              CALL fe_postsolve(vmhd_cmat(imode),cvectr,csln,
     $                          ve,crhs,3_i4,imode,auxv)
            ENDIF
          ENDIF
          SELECT CASE(v_pass)
c-----------------------------------------------------------------------
c         put the fv_vdgv-centered V in work1 for the mhd corrector
c         step.  this is only used for predictor/corrector advection.
c-----------------------------------------------------------------------
          CASE ('v pre')
            DO ibl=1,nbl
              CALL cvector_2D_assign_cvec(cvectr(ibl),ve(ibl),imode)
              CALL vector_add(csln(ibl),cvectr(ibl),v1fac=fv_vdgv,
     $                        v2fac=1-fv_vdgv)
              CALL cvector_assign_cvec2(work1(ibl),csln(ibl),imode)
            ENDDO
c-----------------------------------------------------------------------
c         update ve for this time split.
c-----------------------------------------------------------------------
          CASE ('v cor','vmhd')
            DO ibl=1,nbl
              CALL cvector_assign_cvec2(ve(ibl),csln(ibl),imode)
            ENDDO
c-----------------------------------------------------------------------
c         trap unexpected pass flags.
c-----------------------------------------------------------------------
          CASE DEFAULT
            CALL nim_stop
     $        (TRIM(v_pass)//' pass is inappropriate for adv_v_aniso.')
          END SELECT
        ENDIF convergence
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c     we're done with crhs and it has the same nq as work1, so it is
c     used for temporary storage here.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        SELECT CASE(v_pass)
        CASE('v pre')
          crhs(ibl)=work1(ibl)
          CALL vector_add(crhs(ibl),ve(ibl),1/fv_vdgv,1-1/fv_vdgv)
          CALL extrap_update(crhs(ibl),ibl,v_pass)
        CASE('v cor','vmhd')
          CALL extrap_update(ve(ibl),ibl,v_pass)
        END SELECT
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        SELECT CASE(v_pass)
        CASE('v pre')
          CALL dirichlet_rhs(work1(ibe),seam(ibe),bc_flag,3_i4)
        CASE('v cor','vmhd')
          CALL dirichlet_rhs(ve(ibe),seam(ibe),bc_flag,3_i4)
        END SELECT
      ENDDO
c-----------------------------------------------------------------------
c     if a surface electric field is applied, set the normal component
c     of flow to exb/b**2 for the n=0 mode.
c-----------------------------------------------------------------------
      IF ((loop_volt/=0.OR.i_desired/=0.OR.e_vertical/=0).AND.
     $    norm_flow(1:3)=='exb') THEN
        SELECT CASE(v_pass)
        CASE('v pre')
          CALL surface_exb(work1,fv_vdgv)
        CASE('v cor','vmhd','visc')
          CALL surface_exb(ve,1._r8)
        END SELECT
      ENDIF
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.  also update temporary storage at quad points.
c-----------------------------------------------------------------------
      SELECT CASE(v_pass)
      CASE('v pre')
        CALL regular_ave(work1,3_i4,nmodes,nindex)
        DO ibl=1,nrbl
          CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      CASE('v cor','vmhd')
        CALL regular_ave(ve,3_i4,nmodes,nindex)
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_v_aniso
c-----------------------------------------------------------------------
c     subprogram 13. adv_v_3dn.
c     advance velocity for cases with 3D variations in the mass density
c     or nonlinear computations with anisotropic viscosity.
c     in either case, Fourier coupling appears on the lhs of the
c     equation, these cases require a 3d matrix solution.
c
c     the 3D matrix is never formed; its action on a vector is found by
c     a finite element rhs computation, with the integrand specified in
c     the iter_ky_c3d_solve input.
c-----------------------------------------------------------------------
      SUBROUTINE adv_v_3dn(v_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE pardata
      USE iter_ky_c3d_mod
      USE mpi_nim
      USE pardata
      USE math_tran
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: v_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,ibe,its,iv,ierror,imode,kmode,ig,ng,
     $               adv_it,adv_mx,nqdsc,nbdsc
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err,rhsnrm,nlnrm
      CHARACTER(8) :: seed,bc_flag
      LOGICAL :: do_surf
c-----------------------------------------------------------------------
c     interface block for surface_exb.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE surface_exb(vdum,vcent)
        USE local
        USE fields
        USE seam_storage_mod
        USE global
        USE input
        IMPLICIT NONE

        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: vdum
        REAL(r8), INTENT(IN) :: vcent
        END SUBROUTINE surface_exb
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_v_3dn(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE threed_v_3dn
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_prev_3dn(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE threed_prev_3dn
      END INTERFACE
c-----------------------------------------------------------------------
c     set integer values for the number of discontinuous bases and
c     the number of fields for each.
c-----------------------------------------------------------------------
      IF (poly_divv>=0.AND.nrbl>0) THEN
        nqdsc=rb(1)%auxv%nqty
        nbdsc=rb(1)%auxv%n_int
      ELSE
        nqdsc=0_i4
        nbdsc=0_i4
      ENDIF
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      IF (impladv.AND.maxit_nl>1.AND.
     $    (advect=='V only'.OR.advect=='all'))
     $  THEN
        adv_mx=maxit_nl
      ELSE
        adv_mx=1
      ENDIF
      DO ibl=1,nrbl
        IF (nqdsc>0) THEN
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes,nbdsc,nqdsc)
          work1(ibl)%arrtmp=>rb(ibl)%mwork3%fsi  !  mwork 3 & 4 are
          work4(ibl)%arrtmp=>rb(ibl)%mwork4%fsi  !  allocated for V
          ve(ibl)%arrtmp=>auxv(ibl)%arri
          work5(ibl)%arrtmp=>crhs(ibl)%arrtmp    ! see comments below
        ELSE
          CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes)
        ENDIF
        IF (adv_mx>1)
     $    CALL vector_type_alloc(cvecn(ibl),poly_degree,rb(ibl)%mx,
     $                           rb(ibl)%my,3_i4,nmodes,nbdsc,nqdsc)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
        IF (adv_mx>1)
     $    CALL vector_type_alloc(cvecn(ibl),1_i4,tb(ibl)%mvert,0_i4,
     $                           3_i4,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c-----------------------------------------------------------------------
      IF (flow_bc=='free-slip') THEN
        bc_flag='3vn'
        do_surf=.true.
      ELSE
        bc_flag='all'
        do_surf=.false.
      ENDIF

      IF ( force_recomp .OR. istep>istep_mat .AND. (dt/=dt_old.OR.
     $     b0_changed.OR.p0_changed.OR.n0_changed.OR.nl_changed.OR.
     $     v0_changed) ) THEN
        CALL matrix_create(vmhd_cmat,vmhd_cfac,v_aniso_op,
     $                     dirichlet_comp_op,bc_flag,'vmhd',
     $                     .true.,vmhd_solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes.
c
c     cell-interior data is not eliminated for 3d solves.
c-----------------------------------------------------------------------
      CALL get_rhs(vrhs,crhs,dirichlet_rhs,bc_flag,'cyl_vec',
     $             do_surf,mag_tension)
c-----------------------------------------------------------------------
c     if nonlinear iteration is used for advection, copy the
c     original rhs and find the norm of the dV part, i.e. the continuous
c     bases only.
c-----------------------------------------------------------------------
      IF (adv_mx>1) THEN
        DO ibl=1,nbl
          cvecn(ibl)=crhs(ibl)
          IF (nqdsc>0.AND.ibl<=nrbl) NULLIFY(crhs(ibl)%arrtmp)
        ENDDO
        CALL iter_c3d_err(rhsnrm,nrbl,nbl,poly_degree,crhs,'2-norm')
        IF (nqdsc>0) THEN
          DO ibl=1,nrbl
            crhs(ibl)%arrtmp=>work5(ibl)%arrtmp
          ENDDO
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve,
c     and transfer it to work1 as the initial guess.
c
c     use the old solution as a guess for the discontinuous fields.
c-----------------------------------------------------------------------
      CALL extrap_sln(v_pass,t+dt)
      DO ibl=1,nbl
        CALL extrap_eval_cvec(work1(ibl),ibl,v_pass)
        CALL vector_add(work1(ibl),ve(ibl),v2fac=-1._r8)
        CALL dirichlet_rhs(work1(ibl),seam(ibl),bc_flag,3_i4)
        CALL regular_vec(work1(ibl),seam(ibl),'cyl_vec',3_i4,nmodes,
     $                   nindex)
        IF (nqdsc>0.AND.ibl<=nrbl) THEN
          work1(ibl)%arrtmp(:,:,:,:,:)=auxv(ibl)%arri(:,:,:,:,:)
          auxv(ibl)%arri=0._r8 ! diffusive correction, so 0 old value
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     nonlinear iteration for advection may be done here.
c-----------------------------------------------------------------------
      adv_loop: DO adv_it=1,adv_mx
c-----------------------------------------------------------------------
c       Solve the 3D equation.
c-----------------------------------------------------------------------
        CALL iter_ky_c3d_solve(threed_v_3dn,threed_prev_3dn,
     $                         vmhd_cmat(1)%hermitian,work1,crhs,work4,
     $                         3_i4,poly_degree,nbdsc,nqdsc,nrbl,nbl,
     $                         nmodes,tol,maxit,nmodes_total,
     $                         err,its,seed,old_vec=ve)
        IF (node==0.AND.itflag) CALL iter_out(v_pass,seed,its,err)
        IF (node==0) vmhdits=vmhdits+its
c-----------------------------------------------------------------------
c       test convergence of 3D equation.
c-----------------------------------------------------------------------
        IF (err>tol) THEN
          converged=.false.
          DO ibl=1,nbl
            NULLIFY(ve(ibl)%arrtmp,work1(ibl)%arrtmp,work4(ibl)%arrtmp)
            NULLIFY(work5(ibl)%arrtmp)
          ENDDO
          RETURN
        ENDIF
c-----------------------------------------------------------------------
c       if nonlinear iterations are used for implicit advection,
c       compute the correction to the original rhs.  also assess the
c       norm of the change to the rhs.
c
c       this alters coefficients of the continuous bases of crhs, only,
c       so use work5 to locate the memory for the crhs discontinuous
c       bases, while its pointers are temporarily nullified.  also
c       temporarily nullify the work4 arrtmp pointer, so that nlnrm
c       is the norm for the continuous bases only.  note that the
c       vector_add routine only adds arrtmp arrays when both pointers
c       are associated.
c-----------------------------------------------------------------------
        IF (adv_it<adv_mx) THEN
          CALL regular_ave(work1,3_i4,nmodes,nindex)
          DO ibl=1,nrbl
            CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
            work4(ibl)=crhs(ibl)
            IF (nqdsc>0) NULLIFY(crhs(ibl)%arrtmp,work4(ibl)%arrtmp)
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
            work4(ibl)=crhs(ibl)
          ENDDO
          CALL get_rhs(advect_cor,crhs,dirichlet_rhs,bc_flag,'cyl_vec',
     $                 .false.,no_surf_int)
          DO ibl=1,nbl
            CALL vector_add(crhs(ibl),cvecn(ibl))
            CALL vector_add(work4(ibl),crhs(ibl),v2fac=-1._r8)
          ENDDO
          CALL iter_c3d_err(nlnrm,nrbl,nbl,poly_degree,work4,'2-norm')
          IF (node==0.AND.itflag)
     $      CALL iter_out("Nlin V.gV"," ",adv_it,nlnrm/rhsnrm)
          IF (nqdsc>0) THEN
            DO ibl=1,nrbl
              crhs(ibl)%arrtmp=>work5(ibl)%arrtmp
              work4(ibl)%arrtmp=>rb(ibl)%mwork4%fsi
            ENDDO
          ENDIF
          IF (nlnrm/rhsnrm<=tol_nl) EXIT adv_loop
          CALL vcom_store('advect it') 
        ENDIF
      ENDDO adv_loop
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to avoid
c     accumulating finite cg-solver tolerance errors.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(work1(ibe),seam(ibe),bc_flag,3_i4)
      ENDDO
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.
c-----------------------------------------------------------------------
      CALL regular_ave(work1,3_i4,nmodes,nindex)
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
        IF (adv_mx>1) CALL vector_type_dealloc(cvecn(ibl))
        CALL vector_add(work1(ibl),ve(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     if a surface electric field is applied, set the normal component
c     of flow to exb/b**2 for the n=0 mode.
c-----------------------------------------------------------------------
      IF ((loop_volt/=0.OR.i_desired/=0.OR.e_vertical/=0).AND.
     $    norm_flow(1:3)=='exb') CALL surface_exb(work1,1._r8)
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL extrap_update(work1(ibl),ibl,v_pass)
      ENDDO
c-----------------------------------------------------------------------
c     predictor steps save an average of the new and old fields.
c     only predictor steps need quad-point updates here.
c-----------------------------------------------------------------------
      SELECT CASE(v_pass)
      CASE('v pre')
        DO ibl=1,nrbl
          CALL vector_add(work1(ibl),ve(ibl),fv_vdgv,1._r8-fv_vdgv)
          CALL rblock_qp_update(rb(ibl)%work1,rb(ibl)%qwork1,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_add(work1(ibl),ve(ibl),fv_vdgv,1._r8-fv_vdgv)
          CALL tblock_qp_update(tb(ibl)%work1,tb(ibl)%qwork1,tb(ibl))
        ENDDO
      CASE('v cor','vmhd')
        DO ibl=1,nbl
          ve(ibl)=work1(ibl)
        ENDDO
      CASE DEFAULT
        CALL nim_stop
     $    (TRIM(v_pass)//' pass is inappropriate for adv_v_3dn.')
      END SELECT
c-----------------------------------------------------------------------
c     nullify the arrtmp pointers.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        NULLIFY(ve(ibl)%arrtmp,work1(ibl)%arrtmp,work4(ibl)%arrtmp)
        NULLIFY(work5(ibl)%arrtmp)
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_v_3dn
c-----------------------------------------------------------------------
c     subprogram 14. adv_t_aniso.
c     temperature advance with a 3D operator (all nonlinear cases).  
c     for a predictor step, the result is averaged with the old 
c     temperature and saved in work3.
c
c     this routine has been modified to solve a 3D matrix equation
c     using iter_ky_c3d_solve.  the 3D matrix is never formed; its
c     action on a vector is found by a finite element rhs computation,
c     with the integrand specified in the iter_ky_c3d_solve input.
c-----------------------------------------------------------------------
      SUBROUTINE adv_t_aniso(t_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE pardata
      USE iter_ky_c3d_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: t_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,iv,ierror,imode,kmode,off
      INTEGER(i4), SAVE :: istep_ti_mat=-1
      INTEGER(i4), SAVE :: istep_te_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed
      CHARACTER(32) :: save_flag
      TYPE(cvector_type), DIMENSION(:), POINTER :: temp
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: t_mat
      TYPE(complex_factor_type), DIMENSION(:), POINTER :: t_fac
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: toff_cmat
      LOGICAL :: dirbc
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_t_aniso(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE threed_t_aniso
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_pret_aniso(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE threed_pret_aniso
      END INTERFACE
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     select Te or Ti.
c-----------------------------------------------------------------------
      IF (t_pass(1:2)=='te') THEN
        temp=>tele
        t_mat=>te_cmat
        t_fac=>te_cfac
      ELSE
        temp=>tion
        t_mat=>ti_cmat
        t_fac=>ti_cfac
      ENDIF
      IF (p_model=='adiabat'.OR.p_model=='isothermal'.OR.insulate) THEN
        dirbc=.false.
      ELSE
        dirbc=.true.
      ENDIF
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.b0_changed.OR.kpll_changed.OR.n0_changed.AND.
     $     (continuity=='n=0 only'.OR.continuity=='full').OR.
     $     v0_changed).AND.
     $    (t_pass(1:2)=='ti'.AND.istep>istep_ti_mat.OR.
     $     t_pass(1:2)=='te'.AND.istep>istep_te_mat.OR.kpll_changed).OR.
     $     force_recomp) THEN
        IF (t_pass(1:2)=='ti') THEN
          istep_ti_mat=istep
        ELSE
          istep_te_mat=istep
        ENDIF
        IF (dirbc) THEN
          CALL matrix_create(t_mat,t_fac,t_aniso_op,dirichlet_comp_op,
     $                       'sd','temperature',.true.,temp_solver,
     $                       mass_type='none')
        ELSE
          CALL matrix_create(t_mat,t_fac,t_aniso_op,no_comp_mat_bc,
     $                       ' ','temperature',.true.,temp_solver,
     $                       mass_type='none')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the rhs for the temperature change.
c-----------------------------------------------------------------------
      IF (t_pass(1:2)=='te') THEN
        IF (dirbc) THEN
          CALL get_rhs(terhs,crhs,dirichlet_rhs,'sd','scalar',
     $                 .false.,no_surf_int)
        ELSE
          CALL get_rhs(terhs,crhs,no_rhs_bc,' ','scalar',
     $                 .false.,no_surf_int)
        ENDIF
      ELSE
        IF (dirbc) THEN
          CALL get_rhs(tirhs,crhs,dirichlet_rhs,'sd','scalar',
     $                 .false.,no_surf_int)
        ELSE
          CALL get_rhs(tirhs,crhs,no_rhs_bc,' ','scalar',
     $                 .false.,no_surf_int)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve,
c     and transfer it to work3 as the initial guess.
c
c     the nonlinear correction steps are special cases that use the
c     difference between the first solve and the old solution as the
c     guess, but the old fields need to be reset before the solve.
c-----------------------------------------------------------------------
      IF (t_pass(8:11)/='kart') THEN
        CALL extrap_sln(t_pass,t+dt)
        DO ibl=1,nbl
          CALL extrap_eval_cvec(work3(ibl),ibl,t_pass)
          CALL vector_add(work3(ibl),temp(ibl),v2fac=-1._r8)
          IF (dirbc) CALL dirichlet_rhs(work3(ibl),seam(ibl),'sd',1_i4)
        ENDDO
      ELSE IF (t_pass=='ti cor kart') THEN
        DO ibl=1,nbl
          work3(ibl)=tion(ibl)
          tion(ibl)=ti_old(ibl)
          CALL vector_add(work3(ibl),ti_old(ibl),v2fac=-1._r8)
        ENDDO
      ELSE IF (t_pass=='te cor kart') THEN
        DO ibl=1,nbl
          work3(ibl)=tele(ibl)
          tele(ibl)=te_old(ibl)
          CALL vector_add(work3(ibl),te_old(ibl),v2fac=-1._r8)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     Solve the 3D equation.
c-----------------------------------------------------------------------
      CALL iter_ky_c3d_solve(threed_t_aniso,threed_pret_aniso,
     $                       t_mat(1)%hermitian,work3,crhs,work2,1_i4,
     $                       poly_degree,0_i4,0_i4,nrbl,nbl,nmodes,tol,
     $                       maxit,nmodes_total,err,its,seed,
     $                       old_vec=temp)

      IF (node==0.AND.itflag) CALL iter_out(t_pass,seed,its,err)
      IF (node==0) THEN
        IF (t_pass(1:2)=='te') THEN
          teits=teits+its
        ELSE
          tiits=tiits+its
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     test convergence of 3D equation.
c-----------------------------------------------------------------------
      IF (err>tol) THEN
        converged=.false.
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     combine delta-T and T, then update data for polynomial
c     extrapolation.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(crhs(ibl))
        CALL vector_add(work3(ibl),temp(ibl))
        IF (t_pass(8:11)/='kart')
     $    CALL extrap_update(work3(ibl),ibl,t_pass)
      ENDDO
c-----------------------------------------------------------------------
c     predictor steps save an average of the new and old fields.
c     only predictor steps need quad-point updates here.
c-----------------------------------------------------------------------
      SELECT CASE(t_pass(4:6))
      CASE('pre')
        DO ibl=1,nrbl
          CALL vector_add(work3(ibl),temp(ibl),fp_vdgp,1._r8-fp_vdgp)
          CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qwork3,rb(ibl))
        ENDDO
        DO ibl=nrbl+1,nbl
          CALL vector_add(work3(ibl),temp(ibl),fp_vdgp,1._r8-fp_vdgp)
          CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qwork3,tb(ibl))
        ENDDO
      CASE('cor')
        DO ibl=1,nbl
          temp(ibl)=work3(ibl)
        ENDDO
      CASE DEFAULT
        CALL nim_stop
     $    (TRIM(t_pass)//' pass is inappropriate for adv_t_aniso.')
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_t_aniso
c-----------------------------------------------------------------------
c     subprogram 15. adv_t_nsym.
c     temperature advance with a possibly non-Hermiatian 2D operator. 
c-----------------------------------------------------------------------
      SUBROUTINE adv_t_nsym(t_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE pardata
      USE iter_cg
      USE iter_gmres_c2d
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: t_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,iv,ierror,imode,kmode
      INTEGER(i4), SAVE :: istep_ti_mat=-1
      INTEGER(i4), SAVE :: istep_te_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed
      TYPE(cvector_type), DIMENSION(:), POINTER :: temp,temp_old
      TYPE(complex_matrix_type), DIMENSION(:), POINTER :: t_mat
      TYPE(complex_factor_type), DIMENSION(:), POINTER :: t_fac
      LOGICAL :: dirbc
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(csln(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,1_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     select Te or Ti.
c-----------------------------------------------------------------------
      IF (t_pass(1:2)=='te') THEN
        temp=>tele
        temp_old=>te_old
        t_mat=>te_cmat
        t_fac=>te_cfac
      ELSE
        temp=>tion
        temp_old=>ti_old
        t_mat=>ti_cmat
        t_fac=>ti_cfac
      ENDIF
      IF (p_model=='adiabat'.OR.p_model=='isothermal'.OR.insulate) THEN
        dirbc=.false.
      ELSE
        dirbc=.true.
      ENDIF
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.b0_changed.OR.kpll_changed.OR.n0_changed.AND.
     $     (continuity=='n=0 only'.OR.continuity=='full').OR.
     $     v0_changed).AND.
     $    (t_pass(1:2)=='ti'.AND.istep>istep_ti_mat.OR.
     $     t_pass(1:2)=='te'.AND.istep>istep_te_mat).OR.
     $     force_recomp) THEN
        IF (t_pass(1:2)=='ti') THEN
          istep_ti_mat=istep
        ELSE
          istep_te_mat=istep
        ENDIF
        IF (dirbc) THEN
          CALL matrix_create(t_mat,t_fac,t_aniso_op,dirichlet_comp_op,
     $                       'sd','temperature',.true.,temp_solver,
     $                       mass_type='none')
        ELSE
          CALL matrix_create(t_mat,t_fac,t_aniso_op,no_comp_mat_bc,
     $                       ' ','temperature',.true.,temp_solver,
     $                       mass_type='none')
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the rhs for the temperature change.
c-----------------------------------------------------------------------
      IF (t_pass(1:2)=='te') THEN
        IF (dirbc) THEN
          CALL get_rhs(terhs,crhs,dirichlet_rhs,'sd','scalar',
     $                 .false.,no_surf_int,cmat_elim=te_cmat)
        ELSE
          CALL get_rhs(terhs,crhs,no_rhs_bc,' ','scalar',
     $                 .false.,no_surf_int,cmat_elim=te_cmat)
        ENDIF
      ELSE
        IF (dirbc) THEN
          CALL get_rhs(tirhs,crhs,dirichlet_rhs,'sd','scalar',
     $                 .false.,no_surf_int,cmat_elim=ti_cmat)
        ELSE
          CALL get_rhs(tirhs,crhs,no_rhs_bc,' ','scalar',
     $                 .false.,no_surf_int,cmat_elim=ti_cmat)
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c
c     if this is a recompute step, use the difference between the
c     first result and the old temperature as the guess, then over-write
c     the first results with the old temperature.
c-----------------------------------------------------------------------
      IF (t_pass(8:10)/='rec') THEN
        CALL extrap_sln(t_pass,t+dt)
      ELSE
        DO ibl=1,nbl
          work3(ibl)=temp(ibl)
          CALL vector_add(work3(ibl),temp_old(ibl),v2fac=-1._r8)
          temp(ibl)=temp_old(ibl)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     loop over the modes, which are independent in this version.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       evaluate the product of the matrix and the old temperature.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),temp(ibl),imode)
        ENDDO
        CALL matvec(t_mat(imode),cvectr,csln,1_i4)
c-----------------------------------------------------------------------
c       add the product of the matrix and the old vector to the rhs,
c       then set the initial guess.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
          IF (ASSOCIATED(csln(ibl)%arri).AND.t_mat(imode)%eliminated)
     $      csln(ibl)%arri=0
          CALL vector_add(cvectr(ibl),csln(ibl))
          IF (t_pass(8:10)/='rec') THEN
            CALL extrap_eval_cvec2D(csln(ibl),ibl,imode,t_pass)
          ELSE
            CALL cvector_2D_assign_cvec(csln(ibl),work3(ibl),imode)
          ENDIF
        ENDDO
c-----------------------------------------------------------------------
c       invert the matrix.
c-----------------------------------------------------------------------
        IF (t_mat(imode)%hermitian) THEN
          CALL iter_cg_2d_solve(t_mat(imode),t_fac(imode),csln,cvectr,
     $                          1_i4,tol,maxit,temp_solver,err,its,seed)
        ELSE
          CALL iter_gmr_c2d_solve(t_mat(imode),t_fac(imode),csln,cvectr,
     $                            1_i4,tol,maxit,temp_solver,err,its,
     $                            seed)
        ENDIF
        IF (node==0.AND.itflag) CALL iter_out(t_pass,seed,its,err)
        IF (t_pass(1:2)=='te') THEN
          teits=teits+its
        ELSE
          tiits=tiits+its
        ENDIF
c-----------------------------------------------------------------------
c       convergence result.
c-----------------------------------------------------------------------
        convergence: IF (err>tol) THEN
          converged=.false.
          RETURN
        ELSE convergence
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.
c-----------------------------------------------------------------------
          IF (poly_degree>1.AND.t_mat(imode)%eliminated)
     $      CALL fe_postsolve(t_mat(imode),cvectr,csln,
     $                        temp,crhs,1_i4,imode)
          update: DO ibl=1,nbl
            SELECT CASE(t_pass(4:6))
c-----------------------------------------------------------------------
c           update T for this time split.
c-----------------------------------------------------------------------
            CASE ('cor')
              CALL cvector_assign_cvec2(temp(ibl),csln(ibl),imode)
c-----------------------------------------------------------------------
c           trap unexpected pass flags.
c-----------------------------------------------------------------------
            CASE DEFAULT
              CALL nim_stop
     $          (TRIM(t_pass)//' pass is inappropriate for adv_t_nsym.')
            END SELECT
          ENDDO update
        ENDIF convergence
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      IF (t_pass(8:10)/='rec') THEN
        DO ibl=1,nbl
          CALL extrap_update(temp(ibl),ibl,t_pass)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_t_nsym
c-----------------------------------------------------------------------
c     subprogram 16. adv_nd.
c     number density advance.  for a predictor step, the result is
c     averaged with the old density and saved in work3.  for a
c     corrector step, the result is placed in nd. 
c-----------------------------------------------------------------------
      SUBROUTINE adv_nd(n_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: n_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,i_ri,iv,imode
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,flag
      TYPE(global_matrix_type), POINTER :: n_mat
      TYPE(matrix_factor_type), POINTER :: n_fac
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,1_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,1_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4,
     $                         nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(n_pass,t+dt)
c-----------------------------------------------------------------------
c     create the particle density diffusivity operator.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.n0_changed.AND.nd_floor>0).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        IF (nd_bc=='dirichlet'.AND.nd_diff>0) THEN
          CALL matrix_create(ndiso_mat,ndiso_fac,n_iso_op,dirichlet_op,
     $                       'sd','number density',.true.,solver,
     $                       mass_type='none')
        ELSE
          CALL matrix_create(ndiso_mat,ndiso_fac,n_iso_op,no_mat_bc,
     $                       ' ','number density',.true.,solver,
     $                       mass_type='none')
        ENDIF
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side.
c-----------------------------------------------------------------------
      IF (nd_diff>0.AND.nd_bc=='dirichlet') THEN
        CALL get_rhs(ndrhs,crhs,dirichlet_rhs,'sd','scalar',.false.,
     $               no_surf_int,rmat_elim=ndiso_mat)
      ELSE
        CALL get_rhs(ndrhs,crhs,no_rhs_bc,' ','scalar',.false.,
     $               no_surf_int,rmat_elim=ndiso_mat)
      ENDIF
c-----------------------------------------------------------------------
c     Loop over modes.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       with the generalized bases, the continuity equation no longer
c       relies on the mass matrix.
c-----------------------------------------------------------------------
        n_mat=>ndiso_mat(imode)
        n_fac=>ndiso_fac(imode)
c-----------------------------------------------------------------------
c       compute the real and imaginary parts separately.
c-----------------------------------------------------------------------
        real_imag: DO i_ri=0,1
          IF (keff(imode)==0.AND.i_ri==1) CYCLE mode_loop
          IF (i_ri==0) THEN
            flag='real'
          ELSE
            flag='imag'
          ENDIF
c-----------------------------------------------------------------------
c         evaluate the product of the matrix and the old density.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),nd(ibl),flag,imode)
          ENDDO
          CALL matvec(n_mat,vectr,sln,1_i4)
c-----------------------------------------------------------------------
c         add the product of the matrix and the old vector to the rhs,
c         then set the initial guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),crhs(ibl),flag,
     $                              imode,nqty=1_i4)
            IF (nd_diff>0.AND.ASSOCIATED(sln(ibl)%arri)) sln(ibl)%arri=0
            CALL vector_add(vectr(ibl),sln(ibl))
            CALL extrap_eval_vec(sln(ibl),ibl,imode,n_pass,flag)
          ENDDO
c-----------------------------------------------------------------------
c         call solver.
c-----------------------------------------------------------------------
          CALL iter_cg_2d_solve(n_mat,n_fac,sln,vectr,1_i4,
     $                          tol,maxit,solver,err,its,seed)
          IF (node == 0 .AND. itflag) CALL iter_out(n_pass,seed,its,err)
          ndits=its+ndits
c-----------------------------------------------------------------------
c         save the solution.
c-----------------------------------------------------------------------
          IF (err>tol) THEN
            converged=.false.
            RETURN
          ELSE
c-----------------------------------------------------------------------
c           determine the interior data if eliminated for the matrix
c           solve.
c-----------------------------------------------------------------------
            IF (poly_degree>1.AND.n_mat%eliminated)
     $        CALL fe_postsolve(n_mat,vectr,sln,
     $                          nd,crhs,1_i4,imode,flag)
            update: DO ibl=1,nbl
              SELECT CASE(n_pass)
c-----------------------------------------------------------------------
c             predictor step saves the fn_vdgn-centered result in work3:
c-----------------------------------------------------------------------
              CASE('n pre')
                CALL vector_assign_cvec(vectr(ibl),nd(ibl),flag,imode)
                CALL vector_add(sln(ibl),vectr(ibl),v1fac=fn_vdgn,
     $                          v2fac=1-fn_vdgn)
                CALL cvector_assign_vec(work3(ibl),sln(ibl),flag,imode)
c-----------------------------------------------------------------------
c             corrector step:
c-----------------------------------------------------------------------
              CASE('n cor')
                CALL cvector_assign_vec(nd(ibl),sln(ibl),flag,imode)
c-----------------------------------------------------------------------
c             trap unexpected pass flags.
c-----------------------------------------------------------------------
              CASE DEFAULT
                CALL nim_stop
     $          (TRIM(n_pass)//' pass is inappropriate for adv_nd.')
              END SELECT
            ENDDO update
          ENDIF
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     update data for polynomial extrapolation.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        SELECT CASE(n_pass)
        CASE('n pre')
          crhs(ibl)=work3(ibl)
          CALL vector_add(crhs(ibl),nd(ibl),1/fn_vdgn,1-1/fn_vdgn)
          CALL extrap_update(crhs(ibl),ibl,n_pass)
          IF (ibl<=nrbl) THEN
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qwork3,rb(ibl))
          ELSE
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qwork3,tb(ibl))
          ENDIF
        CASE('n cor')
          CALL extrap_update(nd(ibl),ibl,n_pass)
        END SELECT
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_nd
c-----------------------------------------------------------------------
c     subprogram 17. adv_nd_nsym.
c     number density advance with a non-Hermitian 2D operator.
c-----------------------------------------------------------------------
      SUBROUTINE adv_nd_nsym(n_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_cg
      USE iter_gmres_c2d
      USE extrap_mod
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: n_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,iv,imode,nq
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed
c-----------------------------------------------------------------------
c     set the quantity size to two if using hyper-diffusivity.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        nq=2
      ELSE
        nq=1
      ENDIF
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(csln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,nq)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,nq)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,nq,nmodes)
        IF (nq>1) CALL vector_type_alloc(cvecn(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,nq,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,nq)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,nq)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,nq,
     $                         nmodes)
        IF (nq>1) CALL vector_type_alloc(cvecn(ibl),1_i4,
     $                         tb(ibl)%mvert,0_i4,nq,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve.
c-----------------------------------------------------------------------
      CALL extrap_sln(n_pass,t+dt)
c-----------------------------------------------------------------------
c     create the particle advection operator.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.n0_changed.AND.nd_floor>0).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
          CALL matrix_create(nd_cmat,nd_cfac,cont_op,dirichlet_comp_op,
     $                       'sd sd','particle flow',.true.,solver,
     $                       mass_type='none')
        ELSE
          CALL matrix_create(nd_cmat,nd_cfac,cont_op,no_comp_mat_bc,
     $                       ' ','particle flow',.true.,solver,
     $                       mass_type='none')
        ENDIF
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side.
c-----------------------------------------------------------------------
      IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        CALL get_rhs(ndrhs,crhs,dirichlet_rhs,'sd sd','scalar',.false.,
     $               no_surf_int,cmat_elim=nd_cmat)
      ELSE
        CALL get_rhs(ndrhs,crhs,no_rhs_bc,' ','scalar',.false.,
     $               no_surf_int,cmat_elim=nd_cmat)
      ENDIF
c-----------------------------------------------------------------------
c     Loop over modes.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
c-----------------------------------------------------------------------
c       evaluate the product of the matrix and the old density.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          IF (nq>1) cvectr(ibl)=0._r8
          CALL cvector_2D_assign_cvec(cvectr(ibl),nd(ibl),imode,1_i4)
        ENDDO
        CALL matvec(nd_cmat(imode),cvectr,csln,nq)
c-----------------------------------------------------------------------
c       add the product of the matrix and the old vector to the rhs,
c       then set the initial guess.
c-----------------------------------------------------------------------
        DO ibl=1,nbl
          CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
          IF (ASSOCIATED(csln(ibl)%arri)) csln(ibl)%arri=0
          CALL vector_add(cvectr(ibl),csln(ibl))
          IF (nq>1) csln(ibl)=0._r8
          CALL cvector_2D_assign_cvec(csln(ibl),nd(ibl),imode,1_i4)
        ENDDO
c-----------------------------------------------------------------------
c       call solver.
c-----------------------------------------------------------------------
        CALL iter_gmr_c2d_solve(nd_cmat(imode),nd_cfac(imode),csln,
     $                          cvectr,nq,tol,maxit,solver,err,its,seed)
        IF (node == 0 .AND. itflag) CALL iter_out(n_pass,seed,its,err)
        ndits=its+ndits
c-----------------------------------------------------------------------
c       save the solution.
c-----------------------------------------------------------------------
        IF (err>tol) THEN
          converged=.false.
          RETURN
        ELSE
c-----------------------------------------------------------------------
c         determine the interior data if eliminated for the matrix
c         solve.
c-----------------------------------------------------------------------
          IF (poly_degree>1) THEN
            IF (nq>1) THEN
              DO ibl=1,nbl
                cvecn(ibl)=0._r8
                CALL cvector_assignq_cvec(cvecn(ibl),nd(ibl),1_i4)
              ENDDO
              CALL fe_postsolve(nd_cmat(imode),cvectr,csln,
     $                          cvecn,crhs,nq,imode)
            ELSE
              CALL fe_postsolve(nd_cmat(imode),cvectr,csln,
     $                          nd,crhs,1_i4,imode)
            ENDIF
          ENDIF
          update: DO ibl=1,nbl
            CALL cvector_assign_cvec2(nd(ibl),csln(ibl),imode,1_i4)
          ENDDO update
        ENDIF
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
        IF (nq>1) CALL vector_type_dealloc(cvecn(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_nd_nsym
c-----------------------------------------------------------------------
c     subprogram 18. adv_nd_3dnsym.
c     advance number density with implicit advection and 3D variations
c     in the flow field. Fourier coupling appears on the lhs of the
c     equation, these cases require a 3d matrix solution.
c
c     the 3D matrix is never formed; its action on a vector is found by
c     a finite element rhs computation, with the integrand specified in
c     the iter_ky_c3d_solve input.
c-----------------------------------------------------------------------
      SUBROUTINE adv_nd_3dnsym(n_pass,converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE surface_ints
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE extrap_mod
      USE pardata
      USE iter_ky_c3d_mod
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: n_pass
      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,imode,nq
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed
c-----------------------------------------------------------------------
c     interface block for the external subroutine that calls get_rhs
c     for the matrix-free dot-product computation.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_n_3dnsym(oper,prod,bc_oper)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(INOUT) :: oper,prod
        LOGICAL, INTENT(IN) :: bc_oper
        END SUBROUTINE threed_n_3dnsym
      END INTERFACE
c-----------------------------------------------------------------------
c     interface block for the external subroutine that applies
c     preconditioning operations.
c-----------------------------------------------------------------------
      INTERFACE
        SUBROUTINE threed_pren_3dnsym(resd,zeed,iiter,flex)
        USE vector_type_mod
        USE local
        TYPE(cvector_type), DIMENSION(:), INTENT(IN) :: resd
        TYPE(cvector_type), DIMENSION(:), INTENT(OUT) :: zeed
        INTEGER(i4), INTENT(IN) :: iiter
        LOGICAL, INTENT(IN) :: flex
        END SUBROUTINE threed_pren_3dnsym
      END INTERFACE
c-----------------------------------------------------------------------
c     set the quantity size to two if using hyper-diffusivity.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        nq=2
      ELSE
        nq=1
      ENDIF
c-----------------------------------------------------------------------
c     make space for the rhs computation.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,nq,nmodes)
        IF (nq>1) CALL vector_type_alloc(cvecn(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,nq,nmodes)
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,nq,
     $                         nmodes)
        IF (nq>1) CALL vector_type_alloc(cvecn(ibl),1_i4,
     $                         tb(ibl)%mvert,0_i4,nq,nmodes)
      ENDDO
c-----------------------------------------------------------------------
c     create the 2D operator for each Fourier component used for the
c     preconditioner in the iterations.
c-----------------------------------------------------------------------
      IF ((dt/=dt_old.OR.n0_changed.AND.nd_floor>0.OR.v0_changed).AND.
     $    istep>istep_mat.OR.force_recomp) THEN
        IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
          CALL matrix_create(nd_cmat,nd_cfac,cont_op,dirichlet_comp_op,
     $                       'sd sd','particle flow',.true.,solver,
     $                       mass_type='none')
        ELSE
          CALL matrix_create(nd_cmat,nd_cfac,cont_op,no_comp_mat_bc,
     $                       ' ','particle flow',.true.,solver,
     $                       mass_type='none')
        ENDIF
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     evaluate the right-hand side for all modes.
c
c     cell-interior data is not eliminated for 3d solves.
c-----------------------------------------------------------------------
      IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0)) THEN
        CALL get_rhs(ndrhs,crhs,dirichlet_rhs,'sd sd','scalar',.false.,
     $               no_surf_int)
      ELSE
        CALL get_rhs(ndrhs,crhs,no_rhs_bc,' ','scalar',.false.,
     $               no_surf_int)
      ENDIF
c-----------------------------------------------------------------------
c     use a polynomial fit to guess a solution for the linear solve,
c     and transfer it to work3 as the initial guess.
c
c     the nonlinear correction steps are special cases that use the
c     difference between the first solution and the old solution as
c     the guess, but they need nd reset before the solve.
c-----------------------------------------------------------------------
      IF (n_pass/='n cor dart') THEN
        CALL extrap_sln(n_pass,t+dt)
        DO ibl=1,nbl
          CALL extrap_eval_cvec(work3(ibl),ibl,n_pass)
          CALL vector_add(work3(ibl),nd(ibl),v2fac=-1._r8)
          IF (nd_bc=='dirichlet'.AND.(nd_diff>0.OR.nd_hypd>0))
     $      CALL dirichlet_rhs(work3(ibl),seam(ibl),'sd',1_i4)
        ENDDO
      ELSE
        DO ibl=1,nbl
          work3(ibl)=nd(ibl)
          nd(ibl)=nd_old(ibl)
          CALL vector_add(work3(ibl),nd_old(ibl),v2fac=-1._r8)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     initialize the 2-vector for the guess and solution if hyper-
c     diffusivity is used.
c-----------------------------------------------------------------------
      IF (nq>1) THEN
        DO ibl=1,nbl
          work5(ibl)=0._r8
          cvecn(ibl)=0._r8
          CALL cvector_assignq_cvec(work5(ibl),work3(ibl),1_i4)
          CALL cvector_assignq_cvec(cvecn(ibl),nd(ibl),1_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     Solve the 3D equation.
c-----------------------------------------------------------------------
      IF (nq==1) THEN
        CALL iter_ky_c3d_solve(threed_n_3dnsym,threed_pren_3dnsym,
     $                         nd_cmat(1)%hermitian,work3,crhs,work2,
     $                         1_i4,poly_degree,0_i4,0_i4,nrbl,nbl,
     $                         nmodes,tol,maxit,nmodes_total,err,its,
     $                         seed,old_vec=nd)
      ELSE
        CALL iter_ky_c3d_solve(threed_n_3dnsym,threed_pren_3dnsym,
     $                         .false.,work5,crhs,work6,nq,poly_degree,
     $                         0_i4,0_i4,nrbl,nbl,nmodes,tol,maxit,
     $                         nmodes_total,err,its,seed,old_vec=cvecn)
      ENDIF
      IF (node==0.AND.itflag) CALL iter_out(n_pass,seed,its,err)
      IF (node==0) ndits=ndits+its
c-----------------------------------------------------------------------
c     test convergence of 3D equation.
c-----------------------------------------------------------------------
      IF (err>tol) THEN
        converged=.false.
        RETURN
      ENDIF
c-----------------------------------------------------------------------
c     transfer the solution from the 2-vector if needed.
c-----------------------------------------------------------------------
      IF (nq>1) THEN
        DO ibl=1,nbl
          CALL cvector_assignq_cvec(work3(ibl),work5(ibl),1_i4)
        ENDDO
      ENDIF
c-----------------------------------------------------------------------
c     save the solution.
c-----------------------------------------------------------------------
      SELECT CASE(n_pass(1:5))
c-----------------------------------------------------------------------
c     combine delta-n and n, and update data for polynomial
c     extrapolation.
c-----------------------------------------------------------------------
      CASE('n cor')
        DO ibl=1,nbl
          CALL vector_add(nd(ibl),work3(ibl))
          IF (n_pass/='n cor dart') 
     $      CALL extrap_update(nd(ibl),ibl,n_pass)
          CALL vector_type_dealloc(crhs(ibl))
          IF (nq>1) CALL vector_type_dealloc(cvecn(ibl))
        ENDDO
c-----------------------------------------------------------------------
c     trap unexpected pass flags.
c-----------------------------------------------------------------------
      CASE DEFAULT
        CALL nim_stop
     $    (TRIM(n_pass)//' pass is inappropriate for adv_nd_3dnsym.')
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_nd_3dnsym
c-----------------------------------------------------------------------
c     subprogram 19. project_ndiff.
c     this routine manages a projection of the effective particle
c     source density from artificial particle diffusivity.  the
c     resulting expansion is used to correct momentum and energy
c     errors from this artificial term.
c-----------------------------------------------------------------------
      SUBROUTINE project_ndiff(converged,setave)
      USE local
      USE input
      USE integrands
      USE surface_ints
      USE boundary
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE iter_cg
      USE pardata
      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: converged
      LOGICAL, INTENT(IN) :: setave

      INTEGER(i4) :: ibl,its,imode,ibe,i_ri,nq,ipstart,ipass
      REAL(r8) :: err
      CHARACTER(8) :: seed,comp_flag
c-----------------------------------------------------------------------
c     make space for the linear algebra.  also interpolate nd to the
c     qwork3 quadrature-point space.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,1_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,1_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,1_i4,nmodes)
        CALL rblock_qp_update(rb(ibl)%nd,rb(ibl)%qwork3,rb(ibl))
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,1_i4,
     $                         nmodes)
        CALL tblock_qp_update(tb(ibl)%nd,tb(ibl)%qwork3,tb(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     if hyper particle diffusivity is used, the projection needs to
c     be performed twice.
c-----------------------------------------------------------------------
      IF (nd_hypd>0) THEN
        ipstart=1
      ELSE
        ipstart=2
      ENDIF
      pass_loop: DO ipass=ipstart,2
c-----------------------------------------------------------------------
c       contruct the right side of the projection equation.
c-----------------------------------------------------------------------
        IF (nd_bc=='dirichlet') THEN 
          CALL get_rhs(scal_lapl_rhs,crhs,dirichlet_rhs,'sd','scalar',
     $                 .false.,no_surf_int,rmat_elim=sproj_mat)
        ELSE
          CALL get_rhs(scal_lapl_rhs,crhs,no_rhs_bc,' ','scalar',
     $                 .false.,no_surf_int,rmat_elim=sproj_mat)
        ENDIF
        DO ibl=1,nbl
          work3(ibl)=0._r8
        ENDDO
c-----------------------------------------------------------------------
c       begin the loop over modes and separate real and imaginary parts
c       of the projection.
c-----------------------------------------------------------------------
        mode_loop: DO imode=1,nmodes
          real_imag: DO i_ri=0,1 
            IF (i_ri==0) THEN
              comp_flag='real'
            ELSE
              comp_flag='imag'
            ENDIF
            IF (keff(imode)==0) THEN
              IF (i_ri>0) CYCLE mode_loop
              comp_flag='real'
            ENDIF
c-----------------------------------------------------------------------
c           transfer the current component of the rhs and zero-out sln,
c           which is used as the initial guess.
c-----------------------------------------------------------------------
            DO ibl=1,nbl
              CALL vector_assign_cvec(vectr(ibl),crhs(ibl),
     $                                comp_flag,imode)
              sln(ibl)=0._r8
            ENDDO
c-----------------------------------------------------------------------
c           call solver.
c-----------------------------------------------------------------------
            CALL iter_cg_2d_solve(sproj_mat(imode),sproj_fac(imode),sln,
     $                            vectr,1_i4,tol,maxit,solver,err,
     $                            its,seed)
            IF (node == 0 .AND. itflag)
     $        CALL iter_out('sc prj',seed,its,err)
c-----------------------------------------------------------------------
c           save solution.
c-----------------------------------------------------------------------
            IF (err>tol) THEN
              converged=.false.
              RETURN
            ELSE
c-----------------------------------------------------------------------
c             determine the interior data if eliminated for the matrix
c             solve.
c-----------------------------------------------------------------------
              IF (poly_degree>1)
     $          CALL fe_postsolve(sproj_mat(imode),vectr,sln,
     $                            work3,crhs,1_i4,imode,comp_flag)
c-----------------------------------------------------------------------
c             put the result in work3 for interpolation after solves
c             for this pass are complete.
c-----------------------------------------------------------------------
              DO ibl=1,nbl
                CALL cvector_assign_vec(work3(ibl),sln(ibl),comp_flag,
     $                                  imode)
              ENDDO
            ENDIF
          ENDDO real_imag
        ENDDO mode_loop
c-----------------------------------------------------------------------
c       scale the result in order to have div*(nd_diff*grad(n)) if
c       there is only second-order dissipation or to leave
c       nd_diff*n-nd_hypd*grad**2(n) in qwork3 after the first
c       pass if there is hyper dissipation.
c
c       on the last pass, also create the average of the current and
c       previous result if this follows the first advance of n.
c-PRE
c       numerical analysis shows that centering needs to be the same
c       as the diffusion term in the continuity equation, so use
c       fthc, although hyper-diffn is still problematic.
c-----------------------------------------------------------------------
        IF (ipass==1) THEN  !  first pass with hyper dissipation
          DO ibl=1,nrbl
            CALL vector_add(work3(ibl),nd(ibl),v1fac=-nd_hypd,
     $                      v2fac=nd_diff)
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qwork3,rb(ibl))
          ENDDO
          DO ibl=nrbl+1,nbl
            CALL vector_add(work3(ibl),nd(ibl),v1fac=-nd_hypd,
     $                      v2fac=nd_diff)
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qwork3,tb(ibl))
          ENDDO
        ELSE  !  second or only pass
          DO ibl=1,nrbl
            IF (ipstart==2) CALL vector_mult(work3(ibl),nd_diff)
            IF (setave) rb(ibl)%qndiffa%qpf=rb(ibl)%qndiff%qpf
            CALL rblock_qp_update(rb(ibl)%work3,rb(ibl)%qndiff,rb(ibl))
            IF (setave) rb(ibl)%qndiffa%qpf=
     $        fthc*rb(ibl)%qndiff%qpf+
     $        (1._r8-fthc)*rb(ibl)%qndiffa%qpf
          ENDDO
          DO ibl=nrbl+1,nbl
            IF (ipstart==2) CALL vector_mult(work3(ibl),nd_diff)
            IF (setave) tb(ibl)%qndiffa%qpf=tb(ibl)%qndiff%qpf
            CALL tblock_qp_update(tb(ibl)%work3,tb(ibl)%qndiff,tb(ibl))
            IF (setave) tb(ibl)%qndiffa%qpf=
     $        fthc*tb(ibl)%qndiff%qpf+
     $        (1._r8-fthc)*tb(ibl)%qndiffa%qpf
          ENDDO
          IF (nonlinear.AND.impladv.AND.setave) CALL ndiff_store('ave')
        ENDIF
      ENDDO pass_loop
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE project_ndiff
c-----------------------------------------------------------------------
c     subprogram 20. adv_v_hypv_matv.
c     perform a time-split hyper-viscosity advance, where the operator
c     is the composite of two vector Laplacians.  this version uses a
c     stored matrix for generating the rhs, instead of calling finite-
c     element machinery.
c-----------------------------------------------------------------------
      SUBROUTINE adv_v_hypv_matv(converged,force_recomp)
      USE local
      USE input
      USE physdat
      USE integrands
      USE boundary
      USE matrix_storage_mod
      USE finite_element_mod
      USE computation_pointers
      USE matrix_mod
      USE regularity
      USE global
      USE fields
      USE seam_storage_mod
      USE iter_gmres_c2d
      USE pardata
      USE time
      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: converged
      LOGICAL, INTENT(IN) :: force_recomp

      INTEGER(i4) :: ibl,its,imode,ibe,i_ri,nq,n_int
      INTEGER(i4), SAVE :: istep_mat=-1
      REAL(r8) :: err
      CHARACTER(8) :: seed,bc_flag,comp_flag
c-----------------------------------------------------------------------
c     make space for the linear algebra.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        CALL vector_type_alloc(sln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(vectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(csln(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(cvectr(ibl),poly_degree,
     $                         rb(ibl)%mx,rb(ibl)%my,3_i4)
        CALL vector_type_alloc(crhs(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
        CALL vector_type_alloc(cvecn(ibl),poly_degree,rb(ibl)%mx,
     $                         rb(ibl)%my,3_i4,nmodes)
        work1(ibl)=0._r8
        crhs(ibl)=0._r8
        cvecn(ibl)=0._r8
      ENDDO
      DO ibl=nrbl+1,nbl
        CALL vector_type_alloc(sln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(vectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(csln(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(cvectr(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4)
        CALL vector_type_alloc(crhs(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
        CALL vector_type_alloc(cvecn(ibl),1_i4,tb(ibl)%mvert,0_i4,3_i4,
     $                         nmodes)
        work1(ibl)=0._r8
        crhs(ibl)=0._r8
        cvecn(ibl)=0._r8
      ENDDO
      n_int=(poly_degree-1)**2
c-----------------------------------------------------------------------
c     set the boundary condition flag.
c-----------------------------------------------------------------------
      IF (flow_bc=='free-slip') THEN
        bc_flag='3vn'
      ELSE
        bc_flag='all'
      ENDIF
c-----------------------------------------------------------------------
c     create the separate matrices for each Fourier component with
c     storage order (r,z,+-phi), where +- indicates (real_pol,-imag_phi)
c     or (imag_pol,+real_phi).
c-----------------------------------------------------------------------
      IF (dt/=dt_old.AND.(istep>istep_mat.OR.force_recomp)) THEN
        CALL matrix_create(hypv_cmat,hypv_cfac,hyp_visc_op,
     $                     dirichlet_comp_op,bc_flag,'hypv',.true.,
     $                     solver,mass_type='none')
        istep_mat=istep
      ENDIF
c-----------------------------------------------------------------------
c     begin the loop over modes.  at most, the (real_r,real_z,-imag_phi)
c     and (imag_r,imag_z,real_phi) components couple and each set uses
c     the same operator.  for n=0 phi is uncoupled and all real
c     components are done simultaneously.
c-----------------------------------------------------------------------
      mode_loop: DO imode=1,nmodes
        real_imag: DO i_ri=0,1 
          IF (i_ri==0) THEN
            comp_flag='r12mi3'
          ELSE
            comp_flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loop
            comp_flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         evaluate the product of the matrix for explicit terms and the
c         old velocity field.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            CALL vector_assign_cvec(vectr(ibl),ve(ibl),comp_flag,imode)
          ENDDO
          CALL matvec(hypv_expmat(imode),vectr,sln,3_i4,do_net=.false.)
          IF (i_ri==0) THEN
            DO ibl=1,nbl
              CALL cvector_assign_vec(crhs(ibl),sln(ibl),
     $                                comp_flag,imode)
            ENDDO
          ELSE
            DO ibl=1,nbl
              CALL cvector_assign_vec(cvecn(ibl),sln(ibl),
     $                                comp_flag,imode)
            ENDDO
          ENDIF
c-----------------------------------------------------------------------
c         eliminate interiors for this comp_flag part of this Fourier
c         component.  The complex matrix includes coupling between
c         dV and grad**2(V), so this is not standard and requires the
c         extra cvecn storage.
c-----------------------------------------------------------------------
          IF (poly_degree>1) THEN
            DO ibl=1,nrbl
              csln(ibl)=0._r8
              CALL cvector_2D_assign_vec(csln(ibl),sln(ibl),comp_flag)
            ENDDO

            CALL timer(timestart)
            CALL matelim_presolve(hypv_cmat(imode),csln,cvectr,3_i4)
            CALL timer(timeend)
            time_stcon=time_stcon+timeend-timestart

            IF (i_ri==0) THEN
              DO ibl=1,nrbl
                CALL cvector_assign_cvec2(crhs(ibl),cvectr(ibl),imode)
              ENDDO
            ELSE
              DO ibl=1,nrbl
                CALL cvector_assign_cvec2(cvecn(ibl),cvectr(ibl),imode)
              ENDDO
            ENDIF
          ENDIF
        ENDDO real_imag
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     multiply by the coefficient for the rhs and communicate across
c     block boundaries.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_mult(crhs(ibl),-(0,1)*SQRT(dt*hyp_visc/fhyp_visc))
      ENDDO
      CALL edge_network(3_i4,nmodes,poly_degree-1_i4,.true.)
      DO ibl=1,nbl
        CALL vector_mult(cvecn(ibl),-(0,1)*SQRT(dt*hyp_visc/fhyp_visc))
        CALL edge_load_carr(cvecn(ibl),3_i4,1_i4,nmodes,
     $                      poly_degree-1_i4,seam(ibl))
      ENDDO
      CALL edge_network(3_i4,nmodes,poly_degree-1_i4,.false.)
      DO ibl=1,nbl
        CALL edge_unload_carr(cvecn(ibl),3_i4,1_i4,nmodes,
     $                        poly_degree-1_i4,seam(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     apply boundary and regularity conditions to the rhs vector.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(r0block_list)
        ibe=r0block_list(ibl)
        CALL regular_vec(crhs(ibe),seam(ibe),'cyl_vec',3_i4,nmodes,
     $                   nindex)
        CALL regular_vec(cvecn(ibe),seam(ibe),'cyl_vec',3_i4,nmodes,
     $                   nindex)
      ENDDO
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(crhs(ibe),seam(ibe),bc_flag,3_i4)
        CALL dirichlet_rhs(cvecn(ibe),seam(ibe),bc_flag,3_i4)
      ENDDO
c-----------------------------------------------------------------------
c     start the loop for the solver.
c-----------------------------------------------------------------------
      mode_loops: DO imode=1,nmodes
        real_imags: DO i_ri=0,1 
          IF (i_ri==0) THEN
            comp_flag='r12mi3'
          ELSE
            comp_flag='i12r3'
          ENDIF
          IF (keff(imode)==0) THEN
            IF (i_ri>0) CYCLE mode_loops
            comp_flag='real'
          ENDIF
c-----------------------------------------------------------------------
c         transfer the relevant part of the rhs, then set the initial
c         guess.
c-----------------------------------------------------------------------
          DO ibl=1,nbl
            IF (i_ri==0) THEN
              CALL cvector_2D_assign_cvec(cvectr(ibl),crhs(ibl),imode)
            ELSE
              CALL cvector_2D_assign_cvec(cvectr(ibl),cvecn(ibl),imode)
            ENDIF
            csln(ibl)=0._r8
          ENDDO
c-----------------------------------------------------------------------
c         call solver.
c-----------------------------------------------------------------------
          CALL iter_gmr_c2d_solve(hypv_cmat(imode),hypv_cfac(imode),
     $                            csln,cvectr,3_i4,tol,maxit,solver,
     $                            err,its,seed)
          IF (node==0.AND.itflag) CALL iter_out('hypv',seed,its,err)
          viscits=its+viscits
c-----------------------------------------------------------------------
c         save solution.
c-----------------------------------------------------------------------
          IF (err>tol) THEN
            converged=.false.
            RETURN
          ELSE
c-----------------------------------------------------------------------
c           determine the interior data if eliminated for the matrix
c           solve.
c-----------------------------------------------------------------------
            IF (poly_degree>1) THEN
              IF (i_ri==0) THEN
                CALL fe_postsolve(hypv_cmat(imode),cvectr,csln,
     $                            work1,crhs,3_i4,imode)
              ELSE
                CALL fe_postsolve(hypv_cmat(imode),cvectr,csln,
     $                            work1,cvecn,3_i4,imode)
              ENDIF
            ENDIF
c-----------------------------------------------------------------------
c           update ve.
c-----------------------------------------------------------------------
            DO ibl=1,nbl
              CALL vector_assign_cvec2(sln(ibl),csln(ibl),comp_flag)
              CALL vector_assign_cvec(vectr(ibl),ve(ibl),comp_flag,
     $                                imode)
              CALL vector_add(sln(ibl),vectr(ibl))
              CALL cvector_assign_vec(ve(ibl),sln(ibl),comp_flag,
     $                                imode)
            ENDDO
          ENDIF
        ENDDO real_imags
      ENDDO mode_loops
c-----------------------------------------------------------------------
c     deallocate linear algebra space.
c-----------------------------------------------------------------------
      DO ibl=1,nbl
        CALL vector_type_dealloc(sln(ibl))
        CALL vector_type_dealloc(vectr(ibl))
        CALL vector_type_dealloc(csln(ibl))
        CALL vector_type_dealloc(cvectr(ibl))
        CALL vector_type_dealloc(crhs(ibl))
        CALL vector_type_dealloc(cvecn(ibl))
      ENDDO
c-----------------------------------------------------------------------
c     enforce the external boundary conditions on the solution to
c     avoid accumulating finite tolerance errors.
c-----------------------------------------------------------------------
      DO ibl=1,SIZE(exblock_list)
        ibe=exblock_list(ibl)
        CALL dirichlet_rhs(ve(ibe),seam(ibe),bc_flag,3_i4)
      ENDDO
c-----------------------------------------------------------------------
c     complete regularity conditions on the relation between n=1 r and
c     phi components.
c-----------------------------------------------------------------------
      CALL regular_ave(ve,3_i4,nmodes,nindex)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE adv_v_hypv_matv
