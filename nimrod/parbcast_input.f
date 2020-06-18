c-----------------------------------------------------------------------
c     file parbcast_input.f
c
c     the broadcast_input routine was extracted from the parallel.f
c     file to improve modularity.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 1. broadcast_input
c     broadcast info read by proc 0 out of nimrod.in (namelist reads)
c       to all processors
c     every quantity in input module is broadcast in case it was read
c     IMPORTANT: anytime a variable gets added to input module, it MUST be
c       added to this routine as well
c-----------------------------------------------------------------------

      SUBROUTINE broadcast_input
      USE input
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ierror

c physics specifications

      CALL mpi_bcast(nonlinear,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(lin_nmax,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(zperiod,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ndens,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(eta_model,16)
      CALL mpi_bcast(elecd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(eta_ref_t,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(elecd_max,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(elecd_min,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(elecd_chodura,1,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(f_chodura,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(separate_pe,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(kin_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(iso_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(par_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gyr_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dvac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dexp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(ds_use,8)
      CALL mpi_bcast(xvac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vsink_rate,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(beta,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(pe_frac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(advect,8)
      CALL bcast_str(continuity,16)
      CALL mpi_bcast(nd_diff,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_hypd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(nd_bc,16)
      CALL bcast_str(ohms,8)
      CALL mpi_bcast(be0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(thetab,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(phib,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(lamprof,6)
      CALL mpi_bcast(lam0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rbreak,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(alpha,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(eq_flow,16)
      CALL mpi_bcast(pit_0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(pit_2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(pit_4,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(pres_2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(pres_4,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ve0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(thetav,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(phiv,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(eqflow_width,1,mpi_nim_real,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(glength,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gravity,3,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(loop_volt,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tloopv0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tloopv1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(i_desired,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(loop_rate,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(loop_rate2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(e_vertical,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(t_e_vert0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(t_e_vert1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(init_type,16)
      CALL mpi_bcast(bamp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nx,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ny,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nz,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(zero_bnorm,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(norm_flow,16)
      CALL bcast_str(flow_bc,16)
      CALL mpi_bcast(ferr_amp,nfield_err,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(ferr_phase,nfield_err,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(ferr_n,nfield_err,mpi_nim_int,0,
     $               mpi_comm_world,ierror)

c constants specifications

      CALL mpi_bcast(set_phys_constants,1,mpi_nim_logical,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(chrg_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(zeff_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mi_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(me_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gam_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(kblz_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mu0_input,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(c_input,1,mpi_nim_real,0,mpi_comm_world,ierror)

c other equilibrium specifications.

      CALL mpi_bcast(tor_eqja_fe,1,mpi_nim_logical,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(nedge,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ncoil,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(coil_current,ncoil_max,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(coil_r,ncoil_max,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(coil_z,ncoil_max,mpi_nim_real,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(tanh_byl,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tanh_prl,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tanh_pfrac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tanh_ndl,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tanh_nfrac,1,mpi_nim_real,0,mpi_comm_world,ierror)

c closure specifications

      CALL bcast_str(neoe_flag,8)
      CALL bcast_str(neoi_flag,8)
      CALL bcast_str(neoe_det,8)
      CALL bcast_str(neoi_det,8)
      CALL mpi_bcast(neo_debug,1,mpi_nim_logical,0,
     &     mpi_comm_world,ierror)
      CALL mpi_bcast(ng_ft,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mu_e,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mu_i,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(m_neo,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(n_neo,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(neo_rad,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(p_model,16)
      CALL mpi_bcast(k_perp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_pll,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_perpe,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_plle,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_perpi,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_plli,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_pll_max,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_pll_min,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_pll_ref_t,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(kprp_mnrat,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(k_cross,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tequil_rate,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ohm_heat,1,mpi_nim_logical,0,
     &     mpi_comm_world,ierror)
      CALL mpi_bcast(visc_heat,1,mpi_nim_logical,0,
     &     mpi_comm_world,ierror)
      CALL mpi_bcast(insulate,1,mpi_nim_logical,0,
     &     mpi_comm_world,ierror)
      CALL bcast_str(parvisc_model,16)
      CALL bcast_str(closure_model,12)
      CALL mpi_bcast(tdep_coul_log,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(coulomb_logarithm,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(magfac_ele,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(magfac_ion,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(delta0_ele,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(delta1_ele,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(gamma0_ele,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(gamma1_ele,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(delta0_ion,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(delta1_ion,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(gamma0_ion,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(gamma1_ion,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(tdep_tequil,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)

c numerical specifications

      CALL mpi_bcast(dtm,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dt_initial,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dt_stop,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dt_incr,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(v_cfl,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(mhdadv_alg,8)
      CALL mpi_bcast(nl_cfl_lim,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tmax,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(cpu_tmax,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nstep,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(npc,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fom,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(feta,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(split_resist,1,mpi_nim_logical,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(fvsc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(split_visc,1,mpi_nim_logical,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(fthc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fb_vxb,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fv_vdgv,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fp_vdgp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fn_vdgn,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(divbd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fdivb,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ndivb,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(kdivb_2_limit,1,mpi_nim_real,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(split_divb,1,mpi_nim_logical,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(dt_change_frac,1,mpi_nim_real,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(ave_change_limit,1,mpi_nim_real,0,mpi_comm_world,
     $               ierror)
      CALL mpi_bcast(n_dt_release,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(n_mat_update,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(si_in_v,1,mpi_nim_logical,0,mpi_comm_world,ierror)
      CALL bcast_str(siop_type,8)
      CALL mpi_bcast(si_fac_mhd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mhd_si_iso,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(si_fac_hall,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(si_fac_pres,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(si_fac_j0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(si_fac_nl,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(conform,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(lump_b,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(lump_all,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(ngr,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(integration_formula,8)
      CALL mpi_bcast(poly_degree,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divb,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divb_min,1,mpi_nim_int,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divb_max,1,mpi_nim_int,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(disc_dbd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(hyp_eta,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fhyp_eta,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(split_hypeta,1,mpi_nim_logical,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(hyp_dbd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fhyp_dbd,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(hyp_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fhyp_visc,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fdivv,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(fpvrt,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ddivv,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dpvrt,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divv,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divv_min,1,mpi_nim_int,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divv_max,1,mpi_nim_int,0,
     $               mpi_comm_world,ierror)
      CALL mpi_bcast(poly_divv_auto,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(met_spl,5)
      CALL mpi_bcast(r0dr_weight,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(transfer_eq,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nd_floor,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_exp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_dart_fac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_dart_upw,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(t_dart_upw,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(upw_aniso,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(upw_limit,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_floor_upw,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nd_width_upw,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(t_floor_upw,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(t_width_upw,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(p_computation,16)
      CALL mpi_bcast(fmx_drate,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nfdamp,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nd_nodal_floor,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(ti_nodal_floor,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(te_nodal_floor,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nd_correrr,1,mpi_nim_logical,0,
     $               mpi_comm_world,ierror)

c linear solver specifications

      CALL mpi_bcast(tol,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(maxit,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(solver,8)
      CALL bcast_str(vmhd_solver,8)
      CALL bcast_str(bmhd_solver,8)
      CALL bcast_str(temp_solver,8)
      CALL mpi_bcast(off_diag_fac,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(extrap_order,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nsym_pre_band,1,mpi_nim_int,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nsym_pre_rpass,1,mpi_nim_int,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nsym_pre_rfac,1,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(nsym_pre_rtype,12)
      CALL mpi_bcast(maxit_nl,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(tol_nl,1,mpi_nim_real,0,mpi_comm_world,ierror)

c grid specifications

      CALL bcast_str(gridshape,8)
      CALL bcast_str(eqfile,64)
      CALL bcast_str(periodicity,8)
      CALL mpi_bcast(xmin,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(xmax,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ymin,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ymax,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(xo,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(yo,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mx,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(my,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(mxpie,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(pieflag,8)
      CALL bcast_str(rimflag,8)
      CALL mpi_bcast(firstx,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(firsty,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(skew,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nxbl,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nybl,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(geom,3)
      CALL mpi_bcast(lphi,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(lin_nmodes,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(per_length,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(decompflag,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nlayers,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(qpack,15,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(wpack,15,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(amp,15,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(packbigr,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(dealiase,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(quadxx,4,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(quadyy,4,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(quadarc,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)

c output specifications

      CALL mpi_bcast(detflag,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(itflag,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(nhist,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ihist,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(jhist,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(hist_binary,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(hist_flush,8)
      CALL mpi_bcast(xy_stride,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(xt_stride,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(yt_stride,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(x0fac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(y0fac,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(xdraw_dir,64)
      CALL mpi_bcast(ndxout,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(dx_dir,64)
      CALL mpi_bcast(ndump,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(dump_file,128)
      CALL bcast_str(dump_name,64)
      CALL bcast_str(dump_dir,64)
      CALL mpi_bcast(dump_over,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL bcast_str(reset_file,128)
      CALL mpi_bcast(reset_time,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(reset_step,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(reset_eq_mesh,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)

      RETURN
      END SUBROUTINE broadcast_input
