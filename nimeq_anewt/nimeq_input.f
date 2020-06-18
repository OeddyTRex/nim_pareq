c-----------------------------------------------------------------------
c     file nimeq_input.f
c     this contains a module for the input parameters that are specific
c     to the Grad-Shafranov solver nimeq.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declaration.
c-----------------------------------------------------------------------
      MODULE input_eq
      USE local
      IMPLICIT NONE

      INTEGER(i4), PARAMETER :: neqcoil_max=2000
      INTEGER(i4), PARAMETER :: nvfbcoil_max =1000
      INTEGER(i4), PARAMETER :: nrfbcoil_max =1000

c=======================================================================
c     nimeq parameters.
c=======================================================================
c     the following parameters control the profile of flux-functions
c	for the equilibrium.
c
c     the R*B_phi function (F) is specified according to the f_model
c	selection:
c
c	"linear" :: F(psi)=f0+dfdpsi*psi , where psi is the dimensional
c		poloidal flux function (-disk flux / 2*pi).  for
c		cylindrical configurations it is -1*mu0*J_parallel/B
c		and has units of inverse length.
c
c	"quad_closed" :: F is a quadratic function of the normalized
c		disk flux.  in the closed-flux region, it is
c
c		F = f_open+f1_closed*(1-Y)+4*f2_closed*Y(Y-1)
c
c		where Y=ring_flux/closed_flux,
c		which varies from 0 at the magnetic axis to 1 at the
c		separatrix.  in the open-flux region, F=f_open.
c		the derivative at the separatrix is continuous if
c		f2_closed=0.25*f1_closed.
c
c	"cubic_closed" :: F uses a smooth Hermite cubic basis
c
c		F = f_open+(f_axis-f_open)*(1-3*Y**2+2*Y**3)
c
c		where Y=ring_flux/closed_flux in the closed flux region,
c		and in the open flux region, F=f_open.
c
c	NOTE:	using "cubic_closed" for both f_model and p_model or
c		for one with the other uniform leads to nonlinear
c		convergence problems--possibly due to lack of a unique
c		solution.

      CHARACTER(16) :: f_model="linear"
      REAL(r8) :: dfdpsi=0._r8
      REAL(r8) :: f0=0._r8
      REAL(r8) :: f_open=0._r8
      REAL(r8) :: f_axis=0._r8
      REAL(r8) :: f1_closed=0._r8
      REAL(r8) :: f2_closed=0._r8
c-----------------------------------------------------------------------
c     equilibrium pressure (multiplied by mu0) as a function of poloidal
c	flux is specified according to pres_model.  the selections are:
c
c	"constant" :: mu0*Peq is set to pressure0
c
c	"quad_closed" :: mu0*Peq is a quadratic function of the
c		normalized disk flux.  in the closed-flux region, it is
c
c		mu0*Peq = p_open+p1_closed*(1-Y)+4*p2_closed*Y(Y-1)
c
c		where Y=ring_flux/closed_flux, in the closed flux
c		region.  in the open-flux region, mu0*Peq=p_open.
c		the derivative at the separatrix is continuous if
c		p2_closed=0.25*p1_closed.
c
c	"cubic_closed" :: mu0*Peq is set to 
c
c		mu0*P = p_open+(p_axis-p_open)*(1-3*Y**2+2*Y**3)
c
c		where Y=ring_flux/closed_flux in the closed flux region,
c		and in the open flux region, it's p_open.
c
c	NOTE:	using "cubic_closed" for both f_model and p_model or
c		for one with the other uniform leads to nonlinear
c		convergence problems--possibly due to lack of a unique
c		solution.

      CHARACTER(16) :: pres_model="constant"
      REAL(r8) :: pressure0=0._r8
      REAL(r8) :: p_open=0._r8
      REAL(r8) :: p_axis=0._r8
      REAL(r8) :: p1_closed=0._r8
      REAL(r8) :: p2_closed=0._r8
c-----------------------------------------------------------------------
c     With the quad_closed, cubic_closed, and linear options for
c	f_model, the GS solve can adjust f1_closed, f1_axis, and dfdpsi,
c	respectively to converge on a desired value of plasma current.
c	This option is used if ip_equil is nonzero.

      REAL(r8) :: ip_equil=0._r8
c-----------------------------------------------------------------------
c     F and/or mu0*P profiles are read from the 1d.bin file that
c	is produced by fluxgrid if mfile>0.  however, f_model and/or
c	p_model also need to be set to "read1d" to use the corresponding
c	profile in the GS computation.  if either f0file or p0file,
c	these values replace the data on the magnetic axis that
c	fluxgrid deletes.  otherwise, the splines are extrapolated to
c	the magnetic axis.

      INTEGER(i4) :: mfile=0	! # of equilibrium data nodes in 1d.bin.
      INTEGER(i4) :: muse=1	! nimeq can skip & use every muse point.
      REAL(r8) :: f0file=0._r8	! R*B_phi on the magnetic axis.
      REAL(r8) :: p0file=0._r8	! mu0*P on the magnetic axis.
c-----------------------------------------------------------------------
c     with diverted equilibria, discontinuous P and F derivatives at
c	the separatrix can hamper nonlinear convergence.  the following
c	parameters can be used to smooth these derivatives.  in the
c	range of psin_smth*p < Y < 1, the P or F profiles are evaluated
c	with cubics that smoothly connect their values at psin_smth*p
c	to the open-field values.

      REAL(r8) :: psin_smthpp=1._r8  !  limit for P-prime smoothing
      REAL(r8) :: psin_smthfp=1._r8  !  limit for F-prime smoothing
c-----------------------------------------------------------------------
c     the equilibrium number density may be reset so that
c	P_eq/p0=(n_eq/ndens)**gamma_nimeq .  the default value of 0 for
c	gamma_nimeq leaves the equilibrium number density unchanged
c	from what was in the dump file.  the reference value of p0
c	is the value of pressure at the magnetic axis.

      REAL(r8) :: gamma_nimeq=0._r8
c-----------------------------------------------------------------------
c     if the configuration is topologically rectangular and has a
c	a periodic boundary, the following parameter is used to specify
c	the poloidal flux that threads the periodic gap.

      REAL(r8) :: per_flux=0._r8
c-----------------------------------------------------------------------
c     nimeq can compute the top or bottom half of a vertically
c	symmetric equilibrium, and those halves can be assembled with
c	stitch.  set symm_region to "top" to compute that top half or
c	"bottom" for the bottom half.
c
c	note that the btop_check option does not function with
c	symm_region set to "top" or "bottom."

      CHARACTER(7) :: symm_region="neither"
c-----------------------------------------------------------------------
c     to run nimeq for a free-boundary solve in toroidal geometry,
c	set gs_type to "free", and specify the necessary external
c	field that will provide horizontal and vertical stability
c	without an eddy current distribution along the border of the
c	domain.

      CHARACTER(5) :: gs_type="fixed"  !  or "free"

c     when running free-boundary computations, it is possible to
c	symmetrize the surface-flux updates about z=0 for up-down
c	symmetric configurations.  this is automatic with
c	symm_region set to "top" or "bottom," and it can be imposed
c	without the symm_region option by setting free_symmj to true.

      LOGICAL :: free_symmj=.false.
c-----------------------------------------------------------------------
c     magnetic field from coils outside the domain can be imposed
c	as inhomogeneous conditions in the Grad-Shafranov solve.
c	a distinction from imposing them in nimset is that if nimeq_tx
c	is used and coil_tx is false, the coil field is not transferred
c	to n=0.  the coil fields are then fixed even if resistive walls
c	are used in nimrod. 
c
c	the nimeq coils are only available in toroidal geometry.

      INTEGER(i4) :: neqcoil=0
      REAL(r8), DIMENSION(neqcoil_max) :: eqcoil_i=0._r8   !  (Amps)
      REAL(r8), DIMENSION(neqcoil_max) :: eqcoil_r=0._r8   !  (meters)
      REAL(r8), DIMENSION(neqcoil_max) :: eqcoil_z=0._r8
      LOGICAL :: coil_tx=.false.
c-----------------------------------------------------------------------
c     Free-boundary computations can use sets of coils for feedback
c     during nonlinear iteration. Two independent feedback systems can
c     be used, one for vertical control and one for radial
c     control.
c
c     The feedback system is a generalization of the feedback system
c     described in Jardin's textbook Computational Methods in Plasma
c     physics. The controller calculates the difference in the flux
c     at two observation points specified by r(v)fb_ri and r(v)fb_zi
c
c     The change in the current in the feedback coils is
c     calculated using the formula
c     delta I^n = r(v)fb_amp delta psi^n
c                +r(v)fb_damp (delta psi^n - delta psi(n-1))
c                +r(v)fb_iamp sum_i delta psi^n
c
c    The second term is a derivative like term, and provides dampening
c    which is critical for the stability of the feedback system.
c
c    The third term is an integral like term, but has not yet been
c    thoroughly tested. This term helps provide stability against sudden
c    changes. This could be useful, for example, when the x-point
c    location changes drastically.
c
c    r(v)fbcoil_r and r(v)fbcoil_z specify the locations of the r(v)
c    feedback coils.  r(v)fbcoil_i specifies the relative current in the
c    feedback coils.

c    Preliminary testing shows that rfb_amp=0.1 and rfb_damp = 1.0 works
c    well for gs_center = 0.5.

      INTEGER(i4) :: nrfbcoil=0
      REAL(r8) :: rfb_r1=0._r8,rfb_z1=0._r8
      REAL(r8) :: rfb_r2=0._r8,rfb_z2=0._r8
      REAL(r8), DIMENSION(nrfbcoil_max) :: rfbcoil_i=1._r8   !  (arb)
      REAL(r8), DIMENSION(nrfbcoil_max) :: rfbcoil_r=0._r8   !  (meters)
      REAL(r8), DIMENSION(nrfbcoil_max) :: rfbcoil_z=0._r8
      REAL(r8) :: rfb_amp=1.0
      REAL(r8) :: rfb_damp=0.0
      REAL(r8) :: rfb_iamp=0.0

      INTEGER(i4) :: nvfbcoil=0
      REAL(r8) :: vfb_r1=0._r8,vfb_z1=0._r8
      REAL(r8) :: vfb_r2=0._r8,vfb_z2=0._r8
      REAL(r8), DIMENSION(nvfbcoil_max) :: vfbcoil_i=1._r8 !  (arb)
      REAL(r8), DIMENSION(nvfbcoil_max) :: vfbcoil_r=0._r8 !  (meters)
      REAL(r8), DIMENSION(nvfbcoil_max) :: vfbcoil_z=0._r8
      REAL(r8) :: vfb_amp=1.0
      REAL(r8) :: vfb_damp=0.0
      REAL(r8) :: vfb_iamp=0.0
c-----------------------------------------------------------------------
c     coil-field information can also be read from a binary file, where
c	each record lists (R,Z,I) for one coil.  this is intended for
c	reading the set of current filaments that represent plasma
c	current and are written to nimeq_iint.bin after each free-
c	boundary solve.  [change the name of this file after it is
c	created to avoid over-writing.]

      CHARACTER(128) :: curfil_file="none"  !  default skips the option
c-----------------------------------------------------------------------
c     if the configuration is diverted, then the value of the poloidal
c	flux function does not determine whether a given position is
c	connected to the wall by open magnetic field.  the following
c	parameter is used to specify what algorithm, if any, is used
c	to check whether each integration point is in an open-flux
c	region.

c     when tracing field-lines, blen_check sets the maximum distance
c	(relative to a macroscopic scale) that an integration is run
c	before the field-line is considered closed.  the default works
c	for many cases, but larger values may be needed in some
c	applications.

c     there are now two possible algorithms.  "beq-trace" traces the
c	poloidal part of field lines to distinguish regions that are
c	connected to the wall.  "passive adv" carries out a PDE solve
c	for a passive scalar problem using B-poloidal as the
c	effective velocity field.  the latter can be much faster than
c	the former, but accuracy requires fine spatial resolution.

      CHARACTER(16) :: btop_check="none"  !  or "beq-trace"
					  !  or "passive adv"
      REAL(r8) :: blen_check=3._r8        !  dimensionless

c     the stiff_solver option traces magnetic field with LSODE's BDF
c	solver instead of its Adams solver.  it can be more robust, but
c	it is also significantly (2-3 times) slower.

      LOGICAL :: stiff_solver=.false.
c-----------------------------------------------------------------------
c     this set of parameters affects the computational algorithm
c	for the Grad-Shafranov solver.  eq_iters is the maximum number
c	of nonlinear iterations that the code will attempt in reaching
c	a normalized tolerace of gsh_tol.  note that gsh_tol should not
c	be less than the tolerance of the linear solver (tol), which
c	is specified in nimrod.in.  gscenter is the relaxation parameter
c	for the nonlinear iterations.

      INTEGER(i4) :: eq_iters=80
      REAL(r8) :: gsh_tol=1.e-7_r8
      REAL(r8) :: gscenter=0.9_r8

c     to avoid changing maxit between nimeq and nimrod runs, the
c	following is used as a limit on linear-solver iterations for
c	nimeq only.
c     nimeq now has its own solver (preconditioner) specification.
c	the sequential SuperLU library is recommended but its use
c	may be limited by hardware memory.

      INTEGER(i4) :: linmaxit=200
      CHARACTER(8) :: nimeq_solver='seq_slu' ! 'bl_diaga' is slower but
					     ! uses less memory.

c     linearized F*F' and mu0*R**2*P' terms are included in the
c	linear systems through numerical differencing when dellam
c	is positive.  the numerical differences represent the partial
c	of these nonlinear terms with respect to lambda=psi/R**2.
c	this approximate-Newton approach starts when the relative
c	error is less than linffp.  at larger error, the iteration is
c	the standard fixed-point scheme using the gscenter relaxation
c	parameter.
c
c	when these terms are used and ip_equil is set, error in the
c	plasma current can be weighted with respect to error
c	in solving the GS equation.  ip_weight is used for this
c	purpose.  with its default value of 0, the plasma current
c	is not part of the linearization, even when ip_equil/=0.

      REAL(r8) :: dellam=0._r8
      REAL(r8) :: linffp=5.e-3_r8
      REAL(r8) :: ip_weight=0._r8 ! units are length^2

c     the use_delst_init option applies the poloidal flux computation on
c	equilibrium current density, which is read from the dump file,
c	when initializing the nonlinear GS iterations.  it calls nimeq
c	option 2 at the start of option 1 computations.  this can help
c	reduce the number of nonlinear iterations in the GS solve.

      LOGICAL :: use_delst_init=.false.
c-----------------------------------------------------------------------
c     the following control output other than the nimplot/XDRAW-like
c	plotting.  write_iters is a flag for writing norms of nonlinear
c	residuals during iteration.

      LOGICAL :: write_iters=.true.

c     the nimeq_tx flag is like transfer_eq for nimset--in fact, it uses
c	the same subroutine.  if it's true, and the 'nonlinear' input
c	for nimrod is true, the equilibrium is transferred to the n=0
c	part of the solution.

      LOGICAL :: nimeq_tx=.false.

c     if ndeq_notx>0 this input specifies a uniform value of number
c	density to leave in the equilibrium arrays.  the rest is
c	transferred to the n=0 component of the solution.  if ndeq_notx
c	is 0, all of the number density is left in the equilibrium
c	arrays.

      REAL(r8) :: ndeq_notx=0._r8

c     the rbph_notx paramter specifies a uniform amount of R*Bphi
c	to leave in the equilibrium B array when nimeq_tx is true.

      REAL(r8) :: rbph_notx=0._r8

c     this character string is used as the root name when writing
c	dump files after solving the GS equation.

      CHARACTER(64) :: dumpgs_name="dump_gs"
c-----------------------------------------------------------------------
c     the following input are no longer used but are left in nimeq_input
c	for backward compatibility.

c     the centering for free-boundary flux is increased from 0
c	to gscenter over frcent_iters iterations of the nonlinear
c	solve when gs_type='free' to avoid large initial
c	swings.

      INTEGER(i4) :: frcent_iters=10  !  not used at this time

c     nbcflux is no longer used with surface flux solved simultaneously
c	during free-boundary computations.

      INTEGER(i4) :: nbcflux=2

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE input_eq

c-----------------------------------------------------------------------
c     subprogram 1. read_namelist_eq
c     open and read the namelist input.
c-----------------------------------------------------------------------

      SUBROUTINE read_namelist_eq(infile,echo_in)
      USE local
      USE input_eq
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: infile
      LOGICAL, INTENT(IN) :: echo_in

      INTEGER :: read_stat
c-----------------------------------------------------------------------
c     the following defines the namelist variables for nimeq.  any new
c     input parameters should be added to this list.
c-----------------------------------------------------------------------
      NAMELIST/nimeq_input/pressure0,dfdpsi,f0,gscenter,eq_iters,
     $  gsh_tol,write_iters,dumpgs_name,per_flux,f_model,f_open,f_axis,
     $  f1_closed,f2_closed,pres_model,p_open,p_axis,nimeq_tx,p1_closed,
     $  p2_closed,linmaxit,mfile,f0file,p0file,muse,gamma_nimeq,
     $  symm_region,nimeq_solver,btop_check,blen_check,psin_smthpp,
     $  psin_smthfp,use_delst_init,stiff_solver,eqcoil_i,eqcoil_r,
     $  eqcoil_z,neqcoil,coil_tx,gs_type,frcent_iters,nrfbcoil,nvfbcoil,
     $  rfb_r1,rfb_z1,rfb_r2,rfb_z2,vfb_r1,vfb_z1,vfb_r2,vfb_z2,
     $  ip_equil,rfbcoil_i,rfbcoil_r,rfbcoil_z,vfbcoil_i,vfbcoil_r,
     $  vfbcoil_z,rfb_amp,rfb_damp,rfb_iamp,vfb_amp,vfb_damp,vfb_iamp,
     $  curfil_file,free_symmj,ndeq_notx,dellam,ip_weight,linffp,
     $  rbph_notx
c-----------------------------------------------------------------------
c     open input file.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE=infile,STATUS='OLD',POSITION='REWIND')
c-----------------------------------------------------------------------
c     read namelist.
c-----------------------------------------------------------------------
      READ(UNIT=in_unit,NML=nimeq_input,IOSTAT=read_stat)
      IF (read_stat/=0) CALL nim_stop
     $   ('Error reading nimeq input file.')
c-----------------------------------------------------------------------
c     close input file and terminate.
c-----------------------------------------------------------------------
      CLOSE(in_unit)
      END SUBROUTINE read_namelist_eq

c-----------------------------------------------------------------------
c     subprogram 2. broadcast_eq_input
c     broadcast the variables in input_eq to all processors.
c-----------------------------------------------------------------------
      SUBROUTINE broadcast_eq_input
      USE input_eq
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      INTEGER(i4) :: ierror

      CALL mpi_bcast(pressure0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dfdpsi,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(f0,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gscenter,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(eq_iters,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gsh_tol,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(write_iters,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(dumpgs_name,64)
      CALL mpi_bcast(per_flux,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(f_model,16)
      CALL mpi_bcast(f_open,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(f_axis,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(f1_closed,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(f2_closed,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(pres_model,16)
      CALL mpi_bcast(p_open,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(p_axis,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(p1_closed,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(p2_closed,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(linmaxit,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nimeq_tx,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(mfile,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(f0file,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(p0file,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(muse,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(gamma_nimeq,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL bcast_str(symm_region,7)
      CALL bcast_str(nimeq_solver,8)
      CALL bcast_str(btop_check,16)
      CALL mpi_bcast(blen_check,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(psin_smthpp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(psin_smthfp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(use_delst_init,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(stiff_solver,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(eqcoil_i,neqcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(eqcoil_r,neqcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(eqcoil_z,neqcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(neqcoil,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(coil_tx,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(gs_type,5)
      CALL mpi_bcast(frcent_iters,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nvfbcoil,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfbcoil_i,nvfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(vfbcoil_r,nvfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(vfbcoil_z,nvfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_r1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_z1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_r2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_z2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(nrfbcoil,1,mpi_nim_int,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfbcoil_i,nrfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(rfbcoil_r,nrfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(rfbcoil_z,nrfbcoil_max,mpi_nim_real,0,
     $     mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_r1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_z1,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_r2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_z2,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_amp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_damp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rfb_iamp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_amp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_damp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(vfb_iamp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ip_equil,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(free_symmj,1,mpi_nim_logical,0,
     $     mpi_comm_world,ierror)
      CALL bcast_str(curfil_file,128)
      CALL mpi_bcast(ndeq_notx,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(dellam,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(linffp,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(ip_weight,1,mpi_nim_real,0,mpi_comm_world,ierror)
      CALL mpi_bcast(rbph_notx,1,mpi_nim_real,0,mpi_comm_world,ierror)

      RETURN
      END SUBROUTINE broadcast_eq_input
