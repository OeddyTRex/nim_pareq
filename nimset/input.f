c-----------------------------------------------------------------------
c     file input.f
c     declarations and default values of input parameters.
c     note that new parameters should be added to the lists in module
c     input and to the namelist declarations in subroutine 
c     read_namelist.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  input.
c     2.  read_namelist.
c-----------------------------------------------------------------------
c     subprogram 1. input
c     module containing the type declaration and default values for
c     all input variables.
c-----------------------------------------------------------------------
      MODULE input
      USE local
      IMPLICIT NONE

      INTEGER(i4), PARAMETER :: ncoil_max=200,nfield_err=200

c=======================================================================
c     grid parameters.
c=======================================================================
c     The type of grid is set by gridshape and supporting input:
c       'rect':  rectangular cross section with the periodicity of the
c		 x and y directions also specified as 'none', 'y-dir',
c		 or 'both'.
c		 The dimensions are set by xmin, xmax, ymin and ymax.
c		 This grid may be made nonorthogonal by setting skew/=0.
c		 To get a nonuniform grid, set firstx and firsty to the
c		 first delta-x and delta-y desired, respectively, 
c		 otherwise set these parameters to 0. 
c	'circ':  circular or annular cross section with minimum and
c		 maximum polar radii set by xmin and xmax.
c		 firstx and firsty have some functionality here, too.
c		 for toroidal geometries (see geom below), firstx sets
c		 an offset for the center of the grid from the geometric
c		 center of the cross section to match a Shafranov shift.
c		 firsty stretches the grid in the radial direction from
c		 its center.  Unlike the application to rect gridshapes,
c		 firsty is nondimensional, relative to uniform radial
c		 gridding.
c	'flux':  shape based on equilibrium magnetic field which is
c		 read from the eqfile.
c	'rect_cir': a physically smooth, i.e. circular or elliptical, 
c		 but logically rectangular grid with maximum dimensions
c		 determined by xmin, xmax, ymin, and ymax.  The npc
c		 numerical parameter controls the number of smoothing
c		 passes.
c	'rectquad': produces a logically rectangular, quadrilateral
c		 region with corners set by quadxx(1:4) and quadyy(1:4),
c		 progressing in the counter-clockwise through the four
c		 corners.  The segment between corners 3 and 4 may be
c		 shaped as an arc by setting quadarc to true.

      CHARACTER(8) :: gridshape="circ"
      CHARACTER(8) :: periodicity="none"
      REAL(r8) :: xmin=0
      REAL(r8) :: xmax=1
      REAL(r8) :: ymin=0
      REAL(r8) :: ymax=1
      REAL(r8) :: skew=0
      REAL(r8) :: firstx=0
      REAL(r8) :: firsty=0
      CHARACTER(64) :: eqfile="fluxgrid.dat"
      REAL(r8), DIMENSION(4) :: quadxx=0._r8,quadyy=0._r8
      LOGICAL :: quadarc=.false.
c-----------------------------------------------------------------------
c     Select a periodic linear system ('lin') with length specified by 
c	per_length or a toroidal system ('tor').  For circular grids the
c	major radius is set by xo, and yo controls the vertical offset. 
c	For flux grids, this is determined by the data, and for 
c	rectangular grids, it's set by xmin, xmax, ymin and ymax.  

      CHARACTER(3) :: geom="tor"
      REAL(r8) :: per_length=1
      REAL(r8) :: xo=2
      REAL(r8) :: yo=0
c-----------------------------------------------------------------------
c     The following parameters control the mesh refinement.  The two
c	parameters mx and my set the number of cells in the central
c	part of the grid.  For rectangular grids, mx and my are the
c	number of cells across the x and y directions.  For circular or
c	flux grids, they are the number of radial and poloidal cells,
c	respectively, in the polar region gridded by rectangles
c	(logically rectangular).  This region is annular in shape and
c	extends from a small hole at the magnetic axis to the separatrix
c	(or some fraction of the separatrix, depending on parameter
c	psihigh set in fluxgrid.in).  Note that mesh refinement set by
c	mx and my is independent of the	block decomposition set by nxbl
c	and nybl (later), but each block in the central rblock region
c	must have at least 3 cells in each direction.

      INTEGER(i4) :: mx=32
      INTEGER(i4) :: my=32
c-----------------------------------------------------------------------
c     For linear geometry cases with 1D equilibria and gridshape="circ",
c	one can pack around rational surfaces whose
c	value are specified by qpack(:).  The packing is Gaussian (i.e.,
c	a plot of dr vs. i shows a Gaussian) with amplitude and width 
c	of the Guassian by amp(:) and wpack(:). The amplitude is roughly
c	dr0/dr where dr0 is the value far away from the Gaussian and
c	wpack is the full-width given as a percentage of a = xmax-xmin.
c	For example, if mx=164 and you want to pack an extra 100 points 
c	around the q=1 rational surface within a width of 0.1*a then:
c		qpack(1)=1.
c		wpack(1)=0.1
c		amp(1)=15.625  
c       where amp = dr0/dr = 1/(164-100) / (0.1*1/100) = 15.625
c       A value of zero for any quantity turns off the packing.

      INTEGER(i4), PARAMETER :: npack=15
      REAL(r8), DIMENSION(npack) :: qpack=0._r8
      REAL(r8), DIMENSION(npack) :: wpack=0._r8
      REAL(r8), DIMENSION(npack) :: amp=0._r8

c     Because this type of packing is accomplished easily, the routines
c	are now used for other types of meshes, and the packing
c	spots specified by qpack refer to other independent variables:
c
c     For gridshape = "rect" the packing is in the
c	x-direction if firstx=0.  In these cases, qpack has nothing
c	to do with safety factor, it's just the x-values for packing
c	locations.
c     For gridshape = "circ" and toroidal geometry, packing specified
c	by qpack is with respect to the minor-radius coordinate
c	with origin at the center of the mesh.  Another subtlety
c	for these cases is that packing can be specified separately
c	for the inboard and outboard sides, allowing the user to
c	better follow flux-surface shapes.  The packing transitions
c	smoothly from one to the other over the azimuthal coordinate.
c	If any of the width parameters, wpack, are negative, then
c	for wpack(i)>0, the corresponding amp(i) and location qpack(i)
c	are applied to the outboard side.  For wpack(i)<0, they are
c	applied to the inboard side with width=|wpack(i)|.
c     For gridshape = "flux" and toroidal geometry, the routines can
c	be used to repack a fluxgrid mesh.  Here the independent
c	coordinate for packing is R on the outboard side of the
c	magnetic axis (the location of the cut in the fluxgrid mesh).
c	If the parameter packbigr is false, qpack specifies q-values
c	for packing.  If packbigr is true, qpack specifies R-values
c	on the outboard side, i.e. just the independent coordinate.

      LOGICAL :: packbigr=.false.
c
c     The switch to get this form of packing has been revised (5/06) to
c	be based on amp alone.
c-----------------------------------------------------------------------
c     The block decomposition of the region of rectangular cells is
c	controlled by nxbl and nybl [the number of rblocks in the x and
c	y (or radial and poloidal) directions, respectively].  The
c	distribution of blocks on parallel machines is controlled by
c	decompflag (0 = clumped, 1 = strided, 2 = specified).

c     With decompflag=2, parallel_block_init reads a file that specifies
c	the within-layer rank for each block.  That rank (process #
c	within the layer) for global block id is the single integer on
c	the id-th line of the file.  For plasma regions, the file
c	should be named block_plasma.sup.part, and if there are vacuum
c	regions (vblock_from344 branch), the file for vacuum region 
c	number # should be named block_vacuum#.sup.part.

      INTEGER(i4) :: nxbl=1
      INTEGER(i4) :: nybl=1
      INTEGER(i4) :: decompflag=0
c-----------------------------------------------------------------------
c     The gridding used in the hole at the magnetic axis for circular
c	or flux grids is controlled by this set of input parameters.
c	mxpie sets the radius of the central hole.  The radius is
c	mxpie*(radial dimension of the first cell in the rblock region).

c       Pieflag controls the grid implmenetation of the core/magnetic
c	axis:
c       1) If pieflag is 'rblock', then the core region is a pie of 
c		degenerate rectangles tacked onto the first set of
c		rblocks in the polar region of rectangles.  This
c		increases the number of rectangular cells to the
c		separatrix to mx+1.
c       2) If pieflag is 'tblock0' then the tblock of triangles is 
c		a simple set of tblock triangles (pie slices like
c		the prevous option) with the magnetic axis
c       	as a common point.
c       3) If pieflag is 'tblock1' then the tblock of triangles is read
c       	from a set of input files based on the triangle code
c       	to generate the pie central block.  To use this option
c		one needs to first run fluxgrid (gridshape must be
c		'flux'), then run nimset with pieflag='tblock0',
c		then run the triangle code, then rerun fluxgrid and
c		nimset with pieflag='tblock1'.  The GUI is very helpful
c		with this.

      INTEGER(i4) :: mxpie=1
      CHARACTER(8) :: pieflag="rblock"
c-----------------------------------------------------------------------
c     Gridding from the separatrix to the wall is controlled by this
c	set. If rimflag is set to "none", then the boundary is taken
c	to be the last rblock.  If rimflag is set to "tblock" then
c       the region between the last rblock and the wall is implemented
c	as a set of nbl_rim tblocks.  The tblock elements
c	are defined by a set of input files based on the triangle code.
c	If rimflag is set to "rblock", then the region between the
c	last rblock and the wall is implemented as a set of
c	nbl_rim rblocks.  However, the rblock option has not been
c	implemented.

c     For the tblock option, nbl_rim must be >=2 to get proper seam
c	sequencing, and gridshape must be 'flux'.  In addition, this
c	requires a sequence of steps similar to the pieflag='tblock1'
c	option.  First, run fluxgrid. Second, run nimset with
c	rimflag="none".  Third, run the triangle code.  Fourth, rerun
c	fluxgrid, then rerun nimset with rimflag="tblock".  The GUI is
c	very helpful with this.

c     rim_length is a minimum wall length used in constucting the
c       poly files associated with the outer boundary.  When negative,
c       an internal scheme finds a minimun length by inspecting the
c       rblock boundary and using twice this minimum length, 
c       otherwise rim_length is this minimum value.

      CHARACTER(8) :: rimflag="none"
      INTEGER(i4):: nbl_rim=2
      REAL(r8) :: rim_length=-1.0
c-----------------------------------------------------------------------
c     For nonlinear runs, set the number of cells in the phi-direction
c	through lphi.  The number of cells is 2**lphi, and the number
c       of toroidal modes after dealiasing is 2**lphi/3 + 1 (including
c	n=0) using integer algebra.
c     For linear runs, the number of Fourier modes is set directly with
c	lin_nmodes.  One may specify the maximum toroidal mode number
c	for these cases to get the range
c		lin_nmax-(lin_nmodes-1) <= n <= lin_nmax,
c	otherwise the default gives
c		0 <= n <= lin_nmodes-1

      INTEGER(i4) :: lphi=2
      INTEGER(i4) :: lin_nmodes=1
      INTEGER(i4) :: lin_nmax=0

c     zperiod is the fraction of the torus to simulate (as in BOUT++)
c	zperiod=1 simulates full torus
c	zperiod=n simulates 1/n of the torus
c     nimset just multiplies the wavenumber array keff by this.

      INTEGER(i4) :: zperiod=1

c     The following option allows nimrod to skip the dealiasing in
c	nonlinear computations (when dealiase is set to false).  In this
c	case, there are 2**lphi/2 Fourier components, which still skips
c	the n=2**lphi/2 Fourier component.

      LOGICAL :: dealiase=.true.
c-----------------------------------------------------------------------
c     For parallel computing, set the number of processor-layers into
c     which the Fourier components are distributed.

      INTEGER(i4) :: nlayers=1

c=======================================================================
c     physical constants.
c=======================================================================
c     This set of input parameters allow the user to change physical
c	constants at run-time.  Previously, these constants had be
c	fortran paramters that could only be changes by modifying
c	physdat and recompiling.  These input specifications allow
c	greater flexibility, including running normalized sets of
c	equations, but they also allow inconsistencies between pre-
c	processing and the calculation if the parameters are
c	inadvertently changed between the two steps.

c     The default values are MKS, except the Boltzmann constant is for
c	temperature in eV.  Note that the appearance of factors of c
c	are for MKS and not CGS.  If you really want CGS, the closest
c	possibility may be to consider the computed B as B/c
c	(Gauss-s/cm) and set mu0 to 4*pi/c**2.  Vacuum permittivity
c	is always computed from 1/(c**2*mu0).

c     Note that if the const_input namelist does not appear in the
c	nimrod.in file, the default values appearing in physdat.f
c	are used, not the input defaults below.  This is unique among
c	input parameters.

      REAL(r8) :: chrg_input=1.60217733e-19_r8	!  elementary charge
      REAL(r8) :: zeff_input=1._r8		!  relative ion charge
      REAL(r8) :: mi_input=3.3435860e-27_r8	!  ion mass
      REAL(r8) :: me_input=9.1093898e-31_r8	!  electron mass
      REAL(r8) :: gam_input=5._r8/3._r8		!  ratio of sp. heats
      REAL(r8) :: kblz_input=1.60217733e-19_r8	!  Boltzmann constant
      REAL(r8) :: mu0_input=pi*4.e-7_r8		!  vacuum permeability
      REAL(r8) :: c_input=2.99792458e8_r8	!  speed of light

c     The following is not an input parameter; do not put it in the
c	namelist.

      LOGICAL :: set_phys_constants=.false.

c=======================================================================
c     initialization parameters used in nimset.
c=======================================================================
c     The initial perturbation is controlled by the following set.  bamp
c	is the magnitude of the perturbed B.  nx and ny set the number
c	of half-wavelengths in the x and y (radial and poloidal) 
c	directions, respectively.  init_type selects what is being
c	initialized; though, you are not guaranteed an eigenfunction.
c	Note that "linear b" matches the DEBS perturbation and is only
c	suitable for linear geometry.

      REAL(r8) :: bamp=1.e-5
      INTEGER(i4) :: nx=1
      INTEGER(i4) :: ny=2

c     The parameter nz is used to select a specific Fourier component
c	to initialize. If nz<0, then all components are initialized.
c	Note that nz is the actual Fourier component mode number, so for 
c	zperiod=5, nz should be 0,5,10,15,... and not 0,1,2,3,...

      INTEGER(i4) :: nz=-1

      CHARACTER(16) :: init_type="compr alf" 
c	choose:	"compr alf", "shear alf",
c		"whistler", "linear b", "vertical"

c-----------------------------------------------------------------------
c     In addition to the perturbation specified through init_type, one
c	may add a field error to the magnetic-field solution data.  Note
c	that zero_bnorm must be set to false to avoid having nimrod
c	zero-out the field error along conducting surfaces.  Also, while
c	the distributions are vacuum fields, they are not force-free in
c	the presence of equilibrium current density.

c	Setting ferr_amp(1:nfield_err) to a value other than zero
c	activates this option in nimset.  ferr_n is the integer Fourier
c	index in the periodic coordinate, and ferr_phase (in units of
c	pi) shifts the phase of the perturbation.

c	In toroidal geometry, the field is the gradient of a magnetic
c	potential that only varies in the R-phi plane.  The radius at
c	which the R-component of the perturbation is ferr_amp is R=xo.
c	In linear geometry, the potential only varies in the x-z plane
c	(horizontal-periodic), and the x-component is ferr_amp at
c	xmin and xmax.

      REAL(r8), DIMENSION(nfield_err) :: ferr_amp=0._8
      REAL(r8), DIMENSION(nfield_err) :: ferr_phase=0._8
      INTEGER(i4), DIMENSION(nfield_err) :: ferr_n=1_i4
c-----------------------------------------------------------------------
c     The transfer_eq option allows one to transfer the equilibrium
c	field arrays (from an equilibrium file or an analytic model)
c	into the n=0 solution vector for use as an initial condition.
c	This is done with nimset, and may be done during a dump file
c	reset, but nonlinear must be true.

      LOGICAL :: transfer_eq=.false.

c=======================================================================
c     equilibrium parameters used in nimset.
c     Note that the following parameters also affect how the kernal
c     nimrod works, but they are first used in nimset:
c       beta
c       pe_frac
c       eq_flow
c=======================================================================
c     Shaping parameters for both resistivity and viscosity:  each
c	is the constant specified above times the function,
c	(1+(SQRT(dvac)-1)*(xs/mx)**dexp)**2, where xs is the radial
c	logical coordinate in the collection of polar rblocks, i.e., 
c	xs/mx ~ SQRT(psi).  In blocks outside the separatrix it is dvac.
c	(The character ds_use defined below controls how this profile
c	is used by the kernel.)

      REAL(r8) :: dvac=1
      REAL(r8) :: dexp=20

c     For circular gridshapes, the following allows one to put a 
c	vacuum region from xvac < r < xmax.  The dimension of the
c	current channel is then scaled to xvac instead of xmax.  The
c	default value gives no vacuum region.

      REAL(r8) :: xvac=0
c-----------------------------------------------------------------------
c     Parameters for the 0-th order plasma:  The electron number 
c	density, ndens, is in m**(-3).

      REAL(r8) :: ndens=1.e20
c-----------------------------------------------------------------------
c     Parameters for the 0-th order magnetic field for cases not reading
c	an equilibrium file:  be0 is the magnitude (at the polar origin
c	for circular grids).  For rectangular grids (with lam0=0), the
c	field is uniform, and its direction is given in polar angles
c	thetab and phib.

c	For circular grids in toroidal geometry, be0 is the magnitude
c	at R=xo.  It's direction in the Z-phi plane is controlled by
c	thetab; 0 is a purely toroidal vacuum field and 0.5 is
c	purely vertical.

c	For circular grids in linear geometry, these parameters specify 
c	the J/B profile.  lam0 is the J/B ratio at the magnetic axis
c	(multiplied by plasma radius).  Either rbreak or alpha are used 
c	to control the shape of the profile according to the choice of
c	lamprof.  For paramagnetic equilibria, just specify lam0 and
c	be0, and the parallel current will be determined as if E_axial
c	and resistivity were uniform.

c	These lambda profiles can be applied to rectangular
c	cases, too, if mx is even and the grid spacing is uniform in
c	x.  This produces a sheared slab model, which can be rotated
c	about the x axis by setting thetab/=0.

c       The pressure (and perpendicular current) are specified with the
c       beta, pres_2, and pres_4 parameters described for the pitprs
c       lamprof option below.

      REAL(r8) :: be0=1
      REAL(r8) :: thetab=0
      REAL(r8) :: phib=0.5

      CHARACTER(6) :: lamprof="mbfm  "  ! "mbfm", "alpha", "para"; for
c                     "pitprs", "qspec", "tanh", "sheet", "grbrap", and
c		      "jardel", see below.
      REAL(r8) :: lam0=0
      REAL(r8) :: rbreak=0.2	! break point (in a) for MBFM
      REAL(r8) :: alpha=3.	! lam=lam0*( (a-r)/a )**alpha
c-----------------------------------------------------------------------
c     Set the total plasma pressure based on the equilibrium magnetic
c	field (pres_eq and be_eq, respectively) at the magnetic axis for
c	equilibria not initialized from a fluxgrid file.  In ALL cases,
c	the perturbed pressure is advanced, and momentum is affected by
c	pressure gradients only when beta>0. 

      REAL(r8) :: beta=0

c     Set the electron pressure as a fraction of total.  For cases
c	where the electron pressure is advanced separately
c	(separate_pe=T), pe_frac sets the equilibrium and initial
c	electron pressures.  For cases without a separate advance,
c	electron pressure is pe_frac times the total pressure for all
c	time.  Note: set 0<=pe_frac<=1.
c
c       Note that with separate pressures (now separate temperatures)
c       and temperature equilibration, equilibrium electron and ion
c       temperatures should be the same.  Thus, pe_frac should be
c       zeff/(1+zeff).

      REAL(r8) :: pe_frac=0.5

c     For analytical 1D slab and circular profiles, pres_offset lets
c	one add a constant to the pressure profile.  This can be handy
c	for avoiding P_eq->0 at the wall.  The input value is normalized
c	such that pres_offset*p0 is added, where p0 is the central
c	pressure.

      REAL(r8) :: pres_offset=0.

c     For the same 1D slab and circular profiles, gamma_nimset can be
c	used to specify a number-density profile.  If gamma_nimset is
c	nonzero, the number density is set so that
c
c	  P_eq/p0=(n_eq/ndens)**gamma_nimset
c
c	where p0 is the central pressure.

      REAL(r8) :: gamma_nimset=0.
c-----------------------------------------------------------------------
c     For linear geometry with a circular cross section, one may create
c       nonzero beta equilibria with the pitch & pressure profile or
c       with the q and pressure profile (more convenient for tokamaks).
c       To invoke these profiles, set lamprof="pitprs" or "qspec".
c
c       For lamprof="pitprs":
c       Pitch is defined as (r*Bz/a*Bth), i.e. q*R/a, and it's
c       specified through the input variables pit_0, pit_2, and pit_4:
c
c       pitch(r)=pit_0*(1+pit_2*x**2+pit_4*x**4) ,
c
c       where x=(r-xmin)/(xmax-xmin)
c       The pressure is specified through beta, pres_2, pres_4 as
c
c       p(r)=beta*be0**2/(2*mu0)*(1+pres_2*x**2+pres_4*x**4) .
c
c       For lamprof="qspec":
c               See Holmes et.al., Phys. Fluid 26 (1983) 2569
c       The q profile is specified through the input
c       variables pit_0, pit_2, and pit_4:
c
c       q(r)=pit_0*(1+(x/pit_2)**(2*pit_4))**(1/pit_4) .
c
c       The pressure is specified through beta as
c
c       p(r)=beta*be0**2/(2*mu0)*Integrate[ x/q**3 * (2q - xq'),(x,1,r)]
c
c       Note that xmin must be zero for the qspec option.
c
      REAL(r8) :: pit_0=1
      REAL(r8) :: pit_2=0
      REAL(r8) :: pit_4=0
      REAL(r8) :: pres_2=0
      REAL(r8) :: pres_4=0
c-----------------------------------------------------------------------
c     Two related slab-only (geom=lin and gridshape=rect) equilibria
c	are available.  The magnetic field that is perpendicular to the
c	guide field is proportional to the hyperbolic tangent function
c	in both cases.  Other similarities and the differences are
c	described in the following:
c
c	1) lamprof='sheet' -- the domain is assumed to be symmetric
c	   about x=0, where the magnitude of the guide field and
c	   the background pressure are be0 and beta*be0**2/2*mu0,
c	   respectively.  The magnetic field profile is integrated
c	   numerically from the parallel current profile, and the
c	   prescription of By (shearing part) at walls is through
c		By_inf=be0*SIN(phib)
c	   based on the large guide-field limit.  Also set lam0 as
c		lam0=SIN(phib)*(x_wall)/(equilibrium-length scale).
c	   The polynomial pressure profile can be added with pres_2
c	   and pres_4.
c
c	2) lamprof='tanh' -- both the shearing component, By, and the
c	   pressure profile vary as hyperbolic tangents so that a
c	   diamagnetic drift can be imposed.  Equilibrium fields and
c	   currents are computed analytically, and the domain does
c	   not have to be symmetric about x=0.  Like the 'sheet'
c          profile, By limits to 
c            be0*SIN(phib)
c          for large |x|, which may be beyond the walls).  However, the
c          equilibrium length for By is set directly by parameter
c          tanh_byl.  The equilibrium length for P is tanh_prl.
c          Pressure at x=0 is again set by
c	     beta*be0**2/2*mu0  ,
c          and the parameter tanh_pfrac is the fractional change
c          in pressure over x=0 to x=inf.
c
c	3) lamprof='cos' -- the sheared component varies as a cosine 
c          but the pressure profile and density vary as tanh:
c           B_eq   = B_0(x)*[              phib*sin(x/tanh_byl) s^
c                             +sqrt(1-(phib*sin(x/tanh_byl))^2) g^]
c           g^ = cos(pi*thetab) phi^ + sin(pi*thetab) Z^
c           s^ =-sin(pi*thetab) phi^ + cos(pi*thetab) Z^
c           B_0(x) = B_0(0)*[1-beta*tanh_pfrac*TANH(x/tanh_prl) ]^(1/2)
c           p_0(x) = beta*B_0**2/(2*mu0)*[1+tanh_pfrac*TANH(x/tanh_prl)]
c           n_0(x) =               ndens*[1+tanh_nfrac*TANH(x/tanh_ndl)]
c
c          The guide field and sheared field amplitudes are controlled 
c          by the phib parameter and tanh_byl sets the shear-length.
c          If variation is allowed only in Z (lin_nmodes=1,lin_nmax=0)
c          then k_parallel can be set using thetab to rotate B. 
c          Alternatively, we can keep thetab=0 and control k_par with 
c          lin_nmax=1 and changing per_length.
c          Pressure at x=0 is again set by beta*be0**2/2*mu0 and density
c          is set by ndens. The length scale and fractional change are 
c          set by tanh_prl and tanh_pfrac for pressure, and 
c          tanh_ndl and tanh_nfrac.

      REAL(r8) :: tanh_byl=1._r8	!  eq. By length in m
      REAL(r8) :: tanh_prl=1._r8	!  eq. P length in m
      REAL(r8) :: tanh_pfrac=0.25_r8	!  |P(x=inf)-P(0)|/P(0)
      REAL(r8) :: tanh_ndl=1._r8	!  eq. n-gradient length in m
      REAL(r8) :: tanh_nfrac=0.25_r8	!  |n(x=inf)-n(0)|/n(0)

c     As with other slab equilibria, these two profiles may be rotated
c	about the x axis using the thetab input parameter.
c-----------------------------------------------------------------------
c     The "grbrap" cylindrical profile is from the Gruber and Rappaz
c	book, "Finite Element Methods in Linear Ideal MHD," Springer-
c	Verlag 1985.  The axial field is uniform, and the Btheta
c	and P are set by two parameters, c1 and c2:
c
c	Bth=be0*c1*r/(1+c2**2*r**2)
c	Bz =be0
c       P  =0.5*(be0**2/mu0)*(c1**2/c2**2)
c			    *[1/(1-c2**2*r**2)**2-1/(1-c2**2)**2]
c
c       The un-normalized pitch profile with this equilibrium is
c	q(r)*R=(1+c2**2*r**2)/c1.
c
c	Instead of defining new parameters, the pit_0 parameter is
c	used to set 1/((xmax-xmin)*c1), which is the normalized pitch on
c	axis, and pit_2 is used to set c2; though, it is not directly
c	related to pitch.
c-----------------------------------------------------------------------
c     The "jardel" profile is the cylindrical equilibrium from the
c	DeLucia, Jardin, and Glasser paper, Resistive stability of the
c	cylindrical spheromak, Phys. Fluids 27, 1470 (1984).  The
c	pitch profile is specified as pitch(r)=pit_0*(1-r**2), and the
c	pressure gradient is alpha times the Suydam-marginal value at
c	every point in r.  The equilibrium B and P are computed via
c	integration.  Note that this profile requires xmin=0 and xmax=1.
c-----------------------------------------------------------------------
c     Select an equilibrium flow profile.  This is treated like all
c	other _eq quantities--bicubic and fixed in time.  The choices
c	are:
c	  1) 'none'
c	  2) 'uniform'=>This produces a uniform flow with magnitude
c		ve0, specified in MKS units (m/s).
c		For gridshape='rect', the direction is set by thetav
c		and phiv.  For the other grid shapes, it is in the
c		direction normal to poloidal plane.  When the geometry
c		is toroidal, the flow is rigid rotation.
c	  3) 'gauss'=>This gives a Gaussian flow distribution centered
c		at the axis of polar grids (or the horizontal center 
c		of rectangular grids) with width eqflow_width, relative
c		to the root of the flux function (approximately).
c		The peak of the profile has magnitude ve0, and the
c		direction with different gridshapes is the same as for
c		the uniform option.
c	  4) 'pinch'=>An axial or toroidal electric field is computed
c		from the resistivity and equilibrium current density at
c		the magnetic axis.  ve_eq is then solved from the
c		perpendicular part of E_tor=ve_eqXbe_eq - eta*j_eq.
c		Caution:  E_tor is considered to be ~ 1/R, whether or
c		not this is appropriate for a given equilibrium. 
c	  5) 'eq_file'=>The equilibrium flow is taken from the
c		equilibrium solver file where possible.
c	  6) 'diamagnetic'=>Compute the ion diamagnetic flow from the
c		equilibrium B and J, assuming J_perp is from grad(P)
c		and using P_i=(1-pe_frac)*P.

c     Note that one must set the advect variable to get equilibrium
c	flow in the momentum equation, regardless of the eq_flow
c	selection.

      CHARACTER(16) :: eq_flow='none'
      REAL(r8) :: ve0=0
      REAL(r8) :: thetav=0
      REAL(r8) :: phiv=0.5
      REAL(r8) :: eqflow_width=0.5
c-----------------------------------------------------------------------
c     Set the equilibrium length for the gravitational mode equlibrium.
c	The first vector component of the gravity parameter, declared
c	with physics input, also affects the profile.

      REAL(r8) :: glength=0.	!  in meters

c=======================================================================
c     physics parameters.
c=======================================================================
c     If nonlinear is false, the equilibrium field is held fixed in 
c	time, and electric fields and Lorentz forces use only linear 
c	terms.

      LOGICAL :: nonlinear=.true.
c-----------------------------------------------------------------------
c     Select the terms to appear in the generalized Ohm's law.  Valid
c	choices are: 1) 'mhd', 2) 'hall', 3) 'mhd&hall', or 4) '2fl'
c	Options 2) - 4) now include the electron pressure term
c	term when beta>0, and 4) adds electron inertia.

      CHARACTER(8) :: ohms="mhd"
c-----------------------------------------------------------------------
c     The eta_model input allows the user to select whether to use
c	a temperature (electron) dependent resistivity, proportional
c	to T**-3/2.  Valid input options are:

c	Note that for linear computations, the "eta full" option will
c	include the linear resistivity perturbation according to the
c	linear electron temperature perturbation, while the
c	"eta n=0 only" option will not.  Both will use an unperturbed
c	resistivity profile based on the equilibrium temperature,
c	however.

c	1) "fixed"		skip the T-dependence
c	2) "eta n=0 only"	use only the symmetric part of T
c	3) "eta full"		use 3D eta(T)
c	4) "chodura"		uses the Chodura model--see below

c	For options 2 & 3, the electrical diffusivity has the form

c	D=MAX( MIN( elecd*(T/eta_ref_t)**(-1.5), elecd_max ), elecd_min)

c	The elecd coefficient is the electrical diffusivity in m**2/s
c	(MKS resistivity/mu0).

      CHARACTER(16) :: eta_model="fixed"
      REAL(r8) :: elecd=0
      REAL(r8) :: elecd_max=1000
      REAL(r8) :: elecd_min=0
      REAL(r8) :: eta_ref_t=1
      REAL(r8) :: f_chodura=1
      REAL(r8) :: elecd_chodura=1

c	When eta_model="chodura" the electrical diffusivity is computed
c	from the sum of the "eta full" T**(-1.5) model and the
c	phenomenological Chodura model [Nucl. Fusion 15, 65 (1975)]:

c	elecd_chodura*SQRT(ndens/n)*(1-EXP(-f_chodura*v_drift/v_ia))

c	Where v_drift is the electron drift speed, |J/ne|, and v_ia
c	is the ion-acoustic speed, SQRT(gamma*P/mtot*n), each computed
c	locally in 3D.  At present the drift speed computation is not
c	centered in the time-advance of B, so this contribution is
c	first-order accurate in dt.
c-----------------------------------------------------------------------
c
c     Select the independent calculation of the electron pressure
c
      LOGICAL :: separate_pe=.FALSE. 
c-----------------------------------------------------------------------
c     Select the advection terms to appear in the the momentum equation
c	and the generalized Ohm's law.  Valid choices are: 1) 'none', 
c	2) 'V only', 3) 'J only', or 4) 'all'.

      CHARACTER(8) :: advect="none"
c-----------------------------------------------------------------------
c     Select how the continuity equation for number density is applied.
c	Valid choices are:
c	1) 'none' density perturbations are not evolved.
c	2) 'fix profile' evolve density but use only the equilibrium
c	   density in the velocity advance,
c	3) 'n=0 only' evolve density but use only the equilibrium plus
c	   n=0 part of density in the velocity advance.
c	4) 'full' evolve density and use the full 3D density profile in
c	   the velocity advance.  This is the only self-consistent
c	   option for nonlinear runs, but it's also the most
c	   computationally intensive, since the velocity advance
c	   requires a 3D matrix solve at each step.

      CHARACTER(16) :: continuity="none"
c-----------------------------------------------------------------------
c     Isotropic diffusivity used in the number density evolution
c	(continuity) equation that crudely represents particle transport
c	for long time-scale simulations.  Units are m**2/s.
c     The nd_hypd is a hyper-diffusivity for number density evolution
c	and is only applied to computations with implicit advection,
c	i.e. ohms/=mhd or mhdadv_alg=centered.  Units are m**4/s.

      REAL(r8) :: nd_diff=0
      REAL(r8) :: nd_hypd=0

c     When either is nonzero, it is mathematically sound to choose
c	(homogeneous) Dirichlet conditions.  This input allows the
c	user to choose.

      CHARACTER(16) :: nd_bc="no flux"	! or "dirichlet"
c-----------------------------------------------------------------------
c     Viscous dissipation in the velocity equation:  
c	different viscous stress tensors may now be added together.
c
c	kinematic:     Pi(V) = -rho*kin_visc*grad(V)
c	isotropic:     Pi(V) = -rho*iso_visc*W
c	parallel:      Pi(V) = -rho*3*par_visc/2*(b.W.b)*(bb-I/3)
c	gyroviscosity: Pi(V) = -eta3/2*gyr_visc*
c				((bxW).(3bb+I)-(3bb+I).(Wxb))
c
c	where b is the magnetic direction vector, eta0 is the ion
c	pressure times the ion collision time, and W is the rate of
c	strain tensor,
c
c	  W=grad(V)+grad(V)^T-(2/3)*div(V)
c
c	and kin_visc, iso_visc, and par_visc have units of m**2/s.  
c	the diffusivity shaping specified below is used in both the
c	kinematic and the isotropic coefficients.

      REAL(r8) :: kin_visc=0
      REAL(r8) :: iso_visc=0
      REAL(r8) :: par_visc=0
      REAL(r8) :: gyr_visc=0
c-----------------------------------------------------------------------
c     Specify whether to time advance the tracking marker
c     which determines the plasma-vacuum interface (PVI)

      LOGICAL :: advance_pvi=.FALSE.
c-----------------------------------------------------------------------
c     Shaping parameters for both resistivity and viscosity:  
c	The character variable, ds_use, controls which diffusivities are
c	shaped by this function.  (see dvac and dexp above)

      CHARACTER(8) :: ds_use="both"	!or "elecd","kin_visc","neither"

c     This parameter is used in a sink term for velocity, so that
c	flows in the 'vacuum' region are damped.  When used, the
c	viscosity goes to zero in this region, instead of to a large
c	value, so that the plasma/vacuum interface does not exert a 
c	drag on the plasma.  The coefficient is zero in the core,
c	rising to vsink_rate in the vacuum.  Units are s**-1.  The
c	default is to not use the sink term and to use a large
c	viscosity instead.

      REAL(r8) :: vsink_rate=0		! no longer available.
c-----------------------------------------------------------------------
c     Apply a toroidal electric field to the n=0 'perturbation.'
c	This is used to self-consistently generate the 0-th order
c	part of the field for nonlinear runs when an equilibrium file
c	is not available.  Use be0 to set the vacuum toroidal field at
c	the center of the polar mesh and loop_volt to apply a loop
c	voltage (in Volts).  The 'equilibrium' arrays will remain 0,
c	except for the vacuum toroidal field, and the n=0 'perturbed'
c	part will serve as the equilibrium field - the vacuum toroidal B
c	+ the n=0 perturbation.  tloopv0 & tloopv1 are the
c	beginning and ending times to ramp the applied voltage from 0
c	to the full value of loop_volt.
c	CAUTION:  If the the equilibrium field has current, loop_volt
c	will effectively add to it.  To avoid equilibrium current when
c	gridshape='circ', ensure that lam0=0.

      REAL(r8) :: loop_volt=0
      REAL(r8) :: tloopv0=-1
      REAL(r8) :: tloopv1=0

c     The following allows the voltage to vary in order to achieve a 
c	desired amount of current.  The voltage is adjusted according
c	to V = V_old + dt*loop_rate*(i_desired-i_n0)
c	             - dt*loop_rate2*d(i_n0)/dt      ,
c	where i_n0 is just the perturbed part.  If the equilibrium
c	fields have current, i_desired is interpretted as the difference
c	between the actual desired current and the equilibrium field
c	current.  This is used instead of the voltage ramp whenever
c	i_desired is not 0.  Note that the di/dt part of the feedback
c	is not applied during the first time-step (leading to minor
c	discrepancies if a restart is inserted into a sequence of
c	time-steps), and V_old at the beginning of a simulation is
c	manually set with loop_volt.
c
c	Be aware that the behavior of the current
c	depends on the simulated plasma!

      REAL(r8) :: i_desired=0		!   in Amperes
      REAL(r8) :: loop_rate=1.e-4	!   in Ohms/s
      REAL(r8) :: loop_rate2=0		!   in Ohms

c     The next set applies a vertical electric field, possibly for
c	electrostatic current drive in compact devices.  The
c	t_e_vert0-&-1 are beginning and ending times for ramping the
c       electric field.

      REAL(r8) :: e_vertical=0		!   in Volts/m
      REAL(r8) :: t_e_vert0=-1
      REAL(r8) :: t_e_vert1=0
c-----------------------------------------------------------------------
c     The default behavior provides a normal component of flow equal
c	to the ExB drift when a loop voltage or other tangential
c	electric fields are applied.  However, since there is no surface
c	integral for mass flux, it can cause an unresolved boundary
c	layer in number density if continuity is not none.  Setting
c	norm_flow to "none" will set the normal component of flow to
c	zero regardless of tangential E.

      CHARACTER(16) :: norm_flow="exb"

c     We typically use no-slip boundary conditions, apart from a
c	possible specified surface normal from ExB.  The following
c	allows you to choose free-slip conditions independent of
c	viscosity.  However, no-slip should be used if equilibrium or
c	perturbed magnetic field penetrates the wall.

      CHARACTER(16) :: flow_bc="no-slip"   !   or "free-slip"
c-----------------------------------------------------------------------
c     Homogeneous Dirichlet boundary conditions are applied to the
c	change in the normal component of the magnetic field.  When
c	set to T, the zero_bnorm parameter will zero the entire normal
c	component of B, if it is nonzero in the initial conditions.

      LOGICAL :: zero_bnorm=.true.
c-----------------------------------------------------------------------
c     Add the force from a spatially uniform gravitational acceleration
c	for perturbations from the steady-state mass density.

      REAL(r8), DIMENSION(3) :: gravity=0.	!  (m/s**2)

c=======================================================================
c     other equilibrium  parameters.
c=======================================================================
c     When true, nimrod's initialization recomputes the parallel
c	component of the equilibrium current density at the numerical
c	quadrature points from the curl of Beq.  The perpendicular
c	component is also recompted at the quadrature points
c	to satisfy JeqXBeq=grad(Peq) exactly.  This option should not
c	be used if the force-balance relation includes other effects
c	such as flow.  Despite the name of the parameter, it
c	can be used for linear geometry, too.

      LOGICAL :: tor_eqja_fe=.false.
c-----------------------------------------------------------------------
c     For rectangular grid-shapes in toroidal geometry, and circular
c	or fluxgrid-based meshes in toroidal geometry, add magnetic
c	fields due to symmetric toroidal coils to the equilibrium
c	poloidal magnetic field.  The magnetic field components are
c	based on eqns. 5.40 of Jackson, including second-order terms.

      INTEGER(i4) :: ncoil=0
      REAL(r8), DIMENSION(ncoil_max) :: coil_current=0._r8 ! (Amps)
      REAL(r8), DIMENSION(ncoil_max) :: coil_r=0._r8       ! (meters)
      REAL(r8), DIMENSION(ncoil_max) :: coil_z=0._r8
c-----------------------------------------------------------------------
c     The nedge parameter has been used in certain applications, and it
c	is useful to declare it here.

c     In this version of the code, when nedge is not 0, it is used to
c	indicate a constant level of density to leave in the equilibrium
c	nd_eq data when performing a transfer_eq.  This mirrors nimeq's
c	ndeq_notx input.

      REAL(r8) :: nedge=0._r8
c=======================================================================
c     closure physics parameters.
c=======================================================================
c     The variable neoe_flag controls which closure model is
c       implemented in the Ohm's law. The current choices are "none",
c       "HandS", "Hole", "pp1", and "pp2".
c       The abbreviation "HandS" is for Hirshman and Sigmar
c       and uses the mu_e, mu_i, and ng_ft variables.
c       The "Hole" closure estimates the island width for the mode
c       resonant at q=m_neo/n_neo and does an analytic calculation
c       of the perturbed bootstrap current for that single harmonic.
c       This has only been tested with linear MHD.
c       The pp1 closure is a bootstrap current closure based
c       on making a simplification of the HandS closure where the
c       anisotropic pressure is assumed to be given only by the current
c       term and all equilibrium and gradients on the anisotropic
c       pressure are neglected.
c       The pp2 closure is similar to the pp1 closure, but makes the
c       assumption that the current term in the anisotropic pressure
c       evaluation is based on ideal MHD, j = B_eq x grad p1/beq**2.
c     The variable neoi_flag controls which closure model is
c       implemented in the ion equation. The current choices are "none"
c       and "HandS". The abbreviation "HandS" is for Hirshman and Sigmar
c       and produces poloidal flow damping.  The neoi_flag should always
c       be none when neoe_flag = "Hole".
c
c-TMP Neoclassical options are disabled in the present nimuw version.


      CHARACTER(8) :: neoe_flag="none"
      CHARACTER(8) :: neoi_flag="none"
      LOGICAL :: neo_debug=.false.
      INTEGER(i4) :: ng_ft=10   ! trapezoidal quadrature for trapped
				! particles.
      REAL(r8) :: mu_e=0._r8    ! electon viscous damping frequency
      REAL(r8) :: mu_i=0._r8    ! ion viscous damping frequency
      REAL(r8) :: neo_rad=0.0_r8 ! Smoothing radius near magnetic axis
      INTEGER(i4) :: m_neo = 2  ! poloidal mode number (Hole)
      INTEGER(i4) :: n_neo = 1  ! toroidal mode number (Hole)
      CHARACTER(8) :: neoe_det="all"
      CHARACTER(8) :: neoi_det="all"
c-----------------------------------------------------------------------
c     The variable p_model controls the type of pressure, i.e.
c	temperature, closure for total pressure (separate_pe=F) or
c	electron and ion pressures (separate_pe=T).  The available
c	options are:
c
c       1. 'adiabat' uses advective and compressive terms but no
c		thermal conduction.
c       2. 'isothermal' excludes the compressive term in the
c		temperature equation(s).
c	3. 'isotropic' adds isotropic thermal conduction to the 
c		temperature equation(s) with k_perpi coefficient
c		(and k_perpe for separate_pe=T)
c		This is like the former 'iso' option but renamed to
c		help distinguish the new isothermal option.
c	4. 'aniso1' uses anisotropic thermal conduction in the
c		temperature equation(s) with diffusivities k_perpi
c		and k_plle for single T and k_perpe and k_plli for
c		separate Ti and Te.
c       5. 'aniso_plltdep' nonlinear-only option for running
c               temperature-dependent parallel thermal conductivity,
c               k_plls->MIN( k_plls*(Ts/k_pll_ref_t)**2.5, k_pll_max )
c               where s = e or i (species).  Also, k_plls are
c               not allowed to drop below k_pll_min.
c       6. 'aniso_tdep' nonlinear-only option for running
c               temperature-dependent parallel thermal conductivity
c               (see 6.) and temperature- and mod(B)-dependent
c               perpendicular thermal conductivity with
c               k_perps-> k_perps*(k_pll_ref_t/<Ts>)**0.5*(1/<B>)**2
c               where s = e or i (species), <B> is in Tesla, and
c               <> stands for equilibrium + n=0 only.
c
c               To provide extra thermal diffusivity near boundaries
c               that represent a gap, the k_pll_ref_t/<Ts> ratio
c               has a lower bound of
c
c               1+(kprp_mnrat*(diff_shape-1)/dvac)**2
c
c               where diff_shape is the evalated diffusivity shape
c               profile determined by dvac and dexp.
c
c		Note: also see closure_model below.
c
c       Be aware that all diffusivities assume that the (gamma-1)
c       factor is absorbed into the input parameter, so that k_*
c       specifies a diffusivity coefficient for T.  [With dimensions
c	of m**2/s, all of the kappas and k_ps would be better labeled
c	as chis.]  k_cross is a switch for bXgrad(T) heat fluxes in
c       Braginskii; k_cross=1 gives the correct coefficient including
c       (gamma-1).  These terms are only used with separate_pe=.true.,

        CHARACTER(16) :: p_model="adiabat"
        REAL(r8) :: k_perpe=0
        REAL(r8) :: k_plle=0
        REAL(r8) :: k_perpi=0
        REAL(r8) :: k_plli=0
        REAL(r8) :: k_perp=0
        REAL(r8) :: k_pll=0
        REAL(r8) :: k_pll_max=1.e14
        REAL(r8) :: k_pll_min=0
        REAL(r8) :: k_pll_ref_t=1
        REAL(r8) :: kprp_mnrat=0
        REAL(r8) :: k_cross=0

c     Note that there are also some subtleties with the choice of
c	number density evolution.  The correct form of the temperature
c	equation has conductive heat flux proportional to n
c
c       n*dT/dt=-(gamma-1)*n*T*div(v)+div(n*chi.grad(T))
c
c	where d/dt includes v.grad() and chi is the temperature
c	diffusivity with gamma-1 absorbed.  For p_model='aniso1' the
c	temperature equation has this form and spatial variations in n
c	are used according to the continuity input parameter.

c     Additionally, if anisotropic conduction is used, we assume
c	B_eq.grad(T_eq)=0 and J_eq.grad(n_eq)=0 (separate_pe=T),
c	which are stronger assumptions on the equilibrium fields
c	than the usual steady-state assumption.  If this is not
c	satisfied, consider using transfer_eq to make the equilibrium
c	fields part of the initial condition only.
c-----------------------------------------------------------------------

c     THERMAL CONDUCTION CLOSURE MODEL

c     The variable closure_model determines which closure model is
c       used to calculate thermal transport coefficients when
c       "aniso_tdep".

c     The available options are:
c       1) "std kprp n0"  !  the limited, high-mag described above
c       2) "standard"     !  limited, high-mag w 3D perp conductivity
c       3) "braginskii"   !  Braginskii's magnetizing conductivity in 3D
c       4) "k2"           !  Ji's k2 model

      CHARACTER(12) :: closure_model="std kprp n0"

c     BRAGINSKII AND K2 CLOSURE MODELS:

c     The magnetization coefficients used to calculate
c       chi_perp have the form:  xs=omega_cs*tau_s

c     Unless artifically scaling magnetization coefficients,
c       leave magfac_ele=magfac_ion=1.

      REAL(r8) :: magfac_ele=1._r8
      REAL(r8) :: magfac_ion=1._r8

c     If tdep_coul_log=T, a temperature dependent coulomb logarithm
c       is computed for thermal diffusivity coefficients and
c       the thermal equilibration rate.  Otherwise, a constant value
c       of coulomb_logarithm is used for the calculations.

c     When using tdep_coul_log=T, leave coulomb_logarithm=1.

      LOGICAL :: tdep_coul_log=.false.
      REAL(r8) :: coulomb_logarithm=1._r8

c-----------------------------------------------------------------------

c     STANDARD MODEL

c     Thermal conduction coefficients have the form:
c       kappa_pll,s = k_plls*(Ts/k_pll_ref_t)**2.5
c       kappa_perp,s = k_perps*(k_pll_ref_t/Ts)**0.5/B**2

c-----------------------------------------------------------------------

c     BRAGINSKII MODEL

c     Reference:  Braginskii. "Transport Processes in a Plasma."
c       Reviews of Plasma Physics. 1965.

c     Thermal diffusivity coefficients have the form:
c       chi_pll,s  = Ts*tau_s/ms * (gamma0_s/delta0_s)
c       chi_perp,s = Ts*tau_s/ms *
c           (gamma1_s*xs**2+gamma0_s)/(xs**4+delta1_s*xs**2+delta0_s).
c       where tau_i=SQRT(2)*tau_ii and tau_e=tau_ei

c     Unless artificially scaling thermal diffusivity coefficients,
c       use k_plli=k_perpi=k_plle=k_perpe=1.

c     The default coefficients given below are for a Z = 1 plasma.
c       Ion coefficients for different Z plasma are given in the table
c       on page 25 of the reference.

      REAL(r8) :: gamma1_ele=4.6640_r8
      REAL(r8) :: gamma0_ele=11.920_r8
      REAL(r8) :: delta1_ele=14.790_r8
      REAL(r8) :: delta0_ele=3.7703_r8

      REAL(r8) :: gamma1_ion=2.0000_r8
      REAL(r8) :: gamma0_ion=2.6450_r8
      REAL(r8) :: delta1_ion=2.7000_r8
      REAL(r8) :: delta0_ion=0.6770_r8

c-----------------------------------------------------------------------

c     K2 MODEL

c     Reference:  Ji. "Closure and transport theory for
c       high-collisionality electron-ion plasma." 2011.

c     Thermal diffusivity coefficients have the form:
c       chi_pll,s  = T_s*tau_s/ms * f1(xs,zeta)
c       chi_perp,s = T_s*tau_s/ms * f2(xs,zeta)
c       where tau_s=tau_ss and zeta = sqrt(me/mi*ti/te)/Z

c     Unless artificially scaling thermal diffusivity coefficients,
c       use k_plli=k_perpi=k_plle=k_perpe=1.

c-----------------------------------------------------------------------
c     When separate_pe=.true., tequil_rate sets heat flux density
c	from electrons to ions (and vice versa) as
c
c	Q_i=-Q_e=n_e*k*(Te-Ti)*tequil_rate/(gamma-1)
c
c       and n_s*k*dT_s/dt equation has the additional term, Q_s.
c	Equilibration is not used if tequil_rate=0, the default.

      REAL(r8) :: tequil_rate=0

c     Specify whether to use a temperature-dependent thermal
c       equilibration rate.  When tdep_tequil=T, use tequil_rate=1.
c       Coefficients are calculated within the thermal_equil routine
c       within field_comps.f.  This option only applies to nonlinear
c	computations; linear computations use the fixed tequil_rate.

      LOGICAL :: tdep_tequil=.false.

c-----------------------------------------------------------------------
c     The following parameters determine whether dissipated energy
c	is used as a source in the temperature equations.  If
c	separate_pe=F, they would both add to the single-fluid
c	temperature evolution.  If separate_pe=T, Ohmic heating
c	(eta*J**2) is used in the electron temperature evolution only,
c	and viscous heating -grad(v)^T:Pi(v) is used in
c	the ion temperature evolution only.

      LOGICAL :: ohm_heat=.false.
      LOGICAL :: visc_heat=.false.
c-----------------------------------------------------------------------
c     This is a logical flag to set the flux of heat to zero at the
c	walls intead of Dirichlet boundary conditions on temperature
c	when thermal conduction is used.

      LOGICAL :: insulate=.false.
c-----------------------------------------------------------------------
c     Computation of the parallel viscosity coefficient is controlled
c	by parvisc_model.  Its default behavior uses par_visc as a
c	fixed viscous diffusivity.  When it is set to "plltdep," the
c	Braginskii T**5/2 temperature dependence is used.  Since the
c	coefficient for parallel thermal diffusivity has exactly the
c	same dependence, the physics kernel computes the temperature-
c	dependent thermal diffusivity as described above for option 5
c	of p_model (using k_pll_ref_t, k_pll_max, and k_pll_min) and
c	multiplies the result by par_visc/k_plli to get the viscous
c	diffusivity.  This plltdep option for viscosity may be used
c	with or without the temperature-dependent thermal diffusivity.
c
c	Like the temperature-dependent thermal diffusivities, this
c	option, this option is only for nonlinear computations.

      CHARACTER(16) :: parvisc_model="fixed"	! or "plltdep"
c=======================================================================
c     numerical parameters.
c=======================================================================
c     The following are general run-controlling parameters.
c	The cpu time is checked for each processor during parallel runs,
c	and it's not the total time across processors.  Note that
c	nim_stop is called at the end of the first time step that
c	crosses cpu_tmax, so allow a small margin when setting this
c	value and your batch limit.

      REAL(r8) :: dtm=1.e-8	! maximum time step, seconds
      REAL(r8) :: dt_initial=0	! initial time step if nonzero
      REAL(r8) :: dt_stop=1.e-4 ! code stops if dt<=dt_stop*dtm
      REAL(r8) :: tmax=10.0	! simulation time limit, seconds
      REAL(r8) :: cpu_tmax=1.e9	! approximate cpu time limit, seconds
      INTEGER(i4) :: nstep=2000 ! limit on number of time steps
      INTEGER(i4) :: npc=1	! rect_cir gridshape smoothing passes
c-----------------------------------------------------------------------
c     The next set of parameters influence the finite element
c	discretization. 

c	The degree of the polynomials used in the expansions for the
c	mapping and for the dependent variables are set with
c	poly_degree.

      INTEGER(i4) :: poly_degree=3
c-----------------------------------------------------------------------
c     The following parameter controls whether 1) pressure used in the
c	velocity advance is computed from the nT product at the nodes
c	of the FE expansion and interpolated to the quad points or 2)
c	n and T are interpolated then combined.

      CHARACTER(16) :: p_computation='at quads'  !  or 'at nodes'

c	ngr sets the number of Gaussian quadrature points used in each
c	direction of a rectangular cell.  The number of quad points
c	is set by the function:
c
c	  #_rblock_quad_pts_per_dir=ngr+poly_degree-1
c
c	so that it automatically increases with the polynomial
c	degree of the expansions for the dependent variables.  (see 
c	above).  ngr<2 may produce singular matrices.

c       Also note that 'exact' integration of quadratic nonlinearities,
c       each of polynomial degree p with a basis of degree p, requires
c       at least (3*p+1)/2 quadrature points in each direction with
c       Gaussian integration.  [If nq is the number of quad points,
c       Gaussian is exact for polynomials to degree 2*nq-1.  Also,
c       'exact' here means that variation of the Jacobian within an
c       element does not affect accuracy.]  Considering the relation
c       between nq and ngr, the minimum ngr for exact integration of
c       quadratic nonlinearities is (p+3)/2.

c	The integration_formula lets you choose the type of numerical
c	integration, either Gaussian (no end-point abscissas) or 
c	Lobatto (abscissas located at the GLL nodes if ngr=2, including 
c	end points).  However, degenerate rblocks and rblocks touching
c	the geometric axis always use Gaussian to avoid mass matrices
c	with zeros on the diagonal at end points where the Jacobian is
c	zero.  To get Lobatto in other blocks, set nxbl>1.

c	met_spl selects what form the partial derivatives of the grid
c	have within the rblocks (isoparametric, bilinear, or piecewise
c	constant).

      INTEGER(i4) :: ngr=2
      CHARACTER(5) :: met_spl='iso'	!  'iso', 'liner', or 'pcnst'
      CHARACTER(8) :: integration_formula='gaussian'  !  or 'lobatto'

c     Note:  the conform=F option and lumping options have been
c	disabled permanently.
c
c	When conform is true, the code uses conforming bilinear elements
c	for magnetic field, current density and velocity.

c	lump_b tells the code whether or not to lump the mass matrix
c       associated with dB/dt.

c	lump_all lumps the mass matrix for every equation.

      LOGICAL :: conform=.true. 
      LOGICAL :: lump_b=.false.
      LOGICAL :: lump_all=.false.
c-----------------------------------------------------------------------
c     Time-discretization parameters:

c       fom is now the time-centering for NIMITH's advance.  This is an
c	  old parameter for the original time-stepping that hasn't been
c	  used.  However, NIMITH has a theta-centered advance, and the
c	  original intent of fom is the theta parameter of an implicit
c	  single-step advance with arbitrary centering.
c       feta is the time centering on resistive diffusion.
c       fvsc is the centering on viscous dissipation.
c       fthc is the implicit centering for thermal conduction.

      REAL(r8) :: fom=0.5
      REAL(r8) :: feta=1
      REAL(r8) :: fvsc=1
      REAL(r8) :: fthc=1

c-----------------------------------------------------------------------
c     Time-centering parameters for advective terms:  fb_vxb is the
c	time-centering (fraction from the predictor step used in the
c	corrector step) of the b in vxb of Ohm's law. fv_vdgv, fn_vdgn,
c	and fp_vdgp are similar centerings for v in momentum advection,
c	n in continuity, and p in the pressure equation, respectively.
c	All of these parameters should be <=1 (there are possible
c	exceptions when there is dissipation), but the minimum value
c	depends on v_cfl, which limits the time step to
c	v_cfl*(cell size)/V.

c     1D planar analysis (no dissipation) demonstrates that numerical
c	stability associated with advection is satisfied when
c	1-2*f+f^2*cfl^2<=0, where f is the centering parameter.  Here
c	cfl is an effective cfl [v*dt/(length scale)].  For the Fourier
c	direction, 1/(length scale) -> k_max.  For the poloidal plane,
c	the distributed mass matrix introduces an extra factor of
c	sqrt(3): 1/(length scale) -> sqrt(3)/(cell size).  These
c	factors are now computed within nimrod for the appropriate
c	components.

      REAL(r8) :: fb_vxb=0.54
      REAL(r8) :: fv_vdgv=0.54
      REAL(r8) :: fp_vdgp=0.54
      REAL(r8) :: fn_vdgn=0.54
      REAL(r8) :: v_cfl=0.5

c     The implicit leapfrog algorithm for HMHD requires time-centered
c	implicit advection, and it is always used for ohms/='mhd'.
c       However, the following allows a choice of predictor/corrector
c	advection (with the coefficients defined just above) or 
c	time-centered implicit advection.  This allows us to 
c	reproduce earlier computations with P/C advection, but it leads
c	to additional coding complications and may be removed at some
c	point.

      CHARACTER(8) :: mhdadv_alg='precor'	! or 'centered'

c-----------------------------------------------------------------------
c     The semi-implicit operator in the present formulation should be
c	effective for stabilizing nonlinear activity.   However, in
c       some cases, an additional time-step limitation may be useful.
c       nl_cfl_lim limits the time step so that a cfl condition based
c	on wave speeds computed with nonlinear pressures is satisfied.
c	The default is effectively no limit.

      REAL(r8) :: nl_cfl_lim=1.e10

c-----------------------------------------------------------------------
c     Select where the semi-implicit operators appear.  If si_in_v is
c	set to true, the mhd operator appears in the velocity advance.
c	If not, one operator appears in the magnetic field advance,
c	and another appears in the pressure advance.  The hall advance
c	is separate with its own operator.  At present, the si operator
c	for v is isotropic, and mhd_si_iso only affects the si_in_v=F
c	case.

      LOGICAL :: si_in_v=.true.		! option no longer functions;
					! si_in_v is always true.

c     Nonlinear velocity advances may use either the standard semi-
c	implicit operator that is described in the 2004 JCP paper or a
c	3D version.  The standard operator acts separately on each
c	Fourier component, which is convenient for the linear algebra.
c	The 3D operator more accurately represents toroidal deformations
c	and improves accuracy at large time-step when toroidal variation
c	is substantial.

      CHARACTER(8) :: siop_type='standard'    !   or '3D'

c     Semi-implicit coefficients:  In general, unity gives the minimum
c	coefficient based on linear numerical stability analysis.  
c	Factors>1 are necessary to stabilize nonlinear terms.
c	si_fac_mhd is a coefficient for the semi-implicit
c	operator used for the vxb term in Faraday's law.  si_fac_pres is
c	the coefficient for the pressure advance.  si_fac_j0 is a 
c	coefficient for the parallel and perpendicular j0 terms in the
c       velocity equations (j0Xb_tilde and v.grad(P0)).  si_fac_nl is
c	a coefficient for an isotropic operator in velocity, whose
c	magnitude is based on the nonlinear pressures.

      REAL(r8) :: si_fac_mhd=1
      REAL(r8) :: si_fac_pres=1
      REAL(r8) :: si_fac_j0=1
      REAL(r8) :: si_fac_nl=1.1

c     At present, si_fac_hall is only used in the preconditioner for
c	nonlinear computations with ohms='2fl' and multiplies the 
c	partial(J)/partial(t) term to increase diagonal dominance when
c	its value is greater than unity.  It should not affect the
c	numerical result provided that the solver tolerance is
c	sufficiently small.

      REAL(r8) :: si_fac_hall=2.

c     The parameter mhd_si_iso is used to control the degree of isotropy
c	for the linear semi-implicit mhd operator.  It must satisfy
c	0<=mhd_si_iso<=1.  The operator is fully isotropic
c       when this factor is 1 and fully anisotropic when this factor is
c	zero.  This adds to the nonlinear isotropic operator controlled
c	by si_fac_nl.

      REAL(r8) :: mhd_si_iso=0.
c-----------------------------------------------------------------------
c     Split resistive diffusion from the mhd part of the magnetic field
c	advance?  This may improve numerical stability at larger time
c	steps at the possible expense of the order of convergence in dt.

      LOGICAL :: split_resist=.false.	!  this option is no longer
c					!  functioning.  always = F.

c     Split the viscous forces from and advective & jxb in the velocity
c	advance?  This will improve numerical stability at larger time
c	steps at the possible expense of the order of convergence in dt.

      LOGICAL :: split_visc=.false.
c-----------------------------------------------------------------------
c     Divergence of B cleaning:  divbd is the diffusivity used with
c	the grad(div(B)) term, and fdivb is the time-centering parameter
c	for the implicit solve.  ndivb is the cycle frequency for
c	calling the divergence diffusing equation.  If split_divb is
c	true, the divb cleaning will be done in a separate time split
c	after the rest of the time advance.  Otherwise it is applied
c	in the magnetic-field advance.
c
c	[ndivb is only applicable when split_divb=true; if not,
c	divergence cleaning is applied every time step.]

      REAL(r8) :: divbd=0
      REAL(r8) :: fdivb=1
      INTEGER(i4) :: ndivb=1
      LOGICAL :: split_divb=.false.

c     Limit on the tolerated wavenumber**2.  This should be set to
c	1/L**2, approximately, where L is the length scale of the
c	domain.

      REAL(r8) :: kdivb_2_limit=1

c     This version of nimrod allows a separate discontinuous expansion
c	for cleaning div(B), employing a mixed finite-element method
c	with parabolic correction.  The auxiliary
c	discontinuous scalar is of polynomial degree poly_divb and is
c	allocated if the parameter is set >=0.
c
c	The poly_divb_min and poly_divb_max parameters are used
c	to indicate that the modal expansion for the auxiliary field
c	is incomplete.  Where these limits are used, the bases are
c	the outer product of a set of Legendre polynomials from 0
c	to poly_divb in the x logical coordinate and the Legendre
c	polynomials from poly_divb_min to poly_divb_max in the
c	y logical coordinate.  The union of these 2D polynomials and
c	ones that have x and y swapped defines the incomplete expansion.
c
c	The parabolic correction method solves the system,
c	  dB/dt = -curl(E) + disc_dbd*grad(auxb)
c	   auxb = div(B)
c
c	The method uses fdivb for implicit centering.

      REAL(r8) :: disc_dbd=0._r8   !  in L**2/time
      INTEGER(i4) :: poly_divb=-1  !  set to poly_degree-1 if used
      INTEGER(i4) :: poly_divb_min=-1  !  best left at -1 for a
      INTEGER(i4) :: poly_divb_max=-1  !  complete lower-order expansion
c-----------------------------------------------------------------------
c     When hyp_eta is greater than zero, a numerical hyper-resistivity
c	is applied in the advance of magnetic field.  The coefficient
c	is uniform and time-independent.  It has units of L**4/time.
c	Note that the magnetic advance requires solution of a 6-vector,
c	so it takes much more time per solve.

c	The split_hypeta option allows a choice of time-splitting the
c	hyper-resistive diffusion as a separate solve.

      REAL(r8) :: hyp_eta=0._r8
      REAL(r8) :: fhyp_eta=1._r8  !  temporal centering parameter
      LOGICAL :: split_hypeta=.false.

c     The hyp_dbd is a hyper-diffusivity for divergence cleaning that
c	is similar to the hyper-resistivity, except that the operator
c	is grad(div(grad(div()))).  The split_hypeta input affects this
c	computation, too.

      REAL(r8) :: hyp_dbd=0._r8
      REAL(r8) :: fhyp_dbd=1._r8

c     Note that dissipation from hyper-resistivity and div(B) hyper-
c	diffusivity is added to ohmic heating in nonlinear problems
c	with ohm_heat=.true.  The computation assumes that the two
c	hyper terms have the same centering parameter, so fhyp_dbd
c	is reset in these cases.
c-----------------------------------------------------------------------
c     Time-split hyper-viscosity can be applied to the center-of-mass
c	flow velocity.  The hyp_visc parameter is the hyper diffusivity
c	value [L**4/t] for V.  fhyp_visc is centering parameter for
c	the implicit time-split advance, but there will be first-order
c	errors, regardless of its value.

      REAL(r8) :: hyp_visc=0._r8
      REAL(r8) :: fhyp_visc=1._r8
c-----------------------------------------------------------------------
c     Computations can use auxiliary scalar fields with incomplete
c	discontinuous representations to alter convergence properties.
c	The two fields described here work like the correction method
c	for magnetic divergence described above, but they act on
c	the momentum-density evolution.  One helps stabilize flow
c	divergence, and the other helps stabilize parallel vorticity.
c
c	When used in a computation, the modified system is (loosely)
c
c	   rho*dV/dt = F + sqrt[ddivv*(B^2/mu0+gamma*P)*dt]*grad(auxv)
c	                 + sqrt[dpvrt*dt/mu0]*BXgrad(auxw)
c	   auxv = sqrt[ddivv*(B^2/mu0+gamma*P)*dt]*div(V)
c	   auxw = sqrt[dpvrt*dt/mu0]*B.curl(V)
c
c	where ddivv is a normalized viscosity for divergence error,
c	and dpvrt is a normalized viscosity for vorticity error.
c	This description is not precise, because the implementation is
c	defined directly in weak form and not derived from the strong
c	form shown above.
c
c	The representations need to be incomplete so that they do not
c	duplicate the physics responses to compression
c	and parallel vorticity in numerically resolved evolution.  
c	Thus, these auxiliary fields are expanded in incomplete modal,
c	Legendre polynomial representations.  For each of the two
c	logical coordinates in a 2D element, the polynomials from
c	poly_divv_min to poly_divv_max in one coordinate and from 0
c	to poly_divv in the other (and vice versa).  The same holds
c	for the auxiliary field for parallel-vorticity.
c
c	As described in the Aug. 2013 team meeting presentation on
c	projection, the most useful incomplete expansion has
c	poly_divv=poly_divv_min=poly_divv_max=poly_degree.  If the
c	poly_divv_auto input is set to true, the three separate
c	parameters are set to poly_degree automatically.
c
c	Speeds for the artificial waves can be comparable to
c	characteristic perpendicular and shear waves, respectively.
c
c	fdivv and fpvrt are implicit centering parameters for these
c	terms.

      REAL(r8) :: ddivv=1._r8
      REAL(r8) :: dpvrt=0.1_r8
      REAL(r8) :: fdivv=1.0_r8
      REAL(r8) :: fpvrt=1.0_r8
      INTEGER(i4) :: poly_divv=-1
      INTEGER(i4) :: poly_divv_min=-1
      INTEGER(i4) :: poly_divv_max=-1
      LOGICAL :: poly_divv_auto=.false.

c	NOTE: it is assumed that the two auxiliary fields have the
c	same expansion, so they are combined into a single two-vector,
c	which simplifies coding.
c-----------------------------------------------------------------------
c     The following parameters affect how frequently matrices and their
c	preconditioning factors are computed in nonlinear cases.  first,
c	nearly all of the matrices depend on dt.  thus, dt_change_frac
c	is used to force the time step to decrease by this fraction when
c	CFL conditions are limiting it.  it also prevents increases
c	until a change of dt_change_frac is allowable by the CFL
c	conditions.  this prevents time step from changing every cycle.

c     The n_dt_release parameter is a release from waiting for the CFL-
c	allowable dt to go back up by dt_change_frac.  After 
c	n_dt_release cycles, it is allowed to increase, anyway.  This
c	had been a parameter in subroutine new_dt.

c     When the time-step is allowed to increase, it is limited to
c	changing by a factor of dt_incr.  This helps keep the time-step
c	from repeatedly bumping against the CFL limit with subsequent
c	drops by a fraction of dt_change_frac.

c     The semi-implicit operators used for advancing B and P depend on
c	the equilibrium fields plus the n=0 part of the perturbed
c	solution.  these matrices are recomputed if the fractional
c	change in the sum of the two changes by more than
c	ave_change_limit.

c     For nonlinear computations with slowly changing fields, it is
c	useful to regularly update matrices regardless of the tests
c	on how much the fields have changed.  n_mat_update sets this
c	frequency in number of steps.

	REAL(r8) :: dt_change_frac=0.1
	REAL(r8) :: ave_change_limit=0.01
        REAL(r8) :: dt_incr=1.07
        INTEGER(i4) :: n_dt_release=10
        INTEGER(i4) :: n_mat_update=100000
c-----------------------------------------------------------------------
c     The following parameters enhance local particle diffusion in
c       regions where number density drops below a 'floor' value.  The
c       extra diffusivity is 
c
c         nd_dart_fac*(1/dt)*[cross_section/(mx*my)]*
c           MAX(0._r8,TANH( (nd_floor*ndens-n)/(nd_exp*ndens) ) )
c
c       where the factor in brackets is an average element area.  The
c       computation uses the minimum n over toroidal angle, and the
c       diffusivity is just a function of position in the poloidal
c       plane.
c
c       The input parameters are:
c
c       nd_floor    -- the floor value is nd_floor*ndens
c       nd_exp      -- nd_exp*ndens determines the transition 'width'
c                      of the diffusivity enhancement in terms of n
c       nd_dart_fac -- sets the magnitude of the diffusivity enhancement
c                      in terms of cell area per dt.
c
c       The default nd_floor (=0) skips this enhancement, and it's only
c       active when nd_diff>0.
 
      REAL(r8) :: nd_floor=0._r8
      REAL(r8) :: nd_exp=0.02_r8
      REAL(r8) :: nd_dart_fac=1._r8

c     This is a second artificial diffusivity that reproduces an
c       upwinding-like smoothing.  The coefficient nd_dart_upw is
c       dimensionless, and it is multiplied by
c
c       (jac/dt)*[ (dt*V.grad(n)/n_tot)**2 +
c                  MAX(0,(nd_floor_upw-n_tot)/nd_width_upw)) ] , 
c
c       where jac is the local 2D jacobian for the poloidal plane.  This
c       makes an effective diffusivity of jac/dt that becomes active
c       in places where (dt*V/L_n)**2 is significant and where
c       n_tot drops below nd_floor_upw.  The diffusion is parallel to
c       the flow only with upw_aniso=1 (see below), and it is only
c       applied to nonlinear cases with implicit advection.  Here, the
c       diffusivity is fully 3D.

      REAL(r8) :: nd_dart_upw=0._r8
      REAL(r8) :: nd_floor_upw=0._r8   !  [m**-3 unlike nd_floor]
      REAL(r8) :: nd_width_upw=1._r8   !  [m**-3 unlike nd_exp]
c-----------------------------------------------------------------------
c     Artificial thermal diffusivities for an upwinding-like smoothing
c       are also available for the temperature advance.  The coefficient
c       t_dart_upw is used in the same way as nd_dart_upw for number
c       density.  Here, the diffusivities
c
c       (jac/dt)*[ (dt*V.grad(T)/T_tot)**2 +
c                  MAX(0,(t_floor_upw-T_tot)/t_width_upw)) ] , 
c
c       use the species T and V if separate_pe is true; though, the VV
c       dyad used to create the preconditioner matrix is formed from
c       the COM flow only.  Again, this is only used in nonlinear
c       computations with implicit advection, and the diffusivity is 3D.

      REAL(r8) :: t_dart_upw=0._r8
      REAL(r8) :: t_floor_upw=0._r8   !  in eV
      REAL(r8) :: t_width_upw=1._r8   !  in eV
c-----------------------------------------------------------------------
c     The following parameters alter the upwinding-like diffusivities
c	for number density, temperature, etc.  upw_aniso controls
c	whether the diffusion is purely along the direction of flow.
c	It should have values ranging from 0 to 1 with 0 being
c	isotropic, and 1 being fully anisotropic.  upw_limit sets
c	a limit on (dt*V.grad(n)/n)**2 and (dt*V.grad(T)/T)**2.  Use a
c	value less than 1 and increase nd_dart_upw and t_dart_upw to
c	enlarge the region over which the artificial diffusivity is
c	significant.

      REAL(r8) :: upw_aniso=1.
      REAL(r8) :: upw_limit=1.
c-----------------------------------------------------------------------
c     Either of the artificial particle diffusivities (nd_diff and
c	nd_hypd) leads to violation of conservation of momentum and
c	energy.  The following option subtracts V*S_n from the flow-
c	velocity equation and adds (mV^2/2-k_B*T/gamma-1)*Sn to the
c	temperature equation to correct these errors.  The effective
c	source/sink factor is
c	  Sn=div(nd_diff*grad(n)-nd_hypd*grad(grad^2(n))).

      LOGICAL :: nd_correrr=.false.
c-----------------------------------------------------------------------
c     The following parameters set minimum or "floor" values directly
c	on the values of density and temperature at the nodal locations
c	of the spectral-element/Fourier expansions.  This is only used
c	in nonlinear computations.  Unlike diffusion, this approach is
c	not conservative and should only be used to control the density
c	and/or internal energy where they are negligibly small, such
c	as artificial cold-plasma regions.
c
c	Note that dealiasing, and dropping the largest-possible Fourier
c	coefficient when dealiase=F, make our FFTs non-invertible.
c	Thus, a specified minimum and what will be realized after the 
c	forward FFT will differ.  A test can be used to compare the
c	two (change the check_min parameter in set_nodal_min in
c	utilities.f), but it performs an extra FFT.
c
c       This option is only used when the paramters are set>=0.
 
      REAL(r8) :: nd_nodal_floor=-1._r8  !  in units of 1/m^3
      REAL(r8) :: ti_nodal_floor=-1._r8  !  in units of eV
      REAL(r8) :: te_nodal_floor=-1._r8  !  used when separate_pe=T
c-----------------------------------------------------------------------
c     Damping can be applied directly to the highest Fourier
c	coefficients of a nonlinear computation.  This is intended for
c	computations with dealiase=.false., but it also functions when
c	dealiase is true.  The damping rate for the largest retained
c	Fourier component for all physical fields is specified by
c	fmx_drate.  Components with n > nmax - nfdamp are damped at
c	proportionally larger rates up to fmx_drate:
c
c	  drate(n) = fmx_rate*MAX(0,(n-nmax+nfdamp)/nfdamp)
c
c       Fourier components are then just multiplied
c	F_n->F_n*(1-drate(n)*dt) at the end of each step.

      REAL(r8) :: fmx_drate=0._r8
      INTEGER(i4) :: nfdamp=0  !  default does not apply damping.
c-----------------------------------------------------------------------
c     For simulations that have vertices along r=0, regularity
c	conditions must be satisfied for the different Fourier
c	components.  One of these conditions is that d/dr of any
c	radial or toroidal vector component of n=1 should be 0.  This
c	is not satisfied through usual finite element natural boundary
c	conditions, since the effective surface area is 0.  Instead,
c	a volume constraint integral is added to the matrices to 
c	inhibit the growth of d/dr for these components.  The
c	coefficient r0dr_weight is a factor for the integrand, which
c	is proportional to r**2*(d/dr)**2 along cells touching r=0.
c	The integral is also scaled to diagonal elements of each matrix
c	to which it is added.
c
c	If these components appear to violate the Neumann condition
c	excessively in a particular simulation, increase r0dr_weight.
c	Conversely, if the conditions distort the interior of the
c	solution, decrease r0dr_weight.  In any case, the influence
c	of the constraint integral decreases with increasing radial
c	resolution.

	REAL(r8) :: r0dr_weight=1	!  NO LONGER OPERATIONAL
c=======================================================================
c     algebraic solver parameters.
c=======================================================================
c     The "solver" parameters affect the choice of preconditioner that
c	is used in NIMROD's Krlov-space iterative solves.
c
c	Existing selections include:
c	1) "diagonal"  where the inverted diagonal elements are the
c		sub-matrices consisting of all quantities located at
c		each grid vertex, e.g., the real and imaginary vector
c		components of B at each vertex for the B-field solve.
c	2) "lapack"  a banded serial direct solve over all grid-blocks
c		(rblocks only, tblocks use diagonal) using Lapack
c		routines.  Here, the cg is only used for iterative
c		refinement on 2D matrices. The old "bl_drect" option
c		now reverts to "lapack".
c	3) "bl_ilu_n"  where n is "0" or "1".  this option was not
c		updated for high-order elements and has been removed.
c	4) "bl_diagx" or "bl_diagy" this is a line-Jacobi iteration,
c		where a direct solve is computed along either the x-
c		or y-direction.  It's mainly intended for single-block
c		problems on vector machines.  Another related choice
c		is "bl_diaga", which is an average of the separate
c		x- and y- direction solves.
c	5) "gl_diaga" global line-Jacobi preconditioning, where the
c		lines extend over contiguous, conforming rblocks to the
c		greatest extent possible.  Parallel decomposition
c		swapping is performed with asynchronous mpi calls.
c		This is probably the best large-problem, multi-block
c		preconditioner that we have at this point.
c	6) "seq_slu"  a sparse serial direct solve over all rblocks
c		using the Sequential SuperLU library, which needs to
c		be linked at the load step.  As with "lapack,"
c		the cg is only used for iterative refinement on 2D mats.
c	7) "slu_dist"  a sparse parallel direct solve over all rblocks
c		using the SuperLU_DIST library, which needs to
c		be linked at the load step.  As with "lapack,"
c		the cg is only used for iterative refinement on 2D mats.
c	8) "slu_dstm"  is the same as "slu_dist" except that the library
c		is called through a distributed memory interface and
c		does more communication during the factorization and
c		substitution steps.  it may be somewhat slower while
c		using less memory. 
c	9) "slu_dsta"  is the same as "slu_dsta" with respect to the
c		slu library, but is uses point-to-point communication
c		when setting-up the sparsity pattern to avoid increasing
c		memory requirements with increased parallelization.

      CHARACTER(8) :: solver="diagonal"
 
c     The following parameters allow separate choices for specific
c       equations. Any parameters left at their default value will
c       be set to the value of the solver input parameter.
 
      CHARACTER(8) :: vmhd_solver="none"
      CHARACTER(8) :: bmhd_solver="none"
      CHARACTER(8) :: temp_solver="none"
c-----------------------------------------------------------------------
c     The parameter maxit is the maximum number of CG or GMRES
c	iterations, and tol is the relative convergence tolerance.

c     NOTE memory considerations: The basic resistive-MHD advance with
c	predictor-corrector or no advection uses the conjugate-gradients
c	method to solve the resulting Hermitian matrices.  The amount of
c	memory required is just determined by the mesh in those
c	computations.  However, many options require solution of non-
c	Hermitian matrices, and GMRES is used to solve those matrices.
c	In all of these cases, the solver needs to keep all direction
c	vectors, so the maxit parameter is used to allocate space.
c	Thus, very large values of maxit can cause computations to
c	crash, as a result of memory limitations, when invoking the
c	GMRES solver.  The options that lead to non-Hermitian matrices
c	include mhdadv_alg, ohms, gyr_visc, hyper-diffusion terms,
c	and numerical methods that project onto discontinuous bases.

      INTEGER(i4) :: maxit=50
      REAL(r8) :: tol=1.e-8
c-----------------------------------------------------------------------
c     The existing line-solve based (like bl_diagx
c	and bl_adi_y) preconditioners do not use pivoting and are
c	unstable for matrices with large condition numbers.
c	off_diag_fac is used to improve the condition number (of the
c	factored matrix	only).  The off-diagonal sub-matrices are
c	multiplied by this factor prior to the incomplete factorization.
c	If you run into	problems with the factorization, reduce
c	off_diag_fac.  If you are adventurous, try increasing it.
c	[Line factorizations now have a loop to find the
c	largest factor possible.  off_diag_fac is now used as a lower
c	limit.]

      REAL(r8) :: off_diag_fac=0.9
c-----------------------------------------------------------------------
c     Old solution vectors from each equation are saved to extrapolate
c	guesses.  The extrapolation is based on a polynomial that
c	passes through the old data.  extrap_order is the order of this
c	polynomial.  As its value is increased, there is potentially 
c       more accuracy in the guesses at the expense of more memory.

      INTEGER(i4) :: extrap_order=1
c-----------------------------------------------------------------------
c     The nsym_pre preconditioning options have been removed from nimuw
c	starting after version 3_4_13.  Setting nsym_pre_band>0 will
c	not cause the code to stop, however.
c
c     When nsym_pre_band is greater than zero, nimrod generates
c	bands representing Fourier-component coupling above and below
c	the diagonal for use in block SOR preconditioning during the
c	nonsymmetric 3d solves.  The nsym_pre_rpass sets the number of
c	passes, which use alternating (up and down in F-comp) direction
c	between even and odd numbered iterations.  The relaxation factor
c	is nsym_pre_rfac, and nsym_pre_rtype indicates whether to use
c	Jacobi- or Gauss-Seidel-style passes.

      INTEGER(i4) :: nsym_pre_band=0
      INTEGER(i4) :: nsym_pre_rpass=2
      REAL(r8) :: nsym_pre_rfac=2._r8/3._r8
      CHARACTER(12) :: nsym_pre_rtype="Gauss Seidel"  !  or "Jacobi"

c     Note that the total number of GMRES iterations (the outer loop
c	over these preconditioning passes) is independent of layer
c	decomposition only when using Jacobi-style passes.  When using
c	Gauss-Seidel-style passes, only Fourier couplings within a layer
c	use the most current data; couplings outside the layer use data
c	from the previous iteration, like Jacobi.
c-----------------------------------------------------------------------
c     The following parameters control convergence on nonlinear
c	algebraic systems in the implicit leapfrog algorithm.  In
c	particular, the delta_V.grad(delta_V) term can be iterated in
c	the COM velocity advance, and the delta_Jxdelta_B nonlinear
c	Hall term can be iterated in the two-fluid B-advance.

c	The maxit_nl parameter is the maximum number of iterations
c	allowed, but setting it to 1 (the default) is a switch to
c	avoid nonlinear iteration.  The tol_nl specifies the relative
c	tolerance of the 2-norm of the residual.  Note that we
c	must have tol_nl>>tol.

      INTEGER(i4) :: maxit_nl=1
      REAL(r8) :: tol_nl=1.e-4
c=======================================================================
c     io  parameters.
c=======================================================================
c     These parameters control the restart dumps read and written by
c	nimrod.  The run is initiated by reading dump_file, which was
c	written by nimset or during the course of another nimrod run.
c	It writes dumps into the dump_dir directory with the name prefix
c	dump_name.  ndump controls the dumping frequency in terms of
c	number of time steps, and dump_over tells nimrod how to handle
c	pre-existing dump files with the same name.

      CHARACTER(128) :: dump_file="dump.00000"
      CHARACTER(64) :: dump_dir="."
      CHARACTER(64) :: dump_name="dump"
      INTEGER(i4) :: ndump=10000
      INTEGER(i4) :: dump_over=0 ! 0: overwrite; 1: append; 2: error
c-----------------------------------------------------------------------
c     Parameters for time history plots, created if nhist>0.  This
c	includes the single-point probe, and spatially-integrated
c	energy and div(b) diagnostics.

c	The hist_flush character indicates how frequently the buffers
c	for the history files are flushed to disk by closing and
c	reopening them.  Previously, files were closed after each write
c	to ensure that data was not lost in case of a crash, but this
c	can affect performance.  Options are:
c	  1) "always" -- close after each write (good for debugging)
c	  2) "at dumps" -- close on when writing dumps (a compromise)
c	  3) "at end" -- only close at the end of the run (fastest)

      INTEGER(i4) :: nhist=0	! time-history stride
      LOGICAL :: hist_binary = .true.	! true=binary file, false=text
      INTEGER(i4) :: ihist=0	! i-index on init grid for time hists.
      INTEGER(i4) :: jhist=0	! j-index on init grid for time hists.
      CHARACTER(8) :: hist_flush="at dumps"
c-----------------------------------------------------------------------
c     Parameters for xdraw plots--created if the *_stride 
c	values are >0.  Note:  these should be set to zero for running
c	in parallel, since the associated routines do not have 
c	processor-to-processsor communication.  Intead, use nimplot to
c	create the same files from dump files.

      CHARACTER(64) :: xdraw_dir="."	!   xdraw file directory
      INTEGER(i4) :: xy_stride=0 ! time stride for xy slice
      INTEGER(i4) :: xt_stride=0 ! time stride for xt slice
      INTEGER(i4) :: yt_stride=0 ! time stride for yt slice
      REAL(r8) :: x0fac=0.25    ! relative position of yt slice
      REAL(r8) :: y0fac=0.25    ! relative position of xt slice
c-----------------------------------------------------------------------
c     Dirctory name and time-step cycle frequency for IBM Data 
c	Explorer output.  Note:  set these to zero for running in
c	parallel, too (see above).

      CHARACTER(64) :: dx_dir="."	! DX file directory
      INTEGER(i4) :: ndxout=0		! stride for DX output
c-----------------------------------------------------------------------
c     Debugging output for nimset.

      LOGICAL :: detflag=.false.	! print detflaged ascii output?
c-----------------------------------------------------------------------
c     Long listing from iterative solver.

      LOGICAL :: itflag=.false.
c-----------------------------------------------------------------------
c     If reset_file is not 'none' when running nimset, the equilibrium
c	fields and perturbations are first initialized in the normal
c	way.  Then, nimset looks for reset_file and reads the 
c	perturbed fields contained therein.  If the reset_file has
c	at least as many Fourier components as the new case, all
c	components in the new case will be over-written by those from
c	reset_file.  If reset_file has fewer components (n), only the
c	first n will be over-written.  Note that there is no checking
c	for consistency between the wavenumber arrays of the new and
c	previous cases, and the equilibrium fields are not reset from
c	the reset file.  This is intended to facilitate a change in the
c	number of components for nonlinear runs and for changing the
c	the equilibrium in a linear parameter scan while keeping the
c	last eigenfunction.

c	reset_time and reset_step set the physical time and the
c	numerical time-step index in the new file.

c	setting reset_eq_mesh to true changes the way that a reset
c	works in that the mesh and equilibrium fields are read and
c	retained from the reset file.  the expansion or reduction
c	in the number of Fourier components for the solution fields
c	works as described above.

      CHARACTER(128) :: reset_file="none"
      REAL(r8) :: reset_time=0
      INTEGER(i4) :: reset_step=0
      LOGICAL :: reset_eq_mesh=.false.
c-----------------------------------------------------------------------
c     Some of the diagnostic output needs to represent fields that are
c	discontinuous across element boundaries.  They are written
c	directly to data files for display by external software
c	packages.  The following sets the format of these files.

      CHARACTER(8) :: elout_format="vtk"  !  or "tecplot"

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE input

c-----------------------------------------------------------------------
c     subprogram 2. read_namelist
c     open and read the namelist input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_namelist(infile,echo_in)
      USE local
      USE input
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: infile
      LOGICAL, INTENT(IN) :: echo_in

      INTEGER :: read_stat,nc
      INTEGER(i4) :: number_of_namelists=0,il,inst
      INTEGER(i4), DIMENSION(20) :: namel_order=0
      CHARACTER(14) :: tempfile
      CHARACTER(128) :: ctmp
      LOGICAL :: reading

c=======================================================================
c     the following variables are used by the central repository, also
c     known as "the developer's version."  these variables are not
c     functional in nimuw, but their declarations here allows this
c     routine to read different input formats.
c=======================================================================
c-----------------------------------------------------------------------
c     extra grid variables:
c-----------------------------------------------------------------------
      REAL(r8) :: vac_frac=0
      REAL(r8) :: x1=0,y1=0
c-----------------------------------------------------------------------
c     extra init variables:
c-----------------------------------------------------------------------
      REAL(r8) :: reset_tol=0,rescale_factor=0
      REAL(r8) :: init_x=0,init_y=0
      REAL(r8) :: blob_amp=0,blob_width=0
c-----------------------------------------------------------------------
c     extra equil variables:
c-----------------------------------------------------------------------
      REAL(r8) :: zeff=0    ! local declaration does not affect physdat,
			    ! but see zeff_input.
      REAL(r8) :: npp=0
      REAL(r8) :: te_cent=0,te_edge=0,tepp=0
      REAL(r8) :: divbeqd=0
      CHARACTER(8) :: n_profile='none',pe_profile='none',
     $                ds_function='none'
c-----------------------------------------------------------------------
c     extra physics variables:
c-----------------------------------------------------------------------
      INTEGER(i4) :: hv_power=0,hv_filter=0
      REAL(r8) :: hv_coeff=0,hv_time=0
      LOGICAL :: rt_transfer_eq=.false.
c-----------------------------------------------------------------------
c     extra closure variables:
c-----------------------------------------------------------------------
      REAL(r8) :: gamma_heat=0
      CHARACTER(16) :: heat_source='none',electric_source='none',
     $                 heat_ion_source='none',heat_ele_source='none',
     $                 momentum_source='none'
      CHARACTER(12) :: cel_model='none'
c-----------------------------------------------------------------------
c     extra numerical variables:
c-----------------------------------------------------------------------
      INTEGER(i4) :: nstep_converge=0
      REAL(r8) :: growth_eps=0
c-----------------------------------------------------------------------
c     extra solver variables:
c-----------------------------------------------------------------------
      INTEGER(i4) :: n_gm=0
c-----------------------------------------------------------------------
c     extra output variables:
c-----------------------------------------------------------------------
      CHARACTER(128) :: dump_list='none'
c-----------------------------------------------------------------------
c     all particle variables:
c-----------------------------------------------------------------------
      REAL(r8) :: avu0=0,avw0=0
      REAL(r8) :: mass0=0,charge0=0
      REAL(r8) :: vhmx=0,psipmax=0
      REAL(r8) :: Ppsi0=0,R0=0,betafrac=0
      INTEGER(i4) :: nm=0,restart=0
      CHARACTER(100) :: tpartdr='none'
      LOGICAL :: trace=.false.

c-----------------------------------------------------------------------
c     namelist declarations.  note that some variable appear in more
c     than one namelist for compatibility with the central repository.
c-----------------------------------------------------------------------
c=======================================================================
      NAMELIST/grid_input/gridshape,eqfile,xmin,xmax,ymin,ymax,mx,my,
     $	   mxpie,nxbl,nybl,firstx,firsty,skew,geom,per_length,
     $	   decompflag,periodicity,xo,yo,lphi,lin_nmodes,zperiod,
     $     pieflag,rimflag,nbl_rim,lin_nmax,nlayers,rim_length,
     $     qpack,wpack,amp,vac_frac,x1,y1,poly_degree,dealiase,packbigr,
     $     quadxx,quadyy,quadarc
c=======================================================================
      NAMELIST/const_input/chrg_input,zeff_input,mi_input,me_input,
     $     gam_input,kblz_input,mu0_input,c_input
c=======================================================================
      NAMELIST/init_input/nx,ny,nz,init_type,bamp,reset_file,reset_tol,
     $     rescale_factor,reset_time,reset_step,transfer_eq,init_x,
     $     init_y,blob_amp,blob_width,ferr_amp,ferr_phase,ferr_n
c=======================================================================
      NAMELIST/equil_input/zeff,ndens,nedge,n_profile,npp,be0,
     $     pe_profile,te_cent,te_edge,tepp,thetab,phib,lamprof,lam0,
     $     rbreak,alpha,beta,dvac,dexp,ds_use,ds_function,xvac,eq_flow,
     $     ve0,thetav,phiv,eqflow_width,pit_0,pit_2,pit_4,pres_2,pres_4,
     $     divbeqd,glength,tor_eqja_fe,tanh_byl,tanh_prl,tanh_pfrac,
     $     tanh_ndl,tanh_nfrac,ncoil,coil_current,coil_r,coil_z
c=======================================================================
      NAMELIST/physics_input/ndens,be0,thetab,phib,bamp,nonlinear,
     $	   nx,ny,init_type,ohms,elecd,lamprof,lam0,rbreak,alpha,
     $     beta,kin_visc,dvac,dexp,advect,loop_volt,tloopv0,
     $     tloopv1,eq_flow,ve0,thetav,phiv,eqflow_width,
     $     pit_0,pit_2,pit_4,pres_2,pres_4,ds_use,i_desired,loop_rate,
     $     separate_pe,pe_frac,loop_rate2,e_vertical,t_e_vert0,
     $     t_e_vert1,xvac,vsink_rate,advance_pvi,continuity,nd_diff,
     $     eta_model,eta_ref_t,elecd_max,elecd_min,zero_bnorm,iso_visc,
     $     par_visc,gyr_visc,norm_flow,flow_bc,nd_bc,hv_power,hv_filter,
     $     hv_coeff,hv_time,gravity,rt_transfer_eq,nd_hypd,f_chodura,
     $     elecd_chodura,tanh_byl,tanh_prl,tanh_pfrac,pres_offset,
     $     tanh_ndl,tanh_nfrac,gamma_nimset
c=======================================================================
      NAMELIST/closure_input/neoe_flag,neoi_flag,neo_debug,
     $     mu_e,mu_i,ng_ft,m_neo,n_neo,neo_rad,p_model,k_perp,k_pll,
     $     neoe_det,neoi_det,k_perpe,k_plle,k_perpi,k_plli,tequil_rate,
     $     ohm_heat,visc_heat,k_pll_max,k_pll_ref_t,k_pll_min,
     $     kprp_mnrat,k_cross,insulate,gamma_heat,heat_source,
     $     electric_source,heat_ion_source,heat_ele_source,
     $     momentum_source,cel_model,parvisc_model,tdep_tequil,
     $     tdep_coul_log,coulomb_logarithm,magfac_ele,magfac_ion,
     $     delta0_ele,delta1_ele,gamma0_ele,gamma1_ele,
     $     delta0_ion,delta1_ion,gamma0_ion,gamma1_ion,closure_model
c=======================================================================
      NAMELIST/numerical_input/dtm,tmax,nstep,npc,fom,si_fac_hall,
     $     conform,lump_b,ngr,feta,divbd,fdivb,ndivb,met_spl,
     $     si_fac_pres,si_fac_mhd,fvsc,mhd_si_iso,dt_stop,v_cfl,
     $     split_resist,fb_vxb,fv_vdgv,fp_vdgp,fn_vdgn,kdivb_2_limit,
     $     dt_change_frac,ave_change_limit,dt_incr,split_divb,
     $     nl_cfl_lim,si_in_v,split_visc,cpu_tmax,dt_initial,lump_all,
     $     si_fac_j0,r0dr_weight,si_fac_nl,transfer_eq,poly_degree,fthc,
     $     n_dt_release,mhdadv_alg,nd_floor,nd_exp,nd_dart_fac,
     $     nd_dart_upw,t_dart_upw,upw_aniso,upw_limit,nd_floor_upw,
     $     nd_width_upw,t_floor_upw,t_width_upw,integration_formula,
     $     nstep_converge,growth_eps,n_mat_update,p_computation,
     $     poly_divb,disc_dbd,poly_divv,poly_divv_min,fdivv,fpvrt,
     $     poly_divv_max,ddivv,dpvrt,poly_divb_min,poly_divb_max,
     $     poly_divv_auto,hyp_eta,fhyp_eta,split_hypeta,hyp_dbd,
     $     fhyp_dbd,siop_type,fmx_drate,nfdamp,nd_nodal_floor,
     $     ti_nodal_floor,te_nodal_floor,nd_correrr,hyp_visc,fhyp_visc
c=======================================================================
      NAMELIST/solver_input/tol,maxit,solver,off_diag_fac,extrap_order,
     $     vmhd_solver,bmhd_solver,temp_solver,n_gm,nsym_pre_band,
     $     nsym_pre_rpass,nsym_pre_rfac,nsym_pre_rtype,maxit_nl,tol_nl
c=======================================================================
      NAMELIST/output_input/detflag,nhist,ihist,jhist,hist_binary,ndump,
     $     dump_name,dump_dir,ndxout,dump_over,dump_file,xy_stride,
     $     xt_stride,yt_stride,x0fac,y0fac,xdraw_dir,dx_dir,itflag,
     $     reset_file,reset_time,reset_step,dump_list,elout_format,
     $     reset_eq_mesh,hist_flush
c=======================================================================
      NAMELIST/particle_input/nm,restart,vhmx,avu0,avw0,
     $     mass0,charge0,psipmax,Ppsi0,R0,betafrac,tpartdr,trace

c-----------------------------------------------------------------------
c     open input file.
c     Remove comments from input file and put into temporary file.
c-----------------------------------------------------------------------
      tempfile='tempnimrod.in'
      CALL rmcoment(infile,tempfile)
      OPEN(UNIT=in_unit,FILE=tempfile,STATUS='OLD',POSITION='REWIND')
c-----------------------------------------------------------------------
c     check namelist file for namelist order and number.
c-----------------------------------------------------------------------
      DO
        READ(UNIT=in_unit,FMT='(a)',IOSTAT=read_stat) ctmp 
        IF (read_stat/=0) EXIT
        nc=LEN_TRIM(ctmp)
        IF (nc<1) CYCLE
        ctmp=ADJUSTL(ctmp)
        reading=.false.
        IF (ctmp(1:1)=='&') THEN
          number_of_namelists=number_of_namelists+1
c-----------------------------------------------------------------------
c         trim all but the namelist name.
c-----------------------------------------------------------------------
          DO il=2,nc+1
            IF (ctmp(il:il)/=' ') THEN
              IF (.NOT.reading) inst=il
              reading=.true.
              CYCLE
            ENDIF
            IF (ctmp(il:il)==' '.AND.reading) THEN
              ctmp=ctmp(inst:il-1)
              EXIT
            ENDIF
          ENDDO
          BACKSPACE(in_unit)
c-----------------------------------------------------------------------
c         select and read namelist.
c-----------------------------------------------------------------------
          SELECT CASE(TRIM(ctmp))
          CASE('grid_input')
            READ(UNIT=in_unit,NML=grid_input,IOSTAT=read_stat)
          CASE('const_input')
            set_phys_constants=.true.
            READ(UNIT=in_unit,NML=const_input,IOSTAT=read_stat)
          CASE('init_input')
            READ(UNIT=in_unit,NML=init_input,IOSTAT=read_stat)
          CASE('equil_input')
            READ(UNIT=in_unit,NML=equil_input,IOSTAT=read_stat)
          CASE('physics_input')
            READ(UNIT=in_unit,NML=physics_input,IOSTAT=read_stat)
          CASE('closure_input')
            READ(UNIT=in_unit,NML=closure_input,IOSTAT=read_stat)
          CASE('numerical_input')
            READ(UNIT=in_unit,NML=numerical_input,IOSTAT=read_stat)
          CASE('solver_input')
            READ(UNIT=in_unit,NML=solver_input,IOSTAT=read_stat)
          CASE('output_input')
            READ(UNIT=in_unit,NML=output_input,IOSTAT=read_stat)
          CASE('particle_input')
            READ(UNIT=in_unit,NML=particle_input,IOSTAT=read_stat)
          CASE DEFAULT
            CALL nim_stop
     $        (TRIM(ctmp)//' is an unrecognized namelist.')
          END SELECT
          IF (read_stat/=0) CALL nim_stop
     $      ('Error reading namelist '//TRIM(ctmp)//'.')
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     close input file.
c       Delete it since it is the temporary file
c-----------------------------------------------------------------------
      CLOSE(in_unit,STATUS='DELETE')
c-----------------------------------------------------------------------
c     echo the input parameters to the output file.
c-----------------------------------------------------------------------
      IF (echo_in) THEN
        WRITE(out_unit,'(a,/)') 'grid_input:'
        WRITE(UNIT=out_unit,NML=grid_input)
        IF (set_phys_constants) THEN
          WRITE(out_unit,'(/,a,/)') 'const_input:'
          WRITE(UNIT=out_unit,NML=const_input)
        ENDIF
        WRITE(out_unit,'(/,a,/)') 'init_input:'
        WRITE(UNIT=out_unit,NML=init_input)
        WRITE(out_unit,'(/,a,/)') 'equil_input:'
        WRITE(UNIT=out_unit,NML=equil_input)
        WRITE(out_unit,'(/,a,/)') 'physics_input:'
        WRITE(UNIT=out_unit,NML=physics_input)
        WRITE(out_unit,'(/,a,/)') 'closure_input:'
        WRITE(UNIT=out_unit,NML=closure_input)
        WRITE(out_unit,'(/,a,/)') 'numerical_input:'
        WRITE(UNIT=out_unit,NML=numerical_input)
        WRITE(out_unit,'(/,a,/)') 'solver_input:'
        WRITE(UNIT=out_unit,NML=solver_input)
        WRITE(out_unit,'(/,a,/)') 'output_input:'
        WRITE(UNIT=out_unit,NML=output_input)
        WRITE(out_unit,'(//)')
      ENDIF
c-----------------------------------------------------------------------
c     check for old input specifications, inconsistencies, and set
c     any unset defaults.
c-----------------------------------------------------------------------
      IF (.NOT.conform)
     $  CALL nim_stop('The conform=F option is not available.')
      IF (.NOT.si_in_v)
     $  CALL nim_stop('The si_in_v=F option is not available.')
      IF (lump_b)
     $  CALL nim_stop('The lump_b=T option is not available.')
      IF (ngr<2)
     $  CALL nim_stop('Setting ngr<2 may produce singular matrices.')
      IF (ndxout>0)
     $  CALL nim_stop('DX output is no longer produced directly '//
     $  'from nimrod.  Set ndxout to 0.')
      IF (continuity/='none'.AND.continuity/='fix profile'.AND.
     $    continuity/='n=0 only'.AND.continuity/='full')
     $  CALL nim_stop('Continuity value '//TRIM(continuity)//
     $  ' is not recognized.')
      IF (continuity=='none'.AND.beta>0) CALL nim_stop('Continuity '
     $  //'must not be set to none if beta>0.')
      IF (split_visc.AND.kin_visc<=0) CALL nim_stop('Set split_visc='
     $  //'F when kin_visc=0.')
      IF (split_visc.AND.(iso_visc>0.OR.par_visc>0)) CALL nim_stop
     $  ('Isotropic and parallel viscosities do not work with '
     $  //'split_visc.')
      IF (eta_model/='fixed'.AND.eta_model/='eta n=0 only'.AND.
     $    eta_model/='eta full'.AND.eta_model/='chodura')
     $  CALL nim_stop('Eta_model value '//TRIM(eta_model)//
     $  ' is not recognized.')
      IF (split_resist) CALL nim_stop('The split_resist option is no '
     $  //'longer available.')
      IF (met_spl=='cubic'.OR.met_spl=='bicubic')
     $  CALL nim_stop('Met_spl value '//TRIM(met_spl)//' is not '
     $  //'available with the improved mapping.')
      IF (gyr_visc>0.AND.beta<=0)
     $  CALL nim_stop('Gyroviscosity requires beta>0.')
      IF (par_visc>0.AND.parvisc_model/='fixed'.AND.beta<=0)
     $  CALL nim_stop('T-dependent parallel viscosity requires beta>0.')
c-PRE
      IF (k_cross>0.AND.nonlinear)
     $  CALL nim_stop('Thermal drifts are not yet '
     $  //'ready for nonlinear calculations.')
      IF (k_cross>0.AND.separate_pe.AND..NOT.p_model(1:5)=='aniso')
     $  CALL nim_stop('Thermal drifts require the anisotropic p_model.')
      IF (vsink_rate>0)
     $  CALL nim_stop('Vsink_rate is no longer available.')
      IF (p_model=='aniso_plltdep'.AND.closure_model/='std kprp n0')
     $  CALL nim_stop('Closure_model '//TRIM(closure_model)//
     $                ' requires aniso_tdep p_model.')
      IF (solver(1:6)=='bl_ilu')
     $  CALL nim_stop('The block-ILU preconditioners are no longer'//
     $                 ' available.')

      IF (p_model=='iso') p_model='isotropic'
      IF (k_perpi==0.AND.k_perp>0) k_perpi=k_perp
      IF (k_perpe==0.AND.k_perp>0) k_perpe=k_perp
      IF (k_plli==0.AND.k_pll>0) k_plli=k_pll
      IF (k_plle==0.AND.k_pll>0) k_plle=k_pll
      IF (vmhd_solver=='none') vmhd_solver=solver
      IF (bmhd_solver=='none') bmhd_solver=solver
      IF (temp_solver=='none') temp_solver=solver
      IF (.NOT.nonlinear.AND.(p_model=='aniso_plltdep'.OR.
     $                        p_model=='aniso_tdep')) CALL nim_stop
     $   ('The selected p_model is intended for nonlinear runs only.')
      IF (.NOT.nonlinear.AND.parvisc_model=='plltdep'.AND.
     $                       par_visc>0) CALL nim_stop
     $   ('The selected parvisc_model is for nonlinear runs only.')
      IF (neoi_flag/='none'.OR.neoe_flag/='none')
     $  CALL nim_stop('The neoclassical options are not available'
     $  //' in this nimuw version.')
      IF (solver=='bl_drect') solver='lapack'

      IF (poly_divv_auto) THEN
        poly_divv=poly_degree
        poly_divv_min=poly_degree
        poly_divv_max=poly_degree
      ENDIF
      IF (poly_divv_max<0) poly_divv_max=poly_divv
      IF (poly_divv_max>=0) poly_divv_min=MAX(poly_divv_min,0_i4)
      IF (poly_divb_max<0) poly_divb_max=poly_divb
      IF (poly_divb_max>=0) poly_divb_min=MAX(poly_divb_min,0_i4)
      IF (ohm_heat.AND.nonlinear.AND.hyp_dbd>0) fhyp_dbd=fhyp_eta
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END SUBROUTINE read_namelist
c-----------------------------------------------------------------------
c     subprogram 2. rmcomment
c     routine strips out comments beginning with an exclamation point
c-----------------------------------------------------------------------
      SUBROUTINE rmcoment(fileold,filenew)

      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: fileold,filenew
      CHARACTER(128) :: line
      INTEGER, PARAMETER :: nold=55,nnew=56
      INTEGER cmax, ios
      LOGICAL :: file_exist
c-----------------------------------------------------------------------
c     Open files, but make sure the old one exists first.
c-----------------------------------------------------------------------
      INQUIRE(FILE=fileold,EXIST=file_exist)
      IF(.NOT. file_exist) THEN
         PRINT *,'The file "',fileold,'" could not be found.'
         STOP
      ENDIF

      OPEN(UNIT=nold,FILE=fileold,status="OLD")
      OPEN(UNIT=nnew,FILE=filenew,status='REPLACE')

c-----------------------------------------------------------------------
c     Strip comments.     Note: line lengths limited to 127 characters
c-----------------------------------------------------------------------
      DO
        READ(UNIT=nold,FMT='(a)',IOSTAT=ios) line
        IF (ios /= 0) EXIT
        cmax=1
        DO WHILE(line(cmax:cmax).NE.'!' .AND. cmax .LE. 127)
           cmax=cmax+1
        ENDDO
        IF(cmax .GT. 1) WRITE(nnew,'(a)') line(1:cmax-1)
      ENDDO

c-----------------------------------------------------------------------
c     Close files and exit
c-----------------------------------------------------------------------
      CLOSE(nold)
      CLOSE(nnew)

      RETURN
      END SUBROUTINE rmcoment
