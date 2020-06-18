c-----------------------------------------------------------------------
c     subprogram input
c     handles input.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     subprogram 0. input
c     module declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     input variable declarations.
c-----------------------------------------------------------------------
      MODULE input
      USE local
      IMPLICIT NONE


c-----------------------------------------------------------------------
c    filename gives the path of the data file specifying the equilibrium.  
c     The file may be ascii or binary.  It may be the output of a direct
c     Grad-Shafranov solver, which specifies poloidal flux, R*B_phi, pressure,
c     and safety factor on a 1D grid of flux surfaces, and the poloidal flux
c     psi on a 2D grid in R and Z; or of an inverse solve, which replaces
c     psi(R,Z) with values of R and Z on a 2D grid of psi and theta on flux
c     surfaces.

      CHARACTER(128) :: filename="82205-1.25X"

c-----------------------------------------------------------------------
c     eq_type specifies the type of the input file.  It allows fluxgrid to 
c     determine whether the file is ascii or binary, whether it is direct 
c     or inverse, and how to interpret the data.  
c     Here is a list of allowed types:
c   
c     Name          Description               Direct/Inverse  Ascii/Binary
c     --------------------------------------------------------------------
c     "miller"      Bob Miller's TOQ code       inverse        binary (r4 data)
c     "miller8"     Bob Miller's TOQ code       inverse        binary (r8 data)
c     "millasc"     Bob Miller's TOQ code       inverse        binary (ASCII)
c     "galkin"      Sergei Galkin's code        inverse        binary
c     "chease"      Lausanne CHEASE code        inverse        ascii
c     "chum"        Modified Miller code        inverse        binary
c     "efit"        GA EFIT code                direct         ascii
c     "rsteq"       ORNL RSTEQ code             direct         binary
c     "soloviev"    Soloviev analytical model   direct         binary


      CHARACTER(8) :: eq_type="efit"

c-----------------------------------------------------------------------
c     The indices of the grid are 0:mpsi and 0:mtheta and correspond to
c     mx and my in the nimset/nimrod input file except when extending
c	the grid into the vacuum region.  See extr and mvac below.
c     mpsi is the number
c     of radial cells between the magnetic axis and the last azimuthal
c     grid line (before the separatrix).  mtheta is the number of
c     azimuthal cells.
c
c     Don't touch mxpie unless you really know what you are doing.

      INTEGER(i4) :: mpsi=32   
      INTEGER(i4) :: mtheta=32
      INTEGER(i4) :: mxpie=1    ! radial grid position of pie-shaped tblock
      INTEGER(i4) :: mx=128     ! adapted nodes

c-----------------------------------------------------------------------
c     The grid in the grad_psi direction is controlled by grid_method
c     and pack_method.  
c
c     grid_method controls what to use as the "radial" coordinate in 
c     a "global" type of sense.
c      grid_method="uniform", distributes radial grid points uniformly 
c       in rho=psi_normal.
c      grid_method="sqrt", distributes radial grid points uniformly 
c       in rho=(psi_normal)^1/2  This puts more points near the center
c       as compared to "uniform" method.
c      grid_method="stretch" (suggested by Don Pearlstein) distributes points
c       at both the axis and the edge as compared to the "uniform" method.
c     For inverse codes (such as "read_miller") types, grid_method="sqrt"
c      is recommended.
c
c     pack_method is used to control how to pack around rational surfaces which
c     is most useful for interchange and tearing modes.
c     pack_method="lorentz" and "no_gap" use the monitor function approach 
c     by Alan Glasser, and the "gauss" choice uses a Gaussian distribution
c     For more information on these, see below.
c	Note that for pack_method="lorentz" or "no_gap", grid_method is
c	automatically set to "sqrt"
c      
      CHARACTER(8) :: grid_method="use_gt"   
                                        ! "original", "uniform",
                               		! "stretch", "sqrt",
					! "use_gt" use old style grid_type
      CHARACTER(8) :: pack_method="none"   
                                        ! "none","lorentz","no_gap","gauss"

c-----------------------------------------------------------------------
c     The grid in the theta direction is chosen by angle_method.
c
c     angle_method="jac"  specifies the theta angle by  choosing the
c     Jacobian: Jacobian = R^ipr/|B_P|^ipb .  For a discussion of the
c     types of grid chosen by this method, see the paper:
c	Grimm, Dewar, Manickam JCP 49 (1983) 94
c
c     angle_method="geom" creates the theta contours by forming straight
c     rays coming out of the magnetic axis.  Unlike the delta-W codes,
c     this is probably a good choice for NIMROD.
c     The default is "jac" to make it backword compatible with previous
c     versions
c
c     angle_pack chooses how to distribute the theta mesh in the poloidal
c     direction:  
c       "none" => theta is equally spaced (dtheta=constant) 
c                 This is the default, and is no packing
c       "arclength" => the theta is distributed such that on the last
c               closed flux surface, the arc-lengths are equally spaced
     
      CHARACTER(8) :: angle_method="jac"   
      CHARACTER(10) :: angle_pack="none"   
      INTEGER(i4) :: ipb=1
      INTEGER(i4) :: ipr=0


c-----------------------------------------------------------------------
c     For the RBLOCK region:
c     The input parameter psilow is used to cut off stuff near the axis
c     which can be useful for when the equilibrium codes have problems 
c     near the axis.  A good value is 1e-4.  
c     
c     The input parameter psihigh determines how close to the edge the ODE's
c     are integrated.  It must be <= 1.  For an inverse equilibrium, which
c     never has a separatrix in the computational domain, it may be set to 1.
c     For a direct equilibrium with a separatrix, psihigh determines how close
c     the code approaches the separatrix. For some EFIT files with large
c     currents near the edge, this can be useful but difficult to use
c     correctly since nimrod is sensitive to currents at the edge.  Get some
c     advice on those cases
c
c     Currently, the capability for triangular regions is only contained when
c     reading EFIT (eqdsk) files.  In this case, we have the following choices:
c           psihigh > 0, triangle=.F.      => pure rblock
c           psihigh > 0, triangle=.T.      => rblock and tblock regions
c           psihigh = 0, triangle=.T.      => pure tblock 
c           psihigh = 0, triangle=.F.      => exits.
c
c     The output files are fluxgrid.dat for the rblock regions and nimdsk (a
c     generalized eqdsk file) and poly files.  The poly files are normally
c     meant to be used with the TRIANGLE code to produce the grid.  Then
c     the grid and the equilibrium contained in the nimdsk file are combined
c     in nimset to produce the dump file for nimrod.
c
c
c     iextr and extr allows one to extend the rblock region beyond 
c	that specified by psihigh.  
c
c     extr =< 1. gives the default behavior while extr > 1 
c	extends the rblock region as described below.  Note that only
c	fluxgrid.dat is affected by this option.  All of the other
c	output (graphs and interfaces to other codes are unaffected
c	by this option).
c
c     The specific method distance between the outer plasma boundary and
c       the rblock boundary is specified by extr_type:
c	1) extr_type='conformal'
c		The distance between the wall point and the
c		outer boundary is constant.  The boundary is defined
c		by psihigh.
c	2) extr_type='self_similar'
c		The distance extr=b/a is held constant where b is the
c		distance from the wall point to the origin and a is 
c		the distance from the edge to the origin.  The wall
c		point, the boundary point, and the origin are all 
c		colinear.  The center is chosen using extr_center,
c		and can either be 'geom' for the geometric center (the
c		conventional GATO choice) or 'mag_axis' for the
c		magnetic axis (the conventional MARS choice).
c	For NIMROD's purposes, the 'self-similar' grid is usually superior.
c
c     mvac is the number of points between the outer boundary (as defined
c	by psihigh) and the boundary described by extr.  This normally
c	corresponds to the vacuum region (psihigh ~ 1), but isn't necessarily
c	true.
c     NOTE: To match with nimrod.in: Let mx=mpsi+mvac
c
c     If vac_pack=.TRUE. then packing is done around the boundary of the
c       outer boundary (defined by psihigh) and the vacuum region.
c       The values for the packing are set by amp and wpack in the
c       same manner as Guassian packing.  See that section for details.
c       CURRENTLY NOT RECOMMENDED.
c
c     the firstv option replaces other vacuum-packing schemes.  the
c	first element from the separatrix is the fraction firstv of the
c	distance to the outer boundary, and dr then varies linearly
c	to meet the outer boundary.  this is similar to firstx in nimset.
c
c     wall_length is a minimum wall length used in constucting the
c       poly files associated with the outer boundary.  When negative,
c       an internal scheme finds a minimun length by inspecting the
c       rblock boundary and using twice this minimum length,
c       otherwise wall_length is this minimum value.

      REAL(r8) :: psilow=.01    	! minimum psi value, > 0
      REAL(r8) :: psihigh=.95   	! maximum psi value, < 1
      REAL(r8) :: extr=0.		! extend rblock region
      REAL(r8) :: firstv=0.		! default does not pack
      INTEGER(i4) :: mvac=0     	!  # of radial points in "vacuum"
      CHARACTER(12) :: extr_type='self_similar'	! Wall boundary type
      CHARACTER(12) :: extr_center='geom'	! Used with self-similar
      LOGICAL :: vac_pack=.FALSE.
      LOGICAL :: triangle=.TRUE.	! Not yet functional in fluxgrid

      REAL(r8) :: wall_length=-1.0


c-----------------------------------------------------------------------
c     For newq0 = 0, the code uses the original q profile specified in the
c     equilibrium file.  For newq0 /= 0, it readjusts the q profile to give
c     the specified value of q at the axis, in such a way as to leave the
c     Grad-Shafranov solution invariant.  This can be used to explore the
c     stability of a range of equilibria for each equilibrium file.

      REAL(r8) :: newq0=0    ! desired q on axis; 0 => keep original

c-----------------------------------------------------------------------
c     tol0 gives the tolerance for mapping direct equilibria
c     The smaller the number the longer fluxgrid will take.  
c     10^-6 is OK, if there isn't any radial or poloidal packing
c     10^-8 is a safer choice however.
c
c     interp determines whether cubic spline interpolation is done for
c     some of the files, more detail at the expense of slower operation
c     

      REAL(r8) :: tol0=1.e-8    ! tolerance for mapping direct equilibria
      LOGICAL :: interp=.FALSE.

c-----------------------------------------------------------------------
c     For direct equilibria codes, it is possible to evaluate the 
c     magnetic fields directly from the input without going through
c     the entire mapping.  This is useful for debugging cases where
c     one is not sure if the crappy equilibrium your feeding into NIMROD
c     is fluxgrid's fault or the equilibrium code's fault.
c     
      LOGICAL :: calc_direct=.FALSE.

c-----------------------------------------------------------------------
c     sfac and ndens are not used by fluxgrid in any significant way
c     since equilibria do not depend on S or n; however, one often wants
c     to know the value of S given a density, or vice versa, since
c     these numbers are important for running nimrod. fluxgrid calculates them
c     for you.

      REAL(r8) :: sfac=1e4      ! Lundquist number
      REAL(r8) :: ndens=1e20    ! plasma particle density


c-----------------------------------------------------------------------
c     When packing around the rational surfaces, one has to choose which
c      rational surfaces to pack around.  Here one chooses the toroidal
c      mode number to pack around, and then gives a range of q values
c      to pack around.  For example, two pack around the q=2 and q=3/2
c      surfaces, one can choose nn=2 and qsmin = 1.49 and qsmax=2.01
c      (the q=2 surface corresponds to m=4/n=2).
      INTEGER(i4) :: nn=1       
      REAL(r8) :: qsmin=-HUGE(0) 
      REAL(r8) :: qsmax=HUGE(0) 

c-----------------------------------------------------------------------
c     For pack_method="lorentz" or "no_gap" one can specify the width of
c	the packing via delta0.  

      REAL(r8) :: delta0=.01 

c-----------------------------------------------------------------------
c     For pack_method="gauss", the packing is Gaussian (i.e.,
c       a plot of dr vs. i shows a Gaussian) with the amplitude and width
c       of the Guassian by amp and wpack.  The amplitude is roughly
c       dr0/dr where dr0 is the value far away from the Gaussian and wpack
c       is the full-width given as a percentage of a = xmax-xmin.
c       For example, if mpsi=164 and you want to pack an extra 100 points
c       around the q=1 rational surface within a width of 0.1*a then:
c               nn=1, qsmin = 0.99, qsmax=1.01
c               wpack=0.1
c               amp=15.625 
c        where amp = dr0/dr = 1/(164-100) / (0.1*1/100) = 15.625

      REAL(r8) :: wpack=0._r8
      REAL(r8) :: amp=0._r8

c-----------------------------------------------------------------------
c     For the toroidal current, there are two ways of calculating it:
c	"delstar"  -- direct evaluation of J_t = delstar(psi)
c	"use_gs"   -- use the Grad-Shafranov equation: J_t = -R^2 mu0 p'-FF' 
c	use_gs is recommended because there is less numerical error 
c	  associated with the evaluation of the derivatives.  However,
c	  if you don't think you are reading in a good equilibrium, go for
c	  the more direct evaluation. 

      CHARACTER(8) :: j_t="use_gs"      

c-----------------------------------------------------------------------
c     To calculate the magnetic field in the R,phi,Z coordinates of 
c     nimrod, one requires the Jacobian for the transformation of the
c     flux coordinates to R,Z coordinates.  Fluxgrid has 2 
c     coordinate systems: (rho,eta) and (R,Z).  (rho, eta) tends to work
c     better near the magnetic axis, while sometimes R,Z works better 
c     for the edge in highly shaped plasmas.  The 'rjac' parameter gives
c     the "radius" (in terms of r=SQRT(psi_normal) at which to switch
c     between using rho,eta to calculate the grid, and R,Z coordinate
c     system to calculate the grid.  Basically, if you are having problems
c     with the numerical noise near the edge, switch this puppy to one and
c     see if that helps.

      REAL(r8) :: rjac=0.5

c-----------------------------------------------------------------------
c     A bunch of parameters here for controlling the output files.
c     See the README file here for more info.

      LOGICAL :: stability=.FALSE.
      LOGICAL :: bin_1d=.FALSE.
      LOGICAL :: out_1d=.FALSE.
      LOGICAL :: bin_2d=.FALSE.
      LOGICAL :: out_2d=.FALSE.
      LOGICAL :: out_pies=.FALSE.
      LOGICAL :: out_pest=.FALSE.
      LOGICAL :: out_far=.FALSE.
      LOGICAL :: out_neoclassical=.FALSE.
      LOGICAL :: out_tecplot=.FALSE.

c-----------------------------------------------------------------------
c     These are in here for historical reasons and have been superseded.

      INTEGER(i4) :: grid_type=3        ! 0:original; 1:uniform,
                               		! 2:stretch, 3:sqrt, 4:pack
      CHARACTER(8) :: mon_type="none"   ! "lorentz" or "no_gap"
                                        ! monitor functn type for grid_method ="pack" 

c-----------------------------------------------------------------------
c     namelist declarations. 
c-----------------------------------------------------------------------
      NAMELIST/global_input/filename,mpsi,mtheta,mxpie,psilow,psihigh,
     $     newq0,tol0,j_t,rjac,interp,calc_direct,
     $     eq_type,nn,qsmin,qsmax,delta0,sfac,ndens,
     $     bin_1d,out_1d,bin_2d,out_2d,stability,
     $     out_tecplot,out_neoclassical,out_pest,out_far,out_pies,
     $     grid_method,pack_method,wpack,amp,mx,
     $     angle_method,angle_pack,ipb,ipr,
     $     extr,mvac, extr_type, extr_center,vac_pack,firstv,
     $     triangle,wall_length,
     $     grid_type,mon_type


      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. read_input.
c     reads global input.
c-----------------------------------------------------------------------
      SUBROUTINE read_input

      CHARACTER(14) :: tempfile,infile
c-----------------------------------------------------------------------
c     Strip the comments from the input file and store into a 
c 	temporary file which is more like the normal file
c-----------------------------------------------------------------------
      infile='fluxgrid.in'
      tempfile='tempfg.in'
      CALL rmcoment(infile,tempfile)
      OPEN(UNIT=eq1_unit,FILE=tempfile,STATUS="UNKNOWN")
      READ(UNIT=eq1_unit,NML=global_input)
      CLOSE(UNIT=eq1_unit,STATUS="DELETE")
c-----------------------------------------------------------------------
c     Backward compatability conversions:
c-----------------------------------------------------------------------
      IF(TRIM(grid_method) == 'use_gt')THEN
         SELECT CASE(grid_type)
         CASE(0)
           grid_method = 'original'
         CASE(1)
           grid_method = 'uniform'
         CASE(2)
           grid_method = 'stretch'
         CASE(3)
           grid_method = 'sqrt'
         CASE(4)
           grid_method = 'sqrt'
         CASE default
           grid_method = 'sqrt'
         END SELECT
         IF (mon_type/='none')pack_method=mon_type
      ENDIF
c-----------------------------------------------------------------------
c     extr determines whether we extend to vacuum or not.  We also need
c     mvac > 0
c-----------------------------------------------------------------------
      IF(extr < 1.) mvac=0
      IF(mvac <= 0) extr=0.
c-----------------------------------------------------------------------
c     the input mpsi has been changed in cvs Version 3_0 to be the
c     number of radial cells between the magnetic axis and the last
c     azimuthal grid line.  fluxgrid still uses mpsi as the number of
c     cells between psilow and psihigh, however.
c
c     if "lorentz" or "no_gap" pack_methods are used, mpsi is used to
c     create a temporary sqrt-type grid from which a packed grid with
c     dimension mx will be created.
c-----------------------------------------------------------------------
      mpsi=mpsi-1
      mx=mx-1
c-----------------------------------------------------------------------
      IF(TRIM(pack_method)=='lorentz'.OR.TRIM(pack_method)=='no_gap')
     &                    grid_method='sqrt'
c-----------------------------------------------------------------------
      END SUBROUTINE read_input

c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      END MODULE input



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

      OPEN(UNIT=nold,FILE=fileold,status='OLD')      
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
