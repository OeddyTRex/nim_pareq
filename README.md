# nim_me759

ME759 final project - parallelizing nimeq's field line length calculations

Appended below are the contents of the README from the trunk of nimuw
version 3.5.2 from which this directory has been cloned.

This project will concern only the nimeq portion of the nimrod suite
of codes. 

--------------------------------------------------------------------------

NIMROD CVS revision (nimuw) 3_5_2, 3/23/20

  This README file is intended to help a user get started with the
nimrod package.  It contains a description of the contents of each
subdirectory in the nimrod package, and it contains instructions on
how to install the codes in the package.


  Directory organization is first:

1.  nimrod:  contains the physics kernel code.
  Only sources are stored, so the executable needs to
be made before running (see below).  The physics kernel executable
will be located here.

2.  nimset:  contains the pre-processor code.
  This code reads a flux-grid or generates a simple grid, and it
initializes the fundamental physical quantities.  This information is
then stored in an initial restart dump, which is read by nimrod.  The
nimset executable will be located here after the build is complete.

3.  fluxgrid:  contains a grid generation code.
  This code (grid) reads output from various equilibria codes
(GA's EFIT, for example) to generate a magnetic flux grid.  The output
is used by nimset.

4.  input:  contains sample input files for nimrod.
  These files may be read by the nimset and nimrod.  Sample input
for other codes in the suite are also in this directory.

5.  draw:  contains useful input files for the xdraw postprocessor.
Some of the xdraw files are available directly under the GUI's.

6.  nimplot:  contains programs for postprocessing.  The first program
is nimplot.  It is interactive and is used to create graphics files
from the nimrod dump files.
  The second program is xlog, which can be used to take the logarithm of
data in xdraw binary output files (such as nimhist.bin) and write this
into a new xdraw binary.  It is also interactive.

7.  nimfl:  contains a field line tracing routine that uses nimfl.in
as input (see the input directory for an example).  It works directly
from the nimrod dump files.

8. nimlib: contains commonly used F90 modules such as those defining
basis functions.

9. dump_convert: has code to help convert from older format dump
files to the format for the present version of nimrod.

10. externals: contains libraries to facilitate calling external
mathematical libraries that are written in C and dummy libraries to
satisfy the loader when a library is not needed and not available.

11. make_includes: include files with machine-dependent information
for the make operation.

12. nimbase: this utility program allows a user to change the degree of
the polynomials and/or node spacing of the finite-element basis
functions.  A new dump file is produced.

13. nimeq: this is a Grad-Shafranov solver that uses nimrod's mesh and
basis functions.

14. nimcore: contains modules and libraries that are used by more
than one of the codes but with higher-level module dependencies than
those in nimlib.

15. scripts: contains at least one script that facilitates workflow.


  Installation procedure:

 First, check for external libraries.  Nimrod, and other codes
in the suite can be built with just the basic BLAS library,
but external libraries are needed for full functionality.  Super-
computing centers may have some or all of these libraries
installed for users.  These libraries are:

1. Message-passing interface (MPI) -- MPI is needed for parallel
computation.

2. Odepack -- Fluxgrid, nimfl, and nimeq use odepack for tracing
field-lines and solving other ODEs.  It is available from
https://computation.llnl.gov/casc/odepack/

3. SuperLU, SuperLU_DIST -- Nimrod and nimeq typically use SuperLU
and its parallel, distributed-memory version for preconditioning
Krylov-space solves of algebraic systems.  They are available from
http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ .

4. METIS and ParMETIS -- SuperLU can use these libraries to determine
efficient orderings.  See
http://glaros.dtc.umn.edu/gkhome/metis/metis/overview .

5. LAPACK -- Although it is not recommended, this algebraic solver
for dense systems can be used by nimrod as a preconditioner.  Parts of
SuperLU and SuperLU_DIST may use LAPACK.


 Next, copy the appropriate make_<machine>.inc file from the
make_includes directory to "make.inc" at the top level of the
directory structure.  You may need to modify paths to link the desired
math libraries.  A make.inc template is also available in
make_includes if you wish to start from scratch for a new machine.
Now type "make" at the top level of the directory structure.  This
will run make in each of the subdirectories where executables and
nimrod-specific libraries will reside.  It also creates the rundir
subdirectory with convenient softlinks.  To build individual executables
in the suite, you may type "make nimrodex nimsetex" for example.




  A note for Linux users: when using the Lahey/Fujitsu fortran, you
will need to set the environment variable FORT90L set to "-Wl,-T" to
have binary input and output in standard portable format (big endian).

csh:
setenv FORT90L -Wl,-T

ksh:
FORT90L="-Wl,-T"; export FORT90L

However, we have had difficulty getting environment variables
transferred to nodes of parallel Linux clusters.  Here it may be best
to specify big endian in the command line:

mpirun -np # [...]  ./nimrod -Wl -T  [...]

For the Intel compiler, the environment variable is F_UFMTENDIAN, and
it should be set to "big".

  Alternatively, many compilers now accept a big endian flag as an
option, and that is now the preferred approach.
