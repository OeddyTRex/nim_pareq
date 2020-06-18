c-----------------------------------------------------------------------
c     file input0.f
c     declarations and default values of input parameters.
c     note that new parameters should be added to the lists in module
c     input and to the namelist declarations in subroutine 
c     read_namelist.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  input0.
c     2.  read_namelist0.
c-----------------------------------------------------------------------
c     subprogram 1. input0
c     module containing the type declaration and default values for
c     nimfl input variables.
c-----------------------------------------------------------------------
      MODULE input0
      USE local
      IMPLICIT NONE
c=======================================================================
c     nimfl parameters.
c=======================================================================
      INTEGER(i4) :: nimfl_step=1000 ! number of line segment integrations.
      REAL(r8) :: per_start=0           ! starting periodic coordinate
      REAL(r8) :: frac_step=0.1 ! fraction of characteristic L per step.        
      LOGICAL :: same_color=.FALSE.

      CHARACTER(16) :: plane_type="periodic"    !  or "r" or "z"
      REAL(r8) :: plane_position=0

      LOGICAL :: out_tecplot=.FALSE.
      CHARACTER(8) :: task='poincare'           ! 'poincare','ntm','poloidal'
      CHARACTER(4) :: poincare_positions='qval' ! 'read','qval','cell','axis'
      INTEGER(i4) :: nn = 1     ! Toroidal mode number and below for resonance.
      INTEGER(i4) :: mm = 4     ! Max mode number for poloidal FFT evaluation.
      INTEGER(i4) :: num_surface= 20     ! Number of surfaces for poloidal FFT.
      CHARACTER(8) :: find_axis='default'! 'default','read'
      REAL(r8) :: rmguess=0,zmguess=0
      REAL(r8) :: cross_tol=1.e-4_r8	! plane-crossing tolerance
					! absolute--may need adjusting!
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE input0

c-----------------------------------------------------------------------
c     subprogram 1. read_namelist0
c     open and read the namelist input.
c-----------------------------------------------------------------------

      SUBROUTINE read_namelist0(infile,echo_in)
      USE local
      USE input0
      USE input
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: infile
      LOGICAL, INTENT(IN) :: echo_in

      INTEGER :: read_stat
c=======================================================================
      NAMELIST/nimfl_input/per_length,nimfl_step
     $     ,same_color,plane_type,plane_position,per_start,frac_step
     $     ,out_tecplot,poincare_positions,nn,mm,num_surface,task
     $     ,dump_file
     $     ,find_axis,rmguess,zmguess,cross_tol
c-----------------------------------------------------------------------
c     open input file.
c-----------------------------------------------------------------------
      OPEN(UNIT=in_unit,FILE=infile,STATUS='OLD',POSITION='REWIND')
c-----------------------------------------------------------------------
c     read namelist.
c-----------------------------------------------------------------------
      READ(UNIT=in_unit,NML=nimfl_input,IOSTAT=read_stat)
      IF (read_stat/=0) CALL nim_stop
     $   ('Error reading nimfl input file.')
c-----------------------------------------------------------------------
c     close input file and terminate.
c-----------------------------------------------------------------------
      CLOSE(in_unit)
      END SUBROUTINE read_namelist0
