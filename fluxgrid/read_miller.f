c-----------------------------------------------------------------------
c     module 1. read_miller.
c     reads data from Bob Miller's TOQ equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_miller(gt)
      USE local
      USE inverse
      USE input
      USE global
      USE physdat
      IMPLICIT NONE

      TYPE(global_type), INTENT(OUT) :: gt

      REAL(r4), DIMENSION(:), ALLOCATABLE :: psis_in,fs_in,ps_in,qs_in,
     &                                       mach_in
      REAL(r4), DIMENSION(:,:), ALLOCATABLE :: rg_in,zg_in
c-----------------------------------------------------------------------
c     the reals and integers are assumed to be 32 bit.
c-----------------------------------------------------------------------
      CALL open_bin(eq1_unit,TRIM(filename),"UNKNOWN","REWIND",32_i4)
      READ(eq1_unit)ii%mtau,ii%ma
      ii%mtau=ii%mtau-1
      ii%ma=ii%ma-1
c-----------------------------------------------------------------------
c     allocate and read arrays, close files.
c-----------------------------------------------------------------------
      ALLOCATE(psis_in(0:ii%ma))
      ALLOCATE(fs_in(0:ii%ma))
      ALLOCATE(ps_in(0:ii%ma))
      ALLOCATE(qs_in(0:ii%ma))
      ALLOCATE(rg_in(0:ii%mtau,0:ii%ma))
      ALLOCATE(zg_in(0:ii%mtau,0:ii%ma))
      ALLOCATE(mach_in(0:ii%ma))
      READ(eq1_unit)psis_in
      READ(eq1_unit)fs_in
      READ(eq1_unit)ps_in
      READ(eq1_unit)qs_in
      READ(eq1_unit)rg_in
      READ(eq1_unit)zg_in
      mach_in(:) = 0.
      READ(eq1_unit,END=999) mach_in
      CALL close_bin(eq1_unit,TRIM(filename))
 999  CONTINUE
c-----------------------------------------------------------------------
c     copy and modify 1D arrays.
c-----------------------------------------------------------------------
      ii%psio=psis_in(ii%ma)-psis_in(0)
      ii%ma=ii%ma-1
      CALL spline_alloc(ii%sq_in,ii%ma,4_i4)
      ii%sq_in%xs(:)=(psis_in(1:ii%ma+1)-psis_in(0))/ii%psio
      ii%sq_in%fs(:,1)=fs_in(1:ii%ma+1)
      ii%sq_in%fs(:,2)=ps_in(1:ii%ma+1)*mu0
      ii%sq_in%fs(:,3)=qs_in(1:ii%ma+1)
      ii%sq_in%fs(:,4)=mach_in(1:ii%ma+1)
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      ii%ro=rg_in(0,0)
      ii%zo=zg_in(0,0)
      ALLOCATE(ii%rg(0:ii%mtau,0:ii%ma))
      ALLOCATE(ii%zg(0:ii%mtau,0:ii%ma))
      ii%rg=rg_in(:,1:ii%ma+1)
      ii%zg=zg_in(:,1:ii%ma+1)
c-----------------------------------------------------------------------
c     set flags.
c-----------------------------------------------------------------------
      ii%p1flag=.FALSE.
      ii%symflag=.FALSE.
c-----------------------------------------------------------------------
c     nullify eigenvector pointer, to note that this file type does
c     not contain eigenvector data.
c-----------------------------------------------------------------------
      NULLIFY(ii%eigvec)
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL process_eq(gt)
c-----------------------------------------------------------------------
c     terminate program.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_miller
