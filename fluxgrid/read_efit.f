c-----------------------------------------------------------------------
c     module 1. read_efit.
c     reads data from General Atomic's EFIT equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_efit(gt)
      USE local
      USE direct
      USE input
      USE global
      USE physdat
      IMPLICIT NONE
      
      TYPE(global_type), INTENT(OUT) :: gt

      INTEGER(i4) :: i,j,nw,nh,ia,mr,mz,ma, ir,iz
      INTEGER :: ios
      REAL(r8) :: bcentr,cpasma,rgrid,rmaxis,rzero,ssibry1,
     $     ssibry2,ssimag1,ssimag2,xdim,xdum,zdim,zmaxis,zmid
      INTEGER(i4) :: nbbbs, limitr
      REAL(r8), DIMENSION(:), ALLOCATABLE :: rbbbs,zbbbs
      REAL(r8), DIMENSION(:), ALLOCATABLE :: xlim,ylim
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      OPEN(UNIT=eq1_unit,FILE=TRIM(filename),STATUS='old')
      READ(eq1_unit,'(52x,2i4)')nw,nh
      READ(eq1_unit,'(5e16.9)')xdim,zdim,rzero,rgrid,zmid
      READ(eq1_unit,'(5e16.9)')rmaxis,zmaxis,ssimag1,ssibry1,bcentr
      READ(eq1_unit,'(5e16.9)')cpasma,ssimag2,xdum,rmaxis,xdum
      READ(eq1_unit,'(5e16.9)')zmaxis,xdum,ssibry2,xdum,xdum
      CALL spline_alloc(di%sq_in,nw-1_i4,4_i4)
      READ(eq1_unit,'(5e16.9)')(di%sq_in%fs(i,1),i=0,nw-1)
      READ(eq1_unit,'(5e16.9)')(di%sq_in%fs(i,2),i=0,nw-1)
      READ(eq1_unit,'(5e16.9)')(di%sq_in%fs(i,3),i=0,nw-1)
      READ(eq1_unit,'(5e16.9)')(di%sq_in%fs(i,3),i=0,nw-1)
      di%sq_in%fs(:,4)=0.                                       ! Mach
      CALL bicube_alloc(di%psig,nw-1_i4,nh-1_i4,1_i4)
      READ(eq1_unit,'(5e16.9)')((di%psig%fs(1,i,j),i=0,nw-1),j=0,nh-1)
      READ(eq1_unit,'(5e16.9)',iostat=ios)(di%sq_in%fs(i,3),i=0,nw-1)
c-----------------------------------------------------------------------
c     read and write wall data.
c-----------------------------------------------------------------------
      READ(eq1_unit,'(2i5)',iostat=ios)nbbbs,limitr
      IF(ios.EQ.0_i4)then
        ALLOCATE(rbbbs(nbbbs))
        ALLOCATE(zbbbs(nbbbs))
        ALLOCATE(xlim(limitr))
        ALLOCATE(ylim(limitr))
        READ(eq1_unit,'(5e16.9)') (rbbbs(i),zbbbs(i),i=1,nbbbs)
        READ(eq1_unit,'(5e16.9)') (xlim(i),ylim(i),i=1,limitr)
        OPEN(UNIT=eq2_unit,file='wall.dat',status='replace')
        WRITE(eq2_unit,*)limitr
        DO i=1,limitr
          WRITE(eq2_unit,*)xlim(i),ylim(i)
        ENDDO
        CLOSE(eq2_unit)
      ENDIF
      CLOSE(UNIT=eq1_unit)
c-----------------------------------------------------------------------
c     translate to internal quantities.
c-----------------------------------------------------------------------
      mr=nw-1
      mz=nh-1
      ma=nw-1
      di%psig%mx=mr
      di%psig%my=mz
      di%rmin=rgrid
      di%rmax=rgrid+xdim
      di%zmin=-zdim/2
      di%zmax=zdim/2
      di%psio=ssibry1-ssimag1
      di%sq_in%xs=(/(ia,ia=0,ma)/)/dfloat(ma)
      di%sq_in%fs(:,1)=ABS(di%sq_in%fs(:,1))
      di%sq_in%fs(:,2)=MAX(di%sq_in%fs(:,2)*mu0,0._r8)  !Pressure
c-----------------------------------------------------------------------
c     copy and convert 2D quantities.
c-----------------------------------------------------------------------
      di%ro=rmaxis
      di%zo=zmaxis
      di%psig%fs=ssibry1-di%psig%fs
c-----------------------------------------------------------------------
c     normalize sign convention.
c-----------------------------------------------------------------------
      IF(di%psio < 0)THEN
         di%psio=-di%psio
         di%psig%fs=-di%psig%fs
      ENDIF
c-----------------------------------------------------------------------
c     fit input to cubic and bicubic splines.
c-----------------------------------------------------------------------
      CALL spline_fit(di%sq_in,"extrap")
      dr0=(di%rmax-di%rmin)/di%psig%mx
      dz0=(di%zmax-di%zmin)/di%psig%my
      di%psig%xs=di%rmin+(/(ir,ir=0,di%psig%mx)/)*dr0
      di%psig%ys=di%zmin+(/(iz,iz=0,di%psig%my)/)*dz0
      CALL bicube_fit(di%psig,"extrap","extrap")
c-----------------------------------------------------------------------
c     diagnose 1D input.
c-----------------------------------------------------------------------
      IF(out_1d .OR. bin_1d)THEN
         IF(out_1d)OPEN(UNIT=ascii_unit,FILE='1d.out',STATUS='UNKNOWN')
         IF(bin_1d)CALL open_bin(binary_unit,"1d.bin","UNKNOWN",
     $        "REWIND",32_i4)
         di%sq_in%title=(/' psi  ','  f   ','  p   ','  q   ',' mach '/)
         CALL spline_write(di%sq_in,out_1d,bin_1d,
     $        ascii_unit,binary_unit)
         IF(out_1d)CLOSE(ascii_unit)
         IF(bin_1d)CALL close_bin(binary_unit,"1d.bin")
      ENDIF
c-----------------------------------------------------------------------
c     define positions  
c-----------------------------------------------------------------------
      gt%ro=di%ro
      gt%zo=di%zo
      gt%psio=di%psio
      IF(grid_method == "original")mpsi=di%sq_in%nodes
      CALL position(gt)
c-----------------------------------------------------------------------
c     Handle the rblock (if psihigh > 0) and tblock regions
c-----------------------------------------------------------------------
      IF (psihigh > 0.) THEN 
        CALL process_eq(gt)
        CALL write_piegrid
        IF (triangle) THEN 
           CALL write_rimgrid(gt)
c           CALL hole_tri(gt)
        ENDIF
      ELSE
        IF (triangle) THEN
c            CALL pure_tri
            CALL nim_stop("Finished")
        ELSE
            CALL nim_stop("Must use either rblocks or tblocks")
        ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_efit
