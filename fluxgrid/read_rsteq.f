      SUBROUTINE read_rsteq(gt) 
c-----------------------------------------------------------------------
c     reads data from rsteq equilibrium code.
c-----------------------------------------------------------------------
      USE local
      USE physdat
      USE direct
      USE input
      USE global
      IMPLICIT NONE
      
      TYPE(global_type), INTENT(OUT) :: gt

      INTEGER(i4) :: i,j,nw,nh,ia,mr,mz,ma,ir,iz
      INTEGER(i4) :: nrst
c-----------------------------------------------------------------------
c     read equilibrium data.
c-----------------------------------------------------------------------
      OPEN(UNIT=eq1_unit,FILE=TRIM(filename),STATUS='old',
     $     FORM='UNFORMATTED')
      READ(eq1_unit)nw,nh,nrst
      READ(eq1_unit)di%rmin,di%rmax,di%zmin,di%zmax
      READ(eq1_unit)di%ro,di%zo
      mr=nw-1
      mz=2*nh-3
      ma=nrst-1
      di%psig%mx=mr
      di%psig%my=mz
      CALL bicube_alloc(di%psig,mr,mz,1_i4)
      READ(eq1_unit)((di%psig%fs(1,i,j),i=0,nw-1),j=nh-2,mz)
      DO i=0,nw-1
        DO j=0,nh-3
          di%psig%fs(1,i,j)=di%psig%fs(1,i,mz-j)
        ENDDO
      ENDDO
      CALL spline_alloc(di%sq_in,nrst-1_i4,3_i4)
      READ(eq1_unit)(di%sq_in%fs(i,1),i=0,nrst-1)   ! psi
      di%psio=di%sq_in%fs(0,1)-di%sq_in%fs(nrst-1,1)
      READ(eq1_unit)(di%sq_in%fs(i,1),i=0,nrst-1)   ! F
      READ(eq1_unit)(di%sq_in%fs(i,2),i=0,nrst-1)   ! pressure
      READ(eq1_unit)(di%sq_in%fs(i,3),i=0,nrst-1)   ! q
      di%sq_in%fs(:,4)=0			    ! mach
      CLOSE(UNIT=eq1_unit)
c-----------------------------------------------------------------------
c     translate to internal quantities.
c-----------------------------------------------------------------------
      di%sq_in%xs=(/(ia,ia=0,ma)/)/dfloat(ma)
      di%sq_in%fs(:,1)=ABS(di%sq_in%fs(:,1))
      di%sq_in%fs(:,2)=MAX(di%sq_in%fs(:,2)*mu0,0._r8)
c-----------------------------------------------------------------------
c     normalize sign convention.
c-----------------------------------------------------------------------
      IF(di%psio < 0)THEN
         di%psio=-di%psio
         di%psig%fs=-di%psig%fs
      ENDIF
c-----------------------------------------------------------------------
c     copy axis quantities.
c-----------------------------------------------------------------------
      gt%ro=di%ro
      gt%zo=di%zo
      gt%psio=di%psio
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
      END
