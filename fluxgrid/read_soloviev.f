c-----------------------------------------------------------------------
c     subprogram 1. readeq.
c     reads data from Soloviev equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_soloviev(gt)
      USE local
      USE direct
      USE input
      USE global
      USE physdat
      IMPLICIT NONE

      TYPE(global_type), INTENT(OUT) :: gt

      INTEGER(i4) :: mr,mz,ma,ir,iz
c-----------------------------------------------------------------------
c     read data.
c-----------------------------------------------------------------------
      CALL open_bin(eq1_unit,"sol.dat","OLD","REWIND",64_i4)
      READ(eq1_unit)mr,mz,ma
      READ(eq1_unit)di%rmin,di%rmax,di%zmin,di%zmax
      CALL bicube_alloc(di%psig,mr,mz,1_i4)
      CALL spline_alloc(di%sq_in,ma,4_i4)
      READ(eq1_unit)di%psig%fs,di%sq_in%xs,di%sq_in%fs(:,1:2)
      CALL close_bin(eq1_unit,"sol.dat")
c-----------------------------------------------------------------------
c     modify data.
c-----------------------------------------------------------------------
      di%psio=di%sq_in%xs(ma)
      di%sq_in%xs=di%sq_in%xs/di%psio
      di%sq_in%fs(:,3:4)=0
      di%ro=0
      di%zo=0
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
c     terminate program.
c-----------------------------------------------------------------------
      RETURN
      END
