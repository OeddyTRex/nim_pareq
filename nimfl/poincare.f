      SUBROUTINE poincare
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE start_positions
      USE pardata
      IMPLICIT NONE

      REAL(r8) :: d_ind,x_start,y_start,rmax
      LOGICAL :: any_cross,failure
      INTEGER(i4) :: ix,iy,nlines=0,iline=0,ibl,yend,xst
      INTEGER(i4) :: maxline=1000
      INTEGER(i4) :: istart,ntraces,rem
      TYPE(location_type), POINTER :: p0
      TYPE(cell_type), POINTER :: item
      CHARACTER(16) :: nimfile,binfile
c-----------------------------------------------------------------------
c     open xdraw file for the poincare plot.
c-----------------------------------------------------------------------
      nimfile="nimfl"
      IF (nprocs>1) THEN
        IF (node<10) THEN
          WRITE(nimfile,'(a,i1)') "nimfl"//"00", node
        ELSE IF (node<100) THEN
          WRITE(nimfile,'(a,i2)') "nimfl"//"0", node
        ELSE
          WRITE(nimfile,'(a,i3)') "nimfl", node
        ENDIF
      ENDIF
      binfile=TRIM(nimfile)//".bin"
      CALL open_bin(temp_unit,TRIM(binfile),"UNKNOWN","REWIND",32_i4)
      nimfile=TRIM(nimfile)//".dat"
      IF (out_tecplot)
     $  OPEN(UNIT=tec2d,FILE=TRIM(nimfile),STATUS='UNKNOWN')
c-----------------------------------------------------------------------
c     find rmax.
c-----------------------------------------------------------------------
      rmax=0
      DO ibl=1,nrblc
        rmax=MAX(rmax,MAXVAL(rbc(ibl)%rz%fs(:,:,1)))
      ENDDO
c-----------------------------------------------------------------------
c     set the independent variable increment (phi for toroidal geom,
c     axial length for linear geom).  the parameter frac_step is the
c     fraction of a characteristic length to be integrated at each
c     first-level integration step.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        d_ind=frac_step*nimfl_step*rmax*twopi
      ELSE
        d_ind=frac_step*nimfl_step*per_length
      ENDIF
c-----------------------------------------------------------------------
c     perform integrations.
c-----------------------------------------------------------------------
      ntraces=n_fieldlines/nprocs
      rem=n_fieldlines-nprocs*ntraces
      IF (node<rem) THEN
        ntraces=ntraces+1
        istart=node*ntraces+1
      ELSE
        istart=rem*(ntraces+1)+(node-rem)*ntraces+1
      ENDIF

      DO iline=istart,istart+ntraces-1
        IF (out_tecplot)WRITE(tec2d,*)" ZONE"
        any_cross=.FALSE.
        ibl=ibl_fieldlines(iline)
        x_start=x_fieldlines(iline)
        y_start=y_fieldlines(iline)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        WRITE(nim_wr,'(/,2(a,i3),a,2(x,E12.5))')
     &        'Starting line ',iline,' of ',n_fieldlines,
     &        ' at ',rbc(ibl)%rz%f
        CALL onefl(d_ind,ibl,x_start,y_start,
     &                   per_start,nimfl_step,plane_type,any_cross)
        IF (.NOT.same_color.AND.any_cross) WRITE(temp_unit)
      ENDDO
      IF (same_color) WRITE(temp_unit)
c-----------------------------------------------------------------------
c     close the xdraw file.
c-----------------------------------------------------------------------
      CALL close_bin(temp_unit,TRIM(binfile))
      IF(out_tecplot)CLOSE(UNIT=tec2d)
      END SUBROUTINE poincare
