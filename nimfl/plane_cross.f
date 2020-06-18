c-----------------------------------------------------------------------
c     file plane_cross.f
c     subroutine that checks if a field line crosses a specified plane.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE plane_cross(r,z,per,cross,rc,zc,perc,perc_adj)
      USE local
      USE input
      USE input0
      IMPLICIT NONE

      REAL(r8), DIMENSION(2), INTENT(IN) :: r,z,per
      REAL(r8), INTENT(OUT) :: rc,zc,perc,perc_adj
      LOGICAL, INTENT(OUT) :: cross

      REAL(r8), DIMENSION(2) :: per_adj,delta,per_trunc
      REAL(r8) :: per_dim,trunc,zadj
      INTEGER(i4) :: itrunc,numz
c-----------------------------------------------------------------------
c     initialize.
c-----------------------------------------------------------------------
      cross=.FALSE.
      IF (geom=='tor') THEN
        per_dim=twopi
      ELSE
        per_dim=per_length
      ENDIF
c-----------------------------------------------------------------------
c     select the appropriate type of plane, check for crossing, and
c     if so, save crossing position.
c-----------------------------------------------------------------------
      SELECT CASE(plane_type)
      CASE("r")
        delta=plane_position-r
        IF (PRODUCT(delta)<0) THEN
          cross=.TRUE.
          delta(2)=r(2)-r(1)
          rc=plane_position
          zc=z(1)+delta(1)*(z(2)-z(1))/delta(2)
          perc=per(1)+delta(1)*(per(2)-per(1))/delta(2)
          perc_adj=MODULO(perc,per_dim)
          IF (perc_adj<0) perc_adj=perc_adj+per_dim
        ENDIF
      CASE("z")
        IF (gridshape=='rect'.AND.periodicity=='y-dir') THEN
          numz=(0.5_r8*(z(1)+z(2))-ymin)/(ymax-ymin)
          zadj=numz*(ymax-ymin)+plane_position
          delta=zadj-z
        ELSE
          delta=plane_position-z
          zadj=plane_position
        ENDIF
        IF (PRODUCT(delta)<0) THEN
          cross=.TRUE.
          delta(2)=z(2)-z(1)
          rc=r(1)+delta(1)*(r(2)-r(1))/delta(2)
          zc=zadj
          perc=per(1)+delta(1)*(per(2)-per(1))/delta(2)
          perc_adj=MODULO(perc,per_dim)
          IF (perc_adj<0) perc_adj=perc_adj+per_dim
        ENDIF
      CASE("periodic")
        per_trunc=INT(per/per_dim)
        per_adj=per-per_trunc(1)*per_dim
        IF (per_adj(1)<0) per_adj=per_adj+per_dim
        delta=plane_position-per_adj
        IF (PRODUCT(delta)<0) THEN
          cross=.TRUE.
          delta(2)=per(2)-per(1)
          rc=r(1)+delta(1)*(r(2)-r(1))/delta(2)
          zc=z(1)+delta(1)*(z(2)-z(1))/delta(2)
          perc=plane_position
          perc_adj=plane_position
        ENDIF
        IF (.NOT.cross) THEN
          per_adj=per-per_trunc(2)*per_dim
          IF (per_adj(2)<0) per_adj=per_adj+per_dim
          delta=plane_position-per_adj
          IF (PRODUCT(delta)<0) THEN
            cross=.TRUE.
            delta(2)=per(2)-per(1)
            rc=r(1)+delta(1)*(r(2)-r(1))/delta(2)
            zc=z(1)+delta(1)*(z(2)-z(1))/delta(2)
            perc=plane_position
            perc_adj=plane_position
          ENDIF
        ENDIF
      CASE DEFAULT
        CALL nim_stop("Inappropriate plane_type: "//TRIM(plane_type))
      END SELECT
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE plane_cross
