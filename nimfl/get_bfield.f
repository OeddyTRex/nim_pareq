      SUBROUTINE get_bfield(neq,y_lsode,b_xyz,bigr)
c
c						Intended to compute the
c						magnetic field strength
c						and presumes that x and y
c						are in the block with label ib.
c
      USE input
      USE global
      USE dumpc
      USE cell_type_mod
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: neq
      REAL(r8), DIMENSION(neq+3), INTENT(INOUT) :: y_lsode 
      REAL(r8), DIMENSION(3), INTENT(OUT) :: b_xyz  
      REAL(r8), INTENT(OUT) :: bigr

      REAL(r8) :: jacobian,norm_factor
      REAL(r8) :: x,y,r,z,gcosv,gsinv,angle,rnew,znew
      REAL(r8) :: btol = 1.e-11
      INTEGER(i4) :: ib,kfactor,imodes,isearch,jsearch
      TYPE(location_type):: p0
      TYPE(cell_type), POINTER :: item_pre
      LOGICAL :: failure


c
c						Translate y_lsode to local
c						representaion
      r=y_lsode(1)
      z=y_lsode(2)
      x=y_lsode(neq+1)
      y=y_lsode(neq+2)
      ib=NINT(y_lsode(neq+3))
c
c-TMP                                           In rect geometry, r,z
c                                               may have periodicity
      znew=z
      rnew=r

c     revert to using input and normal correlation of x::r, y::z

      IF(gridshape == 'rect')THEN
        SELECT CASE(periodicity)
        CASE("y-dir")
          IF(z > ymax)znew=ymin+MOD(z-ymin,ymax-ymin)
          IF(z < ymin)znew=ymax+MOD(z-ymin,ymax-ymin)
        CASE("both")
          IF(z > ymax)znew=ymin+MOD(z-ymin,ymax-ymin)
          IF(z < ymin)znew=ymax+MOD(z-ymin,ymax-ymin)
          IF(r > xmax)rnew=xmin+MOD(r-xmin,xmax-xmin)
          IF(r < xmin)rnew=xmax+MOD(r-xmin,xmax-xmin)
        END SELECT
      ENDIF

c
c						Update ib,x,y for r,z
      p0%point(1)=rnew
      p0%point(2)=znew

      isearch = NINT(search_idim*
     &              (rnew-search_rmin)/(search_rmax-search_rmin))
      jsearch = NINT(search_jdim*
     &              (znew-search_zmin)/(search_zmax-search_zmin))
      lucky => search_map(isearch,jsearch)%next
      CALL rz_to_cell(p0,lucky,item_pre,failure,.FALSE.)
c     CALL refine_cell(p0,lucky,failure,x,y,.TRUE.)
      CALL refine_cell(p0,lucky,failure,x,y,.FALSE.)
      ib=lucky%ib
      IF(failure)THEN
c       WRITE(out_unit,*)" # FAILURE at ",p0%point(1:2)
c        CALL nim_stop("get_bfield:ERROR")
c       WRITE(out_unit,*)" # trace off grid at ",p0%point(1:2)
c       WRITE(*,*)" # trace off grid at ",p0%point(1:2)
      ENDIF
      y_lsode(neq+1)=x
      y_lsode(neq+2)=y
      y_lsode(neq+3)=ib

      IF(ib  ==  0)THEN
c
c						no change, try to keep lsode
c						looping at some point.
c       b_xyz=0
        b_xyz(1:2)=0
        b_xyz(3)=1
        IF (geom=='tor') THEN
          bigr=r
        ELSE
          bigr=1
        ENDIF
        RETURN
      ELSE IF(ib <= nrblc)THEN
c
c						Calculate the equilibrium
c						magnetic field from rblock
        call lagr_quad_eval(rbc(ib)%be_eq,x,y,0_i4)
        b_xyz=rbc(ib)%be_eq%f
c
c						Convert toroidal component
c						to cylindrical and set
c						angle.
        angle=y_lsode(3)
        IF (geom=='tor') THEN
          b_xyz(3)=b_xyz(3)/r
          bigr=r
        ELSE
          angle=twopi*angle/per_length
          bigr=1
        ENDIF
c
c
c						Sum the perturbed magnetic
c						fields over the number of
c						Fourier modes.
c	be%f(1,1)-------r-comp, 1st mode (n=0)
c	be%f(2,1)-------z-comp, 1st mode (n=0)
c	be%f(3,1)-------R*phi-comp, 1st mode
c	be%f(1,2)-------r-comp, 2nd mode (n=1)
c	be%f(2,2)-------z-comp, 2nd mode
c	be%f(3,2)-------R*phi-comp, 2nd mode
c	be%f(1,3)-------r-comp, 3rd mode (n=2)
c	
c
c       B_tilde(R,Z,phi)=sum_over_modes[B_tilde_n(R,Z)*exp( i*n*phi)
c                                +conjg(B_tilde_n(R,Z)*exp(-i*n*phi)]
c
        call lagr_quad_eval(rbc(ib)%be,x,y,0_i4)
c
        IF(nonlinear)THEN
c
c						n=0 Fourier component.
c
          b_xyz(:)=b_xyz(:)+rbc(ib)%be%f(:,1)
c
c						n/=0 Fourier components.
c
          DO imodes=2,nmodes
            gcosv=COS(angle*(imodes-1)*zperiod)
            gsinv=SIN(angle*(imodes-1)*zperiod)
            b_xyz(:)=b_xyz(:)+
     &          2*(REAL(rbc(ib)%be%f(:,imodes),r8)*gcosv-
     &             AIMAG(rbc(ib)%be%f(:,imodes)  )*gsinv)
          ENDDO
        ELSE
c
c						Presumes lin_nmodes=1 
          
          gcosv=COS(angle*lin_nmax)
          gsinv=SIN(angle*lin_nmax)
          norm_factor=2._r8
          IF(lin_nmax == 0)norm_factor=1._r8
          b_xyz(:)=b_xyz(:)+norm_factor
     &           *(REAL(rbc(ib)%be%f(:,1),r8)*gcosv-
     &             AIMAG(rbc(ib)%be%f(:,1)  )*gsinv)
        ENDIF
      ELSE
c
c						Calculate the equilibrium
c						magnetic field from tblock
         CALL tri_linear_eval(tbc(lucky%ib)%be_eq,
     &                        tbc(lucky%ib)%tgeom,
     &                   p0%point(1),p0%point(2),lucky%icell,0_i4)
        b_xyz=tbc(ib)%be_eq%f
c
c						Convert toroidal component
c						to cylindrical and set angle.
        angle=y_lsode(3)
        IF (geom=='tor') THEN
          b_xyz(3)=b_xyz(3)/r
          bigr=r
        ELSE
          angle=twopi*angle/per_length
          bigr=1
        ENDIF
c
c						Calculate the equilibrium
c						magnetic field from tblock
        CALL tri_linear_eval(tbc(lucky%ib)%be,
     &                       tbc(lucky%ib)%tgeom,
     &                   p0%point(1),p0%point(2),lucky%icell,0_i4)
c
c
        IF(nonlinear)THEN
c
c						n=0 Fourier component.
c
          b_xyz(:)=b_xyz(:)+rbc(ib)%be%f(:,1)
c
c						n/=0 Fourier components.
c
          DO imodes=2,nmodes
            gcosv=COS(angle*(imodes-1)*zperiod)
            gsinv=SIN(angle*(imodes-1)*zperiod)
            b_xyz(:)=b_xyz(:)+
     &          2*(REAL(tbc(ib)%be%f(:,imodes),r8)*gcosv-
     &             AIMAG(tbc(ib)%be%f(:,imodes)  )*gsinv)
          ENDDO
        ELSE
c
c						Presumes lin_nmodes=1 
          
          gcosv=COS(angle*lin_nmax)
          gsinv=SIN(angle*lin_nmax)
          norm_factor=2._r8
          IF(lin_nmax == 0)norm_factor=1._r8
          kfactor=0
          b_xyz(:)=b_xyz(:)+norm_factor
     &           *(REAL(tbc(ib)%be%f(:,1),r8)*gcosv-
     &             AIMAG(tbc(ib)%be%f(:,1)  )*gsinv)
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE get_bfield
