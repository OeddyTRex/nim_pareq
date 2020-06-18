c-----------------------------------------------------------------------
c     file eq_swap.
c     transfers data from 'equilibrium' fields to the n=0 part of the
c     solution arrays.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. eq_swap.
c     2. replace_coil.
c-----------------------------------------------------------------------
c     subprogram 1.  eq_swap
c     when requested, move the equilibrium fields to the n=0 part of the
c     solution vector to use them as initial conditions.
c
c     this version of eq_swap also transfers the number density, apart
c     from the ndns constant parameter if ndns is positive.  similarly,
c     rbph is the value of R*B_phi left in the equilibrium arrays.
c-----------------------------------------------------------------------
      SUBROUTINE eq_swap(nmodes,keff,geom,eq_flow,ndns,rbph)
      USE local
      USE fields
      USE lagr_quad_mod
      USE physdat
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff
      REAL(r8), INTENT(IN) :: ndns,rbph
      CHARACTER(*), INTENT(IN) :: geom,eq_flow

      INTEGER(i4) :: ibl,imode,ix,iy,ibasis,poly_degree=1
      INTEGER(i4) :: ix0,iy0
      REAL(r8) :: dx,dy,x,y
      COMPLEX(r8), DIMENSION(1) :: lpq,lpe,lte,lti,lnd
      COMPLEX(r8), DIMENSION(3) :: lbq,lvq
      LOGICAL :: n0_found=.false.
c-----------------------------------------------------------------------
c     loop over modes to find the n=0
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (keff(imode)==0) THEN
          n0_found=.true.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT.n0_found) CALL nim_stop("Cannot find n=0 Fourier comp.")
c-----------------------------------------------------------------------
c     determine the degree of the polynomials.
c-----------------------------------------------------------------------
      IF (nrbl>0) poly_degree=rb(1)%be%n_side+1
c-----------------------------------------------------------------------
c     loop over blocks, transfer data, and zero-out eq arrays.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        basis: DO ibasis=1,poly_degree**2
          ix0=rb(ibl)%be%ix0(ibasis)
          iy0=rb(ibl)%be%iy0(ibasis)
          dx=rb(ibl)%be%dx(ibasis)
          dy=rb(ibl)%be%dy(ibasis)
          DO iy=iy0,rb(ibl)%my
            DO ix=ix0,rb(ibl)%mx
              x=ix-ix0+dx
              y=iy-iy0+dy
              CALL lagr_quad_eval(rb(ibl)%be_eq,x,y,0_i4)
              CALL lagr_quad_eval(rb(ibl)%ve_eq,x,y,0_i4)
              CALL lagr_quad_eval(rb(ibl)%pres_eq,x,y,0_i4)
              CALL lagr_quad_eval(rb(ibl)%prese_eq,x,y,0_i4)
              CALL lagr_quad_eval(rb(ibl)%nd_eq,x,y,0_i4)
              lbq=rb(ibl)%be_eq%f
              lbq(3)=lbq(3)-rbph
              lvq=rb(ibl)%ve_eq%f
              lpq=rb(ibl)%pres_eq%f(1)
              lpe=rb(ibl)%prese_eq%f(1)
              lnd=rb(ibl)%nd_eq%f(1)
              lte=lpe/(lnd*kboltz)
              lti=(lpq-lpe)/(lnd*kboltz)
              IF (geom=='tor') THEN
                CALL lagr_quad_eval(rb(ibl)%rz,x,y,0_i4)
                IF (rb(ibl)%rz%f(1)==0) THEN
                  lbq(3)=0
                ELSE
                  lbq(3)=lbq(3)/rb(ibl)%rz%f(1)
                ENDIF
              ENDIF
              CALL lagr_quad_basis_assign_loc(rb(ibl)%be,lbq,ibasis,
     $                                        ix,iy,imode)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%be_eq,
     $                                (/0._r8,0._r8,rbph/),ibasis,ix,iy)
              IF (eq_flow/="none") 
     $          CALL lagr_quad_basis_assign_loc(rb(ibl)%ve,lvq,ibasis,
     $                                          ix,iy,imode)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%pres,lpq,ibasis,
     $                                        ix,iy,imode)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%prese,lpe,ibasis,
     $                                        ix,iy,imode)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%tele,lte,ibasis,
     $                                        ix,iy,imode)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%tion,lti,ibasis,
     $                                        ix,iy,imode)
              IF (ndns>0._r8) THEN
                lnd=lnd-ndns
                CALL lagr_quad_basis_assign_loc(rb(ibl)%nd,lnd,ibasis,
     $                                          ix,iy,imode)
              ENDIF
            ENDDO
          ENDDO
        ENDDO basis
        rb(ibl)%ja_eq=0
        rb(ibl)%ve_eq=0
        rb(ibl)%pres_eq=0
        rb(ibl)%prese_eq=0
        IF (ndns>0._r8) rb(ibl)%nd_eq=ndns
      ENDDO
c-----------------------------------------------------------------------
c-PRE tblocks, linear only for now.
c-----------------------------------------------------------------------
      DO ibl=nrbl+1,nbl
        tb(ibl)%be%fs(1:2,:,:,imode)=tb(ibl)%be_eq%fs(1:2,:,:)
        IF (geom=='tor') THEN
          DO ix=0,tb(ibl)%mvert
            IF (tb(ibl)%tgeom%xs(ix)==0) THEN
              tb(ibl)%be%fs(3,ix,0,imode)=0
            ELSE
              tb(ibl)%be%fs(3,ix,0,imode)=
     $          (tb(ibl)%be_eq%fs(3,ix,0)-rbph)/tb(ibl)%tgeom%xs(ix)
            ENDIF
          ENDDO
        ELSE
          tb(ibl)%be%fs(3,:,:,imode)=(tb(ibl)%be_eq%fs(3,:,:)-rbph)
        ENDIF
        tb(ibl)%be_eq%fs(1:2,:,:)=0
        tb(ibl)%be_eq%fs(3,:,:)=rbph
        tb(ibl)%ja_eq%fs=0
        tb(ibl)%ve%fs(:,:,:,imode)=tb(ibl)%ve_eq%fs
        tb(ibl)%ve_eq%fs=0
        tb(ibl)%pres%fs(:,:,:,imode)=tb(ibl)%pres_eq%fs
        tb(ibl)%pres_eq%fs=0
        tb(ibl)%prese%fs(:,:,:,imode)=tb(ibl)%prese_eq%fs
        tb(ibl)%prese_eq%fs=0
        tb(ibl)%tele%fs(:,:,:,imode)=tb(ibl)%prese%fs(:,:,:,imode)/
     $    (tb(ibl)%nd_eq%fs(:,:,:)*kboltz)
        tb(ibl)%tion%fs(:,:,:,imode)=
     $    (tb(ibl)%pres%fs(:,:,:,imode)-tb(ibl)%prese%fs(:,:,:,imode))/
     $    (tb(ibl)%nd_eq%fs(:,:,:)*kboltz)
        IF (ndns>0._r8) THEN
          tb(ibl)%nd%fs(:,:,:,imode)=tb(ibl)%nd_eq%fs(:,:,:)-ndns
          tb(ibl)%nd_eq%fs(:,:,:)=ndns
        ENDIF
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE eq_swap
c-----------------------------------------------------------------------
c     subprogram 1.  replace_coil
c     this routine adds external coils fields to the equilibrium B
c     and removes them from the n=0 part of the solution.  it is
c     intended for use after an equilibrium transfer, and can be used
c     to keep nimrod from altering the coil currents due to resistive-
c     wall diffusion.
c
c     passing the coil information allows greater flexibility.
c     note that this is intended for toroidal geometry only.
c-----------------------------------------------------------------------
      SUBROUTINE replace_coil(nmodes,keff,ncoil,coilr,coilz,coili)
      USE local
      USE fields
      USE lagr_quad_mod
      USE physdat
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nmodes,ncoil
      REAL(r8), DIMENSION(nmodes), INTENT(IN) :: keff
      REAL(r8), DIMENSION(ncoil), INTENT(IN) :: coilr,coilz,coili

      INTEGER(i4) :: ibl,imode,ix,iy,ibasis,poly_degree=1
      INTEGER(i4) :: ix0,iy0
      REAL(r8) :: dx,dy,x,y
      COMPLEX(r8), DIMENSION(1) :: lpq,lpe,lte,lti,lnd
      COMPLEX(r8), DIMENSION(3) :: lbq,lvq
      LOGICAL :: n0_found=.false.

      REAL(r8), DIMENSION(2) :: xxyy
      REAL(r8), DIMENSION(3) :: bcl
      INTEGER(i4) :: icoil
c-----------------------------------------------------------------------
c     loop over components to find n=0.
c-----------------------------------------------------------------------
      DO imode=1,nmodes
        IF (keff(imode)==0) THEN
          n0_found=.true.
          EXIT
        ENDIF
      ENDDO
      IF (.NOT.n0_found) CALL nim_stop("Cannot find n=0 Fourier comp.")
c-----------------------------------------------------------------------
c     determine the degree of the polynomials.
c-----------------------------------------------------------------------
      IF (nrbl>0) poly_degree=rb(1)%be%n_side+1
c-----------------------------------------------------------------------
c     loop over blocks.
c-----------------------------------------------------------------------
      DO ibl=1,nrbl
        basis: DO ibasis=1,SIZE(rb(ibl)%be%ix0)
          ix0=rb(ibl)%be%ix0(ibasis)
          iy0=rb(ibl)%be%iy0(ibasis)
          dx=rb(ibl)%be%dx(ibasis)
          dy=rb(ibl)%be%dy(ibasis)
          DO iy=iy0,rb(ibl)%my
            DO ix=ix0,rb(ibl)%mx
              x=ix-ix0+dx
              y=iy-iy0+dy
              CALL lagr_quad_eval(rb(ibl)%rz,x,y,1_i4)
              CALL lagr_quad_eval(rb(ibl)%be_eq,x,y,0_i4)
              CALL lagr_quad_eval(rb(ibl)%be,x,y,0_i4)
              bcl=0._r8
              CALL brz_eval(bcl(1:2),rb(ibl)%rz%f(1),rb(ibl)%rz%f(2),
     $                      ncoil,coilr,coilz,coili,mu0)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%be_eq,
     $                                        rb(ibl)%be_eq%f+bcl,
     $                                        ibasis,ix,iy)
              CALL lagr_quad_basis_assign_loc(rb(ibl)%be,
     $                                        rb(ibl)%be%f(:,imode)-bcl,
     $                                        ibasis,ix,iy,imode)
            ENDDO
          ENDDO
        ENDDO basis
      ENDDO
c-----------------------------------------------------------------------
c-PRE tblocks
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE replace_coil
