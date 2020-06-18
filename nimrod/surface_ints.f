c-----------------------------------------------------------------------
c     file surface_ints.f
c     module that includes integrand routines associated with surface
c     integrals.  all subroutines included here must use the same
c     interface block.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  no_surf_int.
c     2.  e_tangential.
c     3.  mag_tension.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE surface_ints

      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE lagr_quad_mod
      USE tri_linear
      USE math_tran
      USE physdat
      USE input
      USE global
      USE fft_mod
      USE generic_evals
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. no_surf_int.
c     dummy routine with the correct interface.
c-----------------------------------------------------------------------
      SUBROUTINE no_surf_int(int,rb,tb,x,y,bigr,norm,
     $                       ijcell,alpha,dxdr,dydr,dxdz,dydz)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: int 
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(tblock_type), INTENT(INOUT) :: tb
      REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(2), INTENT(IN) :: norm
      REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
      INTEGER(i4), DIMENSION(2) :: ijcell

      int=0

      RETURN
      END SUBROUTINE no_surf_int
c-----------------------------------------------------------------------
c     subprogram 2. e_tangential.
c     compute the integrand for tangential electric field.
c-----------------------------------------------------------------------
      SUBROUTINE e_tangential(int,rb,tb,x,y,bigr,norm,
     $                        ijcell,alpha,dxdr,dydr,dxdz,dydz)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: int 
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(tblock_type), INTENT(INOUT) :: tb
      REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(2), INTENT(IN) :: norm
      REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
      INTEGER(i4), DIMENSION(2) :: ijcell

      REAL(r8) :: e_tor
      COMPLEX(r8), DIMENSION(1,SIZE(int,3)) :: tele,ter,tez
      INTEGER(i4) :: imode
c-----------------------------------------------------------------------
c     find toroidal applied E from loop voltage.
c-----------------------------------------------------------------------
      IF (geom=='tor') THEN
        e_tor=volt/(twopi*bigr)
      ELSE
        e_tor=volt/per_length
      ENDIF
c-----------------------------------------------------------------------
c     apply to n=0 only.
c-----------------------------------------------------------------------
      int=0._r8
      mode_loop: DO imode=1,nmodes
        IF (keff(imode)/=0) CYCLE mode_loop
        int(1,:,imode)=-alpha*dt*e_tor*norm(2)
        int(2,:,imode)= alpha*dt*e_tor*norm(1)
        int(3,:,imode)=-alpha*dt*e_vert*norm(1)
      ENDDO mode_loop
c-----------------------------------------------------------------------
c     for two-fluid computations, the electric field in the B-advance
c     is actually E+kboltz*grad(Te)/e, which eliminates the purely
c     electrostatic term.  This term, therefore, needs to be added to
c     the surface integrand.
c-----------------------------------------------------------------------
      IF (ohms/='mhd'.AND.beta>0) THEN
        CALL generic_eval(rb%tele,tb%tele,dxdr,dydr,dxdz,dydz,x,y,
     $                    tb%tgeom,ijcell,tele,ter,tez,1_i4)
        ter=dt*kboltz*(1._r8-meomi)/(elementary_q*(1._r8+meomi))*ter
        tez=dt*kboltz*(1._r8-meomi)/(elementary_q*(1._r8+meomi))*tez
        tele=dt*kboltz*(1._r8-meomi)/(elementary_q*(1._r8+meomi))*tele

        estat_loop: DO imode=1,nmodes
          int(1,:,imode)=int(1,:,imode)-alpha*
     $                   norm(2)*(0,1)*keff(imode)*tele(1,imode)/bigr
          int(2,:,imode)=int(2,:,imode)+alpha*
     $                   norm(1)*(0,1)*keff(imode)*tele(1,imode)/bigr
          int(3,:,imode)=int(3,:,imode)+alpha*
     $                   (ter(1,imode)*norm(2)-tez(1,imode)*norm(1))
        ENDDO estat_loop
      ENDIF
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE e_tangential
c-----------------------------------------------------------------------
c     subprogram 3. mag_tension.
c     compute the integrand for Maxwell stress tension at the 
c     surface, -dS.[BB-0.5*B^2*I].alpha .
c-----------------------------------------------------------------------
      SUBROUTINE mag_tension(int,rb,tb,x,y,bigr,norm,
     $                       ijcell,alpha,dxdr,dydr,dxdz,dydz)

      COMPLEX(r8), DIMENSION(:,:,:), INTENT(OUT) :: int 
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(tblock_type), INTENT(INOUT) :: tb
      REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(2), INTENT(IN) :: norm
      REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
      INTEGER(i4), DIMENSION(2) :: ijcell

      REAL(r8), DIMENSION(3) :: be_eq,dumr
      COMPLEX(r8), DIMENSION(3,SIZE(int,3)) :: be
      COMPLEX(r8), DIMENSION(1,1) :: dumc

      REAL(r8) :: dtmu0,nbeq
      COMPLEX(r8) :: dot,nb
      INTEGER(i4) :: im
c-----------------------------------------------------------------------
c     point-wise evaluations.
c-----------------------------------------------------------------------
      CALL generic_eval(rb%be_eq,tb%be_eq,dxdr,dydr,dxdz,dydz,x,y,
     $                  tb%tgeom,ijcell,be_eq,dumr,dumr,0_i4)
      CALL generic_eval(rb%be,tb%be,dxdr,dydr,dxdz,dydz,x,y,
     $                  tb%tgeom,ijcell,be,dumc,dumc,0_i4)
      IF (geom=='tor') be_eq(3)=be_eq(3)/bigr
      dtmu0=dt/mu0
c-----------------------------------------------------------------------
c     linear terms only so far.
c-----------------------------------------------------------------------
      int=0
      nbeq=norm(1)*be_eq(1)+norm(2)*be_eq(2)
      DO im=1,nmodes
        dot=be_eq(1)*be(1,im)+be_eq(2)*be(2,im)+be_eq(3)*be(3,im)
        nb=norm(1)*be(1,im)+norm(2)*be(2,im)
        int(1,:,im)=dtmu0*(dot*norm(1)-nbeq*be(1,im)-nb*be_eq(1))*alpha
        int(2,:,im)=dtmu0*(dot*norm(2)-nbeq*be(2,im)-nb*be_eq(2))*alpha
        int(3,:,im)=dtmu0*(           -nbeq*be(3,im)-nb*be_eq(3))*alpha
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE mag_tension
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE surface_ints
