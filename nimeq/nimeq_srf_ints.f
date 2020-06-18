c-----------------------------------------------------------------------
c     file nimeq_srf_ints.f
c     module that includes integrand routines associated with surface
c     integrals for nimeq.  all subroutines included here must use the
c     same interface block.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1.  no_surf_int.
c     2.  surface_resid.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE nimeq_srf_ints

      USE local
      USE rblock_type_mod
      USE tblock_type_mod
      USE lagr_quad_mod
      USE tri_linear
      USE math_tran
      USE physdat
      USE input
      USE input_eq
      USE global
      USE generic_evals
      IMPLICIT NONE

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. no_surf_int.
c     dummy routine with the correct interface.
c-----------------------------------------------------------------------
      SUBROUTINE no_surf_int(int,rb,tb,x,y,bigr,norm,
     $                       ijcell,alpha,dxdr,dydr,dxdz,dydz)

      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: int 
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
c     subprogram 2. surface_resid.
c     compute the integrand for the surface computation of the
c     residual integral computations in gsfree.
c-----------------------------------------------------------------------
      SUBROUTINE surface_resid(int,rb,tb,x,y,bigr,norm,
     $                         ijcell,alpha,dxdr,dydr,dxdz,dydz)

      REAL(r8), DIMENSION(:,:), INTENT(OUT) :: int 
      TYPE(rblock_type), INTENT(INOUT) :: rb
      TYPE(tblock_type), INTENT(INOUT) :: tb
      REAL(r8), INTENT(IN) :: x,y,bigr,dxdr,dydr,dxdz,dydz
      REAL(r8), DIMENSION(2), INTENT(IN) :: norm
      REAL(r8), DIMENSION(:), INTENT(IN) :: alpha
      INTEGER(i4), DIMENSION(2) :: ijcell

      REAL(r8), DIMENSION(1) :: lam,lamr,lamz
c-----------------------------------------------------------------------
c     evaluate lambda and its gradient along the surface.
c-----------------------------------------------------------------------
      CALL generic_eval(rb%rwork2,tb%rwork2,dxdr,dydr,dxdz,dydz,x,y,
     $                  tb%tgeom,ijcell,lam,lamr,lamz,1_i4)
c-----------------------------------------------------------------------
c     assemble -alpha*R**2*grad(lambda).n-hat for the first component.
c-----------------------------------------------------------------------
      int(1,:)=-alpha*bigr**2*(lamr(1)*norm(1)+lamz(1)*norm(2))
      int(2,:)=0._r8
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE surface_resid
c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE nimeq_srf_ints
