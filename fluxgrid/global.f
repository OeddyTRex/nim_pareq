      MODULE global
      USE local
      USE spline
      USE bicube
      IMPLICIT NONE

      TYPE :: global_type
      INTEGER(i4) :: ntor,nsing,nn
      REAL(r8) :: amean,aratio,betan,betap1,betap2,betap3,betat,bt0,
     $     crnt,delta1,delta2,kappa,li1,li2,li3,rmean,ro,zo,psio,
     $     q0,qmin,qmax,qa,rs1,rs2,rquot,omega0,growth,volume
      INTEGER(i4), DIMENSION(64) :: msing
      REAL(r8), DIMENSION(64) :: qsing,q1sing,psising
      REAL(r8), DIMENSION(:), POINTER :: eigveco
      REAL(r8), DIMENSION(:,:,:), POINTER :: eigvec
      CHARACTER(7) :: eqtype
      TYPE(spline_type) :: sq,ob,rzsep
      TYPE(bicube_type) :: r2g, twod, vac, dir
      END TYPE global_type

      END MODULE global
c-----------------------------------------------------------------------
c     The dominant information about the equilibrium is stored 
c       in 3 structures calculated in process_eq in inverse.f and 
c       direct.f:
c
c     gt%sq                             Surface quantities
c       gt%sq%xs = Psi_normal
c       gt%sq%fs(:,1) = R B_T (toroidal flux function)
c       gt%sq%fs(:,2) = Pressure  (No flow)
c       gt%sq%fs(:,3) = q (Safety factor)
c       gt%sq%fs(:,4) = Mach number (Normalized toroidal flow)
c
c     gt%r2g                            Mapping information
c	gt%r2g%xs(:) = theta
c	gt%r2g%ys(:) = Psi_normal
c	gt%r2g%fs(1,:,:) = rho
c	gt%r2g%fs(2,:,:) = eta
c
c     gt%twod                           Jacobian and metric elements
c	gt%twod%fs(1,:,:) = R(theta,rho)
c	gt%twod%fs(2,:,:) = Z(theta,rho)
c	gt%twod%fs(3,:,:) = jac(theta,rho)
c	gt%twod%fs(4,:,:) = gpsipsi(theta,rho)
c	gt%twod%fs(5,:,:) = gpsitheta(theta,rho)
c	gt%twod%fs(6,:,:) = gthetatheta(theta,rho)
c
c     Other useful information:
c
c     Type of equilibria: either inverse or direct:
c       eqtype
c
c     Characterization of the plasma calculated in analyze.f:
c      amean,aratio,betan,betap1,betap2,betap3,betat,bt0,
c      crnt,delta1,delta2,kappa,li1,li2,li3,rmean,ro,zo,psio,
c      q0,qmin,qmax,qa,rs1,rs2,rquot,omega0,growth,volume
c
c     Location of singular surfaces calculated in grid.f:
c      msing, qsing,q1sing,psising
c      
c     Boundary information calculated in process_eq in inverse.f
c       and direct.f:
c       ob - R,Z points of outer boundary of rblock region (that is
c            inside separatrix.
c       rzsep - R,Z points of separatrix
c       
c     Other:
c        vac - information for setupping up rblocks in the vacuum region
c        gt%vac%xs=gt%ob%xs
c        gt%vac%ys=(/(REAL(ivac),ivac=0,mvac)/)/REAL(mvac)
c        gt%vac%fs(1,:,:) = R
c        gt%vac%fs(2,:,:) = Z
c        gt%vac%fs(3,:,:) = B_R
c        gt%vac%fs(4,:,:) = B_Z
c        gt%vac%fs(5,:,:) = B_tor
c        gt%vac%fs(6,:,:) = J_R
c        gt%vac%fs(7,:,:) = J_Z
c        gt%vac%fs(8,:,:) = J_tor
c        gt%vac%fs(9,:,:) = Pressure
c        gt%vac%fs(10,:,:) = conc
c        gt%vac%fs(11,:,:) = psi
c
c        dir
c        gt%dir%xs=gt%twod%xs
c        gt%dir%ys=gt%twod%ys
c        gt%dir%fs(1,:,:) = B_R
c        gt%dir%fs(2,:,:) = B_Z
c        gt%dir%fs(3,:,:) = B_tor
c        gt%dir%fs(4,:,:) = J_R
c        gt%dir%fs(5,:,:) = J_Z
c        gt%dir%fs(6,:,:) = J_tor
c        NOTE: Get R,Z values from gt%twod
c-----------------------------------------------------------------------
