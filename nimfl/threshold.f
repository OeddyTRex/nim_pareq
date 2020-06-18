      SUBROUTINE threshold(found)
c
c			Find resonant surfaces and compute NTM threshold
      USE local
      USE io
      USE input0
      USE global
      USE fields
      USE dump
      USE cell_type_mod
      USE dumpc
      USE magnetic_axis
      USE node_type_mod
      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: found

      LOGICAL :: local_debug=.FALSE.
      TYPE(cell_type), POINTER :: item
      INTEGER(i4) :: ibl,iline=0,maxline=1000,maxstep=5000,nline
      INTEGER(i4) :: mvalue,nvalue,i_bisect,max_bisect = 100,i_repeat
      REAL(r8) :: x_start,y_start,rmax
      REAL(r8), PARAMETER :: ind_tol0 = 1.e-10
      REAL(r8) :: qvalue,area,qtest,qvalue_last
      REAL(r8) :: r_min,z_min,r_max,z_max
      REAL(r8) :: x_left,x_right,q_left,q_right
      REAL(r8) :: qmin,qmax,area_last = 0
      REAL(r8), DIMENSION(:,:),ALLOCATABLE :: qprof
      TYPE(location_type) :: p0,p0_last
      REAL(r8) :: q_search_value
      REAL(r8) :: rs_eff,eps_eff,rzero,shear_eff
      REAL(r8) :: pres_eq_last,pres_eq_new,beta_prime,bmag,bigr
      REAL(r8) :: w_nc,w_d
      REAL(r8) :: mu_e_over_nu_e
      REAL(r8) :: deltaprime,determinant
      REAL(r8) :: nu_e,gamma_max,tau_a

      INTEGER(i4), PARAMETER :: neq = 4
      REAL(r8), DIMENSION(neq+3) :: y_lsode
      REAL(r8), DIMENSION(3) :: b0_xyz
      REAL(r8), PARAMETER :: mu0=pi*4.e-7_r8 
      REAL(r8), PARAMETER :: zeff = 1.0_r8
      REAL(r8), PARAMETER :: elementary_q=1.60217733e-19_r8
      REAL(r8), PARAMETER, DIMENSION(2) ::
     $          qs=(/ -elementary_q, elementary_q /),
     $          ms=(/ 9.1093898e-31_r8, 3.3435860e-27_r8 /)



      found = .FALSE.
      NULLIFY(surface_first)
      qmin=HUGE(1.0)
      qmax=-HUGE(1.0)
c
c                                       Calculate Electron collision frequency
      nu_e=3.553e-14*ndens*elecd/ABS(zeff)
c
c                                       Calculate Alfven time
             y_lsode(1)=rmaxis              ! R coordinate
             y_lsode(2)=zmaxis              ! Z coordinate
             y_lsode(3)=0.                  ! PEST angle
             y_lsode(4)=0.                  ! Area integration
             y_lsode(neq+1)=0.              ! x coordinate in logical space
             y_lsode(neq+2)=0.              ! y coordinate in logical space
             y_lsode(neq+3)=1               ! block label
             CALL get_bfield0(neq,y_lsode,b0_xyz,bigr)
             bmag=SQRT(SUM(b0_xyz**2))+TINY(bmag)
             tau_a=(bmag/rmaxis)**2/(mu0*ndens*(ms(1)+ms(2)))
             tau_a=1._r8/tau_a**0.5



c
c					Number of surfaces.
      item => start
      DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO
      nline=iline
      ALLOCATE(qprof(nline,2))

c
c     					Set the coarse mesh
      iline=0
      item => start
      outer : DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        IF(local_debug)WRITE(7,*)"# Loop ",iline
        IF(local_debug)WRITE(8,*)"# Loop ",iline
c
c     					set the start positions
        ibl = item%ib
        x_start = item%p(1,1)
        y_start = item%p(2,1)
        CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
        p0%point(1:2)=rbc(ibl)%rz%f		! (R,Z) coordinate
c
c     					Find qvalue
        CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
        WRITE(nim_wr,*)iline,area,qvalue
        IF(area < area_last)THEN		! How can this be?
          nline = iline - 1
          IF(nline < 1)THEN
            iline=iline-1
            item => item%face(1)%p
            CYCLE
          ENDIF
          EXIT outer
        ENDIF
        IF(area < 0) THEN		! Trap for hitting the wall.
          nline=iline-1
          IF(nline < 1)THEN
            iline=iline-1
            item => item%face(1)%p
            CYCLE
          ENDIF
          WRITE(nim_wr,*)"This n=0 fieldline intersects the wall."
          WRITE(nim_wr,*)"The calculation continues short of the wall."
          EXIT outer
        ENDIF
        area_last = area
        qprof(iline,1)=area
        qprof(iline,2)=qvalue
        IF(qvalue > qmax)qmax=qvalue
        IF(qvalue < qmin)qmin=qvalue
        

c
c						Advance to new start position
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO outer
      WRITE(nim_wr,*)'Range of q values',qmin,' to ',qmax
c-----------------------------------------------------------------
c						Loop through qvalues
      DO nvalue = 1,nn
       mloop: DO mvalue = FLOOR(qmin)*nvalue,CEILING(qmax)*nvalue
          q_search_value = REAL(mvalue,r8)/REAL(nvalue,r8)
          DO i_repeat=1,nvalue-1
            IF(MOD(q_search_value,REAL(i_repeat,r8)) == 0)CYCLE mloop
          ENDDO
c-----------------------------------------------------------------
c
c						Use the coarse grid to
c						bracket a local search.
      iline=0
      item => start
      DO WHILE(ASSOCIATED(item%face(1)%p))
        iline=iline+1
        IF(iline+1 == nline) EXIT

        qtest =(qprof(iline,2)-q_search_value)
     &         *(qprof(iline+1,2)-q_search_value)
c
c						Apply bisection rule to refine
c						for the starting positions.
        IF(qtest < 0 ) THEN
          IF(local_debug)THEN
            WRITE(nim_wr,*)iline,item%id
            WRITE(nim_wr,*)qprof(iline,2),qprof(iline+1,2)
          ENDIF
          ibl = item%ib
          x_left  = item%p(1,1)
          x_right = item%p(1,1)+1.
          q_left  = qprof(iline,2)
          q_right = qprof(iline+1,2)
          x_start = item%p(1,1)+0.5
          y_start = item%p(2,1)
          i_bisect = 0
          DO
            i_bisect = i_bisect+1
            pres_eq_last=pres_eq_new
            area_last=area
            p0_last%point=p0%point
            qvalue_last=qvalue

            CALL lagr_quad_eval(rbc(ibl)%rz,x_start,y_start,0_i4)
            CALL lagr_quad_eval(rbc(ibl)%pres_eq,x_start,y_start,0_i4)
            p0%point(1:2)=rbc(ibl)%rz%f
            pres_eq_new=rbc(ibl)%pres_eq%f(1)
            CALL qcompute(p0,area,qvalue,r_min,z_min,r_max,z_max)
            IF(local_debug)WRITE(nim_wr,*)area,qvalue
            IF(qvalue > q_search_value)THEN
              q_right=qvalue 
              x_right=x_start
            ELSE
              q_left=qvalue 
              x_left=x_start
            ENDIF
c
c						Check convergence
c						compute NTM threshold
            IF(( ((qvalue - q_search_value)<ind_tol0).AND.
     &           (i_bisect > 2))
     &          .OR.(i_bisect > max_bisect))THEN
             found = .TRUE.

             WRITE(nim_wr,*)"NTM information for ",mvalue," / ",nvalue
             rs_eff = SQRT(area/pi)
             rzero = 0.5_r8*(r_max+r_min)
             eps_eff = rs_eff/rzero
             shear_eff = rs_eff/qvalue*(qvalue-qvalue_last)/
     &                 (rs_eff-SQRT(area_last/pi))
c
c					Recall that my mu_e does not have
c					the trap particle fraction
             mu_e_over_nu_e = mu_e*1.46*SQRT(eps_eff)/nu_e

c
c					Evaluate poloidal magnetic field
             y_lsode(1:2)=p0%point(1:2)     ! (R,Z) coordinate
             y_lsode(3)=0.                  ! PEST angle
             y_lsode(4)=0.                  ! Area integration
             y_lsode(neq+1)=x_start         ! x coordinate in logical space
             y_lsode(neq+2)=y_start         ! y coordinate in logical space
             y_lsode(neq+3)=ibl             ! block label
             CALL get_bfield0(neq,y_lsode,b0_xyz,bigr)
             bmag = SQRT(b0_xyz(1)**2+b0_xyz(2)**2+b0_xyz(3)**2)
     &              +TINY(bmag)

             beta_prime =-(pres_eq_new-pres_eq_last)*mu0
     &                /(rs_eff-SQRT(area_last/pi))
     &                *rs_eff*(qvalue/eps_eff)**2/bmag**2
             w_nc = 6.34
     &              *beta_prime/shear_eff
     &              *mu_e_over_nu_e/(1.+mu_e_over_nu_e)
             w_d = rs_eff * 1.8 * SQRT(8.) * (k_perp/k_pll)**0.25/
     &             SQRT(ABS(eps_eff*shear_eff*REAL(nvalue,r8)))
             deltaprime =  -2.*REAL(mvalue,r8)/rs_eff
             gamma_max = (deltaprime + w_nc*0.5/w_d)/0.8227*elecd

             WRITE(nim_wr,*)'rs_eff  (m)     = ',rs_eff
             WRITE(nim_wr,*)'eps_eff         = ',eps_eff
             WRITE(nim_wr,*)'shear_eff       = ',shear_eff
             WRITE(nim_wr,*)'bmag            = ',bmag
             WRITE(nim_wr,*)'beta_prime      = ',beta_prime
             WRITE(nim_wr,*)'mu_e/nu_e       = ',mu_e_over_nu_e
             WRITE(nim_wr,*)'w_nc    (m)     = ',w_nc
             WRITE(nim_wr,*)'w_d     (m)     = ',w_d
             WRITE(nim_wr,*)'gamma_max(s^-1) = ',gamma_max
             WRITE(nim_wr,*)'tau_R local (s) = ',rs_eff**2/elecd
             WRITE(nim_wr,*)'tau_a global(s) = ',tau_a
c
c						NTM analysis
             determinant = (w_nc/deltaprime)**2 - 4.0 * w_d**2
             WRITE(nim_wr,*)'determinant     = ',determinant
             IF(determinant > 0 )THEN
c				Island widths are in meters.
               WRITE(nim_wr,*)'W_min (m) = ',
     &          0.5*(-w_nc/deltaprime-SQRT(determinant))
               WRITE(nim_wr,*)'W_max (m) = ',
     &          0.5*(-w_nc/deltaprime+SQRT(determinant))
             ELSE
               WRITE(nim_wr,*)"STABLE"
             ENDIF
             WRITE(nim_wr,*)			! blank line
             EXIT
            ENDIF
            x_start=0.5*(x_left+x_right)
          ENDDO
        ENDIF
        item => item%face(1)%p
        IF(item%id == 0 ) EXIT
      ENDDO
      ENDDO mloop
      ENDDO
      RETURN

      END SUBROUTINE threshold
