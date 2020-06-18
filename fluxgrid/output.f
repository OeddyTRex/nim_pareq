c-----------------------------------------------------------------------
c     file output.f.
c     Output files describing the equilibrium
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     1. write_out.
c     2. write_2d.
c     3. draw_eig.
c     4. write_tecplot.
c     5. write_neoclassical.
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
c     subprogram 1. write_out.
c     Writes standard ascii and binary output:
c	fluxgrid.out, fluxgrid.bin, stability.bin
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_out(gt,iua,iub)

      USE analyze
      
      TYPE(global_type), INTENT(INOUT) :: gt
      INTEGER, INTENT(IN) :: iua,iub

      INTEGER(i4) :: ipsi,ising
c      REAL(r8), PARAMETER :: mi=3.3435860e-27_r8
      REAL(r8) :: gsnorm
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: diff

c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 10   FORMAT(/3x,"mpsi",1x,"mtheta",2x,"psilow",4x,"psihigh",6x,"tol0",
     $     6x,"ipb",3x,"ipr"//2i6,1p,3e11.3,2i6/)
 20   FORMAT(/4x,"R_m (m)",3x,"Z_m (m)",4x,"R_o (m)",5x,"a (m)",5x,
     $     "A=R_o/a",//1p,5e11.3/)
 30   FORMAT(/4x,"kappa",4x,"delta_top",2x,"delta_bot"//1p,3e11.3/)
 40   FORMAT(/2x,"B_o (T)",5x,"I (MA)",4x,"I_N=I/aB",  //1p,3e11.3/)
 50   FORMAT(/6x,"q0",8x,"qmin",7x,"qmax",8x,"qa",8x,  //1p,6e11.3/)
 60   FORMAT(/5x,"li1",8x,"li2",8x,"li3",5x,"Flux_tot (w)"//1p,4e11.3/)
 70   FORMAT(/4x,"betap1",5x,"betap2",5x,"betap3",5x,"betat",6x,
     $     "betan",//1p,5e11.3/)
 100  FORMAT(/1x,"GIVEN:",9x,"S",7x,"ndens (m^-3)",//11x,1p,2e11.3/)
 110  FORMAT(/6x,"elecd (m^2/s)",2x,"Tau_A (s)",2x,"Tau_R (s)",
     &           //8x,2p,3e11.3/)
 120  FORMAT(/6x,"taua^2/ndens",2x,"taur*elecd",3x,"taua*S^1/3",
     $     //6x,1p,3e12.3/)
 150  FORMAT(/1x,"OR GIVEN:",9x,"S",7x,"Tau_A (s)", //14x,1p,2e11.3/)
 160  FORMAT(/6x,"ndens (m^-3)",3x,"elecd (m^2/s)",3x,"Tau_R (s)",
     &           //4x,3(3x,e11.3)/)
 200   FORMAT(/5x,"i",5x,"m",6x,"q",7x,"dq/dpsi",6x,"psi",8x,"rho"/)
 210   FORMAT(2i6,1p,4e11.3)
c-----------------------------------------------------------------------
c     open ascii output file and write input data.
c-----------------------------------------------------------------------
      WRITE(iua,*)
     $     "Equilibrium: ",TRIM(filename),", TYPE: ", TRIM(eq_type)
      WRITE(iua,fmt='(//1x,"INPUT PARAMETERS:")')
      WRITE(iua,10)mpsi,mtheta,psilow,psihigh,tol0,ipb,ipr
c-----------------------------------------------------------------------
c     Print info on global data.
c     NOTE: In future, need to add qcyl, pressure peaking, GA's S, etc.
c-----------------------------------------------------------------------
      WRITE(iua,fmt='(//1x,"GEOMETRIC QUANTITIES:")')
      WRITE(iua,20)gt%ro,gt%zo,gt%rmean,gt%amean,gt%aratio
      WRITE(iua,30)gt%kappa,gt%delta1,gt%delta2
      WRITE(iua,fmt='(//1x,"GLOBAL PLASMA QUANTITIES:")')
      WRITE(iua,40) gt%bt0,gt%crnt,gt%crnt/(gt%amean*gt%bt0)
      WRITE(iua,fmt='(//1x,"SAFETY FACTOR QUANTITIES:")')
      WRITE(iua,50)gt%q0,gt%qmin,gt%qmax,gt%qa
      WRITE(iua,fmt='(//1x,"CURRENT QUANTITIES:")')
      WRITE(iua,60)gt%li1,gt%li2,gt%li3,gt%psio
      WRITE(iua,fmt='(//1x,"PRESSURE QUANTITIES:")')
      WRITE(iua,70)gt%betap1,gt%betap2,gt%betap3,gt%betat,gt%betan
      WRITE(iua,fmt='(//1x,"MHD QUANTITIES:")')
      WRITE(iua,100)sfac,ndens
      WRITE(iua,110)elecd,taua,taur
      WRITE(iua,120)tauafac,taurfac,taua*sfac**(1/REAL(3,r8))
      WRITE(iua,150)sfac,1.e-7
      WRITE(iua,160)(1.e-7)**2/tauafac,taurfac/(sfac*1.e-7),sfac*1.e-7
c-----------------------------------------------------------------------
c     Print info on singular surfaces.
c-----------------------------------------------------------------------
      WRITE(iua,'(/1x,a,i2,a)')"singular surfaces for n = ",gt%nn,":"
      WRITE(iua,200)
      WRITE(iua,210)(ising,gt%msing(ising),gt%qsing(ising),
     $     gt%q1sing(ising),gt%psising(ising),SQRT(gt%psising(ising)),
     $     ising=1,gt%nsing)
      WRITE(iua,200)


c-----------------------------------------------------------------------
c     write surface quantities.
c	Write everything in gt%sq, plus rho=sqrt(psinorm) and lambdaprof
c-----------------------------------------------------------------------
      diff=ABS(delstr-gsrhs)                         ! Check GS Eq.
      gsnorm=MAXVAL(delstr)

      CALL open_bin(iub,"fluxgrid.bin","UNKNOWN","REWIND",32_i4)
      
      CALL spline_write(gt%sq,.TRUE.,.FALSE.,iua,iub)

      DO ipsi=0,gt%sq%nodes
        WRITE(iub)(/REAL(gt%sq%xs(ipsi),4),
     &              REAL(gt%sq%fs(ipsi,1:gt%sq%nqty),4),
     &              REAL(SQRT(gt%sq%xs(ipsi)),4),
     &              REAL(lambdaprof(ipsi),4),
     &              REAL(MAXVAL(diff(:,ipsi)),4),
     &              REAL(MAXVAL(diff(:,ipsi)/gsnorm),4)/)
      ENDDO

      CALL close_bin(iub,"fluxgrid.bin")
c-----------------------------------------------------------------------
c     Show 2d difference of GS Equation
c-----------------------------------------------------------------------
      CALL open_bin(iub,"gs2d.bin","UNKNOWN","REWIND",32_i4)
      WRITE(iub) INT(1,4),INT(0,4), INT(2,4)
      WRITE(iub)  INT(mtheta,4), INT(mpsi,4)
      WRITE(iub)  REAL(gt%twod%fs(1,:,:),4)
      WRITE(iub)  REAL(gt%twod%fs(2,:,:),4)
      WRITE(iub)  REAL(diff(:,:),4)
      WRITE(iub)  REAL(diff(:,:)/gsnorm,4)
      CALL close_bin(iub,"gs2d.bin")
c-----------------------------------------------------------------------
c     write stability information if asked for
c-----------------------------------------------------------------------
      IF (stability) THEN
       CALL open_bin(iub,"stability.bin","UNKNOWN","REWIND",32_i4)

       DO ipsi=0,mpsi
        WRITE(iub)(/REAL(SQRT(gt%sq%xs(ipsi)),4),
     &     REAL(dideal(ipsi),4),
     &     REAL(dres(ipsi),4),
     &     REAL(dnc(ipsi),4),
     &     REAL(dres(ipsi)/(alphas(ipsi)-hfactor(ipsi)),4),
     &     REAL(dres(ipsi)/(alphas(ipsi)-hfactor(ipsi))+dnc(ipsi),4)/)
       ENDDO

       CALL close_bin(iub,"stability.bin")
      ENDIF


c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_out


c-----------------------------------------------------------------------
c     subprogram 2. write_2d.
c     produces ascii and binary output for R,Z(tau,a).
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_2d(gt,iua,iub)

      USE analyze
      
      TYPE(global_type), INTENT(INOUT) :: gt
      INTEGER, INTENT(IN) :: iua,iub

      INTEGER(i4) :: ipsi,itau
      REAL(r8) :: tau,deta,eta,r2,psi,r,z
c-----------------------------------------------------------------------
c     WRITE formats.
c-----------------------------------------------------------------------
 2010 FORMAT(1x,'ipsi = ',i3,', psi = ',1p,e11.3)
 2015 FORMAT(1x,'ipsi = ',i3,', jpsi = ',i1,', psi = ',1p,e11.3)
 2020 FORMAT(/4x,'it',5x,'tau',9x,'r2',8x,'deta',7x,'eta',9x,'r',10x,
     $     'z'/)
 2030 FORMAT(i6,1p,6e11.3)
c-----------------------------------------------------------------------
c     open output files.
c-----------------------------------------------------------------------
      IF(.NOT.(out_2d .OR. bin_2d))RETURN
      IF(out_2d)OPEN(UNIT=iua,FILE='2d.out',STATUS='UNKNOWN')
      IF(bin_2d)CALL open_bin(iub,"2d.bin","UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     write input data along flux surfaces.
c-----------------------------------------------------------------------
      IF(out_2d .OR. bin_2d)THEN
       IF(out_2d)WRITE(iua,'(1x,a/)')"input data"
       DO ipsi=0,mpsi
          psi=gt%r2g%ys(ipsi)
          IF(out_2d)THEN
             WRITE(iua,2010)ipsi,psi
             WRITE(iua,2020)
          ENDIF
          DO itau=0,mtheta
             tau=gt%r2g%xs(itau)
             r2=gt%r2g%fs(1,itau,ipsi)
             deta=gt%r2g%fs(2,itau,ipsi)
             eta=tau+deta
             r=gt%twod%fs(1,itau,ipsi)
             z=gt%twod%fs(2,itau,ipsi)
             IF(out_2d)WRITE(iua,2030)itau,tau,r2,deta,eta,r,z
             IF(bin_2d)WRITE(iub)(/REAL(tau,4),REAL(psi,4),REAL(r2,4),
     $              REAL(deta,4),REAL(eta,4),REAL(r,4),REAL(z,4)/)
          ENDDO
          IF(out_2d)WRITE(iua,2020)
          IF(bin_2d)WRITE(iub)
       ENDDO
       IF(extr > 1. .AND. mvac > 0) THEN
       DO ipsi=1,mvac
             psi=1. + 0.5*(1.-extr)*(gt%rs2-gt%rs1)*
     $                    (REAL(ipsi)/REAL(mvac))**2
          DO itau=0,mtheta
             r=gt%vac%fs(1,itau,ipsi)
             z=gt%vac%fs(2,itau,ipsi)
             tau=gt%vac%xs(itau)
             deta=gt%r2g%fs(2,itau,mpsi)
             eta=tau+deta
             r2=r**2 + z**2
             IF(bin_2d)WRITE(iub)(/REAL(tau,4),REAL(psi,4),REAL(r2,4),
     $              REAL(deta,4),REAL(eta,4),REAL(r,4),REAL(z,4)/)
          ENDDO
          IF(bin_2d)WRITE(iub)
       ENDDO
       ENDIF
      ENDIF
c-----------------------------------------------------------------------
c     CLOSE output file.
c-----------------------------------------------------------------------
      IF(bin_2d)CALL close_bin(iub,"2d.bin")
c-----------------------------------------------------------------------
c     terminate routine.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_2d



c-----------------------------------------------------------------------
c     subprogram 4. draw_eig.
c     produces 2D binary output for xdraw.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE draw_eig(gt,eig_unit)

      USE analyze

      INTEGER,INTENT(IN) :: eig_unit
      TYPE(global_type), INTENT(INOUT) :: gt

      INTEGER(i4) :: ipsi,itheta
c-----------------------------------------------------------------------
c     open file and fill pointers.
c-----------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(gt%eigvec))RETURN
      CALL open_bin(eig_unit,"eig.bin","UNKNOWN","REWIND",32_i4)
c-----------------------------------------------------------------------
c     write data on axis.
c-----------------------------------------------------------------------
      DO itheta=0,mtheta
         WRITE(eig_unit)
     $        0._4,
     $        REAL(gt%r2g%xs(itheta),4),
     $        REAL(gt%ro,4),REAL(gt%zo,4),
     $        REAL(gt%eigveco,4)
      ENDDO
      WRITE(eig_unit)
c-----------------------------------------------------------------------
c     write data off axis.
c-----------------------------------------------------------------------
      DO ipsi=0,mpsi
         DO itheta=0,mtheta
            WRITE(eig_unit)
     $           REAL(SQRT(gt%r2g%ys(ipsi)),4),
     $           REAL(gt%r2g%xs(itheta),4),
     $           REAL(gt%twod%fs(1:2,itheta,ipsi),4),
     $           REAL(gt%eigvec(ipsi,itheta,:),4)
         ENDDO
         WRITE(eig_unit)
      ENDDO
c-----------------------------------------------------------------------
c     close file.
c-----------------------------------------------------------------------
      CALL close_bin(eig_unit,"eig.bin")
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE draw_eig




c-----------------------------------------------------------------------
c     subprogram 7. write_tecplot.
c     writes output file for Tecplot plotting program
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_tecplot(gt)

      USE analyze

      TYPE(global_type), INTENT(INOUT) :: gt

      INTEGER(i4) :: itheta,ipsi,yspace, nimx, mnode
      REAL(r8) :: tau,psi, bmod,gss,f,fprime,jac,ldelstr
      REAL(r8) :: r,z,br,bz,bt,jr,jz,jt, p,conc
c-----------------------------------------------------------------------
c     Two dimensional data
c-----------------------------------------------------------------------
      nimx=mpsi+mvac
      OPEN(UNIT=tec2d,FILE='2d.dat',STATUS='UNKNOWN')
      WRITE(tec2d,*) 'VARIABLES = "R", "Z", "psi", "tau", "bmod"'
      WRITE(tec2d,*) 'ZONE ',',i=',mtheta+1,' j=',nimx+1,',F=POINT'
      DO ipsi=0,mpsi
         psi=gt%r2g%ys(ipsi)
         DO itheta=0,mtheta
            tau=gt%r2g%xs(itheta)
            r=gt%twod%fs(1,itheta,ipsi)
            z=gt%twod%fs(2,itheta,ipsi)
            f=gt%sq%fs(ipsi,1)
            gss=gt%twod%fs(4,itheta,ipsi)
            bmod=(f**2 + gss)**0.5/r
            WRITE(tec2d,*) r, z, psi, tau, bmod
         ENDDO
      ENDDO
      IF(extr > 1. .AND. mvac > 0) THEN
      DO ipsi=1,mvac
         DO itheta=0,mtheta
           psi=gt%vac%fs(11,itheta,ipsi)
           tau=gt%r2g%xs(itheta)
           r=gt%vac%fs(1,itheta,ipsi)
           z=gt%vac%fs(2,itheta,ipsi)
           br=gt%vac%fs(3,itheta,ipsi)
           bz=gt%vac%fs(4,itheta,ipsi)
           bt=gt%vac%fs(5,itheta,ipsi)
           bmod=(br**2+bz**2+bt**2)**0.5
           WRITE(tec2d,*) r, z, psi, tau, bmod
         ENDDO
      ENDDO
      ENDIF
      CLOSE(UNIT=tec2d)
c-----------------------------------------------------------------------
c     Stability data
c-----------------------------------------------------------------------
      IF (stability) THEN
       OPEN(UNIT=tec1d,FILE='stability.dat',STATUS='UNKNOWN')
       WRITE(tec1d,*) 
     &  "VARIABLES = `r,F',p',q',D_I,D_R, D_n_c,D_t_o_t,f_u_t,Hegna"
       WRITE(tec1d,*) "ZONE i=",mpsi+1, ", F=POINT"

       DO ipsi=0,mpsi
         WRITE(tec1d,*) SQRT(gt%sq%xs(ipsi)), 
     &            gt%sq%fs1(ipsi,1:3)/gt%psio,
     &            dideal(ipsi), dres(ipsi), dnc(ipsi),
     &            dres(ipsi)/(alphas(ipsi)-hfactor(ipsi))+dnc(ipsi),
     &            fcirc(ipsi),alphas(ipsi)-hfactor(ipsi)
       ENDDO

       CLOSE(UNIT=tec1d)
      ENDIF
c-----------------------------------------------------------------------
c     compute and write coordinates and equilibrium magnetic field.
c     This is what fluxgrid.dat contains so good for seeing what nimrod
c	is getting as input
c-----------------------------------------------------------------------
      mnode=mpsi
      IF(extr > 1. .AND. mvac > 0) mnode=mpsi+mvac+1
      OPEN(UNIT=tec2d,FILE="field.dat",STATUS="UNKNOWN")
      WRITE(tec2d,*) 'VARIABLES=R,Z,`y,`q,B_R,B_Z,B_t,J_R,J_Z,J_t,p,con'
      WRITE(tec2d,*) 'ZONE ',',i=',mtheta+1,' j=',mnode,',F=POINT'
      DO ipsi=0,mpsi
       psi=gt%r2g%ys(ipsi)
       f = gt%sq%fs(ipsi,1)					! R B_tor
       fprime = gt%sq%fs1(ipsi,1)/gt%psio			! dF/dpsi
       p = gt%sq%fs(ipsi,2)					! R B_tor
c       pprime = gt%sq%fs1(ipsi,2)/gt%psio			! mu0 pprime
       DO itheta=0,mtheta
           tau=gt%r2g%xs(itheta)
           ldelstr = delstr(itheta,ipsi) 
           IF (j_t == "use_gs") ldelstr= gsrhs(itheta,ipsi)
           r = gt%twod%fs(1,itheta,ipsi)
           z = gt%twod%fs(2,itheta,ipsi)
           jac = gt%twod%fs(3,itheta,ipsi)
           br = gt%twod%fsx(1,itheta,ipsi)/jac
           bz = gt%twod%fsx(2,itheta,ipsi)/jac
           bt = f/r

           jr = -fprime*br/mu0
           jz = -fprime*bz/mu0
           jt = ldelstr/(mu0*r**2)		        !Contravariant

         WRITE(tec2d,*) r,z,psi,tau,br,bz,bt,jr,jz,jt,p,0.
       ENDDO
      ENDDO
      IF(extr > 1. .AND. mvac > 0) THEN
        DO ipsi=1,mvac
         DO itheta=0,mtheta
           psi=gt%vac%fs(11,itheta,ipsi)
           tau=gt%vac%xs(itheta)
           r=gt%vac%fs(1,itheta,ipsi)
           z=gt%vac%fs(2,itheta,ipsi)
           br=gt%vac%fs(3,itheta,ipsi)
           bz=gt%vac%fs(4,itheta,ipsi)
           bt=gt%vac%fs(5,itheta,ipsi)
           jr=gt%vac%fs(6,itheta,ipsi)
           jz=gt%vac%fs(7,itheta,ipsi)
           jt=gt%vac%fs(8,itheta,ipsi)
           p =gt%vac%fs(9,itheta,ipsi)
           conc =gt%vac%fs(10,itheta,ipsi)
         WRITE(tec2d,*) r,z,psi,tau,br,bz,bt,jr,jz,jt,p,conc
         ENDDO
        ENDDO
      ENDIF
      CLOSE(UNIT=tec2d)
c-----------------------------------------------------------------------
c     1D File - Surface quantities 
c-----------------------------------------------------------------------
      OPEN(UNIT=tec1d,FILE='1d.dat',STATUS='UNKNOWN')
       WRITE(tec1d,*)
     &      "VARIABLES=",gt%sq%title(0:gt%sq%nqty),",`r  <J.B/B^2>"
       WRITE(tec1d,*) "ZONE i=",gt%sq%nodes+1, ", F=POINT"

       DO ipsi=0,gt%sq%nodes
         WRITE(tec1d,*) gt%sq%xs(ipsi), gt%sq%fs(ipsi,1:gt%sq%nqty),
     &           SQRT(gt%sq%xs(ipsi)),  lambdaprof(ipsi)
       ENDDO

c-----------------------------------------------------------------------
c     1D File - Text Data (complicated because I want it to look nice)
c-----------------------------------------------------------------------

1700   FORMAT('TEXT X=35 Y=',i3,1x,
     &       'F=HELV-BOLD HU=POINT H=14 T="',a,'"')
1740   FORMAT('TEXT X=8 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',a4,'"')
1750   FORMAT('TEXT X=20 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="=','"')
1760   FORMAT('TEXT X=23 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',f10.2,'"')
1770   FORMAT('TEXT X=55 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',a4,'"')
1780   FORMAT('TEXT X=67 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="=','"')
1790   FORMAT('TEXT X=70 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',f10.2,'"')
1791   FORMAT('TEXT X=70 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',e9.2,'"')
1771   FORMAT('TEXT X=55 Y=',i3,1x, 
     &       'F=HELV-BOLD HU=POINT H=14 T="',a6,'"')

       yspace=5			!Position text data in frame

       WRITE (tec1d,1700) yspace+85,'Global Data'

       WRITE (tec1d,1740) yspace+75,'B_o'
       WRITE (tec1d,1750) yspace+75
       WRITE (tec1d,1760) yspace+75,gt%bt0

       WRITE (tec1d,1740) yspace+70,'I_P'
       WRITE (tec1d,1750) yspace+70
       WRITE (tec1d,1760) yspace+70,gt%crnt

       WRITE (tec1d,1740) yspace+60,'`b  = '
       WRITE (tec1d,1750) yspace+60
       WRITE (tec1d,1760) yspace+60,100*gt%betat

       WRITE (tec1d,1740) yspace+55,'`b_P'
       WRITE (tec1d,1750) yspace+55
       WRITE (tec1d,1760) yspace+55,gt%betap1

       WRITE (tec1d,1740) yspace+50,'`b_N'
       WRITE (tec1d,1750) yspace+50
       WRITE (tec1d,1760) yspace+50,gt%betan

       WRITE (tec1d,1740) yspace+45,'I_N'
       WRITE (tec1d,1750) yspace+45
       WRITE (tec1d,1760) yspace+45,gt%crnt/(gt%amean*gt%bt0)

       WRITE (tec1d,1740) yspace+35,'q_o'
       WRITE (tec1d,1750) yspace+35
       WRITE (tec1d,1760) yspace+35,gt%q0

       WRITE (tec1d,1740) yspace+30,'q_m'
       WRITE (tec1d,1750) yspace+30
       WRITE (tec1d,1760) yspace+30,gt%qmin

       WRITE (tec1d,1740) yspace+25,'q_a'
       WRITE (tec1d,1750) yspace+25
       WRITE (tec1d,1760) yspace+25,gt%qa

       WRITE (tec1d,1740) yspace+20,'l_i'
       WRITE (tec1d,1750) yspace+20
       WRITE (tec1d,1760) yspace+20,gt%li1

       WRITE (tec1d,1770) yspace+75,'Z_o'
       WRITE (tec1d,1780) yspace+75
       WRITE (tec1d,1790) yspace+75, gt%zo

       WRITE (tec1d,1770) yspace+70,'R_o'
       WRITE (tec1d,1780) yspace+70
       WRITE (tec1d,1790) yspace+70, gt%ro

       WRITE (tec1d,1770) yspace+65,'<a>'
       WRITE (tec1d,1780) yspace+65
       WRITE (tec1d,1790) yspace+65,gt%amean

       WRITE (tec1d,1770) yspace+60,'`k'
       WRITE (tec1d,1780) yspace+60
       WRITE (tec1d,1790) yspace+60,gt%kappa

       WRITE (tec1d,1770) yspace+55,'`d'
       WRITE (tec1d,1780) yspace+55
       WRITE (tec1d,1790) yspace+55,gt%delta1

       WRITE (tec1d,1770) yspace+50,'A'
       WRITE (tec1d,1780) yspace+50
       WRITE (tec1d,1790) yspace+50,gt%aratio 

       WRITE (tec1d,1770) yspace+35,'`t_A'
       WRITE (tec1d,1780) yspace+35
       WRITE (tec1d,1790) yspace+35,taua

       WRITE (tec1d,1770) yspace+30,'`t_R'
       WRITE (tec1d,1780) yspace+30
       WRITE (tec1d,1790) yspace+30,taur

       WRITE (tec1d,1770) yspace+25,'S'
       WRITE (tec1d,1780) yspace+25
       WRITE (tec1d,1790) yspace+25,sfac

       WRITE (tec1d,1770) yspace+20,'n'
       WRITE (tec1d,1780) yspace+20
       WRITE (tec1d,1791) yspace+20,ndens

       WRITE (tec1d,1771) yspace+15,'`h/`m_0'
       WRITE (tec1d,1780) yspace+15
       WRITE (tec1d,1790) yspace+15,elecd

      CLOSE(UNIT=tec1d)
c-----------------------------------------------------------------------
c     And routines just aren't finished unless we have the fucking:
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_tecplot


c-----------------------------------------------------------------------
c     subprogram 8. write_neoclassical
c     Analytic calculation of deltaprime and evaluation of pressure
c	contributions to the Rutherford equation at rational surfaces.
c     REFERENCES:
c         C.C. Hegna and J.D. Callen,Phys. Plasmas 1 (1994) 2308.
c         C.C. Hegna, Phys. Plasmas 6 (1999)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE write_neoclassical(gt,iua)

      USE analyze

      USE analyze
      TYPE(global_type), INTENT(INOUT) :: gt
      INTEGER, INTENT(IN) :: iua

      INTEGER(i4) :: itheta,ipsi, num_dp,idp,j
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: mdp,ndp
      REAL(r8) :: f, fprime, psio, gss, ldelstr, r
      REAL(r8) :: density, r_min,r_max
      REAL(r8), DIMENSION(0:mtheta,0:mpsi) :: jdotb, bdotgf,work,bdotb
      LOGICAL :: file_stat
      TYPE(spline_type) :: profiles

      TYPE :: node_type
        REAL(r8) :: psirs,q,qprime,sigmaprime,f,gss,gtt
        REAL(r8) :: lambda,deltaprime,rho,eps,pprime
        REAL(r8) :: w_sat,w_nc,w_ps,w_star1,w_star2
        INTEGER(i4) :: m,n
        TYPE(node_type), POINTER :: next    ! ptr to next node in queue
      END TYPE node_type
      TYPE(node_type), POINTER :: front,rear,rs

c-----------------------------------------------------------------------
c     Check for the right input parameters
c-----------------------------------------------------------------------
      IF(.NOT.((ipb == 0) .AND. (ipr == 2)).AND.
     &  .NOT.(angle_method == 'jac'))
     &  CALL nim_stop("write_neoclassical:  Requires PEST coordinates.")
c-----------------------------------------------------------------------
c     open input file which contains mode list or create
c-----------------------------------------------------------------------
      INQUIRE(FILE="mode.list",EXIST=file_stat)
        IF (.NOT.file_stat)THEN
          OPEN(UNIT=iua,FILE="mode.list",STATUS="NEW")
          num_dp=3
          ALLOCATE(mdp(num_dp),ndp(num_dp))
          mdp(1)=2;  ndp(1)=1
          mdp(2)=3;  ndp(2)=2
          mdp(3)=5;  ndp(3)=4
          WRITE(iua,*)num_dp
          DO idp=1,num_dp
            WRITE(UNIT=iua,fmt='(x,i4,x,i4)') mdp(idp), ndp(idp)
          ENDDO
          CLOSE(UNIT=iua)
        ELSE
          OPEN(UNIT=iua,FILE="mode.list",STATUS="OLD")
          READ(UNIT=iua,FMT=*)num_dp
          ALLOCATE(mdp(num_dp),ndp(num_dp))
          DO idp=1,num_dp
            READ(UNIT=iua,FMT='(x,i4,x,i4)') mdp(idp), ndp(idp)
          ENDDO
          CLOSE(UNIT=iua)
        ENDIF

c-----------------------------------------------------------------------
c     Build linked-list of rational surfaces
c-----------------------------------------------------------------------
      CALL spline_alloc(qprof,mpsi,1_i4)
      qprof%xs(:)=gt%sq%xs(:)
      qprof%fs(:,1)=gt%sq%fs(:,3)
      CALL spline_fit(qprof,"extrap")
      NULLIFY(front,rear)
      DO idp=1,num_dp
          q_search_value=REAL(mdp(idp))/REAL(ndp(idp))
          CALL q_search(qprof%xs(0),qprof%xs(mpsi))
          DO j=1,num_root
             IF (.NOT. ASSOCIATED (front) ) THEN
               ALLOCATE(front)
               rear => front
             ELSE
               ALLOCATE(rear%next)
               rear => rear%next
             ENDIF
             rear%psirs = psi_root(j)
             rear%m = mdp(idp)
             rear%n = ndp(idp)
             rear%q = q_search_value
             CALL spline_eval(gt%sq,rear%psirs,1_i4)
             rear%f = gt%sq%f(1)
             rear%pprime = gt%sq%f1(2)
             rear%qprime = gt%sq%f1(3)
             NULLIFY (rear%next)
          ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     Calculate <J.B>/<B.Grad phi> and an eps for each surface
c-----------------------------------------------------------------------
      CALL spline_alloc(profiles,mpsi,8_i4)
      psio=gt%psio						! total flux
      DO ipsi=0,mpsi
          f = gt%sq%fs(ipsi,1)					! R B_tor
          fprime = gt%sq%fs1(ipsi,1)/psio			! dF/dpsi
          r_min=HUGE(0._r8)
          r_max=TINY(0._r8)
          DO itheta=0,mtheta
            r = gt%twod%fs(1,itheta,ipsi)
            IF(r < r_min) r_min = r
            IF(r > r_max) r_max = r
            gss = gt%twod%fs(4,itheta,ipsi)
            ldelstr= delstr(itheta,ipsi)
            IF (j_t == "use_gs") ldelstr= gsrhs(itheta,ipsi)
            jdotb(itheta,ipsi) = (ldelstr*f - fprime*gss)/r**2
            bdotgf(itheta,ipsi) = f/r**2
            bdotb(itheta,ipsi) = (f**2+ gss)/r**2
          ENDDO
          profiles%fs(ipsi,5)=(r_max/r_min-r_min/r_max)/4.	! eps
      ENDDO
c-----------------------------------------------------------------------
c     Calculation is not based on flux averages, but just chi averages.
c-----------------------------------------------------------------------
c     DO ipsi=0,mpsi
c       DO itheta=0,mtheta
c         jac = 1./gt%twod%fs(3,itheta,ipsi)
c         jdotb(itheta,ipsi) = jdotb(itheta,ipsi) * jac
c         bdotgf(itheta,ipsi) = bdotgf(itheta,ipsi) * jac
c       ENDDO
c     ENDDO
c-----------------------------------------------------------------------
c     Compute flux averages and fit to splines
c-----------------------------------------------------------------------
      profiles%fs(:,7)=dnc(:)
      profiles%fs(:,8)=dres(:)/(alphas(:)-hfactor(:))
      profiles%xs(:) = gt%sq%xs(:)
      work(:,:)=jdotb(:,:)/gt%twod%fs(3,:,:)/bdotb(:,:)
      CALL fluxav(gt,work,profiles%fs(:,1))
      work(:,:)=bdotgf(:,:)/gt%twod%fs(3,:,:)
      CALL fluxav(gt,work,profiles%fs(:,2))
      work(:,:)=gt%twod%fs(4,:,:)/gt%twod%fs(3,:,:)
      CALL fluxav(gt,work,profiles%fs(:,3))   			! gss
      work(:,:)=gt%twod%fs(6,:,:)/gt%twod%fs(3,:,:)
      CALL fluxav(gt,work,profiles%fs(:,4))  		 	! gtt
      work(:,:)=gt%twod%fs(5,:,:)/gt%twod%fs(3,:,:)
      CALL fluxav(gt,work,profiles%fs(:,6)) 		  	! gst
      CALL spline_fit(profiles,"extrap")
      CALL spline_int(profiles)
      profiles%fs(:,2)=profiles%fsi(:,2)/profiles%fsi(mpsi,2)
      CALL spline_fit(profiles,"extrap")
c-----------------------------------------------------------------------
c     Evaluate stability parameters at each rational surface
c     rs = linked list of rational surfaces
c-----------------------------------------------------------------------
      rs => front
      DO WHILE(ASSOCIATED(rs))
        CALL spline_eval(profiles,rs%psirs,1_i4)
        rs%sigmaprime=profiles%f1(1)
        rs%rho=profiles%f(2)
        rs%gss=profiles%f(3)
        rs%gtt=profiles%f(4)
        rs%eps=profiles%f(5)

        rs%lambda=rs%f * rs%q * rs%sigmaprime
     &                 /(2.*REAL(rs%n)*rs%qprime)/(rs%gtt*rs%gss)**0.5 
        rs%lambda=ABS(rs%lambda)
        IF(rs%lambda > 0.9999 )  rs%lambda = 0.9999

        rs%deltaprime=-2.*pi  * REAL(rs%n) * rs%gtt**0.5 * rs%lambda 
     &                /TAN(pi*rs%lambda)

        IF(rs%lambda > 0.5)THEN
           rs%w_sat=2.04*(rs%lambda-0.5) / (REAL(rs%m)*rs%gtt**0.5)
        ELSE
           rs%w_sat=0.
        ENDIF

        rs%w_nc=profiles%f(7)
        rs%w_ps=profiles%f(8)
c       rs%w_nc=4.3 *rs%eps**0.5 *rs%q *rs%pprime *gt%ro**2
c    &           /(rs%qprime*rs%gss)

c                                       Evaluate w_star factor in:
c                                        Collisional limit
        density=1.0e+19
        rs%w_star1=3.275e+13
     &      *(rs%q*rs%pprime/(rs%qprime*mu0))**2
     &      /(rs%gss/gt%ro**2)**2
     &               /density
c
c                                        Collisionless limit
        rs%w_star2=rs%w_star1*rs%eps**1.5


        rs  => rs%next			! Go to next rational surface
      ENDDO
c-----------------------------------------------------------------------
c     open output file.
c-----------------------------------------------------------------------
      OPEN(UNIT=iua,FILE="dp.out",STATUS="UNKNOWN")
      WRITE(iua,*)"Conventional Tearing"
      WRITE(iua,*)" m  n        psi           r         lambda     
     &   deltap    deltap*gxx^0.5    w_sat"

      rs => front
      DO WHILE(ASSOCIATED(rs))
        WRITE(UNIT=iua,fmt='(x,i2,x,i2,7(2x,e12.5))')
     &                    rs%m
     &                    ,rs%n
     &                    ,rs%psirs
     &                    ,rs%rho
     &                    ,rs%lambda
     &                    ,rs%deltaprime
     &                    ,rs%deltaprime*rs%gtt**0.5
     &                    ,rs%w_sat
c    &                    ,rs%qprime/gt%psio
c    &                    ,rs%gss
c    &                    ,rs%pprime/mu0/gt%psio
c    &                    ,rs%gtt
c    &                    ,rs%sigmaprime
c    &                    ,rs%eps
        rs  => rs%next
      ENDDO

      WRITE(iua,*)"Neoclassical Tearing"
      WRITE(UNIT=iua,fmt='(2x,"m",2x,"n",3x,"w_nc-pressure",3x,
     &   "w_ps-pressure",7x,"w_ps+w_nc")')
      rs => front
      DO WHILE(ASSOCIATED(rs))
        WRITE(UNIT=iua,fmt='(x,i2,x,i2,7(4x,e12.5))')
     &                    rs%m
     &                    ,rs%n
     &                    ,rs%w_nc
     &                    ,rs%w_ps
     &                    ,rs%w_ps+rs%w_nc
        rs  => rs%next
      ENDDO

      WRITE(iua,*)"Polarization Tearing with n =",density
      WRITE(UNIT=iua,fmt='(2x,"m",2x,"n",5x,"w-star1",7x,"w-star2")')
      rs => front
      DO WHILE(ASSOCIATED(rs))
        WRITE(UNIT=iua,fmt='(x,i2,x,i2,7(2x,e12.5))')
     &                    rs%m
     &                    ,rs%n
     &                    ,rs%w_star1
     &                    ,rs%w_star2
        rs  => rs%next
      ENDDO
      CLOSE(UNIT=iua)


c-----------------------------------------------------------------------
c     Create secondary input file for the tear code.
c-----------------------------------------------------------------------
      OPEN(UNIT=iua,FILE="tear_input.dat",STATUS="UNKNOWN")
      WRITE(iua,*)mpsi,gt%psio
      DO ipsi=0,mpsi
         WRITE(UNIT=iua,fmt='(9(x,e12.5))')
     $      profiles%xs(ipsi)
     $     ,profiles%fs(ipsi,1)
     $     ,profiles%fs(ipsi,3)
     $     ,profiles%fs(ipsi,4)
     $     ,profiles%fs(ipsi,6)
     $     ,gt%sq%fs(ipsi,1)
     $     ,gt%sq%fs(ipsi,2)
     $     ,gt%sq%fs(ipsi,3)
      ENDDO
      CLOSE(UNIT=iua)

c-----------------------------------------------------------------------
c     terminate routine
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE write_neoclassical
