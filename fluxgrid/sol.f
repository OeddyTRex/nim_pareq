c-----------------------------------------------------------------------
c     program sol.
c     computes Soloviev's analytical equilibrium and writes it to files.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      PROGRAM sol
      USE local
      USE input_sol
      IMPLICIT NONE
      
      INTEGER(i4) :: ir,iz,ia
      REAL(r8) :: rmin,rmax,zmin,zmax,psifac,efac,psis,f0,pfac
      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: psig,rg,zg
      REAL(r8), DIMENSION(:), ALLOCATABLE :: psi,f,p,r,z
c-----------------------------------------------------------------------
c     write formats.
c-----------------------------------------------------------------------
 2010 FORMAT(4x,'mr',4x,'mz',4x,'ma',5x,'rmin',7x,'rmax',7x,'zmax'
     $     //3i6,1p,3e11.3/)
 2020 FORMAT(/4x,'iz =',i3,', z =',1p,e11.3)
 2030 FORMAT(/4x,'ir',6x,'r',9x,'psi'/)
 2040 FORMAT(i6,1p,2e11.3)
 2050 FORMAT(/4x,'ia',5x,'psi',9x,'f',10x,'p'/)
 2060 FORMAT(i6,1p,3e11.3)
c-----------------------------------------------------------------------
c     read input data.
c-----------------------------------------------------------------------
      CALL read_input
c-----------------------------------------------------------------------
c     allocate arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psi(0:ma))
      ALLOCATE(f(0:ma))
      ALLOCATE(p(0:ma))
      ALLOCATE(r(0:mr))
      ALLOCATE(z(0:mz))
      ALLOCATE(rg(0:mr,0:mz))
      ALLOCATE(zg(0:mr,0:mz))
      ALLOCATE(psig(0:mr,0:mz))
c-----------------------------------------------------------------------
c     compute constants.
c-----------------------------------------------------------------------
      f0=r0
      psis=e*f0*a*a/(2*q0*r0)
      psifac=psis/(a*r0)**2
      efac=1/(e*e)
      pfac=2*psis*(e*e+1)/(a*r0*e)**2
      rmin=r0-1.5*a
      rmax=r0+1.5*a
      zmax=1.5*e*a
      zmin=-zmax
      r=rmin+(/(ir,ir=0,mr)/)*(rmax-rmin)/mr
      z=zmin+(/(iz,iz=0,mz)/)*(zmax-zmin)/mz
c-----------------------------------------------------------------------
c     compute 2D data.
c-----------------------------------------------------------------------
      DO iz=0,mz
         DO ir=0,mr
            rg(ir,iz)=r(ir)
            zg(ir,iz)=z(iz)
            psig(ir,iz)=psis
     $           -psifac*(efac*(r(ir)*z(iz))**2+(r(ir)**2-r0**2)**2/4)
         ENDDO
      ENDDO
c-----------------------------------------------------------------------
c     compute 1D data.
c-----------------------------------------------------------------------
      psi=(/(ia,ia=1,ma+1)/)*psis/(ma+1)
      f=f0
      p=pfac*(psis-psi)
c-----------------------------------------------------------------------
c     write ascii data.
c-----------------------------------------------------------------------
      IF(out)THEN
         OPEN(UNIT=eq1_unit,FILE='sol.out')
         WRITE(eq1_unit,2010)mr,mz,ma,rmin,rmax,zmax
         DO iz=0,mz
            WRITE(eq1_unit,2020)iz,z(iz)
            WRITE(eq1_unit,2030)
            WRITE(eq1_unit,2040)(ir,r(ir),psig(ir,iz),ir=0,mr)
            WRITE(eq1_unit,2030)
         ENDDO
         WRITE(eq1_unit,2050)
         WRITE(eq1_unit,2060)(ia,psi(ia),f(ia),p(ia),ia=0,ma)
         WRITE(eq1_unit,2050)
         CLOSE(UNIT=eq1_unit)
      ENDIF
c-----------------------------------------------------------------------
c     write binary data.
c-----------------------------------------------------------------------
      CALL open_bin(eq2_unit,"sol.dat","UNKNOWN","REWIND",64_i4)
      WRITE(eq2_unit)mr,mz,ma
      WRITE(eq2_unit)rmin,rmax,zmin,zmax
      WRITE(eq2_unit)psig,psi,f,p
      CALL close_bin(eq2_unit,"sol.dat")
c-----------------------------------------------------------------------
c     draw contour plot.
c-----------------------------------------------------------------------
      IF(bin)THEN
         CALL open_bin(eq2_unit,"nina.bin","UNKNOWN","REWIND",32_i4)
         WRITE(eq2_unit)1,1
         WRITE(eq2_unit)mr,mz
         WRITE(eq2_unit)REAL(rg,4),REAL(zg,4)
         WRITE(eq2_unit)REAL(psig,4)
         CALL close_bin(eq2_unit,"nina.bin")
      ENDIF
c-----------------------------------------------------------------------
c     terminate program.
c-----------------------------------------------------------------------
      STOP
      END PROGRAM sol
