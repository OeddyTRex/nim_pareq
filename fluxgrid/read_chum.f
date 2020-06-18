c-----------------------------------------------------------------------
c     subprogram 1. read_chum.
c     reads data from Ming Chu's equilibrium and eigenfunction file.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_chum(gt)
      USE local
      USE inverse
      USE input
      USE global
      USE physdat
      IMPLICIT NONE

      TYPE(global_type), INTENT(OUT) :: gt

      LOGICAL, PARAMETER :: transform_flag=.TRUE.
      INTEGER(i4) :: itau,ia,jpsi,itht,ntor,iqty,m
      REAL(4) :: rquot,omega0,growth
      REAL(4), DIMENSION(:), POINTER :: psi,f,p,q
      REAL(4), DIMENSION(:,:), POINTER :: rcc,zcc,delrr,delri,
     $     delzr,delzi,xphir,xphii,bradr,bradi,baxlr,baxli,btorr,btori,
     $     xpsir,xpsii,xchir,xchii,bnrmr,bnrmi,bpolr,bpoli,dprer,dprei
      REAL(r8), DIMENSION(:), POINTER :: eigcono
      REAL(r8), DIMENSION(:,:,:), POINTER :: eigcon
      REAL(r8), DIMENSION(4) :: x,ff1,ff2
      REAL(r8) :: x0,rho,tau
c-----------------------------------------------------------------------
c     open file and read scalars.
c-----------------------------------------------------------------------
      CALL open_bin(eq1_unit,TRIM(filename),"UNKNOWN","REWIND",32_i4)
      READ(eq1_unit)jpsi,itht,ntor
      READ(eq1_unit)rquot,omega0,growth
c-----------------------------------------------------------------------
c     allocate 1D arrays.
c-----------------------------------------------------------------------
      ALLOCATE(psi(jpsi))
      ALLOCATE(f(jpsi))
      ALLOCATE(p(jpsi))
      ALLOCATE(q(jpsi))
c-----------------------------------------------------------------------
c     read 1D arrays.
c-----------------------------------------------------------------------
      READ(eq1_unit)psi
      READ(eq1_unit)f
      READ(eq1_unit)p
      READ(eq1_unit)q
c-----------------------------------------------------------------------
c     allocate 2D arrays.
c-----------------------------------------------------------------------
      ALLOCATE(rcc(itht,jpsi))
      ALLOCATE(zcc(itht,jpsi))
      ALLOCATE(dprer(itht,jpsi))
      ALLOCATE(dprei(itht,jpsi))
      ALLOCATE(delrr(itht,jpsi))
      ALLOCATE(delri(itht,jpsi))
      ALLOCATE(delzr(itht,jpsi))
      ALLOCATE(delzi(itht,jpsi))
      ALLOCATE(xphir(itht,jpsi))
      ALLOCATE(xphii(itht,jpsi))
      ALLOCATE(bradr(itht,jpsi))
      ALLOCATE(bradi(itht,jpsi))
      ALLOCATE(baxlr(itht,jpsi))
      ALLOCATE(baxli(itht,jpsi))
      ALLOCATE(btorr(itht,jpsi))
      ALLOCATE(btori(itht,jpsi))
      ALLOCATE(xpsir(itht,jpsi))
      ALLOCATE(xpsii(itht,jpsi))
      ALLOCATE(xchir(itht,jpsi))
      ALLOCATE(xchii(itht,jpsi))
      ALLOCATE(bnrmr(itht,jpsi))
      ALLOCATE(bnrmi(itht,jpsi))
      ALLOCATE(bpolr(itht,jpsi))
      ALLOCATE(bpoli(itht,jpsi))
c-----------------------------------------------------------------------
c     read 2D arrays and close file
c-----------------------------------------------------------------------
      READ(eq1_unit)rcc
      READ(eq1_unit)zcc
      READ(eq1_unit)dprer
      READ(eq1_unit)dprei
      READ(eq1_unit)delrr
      READ(eq1_unit)delri
      READ(eq1_unit)delzr
      READ(eq1_unit)delzi
      READ(eq1_unit)xphir
      READ(eq1_unit)xphii
      READ(eq1_unit)bradr
      READ(eq1_unit)bradi
      READ(eq1_unit)baxlr
      READ(eq1_unit)baxli
      READ(eq1_unit)btorr
      READ(eq1_unit)btori
      READ(eq1_unit)xpsir
      READ(eq1_unit)xpsii
      READ(eq1_unit)xchir
      READ(eq1_unit)xchii
      READ(eq1_unit)xphir
      READ(eq1_unit)xphii
      READ(eq1_unit)bnrmr
      READ(eq1_unit)bnrmi
      READ(eq1_unit)bpolr
      READ(eq1_unit)bpoli
      READ(eq1_unit)btorr
      READ(eq1_unit)btori
      CALL close_bin(eq1_unit,filename)
c-----------------------------------------------------------------------
c     copy and modify 1D arrays.
c-----------------------------------------------------------------------
      ii%psio=psi(jpsi)
      ii%ma=jpsi-1
      ii%mtau=itht
      CALL spline_alloc(ii%sq_in,ii%ma,4_i4)
      ii%sq_in%xs(:)=psi/ii%psio
      ii%sq_in%fs(:,1)=f
      ii%sq_in%fs(:,2)=mu0*p
      ii%sq_in%fs(:,3)=q
      ii%sq_in%fs(:,4)=0
c-----------------------------------------------------------------------
c     extrapolate coordinates to magnetic axis.
c-----------------------------------------------------------------------
      m=4
      x=psi(1:4)
      x0=0
      DO ia=1,4
         ff1(ia)=SUM(rcc(1:ii%mtau,ia))/ii%mtau
         ff2(ia)=SUM(zcc(1:ii%mtau,ia))/ii%mtau
      ENDDO
      CALL extrap(m,x,ff1,x0,ii%ro)
      CALL extrap(m,x,ff2,x0,ii%zo)
c-----------------------------------------------------------------------
c     copy 2D equilibrium arrays.
c-----------------------------------------------------------------------
      ALLOCATE(ii%rg(0:ii%mtau,0:ii%ma))
      ALLOCATE(ii%zg(0:ii%mtau,0:ii%ma))
      ii%rg(0:ii%mtau-1,:)=rcc
      ii%zg(0:ii%mtau-1,:)=zcc
      ii%rg(ii%mtau,:)=ii%rg(0,:)
      ii%zg(ii%mtau,:)=ii%zg(0,:)
c-----------------------------------------------------------------------
c     set flags.
c-----------------------------------------------------------------------
      ii%p1flag=.FALSE.
      ii%symflag=.FALSE.
c-----------------------------------------------------------------------
c     set up perturbed eigenvectors.
c-----------------------------------------------------------------------
      ALLOCATE(ii%eigvec)
      gt%ntor=ntor
      gt%rquot=rquot
      gt%omega0=omega0
      gt%growth=growth
      CALL bicube_alloc(ii%eigvec,ii%mtau,ii%ma,14_i4)
      ii%eigvec%fs(1,0:ii%mtau-1,:)=dprer
      ii%eigvec%fs(2,0:ii%mtau-1,:)=dprei
      ii%eigvec%fs(3,0:ii%mtau-1,:)=delrr
      ii%eigvec%fs(4,0:ii%mtau-1,:)=delzr
      ii%eigvec%fs(5,0:ii%mtau-1,:)=xphir
      ii%eigvec%fs(6,0:ii%mtau-1,:)=delri
      ii%eigvec%fs(7,0:ii%mtau-1,:)=delzi
      ii%eigvec%fs(8,0:ii%mtau-1,:)=xphii
      ii%eigvec%fs(9,0:ii%mtau-1,:)=bradr
      ii%eigvec%fs(10,0:ii%mtau-1,:)=baxlr
      ii%eigvec%fs(11,0:ii%mtau-1,:)=btorr
      ii%eigvec%fs(12,0:ii%mtau-1,:)=bradi
      ii%eigvec%fs(13,0:ii%mtau-1,:)=baxli
      ii%eigvec%fs(14,0:ii%mtau-1,:)=btori
      ii%eigvec%fs(:,ii%mtau,:)=ii%eigvec%fs(:,0,:)
c-----------------------------------------------------------------------
c     transform from contravariant to cylindrical projections.
c-----------------------------------------------------------------------
      IF(transform_flag)THEN
         ALLOCATE(eigcon(0:ii%mtau,0:ii%ma,8))
         eigcon(0:ii%mtau-1,:,1)=xpsir
         eigcon(0:ii%mtau-1,:,2)=xchir
         eigcon(0:ii%mtau-1,:,3)=xpsii
         eigcon(0:ii%mtau-1,:,4)=xchii
         eigcon(0:ii%mtau-1,:,5)=bnrmr
         eigcon(0:ii%mtau-1,:,6)=bpolr
         eigcon(0:ii%mtau-1,:,7)=bnrmi
         eigcon(0:ii%mtau-1,:,8)=bpoli
         eigcon(ii%mtau,:,:)=eigcon(0,:,:)
         CALL transform(eigcon)
c-----------------------------------------------------------------------
c     extrapolate eigcon to magnetic axis.
c-----------------------------------------------------------------------
         ALLOCATE(eigcono(SIZE(eigcon,3)))
         DO iqty=1,SIZE(eigcon,3)
            DO ia=1,4
               ff1(ia)=SUM(eigcon(1:ii%mtau,ia-1,iqty))/ii%mtau
            ENDDO
            CALL extrap(m,x,ff1,x0,eigcono(iqty))
         ENDDO
c-----------------------------------------------------------------------
c     diagnose contravariant components.
c-----------------------------------------------------------------------
         CALL open_bin(eq1_unit,"eigcon.bin","UNKNOWN","REWIND",32_i4)
         rho=0
         DO itau=0,ii%mtau
            tau=twopi*itau/ii%mtau
            WRITE(eq1_unit)REAL(rho,4),REAL(tau,4),REAL(ii%ro,4),
     $           REAL(ii%zo,4),REAL(eigcono(:),4)
         ENDDO
         WRITE(eq1_unit)
         DO ia=0,ii%ma
            rho=SQRT(ii%sq_in%xs(ia))
            DO itau=0,ii%mtau
               tau=twopi*itau/ii%mtau
               WRITE(eq1_unit)
     $              REAL(rho,4),REAL(tau,4),
     $              REAL(ii%rg(itau,ia),4),
     $              REAL(ii%zg(itau,ia),4),
     $              REAL(eigcon(itau,ia,:),4)
            ENDDO
            WRITE(eq1_unit)
         ENDDO
         CALL close_bin(eq1_unit,"eigcon.bin")
      ENDIF
c-----------------------------------------------------------------------
c     extrapolate eigenvector to magnetic axis.
c-----------------------------------------------------------------------
      ALLOCATE(gt%eigveco(SIZE(ii%eigvec%fs,3)))
      DO iqty=1,SIZE(ii%eigvec%fs,3)
         DO ia=1,4
            ff1(ia)=SUM(ii%eigvec%fs(iqty,1:ii%mtau,ia-1))/ii%mtau
         ENDDO
         CALL extrap(m,x,ff1,x0,gt%eigveco(iqty))
      ENDDO
c-----------------------------------------------------------------------
c     deallocate local arrays.
c-----------------------------------------------------------------------
      DEALLOCATE(psi)
      DEALLOCATE(f)
      DEALLOCATE(p)
      DEALLOCATE(q)
      DEALLOCATE(rcc)
      DEALLOCATE(zcc)
      DEALLOCATE(dprer)
      DEALLOCATE(dprei)
      DEALLOCATE(delrr)
      DEALLOCATE(delri)
      DEALLOCATE(delzr)
      DEALLOCATE(delzi)
      DEALLOCATE(xphir)
      DEALLOCATE(xphii)
      DEALLOCATE(bradr)
      DEALLOCATE(bradi)
      DEALLOCATE(baxlr)
      DEALLOCATE(baxli)
      DEALLOCATE(btorr)
      DEALLOCATE(btori)
      DEALLOCATE(xpsir)
      DEALLOCATE(xpsii)
      DEALLOCATE(xchir)
      DEALLOCATE(xchii)
      DEALLOCATE(bnrmr)
      DEALLOCATE(bnrmi)
      DEALLOCATE(bpolr)
      DEALLOCATE(bpoli)
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL process_eq(gt)
c-----------------------------------------------------------------------
c     terminate program.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 2. extrap.
c     extrapolates function to x=0.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE extrap(m,x,ff,x0,ff0)
      USE local
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: m
      REAL(r8), DIMENSION(m), INTENT(IN) :: x,ff
      REAL(r8), INTENT(IN) :: x0
      REAL(r8), INTENT(OUT) :: ff0
      
      INTEGER(i4) :: i,j
      REAL(r8) :: term
c-----------------------------------------------------------------------
c     computations.
c-----------------------------------------------------------------------
      ff0=0
      DO i=1,m
         term=ff(i)
         DO j=1,m
            IF(j.NE.i)term=term*(x0-x(j))/(x(i)-x(j))
         ENDDO
         ff0=ff0+term
      ENDDO
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END
c-----------------------------------------------------------------------
c     subprogram 3. transform
c     extrapolates function to x=0.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE transform(eigcon)
      USE local
      USE bicube
      USE inverse
      IMPLICIT NONE

      REAL(r8), DIMENSION(0:ii%mtau,0:ii%ma,8) :: eigcon

      INTEGER(i4) :: itau
      TYPE(bicube_type), TARGET :: rz
      REAL(r8), DIMENSION(0:ii%mtau,0:ii%ma) ::
     $     jacfac,psi_r,psi_z,tau_r,tau_z,delpsi,deltau
      REAL(r8), DIMENSION(:,:), POINTER :: r_psi,z_psi,r_tau,z_tau
c-----------------------------------------------------------------------
c     fit coordinates to bicubic splines.
c-----------------------------------------------------------------------
      CALL bicube_alloc(rz,ii%mtau,ii%ma,2_i4)
      rz%xs=(/(itau,itau=0,ii%mtau)/)*twopi/ii%mtau
      rz%ys=ii%sq_in%xs
      rz%fs(1,:,:)=ii%rg
      rz%fs(2,:,:)=ii%zg
      CALL bicube_fit(rz,"periodic","extrap")
c-----------------------------------------------------------------------
c     compute derivatives.
c-----------------------------------------------------------------------
      r_tau => rz%fsx(1,:,:)
      z_tau => rz%fsx(2,:,:)
      r_psi => rz%fsy(1,:,:)
      z_psi => rz%fsy(2,:,:)
c-----------------------------------------------------------------------
c     compute inverse derivatives.
c-----------------------------------------------------------------------
      jacfac=r_tau*z_psi-r_psi*z_tau
      tau_r=z_psi/jacfac
      tau_z=-r_psi/jacfac
      psi_r=-z_tau/jacfac
      psi_z=r_tau/jacfac
c-----------------------------------------------------------------------
c     multiply derivatives by gradients.
c-----------------------------------------------------------------------
      deltau=SQRT(tau_r**2+tau_z**2)
      r_tau=r_tau*deltau
      z_tau=z_tau*deltau
      delpsi=SQRT(psi_r**2+psi_z**2)
      r_psi=r_psi*delpsi
      z_psi=z_psi*delpsi
c-----------------------------------------------------------------------
c     compute new cylindrical components.
c-----------------------------------------------------------------------
      ii%eigvec%fs(3,:,:)=eigcon(:,:,1)*r_psi+eigcon(:,:,2)*r_tau
      ii%eigvec%fs(4,:,:)=eigcon(:,:,1)*z_psi+eigcon(:,:,2)*z_tau
      ii%eigvec%fs(6,:,:)=eigcon(:,:,3)*r_psi+eigcon(:,:,4)*r_tau
      ii%eigvec%fs(7,:,:)=eigcon(:,:,3)*z_psi+eigcon(:,:,4)*z_tau
      ii%eigvec%fs(9,:,:)=eigcon(:,:,5)*r_psi+eigcon(:,:,6)*r_tau
      ii%eigvec%fs(10,:,:)=eigcon(:,:,5)*z_psi+eigcon(:,:,6)*z_tau
      ii%eigvec%fs(12,:,:)=eigcon(:,:,7)*r_psi+eigcon(:,:,8)*r_tau
      ii%eigvec%fs(13,:,:)=eigcon(:,:,7)*z_psi+eigcon(:,:,8)*z_tau
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END
