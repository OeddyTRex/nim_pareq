c-----------------------------------------------------------------------
c     subprogram 1. read_eq_chease.
c     reads data from chease.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_chease(gt)
      USE local
      USE inverse
      USE input
      USE global
      USE physdat
      IMPLICIT NONE

      TYPE(global_type), INTENT(OUT) :: gt

      INTEGER(i4) :: ntnova,npsi1,isym,mtau,itau,ia
      REAL(r8), DIMENSION(5) :: axx
      REAL(r8), DIMENSION(:), POINTER :: zcpr,zcppr,zq,zdq,ztmf,
     $     ztp,zfb,zfbp,zpsi,zpsim
      REAL(r8), DIMENSION(:,:), POINTER :: zrcp,zzcp,zjacm,zjac
c-----------------------------------------------------------------------
c     open file and read sizes.
c-----------------------------------------------------------------------
      CALL open_bin(eq1_unit,TRIM(filename),"UNKNOWN","REWIND",32_i4)
      READ(eq1_unit)ntnova,npsi1,isym	          ! last data is symmettry
      READ(eq1_unit)axx
      IF(isym == 1)CALL nim_stop
     & ("read_chease:  Cannot process symmetric data")
c-----------------------------------------------------------------------
c     allocate local arrays.
c-----------------------------------------------------------------------
      ALLOCATE(zcpr(npsi1-1),zcppr(npsi1),zq(npsi1),zdq(npsi1),
     $     ztmf(npsi1),ztp(npsi1),zfb(npsi1),zfbp(npsi1),zpsi(npsi1),
     $     zpsim(npsi1-1))
      ALLOCATE(zrcp(ntnova+3,npsi1),zzcp(ntnova+3,npsi1),
     $     zjacm(ntnova+3,npsi1),zjac(ntnova+3,npsi1))
c-----------------------------------------------------------------------
c     read arrays and close file.
c-----------------------------------------------------------------------
      READ(eq1_unit)zcpr			! Pressure
      READ(eq1_unit)zcppr			! d Pressure / d psi
      READ(eq1_unit)zq				! q
      READ(eq1_unit)zdq				! d q / d psi
      READ(eq1_unit)ztmf
      READ(eq1_unit)ztp
      READ(eq1_unit)zfb				! F/q
      READ(eq1_unit)zfbp			! d (F/q ) / d psi
      READ(eq1_unit)zpsi
      READ(eq1_unit)zpsim
      READ(eq1_unit)zrcp
      READ(eq1_unit)zzcp
      READ(eq1_unit)zjacm
      READ(eq1_unit)zjac
      CALL close_bin(eq1_unit)
c-----------------------------------------------------------------------
c     copy 1D arrays.
c-----------------------------------------------------------------------
      ii%ma=npsi1-2
      CALL spline_alloc(ii%sq_in,ii%ma,4_i4)
      ii%psio=zpsi(npsi1)-zpsi(1)
      ii%sq_in%xs=(zpsi(2:)-zpsi(1))/ii%psio
      ii%sq_in%fs(:,1)=ztmf(2:)
      ii%sq_in%fs(:,2)=zcppr(2:)
      ii%sq_in%fs(:,3)=zq(2:)
      ii%sq_in%fs(:,4)=0
      ia=ii%ma
      CALL spline_fit(ii%sq_in,"extrap")
c-----------------------------------------------------------------------
c     integrate pressure.
c-----------------------------------------------------------------------
      CALL spline_int(ii%sq_in)
      ii%sq_in%fs(:,2)=(ii%sq_in%fsi(:,2)-ii%sq_in%fsi(ii%ma,2))*ii%psio
      CALL spline_fit(ii%sq_in,"extrap")
c-----------------------------------------------------------------------
c     copy 2D arrays.
c-----------------------------------------------------------------------
      ii%mtau=ntnova
      ii%ro=SUM(zrcp(1:ntnova,1))/ntnova
      ii%zo=SUM(zzcp(1:ntnova,1))/ntnova
      ALLOCATE(ii%rg(0:ii%mtau,0:ii%ma))
      ALLOCATE(ii%zg(0:ii%mtau,0:ii%ma))
c     DO ia=0,ii%ma
c       DO itau=0,ii%mtau
c         ii%rg(itau,ia)=zrcp(itau+1,ia+1)
c         ii%zg(itau,ia)=zzcp(itau+1,ia+1)
c       ENDDO
c     ENDDO
      ii%rg=zrcp(1:ntnova+1,2:)
      ii%zg=zzcp(1:ntnova+1,2:)
      OPEN(UNIT = 7, FILE='temp.dat',STATUS='UNKNOWN')
      DO ia=0,ii%ma
        write(7,*)
        DO itau=0,ii%mtau
          write(7,fmt='(x,e12.5,x,e12.5,x,e12.5)')
     &        REAL(itau)
     &       ,ii%rg(itau,ia)-ii%ro
     &       ,ii%zg(itau,ia)-ii%zo
        ENDDO
      ENDDO
      CLOSE(7)


c     CALL bicube_alloc(rz_in,ma,mtau,2)
c     rz_in%fs(:,:,1)=TRANSPOSE(zrcp(1:ntnova+1,2:))
c     rz_in%fs(:,:,2)=TRANSPOSE(zzcp(1:ntnova+1,2:))
c-----------------------------------------------------------------------
c     deallocate local arrays and process inverse equilibrium.
c-----------------------------------------------------------------------
      DEALLOCATE(zcpr,zcppr,zq,zdq,ztmf,ztp,zfb,zfbp,zpsi,zpsim)
      DEALLOCATE(zrcp,zzcp,zjacm,zjac)
c-----------------------------------------------------------------------
c     set flags.
c-----------------------------------------------------------------------
      ii%p1flag=.FALSE.
      ii%symflag=.FALSE.
c-----------------------------------------------------------------------
c     nullify eigenvector pointer, to note that this file type does
c     not contain eigenvector data.
c-----------------------------------------------------------------------
      NULLIFY(ii%eigvec)
c-----------------------------------------------------------------------
c     process equilibrium.
c-----------------------------------------------------------------------
      CALL process_eq(gt)
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE read_chease
