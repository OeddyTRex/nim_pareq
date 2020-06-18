c-----------------------------------------------------------------------
c     module 1. read_galkin.
c     reads data from Sergei Galkin's equilibrium code.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE read_galkin(gt)
      USE local
      USE inverse
      USE input
      USE global
      USE physdat
      IMPLICIT NONE

      TYPE(global_type), INTENT(OUT) :: gt

      LOGICAL, PARAMETER :: diagnose=.FALSE.
      INTEGER(i4) :: nfeq,n,m,i,j,ia,itau
      REAL(r8) :: dadv,anun,anu0,um,vm
      REAL(r8), DIMENSION(:), ALLOCATABLE :: uk,vk,psi,p,f,q,anu,pp,
     $     pff,a

      REAL(r8), DIMENSION(:,:), ALLOCATABLE :: ro_in
c-----------------------------------------------------------------------
c     read data.  the reals and integers are likely a mix of 64 and 32
c     bit data, so using open_bin will not help porting.
c-----------------------------------------------------------------------
      OPEN(UNIT=eq1_unit,FILE=TRIM(filename),FORM='UNFORMATTED',
     $     STATUS='OLD')
      READ(eq1_unit)nfeq,n,m,dadv,anun,anu0,um,vm
      ALLOCATE(uk(2:m-1))
      ALLOCATE(vk(2:m-1))
      ALLOCATE(p(n))
      ALLOCATE(f(n))
      ALLOCATE(q(n))
      ALLOCATE(anu(n))
      ALLOCATE(pp(n))
      ALLOCATE(pff(n))
      ALLOCATE(a(n))
      ALLOCATE(ro_in(n,2:m-1))
      READ(eq1_unit)uk,vk,ro_in,p,f,q,anu,pp,pff,a
      CLOSE(eq1_unit)
c-----------------------------------------------------------------------
c     integrate psi profile.
c-----------------------------------------------------------------------
      ALLOCATE(psi(n))
      psi(n)=1
      DO i=n,2,-1
         psi(i-1)=psi(i)-anu(i)*(a(i)-a(i-1))
      ENDDO
c-----------------------------------------------------------------------
c     copy and modify 1D arrays.
c-----------------------------------------------------------------------
      CALL spline_alloc(ii%sq_in,n-2_i4,3_i4)
      ii%ma=n-2
      ii%ro=um
      ii%zo=vm
      ii%psio=anun
      ii%sq_in%xs=psi(2:n)
      ii%sq_in%fs(:,1)=f(2:n)
      ii%sq_in%fs(:,2)=p(2:n)
      ii%sq_in%fs(:,3)=q(2:n)
      ii%sq_in%fs(:,4)=0
c-----------------------------------------------------------------------
c     copy and modify 2D arrays.
c-----------------------------------------------------------------------
      ii%mtau=m-2
      ALLOCATE(ii%rg(0:ii%mtau,0:ii%ma))
      ALLOCATE(ii%zg(0:ii%mtau,0:ii%ma))
      DO ia=0,ii%ma
         i=ia+2
         DO itau=0,ii%mtau-1
            j=itau+2
            ii%rg(itau,ia)=um+ro_in(i,j)*(uk(j)-um)
            ii%zg(itau,ia)=vm+ro_in(i,j)*(vk(j)-vm)
         ENDDO
         ii%rg(ii%mtau,ia)=ii%rg(0,ia)
         ii%zg(ii%mtau,ia)=ii%zg(0,ia)
      ENDDO
c-----------------------------------------------------------------------
c     set flags.
c-----------------------------------------------------------------------
      ii%p1flag=.FALSE.
      ii%symflag=.FALSE.
c-----------------------------------------------------------------------
c     diagnose 1D input.
c-----------------------------------------------------------------------
      IF(diagnose)THEN
         OPEN(unit=eq2_unit,file='galkin1.bin',form='unformatted')
         DO i=2,n
            WRITE(eq2_unit)
     $           REAL(p(i),4),REAL(f(i),4),REAL(q(i),4),REAL(anu(i),4),
     $           REAL(pp(i),4),REAL(pff(i),4),REAL(a(i),4)
         ENDDO
         WRITE(eq2_unit)
         CLOSE(unit=eq2_unit)
c-----------------------------------------------------------------------
c     diagnose 2D input.
c-----------------------------------------------------------------------
         OPEN(unit=eq2_unit,file='galkin2.bin',form='unformatted')
         DO ia=0,ii%ma
            DO itau=0,ii%mtau
               WRITE(eq2_unit)
     $              REAL(ii%rg(itau,ia),4),REAL(ii%zg(itau,ia),4)
            ENDDO
            WRITE(eq2_unit)
         ENDDO
         CLOSE(unit=eq2_unit)
      ENDIF
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
      END SUBROUTINE read_galkin
