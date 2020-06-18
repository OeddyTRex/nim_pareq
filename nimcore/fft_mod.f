c-----------------------------------------------------------------------
c     file fft_mod.f
c     module containing fast Fourier transform routines for use with
c     nimrod.  this is a hybrid version having the newer wrapper
c     retrofitted with Zoran Mikic's FFT routines.
c     (CRS and SJP--last revised 8/14/98 to make the ffts and
c     configuration-space operations scale with nlayers.  note that
c     the passed parameters, nf and ng are different than the
c     previous nx and ny.)
c
c     modified 12/10/98 to use a complex work array in cpftv.  CRS
c
c     modified 1/13/99 to handle complex Fourier coefficient arrays, as
c     directly as possible.  CRS
c
c     modified 9/16/99 for separating vector components and Fourier
c     components into separate array indices and having the toroidal
c     index last in the returned real-space arrays.
c
c     the option to perform ffts in nimrod's data structures without
c     dealiasing is being added.  CRS, 5/30/08
c
c     the no-dealiasing option is modified, so that it does not expect
c     the Fourier-coefficient arrays to have the m=2**lphi/2 mode.
c     CRS, 1/01/17
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     code organization.
c-----------------------------------------------------------------------
c     0. fft_mod.
c     1. fft_nim.
c     2. fft2d.
c     3. cpftv.
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      MODULE fft_mod
      USE local

      CONTAINS
c-----------------------------------------------------------------------
c     subprogram 1. fft_nim.
c     performs array packing, and calls the actual fft routine.
c     this version uses a complex array for the passed Fourier
c     components.  the data is three-dimensional, and only the last is
c     transformed.
c
c     the arguments for this subprogram are
c
c	direction:  either 'forward' or 'inverse'.  the former sums
c		functions of configuration space with exp{-inx}, and
c		divides by the number of cells for normalization.  the
c		latter sums Fourier coefficients with exp{+inx}.
c
c	nf:  dimension of the non-transformed direction of the array
c		of Fourier coefficients.
c
c	nr:  dimension of the non-transformed direction of the array
c		of configuration-space data.  If nr<nf, different
c		segments of the non-transformed dimension are allocated
c		to the different processor layers, and the config-space
c		data on a layer does not cover a full grid-block.
c
c	lphi:  log (base 2) of the number of cells of the transformable
c		dimension.  the corresponding number of dealiased modes
c		(for quadratic terms) is 2**lphi/3+1.
c
c	nq:  the number of quantities (vector components) represented
c		in the data field.
c
c	f_coef: complex 3D array (1:nq,1:nf,(N+1)), where N=2**lphi/3,
c		representing the n>=0 Fourier coefficients of nq real
c		functions of configuration space.
c
c	re: 3D array (1:nq,1:nr,1:2**lphi) for the data as
c		a function of configuration space.
c
c       dealiase: if this optional argument is present and set to
c		false, then the maximum Fourier coefficient is
c		N=2**lphi/2 in the fft routine and in the last array
c		dimension of f_coef.  this implies a different layer
c		decomposition for the complex coefficients and for the
c		real data.
c
c     the last two arrays are now assumed-size dummy arrays, so that the
c     non-transformed directions are packed together for better
c     efficiency.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE fft_nim(direction,nf,nr,lphi,nq,f_coef,re,dealiase)
      USE time
      USE mpi_nim
      USE pardata
      IMPLICIT NONE

      CHARACTER(*), INTENT(IN) :: direction
      INTEGER(i4), INTENT(IN) :: nf,nr,lphi,nq
      COMPLEX(r8), DIMENSION(nq,nf,*), INTENT(INOUT) :: f_coef
      REAL(r8), DIMENSION(nq,nr,*), INTENT(INOUT) :: re
      LOGICAL, INTENT(IN), OPTIONAL :: dealiase

      COMPLEX(r8), DIMENSION(:), ALLOCATABLE :: f_ctmp
      COMPLEX(r8), DIMENSION(:,:,:), ALLOCATABLE :: comp
      REAL(r8) :: timestart_fft,timeend_fft
      INTEGER(i4), DIMENSION(:), ALLOCATABLE :: recvcounts,displs,
     $             scounts,sdispls,rcounts,rdispls
      INTEGER(i4) :: isign,nphi,iq,im,im0,idp,idm,ierror,rep,nfl,
     $               nm_tot,il,nl,nm,rem,ip,ifst,ifen,npo,jfc,
     $               nmodes,ifc
c-----------------------------------------------------------------------
c     start timer and allocate arrays.
c-----------------------------------------------------------------------
      CALL timer(timestart_fft)
      nphi=2**lphi
      nm_tot=nphi/3+1
      IF (PRESENT(dealiase)) THEN
        IF (.NOT.dealiase) nm_tot=nphi/2  !  no 2*dx wavelength
      ENDIF
      nl=nprocs/nprocs_layer
      nm=nm_tot/nl
      nfl=nf/nl
      rem=MODULO(nm_tot,nl)
      rep=MODULO(nf,nl)
      nmodes=mode_hi-mode_lo+1
c-----------------------------------------------------------------------
c     check direction.
c-----------------------------------------------------------------------
      IF (direction/='forward'.AND.direction/='inverse') CALL nim_stop
     $  ('fft_nimc does not recognize direction '//direction//'.')
c-----------------------------------------------------------------------
c     pack arrays.
c-----------------------------------------------------------------------
      dir_if: IF (direction=='forward') THEN
c-----------------------------------------------------------------------
c       if lphi<=1, the only component is n=0.  there should only be
c       one processor-layer in these cases.
c-----------------------------------------------------------------------
        IF (lphi<=1) THEN
          f_coef(:,:,1)=re(:,:,1)
          CALL timer(timeend_fft)
          time_fft = time_fft + timeend_fft-timestart_fft
          RETURN
        ENDIF
        isign=1
        ALLOCATE(comp(nq,nr,nm_tot))
      ELSE dir_if
c-----------------------------------------------------------------------
c       if lphi<=1, the only component is n=0.
c-----------------------------------------------------------------------
        IF (lphi<=1) THEN
          DO ip=1,2**lphi
            re(:,:,ip)=f_coef(:,:,1)
          ENDDO
          CALL timer(timeend_fft)
          time_fft = time_fft + timeend_fft-timestart_fft
          RETURN
        ENDIF
c-----------------------------------------------------------------------
c       there are 3 possibilities:  1) Fourier coefficients are located
c       on separate layers but all layers need all of the configuration-
c       space data, 2) coefficients on separate layers and config-space
c       data is broken up into separate portions across the poloidal
c       plane (within a block), 3) there is only one layer.
c
c       1) if different Fourier coefficients are located on different
c       processor layers, collect all before packing.
c       recvcounts = # of values from each proc across layers
c       displs = offsets into comp array for each proc's section
c-----------------------------------------------------------------------
        isign=-1
        ALLOCATE(comp(nq,nr,nm_tot))
        layer_if: IF (nl>1.AND.nr==nf) THEN
          ALLOCATE(recvcounts(0:nl-1),displs(0:nl-1))
          displs(0)=0
          DO il=0,nl-1
            IF (il<rem) THEN
              recvcounts(il)=nq*nf*(nm+1)
            ELSE
              recvcounts(il)=nq*nf*nm
            ENDIF
            IF (il>0) displs(il)=displs(il-1)+recvcounts(il-1)
          ENDDO
          CALL mpi_allgatherv(f_coef,recvcounts(ilayer),mpi_nim_comp,
     $         comp,recvcounts,displs,mpi_nim_comp,comm_mode,ierror)
          DEALLOCATE(recvcounts,displs)
c-----------------------------------------------------------------------
c       2) transpose type communication.  first make a 1D
c       array holding contiguous information for communication and
c       the send&recv displacement and count arrays.
c-----------------------------------------------------------------------
        ELSE IF (nf>nr) THEN layer_if
          ALLOCATE(f_ctmp(nq*nf*nmodes))
          ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),
     $             rdispls(0:nl-1),rcounts(0:nl-1))
          sdispls(0)=0
          rdispls(0)=0
          jfc=1
          DO il=0,nl-1
            ifst=il*nfl+1+MIN(rep,il)
            ifen=(il+1)*nfl+MIN(rep,il+1_i4)
            npo=ifen-ifst+1
            scounts(il)=npo*nq*nmodes
            IF (il<rem) THEN
              rcounts(il)=(nm+1)*nq*nr
            ELSE
              rcounts(il)=nm*nq*nr
            ENDIF
            IF (il>0) THEN
              sdispls(il)=sdispls(il-1)+scounts(il-1)
              rdispls(il)=rdispls(il-1)+rcounts(il-1)
            ENDIF
            DO im=1,nmodes
              DO ifc=ifst,ifen
                f_ctmp(jfc:jfc+nq-1)=f_coef(1:nq,ifc,im)
                jfc=jfc+nq
              ENDDO
            ENDDO
          ENDDO
c-----------------------------------------------------------------------
c         communicate the data.
c-----------------------------------------------------------------------
          CALL mpi_alltoallv(f_ctmp,scounts,sdispls,mpi_nim_comp,
     $         comp,rcounts,rdispls,mpi_nim_comp,comm_mode,ierror)
          DEALLOCATE(f_ctmp,scounts,sdispls,rcounts,rdispls)
        ENDIF layer_if
      ENDIF dir_if
c-----------------------------------------------------------------------
c     perform the transform.
c-----------------------------------------------------------------------
      IF (nl==1) THEN
        CALL fft2d(re,f_coef,nr*nq,lphi,isign,dealiase)
      ELSE
        CALL fft2d(re,comp  ,nr*nq,lphi,isign,dealiase)
      ENDIF
c-----------------------------------------------------------------------
c     put Fourier component data into nimrod form if this is a 
c     forward transform.
c-----------------------------------------------------------------------
      IF (direction=='forward') THEN
c-----------------------------------------------------------------------
c       for multi-layer runs, collect all portions of the block in a 1D
c       array, then sort.  the first part of the if block rearranges
c       from layered configuration-space data (splitting the poloidal
c       directions) to layered Fourier coefficients (splitting
c       components).  the second part goes from unsplit config-space
c       data to split coefficients.
c-----------------------------------------------------------------------
        layer_if2: IF (nf>nr) THEN
          ALLOCATE(f_ctmp(nq*nf*nmodes))
          ALLOCATE(sdispls(0:nl-1),scounts(0:nl-1),
     $             rdispls(0:nl-1),rcounts(0:nl-1))
          rdispls(0)=0
          sdispls(0)=0
          DO il=0,nl-1
            ifst=il*nfl+1+MIN(rep,il)
            ifen=(il+1)*nfl+MIN(rep,il+1_i4)
            npo=ifen-ifst+1
            rcounts(il)=npo*nq*nmodes
            IF (il<rem) THEN
              scounts(il)=(nm+1)*nq*nr
            ELSE
              scounts(il)=nm*nq*nr
            ENDIF
            IF (il>0) THEN
              rdispls(il)=rdispls(il-1)+rcounts(il-1)
              sdispls(il)=sdispls(il-1)+scounts(il-1)
            ENDIF
          ENDDO
c-----------------------------------------------------------------------
c         gather this layer's Fourier components from all block
c         portions in comp.
c-----------------------------------------------------------------------
          CALL mpi_alltoallv(comp,scounts,sdispls,mpi_nim_comp,
     $         f_ctmp,rcounts,rdispls,mpi_nim_comp,comm_mode,ierror)
          DEALLOCATE(scounts,sdispls,rcounts,rdispls)
c-----------------------------------------------------------------------
c         select the data for this Fourier component layer.
c-----------------------------------------------------------------------
          jfc=1
          DO il=0,nl-1
            ifst=il*nfl+1+MIN(rep,il)
            ifen=(il+1)*nfl+MIN(rep,il+1_i4)
            npo =ifen-ifst+1
            DO im=1,nmodes
              DO ifc=ifst,ifen
                f_coef(1:nq,ifc,im)=f_ctmp(jfc:jfc+nq-1)
                jfc=jfc+nq
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(f_ctmp)
          IF (ilayer==0) f_coef(:,:,1)=REAL(f_coef(:,:,1),r8)
c-----------------------------------------------------------------------
c       complete config-space data to component on layers.
c-----------------------------------------------------------------------
        ELSE IF (nl>1) THEN layer_if2
          f_coef(1:nq,1:nf,1:nmodes)=comp(:,:,mode_lo:mode_hi)
          IF (ilayer==0) f_coef(:,:,1)=REAL(f_coef(:,:,1),r8)
        ENDIF layer_if2
      ENDIF
      DEALLOCATE(comp)
c-----------------------------------------------------------------------
c     complete timing.
c-----------------------------------------------------------------------
      CALL timer(timeend_fft)
      time_fft = time_fft + timeend_fft-timestart_fft
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fft_nim

c-----------------------------------------------------------------------
c     subprogram 2. fft2d
c
c     slices data into sections for cray vector machines.  an option
c     for non-vector machines may be desirable.
c
c     ###### CRAY-2 de-aliased version. ######
c
c     Performs the FFT of a two-dimensional array over
c     the second dimension, for all points in the first
c     dimension.  The operations over this first dimension are
c     vectorized.  This routine uses Oscar Buneman's program
c     VCFT.  Load the binary VCFT on the Cray-2.
c     (This binary is available from FILEM 1505 .BUTILITY)
c     WARNING: The CAL VCFT on the Cray-2 requires that MY be
c     between 4 and 21.
c
c     compr and compi have been combined into one array with both
c     real and imaginary parts to reduce data rearranging in nimrod.
c     (CRS, 8/14/98)
c
c     The work array is placed here, so that the loop-splitting
c     can be easily varied by changing the incf paramter.  The
c     indices of the 'imaginary' part of the data have been changed,
c     and the work array is now complex for optimization.
c     (CRS, 12/10/98)
c
c     This version of fft2d uses a complex comp array directly.
c     (CRS, 1/13/99)
c
c     This program can be compiled with CIVIC, CFT2, or CFT77.
c
c     The definition of the arguments is as follows:
c
c       IFLAG:   An integer flag which determines the direction of
c                the FFT.  Fourier analysis is performed when IFLAG=1,
c                and Fourier synthesis (i.e., back to real space) is
c                performed when IFLAG.ne.1.
c
c       NX:      The number of points in the non-FFT direction
c                (over which the vectorization is performed).
c
c       MY:      The power-of-two which gives the number of points
c                in the FFT direction.  Thus the number of points in
c                the FFT direction is NY=2**MY.  Note that MY needs
c                to be between 2 and 10, inclusive.
c
c       RE:      Real array dimensioned NX by NY which contains the
c                data in real space. This array is used as input
c                when IFLAG=1, and is output when IFLAG.ne.1.
c
c       COMP:    Real and imaginary part of complex array dimensioned
c                NX by 2*(NY/3+1) which contains the non-aliased
c                Fourier modes.  This array is used as input when
c                IFLAG.ne.1, and is output when IFLAG=1.
c
c       DEALIASE: When this optional input is present and set to F, the
c                 routine uses all possible Fourier coefficients.  See
c                 below.
c
c     Note that RE is not overwritten in the call with IFLAG=1,
c     and similarly, COMP is not overwritten in the call with
c     IFLAG.ne.1.
c
c     The operation of this routine can be summarized by:
c
c                           NY
c       COMP(I,M) = 1/NY * SUM RE(I,J)*EXP[-2*pi*i*(J-1)*(M-1)/NY]
c                          J=1
c
c                         for I=1,2,...,NX and M=1,2,...,NY/2+1
c                         when IFLAG=1, and
c
c                  NY
c       RE(I,J) = SUM COMP(I,M)*EXP[2*pi*i*(J-1)*(M-1)/NY]
c                 M=1
c
c                         for I=1,2,...,NX and J=1,2,...,NY
c                         when IFLAG.ne.1, where the elements COMP(I,M)
c                         for M=NY/2+2,...,NY are not stored in array
c                         COMP, but obey COMP(I,M)=conj(COMP(I,NY-M+2)).
c
c                         Standard operation for this routine de-aliases
c                         the COMP array such that only the
c                         modes with M=1,2,...,NY/3+1 are nonzero.
c                         However, if the optional dealiase flag is
c                         present and set to false, all (1<=M<=NY/2+1)
c                         are used.
c-----------------------------------------------------------------------
      SUBROUTINE fft2d(re,comp,nx,my,iflag,dealiase)
      IMPLICIT NONE

      INTEGER(i4), INTENT(IN) :: nx,my,iflag
      REAL(r8), DIMENSION(nx,*), INTENT(INOUT) :: re
      COMPLEX(r8), DIMENSION(nx,*), INTENT(INOUT) :: comp
      LOGICAL, INTENT(IN), OPTIONAL :: dealiase

      INTEGER(i4) :: ny,nyd3,n,n2,n2p,i,j,m,k,mm,i1,i2
      INTEGER(i4), PARAMETER :: incf=24,istride=2*incf-2
      COMPLEX(r8), DIMENSION(incf,2**my+1) :: work
      COMPLEX(r8) :: ctmp
      REAL(r8) :: xnrm,xnrm2
c-----------------------------------------------------------------------
c
      ny=2**my
      nyd3=ny/3+1
      if (present(dealiase)) then
        if (.not.dealiase) nyd3=ny/2   !  800 loop is skipped entirely
      endif
      xnrm=1._r8/ny
      xnrm2=xnrm*0.5
c
      if (iflag.eq.1) then
c
c     FORWARD TRANSFORM.
c
      do 400 i=1,nx,istride
      n=MIN(istride,nx-i+1_i4)
      n2=n/2
      n2p=n-n2
c
c     Pack two x-lines at a time.
c
      do 100 j=1,ny
      do 100 k=1,n2
      i1=i-2+2*k
      i2=i1+1
  100 work(k,j)=re(i1,j)+(0,1)*re(i2,j)
c
      if (n2<n2p) then
        do 110 j=1,ny
  110   work(n2p,j)=re(i-1+n,j)
      endif
c
      call cpftv(work,ny,incf,n2p,1_i4,-1_i4)
c
c     Unpack the two x-lines (normalization moved here).
c
cdir$ ivdep
      do 200 k=1,n2
      i1=i-2+2*k
      i2=i1+1
      ctmp=xnrm*work(k,1)
      comp(i1,1)=REAL(ctmp,r8)
  200 comp(i2,1)=AIMAG(ctmp)
c
      do 300 m=2,nyd3
      mm=ny+2-m
cdir$ ivdep
      do 300 k=1,n2
      i1=i-2+2*k
      i2=i1+1
      ctmp=CONJG(work(k,mm))
      comp(i1,m)=xnrm2*(ctmp+work(k,m))
  300 comp(i2,m)=xnrm2*(ctmp-work(k,m))*(0,1)
c
      if (n2<n2p) then
        i1=i-1+n
        comp(i1,1)=xnrm*REAL(work(n2p,1),r8)
        do 310 m=2,nyd3
        mm=ny+2-m
  310   comp(i1,m)=xnrm2*(CONJG(work(n2p,mm))+work(n2p,m))
      endif
c
  400 continue
c
      else
c
c     INVERSE TRANSFORM.
c
      do 950 i=1,nx,istride
      n=MIN(istride,nx-i+1_i4)
      n2=n/2
      n2p=n-n2
c
c     Pack two x-lines at atime.
c
c     m=0
c
      do 500 k=1,n2
      i1=i-2+2*k
      i2=i1+1
  500 work(k,1)=REAL(comp(i1,1),r8)+(0,1)*REAL(comp(i2,1),r8)
c
c     Positive and negative m.
c
      do 600 m=2,nyd3
      mm=ny+2-m
cdir$ ivdep
      do 600 k=1,n2
      i1=i-2+2*k
      i2=i1+1
      work(k,m )=comp(i1,m)+(0,1)*comp(i2,m)
  600 work(k,mm)=CONJG(comp(i1,m))+(0,1)*CONJG(comp(i2,m))
c
      if (n2<n2p) then
        i1=i-1+n
        work(n2p,1)=REAL(comp(i1,1),r8)
        do 610 m=2,nyd3
        mm=ny+2-m
        work(n2p,m )=comp(i1,m)
  610   work(n2p,mm)=CONJG(comp(i1,m))
      endif
c
c     De-aliased m.
c
      do 800 m=nyd3+1,ny-nyd3+1
      do 800 k=1,n2p
  800 work(k,m)=0.
c
      call cpftv(work,ny,incf,n2p,1_i4,1_i4)
c
      do 900 j=1,ny
cdir$ ivdep
      do 900 k=1,n2
      i1=i-2+2*k
      i2=i1+1
      re(i1,j)=work(k,j)
  900 re(i2,j)=AIMAG(work(k,j))
c
      if (n2<n2p) then
        do 910 j=1,ny
        i1=i-1+n
  910   re(i1,j)=work(n2p,j)
      endif
c
  950 continue
c
      end if
c
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE fft2d

c-----------------------------------------------------------------------
c     subprogram 3. cpftv
c     perform transforms.
c
c     This subroutine performs an FFT in one direction of a
c     two-dimensional array, for each row in the non-FFT direction.
c     The loops are vectorized in this perpendicular direction.
c     This routine is a modified version of Langdon's CPFT
c     FFT routine (which is based on Singleton's FFT algorithm).
c     This routine is unnormalized, in the sense that if it is
c     called twice in succession (with the same array), first with
c     ISIGN=1 and then with ISIGN=-1, the array is multiplied by N.
c
c     The meaning of the arguments is as follows:
c     CA    = complex array containing the the complex
c             data to be Fourier-transformed,
c     N     = length of the FFT (must be a power of 2),
c     INCF  = increment between successive elements of R (and I)
c             in the FFT direction when it is regarded
c             as a real one-dimensional array,
c     NV    = number of rows over which to perform the FFT
c             (i.e., in the vector direction),
c     INCV  = increment between successive elements of R (and I)
c             in the perpendicular direction when it is regarded
c             as a real one-dimensional array,
c     ISIGN = sign to be used in the exponential of the FFT
c             (use +1 or -1).
c
c#### Modified by Z. Mikic, SAIC, La Jolla, August 18, 1985.
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     declarations.
c-----------------------------------------------------------------------
      SUBROUTINE cpftv (ca,n,incf,nv,incv,isign)
      IMPLICIT REAL(r8) (a-h,o-z)

      COMPLEX(r8), DIMENSION(:,:), INTENT(INOUT) :: ca
      INTEGER(i4), INTENT(IN) :: n,incf,nv,incv,isign

      INTEGER(i4) :: span,rc
      INTEGER(i4), PARAMETER :: log2nx=15
      REAL(r8), DIMENSION(log2nx), SAVE :: sines=0
      REAL(r8) :: i0,i1,qt,qq
      COMPLEX(r8) :: c0,c1,ct
c-----------------------------------------------------------------------
c
      if(sines(1).eq.1.) go to 1
      sines(1)=1.
      qt=1.
      qt=atan(qt)
      do 2 is=2,log2nx
      qq=sin(qt)
      sines(is)=qq
    2 qt=qt*.5
    1 continue
c
      if(n.eq.1) return
      nvl=1+(nv-1)*incv
c*ZM* For 2-d arrays R and I, set INC=1.
      inc=1
      sgn=isign
      ninc=n*inc
      span=ninc
      it=n/2
      do 1000 is=1,log2nx
c... (2000=recur)
      if(it.eq.1) go to 2000
 1000 it=it/2
c
c  if truncated rather than rounded arithmetic is used,
c  singleton's magnitude correcton should be applied to cos and sin.
 1500 t=sinx+(s*cosx-c*sinx)
      cosx=cosx-(c*cosx+s*sinx)
      sinx=t
c... (3000=repl)
 3000 k1=k0+span
c*ZM* Vector loops introduced all end with 01.
cdir$ ivdep
      do 101 iv=1,nvl,incv
      c0=ca(iv,1+k0)
      c1=ca(iv,1+k1)
      ca(iv,1+k0)=c0+c1
      c0=c0-c1
      ca(iv,1+k1)=(cosx+(0,1)*sinx)*c0
  101 continue
      k0=k1+span
      if(k0.lt.ninc) go to 3000
      k1=k0-ninc
      cosx=-cosx
      k0=span-k1
      if(k1.lt.k0) go to 3000
      k0=k0+inc
      k1=span-k0
      if(k0.lt.k1) go to 1500
 2000 continue
      span=span/2
      k0=0
c... (4000=zero)
 4000 k1=k0+span
cdir$ ivdep
      do 201 iv=1,nvl,incv
      c0=ca(iv,1+k0)
      c1=ca(iv,1+k1)
      ca(iv,1+k0)=c0+c1
      ca(iv,1+k1)=c0-c1
  201 continue
      k0=k1+span
      if(k0.lt.ninc) go to 4000
      if(span.eq.inc) go to 5000
      k0=span/2
 4500 k1=k0+span
cdir$ ivdep
      do 301 iv=1,nvl,incv
      c0=ca(iv,1+k0)
      c1=ca(iv,1+k1)
      ca(iv,1+k0)=c0+c1
      ca(iv,1+k1)=sgn*(0,1)*(c0-c1)
  301 continue
      k0=k1+span
      if(k0.lt.ninc) go to 4500
      k1=inc+inc
      if(span.eq.k1) go to 2000
      c=2.*sines(is)**2
      is=is-1
      sinx=sign(sines(is),sgn)
      s=sinx
      cosx=1.-c
      k0=inc
      go to 3000
c
 5000 n1=ninc-inc
      n2=ninc/2
      ij=0
      ji=0
      rc=0
      if(n2.eq.inc) return
      go to 5020
c... (5010=even)
 5010 ij=n1-ij
      ji=n1-ji
cdir$ ivdep
      do 401 iv=1,nvl,incv
      ct=ca(iv,1+ij)
      ca(iv,1+ij)=ca(iv,1+ji)
      ca(iv,1+ji)=ct
  401 continue
      if(ij.gt.n2) go to 5010
c... (5020=odd)
 5020 ij=ij+inc
      ji=ji+n2
cdir$ ivdep
      do 501 iv=1,nvl,incv
      ct=ca(iv,1+ij)
      ca(iv,1+ij)=ca(iv,1+ji)
      ca(iv,1+ji)=ct
  501 continue
      it=n2
c... (6000=incrv)
 6000 it=it/2
      rc=rc-it
      if(rc.ge.0) go to 6000
      rc=rc+2*it
      ji=rc
      ij=ij+inc
      if(ij.le.ji) go to 5010
      if(ij.lt.n2) go to 5020
c-----------------------------------------------------------------------
c     terminate.
c-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE cpftv

c-----------------------------------------------------------------------
c     close module.
c-----------------------------------------------------------------------
      END MODULE fft_mod
