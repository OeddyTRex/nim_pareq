c-----------------------------------------------------------------------
c
c	FFT routines originally out of fftpack from netlib.
c       Version 4  april 1985, Paul N Swarztrauber, NCAR boulder,colorado 80307
c       8 Dec 2000: modified by T.A. Gianakon for use with NIMROD.
c      CONTAINS
c         0.   rfft_mod
c         1.   rffti     initialize  rfftf and rfftb
c         2.   rfftf     forward transform of a real periodic sequence
c         3.     RADFG
c         4.     RADF5
c         5.     RADF3
c         6.     RADF2
c         7.     RADF4
c         4.   rfftb     backward transform of a real coefficient array
c-----------------------------------------------------------------------
c      SUBROUTINE rffti initializes the array wsave which is used in
c      both rfftf and rfftb. the prime factorization of n together with
c      a tabulation of the trigonometric functions are computed and
c      stored in wsave.
c        input parameter:
c           n       the length of the sequence to be transformed.
c        output parameter
c           wsave   a work array which must be dimensioned at least 2*n+15.
c                   the same work array can be used for both rfftf and rfftb
c                   as long as n remains unchanged. different wsave arrays
c                   are required for different values of n. the contents of
c                   wsave must not be changed between calls of rfftf or rfftb.
c-----------------------------------------------------------------------
c      SUBROUTINE rfftf computes the fourier coefficients of a real
c      perodic sequence (fourier analysis). the transform is defined
c      below at output parameter r.
c        input parameters:
c
c      n       the length of the array r to be transformed.  the method
c              is most efficient when n is a product of small primes.
c              n may change so long as different work arrays are provided
c      
c      r       a real array of length n which contains the sequence
c                    to be transformed
c      
c      wsave   a work array which must be dimensioned at least 2*n+15.
c              in the program that calls rfftf. the wsave array must be
c              initialized by calling subroutine rffti(n,wsave) and a
c              different wsave array must be used for each different
c              value of n. this initialization does not have to be
c              repeated so long as n remains unchanged thus subsequent
c              transforms can be obtained faster than the first.
c              the same wsave array can be used by rfftf and rfftb.
c      
c
c      output parameters
c
c      r       r(1) = the sum from i=1 to i=n of r(i)
c
c              if n is even set l =n/2   , if n is odd set l = (n+1)/2
c      
c                then for k = 2,...,l
c
c                   r(2*k-2) = the sum from i = 1 to i = n of
c
c                        r(i)*cos((k-1)*(i-1)*2*pi/n)
c
c                   r(2*k-1) = the sum from i = 1 to i = n of
c
c                       -r(i)*sin((k-1)*(i-1)*2*pi/n)
c
c              if n is even
c
c                   r(n) = the sum from i = 1 to i = n of
c
c                        (-1)**(i-1)*r(i)
c
c       *****  note
c                   this transform is unnormalized since a call of rfftf
c                   followed by a call of rfftb will multiply the input
c                   sequence by n.
c
c      wsave   contains results which must not be destroyed between
c              calls of rfftf or rfftb.
c      
c-----------------------------------------------------------------------
c      SUBROUTINE rfftb computes the real perodic sequence from its
c      fourier coefficients (fourier synthesis). the transform is defined
c      below at output parameter r.
c        input parameters
c      
c      n       the length of the array r to be transformed.  the method
c              is most efficient when n is a product of small primes.
c              n may change so long as different work arrays are provided
c      
c      r       a real array of length n which contains the sequence
c              to be transformed
c      
c      wsave   a work array which must be dimensioned at least 2*n+15.
c              in the program that calls rfftb. the wsave array must be
c              initialized by calling subroutine rffti(n,wsave) and a
c              different wsave array must be used for each different
c              value of n. this initialization does not have to be
c              repeated so long as n remains unchanged thus subsequent
c              transforms can be obtained faster than the first.
c              the same wsave array can be used by rfftf and rfftb.
c      
c      
c        output parameters
c      
c      r       for n even and for i = 1,...,n
c      
c                   r(i) = r(1)+(-1)**(i-1)*r(n)
c      
c                        plus the sum from k=2 to k=n/2 of
c      
c                         2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c      
c                        -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c      
c              for n odd and for i = 1,...,n
c      
c                   r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of
c      
c                        2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)
c      
c                       -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)
c      
c       *****  note
c                   this transform is unnormalized since a call of rfftf
c                   followed by a call of rfftb will multiply the input
c                   sequence by n.
c      
c      wsave   contains results which must not be destroyed between
c              calls of rfftb or rfftf.
c      

      MODULE rfft_mod
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: nfft
      REAL(r8), DIMENSION(:),ALLOCATABLE  :: WA,CH,C
      INTEGER(i4), DIMENSION(15) :: IFAC

      CONTAINS
c-----------------------------------------------------------------------
      SUBROUTINE RFFTI
      USE local
      IMPLICIT NONE
      INTEGER(i4), DIMENSION(4),PARAMETER :: NTRYH=(/4,2,3,5/)
      REAL(r8) :: ARG,ARGH,ARGLD,FI
      INTEGER(i4) :: I,IB,IDO,II,IP,IPM,IS,J,K1,L1,L2,LD
      INTEGER(i4) :: N,NF,NFM1,NL,NQ,NR,NTRY


      n = nfft
      IF (N .EQ. 1) RETURN
      
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
      ENDDO
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = twopi/REAL(N,r8)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN

      DO K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = REAL(LD,r8)*ARGH
            FI = 0.
            DO II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            ENDDO
            IS = IS+IDO
         ENDDO
         L1 = L2
      ENDDO
      RETURN
      END SUBROUTINE RFFTI


      SUBROUTINE RFFTF 
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IDL1,IDO,IP,IW,IX2,IX3,IX4,K1,KH,L1,L2,N,NA,NF
      N = nfft
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END SUBROUTINE RFFTF


c     SUBROUTINE RFFTB (N,C,CH,WA,IFAC)
c     CALL RFFTB1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      SUBROUTINE RFFTB
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IDL1,IDO,IP,IW,IX2,IX3,IX4,K1,L1,L2,N,NA,NF
      N = nfft
      IF (N .EQ. 1) RETURN
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL RADB2 (IDO,L1,C,CH,WA(IW))
         GO TO 105
  104    CALL RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    CALL RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
         CALL RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    CALL RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    CALL RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      DO 117 I=1,N
         C(I) = CH(I)
  117 CONTINUE
      RETURN
      END SUBROUTINE RFFTB

      END MODULE rfft_mod



      SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      USE local
      IMPLICIT NONE
      REAL(r8) :: AI1,AI2,AR1,AR1H,AR2,AR2H,ARG 
      REAL(r8) :: CI2,CI3,CI4,CI5,CR2,CR3,CR4,CR5 
      REAL(r8) :: DC2,DCP,DI2,DI3,DI4,DI5,DR2,DR3,DR4,DR5,DS2,DSP
      INTEGER(i4) ::  I,IC,IDIJ,IDL1,IDO,IDP2,IK,IP,IPP2,IPPH,IS 
      INTEGER(i4) :: J,J2,JC,K,L,L1,LC,NBD 
      REAL(r8) :: TAUI,TAUR,TI11,TI12,TI2,TI3,TI4,TI5 
      REAL(r8) :: TR11,TR12,TR2,TR3,TR4,TR5,WA1,WA2,WA3,WA4
      REAL(r8), DIMENSION(IDO,IP,L1) :: CC
      REAL(r8), DIMENSION(IDO,L1,IP) :: C1
      REAL(r8), DIMENSION(IDL1,IP) :: C2,CH2
      REAL(r8), DIMENSION(IDO,L1,IP) :: CH
      REAL(r8), DIMENSION(*) :: WA
      ARG = twopi/REAL(IP,r8)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END SUBROUTINE RADFG


      SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      USE local
      IMPLICIT NONE
      REAL(r8) :: CI2,CI3,CI4,CI5,CR2,CR3,CR4,CR5 
      REAL(r8) :: DI2,DI3,DI4,DI5,DR2,DR3,DR4,DR5
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1
      REAL(r8) :: TI2,TI3,TI4,TI5,TR2,TR3,TR4,TR5
      REAL(r8) :: TR11 = 0.309016994374947
      REAL(r8) :: TI11 = 0.951056516295154
      REAL(r8) :: TR12 =-0.809016994374947
      REAL(r8) :: TI12 = 0.587785252292473
      REAL(r8), DIMENSION(*) :: WA1,WA2,WA3,WA4
      REAL(r8),DIMENSION(IDO,L1,5)  :: CC
      REAL(r8),DIMENSION(IDO,5,L1)  :: CH
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE RADF5



      SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)
      USE local
      IMPLICIT NONE
      REAL(r8) :: CI2,CR2,DI2,DI3,DR2,DR3
      INTEGER(i4) ::  I,IC,IDO,IDP2,K,L1
      REAL(r8) :: TI2,TI3,TR2,TR3
      REAL(r8) :: TAUI= 0.866025403784439
      REAL(r8) :: TAUR = -0.5
      REAL(r8), DIMENSION(*) :: WA1,WA2

      REAL(r8),DIMENSION(IDO,3,L1)  :: CH
      REAL(r8),DIMENSION(IDO,L1,3)  :: CC
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE RADF3



      SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1 
      REAL(r8),DIMENSION(IDO,2,L1) :: CH
      REAL(r8),DIMENSION(IDO,L1,2) :: CC
      REAL(r8),DIMENSION(*) :: WA1
      REAL(r8) :: TI2,TR2
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2_i4) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END SUBROUTINE RADF2


      SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      USE local
      IMPLICIT NONE

      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1 
      REAL(r8),DIMENSION(IDO,L1,4) :: CC
      REAL(r8),DIMENSION(IDO,4,L1) :: CH
      REAL(r8),DIMENSION(*) :: WA1,WA2,WA3
      REAL(r8) :: HSQT2=0.7071067811865475
      REAL(r8) :: CI2,CI3,CI4,CR2,CR3,CR4 
      REAL(r8) :: TI1,TI2,TI3,TI4,TR1,TR2,TR3,TR4

      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2_i4) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END SUBROUTINE RADF4


      SUBROUTINE RADBG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDIJ,IDL1,IDO,IDP2,IK,IP,IPP2,IPPH 
      INTEGER(i4) :: IS,J,J2,JC,K,L,L1,LC,NBD
      REAL(r8),DIMENSION(IDO,L1,IP) :: CH
      REAL(r8),DIMENSION(IDO,IP,L1) :: CC
      REAL(r8),DIMENSION(IDO,L1,IP) :: C1
      REAL(r8),DIMENSION(IDL1,IP) :: C2
      REAL(r8),DIMENSION(IDL1,IP) :: CH2
      REAL(r8),DIMENSION(*) :: WA
      REAL(r8) :: AI1,AI2,AR1,AR1H,AR2,AR2H,ARG,DC2,DCP,DS2,DSP 
      ARG = twopi/REAL(IP,r8)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            CH(1,K,J) = CC(IDO,J2-2,K)+CC(IDO,J2-2,K)
            CH(1,K,JC) = CC(1,J2-1,K)+CC(1,J2-1,K)
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               CH(I-1,K,J) = CC(I-1,2*J-1,K)+CC(IC-1,2*J-2,K)
               CH(I-1,K,JC) = CC(I-1,2*J-1,K)-CC(IC-1,2*J-2,K)
               CH(I,K,J) = CC(I,2*J-1,K)-CC(IC,2*J-2,K)
               CH(I,K,JC) = CC(I,2*J-1,K)+CC(IC,2*J-2,K)
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+AR1*CH2(IK,2)
            C2(IK,LC) = AI1*CH2(IK,IP)
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+AR2*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)+AI2*CH2(IK,JC)
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            CH(1,K,J) = C1(1,K,J)-C1(1,K,JC)
            CH(1,K,JC) = C1(1,K,J)+C1(1,K,JC)
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               CH(I-1,K,J) = C1(I-1,K,J)-C1(I,K,JC)
               CH(I-1,K,JC) = C1(I-1,K,J)+C1(I,K,JC)
               CH(I,K,J) = C1(I,K,J)+C1(I-1,K,JC)
               CH(I,K,JC) = C1(I,K,J)-C1(I-1,K,JC)
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            C1(1,K,J) = CH(1,K,J)
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)-WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)+WA(IDIJ)*CH(I-1,K,J)
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END SUBROUTINE RADBG 


      SUBROUTINE RADB5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1 
      REAL(r8),DIMENSION(IDO,5,L1) :: CC          
      REAL(r8),DIMENSION(IDO,L1,5) :: CH
      REAL(r8),DIMENSION(*) :: WA1,WA2,WA3,WA4
      REAL(r8),PARAMETER :: TR11 = 0.309016994374947
      REAL(r8),PARAMETER :: TI11 = 0.951056516295154
      REAL(r8),PARAMETER :: TR12 =-0.809016994374947
      REAL(r8),PARAMETER :: TI12 = 0.587785252292473
      REAL(r8) :: CI2,CI3,CI4,CI5,CR2,CR3,CR4,CR5 
      REAL(r8) :: DI2,DI3,DI4,DI5,DR2,DR3,DR4,DR5
      REAL(r8) :: TI2,TI3,TI4,TI5,TR2,TR3,TR4,TR5
      DO 101 K=1,L1
         TI5 = CC(1,3,K)+CC(1,3,K)
         TI4 = CC(1,5,K)+CC(1,5,K)
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         TR3 = CC(IDO,4,K)+CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI5 = TI11*TI5+TI12*TI4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(1,K,5) = CR2+CI5
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TI5 = CC(I,3,K)+CC(IC,2,K)
            TI2 = CC(I,3,K)-CC(IC,2,K)
            TI4 = CC(I,5,K)+CC(IC,4,K)
            TI3 = CC(I,5,K)-CC(IC,4,K)
            TR5 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            TR4 = CC(I-1,5,K)-CC(IC-1,4,K)
            TR3 = CC(I-1,5,K)+CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4-WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4+WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5-WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5+WA4(I-1)*DR5
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE RADB5 


      SUBROUTINE RADB3 (IDO,L1,CC,CH,WA1,WA2)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1
      REAL(r8),DIMENSION(IDO,3,L1) ::  CC
      REAL(r8),DIMENSION(IDO,L1,3) ::  CH
      REAL(r8),DIMENSION(*) :: WA1,WA2
      REAL(r8),PARAMETER :: TAUR = -.5
      REAL(r8),PARAMETER :: TAUI = 0.866025403784439
      REAL(r8) :: CI2,CI3,CR2,CR3,DI2,DI3,DR2,DR3 
      REAL(r8) :: TI2,TR2
      DO 101 K=1,L1
         TR2 = CC(IDO,2,K)+CC(IDO,2,K)
         CR2 = CC(1,1,K)+TAUR*TR2
         CH(1,K,1) = CC(1,1,K)+TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2-CI3
         CH(1,K,3) = CR2+CI3
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
            CR2 = CC(I-1,1,K)+TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K)+TR2
            TI2 = CC(I,3,K)-CC(IC,2,K)
            CI2 = CC(I,1,K)+TAUR*TI2
            CH(I,K,1) = CC(I,1,K)+TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2-CI3
            DR3 = CR2+CI3
            DI2 = CI2+CR3
            DI3 = CI2-CR3
            CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
      RETURN
      END SUBROUTINE RADB3 


      SUBROUTINE RADB2 (IDO,L1,CC,CH,WA1)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1 
      REAL(r8),DIMENSION(IDO,2,L1) :: CC 
      REAL(r8),DIMENSION(IDO,L1,2) :: CH
      REAL(r8),DIMENSION(*) :: WA1
      REAL(r8) :: TI2,TR2
      DO 101 K=1,L1
         CH(1,K,1) = CC(1,1,K)+CC(IDO,2,K)
         CH(1,K,2) = CC(1,1,K)-CC(IDO,2,K)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CH(I-1,K,1) = CC(I-1,1,K)+CC(IC-1,2,K)
            TR2 = CC(I-1,1,K)-CC(IC-1,2,K)
            CH(I,K,1) = CC(I,1,K)-CC(IC,2,K)
            TI2 = CC(I,1,K)+CC(IC,2,K)
            CH(I-1,K,2) = WA1(I-2)*TR2-WA1(I-1)*TI2
            CH(I,K,2) = WA1(I-2)*TI2+WA1(I-1)*TR2
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2_i4) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(IDO,K,1) = CC(IDO,1,K)+CC(IDO,1,K)
         CH(IDO,K,2) = -(CC(1,2,K)+CC(1,2,K))
  106 CONTINUE
  107 RETURN
      END SUBROUTINE RADB2 


      SUBROUTINE RADB4 (IDO,L1,CC,CH,WA1,WA2,WA3)
      USE local
      IMPLICIT NONE
      INTEGER(i4) :: I,IC,IDO,IDP2,K,L1 
      REAL(r8),DIMENSION(IDO,4,L1)  :: CC
      REAL(r8),DIMENSION(IDO,L1,4)  :: CH
      REAL(r8),DIMENSION(*) :: WA1,WA2,WA3
      REAL(r8) :: CI2,CI3,CI4,CR2,CR3,CR4 
      REAL(r8), PARAMETER :: SQRT2 = 1.414213562373095
      REAL(r8) :: TI1,TI2,TI3,TI4,TR1,TR2,TR3,TR4

      DO 101 K=1,L1
         TR1 = CC(1,1,K)-CC(IDO,4,K)
         TR2 = CC(1,1,K)+CC(IDO,4,K)
         TR3 = CC(IDO,2,K)+CC(IDO,2,K)
         TR4 = CC(1,3,K)+CC(1,3,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,2) = TR1-TR4
         CH(1,K,3) = TR2-TR3
         CH(1,K,4) = TR1+TR4
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TI1 = CC(I,1,K)+CC(IC,4,K)
            TI2 = CC(I,1,K)-CC(IC,4,K)
            TI3 = CC(I,3,K)-CC(IC,2,K)
            TR4 = CC(I,3,K)+CC(IC,2,K)
            TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
            TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
            TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
            TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1-TR4
            CR4 = TR1+TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
            CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
            CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
            CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
            CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
            CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2_i4) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         TI1 = CC(1,2,K)+CC(1,4,K)
         TI2 = CC(1,4,K)-CC(1,2,K)
         TR1 = CC(IDO,1,K)-CC(IDO,3,K)
         TR2 = CC(IDO,1,K)+CC(IDO,3,K)
         CH(IDO,K,1) = TR2+TR2
         CH(IDO,K,2) = SQRT2*(TR1-TI1)
         CH(IDO,K,3) = TI2+TI2
         CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
      END SUBROUTINE RADB4 

