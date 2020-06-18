c-----------------------------------------------------------------------
c     complete elliptic integral expansions multiplied by 2/pi.
c     this code is adapted from Numerical Recipes, Press et al.
c-----------------------------------------------------------------------
      FUNCTION CEL(QQC,PP,AA,BB,tol,maxn) RESULT(CEL_result)
      USE local
      IMPLICIT NONE

      REAL(r8), INTENT(IN) :: QQC,PP,AA,BB,tol
      INTEGER(i4), INTENT(IN) :: maxn
      REAL(r8) :: CEL_result,A,B,P,E,EM,F,Q,G,QC,CA
      INTEGER(i4) :: countn

      CA=SQRT(tol)

      IF(QQC==0._r8) CALL nim_stop('Failure in CEL.')
      QC=ABS(QQC)
      A=AA
      B=BB
      P=PP
      E=QC
      EM=1._r8
      IF(P>0._r8) THEN
        P=SQRT(P)
        B=B/P
      ELSE
        F=QC*QC
        Q=1._r8-F
        G=1._r8-P
        F=F-P
        Q=Q*(B-A*P)
        P=SQRT(F/G)
        A=(A-B)/G
        B=-Q/(G*G*P)+A*P
      ENDIF

      countn=0
      DO
        F=A
        A=A+B/P
        G=E/P
        B=B+F*G
        B=B+B
        P=G+P
        G=EM
        EM=QC+EM
        IF (ABS(G-QC)<G*CA) EXIT
        QC=SQRT(E)
        QC=QC+QC
        E=QC*EM
        countn=countn+1
        IF (countn>maxn) CALL nim_stop("CEL not converging.")
      ENDDO

      CEL_result=0.5_r8*pi*(B+A*EM)/(EM*(EM+P))
      RETURN
      END FUNCTION CEL
