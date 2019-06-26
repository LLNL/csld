      SUBROUTINE SETK02(NA,A,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN THE 1/48 IRREDUCTIBLE WEDGE OF THE BRILLOUIN
C     ZONE FOR A BODY-CENTRED CUBIC LATTICE WITH LATTICE PARAMETER A,
C     SYMMETRY IS OH
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX)
      EQUIVALENCE (IVOL,AVOL)
      PI = 3.141592653589793238D0
      IF(NA.LE.0) GOTO 96
      IF(A.LE.0.0D0) GOTO 98
      NPTK = (NA+1)*(NA+2)*(2*NA+3)/6
      IF(NPTK.GT.NKMAX) STOP '*** <SETK02> NPTK EXCEEDS NKMAX ***'
      NTET = 2*NA**3
      IF(NTET.GT.NTMAX) STOP '*** <SETK02> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      DK=PI/A/NA
      N=2*NA
      NS2=NA
      WRITE(6,100) NPTK,NTET,N*DK,NS2*DK,NS2*DK
      W = 3.0D0/NA**3
      NPTK=0
      I=0
    1 J=0
    2 K=0
    3 NPTK=NPTK+1
C  NPTK = I*(I+1)*(I+2)/6 + J*(J+1)/2 + K + 1                  (I<=N/2)
C  NPTK = NTOT - (N-I+1)*(N-I+2)*(N-I+3)/6 + J*(J+1)/2 + K + 1 (I>=N/2)
C  WHERE    NTOT = (N+2)*(N+3)*(N+4)/24
      WK = W
      IF(I.EQ.J .OR. J.EQ.K .OR. K.EQ.I) WK = WK/2.0D0
      IF(I.EQ.J .AND. J.EQ.K) WK = WK/3.0D0
      IF(K.EQ.0) WK = WK/2.0D0
      IF(J.EQ.0) WK = WK/2.0D0
      IF(I+J.EQ.N) THEN
         WK = WK/2.0D0
         IF(J.EQ.K) WK = W/6.0D0
         IF(K.EQ.NA) WK = W/24.0D0
      ENDIF
      PTK(1,NPTK)=I*DK
      PTK(2,NPTK)=J*DK
      PTK(3,NPTK)=K*DK
      PTK(4,NPTK)=WK
      K=K+1
      IF(K.LE.J) GOTO 3
      J=J+1
      IF(J.LE.I .AND. J.LE.N-I) GOTO 2
      I=I+1
      IF(I.LE.N) GOTO 1
      PTK(4,1) = W/48.0D0
      PTK(4,NPTK) = W/48.0D0
C *** DEFINE THE TETRAHEDRA (NTET = (N**3)/4)
      NTET = 0
      I7=0
      I=0
    4 IX=(I+1)*(I+2)/2
      J=0
    5 K=0
    6 I7=I7+1
      I6=I7+IX
      I2=I6+J+1
      I1=I2+1
      I8=I7+1
      I5=I6+1
      I3=I7+J+1
      I4=I3+1
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I6
      IDEF(3,NTET)=I2
      IDEF(4,NTET)=I1
      IF(K.EQ.J) GOTO 7
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I6
      IDEF(3,NTET)=I5
      IDEF(4,NTET)=I1
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I8
      IDEF(3,NTET)=I5
      IDEF(4,NTET)=I1
    7 IF(J.EQ.I) GOTO 8
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I3
      IDEF(3,NTET)=I2
      IDEF(4,NTET)=I1
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I3
      IDEF(3,NTET)=I4
      IDEF(4,NTET)=I1
      IF(K.EQ.J) GOTO 8
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I8
      IDEF(3,NTET)=I4
      IDEF(4,NTET)=I1
    8 K=K+1
      IF(K.LE.J) GOTO 6
      J=J+1
      IF(J.LE.I) GOTO 5
      I=I+1
      IF(I.LT.NS2) GOTO 4
    9 IX=(N+1-I)*(N+2-I)/2
C     I7=NPTK-(N+1-I)*(N+2-I)*(N+3-I)/6
      I7=NPTK-IX*(N+3-I)/3
      J=0
   10 K=0
   11 I7=I7+1
      I6=I7+IX
      I2=I6+J+1
      I1=I2+1
      I8=I7+1
      I5=I6+1
      I3=I7+J+1
      I4=I3+1
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I3
      IDEF(3,NTET)=I4
      IDEF(4,NTET)=I6
      IF(K.EQ.J) GOTO 12
      NTET=NTET+1
      IDEF(1,NTET)=I7
      IDEF(2,NTET)=I8
      IDEF(3,NTET)=I4
      IDEF(4,NTET)=I6
      NTET=NTET+1
      IDEF(1,NTET)=I8
      IDEF(2,NTET)=I4
      IDEF(3,NTET)=I6
      IDEF(4,NTET)=I5
   12 IF(J.EQ.N-I-1) GOTO 13
      NTET=NTET+1
      IDEF(1,NTET)=I4
      IDEF(2,NTET)=I6
      IDEF(3,NTET)=I2
      IDEF(4,NTET)=I1
      NTET=NTET+1
      IDEF(1,NTET)=I3
      IDEF(2,NTET)=I4
      IDEF(3,NTET)=I6
      IDEF(4,NTET)=I2
      IF(K.EQ.J) GOTO 13
      NTET=NTET+1
      IDEF(1,NTET)=I4
      IDEF(2,NTET)=I6
      IDEF(3,NTET)=I5
      IDEF(4,NTET)=I1
   13 K=K+1
      IF(K.LE.J) GOTO 11
      J=J+1
      IF(J.LE.N-I-1) GOTO 10
      I=I+1
      IF(I.LT.N) GOTO 9
      AVOL=1.0D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   96 WRITE(6,101)
      GOTO 99
   98 WRITE(6,103)
   99 STOP
  100 FORMAT(' SAMPLING THE 48TH PART OF THE B.Z. OF THE',
     ,' BCC LATTICE'/1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' KXMAX =',D11.4,'  KYMAX =',D11.4,'  KZMAX =',D11.4)
  101 FORMAT(' *** <SETK02> NA IS NOT A POSITIVE INTEGER ***')
  103 FORMAT(' *** <SETK02> A IS NOT A POSITIVE NUMBER ***')
      END
