      SUBROUTINE SETK09(NA,NB,NC,A,B,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN ONE-HEIGTH THE RECIPROCAL CELL FOR A
C     SIDE-CENTERED ORTHORHOMBIC LATTICE WITH PARAMETERS A, B, C
C     TYPE CAN BE A, B OR DEPENDING ON WETHER THE (B,C) , (C,A)
C     OR (A,B) FACES ARE CENTERED (NA, NB OR NC NEGATIVE, RESPECTIVELY
C     SYMMETRY IS D2H
C     EXTERNAL ROUTINE : SETK07
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 TYPE
      REAL AVOL
      DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX),AA(3),NN(3)
      EQUIVALENCE (IVOL,AVOL)
      IF(NA*NB*NC.GE.0) GOTO 97
      IF(NA+NB+NC.EQ.-IABS(NA)-IABS(NB)-IABS(NC)) GOTO 97
      IF(A.LE.0.0D0 .OR. B.LE.0.0D0 .OR. C.LE.0.0D0) GOTO 98
C *** INITIALIZATION
      PI = 3.141592653589793238D0
      NN(1) = IABS(NA)
      NN(2) = IABS(NB)
      NN(3) = IABS(NC)
      AA(1) = A
      AA(2) = B
      AA(3) = C
      IF(NA.LT.0) THEN
         K1 = 2
         K2 = 3
         K3 = 1
         TYPE = 'A'
      ENDIF
      IF(NB.LT.0) THEN
         K1 = 3
         K2 = 1
         K3 = 2
         TYPE = 'B'
      ENDIF
      IF(NC.LT.0) THEN
         K1 = 1
         K2 = 2
         K3 = 3
         TYPE = 'C'
      ENDIF
      IF(NN(K1).EQ.NN(K2)) GOTO 1
      IF(NN(K1).LT.NN(K2)) THEN
         K = K1
         K1 = K2
         K2 = K
      ENDIF
      AA(K2) = AA(K2)/2.0D0
      NN(K2) = 2*NN(K2)
      CALL SETK07(NN(1),NN(2),NN(3),AA(1),AA(2),AA(3),PTK,NPTK,IDEF,
     ,            NTET,NKMAX,NTMAX)
      RETURN
    1 NPTK = (NN(K1)+1)*(2*NN(K2)+1)*(NN(K3)+1)
      IF(NPTK.GT.NKMAX) STOP '*** <SETK09> NPTK EXCEEDS NKMAX ***'
      NTET = 12*NN(K1)*NN(K2)*NN(K3)
      IF(NTET.GT.NTMAX) STOP '*** <SETK09> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      AA(K1) = PI/AA(K1)/NN(K1)
      AA(K2) = PI/AA(K2)/NN(K2)
      AA(K3) = PI/AA(K3)/NN(K3)
      NN(K1) = 2*NN(K1)
      NN(K2) = 2*NN(K2)
      WRITE(6,100) TYPE,NPTK,NTET,NN(1)*AA(1),NN(2)*AA(2),NN(3)*AA(3)
      N1 = NN(K1)
      N3 = NN(K3)
      W = 0.5D0/IABS(NA*NB*NC)
      NPTK = 0
      NP1 = N1+1
      NP3 = N3+1
      JMAX = NP1
      DO 3 I=1,NP1
         DO 2 J=1,JMAX
         DO 2 K=1,NP3
C           NPTK = (I-1)*NP3*(2*NP1-I+2)/2 + (J-1)*NP3 + K
            WK = W
            IF(I.EQ.1) WK = WK/2.0D0
            IF(J.EQ.1) WK = WK/2.0D0
            IF(J.EQ.JMAX) WK = WK/2.0D0
            IF(I.EQ.NP1 .OR. J.EQ.NP1) WK = WK/2.0D0
            IF(K.EQ.1 .OR. K.EQ.NP3) WK = WK/2.0D0
            NPTK = NPTK+1
            PTK(K1,NPTK) = (I-1)*AA(K1)
            PTK(K2,NPTK) = (J-1)*AA(K2)
            PTK(K3,NPTK) = (K-1)*AA(K3)
            PTK(4,NPTK) = WK
    2    CONTINUE
         JMAX = JMAX-1
    3 CONTINUE
C *** DEFINE THE TETRAHEDRA
      NTET = 0
      JMAX = N1
      I7 = 0
      DO 8 I=1,N1
         DO 7 J=1,JMAX
            DO 6 K=1,N3
               I7 = I7+1
               I3 = I7+NP3
               I6 = I7+NP3*(JMAX+1)
               I5 = I6+1
               NTET = NTET+1
               IDEF(1,NTET) = I7
               IDEF(2,NTET) = I3
               IDEF(3,NTET) = I6
               IDEF(4,NTET) = I5
               I8 = I7+1
               NTET = NTET+1
               IDEF(1,NTET) = I7
               IDEF(2,NTET) = I8
               IDEF(3,NTET) = I3
               IDEF(4,NTET) = I5
               I4 = I3+1
               NTET = NTET+1
               IDEF(1,NTET) = I8
               IDEF(2,NTET) = I3
               IDEF(3,NTET) = I4
               IDEF(4,NTET) = I5
               IF(J.EQ.JMAX) GOTO 6
               I2 = I6+NP3
               NTET = NTET+1
               IDEF(1,NTET) = I3
               IDEF(2,NTET) = I6
               IDEF(3,NTET) = I5
               IDEF(4,NTET) = I2
               I1 = I2+1
               NTET = NTET+1
               IDEF(1,NTET) = I3
               IDEF(2,NTET) = I4
               IDEF(3,NTET) = I5
               IDEF(4,NTET) = I1
               NTET = NTET+1
               IDEF(1,NTET) = I3
               IDEF(2,NTET) = I5
               IDEF(3,NTET) = I2
               IDEF(4,NTET) = I1
    6       CONTINUE
            I7 = I7+1
    7    CONTINUE
         I7 = I7+NP3
         JMAX = JMAX-1
    8 CONTINUE
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   97 WRITE(6,101)
      STOP
   98 WRITE(6,102)
   99 STOP
  100 FORMAT(' SAMPLING THE 8TH PART OF A RHOMBUS-BASED PRISM ORIENTED',
     ,' ALONG ',A1/1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' KXMAX =',D11.4,'  KYMAX =',D11.4,'  KZMAX =',D11.4)
  101 FORMAT(' *** <SETK09> NA, NB OR NC MUST BE NEGATIVE ***')
  102 FORMAT(' *** <SETK09> A, B OR C IS NOT POSITIVE ***')
      END
