      SUBROUTINE SETK11(NA,NC,A,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN THE 1/16TH THE RECIPROCAL CELL OF A
C     BODY-CENTRED TETRAGONAL LATTICE WITH PARAMETERS A, C
C     SYMMETRY IS D4H
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX)
      EQUIVALENCE (IVOL,AVOL)
      PI = 3.141592653589793238D0
      IF(NA.LE.0.OR.NC.LE.0) GOTO 97
      IF(A.LE.0.0D0 .OR. C.LE.0.0D0) GOTO 98
      IF(NC.LT.NA) THEN
         CALL SETK08(NA,2*NC,A,C/2.0D0,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
         RETURN
      ENDIF
      NPTK = (NC+1)*(NA+1)**2
      IF(NPTK.GT.NKMAX) STOP '*** <SETK11> NPTK EXCEEDS NKMAX ***'
      NTET = 6*NC*NA**2
      IF(NTET.GT.NTMAX) STOP '*** <SETK11> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      AK = PI/A/NA
      CK = PI/C/NC
      N2 = 2*NA
      IF(NC.EQ.NA) WRITE(6,100) NPTK,NTET,N2*AK,NA*AK,N2*CK
      IF(NC.GT.NA) WRITE(6,101) NPTK,NTET,N2*AK,NA*AK,NC*CK
      W = 1.0D0/(NC*NA**2)
      KMAX = NC
      NPTK = 0
      J = 0
    1 IMAX = N2-2*J
      I = 0
    2 IF(NC.EQ.NA) KMAX = N2-I-J
      K = 0
    3 NPTK = NPTK+1
      WK = W
      IF(K.EQ.0) WK = WK/2.0D0
      IF(K.EQ.KMAX) WK = WK/2.0D0
      IF(I.EQ.0) WK = WK/2.0D0
      IF(J.EQ.0) WK = WK/2.0D0
      IF(I+J.EQ.0) WK = WK/2.0D0
      IF(I.EQ.IMAX) WK = WK/2.0D0
      IF(NC.EQ.NA) THEN
         IF(I.EQ.IMAX .AND. K.EQ.KMAX) WK = W/3.0D0
         IF(I.EQ.0 .AND. K.EQ.KMAX) WK = W/6.0D0
         IF(IMAX.EQ.0 .AND. K.EQ.NA) WK = W/8.0D0
         IF(I.EQ.N2) WK = W/24.0D0
         IF(K.EQ.N2) WK = W/48.0D0
      ELSE
         IF(I.EQ.N2) WK = WK/2.0D0
      ENDIF
      PTK(1,NPTK) = (I+J)*AK
      PTK(2,NPTK) = J*AK
      PTK(3,NPTK) = K*CK
      PTK(4,NPTK) = WK
      K = K+1
      IF(K.LE.KMAX) GOTO 3
      I = I+1
      IF(I.LE.IMAX) GOTO 2
      J = J+1
      IF(J.LE.NA) GOTO 1
C *** DEFINE THE TETRAHEDRA
      NTET=0
      I7 = 0
      J = 0
    5 IMAX = N2-2*J
      I = 0
    6 IF(NC.EQ.NA) KMAX = N2-I-J
      K = 0
    7 I7 = I7+1
      I8 = I7+1
      I6 = I8+KMAX
      I5 = I6+1
      I3 = I7+(NC+1)*IMAX
      IF(NC.EQ.NA) I3 = I3-NA+J
      I4 = I3+1
      I2 = I4+KMAX
      I1 = I2+1
      IF(I.GT.0) GOTO 8
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I2
      IF(NC.EQ.NA .AND. K.EQ.KMAX-1) GOTO 11
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I2
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I5
      IDEF(3,NTET) = I2
      IDEF(4,NTET) = I1
      GOTO 10
    8 NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I4
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I4
      IF(NC.EQ.NA .AND. K.EQ.KMAX-1) THEN
         IF(I.EQ.IMAX-1) GOTO 11
         GOTO 9
      ENDIF
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I4
      IF(I.EQ.IMAX-1) GOTO 10
      NTET = NTET+1
      IDEF(1,NTET) = I6
      IDEF(2,NTET) = I5
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I1
      NTET = NTET+1
      IDEF(1,NTET) = I6
      IDEF(2,NTET) = I4
      IDEF(3,NTET) = I2
      IDEF(4,NTET) = I1
    9 NTET = NTET+1
      IDEF(1,NTET) = I6
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I2
   10 K = K+1
      IF(K.LT.KMAX) GOTO 7
   11 I7 = I7+1
      I = I+1
      IF(I.LT.IMAX) GOTO 6
      I7 = I7+KMAX
      IF(NC.GT.NA) I7 = I7+1
      J = J+1
      IF(J.LT.NA) GOTO 5
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   97 WRITE(6,102)
      GOTO 99
   98 WRITE(6,103)
   99 STOP
  100 FORMAT(' SAMPLING THE 1/16TH PART OF A DODECAHEDRON'/
     .1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' KXMAX =',D11.4,'  KYMAX =',D11.4,'  KZMAX =',D11.4)
  101 FORMAT(' SAMPLING THE 16TH PART OF A 45Œ-ROTATED SQUARE-',
     .'BASED PRISM'/1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' KXMAX =',D11.4,'  KYMAX =',D11.4,'  KZMAX =',D11.4)
  102 FORMAT(' *** <SETK11> NA OR NC IS NOT A POSITIVE INTEGER ***')
  103 FORMAT(' *** <SETK11> A AND C MUST BE POSITIVE ***')
      END
