      SUBROUTINE SETK10(NA,NB,NC,A,B,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN ONE-HEIGTH THE RECIPROCAL CELL FOR A
C     BODY-CENTERED ORTHORHOMBIC LATTICE WITH PARAMETERS A, B, C
C     SYMMETRY IS D2H
C     EXTERNAL SOURCES : SETK07, SETK09
      IMPLICIT REAL*8(A-H,O-Z)
      REAL AVOL
      DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX)
      EQUIVALENCE (IVOL,AVOL)
      IF(NA.LE.0 .OR. NB.LE.0 .OR. NC.LE.0) GOTO 97
      IF(A.LE.0.0D0 .OR. B.LE.0.0D0 .OR. C.LE.0.0D0) GOTO 98
      IF(NA.EQ.NB .AND. NB.EQ.NC) GOTO 1
C *** PRALLELIPIPED (SETK07) OR TRIANGULAR PRISM (SETK09)
      IF(NC.LE.NA .AND. NC.LE.NB) THEN
        IF(NA.LT.NC .AND. NB.LT.NC) THEN
          CALL SETK07(NA,NB,2*NC,A,B,C/2.0D0,PTK,NPTK,IDEF,NTET,NKMAX,
     ,               NTMAX)
        ELSE
          IF(NB.EQ.NC) THEN
            CALL SETK09(-NA,NB,NC,A,B,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
          ELSE
            CALL SETK09(NA,-NB,NC,A,B,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
          ENDIF
        ENDIF
      ELSE
        IF(NB.LE.NA) THEN
          IF(NB.EQ.NA) THEN
            CALL SETK09(NA,NB,-NC,A,B,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
          ELSE
            CALL SETK07(NA,2*NB,NC,A,B/2.0D0,C,PTK,NPTK,IDEF,NTET,NKMAX,
     ,               NTMAX)
          ENDIF
        ELSE
          CALL SETK07(2*NA,NB,NC,A/2.0D0,B,C,PTK,NPTK,IDEF,NTET,NKMAX,
     ,               NTMAX)
        ENDIF
      ENDIF
      RETURN
C *** INITIALIZATION
    1 PI = 3.141592653589793238D0
      NPTK = (NA+1)*(4*NA**2+5*NA+2)/2
      IF(NPTK.GT.NKMAX) STOP '*** <SETK10> NPTK EXCEEDS NKMAX ***'
      NTET = 12*NA**3
      IF(NTET.GT.NTMAX) STOP '*** <SETK10> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      AK = PI/A/NA
      BK = PI/B/NA
      CK = PI/C/NA
      N = 2*NA
      WRITE(6,100) NPTK,NTET,N*AK,N*BK,N*CK
      W = 0.5D0/(NA*NB*NC)
      NM1 = N-1
      NP1 = N+1
      NPTK = 0
      I = 0
    2 JMAX = N-I
      J = 0
    3 KMAX = MIN0(N-I,N-J)
      K = 0
      IND7 = NP1*I+J+1
      IDEF(5,IND7) = NPTK+1
    4 NPTK = NPTK+1
      WK = W
      IF(I.EQ.0) WK = WK/2.0D0
      IF(J.EQ.0) WK = WK/2.0D0
      IF(K.EQ.0) WK = WK/2.0D0
      IF(J.EQ.JMAX) WK = WK/2.0D0
      IF(K.EQ.KMAX) WK = WK/2.0D0
      IF(I.EQ.J .AND. K.EQ.KMAX) WK = W/3.0D0
      IF(J.EQ.JMAX .AND. K.EQ.KMAX) WK = W/3.0D0
      IF(I.EQ.N .OR. J.EQ.N .OR. K.EQ.N) WK = W/24.0D0
      IF(I+J+K.EQ.3*NA) WK = W/4.0D0
      PTK(1,NPTK) = I*AK
      PTK(2,NPTK) = J*BK
      PTK(3,NPTK) = K*CK
      PTK(4,NPTK) = WK
      K = K+1
      IF(K.LE.KMAX) GOTO 4
      J = J+1
      IF(J.LE.JMAX) GOTO 3
      I = I+1
      IF(I.LE.N) GOTO 2
C *** DEFINE THE TETRAHEDRA
      NTET = 0
      I = 0
    6 JMAX = NM1-I
      J = 0
    7 IND7 = NP1*I+J+1
      I7 = IDEF(5,IND7)-1
      IJ = I+J
      KMAX = MIN0(NM1-I,NM1-J)
      K = 0
    8 JK = J+K
      KI = K+I
      I7 = I7+1
      I3 = IDEF(5,IND7+1)+K
      I6 = IDEF(5,IND7+NP1)+K
      I2 = IDEF(5,IND7+NP1+1)+K
      I8 = I7+1
      I5 = I6+1
      I4 = I3+1
      I1 = I2+1
      IF(IJ.EQ.NM1) GOTO 9
      IF(JK.EQ.NM1) GOTO 10
      IF(KI.EQ.NM1) GOTO 11
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I3
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I2
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
      GOTO 12
    9 IF(KI.EQ.NM1) GOTO 11
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I5
      IF(JK.EQ.NM1) GOTO 12
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I5
      GOTO 12
   10 NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I2
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I2
      IF(KI.EQ.NM1) GOTO 12
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I2
      GOTO 12
   11 NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I6
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I6
      IF(IJ.EQ.NM1) GOTO 12
      NTET = NTET+1
      IDEF(1,NTET) = I3
      IDEF(2,NTET) = I4
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I2
   12 K = K+1
      IF(K.LE.KMAX) GOTO 8
      J = J+1
      IF(J.LE.JMAX) GOTO 7
      I = I+1
      IF(I.LE.NM1) GOTO 6
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   97 WRITE(6,101)
      STOP
   98 WRITE(6,102)
   99 STOP
  100 FORMAT(' SAMPLING THE 8TH PART OF A DODECAHEDRON'/
     .1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' KXMAX =',D11.4,'  KYMAX =',D11.4,'  KZMAX =',D11.4)
  101 FORMAT(' *** <SETK10> NA, NB OR NC IS NOT POSITIVE ***')
  102 FORMAT(' *** <SETK10> A, B OR C IS NOT POSITIVE ***')
      END
