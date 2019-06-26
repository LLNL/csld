      SUBROUTINE SETK05(NR,AR1,AR2,AR3,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN 1/12TH (NR>0) OR 1/6TH (NR<0) THE PRIMITIVE
C     RECIPROCAL CELL FOR A RHOMBOHEDRAL LATTICE WITH DIRECT PRIMITIVE
C     VECTORS AR1, AR2, AR3 DEFINING THE RHOMBOHEDRON
C     SYMMETRY AT GAMMA IS D3D (NR>0) OR C3I (NR<0)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL AVOL
      DIMENSION PTK(4,NKMAX),IDEF(5,NTMAX),AR1(3),AR2(3),AR3(3),B1(3),
     ,          B2(3),B3(3)
      EQUIVALENCE (IVOL,AVOL)
      PI = 3.141592653589793238D0
      N = IABS(NR)
      IF(N.EQ.0) GOTO 97
      N2 = 2*N
      B1(1) = AR2(2)*AR3(3)-AR2(3)*AR3(2)
      B1(2) = AR2(3)*AR3(1)-AR2(1)*AR3(3)
      B1(3) = AR2(1)*AR3(2)-AR2(2)*AR3(1)
      C = AR1(1)*B1(1)+AR1(2)*B1(2)+AR1(3)*B1(3)
      IF(C.EQ.0.0D0) GOTO 98
      C = PI/C/N
      B1(1) = B1(1)*C
      B1(2) = B1(2)*C
      B1(3) = B1(3)*C
      B2(1) = (AR3(2)*AR1(3)-AR3(3)*AR1(2))*C
      B2(2) = (AR3(3)*AR1(1)-AR3(1)*AR1(3))*C
      B2(3) = (AR3(1)*AR1(2)-AR3(2)*AR1(1))*C
      B3(1) = (AR1(2)*AR2(3)-AR1(3)*AR2(2))*C
      B3(2) = (AR1(3)*AR2(1)-AR1(1)*AR2(3))*C
      B3(3) = (AR1(1)*AR2(2)-AR1(2)*AR2(1))*C
      B3(1) = B1(1)+B2(1)+B3(1)
      B3(2) = B1(2)+B2(2)+B3(2)
      B3(3) = B1(3)+B2(3)+B3(3)
      NTET = 4*N**3
      IF(NR.GT.0) THEN
         JFRAC = 1
         NPTK = N*(N+1)*(4*N+11)/6
      ELSE
         JFRAC = 2
         NPTK = 4*N*(N+1)*(N+2)/3-1
      ENDIF
      NTET = JFRAC*NTET
      IF(NPTK.GT.NKMAX) STOP '*** <SETK05> NPTK EXCEEDS NKMAX ***'
      IF(NTET.GT.NTMAX) STOP '*** <SETK05> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      WRITE(6,100) 12/JFRAC,NPTK,NTET,N*B1(1),N*B1(2),N*B1(3)
      W = 12.0D0/JFRAC/N2**3
      NPTK = 0
      I = 0
    1 NI = I*(I+1)*(6*N-2*I+5)/6
      J = 0
      JMAX = I
      IF(NR.GT.0) JMAX = I/2
      KMAX = N2-I
    2 ICODE = NI+J*(KMAX+1)
      J2 = J
      IF(NR.GT.0) J2 = 2*J
      K = 0
    3 ICODE = ICODE+1
      IF(I+J.EQ.0 .AND. K.LT.N) THEN
         IDEF(5,ICODE) = N-K+1
         GOTO 4
      ENDIF
      IF(I.EQ.N2.AND.(J.EQ.0.OR.J.EQ.N2)) THEN
         IDEF(5,ICODE) = N+1
         GOTO 4
      ENDIF
      WK = W
      IF(J.EQ.0) WK = WK/2.0D0
      IF(J2.EQ.I) WK = WK/2.0D0
      IF(K.EQ.0) WK = WK/2.0D0
      IF(K.EQ.KMAX) WK = WK/2.0D0
      IF(K.EQ.0 .AND. J.EQ.0) WK = W/8.0D0
      IF(NR.GT.0) THEN
         IF(J2.EQ.I) THEN
            IF(K.EQ.0) WK = W*0.3125D0
            IF(K.EQ.KMAX) WK = W*0.1875D0
            IF(KMAX.EQ.0) WK = W/8.0D0
         ENDIF
         IF(I+J.EQ.0) THEN
            WK = W/6.0D0
            IF(K.EQ.N .OR. K.EQ.N2) WK = WK/2.0D0
         ENDIF
      ELSE
         IF(K.EQ.KMAX .AND. J2.EQ.I) WK = W/8.0D0
         IF(I+J.EQ.0) THEN
            WK = W/3.0D0
            IF(K.EQ.N .OR. K.EQ.N2) WK = WK/2.0D0
         ENDIF
      ENDIF
      NPTK = NPTK+1
      PTK(1,NPTK) = I*B1(1)+J*B2(1)+(K-N)*B3(1)
      PTK(2,NPTK) = I*B1(2)+J*B2(2)+(K-N)*B3(2)
      PTK(3,NPTK) = I*B1(3)+J*B2(3)+(K-N)*B3(3)
      PTK(4,NPTK) = WK
      IDEF(5,ICODE) = NPTK
    4 K = K+1
      IF(K.LE.KMAX) GOTO 3
      J = J+1
      IF(J.LE.JMAX) GOTO 2
      I = I+1
      IF(I.LE.N2) GOTO 1
C *** DEFINE THE TETRAHEDRA
      NTET = 0
      I = 0
    5 NI = I*(I+1)*(6*N-2*I+5)/6
      J = 0
      JMAX = I
      IF(NR.GT.0) JMAX = I/2
      KMAX = N2-I
    6 IND7 = NI+J*(KMAX+1)
      K = 0
    7 IND7 = IND7+1
      IND8 = IND7+1
      IND3 = IND7+KMAX+1
      IND4 = IND3+1
      IND6 = IND7+(I+1)*(KMAX+1)-J
      IND5 = IND6+1
      IND2 = IND6+KMAX
      IND1 = IND2+1
      IF(J.EQ.JMAX) GOTO 8
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND8)
      IDEF(3,NTET) = IDEF(5,IND4)
      IDEF(4,NTET) = IDEF(5,IND6)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND3)
      IDEF(3,NTET) = IDEF(5,IND4)
      IDEF(4,NTET) = IDEF(5,IND6)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND3)
      IDEF(2,NTET) = IDEF(5,IND4)
      IDEF(3,NTET) = IDEF(5,IND6)
      IDEF(4,NTET) = IDEF(5,IND2)
      IF(K.EQ.KMAX-1) GOTO 14
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND8)
      IDEF(2,NTET) = IDEF(5,IND4)
      IDEF(3,NTET) = IDEF(5,IND6)
      IDEF(4,NTET) = IDEF(5,IND5)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND4)
      IDEF(2,NTET) = IDEF(5,IND6)
      IDEF(3,NTET) = IDEF(5,IND5)
      IDEF(4,NTET) = IDEF(5,IND1)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND4)
      IDEF(2,NTET) = IDEF(5,IND6)
      IDEF(3,NTET) = IDEF(5,IND2)
      IDEF(4,NTET) = IDEF(5,IND1)
      GOTO 13
    8 IF(NR.GT.0 .AND. I.EQ.2*J) GOTO 10
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND8)
      IDEF(3,NTET) = IDEF(5,IND6)
      IDEF(4,NTET) = IDEF(5,IND2)
      IF(K.EQ.KMAX-1) GOTO 14
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND8)
      IDEF(2,NTET) = IDEF(5,IND6)
      IDEF(3,NTET) = IDEF(5,IND5)
      IDEF(4,NTET) = IDEF(5,IND2)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND8)
      IDEF(2,NTET) = IDEF(5,IND5)
      IDEF(3,NTET) = IDEF(5,IND2)
      IDEF(4,NTET) = IDEF(5,IND1)
      GOTO 13
   10 IND10 = IND2+(I+2)*KMAX-J-1
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND6)
      IDEF(3,NTET) = IDEF(5,IND5)
      IDEF(4,NTET) = IDEF(5,IND10)
      IF(K.EQ.KMAX-2) GOTO 11
      IND11 = IND10+1
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND8)
      IDEF(3,NTET) = IDEF(5,IND5)
      IDEF(4,NTET) = IDEF(5,IND11)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND5)
      IDEF(3,NTET) = IDEF(5,IND10)
      IDEF(4,NTET) = IDEF(5,IND11)
      GOTO 13
   11 IND9 = IND7+2
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND7)
      IDEF(2,NTET) = IDEF(5,IND8)
      IDEF(3,NTET) = IDEF(5,IND5)
      IDEF(4,NTET) = IDEF(5,IND10)
      NTET = NTET+1
      IDEF(1,NTET) = IDEF(5,IND8)
      IDEF(2,NTET) = IDEF(5,IND5)
      IDEF(3,NTET) = IDEF(5,IND9)
      IDEF(4,NTET) = IDEF(5,IND10)
      GOTO 14
   13 K = K+1
      IF(K.LT.KMAX) GOTO 7
   14 J = J+1
      IF(J.LE.JMAX) GOTO 6
      I = I+1
      IF(I.LT.N2) GOTO 5
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 ITET=1,NTET
   15 IDEF(5,ITET) = IVOL
      RETURN
   97 WRITE(6,101)
      GOTO 99
   98 WRITE(6,102)
   99 STOP
  100 FORMAT(' SAMPLING THE',I3,'TH PART OF THE RECIPROCAL',
     .' RHOMBOHEDRAL CELL'/1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' POINT "A" IS AT COORDINATES  ',2(D11.4,' , '),D11.4)
  101 FORMAT(' *** <SETK05> NR IS ZERO ***')
  102 FORMAT(' *** <SETK05> RHOMBOHEDRAL CELL HAS ZERO VOLUME ***')
      END
