      SUBROUTINE SETK06(N1,N2,N3,A1,A2,A3,PTK,NPTK,IDEF,NTET,NKMX,NTMX)
C     SET THE K-POINTS IN ONE HALF OF THE RECIPROCAL CELL FOR A
C     TRICLINIC LATTICE WITH DIRECT LATTICE VECTORS A1, A2, A3
C     SYMMETRY IS C1I
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION PTK(4,NKMX),IDEF(5,NTMX),A1(3),A2(3),A3(3),B(3,3),
     ,          NN(3),KK(3)
      EQUIVALENCE (IVOL,AVOL)
      PI = 3.141592653589793238D0
      IF(N1.LE.0 .OR. N2.LE.0 .OR. N3.LE.0) GOTO 97
C *** GENERATE THE RECIPROCAL BASIS VECTORS
      B(1,1) = A2(2)*A3(3)-A2(3)*A3(2)
      B(1,2) = A2(3)*A3(1)-A2(1)*A3(3)
      B(1,3) = A2(1)*A3(2)-A2(2)*A3(1)
      C = A1(1)*B(1,1)+A1(2)*B(1,2)+A1(3)*B(1,3)
      IF(C.EQ.0.0D0) GOTO 98
      C = PI/C
      B(1,1) = B(1,1)*C/N1
      B(1,2) = B(1,2)*C/N1
      B(1,3) = B(1,3)*C/N1
      B(2,1) = (A3(2)*A1(3)-A3(3)*A1(2))*C/N2
      B(2,2) = (A3(3)*A1(1)-A3(1)*A1(3))*C/N2
      B(2,3) = (A3(1)*A1(2)-A3(2)*A1(1))*C/N2
      B(3,1) = (A1(2)*A2(3)-A1(3)*A2(2))*C/N3
      B(3,2) = (A1(3)*A2(1)-A1(1)*A2(3))*C/N3
      B(3,3) = (A1(1)*A2(2)-A1(2)*A2(1))*C/N3
C *** SORT N1,N2,N3 BY DECREASING ORDER
      NN(1) = N1
      NN(2) = N2
      NN(3) = N3
      KK(1) = 1
      KK(2) = 2
      KK(3) = 3
      DO 2 I=1,2
         DO 1 J=I+1,3
            IF(NN(J).LE.NN(I)) GOTO 1
               N = NN(I)
               NN(I) = NN(J)
               NN(J) = N
               K = KK(I)
               KK(I) = KK(J)
               KK(J) = K
    1    CONTINUE
    2 CONTINUE
      NPTK = 4*N1*N2*N3+2*NN(2)*NN(3)+NN(3)+1
      IF(NPTK.GT.NKMX) STOP '*** <SETK06> NPTK EXCEEDS NKMAX ***'
      NTET = 24*N1*N2*N3
      IF(NTET.GT.NTMX) STOP '*** <SETK06> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      K1 = KK(1)
      K2 = KK(2)
      K3 = KK(3)
      NN1 = NN(1)
      NN2 = NN(2)
      NN3 = NN(3)
      WRITE(6,100) NPTK,NTET,
     ,             NN1*DSQRT(B(K1,1)**2+B(K1,2)**2+B(K1,3)**2),K1,
     ,           2*NN2*DSQRT(B(K2,1)**2+B(K2,2)**2+B(K2,3)**2),K2,
     ,           2*NN3*DSQRT(B(K3,1)**2+B(K3,2)**2+B(K3,3)**2),K3
      W = 1.0D0/(4*N1*N2*N3)
      NPTK = 0
      ICODE = 0
      I2MAX = NN2
      I3MAX = NN3
      I1 = -NN1
    3 I2 = -NN2
    4 I3 = -NN3
    5 ICODE = ICODE+1
      IF(I3.EQ.I3MAX) GOTO 6
      WK = W
      IF(I1.EQ.-NN1) WK = W/2.0D0
      IF(I1.EQ.0) THEN
         IF(I2.EQ.-NN2) WK = W/2.0D0
         IF(I2.EQ.0 .AND. I3.EQ.-NN3) WK = W/2.0D0
      ENDIF
      NPTK = NPTK+1
      PTK(1,NPTK) = I1*B(K1,1)+I2*B(K2,1)+I3*B(K3,1)
      PTK(2,NPTK) = I1*B(K1,2)+I2*B(K2,2)+I3*B(K3,2)
      PTK(3,NPTK) = I1*B(K1,3)+I2*B(K2,3)+I3*B(K3,3)
      PTK(4,NPTK) = WK
      IDEF(5,ICODE) = NPTK
      I3 = I3+1
      GOTO 5
    6 IF(I3.EQ.0) GOTO 9
      IDEF(5,ICODE) = NPTK-2*NN3+1
      I2 = I2+1
      IF(I2.LT.I2MAX) GOTO 4
      IF(I2.EQ.0) THEN
         I3MAX = 0
         GOTO 4
      ENDIF
      N12 = 4*NN2*NN3*(I1+NN1)
      DO 8 I3=1,2*NN3
         ICODE = ICODE+1
         N12 = N12+1
         IDEF(5,ICODE) = N12
    8 CONTINUE
      ICODE = ICODE+1
      IDEF(5,ICODE) = N12-2*NN3+1
      I1 = I1+1
      IF(I1.LT.0) GOTO 3
      IF(I1.EQ.0) THEN
         I2MAX = 0
         GOTO 3
      ENDIF
    9 NPTK = NPTK+1
      PTK(1,NPTK) = 0.0D0
      PTK(2,NPTK) = 0.0D0
      PTK(3,NPTK) = 0.0D0
      PTK(4,NPTK) = W/2.0D0
      IDEF(5,ICODE) = NPTK
      N12 = NPTK
      DO 10 I3=1,NN3
         ICODE = ICODE+1
         N12 = N12-1
         IDEF(5,ICODE) = N12
   10 CONTINUE
      DO 12 I2=1,NN2
         ICODE = ICODE+1
         IDEF(5,ICODE) = N12-2*NN3
         DO 11 I3=1,2*NN3
            ICODE = ICODE+1
            N12 = N12-1
            IDEF(5,ICODE) = N12
   11    CONTINUE
   12 CONTINUE
C *** DEFINE THE TETRAHEDRA
      NTET = 0
      N3P1 = 2*NN3+1
      N12 = N3P1*(2*NN2+1)
      DO 13 I1=1,NN1
      DO 13 I2=1,2*NN2
      IND7 = (I1-1)*N12+N3P1*(I2-1)
      DO 13 I3=1,2*NN3
         IND7 = IND7+1
         IND6 = IND7+N12
         IND2 = IND6+N3P1
         IND1 = IND2+1
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND6)
         IDEF(3,NTET) = IDEF(5,IND2)
         IDEF(4,NTET) = IDEF(5,IND1)
         IND8 = IND7+1
         IND5 = IND6+1
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND6)
         IDEF(3,NTET) = IDEF(5,IND5)
         IDEF(4,NTET) = IDEF(5,IND1)
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND8)
         IDEF(3,NTET) = IDEF(5,IND5)
         IDEF(4,NTET) = IDEF(5,IND1)
         IND3 = IND7+N3P1
         IND4 = IND3+1
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND8)
         IDEF(3,NTET) = IDEF(5,IND4)
         IDEF(4,NTET) = IDEF(5,IND1)
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND3)
         IDEF(3,NTET) = IDEF(5,IND4)
         IDEF(4,NTET) = IDEF(5,IND1)
         NTET = NTET+1
         IDEF(1,NTET) = IDEF(5,IND7)
         IDEF(2,NTET) = IDEF(5,IND3)
         IDEF(3,NTET) = IDEF(5,IND2)
         IDEF(4,NTET) = IDEF(5,IND1)
   13 CONTINUE
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   97 WRITE(6,101)
      GOTO 99
   98 WRITE(6,102)
   99 STOP
  100 FORMAT(' SAMPLING ONE HALF THE RECIPROCAL UNIT CELL'/
     .1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' THE EDGES OF THE SAMPLED PARALLELEPIPED HAVE LENGTHS EQUAL TO'/
     .3(1X,D10.4,' ALONG B',I1,3X))
  101 FORMAT(' *** <SETK06> N1, N2 OR N3 IS NOT POSITIVE ***')
  102 FORMAT(' *** <SETK06> VOLUME OF THE PRIMITIVE CELL IS ZERO ***')
      END
