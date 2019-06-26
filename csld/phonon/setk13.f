      SUBROUTINE SETK13(N1,N2,N3,A1,A2,A3,PTK,NPTK,IDEF,NTET,NKMX,NTMX)
C     SET THE K-POINTS IN ONE FOURTH THE RECIPROCAL CELL FOR A
C     SIMPLE MONOCLINIC LATTICE WITH DIRECT-SPACE TRANSLATION VECTORS
C     A1, A2 AND A3, WITH A3 ALONG THE BINARY AXIS, A1 AND A2
C     PERPENDICULAR TO IT.
C     SYMMETRY IS C2H
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION PTK(4,NKMX),IDEF(5,NTMX),A1(3),A2(3),A3(3),B(3,3)
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
      NPTK = (2*N1*N2+MIN0(N1,N2)+1)*(N3+1)
      IF(NPTK.GT.NKMX) STOP '*** <SETK13> NPTK EXCEEDS NKMAX ***'
      NTET = 12*N1*N2*N3
      IF(NTET.GT.NTMX) STOP '*** <SETK13> NTET EXCEEDS NTMAX ***'
      IF(N2.LE.N1) THEN
         WRITE(6,100) NPTK,NTET,
     ,                N1*DSQRT(B(1,1)**2+B(1,2)**2+B(1,3)**2),
     ,              2*N2*DSQRT(B(2,1)**2+B(2,2)**2+B(2,3)**2),
     ,                N3*DSQRT(B(3,1)**2+B(3,2)**2+B(3,3)**2)
         NN1 = N1
         NN2 = N2
      ELSE
         WRITE(6,100) NPTK,NTET,
     ,              2*N1*DSQRT(B(1,1)**2+B(1,2)**2+B(1,3)**2),
     ,                N2*DSQRT(B(2,1)**2+B(2,2)**2+B(2,3)**2),
     ,                N3*DSQRT(B(3,1)**2+B(3,2)**2+B(3,3)**2)
         NN1 = N2
         NN2 = N1
         DO 1 L=1,3
            C = B(1,L)
            B(1,L) = B(2,L)
            B(2,L) = C
    1    CONTINUE
      ENDIF
C *** SET THE K-POINTS
      N3P1 = N3+1
      N12 = (2*NN2-1)*N3P1
      ICODE = 0
      NPTK = 0
      W = 1.0D0/(2*N1*N2*N3)
      I = 0
    2 J = -NN2
    3 K = 0
    4 ICODE = ICODE+1
      IF(I.GT.0) GOTO 5
      IF(J.LT.0) THEN
         IDEF(5,ICODE) = K+1-N3P1*J
         GOTO 7
      ENDIF
    5 IF(J.EQ.-NN2) THEN
         IDEF(5,ICODE) = NPTK+N12+K+1
         GOTO 7
      ENDIF
      NPTK = NPTK+1
      IDEF(5,ICODE) = NPTK
      WK = W
      IF(I.EQ.0 .AND. (J.EQ.0.OR.J.EQ.NN2)) WK = WK/2.0D0
      IF(I.EQ.NN1) WK = WK/2.0D0
      IF(K.EQ.0 .OR. K.EQ.N3) WK = WK/2.0D0
      DO 6 L=1,3
         PTK(L,NPTK) = I*B(1,L)+J*B(2,L)+K*B(3,L)
    6 CONTINUE
      PTK(4,NPTK) = WK
    7 K = K+1
      IF(K.LE.N3) GOTO 4
      J = J+1
      IF(J.LE.NN2) GOTO 3
      I = I+1
      IF(I.LE.NN1) GOTO 2
C *** DEFINE THE TETRAHEDRA
      N12 = (2*NN2+1)*N3P1
      NTET = 0
      IND7 = 0
      I = 0
   10 J = -NN2
   11 K = 0
   12 IND7 = IND7+1
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
      K = K+1
      IF(K.LT.N3) GOTO 12
      IND7 = IND7+1
      J = J+1
      IF(J.LT.NN2) GOTO 11
      IND7 = IND7+N3P1
      I = I+1
      IF(I.LT.NN1) GOTO 10
      AVOL=1.D0/DFLOAT(NTET)
      DO 15 IT=1,NTET
   15 IDEF(5,IT)=IVOL
      RETURN
   97 WRITE(6,101)
      GOTO 99
   98 WRITE(6,102)
   99 STOP
  100 FORMAT(' SAMPLING ONE FOURTH THE MONOCLINIC RECIPROCAL CELL'/
     ,1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' THE EDGES OF THE SAMPLED PARALLELEPIPED HAVE LENGTHS EQUAL TO'/
     .2D11.4,' IN THE BASAL PLANE AND',D11.4,' ALONG THE BINARY AXIS')
  101 FORMAT(' *** <SETK13> N1, N2 OR N3 IS NOT POSITIVE ***')
  102 FORMAT(' *** <SETK13> VOLUME OF THE PRIMITIVE CELL IS ZERO ***')
      END
