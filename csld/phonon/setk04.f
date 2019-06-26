      SUBROUTINE SETK04(NH,NC,A1,A2,C,PTK,NPTK,IDEF,NTET,NKMAX,NTMAX)
C     SET THE K-POINTS IN THE IRREDUCTIBLE PART OF THE B.Z. OF THE
C     HEXAGONAL LATTICE WITH PRINCIPAL AXIS C ALONG Z AND 120Œ
C     PRIMITIVE VECTORS A1 AND A2 IN THE (X,Y) PLANE
C     WHEN NC>0, THE HEXAGONAL DIVISION IS ASSUMED : 1/24 (NH>0) OR 1/12
C     (NH<0) OF THE HEXAGONAL PRISM WITH KZ>0 IS CONSIDERED
C     WHEN NC<0, THE TRIGONAL DIVISION IS ASSUMED : 1/12 (NH>0) OR 1/6
C     OF THE HEXAGONAL PRISM IS CONSIDERED (-PI/C <= KZ < +PI/C)
C     SYMMETRIES AT GAMMA POINT ARE D6H (NH>0 , NC>0), C6H (NH<0 , NC>0)
C     D3D (NH>0 , NC<0) OR C3I=S6 (NH<0 , NC<0)
C     CAN BE USED FOR 2D HEXAGONAL LATTICE BY SETTING NC = 0 : THERE
C     ARE NO TERHAHEDRA (NTET=0 IS RETURNED), C IS IGNORED AND K-POINTS
C     ARE DEFINED IN THE 1/12TH (NH > 0) OR 1/6TH (NH < 0) IRREDUCTIBLE
C     PART OF THE 1ST BRILLOUIN ZONE (C6V OR C6 SYMMETRY AT GAMMA).
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION A1(2),A2(2),B1(2),B2(2),PTK(4,NKMAX),IDEF(5,NTMAX)
      EQUIVALENCE (IVOL,AVOL)
      NZ = IABS(NC)
      NAS2 = IABS(NH)
      NA = 2*NAS2
      IF(NA.EQ.0) GOTO 96
      IF(NZ.GT.0 .AND. C.LE.0.0D0) GOTO 97
      IF(NC.LT.0) THEN
         JFRAC = 2
         KMIN = -NZ
         KMAX = NZ-1
      ELSE
         JFRAC = 1
         KMIN = 0
         KMAX = NZ
      ENDIF
      IZ = KMAX-KMIN+1
      IF(NH.GT.0) THEN
         NPTK = IZ*(NAS2+1)**2
         NTET = 6*JFRAC*NZ*NAS2**2
         IFRAC = 24
      ELSE
         NPTK = IZ*(NAS2+1)*(NA+1)
         NTET = 3*JFRAC*NZ*NA**2
         IFRAC = 12
      ENDIF
      IF(NPTK.GT.NKMAX) STOP '*** <SETK04> NPTK EXCEEDS NKMAX ***'
      IF(NTET.GT.NTMAX) STOP '*** <SETK04> NTET EXCEEDS NTMAX ***'
C *** SET THE K-POINTS
      PI = 3.141592653589793238D0
      S = A1(1)*A2(2) - A2(1)*A1(2)
      IF(S.EQ.0.0D0) GOTO 98
      S = 2.0D0*PI/3/NA/S
      B1(1) = (2.0D0*A2(2)+A1(2))*S
      B1(2) =-(2.0D0*A2(1)+A1(1))*S
      B2(1) = (A2(2)-A1(2))*S
      B2(2) =-(A2(1)-A1(1))*S
      IF(NZ.GT.0) THEN
         DK3= PI/C/NZ
         W = (IFRAC/JFRAC)/(6.0D0*NZ*NA**2)
      ELSE
         DK3 = 0.0D0
         W = IFRAC/(3.0D0*NA**2)
         IFRAC = IFRAC/2
      ENDIF
      WRITE(6,100) JFRAC,IFRAC,NPTK,NTET,NA*B2(1),NA*B2(2),NZ*DK3
      NPTK=0
      I = 0
      IMAX = NA
      IF(NH.GT.0) IMAX = NAS2
    1 J = 0
      IF(NH.GT.0) J = I
      JMAX = NA-I
    2 K = KMIN
    3 NPTK = NPTK+1
      WK = W
      IF(K.EQ.NZ) WK = WK/2.0D0
      IF(NC.GT.0 .AND. K.EQ.0) WK = WK/2.0D0
      IF(I.EQ.0) WK = WK/2.0D0
      IF(J.EQ.JMAX) WK = WK/2.0D0
      IF(NH.LT.0) THEN
         IF(J.EQ.0) WK = WK/2.0D0
         IF(I+J.EQ.0) WK = WK/1.5D0
      ELSE
         IF(J.EQ.0) WK = WK/3.0D0
         IF(J.EQ.I) WK = WK/2.0D0
      ENDIF
      IF(I.EQ.NA .OR. J.EQ.NA) WK = WK/1.5D0
      PTK(1,NPTK) = I*B1(1)+J*B2(1)
      PTK(2,NPTK) = I*B1(2)+J*B2(2)
      PTK(3,NPTK) = K*DK3
      PTK(4,NPTK) = WK
      K = K+1
      IF(K.LE.KMAX) GOTO 3
      J = J+1
      IF(J.LE.JMAX) GOTO 2
      I = I+1
      IF(I.LE.IMAX) GOTO 1
      IF(NC.EQ.0) RETURN
C *** DEFINE THE TETRAHEDRA
      I7 = 0
      NTET=0
      I = 0
    5 JMAX = NA-I
      IF(NH.GT.0) THEN
         J = I
         II = (NA-2*I)*IZ
      ELSE
         J = 0
         II = (NA+1-I)*IZ
      ENDIF
    6 K = KMIN
    7 I7 = I7+1
      I8 = I7+1
      I6 = I7+II
      I5 = I6+1
      I3 = I7+IZ
      I4 = I3+1
      I2 = I3+II
      I1 = I2+1
      IF(NC.LT.0 .AND. K.EQ.KMAX) THEN
         I8 = I8-2*NZ
         I5 = I5-2*NZ
         I4 = I4-2*NZ
         I1 = I1-2*NZ
      ENDIF
      IF(NH.LT.0 .OR. J.GT.I) GOTO 8
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I2
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I2
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I4
      IDEF(3,NTET) = I2
      IDEF(4,NTET) = I1
      GOTO 9
    8 NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I8
      IDEF(3,NTET) = I3
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I7
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I6
      IDEF(4,NTET) = I5
      NTET = NTET+1
      IDEF(1,NTET) = I8
      IDEF(2,NTET) = I3
      IDEF(3,NTET) = I4
      IDEF(4,NTET) = I5
      IF(J.EQ.JMAX-1) GOTO 9
      NTET = NTET+1
      IDEF(1,NTET) = I3
      IDEF(2,NTET) = I4
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I1
      NTET = NTET+1
      IDEF(1,NTET) = I3
      IDEF(2,NTET) = I6
      IDEF(3,NTET) = I5
      IDEF(4,NTET) = I2
      NTET = NTET+1
      IDEF(1,NTET) = I3
      IDEF(2,NTET) = I5
      IDEF(3,NTET) = I2
      IDEF(4,NTET) = I1
    9 K = K+1
      IF(K.LT.NZ) GOTO 7
      IF(NC.GT.0) I7 = I7+1
      J = J+1
      IF(J.LT.JMAX) GOTO 6
      I7 = I7+IZ
      I = I+1
      IF(I.LT.IMAX) GOTO 5
      AVOL=1.0D0/DFLOAT(NTET)
      DO 10 IT=1,NTET
   10 IDEF(5,IT)=IVOL
      RETURN
   96 WRITE(6,101)
      GOTO 99
   97 WRITE(6,102)
      GOTO 99
   98 WRITE(6,103)
   99 STOP
  100 FORMAT(' SAMPLING THE',I2,'/',I2,'TH PART OF THE B.Z. OF THE',
     .' HEXAGONAL LATTICE'/1X,I5,' K-POINTS',I7,' TETRAHEDRA'/
     .' POINT "H" IS AT COORDINATES  ',2(D11.4,' , '),D11.4)
  101 FORMAT(' *** <SETK04> NH IS ZERO ***')
  102 FORMAT(' *** <SETK04> C IS NOT POSITIVE ***')
  103 FORMAT(' *** <SETK04> PRIMITIVE CELL HAS ZERO VOLUME ***')
      END
