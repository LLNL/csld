      SUBROUTINE PDSTET(CHAR,IDIMC,JDIMC,NPDS,ENER,IDIME,NBAND,IDEF,NTET
     ,                 ,XE,NE,YT,ZT,YP,ZP,IDIMP)
C     COMPUTE PARTIAL DENSITIES AND NUMBERS OF STATES USING THE
C     TETRAHEDRON METHOD ON A MESH OF NE ENERGIES STORED IN XE.
C     ENER(NB,IK) : EIGENENERGY OF THE NTH BAND AT THE K-POINT # IK
C     CHAR(L,NB,IK) : TABLE OF CHARACTERS, L = 1, 2, ... NPDS
C     NBAND : NUMBER OF BANDS
C     NTET  : NUMBER OF TETRAHEDRA (AS GIVEN BY SETK**)
C     ITET  : DEFINES THE TETRAHEDRA (AS GIVEN BY SETK**)
C     IDIMC : FIRST DIMENSION OF CHAR (>= NPDS)
C     JDIMC : SECOND DIMENSION OF CHAR (>= NBAND)
C     IDIME : FIRST DIMENSION OF ENER (>= NBAND)
C     IDIMP : FIRST DIMENSION OF YP AND ZP (>=NPDS)
C     YT(IE)   : TOTAL DENSITY OF STATES AT X(IE), IE = 1 ... NE
C     ZT(IE)   : TOTAL NUMBER OF STATES AT X(IE) , IE = 1 ... NE
C     YP(L,IE) : LTH PARTIAL DENSITY OF STATES AT X(IE), L = 1 ... NPDS
C     ZP(L,IE) : LTH PARTIAL NUMBER OF STATES AT X(IE) , L = 1 ... NPDS
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 AVOL
      DIMENSION CHAR(IDIMC,JDIMC,1),ENER(IDIME,1),XE(NE),YT(NE),ZT(NE),
     ,          YP(IDIMP,NE),ZP(IDIMP,NE),IDEF(5,1)
      DIMENSION C(4),P(4),S(4),IND(4)
      EQUIVALENCE (IVOL,AVOL),(S(1),E1),(S(2),E2),(S(3),E3),(S(4),E4),
     ,            (IND(1),I1),(IND(2),I2),(IND(3),I3),(IND(4),I4)
      DATA EPS/1.0D-17/
C
C *** INITIALIZATION
C
      IF(NE.LE.0) RETURN
      DO 2 IE=1,NE
         YT(IE) = 0.0D0
         ZT(IE) = 0.0D0
         DO 1 IP=1,NPDS
            YP(IP,IE) = 0.0D0
            ZP(IP,IE) = 0.0D0
    1    CONTINUE
    2 CONTINUE
      IF(NBAND.LE.0) RETURN
C
C *** LOOP OVER THE TETRAHEDRA
C
      DO 99 ITET=1,NTET
      I1 = IDEF(1,ITET)
      I2 = IDEF(2,ITET)
      I3 = IDEF(3,ITET)
      I4 = IDEF(4,ITET)
      IVOL = IDEF(5,ITET)
C
C *** LOOP OVER THE BANDS
C
      DO 99 NB=1,NBAND
C *** SITES I1 ... I4 ARE SORTED BY DECREASING ENERGIES
      E1 = ENER(NB,I1)
      E2 = ENER(NB,I2)
      E3 = ENER(NB,I3)
      E4 = ENER(NB,I4)
      DO 4 I=1,3
         SS = S(I)
         II = IND(I)
         J = I
    3    J = J+1
         IF(J.GT.4) GOTO 4
         IF(SS.GE.S(J)) GOTO 3
         S(I) = S(J)
         S(J) = SS
         SS = S(I)
         IND(I) = IND(J)
         IND(J) = II
         II = IND(I)
         GOTO 3
    4 CONTINUE
C *** COMPUTE THE ENERGY DIFFERENCE. WHEN POSITIVE, E1-E4 = ENERGY UNIT
      UNIT = 1.0D0
      IF(E1.GT.E4) UNIT = E1-E4
      E12 = (E1-E2)/UNIT
      E13 = (E1-E3)/UNIT
      E23 = (E2-E3)/UNIT
      E24 = (E2-E4)/UNIT
      E34 = (E3-E4)/UNIT
      FACC = AVOL/UNIT
      FACP = 0.25D0*AVOL
C *** COMPUTE THE WEIGTH FACTORS ON THE ENERGY MESH
      DO 98 IE=1,NE
      EE = XE(IE)
      IF(EE.LT.E4) GOTO 98
      IF(EE.LE.E1) GOTO 5
      P(1) = FACP
      P(2) = FACP
      P(3) = FACP
      P(4) = FACP
      GOTO 96
    5 EE1 = (EE-E1)/UNIT
      IF(DABS(EE1).LT.EPS) EE1 = 0.0D0
      EE2 = (EE-E2)/UNIT
      IF(DABS(EE2).LT.EPS) EE2 = 0.0D0
      EE3 = (EE-E3)/UNIT
      IF(DABS(EE3).LT.EPS) EE3 = 0.0D0
      EE4 = (EE-E4)/UNIT
      IF(DABS(EE4).LT.EPS) EE4 = 0.0D0
      IF(EE.GT.E3) GOTO 10
      IF(E3.EQ.E4) GOTO 6
C     E4 <= EE <= E3
      SURF = (EE4/E34)*(EE4/E24)
      VOLU = EE4*SURF
      C(1) = VOLU
      C(2) = VOLU/E24
      C(3) = VOLU/E34
      C(4) = -SURF*(EE1+EE2/E24+EE3/E34)
      P(1) = C(1)*EE4
      P(2) = C(2)*EE4
      P(3) = C(3)*EE4
      P(4) = VOLU+C(4)*EE4
      GOTO 94
    6 IF(E2.GT.E3) GOTO 98
      IF(E1.EQ.E2) GOTO 7
C     E4 = E3 = E2 = EE < E1
      C(1) = 0.0D0
      C(2) = 1.0D0
      C(3) = 1.0D0
      C(4) = 1.0D0
      P(1) = 0.0D0
      P(2) = 0.0D0
      P(3) = 0.0D0
      P(4) = 0.0D0
      GOTO 94
C     E4 = E3 = E2 = E1 = EE
    7 C(1) = 1.0D+15
      C(2) = 1.0D+15
      C(3) = 1.0D+15
      C(4) = 1.0D+15
      P(1) = 2.0D0
      P(2) = 2.0D0
      P(3) = 2.0D0
      P(4) = 2.0D0
      GOTO 94
   10 IF(EE.GT.E2) GOTO 15
C     E3 < EE <= E2
      SABC = -0.5D0*(EE3/E13)*(EE2/E23)/E24
      SCDA = -0.5D0*(EE1/E13)*(EE4/E24)
      SDAB = -0.5D0*(EE1/E13)*(EE3/E23)
      SBCD = -0.5D0*(EE2/E23)*(EE4/E24)
      VOLU = EE4/E24*EE4
      C(1) = (EE3/E13)*(SABC+SCDA+SDAB) + EE4*(SBCD+SCDA+SDAB)
      C(2) = (EE3/E23)*(SABC+SBCD+SDAB) + (EE4/E24)*(SABC+SBCD+SCDA)
      C(3) =-(EE2/E23)*(SABC+SBCD+SDAB) - (EE1/E13)*(SABC+SCDA+SDAB)
      C(4) =-(EE2/E24)*(SABC+SBCD+SCDA) - EE1*(SBCD+SCDA+SDAB)
      P(1) = EE4*VOLU + EE3*C(1)
      P(2) = (EE4/E24)*VOLU + EE3*C(2)
      P(3) = VOLU + EE3*(SABC+SBCD+SCDA+SDAB+C(3))
      P(4) = VOLU*(1.0D0-EE1-EE2/E24) + EE3*C(4)
      GOTO 94
C     E2 < EE <= E1
   15 SURF = (EE1/E12)*(EE1/E13)
      VOLU =-SURF*EE1
      C(1) = SURF*(EE2/E12+EE3/E13+EE4)
      C(2) = VOLU/E12
      C(3) = VOLU/E13
      C(4) = VOLU
      P(1) = 1.0D0-VOLU+EE1*C(1)
      P(2) = 1.0D0+EE1*C(2)
      P(3) = 1.0D0+EE1*C(3)
      P(4) = 1.0D0+EE1*C(4)
C *** RENORMALIZE THE WEIGHT FACTORS
   94 C(1) = FACC*C(1)
      C(2) = FACC*C(2)
      C(3) = FACC*C(3)
      C(4) = FACC*C(4)
      P(1) = FACP*P(1)
      P(2) = FACP*P(2)
      P(3) = FACP*P(3)
      P(4) = FACP*P(4)
C *** TETRAHEDRON CONTRIBUTION TO THE DENSITIES OF STATES
      YT(IE) = C(1)+C(2)+C(3)+C(4)+YT(IE)
      DO 95 IP=1,NPDS
         YP(IP,IE) = C(1)*CHAR(IP,NB,I1)+C(2)*CHAR(IP,NB,I2) +
     +               C(3)*CHAR(IP,NB,I3)+C(4)*CHAR(IP,NB,I4) + YP(IP,IE)
   95 CONTINUE
C *** TETRAHEDRON CONTRIBUTION TO THE NUMBERS OF STATES
   96 ZT(IE) = P(1)+P(2)+P(3)+P(4)+ZT(IE)
      DO 97 IP=1,NPDS
         ZP(IP,IE) = P(1)*CHAR(IP,NB,I1)+P(2)*CHAR(IP,NB,I2) +
     +               P(3)*CHAR(IP,NB,I3)+P(4)*CHAR(IP,NB,I4) + ZP(IP,IE)
   97 CONTINUE
   98 CONTINUE
   99 CONTINUE
      RETURN
      END
