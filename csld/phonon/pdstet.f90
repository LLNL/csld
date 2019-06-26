!original version
! http://www.unamur.be/sciences/physique/administration/tetpack/pdstet.f
!*==PDSTET.spg  processed by SPAG 6.72Dc at 03:40 on 26 Sep 2016
SUBROUTINE PDSTET(Char,Idimc,Jdimc,Npds,Ener,Idime,Nband,Idef,    &
     & Ntet,Xe,Ne,Yt,Zt,Yp,Zp,Idimp)
  !     COMPUTE PARTIAL DENSITIES AND NUMBERS OF STATES USING THE
  !     TETRAHEDRON METHOD ON A MESH OF NE ENERGIES STORED IN XE.
  !     ENER(NB,IK) : EIGENENERGY OF THE NTH BAND AT THE K-POINT # IK
  !     CHAR(L,NB,IK) : TABLE OF CHARACTERS, L = 1, 2, ... NPDS
  !     NBAND : NUMBER OF BANDS
  !     NTET  : NUMBER OF TETRAHEDRA (AS GIVEN BY SETK**)
  !     ITET  : DEFINES THE TETRAHEDRA (AS GIVEN BY SETK**)
  !     IDIMC : FIRST DIMENSION OF CHAR (>= NPDS)
  !     JDIMC : SECOND DIMENSION OF CHAR (>= NBAND)
  !     IDIME : FIRST DIMENSION OF ENER (>= NBAND)
  !     IDIMP : FIRST DIMENSION OF YP AND ZP (>=NPDS)
  !     YT(IE)   : TOTAL DENSITY OF STATES AT X(IE), IE = 1 ... NE
  !     ZT(IE)   : TOTAL NUMBER OF STATES AT X(IE) , IE = 1 ... NE
  !     YP(L,IE) : LTH PARTIAL DENSITY OF STATES AT X(IE), L = 1 ... NPDS
  !     ZP(L,IE) : LTH PARTIAL NUMBER OF STATES AT X(IE) , L = 1 ... NPDS
  IMPLICIT NONE
  !*--PDSTET20
  !*** Start of declarations inserted by SPAG
  integer, intent(in)   :: Idimc, Jdimc, Npds, Nband, Ntet, Ne, Idime , Idimp
  real(8), intent(in)  :: Char(Idimc,Jdimc,1), Ener(Idime,1), Xe(Ne) 
  integer, intent(in)   :: Idef(5, Ntet)
  real(8), intent(out)   :: Yp(Idimp,Ne) , Yt(Ne) , Zp(Idimp,Ne) , Zt(Ne)

  ! local var
  real(8) c(4),  e1 , e12 , e13 , e2 , e23 , e24 , e3 , e34 ,    &
       & e4 , ee , ee1 , ee2 , ee3 , ee4 , eps , facc , facp
  INTEGER i , i1 , i2 , i3 , i4 , &
       & ie , ii , ind(4) , ip , itet , ivol , j ,   nb 
  real(8) p(4) , s(4) , sabc , sbcd , scda , sdab , ss , surf , unit , volu 
  !*** End of declarations inserted by SPAG
  real(4) avol
  EQUIVALENCE (ivol,avol)
  EQUIVALENCE (s(1),e1)
  EQUIVALENCE (s(2),e2)
  EQUIVALENCE (s(3),e3)
  EQUIVALENCE (s(4),e4)
  EQUIVALENCE (ind(1),i1)
  EQUIVALENCE (ind(2),i2)
  EQUIVALENCE (ind(3),i3)
  EQUIVALENCE (ind(4),i4)
  DATA eps/1.0D-17/
  !
  ! *** INITIALIZATION
  !
  IF ( Ne<=0 ) RETURN
  DO ie = 1 , Ne
     Yt(ie) = 0.0D0
     Zt(ie) = 0.0D0
     DO ip = 1 , Npds
        Yp(ip,ie) = 0.0D0
        Zp(ip,ie) = 0.0D0
     ENDDO
  ENDDO
  IF ( Nband<=0 ) RETURN
  !
  ! *** LOOP OVER THE TETRAHEDRA
  !
  DO itet = 1 , Ntet
     i1 = Idef(1,itet)
     i2 = Idef(2,itet)
     i3 = Idef(3,itet)
     i4 = Idef(4,itet)
     ivol = Idef(5,itet)
     !
     ! *** LOOP OVER THE BANDS
     !
     DO nb = 1 , Nband
        ! *** SITES I1 ... I4 ARE SORTED BY DECREASING ENERGIES
        e1 = Ener(nb,i1)
        e2 = Ener(nb,i2)
        e3 = Ener(nb,i3)
        e4 = Ener(nb,i4)
        DO i = 1 , 3
           ss = s(i)
           ii = ind(i)
           do j = i + 1, 4
              IF ( ss<s(j) ) THEN
                 s(i) = s(j)
                 s(j) = ss
                 ss = s(i)
                 ind(i) = ind(j)
                 ind(j) = ii
                 ii = ind(i)
              ENDIF
           end do
        ENDDO
        ! *** COMPUTE THE ENERGY DIFFERENCE. WHEN POSITIVE, E1-E4 = ENERGY UNIT
        unit = 1.0D0
        IF ( e1>e4 ) unit = e1 - e4
        e12 = (e1-e2)/unit
        e13 = (e1-e3)/unit
        e23 = (e2-e3)/unit
        e24 = (e2-e4)/unit
        e34 = (e3-e4)/unit
        facc = avol/unit
        facp = 0.25D0*avol
        ! *** COMPUTE THE WEIGTH FACTORS ON THE ENERGY MESH
        DO ie = 1 , Ne
           ee = Xe(ie)
           IF ( ee>=e4 ) THEN
              IF ( ee<=e1 ) THEN
                 ee1 = (ee-e1)/unit
                 IF ( DABS(ee1)<eps ) ee1 = 0.0D0
                 ee2 = (ee-e2)/unit
                 IF ( DABS(ee2)<eps ) ee2 = 0.0D0
                 ee3 = (ee-e3)/unit
                 IF ( DABS(ee3)<eps ) ee3 = 0.0D0
                 ee4 = (ee-e4)/unit
                 IF ( DABS(ee4)<eps ) ee4 = 0.0D0
                 IF ( ee>e3 ) THEN
                    IF ( ee>e2 ) THEN
                       !     E2 < EE <= E1
                       surf = (ee1/e12)*(ee1/e13)
                       volu = -surf*ee1
                       c(1) = surf*(ee2/e12+ee3/e13+ee4)
                       c(2) = volu/e12
                       c(3) = volu/e13
                       c(4) = volu
                       p(1) = 1.0D0 - volu + ee1*c(1)
                       p(2) = 1.0D0 + ee1*c(2)
                       p(3) = 1.0D0 + ee1*c(3)
                       p(4) = 1.0D0 + ee1*c(4)
                    ELSE
                       !     E3 < EE <= E2
                       sabc = -0.5D0*(ee3/e13)*(ee2/e23)/e24
                       scda = -0.5D0*(ee1/e13)*(ee4/e24)
                       sdab = -0.5D0*(ee1/e13)*(ee3/e23)
                       sbcd = -0.5D0*(ee2/e23)*(ee4/e24)
                       volu = ee4/e24*ee4
                       c(1) = (ee3/e13)*(sabc+scda+sdab)            &
                            & + ee4*(sbcd+scda+sdab)
                       c(2) = (ee3/e23)*(sabc+sbcd+sdab) + (ee4/e24)&
                            & *(sabc+sbcd+scda)
                       c(3) = -(ee2/e23)*(sabc+sbcd+sdab)           &
                            & - (ee1/e13)*(sabc+scda+sdab)
                       c(4) = -(ee2/e24)*(sabc+sbcd+scda)           &
                            & - ee1*(sbcd+scda+sdab)
                       p(1) = ee4*volu + ee3*c(1)
                       p(2) = (ee4/e24)*volu + ee3*c(2)
                       p(3) = volu + ee3*(sabc+sbcd+scda+sdab+c(3))
                       p(4) = volu*(1.0D0-ee1-ee2/e24) + ee3*c(4)
                    ENDIF
                 ELSEIF ( e3==e4 ) THEN
                    IF ( e2>e3 ) GOTO 20
                    IF ( e1==e2 ) THEN
                       !     E4 = E3 = E2 = E1 = EE
                       c(1) = 1.0D+15
                       c(2) = 1.0D+15
                       c(3) = 1.0D+15
                       c(4) = 1.0D+15
                       p(1) = 2.0D0
                       p(2) = 2.0D0
                       p(3) = 2.0D0
                       p(4) = 2.0D0
                    ELSE
                       !     E4 = E3 = E2 = EE < E1
                       c(1) = 0.0D0
                       c(2) = 1.0D0
                       c(3) = 1.0D0
                       c(4) = 1.0D0
                       p(1) = 0.0D0
                       p(2) = 0.0D0
                       p(3) = 0.0D0
                       p(4) = 0.0D0
                    ENDIF
                 ELSE
                    !     E4 <= EE <= E3
                    surf = (ee4/e34)*(ee4/e24)
                    volu = ee4*surf
                    c(1) = volu
                    c(2) = volu/e24
                    c(3) = volu/e34
                    c(4) = -surf*(ee1+ee2/e24+ee3/e34)
                    p(1) = c(1)*ee4
                    p(2) = c(2)*ee4
                    p(3) = c(3)*ee4
                    p(4) = volu + c(4)*ee4
                 ENDIF
                 ! *** RENORMALIZE THE WEIGHT FACTORS
                 c(1) = facc*c(1)
                 c(2) = facc*c(2)
                 c(3) = facc*c(3)
                 c(4) = facc*c(4)
                 p(1) = facp*p(1)
                 p(2) = facp*p(2)
                 p(3) = facp*p(3)
                 p(4) = facp*p(4)
                 ! *** TETRAHEDRON CONTRIBUTION TO THE DENSITIES OF STATES
                 Yt(ie) = c(1) + c(2) + c(3) + c(4) + Yt(ie)
                 DO ip = 1 , Npds
                    Yp(ip,ie) = c(1)*Char(ip,nb,i1) + c(2)          &
                         & *Char(ip,nb,i2) + c(3)              &
                         & *Char(ip,nb,i3) + c(4)              &
                         & *Char(ip,nb,i4) + Yp(ip,ie)
                 ENDDO
              ELSE
                 p(1) = facp
                 p(2) = facp
                 p(3) = facp
                 p(4) = facp
              ENDIF
              ! *** TETRAHEDRON CONTRIBUTION TO THE NUMBERS OF STATES
              Zt(ie) = p(1) + p(2) + p(3) + p(4) + Zt(ie)
              DO ip = 1 , Npds
                 Zp(ip,ie) = p(1)*Char(ip,nb,i1) + p(2)             &
                      & *Char(ip,nb,i2) + p(3)*Char(ip,nb,i3)  &
                      & + p(4)*Char(ip,nb,i4) + Zp(ip,ie)
              ENDDO
           ENDIF
20      ENDDO
     ENDDO
  ENDDO
END SUBROUTINE PDSTET
