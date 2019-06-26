!original version
! http://www.unamur.be/sciences/physique/administration/tetpack/dostet.f
!*==DOSTET.spg  processed by SPAG 6.72Dc at 03:33 on 26 Sep 2016
SUBROUTINE DOSTET(Ener,Idime,Nband,Idef,Ntet,Xe,Ne,Y,Z)
  !     COMPUTE A DENSITY OF STATES USING THE TETRAHEDRONS METHOD.
  !     XE CONTAINS THE ENERGIES, Y AND Z RETURN THE RELATED DENSITY OF
  !     STATES AND THE INTEGRATED DENSITY OF STATES, RESPECTIVELY.
  IMPLICIT NONE
  !*--DOSTET7
  !*** Start of declarations inserted by SPAG
  integer, intent(in)   :: Idime, Nband, Ntet, Ne
  real(8), intent(in)    :: Ener(Idime,*), Xe(Ne)
  integer, intent(in)   :: Idef(5, Ntet)
  real(8), intent(out)   :: Y(Ne) , Z(Ne)

! local variables
  REAL(8) c(4), cc , e , e1 , e12 , e13 , e14 , e2 , e23 , e24 , e3 ,  &
       & e34 , e4 , ee1 , ee2 , ee3 , ee4 ,  eps , facy
  INTEGER i , ia , ib , ic , id , itet , ivol , ix , j , nb 
  REAL(8) surfac , unite , volume
  !*** End of declarations inserted by SPAG
  REAL(4) avol
  !DIMENSION Ener(Idime,1) , Xe(1) , Y(1) , Z(1) , Idef(5,1) , c(4)
  EQUIVALENCE (ivol,avol)
  EQUIVALENCE (c(1),e1)
  EQUIVALENCE (c(2),e2)
  EQUIVALENCE (c(3),e3)
  EQUIVALENCE (c(4),e4)
  DATA eps/1.0D-05/
  Y = 0.D0
  Z = 0.D0
  !
  !     LOOP OVER THE TETRAHEDRONS
  DO itet = 1 , Ntet
     ia = Idef(1,itet)
     ib = Idef(2,itet)
     ic = Idef(3,itet)
     id = Idef(4,itet)
     ivol = Idef(5,itet)
     !
     !     LOOP OVER THE BANDS
     DO nb = 1 , Nband
        !
        ! *** DEFINE E1, E2, E3, E4, AS THE CORNER ENERGIES ORDERED BY
        ! *** DECREASING SIZE
        c(1) = Ener(nb,ia)
        c(2) = Ener(nb,ib)
        c(3) = Ener(nb,ic)
        c(4) = Ener(nb,id)
        DO i = 1 , 4
           cc = c(i)
           do j = i + 1, 4
              IF ( cc.LT.c(j) ) THEN
                 c(i) = c(j)
                 c(j) = cc
                 cc = c(i)
              ENDIF
           end do
        ENDDO
        unite = 1.0D0
        IF ( e1.GT.e4 ) unite = e1 - e4
        e12 = (e1-e2)/unite
        e13 = (e1-e3)/unite
        e14 = (e1-e4)/unite
        e23 = (e2-e3)/unite
        e24 = (e2-e4)/unite
        e34 = (e3-e4)/unite
        facy = 3.D0*DBLE(avol)/unite
        DO ix = 1 , Ne
           e = Xe(ix)
           surfac = 0.D0
           volume = 1.D0
           IF ( e.LE.e1 ) THEN
              volume = 0.D0
              IF ( e.GE.e4 ) THEN
                 ee1 = (e-e1)/unite
                 IF ( DABS(ee1).LT.eps ) ee1 = 0.D0
                 ee2 = (e-e2)/unite
                 IF ( DABS(ee2).LT.eps ) ee2 = 0.D0
                 ee3 = (e-e3)/unite
                 IF ( DABS(ee3).LT.eps ) ee3 = 0.D0
                 ee4 = (e-e4)/unite
                 IF ( DABS(ee4).LT.eps ) ee4 = 0.D0
                 IF ( e.GT.e3 ) THEN
                    IF ( e.GT.e2 ) THEN
                       ! *** E2.LT.E.AND.E.LE.E1
                       surfac = (ee1/e12)*(ee1/e13)
                       volume = 1.D0 + surfac*ee1
                    ELSE
                       ! *** E3.LT.E.AND.E.LE.E2
                       surfac = -(ee3*ee2/e23+ee4*ee1)/e13/e24
                       volume = (0.5D0*ee3*(2.D0*e13*e34+e13*ee4-e34&
                            & *ee1-2.D0*ee1*ee4+ee3*(ee3-3.D0*ee2)&
                            & /e23)/e13+e34*e34)/e24
                    ENDIF
                    ! *** E4.LE.E.AND.E.LE.E3
                 ELSEIF ( e4.EQ.e3 ) THEN
                    IF ( e3.GE.e2 ) THEN
                       IF ( e2.EQ.e1 ) THEN
                          surfac = 1.0D+15
                          volume = 0.5D0
                       ELSE
                          surfac = 1.D0/e12
                       ENDIF
                    ENDIF
                 ELSE
                    surfac = (ee4/e34)*(ee4/e24)
                    volume = surfac*ee4
                 ENDIF
              ENDIF
           ENDIF
           Y(ix) = Y(ix) + facy*surfac
           Z(ix) = Z(ix) + DBLE(avol)*volume
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE DOSTET

