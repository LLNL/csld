!     original version
!     http://www.unamur.be/sciences/physique/administration/tetpack/setk06.f
!     *==SETK06.spg  processed by SPAG 6.72Dc at 07:20 on 26 Sep 2016
SUBROUTINE SETK06(N1,N2,N3,A1,A2,A3,Ptk,Nptk,Idef,Ntet,Nkmx,Ntmx)
  !     SET THE K-POINTS IN ONE HALF OF THE RECIPROCAL CELL FOR A
  !     TRICLINIC LATTICE WITH DIRECT LATTICE VECTORS A1, A2, A3
  !     SYMMETRY IS C1I
  IMPLICIT NONE
  integer, intent(in) :: N1, N2, N3, Nkmx, Ntmx
  integer, intent(out) :: Nptk, Ntet, Idef(5,Ntmx)
  REAL(8), intent(in) :: A1(3) , A2(3) , A3(3)
  REAL(8), intent(out) :: Ptk(4,Nkmx)

  !local var
  REAL(8)  b(3,3) , c , pi , w , wk
  INTEGER i , i1 , i2 , i2max , i3 , i3max , icode ,  ind1 ,  &
       & ind2 , ind3 , ind4 , ind5 , ind6 , ind7 , ind8 , it ,     &
       & ivol , j , k
  INTEGER k1 , k2 , k3 , kk(3) , n ,  n12 ,  n3p1   ,&
       & nn(3) , nn1 , nn2 , nn3
  REAL(4) avol
  !     DIMENSION Ptk(4,Nkmx) , Idef(5,Ntmx) , A1(3) , A2(3) , A3(3) ,    &
  !     & b(3,3) , nn(3) , kk(3)
  EQUIVALENCE (ivol,avol)
  pi = 3.141592653589793238D0

  IF ( N1<=0 .OR. N2<=0 .OR. N3<=0 ) stop ' *** <SETK06> N1, N2 OR N3 IS NOT POSITIVE ***'
  !     *** GENERATE THE RECIPROCAL BASIS VECTORS
  b(1,1) = A2(2)*A3(3) - A2(3)*A3(2)
  b(1,2) = A2(3)*A3(1) - A2(1)*A3(3)
  b(1,3) = A2(1)*A3(2) - A2(2)*A3(1)
  c = A1(1)*b(1,1) + A1(2)*b(1,2) + A1(3)*b(1,3)
  IF ( abs(c)< 1.0E-8 ) stop ' *** <SETK06> VOLUME OF THE PRIMITIVE CELL IS ZERO ***'
  c = pi/c
  b(1,1) = b(1,1)*c/N1
  b(1,2) = b(1,2)*c/N1
  b(1,3) = b(1,3)*c/N1
  b(2,1) = (A3(2)*A1(3)-A3(3)*A1(2))*c/N2
  b(2,2) = (A3(3)*A1(1)-A3(1)*A1(3))*c/N2
  b(2,3) = (A3(1)*A1(2)-A3(2)*A1(1))*c/N2
  b(3,1) = (A1(2)*A2(3)-A1(3)*A2(2))*c/N3
  b(3,2) = (A1(3)*A2(1)-A1(1)*A2(3))*c/N3
  b(3,3) = (A1(1)*A2(2)-A1(2)*A2(1))*c/N3
  !     *** SORT N1,N2,N3 BY DECREASING ORDER
  nn(1) = N1
  nn(2) = N2
  nn(3) = N3
  kk(1) = 1
  kk(2) = 2
  kk(3) = 3
  DO i = 1 , 2
     DO j = i + 1 , 3
        IF ( nn(j)>nn(i) ) THEN
           n = nn(i)
           nn(i) = nn(j)
           nn(j) = n
           k = kk(i)
           kk(i) = kk(j)
           kk(j) = k
        ENDIF
     ENDDO
  ENDDO
  Nptk = 4*N1*N2*N3 + 2*nn(2)*nn(3) + nn(3) + 1
  IF ( Nptk>Nkmx ) STOP '*** <SETK06> NPTK EXCEEDS NKMAX ***'
  Ntet = 24*N1*N2*N3
  IF ( Ntet>Ntmx ) STOP '*** <SETK06> NTET EXCEEDS NTMAX ***'
  !     *** SET THE K-POINTS
  k1 = kk(1)
  k2 = kk(2)
  k3 = kk(3)
  nn1 = nn(1)
  nn2 = nn(2)
  nn3 = nn(3)
!  WRITE (6,99003) Nptk , Ntet ,                               &
!       & nn1*DSQRT(b(k1,1)**2+b(k1,2)**2+b(k1,3)**2) &
!       & , k1 ,                                      &
!       & 2*nn2*DSQRT(b(k2,1)**2+b(k2,2)**2+b(k2,3)   &
!       & **2) , k2 ,                                 &
!       & 2*nn3*DSQRT(b(k3,1)**2+b(k3,2)**2+b(k3,3)   &
!       & **2) , k3
!99003 FORMAT (' SAMPLING ONE HALF THE RECIPROCAL UNIT CELL'/1X,I5,&
!       &' K-POINTS',I7,' TETRAHEDRA'/                       &
!       &' THE EDGES OF THE SAMPLED PARALLELEPIPED HAVE LENGTHS EQUAL TO'&
!       & /3(1X,D10.4,' ALONG B',I1,3X))
  w = 1.0D0/(4*N1*N2*N3)
  Nptk = 0
  icode = 0
  i2max = nn2
  i3max = nn3
  i1 = -nn1
  i2 = -nn2
  i3 = -nn3
20 icode = icode + 1
  IF ( i3==i3max ) THEN
     IF ( i3/=0 ) THEN
        Idef(5,icode) = Nptk - 2*nn3 + 1
        i2 = i2 + 1
        IF ( i2<i2max ) THEN
           i3 = -nn3
           GOTO 20
        ELSE
           IF ( i2==0 ) THEN
              i3max = 0
              i3 = -nn3
              GOTO 20
           ENDIF
           n12 = 4*nn2*nn3*(i1+nn1)
           DO i3 = 1 , 2*nn3
              icode = icode + 1
              n12 = n12 + 1
              Idef(5,icode) = n12
           ENDDO
           icode = icode + 1
           Idef(5,icode) = n12 - 2*nn3 + 1
           i1 = i1 + 1
           IF ( i1<0 ) THEN
              i2 = -nn2
              i3 = -nn3
              GOTO 20
           ELSEIF ( i1==0 ) THEN
              i2max = 0
              i2 = -nn2
              i3 = -nn3
              GOTO 20
           ENDIF
        ENDIF
     ENDIF
     Nptk = Nptk + 1
     Ptk(1,Nptk) = 0.0D0
     Ptk(2,Nptk) = 0.0D0
     Ptk(3,Nptk) = 0.0D0
     Ptk(4,Nptk) = w/2.0D0
     Idef(5,icode) = Nptk
     n12 = Nptk
     DO i3 = 1 , nn3
        icode = icode + 1
        n12 = n12 - 1
        Idef(5,icode) = n12
     ENDDO
     DO i2 = 1 , nn2
        icode = icode + 1
        Idef(5,icode) = n12 - 2*nn3
        DO i3 = 1 , 2*nn3
           icode = icode + 1
           n12 = n12 - 1
           Idef(5,icode) = n12
        ENDDO
     ENDDO
     !     *** DEFINE THE TETRAHEDRA
     Ntet = 0
     n3p1 = 2*nn3 + 1
     n12 = n3p1*(2*nn2+1)
     DO i1 = 1 , nn1
        DO i2 = 1 , 2*nn2
           ind7 = (i1-1)*n12 + n3p1*(i2-1)
           DO i3 = 1 , 2*nn3
              ind7 = ind7 + 1
              ind6 = ind7 + n12
              ind2 = ind6 + n3p1
              ind1 = ind2 + 1
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind6)
              Idef(3,Ntet) = Idef(5,ind2)
              Idef(4,Ntet) = Idef(5,ind1)
              ind8 = ind7 + 1
              ind5 = ind6 + 1
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind6)
              Idef(3,Ntet) = Idef(5,ind5)
              Idef(4,Ntet) = Idef(5,ind1)
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind8)
              Idef(3,Ntet) = Idef(5,ind5)
              Idef(4,Ntet) = Idef(5,ind1)
              ind3 = ind7 + n3p1
              ind4 = ind3 + 1
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind8)
              Idef(3,Ntet) = Idef(5,ind4)
              Idef(4,Ntet) = Idef(5,ind1)
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind3)
              Idef(3,Ntet) = Idef(5,ind4)
              Idef(4,Ntet) = Idef(5,ind1)
              Ntet = Ntet + 1
              Idef(1,Ntet) = Idef(5,ind7)
              Idef(2,Ntet) = Idef(5,ind3)
              Idef(3,Ntet) = Idef(5,ind2)
              Idef(4,Ntet) = Idef(5,ind1)
           ENDDO
        ENDDO
     ENDDO
     avol = 1.D0/dble(Ntet)
     DO it = 1 , Ntet
        Idef(5,it) = ivol
     ENDDO
     RETURN
  ELSE
     wk = w
     IF ( i1==-nn1 ) wk = w/2.0D0
     IF ( i1==0 ) THEN
        IF ( i2==-nn2 ) wk = w/2.0D0
        IF ( i2==0 .AND. i3==-nn3 ) wk = w/2.0D0
     ENDIF
     Nptk = Nptk + 1
     Ptk(1,Nptk) = i1*b(k1,1) + i2*b(k2,1) + i3*b(k3,1)
     Ptk(2,Nptk) = i1*b(k1,2) + i2*b(k2,2) + i3*b(k3,2)
     Ptk(3,Nptk) = i1*b(k1,3) + i2*b(k2,3) + i3*b(k3,3)
     Ptk(4,Nptk) = wk
     Idef(5,icode) = Nptk
     i3 = i3 + 1
     GOTO 20
  ENDIF


END SUBROUTINE SETK06

