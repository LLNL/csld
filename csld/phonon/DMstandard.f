      subroutine DMstandard(q, DM, natoms, apos)
      implicit none
      
      integer natoms
      real*8 q(3), apos(3,natoms)
      complex*16 DM(3,natoms,3,natoms)

      integer iat1, iat2, j1, j2
      real*8  AposDiff(3), QxR, x(3), tpi
      complex*16 expQxR

      tpi = 8d0*atan(1d0)

c     Conjugate the dynamical matrix, since we read row into column in main.f !!!
      DM = conjg( DM )
      do iat1=1,Natoms
         do iat2=1,Natoms
            AposDiff = apos(:,iat2) - apos(:,iat1)
            QxR = tpi * dot_product(q, AposDiff)
            expQxR = dcmplx( cos(QxR), -sin(QxR) )

            DM(:,iat1,:,iat2) = DM(:,iat1,:,iat2)*expQxR

         end do
      end do

c$$$      print 110, q
c$$$      print 100, DM
c$$$      print *
c$$$ 100  format(6('(',f6.2,',',f6.2,')',2x))
c$$$ 110  format('   q=   ', 3f12.6)

      return
      end
