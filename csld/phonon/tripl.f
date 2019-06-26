      real*8 function tripl(a,b,c)
      implicit none
      real*8 a(3),b(3),c(3)
      tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
     .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
      return
      end
