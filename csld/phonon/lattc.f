      subroutine lattc(rb,qb)
c  sets up the real and reciprocal space lattice vectors
      implicit none
      real*8 rb(3,3), qb(3,3)
      real*8 vol
      integer m,k
      real*8 tripl
      external tripl

      call cross(rb(1,2), rb(1,3), qb(1,1))
      call cross(rb(1,3), rb(1,1), qb(1,2))
      call cross(rb(1,1), rb(1,2), qb(1,3))
      vol = tripl(rb(1,1), rb(1,2), rb(1,3))
      qb = qb * (1.d0/vol)

      return
      end
