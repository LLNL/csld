      subroutine gbox(emax, alat, qb, n1, n2, n3)
c============================================================
c     Find the max. dimensions of the reciprocal-space
c     box which contains all plane waves with Ekin < Emax.
c
c  DATE: Wed Oct 16 14:54:34 MDT 1996
c============================================================
      implicit none

      real*8, intent(IN)                       :: emax, alat
      real*8, dimension(3,3), intent(IN)       :: qb
      integer, intent(OUT)                     :: n1,n2,n3

      real*8 Gmax, twopi, tpiba
      real*8 q11,q12,q13,q22,q23,q33,rvol,range1,range2,range3
      real*8, external                         :: tripl

      twopi = 8d0*datan(1d0)
      tpiba = twopi / alat
      Gmax = sqrt( 2d0 * emax ) / tpiba

      q11 = DOT_PRODUCT( qb(:,1), qb(:,1) )
      q12 = DOT_PRODUCT( qb(:,1), qb(:,2) )
      q13 = DOT_PRODUCT( qb(:,1), qb(:,3) )
      q22 = DOT_PRODUCT( qb(:,2), qb(:,2) )
      q23 = DOT_PRODUCT( qb(:,2), qb(:,3) )
      q33 = DOT_PRODUCT( qb(:,3), qb(:,3) )

      rvol = abs( tripl( qb(:,1), qb(:,2), qb(:,3) ) )

      range1 = rvol / sqrt(q22*q33 - q23**2)
      range2 = rvol / sqrt(q11*q33 - q13**2)
      range3 = rvol / sqrt(q11*q22 - q12**2)

      n1 = int ( Gmax / range1 + 1d0 )
      n2 = int ( Gmax / range2 + 1d0 )
      n3 = int ( Gmax / range3 + 1d0 )

      return
      end
