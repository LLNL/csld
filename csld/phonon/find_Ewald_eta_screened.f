#if defined(CRAY)
#define ERFC derfc
#else
#define ERFC derfc
#endif
      subroutine find_Ewald_eta_screened(rb, errlim, nat,
     $     atp, eps, eta)
c=================================================================
c     Calculates optimal "eta" for Ewald summations.
c
c   INPUT:
c     rb - real-space basis
c     errlim - requested accuracy of summations
c     nat - number of atoms
c     atp - atomic positions
c     eps - high-frequency dielectric tensor
c     
c   OUTPUT:
c     eta - parameter for the Ewald summations
c
c   DATE: Sun Jan  1 12:08:58 PST 2012
c=================================================================
      implicit none

      integer MAXIT
      parameter( MAXIT = 100 )

      integer nat
      real*8 rb(3,3), eta, errlim, atp(3,nat), eps(3,3)

      integer nq1, nq2, nq3, nr1, nr2, nr3, m, n, iter, j
      real*8 dtau(3), rsq, taumax, ekmax, pi, tpi, tpi2,
     $     qb(3,3), x, y, FACT, etaprev

c     LAPACK related variables
      integer lwork
      parameter (lwork=30)
      integer info
      real*8 teps(3,3), evr(3), evi(3), vleft(3,3), vright(3,3),
     $     work(lwork)

c     Complementary error function.
      real*8, external                          :: ERFC

!      print *
      if(errlim > 1d0) then
         print *, 'errlim in find_Ewald_eta is too large = ', errlim
         stop
      end if

      pi = 4d0*atan(1d0)
      tpi = 2d0*pi
      tpi2 = tpi**2

      teps = eps
      call DGEEV('N','V',3,teps,3,evr,evi,vleft,3,vright,3,
     $     work,lwork,info)
      if( info /= 0 ) print *, 'DGEEV error: info=', info

#ifdef DEBUG
      print '(2x,a30)', 'Eigenvalues and eigenvectors of ',
     $     'dielectric matrix:'
      do j = 1, 3
         print '(2f10.4,5x,3f10.4)', evr(j), evi(j), vright(:,j)
      end do
#endif

      call lattc(rb,qb)
      taumax = 0d0
      do m=1,nat
         do n=1,m
            dtau = atp(:,m) - atp(:,n)
            rsq = dot_product(dtau, dtau)
            taumax = max(taumax, rsq)
         end do
      end do

      if( abs(eta) < 1d-5 ) eta = 1d0
      FACT = 2d0
      etaprev = eta

      do iter = 1, MAXIT
         x = 2d0 * eta * sqrt(-log(errlim))
         ekmax = MINVAL(1d0/evr) * x**2 / 2d0
         call gbox(ekmax, 1d0, qb, nq1, nq2, nq3)

         x = sqrt(-log(errlim)) / eta
         y = sqrt(taumax)
         ekmax = MINVAL(evr) * tpi2 * (x+y)**2 / 2d0
         call gbox(ekmax, 1d0, rb, nr1, nr2, nr3)

#ifdef DEBUG
         print *, 'Eta = ', eta
         print *, 'Previous Eta = ', etaprev
         print *, 'Recip-space sums are over ', nq1, nq2, nq3
         print *, 'Real space sums over ', nr1, nr2, nr3
#endif

         if( nq1*nq2*nq3 > 2*nr1*nr2*nr3 ) then
            if( iter == MAXIT ) then
               print *, 'Warning: possibly suboptimal value of',
     $              ' eta in Ewald sums!'
            else
               if( eta <= etaprev ) then
                  x = eta / FACT
               else
                  x = eta - (eta-etaprev) / 3d0
               end if
               etaprev = eta
               eta = x
            end if
         else if ( 2*nq1*nq2*nq3 < nr1*nr2*nr3 ) then
            if( iter == MAXIT ) then
               print *, 'Warning: possibly suboptimal value of',
     $              ' eta in Ewald sums!'
            else
               if( eta >= etaprev ) then
                  x = eta * FACT
               else
                  x = eta + (etaprev-eta) / 3d0
               end if
               etaprev = eta
               eta = x
            end if
         else
#ifdef DEBUG
            print '(1x,a,f12.4)', 'Found an optimal value of eta=', eta
#endif
            exit
         end if
      end do

#ifdef DEBUG
      print '(a,3i3)', ' Real space sums over  ', nr1, nr2, nr3
      print '(a,3i3)', ' Recip-space sums over ', nq1, nq2, nq3
#endif

      return
      end
