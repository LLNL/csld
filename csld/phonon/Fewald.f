#if defined(CRAY)
#define ERFC derfc
#else
#define ERFC derfc
#endif
      subroutine Fewald(rb, eta, errlim, atp, nat, Zh, frc)
c=================================================================
c     Perform the Ewald summations to obtain the electrostatic
c     lattice contribution to the forces.
c
c   INPUT:
c     rb - the real-space basis;
c     eta - parameter for the Ewald summations;
c     errlim - the requested accuracy of summations;
c     atp - the position vectors for the atoms in the basis;
c     nat - number of atoms in the unit cell
c     Zh - charges of atoms;
c
c   OUTPUT:
c     frc - ion-ion forces;
c
c   DATE: Mon Jul 24 12:54:01 METDST 1995
c=================================================================
      implicit none
      integer nat
      real*8 rb(3,3), eta, atp(3,nat), Zh(nat), errlim, frc(3,nat)
      
      integer n1,n2,n3, i1,i2,i3, ifrc, iat, m, n
      real*8 dtau(3), vg(3), vr(3), ekmax, pi, pisqt, gsq, x, omega,
     $     taumax, rsq, rmax, etasq, qb(3,3), tpi, tpi2, y
      real*8 t0, t1

      real*8, external                               :: tripl
      real*8, external                               :: ERFC
      real*8, external                               :: second

#ifdef TIMINGS
      t0 = second()
#endif

      if (errlim > 1d0) stop 'errlim in Fewald() is too large.'
      etasq = eta*eta
      pi = 4d0*atan(1d0)
      pisqt = sqrt(pi)
      tpi = 2d0*pi
      tpi2 = tpi*tpi
      omega = abs(tripl(rb(1,1), rb(1,2), rb(1,3)))
      call lattc(rb,qb)
      frc = 0d0

      x = 2d0 * eta * sqrt(-log(errlim))
      ekmax = x**2 / 2d0

      call gbox(ekmax, 1d0, qb, n1, n2, n3)
#ifdef DEBUG
      write(6,*) 'Fewald: recip-space sums are over ', n1, n2, n3
#endif

      do i1=-n1,n1
         do i2=-n2,n2
            do i3=-n3,n3

               if(i1==0 .and. i2==0 .and. i3==0) cycle

               vg = tpi*(i1*qb(:,1) + i2*qb(:,2) + i3*qb(:,3))
               gsq = vg(1)*vg(1) + vg(2)*vg(2) + vg(3)*vg(3)

               do ifrc=1,nat
                  do iat=1,nat
c$$$                        dtau = atp(:,ifrc) - atp(:,iat)
c$$$                        x = - 2d0 * exp(-gsq/(4d0*etasq))
c$$$     $                       * sin(vg(1)*dtau(1)+vg(2)*dtau(2)
c$$$     $                       + vg(3)*dtau(3)) / gsq
c$$$                        x = x * Zh(ifrc) * Zh(iat)
c$$$                        frc(:,ifrc) = frc(:,ifrc) +  x * vg(:)

                     frc(:,ifrc) = frc(:,ifrc) - 2*Zh(iat)*vg(:)
     $                    * exp(-gsq/(4d0*etasq))
     $                    * sin( vg(1)*(atp(1,ifrc)-atp(1,iat))
     $                          +vg(2)*(atp(2,ifrc)-atp(2,iat))
     $                          +vg(3)*(atp(3,ifrc)-atp(3,iat)))
     $                    / gsq
                  enddo
               enddo

            enddo
         enddo
      enddo

      do ifrc = 1, nat
         frc(:,ifrc) = 4d0 * pi * Zh(ifrc) * frc(:,ifrc) / omega
      end do

      taumax = 0d0
      do m=1,nat
         do n=1,m
            dtau = atp(:,m) - atp(:,n)
            rsq = dot_product(dtau, dtau)
            taumax = max(taumax, rsq)
         enddo
      enddo

      x = sqrt(-log(errlim)) / eta
      y = sqrt(taumax)
      rmax = (x+y)**2
      ekmax = tpi2 * rmax / 2d0
      call gbox(ekmax, 1d0, rb, n1, n2, n3)
#ifdef DEBUG
      write(6,*) 'Fewald: real space sums over ', n1, n2, n3
#endif
      do i1=-n1,n1
         do i2=-n2,n2
            do i3=-n3,n3

               if(i1==0 .and. i2==0 .and. i3==0) cycle
               vr = i1*rb(:,1) + i2*rb(:,2) + i3*rb(:,3)

               do ifrc=1,nat
                  do iat=1,nat
                     dtau = atp(:,ifrc) - atp(:,iat) + vr
                     rsq = sqrt( dtau(1)*dtau(1)
     $                    + dtau(2)*dtau(2) + dtau(3)*dtau(3) )
#ifdef DEBUG
                     if(rsq < 1d-10) then
                        write(6,200) ifrc, iat
                        write(6,210) atp(:,ifrc)
                        write(6,210) atp(:,iat)
                        stop 'Atomic positions coincide in Fewald.'
                     endif
#endif
                     x = -2d0/pisqt * eta * exp(-(rsq*eta)**2)/rsq - 
     $                    ERFC(eta*rsq)/rsq**2
                     x = x * 2d0*Zh(ifrc)*Zh(iat) / rsq
                     frc(:,ifrc) = frc(:,ifrc) + x * dtau(:)

                  enddo
               enddo

            enddo
         enddo
      enddo

c     Term with l=0 in the real space summation.
      do ifrc=1,nat
         do iat=1,nat
            if(ifrc /= iat) then
               dtau = atp(:,ifrc) - atp(:,iat)
               rsq = sqrt( dot_product(dtau,dtau) )
               x = -2d0/pisqt * eta *exp(-(rsq*eta)**2)/rsq - 
     $              ERFC(eta*rsq)/rsq**2
               x = x * 2d0 * Zh(ifrc) * Zh(iat) / rsq
               frc(:,ifrc) = frc(:,ifrc) + x * dtau(:)
            endif
         enddo
      enddo

      frc = - frc * 5d-1

#ifdef TIMINGS
      t1 = second()
      print *, 'CPU time in Fewald:', t1-t0
#endif
      return

 200  format(' Positions for atoms',2i4,' coincide.')
 210  format(3(2x,f16.8))
      end
