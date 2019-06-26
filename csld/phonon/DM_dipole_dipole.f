! #define DEBUG      
      subroutine DM_dipole_dipole(rb, qb, eta_in, errlim, apos, nat,
     $     Zeff, eps, qph, dmat)
c*****************************************************************
c     Find the lattice contribution to the dynamical matrix.
c     Uses Ewald transformation to split the lattice sums into
c     rapidly convergent real and reciprocal space summations.
c     Works with effective charges given as tensors, includes
c     dielectric screening due to electronic polarization.
c
c   INPUT:
c     rb - real-space basis;
c     qb - reciprocal space basis;
c     eta - parameter for the Ewald summations;
c              use find_Ewald_eta_screened to determine optimal value
c     errlim - requested accuracy;
c     nat - number of atoms in the unit cell
c     apos - atomic position in the unit cell;
c     Zeff - Born effective charge tensorse;
c     eps - dielectric matrix;
c     qph - phonon wave vector;
c
c   OUTPUT:
c     dmat - electrostatic contribution to the dynamical matrix,
c            it should be divided by the lattice constant to get
c            the same normalization as the DMs calculated by lrx
c
c     This implementation closely follows Gonze, X., & Lee, C. (1997).
c     "Dynamical matrices, born effective charges, dielectric permittivity
c     tensors, and interatomic force constants from density-functional
c     perturbation theory." Physical Review B, 55(16), 10355--10368.
c
c  DATE: Sun Jan  1 11:45:23 PST 2012
c*****************************************************************
      implicit none

      integer, intent(IN)                    :: nat
      real*8, dimension(3,3), intent(IN)     :: rb, qb, eps
      real*8, intent(IN)                     :: eta_in, errlim
      real*8, dimension(3,nat), intent(IN)   :: apos
      real*8, dimension(3,3,nat), intent(IN) :: Zeff
      real*8, dimension(3), intent(IN)       :: qph
      complex*16, dimension(3,nat,3,nat), intent(OUT) :: dmat

      integer i1, i2, i3, n1, n2, n3, j1, j2, m, n, j
      real*8 pi, tpi, piroot, omega, x, x1, x2, x3, y, etasq, ekmax,
     $     G(3), dtau(3), Gpq(3), g2, Gpq2, rsq, rmod, dr(3), vr(3),
     $     Delta(3), tpi2, taumax, eps_inv(3,3), deteps
      complex*16 cz, cfac, tdm(3,3)
      complex*16, dimension(:,:,:,:), allocatable    :: dmat0
      complex*16, dimension(:), allocatable        :: czz
      real*8 eta

c     LAPACK related variables
      integer lwork
      parameter (lwork=30)
      integer info
      real*8 teps(3,3), evr(3), evi(3), vleft(3,3), vright(3,3),
     $     work(lwork)

      real*8 tripl, derfc
      external tripl, derfc

      eta=eta_in

      deteps = (eps(1,1)*(eps(2,2)*eps(3,3)-eps(2,3)*eps(3,2))
     $     - eps(2,1)*(eps(1,2)*eps(3,3)-eps(1,3)*eps(3,2))
     $     + eps(3,1)*(eps(1,2)*eps(2,3)-eps(2,2)*eps(1,3)) )

      eps_inv(1,1) = eps(2,2)*eps(3,3)-eps(2,3)*eps(3,2)
      eps_inv(2,1) = eps(2,3)*eps(3,1)-eps(2,1)*eps(3,3)
      eps_inv(3,1) = eps(2,1)*eps(3,2)-eps(2,2)*eps(3,1)
      eps_inv(1,2) = eps(1,3)*eps(3,2)-eps(1,2)*eps(3,3)
      eps_inv(2,2) = eps(1,1)*eps(3,3)-eps(1,3)*eps(3,1)
      eps_inv(3,2) = eps(1,2)*eps(3,1)-eps(1,1)*eps(3,2)
      eps_inv(1,3) = eps(1,2)*eps(2,3)-eps(1,3)*eps(2,2)
      eps_inv(2,3) = eps(1,3)*eps(2,1)-eps(1,1)*eps(2,3)
      eps_inv(3,3) = eps(1,1)*eps(2,2)-eps(1,2)*eps(2,1)

      eps_inv = eps_inv / deteps
      
#ifdef DEBUG
      print *, 'input eta=', eta
#ifdef READPARAMETER
      open (unit=99, file='eta.txt', status='old', action='read')
      read(99, *), eta
      close(99)
      print *, 'read eta=', eta
#endif
      print '(a/3(3f12.4/))', ' Dielectric matrix:', eps
      print '(a/3(3f12.4/))', ' Inverse dielectric matrix:', eps_inv
      print '(a,f12.4)', ' Determinant of dielectric matrix =', deteps
      print '(3f12.6)', MATMUL( eps, eps_inv )
#endif

      etasq = eta*eta
      pi = 4d0*atan(1d0)
      piroot = sqrt(pi)
      tpi = 2d0*pi
      tpi2 = tpi*tpi
      omega=abs(tripl(rb(1,1), rb(1,2), rb(1,3)))

      allocate( dmat0(3,nat,3,nat) )

      dmat = (0d0, 0d0)
      dmat0 = (0d0, 0d0)

      teps = eps
!      call LA_DGEEV('N','V',3,teps,3,evr,evi,vleft,3,vright,3,
      call DGEEV('N','V',3,teps,3,evr,evi,vleft,3,vright,3,
     $     work,lwork,info)
      if( info /= 0 ) print *, 'DGEEV error: info=', info
#ifdef DEBUG
      print '(2x,a30)', 'Eigenvalues and eigenvectors:'
      do j = 1, 3
         print '(2f10.4,5x,3f10.4)', evr(j), evi(j), vright(:,j)
      end do
#endif

      x = 2d0 * eta * sqrt(-log(errlim))
      y = sqrt(dot_product(qph,qph))
      ekmax = MINVAL(1d0/evr)*(x+y)**2 / 2d0
      call gbox(ekmax, 1d0, qb, n1, n2, n3)
      
      ALLOCATE( czz(nat) )

#ifdef DEBUG
      print '(a,3i4)', 'DM_dipole_dipole: recip-space sums are over ',
     $     n1, n2, n3
#ifdef READPARAMETER
      open (unit=99, file='nq.txt', status='old', action='read')
      read(99, *), n1, n2, n3
      close(99)
#endif
#endif

      do i1=-n1,n1
         do i2=-n2,n2
            do i3=-n3,n3

               G = i1*qb(:,1) + i2*qb(:,2) + i3*qb(:,3)
               Gpq = G + qph
               G2   = dot_product(G, MATMUL(eps,G) )
               Gpq2 = dot_product(Gpq, MATMUL(eps, Gpq))
               
               do m=1,nat
                  x = tpi*(G(1)*apos(1,m) + G(2)*apos(2,m)
     $                 + G(3)*apos(3,m))
                  czz(m) = dcmplx( cos(x), sin(x) )
               end do

c     Q /= 0 contribution.
               if (Gpq2 > 1d-25) then
                  do m=1,nat
                     do n=1,nat
                        cfac = exp(-tpi2*Gpq2/(4d0*etasq)) / Gpq2
                        do j1=1,3
                           do j2=1,3
                              tdm(j1,j2) = cfac * czz(m) * CONJG(czz(n))
     $                             * Gpq(j1) * Gpq(j2)
                           end do
                        end do
                        dmat(:,m,:,n) = dmat(:,m,:,n) + tdm
                     end do
                  end do
#ifdef DEBUG
               if ((i1==n1).and.(i2==n2).and.(i3==n3)) then
                  print *, 'recip fac=', abs(cfac)
               end if
#endif
               end if

c     Q=0 contribution
               if (G2 > 1d-25) then
                  do m=1,nat
                     do n=1,nat
                        cfac = czz(m) * CONJG(czz(n)) *
     $                       exp(-tpi2*G2/(4d0*etasq)) / G2
                        do j1=1,3
                           do j2=1,3
                              tdm(j1,j2) = cfac * G(j1)*G(j2)
                           enddo
                        enddo
                        dmat0(:,m,:,n) = dmat0(:,m,:,n) + tdm
                     end do
                  end do
               end if

            enddo
         enddo
      enddo
      dmat = (4d0*pi/omega) * dmat
      dmat0 = (4d0*pi/omega) * dmat0
      DEALLOCATE( czz )

!#ifdef DEBUG
!C test hermitian
!      print *, 'after recip sum testing dmat'
!      call diff_hermitian(nat*3, dmat)
!      print *, 'testing dmat0'
!      call diff_hermitian(nat*3, dmat0)
!#endif

c   Debugging section for testing reciprocal space loop by comparing
c   with the recip space loop in DMewald routine from the lrx code
c$$$      do m = 1, nat
c$$$         do n = 1, nat
c$$$            dmat(:,m,:,n) = MATMUL( Zeff(:,:,m),
c$$$     $                              MATMUL( dmat(:,m,:,n),
c$$$     $                                      TRANSPOSE(Zeff(:,:,n))
c$$$     $                              )
c$$$     $                      )
c$$$            dmat0(:,m,:,n) = MATMUL( Zeff(:,:,m),
c$$$     $                              MATMUL( dmat0(:,m,:,n),
c$$$     $                                      TRANSPOSE(Zeff(:,:,n))
c$$$     $                              )
c$$$     $                      )
c$$$         end do
c$$$      end do
c$$$
c$$$      do m = 1, nat
c$$$         do n = 1, nat
c$$$            dmat(:,m,:,m) = dmat(:,m,:,m) - dmat0(:,m,:,n)
c$$$         end do
c$$$      end do
c$$$
c$$$      deallocate( dmat0 )
c$$$      return

      taumax = 0d0
      do m=1,nat
         do n=1,m
            dr = apos(:,m) - apos(:,n)
            rsq = dot_product(dr, MATMUL(eps_inv, dr))
            taumax = max(taumax, rsq)
         enddo
      enddo
      x = sqrt(-log(errlim)) / eta
      y = sqrt(taumax)
      ekmax = MINVAL(evr) * tpi2 * (x+y)**2 / 2d0
      call gbox(ekmax, 1d0, rb, n1, n2, n3)

#ifdef DEBUG
      print '(a,3i4)', 'DM_dipole_dipole: real-space sums over',n1,n2,n3
#ifdef READPARAMETER
      open (unit=99, file='nr.txt', status='old', action='read')
      read(99, *), n1, n2, n3
      close(99)
#endif
#endif


      do i1=-n1,n1
         do i2=-n2,n2
            do i3=-n3,n3

               vr = i1*rb(:,1) + i2*rb(:,2) + i3*rb(:,3)

               do m=1,nat
                  do n=1,nat
                     if( n/=m .or. i1/=0 .or. i2/=0 .or. i3/=0 ) then

                        dr = vr + apos(:,m) - apos(:,n)
                        Delta = MATMUL(eps_inv, dr)
                        rsq = dot_product(Delta, dr)
                        rmod = sqrt(rsq)
                        
c     This line should work for the "standard" definition of dynamical matrix
c     without the additional factor exp[iq(a_i-a_j)].
!                        x = tpi * dot_product(qph, vr)
                        x = tpi * dot_product(qph, dr)
                        cz = dcmplx(cos(x), -sin(x))

                        x1 = derfc(eta*rmod) / (eta*rmod)**3
                        x2 = (2d0/piroot) * exp(-rsq*etasq)
                        x3 = x2 / (eta*rmod)**2

                        do j1=1,3
                           do j2=1,3
                              tdm(j1,j2) = - eps_inv(j1,j2) * (x1+x3)
     $                             + (Delta(j1)*Delta(j2)/rsq)
     $                                 * ( 3d0*(x1+x3) + 2d0*x2)
                           enddo
                        enddo
                        tdm = cz * eta**3 * tdm / sqrt(deteps)

                        dmat(:,m,:,n) = dmat(:,m,:,n) - tdm

                        do j1=1,3
                           do j2=1,3
                              tdm(j1,j2) = - eps_inv(j1,j2) * (x1+x3)
     $                             + (Delta(j1)*Delta(j2)/rsq)
     $                                 * ( 3d0*(x1+x3) + 2d0*x2)
                           enddo
                        enddo
                        tdm = eta**3 * tdm / sqrt(deteps)

                        dmat0(:,m,:,n) = dmat0(:,m,:,n) - tdm

#ifdef DEBUG
               if ((i1==n1).and.(i2==n2).and.(i3==n3) 
     $            .and.(m==1).and.(n==1)) then
             write(*,'(A,4E9.2)') '|fac|=', maxval(abs(tdm)),x1,x2,x3
               end if
               if ((i1==0).and.(i2==0).and.(i3==0) 
     $            .and.(m==1).and.(n==2)) then
       write(*,'(A,5E9.2)') '|fac 1-2|=', maxval(abs(tdm)),x1,x2,x3,rmod
               end if
#endif
                     end if
                  enddo
               enddo

            enddo
         enddo
      enddo

!#ifdef DEBUG
!C test hermitian
!      print *, 'after real space sum testing dmat'
!      call diff_hermitian(nat*3, dmat)
!      print *, 'testing dmat0'
!      call diff_hermitian(nat*3, dmat0)
!#endif

c     Constant term for the limiting contribution. Third term in Eq.(73)
      do m = 1, nat
         dmat(:,m,:,m) = dmat(:,m,:,m) - (4d0/(3d0*piroot)) * eta**3
     $        * eps_inv / sqrt(deteps)
         dmat0(:,m,:,m) = dmat0(:,m,:,m) - (4d0/(3d0*piroot)) * eta**3
     $        * eps_inv / sqrt(deteps)
      end do


! Eq.(72), \hat{C} = linear transform of \bar{C}
! Note the typo in Eq.(72) RHS, should be C_{\alpha' \beta'}
      do m = 1, nat
         do n = 1, nat
            dmat(:,m,:,n) = MATMUL( Zeff(:,:,m),
     $                              MATMUL( dmat(:,m,:,n),
     $                                      TRANSPOSE(Zeff(:,:,n))
     $                              )
     $                      )
            dmat0(:,m,:,n) = MATMUL( Zeff(:,:,m),
     $                              MATMUL( dmat0(:,m,:,n),
     $                                      TRANSPOSE(Zeff(:,:,n))
     $                              )
     $                      )
         end do
      end do

!#ifdef DEBUG
!C test hermitian
!      print *, 'after Z* testing dmat'
!      call diff_hermitian(nat*3, dmat)
!      print *, 'testing dmat0'
!      call diff_hermitian(nat*3, dmat0)
!#endif

c     Enforce acoustic sum rule.
      do m = 1, nat
         do n = 1, nat
            dmat(:,m,:,m) = dmat(:,m,:,m) -  dmat0(:,m,:,n)
!C      dmat(:,m,:,m) = dmat(:,m,:,m)-(dmat0(:,m,:,n)+dmat0(:,n,:,m))/2.0
         end do
      end do

#ifdef DEBUG
C test hermitian
      print *, 'after ASR testing dmat'
      call diff_hermitian(nat*3, dmat)
#endif

#ifdef DEBUG

C test ASR
!      print *, 'DEBUG q=', sqrt(dot_product(qph,qph))
!      if (sqrt(dot_product(qph,qph))<1E-4) then
!         print *, 'testing ASR'
!      do j1=1, 3
!         do j2=1, 3
!            do m=1, nat
!               cfac=(0,0)
!               do n=1, nat
!                  cfac=cfac+ dmat(j1, m, j2, n)
!               end do
!               if (abs(cfac)>1E-7) then
!            write(*,'(A,3I4,2F12.8)'), 'ERROR in ASR m=', j1,j2,m,cfac
!               end if
!            end do
!         end do
!      end do
!      end if
!
!      open (unit = 1, file = "dmat.dat")
!      do n=1,nat
!         do j2=1,3
!            do m=1,nat
!               do j1=1,3
!                  write(1, *) dmat(j1,m,j2,n)
!               end do
!            end do
!         end do
!      end do
!      close(1)

      open(unit = 10, file='dmat.bin',form='unformatted')  ! create a new file, or overwrite an existing on
      write(10) dmat
      close(10) ! close the file
#endif

      deallocate( dmat0 )

      return
      end

#ifdef DEBUG
      subroutine diff_hermitian(size, mat)
      implicit none

      integer, intent(IN)                    :: size
      complex*16, intent(IN)     :: mat(size,size)

      integer i, j, count, maxcount
      real*8 diff

      maxcount=10
      count=0
      do i=1,size
         do j=i+1, size
             diff=abs(mat(i,j)-conjg(mat(j,i)))
             if (diff>1E-14) then
                 write(*,'(A,2I6,3E12.3)') 'Err hermit ', i,j,diff,
     $ real(mat(i,j)-conjg(mat(j,i))),imag(mat(i,j)-conjg(mat(j,i)))
                 count= count+1
                 if (count>maxcount) then
                    write(*,*) 'max non-hermitian elements encountered'
                    return
                 end if 
             end if
         end do
      end do
      return
      end subroutine
#endif
