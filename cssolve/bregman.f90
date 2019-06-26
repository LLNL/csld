MODULE bregman

implicit none 

!
! Public subroutines for compressive sensing using Bregman iterations.
! For large problems, use SplitBregman.
!
! REFERENCES:
!
! The Bregman algorithm is described in W. Yin, S. Osher, D. Goldfarb, and J. Darbon, 
! "Bregman Iterative Algorithms for L1-Minimization with Applications to Compressed Sensing,"
! SIAM J. Img. Sci. 1, 143 (2008).
!
! Fixed-point continuation is described in E. T. Hale, W. Yin, and Y. Zhang, 
! "A fixed-point continuation method for l1-regularized minimization with applications 
! to compressed sensing," CAAM TR07-07, Rice University (2007).
!
! The split Bregman algorithm is described in T. Goldstein and S. Osher,
! "The split Bregman method for L1 regularized problems",
! SIAM Journal on Imaging Sciences Volume 2 Issue 2, Pages 323-343 (2009).
!
! The conjugae gradient algorithm is described in S. Boyd and L. Vandenberghe,
! "Convex Optimization" (Cambridge University Press, 2004).
!
! Coded by Vidvuds Ozolins (UCLA Materials Science & Eng).
! Last modified: August 5, 2012
!
public BregmanFPC
public SplitBregman
!public CGmin
!public Shrink

!
! The same as above, but takes advantage of sparse sensing matrices.
!
public SparseBregmanFPC
public SparseSplitBregman
!public SparseCGmin

! with right preconditioner
public BregmanRPrecond2
public BregmanRPrecond

!
! Private parameters for the fixed-point continuation (FPC) algorithm.
!
! MaxFPCit - maximum number of FPC iterations
! nPrint   - number of FPC its to print out gradients and residuals;
!            specify nPrint > MaxFPCit if output is not desired
! FPCgtol  - tolerance for gradient in FPC iterations; the algorithm will iterate until
!            FPCgtol > ||grad||_infty / mu - 1d0
! FPCxtol  - tolerance for the solution; the algorithm will iterate until
!            FPCxtol > ||u-uprev||_2 / ||uprev||_1
!
integer, private                  :: MaxFPCit = 60000
integer, private                  :: nPrint = 1000
double precision, private         :: FPCgtol = 1d-1
double precision, private         :: FPCxtol = 1d-4

!
! Tolerance factor for the CG algorithm: CG stops if ||grad|| <  CGtol*||deltaU||/||U||
!
double precision, private         :: CGtol = 1d-1

CONTAINS

SUBROUTINE BregmanFPC(nBregIt, tol, tau, mu, N, M, A, b, u) BIND(C,NAME="BregmanFPC")
  ! 
  ! Performs Bregman iterations using fixed-point continuation (FPC)
  ! as the solver for the mixed L1/L2 problem:
  !
  !    u = arg min_u { mu ||u||_1 + 1/2 ||Au-b||_2^2 }
  !
  ! The Bregman algorithm is described in W. Yin, S. Osher, D. Goldfarb, and J. Darbon, 
  ! "Bregman Iterative Algorithms for L1-Minimization with Applications to Compressed Sensing,"
  ! SIAM J. Img. Sci. 1, 143 (2008).
  !
  ! Input parameters:
  !   nBregIt - Number of outer Bregman loops (usually 5)
  !   tol     - Required tolerance for the residual. The algorithm stops when
  !             tol > ||Au-b||_2 / ||b||_2
  !   tau     - step size for the FPC iterations;
  !             tau = MIN(1.999d0, -1.665d0 * M/N + 2.665d0)
  !   mu      - weight for the L1 norm of the solution
  !   N   - number of unknown expansion coefficients [ =length(u) ]
  !   M       - number of measurements [ =length(b) ]
  !   A       - sensing matrix of dimensions (M,N); the maximum eigenvalue of A.AT should be <= 1
  !             this can be enforced by dividing both A and b with lambda=Sqrt(Max(Eigenvalues(A.Transpose(A))))
  !   b       - array with measured signal values of length M 
  !   u       - solution array of length N
  !
  ! Output:
  !   u       - solution
  !

  integer, intent(in) :: nBregIt, N, M
  double precision, intent(in) :: tol, tau, mu
  double precision,intent(in) :: A(M,N), b(M)
  double precision,intent(out) :: u(N)
!f2py intent(hide) :: M, N
!f2py intent(in) :: A, b, nBregIt, tol, tau, mu
!f2py intent(out) :: u

  double precision, allocatable:: uprev(:), grad(:), residual(:), x(:), bp(:)

  integer j, k
  double precision crit1, crit2

  ALLOCATE( grad(N), uprev(N), bp(M), x(M), residual(M) )

  bp = b

  do k = 1, nBregIt

     do j = 1, MaxFPCit
        uprev = u

        x = MATMUL(A, u) - bp
        grad = MATMUL(TRANSPOSE(A), x)
        u = u - tau * grad

        !        print *, ' norm of residual = ', SQRT(dot_product(x,x))
        !        print *, ' norm of gradient = ', SQRT(dot_product(grad,grad))
        !        print *, ' norm of u = ', SQRT(dot_product(u,u))

        call Shrink(u, N, mu*tau)
        !        print *, ' norm of u after shrinkage = ', SQRT(dot_product(u,u))

        crit1 = MAXVAL(ABS(grad)) / mu - 1d0
        crit2 = SQRT(DOT_PRODUCT(u-uprev,u-uprev)) / MAX(SUM(ABS(uprev)),1d0)

        if( MOD(j,nPrint) == 0 ) print '(i6,2f14.6)', j, crit1, crit2
        if( crit1 < FPCgtol .and. crit2 < FPCxtol ) exit
     end do

     if ( j >= MaxFPCit ) then
         print '(2(a,f12.8))', " Unconverged FPC: ||deltaU||/||U|| =", crit2, " ||grad||/mu - 1 = ", crit1 
     end if

     residual = MATMUL(A,u) - b

     crit1 =  SQRT(DOT_PRODUCT(residual,residual)) / SQRT(DOT_PRODUCT(b,b))
     print '(a,i2,a,i6,a,f12.6)', ' Finished Bregman loop #',k,' after ', &
          j,' FPC steps, ||Au-b||/||b|| = ', crit1

     if( crit1 <= tol ) then
        exit
     else
        bp = bp - residual
     end if
  end do

  DEALLOCATE( grad, uprev, bp, x, residual )

  return
end SUBROUTINE BregmanFPC


SUBROUTINE SparseBregmanFPC(nBregIt, tol, tau, mu, N, M, nA, Aval, iAx, b, u)
  ! 
  ! Performs Bregman iterations using fixed-point continuation (FPC)
  ! as the solver for the mixed L1/L2 problem:
  !
  !    u = arg min_u { mu ||u||_1 + 1/2 ||Au-b||_2^2 }
  !
  ! This version uses a sparse sensing matrix A.
  !
  ! The Bregman algorithm is described in W. Yin, S. Osher, D. Goldfarb, and J. Darbon, 
  ! "Bregman Iterative Algorithms for L1-Minimization with Applications to Compressed Sensing,"
  ! SIAM J. Img. Sci. 1, 143 (2008).
  !
  ! Input parameters:
  !   nBregIt - Number of outer Bregman loops (usually 5)
  !   tol     - Required tolerance for the residual. The algorithm stops when
  !             tol > ||Au-b||_2 / ||b||_2
  !   tau     - step size for the FPC iterations;
  !             tau = MIN(1.999d0, -1.665d0 * dble(M) / dble(N) + 2.665d0)
  !   mu      - weight for the L1 norm of the solution
  !   N       - number of unknown expansion coefficients [ =length(u) ]
  !   M       - number of measurements [ =length(b) ]
  !   nA      - number of nonzero values of the sensing matrix A
  !   Aval    - nonzero values of the sensing matrix; the maximum eigenvalue of A.AT should be <= 1; 
  !             this can be enforced by dividing both A and b with lambda=Sqrt(Max(Eigenvalues(A.Transpose(A))))
  !   iAx     - indices of the nonzero elements of A; values of iAx are within the range (1:M,1:N)
  !   b       - array with measured signal values (of length M)
  !   u       - solution array (of length N)
  !
  ! Output:
  !   u       - solution
  !
  integer, intent(in) :: nBregIt, N, M, nA, iAx(nA,2)
  double precision, intent(in) :: tol, tau, mu
  double precision,intent(in) :: Aval(nA), b(M)
  double precision,intent(out) :: u(N)

  double precision, allocatable:: uprev(:), grad(:), residual(:), x(:), bp(:)

  integer i, j, k
  double precision crit1, crit2

  ALLOCATE( grad(N), uprev(N), bp(M), x(M), residual(M) )

  bp = b

  do k = 1, nBregIt

     do j = 1, MaxFPCit
        uprev = u

#ifdef IntelMKL
        !
        ! Use BLAS calls if Intel MKL is available.
        !
        call mkl_dcoogemv('N', M, Aval, iAx(:,1), iAx(:,2), nA, u, x)
        x = x - bp
        call mkl_dcoogemv('N', N, Aval, iAx(:,2), iAx(:,1), nA, x, grad)
#else
        ! 
        ! Manual matrix-vector multiplies
        !
        x = 0d0
        do i = 1, nA
           x(iAx(i,1)) = x(iAx(i,1)) + Aval(i) * u(iAx(i,2))
        end do
        x = x - bp

        grad = 0d0
        do i = 1, nA
           grad(iAx(i,2)) = grad(iAx(i,2)) + Aval(i) * x(iAx(i,1))
        end do
#endif

        u = u - tau * grad

        !        print *, ' norm of residual = ', SQRT(dot_product(x,x))
        !        print *, ' norm of gradient = ', SQRT(dot_product(grad,grad))
        !        print *, ' norm of u = ', SQRT(dot_product(u,u))

        call Shrink(u, N, mu*tau)
        !        print *, ' norm of u after shrinkage = ', SQRT(dot_product(u,u))

        crit1 = MAXVAL(ABS(grad)) / mu - 1d0
        crit2 = SQRT(DOT_PRODUCT(u-uprev,u-uprev)) / MAX(SUM(ABS(uprev)),1d0)

        if( MOD(j,nPrint) == 0 ) print '(i6,2f14.6)', j, crit1, crit2
        if( crit1 < FPCgtol .and. crit2 < FPCxtol ) exit
     end do

     if ( j>= MaxFPCit ) then
        print '(2(a,f12.8))', " Unconverged FPC: ||deltaU||/||U|| =", crit2, " ||grad||/mu - 1 = ", crit1 
     end if

#ifdef IntelMKL
     !
     ! Use BLAS calls if Intel MKL is available.
     !
     call mkl_dcoogemv('N', M, Aval, iAx(:,1), iAx(:,2), nA, u, residual)
#else
     ! 
     ! Manual matrix-vector multiply
     !
     residual = 0d0
     do i = 1, nA
        residual(iAx(i,1)) = residual(iAx(i,1)) + Aval(i) * u(iAx(i,2))
     end do
#endif

     residual = residual - b

     crit1 =  SQRT(DOT_PRODUCT(residual,residual)) / SQRT(DOT_PRODUCT(b,b))
     print '(a,i2,a,i6,a,f12.6)', ' Finished Bregman loop #',k,' after ', &
          j,' FPC steps, ||Au-b||/||b|| = ', crit1

     if( crit1 <= tol ) then
        exit
     else
        bp = bp - residual
     end if
  end do

  DEALLOCATE( grad, uprev, bp, x, residual )

  return
end SUBROUTINE SparseBregmanFPC


SUBROUTINE SplitBregman(MaxIt, tol, mu, lambda, N, M, A, f, u)
  ! 
  ! Performs split Bregman iterations using conjugate gradients (CG) for the
  ! L2 minimizations and shrinkage for L1 minimizations.
  !
  !    u = arg min_u { ||u||_1 + mu/2 ||Au-f||_2^2 }
  !
  ! The algorithm is described in T. Goldstein and S. Osher,
  ! "The split Bregman method for L1 regularized problems",
  ! SIAM Journal on Imaging Sciences Volume 2 Issue 2, Pages 323-343 (2009).
  !
  ! Input parameters:
  !   MaxIt   - Number of outer split Bregman loops
  !   tol     - Required tolerance for the residual. The algorithm stops when
  !             tol > ||Au-f||_2 / ||f||_2
  !   mu      - weight for the L1 norm of the solution
  !   lambda  - weight for the split constraint (affects speed of convergence, not the result)
  !   N       - number of unknown expansion coefficients [ =length(u) ]
  !   M       - number of measurements [ =length(f) ]
  !   A       - sensing matrix A
  !   f       - array with measured signal values (of length M)
  !   u       - solution array (of length N)
  !
  ! Output:
  !   u       - solution
  !
  integer, intent(in)             :: MaxIt, N, M
  double precision, intent(in)    :: tol, mu, lambda
  double precision,intent(in)     :: A(M,N), f(M)
  double precision,intent(out)  :: u(N)
!f2py intent(hide) :: M, N
!f2py intent(in) :: A, f
!f2py intent(out) :: u

  integer k, MaxCGit
  double precision crit1, crit2
  double precision, allocatable:: uprev(:), b(:), d(:), bp(:)

  MaxCGit = MAX(10,N/2)
  crit1 = 1d0
  crit2 = 1d0

  ALLOCATE( uprev(N), b(N), d(N), bp(N) )

  b = 0d0
  d = 0d0

  do k = 1, MaxIt

     uprev = u

     bp = d - b
     call CGmin(N, M, A, f, bp, mu, lambda, MaxCGit, CGtol*crit1, u)

     d = b + u
     call Shrink(d, N, 1d0/lambda)

     crit1 = SQRT(dot_product(u-uprev,u-uprev)/dot_product(u,u))
     crit2 = SQRT(dot_product(u-d,u-d)/dot_product(u,u))
     if (mod(k, 5)==1) then
       print '(a,i3,2(a,f10.6))', ' SplitBregman: it=', k, ', |dU|/|U| = ', crit1, ', |d-U|/|U| = ', crit2
     end if
     if ( crit1 <= tol ) exit

     b = b + u - d

  end do

  if ( crit1 > tol ) then
     print '(2(a,f10.6))', ' Did not reach prescribed accuracy in SplitBregman: |dU|/|U| = ', crit1, &
          ', ||d-U||/||U|| = ', crit2
  end if

  DEALLOCATE( uprev, b, d, bp )
  return

end SUBROUTINE SplitBregman


SUBROUTINE SparseSplitBregman(MaxIt, tol, mu, lambda, N, M, nA, Aval, iAx, f, u)
  ! 
  ! Performs split Bregman iterations using conjugate gradients (CG) for the
  ! L2 minimizations and shrinkage for L1 minimizations.
  !
  !    u = arg min_u { ||u||_1 + mu/2 ||Au-f||_2^2 }
  !
  ! This version uses sparse sensing matrices A.
  !
  ! The algorithm is described in T. Goldstein and S. Osher,
  ! "The split Bregman method for L1 regularized problems",
  ! SIAM Journal on Imaging Sciences Volume 2 Issue 2, Pages 323-343 (2009).
  !
  ! Input parameters:
  !   MaxIt   - Number of outer split Bregman loops
  !   tol     - Required tolerance for the residual. The algorithm stops when
  !             tol > ||Au-f||_2 / ||f||_2
  !   mu      - weight for the L1 norm of the solution
  !   lambda  - weight for the split constraint (affects speed of convergence, not the result)
  !   N       - number of unknown expansion coefficients [ =length(u) ]
  !   M       - number of measurements [ =length(f) ]
  !   nA      - number of nonzero values of the sensing matrix A
  !   Aval    - nonzero values of the sensing matrix; the maximum eigenvalue of A.AT should be <= 1; 
  !             this can be enforced by dividing both A and f with lambda=Sqrt(Max(Eigenvalues(A.Transpose(A))))
  !   iAx     - indices of the nonzero elements of A; values of iAx are within the range (1:M,1:N)
  !   f       - array with measured signal values (of length M)
  !   u       - solution array (of length N)
  !
  ! Output:
  !   u       - solution
  !
  integer, intent(in) :: MaxIt, N, M, nA, iAx(nA,2)
  double precision, intent(in) :: tol, mu, lambda
  double precision,intent(in) :: Aval(nA), f(M)
  double precision,intent(out) :: u(N)

  integer k, MaxCGit
  double precision crit1, crit2
  double precision, allocatable:: uprev(:), b(:), d(:), bp(:)

  MaxCGit = MAX(10,N/5)
  crit1 = 1d0
  crit2 = 1d0

  ALLOCATE( uprev(N), b(N), d(N), bp(N) )
  b = 0d0
  d = 0d0

  do k = 1, MaxIt

     uprev = u

     bp = d - b
     call SparseCGmin(N, M, nA, Aval, iAx, f, bp, mu, lambda, MaxCGit, CGtol*crit1, u)

     d = b + u
     call Shrink(d, N, 1d0/lambda)

     crit1 = SQRT(DOT_PRODUCT(u-uprev,u-uprev)/DOT_PRODUCT(u,u))
     crit2 = SQRT(DOT_PRODUCT(u-d,u-d)/DOT_PRODUCT(u,u))
     print '(a,i3,2(a,f10.6))', ' SplitBregman: it=', k, ', ||deltaU||/||U|| = ', crit1, &
          ', ||d-U||/||U|| = ', crit2
     if ( crit1 <= tol ) exit

     b = b + u - d

  end do

  if ( crit1 > tol ) then
     print '(2(a,f10.6))', ' Did not reach prescribed accuracy in SplitBregman: ||deltaU||/||U|| = ', crit1, &
          ', ||d-U||/||U|| = ', crit2
  end if

  DEALLOCATE( uprev, b, d, bp )
  return

end SUBROUTINE SparseSplitBregman


SUBROUTINE CGmin(N, M, A, f, b, mu, lambda, MaxIt, gtol, u)
  !
  ! Conjugate gradient routine to perform L2-based minimization of 
  !
  !     min_u { mu/2 ||Au-f||_2^2 + lambda/2 ||b-u||_2^2 }
  !
  ! Algorithm is described in S. Boyd and L. Vandenberghe,
  ! "Convex Optimization" (Cambridge University Press, 2004).
  !
  ! Inut parameters:
  !    A    - sensing matrix of dimensions (M,N)
  !    M    - number of measurements
  !    N    - number of expansion coefficients
  !    f    - values of measurements
  !    b    - vector enforcing the split-off L1 constraint
  !    mu   - weight of the L2 constraint on Au=f
  !    lambda - weight of the split-off constraint
  !    MaxIt  - max. number of CG iterations
  !    gtol - tolerance for the gradient; exit if gtol > ||grad||_2
  !    u    - starting guess for the solution
  !
  ! Output parameter:
  !    u    - converged solution
  !
  integer, intent(in) :: MaxIt, N, M
  double precision, intent(in) :: gtol, mu, lambda
  double precision,intent(in) :: A(M,N), f(M), b(N)
  double precision,intent(out) :: u(N)

  integer k
  double precision, allocatable:: p(:), r(:), rp(:), x(:)
  double precision beta, alpha, delta, deltaprev

  ALLOCATE( p(N), r(N), rp(N), x(M) )

  x = MATMUL(A,u) - f
  r = -( mu * MATMUL(TRANSPOSE(A),x) - lambda * (b - u) )

  p = r
  delta = DOT_PRODUCT(r,r)
  !    print '(a,f16.4)', ' initial norm of gradient = ', SQRT(delta)

  do k = 1, MaxIt

     x = MATMUL(A,p)
     rp = mu * MATMUL(TRANSPOSE(A),x) + lambda * p

     alpha = delta / DOT_PRODUCT(p,rp)

     !     print '(a,f16.6)', ' gradient check = ', DOT_PRODUCT(r,r - alpha * rp)

     u = u + alpha * p
     r = r - alpha * rp

     deltaprev = delta
     delta = DOT_PRODUCT(r,r)

     !     print '(a,f16.4)', ' norm of gradient = ', SQRT(delta)

     if( SQRT(delta) < gtol ) then
        !     print '(a,i4,a)', ' Finished after ', k, ' CG iterations.'
        exit
     end if

     beta = delta / deltaprev
     p = beta * p + r

  end do

  if ( SQRT(delta) > gtol ) then
     print '(a,f12.6)', ' Unconverged gradient after CGmin = ', SQRT(delta)
  end if

  DEALLOCATE( p, r, rp, x )

  return
end SUBROUTINE CGmin



SUBROUTINE SparseCGmin(N, M, nA, Aval, iAx, f, b, mu, lambda, MaxIt, gtol, u)
  !
  ! Conjugate gradient routine to perform L2-based minimization of 
  !
  !     min_u { mu/2 ||Au-f||_2^2 + lambda/2 ||b-u||_2^2 }
  !
  ! This version uses a sparse sensing matrix A.
  !
  ! Algorithm is described in S. Boyd and L. Vandenberghe,
  ! "Convex Optimization" (Cambridge University Press, 2004).
  !
  ! Inut parameters:
  !    M    - number of measurements
  !    N    - number of expansion coefficients
  !    nA   - number of nonzero values of the sensing matrix A
  !    Aval - nonzero values of the sensing matrix
  !    iAx    - indices of the nonzero elements of A; values of iAx are within the range (1:M,1:N)
  !    f    - values of measurements
  !    b    - vector enforcing the split-off L1 constraint
  !    mu   - weight of the L2 constraint on Au=f
  !    lambda - weight of the split-off constraint
  !    MaxIt  - max. number of CG iterations
  !    gtol - tolerance for the gradient; exit if gtol > ||grad||_2
  !    u    - starting guess for the solution
  !
  ! Output parameter:
  !    u    - converged solution
  !
  integer, intent(in)            :: MaxIt, N, M, nA, iAx(nA,2)
  double precision, intent(in)   :: gtol, mu, lambda
  double precision,intent(in)    :: Aval(nA), f(M), b(N)
  double precision,intent(out) :: u(N)

  integer i, k
  double precision, allocatable:: p(:), r(:), rp(:), x(:)
  double precision beta, alpha, delta, deltaprev

  ALLOCATE( p(N), r(N), rp(N), x(M) )

#ifdef IntelMKL
  !
  ! Use BLAS calls if Intel MKL is available.
  !
  call mkl_dcoogemv('N', M, Aval, iAx(:,1), iAx(:,2), nA, u, x)
  x = x - f
  call mkl_dcoogemv('N', N, Aval, iAx(:,2), iAx(:,1), nA, x, r)
  r = -( mu * r - lambda * (b - u) )
#else
  !
  ! Manual matrix-vector multiplies.
  !
  x = 0d0
  do i = 1, nA
     x(iAx(i,1)) = x(iAx(i,1)) + Aval(i) * u(iAx(i,2))
  end do
  x = x - f

  r = 0d0
  do i = 1, nA
     r(iAx(i,2)) = r(iAx(i,2)) + Aval(i) * x(iAx(i,1))
  end do
  r = -( mu * r - lambda * (b - u) )
#endif

  p = r
  delta = DOT_PRODUCT(r,r)
  !      print '(a,f16.4)', ' norm of gradient = ', SQRT(delta)

  do k = 1, MaxIt

#ifdef IntelMKL
     !
     ! Use BLAS calls if Intel MKL is available.
     !
     call mkl_dcoogemv('N', M, Aval, iAx(1:nA,1), iAx(1:nA,2), nA, p, x)
     call mkl_dcoogemv('N', N, Aval, iAx(1:nA,2), iAx(1:nA,1), nA, x, rp)
#else
     !
     ! Manual matrix-vector multiplies.
     !
     x = 0d0
     do i = 1, nA
        x(iAx(i,1)) = x(iAx(i,1)) + Aval(i) * p(iAx(i,2))
     end do

     rp = 0d0
     do i = 1, nA
        rp(iAx(i,2)) = rp(iAx(i,2)) + Aval(i) * x(iAx(i,1))
     end do
#endif

     rp = mu * rp + lambda * p

     alpha = delta / DOT_PRODUCT(p,rp)
     !     print '(a,f16.6)', ' gradient check = ', dot_product(r,r - alpha * rp)

     u = u + alpha * p
     r = r - alpha * rp

     deltaprev = delta
     delta = DOT_PRODUCT(r,r)
     !     print '(a,f16.4)', ' norm of gradient = ', SQRT(delta)

     if( SQRT(delta) < gtol ) then
        !     print '(a,i4,a)', ' Finished after ', k, ' CG iterations.'
        exit
     end if

     beta = delta / deltaprev
     p = beta * p + r

  end do

  if ( SQRT(delta) > gtol ) then
     print '(a,f12.6)', ' Unconverged gradient after CGmin = ', SQRT(delta)
  end if

  DEALLOCATE( p, r, rp, x )

  return
end SUBROUTINE SparseCGmin



subroutine Shrink( u, N, alpha )
  !
  ! Defines L1-based shrinkage.
  !
  integer, intent(in)  :: N
  double precision, intent(in) :: alpha
  double precision, intent(out) :: u(N)

  where (u > 0d0) 
     u = MAX(ABS(u) - alpha, 0d0)
  elsewhere
     u = - MAX(ABS(u) - alpha, 0d0)
  end where

  return
end subroutine Shrink



!--------------------------------------------
!Right preconditioner
subroutine RightPreconditioner(M, N, A, mu, lambda, Vm)
  integer, intent(in) :: M, N
  double precision, intent(in) :: mu, lambda
  double precision, intent(in) :: A(M,N)
  double precision, intent(out) :: Vm(N,N)
  
  integer i, LDA, lwork, info, ldvl, ldvr, Debug
  double precision tol, beta
  double precision, allocatable :: wr(:), wi(:), work(:), vl(:,:), vr(:,:)
  double precision, allocatable :: mt(:,:)

  tol=1d-8
  beta=0
  Debug=0
  LDA=N
  lwork=10*N
  ldvl=N
  ldvr=N
  allocate(mt(N,N), wr(N), wi(N), work(max(1,lwork)), vl(ldvl,N), vr(ldvr,N))

  
  !print *, "A:", size(A)
  !mt=mu*matmul(transpose(A), A)
  !much faster
  call dgemm('N', 'N', N, N, M, mu, transpose(A), N, A, M, beta, mt, N)
  !print *, "N:", N

  do i=1, N
     mt(i,i)=mt(i,i)+lambda
  end do
  mt=(mt+transpose(mt))/2d0

  if (Debug > 0) then
     print *, 'mt: '
     print *, shape(mt)
     !print '(4(f8.4,4x))', mt
  end if

  !print *, "---Begin deysv:"
  call dsyev('V', 'L', N, mt, LDA, wr, work, lwork, info)
  !print *, "---End deysv"
  !  call dgeev('N', 'V', N, mt, LDA, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)


  if (Debug >0) then
     print *, 'wr: '
     print '(4(f8.4, 4x))', wr
     print *, 'mt: '
     do i=1, N
        print '(4f16.8)', mt(:,i)
     end do
  end if

  Vm=0d0
  do i=1, N
     Vm(i,i)=1d0/sqrt(wr(i))
  end do
  Vm=matmul(mt, matmul(Vm, transpose(mt)))

  deallocate(wr, wi, work, vl, vr, mt)

  return
end subroutine RightPreconditioner


!--------------------------------------------
!conjugate gradient method for splitbregman

subroutine CGmin2(M, N, A, f, C, b, u, mu, lambda, MaxIt, gtol)
  !
  ! Conjugate gradient routine to perform L2-based minimization of 
  !
  !     min_u { mu/2 ||Au-f||_2^2 + lambda/2 ||b-Cu||_2^2 }
  !
  ! Algorithm is described in S. Boyd and L. Vandenberghe,
  ! "Convex Optimization" (Cambridge University Press, 2004).
  !
  ! Inut parameters:
  !    A    - sensing matrix of dimensions (M,N)
  !    C    - L1 norm matrix of dimensions (N',N)
  !    M    - number of measurements
  !    N    - number of expansion coefficients
  !    N'   - 
  !    f    - values of measurements
  !    b    - vector enforcing the split-off L1 constraint
  !    mu   - weight of the L2 constraint on Au=f
  !    lambda - weight of the split-off constraint
  !    MaxIt  - max. number of CG iterations
  !    gtol - tolerance for the gradient; exit if gtol > ||grad||_2
  !    u    - starting guess for the solution
  !
  ! Output parameter:
  !    u    - converged solution
  !
  integer, intent(in) :: MaxIt, M, N
  double precision, intent(in) :: gtol, mu, lambda
  double precision, intent(in) :: A(M,N), f(M), C(N,N), b(N)
  double precision, intent(out) :: u(N)

  integer k
  double precision, allocatable :: p(:), r(:), rp(:), x(:)
  double precision beta, alpha, delta, deltaprev

  allocate(p(N), r(N), rp(N), x(N))
  p=0d0
  r=0d0
  rp=0d0

!  x=matmul(A,u)-f
  r=-mu*matmul(transpose(A),matmul(A,u)-f)-lambda*matmul(transpose(C),matmul(C,u)-b)
  p=r
  delta=dot_product(r, r)

  do k=1, MaxIt
     rp=mu*matmul(transpose(A),matmul(A, p))+lambda*matmul(transpose(C), matmul(C, p))
     alpha=delta/dot_product(p, rp)

     u=u+alpha*p
     r=r-alpha*rp

     deltaprev=delta
     delta=dot_product(r, r)

     if (sqrt(delta)<gtol) then
        !print *
        !print *, 'Finished: k=', k
        exit
     end if
     
     beta=delta/deltaprev
     p=beta*p+r
  end do

  if (sqrt(delta) > gtol) then
     print '(a,f12.6)', 'Unconverged gradient after CGmin2 = ', sqrt(delta)
  end if
  
  deallocate(p, r, rp, x)

  return

end subroutine CGmin2


!--------------------------------------------
!norm

!--------------------------------------------
!splitbregman2
subroutine splitbregman2(MaxIt, tol, mu, lambda, N, M, A, f, C, u)
  integer, intent(in) :: MaxIt, M, N
  double precision, intent(in) :: tol, mu, lambda
  double precision, intent(in) :: A(M,N), f(M), C(N,N)
  double precision, intent(out) :: u(N)

  integer k, MaxCGit, nPrint, converged
  double precision crit1, crit2, CGtol
  double precision, allocatable :: uprev(:), b(:), d(:)
  MaxCGit=min(100, N)
  CGtol=1d-3
  crit1=1d0
  crit2=1d0
  nPrint=int(MaxIt/10)
  converged=0d0
  
  allocate(uprev(N), b(N), d(N))
  uprev=0d0
  b=0d0
  d=0d0

  do k=1, MaxIt
     uprev=u
     call CGmin2(M, N, A, f, C, d-b, u, mu, lambda, MaxCGit, crit1*CGtol)
     d=b+matmul(C, u)
     call Shrink(d, N, 1/lambda)
     crit1=sqrt(dot_product(uprev-u, uprev-u)/dot_product(u,u))
     crit2=sqrt(dot_product(matmul(C,u)-d, matmul(C,u)-d)/dot_product(u,u))

     if (mod(k,nPrint)==0) then
        !print '(2(f12.6))', "Bregman: ||du||/||u|| = ",crit1, ", ||Cu-d||/||u|| = ",crit2
        print 100, crit1, crit2
     end if
  100 format("Bregman: ||du||/||u|| = ", f12.6, ", ||Cu-d||/||u|| = ", f12.6)
     if (crit1<=tol .and. crit2<=tol) then
        converged=1
        exit
     end if

     b=b+matmul(C,u)-d
  end do
  
  if (converged<=0) then
     print *, "SplitBregman2 stopped with ||du||/||u|| = ", sqrt(dot_product(uprev-u, uprev-u)/dot_product(u,u)) 
     print *, "Finished ",k," split Bregman iteractions in total!"
  end if
  
  deallocate(uprev, b, d)
  
  return
end subroutine splitbregman2


!--------------------------------------------
!Bregman with right preconditioner, with constraint "C"
! hence the suffix "2"
subroutine BregmanRPrecond2(maxIter, M, N, A, b, mu, lambda, tol, C, u)
  integer, intent(in) :: maxIter, M, N
  double precision, intent(in) :: mu, lambda, tol
  double precision, intent(in) :: A(M,N), b(M), C(N,N)
  double precision, intent(out) :: u(N)

  integer i, j
  double precision sum
  double precision, allocatable :: Bm(:,:), Vm(:,:), BmVm(:,:)
  
  allocate(Bm(M,N), Vm(N,N), BmVm(M,N))

  !print *, "-Begin BregmanRPrecond2:"
  
  do i=1, N
     do j=1, N
        sum=sum+C(i,j)*C(i,j)
     end do
  end do
  
  if (sum<=1d-12) then
     Bm=A
  else
     Bm=matmul(A,C)
  end if

  !print *, "--Begin RightPreconditioner:"
  call RightPreconditioner(M, N, A, mu, lambda, Vm)
  !print *, "--End RightPreconditioner."

  BmVm=matmul(Bm, Vm)
  call splitbregman2(maxIter, tol, mu, lambda, N, M, BmVm, b, Vm, u)

  if (sum<=1d-12) then
     u=matmul(Vm, u)
  else
     u=matmul(C, matmul(Vm, u))
  end if

  deallocate(Bm, Vm, BmVm)

  return
end subroutine BregmanRPrecond2

!--------------------------------------------
!Bregman with right preconditioner
subroutine BregmanRPrecond(maxIter, M, N, A, b, mu, lambda, tol, u)
  integer, intent(in) :: maxIter, M, N
  double precision, intent(in) :: mu, lambda, tol
  double precision, intent(in) :: A(M,N), b(M)
  double precision, intent(out) :: u(N)

  integer i, j
  double precision, allocatable :: Vm(:,:), BmVm(:,:)

!f2py intent(hide),depend(A) :: M=shape(M,0), N=shape(M,1)
!f2py intent(in) :: A, b
!f2py intent(out) :: u
  allocate(Vm(N,N), BmVm(M,N))

  ! print *, "--Begin RightPreconditioner:"
  call RightPreconditioner(M, N, A, mu, lambda, Vm)
  !print *, "--End RightPreconditioner."

  BmVm=matmul(A, Vm)
  call splitbregman(maxIter, tol, mu, lambda, N, M, BmVm, b, u)

  u=matmul(Vm, u)

  deallocate(Vm, BmVm)

  return
end subroutine BregmanRPrecond


!--------------------------------------------
!random sampling
subroutine ransam(X, A, B, N, K)
  integer, intent(in) :: N, K
  integer, intent(in) :: X(N)
  integer, intent(inout) :: A(K)
  integer, intent(inout) :: B(N-K)
  
  integer Neff, index, Debug, isin, iB
  integer Xtmp(N)
  double precision val
  integer i, j
  
  if (K>N) then
     print *, 'Pool too small!'
  end if

  Xtmp=X
  if (dot_product(Xtmp, Xtmp) < 1d-12) then
     do i=1, N
        Xtmp(i)=i
     end do
  end if
  A=0d0
  Debug=0d0

! random_seed
  call random_seed()
! check random number generator
  if (Debug>0) then
     do i=1, 10
        call random_number(val)
        print *, 'val:', val
     end do
     print *
  end if

  do i=1, K
     call random_number(val)
     Neff=N-i+1
     index=int(Neff*val)+1

     if (Debug>0) then
        print *, 'val: ', val
        print *, 'index: ', index
        print *, 'Selected: ', X(index)
        print *, '---------------------------'
     end if

     A(i)=Xtmp(index)
     do j=index, N-1
        Xtmp(j)=Xtmp(j+1)
     end do
  end do
  
  if(minval(A)<1 .and. maxval(A)>N) then
     print *, 'Something wrong with random sampling!'
  end if
   
  iB=1
  do i=1, N
     isin=0
     do j=1, K
        if (X(i)==A(j)) then
           isin=1
        end if
     end do
     if (isin==0) then
        B(iB)=X(i)
        iB=iB+1
     end if
  end do

end subroutine ransam


!csfitting Amat(Mtot, Ntot); En(Mtot); Wt(Moto)
subroutine csfit(Amat, En, Wt, Mtot, Ntot, solfinal)
  integer i, j, k
  integer Mtot, Ntot, nline, nout, nin, nmu, nsub
  
  double precision lambda
  !in
  double precision Amat(Mtot,Ntot), En(Mtot), Wt(Mtot)
  double precision, allocatable :: mulist(:), theC(:,:)
  integer, allocatable :: Inlist(:), Outlist(:), X(:)

  integer nsubin, nsubout, index
  double precision mu, relerr, errtmp
  integer, allocatable :: Insublist(:), Outsublist(:)
  double precision, allocatable :: Asub(:,:), Esub(:), Csub(:,:), sol(:)
  double precision, allocatable :: Aout(:,:), Eout(:), Etmp(:)
  !out
  double precision, intent(out) :: solfinal(Ntot)

!f2py intent(hide) :: Mtot, Ntot
!f2py intent(in) :: Amat, En, Wt
!f2py intent(out) :: solfinal


  !main program
  !select training set
  nsub=2
  nout=int(Mtot*0.12)
  nin=Mtot-nout
  allocate(Inlist(nin), Outlist(nout), X(Mtot))
  do i=1, Mtot
     X(i)=i
  end do
  call ransam(X, Inlist, Outlist, Mtot, nin)
  !print *, 'Shape of Inlist: ', shape(Inlist)
  print *, 'Selected for fitting: ', shape(Inlist)
  

  !set mu and lambda
  lambda=3d0
  nmu=5
  allocate(mulist(nmu))
  do i=1, nmu
     mulist(i)=dble(3d0)**dble(i-2d0)
  end do
  print *
  print *, 'lambda: ', lambda
  print *, 'mulist: '
  print '(5(f8.4, 4x))', mulist

  !set theC(selected fitting)
  allocate(theC(Ntot, Ntot))
  theC=0d0
  do i=1, Ntot
     theC(i, i)=1d0
  end do
  print *
  print *, "theC: (Identity Matrix)"
  print *, "Shape: "
  print *, shape(theC)
  print * 

  !start main loop
  !BregmanRPrecond2(M, N, A, b, mu, lambda, C, u)
  nsubin=int(nin*0.6)
  nsubout=nin-nsubin
  allocate(Insublist(nsubin), Outsublist(nsubout))
  allocate(Asub(nsubin, Ntot), Esub(nsubin), Csub(Ntot, Ntot), sol(Ntot))
  allocate(Aout(nout, Ntot), Eout(nout), Etmp(nout))
  Csub=0d0
  sol=0d0
  solfinal=0d0
  errtmp=1d0

!  call ransam(Inlist, Insublist, Outsublist, nin, nsubin) 
  print *
  print *, 'Split Bregman begin...'
  print *
  do i=1, nmu  ! loop over mulist
     mu=mulist(i)
     print 202, mu
 202 format("Mu = ", f8.4, '------------------------------------------------------')
     do j=1, nsub  ! loop over subset list
        
        sol=0d0
        call ransam(Inlist, Insublist, Outsublist, nin, nsubin)
        print *, "Selected forces: ", nsubin
        do k=1, nsubin
           index=Insublist(k)
           Asub(k,:)=Amat(index,:)
           Esub(k)=En(index)
        end do

        call BregmanRPrecond2(200, nsubin, Ntot, Asub, Esub, mu, lambda, 1D-4, Csub, sol)
        do k=1, nout
           index=Outlist(k)
           Aout(k,:)=Amat(index,:)
           Eout(k)=En(index)
        end do
        
        Etmp=matmul(Aout, matmul(theC, sol))
        relerr=sqrt(dot_product(Etmp-Eout, Etmp-Eout)/dot_product(Eout,Eout))
        print *, 'relative error: ', relerr
        
        !select the sol with minimum relerr
        if (relerr < errtmp) then
           errtmp=relerr
           solfinal=sol
        end if
        print *

     end do !end loop over subset

     print *, 'Finished mu: ', mu
     print *

  end do !end loop over mulist

  print *, '------------------------------------------------'
  print *, 'Minimum relative error: ', errtmp
  print *, 'solfinal: '
  print *, 'shape: ', shape(solfinal)
  open(33,file='sol',form='formatted',status='unknown')
  do i=1, Ntot
     write(33, *) solfinal(i)
  end do
  close(33)

end subroutine csfit

end MODULE bregman
