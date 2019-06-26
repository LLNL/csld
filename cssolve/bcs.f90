!!<summary>Wrapper module for calculating finding a solution vector to an
!!underdetermined linear system by minimizing the ell_1 norm with Bayesian
!!constraints. It is assumed that a collection of 'data sets' are available
!!with each one having an 1) experimental measurement and 2) a vector of
!!coefficients from evaluating the set within some basis that matches the
!!constraints of compressive sensing. The collection of data sets' function
!!evaluation vectors is called 'full_pi' and the measurements vector is called
!!'y'.</summary>
module bcs
  use num_types
  use laplace
  public choose_n_random_sets, do_bcs, do_reweighted, do_normal, write_results
  public predict, validate, do_wrapped
contains
  !!<summary>The purpose of this routine is to find a set of structures that are
  !!uncorrelated (in the sense discussed by Candes in his Compressive
  !!Sampling review).</summary>
  !!<comments author="GLWH" date="March 2013">Original routine by LJN March 2013;
  !!Updated to Vidvud's algorithm by GLWH March 2013.</comments>
  !!<comments>The basic algorithm is this:
  !! 1. Start with an empty set of selected data sets (a data set is a vector of
  !!    basis function evaluations for a single CS measurement).
  !! 2. Generate a random vector "pi" on the unit hypersphere.
  !! 3. Orthogonalize "pi" to the subspace spanned by the current set of data sets.
  !! 4. Normalize "pi".
  !! 5. Find the nearest data set to the orthonormalized "pi" using Gus' ArcCos distance.
  !! 6. Add this data set to the set of selected sets.
  !! 7. Construct an orthonormal basis in the space of selected data sets.
  !! 8. Go back to 2 until the desired size of selected set is reached.
  !!</comments>
  !!<parameter name="full_pi" regular="true">The entire matrix of data sets (as defined above).</parameter>
  !!<parameter name="nrandsets" regular="true">The number of random sets to select.</parameter>
  !!<parameter name="selected" regular="true">The indices of the randomly selected data sets (rows in the
  !!@CREF[param.full_pi] matrix.</parameter>
  !!<parameter name="seed_" regular="true">The seed for the random number generator; set using the
  !!system clock if unspecified.</parameter>
  !!<parameter name="outpi_" regular="true">The matrix that includes only the selected rows from the random projection.
  !!</parameter>
  !!<parameter name="y_" regular="true">The vector of experimental measurements if the user would like them to be
  !!selected for @CREF[param.outy_].</parameter>
  !!<parameter name="outy_" regular="true">The vector that includes only the selected measurements from @CREF[param.y_].
  !!</parameter>
  SUBROUTINE choose_n_random_sets(full_pi, nrandsets, selected, seed_, outpi_, y_, outy_)
    real(dp), allocatable :: full_pi(:,:)
    integer, intent(in) :: nrandsets
    integer, intent(inout), pointer :: selected(:)
    integer, intent(inout), optional, allocatable :: seed_(:)
    real(dp), optional :: outpi_(nrandsets, size(full_pi, 2)), y_(size(full_pi, 1)), outy_(nrandsets)

    !!<local name="n, clock">The length of the random seed and system clock time.</local>
    !!<local name="nsets, nbasis">The number of data sets and basis function evaluations
    !!(rows and columns) in @CREF[param.full_pi].</local>
    !!<local name="randset">A randomly data set of basis function evaluations.</local>
    !!<local name="tempset">Variable used for normalizing the data set rows in
    !!@CREF[param.full_pi] for comparison with @CREF[randset].</local>
    !!<local name="closest">*Minimum* distance to the random vector for the current data set so
    !!far when sorting through data sets to find the one closest to the random vector.</local>
    !!<local name="distance">Distance to the random vector for the current data set
    !!when sorting through data sets.</local>
    !!<local name="keepidx">Index of the current data set in the @CREF[param.full_pi]
    !!that is closest.</local>
    !!<local name="inList">Used to hold the result from checking whether a data set is
    !!already in the solution vector of indices @CREF[param.selected].</local>
    integer :: n, clock
    integer :: nsets, nbasis
    real(dp), allocatable :: randset(:), tempset(:)
    real(dp) :: closest, distance ! For sorting the candidate structures
    integer :: i, j, keepIdx, ipi, iy
    integer :: inList(1)
    integer, allocatable :: seed(:)

    if (.not. present(seed_) .or. (.not. allocated(seed_))) then
       call random_seed(size = n)
       allocate(seed(n))
       call system_clock(count = clock)
       seed = clock + 37 * (/ (i -1, i = 1, n)/)
       call random_seed(put = seed)
    else
       call random_seed(put = seed_)
    end if
    
    nsets = size(full_pi, 1)
    nbasis = size(full_pi, 2)
    ipi = 1
    iy = 1

    allocate(randset(nbasis), tempset(nbasis))
    allocate(selected(nrandsets))
    selected = 0

    ! Loop until the desired number of uncorrelated data sets are found.
    do i = 1, nrandsets
       ! Get a random vector of basis function evaluations.
       call random_number(randset)
       ! Rescale the vector so that each number is between -1 and 1
       randset = randset *2 - 1
       ! Before sorting through the list of data sets and determining
       ! which set this vector is "closest" to, orthogonalize it
       ! with respect to all the structures already chosen
       randset = randset/norm(randset)
       
       if (i > 1) then 
          call orthogonalize_to_set_list(randset, full_pi(selected(1:i-1),:))
       end if

       !We are hoping that the minimum distance to the random vector will be less than this
       !it isn't really a risky assumption.
       closest = 1000 
       keepIdx = -1
       do j = 1, nrandsets
          ! Normalize this data set's basis function evaluations' vector.
          tempset(:) = full_pi(j,:)/norm(full_pi(j,:))
          ! Compute the distance (arc length) between this data set's normalized
          ! evaluations' vector to the random vector drawn above. 
          distance = acos(dot_product(tempset,randset))
          ! Is data set j in the list already?
          inList = maxloc(selected, mask=(selected == j))
          ! If not and this distance is smaller than all previous, save this data set.
          if (distance < closest .and. inList(1) == 0) then
             closest = distance
             keepIdx = j
          end if
       end do
       ! Add the data set that was the closest to the random vector in our final list.
       selected(i) = keepIdx
       if (present(outpi_)) then
          outpi_(ipi, :) = full_pi(keepIdx, :)
          ipi = ipi + 1
       end if
       if (present(y_) .and. present(outy_)) then
          outy_(iy) = y_(keepIdx)
          iy = iy + 1
       end if
    end do
  end subroutine choose_n_random_sets

  !!<summary>This is a "modified" Gram-Schmidt orthogonalization routine.</summary>
  !!<comments>Modified Gram-Schmidt routine retrieved January 2015 from
  !!http://terminus.sdsu.edu/SDSU/Math543_s2008/Lectures/07/lecture.pdf; while it
  !!isn't as good as householder, it is really simple.</comments>
  !!<parameter name="A" regular="true">The matrix whose columns represent the vectors
  !!to orthogonalize.</parameter>
  !!<parameter name="Q" regular="true">The orthogonal matrix from the QR decomposition by
  !!classic Gram-Schmidt.</parameter>
  subroutine gram_schmidt(A, Q)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: Q(size(A, 1), size(A, 2))

    !!<local name="qi">The ith unit vector for Q.</local>
    !!<local name="ncol, nrow">The number of columns and rows in A.</local>
    real(dp) :: qi(size(A, 1))
    integer ncol
    integer i, j

    ncol = size(A, 2)
    Q = A
    
    do i = 1, ncol
       qi = Q(:, i)/norm(Q(:, i))
       do j = (i+1), ncol
          Q(:,j) = Q(:,j) - dot_product(qi, Q(:, j))*qi
       end do
    end do
  end subroutine gram_schmidt
  
  !!<summary>This routine takes a random vector on a unit hypersphere and
  !!orthogonalizes it with respect to an existing list of vectors
  !!(themselves not necessarily orthogonal). The output vector is normalized.</summary>
  !!<parameter name="randvec" regular="true">The random vector to orthogonalize.</parameter>
  !!<parameter name="veclist" regular="true">The list of existing vectors to
  !!orthogonalize against.</parameter>
  subroutine orthogonalize_to_set_list(randvec, veclist)
    real(dp), intent(inout) :: randvec(:)
    real(dp), intent(in)    :: veclist(:,:)

    !!<local name="ortholist">The list of orthogonalized data sets in *row* format.</local>
    !!<local name="tr_ortholist">The *transposed* list of orthogonalized data sets
    !!returned by Gram-Schmidt in *column* format.</local>
    real(dp) :: ortholist(size(veclist, 1), size(veclist, 2))
    real(dp) :: tr_ortholist(size(veclist,2), size(veclist, 1))
    integer :: i

    call gram_schmidt(transpose(veclist), tr_ortholist)
    ortholist = transpose(tr_ortholist)

    ! You want to do this projection in a loop (rather than "vectorized"),
    ! updating the vector as you go because it is numerically stable that
    ! way. (See discussions about stability of Gram-Schmidt...same issues
    ! are important here.)
    do i = 1, size(veclist, 1)
       randvec = randvec - dot_product(randvec, ortholist(i,:))*ortholist(i,:)
    enddo
    randvec = randvec/norm(randvec)
  end subroutine orthogonalize_to_set_list

  !!<summary>Returns the value of the weight that should be used for the
  !!weighting matrix in the bcs reweighting scheme.</summary>
  !!<comments author="Conrad W. Rosenbrock" date="11 Jun 2013">
  !!This method is based on Candes "Enhancing sparsity by reweighted ell 1 minimization",
  !!eq 9 and discussion following additional penalty functions are based on the
  !!transfer functions used in machine learning -> neural networks because they
  !!share similarities in the requirements of being bounded, differentiable, etc.    
  !!</comments>
  !!<parameter name="j0">The current value of the j-coefficient for the
  !!corresponding basis</parameter>
  !!<parameter name="epsilon">!A fudge factor for numerical stability.</parameter>
  !!<parameter name="fxn">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  !!<skip enabled="true">This function just has hard-coded analytic functions; there
  !!is nothing tricky since only scalars are involved.</skip>
  real(dp) function reweight_penalty(j0, epsilon, fxn)
    real(dp), intent(in) :: j0 
    real(dp), intent(in) :: epsilon 
    character(len=6), intent(in), optional :: fxn 
 
    select case (fxn)
    case ('logsum')
       reweight_penalty = j0 + epsilon
    case ('logsig')
       reweight_penalty = 4*exp(-epsilon + j0)*(1 + exp(epsilon - j0))**2
    case ('arctan')
       reweight_penalty = j0**2 + epsilon**2
    case ('quarti')
       reweight_penalty = 1+(j0 + 0.7*epsilon)**4
    case ('hexics')
       reweight_penalty = 1+(j0 + 0.57*epsilon)**6
    case ('octics')
       reweight_penalty = 1+(j0 + 0.44*epsilon)**8
    case default
       reweight_penalty = j0 + epsilon
    end select
  end function reweight_penalty
  
  !!<summary>Finds the optimal value for sigma^2 and checks its validity to detect
  !!data that isn't well suited to the BCS routine.</summary>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="nfit" regular="true">The number of data sets allowed to be
  !!used in the fit.</parameter>
  real(dp) function get_sigma2(y, nfit)
    real(dp), allocatable, intent(in) :: y(:)
    integer, intent(in) :: nfit

    !!<local name="lg10">The log10 value of the standard deviation for renormalizing.</local>
    real(dp) :: lg10

    !This estimate for sigma^2 was originally related to the value for \beta suggested by Babacan
    !following the discussion after algorithm 1 in the paper. CWR adjusted it after some
    !experimentation and it gave lower errors for the unit test cases and was more stable for the
    !CEs that were breaking with a negative log error.
    get_sigma2 = (sum( (/( y(i)**2, i=1,nfit )/) )/nfit - (sum( (/( y(i), i=1,nfit )/) )/nfit)**2)*0.01
    ! lg10 = log10(get_sigma2)
    ! get_sigma2 = get_sigma2/(real(10)**ceiling(lg10+1))
  end function get_sigma2
  
  !!<summary>This routine wraps the Bayesian CS routine with an outside loop. Only
  !!regular ell_1 minimization is done (as opposed to p &lt; 1 reweighted minimization).</summary>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi" regular="true">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="sigma2" regular="true">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta" regular="true">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="js" regular="true">The solution vector for the reweighted BCS.</parameter>
  !!<parameter name="error_bars" regular="true">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix.</parameter>
  subroutine do_normal(full_pi, y, sigma2, eta, js, error_bars)
    real(dp), allocatable, intent(in) :: full_pi(:,:), y(:)
    real(dp), intent(in) :: sigma2, eta
    real(dp), intent(out) :: js(size(full_pi, 2)), error_bars(size(full_pi, 2))

    !!<local name="iterator">The laplace iterator to find the solution vector.</local>
    !!<local name="returnsigma2">The value of sigma2 calculated by the laplace
    !!iterator if one was not provided. We don't actually output this anywhere.</local>
    real(dp) returnsigma2
    type(laplace_iterator), pointer :: iterator
    
    allocate(iterator)
    call iterator%initialize(full_pi, y, sigma2, eta)
    call iterator%iterate(js, error_bars, returnsigma2)
    deallocate(iterator)
  end subroutine do_normal

  !!<summary>This routine wraps the Bayesian CS routine with an outside loop. Only
  !!regular ell_1 minimization is done (as opposed to p &lt; 1 reweighted minimization).</summary>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi" regular="true">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="sigma2" regular="true">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta" regular="true">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="js">The solution vector for the reweighted BCS.</parameter>
  !!<parameter name="error_bars">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix.</parameter>
  subroutine do_wrapped(full_pi, y, sigma2, eta, js, error_bars)
    real(dp), intent(in) :: full_pi(:,:), y(:)
    real(dp), intent(in) :: sigma2, eta
    real(dp), intent(out) :: js(size(full_pi, 2)), error_bars(size(full_pi, 2))

    !!<local name="iterator">The laplace iterator to find the solution vector.</local>
    !!<local name="returnsigma2">The value of sigma2 calculated by the laplace
    !!iterator if one was not provided. We don't actually output this anywhere.</local>
    real(dp) returnsigma2
    type(laplace_iterator), pointer :: iterator
    
    allocate(iterator)
    call iterator%initialize(full_pi, y, sigma2, eta)
    call iterator%iterate(js, error_bars, returnsigma2)
    deallocate(iterator)
  end subroutine do_wrapped

  !!<summary>Added by Lance Nelson: 1 Jul 2015 Given a set of Js, this
  !!routine will find which set has the best rms error over the
  !!validation set. </summary> 
  !!<parameter name="trackedJs">list of J vectors for the different solutions.</parameter>-
    !!<parameter name="hold_pi">The correlations of the holdout set to validate against for dealing
  !!with the ringing problem. If unspecified, ringing problems may persist.</parameter>
  !!<parameter name="hold_y">As for @CREF[hold_pi], the corresponding measurements for choosing a
  !!best solution when ringing occurs.</parameter>
  integer function bestSolution(hold_pi, hold_y, trackedJs)
    real(dp), intent(in) :: hold_pi(:,:), hold_y(:)
    real(dp):: trackedJs(:,:)
    !!<local name="rmsErrors">The rms error for each fit in @CREF[trackedJs].</local>
    !!<local name="pred_energies, pred_diffs">The predicted energies for the current fit
    !!and its absolute difference from the actual answers.</local>
    !!<local name="nhold">The number of data sets in the holdout set to validate against.</local>
    real(dp), pointer :: rmsErrors(:)
    real(dp), allocatable :: pred_energies(:)
    integer :: nhold
    integer :: i, l

    nhold = size(hold_pi, 1)
    allocate(pred_energies(nhold))
    allocate(rmsErrors(size(trackedJs, 1)))
    
    do i = 1, size(trackedJs,1)
       pred_energies(:) = (/ (dot_product(hold_pi(l,:), trackedJs(i,:)), l=1,nhold) /)
       rmsErrors(i) = sqrt(sum((pred_energies - hold_y)**2)/nhold)
    enddo
    bestSolution = minloc(rmsErrors, 1)
  end function bestSolution


  !!<summary>Added by Lance Nelson: 1 Jul 2015 Determines whether the
  !!reweighting algorithm is failing to converge.  Sometimes the
  !!algorithm exhibts a ringing behavior where it oscillates between
  !!several different solutions rather than converging.  If this is
  !!the case we want to take action and quickly end the iterative
  !!procedure while saving the best solution.  Best in this case will
  !!be determined by the rms error over the validation set. </summary>
  !!<parameter name="trackedell0s">list of \ell_0 norms for the
  !!previous five(5) iterations.  </parameter>
  logical function isRinging(idx, trackedell0s)
    integer, intent(in) :: trackedell0s(:)
    integer, intent(in) :: idx
    
    !!<local name="subsetell0s">For periodicity 'i', the last 'i' fits tracked.</local>
    !!<local name="ell0perm">Used to compare @CREF[subsetell0s] to the rest of the fits in
    !!the tracking list.</local>
    !!<local name="nTrack">The number of reweight iterations tracked so far.</local>
    !!<local name="maxLength">The (even) number of iterations that will actually be used
    !!to assess the periodicity of the system.</local>
    integer, pointer :: subsetell0s(:)
    integer, pointer :: ell0perm(:)
    integer :: nTrack, maxLength, i

    isRinging = .false.  !Default value for ringing
    if (idx < size(trackedell0s,1)) then  ! If we haven't tracked enough solutions make an accurate assessment, just exit
       return
    endif

    nTrack = size(trackedell0s,1)
    ! To assess any periodicity in the solutions, we need an even length list of solutions.
    ! In the case that the size of trackJs is odd, let's just consider the last n-1 (which is even)
    ! elements
    if (mod(2,nTrack) == 1) then
       maxLength = nTrack - 1
    else
       maxLength = nTrack
    endif
    ! This do loop is to check any order periodicity, starting at 2
    ! and ending at n (or n-1 if n is odd)/2 where n is the number of tracked solutions.
    ! For example, if you are tracking 11 solutions, this 'do' loop will check
    ! periodicity of lengths: 2,3,4,and 5.  If larger periodicites ever needed
    ! to be detected, we would need to increase the number of solutions being
    ! tracked.  This is controlled using the variable "nTrack" found in this function's
    ! calling routine.
    do i = 4, maxLength,2
       allocate(subsetell0s(i), ell0perm(i) )
       subsetell0s = trackedell0s(nTrack-i+1:)  !Go get the last i solutions
       ell0perm = 0
       ell0perm(i/2+1:) = subsetell0s(1:i/2)  ! Shift the list by i/2 (half its length)
       ell0perm(:i/2) = subsetell0s(i/2+1:)
       if (all(ell0perm(:) == subsetell0s(:)) ) then ! See if the two lists are equal
          isRinging = .true.
          deallocate(subsetell0s, ell0perm)
          return
       endif
       deallocate(subsetell0s, ell0perm)          
    enddo
  end function isRinging
  
  !!<summary>This routine wraps the Bayesian CS routine with an outside loop that
  !! applies the p &lt; 1 norm reweighting scheme.</summary>
  !!<comments>See: E. J. Candes, M. Wakin and S. Boyd. Enhancing sparsity by
  !! reweighted \ell_1 minimization. J. Fourier Anal. Appl., 14 877-905.</comments>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi" regular="true">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="sigma2" regular="true">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta" regular="true">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="jcutoff" regular="true">The minimum value a J coefficient has to have
  !!in order to be kept between reweighting iterations.</parameter>
  !!<parameter name="penaltyfxn" regular="true">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  !!<parameter name="js" regular="true">The solution vector for the reweighted BCS.</parameter>
  !!<parameter name="error_bars" regular="true">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix.</parameter>
  !!<parameter name="hold_pi">The correlations of the holdout set to validate against for dealing
  !!with the ringing problem. If unspecified, ringing problems may persist.</parameter>
  !!<parameter name="hold_y">As for @CREF[hold_pi], the corresponding measurements for choosing a
  !!best solution when ringing occurs.</parameter>
  subroutine do_reweighted(full_pi, y, sigma2, eta, jcutoff, penaltyfxn, js, error_bars, hold_pi, hold_y)
    real(dp), allocatable, intent(in) :: full_pi(:,:), y(:)
    real(dp), intent(in) :: eta, jcutoff
    real(dp), intent(inout) :: sigma2
    character(len=6), intent(in) :: penaltyfxn
    real(dp), intent(out) :: js(size(full_pi, 2)), error_bars(size(full_pi, 2))
    real(dp), optional, intent(in) :: hold_pi(:,:), hold_y(:)

    !!<local name="returnsigma2">The value of sigma2 calculated by the laplace
    !!iterator if one was not provided. We don't actually output this anywhere.</local>
    !!<local name="w_pi">The *weighted* @CREF[param.full_pi] matrix.</local>
    !!<local name="weight_matrix">The matrix of adjusted/re-weighted J values.</local>
    !!<local name="nsets, nbasis">The number of data sets (rows) and basis functions
    !!(columns) in the @CREF[param.full_pi] matrix.</local>
    !!<local name="copyJs">A copy of the current js solution vector for manipulation.</local>
    !!<local name="ell0, prevell0">The iteration continues until the ell_0 norm of
    !!the solution vector stops changing significantly. These variables keep track of
    !!the previous and current ell_0 norm.</local>
    !!<local name="iterator">The laplace iterator to find the solution vector.</local>
    !!<local name="i_0">The maximum number of non-negligible coefficient values
    !!in order for the compressive sensing mathematics to be applicable.</local>
    !!<local name="maxidx">The index of the maximum value in the copy of the J
    !!solution vector. Helps determine the value of the reweighting epsilon parameter.</local>
    !!<local name="jmax">The value of the largest j value outside of the range
    !!specifed by 'i_0'.</local>
    !!<local name="epsilon">The reweighting parameter that controls how aggresively the
    !!weight matrix is altered.</local>
    real(dp) :: returnsigma2
    real(dp), allocatable :: w_pi(:,:)
    real(dp), allocatable :: weight_matrix(:,:)
    integer :: nsets, nbasis
    real(dp) :: copyJs(size(full_pi,2))
    integer :: ell0, prevell0
    type(laplace_iterator), pointer :: iterator
    integer  :: i_0, maxidx(1)
    real(dp) :: jmax(1), epsilon
    integer :: i
    !!<local name="trackJs, trackell0">For coping with solution "ringing", tracks the fits
    !!and their ell_0 norms.</local>
    !!<local name="ntrack">The number of previous reweighting fits to check for periodicity.</local>
    real(dp), allocatable :: trackJs(:,:)
    integer,allocatable :: trackell0(:)
    integer :: idx, best, ntrack

    !We start off with the identity matrix for the weight matrix and then
    !update it after each iteration.
    nbasis = size(full_pi, 2)
    nsets = size(full_pi, 1)
    
    allocate(weight_matrix(nbasis,nbasis), w_pi(nsets,nbasis))
    allocate(iterator)
    weight_matrix = 0
    do i = 1, nbasis
       weight_matrix(i,i) = 1
    end do

    !If the value for sigma2 is not right, choose a better one.
    if (sigma2 .lt. 0) sigma2 = get_sigma2(y, nsets)
    
    ntrack = 13
    allocate(trackJs(ntrack, nbasis), trackell0(ntrack))
    trackJs = 0
    trackell0 = 0
    idx = 0
    
    !Here we make a rough estimate of the number of J coefficients that will
    !have non-negligible values (i.e. > 1e-3). We hope since the solution is
    !supposed to be sparse that this condition will hold.
    i_0 = nsets/(4*log(real(nbasis)/real(nsets)))
    prevell0 = -1
    
    do while (.true.)
       w_pi = matmul(full_pi, weight_matrix)
       call iterator%initialize(w_pi, y, sigma2, eta)
       call iterator%iterate(js, error_bars, returnsigma2)

       where(abs(js) .le. jcutoff) js = 0
       ell0 = count(abs(js) > jcutoff)

       trackJs(:nTrack - 1,:) = trackJs(2:nTrack,:)  ! Slide all previous Js down one index
       trackJs(nTrack,:) = js(:) ! Save the current Js at the last index
       trackell0(:nTrack-1) = trackell0(2:nTrack) ! Slide all the previous \ell_0s down one index
       trackell0(nTrack) = ell0 ! Save the current \ell_0 norm.
       if (isRinging(idx, trackell0)) then  ! Check the \ell_0 norms to see if the algorithm is "ringing"
          if (present(hold_pi) .and. present(hold_y)) then
             !Since we have a validation set, we can choose the best predictor of the ringing
             !solutions.
             best = bestSolution(hold_pi, hold_y, trackJs)
          else
             !Choose the one with the smallest ell0 norm since we don't have any validation set to work with.
             best = minloc(trackell0, 1)
          end if
          js(:) = trackJs(best,:)  ! Save the "best" Js.
          exit
       endif
       
       if (abs(prevell0-ell0) .le. 5) then
          exit
       else
          prevell0 = ell0
       end if

       !We should only reweight the js if we are going to do another reweighting
       !run on the data.
       js(:) = matmul(weight_matrix, js)

       !Set the largest i_0 j-values to zero so we can get an estimate of
       !the size of the smallest coefficients. They affect the reweighting
       !matrix and penalty functions.
       copyJs = abs(js)
       do i=1, i_0-1
          maxidx = maxloc(copyJs)
          copyJs(maxidx(1)) = 0
       end do
       jmax = maxval(copyJs)
       if (jmax(1) > 10e-3) then
          epsilon = jmax(1)
       else
          epsilon = 10e-3
       end if

       !Adjust the values of the diagonals on the weighting matrix using the
       !specified penalty function.
       do i=1, nbasis
          !This is the inverse W matrix.
          weight_matrix(i,i) = reweight_penalty(js(i), epsilon, penaltyfxn)
       end do
       
       idx = idx + 1
       call iterator%reset()
    end do

    !Perform cleanup on the iterator object. This will call its finalizer once
    !gfortran implements them. For now we rely on the built-in deallocation.
    deallocate(iterator)    
  end subroutine do_reweighted

  !!<summary>Partitions the collection of data sets into a set for fitting and a set
  !!for validating the fit.</summary>
  !!<parameter name="nfits" regular="true">The number of data sets to use in the fit.</parameter>
  !!<parameter name="nsets" regular="true">The total number of sets available to select from.</parameter>
  !!<parameter name="nholdout" regular="true">The number of data sets to
  !!retain for validation.</parameter>
  !!<parameter name="fitlist" regular="true">A list of the row indices in 'full_pi' that should be used
  !!for fitting.</parameter>
  !!<parameter name="holdlist" regular="true">A list of the row indices in 'full_pi' that should be held
  !!out for validation after the fitting.</parameter>
  !!<parameter name="seed_" regular="true">The seed for the random number generator; set using the
  !!system clock if unspecified.</parameter>
  subroutine partition_holdout_set(nfits, nsets, nholdout, fitlist, holdlist, seed_)
    integer, intent(in) :: nfits, nsets, nholdout
    integer, intent(inout) :: fitlist(nfits, nsets-nholdout), holdlist(nfits, nholdout)
    integer, optional, allocatable :: seed_(:)

    !!<local name="n, clock">The length of the random seed and system clock time.</local>
    !!<local name="list">Holds a temporory list of the remaining items once one has
    !!been randomly selected for inclusion in the fit list.</local>
    !!<local name="indx">The index of the randomly selected set.</local>
    !!<local name="item">Temporary for holding the value of the randomly selected index.</local>
    !!<local name="ifit, iset, i">Iterators over the fits, sets and loop constructs.</local>
    !!<local name="r">Holds the randomly selected value for choosing random fit sets.</local>
    integer :: n, clock
    integer :: ifit, i, list(nsets)
    real(dp) :: r
    integer :: indx, iset, item
    integer, allocatable :: seed(:)

    if (.not. present(seed_) .or. (.not. allocated(seed_))) then
       call random_seed(size = n)
       allocate(seed(n))
       call system_clock(count = clock)
       seed = clock + 37 * (/ (i -1, i = 1, n)/)
       call random_seed(put = seed)
    else
       call random_seed(put = seed_)
    end if
    
    do ifit = 1, nfits
       list = (/(i, i=1, nsets)/)
       do iset = 1, nsets
          call random_number(r)
          ! Pick the index of item at random from the remaining list
          indx = ceiling(r*(nsets - iset + 1))
          item = list(indx) ! Get the value of the list at the selected index
          list = (/list(1:indx-1), list(indx+1:size(list)), item/) ! Update the list of items
       enddo

       fitlist(ifit,:) = list(1:nsets-nholdout)
       holdlist(ifit,:) = list(nsets-nholdout+1:)
    end do
  end subroutine partition_holdout_set

  !!<summary>Returns lowest i/o unit number not in use.</summary>
  !!<parameter name="unit">Out parameter that will contain the lowest i/o number.</parameter>
  !!<skip enabled="true">Only involves scalars and system functions/procedures.</skip>
  integer function newunit(unit) result(n)
    integer, intent(out), optional :: unit
    logical inuse
    integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
    integer, parameter :: nmax=999  ! may be system-dependent
    do n = nmin, nmax
       inquire(unit=n, opened=inuse)
       if (.not. inuse) then
          if (present(unit)) unit=n
          return
       end if
    end do
    stop "newunit ERROR: available unit not found."
  end function newunit

  !!<summary>Returns a value indicating whether the specified file exists.</summary>
  !!<parameter name="filename">The name of the file to check for (in the current directory).</parameter>
  !!<skip enabled="true">Only involves scalars and system functions/procedures.</skip>
  logical function file_exists(filename)
    character(len=*), intent(in) :: filename
    INQUIRE( FILE=filename, EXIST=file_exists) 
  end function file_exists

  !!<summary>Returns the number of values in the specified line assuming
  !!that they are separated by spaces or tabs.</summary>
  !!<parameter name="length">The number of characters in line.</parameter>
  !!<parameter name="line">The string for the line to count values in.</parameter>
  !!<skip enabled="true">Disabled because it is a copy of the fortpy routine called
  !!'fpy_value_count' that has been tested as part of fortpy.</skip>
  integer function value_count(line, length)
    integer, intent(in) :: length
    character(length), intent(in) :: line
    character(2) :: whitespace
    integer           :: success, i, indx, prev = 1, beginning = 1
    real              :: value

    !Initialize the whitespace array. We will cycle through all the characters
    !in the specified line looking for whitespace. Each time we find it, if the
    !character immediately preceding it was not whitespace, we have a value.
    whitespace = '  '
    value_count = 0

    do i = 1, length
       !indx will be zero if the current character is not a whitespace character.
       indx = index(whitespace, line(i:i))
       !The ichar == 9 statement checks for tabs; we used to have it concatenated onto the
       !whitespace array, but there was a bug, so we switched to explicit behavior.
       if ((indx > 0 .or. ichar(line(i:i)) .eq. 9) .and. prev == 0) then
          !We found the first whitespace after the end of a value we want.
          value_count = value_count + 1
       end if

       prev = indx
    end do

    !If the last value on the line ends right before \n, then we wouldn't have
    !picked it up; add an extra one.
    if (indx == 0) value_count = value_count + 1
  end function value_count
  
  !!<summary>Returns the number of lines in the file that aren't comments and
  !!the number of whitespace-separated values on the first non-comment line.</summary>
  !!<parameter name="filename">The name of the file to pass to open.</parameter>
  !!<parameter name="n">The number of characters in 'filename'.</parameter>
  !!<parameter name="commentchar">A single character which, when present at the start
  !!of a line designates it as a comment.</parameter>
  !!<skip enabled="true">Disabled because it is a copy of the fortpy routine called
  !!'fpy_linevalue_count' that has been tested as part of fortpy.</skip>
  subroutine linevalue_count(filename, n, commentchar, nlines, nvalues)
    integer, intent(in) :: n
    character(n), intent(in) :: filename
    character(1), intent(in) :: commentchar
    integer, intent(out) :: nlines, nvalues
    character(len=:), allocatable :: cleaned
    integer :: ioerr, funit
    character(150000) :: line

    !Initialize the value for the result; if we get an error during the read, we
    !end the loop. It can be caused by badly formatted data or the EOF marker.
    nlines = 0
    nvalues = 0

    open(newunit(funit), file=filename, iostat=ioerr)
    if (ioerr == 0) then
       do
          read(funit, "(A)", iostat=ioerr) line
          if (ioerr == 0) then
             cleaned = trim(adjustl(line))
             if (len(cleaned) .gt. 0) then
                if (cleaned(1:1) /= commentchar) then
                   nlines = nlines + 1
                   !We only need to get the number of values present in a line once.
                   !We restrict the file structure to have rectangular arrays.
                   if (nvalues == 0) then
                      nvalues = value_count(cleaned, len(cleaned))
                   end if
                end if
             end if
          else
             exit
          end if
       end do
    end if
    close(funit)
  end subroutine linevalue_count

  !!<summary>Writes the resulting solution vector, error bars and error prediction
  !!to the specified file unit.</summary>
  !!<parameter name="js">The solution vector for the reweighted BCS. There is one row
  !!in this matrix for each fit performed. The values below the cutoff should already
  !!be set to zero.</parameter>
  !!<parameter name="hold_rms_, fit_rms">The RMS error for the holdout set and fitting
  !!set respectively using the coefficients from the BCS ell_1 minimization.</parameter>
  !!<parameter name="hold_err_, fit_err">The absolute error for the holdout set and fitting
  !!set respectively using the coefficients from the BCS ell_1 minimization.</parameter>
  subroutine write_results(js, fit_rms, fit_err, hold_rms_, hold_err_, sigma2_)
    real(dp), intent(in) :: js(:,:), fit_rms(:), fit_err(:)
    real(dp), intent(in), optional :: hold_rms_(:), hold_err_(:), sigma2_(:)

    !!<local name="l0, l1">The ell_0 and ell_1 norms of the solution vectors.</local>
    !!<local name="funit">File unit for writing the results.</local>
    integer :: l0, funit, i
    real(dp) :: l1
    !The s_* character arrays are for formatting the values of the optional parameters
    !so that we can insert blanks if they were not specified.
    character(22) :: s_sigma2, s_fmt
    character(20) :: s_fitrms, s_fiterr, s_holdrms, s_holderr

    open(unit=newunit(funit), file='summary.out', status="replace")
    s_fmt = '(A5, A20, A5, 4A20, A22)'
    write(funit, s_fmt) 'Fit #','|J|_1','|J|_0','RMS(fit)','ERR(fit)','RMS(hold)','ERR(hold)','sigma^2'
    !Overwrite the string format specified for the header with the data item format.
    s_fmt = '(I5, F20.5, I5, 2F20.5, 2A20, A22)'
    
    do i=1, size(js, 1)
       !Handle the writing of the optional arguments; if they weren't specified, they shouldn't show up.
       if (present(hold_err_)) then
          write(s_holderr, '(F20.5)') hold_err_(i)
       else
          s_holderr = ''
       end if
       if (present(hold_rms_)) then
          write(s_holdrms, '(F20.5)') hold_rms_(i)
       else
          s_holdrms = ''
       end if
       if (present(sigma2_)) then
          write(s_sigma2, '(F22.10)') sigma2_(i)
       else
          s_sigma2 = ''
       end if   
       
       l1 = sum(abs(js(i, :)))
       l0 = count(abs(js(i, :)) > 0.0)
       write(funit, s_fmt) i, l1, l0, fit_rms(i), fit_err(i), s_holdrms, s_holderr, s_sigma2
    end do
    close(funit)
  end subroutine write_results

  !!<summary>Predicts measurement values for the specified matrix of function evalutions
  !!and a BCS solution vector.</summary>
  !!<parameter name="pi" regular="true">A matrix of data sets (rows) with each set evaluated at the basis
  !!functions in @CREF[param.js].</parameter>
  !!<parameter name="js" regular="true">The BCS solution vector with coefficients for each basis function
  !!forming the model.</parameter>
  !!<parameter name="prediction" regular="true">The resulting vector of the model's
  !!measurement predictions.</parameter>
  subroutine predict(pi, js, prediction)
    real(dp), allocatable, intent(in) :: pi(:,:)
    real(dp), intent(in) :: js(size(pi, 2))
    real(dp), allocatable, intent(inout) :: prediction(:)

    integer :: k
    if (allocated(prediction)) deallocate(prediction)
    allocate(prediction(size(pi, 1)))
    prediction(:) = (/ (dot_product(pi(k,:), js(:)), k=1, size(pi, 1)) /)
  end subroutine predict

  !!<summary>Predicts measurement values for the specified matrix of function evalutions
  !!and a BCS solution vector, and validates the prediction against a known solution.</summary>
  !!<parameter name="pi" regular="true">A matrix of data sets (rows) with each set evaluated at the basis
  !!functions in @CREF[param.js].</parameter>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="js" regular="true">The BCS solution vector with coefficients for each basis function
  !!forming the model.</parameter>
  !!<parameter name="prediction" regular="true">The resulting vector of the model's
  !!measurement predictions.</parameter>
  !!<parameter name="pred_err" regular="true">The absolute error over the entire
  !!prediction vector.</parameter>
  !!<parameter name="pred_rms" regular="true">The RMS error for the entire prediction vector.</parameter>
  subroutine validate(pi, y, js, prediction, pred_err, pred_rms)
    real(dp), allocatable, intent(in) :: pi(:,:)
    real(dp), allocatable, intent(in) :: y(:)
    real(dp), intent(in) :: js(size(pi, 2))
    real(dp), allocatable, intent(inout) :: prediction(:)
    real(dp), intent(out) :: pred_err, pred_rms

    call predict(pi, js, prediction)
    pred_err = sum(abs(prediction - y))
    pred_rms = sqrt(sum((prediction - y)**2)/size(pi, 1))
  end subroutine validate
  
  !!<summary>Performs BCS for the given sensing matrix and measurement vector.</summary>
  !!<parameter name="y" regular="true">The experimental measurements for each data set.</parameter>
  !!<parameter name="full_pi" regular="true">The full matrix of basis function evaluations for
  !!each data set for which we have measurements.</parameter>
  !!<parameter name="nfits" regular="true">The number of random fits to perform with these args.
  !!Specify a value of 1, for minimal requirements to run.</parameter>
  !!<parameter name="js" regular="true">The solution vector for the reweighted BCS. There is one row
  !!in this matrix for each fit performed.</parameter>
  !!<parameter name="error_bars" regular="true">Error bars for the solution vector @CREF[param.js]; taking
  !!from the diagonal entries of the covariance matrix. There is one row in this matrix
  !!for each fit performed.</parameter>
  !!<parameter name="hold_rms_, fit_rms" regular="true">The RMS error for the holdout set and fitting
  !!set respectively using the coefficients from the BCS ell_1 minimization.</parameter>
  !!<parameter name="hold_err_, fit_err" regular="true">The absolute error for the holdout set and fitting
  !!set respectively using the coefficients from the BCS ell_1 minimization.</parameter>
  !!<parameter name="nholdout_" regular="true">If validation is desired, how many of the data sets in
  !!@CREF[param.full_pi] should be used for validation. If unspecified, then no validation
  !!is performed unless the file 'holdout' is present in the current directory, in which
  !!case the values in that file are used instead.</parameter>
  !!<parameter name="sigma2_">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta_" regular="true">Threshold for stopping the iteration algorithm.</parameter>
  !!<parameter name="jcutoff_" regular="true">The minimum value a J coefficient has to have
  !!in order to be kept between reweighting iterations.</parameter>
  !!<parameter name="penaltyfxn_" regular="true">The name of the penalty function to use. Possible values:
  !!logsum, logsig, arctan, quarti, hexics, octics.</parameter>
  !!<parameter name="reweight_" regular="true">When specified, the reweighted BCS
  !!routine is used, otherwise a normal ell_1 minimization is used.</parameter>
  !!<parameter name="seed" regular="true">The seed for the random number generator; set using the
  !!system clock if unspecified.</parameter>
  subroutine do_bcs(full_pi, y, nfits, js, error_bars, fit_rms, fit_err, hold_rms_, hold_err_, &
                    nholdout_, sigma2_, eta_, jcutoff_, penaltyfxn_, reweight_, seed)
    real(dp), intent(in) :: full_pi(:,:), y(:)
    integer :: nfits
    real(dp), intent(inout) :: js(nfits, size(full_pi, 2)), error_bars(nfits, size(full_pi, 2))
    real(dp), intent(inout) :: fit_rms(nfits), fit_err(nfits)
    integer, intent(inout), optional :: nholdout_
    real(dp), intent(inout), optional :: hold_rms_(nfits), hold_err_(nfits)
    real(dp), intent(inout), optional :: sigma2_(nfits), eta_, jcutoff_
    character(len=6), intent(in), optional :: penaltyfxn_
    logical, intent(in), optional :: reweight_
    integer, optional, allocatable :: seed(:)
    
    !Handle all the optional parameters; these have the same descriptions as their parameter
    !summary tags decorating the subroutine.
    real(dp) :: sigma2, eta, jcutoff
    character(len=6) :: penaltyfxn
    logical :: reweight

    !Here it may seem counter-intuitive to have a temporary vector for 'js' and 'error_bars'
    !for performing multiple fits; we could just write the result directly into the matrix
    !of solution vectors. However, we want the do_reweighted() and do_normal() subroutines
    !to have nice signatures for single fitting use.
    !!<local name="tjs">A single column of j solution vectors for re-using in the loop.</local>
    !!<local name="terror_bars">A single column of error bars for the 'js' vector.</local>
    !!<local name="nsets, nbasis">The number of data sets (rows) and function evalutions (cols)
    !!in the @CREF[param.full_pi] matrix.</local>
    !!<local name="nfitsets">The number of data sets allowed to be used for fitting.
    !!Will be all of them unless holdout sets for validation are specified.</local>
    !!<local name="nlines, nvalues">The number of lines and values in a 2D data file being
    !!read in (for example 'holdout').</local>
    !!<local name="nfixed">The number of data sets specified in the 'fixed' file.</local>
    !!<local name="funit">The next available file unit for reading data files.</local>
    !!<local name="sub_pi, sub_y">A subset of the @CREF[param.full_pi] matrix rows (sub_pi) and
    !!columns (sub_y) respectively for a single 'fitlist'/'holdlist' pair.</local>
    !!<local name="fit_pred, hold_pred">The predicted measurements for the fitting and
    !!holdout sets respectively using the solution vector calculated with BCS.</local>
    real(dp) :: tjs(size(full_pi, 2)), terror_bars(size(full_pi, 2))
    integer, allocatable :: fitlist(:,:), holdlist(:,:), fixed(:)
    integer :: i, j, k, nsets, nbasis, nfitsets, nlines, nvalues, funit
    integer :: nfixed
    real(dp), allocatable :: sub_pi(:,:), sub_y(:), fit_pred(:), hold_pred(:)

    !If the fit needs to be validated afterward, we need to split the data sets in 'full_pi'
    !into sets for fitting and those held out for validation afterward.
    nsets = size(full_pi, 1)
    nbasis = size(full_pi, 2)
    nfitsets = 0

    if (present(nholdout_) .and. nholdout_ .ne. 0) then
       allocate(fitlist(nfits, nsets-nholdout_), holdlist(nfits, nholdout_))
       if (present(seed)) then
          call partition_holdout_set(nfits, nsets, nholdout_, fitlist, holdlist, seed)
       else
          call partition_holdout_set(nfits, nsets, nholdout_, fitlist, holdlist)
       end if
       nfitsets = nsets - nholdout_
    end if

    !The user can also specify the holdout set using a file called 'holdout' in the current
    !directory. If it exists, set up the fitlist and hold list arrays.
    if ((.not. allocated(fitlist)) .and. file_exists('holdout')) then
       call linevalue_count('holdout', 7, '#', nlines, nvalues)
       if (nvalues .gt. nsets .or. nlines .ne. nfits) then
          write(*, *) "Invalid dimensions for 'holdout' file: ", nlines, " x ", nvalues
          write(*, *) "We need as many rows in 'holdout' as requested fits: ", nfits
          write(*, *) "The number of holdout values must be less than total data sets: ", nsets
          stop
       end if
       
       open(newunit(funit), file='holdout', action='read')
       allocate(holdlist(nlines, nvalues))
       read(funit, *)(holdlist(i,:), i=1, nlines)
       close(funit)

       !Now we construct the fit list by taking the complement of those indices that
       !are present in the hold list.
       allocate(fitlist(nlines, nsets-nvalues))
       k = 1
       do i=1, nlines
          do j=1, nsets
             if (.not. any(j .eq. holdlist(i,:))) then
                fitlist(i,k) = j
                k = k+1
             end if
          end do
       end do
       nfitsets = nsets - nvalues
       if (present(nholdout_)) nholdout_ = nvalues
    end if

    !If no holdout set is specified for validation, we use the entire fitting set.
    if (nfitsets .eq. 0) nfitsets = nsets
    
    !One more option that is useful sometimes is to define a collection of data sets that *must*
    !be kept in the fitting set always.
    if (file_exists('fixed')) then
       call linevalue_count('fixed', 5, '#', nlines, nfixed)
       if (nfixed .gt. nsets) then
          write(*, *) "The number of fixed data sets specified exceeds the possible value: ", nsets
          stop
       end if

       allocate(fixed(nfixed))
       open(newunit(funit), file='fixed', action='read')
       read(funit, *) fixed
       close(funit)
    else
       nfixed = 0
    end if

    !Next we set the default values for the optional parameters if they are not present.
    if (present(eta_)) then
       eta = eta_
    else
       eta = 1e-8
    end if
    if (present(jcutoff_)) then
       jcutoff = jcutoff_
    else
       jcutoff = 1e-3
    end if
    if (present(penaltyfxn_)) then
       penaltyfxn = penaltyfxn_
    else
       penaltyfxn = 'arctan'
    end if
    if (present(reweight_)) then
       reweight = reweight_
    else
       reweight = .true.
    end if

    !The size of the arrays for making predictions don't change between fitting
    !iterations, so we can re-use them.
    allocate(fit_pred(nfitsets))
    if (allocated(holdlist)) allocate(hold_pred(size(holdlist, 2)))
    
    do i=1, nfits
       !Because of the possible random selection of fitting sets and validation sets,
       !we need to allocate separately for each fit.
       allocate(sub_pi(nfitsets, nbasis), sub_y(nfitsets))
       do j=1, nfixed
          sub_pi(j,:) = full_pi(fixed(j),:)
          sub_y(j) = y(fixed(j))
       end do
       do j=(nfixed+1), nfitsets
          if (allocated(fitlist)) then
             sub_pi(j,:) = full_pi(fitlist(i, j-nfixed),:)
             sub_y(j) = y(fitlist(i, j-nfixed))
          else
             sub_pi(j,:) = full_pi(j-nfixed,:)
             sub_y(j) = y(j-nfixed)
          end if
       end do

       !Calculate a special sigma2 value for this fitting set and measurement vector.
       if (present(sigma2_)) then
          if (i .lt. size(sigma2_) .and. sigma2_(i) .ne. 0) then
             sigma2 = sigma2_(i)
          else
             sigma2 = get_sigma2(sub_y, nfitsets)
          end if
       else
          sigma2 = get_sigma2(sub_y, nfitsets)
          if (present(sigma2_)) sigma2_(i) = sigma2
       end if
    print *, 'debug i_fit=', i
    print *, 'debug pi',shape(full_pi), shape(y), 'nfits',nfits,'js',shape(js),shape(error_bars),shape(fit_rms), shape(fit_err)
    print *, 'debug hold', shape(hold_rms_), shape(hold_err_),nholdout_,'sigma2',sigma2, 'eta', eta, 'jcutoff',jcutoff_,penaltyfxn_, reweight_

       tjs = 0
       terror_bars = 0
       if (reweight) then
          call do_reweighted(sub_pi, sub_y, sigma2, eta, jcutoff, penaltyfxn, tjs, terror_bars)
       else
          call do_normal(sub_pi, sub_y, sigma2, eta, tjs, terror_bars)
       end if
       
       !Copy the solution vector and error bars across; if this seems wasteful, see
       !the comment above the variable declaration for tjs and terror_bars.
       js(i,:) = tjs(:)
       error_bars(i,:) = terror_bars(:)

       !We have the resulting solution vector and its error bars; at the least we calculate
       !the RMS error for the fitting set; but if a holdout set was also specified, we need
       !to calculate the error for that as well.
       fit_pred(:) = (/ (dot_product(sub_pi(k,:), tjs(:)), k=1, nfitsets) /)
       fit_err(i) = sum(abs(fit_pred - sub_y))
       fit_rms(i) = sqrt(sum((fit_pred - sub_y)**2)/nfitsets)
       deallocate(sub_pi, sub_y)
       
       !Check for the presence of a holdout set; if we have one calculate the RMS and
       !absolute error for that too.
       if (allocated(holdlist) .and. (present(hold_err_) .or. present(hold_rms_))) then
          allocate(sub_pi(size(holdlist, 2), nbasis), sub_y(size(holdlist, 2)))
          do j=1, size(holdlist, 2)
             sub_pi(j, :) = full_pi(holdlist(i, j), :)
             sub_y(j) = y(holdlist(i, j))
          end do
          
          hold_pred(:) = (/ (dot_product(sub_pi(k,:), tjs(:)), k=1, size(holdlist, 2)) /)
          if (present(hold_err_)) hold_err_(i) = sum(abs(hold_pred - sub_y))
          if (present(hold_rms_)) hold_rms_(i) = sqrt(sum((hold_pred - sub_y)**2)/size(holdlist, 2))
          deallocate(sub_pi, sub_y)
       end if
    end do
  end subroutine do_bcs
end module bcs
