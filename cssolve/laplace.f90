!!<summary>Implements a modular version of the fast laplace from Babacan et al.</summary>
module laplace
use matrix_sets
  implicit none
  private
  public norm
  
  !!<usage>Create a laplace iterator object using public type laplace_iterator
  !! - Call iterator.initialize(Phi, y, sigma2, eta)
  !! - Call iterator.iterate(Js, out Js_indices, error_bars, out next_PIs, out returnsigma2).</usage>
  type, public:: laplace_iterator
     !!<member name="Phi">The measurement matrix.</member>
     !!<member name="y">The compressive sensing measurements.</member>
     real(dp), pointer :: Phi(:,:), y(:)
     !!<member name="S,Q,ls,lq">Quantities referenced in Babacan eqns 51, 52.</member>
     real(dp), allocatable :: S(:), Q(:), ls(:), lq(:)
     !!<member name="ninputs">Number of rows in PHI.</member>
     !!<member name="ncorrs">Number of columns in PHI.</member>
     !!<member name="nused">Number of basis functions included in the model.</member>
     integer ninputs, ncorrs, nused 
     !!<member name="sigma2">Initial noise variance (default : std(t)^2/1e2).</member>
     !!<member name="eta">Threshold for stopping the algorithm (default : 1e-8).</member>
     !!<member name="lambda">Parameter controlling the sparsity.</member>
     real(dp) sigma2, eta, lambda
     !!<member name="alpha">The sparse hyperparameters (1/gamma).</member>
     !!<member name="subPhi">A subset of the measurement matrix that includes basis functions
     !!included in the model.</member>
     real(dp), allocatable :: alpha(:), subPhi(:,:)
     !!<member name="indices">A reusable vector of indices computed at each iteration for 
     !!reestimation, add or delete.</member>
     !!<member name="selected">The bases in the measurement matrix that have been added
     !!to the model</member>
     integer, allocatable :: indices(:), selected(:)
     !!<member name="Sigma">The covariance matrix for the model. Diagonals represent 
     !!uncertainty in the solution vector.</member>
     !!<member name="mu">The mean of the multivariate distribution for the model.</member>
     !!<member name="theta">When maximizing the likelihood, theta is the threshold between 
     !!a functional form that has a single maximum and one that has no maximum.</member>
     real(dp), allocatable :: Sigma(:,:), mu(:), theta(:)
     !!<member name="initialized">Specifies whether the iterator has been initialized yet.</member>
     logical :: initialized = .FALSE.

     contains
       private
       procedure, public :: iterate
       procedure, public :: reset => laplace_reset
       procedure, public :: initialize => initialize_laplace
       procedure :: pre_execution
       procedure :: next_alphas
       procedure :: next_igo
       procedure :: max_l
       !final :: finalize_laplace 
       !Finalization not yet implemented in the gfortran compiler. Since we only deallocate
       !in the finalizer, there shouldn't be any memory leaks.
  end type laplace_iterator

contains
  !!<summary>Calculates the norm of the specified vector.</summary>
  !!<parameter name="vector">A real vector with arbitrary dimension.</parameter>
  pure function norm(vector)
    real(dp) norm
    real(dp), intent(in):: vector(:)
    norm = sqrt(dot_product(vector,vector))
  end function norm
  
  !!<summary>Initialize the object, calculate supporting variable values,
  !!do allocations etc.</summary>
  !!<parameter name="Phi">The measurement basis matrix.</parameter>
  !!<parameter name="y">The vector of random measurements to use in the fit.</parameter>
  !!<parameter name="sigma2">Initial noise variance (default : std(t)^2/1e2).</parameter>
  !!<parameter name="eta">Threshold for stopping the algorithm (default : 1e-8).</parameter>
  subroutine initialize_laplace(this, Phi, y, sigma2, eta)
    class(laplace_iterator) :: this
    real(dp), target, intent(in) :: Phi(:,:)
    real(dp), target, intent(in) :: y(:)
    real(dp), intent(in) :: sigma2, eta

    !Associate the pointers to the data sets
    this%Phi => Phi
    this%y => y
    this%sigma2 = sigma2
    this%eta = eta

    !Determine the sizes of the data we are working with
    this%ninputs = size(Phi,1)
    this%ncorrs = size(Phi,2)
    this%nused = 1

    if (.not. allocated(this%S)) then         
       !Allocate the arrays that we will use for all iterations
       associate (ncorrs => this%ncorrs)
         allocate(this%S(ncorrs), this%Q(ncorrs), this%ls(ncorrs), this%lq(ncorrs))
         allocate(this%theta(ncorrs))
       end associate
    end if

    associate (ncorrs => this%ncorrs, nused => this%nused, ninputs => this%ninputs)
      allocate(this%indices(nused), this%selected(nused))
      allocate(this%alpha(nused), this%subPhi(ninputs,nused))
      allocate(this%Sigma(1,1), this%mu(1))
    end associate

    !Prepare everything for iteration
    call this%pre_execution()
  end subroutine initialize_laplace

  !!<summary>Prepares the iterator to perform another iteration in a reweighted ell1
  !!minimization scheme.</summary>
  subroutine laplace_reset(this)
    class(laplace_iterator) :: this

    !Deallocate matrices that need to restart with a size of 1
    deallocate(this%Sigma, this%mu)
    deallocate(this%indices, this%selected)
    deallocate(this%alpha, this%subPhi)
  end subroutine laplace_reset

  !!<summary>Performs iterations of the fast laplace bcs until convergence.</summary>
  !!<parameter name="Js">The  solution vector to the BCS problem.</parameter>
  !!<parameter name="Js_indices">The basis function indices that the @CREF[param.Js] solution
  !!values correspond to.</parameter>
  !!<parameter name="error_bars">Diagonal terms in the covariance matrix Sigma.</parameter>
  !!<parameter name="returnsigma2">If sigma2 was not specified, this is the value that was
  !!calculated from the data.</parameter>
  subroutine iterate(this, Js, error_bars, returnsigma2, Js_indices) 
    class(laplace_iterator) :: this 
    real(dp), intent(inout) :: Js(:)
    real(dp), intent(inout) :: error_bars(:)
    real(dp), intent(inout) :: returnsigma2
    integer, pointer, optional :: Js_indices(:)

    !!<local name="maxit">The maximum number of iterations to perform.</local>
    !!<local name="deleted">Indices/basis functions that have been deleted
    !!from the model.</local>
    !!<local name="nextalphas">Original alphas are the sparse hyperparameters (1/gamma)
    !!in the paper. They are re-estimated every iteration.</local>
    !!<local name="ml">The maximum likelihood estimates for determining which basis
    !!to add to the model next.</local>
    !!<local name="bml">Vector of ml estimates for each iteration.</local>
    !!<local name="which">The special indices in @CREF[this.indices] that will be operated on.</local>
    integer :: maxit = 1000
    integer, allocatable :: deleted(:)
    real(dp), allocatable :: nextalphas(:), ml(:)
    real(dp), allocatable :: bml(:)
    integer i, idx(1), which(1)

    !allocate the local variables
    associate (NCorr => this%ncorrs, NInp => this%ninputs)
      allocate(ml(NCorr))
      allocate(bml(maxit))
      allocate(deleted(0), nextalphas(NCorr))
    end associate
    
    !Initialize the value of bml. It keeps track of the maximum likelihoods for each iteration,
    !so it makes sense to set it to zero since that shows that nothing has happened
    bml = 0  

    do i = 1, maxit
       !First, we need some values of alpha to work with, also
       !we need to calculate theta so we know whether to add, remove, etc.
       !nextalphas also calculates the value of theta (which is buried in the type)
       call this%next_alphas(nextalphas)
       !Calculate the max likelihood estimates for the re-estimate, add and remove cases
       call calc_ml(this, ml, nextalphas, deleted)
       
       !Determine which of the operations from the previous call will max the likelihood
       bml(i) = maxval(ml)
       !Check whether we have reached convergence yet
       if (i > 2) then
          if (abs(bml(i)-bml(i-1)) < abs(bml(i)-bml(1)) * this%eta) exit
       end if
       
       !Find out which index in theta corresponds to the ml.
       !'indices' holds the indices of the basis functions that the j coefficients 
       !of the solution vector point to.
       idx = maxloc(ml)
       which = maxloc(this%indices, mask=this%indices==idx(1))
   
       !Theta tells us when the function has a valid maximum. If lambda is less than theta
       !we have a maximum so adding or reestimating makes sense. Otherwise, we need to
       !remove the basis because it is not maximizing the likelihood
       if (this%theta(idx(1)) > this%lambda) then     
          if (which(1) .ne. 0) then   
             !It is already in the list of included basis functions, re-estimate.
             call reestimate(this, idx, which, nextalphas)
          else  
             ! Add a basis function, Eqns are from Tipping A.2 
             call add_basis(this, idx, nextalphas)
          endif
       else
          if (which(1) .ne. 0 .and. size(this%indices,1) > 1) then          
             !The ml basis for this iteration does not max the overall L. Delete
             !it. Eqns from Tipping A.4
             call delete_basis(this, idx, which, deleted)
          endif
       endif

       which = 0
       call appendto(this%selected, idx(1))
    enddo

    if (present(Js_indices)) then
       !Indices kept track of which bases had been added and deleted etc during the iteration
       allocate(Js_indices(this%nused))
       Js_indices = this%indices
    end if
    
    if (this%nused .ne. size(this%indices,1) ) then
       print *, "Disagreement over how many basis functions we ended up with."
       print *, this%nused, ": number of used Js"
       print *, size(this%indices), ": size of list indices"
    end if
    
    Js = 0
    !Use the mean values of the distributions as the solution vectors
    Js(this%indices) = this%mu

    !Error bars are just the diagonals of the covariance matrix
    error_bars = 0
    error_bars(this%indices) = (/ (sqrt(this%Sigma(i,i)), i =1,size(this%Sigma,1)) /)

    associate (subPhi => this%subPhi, mu => this%mu, NINp => this%ninputs, alpha => this%alpha, &
      Sigma => this%Sigma, y => this%y, Nused => this%nused)
      returnsigma2 = sum((y - matmul(subPhi, mu))**2) / & 
           (NInp - Nused + dot_product(alpha, (/ (Sigma(i,i), i = 1, size(Sigma,1)) /)))
    end associate

    !Clean up the local variables
    deallocate(nextalphas, bml, deleted, ml)
  end subroutine iterate

  !!<summary>Deallocate all of the arrays etc that were used globally by
  !!the laplace object.</summary>
  subroutine finalize_laplace(this)
     type(laplace_iterator) :: this       
     deallocate(this%alpha, this%S, this%Q, this%ls, this%lq)
     deallocate(this%theta, this%Sigma, this%mu, this%indices)
     deallocate(this%selected)
  end subroutine finalize_laplace

  !!<summary>Re-estimates the value of a coefficient for a basis that is
  !!already in the model.</summary>
  !!<parameter name="idx">The index of the the theta value that maximizes the
  !!likelihood. Because of linked array structures, this is also the index of the
  !!alpha value that we need to reestimate the expansion coefficient.</parameter>
  !!<parameter name="which">The special indices in @CREF[this.indices] that will
  !!be operated on.</parameter>
  !!<parameter name="nextalphas">Original alphas are the sparse hyperparameters (1/gamma)
  !!in the paper. They are re-estimated every iteration.</parameter>
  subroutine reestimate(this, idx, which, nextalphas)
    class(laplace_iterator), intent(inout) :: this 
    integer, intent(in) :: idx(1), which(:)
    real(dp), intent(in) :: nextalphas(:)

    !!<local name="Sigmai">The column of values in @CREF[this.Sigma] corresponding
    !!to the basis function that maximizes the log likelihood.</local>
    !!<local name="comm">A term that is common to both eq 35 and 36 Tipping.</local>
    !!<local name="Sigmaii">The value in @CREF[this.Sigma] along the diagonal for the
    !!basis function @CREF[param.which].</local>
    !!<local name="mui, delta, ki, thisalpha">Variables used in the Tipping equations
    !!for reestimating the basis function.</local>
    real(dp) :: Sigmai(this%nused,1), comm(this%ncorrs)
    real(dp) Sigmaii, mui, delta, ki, thisalpha
          
    thisalpha = nextalphas(idx(1))

    associate (Sigma => this%Sigma, mu => this%mu, alpha => this%alpha, Phi => this%Phi, &
      subPhi => this%subPhi, S => this%S, Q => this%Q)
      Sigmaii = Sigma(which(1), which(1))
      mui = mu(which(1))
      Sigmai(:,1)= Sigma(:,which(1))

      !How much did alpha for this basis change since last iteration?
      delta = thisalpha - alpha(which(1))
      !Re-estimation of these comes from Tipping Appendix
      ki = delta / (1 + Sigmaii*delta)
      mu = mu - ki * mui * Sigmai(:,1) !Eq 34 Tipping
      this%Sigma = this%Sigma - ki * matmul(Sigmai, transpose(Sigmai)) !Eq 33 Tipping
      comm = pack(matmul(transpose(Phi), matmul(subPhi,Sigmai)) / this%sigma2, .true.)
      this%S = this%S + ki * comm**2 !Eq 35
      this%Q = this%Q + ki * mui * comm !Eq 36
      !Overwrite the value of alpha for next iteration
      this%alpha(which(1)) = thisalpha
    end associate
  end subroutine reestimate

  !!<summary>Adds a basis function to the model. These eqns all from Tipping A.2.</summary>
  !!<parameter name="idx">The index of the the theta value that maximizes the
  !!likelihood. Because of linked array structures, this is also the index of the
  !!alpha value that we need to add the basis function and perform a first-time evaluation
  !!of the expansion coefficient.</parameter>
  !!<parameter name="nextalphas">Original alphas are the sparse hyperparameters (1/gamma)
  !!in the paper. They are re-estimated every iteration.</parameter>
  subroutine add_basis(this, idx, nextalphas)
    class(laplace_iterator), intent(inout) :: this 
    integer, intent(in) :: idx(1)
    real(dp), intent(in) :: nextalphas(:)

  !!<local name="thisalpha">The alpha value for the basis function being added.</local>
  !!<local name="comm1, subPhii, mui">Term first calculated when updating the Sigma matrix,
  !!but that is also needed in other routines of the calling parent subroutine @CREF[add_basis].
  !!Included as a parameter for sensible code modularization.</local>    
    real(dp) :: comm1(this%nused,1), subPhii(this%ninputs,1)
    real(dp) :: mui, thisalpha
    !----------------------------------------------------------------------

    subPhii(:,1) = this%Phi(:,idx(1))
    thisalpha = nextalphas(idx(1))     
    call add_sigma(this, idx, thisalpha, comm1, subPhii, mui)     
    
    ! Add element to mu
    this%mu = this%mu - mui*pack(comm1, .true.)
    call appendto(this%mu, mui)
    
    !Recalculate the subset matrix
    call add_subPhi(this, subPhii)
   
    ! update list of used Js
    call appendto(this%indices, idx(1))
    this%nused = size(this%indices,1)
    ! Update alphas
    call appendto(this%alpha, thisalpha)   
  end subroutine add_basis

  !!<summary>Recalculates the value of the subset of the measurement matrix for a basis add operation</summary>
  subroutine add_subPhi(this, subPhii)
    class(laplace_iterator), intent(inout) :: this 
    real(dp), intent(in) :: subPhii(:,:)

    !----------------------------------------------------------------------
    !LOCAL VARIABLES BEGIN
    !----------------------------------------------------------------------     
    real(dp), allocatable :: nextsubPhi(:,:)
    !----------------------------------------------------------------------

    associate (Phi => this%Phi, subPhi => this%subPhi, S => this%S, Q => this%Q)
      !updata phi matrix, since the size of subPhi is changing, we need a temp
      !matrix to store the additional components: nextsubPhi
      allocate(nextsubPhi(size(subPhi,1), size(subPhi,2)+1))
      nextsubPhi(:,1:size(subPhi,2)) = subPhi
      nextsubPhi(:,size(subPhi,2)+1) = subPhii(:,1)

      !Now, overwrite the value of the actual subPhi matrix
      call move_alloc(nextsubPhi, this%subPhi)
    end associate                      
  end subroutine

  !!<summary>Recalculates the value of Sigma matrix for a basis add operation.</summary>
  !!<parameter name="idx">The index of the the theta value that maximizes the
  !!likelihood.</parameter>
  !!<parameter name="thisalpha">The alpha value for the basis function being added.</parameter>
  !!<parameter name="comm1, subPhii, mui">Term first calculated when updating the Sigma matrix,
  !!but that is also needed in other routines of the calling parent subroutine @CREF[add_basis].
  !!Included as a parameter for sensible code modularization.</parameter>
  subroutine add_sigma(this, idx, thisalpha, comm1, subPhii, mui)
    class(laplace_iterator), intent(inout) :: this 
    integer, intent(in) :: idx(1)
    real(dp), intent(in) :: thisalpha
    real(dp), intent(inout) :: comm1(:,:), subPhii(:,:) !
    real(dp), intent(out) :: mui

    !----------------------------------------------------------------------
    !LOCAL VARIABLES BEGIN
    !----------------------------------------------------------------------
    !These local variables are used temporarily for calculating the equations in Tipping A.2
    !If you compare those equations with the code below, the variables' usage will be
    !obvious. We don't bother trying to describe them here.
    real(dp) :: off(this%nused)
    real(dp) :: comm2(this%ncorrs), ei(this%ninputs,1)
    real(dp) :: upper(this%nused,this%nused)
    real(dp), allocatable :: nextSigma(:,:)
    integer j, k 
    real(dp) :: Sigmaii
    !----------------------------------------------------------------------
   
    associate (subPhi => this%subPhi, Sigma => this%Sigma, S => this%S, Q => this%Q)
      Sigmaii = 1 / (thisalpha + S(idx(1)))
      mui = Sigmaii * Q(idx(1))
      comm1 = matmul(Sigma, matmul(transpose(subPhi), subPhii)) / this%sigma2
      ei = subPhii - matmul(subPhi, comm1)

      ! Update Sigma, adding a row and a column
      off = pack(-Sigmaii * comm1, .true.)
      upper = Sigma + Sigmaii * matmul(comm1, transpose(comm1))

      !Because we are changing the size of sigma, we need a temp matrix nextSigma
      allocate(nextSigma(size(Sigma,1)+1, size(Sigma,2)+1))

      ! Add row and column to Sigma matrix
      do j=1,size(Sigma,1)
         do k=1,size(Sigma,1)
            nextSigma(j,k) = upper(j,k)
         enddo
         nextSigma(j, size(Sigma,1)+1) = off(j)
      enddo

      do j=1,size(Sigma,1)
         nextSigma(size(Sigma,1) + 1, j) = off(j)
      enddo
      nextSigma(size(Sigma,1) + 1, size(Sigma,2)+1) = Sigmaii

      !Now overwrite the value of the original sigma matrix
      call move_alloc(nextSigma, this%Sigma)
    end associate

    !comm2 is a term that is common to the calculation of both S and Q        
    comm2 = pack(matmul(transpose(this%Phi), ei) / this%sigma2,.true.)
    this%S = this%S - Sigmaii*comm2**2        
    this%Q = this%Q - mui*comm2
  end subroutine

  !!<summary>Removes a basis function from the list</summary>
  !!<parameter name="idx">The index of the the theta value that maximizes the likelihood.</parameter>
  !!<parameter name="which">The special indices in @CREF[this.indices] that will
  !!be operated on.</parameter>
  !!<parameter name="deleted">The list of basis function indices that have been deleted
  !!during the iterations of the algorithm.</parameter>
  subroutine delete_basis(this, idx, which, deleted)
    class(laplace_iterator), intent(inout) :: this 
    integer, intent(in) :: idx(1), which(:)
    integer, allocatable, intent(inout) :: deleted(:)

    !!<local name="Sigmai">The column of values in @CREF[this.Sigma] corresponding
    !!to the basis function that maximizes the log likelihood.</local>
    !!<local name="Sigmaii">The value in @CREF[this.Sigma] along the diagonal for the
    !!basis function @CREF[param.which].</local>
    real(dp) :: Sigmai(this%nused,1)
    real(dp) Sigmaii, mui

    !Recalculate the value of the Sigma matrix
    call appendto(deleted, idx(1))
    call delete_sigma(this, which, Sigmai, Sigmaii, mui)

    !Recalculate the value of mu
    this%mu = this%mu - pack(mui / Sigmaii*Sigmai,.true.)
    !Delete an entry from the mu vector
    call delete(this%mu, which(1))

    !Delete entries from the list of selected bases and corresponding alpha = (1 / gamma)
    call delete(this%indices, which(1))
    this%nused = size(this%indices,1)
    call delete(this%alpha, which(1))

    ! Delete column from subPhi matrix
    call delete_subPhi(this, which)
  end subroutine delete_basis

  !!<summary>Recalculates the subset of the measurement matrix for a basis delete operation</summary>
  subroutine delete_subPhi(this, which)
    class(laplace_iterator), intent(inout) :: this       
    integer, intent(in) :: which(:)

    !----------------------------------------------------------------------
    !LOCAL VARS START
    !----------------------------------------------------------------------
    real(dp), allocatable ::  nextsubPhi(:,:)
    integer dx, j
    !----------------------------------------------------------------------

    associate(subPhi => this%subPhi)
      !Next subPhi is the temp matrix that we need to change the size of the
      !subPhi matrix for the iterator
      allocate(nextsubPhi(size(subPhi,1), size(subPhi,2)-1))
      dx = 1
      do j=1,size(subPhi,2)
         if (j == which(1)) then
            cycle
         else
            nextsubPhi(:,dx) = subPhi(:,j)
            dx = dx + 1
         end if
      end do
    end associate
      
    !Overwrite the value of the subPhi matrix
    call move_alloc(nextsubPhi, this%subPhi)
  end subroutine

  !!<summary>Recalculates the value of the Sigma matrix for a basis delete operation</summary>
  !!<parameter name="which">The special indices in @CREF[this.indices] that will
  !!be operated on.</parameter>
  !!<parameter name="Sigmai">The column of values in @CREF[this.Sigma] corresponding
  !!to the basis function that maximizes the log likelihood.</parameter>
  !!<parameter name="Sigmaii">The value in @CREF[this.Sigma] along the diagonal for the
  !!basis function @CREF[param.which].</parameter>
  !!<parameter name="mui">Term first calculated when updating the Sigma matrix,
  !!but that is also needed in other routines of the calling parent subroutine
  !!@CREF[delete_basis].</parameter>
  subroutine delete_sigma(this, which, Sigmai, Sigmaii, mui)
    class(laplace_iterator), intent(inout) :: this 
    integer, intent(in) :: which(:)
    real(dp), intent(inout) :: Sigmai(:,:)
    real(dp), intent(out) :: Sigmaii, mui

    !----------------------------------------------------------------------
    !LOCAL VARS START
    !----------------------------------------------------------------------
    !!<local name="nextSigma">The temp matrix used to increase the size of Sigma.</local>
    !!<local name="comm">Term used repeatedly in the Tipping equations.</local>
    !!<local name="dx, dy, j, k">Iteration variables.</local>
    real(dp), allocatable :: nextSigma(:,:)
    real(dp) :: comm(this%ncorrs)
    integer dx, dy, j, k
    !----------------------------------------------------------------------

    associate (Sigma => this%Sigma)
      !Recalculate all the values as suggested by Tipping in A.4
      Sigmaii = Sigma(which(1),which(1))
      mui = this%mu(which(1))
      Sigmai(:,1) = Sigma(:,which(1))
      this%Sigma = this%Sigma - matmul(Sigmai, transpose(Sigmai))/Sigmaii

      !Delete a row and column from the Sigma matrix. nextSigma is 
      allocate(nextSigma(size(Sigma,1)-1,size(Sigma,2)-1))
      dx = 1
      do j=1,size(Sigma,1)
         dy = 1
         if (j == which(1)) cycle
         do k=1,size(Sigma,2)
            if (k == which(1)) cycle
            nextSigma(dx,dy) = Sigma(j,k)
            dy = dy + 1
         enddo
         dx = dx + 1
      end do
    end associate

    !Overwrite the value of the Sigma matrix
    call move_alloc(nextSigma, this%Sigma)

    associate(Phi => this%Phi, subPhi => this%subPhi, &
      sigma2 => this%sigma2)
      comm = pack(matmul(transpose(Phi),matmul(subPhi,Sigmai))/sigma2,.true.)
      this%S = this%S + comm**2/Sigmaii
      this%Q = this%Q + mui/Sigmaii*comm
    end associate
  end subroutine

  !!<summary>Calculates the max likelihood estimates for the re-estimate,
  !!add and remove cases.</summary>
  !!<parameter name="ml">The maximum likelihood estimates for determining which basis
  !!to add to the model next.</parameter>
  !!<parameter name="nextalphas">Original alphas are the sparse hyperparameters (1/gamma)
  !!in the paper. They are re-estimated every iteration.</parameter>
  !!<parameter name="deleted">Indices/basis functions that have been deleted
  !!from the model.</parameter>
  subroutine calc_ml(this, ml, nextalphas, deleted)
    class(laplace_iterator), intent(in) :: this 
    real(dp), intent(inout) :: ml(:)
    integer, allocatable, intent(in) :: deleted(:)
    real(dp), intent(in) :: nextalphas(:)
    !----------------------------------------------------------------------
    integer, allocatable :: ire(:), iad(:), add(:), is0(:), ide(:)
    real(dp), allocatable ::  Nalpha(:)
    integer, allocatable :: ig0(:), which(:)
    integer tempvec(1), j

    !allocate the local variables
    allocate(ire(this%ncorrs), iad(this%ncorrs), add(this%ncorrs), is0(this%ncorrs), ide(this%ncorrs))

    !Initialize distribution to large negative value
    ml = -1e10

    !Now look for the indices for reestimation; [C,IA,IB] = intersect(A,B) also returns 
    !index vectors IA and IB such that C = A(IA) and C = B(IB). 
    allocate(ig0(this%ncorrs))
    call this%next_igo(ig0)
    call intersection(ig0, this%indices, ire)

    !We still need to find the index for the intersection values returned in the previous call
    allocate(which(size(ire,1)))
    do j=1,size(ire,1)
       tempvec = maxloc(this%indices, mask=this%indices==ire(j))
       which(j) = tempvec(1)
    enddo     
    
    !Common index between non-zero gamma_i and basis we are working on.
    !According to algorithm, this means re-estimate. Find the indices to
    !work with.
    if (size(ire,1) .ne. 0) then         
       allocate(Nalpha(size(ire,1)))
       Nalpha = nextalphas(ire)
       !Make sure that we have valid values for computing the max likelihood
       if (any(Nalpha / (Nalpha + this%ls(ire)) < 0)) then
          stop 'ERROR:(negative argument to log 1):  You most likely made a very poor choice for sigma2'
       endif
       if (any(this%alpha(which) / (this%alpha(which) + this%ls(ire)) < 0)) stop 'ERROR:(negative argument to log 2)'

       !Compute the max likelihood
       ml(ire) = this%max_l(Nalpha, ire) - this%max_l(this%alpha(which), ire)
       deallocate(Nalpha)
    end if
    deallocate(which)

   !indices for adding; C = setdiff(A,B) for vectors A and B, returns 
   !the values in A that are not in B with no repetitions. This means we
   !will get the indices for bases with theta larger than lambda that
   !have gamma_i = 0.
   call complement(ig0, ire, iad)
   if (size(iad,1) .ne. 0) then
      allocate(Nalpha(size(iad,1)))
      Nalpha = nextalphas(iad)
      if (any(Nalpha/(Nalpha + this%ls(iad)) < 0)) stop 'negative arguement to log 3'
      !Calculate the max likel. again
      ml(iad) = this%max_l(Nalpha, iad)
      call intersection(deleted, iad, add)
      ml(add) = -1e10
      deallocate(Nalpha)
   end if
   
   !Find the indices for deletion
   call complement((/ (j,j=1,this%ncorrs) /), ig0, is0)
   call intersection(is0, this%indices, ide)
   allocate(which(size(ide,1)))
   do j=1,size(ide,1)
      tempvec = maxloc(this%indices, mask=this%indices==ide(j))
      which(j) = tempvec(1)
   enddo

   if (size(ide,1) .ne. 0) then
      if (size(this%indices,1) == 1) then
         ml(ide) = -1e10
      else
         if (any(this%alpha(which)/(this%alpha(which) + this%ls(ide)) < 0)) stop 'negative arguement to log 4'
         !Calculate the max likel. again
         ml(ide) = -this%max_l(this%alpha(which), ide)
      end if        
   end if

   deallocate(which)
   deallocate(ig0)
   deallocate(ire, iad, add, is0, ide)
  end subroutine calc_ml

  !!<summary>Determines the maximum likelihood for the specified index and alphas.</summary>
  !!<parameter name="alpha">The list of alpha values for the current index being processed
  !!in the parent subroutine @CREF[calc_ml].</parameter>
  !!<parameter name="index">The index of the list of @CREF[param.alpha] in the overall
  !!list of alpha values.</parameter>
  function max_l(this, alpha, index)
    class(laplace_iterator) :: this 
    integer, intent(in) :: index(:)
    real(dp), intent(in) :: alpha(:)
    real(dp) :: max_l(size(alpha, 1))

    associate (lq => this%lq, ls => this%ls, lambda => this%lambda)
      max_l = lq(index)**2 / (alpha + ls(index)) + log(alpha / (alpha + ls(index))) - lambda / alpha
    end associate
    return
  end function max_l

  !!<summary>Calculates the next set of indices that satisfy the condition for nonzero gamma_i
  !!(i.e. ig) in eq 50</summary>
  subroutine next_igo(this, ig0)
    class(laplace_iterator) :: this  
    integer, allocatable, intent(inout) :: ig0(:)
    integer, allocatable :: tempig0(:)
    integer dx, numnonzero, k

    ig0 = 0
    dx = 1
    
    associate (theta => this%theta, lambda => this%lambda)
      do k = 1, size(theta, 1)
         if (theta(k) > lambda) then
            ig0(dx) = k
            dx = dx + 1
         endif
      enddo
    end associate

    numnonzero = count(ig0 > 0)
    allocate(tempig0(numnonzero))
    tempig0 = pack(ig0, mask=ig0 > 0)

    !Overwrite the value of ig0 with the new one
    call move_alloc(tempig0, ig0)
  end subroutine next_igo

  !!<summary>Calculates the next set of alpha values to use in the current iteration.</summary>
  subroutine next_alphas(this, nextalphas)
    class(laplace_iterator), intent(inout) :: this
    real(dp), intent(inout) :: nextalphas(:)
    real(dp) :: A(this%ncorrs), B(this%ncorrs), C(this%ncorrs), discriminant(this%ncorrs)

    associate (S => this%S, Q => this%Q, ls => this%ls, lq => this%lq, lambda => this%lambda, &
         indices => this%indices, alpha => this%alpha, theta => this%theta)
      this%nused = size(indices,1)
      ls = S
      lq = Q
      ls(indices) = alpha * S(indices)/abs(alpha-S(indices)) !Eq 53
      lq(indices) = alpha * Q(indices)/abs(alpha-S(indices)) !Eq 54

      !We don't allow the user to specify the lambda. Use the one defined in the paper
      !Using equation (35); initially, lambda = 0, which means that nu=0 by solving
      !eq. (37) with psi being the PolyGamma[0, x] function.
      lambda = 2 *(this%nused-1)/sum(1/alpha)

      !See the numerator in eq 46, we are maximizing that by setting it to 0,
      !use the quadratic eqn.

      !A is the only term that can cause complex roots => eqn failure.
      !So we require s - q^2 < 0 or q^2 > s. See note and plots in Tipping
      !near eq 20.
      A = lambda * ls**2
      B = 2*lambda*ls+ls**2
      C = lambda + ls -lq**2

      !See eq 50, theta is the condition used to decide how to change the
      !basis.
      theta = lq**2-ls
    end associate

    !The paper says which root to keep. Since alpha = Lambda =
    !diag(gamma_i), we can calculate some new ones and then compare it to
    !our theta.
    discriminant = B**2-4*A*C

    !We add this line for numerical stability with the additional penalty functions
    where (discriminant < 0)
       discriminant = 0
    end where

    if (this%lambda /= 0 ) then 
       nextalphas = (2*A)/(-B + sqrt(discriminant)) !Quadratic solution
    else 
       nextalphas = this%ls**2/(this%lq**2 - this%ls) !Linear solution
    endif
  end subroutine next_alphas

  !!<summary>Calculates the initial values of the matrices etc that are required before any iterations
  !!can take place</summary>
  subroutine pre_execution(this)
    class(laplace_iterator) :: this
    real(dp) :: tempmat(this%ncorrs,this%ninputs)
    real(dp) :: Phidoty(this%ncorrs), Phisquared(this%ncorrs), ratio(this%ncorrs), left(this%ncorrs,1)
    real(dp) maxr, hessian
    integer :: i

    !Coefficient values in measurement basis.
    Phidoty = matmul(transpose(this%Phi), this%y)
    tempmat = transpose(this%Phi**2)

    !Returns the sum over each column of phi^2.
    Phisquared = (/(sum(tempmat(i,:)), i=1,this%ncorrs) /)

    !How large are the measured coefficients in phi basis compared to the
    !coefficients of the basis function alone.
    ratio = Phidoty**2/Phisquared
    maxr = maxval(ratio)
    this%indices = maxloc(ratio)
    this%selected = this%indices
          
    !Alpha is Lambda in the paper; Lambda = diag(1/gamma_i)
    !Since we only have a 1x1, the diag amounts to a single value.
    this%alpha(1) = Phisquared(this%indices(1))/(maxr-this%sigma2)
          
    !We selected the basis with the maximum variance to init alpha, so we want
    !our subset of phi to only include that single basis.
    !subPhi is the subset of the basis, Phi in the paper
    this%subPhi(:,1) = (/ (this%Phi(i,this%indices(1)) , i=1,this%ninputs) /)
    !Next 2 lines for eq 25
    hessian = this%alpha(1) + dot_product(this%subPhi(:,1),this%subPhi(:,1))/this%sigma2
    this%Sigma(1,1) = 1/hessian
    !Babcan eq 24
    this%mu(1) = this%Sigma(1,1) * Phidoty(this%indices(1))/this%sigma2

    !Eq 51 and 52, part of 2nd term on RHS
    associate (sigma2 => this%sigma2)
      left = matmul(transpose(this%Phi),this%subPhi)/sigma2
      this%S = Phisquared/sigma2 - this%Sigma(1,1)*left(:,1)**2 !Babacan eq 51
      this%Q = Phidoty/sigma2 - this%Sigma(1,1)*Phidoty(this%indices(1))/sigma2*left(:,1) !Babacan eq 52
    end associate
  end subroutine pre_execution
end module laplace
