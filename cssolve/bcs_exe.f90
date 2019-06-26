!!<summary>Auto-generated unit test for bcs.do_bcs
!!using FORTPY. Generated on 2015-08-25 14:40:26.868183.
!!Compile a driver for the full BCS.</summary>
PROGRAM UNITTEST_do_bcs
  use num_types, only: dp
  use bcs
  use fortpy
  implicit none

  ! Variable 'eta_' was set to be ignored.
  real(dp), allocatable :: js(:,:)
  ! Variable 'sigma2_' was set to be ignored.
  real(dp), allocatable :: fit_rms(:)
  logical :: reweight_
  character(len=6) :: penaltyfxn_
  integer :: nholdout_
  real(dp), allocatable :: full_pi(:, :)
  integer :: nfits
  real(dp), allocatable :: fit_err(:)
  ! Variable 'seed' was set to be ignored.
  real(dp) :: jcutoff_
  real(dp), allocatable :: hold_rms_(:)
  real(dp), allocatable :: y(:)
  real(dp), allocatable :: error_bars(:,:)
  real(dp), allocatable :: hold_err_(:)

  real(fdp) :: fpy_start, fpy_end, fpy_elapsed = 0

  nfits = 3
  nholdout_ = 15
  allocate(fit_rms(nfits))
  fit_rms = 0
  allocate(fit_err(nfits))
  fit_err = 0
  call fpy_read('pi.in', '#', full_pi)
  call fpy_read('y.in', '#', y)
  allocate(js(nfits, size(full_pi, 2)))
  js = 0
  allocate(error_bars(nfits, size(full_pi, 2)))
  error_bars = 0
  jcutoff_ = 1e-3
  penaltyfxn_ = 'arctan'
  reweight_ = .true.
  allocate(hold_rms_(nfits))
  hold_rms_ = 0
  allocate(hold_err_(nfits))
  hold_err_ = 0

  call cpu_time(fpy_start)
  call do_bcs(full_pi, y, nfits, js, error_bars, fit_rms, fit_err, hold_rms_=hold_rms_, &
               hold_err_=hold_err_, nholdout_=nholdout_, jcutoff_=jcutoff_, &
               penaltyfxn_=penaltyfxn_, reweight_=reweight_)
  call cpu_time(fpy_end)
  fpy_elapsed = fpy_elapsed + fpy_end - fpy_start
  call pysave(js, 'js.out')
  call pysave(error_bars, 'ebars.out')
  call pysave(fit_rms, 'rms.out')
  call pysave(fit_err, 'err.out')
  call pysave(hold_rms_, 'holdrms.out')
  call pysave(hold_err_, 'holderr.out')

  call pysave(fpy_elapsed, 'fpytiming.out')
END PROGRAM UNITTEST_do_bcs