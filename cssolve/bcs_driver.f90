MODULE bcs_driver
  use bcs
implicit none

public bcs_fit
public f90wrap_do_wrapped

!integer, public, parameter :: dp=8

contains


subroutine bcs_fit(nfits, nholdout_, reweight_, penaltyfxn_, jcutoff_, sigma2_, eta_, full_pi, y, &
                   Ny, Nj, fit_rms, fit_err, hold_rms_, hold_err_, error_bars, js)
  use num_types, only: dp
  implicit none

  integer Ny, Nj
  integer :: nfits
  real(dp) :: js(nfits,Nj)
  real(dp):: sigma2_(nfits), eta_
  logical :: reweight_
  character(len=6) :: penaltyfxn_
  integer :: nholdout_
  real(dp):: full_pi(Ny, Nj)
  real(dp):: fit_err(nfits), fit_rms(nfits)
  ! Variable 'seed' was set to be ignored.
  real(dp):: jcutoff_
  real(dp):: hold_rms_(nfits), hold_err_(nfits)
  real(dp):: y(Ny)
  real(dp):: error_bars(nfits,Nj)
!f2py intent(in) nfits, nholdout_, reweight_, penaltyfxn_, jcutoff_, full_pi, y
!f2py intent(out) fit_rms, fit_err, hold_rms_, hold_err_, error_bars, js
!f2py integer, intent(hide),depend(full_pi) :: Ny=shape(full_pi,0),Nj=shape(full_pi,1)


  real(dp) :: fpy_start, fpy_end, fpy_elapsed = 0

  fit_rms = 0
  fit_err = 0
  js = 0
  error_bars = 0
  hold_rms_ = 0
  hold_err_ = 0

  call cpu_time(fpy_start)
  !sigma2_=5
  call do_bcs(full_pi, y, nfits, js, error_bars, fit_rms, fit_err, hold_rms_=hold_rms_, &
               hold_err_=hold_err_, nholdout_=nholdout_, jcutoff_=jcutoff_, &
               penaltyfxn_=penaltyfxn_, reweight_=reweight_,sigma2_=sigma2_,eta_=eta_)
  call cpu_time(fpy_end)
  fpy_elapsed = fpy_elapsed + fpy_end - fpy_start

end subroutine bcs_fit


subroutine f90wrap_do_wrapped(full_pi, y, sigma2, eta, js, error_bars, n0, n1)
    implicit none

    real(8), intent(in), dimension(n0,n1) :: full_pi
    real(8), intent(in), dimension(n0) :: y
    real(8), intent(in) :: sigma2
    real(8), intent(in) :: eta
    real(8), intent(out), dimension(n1) :: js
    real(8), intent(out), dimension(n1) :: error_bars
    integer :: n0
    !f2py intent(hide), depend(full_pi) :: n0 = shape(full_pi,0)
    integer :: n1
    !f2py intent(hide), depend(full_pi) :: n1 = shape(full_pi,1)
    call do_wrapped(full_pi=full_pi, y=y, sigma2=sigma2, eta=eta, js=js, &
        error_bars=error_bars)
end subroutine f90wrap_do_wrapped

end module bcs_driver


! program test
!   use bcs
!   implicit none
!   integer, parameter :: Ny=1290, Nj=121,nfits=5
!   real(8) full_pi(Ny,nj), y(Ny), js(nj), error_bars(nj),
! fit_rms(nfits),fit_err(nfits)
  
!   open(10, file='Amat.txt',action='read',status='old')
!   read(10, *) full_pi
!   close(10)
!   open(10, file='y.txt',action='read',status='old')
!   read(10, *) y
!   close(10)
!   call do_bcs(full_pi, y, nfits, js, error_bars, fit_rms, fit_err)
!   print *, 'j=',js
!   print *, 'fit_rms', fit_rms
!   print *, 'fit_err', fit_err
! end program test
