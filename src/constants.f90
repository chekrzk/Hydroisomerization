module constants
  use kinds, only: dp
  implicit none
  private
  public :: gas_constant, pi
  real(dp), parameter :: gas_constant = 8.314462618_dp
  real(dp), parameter :: pi = acos(-1.0_dp)
end module constants
