module arrhenius
  use kinds, only: dp
  use constants, only: gas_constant
  implicit none
  private
  public :: rate_constant_from_log_prefactor
contains
  elemental pure function rate_constant_from_log_prefactor(ln_k0, activation_energy, temperature) result(k)
    real(dp), intent(in) :: ln_k0, activation_energy, temperature
    real(dp) :: k
    k = exp(ln_k0 - activation_energy / (gas_constant * temperature))
  end function rate_constant_from_log_prefactor
end module arrhenius
