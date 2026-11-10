module kinetics
  use kinds, only: dp
  use species, only: ns
  use stoichiometry, only: build_stoichiometric_matrix
  use arrhenius, only: rate_constant_from_log_prefactor
  implicit none
  private
  public :: kinetic_parameters, initialize_kinetics, reaction_rates, species_sources

  type :: kinetic_parameters
    real(dp) :: k_forward(ns) = 0.0_dp
    real(dp) :: k_reverse(ns) = 0.0_dp
    logical :: reversible(ns) = .false.
    real(dp) :: acid_activity = 1.0_dp
    real(dp) :: concentration_reference = 1.0_dp
  end type kinetic_parameters
contains
  pure function initialize_kinetics(temperature, ln_k0_f, e_f, ln_k0_r, e_r, reversible, acid, cref) result(p)
    real(dp), intent(in) :: temperature, ln_k0_f(ns), e_f(ns), ln_k0_r(ns), e_r(ns)
    logical, intent(in) :: reversible(ns)
    real(dp), intent(in) :: acid, cref
    type(kinetic_parameters) :: p
    integer :: j
    p%reversible = reversible
    p%acid_activity = acid
    p%concentration_reference = cref
    do j = 1, ns
      p%k_forward(j) = rate_constant_from_log_prefactor(ln_k0_f(j), e_f(j), temperature)
      if (reversible(j)) p%k_reverse(j) = rate_constant_from_log_prefactor(ln_k0_r(j), e_r(j), temperature)
    end do
    p%reversible(15) = .false.
    p%k_reverse(15) = 0.0_dp
  end function initialize_kinetics

  pure subroutine reaction_rates(c, p, w)
    real(dp), intent(in) :: c(ns)
    type(kinetic_parameters), intent(in) :: p
    real(dp), intent(out) :: w(ns)
    real(dp) :: y(ns), h
    y = c / p%concentration_reference
    h = p%acid_activity

    w(1) = p%k_forward(1)*y(1) - p%k_reverse(1)*y(7)*y(4)
    w(2) = p%k_forward(2)*y(7)*h - p%k_reverse(2)*y(8)
    w(3) = p%k_forward(3)*y(8) - p%k_reverse(3)*y(9)
    w(4) = p%k_forward(4)*y(8) - p%k_reverse(4)*y(10)
    w(5) = p%k_forward(5)*y(9) - p%k_reverse(5)*y(10)
    w(6) = p%k_forward(6)*y(8)
    w(7) = p%k_forward(7)*y(9)
    w(8) = p%k_forward(8)*y(10)
    w(9) = p%k_forward(9)*y(9) - p%k_reverse(9)*y(13)*h
    w(10) = p%k_forward(10)*y(10) - p%k_reverse(10)*y(14)*h
    w(11) = p%k_forward(11)*y(12)
    w(12) = p%k_forward(12)*y(11)*y(4)
    w(13) = p%k_forward(13)*y(13)*y(4) - p%k_reverse(13)*y(2)
    w(14) = p%k_forward(14)*y(14)*y(4) - p%k_reverse(14)*y(3)
    w(15) = p%k_forward(15)*y(15)*y(4)

    ! Convert legacy dy/dminute rates to concentration source rates in SI.
    w = w * p%concentration_reference / 60.0_dp
  end subroutine reaction_rates

  pure subroutine species_sources(c, p, source, w)
    real(dp), intent(in) :: c(ns)
    type(kinetic_parameters), intent(in) :: p
    real(dp), intent(out) :: source(ns)
    real(dp), intent(out), optional :: w(ns)
    real(dp) :: rates(ns), nu(ns,ns)
    nu = build_stoichiometric_matrix()
    call reaction_rates(c, p, rates)
    source = matmul(nu, rates)
    if (present(w)) w = rates
  end subroutine species_sources
end module kinetics
