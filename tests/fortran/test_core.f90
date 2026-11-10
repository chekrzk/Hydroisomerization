program test_core
  use kinds, only: dp
  use constants, only: pi
  use species, only: ns
  use stoichiometry, only: build_stoichiometric_matrix
  use arrhenius, only: rate_constant_from_log_prefactor
  use kinetics, only: kinetic_parameters, initialize_kinetics
  use knudsen, only: effective_knudsen_diffusivity
  use fvm_1d, only: radial_mesh, build_uniform_cylindrical_mesh
  use pellet_1d, only: initialize_pellet_model, pellet_rhs
  use effectiveness, only: cylindrical_effectiveness, cylindrical_first_order_fvm_eta
  implicit none
  real(dp) :: nu(ns,ns), k, deff, rhs(ns*4), state(ns*4)
  real(dp) :: zeros(ns), bulk(ns), diff(ns)
  logical :: reversible(ns)
  type(radial_mesh) :: mesh
  type(kinetic_parameters) :: p

  nu = build_stoichiometric_matrix()
  call assert_close(nu(15,15), -1.0_dp, 0.0_dp, 'R15 consumes C15')
  call assert_close(nu(4,15), -1.0_dp, 0.0_dp, 'R15 consumes H2')
  call assert_close(nu(6,15), 1.0_dp, 0.0_dp, 'R15 produces C6')

  k = rate_constant_from_log_prefactor(45.58_dp, 52.2_dp*4184.0_dp, 593.15_dp)
  call assert_close(k, 3.648461506_dp, 1.0e-9_dp, 'Arrhenius reference value')
  deff = effective_knudsen_diffusivity(2.692307692e-9_dp, 593.15_dp, 0.22644_dp, &
                                       0.435028249_dp, 3.0_dp)
  call assert_close(deff, 6.12946375e-8_dp, 3.0e-7_dp, 'Knudsen reference value')

  mesh = build_uniform_cylindrical_mesh(20, 1.5e-3_dp)
  call assert_close(sum(mesh%volume), pi*(1.5e-3_dp)**2, 1.0e-14_dp, 'FVM volumes conserve cylinder volume')
  call assert_close(mesh%area(0), 0.0_dp, 0.0_dp, 'Center face area is zero')
  call assert_close(cylindrical_effectiveness(0.0_dp), 1.0_dp, 1.0e-14_dp, 'eta(phi=0)=1')
  call assert_close(cylindrical_first_order_fvm_eta(2.0_dp, 200), &
                    cylindrical_effectiveness(2.0_dp), 2.0e-5_dp, &
                    'FVM effectiveness matches analytical cylinder solution')

  zeros = 0.0_dp
  reversible = .false.
  reversible(15) = .true.
  p = initialize_kinetics(593.15_dp, zeros, zeros, zeros, zeros, reversible, 1.0_dp, 1.0_dp)
  call assert_true(.not. p%reversible(15), 'R15 reversibility is forcibly disabled')
  call assert_close(p%k_reverse(15), 0.0_dp, 0.0_dp, 'R15 reverse rate is absent')

  bulk = 2.0_dp
  diff = 1.0e-7_dp
  state = 2.0_dp
  call initialize_pellet_model(4, 1.0e-3_dp, diff, bulk, p, enable_reaction=.false.)
  call pellet_rhs(state, rhs)
  call assert_close(maxval(abs(rhs)), 0.0_dp, 1.0e-14_dp, 'Uniform no-reaction state has zero RHS')
  write(*,'(a)') 'All core tests passed.'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAILED: '//message
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    call assert_true(abs(actual-expected) <= max(tolerance, tolerance*abs(expected)), message)
  end subroutine assert_close
end program test_core
