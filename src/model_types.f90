module model_types
  use kinds, only: dp, i8
  implicit none
  private
  public :: simulation_config, solver_stats

  integer, parameter :: ns_cfg = 15

  type :: simulation_config
    character(len=32) :: run_mode = 'demo'
    character(len=64) :: run_name = 'v1'
    character(len=16) :: simulation_mode = 'transient'
    integer :: nr = 50
    real(dp) :: temperature = 593.15_dp
    real(dp) :: pressure = 3.0e6_dp
    real(dp) :: radius = 1.5e-3_dp
    real(dp) :: porosity = 0.435028249_dp
    real(dp) :: tortuosity = 3.0_dp
    real(dp) :: pore_radius = 2.692307692e-9_dp
    real(dp) :: acid_activity = 1.0_dp
    real(dp) :: concentration_reference = 608.306609_dp
    real(dp) :: t_start = 0.0_dp
    real(dp) :: t_end = 3600.0_dp
    real(dp) :: output_interval = 60.0_dp
    real(dp) :: rtol = 1.0e-9_dp
    real(dp) :: atol = 1.0e-12_dp
    real(dp) :: initial_step = 1.0e-6_dp
    real(dp) :: max_step = 15.0_dp
    real(dp) :: steady_tolerance = 1.0e-10_dp
    real(dp) :: steady_maximum_time = 1.0e7_dp
    integer(i8) :: max_steps = 100000_i8
    real(dp) :: negative_warning = 1.0e-11_dp
    real(dp) :: negative_error = 1.0e-8_dp
    real(dp) :: molar_mass(ns_cfg) = 0.0_dp
    real(dp) :: c_bulk(ns_cfg) = 0.0_dp
    real(dp) :: c_initial(ns_cfg) = 0.0_dp
    real(dp) :: ln_k0_forward(ns_cfg) = 0.0_dp
    real(dp) :: activation_forward(ns_cfg) = 0.0_dp
    real(dp) :: ln_k0_reverse(ns_cfg) = 0.0_dp
    real(dp) :: activation_reverse(ns_cfg) = 0.0_dp
    logical :: reversible(ns_cfg) = .false.
  end type simulation_config

  type :: solver_stats
    integer(i8) :: steps = 0_i8
    integer(i8) :: rhs_evaluations = 0_i8
    integer(i8) :: error_test_failures = 0_i8
    integer(i8) :: nonlinear_iterations = 0_i8
    integer(i8) :: nonlinear_failures = 0_i8
  end type solver_stats
end module model_types
