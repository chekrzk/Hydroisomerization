program hydroiso_v1
  use, intrinsic :: iso_c_binding, only: c_double
  use kinds, only: dp
  use species, only: ns
  use model_types, only: simulation_config, solver_stats
  use config_io, only: read_simulation_config
  use kinetics, only: kinetic_parameters, initialize_kinetics
  use stoichiometry, only: build_stoichiometric_matrix
  use knudsen, only: effective_knudsen_diffusivity
  use pellet_1d, only: initialize_pellet_model, pellet_mesh, volume_average_concentration, &
                       volume_average_rates, surface_fluxes, minimum_concentration, pellet_rhs
  use cvode_backend, only: cvode_solver, create_cvode_solver
  use output_io, only: output_files, open_outputs, write_snapshot, write_statistics, close_outputs
  implicit none

  type(simulation_config) :: cfg
  type(kinetic_parameters) :: kinetic_model
  type(cvode_solver) :: solver
  type(solver_stats) :: stats
  type(output_files) :: files
  character(len=512) :: config_path, output_directory
  real(dp) :: diffusivity(ns), average(ns), average_rates(ns), surface_rates(ns), flux(ns)
  real(dp) :: average_sources(ns), balance_residuals(ns), average_rhs(ns)
  real(dp) :: derivative_average(ns), nu(ns,ns)
  real(dp) :: current_time, target_time, cmin
  real(c_double), allocatable, target :: state(:)
  real(dp), allocatable :: derivative_state(:)
  integer :: i, j, k, narg
  logical :: steady_converged

  narg = command_argument_count()
  if (narg /= 2) then
    write(*,'(a)') 'Usage: hydroiso_v1 <generated_config.nml> <output_directory>'
    error stop 2
  end if
  call get_command_argument(1, config_path)
  call get_command_argument(2, output_directory)
  call read_simulation_config(trim(config_path), cfg)

  if (trim(cfg%run_mode) == 'demo') then
    write(*,'(a)') 'WARNING: RESULTS ARE NOT EXPERIMENTALLY VALIDATED'
  end if
  if (cfg%nr < 2) error stop 'radial_cells must be at least 2.'
  if (cfg%radius <= 0.0_dp) error stop 'pellet radius must be positive.'
  if (cfg%tortuosity <= 0.0_dp) error stop 'tortuosity must be positive.'

  do i = 1, ns
    diffusivity(i) = effective_knudsen_diffusivity(cfg%pore_radius, cfg%temperature, &
                                                   cfg%molar_mass(i), cfg%porosity, cfg%tortuosity)
  end do
  kinetic_model = initialize_kinetics(cfg%temperature, cfg%ln_k0_forward, cfg%activation_forward, &
                                      cfg%ln_k0_reverse, cfg%activation_reverse, cfg%reversible, &
                                      cfg%acid_activity, cfg%concentration_reference)
  call initialize_pellet_model(cfg%nr, cfg%radius, diffusivity, cfg%c_bulk, kinetic_model)
  nu = build_stoichiometric_matrix()

  allocate(state(ns*cfg%nr))
  allocate(derivative_state(ns*cfg%nr))
  do j = 1, cfg%nr
    k = (j-1)*ns
    state(k+1:k+ns) = cfg%c_initial
  end do

  call open_outputs(trim(output_directory), files)
  call calculate_observables()
  call write_current_snapshot(cfg%t_start)
  solver = create_cvode_solver(state, cfg)
  current_time = cfg%t_start
  steady_converged = .false.
  if (trim(cfg%simulation_mode) == 'steady') then
    cfg%t_end = cfg%steady_maximum_time
  end if
  target_time = min(cfg%t_start + cfg%output_interval, cfg%t_end)
  do while (current_time < cfg%t_end - 10.0_dp*epsilon(cfg%t_end))
    call solver%advance(target_time, current_time)
    call validate_nonnegative(current_time)
    call calculate_observables()
    call write_current_snapshot(current_time)
    if (trim(cfg%simulation_mode) == 'steady') then
      if (maxval(abs(derivative_state)) < cfg%steady_tolerance) then
        steady_converged = .true.
        exit
      end if
    end if
    target_time = min(target_time + cfg%output_interval, cfg%t_end)
  end do
  call solver%statistics(stats)
  cmin = minimum_concentration(state)
  call write_statistics(trim(output_directory), stats, cmin)
  call solver%finalize()
  call close_outputs(files)
  write(*,'(a,i0)') 'CVODE steps: ', stats%steps
  write(*,'(a,es12.4)') 'Minimum concentration: ', cmin
  if (trim(cfg%simulation_mode) == 'steady') then
    write(*,'(a,l1)') 'Steady state converged: ', steady_converged
    if (.not. steady_converged) error stop 'Steady-state tolerance was not reached.'
  end if
  write(*,'(a)') 'V1 calculation completed.'
contains
  subroutine write_current_snapshot(time)
    real(dp), intent(in) :: time
    call write_snapshot(files, time, pellet_mesh, state, average, average_rates, surface_rates, flux, &
                        average_sources, balance_residuals)
  end subroutine write_current_snapshot

  subroutine calculate_observables()
    call volume_average_concentration(state, average)
    call volume_average_rates(state, average_rates, surface_rates)
    call surface_fluxes(state, flux)
    average_sources = matmul(nu, average_rates)
    average_rhs = average_sources - (2.0_dp/cfg%radius)*flux
    call pellet_rhs(state, derivative_state)
    call volume_average_concentration(derivative_state, derivative_average)
    balance_residuals = derivative_average-average_rhs
  end subroutine calculate_observables

  subroutine validate_nonnegative(time)
    real(dp), intent(in) :: time
    real(dp) :: minimum_value
    minimum_value = minimum_concentration(state)
    if (minimum_value < -cfg%negative_error) then
      write(*,'(a,es12.4,a,es12.4)') 'ERROR: concentration ', minimum_value, ' at time ', time
      error stop 'Substantial negative concentration; no clamping applied.'
    else if (minimum_value < -cfg%negative_warning) then
      write(*,'(a,es12.4,a,es12.4)') 'WARNING: small negative concentration ', minimum_value, ' at time ', time
    end if
  end subroutine validate_nonnegative
end program hydroiso_v1
