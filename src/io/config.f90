module config_io
  use kinds, only: dp, i8
  use species, only: ns
  use model_types, only: simulation_config
  implicit none
  private
  public :: read_simulation_config
contains
  subroutine read_simulation_config(path, cfg)
    character(len=*), intent(in) :: path
    type(simulation_config), intent(out) :: cfg
    character(len=32) :: run_mode
    character(len=64) :: run_name
    character(len=16) :: simulation_mode
    integer :: nr, unit, ios
    real(dp) :: temperature, pressure, radius, porosity, tortuosity, pore_radius
    real(dp) :: acid_activity, concentration_reference
    real(dp) :: t_start, t_end, output_interval, rtol, atol, initial_step, max_step
    real(dp) :: negative_warning, negative_error
    real(dp) :: steady_tolerance, steady_maximum_time
    integer(i8) :: max_steps
    real(dp) :: molar_mass(ns), c_bulk(ns), c_initial(ns)
    real(dp) :: ln_k0_forward(ns), activation_forward(ns)
    real(dp) :: ln_k0_reverse(ns), activation_reverse(ns)
    logical :: reversible(ns)
    namelist /simulation/ run_mode, run_name, simulation_mode, nr, temperature, pressure, radius, &
      porosity, tortuosity, pore_radius, acid_activity, concentration_reference, &
      t_start, t_end, output_interval, rtol, atol, initial_step, max_step, max_steps, &
      negative_warning, negative_error, steady_tolerance, steady_maximum_time, molar_mass, c_bulk, c_initial, &
      ln_k0_forward, activation_forward, ln_k0_reverse, activation_reverse, reversible

    run_mode = cfg%run_mode; run_name = cfg%run_name; simulation_mode = cfg%simulation_mode; nr = cfg%nr
    temperature = cfg%temperature; pressure = cfg%pressure; radius = cfg%radius
    porosity = cfg%porosity; tortuosity = cfg%tortuosity; pore_radius = cfg%pore_radius
    acid_activity = cfg%acid_activity; concentration_reference = cfg%concentration_reference
    t_start = cfg%t_start; t_end = cfg%t_end; output_interval = cfg%output_interval
    rtol = cfg%rtol; atol = cfg%atol; initial_step = cfg%initial_step; max_step = cfg%max_step
    max_steps = cfg%max_steps; negative_warning = cfg%negative_warning
    negative_error = cfg%negative_error; molar_mass = cfg%molar_mass
    steady_tolerance = cfg%steady_tolerance; steady_maximum_time = cfg%steady_maximum_time
    c_bulk = cfg%c_bulk; c_initial = cfg%c_initial
    ln_k0_forward = cfg%ln_k0_forward; activation_forward = cfg%activation_forward
    ln_k0_reverse = cfg%ln_k0_reverse; activation_reverse = cfg%activation_reverse
    reversible = cfg%reversible

    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open generated Fortran configuration.'
    read(unit, nml=simulation, iostat=ios)
    close(unit)
    if (ios /= 0) then
      write(*,'(a,i0)') 'Failed to parse Fortran namelist, iostat=', ios
      error stop 'Invalid generated Fortran configuration.'
    end if

    cfg%run_mode = run_mode; cfg%run_name = run_name; cfg%simulation_mode = simulation_mode; cfg%nr = nr
    cfg%temperature = temperature; cfg%pressure = pressure; cfg%radius = radius
    cfg%porosity = porosity; cfg%tortuosity = tortuosity; cfg%pore_radius = pore_radius
    cfg%acid_activity = acid_activity; cfg%concentration_reference = concentration_reference
    cfg%t_start = t_start; cfg%t_end = t_end; cfg%output_interval = output_interval
    cfg%rtol = rtol; cfg%atol = atol; cfg%initial_step = initial_step; cfg%max_step = max_step
    cfg%max_steps = max_steps; cfg%negative_warning = negative_warning
    cfg%negative_error = negative_error; cfg%molar_mass = molar_mass
    cfg%steady_tolerance = steady_tolerance; cfg%steady_maximum_time = steady_maximum_time
    cfg%c_bulk = c_bulk; cfg%c_initial = c_initial
    cfg%ln_k0_forward = ln_k0_forward; cfg%activation_forward = activation_forward
    cfg%ln_k0_reverse = ln_k0_reverse; cfg%activation_reverse = activation_reverse
    cfg%reversible = reversible
  end subroutine read_simulation_config
end module config_io
