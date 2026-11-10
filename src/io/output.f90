module output_io
  use kinds, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use species, only: ns, species_names
  use fvm_1d, only: radial_mesh
  use model_types, only: solver_stats
  implicit none
  private
  public :: output_files, open_outputs, write_snapshot, write_statistics, close_outputs

  type :: output_files
    integer :: profiles = -1
    integer :: averages = -1
    integer :: rates = -1
    integer :: fluxes = -1
    integer :: sources = -1
    integer :: balances = -1
  end type output_files
contains
  subroutine open_outputs(directory, files)
    character(len=*), intent(in) :: directory
    type(output_files), intent(out) :: files
    integer :: i
    open(newunit=files%profiles, file=trim(directory)//'/profiles.csv', status='replace', action='write')
    open(newunit=files%averages, file=trim(directory)//'/averages.csv', status='replace', action='write')
    open(newunit=files%rates, file=trim(directory)//'/reaction_rates.csv', status='replace', action='write')
    open(newunit=files%fluxes, file=trim(directory)//'/surface_fluxes.csv', status='replace', action='write')
    open(newunit=files%sources, file=trim(directory)//'/average_sources.csv', status='replace', action='write')
    open(newunit=files%balances, file=trim(directory)//'/balance_residuals.csv', status='replace', action='write')
    write(files%profiles,'(a)',advance='no') 'time_s,cell,r_m'
    write(files%averages,'(a)',advance='no') 'time_s'
    write(files%rates,'(a)',advance='no') 'time_s'
    write(files%fluxes,'(a)',advance='no') 'time_s'
    write(files%sources,'(a)',advance='no') 'time_s'
    write(files%balances,'(a)',advance='no') 'time_s'
    do i = 1, ns
      write(files%profiles,'(a)',advance='no') ','//trim(species_names(i))
      write(files%averages,'(a)',advance='no') ','//trim(species_names(i))
      write(files%rates,'(a,i0)',advance='no') ',R',i
      write(files%rates,'(a,i0)',advance='no') ',eta_R',i
      write(files%fluxes,'(a)',advance='no') ','//trim(species_names(i))
      write(files%sources,'(a)',advance='no') ','//trim(species_names(i))
      write(files%balances,'(a)',advance='no') ','//trim(species_names(i))
    end do
    write(files%profiles,*)
    write(files%averages,*)
    write(files%rates,*)
    write(files%fluxes,*)
    write(files%sources,*)
    write(files%balances,*)
  end subroutine open_outputs

  subroutine write_snapshot(files, time, mesh, c_flat, average, avg_rates, surface_rates, flux, sources, residuals)
    type(output_files), intent(in) :: files
    real(dp), intent(in) :: time, c_flat(:), average(ns), avg_rates(ns), surface_rates(ns), flux(ns)
    real(dp), intent(in) :: sources(ns), residuals(ns)
    type(radial_mesh), intent(in) :: mesh
    real(dp) :: eta
    integer :: i, j, k
    do j = 1, mesh%n
      k = (j-1)*ns
      write(files%profiles,'(es24.16e3,a,i0,a,es24.16e3)',advance='no') time, ',', j, ',', mesh%center(j)
      do i = 1, ns
        write(files%profiles,'(a,es24.16e3)',advance='no') ',', c_flat(k+i)
      end do
      write(files%profiles,*)
    end do
    write(files%averages,'(es24.16e3)',advance='no') time
    write(files%rates,'(es24.16e3)',advance='no') time
    write(files%fluxes,'(es24.16e3)',advance='no') time
    write(files%sources,'(es24.16e3)',advance='no') time
    write(files%balances,'(es24.16e3)',advance='no') time
    do i = 1, ns
      write(files%averages,'(a,es24.16e3)',advance='no') ',', average(i)
      if (abs(surface_rates(i)) > tiny(1.0_dp)) then
        eta = avg_rates(i)/surface_rates(i)
      else
        eta = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      write(files%rates,'(a,es24.16e3,a,es24.16e3)',advance='no') ',', avg_rates(i), ',', eta
      write(files%fluxes,'(a,es24.16e3)',advance='no') ',', flux(i)
      write(files%sources,'(a,es24.16e3)',advance='no') ',', sources(i)
      write(files%balances,'(a,es24.16e3)',advance='no') ',', residuals(i)
    end do
    write(files%averages,*)
    write(files%rates,*)
    write(files%fluxes,*)
    write(files%sources,*)
    write(files%balances,*)
  end subroutine write_snapshot

  subroutine write_statistics(directory, stats, min_concentration)
    character(len=*), intent(in) :: directory
    type(solver_stats), intent(in) :: stats
    real(dp), intent(in) :: min_concentration
    integer :: unit
    open(newunit=unit, file=trim(directory)//'/solver_statistics.csv', status='replace', action='write')
    write(unit,'(a)') 'metric,value'
    write(unit,'(a,i0)') 'steps,', stats%steps
    write(unit,'(a,i0)') 'rhs_evaluations,', stats%rhs_evaluations
    write(unit,'(a,i0)') 'error_test_failures,', stats%error_test_failures
    write(unit,'(a,i0)') 'nonlinear_iterations,', stats%nonlinear_iterations
    write(unit,'(a,i0)') 'nonlinear_failures,', stats%nonlinear_failures
    write(unit,'(a,es24.16e3)') 'minimum_concentration,', min_concentration
    close(unit)
  end subroutine write_statistics

  subroutine close_outputs(files)
    type(output_files), intent(inout) :: files
    if (files%profiles /= -1) close(files%profiles)
    if (files%averages /= -1) close(files%averages)
    if (files%rates /= -1) close(files%rates)
    if (files%fluxes /= -1) close(files%fluxes)
    if (files%sources /= -1) close(files%sources)
    if (files%balances /= -1) close(files%balances)
  end subroutine close_outputs
end module output_io
