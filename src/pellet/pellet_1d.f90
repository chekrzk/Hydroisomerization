module pellet_1d
  use kinds, only: dp
  use species, only: ns
  use kinetics, only: kinetic_parameters, species_sources, reaction_rates
  use fvm_1d, only: radial_mesh, build_uniform_cylindrical_mesh
  implicit none
  private
  public :: initialize_pellet_model, pellet_rhs, pellet_size, pellet_mesh, &
            volume_average_concentration, volume_average_rates, surface_fluxes, &
            minimum_concentration

  type(radial_mesh), save :: pellet_mesh
  type(kinetic_parameters), save :: kinetics_model
  real(dp), save :: deff(ns) = 0.0_dp
  real(dp), save :: bulk(ns) = 0.0_dp
  logical, save :: reaction_enabled = .true.
contains
  subroutine initialize_pellet_model(nr, radius, diffusivity, c_bulk, parameters, enable_reaction)
    integer, intent(in) :: nr
    real(dp), intent(in) :: radius, diffusivity(ns), c_bulk(ns)
    type(kinetic_parameters), intent(in) :: parameters
    logical, intent(in), optional :: enable_reaction
    pellet_mesh = build_uniform_cylindrical_mesh(nr, radius)
    deff = diffusivity
    bulk = c_bulk
    kinetics_model = parameters
    reaction_enabled = .true.
    if (present(enable_reaction)) reaction_enabled = enable_reaction
  end subroutine initialize_pellet_model

  pure integer function pellet_size() result(n)
    n = pellet_mesh%n * ns
  end function pellet_size

  subroutine pellet_rhs(c_flat, dc_flat)
    real(dp), intent(in) :: c_flat(:)
    real(dp), intent(out) :: dc_flat(:)
    real(dp) :: source(ns), c_cell(ns), flux_in, flux_out
    integer :: i, j, k
    dc_flat = 0.0_dp
    do j = 1, pellet_mesh%n
      k = (j-1)*ns
      c_cell = c_flat(k+1:k+ns)
      if (reaction_enabled) then
        call species_sources(c_cell, kinetics_model, source)
      else
        source = 0.0_dp
      end if
      do i = 1, ns
        if (j == 1) then
          flux_in = 0.0_dp
        else
          flux_in = deff(i) * pellet_mesh%area(j-1) * &
                    (c_flat(k-ns+i) - c_flat(k+i)) / pellet_mesh%dr
        end if
        if (j == pellet_mesh%n) then
          flux_out = deff(i) * pellet_mesh%area(j) * &
                     (bulk(i) - c_flat(k+i)) / (0.5_dp*pellet_mesh%dr)
        else
          flux_out = deff(i) * pellet_mesh%area(j) * &
                     (c_flat(k+ns+i) - c_flat(k+i)) / pellet_mesh%dr
        end if
        dc_flat(k+i) = (flux_in + flux_out) / pellet_mesh%volume(j) + source(i)
      end do
    end do
  end subroutine pellet_rhs

  subroutine volume_average_concentration(c_flat, average)
    real(dp), intent(in) :: c_flat(:)
    real(dp), intent(out) :: average(ns)
    real(dp) :: total_volume
    integer :: j, k
    average = 0.0_dp
    total_volume = sum(pellet_mesh%volume)
    do j = 1, pellet_mesh%n
      k = (j-1)*ns
      average = average + pellet_mesh%volume(j)*c_flat(k+1:k+ns)
    end do
    average = average/total_volume
  end subroutine volume_average_concentration

  subroutine volume_average_rates(c_flat, average_rates, surface_rates)
    real(dp), intent(in) :: c_flat(:)
    real(dp), intent(out) :: average_rates(ns), surface_rates(ns)
    real(dp) :: rates(ns), total_volume
    integer :: j, k
    average_rates = 0.0_dp
    total_volume = sum(pellet_mesh%volume)
    do j = 1, pellet_mesh%n
      k = (j-1)*ns
      call reaction_rates(c_flat(k+1:k+ns), kinetics_model, rates)
      average_rates = average_rates + pellet_mesh%volume(j)*rates
    end do
    average_rates = average_rates/total_volume
    call reaction_rates(bulk, kinetics_model, surface_rates)
  end subroutine volume_average_rates

  subroutine surface_fluxes(c_flat, flux)
    real(dp), intent(in) :: c_flat(:)
    real(dp), intent(out) :: flux(ns)
    integer :: i, k
    k = (pellet_mesh%n-1)*ns
    do i = 1, ns
      flux(i) = -deff(i)*(bulk(i)-c_flat(k+i))/(0.5_dp*pellet_mesh%dr)
    end do
  end subroutine surface_fluxes

  pure function minimum_concentration(c_flat) result(cmin)
    real(dp), intent(in) :: c_flat(:)
    real(dp) :: cmin
    cmin = minval(c_flat)
  end function minimum_concentration
end module pellet_1d
