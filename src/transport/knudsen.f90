module knudsen
  use kinds, only: dp
  use constants, only: gas_constant, pi
  implicit none
  private
  public :: knudsen_diffusivity, effective_knudsen_diffusivity
contains
  elemental pure function knudsen_diffusivity(pore_radius, temperature, molar_mass) result(dk)
    real(dp), intent(in) :: pore_radius, temperature, molar_mass
    real(dp) :: dk
    dk = (2.0_dp/3.0_dp) * pore_radius * sqrt(8.0_dp*gas_constant*temperature/(pi*molar_mass))
  end function knudsen_diffusivity

  elemental pure function effective_knudsen_diffusivity(pore_radius, temperature, molar_mass, porosity, tortuosity) result(deff)
    real(dp), intent(in) :: pore_radius, temperature, molar_mass, porosity, tortuosity
    real(dp) :: deff
    deff = porosity / tortuosity * knudsen_diffusivity(pore_radius, temperature, molar_mass)
  end function effective_knudsen_diffusivity
end module knudsen
