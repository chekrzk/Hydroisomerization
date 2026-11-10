module species
  implicit none
  private
  public :: ns, species_names
  integer, parameter :: ns = 15
  character(len=24), parameter :: species_names(ns) = [character(len=24) :: &
    'n-C16H34', 'mono-iso-C16H34', 'multi-iso-C16H34', 'H2', &
    'C2H6-C5H12', 'C6H14-C14H30', 'n-C16H32', 'n-C16H33+', &
    'mono-iso-C16H33+', 'multi-iso-C16H33+', 'C2H4-C5H10', &
    'C6H13-C14H29+', 'mono-iso-C16H32', 'multi-iso-C16H32', &
    'C6H12-C14H28']
end module species
