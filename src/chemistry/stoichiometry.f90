module stoichiometry
  use kinds, only: dp
  use species, only: ns
  implicit none
  private
  public :: build_stoichiometric_matrix
contains
  pure function build_stoichiometric_matrix() result(nu)
    real(dp) :: nu(ns, ns)
    nu = 0.0_dp

    nu(1,1) = -1.0_dp; nu(7,1) = 1.0_dp; nu(4,1) = 1.0_dp
    nu(7,2) = -1.0_dp; nu(8,2) = 1.0_dp
    nu(8,3) = -1.0_dp; nu(9,3) = 1.0_dp
    nu(8,4) = -1.0_dp; nu(10,4) = 1.0_dp
    nu(9,5) = -1.0_dp; nu(10,5) = 1.0_dp
    nu(8,6) = -1.0_dp; nu(11,6) = 1.0_dp; nu(12,6) = 1.0_dp
    nu(9,7) = -1.0_dp; nu(11,7) = 1.0_dp; nu(12,7) = 1.0_dp
    nu(10,8) = -1.0_dp; nu(11,8) = 1.0_dp; nu(12,8) = 1.0_dp
    nu(9,9) = -1.0_dp; nu(13,9) = 1.0_dp
    nu(10,10) = -1.0_dp; nu(14,10) = 1.0_dp
    nu(12,11) = -1.0_dp; nu(15,11) = 1.0_dp
    nu(11,12) = -1.0_dp; nu(4,12) = -1.0_dp; nu(5,12) = 1.0_dp
    nu(13,13) = -1.0_dp; nu(4,13) = -1.0_dp; nu(2,13) = 1.0_dp
    nu(14,14) = -1.0_dp; nu(4,14) = -1.0_dp; nu(3,14) = 1.0_dp
    nu(15,15) = -1.0_dp; nu(4,15) = -1.0_dp; nu(6,15) = 1.0_dp
  end function build_stoichiometric_matrix
end module stoichiometry
