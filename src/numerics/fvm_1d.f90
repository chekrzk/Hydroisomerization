module fvm_1d
  use kinds, only: dp
  use constants, only: pi
  implicit none
  private
  public :: radial_mesh, build_uniform_cylindrical_mesh

  type :: radial_mesh
    integer :: n = 0
    real(dp) :: radius = 0.0_dp
    real(dp) :: dr = 0.0_dp
    real(dp), allocatable :: face(:)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: volume(:)
    real(dp), allocatable :: area(:)
  end type radial_mesh
contains
  function build_uniform_cylindrical_mesh(n, radius) result(mesh)
    integer, intent(in) :: n
    real(dp), intent(in) :: radius
    type(radial_mesh) :: mesh
    integer :: j
    mesh%n = n
    mesh%radius = radius
    mesh%dr = radius / real(n, dp)
    allocate(mesh%face(0:n), mesh%center(n), mesh%volume(n), mesh%area(0:n))
    do j = 0, n
      mesh%face(j) = real(j, dp) * mesh%dr
      mesh%area(j) = 2.0_dp*pi*mesh%face(j)
    end do
    do j = 1, n
      mesh%center(j) = 0.5_dp*(mesh%face(j-1) + mesh%face(j))
      mesh%volume(j) = pi*(mesh%face(j)**2 - mesh%face(j-1)**2)
    end do
  end function build_uniform_cylindrical_mesh
end module fvm_1d
