module solver_api
  use kinds, only: dp
  use model_types, only: solver_stats
  implicit none
  private
  public :: stiff_solver

  type, abstract :: stiff_solver
  contains
    procedure(advance_interface), deferred :: advance
    procedure(statistics_interface), deferred :: statistics
    procedure(finalize_interface), deferred :: finalize
  end type stiff_solver

  abstract interface
    subroutine advance_interface(self, target_time, reached_time)
      import :: stiff_solver, dp
      class(stiff_solver), intent(inout) :: self
      real(dp), intent(in) :: target_time
      real(dp), intent(out) :: reached_time
    end subroutine advance_interface
    subroutine statistics_interface(self, stats)
      import :: stiff_solver, solver_stats
      class(stiff_solver), intent(in) :: self
      type(solver_stats), intent(out) :: stats
    end subroutine statistics_interface
    subroutine finalize_interface(self)
      import :: stiff_solver
      class(stiff_solver), intent(inout) :: self
    end subroutine finalize_interface
  end interface
end module solver_api
