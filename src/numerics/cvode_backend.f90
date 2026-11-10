module cvode_backend
  use, intrinsic :: iso_c_binding
  use kinds, only: dp, i8
  use model_types, only: simulation_config, solver_stats
  use solver_api, only: stiff_solver
  use pellet_1d, only: pellet_rhs
  use fsundials_core_mod
  use fcvode_mod
  use fnvector_serial_mod
  use fsunmatrix_dense_mod
  use fsunlinsol_dense_mod
  implicit none
  private
  public :: cvode_solver, create_cvode_solver

  type, extends(stiff_solver) :: cvode_solver
    type(c_ptr) :: context = c_null_ptr
    type(c_ptr) :: memory = c_null_ptr
    type(N_Vector), pointer :: vector => null()
    type(SUNMatrix), pointer :: matrix => null()
    type(SUNLinearSolver), pointer :: linear_solver => null()
  contains
    procedure :: advance => cvode_advance
    procedure :: statistics => cvode_statistics
    procedure :: finalize => cvode_finalize
  end type cvode_solver
contains
  function create_cvode_solver(state, cfg) result(solver)
    real(c_double), target, intent(inout) :: state(:)
    type(simulation_config), intent(in) :: cfg
    type(cvode_solver) :: solver
    integer(c_int) :: ierr
    integer(c_int64_t) :: neq

    neq = int(size(state), c_int64_t)
    ierr = FSUNContext_Create(SUN_COMM_NULL, solver%context)
    call check_flag(ierr, 'FSUNContext_Create')
    solver%vector => FN_VMake_Serial(neq, state, solver%context)
    if (.not. associated(solver%vector)) error stop 'FN_VMake_Serial failed.'
    solver%matrix => FSUNDenseMatrix(neq, neq, solver%context)
    if (.not. associated(solver%matrix)) error stop 'FSUNDenseMatrix failed.'
    solver%linear_solver => FSUNLinSol_Dense(solver%vector, solver%matrix, solver%context)
    if (.not. associated(solver%linear_solver)) error stop 'FSUNLinSol_Dense failed.'
    solver%memory = FCVodeCreate(CV_BDF, solver%context)
    if (.not. c_associated(solver%memory)) error stop 'FCVodeCreate failed.'
    ierr = FCVodeInit(solver%memory, c_funloc(rhs_callback), cfg%t_start, solver%vector)
    call check_flag(ierr, 'FCVodeInit')
    ierr = FCVodeSStolerances(solver%memory, cfg%rtol, cfg%atol)
    call check_flag(ierr, 'FCVodeSStolerances')
    ierr = FCVodeSetLinearSolver(solver%memory, solver%linear_solver, solver%matrix)
    call check_flag(ierr, 'FCVodeSetLinearSolver')
    ierr = FCVodeSetMaxNumSteps(solver%memory, int(cfg%max_steps, c_long))
    call check_flag(ierr, 'FCVodeSetMaxNumSteps')
    ierr = FCVodeSetInitStep(solver%memory, cfg%initial_step)
    call check_flag(ierr, 'FCVodeSetInitStep')
    ierr = FCVodeSetMaxStep(solver%memory, cfg%max_step)
    call check_flag(ierr, 'FCVodeSetMaxStep')
  end function create_cvode_solver

  subroutine cvode_advance(self, target_time, reached_time)
    class(cvode_solver), intent(inout) :: self
    real(dp), intent(in) :: target_time
    real(dp), intent(out) :: reached_time
    real(c_double) :: tret(1)
    integer(c_int) :: ierr
    ierr = FCVode(self%memory, target_time, self%vector, tret, CV_NORMAL)
    call check_flag(ierr, 'FCVode')
    reached_time = tret(1)
  end subroutine cvode_advance

  subroutine cvode_statistics(self, stats)
    class(cvode_solver), intent(in) :: self
    type(solver_stats), intent(out) :: stats
    integer(c_long) :: value(1)
    integer(c_int) :: ierr
    stats = solver_stats()
    ierr = FCVodeGetNumSteps(self%memory, value)
    if (ierr == 0) stats%steps = int(value(1), i8)
    ierr = FCVodeGetNumRhsEvals(self%memory, value)
    if (ierr == 0) stats%rhs_evaluations = int(value(1), i8)
    ierr = FCVodeGetNumErrTestFails(self%memory, value)
    if (ierr == 0) stats%error_test_failures = int(value(1), i8)
    ierr = FCVodeGetNumNonlinSolvIters(self%memory, value)
    if (ierr == 0) stats%nonlinear_iterations = int(value(1), i8)
    ierr = FCVodeGetNumNonlinSolvConvFails(self%memory, value)
    if (ierr == 0) stats%nonlinear_failures = int(value(1), i8)
  end subroutine cvode_statistics

  subroutine cvode_finalize(self)
    class(cvode_solver), intent(inout) :: self
    integer(c_int) :: ierr
    if (c_associated(self%memory)) call FCVodeFree(self%memory)
    if (associated(self%linear_solver)) ierr = FSUNLinSolFree(self%linear_solver)
    if (associated(self%matrix)) call FSUNMatDestroy(self%matrix)
    if (associated(self%vector)) call FN_VDestroy(self%vector)
    if (c_associated(self%context)) ierr = FSUNContext_Free(self%context)
    nullify(self%linear_solver, self%matrix, self%vector)
  end subroutine cvode_finalize

  integer(c_int) function rhs_callback(tn, sunvec_y, sunvec_f, user_data) result(ierr) bind(C)
    real(c_double), value :: tn
    type(N_Vector) :: sunvec_y, sunvec_f
    type(c_ptr), value :: user_data
    real(c_double), pointer :: y(:), f(:)
    y => FN_VGetArrayPointer(sunvec_y)
    f => FN_VGetArrayPointer(sunvec_f)
    call pellet_rhs(y, f)
    ierr = 0_c_int
  end function rhs_callback

  subroutine check_flag(flag, name)
    integer(c_int), intent(in) :: flag
    character(len=*), intent(in) :: name
    if (flag /= 0_c_int) then
      write(*,'(a,a,a,i0)') 'SUNDIALS call ', trim(name), ' failed with flag ', flag
      error stop 'CVODE failure.'
    end if
  end subroutine check_flag
end module cvode_backend
