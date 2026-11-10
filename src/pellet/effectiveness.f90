module effectiveness
  use kinds, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: cylindrical_effectiveness, numerical_effectiveness, cylindrical_first_order_fvm_eta
contains
  elemental function cylindrical_effectiveness(phi) result(eta)
    real(dp), intent(in) :: phi
    real(dp) :: eta, x
    x = abs(phi)
    if (abs(phi) < 1.0e-6_dp) then
      eta = 1.0_dp - phi*phi/8.0_dp
    else
      eta = 2.0_dp*modified_bessel_ratio_i1_i0(x)/x
    end if
  end function cylindrical_effectiveness

  elemental pure function modified_bessel_ratio_i1_i0(x) result(ratio)
    real(dp), intent(in) :: x
    real(dp) :: ratio, sum0, sum1, term0, term1, q
    integer :: k
    if (x > 20.0_dp) then
      ! Asymptotic ratio; avoids overflow of I0 and I1 for large Thiele moduli.
      ratio = 1.0_dp - 0.5_dp/x - 0.125_dp/(x*x)
      return
    end if
    q = 0.25_dp*x*x
    term0 = 1.0_dp
    sum0 = term0
    term1 = 0.5_dp*x
    sum1 = term1
    do k = 1, 200
      term0 = term0*q/real(k*k, dp)
      term1 = term1*q/real(k*(k+1), dp)
      sum0 = sum0 + term0
      sum1 = sum1 + term1
      if (abs(term0) <= epsilon(1.0_dp)*abs(sum0) .and. &
          abs(term1) <= epsilon(1.0_dp)*abs(sum1)) exit
    end do
    ratio = sum1/sum0
  end function modified_bessel_ratio_i1_i0

  elemental function numerical_effectiveness(average_rate, surface_rate) result(eta)
    real(dp), intent(in) :: average_rate, surface_rate
    real(dp) :: eta
    if (abs(surface_rate) <= tiny(1.0_dp)) then
      eta = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      eta = average_rate/surface_rate
    end if
  end function numerical_effectiveness

  function cylindrical_first_order_fvm_eta(phi, n) result(eta)
    use constants, only: pi
    real(dp), intent(in) :: phi
    integer, intent(in) :: n
    real(dp) :: eta
    real(dp), allocatable :: lower(:), diagonal(:), upper(:), rhs(:), c(:), volume(:)
    real(dp) :: dr, r_left, r_right, area_left, area_right, g_left, g_right, factor
    integer :: j
    allocate(lower(n), diagonal(n), upper(n), rhs(n), c(n), volume(n))
    dr = 1.0_dp/real(n,dp)
    lower = 0.0_dp; upper = 0.0_dp; rhs = 0.0_dp
    do j = 1, n
      r_left = real(j-1,dp)*dr
      r_right = real(j,dp)*dr
      area_left = 2.0_dp*pi*r_left
      area_right = 2.0_dp*pi*r_right
      volume(j) = pi*(r_right*r_right-r_left*r_left)
      g_left = merge(area_left/dr, 0.0_dp, j > 1)
      if (j < n) then
        g_right = area_right/dr
      else
        g_right = area_right/(0.5_dp*dr)
        rhs(j) = -g_right
      end if
      diagonal(j) = -(g_left+g_right+phi*phi*volume(j))
      if (j > 1) lower(j) = g_left
      if (j < n) upper(j) = g_right
    end do
    do j = 2, n
      factor = lower(j)/diagonal(j-1)
      diagonal(j) = diagonal(j)-factor*upper(j-1)
      rhs(j) = rhs(j)-factor*rhs(j-1)
    end do
    c(n) = rhs(n)/diagonal(n)
    do j = n-1, 1, -1
      c(j) = (rhs(j)-upper(j)*c(j+1))/diagonal(j)
    end do
    eta = sum(c*volume)/sum(volume)
  end function cylindrical_first_order_fvm_eta
end module effectiveness
