module weno_mod

  use const_mod

  implicit none

  private

  public weno3
  public weno5

contains

  pure real(r8) function weno3(dir, f) result(res)

    integer , intent(in) :: dir   ! Which side to reconstruct (-1: i-1/2, i: i+1/2)
    real(r8), intent(in) :: f(3)

    real(r8), parameter :: c11 = -0.5_r8, c12 =  1.5_r8
    real(r8), parameter :: c21 =  0.5_r8, c22 =  0.5_r8
    real(r8), parameter :: linear_coef(2) = [1.0_r8 / 3.0_r8, 2.0_r8 / 3.0_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    integer im1, i, ip1
    real(r8) fs(2), beta(2), alpha(2)

    select case (dir)
    case (-1) ! i-1/2
      im1 = 3; i = 2; ip1 = 1
    case ( 1) ! i+1/2
      im1 = 1; i = 2; ip1 = 3
    end select

    ! Calculate reconstructed value at each small stencil.
    fs(1) = c11 * f(im1) + c12 * f(i  )
    fs(2) = c21 * f(i  ) + c22 * f(ip1)

    ! Calculate smoothness indicators beta.
    beta(1) = (f(i  ) - f(im1))**2
    beta(2) = (f(ip1) - f(i  ))**2

    ! Calculate nonlinear weights.
    alpha = linear_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno3

  pure real(r8) function weno5(dir, f) result(res)

    integer , intent(in) :: dir   ! Which side to reconstruct (-1: i-1/2, i: i+1/2)
    real(r8), intent(in) :: f(5)

    real(r8), parameter :: c11 =  1.0_r8 / 3.0_r8, c12 = -7.0_r8 / 6.0_r8, c13 = 11.0_r8 / 6.0_r8
    real(r8), parameter :: c21 = -1.0_r8 / 6.0_r8, c22 =  5.0_r8 / 6.0_r8, c23 =  1.0_r8 / 3.0_r8
    real(r8), parameter :: c31 =  1.0_r8 / 3.0_r8, c32 =  5.0_r8 / 6.0_r8, c33 = -1.0_r8 / 6.0_r8
    real(r8), parameter :: b1 = 13.0_r8 / 12.0_r8, b2 = 0.25_r8
    real(r8), parameter :: linear_coef(3) = [1.0_r8 / 10.0_r8, 3.0_r8 / 5.0_r8, 3.0_r8 / 10.0_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    integer im2, im1, i, ip1, ip2
    real(r8) fs(3), beta(3), alpha(3)

    select case (dir)
    case (-1) ! i-1/2
      im2 = 5; im1 = 4; i = 3; ip1 = 2; ip2 = 1
    case ( 1) ! i+1/2
      im2 = 1; im1 = 2; i = 3; ip1 = 4; ip2 = 5
    end select

    ! Calculate reconstructed value at each small stencil.
    fs(1) = c11 * f(im2) + c12 * f(im1) + c13 * f(i  )
    fs(2) = c21 * f(im1) + c22 * f(i  ) + c23 * f(ip1)
    fs(3) = c31 * f(i  ) + c32 * f(ip1) + c33 * f(ip2)

    ! Calculate smoothness indicators beta.
    beta(1) = b1 * (f(im2) - 2 * f(im1) + f(i  ))**2 + b2 * (    f(im2) - 4 * f(im1) + 3 * f(i  ))**2
    beta(2) = b1 * (f(im1) - 2 * f(i  ) + f(ip1))**2 + b2 * (    f(im1)          -         f(ip1))**2
    beta(3) = b1 * (f(i  ) - 2 * f(ip1) + f(ip2))**2 + b2 * (3 * f(i  ) - 4 * f(ip1) +     f(ip2))**2

    ! Calculate nonlinear weights.
    alpha = linear_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno5

end module weno_mod
