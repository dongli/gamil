module swm_mod

  use const_mod
  use kinds_mod
  use state_mod
  use static_mod
  use latlon_parallel_mod

  implicit none

  private

  public swm_wave_speed_x
  public swm_wave_speed_y
  public swm_flux_x
  public swm_flux_y
  public swm_sources
  public swm_raw_to_conservative

contains

  pure real(r8) function swm_wave_speed_x(iG, J, q) result(res)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: q(:)

    real(r8) u, h

    u = q(2) / q(1)  ! NOTE: u is contravariant wind component.
    h = q(1) / J ! NOTE: q(1) contains surface height which is added in before_recon.

    res = abs(u) + sqrt(iG(1,1) * g * h)

  end function swm_wave_speed_x

  pure real(r8) function swm_wave_speed_y(iG, J, q) result(res)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: q(:)

    real(r8) v, h

    v = q(3) / q(1)  ! NOTE: v is contravariant wind component.
    h = q(1) / J ! NOTE: q(1) contains surface height which is added in before_recon.

    res = abs(v) + sqrt(iG(2,2) * g * h)

  end function swm_wave_speed_y

  pure function swm_flux_x(iG, J, q)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: q(:)
    real(r8) swm_flux_x(size(q))

    real(r8) u, h

    u = q(2) / q(1) ! NOTE: u is contravariant wind component.
    h = q(1) / J    ! NOTE: q(1) contains surface height which is added in before_recon.

    swm_flux_x(1) = q(2)
    swm_flux_x(2) = q(2) * u + J * iG(1,1) * g * h**2 / 2.0_r8
    swm_flux_x(3) = q(3) * u

  end function swm_flux_x

  pure function swm_flux_y(iG, J, q)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: q(:)
    real(r8) swm_flux_y(size(q))

    real(r8) v, h

    v = q(3) / q(1) ! NOTE: v is contravariant wind component.
    h = q(1) / J    ! NOTE: q(1) contains surface height which is added in before_recon.

    swm_flux_y(1) = q(3)
    swm_flux_y(2) = q(2) * v
    swm_flux_y(3) = q(3) * v + J * iG(2,2) * g * h**2 / 2.0_r8

  end function swm_flux_y

  pure function swm_sources(iG, J, CS, f, q, zs, dhdx, dhdy)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: CS(3,3,3)
    real(r8), intent(in) :: f
    real(r8), intent(in) :: q(:)
    real(r8), intent(in) :: zs
    real(r8), intent(in) :: dhdx
    real(r8), intent(in) :: dhdy
    real(r8) swm_sources(3)

    real(r8) h, u, v, J_iG11, J_iG22, T11, T12, T22

    h = q(1) / J + zs
    u = q(2) / q(1) ! NOTE: u is contravariant wind component.
    v = q(3) / q(1) ! NOTE: v is contravariant wind component.

    J_iG11 = J * iG(1,1)
    J_iG22 = J * iG(2,2)

    T11 = q(2) * u + J_iG11 * g * h**2 / 2.0_r8
    T12 = q(2) * v
    T22 = q(3) * v + J_iG22 * g * h**2 / 2.0_r8

    ! Coriolis source terms
    swm_sources(2) = &
      - CS(1,1,1) * T11 - 2 * CS(1,2,1) * T12 & ! Metric source term
      + J_iG11 * q(1) * f * v                 & ! Coriolis source term
      + J_iG11 * g * zs * dhdx                  ! Topographic source term
    swm_sources(3) = &
      - CS(1,1,2) * T11 - 2 * CS(1,2,2) * T12 & ! Metric source term
      - J_iG22 * q(1) * f * u                 & ! Coriolis source term
      + J_iG22 * g * zs * dhdy                  ! Topographic source term

  end function swm_sources

  subroutine swm_raw_to_conservative(state, static)

    ! Only called at run beginning.

    type(state_type ), intent(inout) :: state
    type(static_type), intent(in   ) :: static
    
    integer i, j, k

    k = 1
    associate (mesh => state%mesh)
    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        state%q(i,j,k,1) = mesh%J(1,i,j) * (state%h(i,j,k) - static%zs(i,j))
        state%q(i,j,k,2) = state%q(i,j,k,1) * state%uc(i,j,k)
        state%q(i,k,k,3) = state%q(i,j,k,1) * state%vc(i,j,k)
      end do
    end do
    end associate

    call fill_halo(state%array)

  end subroutine swm_raw_to_conservative

end module swm_mod
