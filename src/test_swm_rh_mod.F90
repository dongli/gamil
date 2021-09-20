module test_swm_rh_mod

  use flogger
  use kinds_mod
  use const_mod
  use gamil_params_mod
  use dycore_mod
  use swm_mod
  use latlon_process_mod
  use latlon_parallel_mod

  implicit none

  private

  public test_swm_rh_set_ic

  real(r8), parameter :: R   = 4.0d0
  real(r8), parameter :: omg = 7.848d-6
  real(r8), parameter :: z0  = 8.0d3

contains

  subroutine test_swm_rh_set_ic(dycore)

    type(dycore_type), intent(inout), target :: dycore

    real(r8) lon, lat, cos_lat, sin_lat
    real(r8) a, b, c
    integer i, j

    associate (mesh => dycore%mesh, state => dycore%state(1), static => dycore%static)
    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        lon = mesh%lon(1,i,j)
        lat = mesh%lat(1,i,j)
        cos_lat = cos(lat)
        sin_lat = sin(lat)

        static%zs(i,j) = 0

        a = cos_lat
        b = R * cos_lat**(R - 1) * sin_lat**2 * cos(R * lon)
        c = cos_lat**(R + 1) * cos(R * lon)
        state%u(i,j,1) = radius * omg * (a + b - c)

        a = R * cos_lat**(R - 1) * sin_lat * sin(R * lon)
        state%v(i,j,1) = - radius * omg * a

        a = 0.5 * omg * (2 * omega + omg) * cos_lat**2 + &
          0.25 * omg**2 * ((R + 1) * cos_lat**(2 * R + 2) + (2 * R**2 - R - 2) * cos_lat**(2 * R) - 2 * R**2 * cos_lat**(2 * R - 2))
        b = 2 * (omega + omg) * omg * cos_lat**R * &
          (R**2 + 2 * R + 2 - (R + 1)**2 * cos_lat**2) / (R + 1) / (R + 2)
        c = 0.25 * omg**2 * cos_lat**(2 * R) * ((R + 1) * cos_lat**2 - R - 2)
        state%h(i,j,1) = z0 + radius**2 * (a + b * cos(R * lon) + c * cos(2 * R * lon))
      end do
    end do
    call fill_halo(static%array)
    call dycore%calc_contravar_wind(1)
    call swm_raw_to_conservative(state, static)
    end associate

 
  end subroutine test_swm_rh_set_ic

end module test_swm_rh_mod
