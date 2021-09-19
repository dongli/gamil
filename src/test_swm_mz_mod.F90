module test_swm_mz_mod

  use flogger
  use kinds_mod
  use const_mod
  use dycore_mod
  use latlon_process_mod
  use latlon_parallel_mod

  implicit none

  private

  public test_swm_mz_set_ic

contains

  subroutine test_swm_mz_set_ic(dycore)

    type(dycore_type), intent(inout), target :: dycore

    real(r8) lon, lat, cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    real(r8) alpha, u0, z0, lon0, lat0, zs0, R
    integer i, j, k

    alpha  = 0.0
    u0     = 20.0
    z0     = 5960.0
    lon0   = pi * 1.5
    lat0   = pi / 6.0
    zs0    = 2000.0
    R      = pi / 9.0
    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    associate (mesh => dycore%mesh, state => dycore%state(1), static => dycore%static)
    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        lon = mesh%lon(1,i,j)
        lat = mesh%lat(1,i,j)
        dlon = abs(lon - lon0)
        dlon = min(dlon, 2 * pi - dlon)
        d = min(R, sqrt(dlon**2 + (lat - lat0)**2))
        cos_lat = cos(lat)
        sin_lat = sin(lat)
        cos_lon = cos(lon)
        sin_lon = sin(lon)

        static%zs(i,j)   = zs0 * (1.0 - d / R)
        state%u(i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
        state%v(i,j,1) = - u0 * sin_lon * sin_alpha
        state%h(i,j,1) = z0 - (mesh%r * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do
    call fill_halo(static%array)
    call dycore%calc_contra_wind(1)
    call raw_to_cons(state, static)
    end associate

  end subroutine test_swm_mz_set_ic

end module test_swm_mz_mod
