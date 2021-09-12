module test_swm_mz_mod

  use flogger
  use kinds_mod
  use const_mod
  use dycore_mod

  implicit none

  private

  public test_swm_mz_set_ic

contains

  subroutine test_swm_mz_set_ic(dycore)

    type(dycore_type), intent(inout), target :: dycore

    real(r8) lon, lat, cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    real(r8) alpha, u0, gz0, lon0, lat0, gzs0, R
    integer ids, ide, jds, jde, pc
    integer i, j, k

    alpha  = 0.0
    u0     = 20.0
    gz0    = 5960.0 * g
    lon0   = pi * 1.5
    lat0   = pi / 6.0
    gzs0   = 2000.0 * g
    R      = pi / 9.0
    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    call dycore%mesh%get_params(ids=ids, ide=ide, jds=jds, jde=jde, pc=pc)

    do j = jds, jde
      do i = ids, ide
        lon = dycore%mesh%lon(pc,i,j)
        lat = dycore%mesh%lat(pc,i,j)
        dlon = abs(lon - lon0)
        dlon = min(dlon, 2 * pi - dlon)
        d = min(R, sqrt(dlon**2 + (lat - lat0)**2))
        cos_lat = cos(lat)
        sin_lat = sin(lat)
        cos_lon = cos(lon)
        sin_lon = sin(lon)

        dycore%static%gzs(i,j) = gzs0 * (1.0 - d / R)
        dycore%state(1)%u (i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
        dycore%state(1)%v (i,j,1) = - u0 * sin_lon * sin_alpha
        dycore%state(1)%gz(i,j,1) = gz0 - (dycore%mesh%radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do

  end subroutine test_swm_mz_set_ic

end module test_swm_mz_mod
