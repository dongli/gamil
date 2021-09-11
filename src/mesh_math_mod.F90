module mesh_math_mod

  use mesh_const_mod
  use gauss_quad_mod

  implicit none

  private

  public cart_to_latlon
  public latlon_to_cart
  public rotate
  public rotate_back
  public gaussian_legendre

contains

  subroutine cart_to_latlon(x, y, z, lon, lat)

    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: z
    real(8), intent(out) :: lon
    real(8), intent(out) :: lat

    lon = atan2(y, x)
    lat = asin(z) ! Assume x**2 + y**2 + z**2 = 1.
    if (lon < 0.0) lon = lon + pi2
    if (lon > pi2) lon = lon - pi2

  end subroutine cart_to_latlon

  subroutine latlon_to_cart(lon, lat, x, y, z)

    real(8), intent(in) :: lon
    real(8), intent(in) :: lat
    real(8), intent(out) :: x
    real(8), intent(out) :: y
    real(8), intent(out) :: z

    real(8) cos_lat

    cos_lat = cos(lat)
    x = cos_lat * cos(lon)
    y = cos_lat * sin(lon)
    z = sin(lat)

  end subroutine latlon_to_cart

  subroutine rotate(lon_np, lat_np, lon, lat, lon_rot, lat_rot)

    real(8), intent(in) :: lon_np
    real(8), intent(in) :: lat_np
    real(8), intent(in) :: lon
    real(8), intent(in) :: lat
    real(8), intent(out) :: lon_rot
    real(8), intent(out) :: lat_rot

    real(8) cos_lat_np, sin_lat_np, cos_lat, sin_lat, cos_dlon, sin_dlon
    real(8) tmp1, tmp2, tmp3

    cos_lat_np = cos(lat_np)
    sin_lat_np = sin(lat_np)
    cos_lat    = cos(lat   )
    sin_lat    = sin(lat   )
    cos_dlon   = cos(lon - lon_np)
    sin_dlon   = sin(lon - lon_np)

    tmp1 = cos_lat * sin_dlon
    tmp2 = cos_lat * sin_lat_np * cos_dlon - cos_lat_np * sin_lat
    lon_rot = atan2(tmp1, tmp2)
    if (lon_rot < 0.0) lon_rot = pi2 + lon_rot

    tmp1 = sin_lat * sin_lat_np
    tmp2 = cos_lat * cos_lat_np * cos_dlon
    tmp3 = tmp1 + tmp2
    tmp3 = min(1.0, max(-1.0, tmp3))
    lat_rot = asin(tmp3)

  end subroutine rotate

  subroutine rotate_back(lon_np, lat_np, lon_rot, lat_rot, lon, lat)

    real(8), intent(in) :: lon_np
    real(8), intent(in) :: lat_np
    real(8), intent(in) :: lon_rot
    real(8), intent(in) :: lat_rot
    real(8), intent(out) :: lon
    real(8), intent(out) :: lat

    real(8) cos_lat_np, sin_lat_np, cos_lon_rot, sin_lon_rot, cos_lat_rot, sin_lat_rot
    real(8) tmp1, tmp2, tmp3

    cos_lat_np = cos(lat_np)
    sin_lat_np = sin(lat_np)
    cos_lon_rot = cos(lon_rot)
    sin_lon_rot = sin(lon_rot)
    cos_lat_rot = cos(lat_rot)
    sin_lat_rot = sin(lat_rot)

    tmp1 = cos_lat_rot * sin_lon_rot
    tmp2 = sin_lat_rot * cos_lat_np + cos_lat_rot * cos_lon_rot * sin_lat_np
    lon = lon_np + atan2(tmp1, tmp2)
    if (lon > pi2) lon = lon - pi2
    if (lon < 0.0) lon = lon + pi2

    tmp1 = sin_lat_rot * sin_lat_np
    tmp2 = cos_lat_rot * cos_lat_np * cos_lon_rot
    tmp3 = tmp1 - tmp2
    tmp3 = min(1.0, max(-1.0, tmp3))
    lat = asin(tmp3)

  end subroutine rotate_back

end module mesh_math_mod