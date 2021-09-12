module mesh_const_mod

  use kinds_mod

  implicit none

  real(r8), parameter :: pi     = atan(1.0d0) * 4.0d0
  real(r8), parameter :: pi2    = pi * 2.0d0
  real(r8), parameter :: pi0p5  = pi * 0.5d0
  real(r8), parameter :: pi0p25 = pi * 0.25d0
  real(r8), parameter :: rad    = pi / 180.0d0
  real(r8), parameter :: deg    = 180.0d0 / pi

  ! Domain types
  integer, parameter ::   global_domain = 1
  integer, parameter :: regional_domain = 2

  ! Orient
  integer, parameter :: left   = 1
  integer, parameter :: right  = 2
  integer, parameter :: bottom = 3
  integer, parameter :: top    = 4

end module mesh_const_mod
