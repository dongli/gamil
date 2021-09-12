module const_mod

  use mesh_const_mod

  implicit none

  real(r8) :: omega   = 2.0_r8 * pi / 86400.0_r8  ! s-1
  real(r8) :: g       = 9.80616_r8                ! m2 s-2
  real(r8) :: eps     = epsilon(1.0_r8)
  real(r8) :: inf     = huge(1.0_r8)

end module const_mod
