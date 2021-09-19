module dycore_interfaces_mod

  use kinds_mod
  use state_mod
  use static_mod

  implicit none

  interface
    pure real(r8) function wave_speed_interface(iG, J, q)
      import r8
      real(r8), intent(in) :: iG(:,:)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: q(:)
    end function wave_speed_interface
    pure function flux_interface(iG, J, q)
      import r8
      real(r8), intent(in) :: iG(3,3)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: q(3)
      real(r8) flux_interface(3)
    end function flux_interface
    pure function sources_interface(iG, J, f, CS, q, zs, dhdx, dhdy)
      import r8
      real(r8), intent(in) :: iG(3,3)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: CS(3,3,3)
      real(r8), intent(in) :: f
      real(r8), intent(in) :: q(3)
      real(r8), intent(in) :: zs
      real(r8), intent(in) :: dhdx
      real(r8), intent(in) :: dhdy
      real(r8) sources_interface(3)
    end function sources_interface
    subroutine raw_to_cons_interface(state, static)
      import state_type, static_type
      type(state_type ), intent(inout) :: state
      type(static_type), intent(in   ) :: static
    end subroutine raw_to_cons_interface
    subroutine before_recon_interface(state, static)
      import state_type, static_type
      type(state_type ), intent(inout) :: state
      type(static_type), intent(in   ) :: static
    end subroutine before_recon_interface
    subroutine after_flux_interface(state, static)
      import state_type, static_type
      type(state_type ), intent(inout) :: state
      type(static_type), intent(in   ) :: static
    end subroutine after_flux_interface
  end interface

  procedure(wave_speed_interface  ), pointer :: wave_speed_x
  procedure(wave_speed_interface  ), pointer :: wave_speed_y
  procedure(flux_interface        ), pointer :: flux_x
  procedure(flux_interface        ), pointer :: flux_y
  procedure(sources_interface     ), pointer :: sources
  procedure(raw_to_cons_interface ), pointer :: raw_to_cons
  procedure(before_recon_interface), pointer :: before_recon
  procedure(after_flux_interface  ), pointer :: after_flux

end module dycore_interfaces_mod
