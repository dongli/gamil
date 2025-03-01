module fv_mod

  use kinds_mod
  use const_mod
  use gamil_params_mod
  use state_mod
  use static_mod
  use weno_mod

  implicit none

  private

  public fv_init
  public reconstruct
  public riemann_solver

  interface
    pure real(r8) function wave_speed_interface(iG, J, q)
      import r8
      real(r8), intent(in) :: iG(3,3)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: q(:)
    end function wave_speed_interface
    pure function flux_interface(iG, J, q)
      import r8
      real(r8), intent(in) :: iG(3,3)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: q(:)
      real(r8) flux_interface(size(q))
    end function flux_interface
    pure real(r8) function reconstruct_interface(dir, q)
      import r8
      integer , intent(in) :: dir
      real(r8), intent(in) :: q(:)
    end function reconstruct_interface
    pure function riemann_solver_interface(iG, J, ql, qr, wave_speed, flux)
      import r8, wave_speed_interface, flux_interface
      real(r8), intent(in) :: iG(3,3)
      real(r8), intent(in) :: J
      real(r8), intent(in) :: ql(:)
      real(r8), intent(in) :: qr(:)
      procedure(wave_speed_interface), intent(in), pointer :: wave_speed
      procedure(flux_interface), intent(in), pointer :: flux
      real(r8) riemann_solver_interface(size(ql))
    end function riemann_solver_interface
  end interface

  procedure(reconstruct_interface   ), pointer :: reconstruct    => null()
  procedure(riemann_solver_interface), pointer :: riemann_solver => null()

contains

  subroutine fv_init(recon_h_scheme, riemann_solver_type)

    character(*), intent(in) :: recon_h_scheme
    character(*), intent(in) :: riemann_solver_type

    select case (recon_h_scheme)
    case ('weno5')
      reconstruct => weno5
    end select
    select case (riemann_solver_type)
    case ('llf')
      riemann_solver => riemann_solver_llf
    end select

  end subroutine fv_init

  pure function riemann_solver_llf(iG, J, ql, qr, wave_speed, flux)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: ql(:)
    real(r8), intent(in) :: qr(:)
    procedure(wave_speed_interface), intent(in), pointer :: wave_speed
    procedure(flux_interface), intent(in), pointer :: flux
    real(r8) riemann_solver_llf(size(ql))

    real(r8) cl, cr, fl(size(ql)), fr(size(qr))
    integer i

    cl = wave_speed(iG, J, ql)
    cr = wave_speed(iG, J, qr)

    fl = flux(iG, J, ql)
    fr = flux(iG, J, qr)

    riemann_solver_llf = 0.5_r8 * (fl + fr - max(cl, cr) * (qr - ql))

  end function riemann_solver_llf

end module fv_mod
