module fv_mod

  use kinds_mod
  use const_mod
  use gamil_params_mod
  use dycore_interfaces_mod
  use weno_mod

  implicit none

  private

  public recon
  public flux_llf

contains

  subroutine recon(state)

    type(state_type), intent(inout) :: state

    integer i, j, k, l

    associate (mesh => state%mesh)
    select case (recon_h_scheme)
    case ('weno5')
      !         qt
      !     ____x____
      !    |         |
      !    |         |
      ! ql x    o    x qr
      !    |         |
      !    |____x____|
      !         qb
      do l = 1, state%nvar
        do k = mesh%kds, mesh%kde
          do j = mesh%jds, mesh%jde
            do i = mesh%ids, mesh%ide
              state%ql(i,j,k,l) = weno5(-1, state%q(i-2:i+2,j,k,l))
              state%qr(i,j,k,l) = weno5( 1, state%q(i-2:i+2,j,k,l))
              state%qb(i,j,k,l) = weno5(-1, state%q(i,j-2:j+2,k,l))
              state%qt(i,j,k,l) = weno5( 1, state%q(i,j-2:j+2,k,l))
            end do
          end do
        end do
      end do
    end select
    end associate

  end subroutine recon

  pure function flux_llf(iG, J, ql, qr, flux)

    real(r8), intent(in) :: iG(3,3)
    real(r8), intent(in) :: J
    real(r8), intent(in) :: ql(:)
    real(r8), intent(in) :: qr(:)
    procedure(flux_interface), pointer :: flux
    real(r8) flux_llf(size(ql))

    real(r8) cl, cr, fl(size(ql)), fr(size(qr))
    integer i

    cl = wave_speed_x(iG, J, ql)
    cr = wave_speed_x(iG, J, qr)

    fl = flux(iG, J, ql)
    fr = flux(iG, J, qr)

    flux_llf = 0.5_r8 * (fl + fr - max(cl, cr) * (qr - ql))

  end function flux_llf

end module fv_mod
