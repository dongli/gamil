module gamil_params_mod

  use kinds_mod

  implicit none

  character(30) :: case_name      = ''
  character(10) :: model_type     = 'swm'
  character(30) :: recon_h_scheme = 'weno5'
  character(30) :: recon_v_scheme = 'weno5'

  integer :: nx = 0
  integer :: ny = 0
  integer :: nz = 1

  real(r8) :: reduce_factors(100)

  real(r8) :: radius = 6.37122e6_r8

  namelist /gamil_control/ &
    case_name            , &
    model_type           , &
    nx, ny, nz           , &
    reduce_factors

contains

  subroutine gamil_params_init(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gamil_control)
    close(10)

  end subroutine gamil_params_init

end module gamil_params_mod
