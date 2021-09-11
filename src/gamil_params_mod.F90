module gamil_params_mod

  use kinds_mod

  implicit none

  integer :: nx = 0
  integer :: ny = 0
  integer :: nz = 1

  real(r8) :: radius = 6.37122e6_r8

  namelist /gamil_control/ &
    nx                   , &
    ny                   , &
    nz

contains

  subroutine gamil_params_init(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gamil_control)
    close(10)

    write(*, *) 'nx = ', nx
    write(*, *) 'ny = ', ny
    write(*, *) 'nz = ', nz

  end subroutine gamil_params_init

end module gamil_params_mod
