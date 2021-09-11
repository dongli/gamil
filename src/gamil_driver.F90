program gamil_driver

  use gamil_mod

  implicit none

  character(256) namelist

  call get_command_argument(1, namelist)

  call gamil_init(namelist)

  call gamil_final()

end program gamil_driver
