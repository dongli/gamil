program gamil_driver

  use gamil_mod
  use test_swm_mz_mod
  use test_swm_rh_mod

  implicit none

  character(256) namelist

  call get_command_argument(1, namelist)

  call gamil_init(namelist)

  call test_swm_mz_set_ic(dycore)

  call gamil_run()

  call gamil_final()

end program gamil_driver
