module gamil_mod

  use latlon_process_mod
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use gamil_params_mod
  use dycore_mod
  use history_mod

  implicit none

  private

  public dycore

  public gamil_init
  public gamil_run
  public gamil_final

  type(mesh_type), allocatable :: mesh
  type(dycore_type), allocatable :: dycore

contains

  subroutine gamil_init(namelist)

    character(*), intent(in) :: namelist

    allocate(proc, mesh, dycore)

    call fiona_init()
    call gamil_params_init(namelist)
    call proc%init()
    call proc%decomp_domain(mesh, nx, ny, nz, hwx=2, hwy=2, neq=1, r=radius)
    call dycore%init(mesh, model_type)
    call history_init(case_name, dycore)

  end subroutine gamil_init

  subroutine gamil_run()

    call history_write(dycore, 1)

  end subroutine gamil_run

  subroutine gamil_final()

    deallocate(proc, mesh, dycore)

  end subroutine gamil_final

end module gamil_mod
