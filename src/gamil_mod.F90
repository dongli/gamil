module gamil_mod

  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use gamil_params_mod

  implicit none

  private

  public gamil_init
  public gamil_run
  public gamil_final

  type(mesh_type) mesh

contains

  subroutine gamil_init(namelist)

    character(*), intent(in) :: namelist

    call fiona_init()
    call gamil_params_init(namelist)
    call mesh%init(nx, ny, nz, hwh=2, neq=1, radius=radius)
    call mesh%write('mesh.nc')

  end subroutine gamil_init

  subroutine gamil_run()

  end subroutine gamil_run

  subroutine gamil_final()

  end subroutine gamil_final

end module gamil_mod
