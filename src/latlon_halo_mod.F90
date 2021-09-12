module latlon_halo_mod

  use mpi
  use latlon_mesh_mod

  implicit none

  type latlon_halo_type
    logical :: initialized = .true.
  contains
    procedure :: init    => latlon_halo_init
    procedure :: clear   => latlon_halo_clear
    final :: latlon_halo_final
  end type latlon_halo_type

contains

  subroutine latlon_halo_init(this, mesh, orient, dtype, host_id, proc_id)

    class(latlon_halo_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: mesh
    integer, intent(in) :: orient
    integer, intent(in) :: dtype
    integer, intent(in) :: host_id
    integer, intent(in) :: proc_id

    call this%clear()

  end subroutine latlon_halo_init

  subroutine latlon_halo_clear(this)

    class(latlon_halo_type), intent(inout) :: this

  end subroutine latlon_halo_clear

  subroutine latlon_halo_final(this)

    type(latlon_halo_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_halo_final

end module latlon_halo_mod
