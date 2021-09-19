module static_mod

  use kinds_mod
  use latlon_mesh_mod
  use latlon_array_mod

  implicit none

  type static_type
    logical :: initialized = .false.
    type(latlon_array_type), allocatable :: array
    real(r8), pointer, dimension(:,:) :: zs
  contains
    procedure :: init     => static_init
    procedure :: clear    => static_clear
    final :: static_final
  end type static_type

contains

  subroutine static_init(this, mesh)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh

    call this%clear()

    allocate(this%array)

    call this%array%init(mesh)
    call this%array%add_var('zs', 'Surface height', 'm', loc='CA', with_halo='T', fill_halo='T', output='T', only_2d=.true.)
    call this%array%allocate_arrays()

    call this%array%get_array(this%zs, var_name='zs', loc='CA')

    this%initialized = .true.

  end subroutine static_init

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%array)) deallocate(this%array)

    this%initialized = .false.

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod
