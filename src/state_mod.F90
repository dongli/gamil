module state_mod

  use kinds_mod
  use latlon_array_mod

  implicit none

  type state_type
    logical :: initialized = .false.
    type(latlon_array_type), allocatable :: array
    real(r8), pointer, dimension(:,:,:) :: gz
    real(r8), pointer, dimension(:,:,:) :: u
    real(r8), pointer, dimension(:,:,:) :: v
  contains
    procedure :: init  => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh)

    class(state_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in) :: mesh

    call this%clear()

    allocate(this%array)

    call this%array%init(mesh)
    call this%array%add_var('gz', 'Geopotential'             , 'm2 s-2', loc='CA', output='T')
    call this%array%add_var('u' , 'Zonal wind component'     , 'm s-1' , loc='CA', output='T')
    call this%array%add_var('v' , 'Meridional wind component', 'm s-1' , loc='CA', output='T')
    call this%array%allocate_arrays()

    call this%array%get_array(this%gz, var_name='gz', loc='CA')
    call this%array%get_array(this%u , var_name='u' , loc='CA')
    call this%array%get_array(this%v , var_name='v' , loc='CA')

    this%initialized = .true.

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%array)) deallocate(this%array)

    this%initialized = .false.

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
