module dycore_mod

  use latlon_mesh_mod
  use static_mod
  use state_mod

  implicit none

  private

  public dycore_type

  type dycore_type
    logical :: initialized = .false.
    type(latlon_mesh_type), pointer :: mesh => null()
    type(static_type), allocatable :: static
    type(state_type), allocatable :: state(:)
  contains
    procedure :: init   => dycore_init
    procedure :: clear  => dycore_clear
    final :: dycore_final
  end type dycore_type

contains

  subroutine dycore_init(this, mesh)

    class(dycore_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh

    integer i

    call this%clear()

    this%mesh => mesh

    allocate(this%static)
    call this%static%init(mesh)

    allocate(this%state(1:1))
    do i = 1, size(this%state)
      call this%state(i)%init(mesh)
    end do

    this%initialized = .true.

  end subroutine dycore_init

  subroutine dycore_clear(this)

    class(dycore_type), intent(inout) :: this

    nullify(this%mesh)

    if (allocated(this%static)) deallocate(this%static)
    if (allocated(this%state )) deallocate(this%state )

    this%initialized = .false.

  end subroutine dycore_clear

  subroutine dycore_final(this)

    type(dycore_type), intent(inout) :: this

    call this%clear()

  end subroutine dycore_final

end module dycore_mod
