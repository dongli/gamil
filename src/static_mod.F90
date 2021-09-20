module static_mod

  use const_mod
  use kinds_mod
  use latlon_mesh_mod
  use latlon_array_mod

  implicit none

  type static_type
    logical :: initialized = .false.
    type(latlon_array_type), allocatable :: array
    real(r8), pointer, dimension(:,:) :: zs
    real(r8), pointer, dimension(:,:) :: f1, f2 ! Coriolis coefficent
  contains
    procedure :: init     => static_init
    procedure :: clear    => static_clear
    final :: static_final
  end type static_type

contains

  subroutine static_init(this, mesh)

    class(static_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh

    integer i, j

    call this%clear()

    allocate(this%array)

    call this%array%init(mesh)
    call this%array%add_var('zs', 'Surface height', 'm', loc='CA', with_halo='T', fill_halo='T', output='T', only_2d=.true.)
    call this%array%add_var('f1', loc='C')
    call this%array%add_var('f2', loc='C')
    call this%array%allocate_arrays()

    call this%array%get_array(this%zs, var_name='zs', loc='CA')
    call this%array%get_array(this%f1, var_name='f1', loc='C' )
    call this%array%get_array(this%f2, var_name='f2', loc='C' )

    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        this%f1(i,j) = 2 * omega * sin(mesh%lat(1,i,j))
        this%f2(i,j) = 2 * omega * cos(mesh%lat(1,i,j))
      end do
    end do

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
