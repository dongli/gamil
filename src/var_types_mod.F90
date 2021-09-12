module var_types_mod

  implicit none

  integer, parameter :: name_len = 10

  type var_info_type
    character(name_len) name
    character(100) long_name
    character(30 ) units
    character(2  ) loc        ! Variable location:
                              ! - C : Cell
                              ! - CA: Cell average
                              ! - CQ: Cell quadrature
                              ! - L : Left edge
                              ! - LQ: Left edge quadrature
                              ! - R : Right edge
                              ! - RQ: Right edge quadrature
                              ! - T : Top edge
                              ! - TQ: Top edge quadrature
                              ! - B : Bottom edge
                              ! - BQ: Bottom edge quadrature
    character(30) :: tag = '' ! Used to tag variables that can be used together, e.g., tag='fcst', tag='tracer'
    logical :: with_halo = .false.
    logical :: fill_halo = .false.
    logical :: output    = .false.
    logical :: only_2d   = .false.
    integer :: array_idx = 0
  end type var_info_type

  type var_stack_type
    integer :: size         = 0
    type(var_info_type), allocatable :: var_info(:)
  contains
    procedure :: init    => var_stack_init
    procedure :: clear   => var_stack_clear
    procedure :: append  => var_stack_append
    procedure :: reorder => var_stack_reorder
    procedure :: name_at => var_stack_name_at
    final :: var_stack_final
  end type var_stack_type

contains

  subroutine var_stack_init(this, max_var)

    class(var_stack_type), intent(inout) :: this
    integer, intent(in), optional :: max_var

    call this%clear()

    allocate(this%var_info(merge(max_var, 100, present(max_var))))

  end subroutine var_stack_init

  subroutine var_stack_clear(this)

    class(var_stack_type), intent(inout) :: this

    if (allocated(this%var_info)) deallocate(this%var_info)

  end subroutine var_stack_clear

  subroutine var_stack_append(this, name, long_name, units, loc, tag, with_halo, fill_halo, output, only_2d)

    class(var_stack_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: long_name
    character(*), intent(in) :: units
    character(*), intent(in) :: loc
    character(*), intent(in), optional :: tag
    character(*), intent(in), optional :: with_halo
    character(*), intent(in), optional :: fill_halo
    character(*), intent(in), optional :: output
    logical, intent(in), optional :: only_2d

    type(var_info_type), allocatable :: tmp(:)
    integer i, loc_is, loc_ie, halo_i, n

    n = 1
    do i = 1, len_trim(loc)
      if (loc(i:i) == ':') n = n + 1
    end do

    this%size = this%size + n

    if (this%size > size(this%var_info)) then
      ! Enlarge var_info array.
      allocate(tmp(this%size + 100))
      do i = 1, size(this%var_info)
        tmp(i) = this%var_info(i)
      end do
      if (allocated(this%var_info)) deallocate(this%var_info)
      allocate(this%var_info(size(tmp)), source=tmp)
      deallocate(tmp)
    end if
    loc_is = 1
    halo_i = 1
    do i = this%size - n + 1, this%size
      this%var_info(i)%name      = name
      this%var_info(i)%long_name = long_name
      this%var_info(i)%units     = units
      this%var_info(i)%only_2d   = merge(only_2d, .false., present(only_2d))
      loc_ie = loc_is + index(loc(loc_is:len_trim(loc)), ':') - 2
      loc_ie = merge(len_trim(loc), loc_ie, loc_ie < loc_is)
      this%var_info(i)%loc       = loc(loc_is:loc_ie)
      if (present(tag)) then
        this%var_info(i)%tag     = tag
      end if
      if (present(with_halo)) then
        this%var_info(i)%with_halo = with_halo(halo_i:halo_i) == 'T'
      end if
      if (present(fill_halo)) then
        this%var_info(i)%fill_halo = fill_halo(halo_i:halo_i) == 'T'
      end if
      if (present(output)) then
        this%var_info(i)%output = output(halo_i:halo_i) == 'T'
      end if
      loc_is = loc_ie + 2
      halo_i = halo_i + 2
    end do

  end subroutine var_stack_append

  subroutine var_stack_reorder(this)

    class(var_stack_type), intent(inout) :: this

    type(var_info_type) tmp(this%size)
    character(2) :: loc(11) = ['C ', 'CA', 'CQ', 'L ', 'LQ', 'R ', 'RQ', 'T ', 'TQ', 'B ', 'BQ']
    integer i, j, k

    k = 1
    do j = 1, size(loc)
      do i = 1, this%size
        if (this%var_info(i)%loc == loc(j)) then
          tmp(k) = this%var_info(i); k = k + 1
        end if
      end do
    end do
    if (allocated(this%var_info)) deallocate(this%var_info)
    allocate(this%var_info(size(tmp)), source=tmp)

  end subroutine var_stack_reorder

  character(name_len) function var_stack_name_at(this, idx) result(res)

    class(var_stack_type), intent(in) :: this
    integer, intent(in) :: idx

    res = this%var_info(idx)%name

  end function var_stack_name_at

  subroutine var_stack_final(this)

    type(var_stack_type), intent(inout) :: this

    call this%clear()

  end subroutine var_stack_final

end module var_types_mod
