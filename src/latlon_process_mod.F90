module latlon_process_mod

  use mpi
  use flogger
  use mesh_const_mod

  implicit none

  private

  public proc

  type process_neighbor_type
    integer :: id             = MPI_PROC_NULL
    integer :: cart_id        = MPI_PROC_NULL
    integer :: orient         = 0
    integer :: hwx            = 0
  end type process_neighbor_type

  type latlon_process_type
    logical :: initialized    = .false.
    logical :: mpi_master     = .true.
    integer :: comm           = MPI_COMM_NULL
    integer :: group          = MPI_GROUP_NULL
    integer :: np             = 0
    integer :: id             = MPI_PROC_NULL
    integer :: cart_comm      = MPI_COMM_NULL
    integer :: cart_group     = MPI_GROUP_NULL
    integer :: cart_id        = MPI_PROC_NULL
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes
  contains
    procedure :: init           => latlon_process_init
    procedure :: decomp_domain  => latlon_process_decomp_domain
    procedure :: is_root        => latlon_process_is_root
    final :: latlon_process_final
  end type latlon_process_type

  type(latlon_process_type), allocatable :: proc

contains

  subroutine latlon_process_init(this, comm, npx, npy, periods)

    class(latlon_process_type), intent(inout) :: this
    integer, intent(in), optional :: comm
    integer, intent(in), optional :: npx
    integer, intent(in), optional :: npy
    integer, intent(in), optional :: periods(2)

    integer ierr, tmp_comm, tmp_id(1), i, j

    if (present(comm)) then
      this%comm = comm
      this%mpi_master = .false.
    else
      call MPI_INIT(ierr)
      this%comm = MPI_COMM_WORLD
      call log_notice('Initialize MPI calling.')
    end if
    call MPI_COMM_GROUP(this%comm, this%group, ierr)
    call MPI_COMM_SIZE(this%comm, this%np, ierr)
    call MPI_COMM_RANK(this%comm, this%id, ierr)

    if (merge(npx * npy == this%np, .false., present(npx) .and. present(npy))) then
      this%cart_dims = [npx, npy]
    else
      this%cart_dims = [1, this%np]
    end if

    call MPI_COMM_SPLIT(this%comm, 0, this%id, tmp_comm, ierr)
    if (present(periods)) then
      call MPI_CART_CREATE(tmp_comm, 2, this%cart_dims, periods, .true., this%cart_comm, ierr)
    else
      call MPI_CART_CREATE(tmp_comm, 2, this%cart_dims, [.true., .false.], .true., this%cart_comm, ierr)
    end if
    call MPI_COMM_GROUP(this%cart_comm, this%cart_group, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_RANK(this%cart_comm, this%cart_id, ierr)
    call MPI_CART_COORDS(this%cart_comm, this%cart_id, 2, this%cart_coords, ierr)

    ! Set neighborhood of the process.
    if (allocated(this%ngb)) deallocate(this%ngb)
    allocate(this%ngb(4))
    call MPI_CART_SHIFT(this%cart_comm, 0, 1, this%ngb(left)%cart_id, this%ngb(right)%cart_id, ierr)
    call MPI_CART_SHIFT(this%cart_comm, 1, 1, this%ngb(bottom)%cart_id, this%ngb(top)%cart_id, ierr)

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, size(this%ngb)
      if (this%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(this%cart_group, 1, [this%ngb(i)%cart_id], this%group, tmp_id, ierr)
        this%ngb(i)%id = tmp_id(1)
      end if
    end do

  end subroutine latlon_process_init

  subroutine latlon_process_decomp_domain(this)

    class(latlon_process_type), intent(inout) :: this



  end subroutine latlon_process_decomp_domain

  pure logical function latlon_process_is_root(this) result(res)

    class(latlon_process_type), intent(in) :: this

    res = this%id == 0

  end function latlon_process_is_root

  subroutine latlon_process_final(this)

    type(latlon_process_type), intent(inout) :: this

    integer ierr

    if (this%mpi_master) then
      call MPI_FINALIZE(ierr)
      call log_notice('Finalize MPI calling.')
    end if

  end subroutine latlon_process_final

  subroutine round_robin(dim, coord, n, is, ie)

    integer, intent(in) :: dim
    integer, intent(in) :: coord
    integer, intent(inout) :: n
    integer, intent(out) :: is ! Start from 1.
    integer, intent(out) :: ie ! Start from 1.

    integer res_n, tmp_n, i

    res_n = mod(n, dim)
    is = 1
    do i = 0, coord - 1
      if (res_n /= 0 .and. i < res_n) then
        tmp_n = n / dim + 1
      else
        tmp_n = n / dim
      end if
      is = is + tmp_n
    end do
    if (res_n /= 0 .and. coord < res_n) then
      n = n / dim + 1
    else
      n = n / dim
    end if
    ie = is + n - 1

  end subroutine round_robin

end module latlon_process_mod
