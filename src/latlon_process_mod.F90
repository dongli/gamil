module latlon_process_mod

  use mpi
  use flogger
  use mesh_const_mod
  use latlon_mesh_mod

  implicit none

  private

  public proc

  type process_neighbor_type
    integer :: id             = MPI_PROC_NULL
    integer :: cart_id        = MPI_PROC_NULL
    integer :: orient         = 0
    integer :: hwx            = 0
    integer :: hwy            = 0
    integer :: is             = 0
    integer :: ie             = 0
    integer :: js             = 0
    integer :: je             = 0
  contains
    procedure :: init => process_neighbor_init
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
    integer :: is, ie, nx
    integer :: js, je, ny
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
    integer, intent(in), optional :: npx        ! Number of processes along x-axis
    integer, intent(in), optional :: npy        ! Number of processes along y-axis
    integer, intent(in), optional :: periods(2)

    integer ierr, tmp_comm

    if (present(comm)) then
      this%comm = comm
      this%mpi_master = .false.
    else
      call MPI_INIT(ierr)
      this%comm = MPI_COMM_WORLD
    end if
    call MPI_COMM_GROUP(this%comm, this%group, ierr)
    call MPI_COMM_SIZE(this%comm, this%np, ierr)
    call MPI_COMM_RANK(this%comm, this%id, ierr)

    call log_notice('Initialize MPI environment.', pid=this%id)

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

  end subroutine latlon_process_init

  subroutine latlon_process_decomp_domain(this, mesh, nx, ny, nz, hwx, hwy, nw, neq, r)

    class(latlon_process_type), intent(inout) :: this
    type(latlon_mesh_type), intent(inout) :: mesh
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer, intent(in) :: hwx
    integer, intent(in) :: hwy
    integer, intent(in) :: nw
    integer, intent(in) :: neq
    real(r8), intent(in) :: r

    integer ierr, cart_coords(2), tmp_id(1), i

    ! Set neighborhood of the process.
    if (allocated(this%ngb)) deallocate(this%ngb)
    allocate(this%ngb(8))
    call MPI_CART_SHIFT(this%cart_comm, 0, 1, this%ngb(left)%cart_id, this%ngb(right)%cart_id, ierr)
    call MPI_CART_SHIFT(this%cart_comm, 1, 1, this%ngb(bottom)%cart_id, this%ngb(top)%cart_id, ierr)
    ! West-south neighbor
    cart_coords = [proc%cart_coords(1)-1,proc%cart_coords(2)-1]
    if (cart_coords(1) == -1) cart_coords(1) = proc%cart_dims(1) - 1
    if (cart_coords(2) /= -1) then
      call MPI_CART_RANK(proc%cart_comm, cart_coords, proc%ngb(left_bottom)%cart_id, ierr)
    end if
    ! West-north neighbor
    cart_coords = [proc%cart_coords(1)-1,proc%cart_coords(2)+1]
    if (cart_coords(1) == -1) cart_coords(1) = proc%cart_dims(1) - 1
    if (cart_coords(2) /= proc%cart_dims(2)) then
      call MPI_CART_RANK(proc%cart_comm, cart_coords, proc%ngb(left_top)%cart_id, ierr)
    end if
    ! East-south neighbor
    cart_coords = [proc%cart_coords(1)+1,proc%cart_coords(2)-1]
    if (cart_coords(1) == proc%cart_dims(1)) cart_coords(1) = 0
    if (cart_coords(2) /= -1) then
      call MPI_CART_RANK(proc%cart_comm, cart_coords, proc%ngb(right_bottom)%cart_id, ierr)
    end if
    ! East-north neighbor
    cart_coords = [proc%cart_coords(1)+1,proc%cart_coords(2)+1]
    if (cart_coords(1) == proc%cart_dims(1)) cart_coords(1) = 0
    if (cart_coords(2) /= proc%cart_dims(2)) then
      call MPI_CART_RANK(proc%cart_comm, cart_coords, proc%ngb(right_top)%cart_id, ierr)
    end if

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, size(this%ngb)
      if (this%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(this%cart_group, 1, [this%ngb(i)%cart_id], this%group, tmp_id, ierr)
        this%ngb(i)%id = tmp_id(1)
      end if
    end do

    call round_robin(proc%cart_dims(1), proc%cart_coords(1), nx, proc%nx, proc%is, proc%ie)
    call round_robin(proc%cart_dims(2), proc%cart_coords(2), ny, proc%ny, proc%js, proc%je)

    ! NOTE: We do not include the Poles as a cell.
    call mesh%init(nx, ny, nz, dx=pi2/nx, dy=pi/ny, hwx=hwx, hwy=hwy, nw=nw, neq=neq, r=r, &
                   ids=proc%is, ide=proc%ie, jds=proc%js, jde=proc%je)

    call this%ngb(left        )%init(left        , hwx, hwy, mesh%ims    , mesh%ids - 1, mesh%jds    , mesh%jde    )
    call this%ngb(right       )%init(right       , hwx, hwy, mesh%ide + 1, mesh%ime    , mesh%jds    , mesh%jde    )
    call this%ngb(bottom      )%init(bottom      , hwx, hwy, mesh%ids    , mesh%ide    , mesh%jms    , mesh%jds - 1)
    call this%ngb(top         )%init(top         , hwx, hwy, mesh%ids    , mesh%ide    , mesh%jde + 1, mesh%jme    )
    call this%ngb(left_bottom )%init(left_bottom , hwx, hwy, mesh%ims    , mesh%ids - 1, mesh%jms    , mesh%jds - 1)
    call this%ngb(left_top    )%init(left_top    , hwx, hwy, mesh%ims    , mesh%ids - 1, mesh%jde + 1, mesh%jme    )
    call this%ngb(right_bottom)%init(right_bottom, hwx, hwy, mesh%ide + 1, mesh%ime    , mesh%jms    , mesh%jds - 1)
    call this%ngb(right_top   )%init(right_top   , hwx, hwy, mesh%ide + 1, mesh%ime    , mesh%jde + 1, mesh%jme    )

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
      call log_notice('Finalize MPI environment.', pid=this%id)
    end if

  end subroutine latlon_process_final

  subroutine round_robin(dim, coord, n0, n, is, ie)

    integer, intent(in) :: dim
    integer, intent(in) :: coord
    integer, intent(in) :: n0
    integer, intent(out) :: n
    integer, intent(out) :: is ! Start from 1.
    integer, intent(out) :: ie ! Start from 1.

    integer res_n, tmp_n, i

    res_n = mod(n0, dim)
    is = 1
    do i = 0, coord - 1
      if (res_n /= 0 .and. i < res_n) then
        tmp_n = n0 / dim + 1
      else
        tmp_n = n0 / dim
      end if
      is = is + tmp_n
    end do
    if (res_n /= 0 .and. coord < res_n) then
      n = n0 / dim + 1
    else
      n = n0 / dim
    end if
    ie = is + n - 1

  end subroutine round_robin

  subroutine process_neighbor_init(this, orient, hwx, hwy, is, ie, js, je)

    class(process_neighbor_type), intent(inout) :: this
    integer, intent(in) :: orient
    integer, intent(in) :: hwx
    integer, intent(in) :: hwy
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: js
    integer, intent(in) :: je

    this%orient = orient
    this%hwx    = hwx
    this%hwy    = hwy
    this%is     = is
    this%ie     = ie
    this%js     = js
    this%je     = je

  end subroutine process_neighbor_init

end module latlon_process_mod
