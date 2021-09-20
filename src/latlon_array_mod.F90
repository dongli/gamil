module latlon_array_mod

  use flogger
  use fiona
  use kinds_mod
  use latlon_process_mod
  use latlon_mesh_mod
  use latlon_halo_mod

  type latlon_array_type
    logical :: initialized = .false.
    type(var_stack_type) var_stack
    type(latlon_mesh_type), pointer :: mesh => null()
    type(latlon_halo_type), allocatable :: halo(:)
    ! Arrays for storing user data
    ! Dimension columns:
    ! 1. point location
    ! 2. x-axis
    ! 3. y-axis
    ! 4. vertical level
    ! 5. variable
    ! Subscript:
    ! c. cell
    ! l. left cell edge quadtrature
    ! r. right cell edge quadtrature
    ! t. top cell edge quadtrature
    ! b. bottom cell edge quadtrature
    ! q. cell quadtrature
    ! h. with halo elements
    !                                1 2 3 4 5
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_c_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_ca_h
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_cq_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_c
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_ca
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_cq
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_l_h
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_lq_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_l
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_lq
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_r_h
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_rq_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_r
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_rq
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_t_h
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_tq_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_t
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_tq
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_b_h
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_bq_h
    real(r8), allocatable, dimension(  :,:,  :) :: a_2d_b
    real(r8), allocatable, dimension(:,:,:,  :) :: a_2d_bq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_c_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_ca_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_cq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_c
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_ca
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_cq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_l_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_lq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_l
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_lq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_r_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_rq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_r
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_rq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_t_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_tq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_t
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_tq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_b_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_bq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_b
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_bq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_u_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_uq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_u
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_uq
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_d_h
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_dq_h
    real(r8), allocatable, dimension(  :,:,:,:) :: a_3d_d
    real(r8), allocatable, dimension(:,:,:,:,:) :: a_3d_dq
  contains
    procedure :: init            => latlon_array_init
    procedure :: clear           => latlon_array_clear
    procedure :: add_var         => latlon_array_add_var
    procedure :: var_idx         => latlon_array_var_idx
    procedure :: allocate_arrays => latlon_array_allocate_arrays
    generic   :: get_array       => latlon_array_get_array_2d_1, &
                                    latlon_array_get_array_2d_2, &
                                    latlon_array_get_array_3d_1, &
                                    latlon_array_get_array_3d_2, &
                                    latlon_array_get_array_4d_1, &
                                    latlon_array_get_array_4d_2, &
                                    latlon_array_get_array_5d
    procedure, private :: latlon_array_get_array_2d_1
    procedure, private :: latlon_array_get_array_2d_2
    procedure, private :: latlon_array_get_array_3d_1
    procedure, private :: latlon_array_get_array_3d_2
    procedure, private :: latlon_array_get_array_4d_1
    procedure, private :: latlon_array_get_array_4d_2
    procedure, private :: latlon_array_get_array_5d
    procedure :: create_dataset  => latlon_array_create_dataset
    procedure :: append_dataset  => latlon_array_append_dataset
    procedure :: write           => latlon_array_write
    final :: latlon_array_final
  end type latlon_array_type

contains

  subroutine latlon_array_init(this, mesh)

    class(latlon_array_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh

    integer i, dtype

    call this%clear()

    this%mesh => mesh
    call this%var_stack%init()

    ! Initialize halo objects.
    select case (r8)
    case (4)
      dtype = MPI_REAL
    case (8)
      dtype = MPI_DOUBLE
    end select
    allocate(this%halo(size(proc%ngb)))
    do i = 1, size(this%halo)
      call this%halo(i)%init(mesh, proc%ngb(i)%orient, dtype, host_id=proc%id, proc_id=proc%ngb(i)%id)
    end do

    this%initialized = .true.

  end subroutine latlon_array_init

  subroutine latlon_array_clear(this)

    class(latlon_array_type), intent(inout) :: this

    call this%var_stack%clear()

    if (allocated(this%halo)) deallocate(this%halo)

    if (allocated(this%a_2d_c_h )) deallocate(this%a_2d_c_h )
    if (allocated(this%a_2d_c   )) deallocate(this%a_2d_c   )
    if (allocated(this%a_2d_ca_h)) deallocate(this%a_2d_ca_h)
    if (allocated(this%a_2d_ca  )) deallocate(this%a_2d_ca  )
    if (allocated(this%a_2d_cq_h)) deallocate(this%a_2d_cq_h)
    if (allocated(this%a_2d_cq  )) deallocate(this%a_2d_cq  )
    if (allocated(this%a_3d_c_h )) deallocate(this%a_3d_c_h )
    if (allocated(this%a_3d_c   )) deallocate(this%a_3d_c   )
    if (allocated(this%a_3d_ca_h)) deallocate(this%a_3d_ca_h)
    if (allocated(this%a_3d_ca  )) deallocate(this%a_3d_ca  )
    if (allocated(this%a_3d_cq_h)) deallocate(this%a_3d_cq_h)
    if (allocated(this%a_3d_cq  )) deallocate(this%a_3d_cq  )
    if (allocated(this%a_2d_l_h )) deallocate(this%a_2d_l_h )
    if (allocated(this%a_2d_l   )) deallocate(this%a_2d_l   )
    if (allocated(this%a_2d_lq_h)) deallocate(this%a_2d_lq_h)
    if (allocated(this%a_2d_lq  )) deallocate(this%a_2d_lq  )
    if (allocated(this%a_3d_l_h )) deallocate(this%a_3d_l_h )
    if (allocated(this%a_3d_l   )) deallocate(this%a_3d_l   )
    if (allocated(this%a_3d_lq_h)) deallocate(this%a_3d_lq_h)
    if (allocated(this%a_3d_lq  )) deallocate(this%a_3d_lq  )
    if (allocated(this%a_2d_r_h )) deallocate(this%a_2d_r_h )
    if (allocated(this%a_2d_r   )) deallocate(this%a_2d_r   )
    if (allocated(this%a_2d_rq_h)) deallocate(this%a_2d_rq_h)
    if (allocated(this%a_2d_rq  )) deallocate(this%a_2d_rq  )
    if (allocated(this%a_3d_r_h )) deallocate(this%a_3d_r_h )
    if (allocated(this%a_3d_r   )) deallocate(this%a_3d_r   )
    if (allocated(this%a_3d_rq_h)) deallocate(this%a_3d_rq_h)
    if (allocated(this%a_3d_rq  )) deallocate(this%a_3d_rq  )
    if (allocated(this%a_2d_t_h )) deallocate(this%a_2d_t_h )
    if (allocated(this%a_2d_t   )) deallocate(this%a_2d_t   )
    if (allocated(this%a_2d_tq_h)) deallocate(this%a_2d_tq_h)
    if (allocated(this%a_2d_tq  )) deallocate(this%a_2d_tq  )
    if (allocated(this%a_3d_t_h )) deallocate(this%a_3d_t_h )
    if (allocated(this%a_3d_t   )) deallocate(this%a_3d_t   )
    if (allocated(this%a_3d_tq_h)) deallocate(this%a_3d_tq_h)
    if (allocated(this%a_3d_tq  )) deallocate(this%a_3d_tq  )
    if (allocated(this%a_2d_b_h )) deallocate(this%a_2d_b_h )
    if (allocated(this%a_2d_b   )) deallocate(this%a_2d_b   )
    if (allocated(this%a_2d_bq_h)) deallocate(this%a_2d_bq_h)
    if (allocated(this%a_2d_bq  )) deallocate(this%a_2d_bq  )
    if (allocated(this%a_3d_b_h )) deallocate(this%a_3d_b_h )
    if (allocated(this%a_3d_b   )) deallocate(this%a_3d_b   )
    if (allocated(this%a_3d_bq_h)) deallocate(this%a_3d_bq_h)
    if (allocated(this%a_3d_bq  )) deallocate(this%a_3d_bq  )
    if (allocated(this%a_3d_u_h )) deallocate(this%a_3d_u_h )
    if (allocated(this%a_3d_u   )) deallocate(this%a_3d_u   )
    if (allocated(this%a_3d_uq_h)) deallocate(this%a_3d_uq_h)
    if (allocated(this%a_3d_uq  )) deallocate(this%a_3d_uq  )
    if (allocated(this%a_3d_d_h )) deallocate(this%a_3d_d_h )
    if (allocated(this%a_3d_d   )) deallocate(this%a_3d_d   )
    if (allocated(this%a_3d_dq_h)) deallocate(this%a_3d_dq_h)
    if (allocated(this%a_3d_dq  )) deallocate(this%a_3d_dq  )

    this%initialized = .false.

  end subroutine latlon_array_clear

  subroutine latlon_array_add_var(this, name, long_name, units, loc, tag, with_halo, fill_halo, output, only_2d)

    class(latlon_array_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units
    character(*), intent(in), optional :: loc
    character(*), intent(in), optional :: tag
    character(*), intent(in), optional :: output
    character(*), intent(in), optional :: with_halo
    character(*), intent(in), optional :: fill_halo
    logical, intent(in), optional :: only_2d

    call this%var_stack%append(name, long_name, units, loc, tag, with_halo, fill_halo, output, only_2d)

  end subroutine latlon_array_add_var

  integer function latlon_array_var_idx(this, name, loc) result(res)

    class(latlon_array_type), intent(in) :: this
    character(*), intent(in) :: name
    character(*), intent(in) :: loc

    do res = 1, this%var_stack%size
      if (this%var_stack%var_info(res)%name == name .and. this%var_stack%var_info(res)%loc == loc) return
    end do
    res = 0

  end function latlon_array_var_idx

  subroutine latlon_array_allocate_arrays(this, fill_value)

    class(latlon_array_type), intent(inout) :: this
    real(r8), intent(in), optional :: fill_value

    integer i
    integer n_2d_c_h, n_2d_c, n_3d_c_h, n_3d_c, n_2d_cq_h, n_2d_cq, n_3d_cq_h, n_3d_cq
    integer n_2d_l_h, n_2d_l, n_3d_l_h, n_3d_l, n_2d_lq_h, n_2d_lq, n_3d_lq_h, n_3d_lq
    integer n_2d_r_h, n_2d_r, n_3d_r_h, n_3d_r, n_2d_rq_h, n_2d_rq, n_3d_rq_h, n_3d_rq
    integer n_2d_t_h, n_2d_t, n_3d_t_h, n_3d_t, n_2d_tq_h, n_2d_tq, n_3d_tq_h, n_3d_tq
    integer n_2d_b_h, n_2d_b, n_3d_b_h, n_3d_b, n_2d_bq_h, n_2d_bq, n_3d_bq_h, n_3d_bq
    integer                                     n_2d_ca_h, n_2d_ca, n_3d_ca_h, n_3d_ca
    integer n_3d_u_h, n_3d_u, n_3d_uq_h, n_3d_uq
    integer n_3d_d_h, n_3d_d, n_3d_dq_h, n_3d_dq
    integer ims, ime, ids, ide, jms, jme, jds, jde, kms, kme
    integer pqs, pqe, pes, pee
    real(r8) fill_value_

    call this%var_stack%reorder()

    fill_value_ = merge(fill_value, 0.0_r8, present(fill_value))

    n_2d_c_h = 0; n_2d_c = 0; n_3d_c_h = 0; n_3d_c = 0; n_2d_cq_h = 0; n_2d_cq = 0; n_3d_cq_h = 0; n_3d_cq = 0
    n_2d_l_h = 0; n_2d_l = 0; n_3d_l_h = 0; n_3d_l = 0; n_2d_lq_h = 0; n_2d_lq = 0; n_3d_lq_h = 0; n_3d_lq = 0
    n_2d_r_h = 0; n_2d_r = 0; n_3d_r_h = 0; n_3d_r = 0; n_2d_rq_h = 0; n_2d_rq = 0; n_3d_rq_h = 0; n_3d_rq = 0
    n_2d_t_h = 0; n_2d_t = 0; n_3d_t_h = 0; n_3d_t = 0; n_2d_tq_h = 0; n_2d_tq = 0; n_3d_tq_h = 0; n_3d_tq = 0
    n_2d_b_h = 0; n_2d_b = 0; n_3d_b_h = 0; n_3d_b = 0; n_2d_bq_h = 0; n_2d_bq = 0; n_3d_bq_h = 0; n_3d_bq = 0
                                                        n_2d_ca_h = 0; n_2d_ca = 0; n_3d_ca_h = 0; n_3d_ca = 0
    n_3d_u_h = 0; n_3d_u = 0; n_3d_uq_h = 0; n_3d_uq = 0
    n_3d_d_h = 0; n_3d_d = 0; n_3d_dq_h = 0; n_3d_dq = 0
    do i = 1, this%var_stack%size
      select case (this%var_stack%var_info(i)%loc)
      case ('C')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_c_h = n_2d_c_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_c_h
          else
            n_3d_c_h = n_3d_c_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_c_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_c = n_2d_c + 1
            this%var_stack%var_info(i)%array_idx = n_2d_c
          else
            n_3d_c = n_3d_c + 1
            this%var_stack%var_info(i)%array_idx = n_3d_c
          end if
        end if
      case ('CA')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_ca_h = n_2d_ca_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_ca_h
          else
            n_3d_ca_h = n_3d_ca_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_ca_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_ca = n_2d_ca + 1
            this%var_stack%var_info(i)%array_idx = n_2d_ca
          else
            n_3d_ca = n_3d_ca + 1
            this%var_stack%var_info(i)%array_idx = n_3d_ca
          end if
        end if
      case ('CQ')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_cq_h = n_2d_cq_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_cq_h
          else
            n_3d_cq_h = n_3d_cq_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_cq_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_cq = n_2d_cq + 1
            this%var_stack%var_info(i)%array_idx = n_2d_cq
          else
            n_3d_cq = n_3d_cq + 1
            this%var_stack%var_info(i)%array_idx = n_3d_cq
          end if
        end if
      case ('L')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_l_h = n_2d_l_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_l_h
          else
            n_3d_l_h = n_3d_l_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_l_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_l = n_2d_l + 1
            this%var_stack%var_info(i)%array_idx = n_2d_l
          else
            n_3d_l = n_3d_l + 1
            this%var_stack%var_info(i)%array_idx = n_3d_l
          end if
        end if
      case ('LQ')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_lq_h = n_2d_lq_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_lq_h
          else
            n_3d_lq_h = n_3d_lq_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_lq_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_lq = n_2d_lq + 1
            this%var_stack%var_info(i)%array_idx = n_2d_lq
          else
            n_3d_lq = n_3d_lq + 1
            this%var_stack%var_info(i)%array_idx = n_3d_lq
          end if
        end if
      case ('R')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_r_h = n_2d_r_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_r_h
          else
            n_3d_r_h = n_3d_r_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_r_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_r = n_2d_r + 1
            this%var_stack%var_info(i)%array_idx = n_2d_r
          else
            n_3d_r = n_3d_r + 1
            this%var_stack%var_info(i)%array_idx = n_3d_r
          end if
        end if
      case ('RQ')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_rq_h = n_2d_rq_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_rq_h
          else
            n_3d_rq_h = n_3d_rq_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_rq_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_rq = n_2d_rq + 1
            this%var_stack%var_info(i)%array_idx = n_2d_rq
          else
            n_3d_rq = n_3d_rq + 1
            this%var_stack%var_info(i)%array_idx = n_3d_rq
          end if
        end if
      case ('T')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_t_h = n_2d_t_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_t_h
          else
            n_3d_t_h = n_3d_t_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_t_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_t = n_2d_t + 1
            this%var_stack%var_info(i)%array_idx = n_2d_t
          else
            n_3d_t = n_3d_t + 1
            this%var_stack%var_info(i)%array_idx = n_3d_t
          end if
        end if
      case ('TQ')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_tq_h = n_2d_tq_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_tq_h
          else
            n_3d_tq_h = n_3d_tq_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_tq_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_tq = n_2d_tq + 1
            this%var_stack%var_info(i)%array_idx = n_2d_tq
          else
            n_3d_tq = n_3d_tq + 1
            this%var_stack%var_info(i)%array_idx = n_3d_tq
          end if
        end if
      case ('B')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_b_h = n_2d_b_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_b_h
          else
            n_3d_b_h = n_3d_b_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_b_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_b = n_2d_b + 1
            this%var_stack%var_info(i)%array_idx = n_2d_b
          else
            n_3d_b = n_3d_b + 1
            this%var_stack%var_info(i)%array_idx = n_3d_b
          end if
        end if
      case ('BQ')
        if (this%var_stack%var_info(i)%with_halo) then
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_bq_h = n_2d_bq_h + 1
            this%var_stack%var_info(i)%array_idx = n_2d_bq_h
          else
            n_3d_bq_h = n_3d_bq_h + 1
            this%var_stack%var_info(i)%array_idx = n_3d_bq_h
          end if
        else
          if (this%var_stack%var_info(i)%only_2d) then
            n_2d_bq = n_2d_bq + 1
            this%var_stack%var_info(i)%array_idx = n_2d_bq
          else
            n_3d_bq = n_3d_bq + 1
            this%var_stack%var_info(i)%array_idx = n_3d_bq
          end if
        end if
      case ('U')
        if (this%var_stack%var_info(i)%with_halo) then
          n_3d_u_h = n_3d_u_h + 1
          this%var_stack%var_info(i)%array_idx = n_3d_u_h
        else
          n_3d_u = n_3d_u + 1
          this%var_stack%var_info(i)%array_idx = n_3d_u
        end if
      case ('UQ')
        if (this%var_stack%var_info(i)%with_halo) then
          n_3d_uq_h = n_3d_uq_h + 1
          this%var_stack%var_info(i)%array_idx = n_3d_uq_h
        else
          n_3d_uq = n_3d_uq + 1
          this%var_stack%var_info(i)%array_idx = n_3d_uq
        end if
      case ('D')
        if (this%var_stack%var_info(i)%with_halo) then
          n_3d_d_h = n_3d_d_h + 1
          this%var_stack%var_info(i)%array_idx = n_3d_d_h
        else
          n_3d_d = n_3d_d + 1
          this%var_stack%var_info(i)%array_idx = n_3d_d
        end if
      case ('DQ')
        if (this%var_stack%var_info(i)%with_halo) then
          n_3d_dq_h = n_3d_dq_h + 1
          this%var_stack%var_info(i)%array_idx = n_3d_dq_h
        else
          n_3d_dq = n_3d_dq + 1
          this%var_stack%var_info(i)%array_idx = n_3d_dq
        end if
      end select
    end do

    ims = this%mesh%ims; ime = this%mesh%ime; ids = this%mesh%ids; ide = this%mesh%ide
    jms = this%mesh%jms; jme = this%mesh%jme; jds = this%mesh%jds; jde = this%mesh%jde
    kms = this%mesh%kms; kme = this%mesh%kme
    pqs = this%mesh%pqs; pqe = this%mesh%pqe

    ! Center
    if (n_2d_c_h > 0) then
      allocate(this%a_2d_c_h(ims:ime,jms:jme,n_2d_c_h))
      this%a_2d_c_h = fill_value_
    end if
    if (n_3d_c_h > 0) then
      allocate(this%a_3d_c_h(ims:ime,jms:jme,kms:kme,n_3d_c_h))
      this%a_3d_c_h = fill_value_
    end if
    if (n_2d_c > 0) then
      allocate(this%a_2d_c(ids:ide,jds:jde,n_2d_c))
      this%a_2d_c = fill_value_
    end if
    if (n_3d_c > 0) then
      allocate(this%a_3d_c(ids:ide,jds:jde,kms:kme,n_3d_c))
      this%a_3d_c = fill_value_
    end if
    if (n_2d_ca_h > 0) then
      allocate(this%a_2d_ca_h(ims:ime,jms:jme,n_2d_ca_h))
      this%a_2d_ca_h = fill_value_
    end if
    if (n_3d_ca_h > 0) then
      allocate(this%a_3d_ca_h(ims:ime,jms:jme,kms:kme,n_3d_ca_h))
      this%a_3d_ca_h = fill_value_
    end if
    if (n_2d_ca > 0) then
      allocate(this%a_2d_ca(ids:ide,jds:jde,n_2d_ca))
      this%a_2d_ca = fill_value_
    end if
    if (n_3d_ca > 0) then
      allocate(this%a_3d_ca(ids:ide,jds:jde,kms:kme,n_3d_ca))
      this%a_3d_ca = fill_value_
    end if
    if (n_2d_cq_h > 0) then
      allocate(this%a_2d_cq_h(pqs:pqe,ims:ime,jms:jme,n_2d_cq_h))
      this%a_2d_cq_h = fill_value_
    end if
    if (n_3d_cq_h > 0) then
      allocate(this%a_3d_cq_h(pqs:pqe,ims:ime,jms:jme,kms:kme,n_3d_cq_h))
      this%a_3d_cq_h = fill_value_
    end if
    if (n_2d_cq > 0) then
      allocate(this%a_2d_cq(pqs:pqe,ids:ide,jds:jde,n_2d_cq))
      this%a_2d_cq = fill_value_
    end if
    if (n_3d_cq > 0) then
      allocate(this%a_3d_cq(pqs:pqe,ids:ide,jds:jde,kms:kme,n_3d_cq))
      this%a_3d_cq = fill_value_
    end if
    ! Left edge
    if (n_2d_l_h > 0) then
      allocate(this%a_2d_l_h(ims:ime,jms:jme,n_2d_l_h))
      this%a_2d_l_h = fill_value_
    end if
    if (n_3d_l_h > 0) then
      allocate(this%a_3d_l_h(ims:ime,jms:jme,kms:kme,n_3d_l_h))
      this%a_3d_l_h = fill_value_
    end if
    if (n_2d_l > 0) then
      allocate(this%a_2d_l(ids:ide+1,jds:jde,n_2d_l))
      this%a_2d_l = fill_value_
    end if
    if (n_3d_l > 0) then
      allocate(this%a_3d_l(ids:ide+1,jds:jde,kms:kme,n_3d_l))
      this%a_3d_l = fill_value_
    end if
    pes = this%mesh%pes(left)
    pee = this%mesh%pee(left)
    if (n_2d_lq_h > 0) then
      allocate(this%a_2d_lq_h(pes:pee,ims:ime,jms:jme,n_2d_lq_h))
      this%a_2d_lq_h = fill_value_
    end if
    if (n_3d_lq_h > 0) then
      allocate(this%a_3d_lq_h(pes:pee,ims:ime,jms:jme,kms:kme,n_3d_lq_h))
      this%a_3d_lq_h = fill_value_
    end if
    if (n_2d_lq > 0) then
      allocate(this%a_2d_lq(pes:pee,ids:ide+1,jds:jde,n_2d_lq))
      this%a_2d_lq = fill_value_
    end if
    if (n_3d_lq > 0) then
      allocate(this%a_3d_lq(pes:pee,ids:ide+1,jds:jde,kms:kme,n_3d_lq))
      this%a_3d_lq = fill_value_
    end if
    ! Right edge
    if (n_2d_r_h > 0) then
      allocate(this%a_2d_r_h(ims:ime,jms:jme,n_2d_r_h))
      this%a_2d_r_h = fill_value_
    end if
    if (n_3d_r_h > 0) then
      allocate(this%a_3d_r_h(ims:ime,jms:jme,kms:kme,n_3d_r_h))
      this%a_3d_r_h = fill_value_
    end if
    if (n_2d_r > 0) then
      allocate(this%a_2d_r(ids:ide+1,jds:jde,n_2d_r))
      this%a_2d_r = fill_value_
    end if
    if (n_3d_r > 0) then
      allocate(this%a_3d_r(ids:ide+1,jds:jde,kms:kme,n_3d_r))
      this%a_3d_r = fill_value_
    end if
    pes = this%mesh%pes(right)
    pee = this%mesh%pee(right)
    if (n_2d_rq_h > 0) then
      allocate(this%a_2d_rq_h(pes:pee,ims:ime,jms:jme,n_2d_rq_h))
      this%a_2d_rq_h = fill_value_
    end if
    if (n_3d_rq_h > 0) then
      allocate(this%a_3d_rq_h(pes:pee,ims:ime,jms:jme,kms:kme,n_3d_rq_h))
      this%a_3d_rq_h = fill_value_
    end if
    if (n_2d_rq > 0) then
      allocate(this%a_2d_rq(pes:pee,ids:ide+1,jds:jde,n_2d_rq))
      this%a_2d_rq = fill_value_
    end if
    if (n_3d_rq > 0) then
      allocate(this%a_3d_rq(pes:pee,ids:ide+1,jds:jde,kms:kme,n_3d_rq))
      this%a_3d_rq = fill_value_
    end if
    ! Top edge quadtrature
    if (n_2d_t_h > 0) then
      allocate(this%a_2d_t_h(ims:ime,jms:jme,n_2d_t_h))
      this%a_2d_t_h = fill_value_
    end if
    if (n_3d_t_h > 0) then
      allocate(this%a_3d_t_h(ims:ime,jms:jme,kms:kme,n_3d_t_h))
      this%a_3d_t_h = fill_value_
    end if
    if (n_2d_t > 0) then
      allocate(this%a_2d_t(ids:ide,jds:jde+1,n_2d_t))
      this%a_2d_t = fill_value_
    end if
    if (n_3d_t > 0) then
      allocate(this%a_3d_t(ids:ide,jds:jde+1,kms:kme,n_3d_t))
      this%a_3d_t = fill_value_
    end if
    pes = this%mesh%pes(top)
    pee = this%mesh%pee(top)
    if (n_2d_tq_h > 0) then
      allocate(this%a_2d_tq_h(pes:pee,ims:ime,jms:jme,n_2d_tq_h))
      this%a_2d_tq_h = fill_value_
    end if
    if (n_3d_tq_h > 0) then
      allocate(this%a_3d_tq_h(pes:pee,ims:ime,jms:jme,kms:kme,n_3d_tq_h))
      this%a_3d_tq_h = fill_value_
    end if
    if (n_2d_tq > 0) then
      allocate(this%a_2d_tq(pes:pee,ids:ide,jds:jde+1,n_2d_tq))
      this%a_2d_tq = fill_value_
    end if
    if (n_3d_tq > 0) then
      allocate(this%a_3d_tq(pes:pee,ids:ide,jds:jde+1,kms:kme,n_3d_tq))
      this%a_3d_tq = fill_value_
    end if
    ! Bottom cell edge quadtrature
    if (n_2d_b_h > 0) then
      allocate(this%a_2d_b_h(ims:ime,jms:jme,n_2d_b_h))
      this%a_2d_b_h = fill_value_
    end if
    if (n_3d_b_h > 0) then
      allocate(this%a_3d_b_h(ims:ime,jms:jme,kms:kme,n_3d_b_h))
      this%a_3d_b_h = fill_value_
    end if
    if (n_2d_b > 0) then
      allocate(this%a_2d_b(ids:ide,jds:jde+1,n_2d_b))
      this%a_2d_b = fill_value_
    end if
    if (n_3d_b > 0) then
      allocate(this%a_3d_b(ids:ide,jds:jde+1,kms:kme,n_3d_b))
      this%a_3d_b = fill_value_
    end if
    pes = this%mesh%pes(bottom)
    pee = this%mesh%pee(bottom)
    if (n_2d_bq_h > 0) then
      allocate(this%a_2d_bq_h(pes:pee,ims:ime,jms:jme,n_2d_bq_h))
      this%a_2d_bq_h = fill_value_
    end if
    if (n_3d_bq_h > 0) then
      allocate(this%a_3d_bq_h(pes:pee,ims:ime,jms:jme,kms:kme,n_3d_bq_h))
      this%a_3d_bq_h = fill_value_
    end if
    if (n_2d_bq > 0) then
      allocate(this%a_2d_bq(pes:pee,ids:ide,jds:jde+1,n_2d_bq))
      this%a_2d_bq = fill_value_
    end if
    if (n_3d_bq > 0) then
      allocate(this%a_3d_bq(pes:pee,ids:ide,jds:jde+1,kms:kme,n_3d_bq))
      this%a_3d_bq = fill_value_
    end if
    ! Up cell edge quadtrature
    if (n_3d_u_h > 0) then
      allocate(this%a_3d_u_h(ims:ime,jms:jme,kms:kme+1,n_3d_u_h))
      this%a_3d_u_h = fill_value_
    end if
    if (n_3d_u > 0) then
      allocate(this%a_3d_u(ids:ide,jds:jde,kms:kme+1,n_3d_u))
      this%a_3d_u = fill_value_
    end if
    if (n_3d_uq_h > 0) then
      allocate(this%a_3d_uq_h(this%mesh%neqv,ims:ime,jms:jme,kms:kme+1,n_3d_uq_h))
      this%a_3d_uq_h = fill_value_
    end if
    if (n_3d_uq > 0) then
      allocate(this%a_3d_uq(this%mesh%neqv,ids:ide,jds:jde,kms:kme+1,n_3d_uq))
      this%a_3d_uq = fill_value_
    end if
    ! Down cell edge quadtrature
    if (n_3d_d_h > 0) then
      allocate(this%a_3d_d_h(ims:ime,jms:jme,kms:kme+1,n_3d_d_h))
      this%a_3d_d_h = fill_value_
    end if
    if (n_3d_d > 0) then
      allocate(this%a_3d_d(ids:ide,jds:jde,kms:kme+1,n_3d_d))
      this%a_3d_d = fill_value_
    end if
    if (n_3d_dq_h > 0) then
      allocate(this%a_3d_dq_h(this%mesh%neqv,ims:ime,jms:jme,kms:kme+1,n_3d_dq_h))
      this%a_3d_dq_h = fill_value_
    end if
    if (n_3d_dq > 0) then
      allocate(this%a_3d_dq(this%mesh%neqv,ids:ide,jds:jde,kms:kme+1,n_3d_dq))
      this%a_3d_dq = fill_value_
    end if

  end subroutine latlon_array_allocate_arrays

  subroutine latlon_array_get_array_2d_1(this, ptr, var_info)

    class(latlon_array_type), intent(in), target :: this
    real(r8), intent(out), pointer :: ptr(:,:)
    type(var_info_type), intent(in) :: var_info

    select case (var_info%loc)
    case ('C')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_c_h, 1):ubound(this%a_2d_c_h, 1), &
            lbound(this%a_2d_c_h, 2):ubound(this%a_2d_c_h, 2)) => this%a_2d_c_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_c  , 1):ubound(this%a_2d_c  , 1), &
            lbound(this%a_2d_c  , 2):ubound(this%a_2d_c  , 2)) => this%a_2d_c  (:,:,var_info%array_idx)
      end if
    case ('CA')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_ca_h, 1):ubound(this%a_2d_ca_h, 1), &
            lbound(this%a_2d_ca_h, 2):ubound(this%a_2d_ca_h, 2)) => this%a_2d_ca_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_ca  , 1):ubound(this%a_2d_ca  , 1), &
            lbound(this%a_2d_ca  , 2):ubound(this%a_2d_ca  , 2)) => this%a_2d_ca  (:,:,var_info%array_idx)
      end if
    case ('L')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_l_h, 1):ubound(this%a_2d_l_h, 1), &
            lbound(this%a_2d_l_h, 2):ubound(this%a_2d_l_h, 2)) => this%a_2d_l_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_l  , 1):ubound(this%a_2d_l  , 1), &
            lbound(this%a_2d_l  , 2):ubound(this%a_2d_l  , 2)) => this%a_2d_l  (:,:,var_info%array_idx)
      end if
    case ('R')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_r_h, 1):ubound(this%a_2d_r_h, 1), &
            lbound(this%a_2d_r_h, 2):ubound(this%a_2d_r_h, 2)) => this%a_2d_r_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_r  , 1):ubound(this%a_2d_r  , 1), &
            lbound(this%a_2d_r  , 2):ubound(this%a_2d_r  , 2)) => this%a_2d_r  (:,:,var_info%array_idx)
      end if
    case ('T')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_t_h, 1):ubound(this%a_2d_t_h, 1), &
            lbound(this%a_2d_t_h, 2):ubound(this%a_2d_t_h, 2)) => this%a_2d_t_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_t  , 1):ubound(this%a_2d_t  , 1), &
            lbound(this%a_2d_t  , 2):ubound(this%a_2d_t  , 2)) => this%a_2d_t  (:,:,var_info%array_idx)
      end if
    case ('B')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_b_h, 1):ubound(this%a_2d_b_h, 1), &
            lbound(this%a_2d_b_h, 2):ubound(this%a_2d_b_h, 2)) => this%a_2d_b_h(:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_b  , 1):ubound(this%a_2d_b  , 1), &
            lbound(this%a_2d_b  , 2):ubound(this%a_2d_b  , 2)) => this%a_2d_b  (:,:,var_info%array_idx)
      end if
    case default
      nullify(ptr)
    end select

  end subroutine latlon_array_get_array_2d_1

  subroutine latlon_array_get_array_2d_2(this, ptr, var_name, loc)

    class(latlon_array_type), intent(in) :: this
    real(r8), intent(out), pointer :: ptr(:,:)
    character(*), intent(in) :: var_name
    character(*), intent(in) :: loc

    integer i

    if (.not. this%initialized) then
      nullify(ptr)
      return
    end if
    i = this%var_idx(var_name, loc)
    call this%get_array(ptr, this%var_stack%var_info(i))
    if (.not. associated(ptr)) call log_error('Unable to get array "' // trim(var_name) // '"!', __FILE__, __LINE__, proc%id)

  end subroutine latlon_array_get_array_2d_2

  subroutine latlon_array_get_array_3d_1(this, ptr, var_info)

    class(latlon_array_type), intent(in), target :: this
    real(r8), intent(out), pointer :: ptr(:,:,:)
    type(var_info_type), intent(in) :: var_info

    select case (var_info%loc)
    case ('C')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_c_h, 1):ubound(this%a_3d_c_h, 1), &
            lbound(this%a_3d_c_h, 2):ubound(this%a_3d_c_h, 2), &
            lbound(this%a_3d_c_h, 3):ubound(this%a_3d_c_h, 3)) => this%a_3d_c_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_c  , 1):ubound(this%a_3d_c  , 1), &
            lbound(this%a_3d_c  , 2):ubound(this%a_3d_c  , 2), &
            lbound(this%a_3d_c  , 3):ubound(this%a_3d_c  , 3)) => this%a_3d_c  (:,:,:,var_info%array_idx)
      end if
    case ('CA')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_ca_h, 1):ubound(this%a_3d_ca_h, 1), &
            lbound(this%a_3d_ca_h, 2):ubound(this%a_3d_ca_h, 2), &
            lbound(this%a_3d_ca_h, 3):ubound(this%a_3d_ca_h, 3)) => this%a_3d_ca_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_ca  , 1):ubound(this%a_3d_ca  , 1), &
            lbound(this%a_3d_ca  , 2):ubound(this%a_3d_ca  , 2), &
            lbound(this%a_3d_ca  , 3):ubound(this%a_3d_ca  , 3)) => this%a_3d_ca  (:,:,:,var_info%array_idx)
      end if
    case ('CQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_cq_h, 1):ubound(this%a_2d_cq_h, 1), &
            lbound(this%a_2d_cq_h, 2):ubound(this%a_2d_cq_h, 2), &
            lbound(this%a_2d_cq_h, 3):ubound(this%a_2d_cq_h, 3)) => this%a_2d_cq_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_cq  , 1):ubound(this%a_2d_cq  , 1), &
            lbound(this%a_2d_cq  , 2):ubound(this%a_2d_cq  , 2), &
            lbound(this%a_2d_cq  , 3):ubound(this%a_2d_cq  , 3)) => this%a_2d_cq  (:,:,:,var_info%array_idx)
      end if
    case ('L')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_l_h, 1):ubound(this%a_3d_l_h, 1), &
            lbound(this%a_3d_l_h, 2):ubound(this%a_3d_l_h, 2), &
            lbound(this%a_3d_l_h, 3):ubound(this%a_3d_l_h, 3)) => this%a_3d_l_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_l  , 1):ubound(this%a_3d_l  , 1), &
            lbound(this%a_3d_l  , 2):ubound(this%a_3d_l  , 2), &
            lbound(this%a_3d_l  , 3):ubound(this%a_3d_l  , 3)) => this%a_3d_l  (:,:,:,var_info%array_idx)
      end if
    case ('LQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_lq_h, 1):ubound(this%a_2d_lq_h, 1), &
            lbound(this%a_2d_lq_h, 2):ubound(this%a_2d_lq_h, 2), &
            lbound(this%a_2d_lq_h, 3):ubound(this%a_2d_lq_h, 3)) => this%a_2d_lq_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_lq  , 1):ubound(this%a_2d_lq  , 1), &
            lbound(this%a_2d_lq  , 2):ubound(this%a_2d_lq  , 2), &
            lbound(this%a_2d_lq  , 3):ubound(this%a_2d_lq  , 3)) => this%a_2d_lq  (:,:,:,var_info%array_idx)
      end if
    case ('R')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_r_h, 1):ubound(this%a_3d_r_h, 1), &
            lbound(this%a_3d_r_h, 2):ubound(this%a_3d_r_h, 2), &
            lbound(this%a_3d_r_h, 3):ubound(this%a_3d_r_h, 3)) => this%a_3d_r_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_r  , 1):ubound(this%a_3d_r  , 1), &
            lbound(this%a_3d_r  , 2):ubound(this%a_3d_r  , 2), &
            lbound(this%a_3d_r  , 3):ubound(this%a_3d_r  , 3)) => this%a_3d_r  (:,:,:,var_info%array_idx)
      end if
    case ('RQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_rq_h, 1):ubound(this%a_2d_rq_h, 1), &
            lbound(this%a_2d_rq_h, 2):ubound(this%a_2d_rq_h, 2), &
            lbound(this%a_2d_rq_h, 3):ubound(this%a_2d_rq_h, 3)) => this%a_2d_rq_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_rq  , 1):ubound(this%a_2d_rq  , 1), &
            lbound(this%a_2d_rq  , 2):ubound(this%a_2d_rq  , 2), &
            lbound(this%a_2d_rq  , 3):ubound(this%a_2d_rq  , 3)) => this%a_2d_rq  (:,:,:,var_info%array_idx)
      end if
    case ('T')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_t_h, 1):ubound(this%a_3d_t_h, 1), &
            lbound(this%a_3d_t_h, 2):ubound(this%a_3d_t_h, 2), &
            lbound(this%a_3d_t_h, 3):ubound(this%a_3d_t_h, 3)) => this%a_3d_t_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_t  , 1):ubound(this%a_3d_t  , 1), &
            lbound(this%a_3d_t  , 2):ubound(this%a_3d_t  , 2), &
            lbound(this%a_3d_t  , 3):ubound(this%a_3d_t  , 3)) => this%a_3d_t  (:,:,:,var_info%array_idx)
      end if
    case ('TQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_tq_h, 1):ubound(this%a_2d_tq_h, 1), &
            lbound(this%a_2d_tq_h, 2):ubound(this%a_2d_tq_h, 2), &
            lbound(this%a_2d_tq_h, 3):ubound(this%a_2d_tq_h, 3)) => this%a_2d_tq_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_tq  , 1):ubound(this%a_2d_tq  , 1), &
            lbound(this%a_2d_tq  , 2):ubound(this%a_2d_tq  , 2), &
            lbound(this%a_2d_tq  , 3):ubound(this%a_2d_tq  , 3)) => this%a_2d_tq  (:,:,:,var_info%array_idx)
      end if
    case ('B')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_b_h, 1):ubound(this%a_3d_b_h, 1), &
            lbound(this%a_3d_b_h, 2):ubound(this%a_3d_b_h, 2), &
            lbound(this%a_3d_b_h, 3):ubound(this%a_3d_b_h, 3)) => this%a_3d_b_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_b  , 1):ubound(this%a_3d_b  , 1), &
            lbound(this%a_3d_b  , 2):ubound(this%a_3d_b  , 2), &
            lbound(this%a_3d_b  , 3):ubound(this%a_3d_b  , 3)) => this%a_3d_b  (:,:,:,var_info%array_idx)
      end if
    case ('BQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_2d_bq_h, 1):ubound(this%a_2d_bq_h, 1), &
            lbound(this%a_2d_bq_h, 2):ubound(this%a_2d_bq_h, 2), &
            lbound(this%a_2d_bq_h, 3):ubound(this%a_2d_bq_h, 3)) => this%a_2d_bq_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_2d_bq  , 1):ubound(this%a_2d_bq  , 1), &
            lbound(this%a_2d_bq  , 2):ubound(this%a_2d_bq  , 2), &
            lbound(this%a_2d_bq  , 3):ubound(this%a_2d_bq  , 3)) => this%a_2d_bq  (:,:,:,var_info%array_idx)
      end if
    case ('U')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_u_h, 1):ubound(this%a_3d_u_h, 1), &
            lbound(this%a_3d_u_h, 2):ubound(this%a_3d_u_h, 2), &
            lbound(this%a_3d_u_h, 3):ubound(this%a_3d_u_h, 3)) => this%a_3d_u_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_u  , 1):ubound(this%a_3d_u  , 1), &
            lbound(this%a_3d_u  , 2):ubound(this%a_3d_u  , 2), &
            lbound(this%a_3d_u  , 3):ubound(this%a_3d_u  , 3)) => this%a_3d_u  (:,:,:,var_info%array_idx)
      end if
    case ('D')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_d_h, 1):ubound(this%a_3d_d_h, 1), &
            lbound(this%a_3d_d_h, 2):ubound(this%a_3d_d_h, 2), &
            lbound(this%a_3d_d_h, 3):ubound(this%a_3d_d_h, 3)) => this%a_3d_d_h(:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_d  , 1):ubound(this%a_3d_d  , 1), &
            lbound(this%a_3d_d  , 2):ubound(this%a_3d_d  , 2), &
            lbound(this%a_3d_d  , 3):ubound(this%a_3d_d  , 3)) => this%a_3d_d  (:,:,:,var_info%array_idx)
      end if
    case default
      nullify(ptr)
    end select

  end subroutine latlon_array_get_array_3d_1

  subroutine latlon_array_get_array_3d_2(this, ptr, var_name, loc)

    class(latlon_array_type), intent(in) :: this
    real(r8), intent(out), pointer :: ptr(:,:,:)
    character(*), intent(in) :: var_name
    character(*), intent(in) :: loc

    integer i

    if (.not. this%initialized) then
      nullify(ptr)
      return
    end if
    i = this%var_idx(var_name, loc)
    if (i == 0) then
      nullify(ptr)
      return
    end if
    call this%get_array(ptr, this%var_stack%var_info(i))

  end subroutine latlon_array_get_array_3d_2

  subroutine latlon_array_get_array_4d_1(this, ptr, var_info)

    class(latlon_array_type), intent(in), target :: this
    real(r8), intent(out), pointer :: ptr(:,:,:,:)
    type(var_info_type), intent(in) :: var_info

    select case (var_info%loc)
    case ('CQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_cq_h, 1):ubound(this%a_3d_cq_h, 1), &
            lbound(this%a_3d_cq_h, 2):ubound(this%a_3d_cq_h, 2), &
            lbound(this%a_3d_cq_h, 3):ubound(this%a_3d_cq_h, 3), &
            lbound(this%a_3d_cq_h, 4):ubound(this%a_3d_cq_h, 4)) => this%a_3d_cq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_cq  , 1):ubound(this%a_3d_cq  , 1), &
            lbound(this%a_3d_cq  , 2):ubound(this%a_3d_cq  , 2), &
            lbound(this%a_3d_cq  , 3):ubound(this%a_3d_cq  , 3), &
            lbound(this%a_3d_cq  , 4):ubound(this%a_3d_cq  , 4)) => this%a_3d_cq  (:,:,:,:,var_info%array_idx)
      end if
    case ('LQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_lq_h, 1):ubound(this%a_3d_lq_h, 1), &
            lbound(this%a_3d_lq_h, 2):ubound(this%a_3d_lq_h, 2), &
            lbound(this%a_3d_lq_h, 3):ubound(this%a_3d_lq_h, 3), &
            lbound(this%a_3d_lq_h, 4):ubound(this%a_3d_lq_h, 4)) => this%a_3d_lq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_lq  , 1):ubound(this%a_3d_lq  , 1), &
            lbound(this%a_3d_lq  , 2):ubound(this%a_3d_lq  , 2), &
            lbound(this%a_3d_lq  , 3):ubound(this%a_3d_lq  , 3), &
            lbound(this%a_3d_lq  , 4):ubound(this%a_3d_lq  , 4)) => this%a_3d_lq  (:,:,:,:,var_info%array_idx)
      end if
    case ('RQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_rq_h, 1):ubound(this%a_3d_rq_h, 1), &
            lbound(this%a_3d_rq_h, 2):ubound(this%a_3d_rq_h, 2), &
            lbound(this%a_3d_rq_h, 3):ubound(this%a_3d_rq_h, 3), &
            lbound(this%a_3d_rq_h, 4):ubound(this%a_3d_rq_h, 4)) => this%a_3d_rq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_rq  , 1):ubound(this%a_3d_rq  , 1), &
            lbound(this%a_3d_rq  , 2):ubound(this%a_3d_rq  , 2), &
            lbound(this%a_3d_rq  , 3):ubound(this%a_3d_rq  , 3), &
            lbound(this%a_3d_rq  , 4):ubound(this%a_3d_rq  , 4)) => this%a_3d_rq  (:,:,:,:,var_info%array_idx)
      end if
    case ('TQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_tq_h, 1):ubound(this%a_3d_tq_h, 1), &
            lbound(this%a_3d_tq_h, 2):ubound(this%a_3d_tq_h, 2), &
            lbound(this%a_3d_tq_h, 3):ubound(this%a_3d_tq_h, 3), &
            lbound(this%a_3d_tq_h, 4):ubound(this%a_3d_tq_h, 4)) => this%a_3d_tq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_tq  , 1):ubound(this%a_3d_tq  , 1), &
            lbound(this%a_3d_tq  , 2):ubound(this%a_3d_tq  , 2), &
            lbound(this%a_3d_tq  , 3):ubound(this%a_3d_tq  , 3), &
            lbound(this%a_3d_tq  , 4):ubound(this%a_3d_tq  , 4)) => this%a_3d_tq  (:,:,:,:,var_info%array_idx)
      end if
    case ('BQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_bq_h, 1):ubound(this%a_3d_bq_h, 1), &
            lbound(this%a_3d_bq_h, 2):ubound(this%a_3d_bq_h, 2), &
            lbound(this%a_3d_bq_h, 3):ubound(this%a_3d_bq_h, 3), &
            lbound(this%a_3d_bq_h, 4):ubound(this%a_3d_bq_h, 4)) => this%a_3d_bq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_bq  , 1):ubound(this%a_3d_bq  , 1), &
            lbound(this%a_3d_bq  , 2):ubound(this%a_3d_bq  , 2), &
            lbound(this%a_3d_bq  , 3):ubound(this%a_3d_bq  , 3), &
            lbound(this%a_3d_bq  , 4):ubound(this%a_3d_bq  , 4)) => this%a_3d_bq  (:,:,:,:,var_info%array_idx)
      end if
    case ('UQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_uq_h, 1):ubound(this%a_3d_uq_h, 1), &
            lbound(this%a_3d_uq_h, 2):ubound(this%a_3d_uq_h, 2), &
            lbound(this%a_3d_uq_h, 3):ubound(this%a_3d_uq_h, 3), &
            lbound(this%a_3d_uq_h, 4):ubound(this%a_3d_uq_h, 4)) => this%a_3d_uq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_uq  , 1):ubound(this%a_3d_uq  , 1), &
            lbound(this%a_3d_uq  , 2):ubound(this%a_3d_uq  , 2), &
            lbound(this%a_3d_uq  , 3):ubound(this%a_3d_uq  , 3), &
            lbound(this%a_3d_uq  , 4):ubound(this%a_3d_uq  , 4)) => this%a_3d_uq  (:,:,:,:,var_info%array_idx)
      end if
    case ('DQ')
      if (var_info%with_halo) then
        ptr(lbound(this%a_3d_dq_h, 1):ubound(this%a_3d_dq_h, 1), &
            lbound(this%a_3d_dq_h, 2):ubound(this%a_3d_dq_h, 2), &
            lbound(this%a_3d_dq_h, 3):ubound(this%a_3d_dq_h, 3), &
            lbound(this%a_3d_dq_h, 4):ubound(this%a_3d_dq_h, 4)) => this%a_3d_dq_h(:,:,:,:,var_info%array_idx)
      else
        ptr(lbound(this%a_3d_dq  , 1):ubound(this%a_3d_dq  , 1), &
            lbound(this%a_3d_dq  , 2):ubound(this%a_3d_dq  , 2), &
            lbound(this%a_3d_dq  , 3):ubound(this%a_3d_dq  , 3), &
            lbound(this%a_3d_dq  , 4):ubound(this%a_3d_dq  , 4)) => this%a_3d_dq  (:,:,:,:,var_info%array_idx)
      end if
    case default
      nullify(ptr)
    end select

  end subroutine latlon_array_get_array_4d_1

  subroutine latlon_array_get_array_4d_2(this, ptr, loc, tag)

    class(latlon_array_type), intent(in), target :: this
    real(r8), intent(out), pointer :: ptr(:,:,:,:)
    character(*), intent(in) :: loc
    character(*), intent(in) :: tag

    integer ivs, ive, loc_is
    logical with_halo

    if (.not. this%initialized) then
      nullify(ptr)
      return
    end if

    ivs = 0; ive = 0; loc_is = 0
    do i = 1, this%var_stack%size
      associate (info => this%var_stack%var_info(i))
      if (info%loc == loc .and. info%tag == tag) then
        if (loc_is == 0) loc_is = i
        if (ivs == 0) then
          ivs = i
          with_halo = info%with_halo
        end if
      else if (ivs /= 0) then
        ive = i - 1
        exit
      end if
      end associate
    end do
    if (ivs == 0) return
    if (ive == 0) ive = i - 1 ! Matched variables are at tail.
    ivs = ivs - loc_is + 1
    ive = ive - loc_is + 1

    select case (loc)
    case ('C')
      if (with_halo) then
        ptr(lbound(this%a_3d_c_h, 1):ubound(this%a_3d_c_h, 1), &
            lbound(this%a_3d_c_h, 2):ubound(this%a_3d_c_h, 2), &
            lbound(this%a_3d_c_h, 3):ubound(this%a_3d_c_h, 3), &
            ivs:ive) => this%a_3d_c_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_c  , 1):ubound(this%a_3d_c  , 1), &
            lbound(this%a_3d_c  , 2):ubound(this%a_3d_c  , 2), &
            lbound(this%a_3d_c  , 3):ubound(this%a_3d_c  , 3), &
            ivs:ive) => this%a_3d_c  (:,:,:,ivs:ive)
      end if
    case ('CA')
      if (with_halo) then
        ptr(lbound(this%a_3d_ca_h, 1):ubound(this%a_3d_ca_h, 1), &
            lbound(this%a_3d_ca_h, 2):ubound(this%a_3d_ca_h, 2), &
            lbound(this%a_3d_ca_h, 3):ubound(this%a_3d_ca_h, 3), &
            ivs:ive) => this%a_3d_ca_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_ca  , 1):ubound(this%a_3d_ca  , 1), &
            lbound(this%a_3d_ca  , 2):ubound(this%a_3d_ca  , 2), &
            lbound(this%a_3d_ca  , 3):ubound(this%a_3d_ca  , 3), &
            ivs:ive) => this%a_3d_ca  (:,:,:,ivs:ive)
      end if
    case ('L')
      if (with_halo) then
        ptr(lbound(this%a_3d_l_h, 1):ubound(this%a_3d_l_h, 1), &
            lbound(this%a_3d_l_h, 2):ubound(this%a_3d_l_h, 2), &
            lbound(this%a_3d_l_h, 3):ubound(this%a_3d_l_h, 3), &
            ivs:ive) => this%a_3d_l_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_l  , 1):ubound(this%a_3d_l  , 1), &
            lbound(this%a_3d_l  , 2):ubound(this%a_3d_l  , 2), &
            lbound(this%a_3d_l  , 3):ubound(this%a_3d_l  , 3), &
            ivs:ive) => this%a_3d_l  (:,:,:,ivs:ive)
      end if
    case ('R')
      if (with_halo) then
        ptr(lbound(this%a_3d_r_h, 1):ubound(this%a_3d_r_h, 1), &
            lbound(this%a_3d_r_h, 2):ubound(this%a_3d_r_h, 2), &
            lbound(this%a_3d_r_h, 3):ubound(this%a_3d_r_h, 3), &
            ivs:ive) => this%a_3d_r_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_r  , 1):ubound(this%a_3d_r  , 1), &
            lbound(this%a_3d_r  , 2):ubound(this%a_3d_r  , 2), &
            lbound(this%a_3d_r  , 3):ubound(this%a_3d_r  , 3), &
            ivs:ive) => this%a_3d_r  (:,:,:,ivs:ive)
      end if
    case ('T')
      if (with_halo) then
        ptr(lbound(this%a_3d_t_h, 1):ubound(this%a_3d_t_h, 1), &
            lbound(this%a_3d_t_h, 2):ubound(this%a_3d_t_h, 2), &
            lbound(this%a_3d_t_h, 3):ubound(this%a_3d_t_h, 3), &
            ivs:ive) => this%a_3d_t_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_t  , 1):ubound(this%a_3d_t  , 1), &
            lbound(this%a_3d_t  , 2):ubound(this%a_3d_t  , 2), &
            lbound(this%a_3d_t  , 3):ubound(this%a_3d_t  , 3), &
            ivs:ive) => this%a_3d_t  (:,:,:,ivs:ive)
      end if
    case ('B')
      if (with_halo) then
        ptr(lbound(this%a_3d_b_h, 1):ubound(this%a_3d_b_h, 1), &
            lbound(this%a_3d_b_h, 2):ubound(this%a_3d_b_h, 2), &
            lbound(this%a_3d_b_h, 3):ubound(this%a_3d_b_h, 3), &
            ivs:ive) => this%a_3d_b_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_b  , 1):ubound(this%a_3d_b  , 1), &
            lbound(this%a_3d_b  , 2):ubound(this%a_3d_b  , 2), &
            lbound(this%a_3d_b  , 3):ubound(this%a_3d_b  , 3), &
            ivs:ive) => this%a_3d_b  (:,:,:,ivs:ive)
      end if
    case ('U')
      if (with_halo) then
        ptr(lbound(this%a_3d_u_h, 1):ubound(this%a_3d_u_h, 1), &
            lbound(this%a_3d_u_h, 2):ubound(this%a_3d_u_h, 2), &
            lbound(this%a_3d_u_h, 3):ubound(this%a_3d_u_h, 3), &
            ivs:ive) => this%a_3d_u_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_u  , 1):ubound(this%a_3d_u  , 1), &
            lbound(this%a_3d_u  , 2):ubound(this%a_3d_u  , 2), &
            lbound(this%a_3d_u  , 3):ubound(this%a_3d_u  , 3), &
            ivs:ive) => this%a_3d_u  (:,:,:,ivs:ive)
      end if
    case ('D')
      if (with_halo) then
        ptr(lbound(this%a_3d_d_h, 1):ubound(this%a_3d_d_h, 1), &
            lbound(this%a_3d_d_h, 2):ubound(this%a_3d_d_h, 2), &
            lbound(this%a_3d_d_h, 3):ubound(this%a_3d_d_h, 3), &
            ivs:ive) => this%a_3d_d_h(:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_d  , 1):ubound(this%a_3d_d  , 1), &
            lbound(this%a_3d_d  , 2):ubound(this%a_3d_d  , 2), &
            lbound(this%a_3d_d  , 3):ubound(this%a_3d_d  , 3), &
            ivs:ive) => this%a_3d_d  (:,:,:,ivs:ive)
      end if
    case default
      nullify(ptr)
    end select

  end subroutine latlon_array_get_array_4d_2

  subroutine latlon_array_get_array_5d(this, ptr, loc, tag)

    class(latlon_array_type), intent(in), target :: this
    real(r8), intent(out), pointer :: ptr(:,:,:,:,:)
    character(*), intent(in) :: loc
    character(*), intent(in) :: tag

    integer ivs, ive, loc_is
    logical with_halo

    if (.not. this%initialized) then
      nullify(ptr)
      return
    end if

    ivs = 0; ive = 0; loc_is = 0
    do i = 1, this%var_stack%size
      associate (info => this%var_stack%var_info(i))
      if (info%loc == loc .and. info%tag == tag) then
        if (loc_is == 0) loc_is = i
        if (ivs == 0) then
          ivs = i
          with_halo = info%with_halo
        end if
      else if (ivs /= 0) then
        ive = i - 1
        exit
      end if
      end associate
    end do
    if (ivs == 0) return
    if (ive == 0) ive = i - 1 ! Matched variables are at tail.
    ivs = ivs - loc_is + 1
    ive = ive - loc_is + 1

    select case (loc)
    case ('LQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_lq_h, 1):ubound(this%a_3d_lq_h, 1), &
            lbound(this%a_3d_lq_h, 2):ubound(this%a_3d_lq_h, 2), &
            lbound(this%a_3d_lq_h, 3):ubound(this%a_3d_lq_h, 3), &
            lbound(this%a_3d_lq_h, 4):ubound(this%a_3d_lq_h, 4), &
            ivs:ive) => this%a_3d_lq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_lq  , 1):ubound(this%a_3d_lq  , 1), &
            lbound(this%a_3d_lq  , 2):ubound(this%a_3d_lq  , 2), &
            lbound(this%a_3d_lq  , 3):ubound(this%a_3d_lq  , 3), &
            lbound(this%a_3d_lq  , 4):ubound(this%a_3d_lq  , 4), &
            ivs:ive) => this%a_3d_lq  (:,:,:,:,ivs:ive)
      end if
    case ('RQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_rq_h, 1):ubound(this%a_3d_rq_h, 1), &
            lbound(this%a_3d_rq_h, 2):ubound(this%a_3d_rq_h, 2), &
            lbound(this%a_3d_rq_h, 3):ubound(this%a_3d_rq_h, 3), &
            lbound(this%a_3d_rq_h, 4):ubound(this%a_3d_rq_h, 4), &
            ivs:ive) => this%a_3d_rq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_rq  , 1):ubound(this%a_3d_rq  , 1), &
            lbound(this%a_3d_rq  , 2):ubound(this%a_3d_rq  , 2), &
            lbound(this%a_3d_rq  , 3):ubound(this%a_3d_rq  , 3), &
            lbound(this%a_3d_rq  , 4):ubound(this%a_3d_rq  , 4), &
            ivs:ive) => this%a_3d_rq  (:,:,:,:,ivs:ive)
      end if
    case ('TQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_tq_h, 1):ubound(this%a_3d_tq_h, 1), &
            lbound(this%a_3d_tq_h, 2):ubound(this%a_3d_tq_h, 2), &
            lbound(this%a_3d_tq_h, 3):ubound(this%a_3d_tq_h, 3), &
            lbound(this%a_3d_tq_h, 4):ubound(this%a_3d_tq_h, 4), &
            ivs:ive) => this%a_3d_tq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_tq  , 1):ubound(this%a_3d_tq  , 1), &
            lbound(this%a_3d_tq  , 2):ubound(this%a_3d_tq  , 2), &
            lbound(this%a_3d_tq  , 3):ubound(this%a_3d_tq  , 3), &
            lbound(this%a_3d_tq  , 4):ubound(this%a_3d_tq  , 4), &
            ivs:ive) => this%a_3d_tq  (:,:,:,:,ivs:ive)
      end if
    case ('BQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_bq_h, 1):ubound(this%a_3d_bq_h, 1), &
            lbound(this%a_3d_bq_h, 2):ubound(this%a_3d_bq_h, 2), &
            lbound(this%a_3d_bq_h, 3):ubound(this%a_3d_bq_h, 3), &
            lbound(this%a_3d_bq_h, 4):ubound(this%a_3d_bq_h, 4), &
            ivs:ive) => this%a_3d_bq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_bq  , 1):ubound(this%a_3d_bq  , 1), &
            lbound(this%a_3d_bq  , 2):ubound(this%a_3d_bq  , 2), &
            lbound(this%a_3d_bq  , 3):ubound(this%a_3d_bq  , 3), &
            lbound(this%a_3d_bq  , 4):ubound(this%a_3d_bq  , 4), &
            ivs:ive) => this%a_3d_bq  (:,:,:,:,ivs:ive)
      end if
    case ('UQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_uq_h, 1):ubound(this%a_3d_uq_h, 1), &
            lbound(this%a_3d_uq_h, 2):ubound(this%a_3d_uq_h, 2), &
            lbound(this%a_3d_uq_h, 3):ubound(this%a_3d_uq_h, 3), &
            lbound(this%a_3d_uq_h, 4):ubound(this%a_3d_uq_h, 4), &
            ivs:ive) => this%a_3d_uq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_uq  , 1):ubound(this%a_3d_uq  , 1), &
            lbound(this%a_3d_uq  , 2):ubound(this%a_3d_uq  , 2), &
            lbound(this%a_3d_uq  , 3):ubound(this%a_3d_uq  , 3), &
            lbound(this%a_3d_uq  , 4):ubound(this%a_3d_uq  , 4), &
            ivs:ive) => this%a_3d_uq  (:,:,:,:,ivs:ive)
      end if
    case ('DQ')
      if (with_halo) then
        ptr(lbound(this%a_3d_dq_h, 1):ubound(this%a_3d_dq_h, 1), &
            lbound(this%a_3d_dq_h, 2):ubound(this%a_3d_dq_h, 2), &
            lbound(this%a_3d_dq_h, 3):ubound(this%a_3d_dq_h, 3), &
            lbound(this%a_3d_dq_h, 4):ubound(this%a_3d_dq_h, 4), &
            ivs:ive) => this%a_3d_dq_h(:,:,:,:,ivs:ive)
      else
        ptr(lbound(this%a_3d_dq  , 1):ubound(this%a_3d_dq  , 1), &
            lbound(this%a_3d_dq  , 2):ubound(this%a_3d_dq  , 2), &
            lbound(this%a_3d_dq  , 3):ubound(this%a_3d_dq  , 3), &
            lbound(this%a_3d_dq  , 4):ubound(this%a_3d_dq  , 4), &
            ivs:ive) => this%a_3d_dq  (:,:,:,:,ivs:ive)
      end if
    case default
      nullify(ptr)
    end select

  end subroutine latlon_array_get_array_5d

  subroutine latlon_array_create_dataset(this, dataset_name, desc, file_path, file_prefix)

    class(latlon_array_type), intent(in) :: this
    character(*), intent(in) :: dataset_name
    character(*), intent(in), optional :: desc
    character(*), intent(in), optional :: file_path
    character(*), intent(in), optional :: file_prefix

    integer nx, ny, nz, neq

    call this%mesh%get_params(nx=nx, ny=ny, nz=nz, neq=neq)

    call fiona_init()

    if (.not. fiona_has_dataset(dataset_name)) then
      call fiona_create_dataset(dataset_name, desc=desc, file_path=file_path, file_prefix=file_prefix, mpi_comm=proc%comm)
      call fiona_add_dim(dataset_name, 'x', 'Domain x coordinate', '', nx, decomp=.true.)
      call fiona_add_dim(dataset_name, 'y', 'Domain y coordinate', '', ny, decomp=.true.)
      call fiona_add_dim(dataset_name, 'z', 'Domain z coordinate', '', nz, decomp=.true.)
      call fiona_add_dim(dataset_name, 'time', add_var=.true.)
      call fiona_add_var(dataset_name, 'lon', 'Longitude', 'degrees_east' , ['x', 'y'], 'r4')
      call fiona_add_var(dataset_name, 'lat', 'Latitude' , 'degrees_north', ['x', 'y'], 'r4')
      call fiona_add_var(dataset_name,   'J', 'Jacobian sqrt(det(G))', '-', ['x', 'y'], 'r4')
    end if
    call this%append_dataset(dataset_name)

  end subroutine latlon_array_create_dataset

  subroutine latlon_array_append_dataset(this, dataset_name)

    class(latlon_array_type), intent(in) :: this
    character(*), intent(in) :: dataset_name

    character(4) cell_3d_dims(4), cell_2d_dims(3)
    integer i

    cell_3d_dims(1) = 'x'
    cell_3d_dims(2) = 'y'
    cell_3d_dims(3) = 'z'
    cell_3d_dims(4) = 'time'
    cell_2d_dims(1) = 'x'
    cell_2d_dims(2) = 'y'
    cell_2d_dims(3) = 'time'

    do i = 1, this%var_stack%size
      associate (info => this%var_stack%var_info(i))
      if (info%output) then
        if (fiona_has_var(dataset_name, info%name)) cycle
        if (proc%is_root()) call log_notice('Define output for variable ' // trim(info%name) // '.', pid=proc%id)
        if (info%only_2d) then
          select case (info%loc)
          case ('C', 'CA')
            call fiona_add_var(dataset_name, info%name, info%long_name, info%units, cell_2d_dims, 'r4')
            call fiona_add_att(dataset_name, info%name, 'coordinates', 'lon lat')
          end select
        else
          select case (info%loc)
          case ('C', 'CA')
            call fiona_add_var(dataset_name, info%name, info%long_name, info%units, cell_3d_dims, 'r4')
            call fiona_add_att(dataset_name, info%name, 'coordinates', 'lon lat z')
          end select
        end if
      end if
      end associate
    end do

  end subroutine latlon_array_append_dataset

  subroutine latlon_array_write(this, dataset_name, time_in_seconds, new_file, tag)

    class(latlon_array_type), intent(in), target :: this
    character(*), intent(in) :: dataset_name
    real(8), intent(in), optional :: time_in_seconds
    logical, intent(in), optional :: new_file
    character(*), intent(in), optional :: tag

    real(r8), pointer :: ptr_2d(:,:)
    real(r8), pointer :: ptr_3d(:,:,:)
    integer i, ids, ide, jds, jde, kds, kde, pc, nx, ny, nz

    call this%mesh%get_params(ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, pc=pc)

    nx = ide - ids + 1
    ny = jde - jds + 1
    nz = kde - kds + 1

    call fiona_start_output(dataset_name, time_in_seconds, new_file, tag)
    do i = 1, this%var_stack%size
      associate (info => this%var_stack%var_info(i))
      if (info%output) then
        if (info%only_2d) then
          call this%get_array(ptr_2d, info)
          call fiona_output(dataset_name, info%name, ptr_2d(ids:ide,jds:jde), &
                            start=[ids,jds], count=[nx,ny])
        else
          call this%get_array(ptr_3d, info)
          call fiona_output(dataset_name, info%name, ptr_3d(ids:ide,jds:jde,kds:kde), &
                            start=[ids,jds,kds], count=[nx,ny,nz])
        end if
      end if
      end associate
    end do
    call fiona_output(dataset_name, 'lon', this%mesh%lon(pc,ids:ide,jds:jde) * deg, start=[ids,jds], count=[nx,ny])
    call fiona_output(dataset_name, 'lat', this%mesh%lat(pc,ids:ide,jds:jde) * deg, start=[ids,jds], count=[nx,ny])
    call fiona_output(dataset_name,   'J', this%mesh% J (pc,ids:ide,jds:jde)      , start=[ids,jds], count=[nx,ny])
    call fiona_end_output(dataset_name)

  end subroutine latlon_array_write

  subroutine latlon_array_final(this)

    type(latlon_array_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_array_final

end module latlon_array_mod
