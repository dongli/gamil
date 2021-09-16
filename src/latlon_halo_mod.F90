module latlon_halo_mod

  use mpi
  use latlon_mesh_mod

  implicit none

  type latlon_halo_type
    logical :: initialized = .false.
    integer :: host_id = MPI_PROC_NULL
    integer :: proc_id = MPI_PROC_NULL
    integer :: orient  = 0
    integer :: dtype   = 0
    integer :: send_type_2d = MPI_DATATYPE_NULL
    integer :: recv_type_2d = MPI_DATATYPE_NULL
    integer :: send_type_3d = MPI_DATATYPE_NULL
    integer :: recv_type_3d = MPI_DATATYPE_NULL
    integer :: send_array_size    (3) =  0
    integer :: send_subarray_start(3) = -1
    integer :: send_subarray_size (3) =  0
    integer :: recv_array_size    (3) =  0
    integer :: recv_subarray_start(3) = -1
    integer :: recv_subarray_size (3) =  0
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

    integer ims, ime, ids, ide, jms, jme, jds, jde, kms, kme, hwx, hwy, nx, ny, nz
    integer ierr

    call this%clear()

    if (proc_id == MPI_PROC_NULL) return

    this%orient  = orient
    this%dtype   = dtype
    this%host_id = host_id
    this%proc_id = proc_id

    call mesh%get_params(ims=ims, ime=ime, ids=ids, ide=ide, &
                         jms=jms, jme=jme, jds=jds, jde=jde, &
                         kms=kms, kme=kme, hwx=hwx, hwy=hwy)

    nx = ide - ids + 1
    ny = jde - jds + 1
    nz = kme - kms + 1

    this%send_array_size = [ime-ims+1,jme-jms+1,kme-kms+1]
    this%recv_array_size = this%send_array_size

    ! NOTE: MPI array index starts from zero.
    select case (orient)
    case (left)
      this%send_subarray_start = [hwx         ,hwy          ,0 ]
      this%recv_subarray_start = [0           ,hwy          ,0 ]
      this%send_subarray_size  = [hwx         ,ny           ,nz]
      this%recv_subarray_size  = [hwx          ,ny           ,nz]
    case (left_bottom)
      this%send_subarray_start = [hwx          ,hwy          ,0 ]
      this%recv_subarray_start = [0            ,0            ,0 ]
      this%send_subarray_size  = [hwx          ,hwy          ,nz]
      this%recv_subarray_size  = [hwx          ,hwy          ,nz]
    case (left_top)
      this%send_subarray_start = [hwx          ,jde-hwy+1-jms,0 ]
      this%recv_subarray_start = [0            ,jde    +1-jms,0 ]
      this%send_subarray_size  = [hwx          ,hwy          ,nz]
      this%recv_subarray_size  = [hwx          ,hwy          ,nz]
    case (right)
      this%send_subarray_start = [ide-hwx+1-ims,hwy          ,0 ]
      this%recv_subarray_start = [ide    +1-ims,hwy          ,0 ]
      this%send_subarray_size  = [hwx          ,ny           ,nz]
      this%recv_subarray_size  = [hwx          ,ny           ,nz]
    case (right_bottom)
      this%send_subarray_start = [ide-hwx+1-ims,hwy          ,0 ]
      this%recv_subarray_start = [ide    +1-ims,0            ,0 ]
      this%send_subarray_size  = [hwx          ,hwy          ,nz]
      this%recv_subarray_size  = [hwx          ,hwy          ,nz]
    case (right_top)
      this%send_subarray_start = [ide-hwx+1-ims,jde-hwy+1-jms,0 ]
      this%recv_subarray_start = [ide    +1-ims,jde    +1-jms,0 ]
      this%send_subarray_size  = [hwx          ,hwy          ,nz]
      this%recv_subarray_size  = [hwx          ,hwy          ,nz]
    case (top)
      this%send_subarray_start = [hwx          ,jde-hwy+1-jms,0 ]
      this%recv_subarray_start = [hwx          ,jde    +1-jms,0 ]
      this%send_subarray_size  = [nx           ,hwy          ,nz]
      this%recv_subarray_size  = [nx           ,hwy          ,nz]
    case (bottom)
      this%send_subarray_start = [hwx          ,hwy          ,0 ]
      this%recv_subarray_start = [hwx          ,0            ,0 ]
      this%send_subarray_size  = [nx           ,hwy          ,nz]
      this%recv_subarray_size  = [nx           ,hwy          ,nz]
    end select

    call MPI_TYPE_CREATE_SUBARRAY(3, this%send_array_size, this%send_subarray_size, &
                                  this%send_subarray_start, MPI_ORDER_FORTRAN, dtype, &
                                  this%send_type_3d, ierr)
    call MPI_TYPE_COMMIT(this%send_type_3d, ierr)
    call MPI_TYPE_CREATE_SUBARRAY(3, this%recv_array_size, this%recv_subarray_size, &
                                  this%recv_subarray_start, MPI_ORDER_FORTRAN, dtype, &
                                  this%recv_type_3d, ierr)
    call MPI_TYPE_COMMIT(this%recv_type_3d, ierr)

    call MPI_TYPE_CREATE_SUBARRAY(2, this%send_array_size, this%send_subarray_size, &
                                  this%send_subarray_start, MPI_ORDER_FORTRAN, dtype, &
                                  this%send_type_2d, ierr)
    call MPI_TYPE_COMMIT(this%send_type_2d, ierr)
    call MPI_TYPE_CREATE_SUBARRAY(2, this%recv_array_size, this%recv_subarray_size, &
                                  this%recv_subarray_start, MPI_ORDER_FORTRAN, dtype, &
                                  this%recv_type_2d, ierr)
    call MPI_TYPE_COMMIT(this%recv_type_2d, ierr)

    this%initialized = .true.

  end subroutine latlon_halo_init

  subroutine latlon_halo_clear(this)

    class(latlon_halo_type), intent(inout) :: this

    integer ierr

    if (this%send_type_2d /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_2d, ierr)
    if (this%recv_type_2d /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_2d, ierr)
    if (this%send_type_3d /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_3d, ierr)
    if (this%recv_type_3d /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_3d, ierr)

    this%host_id             = MPI_PROC_NULL
    this%proc_id             = MPI_PROC_NULL
    this%orient              =  0
    this%dtype               =  0
    this%send_array_size     =  0
    this%send_subarray_start = -1
    this%send_subarray_size  =  0
    this%recv_array_size     =  0
    this%recv_subarray_start = -1
    this%recv_subarray_size  =  0

    this%initialized = .false.

  end subroutine latlon_halo_clear

  subroutine latlon_halo_final(this)

    type(latlon_halo_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_halo_final

end module latlon_halo_mod
