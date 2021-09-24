module latlon_parallel_mod

  use mpi
  use latlon_halo_mod
  use latlon_array_mod

  implicit none

  private

  public fill_halo
  public global_min
  public global_max
  public global_sum

  interface global_min
    module procedure global_min_3d_r4
    module procedure global_min_3d_r8
  end interface global_min

  interface global_max
    module procedure global_max_3d_r4
    module procedure global_max_3d_r8
  end interface global_max

  interface global_sum
    module procedure global_sum_0d_r8
  end interface global_sum

contains

  subroutine fill_halo(array)

    type(latlon_array_type), intent(inout) :: array

    integer io, iv

    do io = 1, size(array%halo)
      do iv = 1, array%var_stack%size
        call send_halo(array, array%var_stack%var_info(iv), array%halo(io))
      end do
    end do
    ! Recv data.
    do io = 1, size(array%halo)
      do iv = 1, array%var_stack%size
        call recv_halo(array, array%var_stack%var_info(iv), array%halo(io))
      end do
    end do

  end subroutine fill_halo

  subroutine send_halo(array, var_info, halo)

    type(latlon_array_type), intent(inout) :: array
    type(var_info_type), intent(in) :: var_info
    type(latlon_halo_type), intent(in) :: halo

    integer i, tag, ierr

    if (.not. var_info%fill_halo .or. .not. halo%initialized) return

    i = var_info%array_idx
    tag = 10 * i + halo%orient

    select case (var_info%loc)
    case ('C')
      if (var_info%only_2d) then
        call MPI_SEND(array%a_2d_c_h(:,:,i), i, halo%send_type_2d, halo%proc_id, tag, &
                      proc%comm, ierr)
      else
        call MPI_SEND(array%a_3d_c_h(:,:,:,i), i, halo%send_type_3d, halo%proc_id, tag, &
                      proc%comm, ierr)
      end if
    case ('CA')
      if (var_info%only_2d) then
        call MPI_SEND(array%a_2d_ca_h(:,:,i), i, halo%send_type_2d, halo%proc_id, tag, &
                      proc%comm, ierr)
      else
        call MPI_SEND(array%a_3d_ca_h(:,:,:,i), i, halo%send_type_3d, halo%proc_id, tag, &
                      proc%comm, ierr)
      end if
    case default
      stop 'Not support halo filling except for cell arrays!'
    end select

  end subroutine send_halo

  subroutine recv_halo(array, var_info, halo)

    type(latlon_array_type), intent(inout), target :: array
    type(var_info_type), intent(in) :: var_info
    type(latlon_halo_type), intent(in) :: halo

    integer i, tag, ierr

    if (.not. var_info%fill_halo .or. .not. halo%initialized) return

    i = var_info%array_idx
    tag = 10 * i + halo%pair_orient

    select case (var_info%loc)
    case ('C')
      if (var_info%only_2d) then
        call MPI_RECV(array%a_2d_c_h(:,:,i), i, halo%recv_type_2d, halo%proc_id, tag, &
                      proc%comm, MPI_STATUS_IGNORE, ierr)
      else
        call MPI_RECV(array%a_3d_c_h(:,:,:,i), i, halo%recv_type_3d, halo%proc_id, tag, &
                      proc%comm, MPI_STATUS_IGNORE, ierr)
      end if
    case ('CA')
      if (var_info%only_2d) then
        call MPI_RECV(array%a_2d_ca_h(:,:,i), i, halo%recv_type_2d, halo%proc_id, tag, &
                      proc%comm, MPI_STATUS_IGNORE, ierr)
      else
        call MPI_RECV(array%a_3d_ca_h(:,:,:,i), i, halo%recv_type_3d, halo%proc_id, tag, &
                      proc%comm, MPI_STATUS_IGNORE, ierr)
      end if
    case default
      stop 'Not support halo filling except for cell arrays!'
    end select

  end subroutine recv_halo

  real(4) function global_min_3d_r4(x) result(res)

    real(4), intent(in) :: x(:,:,:)

    integer ierr

    call MPI_ALLREDUCE(minval(x), res, 1, MPI_REAL, MPI_MIN, proc%comm, ierr)

  end function global_min_3d_r4

  real(8) function global_min_3d_r8(x) result(res)

    real(8), intent(in) :: x(:,:,:)

    integer ierr

    call MPI_ALLREDUCE(minval(x), res, 1, MPI_DOUBLE, MPI_MIN, proc%comm, ierr)

  end function global_min_3d_r8

  real(4) function global_max_3d_r4(x) result(res)

    real(4), intent(in) :: x(:,:,:)

    integer ierr

    call MPI_ALLREDUCE(maxval(x), res, 1, MPI_REAL, MPI_MAX, proc%comm, ierr)

  end function global_max_3d_r4

  real(8) function global_max_3d_r8(x) result(res)

    real(8), intent(in) :: x(:,:,:)

    integer ierr

    call MPI_ALLREDUCE(maxval(x), res, 1, MPI_DOUBLE, MPI_MAX, proc%comm, ierr)

  end function global_max_3d_r8

  real(8) function global_sum_0d_r8(x) result(res)

    real(8), intent(in) :: x

    integer ierr

    call MPI_ALLREDUCE(x, res, 1, MPI_DOUBLE, MPI_SUM, proc%comm, ierr)

  end function global_sum_0d_r8

end module latlon_parallel_mod
