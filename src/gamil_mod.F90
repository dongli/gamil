module gamil_mod

  use flogger
  use string
  use latlon_process_mod
  use latlon_mesh_mod, mesh_type => latlon_mesh_type
  use gamil_params_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use dycore_mod
  use history_mod

  implicit none

  private

  public dycore

  public gamil_init
  public gamil_run
  public gamil_final

  type(mesh_type), allocatable :: mesh
  type(dycore_type), allocatable :: dycore

contains

  subroutine gamil_init(namelist)

    character(*), intent(in) :: namelist

    allocate(proc, mesh, dycore)

    call log_init()
    call fiona_init()
    call gamil_params_init(namelist)
    call time_init()
    call proc%init()
    call proc%decomp_domain(mesh, nx, ny, nz, hwx=2, hwy=2, nw=1, neq=1, r=radius)
    call dycore%init(mesh, model_type, recon_h_scheme, riemann_solver_type, time_scheme)
    call history_init(case_name, history_interval, dycore)

  end subroutine gamil_init

  subroutine gamil_run()

    call dycore%calc_swm_tend(old)
    call dycore%diag(old)
    call output(old)
    stop

    do while (.not. time_is_finished())
      call dycore%run(dt_in_seconds, old, new)
      call time_advance(dt_in_seconds)
      call dycore%diag(new)
      call output(new)
      stop 111
    end do

  end subroutine gamil_run

  subroutine output(itime)

    integer, intent(in) :: itime

    real(8), save :: time1 = 0, time2

    call log_print_diag(curr_time_str, pid=proc%id)
    if (.true. .or. time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.', pid=proc%id)
        time1 = time2
      end if
      call history_write(dycore, itime)
    end if

  end subroutine output

  subroutine gamil_final()

    deallocate(proc, mesh, dycore)

  end subroutine gamil_final

end module gamil_mod
