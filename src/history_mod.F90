module history_mod

  use flogger
  use string
  use time_mod
  use dycore_mod
  use latlon_process_mod

  implicit none

  private

  public history_init
  public history_write

contains

  subroutine history_init(case_name, history_interval, dycore)

    character(*), intent(in) :: case_name
    character(*), intent(in) :: history_interval(:)
    type(dycore_type), intent(in) :: dycore

    character(10) time_value, time_units
    real(8) seconds, months

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!', pid=proc%id)
    end select

    call time_add_alert('history_write', seconds=seconds)

    call dycore%static%array%create_dataset('h0', file_prefix=case_name)
    call dycore%state(1)%array%append_dataset('h0')

  end subroutine history_init

  subroutine history_write(dycore, itime)

    type(dycore_type), intent(in) :: dycore
    integer, intent(in) :: itime

    call dycore%static%array%write('h0', time_in_seconds=0.0d0, new_file=time_step==0)
    call dycore%state(itime)%array%write('h0', time_in_seconds=0.0d0, new_file=.false.)

    call log_notice('Write history file.', pid=proc%id)

  end subroutine history_write

end module history_mod
