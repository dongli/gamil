module history_mod

  use dycore_mod

  implicit none

  private

  public history_init
  public history_write

contains

  subroutine history_init(case_name, dycore)

    character(*), intent(in) :: case_name
    type(dycore_type), intent(in) :: dycore

    call dycore%static%array%create_dataset('h0', file_prefix=case_name)
    call dycore%state(1)%array%append_dataset('h0')

  end subroutine history_init

  subroutine history_write(dycore, itime)

    type(dycore_type), intent(in) :: dycore
    integer, intent(in) :: itime

    call dycore%static%array%write('h0', time_in_seconds=0.0d0, new_file=.true.)
    call dycore%state(itime)%array%write('h0', time_in_seconds=0.0d0, new_file=.false.)

  end subroutine history_write

end module history_mod
