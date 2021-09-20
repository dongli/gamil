module time_schemes_mod

  use kinds_mod
  use dycore_mod

  implicit none

  private

  public tvd3

contains

  subroutine tvd3(dt, old, new, dycore)

    real(r8), intent(in) :: dt
    integer , intent(in) :: old
    integer , intent(in) :: new
    type(dycore_type), intent(inout) :: dycore

    call dycore%calc_swm_tend(old)
    call dycore%update_state(dt, old, new)



  end subroutine tvd3

end module time_schemes_mod
