module dycore_mod

  use latlon_mesh_mod
  use static_mod
  use state_mod
  use swm_mod
  use dycore_interfaces_mod
  use fv_mod

  implicit none

  private

  public dycore_type
  public wave_speed_x
  public wave_speed_y
  public flux_x
  public flux_y
  public raw_to_cons

  type dycore_type
    logical :: initialized = .false.
    type(latlon_mesh_type), pointer :: mesh => null()
    type(static_type), allocatable :: static
    type(state_type), allocatable :: state(:)
  contains
    procedure :: init             => dycore_init
    procedure :: calc_contra_wind => dycore_calc_contra_wind
    procedure :: calc_sph_wind    => dycore_calc_sph_wind
    procedure :: calc_tend        => dycore_calc_tend
    procedure :: clear            => dycore_clear
    final :: dycore_final
  end type dycore_type

contains

  subroutine dycore_init(this, mesh, model_type)

    class(dycore_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: model_type

    integer i

    call this%clear()

    this%mesh => mesh

    allocate(this%static)
    call this%static%init(mesh)

    allocate(this%state(1:1))
    do i = 1, size(this%state)
      call this%state(i)%init(mesh, model_type)
    end do

    select case (model_type)
    case ('swm')
      wave_speed_x => swm_wave_speed_x
      wave_speed_y => swm_wave_speed_y
      flux_x       => swm_flux_x
      flux_y       => swm_flux_y
      sources      => swm_sources
      raw_to_cons  => swm_raw_to_cons
      before_recon => swm_before_recon
      after_flux   => swm_after_flux
    end select

    this%initialized = .true.

  end subroutine dycore_init

  subroutine dycore_calc_contra_wind(this, itime)

    class(dycore_type), intent(in) :: this
    integer, intent(in) :: itime

    integer i, j, k

    associate (mesh => this%mesh, state => this%state(itime))
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          state%uc(i,j,k) = state%u(i,j,k) * sqrt(mesh%iG(1,1,1,i,j))
          state%vc(i,j,k) = state%v(i,j,k) * sqrt(mesh%iG(2,2,1,i,j))
        end do
      end do
    end do
    end associate

  end subroutine dycore_calc_contra_wind

  subroutine dycore_calc_sph_wind(this, itime)

    class(dycore_type), intent(in) :: this
    integer, intent(in) :: itime

    integer i, j, k

    associate (mesh => this%mesh, state => this%state(itime))
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          state%u(i,j,k) = state%uc(i,j,k) / sqrt(mesh%iG(1,1,1,i,j))
          state%v(i,j,k) = state%vc(i,j,k) / sqrt(mesh%iG(2,2,1,i,j))
        end do
      end do
    end do
    end associate

  end subroutine dycore_calc_sph_wind

  subroutine dycore_calc_tend(this, itime)

    class(dycore_type), intent(in) :: this
    integer, intent(in) :: itime

    integer i, j, k

    associate (mesh => this%mesh, state => this%state(itime))
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
        end do
      end do
    end do
    end associate

  end subroutine dycore_calc_tend

  subroutine dycore_clear(this)

    class(dycore_type), intent(inout) :: this

    nullify(this%mesh)

    if (allocated(this%static)) deallocate(this%static)
    if (allocated(this%state )) deallocate(this%state )

    nullify(wave_speed_x)
    nullify(wave_speed_y)

    this%initialized = .false.

  end subroutine dycore_clear

  subroutine dycore_final(this)

    type(dycore_type), intent(inout) :: this

    call this%clear()

  end subroutine dycore_final

end module dycore_mod
