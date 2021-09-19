module dycore_mod

  use latlon_mesh_mod
  use static_mod
  use state_mod
  use swm_mod
  use fv_mod

  implicit none

  private

  public dycore_type

  type dycore_type
    logical :: initialized = .false.
    type(latlon_mesh_type), pointer :: mesh => null()
    type(static_type), allocatable :: static
    type(state_type), allocatable :: state(:)
  contains
    procedure :: init                    => dycore_init
    procedure :: calc_contravariant_wind => dycore_calc_contravariant_wind
    procedure :: calc_spherical_wind     => dycore_calc_spherical_wind
    procedure :: calc_swm_tendencies     => dycore_calc_swm_tendencies
    procedure :: clear                   => dycore_clear
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

    this%initialized = .true.

  end subroutine dycore_init

  subroutine dycore_calc_contravariant_wind(this, itime)

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

  end subroutine dycore_calc_contravariant_wind

  subroutine dycore_calc_spherical_wind(this, itime)

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

  end subroutine dycore_calc_spherical_wind

  subroutine dycore_calc_swm_tendencies(this, itime)

    class(dycore_type), intent(inout) :: this
    integer, intent(in) :: itime

    real(r8) flux(this%state(itime)%nvar)
    real(r8) ql  (this%state(itime)%nvar)
    real(r8) qr  (this%state(itime)%nvar)
    integer i, j, k, e

    associate (mesh => this%mesh, state => this%state(itime))
    call swm_before_reconstruct(state, this%static)
    call reconstruct(state)
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          do e = mesh%pes(left), mesh%pee(left)
            ql = state%qr(i-1,j,k,:)
            qr = state%ql(i  ,j,k,:)
            flux = riemann_solver(mesh%iG(:,:,e,i,j), mesh%J(e,i,j), ql, qr, swm_wave_speed_x, swm_flux_x)
          end do
        end do
      end do
    end do
    end associate

  end subroutine dycore_calc_swm_tendencies

  subroutine dycore_clear(this)

    class(dycore_type), intent(inout) :: this

    nullify(this%mesh)

    if (allocated(this%static)) deallocate(this%static)
    if (allocated(this%state )) deallocate(this%state )

    this%initialized = .false.

  end subroutine dycore_clear

  subroutine dycore_final(this)

    type(dycore_type), intent(inout) :: this

    call this%clear()

  end subroutine dycore_final

end module dycore_mod
