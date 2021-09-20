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
    procedure :: init                => dycore_init
    procedure :: calc_contravar_wind => dycore_calc_contravar_wind
    procedure :: calc_spherical_wind => dycore_calc_spherical_wind
    procedure :: calc_swm_tend       => dycore_calc_swm_tend
    procedure :: update_state        => dycore_update_state
    procedure :: clear               => dycore_clear
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

  subroutine dycore_calc_contravar_wind(this, itime)

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

  end subroutine dycore_calc_contravar_wind

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

  subroutine dycore_calc_swm_tend(this, itime)

    class(dycore_type), intent(inout) :: this
    integer, intent(in) :: itime

    real(r8) res(this%state(itime)%nvar)
    real(r8) q  (this%state(itime)%nvar)
    real(r8) ql (this%state(itime)%nvar)
    real(r8) qr (this%state(itime)%nvar)
    integer i, j, k, l, e

    k = 1 ! Shallow-water model only has one level.
    associate (mesh => this%mesh, static => this%static, state => this%state(itime))
    ! Change depth to height for reconstruction.
    do j = mesh%jms, mesh%jme
      do i = mesh%ims, mesh%ime
        state%q(i,j,k,1) = state%q(i,j,k,1) + mesh%J(1,i,j) * static%zs(i,j)
      end do
    end do
    ! Calculate height gradient at cell center.
    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        state%dhdx(i,j,k) = (state%h(i-2,j,k) - 8 * state%h(i-1,j,k) + 8 * state%h(i+1,j,k) - state%h(i+2,j,k)) / (12 * mesh%dx)
        state%dhdy(i,j,k) = (state%h(i,j-2,k) - 8 * state%h(i,j-1,k) + 8 * state%h(i,j+1,k) - state%h(i,j+2,k)) / (12 * mesh%dy)
      end do
    end do
    ! Handle poles.
    if (mesh%jds == 1) then
      j = 1
      do i = mesh%ids, mesh%ide
        state%dhdy(i,j,k) = (4 * state%h(i,j+1,k) - state%h(i,j+2,k) - 3 * state%h(i,j,k)) / (2 * mesh%dy)
      end do
    end if
    if (mesh%jde == mesh%ny) then
      j = mesh%ny
      do i = mesh%ids, mesh%ide
        state%dhdy(i,j,k) = (3 * state%h(i,j,k) + state%h(i,j-2,k) - 4 * state%h(i,j-1,k)) / (2 * mesh%dy)
      end do
    end if
    ! Reconstruct forecast variables to cell edges.
    do l = 1, state%nvar
      do j = mesh%jds - 1, mesh%jde + 1
        do i = mesh%ids - 1, mesh%ide + 1
          ! FIXME: Here we assume there is only one edge quadrature point.
          state%ql(1,i,j,k,l) = reconstruct(-1, state%q(i-2:i+2,j,k,l))
          state%qr(1,i,j,k,l) = reconstruct( 1, state%q(i-2:i+2,j,k,l))
          state%qb(1,i,j,k,l) = reconstruct(-1, state%q(i,j-2:j+2,k,l))
          state%qt(1,i,j,k,l) = reconstruct( 1, state%q(i,j-2:j+2,k,l))
        end do
      end do
    end do
    ! Calculate numerical flux at cell edges using Riemann solver.
    do j = mesh%jds - 1, mesh%jde + 1
      do i = mesh%ids - 1, mesh%ide + 1
        do e = mesh%pes(left), mesh%pee(left)
          ql  = state%qr(e,i-1,j,k,:)
          qr  = state%ql(e,i  ,j,k,:)
          res = riemann_solver(mesh%iG(:,:,e,i,j), mesh%J(e,i,j), ql, qr, &
                               swm_wave_speed_x, swm_flux_x)
          state%fx(e,i,j,k,:) = res
        end do
        do e = mesh%pes(bottom), mesh%pee(bottom)
          ql  = state%qt(e,i,j-1,k,:)
          qr  = state%qb(e,i,j  ,k,:)
          res = riemann_solver(mesh%iG(:,:,e,i,j), mesh%J(e,i,j), ql, qr, &
                               swm_wave_speed_y, swm_flux_y)
          state%fy(e,i,j,k,:) = res
        end do
      end do
    end do
    ! Change height back to depth.
    do j = mesh%jms, mesh%jme
      do i = mesh%ims, mesh%ime
        state%q(i,j,k,1) = state%q(i,j,k,1) - mesh%J(1,i,j) * static%zs(i,j)
      end do
    end do
    ! Calculate sources and set tendencies.
    do j = mesh%jds, mesh%jde + 1
      do i = mesh%ids, mesh%ide + 1
        q   = state%q(i,j,k,:)
        res = swm_sources(mesh%iG(:,:,1,i,j), mesh%J(1,i,j), mesh%CS(:,:,:,i,j), static%f1(i,j), &
                          q, static%zs(i,j), state%dhdx(i,j,k), state%dhdy(i,j,k))
        state%dqdt(i,j,k,:) = - (state%fx(1,i+1,j,k,:) - state%fx(1,i,j,k,:)) / mesh%dx &
                              - (state%fy(1,i,j+1,k,:) - state%fy(1,i,j,k,:)) / mesh%dy &
                              + res
      end do
    end do
    end associate

  end subroutine dycore_calc_swm_tend

  subroutine dycore_update_state(this, dt, old, new)

    class(dycore_type), intent(inout) :: this
    real(r8), intent(in) :: dt
    integer , intent(in) :: old
    integer , intent(in) :: new

    integer i, j, k

    associate (mesh => this%mesh, old_state => this%state(old), new_state => this%state(new))
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          new_state%q(i,j,k,:) = old_state%q(i,j,k,:) + dt * old_state%dqdt(i,j,k,:)
        end do
      end do
    end do
    end associate

  end subroutine dycore_update_state

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
