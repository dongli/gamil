module dycore_mod

  use latlon_mesh_mod
  use latlon_parallel_mod
  use static_mod
  use state_mod
  use swm_mod
  use fv_mod

  use weno_mod

  implicit none

  private

  public dycore_type

  type dycore_type
    logical :: initialized = .false.
    character(10) model_type
    type(latlon_mesh_type), pointer :: mesh => null()
    type(static_type), allocatable :: static
    type(state_type), allocatable :: state(:)
  contains
    procedure :: init                => dycore_init
    procedure :: run                 => dycore_run
    procedure :: diag                => dycore_diag
    procedure :: calc_spherical_wind => dycore_calc_spherical_wind
    procedure :: calc_swm_tend       => dycore_calc_swm_tend
    procedure :: update_state        => dycore_update_state
    procedure :: clear               => dycore_clear
    final :: dycore_final
  end type dycore_type

  interface
    subroutine time_integrate_interface(dt, old, new, dycore)
      import r8, dycore_type
      real(r8), intent(in) :: dt
      integer , intent(in) :: old
      integer , intent(in) :: new
      type(dycore_type), intent(inout) :: dycore
    end subroutine time_integrate_interface
  end interface

  procedure(time_integrate_interface), pointer :: time_integrate

contains

  subroutine dycore_init(this, mesh, model_type, recon_h_scheme, riemann_solver_type, time_scheme)

    class(dycore_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: model_type
    character(*), intent(in) :: recon_h_scheme
    character(*), intent(in) :: riemann_solver_type
    character(*), intent(in) :: time_scheme

    integer i

    call this%clear()

    call fv_init(recon_h_scheme, riemann_solver_type)

    this%model_type = model_type
    this%mesh => mesh

    allocate(this%static)
    call this%static%init(mesh)

    select case (time_scheme)
    case ('tvd3')
      allocate(this%state(-1:2))
      time_integrate => tvd3
    end select
    do i = lbound(this%state, 1), ubound(this%state, 1)
      call this%state(i)%init(mesh, model_type)
    end do

    this%initialized = .true.

  end subroutine dycore_init

  subroutine dycore_run(this, dt, old, new)

    class(dycore_type), intent(inout) :: this
    real(r8), intent(in) :: dt
    integer , intent(in) :: old
    integer , intent(in) :: new

    call time_integrate(dt, old, new, this)

  end subroutine dycore_run

  subroutine dycore_diag(this, itime)

    class(dycore_type), intent(in) :: this
    integer, intent(in) :: itime

    real(r8) total_mass
    integer i, j, k

    total_mass = 0
    associate (mesh => this%mesh, state => this%state(itime))
    select case (this%model_type)
    case ('swm')
      k = 1
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          total_mass = total_mass + state%h(i,j,k) - this%static%zs(i,j)
        end do
      end do
    end select
    end associate

    total_mass = global_sum(total_mass)

    call log_add_diag('tm', total_mass)

  end subroutine dycore_diag

  subroutine dycore_calc_spherical_wind(this, itime)

    class(dycore_type), intent(in) :: this
    integer, intent(in) :: itime

    integer i, j, k

    associate (mesh => this%mesh, state => this%state(itime))
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          state%u(i,j,k) = state%q(i,j,k,2) / mesh%J(1,i,j) / sqrt(mesh%iG(1,1,1,i,j))
          state%v(i,j,k) = state%q(i,j,k,3) / mesh%J(1,i,j) / sqrt(mesh%iG(2,2,1,i,j))
        end do
      end do
    end do
    end associate

  end subroutine dycore_calc_spherical_wind

  subroutine dycore_calc_swm_tend(this, itime)

    class(dycore_type), intent(inout) :: this
    integer, intent(in) :: itime

    real(r8) q (this%state(itime)%nvar)
    real(r8) ql(this%state(itime)%nvar)
    real(r8) qr(this%state(itime)%nvar)
    integer i, j, k, l, pc

    k  = 1 ! Shallow-water model only has one level.
    pc = 1 ! Cell center point index

    associate (mesh => this%mesh, static => this%static, state => this%state(itime))
    ! Change depth to height for reconstruction.
    do j = mesh%jms, mesh%jme
      do i = mesh%ims, mesh%ime
        state%q(i,j,k,1) = state%q(i,j,k,1) + mesh%J(pc,i,j) * static%zs(i,j)
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
        ql  = state%qr(1,i-1,j,k,:)
        qr  = state%ql(1,i  ,j,k,:)
        state%fx(1,i,j,k,:) = riemann_solver(mesh%iG(:,:,mesh%pes(left),i,j), &
                                             mesh% J(    mesh%pes(left),i,j), &
                                             ql, qr, swm_wave_speed_x, swm_flux_x)
        ql  = state%qt(1,i,j-1,k,:)
        qr  = state%qb(1,i,j  ,k,:)
        state%fy(1,i,j,k,:) = riemann_solver(mesh%iG(:,:,mesh%pes(bottom),i,j), &
                                             mesh% J(    mesh%pes(bottom),i,j), &
                                             ql, qr, swm_wave_speed_y, swm_flux_y)
      end do
    end do
    if (mesh%jds == 1) then
      j = 1
      state%fy(:,:,j,:,:) = 0
    end if
    if (mesh%jde == mesh%ny) then
      j = mesh%jde + 1
      state%fy(:,:,j,:,:) = 0
    end if
    ! Change height back to depth.
    do j = mesh%jms, mesh%jme
      do i = mesh%ims, mesh%ime
        state%q(i,j,k,1) = state%q(i,j,k,1) - mesh%J(pc,i,j) * static%zs(i,j)
      end do
    end do
    ! Calculate sources and set tendencies.
    do j = mesh%jds, mesh%jde
      do i = mesh%ids, mesh%ide
        state%dqdt(i,j,k,:) = &
          - (state%fx(1,i+1,j,k,:) - state%fx(1,i,j,k,:)) / mesh%dx &
          - (state%fy(1,i,j+1,k,:) - state%fy(1,i,j,k,:)) / mesh%dy &
          + swm_sources(mesh%iG(:,:,pc,i,j), mesh%J(pc,i,j), mesh%CS(:,:,:,i,j), &
                        static%f1(i,j), state%q(i,j,k,:), static%zs(i,j), &
                        state%dhdx(i,j,k), state%dhdy(i,j,k))
      end do
    end do
    end associate

  end subroutine dycore_calc_swm_tend

  subroutine dycore_update_state(this, dt, dqdt_idx, old_q_idx, new_q_idx)

    class(dycore_type), intent(inout) :: this
    real(r8), intent(in) :: dt
    integer , intent(in) :: dqdt_idx
    integer , intent(in) :: old_q_idx
    integer , intent(in) :: new_q_idx

    integer i, j, k, pc

    pc = 1

    associate (mesh => this%mesh, dqdt => this%state(dqdt_idx)%dqdt, &
               old_q => this%state(old_q_idx)%q, new_q => this%state(new_q_idx)%q)
    do k = mesh%kds, mesh%kde
      do j = mesh%jds, mesh%jde
        do i = mesh%ids, mesh%ide
          new_q(i,j,k,:) = old_q(i,j,k,:) + dt * dqdt(i,j,k,:)
        end do
      end do
    end do
    end associate

    call swm_conservative_to_raw(this%state(new_q_idx), this%static)

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

  subroutine tvd3(dt, old, new, dycore)

    real(r8), intent(in) :: dt
    integer , intent(in) :: old
    integer , intent(in) :: new
    type(dycore_type), intent(inout) :: dycore

    integer, parameter :: s1 =  0
    integer, parameter :: s2 = -1
    real(r8), parameter :: c11 =  1.0_r8 / 2.0_r8, c12 = 1.0_r8 / 4.0_r8, c13 = 1.0_r8 / 4.0_r8
    real(r8), parameter :: c21 = -1.0_r8 / 3.0_r8, c22 = 2.0_r8 / 3.0_r8, c23 = 2.0_r8 / 3.0_r8

    associate (state => dycore%state)
    call dycore%calc_swm_tend(old)
    call dycore%update_state(dt, old, old, s1)

    call dycore%calc_swm_tend(s1)
    call dycore%update_state(dt, s1, old, s2)
    state(s2 )%q = c11 * state(old)%q + c12 * state(s1)%q + c13 * state(s2)%q

    call dycore%calc_swm_tend(s2)
    call dycore%update_state(dt, s2, old, s1)
    state(new)%q = c21 * state(old)%q + c22 * state(s2)%q + c23 * state(s1)%q
    end associate

  end subroutine tvd3

end module dycore_mod
