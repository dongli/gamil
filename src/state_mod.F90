module state_mod

  use const_mod
  use kinds_mod
  use latlon_array_mod
  use latlon_parallel_mod

  implicit none

  type state_type
    logical :: initialized = .false.
    type(latlon_mesh_type), pointer :: mesh => null()
    type(latlon_array_type), allocatable :: array
    integer :: nvar = 0
    real(r8), pointer, dimension(:,:,:) :: h
    real(r8), pointer, dimension(:,:,:) :: u
    real(r8), pointer, dimension(:,:,:) :: v
    real(r8), pointer, dimension(:,:,:) :: uc
    real(r8), pointer, dimension(:,:,:) :: vc
    real(r8), pointer, dimension(:,:,:,:) :: q
    real(r8), pointer, dimension(:,:,:,:) :: ql
    real(r8), pointer, dimension(:,:,:,:) :: qr
    real(r8), pointer, dimension(:,:,:,:) :: qt
    real(r8), pointer, dimension(:,:,:,:) :: qb
    real(r8), pointer, dimension(:,:,:,:) :: fx
    real(r8), pointer, dimension(:,:,:,:) :: fy
    real(r8), pointer, dimension(:,:,:,:) :: dqdt
  contains
    procedure :: init  => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh, model_type)

    class(state_type), intent(inout) :: this
    type(latlon_mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: model_type

    call this%clear()

    allocate(this%array)

    this%mesh => mesh

    call this%array%init(mesh)

    select case (model_type)
    case ('swm')
      ! Raw variables
      call this%array%add_var('h' , 'Height'                   , 'm'    , loc='CA', output='T')
      call this%array%add_var('u' , 'Zonal wind component'     , 'm s-1', loc='CA', output='T')
      call this%array%add_var('v' , 'Meridional wind component', 'm s-1', loc='CA', output='T')
      call this%array%add_var('uc', 'Contravariant u'          , 's-1'  , loc='CA')
      call this%array%add_var('uv', 'Contravariant u'          , 's-1'  , loc='CA')

      ! Conservative variables
      this%nvar = 3
      call this%array%add_var('q1', loc='CA:LQ:RQ:TQ:BQ', with_halo='T:T:T:T:T', fill_halo='T:F:F:F:F', tag='fcst')
      call this%array%add_var('q2', loc='CA:LQ:RQ:TQ:BQ', with_halo='T:T:T:T:T', fill_halo='T:F:F:F:F', tag='fcst')
      call this%array%add_var('q2', loc='CA:LQ:RQ:TQ:BQ', with_halo='T:T:T:T:T', fill_halo='T:F:F:F:F', tag='fcst')
      call this%array%add_var('dq1dt', loc='CA', tag='tend')
      call this%array%add_var('dq2dt', loc='CA', tag='tend')
      call this%array%add_var('dq3dt', loc='CA', tag='tend')
    end select

    call this%array%add_var('fx1', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fx')
    call this%array%add_var('fx2', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fx')
    call this%array%add_var('fx3', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fx')
    call this%array%add_var('fy1', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fy')
    call this%array%add_var('fy2', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fy')
    call this%array%add_var('fy3', loc='LQ:BQ', with_halo='T:T', fill_halo='F:F', tag='fy')

    call this%array%allocate_arrays()

    select case (model_type)
    case ('swm')
      call this%array%get_array(this%h , var_name='h' , loc='CA')
      call this%array%get_array(this%u , var_name='u' , loc='CA')
      call this%array%get_array(this%v , var_name='v' , loc='CA')
      call this%array%get_array(this%uc, var_name='uc', loc='CA')
      call this%array%get_array(this%vc, var_name='vc', loc='CA')
    end select

    call this%array%get_array(this%q , tag='fcst', loc='CA')
    call this%array%get_array(this%ql, tag='fcst', loc='LQ')
    call this%array%get_array(this%qr, tag='fcst', loc='RQ')
    call this%array%get_array(this%qt, tag='fcst', loc='TQ')
    call this%array%get_array(this%qb, tag='fcst', loc='BQ')

    call this%array%get_array(this%dqdt, tag='tend', loc='CA')

    call this%array%get_array(this%fx, tag='fx'  , loc='LQ')
    call this%array%get_array(this%fy, tag='fy'  , loc='BQ')

    this%q = inf

    this%initialized = .true.

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    nullify(this%mesh)
    this%nvar = 0
    if (allocated(this%array)) deallocate(this%array)

    this%initialized = .false.

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
