module latlon_mesh_mod

  use fiona
  use kinds_mod
  use var_types_mod
  use mesh_const_mod
  use mesh_math_mod

  implicit none

  type latlon_mesh_type
    logical :: initialized = .false.
    integer :: domain_type = global_domain
    integer :: nx  = 0                                  ! Cell  size  along x-axis
    integer :: ny  = 0                                  ! Cell  size  along y-axis
    integer :: nz  = 0                                  ! Cell  size  along z-axis
    integer :: ims = 0, ime = 0, ids = 0, ide = 0       ! Index range along x-axis
    integer :: jms = 0, jme = 0, jds = 0, jde = 0       ! Index range along y-axis
    integer :: kms = 1, kme = 1, kds = 1, kde = 1       ! Index range along z-axis
    integer :: hwx = 0                                  ! Halo  width along x-axis
    integer :: hwy = 0                                  ! Halo  width along y-axis
    integer :: hwz = 0                                  ! Halo  width along z-axis
    integer :: nw  = 0                                  ! Extra halo width for easing loop writing
    ! Point type and their indicators
    integer :: npt  = 0                                 ! Number of points in each cell (e.g., cell center, vertex)
    integer :: ncq  = 0                                 ! Number of quadrature points in cell
    integer :: neq  = 1                                 ! Number of quadrature points along each horizontal edge
    integer :: neqv = 1                                 ! Number of quadrature points along each vertical edge
    integer :: pc   = 1                                 ! Point index for cell center
    integer :: pvs  = 2                                 ! Point start index for vertices
    integer :: pve  = 5                                 ! Point end index for vertices
    integer :: pes(4)                                   ! Point start index for quadrature points on edges
    integer :: pee(4)                                   ! Point end index for quadrature points on edges
    integer :: pqs                                      ! Point start index for quadrature points in cell
    integer :: pqe                                      ! Point end index for quadrature points in cell
    real(r8) :: r                                       ! Sphere radius
    real(r8) :: xmin = 0, xmax = 0                      ! Longitude range
    real(r8) :: ymin = 0, ymax = 0                      ! Latitude range
    real(r8) :: dx, dy
    real(r8), allocatable, dimension(:        ) :: xeq  ! Gaussian-Legendre quadrature points
    real(r8), allocatable, dimension(:        ) :: weq  ! Gaussian-Legendre quadrature weights
    real(r8), allocatable, dimension(    :,:,:) :: lon
    real(r8), allocatable, dimension(    :,:,:) :: lat
    real(r8), allocatable, dimension(:,:,:,:,:) ::  G
    real(r8), allocatable, dimension(:,:,:,:,:) :: iG
    real(r8), allocatable, dimension(    :,:,:) ::  J   ! Jacobian (e.g., sqrt(det(G)))
    real(r8), allocatable, dimension(:,:,:,:,:) :: CS   ! Christoffel symbols
  contains
    procedure :: init        => latlon_mesh_init
    procedure :: set_metrics => latlon_mesh_set_metrics
    procedure :: get_params  => latlon_mesh_get_params
    procedure :: clear       => latlon_mesh_clear
    procedure :: write       => latlon_mesh_write
    final :: latlon_mesh_final
  end type latlon_mesh_type

contains

  subroutine latlon_mesh_init(this, nx, ny, nz, dx, dy, rlon0, rlat0, xmin, xmax, ymin, ymax, hwx, hwy, hwz, nw, neq, r, &
                              ids, ide, jds, jde)

    class(latlon_mesh_type), intent(inout) :: this
    integer , intent(in)           :: nx
    integer , intent(in)           :: ny
    integer , intent(in)           :: nz
    real(r8), intent(in), optional :: dx
    real(r8), intent(in), optional :: dy
    real(r8), intent(in), optional :: rlon0
    real(r8), intent(in), optional :: rlat0
    real(r8), intent(in), optional :: xmin
    real(r8), intent(in), optional :: xmax
    real(r8), intent(in), optional :: ymin
    real(r8), intent(in), optional :: ymax
    integer , intent(in), optional :: hwx
    integer , intent(in), optional :: hwy
    integer , intent(in), optional :: hwz
    integer , intent(in), optional :: nw
    integer , intent(in), optional :: neq
    real(r8), intent(in), optional :: r
    integer , intent(in), optional :: ids
    integer , intent(in), optional :: ide
    integer , intent(in), optional :: jds
    integer , intent(in), optional :: jde

    integer i, j, p, q1, q2

    call this%clear()

    this%nx = nx
    this%ny = ny
    this%nz = nz

    if (present(hwx)) this%hwx = hwx
    if (present(hwy)) this%hwy = hwy
    if (present(hwz)) this%hwz = hwz
    if (present(nw )) this%nw  = nw
    if (present(neq)) this%neq = neq
    if (present(r  )) this%r   = r

    this%ids = merge(ids, 1 , present(ids))
    this%ide = merge(ide, nx, present(ide))
    this%ims = this%ids - this%hwx - this%nw
    this%ime = this%ide + this%hwx + this%nw
    this%jds = merge(jds, 1 , present(jds))
    this%jde = merge(jde, ny, present(jde))
    this%jms = this%jds - this%hwy - this%nw
    this%jme = this%jde + this%hwy + this%nw
    this%kds = 1
    this%kde = nz
    this%kms = this%kds - this%hwz
    this%kme = this%kde + this%hwz
    this%npt = 5 + 4 * this%neq + this%neq**2
    this%ncq = this%neq**2

    if (present(xmin) .and. present(xmax)) then
      this%xmin = xmin
      this%xmax = xmax
      this%dx   = (this%xmax - this%xmin) / (this%ide - this%ids)
    else if (present(dx)) then
      this%xmin = 0.5_r8 * dx + dx * (this%ids - 1) - 0.5_r8 * dx
      this%xmax = 0.5_r8 * dx + dx * (this%ide - 1) + 0.5_r8 * dx
      this%dx   = dx
    end if

    if (present(ymin) .and. present(ymax)) then
      this%ymin = ymin
      this%ymax = ymax
      this%dy   = (this%ymax - this%ymin) / (this%jde - this%jds)
    else if (present(dy)) then
      this%ymin = -pi0p5 + 0.5_r8 * dy + dy * (this%jds - 1) - 0.5_r8 * dy
      this%ymax = -pi0p5 + 0.5_r8 * dy + dy * (this%jde - 1) + 0.5_r8 * dy
      this%dy   = dy
    end if

    allocate(this%xeq(this%neq))
    allocate(this%weq(this%neq))
    allocate(this%lon(    this%npt,this%ims:this%ime,this%jms:this%jme))
    allocate(this%lat(    this%npt,this%ims:this%ime,this%jms:this%jme))
    allocate(this%G  (3,3,this%npt,this%ims:this%ime,this%jms:this%jme))
    allocate(this%iG (3,3,this%npt,this%ims:this%ime,this%jms:this%jme))
    allocate(this%J  (    this%npt,this%ims:this%ime,this%jms:this%jme))
    allocate(this%CS (3,3,3       ,this%ims:this%ime,this%jms:this%jme))

    if (this%neq == 1) then
      this%xeq = 0.5_r8
      this%weq = 1.0_r8
    else
      call gaussian_legendre(this%neq, this%xeq, this%weq)
      this%xeq = (this%xeq + 1.0_r8) * 0.5_r8
      this%weq = this%weq / 2.0_r8
    end if

    ! Rotate coordinate
    if (present(rlon0) .and. present(rlat0)) then

    else
      do j = this%jms, this%jme
        do i = this%ims, this%ime
          ! Cell center
          !   o---------------o
          !   |               |
          !   |               |
          !   |               |
          !   |       1       |
          !   |               |
          !   |               |
          !   |               |
          !   o---------------o
          !
          this%pc = 1
          this%lon(1,i,j) = this%xmin + (i - this%ids + 0.5_r8) * this%dx
          this%lat(1,i,j) = this%ymin + (j - this%jds + 0.5_r8) * this%dy
          ! Vertices
          ! 5 o---------------o 4
          !   |               |
          !   |               |
          !   |               |
          !   |               |
          !   |               |
          !   |               |
          !   |               |
          ! 2 o---------------o 3
          !
          this%pvs = 2
          this%pve = 5
          this%lon(2,i,j) = this%lon(1,i,j) - 0.5_r8 * this%dx
          this%lat(2,i,j) = this%lat(1,i,j) - 0.5_r8 * this%dy
          this%lon(3,i,j) = this%lon(1,i,j) + 0.5_r8 * this%dx
          this%lat(3,i,j) = this%lat(1,i,j) - 0.5_r8 * this%dy
          this%lon(4,i,j) = this%lon(1,i,j) + 0.5_r8 * this%dx
          this%lat(4,i,j) = this%lat(1,i,j) + 0.5_r8 * this%dy
          this%lon(5,i,j) = this%lon(1,i,j) - 0.5_r8 * this%dx
          this%lat(5,i,j) = this%lat(1,i,j) + 0.5_r8 * this%dy
          ! Edge quadrature with 2 points as an example
          !     10        11
          !   o--x---------x--o
          !   |               |
          !13 x               x 9
          !   |               |
          !   |               |
          !   |               |
          !12 x               x 8
          !   |               |
          !   o--x---------x--o
          !      6         7
          p = 6
          ! Bottom edge points
          this%pes(bottom) = p
          do q1 = 1, this%neq
            this%lon(p,i,j) = this%lon(2,i,j) + this%dx * this%xeq(q1)
            this%lat(p,i,j) = this%lat(2,i,j)
            p = p + 1
          end do
          this%pee(bottom) = p - 1
          ! Right edge points
          this%pes(right) = p
          do q1 = 1, this%neq
            this%lon(p,i,j) = this%lon(3,i,j) 
            this%lat(p,i,j) = this%lat(3,i,j) + this%dy * this%xeq(q1)
            p = p + 1
          end do
          this%pee(right) = p - 1
          ! Top edge points
          this%pes(top) = p
          do q1 = 1, this%neq
            this%lon(p,i,j) = this%lon(5,i,j) + this%dx * this%xeq(q1)
            this%lat(p,i,j) = this%lat(5,i,j)
            p = p + 1
          end do
          this%pee(top) = p - 1
          ! Left edge points
          this%pes(left) = p
          do q1 = 1, this%neq
            this%lon(p,i,j) = this%lon(2,i,j)
            this%lat(p,i,j) = this%lat(2,i,j) + this%dy * this%xeq(q1)
            p = p + 1
          end do
          this%pee(left) = p - 1
          ! Cell quadrature with 4 points as an example
          !   o---------------o
          !   | 16        17  |
          !   |  x         x  |
          !   |               |
          !   |               |
          !   |               |
          !   |  x         x  |
          !   | 14        15  |
          !   o---------------o
          this%pqs = p
          do q2 = 1, this%neq
            do q1 = 1, this%neq
              this%lon(p,i,j) = this%lon(2,i,j) + this%dx * this%xeq(q2)
              this%lat(p,i,j) = this%lat(2,i,j) + this%dy * this%xeq(q1)
              p = p + 1
            end do
          end do
          this%pqe = p - 1
        end do
      end do
    end if

    call this%set_metrics(this%r)

  end subroutine latlon_mesh_init

  subroutine latlon_mesh_set_metrics(this, r)

    class(latlon_mesh_type), intent(inout) :: this
    real(r8), intent(in) :: r

    real(r8) lon, lat, cos_lat, sin_lat
    integer i, j, p

    do j = this%jms, this%jme
      do i = this%ims, this%ime
        do p = 1, this%npt
          lon = this%lon(p,i,j)
          lat = this%lat(p,i,j)
          cos_lat = cos(lat); cos_lat = merge(0.0_r8, cos_lat, abs(cos_lat) < 1.0e-16)
          sin_lat = sin(lat); sin_lat = merge(0.0_r8, sin_lat, abs(sin_lat) < 1.0e-16)
          this% G(1,1,p,i,j) = (r * cos_lat)**2
          this% G(2,2,p,i,j) = r**2
          this% G(3,3,p,i,j) = 1.0_r8
          this%iG(1,1,p,i,j) = 1.0_r8 / this%G(1,1,p,i,j)
          this%iG(2,2,p,i,j) = 1.0_r8 / this%G(2,2,p,i,j)
          this%iG(3,3,p,i,j) = 1.0_r8 / this%G(3,3,p,i,j)
          this% J(    p,i,j) = r**2 * cos_lat
          if (p /= 1) cycle
          ! FIXME: Christoffel symbols are only needed on cell centers?
          this%CS(1,2,1,i,j) = -tan(lat)
          this%CS(2,1,1,i,j) = -tan(lat)
          this%CS(1,1,2,i,j) =  sin(lat) * cos_lat
        end do
      end do
    end do

  end subroutine latlon_mesh_set_metrics

  subroutine latlon_mesh_get_params(this, ims, ime, ids, ide, nx, &
                                          jms, jme, jds, jde, ny, &
                                          kms, kme, kds, kde, nz, &
                                          hwx, hwy, hwz, pc ,     &
                                          neq, ncq)

    class(latlon_mesh_type), intent(in) :: this
    integer, intent(out), optional :: ims, ime, ids, ide, nx
    integer, intent(out), optional :: jms, jme, jds, jde, ny
    integer, intent(out), optional :: kms, kme, kds, kde, nz
    integer, intent(out), optional :: hwx, hwy, hwz, pc , neq, ncq

    if (present(ims)) ims = this%ims
    if (present(ime)) ime = this%ime
    if (present(ids)) ids = this%ids
    if (present(ide)) ide = this%ide
    if (present(nx )) nx  = this%nx
    if (present(jms)) jms = this%jms
    if (present(jme)) jme = this%jme
    if (present(jds)) jds = this%jds
    if (present(jde)) jde = this%jde
    if (present(ny )) ny  = this%ny
    if (present(kms)) kms = this%kms
    if (present(kme)) kme = this%kme
    if (present(kds)) kds = this%kds
    if (present(kde)) kde = this%kde
    if (present(nz )) nz  = this%nz
    if (present(hwx)) hwx = this%hwx
    if (present(hwy)) hwy = this%hwy
    if (present(hwz)) hwz = this%hwz
    if (present(pc )) pc  = this%pc
    if (present(neq)) neq = this%neq
    if (present(ncq)) ncq = this%ncq

  end subroutine latlon_mesh_get_params

  subroutine latlon_mesh_clear(this)

    class(latlon_mesh_type), intent(inout) :: this

    if (allocated(this%xeq)) deallocate(this%xeq)
    if (allocated(this%weq)) deallocate(this%weq)
    if (allocated(this%lon)) deallocate(this%lon)
    if (allocated(this%lat)) deallocate(this%lat)
    if (allocated(this%G  )) deallocate(this%G  )
    if (allocated(this%iG )) deallocate(this%iG )
    if (allocated(this%J  )) deallocate(this%J  )
    if (allocated(this%CS )) deallocate(this%CS )

  end subroutine latlon_mesh_clear

  subroutine latlon_mesh_write(this, file_path)

    class(latlon_mesh_type), intent(inout) :: this
    character(*), intent(in) :: file_path

    call fiona_create_dataset('mesh', file_path=file_path)
    call fiona_add_dim('mesh', 'lon', 'Longitude of cell center', size=this%nx)
    call fiona_add_dim('mesh', 'lat', 'Latitude of cell center' , size=this%ny)
    call fiona_add_dim('mesh', 'three', size=3)
    call fiona_add_var('mesh', 'xlon', 'Longitude of cell center', units='degree_east', dim_names=['lon','lat'], data_type='r8')
    call fiona_add_var('mesh', 'xlat', 'Latitude of cell center' , units='degree_east', dim_names=['lon','lat'], data_type='r8')
    call fiona_add_var('mesh', 'J', 'Jacobian of cell center', units='1', dim_names=['lon','lat'], data_type='r8')
    call fiona_add_att('mesh', 'J', 'coordinates', 'xlon xlat')
    call fiona_start_output('mesh')
    call fiona_output('mesh', 'xlon', this%lon(1,this%ids:this%ide,this%jds:this%jde) * deg)
    call fiona_output('mesh', 'xlat', this%lat(1,this%ids:this%ide,this%jds:this%jde) * deg)
    call fiona_output('mesh', 'J'   , this%J  (1,this%ids:this%ide,this%jds:this%jde))
    call fiona_end_output('mesh')

  end subroutine latlon_mesh_write

  subroutine latlon_mesh_final(this)

    type(latlon_mesh_type), intent(inout) :: this

    call this%clear()

  end subroutine latlon_mesh_final

end module latlon_mesh_mod
