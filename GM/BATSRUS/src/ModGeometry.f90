!^CFG COPYRIGHT UM
module ModGeometry

  use ModSize
  use ModMain,       ONLY: UseBody2, body2_, ExtraBc_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  SAVE

  character(len=20):: TypeGeometry = 'cartesian'

  !\
  ! Geometry parameters.
  !/
  real :: x1, x2, y1, y2, z1, z2
  real :: DomainVolume = -1.0
  real :: xyzStart(3)
  real :: xyzStart_BLK(3,MaxBlock)
  real :: XyzMin_D(3)
  real :: XyzMax_D(3)
  real :: RadiusMin = -1.0, RadiusMax = -1.0

  ! Mirror symmetry in the 3 directions
  integer:: nMirror_D(3) = 1

  !\
  ! Other block solution and geometry parameters.
  !/
  real :: minDXvalue, maxDXvalue

  real, dimension(MaxBlock) :: dx_BLK, dy_BLK, dz_BLK, Rmin_BLK
  real, dimension(MaxBlock) :: Rmin2_BLK
  real, dimension(MaxBlock) :: fAx_BLK, fAy_BLK, fAz_BLK


  ! Variables describing cells inside boundaries
  !true when at least one cell in the block (including ghost cells) is not true
  logical :: body_BLK(MaxBlock)

  ! true when all cells in block (not including ghost cells) are true_cells 
  logical :: true_BLK(MaxBlock)

  ! true cells are cells that are not inside a body
  logical, allocatable :: true_cell(:,:,:,:)
  logical,dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn,body2_:Top_):: &
       IsBoundaryCell_GI
  logical,dimension(body2_:Top_,MaxBlock):: IsBoundaryBlock_IB 
  integer :: MinBoundary=Top_, MaxBoundary=body2_                    
  logical :: far_field_BCs_BLK(MaxBlock)

  ! Block cell coordinates
  real, allocatable :: x_BLK(:,:,:,:)
  real, allocatable :: y_BLK(:,:,:,:)
  real, allocatable :: z_BLK(:,:,:,:)
  real, allocatable :: R_BLK(:,:,:,:)
  real, allocatable :: R2_BLK(:,:,:,:)

  ! ADDED FOR general r grid in spherical geometry!
  ! Main Idea is to have a tabulated function that maps
  ! a general coordinate to log(r). This way, r resolution can 
  ! be arbitrarily defined. Gen. coord is taken to be 0 to 1.0
  ! but will be linearly extrapolated outside of this domain
  real, allocatable :: LogRGen_I(:)
  character(len=100):: NameGridFile=''

contains
  !============================================================================
  subroutine init_mod_geometry

    if(allocated(true_cell)) return
    allocate(true_cell(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(x_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(y_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(z_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(R_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(UseBody2) allocate(R2_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_geometry allocated arrays'
    end if

  end subroutine init_mod_geometry
  !============================================================================
  subroutine clean_mod_geometry

    if(.not.allocated(true_cell)) RETURN

    deallocate(true_cell)
    if(allocated(x_BLK))     deallocate(x_BLK)
    if(allocated(y_BLK))     deallocate(y_BLK)
    if(allocated(z_BLK))     deallocate(z_BLK)
    if(allocated(R_BLK))     deallocate(R_BLK)
    if(allocated(R2_BLK))    deallocate(R2_BLK)
    if(allocated(LogRGen_I)) deallocate(LogRGen_I)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_geometry deallocated arrays'
    end if

  end subroutine clean_mod_geometry
  !===========================================================================
  subroutine read_gen_radial_grid(NameFile)

    use ModIoUnit, ONLY: UnitTmp_

    character(len=*), intent(in) :: NameFile

    integer :: i, iError, nGrid
    real :: LogR

    character(len=*), parameter:: NameSub = 'read_gen_radial_grid'
    !------------------------------------------------------------------------
    ! This function is for reading in the general grid function.
    ! For simplicity the files Must be in a fixed format that is read by this
    ! routine.

    ! The formate is a single column with number of points in the 
    ! first row, and the folling rows the log(r) value from gen = 0 to 1

    open(UnitTmp_, FILE=NameFile, STATUS='old', IOSTAT=iError)
    if(iError /= 0) call CON_stop(NameSub// &
         ' could not open grid file = ' // trim(NameFile))

    ! read in nGrid 
    read(UnitTmp_,*,iostat=iError) nGrid

    if(iError /= 0) call CON_stop(NameSub// &
         ' could not read nGrid from file ' // trim(NameFile))

    ! Allocate LogRGen_I
    if(allocated(LogRGen_I)) deallocate(LogRGen_I)
    allocate(LogRGen_I(nGrid))
    LogRGen_I = 0.0

    ! read in tabulated log(radius)
    do i=1, nGrid
       read(UnitTmp_, *, IOSTAT=iError) LogR
       if(iError /= 0) then
          write(*,*) NameSub,': ERROR at i=',i
          call CON_stop(NameSub//' could not read logR from '//trim(NameFile))
       end if
       LogRGen_I(i) = LogR
    enddo

    close(UnitTmp_)

  end subroutine read_gen_radial_grid
  !===========================================================================
  subroutine set_gen_radial_grid

    !Set a single elemene

    if(allocated(LogRGen_I)) deallocate(LogRGen_I)
    allocate(LogRGen_I(1))
    LogRGen_I = 0.0 

  end subroutine set_gen_radial_grid

end module ModGeometry
