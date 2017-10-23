!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModGeometry

  use ModSize
  use ModMain,       ONLY: UseBody2, ExtraBc_, SolidBc_
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc
  use BATL_grid,     ONLY: Xyz_DGB, CellSize_DB

  implicit none
  SAVE

  ! Grid geometry description
  character(len=20):: TypeGeometry = 'cartesian'

  ! Cartesian box limits that may be cut out of a non-Cartesian grid
  real :: x1, x2, y1, y2, z1, z2

  ! Total volume of used (true) cells
  real :: DomainVolume = -1.0

  ! Coordinates of 1,1,1 cell. Coord111_DB would be a good name.
  ! Same as BATL_lib::CoordMin_DB + 0.5*CellSize_DB
  real :: XyzStart_BLK(3,MaxBlock)

  ! Obsolete variables. Same as BATL_lib::CoordMin_D and CoordMax_D
  real :: XyzMin_D(3), XyzMax_D(3)

  ! Coodinate limits in true radius and degrees (used in restart.H)
  real :: CoordDimMin_D(3) = 0.0, CoordDimMax_D(3) = 0.0
  real :: RadiusMin = -1.0, RadiusMax = -1.0
  
  ! Mirror symmetry in the 3 directions
  integer:: nMirror_D(3) = 1

  ! Smallest and largest cell sizes in the first dimension
  real:: CellSize1Min, CellSize1Max

  ! Smallest and largest cell sizes in either the first direction
  ! or the Phi direction in degrees if IsLogRadius or IsGenRadius is true
  real:: CellSizeMin, CellSizeMax

  ! Variables describing cells inside boundaries
  !true when at least one cell in the block (including ghost cells) is not true
  logical :: body_BLK(MaxBlock)

  ! true when all cells in block (not including ghost cells) are true_cells 
  logical :: true_BLK(MaxBlock)

  ! true cells are cells that are not inside a body
  logical, allocatable :: true_cell(:,:,:,:)

  ! Number of true cells (collected for processor 0)
  integer :: nTrueCells = -1

  ! True for blocks next to the cell based boundaries
  logical :: far_field_BCs_BLK(MaxBlock)

  ! Radial distance from origin and second body
  real, allocatable :: R_BLK(:,:,:,:)
  real, allocatable :: R2_BLK(:,:,:,:)

  ! Smallest value of r_BLK and r2_BLK within a block
  real:: Rmin_BLK(MaxBlock), Rmin2_BLK(MaxBlock)

  ! ADDED FOR general r grid in spherical geometry!
  ! Main Idea is to have a tabulated function that maps
  ! a general coordinate to log(r). This way, r resolution can 
  ! be arbitrarily defined. Gen. coord is taken to be 0 to 1.0
  ! but will be linearly extrapolated outside of this domain
  real, allocatable :: LogRGen_I(:)
  character(len=100):: NameGridFile=''

  ! Jacobian matrix for general grid: Dgencoord/Dcartesian
  ! This can be set with call set_block_jacobian_cell(iBlock)
  real, public :: DgenDxyz_DDC(MaxDim,MaxDim,nI,nJ,nK)

contains
  !============================================================================
  subroutine init_mod_geometry

    if(allocated(true_cell)) return
    allocate(true_cell(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
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
    use ModUtilities, ONLY: open_file, close_file

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

    call open_file(FILE=NameFile, STATUS='old')

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

    call close_file

  end subroutine read_gen_radial_grid
  !===========================================================================
  subroutine set_gen_radial_grid

    ! Set a single element array when generalized radial coordinate
    ! is not used, but still should be set.

    if(allocated(LogRGen_I)) deallocate(LogRGen_I)
    allocate(LogRGen_I(1))
    LogRGen_I = 0.0 

  end subroutine set_gen_radial_grid

  !============================================================================

  subroutine set_block_jacobian_cell(iBlock)

    use BATL_lib, ONLY: nDim
    use ModNumConst, ONLY: i_DD
    use ModCoordTransform, ONLY: inverse_matrix

    integer, intent(in):: iBlock
    real:: InvDx1Half, InvDx2Half, InvDx3Half
    real:: DxyzDgen_DD(MaxDim, MaxDim)
    integer:: i,j,k
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='set_block_jacobian'
    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Calculate the dCartesian/dGencoord matrix
    
    InvDx1Half = 0.5/CellSize_DB(1,iBlock)
    InvDx2Half = 0.5/CellSize_DB(2,iBlock)
    InvDx3Half = 0.5/CellSize_DB(3,iBlock)

    ! Get the dCartesian/dGencoord matrix with finite differences
    do k=1,nK; do j=1,nJ; do i=1,nI
       DxyzDgen_DD(:,1) = InvDx1Half &
            *(Xyz_DGB(:,i+1,j,k,iBlock) - Xyz_DGB(:,i-1,j,k,iBlock))

       DxyzDgen_DD(:,2) = InvDx2Half &
            *(Xyz_DGB(:,i,j+1,k,iBlock) - Xyz_DGB(:,i,j-1,k,iBlock))

       if(nDim==3)then
          DxyzDgen_DD(:,3) = InvDx3Half &
               *(Xyz_DGB(:,i,j,k+1,iBlock) - Xyz_DGB(:,i,j,k-1,iBlock))
       else
          DxyzDgen_DD(:,3) = i_DD(:,3)
       end if

       DgenDxyz_DDC(:,:,i,j,k) = &
            inverse_matrix(DxyzDgen_DD, DoIgnoreSingular=.true.)
    end do; end do; end do

  end subroutine set_block_jacobian_cell
 
  !============================================================================

  subroutine count_true_cells

    use BATL_lib, ONLY: nI, nJ, nK, nBlock, Unused_B, iComm
    use ModMpi

    integer :: iBlock, iError
    !----------------------------------------------------------------------------

    nTrueCells=0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       nTrueCells = nTrueCells + count(true_cell(1:nI,1:nJ,1:nK,iBlock))
    end do
    call MPI_reduce_integer_scalar(nTrueCells, MPI_SUM, 0, iComm, iError)

  end subroutine count_true_cells

end module ModGeometry
