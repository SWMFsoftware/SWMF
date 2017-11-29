!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_region

  ! Define an array of geometric regions (areas)
  ! that can be used for AMR and other criteria.
  ! For AMR regions defined by the #GRIDLEVEL/RESOLUTION commands
  ! also store the required grid resolution/level.
  ! Also read the initial grid refinement resolution/level.

  use BATL_mpi,  ONLY: iProc
  use BATL_size, ONLY: nDim, Dim2_, Dim3_, j0_, nJp1_, k0_, nKp1_, &
       nI, nJ, nK, nIJK, MinI, MinJ, MinK, MaxI, MaxJ, MaxK, MaxIJK,&
       nINode, nJNode, nKNode

  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  private ! except

  public:: read_region_param
  public:: init_region
  public:: clean_region
  public:: i_signed_region
  public:: get_region_indexes
  public:: block_inside_regions
  public:: points_inside_region

  ! Maximum number of geometric areas
  integer, public, parameter :: MaxArea = 100

  ! Number of areas defined
  integer, public :: nArea = 0

  ! Number of geometrical criteria from #GRIDLEVEL/RESOLUTION commands
  integer, public :: nCritGrid    = 0

  ! we have read in new data from PARAM.in
  logical, public :: IsNewGeoParam =.false.

  ! Lengths of strings
  integer, parameter:: lNameArea = 20, lNameRegion = 30

  type, public :: AreaType
     character(lNameRegion):: NameRegion
     character(lNameArea)  :: NameShape
     real                  :: Resolution
     integer               :: Level
     real                  :: Center_D(nDim)
     real                  :: Size_D(nDim)
     real                  :: Radius1
     real                  :: Taper
     logical               :: DoRotate
     real, allocatable     :: Rotate_DD(:,:)
     logical               :: IsSimple
  end type AreaType

  ! List of all geometric regions defined
  type(AreaType), target, public :: Area_I(1:MaxArea)

  ! the area being processed (for easier access)
  type(AreaType), pointer, public:: Area

  ! These are set by #GRIDLEVEL/RESOLTION ... initial command
  ! Values are corrected and used in init_region
  ! The number of levels is public so that BATSRUS can do the
  ! initial refinement loop.
  integer, public :: nInitialAmrLevel  = 0

  ! Local variables
  real            :: InitialResolution = -1.0

  ! Name of shape being processed
  character(lNameArea) :: NameShape

  ! index for the parallel and one or two perpendicular directions
  integer:: iPar
  integer, allocatable:: iPerp_I(:)
  real,    allocatable:: SlopePerp_D(:)

  ! Allocatable storage
  logical, allocatable:: IsInsideOld_I(:)
  real,    allocatable:: ValueOld_I(:), Xyz_DI(:,:), Coord_DI(:,:)
  real,    allocatable:: Norm_DI(:,:)
  real,    allocatable:: Corner_DI(:,:), CornerNorm_DI(:,:)

contains
  !============================================================================
  subroutine init_region

    use BATL_geometry,     ONLY: Phi_, Theta_, IsCartesianGrid, &
         IsLogRadius, IsGenRadius, radius_to_gen, CellSizeRoot
    use ModNumConst,       ONLY: cDegToRad

    integer :: iGeo
    real:: rMin, rMax

    logical, parameter :: DoTest = .false.
    character(len=*), parameter:: NameSub = 'init_region'
    !--------------------------------------------------------------------------

    ! Set initial resolutions (this depends on domain size set in BATL_grid)
    if(InitialResolution > 0.0) nInitialAmrLevel &
         = ceiling(log(CellSizeRoot/InitialResolution)/log(2.0) - 0.01)

    if(DoTest)write(*,*) NameSub, &
         ' CellSizeRoot, InitialResolution, nInitialAmrLevel=', &
         CellSizeRoot, InitialResolution, nInitialAmrLevel

    ! Set IsSimple for all the areas
    do iGeo = 1, nArea
       Area => Area_I(iGeo)
       Area%IsSimple = .false.
       ! generalized coordinates
       if(Area%NameShape == "brick_gen")then
          Area%IsSimple = .true.
          CYCLE
       end if
       if(.not.IsCartesianGrid)CYCLE
       select case(Area%NameShape(1:4))
       case('sphe','shel','bric','cyli','ring')
          Area%IsSimple = .not.Area%DoRotate
       end select
    end do

    ! Fix generalized coordinates for "brick_gen" and rename to "brick_coord"
    do iGeo = 1, nArea
       Area => Area_I(iGeo)
       if( Area % NameShape /= "brick_gen") CYCLE

       ! Rename area to indicate that it has been processed
       Area%NameShape = "brick_coord"

       ! Recalculate area center and size in generalized coordinates
       if(IsLogRadius .or.IsGenRadius)then
          rMin = max(1e-30, Area % Center_D(1) - Area % Size_D(1))
          rMax = max(1e-30, Area % Center_D(1) + Area % Size_D(1))
          if(IsLogRadius) rMin = log(rMin)
          if(IsLogRadius) rMax = log(rMax)
          if(IsGenRadius) call radius_to_gen(rMin)
          if(IsGenRadius) call radius_to_gen(rMax)
          Area % Center_D(1) = 0.5*(rMax + rMin)
          Area % Size_D(1)   = 0.5*(rMax - rMin)
       endif

       ! Convert degrees to radians
       if(Phi_ > 1)   Area%Center_D(Phi_)  = Area%Center_D(Phi_)*cDegToRad
       if(Phi_ > 1)   Area%Size_D(Phi_)    = Area%Size_D(Phi_)*cDegToRad
       if(Theta_ > 1) Area%Center_D(Theta_)= Area%Center_D(Theta_)*cDegToRad
       if(Theta_ > 1) Area%Size_D(Theta_)  = Area%Size_D(Theta_)*cDegToRad
    end do

    ! Loop over the areas defined by #GRIDLEVEL/RESOLUTION commands
    ! That do not have a name
    do iGeo = 1, nArea
       Area => Area_I(iGeo)

       ! Exclude named areas defined by #REGION / #AMRREGION commands
       if( Area%NameRegion /= "NULL") CYCLE

       ! Set level and resolution based on the read value
       if(Area%Level < 0)then
          ! Set resolution corresponding to this level
          Area%Resolution = CellSizeRoot * 2.0**Area%Level
       else
          ! Set AMR level (note that negative value is used)
          Area%Level = &
               -ceiling(log(CellSizeRoot/Area%Resolution)/log(2.0) - 0.01)
       end if

    end do

    if(DoTest .and. iProc == 0) then
       do iGeo = 1, nArea
          call show_region(iGeo, NameSub)
          write(*,*) "-----------------------"
       end do
    end if

    IsNewGeoParam =.false.

  end subroutine init_region
  !============================================================================

  subroutine clean_region

    !--------------------------------------------------------------------------
    if(allocated(iPerp_I))        deallocate(iPerp_I)
    if(allocated(SlopePerp_D))    deallocate(SlopePerp_D)
    if(allocated(IsInsideOld_I))  deallocate(IsInsideOld_I)
    if(allocated(ValueOld_I))     deallocate(ValueOld_I)
    if(allocated(Xyz_DI))         deallocate(Xyz_DI)
    if(allocated(Coord_DI))       deallocate(Coord_DI)
    if(allocated(Norm_DI))        deallocate(Norm_DI)
    if(allocated(Corner_DI))      deallocate(Corner_DI)
    if(allocated(CornerNorm_DI))  deallocate(CornerNorm_DI)

    nArea                 = 0
    IsNewGeoParam         = .false.
    InitialResolution     = -1.0
    nInitialAmrLevel      = 0

  end subroutine clean_region
  !============================================================================
  subroutine show_region(iRegion, String)

    integer, intent(in):: iRegion
    character(len=*), intent(in):: String

    type(AreaType), pointer:: Area1
    !--------------------------------------------------------------------------
    Area1 => Area_I(iRegion)
    write(*,*) "show_region for ", String,' iRegion=', iRegion
    write(*,*) "Region name      :: ", Area1%NameRegion
    write(*,*) "Shape name       :: ", Area1%NameShape
    write(*,*) "Region resolution:: ", Area1%Resolution
    write(*,*) "Region level     :: ", Area1%Level
    write(*,*) "Region Center_D  :: ", Area1%Center_D
    write(*,*) "Region Size_D    :: ", Area1%Size_D
    write(*,*) "Region Radius1   :: ", Area1%Radius1
    write(*,*) "Region Taper     :: ", Area1%Taper
    write(*,*) "Region DoRotate  :: ", Area1%DoRotate
    if(Area1%DoRotate) write(*,*) &
         "Region Rotate_DD :: ",       Area1%Rotate_DD
    write(*,*) "Region IsSimple  :: ", Area1%IsSimple
  end subroutine show_region
  !============================================================================

  integer function i_signed_region(NameRegionIn)

    ! Return the index of the region matching NameRegionIn.
    ! NameRegionIn may start with a '+' or '-' character.

    character(len=*), intent(in):: NameRegionIn

    character(lNameRegion):: NameRegion
    integer:: iArea, iSign

    character(len=*), parameter:: NameSub = 'i_signed_region'
    !--------------------------------------------------------------------------
    if(NameRegionIn(1:1) == "-") then
       iSign = -1
       NameRegion = NameRegionIn(2:)
    elseif (NameRegionIn(1:1) == "+") then
       NameRegion = NameRegionIn(2:)
       iSign = 1
    else
       NameRegion = NameRegionIn
       iSign = 1
    end if

    do iArea = 1, nArea
       if(NameRegion == Area_I(iArea)%NameRegion) then
          i_signed_region = iSign*iArea
          RETURN
       end if
    end do

    call CON_stop(NameSub//' could not find NameRegion='//NameRegionIn)

  end function i_signed_region
  !============================================================================
  subroutine get_region_indexes(StringRegion, iRegion_I)

    use ModUtilities, ONLY: split_string

    ! Create an array of signed region indexes from the region names
    ! listed in StringRegion. The names can have an optional "+" or "-"
    ! sign in front of them, which will be stored as the sign of the index.

    character(len=*), intent(in):: StringRegion
    integer, allocatable, intent(inout):: iRegion_I(:)

    integer, parameter:: MaxName = 100
    integer:: nName, iName
    character(lNameRegion):: NameRegion_I(MaxName)
    !--------------------------------------------------------------------------
    if(allocated(iRegion_I)) deallocate(iRegion_I)
    if(StringRegion == 'none') RETURN

    call split_string(StringRegion, MaxName, NameRegion_I, nName)
    allocate(iRegion_I(nName))
    do iName = 1, nName
       iRegion_I(iName) = i_signed_region(NameRegion_I(iName))
    end do

  end subroutine get_region_indexes
  !============================================================================
  subroutine set_i_par_perp

    ! Set the indexes parallel and perpendicular to the symmetry axis
    ! of various shapes

    integer:: l
    character:: CharLast

    !--------------------------------------------------------------------------
    iPar = 0 ! default is that there is no symmetry axis
    if(nDim == 1) RETURN

    if(.not.allocated(iPerp_I)) allocate(iPerp_I(nDim-1))

    l = len_trim(NameShape)
    CharLast = NameShape(l:l)

    select case(CharLast)
    case('x')
       iPar = 1; iPerp_I(1) = 2
    case('y')
       iPar = 2; iPerp_I(1) = 1
    case('z')
       iPar = 3; iPerp_I(1) = 1
    end select
    if(nDim == 3) iPerp_I(2) = 6 - iPar - iPerp_I(1)

  end subroutine set_i_par_perp
  !============================================================================

  subroutine read_region_param(NameCommand, UseStrictIn)

    use ModReadParam,      ONLY: read_var, lStringLine
    use ModNumConst,       ONLY: cDegToRad
    use ModCoordTransform, ONLY: &
         rot_matrix, rot_matrix_x, rot_matrix_y, rot_matrix_z

    character(len=*),  intent(inout) :: NameCommand
    logical, optional, intent(in)    :: UseStrictIn

    logical :: UseStrict

    character(len=lStringLine):: StringShape = 'all'
    real    :: RadiusArea    =0.0
    logical :: DoReadAreaCenter = .false.
    real    :: XyzStartArea_D(nDim)=0.0, XyzEndArea_D(nDim)=0.0
    real    :: xRotateArea   = 0.0
    real    :: yRotateArea   = 0.0
    real    :: zRotateArea   = 0.0
    logical :: DoTaperArea   = .false.
    logical :: DoStretch     = .false.

    character(lNameRegion):: NameRegion
    integer :: i, iDim

    real    :: AreaResolution  = 0.0
    integer :: nLevelArea = 0

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'read_region_param'
    !--------------------------------------------------------------------------

    UseStrict = .true.
    if(present(UseStrictIn)) UseStrict = UseStrictIn

    ! Default values
    NameRegion     = 'NULL'
    AreaResolution = 0.0
    nLevelArea     = 0
    if(index(NameCommand, "REGION") > 0) then
       ! Read name of region for #REGION (or #AMRREGION) command
       call read_var('NameRegion', NameRegion, IsLowerCase=.true.)
    elseif(index(NameCommand,"RESOLUTION") > 0)then
       ! Read grid resolution for #RESOLUTION command
       call read_var('AreaResolution', AreaResolution)
    else
       ! Read grid level for #GRIDLEVEL command
       call read_var('nLevelArea', nLevelArea)
       ! Requesting 0 level refinement is meaningless
       if(nLevelArea <= 0) RETURN
    end if
    call read_var('StringShape', StringShape, IsLowerCase=.true.)
    StringShape = adjustl(StringShape)

    ! 'init' or 'initial' means that the initial resolution is set,
    ! and no area is created.
    if(StringShape(1:4) == 'init')then
       if(AreaResolution > 0)then
          InitialResolution = AreaResolution
       else
          nInitialAmrLevel = nLevelArea
       end if

       if(DoTest)write(*,*) NameSub, &
            ' setting InitialResolution, nInitialAmrLevel =', &
            InitialResolution, nInitialAmrLevel
       RETURN
    end if

    nArea = nArea + 1
    if(nArea > MaxArea)then
       if(UseStrict) &
            call CON_stop(NameSub//' ERROR: Too many grid areas were defined')
       if(iProc == 0)then
          write(*,*) NameSub," Too many regions! MaxArea = ", MaxArea
          write(*,*) NameSub," ignoring command ", NameCommand, &
               " with StringShape=", trim(StringShape)
       end if
       nArea = MaxArea
       RETURN
    end if

    Area => Area_I(nArea)

    ! Regions without a name are produced by #GRIDLEVEL/RESOLUTION commands
    if(NameRegion == "NULL")then
       nCritGrid = nCritGrid + 1
       IsNewGeoParam =.true.
    end if

    ! Store the information read above
    Area%NameRegion = NameRegion
    Area%Resolution = AreaResolution
    Area%Level      = -nLevelArea    ! Level is stored as a negative integer

    ! Now process the information read into StringShape

    ! Set defaults
    Area%Center_D = 0.0
    Area%Size_D   = 1.0
    Area%Radius1  = 0.0
    Area%Taper    = 0.0
    Area%DoRotate = .false.

    ! Check for the word rotated in the name
    i = index(StringShape, 'rotated')
    Area%DoRotate = (i > 0 .and. nDim > 1)
    if(i > 0) StringShape = StringShape(1:i-1)//StringShape(i+7:)

    ! check for the word stretched in the name
    i = index(StringShape, 'stretched')
    DoStretch = i > 0
    if(i > 0) StringShape = StringShape(:i-1)//StringShape(i+9:)

    ! check for the word tapered in the name
    i = index(StringShape,'tapered')
    DoTaperArea = i > 0

    if(i > 0) StringShape = StringShape(:i-1)//StringShape(i+7:)

    ! Extract character '0' from the name
    i = index(StringShape,'0')
    if(i > 0) StringShape = StringShape(:i-1)//StringShape(i+1:)

    DoReadAreaCenter = (i < 1 .and. StringShape(1:3) /= 'box')

    ! Set name of the shape after all options were removed from StringShape
    NameShape = StringShape

    ! Store NameShape
    Area%NameShape = NameShape

    ! The 'all' types does not need any parameters
    if(NameShape == 'all' .or. NameShape(1:4) == 'user') RETURN

    ! Read center of area if needed
    if(DoReadAreaCenter)then
       do iDim = 1, nDim
          call read_var("Position", Area%Center_D(iDim))
       end do
    endif

    ! Read shape parameters
    call set_i_par_perp
    select case(NameShape)
    case("box", "box_gen")
       do iDim = 1, nDim
          call read_var("XyzMinBox", XyzStartArea_D(iDim))
       end do
       do iDim = 1, nDim
          call read_var("XyzMaxBox", XyzEndArea_D(iDim))
       end do
       ! Convert to center and size information
       Area%Center_D= 0.5*   (XyzStartArea_D + XyzEndArea_D)
       Area%Size_D  = 0.5*abs(XyzEndArea_D - XyzStartArea_D)

       ! Overwrite name with brick
       if(NameShape == "box_gen")then
          Area%NameShape = "brick_gen"
       else
          Area%NameShape = "brick"
       end if

    case("brick", "brick_gen")
       do iDim = 1, nDim
          call read_var("Size", Area%Size_D(iDim))
       end do

       ! Area size is measured from the center: half of brick size
       Area%Size_D = 0.5*Area%Size_D

    case("sphere", "shell")
       if(NameShape=="shell") call read_var("RadiusInner", Area%Radius1)
       call read_var("Radius", RadiusArea)
       Area%Size_D  = RadiusArea
       Area%Radius1 = Area%Radius1/RadiusArea
       if(DoStretch)then
          if(nDim>1) call read_var("RadiusY", Area%Size_D(Dim2_))
          if(nDim>2) call read_var("RadiusZ", Area%Size_D(Dim3_))
       end if

    case("cylinderx", "cylindery", "cylinderz", "ringx", "ringy", "ringz")
       call read_var("Length", Area%Size_D(iPar))
       Area%Size_D(iPar) = 0.5*Area%Size_D(iPar)
       if(NameShape(1:4)=="ring") call read_var("RadiusInner", Area%Radius1)
       call read_var("Radius", RadiusArea)
       Area%Size_D(iPerp_I) = RadiusArea
       Area%Radius1         = Area%Radius1/RadiusArea
       if(DoStretch .and. nDim==3) &
            call read_var("Radius2",Area%Size_D(iPerp_I(2)))

    case("conex", "coney", "conez", "funnelx", "funnely", "funnelz")
       call read_var("Height", Area%Size_D(iPar))
       if(NameShape(1:4) == "1.0")then
          ! A 1.0 is a funnel with 0 radius on one end
          Area%NameShape = "funnel" // NameShape(5:5)
       else
          call read_var("RadiusSmall", Area%Radius1)
       end if
       call read_var("Radius", RadiusArea)
       Area%Radius1 = Area%Radius1/RadiusArea    ! Normalized small radius
       Area%Size_D(iPerp_I) = RadiusArea
       if(DoStretch .and. nDim==3) &
            call read_var("RadiusPerp", Area%Size_D(iPerp_I(2)))

    case("doubleconex", "doubleconey", "doubleconez")
       call read_var("Height", Area%Size_D(iPar))
       Area%Size_D(iPar) = 0.5*Area%Size_D(iPar)
       call read_var("Radius", RadiusArea)
       Area%Size_D(iPerp_I) = RadiusArea
       if(DoStretch .and. nDim==3) &
            call read_var("RadiusPerp", Area%Size_D(iPerp_I(2)))

    case("paraboloidx", "paraboloidy", "paraboloidz")
       call read_var("Height", Area%Size_D(iPar))
       call read_var("Radius", RadiusArea)
       Area%Size_D(iPerp_I) = RadiusArea
       if(DoStretch .and. nDim==3) &
            call read_var("RadiusPerp", Area%Size_D(iPerp_I(2)))

    case default
       if(UseStrict) call CON_stop(NameSub// &
            ' ERROR: unknown NameShape='//trim(Area%NameShape))

       if(iProc == 0) &
            write(*,*) NameSub//' WARNING: unknown NameShape=' // &
            trim(Area%NameShape) // ', ignoring command ' // &
            trim(NameCommand)

       nArea = nArea - 1
       RETURN
    end select

    ! Read the tapering width if required
    if(DoTaperArea) call read_var('Taper', Area % Taper)

    ! Read rotation angles and create rotation matrix
    if(Area%DoRotate)then
       if(.not.allocated(Area%Rotate_DD))allocate(Area%Rotate_DD(nDim,nDim))

       if(nDim == 2)then
          call read_var('zRotate', zRotateArea)
          Area%Rotate_DD = rot_matrix(-zRotateArea*cDegToRad)
       elseif(nDim == 3)then
          ! Read 3 angles for the rotation matrix in degrees
          call read_var('xRotate', xRotateArea)
          call read_var('yRotate', yRotateArea)
          call read_var('zRotate', zRotateArea)

          ! Rotation matrix rotates around X, Y and Z axes in this order
          Area%Rotate_DD = matmul( matmul( &
               rot_matrix_z(-zRotateArea*cDegToRad), &
               rot_matrix_y(-yRotateArea*cDegToRad)),&
               rot_matrix_x(-xRotateArea*cDegToRad))
       end if

    end if

  end subroutine read_region_param
  !============================================================================

  subroutine block_inside_regions(iRegion_I, iBlock, nValue, StringLocation, &
       IsInside, IsInside_I, Value_I)

    use BATL_geometry, ONLY: IsCartesianGrid
    use ModUtilities,  ONLY: lower_case

    ! Interface for the external user routine
    interface
       subroutine user_specify_region(iArea, iBlock, nValue, NameLocation, &
            IsInside, IsInside_I, Value_I)

         implicit none

         integer,   intent(in):: iArea        ! area index in BATL_region
         integer,   intent(in):: iBlock       ! block index
         integer,   intent(in):: nValue       ! number of output values
         character, intent(in):: NameLocation ! c, g, x, y, z, or n

         logical, optional, intent(out) :: IsInside
         logical, optional, intent(out) :: IsInside_I(nValue)
         real,    optional, intent(out) :: Value_I(nValue)

       end subroutine user_specify_region
    end interface

    ! Check the intersection of block iBlock with one or more regions
    ! indexed by the iRegion_I array. Positive region index means
    ! an "OR", while negative region index means "AND NOT".
    ! If the first region index is negative, then it is relative to
    ! an "all" condition. The regions are evaluated in the order
    ! they are listed in the iRegion_I array.
    !
    ! If IsInside is present and neither IsInside_I/Value_I are present,
    ! then IsInside is set to true if the block intersects
    ! any of the + regions and completely avoids all the - regions.
    !
    ! If IsInside_I is present, it is set true for each point of the block
    ! (defined by StringLocation) that is inside any of the +region(s)
    ! and outside all of the -regions.
    ! If IsInside is also present, it is set to IsInside = any(IsInside_I).
    !
    ! If Value_I is present, then it is set to 1 inside any of the +region(s)
    ! and outside all of the -regions(s). In the tapering region the
    ! value gradually decreases to 0 and the rest is set to 0 value.
    ! This is done for each point of the block defined by StringLocation.
    ! If IsInside is also present, it is set to IsInside = any(Value_I > 0).
    !
    ! The StringLocation argument tells which points of the block are to be
    ! evaluated for IsInside_I and/or Value_I. Possible values are
    ! "cells", "ghosts", "xfaces", "yfaces", "zfaces", "faces", "nodes"
    ! Only the first character is significant and it is case insensitive.

    integer, intent(in):: iRegion_I(:)
    integer, intent(in):: iBlock
    integer, intent(in):: nValue

    character(*), optional, intent(in) :: StringLocation
    logical,      optional, intent(out):: IsInside
    logical,      optional, intent(out):: IsInside_I(nValue)
    real,         optional, intent(out):: Value_I(nValue)

    integer:: iRegion, nRegion, iArea, iSign
    integer:: nPoint
    character:: NameLocation

    logical:: DoBlock, DoMask, DoValue, DoBlockOnly

    logical             :: IsInsideOld

    logical:: DoSetCorner, DoSetCoord, DoSetXyz

    ! logical:: DoTest
    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'block_inside_regions'
    !--------------------------------------------------------------------------
    ! DoTest = iBlock == 136

    DoBlock  = present(IsInside)
    DoMask   = present(IsInside_I)
    DoValue  = present(Value_I)
    DoBlockOnly = DoBlock .and. .not.(DoMask .or. DoValue)

    if(.not.(DoBlock .or. DoMask .or. DoValue)) call CON_stop(NameSub// &
         ': no output argument is present')

    nRegion = size(iRegion_I)
    if(nRegion < 1) call CON_stop(NameSub//': empty region index array')

    if(DoTest)write(*,*) NameSub, &
         ' nRegion, iBlock, nValue, DoBlock, DoMask, DoValue=', &
         nRegion, iBlock, nValue, DoBlock, DoMask, DoValue

    ! Default number of points is the same as the number of returned values
    nPoint = nValue

    ! Allocate storage if multiple regions are used
    if(DoMask  .and. nRegion > 1)then
       if(allocated(IsInsideOld_I))then
          if(size(IsInsideOld_I) /= nPoint) deallocate(IsInsideOld_I)
       end if
       if(.not.allocated(IsInsideOld_I)) allocate(IsInsideOld_I(nPoint))
    end if
    if(DoValue .and. nRegion > 1)then
       if(allocated(ValueOld_I))then
          if(size(ValueOld_I) /= nPoint) deallocate(ValueOld_I)
       end if
       if(.not.allocated(ValueOld_I)) allocate(ValueOld_I(nPoint))
    end if

    ! Default location is nodes
    NameLocation = "n"
    if(present(StringLocation)) NameLocation = StringLocation(1:1)
    call lower_case(NameLocation)

    ! Check nodes of block by default (optimize for corners???)
    if(nPoint == 0) nPoint = nINode*nJNode*nKNode

    ! New block requires setting positions once
    DoSetCorner= .true.
    DoSetCoord = .true.
    DoSetXyz   = .true.

    if(DoTest)write(*,*) NameSub,' NameLocation, nPoint=', NameLocation, nPoint

    ! Evaluate all regions
    do iRegion = 1, nRegion

       ! Store results obtained from previous regions
       if(iRegion > 1)then
          if(DoBlock) IsInsideOld   = IsInside
          if(DoMask)  IsInsideOld_I = IsInside_I
          if(DoValue) ValueOld_I    = Value_I
       end if

       ! Get area index and evaluate it
       iArea = abs(iRegion_I(iRegion))
       iSign = sign(1, iRegion_I(iRegion))
       Area => Area_I(iArea)
       NameShape = Area%NameShape

       if(DoTest)write(*,*) NameSub,' iArea,iSign,NameShape=', &
            iArea, iSign, NameShape

       if(NameShape == 'all')then
          ! Set everything inside
          if(DoTest)write(*,*) NameSub,' set for all'
          if(DoBlock) IsInside = .true.
          if(DoMask)  IsInside_I = .true.
          if(DoValue) Value_I = 1.0
       elseif(NameShape(1:4) == 'user')then
          if(DoTest)write(*,*) NameSub,' call user_region'
          call user_specify_region(iArea, iBlock, nValue, NameLocation, &
               IsInside, IsInside_I, Value_I)
       elseif(DoBlockOnly .and. Area%IsSimple)then
          ! Use corners of the block to set IsInside
          if(DoTest)write(*,*) NameSub,' call is_block_inside'
          IsInside = is_block_inside(iBlock, Area, DoSetCorner)
       elseif(NameShape == 'brick_coord' .and. .not.IsCartesianGrid)then
          ! Use generalized coordinates for brick_coord
          if(DoTest)write(*,*) NameSub,' DoSetCoord=', DoSetCoord
          if(DoSetCoord) call set_coord
          call points_inside_region(nPoint, Coord_DI, Area, &
               IsInside, IsInside_I, Value_I)
       else
          ! The generic case uses Cartesian point  coordinates
          if(DoTest)write(*,*) NameSub,' DoSetXyz=', DoSetXyz
          if(DoSetXyz) call set_xyz
          call points_inside_region(nPoint, Xyz_DI, Area, &
               IsInside, IsInside_I, Value_I)
       end if

       if(DoTest.and.iRegion>1)then
          if(DoBlock)write(*,*) NameSub,' IsInsideOld, IsInside=', &
               IsInsideOld, IsInside
          if(DoValue)write(*,*) NameSub,' min, maxval ValueOld =', &
               minval(ValueOld_I), maxval(ValueOld_I)
          if(DoValue)write(*,*) NameSub,' min, maxval Value    =', &
               minval(Value_I), maxval(Value_I)
       end if
       ! Negate results for negative sign
       if(iRegion == 1)then
          if(iSign < 0)then
             if(DoBlock) IsInside   = .not. IsInside
             if(DoMask)  IsInside_I = .not. IsInside_I
             if(DoValue) Value_I    = 1 - Value_I
          end if
       else
          ! Combine last region info with previous
          if(iSign < 0)then
             if(DoBlock) IsInside   = IsInsideOld   .and. .not. IsInside
             if(DoMask)  IsInside_I = IsInsideOld_I .and. .not. IsInside_I
             if(DoValue) Value_I    = min(ValueOld_I, 1 - Value_I)
          else
             if(DoBlock) IsInside   = IsInsideOld   .or. IsInside
             if(DoMask)  IsInside_I = IsInsideOld_I .or. IsInside_I
             if(DoValue) Value_I    = max(ValueOld_I, Value_I)
          end if
       end if

       if(DoTest)then
          if(DoBlock)write(*,*) NameSub,' IsInside=', IsInside
          if(DoValue)write(*,*) NameSub,' min, max Value =', &
               minval(Value_I), maxval(Value_I)
       end if
    end do

    ! Correct IsInside in case an array of values is also returned
    if(DoBlock)then
       if(DoMask)then
          IsInside = any(IsInside_I)
       elseif(DoValue)then
          IsInside = any(Value_I > 0)
       end if
    endif

    if(DoTest)then
       write(*,*) NameSub,' finished with'
       if(DoBlock) write(*,*)' IsInside       =', IsInside
       if(DoMask)  write(*,*)' any(IsInside_I)=', any(IsInside_I)
       if(DoValue) write(*,*)' min, max Value =', &
            minval(Value_I), maxval(Value_I)
    end if

  contains
    !==========================================================================

    subroutine set_coord

      ! Set point positions in generalized coordinates
      use BATL_grid, ONLY: CoordMin_DB, CellSize_DB

      integer:: i, j, k, n
      real:: Coord_D(3), CellSize_D(3), CoordMinBlock_D(3)
      real:: CoordFace1, CoordFace2, CoordFace3
      real:: CoordCell1, CoordCell2, CoordCell3

      ! Allocate Coord array if new or size changed
      !------------------------------------------------------------------------
      if(allocated(Coord_DI))then
         if(size(Coord_DI) /= nDim*nPoint) deallocate(Coord_DI)
      end if
      if(.not.allocated(Coord_DI)) allocate(Coord_DI(nDim,nPoint))

      CellSize_D = CellSize_DB(:,iBlock)
      CoordMinBlock_D = CoordMin_DB(:,iBlock)

      n = 0
      select case(NameLocation)
      case('c')
         if(nPoint /= nIJK) call CON_stop(NameSub// &
              ': incorrect number of points for cell centers')
         do k = 1, nK
            Coord_D(3) = CoordMinBlock_D(3) + (k-0.5)*CellSize_D(3)
            do j = 1, nJ
               Coord_D(2) = CoordMinBlock_D(2) + (j-0.5)*CellSize_D(2)
               do i = 1, nI
                  Coord_D(1) = CoordMinBlock_D(1) + (i-0.5)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do
      case('g')
         if(nPoint /= MaxIJK) call CON_stop(NameSub// &
              ': incorrect number of points for cell centers with ghosts')

         do k = MinK, MaxK
            Coord_D(3) = CoordMinBlock_D(3) + (k-0.5)*CellSize_D(3)
            do j = MinJ, MaxJ
               Coord_D(2) = CoordMinBlock_D(2) + (j-0.5)*CellSize_D(2)
               do i = MinI, MaxI
                  Coord_D(1) = CoordMinBlock_D(1) + (i-0.5)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do
      case('x')
         if(nPoint /= (nI+1)*nJ*nK) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')
         do k = 1, nK
            Coord_D(3) = CoordMinBlock_D(3) + (k-0.5)*CellSize_D(3)
            do j = 1, nJ
               Coord_D(2) = CoordMinBlock_D(2) + (j-0.5)*CellSize_D(2)
               do i = 1, nI+1
                  Coord_D(1) = CoordMinBlock_D(1) + (i-1)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do

      case('y')
         if(nPoint /= nI*(nJ+1)*nK) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')

         do k = 1, nK
            Coord_D(3) = CoordMinBlock_D(3) + (k-0.5)*CellSize_D(3)
            do j = 1, nJ+1
               Coord_D(2) = CoordMinBlock_D(2) + (j-1)*CellSize_D(2)
               do i = 1, nI
                  Coord_D(1) = CoordMinBlock_D(1) + (i-0.5)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do

      case('z')
         if(nPoint /= nI*nJ*(nK+1)) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')

         do k = 1, nK+1
            Coord_D(3) = CoordMinBlock_D(3) + (k-1)*CellSize_D(3)
            do j = 1, nJ
               Coord_D(2) = CoordMinBlock_D(2) + (j-0.5)*CellSize_D(2)
               do i = 1, nI
                  Coord_D(1) = CoordMinBlock_D(1) + (i-0.5)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do

      case('f')
         if(nPoint /= nDim*nINode*nJNode*nKNode) call CON_stop(NameSub// &
              ': incorrect number of points for faces')
         n = 0
         do k = 1, nKNode
            CoordFace3 = CoordMinBlock_D(3) + (k-1)*CellSize_D(3)
            CoordCell3 = CoordMinBlock_D(3) + (min(k,nK)-0.5)*CellSize_D(3)
            do j = 1, nJNode
               CoordFace2 = CoordMinBlock_D(2) + (j-1)*CellSize_D(2)
               CoordCell2 = CoordMinBlock_D(2) + (min(j,nJ)-0.5)*CellSize_D(2)
               do i = 1, nINode
                  CoordFace1 = CoordMinBlock_D(1) + (i-1)*CellSize_D(1)
                  CoordCell1 = CoordMinBlock_D(1) + (min(i,nI)-0.5)*CellSize_D(1)

                  n = n + 1
                  Coord_D = (/ CoordFace1, CoordCell2, CoordCell3 /)
                  Coord_DI(:,n) = Coord_D(1:nDim)

                  if(nDim == 1) CYCLE
                  n = n + 1
                  Coord_D = (/ CoordCell1, CoordFace2, CoordCell3 /)
                  Coord_DI(:,n) = Coord_D(1:nDim)

                  if(nDim == 2) CYCLE
                  n = n + 1
                  Coord_D = (/ CoordCell1, CoordCell2, CoordFace3 /)
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do

      case('n')
         if(nPoint /= nINode*nJNode*nKNode) call CON_stop(NameSub// &
              ': incorrect number of points for nodes')

         do k = 1, nKNode
            Coord_D(3) = CoordMinBlock_D(3) + (k-1)*CellSize_D(3)
            do j = 1, nJNode
               Coord_D(2) = CoordMinBlock_D(2) + (j-1)*CellSize_D(2)
               do i = 1, nINode
                  Coord_D(1) = CoordMinBlock_D(1) + (i-1)*CellSize_D(1)
                  n = n + 1
                  Coord_DI(:,n) = Coord_D(1:nDim)
               end do
            end do
         end do

      case default
         write(*,*) NameSub,': StringLocation=', StringLocation
         call CON_stop(NameSub//' incorrect value for StringLocation')
      end select

      DoSetCoord = .false.

    end subroutine set_coord
    !==========================================================================
    subroutine set_xyz

      ! Set point positions in Cartesian coordinate

      use BATL_grid, ONLY: Xyz_DGB, Xyz_DNB

      integer:: i, j, k, iC, jC, kC, n

      ! Allocate Xyz array if new or size changed
      !------------------------------------------------------------------------
      if(allocated(Xyz_DI))then
         if(size(Xyz_DI) /= nDim*nPoint) deallocate(Xyz_DI)
      end if
      if(.not.allocated(Xyz_DI)) allocate(Xyz_DI(nDim,nPoint))

      select case(NameLocation)
      case('c')
         if(nPoint /= nIJK) call CON_stop(NameSub// &
              ': incorrect number of points for cell centers')

         Xyz_DI = reshape(Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock), &
              (/nDim, nPoint/))

      case('g')
         if(nPoint /= MaxIJK) call CON_stop(NameSub// &
              ': incorrect number of points for cell centers with ghosts')

         Xyz_DI = reshape(Xyz_DGB(1:nDim,:,:,:,iBlock), (/nDim, nPoint/))

      case('x')
         if(nPoint /= (nI+1)*nJ*nK) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')

         Xyz_DI = reshape( &
              0.5*(Xyz_DGB(1:nDim,0:nI  ,1:nJ,1:nK,iBlock) &
              +    Xyz_DGB(1:nDim,1:nI+1,1:nJ,1:nK,iBlock)), (/nDim, nPoint/))

      case('y')
         if(nPoint /= nI*(nJ+1)*nK) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')

         Xyz_DI = reshape( &
              0.5*(Xyz_DGB(1:nDim,1:nI,j0_:nJ, 1:nK,iBlock) &
              +    Xyz_DGB(1:nDim,1:nI,1:nJp1_,1:nK,iBlock)), (/nDim, nPoint/))

      case('z')
         if(nPoint /= nI*nJ*(nK+1)) call CON_stop(NameSub// &
              ': incorrect number of points for X faces')

         Xyz_DI = reshape( &
              0.5*(Xyz_DGB(1:nDim,1:nI,1:nJ,k0_:nK ,iBlock) &
              +    Xyz_DGB(1:nDim,1:nI,1:nJ,1:nKp1_,iBlock)), (/nDim, nPoint/))

      case('f')
         if(nPoint /= nDim*nINode*nJNode*nKNode) call CON_stop(NameSub// &
              ': incorrect number of points for faces')
         n = 0
         do k = 1, nKNode
            kC = min(k,nK)
            do j = 1, nJNode
               jC = min(j,nJ)
               do i = 1, nINode
                  iC = min(i,nI)
                  n = n + 1
                  Xyz_DI(:,n) = 0.5*(Xyz_DGB(1:nDim,i-1,jC,kC,iBlock) &
                       +             Xyz_DGB(1:nDim,i  ,jC,kC,iBlock))
                  if(nDim == 1) CYCLE
                  n = n + 1
                  Xyz_DI(:,n) = 0.5*(Xyz_DGB(1:nDim,iC,j-1,kC,iBlock) &
                       +             Xyz_DGB(1:nDim,iC,j  ,kC,iBlock))
                  if(nDim == 2) CYCLE
                  n = n + 1
                  Xyz_DI(:,n) = 0.5*(Xyz_DGB(1:nDim,iC,jC,k-1,iBlock) &
                       +             Xyz_DGB(1:nDim,iC,jC,k  ,iBlock))
               end do
            end do
         end do
      case('n')
         if(nPoint /= nINode*nJNode*nKNode) call CON_stop(NameSub// &
              ': incorrect number of points for nodes')

         Xyz_DI = reshape(Xyz_DNB(1:nDim,:,:,:,iBlock), (/nDim, nPoint/))

      case default
         write(*,*) NameSub,': StringLocation=', StringLocation
         call CON_stop(NameSub//' incorrect value for StringLocation')
      end select

      DoSetXyz = .false.

    end subroutine set_xyz
    !==========================================================================

  end subroutine block_inside_regions
  !============================================================================
  logical function is_block_inside(iBlock, Area, DoSetCorner)

    use ModNumConst, ONLY: cTiny
    use BATL_grid, ONLY: CoordMin_DB, CoordMax_DB

    ! Return true if block intersects the Area.
    ! This algorithm based on a very fast evaluation of the block corners
    ! for the symmetric shapes:
    ! brick, sphere, shell, cylinder, ring when they are not rotated.

    integer,        intent(in):: iBlock
    type(AreaType), intent(in):: Area
    logical, optional, intent(inout):: DoSetCorner

    logical:: DoSet
    integer, parameter :: nCorner = 2**nDim
    real     :: DistMin_D(nDim), DistMax_D(nDim), Radius1Sqr
    integer  :: iDim, iCorner

    logical, parameter:: DoTest = .false.

    character(len=*), parameter:: NameSub = 'is_block_inside'
    !--------------------------------------------------------------------------
    DoSet = .true.
    if(present(DoSetCorner))then
       DoSet = DoSetCorner
       DoSetCorner = .false.
    end if

    if(DoSet)then
       if(.not.allocated(Corner_DI)) &
            allocate(Corner_DI(nDim,nCorner), CornerNorm_DI(nDim,nCorner))

       ! Block corner coordinates
       Corner_DI(:,1)       = CoordMin_DB(1:nDim,iBlock)
       Corner_DI(:,nCorner) = CoordMax_DB(1:nDim,iBlock)
       if(nDim==2)then
          Corner_DI(:,2) = (/ CoordMax_DB(1,iBlock), CoordMin_DB(2,iBlock) /)
          Corner_DI(:,3) = (/ CoordMin_DB(1,iBlock), CoordMax_DB(2,iBlock) /)
       else if(nDim == 3)then
          Corner_DI(:,2) = (/ CoordMax_DB(1,iBlock), &
               CoordMin_DB(2,iBlock), CoordMin_DB(3,iBlock) /)
          Corner_DI(:,3) = (/ CoordMin_DB(1,iBlock), &
               CoordMax_DB(2,iBlock), CoordMin_DB(3,iBlock) /)
          Corner_DI(:,4) = (/ CoordMax_DB(1,iBlock), &
               CoordMax_DB(2,iBlock), CoordMin_DB(3,iBlock) /)
          Corner_DI(:,5) = (/ CoordMin_DB(1,iBlock), &
               CoordMin_DB(2,iBlock), CoordMax_DB(3,iBlock) /)
          Corner_DI(:,6) = (/ CoordMax_DB(1,iBlock), &
               CoordMin_DB(2,iBlock), CoordMax_DB(3,iBlock) /)
          Corner_DI(:,7) = (/ CoordMin_DB(1,iBlock), &
               CoordMax_DB(2,iBlock), CoordMax_DB(3,iBlock) /)
       end if
    end if

    ! Shift corner coordinates to the center of area and normalize to size
    do iCorner = 1, nCorner
       CornerNorm_DI(:,iCorner) = &
            (Corner_DI(:,iCorner) - Area%Center_D) / Area%Size_D
    end do

    ! Calculate maximum and minimum distances in all 3 directions
    ! Avoid rounding errors if possible
    do iDim = 1, nDim
       DistMax_D(iDim) = (1 - cTiny)*maxval(abs(CornerNorm_DI(iDim,:)))

       if(maxval(CornerNorm_DI(iDim,:))*minval(CornerNorm_DI(iDim,:)) <= 0)then
          ! The block covers the center point in this dimension
          DistMin_D(iDim) = 0
       else
          ! Select the point that is closer in this dimension
          DistMin_D(iDim) = (1 + cTiny)*minval(abs(CornerNorm_DI(iDim,:)))
       end if
    end do

    ! Set global shape name
    NameShape = Area%NameShape

    ! This occurs multiple times
    Radius1Sqr = (Area%Radius1)**2

    ! Check if this area is intersecting with the block
    call set_i_par_perp
    select case(NameShape)
    case('brick', 'brick_coord')
       is_block_inside = all(DistMin_D < 1)
    case('sphere')
       is_block_inside = sum(DistMin_D**2) < 1
    case('shell')
       ! Check if block intersects with the enclosing sphere
       ! but it is not fully inside the inner sphere
       is_block_inside = sum(DistMin_D**2) < 1 &
            .and. sum(DistMax_D**2) > Radius1Sqr
    case('cylinderx','cylindery','cylinderz')
       is_block_inside = DistMin_D(iPar) < 1 &
            .and. sum(DistMin_D(iPerp_I)**2) < 1
    case('ringx', 'ringy', 'ringz')
       ! Check if block intersects with the enclosing cylinder
       ! but it is not fully inside the inner cylinder
       is_block_inside = DistMin_D(iPar) < 1 &
            .and. sum(DistMin_D(iPerp_I)**2) < 1 &
            .and. sum(DistMax_D(iPerp_I)**2) > Radius1Sqr
    case default
       call CON_stop(NameSub //' ERROR: invalid NameShape = '//NameShape)
    end select

    if(DoTest)then
       write(*,*) NameSub,' iBlock, DistMin_D, inside=', &
            iBlock, DistMin_D, is_block_inside
    end if

  end function is_block_inside
  !============================================================================
  subroutine points_inside_region(&
       nPoint, Xyz_DI, Area, IsInside, IsInside_I, Value_I)

    use BATL_geometry, ONLY: IsCartesianGrid, IsPeriodic_D, DomainSize_D

    ! Check if the points listed in Xyz_DI are inside the Area

    ! The optional logical IsInside is set to true if any of the points
    ! are inside.
    !
    ! The optional logical array IsInside_I is set to true for points inside.
    !
    ! The optional real array Value_I is set to
    ! 0 if the point is fully outside the Area
    ! 1 if the point is fully inside the Area
    ! 0 to 1 if the point is in the taper region of the area (if any)

    integer,       intent(in)   :: nPoint
    real,          intent(inout):: Xyz_DI(nDim,nPoint) ! Gets normalized
    type(AreaType),intent(in)   :: Area
    logical, optional, intent(out):: IsInside
    logical, optional, intent(out):: IsInside_I(nPoint)
    real,    optional, intent(out):: Value_I(nPoint)

    integer  :: iPoint, iDim

    logical:: DoBlock, DoMask, DoValue, DoBlockOnly

    logical:: DoTaper, DoPeriodic
    real:: Xyz_D(nDim), Size_D(nDim)
    real:: Radius, RadiusSqr, Dist1, Dist2
    real:: Radius1, Radius1Sqr, Slope1
    real:: Taper, TaperFactor_D(nDim), TaperFactor1_D(nDim)
    real:: Slope_D(nDim)

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'points_inside_region'
    !--------------------------------------------------------------------------
    DoBlock  = present(IsInside)
    DoMask   = present(IsInside_I)
    DoValue  = present(Value_I)
    DoBlockOnly = DoBlock .and. .not. (DoMask .or. DoValue)

    if(DoTest) write(*,*) NameSub, &
         ': nPoint, DoBlock, DoMask, DoValue, DoBlockOnly=',&
         nPoint, DoBlock, DoMask, DoValue, DoBlockOnly

    NameShape = Area%NameShape

    if(DoTest) write(*,*) NameSub,' NameShape=', NameShape

    ! Treat special cases first
    if(NameShape == 'all')then
       if(DoBlock) IsInside   = .true.
       if(DoMask)  IsInside_I = .true.
       if(DoValue) Value_I    = 1.0
       RETURN
    end if

    ! Set iPar and iPerp_I indexes based on NameShape
    call set_i_par_perp

    ! Allocate normalized (to size of shape) coordinates
    if(allocated(Norm_DI))then
       if(size(Norm_DI) /= nDim*nPoint) deallocate(Norm_DI)
    endif
    if(.not. allocated(Norm_DI)) allocate(Norm_DI(nDim,nPoint))

    ! This occurs multiple times
    Radius1    = Area % Radius1
    Radius1Sqr = Radius1**2
    Slope1     = 1 - Radius1

    ! Use local variable for size
    Size_D = Area % Size_D

    ! The distance from the center will be normalized by Size_D,
    ! so the tapering length needs to be scaled accordingly
    Taper = Area % Taper
    DoTaper = DoValue .and. Taper > 0.0
    if(DoTaper)then
       Slope_D = abs(Size_D) / Taper
       TaperFactor_D = Slope_D/(Slope_D + 1)
       ! For the smaller radius of shell and ring
       if(Radius1 > 0.0) &
            TaperFactor1_D = Slope_D/max(Radius1*Slope_D - 1, 1e-30)
       ! Modified tapering slope for the sides of funnel, 1.0 and doublecone
       if(.not.allocated(SlopePerp_D)) allocate(SlopePerp_D(nDim-1))
       if(iPar > 0) SlopePerp_D = Slope_D(iPerp_I) &
            /sqrt(1 + (Slope1*Size_D(iPerp_I)/Size_D(iPar))**2)
    endif

    ! For now periodic boundaries are checked for
    ! Cartesian coordinates and for brick_coord.
    DoPeriodic = any(IsPeriodic_D(1:nDim)) .and. &
         (IsCartesianGrid .or. NameShape == 'brick_coord')

    ! Transform point coordinates for easy comparison
    do iPoint = 1, nPoint

       ! Shift to area center
       Xyz_D = Xyz_DI(:,iPoint) - Area % Center_D

       if(DoPeriodic)then
          do iDim = 1, nDim
             if(IsPeriodic_D(iDim))then
                if(Xyz_D(iDim) > 0.5*DomainSize_D(iDim)) &
                     Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)
                if(Xyz_D(iDim) < -0.5*DomainSize_D(iDim)) &
                     Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
             end if
          end do
       end if

       ! Rotate into area coordinates
       if(Area % DoRotate) Xyz_D = matmul(Area % Rotate_DD, Xyz_D)

       ! Rescale coordinates to the size of the area in all directions
       Norm_DI(:,iPoint) = Xyz_D / Size_D
    end do

    ! Initialize values for outside
    if(DoBlock) IsInside   = .false.
    if(DoMask)  IsInside_I = .false.
    if(DoValue) Value_I    = 0

    ! Set value for the shape
    select case(NameShape)
    case('brick', 'brick_coord')
       Norm_DI = abs(Norm_DI)
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(maxval(Norm_DI(:,iPoint)) >= 1) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - maxval( Slope_D*(Norm_DI(:,iPoint)-1) )))
          end do
       end if
    case('sphere')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(sum(Norm_DI(:,iPoint)**2) >= 1) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Dist1 = sum(Norm_DI(:,iPoint)**2)
             if(Dist1 < 1)then
                Value_I(iPoint) = 1.0 ! inside
                CYCLE
             end if
             Dist2 = sum((TaperFactor_D*Norm_DI(:,iPoint))**2)
             if(Dist2 >= 1) CYCLE     ! outside
             ! Use a roughly linear function between the ellipsoids
             Dist1 = sqrt(Dist1) - 1
             Dist2 = 1 - sqrt(Dist2)
             Value_I(iPoint) = Dist2/(Dist1 + Dist2)
          end do
       endif
    case('shell')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             Dist1 = sum(Norm_DI(:,iPoint)**2)
             if(Dist1 >= 1 .or. Dist1 <= Radius1Sqr)CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Dist1 = sum(Norm_DI(:,iPoint)**2)
             if(Dist1 >= 1)then
                Dist2 = sum((TaperFactor_D*Norm_DI(:,iPoint))**2)
                if(Dist2 >= 1) CYCLE     ! outside outer radius
                ! Use a roughly linear function between the ellipsoids
                Dist1 = sqrt(Dist1) - 1
                Dist2 = 1 - sqrt(Dist2)
                Value_I(iPoint) = Dist2/(Dist1 + Dist2)
             elseif(Dist1 <= Radius1Sqr)then
                Dist2 = sum((TaperFactor1_D*Norm_DI(:,iPoint))**2)
                if(Dist2 <= 1) CYCLE     ! inside inner radius
                Dist1 = 1 - sqrt(Dist1)/Radius1
                Dist2 = sqrt(Dist2) - 1
                Value_I(iPoint) = Dist2/(Dist1 + Dist2)
             else
                Value_I(iPoint) = 1.0 ! inside ring
             end if
          end do
       end if
    case('cylinderx', 'cylindery', 'cylinderz')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(abs(Norm_DI(iPar,iPoint)) >= 1) CYCLE
             if(sum(Norm_DI(iPerp_I,iPoint)**2) >= 1) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             ! Handle parallel direction first
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - Slope_D(iPar)*(abs(Norm_DI(iPar,iPoint))-1)))
             if(Value_I(iPoint) == 0) CYCLE          ! outside parallel to axis
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 <= 1) CYCLE                    ! inside perpendicular
             Dist2= sum((TaperFactor_D(iPerp_I)*Norm_DI(iPerp_I,iPoint))**2)
             if(Dist2 >= 1)then                      ! outside perpendicular
                Value_I(iPoint) = 0.0
                CYCLE
             end if
             ! Use a roughly linear function between the elliptic cylinders
             Dist1 = sqrt(Dist1) - 1
             Dist2 = 1 - sqrt(Dist2)
             Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
          end do
       end if
    case('ringx','ringy','ringz')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(abs(Norm_DI(iPar,iPoint)) > 1) CYCLE
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 >= 1 .or. Dist1 <= Radius1Sqr)CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - Slope_D(iPar)*(abs(Norm_DI(iPar,iPoint))-1)))
             if(Value_I(iPoint) == 0) CYCLE          ! outside parallel to axis
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 > 1)then
                Dist2= sum((TaperFactor_D(iPerp_I)*Norm_DI(iPerp_I,iPoint))**2)
                if(Dist2 >= 1)then                   ! outside outer radius
                   Value_I(iPoint) = 0.0
                   CYCLE
                end if
                ! Use a roughly linear function between the ellipsoids
                Dist1 = sqrt(Dist1) - 1
                Dist2 = 1 - sqrt(Dist2)
                Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
             elseif(Dist1 < Radius1Sqr)then
                Dist2=sum((TaperFactor1_D(iPerp_I)*Norm_DI(iPerp_I,iPoint))**2)
                if(Dist2 <= 1)then                   ! inside inner radius
                   Value_I(iPoint) = 0.0
                   CYCLE
                end if
                Dist1 = 1 - sqrt(Dist1)/Radius1
                Dist2 = sqrt(Dist2) - 1
                Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
             end if
          end do
       end if
    case('funnelx', 'funnely', 'funnelz')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(abs(Norm_DI(iPar,iPoint) - 0.5) >= 0.5) CYCLE
             Radius = Radius1 + Slope1*Norm_DI(iPar,iPoint)
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 >= Radius**2) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - Slope_D(iPar)*(abs(Norm_DI(iPar,iPoint)-0.5)-0.5)))
             if(Value_I(iPoint) == 0) CYCLE          ! outside parallel to axis
             ! Radius of funnel at this cross section (can be negative!)
             Radius = Radius1 + Slope1*Norm_DI(iPar,iPoint)
             Dist1  = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 <= max(0.0, Radius)**2) CYCLE  ! inside perpendicular
             ! Calculate distance from tapered funnel
             ! Use taper factor specific for this particular funnel radius
             Dist2 = sum( (SlopePerp_D/(Radius*SlopePerp_D + 1) &
                  *Norm_DI(iPerp_I,iPoint))**2 )
             if(Dist2 >= 1)then              ! outside perpendicular
                Value_I(iPoint) = 0.0
                CYCLE
             end if
             ! Use a roughly linear function between the elliptic cones
             Dist1 = sqrt(Dist1) - Radius
             ! Dist2 is multiplied by r+taper so it comparable to Dist1
             Dist2 = (1 - sqrt(Dist2))*(Radius + 1/SlopePerp_D(1))
             Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
          end do
       end if
    case('paraboloidx', 'paraboloidy', 'paraboloidz')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(abs(Norm_DI(iPar,iPoint)-0.5) >= 0.5) CYCLE ! outside parallel
             if(sum(Norm_DI(iPerp_I,iPoint)**2) &           ! outside
                  >= Norm_DI(iPar,iPoint)) CYCLE            !   perpendicular
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - Slope_D(iPar)*(abs(Norm_DI(iPar,iPoint)-0.5)-0.5)))
             if(Value_I(iPoint) == 0) CYCLE          ! outside parallel to axis
             Dist1  = sum(Norm_DI(iPerp_I,iPoint)**2)
             ! Radius squared of paraboloid is the parallel coordinate
             if(Dist1 <= Norm_DI(iPar,iPoint)) CYCLE ! inside perpendicular
             ! Calculate distance from tapered paraboloid
             RadiusSqr = (Slope_D(iPar)*Norm_DI(iPar,iPoint) + 1) &
                  / (Slope_D(iPar) + 1)
             Dist2 = sum((TaperFactor_D(iPerp_I)*Norm_DI(iPerp_I,iPoint))**2)
             if(Dist2 >= RadiusSqr)then              ! outside perpendicular
                Value_I(iPoint) = 0.0
                CYCLE
             end if
             ! Use a roughly linear function between the elliptic cones
             Dist1 = sqrt(Dist1) - sqrt(max(0.0, Norm_DI(iPar,iPoint)))
             Dist2 = sqrt(max(0.0, RadiusSqr)) - sqrt(Dist2)
             Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
          end do
       end if
    case('doubleconex', 'doubleconey', 'doubleconez')
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(abs(Norm_DI(iPar,iPoint)) >= 1) CYCLE ! outside parallel
             if(sum(Norm_DI(iPerp_I,iPoint)**2) &     ! outside perpendicular
                  >= Norm_DI(iPar,iPoint)**2) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Value_I(iPoint) = max(0.0, min(1.0, &
                  1 - Slope_D(iPar)*(abs(Norm_DI(iPar,iPoint)) - 1)))
             if(Value_I(iPoint) == 0) CYCLE
             ! Radius of double 1.0 at this cross section
             Radius = abs(Norm_DI(iPar,iPoint))
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 < Radius**2) CYCLE             ! inside perpendicular
             ! Use taper factor specific for this particular 1.0 radius
             Dist2 = sum( (SlopePerp_D/(Radius*SlopePerp_D + 1) &
                  *Norm_DI(iPerp_I,iPoint))**2 )
             if(Dist2 >= 1)then                    ! outside perpendicular
                Value_I(iPoint) = 0.0
                CYCLE
             end if
             ! Use a roughly linear function between the elliptic cylinders
             Dist1 = sqrt(Dist1) - Radius
             Dist2 = (1 - sqrt(Dist2))*(Radius + 1/SlopePerp_D(1))
             Value_I(iPoint) = Dist2/(Dist1 + Dist2) * Value_I(iPoint)
          end do
       end if

    case default
       if(NameShape(1:4) == 'user')then
          call CON_stop(NameSub// &
               ' ERROR: not yet implemented for NameShape='//trim(NameShape))
          ! call user_region(....)
       else
          call CON_stop(NameSub//' ERROR: unknown NameShape='//trim(NameShape))
       end if
    end select

    ! Get logical information if not yet calculated
    ! The tapered area counts as "inside".
    if(DoTaper)then
       if(DoMask)  IsInside_I =     Value_I > 0
       if(DoBlock) IsInside   = any(Value_I > 0)
    end if

    if(DoTest)then
       write(*,*) NameSub,' finished with'
       if(DoBlock) write(*,*)' IsInside       =', IsInside
       if(DoMask)  write(*,*)' any(IsInside_I)=', any(IsInside_I)
       if(DoValue) write(*,*)' maxval(Value_I)=', maxval(Value_I)
    end if

  end subroutine points_inside_region
  !============================================================================

end module BATL_region
!==============================================================================
