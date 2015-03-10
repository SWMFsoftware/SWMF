!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_amr_geometry

  use BATL_mpi, ONLY: iProc
  use BATL_size, ONLY: nI, nJ, nK

  implicit none

  SAVE

  private ! except

  public read_amr_geometry
  public init_amr_geometry
  public apply_amr_geometry
  public clean_amr_geometry

  ! compatible with BATSRUS way of doing AMR
  logical, public :: IsBatsrusAmr = .true. 

  ! Masking criteria for use
  ! This belongs more in batl_amr_criteria but also
  ! needed here so it is declared here.
  logical, public, allocatable :: UseCrit_IB(:,:)
  logical, public, allocatable :: UseCrit_B(:)

  ! We have only two values to decide the criteria
  ! for geometrick refinment cell size (1) , and AMR level (2)
  integer, parameter, public :: nGeoCrit = 2

  ! number of geometrical criteria actually used
  integer, public :: nCritGeoUsed = 0
  integer, public :: nCritGeoBackword = 0
  !Number of geometric criterias from AMRCRITERIALEVEL/RESOLUTION
  integer, public :: nCritDxLevel = 0

  ! we have read in new data from PARAM.in
  logical, public :: IsNewGeoParam =.false.

  integer, public, parameter :: MaxArea = 100, lNameArea = 20

  type, public :: AreaType 
     character(len=lNameArea) :: NameRegion
     character(len=lNameArea) :: Name
     real                     :: Resolution
     integer                  :: Level
     real, dimension(3)       :: Center_D,  Size_D
     real                     :: Radius1
     logical                  :: DoRotate
     real, dimension(3,3)     :: Rotate_DD
  end type AreaType

  ! index zero contains the 'all' area
  type(AreaType), public :: AreaGeo_I(0:MaxArea)

  ! Choosing with blocks we want to refine is based on a list of criteria 
  ! and a set of upper (refine)  and lower (coarsen) limits. The criteria 
  ! can be external or calculated internally by estimating the 
  ! numerical errors (calc_error_amr_criteria) based on the state variables.

  ! Storing All the AMR criteria
  integer, public          :: nAmrCrit = 0
  integer,public           :: nAmrCritUsed = 0
  real, public, allocatable:: AmrCrit_IB(:,:)

  real,    public, allocatable:: CoarsenCritAll_I(:), RefineCritAll_I(:)
  integer, public, allocatable:: iVarCritAll_I(:),iResolutionLimit_I(:)
  real,    public, allocatable:: ResolutionLimit_I(:)
  !type(AreaType), public, allocatable :: AreaAll_I(:)
  integer,    public, allocatable:: iAreaIdx_II(:,:)
  ! Storing names of areas for each criteria given by #AMRCRITERIA.....
  integer, public, allocatable :: nAreaPerCritAll_I(:)

  ! Local variables

  real    :: InitialResolution = -1.0
  integer :: initial_refine_levels = 0
  real    :: AreaResolution=0.0
  integer :: nLevelArea=0

contains
  !============================================================================
  subroutine init_amr_geometry

    use BATL_grid,         ONLY: CoordMin_D,CoordMax_D
    use BATL_geometry,     ONLY: x_
    use BATL_tree,         ONLY: nRoot_D

    integer :: iGeo, iVar,iGeoAll
    logical :: DoTestMe = .false.
    !-------------------------------------------------------------------------

    ! Fix resolutions (this depends on domain size set above)
    if(InitialResolution > 0.0) initial_refine_levels = nint( &
         alog(((CoordMax_D(x_) - CoordMin_D(x_)) / (nRoot_D(x_) * nI))  &
         / InitialResolution) / alog(2.0) )

    do iGeo = 1, nCritGeoUsed
       AreaResolution = AreaGeo_I(iGeo) % Resolution

       if(AreaResolution <= 0.0) then
          ! Convert back to integer area
          nLevelArea = ceiling(abs(AreaResolution))
          AreaGeo_I(iGeo) % Level = sign(nLevelArea,nint(AreaResolution))
          ! Set actual resolution
          AreaGeo_I(iGeo) % Resolution = (CoordMax_D(x_)-CoordMin_D(x_)) &
               / (nRoot_D(x_) * nI * 2.0**nLevelArea)
       else
          AreaGeo_I(iGeo) % Level = ceiling( &
               alog(((CoordMax_D(x_)-CoordMin_D(x_)) / (nRoot_D(x_) * nI))  &
               / AreaResolution) / alog(2.0) )
       end if

    end do

    iGeo =0
    do iGeoAll = 1, nCritGeoUsed
       iVar = 1
       !exclude areas given by #AREAREGION
       if( AreaGeo_I(iGeoAll) % NameRegion .ne. "NULL") CYCLE
       iGeo = iGeo+1
       if(AreaGeo_I(iGeoAll) % Level  < 0) then
          RefineCritAll_I(nAmrCritUsed+iGeo)  = AreaGeo_I(iGeoAll)%Level
          CoarsenCritAll_I(nAmrCritUsed+iGeo) = &
               RefineCritAll_I(nAmrCritUsed+iGeo)-1
          iVar = 2
          ResolutionLimit_I(nAmrCritUsed+iGeo) = &
               -abs(AreaGeo_I(iGeo)%Level)
       else
          RefineCritAll_I(nAmrCritUsed+iGeo)  = AreaGeo_I(iGeoAll)%Resolution 
          CoarsenCritAll_I(nAmrCritUsed+iGeo) = AreaGeo_I(iGeoAll)%Resolution/2.0
          iVar = 1
          ResolutionLimit_I(nAmrCritUsed+iGeo) = &
               abs(AreaGeo_I(iGeo)%Resolution)
       end if

       ! All geometrical criteias is a comparson to gird values
       if(nAmrCritUsed  > 0) then
          iVarCritAll_I(nAmrCritUsed+iGeo) = &
               maxval(iVarCritAll_I(1:nAmrCritUsed-nCritDxLevel))+iVar
          iResolutionLimit_I(nAmrCritUsed+iGeo) = &
               maxval(iVarCritAll_I(1:nAmrCritUsed-nCritDxLevel))+iVar
       else
          iVarCritAll_I(iGeo) = iVar   
          iResolutionLimit_I(nAmrCritUsed+iGeo) = iVar
       end if

    end do


    ! Add geometric info for BATSRUS amr type of params
    iGeo = 1
    do iGeoAll = 1, nCritGeoUsed
       if( AreaGeo_I(iGeoAll) % NameRegion .ne. "NULL") CYCLE
       iAreaIdx_II(1,nAmrCritUsed+iGeo) = iGeoAll
       nAreaPerCritAll_I(nAmrCritUsed+iGeo) = 1
       iGeo = iGeo+1
    end do


    if(DoTestMe .and. iProc == 0) then
       write(*,*) "START init_amr_geometry "
       write(*,*) "nCritGeoUsed     = ",nCritGeoUsed
       do iGeoAll = 0, nCritGeoUsed
          write(*,*) "--------------------"
          write(*,*) "Area namesRegion :: ", AreaGeo_I(iGeoAll)%NameRegion
          write(*,*) "Area names       :: ", AreaGeo_I(iGeoAll)%Name
          write(*,*) "Area resolution  :: ", AreaGeo_I(iGeoAll)%Resolution
          write(*,*) "Area level       :: ", AreaGeo_I(iGeoAll)%Level
          write(*,*) "Area Center_D    :: ", AreaGeo_I(iGeoAll)%Center_D(1:3)
          write(*,*) "Area Size_D      :: ", AreaGeo_I(iGeoAll)%Size_D(1:3)
          write(*,*) "Area Radius1     :: ", AreaGeo_I(iGeoAll)%Radius1
          write(*,*) "Area DoRotate    :: ", AreaGeo_I(iGeoAll)%DoRotate
          write(*,*) "Area Rotate_DD   :: ", AreaGeo_I(iGeoAll)%Rotate_DD(1:3,1:3)
          write(*,*) "--------------------"
       end do
       write(*,*) "END init_amr_geometry "
    end if

    IsNewGeoParam =.false.

  end subroutine init_amr_geometry

  !============================================================================
  subroutine apply_amr_geometry(iBlock, Area, UseBlock,DoCalcCritIn,&
       user_amr_geometry)

    !DESCRIPTION:
    ! Set DoRefine_B to .true. for blocks touching the predefined areas
    ! if the area has a finer resolution than the block

    ! Set UseCrit_IB to true  for block touching the predefined areas
    ! if the area has a finer resolution than the block
    ! Also return max cell size for the block in CritGeo_IB

    ! WARNING
    ! user_specify_refinement,     not suported
    ! currentsheet,                not suported

    use BATL_geometry, ONLY: gen_to_radius,r_, Phi_, Theta_,TypeGeometry,&
         IsCartesianGrid
    use BATL_grid,     ONLY: Xyz_DNB
    use BATL_size,     ONLY: nINode,nJNode,nKNode,nDim
    use ModNumConst,   ONLY: cTiny, cRadToDeg
    use BATL_grid,     ONLY: CoordMin_DB, CoordMax_DB
    implicit none

    interface
       subroutine user_amr_geometry(iBlock, iArea, DoRefine)
         integer, intent(in) :: iBlock, iArea
         logical,intent(out) :: DoRefine
       end subroutine user_amr_geometry
    end interface
    optional :: user_amr_geometry

    integer, intent(in)       :: iBlock
    type(AreaType),intent(in) :: Area
    logical, intent(out)      :: UseBlock
    logical, intent(in), optional :: DoCalcCritIn

    !LOCAL VARIABLES:
    character(len=lNameArea) :: NameArea

    logical :: DoRefine, IsSpecialArea
    integer, parameter :: nCorner = 8
    real     :: Corner_DI(3, nCorner)
    real     :: DistMin_D(3), DistMax_D(3), Radius1Sqr
    integer  :: i, j, k, iDim, iCorner

    ! These variables are needed for generalized coordinates
    real :: Xyz_D(3)
    real, dimension(nINode,nJNode,nKNode):: x_N, y_N, z_N, R2_N

    character(len=*), parameter:: NameSub = 'apply_amr_geometry'
    logical, parameter:: DoTestBlock = .false.

    logical :: DoCalcCrit
    !--------------------------------------------------------------------------
    !if(nCritGeoUsed <= 0) RETURN

    DoCalcCrit = .true.
    if(present(DoCalcCritIn)) DoCalcCrit =  DoCalcCritIn


    if(DoCalcCrit) call calc_crit(iBlock)

    UseBlock = .false.

    ! Blocks outer corners
    Corner_DI(:,1) = CoordMin_DB(:,iBlock)
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
    Corner_DI(:,8) = CoordMax_DB(:,iBlock)

    NameArea = Area % Name

    if(DoTestBlock) write(*,*) NameSub,'Name,Resolution=',&
         ' ',trim(NameArea),' ',Area % Resolution

    ! Treat special cases first
    IsSpecialArea = .true.
    select case(NameArea)
    case('all','initial')
       UseBlock = .true.
    case('usergeo')
       if(present(user_amr_geometry)) then
          call user_amr_geometry(iBlock, -1, UseBlock)
       else
          call CON_stop(NameSub //' ERROR: need user_amr_geometry')
       end if
    case default
       IsSpecialArea = .false.
    end select

    if(IsSpecialArea) RETURN

    ! Check if it is a brick defined in generalized coordinates
    if(NameArea == "brick_gen" .and. .not.IsCartesianGrid)then

!!! TO BE IMPROVED !!!

       ! Convert angles to degrees and ln(r)/genr to r
       ! Make sure that phi=360 is not 0 but really 360
       select case(TypeGeometry)
       case('rlonlat')
          Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
          Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg

       case('rlonlat_lnr')
          Corner_DI(r_,:)     = exp(Corner_DI(r_,:))
          Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
          Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg

       case('rlonlat_genr')
          do iCorner = 1, nCorner
             call gen_to_radius(Corner_DI(r_,iCorner))
          end do
          if(DoTestBlock) write(*,*) "gen_to_radius  Corner_DI = ",Corner_DI
          Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
          Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg

       case('cylindrical')
          Corner_DI(Phi_,:) = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
       end select
    end if

    if(DoTestBlock) write(*,*) "post Corner_DI = ",Corner_DI

    ! Shift corner coordinates to the center of area
    do iCorner = 1, nCorner
       Corner_DI(1:nDim,iCorner) = &
            Corner_DI(1:nDim,iCorner) - Area % Center_D(1:nDim)
    end do

    ! Rotate corners into the orientation of the area if required
    if(Area % DoRotate) Corner_DI = matmul(Area % Rotate_DD, Corner_DI)

    ! Normalize coordinates to the size of the area in all 3 directions
    do iCorner = 1, nCorner
       Corner_DI(:,iCorner) = Corner_DI(:,iCorner) / Area % Size_D
    end do

    ! Calculate maximum and minimum distances in all 3 directions
    ! Avoid rounding errors if possible
    do iDim = 1, 3
       DistMax_D(iDim) = (1-cTiny)*maxval(abs(Corner_DI(iDim,:)))

       if( maxval(Corner_DI(iDim,:))*minval(Corner_DI(iDim,:)) <= 0.0)then
          ! The block covers the center point in this dimension
          DistMin_D(iDim) = 0.0
       else
          ! Select the point that is closer in this dimension
          DistMin_D(iDim) = (1+cTiny)*minval(abs(Corner_DI(iDim,:)))
       end if
    end do

    ! This occurs multiple times
    Radius1Sqr = (Area % Radius1)**2

    if(DoTestBlock)then
       write(*,*) NameSub,' DistMin_D=',DistMin_D
       write(*,*) NameSub,' DistMax_D=',DistMax_D
    end if

    ! Check if this area is intersecting with the block
    select case( NameArea)
    case('brick', 'brick_gen')
       DoRefine = all( DistMin_D < 1.0 )
    case('sphere')
       DoRefine = sum(DistMin_D**2) < 1.0
    case('shell')
       ! Check if block intersects with the enclosing sphere
       ! but it is not fully inside the inner sphere
       DoRefine = sum(DistMin_D**2)<1.0 .and. sum(DistMax_D**2)>Radius1Sqr
    case('cylinderx')
       DoRefine = DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0
    case('cylindery')
       DoRefine = DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0
    case('cylinderz')
       DoRefine = DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0
    case('ringx')
       ! Check if block intersects with the enclosing cylinder
       ! but it is not fully inside the inner cylinder
       DoRefine = DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0 &
            .and. sum(DistMax_D(2:3)**2) > Radius1Sqr 
    case('ringy')
       DoRefine = DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0 &
            .and. sum(DistMax_D(1:3:2)**2) > Radius1Sqr
    case('ringz')
       DoRefine = DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0 &
            .and. sum(DistMax_D(1:2)**2) > Radius1Sqr 

    case default
       call CON_stop(NameSub //' ERROR: Unknown NameArea = '//NameArea)

    end select

    if(DoTestBlock)write(*,*) NameSub,' DoRefine (from corners)=',DoRefine

    if(DoRefine) UseBlock = .true.
    
    if(NameArea == 'brick_gen' .or. IsCartesianGrid) RETURN

    ! Non-cartesian case
    UseBlock = .true.

    ! Check all nodes of the block
    do k=1,nKNode; do j=1,nJNode; do i=1,nINode

       ! Shift to area center
       Xyz_D = Xyz_DNB(:,i,j,k,iBlock) - Area % Center_D

       ! Rotate into area coordinates
       if(Area % DoRotate) &
            Xyz_D = matmul(Area % Rotate_DD, Xyz_D)

       ! Rescale coordinates to the size of the area in all directions
       Xyz_D = Xyz_D / Area % Size_D

       ! We only need the absolute values of the coordinates
       x_N(i,j,k) = abs(Xyz_D(1))
       y_N(i,j,k) = abs(Xyz_D(2))
       z_N(i,j,k) = abs(Xyz_D(3))
    end do; end do; end do

    select case( NameArea)
    case('brick')
       if( any( x_N<1.0 .and. y_N<1.0 .and. z_N<1.0) ) RETURN
    case('sphere')
       if( any(x_N**2 + y_N**2 + z_N**2 < 1.0) ) RETURN
    case('shell')
       R2_N = x_N**2 + y_N**2 + z_N**2
       if( any(R2_N < 1.0 .and. R2_N > Radius1Sqr )) RETURN
    case('cylinderx')
       if( any(x_N < 1.0 .and. y_N**2+z_N**2 < 1.0 ) ) RETURN
    case('cylindery')
       if( any(y_N < 1.0 .and. x_N**2+z_N**2 < 1.0 ) ) RETURN
    case('cylinderz')
       if( any(z_N < 1.0 .and. x_N**2+y_N**2 < 1.0 ) ) RETURN
    case('ringx')
       R2_N = y_N**2+z_N**2
       if(any(x_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) RETURN
    case('ringy')
       R2_N = x_N**2+z_N**2
       if(any(y_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) RETURN
    case('ringz')
       R2_N = x_N**2+y_N**2
       if(any(z_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) RETURN
    case default
       call CON_stop(NameSub //' ERROR: Unknown NameArea = '//NameArea)
    end select

    ! The block is not inside
    UseBlock = .false.

    if(DoTestBlock) write(*,*)NameSub,' DoRefine=false'

    if(DoTestBlock) write(*,*)NameSub,' DoRefine final=',UseBlock


  end subroutine apply_amr_geometry

  !============================================================================
  subroutine calc_crit(iBlock)

    use BATL_grid,     ONLY: Xyz_DNB,CellSize_DB
    use BATL_size,     ONLY: nINode,nJNode,nKNode,nDim
    use BATL_geometry, ONLY: IsCartesian
    use BATL_tree,     ONLY: iNode_B, iTree_IA, Level_

    integer, intent(in) :: iBlock

    real    :: MaxLength
    integer :: i,j,k
    !-------------------------------------------------------------------------

    if(IsCartesian) then
       MaxLength = maxval(CellSize_DB(1:nDim,iBlock))
    else
       ! Find the longest cell edge in the block
       MaxLength = 0.0

       ! Get maximum length squared of i-edges
       do k = 1, nKNode; do j = 1, nJNode; do i = 2, nINode
          MaxLength = max(MaxLength,&
               sum((Xyz_DNB(:,i  ,j,k,iBlock) &
               -    Xyz_DNB(:,i-1,j,k,iBlock))**2))
       end do; end do ; end do

       if(nDim >1) then
          ! Get maximum length squared of j-edges
          do k = 1, nKNode; do j = 2, nJNode; do i = 1, nINode
             MaxLength = max(MaxLength, &
                  sum((Xyz_DNB(:,i,j  ,k,iBlock) &
                  -    Xyz_DNB(:,i,j-1,k,iBlock))**2))
          end do; end do ; end do
       end if

       if(nDim >2) then
          ! Get maximum length squared of k-edges
          do k = 2, nKNode; do j = 1, nJNode; do i = 1, nINode
             MaxLength = max(MaxLength,               &
                  sum((Xyz_DNB(:,i,j,k  ,iBlock)      &
                  -    Xyz_DNB(:,i,j,k-1,iBlock))**2))
          end do; end do ; end do
       end if

       ! Get maximum length
       MaxLength = sqrt(MaxLength)

    end if

    AmrCrit_IB(nAmrCrit-nGeoCrit+1:nAmrCrit,iBlock) = &
         (/ MaxLength, -real(iTree_IA(Level_,iNode_B(iBlock))) /)

  end subroutine calc_crit

  !============================================================================
  subroutine read_amr_geometry(NameCommand, UseStrictIn, &
       InitLevelInOut, InitResInOut)

    use ModReadParam ,     ONLY: read_var, lStringLine
    use ModNumConst,       ONLY: cUnit_DD,cDegToRad
    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z

    character(len=*),  intent(inout) :: NameCommand
    logical, optional, intent(in)    :: UseStrictIn
    integer, optional, intent(inout) :: InitLevelInOut 
    real,    optional, intent(inout) :: InitResInOut 

    character(len=lStringLine):: NameArea='all'
    real    :: RadiusArea=0.0
    logical :: DoReadAreaCenter=.false.
    real    :: XyzStartArea_D(3)=0.0, XyzEndArea_D(3)=0.0
    real    :: xRotateArea=0., yRotateArea=0., zRotateArea=0.
    logical :: DoStretchArea = .false.
    real    :: xStretchArea=1.0, yStretchArea=1.0, zStretchArea=1.0
    character(len=lNameArea) :: NameRegion
    integer :: i
    logical :: UseStrict

    character (len=17) :: NameSub='read_amr_geometry'
    !-------------------------------------------------------------------------

    isNewGeoParam =.true.

    UseStrict = .true.
    if(present(UseStrictIn)) UseStrict = UseStrictIn
    NameRegion ="NULL"
    if(NameCommand == "#AMRREGION") then
       call read_var('NameRegion',NameRegion, IsLowerCase=.true.)
       AreaResolution = 0
    else
       ! #GRIDLEVEL/RESOLUTION
       if(index(NameCommand,"RESOLUTION") > 0)then
          call read_var('AreaResolution', AreaResolution) 
       else
          call read_var('nLevelArea', nLevelArea)
          ! Store level as a negative integer resolution.
          AreaResolution = -nLevelArea
       end if
    end if
    call read_var('TypeRegion', NameArea, IsLowerCase=.true.)

    ! Remove leading spaces
    NameArea = adjustl(NameArea)

    if(NameArea(1:4) == 'init')then
       ! 'init' or 'initial' means that the initial resolution is set,
       ! and no area is created. 
       if(AreaResolution > 0)then
          InitialResolution = AreaResolution
       else
          initial_refine_levels = nLevelArea
       end if
       if(present(InitLevelInOut)) InitLevelInOut = initial_refine_levels
       if(present(InitResInOut)) InitResInOut = AreaResolution
       RETURN
    end if

    nCritGeoUsed = nCritGeoUsed + 1
    if(NameRegion == "NULL") nCritGeoBackword = nCritGeoBackword + 1

    if(nCritGeoUsed > MaxArea)then
       if(UseStrict) &
            call CON_stop(NameSub//' ERROR: Too many grid areas were defined')
       if(iProc == 0)then
          write(*,*) NameSub," nCritGeoUsed = ",nCritGeoUsed
          write(*,*) NameSub," WARNING: Too many grid areas were defined"
          write(*,*) NameSub," ignoring command ",NameCommand
       end if
       nCritGeoUsed = MaxArea
    end if

    AreaGeo_I(nCritGeoUsed) % Resolution = AreaResolution
    AreaGeo_I(nCritGeoUsed) % NameRegion = NameRegion
    ! Set the default center to be the origin, 
    ! the size and radii to be 1, and no rotation
    AreaGeo_I(nCritGeoUsed)%Center_D  = 0.0
    AreaGeo_I(nCritGeoUsed)%Size_D    = 1.0
    AreaGeo_I(nCritGeoUsed)%Radius1   = 1.0
    AreaGeo_I(nCritGeoUsed)%DoRotate  = .false.
    AreaGeo_I(nCritGeoUsed)%Rotate_DD = cUnit_DD

    ! Check for the word rotated in the name
    i = index(NameArea,'rotated')
    AreaGeo_I(nCritGeoUsed)%DoRotate = i > 0
    if(i>0) NameArea = NameArea(1:i-1)//NameArea(i+7:len(NameArea))

    ! check for the word stretched in the name
    i = index(NameArea,'stretched')
    DoStretchArea = i > 0
    if(i>0) NameArea = NameArea(1:i-1)//NameArea(i+9:len(NameArea))

    ! Extract character '0' from the name
    i = index(NameArea,'0')
    if(i>0) NameArea = NameArea(1:i-1)//NameArea(i+1:len(NameArea))

    DoReadAreaCenter = (i < 1 .and. NameArea(1:3) /= 'box')

    ! Store name
    AreaGeo_I(nCritGeoUsed)%Name = NameArea

    ! These types do not need any more parameters
    select case(NameArea)
    case('all')
       return
    case("user")
       NameArea ='usergeo'
       AreaGeo_I(nCritGeoUsed)%Name = NameArea
       return
    case('currentsheet')
       if(iProc == 0) then
          write(*,*) "OLD PARAM.in FILE!!!!!"
          write(*,*) ""
          write(*,*) &
               "Use one of the following commands for current sheet refinement:"
          write(*,*) ""
          write(*,*) ""
          write(*,*) "#AMRCRITERIARESOLUTION"
          write(*,*) "1                       nCriteria "
          write(*,*) "currentsheet            TypeCriteria"
          write(*,*) "0.5                     CoarsenLimit"
          write(*,*) "0.5                     RefineLimit"
          write(*,*) "0.2                     MinResolution"
          write(*,*) ""
          write(*,*) "or"
          write(*,*) ""
          write(*,*) "#AMRCRITERIALEVEL"
          write(*,*) "1                       nCriteria"
          write(*,*) "currentsheet            TypeCriteria"
          write(*,*) "0.5                     CoarsenLimit"
          write(*,*) "0.5                     RefineLimit"
          write(*,*) "5                       MaxLevel"
          write(*,*) ""
          write(*,*) &
               "with the desired level/resolution values and additional criteria" 
       end if
       call CON_stop(NameSub//': fix PARAM.in file')
    end select

    ! Read center of area if needed
    if(DoReadAreaCenter)then
       call read_var("xCenter",AreaGeo_I(nCritGeoUsed)%Center_D(1))
       call read_var("yCenter",AreaGeo_I(nCritGeoUsed)%Center_D(2))
       call read_var("zCenter",AreaGeo_I(nCritGeoUsed)%Center_D(3))
    endif

    select case(NameArea)
    case("box", "box_gen")
       call read_var("xMinBox",XyzStartArea_D(1))
       call read_var("yMinBox",XyzStartArea_D(2))
       call read_var("zMinBox",XyzStartArea_D(3))
       call read_var("xMaxBox",XyzEndArea_D(1))
       call read_var("yMaxBox",XyzEndArea_D(2))
       call read_var("zMaxBox",XyzEndArea_D(3))
       ! Convert to center and size information
       AreaGeo_I(nCritGeoUsed)%Center_D= 0.5*   (XyzStartArea_D + XyzEndArea_D)
       AreaGeo_I(nCritGeoUsed)%Size_D  = 0.5*abs(XyzEndArea_D - XyzStartArea_D)

       ! Overwrite name with brick
       if(NameArea == "box_gen")then
          AreaGeo_I(nCritGeoUsed)%Name = "brick_gen"
       else
          AreaGeo_I(nCritGeoUsed)%Name = "brick"
       end if
    case("brick", "brick_gen")
       call read_var("xSize", AreaGeo_I(nCritGeoUsed)%Size_D(1))
       call read_var("ySize", AreaGeo_I(nCritGeoUsed)%Size_D(2))
       call read_var("zSize", AreaGeo_I(nCritGeoUsed)%Size_D(3))

       ! Area size is measured from the center: half of brick size
       AreaGeo_I(nCritGeoUsed)%Size_D   = 0.5*AreaGeo_I(nCritGeoUsed)%Size_D

    case("shell")
       call read_area_radii
       AreaGeo_I(nCritGeoUsed)%Size_D = RadiusArea

    case("sphere")
       call read_var("Radius", RadiusArea)
       AreaGeo_I(nCritGeoUsed)%Size_D   = RadiusArea

    case("ringx")
       call read_var("Height", AreaGeo_I(nCritGeoUsed)%Size_D(1))
       AreaGeo_I(nCritGeoUsed)%Size_D(1)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(1)
       call read_area_radii
       AreaGeo_I(nCritGeoUsed)%Size_D(2:3)   = RadiusArea

    case("ringy")
       call read_var("Height", AreaGeo_I(nCritGeoUsed)%Size_D(2))
       AreaGeo_I(nCritGeoUsed)%Size_D(2)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(2)
       call read_area_radii
       AreaGeo_I(nCritGeoUsed)%Size_D(1:3:2) = RadiusArea

    case("ringz")
       call read_var("Height", AreaGeo_I(nCritGeoUsed)%Size_D(3))
       AreaGeo_I(nCritGeoUsed)%Size_D(3)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(3)
       call read_area_radii
       AreaGeo_I(nCritGeoUsed)%Size_D(1:2)   = RadiusArea

    case("cylinderx")
       call read_var("Length", AreaGeo_I(nCritGeoUsed)%Size_D(1))
       AreaGeo_I(nCritGeoUsed)%Size_D(1)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(1)
       call read_var("Radius", RadiusArea)
       AreaGeo_I(nCritGeoUsed)%Size_D(2:3)   = RadiusArea

    case("cylindery")
       call read_var("Length", AreaGeo_I(nCritGeoUsed)%Size_D(2))
       AreaGeo_I(nCritGeoUsed)%Size_D(2)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(2)
       call read_var("Radius", RadiusArea)
       AreaGeo_I(nCritGeoUsed)%Size_D(1:3:2) = RadiusArea

    case("cylinderz")
       call read_var("Length", AreaGeo_I(nCritGeoUsed)%Size_D(3))
       AreaGeo_I(nCritGeoUsed)%Size_D(3)     = &
            0.5*AreaGeo_I(nCritGeoUsed)%Size_D(3)
       call read_var("Radius", RadiusArea)
       AreaGeo_I(nCritGeoUsed)%Size_D(1:2)   = RadiusArea

    case default
       if(UseStrict) call CON_stop(NameSub// &
            ' ERROR: unknown Name='//trim(AreaGeo_I(nCritGeoUsed)%Name))

       if(iProc == 0) &
            write(*,*) NameSub//' WARNING: unknown NameArea=' // &
            trim(AreaGeo_I(nCritGeoUsed)%Name) // ', ignoring command ' // &
            trim(NameCommand)

       nCritGeoUsed = nCritGeoUsed - 1
       RETURN
    end select

    if(DoStretchArea)then
       ! Read 3 stretch factors for the sizes
       call read_var('xStretch', xStretchArea)
       call read_var('yStretch', yStretchArea)
       call read_var('zStretch', zStretchArea)

       ! Stretch the x, y, z sizes
       AreaGeo_I(nCritGeoUsed)%Size_D(1) = &
            AreaGeo_I(nCritGeoUsed)%Size_D(1)*xStretchArea
       AreaGeo_I(nCritGeoUsed)%Size_D(2) = &
            AreaGeo_I(nCritGeoUsed)%Size_D(2)*yStretchArea
       AreaGeo_I(nCritGeoUsed)%Size_D(3) = &
            AreaGeo_I(nCritGeoUsed)%Size_D(3)*zStretchArea

       DoStretchArea = .false.
    end if

    if(AreaGeo_I(nCritGeoUsed) % DoRotate)then

       ! Read 3 angles for the rotation matrix in degrees
       call read_var('xRotate', xRotateArea)
       call read_var('yRotate', yRotateArea)
       call read_var('zRotate', zRotateArea)

       ! Rotation matrix rotates around X, Y and Z axes in this order
       AreaGeo_I(nCritGeoUsed)%Rotate_DD = matmul( matmul( &
            rot_matrix_z(-zRotateArea*cDegToRad), &
            rot_matrix_y(-yRotateArea*cDegToRad)),&
            rot_matrix_x(-xRotateArea*cDegToRad))

    end if

  contains
    !==========================================================================
    subroutine read_area_radii

      real :: Radius1, Radius2

      ! Read inner and outer radii for areas shell and ring
      call read_var("Radius1", Radius1)
      call read_var("Radius2", Radius2)

      ! Set outer size of the area
      RadiusArea = max(Radius1, Radius2)
      ! Normalize inner radius to the outer size
      AreaGeo_I(nCritGeoUsed)%Radius1  = min(Radius1, Radius2) / RadiusArea

    end subroutine read_area_radii

  end subroutine read_amr_geometry
  !===========================================================================
  subroutine clean_amr_geometry

    if(allocated(UseCrit_IB)) deallocate(UseCrit_IB)
    if(allocated(UseCrit_B)) deallocate(UseCrit_B)
    if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
    if(allocated(CoarsenCritAll_I)) deallocate(CoarsenCritAll_I)
    if(allocated(RefineCritAll_I)) deallocate(RefineCritAll_I)
    if(allocated(ResolutionLimit_I)) deallocate(ResolutionLimit_I)
    if(allocated(iVarCritAll_I)) deallocate(iVarCritAll_I)
    !if(allocated(AreaAll_I)) deallocate(AreaAll_I)

    !nCritGeoUsed          = 0
    IsNewGeoParam         = .true.
    InitialResolution     = -1.0 
    initial_refine_levels = 0
    AreaResolution        = 0.0
    nAmrCrit              = 0
    nAmrCritUsed          = 0

  end subroutine clean_amr_geometry

end module BATL_amr_geometry
