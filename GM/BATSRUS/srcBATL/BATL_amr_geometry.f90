!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_amr_geometry

  use BATL_mpi,  ONLY: iProc
  use BATL_size, ONLY: nI, nJ, nK, nDim, j0_, nJp1_, k0_, nKp1_
  use BATL_geometry, ONLY: Dim2_, Dim3_, IsPeriodic_D, DomainSize_D

  implicit none

  SAVE

  private ! except

  public:: read_amr_geometry
  public:: init_amr_geometry
  public:: apply_amr_geometry
  public:: clean_amr_geometry
  public:: block_inside_region
  public:: region_value  ! return a 0 to 1 real number

  ! compatible with BATSRUS way of doing AMR
  logical, public :: IsBatsrusAmr = .true. 

  ! Masking criteria for use
  ! This belongs more in batl_amr_criteria but also
  ! needed here so it is declared here.
  logical, public, allocatable :: UseCrit_IB(:,:)
  logical, public, allocatable :: UseCrit_B(:)

  ! We have only two values to decide the criteria
  ! for geometric refinement: cell size (1) and AMR level (2)
  integer, parameter, public :: nGeoCrit = 2

  ! Number of geometrical criteria actually used
  integer, public :: nCritGeoUsed = 0

  ! Number of geometrical criteria from #GRIDLEVEL/RESOLUTION commands
  integer, public :: nCritGrid    = 0

  ! Number of geometric criteria from #AMRCRITERIALEVEL/RESOLUTION commands
  integer, public :: nCritDxLevel = 0

  ! we have read in new data from PARAM.in
  logical, public :: IsNewGeoParam =.false.

  integer, public, parameter :: MaxArea = 100, lNameArea = 20

  type, public :: AreaType 
     character(lNameArea) :: NameRegion
     character(lNameArea) :: NameShape
     real                 :: Resolution
     integer              :: Level
     real                 :: Center_D(nDim)
     real                 :: Size_D(nDim)
     real                 :: Radius1
     real                 :: Taper
     logical              :: DoRotate
     real, allocatable    :: Rotate_DD(:,:)
  end type AreaType

  ! index zero contains the 'all' area
  type(AreaType), target, public :: AreaGeo_I(0:MaxArea)

  ! Choosing which blocks we want to refine is based on a list of criteria 
  ! and a set of upper (refine)  and lower (coarsen) limits. The criteria 
  ! can be external or calculated internally by estimating the 
  ! numerical errors (calc_error_amr_criteria) based on the state variables.

  ! Storing All the AMR criteria
  integer, public:: nAmrCrit = 0
  integer, public:: nAmrCritUsed = 0
  real, public, allocatable:: AmrCrit_IB(:,:)

  real,    public, allocatable:: CoarsenCritAll_I(:), RefineCritAll_I(:)
  integer, public, allocatable:: iVarCritAll_I(:),iResolutionLimit_I(:)
  real,    public, allocatable:: ResolutionLimit_I(:)
  !type(AreaType), public, allocatable :: AreaAll_I(:)
  integer, public, allocatable:: iAreaIdx_II(:,:)

  ! Storing names of areas for each criteria given by #AMRCRITERIA.....
  integer, public, allocatable :: nAreaPerCritAll_I(:)

  ! Local variables

  real    :: InitialResolution = -1.0
  integer :: initial_refine_levels = 0
  real    :: AreaResolution  = 0.0
  integer :: nLevelArea = 0

  ! Name of area being processed
  character(len=lNameArea) :: NameShape

  ! index for the parallel and two perpendicular directions
  integer:: iPar
  integer, allocatable:: iPerp_I(:)

contains
  !============================================================================
  subroutine init_amr_geometry

    use BATL_grid,         ONLY: CoordMin_D,CoordMax_D
    use BATL_geometry,     ONLY: x_
    use BATL_tree,         ONLY: nRoot_D

    integer :: iGeo, iVar, iGeoAll
    logical, parameter :: DoTest=.false.
    character(len=*), parameter:: NameSub = 'init_amr_geometry'
    !-------------------------------------------------------------------------
    ! Fix resolutions (this depends on domain size set in BATL_grid)
    if(InitialResolution > 0.0) initial_refine_levels = nint( &
         alog(((CoordMax_D(x_) - CoordMin_D(x_)) / (nRoot_D(x_) * nI))  &
         / InitialResolution) / alog(2.0) )

    do iGeo = 1, nCritGeoUsed
       AreaResolution = AreaGeo_I(iGeo) % Resolution

       if(AreaResolution <= 0.0) then
          ! Negative value is understood as a "grid level"
          nLevelArea = ceiling(abs(AreaResolution))
          AreaGeo_I(iGeo) % Level = sign(nLevelArea, nint(AreaResolution))
          ! Set actual resolution
          AreaGeo_I(iGeo) % Resolution = (CoordMax_D(x_)-CoordMin_D(x_)) &
               / (nRoot_D(x_) * nI * 2.0**nLevelArea)
       else
          AreaGeo_I(iGeo) % Level = ceiling( &
               alog(((CoordMax_D(x_)-CoordMin_D(x_)) / (nRoot_D(x_) * nI))  &
               / AreaResolution) / alog(2.0) )
       end if

    end do

    iGeo = 0
    do iGeoAll = 1, nCritGeoUsed
       iVar = 1
       ! exclude named areas
       if( AreaGeo_I(iGeoAll) % NameRegion /= "NULL") CYCLE
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

       ! All geometrical criteia is a comparison to grid values
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

    ! Add geometric info for BATSRUS AMR type of params
    iGeo = 1
    do iGeoAll = 1, nCritGeoUsed
       if( AreaGeo_I(iGeoAll) % NameRegion /= "NULL") CYCLE
       iAreaIdx_II(1,nAmrCritUsed+iGeo) = iGeoAll
       nAreaPerCritAll_I(nAmrCritUsed+iGeo) = 1
       iGeo = iGeo + 1
    end do

    if(DoTest .and. iProc == 0) then
       do iGeoAll = 0, nCritGeoUsed
          write(*,*) NameSub,' iGeoAll=', iGeoAll
          write(*,*) "Region name      :: ", AreaGeo_I(iGeoAll)%NameRegion
          write(*,*) "Shape name       :: ", AreaGeo_I(iGeoAll)%NameShape
          write(*,*) "Region resolution:: ", AreaGeo_I(iGeoAll)%Resolution
          write(*,*) "Region level     :: ", AreaGeo_I(iGeoAll)%Level
          write(*,*) "Region Center_D  :: ", AreaGeo_I(iGeoAll)%Center_D
          write(*,*) "Region Size_D    :: ", AreaGeo_I(iGeoAll)%Size_D
          write(*,*) "Region Radius1   :: ", AreaGeo_I(iGeoAll)%Radius1
          write(*,*) "Area Taper       :: ", AreaGeo_I(iGeoAll)%Taper
          write(*,*) "Area DoRotate    :: ", AreaGeo_I(iGeoAll)%DoRotate
          if(AreaGeo_I(iGeoAll)%DoRotate) write(*,*) &
               "Area Rotate_DD :: ", AreaGeo_I(iGeoAll)%Rotate_DD
          write(*,*) "--------------------"
       end do
    end if

    IsNewGeoParam =.false.

  end subroutine init_amr_geometry

  !===========================================================================
  subroutine clean_amr_geometry

    if(allocated(UseCrit_IB))       deallocate(UseCrit_IB)
    if(allocated(UseCrit_B))        deallocate(UseCrit_B)
    if(allocated(AmrCrit_IB))       deallocate(AmrCrit_IB)
    if(allocated(CoarsenCritAll_I)) deallocate(CoarsenCritAll_I)
    if(allocated(RefineCritAll_I))  deallocate(RefineCritAll_I)
    if(allocated(ResolutionLimit_I))deallocate(ResolutionLimit_I)
    if(allocated(iVarCritAll_I))    deallocate(iVarCritAll_I)
    !if(allocated(AreaAll_I))       deallocate(AreaAll_I)

    ! nCritGeoUsed        = 0
    IsNewGeoParam         = .true.
    InitialResolution     = -1.0 
    initial_refine_levels = 0
    AreaResolution        = 0.0
    nAmrCrit              = 0
    nAmrCritUsed          = 0

  end subroutine clean_amr_geometry
  !============================================================================
  subroutine set_i_par_perp

    ! Set the indexes parallel and perpendicular to the symmetry axis
    ! of various shapes

    integer:: l
    character:: CharLast
    !------------------------------------------------------------------------
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
  subroutine apply_amr_geometry(iBlock, Area, UseBlock, DoCalcCritIn,&
       user_amr_geometry)

    !DESCRIPTION:
    ! Set UseBlock to .true. if block iBlock touches the region
    ! described by Area.
    ! if the area has a finer resolution than the block

    ! Set UseCrit_IB to true  for block touching the predefined areas
    ! if the area has a finer resolution than the block
    ! Also return max cell size for the block in CritGeo_IB

    ! WARNING
    ! user_specify_refinement,     not suported
    ! currentsheet,                not suported

    use BATL_geometry, ONLY: gen_to_radius,r_, Phi_, Theta_,  &
         IsCartesianGrid, IsLogRadius, IsGenRadius, IsRLonLat, &
         IsCylindrical
    use BATL_grid,     ONLY: Xyz_DNB
    use BATL_size,     ONLY: nINode, nJNode, nKNode, nDim
    use ModNumConst,   ONLY: cTiny, cRadToDeg
    use BATL_grid,     ONLY: CoordMin_DB, CoordMax_DB

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

    ! LOCAL VARIABLES:
    logical :: DoRefine
    integer, parameter :: nCorner = 2**nDim
    real, allocatable, save:: Corner_DI(:,:)
    real     :: DistMin_D(nDim), DistMax_D(nDim), Radius1Sqr
    integer  :: i, j, k, iDim, iCorner

    ! These variables are needed for generalized coordinates
    real :: Xyz_D(nDim)
    real, dimension(nINode,nJNode,nKNode):: x_N, y_N, z_N, R2_N

    character(len=*), parameter:: NameSub = 'apply_amr_geometry'
    logical, parameter:: DoTest = .false.

    logical :: DoCalcCrit
    !--------------------------------------------------------------------------
    !if(nCritGeoUsed <= 0) RETURN

    DoCalcCrit = .true.
    if(present(DoCalcCritIn)) DoCalcCrit =  DoCalcCritIn

    if(DoCalcCrit) call calc_crit(iBlock)

    UseBlock = .false.

    NameShape = Area % NameShape

    if(DoTest) write(*,*) NameSub,'NameShape, Resolution=',&
         ' ',trim(NameShape),' ',Area % Resolution

    ! Treat special cases first
    select case(NameShape)
    case('all', 'initial')
       UseBlock = .true.
       RETURN
    case('usergeo')
       if(present(user_amr_geometry)) then
          call user_amr_geometry(iBlock, -1, UseBlock)
          RETURN
       else
          call CON_stop(NameSub //' ERROR: need user_amr_geometry')
       end if
    end select

    if(.not.allocated(Corner_DI)) allocate(Corner_DI(nDim,nCorner))

    ! Blocks corner coordinates
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

    ! Check if it is a brick defined in generalized coordinates
    if(NameShape == "brick_gen" .and. .not.IsCartesianGrid)then

!!! TO BE IMPROVED !!!

       ! Convert angles to degrees and ln(r)/genr to r
       ! Make sure that phi=360 is not 0 but really 360
       if(IsLogRadius) Corner_DI(r_,:)     = exp(Corner_DI(r_,:))
       if(IsGenRadius) then
          do iCorner = 1, nCorner
             call gen_to_radius(Corner_DI(r_,iCorner))
          end do
       endif
       if(IsRLonLat)then
          Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
          Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg
       else if(IsCylindrical)then
          Corner_DI(Phi_,:) = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
          where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
          where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
       endif       
    end if

    if(DoTest) write(*,*) "post Corner_DI = ",Corner_DI

    ! Shift corner coordinates to the center of area
    do iCorner = 1, nCorner
       Corner_DI(:,iCorner) = &
            Corner_DI(:,iCorner) - Area % Center_D
    end do

    ! Rotate corners into the orientation of the area if required
    if(Area % DoRotate) &
         Corner_DI = matmul(Area % Rotate_DD, Corner_DI(1:nDim,:))

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

    if(DoTest)then
       write(*,*) NameSub,' DistMin_D=',DistMin_D
       write(*,*) NameSub,' DistMax_D=',DistMax_D
    end if

    ! Check if this area is intersecting with the block
    call set_i_par_perp
    select case(NameShape)
    case('brick', 'brick_gen')
       DoRefine = all(DistMin_D < 1)
    case('sphere')
       DoRefine = sum(DistMin_D**2) < 1
    case('shell')
       ! Check if block intersects with the enclosing sphere
       ! but it is not fully inside the inner sphere
       DoRefine = sum(DistMin_D**2) < 1 .and. sum(DistMax_D**2)>Radius1Sqr
    case('cylinderx','cylindery','cylinderz')
       DoRefine = DistMin_D(iPar) < 1 .and. sum(DistMin_D(iPerp_I)**2) < 1
    case('ringx', 'ringy', 'ringz')
       ! Check if block intersects with the enclosing cylinder
       ! but it is not fully inside the inner cylinder
       DoRefine = DistMin_D(iPar) < 1 .and. sum(DistMin_D(iPerp_I)**2) < 1 &
            .and. sum(DistMax_D(iPerp_I)**2) > Radius1Sqr 
    case default
       call CON_stop(NameSub //' ERROR: Unknown NameShape = '//NameShape)

    end select

    if(DoTest)write(*,*) NameSub,' DoRefine (from corners)=',DoRefine

    if(DoRefine) UseBlock = .true.
    
    if(NameShape == 'brick_gen' .or. IsCartesianGrid) RETURN

    ! Non-cartesian case
    UseBlock = .true.

    ! Check all nodes of the block
    do k=1, nKNode; do j=1, nJNode; do i=1, nINode

       ! Shift to area center
       Xyz_D = Xyz_DNB(1:nDim,i,j,k,iBlock) - Area % Center_D

       ! Rotate into area coordinates
       if(Area % DoRotate) Xyz_D = matmul(Area % Rotate_DD, Xyz_D)

       ! Rescale coordinates to the size of the area in all directions
       Xyz_D = Xyz_D / Area % Size_D

       ! We only need the absolute values of the coordinates
       x_N(i,j,k) = abs(Xyz_D(1))
       if(nDim>1) y_N(i,j,k) = abs(Xyz_D(Dim2_))
       if(nDim>2) z_N(i,j,k) = abs(Xyz_D(Dim3_))
    end do; end do; end do
    if(nDim < 2) y_N = 0.0
    if(nDim < 3) z_N = 0.0

    select case( NameShape)
    case('brick')
       if(any( x_N<1.0 .and. y_N<1.0 .and. z_N<1.0) ) RETURN
    case('sphere')
       if(any(x_N**2 + y_N**2 + z_N**2 < 1.0) ) RETURN
    case('shell')
       R2_N = x_N**2 + y_N**2 + z_N**2
       if(any(R2_N < 1.0 .and. R2_N > Radius1Sqr )) RETURN
    case('cylinderx')
       if(any(x_N < 1.0 .and. y_N**2+z_N**2 < 1.0 ) ) RETURN
    case('cylindery')
       if(any(y_N < 1.0 .and. x_N**2+z_N**2 < 1.0 ) ) RETURN
    case('cylinderz')
       if(any(z_N < 1.0 .and. x_N**2+y_N**2 < 1.0 ) ) RETURN
    case('ringx')
       R2_N = y_N**2+z_N**2
       if(any(x_N < 1.0 .and. R2_N < 1.0 .and. R2_N > Radius1Sqr)) RETURN
    case('ringy')
       R2_N = x_N**2+z_N**2
       if(any(y_N < 1.0 .and. R2_N < 1.0 .and. R2_N > Radius1Sqr)) RETURN
    case('ringz')
       R2_N = x_N**2+y_N**2
       if(any(z_N < 1.0 .and. R2_N < 1.0 .and. R2_N > Radius1Sqr)) RETURN
    case default
       call CON_stop(NameSub //' ERROR: Unknown NameShape = '//NameShape)
    end select

    ! The block is not inside
    UseBlock = .false.

    if(DoTest) write(*,*)NameSub,' DoRefine final=',UseBlock


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

    use ModReadParam,      ONLY: read_var, lStringLine
    use ModNumConst,       ONLY: cDegToRad
    use ModCoordTransform, ONLY: &
         rot_matrix, rot_matrix_x, rot_matrix_y, rot_matrix_z

    character(len=*),  intent(inout) :: NameCommand
    logical, optional, intent(in)    :: UseStrictIn
    integer, optional, intent(inout) :: InitLevelInOut 
    real,    optional, intent(inout) :: InitResInOut 

    character(len=lStringLine):: StringShape = 'all'
    real    :: RadiusArea    =0.0
    logical :: DoReadAreaCenter = .false.
    real    :: XyzStartArea_D(nDim)=0.0, XyzEndArea_D(nDim)=0.0
    real    :: xRotateArea   = 0.0
    real    :: yRotateArea   = 0.0
    real    :: zRotateArea   = 0.0
    logical :: DoTaperArea   = .false.
    logical :: DoStretch     = .false.

    character(len=lNameArea):: NameRegion
    integer :: i, iDim
    logical :: UseStrict

    type(AreaType), pointer:: Area

    character (len=17) :: NameSub='read_amr_geometry'
    !-------------------------------------------------------------------------

    IsNewGeoParam =.true.

    UseStrict = .true.
    if(present(UseStrictIn)) UseStrict = UseStrictIn

    NameRegion ="NULL"
    if(index(NameCommand, "REGION") > 0) then
       ! Read name of region for #REGION (or #AMRREGION) command
       call read_var('NameRegion',NameRegion, IsLowerCase=.true.)
       AreaResolution = 0
    else
       ! Read grid level/resolution for #GRIDLEVEL/RESOLUTION command
       if(index(NameCommand,"RESOLUTION") > 0)then
          call read_var('AreaResolution', AreaResolution) 
       else
          call read_var('nLevelArea', nLevelArea)
          ! Store level as a negative integer resolution.
          AreaResolution = -nLevelArea
       end if
    end if
    call read_var('TypeRegion', StringShape, IsLowerCase=.true.)
    NameShape = adjustl(StringShape)

    if(StringShape(1:4) == 'init')then
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
    Area => AreaGeo_I(nCritGeoUsed)

    if(NameRegion == "NULL") nCritGrid = nCritGrid + 1

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

    Area % Resolution = AreaResolution
    Area % NameRegion = NameRegion
    ! Set the default center to be the origin, 
    ! size 1, smaller radius 0, no tapering, no rotation
    Area%Center_D  = 0.0
    Area%Size_D    = 1.0
    Area%Radius1   = 0.0
    Area%Taper     = 0.0
    Area%DoRotate  = .false.

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
    if(i>0) StringShape = StringShape(:i-1)//StringShape(i+1:)

    DoReadAreaCenter = (i < 1 .and. StringShape(1:3) /= 'box')

    ! Store name after all the options were removed from StringShape
    NameShape = StringShape

    if(NameShape == 'user') NameShape = 'usergeo'

    ! Store NameShape
    Area % NameShape = NameShape

    ! The 'all' types does not need any parameters
    if(NameShape == 'all') RETURN

    ! Read center of area if needed
    if(DoReadAreaCenter)then
       do iDim = 1, nDim
          call read_var("Center", Area%Center_D(iDim))
       end do
    endif

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
       if(NameShape(1:4) == "cone")then
          ! A cone is a funnel with 0 radius on one end
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

       nCritGeoUsed = nCritGeoUsed - 1
       RETURN
    end select

    ! Read the tapering width if required
    if(DoTaperArea) call read_var('Taper', Area % Taper)

    if(Area % DoRotate)then
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

  end subroutine read_amr_geometry

  !============================================================================
  subroutine block_inside_region(iRegion_I, iBlock, nValue, StringLocation, &
       IsInside, IsInside_I, Value_I)

    use BATL_size, ONLY: nINode, nJNode, nKNode, nIJK, MaxIJK_D, MinIJK_D
    use BATL_grid, ONLY: Xyz_DGB, Xyz_DNB

    ! Check the intersection of block iBlock with one or more regions
    ! indexed by the iRegion_I array. Positive region index means 
    ! an "OR", while negative region index means "AND NOT".
    ! If the first region index is negative, then it is relative to
    ! an "all" condition.
    !
    ! The IsInside is present it is set to true if the block intersects 
    ! any of the + regions and avoids all the - regions.
    ! 
    ! If IsInside_I is present, it is set true for each point of the block
    ! (defined by StringLocation) that is inside the region(s).
    !
    ! If Value_I is present, then it is set to 1 inside the region(s), 
    ! 0 outside the regions, and between 0 and 1 inside the tapering
    ! for each point of the block defined by StringLocation.
    ! 
    ! The StringLocation argument tells which points of the block are to be 
    ! evaluated for IsInside_I and/or Value_I. Possible values are 
    ! "cell", "ghost", "xface", "yface", "zface", "node"
    ! Only the first character is significant and it is case insensitive.

    integer, intent(in):: iRegion_I(:)
    integer, intent(in):: iBlock
    integer, intent(in):: nValue

    character(*), optional, intent(in) :: StringLocation
    logical,      optional, intent(out):: IsInside
    logical,      optional, intent(out):: IsInside_I(nValue)
    real,         optional, intent(out):: Value_I(nValue)

    integer:: iRegion, nRegion, iArea
    integer:: nPoint
    character:: NameLocation

    logical:: DoBlock, DoMask, DoValue

    logical             :: IsInsideOld
    logical, allocatable:: IsInsideOld_I(:)
    real,    allocatable:: ValueOld_I(:)

    real, allocatable, save:: Xyz_DI(:,:)

    logical, parameter:: DoTest=.false.

    character(len=*), parameter:: NameSub='region_block'
    !------------------------------------------------------------------------
    DoBlock  = present(IsInside)
    DoMask   = present(IsInside_I)
    DoValue  = present(Value_I)

    if(.not.(DoBlock .or. DoMask .or. DoValue)) call CON_stop(NameSub// &
         ': no output argument is present')

    nRegion = size(iRegion_I)
    if(nRegion < 1) call CON_stop(NameSub//': empty region index array')

    nPoint = nValue
    if(DoMask  .and. nRegion > 1) allocate(IsInsideOld_I(nPoint))
    if(DoValue .and. nRegion > 1) allocate(ValueOld_I(nPoint))

    ! Default location is nodes
    NameLocation = "n"
    if(present(StringLocation)) NameLocation = StringLocation(1:1)

    ! Check nodes of block by default (optimize for corners???)
    if(nPoint == 0) nPoint = nINode*nJNode*nKNode

    ! Allocate coordinate array
    allocate(Xyz_DI(nDim,nPoint))

    select case(NameLocation)
    case('c', 'C')
       if(nPoint /= nIJK) call CON_stop(NameSub// &
            ': incorrect number of points for cell centers')

       Xyz_DI = reshape(Xyz_DGB(1:nDim,1:nI,1:nJ,1:nK,iBlock), &
            (/nDim, nPoint/))

    case('g', 'G')
       if(nPoint /= product(MaxIJK_D-MinIJK_D+1)) call CON_stop(NameSub// &
            ': incorrect number of points for cell centers including ghosts')

       Xyz_DI = reshape(Xyz_DGB(:,:,:,:,iBlock), (/nDim, nPoint/))

    case('x', 'X')
       if(nPoint /= (nI+1)*nJ*nK) call CON_stop(NameSub// &
            ': incorrect number of points for X faces')

       Xyz_DI = reshape( &
            0.5*(Xyz_DGB(1:nDim,0:nI  ,1:nJ,1:nK,iBlock) &
            +    Xyz_DGB(1:nDim,1:nI+1,1:nJ,1:nK,iBlock)), (/nDim, nPoint/))

    case('y', 'Y')
       if(nPoint /= nI*(nJ+1)*nK) call CON_stop(NameSub// &
            ': incorrect number of points for X faces')

       Xyz_DI = reshape( &
            0.5*(Xyz_DGB(1:nDim,1:nI,j0_:nJ, 1:nK,iBlock) &
            +    Xyz_DGB(1:nDim,1:nI,1:nJp1_,1:nK,iBlock)), (/nDim, nPoint/))

    case('z', 'Z')
       if(nPoint /= nI*nJ*(nK+1)) call CON_stop(NameSub// &
            ': incorrect number of points for X faces')

       Xyz_DI = reshape( &
            0.5*(Xyz_DGB(1:nDim,1:nI,1:nJ,k0_:nK ,iBlock) &
            +    Xyz_DGB(1:nDim,1:nI,1:nJ,1:nKp1_,iBlock)), (/nDim, nPoint/))

    case('n', 'N')
       if(nPoint /= nINode*nJNode*nKNode) call CON_stop(NameSub// &
            ': incorrect number of points for nodes')

       Xyz_DI = reshape(Xyz_DNB(1:nDim,:,:,:,iBlock), (/nDim, nPoint/))

    case default
       write(*,*) NameSub,': StringLocation=', StringLocation
       call CON_stop(NameSub//' incorrect value for StringLocation')
    end select

    if(DoTest)write(*,*) NameSub, 'iBlock, nPoint, NameLoc, size(Xyz_DI) = ', &
         iBlock, nPoint, NameLocation, size(Xyz_DI)

    ! Evaluate all regions
    do iRegion = 1, nRegion

       ! Store results obtained from previous regions
       if(iRegion > 0)then
          if(DoBlock) IsInsideOld   = IsInside
          if(DoMask)  IsInsideOld_I = IsInside_I
          if(DoValue) ValueOld_I    = Value_I
       end if

       ! Get area index and evaluate it
       iArea = iRegion_I(iRegion)
       call region_value(nPoint, Xyz_DI, AreaGeo_I(abs(iArea)), &
            IsInside, IsInside_I, Value_I)

       if(DoTest)write(*,*)'maxval(Value_I)=', maxval(Value_I)

       if(iRegion == 1)then
          if(iArea < 0)then
             if(DoBlock) IsInside   = .not. IsInside
             if(DoMask)  IsInside_I = .not. IsInside_I
             if(DoValue) Value_I    = 1 - Value_I
          end if
       else
          ! Combine last region info with previous
          if(iArea < 0)then
             if(DoBlock) IsInside   = IsInsideOld   .and. .not. IsInside
             if(DoMask)  IsInside_I = IsInsideOld_I .and. .not. IsInside_I
             if(DoValue) Value_I    = min(ValueOld_I, 1 - Value_I)
          else
             if(DoBlock) IsInside   = IsInsideOld   .or. IsInside
             if(DoMask)  IsInside_I = IsInsideOld_I .or. IsInside_I
             if(DoValue) Value_I    = max(ValueOld_I, Value_I)
          end if
       end if
    end do

    deallocate(Xyz_DI)
    if(allocated(IsInsideOld_I)) deallocate(IsInsideOld_I)
    if(allocated(ValueOld_I))    deallocate(ValueOld_I)

  end subroutine block_inside_region
  !============================================================================

  subroutine region_value(nPoint, Xyz_DI, Area, IsInside, IsInside_I, Value_I)

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

    real, allocatable, save:: Norm_DI(:,:)

    integer  :: iPoint, iDim

    logical:: DoBlock, DoMask, DoValue, DoBlockOnly
    
    logical:: DoTaper
    real:: Xyz_D(nDim), Size_D(nDim)
    real:: Radius, RadiusSqr, Dist1, Dist2
    real:: Radius1, Radius1Sqr, Slope1
    real:: Taper, TaperFactor_D(nDim), TaperFactor1_D(nDim)
    real:: Slope_D(nDim)
    real, allocatable, save:: SlopePerp_D(:)

    logical, parameter:: DoTest = .false.
    character(len=*), parameter:: NameSub = 'region_value'
    !--------------------------------------------------------------------------
    DoBlock  = present(IsInside)
    DoMask   = present(IsInside_I)
    DoValue  = present(Value_I)
    DoBlockOnly = DoBlock .and. .not. (DoMask .or. DoValue)

    if(DoTest) write(*,*) NameSub, &
         ': nPoint, DoBlock, DoMask, DoValue, DoBlockOnly=',&
         nPoint, DoBlock, DoMask, DoValue, DoBlockOnly

    NameShape = Area % NameShape

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

    ! Check if it is a brick defined in generalized coordinates
!    if(NameShape == "brick_gen" .and. .not.IsCartesianGrid)then
!

    if(allocated(Norm_DI))then
       if(size(Norm_DI) < nDim*nPoint) deallocate(Norm_DI)
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
       ! Modified tapering slope for the sides of funnel, cone and doublecone
       if(.not.allocated(SlopePerp_D)) allocate(SlopePerp_D(nDim-1))
       if(iPar > 0) SlopePerp_D = Slope_D(iPerp_I) &
            /sqrt(1 + (Slope1*Size_D(iPerp_I)/Size_D(iPar))**2)
    endif

    ! Transform point coordinates for easy comparison
    do iPoint = 1, nPoint

       ! Shift to area center
       Xyz_D = Xyz_DI(:,iPoint) - Area % Center_D

       do iDim = 1, nDim
          if(IsPeriodic_D(iDim))then
             if(Xyz_D(iDim) > 0.5*DomainSize_D(iDim)) &
                  Xyz_D(iDim) = Xyz_D(iDim) - DomainSize_D(iDim)
             if(Xyz_D(iDim) < -0.5*DomainSize_D(iDim)) &
                  Xyz_D(iDim) = Xyz_D(iDim) + DomainSize_D(iDim)
          end if
       end do

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
       !case('usergeo')
       ! call user_region(nPoint, Norm_DI, Area, Value_I)
    case('brick')
       Norm_DI = abs(Norm_DI)
       if(.not.DoTaper)then
          do iPoint = 1, nPoint
             if(maxval(Norm_DI(:,iPoint)) > 1) CYCLE
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
             if(sum(Norm_DI(:,iPoint)**2) > 1) CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Dist1 = sum(Norm_DI(:,iPoint)**2)
             if(Dist1 <= 1)then
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
             if(Dist1 > 1 .or. Dist1 < Radius1Sqr)CYCLE
             if(DoBlock) IsInside = .true.
             if(DoBlockOnly) EXIT
             if(DoMask)  IsInside_I(iPoint) = .true.
             if(DoValue) Value_I(iPoint) = 1
          end do
       else
          do iPoint = 1, nPoint
             Dist1 = sum(Norm_DI(:,iPoint)**2)
             if(Dist1 > 1)then
                Dist2 = sum((TaperFactor_D*Norm_DI(:,iPoint))**2)
                if(Dist2 >= 1) CYCLE     ! outside outer radius
                ! Use a roughly linear function between the ellipsoids
                Dist1 = sqrt(Dist1) - 1
                Dist2 = 1 - sqrt(Dist2)
                Value_I(iPoint) = Dist2/(Dist1 + Dist2)
             elseif(Dist1 < Radius1Sqr)then
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
             if(abs(Norm_DI(iPar,iPoint)) > 1) CYCLE
             if(sum(Norm_DI(iPerp_I,iPoint)**2) > 1) CYCLE
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
             if(Dist1 > 1 .or. Dist1 < Radius1Sqr)CYCLE
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
             if(abs(Norm_DI(iPar,iPoint) - 0.5) > 0.5) CYCLE
             Radius = Radius1 + Slope1*Norm_DI(iPar,iPoint)
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 > Radius**2) CYCLE
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
             if(abs(Norm_DI(iPar,iPoint)-0.5) > 0.5) CYCLE
             if(sum(Norm_DI(iPerp_I,iPoint)**2) &
                  > Norm_DI(iPar,iPoint)) CYCLE
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
             if(abs(Norm_DI(iPar,iPoint)) > 1) CYCLE
             if(sum(Norm_DI(iPerp_I,iPoint)**2) &
                  > Norm_DI(iPar,iPoint)**2) CYCLE
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
             ! Radius of double cone at this cross section
             Radius = abs(Norm_DI(iPar,iPoint))
             Dist1 = sum(Norm_DI(iPerp_I,iPoint)**2)
             if(Dist1 < Radius**2) CYCLE             ! inside perpendicular
             ! Use taper factor specific for this particular cone radius
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
       call CON_stop(NameSub//' ERROR: unknown NameShape='//trim(NameShape))
    end select

    ! Get logical information if not yet calculated
    if(DoTaper)then
       if(DoMask)  IsInside_I =     Value_I == 1
       if(DoBlock) IsInside   = any(Value_I == 1)
    end if

    if(DoTest .and. DoValue) &
         write(*,*) NameSub, ': maxval(Value_I)=', maxval(Value_I)

  end subroutine region_value

end module BATL_amr_geometry
