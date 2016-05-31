module ModGrid

  use ModSize, ONLY: &
       nDim, nVar, nLat, nLon, nNode, &
       RMin, iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_, R_, Lat_, Lon_

  implicit none

  SAVE

  private ! except

  public:: set_grid_param, init_grid, get_node_indexes
  public:: iComm, iProc, nProc, nBlock, Proc_, Block_, Begin_, End_
  public:: LatMin, LatMax, LonMin, LonMax
  public:: iGrid_IA, iNode_II, iNode_B, State_VIB, CoordOrigin_DA

  !\
  ! MPI information
  !----------------------------------------------------------------------------
  integer:: iComm = -1
  integer:: iProc = -1
  integer:: nProc = -1
  !/
  !\
  ! Grid info
  ! Containers for coordinates and data
  !----------------------------------------------------------------------------
  ! Size of angular grid, latitude and longitude, at origin surface R=RMin
  real:: LatMin, LatMax, DLat
  real:: LonMin, LonMax, DLon
  ! Mark that grid has been set
  logical:: IsSetGrid = .false.
  !----------------------------------------------------------------------------
  ! Said angular grids itself; each field line is identified by latitude
  ! and longitude of the origin point at surface R=RMin as it is set
  ! at the beginning of simulation;
  ! 1st index - three HGI coordinates (R is added for compliteness)
  ! 2nd index - node number (equivalent to line number)
  real,    allocatable:: CoordOrigin_DA(:,:)
  !----------------------------------------------------------------------------
  ! Node number based on the field line identified by 2 angular grid indices,
  ! latitude and longitude;
  ! 1st index - latitude index
  ! 2nd index - longitude index
  integer, allocatable:: iNode_II(:,:)
  !----------------------------------------------------------------------------
  ! Number of blocks on this processor
  integer:: nBlock
  !----------------------------------------------------------------------------
  ! Node number based on the local block number
  ! 1st index - block number
  integer, allocatable:: iNode_B(:)
  !----------------------------------------------------------------------------
  ! Various house-keeping information about the node/line;
  ! 1st index - identification of info field
  ! 2nd index - node number
  integer, allocatable:: iGrid_IA(:,:)
  !----------------------------------------------------------------------------
  ! Number of info fields per node and their identifications
  integer, parameter:: nNodeIndexes = 4
  integer, parameter:: &
       Proc_  = 1, & ! Processor that has this line/node
       Block_ = 2, & ! Block that has this line/node
       Begin_ = 3, & ! Index of the 1st particle on this line/node
       End_   = 4    ! Index of the last particle on this line/node
  !----------------------------------------------------------------------------
  ! State vector;
  ! 1st index - identification of variable
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real,    allocatable:: State_VIB(:,:,:)
  !/

contains
  
  subroutine set_grid_param
    use ModReadParam, ONLY: read_var
    character(len=*), parameter:: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    call read_var('LatMin', LatMin)
    call read_var('LatMax', LatMax)
    if(LatMax <= LatMin)&
         call CON_stop('Origin surface grid is inconsistent:'//NameSub)
    DLat = (LatMax - LatMin) / nLat

    call read_var('LonMin', LonMin)
    call read_var('LonMax', LonMax)
    if(LonMax <= LonMin)&
         call CON_stop('Origin surface grid is inconsistent:'//NameSub)
    DLon = (LonMax - LonMin) / nLon
    IsSetGrid = .true.
  end subroutine set_grid_param

  !============================================================================

  subroutine init_grid
    ! allocate the grid used in this model
    use ModUtilities, ONLY: check_allocate
    integer:: iError
    integer:: iLat, iLon, iNode, iBlock, iProcNode
    character(LEN=*),parameter:: NameSub='SP:init_grid'
    !--------------------------------------------------------------------------
    !\
    ! Check if everything's ready for initialization
    !/
    if(.not.IsSetGrid)&
         call CON_stop(NameSub//': grid is not set in PARAM.in file')

    !\
    ! distribute nodes between processors
    !/
    if(nNode < nProc)&
         call CON_stop('There are less processor than field lines:'//NameSub)
    nBlock = ((iProc+1)*nNode) / nProc - (iProc*nNode) / nProc
    !\
    ! check consistency
    !/
    if(nLat <= 0 .or. nLon <= 0)&
         call CON_stop('Origin surface grid is invalid:'//NameSub)
    if(iParticleMin > 0 .or. iParticleMax < 0)&
         call CON_stop('Origin surface is not included:'//NameSub)
    !\
    ! allocate data and grid containers
    !/
    allocate(iNode_II(nLat, nLon), stat=iError)
    call check_allocate(iError, NameSub//'iNode_II')
    allocate(iNode_B(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iNode_B')
    allocate(iGrid_IA(nNodeIndexes, nNode), stat=iError)
    call check_allocate(iError, NameSub//'iGrid_IA')
    allocate(CoordOrigin_DA(nDim, nNode), stat=iError)
    call check_allocate(iError, NameSub//'CoordOrigin_DA')
    allocate(State_VIB(nVar,iParticleMin:iParticleMax,nBlock), stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    !\
    ! fill grid containers
    !/
    iBlock = 1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iLon + nLon * (iLat-1)
          iNode_II(iLat, iLon) = iNode
          iProcNode = ceiling(real(iNode*nProc)/nNode) - 1
          if(iProcNode==iProc)then
             iNode_B(iBlock) = iNode
          end if
          iGrid_IA(Proc_, iNode) = iProcNode
          iGrid_IA(Block_,iNode) = iBlock
          iGrid_IA(Begin_,iNode) = 0
          iGrid_IA(End_,  iNode) = 0
          if(iNode == ((iProcNode+1)*nNode)/nProc)then
             iBlock = 1
          else
             iBlock = iBlock + 1
          end if
       end do
    end do
    !\
    ! fill data containers
    !/
    State_VIB = -1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iNode_II(iLat, iLon)
          CoordOrigin_DA(:, iNode) = &
               (/RMin, LatMin + (iLat-0.5)*DLat, LonMin + (iLon-0.5)*DLon/)
          iBlock = iGrid_IA(Block_, iNode)
          State_VIB(1:nDim,0,iBlock) = CoordOrigin_DA(:,iNode)
       end do
    end do
  end subroutine init_grid

  !============================================================================

  subroutine get_node_indexes(iNodeIn, iLonOut, iLatOut)
    ! return angular grid's indexes corresponding to this node
    integer, intent(in) :: iNodeIn
    integer, intent(out):: iLonOut
    integer, intent(out):: iLatOut
    !---------------------------------------------------------
    iLatOut = 1 + (iNodeIn-1) / nLon
    iLonOut = iNodeIn - nLon * (iLatOut-1)
  end subroutine get_node_indexes

  !============================================================================

  subroutine get_cell(CoordIn_D, iCellOut_D)
    real,    intent(in) :: CoordIn_D(nDim)
    integer, intent(out):: iCellOut_D(nDim)
    !--------------------------------------------------------------------------
    iCellOut_D(Particle_)  = nint( CoordIn_D(Particle_))
    iCellOut_D(OriginLat_) = nint((CoordIn_D(OriginLat_)-LatMin)/DLat + 0.5)
    iCellOut_D(OriginLon_) = nint((CoordIn_D(OriginLon_)-LonMin)/DLon + 0.5)
  end subroutine get_cell

  !============================================================================
  
  subroutine get_node(CoordIn_D, iNodeOut)
    real,    intent(in) :: CoordIn_D(nDim)
    integer, intent(out):: iNodeOut
    ! angular grid indices
    integer:: iLat, iLon
    !--------------------------------------------------------------------------
    iLat = nint((CoordIn_D(OriginLat_)-LatMin)/DLat + 0.5)
    iLon = nint((CoordIn_D(OriginLon_)-LonMin)/DLon + 0.5)
    iNodeOut = iNode_II(iLat, iLon)
  end subroutine get_node

  !============================================================================

  subroutine convert_to_hgi(CoordIn_D, CoordOut_D)
    real, intent(in) :: CoordIn_D(nDim)
    real, intent(out):: CoordOut_D(nDim)

    integer:: iBlock, iCell_D(nDim)
    !--------------------------------------------------------------------------
    call get_cell(CoordIn_D, iCell_D)
    iBlock = &
         iGrid_IA(Block_, iNode_II(iCell_D(OriginLat_), iCell_D(OriginLon_)))
    CoordOut_D((/R_, Lat_, Lon_/)) = &
         State_VIB((/R_,Lat_,Lon_/), iCell_D(Particle_), iBlock)
  end subroutine convert_to_hgi

  !============================================================================


end module ModGrid
