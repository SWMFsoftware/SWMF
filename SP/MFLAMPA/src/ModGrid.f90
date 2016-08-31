module ModGrid

  use ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nMomentumBin, nPitchAngleBin, &
       ROrigin, iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_

  implicit none

  SAVE

  private ! except

  public:: set_grid_param, init_grid, get_node_indexes, distance_to_next
  public:: fix_grid_consistency
  public:: iComm, iProc, nProc, nBlock, Proc_, Block_
  public:: LatMin, LatMax, LonMin, LonMax
  public:: iGridGlobal_IA, iGridLocal_IB, iNode_II, iNode_B
  public:: State_VIB, Distribution_IIIB, MomentumScale_I
  public:: Begin_, End_, Shock_, ShockOld_
  public:: nVar, R_, Lon_, Lat_, D_
  public:: Rho_, T_, Ux_,Uy_,Uz_,U_, Bx_,By_,Bz_,B_, RhoOld_, BOld_  
  public:: NameVar_V

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
  ! Size of angular grid, latitude and longitude, at origin surface R=ROrigin
  real:: LatMin, LatMax, DLat
  real:: LonMin, LonMax, DLon
  ! Mark that grid has been set
  logical:: IsSetGrid = .false.
  !----------------------------------------------------------------------------
  ! Said angular grids itself; each field line is identified by latitude
  ! and longitude of the origin point at surface R=ROrigin as it is set
  ! at the beginning of simulation;
  ! 1st index - three spherical coordinates (R is added for completeness)
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
  ! 2nd index - node number / block number
  integer, allocatable:: iGridGlobal_IA(:,:)
  integer, allocatable:: iGridLocal_IB(:,:)
  !----------------------------------------------------------------------------
  ! Number of info fields per node/block and their identifications
  integer, parameter:: nNodeIndexes = 2
  integer, parameter:: &
       Proc_  = 1, & ! Processor that has this line/node
       Block_ = 2    ! Block that has this line/node
  integer, parameter:: nBlockIndexes = 4
  integer, parameter:: &
       Begin_   = 1, & ! Index of the 1st particle on this line/node
       End_     = 2, & ! Index of the last particle on this line/node
       Shock_   = 3, & ! Current location of a shock wave
       ShockOld_= 4    ! Old location of a shock wave
  !----------------------------------------------------------------------------
  ! State vector;
  ! 1st index - identification of variable
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, allocatable:: State_VIB(:,:,:)
  !----------------------------------------------------------------------------
  ! Number of variables in the state vector and their identifications
  integer, parameter:: nVar = 16
  integer, parameter:: &
       R_     = 1, & ! Radial coordinate
       Lon_   = 2, & ! Longitude
       Lat_   = 3, & ! Latitude
       D_     = 5, & ! Distance to the next particle
       ! Current values
       Rho_   = 5, & ! Background plasma density
       T_     = 6, & ! Background temperature
       Ux_    = 7, &
       Uy_    = 8, &
       Uz_    = 9, &
       U_     =10, &
       Bx_    =11, & !
       By_    =12, & ! Background magnetic field
       Bz_    =13, & !
       B_     =14, & ! Magnitude of magnetic field
       ! Old values
       RhoOld_=15, & ! Background plasma density
       BOld_  =16 ! Magnitude of magnetic field
  ! variable names
  character(len=10), parameter:: NameVar_V(nVar) = (/&
       'R     ', &
       'Lon   ', &
       'Lat   ', &
       'D     ', &
       'Rho   ', &
       'T     ', &
       'Ux    ', &
       'Uy    ', &
       'Uz    ', &
       'U     ', &
       'Bx    ', &
       'By    ', &
       'Bz    ', &
       'B     ', &
       'RhoOld', &
       'BOld  '/)
  !----------------------------------------------------------------------------
  ! Distribution vector;
  ! Number of bins in the distribution is set in ModSize
  ! 1st index - log(momentum) bin
  ! 2nd index - cos(pitch-angle) bin
  ! 3rd index - particle index along the field line
  ! 4th index - local block number
  real, allocatable:: Distribution_IIIB(:,:,:,:)
  ! scale with respect to log(Momentum)
  real:: MomentumScale_I(nMomentumBin)
  !/

contains
  
  subroutine set_grid_param
    use ModReadParam, ONLY: read_var
    use ModNumConst, ONLY: cPi
    character(len=*), parameter:: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    call read_var('LonMin', LonMin)
    call read_var('LonMax', LonMax)
    if(LonMax <= LonMin)&
         call CON_stop('Origin surface grid is inconsistent:'//NameSub)
    ! convert angels from degrees to radians
    LonMax = LonMax * cPi / 180
    LonMin = LonMin * cPi / 180
    ! angular grid's step
    DLon = (LonMax - LonMin) / nLon

    call read_var('LatMin', LatMin)
    call read_var('LatMax', LatMax)
    if(LatMax <= LatMin)&
         call CON_stop('Origin surface grid is inconsistent:'//NameSub)
    ! convert angels from degrees to radians
    LatMax = LatMax * cPi / 180
    LatMin = LatMin * cPi / 180
    ! angular grid's step
    DLat = (LatMax - LatMin) / nLat

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
         call CON_stop('There are more processors than field lines:'//NameSub)
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
    allocate(iNode_II(nLon, nLat), stat=iError)
    call check_allocate(iError, NameSub//'iNode_II')
    allocate(iNode_B(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iNode_B')
    allocate(iGridGlobal_IA(nNodeIndexes, nNode), stat=iError)
    call check_allocate(iError, NameSub//'iGridGlobal_IA')
    allocate(iGridLocal_IB(nBlockIndexes, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iGridLocal_IB')
    allocate(CoordOrigin_DA(nDim, nNode), stat=iError)
    call check_allocate(iError, NameSub//'CoordOrigin_DA')
    allocate(State_VIB(nVar,iParticleMin:iParticleMax,nBlock), stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    allocate(Distribution_IIIB(&
         nMomentumBin,nPitchAngleBin,iParticleMin:iParticleMax,nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'Distribution_IIIB')
    !\
    ! fill grid containers
    !/
    iBlock = 1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iLon + nLon * (iLat-1)
          iNode_II(iLon, iLat) = iNode
          iProcNode = ceiling(real(iNode*nProc)/nNode) - 1
          if(iProcNode==iProc)then
             iNode_B(iBlock) = iNode
          end if
          iGridGlobal_IA(Proc_,   iNode)  = iProcNode
          iGridGlobal_IA(Block_,  iNode)  = iBlock
          iGridLocal_IB(Begin_,   iBlock) = 0
          iGridLocal_IB(End_,     iBlock) = 0
          iGridLocal_IB(Shock_,   iBlock) = iParticleMin - 1
          iGridLocal_IB(ShockOld_,iBlock) = iParticleMin - 1
          if(iNode == ((iProcNode+1)*nNode)/nProc)then
             iBlock = 1
          else
             iBlock = iBlock + 1
          end if
       end do
    end do
    !\
    ! reset and fill data containers
    !/
    Distribution_IIIB = tiny(1.0)
    State_VIB = -1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iNode_II(iLon, iLat)
          CoordOrigin_DA(:, iNode) = &
               (/ROrigin, LonMin + (iLon-0.5)*DLon, LatMin + (iLat-0.5)*DLat/)
          iBlock = iGridGlobal_IA(Block_, iNode)
          if(iProc == iGridGlobal_IA(Proc_, iNode))&
               State_VIB(1:nDim,0,iBlock) = CoordOrigin_DA(:,iNode)
       end do
    end do
  end subroutine init_grid

  !============================================================================

  subroutine fix_grid_consistency
    ! recompute some values (magnitudes of plasma velocity and magnetic field)
    ! so they are consistent with components for all lines
    integer:: iBlock, iParticle, iBegin, iEnd
    !--------------------------------------------------------------------------
    do iBlock = 1, nBlock
       iBegin = iGridLocal_IB(Begin_,iBlock)
       iEnd   = iGridLocal_IB(End_,  iBlock)
       do iParticle = iBegin, iEnd
          ! plasma speed
          State_VIB(U_,iParticle, iBlock) = &
               sqrt(sum(State_VIB(Ux_:Uz_,iParticle,iBlock)**2))
          ! magnetic field
          State_VIB(B_,iParticle, iBlock) = &
               sqrt(sum(State_VIB(Bx_:Bz_,iParticle,iBlock)**2))

          ! distances between particles
          if(iParticle < iGridLocal_IB(End_,  iBlock))&
               State_VIB(D_, iParticle, iBlock) = &
               distance_to_next(iParticle, iBlock)
       end do
       ! location of shock
       if(iGridLocal_IB(ShockOld_, iBlock) < iParticleMin)&
            iGridLocal_IB(ShockOld_, iBlock)= iBegin
       if(iGridLocal_IB(Shock_, iBlock) < iParticleMin)&
            iGridLocal_IB(Shock_, iBlock)   = iBegin
    end do
  end subroutine fix_grid_consistency

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
  function distance_to_next(iParticle, iBlock) result(Distance)
    ! the function returns distance to the next particle measured in Rs;
    ! formula for distance between 2 points in rlonlat system:
    !  Distance = R1**2 + R2**2 - 
    !     2*R1*R2*(Cos(Lat1)*Cos(Lat2) * Cos(Lon1-Lon2) + Sin(Lat1)*Sin(Lat2))
    ! NOTE: function doesn't check whether iParticle is last on the field line
    integer, intent(in):: iParticle
    integer, intent(in):: iBlock
    real               :: Distance

    real:: CosLat1xCosLat2, SinLat1xSinLat2, CosLon1MinusLon2
    real:: CosLat1PlusLat2, CosLat1MinusLat2
    !--------------------------------------------------------------------
    CosLat1PlusLat2 = cos(&
         State_VIB(Lat_,iParticle,iBlock) + State_VIB(Lat_,iParticle+1,iBlock))
    CosLat1MinusLat2 = cos(&
         State_VIB(Lat_,iParticle,iBlock) - State_VIB(Lat_,iParticle+1,iBlock))
    CosLon1MinusLon2 = cos(&
         State_VIB(Lon_,iParticle,iBlock) - State_VIB(Lon_,iParticle+1,iBlock))
    CosLat1xCosLat2 = 0.5 * (CosLat1MinusLat2 + CosLat1PlusLat2)
    SinLat1xSinLat2 = 0.5 * (CosLat1MinusLat2 - CosLat1PlusLat2)
    Distance = &
         State_VIB(R_,iParticle,  iBlock)**2 + &
         State_VIB(R_,iParticle+1,iBlock)**2 - &
         2*State_VIB(R_,iParticle,iBlock)*State_VIB(R_,iParticle+1,iBlock) * &
         (CosLat1xCosLat2 * CosLon1MinusLon2 + SinLat1xSinLat2)
  end function distance_to_next
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
    iNodeOut = iNode_II(iLon, iLat)
  end subroutine get_node

  !============================================================================

  subroutine convert_to_hgi(CoordIn_D, CoordOut_D)
    real, intent(in) :: CoordIn_D(nDim)
    real, intent(out):: CoordOut_D(nDim)

    integer:: iBlock, iCell_D(nDim)
    !--------------------------------------------------------------------------
    call get_cell(CoordIn_D, iCell_D)
    iBlock = &
         iGridGlobal_IA(Block_, iNode_II(iCell_D(OriginLon_), iCell_D(OriginLat_)))
    CoordOut_D((/R_, Lat_, Lon_/)) = &
         State_VIB((/R_,Lat_,Lon_/), iCell_D(Particle_), iBlock)
  end subroutine convert_to_hgi

  !============================================================================


end module ModGrid
