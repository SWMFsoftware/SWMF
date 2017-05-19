module SP_ModGrid

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nMomentumBin, nPitchAngleBin, &
       iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_

  implicit none

  SAVE

  private ! except

  public:: set_grid_param, init_grid, get_node_indexes, distance_to_next
  public:: fix_grid_consistency
  public:: iComm, iProc, nProc, nBlock, Proc_, Block_
  public:: LatMin, LatMax, LonMin, LonMax, RMin, RSc, RMax, ROrigin
  public:: iGridGlobal_IA, iGridLocal_IB, iNode_II, iNode_B
  public:: CoordMin_DI
  public:: State_VIB, Distribution_IIB
  public:: MomentumScale_I, LogMomentumScale_I, EnergyScale_I, LogEnergyScale_I
  public:: DMomentumOverDEnergy_I
  public:: Begin_, End_, Shock_, ShockOld_
  public:: nVar, X_, Y_, Z_, D_, S_
  public:: Rho_, T_, Ux_,Uy_,Uz_,U_, Bx_,By_,Bz_,B_, RhoOld_, BOld_, EFlux_
  public:: NameVar_V
  public:: TypeCoordSystem

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
  ! Starting position of field lines in Rs
  real:: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin surface R=ROrigin
  real:: LatMin, LatMax, DLat
  real:: LonMin, LonMax, DLon
  ! Lower boundary of the domain in Rs
  real:: RMin=-1.
  ! Boundary of the solar corona in Rs
  real:: RSc =-1.
  ! Upper boundary of the domain in Rs
  real:: RMax=-1.
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
  ! Foot-prints of the traced field lines on the surface R=RMin;
  ! 1st index - three spherical coordinates (R is added for completeness)
  ! 2nd index - block number
  real, allocatable:: CoordMin_DI(:,:)
  ! the initial length of segment 1-2: to control its as new particles
  ! are appended to the beginnings of lines
  real, allocatable:: Length_I(:)
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
  integer, parameter:: nVar = 18
  integer, parameter:: &
       X_     = 1, & ! 
       Y_     = 2, & ! Cartesian coordinates
       Z_     = 3, & ! 
       D_     = 4, & ! Distance to the next particle
       S_     = 5, & ! Distance from the beginning of the line
       ! Current values
       Rho_   = 6, & ! Background plasma density
       T_     = 7, & ! Background temperature
       Ux_    = 8, &
       Uy_    = 9, &
       Uz_    =10, &
       U_     =11, &
       Bx_    =12, & !
       By_    =13, & ! Background magnetic field
       Bz_    =14, & !
       B_     =15, & ! Magnitude of magnetic field
       ! Old values
       RhoOld_=16, & ! Background plasma density
       BOld_  =17, & ! Magnitude of magnetic field
       ! derived values
       EFlux_ =18    ! Integrated particle flux
  ! variable names
  character(len=10), parameter:: NameVar_V(nVar) = (/&
       'X     ', &
       'Y     ', &
       'Z     ', &
       'D     ', &
       'S     ', &
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
       'BOld  ', &
       'EFlux '  /)
  !----------------------------------------------------------------------------
  ! Distribution vector;
  ! Number of bins in the distribution is set in ModSize
  ! 1st index - log(momentum) bin
  ! 2nd index - particle index along the field line
  ! 4th index - local block number
  real, allocatable:: Distribution_IIB(:,:,:)
  ! scale with respect to Momentum and log(Momentum)
  real, target:: MomentumScale_I(nMomentumBin)
  real, target:: LogMomentumScale_I(nMomentumBin)
  real, target:: EnergyScale_I(nMomentumBin)
  real, target:: LogEnergyScale_I(nMomentumBin)
  real, target:: DMomentumOverDEnergy_I(nMomentumBin)
  !----------------------------------------------------------------------------
  ! Coordinate system and geometry
  character(len=3) :: TypeCoordSystem = 'HGI'
  !/

contains
  
  subroutine set_grid_param
    use ModReadParam, ONLY: read_var
    use ModNumConst, ONLY: cPi
    character(len=*), parameter:: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    call read_var('ROrigin', ROrigin)
    call read_var('RSc', RSc)
    call read_var('RMax',RMax)
    if(ROrigin < 0.0 .or. RSc < 0.0 .or. RMax < 0.0)&
         call CON_stop(NameSub//&
         ': all values ROrigin, RSc, RMax msut be set to positive values')
    if(RMax < RSc .or. RMax < ROrigin)&
         call CON_stop(NameSub//&
         ': value of RMax is inconsistent with ROrigin or RSc')
    call read_var('LonMin', LonMin)
    call read_var('LonMax', LonMax)
    if(LonMax <= LonMin)&
         call CON_stop(NameSub//': Origin surface grid is inconsistent')
    ! convert angels from degrees to radians
    LonMax = LonMax * cPi / 180
    LonMin = LonMin * cPi / 180
    ! angular grid's step
    DLon = (LonMax - LonMin) / nLon

    call read_var('LatMin', LatMin)
    call read_var('LatMax', LatMax)
    if(LatMax <= LatMin)&
         call CON_stop(NameSub//': Origin surface grid is inconsistent')
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
    use ModCoordTransform, ONLY: rlonlat_to_xyz
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
         call CON_stop(NameSub//': There are more processors than field lines')
    nBlock = ((iProc+1)*nNode) / nProc - (iProc*nNode) / nProc
    !\
    ! check consistency
    !/
    if(nLat <= 0 .or. nLon <= 0)&
         call CON_stop(NameSub//': Origin surface grid is invalid')
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
    allocate(CoordMin_DI(nDim, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'CoordMin_DI')
    allocate(Length_I(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'Length_I')
    allocate(State_VIB(nVar,iParticleMin:iParticleMax,nBlock), stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    allocate(Distribution_IIB(&
         nMomentumBin,iParticleMin:iParticleMax,nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'Distribution_IIB')
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
             iGridLocal_IB(Begin_,   iBlock) = 1
             iGridLocal_IB(End_,     iBlock) = 1
             iGridLocal_IB(Shock_,   iBlock) = iParticleMin - 1
             iGridLocal_IB(ShockOld_,iBlock) = iParticleMin - 1
          end if
          iGridGlobal_IA(Proc_,   iNode)  = iProcNode
          iGridGlobal_IA(Block_,  iNode)  = iBlock
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
    Distribution_IIB = tiny(1.0)
    State_VIB = -1
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iNode_II(iLon, iLat)
          CoordOrigin_DA(:, iNode) = &
               (/ROrigin, LonMin + (iLon-0.5)*DLon, LatMin + (iLat-0.5)*DLat/)
          iBlock = iGridGlobal_IA(Block_, iNode)
          if(iProc == iGridGlobal_IA(Proc_, iNode))&
               call rlonlat_to_xyz(&
               CoordOrigin_DA(:,iNode), State_VIB(X_:Z_,1,iBlock))
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
          ! if particle has left the domain -> cut the rest of the line
          if(sum(State_VIB(X_:Z_, iParticle, iBlock)**2) > RMax**2)then
             iGridLocal_IB(End_,  iBlock) = iParticle - 1
             EXIT
          end if
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
          
          ! distance from the beginning of the line
          if(iParticle == iGridLocal_IB(Begin_,  iBlock))then
             State_VIB(S_, iParticle, iBlock) = 0.0
          else
             State_VIB(S_, iParticle, iBlock) = &
                  State_VIB(S_, iParticle-1, iBlock) + &
                  State_VIB(D_, iParticle-1, iBlock)
          end if
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
    !  Distance**2 = R1**2 + R2**2 - 
    !     2*R1*R2*(Cos(Lat1)*Cos(Lat2) * Cos(Lon1-Lon2) + Sin(Lat1)*Sin(Lat2))
    ! NOTE: function doesn't check whether iParticle is last on the field line
    integer, intent(in):: iParticle
    integer, intent(in):: iBlock
    real               :: Distance
    !--------------------------------------------------------------------
    Distance = sqrt(sum((&
         State_VIB(X_:Z_, iParticle,   iBlock) - &
         State_VIB(X_:Z_, iParticle+1, iBlock))**2))
  end function distance_to_next

  !============================================================================

  subroutine append_particles
    !appends a new particle at the beginning of lines if necessary
    integer:: iBlock
    real:: DistanceToMin, Alpha

    character(len=*), parameter:: NameSub = 'append_particles'
    !--------------------------------------------------------------------
    do iBlock = 1, nBlock
       ! check if the beginning of the line moved far enough from its 
       ! footprint on the solar surface
       DistanceToMin = sqrt(sum((&
            State_VIB(X_:Z_, 1, iBlock) - CoordMin_DI(:,iBlock))**2))
       ! skip the line if it's still close to the Sun
       if(DistanceToMin <= Length_I(iBlock)) CYCLE
       
       ! append a new particle
       !-----------------------
       ! check if have enough space
       if(nParticle == iGridLocal_IB(End_, iBlock))&
            call CON_Stop(NameSub//&
            ': not enough memory allocated to append a new particle')
       ! shift the grid:
       State_VIB(     :,2:iGridLocal_IB(End_, iBlock)+1, iBlock) = &
            State_VIB(:,1:iGridLocal_IB(End_, iBlock),   iBlock)
       Distribution_IIB(     :,2:iGridLocal_IB(End_, iBlock)+1, iBlock) = &
            Distribution_IIB(:,1:iGridLocal_IB(End_, iBlock),   iBlock)
       iGridLocal_IB(End_, iBlock) = iGridLocal_IB(End_, iBlock) + 1
       ! compute new coordinates
       ! TO BE CHANGED: NOW NEW PARTICLE IS ON THE LINE CONNECTING 
       ! THE CURRENT BEGINNING WITH MIN 
       Alpha = Length_I(iBlock) / DistanceToMin
       State_VIB(X_:Z_, 1, iBlock) = Alpha * CoordMin_DI(:,iBlock) + &
            (1 - Alpha) * State_VIB(X_:Z_, 2, iBlock)
    end do
  end subroutine append_particles

end module SP_ModGrid
