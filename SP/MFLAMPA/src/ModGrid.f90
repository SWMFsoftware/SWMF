module SP_ModGrid

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nMomentumBin, nPitchAngleBin, &
       iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_

  implicit none

  SAVE

  private ! except

  public:: set_grid_param, init_grid, get_node_indexes, distance_to_next
  public:: fix_grid_consistency, append_particles
  public:: iComm, iProc, nProc, nBlock, Proc_, Block_, nBlockIndexes, nBlockParam
  public:: LatMin, LatMax, LonMin, LonMax
  public:: RMin, RBufferMin, RBufferMax, RMax, ROrigin
  public:: iGridGlobal_IA, iGridLocal_IB, ParamLocal_IB, iNode_II, iNode_B
  public:: State_VIB, Distribution_IIB
  public:: MomentumScale_I, LogMomentumScale_I, EnergyScale_I, LogEnergyScale_I
  public:: DMomentumOverDEnergy_I
  public:: Begin_, End_, Shock_, ShockOld_, XMin_, YMin_, ZMin_, Length_
  public:: nVar, nVarRead,  X_, Y_, Z_, D_, S_, LagrID_, Offset_
  public:: Rho_,T_, Ux_,Uy_,Uz_,U_,DLogRho_, Bx_,By_,Bz_,B_, RhoOld_,BOld_
  public:: EFlux_, Flux0_, Flux1_, Flux2_, Flux3_, Flux4_, Flux5_, Flux6_
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
  ! Upper boundary of the domain in Rs
  real:: RMax=-1.
  ! Boundaries of the buffer layer between SC and IH Rs
  real:: RBufferMin=-1.
  real:: RBufferMax=-1.
  ! Mark that grid or lines' origin have been set
  logical:: IsSetGrid   = .false.
  logical:: IsSetOrigin = .false.
  !----------------------------------------------------------------------------
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
  real,    allocatable:: ParamLocal_IB(:,:)
  !----------------------------------------------------------------------------
  ! Number of info fields per node/block and their identifications
  integer, parameter:: nNodeIndexes = 2
  integer, parameter:: &
       Proc_  = 1, & ! Processor that has this line/node
       Block_ = 2    ! Block that has this line/node
  integer, parameter:: nBlockIndexes = 5
  integer, parameter:: &
       Begin_   = 1, & ! Index of the 1st particle on this line/node
       End_     = 2, & ! Index of the last particle on this line/node
       Shock_   = 3, & ! Current location of a shock wave
       ShockOld_= 4, & ! Old location of a shock wave
       Offset_  = 5    ! To account for the dymaical grid distinction 
                       ! from that updated in the other components
  
  integer, parameter:: nBlockParam = 4
  integer, parameter:: &
       XMin_   = 1, & ! 
       YMin_   = 2, & ! Foot-prints of the field lines on the surface R=RMin;
       ZMin_   = 3, & ! 
       Length_ = 4    ! init length of segment 1-2: control for new particles
                      ! being appended to the beginnings of lines
  !----------------------------------------------------------------------------
  ! State vector;
  ! 1st index - identification of variable
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, allocatable:: State_VIB(:,:,:)
  !----------------------------------------------------------------------------
  ! Number of variables in the state vector and their identifications
  integer, parameter:: nVar     = 27
  integer, parameter:: nVarRead = 12
  integer, parameter:: &
       !\
       !-- The following variables MUST be in CONTIGUOUS  order --------------
       !-- used in subroutines read_mh_data, write_restart, read_restart -----
       !-- DO NOT CHANGE WITHOUT CAREFULL CONSIDERATION !!! ------------------
       LagrID_ = 1, & ! Lagrangian id
       X_      = 2, & ! 
       Y_      = 3, & ! Cartesian coordinates
       Z_      = 4, & ! 
       Rho_    = 5, & ! Background plasma density
       T_      = 6, & ! Background temperature
       Ux_     = 7, & !
       Uy_     = 8, & ! Background plasma bulk velocity
       Uz_     = 9, & !
       Bx_     =10, & !
       By_     =11, & ! Background magnetic field
       Bz_     =12, & !
       !-----------------------------------------------------------------------
       D_      =13, & ! Distance to the next particle
       S_      =14, & ! Distance from the beginning of the line
       U_      =15, & ! Magnitude of plasma bulk velocity
       B_      =16, & ! Magnitude of magnetic field
       DLogRho_=17, & ! Dln(Rho), i.e. -div(U) * Dt
       RhoOld_ =18, & ! Background plasma density
       BOld_   =19, & ! Magnitude of magnetic field
       Flux0_  =20, & ! Total integral (simulated) particle flux
       Flux1_  =21, & ! Integral particle flux >  5 MeV (GOES Channel 1)
       Flux2_  =22, & ! Integral particle flux > 10 MeV (GOES Channel 2)
       Flux3_  =23, & ! Integral particle flux > 30 MeV (GOES Channel 3)
       Flux4_  =24, & ! Integral particle flux > 50 MeV (GOES Channel 4)
       Flux5_  =25, & ! Integral particle flux > 60 MeV (GOES Channel 5)
       Flux6_  =26, & ! Integral particle flux >100 MeV (GOES Channel 6)
       EFlux_  =27    ! Total integral energy flux

  ! variable names
  character(len=10), parameter:: NameVar_V(nVar) = (/&
       'LagrID    ', &
       'X         ', &
       'Y         ', &
       'Z         ', &
       'Rho       ', &
       'T         ', &
       'Ux        ', &
       'Uy        ', &
       'Uz        ', &
       'Bx        ', &
       'By        ', &
       'Bz        ', &
       'D         ', &
       'S         ', &
       'U         ', &
       'B         ', &
       'DLogRho   ', &
       'RhoOld    ', &
       'BOld      ', &
       'Flux_Total', &
       'Flux_GOES1', &
       'Flux_GOES2', &
       'Flux_GOES3', &
       'Flux_GOES4', &
       'Flux_GOES5', &
       'Flux_GOES6', &
       'EFlux     '  /)
  !----------------------------------------------------------------------------
  ! Distribution vector;
  ! Number of bins in the distribution is set in ModSize
  ! 1st index - log(momentum) bin
  ! 2nd index - particle index along the field line
  ! 4th index - local block number
  real, allocatable:: Distribution_IIB(:,:,:)
  ! scale with respect to Momentum and log(Momentum)
  real:: MomentumScale_I(nMomentumBin)
  real:: LogMomentumScale_I(nMomentumBin)
  real:: EnergyScale_I(nMomentumBin)
  real:: LogEnergyScale_I(nMomentumBin)
  real:: DMomentumOverDEnergy_I(nMomentumBin)
  !----------------------------------------------------------------------------
  ! Coordinate system and geometry
  character(len=3) :: TypeCoordSystem = 'HGR'
  !/

contains
  
  subroutine set_grid_param(TypeAction)
    use ModReadParam, ONLY: read_var
    use ModNumConst, ONLY: cPi
    character (len=*), intent(in):: TypeAction ! What to do  
    character(len=*), parameter:: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    select case(TypeAction)
    case('#ORIGIN')
       call read_var('ROrigin', ROrigin)
       call read_var('LonMin', LonMin)
       call read_var('LatMin', LatMin)
       call read_var('LonMax', LonMax)
       call read_var('LatMax', LatMax)

       ! check consistency
       if(LonMax <= LonMin .or. LatMax <= LatMin)&
            call CON_stop(NameSub//': Origin surface grid is inconsistent')
       if(ROrigin < 0.0)&
            call CON_stop(NameSub//&
            ': ROrigin, if set, must have a positive values')
       if(any((/RBufferMin, RBufferMax, RMax/) < ROrigin) .and. IsSetGrid)&
            call CON_stop(NameSub//&
            ': inconsistent values of ROrigin, RBufferMin, RBufferMax, RMax')

       ! convert angels from degrees to radians
       LonMax = LonMax * cPi / 180
       LonMin = LonMin * cPi / 180
       ! angular grid's step
       DLon = (LonMax - LonMin) / nLon

       ! convert angels from degrees to radians
       LatMax = LatMax * cPi / 180
       LatMin = LatMin * cPi / 180
       ! angular grid's step
       DLat = (LatMax - LatMin) / nLat

       IsSetOrigin = .true.
    case('#GRID')
       call read_var('RMin',RMin)
       call read_var('RBufferMin', RBufferMin)
       call read_var('RBufferMax', RBufferMax)
       call read_var('RMax',RMax)

       ! check consistency
       if(RBufferMin < 0.0 .or.RBufferMax < 0.0 .or.RMax < 0.0)&
            call CON_stop(NameSub//&
            ': RBufferMin, RBufferMax, RMax must be set to positive values')
       if(any((/RMax, RBufferMax/) < RBufferMin) .or. RMax < RBufferMax .or. &
            any((/RBufferMin, RBufferMax, RMax/) < ROrigin) .and. IsSetOrigin)&
            call CON_stop(NameSub//&
            ': inconsistent values of ROrigin, RBufferMin, RBufferMax, RMax')

       IsSetGrid = .true.
    end select
  end subroutine set_grid_param

  !============================================================================

  subroutine init_grid(DoReadInput)
    ! allocate the grid used in this model
    use ModUtilities, ONLY: check_allocate
    use ModCoordTransform, ONLY: rlonlat_to_xyz

    logical, intent(in):: DoReadInput
    integer:: iError
    integer:: iLat, iLon, iNode, iBlock, iProcNode
    character(LEN=*),parameter:: NameSub='SP:init_grid'
    !--------------------------------------------------------------------------
    !\
    ! Check if everything's ready for initialization
    !/
    if(.not.IsSetGrid)&
         call CON_stop(NameSub//': grid is not set in PARAM.in file')
    if(.not.IsSetOrigin .and. .not. DoReadInput)&
         call CON_stop(NameSub//": neither lines' origin is set, "//&
         "nor input files are provided; change PARAM.in file!!!")
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
    allocate(ParamLocal_IB(nBlockParam, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'ParamLocal_IB')
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
             iGridLocal_IB(Offset_,  iBlock) = 0
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

    if(DoReadInput)then
       if(IsSetOrigin)&
            write(*,*)NameSub//": input files are provided, "//&
            "but lines' origin is set in PARAM.in. "//&
            "The latter will be IGNORED!!!"
       RETURN
    end if
    
    do iLat = 1, nLat
       do iLon = 1, nLon
          iNode = iNode_II(iLon, iLat)
          iBlock = iGridGlobal_IA(Block_, iNode)
          if(iProc == iGridGlobal_IA(Proc_, iNode))then
             call rlonlat_to_xyz(&
                  (/ROrigin, LonMin+(iLon-0.5)*DLon, LatMin+(iLat-0.5)*DLat/),&
                  State_VIB(X_:Z_,1,iBlock))
             State_VIB(LagrID_,1,iBlock) = 1
          end if
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

          ! divergence of plasma velocity
          State_VIB(DLogRho_,iParticle,iBlock) = log(&
               State_VIB(Rho_,iParticle,iBlock) / &
               State_VIB(RhoOld_,iParticle,iBlock))

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
    real, parameter:: cTol = 1E-06

    character(len=*), parameter:: NameSub = 'append_particles'
    !--------------------------------------------------------------------
    do iBlock = 1, nBlock
       ! check if the beginning of the line moved far enough from its 
       ! footprint on the solar surface
       DistanceToMin = sqrt(sum((&
            State_VIB(X_:Z_,1,iBlock) - ParamLocal_IB(XMin_:ZMin_,iBlock))**2))
       ! skip the line if it's still close to the Sun
       if(DistanceToMin * (1.0 + cTol) < ParamLocal_IB(Length_, iBlock)) CYCLE
       
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
       iGridLocal_IB(End_, iBlock)  = iGridLocal_IB(End_, iBlock) + 1
       !Particles ID as handled by other components keep unchanged
       !while their order numbers in SP are increased by 1. Therefore,
       iGridLocal_IB(Offset_, iBlock)  = iGridLocal_IB(Offset_, iBlock) + 1
       ! put the new particle just above the lower boundary
       State_VIB(X_:Z_,  1, iBlock) = &
            ParamLocal_IB(XMin_:ZMin_, iBlock) * (1.0 + cTol)
       State_VIB(LagrID_,1, iBlock) = State_VIB(LagrID_, 2, iBlock) - 1.0
       ! for old values of background parameters use extrapolation
       Alpha = DistanceToMin / (DistanceToMin + State_VIB(D_, 2, iBlock))
       State_VIB((/RhoOld_, BOld_/), 1, iBlock) = &
            (Alpha+1) * State_VIB((/RhoOld_, BOld_/), 2, iBlock) &
            -Alpha    * State_VIB((/RhoOld_, BOld_/), 3, iBlock)
       
    end do
  end subroutine append_particles

end module SP_ModGrid
