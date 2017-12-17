!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModGrid
  !Multi-line grid, D.Borovikov & I.Sokolov, Dec,17, 2017.
  !Further revision history:
  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, nMomentumBin, nPitchAngleBin, &
       nParticleMax, Particle_, OriginLat_, OriginLon_

  implicit none

  SAVE

  private ! except
  !Public members:
  public:: set_grid_param  !read parameters related to grid 
  public:: init_grid       !Initialize arrays on the grid
  public:: copy_old_state  !save old arrays before getting new ones  
  public:: get_other_state_var !Auxiliary components of state vector 
 
  ! Coordinate system and geometry
  character(len=3), public :: TypeCoordSystem = 'HGR'
  !\
  ! Grid info
  !/
  !\
  ! Grid size, boundaries, coordinates
  ! Starting position of field lines in Rs
  real, public :: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin 
  ! surface R=ROrigin
  real, public ::  LonMin, LonMax, LatMin, LatMax
  !Sell size on the origin surface, per line
  real         ::  DLon, DLat
  ! Lower/Upper boundary of the domain in Rs
  real, public :: RMin=-1.0, RMax = -1.0
  ! Boundaries of the buffer layer between SC and IH Rs
  real, public :: RBufferMin=-1.0, RBufferMax=-1.0
  ! Number of blocks on this processor
  integer, public              :: nBlock
  ! Number of particles per block (line):
  integer, public, allocatable :: nParticle_B(:)
  !/
  !\
  ! Node number based on the field line identified by 
  ! 2 angular grid indices, latitude and longitude;
  ! 1st index - longitude index
  ! 2nd index - latitude index
  integer, public, allocatable:: iNode_II(:,:)
  !inverse function:get_node_indexes(iNodeIn,iLonOut,iLatOut)
  public :: get_node_indexes
  !/
  !\
  ! Node number based on the local block number
  ! 1st index - block number
  integer, public, allocatable:: iNode_B(:)
  !/
  !\
  ! Various house-keeping information about the node/line;
  ! 1st index - identification of info field
  ! 2nd index - node number / block number
  !/
  !\
  !Proc_ and Block_ number, for a given node:
  integer, public, parameter:: &
       Proc_  = 1, & ! Processor that has this line/node
       Block_ = 2    ! Block that has this line/node
    integer, public, allocatable:: iGridGlobal_IA(:,:)
  !/
  !\
  ! Array for current and present location of shock wave
  integer, public, parameter:: nShockParam = 2,  &
       Shock_   = 1, & ! Current location of a shock wave
       ShockOld_= 2    ! Old location of a shock wave
  integer, public, allocatable:: iShock_IB(:,:)
  !/ 
  !\
  ! Information about the magnetic field line foot point:
  ! the Lagrangian (0) and Cartesian (1:3) coordinates, and
  integer, public, parameter :: &! init length of segment 1-2: control 
       Length_ = 4               ! appending  new particles 
  real, public, allocatable:: FootPoint_VB(:,:)
  !/
  !\
  ! State vector;
  ! 1st index - identification of variable
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, public, allocatable:: State_VIB(:,:,:)
  real, public, allocatable:: Flux_VIB( :,:,:)
  ! Number of variables in the state vector and their identifications
  integer, public, parameter :: nVarRead = 13, nVar = 29, FluxMax_ = 29,&
       !\
       LagrID_ = 0, & ! Lagrangian id           ^saved/     ^set to 0 in
       X_      = 1, & !                         |read in    |copy_old_state
       Y_      = 2, & ! Cartesian coordinates   |restart    |
       Z_      = 3, & !                         v/          |saved to 
       Rho_    = 4, & ! Background plasma density           |mhd1
       T_      = 5, & ! Background temperature              | 
       Ux_     = 6, & !                                     |may be
       Uy_     = 7, & ! Background plasma bulk velocity     |read from
       Uz_     = 8, & !                                     |mhd1
       Bx_     = 9, & !                                     |or
       By_     =10, & ! Background magnetic field           |received 
       Bz_     =11, & !                                     |from
       Wave1_  =12, & !\                                    |coupler
       Wave2_  =13, & ! Alfven wave turbulence              v
       !-----------------------------------------------------------------------
       D_      =14, & ! Distance to the next particle   ^
       S_      =15, & ! Distance from the foot point    |derivatives from
       R_      =16, & ! Heliocentric distance           |the mhd1, set in
       U_      =17, & ! Plasma speed                    |get_other_state_var
       B_      =18, & ! Magnitude of magnetic field     v
       DLogRho_=19, & ! Dln(Rho), i.e. -div(U) * Dt
       RhoOld_ =20, & ! Background plasma density       !\copy_old_state
       BOld_   =21, & ! Magnitude of magnetic field     !/
       Flux0_  =22, & ! Total integral (simulated) particle flux
       Flux1_  =23, & ! Integral particle flux >  5 MeV (GOES Channel 1)
       Flux2_  =24, & ! Integral particle flux > 10 MeV (GOES Channel 2)
       Flux3_  =25, & ! Integral particle flux > 30 MeV (GOES Channel 3)
       Flux4_  =26, & ! Integral particle flux > 50 MeV (GOES Channel 4)
       Flux5_  =27, & ! Integral particle flux > 60 MeV (GOES Channel 5)
       Flux6_  =28, & ! Integral particle flux >100 MeV (GOES Channel 6)
       EFlux_  =29    ! Total integral energy flux
  !/
  !\
  ! variable names
  character(len=10), public, parameter:: NameVar_V(LagrID_:EFlux_) = (/&
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
       'Wave1     ', &
       'Wave2     ', &
       'D         ', &
       'S         ', &
       'R         ', &
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
  !/
  !\
  ! Velosity Distribution Function (VDF) 
  ! Number of bins in the distribution is set in ModSize
  ! 1st index - log(momentum) bin
  ! 2nd index - particle index along the field line
  ! 3rd index - local block number
  real, public, allocatable:: Distribution_IIB(:,:,:)
  !/
  !Misc
  ! Mark that grid or lines' origin have been set
  logical:: IsSetGrid   = .false.
  logical:: IsSetOrigin = .false.
contains  
  subroutine set_grid_param(NameCommand)
    use ModReadParam, ONLY: read_var
    use ModNumConst, ONLY : cDegToRad
    character (len=*), intent(in):: NameCommand ! From PARAM.in  
    character(len=*), parameter  :: NameSub = 'SP:set_grid_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
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
       LonMax = LonMax*cDegToRad
       LonMin = LonMin*cDegToRad
       ! angular grid's step
       DLon = (LonMax - LonMin) / nLon

       ! convert angels from degrees to radians
       LatMax = LatMax*cDegToRad
       LatMin = LatMin*cDegToRad
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

  !==========================================================================
  subroutine init_grid(DoReadInput)
    ! allocate the grid used in this model
    use ModUtilities,      ONLY: check_allocate
    use ModCoordTransform, ONLY: rlonlat_to_xyz
    use SP_ModProc,        ONLY: iProc, nProc

    logical, intent(in):: DoReadInput
    integer:: iError
    integer:: iLat, iLon, iNode, iBlock, iProcNode, iParticle
    character(LEN=*),parameter:: NameSub='SP:init_grid'
    !-------------------------------------------------------------------------
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
    allocate(iGridGlobal_IA(Proc_:Block_, nNode), stat=iError)
    call check_allocate(iError, NameSub//'iGridGlobal_IA')
    allocate(nParticle_B(nBlock), stat=iError)
    call check_allocate(iError, NameSub//'nParticle_B')
    allocate(iShock_IB(nShockParam, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'iShock_IB')
    allocate(FootPoint_VB(LagrID_:Length_, nBlock), stat=iError)
    call check_allocate(iError, NameSub//'FootPoint_VB')
    allocate(State_VIB(LagrID_:nVar,1:nParticleMax,nBlock), &
         stat=iError)
    call check_allocate(iError, NameSub//'State_VIB')
    allocate(Flux_VIB(Flux0_:FluxMax_,1:nParticleMax,nBlock), &
         stat=iError); call check_allocate(iError, 'Flux_VIB')
    allocate(Distribution_IIB(&
         nMomentumBin,1:nParticleMax,nBlock), &
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
             nParticle_B(     iBlock) = 1
             iShock_IB(:, iBlock) = 0
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
    State_VIB = -1; Flux_VIB = -1; FootPoint_VB = -1
    
    !\
    ! reset lagrangian ids
    !/
    do iParticle = 1, nParticleMax
       State_VIB(LagrID_, iParticle, 1:nBlock) = real(iParticle)
    end do

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
          end if
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
  !===================
  subroutine copy_old_state
    ! copy current state to old state for all field lines
    integer:: iEnd, iBlock
    !-----------------------------------------------------------------------
    do iBlock = 1, nBlock
       iEnd   = nParticle_B(  iBlock)
       State_VIB((/RhoOld_,BOld_/), 1:iEnd, iBlock) = &
            State_VIB((/Rho_,B_/),  1:iEnd, iBlock)
       ! reset variables to be read from file or received via coupler
       State_VIB(1:nVarRead,1:iEnd, iBlock) = 0.0
    end do
  end subroutine copy_old_state
  !===================
  subroutine get_other_state_var
    integer:: iBlock, iParticle, iEnd
    !----------------------------------------------------------------------
    do iBlock = 1, nBlock
       iEnd   = nParticle_B(  iBlock)
       do iParticle = 1, iEnd
          ! plasma speed
          State_VIB(U_,iParticle, iBlock) = &
               sqrt(sum(State_VIB(Ux_:Uz_,iParticle,iBlock)**2))
          ! magnetic field
          State_VIB(B_,iParticle, iBlock) = &
               sqrt(sum(State_VIB(Bx_:Bz_,iParticle,iBlock)**2))
          ! distances between particles
          if(iParticle /=nParticle_B(iBlock))&
               State_VIB(D_, iParticle, iBlock) = sqrt(sum((&
               State_VIB(X_:Z_, iParticle    , iBlock) - &
               State_VIB(X_:Z_, iParticle + 1, iBlock))**2))
          ! distance from the beginning of the line
          if(iParticle == 1)then
             State_VIB(S_, iParticle, iBlock) = 0.0
          else
             State_VIB(S_, iParticle, iBlock) = &
                  State_VIB(S_, iParticle-1, iBlock) + &
                  State_VIB(D_, iParticle-1, iBlock)
          end if
          !Heliocentric Distance
          State_VIB(R_, iParticle, iBlock) = &
               sqrt(sum(State_VIB(X_:Z_, iParticle, iBlock)**2, 1))
       end do
    end do
  end subroutine get_other_state_var
end module SP_ModGrid
