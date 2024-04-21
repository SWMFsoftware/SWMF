!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModGrid

  ! Multi-line grid, D.Borovikov & I.Sokolov, Dec,17, 2017.
  ! Dec.23 2017: exclude fluxes from the state vector.
  ! Dec.25 2017: standard init and read_param
  ! Dec.25 2017: rename nVarRead=>nMhData, add NoShock_ param.

#ifdef OPENACC
  use ModUtilities, ONLY: norm2
#endif
  use ModUtilities, ONLY: CON_stop
  use PT_ModSize,  ONLY: nVertexMax
  use PT_ModProc,  ONLY: iProc
  use ModNumConst, ONLY: cTwoPi, cPi

  implicit none
  SAVE

  PRIVATE ! except
  ! Public members:
  public:: read_param          ! read parameters related to grid
  public:: init                ! Initialize arrays on the grid
  public:: init_stand_alone    ! Initialize arrays on the grid
  public:: copy_old_state      ! save old arrays before getting new ones
  public:: get_other_state_var ! Auxiliary components of state vector
  public:: get_shock_location  ! finds shock location on all lines
  public:: search_line         ! find particle index corresponding to radius

  ! Coordinate system and geometry
  character(len=3), public :: TypeCoordSystem = 'HGR'
  !
  ! Grid info
  !
  ! Angular grid at origin surface
  !
  integer, public :: nLon  = 4
  integer, public :: nLat  = 4
  !
  ! Total number of magnetic field lines on all PEs
  ! (just a product of nLat * nLon)
  integer, public :: nLineAll = 16
  !
  ! All nodes are enumerated. The last node number on the previous proc
  ! (iProc-1)
  ! equals (iProc*nLineAll)/nProc. Store this:
  !
  integer, public :: iLineAll0
  !
  ! The nodes on a given PE have node numbers ranging from iLineAll0 +1 to
  ! iNodeLast =((iProc + 1)*nLineAll)/nProc. The iLine index to enumerate
  ! lines on a given proc ranges from 1 to iNodeLast.
  ! nLine = nNodeLast - iLineAll0 is the number of
  ! lines (blocks) on this processor. For iLine=1:nLine
  ! iLineAll = iLineAll0+1:iNodeLast
  !
  integer, public :: nLine
  !
  ! Number of particles (vertexes, Lagrangian meshes) per line (line):
  integer, public,     pointer :: nVertex_B(:)
  !
  ! Function converting line number to lon-lat location of the line
  !
  public :: iBlock_to_lon_lat

  !
  ! Array for current and present location of shock wave
  !
  integer, public, parameter:: nShockParam = 2,  &
       Shock_   = 1, & ! Current location of a shock wave
       ShockOld_= 2    ! Old location of a shock wave
  integer, public, allocatable:: iShock_IB(:,:)
  integer, public, parameter:: NoShock_ = 1
  !
  ! Information about the magnetic field line foot point:
  ! the Lagrangian (0) and Cartesian (1:3) coordinates, and
  integer, public, parameter :: &! init length of segment 1-2:
       Length_ = 4               ! control appending  new particles
  real, public, pointer :: FootPoint_VB(:,:)
  !
  ! Logical to mark unusable lines
  !
  logical, public, pointer :: Used_B(:)
  !
  ! MHD state vector;
  ! 1st index - identification of variable (LagrID_:Wave2_)
  ! 2nd index - particle index along the field line
  ! 3rd index - local line number
  !
  real, public, pointer     :: MhData_VIB(:,:,:)
  !
  ! Aux state vector;
  ! 1st index - identification of variable (D_:BOld_)
  ! 2nd index - particle index along the field line
  ! 3rd index - local line number
  !
  real, public, pointer     :: State_VIB(:,:,:)
  !
  ! Number of variables in the state vector and the identifications
  !
  integer, public, parameter :: nMhData = 13, nVar = 21,          &
       !
       LagrID_     = 0, & ! Lagrangian id           ^saved/   ^set to 0
       X_          = 1, & !                         |read in  |in copy_
       Y_          = 2, & ! Cartesian coordinates   |restart  |old_stat
       Z_          = 3, & !                         v/        |saved to
       Rho_        = 4, & ! Background plasma density         |mhd1
       T_          = 5, & ! Background temperature            |
       Ux_         = 6, & !                                   |may be
       Uy_         = 7, & ! Background plasma bulk velocity   |read from
       Uz_         = 8, & !                                   |mhd1
       Bx_         = 9, & !                                   |or
       By_         =10, & ! Background magnetic field         |received
       Bz_         =11, & !                                   |from
       Wave1_      =12, & !                                   |coupler
       Wave2_      =13, & !-Alfven wave turbulence            v
       !-
       R_          =14, & ! Heliocentric distance          ^derived from
       D_          =15, & ! Distance to the next particle  |MHdata in
       S_          =16, & ! Distance from the foot point   |get_other_
       U_          =17, & ! Plasma speed along field line  |state_var
       B_          =18, & ! Magnitude of magnetic field    v
       RhoOld_     =19, & ! Background plasma density      ! copy_
       UOld_       =20, & ! Background plasma bulk speed   ! old_
       BOld_       =21    ! Magnitude of magnetic field    ! state
  !
  ! variable names
  !
  character(len=10), public, parameter:: NameVar_V(LagrID_:nVar)&
    = ['LagrID    ', &
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
       'R         ', &
       'D         ', &
       'S         ', &
       'U         ', &
       'B         ', &
       'RhoOld    ', &
       'UOld      ', &
       'BOld      ' ]
  !
  logical:: DoInit = .true.
  ! Shock algorithm parameters:
  real,    public, parameter :: dLogRhoThreshold = 0.01
  integer, public, parameter :: nWidth = 50
contains
  !============================================================================
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    ! Misc
    integer :: nParticleCheck, nLonCheck, nLatCheck
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#CHECKGRIDSIZE')
       call read_var('nVertexMax',nParticleCheck)
       call read_var('nLon',     nLonCheck)
       call read_var('nLat',     nLatCheck)
       if(iProc==0.and.any(&
            [nLon,     nLat] /= &
            [nLonCheck,nLatCheck])&
            )write(*,'(a,2I5)') 'nLon,nLat are reset to ',&
            nLonCheck, nLatCheck
       nLon = nLonCheck
       nLat = nLatCheck
       nLineAll = nLon*nLat
       if(nParticleCheck > nVertexMax)then
          if(iProc==0)write(*,*)&
               'nVertexMax is too small, use ./Config.pl -g=',nParticleCheck
          call CON_stop('Code stopped')
       end if
    case('#COORDSYSTEM','#COORDINATESYSTEM')
       call read_var('TypeCoordSystem', TypeCoordSystem, &
            IsUpperCase=.true.)
    case('#GRIDNODE')
       call read_var('nLat',  nLat)
       call read_var('nLon',  nLon)
       nLineAll = nLat * nLon
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select

  end subroutine read_param
  !============================================================================
  subroutine init

    ! allocate the grid used in this model
    use ModUtilities,      ONLY: check_allocate
    use PT_ModProc,        ONLY: nProc

    integer:: iError
    integer:: iNodeLast

    character(len=*), parameter:: NameSub = 'init'
    !--------------------------------------------------------------------------
    if(.not.DoInit)RETURN
    DoInit = .false.
    !
    ! distribute nodes between processors
    !
    if(nLineAll < nProc)call CON_stop(NameSub//&
         ': There are more processors than field lines')
    iLineAll0 = (iProc*nLineAll)/nProc
    iNodeLast = ((iProc+1)*nLineAll)/nProc
    nLine = iNodeLast - iLineAll0
    !
    ! check consistency
    !
    if(nLat <= 0 .or. nLon <= 0)&
         call CON_stop(NameSub//': Origin surface grid is invalid')
    !
    ! allocate data and grid containers
    !
    allocate(iShock_IB(nShockParam, nLine), stat=iError)
    call check_allocate(iError, NameSub//'iShock_IB')
    iShock_IB = NoShock_

  end subroutine init
  !============================================================================
  subroutine init_stand_alone

    ! allocate the grid used in this model
    use ModUtilities,      ONLY: check_allocate
    integer :: iVertex, iError
    character(len=*), parameter:: NameSub = 'init_stand_alone'
    !--------------------------------------------------------------------------
    ! Allocate here if stand alone
    allocate(MhData_VIB(LagrID_:nMhData, 1:nVertexMax, nLine))
    !
    MhData_VIB(1:nMhData,:,:) = 0.0
    !
    ! reset lagrangian ids
    !
    do iVertex = 1, nVertexMax
       MhData_VIB(LagrID_, iVertex, 1:nLine) = real(iVertex)
    end do
    ! Allocate auxiliary State vector
    allocate(State_VIB(nMhData+1:nVar, 1:nVertexMax, nLine))
    State_VIB = -1
    allocate(nVertex_B(nLine))
    nVertex_B = 0
    allocate(FootPoint_VB(LagrID_:Length_, nLine))
    FootPoint_VB = -1
    allocate(Used_B(nLine)); Used_B = .true.

  end subroutine init_stand_alone
  !============================================================================
  subroutine iblock_to_lon_lat(iBlockIn, iLonOut, iLatOut)

    ! return angular grid's indexes corresponding to this line
    integer, intent(in) :: iBlockIn
    integer, intent(out):: iLonOut
    integer, intent(out):: iLatOut

    integer :: iLineAll
    !--------------------------------------------------------------------------
    !
    ! Get node number from line number
    !
    iLineAll = iBlockIn + iLineAll0
    iLatOut = 1 + (iLineAll - 1)/nLon
    iLonOut = iLineAll - nLon*(iLatOut - 1)

  end subroutine iblock_to_lon_lat
  !============================================================================
  subroutine copy_old_state

    ! copy current state to old state for all field lines
    integer:: iEnd, iLine
    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       if(.not.Used_B(iLine))CYCLE
       iEnd   = nVertex_B(  iLine)
       iShock_IB(ShockOld_,iLine) = iShock_IB(Shock_, iLine)
       State_VIB(RhoOld_, 1:iEnd, iLine) = &
            MhData_VIB(Rho_,  1:iEnd, iLine)
       State_VIB(UOld_, 1:iEnd, iLine) = &
            State_VIB(U_,  1:iEnd, iLine)
       State_VIB(BOld_, 1:iEnd, iLine) = &
            State_VIB(B_,  1:iEnd, iLine)
       ! reset variables read from file or received via coupler
       MhData_VIB(1:nMhData, 1:iEnd, iLine) = 0.0
    end do

  end subroutine copy_old_state
  !============================================================================
  subroutine get_other_state_var

    integer:: iLine, iVertex, iEnd
    integer:: iAux1, iAux2
    real   :: XyzAux1_D(x_:z_), XyzAux2_D(x_:z_)
    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       if(.not.Used_B(iLine))CYCLE
       iEnd   = nVertex_B(  iLine)
       do iVertex = 1, iEnd
          ! magnetic field
          State_VIB(B_,iVertex, iLine) = &
               norm2(MhData_VIB(Bx_:Bz_,iVertex,iLine))
          ! plasma speed (negative if velocity is antiparallel to
          ! magnetic field)
          State_VIB(U_,iVertex, iLine) = &
               sum(MhData_VIB(Ux_:Uz_,iVertex,iLine)*&
               MhData_VIB(Bx_:Bz_,iVertex,iLine))/   &
               State_VIB(B_,iVertex, iLine)
          ! distances between particles
          ! if(.not.DoSmooth)then
          if(iVertex /=nVertex_B(iLine))&
               State_VIB(D_, iVertex, iLine) = norm2(&
               MhData_VIB(X_:Z_, iVertex    , iLine) - &
               MhData_VIB(X_:Z_, iVertex + 1, iLine))
          ! else
          ! smoothing is done by groups:
          ! nSmooth particles are aggeregated into single effective one,
          ! find length increment between effective particles are used
          ! to find length increment between regular particles
          ! iAux1 = nSmooth * max(1, min(&
          !     iVertex/nSmooth,nVertex_B(iLine)/nSmooth-1))
          ! iAux2 = iAux1 + nSmooth
          ! XyzAux1_D = sum(MhData_VIB(&
          !     X_:Z_,iAux1-nSmooth+1:iAux1,iLine),DIM=2)/nSmooth
          ! XyzAux2_D = sum(MhData_VIB(&
          !     X_:Z_,iAux2-nSmooth+1:iAux2,iLine),DIM=2)/nSmooth
          ! State_VIB(D_, iVertex, iLine) = &
          !     sqrt(sum((XyzAux2_D - XyzAux1_D)**2)) / nSmooth
          ! end if
          ! distance from the beginning of the line
          if(iVertex == 1)then
             State_VIB(S_, iVertex, iLine) = 0.0
          else
             State_VIB(S_, iVertex, iLine) = &
                  State_VIB(S_, iVertex-1, iLine) + &
                  State_VIB(D_, iVertex-1, iLine)
          end if
          ! Heliocentric Distance
          State_VIB(R_, iVertex, iLine) = &
               norm2(MhData_VIB(X_:Z_, iVertex, iLine))
       end do
    end do

  end subroutine get_other_state_var
  !============================================================================
  subroutine get_shock_location
    use PT_ModSize, ONLY: nVertexMax
    ! find location of a shock wave on a given line (line)
    ! shock front is assumed to be location of max log(Rho/RhoOld)
    real :: dLogRho_I(1:nVertexMax)
    ! Threshold for shock tracing
    real, parameter :: DLogRhoThreshold = 0.01
    ! Do not search too close to the Sun
    real, parameter :: RShockMin = 1.20  ! *RSun
    integer         :: iShockMin
    ! Do not search too close to the heliosphere boundary
    integer:: iShockMax
    ! Misc
    integer:: iShockCandidate
    integer:: iLine, iVertex
    !--------------------------------------------------------------------------
    do iLine = 1, nLine
       if(.not.Used_B(iLine))then
          iShock_IB(Shock_,iLine) = NoShock_
          CYCLE
       end if
       ! shock front is assumed to be location of max log(Rho/RhoOld);
       do iVertex = 1, nVertex_B(  iLine)
          ! divergence of plasma velocity
          dLogRho_I(iVertex) = log(MhData_VIB(Rho_,iVertex,iLine)/&
               State_VIB(RhoOld_,iVertex,iLine))
       end do
       ! shock never moves back
       iShockMin = max(iShock_IB(ShockOld_, iLine), 1 + nWidth )
       iShockMax = nVertex_B(iLine) - nWidth - 1
       iShockCandidate = iShockMin - 1 + maxloc(&
            DLogRho_I(   iShockMin:iShockMax),1, MASK = &
            State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
            DLogRho_I(   iShockMin:iShockMax)       > DLogRhoThreshold)
       if(iShockCandidate >= iShockMin)&
            iShock_IB(Shock_, iLine) = iShockCandidate
    end do
  end subroutine get_shock_location
  !============================================================================
  subroutine search_line(iLine, Radius, iParticleOut, IsFound, Weight)

    ! performs search along given line
    ! for FIRST location ABOVE given heliocentric radius;
    ! if found, IsFound is set to .true. (.false. otherwise)
    ! and iParticleOut is set to index of particle just above Radius
    integer,      intent(in) :: iLine ! line/line index
    real,         intent(in) :: Radius ! heliocentric distance to find
    integer,      intent(out):: iParticleOut ! result: index
    logical,      intent(out):: IsFound! result: whether search was successful
    real,optional,intent(out):: Weight ! interpolation weight for output index

    integer:: iVertex ! loop variable
    !--------------------------------------------------------------------------
    ! check whether line reaches given radial distance
    if(State_VIB(R_, nVertex_B(iLine), iLine) < Radius)then
       ! mark failure to find location
       IsFound = .false.
       iParticleOut = -1
       RETURN
    end if

    ! line reaches given radial distance
    IsFound = .true.
    ! find index of first particle above Radius
    do iVertex = 1, nVertex_B(iLine)
       if(State_VIB(R_, iVertex, iLine) > Radius)then
          iParticleOut = iVertex
          EXIT
       end if
    end do

    ! get interpolation weight is necessary
    if(present(Weight))then
       if(iParticleOut > 1)then
          Weight = (Radius - State_VIB(R_,iParticleOut-1, iLine)) / (&
               State_VIB(R_,iParticleOut,  iLine) - &
               State_VIB(R_,iParticleOut-1,iLine))
       else
          Weight = 1.0
       end if
    end if

  end subroutine search_line
  !============================================================================
end module PT_ModGrid
!==============================================================================
