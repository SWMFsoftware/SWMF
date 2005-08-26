!^CFG COPYRIGHT UM
Module ModMain
  !\
  ! Get values of nBLK, nI, nJ, nK, and gcn from ModSize
  !/
  use ModKind
  use ModSize

  implicit none

  save

  !\
  ! Version of BATSRUS
  !/
  real :: CodeVersion = 7.73, CodeVersionRead = -1.0

  !\
  ! Standalone and component information
  !/
  ! In stand alone mode this variable is set to true
  logical             :: IsStandAlone = .false.

  ! In the SWMF the BATSRUS may run as GM, SC or IH component
  character (len=2)   :: NameThisComp=''

  !\
  ! Precision of reals
  !/
  ! Number of bytes in the default real number (precision)
  ! nByteReal = 4 if range(1.0) is smaller than range(1.D0), otherwise 8
  ! This is a standard F90 initialization expression but may give warnings:
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0

  !\
  ! Problem definition.
  !/
  integer :: problem_type
  integer, parameter :: &
       problem_uniform    =1, &
       problem_shocktube  =2, &
       problem_heliosphere=3, &
       problem_comet      =5, &
       problem_rotation   =6, &
       problem_diffusion  =7, &
       problem_earth      =11,&
       problem_saturn     =12,&
       problem_jupiter    =13,&
       problem_venus      =14,&
       problem_mars       =15,&
       problem_cylinder   =21,&
       problem_sphere     =22,&
       problem_arcade     =25,&
       problem_cme        =26,&
       problem_globalhelio=28,&
       problem_dissipation=30
  character (len=80) :: StringProblemType_I(30)
  data StringProblemType_I / &
       'MHD Uniform Flow Problem', &
       'Shock Tube Initial Value Problem', &
       'Solar Wind and Inner Heliosphere Problem', &
       ' 4-Nameless Problem', &
       'Mass-Loaded Comet Problem', &
       'Rotation Test Case Problem', &
       'Magnetic Diffusion Test Case Problem', &
       ' 8-Nameless Problem', &
       ' 9-Nameless Problem', &
       '10-Nameless Problem', &
       'Earth Magnetosphere Problem', &
       'Saturn Magnetosphere Problem', &
       'Jupiter Magnetosphere Problem', &
       'Venus Ionosphere Problem', &
       'Mars Ionosphere Problem', &
       '16-Nameless Problem', &
       '17-Nameless Problem', &
       '18-Nameless Problem', &
       '19-Nameless Problem', &
       '20-Nameless Problem', &
       'Conducting Cylinder (2-D), MHD Shock Problem', &
       'Conducting Sphere (3-D), MHD Shock Problem', &
       '23-Nameless Problem', &
       '24-Nameless Problem', &
       'Arcade Eruption', &
       'CME Initiated by Magnetic Flux Rope', &
       '27-Nameless Problem', &
       'globalhelio', &
       '29-Nameless Problem', &
       'dissipation' /

  !\
  ! Named indexes for directions
  !/
  integer, parameter :: x_ = 1,  y_   = 2,  z_     = 3
  integer, parameter :: R_ = x_, Phi_ = y_, Theta_ = z_

  !\
  ! Time stepping parameters and values.
  !/
  integer :: n_step, nOrder, nStage, iteration_number=0
  real :: dt, DtFixed, DtFixedOrig, DtFixedDim, cfl, dt_BLK(nBLK)
  logical :: time_accurate=.true.,           &
       boris_correction, &                   !^CFG IF BORISCORR       
       UseBorisSimple,   &                   !^CFG IF SIMPLEBORIS
       time_loop=.false.

  logical :: UseDtFixed
  logical :: UsePointImplicit                !^CFG IF POINTIMPLICIT
  logical :: UseImplicit,  UsePartLocal      !^CFG IF IMPLICIT

  !\
  ! Model Coupling variables
  !/
  
  ! Send true MHD solution from IH to GM or values from an IMF file
  logical :: DoSendMHD   = .true.

  ! Dimensions of the buffer grid between SC and IH
  integer :: nPhiBuff = 90,   nThetaBuff = 45
  real    :: rBuffMin = 19.0, rBuffMax = 21.0 

  logical :: UseIonosphere=.false.

  logical :: UseRaytrace=.false.                           !^CFG IF RAYTRACE

  logical :: UseIM = .false.                               !^CFG IF RCM
  logical :: DoCoupleImPressure = .true.                   !^CFG IF RCM
  logical :: DoCoupleImDensity  = .false.                  !^CFG IF RCM
  real    :: TauCoupleIM = 20.0                            !^CFG IF RCM

  !\
  ! Parameters for the B0 field
  !/
  real    :: Dt_UpdateB0  = 0.0001
  logical :: DoUpdateB0   = .true.
  logical :: DoSplitDb0Dt = .true.

  !\
  ! Inner and outer boundaries
  !/
  ! Indexes for boundaries
  integer, parameter :: body1_   = -1
  integer, parameter :: body2_   = -2
  integer, parameter :: ExtraBc_ =  0

  ! Inner and outer boundary conditions
  character (len=20) :: TypeBc_I(body2_:top_) = 'float'

  ! Solar Wind Input Parameters
  logical :: UseUpstreamInputFile
  real :: Satellite_Y_Pos = 0.0, Satellite_Z_Pos = 0.0

  ! Logicals for bodies
  logical :: Body1    = .false.
  logical :: UseBody2 = .false.         !^CFG IF SECONDBODY

  !\
  ! Block AMR grid parameters
  !/

  ! Identifiers for the grid and decomposition changes
  integer :: iNewGrid=0, iNewDecomposition=0

  ! Logigal array for the blocks used (=.false.) on the given processor
  logical, dimension(nBLK) :: UnusedBLK

  ! nBlock is a maximum block index number used on the given processor
  ! nBlockMax is a maximum of nBlock over all the processors, 
  ! nBlockALL is total number of blocks used.
  
  integer :: nBlockMax, nBlock, nBlockALL
  integer :: nBlockExplALL, nBlockImplALL        !^CFG IF IMPLICIT

  ! The index of the block currently being dealt with.
  integer :: GlobalBLK

  ! Unique index of blocks for all the processors
  integer, dimension(nBLK) :: global_block_number

  ! Parameters for block location among eight subcubes.
  ! T=top, B=bottom, S=south, N=north, E=east, W=west
  integer, parameter :: &
       LOC_TSE=1, &
       LOC_TSW=2, &
       LOC_BSW=3, &
       LOC_BSE=4, &
       LOC_BNE=5, &
       LOC_BNW=6, &
       LOC_TNW=7, &
       LOC_TNE=8

  ! Variable for setting AMR levels
  logical :: DoSetLevels = .false.

  ! Index limits for the cell faces (needed for the constrained trasnsport)
  integer, parameter :: nIFace=nI+1 
  integer, parameter :: nJFace=nJ+1
  integer, parameter :: nKFace=nK+1  

  integer:: &
       jMinFaceX=1, jMaxFaceX=nJ,&
       kMinFaceX=1, kMaxFaceX=nK,&
       iMinFaceY=1, iMaxFaceY=nI,&
       kMinFaceY=1, kMaxFaceY=nK,&
       iMinFaceZ=1, iMaxFaceZ=nI,&
       jMinFaceZ=1, jMaxFaceZ=nJ

  !\
  ! How to deal with div B = 0
  !/
  logical :: UseDivbSource    = .true.
  logical :: UseDivbDiffusion = .false. !^CFG IF DIVBDIFFUSE
  logical :: UseProjection    = .false. !^CFG IF PROJECTION
  logical :: UseConstrainB    = .false. !^CFG IF CONSTRAINB
  logical :: UseB0Source

  !^CFG IF DIVBDIFFUSE BEGIN
  ! Default value for unsplit div B diffusion coefficient.
  real :: Divb_DiffCoeff=1.0/6.0
  !^CFG END DIVBDIFFUSE

  !\
  ! More numerical scheme parameters
  !/
  ! Choice of limiter
  character (len=6) :: Limiter_Type = 'minmod'
  real              :: BetaLimiter

  ! Prolongation order (1 or 2) and type ('central','lr','minmod','central2'..)
  integer           :: prolong_order = 1
  character(len=10) :: prolong_type  = 'lr'

  ! Message passing mode ('dir', 'opt', 'all', 'allopt' ...)
  character(len=10) :: optimize_message_pass = 'dir'

  !\
  ! Source terms
  !/

  ! Logical for mass loading
  logical :: UseMassLoading=.false.
  logical :: AccelerateMassLoading=.false.

  !^CFG IF DISSFLUX BEGIN
  ! Logical for adding heat conduction
  logical:: UseHeatFlux=.false.
  logical:: UseSpitzerForm=.false.

  ! Logical for adding resistivity
  logical:: UseResistFlux=.false.
  logical:: UseAnomResist=.false.
  !^CFG END DISSFLUX

  ! Logical and type for gravity
  logical :: UseGravity
  integer :: GravityDir

  !\
  ! Rotation and coordinate system
  !/

  ! Logical for rotating inner boundary
  logical          :: UseRotatingBc = .false.

  ! Coordinate system
  character(len=3) :: TypeCoordSystem = 'GSM'

  ! Rotating frame or (at least approximately) inertial frame
  logical :: UseRotatingFrame = .false. 

  ! Transformation from HGC to HGI systems can be done in a simulation
  logical :: DoTransformToHgi = .false.

  !\
  ! Variables for debugging. 
  !/

  ! Shall we be strict about errors in the input parameter file
  logical :: UseStrict=.true.

  ! Verbosity level (0, 1, 10, 100)
  integer :: lVerbose = 1

  ! A space separated list of words, typically names of subroutines.
  character (len=79) :: test_string=''

  ! Location for test
  integer :: Itest=1, Jtest=1, Ktest=1, BLKtest=1, PROCtest=0, ITERtest=-1
  integer :: VARtest=1, DIMtest=1
  real    :: Xtest, Ytest, Ztest, Ttest
  real    :: Xtest_mod, Ytest_mod, Ztest_mod
  logical :: UseTestCell=.false., coord_test=.false.

  ! Debug logicals
  logical, parameter :: debug1=.false., debug2=.false.
  logical :: okdebug=.false., ShowGhostCells=.true.

  !\
  ! Time and timing variables
  !/
  real :: Time_Simulation = 0.0

  ! This is the same default value as in the SWMF
  integer, dimension(7) :: iStartTime_I = (/2000,3,21,10,45,0,0/)
  real(Real8_)          :: StartTime

  ! Timing variables
  logical:: UseTiming = .true.
  integer:: dn_timing = -2

  ! Simulation time in hour,minute,second format
  character (LEN=8) :: StringTimeH4M2S2

  !\
  ! Stopping conditions. These variables are only used in stand alone mode.
  !/
  real    :: t_Max = -1.0, cputime_max = -1.0
  integer :: nIter = -1
  logical :: Check_Stopfile = .true.
  logical :: IsLastRead = .false.

  !\
  ! Controling the use of the features implemented in user files
  !/
  logical:: UseUserInnerBcs          = .false.
  logical:: UseUserSource            = .false.
  logical:: UseExtraBoundary         = .false.
  logical:: DoFixExtraBoundary       = .false.
  logical:: DoFixOuterBoundary       = .false.  !^CFG IF FACEOUTERBC
  logical:: UseUserPerturbation      = .false.
  logical:: UseUserOuterBcs          = .false.
  logical:: UseUserICs               = .false.
  logical:: UseUserSpecifyRefinement = .false.
  logical:: UseUserLogFiles          = .false.
  logical:: UseUserWritePlot         = .false.
  logical:: UseUserAMR               = .false.
  logical:: UseUserEchoInput         = .false.
  logical:: UseUserB0                = .false.
  logical:: UseUserHeating           = .false.
  logical:: UseUserSetPhysConst      = .false.
  logical:: UseUserUpdateStates      = .false.

end module ModMain
