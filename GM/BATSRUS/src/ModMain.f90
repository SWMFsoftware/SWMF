!^CFG COPYRIGHT UM
Module ModMain
  !\
  ! Get values of nBLK, nI, nJ, nK, and gcn from ModSize
  !/
  use ModKind
  use ModSize

  implicit none

  save

  real :: CodeVersion=7.73, CodeVersionRead=-1.

  ! BATSRUS may run as GM or IH component. This will be set by CON.
  character (len=2)   :: NameThisComp=''

  ! In stand alone mode this variable is set to true
  logical             :: IsStandAlone=.false.

  ! We can use the old algorithms to calculate axes or the new one
  logical             :: UseNewAxes  =.true.

  ! We can use the old parameter format (parameters depend on time_accurate)
  ! or the new format compatible with SWMF
  logical             :: UseNewParam = .true.

  ! Number of bytes in the default real number (precision)
  ! nByteReal = 4 if range(1.0) is smaller than range(1.D0), otherwise 8
  ! This is standard F90 initialization expression but may give warnings:
  integer, parameter :: nByteReal = 4 + (1.00000000041 - 1.0)*10000000000.0

  !\
  ! Block parameters: nBlock is a maximum BLK number,  which is
  ! used by the given processor, nBlockMax is a maximum of nBlock
  ! over all the processors, nBlockALL is total number of blocks used.
  !/
  integer :: nBlockMax, nBlock, nBlockALL
  integer :: nBlockExplALL, nBlockImplALL        !^CFG IF IMPLICIT
  integer :: globalBLK

  ! Integer identifiers for the grid and the decomposition
  integer :: iNewGrid=0, iNewDecomposition=0
!^CFG IF NOT SIMPLE BEGIN
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
!^CFG END SIMPLE

  ! Named indexes for directions, number of dimensions
  integer, parameter:: x_=1, y_=2, z_=3

  integer, parameter:: R_=x_,Phi_=y_,Theta_=z_

  !\
  ! Face array parameters.
  !/
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
  ! Time stepping parameters and values.
  !/
  integer :: n_step, nOrder, nStage, iteration_number=0
  real :: dt, DtFixed, DtFixedDim, cfl, dt_BLK(nBLK)
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
  
  !\
  ! Dimensions of the buffer grid
  !/
  integer::  nPhiBuff=90,nThetaBuff=45
  real   ::RBuffMin=19.0,RBuffMax=21.0 

  logical :: UseIonosphere=.false.

  logical :: UseRaytrace=.false.                           !^CFG IF RAYTRACE

  logical :: UseIM = .false.                               !^CFG IF RCM
  real    :: tauCoupleIM=0.01                              !^CFG IF RCM

  real    :: ThetaTiltDeg = 0.0
  real    :: dt_UpdateB0  = 0.0001
  logical :: DoUpdateB0   = .true.
  logical :: DoSplitDb0Dt = .true.

  logical, dimension(nBLK) :: unusedBLK

  integer, dimension(nBLK) :: global_block_number

 ! Logical for bodies
  logical :: body1=.false.
  logical :: UseBody2=.false.         !^CFG IF SECONDBODY
    ! Integers for bodies
  integer, parameter :: body1_=-1
  integer, parameter :: body2_=-2    
  integer, parameter :: ExtraBc_=0   

  ! Inner and Outer boundary conditions
  character (len=20) :: TypeBc_I(body2_:top_)='float'

  ! Parameters for block location among eight subcubes.
  integer, parameter :: &
       LOC_TSE=1, &
       LOC_TSW=2, &
       LOC_BSW=3, &
       LOC_BSE=4, &
       LOC_BNE=5, &
       LOC_BNW=6, &
       LOC_TNW=7, &
       LOC_TNE=8

  ! How to deal with div B = 0
  logical :: UseDivbSource   =.true.
  logical :: UseDivbDiffusion=.false. !^CFG IF DIVBDIFFUSE
  logical :: UseProjection   =.false. !^CFG IF PROJECTION
  logical :: UseConstrainB   =.false. !^CFG IF CONSTRAINB
  logical :: UseB0Source

  !^CFG IF DIVBDIFFUSE BEGIN
  ! Default value for diffusion coefficient. Should be stable for CFL<1.
  ! Set 11/7/2000 - criteria defined/chosen by Gabor Toth
  real :: divb_diffcoeff=1.0/6.0
  !^CFG END DIVBDIFFUSE
  ! Choice of limiter
  character*6 :: limiter_type='minmod'
  real :: v_limiter_beta_param

  ! Prolongation order (1 or 2) and type ('central','lr','minmod','central2'..)
  integer :: prolong_order=1
  character*10 :: prolong_type='lr'

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

  ! Logical for rotating frame
  logical :: UseRotatingFrame

  ! Logical for corotation
  logical :: UseCorotation=.false.
  character(len=3) :: TypeCoordSystem='GSM'

  ! Logical for inertial frame
  logical :: UseInertial=.true. 

  ! String variable for debugging. A space separated list of words,
  ! typically names of subroutines. Listing "QUIET" suppresses progress reports

  character (len=79) :: test_string=''

  ! Location for test

  integer :: Itest=1, Jtest=1, Ktest=1, BLKtest=1, PROCtest=0, ITERtest=-1
  integer :: VARtest=1, DIMtest=1
  real    :: Xtest, Ytest, Ztest, Ttest
  real    :: Xtest_mod, Ytest_mod, Ztest_mod
  logical :: UseTestCell=.false., coord_test=.false.

  !\
  ! Debug parameters. 
  !/
  integer :: lVerbose = 1

  logical, parameter :: debug1=.false., debug2=.false.

  logical :: okdebug=.false., ShowGhostCells=.true.

  ! Time variables
  real :: Time_Simulation = 0.0

  ! This is the same default value as in SWMF
  integer, dimension(7) :: iStartTime_I=(/2000,3,21,10,45,0,0/)
  real(Real8_)          :: StartTime

  ! Optimization switches
  character*10 :: optimize_message_pass='dir'

  ! Solar Wind Input Parameters

  logical :: UseUpstreamInputFile
  real :: Satellite_Y_Pos, Satellite_Z_Pos

  ! Timing variables
  logical:: UseTiming = .true.
  integer:: dn_timing = -2

  ! Variables for controling the use of user files
  logical:: UseUserInnerBcs=.false.
  logical:: UseUserSource =.false.
  logical:: UseExtraBoundary = .false.
  logical:: DoFixExtraBoundary =.false.
  logical:: DoFixOuterBoundary =.false.  !^CFG IF FACEOUTERBC
  logical:: UseUserPerturbation =.false.
  logical:: UseUserOuterBcs=.false.
  logical:: UseUserICs=.false.
  logical:: UseUserSpecifyRefinement=.false.
  logical:: UseUserLogFiles=.false.
  logical:: UseUserWritePlot=.false.
  logical:: UseUserAMR=.false.
  logical:: UseUserEchoInput=.false.
  logical:: UseUserB0=.false.
  logical:: UseUserHeating=.false.
  logical:: UseUserSetPhysConst=.false.
  logical:: UseUserUpdateStates=.false.

  character (LEN=8) :: StringTimeH4M2S2

  logical :: DoSetLevels = .false.
  logical :: DoSendMHD   = .false.

  ! Shall we be strict about errors in the input parameter file
  logical :: UseStrict=.true.

  ! The following  variables are only used in stand alone mode
  real    :: t_max = -1.0, cputime_max = -1.0
  integer :: nIter = -1
  logical :: check_stopfile = .true., IsLastRead = .false.

end module ModMain
