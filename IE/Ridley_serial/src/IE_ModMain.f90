!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module IE_ModMain

  use ModConst

  implicit none
  save

  integer, parameter :: iono_init=1, iono_fac=2, iono_read=3, iono_save=4, &
       iono_save_restart=5, iono_solve=6

  !\
  ! Ionospheric Model Parameters
  !/
  logical :: UseFakeRegion2=.false.
  integer :: conductance_model
  logical :: UseFullCurrent
  real    :: f107_flux, PolarCapPedConductance, StarLightPedConductance
  real    :: Hall_to_Ped_Ratio

  !\
  ! Conductance model 4 factors
  !/
  real :: OvalWidthFactor=1., OvalStrengthFactor=1., CondFactor=1.7


  !\
  ! Time variables obtained from CON_physics
  !/
  logical               :: time_accurate
  integer, dimension(7) :: Time_Array
  real                  :: Time_Simulation

  !\
  ! Counter for number of solves (like nStep in GM)
  !/
  integer :: nSolve = 0

  !\
  ! Krylov solver parameters
  !/
  character(len=10):: NameSolver='bicgstab' ! Name of krylov solver
  logical :: UsePreconditioner = .true.     ! Use preconditioner
  logical :: UseInitialGuess   = .true.     ! Use previous solution as initial guess
  real    :: Tolerance         = 1.e-2      ! Solution accuracy: abs norm of residual
  integer :: MaxIteration      = 200        ! Maximum number of Krylov iterations

  !\
  ! Logical which tells if there is any new information to use
  !/
  logical :: IsNewInput = .false.

  !\
  ! Character string selecting the potential sent to the IM module
  ! Possible values are 'north', 'south', 'average', 'cpcpmin'
  !/
  character (len=7) :: TypeImCouple = 'north   '
  !\
  ! Fraction of IM current to consider
  !/
  real :: FractionImJr = 1.0

  !\
  ! Parameter for coupling the UA current and 
  ! latitude boundary for currents from GM and for calculating Phi
  !/
  logical :: DoCoupleUaCurrent = .false.
  real    :: LatBoundaryGm, LatBoundary = 90.0 * cDegToRad

  !\
  ! Dipole parameters obtained from CON_physics
  !/
  real :: ThetaTilt, SinThetaTilt, CosThetaTilt

  !\
  ! Save logfile?
  !/
  logical :: DoSaveLogfile = .true.

  integer :: iDebugLevel = 0

end module IE_ModMain
