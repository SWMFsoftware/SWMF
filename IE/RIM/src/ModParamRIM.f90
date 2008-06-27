
module ModParamRIM

  use ModNumConst, only : cDegToRad

  implicit none

  logical :: DoSolve = .false.
  logical :: DoFold  = .false.
  real    :: HighLatBoundary = 90.0 * cDegToRad
  real    :: LowLatBoundary  = 55.0 * cDegToRad

  logical :: UseGMCurrents = .true.
  logical :: UseIMCurrents = .false.
  logical :: UseUACurrents = .false.
  logical :: UseUAConductances = .false.

  character (len=100) :: NameEFieldModel="weimer96"
  character (len=100) :: NameAuroralModel="ihp"
  character (len=100) :: NameSolarModel="mb"

  logical :: UseStaticIMF=.true.

  integer :: iConductanceModel=5
  real    :: f107flux=150.
  real    :: StarlightPedConductance=1.
  real    :: PolarCapPedConductance=0.25

  character (len=7) :: TypeImCouple = 'north'

  !\
  ! Krylov solver (GMRES) parameters
  !/
  logical :: UsePreconditioner = .true.! Use preconditioner
  logical :: UseInitialGuess = .true.  ! Use previous solution as initial guess
  real    :: Tolerance = 1.e-2        ! Solution accuracy: 2nd norm of residual
  integer :: MaxIteration = 500       ! Maximum number of Krylov iterations

  integer :: iDebugLevel=0

  logical :: DoSaveLogfile=.true.

  integer, parameter :: SolveAcrossEquator_  = 1
  integer, parameter :: SolveWithOutEquator_ = 2
  integer, parameter :: SolveWithFold_       = 3

  integer :: SolveType = 1

  logical :: DoTouchNorthPole = .false.
  logical :: DoTouchSouthPole = .false.

  logical :: DoPrecond = .false.

end module ModParamRIM

