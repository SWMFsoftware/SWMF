
module ModParamRIM
  implicit none

  logical :: DoSolve=.true.
  logical :: DoFold=.false.
  real :: HighLatBoundary=90., LowLatBoundary=55.

  logical :: UseGMCurrents=.true.
  logical :: UseIMCurrents=.false.
  logical :: UseUACurrents=.false.

  character (len=100) :: NameEFieldModel="weimer96"
  character (len=100) :: NameAuroralModel="ihp"
  character (len=100) :: NameSolarModel="mb"

  logical :: UseStaticIMF=.true.
  real :: IMFBx, IMFBy, IMFBz, SWVx

  integer :: iConductanceModel=5
  real :: f107flux=150., StarlightPedConductance=1., PolarCapPedConductance=0.25

  character (len=7) :: TypeImCouple = 'north'

  !\
  ! Krylov solver (GMRES) parameters
  !/
  logical :: UsePreconditioner = .true. ! Use preconditioner
  logical :: UseInitialGuess = .true.   ! Use previous solution as initial guess
  real    :: Tolerance = 1.e-2          ! Solution accuracy: 2nd norm of residual
  integer :: MaxIteration = 200         ! Maximum number of Krylov iterations

  integer :: iDebugLevel=0

  logical :: DoSaveLogfile=.true.



end module ModParamRIM

