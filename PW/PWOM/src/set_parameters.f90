subroutine PW_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModPwom
  implicit none

  character (len=*), intent(in) :: NameAction
  character (len=*), parameter :: NameSub = 'PW_set_parameters'

  !****************************************************************************
  ! This subroutine gets the inputs for PWOM
  !****************************************************************************

  real:: ddt1, xxx
  integer:: ns
  !---------------------------------------------------------------------------
  write(iUnitOut,*) NameSub,': called with action=',NameAction

  iUnitInput = UnitTmp_
  open(iUnitInput,       FILE=NameInput)

  iUnitOutput = io_unit_new()
  open(UNIT=iUnitOutput, FILE=NameOutput)

  READ(iUnitInput,*) TMAX
  WRITE(iUnitOutput,*) TMAX

  READ(iUnitInput,*) DToutput
  WRITE(iUnitOutput,*) DToutput

  READ(iUnitInput,*) TypeSolver
  WRITE(iUnitOutput,*) TypeSolver

  READ(iUnitInput,*) IsImplicit
  WRITE(iUnitOutput,*) IsImplicit

  read(iUnitInput,*) IsRestart
  if (IsRestart) then
     write(*,*) 'Is Restart', IsRestart
  endif
  read(iUnitInput,*) IsVariableDt
  if (IsVariableDt) then
     write(*,*) 'IsVariableDT', IsVariableDt
  endif

  READ(iUnitInput,*)   DTpolarwind
  WRITE(iUnitOutput,*) DTpolarwind

  READ(iUnitInput,*)   IsMoveFluxTube 
  READ(iUnitInput,*)   IsUseJr
  READ(iUnitInput,*)   IsCentrifugal

  CLOSE(UNIT=iUnitOutput)
  CLOSE(UNIT=iUnitInput)

end subroutine PW_set_parameters
