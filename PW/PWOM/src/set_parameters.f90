subroutine PW_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModPwom
  use ModReadParam
  implicit none

  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=*), parameter  :: NameSub = 'PW_set_parameters'

  !****************************************************************************
  ! This subroutine gets the inputs for PWOM
  !****************************************************************************

  real:: ddt1, xxx
  integer:: ns
  !---------------------------------------------------------------------------
  NameOutput  = 'log.out'
  iUnitOutput = io_unit_new()
  open(UNIT=iUnitOutput, FILE=NameOutput)
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case('#POLARWIND')
        call read_var('Tmax',Tmax)
        call read_var('DToutput',DToutput)
        call read_var('TypeSolver',TypeSolver)
        call read_var('IsImplicit',IsImplicit)
        call read_var('IsRestart',IsRestart)
        call read_var('IsVariableDt',IsVariableDt)
        call read_var('DtPolarWind',DTpolarwind)
        call read_var('IsMoveFluxTube',IsMoveFluxTube)
        call read_var('IsUseJr',IsUseJr)
        call read_var('IsCentrifugal',IsCentrifugal)
     endselect
  enddo

  write(iUnitOutput,*) TMAX
  write(iUnitOutput,*) DToutput
  write(iUnitOutput,*) TypeSolver
  write(iUnitOutput,*) IsImplicit
  write(iUnitOutput,*) DTpolarwind
  write(iUnitOutput,*) TypeSolver

end subroutine PW_set_parameters
