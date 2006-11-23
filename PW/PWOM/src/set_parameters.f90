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
  NameOutput  = 'PW/log.out'
  iUnitOutput = io_unit_new()
  open(UNIT=iUnitOutput, FILE=NameOutput)
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case('#STOP')
        if(IsStandAlone)then
           call read_var('Tmax',Tmax)
        else
           write(*,*)'PWOM WARNING: #STOP command is ignored in the framework'
        end if
     case('#SAVEPLOT')
        call read_var('DtSavePlot',DtOutput)
     case('#SCHEME')
        call read_var('TypeSolver',TypeSolver)
        call read_var('IsImplicit',IsImplicit)
        call read_var('DtVertical',DtPolarwind)
        
     case('#RESTART')
        call read_var('IsRestart',IsRestart)
     case('#MOTION')
        call read_var('IsMoveFluxTube',IsMoveFluxTube)
     case('#FAC')
        call read_var('IsUseJr',IsUseJr)
     case('#ROTATION')
        call read_var('IsCentrifugal',IsCentrifugal)
     case('#TIMESTEP')
        call read_var('DtMax',DtMax)
     case('#VERTICALGRID')
        call read_var('nPoints',nDim)
     case('#FIELDLINE')
        call read_var('nTotalLine',nTotalLine)
     endselect
  enddo

  Dt = DtMax
  write(iUnitOutput,*) tMax
  write(iUnitOutput,*) DToutput
  write(iUnitOutput,*) TypeSolver
  write(iUnitOutput,*) IsImplicit
  write(iUnitOutput,*) DTpolarwind
  write(iUnitOutput,*) TypeSolver

end subroutine PW_set_parameters
