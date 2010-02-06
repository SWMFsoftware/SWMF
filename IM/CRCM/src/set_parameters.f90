subroutine CRCM_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModCrcmInitialize, ONLY: IsRestart
  use ModCrcmPlot,       ONLY: DtOutput, DoSavePlot
  implicit none

  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=*), parameter  :: NameSub = 'CRCM_set_parameters'

  !\
  ! Description:
  ! This subroutine gets the inputs for CRCM
  !/

  !---------------------------------------------------------------------------
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     
     case('#SAVEPLOT')
        call read_var('DtSavePlot',DtOutput)
        DoSavePlot = .true.
        
     case('#RESTART')
        call read_var('IsRestart',IsRestart) !T:Continuous run
                                             !F:Initial run
     end select
  enddo

end subroutine CRCM_set_parameters
