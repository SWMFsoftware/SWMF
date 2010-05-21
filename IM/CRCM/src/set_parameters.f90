subroutine CRCM_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModCrcmInitialize, ONLY: IsEmptyInitial
  use ModCrcmPlot,       ONLY: DtOutput, DoSavePlot
  use ModFieldTrace,     ONLY: UseEllipse
  use ModCrcm,           ONLY: UseMcLimiter, BetaLimiter
  use ModCrcmRestart,    ONLY: IsRestart
  implicit none

  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=7)             :: TypeBoundary
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

     case('#INITIALF2')
        call read_var('IsEmptyInitial',IsEmptyInitial)

     case('#TYPEBOUNDARY')
        call read_var('TypeBoundary',TypeBoundary)
        if(TypeBoundary == 'Ellipse') then
           UseEllipse = .true.
        else
           UseEllipse = .false.
        endif
        
     case('#RESTART')
        call read_var('IsRestart',IsRestart) !T:Continuous run
 
                                            !F:Initial run
     case('#LIMITER')
        call read_var('UseMcLimiter', UseMcLimiter)
        if(UseMcLimiter) call read_var('BetaLimiter', BetaLimiter)
     
     end select
  enddo

end subroutine CRCM_set_parameters
