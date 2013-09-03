subroutine CRCM_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModCrcmInitialize, ONLY: IsEmptyInitial,IsDataInitial, IsGmInitial
  use ModCrcmPlot,       ONLY: DtOutput, DoSavePlot, DoSaveFlux, DoSaveLog,&
                               UseSeparatePlotFiles, DtLogOut
  use ModFieldTrace,     ONLY: UseEllipse
  use ModCrcm,           ONLY: UseMcLimiter, BetaLimiter, time, Pmin
  use ModCrcmRestart,    ONLY: IsRestart
  use ModCrcmPlanet,     ONLY: nspec
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
        call read_var('DoSaveFlux',DoSaveFlux)
        ! If saving flux then decide if it should be just one file or many
        if (DoSaveFlux) then
           call read_var('UseSeparatePlotFiles',UseSeparatePlotFiles)
        endif
        
        DoSavePlot = .true.

     case('#SAVELOG')
        call read_var('DtLogOut',DtLogOut)
        DoSaveLog = .true.

     case('#INITIALF2')
        call read_var('IsEmptyInitial',IsEmptyInitial)
        call read_var('IsGmInitial',   IsGmInitial)
        call read_var('IsDataInitial', IsDataInitial)
        
        !IsDataInitial only works with EarthHO or EarthH configurations
        if (nspec > 3 .and. IsDataInitial) &
             call CON_STOP('IsDataInitial only works with EarthHO or EarthH')  
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

     ! minimum pressure in nPa passed to GM
     case('#MINIMUMPRESSURETOGM')   
        call read_var('MinimumPressureToGM', Pmin)
        
     case('#TIMESIMULATION')
        call read_var('TimeSimulation',time)
     
     end select
  enddo

end subroutine CRCM_set_parameters
