subroutine PW_set_parameters(NameAction)

  use ModIoUnit, ONLY: UnitTmp_, io_unit_new
  use ModPwom
  use ModReadParam
  use ModCommonVariables, ONLY: F107,F107A,AP,UseStaticAtmosphere,DrBnd,&
                                UsePhotoElectronHeatFlux,UseAuroralHeatFlux, &
                                UseCuspHeatFlux
  use ModPwTime
  implicit none
  

  character (len=100) :: cLine
  character (len=100) :: cTempLine
  character (len=100), dimension(100) :: cTempLines
  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=*), parameter  :: NameSub = 'PW_set_parameters'
  real :: Vx, Bx, Bz, By, HPI
  
  !****************************************************************************
  ! This subroutine gets the inputs for PWOM
  !****************************************************************************

  real:: ddt1, xxx
  integer:: ns,iDate,iError

  !---------------------------------------------------------------------------
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case('#STOP')
        if(IsStandAlone)then
           call read_var('Tmax',Tmax)
           call read_var('MaxStep',MaxStep)
        else
           write(*,*)'PWOM WARNING: #STOP command is ignored in the framework'
        end if
     case('#STARTTIME')
        if(IsStandAlone)then
           !read in iYear,iMonth,iDay,iHour,iMinute,iSecond into iStartTime
           do iDate=1,6
              call read_var('iStartTime',iStartTime(iDate))
           enddo
           
        else
           write(*,*)'PWOM WARNING: #STARTTIME command is ignored in the framework'
        end if
     case('#STATICATMOSPHERE')
        call read_var('UseStaticAtmosphere', UseStaticAtmosphere)
     case('#MSISPARAM')
        call read_var('F107' ,F107)
        call read_var('F107A',F107A)
        call read_var('AP(1)',AP(1))
        call read_var('AP(2)',AP(2))
        call read_var('AP(3)',AP(3))
        call read_var('AP(4)',AP(4))
        call read_var('AP(5)',AP(5))
        call read_var('AP(6)',AP(6))
        call read_var('AP(7)',AP(7))
     case('#TIMEACCURATE')
        call read_var('DoTimeAccurate',DoTimeAccurate)
     case('#SAVEPLOT')
        call read_var('DtSavePlot',DtOutput)
        call read_var('DnSavePlot',DnOutput)
        call read_var('SaveFirst',DoSavePlot)
     case('#SAVEPLOTELECTRODYNAMICS')
        call read_var('DoPlotElectrodynamics',DoPlotElectrodynamics)
        call read_var('DtPlotElectrodynamics',DtPlotElectrodynamics)
     case('#SCHEME')
        call read_var('TypeSolver',TypeSolver)
        call read_var('DtVertical',DtVertical)
        call read_var('IsFullyImplicit'   ,IsFullyImplicit)
        if(IsFullyImplicit)then
           IsPointImplicit = .false.
           IsPointImplicitAll = .false.
        else
           call read_var('IsPointImplicit'   ,IsPointImplicit)
           call read_var('IsPointImplicitAll',IsPointImplicitAll)
        end if
     case('#VARIABLEDT')
        call read_var('IsVariableDt',IsVariableDt)
     case('#DIFFUSION')
        call read_var('TypeDiffusion',TypeDiffusion)
     case('#LIMITER')
        call read_var('LimiterBeta',BetaIn)
        Beta = BetaIn
     case('#RESTART')
        call read_var('IsRestart',IsRestart)
     case('#MOTION')
        call read_var('DoMoveLine',DoMoveLine)
     case('#FAC')
        call read_var('UseJr',UseJr)
     case('#AURORA')
        call read_var('UseAurora',UseAurora)
     case('#JOULEHEATING')
        call read_var('UseJouleHeating',UseJouleHeating)
     case('#ROTATION')
        call read_var('UseCentrifugal',UseCentrifugal)
     case('#TIMESTEP')
        call read_var('DtHorizontal',DtHorizontal)
        DtHorizontalOrig = DtHorizontal
     case('#VERTICALGRID')
        call read_var('nPoints',nAlt)
        call read_var('DeltaR',DrBnd)
     case('#FIELDLINE')
        call read_var('nTotalLine',nTotalLine)
     case('#LOG')
        call read_var('WriteLog',nLog) ! nLog=-1 write for all lines
                                       ! nLog= 0 write for no lines
                                       ! nLog= 1..nTotalLine, write one line
     case('#TEST')
        call read_var('StringTest',StringTest)
        call read_var('iProcTest', iProcTest)
        call read_var('iLinetest', iLineTest)
     case('#HEAT')
        call read_var('UseIonHeat',UseIonHeat)
        call read_var('UseEleHeat',UseEleHeat)
        if (UseEleHeat) then
           call read_var('UseExplicitHeat',UseExplicitHeat)
        else
           UseExplicitHeat = .false.
        endif
        
     case('#HEATFLUX')
        call read_var('UsePhotoElectronHeatFlux',UsePhotoElectronHeatFlux)
        call read_var('UseAuroralHeatFlux',UseAuroralHeatFlux)
        call read_var('UseCuspHeatFlux',UseCuspHeatFlux)
     case ("#MHD_INDICES")
        cTempLines(1) = NameCommand
        call read_var('UpstreamFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_MHDIMF_Indices(iError)
     case ("#SOLARWIND")
        call read_var('bx',bx)
        call read_var('by',by)
        call read_var('bz',bz)
        call read_var('vx',vx)
        call IO_set_imf_by_single(by)
        call IO_set_imf_bz_single(bz)
        call IO_set_sw_v_single(abs(vx))
        
     case ("#HPI")
        call read_var('HemisphericPower', HPI)
        call IO_set_hpi_single(HPI)
        
     case ("#NOAAHPI_INDICES")
        cTempLines(1) = "#NOAAHPI_INDICES"
        call read_var('NameHpiFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        call IO_set_inputs(cTempLines)
        call read_NOAAHPI_Indices(iError)
        if (iError /= 0) then
           write(*,*) "PW_ERROR: read hpi indices was NOT successful"
        endif

     case ("#NGDC_INDICES")
        cTempLines(1) = "#NGDC_INDICES"
        call read_var('NameNgdcFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_NGDC_Indices(iError)
        UseIndicies = .true.
    
        if (iError /= 0) then
           write(*,*) "PW_ERROR: read indices was NOT successful"
        endif
        
     endselect
  enddo
  
  
end subroutine PW_set_parameters

