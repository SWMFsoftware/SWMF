subroutine CIMI_set_parameters(NameAction)

  use ModIoUnit,         ONLY: UnitTmp_, io_unit_new
  use ModReadParam
  use ModCimiInitialize, ONLY: IsEmptyInitial,IsDataInitial,IsRBSPData, &
       IsGmInitial
  use ModCimiPlot,       ONLY: DtOutput, DoSavePlot, DoSaveFlux, DoSaveDrifts,&
       DoSaveLog, UseSeparatePlotFiles, DtLogOut,DoSavePSD
  use ModCimiTrace,     ONLY: UseEllipse, UseSmooth, UseCorotation, &
       UsePotential, SmoothWindow, imod, iLatTest, iLonTest, rb
  use ModCimi,           ONLY: UseMcLimiter, BetaLimiter, time, Pmin,&
       IsStandAlone, UseStrongDiff, UseDecay, DecayTimescale, dt, dtmax
  use ModCimiRestart,    ONLY: IsRestart,DtSaveRestart
  use ModCimiPlanet,     ONLY: nspec
  use ModImTime,         ONLY: iStartTime_I, TimeMax
  use ModCimiBoundary,   ONLY: UseBoundaryEbihara,UseYoungEtAl
  use ModIeCimi,         ONLY: UseWeimer
  use ModPrerunField,    ONLY: DoWritePrerun, UsePrerun, DtRead
  use ModGmCIMI,         ONLY: UseGm
  use ModWaveDiff,       ONLY: UseWaveDiffusion,UseHiss,UseChorus,UseChorusUB, &
       DiffStartT,HissWavesD, ChorusWavesD,ChorusUpperBandD, &
       testDiff_aa, testDiff_EE, testDiff_aE, &
       NameAeFile, read_ae_wdc_kyoto
  use DensityTemp, 	 ONLY: densityP
  use ModImSat,          ONLY: DtSatOut, DoWritePrerunSat, UsePrerunSat, &
       DtReadSat, DoWriteSats, ReadRestartSat
  use ModCimiGrid
  
  implicit none

  logical 			:: DoEcho=.false.
  character (len=100)           :: NameCommand
  character (len=*), intent(in) :: NameAction
  character (len=7)             :: TypeBoundary
  character (len=3)             :: NameModel
  character (len=*), parameter  :: NameSub = 'CIMI_set_parameters'

  character (len=100) :: cTempLine
  character (len=100), dimension(100) :: cTempLines
  
  integer :: iError, iDate
  real :: DensitySW, VelSW, BxSW, BySW,BzSW
  !\
  ! Description:
  ! This subroutine gets the inputs for CIMI
  !/

  !---------------------------------------------------------------------------
  
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE

     select case(NameCommand)

     case("#ECHO")
!        call check_stand_alone
        call read_var('DoEcho', DoEcho)
        if(iProc==0)call read_echo_set(DoEcho)

     case('#STOP')
        if(IsStandAlone)then
           call read_var('TimeMax',TimeMax)
        else
           write(*,*)'IM WARNING: #STOP command is ignored in the framework'
        end if

     case('#STARTTIME')
        if(IsStandAlone)then
           !read in iYear,iMonth,iDay,iHour,iMinute,iSecond into iStartTime
           do iDate=1,6
              call read_var('iStartTime',iStartTime_I(iDate))
           enddo
           
        else
           write(*,*)'IM WARNING:#STARTTIME command is ignored in the framework'
        end if

     case ("#NGDC_INDICES")
        cTempLines(1) = NameCommand
        call read_var('NameNgdcFile', cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_NGDC_Indices(iError)
        
        if (iError /= 0) then 
           write(*,*) "read indices was NOT successful (NOAA file)"
        endif

     case ("#MHD_INDICES")
        cTempLines(1) = NameCommand
        call read_var('UpstreamFile',cTempLine)
        cTempLines(2) = cTempLine
        cTempLines(3) = " "
        cTempLines(4) = "#END"
        
        call IO_set_inputs(cTempLines)
        call read_MHDIMF_Indices(iError)

     case ("#SOLARWIND")
        call read_var('n',DensitySW)
        call read_var('vx',VelSW)
        call read_var('bx',BxSW)
        call read_var('by',BySW)
        call read_var('bz',BzSW)
        
        call IO_set_imf_bx_single(BxSW)
        call IO_set_imf_by_single(BySW)
        call IO_set_imf_bz_single(BzSW)
        call IO_set_sw_v_single(abs(VelSW))
        call IO_set_sw_n_single(abs(DensitySW))

     case('#SMOOTH')
        !do you want to smooth the solar wind for Tsy model
        call read_var('UseSmooth',UseSmooth)
        call read_var('SmoothWindow',SmoothWindow)

     case('#BMODEL')
        call read_var('NameModel',NameModel)!t96,t04,MHD,Dip
        !call read_var('UseFixedB',UseFixedB)!T=fixed B config or 
                                            !F=changing B config
        if (NameModel == 'Dip') then
           iMod=0
           ! Checks the next two lines to see if cor and/or pot
           ! variables should be set to 0.
           call read_var('UseCorotation',UseCorotation)
           call read_var('UsePotential',UsePotential)
        elseif (NameModel == 't96') then
           iMod=1
        elseif(NameModel == 't04') then
           iMod=2
        elseif(NameModel == 'MHD')then
           iMod=3
           UseGm=.true.
        else
           call con_stop('Error: Model not found') 
        endif

     case('#IEMODEL')
        call read_var('UseWeimer',UseWeimer)

     case('#SAVEPLOT')
        call read_var('DtSavePlot',DtOutput)
        call read_var('DoSaveFlux',DoSaveFlux)
        call read_var('DoSaveDrifts',DoSaveDrifts)
        call read_var('DoSavePSD',DoSavePSD)
        ! If saving flux then decide if it should be just one file or many
        if (DoSaveFlux .or. DoSaveDrifts .or. DoSavePSD) then
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
        if (IsDataInitial) &
             call read_var('IsRBSPData', IsRBSPData)
        
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
        
     case('#PLASMASHEET')
        call read_var('UseYoungEtAl',UseYoungEtAl)
        if(IsStandAlone)then
           ! When Standalone specify boundary type
           call read_var('UseBoundaryEbihara',UseBoundaryEbihara)
        endif
        
     case('#RESTART')
        call read_var('IsRestart',IsRestart) !T:Continuous run, F: Initial Run
        if(IsRestart) call read_var('ReadRestartSat',ReadRestartSat)
        
     case('#SAVERESTART')
        ! when in standalone mode read restart save frequency
        if(IsStandAlone) then
           call read_var('DtSaveRestart',DtSaveRestart)
        else
           call CON_STOP('ERROR: Restart save frequency only set in standalone')
        end if

     case('#LIMITER')
        call read_var('UseMcLimiter', UseMcLimiter)
        if(UseMcLimiter) call read_var('BetaLimiter', BetaLimiter)

     ! minimum pressure in nPa passed to GM
     case('#MINIMUMPRESSURETOGM')   
        call read_var('MinimumPressureToGM', Pmin)
        
     case('#TIMESIMULATION')
        call read_var('TimeSimulation',time)

     case('#PRERUNSAT')
        call read_var('DoWritePrerunSat',DoWritePrerunSat)
        if(.not.DoWritePrerunSat) call read_var('UsePrerunSat',UsePrerunSat)
        if(UsePrerunSat) then
           call read_var('DtReadSat',   DtReadSat)
           DtSatout = DtReadSat
           DoWriteSats = .true.
        end if
        
     case('#PRERUNFIELD')
        call read_var('DoWritePrerun',DoWritePrerun)
        if(.not.DoWritePrerun) call read_var('UsePrerun',UsePrerun)
        if(UsePrerun)          call read_var('DtRead',   DtRead)
        if(UsePrerun .or. DoWritePrerun) UseGm=.true.

     case('#DECAY')
        call read_var('UseDecay',UseDecay)
        if ( UseDecay ) &
             call read_var('DecayTimescale in seconds', DecayTimescale)
        
     case('#STRONGDIFFUSION')
        call read_var('UseStrongDiff',UseStrongDiff)
        
     case('#WAVEDIFFUSION')
        call read_var('UseWaveDiffusion',UseWaveDiffusion)
        if(UseWaveDiffusion) then
           call read_var('DiffStartT',  DiffStartT)
           call read_var('UseHiss',  UseHiss)
           call read_var('UseChorus',  UseChorus)
           call read_var('UseChorusUB',  UseChorusUB)
           call read_var('HissWavesD', HissWavesD)
           call read_var('ChorusWavesD', ChorusWavesD)
           call read_var('ChorusUpperBandD', ChorusUpperBandD)
           call read_var('NameAeFile',NameAeFile)
           call read_ae_wdc_kyoto(iError)

           if (iError /= 0) then
              write(*,*) "read AE index was not successful "//&
                   "(WDC Kyoto file; IAGA2002 Format)"
           end if

        end if

     case('#ENERGYGRID')
        call read_var('MinIonEnergy (in keV), MinElectronEnergy x10',&
             MinIonEnergy)
        call read_var('MaxIonEnergy (in keV), MaxElectronEnergy x10',&
             MaxIonEnergy)
        
     case('#IMTIMESTEP')
        call read_var('IMDeltaT [s]',dt)
        call read_var('IMDeltaTMax [s]',dtmax)

     case('#TESTFIELDLINE')
        call read_var('iLatTest',iLatTest)
        call read_var('iLonTest',iLonTest)
        
     case('#PLASMAPAUSEDENSITY')
        call read_var('DensityP [m^3]',densityP)
        
     case('#SETRB')
        call read_var('rb [R_E]', rb)
        
     end select
     
  enddo

end subroutine CIMI_set_parameters
