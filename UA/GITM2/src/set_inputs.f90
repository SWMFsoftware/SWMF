
subroutine set_inputs

  use ModInputs
  use ModSizeGitm
  use ModGITM, only: iProc
  use ModTime
  use ModPlanet
  use ModSatellites

  implicit none

  integer, external :: bad_outputtype
  integer, external :: jday

  integer, dimension(7) :: iEndTime

  logical :: IsDone
  integer :: iDebugProc=0, n
  character (len=iCharLen_) :: cLine
  integer :: iLine, iSpecies, iSat
  integer :: i, iError, iOutputTypes
  integer, dimension(7) :: iTimeEnd

  character (len=iCharLen_)                 :: cTempLine
  character (len=iCharLen_)                 :: sIonChemistry, sNeutralChemistry
  character (len=iCharLen_), dimension(100) :: cTempLines

  real :: Vx, Bx, Bz, By, Kp, HemisphericPower, tsim_temp
  real*8 :: DTime

  call report("set_inputs",1)

  iError = 0
  IsDone = .false.
  iLine  = 1

  do while (.not. IsDone)

     cLine = cInputText(iLine)

     if (cLine(1:1) == "#") then

        ! Remove anything after a space or TAB
        i=index(cLine,' '); if(i>0)cLine(i:len(cLine))=' '
        i=index(cLine,char(9)); if(i>0)cLine(i:len(cLine))=' '

        if (iDebugLevel > 3) write(*,*) "====> cLine : ",cLine(1:40)

        select case (cLine)

        case ("#TIMESTART","#STARTTIME")

           if (IsFramework) then

              if (iDebugLevel > 0) then
                 write(*,*) "----------------------------------"
                 write(*,*) "-   UAM trying to set STARTTIME  -"
                 write(*,*) "-          Ignoring              -"
                 write(*,*) "----------------------------------"
              endif

           else

              iStartTime = 0
              do i=1,6
                 call read_in_int(iStartTime(i), iError)
              enddo

              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #TIMESTART:'
                 write(*,*) '#TIMESTART'
                 write(*,*) 'iYear    (integer)'
                 write(*,*) 'iMonth   (integer)'
                 write(*,*) 'iDay     (integer)'
                 write(*,*) 'iHour    (integer)'
                 write(*,*) 'iMinute  (integer)'
                 write(*,*) 'iSecond  (integer)'
                 IsDone = .true.
              else

                 iTimeArray = iStartTime
                 call time_int_to_real(iStartTime, CurrentTime)
                 StartTime = CurrentTime
                 if (tSimulation > 0) then
                    CurrentTime = CurrentTime + tSimulation
                    call time_real_to_int(CurrentTime, iTimeArray)
                 endif

                 call fix_vernal_time

              endif

           endif

        case ("#TIMEEND","#ENDTIME")

           if (IsFramework) then
              if (iDebugLevel > 0) then
                 write(*,*) "--------------------------------"
                 write(*,*) "-   UAM trying to set ENDTIME  -"
                 write(*,*) "-          Ignoring            -"
                 write(*,*) "--------------------------------"
              endif
           else

              iTimeEnd = 0
              do i=1,6
                 call read_in_int(iTimeEnd(i), iError)
              enddo

              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #TIMEEND:'
                 write(*,*) '#TIMEEND'
                 write(*,*) 'iYear    (integer)'
                 write(*,*) 'iMonth   (integer)'
                 write(*,*) 'iDay     (integer)'
                 write(*,*) 'iHour    (integer)'
                 write(*,*) 'iMinute  (integer)'
                 write(*,*) 'iSecond  (integer)'
                 IsDone = .true.
              else
                 call time_int_to_real(iTimeEnd, EndTime)
              endif

           endif

        case ("#ISTEP")
           call read_in_int(iStep, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ISTEP:'
              write(*,*) '#ISTEP'
              write(*,*) 'iStep     (integer)'
              write(*,*) 'This is typically only specified in a'
              write(*,*) 'restart header.'
              IsDone = .true.
           endif

        case ("#CPUTIMEMAX")
           call read_in_real(CPUTIMEMAX, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #CPUTIMEMAX:'
              write(*,*) '#CPUTIMEMAX'
              write(*,*) 'CPUTimeMax    (real)'
           endif

        case ("#STATISTICALMODELSONLY")
           call read_in_logical(UseStatisticalModelsOnly, iError)
           if (UseStatisticalModelsOnly .and. iError == 0) &
                call read_in_real(DtStatisticalModels,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #STATISTICALMODELSONLY:'
              write(*,*) '#STATISTICALMODELSONLY'
              write(*,*) 'UseStatisticalModelsOnly    (logical)'
              write(*,*) 'DtStatisticalModels         (real)'
           endif

        case ("#TSIMULATION")

           if (IsFramework) then
              if (iDebugLevel > 0) then
                 write(*,*) "------------------------------------"
                 write(*,*) "-   UAM trying to set tsimulation  -"
                 write(*,*) "-          Ignoring                -"
                 write(*,*) "------------------------------------"
              endif
           else
              call read_in_real(tSim_temp, iError)
              tSimulation = tSim_temp
              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #TSIMULATION:'
                 write(*,*) '#TSIMULATION'
                 write(*,*) 'tsimulation    (real)'
                 write(*,*) 'This is typically only specified in a'
                 write(*,*) 'restart header.'
                 IsDone = .true.
              else
                 CurrentTime = CurrentTime + tsimulation
                 call time_real_to_int(CurrentTime, iTimeArray)
              endif
           endif

        case ("#F107")

           call read_in_real(f107, iError)
           call read_in_real(f107a, iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #F107:'
              write(*,*) '#F107'
              write(*,*) 'f107  (real)'
              write(*,*) 'f107a (real - 81 day average of f107)'
              f107  = 150.0
              f107a = 150.0
              IsDone = .true.
           endif

           call IO_set_f107_single(f107)
           call IO_set_f107a_single(f107a)

        case ("#INITIAL")

           call read_in_logical(UseMSIS, iError)
           call read_in_logical(UseIRI, iError)
           if (.not.UseMSIS) then
              call read_in_real(TempMin,    iError)
              call read_in_real(TempMax,    iError)
              call read_in_real(TempHeight, iError)
              call read_in_real(TempWidth,  iError)
              do i=1,nSpecies
                 call read_in_real(LogNS0(i), iError)
                 if (iError == 0) LogNS0(i) = alog(LogNS0(i))
              enddo
              LogRho0 = 0.0
              do iSpecies = 1, nSpecies
                 LogRho0 = LogRho0 + &
                      exp(LogNS0(iSpecies)) * Mass(iSpecies)
              enddo
              LogRho0 = alog(LogRho0)
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #INITIAL:'
              write(*,*) '#INITIAL'
              write(*,*) 'UseMSIS        (logical)'
              write(*,*) 'UseIRI         (logical)'
              write(*,*) 'If UseMSIS is .false. then :'
              write(*,*) 'TempMin        (real, bottom temperature)'
              write(*,*) 'TempMax        (real, top initial temperature)'
              write(*,*) 'TempHeight     (real, Height of the middle of temp gradient)'
              write(*,*) 'TempWidth      (real, Width of the temperature gradient)'
              do i=1,nSpecies
                 write(*,*) 'Bottom N Density (real), species',i
              enddo
           endif

        case ("#HPI")

           call read_in_real(HemisphericPower, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #HPI:'
              write(*,*) '#HPI'
              write(*,*) 'HemisphericPower  (real)'
              IsDone = .true.
           else
              call IO_set_hpi_single(HemisphericPower)
           endif

        case ("#KP")
           call read_in_real(kp, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #KP:'
              write(*,*) '#KP'
              write(*,*) 'kp  (real)'
              IsDone = .true.
           else
              call IO_set_kp_single(kp)
           endif

        case ("#CFL")
           call read_in_real(cfl, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #cfl:'
              write(*,*) '#CFL'
              write(*,*) 'cfl  (real)'
              IsDone = .true.
           endif

        case ("#SOLARWIND")
           call read_in_real(bx, iError)
           call read_in_real(by, iError)
           call read_in_real(bz, iError)
           call read_in_real(vx, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SOLARWIND:'
              write(*,*) '#SOLARWIND'
              write(*,*) 'bx  (real)'
              write(*,*) 'by  (real)'
              write(*,*) 'bz  (real)'
              write(*,*) 'vx  (real)'
              IsDone = .true.
           else
              call IO_set_imf_by_single(by)
              call IO_set_imf_bz_single(bz)
              call IO_set_sw_v_single(abs(vx))
           endif

        case ("#AMIEFILES")
           call read_in_string(cAMIEFileNorth, iError)
           call read_in_string(cAMIEFileSouth, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #AMIEFILES:'
              write(*,*) '#AMIEFILES'
              write(*,*) 'cAMIEFileNorth  (string)'
              write(*,*) 'cAMIEFileSouth  (string)'
              IsDone = .true.
           endif

        case ("#LIMITER")
           call read_in_string(TypeLimiter, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LIMITER:'
              write(*,*) '#LIMITER'
              write(*,*) 'TypeLimiter  (string)'
              IsDone = .true.
           endif

           call read_in_real(BetaLimiter, iError)
           if (iError /= 0) then
              BetaLimiter = 1.6
              if (TypeLimiter == "minmod") BetaLimiter = 1.0
              if (iProc == 0) then
                 write(*,*) "You can now set the limiter value yourself!"
                 write(*,*) '#LIMITER'
                 write(*,*) 'TypeLimiter  (string)'
                 write(*,*) 'BetaLimiter  (real between 1.0-minmod and 2.0-mc)'
              endif
              iError = 0
           else
              if (BetaLimiter < 1.0) BetaLimiter = 1.0
              if (BetaLimiter > 2.0) BetaLimiter = 2.0
              if (iDebugLevel > 2) &
                   write(*,*) "===>Beta Limiter set to ",BetaLimiter
           endif

        case ("#DEBUG")
           call read_in_int(iDebugLevel, iError)
           call read_in_int(iDebugProc, iError)
           call read_in_real(DtReport, iError)
           call read_in_logical(UseBarriers, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DEBUG:'
              write(*,*) '#DEBUG'
              write(*,*) 'iDebugLevel (integer)'
              write(*,*) 'iDebugProc  (integer)'
              write(*,*) 'DtReport    (real)'
              write(*,*) 'UseBarriers (logical)'
              IsDone = .true.
           endif

        case ("#THERMO")
           call read_in_logical(UseSolarHeating, iError)
           call read_in_logical(UseJouleHeating, iError)
           call read_in_logical(UseAuroralHeating, iError)
           call read_in_logical(UseNOCooling, iError)
           call read_in_logical(UseOCooling, iError)
           call read_in_logical(UseConduction, iError)
           call read_in_logical(UseTurbulentCond, iError)
           call read_in_logical(UseUpdatedTurbulentCond, iError)
           call read_in_real(EddyScaling, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMO:'
              write(*,*) '#THERMO'
              write(*,*) "UseSolarHeating   (logical)"
              write(*,*) "UseJouleHeating   (logical)"
              write(*,*) "UseAuroralHeating (logical)"
              write(*,*) "UseNOCooling      (logical)"
              write(*,*) "UseOCooling       (logical)"
              write(*,*) "UseConduction     (logical)"
              IsDone = .true.
           endif

        case ("#THERMALDIFFUSION")
           call read_in_real(KappaTemp0, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMALDIFFUSION:'
              write(*,*) '#THERMALDIFFUSION'
              write(*,*) "KappaTemp0    (thermal conductivity, real)"
           endif

        case ("#VERTICALSOURCES")
           call read_in_logical(UseEddyInSolver, iError)
           call read_in_logical(UseNeutralFrictionInSolver, iError)
           call read_in_real(MaximumVerticalVelocity, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #VERTICALSOURCES:'
              write(*,*) '#VERTICALSOURCES'
              write(*,*) "UseEddyInSolver              (logical)"
              write(*,*) "UseNeutralFrictionInSolver   (logical)"
              write(*,*) "MaximumVerticalVelocity      (real)"
           endif

        case ("#DIFFUSION")
           call read_in_logical(UseDiffusion, iError)
           if (UseDiffusion .and. iError == 0) then
              call read_in_real(EddyDiffusionCoef,iError)
              call read_in_real(EddyDiffusionPressure0,iError)
              call read_in_real(EddyDiffusionPressure1,iError)
           endif

           if (EddyDiffusionPressure0 < EddyDiffusionPressure1) then
              write(*,*) "If you use eddy diffusion, you must specify two pressure"
              write(*,*) "levels - under the first, the eddy diffusion is constant."
              write(*,*) "Between the first and the second, there is a linear drop-off."
              write(*,*) "Therefore The first pressure must be larger than the second!"
              iError = 1
           endif
             
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #DIFFUSION:'
              write(*,*) '#DIFFUSION'
              write(*,*) "UseDiffusion (logical)"
              write(*,*) "EddyDiffusionCoef (real)"
              write(*,*) "EddyDiffusionPressure0 (real)"
              write(*,*) "EddyDiffusionPressure1 (real)"
           endif

        case ("#FORCING")
           call read_in_logical(UsePressureGradient, iError)
           call read_in_logical(UseIonDrag, iError)
           call read_in_logical(UseNeutralFriction, iError)
           call read_in_logical(UseViscosity, iError)
           call read_in_logical(UseCoriolis, iError)
           call read_in_logical(UseGravity, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #THERMO:'
              write(*,*) '#FORCING'
              write(*,*) "UsePressureGradient (logical)"
              write(*,*) "UseIonDrag          (logical)"
              write(*,*) "UseNeutralFriction  (logical)"
              write(*,*) "UseViscosity        (logical)"
              write(*,*) "UseCoriolis         (logical)"
              write(*,*) "UseGravity          (logical)"
              IsDone = .true.
           endif

        case ("#IONFORCING")
           call read_in_logical(UseExB, iError)
           call read_in_logical(UseIonPressureGradient, iError)
           call read_in_logical(UseIonGravity, iError)
           call read_in_logical(UseNeutralDrag, iError)
           call read_in_logical(UseDynamo, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #IONFORCING:'
              write(*,*) '#IONFORCING'
              write(*,*) "UseExB                 (logical)"
              write(*,*) "UseIonPressureGradient (logical)"
              write(*,*) "UseIonGravity          (logical)"
              write(*,*) "UseNeutralDrag         (logical)"
              write(*,*) "UseDynamo              (logical)"
              IsDone = .true.
           endif

        case ("#CHEMISTRY")

           call read_in_logical(UseIonChemistry, iError)
           call read_in_logical(UseIonAdvection, iError)
           call read_in_logical(UseNeutralChemistry, iError)

!           call read_in_string(sNeutralChemistry, iError)
!           call read_in_string(sIonChemistry, iError)
!
!           iInputIonChemType = -1
!           iInputNeutralChemType = -1
!
!           do i = 1, nChemTypes_
!              if (sNeutralChemistry == sChemType(i)) iInputNeutralChemType = i
!              if (sIonChemistry == sChemType(i)) iInputIonChemType = i
!           enddo
!
!           if (iInputNeutralChemType < 1) then
!              write(*,*) "Error in #CHEMISTRY for Neutrals"
!              write(*,*) "Input type : ", sNeutralChemistry
!              write(*,*) "Acceptable types :"
!              do i = 1, nChemTypes_
!                 write(*,*) sChemType(i)
!              enddo
!              IsDone = .true.
!           endif
!
!           if (iInputIonChemType < 1) then
!              write(*,*) "Error in #CHEMISTRY for Ions"
!              write(*,*) "Input type : ", sIonChemistry
!              write(*,*) "Acceptable types :"
!              do i = 1, nChemTypes_
!                 write(*,*) sChemType(i)
!              enddo
!              IsDone = .true.
!           endif

        case ("#APEX")

           if (IsFramework .and. UseApex) then
              if (iDebugLevel >= 0) then
                 write(*,*) "---------------------------------------"
                 write(*,*) "-  While using Framework, you can not -"
                 write(*,*) "-   use APEX coordinates, sorry.      -"
                 write(*,*) "-          Ignoring                   -"
                 write(*,*) "---------------------------------------"
              endif
              UseApex = .false.
           else

              call read_in_logical(UseApex, iError)
              
              if (iError /= 0) then
                 write(*,*) 'Incorrect format for #APEX:'
                 write(*,*) '#APEX'
                 write(*,*) 'UseApex (logical)'
                 write(*,*) '        Sets whether to use a realistic magnetic'
                 write(*,*) '        field (T) or a dipole (F)'
                 IsDone = .true.
              endif

           endif


        case ("#GRID")
           if (nLats == 1 .and. nLons == 1) Is1D = .true.
           call read_in_int(nBlocksLon, iError)
           call read_in_int(nBlocksLat, iError)
           call read_in_real(LatStart, iError)
           call read_in_real(LatEnd, iError)
           call read_in_real(LonStart, iError)

           if (nLats > 1 .and. LatEnd-LatStart < 1) iError=1

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #GRID:'
              write(*,*) '#GRID'
              write(*,*) 'nBlocksLon   (integer)'
              write(*,*) 'nBlocksLat   (integer)'
              write(*,*) 'LatStart     (real)'
              write(*,*) 'LatEnd       (real)'
              write(*,*) 'LonStart     (real)'
              write(*,*) 'If LatStart and LatEnd are set to < -90 and'
              write(*,*) '> 90, respectively, then GITM does a whole'
              write(*,*) 'sphere.  If not, it models between the two.'
              write(*,*) 'If you want to do 1-D, set nLons=1, nLats=1 in'
              write(*,*) 'ModSizeGitm.f90, then recompile, then set LatStart'
              write(*,*) 'and LonStart to the point on the Globe you want'
              write(*,*) 'to model.'
              IsDone = .true.
           else

              if (LatStart <= -90.0 .and. LatEnd >= 90.0) then
                 IsFullSphere = .true.
              else
                 IsFullSphere = .false.
                 LatStart = LatStart * pi / 180.0
                 LatEnd   = LatEnd * pi / 180.0
                 LonStart = LonStart * pi / 180.0
              endif

           endif

        case ("#STRETCH")
           call read_in_real(ConcentrationLatitude, iError)
           call read_in_real(StretchingPercentage, iError)
           call read_in_real(StretchingFactor, iError)

           if (iError /= 0) then
              write(*,*) 'Incorrect format for #STRETCH:'
              write(*,*) '#STRETCH'
              write(*,*) 'ConcentrationLatitude (real, degrees)'
              write(*,*) 'StretchingPercentage  (real, 0-1)'
              write(*,*) 'StretchingFactor      (real)'
              write(*,*) 'Example (no stretching):'
              write(*,*) '#STRETCH'
              write(*,*) '65.0 ! location of minimum grid spacing'
              write(*,*) '0.0	 ! Amount of stretch 0 (none) to 1 (lots)'
              write(*,*) '1.0  '// &
                   '! More control of stretch ( > 1 stretch less '//&
                   '< 1 stretch more)'
              IsDone = .true.
           endif

        case("#TOPOGRAPHY")
           call read_in_logical(UseTopography, iError)
           if(iError /= 0) then
              write(*,*) 'Incorrect format for #TOPOGRAPHY:'
              write(*,*) '#TOPOGRAPHY'
              write(*,*) 'UseTopography (logical)'
              IsDone = .true.
           endif

        case ("#RESTART")
           call read_in_logical(DoRestart, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #RESTART:'
              write(*,*) '#RESTART'
              write(*,*) 'DoRestart (logical)'
              IsDone = .true.
           endif

        case ("#SAVEPLOTS", "#SAVEPLOT")
           call read_in_real(DtRestart, iError)
           call read_in_int(nOutputTypes, iError)
           if (nOutputTypes > nMaxOutputTypes) then
              write(*,*) "Sorry, nOutputTypes is Limited to ",nMaxOutputTypes
              iError = 1
           else
              do iOutputTypes=1,nOutputTypes
                 call read_in_string(OutputType(iOutputTypes), iError)
                 call read_in_real(DtPlot(iOutputTypes), iError)
              enddo
              iError = bad_outputtype()
              if (iError /= 0) then
                 write(*,*) "Error in #SAVEPLOTS"
                 write(*,*) "Bad output type : ",OutputType(iError)
              endif
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SAVEPLOT'
              write(*,*) '#SAVEPLOT'
              write(*,*) 'DtRestart (real, seconds)'
              write(*,*) 'nOutputTypes  (integer)'
              write(*,*) 'Outputtype (string, 3D, 2D, ION, NEUTRAL, ...)'
              write(*,*) 'DtPlot    (real, seconds)'
              IsDone = .true.
           endif

        case ("#SATELLITES")
           call read_in_int(nSats, iError)
           if (nSats > nMaxSats) then
              iError = 1
              write(*,*) "Too many satellites requested!!"
           else
              do iSat=1,nSats
                 call read_in_string(cSatFileName(iSat), iError)
                 call read_in_real(SatDtPlot(iSat), iError)
                 iSatCurrentIndex(iSat) = 0
              enddo
           endif
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #SATELLITES'
              write(*,*) '#SATELLITES'
              write(*,*) 'nSats     (integer - max = ',nMaxSats,')'
              write(*,*) 'SatFile1  (string)'
              write(*,*) 'DtPlot1   (real, seconds)'
              write(*,*) 'etc...'
              IsDone = .true.
           else
              call read_satellites(iError)
              if (iError /= 0) IsDone = .true.
           endif

        case ("#ELECTRODYNAMICS")
           call read_in_real(dTPotential, iError)
           call read_in_real(dTAurora, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ELECTRODYNAMICS'
              write(*,*) '#ELECTRODYNAMICS'
              write(*,*) 'DtPotential (real, seconds)'
              write(*,*) 'DtAurora    (real, seconds)'
              IsDone = .true.
           endif

        case ("#LTERadiation")
           call read_in_real(DtLTERadiation, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LTERadiation:'
              write(*,*) '#LTERadiation'
              write(*,*) 'DtLTERadiation (real)'
              IsDone = .true.
           endif

        case ("#IONPRECIPITATION")
           call read_in_logical(UseIonPrecipitation, iError)
           call read_in_string(IonIonizationFilename, iError)
           call read_in_string(IonHeatingRateFilename, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #IONPRECIPITATION'
              write(*,*) '#IONPRECIPITATION'
              write(*,*) 'UseIonPrecipitation     (logical)'
              write(*,*) 'IonIonizationFilename   (string)'
              write(*,*) 'IonHeatingRateFilename  (string)'
              IsDone = .true.
           endif

        case ("#LOGFILE")
           call read_in_real(dTLogFile, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #LOGFILE'
              write(*,*) '#LOGFILE'
              write(*,*) 'DtLogFile   (real, seconds)'
              IsDone = .true.
           endif

        case ("#ALTITUDE")
           call read_in_real(AltMin, iError)
           call read_in_real(AltMax, iError)
           call read_in_logical(UseStretchedAltitude, iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #ALTITUDE'
              write(*,*) '#ALTITUDE'
              write(*,*) 'AltMin                (real, km)'
              write(*,*) 'AltMax                (real, km)'
              write(*,*) 'UseStretchedAltitude  (logical)'
           else
              AltMin = AltMin * 1000.0
              AltMax = AltMax * 1000.0
           endif

        case ("#EUV_DATA")
           call read_in_logical(UseEUVData, iError)
           call read_in_string(cEUVFile, iError)
           
           if (UseEUVData) call Set_Euv(cEUVFile,iError)
           if (iError /= 0) then
              write(*,*) 'Incorrect format for #EUV_DATA'
              write(*,*) '#EUV_DATA'
              write(*,*) 'UseEUVData            (logical)'
              write(*,*) 'cEUVFile              (string)'
           endif

        case ("#GLOW")
           call read_in_logical(UseGlow, iError) 
           call read_in_real(dTGlow, iError)

        case ("#MHD_INDICES")

           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_MHDIMF_Indices(iError)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#NGDC_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_NGDC_Indices(iError)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#NOAAHPI_INDICES")
           cTempLines(1) = cLine
           call read_in_string(cTempLine, iError)
           cTempLines(2) = cTempLine
           cTempLines(3) = " "
           cTempLines(4) = "#END"

           call IO_set_inputs(cTempLines)
           call read_NOAAHPI_Indices(iError)

           if (iError /= 0) then 
              write(*,*) "read indices was NOT successful"
              IsDone = .true.
           else
              UseVariableInputs = .true.
           endif

        case ("#END")
           IsDone = .true.

        end select

        if (iError /= 0) IsDone = .true.

    endif

    iLine = iLine + 1

    if (iLine >= nInputLines) IsDone = .true.

  enddo

  if (iDebugProc >= 0 .and. iProc /= iDebugProc) then
     iDebugLevel = -1
  endif

  if (iError /= 0) then
     call stop_gitm("Must Stop!!")
  endif

!  KappaTemp0 = 3.6e-4

contains

  subroutine read_in_int(variable, iError)
    integer, intent(out) :: variable
    integer :: iError
    if (iError == 0) then 
       iline = iline + 1
       read(cInputText(iline),*,iostat=iError) variable
    endif
  end subroutine read_in_int

  subroutine read_in_logical(variable, iError)
    logical, intent(out) :: variable
    integer :: iError
    if (iError == 0) then
       iline = iline + 1
       read(cInputText(iline),*,iostat=iError) variable
    endif
  end subroutine read_in_logical

  subroutine read_in_string(variable, iError)
    character (len=iCharLen_), intent(out) :: variable
    integer :: iError
    if (iError == 0) then 
       iline = iline + 1
       variable = cInputText(iline)
       ! Remove anything after a space or TAB
       i=index(variable,' '); if(i>0)variable(i:len(cLine))=' '
       i=index(variable,char(9)); if(i>0)variable(i:len(cLine))=' '
    endif
  end subroutine read_in_string

  subroutine read_in_real(variable, iError)
    real, intent(out) :: variable
    integer :: iError
    if (iError == 0) then
       iline = iline + 1
       read(cInputText(iline),*, iostat=iError) variable
    endif
  end subroutine read_in_real

end subroutine set_inputs


subroutine fix_vernal_time

  use ModTime
  use ModPlanet

  implicit none

  real*8 :: DTime
  integer :: n
  integer, external :: jday

  DTime = CurrentTime - VernalTime

  if (DTime < 0.0) then
     n = 0
     do while (CurrentTime < VernalTime)
        VernalTime = VernalTime-int(DaysPerYear)*Rotation_Period
        n = n + 1
        if (floor(n*(DaysPerYear - int(DaysPerYear))) > 0) then
           VernalTime = VernalTime - &
                floor(n*(DaysPerYear - int(DaysPerYear))) * &
                Rotation_Period
           n = 0
        endif
     enddo
     DTime = CurrentTime - VernalTime
  else
     n = 0
     do while (DTime > SecondsPerYear)
        VernalTime = VernalTime+int(DaysPerYear)*Rotation_Period
        n = n + 1
        if (floor(n*(DaysPerYear - int(DaysPerYear))) > 0) then
           VernalTime = VernalTime + &
                floor(n*(DaysPerYear - int(DaysPerYear))) * &
                Rotation_Period
           n = 0
        endif
        DTime = CurrentTime - VernalTime
     enddo
  endif
  iDay  = DTime / Rotation_Period
  uTime = (DTime / Rotation_Period - iDay) * Rotation_Period
  iJulianDay = jday( &
       iTimeArray(1), iTimeArray(2), iTimeArray(3)) 

end subroutine fix_vernal_time
