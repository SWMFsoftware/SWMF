!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_set_plot_var,               &
       IMPLEMENTED7 => user_material_properties,        &
       IMPLEMENTED8 => user_init_session

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'multi-group radiation diffusion tests'

  character(len=20) :: TypeProblem

  real :: EphotonMin = 0.1 ! eV
  real :: EphotonMax = 2e4 ! eV
  real :: TeFinalSi, TeFinal, TeInit

contains

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    integer, intent(in):: iStage, iBlock

    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    ! No call to update_states_MHD to nullify the effect of the hydro solver
    ! call update_states_MHD(iStage,iBlock)

  end subroutine user_update_states

  !============================================================================

  subroutine user_read_inputs

    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModMain,      ONLY: lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_line, read_command, read_var

    character (len=100) :: NameCommand
    !--------------------------------------------------------------------------
    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case("#TYPEPROBLEM")
          call read_var('TypeProblem',TypeProblem)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

  end subroutine user_read_inputs

  !============================================================================

  subroutine user_init_session

    use CRASH_ModMultiGroup, ONLY: TradMin, EradMin, set_multigroup
    use ModConst,            ONLY: cHPlanckEV, cRadiation
    use ModVarIndexes,       ONLY: nWave

    real ::  FreqMinSi, FreqMaxSi

    character (len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    select case(TypeProblem)
    case('lightfront')
    case('infinitemedium')
    case('planckian')
       if(nWave > 1)then
          TradMin = 1e-4; EradMin = cRadiation*TradMin**4
          ! Reset the minimum photon energy to be 0.1 eV
          FreqMinSi = EphotonMin/cHPlanckEV
          ! Reset the maximum photon energy to be 20 keV
          FreqMaxSi = EphotonMax/cHPlanckEV
          call set_multigroup(nWave, FreqMinSi, FreqMaxSi)
       end if
    end select

  end subroutine user_init_session

  !============================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKEVToK
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK
    use ModPhysics,    ONLY: cRadiationNo, inv_gm1, g, No2Si_V, Si2No_V, &
         UnitTemperature_
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, ExtraEint_, p_, &
         nWave, WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    real :: Rho, Temperature, Pressure, Erad, Trad, ExtraEint

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(TypeProblem)
    case('lightfront')
       Rho = 1.0
       Temperature = 1.0e-16
       Erad = 1.0e-12
       Pressure = Rho*Temperature
       ExtraEint = 0.0
    case('infinitemedium')
       ! initial zero radiation, at time=infinity Erad(final)=a*te(final)**4,
       ! where a is the radiation constant. Set Cv = 4*a*Te**3,
       ! so that a*te(final)**4 + Erad(final)= a*te(initial)**4
       ! or Te(initial)**4 = 2*Te(final)**4
       TeFinalSi = cKEVToK
       TeFinal = TeFinalSi*Si2No_V(UnitTemperature_)
       Rho = 1.0 ! fake, it is not used
       TeInit = sqrt(sqrt(1.0 + nWave))*TeFinal
       Temperature = TeInit
       Erad = 0.0
       Pressure = Rho*Temperature
       ExtraEint = cRadiationNo*Temperature**4 - inv_gm1*Pressure
    case('planckian')
       ! initial zero radiation and infinite heat capacity and given
       ! electron temperature Te.
       ! at time infinity: Erad(final)= a*Te**4
       TeFinalSi = cKEVToK
       TeFinal = TeFinalSi*Si2No_V(UnitTemperature_)
       Rho = 1.0e18 ! this ensure large heat capacity
       Temperature = TeFinal
       Erad = 0.0
       Pressure = Rho*Temperature
       ExtraEint = 0.0
    end select

    do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
       State_VGB(Rho_,i,j,k,iBlock) = Rho
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = Erad
       State_VGB(ExtraEint_,i,j,k,iBlock) = ExtraEint
       State_VGB(p_,i,j,k,iBlock) = Pressure
    end do; end do; end do


  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: StateSemi_VGB, iTrImplFirst, iTrImplLast, &
         UseSplitSemiImplicit, iVarSemi
    use ModMain,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: cRadiationNo
    use ModVarIndexes, ONLY: nWave, WaveFirst_, WaveLast_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    IsFound = .true.

    select case(iSide)
    case(1)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,0,:,:,iBlock)  = State_VGB(:,1,:,:,iBlock)
          State_VGB(:,-1,:,:,iBlock) = State_VGB(:,1,:,:,iBlock)
          ! bin 1 starts on the left
          State_VGB(WaveFirst_,-1:0,:,:,iBlock) = 1.0
       case('usersemi')
          ! bin 1 starts on the left
          if(nWave == 1)then
             ! set Erad
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,0,j,k,iBlock) = 1.0
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,0,1:nJ,1:nK,iBlock) = 1.0
             else
                StateSemi_VGB(1,0,1:nJ,1:nK,iBlock) = &
                     StateSemi_VGB(1,1,1:nJ,1:nK,iBlock)
             end if
          else
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,0,j,k,iBlock) = 1.0
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,0,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst+1:iTrImplLast,1,j,k,iBlock)
             end do; end do
          end if
       end select
    case(2)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          ! bin 2 starts on the right
          if(nWave > 1) &
               State_VGB(WaveFirst_+1:WaveLast_,nI+1:nI+2,:,:,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,nI+1,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,nI,j,k,iBlock)
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,nI+1,1:nJ,1:nK,iBlock) = &
                     StateSemi_VGB(1,nI,1:nJ,1:nK,iBlock)
             else
                StateSemi_VGB(1,nI+1,1:nJ,1:nK,iBlock) = 1.0
             end if
          else
             ! bin 2 starts on the right
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,nI+1,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,nI,j,k,iBlock)
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,nI+1,j,k,iBlock) = 1.0
             end do; end do
          end if
       end select
    case(3)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,0,:,iBlock)  = State_VGB(:,:,1,:,iBlock)
          State_VGB(:,:,-1,:,iBlock)  = State_VGB(:,:,1,:,iBlock)
          ! bin 1 starts on the left
          State_VGB(WaveFirst_,:,-1:0,:,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,0,k,iBlock) = 1.0
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,1:nI,0,1:nK,iBlock) = 1.0
             else
                StateSemi_VGB(1,1:nI,0,1:nK,iBlock) = &
                     StateSemi_VGB(1,1:nI,1,1:nK,iBlock)
             end if
          else
             ! bin 1 starts on the left
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,0,k,iBlock) = 1.0
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,0,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,1,k,iBlock)
             end do; end do
          end if
       end select
    case(4)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,nJ+1,:,iBlock) = State_VGB(:,:,nJ,:,iBlock)
          State_VGB(:,:,nJ+2,:,iBlock) = State_VGB(:,:,nJ,:,iBlock)
          ! bin 2 starts on the right
          if(nWave > 1) &
               State_VGB(WaveFirst_+1:WaveLast_,:,nJ+1:nJ+2,:,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,nJ+1,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,nJ,k,iBlock)
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,1:nI,nJ+1,1:nK,iBlock) = &
                     StateSemi_VGB(1,1:nI,nJ,1:nK,iBlock)
             else
                StateSemi_VGB(1,1:nI,nJ+1,1:nK,iBlock) = 1.0
             end if
          else
             ! bin 2 starts on the right
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,nJ+1,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,nJ,k,iBlock)
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,nJ+1,k,iBlock) = 1.0
             end do; end do
          end if
       end select
    case(5)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,:,0,iBlock)  = State_VGB(:,:,:,1,iBlock)
          State_VGB(:,:,:,-1,iBlock)  = State_VGB(:,:,:,1,iBlock)
          ! bin 1 starts on the left
          State_VGB(WaveFirst_,:,:,-1:0,iBlock) = 1.0
       case('usersemi')
          ! bin 1 starts on the left
          if(nWave == 1)then
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,0,iBlock) = 1.0
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,1:nI,1:nJ,0,iBlock) = 1.0
             else
                StateSemi_VGB(1,1:nI,1:nJ,0,iBlock) = &
                     StateSemi_VGB(1,1:nI,1:nJ,1,iBlock)
             end if
          else
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,0,iBlock) = 1.0
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,j,0,iBlock) = &
                     StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,j,1,iBlock)
             end do; end do
          end if
       end select
    case(6)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,:,nK+1,iBlock) = State_VGB(:,:,:,nK,iBlock)
          State_VGB(:,:,:,nK+2,iBlock) = State_VGB(:,:,:,nK,iBlock)
          ! bin 2 starts on the right
          if(nWave > 1) &
               State_VGB(WaveFirst_+1:WaveLast_,:,:,nK+1:nK+2,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,nK+1,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,j,nK,iBlock)
             end do; end do
          elseif(UseSplitSemiImplicit)then
             if(iVarSemi == iTrImplFirst)then
                StateSemi_VGB(1,1:nI,1:nJ,nK+1,iBlock) = &
                     StateSemi_VGB(1,1:nI,1:nJ,nK,iBlock)
             else
                StateSemi_VGB(1,1:nI,1:nJ,nK+1,iBlock) = 1.0
             end if
          else
             ! bin 2 starts on the right
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,nK+1,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,j,nK,iBlock)
                StateSemi_VGB(iTrImplFirst+1:iTrImplLast,i,j,nK+1,iBlock) = 1.0
             end do; end do
          end if
       end select
    case default
       write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
       call stop_mpi(NameSub)
    end select

  end subroutine user_set_outerbcs

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use CRASH_ModMultiGroup, ONLY: get_planck_g_from_temperature
    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKtoKev, cKEV
    use ModGeometry,   ONLY: dy_BLK, dz_BLK
    use ModIo,         ONLY: NamePlotDir
    use ModIoUnit,     ONLY: io_unit_new
    use ModMain,       ONLY: nI, nJ, nK, Time_Simulation, n_step
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, UnitEnergyDens_, &
         cRadiationNo, inv_gm1, Si2No_V
    use ModVarIndexes, ONLY: Rho_, p_, nWave, WaveFirst_, WaveLast_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k, iVar, iWave
    real :: DelLogEphoton, Coord_W(nWave), PlanckSi, EinternalSi, TeSi
    integer :: iUnit, iError
    character(len=10) :: NameWave, NameFormat, TypeStatus
    character(len=100):: NameFile = 'planckian.outs'
    character(len=10) :: TypePosition = 'rewind'

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('planck')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          PlotVar_G(i,j,k) = cRadiationNo &
               *(State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**4
       end do; end do; end do
    case('bfinal')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          PlotVar_G(i,j,k) = cRadiationNo*TeFinal**4
       end do; end do; end do
    case('tkev')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          PlotVar_G(i,j,k) = &
               State_VGB(p_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock) &
               *No2Si_V(UnitTemperature_)*cKToKev
       end do; end do; end do
    case('trkev')
       NameIdlUnit = 'KeV'
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          PlotVar_G(i,j,k) = sqrt(sqrt(sum( &
               State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))/cRadiationNo)) &
               *No2Si_V(UnitTemperature_)*cKToKev
       end do; end do; end do
    case('etotal')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               EinternalOut = EinternalSi)
          PlotVar_G(i,j,k) = sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))&
               + EinternalSi*Si2No_V(UnitEnergyDens_)
       end do; end do; end do
    case default
       IsFound = .false.
    end select

    if(TypeProblem == 'planckian' .and. NameVar == 'planckian')then
       IsFound = .true.

       DelLogEphoton = (log(EphotonMax) - log(EphotonMin))/nWave
       do iWave = 1, nWave
          Coord_W(iWave) = exp(log(EphotonMin) + (iWave - 0.5)*DelLogEphoton)
       end do

       iUnit = io_unit_new()

       TypeStatus = 'replace'
       if(TypePosition == 'append')TypeStatus = 'unknown'
       open(iUnit, file=trim(NamePlotDir)//NameFile, &
            position = TypePosition, status=TypeStatus, iostat=iError)
       if(iError /= 0)call CON_stop(NameSub // &
            ' could not open ascii file=' // trim(NamePlotDir)//NameFile)

       write(iUnit, "(a)")             'group energy'
       write(iUnit, "(i7,es13.5,3i3)") n_step, Time_Simulation, 1, 1, 2
       write(iUnit, "(3i8)")           nWave
       write(iUnit, "(100es13.5)")     0.0
       write(iUnit, "(a)")             'Ephoton[eV] Egroup[eV] Bgroup[eV]'

       do iWave = 1, nWave
          call get_planck_g_from_temperature(iWave, TeFinalSi, PlanckSi)
          iVar = WaveFirst_ + iWave - 1
          write(iUnit, "(100es18.10)") Coord_W(iWave), &
               State_VGB(iVar,1,1,1,1)*No2Si_V(UnitEnergyDens_)*cKEV, &
               PlanckSi*cKEV
       end do

       close(iUnit)

       TypePosition = 'append'
    end if

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

    ! The State_V vector is in normalized units

    use CRASH_ModMultiGroup, ONLY: get_planck_g_from_temperature
    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_, cRadiationNo, Clight
    use ModVarIndexes, ONLY: nVar, Rho_, p_, nWave, WaveFirst_, WaveLast_

    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in):: i, j, k, iBlock, iDir  ! cell/face index
    real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
    real, optional, intent(in)  :: TeIn                    ! [K]
    real, optional, intent(out) :: NatomicOut              ! [1/m^3]
    real, optional, intent(out) :: EinternalOut            ! [J/m^3]
    real, optional, intent(out) :: TeOut                   ! [K]
    real, optional, intent(out) :: PressureOut             ! [Pa]
    real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: GammaOut                ! dimensionless
    real, optional, intent(out) :: HeatCondOut             ! [J/(m*K*s)]
    real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]

    integer :: iVar, iWave
    real :: Rho, Pressure, Te
    real :: RhoSi, pSi, TeSi
    real :: PlanckSi, EgSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalIn))then
       select case(TypeProblem)
       case('infinitemedium')
          Te = sqrt(sqrt(EinternalIn*Si2No_V(UnitEnergyDens_)/cRadiationNo))
          Pressure = Rho*Te
          pSi = Pressure*No2Si_V(UnitP_)
       case default
          pSi = EinternalIn*gm1
          Pressure = pSi*Si2No_V(UnitP_)
          Te = Pressure/Rho
       end select
       TeSi = Te*No2Si_V(UnitTemperature_)
    elseif(present(TeIn))then
       TeSi = TeIn
       Te = TeSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Te
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Te = Pressure/Rho
       TeSi = Te*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalOut))then
       select case(TypeProblem)
       case('infinitemedium')
          EinternalOut = cRadiationNo*Te**4*No2Si_V(UnitEnergyDens_)
       case default
          EinternalOut = pSi*inv_gm1
       end select
    end if
    if(present(TeOut)) TeOut = TeSi
    if(present(PressureOut)) PressureOut = pSi

    if(present(CvOut))then
       select case(TypeProblem)
       case('lightfront','planckian')
          CvOut = inv_gm1*Rho &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)
       case('infinitemedium')
          CvOut = 4.0*cRadiationNo*Te**3 &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)
       end select
    end if

    if(present(HeatCondOut)) HeatCondOut = 0.0
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0

    select case(TypeProblem)
    case('lightfront')
       if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = 0.0

       if(present(OpacityRosselandOut_W)) &
            OpacityRosselandOut_W = 1.0e-16/No2Si_V(UnitX_)

       if(present(PlanckOut_W)) PlanckOut_W = cRadiationNo*Te**4 &
            *No2Si_V(UnitEnergyDens_)

    case('infinitemedium')
       if(present(OpacityPlanckOut_W))OpacityPlanckOut_W = 2.0/No2Si_V(UnitX_)

       if(present(OpacityRosselandOut_W)) &
            OpacityRosselandOut_W = 1.0e20/No2Si_V(UnitX_)

       if(present(PlanckOut_W)) PlanckOut_W = cRadiationNo*Te**4 &
            *No2Si_V(UnitEnergyDens_)

    case('planckian')
       if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = &
            2.0/No2Si_V(UnitX_)

       if(present(OpacityRosselandOut_W)) &
            OpacityRosselandOut_W = 1.0e20/No2Si_V(UnitX_)

       if(nWave == 1)then
          if(present(PlanckOut_W)) PlanckOut_W = cRadiationNo*Te**4 &
               *No2Si_V(UnitEnergyDens_)
       else

          if(present(PlanckOut_W))then
             do iWave = 1, nWave
                call get_planck_g_from_temperature(iWave, TeSi, PlanckSi)
                PlanckOut_W(iWave) = PlanckSi
             end do
          end if
       end if

    end select

  end subroutine user_material_properties

end module ModUser
