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
       if(nWave > 1)then
          TradMin = 1.0; EradMin = cRadiation*TradMin**4
          ! Reset the minimum photon energy to be 0.1 eV
          FreqMinSi = 0.1/cHPlanckEV
          ! Reset the maximum photon energy to be 10 keV
          FreqMaxSi = 10000.0/cHPlanckEV
          call set_multigroup(nWave, FreqMinSi, FreqMaxSi)
       end if
    end select

  end subroutine user_init_session

  !============================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKEVToK
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK
    use ModPhysics,    ONLY: cRadiationNo, g, No2Si_V, Si2No_V, &
         UnitTemperature_
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, ExtraEint_, p_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    real :: Rho, Temperature, Pressure, Erad
    real :: TeFinal

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(TypeProblem)
    case('lightfront')
       Rho = 1.0
       Temperature = 1.0e-16
       Erad = 1.0e-12
    case('infinitemedium')
       ! initial zero radiation, at time=infinity Erad(final)=a*te(final)**4
       ! so that inv_gm1*rho*te(final)+a*te(final)**4 = inv_gm1*rho*te(initial)
       TeFinal = cKEVToK*Si2No_V(UnitTemperature_)
       Rho = cRadiationNo*TeFinal**3
       Temperature = (Rho*TeFinal + (g - 1)*cRadiationNo*TeFinal**4)/Rho
       Erad = 0.0
    end select

    Pressure = Rho*Temperature

    do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
       State_VGB(Rho_,i,j,k,iBlock) = Rho
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = Erad
       State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0
       State_VGB(p_,i,j,k,iBlock) = Pressure
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModImplicit,   ONLY: StateSemi_VGB, iTrImplFirst, iTrImplLast
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
          State_VGB(:,-1,:,:,iBlock)  = State_VGB(:,1,:,:,iBlock)
          ! bin 1 starts on the left
          State_VGB(WaveFirst_,-1:0,:,:,iBlock) = 1.0
       case('usersemi')
          ! bin 1 starts on the left
          if(nWave == 1)then
             ! set Erad
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,0,j,k,iBlock) = 1.0
             end do; end do
          else
             ! set group temperature
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,0,j,k,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
                StateSemi_VGB(iTrImplLast,0,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplLast,1,j,k,iBlock)
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
           if(nWave > 1) State_VGB(WaveLast_,nI+1:nI+2,:,:,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,nI+1,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,nI,j,k,iBlock)
             end do; end do
          else
             ! bin 2 starts on the right
             do k = 1, nK; do j = 1, nJ
                StateSemi_VGB(iTrImplFirst,nI+1,j,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,nI,j,k,iBlock)
                StateSemi_VGB(iTrImplLast,nI+1,j,k,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
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
          else
             ! bin 1 starts on the left
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,0,k,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
                StateSemi_VGB(iTrImplLast,i,0,k,iBlock) = &
                     StateSemi_VGB(iTrImplLast,i,1,k,iBlock)
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
          if(nWave > 1) State_VGB(WaveLast_,:,nJ+1:nJ+2,:,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,nJ+1,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,nJ,k,iBlock)
             end do; end do
          else
             ! bin 2 starts on the right
             do k = 1, nK; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,nJ+1,k,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,nJ,k,iBlock)
                StateSemi_VGB(iTrImplLast,i,nJ+1,k,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
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
          else
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,0,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
                StateSemi_VGB(iTrImplLast,i,j,0,iBlock) = &
                     StateSemi_VGB(iTrImplLast,i,j,1,iBlock)
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
          if(nWave > 1) State_VGB(WaveLast_,:,:,nK+1:nK+2,iBlock) = 1.0
       case('usersemi')
          if(nWave == 1)then
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,nK+1,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,j,nK,iBlock)
             end do; end do
          else
             ! bin 2 starts on the right
             do j = 1, nJ; do i = 1, nI
                StateSemi_VGB(iTrImplFirst,i,j,nK+1,iBlock) = &
                     StateSemi_VGB(iTrImplFirst,i,j,nK,iBlock)
                StateSemi_VGB(iTrImplLast,i,j,nK+1,iBlock) = &
                     sqrt(sqrt(1.0/cRadiationNo))
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

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKtoKev
    use ModGeometry,   ONLY: dy_BLK, dz_BLK
    use ModMain,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, cRadiationNo
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
    character(len=10) :: NameWave, NameFormat

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
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
    case('dy')
       PlotVar_G(:,:,:) = dy_BLK(iBlock)
    case('dz')
       PlotVar_G(:,:,:) = dz_BLK(iBlock)
    case default
       IsFound = .false.
    end select

    if(nWave < 10)then
       NameFormat = "(a,i1)"
    else
       NameFormat = "(a,i2.2)"
    end if

    do iWave = 1, nWave
       write(NameWave, NameFormat) 'erad', iWave
       if(NameVar == NameWave)then
          iVar = WaveFirst_ + iWave -1
          do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
             PlotVar_G(i,j,k) = State_VGB(iVar,i,j,k,iBlock)
          end do; end do; end do
          IsFound = .true.
          EXIT
       end if
    end do

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, &
       PlanckOut_W, CgTeOut_W, CgTgOut_W, TgOut_W)

    ! The State_V vector is in normalized units

    use CRASH_ModMultiGroup, ONLY: get_energy_g_from_temperature, &
         get_temperature_from_energy_g
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
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    !  Derivative of group Planckian by electron temperature
    !  Group specific heat of the radiation
    !  Group radiation temperature
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]
    real, optional, intent(out) :: CgTeOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: CgTgOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: TgOut_W(nWave)          ! [K]

    integer :: iVar, iWave
    real :: Rho, Pressure, Te, Tg_W(nWave)
    real :: RhoSi, pSi, TeSi
    real :: PlanckSi, CgTeSi, EgSi, TgSi, CgTgSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalIn))then
       pSi = EinternalIn*gm1
       Pressure = pSi*Si2No_V(UnitP_)
       Te = Pressure/Rho
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

    if(present(EinternalOut)) EinternalOut = pSi*inv_gm1
    if(present(TeOut)) TeOut = TeSi
    if(present(PressureOut)) PressureOut = pSi

    if(present(CvOut)) CvOut = inv_gm1*Rho &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(HeatCondOut)) HeatCondOut = 0.0
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0

    select case(TypeProblem)
    case('lightfront')
       if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = 0.0

       if(present(OpacityRosselandOut_W)) &
            OpacityRosselandOut_W = 1.0e-16/No2Si_V(UnitX_)

       Tg_W = sqrt(sqrt(State_V(WaveFirst_:WaveLast_)/cRadiationNo))

       if(present(PlanckOut_W)) PlanckOut_W = cRadiationNo*Te**4 &
            *No2Si_V(UnitEnergyDens_)

       if(present(CgTeOut_W)) CgTeOut_W = &
            cRadiationNo*(Te+Tg_W)*(Te**2+Tg_W**2) &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

       if(present(TgOut_W)) TgOut_W = Tg_W*No2Si_V(UnitTemperature_)

       if(present(CgTgOut_W)) CgTgOut_W = 4.0*cRadiationNo*Tg_W**3 &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    case('infinitemedium')
       if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = &
            1.0/No2Si_V(UnitX_)

       if(present(OpacityRosselandOut_W)) &
            OpacityRosselandOut_W = 1.0e20/No2Si_V(UnitX_)

       if(nWave == 1)then
          if(present(PlanckOut_W)) PlanckOut_W = 0.0
          if(present(CgTeOut_W)) CgTeOut_W = 0.0
          if(present(TgOut_W)) TgOut_W = 0.0
          if(present(CgTgOut_W)) CgTgOut_W = 0.0
       else

          if(present(PlanckOut_W) .or. present(CgTeOut_W))then
             do iWave = 1, nWave
                call get_energy_g_from_temperature( &
                     iWave, TeSi, EgSI=PlanckSi, CgSI=CgTeSi)

                if(present(PlanckOut_W)) PlanckOut_W(iWave) = PlanckSi
                if(present(CgTeOut_W)) CgTeOut_W(iWave) = CgTeSi
             end do
          end if

          if(present(TgOut_W) .or. present(CgTgOut_W))then
             do iWave = 1, nWave
                iVar = WaveFirst_ + iWave - 1
                EgSi = State_V(iVar)*No2Si_V(UnitEnergyDens_)
                call get_temperature_from_energy_g(iWave, EgSi, &
                     TgSIOut=TgSi, CgSIOut=CgTgSi)

                if(present(TgOut_W)) TgOut_W(iWave) = TgSi
                if(present(CgTgOut_W)) CgTgOut_W(iWave) = CgTgSi
             end do
          end if

       end if

    end select

  end subroutine user_material_properties

end module ModUser
