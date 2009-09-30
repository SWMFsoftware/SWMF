!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_set_plot_var,               &
       IMPLEMENTED7 => user_material_properties

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

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKevToK
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK
    use ModPhysics,    ONLY: inv_gm1, cRadiationNo, Si2No_V, UnitTemperature_
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, ExtraEint_, p_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    real :: Temperature, Pressure

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    Temperature = 1.0e-16
    Pressure = Temperature
    do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
       State_VGB(Rho_,i,j,k,iBlock) = 1.0
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 1.0d-12
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
    use ModVarIndexes, ONLY: WaveFirst_, WaveLast_

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
          ! bin 1 is starts on the left
          State_VGB(WaveFirst_,-1:0,:,:,iBlock) = 1.0
       case('usersemi')
          ! bin 1 is starts on the left
          do k = 1, nK; do j = 1, nJ
             StateSemi_VGB(iTrImplFirst,0,j,k,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
             StateSemi_VGB(iTrImplLast,0,j,k,iBlock) = &
                  StateSemi_VGB(iTrImplLast,1,j,k,iBlock)
          end do; end do
       end select
    case(2)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          ! bin 2 starts on the right
          State_VGB(WaveLast_,nI+1:nI+2,:,:,iBlock) = 1.0
       case('usersemi')
          ! bin 2 starts on the right
          do k = 1, nK; do j = 1, nJ
             StateSemi_VGB(iTrImplFirst,nI+1,j,k,iBlock) = &
                  StateSemi_VGB(iTrImplFirst,nI,j,k,iBlock)
             StateSemi_VGB(iTrImplLast,nI+1,j,k,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
          end do; end do
       end select
    case(3)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,0,:,iBlock)  = State_VGB(:,:,1,:,iBlock)
          State_VGB(:,:,-1,:,iBlock)  = State_VGB(:,:,1,:,iBlock)
          ! bin 1 is starts on the left
          State_VGB(WaveFirst_,:,-1:0,:,iBlock) = 1.0
       case('usersemi')
          ! bin 1 is starts on the left
          do k = 1, nK; do i = 1, nI
             StateSemi_VGB(iTrImplFirst,i,0,k,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
             StateSemi_VGB(iTrImplLast,i,0,k,iBlock) = &
                  StateSemi_VGB(iTrImplLast,i,1,k,iBlock)
          end do; end do
       end select
    case(4)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,nJ+1,:,iBlock) = State_VGB(:,:,nJ,:,iBlock)
          State_VGB(:,:,nJ+2,:,iBlock) = State_VGB(:,:,nJ,:,iBlock)
          ! bin 2 starts on the right
          State_VGB(WaveLast_,:,nJ+1:nJ+2,:,iBlock) = 1.0
       case('usersemi')
          ! bin 2 starts on the right
          do k = 1, nK; do i = 1, nI
             StateSemi_VGB(iTrImplFirst,i,nJ+1,k,iBlock) = &
                  StateSemi_VGB(iTrImplFirst,i,nJ,k,iBlock)
             StateSemi_VGB(iTrImplLast,i,nJ+1,k,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
          end do; end do
       end select
    case(5)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,:,0,iBlock)  = State_VGB(:,:,:,1,iBlock)
          State_VGB(:,:,:,-1,iBlock)  = State_VGB(:,:,:,1,iBlock)
          ! bin 1 is starts on the left
          State_VGB(WaveFirst_,:,:,-1:0,iBlock) = 1.0
       case('usersemi')
          ! bin 1 is starts on the left
          do j = 1, nJ; do i = 1, nI
             StateSemi_VGB(iTrImplFirst,i,j,0,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
             StateSemi_VGB(iTrImplLast,i,j,0,iBlock) = &
                  StateSemi_VGB(iTrImplLast,i,j,1,iBlock)
          end do; end do
       end select
    case(6)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,:,:,nK+1,iBlock) = State_VGB(:,:,:,nK,iBlock)
          State_VGB(:,:,:,nK+2,iBlock) = State_VGB(:,:,:,nK,iBlock)
          ! bin 2 starts on the right
          State_VGB(WaveLast_,:,:,nK+1:nK+2,iBlock) = 1.0
       case('usersemi')
          ! bin 2 starts on the right
          do j = 1, nJ; do i = 1, nI
             StateSemi_VGB(iTrImplFirst,i,j,nK+1,iBlock) = &
                  StateSemi_VGB(iTrImplFirst,i,j,nK,iBlock)
             StateSemi_VGB(iTrImplLast,i,j,nK+1,iBlock) = &
                  sqrt(sqrt(1.0/cRadiationNo))
          end do; end do
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
    use ModConst,      ONLY: cKToEv
    use ModMain,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, UnitEnergyDens_
    use ModVarIndexes, ONLY: p_, Rho_, WaveFirst_, WaveLast_

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

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('erad1')
!       NameIdlUnit = 'erg/cm3'
       PlotVar_G = &
            State_VGB(WaveFirst_,:,:,:,iBlock) !*No2Si_V(UnitEnergyDens_)
    case('erad2')
!       NameIdlUnit = 'erg/cm3'
       PlotVar_G = &
            State_VGB(WaveLast_,:,:,:,iBlock) !*No2Si_V(UnitEnergyDens_)
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, &
       PlanckOut_W, CgTeOut_W, CgTgOut_W, TgOut_W)

    ! The State_V vector is in normalized units

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

    real :: Rho, Pressure, Te, Tg_W(nWave)
    real :: RhoSi, pSi, TeSi

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

    if(present(OpacityPlanckOut_W)) OpacityPlanckOut_W = 0.0

    if(present(OpacityRosselandOut_W)) &
         OpacityRosselandOut_W = 1.0e-16/No2Si_V(UnitX_)

    if(present(HeatCondOut)) HeatCondOut = 0.0
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0

    Tg_W = sqrt(sqrt(State_V(WaveFirst_:WaveLast_)/cRadiationNo))

    if(present(PlanckOut_W)) PlanckOut_W = cRadiationNo*Te**4 &
         *No2Si_V(UnitEnergyDens_)

    if(present(CgTeOut_W)) CgTeOut_W = cRadiationNo*(Te+Tg_W)*(Te**2+Tg_W**2) &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(TgOut_W)) TgOut_W = Tg_W*No2Si_V(UnitTemperature_)

    if(present(CgTgOut_W)) CgTgOut_W = 4.0*cRadiationNo*Tg_W**3 &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

  end subroutine user_material_properties

end module ModUser
