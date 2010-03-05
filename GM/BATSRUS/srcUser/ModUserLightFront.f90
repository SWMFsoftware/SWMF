!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_init_session,               &
       IMPLEMENTED2 => user_read_inputs,                &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_update_states,              &
       IMPLEMENTED6 => user_set_plot_var,               &
       IMPLEMENTED7 => user_material_properties

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'Marshak Wave'

  character(len=20) :: TypeProblem

  real :: TradBc, Density, SpecificOpacity, Epsilon

contains

  !============================================================================

  subroutine user_init_session

    use ModConst,   ONLY: cKevToK
    use ModPhysics, ONLY: Si2No_V, UnitTemperature_, UnitRho_, UnitX_

    character (len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    TradBc = cKevToK*Si2No_V(UnitTemperature_)                ! 1 keV
    Density = 1000.0*Si2No_V(UnitRho_)                        ! 1 g/cm^3
    SpecificOpacity = 1.0e-18/(Si2No_V(UnitRho_)*Si2No_V(UnitX_)) ! 1 cm^2/g
    Epsilon = 1.0

  end subroutine user_init_session

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

!esm [
    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk
    use ModConst,      ONLY: cKToEv
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, cRadiationNo
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Erad_, ExtraEint_, p_
!esm ]
    integer, intent(in) :: iStage, iBlock
!esm [
    integer ::  i,j,k
!esm ]
    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    ! No call to update_states_MHD to nullify the effect of the hydro solver
    ! call update_states_MHD(iStage,iBlock)
!esm [
!    do i = -1, nI+2
!       write(6,*) 'i, tmat = ', i, &
!            (State_VGB(p_,i,nJ/2,nK/2,iBlock)/ &
!            State_VGB(Rho_,i,nJ/2,nK/2,iBlock)) &
!            *No2Si_V(UnitTemperature_)*cKToEv
!       write(6,*) 'i, erad = ', i, &
!            x_Blk(i,1,1,iBlock),   &
!            State_VGB(Erad_,i,nJ/2,nK/2,iBlock)
!    enddo
!esm ]
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
       case("#LIGHTFRONTTEST")
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
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Erad_, ExtraEint_, p_

    integer :: i, j, k, iBlock
    real :: Temperature, Pressure

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    select case(TypeProblem)
    case('cartesian')
!!       Temperature = 1.0e-5*cKevToK*Si2No_V(UnitTemperature_)
       Temperature = 1.0e-16
!       Density     = 1.0e-16
       Pressure = Density*Temperature
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          State_VGB(Rho_,i,j,k,iBlock) = Density
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
!!!          State_VGB(Erad_,i,j,k,iBlock) = cRadiationNo*Temperature**4
!!!          State_VGB(Erad_,i,j,k,iBlock) = 1.0d-06*cRadiationNo*Temperature**4
          State_VGB(Erad_,i,j,k,iBlock) = 1.0d-06
!!!          State_VGB(ExtraEint_,i,j,k,iBlock) = &
!!!               cRadiationNo*Temperature**4 - inv_gm1*Pressure
          State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0

          State_VGB(p_,i,j,k,iBlock) = Pressure
       end do; end do; end do

    case default
       call stop_mpi(NameSub//' : undefined problem type='//TypeProblem)
    end select

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,  ONLY: State_VGB
    use ModGeometry, ONLY: Dx_Blk
    use ModImplicit, ONLY: StateSemi_VGB, TypeSemiImplicit
    use ModMain,     ONLY: nJ, nK
    use ModPhysics,  ONLY: No2Si_V, Si2No_V, UnitEnergyDens_, &
         cRadiationNo

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: j, k, iEradImpl
    real :: Coef

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------
    if(iSide /= 1)then
       write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
       call stop_mpi(NameSub)
    end if

    IsFound = .true.

    select case(TypeBc)
    case('user')

       State_VGB(:,0,:,:,iBlock)  = 1.0
       State_VGB(:,-1,:,:,iBlock) = 1.0
       State_VGB(:,1,:,:,iBlock) = 1.0


    case('usersemi')

       ! Marshak boundary conditions

       select case(TypeSemiImplicit)
       case('radiation')
          iEradImpl = 1
       case('radcond')
          iEradImpl = 2
       end select

       Coef = 2.0/(3.0*SpecificOpacity*Density*Dx_Blk(iBlock))
       do k = 1, nK; do j = 1, nJ
          StateSemi_VGB(iEradImpl,0:1,j,k,iBlock) = 1.0 !&
!               (cRadiationNo*TradBc**4 - StateSemi_VGB(iEradImpl,1,j,k,iBlock)&
!               *(0.5 - Coef) )/(0.5 + Coef)
       end do; end do
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
    use ModVarIndexes, ONLY: p_, Rho_, Erad_

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
    case('tmat')
       NameIdlUnit = 'eV'
       PlotVar_G = (State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)) &
            *No2Si_V(UnitTemperature_)*cKToEv
    case('erad')
       NameIdlUnit = 'erg/cm3'
       PlotVar_G = &
            State_VGB(Erad_,:,:,:,iBlock)*No2Si_V(UnitEnergyDens_)
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

    ! The State_V vector is in normalized units

    use ModAdvance,    ONLY: nWave
    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_, cRadiationNo, Clight
    use ModVarIndexes, ONLY: nVar, Rho_, p_

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

    real :: Rho, Pressure, Temperature
    real :: RhoSi, pSi, TemperatureSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalIn))then
       pSi = EinternalIn*gm1
       Pressure = pSi*Si2No_V(UnitP_)
       Temperature = Pressure/Rho
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
    elseif(present(TeIn))then
       TemperatureSi = TeIn
       Temperature = TemperatureSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Temperature
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Temperature = Pressure/Rho
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalOut)) EinternalOut = pSi*inv_gm1
    if(present(TeOut)) TeOut = TemperatureSi
    if(present(PressureOut)) PressureOut = pSi

    if(present(CvOut)) CvOut = inv_gm1*Rho &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(OpacityPlanckOut_W)) &
         OpacityPlanckOut_W = 0.0 !SpecificOpacity*Rho/No2Si_V(UnitX_)
    if(present(OpacityRosselandOut_W)) &
         OpacityRosselandOut_W = SpecificOpacity*Rho/No2Si_V(UnitX_)
    if(present(HeatCondOut)) HeatCondOut = 0.0
    if(present(TeTiRelaxOut)) TeTiRelaxOut = 0.0

  end subroutine user_material_properties

end module ModUser
