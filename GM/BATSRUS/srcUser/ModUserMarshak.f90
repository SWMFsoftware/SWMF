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
    SpecificOpacity = 0.1/(Si2No_V(UnitRho_)*Si2No_V(UnitX_)) ! 1 cm^2/g
    Epsilon = 1.0

  end subroutine user_init_session

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    integer, intent(in) :: iStage, iBlock

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
       case("#MARSHAKWAVETEST")
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
    case('suolson')
       Temperature = 1.0e-5*cKevToK*Si2No_V(UnitTemperature_)
       Pressure = Density*Temperature
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          State_VGB(Rho_,i,j,k,iBlock) = Density
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          State_VGB(Erad_,i,j,k,iBlock) = cRadiationNo*Temperature**4
          State_VGB(ExtraEint_,i,j,k,iBlock) = &
               cRadiationNo*Temperature**4 - inv_gm1*Pressure
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
    use ModImplicit, ONLY: StateSemi_VGB, TypeSemiImplicit, iEradImpl
    use ModMain,     ONLY: nI, nJ, nK
    use ModPhysics,  ONLY: cRadiationNo

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    real :: Coef

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    select case(iSide)
    case(1)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,0,:,:,iBlock)  = State_VGB(:,1,:,:,iBlock)
          State_VGB(:,-1,:,:,iBlock) = State_VGB(:,1,:,:,iBlock)
       case('usersemi')
          ! Marshak boundary conditions
          Coef = 2.0/(3.0*SpecificOpacity*Density*Dx_Blk(iBlock))
          StateSemi_VGB(iEradImpl,0,:,:,iBlock) = &
               (cRadiationNo*TradBc**4 + StateSemi_VGB(iEradImpl,1,:,:,iBlock)&
               *(Coef - 0.5) )/(Coef + 0.5)
       case('usersemilinear')
          Coef = 2.0/(3.0*SpecificOpacity*Density*Dx_Blk(iBlock))
          StateSemi_VGB(iEradImpl,0,:,:,iBlock) = &
               StateSemi_VGB(iEradImpl,1,:,:,iBlock)*(Coef - 0.5)/(Coef + 0.5)
       end select
    case(2)
       select case(TypeBc)
       case('user')
          ! float, just for the sake of having filled in ghost cells
          State_VGB(:,nI+1,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
          State_VGB(:,nI+2,:,:,iBlock) = State_VGB(:,nI,:,:,iBlock)
       case('usersemi','usersemilinear')
          ! zero radiation influx boundary conditions
          Coef = 2.0/(3.0*SpecificOpacity*Density*Dx_Blk(iBlock))
          StateSemi_VGB(iEradImpl,nI+1,:,:,iBlock) = &
               StateSemi_VGB(iEradImpl,nI,:,:,iBlock)*(Coef - 0.5)/(Coef + 0.5)
       end select
    case default
       write(*,*) NameSub//' : user boundary not defined at iSide = ', iSide
       call stop_mpi(NameSub)
    end select

    IsFound = .true.

  end subroutine user_set_outerbcs

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModConst,      ONLY: cKToEv
    use ModMain,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, cRadiationNo
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
    case('trad')
       NameIdlUnit = 'eV'
       PlotVar_G = &
            sqrt(sqrt(abs(State_VGB(Erad_,:,:,:,iBlock))/cRadiationNo))&
            *No2Si_V(UnitTemperature_)*cKToEv
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalSiIn, TeSiIn, &
       EinternalSiOut, TeSiOut, PeSiOut, EeSiOut, PressureSiOut, &
       CvSiOut, CveSiOut, GammaOut, HeatCondSiOut, TeTiRelaxSiOut, &
       AbsorptionOpacitySiOut_I, DiffusionOpacitySiOut_I)

    ! The State_V vector is in normalized units

    use ModAdvance,    ONLY: nOpacity
    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_, cRadiationNo, Clight
    use ModVarIndexes, ONLY: nVar, Rho_, p_

    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in):: i, j, k, iBlock, iDir    ! cell/face index
    real, optional, intent(in)  :: EinternalSiIn             ! [J/m^3]
    real, optional, intent(in)  :: TeSiIn                    ! [K]
    real, optional, intent(out) :: EinternalSiOut            ! [J/m^3]
    real, optional, intent(out) :: TeSiOut                   ! [K]
    real, optional, intent(out) :: PeSiOut                   ! [Pa]
    real, optional, intent(out) :: EeSiOut                   ! [J/m^3]
    real, optional, intent(out) :: PressureSiOut             ! [Pa]
    real, optional, intent(out) :: CvSiOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: CveSiOut                  ! [J/(K*m^3)]
    real, optional, intent(out) :: GammaOut                  ! dimensionless
    real, optional, intent(out) :: HeatCondSiOut             ! [J/(m*K*s)]
    real, optional, intent(out) :: TeTiRelaxSiOut            ! [1/s]
    real, optional, intent(out) :: &
         AbsorptionOpacitySiOut_I(nOpacity)                  ! [1/m]
    real, optional, intent(out) :: &
         DiffusionOpacitySiOut_I(nOpacity)                   ! [1/m]

    real :: Rho, Pressure, Temperature
    real :: RhoSi, pSi, TemperatureSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalSiIn))then
       Temperature = sqrt(sqrt(EinternalSiIn*Si2No_V(UnitEnergyDens_) &
            /cRadiationNo))
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
       Pressure = Rho*Temperature
       pSi = Pressure*No2Si_V(UnitP_)
    elseif(present(TeSiIn))then
       TemperatureSi = TeSiIn
       Temperature = TemperatureSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Temperature
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Temperature = Pressure/Rho
       TemperatureSi = Temperature*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalSiOut)) EinternalSiOut = &
         cRadiationNo*Temperature**4*No2Si_V(UnitEnergyDens_)

    if(present(TeSiOut)) TeSiOut = TemperatureSi
    if(present(PressureSiOut)) PressureSiOut = pSi

    if(present(CvSiOut)) CvSiOut = 4.0*cRadiationNo*Temperature**3/Epsilon &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(AbsorptionOpacitySiOut_I)) &
         AbsorptionOpacitySiOut_I = SpecificOpacity*Rho/No2Si_V(UnitX_)
    if(present(DiffusionOpacitySiOut_I)) &
         DiffusionOpacitySiOut_I = SpecificOpacity*Rho/No2Si_V(UnitX_)

    if(present(HeatCondSiOut)) HeatCondSiOut = 0.0
    if(present(CveSiOut)) CveSiOut = 0.0
    if(present(TeTiRelaxSiOut)) TeTiRelaxSiOut = 0.0

  end subroutine user_material_properties

end module ModUser
