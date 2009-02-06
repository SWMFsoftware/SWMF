!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_init_session,               &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_update_states,              &
       IMPLEMENTED4 => user_set_plot_var,               &
       IMPLEMENTED5 => user_material_properties

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Mhd with heat conduction'

  real, parameter :: HeatConductionCoef = 10.0
  real, parameter :: AmplitudeTemperature = 1e3
  real, parameter :: T0 = 1.0

contains

  !============================================================================

  subroutine user_init_session

    use ModGeometry, ONLY: x1, x2, XyzMin_D, XyzMax_D
    use ModMain,     ONLY: nI
    use ModParallel, ONLY: proc_dims

    character(len=*), parameter :: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------

    if(x2 /= -x1) &
       call stop_mpi(NameSub//' : xMin should be -xMax in the #GRID command')

    ! Shift grid to the right by half a grid cell size,
    ! so that x=0 is at a cell center
    x1 = x1 - x1/(2*nI*proc_dims(1))
    x2 = x2 + x2/(2*nI*proc_dims(2))

    XyzMin_D(1) = x1
    XyzMax_D(1) = x2

  end subroutine user_init_session

  !============================================================================

  subroutine user_update_states(iStage, iBlock)

    use ModAdvance,  ONLY: VdtFace_x, VdtFace_y, VdtFace_z
    use ModGeometry, ONLY: fAx_Blk

    integer, intent(in) :: iStage, iBlock

    real, parameter :: Speed = 10.

    character(len=*), parameter :: NameSub = 'user_update_states'
    !--------------------------------------------------------------------------

    ! No call to update_states_MHD to nullify the effect of the hydro solver
    ! call update_states_MHD(iStage,iBlock)

    ! the following results in a time step of Dx/Speed
    VdtFace_x = fAx_Blk(iBlock)*Speed
    VdtFace_y = 0.0
    VdtFace_z = 0.0

  end subroutine user_update_states

  !============================================================================

  subroutine user_set_ics

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, Dx_Blk
    use ModMain,       ONLY: GlobalBlk, nI, nJ, nK, x_, Time_Simulation
    use ModNumConst,   ONLY: cPi
    use ModPhysics,    ONLY: gm1
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, ExtraEint_, p_

    integer :: iBlock, i, j, k
    real :: x, Dx, Temperature, Spread

    character(len=*), parameter :: NameSub = "user_set_ics"
    !--------------------------------------------------------------------------
    iBlock = GlobalBlk

    ! temperature spike at x=0, otherwise set the temperature to zero
    do i = -1, nI+2
       if(Time_Simulation == 0.0)then
          if(x_Blk(i,1,1,iBlock) > -0.5*Dx_Blk(iBlock) .and. &
               x_Blk(i,1,1,iBlock) < 0.5*Dx_Blk(iBlock)) then
             Temperature = T0 + AmplitudeTemperature/Dx_Blk(iBlock)
          else
             Temperature = T0
          end if
       else
          Spread = 4.0*HeatConductionCoef*Time_Simulation
          Temperature = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
               *exp(-x_Blk(i,1,1,iBlock)**2/Spread)
       end if
       do k = -1, nK+2; do j = -1, nJ+2
          State_VGB(Rho_,i,j,k,iBlock) = 1.0
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          State_VGB(ExtraEint_,i,j,k,iBlock) = 0.0
          State_VGB(p_,i,j,k,iBlock) = Temperature
       end do; end do
    end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_Blk, Dx_Blk
    use ModMain,       ONLY: nI, nJ, nK, Time_Simulation
    use ModNumConst,   ONLY: cPi
    use ModVarIndexes, ONLY: p_, Rho_

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

    real :: Spread

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('te')
       PlotVar_G = State_VGB(p_,:,:,:,iBlock)/State_VGB(Rho_,:,:,:,iBlock)
    case('te0')
       if(Time_Simulation == 0.0)then
          PlotVar_G = State_VGB(p_,:,:,:,iBlock) &
               /State_VGB(Rho_,:,:,:,iBlock)
       else
          Spread = 4.0*HeatConductionCoef*Time_Simulation
          PlotVar_G = T0 + AmplitudeTemperature/(sqrt(cPi*Spread)) &
               *exp(-x_Blk(:,:,:,iBlock)**2/Spread)
       end if
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, EinternalSiIn, &
       TeSiIn, EinternalSiOut, TeSiOut, PressureSiOut, CvSiOut, &
       AbsorptionOpacitySiOut, RosselandMeanOpacitySiOut, &
       HeatConductionCoefSiOut)

    ! The State_V vector is in normalized units

    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_
    use ModVarIndexes, ONLY: nVar, Rho_, p_

    real, intent(in) :: State_V(nVar)
    real, optional, intent(in)  :: EinternalSiIn             ! [J/m^3]
    real, optional, intent(in)  :: TeSiIn                    ! [K]
    real, optional, intent(out) :: EinternalSiOut            ! [J/m^3]
    real, optional, intent(out) :: TeSiOut                   ! [K]
    real, optional, intent(out) :: AbsorptionOpacitySiOut    ! [1/m]
    real, optional, intent(out) :: RosselandMeanOpacitySiOut ! [1/m]
    real, optional, intent(out) :: CvSiOut                   ! [J/(K*m^3)]
    real, optional, intent(out) :: PressureSiOut             ! [Pa]
    real, optional, intent(out) :: HeatConductionCoefSiOut   ! [Jm^2/(Ks)]

    real :: Rho, Pressure, Tmat
    real :: RhoSi, pSi, TmatSi

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(EinternalSiIn))then
       pSi = EinternalSiIn*gm1
       Pressure = pSi*Si2No_V(UnitP_)
       Tmat = Pressure/Rho
       TmatSi = Tmat*No2Si_V(UnitTemperature_)
    elseif(present(TeSiIn))then
       TmatSi = TeSiIn
       Tmat = TmatSi*Si2No_V(UnitTemperature_)
       Pressure = Rho*Tmat
       pSi = Pressure*No2Si_V(UnitP_)
    else
       Pressure = State_V(p_)
       pSi = Pressure*No2Si_V(UnitP_)
       Tmat = Pressure/Rho
       TmatSi = Tmat*No2Si_V(UnitTemperature_)
    end if

    if(present(EinternalSiOut)) EinternalSiOut = pSi*inv_gm1
    if(present(TeSiOut)) TeSiOut = TmatSi
    if(present(PressureSiOut)) PressureSiOut = pSi

    if(present(CvSiOut)) CvSiOut = inv_gm1*Rho &
         *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_)

    if(present(HeatConductionCoefSiOut)) &
       HeatConductionCoefSiOut = HeatConductionCoef &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
            *No2Si_V(UnitU_)*No2Si_V(UnitX_)

    if(present(AbsorptionOpacitySiOut)) AbsorptionOpacitySiOut = 0.0
    if(present(RosselandMeanOpacitySiOut)) RosselandMeanOpacitySiOut = 0.0

  end subroutine user_material_properties

end module ModUser
