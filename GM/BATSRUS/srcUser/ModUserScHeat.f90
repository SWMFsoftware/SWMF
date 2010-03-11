!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModReadParam, ONLY: lStringLine
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_get_log_var,                &
       IMPLEMENTED5 => user_calc_sources,               &
       IMPLEMENTED6 => user_specify_refinement,         &
       IMPLEMENTED7 => user_set_plot_var,               &
       IMPLEMENTED8 => user_set_outerbcs,               &
       IMPLEMENTED9 => user_face_bcs,                   &
       IMPLEMENTED10=> user_set_boundary_cells,         &
       IMPLEMENTED11=> user_set_resistivity

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Global Corona'

  logical ::  UseWaveDissipation = .false.
  real :: DissipationScaleFactor  ! unit = m*T^0.5

  real :: TeFraction, TiFraction
  real :: EtaPerpSi

  character(len=lStringLine) :: TypeCoronalHeating

contains

  !============================================================================

  subroutine user_read_inputs

    use ModMain,           ONLY: UseUserInitSession, lVerbose
    use ModProcMH,         ONLY: iProc
    use ModReadParam,      ONLY: read_line, read_command, read_var
    use ModIO,             ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand
    character(len=*), parameter :: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    UseUserInitSession = .true.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input SOLAR CORONA starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)

       case("#WSACOEFF")
          call read_wsa_coeff

       case("#CORONALHEATING")
          call read_var('TypeCoronalHeating', TypeCoronalHeating)
          select case(TypeCoronalHeating)
          case('wavedissipation')
             UseWaveDissipation = .true.
             call read_var('DissipationScaleFactor', DissipationScaleFactor)
          case default
             call stop_mpi(NameSub//': unknown TypeCoronalHeating = ' &
                  //TypeCoronalHeating)
          end select

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input SOLAR CORONA ends'
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

    use ModAdvance,     ONLY: UseElectronPressure
    use ModConst,       ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
         cElectronMass
    use ModIO,          ONLY: write_prefix, iUnitOut
    use ModMultiFluid,  ONLY: MassIon_I
    use ModNumConst,    ONLY: cTwoPi
    use ModPhysics,     ONLY: ElectronTemperatureRatio, AverageIonCharge
    use ModProcMH,      ONLY: iProc
    use ModWaves,       ONLY: UseWavePressure, UseAlfvenWaves


    real, parameter :: CoulombLog = 20.0
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    UseAlfvenWaves = .true.
    UseWavePressure = .true.

    ! TeFraction is used for ideal EOS:
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TiFraction = MassIon_I(1)
       TeFraction = MassIon_I(1)/AverageIonCharge
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       TeFraction = TiFraction*ElectronTemperatureRatio
    end if

    ! perpendicular resistivity, used for temperature relaxation
    ! Note EtaPerpSi is divided by cMu.
    EtaPerpSi = sqrt(cElectronMass)*CoulombLog &
         *(cElectronCharge*cLightSpeed)**2/(3*(cTwoPi*cBoltzmann)**1.5*cEps)

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session

  !============================================================================

  subroutine get_plasma_parameters_base(x_D, RhoBase, Tbase)

    ! This subroutine computes the base values for mass density and temperature

    use ModPhysics,          ONLY: BodyRho_I, Si2No_V, UnitTemperature_
    use ModExpansionFactors, ONLY: Umin, CoronalT0Dim

    real, intent(in) :: x_D(3)
    real, intent(out):: RhoBase, Tbase

    real :: Ufinal ! The solar wind speed at the far end of the Parker spiral,
                   ! which originates from the given point.
    real :: Uratio ! The coronal based values for temperature and density
                   ! are scaled as functions of UFinal/UMin ratio.
    real :: Runit_D(3)
    !--------------------------------------------------------------------------

    Runit_D = x_D/sqrt(sum(x_D**2))

    call get_bernoulli_integral(Runit_D(1), Runit_D(2), Runit_D(3), Ufinal)
    Uratio = Ufinal/Umin

    ! This is the temperature variation.
    ! In coronal holes the temperature is reduced 
    Tbase = CoronalT0Dim*Si2No_V(UnitTemperature_) / min(Uratio, 1.5)

    ! This is the density variation
    RhoBase = BodyRho_I(1)/URatio

  end subroutine get_plasma_parameters_base

  !============================================================================

  subroutine get_total_wave_energy(x, y, z, VAlfvenSi, WaveEnergyDensSi)

    ! Provides the distribution of the total Alfven wave energy density
    ! at the coronal base complying with the WSA semi-empirical model
    ! Borrowed from ModExpansionFactors
    ! No nonzero heat conduction contribution at base yet

    use ModExpansionFactors
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: g, inv_gm1, AverageIonCharge

    real, intent(in) :: x, y, z, VAlfvenSi !VAlfven should be in m/s
    real, intent(out):: WaveEnergyDensSi

    real :: r, Uf, ExpansionFactorInv
    real, parameter :: RhoVAt1AU = 5.40e-15 !kg/(m2*s)
    real, parameter :: AreaRatio = (cAU/rSun)**2
    real, parameter :: RhoV =  AreaRatio * RhoVAt1AU
    real, parameter :: VAlfvenMin = 1.0e5   !100 km/s

    real :: HeatFluxSi
    !--------------------------------------------------------------------------

    !\
    ! Calculate cell-centered spherical coordinates::
    r = sqrt(x**2 + y**2 + z**2)
    !\
    ! Avoid calculating inside a critical radius = 0.5*Rsun
    !/
    if (r < max(Ro_PFSSM-dR*nRExt,0.90*Ro_PFSSM)) then 
       WaveEnergyDensSi = 0.0
       RETURN
    end if

    !v_\infty from WSA model:
    call get_bernoulli_integral(x, y, z, Uf)

    !An expansion factor
    call get_interpolated(ExpansionFactorInv_N, x, y, z, ExpansionFactorInv)

    HeatFluxSi = 0.0

    WaveEnergyDensSi = (RhoV*(0.5*Uf**2 + cSunGravitySi - g*inv_gm1*&
         cBoltzmann/(cProtonMass*MassIon_I(1))*(1.0+AverageIonCharge) &
         *CoronalT0Dim/min(Uf/UMin, 1.5) ) & !This is a modulated Tc
         - ExpansionFactorInv*HeatFluxSi) &
         /max(abs(VAlfvenSi)*ExpansionFactorInv, VAlfvenMin)

  end subroutine get_total_wave_energy

  !============================================================================

  subroutine user_set_ics

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk, r_Blk, true_cell
    use ModMain,       ONLY: nI, nJ, nK, globalBLK
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         BodyRho_I, BodyTDim_I, UnitU_, AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_, Pe_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    integer :: IterCount
    real :: x, y, z, r, RhoBase, Rho, NumDensIon, NumDensElectron
    real :: Tcorona, Tbase, Temperature
    real :: Ur, Ur0, Ur1, del, Ubase, rTransonic, Uescape, Usound

    real, parameter :: Epsilon = 1.0e-6
    !--------------------------------------------------------------------------

    iBlock = globalBLK

    ! Initially, the electron and ion temperature are at 1.5e6(K) in the corona
    Tcorona = 1.5e6*Si2No_V(UnitTemperature_)

    ! normalize with isothermal sound speed.
    Usound = sqrt(Tcorona*(1.0+AverageIonCharge)/MassIon_I(1))
    Uescape = sqrt(-GBody*2.0)/Usound

    !\
    ! Initialize MHD wind with Parker's solution
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    rTransonic = 0.25*Uescape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    Ubase = rTransonic**2*exp(1.5 - 2.0*rTransonic)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)
       z = z_BLK(i,j,k,iBlock)
       r = r_BLK(i,j,k,iBlock)

       ! The electron and ion temperature are initially equal to
       ! the temperature at the base.
       call get_plasma_parameters_base((/x,y,z/), RhoBase, Tbase)

       if(r > rTransonic)then
          !\
          ! Inside supersonic region
          !/
          Ur0 = 1.0
          IterCount = 0
          do
             IterCount = IterCount + 1
             Ur1 = sqrt(Uescape**2/r - 3.0 + 2.0*log(16.0*Ur0*r**2/Uescape**4))
             del = abs(Ur1 - Ur0)
             if(del < Epsilon)then
                Ur = Ur1
                EXIT
             elseif(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end do
       else
          !\
          ! Inside subsonic region
          !/
          Ur0 = 1.0
          IterCount = 0
          do
             IterCount = IterCount + 1
             Ur1 = (Uescape**2/(4.0*r))**2 &
                  *exp(0.5*(Ur0**2 + 3.0 - Uescape**2/r))
             del = abs(Ur1 - Ur0)
             if(del < Epsilon)then
                Ur = Ur1
                EXIT
             elseif(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end do
       end if

       Rho = rBody**2*RhoBase*Ubase/(r**2*Ur)
       State_VGB(Rho_,i,j,k,iBlock) = Rho
       State_VGB(RhoUx_,i,j,k,iBlock) = Rho*Ur*x/r *Usound
       State_VGB(RhoUy_,i,j,k,iBlock) = Rho*Ur*y/r *Usound
       State_VGB(RhoUz_,i,j,k,iBlock) = Rho*Ur*z/r *Usound
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 1.0e-12 !0.0
       NumDensIon = Rho/MassIon_I(1)
       NumDensElectron = NumDensIon*AverageIonCharge

       if(UseElectronPressure)then
          State_VGB(p_,i,j,k,iBlock) = NumDensIon*Tbase
          State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*Tbase
       else
          State_VGB(p_,i,j,k,iBlock) = &
               (NumDensIon + NumDensElectron)*Tbase
       end if
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_calc_sources

    use ModAdvance,        ONLY: State_VGB, Source_VC, UseElectronPressure, &
         time_BLK, B0_DGB
    use ModGeometry,       ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModMain,           ONLY: nI, nJ, nK, GlobalBlk, Cfl, UseB0
    use ModPhysics,        ONLY: gm1, inv_gm1, rBody, Si2No_V, No2Si_V, &
         UnitT_, UnitB_, UnitX_
    use ModVarIndexes,     ONLY: Rho_, Bx_, Bz_, Energy_, p_, Pe_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock, iWave
    real :: CoronalHeating, RadiativeCooling
    real :: DissipationLength, WaveHeating, FullB_D(3), FullBSi

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------

    iBlock = globalBlk

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(r_BLK(i,j,k,iBlock) < rBody) CYCLE

       if(UseWaveDissipation)then
          if(UseB0)then
             FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
          else
             FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if
          FullBSi = sqrt(sum(FullB_D**2))*No2Si_V(UnitB_)
          DissipationLength = DissipationScaleFactor &
               /sqrt(FullBSi)*Si2No_V(UnitX_)
          CoronalHeating = 0.0
          do iWave = WaveFirst_, WaveLast_
             WaveHeating = State_VGB(iWave,i,j,k,iBlock)**1.5 &
                  /sqrt(State_VGB(Rho_,i,j,k,iBlock))/DissipationLength
             Source_VC(iWave,i,j,k) = Source_VC(iWave,i,j,k) - WaveHeating
             CoronalHeating = CoronalHeating + WaveHeating
          end do
       else
          CoronalHeating = 0.0
       end if
       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + CoronalHeating

       CYCLE

       call get_radiative_cooling( &
            State_VGB(:,i,j,k,iBlock), RadiativeCooling)

       if(UseElectronPressure)then
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) + gm1*RadiativeCooling
       else
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k)+RadiativeCooling
       end if
    end do; end do; end do

  end subroutine user_calc_sources

  !============================================================================

  subroutine get_radiative_cooling(State_V, RadiativeCooling)

    use ModAdvance,    ONLY: UseElectronPressure
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitT_, UnitN_, &
         UnitEnergyDens_, UnitTemperature_, AverageIonCharge
    use ModVarIndexes, ONLY: nVar, Rho_, p_, Pe_

    real, intent(in) :: State_V(nVar)
    real, intent(out):: RadiativeCooling

    real :: Te, TeSi, CoolingFunctionCgs
    real :: NumDensIonCgs, NumDensElectronCgs
    !--------------------------------------------------------------------------

    if(UseElectronPressure)then
       Te = TeFraction*State_V(Pe_)/State_V(Rho_)
    else
       Te = TeFraction*State_V(p_)/State_V(Rho_)
    end if
    TeSi =Te*No2Si_V(UnitTemperature_)

    ! CGS is used to avoid insane numbers
    call get_cooling_function(TeSi, CoolingFunctionCgs)
    NumDensIonCgs = (State_V(Rho_)/MassIon_I(1))*No2Si_V(UnitN_)*1.0e-6
    NumDensElectronCgs = AverageIonCharge*NumDensIonCgs

    RadiativeCooling = -NumDensIonCgs*NumDensElectronCgs*CoolingFunctionCgs &
         *0.1*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

  contains

    subroutine get_cooling_function(TeSi, CoolingFunctionCgs)

      ! Based on Rosner et al. (1978) and Peres et al. (1982)
      ! Need to be replaced by Chianti tables

      real, intent(in) :: TeSi
      real, intent(out):: CoolingFunctionCgs
      !------------------------------------------------------------------------

      if(TeSi <= 8e3)then
         CoolingFunctionCgs = (1.0606e-6*TeSi)**11.7
      elseif(TeSi <= 2e4)then
         CoolingFunctionCgs = (1.397e-8*TeSi)**6.15
      elseif(TeSi <= 10**4.575)then
         CoolingFunctionCgs = 10**(-21.85)
      elseif(TeSi <= 10**4.9)then
         CoolingFunctionCgs = 10**(-31.0)*TeSi**2
      elseif(TeSi <= 10**5.4)then
         CoolingFunctionCgs = 10**(-21.2)
      elseif(TeSi <= 10**5.77)then
         CoolingFunctionCgs = 10**(-10.4)/TeSi**2
      elseif(TeSi <= 10**6.315)then
         CoolingFunctionCgs = 10**(-21.94)
      elseif(TeSi <= 10**7.60457)then
         CoolingFunctionCgs = 10**(-17.73)/TeSi**(2.0/3.0)
      else
         CoolingFunctionCgs = 10**(-26.6)*sqrt(TeSi)
      end if

    end subroutine get_cooling_function

  end subroutine get_radiative_cooling

  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, B0_DGB, UseElectronPressure
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK, nBlock, x_, y_, z_
    use ModPhysics,    ONLY: inv_gm1, No2Io_V, UnitEnergydens_, UnitX_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_, p_, Pe_

    real, intent(out) :: VarValue
    character(len=10), intent(in) :: TypeVar 
    real, optional, intent(in) :: Radius

    integer :: iBlock
    real :: unit_energy
    real, external :: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = No2Io_V(UnitEnergydens_)*No2Io_V(UnitX_)**3
    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)
    case('eint')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          if(UseElectronPressure)then
             tmp1_BLK(:,:,:,iBlock) = &
                  State_VGB(p_,:,:,:,iBlock) + State_VGB(Pe_,:,:,:,iBlock)
          else
             tmp1_BLK(:,:,:,iBlock) = State_VGB(p_,:,:,:,iBlock)
          end if
       end do
       VarValue = unit_energy*inv_gm1*integrate_BLK(1,tmp1_BLK)

    case('emag')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = & 
               ( B0_DGB(x_,:,:,:,iBlock) + State_VGB(Bx_,:,:,:,iBlock))**2 &
               +(B0_DGB(y_,:,:,:,iBlock) + State_VGB(By_,:,:,:,iBlock))**2 &
               +(B0_DGB(z_,:,:,:,iBlock) + State_VGB(Bz_,:,:,:,iBlock))**2
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)

    case('vol')
       tmp1_BLK(:,:,:,iBlock) = 1.0
       VarValue = integrate_BLK(1,tmp1_BLK)

    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select

  end subroutine user_get_log_var

  !============================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_
    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: Rho_, p_, Pe_

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

    integer :: i, j, k

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    IsFound = .true.

    select case(NameVar)
    case('te')
       NameIdlUnit = 'K'
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          if(UseElectronPressure)then
             PlotVar_G(i,j,k) = TeFraction*State_VGB(Pe_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
          else
             PlotVar_G(i,j,k) = TeFraction*State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
          end if
       end do; end do; end do
    case('ti')
       NameIdlUnit = 'K'
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          PlotVar_G(i,j,k) = TiFraction*State_VGB(p_,i,j,k,iBlock) &
               /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
       end do; end do; end do
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, B0_DGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, far_field_BCs_BLK
    use ModNumConst, ONLY: cTiny
    use ModMain,     ONLY: x_, y_, z_, time_loop

    integer, intent(in) :: iBlock, iArea
    logical,intent(out) :: DoRefine

    real :: rDotB_G(1:nI,1:nJ,0:nK+1)
    integer :: i, j, k
    character (len=*), parameter :: NameSub = 'user_specify_refinement'
    !--------------------------------------------------------------------------

    if(.not.time_loop)then
       DoRefine = .false.

       RETURN
    end if

    ! Calculate r.B in all physical cells and ghost cells 
    ! in the Z/Theta direction to find current sheet 
    ! passing between blocks
    do k=0, nK+1; do j=1, nJ; do i=1, nI
       rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
            * (B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock)) &
            +           y_BLK(i,j,k,iBlock)   &
            * (B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock)) &
            +           z_BLK(i,j,k,iBlock)   &
            * (B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))
    end do; end do; end do;

    DoRefine = maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny

  end subroutine user_specify_refinement

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    ! Fill one layer of ghost cells with the temperature for heat conduction

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk
    use ModMain,       ONLY: x_, y_, z_, nJ, nK, East_
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, p_, Pe_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k
    real :: x_D(3)
    real :: RhoBase, Tbase, NumDensIon, NumDensElectron

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    if(iSide /= East_) call stop_mpi('Wrong iSide in user_set_outerBCs')

    IsFound = .true.

    do k = -1, nK+2; do j = -1, nJ+2
       x_D(x_) = 0.5*sum(x_Blk(0:1,j,k,iBlock))
       x_D(y_) = 0.5*sum(y_Blk(0:1,j,k,iBlock))
       x_D(z_) = 0.5*sum(z_Blk(0:1,j,k,iBlock))

       call get_plasma_parameters_base(x_D, RhoBase, Tbase)

       ! Fixed density
       State_VGB(Rho_,0,j,k,iBlock) = &
            2.0*RhoBase - State_VGB(Rho_,1,j,k,iBlock)
       State_VGB(Rho_,-1,j,k,iBlock) = State_VGB(Rho_,0,j,k,iBlock)

       ! fixed electron and ion temperature
       do i = -1, 0
          NumDensIon = State_VGB(Rho_,i,j,k,iBlock)/MassIon_I(1)
          NumDensElectron = NumDensIon*AverageIonCharge
          if(UseElectronPressure)then
             State_VGB(p_,i,j,k,iBlock) = NumDensIon*Tbase
             State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*Tbase
          else
             State_VGB(p_,i,j,k,iBlock) = (NumDensIon + NumDensElectron)*Tbase
          end if
       end do

    end do; end do

  end subroutine user_set_outerbcs

  !============================================================================

  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance,     ONLY: State_VGB, UseElectronPressure
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: OmegaBody, BodyRho_I, BodyTDim_I, No2Si_V, &
         UnitTemperature_, Si2No_V, AverageIonCharge, UnitU_, UnitEnergyDens_
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_, &
         WaveFirst_, WaveLast_, Pe_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: RhoBase, NumDensIon, NumDensElectron, Tbase, FullBr
    real :: Runit_D(3), U_D(3)
    real :: B1_D(3), B1t_D(3), B1r_D(3), FullB_D(3)
    real :: UalfvenSi, EwaveSi, Ewave
    !--------------------------------------------------------------------------

    Runit_D = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    U_D   = VarsTrueFace_V(Ux_:Uz_)
    B1_D  = VarsTrueFace_V(Bx_:Bz_)
    B1r_D = dot_product(Runit_D, B1_D)*Runit_D
    B1t_D = B1_D - B1r_D

    VarsGhostFace_V(Ux_:Uz_) = -U_D
    VarsGhostFace_V(Bx_:Bz_) = B1t_D !- B1r_D

    FullB_D = B0Face_D + B1t_D
    FullBr = dot_product(Runit_D, FullB_D)

    call get_plasma_parameters_base(FaceCoords_D, RhoBase, Tbase)

    VarsGhostFace_V(Rho_) =  2.0*RhoBase - VarsTrueFace_V(Rho_)
    NumDensIon = VarsGhostFace_V(Rho_)/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    if(UseElectronPressure)then
       VarsGhostFace_V(p_) = NumDensIon*Tbase
       VarsGhostFace_V(Pe_) = NumDensElectron*Tbase
    else
       VarsGhostFace_V(p_) = (NumDensIon + NumDensElectron)*Tbase
    end if

    ! Set Alfven waves energy density based on Bernoulli function
    UalfvenSi = (FullBr/sqrt(RhoBase))*No2Si_V(UnitU_)
    call get_total_wave_energy( &
         FaceCoords_D(x_), &
         FaceCoords_D(y_), &
         FaceCoords_D(z_), &
         UalfvenSi, EwaveSi)
    Ewave = EwaveSi*Si2No_V(UnitEnergyDens_)

    if(UalfvenSi > 0.0)then
       VarsGhostFace_V(WaveFirst_) = Ewave
       VarsGhostFace_V(WaveLast_) = VarsTrueFace_V(WaveLast_)
    else
       VarsGhostFace_V(WaveFirst_) = VarsTrueFace_V(WaveFirst_)
       VarsGhostFace_V(WaveLast_) = Ewave
    end if

    !\
    ! Apply corotation if needed
    !/
    if(.not.UseRotatingFrame)then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2.0*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2.0*OmegaBody*FaceCoords_D(x_)
    end if

  end subroutine user_face_bcs

  !============================================================================

  subroutine user_set_boundary_cells(iBLK)

    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, r_Blk
    use ModBoundaryCells, ONLY: SaveBoundaryCells
    use ModPhysics,       ONLY: rBody

    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = r_Blk(:,:,:,iBLK) < rBody

    if(SaveBoundaryCells) RETURN
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

  !============================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: nI, nJ, nK
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitTemperature_, UnitX_, UnitT_
    use ModVarIndexes, ONLY: Rho_, Pe_

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(-1:nI+2,-1:nJ+2,-1:nK+2)

    integer :: i, j, k
    real :: Te, TeSi

    character (len=*), parameter :: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------

    do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
       Te = TeFraction*State_VGB(Pe_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       TeSi = Te*No2Si_V(UnitTemperature_)

       Eta_G(i,j,k) = EtaPerpSi/TeSi**1.5 *Si2No_V(UnitX_)**2/Si2No_V(UnitT_)
    end do; end do; end do

  end subroutine user_set_resistivity

end module ModUser
