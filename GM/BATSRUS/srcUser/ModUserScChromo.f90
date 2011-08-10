!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain, ONLY: nI, nJ,nK
  use ModSize, ONLY: MaxBlock
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
       IMPLEMENTED11=> user_set_resistivity,            &
       IMPLEMENTED12=> user_update_states,              &
       IMPLEMENTED13=> user_initial_perturbation

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Chromosphere - solar wind model with Alfven waves - Oran, van der Holst'


  ! Input parameters for chromospheric inner BC's
  logical :: UseChromoBc = .true.
  real    :: WaveDeltaU = 0.0
  real    :: nChromoSi = 0.0, tChromoSi = 0.0
  real    :: nChromo = 0.0, RhoChromo = 0.0, tChromo = 0.0

  ! Input parametes for coronal inner BC's
  logical :: UseScaledWaveBcs = .false.
  real    :: WaveEnergyFactor = 0.0, WaveEnergyPower = 1./3.
  real    :: ExpFactorInvMin = 1e-3

  ! Input parameters controling wave dissipation
  logical :: UseWaveDissipation = .false.
  real    :: DissipationScaleFactorSi  ! unit = m*T^0.5
  real    :: DissipationScaleFactor
  real    :: LminIo , LmaxIo, WaveRatio

  ! variables for magnetic (unsigned flux) heating
  logical :: DoMagneticHeating = .false.
  real    :: TempLimitSi, HeatFactorCoeff = 1.0

  ! variables for Parker initial condition
  real    :: nCoronaSi = 0.0, tCoronaSi = 0.0

  ! Global arrays for plotting
  real       :: WaveDissip_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)      = 1.e-30
  real       :: WaveDissipPlus_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)  = 1.e-30
  real       :: WaveDissipMinus_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock) = 1.e-30
  real       :: DissipLength_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)    = 1.e-30
  integer    :: IsCounter_GB(-1:nI+2,-1:nJ+2,-1:nK+2,MaxBlock)       = 0
  
  ! Input parameters for two-temperature effects
  real :: TeFraction, TiFraction
  real :: QeByQtotal = 0.0
  real :: EtaPerpSi

contains 
 !============================================================================
  subroutine user_read_inputs

    use ModMain,      ONLY: UseUserInitSession, lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand
    character(len=*), parameter :: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    UseUserInitSession = .true.

    if(iProc == 0 .and. lVerbose > 0)then
       call write_prefix;
       write(iUnitOut,*)'User read_input CHROMOSPHERE-CORONA starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)

       case("#CHROMOBC")
          call read_var('UseChromoBc',UseChromoBc)
          if(UseChromoBc) then
             call read_var('WaveDeltaU', WaveDeltaU)
             call read_var('nChromoSi', nChromoSi)
             call read_var('tChromoSi', tChromoSi)
          end if

       case('#MAGHEATING')
          ! Heating related to unsigned magnetic flux, Abbet
          call read_var('DoMagneticHeating',DoMagneticHeating)
          call read_var('TempLimitSi',      TempLimitSi)
          call read_var('HeatFactorCoeff',  HeatFactorCoeff)

       case("#PARKERIC")
          call read_var('nCoronaSi', nCoronaSi)
          call read_var('tCoronaSi', tCoronaSi)

       case("#WAVEDISSIPATION")
          call read_var('UseWaveDissipation', UseWaveDissipation)
          if(UseWaveDissipation) then 
             call read_var('LminIo', LminIo)
             call read_var('LmaxIo', LmaxIo)
             call read_var('WaveRatio', WaveRatio)
          end if

       case("#WAVEBOUNDARY")
          call read_var('UseScaledWaveBcs', UseScaledWaveBcs)
          if (UseScaledWaveBcs) then
             if (UseChromoBc) call CON_stop('ERROR in '//NameSub//'- conflicting wave BCs')
             call read_var('WaveEnergyFactor', WaveEnergyFactor)
             call read_var('WaveEnergyPower',  WaveEnergyPower)
             call read_var('ExpFactorInvMin',  ExpFactorInvMin)
          end if
         
       case('#ELECTRONHEATING')
          ! Steven Cranmer his electron heating fraction (Apj 2009)
          call read_var('QeByQtotal', QeByQtotal)

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

    use ModMain,           ONLY: UseB0, UseMagnetogram
    use ModProcMH,         ONLY: iProc
    use ModReadParam,      ONLY: i_session_read
    use ModIO,             ONLY: write_prefix, iUnitOut
    use ModWaves,          ONLY: UseWavePressure, UseAlfvenWaves
    use ModAdvance,        ONLY: UseElectronPressure
    use ModCoronalHeating, ONLY: get_coronal_heat_factor, HeatFactor
    use ModMultiFluid,     ONLY: MassIon_I
    use ModConst,          ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
                                 cElectronMass
    use ModNumConst,       ONLY: cTwoPi
    use ModPhysics,        ONLY: ElectronTemperatureRatio, AverageIonCharge, &
                                 Si2No_V, UnitTemperature_, UnitN_

    real            :: HeatCondParSi
    real, parameter :: CoulombLog = 20.0
    character (len=*),parameter :: NameSub = 'uset_init_session'
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    UseAlfvenWaves = .true.
    UseWavePressure = .true.

    if(UseChromoBc) then
       ! convert to normalized units
       nChromo = nChromoSi*Si2No_V(UnitN_)
       RhoChromo = nChromo*MassIon_I(1)
       tChromo = tChromoSi*Si2No_V(UnitTemperature_)
    end if

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


    if(.not.UseB0)UseMagnetogram=.false.

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session
  !============================================================================
  subroutine user_initial_perturbation

    use ModCoronalHeating,  ONLY: get_coronal_heat_factor, HeatFactor
    use ModProcMH,          ONLY: iProc

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    ! ---------------------------------------------------------------
    if (DoMagneticHeating) then
       call get_coronal_heat_factor
       if (iProc == 0) write(*,*) 'Coronal magnetic heating set: ',HeatFactor
    end if

  end subroutine user_initial_perturbation
  !============================================================================
  subroutine get_plasma_parameters_base(x_D, RhoBase, Tbase)

    ! This subroutine computes the base values for mass density and temperature
    ! OPTION 1:
    ! Coronal boundary conditions: use LDEM moments (from EUV Tomography)
    !
    ! OPTION 2:
    ! Chromospheric boundary conditions - uniform T and Rho from PARAM.in

    use ModLdem,             ONLY: get_ldem_moments, UseLdem
    use ModPhysics,          ONLY: BodyRho_I, Si2No_V, UnitTemperature_, &
                                   UnitN_, AverageIonCharge
    use ModMultiFluid,       ONLY: MassIon_I

    real, intent(in) :: x_D(3)
    real, intent(out):: RhoBase, Tbase

    ! Electron density and temperature from ldem moments (cm^-3, MK respectively)
    real :: Ne, Te

    character (len=*), parameter :: NameSub = 'get_plasma_parameters_base'
    !--------------------------------------------------------------------------   
    if(UseLdem)then
       call get_ldem_moments(x_D, Ne, Te)
       RhoBase = Ne*Si2No_V(UnitN_)*MassIon_I(1)/AverageIonCharge
       Tbase = Te*Si2No_V(UnitTemperature_)

       RETURN
    end if

    if(UseChromoBc) then
       RhoBase = RhoChromo
       TBase   = tChromo
    end if

  end subroutine get_plasma_parameters_base
  !============================================================================
  subroutine get_wave_energy_base(Btot, Ewave)

    ! Provides the distribution of the total Alfven wave energy density
    ! at the lower boundary

    use ModPhysics,    ONLY: Si2No_V, UnitU_

    real, intent(in) :: Btot
    real, intent(out):: Ewave

    real :: RhoBase, Tbase

    character (len=*), parameter :: NameSub = 'get_wave_energy_base'
    !--------------------------------------------------------------------------
    Ewave = 0.0

    if (UseChromoBc) then
       
       Ewave = RhoChromo*(WaveDeltaU*Si2No_V(UnitU_))**2

    else if(UseScaledWaveBcs)then
      
       Ewave = WaveEnergyFactor*(Btot)**WaveEnergyPower

    end if

  end subroutine get_wave_energy_base
  !============================================================================
  subroutine user_set_ics

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, B0_DGB
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk, r_Blk, true_cell
    use ModMain,       ONLY: nI, nJ, nK, globalBLK, UseB0, unusedBLK
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         BodyRho_I, BodyP_I, BodyNDim_I, UnitU_, UnitN_, AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_, Pe_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    real :: x, y, z, r, Rho, NumDensIon, NumDensElectron
    real :: RhoBase, Tbase, Ubase, B_D(3), r_D(3), Br
    ! variables for iterative Parker solution
    integer :: IterCount
    real :: Ur, Ur0, Ur1, del, rTransonic, Uescape, Usound

    real, parameter :: Epsilon = 1.0e-6
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------

    iBlock = globalBLK

    ! Initially, density, electron and ion temperature are at coronal
    ! values starting from just above the boundary
    RhoBase = nCoronaSi*Si2No_V(UnitN_)*MassIon_I(1)
    Tbase   = tCoronaSi*Si2No_V(UnitTemperature_)

    ! normalize with isothermal sound speed.
    Usound = sqrt(Tbase*(1.0+AverageIonCharge)/MassIon_I(1))
    Uescape = sqrt(-GBody*2.0)/Usound

    !\
    ! Initialize MHD wind with Parker's solution
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    rTransonic = 0.25*Uescape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    Ubase = rTransonic**2*exp(1.5 - 2.0*rTransonic)

    do k = 1, nK ; do j = 1, nJ ; do i = 1, nI
       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)
       z = z_BLK(i,j,k,iBlock)
       r = r_BLK(i,j,k,iBlock)
       r_D = (/x,y,z/)

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
                call CON_stop('PARKER > 1000 it.')
             end if
          end do
       end if

       ! Set chromospheric initial condition inside the body, if needed
       if(UseChromoBc .and. r <= rBody)then

          State_VGB(rho_,i,j,k,iBlock) = RhoChromo
          State_VGB(p_,  i,j,k,iBlock) = 2*tChromo*nChromoSi*Si2No_V(UnitN_)
       else

          Rho = rBody**2*RhoBase*Ubase/(r**2*Ur)
          State_VGB(Rho_,i,j,k,iBlock) = Rho

          NumDensIon = Rho/MassIon_I(1)
          NumDensElectron = NumDensIon*AverageIonCharge
          
          if(UseElectronPressure)then
             State_VGB(p_,i,j,k,iBlock) = NumDensIon*Tbase
             State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*Tbase
          else
             State_VGB(p_,i,j,k,iBlock) = &
                  (NumDensIon + NumDensElectron)*Tbase
          end if
       end if

       State_VGB(RhoUx_,i,j,k,iBlock) = Rho*Ur*x/r *Usound
       State_VGB(RhoUy_,i,j,k,iBlock) = Rho*Ur*y/r *Usound
       State_VGB(RhoUz_,i,j,k,iBlock) = Rho*Ur*z/r *Usound

       if(UseB0)then
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
       else
          call get_coronal_b0(x, y, z, B_D)
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = B_D
       end if

       if (UseChromoBc) then
          Br = sum(B0_DGB(1:3,i,j,k,iBlock)*r_D)
          if (Br >= 0.0) then
             State_VGB(WaveFirst_,i,j,k,iBlock) =  &
                  State_VGB(rho_,i,j,k,iBlock)*(WaveDeltaU*Si2No_V(UnitU_))**2
             State_VGB(WaveLast_,i,j,k,iBlock) = 1e-12
          else
             State_VGB(WaveLast_,i,j,k,iBlock) =  &
                  State_VGB(rho_,i,j,k,iBlock)*(WaveDeltaU*Si2No_V(UnitU_))**2
             State_VGB(WaveFirst_,i,j,k,iBlock) = 1e-12
          end if
       else
          State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 1e-12
       end if

    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
  subroutine user_update_states(iStage, iBlock)

    integer, intent(in) :: iStage, iBlock
    !--------------------------------------------------------------------------
    call update_states_MHD(iStage, iBlock)

  end subroutine user_update_states
  !============================================================================
  subroutine user_calc_sources

    use ModAdvance,        ONLY: State_VGB, Source_VC, UseElectronPressure, &
                                 B0_DGB
    use ModGeometry,       ONLY: r_BLK
    use ModMain,           ONLY: nI, nJ, nK, GlobalBlk, UseB0
    use ModPhysics,        ONLY: gm1, rBody, BodyRho_I, Si2No_V, UnitX_,&
                                 UnitB_,  UnitTemperature_, No2Si_V
    use ModVarIndexes,     ONLY: Rho_, Bx_, Bz_, Energy_, p_, Pe_, &
                                 WaveFirst_, WaveLast_
    use ModMultifluid,     ONLY: MassIon_I
    use ModCoronalHeating, ONLY: HeatFactor

    integer :: i, j, k, iBlock, iWave
    real    :: WaveEnergyPlus, WaveEnergyMinus, TemperatureSi, FullB_D(3), FullB
    real    :: CoronalHeating, RadiativeCooling, MagneticHeating

    ! varaibles for wave dissipation
    real    :: LocalDissipationFactor, Lmin, Lmax
    real    :: WaveDissipationPlus, WaveDissipationMinus, CounterWaveDissipRate

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
          FullB = sqrt(sum(FullB_D**2))

          WaveEnergyPlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
          WaveEnergyMinus = State_VGB(WaveLast_,i,j,k,iBlock)

          !\                                                             
          ! Calculate Local Dissipation factor, dimensions: length/sqrt(density/B)  
          !/ 
          ! Local dissipation factor: 
          ! Assumed to vary depending on density and magnetic field

          Lmin = LminIo*Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
          Lmax = LmaxIo*Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
          LocalDissipationFactor = 1./(Lmax * &
               sqrt(State_VGB(Rho_,i,j,k,iBlock)/FullB))

          !\                                                                        
          ! Calculate coronal heating due to wave dissipation
          !/                                                                                       
          CoronalHeating =  0.0

          if(WaveEnergyPlus <= WaveEnergyMinus/WaveRatio .or. &
             WaveEnergyMinus <= WaveEnergyPlus/WaveRatio )  then
             ! Dissipate waves in "coronal hole"
             WaveDissipationPlus = LocalDissipationFactor * &
                  WaveEnergyPlus**(3./2.)
             WaveDissipationMinus = LocalDissipationFactor * &
                  WaveEnergyMinus**(3./2.)
             DissipLength_GB(i,j,k,iBlock) = Lmax/sqrt(FullB)
             IsCounter_GB(i,j,k,iBlock) = 0
          else
             ! Dissipate counter propagating waves
             CounterWaveDissipRate = LocalDissipationFactor *(Lmax/Lmin)* &
                  sqrt(2.*WaveEnergyPlus*WaveEnergyMinus/&
                  (WaveEnergyPlus + WaveEnergyMinus))

             WaveDissipationPlus  = CounterWaveDissipRate* WaveEnergyPlus
             WaveDissipationMinus = CounterWaveDissipRate* WaveEnergyMinus
             DissipLength_GB(i,j,k,iBlock) = Lmin/sqrt(FullB)
             IsCounter_GB(i,j,k,iBlock) = 1

          end if

          CoronalHeating = WaveDissipationPlus + WaveDissipationMinus

          ! Remove wave energy from state variables
          Source_VC(WaveFirst_,i,j,k) = Source_VC(WaveFirst_,i,j,k) &
               - WaveDissipationPlus
          Source_VC(WaveLast_,i,j,k) = Source_VC(WaveLast_,i,j,k) &
               - WaveDissipationMinus

          ! save for plotting
          WaveDissip_GB(i,j,k,iBlock) = CoronalHeating
          WaveDissipPlus_GB(i,j,k,iBlock)  = WaveDissipationPlus
          WaveDissipMinus_GB(i,j,k,iBlock) = WaveDissipationMinus
       end if
      
       ! Add magnetic heating if desired (due to Abbett 2007)
       if (DoMagneticHeating) then
          TemperatureSi = No2Si_V(UnitTemperature_)* MassIon_I(1)* &
               State_VGB(p_,i,j,k,iBlock)/State_VGB(rho_,i,j,k,iBlock)
          MagneticHeating = 1e-30
          if (TemperatureSi <= TempLimitSi .and. &
              r_BLK(i,j,k,iBlock) < 2.*rBody) then

             MagneticHeating = HeatFactorCoeff*HeatFactor*FullB
             CoronalHeating = CoronalHeating + MagneticHeating
             ! for plotting
             ! MagHeating_GB(i,j,k,iBlock) = MagneticHeating
          end if
       end if

       if(UseElectronPressure)then
          Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
               + gm1*(1.0 - QeByQtotal)*CoronalHeating
          Source_VC(Pe_,i,j,k) = Source_VC(Pe_,i,j,k) &
               + gm1*QeByQtotal*CoronalHeating
       else
          Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) + gm1*CoronalHeating
       end if

    end do; end do; end do

  end subroutine user_calc_sources
  !============================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, B0_DGB, UseElectronPressure
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK, nBlock, x_, y_, z_, UseB0
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
          if(UseB0)then
             tmp1_BLK(:,:,:,iBlock) = & 
                  ( B0_DGB(x_,:,:,:,iBlock) + State_VGB(Bx_,:,:,:,iBlock))**2 &
                  +(B0_DGB(y_,:,:,:,iBlock) + State_VGB(By_,:,:,:,iBlock))**2 &
                  +(B0_DGB(z_,:,:,:,iBlock) + State_VGB(Bz_,:,:,:,iBlock))**2
          else
             tmp1_BLK(:,:,:,iBlock) = State_VGB(Bx_,:,:,:,iBlock)**2 &
                  + State_VGB(By_,:,:,:,iBlock)**2 &
                  + State_VGB(Bz_,:,:,:,iBlock)**2
          end if
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

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, B0_DGB
    use ModMain,       ONLY: UseB0
    use ModNumConst,   ONLY: cTolerance
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, UnitEnergyDens_, UnitX_
    use ModSize,       ONLY: nI, nJ, nK
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Bx_, Bz_, WaveFirst_, WaveLast_
    use ModRadiativeCooling, ONLY: RadCooling_CB

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
    real :: FullB, FullB2, p, Ewave

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    IsFound = .true.

    select case(NameVar)

    case('disstot')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          PlotVar_G(i,j,k) = max(1e-30,WaveDissip_GB(i,j,k,iBlock) * &
               No2Si_V(UnitEnergyDens_))

       end do; end do ; end do


    case('dissplus')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          PlotVar_G(i,j,k) = max(1e-30,WaveDissipPlus_GB(i,j,k,iBlock) * &
               No2Si_V(UnitEnergyDens_))

       end do; end do ; end do


    case('dissminus')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          PlotVar_G(i,j,k) = max(1e-30,WaveDissipMinus_GB(i,j,k,iBlock) * &
               No2Si_V(UnitEnergyDens_))

       end do; end do ; end do


    case('disslen')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          PlotVar_G(i,j,k) = max(1e-30,DissipLength_GB(i,j,k,iBlock) * &
               No2Si_V(UnitX_))

       end do; end do ; end do

    case('iscounter')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          PlotVar_G(i,j,k) = IsCounter_GB(i,j,k,iBlock)

       end do; end do ; end do

    case('radcool')

       do k = 1, nK; do j = 1, nJ; do i = 1, nI

          PlotVar_G(i,j,k) = max(1e-30,-RadCooling_CB(i,j,k,iBlock) * &
               No2Si_V(UnitEnergyDens_))

       end do; end do ; end do

    case('magheat')

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2

          !PlotVar_G(i,j,k) = max(1e-30,MagHeating_GB(i,j,k,iBlock) * &
          !     No2Si_V(UnitEnergyDens_))

       end do; end do ; end do

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

    case('deltabperb')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          if(UseB0)then
             FullB = sqrt(sum((B0_DGB(:,i,j,k,iBlock) &
                  + State_VGB(Bx_:Bz_,i,j,k,iBlock))**2))
          else
             FullB = sqrt(sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2))
          end if
          Ewave = sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
          PlotVar_G(i,j,k) = sqrt(Ewave)/max(FullB, cTolerance)
       end do; end do; end do

    case('beta')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          if(UseB0)then
             FullB2 = sum((B0_DGB(:,i,j,k,iBlock) &
                  + State_VGB(Bx_:Bz_,i,j,k,iBlock))**2)
          else
             FullB2 = sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)
          end if
          if(UseElectronPressure)then
             p = State_VGB(p_,i,j,k,iBlock) + State_VGB(Pe_,i,j,k,iBlock)
          else
             p = State_VGB(p_,i,j,k,iBlock)
          end if
          PlotVar_G(i,j,k) = 2.0*p/max(FullB2, cTolerance)
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
    use ModMain,     ONLY: x_, y_, z_, time_loop, UseB0

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
    if(UseB0)then
       do k=0, nK+1; do j=1, nJ; do i=1, nI
          rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
               * (B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock)) &
               +           y_BLK(i,j,k,iBlock)   &
               * (B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock)) &
               +           z_BLK(i,j,k,iBlock)   &
               * (B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))
       end do; end do; end do
    else
       do k=0, nK+1; do j=1, nJ; do i=1, nI
          rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
               * State_VGB(Bx_,i,j,k,iBlock) &
               +           y_BLK(i,j,k,iBlock)   &
               * State_VGB(By_,i,j,k,iBlock) &
               +           z_BLK(i,j,k,iBlock)   &
               * State_VGB(Bz_,i,j,k,iBlock)
       end do; end do; end do
    end if

    DoRefine = maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny

  end subroutine user_specify_refinement
  !============================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    ! Fill one layer of ghost cells with the temperature for heat conduction

    use ModAdvance,    ONLY: State_VGB
    use ModMain,       ONLY: East_
    use ModVarIndexes, ONLY: Rho_, p_
 
    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------

    if(iSide /= East_) call CON_stop('Wrong iSide in user_set_outerBCs')

    IsFound = .true.

    if (UseChromoBc) then

       State_VGB(Rho_,-1:0,:,:,iBlock) = RhoChromo
       State_VGB(p_,-1:0,:,:,iBlock) = 2.*nChromo*tChromo
       
    else
       call CON_stop('User boundary conditions are not specified')
    end if

  end subroutine user_set_outerbcs
  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance,     ONLY: State_VGB, UseElectronPressure
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame, GlobalBLK, UseB0, &
         UseHeatConduction
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: OmegaBody, AverageIonCharge, BodyRho_I, &
                              BodyTDim_I, Si2No_V, UnitTemperature_, UnitN_
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_, &
         WaveFirst_, WaveLast_, Pe_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: RhoBase, NumDensIon, NumDensElectron, Tbase, FullBr, B0tot
    real :: Runit_D(3), U_D(3)
    real :: B1_D(3), B1t_D(3), B1r_D(3), FullB_D(3), B0_D(3)
    real :: Ewave

    character (len=*), parameter :: NameSub = 'user_face_bcs'
    !--------------------------------------------------------------------------

    Runit_D = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    U_D   = VarsTrueFace_V(Ux_:Uz_)
    VarsGhostFace_V(Ux_:Uz_) = -U_D

    if(UseB0)then
       B1_D  = VarsTrueFace_V(Bx_:Bz_)
       B1r_D = sum(Runit_D*B1_D)*Runit_D
       B1t_D = B1_D - B1r_D
       VarsGhostFace_V(Bx_:Bz_) = B1t_D !- B1r_D
       FullB_D = B0Face_D + B1t_D
       B0tot = sqrt(sum(B0Face_D**2))
    else
       call get_coronal_b0(FaceCoords_D(x_), FaceCoords_D(y_), &
            FaceCoords_D(z_), B0_D)
       B1_D  = VarsTrueFace_V(Bx_:Bz_) - B0_D
       B1r_D = sum(Runit_D*B1_D)*Runit_D
       B1t_D = B1_D - B1r_D
       VarsGhostFace_V(Bx_:Bz_) = B1t_D + B0_D
       FullB_D = VarsGhostFace_V(Bx_:Bz_)
       B0tot = sqrt(sum(B0_D**2))
 
    end if
    FullBr = sum(Runit_D*FullB_D)

    ! Set Alfven waves energy density
    call get_wave_energy_base(B0tot, Ewave)

    if(FullBr > 0.0)then
       VarsGhostFace_V(WaveFirst_) = Ewave
       VarsGhostFace_V(WaveLast_) = 1.e-12
    else
       VarsGhostFace_V(WaveFirst_) = 1e-12 !VarsTrueFace_V(WaveFirst_)
       VarsGhostFace_V(WaveLast_) = Ewave
    end if

    call get_plasma_parameters_base(FaceCoords_D, RhoBase, Tbase)
   
    VarsGhostFace_V(Rho_) =  2.0*RhoBase - VarsTrueFace_V(Rho_)
    NumDensIon = VarsGhostFace_V(Rho_)/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    if(UseElectronPressure)then
       VarsGhostFace_V(p_) = max(NumDensIon*Tbase, VarsTrueFace_V(p_))
       VarsGhostFace_V(Pe_) = NumDensElectron*Tbase
    else
       if(UseHeatConduction)then
          VarsGhostFace_V(p_) = (NumDensIon + NumDensElectron)*Tbase
       else
          VarsGhostFace_V(p_) = &
               max((NumDensIon + NumDensElectron)*Tbase, VarsTrueFace_V(p_))
       end if
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
