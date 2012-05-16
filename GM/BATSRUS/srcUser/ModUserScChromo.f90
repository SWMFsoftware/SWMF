!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain, ONLY: nI, nJ,nK

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
       IMPLEMENTED11=> user_set_resistivity,            &
       IMPLEMENTED12=> user_update_states,              &
       IMPLEMENTED13=> user_initial_perturbation

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Chromosphere to solar wind model with Alfven waves - Oran, van der Holst'

  ! Input parameters for chromospheric inner BC's
  logical :: UseChromoBc = .true., UseUparBc = .false., &
             UseExtrapolatedEwave = .true.
  real    :: WaveDeltaU = 0.0
  real    :: nChromoSi = 0.0, tChromoSi = 0.0
  real    :: nChromo = 0.0, RhoChromo = 0.0, tChromo = 0.0

  ! Input parametes for coronal inner BC's
  logical :: UseScaledWaveBcs = .false.
  real    :: WaveEnergyFactor = 0.0, WaveEnergyPower = 1./3.
  real    :: ExpFactorInvMin = 1e-3

  ! Input parameters controling wave dissipation
  logical :: UseWaveDissipation = .false.
  real    :: DissipationScaleFactorSi = 0.0  ! unit = m*T^0.5
  real    :: DissipationScaleFactor = 0.0
  real    :: LkolIo = 0.0, Lkol, Lcp, Lratio = 0.0

  ! variables for magnetic (unsigned flux) heating
  logical :: DoMagneticHeating = .false.
  real    :: TempLimitSi = 0.0, HeatFactorCoeff = 1.0

  ! variables for Parker initial condition
  real    :: nCoronaSi = 0.0, tCoronaSi = 0.0
  real    :: DipoleTiltDeg = 0.0

  ! Current sheet refinement 
  real,allocatable    :: rRefine_I(:)
  logical :: DoRefineGradualSheet = .false.
  integer :: nLevelGradual = 0

  ! Input parameters for two-temperature effects
  real    :: TeFraction, TiFraction
  real    :: QeByQtotal = 0.0
  real    :: EtaPerpSi

contains 
 !============================================================================
  subroutine user_read_inputs

    use ModMain,      ONLY: UseUserInitSession, lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModPhysics,   ONLY: DipoleStrengthSi, rBody

    integer :: iLevel
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
          
       ! This commans is used when the inner boundary is the chromosphere   
       case("#CHROMOBC")
          call read_var('UseChromoBc',UseChromoBc)
          if(UseChromoBc) then
             call read_var('WaveDeltaU', WaveDeltaU)
             call read_var('nChromoSi', nChromoSi)
             call read_var('tChromoSi', tChromoSi)
             call read_var('UseUparBc',UseUparBc)
             call read_var('UseExtrapolatedEwave',UseExtrapolatedEwave)
          end if

       ! This command is used when the inner boundary is the coronal base
       ! Replaces WSA - based wave boundary conditions
       case("#CORONAWAVEBC")
          call read_var('UseScaledWaveBcs', UseScaledWaveBcs)
          if (UseScaledWaveBcs) then
             if (UseChromoBc) call CON_stop('ERROR in '//NameSub//'- conflicting wave BCs')
             call read_var('WaveEnergyFactor', WaveEnergyFactor)
             call read_var('WaveEnergyPower',  WaveEnergyPower)
             call read_var('ExpFactorInvMin',  ExpFactorInvMin)
          end if

       case('#SOLARDIPOLE')
          call read_var('DipoleStrengthSi',DipoleStrengthSi)
          call read_var('DipoleTiltDeg',DipoleTiltDeg)

       case("#WAVEDISSIPATION")
          call read_var('UseWaveDissipation', UseWaveDissipation)
          if(UseWaveDissipation) then 
             call read_var('LkolIo', LkolIo)
             call read_var('Lratio', Lratio)
          end if

       case('#MAGHEATING')
          ! Heating related to unsigned magnetic flux, Abbet
          call read_var('DoMagneticHeating',DoMagneticHeating)
          call read_var('TempLimitSi',      TempLimitSi)
          call read_var('HeatFactorCoeff',  HeatFactorCoeff)

       case("#PARKERIC")
          call read_var('nCoronaSi', nCoronaSi)
          call read_var('tCoronaSi', tCoronaSi)

       case("#REFINEGRADUALSHEET")
          call read_var('DoRefineGradualSheet', DoRefineGradualSheet)
          call read_var('nLevelGradual'     , nLevelGradual)
          if (DoRefineGradualSheet .and. nLevelGradual > 0) then
             allocate(rRefine_I(0:nLevelGradual))
             rRefine_I(0) = rBody
             do iLevel = 1, nLevelGradual
                call read_var('rRefine', rRefine_I(iLevel))
             end do
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

    use ModMain,           ONLY: UseB0, UseMagnetogram, UseUserPerturbation
    use ModProcMH,         ONLY: iProc
    use ModReadParam,      ONLY: i_session_read
    use ModIO,             ONLY: write_prefix, iUnitOut
    use ModWaves,          ONLY: UseWavePressure, UseAlfvenWaves
    use ModAdvance,        ONLY: UseElectronPressure
    use ModMultiFluid,     ONLY: MassIon_I
    use ModConst,          ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
                                 cElectronMass
    use ModNumConst,       ONLY: cTwoPi, cDegToRad
    use ModPhysics,        ONLY: ElectronTemperatureRatio, AverageIonCharge, &
                                 Si2No_V, UnitTemperature_, UnitN_, UnitX_, UnitB_, &
                                 SinThetaTilt, CosThetaTilt

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

    if(UseWaveDissipation)then
       ! Reference dissipation length: 
       ! Lcp - used for counter propogating wave dissipation
       ! Lkol - used for Kolmogorov dissipation
       Lkol = LkolIo*Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
       if (Lratio > 0.0)  Lcp = Lkol/Lratio
    end if

    if (DoMagneticHeating .and. .not. UseUserPerturbation) &
         call CON_stop('UseUserInitPerturbation must be set to T'// &
         ' when using magnetic heating. Correct PARAM.in file')

    if (.not. UseMagnetogram) then
       SinThetaTilt = sin(cDegToRad*DipoleTiltDeg)
       CosThetaTilt = cos(cDegToRad*DipoleTiltDeg)
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
    use ModPhysics,          ONLY: Si2No_V, UnitTemperature_, &
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
  subroutine get_wave_energy_base(FullB, Ewave)

    ! Provides the distribution of the total Alfven wave energy density
    ! at the lower boundary

    use ModPhysics,    ONLY: Si2No_V, UnitU_

    real, intent(in) :: FullB
    real, intent(out):: Ewave

    character (len=*), parameter :: NameSub = 'get_wave_energy_base'
    !--------------------------------------------------------------------------
    Ewave = 0.0

    if (UseChromoBc) then
       
       Ewave = RhoChromo*(WaveDeltaU*Si2No_V(UnitU_))**2

    else if(UseScaledWaveBcs)then
      
       Ewave = WaveEnergyFactor*(FullB)**WaveEnergyPower

    end if

  end subroutine get_wave_energy_base
  !============================================================================
  subroutine user_set_ics

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, B0_DGB
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk, r_Blk, true_cell
    use ModMain,       ONLY: globalBLK, UseB0, unusedBLK
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         BodyRho_I, BodyP_I, BodyNDim_I, UnitU_, UnitN_, AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_, Pe_, &
         WaveFirst_, WaveLast_, Hyp_

    integer :: i, j, k, iBlock
    real :: x, y, z, r, Rho, NumDensIon, NumDensElectron, Temperature
    real :: RhoCorona, tCorona, uCorona, rCorona, TemperatureGradient
    real :: B_D(3), r_D(3), Br
    ! variables for iterative Parker solution
    integer :: IterCount
    real :: Ur, Ur0, Ur1, del, rTransonic, Uescape, Usound

    real, parameter :: Epsilon = 1.0e-6
    character (len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    ! Initially, density, electron and ion temperature are at coronal
    ! values starting from just above the boundary
    RhoCorona = nCoronaSi*Si2No_V(UnitN_)*MassIon_I(1)
    tCorona   = tCoronaSi*Si2No_V(UnitTemperature_)

    ! normalize with isothermal sound speed.
    Usound = sqrt(tCorona*(1.0+AverageIonCharge)/MassIon_I(1))
    Uescape = sqrt(-GBody*2.0)/Usound

    !\
    ! Initialize MHD wind with Parker's solution
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    rTransonic = 0.25*Uescape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    uCorona = rTransonic**2*exp(1.5 - 2.0*rTransonic)

    do k = -1, nK+2 ; do j = -1, nJ+2 ; do i = -1, nI+2
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
       if(UseChromoBc) then
          if ( r <= rBody)then
             Rho = RhoChromo
             Temperature = tChromo
          else
             ! Density jumps to coronal values right outside body,
             ! then set according to constant mass flux
             Rho = rBody**2*RhoCorona*uCorona/(r**2*Ur)
          
             ! Temperature increases linearily from chromospheric 
             ! to coronal values at r=rCorona (for stability)
             rCorona = 1.1
             if (r < rCorona) then
                TemperatureGradient = (tCorona - tChromo)/(rCorona - rBody)    
                Temperature = &
                     (tChromo + TemperatureGradient*(r - rBody))
             else
                Temperature = tCorona
             end if
          end if
       else
          Rho = rBody**2*RhoCorona*uCorona/(r**2*Ur)
          Temperature = tCorona
       end if
      
       NumDensIon = Rho/MassIon_I(1)
       NumDensElectron = NumDensIon*AverageIonCharge

       if(UseElectronPressure)then
          State_VGB(p_,i,j,k,iBlock) = NumDensIon*Temperature
          State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*Temperature
       else
          State_VGB(p_,i,j,k,iBlock) = &
               (NumDensIon + NumDensElectron)*Temperature
       end if
       if (r > 1.5) then
          ! NOTE: If you wish to recover the Parker solution, remove this case.

          ! "Vacuum cleaner" initial condition to speed up convergence.
          ! Basically, it makes the initial solution lighter and easier to
          ! push out of the domain. Note, this initial solution is not
          ! self consistent, but this should not affect the final state.
          State_VGB(Rho_,i,j,k,iBlock) = Rho/2
       else
          State_VGB(Rho_,i,j,k,iBlock) = Rho
       end if

       State_VGB(RhoUx_,i,j,k,iBlock) = Rho*Ur*x/r *Usound
       State_VGB(RhoUy_,i,j,k,iBlock) = Rho*Ur*y/r *Usound
       State_VGB(RhoUz_,i,j,k,iBlock) = Rho*Ur*z/r *Usound

       if(UseB0)then
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
          Br = sum(B0_DGB(1:3,i,j,k,iBlock)*r_D)
       else
          call get_coronal_b0(x, y, z, B_D)
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = B_D
 	  Br = sum(B_D*r_D)
       end if
       if(Hyp_>1) State_VGB(Hyp_,i,j,k,iBlock) = 0.0

       if (UseChromoBc) then
          if (Br >= 0.0) then
             State_VGB(WaveFirst_,i,j,k,iBlock) =  &
                  State_VGB(rho_,i,j,k,iBlock)*(WaveDeltaU*Si2No_V(UnitU_))**2
             State_VGB(WaveLast_,i,j,k,iBlock) = 1e-30
          else
             State_VGB(WaveLast_,i,j,k,iBlock) =  &
                  State_VGB(rho_,i,j,k,iBlock)*(WaveDeltaU*Si2No_V(UnitU_))**2
             State_VGB(WaveFirst_,i,j,k,iBlock) = 1e-30
          end if
       else
          State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 1e-30
       end if

    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
  subroutine user_update_states(iStage, iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModVarIndexes, ONLY: nVar

    integer, intent(in) :: iStage, iBlock

    !--------------------------------------------------------------------------
    call update_states_MHD(iStage, iBlock)
    
  end subroutine user_update_states
  !============================================================================
  subroutine write_ghost_and_boundary

    ! output ghost cells and true cells values for blocks that cross the
    ! inner boundary, in a meridional plane close to the x=0 plane

    ! ------ SPHERICAL GRID ONLY ------------

    use ModProcMH,     ONLY: iProc
    use ModIoUnit,     ONLY: UNITTMP_
    use ModVarIndexes, ONLY: WaveFirst_, WaveLast_, rho_, p_, Bx_, Bz_
    use ModMain,       ONLY: iteration_number, unusedBLK, nBlock
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModAdvance,    ONLY: State_VGB

    integer            :: iBLK, i ,j ,k
    real               :: x, y, z, r, x_D(3)
    real               :: IPlus, IMinus, B, Br, rho, p
    character(len=40)  :: FileNameTec 
    character(len=11)  :: NameStage
    character(len=7)   :: NameProc

    character(len=*),parameter   :: NameSub='write_ghost_and_boundary'
    !-------------------------------------------------------------------
    write(NameStage,'(i5.5)') iteration_number
    write(NameProc, '(i4.4)') iProc
    FileNameTec = &
         'SC/IO2/Boundary_n_'//trim(NameStage)//'_p'//trim(NameProc)//'.txt'
    open(UNITTMP_, file=FileNameTec, form='formatted', access='sequential', &
         status = 'replace')
    write(*,*) 'Writing ', FileNameTec
      
    ! choose a single plane (more conditions inside loop)
    j = 4
    k = 1
    do iBLK = 1, nBlock
       ! Only choose blocks that cross the boundary
       if ( minval(r_BLK(:,:,:,iBLK)) <= 1. .and. &
            maxval(r_BLK(:,:,:,iBLK)) > 1.) then

          ! Advance in the r direction
          do i=-1,nI+2

             x=x_BLK(i,j,k,iBLK)
             ! Extract data in the x=0 plane only
             if (abs(X) > 0.1) CYCLE
             
             y = y_BLK(i,j,k,iBLK)
             z = z_BLK(i,j,k,iBLK)
             r = r_BLK(i,j,k,iBLK)
             x_D = (/x,y,z/)
             
             Iplus = State_VGB(WaveFirst_,i,j,k,iBLK)
             Iminus = State_VGB(WaveLast_,i,j,k,iBLK)
             B = sum(State_VGB(Bx_:Bz_,i,j,k,iBLK)**2)
             Br = sum(x_D*State_VGB(Bx_:Bz_,i,j,k,iBLK))
             rho = State_VGB(rho_,i,j,k,iBLK)
             p = State_VGB(p_,i,j,k,iBLK)

             ! write to file
             write(UNITTMP_, '(17e12.5)') real(i),real(j),real(k),&
                  real(iBLK), real(iProc),&
                  x, y, z, r, Iplus, Iminus, B, Br,  rho, p
          end do

       end if
    end do
    close(UNITTMP_)

  end subroutine write_ghost_and_boundary
  !============================================================================
  subroutine write_b0

    ! output b0 components in all cells

    use ModProcMH,     ONLY: iProc
    use ModIoUnit,     ONLY: UNITTMP_
    use ModMain,       ONLY: iteration_number, unusedBLK, nBlock
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModAdvance,    ONLY: B0_DGB

    integer            :: iBlock, i ,j ,k
    real               :: x, y, z, r
    real               :: B0_x, B0_y, B0_z
    character(len=40)  :: FileNameTec 
    character(len=11)  :: NameStage
    character(len=7)   :: NameProc

    character(len=*),parameter   :: NameSub='write_b0'
    !-------------------------------------------------------------------
    write(NameStage,'(i5.5)') iteration_number
    write(NameProc, '(i4.4)') iProc
    FileNameTec = &
         'SC/IO2/B0_n_'//trim(NameStage)//'_p'//trim(NameProc)//'.txt'
    open(UNITTMP_, file=FileNameTec, form='formatted', access='sequential', &
         status = 'replace')
    write(*,*) 'Writing ', FileNameTec
      
    ! choose a single plane (more conditions inside loop)
    j = 4
    do iBlock = 1, nBlock
       ! Advance in the r direction
       do k=-1,nK+2 ; do i=-1,nI+2

          x=x_BLK(i,j,k,iBlock)
          ! Extract data in the x=0 plane only
          if (abs(X) > 0.1) CYCLE
             
          y = y_BLK(i,j,k,iBlock)
          z = z_BLK(i,j,k,iBlock)
          r = r_BLK(i,j,k,iBlock)
             
          B0_x = B0_DGB(1,i,j,k,iBlock)
          B0_y = B0_DGB(2,i,j,k,iBlock)
          B0_z = B0_DGB(3,i,j,k,iBlock)

          ! write to file
          write(UNITTMP_, '(17e12.5)') real(i),real(j),real(k),&
               real(iBlock), real(iProc),&
               x, y, z, r, B0_x, B0_y, B0_z
       end do; end do
    end do
    close(UNITTMP_)

  end subroutine write_b0
  !============================================================================
  subroutine user_calc_sources

    use ModAdvance,        ONLY: State_VGB, Source_VC, UseElectronPressure, &
         B0_DGB, VdtFace_x, VdtFace_y, VdtFace_z
    use ModGeometry,       ONLY: r_BLK, vInv_CB
    use ModMain,           ONLY: GlobalBlk
    use ModPhysics,        ONLY: gm1, rBody, UnitTemperature_, No2Si_V
    use ModVarIndexes,     ONLY: Rho_, Energy_, p_, Pe_, WaveFirst_, WaveLast_
    use ModMultifluid,     ONLY: MassIon_I
    use ModCoronalHeating, ONLY: HeatFactor

    integer :: i, j, k, iBlock
    real    :: TemperatureSi, FullB
    real    :: CoronalHeating, RadiativeCooling, MagneticHeating
    real    :: DtInvWavePlus, DtInvWaveMinus, Vdt_Source, Vdt

    ! varaibles for wave dissipation
    real    :: WaveDissipationPlus, WaveDissipationMinus

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !------------------------------------------------------------------------- 
    iBlock = globalBlk

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(r_BLK(i,j,k,iBlock) < rBody) CYCLE

       !\
       ! Calculate coronal heating due to wave dissipation
       !/         
       if(UseWaveDissipation)then
          CoronalHeating =  0.0
          call calc_alfven_wave_dissipation(i, j, k, iBlock, &
               WaveDissipationPlus,  &
               WaveDissipationMinus) 

          CoronalHeating = WaveDissipationPlus + WaveDissipationMinus

          ! Remove wave energy from state variables
          Source_VC(WaveFirst_,i,j,k) = Source_VC(WaveFirst_,i,j,k) &
               - WaveDissipationPlus
          Source_VC(WaveLast_,i,j,k) = Source_VC(WaveLast_,i,j,k) &
               - WaveDissipationMinus
       end if
      

       ! Limit time step according to size of source terms:
       DtInvWavePlus = WaveDissipationPlus/ &
            max(State_VGB(WaveFirst_,i,j,k,iBlock),1e-30)
       DtInvWaveMinus = WaveDissipationMinus/ &
            max(State_VGB(WaveLast_,i,j,k,iBlock),1e-30)
       
       Vdt_Source = 2.0*max(DtInvWavePlus, DtInvWaveMinus)/&
            vInv_CB(i,j,k,iBlock)

       Vdt = min(VdtFace_x(i,j,k),VdtFace_y(i,j,k),VdtFace_z(i,j,k))

       if(Vdt_Source > Vdt)then
          ! The following prevents the wave energy from becoming negative
          ! due to too large source terms. This should be cell-centered,
          ! since the source are, but we add
          ! them for now to the left face of VdtFace.
          VdtFace_x(i,j,k) = max(VdtFace_x(i,j,k), Vdt_Source)
          VdtFace_y(i,j,k) = max(VdtFace_y(i,j,k), Vdt_Source)
          VdtFace_z(i,j,k) = max(VdtFace_z(i,j,k), Vdt_Source)
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
  subroutine calc_alfven_wave_dissipation(i,j,k,iBlock, &
       WaveDissipationPlus, WaveDissipationMinus, &
       DissipationLengthKol, DissipationLengthCp )

    use ModAdvance,        ONLY: State_VGB, B0_DGB
    use ModMain,           ONLY: UseB0
    use ModPhysics,        ONLY: Si2No_V, No2Si_V, UnitX_, UnitB_
    use ModVarIndexes,     ONLY: Rho_, Bx_, Bz_, WaveFirst_, WaveLast_

    integer,intent(in)        :: i, j, k, iBlock
    real,intent(out)          :: WaveDissipationPlus, WaveDissipationMinus
    real,intent(out),optional :: DissipationLengthKol, DissipationLengthCp
    
    real    :: WaveEnergyPlus, WaveEnergyMinus, FullB_D(3), FullB
    real    :: WavePressure, LocalDissipationFactor

    character (len=*), parameter :: NameSub = 'calc_alfven_wave_dissipation'
    !------------------------------------------------------------------------- 
    if(UseB0)then
       FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    FullB = sqrt(sum(FullB_D**2))
    
    ! Ouput dissipation lengths in units of length for plotting
    if(present(DissipationLengthKol)) &
         DissipationLengthKol = Lkol/sqrt(FullB)
    if(present(DissipationLengthCp)) &
         DissipationLengthCp = DissipationLengthKol/Lratio

    ! Local Dissipation factor, dimensions length/sqrt(density/B)  
    LocalDissipationFactor = &
         sqrt(FullB/State_VGB(Rho_,i,j,k,iBlock))/Lkol

    WaveEnergyPlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
    WaveEnergyMinus = State_VGB(WaveLast_,i,j,k,iBlock)
    WavePressure = (WaveEnergyPlus + WaveEnergyMinus)/2.
    
    ! Wave dissipation:
    ! This is the sum of a Kolmogorov dissipation of each wave polarity by itself,
    ! and dissipation due to counter propagating waves
    ! Note the use of Lratio > 1, which in effect makes the dissipation length of counter
    ! propagating waves smaller than the Kolmogorov length.
    
    WaveDissipationPlus =  LocalDissipationFactor * &
         (1. + Lratio*sqrt(WaveEnergyMinus/WavePressure)) * &
         WaveEnergyPlus**1.5
    
    WaveDissipationMinus =  LocalDissipationFactor * &
         (1. +  Lratio*sqrt(WaveEnergyPlus/WavePressure)) * &
         WaveEnergyMinus**1.5
    
  end subroutine calc_alfven_wave_dissipation
  ! ===================================================================================
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
    case('bcs')
       call write_ghost_and_boundary

    case('b0cells')
       call write_b0

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

    case('b')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          if(UseB0)then
             tmp1_BLK(:,:,:,iBlock) = &
                  sqrt(sum((B0_DGB(:,:,:,:,iBlock) + State_VGB(Bx_:Bz_,:,:,:,iBlock))**2))
          else
             tmp1_BLK(:,:,:,iBlock) = sqrt(sum(State_VGB(Bx_:Bz_,:,:,:,iBlock)**2))
          end if
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('b0')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          if(UseB0)then
             tmp1_BLK(:,:,:,iBlock) = & 
                  sqrt(sum(B0_DGB(:,:,:,:,iBlock)**2)) 
          else
             tmp1_BLK(:,:,:,iBlock) = -7777.
          end if
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('b1')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = &
               sqrt(sum(State_VGB(Bx_:Bz_,:,:,:,iBlock)**2))
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('bx0')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = B0_DGB(1,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('by0')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = B0_DGB(2,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('bz0')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = B0_DGB(3,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('bx1')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = State_VGB(Bx_,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('by1')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = State_VGB(By_,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

    case('bz1')
       do iBlock = 1, nBlock
          if(unusedBLK(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = State_VGB(Bz_,:,:,:,iBlock)
       end do
       VarValue = integrate_BLK(1,tmp1_BLK)

       
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

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, B0_DGB, StateOld_VCB
    use ModMain,       ONLY: UseB0
    use ModNumConst,   ONLY: cTolerance
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, UnitEnergyDens_, UnitX_
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Bx_, By_, Bz_, RhoUx_, RhoUy_, RhoUz_, &
         WaveFirst_, WaveLast_

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
    real    :: WaveDissipPlus,WaveDissipMinus,LengthKol,LengthCp

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------
    IsFound = .true.

    select case(NameVar)

    case('disstot','dissplus','dissminus','lkol','lcp')
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          call calc_alfven_wave_dissipation(i,j,k,iBlock,&
               WaveDissipPlus, WaveDissipMinus,LengthKol,LengthCp)

          select case(NameVar)
          case('dissplus')
             PlotVar_G(i,j,k) = max(1e-30,WaveDissipPlus* &
                  No2Si_V(UnitEnergyDens_))
             NameIdlUnit = '[J/m3]'
             NameTecUnit = '[J/m3]'

          case('dissminus')
             PlotVar_G(i,j,k) = max(1e-30,WaveDissipMinus* &
                  No2Si_V(UnitEnergyDens_))
             NameIdlUnit = '[J/m3]'
             NameTecUnit = '[J/m3]'

          case('disstot')
             PlotVar_G(i,j,k) = max(1e-30,&
                  (WaveDissipPlus + WaveDissipMinus)* &
                  No2Si_V(UnitEnergyDens_))
             NameIdlUnit = '[J/m3]'
             NameTecUnit = '[J/m3]'

          case('lkol')
             PlotVar_G(i,j,k) = LengthKol * No2Si_V(UnitX_)
             NameIdlUnit = '[m]'
             NameTecUnit = '[m]'
             NameTecVar  = 'Lkol'

          case('lcp')
             PlotVar_G(i,j,k) = LengthKol * No2Si_V(UnitX_)
             NameIdlUnit = '[m]'
             NameTecUnit = '[m]'
             NameTecVar  = 'Lcp'
          end select

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

    case('rhoerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(rho_,i,j,k,iBlock) &
               -StateOld_VCB(rho_,i,j,k,iBlock))/State_VGB(rho_,i,j,k,iBlock)
       end do; end do; end do

    case('perr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(p_,i,j,k,iBlock) &
               -StateOld_VCB(p_,i,j,k,iBlock))/State_VGB(p_,i,j,k,iBlock)
       end do; end do; end do

    case('mxerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(rhoUx_,i,j,k,iBlock) &
               -StateOld_VCB(rhoUx_,i,j,k,iBlock))/State_VGB(rhoUx_,i,j,k,iBlock)
       end do; end do; end do

    case('myerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(rhoUy_,i,j,k,iBlock) &
               -StateOld_VCB(rhoUy_,i,j,k,iBlock))/State_VGB(rhoUy_,i,j,k,iBlock)
       end do; end do; end do

    case('mzerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(rhoUz_,i,j,k,iBlock) &
               -StateOld_VCB(rhoUz_,i,j,k,iBlock))/State_VGB(rhoUz_,i,j,k,iBlock)
       end do; end do; end do

    case('bxerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(bx_,i,j,k,iBlock) &
               -StateOld_VCB(bx_,i,j,k,iBlock))/State_VGB(bx_,i,j,k,iBlock)
       end do; end do; end do

    case('byerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(by_,i,j,k,iBlock) &
               -StateOld_VCB(by_,i,j,k,iBlock))/State_VGB(by_,i,j,k,iBlock)
       end do; end do; end do

    case('bzerr')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(bz_,i,j,k,iBlock) &
               -StateOld_VCB(bz_,i,j,k,iBlock))/State_VGB(bz_,i,j,k,iBlock)
       end do; end do; end do

    case('i01err')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(WaveFirst_,i,j,k,iBlock) &
               -StateOld_VCB(WaveFirst_,i,j,k,iBlock))/State_VGB(WaveFirst_,i,j,k,iBlock)
       end do; end do; end do

    case('i02err')
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = (State_VGB(WaveLast_,i,j,k,iBlock) &
               -StateOld_VCB(WaveLast_,i,j,k,iBlock))/State_VGB(WaveLast_,i,j,k,iBlock)
       end do; end do; end do
       
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var
  !============================================================================
  subroutine user_specify_refinement(iBlock, iArea, DoRefine)

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, B0_DGB, SignB_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, far_field_BCs_BLK
    use ModNumConst, ONLY: cTiny
    use ModMain,     ONLY: x_, y_, z_, time_loop, UseB0, iNewGrid, &
                           DoThinCurrentSheet
    use ModPhysics,  ONLY: rBody

    integer, intent(in) :: iBlock, iArea
    logical,intent(out) :: DoRefine

    real :: rDotB_G(1:nI,1:nJ,0:nK+1), r_G(1:nI,1:nJ,1:nK), rMin=0.0
    integer       :: i, j, k
    integer, save :: iLevel = 0, iLastGrid = -1
    character (len=*), parameter :: NameSub = 'user_specify_refinement'
    !--------------------------------------------------------------------------

    if(.not.time_loop)then
       DoRefine = .false.

       RETURN
    end if
    
    ! Check if current sheet passess through this block
    ! (i.e. if B.r changes sign)
    if(SignB_>1 .and. DoThinCurrentSheet)then
       DoRefine = &
            maxval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))>0.0.and. &
            minval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))<0.0
    else

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

    end if
  
    ! Do gradual refinement along R if required
    if(DoRefineGradualSheet) then
       if(iNewGrid > iLastGrid ) then
          iLevel = iLevel + 1
          iLastGrid = iNewGrid
       end if

       if (iLevel > nLevelGradual) then
          DoRefine = .false.
          RETURN
       end if

       if(rRefine_I(iLevel) > rRefine_I(iLevel-1)) then
          do k=1,nK ; do j=1,nJ ; do i=1,nI
             r_G(i,j,k) = sqrt(x_BLK(i,j,k,iBlock)**2 + &
                  y_BLK(i,j,k,iBlock)**2 + &
                  z_BLK(i,j,k,iBlock)**2)
          end do; end do; end do
          rMin = minval(r_G)
          if (rMin < rRefine_I(iLevel)) DoRefine = .false.
       end if
    end if

  end subroutine user_specify_refinement
  !============================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    ! Fill ghost cells inside body for spherical grid - this subroutine only 
    ! modifies ghost cells in the r direction
    
    use ModAdvance,    ONLY: State_VGB, B0_DGB, UseElectronPressure
    use ModMain,       ONLY: East_, UseB0, UseUSerInnerBcs
    use ModGeometry,   ONLY: TypeGeometry, x_BLK, y_BLK, z_BLK, r_BLK
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, WaveFirst_, WaveLast_, &
                             Bx_, Bz_, Ux_, Uz_
    use ModMultiFluid, ONLY: MassIon_I
    use ModImplicit,   ONLY: StateSemi_VGB, iTeImpl
    use ModPhysics,    ONLY: AverageIonCharge

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i,j,k
    real    :: x, y, z, r, r_D(3), rUnit_DG(3,3), bUnit_DG(3,3), U_DG(3,3)
    real    :: FullBr, FullB, FullB_D(3), Br1_D(3), Bt1_D(3),Ewave
    real    :: NumDensIon, NumDensElectron
    real :: Runit_D(3)

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !--------------------------------------------------------------------------
    if(iSide /= East_ .or. TypeGeometry(1:9) /='spherical') &
         call CON_stop('Wrong iSide in user_set_outerBCs')

    IsFound = .true.

    if(TypeBc == 'usersemi')then
       StateSemi_VGB(iTeImpl,0,:,:,iBlock) = tChromo
       RETURN
    elseif(TypeBc == 'usersemilinear')then
       RETURN
    end if

    if (UseChromoBc) then

       ! Extrapolate density around fixed RhoChromo value
       State_VGB(Rho_,0,:,:,iBlock) = &
            2.*RhoChromo - State_VGB(rho_,1,:,:,iBlock)
       State_VGB(Rho_,-1,:,:,iBlock) = &
            2.*State_VGB(Rho_,0,:,:,iBlock)  - State_VGB(rho_,1,:,:,iBlock)

       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, 0
          NumDensIon = State_VGB(Rho_,i,j,k,iBlock)/MassIon_I(1)
          NumDensElectron = NumDensIon*AverageIonCharge
          if(UseElectronPressure)then
             State_VGB(p_,i,j,k,iBlock) = NumDensIon*tChromo
             State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*tChromo
          else
             State_VGB(p_,i,j,k,iBlock) = &
                  (NumDensIon + NumDensElectron)*tChromo
          end if
       end do; end do; end do

       do k = -1, nK+2; do j = -1, nJ+2
          Runit_D = (/ x_BLK(1,j,k,iBlock), y_BLK(1,j,k,iBlock), &
               z_BLK(1,j,k,iBlock) /) / r_BLK(1,j,k,iBlock)

          Br1_D = sum(State_VGB(Bx_:Bz_,1,j,k,iBlock)*Runit_D)*Runit_D
          Bt1_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) - Br1_D

          do i = -1, 0
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = Bt1_D
          end do
          
       end do; end do

    end if

    ! Safety
    if (UseUserInnerBcs) RETURN

    do k = -1,nK+2 ; do j = -1,nJ+2

       ! Update B1 in ghost cells (r direction only
       ! Reflect normal component, float tangential component 
       do i= 1,0,-1
          x = x_BLK(i,j,k,iBlock)
          y = y_BLK(i,j,k,iBlock)
          z = z_BLK(i,j,k,iBlock)
          r = r_BLK(i,j,k,iBlock)
          r_D = (/x,y,z/)
          rUnit_DG(:,i) = r_D/r             
          
          Br1_D = rUnit_DG(:,i) * &
               sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)*rUnit_DG(:,i))
          Bt1_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) - Br1_D
          State_VGB(Bx_:Bz_, i-1,j,k,iBlock) = Bt1_D - Br1_D
       end do

       ! Full magnetic field
       do i = -1,1
          FullB_D = B0_DGB(:,i,j,k,iBlock) + &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)
          FullB = sqrt(sum(FullB_D**2))
          FullBr = sum(rUnit_DG(:,i)*FullB_D)
          bUnit_DG(:,i) = FullB_D / FullB
       end do

       !\
       ! Velocity
       !/
       do i = 0,-1,-1
          if (UseUparBc) then
             ! Conserve momentum along field lines
             U_DG(:,i+1) = State_VGB(Ux_:Uz_,i+1,j,k,iBlock)
             U_DG(:,i) = &
                  (State_VGB(Rho_,i+1,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))* &
                  sum(U_DG(:,i+1)*bUnit_DG(:,i+1))*bUnit_DG(:,i)
             State_VGB(Ux_:Uz_,i,j,k,iBlock) = U_DG(:,i)
          else
             ! Reflect U - same as van der Holst, 2010, Cohen 2007
             State_VGB(Ux_:Uz_,i,j,k,iBlock) = -U_DG(:,i+1)
          end if
       end do

       ! set wave energy in ghost cells in the r direction
       call get_wave_energy_base(FullB, Ewave)
          
       if (UseExtrapolatedEwave) then
          if(FullBr > 0.0)then
             State_VGB(WaveFirst_,0,j,k,iBlock) = &
                  2.*Ewave - State_VGB(WaveFirst_,1,j,k,iBlock)
             State_VGB(WaveFirst_,-1,j,k,iBlock) = &
                  2.*State_VGB(WaveFirst_,0,j,k,iBlock) - &
                  State_VGB(WaveFirst_,1,j,k,iBlock)
             
             State_VGB(WaveLast_,0,j,k,iBlock) = &
                  2.*State_VGB(WaveLast_,1,j,k,iBlock) - &
                  State_VGB(WaveLast_,2,j,k,iBlock)
             State_VGB(WaveLast_,-1,j,k,iBlock) = &
                  2.*State_VGB(WaveLast_,0,j,k,iBlock) - & 
                  State_VGB(WaveLast_,1,j,k,iBlock)
          else
             State_VGB(WaveLast_,0,j,k,iBlock) = &
                  2.*Ewave - State_VGB(WaveLast_,1,j,k,iBlock)
             State_VGB(WaveLast_,-1,j,k,iBlock) = &
                  2.*State_VGB(WaveLast_,0,j,k,iBlock) - & 
                  State_VGB(WaveLast_,1,j,k,iBlock)
             
             State_VGB(WaveFirst_,0,j,k,iBlock) = &
                  2.*State_VGB(WaveFirst_,1,j,k,iBlock) - &
                  State_VGB(WaveFirst_,2,j,k,iBlock)
             State_VGB(WaveFirst_,-1,j,k,iBlock) = &
                  2.*State_VGB(WaveFirst_,0,j,k,iBlock) - & 
                  State_VGB(WaveFirst_,1,j,k,iBlock)
          end if
       else
          ! Use fixed wave energy
          if(FullBr > 0.0)then
             State_VGB(WaveFirst_,-1:0,j,k,iBlock) = Ewave          
             State_VGB(WaveLast_,-1:0,j,k,iBlock) = 0.0 
          else
             State_VGB(WaveFirst_,-1:0,j,k,iBlock) = Ewave          
             State_VGB(WaveLast_,-1:0,j,k,iBlock) = 0.0 
          end if
       end if
    end do; end do

  end subroutine user_set_outerbcs
  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance,     ONLY: State_VGB, UseElectronPressure
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame, GlobalBLK, UseB0
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: OmegaBody, AverageIonCharge, BodyRho_I, &
                              BodyTDim_I, Si2No_V, UnitTemperature_, UnitN_
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_, &
         WaveFirst_, WaveLast_, Pe_,Hyp_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: RhoBase, NumDensIon, NumDensElectron, Tbase, FullBr, EwaveBase
    real,dimension(3) :: U_D, B0_D, B1_D, B1t_D, B1r_D
    real,dimension(3) :: bUnitGhost_D, bUnitTrue_D, rUnit_D
    real,dimension(3) :: FullBGhost_D, FullBTrue_D
    real              :: RhoTrue, RhoGhost, Ur
   
    character (len=*), parameter :: NameSub = 'user_face_bcs'
    !--------------------------------------------------------------------------
    ! Check that extrapolated wave inner BC's were not chosen - can only be performed
    ! while using ghost-cells based BC's (see set_outerbc for a spherical grid)
    if (UseExtrapolatedEwave) call CON_stop('UseExtrapolatedEwave=T but you are using face_bcs')

    rUnit_D = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    !\
    ! Magnetic field
    !/
    ! Float tangential B, reflect radial B
    if(UseB0)then
       B1_D  = VarsTrueFace_V(Bx_:Bz_)
       B1r_D = sum(rUnit_D*B1_D)*rUnit_D
       B1t_D = B1_D - B1r_D
       VarsGhostFace_V(Bx_:Bz_) = B1t_D
       FullBGhost_D = B0Face_D + VarsGhostFace_V(Bx_:Bz_)
       FullBTrue_D  = B0Face_D + VarsTrueFace_V(Bx_:Bz_)
    else
       call get_coronal_b0(FaceCoords_D(x_), FaceCoords_D(y_), &
            FaceCoords_D(z_), B0_D)
       B1_D  = VarsTrueFace_V(Bx_:Bz_) - B0_D 
       B1r_D = sum(rUnit_D*B1_D)*rUnit_D
       B1t_D = B1_D - B1r_D
       VarsGhostFace_V(Bx_:Bz_) = B1t_D + B0_D
       FullBGhost_D = VarsGhostFace_V(Bx_:Bz_)
       FullBTrue_D  = VarsTrueFace_V(Bx_:Bz_)
    end if
    FullBr = sum(FullBGhost_D*rUnit_D)

    ! Density
    call get_plasma_parameters_base(FaceCoords_D, RhoBase, Tbase)

    VarsGhostFace_V(Rho_) =  2.0*RhoBase - VarsTrueFace_V(Rho_)

    RhoTrue = VarsTrueFace_V(Rho_)
    RhoGhost = VarsGhostFace_V(Rho_)

    !\
    ! Velocity
    !/
    U_D   = VarsTrueFace_V(Ux_:Uz_)
    if (UseUparBc) then
       ! Conserve momentum along field lines
       bUnitGhost_D = FullBGhost_D/sqrt(sum(FullBGhost_D**2))
       bUnitTrue_D = FullBTrue_D/sqrt(sum(FullBTrue_D**2))

       VarsGhostFace_V(Ux_:Uz_) = RhoTrue/RhoGhost* &
            sum(U_D*bUnitTrue_D)*bUnitGhost_D       

       ! Apply corotation if needed
       if(.not.UseRotatingFrame)then
          VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
               - OmegaBody*FaceCoords_D(y_)
          VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
               + OmegaBody*FaceCoords_D(x_)
       end if
    else
       ! Reflect U - same as van der Holst, 2010, Cohen 2007
       VarsGhostFace_V(Ux_:Uz_) = -U_D
       ! Apply corotation if needed
       if(.not.UseRotatingFrame)then
          VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
               - 2.0*OmegaBody*FaceCoords_D(y_)
          VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
               + 2.0*OmegaBody*FaceCoords_D(x_)
       end if
    end if

    !\
    ! Fixed wave BC's
    !/
    call get_wave_energy_base(abs(FullBr), EwaveBase)

    if (FullBr > 0. ) then
       VarsGhostFace_V(WaveFirst_) = EwaveBase*RhoGhost/RhoBase
       VarsGhostFace_V(WaveLast_) = 0.0
    else
       VarsGhostFace_V(WaveFirst_) = 0.0
       VarsGhostFace_V(WaveLast_) = EwaveBase*RhoGhost/RhoBase
    end if
  
    !\
    ! Pressure
    !/
    NumDensIon = VarsGhostFace_V(Rho_)/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    if(UseElectronPressure)then
       VarsGhostFace_V(p_) = NumDensIon*Tbase
       VarsGhostFace_V(Pe_) = NumDensElectron*Tbase
    else
       VarsGhostFace_V(p_) = (NumDensIon + NumDensElectron)*Tbase
    end if

    if(Hyp_>1) VarsGhostFace_V(Hyp_) = 0.0

  end subroutine user_face_bcs
  !============================================================================
  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModAdvance,    ONLY: State_VGB
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
