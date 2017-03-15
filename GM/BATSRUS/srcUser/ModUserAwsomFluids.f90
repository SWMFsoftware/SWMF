!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:bartvand@umich.edu  expires:12/31/2099
!==============================================================================
module ModUser

  use ModVarIndexes, ONLY: IonFirst_, nFluid
  use ModMultiFluid, ONLY: nIonFluid
  use ModMain, ONLY: nI, nJ,nK
  use ModCoronalHeating, ONLY: PoyntingFluxPerB
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_get_log_var,                &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_set_cell_boundary,          &
       IMPLEMENTED7 => user_set_resistivity,            &
       IMPLEMENTED8 => user_calc_sources,               &
       IMPLEMENTED9 => user_init_point_implicit,        &
       IMPLEMENTED10=> user_get_b0

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'AWSoM multi-fluid model'

  ! Input parameters for chromospheric inner BC's
  real    :: NchromoSi_I(IonFirst_:nFluid) = 2e17, TchromoSi = 5e4
  real    :: Nchromo_I(IonFirst_:nFluid), Tchromo

  ! variables for Parker initial condition
  real    :: nCoronaSi = 1.5e14, tCoronaSi = 1.5e6

  ! Input parameters for two-temperature effects
  real    :: TeFraction, TiFraction
  real    :: EtaPerpSi

  real    :: Mass_I(0:nIonFluid), Charge_I(0:nIonFluid)
  real    :: ReducedMass_II(0:nIonFluid,0:nIonFluid)
  real    :: CollisionCoef_II(0:nIonFluid,0:nIonFluid)

  ! variables for polar jet application              
  ! Dipole under surface                                                   
  real    :: UserDipoleStrength, UserDipoleStrengthSi
  real    :: UserDipoleDepth
  real    :: UserDipoleLatitude
  real    :: UserDipoleLongitude
  real    :: UserDipoleAxisLatitude
  real    :: UserDipoleAxisLongitude

  ! Rotating boundary condition                                             
  real    :: TbeginJet, TendJet
  logical :: IsRamping
  logical :: IsPolarDipole, IsJetBC
  real    :: DistMinJet
  real    :: DistMaxJet
  real    :: LocationMaxJet
  real    :: UmaxJet
  real    :: ProfileExponentJet
  real    :: LinearCoeffJet 
  real    :: PowerCoeffJet 

  ! Multiply the collisionfrequencies by CollisionFactor
  real :: CollisionFactor = 1.0

contains 
  !============================================================================
  subroutine user_read_inputs

    use ModMain,       ONLY: UseUserInitSession, lVerbose
    use ModProcMH,     ONLY: iProc
    use ModReadParam,  ONLY: read_line, read_command, read_var
    use ModIO,         ONLY: write_prefix, write_myname, iUnitOut

    integer :: iFluid
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
          do iFluid = IonFirst_, nFluid
             call read_var('NchromoSi', NchromoSi_I(iFluid))
          end do
          call read_var('TchromoSi', TchromoSi)

       case("#PARKERIC")
          call read_var('nCoronaSi', nCoronaSi)
          call read_var('tCoronaSi', tCoronaSi)

       case("#POLARJETDIPOLE")
          call read_var('UserDipoleStrengthSi', UserDipoleStrengthSi)
          call read_var('UserDipoleDepth', UserDipoleDepth)
          call read_var('UserDipoleLatitude',UserDipoleLatitude)
          call read_var('UserDipoleLongitude',UserDipoleLongitude)
          call read_var('UserDipoleAxisLatitude',UserDipoleAxisLatitude)
          call read_var('UserDipoleAxisLongitude',UserDipoleAxisLongitude)
          IsPolarDipole = .true.

       case("#POLARJETBC")
          call read_var('IsRamping', IsRamping)
          if(IsRamping)then
             call read_var('TbeginJet', TBeginJet)
             call read_var('TendJet',   TEndJet)
          endif
          call read_var('DistMinJet',DistMinJet)
          call read_var('DistMaxJet',DistMaxJet)
          call read_var('ProfileExponentJet',ProfileExponentJet)
          call read_var('UmaxJet',UmaxJet)
          IsJetBC = .true.

       case("#COLLISIONFACTOR")
          call read_var('CollisionFactor', CollisionFactor)

       case('#USERINPUTEND')
          if(iProc == 0 .and. lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input SOLAR CORONA ends'
          endif
          EXIT

       case default
          if(iProc == 0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs.'
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

    use ModProcMH,     ONLY: iProc
    use ModIO,         ONLY: write_prefix, iUnitOut
    use ModWaves,      ONLY: UseWavePressure, UseAlfvenWaves
    use ModAdvance,    ONLY: UseElectronPressure
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I
    use ModConst,      ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
         cElectronMass, cProtonMass
    use ModNumConst,   ONLY: cTwoPi, cDegToRad
    use ModPhysics,    ONLY: ElectronTemperatureRatio, AverageIonCharge, &
         Si2No_V, UnitTemperature_, UnitN_, UnitX_, No2Si_V, UnitT_, UnitB_, &
         UnitU_
    !UnitB_ and Unit_U are for jet only
    integer :: iIon, jIon
    real, parameter :: CoulombLog = 20.0

    character (len=*),parameter :: NameSub = 'user_init_session'

    !--------------------------------------------------------------------------

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    UseAlfvenWaves = .true.
    UseWavePressure = .true.

    ! convert to normalized units
    Nchromo_I = NchromoSi_I*Si2No_V(UnitN_)
    Tchromo = TchromoSi*Si2No_V(UnitTemperature_)

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

    Mass_I(0) = cElectronMass/cProtonMass
    Mass_I(1:) = MassIon_I
    Charge_I(0) = 1.0
    Charge_I(1:) = ChargeIon_I

    ! Coefficient for effective ion-ion collision frequencies
    do jIon = 0, nIonFluid
       do iIon = 0, nIonFluid
          ReducedMass_II(iIon,jIon) = Mass_I(iIon)*Mass_I(jIon) &
               /(Mass_I(iIon) + Mass_I(jIon))

          CollisionCoef_II(iIon,jIon) = CollisionFactor*CoulombLog &
               *sqrt(ReducedMass_II(iIon,jIon)/cProtonMass)/Mass_I(iIon) &
               *(Charge_I(iIon)*Charge_I(jIon)*cElectronCharge**2/cEps)**2 &
               /(3*(cTwoPi*cBoltzmann)**1.5)
       end do
    end do
    ! To obtain the effective ion-ion collision frequencies, the
    ! coefficients still need to be multiplied by Nion(jIon)/reducedTemp**1.5.
    ! Here, we already take care of the units.
    CollisionCoef_II = CollisionCoef_II &
         *(1/Si2No_V(UnitT_))*No2Si_V(UnitN_)/No2Si_V(UnitTemperature_)**1.5

    ! dipole (jet) parameter converted to normalized units                  
    if(IsPolarDipole)then
       UserDipoleStrength = UserDipoleStrengthSi*Si2No_V(UnitB_)
       UserDipoleLatitude = UserDipoleLatitude*cDegToRad
       UserDipoleLongitude = UserDipoleLongitude*cDegToRad
       UserDipoleAxisLatitude = UserDipoleAxisLatitude*cDegToRad
       UserDipoleAxisLongitude = UserDipoleAxisLongitude*cDegToRad
    endif

    ! the rotation's velocity profile is the following:
    ! u(z) = Ax-Bx^C
    ! LinearCoeff =     A = umax/xmax * rmax^(C-1)/(rmax^(C-1)-xmax^(C-1)) 
    ! PowerCoeff  =     B = umax/xmax             /(rmax^(C-1)-xmax^(C-1)) 
    !                   B = A * rmax^(C-1) 
    ! ProfileExponent = C = 5.14
    ! A = 3600.
    ! B = 164.^1./5.14 = 2.4226857e11
    ! xmax = rmax*C^(1/1-C) = LocationMaxJet
    ! for:
    !     DistMminJet = 2.13e-3 Rs = rmin
    !     DistMaxJet = 1.28582e-2 Rs = rmax
    !     LocationMaxJet = 8.658e-3 Rs = xmax
    !     UmaxJet = 25.1067 km/s = umax
    if(IsJetBC) then

       LocationMaxJet = DistMaxJet*ProfileExponentJet** &
            (1./(1.-ProfileExponentJet))

       UmaxJet = UmaxJet * Si2No_V(UnitU_) * 1e3 

       LinearCoeffJet = UmaxJet/LocationMaxJet &
            * DistMaxJet**(ProfileExponentJet-1.) &
            /(DistMaxJet**(ProfileExponentJet-1.) &
            - LocationMaxJet**(ProfileExponentJet-1.))

       PowerCoeffJet = LinearCoeffJet* DistMaxJet** &
            (1.-ProfileExponentJet)

    endif

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session

  !============================================================================

  subroutine user_set_ics(iBlock)

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, UseAnisoPressure
    use ModB0,         ONLY: B0_DGB
    use ModGeometry,   ONLY: Xyz_DGB, r_Blk
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, UnitN_
    use ModVarIndexes, ONLY: Rho_, Bx_, Bz_, p_, Pe_, WaveFirst_, WaveLast_, &
         Ew_
    use ModMultiFluid, ONLY: iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I, &
         MassFluid_I, iRhoIon_I, iPIon_I, MassIon_I, ChargeIon_I, IonLast_, &
         UseMultiIon, iPparIon_I, IsIon_I
    use ModWaves, ONLY: UseWavePressureLtd

    integer, intent(in) :: iBlock

    integer :: i, j, k, iFluid
    real :: x, y, z, r, Rho, NumDens_I(IonFirst_:nFluid)
    real :: RhoCorona, tCorona, uCorona
    real :: r_D(3), Br
    ! variables for iterative Parker solution
    integer :: IterCount
    real :: Ur, Ur0, Ur1, del, rTransonic, Uescape, Usound

    real, parameter :: Epsilon = 1.0e-6
    character (len=*), parameter :: NameSub = 'user_set_ics'

    !--------------------------------------------------------------------------

    ! Initially, density, electron and ion temperature are at coronal
    ! values starting from just above the boundary
    RhoCorona = nCoronaSi*Si2No_V(UnitN_)*MassIon_I(1)
    tCorona   = tCoronaSi*Si2No_V(UnitTemperature_)

    ! normalize with isothermal sound speed.
    Usound = sqrt(tCorona*(1.0+ChargeIon_I(1))/MassIon_I(1))
    Uescape = sqrt(-GBody*2.0)/Usound

    !\
    ! Initialize MHD wind with Parker's solution
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    rTransonic = 0.25*Uescape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    uCorona = rTransonic**2*exp(1.5 - 2.0*rTransonic)

    do k = MinK, MaxK ; do j = MinJ, MaxJ ; do i = MinI, MaxI
       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)
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

       Rho = rBody**2*RhoCorona*uCorona/(r**2*Ur)

       NumDens_I = Rho/MassIon_I(1)*Nchromo_I/Nchromo_I(IonFirst_)

       do iFluid = IonFirst_, nFluid
          State_VGB(iRho_I(iFluid),i,j,k,iBlock) = &
               NumDens_I(iFluid)*MassFluid_I(iFluid)
       end do

       State_VGB(iP_I(IonFirst_:),i,j,k,iBlock) = NumDens_I*Tcorona
       if(UseAnisoPressure)then
          do iFluid = IonFirst_, IonLast_
             State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  State_VGB(iP_I(iFluid),i,j,k,iBlock)
          end do
       end if
       State_VGB(Pe_,i,j,k,iBlock) = &
            sum(ChargeIon_I*State_VGB(iPIon_I,i,j,k,iBlock))

       State_VGB(iRhoUx_I,i,j,k,iBlock) = &
            State_VGB(iRho_I,i,j,k,iBlock)*Ur*x/r *Usound
       State_VGB(iRhoUy_I,i,j,k,iBlock) = &
            State_VGB(iRho_I,i,j,k,iBlock)*Ur*y/r *Usound
       State_VGB(iRhoUz_I,i,j,k,iBlock) = &
            State_VGB(iRho_I,i,j,k,iBlock)*Ur*z/r *Usound

       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
       Br = sum(B0_DGB(1:3,i,j,k,iBlock)*r_D)

       if (Br >= 0.0) then
          State_VGB(WaveFirst_,i,j,k,iBlock) =  &
               PoyntingFluxPerB*sqrt(State_VGB(iRho_I(IonFirst_),i,j,k,iBlock))
          State_VGB(WaveLast_,i,j,k,iBlock) = &
               1e-4*State_VGB(WaveFirst_,i,j,k,iBlock)
       else
          State_VGB(WaveLast_,i,j,k,iBlock) =  &
               PoyntingFluxPerB*sqrt(State_VGB(iRho_I(IonFirst_),i,j,k,iBlock))
          State_VGB(WaveFirst_,i,j,k,iBlock) = &
               1e-4*State_VGB(WaveLast_,i,j,k,iBlock)
       end if
       if(UseWavePressureLtd) State_VGB(Ew_,i,j,k,iBlock) = &
            sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModAdvance,    ONLY: State_VGB, tmp1_BLK
    use ModB0,         ONLY: B0_DGB
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: Unused_B, nBlock, x_, y_, z_, UseB0
    use ModPhysics,    ONLY: InvGammaMinus1, No2Io_V, No2Si_V, &
         UnitEnergydens_, UnitX_, UnitRho_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_, Pe_, iP_I, Ew_, &
                             rho_, rhoUx_, rhoUy_, rhoUz_
    use ModGeometry,   ONLY: R_BLK

    real, intent(out) :: VarValue
    character(len=10), intent(in) :: TypeVar 
    real, optional, intent(in) :: Radius

    integer :: i, j, k, iBlock
    real :: unit_energy, unit_mass
    real, external :: integrate_BLK
    !--------------------------------------------------------------------------
    unit_energy = No2Io_V(UnitEnergydens_)*No2Io_V(UnitX_)**3
    unit_mass   = 1.0e3*No2Si_V(UnitRho_)*No2Si_V(UnitX_)**3

    !\
    ! Define log variable to be saved::
    !/
    select case(TypeVar)

    case('eint')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             tmp1_BLK(i,j,k,iBlock) = &
                  sum(State_VGB(iP_I(IonFirst_:),i,j,k,iBlock)) &
                  + State_VGB(Pe_,i,j,k,iBlock)
          end do; end do; end do
       end do
       VarValue = unit_energy*InvGammaMinus1*integrate_BLK(1,tmp1_BLK)

    case('emag')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
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

    case('ekin')
       do iBlock=1,nBlock
          if (Unused_B(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = &
               (State_VGB(rhoUx_,:,:,:,iBlock)**2 +&
               State_VGB(rhoUy_,:,:,:,iBlock)**2 +&
               State_VGB(rhoUz_,:,:,:,iBlock)**2)/&
               State_VGB(rho_  ,:,:,:,iBlock)
       end do
       VarValue = unit_energy*0.5*integrate_BLK(1,tmp1_BLK)

    case('ew')
       do iBlock=1,nBlock
          if (Unused_B(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = State_VGB(Ew_,:,:,:,iBlock)
       end do
       VarValue = unit_energy*integrate_BLK(1,tmp1_BLK)

    case('mass')
       do iBlock=1,nBlock
          if (Unused_B(iBlock)) CYCLE
          tmp1_BLK(:,:,:,iBlock) = &
               State_VGB(rho_,:,:,:,iBlock)/R_BLK(:,:,:,iBlock)
       end do
       VarValue = unit_mass*integrate_BLK(1,tmp1_BLK)

    case('vol')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          tmp1_BLK(:,:,:,iBlock) = 1.0
       end do
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

    use ModAdvance,        ONLY: UseElectronPressure, UseAnisoPressure, &
         State_VGB, Source_VC
    use ModB0,             ONLY: set_b0_face
    use ModChromosphere,   ONLY: DoExtendTransitionRegion, extension_factor, &
         TeSi_C, get_tesi_c
    use ModCoronalHeating, ONLY: IsNewBlockAlfven, CoronalHeating_C, &
         apportion_coronal_heating, get_block_heating, get_wave_reflection
    use ModFaceValue,      ONLY: calc_face_value
    use ModMultiFluid,     ONLY: IonLast_
    use ModPhysics,        ONLY: No2Si_V, UnitT_, UnitEnergyDens_
    use ModVarIndexes,     ONLY: WaveFirst_, WaveLast_

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k
    real :: QPerQtotal_I(IonFirst_:IonLast_)
    real :: QparPerQtotal_I(IonFirst_:IonLast_)
    real :: QePerQtotal

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------

    IsFound = .true.

    select case(NameVar)
    case('refl')
       Source_VC(WaveFirst_:WaveLast_,:,:,:) = 0.0
       call set_b0_face(iBlock)
       call calc_face_value(.false., iBlock)
       IsNewBlockAlfven = .true.
       call get_wave_reflection(iBlock)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = Source_VC(WaveLast_,i,j,k) &
               /sqrt(State_VGB(WaveFirst_,i,j,k,iBlock) &
               *     State_VGB(WaveLast_,i,j,k,iBlock))/No2Si_V(UnitT_)
          Source_VC(WaveFirst_:WaveLast_,i,j,k) = 0.0
       end do; end do; end do
       NameIdlUnit = '1/s'
       NameTecUnit = '1/s'

    case('qheat')
       call get_block_heating(iBlock)
       if(DoExtendTransitionRegion) call get_tesi_c(iBlock, TeSi_C)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(DoExtendTransitionRegion) CoronalHeating_C(i,j,k) = &
               CoronalHeating_C(i,j,k)/extension_factor(TeSi_C(i,j,k))
          PlotVar_G(i,j,k) = CoronalHeating_C(i,j,k) &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitT_)
       end do; end do; end do
       NameIdlUnit = 'J/m^3/s'
       NameTecUnit = 'J/m^3/s'

    case('qebyq', 'qparbyq')
       ! Not yet generalized to multi-fluid
       if(UseElectronPressure)then
          call set_b0_face(iBlock)
          call calc_face_value(.false., iBlock)
          call get_block_heating(iBlock)
          if(DoExtendTransitionRegion) call get_tesi_c(iBlock, TeSi_C)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(DoExtendTransitionRegion) CoronalHeating_C(i,j,k) = &
                  CoronalHeating_C(i,j,k)/extension_factor(TeSi_C(i,j,k))
             call apportion_coronal_heating(i, j, k, iBlock, &
                  CoronalHeating_C(i,j,k), QPerQtotal_I, QparPerQtotal_I, &
                  QePerQtotal)
             select case(NameVar)
             case('qebyq')
                PlotVar_G(i,j,k) = QePerQtotal
             case('qparbyq')
                if(UseAnisoPressure) &
                     PlotVar_G(i,j,k) = QparPerQtotal_I(IonFirst_)
             end select
          end do; end do; end do
       else
          PlotVar_G(i,j,k) = 0.0
       end if
       NameIdlUnit = 'J/m^3/s'
       NameTecUnit = 'J/m^3/s'

    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    ! Fill ghost cells inside body for spherical grid - this subroutine only 
    ! modifies ghost cells in the r direction

    use BATL_lib,      ONLY: CellSize_DB, x_, y_, z_
    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, UseAnisoPressure
    use ModB0,         ONLY: B0_DGB
    use ModGeometry,   ONLY: TypeGeometry, Xyz_DGB, r_BLK
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless
    use ModPhysics,    ONLY: InvGammaMinus1, GBody, rBody, &
         UnitU_ ! for jet only
    use ModVarIndexes, ONLY: Pe_, Bx_, Bz_, WaveFirst_, WaveLast_, Ew_, &
         Ehot_, p_
    use ModMultiFluid, ONLY: MassIon_I, iRhoIon_I, ChargeIon_I, IonLast_, &
         iRho_I, MassFluid_I, iUx_I, iUz_I, iRhoUx_I, iRhoUz_I, iPIon_I, &
         iP_I, iPparIon_I, IsIon_I
    use ModImplicit,   ONLY: StateSemi_VGB, iTeImpl
    use ModWaves,      ONLY: UseWavePressureLtd
    ! Below is for jet only
    use ModMain,       ONLY: time_simulation, time_accurate
    use ModConst,      ONLY: cPi
    use ModCoordTransform, ONLY: rlonlat_to_xyz, cross_product
    ! Above is for jet only                    
    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: Minor_, Major_
    integer :: i, j, k, iFluid, iRho, iRhoUx, iRhoUz, iP
    real    :: Br1_D(3), Bt1_D(3), Runit_D(3)
    real    :: FullB_D(3), SignBr
    real    :: Gamma
    real    :: U, U_D(3), Bdir_D(3)

    ! Below is for jet only  
    real    :: Framp
    real    :: Xyz_D(3), JetCenter_D(3)
    real    :: CrossProduct_D(3), Urot_D(3)
    real    :: DistanceJet
    real    :: Ucoeff = 0.0
    real    :: r
    ! Above is for jet only                    

    character (len=*), parameter :: NameSub = 'user_set_cell_boundary'

    !--------------------------------------------------------------------------

    if(iSide /= 1 .or. TypeGeometry(1:9) /='spherical') &
         call CON_stop('Wrong iSide in user_set_cell_boundary')

    IsFound = .true.

    select case(TypeBc)
    case('usersemi','user_semi')
       IsFound = .true.
       StateSemi_VGB(iTeImpl,0,:,:,iBlock) = Tchromo
       RETURN
    case('usersemilinear','user_semilinear')
       IsFound = .true.
       ! Value was already set to zero in ModCellBoundary     
       RETURN
       ! jet BC                          
    case('user')
       IsFound = .true.
    case default
       IsFound = .false.
       RETURN
    end select

    ! This part is for jet only: rotation starts immediately or gradually
    if(IsJetBC)then
       ! Check if time accurate is set.              
       if(time_accurate .and. IsRamping)then
          if(time_simulation<TbeginJet)then
             Framp = 0.0
          elseif(time_simulation>TendJet)then
             Framp = 1.0
          else
             Framp = 0.5 * (1 - cos(2*cPi*(time_simulation-TbeginJet) / &
                  (TendJet-TbeginJet)))
          endif
       elseif(.not. IsRamping)then
          Framp = 1.0
       else
          Framp = 0.0
       endif
    endif
    ! This part for jet only ends here.

    do k = MinK, MaxK; do j = MinJ, MaxJ

       Runit_D = Xyz_DGB(:,1,j,k,iBlock) / r_BLK(1,j,k,iBlock)

       Br1_D = sum(State_VGB(Bx_:Bz_,1,j,k,iBlock)*Runit_D)*Runit_D
       Bt1_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) - Br1_D

       ! Set B1r=0, and B1theta = B1theta(1) and B1phi = B1phi(1)
       do i = MinI, 0
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = Bt1_D
       end do

       do iFluid = IonFirst_, nFluid
          iRho = iRho_I(iFluid)

          do i = MinI, 0
             ! exponential scaleheight
             State_VGB(iRho,i,j,k,iBlock) = &
                  Nchromo_I(iFluid)*MassFluid_I(iFluid)*exp(-GBody/rBody &
                  *MassFluid_I(iFluid)/Tchromo &
                  *(rBody/r_BLK(i,j,k,iBlock) - 1.0))

             ! Fix ion temperature T_s
             State_VGB(iP_I(iFluid),i,j,k,iBlock) = Tchromo &
                  *State_VGB(iRho,i,j,k,iBlock)/MassFluid_I(iFluid)

             if(UseAnisoPressure .and. IsIon_I(iFluid)) &
                  State_VGB(iPparIon_I(iFluid),i,j,k,iBlock) = &
                  State_VGB(iP_I(iFluid),i,j,k,iBlock)
          end do
       end do

       FullB_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) + B0_DGB(:,1,j,k,iBlock)
       SignBr = sign(1.0, sum(Xyz_DGB(:,1,j,k,iBlock)*FullB_D))
       if(SignBr < 0.0)then
          Major_ = WaveLast_
          Minor_ = WaveFirst_
       else
          Major_ = WaveFirst_
          Minor_ = WaveLast_
       end if

       do i = MinI, 0
          ! Te = T_s and ne = sum(q_s n_s)
          if(UseElectronPressure) State_VGB(Pe_,i,j,k,iBlock) = &
               sum(ChargeIon_I*State_VGB(iPIon_I,i,j,k,iBlock))

          ! Outgoing wave energy
          State_VGB(Major_,i,j,k,iBlock) = PoyntingFluxPerB &
               *sqrt(State_VGB(iRho,i,j,k,iBlock))

          ! Ingoing wave energy
          State_VGB(Minor_,i,j,k,iBlock) = 0.0

          if(UseWavePressureLtd) State_VGB(Ew_,i,j,k,iBlock) &
               = sum(State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock))
       end do

       ! At the inner boundary this seems to be unnecessary... 
       ! Ehot=0 may be sufficient?!
       if(Ehot_ > 1)then
          if(UseHeatFluxCollisionless)then
             iP = p_; if(UseElectronPressure) iP = Pe_
             do i = MinI, 0
                call get_gamma_collisionless(Xyz_DGB(:,i,j,k,iBlock), Gamma)
                State_VGB(Ehot_,i,j,k,iBlock) = &
                     State_VGB(iP,i,j,k,iBlock)*(1.0/(Gamma-1) -InvGammaMinus1)
             end do
          else
             State_VGB(Ehot_,MinI:0,j,k,iBlock) = 0.0
          end if
       end if

    end do; end do

    ! BEGINNING of rotating BC (jet) part
    if(IsJetBC)then
       do k = MinK, MaxK; do j = MinJ, MaxJ

          ! Get total B direction at surface
          FullB_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) &
               + 0.5*(B0_DGB(:,0,j,k,iBlock) + B0_DGB(:,1,j,k,iBlock))
          Bdir_D = FullB_D/sqrt(max(sum(FullB_D**2), 1e-30))

          Xyz_D = Xyz_DGB(:,1,j,k,iBlock)
          r  = sqrt(sum(Xyz_D**2))

          ! Center is in the direction of the jet's rotation vector 
          call rlonlat_to_xyz( &
               (/r, UserDipoleLongitude, UserDipoleLatitude/), &
               JetCenter_D &
               )
          DistanceJet = sqrt(sum((Xyz_D-JetCenter_D)**2))

          ! Calculate the rotation speed profile and velocity vector
          if(DistanceJet < DistMinJet .or. DistanceJet > DistMaxJet)then
             Urot_D = 0.0
          else
             Ucoeff = LinearCoeffJet*DistanceJet &
                  - PowerCoeffJet*DistanceJet**ProfileExponentJet

             CrossProduct_D = cross_product(JetCenter_D,Xyz_D)
             Urot_D = Ucoeff * CrossProduct_D / sqrt(sum(CrossProduct_D**2))
          endif

          do iFluid = IonFirst_, IonLast_
             iRho = iRho_I(iFluid)
             iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)

             ! BC with reflected velocity for each cell 
             ! (relative to the surface)
             do i = MinI, 0

                ! Calculate velocity in cell 1-i                           
                u_D = State_VGB(iRhoUx:iRhoUz,1-i,j,k,iBlock) &
                     /State_VGB(iRho,1-i,j,k,iBlock)
                u_D = u_D - Urot_D

                ! u_ghost = 2*u_par - u_cell                  
                u_D = 2*(sum(u_D*Bdir_D))*Bdir_D - u_D
                u_D = u_D + Urot_D

                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     u_D*State_VGB(iRho,i,j,k,iBlock)
             end do
          end do
       end do;end do
       ! END of rotating BC (jet) part

    else
       do iFluid = IonFirst_, IonLast_
          iRho = iRho_I(iFluid)
          iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)

          do k = MinK, MaxK; do j = MinJ, MaxJ
             FullB_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) &
                  + 0.5*(B0_DGB(:,0,j,k,iBlock) + B0_DGB(:,1,j,k,iBlock))
             Bdir_D = FullB_D/sqrt(max(sum(FullB_D**2), 1e-30))

             ! Copy field-aligned velocity component. 
             ! Reflect the other components
             do i = MinI, 0
                U_D = State_VGB(iRhoUx:iRhoUz,1-i,j,k,iBlock) &
                     /State_VGB(iRho,1-i,j,k,iBlock)
                U   = sum(U_D*Bdir_D); U_D = U_D - U*Bdir_D
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     (U*Bdir_D - U_D)*State_VGB(iRho,i,j,k,iBlock)
             end do
          end do; end do
       end do
    endif

  end subroutine user_set_cell_boundary

  !============================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitTemperature_, UnitX_, UnitT_
    use ModVarIndexes, ONLY: Rho_, Pe_

    integer, intent(in) :: iBlock
    real,    intent(out):: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    integer :: i, j, k
    real :: Te, TeSi

    character (len=*), parameter :: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------

    do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
       Te = TeFraction*State_VGB(Pe_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       TeSi = Te*No2Si_V(UnitTemperature_)

       Eta_G(i,j,k) = EtaPerpSi/TeSi**1.5 *Si2No_V(UnitX_)**2/Si2No_V(UnitT_)
    end do; end do; end do

  end subroutine user_set_resistivity

  !============================================================================

  subroutine user_calc_sources(iBlock)

    use ModAdvance,    ONLY: State_VGB, Source_VC, UseElectronPressure, &
         UseAnisoPressure
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I, iRhoIon_I, iRhoUxIon_I, &
         iRhoUyIon_I, iRhoUzIon_I, iPIon_I, iPparIon_I
    use ModPhysics,    ONLY: GammaMinus1, InvGammaMinus1
    use ModPointImplicit, ONLY: UsePointImplicit, IsPointImplSource
    use ModVarIndexes, ONLY: nVar, Energy_, Pe_, Bx_, Bz_
    use ModConst,      ONLY: cElectronMass, cProtonMass
    use ModB0,         ONLY: UseB0, B0_DGB

    integer, intent(in) :: iBlock

    integer :: i, j, k
    integer :: iRhoUx, iRhoUz, iP, iEnergy, iIon, jIon, iIonFirst, iPpar
    real :: ReducedTemp, CollisionRate, Coef
    real :: Du2, Phi, Psi, RelativeDrift
    real, dimension(nIonFluid) :: RhoIon_I, ChargeDensIon_I, PparIon_I
    real, dimension(0:nIonFluid) :: Ux_I, Uy_I, Uz_I, P_I, N_I, T_I
    real :: U_D(3), Du_D(3), Me_D(3), B_D(3)
    real :: State_V(nVar), Source_V(nVar+nFluid)

    character(len=*), parameter :: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------
    ! Do not provide explicit source term when point-implicit scheme is used
    ! IsPointImplSource is true only when called from ModPointImplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iIonFirst = 1
    if(UseElectronPressure) iIonFirst = 0

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       State_V = State_VGB(:,i,j,k,iBlock)

       RhoIon_I = State_V(iRhoIon_I)
       Ux_I(1:) = State_V(iRhoUxIon_I)/RhoIon_I
       Uy_I(1:) = State_V(iRhoUyIon_I)/RhoIon_I
       Uz_I(1:) = State_V(iRhoUzIon_I)/RhoIon_I
       P_I(1:) = State_V(iPIon_I)
       N_I(1:) = RhoIon_I/MassIon_I
       T_I(1:) = P_I(1:)/N_I(1:)

       if(UseElectronPressure)then
          ChargeDensIon_I = ChargeIon_I*N_I(1:)
          P_I(0) = State_V(Pe_)
          N_I(0) = sum(ChargeDensIon_I)
          T_I(0) = P_I(0)/N_I(0)
          Ux_I(0) = sum(ChargeDensIon_I*Ux_I(1:))/N_I(0)
          Uy_I(0) = sum(ChargeDensIon_I*Uy_I(1:))/N_I(0)
          Uz_I(0) = sum(ChargeDensIon_I*Uz_I(1:))/N_I(0)
       end if

       if(UseAnisoPressure) PparIon_I = State_V(iPparIon_I)

       Source_V = 0.0

       do iIon = iIonFirst, nIonFluid

          if(iIon == 0)then
             Me_D = 0.0
             iP = Pe_
          else
             iRhoUx = iRhoUxIon_I(iIon); iRhoUz = iRhoUzIon_I(iIon)
             iP = iPIon_I(iIon)
             if(UseAnisoPressure) iPpar = iPparIon_I(iIon)
          end if

          do jIon = iIonFirst, nIonFluid
             if(iIon == jIon .and. .not.(UseAnisoPressure.and.iIon>0)) CYCLE

             ReducedTemp = (Mass_I(jIon)*T_I(iIon)+Mass_I(iIon)*T_I(jIon)) &
                  /(Mass_I(iIon) + Mass_I(jIon))

             ! Turbulence modifies the collision rate, but we do not
             ! incorporate that here
             CollisionRate = CollisionCoef_II(iIon,jIon) &
                  *N_I(jIon)/(ReducedTemp*sqrt(ReducedTemp))

             if(UseAnisoPressure .and. iIon>0)then
                Source_V(iPpar) = Source_V(iPpar) &
                     + (P_I(iIon) - PparIon_I(iIon))*CollisionRate &
                     *Mass_I(iIon)/ReducedMass_II(iIon,jIon)

                if(iIon == jIon) CYCLE
             end if

             Du_D = (/ Ux_I(jIon) - Ux_I(iIon), Uy_I(jIon) - Uy_I(iIon), &
                  Uz_I(jIon) - Uz_I(iIon) /)

             Du2 = sum(Du_D**2)

             if(UseAnisoPressure)then
                Psi = 1.0
                Phi = 1.0
             else
                RelativeDrift = &
                     sqrt(Du2/(2.0*ReducedTemp/ReducedMass_II(iIon,jIon)))

                ! Velocity dependent correction factors
                ! The Phi correction factor is simplified (Nakada, 1970)
                Psi = exp(-RelativeDrift**2)
                Phi = 1.0/(1.0 + 0.74*RelativeDrift**3)
             end if

             Coef = N_I(iIon)*ReducedMass_II(iIon,jIon)*CollisionRate

             if(iIon == 0)then
                Me_D = Me_D + Mass_I(0)*N_I(0)*CollisionRate*Du_D*Phi
             else
                Source_V(iRhoUx:iRhoUz) = Source_V(iRhoUx:iRhoUz) &
                     + RhoIon_I(iIon)*CollisionRate*Du_D*Phi

                if(UseAnisoPressure)then
                   if(UseB0) then
                      B_D = B0_DGB(:,i,j,k,iBlock) &
                           + State_VGB(Bx_:Bz_,i,j,k,iBlock)
                   else
                      B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
                   end if
                   B_D = B_D/max(sqrt(sum(B_D**2)), 1e-15)

                   Source_V(iPpar) = Source_V(iPpar) &
                        + Coef*(2.0*(T_I(jIon) - T_I(iIon))/Mass_I(jIon) &
                        + (1.0/3.0)*Du2 + (sum(Du_D*B_D))**2)
                end if
             end if

             Source_V(iP) = Source_V(iP) &
                  + Coef*(2.0*(T_I(jIon) - T_I(iIon))/Mass_I(jIon)*Psi &
                  + (2.0/3.0)*Du2*Phi)
          end do
       end do

       do iIon = 1, nIonFluid
          iRhoUx = iRhoUxIon_I(iIon); iRhoUz = iRhoUzIon_I(iIon)
          iP = iPIon_I(iIon)
          iEnergy = Energy_ + IonFirst_ - 2 + iIon

          if(UseElectronPressure) Source_V(iRhoUx:iRhoUz) &
               = Source_V(iRhoUx:iRhoUz) + ChargeDensIon_I(iIon)/N_I(0)*Me_D

          U_D = (/ Ux_I(iIon), Uy_I(iIon), Uz_I(iIon) /)
          Source_V(iEnergy) = Source_V(iEnergy) + InvGammaMinus1*Source_V(iP) &
               + sum(U_D*Source_V(iRhoUx:iRhoUz))
       end do

       Source_VC(:,i,j,k) = Source_VC(:,i,j,k) + Source_V

    end do; end do; end do

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModAdvance, ONLY: UseElectronPressure, UseAnisoPressure
    use ModMultiFluid, ONLY: iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I, iPIon_I, &
         iPparIon_I, IonFirst_, IonLast_
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    use ModVarIndexes, ONLY: nVar, Pe_

    logical :: IsPointImpl_V(nVar)
    integer :: iVar, iPointImplVar, nPointImplVar, iFluid

    character(len=*), parameter :: NameSub = 'user_init_point_implicit'
    !--------------------------------------------------------------------------

    IsPointImpl_V = .false.

    ! All ion momenta and pressures are implicit
    IsPointImpl_V(iRhoUxIon_I) = .true.
    IsPointImpl_V(iRhoUyIon_I) = .true.
    IsPointImpl_V(iRhoUzIon_I) = .true.
    IsPointImpl_V(iPIon_I)     = .true.
    if(UseElectronPressure) IsPointImpl_V(Pe_) = .true.
    if(UseAnisoPressure)then
       do iFluid = IonFirst_, IonLast_
          IsPointImpl_V(iPparIon_I(iFluid)) = .true.
       end do
    end if
    nPointImplVar = count(IsPointImpl_V)

    allocate(iVarPointImpl_I(nPointImplVar))

    iPointImplVar = 0
    do iVar = 1, nVar
       if(.not. IsPointImpl_V(iVar)) CYCLE
       iPointImplVar = iPointImplVar + 1
       iVarPointImpl_I(iPointImplVar) = iVar
    end do

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.
    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

  !===================================================================== 

  subroutine user_get_b0(x, y, z, B0_D)

    use ModCoordTransform, ONLY: rlonlat_to_xyz
    use ModGeometry, ONLY: RadiusMin

    real, intent(in)   :: x, y, z
    real, intent(inout):: B0_D(3)

    real :: Xyz_D(3), Dp, rInv, r2Inv, r3Inv, Dipole_D(3)
    real :: UserDipoleAxis_D(3)

    character(len=*), parameter :: NameSub = 'user_get_b0'

    !-------------------------------------------------------------------    

    ! Center of dipole shifted by UserDipoleDepth below RadiusMin
    call rlonlat_to_xyz( &
         (/RadiusMin-UserDipoleDepth, &
         UserDipoleLongitude, UserDipoleLatitude/), &
         Xyz_D &
         )
    ! Normalize with depth
    Xyz_D = ((/x, y, z/) - Xyz_D)/UserDipoleDepth

    call rlonlat_to_xyz(&
         (/1.,UserDipoleAxisLongitude, UserDipoleAxisLatitude/), &
         UserDipoleAxis_D)

    ! Determine radial distance and powers of it                            
    rInv  = 1.0/sqrt(sum(Xyz_D**2))
    r2Inv = rInv**2
    r3Inv = rInv*r2Inv

    ! Compute dipole moment of the intrinsic magnetic field B0.             
    Dipole_D = UserDipoleStrength * UserDipoleAxis_D 

    Dp = 3*sum(Dipole_D*Xyz_D)*r2Inv 

    B0_D = B0_D + (Dp*Xyz_D - Dipole_D) * r3Inv
  end subroutine user_get_b0

  !=====================================================================    

end module ModUser
