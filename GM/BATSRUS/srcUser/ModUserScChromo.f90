!This code is a copyright protected software (c) 2002- University of Michigan
!==============================================================================
module ModUser

  use ModMain, ONLY: nI, nJ,nK
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_get_log_var,                &
       IMPLEMENTED5 => user_set_plot_var,               &
       IMPLEMENTED6 => user_set_cell_boundary,          &
       IMPLEMENTED7 => user_set_face_boundary,          &
       IMPLEMENTED8 => user_set_resistivity,            &
       IMPLEMENTED9 => user_initial_perturbation

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Chromosphere to solar wind model with Alfven waves - Oran, van der Holst'

  ! Input parameters for chromospheric inner BC's
  real    :: nChromoSi = 2e17, tChromoSi = 5e4
  real    :: nChromo, RhoChromo, tChromo
  real    :: PoyntingFluxPerBSi = 1e6, PoyntingFluxPerB
  logical :: UseUparBc = .false.

  ! variables for Parker initial condition
  real    :: nCoronaSi = 1.5e14, tCoronaSi = 1.5e6

  ! Dipole test
  real    :: DipoleTiltDeg = 0.0

  ! Input parameters for two-temperature effects
  real    :: TeFraction, TiFraction
  real    :: EtaPerpSi

contains 
  !============================================================================
  subroutine user_read_inputs

    use ModMain,      ONLY: UseUserInitSession, lVerbose
    use ModProcMH,    ONLY: iProc
    use ModReadParam, ONLY: read_line, read_command, read_var
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut
    use ModPhysics,   ONLY: DipoleStrengthSi

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
          call read_var('nChromoSi', nChromoSi)
          call read_var('tChromoSi', tChromoSi)
          call read_var('PoyntingFluxPerBSi', PoyntingFluxPerBSi)

       case("#LINETIEDBC")
          call read_var('UseUparBc', UseUparBc)

       case('#SOLARDIPOLE')
          call read_var('DipoleStrengthSi',DipoleStrengthSi)
          call read_var('DipoleTiltDeg',DipoleTiltDeg)

       case("#PARKERIC")
          call read_var('nCoronaSi', nCoronaSi)
          call read_var('tCoronaSi', tCoronaSi)

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

    use EEE_ModCommonVariables, ONLY: UseCme
    use EEE_ModMain,   ONLY: EEE_initialize
    use ModMain,       ONLY: UseMagnetogram
    use ModProcMH,     ONLY: iProc
    use ModIO,         ONLY: write_prefix, iUnitOut
    use ModWaves,      ONLY: UseWavePressure, UseAlfvenWaves
    use ModAdvance,    ONLY: UseElectronPressure
    use ModMultiFluid, ONLY: MassIon_I
    use ModConst,      ONLY: cElectronCharge, cLightSpeed, cBoltzmann, cEps, &
         cElectronMass
    use ModNumConst,   ONLY: cTwoPi, cDegToRad
    use ModPhysics,    ONLY: ElectronTemperatureRatio, AverageIonCharge, &
         Si2No_V, UnitTemperature_, UnitN_, UnitX_, UnitB_, UnitU_, &
         UnitEnergyDens_, SinThetaTilt, CosThetaTilt, BodyNDim_I, BodyTDim_I, g

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
    nChromo = nChromoSi*Si2No_V(UnitN_)
    RhoChromo = nChromo*MassIon_I(1)
    tChromo = tChromoSi*Si2No_V(UnitTemperature_)

    PoyntingFluxPerB = PoyntingFluxPerBSi &
         *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitU_)/Si2No_V(UnitB_)

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

    if(UseCme) call EEE_initialize(BodyNDim_I(1), BodyTDim_I(1), g)

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session
  !============================================================================
  subroutine user_set_ics(iBlock)

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, B0_DGB, UseElectronPressure, &
         UseAnisoPressure
    use ModCoronalHeating, ONLY: UseTurbulentCascade
    use ModGeometry,   ONLY: Xyz_DGB, r_Blk
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         UnitN_, AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_, Pe_, &
         Ppar_, WaveFirst_, WaveLast_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: x, y, z, r, Rho, NumDensIon, NumDensElectron
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

    do k = MinK,MaxK ; do j = MinJ,MaxJ ; do i = MinI,MaxI
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

       NumDensIon = Rho/MassIon_I(1)
       NumDensElectron = NumDensIon*AverageIonCharge

       if(UseElectronPressure)then
          State_VGB(p_,i,j,k,iBlock) = NumDensIon*tCorona
          State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*tCorona
          if(UseAnisoPressure) &
               State_VGB(Ppar_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock)
       else
          State_VGB(p_,i,j,k,iBlock) = (NumDensIon + NumDensElectron)*tCorona
       end if
       State_VGB(Rho_,i,j,k,iBlock) = Rho

       State_VGB(RhoUx_,i,j,k,iBlock) = Rho*Ur*x/r *Usound
       State_VGB(RhoUy_,i,j,k,iBlock) = Rho*Ur*y/r *Usound
       State_VGB(RhoUz_,i,j,k,iBlock) = Rho*Ur*z/r *Usound

       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
       Br = sum(B0_DGB(1:3,i,j,k,iBlock)*r_D)

       if (Br >= 0.0) then
          State_VGB(WaveFirst_,i,j,k,iBlock) =  &
               PoyntingFluxPerB*sqrt(State_VGB(rho_,i,j,k,iBlock))
          if(UseTurbulentCascade)then
             State_VGB(WaveLast_,i,j,k,iBlock) = &
                  1e-4*State_VGB(WaveFirst_,i,j,k,iBlock)
          else
             State_VGB(WaveLast_,i,j,k,iBlock) = 1e-30
          end if
       else
          State_VGB(WaveLast_,i,j,k,iBlock) =  &
               PoyntingFluxPerB*sqrt(State_VGB(rho_,i,j,k,iBlock))
          if(UseTurbulentCascade)then
             State_VGB(WaveFirst_,i,j,k,iBlock) = &
                  1e-4*State_VGB(WaveLast_,i,j,k,iBlock)
          else
             State_VGB(WaveFirst_,i,j,k,iBlock) = 1e-30
          end if
       end if

    end do; end do; end do

  end subroutine user_set_ics
  !============================================================================
  subroutine write_ghost_and_boundary

    ! output ghost cells and true cells values for blocks that cross the
    ! inner boundary, in a meridional plane close to the x=0 plane

    ! ------ SPHERICAL GRID ONLY ------------

    use ModProcMH,     ONLY: iProc
    use ModIoUnit,     ONLY: UNITTMP_
    use ModVarIndexes, ONLY: WaveFirst_, WaveLast_, rho_, p_, Bx_, Bz_
    use ModMain,       ONLY: iteration_number, nBlock
    use ModGeometry,   ONLY: Xyz_DGB, r_BLK
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
          do i=MinI,MaxI

             x=Xyz_DGB(x_,i,j,k,iBLK)
             ! Extract data in the x=0 plane only
             if (abs(X) > 0.1) CYCLE

             y = Xyz_DGB(y_,i,j,k,iBLK)
             z = Xyz_DGB(z_,i,j,k,iBLK)
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
    use ModMain,       ONLY: iteration_number, nBlock
    use ModGeometry,   ONLY: Xyz_DGB, r_BLK
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
       do k=MinK,MaxK ; do i=MinI,MaxI

          x=Xyz_DGB(x_,i,j,k,iBlock)
          ! Extract data in the x=0 plane only
          if (abs(X) > 0.1) CYCLE

          y = Xyz_DGB(y_,i,j,k,iBlock)
          z = Xyz_DGB(z_,i,j,k,iBlock)
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
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, B0_DGB, UseElectronPressure
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: Unused_B, nBlock, x_, y_, z_, UseB0
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
          if(Unused_B(iBlock)) CYCLE
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

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure, &
         UseAnisoPressure, B0_DGB, StateOld_VCB, Source_VC
    use ModChromosphere, ONLY: DoExtendTransitionRegion, extension_factor, &
         get_tesi_c, TeSi_C
    use ModConst,      ONLY: rSun
    use ModCoronalHeating, ONLY: get_block_heating, CoronalHeating_C, &
         apportion_coronal_heating, IsNewBlockAlfven, get_wave_reflection
    use ModMain,       ONLY: UseB0, UseRotatingFrame
    use ModPhysics,    ONLY: No2Si_V, UnitTemperature_, UnitEnergyDens_, &
         UnitX_, UnitU_, UnitB_, UnitT_, OmegaBody
    use ModRadiativeCooling, ONLY: RadCooling_C, get_radiative_cooling
    use ModGeometry,   ONLY: Xyz_DGB
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Bx_, By_, Bz_, RhoUx_, RhoUy_, &
         RhoUz_, WaveFirst_, WaveLast_
    use ModCoordTransform,  ONLY: xyz_to_sph
    use ModFaceValue, ONLY: calc_face_value
    use ModB0, ONLY: set_b0_face

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k
    real    :: U_D(3), B_D(3), r, phi, theta
    real    :: sintheta, sinphi, costheta, cosphi
    real    :: QeFraction, QparFraction

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !--------------------------------------------------------------------------
    IsFound = .true.

    select case(NameVar)
    case('te')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
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
       NameTecUnit = '[K]'
       do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
          PlotVar_G(i,j,k) = TiFraction*State_VGB(p_,i,j,k,iBlock) &
               /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
       end do; end do; end do

    case('qrad')
       call get_tesi_c(iBlock, TeSi_C)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          call get_radiative_cooling(i, j, k, iBlock, TeSi_C(i,j,k), &
               RadCooling_C(i,j,k))
          PlotVar_G(i,j,k) = RadCooling_C(i,j,k) &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitT_)
       end do; end do; end do
       NameIdlUnit = 'J/m^3/s'
       NameTecUnit = 'J/m^3/s'

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
       ! some of the heating terms need face values
       call set_b0_face(iBlock)
       call calc_face_value(.false., iBlock)
       call get_block_heating(iBlock)
       call get_tesi_c(iBlock, TeSi_C)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          PlotVar_G(i,j,k) = CoronalHeating_C(i,j,k) &
               /extension_factor(TeSi_C(i,j,k)) &
               *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitT_)
       end do; end do; end do
       NameIdlUnit = 'J/m^3/s'
       NameTecUnit = 'J/m^3/s'

    case('qebyq', 'qparbyq')
       if(UseElectronPressure)then
          call set_b0_face(iBlock)
          call calc_face_value(.false., iBlock)
          call get_block_heating(iBlock)
          if(DoExtendTransitionRegion) call get_tesi_c(iBlock, TeSi_C)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(DoExtendTransitionRegion) CoronalHeating_C(i,j,k) = &
                  CoronalHeating_C(i,j,k)/extension_factor(TeSi_C(i,j,k))
             call apportion_coronal_heating(i, j, k, iBlock, &
                  CoronalHeating_C(i,j,k), QeFraction, QparFraction)
             select case(NameVar)
             case('qebyq')
                PlotVar_G(i,j,k) = QeFraction
             case('qparbyq')
                if(UseAnisoPressure) PlotVar_G(i,j,k) = QparFraction
             end select
          end do; end do; end do
       end if
       NameIdlUnit = 'J/m^3/s'
       NameTecUnit = 'J/m^3/s'

       ! Vector components in spherical coordinates
    case('u_r','uphi','utheta')
       do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
          !calculate angles
          call xyz_to_sph(Xyz_DGB(x_,i,j,k,iBlock), &
               Xyz_DGB(y_,i,j,k,iBlock), &
               Xyz_DGB(z_,i,j,k,iBlock), r, theta, phi)
          sinphi = sin(phi)
          cosphi = cos(phi)
          costheta = cos(theta)
          sintheta = sin(theta)

          ! tranform to non rotating frame if needed
          U_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)/ &
               State_VGB(Rho_,i,j,k,iBlock)
          if(UseRotatingFrame) then
             U_D(x_) = U_D(x_) - OmegaBody*Xyz_DGB(y_,i,j,k,iBlock)
             U_D(y_) = U_D(y_) + OmegaBody*Xyz_DGB(x_,i,j,k,iBlock)
          end if
          U_D = 1e-3*No2Si_V(UnitU_)*U_D

          select case(NameVar)
          case('u_r')
             PlotVar_G(i,j,k) = U_D(x_)*sintheta*cosphi + &
                  U_D(y_)*sintheta*sinphi + &
                  U_D(z_)*costheta

          case('utheta')
             PlotVar_G(i,j,k) = U_D(x_)*costheta*cosphi + &
                  U_D(y_)*costheta*sinphi - &
                  U_D(z_)*sintheta

          case('uphi')
             PlotVar_G(i,j,k) =-U_D(x_)*sinphi + U_D(y_)*cosphi
          end select
       end do; end do; end do 
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'

    case('b_r','bphi','btheta')
       do k = MinK,MaxK; do j = MinJ,MaxJ; do i = MinI,MaxI
          call xyz_to_sph(Xyz_DGB(x_,i,j,k,iBlock), &
               Xyz_DGB(y_,i,j,k,iBlock), &
               Xyz_DGB(z_,i,j,k,iBlock), r, theta, phi)
          sinphi = sin(phi)
          cosphi = cos(phi)
          costheta = cos(theta)
          sintheta = sin(theta)

          if(UseB0) then
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + &
                  B0_DGB(:,i,j,k,iBlock)
          else
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if
          B_D = 1e4*No2Si_V(UnitB_)*B_D

          select case(NameVar)
          case('b_r')
             PlotVar_G(i,j,k) = B_D(x_)*sintheta*cosphi + &
                  B_D(y_)*sintheta*sinphi + &
                  B_D(z_)*costheta

          case('btheta')
             PlotVar_G(i,j,k) = B_D(x_)*costheta*cosphi + &
                  B_D(y_)*costheta*sinphi - &
                  B_D(z_)*sintheta

          case('bphi')
             PlotVAr_G(i,j,k) =-B_D(x_)*sinphi + B_D(y_)*cosphi
          end select
       end do; end do; end do 
       NameIdlUnit = 'G'
       NameTecUnit = '[G]'

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
  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    ! Fill ghost cells inside body for spherical grid - this subroutine only 
    ! modifies ghost cells in the r direction

    use EEE_ModCommonVariables, ONLY: UseCme
    use EEE_ModMain,            ONLY: EEE_get_state_BC
    use ModAdvance,    ONLY: State_VGB, B0_DGB, UseElectronPressure
    use ModGeometry,   ONLY: TypeGeometry, Xyz_DGB, r_BLK
    use ModVarIndexes, ONLY: Rho_, p_, Pe_, Bx_, Bz_
    use ModMultiFluid, ONLY: MassIon_I
    use ModImplicit,   ONLY: StateSemi_VGB, iTeImpl
    use ModPhysics,    ONLY: AverageIonCharge, UnitRho_, UnitB_, UnitP_, &
         UnitU_, Si2No_V
    use ModMain,       ONLY: n_step, iteration_number, time_simulation

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i, j, k
    real    :: Br1_D(3), Bt1_D(3)
    real    :: NumDensIon, NumDensElectron
    real    :: Runit_D(3)
    real    :: RhoCme, Ucme_D(3), Bcme_D(3), pCme, BrCme, BrCme_D(3)

    character (len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !--------------------------------------------------------------------------
    if(iSide /= 1 .or. TypeGeometry(1:9) /='spherical') &
         call CON_stop('Wrong iSide in user_set_cell_boundary')

    IsFound = .true.

    if(TypeBc == 'usersemi')then
       StateSemi_VGB(iTeImpl,0,:,:,iBlock) = tChromo
       RETURN
    elseif(TypeBc == 'usersemilinear')then
       RETURN
    end if

    ! The electron heat conduction requires the electron temperature
    ! in the ghost cells
    NumDensIon = RhoChromo/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = -1, 0
       State_VGB(Rho_,i,j,k,iBlock) = RhoChromo
       if(UseElectronPressure)then
          State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*tChromo
       else
          State_VGB(p_,i,j,k,iBlock) = (NumDensIon + NumDensElectron)*tChromo
       end if
    end do; end do; end do

    ! The following is only needed for the semi-implicit heat conduction,
    ! which averages the cell centered heat conduction coefficient towards
    ! the face
    do k = MinK,MaxK; do j = MinJ,MaxJ
       Runit_D = Xyz_DGB(:,1,j,k,iBlock) / r_BLK(1,j,k,iBlock)

       Br1_D = sum(State_VGB(Bx_:Bz_,1,j,k,iBlock)*Runit_D)*Runit_D
       Bt1_D = State_VGB(Bx_:Bz_,1,j,k,iBlock) - Br1_D

       do i = -1, 0
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = Bt1_D
       end do

    end do; end do

    if(UseCme)then
       do k = MinK, MaxK; do j = MinJ, MaxJ
          Runit_D = Xyz_DGB(:,1,j,k,iBlock) / r_BLK(1,j,k,iBlock)

          call EEE_get_state_BC(Runit_D, RhoCme, Ucme_D, Bcme_D, pCme, &
               time_simulation, n_step, iteration_number)

          RhoCme = RhoCme*Si2No_V(UnitRho_)
          Bcme_D = Bcme_D*Si2No_V(UnitB_)
          pCme   = pCme*Si2No_V(UnitP_)

          BrCme   = sum(Runit_D*Bcme_D)
          BrCme_D = BrCme*Runit_D

          do i = -1, 0
             State_VGB(Rho_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)+RhoCme
             if(UseElectronPressure)then
                State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
                     + 0.5*pCme
             else
                State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) + pCme
             end if
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock) + BrCme_D
          end do
       end do; end do
    end if

  end subroutine user_set_cell_boundary
  !============================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use EEE_ModCommonVariables, ONLY: UseCme
    use EEE_ModMain,            ONLY: EEE_get_state_BC
    use ModAdvance,      ONLY: UseElectronPressure, UseAnisoPressure
    use ModFaceBoundary, ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D, TimeBc
    use ModMain,         ONLY: x_, y_, UseRotatingFrame, n_step, &
         iteration_number
    use ModMultiFluid,   ONLY: MassIon_I
    use ModPhysics,      ONLY: OmegaBody, AverageIonCharge, UnitRho_, &
         UnitP_, UnitB_, UnitU_, Si2No_V, inv_gm1
    use ModVarIndexes,   ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, Bz_, p_, &
         WaveFirst_, WaveLast_, Pe_, Ppar_, Hyp_, Ehot_
    use ModConst,        ONLY: cProtonMass
    use ModHeatFluxCollisionless, ONLY: UseHeatFluxCollisionless, &
         get_gamma_collisionless

    real, intent(out) :: VarsGhostFace_V(nVar)

    integer :: iP
    real :: NumDensIon, NumDensElectron, FullBr, Ewave, Pressure, Temperature
    real :: Gamma
    real,dimension(3) :: U_D, B1_D, B1t_D, B1r_D, rUnit_D

    ! Line-tied related variables
    real              :: RhoTrue, RhoGhost
    real,dimension(3) :: bUnitGhost_D, bUnitTrue_D
    real,dimension(3) :: FullBGhost_D, FullBTrue_D

    ! CME related variables
    real :: RhoCme, Ucme_D(3), Bcme_D(3), pCme
    real :: BrCme, BrCme_D(3), UrCme, UrCme_D(3), UtCme_D(3)
   
    character (len=*), parameter :: NameSub = 'user_set_face_boundary'
    !--------------------------------------------------------------------------

    rUnit_D = FaceCoords_D/sqrt(sum(FaceCoords_D**2))

    B1_D  = VarsTrueFace_V(Bx_:Bz_)
    B1r_D = sum(rUnit_D*B1_D)*rUnit_D
    B1t_D = B1_D - B1r_D
    VarsGhostFace_V(Bx_:Bz_) = B1t_D

    ! Fix density
    VarsGhostFace_V(Rho_) = RhoChromo

    if (UseUparBc) then
       ! Use line-tied boundary conditions
       U_D = VarsTrueFace_V(Ux_:Uz_)

       RhoTrue = VarsTrueFace_V(Rho_)
       RhoGhost = VarsGhostFace_V(Rho_)
       FullBGhost_D = B0Face_D + VarsGhostFace_V(Bx_:Bz_)
       FullBTrue_D  = B0Face_D + VarsTrueFace_V(Bx_:Bz_)
 
       bUnitGhost_D = FullBGhost_D/sqrt(max(1e-30,sum(FullBGhost_D**2)))
       bUnitTrue_D = FullBTrue_D/sqrt(max(1e-30,sum(FullBTrue_D**2)))

       ! extrapolate field-aligned velocity component
       VarsGhostFace_V(Ux_:Uz_) = RhoTrue/RhoGhost* &
            sum(U_D*bUnitTrue_D)*bUnitGhost_D
    else
       ! zero velocity at inner boundary
       VarsGhostFace_V(Ux_:Uz_) = -VarsTrueFace_V(Ux_:Uz_)
    end if

    ! Apply corotation if needed
    if(.not.UseRotatingFrame)then                                            
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_)-2*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_)+2*OmegaBody*FaceCoords_D(x_)
    end if   

    Temperature = tChromo

    if(UseCme)then
       call EEE_get_state_BC(Runit_D, RhoCme, Ucme_D, Bcme_D, pCme, TimeBc, &
            n_step, iteration_number)

       RhoCme = RhoCme*Si2No_V(UnitRho_)
       Ucme_D = Ucme_D*Si2No_V(UnitU_)
       Bcme_D = Bcme_D*Si2No_V(UnitB_)
       pCme   = pCme*Si2No_V(UnitP_)

       ! Add CME density
       VarsGhostFace_V(Rho_) = VarsGhostFace_V(Rho_) + RhoCme

       ! Fix the normal component of the CME field to BrCme_D at the Sun
       BrCme   = sum(Runit_D*Bcme_D)
       BrCme_D = BrCme*Runit_D
       VarsGhostFace_V(Bx_:Bz_) = VarsGhostFace_V(Bx_:Bz_) + BrCme_D

       ! Fix the tangential components of the CME velocity at the Sun
       UrCme   = sum(Runit_D*Ucme_D)
       UrCme_D = UrCme*Runit_D
       UtCme_D = UCme_D - UrCme_D
       VarsGhostFace_V(Ux_:Uz_) = VarsGhostFace_V(Ux_:Uz_) + 2*UtCme_D

       Pressure = RhoChromo/MassIon_I(1)*(1 + AverageIonCharge)*tChromo
       Temperature = (Pressure + pCme) &
            / (VarsGhostFace_V(Rho_)/MassIon_I(1)*(1 + AverageIonCharge))
    end if
 
    FullBr = sum((B0Face_D + VarsGhostFace_V(Bx_:Bz_))*rUnit_D)

    ! Ewave \propto sqrt(rho) for U << Ualfven
    Ewave = PoyntingFluxPerB*sqrt(VarsGhostFace_V(Rho_))
    if (FullBr > 0. ) then
       VarsGhostFace_V(WaveFirst_) = Ewave
       VarsGhostFace_V(WaveLast_) = 0.0
    else
       VarsGhostFace_V(WaveFirst_) = 0.0
       VarsGhostFace_V(WaveLast_) = Ewave
    end if

    ! Fix temperature
    NumDensIon = VarsGhostFace_V(Rho_)/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    if(UseElectronPressure)then
       VarsGhostFace_V(p_) = NumDensIon*Temperature
       VarsGhostFace_V(Pe_) = NumDensElectron*Temperature
       if(UseAnisoPressure) VarsGhostFace_V(Ppar_) = VarsGhostFace_V(p_)
    else
       VarsGhostFace_V(p_) = (NumDensIon + NumDensElectron)*Temperature
    end if

    if(Hyp_>1) VarsGhostFace_V(Hyp_) = VarsTrueFace_V(Hyp_)

    if(Ehot_ > 1)then
       if(UseHeatFluxCollisionless)then
          call get_gamma_collisionless(FaceCoords_D, Gamma)
          iP = p_; if(UseElectronPressure) iP = Pe_
          VarsGhostFace_V(Ehot_) = &
               VarsGhostFace_V(iP)*(1.0/(Gamma - 1) - inv_gm1)
       else
          VarsGhostFace_V(Ehot_) = 0.0
       end if
    end if

  end subroutine user_set_face_boundary
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
  subroutine user_initial_perturbation

    use ModProcMH,    ONLY: iProc
    use EEE_ModMain,  ONLY: EEE_get_state_init
    use ModMain, ONLY: nI, nJ, nK, nBLK, unused_B, n_step, iteration_number
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB, UseElectronPressure
    use ModPhysics,   ONLY: Si2No_V, UnitU_, UnitRho_, UnitP_, UnitB_
    use ModGeometry,  ONLY: Xyz_DGB
    use ModEnergy,    ONLY: calc_energy_cell
    use BATL_lib,     ONLY: nDim, MaxDim

    integer :: i, j, k, iBlock
    logical :: oktest, oktest_me
    real :: x_D(nDim), Rho, B_D(MaxDim), p

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    ! -------------------------------------------------------------------------
    call set_oktest('user_initial_perturbation',oktest,oktest_me)

    do iBlock = 1, nBLK
       if(unused_B(iBlock))CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI

          x_D = Xyz_DGB(:,i,j,k,iBlock)

          call EEE_get_state_init(x_D, Rho, B_D, p, n_step, iteration_number)

          Rho = Rho*Si2No_V(UnitRho_)
          B_D = B_D*Si2No_V(UnitB_)
          p = p*Si2No_V(UnitP_)

          !\
          ! Add the eruptive event state to the solar wind
          !/
          if(State_VGB(Rho_,i,j,k,iBlock) + Rho < 0.25*State_VGB(Rho_,i,j,k,iBlock))then
             State_VGB(Rho_,i,j,k,iBlock) = 0.25*State_VGB(Rho_,i,j,k,iBlock)
          else
             State_VGB(Rho_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) + Rho
          endif

          State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock) + B_D

          if(UseElectronPressure)then
             if(State_VGB(Pe_,i,j,k,iBlock) + 0.5*p < 0.25*State_VGB(Pe_,i,j,k,iBlock))then
                State_VGB(Pe_,i,j,k,iBlock) = 0.25*State_VGB(Pe_,i,j,k,iBlock)
             else
                State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) + 0.5*p
             endif

             if(State_VGB(p_,i,j,k,iBlock)  + 0.5*p < 0.25*State_VGB(p_,i,j,k,iBlock))then
                State_VGB(p_,i,j,k,iBlock) = 0.25*State_VGB(p_,i,j,k,iBlock)
             else
                State_VGB(p_,i,j,k,iBlock)  = State_VGB(p_,i,j,k,iBlock)  + 0.5*p
             endif

          else
             if(State_VGB(p_,i,j,k,iBlock) + p < 0.25*State_VGB(p_,i,j,k,iBlock))then
                State_VGB(p_,i,j,k,iBlock) = 0.25*State_VGB(p_,i,j,k,iBlock)
             else
                State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) + p
             endif
          end if
       end do; end do; end do

       call calc_energy_cell(iBlock)

    end do

  end subroutine user_initial_perturbation

end module ModUser
