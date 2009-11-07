!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModReadParam, ONLY: lStringLine
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_face_bcs,                   &
       IMPLEMENTED5 => user_get_log_var,                &
       IMPLEMENTED6 => user_get_b0,                     &
       IMPLEMENTED7 => user_calc_sources,               &
       IMPLEMENTED8 => user_update_states,              &
       IMPLEMENTED9 => user_specify_refinement,         &
       IMPLEMENTED10=> user_set_boundary_cells,         &
       IMPLEMENTED11=> user_set_plot_var

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Global Corona'

  real :: BcAlfvenWavePressure, BcAlfvenWavePressureCgs = 8.36e-2
  real :: HeatingAmplitude, HeatingAmplitudeCgs = 6.07e-7
  real :: DecayLength = 0.7
  logical :: UseExponentialHeating = .false.

  real :: TeFraction

  character(len=lStringLine) :: NameModel, TypeCoronalHeating

contains

  !============================================================================

  subroutine user_read_inputs

    use ModCoronalHeating, ONLY: UseUnsignedFluxModel
    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMagnetogram, ONLY: set_parameters_magnetogram

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
       case("#PFSSM")
          call read_var('UseUserB0', UseUserB0)
          if(UseUserB0)then
             call set_parameters_magnetogram
             call read_var('dt_UpdateB0', dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          end if

       case("#ALFVENWAVE")
          call read_var('BcAlfvenWavePressureCgs', BcAlfvenWavePressureCgs)

       case("#CORONALHEATING")
          call read_var('TypeCoronalHeating', TypeCoronalHeating)
          select case(TypeCoronalHeating)
          case('exponential')
             UseUnsignedFluxModel = .false.
             UseExponentialHeating = .true.
             call read_var('DecayLength', DecayLength)
             call read_var('HeatingAmplitudeCgs', HeatingAmplitudeCgs)
          case('unsignedflux')
             UseUnsignedFluxModel = .true.
             UseExponentialHeating = .false.
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
    use ModIO,          ONLY: write_prefix, iUnitOut,NamePlotDir
    use ModMagnetogram, ONLY: read_magnetogram_file
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: Si2No_V, UnitP_, UnitEnergyDens_, UnitT_, &
         ElectronTemperatureRatio, AverageIonCharge
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
    use ModWaves,       ONLY: UseWavePressure, UseAlfvenSpeed

    !--------------------------------------------------------------------------
    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session:'
       call write_prefix; write(iUnitOut,*) ''
    end if

    if(i_line_command("#PFSSM", iSessionIn = 1) < 0)then
       write(*,*) 'In session 1, a magnetogram file has to be read via #PFSSM'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if
    if(i_line_command("#PFSSM") > 0)then
       call read_magnetogram_file(NamePlotDir)
    end if

    if(TypeCoronalHeating == 'exponential')then
        HeatingAmplitude =  HeatingAmplitudeCgs*0.1 &
             *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
     end if

    if(BcAlfvenWavePressureCgs < 0.0)then
       write(*,*) 'Negative Alfven wave pressure not allowed'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if
    BcAlfvenWavePressure = BcAlfvenWavePressureCgs*0.1*Si2No_V(UnitP_)

    UseAlfvenSpeed = .true.
    UseWavePressure = .true.

    ! TeFraction is used for ideal EOS:
    if(UseElectronPressure)then
       ! Pe = ne*Te (dimensionless) and n=rho/ionmass
       ! so that Pe = ne/n *n*Te = (ne/n)*(rho/ionmass)*Te
       ! TeFraction is defined such that Te = Pe/rho * TeFraction
       TeFraction = MassIon_I(1)/AverageIonCharge
    else
       ! p = n*T + ne*Te (dimensionless) and n=rho/ionmass
       ! so that p=rho/massion *T*(1+ne/n Te/T)
       ! TeFraction is defined such that Te = p/rho * TeFraction
       TeFraction = MassIon_I(1)*ElectronTemperatureRatio &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
    end if

    if(iProc == 0)then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session

  !============================================================================

  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance,     ONLY: State_VGB, UseElectronPressure
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: OmegaBody, BodyRho_I, BodyTDim_I, &
         UnitTemperature_, Si2No_V, AverageIonCharge
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_, &
         WaveFirst_, WaveLast_, Pe_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: Density, NumDensIon, NumDensElectron, Tbase, FullBr
    real :: Runit_D(3), U_D(3)
    real :: B1_D(3), B1t_D(3), B1r_D(3), FullB_D(3)
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

    if(FullBr > 0.0)then
       VarsGhostFace_V(WaveFirst_) = 2.0*BcAlfvenWavePressure
       VarsGhostFace_V(WaveLast_) = VarsTrueFace_V(WaveLast_)
    else
       VarsGhostFace_V(WaveFirst_) = VarsTrueFace_V(WaveFirst_)
       VarsGhostFace_V(WaveLast_) = 2.0*BcAlfvenWavePressure
    end if

    Density = BodyRho_I(1)
    ! The electron and ion temperature are equal to Tbase at the coronal base
    Tbase = BodyTDim_I(1)*Si2No_V(UnitTemperature_)

!!!    VarsGhostFace_V(Rho_) = &
!!!         max(-VarsTrueFace_V(Rho_) + 2.0*Density, VarsTrueFace_V(Rho_))
!!!    VarsGhostFace_V(p_) = &
!!!         max(VarsGhostFace_V(Rho_)*Temperature, VarsTrueFace_V(p_))

    VarsGhostFace_V(Rho_) =  2.0*Density - VarsTrueFace_V(Rho_)
    NumDensIon = VarsGhostFace_V(Rho_)/MassIon_I(1)
    NumDensElectron = NumDensIon*AverageIonCharge
    if(UseElectronPressure)then
       VarsGhostFace_V(p_) = NumDensIon*Tbase
       VarsGhostFace_V(Pe_) = NumDensElectron*Tbase
    else
       VarsGhostFace_V(p_) = (NumDensIon + NumDensElectron)*Tbase
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

  subroutine user_set_ics

    ! The isothermal parker wind solution is used as initial condition

    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk, r_Blk, true_cell
    use ModMain,       ONLY: nI, nJ, nK, globalBLK
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         BodyRho_I, BodyTDim_I, No2Si_V, UnitU_, AverageIonCharge
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_, Pe_, &
         WaveFirst_, WaveLast_

    integer :: i, j, k, iBlock
    integer :: IterCount
    real :: x, y, z, r, RhoBase, Rho, NumDensIon, NumDensElectron
    real :: Tbase, Tcorona, Temperature
    real :: Ur, Ur0, Ur1, del, Ubase, rTransonic, Uescape, Usound

    real, parameter :: Epsilon = 1.0e-6
    !--------------------------------------------------------------------------

    iBlock = globalBLK

    RhoBase = BodyRho_I(1)
    ! The electron and ion temperature are equal to Tbase at the coronal base
    Tbase = BodyTDim_I(1)*Si2No_V(UnitTemperature_)
    ! Initially, the electron and ion temperature are at 1.5e6(K) in the corona
    Tcorona = 1.5e6*Si2No_V(UnitTemperature_)

    Usound = sqrt(Tcorona*(1.0+AverageIonCharge)/MassIon_I(1))
    Uescape = sqrt(-GBody*2.0/rBody)/Usound

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
       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 0.0
       NumDensIon = Rho/MassIon_I(1)
       NumDensElectron = NumDensIon*AverageIonCharge
       if(true_cell(i,j,k,iBlock))then
          Temperature = Tcorona
       else
          Temperature = Tbase
       end if
       if(UseElectronPressure)then
          State_VGB(p_,i,j,k,iBlock) = NumDensIon*Temperature
          State_VGB(Pe_,i,j,k,iBlock) = NumDensElectron*Temperature
       else
          State_VGB(p_,i,j,k,iBlock) = &
               (NumDensIon + NumDensElectron)*Temperature
       end if
    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_get_b0(x, y, z, B0_D)

    use ModPhysics,     ONLY: Si2No_V, UnitB_
    use ModMagnetogram, ONLY: get_magnetogram_field

    real, intent(in) :: x, y, z
    real, intent(out):: B0_D(3)
    !--------------------------------------------------------------------------

    call get_magnetogram_field(x, y, z, B0_D)
    B0_D = B0_D*Si2No_V(UnitB_)

  end subroutine user_get_b0

  !============================================================================

  subroutine user_calc_sources

    use ModAdvance,        ONLY: State_VGB, Source_VC, UseElectronPressure
    use ModCoronalHeating, ONLY: UseUnsignedFluxModel, get_coronal_heating
    use ModGeometry,       ONLY: r_BLK
    use ModMain,           ONLY: nI, nJ, nK, GlobalBlk
    use ModPhysics,        ONLY: gm1
    use ModVarIndexes,     ONLY: Rho_, Energy_, Pe_

    integer :: i, j, k, iBlock
    real :: CoronalHeating, RadiativeCooling

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------

    iBlock = globalBlk

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(UseUnsignedFluxModel)then
          call get_coronal_heating(i, j, k, iBlock, CoronalHeating)
       elseif(UseExponentialHeating)then
          CoronalHeating = HeatingAmplitude &
               *exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLength)
       else
          CoronalHeating = 0.0
       end if

       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + CoronalHeating

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

    real :: Te, TeSi, TeFraction, CoolingFunctionCgs
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

  subroutine user_update_states(iStage, iBlock)

    integer, intent(in) :: iStage, iBlock
    !--------------------------------------------------------------------------

    call update_states_MHD(iStage, iBlock)

  end subroutine user_update_states

  !============================================================================

  subroutine user_get_log_var(VarValue,TypeVar,Radius)

    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK, nBlock, x_, y_, z_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_, p_ 
    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, B0_DGB
    use ModPhysics,    ONLY: inv_gm1, No2Io_V, UnitEnergydens_, UnitX_

    real, intent(out) :: VarValue
    character (LEN=10), intent(in) :: TypeVar 
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
          tmp1_BLK(:,:,:,iBlock) = State_VGB(P_,:,:,:,iBlock)
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
    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: AverageIonCharge, ElectronTemperatureRatio, &
         No2Si_V, UnitTemperature_
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
    real :: TiFraction

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
       TiFraction = MassIon_I(1) &
            /(1 + AverageIonCharge*ElectronTemperatureRatio)
       do k = -1, nK+2; do j = -1, nJ+2; do i = -1, nI+2
          if(UseElectronPressure)then
             PlotVar_G(i,j,k) = MassIon_I(1)*State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
          else
             PlotVar_G(i,j,k) = TiFraction*State_VGB(p_,i,j,k,iBlock) &
                  /State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
          end if
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

  subroutine user_set_boundary_cells(iBLK)

    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, r_Blk
    use ModBoundaryCells, ONLY: SaveBoundaryCells
    use ModPhysics,       ONLY: rBody

    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = r_Blk(:,:,:,iBLK) < rBody

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

end module ModUser

