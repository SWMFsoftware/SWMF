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
       IMPLEMENTED11=> user_set_outerbcs


  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Low Corona / TransRegion'
  ! quantitative parameters for exponential heating model
  real :: HeatingAmplitude, HeatingAmplitudeCgs = 6.07e-7
  real :: DecayLengthExp = 0.7  ! in Solar Radii units
  logical :: UseExponentialHeating = .false.

  ! ratio of electron Temperature / Total Temperature
  ! given by the formula P = ne * k_B*Te + n_i k_B*T_i
  real :: TeFraction

  character(len=lStringLine) :: NameModel, TypeCoronalHeating

  ! additional variables for TR boundary / heating models
  character(len=lStringLine) :: TypeTRBoundary
  
  logical :: DoChromoBC  = .false.
  logical :: DoCoronalBC = .false.
  logical :: DoREBModel  = .false.

  ! boundary condition variables to be read from PARAM.in
  real :: BoundaryTeSi = 1.5E+6
  real :: BoundaryNeCgs = 3.0E+8
  ! dimensionless values of boundary electron temp Te and mass density
  real :: BoundaryTe
  real :: BoundaryRho

  integer :: iTableRadCool = -1

  ! Same use as the variables in ModHeatConduction
  ! but need them here as well
  logical :: DoModHeatConduction = .false.
  real :: TeModSi = 3.0E+5
  real :: DeltaTeModSi = 1E+4
  real :: TeMod, DeltaTeMod

contains

  !============================================================================

  subroutine user_read_inputs

    use ModCoronalHeating, ONLY: UseUnsignedFluxModel
    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMagnetogram, ONLY: set_parameters_magnetogram
    use ModPhysics,     ONLY: SW_T_dim, SW_n_dim

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

!       case("#ALFVENWAVE")
!          call read_var('BcAlfvenWavePressureCgs', BcAlfvenWavePressureCgs)

       case("#CORONALHEATING")
          call read_var('TypeCoronalHeating', TypeCoronalHeating)
          select case(TypeCoronalHeating)
          case('exponential')
             UseUnsignedFluxModel = .false.
             UseExponentialHeating = .true.
             call read_var('DecayLengthExp', DecayLengthExp)
             call read_var('HeatingAmplitudeCgs', HeatingAmplitudeCgs)
          case('unsignedflux')
             UseUnsignedFluxModel = .true.
             UseExponentialHeating = .false.
          case default
             call stop_mpi(NameSub//': unknown TypeCoronalHeating = ' &
                  //TypeCoronalHeating)
       end select

       ! for compatibility need to read in modified heat values twice:
       ! Once here and once in the #PARALLELCONDUCTION call.
       ! Problem is ModHeatConduction needs user_material_properties
       ! so can't have ModUser needing variables from ModHeatConduction
       case("#MODIFYHEAT")
          call read_var('DoModHeatConduction',DoModHeatConduction)
          call read_var('TeModSi',TeModSi)
          call read_var('DeltaTeModSi',DeltaTeModSi)

       case("#TRBOUNDARY")
          call read_var('TypeTRBoundary', TypeTRBoundary)
          call read_var('BoundaryNeCgs',BoundaryNeCgs)
          call read_var('BoundaryTeSi',BoundaryTeSi)
          select case(TypeTRBoundary)
          case('chromo')
             DoChromoBC  = .true.
             DoCoronalBC = .false.
             DoREBModel  = .false.
          case('coronal')
             DoChromoBC  = .false.
             DoCoronalBC = .true.
             DoREBModel  = .false.
          case('reb','REB')
             DoChromoBC  = .false.
             DoCoronalBC = .false.
             DoREBModel  = .true.
          case default
             call stop_mpi(NameSub//': unknown TypeTRBoundary = ' &
                  //TypeTRBoundary)
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

    use ModAdvance,     ONLY: WaveFirst_, WaveLast_
    use ModIO,          ONLY: write_prefix, iUnitOut,NamePlotDir
    use ModMagnetogram, ONLY: read_magnetogram_file
    use ModMultiFluid,  ONLY: MassIon_I
    use ModPhysics,     ONLY: Si2No_V, UnitP_, UnitEnergyDens_, UnitT_, &
         ElectronTemperatureRatio, AverageIonCharge, UnitTemperature_, &
         UnitRho_, No2Si_V, UnitN_
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
    use ModWaves
    use ModMain,        ONLY: optimize_message_pass
    use ModLookupTable, ONLY: i_lookup_table

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

    if(i_line_command("#TRBOUNDARY", iSessionIn = 1) < 0)then
       write(*,*) 'In session 1, need to specify a BC with #TRBOUNDARY!'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if

    if(optimize_message_pass/='all') then
       write(*,*) 'For Heat Conduction need message pass = all with ',&
                   '#MESSAGEPASS!'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if

    if(TypeCoronalHeating == 'exponential')then
        HeatingAmplitude =  HeatingAmplitudeCgs*0.1 &
             *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
    end if

    iTableRadCool = i_lookup_table('radcool')

    if(i_line_command("#PLASMA", iSessionIn = 1) < 0)then
       write(*,*) 'Need to set electron temp ration with #PLASMA!'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if

    ! ratio of Te / (Te + Tp)
    ! given by the formula P = ne * k_B*Te + n_i k_B*T_i
    TeFraction = MassIon_I(1)*ElectronTemperatureRatio &
         /(1 + AverageIonCharge*ElectronTemperatureRatio)

    if(DoModHeatConduction) then
       TeMod = TeModSi * Si2No_V(UnitTemperature_)
       DeltaTeMod = DeltaTeModSi * Si2No_V(UnitTemperature_)
    end if

    ! calc normalized values of BC Te and Ne
    ! note, implicitly assuming Ne = Ni here
    BoundaryTe = BoundaryTeSi * Si2No_V(UnitTemperature_)
    BoundaryRho = BoundaryNeCgs * 1.0E+6 * Si2No_V(UnitN_)

    if(iProc == 0)then
       if(iTableRadCool>0) then
          call write_prefix;  write(iUnitOut,*) 'Using Tabulated RadCooling'
       else
          call write_prefix;  write(iUnitOut,*) 'Using Analytic Fit to ', &
                                 'RadCooling'
       end if
       
       call write_prefix; write(iUnitOut,*) 'TeFraction = ',TeFraction
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'user_init_session finished'
       call write_prefix; write(iUnitOut,*) ''
    end if

  end subroutine user_init_session

  !============================================================================

  subroutine user_face_bcs(VarsGhostFace_V)

    use ModAdvance,     ONLY: State_VGB, WaveFirst_, WaveLast_
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame
    use ModNumConst,    ONLY: cTolerance
    use ModPhysics,     ONLY: OmegaBody, BodyRho_I, BodyTDim_I, &
         UnitTemperature_, Si2No_V
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: Density, Temperature, FullBr
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

    

    Density = BoundaryRho
!    WILL RE-Implement REB MODEL SOON!
!    if(DoREBModel) call calc_REB_density
    Temperature = BoundaryTe/TeFraction


    VarsGhostFace_V(Rho_) =  2.0*Density - VarsTrueFace_V(Rho_)
    VarsGhostFace_V(p_) = VarsGhostFace_V(Rho_)*Temperature

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

    use ModAdvance,    ONLY: State_VGB, WaveFirst_, WaveLast_
    use ModGeometry,   ONLY: x_Blk, y_Blk, z_Blk, r_Blk, true_cell
    use ModMain,       ONLY: nI, nJ, nK, globalBLK
    use ModPhysics,    ONLY: Si2No_V, UnitTemperature_, rBody, GBody, &
         BodyRho_I, BodyTDim_I, No2Si_V, UnitU_, UnitN_
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, Bz_, p_

    integer :: i, j, k, iBlock
   
    ! Variables for IC using Parker Isothermal Wind Solution
    ! Parker part copied from ModUserScHeat.f90 on 9/20/2009
    integer :: IterCount
    real :: x, y, z, r, RhoBase, Tbase, Rho, rBase
    real :: Ur, Ur0, Ur1, del, Ubase, rTransonic, Uescape, Usound, T0
    real, parameter :: Epsilon = 1.0e-6

    ! variables for adding in TR like icond if needed
    ! this is from 2 lines that peicewise fit an old relaxed solution
    ! not really physical! just gets you close
    ! also, if solutions change significantly, will need to update
    real, parameter :: rTransRegion = 1.03       ! widened w/ modify heatflux approximation
    real, parameter :: PowerTe = 341.81783
    real, parameter :: PowerNe = -745.56192
    real, parameter :: NeEmpirical = 6.0E+8
    real, parameter :: TeEmpirical = 6.0E+5
    real, parameter :: rEmpirical = 1.01
    ! TeBaseCGS and NeBaseCGS and rBase are values that match the coronal
    ! Solution at rTransRegion
    real :: TeBaseCgs, NeBaseCgs, NeValue, TeValue
     
    !--------------------------------------------------------------------------
    ! This icond implements the radially symmetric parker isothermal wind solution 
    ! for the corona. If using a chromospheric condition, extend base of wind solution
    ! to r=rTransRegion, and put in approximate transition region below 
    ! 
    ! NOTE this coronal icond does NOT use T value from the PARAM.in
    ! #BODY command, instead it sets the coronal Te to 1.5 MK everywhere.
    ! It DOES however, set the base density of the coronal part according to 
    ! BodyRho_I(1)

    iBlock = globalBLK

    T0 = 3.0e6*Si2No_V(UnitTemperature_)
    if(DoChromoBC) then 
        rBase = rTransRegion
    else
        rBase = rBody
    endif
    RhoBase = rBase**2 * BodyRho_I(1)
    NeBaseCGS = RhoBase * No2Si_V(UnitN_) * 1.0E-6
    TeBaseCGS = TeFraction * T0 * No2Si_V(UnitTemperature_)

    Usound = sqrt(T0)
    Uescape = sqrt(-GBody*2.0/rBase)/Usound

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
       r = max(r_BLK(i,j,k,iBlock),1.0)
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       if (r >= rBase) then 

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
   
          Rho = RhoBase*Ubase/(r**2*Ur)
          State_VGB(Rho_,i,j,k,iBlock) = Rho
          State_VGB(RhoUx_,i,j,k,iBlock) = Rho*Ur*x/r *Usound
          State_VGB(RhoUy_,i,j,k,iBlock) = Rho*Ur*y/r *Usound
          State_VGB(RhoUz_,i,j,k,iBlock) = Rho*Ur*z/r *Usound
          State_VGB(P_,i,j,k,iBlock) = Rho*T0

       else ! below parker part, apply TR fix

          if (r <= rEmpirical) then
             ! power law approx to what I saw in early run from r=1.00 to r=1.01
             ! exponents values are ridiculous because of tiny range in r
             ! DO NOT APPLY TO OTHER LIMITS!
             NeValue = BoundaryNeCgs*(r/1.0)**PowerNe
             TeValue = BoundaryTeSi*(r/1.0)**PowerTe
          else
             ! now do linear span from rEmpirical to rTransregion
             NeValue = (NeBaseCgs - NeEmpirical) / (rTransRegion - rEmpirical) &
                        * (r-rEmpirical) + NeEmpirical
             TeValue = (TeBaseCgs - TeEmpirical) / (rTransRegion - rEmpirical) &
                        * (r-rEmpirical) + TeEmpirical
          endif

          State_VGB(Rho_,i,j,k,iBlock) = NeValue * Si2No_V(UnitN_) * 1.0E+6
          State_VGB(P_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) * TeValue &
                                        * Si2No_V(UnitTemperature_) / TeFraction                   
          State_VGB(RhoUx_,i,j,k,iBlock) = 0.0
          State_VGB(RhoUy_,i,j,k,iBlock) = 0.0
          State_VGB(RhoUz_,i,j,k,iBlock) = 0.0

       end if
         
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
!       State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = 0.0

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

    use ModAdvance,        ONLY: State_VGB, Source_VC, Rho_, p_, Energy_, &
         UseNonConservative
    use ModCoronalHeating, ONLY: UseUnsignedFluxModel, get_coronal_heating
    use ModGeometry,       ONLY: r_BLK
    use ModMain,           ONLY: nI, nJ, nK, GlobalBlk
    use ModPhysics,        ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         inv_gm1

    integer :: i, j, k, iBlock
    real :: CoronalHeating, RadiativeCooling, EinternalSource, Cv

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------

    iBlock = globalBlk

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(UseUnsignedFluxModel)then
          call get_coronal_heating(i, j, k, iBlock, CoronalHeating)
       elseif(UseExponentialHeating)then
          CoronalHeating = HeatingAmplitude &
               *exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLengthExp)
       else
          CoronalHeating = 0.0
       end if

       call get_radiative_cooling( &
            State_VGB(:,i,j,k,iBlock), RadiativeCooling)

       EinternalSource = CoronalHeating + RadiativeCooling

       if(UseNonConservative)then
          Cv = inv_gm1*State_VGB(Rho_,i,j,k,iBlock)/TeFraction

          Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
               + EinternalSource*State_VGB(Rho_,i,j,k,iBlock)/Cv
       else
          Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + EinternalSource
       end if
    end do; end do; end do

  end subroutine user_calc_sources

  !============================================================================

  subroutine get_radiative_cooling(State_V, RadiativeCooling)

    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitT_, UnitN_, &
         UnitEnergyDens_, UnitTemperature_
    use ModVarIndexes, ONLY: nVar, Rho_, p_
    use ModLookupTable, ONLY: interpolate_lookup_table

    real, intent(in) :: State_V(nVar)
    real, intent(out):: RadiativeCooling

    real :: Te, TeSi, CoolingFunctionCgs, NumberDensCgs
    real :: Log10TeSi, Log10NeCgs, CoolingTableOut(2)
    real, parameter :: RadNorm = 1.0E+22
    real, parameter :: TeModMinSi = 2.0E+4
    real :: TeFactor, TeModMin, FractionSpitzer
    !--------------------------------------------------------------------------

    Te = TeFraction*State_V(p_)/State_V(Rho_)
    TeSi =Te*No2Si_V(UnitTemperature_)

    ! currently proton plasma only
    NumberDensCgs = State_V(Rho_)*No2Si_V(UnitN_)*1.0e-6
    if(iTableRadCool>0) then
       Log10TeSi = log10(TeSi)
       Log10NeCgs = log10(NumberDensCgs)
       ! at the moment, radC not a function of Ne, but need a dummy 2nd
       ! index, and might want to include Ne dependence in table later.
       ! Also lookuptables need 2 output values, but we need only 1 in this case
       ! so use dummy variable in 2nd var column.
       ! *** also table variable should be normalized to radloss_cgs * 10E+22
       ! since don't want to deal with such tiny numbers 
       call interpolate_lookup_table(iTableRadCool,Log10TeSi,Log10NeCgs, &
             CoolingTableOut, DoExtrapolate = .true.)
       CoolingFunctionCgs = CoolingTableOut(1) / RadNorm
    else
       call get_cooling_function_fit(TeSi, CoolingFunctionCgs)
    end if

    RadiativeCooling = -NumberDensCgs**2*CoolingFunctionCgs &
         *0.1*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

    ! include multiplicative factors to make up for extention of
    ! perpendicular heating at low temperatures (as per Abbet 2007).
    ! Need this to strech transition region to larger scales
    ! Also, need radcool modification to become const below TeModMin
    if(DoModHeatConduction) then
       TeModMin = TeModMinSi * Si2No_V(UnitTemperature_)
       FractionSpitzer = 0.5*(1.0+tanh((Te-TeMod)/DeltaTeMod))
       TeFactor = (FractionSpitzer * Te**2.5 &
                  + (1.0 - FractionSpitzer)*TeMod**2.5)
       if(Te >= TeModMin) then
          RadiativeCooling = RadiativeCooling * Te**2.5 / TeFactor
       else
          RadiativeCooling = RadiativeCooling * (TeModMin / TeMod)**2.5
       endif
    end if

  contains

    subroutine get_cooling_function_fit(TeSi, CoolingFunctionCgs)

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

    end subroutine get_cooling_function_fit

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
  
  !=====================================================================
  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)
    
    use ModAdvance,  ONLY: Rho_, P_, State_VGB
    use ModGeometry, ONLY: TypeGeometry
    use ModSize,     ONLY: East_

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: NameSub = 'user_set_outerbcs'
    !-------------------------------------------------------------------
    
    ! This routine used only for setting the inner r ghost cells for
    ! spherical geometry. Need to fix the temperature to the boundary
    ! temperature (which is NOT necessarily the BODY normalization values)
    ! for the heat conduction calculation. The face_gradient calculation
    ! uses ghost cells! If Face gradient was checking values other than
    ! P/rho, would need to set those as well!

    if(iSide==East_) then
       State_VGB(Rho_,-1:0,:,:,iBlock) = BoundaryRho
       State_VGB(P_  ,-1:0,:,:,iBlock) = BoundaryRho * BoundaryTe/TeFraction
    else
       call stop_mpi('For TR Model ONLY East_ (low R) user boundary can be used')
    endif

    IsFound = .true.
  end subroutine user_set_outerbcs

end module ModUser

