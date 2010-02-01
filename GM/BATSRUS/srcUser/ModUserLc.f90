!^CFG COPYRIGHT UM
!==============================================================================
module ModUser
  use ModMain,      ONLY: nBLK, nI, nJ, nK
  use ModReadParam, ONLY: lStringLine
  use ModCoronalHeating,ONLY: CoronalHeating_C,get_cell_heating
  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_init_session,               &
       IMPLEMENTED3 => user_set_ics,                    &
       IMPLEMENTED4 => user_initial_perturbation,       &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_get_log_var,                &
       IMPLEMENTED7 => user_get_b0,                     &
       IMPLEMENTED8 => user_calc_sources,               &
       IMPLEMENTED9 => user_update_states,              &
       IMPLEMENTED10=> user_specify_refinement,         &
       IMPLEMENTED11=> user_set_boundary_cells,         &
       IMPLEMENTED12=> user_set_outerbcs,               &
       IMPLEMENTED13=> user_set_plot_var,               &
       IMPLEMENTED14=> user_material_properties


  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Low Corona / TransRegion'


  ! ratio of electron Temperature / Total Temperature
  ! given by the formula P = ne * k_B*Te + n_i k_B*T_i
  real :: TeFraction


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
  real :: HeatCondParSi, HeatCondPar

  ! Variables for the REB model
  logical :: IsNewBlockTeCalc(nBLK) = .true.
  ! cell centered electron temperature for entire block
  ! put here so not always re-computed during boundary calculation
  real :: Te_G(-1:nI+2,-1:nJ+2,-1:nK+2) 

contains

  !============================================================================

  subroutine user_read_inputs

    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: read_line, read_command, read_var
    use ModIO,          ONLY: write_prefix, write_myname, iUnitOut
    use ModMagnetogram, ONLY: set_parameters_magnetogram
    use ModPhysics,     ONLY: SW_T_dim, SW_n_dim
    use ModCoronalHeating, ONLY: UseUnsignedFluxModel, DecayLength, UseCranmerHeating,&
                                 read_corona_heating, UseExponentialHeating, &
                                 read_active_region_heating,read_longscale_heating,&
                                 DoOpenClosedHeat, WsaT0

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

       case("#CORONALHEATING")
          call read_corona_heating
          if(UseExponentialHeating)call read_active_region_heating
          
       case("#OPENCLOSEDHEAT")
          call read_var('DoOpenClosedHeat', DoOpenClosedHeat)
          if(DoOpenClosedHeat) call read_var('WsaT0', WsaT0)

       case("#LONGSCALEHEAT")
         call read_longscale_heating

       case("#MODIFYHEATCONDUCTION")
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
         UnitRho_, No2Si_V, UnitN_, UnitU_, UnitX_
    use ModProcMH,      ONLY: iProc
    use ModReadParam,   ONLY: i_line_command
    use ModWaves
    use ModMain,        ONLY: optimize_message_pass
    use ModLookupTable, ONLY: i_lookup_table
    use ModConst,       ONLY: cBoltzmann, cElectronMass, cProtonMass, &
         cEps, cElectronCharge, cTwoPi

    real, parameter:: CoulombLog = 20.0

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
    
    iTableRadCool = i_lookup_table('radcool')

    if(i_line_command("#PLASMA", iSessionIn = 1) < 0)then
       write(*,*) 'Need to set electron temp ration with #PLASMA!'
       call stop_mpi('ERROR: Correct PARAM.in!')
    end if

    ! ratio of Te / (Te + Tp)
    ! given by the formula P = ne * k_B*Te + n_i k_B*T_i
    TeFraction = MassIon_I(1)*ElectronTemperatureRatio &
         /(1 + AverageIonCharge*ElectronTemperatureRatio)

    ! *** HEAT CONDUCTION COEFF calc from ModHeatConduction
    ! electron heat conduct coefficient for single charged ions
    ! = 9.2e-12 W/(m*K^(7/2))
    HeatCondParSi = 3.2*3.0*cTwoPi/CoulombLog &
         *sqrt(cTwoPi*cBoltzmann/cElectronMass)*cBoltzmann &
         *((cEps/cElectronCharge)*(cBoltzmann/cElectronCharge))**2

    ! unit HeatCondParSi is W/(m*K^(7/2))
    HeatCondPar = HeatCondParSi &
         *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitTemperature_)**3.5 &
         *Si2No_V(UnitU_)*Si2No_V(UnitX_)

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
       call write_prefix; write(iUnitOut,*) ''
       
       if(DoModHeatConduction)then
          call write_prefix;  write(iUnitOut,*) 'using Modified Heat Conduction'
          call write_prefix;  write(iUnitOut,*) 'TeModSi      = ', TeModSi
          call write_prefix;  write(iUnitOut,*) 'DeltaTeModSi = ', DeltaTeModSi
          call write_prefix;  write(iUnitOut,*) ''
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
    use ModFaceBc,      ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D, &
         iSide, iFace, jFace, kFace
    use ModMain,        ONLY: x_, y_, z_, UseRotatingFrame, GlobalBLK, &
         nI, nJ, nK, East_, West_, South_, North_, Bot_, Top_
    use ModNumConst,    ONLY: cTolerance
    use ModPhysics,     ONLY: OmegaBody, BodyRho_I, BodyTDim_I, &
         UnitTemperature_, Si2No_V, No2Si_V, UnitN_, UnitEnergyDens_, &
         UnitT_, UnitX_
    use ModVarIndexes,  ONLY: nVar, Rho_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, p_
    use ModFaceGradient, ONLY: get_face_gradient
        

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
    if(DoREBModel) Density = calc_reb_density()

    Temperature = BoundaryTe/TeFraction


!    VarsGhostFace_V(Rho_) =  2.0*Density - VarsTrueFace_V(Rho_)
!    VarsGhostFace_V(p_) = VarsGhostFace_V(Rho_)*Temperature
    VarsGhostFace_V(Rho_) =  Density
    VarsGhostFace_V(p_) = Density*Temperature

    !\
    ! Apply corotation if needed
    !/
    if(.not.UseRotatingFrame)then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2.0*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2.0*OmegaBody*FaceCoords_D(x_)
    end if

  contains
    !==========================================================================
    real function calc_reb_density()

      ! function to return the density given by the Radiative Energy Balance Model
      ! (REB) for the Transition region. Originally given in Withbroe 1988, this
      ! uses eq from Lionell 2001. NO enthalpy flux correction in this
      ! implementation.

      real :: FaceGrad_D(3),GradTeSi_D(3)

      ! Here Rad integral is integral of lossfunction*T^(1/2) from T=10,000 to
      ! 500,000. Use same approximate loss function used in BATS to calculate
      ! This is in SI units [J m^3 K^(3/2)]
      real :: RadIntegralSi = 1.009E-26

      ! Left and right cell centered heating
      real :: CoronalHeatingLeft, CoronalHeatingRight, CoronalHeating

      ! Condensed terms in the REB equation
      real :: qCondSi, qHeatSi

      integer :: iBlock, iDir=0

      !--------------------------------------------------------------------------

      iBlock = GlobalBLK

      ! need to get direction for face gradient calc
      ! also put left cell centered heating call here (since index depends on
      ! the direction)
      if(iSide==East_ .or. iSide==West_) then 
         iDir = x_
         call get_cell_heating(iFace-1, jFace, kFace, iBlock, CoronalHeatingLeft)
      elseif(iSide==South_ .or. iSide==North_) then 
         iDir = y_
         call get_cell_heating(iFace, jFace-1, kFace, iBlock, CoronalHeatingLeft)
      elseif(iSide==Bot_ .or. iSide==Top_) then
         iDir = z_
         call get_cell_heating(iFace, jFace, kFace-1, iBlock, CoronalHeatingLeft)
      else
         call stop_mpi('REB model got bad face direction')
      endif

      call get_cell_heating(iFace, jFace, kFace, iBlock, CoronalHeatingRight)

      CoronalHeating = 0.5 * (CoronalHeatingLeft + CoronalHeatingRight)

      ! term based on coronal heating into trans region (calc face centered avg)
      qHeatSi = (2.0/7.0) * CoronalHeating * BoundaryTeSi**1.5 &
                 * No2Si_V(UnitEnergyDens_) / No2Si_V(UnitT_)

      ! now calculate the contribution due to heat conduction into the boundary
      if(IsNewBlockTeCalc(iBlock)) Te_G = State_VGB(P_,:,:,:,iBlock) / &
           State_VGB(Rho_,:,:,:,iBlock) * TeFraction

      call get_face_gradient(iDir, iFace, jFace, kFace, iBlock, &
           IsNewBlockTeCalc(iBlock), Te_G, FaceGrad_D)

      qCondSi = 0.5 * HeatCondParSi * BoundaryTeSi**3 * sum(FaceGrad_D**2) &
           * (No2Si_V(UnitTemperature_) / No2Si_V(UnitX_))**2

      ! put the terms together and calculate the REB density
      calc_reb_density = sqrt((qCondSi + qHeatSi) / RadIntegralSi) &
           * Si2No_V(UnitN_)

    end function calc_reb_density


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
    ! this is from 2 lines that piecewise fit an old relaxed solution
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
    rBase = rBody
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
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0

    end do; end do; end do

  end subroutine user_set_ics

  !============================================================================

  subroutine user_initial_perturbation
    use ModMain, ONLY: nI ,nJ , nK, nBLK, unusedBLK, x_, y_, z_
    use ModProcMH,    ONLY: iProc, iComm
    use ModPhysics,   ONLY: No2Si_V, UnitX_, UnitEnergyDens_, UnitT_, rBody
    use ModGeometry,   ONLY: vInv_CB
    use ModCoronalHeating, ONLY: TotalCoronalHeatingCgs, &
         UseUnsignedFluxModel, get_coronal_heat_factor
    use ModProcMH,      ONLY: iProc
    use ModIO,          ONLY: write_prefix, iUnitOut
    use ModMpi
    use ModCoronalHeating,ONLY:UseExponentialHeating,&
         DecayLengthExp,HeatingAmplitudeCGS
    implicit none

    integer :: i, j, k, iBlock, iError
    logical :: oktest, oktest_me

    real :: TotalHeatingProc, TotalHeating, TotalHeatingCgs, CoronalHeating
    real :: TotalHeatingModel = 0.0
    !--------------------------------------------------------------------------
    call set_oktest('user_initial_perturbation',oktest,oktest_me)

    ! Calculate the total power into the Computational Domain, loop over
    ! every cell, add the heating, and after block loop, MPI_reduce

    ! Do this because want to be able to generalize models, which can depend on 
    ! topology of domain --> total heating not always known beforehand 

    ! if using open closed heating initialize auxilary WSA grid
    if(DoOpenClosedHeat) call set_empirical_model(trim('WSA'),WsaT0)
    ! Need to initialize unsigned flux model first
    if(UseUnsignedFluxModel) call get_coronal_heat_factor

    TotalHeatingProc = 0.0

    
    do iBlock=1,nBLK
       if(unusedBLK(iBlock))CYCLE
       do k=1,nK; do j=1,nJ; do i=1,nI

          ! Calc heating (Energy/Volume/Time) for the cell 
          call get_cell_heating(i, j, k, iBlock, CoronalHeating)
          
          ! Multiply by cell volume and add to sum
          TotalHeatingProc = TotalHeatingProc + CoronalHeating &
                          / vInv_CB(i, j, k, iBlock)

       end do; end do; end do

    end do

    ! now collect sum over procs
    call MPI_allreduce(TotalHeatingProc, TotalHeating, 1, &
                       MPI_REAL, MPI_SUM, iComm, iError)

    ! Convert total into CGS units
    TotalHeatingCgs = TotalHeating * No2Si_V(UnitEnergyDens_) * 10.0 &
                      / No2Si_V(UnitT_) * (No2Si_V(UnitX_) * 100.0)**3

    ! now compute the total heating of the main models alone
    ! if it were applied to entire domain (to check consistency)
    if(UseUnsignedFluxModel) TotalHeatingModel = TotalCoronalHeatingCgs

    if(UseExponentialHeating) then
       TotalHeatingModel = HeatingAmplitudeCgs * 4.0 * 3.1415927 &
            * (DecayLengthExp*rBody**2 + 2.0*DecayLengthExp**2 * rBody &
            + 2.0*DecayLengthExp**3) * (No2Si_V(UnitX_)*100.0)**3
    end if

    if(iProc==0) then
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '----- START Coronal Heating #s -----------'
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'Total Heat of uniform single model'&
             //' (ergs / s) = ', TotalHeatingModel
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) 'Total Heat into corona (ergs / s) = ',&
             TotalHeatingCgs
       call write_prefix; write(iUnitOut,*) ''
       call write_prefix; write(iUnitOut,*) '------- END Coronal Heating #s -----------'
       call write_prefix; write(iUnitOut,*) ''
       write(*,*) ''
    end if

  end subroutine user_initial_perturbation

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
         VdtFace_x, VdtFace_y, VdtFace_z
    use ModGeometry,       ONLY: r_BLK, vInv_CB
    use ModMain,           ONLY: nI, nJ, nK, GlobalBlk
    use ModPhysics,        ONLY: Si2No_V, UnitEnergyDens_, UnitTemperature_, &
         inv_gm1

    integer :: i, j, k, iBlock
    real :: CoronalHeating, RadiativeCooling, EinternalSource

    ! variables for checking timestep control
    logical, parameter :: DoCalcTime = .true.
    real :: TimeInvRad, TimeInvHeat, Einternal, Vdt_MaxSource

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    !--------------------------------------------------------------------------

    iBlock = globalBlk

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       call get_radiative_cooling(i, j, k, iBlock, RadiativeCooling)

       EinternalSource = RadiativeCooling

       Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) + EinternalSource

       ! Add this in for tentative timestep control from large source terms
       ! need this because radiative loss term becomes INSANELY large at
       ! Chromospheric densities
       if(DoCalcTime)then
          Einternal = inv_gm1 * State_VGB(P_,i,j,k,iBlock)
          TimeInvRad  = abs(RadiativeCooling / Einternal)
          TimeInvHeat = abs(CoronalHeating_C(i,j,k)   / Einternal)
   
          Vdt_MaxSource = (TimeInvRad + TimeInvHeat) / vInv_CB(i,j,k,iBlock)
   
          !**** NOTE This Is a CELL CENTERED TIMESCALE since sources are cell
          ! centered. For now, add to lefthand VdtFace, knowing that calc timestep 
          ! looks at MAX of VdtFaces on all sides
          ! (however cells at the edge of the block will only see one neighbor...) 
          VdtFace_x(i,j,k) = VdtFace_x(i,j,k) + 2.0 * Vdt_maxsource
          VdtFace_y(i,j,k) = VdtFace_y(i,j,k) + 2.0 * Vdt_maxsource
          VdtFace_z(i,j,k) = VdtFace_z(i,j,k) + 2.0 * Vdt_maxsource
       end if 

    end do; end do; end do

  end subroutine user_calc_sources

  !============================================================================

  subroutine get_radiative_cooling(i, j, k, iBlock, RadiativeCooling)

    use ModMultiFluid, ONLY: MassIon_I
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitT_, UnitN_, &
         UnitEnergyDens_, UnitTemperature_
    use ModVarIndexes, ONLY: nVar, Rho_, P_
    use ModLookupTable, ONLY: interpolate_lookup_table
    use ModAdvance,    ONLY: State_VGB

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out):: RadiativeCooling

    real :: Te, TeSi, CoolingFunctionCgs, NumberDensCgs
    real :: Log10TeSi, Log10NeCgs, CoolingTableOut(2)
    real, parameter :: RadNorm = 1.0E+22
    real, parameter :: TeModMinSi = 2.0E+4
    real :: TeFactor, TeModMin, FractionSpitzer
    !--------------------------------------------------------------------------

    Te = TeFraction * State_VGB(P_, i, j, k, iBlock) &
         / State_VGB(Rho_, i, j, k, iBlock)
    TeSi =Te*No2Si_V(UnitTemperature_)

    !if(TeSi<=2*TeModMinSi) then
    !   RadiativeCooling = 0.0
    !   RETURN
    !endif

    ! currently proton plasma only
    NumberDensCgs = State_VGB(Rho_, i, j, k, iBlock) * No2Si_V(UnitN_) * 1.0e-6
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

    ! Avoid extrapolation past zero
    CoolingFunctionCgs = max(CoolingFunctionCgs,0.0)

    RadiativeCooling = -NumberDensCgs**2*CoolingFunctionCgs &
         *0.1*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

    ! include multiplicative factors to make up for extention of
    ! perpendicular heating at low temperatures (as per Abbett 2007).
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

    ! REB model calls face gradient calculation, reset block logical
    ! so that the Te block will be re-calculated next pass
    if(DoREBModel) IsNewBlockTeCalc(iBlock) = .true.

  end subroutine user_update_states

  !============================================================================

  subroutine user_get_log_var(VarValue,TypeVar,Radius)

    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK, nBlock, x_, y_, z_
    use ModVarIndexes, ONLY: Bx_, By_, Bz_, p_ 
    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, B0_DGB
    use ModPhysics,    ONLY: inv_gm1, No2Io_V, UnitEnergydens_, UnitX_, &
         UnitT_, No2Si_V
    use ModCoronalHeating, ONLY: HeatFactor,HeatNormalization
    use ModProcMH,     ONLY: nProc

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

    case('psi')
       VarValue = HeatFactor * No2Si_V(UnitEnergyDens_) / No2Si_V(UnitT_) &
            * 10.0 / nProc * HeatNormalization
       
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
    ! uses ghost cells! If face gradient was checking values other than
    ! P/rho, would need to set those as well!

    if(iSide==East_) then
       State_VGB(Rho_,-1:0,:,:,iBlock) = BoundaryRho
       State_VGB(P_  ,-1:0,:,:,iBlock) = BoundaryRho * BoundaryTe/TeFraction
    else
       call stop_mpi('For TR Model ONLY East_ (low R) user boundary can be used')
    endif

    IsFound = .true.
  end subroutine user_set_outerbcs
  
 !===========================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize,    ONLY: nI, nJ, nK
    use ModPhysics, ONLY: No2Si_V, UnitT_,UnitEnergyDens_
    use ModAdvance,  ONLY: State_VGB

    implicit none

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
    real                         :: UnitEnergyDensPerTime, CoronalHeating
    real                         :: RadiativeCooling
    integer                      :: i, j, k
    !-------------------------------------------------------------------    
    !UsePlotVarBody = .true. 
    !PlotVarBody = 0.0 
    IsFound=.true.

    UnitEnergyDensPerTime = 10.0 * No2Si_V(UnitEnergydens_) / No2Si_V(UnitT_)
    !\                                                                              
    ! Define plot variable to be saved::
    !/ 
    !
    select case(NameVar)
       !Allways use lower case !!

    case('qheat')
       do k=-1,nK+2 ; do j=-1,nJ+2 ; do i=-1,nI+2
          call get_cell_heating(i, j, k, iBlock, CoronalHeating)
          PlotVar_G(i,j,k) = CoronalHeating
       end do ; end do ; end do
       PlotVar_G= PlotVar_G * UnitEnergyDensPerTime
       NameTecVar = 'qH'
       NameTecUnit = '[erg/cm^3/s]'
       NameIdlUnit = '[erg/cm^3/s]'

    case('qrad')
       do k=-1,nK+2 ; do j=-1,nJ+2 ; do i=-1,nI+2
          call get_radiative_cooling(i, j, k, iBlock, RadiativeCooling)
          PlotVar_G(i,j,k) = RadiativeCooling
       end do ; end do ; end do
       PlotVar_G= PlotVar_G * UnitEnergyDensPerTime
       NameTecVar = 'qR'
       NameTecUnit = '[erg/cm^3/s]'
       NameIdlUnit = '[erg/cm^3/s]'

    case default
       IsFound= .false.
    end select
  end subroutine user_set_plot_var

  !============================================================================

  subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
       EinternalIn, TeIn, NatomicOut, &
       EinternalOut, TeOut, PressureOut, &
       CvOut, GammaOut, HeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, &
       PlanckOut_W, CgTeOut_W, CgTgOut_W, TgOut_W)

    ! The State_V vector is in normalized units

    use ModAdvance,    ONLY: nWave, UseElectronPressure
    use ModConst,      ONLY: cBoltzmann
    use ModPhysics,    ONLY: gm1, inv_gm1, No2Si_V, Si2No_V, &
         UnitRho_, UnitP_, UnitEnergyDens_, UnitTemperature_, &
         UnitX_, UnitT_, UnitU_, UnitN_, cRadiationNo, g, Clight
    use ModVarIndexes, ONLY: nVar, Rho_, p_, Pe_, ExtraEint_

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
    real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
    real, optional, intent(out) :: &
         OpacityPlanckOut_W(nWave)                         ! [1/m]
    real, optional, intent(out) :: &
         OpacityRosselandOut_W(nWave)                      ! [1/m]

    ! Multi-group specific interface. The variables are respectively:
    !  Group Planckian spectral energy density
    !  Derivative of group Planckian by electron temperature
    !  Group specific heat of the radiation
    !  Group radiation temperature
    real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]
    real, optional, intent(out) :: CgTeOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: CgTgOut_W(nWave)        ! [J/(m^3*K)]
    real, optional, intent(out) :: TgOut_W(nWave)          ! [K]

    real :: Rho, Pressure, Te, Ti
    real :: RhoSi, pSi, TeSi
    real :: HeatCond
    real :: FractionSpitzer

    character(len=*), parameter :: NameSub = 'user_material_properties'
    !--------------------------------------------------------------------------

    Rho = State_V(Rho_)
    RhoSi = Rho*No2Si_V(Rho_)

    if(present(TeIn))then
       TeSi = TeIn
       Te = TeSi*Si2No_V(UnitTemperature_)
    endif

    if(present(TeOut)) then
       Te = TeFraction * State_V(P_) &
            / State_V(Rho_)
       TeSi =Te*No2Si_V(UnitTemperature_)
       TeOut = TeSi
    endif

    if(present(HeatCondOut))then
       if(DoModHeatConduction) then
          ! Artificial modified heat conduction for a smoother transition
          ! region, Linker et al. (2001)
          FractionSpitzer = 0.5*(1.0+tanh((Te-TeMod)/DeltaTeMod))
          HeatCond = HeatCondPar*(FractionSpitzer*Te**2.5 &
               + (1.0 - FractionSpitzer)*TeMod**2.5)
       else
          HeatCond = HeatCondPar
       endif
       
       HeatCondOut = HeatCond &
            *No2Si_V(UnitEnergyDens_)/No2Si_V(UnitTemperature_) &
            *No2Si_V(UnitU_)*No2Si_V(UnitX_)
    end if

  end subroutine user_material_properties

  !============================================================================


end module ModUser

