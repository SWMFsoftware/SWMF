!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!=========================================! Master module!======================
module ModCoronalHeating

  use BATL_lib, ONLY: &
       test_start, test_stop
  use ModMain,       ONLY: nI, nJ, nK
  use ModReadParam,  ONLY: lStringLine
  use ModVarIndexes, ONLY: WaveFirst_, WaveLast_
  use ModMultiFluid, ONLY: IonFirst_, IonLast_
  use ModExpansionFactors, ONLY: get_bernoulli_integral

  implicit none
  SAVE

  PRIVATE  ! except

  ! The Poynting flux to magnetic field ratio (one of the input parameters
  ! in SI unins and diminsionless:
  real, public :: PoyntingFluxPerBSi = 1.0e6, PoyntingFluxPerB
  real, public :: MaxImbalance = 2.0, MaxImbalance2 = 4.0

  logical, public :: UseCoronalHeating = .false.
  character(len=lStringLine) :: NameModel, TypeCoronalHeating

  ! Exponential Model ---------
  ! Variables and parameters for various heating models

  ! quantitative parameters for exponential heating model

  real :: HeatingAmplitude, HeatingAmplitudeCgs = 6.07e-7
  real :: DecayLengthExp = 0.7  ! in Solar Radii units
  logical :: UseExponentialHeating = .false.

  ! parameters for high-B transition (in Gauss)
  ! idea is to grossly approx 'Active Region' heating values
  logical :: UseArComponent = .false.
  real :: ArHeatB0 = 30.0
  real :: DeltaArHeatB0 = 5.0
  real :: ArHeatFactorCgs = 4.03E-05  ! cgs energy density = [ergs cm-3 s-1]

  ! open closed heat is if you want to have different heating
  ! between open and closed field lines. The routines for the WSA
  ! model in the ExpansionFactors module essentially do this, so will
  ! need to call them
  logical, public :: DoOpenClosedHeat = .false.

  ! Abbett's model -------------

  ! Normalization constant for Abbett Model
  real :: HeatNormalization = 1.0

  ! Alfven wave dissipation
  logical,public :: UseAlfvenWaveDissipation = .false.
  real,   public :: LperpTimesSqrtBSi = 7.5e4 ! m T^(1/2)
  real,   public :: LperpTimesSqrtB
  real    :: Crefl = 0.04

  logical,public :: UseTurbulentCascade = .false.
  logical,public :: UseWaveReflection = .true.

  logical,public :: IsNewBlockAlfven = .true.

  ! long scale height heating (Ch = Coronal Hole)
  logical :: DoChHeat = .false.
  real :: HeatChCgs = 5.0e-7
  real :: DecayLengthCh = 0.7

  ! Arrays for the calculated heat function and dissipated wave energy
  real,public :: CoronalHeating_C(1:nI,1:nJ,1:nK)
  real,public :: WaveDissipation_VC(WaveFirst_:WaveLast_,1:nI,1:nJ,1:nK)

  character(len=lStringLine) :: TypeHeatPartitioning

  ! Switch whether to use uniform heat partition
  logical :: UseUniformHeatPartition = .true.
  real :: QionRatio_I(IonFirst_:IonLast_) = 0.6
  real :: QionParRatio_I(IonFirst_:IonLast_) = 0.0
  real, public :: QeRatio = 0.4

  ! Dimensionless parameters for stochastic heating
  logical :: UseStochasticHeating = .false.
  real :: StochasticExponent = 0.34
  real :: StochasticAmplitude = 0.18

  logical :: DoInit = .true.

  ! The power spectrum of magnetic energy fluctuations in the solar wind
  ! follows a Kolmogorov spectrum (\vardelta B)^(-5/3) in the inertial range
  ! of the turbulence. Podesta et al. (2006) found based on in situ
  ! observations  at 1 AU a spectral index of 3/2 for kinetic energy
  ! fluctuations.
  ! If the logical UseKolmogorov is .false., we assume a 3/2 spectral
  ! index (for the entire solar corona and inner heliosphere) in the
  ! stochastic heating mechanism, otherwise we assume it to be 5/3
  logical :: UseKolmogorov = .false.

  public :: get_coronal_heat_factor
  public :: get_coronal_heating
  public :: get_cell_heating
  public :: get_block_heating
  public :: apportion_coronal_heating
  public :: get_wave_reflection
  public :: init_coronal_heating
  public :: read_corona_heating

  ! Bill Abbet's model, if .true.
  logical, public :: UseUnsignedFluxModel = .false.

  ! Normalized value of Heating constant
  real :: HeatFactor = 0.0

  ! Cgs value of total power input from coronal heating
  real :: TotalCoronalHeatingCgs = 1.0e+28

  ! Exponential Scale height to truncate heating function
  real :: DecayLength = 1.0
  real, public :: DtUpdateFlux = -1.0
  real, public :: UnsignedFluxHeight = -99999.0

contains
  !============================================================================

  subroutine get_coronal_heat_factor

    use ModAdvance,     ONLY: State_VGB, Bz_
    use ModGeometry,    ONLY: true_BLK, true_cell, TypeGeometry
    use ModMagnetogram, ONLY: nTheta, nPhi, dSinTheta, dPhi, &
         get_magnetogram_field
    use ModMain,        ONLY: nI, nJ, nK, nBlock, Unused_B, Time_Simulation,z_
    use ModMpi,         ONLY: MPI_REAL, MPI_SUM
    use ModNumConst,    ONLY: cHalfPi
    use ModPhysics,     ONLY: Si2No_V, No2Si_V, UnitX_, UnitT_, &
         UnitEnergyDens_, rBody
    use ModProcMH,      ONLY: nProc, iComm
    use BATL_lib,       ONLY: CellFace_DB, CellVolume_GB

    integer :: i, j, k, iBlock
    integer :: iTheta, iPhi, iError
    real :: UnsignedFluxCgs, dAreaCgs
    real :: HeatFunction, HeatFunctionVolume, HeatFunctionVolumePe
    real :: x, y, z, Theta, Phi, SinTheta, CosTheta, SinPhi, CosPhi
    real :: B0_D(3), BrSi, BrCgs, SumUnsignedBrCgs
    real :: BzCgs_II(1:nI,1:nJ), SumUnsignedBzCgs, UnsignedFluxCgsPe
    real    :: TotalCoronalHeating = -1.0, TimeUpdateLast = -1.0
    logical :: DoFirst = .true.

    real, parameter :: HeatExponent = 1.1488, HeatCoef = 89.4

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_coronal_heat_factor'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoFirst .and. DtUpdateFlux <= 0.0)then

       ! uniform cell area on sphere
       dAreaCgs = rBody**2*dSinTheta*dPhi*No2Si_V(UnitX_)**2*1e4
       SumUnsignedBrCgs = 0.0

       do iTheta = 0, nTheta

          Theta = cHalfPi - asin((real(iTheta) + 0.5)*dSinTheta - 1.0)
          SinTheta = sin(Theta)
          CosTheta = cos(Theta)
          do iPhi = 1, nPhi
             Phi=(real(iPhi)-0.5)*dPhi
             SinPhi = sin(Phi)
             CosPhi = cos(Phi)

             x = rBody*SinTheta*CosPhi
             y = rBody*SinTheta*SinPhi
             z = rBody*CosTheta

             call get_magnetogram_field(x, y, z, B0_D)
             BrSi = (x*B0_D(1) + y*B0_D(2) + z*B0_D(3))/rBody
             BrCgs = BrSi*1e4
             SumUnsignedBrCgs = SumUnsignedBrCgs + abs(BrCgs)

          end do
       end do

       UnsignedFluxCgs = SumUnsignedBrCgs*dAreaCgs

       TotalCoronalHeatingCgs = HeatCoef*UnsignedFluxCgs**HeatExponent

       TotalCoronalHeating = TotalCoronalHeatingCgs*1e-7 &
            *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

       DoFirst = .false.

    elseif( DtUpdateFlux > 0.0 .and. &
         Time_Simulation - TimeUpdateLast > DtUpdateFlux ) then

       UnsignedFluxCgs = 0.0
       if(TypeGeometry == 'spherical')then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             if(true_BLK(iBlock)) then
                call get_photosphere_unsignedflux(iBlock, UnsignedFluxCgs)
             end if
          end do
       elseif(TypeGeometry == 'cartesian')then
          do iBlock = 1, nBlock
             if(Unused_B(iBlock)) CYCLE
             if(true_BLK(iBlock)) then
                dAreaCgs = CellFace_DB(z_,iBlock)*No2Si_V(UnitX_)**2*1e4

                call get_photosphere_field(iBlock, &
                     State_VGB(Bz_,1:nI,1:nJ,0:nK+1,iBlock), BzCgs_II)

                SumUnsignedBzCgs = sum(abs(BzCgs_II))
                UnsignedFluxCgs = UnsignedFluxCgs +  SumUnsignedBzCgs*dAreaCgs
             end if
          end do
       else
          call CON_stop(NameSub//': '//TypeGeometry// &
              ' geometry is not yet implemented')
       end if
       if(nProc>1)then
          UnsignedFluxCgsPe = UnsignedFluxCgs
          call MPI_allreduce(UnsignedFluxCgsPe, UnsignedFluxCgs, 1, &
               MPI_REAL, MPI_SUM, iComm, iError)
       end if
       TotalCoronalHeatingCgs = HeatCoef*UnsignedFluxCgs**HeatExponent

       TotalCoronalHeating = TotalCoronalHeatingCgs*1e-7 &
            *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)
       TimeUpdateLast = Time_Simulation

    end if

    HeatFunctionVolume = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(true_BLK(iBlock)) then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call get_heat_function(i, j, k, iBlock, HeatFunction)
             HeatFunctionVolume = HeatFunctionVolume &
                  + HeatFunction*CellVolume_GB(i,j,k,iBlock)
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(true_cell(i,j,k,iBlock))then
                call get_heat_function(i, j, k, iBlock, HeatFunction)
                HeatFunctionVolume = HeatFunctionVolume &
                     + HeatFunction*CellVolume_GB(i,j,k,iBlock)
             end if
          end do; end do; end do
       end if
    end do

    if(nProc>1)then
       HeatFunctionVolumePe = HeatFunctionVolume
       call MPI_allreduce(HeatFunctionVolumePe, HeatFunctionVolume, 1, &
            MPI_REAL, MPI_SUM, iComm, iError)
    end if

    HeatFactor = TotalCoronalHeating/HeatFunctionVolume

    call test_stop(NameSub, DoTest)
  end subroutine get_coronal_heat_factor
  !============================================================================

  subroutine get_coronal_heating(i, j, k, iBlock, CoronalHeating)

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating

    real :: HeatFunction
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_coronal_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    call get_heat_function(i, j, k, iBlock, HeatFunction)

    CoronalHeating = HeatFactor*HeatFunction

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_coronal_heating
  !============================================================================

  subroutine get_heat_function(i, j, k, iBlock, HeatFunction)

    use ModMain, ONLY: UseB0, z_
    use ModAdvance, ONLY: State_VGB, Bx_, Bz_
    use ModB0, ONLY: B0_DGB
    use ModGeometry, ONLY: r_BLK
    use BATL_lib, ONLY: Xyz_DGB, IsCartesian

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: HeatFunction

    real :: Bmagnitude, B_D(3)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_heat_function'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(UseB0) then
       B_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if

    Bmagnitude = sqrt(sum(B_D**2))

    if(DtUpdateFlux <= 0.0)then
       HeatFunction = Bmagnitude*exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLength)
    else
       if(IsCartesian)then
          if(Xyz_DGB(z_,i,j,k,iBlock)<UnsignedFluxHeight)then
             HeatFunction = 0.0
          else
             HeatFunction = Bmagnitude
          end if
       else
          if(r_BLK(i,j,k,iBlock)<UnsignedFluxHeight)then
             HeatFunction = 0.0
          else
             HeatFunction = Bmagnitude &
                  *exp(-(r_BLK(i,j,k,iBlock)-1.0)/DecayLength)
          end if
       end if
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_heat_function
  !============================================================================

  subroutine get_photosphere_field(iBlock, Bz_V, BzCgs_II)

    use ModMain,      ONLY: nI, nJ, nK, z_
    use ModInterpolate, ONLY: find_cell
    use ModPhysics,   ONLY: No2Si_V, UnitB_
    use BATL_lib,     ONLY: CoordMin_DB, CoordMax_DB, CellSize_DB

    integer, intent(in) :: iBlock
    real, intent(in)    :: Bz_V(1:nI, 1:nJ, 0:nK+1)
    real, intent(out)   :: BzCgs_II(1:nI, 1:nJ)
    real :: MinZ, MaxZ, DxLeft, z
    integer :: iLeft
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_photosphere_field'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    MinZ = CoordMin_DB(z_,iBlock)
    MaxZ = CoordMax_DB(z_,iBlock)

    BzCgs_II = 0.0
    if((UnsignedFluxHeight > MaxZ) .or. (UnsignedFluxHeight < MinZ)) RETURN

    z = (UnsignedFluxHeight - MinZ)/CellSize_DB(z_,iBlock) + 0.5
    call find_cell(0, nK+1, z, iLeft, DxLeft)

    BzCgs_II = ((1.0 - DxLeft)*Bz_V(1:nI, 1:nJ, iLeft) + &
         DxLeft*Bz_V(1:nI, 1:nJ, iLeft+1))*No2Si_V(UnitB_)*1e4

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_photosphere_field
  !============================================================================
  subroutine get_photosphere_unsignedflux(iBlock, UnsignedFluxCgs)

    use ModAdvance,     ONLY: State_VGB
    use ModGeometry,    ONLY: r_BLK
    use ModMain,        ONLY: nJ, nK, r_
    use ModInterpolate, ONLY: find_cell
    use ModPhysics,     ONLY: No2Si_V, UnitB_, UnitX_
    use ModVarIndexes,  ONLY: Bx_, Bz_
    use BATL_lib,       ONLY: CoordMin_DB, CoordMax_DB, CellSize_DB, &
         CellFace_DFB, Xyz_DGB

    integer, intent(in) :: iBlock
    real, intent(inout) :: UnsignedFluxCgs

    real :: MinR, MaxR, r, DrLeft, BrLeft, BrRight, BrCgs, DrL, dAreaCgs
    integer :: iLeft, j, k, iL
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_photosphere_unsignedflux'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    MinR = CoordMin_DB(r_,iBlock)
    MaxR = CoordMax_DB(r_,iBlock)

    if((UnsignedFluxHeight > MaxR) .or. (UnsignedFluxHeight < MinR)) RETURN

    ! Cells used to interpolate Br
    r = (UnsignedFluxHeight - MinR)/CellSize_DB(r_,iBlock) + 0.5
    call find_cell(0, nI+1, r, iLeft, DrLeft)

    ! Cells used to interpolate face area
    r = r + 0.5
    call find_cell(0, nI+1, r, iL, DrL)
    ! Make sure that CellFace_DFB is never out of index range
    if(iL == 0)then
       iL = 1; DrL = 0.0
    else if(iL == nI+1)then
       iL = nI; DrL = 1.0
    end if

    do k = 1, nK; do j = 1, nJ
       BrLeft = sum(Xyz_DGB(:,iLeft,j,k,iBlock) &
            *State_VGB(Bx_:Bz_,iLeft,j,k,iBlock))/r_BLK(iLeft,j,k,iBlock)
       BrRight = sum(Xyz_DGB(:,iLeft+1,j,k,iBlock) &
            *State_VGB(Bx_:Bz_,iLeft+1,j,k,iBlock))/r_BLK(iLeft+1,j,k,iBlock)

       BrCgs = ((1.0 - DrLeft)*BrLeft + DrLeft*BrRight)*No2Si_V(UnitB_)*1e4

       dAreaCgs = ((1.0-DrL)*CellFace_DFB(r_,iL,j,k,iBlock) &
            +            DrL*CellFace_DFB(r_,iL+1,j,k,iBlock)) &
            *No2Si_V(UnitX_)**2*1e4

       UnsignedFluxCgs = UnsignedFluxCgs + abs(BrCgs)*dAreaCgs
    end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_photosphere_unsignedflux
  !============================================================================
  subroutine read_corona_heating(NameCommand)

    use ModAdvance,    ONLY: UseAnisoPressure
    use ModReadParam,  ONLY: read_var

    integer :: iFluid

    character(len=*), intent(in):: NameCommand
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_corona_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case("#CORONALHEATING")
       call read_var('TypeCoronalHeating', TypeCoronalHeating)

       ! Initialize logicals
       UseCoronalHeating = .true.
       UseUnsignedFluxModel = .false.
       UseExponentialHeating= .false.
       UseAlfvenWaveDissipation = .false.
       UseTurbulentCascade = .false.
       select case(TypeCoronalHeating)
       case('F','none')
          UseCoronalHeating = .false.
       case('exponential')
          UseExponentialHeating = .true.
          call read_var('DecayLengthExp', DecayLengthExp)
          call read_var('HeatingAmplitudeCgs', HeatingAmplitudeCgs)

       case('unsignedflux','Abbett')
          UseUnsignedFluxModel = .true.
          call read_var('DecayLength', DecayLength)
          call read_var('HeatNormalization', HeatNormalization)
       case('alfvenwavedissipation')
          UseAlfvenWaveDissipation = .true.
          call read_var('LperpTimesSqrtBSi', LperpTimesSqrtBSi)
          call read_var('Crefl', Crefl)
       case('turbulentcascade')
          UseAlfvenWaveDissipation = .true.
          UseTurbulentCascade = .true.
          call read_var('UseWaveReflection', UseWaveReflection)
          call read_var('LperpTimesSqrtBSi', LperpTimesSqrtBSi)
       case default
          call stop_mpi('Read_corona_heating: unknown TypeCoronalHeating = ' &
               // TypeCoronalHeating)
       end select
    case('#LIMITIMBALANCE')
       call read_var('MaxImbalance',MaxImbalance)
       MaxImbalance2 = MaxImbalance**2

    case("#POYNTINGFLUX")
       call read_var('PoyntingFluxPerBSi', PoyntingFluxPerBSi)

    case("#ACTIVEREGIONHEATING")
       call read_var('UseArComponent', UseArComponent)
       if(UseArComponent) then
          call read_var('ArHeatFactorCgs', ArHeatFactorCgs)
          call read_var('ArHeatB0', ArHeatB0)
          call read_var('DeltaArHeatB0', DeltaArHeatB0)
       endif

    case("#LONGSCALEHEATING")
       call read_var('DoChHeat', DoChHeat)
       if(DoChHeat)then
          call read_var('HeatChCgs', HeatChCgs)
          call read_var('DecayLengthCh', DecayLengthCh)
       end if

    case("#HEATPARTITIONING")
       UseUniformHeatPartition = .false.
       UseStochasticHeating = .false.
       call read_var('TypeHeatPartitioning', TypeHeatPartitioning)
       select case(TypeHeatPartitioning)
       case('uniform')
          UseUniformHeatPartition = .true.
          do iFluid = IonFirst_, IonLast_
             call read_var('QionRatio', QionRatio_I(iFluid))
          end do
          if(UseAnisoPressure)then
             do iFluid = IonFirst_, IonLast_
                call read_var('QionParRatio', QionParRatio_I(iFluid))
             end do
          end if
          QeRatio = 1.0 - sum(QionRatio_I)
       case('stochasticheating')
          UseStochasticHeating = .true.
          call read_var('StochasticExponent', StochasticExponent)
          call read_var('StochasticAmplitude', StochasticAmplitude)
       case default
          call stop_mpi('Read_corona_heating: unknown TypeHeatPartitioning = '&
               // TypeHeatPartitioning)
       end select

    case("#KOLMOGOROV")
       call read_var('UseKolmogorov', UseKolmogorov)

    case default
       call stop_mpi('Read_corona_heating: unknown command = ' &
            // NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_corona_heating
  !============================================================================
  subroutine init_coronal_heating
    use ModPhysics, ONLY: Si2No_V, UnitEnergyDens_, UnitT_, UnitB_, UnitX_, &
         UnitU_

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_coronal_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    if(.not.DoInit) RETURN
    DoInit = .false.

    if(UseExponentialHeating)then
       HeatingAmplitude =  HeatingAmplitudeCgs*0.1 &
            *Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
    end if

    if(UseAlfvenWaveDissipation)then
       LperpTimesSqrtB = LperpTimesSqrtBSi &
            *Si2No_V(UnitX_)*sqrt(Si2No_V(UnitB_))
    end if

    PoyntingFluxPerB = PoyntingFluxPerBSi &
         *Si2No_V(UnitEnergyDens_)*Si2No_V(UnitU_)/Si2No_V(UnitB_)

    call test_stop(NameSub, DoTest)
  end subroutine init_coronal_heating
  !============================================================================
  subroutine get_cell_heating(i, j, k, iBlock, CoronalHeating)

    use ModGeometry,       ONLY: r_BLK
    use ModPhysics,        ONLY: Si2No_V, UnitEnergyDens_, UnitT_, &
         No2Io_V, UnitB_
    use ModExpansionFactors, ONLY: UMin
    use ModMain,       ONLY: x_, z_, UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModAdvance,    ONLY: State_VGB
    use ModB0,         ONLY: B0_DGB
    use BATL_lib,      ONLY: Xyz_DGB

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CoronalHeating

    real :: HeatCh

    ! parameters for open/closed. This uses WSA model to determine if
    ! cell is in an 'open' or 'closed' field region.
    !
    ! ** NOTE ** WSA does field line tracing on an auxiliary grid.
    ! should really be using the computational domain, but global
    ! feild line tracing for this purpose is not easily implemented
    real :: Ufinal
    real :: UminIfOpen = 290.0
    real :: UmaxIfOpen = 300.0
    real :: Weight
    real :: B_D(3)

    ! local variables for ArHeating (Active Region Heating)
    real :: FractionB, Bcell

    real :: WaveDissipation_V(WaveFirst_:WaveLast_)
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_cell_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseB0)then
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(x_:z_,i,j,k,iBlock)
    else
       B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if

    if(UseAlfvenWaveDissipation)then

       if(UseTurbulentCascade)then
          call turbulent_cascade(i, j, k, iBlock, WaveDissipation_V, &
               CoronalHeating)
       else
          call calc_alfven_wave_dissipation(i, j, k, iBlock, &
               WaveDissipation_V, CoronalHeating)
       end if

    elseif(UseUnsignedFluxModel)then

       call get_coronal_heating(i, j, k, iBlock, CoronalHeating)
       CoronalHeating = CoronalHeating * HeatNormalization

    elseif(UseExponentialHeating)then

       CoronalHeating = HeatingAmplitude &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthExp)
    else
       CoronalHeating = 0.0
    end if

    if(DoOpenClosedHeat)then
       ! If field is less than 1.05 times the minimum speed, mark as closed
       ! Interpolate between 1.05 and 1.10 for smoothness
       UminIfOpen = UMin*1.05
       UmaxIfOpen = UMin*1.1
       call get_bernoulli_integral( &
            Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
       if(UFinal <= UminIfOpen) then
          Weight = 0.0
       else if (UFinal >= UmaxIfOpen) then
          Weight = 1.0
       else
          Weight = (UFinal - UminIfOpen)/(UmaxIfOpen - UminIfOpen)
       end if

       CoronalHeating = (1.0 - Weight) * CoronalHeating
    end if

    if(DoChHeat) then
       HeatCh = HeatChCgs * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
       CoronalHeating = CoronalHeating + HeatCh &
            *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthCh)
    end if

    if(UseExponentialHeating.and.UseArComponent) then

       Bcell = No2Io_V(UnitB_) * sqrt( sum( B_D**2 ) )

       FractionB = 0.5*(1.0+tanh((Bcell - ArHeatB0)/DeltaArHeatB0))
       CoronalHeating = max(CoronalHeating, &
            FractionB * ArHeatFactorCgs * Bcell &
            * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_))

    endif

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_cell_heating
  !============================================================================

  subroutine get_block_heating(iBlock)

    use ModGeometry,       ONLY: r_BLK
    use ModPhysics,        ONLY: Si2No_V, UnitEnergyDens_, UnitT_, &
         No2Io_V, UnitB_
    use ModExpansionFactors, ONLY: UMin
    use ModMain,       ONLY: x_, z_, UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModAdvance,    ONLY: State_VGB
    use ModB0,         ONLY: B0_DGB
    use BATL_lib,      ONLY: Xyz_DGB

    integer, intent(in) :: iBlock

    integer             :: i, j, k
    real :: HeatCh

    ! parameters for open/closed. This uses WSA model to determine if
    ! cell is in an 'open' or 'closed' field region.
    !
    ! ** NOTE ** WSA does field line tracing on an auxiliary grid.
    ! should really be using the computational domain, but global
    ! feild line tracing for this purpose is not easily implemented
    real :: Ufinal
    real :: UminIfOpen = 290.0
    real :: UmaxIfOpen = 300.0
    real :: Weight
    real :: B_D(3)

    ! local variables for ArHeating (Active Region Heating)
    real :: FractionB, Bcell

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_block_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseAlfvenWaveDissipation)then

       if(UseTurbulentCascade)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call turbulent_cascade(i, j, k, iBlock, &
                  WaveDissipation_VC(:,i,j,k), CoronalHeating_C(i,j,k))
          end do; end do; end do
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             call calc_alfven_wave_dissipation(i, j, k, iBlock, &
                  WaveDissipation_VC(:,i,j,k), CoronalHeating_C(i,j,k))
          end do; end do; end do
       end if

    elseif(UseUnsignedFluxModel)then

       do k=1,nK;do j=1,nJ; do i=1,nI

          call get_coronal_heating(i, j, k, iBlock, CoronalHeating_C(i,j,k))
          CoronalHeating_C(i,j,k) = CoronalHeating_C(i,j,k) * HeatNormalization

       end do; end do; end do

    elseif(UseExponentialHeating)then
       do k=1,nK;do j=1,nJ; do i=1,nI

          CoronalHeating_C(i,j,k) = HeatingAmplitude &
               *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthExp)

       end do; end do; end do
    else
       CoronalHeating_C = 0.0
    end if

    if(DoOpenClosedHeat)then
       ! If field is less than 1.05 times the minimum speed, mark as closed
       ! Interpolate between 1.05 and 1.10 for smoothness
       UminIfOpen = UMin*1.05
       UmaxIfOpen = UMin*1.1
       do k=1,nK; do j=1,nJ; do i=1,nI
          call get_bernoulli_integral( &
               Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock), UFinal)
          if(UFinal <= UminIfOpen) then
             Weight = 0.0
          else if (UFinal >= UmaxIfOpen) then
             Weight = 1.0
          else
             Weight = (UFinal - UminIfOpen)/(UmaxIfOpen - UminIfOpen)
          end if

          CoronalHeating_C(i,j,k) = (1.0 - Weight) * CoronalHeating_C(i,j,k)
       end do; end do; end do
    end if

    if(DoChHeat) then
       HeatCh = HeatChCgs * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)
       do k=1,nK; do j=1,nJ; do i=1,nI
          CoronalHeating_C(i,j,k) = CoronalHeating_C(i,j,k) + HeatCh &
               *exp(- max(r_BLK(i,j,k,iBlock) - 1.0, 0.0) / DecayLengthCh)
       end do; end do; end do
    end if

    if(UseExponentialHeating.and.UseArComponent) then
       do k=1,nK; do j=1,nJ; do i=1,nI

          if(UseB0)then
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(x_:z_,i,j,k,iBlock)
          else
             B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          end if

          Bcell = No2Io_V(UnitB_) * sqrt( sum( B_D**2 ) )

          FractionB = 0.5*(1.0+tanh((Bcell - ArHeatB0)/DeltaArHeatB0))
          CoronalHeating_C(i,j,k) = max(CoronalHeating_C(i,j,k), &
               FractionB * ArHeatFactorCgs * Bcell &
               * 0.1 * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_))
       end do; end do; end do
    endif

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_block_heating
  !============================================================================

  subroutine calc_alfven_wave_dissipation(i, j, k, iBlock, WaveDissipation_V, &
       CoronalHeating)

    use ModAdvance, ONLY: State_VGB
    use ModB0,      ONLY: B0_DGB
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Rho_, Bx_, Bz_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out)   :: WaveDissipation_V(WaveFirst_:WaveLast_), &
         CoronalHeating

    real :: EwavePlus, EwaveMinus, FullB_D(3), FullB, Coef
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'calc_alfven_wave_dissipation'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseB0)then
       FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    FullB = sqrt(sum(FullB_D**2))

    Coef = 2.0*sqrt(FullB/State_VGB(Rho_,i,j,k,iBlock))/LperpTimesSqrtB

    EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
    EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

    WaveDissipation_V(WaveFirst_) = Coef*EwavePlus &
         *sqrt(max(EwaveMinus,Crefl**2*EwavePlus))

    WaveDissipation_V(WaveLast_) = Coef*EwaveMinus &
         *sqrt(max(EwavePlus,Crefl**2*EwaveMinus))

    CoronalHeating = sum(WaveDissipation_V)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine calc_alfven_wave_dissipation
  !============================================================================

  subroutine turbulent_cascade(i, j, k, iBlock, WaveDissipation_V, &
       CoronalHeating)

    use ModAdvance, ONLY: State_VGB
    use ModB0, ONLY: B0_DGB
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModMultiFluid, ONLY: iRho_I, IonFirst_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out)   :: WaveDissipation_V(WaveFirst_:WaveLast_), &
         CoronalHeating

    real :: FullB_D(3), FullB, Coef
    real :: EwavePlus, EwaveMinus

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'turbulent_cascade'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Low-frequency cascade due to small-scale nonlinearities

    if(UseB0)then
       FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
    else
       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    FullB = sqrt(sum(FullB_D**2))

    Coef = 2.0*sqrt(FullB/State_VGB(iRho_I(IonFirst_),i,j,k,iBlock)) &
         /LperpTimesSqrtB

    EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
    EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

    WaveDissipation_V(WaveFirst_) = Coef*sqrt(EwaveMinus)*EwavePlus
    WaveDissipation_V(WaveLast_) = Coef*sqrt(EwavePlus)*EwaveMinus

    CoronalHeating = sum(WaveDissipation_V)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine turbulent_cascade
  !============================================================================

  subroutine get_wave_reflection(iBlock)

    use BATL_size, ONLY: nDim, nI, nJ, nK
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModB0, ONLY: B0_DGB
    use ModChromosphere,  ONLY: DoExtendTransitionRegion, extension_factor, &
         get_tesi_c, TeSi_C
    use ModGeometry, ONLY: true_cell
    use ModMain, ONLY: UseB0
    use ModVarIndexes, ONLY: Bx_, Bz_
    use ModMultiFluid, ONLY: iRho_I, IonFirst_

    integer, intent(in) :: iBlock

    integer :: i, j, k
    real :: GradLogAlfven_D(nDim), CurlU_D(3), b_D(3)
    real :: FullB_D(3), FullB, Rho, DissipationRateMax, ReflectionRate
    real :: EwavePlus, EwaveMinus

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_wave_reflection'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(DoExtendTransitionRegion) call get_tesi_c(iBlock, TeSi_C)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       call get_grad_log_alfven_speed(i, j, k, iBlock, GradLogAlfven_D)

       call get_curl_u(i, j, k, iBlock, CurlU_D)

       if(UseB0)then
          FullB_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
       else
          FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end if
       FullB = sqrt(sum(FullB_D**2))
       b_D = FullB_D/max(1e-15, FullB)

       Rho = State_VGB(iRho_I(IonFirst_),i,j,k,iBlock)

       EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
       EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)

       DissipationRateMax = &
            2.0*sqrt(max(EwavePlus,EwaveMinus)*FullB/Rho)/LperpTimesSqrtB

       if(DoExtendTransitionRegion) DissipationRateMax = &
            DissipationRateMax/extension_factor(TeSi_C(i,j,k))

       ! Reflection rate driven by Alfven speed gradient and
       ! vorticity along the field lines
       ! if(R_BLK(i,j,k,iBlock)<2.5)then
       ReflectionRate = sqrt( (sum(b_D*CurlU_D))**2 &
            + (sum(FullB_D(:nDim)*GradLogAlfven_D))**2/Rho )
       ! else
       !   ReflectionRate = sqrt( (sum(b_D*CurlU_D))**2 &
       !        + min((sum(FullB_D(1:nDim)*GradLogAlfven_D))**2,&
       !        (sum(FullB_D(1:nDim)*Xyz_DGB(1:nDim,i,j,k,iBlock))/&
       !        R_BLK(i,j,k,iBlock)**2)**2)/Rho )
       ! end if
       ! Clip the reflection rate from above with maximum dissipation rate
       ReflectionRate = min(ReflectionRate, DissipationRateMax)

       ! No reflection when turbulence is balanced (waves are then
       ! assumed to be uncorrelated)
       if(MaxImbalance2*EwaveMinus <= EwavePlus)then
          ReflectionRate = ReflectionRate*&
               (1.0 - MaxImbalance*sqrt(EwaveMinus/EwavePlus))
       elseif(MaxImbalance2*EwavePlus <= EwaveMinus)then
          ReflectionRate = ReflectionRate*&
               (MaxImbalance*sqrt(EwavePlus/EwaveMinus)-1.0)
       else
          ReflectionRate = 0.0
       end if

       Source_VC(WaveFirst_,i,j,k) = Source_VC(WaveFirst_,i,j,k) &
            - ReflectionRate*sqrt(EwavePlus*EwaveMinus)
       Source_VC(WaveLast_,i,j,k) = Source_VC(WaveLast_,i,j,k) &
            + ReflectionRate*sqrt(EwavePlus*EwaveMinus)

    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_wave_reflection
  !============================================================================

  subroutine get_grad_log_alfven_speed(i, j, k, iBlock, GradLogAlfven_D)

    use BATL_lib, ONLY: IsCartesianGrid, &
         CellSize_DB, FaceNormal_DDFB, CellVolume_GB, &
         x_, y_, z_, Dim1_, Dim2_, Dim3_
    use BATL_size, ONLY: nDim, nI, j0_, nJp1_, k0_, nKp1_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: GradLogAlfven_D(nDim)

    real, save :: LogAlfven_FD(0:nI+1,j0_:nJp1_,k0_:nKp1_,nDim)
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_grad_log_alfven_speed'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(IsNewBlockAlfven)then
       call get_log_alfven_speed

       IsNewBlockAlfven = .false.
    end if

    if(IsCartesianGrid)then
       GradLogAlfven_D(Dim1_) = 1.0/CellSize_DB(x_,iBlock) &
            *(LogAlfven_FD(i+1,j,k,Dim1_) - LogAlfven_FD(i,j,k,Dim1_))
       if(nJ > 1) GradLogAlfven_D(Dim2_) = 1.0/CellSize_DB(y_,iBlock) &
            *(LogAlfven_FD(i,j+1,k,Dim2_) - LogAlfven_FD(i,j,k,Dim2_))
       if(nK > 1) GradLogAlfven_D(Dim3_) = 1.0/CellSize_DB(z_,iBlock) &
            *(LogAlfven_FD(i,j,k+1,Dim3_) - LogAlfven_FD(i,j,k,Dim3_))
    else
       GradLogAlfven_D = &
            LogAlfven_FD(i+1,j,k,Dim1_)*FaceNormal_DDFB(:,Dim1_,i+1,j,k,iBlock) &
            - LogAlfven_FD(i,j,k,Dim1_)*FaceNormal_DDFB(:,Dim1_,i,j,k,iBlock)
       if(nJ > 1) GradLogAlfven_D = GradLogAlfven_D + &
            LogAlfven_FD(i,j+1,k,Dim2_)*FaceNormal_DDFB(:,Dim2_,i,j+1,k,iBlock) &
            - LogAlfven_FD(i,j,k,Dim2_)*FaceNormal_DDFB(:,Dim2_,i,j,k,iBlock)
       if(nK > 1) GradLogAlfven_D = GradLogAlfven_D + &
            LogAlfven_FD(i,j,k+1,Dim3_)*FaceNormal_DDFB(:,Dim3_,i,j,k+1,iBlock) &
            - LogAlfven_FD(i,j,k,Dim3_)*FaceNormal_DDFB(:,Dim3_,i,j,k,iBlock)

       GradLogAlfven_D = GradLogAlfven_D/CellVolume_GB(i,j,k,iBlock)
    end if

    call test_stop(NameSub, DoTest, iBlock)
  contains
    !==========================================================================

    subroutine get_log_alfven_speed

      use ModAdvance, ONLY: &
           LeftState_VX,  LeftState_VY,  LeftState_VZ,  &
           RightState_VX, RightState_VY, RightState_VZ
      use ModB0, ONLY: B0_DX, B0_DY, B0_DZ
      use ModMain, ONLY: UseB0
      use ModVarIndexes, ONLY: Bx_, Bz_
      use ModMultiFluid, ONLY: iRho_I, IonFirst_

      integer :: i, j, k
      real :: Rho, FullB_D(3)
      !------------------------------------------------------------------------
      do k = 1, nK; do j = 1, nJ; do i = 1, nI+1
         FullB_D = 0.5*(LeftState_VX(Bx_:Bz_,i,j,k) &
              + RightState_VX(Bx_:Bz_,i,j,k))
         if(UseB0) FullB_D = FullB_D + B0_DX(:,i,j,k)
         Rho = 0.5*(LeftState_VX(iRho_I(IonFirst_),i,j,k) &
              +     RightState_VX(iRho_I(IonFirst_),i,j,k))
         LogAlfven_FD(i,j,k,x_) = log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
      end do; end do; end do

      if(nJ > 1)then
         do k = 1, nK; do j = 1, nJ+1; do i = 1, nI
            FullB_D = 0.5*(LeftState_VY(Bx_:Bz_,i,j,k) &
                 + RightState_VY(Bx_:Bz_,i,j,k))
            if(UseB0) FullB_D = FullB_D + B0_DY(:,i,j,k)
            Rho = 0.5*(LeftState_VY(iRho_I(IonFirst_),i,j,k) &
                 +     RightState_VY(iRho_I(IonFirst_),i,j,k))
            LogAlfven_FD(i,j,k,Dim2_) = &
                 log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
         end do; end do; end do
      end if

      if(nK > 1)then
         do k = 1, nK+1; do j = 1, nJ; do i = 1, nI
            FullB_D = 0.5*(LeftState_VZ(Bx_:Bz_,i,j,k) &
                 + RightState_VZ(Bx_:Bz_,i,j,k))
            if(UseB0) FullB_D = FullB_D + B0_DZ(:,i,j,k)
            Rho = 0.5*(LeftState_VZ(iRho_I(IonFirst_),i,j,k) &
                 +     RightState_VZ(iRho_I(IonFirst_),i,j,k))
            LogAlfven_FD(i,j,k,Dim3_) = &
                 log(sqrt(max(sum(FullB_D**2), 1e-30)/Rho))
         end do; end do; end do
      end if

    end subroutine get_log_alfven_speed
    !==========================================================================

  end subroutine get_grad_log_alfven_speed
  !============================================================================

  subroutine get_curl_u(i, j, k, iBlock, CurlU_D)

    use BATL_lib, ONLY: IsCartesianGrid, CellSize_DB, FaceNormal_DDFB, &
         CellVolume_GB, x_, y_, z_
    use ModAdvance, ONLY: &
         LeftState_VX,  LeftState_VY,  LeftState_VZ,  &
         RightState_VX, RightState_VY, RightState_VZ
    use ModCoordTransform, ONLY: cross_product
    use ModSize, ONLY: MaxDim
    use ModMultiFluid, ONLY: iUx_I, iUy_I, iUz_I, IonFirst_

    integer, intent(in) :: i, j, k, iBlock
    real, intent(out) :: CurlU_D(MaxDim)

    real :: DxInvHalf, DyInvHalf, DzInvHalf

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'get_curl_u'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(IsCartesianGrid)then
       DxInvHalf = 0.5/CellSize_DB(x_,iBlock)
       DyInvHalf = 0.5/CellSize_DB(y_,iBlock)
       DzInvHalf = 0.5/CellSize_DB(z_,iBlock)

       CurlU_D(x_) = &
            DyInvHalf*(LeftState_VY(iUz_I(IonFirst_),i,j+1,k)  &
            +          RightState_VY(iUz_I(IonFirst_),i,j+1,k) &
            -          LeftState_VY(iUz_I(IonFirst_),i,j,k)    &
            -          RightState_VY(iUz_I(IonFirst_),i,j,k))  &
            - DzInvHalf*(LeftState_VZ(iUy_I(IonFirst_),i,j,k+1)  &
            +            RightState_VZ(iUy_I(IonFirst_),i,j,k+1) &
            -            LeftState_VZ(iUy_I(IonFirst_),i,j,k)    &
            -            RightState_VZ(iUy_I(IonFirst_),i,j,k))

       CurlU_D(y_) = &
            DzInvHalf*(LeftState_VZ(iUx_I(IonFirst_),i,j,k+1)  &
            +          RightState_VZ(iUx_I(IonFirst_),i,j,k+1) &
            -          LeftState_VZ(iUx_I(IonFirst_),i,j,k)    &
            -          RightState_VZ(iUx_I(IonFirst_),i,j,k))  &
            - DxInvHalf*(LeftState_VX(iUz_I(IonFirst_),i+1,j,k)  &
            +            RightState_VX(iUz_I(IonFirst_),i+1,j,k) &
            -            LeftState_VX(iUz_I(IonFirst_),i,j,k)    &
            -            RightState_VX(iUz_I(IonFirst_),i,j,k))

       CurlU_D(z_) = &
            DxInvHalf*(LeftState_VX(iUy_I(IonFirst_),i+1,j,k)  &
            +          RightState_VX(iUy_I(IonFirst_),i+1,j,k) &
            -          LeftState_VX(iUy_I(IonFirst_),i,j,k)    &
            -          RightState_VX(iUy_I(IonFirst_),i,j,k))  &
            - DyInvHalf*(LeftState_VY(iUx_I(IonFirst_),i,j+1,k)  &
            +            RightState_VY(iUx_I(IonFirst_),i,j+1,k) &
            -            LeftState_VY(iUx_I(IonFirst_),i,j,k)    &
            -            RightState_VY(iUx_I(IonFirst_),i,j,k))
    else
       CurlU_D(:) = &
            + cross_product( FaceNormal_DDFB(:,1,i+1,j,k,iBlock),       &
            LeftState_VX(iUx_I(IonFirst_):iUz_I(IonFirst_),i+1,j,k)     &
            + RightState_VX(iUx_I(IonFirst_):iUz_I(IonFirst_),i+1,j,k)) &
            - cross_product( FaceNormal_DDFB(:,1,i  ,j,k,iBlock),       &
            LeftState_VX(iUx_I(IonFirst_):iUz_I(IonFirst_),i  ,j,k)     &
            + RightState_VX(iUx_I(IonFirst_):iUz_I(IonFirst_),i  ,j,k)) &
            + cross_product( FaceNormal_DDFB(:,2,i,j+1,k,iBlock),       &
            LeftState_VY(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j+1,k)     &
            + RightState_VY(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j+1,k)) &
            - cross_product( FaceNormal_DDFB(:,2,i,j  ,k,iBlock),       &
            LeftState_VY(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j  ,k)     &
            + RightState_VY(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j  ,k)) &
            + cross_product( FaceNormal_DDFB(:,3,i,j,k+1,iBlock),       &
            LeftState_VZ(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j,k+1)     &
            + RightState_VZ(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j,k+1)) &
            - cross_product( FaceNormal_DDFB(:,3,i,j,k  ,iBlock),       &
            LeftState_VZ(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j,k  )     &
            + RightState_VZ(iUx_I(IonFirst_):iUz_I(IonFirst_),i,j,k  ))

       CurlU_D(:) = 0.5*CurlU_D(:)/CellVolume_GB(i,j,k,iBlock)
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine get_curl_u
  !============================================================================

  subroutine apportion_coronal_heating(i, j, k, iBlock, &
       CoronalHeating, QPerQtotal_I, QparPerQtotal_I, QePerQtotal)

    ! Apportion the coronal heating to the electrons and protons based on
    ! how the Alfven waves dissipate at length scales << Lperp

    use ModMain, ONLY: UseB0
    use ModPhysics, ONLY: IonMassPerCharge
    use ModAdvance, ONLY: State_VGB, UseAnisoPressure, &
         Bx_, Bz_, Pe_
    use ModB0, ONLY: B0_DGB
    use ModChromosphere,  ONLY: DoExtendTransitionRegion, extension_factor, &
         TeSi_C
    use ModMultiFluid, ONLY: ChargeIon_I, MassIon_I, UseMultiIon, &
         nIonFluid, IonFirst_, iRhoIon_I, iRhoUxIon_I, iRhoUzIon_I, iPIon_I, &
         iPparIon_I

    integer, intent(in) :: i, j, k, iBlock
    real, intent(in) :: CoronalHeating
    real, intent(out) :: QPerQtotal_I(nIonFluid), &
         QparPerQtotal_I(nIonFluid), QePerQtotal

    integer :: iIon
    real :: Qtotal, Udiff_D(3), Upar, Valfven, Vperp
    real :: Ne, B_D(3), B, B2, InvGyroRadius, AlfvenRatio, DeltaU, Epsilon
    real :: TeByTp, BetaElectron, BetaProton, Pperp, LperpInvGyroRad
    real :: Ewave, EwavePlus, EwaveMinus, EkinCascade
    real :: DampingElectron, DampingPar_I(nIonFluid) = 0.0
    real :: DampingPerp_I(nIonFluid), DampingTotal

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'apportion_coronal_heating'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseStochasticHeating)then
       ! Damping rates and wave energy partition based on Chandran et al.[2011]

       if(DoExtendTransitionRegion)then
          Qtotal = CoronalHeating*extension_factor(TeSi_C(i,j,k))
       else
          Qtotal = CoronalHeating
       end if

       TeByTp = State_VGB(Pe_,i,j,k,iBlock) &
            /max(State_VGB(iPIon_I(1),i,j,k,iBlock), 1e-15)

       if(UseMultiIon)then
          Ne = sum(State_VGB(iRhoIon_I,i,j,k,iBlock)*ChargeIon_I/MassIon_I)
          TeByTp = TeByTp*State_VGB(iRhoIon_I(1),i,j,k,iBlock)/Ne
       end if

       if(UseB0) then
          B_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
       else
          B_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end if
       B2 = max(sum(B_D**2), 1e-30)
       B = sqrt(B2)

       BetaElectron = 2.0*State_VGB(Pe_,i,j,k,iBlock)/B2
       BetaProton = 2.0*State_VGB(iPIon_I(1),i,j,k,iBlock)/B2

       ! Linear Landau damping and transit-time damping of kinetic Alfven
       ! waves contributes to electron and parallel proton heating
       DampingElectron = 0.01*sqrt(TeByTp/max(BetaProton, 1.0e-8)) &
            *(1.0 + 0.17*BetaProton**1.3)/(1.0 +(2800.0*BetaElectron)**(-1.25))
       DampingPar_I(1) = 0.08*sqrt(sqrt(TeByTp))*BetaProton**0.7 &
            *exp(-1.3/max(BetaProton, 1.0e-8))

       EwavePlus  = State_VGB(WaveFirst_,i,j,k,iBlock)
       EwaveMinus = State_VGB(WaveLast_,i,j,k,iBlock)
       Ewave = EwavePlus + EwaveMinus

       do iIon = 1, nIonFluid

          if(UseAnisoPressure)then
             Pperp = 0.5*(3*State_VGB(iPIon_I(iIon),i,j,k,iBlock) &
                  - State_VGB(iPparIon_I(IonFirst_-1+iIon),i,j,k,iBlock))
          else
             Pperp = State_VGB(iPIon_I(iIon),i,j,k,iBlock)
          end if

          ! Perpendicular ion thermal speed
          Vperp = sqrt(2.0*Pperp/State_VGB(iRhoIon_I(iIon),i,j,k,iBlock))

          InvGyroRadius = &
               B/(IonMassPerCharge*MassIon_I(iIon)/ChargeIon_I(iIon))/Vperp
          LperpInvGyroRad = InvGyroRadius*LperpTimesSqrtB/sqrt(B)

          if(UseKolmogorov)then
             EkinCascade = 1.0/LperpInvGyroRad**(2.0/3.0)
          else
             EkinCascade = 1.0/sqrt(LperpInvGyroRad)
          end if

          if(iIon == 1)then
             DeltaU = sqrt(Ewave/State_VGB(iRhoIon_I(1),i,j,k,iBlock) &
                  *EkinCascade)
          else
             ! difference bulk speed between ions and protons
             Udiff_D = &
                  State_VGB(iRhoUxIon_I(iIon):iRhoUzIon_I(iIon),i,j,k,iBlock) &
                  /State_VGB(iRhoIon_I(iIon),i,j,k,iBlock) &
                  -State_VGB(iRhoUxIon_I(1):iRhoUzIon_I(1),i,j,k,iBlock) &
                  /State_VGB(iRhoIon_I(1),i,j,k,iBlock)
             Upar = sum(Udiff_D*B_D)/B
             Valfven = B/sqrt(State_VGB(iRhoIon_I(1),i,j,k,iBlock))

             ! Alfven ratio at ion gyro scale under the assumption that the
             ! kinetic and magnetic fluctuation energies at the Lperp scale are
             ! near equipartition (with a full wave reflection in the
             ! transition region this would not be true)
             if(UseKolmogorov)then
                DeltaU = sqrt( (EwavePlus*(1 - Upar/Valfven)**2 &
                     + EwaveMinus*(1 + Upar/Valfven)**2) &
                     *EkinCascade/State_VGB(iRhoIon_I(1),i,j,k,iBlock) )
             else
                AlfvenRatio = LperpInvGyroRad**(1.0/6.0)

                DeltaU = sqrt(max((AlfvenRatio+(Upar/Valfven)**2)/AlfvenRatio &
                     *Ewave - 2.0*(Upar/Valfven)*(EwavePlus - EwaveMinus) &
                     /sqrt(AlfvenRatio), 0.0) &
                     *EkinCascade/State_VGB(iRhoIon_I(1),i,j,k,iBlock))
             end if
          end if

          Epsilon = DeltaU/Vperp

          ! Nonlinear damping/stochastic heating to perpendicular ion heating
          DampingPerp_I(iIon) = StochasticAmplitude*InvGyroRadius &
               *State_VGB(iRhoIon_I(iIon),i,j,k,iBlock)*DeltaU**3 &
               *exp(-StochasticExponent/max(Epsilon,1e-15)) &
               /max(Qtotal, 1e-30)
       end do

       ! The 1+ is due to the fraction of the cascade power that succeeds
       ! to cascade to the smallest scale (<< proton gyroradius),
       ! where the dissipation is via interactions with the electrons
       DampingTotal = 1.0 + DampingElectron &
            + sum(DampingPar_I) + sum(DampingPerp_I)

       QPerQtotal_I = (DampingPar_I + DampingPerp_I)/Dampingtotal
       QparPerQtotal_I = DampingPar_I/Dampingtotal
       QePerQtotal = (1.0 + DampingElectron)/DampingTotal

    elseif(UseUniformHeatPartition)then
       QPerQtotal_I = QionRatio_I
       QparPerQtotal_I = QionParRatio_I
       QePerQtotal = QeRatio

    else
       call stop_mpi(NameSub//' Unknown heat partitioning')
    end if

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine apportion_coronal_heating
  !============================================================================

end module ModCoronalHeating
!==============================================================================
