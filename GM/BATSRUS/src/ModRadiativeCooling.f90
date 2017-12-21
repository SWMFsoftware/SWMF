!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModRadiativeCooling

  use BATL_lib, ONLY: &
       test_start, test_stop, StringTest
  use ModChromosphere
  use ModSize
  implicit none
  save

  logical:: UseRadCooling      =.false.
  logical:: UseRadCoolingTable =.false.
  integer,private:: iTableRadCool      =-1

  real :: AuxTeSi
  !\
  ! Recommended usage of get_radiative_cooling function
  !/
  ! use ModRadiativeCooling
  ! call user_material_properties(i,j,k,iBlock,TeOut=AuxTeSi)
  ! call get_radiative_cooling(State_VGB(i, j, k, iBlock,AuxTeSi, RadCooling)
  !
  ! The dimensionless cooling function WITH NEGATIVE SIGN is
  ! in RadCooling

  real :: RadCooling_C(1:nI,1:nJ,1:nK)

  real, parameter :: Cgs2SiEnergyDens = &
       1.0e-7&   ! erg = 1e-7 J
       /1.0e-6    ! cm3 = 1e-6 m3
contains
  !============================================================================
  subroutine read_modified_cooling
    use ModReadParam
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_modified_cooling'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    call read_var('DoExtendTransitionRegion',DoExtendTransitionRegion)
    if(DoExtendTransitionRegion)then

       call read_var('TeModSi',TeModSi)
       call read_var('DeltaTeModSi',DeltaTeModSi)
    else
       call read_var('TeTransitionRegionTopSi',TeTransitionRegionTopSi)
    end if
    call test_stop(NameSub, DoTest)
  end subroutine read_modified_cooling
  !============================================================================

  subroutine check_cooling_param
    use ModLookupTable, ONLY: i_lookup_table
    !--------------------------------------------------------------------------

    iTableRadCool = i_lookup_table('radcool')
    UseRadCoolingTable = iTableRadCool>0

  end subroutine check_cooling_param
  !============================================================================
  subroutine get_radiative_cooling(i, j, k, iBlock, TeSiIn, RadiativeCooling, iError)

    use ModPhysics,    ONLY: No2Si_V, UnitN_
    use ModVarIndexes, ONLY: Rho_
    use ModAdvance,    ONLY: State_VGB
    use ModMultiFluid, ONLY: UseMultiIon, ChargeIon_I, iRhoIon_I, MassIon_I

    integer, intent(in) :: i, j, k, iBlock
    real,    intent(in) :: TeSiIn
    real,   intent(out) :: RadiativeCooling
    integer,intent(out),optional::iError
    real :: NumberDensCgs
    character(len=*), parameter:: NameSub = 'get_radiative_cooling'
    !--------------------------------------------------------------------------
    if(present(iError))iError=0
    if(UseMultiIon)then
       NumberDensCgs = &
            sum(ChargeIon_I*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I) &
            *No2Si_V(UnitN_)*1e-6
    else
       NumberDensCgs = State_VGB(Rho_,i,j,k,iBlock)*No2Si_V(UnitN_)*1.0e-6
    end if

    RadiativeCooling = - radiative_cooling(TeSiIn, NumberDensCgs,iError)

    ! include multiplicative factors to make up for extention of
    ! perpendicular heating at low temperatures (as per Abbett 2007).
    ! Need this to strech transition region to larger scales
    ! Also, need radcool modification to become const below TeModMin
    if(DoExtendTransitionRegion) then
       if(TeSiIn >= TeChromosphereSi) then
          RadiativeCooling = RadiativeCooling / extension_factor(TeSiIn)
       else
          RadiativeCooling = RadiativeCooling * (TeChromosphereSi/TeModSi)**2.5
       endif
    end if
  end subroutine get_radiative_cooling
  !============================================================================
  real function radiative_cooling(TeSiIn, NumberDensCgs,iError)
    use ModPhysics,       ONLY: Si2No_V, UnitT_, UnitEnergyDens_
    use ModLookupTable,   ONLY: interpolate_lookup_table
    use ModMultiFluid,    ONLY: UseMultiIon

    ! Imput - dimensional, output: dimensionless
    real, intent(in):: TeSiIn, NumberDensCgs
    integer,optional,intent(out)::iError

    real :: CoolingFunctionCgs
    real :: CoolingTableOut_I(1)
    real, parameter :: RadNorm = 1.0E+22

    character(len=*), parameter:: NameSub = 'radiative_cooling'
    !--------------------------------------------------------------------------

    if(UseRadCoolingTable) then
       ! at the moment, radC not a function of Ne, but need a dummy 2nd
       ! index, and might want to include Ne dependence in table later.
       ! Table variable should be normalized to radloss_cgs * 10E+22
       ! since we don't want to deal with such tiny numbers
       if(index(StringTest, NameSub)>0)then
          if(TeSiIn<1.0e4.or.TeSiIn>1.0e8.or.NumberDensCgs<1.0e2.or.&
               NumberDensCgs<1.0e2.or.NumberDensCgs>1.0e14)then
             write(*,*)'TeSiIn, NumberDensCgs=',TeSiIn, NumberDensCgs
             if(present(iError))iError=1
             if(TeSiIn<0.0)call stop_mpi(&
                  'Negative temperature in calculating rad. cooling')
          end if
       end if
       call interpolate_lookup_table(iTableRadCool, max(TeSiIn,1.0e4),&
            CoolingTableOut_I, DoExtrapolate = .false.)
       CoolingFunctionCgs = CoolingTableOut_I(1) / RadNorm
       ! The presently stored cooling function assumes that it will be
       ! multiplied by Nelectron*Nhydrogen. The following turns this
       ! Nelectron*Nhydrogen into Nelectron**2. There is some slight
       ! temperature dependence in this multiplicative coefficient below
       ! 100,000 K, but we ignore that here.
       if(UseMultiIon) CoolingFunctionCgs = CoolingFunctionCgs*0.83
    else
       call get_cooling_function_fit(TeSiIn, CoolingFunctionCgs)
    end if

    ! Avoid extrapolation past zero
    CoolingFunctionCgs = max(CoolingFunctionCgs, 0.0)

    radiative_cooling = NumberDensCgs**2*CoolingFunctionCgs &
         * Cgs2SiEnergyDens * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

  end function radiative_cooling
  !============================================================================
  subroutine get_cooling_function_fit(TeSi, CoolingFunctionCgs)

    ! Based on Rosner et al. (1978) and Peres et al. (1982)
    ! Need to be replaced by Chianti tables

    real, intent(in) :: TeSi
    real, intent(out):: CoolingFunctionCgs

    character(len=*), parameter:: NameSub = 'get_cooling_function_fit'
    !--------------------------------------------------------------------------
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
  !============================================================================
  real function cooling_function_integral_si(TeTransitionRegionSi)

    real,intent(in):: TeTransitionRegionSi

    integer,parameter:: nStep = 1000
    integer:: iStep
    real:: DeltaTeSi, IntegralCgs, IntegralSi
    real, dimension(0:nStep):: X_I, Y_I, Simpson_I

    ! Calculate the integral, \int{\sqrt{T_e}\Lambda{T_e}dT_e}:
    !
    !--------------------------------------------------------------------------
    DeltaTeSi = (TeTransitionRegionSi - TeChromosphereSi) / nStep

    X_I(0) = TeChromosphereSi
    call get_cooling_function_fit(X_I(0), Y_I(0))
    Simpson_I(0) = 2.0

    do iStep = 1,nStep
       X_I(iStep) = X_I(iStep - 1) + DeltaTeSi
       call get_cooling_function_fit(X_I(iStep), Y_I(iStep))
       Simpson_I(iStep) = 6.0 - Simpson_I(iStep - 1)
    end do
    Simpson_I(0) = 1.0; Simpson_I(nStep) = 1.0
    Y_I = Y_I * sqrt(X_I)
    IntegralCgs = sum(Simpson_I*Y_I) * DeltaTeSi / 3.0
    !\
    ! Transform to SI. Account for the transformation coefficient for n_e^2
    !/
    IntegralSi = IntegralCgs * Cgs2SiEnergyDens * 1.0e-12
    cooling_function_integral_si = IntegralSi
  end function cooling_function_integral_si
  !============================================================================
  subroutine add_chromosphere_heating(TeSi_C,iBlock)
    use ModGeometry, ONLY: r_BLK
    use ModConst,    ONLY: mSun, rSun, cProtonMass, cGravitation, cBoltzmann
    use ModPhysics,  ONLY: UnitX_, Si2No_V,UseStar,RadiusStar,MassStar
    use ModCoronalHeating, ONLY: CoronalHeating_C

    real,    intent(in):: TeSi_C(1:nI,1:nJ,1:nK)
    integer, intent(in):: iBlock

    real, parameter:: cGravityAcceleration = cGravitation * mSun / rSun**2
    real, parameter:: cBarometricScalePerT = cBoltzmann/&
         (cProtonMass * cGravityAcceleration)

    real:: HeightSi_C(1:nI,1:nJ,1:nK), BarometricScaleSi, Amplitude
    character(len=*), parameter:: NameSub = 'add_chromosphere_heating'
    !--------------------------------------------------------------------------
    HeightSi_C = (r_BLK(1:nI,1:nJ,1:nK,iBlock) - 1) * Si2No_V(UnitX_)
    if(UseStar)then
       BarometricScaleSi = TeChromosphereSi *  cBarometricScalePerT*RadiusStar**2/MassStar
    else
       BarometricScaleSi = TeChromosphereSi *  cBarometricScalePerT
    endif
    Amplitude = radiative_cooling(TeChromosphereSi, NumberDensChromosphereCgs)

    where(HeightSi_C < 10.0 * BarometricScaleSi&
         .and.TeSi_C < 1.1 * TeChromosphereSi)&
         CoronalHeating_C = CoronalHeating_C + &
         Amplitude * exp(-HeightSi_C/BarometricScaleSi)

  end subroutine add_chromosphere_heating
  !============================================================================
  subroutine calc_reb_density(iSide, iFace, jFace, kFace, iBlock,&
       IsNewBlock, TotalFaceB_D,                                 &
       Te_G, TeTRTopSiIn, TeChromoSiIn, RadIntegralSiIn, DensityReb)
    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModConst, ONLY: kappa_0_e
    use ModCoronalHeating, ONLY: get_cell_heating
    use ModPhysics, ONLY: Si2No_V, No2Si_V, UnitX_, UnitT_,&
                          UnitEnergyDens_, UnitN_, UnitTemperature_
    use ModFaceGradient, ONLY: get_face_gradient
    ! function to return the density given by the Radiative Energy Balance Model
    ! (REB) for the Transition region. Originally given in Withbroe 1988, this
    ! uses eq from Lionell 2001. NO enthalpy flux correction in this
    ! implementation.
    integer, intent(in):: iSide,iFace,jFace,kFace,iBlock

    logical, intent(inout):: IsNewBlock

    real, intent(in):: TotalFaceB_D(3)

    real,intent(inout) :: Te_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK)

    ! Temperature value in K chosen as the temperature at the top of
    ! the transition region
    real, intent(in),optional :: TeTRTopSiIn
    ! The same for
    ! the bottom of the transition region=the top of the chromosphere
    real, intent(in),optional :: TeChromoSiIn
    ! The value of integlar of the radiation loss function time \sqrt{T}
    ! from TeChromoSi till TeTRTopSi
    real, intent(in),optional :: RadIntegralSiIn
    real, intent(out) :: DensityReb
    real :: TeTRTopSi     = 1.0e4
    real :: TeChromoSi    = 5.0e5
    ! Here Rad integral is integral of lossfunction*T^(1/2) from T=10,000 to
    ! 500,000. Use same approximate loss function used in BATS to calculate
    ! This is in SI units [J m^3 K^(3/2)]
    real :: RadIntegralSi = 1.009E-26

    real :: FaceGrad_D(3)
    real :: TotalFaceBunit_D(3)

    ! Left and right cell centered heating
    real :: CoronalHeatingLeft, CoronalHeatingRight, CoronalHeating

    ! Condensed terms in the REB equation
    real :: qCondSi, qHeatSi

    integer :: iDir=0

    character(len=*), parameter:: NameSub = 'calc_reb_density'
    !--------------------------------------------------------------------------
    if(present(TeChromoSiIn))then
       TeChromoSi = TeChromoSiIn
    else
       TeChromoSi = 1.0e4
    end if

    if(present(TeTRTopSiIn))then
       TeTRTopSi = TeTRTopSiIn
    else
       TeTRTopSi = 5.0e5
    end if

    if(present(RadIntegralSiIn))then
       RadIntegralSi = RadIntegralSiIn
    else
       ! Here Rad integral is integral of lossfunction*T^(1/2) from T=10,000 to
       ! 500,000. Use same approximate loss function used in BATS to calculate
       ! This is in SI units [J m^3 K^(3/2)]
       RadIntegralSi = 1.009E-26
    end if

    ! need to get direction for face gradient calc
    ! also put left cell centered heating call here (since index depends on
    ! the direction)
    if(iSide==1 .or. iSide==2) then
       iDir = x_
       call get_cell_heating(iFace-1, jFace, kFace, iBlock, CoronalHeatingLeft)
    elseif(iSide==3 .or. iSide==4) then
       iDir = y_
       call get_cell_heating(iFace, jFace-1, kFace, iBlock, CoronalHeatingLeft)
    elseif(iSide==5 .or. iSide==6) then
       iDir = z_
       call get_cell_heating(iFace, jFace, kFace-1, iBlock, CoronalHeatingLeft)
    else
       call stop_mpi('REB model got bad face direction')
    endif

    call get_cell_heating(iFace, jFace, kFace, iBlock, CoronalHeatingRight)

    CoronalHeating = 0.5 * (CoronalHeatingLeft + CoronalHeatingRight)

    ! term based on coronal heating into trans region (calc face centered avg)
    qHeatSi = (2.0/7.0) * CoronalHeating * TeTRTopSi**1.5 &
         * No2Si_V(UnitEnergyDens_) / No2Si_V(UnitT_)

    call get_face_gradient(iDir, iFace, jFace, kFace, iBlock, &
         IsNewBlock, Te_G, FaceGrad_D)

    ! calculate the unit vector of the total magnetic field
    TotalFaceBunit_D = TotalFaceB_D / sqrt(sum(TotalFaceB_D**2))

    ! calculate the heat conduction term in the REB numerator
    qCondSi = 0.5 * kappa_0_e(20.) * TeTRTopSi**3 &
         * sum(FaceGrad_D*TotalFaceBunit_D)**2 &
         * (No2Si_V(UnitTemperature_) / No2Si_V(UnitX_))**2

    ! put the terms together and calculate the REB density
    DensityReb = sqrt((qCondSi + qHeatSi) / RadIntegralSi) &
         * Si2No_V(UnitN_)

  end subroutine calc_reb_density
  !============================================================================
end module ModRadiativeCooling
!==============================================================================
