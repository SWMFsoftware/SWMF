!^CFG COPYRIGHT UM
!============================================================================
module ModChromosphere
  !Here all parameters relating to chromosphere and tansition region are
  !collected
  
  implicit none
  save
  
  !The use of short-scale exponential heat function with 
  !the decay length = (30 m/K)*TeCromosphere SI
  logical:: UseChromosphereHeating    = .false. 
  real   :: NumberDensChromosphereCgs = 1.0e+12 ![cm^{-3}  
  real   :: TeChromosphereSi = 1.0e4            ![K]

  !\
  ! TRANSITION REGION
  !/
  ! 
  
  logical :: DoExtendTransitionRegion = .false.

  !The following variables are meaningful if
  !DoExtendTransitionRegion=.true.

  real :: TeModSi = 3.0E+5                !K
  real :: DeltaTeModSi = 1E+4             !K

  !The following variable is meaningful if 
  !DoExtendTransitionRegion = .false. . Al long as
  !the unextended transition region cannot be resolved
  !we set the 'corona base temperature' equal to the
  !temperature at the top of the transition region and 
  !use the integral ralationship across the transition
  !region to find the number density at the top of the
  !transition region

  real :: TeTransitionRegionTopSi = 4.0e+5 ![K]

contains
  !================================

  subroutine read_chromosphere
    use ModReadParam, ONLY: read_var
    !-------------------------------
    call read_var('UseChromosphereHeating'   , UseChromosphereHeating)
    call read_var('NumberDensChromosphereCgs', NumberDensChromosphereCgs)
    call read_var('TeChromosphereSi',          TeChromosphereSi      )
  end subroutine read_chromosphere

  !================================

  real function extension_factor(TeSi)
    real, intent(in) :: TeSi    !Dimensionless
    
    real :: FractionSpitzer
    !--------------------------------
    FractionSpitzer = 0.5*(1.0+tanh((TeSi-TeModSi)/DeltaTeModSi))
    
    extension_factor = FractionSpitzer + &
         (1.0 - FractionSpitzer)*(TeModSi/TeSi)**2.5
  end function extension_factor
  
end module ModChromosphere
!=========================
module ModRadiativeCooling
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
                               1.0e-7&   !erg = 1e-7 J
                              /1.0e-6    !cm3 = 1e-6 m3 
contains
  !==============================
  subroutine read_modified_cooling
    use ModReadParam
    !---------------
    call read_var('DoExtendTransitionRegion',DoExtendTransitionRegion)
    if(DoExtendTransitionRegion)then

       call read_var('TeModSi',TeModSi)
       call read_var('DeltaTeModSi',DeltaTeModSi)
    else
       call read_var('TeTransitionRegionTopSi',TeTransitionRegionTopSi)
    end if
  end subroutine read_modified_cooling
  !==============================

  subroutine check_cooling_param
    use ModPhysics,ONLY: Si2No_V, UnitTemperature_
    use ModLookupTable, ONLY: i_lookup_table
    !-----------------------

    iTableRadCool = i_lookup_table('radcool')
    UseRadCoolingTable = iTableRadCool>0

  end subroutine check_cooling_param
  !==============================
  subroutine get_radiative_cooling(i, j, k, iBlock, TeSiIn, RadiativeCooling)

    use ModPhysics,       ONLY: Si2No_V, No2Si_V, UnitN_
    use ModVarIndexes,    ONLY: Rho_
    use ModAdvance,       ONLY: State_VGB

    integer, intent(in) :: i, j, k, iBlock
    real,    intent(in) :: TeSiIn

    real,   intent(out) :: RadiativeCooling


    real :: NumberDensCgs

    !--------------------------------------------------------------------------



    ! currently proton plasma only
    NumberDensCgs = State_VGB(Rho_, i, j, k, iBlock) * No2Si_V(UnitN_) * 1.0e-6

    RadiativeCooling = - radiative_cooling(TeSiIn, NumberDensCgs)

    ! include multiplicative factors to make up for extention of
    ! perpendicular heating at low temperatures (as per Abbett 2007).
    ! Need this to strech transition region to larger scales
    ! Also, need radcool modification to become const below TeModMin
    if(DoExtendTransitionRegion) then
       if(TeSiIn >= TeChromosphereSi) then
          RadiativeCooling = RadiativeCooling / extension_factor(TeSiIn)
       else
          RadiativeCooling = RadiativeCooling * (TeChromosphereSi / TeModSi)**2.5
       endif
    end if
  end subroutine get_radiative_cooling
  !============================================================
  real function radiative_cooling(TeSiIn, NumberDensCgs)
    use ModMultiFluid,    ONLY: MassIon_I
    use ModPhysics,       ONLY: Si2No_V, UnitT_, &
         UnitEnergyDens_
    use ModLookupTable,   ONLY: interpolate_lookup_table

    !Imput - dimensional, output: dimensionless
    real, intent(in):: TeSiIn, NumberDensCgs
 
    real :: CoolingFunctionCgs
    real :: Log10TeSi, Log10NeCgs, CoolingTableOut(2)
    real, parameter :: RadNorm = 1.0E+22
    !-----------------------------------
    if(UseRadCoolingTable) then
       Log10TeSi = log10(TeSiIn)
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
       call get_cooling_function_fit(TeSiIn, CoolingFunctionCgs)
    end if

    ! Avoid extrapolation past zero
    CoolingFunctionCgs = max(CoolingFunctionCgs,0.0)

    radiative_cooling = NumberDensCgs**2*CoolingFunctionCgs &
         * Cgs2SiEnergyDens * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

  end function radiative_cooling
!===============================
  subroutine get_cooling_function_fit(TeSi, CoolingFunctionCgs)
    
    ! Based on Rosner et al. (1978) and Peres et al. (1982)
    ! Need to be replaced by Chianti tables
    
    real, intent(in) :: TeSi
    real, intent(out):: CoolingFunctionCgs
    !-----------------------------------------------------------
    
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
  !=============================================================
  real function cooling_function_integral_si(TeTransitionRegionSi)
    
    real,intent(in):: TeTransitionRegionSi
    
    integer,parameter:: nStep = 1000
    integer:: iStep
    real:: DeltaTeSi, IntegralCgs, IntegralSi
    real, dimension(0:nStep):: X_I, Y_I, Simpson_I
    !-----------------
    !Calculate the integral, \int{\sqrt{T_e}\Lambda{T_e}dT_e}:
    !
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
    !Transform to SI. Account for the transformation coefficient for n_e^2
    !/
    IntegralSi = IntegralCgs * Cgs2SiEnergyDens * 1.0e-12 
    cooling_function_integral_si = IntegralSi
  end function cooling_function_integral_si

end module ModRadiativeCooling
