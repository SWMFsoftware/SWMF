!^CFG COPYRIGHT UM
!============================================================================
module ModRadiativeCooling
  use ModSize
  implicit none
  save

  logical:: UseRadCooling      =.false.
  logical:: UseRadCoolingTable =.false.
  integer,private:: iTableRadCool      =-1

  ! Same use as the variables in ModHeatConduction
  ! but need them here as well
  logical :: DoModHeatConduction = .false.

  real :: TeModSi = 3.0E+5                !K
  real :: DeltaTeModSi = 1E+4             !K
  real :: TeMod=0.0, DeltaTeMod=1.0

  real, parameter :: TeModMinSi = 2.0E+4  !K
  real ::  TeModMin

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
contains
  !==============================
  subroutine read_modified_cooling
    use ModReadParam
    !---------------
    call read_var('DoModHeatConduction',DoModHeatConduction)
    if(.not.DoModHeatConduction)return

    call read_var('TeModSi',TeModSi)
    call read_var('DeltaTeModSi',DeltaTeModSi)
  end subroutine read_modified_cooling
  !==============================

  subroutine check_cooling_param
    use ModPhysics,ONLY: Si2No_V, UnitTemperature_
    use ModLookupTable, ONLY: i_lookup_table
    !-----------------------

    iTableRadCool = i_lookup_table('radcool')
    UseRadCoolingTable = iTableRadCool>0

    TeMod      = TeModSi      * Si2No_V(UnitTemperature_)
    DeltaTeMod = DeltaTeModSi * Si2No_V(UnitTemperature_)
    TeModMin   = TeModMinSi   * Si2No_V(UnitTemperature_)

  end subroutine check_cooling_param
  !==============================
  subroutine get_radiative_cooling(i, j, k, iBlock, TeSiIn, RadiativeCooling)

    use ModMultiFluid,    ONLY: MassIon_I
    use ModPhysics,       ONLY: No2Si_V, Si2No_V, UnitT_, UnitN_, &
         UnitEnergyDens_, UnitTemperature_
    use ModVarIndexes,    ONLY: nVar, Rho_, P_
    use ModLookupTable,   ONLY: interpolate_lookup_table
    use ModAdvance,       ONLY: State_VGB

    integer, intent(in) :: i, j, k, iBlock
    real,    intent(in) :: TeSiIn

    real,   intent(out) :: RadiativeCooling


    real :: Te, CoolingFunctionCgs, NumberDensCgs
    real :: Log10TeSi, Log10NeCgs, CoolingTableOut(2)
    real, parameter :: RadNorm = 1.0E+22
    real, parameter :: Cgs2SiEnergyDens = &
                               1.0e-7&   !erg = 1e-7 J
                              /1.0e-6    !cm3 = 1e-6 m3 
    real :: TeFactor, FractionSpitzer
    !--------------------------------------------------------------------------



    ! currently proton plasma only
    NumberDensCgs = State_VGB(Rho_, i, j, k, iBlock) * No2Si_V(UnitN_) * 1.0e-6

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

    RadiativeCooling = - NumberDensCgs**2*CoolingFunctionCgs &
         * Cgs2SiEnergyDens * Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)

    ! include multiplicative factors to make up for extention of
    ! perpendicular heating at low temperatures (as per Abbett 2007).
    ! Need this to strech transition region to larger scales
    ! Also, need radcool modification to become const below TeModMin
    if(DoModHeatConduction) then
       Te = TeSiIn * Si2No_V(UnitTemperature_)

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
end module ModRadiativeCooling
