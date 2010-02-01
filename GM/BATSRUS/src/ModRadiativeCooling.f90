!^CFG COPYRIGHT UM
!============================================================================
module ModRadiativeCooling
  implicit none
  save
  
  logical:: UseRadCooling      =.false.
  logical:: UseRadCoolingTable =.false.
  integer:: iTableRadCool      =-1

  ! Same use as the variables in ModHeatConduction
  ! but need them here as well
  logical :: DoModHeatConduction = .false.

  real :: TeModSi = 3.0E+5
  real :: DeltaTeModSi = 1E+4
  real :: TeMod, DeltaTeMod
contains
  subroutine get_radiative_cooling(i, j, k, iBlock, RadiativeCooling)
    
    use ModMultiFluid,    ONLY: MassIon_I
    use ModPhysics,       ONLY: No2Si_V, Si2No_V, UnitT_, UnitN_, &
                                UnitEnergyDens_, UnitTemperature_
    use ModVarIndexes,    ONLY: nVar, Rho_, P_
    use ModLookupTable,   ONLY: interpolate_lookup_table
    use ModAdvance,       ONLY: State_VGB
    use ModUser,          ONLY: user_material_properties
    integer, intent(in) :: i, j, k, iBlock
    real, intent(out):: RadiativeCooling
    
    real :: Te, TeSi, CoolingFunctionCgs, NumberDensCgs
    real :: Log10TeSi, Log10NeCgs, CoolingTableOut(2)
    real, parameter :: RadNorm = 1.0E+22
    real, parameter :: TeModMinSi = 2.0E+4
    real :: TeFactor, TeModMin, FractionSpitzer
    !--------------------------------------------------------------------------
    
    call user_material_properties(State_VGB(:, i, j, k, iBlock),TeOut=TeSI)
    
    ! currently proton plasma only
    NumberDensCgs = State_VGB(Rho_, i, j, k, iBlock) * No2Si_V(UnitN_) * 1.0e-6
    if(UseRadCoolingTable) then
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
    
    RadiativeCooling = - NumberDensCgs**2*CoolingFunctionCgs &
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
end module ModRadiativeCooling
