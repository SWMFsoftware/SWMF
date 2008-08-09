!^CFG COPYRIGHT UM
!=====================Equation Of State (EOS)===========================!!!!
!
! Thermodynamical variables and other notations
!         \rho, Rho - the mass density
!         {\cal E}, E - internal energy of the unit of mass
!         e, i - electron, ion
!        V, vol - volume or volumetric
!        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
!             derivative at constant volume
!        T_{e,i}, Te,Ti - electron and ion temperature
!        iMaterial - integer variable, a signature of the material:
!        iMaterial=0 - xenon
!        iMaterial=1 - beryllium         
module ModEos
  use ModStatSum
  use ModAtomicMass
  implicit none
  include 'CRASH_definitions.h'
contains
  subroutine eos(UDensityTotal,& !Input total energy density SI,[J/m^3]
                 Rho,          & !Input mass density, SI [kg/m^3] 
                 iMaterial,    & !Input: sort of material
                 TeOut,        & !Output, OPTIONAL, temperature SI[K]
                 PTotalOut,    & !Output, OPTIONAL, pressure, SI [Pa]
                 GammaOut,     & !Output, OPTIONAL, polytropic index
                 Energy0Out)     !Output   (E-P/(\gamma-1))/\rho
    real,intent(in)::UDensityTotal,Rho
    integer,intent(in)::iMaterial
    real,optional,intent(out)::TeOut,PTotalOut,GammaOut,Energy0Out
    real::NAtomic, UPerAtom
    logical::IsDegenerated
    !----------------------------------------------------------------------!
    call set_element(nZ_I(iMaterial))

    !Get the atomic concentration
    NAtomic=Rho/(cAtomicMass*cAtomicMass_I(nZ_I(iMaterial)))

    !Get an energy per the atomic cell, express in eV
    UPerAtom = UDensityTotal / (cEV * NAtomic)
    
    !Find temperature from dentity and internal energy
    call set_temperature(UPerAtom,NAtomic,IsDegenerated)

    if(IsDegenerated)call CON_stop(&
         'No EOS for Fermi=degenerated state')
    if(present(TeOut))TeOut = Te*cEVToK
    if(present(PTotalOut))PTotalOut = pressure()
    if(present(GammaOut))then
       call get_termodyn_derivatives(GammaOut=GammaOut)
       if(present(Energy0Out))&
          Energy0Out = (UDensityTotal - pressure()/(GammaOut-1.0))/Rho
    end if
  end subroutine eos
end module ModEos
