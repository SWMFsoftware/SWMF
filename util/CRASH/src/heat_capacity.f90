!^CFG COPYRIGHT UM
!=====================Heat capacities, at constant volume===============!!!!
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
!===================== Electron part ===================================!!!!
!       
!=====================  Heat capacity per volume:=======================!!!!
!               Rho\left(\frac{\partial E}{\partial T_e}\right)_{V,T_i} !!!!
subroutine get_e_heat_capacity_per_vol(HeatCapacity,Te,Rho,iMaterial)
  implicit none
  include 'CRASH_definitions.h'
  real,intent(out)::HeatCapacity  !In SI: J/(K m^3)
  real,intent(in):: Te            !Electron temperature in SI: K
  real,intent(in):: Rho           !Mass density in SI: kg/m^3
  integer,intent(in)::iMaterial   !iMaterial=0 - xenon,iMaterial=1 - beryllium
  !-----------------------------------------------------------------------!
  real::HeatCapacityPerMass
  !-------------------------!

  call get_e_heat_capacity_per_mass(HeatCapacityPerMass,Te,Rho,iMaterial)
  HeatCapacity = HeatCapacityPerMass * Rho

end subroutine get_e_heat_capacity_per_vol
!========================================================================!!!!
!
!=====================  Heat capacity per mass:  ========================!!!!
!                \left(\frac{\partial E}{\partial T_e}\right)_{V,T_i}    !!!!
subroutine get_e_heat_capacity_per_mass(HeatCapacity,Te,Rho,iMaterial)      !!!!

  use ModConst

  implicit none
  include 'CRASH_definitions.h'
  real,intent(out)::HeatCapacity  !In SI: J/(K m^3)
  real,intent(in):: Te            !Electron temperature in SI: K
  real,intent(in):: Rho           !Mass density in SI: kg/m^3
  integer,intent(in)::iMaterial   !iMaterial=0 - xenon, iMaterial=1 - beryllium
  !----------------------------------------------------------------------!
  real,parameter :: cAtomicToMass = cBoltzmann / cAtomicMass 
  real,parameter,dimension(0:1) :: AtomicMass_I=(/131.29, &
                                                  9.012 /)
  integer,parameter,dimension(0:1) :: nZ_I=(/54, &
                                              4/)
  real::HeatCapacityPerAtom
  !----------------------------------------------------------------------!
  !        Version for fully ionized plasma                              !
  !         Only a stab, should be 1.5*Z(n,Te)+(1.5*Te +I(Z))dZ/dTe      !
  HeatCapacityPerAtom = 1.50*real(nZ_I(iMaterial))
  !----------------------------------------------------------------------!
  HeatCapacity = HeatCapacityPerAtom * cAtomicToMass / &
                     AtomicMass_I(iMaterial)
end subroutine get_e_heat_capacity_per_mass
!========================================================================!
