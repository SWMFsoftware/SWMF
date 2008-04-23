!^CFG COPYRIGHT UM
!!!=====================Heat capacities, at constant volume==============!!!!
!!
!! Thermodynamical variables and other notations
!!         \rho, Rho - the mass density
!!         {\cal E}, E - internal energy of the unit of mass
!!         e, i - electron, ion
!!        V, vol - volume or volumetric
!!        \left(\frac{\partial...}{\partial...}\right)_V - thermodynamical
!!             derivative at constant volume
!!        T_{e,i}, Te,Ti - electron and ion temperature
!!        iLevel - integer variable, a signature of the material:
!!        iLevel=0 - xenon
!!        ilevel=1 - berillium         
!!!===================== Electron part ==================================!!!!
!!       
!!!====================  Heat capacity per volume:=======================!!!!
!!               Rho\left(\frac{\partial E}{\partial T_e}\right)_{V,T_i} !!!!
subroutine get_e_heat_capacity_per_vol(HeatCapacity,Te,Rho,iLevel)
  implicit none
  real,intent(out)::HeatCapacity  !In SI: J/(K m^3)
  real,intent(in):: Te            !Electron temperature in SI: K
  real,intent(in):: Rho           !Mass density in SI: kg/m^3
  integer,intent(in)::iLevel      !iLevel=0 - xenon, ilevel=1 - berillium
  !-----------------------------------------------------------------------!
  real::HeatCapacityPerMass
  !-------------------------!

  call get_e_heat_capacity_per_mass(HeatCapacityPerMass,Te,Rho,iLevel)
  HeatCapacity = HeatCapacityPerMass * Rho

end subroutine get_e_heat_capacity_per_vol
!!=======================================================================!!!!
!!
!!!====================  Heat capacity per mass:  =======================!!!!
!!               \left(\frac{\partial E}{\partial T_e}\right)_{V,T_i}    !!!!
subroutine get_e_heat_capacity_per_mass(HeatCapacity,Te,Rho,iLevel)      !!!!
  use ModConst
  implicit none
  real,intent(out)::HeatCapacity  !In SI: J/(K m^3)
  real,intent(in):: Te            !Electron temperature in SI: K
  real,intent(in):: Rho           !Mass density in SI: kg/m^3
  integer,intent(in)::iLevel      !iLevel=0 - xenon, ilevel=1 - berillium
  !----------------------------------------------------------------------!
  real,parameter :: cAtomicToMass = cBoltzmann / cAtomicMass 
  real,parameter,dimension(0:1) :: AtomicMass_I=(/131.29, &
                                                  9.012 /)
  integer,parameter,dimension(0:1) :: nZ_I=(/54, &
                                              4/)
  real::HeatCapacityPerAtom
  !----------------------------------------------------------------------!
  !        Version for fully ionized plasma                              !
  HeatCapacityPerAtom = 1.50*real(nZ_I(iLevel))
  !----------------------------------------------------------------------!
  HeatCapacity = HeatCapacityPerAtom * cAtomicToMass / &
                     AtomicMass_I(iLevel)
end subroutine get_e_heat_capacity_per_mass
!============================================================================!
