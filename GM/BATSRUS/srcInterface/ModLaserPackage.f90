!=================================THE LASER PACKAGE STARTS HERE=================!
!Our intent is to make ray_path routine re-usable,
!which is a reason to put the package here.
!===============================================================================!
module ModLaserPulse
  implicit none
  !The laser pulse is assumed to have the linear raise front
  !and linear decay with the constant irradiance within the
  !pulse,  at  t_raise < t < tPulse - tDecay
  real:: tRaise = 1.0e-10 ![s]
  real:: tDecay = 1.0e-10 ![s]
  real:: tPulse = 1.1e-9  ![s]
  real:: IrradianceSi = 3.8e12 ![J/s]
contains
  real function irradiance_t(TimeSi)
    real,intent(in)::TimeSi
    !----------------------
    irradiance_t = IrradianceSi *&
         max(0.0, min(1.0,       &
         TimeSi/tRaise,          &
         (tPulse -TimeSi)/tDecay))
         
  end function irradiance_t
end module ModLaserPulse
!===========================
subroutine read_laser_pulse_param
  use ModReadparam
  use ModLaserPulse
  implicit none
  !-------------
  !TO DO !!!
end subroutine read_laser_pulse_param
!===========================
module ModLaserPackage
  use ModImplicit, ONLY: ResExpl_VCB, iTeImpl
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  use ModRadioWaveRaytracing, ONLY: ray_path
  implicit none

  real:: Irradiance
  !------------
contains
  subroutine update_impl_energy_source
  end subroutine update_impl_energy_source
end module ModLaserPackage

subroutine add_laser_energy_deposition
  !This routine should add the laser energy deposition to the "Explicit residual"
  !within the semi-implicit scheme. Effectively, in this way the heat conduction
  !equation, which otherwise is solved within the semi-implicit scheme WITH ZERO
  !right-hand-side, is solved with the right-hand-side including the laser energy
  !deposition

  !use ModImplicit, ONLY: ResExpl_VCB, iTeImpl
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  use ModPhysics, ONLY: Si2No_V, UnitEnergydens_, UnitX_, UnitT_
  use ModLaserPackage, ONLY: Irradiance, update_impl_energy_source
  use ModMain, ONLY: Time_Simulation
  use ModLaserPulse, ONLY: irradiance_t
  implicit none
  !----------------
  Irradiance = irradiance_t(Time_Simulation) * &
       Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  call update_impl_energy_source
  
end subroutine add_laser_energy_deposition
!==========================================
