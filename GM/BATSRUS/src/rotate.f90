!^CFG COPYRIGHT UM
subroutine calc_corotation_velocities(iter,time_now,Xyz_D,uRot_D)
  !-------------------------------------------------------------------------
  !\
  ! This routine calculates cartesian corotation velocity uRot_D as a
  ! function of the cartesian coordinates Xyz_D
  !/
  !-------------------------------------------------------------------------

  use CON_axes,          ONLY: get_axes
  use ModCompatibility,  ONLY: &
       compatible_rotation => calc_corotation_velocities
  use ModCoordTransform, ONLY: cross_product
  use ModMain,           ONLY: Time_Simulation, TypeCoordSystem, UseNewAxes
  use ModPhysics,        ONLY: OmegaBody
  use ModNumConst
  implicit none

  integer, intent(in) :: iter
  integer, intent(in) :: time_now

  real, intent(in) :: Xyz_D(3)
  real, intent(out):: uRot_D(3)

  real    :: Omega_D(3)
  logical :: IsUninitialized = .true.

  !------------------------------------------------------------------------
  if(.not.UseNewAxes)then
     ! The user decided to call an overcomplicated, incorrect and inefficient 
     ! algorithm. Enjoy...
     call compatible_rotation(iter,time_now,Xyz_D,uRot_D)
     RETURN
  end if

  select case(TypeCoordSystem)
  case('HGI')
     ! In the HGI system the Solar angular velocity vector points towards +Z
     Omega_D = (/ 0., 0., OmegaBody /)
  case('GSE')
     if(IsUninitialized)then
        call get_axes(Time_Simulation,RotAxisGseOut_D=Omega_D)
        Omega_D = OmegaBody * Omega_D
        IsUninitialized = .false.
     end if
  case('GSM')
     ! GSM system, Omega_D may be changing
     call get_axes(Time_Simulation,RotAxisGsmOut_D=Omega_D)
     Omega_D = OmegaBody*Omega_D
  end select

  ! The corotation velocity is u = Omega x R

  uRot_D = cross_product(Omega_D, Xyz_D)

end subroutine calc_corotation_velocities
