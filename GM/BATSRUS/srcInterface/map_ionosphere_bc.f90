!^CFG COPYRIGHT UM


subroutine Get_MagneticField_Orient(qx,qb)

  implicit none

  ! Obtain normalized bb field at true location qx and put it into qb

  real, intent(in) :: qx(3)
  real, intent(out):: qb(3)

  ! Get B0

  call get_b0(qx(1),qx(2),qx(3),qb)

  ! Normalize
  qb=qb/sqrt(sum(qb**2))

end subroutine Get_MagneticField_Orient

!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------

subroutine Get_Mapping_Point(XyzIn_D, DirectionIn, XyzOut_D)

  ! Given the XyzIn_D coordinates in the magnetosphere system, 
  ! find the mapping point along the magnetic field lines for
  ! the hemisphere signaled by HemisphereIn (+1.0 for north, -1.0 for south).
  ! The coordinates of the mapping point will be in XyzOut_D in the
  ! ionosphere coordinate system. If HemisphereIn differs from the
  ! hemisphere where XyzIn_D is, then XyzOut_D is set to 0.0.

  use ModMain, ONLY: TypeCoordSystem, Time_Simulation
  use ModMappingParam, ONLY: rMap
  use CON_physics, ONLY: map_planet_field, transform_matrix

  implicit none

  real, intent(in)  :: XyzIn_D(3)
  real, intent(in)  :: DirectionIn
  real, intent(out) :: XyzOut_D(3)

  real :: XyzSmg_D(3) ! Coordinates is solar magnetic system

  integer    :: iHemisphere
  !----------------------------------------------------------------------------


  ! Transform input coordinates to SMG system (used by IE)
  XyzSmg_D = matmul(transform_matrix(Time_Simulation,TypeCoordSystem,'SMG'),&
       XyzIn_D)

  ! Check if the initial point is on the required hemisphere
  if( XyzSmg_D(3)*DirectionIn < 0.0) then
     ! We are on the wrong hemisphere, return 0. location
     XyzOut_D = 0.0
     RETURN
  end if

  ! Get the mapping point
  call map_planet_field(Time_Simulation,XyzSmg_D,'SMG NORM',rMap, &
       XyzOut_D,iHemisphere)

end subroutine Get_Mapping_Point
!==============================================================================

!-------------------------------------------------------------------------
! Get_Mapping_Information
!
subroutine GM_get_mapping_param_for_ie(rCurrentsGm,rIonosphere)
  use ModPhysics,ONLY:rCurrents
  use ModMappingParam,ONLY:rMap
  implicit none
  real,intent(out)::rCurrentsGM,rIonosphere
  rCurrentsGm=rCurrents;rIonosphere=rMap
end subroutine GM_get_mapping_param_for_ie
