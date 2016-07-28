!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
! Calculate the centrifugal acceleration dV/dt = V_{\perp}*db/dt = 
! V_{\perp} * (db/dt +V_{\|}*db/ds + V_{\perp}*\nabla b
! For details see Horwitz et al,[1994]

subroutine calc_centrifugal(nAlt, uIon_C,aCentrifugal_C)
  use ModCommonVariables, ONLY: SmLat, RAD,uJoule2, rPlanet
  use ModNumConst,        ONLY: cDegToRad
  implicit none
  
  integer, intent(in) :: nAlt
  real,    intent(in) :: uIon_C(nAlt)
  real,    intent(out):: aCentrifugal_C(nAlt)
  
  integer :: iAlt
  real    :: Lshell, uConvection
  real,allocatable :: uConvection_C(:), SmLat_C(:)
  !----------------------------------------------------------------------------
  
  if (.not. allocated(uConvection_C)) allocate(uConvection_C(nAlt))
  if (.not. allocated(SmLat_C)) allocate(SmLat_C(nAlt))

  Lshell = 1.0/(cos(SmLat*cDegToRad))**2.0
  uConvection = sqrt(uJoule2 * 1.0e4)

  do iAlt=1,nAlt
     uConvection_C(iAlt) = uConvection * (RAD(iAlt)/RAD(1))**1.5
     SmLat_C (iAlt) = acos(sqrt(RAD(iAlt)/(Lshell*rPlanet)))
     
     ! Set curvature term V_{\|}*db/ds
     aCentrifugal_C(iAlt) = uConvection_C(iAlt) * uIon_C(iAlt) &
          * 3.0 * ( 1 + (sin(SmLat_C(iAlt)))**2.0 ) * cos(SmLat_C(iAlt)) &
          / (RAD(iAlt) * ( 1 + 3.0*(sin(SmLat_C(iAlt)))**2.0 )**1.5)
     
     ! Set curvature term V_{\perp}*\dot \nabla b
     aCentrifugal_C(iAlt) = aCentrifugal_C(iAlt) + uConvection_C(iAlt)**2.0 &
          * 6.0 * ( 1 + (sin(SmLat_C(iAlt)))**2.0 ) * sin(SmLat_C(iAlt)) &
          / (RAD(iAlt) * ( 1 + 3.0*(sin(SmLat_C(iAlt)))**2.0 )**1.5)
  enddo

  if (allocated(uConvection_C)) deallocate(uConvection_C)
  if (allocated(SmLat_C)) deallocate(SmLat_C)
end subroutine calc_centrifugal
