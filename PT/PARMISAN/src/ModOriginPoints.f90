!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModOriginPoints
  use ModNumConst,       ONLY: cDegToRad
  implicit none
  ! Grid size, boundaries, coordinates
  ! Starting position of field lines in Rs
  real         :: ROrigin = 2.5
  ! Size of angular grid, latitude and longitude, at origin
  ! surface R=ROrigin
  real         :: LonMin =   0.0*cDegToRad
  real         :: LonMax = 360.0*cDegToRad
  real         :: LatMin = -70.0*cDegToRad
  real         :: LatMax =  70.0*cDegToRad
contains
  !============================================================================
  subroutine read_param
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var
    use CON_axes,     ONLY: dLongitudeHgrDeg

    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    call read_var('ROrigin', ROrigin)
    call read_var('LonMin' , LonMin )
    call read_var('LatMin' , LatMin )
    call read_var('LonMax' , LonMax )
    call read_var('LatMax' , LatMax )

    ! Correct for Rotation
    LonMax = LonMax - dLongitudeHgrDeg
    LonMin = LonMin - dLongitudeHgrDeg

    !
    ! convert angels from degrees to radians

    LonMax = LonMax*cDegToRad
    LonMin = LonMin*cDegToRad

    ! convert angels from degrees to radians

    LatMax = LatMax*cDegToRad
    LatMin = LatMin*cDegToRad
  end subroutine read_param
  !============================================================================
end module PT_ModOriginPoints
!==============================================================================
