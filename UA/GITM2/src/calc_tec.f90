!
!-----------------------------------------------------------------------------
! $Id$
!
! calc_tec
!
! Author: Angeline G. Burrell, UMichigan, December 2012
!
! Comments: Routines to compute the total electron content (TEC) from
!           the GITM electron density profiles.  TEC is defined as the integral
!           of the electron density from the receiver to the satellite.
!           VTEC (vertical TEC) integrates straight up from the ground up to
!           the GPS satellite orbit height of ~20,200 km.  As GITM reaches up
!           to ~600 km, the VTEC computed here should be less than the VTEC 
!           measured by ground-based receivers.  Extending the topside 
!           ionosphere and protonosphere would make up for much of the 
!           difference.  The plasmapause (above 2,000 km) typically contributes
!           less than 5% to the VTEC.  TEC is given in TECU, where 
!           1 TECU = 10^16 m^-2
!
!           Numerical integration is performed using the Trapazoidal Rule
!
! Contains: calc_single_vtec <- computes the VTEC at a specified latitude
!                               and longitude, using the latitude, longitude
!                               and block indeces
!           calc_vtec <- computes the VTEC at all locations
!-----------------------------------------------------------------------------

subroutine calc_single_vtec(iLon, iLat, iBlock, single_vtec)

  use ModGITM

  implicit none

  integer, intent(in) :: iLon, iLat, iBlock
  real, intent(out) :: single_vtec

  integer :: m, i
  real :: height, density

  ! Perform a simple numerical integration, using the trapazoidal rule.  VTEC
  ! is in TECU while electron density is in m^-3

  m           = nAlts - 1
  single_vtec = 0.0

  do i = 1, m
     ! Calculate height incriment and average density
     height  = Altitude_GB(iLon,iLat,i+1,iBlock)-Altitude_GB(iLon,iLat,i,iBlock)
     density = 0.5 * (IDensityS(iLon,iLat,i,ie_,iBlock) &
          + IDensityS(iLon,iLat,i+1,ie_,iBlock))

     ! Sum successive incrimentations of height * density
     single_vtec = single_vtec + height * density
  enddo

  ! Convert from SI units to TEC units

  single_vtec = single_vtec * (10.0**(-16.0))
end subroutine calc_single_vtec


subroutine calc_vtec(iBlock)

  use ModGITM

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat

  interface
     subroutine calc_single_vtec(iLon, iLat, iBlock, single_vtec)
       integer, intent(in) :: iLon, iLat, iBlock
       real, intent(out) :: single_vtec
     end subroutine calc_single_vtec
  end interface

  ! Perform a simple numerical integration, summing the electron density
  ! and multiplying it by the altitude range.  VTEC is in TECU while
  ! electron density is in m^-3

  do iLon=-1,nLons+2
     do iLat=-1,nLats+2
        call calc_single_vtec(iLon, iLat, iBlock, VTEC(iLon,iLat,iBlock))
     enddo
  enddo

  return
end subroutine calc_vtec
