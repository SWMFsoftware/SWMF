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
!           Numerical integration is performed using Simpson's Rule:
!           S = (h/3) * [ f(a) + f(b) + 2 Sum[i=1,n-1](f(x_2i))
!                         + 4 Sum[i=1,n](f(x_{2i-1})) ]
!
! References: http://www.salihnet.freeservers.com/engineering/
!                     fortran_codes/unequal_simps.html
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

  integer :: m, n, i
  real :: sum1, sum2, sum3
  real, dimension(nAlts) :: height

  ! Perform a simple numerical integration, using Simpson's rule for unequally
  ! spaced data.  VTEC is in TECU while electron density is in m^-3

  ! Evaluate the segment width

  m = nAlts - 1

  do i = 1, m
     height(i) = Altitude_GB(iLon,iLat,i+1,iBlock) &
          - Altitude_GB(iLon,iLat,i,iBlock)
  enddo

  n = nAlts + 1

  do i = nAlts, n
     height(i) = 0.0
  enddo

  sum1 = 0.0
  sum2 = 0.0
  sum3 = 0.0
  i    = 1

  ! Loop for integration

  do while(i < nAlts)
     if((height(i) == height(i+1)).and.(height(i)==height(i+2))) then
        ! Simpson's 3/8 rule
        sum1 = sum1 + (3.0 * height(i) * (IDensityS(iLon,iLat,i,ie_,iBlock) &
             + 3.0 * (IDensityS(iLon,iLat,i+1,ie_,iBlock) &
             + IDensityS(iLon,iLat,i+2,ie_,iBlock)) &
             + IDensityS(iLon,iLat,i+3,ie_,iBlock))) / 8.0
        i = i + 3
     elseif(height(i) == height(i+1)) then
        ! Simpson's 1/3 rule
        sum2 = sum2 + (2.0 * height(i) * (IDensityS(iLon,iLat,i,ie_,iBlock) &
             + 4.0 * IDensityS(iLon,iLat,i+1,ie_,iBlock) &
             + IDensityS(iLon,iLat,i+2,ie_,iBlock))) / 6.0
        i = i + 2
     elseif(height(i).ne.height(i+1)) then
        ! Trapezoidal rule
        sum3 = sum3 + height(i) * (IDensityS(iLon,iLat,i,ie_,iBlock) &
             + IDensityS(iLon,iLat,i+1,ie_,iBlock)) / 2.0
        i = i + 1
     endif
  enddo

  single_vtec = (sum1 + sum2 + sum3) * (10.0**(-16.0))

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
