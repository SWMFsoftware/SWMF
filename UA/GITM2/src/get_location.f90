!------------------------------------------------------------------------------
! $Id$
!
! Author: Angeline G. Burrell, UMichigan, Jan 2013
!
! LocationIndex: A routine to retireve the longitude, latitude, and block
!                indeces for a specified location.  Shamelessly stolen from
!                another place in the GITM code and put in a subroutine so that
!                it can be used in multiple places.  Exit statements were
!                added to prevent additional cycling through do-loops.
!
! Inputs: LonFind = Desired longitude
!         LatFind = Desired latitude
!
! Outputs: iiBlock = Block index containing the desired location
!          iiLon   = Longitude index for LonFind
!          iiLat   = Latitude index for LatFind
!          rLon    = Longitude interpolation scaling factor
!          rLat    = Latitude interpolation scaling factor
!------------------------------------------------------------------------------

subroutine LocationIndex(LonFind, LatFind, iiBlock, iiLon, iiLat, rLon, rLat)

  use ModGITM

  real, intent(in) :: LonFind, LatFind
  integer, intent(out) :: iiBlock, iiLon, iiLat
  real, intent(out) :: rLon, rLat

  integer iBlock, iLon, iLat

  iiBlock = -1
  iiLon   = -1
  iiLat   = -1

  do iBlock=1,nBlocks

     if((Longitude(0,iBlock)+Longitude(1,iBlock))/2 <=LonFind .and. &
          (Longitude(nLons,iBlock)+Longitude(nLons+1,iBlock))/2 >LonFind) then

        if((Latitude(0,iBlock)+Latitude(1,iBlock))/2 <=LatFind .and. &
             (Latitude(nLats,iBlock)+Latitude(nLats+1,iBlock))/2 >LatFind) then

           iiBlock = iBlock

           do iLon = 0,nLons
              if(Longitude(iLon,iBlock) <= LonFind .and. &
                   Longitude(iLon+1,iBlock) > LonFind) then
                 iiLon = iLon
                 rLon = 1.0 - (LonFind - Longitude(iLon,iBlock)) / &
                      (Longitude(iLon+1,iBlock) - Longitude(iLon,iBlock))
                 exit
              endif
           enddo

           do iLat = 0,nLats
              if(Latitude(iLat,iBlock) <= LatFind .and. &
                   Latitude(iLat+1,iBlock) > LatFind) then
                 iiLat = iLat
                 rLat = 1.0 - (LatFind - Latitude(iLat,iBlock)) / &
                      (Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock))
                 exit
              endif
           enddo

           if(iiLon >= 0 .and. iiLat >= 0) then
              exit
           end if
        end if
     end if
  end do

end subroutine LocationIndex
