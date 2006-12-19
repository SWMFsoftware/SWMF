
subroutine UAM_Gradient(InArray, OutArray, iBlock)

  ! This routine calculates the gradient of a scalar quantity.
  ! It assumes that it is working with a single block, but needs the
  ! variable "iBlock" because it has to figure out where it is in the
  ! grid.  This can be generalized by feeding in the specific grid
  ! for the given variable.  It is also assumed that the "InArray"
  ! and "OutArray" have the ghostcells defined.


  use ModSizeGitm, only : nLons, nLats, nAlts
  use ModGITM, only : Latitude, Longitude, RadialDistance, &
                     iEast_, iNorth_, iUp_

  implicit none

  integer, intent(in) :: iBlock
  real, dimension(-1:nLons+2,-1:nLats+2, -1:nAlts+2), intent(in)  :: InArray
  real, dimension(nLons,nLats,nAlts, 3), intent(out) :: OutArray

  real :: InArray1D(0:nAlts+1), OutArray1D(nAlts)

  ! These are delta altitude variables for a nonuniform grid:
  real, dimension(nAlts) :: drm, drp, drmodrp, drr2, dr, bottom

  ! These are delta theta (i.e. latitude) variables for a nonuniform grid:
  real, dimension(nLats) :: dtm, dtp, dtmodtp, dtr2, dt

  real :: maxi

  integer :: iLat, iLon, iAlt

  OutArray = 0.0

  if (nAlts > 1) then

     do iAlt=1,nAlts
        drp(iAlt) = RadialDistance(iAlt+1) - RadialDistance(iAlt)
        drm(iAlt) = RadialDistance(iAlt)   - RadialDistance(iAlt-1)
        drmodrp(iAlt) = drm(iAlt)/drp(iAlt)
        drr2(iAlt) = drmodrp(iAlt) * drmodrp(iAlt)
     enddo

  endif

  if (nLats > 1) then

     do iLat = 1, nLats
        dtm(iLat) = Latitude(iLat,iBlock) - Latitude(iLat-1,iBlock)
        dtp(iLat) = Latitude(iLat+1,iBlock) - Latitude(iLat,iBlock)
        dtmodtp(iLat) = dtm(iLat) / dtp(iLat)
        dtr2(iLat) = dtmodtp(iLat) * dtmodtp(iLat)
     enddo

  endif

  !\
  ! East First
  !/

  if (nLons > 1) then

     do iLon = 1, nLons
        do iLat = 1, nLats
           ! this is 80 degrees...
           maxi = max(cos(Latitude(iLat,iBlock)),0.17)  
           do iAlt = 1, nAlts
              OutArray(iLon,iLat,iAlt, iEast_) = &
                   ((InArray(iLon+1,iLat,iAlt) - InArray(iLon-1,iLat,iAlt)) / &
                   (Longitude(iLon+1,iBlock) - Longitude(iLon-1,iBlock)))/ &
                   (maxi*RadialDistance(iAlt))
           enddo
        enddo
     enddo

  endif

  !\
  ! North Second
  !/

  if (nLats > 1) then

     do iLon = 1, nLons
        do iLat = 1, nLats
           do iAlt = 1, nAlts
              OutArray(iLon, iLat, iAlt, iNorth_) = &
                   ((dtr2(iLat)*InArray(iLon,iLat+1,iAlt) - &
                                InArray(iLon,iLat-1,iAlt) - &
                   (dtr2(iLat)-1)*InArray(iLon,iLat,iAlt)) / &
                   (dtr2(iLat)*dtp(iLat) + dtm(iLat))) / &
                   RadialDistance(iAlt)
           enddo
        enddo
     enddo

  endif

  !\
  ! Up Third
  !/

  if (nAlts > 1) then

     bottom = drr2*drp + drm

     do iLon = 1, nLons
        do iLat = 1, nLats
           InArray1d = InArray(iLon,iLat,0:nAlts+1)
           OutArray1D = (drr2*InArray1D(2:nAlts+1) - &
                InArray1D(0:nAlts-1) - &
                (drr2-1)*InArray1D(1:nAlts)) / bottom
           OutArray(iLon,iLat,:,iUp_) = OutArray1D
        enddo
     enddo

  endif

end subroutine UAM_Gradient
