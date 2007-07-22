
subroutine UAM_Gradient(Array_G, Gradient_CD, iBlock)

  ! This routine calculates the gradients of a scalar quantity 
  ! in spherical coordinates. It is working with a single block, but needs the
  ! variable "iBlock" because it has to figure out where it is in the
  ! grid. It is assumed that in Array_G the ghostcells are defined.


  use ModSizeGitm, only : nLons, nLats, nAlts
  use ModGITM, only : Longitude, Latitude, CosLatitude, &
       GradLonP_CB, GradLon0_CB, GradLonM_CB, &
       GradLatP_CB, GradLat0_CB, GradLatM_CB, &
       RadialDistance_GB, InvRadialDistance_GB, iEast_, iNorth_, iUp_

  implicit none

  integer, intent(in) :: iBlock
  real, dimension(-1:nLons+2,-1:nLats+2, -1:nAlts+2), intent(in)  :: Array_G
  real, dimension(nLons,nLats,nAlts, 3), intent(out) :: Gradient_CD

  ! These are delta altitude variables for a nonuniform grid:
  real :: Drp, Drm, DrmOverDrp, Drr2, Bottom


  real :: InvCosLat

  integer :: iLat, iLon, iAlt
  !---------------------------------------------------------------------------
  Gradient_CD = 0.0

  !\
  ! East First
  !/

  if (nLons > 1) then

     do iLat = 1, nLats
        !!! angle is limited to 80 degrees. Why ???
        InvCosLat = 1./max(CosLatitude(iLat,iBlock),0.17)  
        do iAlt = 1, nAlts
           do iLon = 1, nLons
              Gradient_CD(iLon, iLat, iAlt, iEast_) = &
                   (GradLonM_CB(iLon, iBlock)*Array_G(iLon-1,iLat,iAlt) &
                   +GradLon0_CB(iLon, iBlock)*Array_G(iLon  ,iLat,iAlt) &
                   +GradLonP_CB(iLon, iBlock)*Array_G(iLon+1,iLat,iAlt) &
                   )*InvRadialDistance_GB(iLon, iLat, iAlt, iBlock)*InvCosLat
           enddo
        enddo
     enddo

  endif

  !\
  ! North Second
  !/

  if (nLats > 1) then
     do iAlt = 1, nAlts
        do iLat = 1, nLats
           do iLon = 1, nLons
              Gradient_CD(iLon, iLat, iAlt, iNorth_) = &
                   (GradLatM_CB(iLat, iBlock)*Array_G(iLon,iLat-1,iAlt) &
                   +GradLat0_CB(iLat, iBlock)*Array_G(iLon,iLat  ,iAlt) &
                   +GradLatP_CB(iLat, iBlock)*Array_G(iLon,iLat+1,iAlt) &
                   )*InvRadialDistance_GB(iLon, iLat, iAlt, iBlock)
           enddo
        enddo
     enddo

  endif

  !\
  ! Up Third
  !/

  if (nAlts > 1) then

     do iAlt = 1, nAlts
        do iLat = 1, nLats
           do iLon = 1, nLons
              ! Second order gradient for non-uniform grid
              Drp        = RadialDistance_GB(iLon, iLat, iAlt+1, iBlock) &
                   -       RadialDistance_GB(iLon, iLat, iAlt  , iBlock)
              Drm        = RadialDistance_GB(iLon, iLat, iAlt  , iBlock) &
                   -       RadialDistance_GB(iLon, iLat, iAlt-1, iBlock)
              DrmOverDrp = Drm / Drp
              Drr2       = DrmOverDrp**2
              Bottom     = Drm*(1+DrmOverDrp)

              Gradient_CD(iLon,iLat,iAlt,iUp_) = &
                   (     Drr2*Array_G(iLon,iLat,iAlt+1) &
                   -          Array_G(iLon,iLat,iAlt-1) &
                   - (Drr2-1)*Array_G(iLon,iLat,iAlt)) / Bottom
              
           enddo
        enddo
     enddo

  endif

end subroutine UAM_Gradient
