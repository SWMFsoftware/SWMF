
subroutine chapman_integrals(iBlock)

  use ModGITM
  use ModEUV
  use ModPlanet
  use ModConstants
  implicit none

  integer, intent(in) :: iBlock

  integer :: iLon, iLat, iAlt, iSpecies
  real :: ScaleHeightS

  Chapman = 0.0

  call report("chapman",2)
  call start_timing("chapman_integrals")

  do iLon = 1, nLons

     do iLat = 1, nLats

        if (SZA(iLon,iLat,iBlock) < TwoPi/4 .or. SZA(iLon,iLat,iBlock) > 3*TwoPi/4) then

           do iSpecies = 1, nSpecies

              iAlt = nAlts

              ScaleHeightS = &
                   Temperature(iLon,iLat,iAlt,iBlock)*TempUnit(iLon,iLat,iAlt) * &
                   Boltzmanns_Constant / (-Gravity(iAlt) * Mass(iSpecies))

              Chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
                   NDensityS(iLon,iLat,iAlt,iSpecies, iBlock) * ScaleHeightS

              do iAlt = nAlts-1, 1, -1
                 Chapman(iLon,iLat,iAlt,iSpecies,iBlock) = &
                      Chapman(iLon,iLat,iAlt+1,iSpecies,iBlock) + &
                      NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)*&
                      dAlt(iAlt)*sqrt(1+tan(SZA(iLon,iLat,iBlock))**2)
              enddo

           enddo

        else

           chapman(iLon, iLat, :, :,iBlock) = 1.0e30

        endif

     enddo
  enddo

  call end_timing("chapman_integrals")

end subroutine chapman_integrals
