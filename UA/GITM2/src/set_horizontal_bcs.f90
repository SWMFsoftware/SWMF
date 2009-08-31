
subroutine set_horizontal_bcs(iBlock)

  use ModSizeGitm
  use ModPlanet, only : nSpecies, nIons
  use ModInputs, only : LatStart, LatEnd
  use ModGITM

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iIon, iSpecies, iLat

  call report("set_horizontal_bcs",2)

  if ( LatStart > -pi/2.0 .and. &
       minval(Latitude(:,iBlock)) < LatStart) then

     do iAlt = -1, nAlts+2
        do iLat = 0,-1,-1

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat+1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat+1,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = Rho(:,iLat+1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= Temperature(:,iLat+1,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat+1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat+1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat+1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat+1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat+1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

  if (LatEnd < pi/2.0 .and. &
       maxval(Latitude(:,iBlocK)) > LatEnd) then

     do iAlt = -1, nAlts+2
        do iLat = nLats+1, nLats+2

           do iSpecies = 1, nSpecies
              VerticalVelocity(:,iLat,iAlt,iSpecies,iBlock) = &
                   VerticalVelocity(:,iLat-1,iAlt,iSpecies,iBlock)
              nDensityS(:,iLat,iAlt,iSpecies,iBlock) = &
                   nDensityS(:,iLat-1,iAlt,iSpecies,iBlock)
           enddo

           Rho(:, iLat,iAlt,iBlock) = Rho(:,iLat-1,iAlt,iBlock)
           Temperature( :, iLat,iAlt,iBlock)= Temperature(:,iLat-1,iAlt,iBlock)
           iTemperature(:, iLat,iAlt,iBlock)=iTemperature(:,iLat-1,iAlt,iBlock)
           eTemperature(:, iLat,iAlt,iBlock)=eTemperature(:,iLat-1,iAlt,iBlock)

           do iIon = 1, nIons
              IDensityS(:,iLat,iAlt,iIon,iBlock) = &
                   IDensityS(:,iLat-1,iAlt,iIon,iBlock)
           enddo

           Velocity(:,iLat,iAlt,iNorth_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iNorth_,iBlock)
           Velocity(:,iLat,iAlt,iEast_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iEast_,iBlock)
           Velocity(:,iLat,iAlt,iUp_,iBlock) = &
                Velocity(:,iLat-1,iAlt,iUp_,iBlock)

           iVelocity(:,iLat,iAlt,iNorth_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iNorth_,iBlock)
           iVelocity(:,iLat,iAlt,iEast_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iEast_,iBlock)
           iVelocity(:,iLat,iAlt,iUp_,iBlock) = &
                iVelocity(:,iLat-1,iAlt,iUp_,iBlock)

        enddo
     enddo
  endif

end subroutine set_horizontal_bcs
