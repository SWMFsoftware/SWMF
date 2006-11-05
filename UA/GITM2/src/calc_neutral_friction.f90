
subroutine calc_neutral_friction(iBlock)

  use ModGITM
  use ModSources
  use ModPlanet, only: Diff0, DiffExp
  use ModInputs, only:f107a

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon
  real :: RealTemp(-1:nAlts+2), de(4,4)

  call report("calc_neutral_friction",2)

  NeutralFriction = 0.0

  do iLon = 1, nLons
     do iLat = 1, nLats

        RealTemp = Temperature(iLon, iLat, :, iBlock) * &
             TempUnit(iLon, iLat, :)

        do iAlt = 1, nAlts

           call calc_friction

        enddo
     enddo
  enddo

contains

  subroutine calc_friction

    integer :: iSpecies, jSpecies
    real :: CoefMatrix(nSpecies, nSpecies), kTOverM
    real :: Matrix(nSpecies, nSpecies)
    real :: Vel(nSpecies), Parity, Fraction=1.0, n0 = 1.28E+19
    integer :: iPivot(nSpecies)

    Vel = VerticalVelocity(iLon,iLat,iAlt,:,iBlock)

    CoefMatrix = 0.0

    n0 = 1.28e19

    do iSpecies = 1, nSpecies
       kTOverM = Boltzmanns_Constant * RealTemp(iAlt) / Mass(iSpecies)
       do jSpecies = 1, nSpecies
          if (jSpecies == iSpecies) cycle

          CoefMatrix(iSpecies, jSpecies) = &
               kTOverM * &
               NDensityS(iLon, iLat, iAlt, jSpecies, iBlock) / &
               (Diff0(iSpecies, jSpecies) &
               * (RealTemp(iAlt) / RealTemp(1)) &
               ** DiffExp(iSpecies, jSpecies) &
               * N0)

       enddo
    enddo

    Fraction = 1.4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CoefMatrix = Fraction * CoefMatrix

    Matrix = -dt*CoefMatrix

    do iSpecies = 1, nSpecies
       Matrix(iSpecies,iSpecies) = &
            1.0 + dt*(sum(CoefMatrix(iSpecies,:)))
    enddo

    call ludcmp(Matrix, nSpecies, nSpecies, iPivot, Parity)
    call lubksb(Matrix, nSpecies, nSpecies, iPivot, Vel)

    NeutralFriction(iLon, iLat, iAlt, :) = &
         Vel - VerticalVelocity(iLon,iLat,iAlt,:,iBlock)

  end subroutine calc_friction

end subroutine calc_neutral_friction
