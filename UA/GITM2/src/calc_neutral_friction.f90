
subroutine calc_neutral_friction(nVel)

  use ModGITM
  use ModSources
  use ModPlanet, only: Diff0, DiffExp
  use ModInputs, only: UseNeutralFriction
  use ModVertical, only: Temp, NDensityS_1D

  implicit none

  real,intent(inout) :: nVel(1:nAlts,1:nSpecies)

  integer :: iSpecies, jSpecies
  real :: CoefMatrix(nSpecies, nSpecies), kTOverM
  real :: Matrix(nSpecies, nSpecies)
  real :: Vel(nSpecies), Parity, Fraction=1.0, n0 = 1.28E+19
  integer :: iPivot(nSpecies)

  integer :: iAlt

  call report("calc_neutral_friction",2)

  if (.not.UseNeutralFriction) return

  do iAlt = 1, nAlts

    Vel = nVel(iAlt,:)

    CoefMatrix = 0.0

    n0 = 1.28e19

    do iSpecies = 1, nSpecies
       kTOverM = Boltzmanns_Constant * Temp(iAlt) / Mass(iSpecies)
       do jSpecies = 1, nSpecies
          if (jSpecies == iSpecies) cycle

          CoefMatrix(iSpecies, jSpecies) = &
               kTOverM * &
               NDensityS_1D(iAlt, jSpecies) / &
               (Diff0(iSpecies, jSpecies) &
               * (Temp(iAlt) / Temp(1)) &
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

    nVel(iAlt, :) = Vel

 enddo

end subroutine calc_neutral_friction
