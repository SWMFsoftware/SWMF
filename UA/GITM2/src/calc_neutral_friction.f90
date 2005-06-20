
subroutine calc_neutral_friction(iBlock)

  use ModGITM
  use ModSources
  use ModInputs, only:f107a

  implicit none

  integer, intent(in) :: iBlock

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  real, parameter, dimension(4, 4) :: Diff0 = 1.0e4 * reshape( (/ &
       ! 0      02     N2      N     NO
       !---------------------------------+
       0.00,  0.260, 0.260, 0.300, &            ! O
       0.26,  0.000, 0.181, 0.220, &            ! O2
       0.26,  0.181, 0.000, 0.220, &            ! N2
       0.30,  0.220, 0.220, 0.000 /), (/4,4/) )  ! N

!  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
!  real, parameter, dimension(4, 4) :: Diff0 = 1.0e4 * reshape( (/ &
!       ! 0      02     N2      N     NO
!       !---------------------------------+
!       0.00,  0.260, 0.260, 0.260, &            ! O
!       0.26,  0.000, 0.181, 0.181, &            ! O2
!       0.26,  0.181, 0.000, 0.181, &            ! N2
!       0.26,  0.181, 0.181, 0.000 /), (/4,4/) )  ! N

  ! These are the exponents
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       ! 0      02     N2
       !---------------------------------+
       0.00,  0.75,  0.75, 0.75, &             ! O
       0.75,  0.00,  0.75, 0.75, &             ! O2
       0.75,  0.75,  0.00, 0.75, &             ! N2
       0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N

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

    Fraction = 1.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CoefMatrix = Fraction * CoefMatrix

!    if (iAlt == 1) write(*,*) "low : ",CoefMatrix
!    if (iAlt == 30) write(*,*) "high : ",CoefMatrix
!    if (iAlt == 30) stop

!    Matrix = -dt*CoefMatrix
!
!    do iSpecies = 1, nSpecies
!       Matrix(iSpecies,iSpecies) = &
!            1.0 + dt*sum(CoefMatrix(iSpecies,:))
!    enddo

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
