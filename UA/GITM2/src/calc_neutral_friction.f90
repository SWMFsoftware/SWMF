
subroutine calc_neutral_friction(iBlock)

  use ModGITM
  use ModSources

  implicit none

  integer, intent(in) :: iBlock

  ! These are the numerical coefficients in Table 1 in m^2 instead of cm^2
  real, parameter, dimension(4, 4) :: Diff0 = 0.0001 * reshape( (/ &
       ! 0      02     N2      N     NO
       !---------------------------------+
       0.00,  0.260, 0.260, 0.260, &            ! O
       0.26,  0.000, 0.181, 0.181, &            ! O2
       0.26,  0.181, 0.000, 0.181, &            ! N2
       0.26,  0.181, 0.181, 0.000 /), (/4,4/) )  ! N

  ! These are the exponents
  real, parameter, dimension(4, 4) :: DiffExp = reshape( (/ &
       ! 0      02     N2
       !---------------------------------+
       0.00,  0.75,  0.75, 0.75, &             ! O
       0.75,  0.00,  0.75, 0.75, &             ! O2
       0.75,  0.75,  0.00, 0.75, &             ! N2
       0.75,  0.75,  0.75, 0.00 /), (/4,4/) )  ! N

  integer :: iAlt, iLat, iLon

  NeutralFriction = 0.0

  do iAlt = 1, nAlts
     do iLat = 1, nLats
        do iLon = 1, nLons

           call calc_friction

        enddo
     enddo
  enddo

contains

  subroutine calc_friction

    integer :: iSpecies, jSpecies
    real :: CoefMatrix(nSpecies, nSpecies), kTOverM
    real :: Matrix(nSpecies, nSpecies)
    real :: Vel(nSpecies), Parity, Fraction=1.0
    integer :: iPivot(nSpecies)

    Vel = VerticalVelocity(iLon,iLat,iAlt,:,iBlock)

    do iSpecies = 1, nSpecies
       kTOverM = Boltzmanns_Constant * &
            Temperature(iLon, iLat, iAlt,iBlock) * TempUnit / &
            Mass(iSpecies)
       do jSpecies = 1, nSpecies
          if (jSpecies == iSpecies) cycle
          CoefMatrix(iSpecies, jSpecies) = &
               kTOverM * &
               NDensityS(iLon, iLat, iAlt, jSpecies, iBlock) / &
               (Diff0(iSpecies, jSpecies) &
               * (Temperature(iLon, iLat, iAlt,iBlock) &
               / Temperature(iLon, iLat, 1, iBlock)) &
               ** DiffExp(iSpecies, jSpecies) &
               * NDensity(iLon,iLat,1,iBlock))
!          CoefMatrix(iSpecies, jSpecies) = 1.0e10 * kTOverM * &
!               NDensityS(iLon, iLat, iAlt, jSpecies, iBlock)/ &
!               NDensity(iLon, iLat, iAlt, iBlock)
       enddo
    enddo

    Fraction = 1.0e-8 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CoefMatrix = Fraction * CoefMatrix

!    if (iAlt == 1) write(*,*) "low : ",CoefMatrix
!    if (iAlt == 30) write(*,*) "high : ",CoefMatrix
!    if (iAlt == 30) stop

    Matrix = -dt*CoefMatrix

    do iSpecies = 1, nSpecies
       Matrix(iSpecies,iSpecies) = &
            1.0 + dt*sum(CoefMatrix(iSpecies,:))
    enddo

!    if (iAlt == 1) then
!       write(*,*) "iAlt : ",iAlt
!       write(*,*) "Matrix : ", Matrix(1,2), KtOverM, &
!            Diff0(1, 2), &
!            (Temperature(iLon, iLat, iAlt,iBlock) &
!            / Temperature(iLon, iLat, 1, iBlock)) &
!            ** DiffExp(1, 2), &
!            NDensity(iLon,iLat,1,iBlock)
!    endif

!    if (iAlt == 30) write(*,*) "BVel : ", Vel, &
!         sum(mass(1:nSpecies)*vel)/sum(mass(1:nSpecies))

    call ludcmp(Matrix, nSpecies, nSpecies, iPivot, Parity)
    call lubksb(Matrix, nSpecies, nSpecies, iPivot, Vel)

!    Vel = &
!        sum(mass(1:nSpecies)*NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock)*vel)&
!         /sum(mass(1:nSpecies)*NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))

!    if (iAlt == 30) write(*,*) "AVel : ", Vel, &
!         sum(mass(1:nSpecies)*vel)/sum(mass(1:nSpecies))

!    Vel = 0.0

    NeutralFriction(iLon, iLat, iAlt, :) = &
         Vel - VerticalVelocity(iLon,iLat,iAlt,:,iBlock)

!    do iSpecies = 1, nSpecies
!!!!       NeutralFriction(iLon, iLat, iAlt, iSpecies) = &
!!!!            sum(Vel*Matrix(iSpecies,:))
!       if (iAlt ==7) write(*,*) "NF : ", iSpecies, iAlt, &
!            NeutralFriction(iLon, iLat, iAlt, iSpecies), &
!            Matrix(iSpecies,:), &
!            VerticalVelocity(iLon,iLat,iAlt,iSpecies,iBlock), &
!            Vel(iSpecies)
!    enddo

  end subroutine calc_friction

end subroutine calc_neutral_friction
