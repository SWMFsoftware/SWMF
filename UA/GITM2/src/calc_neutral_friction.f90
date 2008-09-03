
subroutine calc_neutral_friction(oVel, EddyCoef_1d, NDensity_1d, NDensityS_1d, &
                                 EddyCoefRatio_1d, Temp, Gravity_1d )

  use ModGITM
  use ModSources
  use ModPlanet, only: Diff0, DiffExp, IsEarth
  use ModInputs, only: UseNeutralFriction

  implicit none

  real,intent(inout) :: oVel(1:nAlts,1:nSpecies)
  real,intent(in) :: EddyCoef_1d(1:nAlts)
  real,intent(in) :: NDensity_1d(1:nAlts)
  real,intent(in) :: NDensityS_1d(1:nAlts,1:nSpecies)
  real,intent(inout) :: EddyCoefRatio_1d(1:nAlts,1:nSpecies)
  real,intent(in) :: Temp(1:nAlts)
  real,intent(in) :: Gravity_1d(1:nAlts)

  integer :: iSpecies, jSpecies
  real :: CoefMatrix(nSpecies, nSpecies), kTOverM
  real :: Matrix(nSpecies, nSpecies)
  real :: Vel(nSpecies), Parity, Fraction=1.0, N0 = 1.28E+19
  integer :: iPivot(nSpecies)

! Added by Jared 11-29-2007
  real :: TempDij
  real :: InvDij(nSpecies)
  real :: Dij(nSpecies)
  real :: denscale
  real :: mscale
  real :: mms
  real :: mmwos(nSpecies)

  integer :: iAlt

  call report("calc_neutral_friction",4)

  if (.not.UseNeutralFriction) return

    EddyCoefRatio_1d(1:nAlts,1:nSpecies) = 0.0

  do iAlt = 1, nAlts

    Vel = oVel(iAlt,1:nSpecies)
    CoefMatrix = 0.0

    mms = 0.0
    mmwos = 0.0
    InvDij = 0.0


! \
! Below, mms = mean molecular weight
! at the Altitude(iAlt)
!
! In contrast, mmwos = Mean Molecular Weight
! excluding the species (iSpecies)

    do iSpecies = 1, nSpecies
         mms = mms + &
           Mass(iSpecies)*NDensityS_1d(iAlt,iSpecies)/NDensity_1d(iAlt)

         do jSpecies = 1, nSpecies

          if(jSpecies == iSpecies) cycle

         mmwos(iSpecies) = mmwos(iSpecies) +  & 
          Mass(jSpecies)*NDensityS_1d(iAlt,jSpecies)/&
       (NDensity_1d(iAlt) - NDensityS_1d(iAlt,iSpecies) )

            
         enddo ! jSpecies 

    enddo ! iSpecies 



    do iSpecies = 1, nSpecies


       InvDij(iSpecies) = 0.0


       kTOverM = Boltzmanns_Constant * Temp(iAlt) / Mass(iSpecies)


       mscale = (1.0 - Mass(iSpecies)/mmwos(iSpecies)) * &
        (NDensityS_1d(iAlt,iSpecies)/NDensity_1d(iAlt) )

       denscale = 1.0/ &
        (NDensity_1d(iAlt) - NDensityS_1d(iAlt,iSpecies))


       do jSpecies = 1, nSpecies
          if (jSpecies == iSpecies) cycle


! \
! Please note that TempDij is the Dij binary coefficients
! Based upon the formulation by Banks and Kokarts.
! These coefficients demand that 
! (1) NDensity be in cm^-3 (hence the 1.0e-06) factor below
! (2) Additionally, the Dij's are in cm^2/s, thus the 1.0e-04 factor

!         TempDij = (1.0e-04)*&              ! Scales the Dij from cm^2/s -> m^2/s
!           (   Diff0(iSpecies,jSpecies)*( Temp(iAlt)**DiffExp(iSpecies,jSpecies) )   ) / &
!           (   NDensity_1d(iAlt)*(1.0e-06) )     ! Converts to #/cm^-3


! Old CoefMatrix

          if (IsEarth) then 

! TempDij is the old Dij coefficient according to Ridley [2006]

             TempDij = Diff0(iSpecies, jSpecies)*&
                      (N0/NDensity_1d(iAlt))*&
                      (Temp(iAlt) / Temp(1))** DiffExp(iSpecies, jSpecies) 

             CoefMatrix(iSpecies, jSpecies) = &
                  kTOverM * NDensityS_1d(iAlt, jSpecies) / &
                   ( (TempDij/(1.0 - mscale)) * (NDensity_1d(iAlt) - NDensityS_1d(iAlt,iSpecies) ) + &
                         NDensity_1d(iAlt)*EddyCoef_1d(iAlt)  )

             InvDij(iSpecies) = InvDij(iSpecies) + &
                               NDensityS_1d(iAlt, jSpecies) / &
                              ( ( TempDij/(1.0 - mscale)) * (NDensity_1d(iAlt) - NDensityS_1d(iAlt,iSpecies) ) + &
                                    NDensity_1d(iAlt)*EddyCoef_1d(iAlt)  )

! Now we add the turbulent collision frequencies
! to the molecular collision frequencies
!
! V_{turb} = (kT/m)*[N_s/( N * Kzz) ]  From Boqueho and Blelly (2005)


!             CoefMatrix(iSpecies, jSpecies) = &
!                  kTOverM * (NDensityS_1d(iAlt, jSpecies) / NDensity_1d(iAlt)) *&
!                   1.0/( TempDij  + EddyCoef_1d(iAlt)  )

!             CoefMatrix(iSpecies, jSpecies) = &
!                  kTOverM * (NDensityS_1d(iAlt, jSpecies) / (NDensity_1d(iAlt) - NDensityS_1d(iAlt,iSpecies) ) ) *&
!                   1.0/( TempDij  + EddyCoef_1d(iAlt)  )


! InvDij is the inverse of the total "diffusion coefficient" or the 
! total collision frequency (whatever you prefer)

           else
              TempDij = (1.0e-04)*&              ! Scales the Dij from cm^2/s -> m^2/s
                (   Diff0(iSpecies,jSpecies)*( Temp(iAlt)**DiffExp(iSpecies,jSpecies) )   ) / &
                (   NDensity_1d(iAlt)*(1.0e-06) )     ! Converts to #/cm^-3

            !   CoefMatrix(iSpecies, jSpecies) = &
            !        kTOverM *  &
            !        NDensityS_1d(iAlt, jSpecies) / &
            !        TempDij

               CoefMatrix(iSpecies, jSpecies) = &
                    kTOverM * denscale * &
                    NDensityS_1d(iAlt, jSpecies) / &
                    ( TempDij/(1.0 - mscale) + EddyCoef_1d(iAlt))

            !   InvDij(iSpecies) = InvDij(iSpecies) + &
            !      NDensityS_1d(iAlt, jSpecies)/TempDij

               InvDij(iSpecies) = InvDij(iSpecies) + &
                  denscale*NDensityS_1d(iAlt, jSpecies)/ &
                 ( TempDij/(1.0 - mscale) + EddyCoef_1d(iAlt) )

           endif




       enddo

         EddyCoefRatio_1d(iAlt,iSpecies) =  &
             -1.0*Dt*( EddyCoef_1d(iAlt)*InvDij(iSpecies) )*(1.0 - mms/Mass(iSpecies))*Gravity_1d(iAlt)

         


    enddo

!    Vel(1:nSpecies) = Vel(1:nSpecies) + EddyCoefRatio_1d(iAlt,iSpecies)

    Fraction = 1.4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CoefMatrix = Fraction * CoefMatrix

    Matrix = -dt*CoefMatrix

    do iSpecies = 1, nSpecies
       Matrix(iSpecies,iSpecies) = &
            1.0 + dt*(sum(CoefMatrix(iSpecies,:)))
    enddo

    call ludcmp(Matrix, nSpecies, nSpecies, iPivot, Parity)
    call lubksb(Matrix, nSpecies, nSpecies, iPivot, Vel)

    oVel(iAlt, 1:nSpecies) = Vel(1:nSpecies) + EddyCoefRatio_1d(iAlt,1:nSpecies) 

 enddo

end subroutine calc_neutral_friction
