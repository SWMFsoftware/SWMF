
subroutine init_msis(iBlock)

  use ModPlanet
  use ModGITM
  use ModEUV
  use ModInputs, only: iDebugLevel
  use ModSources, only: EHeatingRate

  implicit none

  integer, intent(in) :: iBlock
  real , Dimension(-1:nAlts + 2) :: newalt
  real , Dimension(-1:nAlts + 2) :: Heffalt, HeffPH
  real , Dimension(-1:nAlts + 2) :: eheatalt, eheat
  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InNDensityS 
  real , Dimension(-1:nAlts + 2,nIons - 1) :: InIDensityS

!\
! Below are the Magnetospheric Forcing Terms

  real , Dimension(-1:nAlts + 2,nSpeciesTotal) :: InMagDiss 
  real , Dimension(-1:nAlts + 2,nIons - 1) :: InMagIon

  integer :: iiLon,iiLat,iiAlt
  integer :: iLon,iLat,iAlt, iSpecies, iIon

  write(*,*) '==> Now Initializing Titan Background Composition', iBlock   

!\
! Code for diagnostic purposes only ---------------------------- 
!/
  if(iDebugLevel > 4) then
  write(*,*) '==> nLats ==', nLats   
  write(*,*) '==> nLons ==', nLons   
  write(*,*) '==> nAlts ==', nAlts   
  write(*,*) '==> nSpecies ==', nSpecies   
  write(*,*) '==> nBlocks ==', nBlocks   

!       write(*,*) 'Sequence of Altitudes' 
!    do iAlt = -1,nAlts + 2
!       write(*,*) Altitude(iAlt)
!    enddo

  endif
!\
! Code for diagnostic purposes only ---------------------------- 
!/


!\
! New File with all 12 Neutral Species in it---------------------------------------+
!/
  open(UNIT = 25, FILE = 'DataIn/Comp_Lcv_GITM_25km.txt', STATUS='OLD',ACTION = 'READ')
!  open(UNIT = 25, FILE = 'DataIn/Comp_Lcv_GITM_10km.txt', STATUS='OLD',ACTION = 'READ')

!\
! Old File with only N2, CH4, H2, and HCN in it!------------------------------------+
!  open(UNIT = 25, FILE = 'DataIn/Leb_GITM_25km.txt', STATUS='OLD',ACTION = 'READ')
!  open(UNIT = 25, FILE = 'DataIn/Leb_GITM.txt', STATUS='OLD',ACTION = 'READ')
!  open(UNIT = 25, FILE = 'DataIn/Comp_Lcv_GITM_25km.txt', STATUS='OLD',ACTION = 'READ')
!  open(UNIT = 25, FILE = 'DataIn/Comp_Ta_GITM_25km.txt', STATUS='OLD',ACTION = 'READ')
!  open(UNIT = 25, FILE = 'DataIn/Comp_Ta_GITM_10km.txt', STATUS='OLD',ACTION = 'READ')
!110 FORMAT(F6.1,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3)

  iiLon = 1
  iiLat = 1

111 FORMAT(F6.1,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
              1X,ES10.3)

!111 FORMAT(F6.1,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3)

  do iiAlt = -1,nAlts + 2 
        read(25,111) &
      	newalt(iiAlt), &
       	InTemp(iiLon,iiLat,iiAlt,iBlock), &
!
       	InNDensityS(iiAlt,iH_), &
       	InNDensityS(iiAlt,iH2_), &
       	InNDensityS(iiAlt,iCH_), &
       	InNDensityS(iiAlt,i1CH2_), &
       	InNDensityS(iiAlt,i3CH2_), &
       	InNDensityS(iiAlt,iCH3_), &
       	InNDensityS(iiAlt,iCH4_), &
       	InNDensityS(iiAlt,iC2H4_), &
       	InNDensityS(iiAlt,iN4S_), &
       	InNDensityS(iiAlt,iN2_), &
       	InNDensityS(iiAlt,iHCN_), &
       	InNDensityS(iiAlt,iH2CN_)
!
!       	InIDensityS(iiAlt,iN2P_), &
!       	InIDensityS(iiAlt,iNP_), &
!       	InIDensityS(iiAlt,iCH3P_), &
!       	InIDensityS(iiAlt,iH2CNP_), &
!       	InIDensityS(iiAlt,iC2H5P_)
  end do

  close(Unit = 25)

131 FORMAT(F6.1,1X,ES10.3)
  open(UNIT = 27, FILE = 'DataIn/HCN_25km.txt', STATUS='OLD',ACTION = 'READ')
  do iiAlt = -1,nAlts + 2 
        read(27,131) &
      	newalt(iiAlt), &
       	InNDensityS(iiAlt,iHCN_)
  end do
  close(Unit = 27)


!\
! Below we read in the Eddy Diffusion Profile and the 
! Electron Temperature
!---------------------------------------------------------+

125 FORMAT(1X,F9.4,1X,ES10.3)

!  open(UNIT = 26, FILE = 'DataIn/TeKe_GITM_10km.txt', STATUS='OLD',ACTION = 'READ')
  open(UNIT = 26, FILE = 'DataIn/TeKe_GITM_25km.txt', STATUS='OLD',ACTION = 'READ')

  do iiAlt = -1,nAlts + 2 
        read(26,125) &
        eTemperature(iiLon,iiLat,iiAlt,iBlock),&
        Ke1D(iiAlt) 
  enddo

  close(Unit = 26)

!\
! Below we read in the latest heating efficiency 
!---------------------------------------------------------+

135 FORMAT(F6.1,1X,ES10.3,1X,ES10.3)

     write(*,*) 'Reading New Heating Efficiency File'

  open(UNIT = 35,FILE = 'DataIN/Heff_25km.txt',STATUS='OLD',ACTION='READ')

  do iiAlt = -1,nAlts+2
     read(35,135) &
     Heffalt(iiAlt),&
    Heff3D(iiAlt,iBlock), &
     HeffPH(iiAlt)
  enddo

  close(Unit = 35)

!---------------------------------------------------------+
!
!/



!\
! Below we read in the latest heating efficiency 
!---------------------------------------------------------+

!145 FORMAT(F8.3,1X,F10.5)
!
!     write(*,*) 'Reading New Electron Heating Rates File'
!
!  !open(UNIT = 36,FILE = 'DataIN/GitmEheatp.txt',STATUS='OLD',ACTION='READ')
!  open(UNIT = 36,FILE = 'DataIN/GitmEheatb.txt',STATUS='OLD',ACTION='READ')
!
!  do iiAlt = -1,nAlts+2
!     read(36,145) &
!     eheatalt(iiAlt),&
!     eheat(iiAlt)
!  enddo

!\
! eheat is in units of eV/cm^3/s
! thus we want to multiply by element of charge (1.6022e-19) and multiply by 10^6
! /

!  close(Unit = 36)



!---------------------------------------------------------+
!
!/


115 FORMAT(55(1X,ES15.3))

!  open(UNIT = 27, FILE = 'DataIn/EuvFluxes_10km.txt', STATUS='OLD',ACTION = 'READ')
!
!    do iiAlt = -1,nAlts + 2 
!        read(27,115) &
!            SolarEuv1D(iiAlt,1:55,iBlock)
!    read(27,115) SolarEuv1D(iiAlt,1,iBlock),SolarEuv1D(iiAlt,2,iBlock),SolarEuv1D(iiAlt,3,iBlock),SolarEuv1D(iiAlt,4,iBlock),  &
!                 SolarEuv1D(iiAlt,5,iBlock), &
!                 SolarEuv1D(iiAlt,6,iBlock),SolarEuv1D(iiAlt,7,iBlock),SolarEuv1D(iiAlt,8,iBlock),SolarEuv1D(iiAlt,9,iBlock), &
!                 SolarEuv1D(iiAlt,10,iBlock),&
!                 SolarEuv1D(iiAlt,11,iBlock),SolarEuv1D(iiAlt,12,iBlock),SolarEuv1D(iiAlt,13,iBlock),SolarEuv1D(iiAlt,14,iBlock), &
!                 SolarEuv1D(iiAlt,15,iBlock),&
!                 SolarEuv1D(iiAlt,16,iBlock),SolarEuv1D(iiAlt,17,iBlock),SolarEuv1D(iiAlt,18,iBlock),SolarEuv1D(iiAlt,19,iBlock), &
!                 SolarEuv1D(iiAlt,20,iBlock),&
!                 SolarEuv1D(iiAlt,21,iBlock),SolarEuv1D(iiAlt,22,iBlock),SolarEuv1D(iiAlt,23,iBlock),SolarEuv1D(iiAlt,24,iBlock), &
!                 SolarEuv1D(iiAlt,25,iBlock),&
!                 SolarEuv1D(iiAlt,26,iBlock),SolarEuv1D(iiAlt,27,iBlock),SolarEuv1D(iiAlt,28,iBlock),SolarEuv1D(iiAlt,29,iBlock), &
!                 SolarEuv1D(iiAlt,30,iBlock),&
!                 SolarEuv1D(iiAlt,31,iBlock),SolarEuv1D(iiAlt,32,iBlock),SolarEuv1D(iiAlt,33,iBlock),SolarEuv1D(iiAlt,34,iBlock), &
!                 SolarEuv1D(iiAlt,35,iBlock),&
!                 SolarEuv1D(iiAlt,36,iBlock),SolarEuv1D(iiAlt,37,iBlock),SolarEuv1D(iiAlt,38,iBlock),SolarEuv1D(iiAlt,39,iBlock), &
!                 SolarEuv1D(iiAlt,40,iBlock),&
!                 SolarEuv1D(iiAlt,41,iBlock),SolarEuv1D(iiAlt,42,iBlock),SolarEuv1D(iiAlt,43,iBlock),SolarEuv1D(iiAlt,44,iBlock), &
!                 SolarEuv1D(iiAlt,45,iBlock),&
!                 SolarEuv1D(iiAlt,46,iBlock),SolarEuv1D(iiAlt,47,iBlock),SolarEuv1D(iiAlt,48,iBlock),SolarEuv1D(iiAlt,49,iBlock), &
!                 SolarEuv1D(iiAlt,50,iBlock),&
!                 SolarEuv1D(iiAlt,51,iBlock),SolarEuv1D(iiAlt,52,iBlock),SolarEuv1D(iiAlt,53,iBlock),SolarEuv1D(iiAlt,54,iBlock), &
!                 SolarEuv1D(iiAlt,55,iBlock)
!    end do
!
!  close(Unit = 27)
!
!  open(UNIT = 28, FILE = '1dEuvFluxes.txt', STATUS='NEW',ACTION = 'WRITE')
!
!    do iiAlt = -1,nAlts + 2 
!
!        write(28,115) &
!            SolarEuv1D(iiAlt,1:55,iBlock)
!
!              write(28,115)  &
!                 SolarEuv1D(iiAlt,1,iBlock),SolarEuv1D(iiAlt,2,iBlock),SolarEuv1D(iiAlt,3,iBlock),SolarEuv1D(iiAlt,4,iBlock),  &
!                 SolarEuv1D(iiAlt,5,iBlock), &
!                 SolarEuv1D(iiAlt,6,iBlock),SolarEuv1D(iiAlt,7,iBlock),SolarEuv1D(iiAlt,8,iBlock),SolarEuv1D(iiAlt,9,iBlock), &
!                 SolarEuv1D(iiAlt,10,iBlock),&
!                 SolarEuv1D(iiAlt,11,iBlock),SolarEuv1D(iiAlt,12,iBlock),SolarEuv1D(iiAlt,13,iBlock),SolarEuv1D(iiAlt,14,iBlock), &
!                 SolarEuv1D(iiAlt,15,iBlock),&
!                 SolarEuv1D(iiAlt,16,iBlock),SolarEuv1D(iiAlt,17,iBlock),SolarEuv1D(iiAlt,18,iBlock),SolarEuv1D(iiAlt,19,iBlock), &
!                 SolarEuv1D(iiAlt,20,iBlock),&
!                 SolarEuv1D(iiAlt,21,iBlock),SolarEuv1D(iiAlt,22,iBlock),SolarEuv1D(iiAlt,23,iBlock),SolarEuv1D(iiAlt,24,iBlock), &
!                 SolarEuv1D(iiAlt,25,iBlock),&
!                 SolarEuv1D(iiAlt,26,iBlock),SolarEuv1D(iiAlt,27,iBlock),SolarEuv1D(iiAlt,28,iBlock),SolarEuv1D(iiAlt,29,iBlock), &
!                 SolarEuv1D(iiAlt,30,iBlock),&
!                 SolarEuv1D(iiAlt,31,iBlock),SolarEuv1D(iiAlt,32,iBlock),SolarEuv1D(iiAlt,33,iBlock),SolarEuv1D(iiAlt,34,iBlock), &
!                 SolarEuv1D(iiAlt,35,iBlock),&
!                 SolarEuv1D(iiAlt,36,iBlock),SolarEuv1D(iiAlt,37,iBlock),SolarEuv1D(iiAlt,38,iBlock),SolarEuv1D(iiAlt,39,iBlock), &
!                 SolarEuv1D(iiAlt,40,iBlock),&
!                 SolarEuv1D(iiAlt,41,iBlock),SolarEuv1D(iiAlt,42,iBlock),SolarEuv1D(iiAlt,43,iBlock),SolarEuv1D(iiAlt,44,iBlock), &
!                 SolarEuv1D(iiAlt,45,iBlock),&
!                 SolarEuv1D(iiAlt,46,iBlock),SolarEuv1D(iiAlt,47,iBlock),SolarEuv1D(iiAlt,48,iBlock),SolarEuv1D(iiAlt,49,iBlock), &
!                 SolarEuv1D(iiAlt,50,iBlock),&
!                 SolarEuv1D(iiAlt,51,iBlock),SolarEuv1D(iiAlt,52,iBlock),SolarEuv1D(iiAlt,53,iBlock),SolarEuv1D(iiAlt,54,iBlock), &
!                 SolarEuv1D(iiAlt,55,iBlock)
!    end do
!
!  close(Unit = 28)

                     InMagDiss(-1:nAlts+2, 1:nSpeciesTotal) = 0.0
                     InMagIon(-1:nAlts+2, 1:nIons - 1) = 0.0

     write(*,*) 'Opening Magnetospheric Forcing File'

open(UNIT = 55, FILE = 'DataIn/MagForcing_12SLT.txt', STATUS='OLD',ACTION = 'READ')

155 FORMAT(     F6.1, &
                1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3, &
                1X,ES10.3, 1X,ES10.3, 1X,ES10.3, 1X,ES10.3 )

  do iiAlt = 1, nAlts
       read(55,155)   eheatalt(iiAlt), &
                     InMagDiss(iiAlt,iN4S_), &
                     InMagDiss(iiAlt,iCH_), &
                     InMagDiss(iiAlt,i3CH2_), &
                     InMagDiss(iiAlt,i1CH2_), &
                     InMagDiss(iiAlt,iCH3_), &
                      InMagIon(iiAlt,iNP_), &
                      InMagIon(iiAlt,iN2P_), &
                      InMagIon(iiAlt,iCH3P_), &
                         eheat(iiAlt)
  enddo
 close(Unit = 55)

     write(*,*) 'End Magnetospheric Forcing File'

!     write(*,*) 'Opening Magnetospheric Output File'
!
!open(UNIT = 56, FILE = 'DataIn/MagForcing.out', STATUS='NEW',ACTION = 'WRITE')
!
!  do iiAlt = 1, nAlts
!       write(56,155)  eheatalt(iiAlt), &
!                     InMagDiss(iiAlt,iN4S_), &
!                     InMagDiss(iiAlt,iCH_), &
!                     InMagDiss(iiAlt,i3CH2_), &
!                     InMagDiss(iiAlt,i1CH2_), &
!                     InMagDiss(iiAlt,iCH3_), &
!                      InMagIon(iiAlt,iNP_), &
!                      InMagIon(iiAlt,iN2P_), &
!                      InMagIon(iiAlt,iCH3P_), &
!                         eheat(iiAlt)
!  enddo
!  close(Unit = 56)
!
!     write(*,*) 'End Magnetospheric Output'




!\
! Scale Densities to m^-3 from cm^-3 multiplying by 10^6
!/
!  NDensityS(iiLon,iiLat,-1:nAlts+2,1:nSpeciesTotal,iBlock) = &
!  NDensityS(iiLon,iiLat,-1:nAlts+2,1:nSpeciesTotal,iBlock) * &
!  1.0E06


!\
! Set the Fixed background
!/

    Ke1D(:) = Ke1D(:)*(1.0e-04)

!\
! Initializes the Planet with the same Chemistry as Above 
!/

    do iLat = -1, nLats + 2
      do iLon = -1, nLons + 2

       Temperature(iLon,iLat,:,iBlock) =  &
     Temperature(iiLat,iiLon,:,iBlock)

       eTemperature(iLon,iLat,:,iBlock) =  &
     eTemperature(iiLat,iiLon,:,iBlock)

       NDensityS(iLon,iLat,-1:nAlts + 2,1:nSpeciesTotal,iBlock) =  0.0

              do iSpecies = 1, nSpeciesTotal

       NDensityS(iLon,iLat,-1:nAlts + 2,iSpecies,iBlock) =  &
             InNDensityS(-1:nAlts + 2,iSpecies)*(1.0e+06)       !Converts from cm^-3 to m^-3

             enddo ! end inner ispecies loop

                IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = 0.0 

              do iIon = 1, nIons-1


                IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) =  1.0


                IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) = &
                   IDensityS(iLon,iLat,-1:nAlts + 2,ie_,iBlock) + &
                   IDensityS(iLon,iLat,-1:nAlts + 2,iIon,iBlock) 

             enddo ! end inner ispecies loop




       enddo! end iLon loop
     enddo ! end iLat loop

!\
! These next two lines allow me to set HCN to a factor of its input value
! without causing any bizarre problems from before
!
    FixHCN(:,:,-1:nAlts+2,iBlock) = NDensityS(:,:,-1:nAlts+2,iHCN_,iBlock)
    NDensityS(:,:,-1:nAlts+2,iHCN_,iBlock) = (1.00)*FixHCN(:,:,-1:nAlts+2,iBlock) 

    FixN2(:,:,-1:nAlts+2,iBlock) = NDensityS(:,:,-1:nAlts+2,iN2_,iBlock)
    FixCH4(:,:,-1:nAlts+2,iBlock) = NDensityS(:,:,-1:nAlts+2,iCH4_,iBlock)
    FixH2(:,:,-1:nAlts+2,iBlock) = NDensityS(:,:,-1:nAlts+2,iH2_,iBlock)
!/
!


!\
! Diagnostic Outputs
!/
!112 FORMAT(F6.1,1X,F8.4,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3,1X,ES10.3, &
!              1X,ES10.3,1X,ES10.3)
!
!  open(UNIT = 26, FILE = 'gimpneutrals.txt', STATUS='NEW',ACTION = 'WRITE')
!
!  do iiAlt = -1,nAlts + 2 
!        write(26,112) &
!      	newalt(iiAlt), &
!       	Temperature(iiLon,iiLat,iiAlt,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,i1CH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,i3CH2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH3_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iCH4_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iC2H4_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iN4S_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iN2_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iHCN_,iBlock), &
!       	NDensityS(iiLon,iiLat,iiAlt,iH2CN_,iBlock), &
!        FixHCN(iiLon,iiLat,iiAlt,iBlock) 
!  end do
!
!  close(Unit = 26)




!\
! Calculating MeanMajorMass -----------------------------
!/

!\
! Initialize MeanMajorMass to 0.0
!/
  MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = 0.0
  MeanIonMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = 0.0


! Calculate MeanMajorMass -----------------------------
! Calculate TempUnit -----------------------------

  MagDissRateS(1:nLons,1:nLats,1:nAlts,1:nSpeciesTotal,iBlock) = 0.0
  MagIonRateS(1:nLons,1:nLats,1:nAlts,1:nIons - 1,iBlock) = 0.0

    do iLat = 1,nLats 
     do iLon = 1,nLons 
      do iAlt = 1,nAlts 

       EHeatingRate(iLon,iLat,iAlt,iBlock) = &
           eheat(iAlt)*(0.05)*Heff3D(iAlt,iBlock)

         do iSpecies = 1,nSpeciesTotal
           MagDissRateS(iLon,iLat,iAlt,iSpecies,iBlock) = &
                     InMagDiss(iAlt,iSpecies)*(0.05)
         enddo

         do iIon = 1,nIons - 1
           MagIonRateS(iLon,iLat,iAlt,iIon, iBlock) = &
                     InMagIon(iAlt,iIon)*(0.05)
         enddo

      enddo
     enddo
    enddo

    do iLat = -1,nLats + 2
     do iLon = -1,nLons + 2
      do iAlt = -1,nAlts + 2

! should convert eheat to J/m^3/s

!       EHeatingRate(iLon,iLat,iAlt,iBlock) = 0.0 

!       EHeatingRate(iLon,iLat,iAlt,iBlock) = &
!           eheat(iAlt)*(0.05)*Heff3D(iAlt,iBlock)

           NDensity(iLon,iLat,iAlt,iBlock) = 0.0

       do iSpecies = 1,nSpeciesTotal
           NDensity(iLon,iLat,iAlt,iBlock) = &
              NDensity(iLon,iLat,iAlt,iBlock) + &
              NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)
       enddo

       do iSpecies = 1,nSpeciesTotal
           MeanMajorMass(iLon,iLat,iAlt) = &
                 MeanMajorMass(iLon,iLat,iAlt) + &
           Mass(iSpecies)*NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)/ &
           NDensity(iLon,iLat,iAlt,iBlock) 
       enddo

       do iIon = 1,nIons - 1
           MeanIonMass(iLon,iLat,iAlt) = &
                 MeanIonMass(iLon,iLat,iAlt) + &
           MassI(iIon)*IDensityS(iLon,iLat,iAlt,iIon,iBlock)/ &
           IDensityS(iLon,iLat,iAlt,ie_,iBlock) 
       enddo


       enddo
      enddo
     enddo

!\
! This sets Ar at the bottom of the model
!
    do iLon = -1,nLons+2
     do iLat = -1,nLats+2

       NDensityS(iLon,iLat,-1:1,iAr_,iBlock) =  &
         NDensity(iLon,iLat,-1:1,iBlock)*(4.32e-05)

       NDensityS(iLon,iLat,2:nAlts+2,iAr_,iBlock) = 1.0e3

      enddo
     enddo


  TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2) = &
             MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)/&
             Boltzmanns_Constant

  if(iDebugLevel > 4) then
  write(*,*) '=====> Values for INIT_MSIS.Titan.f90'
  write(*,*) '=====> Before TEMP UNIT'
  write(*,*) '==> Maxval for Temperature is:', &
              MAXVAL(Temperature(:,:,:,iBlock))   
  write(*,*) '==> Minval for Temperature is:', &
              MINVAL(Temperature(:,:,:,iBlock))   
  endif

!\
! Initialize Rho to 0.0
!/
  Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = 0.0

  Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
  Temperature(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) / &
     TempUnit(-1:nLons+2,-1:nLats+2,-1:nAlts+2)

  Rho(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock) = &
      MeanMajorMass(-1:nLons+2,-1:nLats+2,-1:nAlts+2)* &
      NDensity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iBlock)

  if(iDebugLevel > 4) then
  write(*,*) '=====> Values for INIT_MSIS.Titan.f90'
  write(*,*) '=====> After TEMP UNIT' 
  write(*,*) '==> Maxval for Temperature is:', &
              MAXVAL(Temperature(:,:,:,iBlock))   
  write(*,*) '==> Minval for Temperature is:', &
              MINVAL(Temperature(:,:,:,iBlock))   
  write(*,*) '==> Maxval for Rho is:', &
              MAXVAL(Rho(:,:,:,iBlock))   
  write(*,*) '==> Minval for Rho is:', &
              MINVAL(Rho(:,:,:,iBlock))   
  endif

  if(iDebugLevel > 4) then
  write(*,*) '=====> Checking At the END of INIT_MSIS.Titan.f90'
  write(*,*) '=====> After TEMP UNIT' 

  write(*,*) '==> Maxval for NHCN is:', &
              MAXVAL(NDensityS(:,:,:,iHCN_,iBlock))   
  write(*,*) '==> Minval for NHCN is:', &
              MINVAL(NDensityS(:,:,:,iHCN_,iBlock))   
  write(*,*) '==> Maxval for FIXHCN is:', &
              MAXVAL(FixHCN(:,:,:,iBlock))   
  write(*,*) '==> Minval for FIXHCN is:', &
              MINVAL(FixHCN(:,:,:,iBlock))   

  write(*,*) '==> Maxval for N2 is:', &
              MAXVAL(NDensityS(:,:,:,iN2_,iBlock))   
  write(*,*) '==> Minval for N2 is:', &
              MINVAL(NDensityS(:,:,:,iN2_,iBlock))   
  endif

  if(iDebugLevel > 4) then
  write(*,*) 'Ion Density Check'
  write(*,*) '==> Maxval for IDensityS is:', &
              MAXVAL(IDensityS(:,:,:,:,iBlock))   
  write(*,*) '==> Minval for IDensityS is:', &
              MINVAL(IDensityS(:,:,:,:,iBlock))   

  write(*,*) '==> Maxval for IDensityS(iN2P) is:', &
              MAXVAL(IDensityS(:,:,:,iN2P_,iBlock))   
  write(*,*) '==> Minval for IDensityS(iN2P_) is:', &
              MINVAL(IDensityS(:,:,:,iN2P_,iBlock))   

  write(*,*) '==> Maxval for IDensityS(iNP) is:', &
              MAXVAL(IDensityS(:,:,:,iNP_,iBlock))   
  write(*,*) '==> Minval for IDensityS(iNP_) is:', &
              MINVAL(IDensityS(:,:,:,iNP_,iBlock))   

  write(*,*) '==> Maxval for IDensityS(iH2CNP) is:', &
              MAXVAL(IDensityS(:,:,:,iH2CNP_,iBlock))   
  write(*,*) '==> Minval for IDensityS(H2CNP_) is:', &
              MINVAL(IDensityS(:,:,:,iH2CNP_,iBlock))   

  write(*,*) '==> Maxval for IDensityS(iC2H5P_) is:', &
              MAXVAL(IDensityS(:,:,:,iC2H5P_,iBlock))   
  write(*,*) '==> Minval for IDensityS(C2H5P_) is:', &
              MINVAL(IDensityS(:,:,:,iC2H5P_,iBlock))   

  write(*,*) '==> Maxval for IDensityS(iCH3P_) is:', &
              MAXVAL(IDensityS(:,:,:,iCH3P_,iBlock))   
  write(*,*) '==> Minval for IDensityS(iCH3P_) is:', &
              MINVAL(IDensityS(:,:,:,iCH3P_,iBlock))   

  write(*,*) '==> Maxval for IDensityS(ie_) is:', &
              MAXVAL(IDensityS(:,:,:,ie_,iBlock))   
  write(*,*) '==> Minval for IDensityS(ie_) is:', &
              MINVAL(IDensityS(:,:,:,ie_,iBlock))   


  write(*,*) '==> Maxval for MeanIonMass is:', &
              MAXVAL(MeanIonMass(:,:,:))   
  write(*,*) '==> Minval for MeanIonMass is:', &
              MINVAL(MeanIonMass(:,:,:))   
  endif

end subroutine init_msis

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS, Temp, LogRho)

  write(*,*) "You can not use MSIS with any planet except Earth!!!"
  write(*,*) "If you ARE running Earth, then make the code again, using"
  write(*,*) "configure Earth ; make"
  call stop_gitm("I can not continue...")

end subroutine msis_bcs
