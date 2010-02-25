subroutine UA_fill_electrodynamics(UAr2_fac, UAr2_ped, UAr2_hal, UAr2_lats, UAr2_mlts)
  use ModElectrodynamics

  implicit none

  real, dimension(nMagLons+1,nMagLats), intent(out) :: &
       UAr2_fac, UAr2_ped, UAr2_hal, UAr2_lats, UAr2_mlts
  !----------

  UAr2_Fac  = DivJuAltMC
  UAr2_Ped  = SigmaPedersenMC
  UAr2_Hal  = SigmaHallMC
  UAr2_lats = MagLatMC
  UAr2_mlts = MagLocTimeMC

end subroutine UA_fill_electrodynamics

!\

subroutine UA_calc_electrodynamics(UAi_nMLTs, UAi_nLats)
  use ModGITM
  use ModInputs, only: iDebugLevel, iStartTime, UseApex, UseDynamo, Is1D
  use ModConstants
  use ModElectrodynamics
  use ModLinearSolver
  use ModMPI
  use ModTime
  use ModMagTrace

  implicit none

  integer, intent(out) :: UAi_nMLTs, UAi_nLats

  integer, external :: jday

  integer :: i,j,k,bs, iError, iBlock, iDir, iLon, iLat, iAlt, ip, im

  real :: zzV
  real :: GeoLat, GeoLon, GeoAlt, xAlt, len
  real :: xmag, ymag, zmag, bmag, signz, magpot, lShell
  real :: mlatMC, mltMC, jul, shl, spl, length, kdpm_s, kdlm_s
  real :: kpm_s, klm_s, xstretch, ystretch
  real :: sinIm, spp, sll, scc, be3

  real :: q2, dl, cD

  real :: magloctime_local(0:nLons+1,0:nLats+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       JuDotB, ue1, ue2, ed1, ed2, sigma0

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocksMax) :: &
       sigmap_d1d1_d, sigmap_d2d2_d, sigmap_d1d2_d, sigmah, sigmap, DivJu, kmp, kml, je1, je2

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: Gradient_GC

  real :: aLat, aLon, gLat, gLon, Date, sLat, sLon, gLatMC, gLonMC

  real :: residual, oldresidual, a, tmp

  logical :: IsDone, IsFirstTime = .true., DoTestMe

  integer :: iLm, iLp, jLat, iI, MaxIteration, nIteration
  integer :: nX

  external :: matvec_gitm
  !----------

  if (DipoleStrength == 0) return

  if (.not. UseDynamo .or. Is1D) return

  call report("UA_calc_electrodynamics",1)
  call start_timing("calc_electrodyn")

  if (IsFirstTime) then
     IsFirstTime = .false.

     allocate(DivJuAltMC(nMagLons+1,nMagLats), &
          SigmaHallMC(nMagLons+1,nMagLats), &
          SigmaPedersenMC(nMagLons+1,nMagLats), &
          SigmaLLMC(nMagLons+1,nMagLats), &
          SigmaPPMC(nMagLons+1,nMagLats), &
          SigmaHHMC(nMagLons+1,nMagLats), &
          SigmaCCMC(nMagLons+1,nMagLats), &
          SigmaLPMC(nMagLons+1,nMagLats), &
          SigmaPLMC(nMagLons+1,nMagLats), &
          KpmMC(nMagLons+1,nMagLats), &
          KlmMC(nMagLons+1,nMagLats), &
          KDpmMC(nMagLons+1,nMagLats), &
          KDlmMC(nMagLons+1,nMagLats), &
          MagBufferMC(nMagLons+1,nMagLats), &
          LengthMC(nMagLons+1,nMagLats), &
          GeoLatMC(nMagLons+1,nMagLats), &
          GeoLonMC(nMagLons+1,nMagLats), &
          MagLocTimeMC(nMagLons+1,nMagLats), &
          MagLonMC(nMagLons+1,nMagLats), &
          MagLatMC(nMagLons+1,nMagLats), &
          solver_a_mc(nMagLons+1,nMagLats), &
          solver_b_mc(nMagLons+1,nMagLats), &
          solver_c_mc(nMagLons+1,nMagLats), &
          solver_d_mc(nMagLons+1,nMagLats), &
          solver_e_mc(nMagLons+1,nMagLats), &
          solver_s_mc(nMagLons+1,nMagLats), &
          deltalmc(nMagLons+1,nMagLats), &
          deltapmc(nMagLons+1,nMagLats),  &
          dSigmaLLdlMC(nMagLons+1,nMagLats),  &
          dSigmaLPdlMC(nMagLons+1,nMagLats),  &
          dSigmaPLdpMC(nMagLons+1,nMagLats),  &
          dSigmaPPdpMC(nMagLons+1,nMagLats), &
          dKDpmdpMC(nMagLons+1,nMagLats), &
          dKDlmdlMC(nMagLons+1,nMagLats), &
          dKpmdpMC(nMagLons+1,nMagLats), &
          dKlmdlMC(nMagLons+1,nMagLats), &
          DynamoPotentialMC(nMagLons+1,nMagLats), &
          OldPotMC(nMagLons+1,nMagLats), &
          stat = iError)

     DivJuAltMC = 0.0
     SigmaHallMC = 0.0
     SigmaPedersenMC = 0.0
     SigmaLLMC = 0.0
     SigmaPPMC = 0.0
     SigmaHHMC = 0.0
     SigmaCCMC = 0.0
     SigmaLPMC = 0.0
     SigmaPLMC = 0.0
     KDpmMC = 0.0
     KDlmMC = 0.0
     KpmMC = 0.0
     KlmMC = 0.0
     MagBufferMC = 0.0
     LengthMC = 0.0
     GeoLatMC = 0.0
     GeoLonMC = 0.0
     MagLocTimeMC = 0.0
     MagLonMC = 0.0
     MagLatMC = 0.0
     solver_a_mc = 0.0
     solver_b_mc = 0.0
     solver_c_mc = 0.0
     solver_d_mc = 0.0
     solver_e_mc = 0.0
     solver_s_mc = 0.0
     deltalmc = 0.0
     deltapmc = 0.0
     dSigmaLLdlMC = 0.0
     dSigmaLPdlMC = 0.0
     dSigmaPLdpMC = 0.0
     dSigmaPPdpMC = 0.0
     dKDpmdpMC = 0.0
     dKDlmdlMC = 0.0
     dKpmdpMC = 0.0
     dKlmdlMC = 0.0
     DynamoPotentialMC = 0.0
     OldPotMC = 0.0

     if (iError /= 0) call CON_stop("Error allocating array DivJuAltMC")

     date = iStartTime(1) + float(iJulianDay)/float(jday(iStartTime(1),12,31))

     do i=1,nMagLons+1
        if (iDebugLevel > 1) &
             write(*,*) "==> Calculating Apex->Geo", i, nMagLons+1
        do j=1,nMagLats

           MagLatMC(i,j) = 180.0 * (float(j)-0.5) / float(nMagLats) - 90.0
           MagLonMC(i,j) = 360.0 *  float(i-1)    / float(nMagLons)

           ! we could also use the stretch_grid function....
!MAGSTRETCH       xstretch =(float(j)-0.5)/float(nMagLats/2) - 1.0
!MAGSTRETCH       call stretch_grid(xstretch,ystretch)
!MAGSTRETCH       MagLatMC(i,j) = ystretch * 90.0

           if (UseApex) then 

              aLat = MagLatMC(i,j)
              aLon = MagLonMC(i,j)

              if (i > 1) then
                 sLat = GeoLatMC(i-1,j)
                 sLon = GeoLonMC(i-1,j)
              elseif (j > 1) then
                 sLat = GeoLatMC(i,j-1)
                 sLon = GeoLonMC(i,j-1)
              else
                 sLat = -100.0
              endif

              call apex_to_geo(date, aLat, aLon, AltMinIono, gLat, gLon, sLat, sLon)

              GeoLatMC(i,j) = gLat*pi/180.0
              GeoLonMC(i,j) = gLon*pi/180.0

           else

              GeoLatMC(i,j) = MagLatMC(i,j)*pi/180.0
              GeoLonMC(i,j) = MagLonMC(i,j)*pi/180.0

           endif

        enddo
     enddo

     do i=1,nMagLons+1
        do j=2,nMagLats-1
           deltalmc(i,j) = MagLatMC(i,j+1) - MagLatMC(i,j-1) 
        enddo
        deltalmc(i,1) = 2*(MagLatMC(i,2) - MagLatMC(i,1))
        deltalmc(i,nMagLats) = 2*(MagLatMC(i,nMagLats) - MagLatMC(i,nMagLats-1))
     enddo

     do j=1,nMagLats
        do i=2,nMagLons
           deltapmc(i,j) = MagLonMC(i+1,j) - MagLonMC(i-1,j) 
        enddo
        deltapmc(1,j) = 2*(MagLonMC(2,j) - MagLonMC(1,j))
        deltapmc(nMagLons+1,j) = deltapmc(1,j)
     enddo

     deltapmc = deltapmc * Pi / 180.0 / 2.0
     deltalmc = deltalmc * Pi / 180.0 / 2.0

  endif

  if(UseApex .and. IsEarth) then
     do i=1,nMagLons+1
        do j=1,nMagLats
           call magloctm(MagLonMC(i,j), SubsolarLatitude, SubsolarLongitude, &
                MagneticPoleColat, MagneticPoleLon, mltMC)
           MagLocTimeMC(i,j) = mod(mltMC+24.0,24.0)
        enddo
     enddo
  end if

  if (.not. UseDynamo .or. Is1D) return

  !\
  ! Magnetic grid is defined as:
  ! MLT is in hours = 0 - 24
  ! Lat is in degrees = -90 - 90
  !/

  DivJuAltMC      = -1.0e32
  SigmaHallMC     = 0.0
  SigmaPedersenMC = 0.0
  LengthMC        = -1.0e32
  KDlmMC = -1.0e32
  KDpmMC = -1.0e32
  KlmMC = -1.0e32
  KpmMC = -1.0e32
  SigmaLPMC = -1.0e32
  SigmaPLMC = -1.0e32
  DivJuAltMC = -1.0e32

  UAi_nLats = nMagLats
  UAi_nMlts = nMagLons+1

  q2 = Element_Charge * Element_Charge

  HallFieldLine     = 0.
  PedersenFieldLine = 0.
  DivJuFieldLine    = 0.
  LengthFieldLine   = 0.
  SigmaPP = 0.
  SigmaLL = 0.
  SigmaCC = 0.
  KDpm = 0.
  KDlm = 0.
  Kpm  = 0.
  Klm  = 0.
  DivJu    = 0.
  DivJuAlt = 0.

  do iBlock = 1, nBlocks

     call calc_physics(iBlock)
     call calc_rates(iBlock)
     call calc_collisions(iBlock)
     call get_potential(iBlock)
     call calc_efield(iBlock)

     e_density = IDensityS(:,:,:,ie_,iBlock)

     Vi = Collisions(:,:,:,iVIN_)
     Ve = Collisions(:,:,:,iVEN_) ! + Collisions(:,:,:,iVEI_)

     MeVen = Mass_Electron * Collisions(:,:,:,iVEN_)
     MeVei = Mass_Electron * Collisions(:,:,:,iVEI_)
     MiVin = MeanIonMass * Collisions(:,:,:,iVIN_)

     VeOe = Ve**2 + e_gyro**2
     ViOi = Vi**2 + i_gyro**2

     sigma0 = q2 * E_Density / (1.0/MeVen + 1.0/MiVin)

     sigmap(:,:,:,iBlock) = ((1.0/MeVen) * (Ve*Ve/VeOe) + &
          (1.0/MiVin) * (Vi*Vi/ViOi)) * E_Density * q2

     sigmah(:,:,:,iBlock) = ((1.0/MeVen) * (Ve*e_gyro/VeOe) - &
          (1.0/MiVin) * (Vi*i_gyro/ViOi)) * E_Density * q2

     PedersenConductance(:,:,iBlock) = 0.0
     HallConductance(:,:,iBlock)     = 0.0

     do k=-1,nAlts+2
        PedersenConductance(:,:,iBlock) = PedersenConductance(:,:,iBlock) + &
             sigmap(:,:,k,iBlock) * dAlt_GB(:,:,k,iBlock)
        HallConductance(:,:,iBlock)     = HallConductance(:,:,iBlock)     + &
             sigmah(:,:,k,iBlock) * dAlt_GB(:,:,k,iBlock)
     enddo

     call report("Starting Conductances",2)
     do k=-1,nAlts+2
        do j=-1,nLats+2
           do i=-1,nLons+2

              sigmap_d1d1_d(i,j,k,iBlock) = sigmap(i,j,k,iBlock) * &
                   sum(b0_d1(i,j,k,:,iBlock)**2) / b0_cD(i,j,k,iBlock)

              sigmap_d2d2_d(i,j,k,iBlock) = sigmap(i,j,k,iBlock) * &
                   sum(b0_d2(i,j,k,:,iBlock)**2) / b0_cD(i,j,k,iBlock)

              sigmap_d1d2_d(i,j,k,iBlock) = sigmap(i,j,k,iBlock) * &
                   sum(b0_d1(i,j,k,:,iBlock)*b0_d2(i,j,k,:,iBlock)) / b0_cD(i,j,k,iBlock)

              ! Don't know why Richmond names these Ue1 and Ue2, when the
              ! paper describes them as being Uei = di . U
              ue1(i,j,k) = sum(Velocity(i,j,k,:,iBlock) * b0_d1(i,j,k,:,iBlock))

              ue2(i,j,k) = sum(Velocity(i,j,k,:,iBlock) * b0_d2(i,j,k,:,iBlock))

              Ed1(i,j,k) = sum(EField(i,j,k,:) * b0_d1(i,j,k,:,iBlock)) / &
                   sqrt(sum(b0_d1(i,j,k,:,iBlock)**2))

              Ed2(i,j,k) = sum(EField(i,j,k,:) * b0_d2(i,j,k,:,iBlock)) / &
                   sqrt(sum(b0_d2(i,j,k,:,iBlock)**2))

              ! These are from eqns 5.19 and 5.20
              kmp(i,j,k,iBlock) = ue2(i,j,k) * sigmap_d1d1_d(i,j,k,iBlock) + &
                   (sigmah(i,j,k,iBlock) - sigmap_d1d2_d(i,j,k,iBlock)) * ue1(i,j,k)

              kml(i,j,k,iBlock) = (sigmah(i,j,k,iBlock) + sigmap_d1d2_d(i,j,k,iBlock)) * ue2(i,j,k) - &
                   ue1(i,j,k) * sigmap_d2d2_d(i,j,k,iBlock)

              ! The Capital D is removed from the sigmah, since the d1d?_d are /D...
              je1(i,j,k,iBlock) = &
                   sigmap_d1d1_d(i,j,k,iBlock) &
                   *(Ed1(i,j,k) + ue2(i,j,k)*b0_be3(i,j,k,iBlock)) + &
                   (sigmap_d1d2_d(i,j,k,iBlock) - sigmah(i,j,k,iBlock)) &
                   *(Ed2(i,j,k) + ue1(i,j,k)*b0_be3(i,j,k,iBlock))

              je2(i,j,k,iBlock) = &
                   sigmap_d2d2_d(i,j,k,iBlock) &
                   *(Ed2(i,j,k) - ue1(i,j,k)*b0_be3(i,j,k,iBlock)) + &
                   (sigmap_d1d2_d(i,j,k,iBlock) + sigmah(i,j,k,iBlock)) &
                   *(Ed1(i,j,k) + ue2(i,j,k)*b0_be3(i,j,k,iBlock))

              SigmaR(i,j,k,iEast_,iEast_)   =  sigmap(i,j,k,iBlock)
              SigmaR(i,j,k,iEast_,iNorth_)  = -sigmah(i,j,k,iBlock) * sin(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iEast_,iUp_)     =  sigmah(i,j,k,iBlock) * cos(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iNorth_,iEast_)  = -SigmaR(i,j,k,iEast_,iNorth_)
              SigmaR(i,j,k,iNorth_,iNorth_) = &
                   sigmap(i,j,k,iBlock)*sin(DipAngle(i,j,k,iBlock))**2+&
                   sigma0(i,j,k)*cos(DipAngle(i,j,k,iBlock))**2
              SigmaR(i,j,k,iNorth_,iUp_)    = &
                   (sigma0(i,j,k) - sigmap(i,j,k,iBlock)) * &
                   sin(DipAngle(i,j,k,iBlock)) * cos(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iUp_,iEast_)     = -SigmaR(i,j,k,iEast_,iUp_)
              SigmaR(i,j,k,iUp_,iNorth_)    =  SigmaR(i,j,k,iNorth_,iUp_)
              SigmaR(i,j,k,iUp_,iUp_)       = &
                   sigmap(i,j,k,iBlock)*cos(DipAngle(i,j,k,iBlock))**2+&
                   sigma0(i,j,k)*sin(DipAngle(i,j,k,iBlock))**2
           enddo
        enddo
     enddo

     UxB(:,:,:,iEast_)  =  &
          Velocity(:,:,:,iNorth_,iBlock)* B0(:,:,:,iUp_,iBlock)    - &
          Velocity(:,:,:,iUp_,iBlock)   * B0(:,:,:,iNorth_,iBlock)

     UxB(:,:,:,iNorth_) = -( &
          Velocity(:,:,:,iEast_,iBlock) * B0(:,:,:,iUp_,iBlock)    - &
          Velocity(:,:,:,iUp_,iBlock)   * B0(:,:,:,iEast_,iBlock))

     UxB(:,:,:,iUp_)    =  &
          Velocity(:,:,:,iEast_,iBlock) * B0(:,:,:,iNorth_,iBlock) - &
          Velocity(:,:,:,iNorth_,iBlock)* B0(:,:,:,iEast_,iBlock)

     Ju(:,:,:,iEast_)  = &
          SigmaR(:,:,:,iEast_,iEast_)   * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iEast_,iNorth_)  * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iEast_,iUp_)     * UxB(:,:,:,iUp_)
  
     Ju(:,:,:,iNorth_) = &
          SigmaR(:,:,:,iNorth_,iEast_)  * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iNorth_,iNorth_) * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iNorth_,iUp_)    * UxB(:,:,:,iUp_)
  
     Ju(:,:,:,iUp_)    = &
          SigmaR(:,:,:,iUp_,iEast_)     * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iUp_,iNorth_)    * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iUp_,iUp_)       * UxB(:,:,:,iUp_)

     ! We want to take the divergence of Ju perpendicular to the
     ! magnetic field line.  

     do iAlt = -1, nAlts+2
        do iLat = -1, nLats+2
           do iLon = -1, nLons+2
              JuDotB(iLon, iLat, iAlt) = sum( Ju(iLon,iLat,iAlt,:)* &
                   B0(iLon,iLat,iAlt,:,iBlock)/B0(iLon,iLat,iAlt,iMag_,iBlock) )
           enddo
        enddo
     enddo

     do iDir = 1, 3
        call UAM_Gradient_GC(Ju(:,:,:,iDir), Gradient_GC, iBlock)
        DivJu(:,:,1:nAlts,iBlock) = DivJu(:,:,1:nAlts,iBlock) + Gradient_GC(:,:,1:nAlts,iDir)

        call UAM_Gradient_GC(JuDotB(:,:,:), Gradient_GC, iBlock)
        DivJu(:,:,1:nAlts,iBlock) = DivJu(:,:,1:nAlts,iBlock) - Gradient_GC(:,:,1:nAlts,iDir) * &
             B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)
     enddo

     do k=1,nAlts
        DivJuAlt(:,:,iBlock) = DivJuAlt(:,:,iBlock) + DivJu(:,:,k,iBlock)*dAlt_GB(:,:,k,iBlock)
     enddo

  enddo  ! block loop

  call report("Starting Magfield Traces",2)

  call MMT_Integrate(sigmap_d1d1_d,SigmaPP)
  call MMT_Integrate(sigmap_d2d2_d,SigmaLL)
  call MMT_Integrate(sigmap_d1d2_d,SigmaCC)
  call MMT_Integrate(sigmah,HallFieldLine)
  call MMT_Integrate(sigmap,PedersenFieldLine)
  call MMT_Integrate(DivJu,DivJuFieldLine)
  call MMT_Integrate(kmp,KDpm)
  call MMT_Integrate(kml,KDlm)
  call MMT_Integrate(je1,Kpm)
  call MMT_Integrate(je2,Klm)

  do iBlock = 1, nBlocks
     call calc_mltlocal

     do i=1,nMagLons
        do j=1,nMagLats

           mlatMC = MagLatMC(i,j)
           mltMC  = MagLocTimeMC(i,j)
           gLatMC = GeoLatMC(i,j)
           gLonMC = GeoLonMC(i,j)

           call find_mag_point(jul, shl, spl, length, spp, sll, scc, &
                kdpm_s, kdlm_s, be3, kpm_s, klm_s)

           !write(*,*) "find_mag_point : ", i,j,mLatMC, mltMC, length

           if (length > 0) then
              DivJuAltMC(i,j)      = jul
              SigmaHallMC(i,j)     = shl
              SigmaPedersenMC(i,j) = spl
              LengthMC(i,j)        = length

              sinim = abs(2.0 * sin(MagLatMC(i,j)*pi/180) / &
                   sqrt(4.0 - 3.0 * cos(MagLatMC(i,j)*pi/180)))

              SigmaPPMC(i,j) = spp * sinim
              SigmaLLMC(i,j) = sll / sinim
              SigmaHHMC(i,j) = shl
              SigmaCCMC(i,j) = scc
              if (MagLatMC(i,j) > 0.0) then
                 SigmaPLMC(i,j) = + (shl - scc)
                 SigmaLPMC(i,j) = - (shl + scc)
              else
                 SigmaPLMC(i,j) = - (shl - scc)
                 SigmaLPMC(i,j) = + (shl + scc)
              endif

              KDlmMC(i,j) = -sign(1.0,MagLatMC(i,j)) * kdlm_s * be3
              KDpmMC(i,j) = kdpm_s * be3 * abs(sinim)

              KlmMC(i,j)  = -sign(1.0,MagLatMC(i,j)) * klm_s
              KpmMC(i,j)  = kpm_s * abs(sinim)
           endif

        enddo
     enddo

  enddo  ! block loop
  
  if (iDebugLevel > 2) write(*,*) "===> Beginning Sum of Electrodynamics"

  DivJuAltMC(nMagLons+1,:) = DivJuAltMC(1,:)
  SigmaHallMC(nMagLons+1,:) = SigmaHallMC(1,:)
  SigmaPedersenMC(nMagLons+1,:) = SigmaPedersenMC(1,:)
  LengthMC(nMagLons+1,:) = LengthMC(1,:)

  SigmaPPMC(nMagLons+1,:) = SigmaPPMC(1,:)
  SigmaLLMC(nMagLons+1,:) = SigmaLLMC(1,:)
  SigmaHHMC(nMagLons+1,:) = SigmaHHMC(1,:)
  SigmaCCMC(nMagLons+1,:) = SigmaCCMC(1,:)
  SigmaPLMC(nMagLons+1,:) = SigmaPLMC(1,:)
  SigmaLPMC(nMagLons+1,:) = SigmaLPMC(1,:)

  KDlmMC(nMagLons+1,:) = KDlmMC(1,:)
  KDpmMC(nMagLons+1,:) = KDpmMC(1,:)
  
  KlmMC(nMagLons+1,:) = KlmMC(1,:)
  KpmMC(nMagLons+1,:) = KpmMC(1,:)
  
  bs = nMagLats * (nMagLons+1)

  MagBufferMC = DivJuAltMC
  call MPI_AllREDUCE(MagBufferMC, DivJuAltMC,      bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaHallMC
  call MPI_AllREDUCE(MagBufferMC, SigmaHallMC,     bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPedersenMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPedersenMC, bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = LengthMC
  call MPI_AllREDUCE(MagBufferMC, LengthMC,        bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPPMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPPMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaLLMC
  call MPI_AllREDUCE(MagBufferMC, SigmaLLMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaHHMC
  call MPI_AllREDUCE(MagBufferMC, SigmaHHMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaCCMC
  call MPI_AllREDUCE(MagBufferMC, SigmaCCMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaLPMC
  call MPI_AllREDUCE(MagBufferMC, SigmaLPMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPLMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPLMC,       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KDlmMC
  call MPI_AllREDUCE(MagBufferMC, KDlmMC,          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KDpmMC
  call MPI_AllREDUCE(MagBufferMC, KDpmMC,          bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KlmMC
  call MPI_AllREDUCE(MagBufferMC, KlmMC,           bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KpmMC
  call MPI_AllREDUCE(MagBufferMC, KpmMC,           bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  do i=1,nMagLons+1
     do j=2,nMagLats-1
        dSigmaLLdlMC(i,j) = 0.5*(SigmaLLMC(i,j+1) - SigmaLLMC(i,j-1))/deltalmc(i,j)
        dSigmaLPdlMC(i,j) = 0.5*(SigmaLPMC(i,j+1) - SigmaLPMC(i,j-1))/deltalmc(i,j)
        dkdlmdlMC(i,j) = cos(MagLatMC(i,j)*pi/180) * &
             0.5*(KDlmMC(i,j+1) - KDlmMC(i,j-1)) / deltalmc(i,j)
        dkdlmdlMC(i,j) = dkdlmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KDlmMC(i,j)
        dklmdlMC(i,j) = cos(MagLatMC(i,j)*pi/180) * &
             0.5*(KlmMC(i,j+1) - KlmMC(i,j-1)) / deltalmc(i,j)
        dklmdlMC(i,j) = dklmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KlmMC(i,j)
     enddo

     dSigmaLLdlMC(i,1) = (SigmaLLMC(i,2) - SigmaLLMC(i,1)) / deltalmc(i,1)
     dSigmaLLdlMC(i,nMagLats) = &
          (SigmaLLMC(i,nMagLats) - SigmaLLMC(i,nMagLats-1)) / deltalmc(i,nMagLats)

     dSigmaLPdlMC(i,1) = (SigmaLPMC(i,2) - SigmaLPMC(i,1)) / deltalmc(i,1)
     dSigmaLPdlMC(i,nMagLats) = &
          (SigmaLPMC(i,nMagLats) - SigmaLPMC(i,nMagLats-1)) / deltalmc(i,nMagLats)

     dkdlmdlMC(i,1) = cos(MagLatMC(i,1)*pi/180) * &
          (KDlmMC(i,2) -KDlmMC(i,1)) / deltalmc(i,1)
     dkdlmdlMC(i,nMagLats) = cos(MagLatMC(i,nMagLats)*pi/180) * &
          (KDlmMC(i,nMagLats) - KDlmMC(i,nMagLats-1)) / deltalmc(i,nMagLats)
     j = 1
     dkdlmdlMC(i,j) = dkdlmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KDlmMC(i,j)
     j = nMagLats
     dkdlmdlMC(i,j) = dkdlmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KDlmMC(i,j)

     dklmdlMC(i,1) = cos(MagLatMC(i,1)*pi/180) * &
          (KlmMC(i,2) -KlmMC(i,1)) / deltalmc(i,1)
     dklmdlMC(i,nMagLats) = cos(MagLatMC(i,nMagLats)*pi/180) * &
          (KlmMC(i,nMagLats) - KlmMC(i,nMagLats-1)) / deltalmc(i,nMagLats)
     j = 1
     dklmdlMC(i,j) = dklmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KlmMC(i,j)
     j = nMagLats
     dklmdlMC(i,j) = dklmdlMC(i,j) - sin(MagLatMC(i,j)*pi/180) * KlmMC(i,j)

  enddo

  do j=1,nMagLats
     do i=2,nMagLons
        dSigmaPLdpMC(i,j) = 0.5*(SigmaPLMC(i+1,j) - SigmaPLMC(i-1,j)) / deltapmc(i,j)
        dSigmaPPdpMC(i,j) = 0.5*(SigmaPPMC(i+1,j) - SigmaPPMC(i-1,j)) / deltapmc(i,j)
        dKDpmdpMC(i,j)    = 0.5*(KDpmMC(i+1,j)    - KDpmMC(i-1,j)   ) / deltapmc(i,j)
        dKpmdpMC(i,j)     = 0.5*(KpmMC(i+1,j)     - KpmMC(i-1,j)    ) / deltapmc(i,j)
     enddo
     dSigmaPLdpMC(1,j)          = (SigmaPLMC(2,j) - SigmaPLMC(1,j))/deltapmc(1,j)
     dSigmaPLdpMC(nMagLons+1,j) = dSigmaPLdpMC(1,j)
     dSigmaPPdpMC(1,j)          = (SigmaPPMC(2,j) - SigmaPPMC(1,j))/deltapmc(1,j)
     dSigmaPPdpMC(nMagLons+1,j) = dSigmaPPdpMC(1,j)
     dKDpmdpMC(1,j)             = (KDpmMC(2,j) - KDpmMC(1,j))/deltapmc(1,j)
     dKDpmdpMC(nMagLons+1,j)    = dKDpmdpMC(1,j)
     dKpmdpMC(1,j)              = (KpmMC(2,j) - KpmMC(1,j))/deltapmc(1,j)
     dKpmdpMC(nMagLons+1,j)     = dKpmdpMC(1,j)
  enddo

  solver_a_mc = 4 * deltalmc**2 * sigmappmc / cos(MagLatMC*pi/180)
  solver_b_mc = 4 * deltapmc**2 * cos(MagLatMC*pi/180) * sigmallmc
  solver_c_mc = deltalmc * deltapmc * (SigmaPLmc + SigmaLPmc)  
     
  solver_d_mc = 2.0 * deltalmc * deltapmc**2 * ( &
       sign(1.0, MagLatMC) * ( dSigmaPLdpMC + cos(MagLatMC*pi/180) * dSigmaLLdlMC)&
       - sin(MagLatMC*pi/180) * sigmallmc)
  solver_e_mc = 2.0 * sign(1.0, MagLatMC) * deltalmc**2 * deltapmc * ( &
       dSigmaPPdpMC / cos(MagLatMC*pi/180) + dSigmaLPdlMC)
!  solver_s_mc =  4 * deltalmc**2 * deltapmc**2 * (RBody) * &
!       (dkdlmdlMC + dKDpmdpMC - dklmdlMC - dkpmdpMC)
  solver_s_mc = 4 * deltalmc**2 * deltapmc**2 * (RBody) * &
       (dkdlmdlMC + dKDpmdpMC)
     
!  solver_d_mc(:,42:49) = 0.0
!  solver_e_mc(:,42:49) = 0.0
     
!  OldPotMc = 0.0
  DynamoPotentialMC = 0.0

  ! Fill in the diagonal vectors
  iI = 0

  nX = (nMagLats-2) * (nMagLons)

  allocate( x(nX), y(nX), rhs(nX), b(nX), &
       d_I(nX), e_I(nX), e1_I(nX), f_I(nX), f1_I(nX) )

  do iLat=2,nMagLats-1
     do iLon=1,nMagLons

        iI = iI + 1

        ! Right hand side
        b(iI)    = solver_s_mc(iLon, iLat)

        ! Initial Guess
        x(iI)    = DynamoPotentialMC(iLon, iLat)
     
        ! i,j
        d_I(iI)  = -2*(solver_a_mc(iLon, iLat)+solver_b_mc(iLon, iLat))

        ! ilon-1, ilat
        e_I(iI)  = solver_a_mc(iLon, iLat)-solver_e_mc(iLon, iLat)
  
        ! ilon+1, ilat
        f_I(iI)  = solver_a_mc(iLon, iLat)+solver_e_mc(iLon, iLat)
       
        ! ilon, ilat-1
        e1_I(iI) = solver_b_mc(iLon, iLat)-solver_d_mc(iLon, iLat)
       
        ! ilon, ilat+1
        f1_I(iI) = solver_b_mc(iLon, iLat)+solver_d_mc(iLon, iLat)
  
        if (iLat == 2)          e1_I(iI) = 0.0
        if (iLat == nMagLats-1) f1_I(iI) = 0.0
        if (iLon == 1)          e_I(iI)  = 0.0
        if (iLon == nMagLons)   f_I(iI)  = 0.0

     end do
  end do

  Rhs = b

!!!  write(*,*) "prehepta"
!!!  ! A -> LU
!!!
!!!  write(*,*) "pre : ", &
!!!       sum(b),sum(abs(b)),sum(x),sum(d_I),sum(e_I),sum(f_I),&
!!!       sum(e1_I),sum(f1_I)
!!!
!!!  call prehepta(nX,1,nMagLons,nX,-0.5,d_I,e_I,f_I,e1_I,f1_I)
!!!
!!!  ! Left side preconditioning: U^{-1}.L^{-1}.A.x = U^{-1}.L^{-1}.rhs
!!!
!!!  ! rhs'=U^{-1}.L^{-1}.rhs
!!!  write(*,*) "Lhepta"
!!!  call Lhepta(       nX,1,nMagLons,nX,b,d_I,e_I,e1_I)
!!!  write(*,*) "Uhepta"
!!!  call Uhepta(.true.,nX,1,nMagLons,nX,b,    f_I,f1_I)

  MaxIteration = 500
  nIteration = 0
  iError = 0
  if (iDebugLevel > 2) then
     DoTestMe = .true.
  else
     DoTestMe = .false.
  endif

  Residual = 0.01
  call gmres(matvec_gitm,b,x,.true.,nX,&
       MaxIteration,Residual,'rel',nIteration,iError,DoTestMe)

  if (iDebugLevel > 0) &
       write(*,*) "=> gmres : ",MaxIteration,Residual, nIteration, iError

  iI = 0
  do iLat=2,nMagLats-1
     do iLon=1,nMagLons
        iI = iI + 1
        DynamoPotentialMC(iLon, iLat) = x(iI)
     enddo
  enddo

!  DynamoPotentialMC(1,:) = DynamoPotentialMC(nMagLons,:)
  DynamoPotentialMC(nMagLons+1,:) = DynamoPotentialMC(1,:)

  DynamoPotentialMC(:,1) = 0.0
  DynamoPotentialMC(:,nMagLats) = 0.0


  if (iDebugLevel > 0) &
       write(*,*) "=> CPCP of Dynamo : ", &
       maxval(dynamopotentialmc)-minval(dynamopotentialmc)

  where (SigmaHallMC < 0.1) SigmaHallMC = 0.1
  where (SigmaPedersenMC < 0.1) SigmaPedersenMC = 0.1

  if(allocated(b)) deallocate(x, y, b, rhs, d_I, e_I, f_I, e1_I, f1_I)

  call end_timing("calc_electrodyn")

contains

  !\
  ! ------------------------------------------------------------------
  !/

  subroutine calc_mltlocal

    magloctime_local = mod(MLT(0:nLons+1,0:nLats+1,1)+24.0,24.0)
    do iLat = 0,nLats+1
       do iLon = 0, nLons
          dl = magloctime_local(iLon+1,iLat)-magloctime_local(iLon,iLat)
          if (dl < 0) then
             if (iLon == 0) then
                magloctime_local(iLon,iLat) = &
                     magloctime_local(iLon,iLat)-24.0
             else
                if (iLon == nLons) then
                   magloctime_local(iLon+1,iLat) = &
                        magloctime_local(iLon+1,iLat)+24.0
                endif
             endif
          endif
       enddo
    enddo

  end subroutine calc_mltlocal


  !\
  ! ------------------------------------------------------------------
  !/

  !\
  ! Since this is a contains subroutine, we don't need to pass in
  ! mlatMC or mltMC or iBlock.
  !/

  subroutine find_mag_point(juline, shline, spline, length, &
       sppline, sllline, sccline, kdpline, kdlline, be3, &
       kpline, klline)

    real, intent(out) :: juline, shline, spline, &
         sppline, sllline, sccline, kdpline, kdlline, &
         be3, kpline, klline

    integer :: ii, jj

    real :: dip, dec, length, mfac, lfac, gmlt, gmlt2, cdip, sdip, mt

    logical :: IsFound

    juline = 0.0
    shline = 0.0
    spline = 0.0
    length = 0.0

    sppline=0.
    sllline=0.
    sccline=0.
    kdpline=0.
    kdlline=0.
    be3=0.
    kpline=0.
    klline=0.

    if (gLatMC > Latitude(nLats+1,iBlock)) return
    if (gLatMC < Latitude(0,iBlock)) return

    if (gLonMC > Longitude(nLons+1,iBlock)) return
    if (gLonMC < Longitude(0,iBlock)) return

    IsFound = .false.
    jj = -1
    do while (.not.IsFound)
       jj = jj + 1
       if ( gLatMC >= Latitude(jj  ,iBlock) .and. &
            gLatMC <  Latitude(jj+1,iBlock)) IsFound = .true.
    enddo

    IsFound = .false.
    ii = -1
    do while (.not.IsFound)
       ii = ii + 1
       if ( gLonMC >= Longitude(ii  ,iBlock) .and. &
            gLonMC <  Longitude(ii+1,iBlock)) IsFound = .true.
    enddo

    mfac = (gLonMC - Longitude(ii  ,iBlock)) / &
         (Longitude(ii+1,iBlock)-Longitude(ii  ,iBlock))

    lfac = (gLatMC - Latitude(jj  ,iBlock)) / &
         (Latitude(jj+1,iBlock)-Latitude(jj  ,iBlock))

    call interpolate_localB(HallFieldLine,     mfac, lfac, ii, jj, shline)
    call interpolate_localB(PedersenFieldLine, mfac, lfac, ii, jj, spline)
    call interpolate_localB(DivJuFieldLine,    mfac, lfac, ii, jj, juline)
    call interpolate_localB(LengthFieldLine,   mfac, lfac, ii, jj, length)

    call interpolate_localB(SigmaPP, mfac, lfac, ii, jj, sppline)
    call interpolate_localB(SigmaLL, mfac, lfac, ii, jj, sllline)
    call interpolate_localB(SigmaCC, mfac, lfac, ii, jj, sccline)

    call interpolate_localB(KDpm, mfac, lfac, ii, jj, kdpline)
    call interpolate_localB(KDlm, mfac, lfac, ii, jj, kdlline)
    
    call interpolate_localB(Kpm, mfac, lfac, ii, jj, kpline)
    call interpolate_localB(Klm, mfac, lfac, ii, jj, klline)

    call interpolate_local(b0_be3(:,:,1,iBlock), mfac, lfac, ii, jj, be3)

    if (sppline > 1000.) write(*,*) "sppline : ", &
         mfac, lfac, ii, jj, sppline, &
         SigmaPP(ii,jj  ,iBlock), SigmaPP(ii+1,jj  ,iBlock), &
         SigmaPP(ii,jj+1,iBlock), SigmaPP(ii+1,jj+1,iBlock)

  end subroutine find_mag_point

  subroutine interpolate_local(VarToInter, mfac, lfac, ii, jj, Output)

    real, intent(in)    :: VarToInter(-1:nLons+2,-1:nLats+2), mfac, lfac
    integer, intent(in) :: ii, jj
    real, intent(out)   :: Output
    !----------

    Output = (1.0-mfac)*(1.0-lfac)*VarToInter(ii  ,jj  ) + &
             (    mfac)*(1.0-lfac)*VarToInter(ii+1,jj  ) + &
             (1.0-mfac)*(    lfac)*VarToInter(ii  ,jj+1) + &
             (    mfac)*(    lfac)*VarToInter(ii+1,jj+1)

  end subroutine interpolate_local

  subroutine interpolate_localB(VarToInter, mfac, lfac, ii, jj, Output)

    real, intent(in)    :: VarToInter(-1:nLons+2,-1:nLats+2,nBlocksMax), mfac, lfac
    integer, intent(in) :: ii, jj
    real, intent(out)   :: Output
    !----------

    Output = (1.0-mfac)*(1.0-lfac)*VarToInter(ii  ,jj  ,iBlock) + &
             (    mfac)*(1.0-lfac)*VarToInter(ii+1,jj  ,iBlock) + &
             (1.0-mfac)*(    lfac)*VarToInter(ii  ,jj+1,iBlock) + &
             (    mfac)*(    lfac)*VarToInter(ii+1,jj+1,iBlock)

  end subroutine interpolate_localB

end subroutine UA_calc_electrodynamics

!============================================================================
subroutine matvec_gitm(x_I, y_I, n)

  use ModElectrodynamics
  use ModLinearsolver, ONLY: Uhepta, Lhepta

  implicit none

  ! Calculate y = A.x where A is the (pentadiagonal) matrix

  integer, intent(in) :: n          ! number of unknowns
  real, intent(in) :: x_I(n)        ! vector of unknowns
  real, intent(out):: y_I(n)        ! y = A.x

  integer :: iLat, iLon, i, iLm, iLp
  real :: x_G(nMagLons+1, nMagLats) ! 2D array with ghost cells
  !----------

  ! Put 1D vector into 2D solution
  i = 0;
  do iLat = 2, nMagLats-1
     do iLon = 1, nMagLons
        i = i+1
        x_G(iLon, iLat) = x_I(i)
     enddo
  enddo

  x_G(:,1) = 0.0
  x_G(:,nMagLats) = 0.0

  ! Apply periodic boundary conditions in Psi direction
  x_G(nMagLons+1,:) = x_G(1,:)

  i = 0;
  do iLat = 2, nMagLats-1
     do iLon = 1, nMagLons
        i = i+1

        iLm = iLon-1
        if (iLm == 0) iLm = nMagLons
        iLp = iLon+1
        if (iLp == nMagLons+1) iLp = 0

           y_I(i) = &
             -(2*solver_a_mc(iLon, iLat)+2*solver_b_mc(iLon, iLat)) * &
             x_G(iLon   , iLat  ) + &
             ( solver_a_mc(iLon, iLat)+solver_e_mc(iLon, iLat)) * &
             x_G(iLp , iLat  ) + &
             ( solver_b_mc(iLon, iLat)+solver_d_mc(iLon, iLat)) * &
             x_G(iLon  , iLat+1) + &
             ( solver_c_mc(iLon, iLat)                        ) * &
             x_G(iLp, iLat+1) + &
             (-solver_c_mc(iLon, iLat)                        ) * &
             x_G(iLp, iLat-1) + &
             ( solver_a_mc(iLon, iLat)-solver_e_mc(iLon, iLat)) * &
             x_G(iLm, iLat  ) + &
             ( solver_b_mc(iLon, iLat)-solver_d_mc(iLon, iLat)) * &
             x_G(iLon  , iLat-1)    + &
             (-solver_c_mc(iLon, iLat)                        ) * &
             x_G(iLm, iLat+1) + &
             ( solver_c_mc(iLon, iLat)                        ) * &
             x_G(iLm, iLat-1)

     end do
  end do

!!!  ! Preconditioning: y'= U^{-1}.L^{-1}.y
!!!  call Lhepta(       n,1,nMagLons,n,y_I,d_I,e_I,e1_I)
!!!  call Uhepta(.true.,n,1,nMagLons,n,y_I,    f_I,f1_I)

end subroutine matvec_gitm

