subroutine UA_fill_electrodynamics(UAr2_fac, UAr2_ped, UAr2_hal, &
            UAr2_lats, UAr2_mlts)

  use ModElectrodynamics

  implicit none

  real, dimension(nMagLons+1,nMagLats), intent(out) :: &
       UAr2_fac, UAr2_ped, UAr2_hal, UAr2_lats, UAr2_mlts

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

  implicit none

  integer, intent(out) :: UAi_nMLTs, UAi_nLats

  integer, external :: jday

  integer :: i,j,k,bs, iError, iBlock, iDir, iLon, iLat, iAlt, ip, im

  real :: GeoLat, GeoLon, GeoAlt, xAlt, len, ped, hal
  real :: sp_d1d1_d, sp_d2d2_d, sp_d1d2_d, sh
  real :: xmag, ymag, zmag, bmag, signz, magpot, lShell
  real :: mlatMC, mltMC, jul, shl, spl, length, kdpm_s, kdlm_s, je1_s, je2_s
  real :: kpm_s, klm_s, xstretch, ystretch
  real :: sinIm, spp, sll, shh, scc, sccline, sppline, sllline, shhline, be3

  real :: q2, dju,dl, cD

  real :: magloctime_local(0:nLons+1,0:nLats+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       JuDotB, sigmap_d1d1_d, sigmap_d2d2_d, sigmap_d1d2_d, sigmah, &
       ue1, ue2, kmp, kml, je1, je2, ed1, ed2

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       Gradient_GC

  real :: aLat, aLon, gLat, gLon, Date, sLat, sLon, gLatMC, gLonMC

  real :: residual, oldresidual, a, tmp

  logical :: IsDone, IsFirstTime = .true., DoTestMe

  integer :: iLm, iLp, jLat, iI, MaxIteration, nIteration
  integer :: nX

  external :: matvec_gitm

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

     if (iError /= 0) then
        call CON_stop("Error allocating array DivJuAltMC")
     endif

     date = iStartTime(1) + float(iJulianDay)/float(jday(iStartTime(1),12,31))

     do i=1,nMagLons+1
        if (iDebugLevel > 1) &
             write(*,*) "==> Calculating Apex->Geo", i, nMagLons+1
        do j=1,nMagLats

           MagLatMC(i,j)     = 180.0 * (float(j)-0.5) / float(nMagLats) - 90.0
           MagLonMC(i,j)     = 360.0 * float(i-1) / float(nMagLons)

           ! we could also use the stretch_grid function....
!MAGSTRECH       xstretch =(float(j)-0.5)/float(nMagLats/2) - 1.0
!MAGSTRECH       call stretch_grid(xstretch,ystretch)
!MAGSTRECH       MagLatMC(i,j) = ystretch * 90.0

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

              call apex_to_geo(date, aLat, aLon, AltMinIono, &
                   gLat, gLon, sLat, sLon)

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

           call magloctm(MagLonMC(i,j),SubsolarLatitude,   &
                SubsolarLongitude,  &
                MagneticPoleColat, &
                MagneticPoleLon,mltMC)

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

     Sigma_0 = q2 * E_Density / (1.0/MeVen + 1.0/MiVin)

     Sigma_Pedersen = ((1.0/MeVen) * (Ve*Ve/VeOe) + &
          (1.0/MiVin) * (Vi*Vi/ViOi)) * E_Density * q2

     Sigma_Hall = ((1.0/MeVen) * (Ve*e_gyro/VeOe) - &
          (1.0/MiVin) * (Vi*i_gyro/ViOi)) * E_Density * q2

     PedersenConductance(:,:,iBlock) = 0.0
     HallConductance(:,:,iBlock)     = 0.0

     do k=-1,nAlts+2
        PedersenConductance(:,:,iBlock) = PedersenConductance(:,:,iBlock) + &
             Sigma_Pedersen(:,:,k)*dAlt_GB(:,:,k,iBlock)
        HallConductance(:,:,iBlock)     = HallConductance(:,:,iBlock)     + &
             Sigma_Hall(:,:,k)    *dAlt_GB(:,:,k,iBlock)
     enddo

     call report("Starting Conductances",2)
     do k=-1,nAlts+2
        do j=-1,nLats+2
           do i=-1,nLons+2

              sigmap_d1d1_d(i,j,k) = &
                   Sigma_Pedersen(i,j,k) * &
                   sum(b0_d1(i,j,k,:,iBlock)**2) / &
                   b0_cD(i,j,k,iBlock)

              sigmap_d2d2_d(i,j,k) = &
                   Sigma_Pedersen(i,j,k) * &
                   sum(b0_d2(i,j,k,:,iBlock)**2) / &
                   b0_cD(i,j,k,iBlock)

              sigmah(i,j,k) = Sigma_Hall(i,j,k)
              sigmap_d1d2_d(i,j,k) = &
                   Sigma_Pedersen(i,j,k) * &
                   sum(b0_d1(i,j,k,:,iBlock) * &
                       b0_d2(i,j,k,:,iBlock))/&
                   b0_cD(i,j,k,iBlock)

              ! Don't know why Richmond names these Ue1 and Ue2, when the
              ! paper describes them as being Uei = di . U
              ue1(i,j,k) = sum( &
                   Velocity(i,j,k,:,iBlock) * b0_d1(i,j,k,:,iBlock))

              ue2(i,j,k) = sum( &
                   Velocity(i,j,k,:,iBlock) * b0_d2(i,j,k,:,iBlock))

              Ed1(i,j,k) = sum( &
                   EField(i,j,k,:) * b0_d1(i,j,k,:,iBlock)) / &
                   sqrt(sum(b0_d1(i,j,k,:,iBlock)**2))

              Ed2(i,j,k) = sum( &
                   EField(i,j,k,:) * b0_d2(i,j,k,:,iBlock)) / &
                   sqrt(sum(b0_d2(i,j,k,:,iBlock)**2))

              ! These are from eqns 5.19 and 5.20
              kmp(i,j,k) = ue2(i,j,k) * sigmap_d1d1_d(i,j,k) + &
                   (sigmah(i,j,k) - sigmap_d1d2_d(i,j,k)) * ue1(i,j,k)

              kml(i,j,k) = (sigmah(i,j,k) + sigmap_d1d2_d(i,j,k)) * ue2(i,j,k) - &
                   ue1(i,j,k) * sigmap_d1d1_d(i,j,k)

              ! The Capital D is removed from the sigmah, since the d1d?_d are /D...
              je1(i,j,k) = sigmap_d1d1_d(i,j,k)*(Ed1(i,j,k) + ue2(i,j,k)*b0_be3(i,j,k,iBlock)) + &
                   (sigmap_d1d2_d(i,j,k) - sigmah(i,j,k)) * &
                   (Ed2(i,j,k) + ue1(i,j,k)*b0_be3(i,j,k,iBlock))

              je2(i,j,k) = sigmap_d2d2_d(i,j,k)*(Ed2(i,j,k) - ue1(i,j,k)*b0_be3(i,j,k,iBlock)) + &
                   (sigmap_d1d2_d(i,j,k) + sigmah(i,j,k)) * &
                   (Ed1(i,j,k) + ue2(i,j,k)*b0_be3(i,j,k,iBlock))

              SigmaR(i,j,k,iEast_,iEast_)   = &
                   Sigma_Pedersen(i,j,k)
              SigmaR(i,j,k,iEast_,iNorth_)  = &
                   -Sigma_Hall(i,j,k) * sin(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iEast_,iUp_)     = &
                   Sigma_Hall(i,j,k) * cos(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iNorth_,iEast_)  = &
                   -SigmaR(i,j,k,iEast_,iNorth_)
              SigmaR(i,j,k,iNorth_,iNorth_) = &
                   Sigma_Pedersen(i,j,k)*sin(DipAngle(i,j,k,iBlock))**2+&
                   Sigma_0(i,j,k)*cos(DipAngle(i,j,k,iBlock))**2
              SigmaR(i,j,k,iNorth_,iUp_)    = &
                   (Sigma_0(i,j,k) - Sigma_Pedersen(i,j,k)) * &
                   sin(DipAngle(i,j,k,iBlock)) * cos(DipAngle(i,j,k,iBlock))
              SigmaR(i,j,k,iUp_,iEast_)     = &
                   -SigmaR(i,j,k,iEast_,iUp_)
              SigmaR(i,j,k,iUp_,iNorth_)    = &
                   SigmaR(i,j,k,iNorth_,iUp_)
              SigmaR(i,j,k,iUp_,iUp_)       = &
                   Sigma_Pedersen(i,j,k)*cos(DipAngle(i,j,k,iBlock))**2+&
                   Sigma_0(i,j,k)*sin(DipAngle(i,j,k,iBlock))**2
           enddo
        enddo
     enddo

     UxB(:,:,:,iEast_)  =  &
          Velocity(:,:,:,iNorth_,iBlock)*B0(:,:,:,iUp_,iBlock)    - &
          Velocity(:,:,:,iUp_,iBlock)   *B0(:,:,:,iNorth_,iBlock)

     UxB(:,:,:,iNorth_) = -( &
          Velocity(:,:,:,iEast_,iBlock) *B0(:,:,:,iUp_,iBlock)    - &
          Velocity(:,:,:,iUp_,iBlock)   *B0(:,:,:,iEast_,iBlock))

     UxB(:,:,:,iUp_)    =  &
          Velocity(:,:,:,iEast_,iBlock) *B0(:,:,:,iNorth_,iBlock) - &
          Velocity(:,:,:,iNorth_,iBlock)*B0(:,:,:,iEast_,iBlock)

     Ju(:,:,:,iEast_)  = &
          SigmaR(:,:,:,iEast_,iEast_)  * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iEast_,iNorth_) * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iEast_,iUp_)    * UxB(:,:,:,iUp_)
  
     Ju(:,:,:,iNorth_) = &
          SigmaR(:,:,:,iNorth_,iEast_)  * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iNorth_,iNorth_) * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iNorth_,iUp_)    * UxB(:,:,:,iUp_)
  
     Ju(:,:,:,iUp_)    = &
          SigmaR(:,:,:,iUp_,iEast_)  * UxB(:,:,:,iEast_) + &
          SigmaR(:,:,:,iUp_,iNorth_) * UxB(:,:,:,iNorth_) + &
          SigmaR(:,:,:,iUp_,iUp_)    * UxB(:,:,:,iUp_)

     ! We want to take the divergence of Ju perpendicular to the
     ! magnetic field line.  

     do iAlt = -1, nAlts+2
        do iLat = -1, nLats+2
           do iLon = -1, nLons+2

              JuDotB(iLon, iLat, iAlt) = sum( &
                   Ju(iLon,iLat,iAlt,:)* &
                   B0(iLon,iLat,iAlt,:,iBlock)/B0(iLon,iLat,iAlt,iMag_,iBlock))

           enddo
        enddo
     enddo

     do iDir = 1, 3

        call UAM_Gradient_GC(Ju(:,:,:,iDir), Gradient_GC, iBlock)
        DivJu(:,:,:) = DivJu(:,:,:) + Gradient_GC(:,:,1:nAlts,iDir)

        call UAM_Gradient_GC(JuDotB(:,:,:), Gradient_GC, iBlock)
        DivJu(:,:,:) = DivJu(:,:,:) - Gradient_GC(:,:,1:nAlts,iDir) * &
             B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)

     enddo

     DivJuAlt = 0.0
     do k=1,nAlts
        DivJuAlt(:,:) = DivJuAlt(:,:) + DivJu(:,:,k)*dAlt_GB(:,:,k,iBlock)
     enddo

     PedersenFieldLine = 0.0
     HallFieldLine     = 0.0
     LengthFieldLine   = 0.0
     DivJuFieldLine    = 0.0

     SigmaPP = 0.0
     SigmaLL = 0.0
     SigmaHH = 0.0
     SigmaCC = 0.0

     KDpm = 0.0
     KDlm = 0.0

     Kpm = 0.0
     Klm = 0.0

     call report("Starting Magfield Traces",2)

     do iLon = -1, nLons+2
        do iLat = -1, nLats+2

           GeoLat = Latitude(iLat, iBlock)
           GeoLon = Longitude(iLon,iBlock)

           if (GeoLat > pi) then
              GeoLat = 2.0*pi - GeoLat
              GeoLon = mod(GeoLon + pi,twopi)
           endif
              
           if (GeoLat < -pi) then
              GeoLat = -pi - GeoLat
              GeoLon = mod(GeoLon + pi,twopi)
           endif

           GeoAlt = Altitude_GB(iLon,iLat,1,iBlock)
           IsDone = .false.
           len = 250.0
           xAlt = 1.0
           iAlt = 1

           CALL get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0, &
                   XMAG,YMAG,ZMAG)
           signz = sign(1.0,zmag)

           do while (.not. IsDone)

              sp_d1d1_d = &
                          xAlt  * sigmap_d1d1_d(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * sigmap_d1d1_d(iLon, iLat, iAlt+1)

              sp_d2d2_d = &
                          xAlt  * sigmap_d2d2_d(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * sigmap_d2d2_d(iLon, iLat, iAlt+1)

              sp_d1d2_d = &
                          xAlt  * sigmap_d1d2_d(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * sigmap_d1d2_d(iLon, iLat, iAlt+1)

              sh        = &
                          xAlt  * sigmah(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * sigmah(iLon, iLat, iAlt+1)

              kdpm_s     = &
                          xAlt  * kmp(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * kmp(iLon, iLat, iAlt+1)

              kdlm_s     = &
                          xAlt  * kml(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * kml(iLon, iLat, iAlt+1)

              je1_s     = &
                          xAlt  * je1(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * je1(iLon, iLat, iAlt+1)

              je2_s     = &
                          xAlt  * je2(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * je2(iLon, iLat, iAlt+1)

              ped = &
                          xAlt  * Sigma_Pedersen(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * Sigma_Pedersen(iLon, iLat, iAlt+1)
              hal = &
                        xAlt  * Sigma_Hall(iLon, iLat, iAlt) + &
                   (1.0-xAlt) * Sigma_Hall(iLon, iLat, iAlt+1)
              dju = &
                        xAlt  * DivJu(iLon, iLat, iAlt) + &
                   (1.0-xAlt) * DivJu(iLon, iLat, iAlt+1)

              SigmaPP(iLon, iLat) = SigmaPP(iLon, iLat) + len * sp_d1d1_d

              if(SigmaPP(iLon, iLat) > 1000.) write(*,*) "integrating :",&
                   iLon, iLat, iAlt, SigmaPP(iLon, iLat), len, sp_d1d1_d, &
                   Sigma_Pedersen(iLon, iLat, iAlt), &
                   b0_d1(iLon,iLat,iAlt,1:3,iBlock), &
                   b0_cD(iLon,iLat,iAlt,iBlock), iBlock, iProc

              if(SigmaPP(iLon, iLat) > 1000.) write(*,*) "neighbor :",&
                   iLon, iLat, iAlt, SigmaPP(iLon, iLat), len, sp_d1d1_d, &
                   Sigma_Pedersen(iLon-1, iLat, iAlt), &
                   b0_d1(iLon-1,iLat,iAlt,1:3,iBlock), &
                   b0_cD(iLon-1,iLat,iAlt,iBlock), iBlock, iProc

              SigmaLL(iLon, iLat) = SigmaLL(iLon, iLat) + len * sp_d2d2_d
              SigmaHH(iLon, iLat) = SigmaHH(iLon, iLat) + len * sh
              SigmaCC(iLon, iLat) = SigmaCC(iLon, iLat) + len * sp_d1d2_d

              KDpm(iLon, iLat) = KDpm(iLon, iLat) + len * kdpm_s
              KDlm(iLon, iLat) = KDlm(iLon, iLat) + len * kdlm_s

              Kpm(iLon, iLat) = Kpm(iLon, iLat) + len * je1_s
              Klm(iLon, iLat) = Klm(iLon, iLat) + len * je2_s

              PedersenFieldLine(iLon, iLat) = &
                   PedersenFieldLine(iLon, iLat) + len * ped
              HallFieldLine(iLon, iLat) = &
                   HallFieldLine(iLon, iLat) + len * hal
              DivJuFieldLine(iLon, iLat) = &
                   DivJuFieldLine(iLon, iLat) + len * dju

              LengthFieldLine(iLon, iLat) = &
                   LengthFieldLine(iLon, iLat) + len

              CALL get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0,&
                   XMAG,YMAG,ZMAG)

              if (sign(1.0,zmag)*signz < 0) then
                 IsDone = .true.
              else
                 bmag = sqrt(xmag*xmag + ymag*ymag + zmag*zmag)
                 GeoAlt = GeoAlt + abs(zmag)/bmag * len
                 if (GeoAlt > Altitude_GB(iLon,iLat,nAlts,iBlock)) then
                    IsDone = .true.
                 else
                    if (GeoAlt > Altitude_GB(iLon,iLat,iAlt+1,iBlock)) &
                         iAlt = iAlt+1
                    xAlt = (GeoAlt - Altitude_GB(iLon,iLat,iAlt,iBlock)) / &
                         ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) &
                         - Altitude_GB(iLon,iLat,iAlt  ,iBlock))
                    GeoLat = GeoLat + signz*xmag/bmag * len/(RBody + GeoAlt)*pi
                    GeoLon = GeoLon + &
                         signz*ymag/bmag * len/(RBody + GeoAlt)*pi/cos(GeoLon)

                    if (GeoLat > pi) then
                       GeoLat = 2.0*pi - GeoLat
                       GeoLon = mod(GeoLon + pi,twopi)
                    endif
              
                    if (GeoLat < -pi) then
                       GeoLat = -pi - GeoLat
                       GeoLon = mod(GeoLon + pi,twopi)
                    endif

                 endif
              endif

           enddo

        enddo
     enddo

     call calc_mltlocal

     do i=1,nMagLons
        do j=1,nMagLats

           mlatMC = MagLatMC(i,j)
           mltMC  = MagLocTimeMC(i,j)
           gLatMC = GeoLatMC(i,j)
           gLonMC = GeoLonMC(i,j)

           call find_mag_point(jul, shl, spl, length, spp, sll, shh, scc, &
                kdpm_s, kdlm_s, be3, kpm_s, klm_s)

           !write(*,*) "find_mag_points : ", i,j,mLatMC, mltMC, length

           if (length > 0) then
              DivJuAltMC(i,j)      = jul
              SigmaHallMC(i,j)     = shl
              SigmaPedersenMC(i,j) = spl
              LengthMC(i,j)        = length

              sinim = abs(2.0 * sin(MagLatMC(i,j)*pi/180) / &
                   sqrt(4.0 - 3.0 * cos(MagLatMC(i,j)*pi/180)))

              SigmaPPMC(i,j) = spp * sinim
              SigmaLLMC(i,j) = sll / sinim
              SigmaHHMC(i,j) = shh
              SigmaCCMC(i,j) = scc
              if (MagLatMC(i,j) > 0.0) then
                 SigmaPLMC(i,j) = + (shh - scc)
                 SigmaLPMC(i,j) = - (shh + scc)
              else
                 SigmaPLMC(i,j) = - (shh - scc)
                 SigmaLPMC(i,j) = + (shh + scc)
              endif

              KDlmMC(i,j) = -sign(1.0,MagLatMC(i,j)) * kdlm_s * be3
              KDpmMC(i,j) = kdpm_s * be3 * abs(sinim)

              KlmMC(i,j) = -sign(1.0,MagLatMC(i,j)) * klm_s
              KpmMC(i,j) = kpm_s * abs(sinim)

           endif

        enddo
     enddo

  enddo
  
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
  call MPI_AllREDUCE(MagBufferMC, DivJuAltMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaHallMC
  call MPI_AllREDUCE(MagBufferMC, SigmaHallMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPedersenMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPedersenMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = LengthMC
  call MPI_AllREDUCE(MagBufferMC, LengthMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPPMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPPMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaLLMC
  call MPI_AllREDUCE(MagBufferMC, SigmaLLMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaHHMC
  call MPI_AllREDUCE(MagBufferMC, SigmaHHMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaCCMC
  call MPI_AllREDUCE(MagBufferMC, SigmaCCMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaLPMC
  call MPI_AllREDUCE(MagBufferMC, SigmaLPMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = SigmaPLMC
  call MPI_AllREDUCE(MagBufferMC, SigmaPLMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KDlmMC
  call MPI_AllREDUCE(MagBufferMC, KDlmMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KDpmMC
  call MPI_AllREDUCE(MagBufferMC, KDpmMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KlmMC
  call MPI_AllREDUCE(MagBufferMC, KlmMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KpmMC
  call MPI_AllREDUCE(MagBufferMC, KpmMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

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
        dSigmaPLdpMC(i,j) = 0.5*(SigmaPLMC(i+1,j) - SigmaPLMC(i-1,j))/deltapmc(i,j)
        dSigmaPPdpMC(i,j) = 0.5*(SigmaPPMC(i+1,j) - SigmaPPMC(i-1,j))/deltapmc(i,j)
        dKDpmdpMC(i,j) = 0.5*(KDpmMC(i+1,j) - KDpmMC(i-1,j)) / deltapmc(i,j)
        dKpmdpMC(i,j)  = 0.5*(KpmMC(i+1,j)  - KpmMC(i-1,j) ) / deltapmc(i,j)
     enddo
     dSigmaPLdpMC(1,j) = (SigmaPLMC(2,j) - SigmaPLMC(1,j))/deltapmc(1,j)
     dSigmaPLdpMC(nMagLons+1,j) = dSigmaPLdpMC(1,j)
     dSigmaPPdpMC(1,j) = (SigmaPPMC(2,j) - SigmaPPMC(1,j))/deltapmc(1,j)
     dSigmaPPdpMC(nMagLons+1,j) = dSigmaPPdpMC(1,j)
     dKDpmdpMC(1,j) = (KDpmMC(2,j) - KDpmMC(1,j))/deltapmc(1,j)
     dKDpmdpMC(nMagLons+1,j) = dKDpmdpMC(1,j)
     dKpmdpMC(1,j) = (KpmMC(2,j) - KpmMC(1,j))/deltapmc(1,j)
     dKpmdpMC(nMagLons+1,j) = dKpmdpMC(1,j)
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
       sppline, sllline, shhline, sccline, kdpline, kdlline, be3, &
       kpline, klline)

    real, intent(out) :: juline, shline, spline, &
         sppline, sllline, shhline, sccline, kdpline, kdlline, &
         be3, kpline, klline

    integer :: ii, jj

    real :: dip, dec, length, mfac, lfac, gmlt, gmlt2, cdip, sdip, mt

    logical :: IsFound

    juline = 0.0
    shline = 0.0
    spline = 0.0
    length = 0.0

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

    call interpolate_local(HallFieldLine,     mfac, lfac, ii, jj, shline)
    call interpolate_local(PedersenFieldLine, mfac, lfac, ii, jj, spline)
    call interpolate_local(DivJuFieldLine,    mfac, lfac, ii, jj, juline)
    call interpolate_local(LengthFieldLine,   mfac, lfac, ii, jj, length)

    call interpolate_local(SigmaPP, mfac, lfac, ii, jj, sppline)
    call interpolate_local(SigmaLL, mfac, lfac, ii, jj, sllline)
    call interpolate_local(SigmaHH, mfac, lfac, ii, jj, shhline)
    call interpolate_local(SigmaCC, mfac, lfac, ii, jj, sccline)

    call interpolate_local(KDpm, mfac, lfac, ii, jj, kdpline)
    call interpolate_local(KDlm, mfac, lfac, ii, jj, kdlline)
    
    call interpolate_local(Kpm, mfac, lfac, ii, jj, kpline)
    call interpolate_local(Klm, mfac, lfac, ii, jj, klline)

    call interpolate_local(b0_be3(:,:,1,iBlock), mfac, lfac, ii, jj, be3)

    if (sppline > 1000.) write(*,*) "sppline : ", &
         mfac, lfac, ii, jj, sppline, &
         SigmaPP(ii,jj), SigmaPP(ii+1, jj), SigmaPP(jj+1,ii), SigmaPP(ii+1,jj+1)

  end subroutine find_mag_point

  !\
  ! Since this is a contains subroutine, we don't need to pass in
  ! mlatMC or mltMC or iBlock.
  !/

  subroutine find_mag_point_old(juline, shline, spline, length, &
       sppline, sllline, shhline, sccline, kdpline, kdlline)

    real, intent(out) :: juline, shline, spline, &
         sppline, sllline, shhline, sccline, kdpline, kdlline

    integer :: ii, jj

    real :: dip, dec, length, mfac, lfac, gmlt, gmlt2, cdip, sdip, mt

    logical :: IsFound

    juline = 0.0
    shline = 0.0
    spline = 0.0
    length = 0.0

    mltMC = mod(mltMC+24.0,24.0)

    gmlt  = maxval(magloctime_local)
    gmlt2 = minval(magloctime_local)

!    if (gmlt2 < gmlt .and. mltMC < gmlt2) then
!       gmlt  = gmlt - 24.0
!    else
!       if (gmlt2 < gmlt .and. mltMC > gmlt)  gmlt2 = gmlt2 + 24.0
!    endif
!
!    if (mltMC > gmlt) return
!    if (mltMC < gmlt2) return

    if (mlatMC > maxval(MLatitude(0:nLons+1,0:nLats+1,k,iBlock))) return
    if (mlatMC < minval(MLatitude(0:nLons+1,0:nLats+1,k,iBlock))) return

    ii = 0
    do while (ii <= nLons)
       jj = 0
       do while (jj <= nLats-1)

          gmlt  = magloctime_local(ii  ,jj)
          gmlt2 = magloctime_local(ii+1,jj)

          if (gmlt2 < gmlt .and. mltMC < gmlt2) then
             gmlt  = gmlt - 24.0
          else
             if (gmlt2 < gmlt .and. mltMC > gmlt)  gmlt2 = gmlt2 + 24.0
          endif

          IsFound = .false.
          
          if ( mltMC  >= gmlt  .and. &
               mltMC  <=  gmlt2 .and. &
               mlatMC >= MLatitude( ii  ,jj  ,k, iBlock) .and. &
               mlatMC <= MLatitude( ii  ,jj+1,k, iBlock)) then
             if ( mltMC  >= gmlt  .and. &
                  mltMC  <  gmlt2 .and. &
                  mlatMC >= MLatitude( ii  ,jj  ,k, iBlock) .and. &
                  mlatMC <  MLatitude( ii  ,jj+1,k, iBlock)) then
                IsFound = .true.
             else
                if (ii == nLons .and. mltMC ==  gmlt2 .and. &
                     mlatMC >= MLatitude( ii  ,jj  ,k, iBlock) .and. &
                     mlatMC <  MLatitude( ii  ,jj+1,k, iBlock)) &
                     IsFound = .true.
                
                if (jj == nLats .and. &
                     mlatMC == MLatitude( ii  ,jj+1,k, iBlock) .and. &
                     mltMC  >= gmlt  .and. &
                     mltMC  <  gmlt2) &
                     IsFound = .true.

                if (ii == nLons .and. mltMC ==  gmlt2 .and. &
                     jj == nLats .and. &
                     mlatMC ==  MLatitude( ii  ,jj+1,k, iBlock)) &
                     IsFound = .true.

             endif
          endif

!          if (ii == 0 .or. jj == 0) IsFound = .false.

          if (IsFound) then

             mfac = (mltMC - gmlt) / (gmlt2-gmlt)

             lfac = (mlatMC - MLatitude(ii,jj,k,iBlock)) / &
                    (MLatitude(ii,jj+1,k,iBlock) - &
                     MLatitude(ii,jj  ,k,iBlock))

             call interpolate_local(HallFieldLine,     mfac, lfac, ii, jj, shline)
             call interpolate_local(PedersenFieldLine, mfac, lfac, ii, jj, spline)
             call interpolate_local(DivJuFieldLine,    mfac, lfac, ii, jj, juline)
             call interpolate_local(LengthFieldLine,   mfac, lfac, ii, jj, length)

             call interpolate_local(SigmaPP, mfac, lfac, ii, jj, sppline)
             call interpolate_local(SigmaLL, mfac, lfac, ii, jj, sllline)
             call interpolate_local(SigmaHH, mfac, lfac, ii, jj, shhline)
             call interpolate_local(SigmaCC, mfac, lfac, ii, jj, sccline)

             call interpolate_local(KDpm, mfac, lfac, ii, jj, kdpline)
             call interpolate_local(KDlm, mfac, lfac, ii, jj, kdlline)

!!             shline = (1.0-mfac)*(1.0-lfac)*HallFieldLine(ii  ,jj  ) + &
!!                      (    mfac)*(1.0-lfac)*HallFieldLine(ii+1,jj  ) + &
!!                      (1.0-mfac)*(    lfac)*HallFieldLine(ii  ,jj+1) + &
!!                      (    mfac)*(    lfac)*HallFieldLine(ii+1,jj+1)
!!
!!             spline = (1.0-mfac)*(1.0-lfac)*PedersenFieldLine(ii  ,jj  ) + &
!!                      (    mfac)*(1.0-lfac)*PedersenFieldLine(ii+1,jj  ) + &
!!                      (1.0-mfac)*(    lfac)*PedersenFieldLine(ii  ,jj+1) + &
!!                      (    mfac)*(    lfac)*PedersenFieldLine(ii+1,jj+1)
!!
!!             juline = (1.0-mfac)*(1.0-lfac)*DivJuFieldLine(ii  ,jj  ) + &
!!                      (    mfac)*(1.0-lfac)*DivJuFieldLine(ii+1,jj  ) + &
!!                      (1.0-mfac)*(    lfac)*DivJuFieldLine(ii  ,jj+1) + &
!!                      (    mfac)*(    lfac)*DivJuFieldLine(ii+1,jj+1)
!!
!!             length = (1.0-mfac)*(1.0-lfac)*LengthFieldLine(ii  ,jj  ) + &
!!                      (    mfac)*(1.0-lfac)*LengthFieldLine(ii+1,jj  ) + &
!!                      (1.0-mfac)*(    lfac)*LengthFieldLine(ii  ,jj+1) + &
!!                      (    mfac)*(    lfac)*LengthFieldLine(ii+1,jj+1)

!             dip    = (1.0-mfac)*(1.0-lfac)*DipAngle(ii  ,jj  ,k,iBlock) + &
!                      (    mfac)*(1.0-lfac)*DipAngle(ii+1,jj  ,k,iBlock) + &
!                      (1.0-mfac)*(    lfac)*DipAngle(ii  ,jj+1,k,iBlock) + &
!                      (    mfac)*(    lfac)*DipAngle(ii+1,jj+1,k,iBlock)
!
!             dec    = (1.0-mfac)*(1.0-lfac)*DecAngle(ii  ,jj  ,k,iBlock) + &
!                      (    mfac)*(1.0-lfac)*DecAngle(ii+1,jj  ,k,iBlock) + &
!                      (1.0-mfac)*(    lfac)*DecAngle(ii  ,jj+1,k,iBlock) + &
!                      (    mfac)*(    lfac)*DecAngle(ii+1,jj+1,k,iBlock)

!             ! Lets get the length of the field line in the region:
!
!!             mt = mlatMC*pi/180.0
!
!!             cdip = 2.0 * sin(mt) / sqrt(1.0 + 3.0*sin(mt)**2)
!!             sdip = sqrt(1.0 - cdip*cdip)
!
!             sdip = sin(dip)
!             cdip = abs(cos(dip))
!
!             ! Length in altitude
!             length = dAlt_GB(ii,jj,k,iBlock) * sdip
!
!             ! Length in latitude
!             length = length + &
!                  (Latitude(jj+1,iBlock) - Latitude(jj,iBlock)) * &
!                  RadialDistance(k) * cdip

!             ! Length in Longitude
!             length = length * abs(cos(dec)) + &
!                  (Phi(ii+1,jj,k,iBlock) - Phi(ii,jj,k,iBlock)) * &
!                  radial_distance(ii,jj,k) * &
!                  abs(cos(Theta(ii,jj,k,iBlock))) * abs(sin(dec))

!             spline = 1.0
!             shline = 1.0
!             juline = 1.0

             shline = shline
             spline = spline
             juline = juline

!             shline = shline * length
!             spline = spline * length
!             juline = juline * length

             ii = nLons
             jj = nLats

          endif

          jj=jj+1
       enddo
       ii=ii+1
    enddo

  end subroutine find_mag_point_old

  subroutine interpolate_local(VarToInter, mfac, lfac, ii, jj, Output)

    real, intent(in)    :: VarToInter(-1:nLons+2,-1:nLats+2), mfac, lfac
    integer, intent(in) :: ii, jj
    real, intent(out)   :: Output
      
    Output = (1.0-mfac)*(1.0-lfac)*VarToInter(ii  ,jj  ) + &
             (    mfac)*(1.0-lfac)*VarToInter(ii+1,jj  ) + &
             (1.0-mfac)*(    lfac)*VarToInter(ii  ,jj+1) + &
             (    mfac)*(    lfac)*VarToInter(ii+1,jj+1)

  end subroutine interpolate_local

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
  !-------------------------------------------------------------------------

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

