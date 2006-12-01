
!\
! --------------------------------------------------------------------
! --------------------------------------------------------------------
!/

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
! --------------------------------------------------------------------
! --------------------------------------------------------------------
!/

subroutine UA_calc_electrodynamics(UAi_nMLTs, UAi_nLats)

  use ModGITM
  use ModInputs, only: iDebugLevel
  use ModElectrodynamics
  use ModMPI

  implicit none

  integer, intent(out) :: UAi_nMLTs, UAi_nLats

  integer :: i,j,k,bs, iError, iBlock, iDir, iLon, iLat, iAlt, ip, im

  real :: GeoLat, GeoLon, GeoAlt, xAlt, len, ped, hal
  real :: sp_d1d1_d, sp_d2d2_d, sp_d1d2_d, sh
  real :: xmag, ymag, zmag, bmag, signz, alon, alat
  real :: mlatMC, mltMC, jul, shl, spl, length, kdmp_s, kdml_s
  real :: sinIm, spp, sll, shh, scc, sccline, sppline, sllline, shhline

  real :: q2, dju,dl, cD

  real :: magloctime_local(0:nLons+1,0:nLats+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       JuDotB, sigmap_d1d1_d, sigmap_d2d2_d, sigmap_d1d2_d, sigmah, &
       ue1, ue2, kmp, kml

  logical :: IsDone, IsFirstTime = .true.

  if (DipoleStrength == 0) return

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
          KDmpMC(nMagLons+1,nMagLats), &
          KDmlMC(nMagLons+1,nMagLats), &
          MagBufferMC(nMagLons+1,nMagLats), &
          LengthMC(nMagLons+1,nMagLats), &
          MagLocTimeMC(nMagLons+1,nMagLats), &
          MagLatMC(nMagLons+1,nMagLats), stat = iError)
     if (iError /= 0) then
        call CON_stop("Error allocating array DivJuAltMC")
     endif

     do i=1,nMagLons+1
        do j=1,nMagLats
           MagLocTimeMC(i,j) = 24.0 * float(i-1) / float(nMagLons)
           MagLatMC(i,j)     = 180.0 * (float(j)-0.5) / float(nMagLats) - 90.0
           ! we could also us the stretch_grid function....
        enddo
     enddo

  endif

  !\
  ! Magnetic grid is defined as:
  ! MLT is in hours = 0 - 24
  ! Lat is in degrees = -90 - 90
  !/

  DivJuAltMC      = -1.0e32
  SigmaHallMC     = 0.0
  SigmaPedersenMC = 0.0
  LengthMC        = -1.0e32

  UAi_nLats = nMagLats
  UAi_nMlts = nMagLons+1

  q2 = Element_Charge * Element_Charge

  do iBlock = 1, nBlocks

     call calc_physics(iBlock)
     call calc_rates(iBlock)

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

     do k=1,nAlts
        PedersenConductance(:,:,iBlock) = PedersenConductance(:,:,iBlock) + &
             Sigma_Pedersen(:,:,k)*dAlt(k)
        HallConductance(:,:,iBlock)     = HallConductance(:,:,iBlock)     + &
             Sigma_Hall(:,:,k)    *dAlt(k)
     enddo

     ! We need to calculate the rotation matrix.  
     ! I am simply using the formula on
     ! page 89 of Kelley, but this really only works for a magnetic field 
     ! with east-west (geographic) component.  Since this is not true
     ! for us, we need to come up with a more complicated formulation.

     call report("Starting Conductances",2)
     do k=1,nAlts
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

              ! These are from eqns 5.19 and 5.20
              kmp(i,j,k) = ue2(i,j,k) * sigmap_d1d1_d(i,j,k) + &
                   (sigmah(i,j,k) - sigmap_d1d2_d(i,j,k)) * ue1(i,j,k)

              kml(i,j,k) = (sigmah(i,j,k) + sigmap_d1d2_d(i,j,k)) * ue2(i,j,k) - &
                   ue1(i,j,k) * sigmap_d1d1_d(i,j,k)

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
     ! magnetic field line.  So, let us first subtract out the 
     ! component of Ju that is along the field line.


     JuDotB = 0.0
     do iDir = 1, 3
        JuDotB = JuDotB + &
             Ju(:,:,:,iDir) * B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)
     enddo

     do iDir = 1, 3
        Ju(:,:,:,iDir) = Ju(:,:,:,iDir) - &
             JuDotB * B0(:,:,:,iDir,iBlock)/B0(:,:,:,iMag_,iBlock)
     enddo

     call divergence(Ju, DivJu, iBlock)

     DivJuAlt = 0.0
     do k=1,nAlts
        DivJuAlt(:,:) = DivJuAlt(:,:) + DivJu(:,:,k)*dAlt(k)
     enddo

     PedersenFieldLine = 0.0
     HallFieldLine     = 0.0
     LengthFieldLine   = 0.0
     DivJuFieldLine    = 0.0

     SigmaPP = 0.0
     SigmaLL = 0.0
     SigmaHH = 0.0
     SigmaCC = 0.0

     KDmp = 0.0
     KDml = 0.0

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

           GeoAlt = Altitude(1)
           IsDone = .false.
           len = 100.0
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

              kdmp_s     = &
                          xAlt  * kmp(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * kmp(iLon, iLat, iAlt+1)

              kdml_s     = &
                          xAlt  * kml(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * kml(iLon, iLat, iAlt+1)


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
              SigmaLL(iLon, iLat) = SigmaLL(iLon, iLat) + len * sp_d2d2_d
              SigmaHH(iLon, iLat) = SigmaHH(iLon, iLat) + len * sh
              SigmaCC(iLon, iLat) = SigmaCC(iLon, iLat) + len * sp_d1d2_d

              KDmp(iLon, iLat) = KDmp(iLon, iLat) + len * kdmp_s
              KDml(iLon, iLat) = KDml(iLon, iLat) + len * kdml_s

              PedersenFieldLine(iLon, iLat) = &
                   PedersenFieldLine(iLon, iLat) + len * ped
              HallFieldLine(iLon, iLat) = &
                   HallFieldLine(iLon, iLat) + len * hal
              DivJuFieldLine(iLon, iLat) = &
                   DivJuFieldLine(iLon, iLat) + len * dju

              LengthFieldLine(iLon, iLat) = &
                   LengthFieldLine(iLon, iLat) + len

              CALL get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoALT/1000.0, &
                   XMAG,YMAG,ZMAG)

              if (sign(1.0,zmag)*signz < 0) then
                 IsDone = .true.
              else
                 bmag = sqrt(xmag*xmag + ymag*ymag + zmag*zmag)
                 GeoAlt = GeoAlt + abs(zmag)/bmag * len
                 if (GeoAlt > Altitude(nAlts)) then
                    IsDone = .true.
                 else
                    if (GeoAlt > Altitude(iAlt+1)) iAlt = iAlt+1
                    xAlt = (GeoAlt - Altitude(iAlt)) / &
                         (Altitude(iAlt+1) - Altitude(iAlt))
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

     do k=1,1

        call calc_mltlocal

        do i=1,nMagLons
           do j=1,nMagLats

              mlatMC = MagLatMC(i,j)
              mltMC  = MagLocTimeMC(i,j)

              call find_mag_point(jul, shl, spl, length, spp, sll, shh, scc, &
              kdmp_s, kdml_s)

!              write(*,*) "find_mag_points : ", i,j,mLatMC, mltMC, length
              if (length > 0) then
                 DivJuAltMC(i,j)      = jul
                 SigmaHallMC(i,j)     = shl
                 SigmaPedersenMC(i,j) = spl
                 LengthMC(i,j)        = length

                 sinim = abs(2.0 * sin(MagLatMC(i,j)*pi/180) * &
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

                 KDmlMC(i,j) = kdml_s
                 KDmpMC(i,j) = kdmp_s

              endif

           enddo
        enddo

     enddo

  enddo
  
  if (iDebugLevel > 2) write(*,*) "Beginning Sum of Electrodynamics"

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

  KDmlMC(nMagLons+1,:) = KDmlMC(1,:)
  KDmpMC(nMagLons+1,:) = KDmpMC(1,:)
  
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

  MagBufferMC = KDmlMC
  call MPI_AllREDUCE(MagBufferMC, KDmlMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  MagBufferMC = KDmpMC
  call MPI_AllREDUCE(MagBufferMC, KDmpMC,  &
       bs, MPI_REAL, MPI_MAX, iCommGITM, iError)

  ! We may have missed some points in MLT, so let's check

  do j=3, nMagLats-3

     i = 1
     if (LengthMC(i,j) < 0.0) then
        ip = i+1
        do while (ip < nMagLons .and. LengthMC(ip,j) < 0.0) 
           ip = ip+1
        enddo
        im = nMagLons+1
        do while (im > 2 .and. LengthMC(im,j) < 0.0) 
           im = im-1
        enddo

        if (ip >= nMagLons .or. im <= 2 .or. &
             LengthMC(ip,j) < 0.0 .or. LengthMC(im,j) < 0.0) then
           if (iDebugLevel > -1) then
              write(*,*) "j=5,nMagLats-5 loop"
              write(*,*) "Problem with electrodynamics. A field-line length"
              write(*,*) "is less than 0.0. ip, im, j, LengthMC(i,j) :",&
                   ip,im,j,LengthMC(i,j)
           endif
           call stop_gitm("Can't continue")
        else

           ! These are the first MLTs (i.e. MLT = 0)
           SigmaHallMC(i,j) = (SigmaHallMC(ip,j)+SigmaHallMC(im,j))/2.0
           SigmaPedersenMC(i,j) = &
                (SigmaPedersenMC(ip,j)+SigmaPedersenMC(im,j))/2.0
           LengthMC(i,j) = (LengthMC(ip,j)+LengthMC(im,j))/2.0
           DivJuAltMC(i,j) = (DivJuAltMC(ip,j)+DivJuAltMC(im,j))/2.0

           SigmaPPMC(i,j) = (SigmaPPMC(ip,j)+SigmaPPMC(im,j))/2.0
           SigmaLLMC(i,j) = (SigmaLLMC(ip,j)+SigmaLLMC(im,j))/2.0
           SigmaHHMC(i,j) = (SigmaHHMC(ip,j)+SigmaHHMC(im,j))/2.0
           SigmaCCMC(i,j) = (SigmaCCMC(ip,j)+SigmaCCMC(im,j))/2.0
           SigmaPLMC(i,j) = (SigmaPLMC(ip,j)+SigmaPLMC(im,j))/2.0
           SigmaLPMC(i,j) = (SigmaLPMC(ip,j)+SigmaLPMC(im,j))/2.0

           KDmlMC(i,j) = (KDmlMC(ip,j)+KDmlMC(im,j))/2.0
           KDmpMC(i,j) = (KDmpMC(ip,j)+KDmpMC(im,j))/2.0

           ! These are the Lats MLTs (i.e. MLT = 24 = 0)
           SigmaHallMC(nMagLons,j)     = SigmaHallMC(i,j)
           SigmaPedersenMC(nMagLons,j) = SigmaPedersenMC(i,j)
           LengthMC(nMagLons,j)        = LengthMC(i,j)
           DivJuAltMC(nMagLons,j)      = DivJuAltMC(i,j)
           SigmaPPMC(nMagLons,j)       = SigmaPPMC(i,j)
           SigmaLLMC(nMagLons,j)       = SigmaLLMC(i,j)
           SigmaHHMC(nMagLons,j)       = SigmaHHMC(i,j)
           SigmaCCMC(nMagLons,j)       = SigmaCCMC(i,j)
           SigmaLPMC(nMagLons,j)       = SigmaLPMC(i,j)
           SigmaPLMC(nMagLons,j)       = SigmaPLMC(i,j)

           KDmlMC(nMagLons,j)       = KDmlMC(i,j)
           KDmpMC(nMagLons,j)       = KDmpMC(i,j)

           if (iDebugLevel > 2) then
              write(*,*) "Bad Field line length found in calc_electrodynamics."
              write(*,*) "Correcting.  This is because of a bad MLT.", i,j
           endif

        endif
     endif

     do i=2,nMagLons+1
        if (LengthMC(i,j) < 0.0) then

           ! Set im (i-minus) to be the last cell, since we are moving from
           ! minus to plus.  We don't need to search backwards since we 
           ! have filled all of those.
           im = i-1

           ! Set ip (i-plus) to be the next cell, and search forward to
           ! find the next point in which there is a good value
           ip = i+1
           if (ip >= nMagLons+1) ip = 1
           do while (ip <= nMagLons .and. LengthMC(ip,j) < 0.0) 
              ip = ip+1
           enddo

           ! We either found a good point or didn't.... Let's check...
           if (LengthMC(ip,j) < 0.0 .or. LengthMC(im,j) < 0.0) then
              if (iDebugLevel > -1) then
                 write(*,*) "i=2,nMagLons+1 loop"
                 write(*,*) "Problem with electrodynamics. A field-line length"
                 write(*,*) "is less than 0.0. ip, im, j, LengthMC(i,j) :",&
                      im,i,ip,j,nMagLons, &
                      LengthMC(i,j), LengthMC(im,j), LengthMC(ip,j)
              endif
              call stop_gitm("Can't continue")
           else

              ! if we did find a good point, average to get a good value.

              SigmaHallMC(i,j) = (SigmaHallMC(ip,j)+SigmaHallMC(im,j))/2.0
              SigmaPedersenMC(i,j) = &
                   (SigmaPedersenMC(ip,j)+SigmaPedersenMC(im,j))/2.0
              LengthMC(i,j) = (LengthMC(ip,j)+LengthMC(im,j))/2.0
              DivJuAltMC(i,j) = (DivJuAltMC(ip,j)+DivJuAltMC(im,j))/2.0

              SigmaPPMC(i,j) = (SigmaPPMC(ip,j)+SigmaPPMC(im,j))/2.0
              SigmaLLMC(i,j) = (SigmaLLMC(ip,j)+SigmaLLMC(im,j))/2.0
              SigmaHHMC(i,j) = (SigmaHHMC(ip,j)+SigmaHHMC(im,j))/2.0
              SigmaCCMC(i,j) = (SigmaCCMC(ip,j)+SigmaCCMC(im,j))/2.0
              SigmaPLMC(i,j) = (SigmaPLMC(ip,j)+SigmaPLMC(im,j))/2.0
              SigmaLPMC(i,j) = (SigmaLPMC(ip,j)+SigmaLPMC(im,j))/2.0

              if (iDebugLevel > 2) then
                 write(*,*) &
                      "Bad Field line length found in calc_electrodynamics."
                 write(*,*) "Correcting.  This is because of a bad MLT.", i,j
              endif
           endif
        endif
     enddo
  enddo

  ! There is the possibility that there may be some points missing
  ! near the magnetic pole.  This is due to the grid.  So, let's find
  ! those places and fill them in.

  do j=nMagLats-5, nMagLats
     do i=1,nMagLons+1

        if (LengthMC(i,j) < 0.0) then

           SigmaHallMC(i,j) = sum(SigmaHallMC(:,j-1))/(nMagLons+1)
           SigmaPedersenMC(i,j) = sum(SigmaPedersenMC(:,j-1))/(nMagLons+1)
           LengthMC(i,j) = sum(LengthMC(:,j-1))/(nMagLons+1)
           DivJuAltMC(i,j) = sum(DivJuAltMC(:,j-1))/(nMagLons+1)

           SigmaPPMC(i,j) = sum(SigmaPPMC(:,j-1))/(nMagLons+1)
           SigmaLLMC(i,j) = sum(SigmaLLMC(:,j-1))/(nMagLons+1)
           SigmaHHMC(i,j) = sum(SigmaHHMC(:,j-1))/(nMagLons+1)
           SigmaCCMC(i,j) = sum(SigmaCCMC(:,j-1))/(nMagLons+1)
           SigmaPLMC(i,j) = sum(SigmaPLMC(:,j-1))/(nMagLons+1)
           SigmaLPMC(i,j) = sum(SigmaLPMC(:,j-1))/(nMagLons+1)

           if (iDebugLevel > 2) then
              write(*,*) "Bad Field line length found in calc_electrodynamics."
              write(*,*) "Correcting.  This is because of a bad Latitude.", i,j
              write(*,*) LengthMC(i,j)
           endif

        endif
     enddo
  enddo
  do j=5, 1, -1
     do i=1,nMagLons+1
        if (LengthMC(i,j) < 0.0) then

           SigmaHallMC(i,j) = sum(SigmaHallMC(:,j+1))/(nMagLons+1)
           SigmaPedersenMC(i,j) = sum(SigmaPedersenMC(:,j+1))/(nMagLons+1)
           LengthMC(i,j) = sum(LengthMC(:,j+1))/(nMagLons+1)
           DivJuAltMC(i,j) = sum(DivJuAltMC(:,j+1))/(nMagLons+1)


           SigmaPPMC(i,j) = sum(SigmaPPMC(:,j+1))/(nMagLons+1)
           SigmaLLMC(i,j) = sum(SigmaLLMC(:,j+1))/(nMagLons+1)
           SigmaHHMC(i,j) = sum(SigmaHHMC(:,j+1))/(nMagLons+1)
           SigmaCCMC(i,j) = sum(SigmaCCMC(:,j+1))/(nMagLons+1)
           SigmaPLMC(i,j) = sum(SigmaPLMC(:,j+1))/(nMagLons+1)
           SigmaLPMC(i,j) = sum(SigmaLPMC(:,j+1))/(nMagLons+1)

           if (iDebugLevel > 2) then
              write(*,*) "Bad Field line length found in calc_electrodynamics."
              write(*,*) "Correcting.  This is because of a bad Latitude.", i,j
           endif

        endif
     enddo
  enddo

  do iLat = 1, nMagLats
     if (LengthMC(10,iLat) < 0.0) then
        write(*,*) "iLat=1,nMagLats loop"
        write(*,*) iLat, &
             LengthMC(10,iLat),MagLatMC(10,iLat)
        write(*,*) "Problem with electrodynamics. A field-line length"
        write(*,*) "is less than 0.0"
        call stop_gitm("Can't continue")
     endif
  enddo

  where (SigmaHallMC < 0.1) SigmaHallMC = 0.1
  where (SigmaPedersenMC < 0.1) SigmaPedersenMC = 0.1

  call end_timing("calc_electrodyn")

contains

  !\
  ! ------------------------------------------------------------------
  !/

  subroutine calc_mltlocal

    magloctime_local = mod(MLT(0:nLons+1,0:nLats+1,k)+24.0,24.0)
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

             call interpolate_local(KDmp, mfac, lfac, ii, jj, kdpline)
             call interpolate_local(KDml, mfac, lfac, ii, jj, kdlline)

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
!             length = dAlt(k) * sdip
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

  end subroutine find_mag_point

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
