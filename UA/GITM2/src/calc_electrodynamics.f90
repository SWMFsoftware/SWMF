
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

  integer :: i,j,k,bs, iError, iBlock, iDir, iLon, iLat, iAlt

  real :: GeoLat, GeoLon, GeoAlt, xAlt, len, ped, hal
  real :: xmag, ymag, zmag, bmag, signz, alon, alat
  real :: mlatMC, mltMC, jul, shl, spl,length

  real :: e2, dju,dl

  real :: magloctime_local(0:nLons+1,0:nLats+1)

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2) :: &
       e_density, Vi, Ve, MeVen, MeVei, MiVin, VeOe, ViOi, &
       JuDotB

  logical :: IsDone, IsFirstTime = .true.

  call report("UA_calc_electrodynamics",1)
  call start_timing("calc_electrodyn")

  if (IsFirstTime) then

     IsFirstTime = .false.

     allocate(DivJuAltMC(nMagLons+1,nMagLats), &
          SigmaHallMC(nMagLons+1,nMagLats), &
          SigmaPedersenMC(nMagLons+1,nMagLats), &
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

  e2 = Element_Charge * Element_Charge

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

     Sigma_0 = e2 * E_Density / (1.0/MeVen + 1.0/MiVin)

     Sigma_Pedersen = ((1.0/MeVen) * (Ve*Ve/VeOe) + &
          (1.0/MiVin) * (Vi*Vi/ViOi)) * E_Density * e2

     Sigma_Hall = ((1.0/MeVen) * (Ve*e_gyro/VeOe) - &
          (1.0/MiVin) * (Vi*i_gyro/ViOi)) * E_Density * e2

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

     do k=1,nAlts
        do j=-1,nLats+2
           do i=-1,nLons+2
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
           len = 500.0
           xAlt = 1.0
           iAlt = 1

           call get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoAlt/1000.0,&
                alat,alon,xmag,ymag,zmag)
           signz = sign(1.0,zmag)

           do while (.not. IsDone)

              ped = &
                          xAlt  * Sigma_Pedersen(iLon, iLat, iAlt) + &
                   (1.0 - xAlt) * Sigma_Pedersen(iLon, iLat, iAlt+1)
              hal = &
                        xAlt  * Sigma_Hall(iLon, iLat, iAlt) + &
                   (1.0-xAlt) * Sigma_Hall(iLon, iLat, iAlt+1)
              dju = &
                        xAlt  * DivJu(iLon, iLat, iAlt) + &
                   (1.0-xAlt) * DivJu(iLon, iLat, iAlt+1)

              PedersenFieldLine(iLon, iLat) = &
                   PedersenFieldLine(iLon, iLat) + len * ped
              HallFieldLine(iLon, iLat) = &
                   HallFieldLine(iLon, iLat) + len * hal
              DivJuFieldLine(iLon, iLat) = &
                   DivJuFieldLine(iLon, iLat) + len * dju

              LengthFieldLine(iLon, iLat) = &
                   LengthFieldLine(iLon, iLat) + len

              call get_magfield(GeoLat*180.0/pi,GeoLon*180.0/pi,GeoAlt/1000.0,&
                   alat,alon,xmag,ymag,zmag)

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

!           if (iLon == 2) write(*,*) "Divjuline : ",&
!                divjuFieldLine(iLon,iLat), divjualt(iLon,iLat), &
!                HallFieldLine(iLon,iLat), HallConductance(iLon,iLat,iBlock), &
!                LengthFieldLine(iLon,iLat), Altitude(nAlts) - Altitude(1)


        enddo
     enddo

     do k=1,1

        call calc_mltlocal

        do i=1,nMagLons
           do j=1,nMagLats

              mlatMC = MagLatMC(i,j)
              mltMC  = MagLocTimeMC(i,j)

              call find_mag_point(jul, shl, spl, length)

              if (length > 0) then
                 DivJuAltMC(i,j)      = jul
                 SigmaHallMC(i,j)     = shl
                 SigmaPedersenMC(i,j) = spl
                 LengthMC(i,j)        = length
              endif

           enddo
        enddo

     enddo

  enddo
  
  if (iDebugLevel > 2) write(*,*) "Beginning Sum of Electrodynamics"

!  DivJuAltMC(1,:) = (DivJuAltMC(2,:) + DivJuAltMC(nMagLons,:))/2.0
!  LengthMC(1,:)   = (LengthMC(2,:) + LengthMC(nMagLons,:))/2.0
!  SigmaHallMC(1,:) = (SigmaHallMC(2,:) + SigmaHallMC(nMagLons,:))/2.0
!  SigmaPedersenMC(1,:) = (SigmaPedersenMC(2,:) + &
!       SigmaPedersenMC(nMagLons,:))/2.0
  
  DivJuAltMC(nMagLons+1,:) = DivJuAltMC(1,:)
  SigmaHallMC(nMagLons+1,:) = SigmaHallMC(1,:)
  SigmaPedersenMC(nMagLons+1,:) = SigmaPedersenMC(1,:)
  LengthMC(nMagLons+1,:) = LengthMC(1,:)
  
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

  ! There is the possibility that there may be some points missing
  ! near the magnetic pole.  This is due to the grid.  So, let's find
  ! those places and fill them in.

  do i=1,nMagLons+1
     do j=nMagLats-5, nMagLats
        if (LengthMC(i,j) < 0.0) then
           SigmaHallMC(i,j) = sum(SigmaHallMC(:,j-1))/(nMagLons+1)
           SigmaPedersenMC(i,j) = sum(SigmaPedersenMC(:,j-1))/(nMagLons+1)
           LengthMC(i,j) = sum(LengthMC(:,j-1))/(nMagLons+1)
           DivJuAltMC(i,j) = sum(DivJuAltMC(:,j-1))/(nMagLons+1)
        endif
     enddo
     do j=5, 1, -1
        if (LengthMC(i,j) < 0.0) then
           SigmaHallMC(i,j) = sum(SigmaHallMC(:,j+1))/(nMagLons+1)
           SigmaPedersenMC(i,j) = sum(SigmaPedersenMC(:,j+1))/(nMagLons+1)
           LengthMC(i,j) = sum(LengthMC(:,j+1))/(nMagLons+1)
           DivJuAltMC(i,j) = sum(DivJuAltMC(:,j+1))/(nMagLons+1)
        endif
     enddo
  enddo

  do iLat = 1, nMagLats
     if (LengthMC(10,iLat) < 1.0) then
        write(*,*) iLat, &
             LengthMC(10,iLat),MagLatMC(10,iLat)
!        do i = 0, nLats+1
!           write(*,*) i,Mlatitude(1,i,1,1),Mlongitude(1,i,1,1)
!        enddo
        write(*,*) "Problem with electrodynamics. A field-line length"
        write(*,*) "is less than 0.0"
        call stop_gitm("Can't continue")
     endif
  enddo

  where (SigmaHallMC < 0.1) SigmaHallMC = 0.1
  where (SigmaPedersenMC < 0.1) SigmaPedersenMC = 0.1

  call start_timing("calc_electrodyn")

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

  subroutine find_mag_point(juline, shline, spline,length)

    real, intent(out) :: juline, shline, spline

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
          if (IsFound) then

             mfac = (mltMC - gmlt) / (gmlt2-gmlt)

             lfac = (mlatMC - MLatitude(ii,jj,k,iBlock)) / &
                    (MLatitude(ii,jj+1,k,iBlock) - &
                     MLatitude(ii,jj  ,k,iBlock))

             shline = (1.0-mfac)*(1.0-lfac)*HallFieldLine(ii  ,jj  ) + &
                      (    mfac)*(1.0-lfac)*HallFieldLine(ii+1,jj  ) + &
                      (1.0-mfac)*(    lfac)*HallFieldLine(ii  ,jj+1) + &
                      (    mfac)*(    lfac)*HallFieldLine(ii+1,jj+1)

!             shline = (1.0-mfac)*(1.0-lfac)*Sigma_Hall(ii  ,jj  ,k) + &
!                      (    mfac)*(1.0-lfac)*Sigma_Hall(ii+1,jj  ,k) + &
!                      (1.0-mfac)*(    lfac)*Sigma_Hall(ii  ,jj+1,k) + &
!                      (    mfac)*(    lfac)*Sigma_Hall(ii+1,jj+1,k)

             spline = (1.0-mfac)*(1.0-lfac)*PedersenFieldLine(ii  ,jj  ) + &
                      (    mfac)*(1.0-lfac)*PedersenFieldLine(ii+1,jj  ) + &
                      (1.0-mfac)*(    lfac)*PedersenFieldLine(ii  ,jj+1) + &
                      (    mfac)*(    lfac)*PedersenFieldLine(ii+1,jj+1)

!             spline = (1.0-mfac)*(1.0-lfac)*Sigma_Pedersen(ii  ,jj  ,k) + &
!                      (    mfac)*(1.0-lfac)*Sigma_Pedersen(ii+1,jj  ,k) + &
!                      (1.0-mfac)*(    lfac)*Sigma_Pedersen(ii  ,jj+1,k) + &
!                      (    mfac)*(    lfac)*Sigma_Pedersen(ii+1,jj+1,k)

!             juline = (1.0-mfac)*(1.0-lfac)*DivJu(ii  ,jj  ,k) + &
!                      (    mfac)*(1.0-lfac)*DivJu(ii+1,jj  ,k) + &
!                      (1.0-mfac)*(    lfac)*DivJu(ii  ,jj+1,k) + &
!                      (    mfac)*(    lfac)*DivJu(ii+1,jj+1,k)

             juline = (1.0-mfac)*(1.0-lfac)*DivJuFieldLine(ii  ,jj  ) + &
                      (    mfac)*(1.0-lfac)*DivJuFieldLine(ii+1,jj  ) + &
                      (1.0-mfac)*(    lfac)*DivJuFieldLine(ii  ,jj+1) + &
                      (    mfac)*(    lfac)*DivJuFieldLine(ii+1,jj+1)

             length = (1.0-mfac)*(1.0-lfac)*LengthFieldLine(ii  ,jj  ) + &
                      (    mfac)*(1.0-lfac)*LengthFieldLine(ii+1,jj  ) + &
                      (1.0-mfac)*(    lfac)*LengthFieldLine(ii  ,jj+1) + &
                      (    mfac)*(    lfac)*LengthFieldLine(ii+1,jj+1)

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

end subroutine UA_calc_electrodynamics
