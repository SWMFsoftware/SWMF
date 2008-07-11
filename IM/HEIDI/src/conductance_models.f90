
subroutine read_cond_model_coef

  use ModFiles
  use ModAMIE
  use ModConductance

  implicit none

  !  NOAA HPI model:  LONGMX=30,LATMX=20,NDX=10,STEPLAT=2.
  !  Kroelh PEM model:  LONGMX=24,LATMX=40,NDX=9,STEPLAT=1.

  real, dimension(0:lonmdx,0:latmdx,mndx,4) :: hpaf
  real    :: scale
  integer :: iunit, ierr, n, ilat, ilon, maxflx, ilatsave

  character (len=80) :: char80

  scale = 0.001

  iunit = lun_conductivity

  if (index(cond_model,'ihp') > 0) then

     if (DebugLevel > 1) &
          write(*,*) 'Reading IHP background conductance model'
     call merge_str(models_dir, ihp_file)
     
     open(iunit, file = ihp_file, status='old', iostat = ierr)
     if (ierr /= 0) then
        write(6,*) 'Error opening file :',ihp_file
     endif

     longmx = 30
     latmx = 20
     ndx = 10
     steplat = 2.

     do n=1,4
        read(iunit,*) char80
     enddo

     read (iunit,*) char80
     do n=1,ndx
        do ilat=latmx,0,-1
           read (iunit,"(15f7.0)") (halar(ilon,ilat,n),ilon=0,14)
           read (iunit,"(15f7.0)") (halar(ilon,ilat,n),ilon=15,29)
        enddo
     enddo
     halar(longmx,:,:) = halar(0,:,:)
     halar = scale * halar

     read (iunit,*) char80
     do n=1,ndx
        do ilat=latmx,0,-1
           read (iunit,"(15f7.0)") (pedar(ilon,ilat,n),ilon=0,14)
           read (iunit,"(15f7.0)") (pedar(ilon,ilat,n),ilon=15,29)
        enddo
     enddo
     pedar(longmx,:,:) = pedar(0,:,:)
     pedar = scale * pedar

     read (iunit,*) char80
     do n=1,ndx
        do ilat=latmx,0,-1
           read (iunit,"(15f7.0)") (avkar(ilon,ilat,n),ilon=0,14)
           read (iunit,"(15f7.0)") (avkar(ilon,ilat,n),ilon=15,29)
           do ilon=0,29
              if (avkar(ilon,ilat,n) == 2855) &
                   avkar(ilon,ilat,n) = ConductanceBackground(avee_)/scale
           enddo
        enddo
     enddo
     avkar(longmx,:,:) = avkar(0,:,:)
     avkar = scale * avkar

     read (iunit,*) char80
     do n=1,ndx
        do ilat=latmx,0,-1
           read (iunit,"(15f7.0)") (efxar(ilon,ilat,n),ilon=0,14)
           read (iunit,"(15f7.0)") (efxar(ilon,ilat,n),ilon=15,29)
        enddo
!        if (UseExperimentalCode) then
!           do ilon=0,29
!              maxflx = 0
!              do ilat=latmx,0,-1
!                 if (efxar(ilon,ilat,n) > maxflx) then
!                    maxflx = efxar(ilon,ilat,n)
!                    ilatsave = ilat
!                 endif
!              enddo
!              efxar(ilon,ilatsave,n) = efxar(ilon,ilatsave,n)*2
!           enddo
!        endif
     enddo
     efxar(longmx,:,:) = efxar(0,:,:)
     efxar = scale * efxar

     close(iunit)

  else

     if (DebugLevel > 1) &
          write(*,*) 'Reading PEM background conductance model'

     longmx = 24
     latmx = 40
     ndx =  9
     steplat = 1.

     call merge_str(models_dir, pem_file)

     open(iunit, file = pem_file, status='old', iostat = ierr)
     if (ierr /= 0) then
        write(6,*) 'Error opening file :',pem_file
     endif

     do n=1,4
        read(iunit,*) char80
     enddo

     do n=1,ndx
        read (iunit,"(a80)") char80
        do ilat=latmx,0,-1
           read (iunit,"(24f6.0)") (halar(ilon,ilat,n),ilon=0,23)
        enddo
     enddo
     halar(longmx,:,:) = halar(0,:,:)
     halar = scale * halar

     do n=1,ndx
        read (iunit,"(a80)") char80
        do ilat=latmx,0,-1
           read (iunit,"(24f6.0)") (pedar(ilon,ilat,n),ilon=0,23)
        enddo
     enddo
     pedar(longmx,:,:) = pedar(0,:,:)
     pedar = scale * pedar

     do n=1,ndx
        read (iunit,"(a80)") char80
        do ilat=latmx,0,-1
           read (iunit,"(24f6.0)") (avkar(ilon,ilat,n),ilon=0,23)
        enddo
     enddo
     avkar(longmx,:,:) = avkar(0,:,:)
     avkar = scale * avkar

     do n=1,ndx
        read (iunit,"(a80)") char80
        do ilat=latmx,0,-1
           read (iunit,"(24f6.0)") (efxar(ilon,ilat,n),ilon=0,23)
        enddo
     enddo
     efxar(longmx,:,:) = efxar(0,:,:)
     efxar = scale * efxar

     close(iunit)

  endif

  ! get minimum or average long. value at lowest lat.

  do n=1,ndx

     halmin(n) = halar(1,latmx,n)
     pedmin(n) = pedar(1,latmx,n)
     avk50(n) = avkar(1,latmx,n)
     efx50(n) = efxar(1,latmx,n)

     do ilon=2,longmx
        halmin(n) = amin1(halmin(n),halar(ilon,latmx,n))
        pedmin(n) = amin1(pedmin(n),pedar(ilon,latmx,n))
        avk50(n) = avk50(n) + avkar(1,latmx,n)
        efx50(n) = efx50(n) + efxar(1,latmx,n)
     enddo

     avk50(n) = avk50(n) / float(longmx)
     efx50(n) = efx50(n) / float(longmx)

  enddo

  return

end subroutine read_cond_model_coef

subroutine get_auroral_conductance(alatd, amlt, hpi, ped, hal, avkev, eflx)

  Use ModConductance
  Use ModAMIE

  implicit none

  real, intent(in) :: alatd, amlt, hpi
  real, intent(out) :: ped, hal, avkev, eflx

  real    :: dx, tl, x, y
  integer :: i, j, lon

  dx = 24./float(longmx)
  i = hpi + 0.49999
  if (i < 1) i = 1
  if (i > ndx) i = ndx

  y = (90.-abs(alatd))/steplat

  if (y > float(latmx)) then

     ped = pedmin(i)
     hal = halmin(i)
     avkev = avk50(i)
     eflx = efx50(i)

  else

     j = y
     if (j < 0) j = 0
     if (j > latmx-1) j = latmx-1

     y = y - j

     tl = amlt
     if (tl < 0.) tl = tl + 24.

     lon = tl*longmx/24.

     x = tl - 24.*lon/float(longmx)
     x = x / dx

     if (lon > longmx-1) lon = lon - longmx
     if (lon < 0.) lon = lon + longmx

     ped = x*(y*pedar(lon+1,j+1,i) + (1.-y)*pedar(lon+1,j,i)) + &
          (1.-x)*(y*pedar(lon,j+1,i) + (1.-y)*pedar(lon,j,i))
     hal = x*(y*halar(lon+1,j+1,i) + (1.-y)*halar(lon+1,j,i)) + &
          (1.-x)*(y*halar(lon,j+1,i) + (1.-y)*halar(lon,j,i))
     avkev = x*(y*avkar(lon+1,j+1,i) + (1.-y)*avkar(lon+1,j,i)) + &
          (1.-x)*(y*avkar(lon,j+1,i) + (1.-y)*avkar(lon,j,i))
     eflx  = x*(y*efxar(lon+1,j+1,i) + (1.-y)*efxar(lon+1,j,i)) + &
          (1.-x)*(y*efxar(lon,j+1,i) + (1.-y)*efxar(lon,j,i))

     if (ped   < ConductanceBackground(pedersen_)) ped   = ConductanceBackground(pedersen_)
     if (hal   < ConductanceBackground(hall_))     hal   = ConductanceBackground(hall_)
     if (avkev < ConductanceBackground(avee_))     avkev = ConductanceBackground(avee_)
     if (eflx  < ConductanceBackground(eflux_))    eflx  = ConductanceBackground(eflux_)

  endif

  return

end subroutine get_auroral_conductance


subroutine get_old_solar_cond(alatd, amlt, sigma_ped, sigma_hall, sza)

  Use ModAMIE
  use ModConstants

  implicit none

!  Simple routine to calculate the solar zenith angle and assign
!   a Pedersen and Hall conductance (sgp,sgh) due to UV radiation
!  A. Richmond, HAO/NCAR

  real, intent(in)  :: alatd, amlt
  real, intent(out) :: sigma_ped, sigma_hall, sza

  real :: clatrad, sz, cz, cosk, glat, glon
  real :: alonx, fcosk, bbp, bbh
  integer :: ierr

  clatrad = (90.-abs(alatd))/rad
  sz      = sin(clatrad)
  cz      = cos(clatrad)

  cosk=cz*cos_subsol_mag_colat + &
       sz*cos(amlt*pi/12.-pi)*sin_subsol_mag_colat

  !  find alon for mlt
  call mlt2alon (amlt, subsolar_lat, subsolar_lon,   &
       pole_colat, pole_east_lon, alonx)

  !  find glat,glon
  call apxq2g (alatd,alonx,120.,wkap,glat,glon,ierr)

  !  find cos of solar zenith angle
  call cossza (glat,glon,subsolar_lat,subsolar_lon,cosk)

  sza = acos(cosk) * rad

  fcosk=0.03+exp(1.803*tanh(3.833*cosk)+0.5*cosk-2.332)
  fcosk=fcosk+0.03

  bbp=sqrt(1.0-0.99524*sz*sz)*(1.0+0.3*cz*cz)
  bbh=sqrt(1.0-0.01504*(1.0-cz)-0.97986*sz*sz)*(1.0+0.5*cz*cz)

  sigma_ped  = 12.5 * fcosk / bbp
  sigma_hall = 17.0 * fcosk / bbh

  return

end subroutine get_old_solar_cond

subroutine get_new_solar_cond( f107, alatd, rmlt, sigma_ped, sigma_hall, sza)

  use ModAMIE
  use ModConstants

  implicit none

  real, intent(in)  :: f107, alatd, rmlt
  real, intent(out) :: sigma_ped, sigma_hall, sza

!  solpdhl calculates the solar component of the hall and pedersen
!   conductivities in mhos, sigma_hall and sigma_ped.
!  input:
!    f107 = observed daily 10.7 cm flux
!    alatd = apex latitude in degrees
!    colatm = magnetic colatitude in radians.
!    rmlt = magnetic local time (from apex lon) in hours
!    rmlon0 = longitude from midnight in radians.
!
!  internal variables:
!    rmlt = local magnetic time in hours.
!    sza = solar zenith angle in deg
!    szan = solar zenith angle in deg at noon
!    szam = solar zenith angle in deg at midnight
!
!    have an am and pm curve.  this can lead to discontinuities at noon
!  (when sun is low) or at midnight (when sun is high).  discontinuities
!  are eliminated by using am curve only when sun is low and pm curve
!  only when sun is high.
!
!  formulas by h. kroehl (noaa, boulder), 4/89, based on data from
!   chatanika (corrected magnetic latitude = 65.3 deg).
!  formulas revised by a. richmond (ncar, boulder), 6/89, to account
!   for the change in the magnitude of b over the globe.  formulas
!   were normalized for the magnetic latitude of chatanika.
!  implemented by b. emery (ncar, boulder), 6/89.
!***********************************************************************

  real, parameter :: pwr = 2.0/3.0

  real :: szan, szam, csza, cszam, cszan, csza65
  real :: scol, ccol, colatm, usemlt, alonx, glon, glat, alon0, alon12
  real :: sigp65, pmid, sigp00, hcalc, pcalc, sigh65
  real :: bdipol, bnorm, bbp, bbh, bbpn, bbhn, pcalcb, hcalcb

  integer :: ierr

  colatm = (90.-abs(alatd))*dtor
  ccol = cos(colatm)
  scol = sin(colatm)

  csza = ccol*cos_subsol_mag_colat + &
         scol*cos(rmlt*pi/12.-pi)*sin_subsol_mag_colat

  !  find alon for mlt
  call mlt2alon (rmlt,subsolar_lat,subsolar_lon, &
       pole_colat,pole_east_lon,alonx)

  !  find glat,glon
  call apxq2g (alatd,alonx,120.,wkap,glat,glon,ierr)

  !  find cos of solar zenith angle
  call cossza (glat,glon,subsolar_lat,subsolar_lon,csza)
  sza = acos(csza) * rad

  !  find solar zenith angle in deg at noon (rmlon0=pi) and midnight (=0)

  cszan = ccol*cos_subsol_mag_colat + scol*sin_subsol_mag_colat
  call mlt2alon (12.,subsolar_lat,subsolar_lon, &
       pole_colat, pole_east_lon,alon12)
  call apxq2g (alatd,alon12,120.,wkap,glat,glon,ierr)
  call cossza (glat,glon,subsolar_lat,subsolar_lon,cszan)
  szan = acos(cszan) * rad

  cszam = ccol*cos_subsol_mag_colat - scol*sin_subsol_mag_colat
  call mlt2alon (0.,subsolar_lat,subsolar_lon, &
       pole_colat, pole_east_lon,alon0)
  call apxq2g (alatd,alon0,120.,wkap,glat,glon,ierr)

  call cossza (glat,glon,subsolar_lat,subsolar_lon,cszam)
  szam = acos(cszam) * rad

  ! if sza > 65 deg. at noon use am curve for pm too.

  if (szan > 65.0) then

     usemlt = 11.0

     !          check conductance at midnight--if > 0.4 mho
     !                               use pm curve for entire day.

  else

     csza65 = cos(65.0*dtor)
     sigp65 = 0.5 * (csza65**pwr) * (f107**pwr)
     if (szam.le.100.0) then
        pmid = sigp65 - 0.22 * (szam-65.0)
     else
        sigp00 = sigp65 - 0.22 * 35.0
        pmid = sigp00 - 0.13 * (szam-100.0)
     endif
     if (pmid.gt.0.4) then
        usemlt = 13.00
     else
        usemlt = rmlt
     endif

  endif

  !                       sza <= 65 deg. has power relationship.

  if (sza < 65.0) then

     hcalc = 1.8 * csza * sqrt(f107)
     pcalc = 0.5 * (csza**pwr) * (f107**pwr)

  else

     !          sza > 65 deg. has single line segment to min for hall.
     !          sza > 65 deg. has straight line segments to min for ped.

     csza65 = cos(65.0*dtor)
     sigh65 = 1.8 * csza65 * sqrt(f107)
     hcalc = sigh65 - 0.27 * (sza-65.0)

     if (hcalc < 0.8) hcalc = 0.8

     sigp65 = 0.5 * (csza65**pwr) * (f107**pwr)
     if (usemlt.lt.12.0) then

        !                          am has single line segment to minimum.
        pcalc = sigp65 - 0.26 * (sza-65.0)

     else

        !                          pm changes slope at sza = 100 deg.
        if (sza.le.100.0) then

        !                          pm, 65 < sza <= 100.
           pcalc = sigp65 - 0.22 * (sza-65.0)
               
        else

        !                           pm, sza > 100.

           sigp00 = sigp65 - 0.22 * 35.0
           pcalc = sigp00 - 0.13 * (sza-100.0)

        endif

     endif

     if (pcalc < 0.4) pcalc = 0.4

  endif

!
!  add in the effect of the changing magnitude of b and normalize at
!   the magnetic latitude of chatanika (65.3 deg, 24.7 colat).
!
!  rasmussen and schunk (jgr, 1988, p. xxxx) propose a variation of
!   1/b when examining arecibo and chatanika data.
!  dipole field b = beq sqrt(1+3cos(colat)**2) * (re/(re+h))**3
!   where beq=0.282 gauss, re=6371 km, h=110 km.  b(chatanika) =0.499

  bdipol = 0.282 * sqrt(1.+3.*ccol**2) * 0.94994
  pcalcb = pcalc * 0.499 / bdipol
  hcalcb = hcalc * 0.499 / bdipol
  bnorm = 0.499 / bdipol

!  richmond has been using another variation of 1/bbp and 1/bbh where
!   the variation is stronger for hall than for pedersen conductance.
!   this is based on fig. 9. in kamide and matsushita (jgr, 1979,
!   p. 4083) which describes model results.
!   bbp(chatanika) = 1.134, bbh(chatanika) = 1.285

  bbp = sqrt(1.0-0.99524*scol*scol) * (1.0+0.3*ccol*ccol)
  bbh = sqrt(1.0-0.01504*(1.0-ccol)-0.97986*scol*scol) * &
       (1.0+0.5*ccol*ccol)
  sigma_ped  = pcalc * 1.134 / bbp
  sigma_hall = hcalc * 1.285 / bbh
  bbpn = 1.134 / bbp
  bbhn = 1.285 / bbh

end subroutine get_new_solar_cond


subroutine calc_background_conductivity

  use ModAMIE
  use ModIndices
  use ModConductance
  use ModConstants

  implicit none

  integer :: i,j,n
  real    :: sza, alatd

  n = iteration_number

  HPI_Model(n) = 0.0

  do i = 0, ithmx

     alatd = 90.0 - clatd(i)

     do j = 1, lonmx

        call get_auroral_conductance(alatd, mlt(j), Indices(n, hpi_norm_), &
             Model_Auroral_Pedersen(j,i), Model_Auroral_Hall(j,i),       &
             Model_Average_Energy(j,i), Model_Energy_Flux(j,i))

        if (solar_model == 'old') then
           call get_old_solar_cond(alatd*ihsoln, mlt(j), &
                Model_Solar_Pedersen(j,i), Model_Solar_Hall(j,i), sza)
        else
           call get_new_solar_cond( Indices(n,f107_), alatd*ihsoln, mlt(j), &
                Model_Solar_Pedersen(j,i), Model_Solar_Hall(j,i), sza)
        endif

        Model_Total_Pedersen(j,i) = sqrt(Model_Auroral_Pedersen(j,i)**2 + &
                                         Model_Solar_Pedersen(j,i)**2)

        Model_Total_Hall(j,i) = sqrt(Model_Auroral_Hall(j,i)**2 +         &
                                     Model_Solar_Hall(j,i)**2)

        HPI_Model(n) = HPI_Model(n) + Model_Energy_Flux(j,i) * 1.0e-3 * &
             st(i) * Ri * Ri * dlon * dth * 1.0e-9

     enddo

     Model_Average_Energy(0,i)   = Model_Average_Energy(lonmx,i)
     Model_Energy_Flux(0,i)      = Model_Energy_Flux(lonmx,i)
     Model_Auroral_Pedersen(0,i) = Model_Auroral_Pedersen(lonmx,i)
     Model_Auroral_Hall(0,i)     = Model_Auroral_Hall(lonmx,i)
     Model_Solar_Pedersen(0,i)   = Model_Solar_Pedersen(lonmx,i)
     Model_Solar_hall(0,i)       = Model_Solar_Hall(lonmx,i)
     Model_Total_Pedersen(0,i)   = Model_Total_Pedersen(lonmx,i)
     Model_Total_Hall(0,i)       = Model_Total_Hall(lonmx,i)

  enddo

end subroutine calc_background_conductivity

