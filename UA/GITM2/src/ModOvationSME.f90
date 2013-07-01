
module ModOvationSME

  implicit none

  integer, parameter :: nMlts = 96
  integer, parameter :: nMLats = 80

  ! je = energy flux
  ! jn = number flux

  real, dimension(5, nMlts, nMLats) :: &
       beta_je_diff, beta_je_mono, beta_je_wave, beta_je_iono, &
       beta_jn_diff, beta_jn_mono, beta_jn_wave, beta_jn_iono

  real, dimension(nMlts, nMLats) :: &
       NumberFluxDiff, EnergyFluxDiff, &
       NumberFluxMono, EnergyFluxMono, &
       NumberFluxWave, EnergyFluxWave, &
       NumberFluxIons, EnergyFluxIons, &
       Area

  character (len=*), parameter :: cDirectory = "/UA/DataIn/Aurora/"

contains

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine run_ovationsme(GitmStartTime, GitmCurrentTime, iBlock)

    use ModSizeGITM
    use ModGITM
    use ModInputs

    use ModKind, ONLY: Real8_
    use ModIndicesInterfaces

    real(Real8_), intent(in) :: GitmStartTime, GitmCurrentTime
    integer, intent(in) :: iBlock

    logical :: IsFirstTime=.true.

    real :: t, onsettimePast, onsettimeFuture
    real :: DeltaTone, DeltaTtwo, sme

    integer :: iLon, iLat, iMlat, iMlt, iMlat2
    real :: numflux, hps, hpn

    real, parameter :: Q = 1.602e-12 ! eV to ergs

    integer :: iError

    call start_timing("run_ovationsme")

    t = GitmCurrentTime - GitmStartTime

    call get_onsetut(GitmCurrentTime, 1, onsettimePast, iError)
    call get_onsetut(GitmCurrentTime, 2, onsettimeFuture, iError)
    call get_ae(GitmCurrentTime, sme, iError)

    DeltaTone = t - onsettimePast
    DeltaTtwo = onsettimeFuture - t

    call calc_flux(beta_je_diff, sme, DeltaTone, DeltaTtwo, EnergyFluxDiff)
    call calc_flux(beta_je_mono, sme, DeltaTone, DeltaTtwo, EnergyFluxMono)
    call calc_flux(beta_je_wave, sme, DeltaTone, DeltaTtwo, EnergyFluxWave)
    call calc_flux(beta_je_iono, sme, DeltaTone, DeltaTtwo, EnergyFluxIons)

    EnergyFluxDiff = EnergyFluxDiff * Q * &
         beta_je_diff(5,:,:) / (beta_je_iono(5,:,:)+0.01)
    EnergyFluxMono = EnergyFluxMono * Q * &
         beta_je_mono(5,:,:) / (beta_je_iono(5,:,:)+0.01)
    EnergyFluxWave = EnergyFluxWave * Q * &
         beta_je_wave(5,:,:) / (beta_je_iono(5,:,:)+0.01)
    EnergyFluxIons = EnergyFluxIons * Q * &
         beta_je_iono(5,:,:) / (beta_je_iono(5,:,:)+0.01)

    call calc_flux(beta_jn_diff, sme, DeltaTone, DeltaTtwo, NumberFluxDiff)
    call calc_flux(beta_jn_mono, sme, DeltaTone, DeltaTtwo, NumberFluxMono)
    call calc_flux(beta_jn_wave, sme, DeltaTone, DeltaTtwo, NumberFluxWave)
    call calc_flux(beta_jn_iono, sme, DeltaTone, DeltaTtwo, NumberFluxIons)

    NumberFluxDiff = NumberFluxDiff * &
         beta_jn_diff(5,:,:) / (beta_jn_iono(5,:,:)+0.01)
    NumberFluxMono = NumberFluxMono * &
         beta_jn_mono(5,:,:) / (beta_jn_iono(5,:,:)+0.01)
    NumberFluxWave = NumberFluxWave * &
         beta_jn_wave(5,:,:) / (beta_jn_iono(5,:,:)+0.01)
    NumberFluxIons = NumberFluxIons * &
         beta_jn_iono(5,:,:) / (beta_jn_iono(5,:,:)+0.01)

    ElectronEnergyFlux = 0.01
    ElectronAverageEnergy = 0.1

    do iLon = -1, nLons+2
       do iLat = -1, nLats+2

!          write(*,*) iLon, iLat, MLatitude(iLon, iLat, nAlts+1, iBlock), &
!               mod(MLT(iLon, iLat, nAlts+1)+24.0,24.0)

          if (abs(MLatitude(iLon, iLat, nAlts+1, iBlock)) > 50.0) then

             iMlat = floor(abs(MLatitude(iLon, iLat, nAlts+1, iBlock)) - 50.0)*2
!             iMlat = min(max(iMlat,1),nMLats/2)
!             if (MLatitude(iLon, iLat, nAlts+1, iBlock) > 0.0) &
!                  iMlat = iMlat + nMlats/2
             iMlt = mod(floor(mod(MLT(iLon, iLat, nAlts+1)+24.0,24.0)*4),nMlts)
             if (iMlt == 0) iMlt = nMlts

!             if (iMlat > nMlats/2) then
!                iMlat2 = iMlat-nMlats/2
!             else
!                iMlat2 = iMlat+nMlats/2
!             endif
!
!             iMlat2 = min(max(iMlat2,1),nMlats)

             ! Diffuse Energy Flux

!             if (EnergyFluxDiff(iMlt, iMlat)==0) then
!
!                ! Add North and South together
!                ElectronEnergyFlux(iLon, iLat) = &
!                     EnergyFluxDiff(iMlt, iMlat) + EnergyFluxDiff(iMlt, iMlat2)
!
!                ! If there are values in both hemisphere, then divide by 2
!                if ( EnergyFluxDiff(iMlt, iMlat) * &
!                     EnergyFluxDiff(iMlt, iMlat2) /= 0) &
!                     ElectronEnergyFlux(iLon, iLat) = &
!                     ElectronEnergyFlux(iLon, iLat)/2.0
!             else
                ElectronEnergyFlux(iLon, iLat) = &
                     EnergyFluxDiff(iMlt, iMlat)
!             endif


             ! Diffuse Number Flux

!             if (NumberFluxDiff(iMlt, iMlat)==0) then
!
!                ! Add North and South together
!                numflux = &
!                     NumberFluxDiff(iMlt, iMlat) + NumberFluxDiff(iMlt, iMlat2)
!
!                ! If there are values in both hemisphere, then divide by 2
!                if ( NumberFluxDiff(iMlt, iMlat) * &
!                     NumberFluxDiff(iMlt, iMlat2) /= 0) &
!                     numflux = numflux/2
!
!             else
                numflux = NumberFluxDiff(iMlt, iMlat)
!             endif

             if (numflux /= 0) then
                ElectronAverageEnergy(iLon,iLat) = &
                     ElectronEnergyFlux(iLon, iLat)/numflux * &
                     6.242e11 / 1000.0 ! ergs -> keV
             endif

          endif

!          write(*,*) iLon, iLat, &
!               ElectronEnergyFlux(iLon, iLat), ElectronAverageEnergy(iLon,iLat)

       enddo

    enddo

    call end_timing("run_ovationsme")

  end subroutine run_ovationsme

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine calc_flux(beta, sme, delta_t_one, delta_t_two, flux)

    real, intent(in)  :: beta(5, nMlts, nMLats)
    real, intent(in)  :: sme, delta_t_one, delta_t_two
    real, intent(out) :: flux(nMlts, nMLats)

    flux = (beta(1,:,:) + &
            beta(2,:,:) * sme + &
            beta(3,:,:) * delta_t_one + &
            beta(4,:,:) * delta_t_two)

  end subroutine calc_flux

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine read_ovationsm_files

    call read_single_file(&
         "Diff_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_je_diff)
    call read_single_file(&
         "Mono_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_je_mono)
    call read_single_file(&
         "Wave_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_je_wave)
    call read_single_file(&
         "Iono_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_je_iono)

    call read_single_file(&
         "Diff_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_jn_diff)
    call read_single_file(&
         "Mono_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_jn_mono)
    call read_single_file(&
         "Wave_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_jn_wave)
    call read_single_file(&
         "Iono_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt", &
         beta_jn_iono)

  end subroutine read_ovationsm_files

  ! ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------

  subroutine read_single_file(infile, VarToRead)

    character (len=*), intent(in) :: infile
    real, intent(out) :: VarToRead(5, nMlts, nMLats)

    integer :: iInputUnit_ = 31

    character (len=150) :: cFile
    logical :: IsThere

    integer :: i,j, npnts, iError
    real    :: mlt_set, mlat_set, B0, Bsme,Bt1,Bt2,Rsq_e

    cFile = "."//cDirectory//infile

    inquire(file=cFile,EXIST=IsThere)
    if (.not.IsThere) then
       write(*,*) cFile//" cannot be found by read_ovationsm_files"
!         call stop_gitm(cFile//" cannot be found by read_ovationsm_files")
    endif

    iError = 0

    open(iInputUnit_,file=cFile,status="old", iostat=iError)

    do while (iError == 0)
       
       read(iInputUnit_,*,iostat=iError) i,j,npnts, mlt_set, mlat_set, B0
       read(iInputUnit_,*,iostat=iError) Bsme,Bt1,Bt2,Rsq_e

       if (iError == 0) then

          if (npnts < 160) then
             B0 = 0.0
             Bsme = 0.0
             Bt1 = 0.0
             Bt2 = 0.0
          endif

          VarToRead(1,i+1,j+1) = B0 
          VarToRead(2,i+1,j+1) = Bsme
          VarToRead(3,i+1,j+1) = Bt1
          VarToRead(4,i+1,j+1) = Bt2
          VarToRead(5,i+1,j+1) = npnts

       endif
       
    enddo

    close(iInputUnit_)

  end subroutine read_single_file

end module ModOvationSME
