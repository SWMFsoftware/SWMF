Module ModAurora
  implicit none
  save
  
  real :: Theta0,dtheta=0.174532925199 !10 degrees
  real,parameter :: E0=30.0e-3 !ergs/cm**2
  real :: Emax

  !\
  ! Stuff for auroral energy deposition and ionization
  !/
  
  real, dimension(:), allocatable :: &
       ED_grid, ED_Energies, ED_Flux, ED_Ion, ED_Heating
  integer :: ED_N_Energies, ED_N_Alts
  real,    allocatable :: ED_Interpolation_Weight(:)
  integer, allocatable :: ED_Interpolation_Index (:)
  real,    allocatable :: AuroralIonRateO_C(:),HeatingRate_C(:)


  integer,parameter :: nSpecies = 3
contains
  !============================================================================
  
  real function get_etop(GmLat,GmLon)
    !returns topside electron energy input
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    real :: eTopMax
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    eTopMax = 10.0*E0
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_etop = &
            (1 - 2.0*DeltaPhi/cPi)*(eTopMax-(eTopMax-E0)*DeltaTheta/dtheta)
       
    else
       get_etop = E0
    endif
    
  end function get_etop
  
  !============================================================================
  !============================================================================
  
  real function get_eflux(GmLat,GmLon)
    !returns topside electron energy input
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_eflux = (1 - 2.0*DeltaPhi/cPi)*(Emax-(Emax-E0)*DeltaTheta/dtheta)
       
    else
       get_eflux = E0
    endif
    
  end function get_eflux
  
  !============================================================================
  !============================================================================
  
  real function get_eAverageE(GmLat,GmLon)
    !returns average precipitating energy
    use ModNumConst, ONLY:cPi,cDegToRad
    
    real,intent(in) :: GmLat,GmLon
    real :: Theta,Phi,DeltaTheta,DeltaPhi
    real,parameter:: AveE0=.2, AveEmax=.75
    !--------------------------------------------------------------------------
    Theta=cPi/2-GmLat*cDegToRad
    Phi = GmLon*cDegToRad
    
    ! If point is inside aurora return auroral etop value, else return default
    
    if(Theta > Theta0-dtheta/2 .and. Theta < Theta0+dtheta/2 &
         .and. phi > 0.5*cPi .and. phi < 1.5*cPi) then
       DeltaTheta = abs(Theta-Theta0)
       DeltaPhi = abs(Phi - cPi)
       get_eAverageE = &
            (1 - 2.0*DeltaPhi/cPi)*(AveEmax-(AveEmax-AveE0)*DeltaTheta/dtheta)
       
    else
       get_eAverageE = AveE0
    endif
    
  end function get_eAverageE
  
  !============================================================================
  
  subroutine set_Emax(Ap)
    real,intent(in) :: Ap
    !--------------------------------------------------------------------------
    
    Emax=1.0+7.0*Ap/207.0
  end subroutine set_Emax
  
  !============================================================================
  
  subroutine set_theta0(nPhi,nTheta,uExBphi_C,uExBtheta_C,Theta_G)
    use ModPWOM, ONLY: DoMoveLine
    integer,intent(in) :: nPhi,nTheta
    real   ,intent(in) :: uExBphi_C(nPhi,nTheta),uExBtheta_C(nPhi,nTheta)
    real   ,intent(in) :: Theta_G(0:nPhi+1,0:nTheta+1)
    real,parameter :: AlphaMax = 0.5            !30 degrees  
    real,parameter :: ThetaMin = 0.174532925199 !10 degrees  
    integer:: iPhi,iTheta
    real :: Alpha
    
    !--------------------------------------------------------------------------
    iPhi = nPhi/2
    do iTheta=nTheta,1,-1
       Theta0=Theta_G(iPhi,iTheta)
       if (.not. DoMoveLine) return
       Alpha = atan(uExBtheta_C(iPhi,iTheta)/uExBphi_C(iPhi,iTheta))
       
       if (Alpha >= AlphaMax .or. Theta0 <= ThetaMin ) then
          return
       endif
    enddo
  end subroutine set_theta0
  
  !============================================================================
  subroutine set_aurora
    use ModNumConst,ONLY:cPi,cHalfPi,cRadToDeg
    Use ModPWOM,    ONLY: ElectronAverageEnergy_C,ElectronEnergyFlux_C,&
                          UseAurora,nLine, PhiLine_I,ThetaLine_I,Time
    use ModPwTime,  ONLY: CurrentTime,StartTime
    use ModIndicesInterfaces, only: get_HPI

    real, allocatable ::MLT_II(:,:), MLatitude_II(:,:)
    integer :: iLine, iLat, iLon, nLat, nLon, iError
    real    :: temp
    !--------------------------------------------------------------------------
    
    ! Allocated and set reduced grid of fieldline foot points
    if (.not.allocated(MLT_II)) &
         allocate(MLT_II(nLine,nLine),MLatitude_II(nLine,nLine))
    if (.not.allocated(ElectronAverageEnergy_C)) &
         allocate(ElectronAverageEnergy_C(nLine,nLine),ElectronEnergyFlux_C(nLine,nLine))
    
    nLat = nLine
    nLon = nLine
    do iLine = 1,nLine 
       MLT_II(iLine,1:nLat)=mod(PhiLine_I(iLine)*12.0/cPi+12.0,24.0)
       MLatitude_II(1:nLon,iLine)=(cHalfPi-ThetaLine_I(iLine))*cRadToDeg
    enddo

    call IO_SetnMLTs(nLon)
    call IO_SetnLats(nLat)
    CurrentTime=StartTime+Time
    call IO_SetTime(CurrentTime)
    
    call IO_SetGrid(                    &
         MLT_II(1:nLon,1:nLat), &
         MLatitude_II(1:nLon,1:nLat), iError)
    if (iError /= 0) then
       write(*,*) 'Error in set_aurora (IO_SetGrid)'
    endif
    
    ! get HPI index
    call get_HPI(CurrentTime, temp, iError)
    call IO_SetHPI(temp)
    if (iError /= 0) then
       call con_stop("PW_Error in get_hpi called from set_aurora.f90")
    endif
    
    
    ! get the average electron energy and make sure that it falls in 
    ! a reasonable range
    call IO_GetAveE(ElectronAverageEnergy_C, iError)

    if (iError /= 0) then
       write(*,*) "Error in set_aurora (IO_GetAveE):"
       write(*,*) iError
       call con_stop("Stopping in set_aurora")
    endif
    
    do iLat=1,nLat
       do iLon=1,nLon
          if (ElectronAverageEnergy_C(iLon,iLat) < 0.0) then
             ElectronAverageEnergy_C(iLon,iLat) = 0.1
             write(*,*) "i,j Negative : ",iLon,iLat,&
                  ElectronAverageEnergy_C(iLon,iLat)
          endif
          if (ElectronAverageEnergy_C(iLon,iLat) > 100.0) then
             write(*,*) "i,j Positive : ",iLon,iLat,&
                  ElectronAverageEnergy_C(iLon,iLat)
             ElectronAverageEnergy_C(iLon,iLat) = 0.1
          endif
       enddo
    enddo
    
    ! get the electron energy flux
    call IO_GetEFlux(ElectronEnergyFlux_C, iError)
    if (iError /= 0) then
       write(*,*) "Error in set_aurora (IO_GetEFlux):"
       write(*,*) iError
       call con_stop("Stopping in set_aurora")
    endif
    
  end subroutine set_aurora
  !============================================================================
  subroutine set_auroral_rates
    use ModCommonVariables, ONLY: nAlts => nDim
    use ModPWOM, ONLY: nLine,ElectronAverageEnergy_C,ElectronEnergyFlux_C
    use ModNumConst, ONLY:cPi

    real :: alat, hpi, ped, hal, av_kev, eflx_ergs, a,b, maxi
    real :: Factor,temp_ED, avee, eflux, p
    integer :: i, j, k, n, iError, iED, nLat, nLon
    logical :: IsDone
    
!    real, dimension(nLon,nLats,nAlts) :: temp, AuroralBulkIonRate, &
!         IonPrecipitationBulkIonRate, IonPrecipitationHeatingRate
    
    logical :: IsFirstTime = .true.
    
    nLat=nLine
    nLon=nLine
    
    do i=1,nLat
       do j=1,nLon
          
          eflx_ergs = ElectronEnergyFlux_C(j,i) !/ (1.0e-7 * 100.0 * 100.0)
          av_kev    = ElectronAverageEnergy_C(j,i)
          
          if (eflx_ergs > 0.02) then
             
             avee = av_kev * 1000.0
             eflux = eflx_ergs * 6.242e11
             
             a= sqrt(27.0/(2.0*cPi)) * eflux /(avee**2.5)
             
             do n=1,ED_N_Energies-1
                ED_flux(n) = &
                     a*sqrt(ed_energies(n))*exp(-1.5*ed_energies(n)/avee)/cPi
             enddo
             
             call R_ELEC_EDEP (ED_Flux, 30, ED_Energies, 3, ED_Ion, 13)!Ionizations/cc/s
             call R_ELEC_EDEP (ED_Flux, 30, ED_Energies, 3, ED_Heating, 18)!erg/cc/s
          endif
       enddo
    enddo
    
  end subroutine set_auroral_rates
  
  !============================================================================
  subroutine get_aurora(nAlt,Alt_C,NDensity_CI,Pressure_C)
    use ModCommonPlanet, ONLY:O_,O2_,N2_
    use ModPWOM, ONLY: iLine,ElectronAverageEnergy_C,ElectronEnergyFlux_C
    integer, intent(in)  :: nAlt
    real,    intent(in)  :: NDensity_CI(nAlt,nSpecies),Alt_C(nAlt),Pressure_C(nAlt)
    real,    allocatable :: temp(:), AuroralBulkIonRate(:), &
                IonPrecipitationBulkIonRate(:), IonPrecipitationHeatingRate(:)
    integer :: iLat, iLon, k, n, iED, iTop
    real    :: p, eflx_ergs, av_kev, ScaleHeight = 50.0e5 !50km
    logical :: IsDone, IsTopFirst, IsTop
    !--------------------------------------------------------------------------
    
    if (.not.allocated(temp))then
       allocate(temp(nAlt), AuroralBulkIonRate(nAlt), &
         IonPrecipitationBulkIonRate(nAlt),IonPrecipitationHeatingRate(nAlt))
    endif

    if (.not.allocated(ED_Interpolation_Weight))then
       allocate(ED_Interpolation_Weight(nAlt),ED_Interpolation_Index(nAlt),&
            AuroralIonRateO_C(nAlt),HeatingRate_C(nAlt))
    endif
    
    iLat=iLine
    iLon=iLine

    AuroralBulkIonRate (:) = 0.0
    HeatingRate_C (:) = 0.0 
    
    eflx_ergs = ElectronEnergyFlux_C(iLon,iLat) !/ (1.0e-7 * 100.0 * 100.0)
    av_kev    = ElectronAverageEnergy_C(iLon,iLat)
    
    if (eflx_ergs > 0.02) then
       iED = 1
       IsTopFirst=.true.
       do k = 1, nAlt
          p = alog(Pressure_C(k))

          
          IsDone = .false.
          IsTop  = .false.
          do while (.not.IsDone)
             !write(*,*) k, p, ED_grid(iED)
             if (ED_grid(iED) >= p .and. ED_grid(iED+1) <= p) then
                IsDone = .true.
                ED_Interpolation_Index(k) = iED
                ED_Interpolation_Weight(k) = (ED_grid(iED) - p) /  &
                     (ED_grid(iED) - ED_grid(iED+1))
             else
                if (iED == ED_N_Alts-1) then
                   IsDone = .true.
                   if (IsTopFirst) then !save altitude index at top
                      iTop = k
                      IsTopFirst=.false.
                   endif
                   ED_Interpolation_Weight(k) = 1.0
                   ED_Interpolation_Index(k) = ED_N_Alts - 1
                   IsTop = .true.
                else
                   iED = iED + 1
                endif
             endif
          enddo
          
          n = ED_Interpolation_Index(k)
          AuroralBulkIonRate(k) = ED_Ion(n) - &
               (ED_Ion(n) - ED_Ion(n+1))*ED_Interpolation_Weight(k)
          HeatingRate_C(k) = ED_Heating(n) - &
               (ED_Heating(n)-ED_Heating(n+1))*ED_Interpolation_Weight(k)
          ! if you are past top of Energy deposition model, drop contributions with scaleheight approx
          if (IsTop) then
             AuroralBulkIonRate(k) = AuroralBulkIonRate(k) * exp((Alt_C(iTop) - Alt_C(k))/ScaleHeight)
             HeatingRate_C(k)      = HeatingRate_C(k)      * exp((Alt_C(iTop) - Alt_C(k))/ScaleHeight)
          endif
       enddo
    endif
    
!    do k=1,ED_N_Alts
!       write(*,*) ED_ion(k)
!    enddo
    
    temp(1:nAlt) = 0.92 * NDensity_CI(1:nAlt,N2_) + &
         1.00 * NDensity_CI(1:nAlt,O2_) + &
         0.56 * NDensity_CI(1:nAlt,O_)
    
    AuroralIonRateO_C(1:nAlt)  = &
         0.56*AuroralBulkIonRate(1:nAlt)*&
         NDensity_CI(1:nAlt,O_)/temp(1:nAlt)

!    AuroralIonRateS(:,:,:,iO2_) = &
!         1.00*AuroralBulkIonRate*&
!         NDensity_CI(1:nLon,1:nLat,1:nAlt,iO2_)/temp
!    AuroralIonRateS(:,:,:,iN2_) = &
!         0.92*AuroralBulkIonRate*&
!         NDensity_CI(1:nLon,1:nLat,1:nAlt,iN2_)/temp

  end subroutine get_aurora
  !=============================================================================
  subroutine init_aurora
    
    implicit none
    
    integer :: ierr, iAlt, i
    real :: a
    logical :: IsDone
    
    !--------------------------------------
    ! Start doing energy deposition stuff
    !--------------------------------------
    
    call ED_Init(ierr)
    if (ierr /= 0) then
       call con_stop("Error in initilizing the energy deposition tables")
    endif
    
    call ED_Get_Grid_Size(ED_N_Alts)
    allocate(ED_grid(ED_N_Alts), stat=ierr)
    if (ierr /= 0) then
       call con_stop("Error allocating array ED_grid")
    endif
    
    allocate(ED_Ion(ED_N_Alts), stat=ierr)
    if (ierr /= 0) then
       call con_stop("Error allocating array ED_Ion")
    endif
    
    allocate(ED_Heating(ED_N_Alts), stat=ierr)
    if (ierr /= 0) then
       call con_stop("Error allocating array ED_Heating")
    endif
    
    call ED_Get_Number_of_Energies(ED_N_Energies)
    allocate(ED_Energies(ED_N_Energies), stat=ierr)
    if (ierr /= 0) then
       call con_stop("Error allocating array ED_Energies")
    endif
    
    allocate(ED_Flux(ED_N_Energies), stat=ierr)
    if (ierr /= 0) then
       call con_stop("Error allocating array ED_Flux")
    endif
    
    call ED_Get_Grid(ED_grid, .true., ierr)
    
    do iAlt = 1, ED_N_Alts
       ED_grid(iAlt) = alog(ED_grid(iAlt))
    enddo
    
    call ED_Get_Energies(ED_Energies)
    
  end subroutine init_aurora
  
end Module ModAurora
