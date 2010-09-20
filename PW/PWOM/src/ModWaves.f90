Module ModPwWaves
  implicit none
  
  private ! except
  public :: calc_wave_acceleration
  public :: wave_init
  real,    public, allocatable :: WaveAcceleration_C(:),Bwave_C(:)
  logical, public :: UseWaveAcceleration=.false.

contains
  !=============================================================================
  subroutine wave_init(nAlt)
    integer, intent(in) :: nAlt
    if(.not. allocated(WaveAcceleration_C)) allocate(WaveAcceleration_C(nAlt))
    WaveAcceleration_C(:)=0.0

    if (.not. allocated(Bwave_C)) allocate(Bwave_C(nAlt))
    
  end subroutine wave_init
  !=============================================================================
  subroutine calc_wave_acceleration
    use ModCommonVariables, ONLY: nAlt=>nDim,State_GV,iRho_I,nIon,DrBnd
    use ModConst, ONLY: cMu

    real, parameter :: cGperCm3toKgperM3 = 1.0e3,cM3toCm3 = 1.0e6, &
         cCmToM=1.0e-2, cMtoCm=1.0e2
    real, allocatable:: Rho_G(:), dlnRhodS_C(:)
    real :: bwave=10.0e-9
    integer :: iAlt,iIon,iAlt1,iAlt2,nWindow=10
    real :: dnAlt
    !---------------------------------------------------------------------------
    if (.not.allocated(Rho_G)) allocate(Rho_G(0:nAlt+1))
    if (.not.allocated(dlnRhodS_C)) allocate(dlnRhodS_C(nAlt))

    Rho_G(:)=0.0

    do iIon=1,nIon
       Rho_G(:)=Rho_G(:)+cGperCm3toKgperM3*State_GV(0:nAlt+1,iRho_I(iIon))
    end do
    
    do iAlt=1,nAlt
       ! Central diff
       !dlnRhodS_C(iAlt)=(log(Rho_G(iAlt+1))-log(Rho_G(iAlt-1))) & 
       !     / (2.0*DrBnd*cCmToM)
       
       ! upwind diff
       !dlnRhodS_C(iAlt)=min((log(Rho_G(iAlt))-log(Rho_G(iAlt-1))) & 
       !     / (DrBnd*cCmToM),(log(Rho_G(iAlt+1))-log(Rho_G(iAlt))) & 
       !     / (DrBnd*cCmToM))
       
       !derivative evaluated over large alt range
       iAlt1=min(iAlt+nWindow,nAlt)
       iAlt2=max(iAlt-nWindow,1)
       dnAlt= real(iAlt1-iAlt2)
       dlnRhodS_C(iAlt)=(log(Rho_G(iAlt1))-log(Rho_G(iAlt2))) & 
            / (dnAlt*DrBnd*cCmToM)

    end do
    
    call get_wave_b(nAlt,Rho_G(1:nAlt))

    do iAlt=1,nAlt
!       WaveAcceleration_C(iAlt) = &
!            -1.0*bwave**2*cMtoCm*dlnRhodS_C(iAlt)/Rho_G(iAlt)/cMu/8.0

       WaveAcceleration_C(iAlt) = &
            -1.0*Bwave_C(iAlt)**2*cMtoCm*dlnRhodS_C(iAlt)/Rho_G(iAlt)/cMu/8.0
       
    enddo
!    write(*,*) 'max(WaveAcceleration_C)',maxval(WaveAcceleration_C)
!    write(*,*) 'min(WaveAcceleration_C)',minval(WaveAcceleration_C)
!    do iAlt=1,nAlt
!       write(*,*) iAlt, WaveAcceleration_C(iAlt)
!       write(*,*) iAlt, Rho_G(iAlt)
!       write(*,*) iAlt, -1.0*bwave**2*cMtoCm*dlnRhodS_C(iAlt)/Rho_G(iAlt)/cMu/8.0
!    enddo
    deallocate(Rho_G,dlnRhodS_C)
  end subroutine calc_wave_acceleration

  !============================================================================
  subroutine get_wave_b(nAlt,Rho_C)
    
    
    use ModCommonVariables, ONLY: GmLat, RAD,uJoule2
    use ModPlanetConst,     ONLY: Earth_,DipoleStrengthPlanet_I,rPlanet_I
    use ModConst, ONLY: cMu, cLightSpeed
    use ModNumConst,        ONLY: cDegToRad
    integer, intent(in):: nAlt
    real, intent(in)   :: Rho_C(nAlt)
    real,allocatable   :: GmLat_C(:),E_C(:)
    real, parameter :: cCmToM=1.0e-2
    real    :: Lshell, uConvection,rPlanet, dipmom, E0, Bfield, AlfvenSpeed
    integer :: iAlt
    
    !--------------------------------------------------------------------------

    rPlanet = rPlanet_I(Earth_)                            ! planet's radius (m)
    dipmom  = abs(DipoleStrengthPlanet_I(Earth_)*rPlanet**3)  ! planet's dipole 
    
    
    if (.not. allocated(GmLat_C)) allocate(GmLat_C(nAlt))
    if (.not. allocated(E_C)) allocate(E_C(nAlt))
    
    Lshell = 1.0/(cos(GmLat*cDegToRad))**2.0
!    uConvection = sqrt(uJoule2) ! m/s
    uConvection = 3200.0 ! m/s
    E0 = uConvection*dipmom*sqrt(1+3.0*(sin(GmLat))**2.0)/(RAD(1)*cCmToM)**3.0

    do iAlt=1,nAlt
       E_C(iAlt) = E0 * (RAD(1)/RAD(iAlt))**1.5
       GmLat_C (iAlt) = acos(sqrt(RAD(iAlt)*cCmToM/(Lshell*rPlanet)))
       Bfield = &
            dipmom*sqrt(1+3.0*(sin(GmLat_C(iAlt)))**2.0)/(RAD(iAlt)*cCmToM)**3.0
       AlfvenSpeed = min(Bfield/sqrt(cMu*Rho_C(iAlt)),cLightSpeed)
       Bwave_C(iAlt) = E_C(iAlt)/AlfvenSpeed
       
       !write(*,*) 'iAlt,Bwave_C(iAlt)*1.0e9',iAlt,Bwave_C(iAlt)*1.0e9
    end do

  end subroutine get_wave_b

end Module ModPwWaves
