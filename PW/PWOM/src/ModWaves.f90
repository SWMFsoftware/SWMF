!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
    real, allocatable:: Rho_G(:), dRhoDr_C(:)
    real :: Ewave=250.0e-3 !V/m
    integer :: iAlt,iIon,iAlt1,iAlt2,nWindow=10
    real :: dnAlt
    real,parameter :: AltRef = 4000.0e3 !reference altitude in meters
    integer, parameter :: iAltRef=190
    real :: RhoRef ! rho at reference altitude [kg m^-3]
    real :: dRhoDrCentral,dRhoDrUpwind,B0ref
    !---------------------------------------------------------------------------
    if (.not.allocated(Rho_G)) allocate(Rho_G(0:nAlt+1))
    if (.not.allocated(dRhoDr_C)) allocate(dRhoDr_C(nAlt))

    Rho_G(:)=0.0

    ! fill Rho on altitude grid and convert to kg/m^3
    do iIon=1,nIon
       Rho_G(:)=Rho_G(:)+cGperCm3toKgperM3*State_GV(0:nAlt+1,iRho_I(iIon))
    end do
    
    ! Set density at reference altitude (assuming 8000km for reference)
    RhoRef=Rho_G(iAltRef)

    ! Get the background B in Tesla at the wave reference altitude
    call get_b0(AltRef,B0ref)

    ! loop over altitude and set the force
    do iAlt=1,nAlt
       ! Calculate the central and upwind derivative and then choose the min.
       ! Central diff
       dRhoDrCentral=(Rho_G(iAlt+1)-Rho_G(iAlt-1)) & 
            / (2.0*DrBnd*cCmToM)
       
       ! upwind diff
       dRhoDrUpwind=min((Rho_G(iAlt)-Rho_G(iAlt-1)) & 
            / (DrBnd*cCmToM),(Rho_G(iAlt+1)-Rho_G(iAlt)) & 
            / (DrBnd*cCmToM))
       
       dRhoDr_C(iAlt) = min(dRhoDrCentral,dRhoDrUpwind)
       
!       
!       !derivative evaluated over large alt range
!       iAlt1=min(iAlt+nWindow,nAlt)
!       iAlt2=max(iAlt-nWindow,1)
!       dnAlt= real(iAlt1-iAlt2)
!       dRhoDr_C(iAlt)=(Rho_G(iAlt1)-Rho_G(iAlt2)) & 
!            / (dnAlt*DrBnd*cCmToM)

       !get the wave acceleration in cm/s^2
       WaveAcceleration_C(iAlt) = &
            -0.125*(Ewave**2)/(B0ref**2)*sqrt(RhoRef)/(Rho_G(iAlt)**1.5) &
            * dRhoDr_C(iAlt)*cMtoCm

!       write(*,*) 'max(WaveAcceleration_C)',maxval(WaveAcceleration_C)
!       write(*,*) 'min(WaveAcceleration_C)',minval(WaveAcceleration_C)
       
    enddo

!    do iAlt=1,nAlt
!       write(*,*) iAlt, WaveAcceleration_C(iAlt)
!       write(*,*) iAlt, Rho_G(iAlt)
!       write(*,*) iAlt, -1.0*bwave**2*cMtoCm*dRhoDr_C(iAlt)/Rho_G(iAlt)/cMu/8.0
!    enddo
    deallocate(Rho_G,dRhoDr_C)
  end subroutine calc_wave_acceleration

  !============================================================================
  subroutine get_b0(AltRef, B0ref)
    
    use ModCommonVariables, ONLY: SmLat
    use ModPlanetConst,     ONLY: Earth_,DipoleStrengthPlanet_I,rPlanet_I
    use ModNumConst,        ONLY: cDegToRad
    real, intent(in):: AltRef !incomming reference alt [m]
    real, intent(out)   :: B0ref

    real, parameter :: cCmToM=1.0e-2
    real    :: Lshell, rPlanet, dipmom, Lat
    real    :: rRef ! reference radius
    !--------------------------------------------------------------------------

    rPlanet = rPlanet_I(Earth_)                            ! planet's radius (m)
    dipmom  = abs(DipoleStrengthPlanet_I(Earth_)*rPlanet**3)  ! planet's dipole 
    
    !set the reference radius
    rRef = rPlanet+AltRef
    ! find corresponding l-shell
    Lshell = 1.0/(cos(SmLat*cDegToRad))**2.0
    
    ! find corresponding latitude for location on l-shell
    Lat = acos(sqrt(rRef*cCmToM/(Lshell*rPlanet)))
    
    ! get the magnetic field of the reference altitude
    B0ref = &
         dipmom*sqrt(1+3.0*(sin(Lat))**2.0)/(rRef)**3.0

  end subroutine get_b0

end Module ModPwWaves
