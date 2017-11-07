Module ModSeGrid
  implicit none
  
  private !except
  ! should the code write in a verbose manner?
  logical, public :: IsVerbose=.false.

  real,public :: Time
 
  !E Grid
  integer, public :: nEnergy

  real, allocatable, public :: DeltaE_I(:),EnergyGrid_I(:)
  real,   public :: EnergyMin, EnergyMax, DeltaE
!  real, allocatable,public  :: KineticEnergy_IIC(:,:,:)

  !Alt Grid
  integer, public :: nAlt
  real, allocatable, public :: Alt_C(:) !altitude grid in cm

  !Alt Grid Extended above two stream calculation by mapping.
  ! should be nAlt+some points
  integer, public :: nAltExtended
  real, allocatable, public :: AltExtended_C(:) !altitude grid in cm
  real,public :: AltPwUpper=8000.0e5 !alt of PW upper boundary in CM

  real, allocatable, public :: DeltaPot_C(:),Efield_C(:)

  ! information on how global line number (as opposed to line number on 
  ! given proc)
  integer, public :: iLineGlobal
  
  ! public methods
  public :: allocate_grid_arrays
  public :: set_egrid
  public :: set_altgrid
  public :: BINNUM
  public :: calc_potential
  real, public :: rPlanetCM


contains
  !=============================================================================
  ! Subroutine Based on EGRID from GLOW sets up electron energy grid
  !
  ! Stan Solomon, 1/92
  !
  ! This software is part of the GLOW model.  Use is governed by the Open Source
  ! Academic Research License Agreement contained in the file glowlicense.txt.
  ! For more information see the file glow.txt.
  !
  SUBROUTINE set_egrid
    use ModPlanetConst, ONLY: Planet_, rPlanet_I, NamePlanet_I
    implicit none

    ! Local:
    Integer iEnergy, i0,i1,i2
    real A,B,B1,C,C2

    if (NamePlanet_I(Planet_).EQ.'JUPITER') then
       A = 2.     ! use 1.-2.
       B = 16.    ! use 16.-20.
       C = 0.01   ! use 0.01-0.002
       
       i0 = A/(60./740.)
       B1 = B/(60./740.)*2.
       i1=i0+B1*(B-A)/B
       C2 = 740
       i2=i1+C2*(60.-B)/60.
       ! nEnergy = i2 + 0.05/C*50  for jupiter EUV
       ! nEnergy = i2 + 0.05/C*110 for jupiter precip

       DO iEnergy=1,nEnergy
          IF (iEnergy .LE. i0) THEN
             EnergyGrid_I(iEnergy) = A/real(i0) * REAL(iEnergy)
          ELSE IF (iEnergy .LE. i1) THEN
             EnergyGrid_I(iEnergy) = B/real(B1) * REAL(iEnergy-i0+A/(B/B1))
          ELSE IF (iEnergy .LE. i2) THEN
             EnergyGrid_I(iEnergy) = 60./real(C2) * REAL(iEnergy-i1+B/(60./C2))
          ELSE
             EnergyGrid_I(iEnergy) = EXP(0.05*81.887)*EXP(C*REAL(iEnergy-i2))
             if (IsVerbose) write(*,*) iEnergy,EnergyGrid_I(iEnergy)
          ENDIF
       end do
    else   
       DO iEnergy=1,nEnergy
          IF (iEnergy .LE. 21) THEN
             EnergyGrid_I(iEnergy) = 0.5 * REAL(iEnergy)
          ELSE
             EnergyGrid_I(iEnergy) = EXP (0.05 * REAL(iEnergy+26))
          ENDIF
       end do
    end if
    
    if (IsVerbose) write(*,*) 'Energy Grid: ',EnergyGrid_I(:)

    DeltaE_I(1) = EnergyGrid_I(1)
    DO iEnergy=2,nEnergy
       DeltaE_I(iEnergy) = EnergyGrid_I(iEnergy)-EnergyGrid_I(iEnergy-1)
    End Do
    
    DO iEnergy=1,nEnergy
       EnergyGrid_I(iEnergy) = EnergyGrid_I(iEnergy) - DeltaE_I(iEnergy)/2.0
    End Do
    
    EnergyMin=EnergyGrid_I(1)
    EnergyMax=EnergyGrid_I(nEnergy)
    
    if (IsVerbose) write (*,*) sum(DeltaE_I(:)),&
         EnergyGrid_I(nEnergy)+DeltaE_I(nEnergy)/2.0,'energies!'
    
  End Subroutine set_egrid

  !=============================================================================
  ! Set Altitude Grid

  SUBROUTINE set_altgrid
    use ModPlanetConst, ONLY: Planet_, rPlanet_I, NamePlanet_I
    implicit none

    ! Local:
    real,parameter :: cKmToCm=1.0e5, cMtoCM = 1.0e2
    Integer ::iAlt
    real :: dAltExtended,dr
    real :: AltKM_C(120)
    DATA AltKM_C/     80., 81., 82., 83., 84., 85., 86., 87., 88., 89., &
         90., 91., 92., 93., 94., 95., 96., 97., 98., 99., &
         100.,101.,102.,103.,104.,105.,106.,107.,108.,109., &
         110.,111.5,113.,114.5,116.,118.,120.,122.,124.,126., &
         128.,130.,132.,134.,136.,138.,140.,142.,144.,146., &
         148.,150.,153.,156.,159.,162.,165.,168.,172.,176., &
         180.,185.,190.,195.,200.,205.,211.,217.,223.,230., &
         237.,244.,252.,260.,268.,276.,284.,292.,300.,309., &
         318.,327.,336.,345.,355.,365.,375.,385.,395.,406., &
         417.,428.,440.,453.,467.,482.,498.,515.,533.,551., &
         570.,590.,610.,630.,650.,670.,690.,710.,730.,750., &
         770.,790.,810.,830.,850.,870.,890.,910.,930.,950./
    
    if (NamePlanet_I(Planet_).EQ.'JUPITER') then

       Alt_C(1:120) = (AltKM_C/4.0+230.) * cKmToCm
       dr = 5.0

       do iAlt = 1,nAlt-120
          Alt_C(120+iAlt) = Alt_C(119+iAlt) + dr*cKmToCm
          dr = dr + iAlt/1280.0
!          write(*,*) iAlt, Alt_C(120+iAlt), dr
       end do
    else
       Alt_C(1:120) = AltKM_C * cKmToCm
    endif

    rPlanetCM=rPlanet_I(Planet_)*cMtoCM

    !add in points in extended grid above 2 stream calc. Fill in equally spaced 
    ! up to PWOM boundary
    
    !first set AltExtended_C to Alt in two stream region
    AltExtended_C(1:nAlt)=Alt_C(1:nAlt)
    
    !find grid spacing in extended region
    dAltExtended=(AltPwUpper-Alt_C(nAlt))/ real(nAltExtended-nAlt)
    
    !fill in Altitude grid in extended region
    do iAlt=nAlt+1,nAltExtended
       AltExtended_C(iAlt)=AltExtended_C(iAlt-1)+dAltExtended
    enddo

  End Subroutine set_altgrid
  


  !=============================================================================
  ! ------------------------------------------------------------------ ** 
  ! This function finds the energy grid number of the input energy E.
  !  
  INTEGER FUNCTION BINNUM(E)
    
!    use ModSeGrid,      ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
!         DeltaE_I,EnergyGrid_I, EnergyMin
    
    real, intent(in) :: E
    integer :: N, i
    logical :: IsBinFound
    !---------------------------------------------------------------------------
    IsBinFound=.false.
    i=1
    DO WHILE ((.not.IsBinFound).AND.(i.LE.nEnergy))
       IF (E.LE.EnergyGrid_I(i)+DeltaE_I(i)/2) IsBinFound=.true.
       i=i+1
    END DO
    IF (IsBinFound) THEN
       BINNUM=i-1
    ELSE
       BINNUM=nEnergy+1
    END IF
    IF (E.LE.EnergyMin) BINNUM=0
    RETURN
  END FUNCTION BINNUM


  !============================================================================
  subroutine calc_potential
    use ModConst,only: cElectronCharge ! in Coulombs 
    integer :: iAlt
    real    :: Pot !the electric potential difference from the equator [Volts]
    real,parameter :: cCmToM = 1.0e-2, cJoulesToeV=6.24150934e18
    !--------------------------------------------------------------------------
    

    !fill the DeltaPot_C array. The reference point (zero potential) is at the 
    ! top of the ionosphere. Ignore E|| effects in iono. note that you 
    ! want the total energy array to be positive so the reference point for 
    ! the potential energy should be at minimum  potential energy location. 
    ! since E|| =-dPhi/ds, DeltaPhi = int(- E||) from s_iono to s.
    DeltaPot_C(:)=0.0
    do iAlt=nAlt+1,nAltExtended
       call midpnt_int(Pot,-1.0*Efield_C(:),&
            AltExtended_C(:)*cCmToM,nAlt,iAlt,nAltExtended,1)

       ! from the potential change
       DeltaPot_C(iAlt) = -cElectronCharge*Pot*cJoulesToeV
       
    end do

!    do iAlt=1,nAltExtended
!       write(*,*) 'AltExtended_C(iAlt)*1e-5,DeltaPot_C(iAlt) ',&
!            AltExtended_C(iAlt)*1e-5,DeltaPot_C(iAlt) 
!    enddo
    
              
    
  end subroutine calc_potential

  
  subroutine allocate_grid_arrays
    use ModPlanetConst, ONLY: Planet_, rPlanet_I, NamePlanet_I
    
    if (NamePlanet_I(Planet_).EQ.'JUPITER') then
       nEnergy = 1463 ! use 1163 w/o precip; 1463 w/ precip
       nAlt    = 385
       nAltExtended = 480
    else
       nEnergy = 190
       nAlt    = 120
       nAltExtended = 220
    endif
    
    if(.not.allocated(DeltaE_I))      allocate(DeltaE_I(nEnergy))
    if(.not.allocated(EnergyGrid_I))  allocate(EnergyGrid_I(nEnergy))
    if(.not.allocated(Alt_C))         allocate(Alt_C(nAlt))
    if(.not.allocated(AltExtended_C)) allocate(AltExtended_C(nAltExtended))
    if(.not.allocated(DeltaPot_C))    allocate(DeltaPot_C(nAltExtended))
    if(.not.allocated(Efield_C))      allocate(Efield_C(nAltExtended))
    
!    if (.not.allocated(iLineGlobal_I))  allocate(iLineGlobal_I(nLine))
  end subroutine allocate_grid_arrays

  SUBROUTINE midpnt_int(sum,f,x,a,b,N,ii)
    IMPLICIT NONE
    INTEGER a,b,N,j,ii
    REAL f(N),x(N),sum
    
    sum=0.
    if ((b-a.LT.0).OR.(b.GT.N).OR.(a.LE.0)) RETURN
    if (ii.EQ.1) then
       do  j=a,b-1
          sum=sum+(x(j+1)-x(j))*(f(j)+f(j+1))*0.5
       enddo
    else      ! ii=2
       do j=a,b
          sum=sum+x(j)*f(j)
       enddo
    END IF
    RETURN
  END SUBROUTINE midpnt_int

 
end Module ModSeGrid
