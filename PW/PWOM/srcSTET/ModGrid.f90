Module ModSeGrid
  implicit none
  
  private !except
  ! should the code write in a verbose manner?
  logical, public :: IsVerbose=.false.

  !field line grid
  real, public, allocatable    :: FieldLineGrid_IC(:,:) 
  
  integer, public,allocatable  :: MaxAlt_IC(:,:) !max alt for given energy

  ! Grid spacings in each zone in the ionosphere
  real, public                 :: DrIono1,DrIono2,DrIono3,DrIono4

  ! Index Size of Each ionospheric zone
  integer, public :: nIono1,nIono2,nIono3,nIono4

  ! Index Size of total ionosphere
!  integer, public, allocatable :: nIono
  
  ! Altitude of Ionospheric and Plasmaspheric base 
  real, public :: BaseAltIono, BaseAltPlas

  
  ! Angle Grid
!  real, public, allocatable :: AngleGrid_ICI(:,:,:) 
  real, public, allocatable :: EqAngleGrid_IG(:,:) 
 
  integer, public,allocatable ::nTheta_II(:,:)! PA steps added for each zone on 
                                              !  each line
  real, public, allocatable   :: ThetaZone_II(:,:)! PA width for each zone 
                                                  !  on each line
  real, public,allocatable    :: dTheta_II(:,:)   ! PA spacing in each zone 
                                                  !  on each line
  integer, public,allocatable :: nThetaAlt_II(:,:)! number of PA points at 
                                                  !  each alt
  integer, public,allocatable :: nThetaAlt_IIC(:,:,:)! same as nThetaAlt_II but 
                                                  ! energy dependent for case 
                                                  ! with potential

  real, public,allocatable    :: dThetaEnd_II(:,:)   ! PA spacing at edge 
                                                     !of PA range 
  real, public,allocatable    :: dThetaEnd_III(:,:,:) !same as dThetaEnd but 
                                                      !Energy dependent for 
                                                      !case with potential
  
  real, public,allocatable    :: mu_III(:,:,:)   !Cosine of local PA at each alt
  real, public,allocatable    :: mu_IIIC(:,:,:,:)!Cosine of local PA at each alt
                                                 ! in case of potential 
                                                 ! (Energy dependent)
  !altitude of reference altitude for 1st invarient for each line and energy
  integer, public,allocatable :: nMu0RefAlt_II(:,:) 
  
  !B Grid
  real, public, allocatable :: Bfield_IC(:,:) ! Bfield in G at each alt step and
                                              ! for each line
  real, public, allocatable :: BFieldIono_I(:)! B at ionosphere for each line
  real, public,allocatable  :: BFieldEq_I(:)  ! B at equator for each line
  real, public,allocatable  :: BField0_II(:,:) ! B at mu0 ref for each line and 
                                              ! total energy

  real, public, allocatable :: Lshell_I(:) !Lshell foreach line

  integer, public :: nAngle=135 ! number of points in equatorial angle
!  integer, public :: nAngle=235 ! number of points in equatorial angle
  integer, public :: nPoint=200 ! total number of points on grid
  integer, public :: nIono =34  ! number of points in each ionosphere
  integer, public :: nPlas =132  ! number of points in plasmasphere
  integer, public :: nPlasHalf   ! number of points from BasePlas to equator
  integer, public :: nTop   ! number of points to top (in open) 
                         !   or equatorial (in closed)

  integer, public :: nLine=1! number of field lines
  integer, public :: nZone=4! number of pitchangle zones
  integer, public :: MaxTheta

  !E Grid
  character(len=10),public :: TypeGridE
  integer,public :: nEnergy
  real, allocatable, public :: DeltaE_I(:),EnergyGrid_I(:)
  real,   public :: EnergyMin, EnergyMax, DeltaE
  real, allocatable,public  :: KineticEnergy_IIC(:,:,:)

  ! Potential 
  real, public, allocatable :: Efield_IC(:,:) ! in volts/m
  real, public, allocatable :: DeltaPot_IC(:,:) !delta potential energy in eV
  real, public, allocatable :: DeltaPot0_II(:,:) !delta potential energy between
                                                 ! mu0 ref and pot ref altitutes
                                                 ! for each line and total 
                                                 ! energy
!  real,parameter :: KineticEnergyMin = 0.25 ! eV, the min kinetic energy allowed
  real,parameter :: KineticEnergyMin = 1.0 ! eV, the min kinetic energy allowed
  integer,public,allocatable:: MinEnergy_IC(:,:) ! minimum index in region of 
                                                 ! existence for energy

  logical,public :: DoIncludePotential=.false.  ! should we include potential 
                                                ! in calculation.
  
  ! information on how global line number (as opposed to line number on 
  ! given proc)
  integer, public,allocatable :: iLineGlobal_I(:)

  ! If we include extra points in PW overlap region (only when coupling PWOM)
  logical,public :: UsePwRegion = .false.
  integer,public :: nPwRegion=50 ! points in overlap region above ionosphere


  ! public methods
  public :: allocate_grid_arrays
  public :: init_se_grid
  public :: se_grid_test
  public :: create_se_test_grid
  public :: BINNUM
  public :: plot_grid_pot
  public :: set_grid_pot
  public :: update_grid

  real, public :: rPlanetCM
contains
  !============================================================================
  subroutine init_se_grid(iLine)
    use ModPlanetConst, ONLY: Planet_, rPlanet_I, DipoleStrengthPlanet_I
    use ModNumConst,    ONLY: cPi
    integer, intent(in) :: iLine
    real, parameter :: cMtoCM = 1.0e2
    real,parameter :: cTeslaToGauss=10000.0 ! convert tesla to gauss
    real    :: Biono, Beq, MLAT1, PhiBasePlas, QO, SphiO
    !--------------------------------------------------------------------------

    nIono=nIono1+nIono2+nIono3+nIono4
    BaseAltPlas = &
         DrIono1*nIono1+DrIono2*nIono2+DrIono3*nIono3+DrIono4*nIono4
    BaseAltIono = 9e6 !cm
    
    nPlas=nPoint-2*nIono
    nTop = nPoint/2
    nPlasHalf=nPlas/2
    
    rPlanetCM=rPlanet_I(Planet_)*cMtoCM
    ! Set field line info
    PhiBasePlas=&
         ACOS(SQRT((rPlanetCM+BaseAltPlas)/(Lshell_I(iLine)*rPlanetCM)))
    SphiO=SIN(PhiBasePlas)
    QO=SQRT(1+3*SphiO**2)
    MLAT1=PhiBasePlas*180./cPi
!    Beq   = 0.31/Lshell_I(iLine)**3
    Beq   = abs(DipoleStrengthPlanet_I(Planet_))/Lshell_I(iLine)**3*cTeslaToGauss
!    Biono = 0.31*QO/(Lshell_I(iLine)**3*(1-SphiO**2)**3)
    Biono = abs(DipoleStrengthPlanet_I(Planet_))*QO&
         /(Lshell_I(iLine)**3*(1-SphiO**2)**3)*cTeslaToGauss
    BFieldIono_I(iLine) = Biono
    BFieldEq_I(iLine)   = Beq
    if(IsVerbose) write(*,*) 'calc bfield and s grid for iLine = ', iLine
    call calc_bfield_sgrid(iLine,Biono,PhiBasePlas)
    if(IsVerbose) write(*,*) 'calc equatorial PA grid', iLine
    call calc_equatorial_pitchangle(iLine,Beq,Biono)
    if(IsVerbose) write(*,*) 'calc mu'
    call calc_mu(iLine,Beq)
    if(IsVerbose) write(*,*) 'calc energy grid', iLine
    call calc_energy_grid
    
  end subroutine init_se_grid

  !============================================================================
  subroutine calc_bfield_sgrid(iLine, Biono,PhiBasePlas)
    use ModNumConst,    ONLY: cDegToRad,cRadToDeg
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I
    integer, intent(in) :: iLine
    real   , intent(in) :: Biono,PhiBasePlas
    integer :: iAlt
    real :: Phi, LastPhi
    real :: DeltaS, dPhi
    !length of iono regions
    real :: LengthIono1,LengthIono2,LengthIono3,LengthIono4 

    !variables for PW overlap region
    real :: PhiTopPw, dPhiPw
    real    :: TopAltPw=8000.0e5
    !--------------------------------------------------------------------------

!    if(.not.allocated(Bfield_IC)) allocate(Bfield_IC(nLine,nPoint))


    ! set dPhi
    if (UsePwRegion) then
       !choose top of Pw grid based on planet
       select case(NamePlanet_I(Planet_))
       case('EARTH')
          TopAltPw=8000.0e5
       case('JUPITER')
          TopAltPw=60000.0e5
       end select


       PhiTopPw=&
            ACOS(SQRT((rPlanetCM+TopAltPw)/(Lshell_I(iLine)*rPlanetCM)))
       dPhi=PhiTopPw/(nPlasHalf-nPwRegion)
       dPhiPw=(PhiBasePlas-PhiTopPW)/nPwRegion
    else
       dPhi=PhiBasePlas/nPlasHalf
       PhiTopPw=0.0
       nPwRegion=0
       dPhiPw=0.0
    endif
!    write(*,*) 'nPlasHalf,dPhi*cRadToDeg',nPlasHalf,dPhi*cRadToDeg

    !set Length of iono regions
    LengthIono1 = DrIono1*(nIono1-1)
    LengthIono2 = DrIono2*(nIono2)
    LengthIono3 = DrIono3*(nIono3)
    LengthIono4 = DrIono4*(nIono4)

    
    do iAlt = 1, nPoint
       if (iAlt <= nIono1) then
          !set alt zone 1 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = BaseAltIono+(iAlt-1)*DrIono1
       elseif(iAlt <= nIono1+nIono2) then
          !set alt zone 2 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+(iAlt-nIono1)*DrIono2
       elseif(iAlt <=nIono1+nIono2+nIono3) then
          !set alt zone 3 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               BaseAltIono+LengthIono1+LengthIono2 &
               +(iAlt-nIono1-nIono2)*DrIono3
       elseif(iAlt <=nIono1+nIono2+nIono3+nIono4) then
          !set alt zone 4 of N. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
                  BaseAltIono+LengthIono1+LengthIono2+LengthIono3&
                  +(iAlt-nIono1-nIono2-nIono3)*DrIono4

       elseif(iAlt <=nIono1+nIono2+nIono3+nIono4+nPwRegion.and.UsePwRegion) then
          !set alt zone in PW region above ionosphere
          Phi = abs (PhiBasePlas-(iAlt-nIono)*dPhiPW)
          LastPhi = abs (PhiBasePlas-(iAlt-nIono-1)*dPhiPW)
          ! Get B and Delta s at each point
          call get_b_deltaS_point(Phi,LastPhi,DeltaS,Lshell_I(iLine),&
               Bfield_IC(iLine,iAlt))
          FieldLineGrid_IC(iLine,iAlt)=FieldLineGrid_IC(iLine,iAlt-1)+DeltaS

       elseif(iAlt < nPoint-nIono-nPwRegion) then
          ! set plasmaspheric region
          ! Get polar angle at fieldline point
          if(UsePwRegion) then
             Phi = abs (PhiTopPw-(iAlt-nIono-nPwRegion)*dPhi)
             LastPhi = abs (PhiTopPw-(iAlt-nIono-nPwRegion-1)*dPhi)
          else
             Phi = abs (PhiBasePlas-(iAlt-nIono)*dPhi)
             LastPhi = abs (PhiBasePlas-(iAlt-nIono-1)*dPhi)
          endif
          ! Get B and Delta s at each point
          call get_b_deltaS_point(Phi,LastPhi,DeltaS,Lshell_I(iLine),Bfield_IC(iLine,iAlt))
          FieldLineGrid_IC(iLine,iAlt)=FieldLineGrid_IC(iLine,iAlt-1)+DeltaS

       elseif(iAlt < nPoint-nIono .and. UsePwRegion) then
          !set alt zone in PW region above ionosphere
          Phi = abs (PhiBasePlas-(iAlt-nIono)*dPhiPW)
          LastPhi = abs (PhiBasePlas-(iAlt-nIono-1)*dPhiPW)
          ! Get B and Delta s at each point
          call get_b_deltaS_point(Phi,LastPhi,DeltaS,Lshell_I(iLine),Bfield_IC(iLine,iAlt))
          FieldLineGrid_IC(iLine,iAlt)=FieldLineGrid_IC(iLine,iAlt-1)+DeltaS

       elseif(iAlt < nPoint-nIono1-nIono2-nIono3) then
          !set alt zone 4 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono-1)&
               +(iAlt-nPlas-nIono+1)*DrIono4
       elseif(iAlt < nPoint-nIono1-nIono2) then
          !set alt zone 3 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-nIono3-1)&
               +(iAlt-nIono4-nPlas-nIono+1)*DrIono3
       elseif(iAlt < nPoint-nIono1)then
          !set alt zone 2 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-nIono2-1)&
               +(iAlt-nIono4-nIono3-nPlas-nIono+1)*DrIono2
       elseif(iAlt <= nPoint)then
          !set alt zone 1 of S. ionosphere
          FieldLineGrid_IC(iLine,iAlt) = &
               FieldLineGrid_IC(iLine,nPoint-nIono1-1)&
               +(iAlt-nIono4-nIono3-nIono2-nPlas-nIono+1)*DrIono1
       endif
       
!       write(*,*) 'iAlt,FieldLineGrid_IC(iLine,iAlt)',iAlt,FieldLineGrid_IC(iLine,iAlt)/1e5
    end do

    ! Set magnetic field in ionosphere to value at top of ionosphere
    Bfield_IC(iLine,1:nIono)= Biono
    Bfield_IC(iLine,nPoint-nIono:nPoint)= Biono
!    write(*,*) 'Biono',Biono
  end subroutine calc_bfield_sgrid
  !============================================================================
  
  subroutine calc_equatorial_pitchangle(iLine,Beq,Biono)
    use ModNumConst,    ONLY: cPi
    integer, intent(in) :: iLine
    real,    intent(in) :: Beq, Biono
    integer :: iAlt, iAngle
    real    :: LocalThetaMax

    integer, allocatable :: tmp_array(:)
    !--------------------------------------------------------------------------
    
    
    ! Set Maximum number of points in equatorial PA
    MaxTheta = sum(nTheta_II(iLine,:))
    
    if (.not.allocated(tmp_array)) allocate(tmp_array(0:MaxTheta))
    
    ! Set the theta range for each zone and each line
    !  (currently assumes 4 zones)
    ThetaZone_II(iLine,1) = asin(sqrt(Beq/Biono)) !zone 1 is the loss cone
    ThetaZone_II(iLine,2) = ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,3) = 1.2-2.0*ThetaZone_II(iLine,1)
    ThetaZone_II(iLine,4) = 0.5*cPi-1.2
    
    ! Get DeltaTheta for each zone
    dTheta_II(iLine,:) = ThetaZone_II(iLine,:)/nTheta_II(iLine,:)
    
    ! Fill in compuational region of equatorial PA grid
    do iAngle=1,MaxTheta
       if (iAngle <= nTheta_II(iLine,1)) then 
          EqAngleGrid_IG(iLine,iAngle) = dTheta_II(iLine,1) * iAngle
       elseif (iAngle <= sum(nTheta_II(iLine,1:2))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,2) * (iAngle-nTheta_II(iLine,1)) &
               + ThetaZone_II(iLine,1)
       elseif (iAngle <= sum(nTheta_II(iLine,1:3))) then
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,3) * (iAngle-sum(nTheta_II(iLine,1:2))) &
               + sum(ThetaZone_II(iLine,1:2))
       else 
          EqAngleGrid_IG(iLine,iAngle) = &
               dTheta_II(iLine,4) * (iAngle-sum(nTheta_II(iLine,1:3))) &
               + sum(ThetaZone_II(iLine,1:3))
       end if
 
   end do
    
    ! Fill in Ghost Cells of Equatorial Angle Grid
    EqAngleGrid_IG(iLine,0) = 0.0
    EqAngleGrid_IG(iLine,MaxTheta) = 0.5 * cPi
    
    ! Fill in nThetaAlt array (Number of points in theta for each alt)
    do iAlt = 1, nPoint
       ! find equatorial PA corresponding to local Pi/2  
       LocalThetaMax = asin(sqrt(min(1.0,Beq/Bfield_IC(iLine,iAlt))))
       
       ! From LocalThetaMax and EqAngleGrid find nThetaAlt
       Tmp_array = 0
       where(EqAngleGrid_IG(iLine,:) <= LocalThetaMax)
          Tmp_array = 1 
       end where
       nThetaAlt_II(iLine,iAlt) = max(sum(Tmp_array)-1,0)
       
       ! Set dThetaEnd, the delta theta at the end of the PA grid
       dThetaEnd_II(iLine,iAlt)=EqAngleGrid_IG(iLine,nThetaAlt_II(iLine,iAlt)) &
            - EqAngleGrid_IG(iLine,nThetaAlt_II(iLine,iAlt)-1)
    enddo
    

  end subroutine calc_equatorial_pitchangle

  !============================================================================
  ! subroutine to calculate the array, mu_III, that holds cos(local PA)
  subroutine calc_mu(iLine,Beq)
    integer, intent(in) :: iLine
    real   , intent(in) :: Beq
    
    real    :: coef
    integer :: iAngle, iAlt
    !---------------------------------------------------------------------------
    
    ALONG_LINE: do iAlt= 1, nPoint
       coef = min(1.0,Beq/Bfield_IC(iLine,iAlt))
       
       ! set mu at 0 and 90 local PA
       mu_III(iLine,0,iAlt) = 1.0
       mu_III(iLine,nThetaAlt_II(iLine,iAlt),iAlt) = 0.0
       ! set mu at angles between 0 and 90 
      ANGLE: do iAngle=1,nThetaAlt_II(iLine,iAlt)-1
          mu_III(iLine,iAngle,iAlt) = &
               sqrt(1.0-1.0/coef*(1.0-(cos(EqAngleGrid_IG(iLine,iAngle)))**2.0))
       end do ANGLE
    end do ALONG_LINE
  end subroutine calc_mu
  !============================================================================
  ! subroutine to calculate the array, mu_III, that holds cos(local PA) in the 
  ! case of the field aligned potential. Note that reference altitudes must be 
  ! set first
  subroutine calc_mu_pot(iLine)
    integer, intent(in) :: iLine

    real    :: eta
    integer :: iAngle, iAlt, iEnergy
    !---------------------------------------------------------------------------
    !now pitchangle is energy dependent
    !write(*,*)'MaxAlt_IC(iLine,1)',MaxAlt_IC(iLine,1)
    !call con_stop('')
    ENERGY: do iEnergy=1,nEnergy
       ALONG_LINE: do iAlt= 1, nPoint
          ! check if we are outside of altitude range if we are then cycle the 
          ! loop along the line
          if (iAlt>MaxAlt_IC(iLine,iEnergy) &
               .and. iAlt<nPoint-MaxAlt_IC(iLine,iEnergy)) cycle ALONG_LINE

          eta = min(1.0,(Bfield0_II(iLine,iEnergy)/Bfield_IC(iLine,iAlt)) &
               *(EnergyGrid_I(iEnergy)-DeltaPot_IC(iLine,iAlt))   &
               /(EnergyGrid_I(iEnergy)-DeltaPot0_II(iLine,iEnergy)))

          ! set mu at 0 and 90 local PA
          mu_IIIC(iLine,0,iEnergy,iAlt) = 1.0
          mu_IIIC(iLine,nThetaAlt_IIC(iLine,iEnergy,iAlt),iEnergy,iAlt) = 0.0

          ! set mu at angles between 0 and 90 
          ANGLE: do iAngle=1,nThetaAlt_IIC(iLine,iEnergy,iAlt)-1
             ! calculation of mu from mu0 in case of potential. formula from 
             ! Liemohn et al., 1997
             mu_IIIC(iLine,iAngle,iEnergy,iAlt) = &
                  sqrt(1.0-1.0/eta*(1.0-&
                  (cos(EqAngleGrid_IG(iLine,iAngle)))**2.0))
          end do ANGLE
       end do ALONG_LINE
    end do ENERGY
  end subroutine calc_mu_pot
  !============================================================================
  !============================================================================
  ! Subroutine SpaceVar calculates h and Bfield for a given spatial step.
  !  VARIABLE DESCRIPTIONS
  !      Q       Used in distance calculations, SQRT(1+3*(SIN(phi))**2)
  !      lQ      Same as Q, but for lphi
  !      phi     Azimuthal angle for the base of the plasmasphere; rad
  !      L       Equatorial distance in Earth radii of flux tube
  !      Re      Radius of the Earth; 8.378E8 cm
  !      lphi    Azimuthal angle for the previous step; radians
  !      h`       Size of this step; cm
  !      B       Magnetic field strength for this step; G
  !
  SUBROUTINE get_b_deltaS_point(phi,lphi,h,L,B)
    use ModPlanetConst, ONLY: Planet_,DipoleStrengthPlanet_I
    REAL phi,lphi,h,L,B,Q,lQ
    real,parameter :: cTeslaToGauss=10000.0 ! convert tesla to gauss
    !---------------------------------------------------------------------------
    
    Q=SQRT(1+3*(SIN(phi))**2)
    lQ=SQRT(1+3*(SIN(lphi))**2)
    h=ABS(.5*L*rPlanetCM*(lQ*SIN(lphi)-Q*SIN(phi)+1/SQRT(3.)*LOG((SQRT(3.) &
         *SIN(lphi)+lQ)/(SQRT(3.)*SIN(phi)+Q))))
!    B=0.31*Q/(L**3*(COS(phi))**6)
    B=abs(DipoleStrengthPlanet_I(Planet_))*Q/(L**3*(COS(phi))**6)*cTeslaToGauss
    RETURN
  END SUBROUTINE get_b_deltaS_point
  !============================================================================
  subroutine calc_energy_grid
    integer :: iEnergy
    real    :: EnergySide1,EnergySide2

    select case(TypeGridE)
    case('ConstDE')
       EnergyMin=EnergyMax-nEnergy*DeltaE
       do iEnergy = 1, nEnergy
          !Constant DeltaE grid
          EnergyGrid_I(iEnergy)=EnergyMin+DeltaE*(iEnergy-0.5)
          DeltaE_I(iEnergy)=DeltaE
       end do
       DeltaE_I(nEnergy+1)=DeltaE
    case('ConstLogDE')
       ! Constant log(DeltaE) grid
       EnergySide2=alog(EnergyMax)+0.5*DeltaE-DeltaE*nEnergy
       EnergyMin=exp(EnergySide2)
       do iEnergy =1, nEnergy
          EnergySide1=EnergySide2
          EnergySide2=EnergySide1+DeltaE
          DeltaE_I(iEnergy) = exp(EnergySide2)-exp(EnergySide1)
          EnergyGrid_I(iEnergy) = exp(EnergySide1)+0.5*DeltaE_I(iEnergy)
       end do
       Energymax=exp(EnergySide2)
    case('Staggered')
       !Staggered dE grid. Two regions with 2eV and 4eV spacing
       EnergyMin=1.0			! Input Emax is Emax
       EnergyGrid_I(1)=1.5
       DeltaE_I(1)=1.0
       do iEnergy = 2, 20	
          EnergyGrid_I(iEnergy)=REAL(2*iEnergy)-1.0	
	  DeltaE_I(iEnergy)=2.0
       end do
       
       do iEnergy=21,21+nEnergy
          EnergyGrid_I(iEnergy)=40.+(0.5+iEnergy-21)*DeltaE
	  DeltaE_I(iEnergy)=DeltaE
       enddo
       EnergyMax=EnergyGrid_I(nEnergy)+0.5*DeltaE_I(nEnergy)
       
    case('Geometric')
       EnergyMin=2.			! Given unless input file is changed
       DeltaE_I(1)=DeltaE
       EnergyGrid_I(1)=EnergyMin+0.5*(1)
       do iEnergy=2,nEnergy
          ! Input Emax is dE growth factor
	  DeltaE_I(iEnergy)=EnergyMax*DeltaE_I(iEnergy-1)
          EnergyGrid_I(iEnergy)=&
               EnergyGrid_I(iEnergy-1)+.5*(DeltaE_I(iEnergy-1)+DeltaE_I(iEnergy))
       end do
       EnergyMax=EnergyGrid_I(nEnergy)+.5*DeltaE_I(nEnergy)
    end select
 
  end subroutine calc_energy_grid

  !=============================================================================
  ! ------------------------------------------------------------------ ** 
  ! This function finds the energy grid number of the input energy E.
  !  VARIABLE DESCRIPTIONS
  !      Elen    Array size for energy step variables, >= Jo
  !      ener    Energy array; eV
  !      Emin    Minumum energy; eV
  !      del     Energy step array; eV 
  !      Jo      Number of energy steps
  !      N       Flag indicating whether the energy is in the range
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

  !=============================================================================
  subroutine calc_potential(iLine)
    use ModMath, only: midpnt_int
    use ModConst,only: cElectronCharge ! in Coulombs 
    integer, intent(in) :: iLine
    integer :: iAlt
    real    :: Pot !the electric potential difference from the equator [Volts]
    real,parameter :: cCmToM = 1.0e-2, cJoulesToeV=6.24150934e18
    !--------------------------------------------------------------------------
    

    !fill the DeltaPot_IC array. The reference point (zero potential) is at the 
    ! top of the ionosphere. Ignore E|| effects in iono. note that you 
    ! want the total energy array to be positive so the reference point for 
    ! the potential energy should be at minimum  potential energy location. 
    ! since E|| =-dPhi/ds, DeltaPhi = int(- E||) from s_iono to s.
    DeltaPot_IC(iLine,:)=0.0
    do iAlt=nIono+1,nTop
       call midpnt_int(Pot,-1.0*Efield_IC(iLine,:),&
            FieldLineGrid_IC(iLine,:)*cCmToM,nIono,iAlt,nPoint,1)
       ! from the potential change
       DeltaPot_IC(iLine,iAlt) = -cElectronCharge*Pot*cJoulesToeV

       ! assume that other hemisphere has identical potential structure
       DeltaPot_IC(iLine,nPoint-iAlt) = DeltaPot_IC(iLine,iAlt) 
    end do
    
  end subroutine calc_potential

  !=============================================================================
  subroutine set_energy_bounds(iLine)
    integer, intent(in) :: iLine
    integer, allocatable :: tmp_array(:)
    integer :: iPoint
    ! --------------------------------------------------------------------------
    
    !allocate temporary array to help determine bounds
    if (.not.allocated(tmp_array)) allocate(tmp_array(nEnergy))

    !For each altitude we need to restrict the range of the total energy 
    ! array we consider when a potential is included, whenever 
    ! |Potential Energy|>TotalEnergy-minKinetic Energy 
    ! there is not enough KE to overcome the potential so it will not be in 
    ! our region of exisitence
    
    do iPoint=1,nPoint
       Tmp_array = 0
       where(EnergyGrid_I(:) < abs(DeltaPot_IC(iLine,iPoint))+KineticEnergyMin)
          Tmp_array = 1 
       end where
       MinEnergy_IC(iLine,iPoint) = max(sum(Tmp_array)+1,1)
    end do
    
    !deallocate array to save memory
    deallocate (tmp_array)
  end subroutine set_energy_bounds
  !=============================================================================
  subroutine set_alt_bounds(iLine)
    integer, intent(in) :: iLine
    integer, allocatable :: tmp_array(:)
    integer :: iEnergy
    ! --------------------------------------------------------------------------
    
    !allocate temporary array to help determine bounds
    if (.not.allocated(tmp_array)) allocate(tmp_array(nTop))

    !For each total energy we need to restrict the range of the altitude  
    ! array we consider when a potential is included, whenever 
    ! |Potential Energy|>TotalEnergy-minKinetic Energy 
    ! there is not enough KE present so it will not be in 
    ! our region of exisitence
    
    do iEnergy=1,nEnergy
       Tmp_array = 0
       where(EnergyGrid_I(iEnergy) > abs(DeltaPot_IC(iLine,1:nTop))&
            +KineticEnergyMin)
          Tmp_array = 1 
       end where
       MaxAlt_IC(iLine,iEnergy) = max(sum(Tmp_array),1)
       
!       write(*,*) 'iEnergy,MaxAlt_IC(iLine,iEnergy)',iEnergy,MaxAlt_IC(iLine,iEnergy)
    end do
!    call con_stop('')
    
    !deallocate array to save memory
    deallocate (tmp_array)
  end subroutine set_alt_bounds
 
 !=============================================================================
  ! This subroutine finds the reference altitude for the first invariant 
  ! for each energy by minimizing [(Total Energy - Potential Energy)/B(s)]^-1
  subroutine locate_reference_alt_for_mu0(iLine)
    integer, intent(in) :: iLine
    integer :: iEnergy, iPoint, nMu0Alt
    logical :: IsFoundRefAlt
    real    :: eta
    !---------------------------------------------------------------------------
    
    
    ! for each energy find the reference altitude by fixing the reflection 
    ! altitude as the reference altitude and seeing the pi/2 local pitch angle 
    ! at each at altitude can map into that reference. If not, then lower the 
    ! reference altitude by one and try again. 
    ENERGY_LOOP: do iEnergy = 1,nEnergy
       !Initial guess for the reference altitude is the reflection altitude
       nMu0Alt = MaxAlt_IC(iLine,iEnergy)

       !test the initial guess and move down to next altitude if it doesnt work
       IsFoundRefAlt = .false.
       FIND_REF_ALT: do while (.not.IsFoundRefAlt)
          !set the minimum B
          BField0_II(iLine,iEnergy) = Bfield_IC(iLine,nMu0Alt)
          
          !set potential difference to mu0 reference altitude
          DeltaPot0_II(iLine,iEnergy) = DeltaPot_IC(iLine,nMu0Alt)
          
          !for each altitude, the pi/2 local PA maps to the reference PA using
          ! the following relation: mu0=sqrt(1-eta) where eta is defined by
          ! eta = (1-PotEnergy(s))/(1-PotEnergy0) * B0/B(s)
          ALT_LOOP: do iPoint = 1, MaxAlt_IC(iLine,iEnergy)
             !start by calculating eta for each altitude
             eta = ((EnergyGrid_I(iEnergy)-DeltaPot_IC(iLine,iPoint))&
                  /((EnergyGrid_I(iEnergy)-DeltaPot_IC(iLine,nMu0Alt))))&
                  * (Bfield_IC(iLine,nMu0Alt)/Bfield_IC(iLine,iPoint))
             
             ! if eta exceeds 1 than pi/2 local PA cannot map to mu0 grid so 
             ! then this choice of reference altitude does not work
             if (eta > 1.0) then
                ! set the new reference altitude to one less than before and 
                !  cycle the while loop
                nMu0Alt=nMu0Alt-1
                CYCLE FIND_REF_ALT
             endif
          enddo ALT_LOOP
          
          ! If eta<1 for all altitude (which it must be to reach this point)
          !  then the mu0 reference altitude is found and the stop condition
          !  can be set. 
          
          nMu0RefAlt_II(iLine,iEnergy) = nMu0Alt
          IsFoundRefAlt =.true.
       end do FIND_REF_ALT
       !write(*,*) iEnergy,BField0_II(iLine,iEnergy),nMu0Alt,MaxAlt_IC(iLine,iEnergy)
    enddo ENERGY_LOOP
    !call con_stop('')

  end subroutine locate_reference_alt_for_mu0

  !=============================================================================
  subroutine find_angle_boundary(iLine)
    integer, intent(in) :: iLine
    integer :: iPoint,iEnergy
    real    :: eta, LocalThetaMax
    integer, allocatable :: tmp_array(:)
    !--------------------------------------------------------------------------
    
    if (.not.allocated(tmp_array)) allocate(tmp_array(0:MaxTheta))
    
    ! In the presence of an electric potential, find the number of points in 
    !  PA to consider at each altitude
    do iEnergy =1,nEnergy
       do iPoint=1,MaxAlt_IC(iLine,iEnergy)
          eta = (Bfield0_II(iLine,iEnergy)/Bfield_IC(iLine,iPoint)) &
               *(EnergyGrid_I(iEnergy)-DeltaPot_IC(iLine,iPoint))   &
               /(EnergyGrid_I(iEnergy)-DeltaPot0_II(iLine,iEnergy))

          ! Find the pitch angle at the reference altitude that corresponds to 
          !  local PA Pi/2. 
          LocalThetaMax = acos(sqrt(1.0-min(1.0,eta)))
          
          ! From LocalThetaMax and EqAngleGrid find nThetaAlt
          tmp_array = 0
          ! note a small tolerance added to the inequality test
          where(EqAngleGrid_IG(iLine,:) <= LocalThetaMax+.0000001)
             tmp_array = 1 
          end where
          nThetaAlt_IIC(iLine,iEnergy,iPoint) = max(sum(Tmp_array)-1,0)
          if(nThetaAlt_IIC(iLine,iEnergy,iPoint)==136) then
             write(*,*) iEnergy,iPoint,LocalThetaMax
             call con_stop('')
          endif
          ! Set dThetaEnd, the delta theta at the end of the PA grid
          dThetaEnd_III(iLine,iEnergy,iPoint)=&
               EqAngleGrid_IG(iLine,nThetaAlt_IIC(iLine,iEnergy,iPoint)) &
               - EqAngleGrid_IG(iLine,nThetaAlt_IIC(iLine,iEnergy,iPoint)-1)
          
          ! fill in conjugate hemisphere assuming symmetry
          nThetaAlt_IIC(iLine,iEnergy,nPoint-iPoint) = &
               nThetaAlt_IIC(iLine,iEnergy,iPoint)
          
          dThetaEnd_III(iLine,iEnergy,nPoint-iPoint) = &
               dThetaEnd_III(iLine,iEnergy,iPoint)
     
       enddo
       !fill in last point in conjugate hemisphere
       nThetaAlt_IIC(iLine,iEnergy,nPoint) = &
            nThetaAlt_IIC(iLine,iEnergy,1)
       dThetaEnd_III(iLine,iEnergy,nPoint) = &
            dThetaEnd_III(iLine,iEnergy,1)
    enddo
    
    deallocate(tmp_array)
  end subroutine find_angle_boundary
  !=============================================================================
  ! Fill the kinetic energy array
  subroutine calc_kinetic_energy(iLine)
    integer, intent(in) :: iLine
    integer :: iPoint,iEnergy
    !---------------------------------------------------------------------------
    
    !find the kinetic energy by subtracting the potential energy from the 
    ! total energy
    do iEnergy =1,nEnergy
       do iPoint=1,MaxAlt_IC(iLine,iEnergy)
          KineticEnergy_IIC(iLine,iEnergy,iPoint) = &
               EnergyGrid_I(iEnergy) - DeltaPot_IC(iLine,iPoint)
       enddo
       do iPoint=nPoint-MaxAlt_IC(iLine,iEnergy),nPoint
          KineticEnergy_IIC(iLine,iEnergy,iPoint) = &
               EnergyGrid_I(iEnergy) - DeltaPot_IC(iLine,iPoint)
       enddo
    enddo
  end subroutine calc_kinetic_energy
  !=============================================================================
  subroutine set_grid_pot(iLine)
    integer, intent(in) :: iLine
    !when including a potential set up energy dependent grid for calculation
    ! in new variables
    
    DoIncludePotential=.true.
    
    ! calculate the potential energy
    call calc_potential(iLine)
    
    ! set the bounds on the altitude
    call set_alt_bounds(iLine)
    
    ! get mu0 reference altitude
    call locate_reference_alt_for_mu0(iLine)
    
    ! set the pitchangle range for each energy and altitude 
    call find_angle_boundary(iLine)

    ! find the kinetic energy for each altitude and PA
    call calc_kinetic_energy(iLine)
    
    ! set local mu from mu0 grid
    call calc_mu_pot(iLine)
  end subroutine set_grid_pot
  !=============================================================================

  subroutine allocate_grid_arrays
    
    ! Allocate PA related grid
    if(.not.allocated(nTheta_II))       allocate(nTheta_II(nLine,nZone))
    if(.not.allocated(ThetaZone_II))    allocate(ThetaZone_II(nLine,nZone))
    if(.not.allocated(dTheta_II))       allocate(dTheta_II(nLine,nZone))
    if(.not.allocated(dThetaEnd_II))    allocate(dThetaEnd_II(nLine,nPoint))
    if(.not.allocated(nThetaAlt_II))    allocate(nThetaAlt_II(nLine,nPoint))
    if(.not.allocated(EqAngleGrid_IG))  allocate(EqAngleGrid_IG(nLine,0:nAngle))
    if(.not.allocated(Bfield_IC))       allocate(Bfield_IC(nLine,nPoint))
    if(.not.allocated(BFieldIono_I))    allocate(BFieldIono_I(nLine))
    if(.not.allocated(BFieldEq_I))      allocate(BFieldEq_I(nLine))
    if(.not.allocated(Lshell_I))        allocate(Lshell_I(nLine))
    if(.not.allocated(FieldLineGrid_IC))allocate(FieldLineGrid_IC(nLine,nPoint))
    if(.not.allocated(DeltaE_I))        allocate(DeltaE_I(nEnergy+1))
    if(.not.allocated(EnergyGrid_I))    allocate(EnergyGrid_I(nEnergy))
    if(.not.allocated(mu_III))          allocate(mu_III(nLine,0:nAngle,nPoint))

    ! variables needed for case with potential. to increase efficiency these 
    ! may be allocated separately in the future 

    if(.not.allocated(Efield_IC))       allocate(Efield_IC(nLine,nPoint))
    if(.not.allocated(DeltaPot_IC))     allocate(DeltaPot_IC(nLine,nPoint))
    if(.not.allocated(KineticEnergy_IIC))allocate(KineticEnergy_IIC(nLine,nEnergy,nPoint))
    if(.not.allocated(MinEnergy_IC))    allocate(MinEnergy_IC(nLine,nPoint))
    if (.not.allocated(nMu0RefAlt_II))  allocate(nMu0RefAlt_II(nLine,nEnergy))
    if (.not.allocated(DeltaPot0_II))   allocate(DeltaPot0_II(nLine,nEnergy))
    if (.not.allocated(Bfield0_II))     allocate(Bfield0_II(nLine,nEnergy))
    if (.not.allocated(MaxAlt_IC))      allocate(MaxAlt_IC(nLine,nEnergy))
    if (.not.allocated(nThetaAlt_IIC))  allocate(nThetaAlt_IIC(nLine,nEnergy,nPoint))
    if (.not.allocated(dThetaEnd_III))  allocate(dThetaEnd_III(nLine,nEnergy,nPoint))
    if (.not.allocated(mu_IIIC))        allocate(mu_IIIC(nLine,0:nAngle,nEnergy,nPoint))
    if (.not.allocated(iLineGlobal_I))  allocate(iLineGlobal_I(nLine))
  end subroutine allocate_grid_arrays

  !============================================================================
  ! save grid plot for verification
  subroutine plot_grid
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg
    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    integer, parameter :: nDim =2, nVar=2, S_=1, PA_=2,B_=1
    character(len=100),parameter :: NamePlotVar='S PA B PA g r'
    character(len=12) :: NamePlot='GridPlot.out'
    character(len=*),parameter :: NameHeader='Grid output'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iLine,iPoint
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,nAngle), PlotState_IIV(nPoint,nAngle,nVar))

    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
          do iPoint=1,nPoint
             do iAngle=1,nAngle
!             Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)/6375.0e5
             Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)/rPlanetCM
             Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             if (iAngle <= nThetaAlt_II(iLine,iPoint) )then
                PlotState_IIV(iPoint,iAngle,B_) = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo
       enddo
            
       !Plot grid for given line
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= 1,TimeIn=1.0,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
    end do
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_grid


  !============================================================================
  ! Plot region of existance for a particular total energy in the case  
  ! when an electric potential is included
  subroutine plot_grid_pot(iLine,nStep,iEnergy,time)
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst,ONLY: Planet_, rPlanet_I
    integer, intent(in) :: iLine, nStep, iEnergy
    real,    intent(in) :: time

    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !real, parameter     :: rEarthCM = 6375.0e5
    !grid parameters
    integer, parameter :: nDim =2, nVar=2, S_=1, PA_=2,B_=1
    real, parameter :: cMtoCM = 1.0e2

    !set the corresponding energy channels for E1-E5
    character(len=100),parameter :: NamePlotVar='S PA B PA g r'
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='Region of existance'
    character(len=5) :: TypePlot='ascii'
    integer :: iAngle,iAngleDn,iPoint,iIono,iPlas
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nPoint,2*nAngle),PlotState_IIV(nPoint,2*nAngle,nVar))

    rPlanetCM=rPlanet_I(Planet_)*cMtoCM
!    do iLine=1,nLine
       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iPoint=1,nPoint
          do iAngle=1,2*nAngle
             !set coord based on up or down region
             if(iAngle<nAngle) then
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= EqAngleGrid_IG(iLine,iAngle)
             else
                iAngleDn = 2*nAngle-iAngle
                Coord_DII(S_,iPoint,iAngle) = FieldLineGrid_IC(iLine,iPoint)&
                     /rPlanetCM
                Coord_DII(PA_,iPoint,iAngle)= &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
             endif
             
             !set plotstate based on up or down region  
             if (iAngle <= nThetaAlt_IIC(iLine,iEnergy,iPoint))then
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = EqAngleGrid_IG(iLine,iAngle)
             elseif(iAngle > 2*nAngle-nThetaAlt_IIC(iLine,iEnergy,iPoint)) then
                iAngleDn = 2*nAngle-iAngle
                PlotState_IIV(iPoint,iAngle,B_)  = Bfield_IC(iLine,iPoint)
                PlotState_IIV(iPoint,iAngle,PA_) = &
                     cPi-EqAngleGrid_IG(iLine,iAngleDn)
             else
                PlotState_IIV(iPoint,iAngle,:)=0.0
             endif
          enddo
       enddo
       
      ! set name for plotfile
       write(NamePlot,"(a,i4.4,a,i4.4,a)") 'RegionOfExist_Enum_',iEnergy,'_line_',iLine,'.out'
 
       !Plot grid for given line
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_DII=Coord_DII,                &
            VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/))
 !   end do
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_grid_pot

  !=============================================================================
  ! update the grid for a given line when it's position moves or the 
  ! potential changes
  subroutine update_grid(iLine,Coord_D,DoOnlySpatial)
    use ModNumConst,    ONLY: cDegToRad
    implicit none
    integer, intent(in) :: iLine ! line whose grid is being updated
    real,    intent(in) :: Coord_D(2) ! foot point Lat and Lon of line
    ! When DoOnlySpatial is included then leave the potential update to later
    ! only set the spatial grid. This is important for PW coupling
    logical, optional,intent(in) :: DoOnlySpatial
    logical :: DoOnlySpatial1
    integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID
    !---------------------------------------------------------------------------
    if(.not.present(DoOnlySpatial)) then 
       DoOnlySpatial1=.false.
    else
       DoOnlySpatial1=DoOnlySpatial
    endif
    ! Set the Lshell for each line based on the input latitude
    Lshell_I(iLine) = (cos(Coord_D(Lat_)*cDegToRad))**-2.0
        
    if(IsVerbose) write(*,*) 'calling init_se_grid'
    call init_se_grid(iLine)
    
    ! when including a potential a new grid is needed
    if (DoIncludePotential .and. .not.DoOnlySpatial1) then
       call set_grid_pot(iLine)
    endif
  end subroutine update_grid
  !=============================================================================
  

  !============================================================================
  ! Test Grid: This is a routine that creates the test grid for all unit tests 
  ! (since all other modules rely on the grid)
  subroutine create_se_test_grid
    use ModPlanetConst, ONLY: Planet_, NamePlanet_I

    select case(NamePlanet_I(Planet_))

    case('EARTH')

    DrIono1 = 1e6
    nIono1  = 12
    DrIono2 = 2e6
    nIono2  = 10
    DrIono3 = 3e6
    nIono3  = 10
    DrIono4 = 5e6
    nIono4  = 2    

    case('JUPITER')

    DrIono1 = 5.0*1e6
    nIono1  = 12
    DrIono2 = 5.0*2e6
    nIono2  = 10
    DrIono3 = 5.0*3e6
    nIono3  = 10
    DrIono4 = 5.0*5e6
    nIono4  = 2    


    end select



    nIono=nIono1+nIono2+nIono3+nIono4
    ! set the energy parameters for the energy grid
    TypeGridE = 'ConstDE'
    nEnergy=800
!    nEnergy=94
    EnergyMax=100.5
!    nEnergy=2000
!    EnergyMax=2000.5
    DeltaE = 0.125

    ! Allocated the grid arrays and populate the bfield, sgrid, and PA grid  
    if(IsVerbose) write(*,*) 'allocating arrays'
    call allocate_grid_arrays

    Lshell_I(1)=10.0
    nTheta_II(1,1)=5
    nTheta_II(1,2)=20
    nTheta_II(1,3)=90
    nTheta_II(1,4)=20


    nAngle = sum(nTheta_II(1,:))

    if(IsVerbose) write(*,*) 'calling init_se_grid'
    call init_se_grid(1)
    
  end subroutine create_se_test_grid
  !============================================================================
  ! UNIT test for SE grid
  subroutine se_grid_test
    call create_se_test_grid
    write(*,*) 'calling plot_grid'
    call plot_grid
  end subroutine se_grid_test

end Module ModSeGrid
