Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_C(:),eThermalTemp_C(:)

  real, public :: mLat, mLon, gLat, gLon


  real, public :: UT = 43200.0 ! default at noon
  integer,public      :: Idate=97046 !day in the form of YYDDD

  ! when the dipole and rotation axis are aligned
  logical, public      :: DoAlignDipoleRot = .false.
  
  ! need to know if we are coupled to PWOM when setting background
  logical, public :: DoUsePWOM = .false.

  ! dip inclination angle in radians
  real, public :: dip

  ! variables for thinning the topside ionosophere 
  real :: facn=-1, fact=0

  ! exponent for extending solution above IRI or PWOM solution 
!  real,public :: ZEP=1

  ! Neutral Atmosphere arrays and variables
  integer, public :: nNeutralSpecies
  real, allocatable,public  :: NeutralDens_IC(:,:)
  real, allocatable  :: NeutralTemp_C(:)
  !earth
  integer,parameter  :: O_=1, O2_=2, N2_=3 
  !jupiter
  integer,parameter  :: H2_=1, He_=2, H_=3, CH4_=4, T_=5 
  
  ! Arrays that hold the photo electron production spectrum 
  real, allocatable,public :: ePhotoProdSpec_IC(:,:)
  real, allocatable,public :: PhotoIonRate_IC(:,:)

  !number of ions for ionization rate array
  integer :: nIons

  public :: allocate_background_arrays
  public :: fill_thermal_plasma_empirical
  public :: set_footpoint_locations
  public :: get_neutrals_and_pe_spectrum
  public :: plot_background
  public :: plot_ephoto_prod
  public :: plot_ionization_rate
contains
  !subroutines to fill in the neutral atmosphere and thermal plasma
  !=============================================================================
  subroutine set_footpoint_locations
    !---------------------------------------------------------------------------

    !  Find the geographic coordinates for our geomagnetic coordinates iono1
    CALL GEOMAG(1,gLon,gLat,mLon,mLat)
    IF (DoAlignDipoleRot) THEN
       gLat=mLat      ! Use these two lines if you want the
       gLon=mLon      ! magnetic and geographic poles aligned
    END IF
    
  end subroutine set_footpoint_locations
  !=============================================================================
  subroutine fill_thermal_plasma_empirical(F107,F107A,t)
    use ModSeGrid, only: nAlt,Alt_C,IsVerbose
    use ModPlanetConst, only: Planet_, NamePlanet_I

    real   , intent(in) :: F107, F107A, t
    real    :: factor
    integer :: iAlt
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !  Use IRI to fill the thermal plasma
       if(IsVerbose) write(*,*) 'calling get_iri'
       CALL get_iri(F107A)
       if(IsVerbose) write(*,*) 'finish get_iri'

    case('JUPITER')
       do iAlt=1,nAlt
          eThermalDensity_C(iAlt) = &
               Ne_Jupiter(Alt_C(iAlt))
          eThermalTemp_C(iAlt) = &
               eTemp_Jupiter(Alt_C(iAlt))
       end do
    end select
    
    
  end subroutine fill_thermal_plasma_empirical
  !=============================================================================
  !*  Subroutine get_iri calls IRI-90 for each ionosphere.
  SUBROUTINE get_iri(F107A)
    use ModSeGrid, only: nAlt,Alt_C,IsVerbose
    use ModNumConst,    only: cDegToRad    
    real,    intent(in) :: F107A

    real    :: RZ12 ! parameter in MSIS equal to negative of F107A
    integer :: MMDD ! -day of year
    integer :: i,j  !Generic loop indices
    real    :: STL1, STL2 ! Solar Local Time in iono 1 or 2
    
    ! IRI has output for a number of variables for a provided alt array
    real,allocatable   :: IriOutput_VC(:,:)
    integer, parameter :: nIriOutputs = 11
    real               :: OARR(30) !additional IRI outputs
    integer :: JMAG
    LOGICAL JF(12)
 
    ! unit conversion parameters
    real, parameter :: KtoeV=8.6149E-5, PerM3toPerCm3=1.E-6
    !---------------------------------------------------------------------------

    ! Set IRI output array
    if(.not.allocated(IriOutput_VC)) allocate(IriOutput_VC(nIriOutputs,nAlt))

    !  Set up IRI inputs
    do i=1,12
       JF(i)=.TRUE.
    end do
    JF(4)=.FALSE.
    JF(5)=.FALSE.
    JMAG=0
    RZ12=-F107A
    MMDD=-(Idate-(Idate/1000)*1000)
    
    !  Calculate the local solar time
    STL1=(UT/240+gLon)/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    if(IsVerbose) write(*,*) 'calling iri'
    CALL IRI90(JF,JMAG,gLat,gLon,RZ12,MMDD,STL1, &
         Alt_C(1:nAlt)/1e5,nAlt,'PW/IRI_DATA/ ',IriOutput_VC,OARR)
    !save dip inclination angle in radians from IRI output
    dip=OARR(25)*cDegToRad
    if(IsVerbose) write(*,*) 'finish iri'
    do i=nAlt,1,-1
       eThermalDensity_C(i)=IriOutput_VC(1,i)*PerM3toPerCm3 
       IF (IRIOUTPUT_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_C(i)=IriOutput_VC(4,i)*KtoeV
    end do
    
  end SUBROUTINE get_iri

  !============================================================================
  ! subroutine that fills the neutral atmosphere and PE production spectrum
  subroutine get_neutrals_and_pe_spectrum(F107,F107A,AP)
    use ModSeGrid,      only: nAlt,nEnergy,Alt_C,IsVerbose
    use ModSeProduction,only: RCOLUM,RCOLUM_ABOVE,EPHOTO,SOLZEN,SSFLUX,&
                              init_production
    use ModSeCross,     only: EXSECT,cross_jupiter
!    use ModSeCross,     only: CROSS,cross_jupiter
    use EUA_ModMsis90,  only: GTD6,TSELEC
    use ModNumConst,    only: cDegToRad,cRadToDeg
    use ModPlanetConst, only: Planet_, NamePlanet_I

    real   , intent(in) :: F107, F107A,AP(7)
    
    integer :: iAlt, iEnergy ! loop variables
    
    real    :: STL ! Solar Local Time 
    ! production variables
    real    :: SZA 
!    real,allocatable :: ColumnDens_IC(:,:)

    real, parameter :: cCmToKm=1.0e-5

    !MSIS variables
    integer,parameter :: msisO_=2, msisO2_=4, msisN2_=3 
    real :: SW(25),DN(8),TN(2)
    DATA sw/8*1.,-1.,16*1./
    
    character(len=100) :: NeutralFile

    logical,save :: IsFirstCall = .true.
    !--------------------------------------------------------------------------

    !set the solar flux
    CALL SSFLUX(0,F107,F107A,0.,0.,0.,0.,1.)


    ! on first call initialize the production parameters
    if(IsFirstCall) then
       call init_production
       IsFirstCall = .false.
    endif

    !  Calculate the local solar time
    STL=(UT/240+gLon)/15
    IF (STL.LT.0.) STL=STL+24.
    IF (STL.GT.24.) STL=STL-24.
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !  Call MSIS to get the neutral densities and temperature
       CALL TSELEC(SW)
       
       do  iAlt=1,nAlt
          CALL GTD6(Idate,UT,Alt_C(iAlt)/1e5, &
               gLat,gLon,STL,F107A,F107,AP,48,DN,TN)
          NeutralDens_IC(O_ ,iAlt)=DN(msisO_)
          NeutralDens_IC(O2_,iAlt)=DN(msisO2_)
          NeutralDens_IC(N2_,iAlt)=DN(msisN2_)
          NeutralTemp_C (iAlt)=TN(2)
       end do

       ! Calculate the solar zenith angle
       CALL SOLZEN(Idate,UT,gLat,gLon,SZA)
       SZA=SZA*cDegToRad
    case('JUPITER')
       ! interpolate from AtmosArray(1,:) to FieldLineGrid_IC(iLine,iAlt)
       NeutralFile = 'PW/JGITM-1D-atmos.dat'
       CALL get_jupiter_atmos(NeutralFile,nNeutralSpecies, &
            NeutralDens_IC(:,1:nAlt),NeutralTemp_C(1:nAlt))

       ! Calculate the solar zenith angle
       SZA=acos(cos(gLat*cDegToRad)*cos(gLon*cDegToRad))
    end select
    if (IsVerbose) write(*,*) 'SZA ',SZA*cRadToDeg
    
    !  Set the slant path column densities for O,O2 and N2
    CALL RCOLUM(SZA,Alt_C(1:nAlt), &
         NeutralDens_IC(:,:),NeutralTemp_C(:),nAlt)    

    ! for starlight input get the column density above 
    CALL RCOLUM_ABOVE(Alt_C(1:nAlt), &
         NeutralDens_IC(:,:),NeutralTemp_C(:),nAlt)    
    
    !  Calculate the photoelectron production spectrum
    CALL EPHOTO(NeutralDens_IC(:,:),ePhotoProdSpec_IC(:,:),&
         nAlt,0,SZA*cRadToDeg,Alt_C(1:nAlt)*cCmToKm,&
         nIons,PhotoIonRate_IC(:,:))

    ! set the cross sections (perhaps this should only be called once?)
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       call EXSECT
    case('JUPITER')
       call cross_jupiter(nNeutralSpecies)
    end select
    


  end subroutine get_neutrals_and_pe_spectrum
  !============================================================================
  
  subroutine plot_background(nStep,time)
    use ModPlanetConst, only: Planet_, NamePlanet_I
    use ModSeGrid,     ONLY: Alt_C,nAlt,rPlanetCM
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg

    integer, intent(in) :: nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    integer, parameter :: nDim =1, eDens_=1,eTemp_=2, &
                          VarO_=3,VarO2_=4,VarN2_=5
    integer :: nVar,iVar
    character(len=100) :: NamePlotVar
    character(len=100) :: NamePlot
    character(len=*),parameter :: NameHeader='background output'
    character(len=5) :: TypePlot='ascii'
    integer :: iPoint
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar='S ne te nO nO2 nN2 g r'
       nVar = 5                 ! neutral species + 2
    case('JUPITER')
       NamePlotVar='S ne te nH2 nHe nH nCH4 g r'
       nVar = 6                 ! neutral species + 2
    end select


    allocate(Coord_I(nAlt), PlotState_IV(nAlt,nVar))
    
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    !Set Coordinates along field line and PA
    do iPoint=1,nAlt
       Coord_I(iPoint) = Alt_C(iPoint)/1.0e5
       PlotState_IV(iPoint,eDens_) = eThermalDensity_C(iPoint)
       PlotState_IV(iPoint,eTemp_) = eThermalTemp_C(iPoint)
       do iVar=3,nVar
          PlotState_IV(iPoint,iVar) = &
               NeutralDens_IC(iVar-2,iPoint)
       enddo
    enddo
    
    ! set name for plotfile
!    write(NamePlot,"(a,i4.4,a)") 'background_iLine',iLine,'.out'
    write(NamePlot,"(a,i4.4,a)") 'background.out'
    
    !Plot grid for given line. Overwrite old results on firstcall
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) 
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) 
    end if
     
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_background


  !============================================================================
  ! save state plot for verification
  subroutine plot_ephoto_prod(nStep,time)
    use ModSeGrid,     ONLY: Alt_C,nEnergy, nAlt, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, nVar=1, E_=1, S_=2
    integer, parameter :: spec_=1

    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] eProd[cm-3eV-1s-1sr-1] g r'
    character(len=*),parameter :: NameHeader='ePhoto Production Spectrum output'
    character(len=5) :: TypePlot='ascii'
    integer :: iEnergy,iAlt
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nEnergy,nAlt),PlotState_IIV(nEnergy,nAlt,nVar))


       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iEnergy=1,nEnergy
          do iAlt=1,nAlt
             Coord_DII(E_,iEnergy,iAlt) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iAlt) = Alt_C(iAlt)/1e5
             PlotState_IIV(iEnergy,iAlt,spec_)  = &
                  ePhotoProdSpec_IC(iEnergy,iAlt)&
                  /4.0/cPi/DeltaE_I(iEnergy)
          enddo
       enddo

       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'ephotoprod.out'
       
       !Plot grid for given line
       if(IsFirstCall) then
          call save_plot_file(NamePlot, TypePositionIn='rewind', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/)) !***
          IsFirstCall = .false.
       else
          call save_plot_file(NamePlot, TypePositionIn='append', &
               TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
               NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
               nDimIn=nDim,CoordIn_DII=Coord_DII,                &
               VarIn_IIV = PlotState_IIV, ParamIn_I = (/1.6, 1.0/)) !***
       endif
    
    deallocate(Coord_DII, PlotState_IIV)
  end subroutine plot_ephoto_prod
  !============================================================================
  ! save state plot for verification
  subroutine plot_ionization_rate(nStep,time)
    use ModSeGrid,     ONLY: Alt_C, nAlt, nEnergy, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst, only: Planet_, NamePlanet_I

    integer, intent(in) :: nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_I(:), PlotState_IV(:,:)
    !grid parameters
    integer, parameter :: nDim =1,S_=1
    integer :: nVar

    !Jupiter
    integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
         CH3plus_=5,CH2plus_=6,CHplus_=7
    
    !Earth
    integer, parameter :: Oplus_=1
    
    character(len=100),parameter :: NamePlotVarEarth=&
         'Alt[km] O+[cm-3s-1] g r'
    character(len=100),parameter :: NamePlotVarJupiter=&
         'Alt[km] H2+[cm-3s-1] He+[cm-3s-1] H+[cm-3s-1] CH4+[cm-3s-1] CH3+[cm-3s-1] CH2+[cm-3s-1] CH+[cm-3s-1] g r'

    character(len=100) :: NamePlotVar
    character(len=*),parameter :: NameHeader='Photoionization Rates'
    character(len=5) :: TypePlot='ascii'
    integer :: iIon,iAlt
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------

    nVar=nIons!+2
    allocate(Coord_I(nAlt),PlotState_IV(nAlt,nVar))
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar=NamePlotVarEarth
    case('JUPITER')
       NamePlotVar=NamePlotVarJupiter
    end select
       

    !Set Coordinates along field line and PA
    do iAlt=1,nAlt
       Coord_I(iAlt) = Alt_C(iAlt)/1e5

       select case(NamePlanet_I(Planet_))
       case('EARTH')
          PlotState_IV(iAlt,Oplus_)  = &
               PhotoIonRate_IC(Oplus_,iAlt)
       case('JUPITER')
          PlotState_IV(iAlt,H2plus_)  = &
               PhotoIonRate_IC(H2plus_,iAlt)
          PlotState_IV(iAlt,Heplus_)  = &
               PhotoIonRate_IC(Heplus_,iAlt)
          PlotState_IV(iAlt,Hplus_)  = &
               PhotoIonRate_IC(Hplus_,iAlt)
          PlotState_IV(iAlt,CH4plus_)  = &
               PhotoIonRate_IC(CH4plus_,iAlt)
          PlotState_IV(iAlt,CH3plus_)  = &
               PhotoIonRate_IC(CH3plus_,iAlt)
          PlotState_IV(iAlt,CH2plus_)  = &
               PhotoIonRate_IC(CH2plus_,iAlt)
          PlotState_IV(iAlt,CHplus_)  = &
               PhotoIonRate_IC(CHplus_,iAlt)
!          PlotState_IV(iAlt,8)  = &
!               NeutralDens_IC(iLine,H2_,iAlt)
!          PlotState_IV(iAlt,9)  = &
!               NeutralDens_IC(iLine,H_,iAlt)
       end select
       
       
    enddo

    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'PhotoIonization.out'
    
    !Plot grid for given line
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn=nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
    endif
    
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_ionization_rate
  
  !=============================================================================
  subroutine allocate_background_arrays
    use ModSeGrid,     ONLY: nAlt, nEnergy
    use ModPlanetConst, only: Planet_, NamePlanet_I
    !---------------------------------------------------------------------------

    select case(NamePlanet_I(Planet_))
    case('EARTH')
       nNeutralSpecies=3
       nIons=1
    case('JUPITER')
       nNeutralSpecies=4
       nIons=7
    end select

    if(.not.allocated(eThermalDensity_C)) &
         allocate(eThermalDensity_C(nAlt))
    if(.not.allocated(eThermalTemp_C)) &
         allocate(eThermalTemp_C(nAlt))

    if(.not.allocated(NeutralDens_IC)) &
         allocate(NeutralDens_IC(nNeutralSpecies,nAlt))

    if(.not.allocated(NeutralTemp_C)) &
         allocate(NeutralTemp_C(nAlt))

    if(.not.allocated(ePhotoProdSpec_IC)) &
         allocate(ePhotoProdSpec_IC(nEnergy,nAlt))

    if(.not.allocated(PhotoIonRate_IC)) &
         allocate(PhotoIonRate_IC(nIons,nAlt))

  end subroutine allocate_background_arrays

  !============================================================================
  
  subroutine get_jupiter_atmos(DatafileName,nSpecies,NeutralDens_IC,NeutralTemp_C)
    use ModInterpolate, ONLY: linear
    use ModSeGrid, ONLY : Alt_C, nAlt
    use ModIoUnit, ONLY : UnitTmp_
    character (len=100), intent(in) :: DatafileName ! 'JGITM-1D-atmos.dat'
    integer,            intent(in)  :: nSpecies     ! 4
    real,               intent(out) :: NeutralDens_IC(nSpecies,nAlt)
    real,               intent(out) :: NeutralTemp_C(nAlt)
    
    integer, parameter :: nAltGrid= 10000
    integer            :: iAlt
    character (len=189) :: line1
    real :: AtmosArray(5+nSpecies,nAltGrid)
    
    write(*,*) 'starting get_jupiter_atmos'
    open(UnitTmp_,FILE=DatafileName,STATUS='OLD')
    
    read(UnitTmp_,'(a)') line1
    read(UnitTmp_,*) AtmosArray
    
    close(UnitTmp_)
    
    print *,line1
    print *, AtmosArray(1,1),AtmosArray(2,1),AtmosArray(9,1)
    print *, AtmosArray(1,10000),AtmosArray(2,10000),AtmosArray(9,10000)
    write(*,*) 'starting get_jupiter_atmos linear'
    do iAlt=1,nAlt
       NeutralDens_IC(1,iAlt) = linear(AtmosArray(5,:), &        ! H2
            1,nAltGrid,Alt_C(iAlt)/1e5,AtmosArray(1,:))
       NeutralDens_IC(2,iAlt) = linear(AtmosArray(6,:), &        ! He
            1,nAltGrid,Alt_C(iAlt)/1e5,AtmosArray(1,:))
       NeutralDens_IC(3,iAlt) = linear(AtmosArray(7,:), &        ! H
            1,nAltGrid,Alt_C(iAlt)/1e5,AtmosArray(1,:))
       NeutralDens_IC(4,iAlt) = linear(AtmosArray(8,:), &        ! CH4
            1,nAltGrid,Alt_C(iAlt)/1e5,AtmosArray(1,:))
       NeutralTemp_C(iAlt) = linear(AtmosArray(2,:), &        ! Temp
            1,nAltGrid,Alt_C(iAlt)/1e5,AtmosArray(1,:))
    end do
    write(*,*) 'done get_jupiter_atmos'
  end subroutine get_jupiter_atmos
  !=============================================================================
  function Ne_Jupiter(z)
    use ModSeGrid, ONLY: rPlanetCM
    real, intent(in) :: z
    ! fit based on Kitamura, [2011]
    ! parameters for fit choosen to match data from Yelle and Miller, [2004]
    real, parameter :: n600 = 3.5e5
    real, parameter :: n3500 = 3.e3
    real :: h600, alphav
    real :: r, eqn1, eqn2, topside, ne_cutoff
    real :: Ne_Jupiter
    !-------------------------------------------------------------------------
    h600 = 175.
    alphav = 70.
    
    r = 1. + z/rPlanetCM
    eqn1 = n600*exp(1.01*(600.-z)/(r*h600))
    eqn2 = n3500*(r/1.05)**(-alphav)
    
    topside = eqn1 + eqn2
    ! cutoff at z = 600 km
    ne_cutoff = n600 + n3500*((1.+600./rPlanetCM)/1.05)**(-alphav)
    
    Ne_Jupiter = min(topside,z*ne_cutoff/600.)
    
    return
    
  end function Ne_Jupiter
  !=============================================================================
  function eTemp_Jupiter(z)
    use ModConst, ONLY: cBoltzmann
    real, intent(in)  :: z
    real :: eTemp_Jupiter
    real :: eTemp_func
    real, parameter :: cJtoeV = 6.242e+18

    ! fit choosen to match data from Yelle and Miller, [2004]
    eTemp_func = 25.*sqrt(max(1.e-10,z-300.))+200.
    eTemp_Jupiter = min(eTemp_func,900.) *cBoltzmann* cJtoeV   
       
    return

  end function eTemp_Jupiter
end Module ModSeBackground
