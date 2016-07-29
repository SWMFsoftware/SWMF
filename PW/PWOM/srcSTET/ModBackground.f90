Module ModSeBackground
  implicit none

  private !except

  real, public,allocatable :: eThermalDensity_IC(:,:),eThermalTemp_IC(:,:)

  real, public,allocatable :: mLat_I(:), mLon_I(:), gLat1_I(:), gLon1_I(:), &
                              gLat2_I(:), gLon2_I(:)


  real, public :: UT = 43200.0 ! default at noon
  integer,public      :: Idate=97046 !day in the form of YYDDD

  ! when the dipole and rotation axis are aligned
  logical, public      :: DoAlignDipoleRot = .false.
  
  ! need to know if we are coupled to PWOM when setting background
  logical, public :: DoUsePWOM = .false.


  ! variables for thinning the topside ionosophere 
  real :: facn=-1, fact=0

  ! exponent for extending solution above IRI or PWOM solution 
  real,public :: ZEP=1

  ! Neutral Atmosphere arrays and variables
  integer, public :: nNeutralSpecies
  real, allocatable,public  :: NeutralDens1_IIC(:,:,:),NeutralDens2_IIC(:,:,:)
  real, allocatable  :: NeutralTemp1_IC(:,:),NeutralTemp2_IC(:,:)
  !earth
  integer,parameter  :: O_=1, O2_=2, N2_=3 
  !jupiter
  integer,parameter  :: H2_=1, He_=2, H_=3, CH4_=4, T_=5 
  
  ! Logical variables for if we are calculating photo e spectrum in iono 1 or 2
  logical :: DoCalcPeIono1=.true., DoCalcPeIono2=.true.

  ! Arrays that hold the photo electron production spectrum in iono 1 or 2
  real, allocatable,public :: ePhotoProdSpec1_IIC(:,:,:),ePhotoProdSpec2_IIC(:,:,:)
  real, allocatable,public :: PhotoIonRate1_IIC(:,:,:),PhotoIonRate2_IIC(:,:,:)

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
  subroutine set_footpoint_locations(iLine)
    integer, intent(in) :: iLine
    !---------------------------------------------------------------------------

    !  Find the geographic coordinates for our geomagnetic coordinates iono1
    CALL GEOMAG(1,gLon1_I(iLine),gLat1_I(iLine),mLon_I(iLine),mLat_I(iLine))
    IF (DoAlignDipoleRot) THEN
       gLat1_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want the
       gLon1_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligned
    END IF
    
    !  Find the geographic coordinates for our geomagnetic coordinates iono2
    CALL GEOMAG(1,gLon2_I(iLine),gLat2_I(iLine),mLon_I(iLine),-mLat_I(iLine))
    IF (DoAlignDipoleRot) THEN
       gLat2_I(iLine)=mLat_I(iLine)      ! Use these two lines if you want t
       gLon2_I(iLine)=mLon_I(iLine)      ! magnetic and geographic poles aligne
    END IF
  end subroutine set_footpoint_locations
  !=============================================================================
  subroutine fill_thermal_plasma_empirical(iLine,F107,F107A,t)
    use ModSeGrid, only: nIono1,nIono2,nIono,nPlas, nPoint, &
                         FieldLineGrid_IC,Bfield_IC,IsVerbose
    use ModPlanetConst, only: Planet_, NamePlanet_I

    integer, intent(in) :: iLine
    real   , intent(in) :: F107, F107A, t
    integer :: iIono,iIono2,iPlas,nTopIono1,nTopIono2
    real    :: factor

    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !  Use IRI to fill the thermal plasma
       if(IsVerbose) write(*,*) 'calling get_iri'
       CALL get_iri(iLine,F107A)
       if(IsVerbose) write(*,*) 'finish get_iri'
       ! The next few lines are for thinning the topside ionosphere densities
       IF (ABS(facn).GT.0.) THEN
          do iIono=nIono1+nIono2+1,nIono
             iIono2=nPoint-iIono+1
             factor= (FieldLineGrid_IC(iLine,iIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))&
                  /(FieldLineGrid_IC(iLine,nIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))
             factor=10**(ALOG10(ABS(facn))*factor)
             eThermalDensity_IC(iLine,iIono)=eThermalDensity_IC(iLine,iIono) &
                  /factor
             eThermalDensity_IC(iLine,iIono2)=&
                  eThermalDensity_IC(iLine,iIono2)/factor
          end do
       END IF
       
       ! Alter the topside ionospheric thermal temperatures
       IF (fact.NE.0.) THEN
          do iIono=nIono1+nIono2+1,nIono
             iIono2=nPoint-iIono+1
             factor= (FieldLineGrid_IC(iLine,iIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))&
                  /(FieldLineGrid_IC(iLine,nIono) &
                  - FieldLineGrid_IC(iLine,nIono1+nIono2))
             
             eThermalTemp_IC(iLine,iIono) = eThermalTemp_IC(iLine,iIono) &
                  + factor*(fact-eThermalTemp_IC(iLine,iIono))
             eThermalTemp_IC(iLine,iIono2) = eThermalTemp_IC(iLine,iIono2) &
                  + factor*(fact-eThermalTemp_IC(iLine,iIono2))
          end do
       END IF
       ! End IRI setup

    case('JUPITER')
       do iIono=1,nIono
          !fill in iono1
          eThermalDensity_IC(iLine,iIono) = &
               Ne_Jupiter(FieldLineGrid_IC(iLine,iIono))
          eThermalTemp_IC(iLine,iIono) = &
               eTemp_Jupiter(FieldLineGrid_IC(iLine,iIono))
          !fill in iono2
          iIono2=nPoint-iIono+1
          eThermalDensity_IC(iLine,iIono2) = eThermalDensity_IC(iLine,iIono)
          eThermalTemp_IC(iLine,iIono2)    = eThermalTemp_IC(iLine,iIono)
       end do
    end select
    
    ! Fill in the plasmaspheric thermal densities
    nTopIono1 = nIono		! Simplifying notation, not en. index
    nTopIono2 = nIono+nPlas+1	! Simplifying notation, not angle index
    
    do iPlas=nTopIono1+1,nTopIono2-1
       factor=&
            ( FieldLineGrid_IC(iLine,iPlas) &
            - FieldLineGrid_IC(iLine,nTopIono1) ) &
            / ( FieldLineGrid_IC(iLine,nTopIono2) &
            - FieldLineGrid_IC(iLine,nTopIono1) )
       eThermalDensity_IC(iLine,iPlas) = &
            ((1.-factor)*eThermalDensity_IC(iLine,nTopIono1) + &
            factor*eThermalDensity_IC(iLine,nTopIono2))        &
            * (Bfield_IC(iLine,iPlas)/Bfield_IC(iLine,nTopIono1))**ZEP
    end do
    
    
    ! Fill in the plasmaspheric thermal temperatures
    do iPlas=nTopIono1+1,nTopIono2-1
       factor=&
            ( FieldLineGrid_IC(iLine,iPlas) &
            - FieldLineGrid_IC(iLine,nTopIono1) ) &
            / ( FieldLineGrid_IC(iLine,nTopIono2)&
            - FieldLineGrid_IC(iLine,nTopIono1) )
       eThermalTemp_IC(iLine,iPlas)= & 
            ((1.-factor)*eThermalTemp_IC(iLine,nTopIono1) &
            + factor*eThermalTemp_IC(iLine,nTopIono2))
    end do
    
  end subroutine fill_thermal_plasma_empirical
  !=============================================================================
  !*  Subroutine get_iri calls IRI-90 for each ionosphere.
  SUBROUTINE get_iri(iLine,F107A)
    use ModSeGrid, only: nIono,nPoint,FieldLineGrid_IC,IsVerbose
    use EUA_ModIri90, only: iri90

    real,    intent(in) :: F107A
    integer, intent(in) :: iLine

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
    if(.not.allocated(IriOutput_VC)) allocate(IriOutput_VC(nIriOutputs,nIono))

    !  Set up IRI inputs
    do i=1,12
       JF(i)=.TRUE.
    end do
    JF(4)=.FALSE.
    JF(5)=.FALSE.
    JF(12)=.FALSE.
    JMAG=0
    RZ12=-F107A
    MMDD=-(Idate-(Idate/1000)*1000)
    
    !\
    ! Work on first ionosphere
    !/
 
    !  Calculate the local solar time
    STL1=(UT/240+gLon1_I(iLine))/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    if(IsVerbose) write(*,*) 'calling iri'
    CALL IRI90(JF,JMAG,gLat1_I(iLine),gLon1_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,1:nIono)/1e5,nIono, &
         'PW/IRI_DATA/ccir.cofcnts', &
         'PW/IRI_DATA/ursi.cofcnts', IriOutput_VC,OARR, 0)
    if(IsVerbose) write(*,*) 'finish iri'
    do i=nIono,1,-1
       eThermalDensity_IC(iLine,i)=IriOutput_VC(1,i)*PerM3toPerCm3 
       IF (IRIOUTPUT_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,i)=IriOutput_VC(4,i)*KtoeV
    end do
    
    !\
    ! Work on second ionosphere
    !/
    
    !  Calculate the local solar time
    STL2=(UT/240+gLon2_I(iLine))/15
    IF (STL2.LT.0.) STL2=STL2+24.
    IF (STL2.GT.24.) STL2=STL2-24.
    
    !  Call IRI for the second ionosphere
    CALL IRI90(JF,JMAG,gLat2_I(iLine),gLon2_I(iLine),RZ12,MMDD,STL1, &
         FieldLineGrid_IC(iLine,1:nIono)/1e5,nIono, &
         'PW/IRI_DATA/ccir.cofcnts', &
         'PW/IRI_DATA/ursi.cofcnts', IriOutput_VC, OARR, 0)
    
    do i=nIono,1,-1
       j=nPoint-i+1
       eThermalDensity_IC(iLine,j)=IriOutput_VC(1,i)*PerM3toPerCm3
       IF (IriOutput_VC(4,i).LT.0.) IriOutput_VC(4,i)=IriOutput_VC(4,i+1)
       eThermalTemp_IC(iLine,j)=IriOutput_VC(4,i)*KtoeV
    end do
    
  end SUBROUTINE get_iri

  !============================================================================
  ! subroutine that fills the neutral atmosphere and PE production spectrum
  subroutine get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP)
    use ModSeGrid,      only: nIono,nEnergy,nPoint,FieldLineGrid_IC,IsVerbose
    use ModSeProduction,only: RCOLUM,RCOLUM_ABOVE,ESPEC,SOLZEN,SSFLUX,&
                              init_production
    use ModSeCross,     only: cross,cross_jupiter
    use EUA_ModMsis90,  only: GTD6,TSELEC
    use ModNumConst,    only: cDegToRad,cRadToDeg
    use ModPlanetConst, only: Planet_, NamePlanet_I

    integer, intent(in) :: iLine
    real   , intent(in) :: F107, F107A,AP(7)
    
    integer :: iIono, iEnergy ! loop variables
    
    real    :: STL1, STL2 ! Solar Local Time in iono 1 or 2
    ! production variables
    real    :: SZA1,SZA2
    real,allocatable :: ColumnDens_IC(:,:)

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


    !\
    ! Work on ionosphere 1
    !/
    !  Calculate the local solar time
    STL1=(UT/240+gLon1_I(iLine))/15
    IF (STL1.LT.0.) STL1=STL1+24.
    IF (STL1.GT.24.) STL1=STL1-24.
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !  Call MSIS to get the neutral densities and temperature
       CALL TSELEC(SW)
       
       do  iIono=1,nIono
          CALL GTD6(Idate,UT,FieldLineGrid_IC(iLine,iIono)/1e5, &
               gLat1_I(iLine),gLon1_I(iLine),STL1,F107A,F107,AP,48,DN,TN)
          NeutralDens1_IIC(iLine,O_ ,iIono)=DN(msisO_)
          NeutralDens1_IIC(iLine,O2_,iIono)=DN(msisO2_)
          NeutralDens1_IIC(iLine,N2_,iIono)=DN(msisN2_)
          NeutralTemp1_IC (iLine,iIono)=TN(2)
       end do

       ! Calculate the solar zenith angle
       CALL SOLZEN(Idate,UT,gLat1_I(iLine),gLon1_I(iLine),SZA1)
       SZA1=SZA1*cDegToRad
    case('JUPITER')
       ! interpolate from AtmosArray(1,:) to FieldLineGrid_IC(iLine,iAlt)
       NeutralFile = 'PW/JGITM-1D-atmos.dat'
       CALL get_jupiter_atmos(NeutralFile,nNeutralSpecies, iLine, &
            NeutralDens1_IIC(iLine,:,1:nIono),NeutralTemp1_IC(iLine,1:nIono))

       ! Calculate the solar zenith angle
       SZA1=acos(cos(gLat1_I(iLine)*cDegToRad)*cos(gLon1_I(iLine)*cDegToRad))
    end select
    if (IsVerbose) write(*,*) 'SZA1 ',SZA1*cRadToDeg
    
    !  Set the slant path column densities for O,O2 and N2
    CALL RCOLUM(SZA1,FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens1_IIC(iLine,:,:),NeutralTemp1_IC(iLine,:),nIono)    

    ! for starlight input get the column density above 
    CALL RCOLUM_ABOVE(FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens1_IIC(iLine,:,:),NeutralTemp1_IC(iLine,:),nIono)    
    
    !  Calculate the photoelectron production spectrum
!    IF ((SZA1.LT.2.).AND.(DoCalcPeIono1)) THEN
    IF (DoCalcPeIono1) THEN
       CALL ESPEC(NeutralDens1_IIC(iLine,:,:),ePhotoProdSpec1_IIC(iLine,:,:),&
            nIono,0,SZA1*cRadToDeg,FieldLineGrid_IC(iLine,1:nIono)*cCmToKm,&
            nIons,PhotoIonRate1_IIC(iLine,:,:))
    ELSE
       do iIono=1,nIono
          do iEnergy=1,nEnergy
             ePhotoProdSpec1_IIC(iLine,iEnergy,iIono)=0.
	  end do
       end do
    END IF
!    if(IsVerbose) write (*,*) PhotoIonRate1_IIC(1,10:15,10:15)

    !\
    ! Work on ionosphere 2
    !/
    !  Calculate the local solar time
    STL2=(UT/240+gLon2_I(iLine))/15
    IF (STL2.LT.0.) STL2=STL2+24.
    IF (STL2.GT.24.) STL2=STL2-24.
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       !  Call MSIS to get the neutral densities and temperature
       CALL TSELEC(SW)
       do  iIono=1,nIono
          CALL GTD6(Idate,UT,FieldLineGrid_IC(iLine,iIono)/1e5, &
               gLat2_I(iLine),gLon2_I(iLine),STL2,F107A,F107,AP,48,DN,TN)
          NeutralDens2_IIC(iLine,O_ ,iIono)=DN(msisO_)
          NeutralDens2_IIC(iLine,O2_,iIono)=DN(msisO2_)
          NeutralDens2_IIC(iLine,N2_,iIono)=DN(msisN2_)
          NeutralTemp2_IC (iLine,iIono)=TN(2)
       end do

       ! Calculate the solar zenith angle
       CALL SOLZEN(Idate,UT,gLat2_I(iLine),gLon2_I(iLine),SZA2)
       SZA2=SZA2*cDegToRad
       
    case('JUPITER')
       NeutralDens2_IIC(iLine,:,:)  = NeutralDens1_IIC(iLine,:,:)
       NeutralTemp2_IC(iLine,:) = NeutralTemp1_IC(iLine,:)

       ! Calculate the solar zenith angle
       SZA2=acos(cos(gLat2_I(iLine)*cDegToRad)*cos(gLon2_I(iLine)*cDegToRad))
    end select

    if (IsVerbose) write(*,*) 'SZA2 ',SZA2*cRadToDeg


    !  Set the slant path column densities for O,O2 and N2
    CALL RCOLUM(SZA2,FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens2_IIC(iLine,:,:),NeutralTemp2_IC(iLine,:),nIono)    

    ! for starlight input get the column density above 
    CALL RCOLUM_ABOVE(FieldLineGrid_IC(iLine,1:nIono), &
         NeutralDens1_IIC(iLine,:,:),NeutralTemp1_IC(iLine,:),nIono)    


    !  Calculate the photoelectron production spectrum
    IF (DoCalcPeIono2) THEN
       CALL ESPEC(NeutralDens2_IIC(iLine,:,:),ePhotoProdSpec2_IIC(iLine,:,:),&
            nIono,nPoint,SZA2*cRadToDeg,FieldLineGrid_IC(iLine,1:nIono)*cCmToKm,&
            nIons,PhotoIonRate2_IIC(iLine,:,:))
    ELSE
       do iIono=1,nIono
          do iEnergy=1,nEnergy
             ePhotoProdSpec2_IIC(iLine,iEnergy,iIono)=0.
	  end do
       end do
    END IF
!    if(IsVerbose) write (*,*) 'PhotoIon2: ',PhotoIonRate2_IIC(1,10:15,10:15)


    ! set the cross sections (perhaps this should only be called once?)
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       call cross
    case('JUPITER')
       call cross_jupiter(nNeutralSpecies)
    end select
    


  end subroutine get_neutrals_and_pe_spectrum
  !============================================================================
  
  subroutine plot_background(iLine,nStep,time)
    use ModPlanetConst, only: Planet_, NamePlanet_I
    use ModSeGrid,     ONLY: FieldLineGrid_IC,nLine,nPoint,nIono,rPlanetCM
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg

    integer, intent(in) :: iLine,nStep
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


    allocate(Coord_I(nPoint), PlotState_IV(nPoint,nVar))
    
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    !Set Coordinates along field line and PA
    do iPoint=1,nPoint
       Coord_I(iPoint) = FieldLineGrid_IC(iLine,iPoint)/rPlanetCM
       PlotState_IV(iPoint,eDens_) = eThermalDensity_IC(iLine,iPoint)
       PlotState_IV(iPoint,eTemp_) = eThermalTemp_IC(iLine,iPoint)
       do iVar=3,nVar
          if (iPoint <= nIono) then
!          PlotState_IV(iPoint,VarO_)  = NeutralDens1_IIC(iLine,O_,iPoint)
!          PlotState_IV(iPoint,VarO2_) = NeutralDens1_IIC(iLine,O2_,iPoint)
!          PlotState_IV(iPoint,VarN2_) = NeutralDens1_IIC(iLine,N2_,iPoint)
             PlotState_IV(iPoint,iVar) = &
                  NeutralDens1_IIC(iLine,iVar-2,iPoint)
          elseif(iPoint>nPoint-nIono) then
!          PlotState_IV(iPoint,VarO_)  = &
!               NeutralDens2_IIC(iLine,O_,nPoint-iPoint+1)
!          PlotState_IV(iPoint,VarO2_) = &
!               NeutralDens2_IIC(iLine,O2_,nPoint-iPoint+1)
!          PlotState_IV(iPoint,VarN2_) = &
!               NeutralDens2_IIC(iLine,N2_,nPoint-iPoint+1)

             PlotState_IV(iPoint,iVar)  = &
                  NeutralDens2_IIC(iLine,iVar-2,nPoint-iPoint+1)
          else
          !neutral atmosphere not considered in plasmaphere
             PlotState_IV(iPoint,iVar)  = 0.0
!          PlotState_IV(iPoint,VarO2_) = 0.0
!          PlotState_IV(iPoint,VarN2_) = 0.0
          endif
       enddo
    enddo
    
    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'background_iLine',iLine,'.out'
    
    !Plot grid for given line. Overwrite old results on firstcall
    if(IsFirstCall) then
       call save_plot_file(NamePlot, TypePositionIn='rewind', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
       IsFirstCall = .false.
    else
       call save_plot_file(NamePlot, TypePositionIn='append', &
            TypeFileIn=TypePlot,StringHeaderIn = NameHeader,  &
            NameVarIn = NamePlotVar, nStepIn= nStep,TimeIn=time,     &
            nDimIn=nDim,CoordIn_I=Coord_I,                &
            VarIn_IV = PlotState_IV, ParamIn_I = (/1.6, 1.0/)) !***
    end if
     
    deallocate(Coord_I, PlotState_IV)
  end subroutine plot_background


  !============================================================================
  !============================================================================
  ! save state plot for verification
  subroutine plot_ephoto_prod(iLine,nStep,time)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi

    integer, intent(in) :: iLine, nStep
    real,    intent(in) :: time

    real, allocatable   :: Coord_DII(:,:,:), PlotState_IIV(:,:,:)
    !grid parameters
    integer, parameter :: nDim =2, nVar=2, E_=1, S_=2
    integer, parameter :: spec1_=1, spec2_=2

    character(len=100),parameter :: NamePlotVar='E[eV] Alt[km] eProd1[cm-3eV-1s-1sr-1] eProd2[cm-3eV-1s-1sr-1] g r'
    character(len=*),parameter :: NameHeader='ePhoto Production Spectrum output'
    character(len=5) :: TypePlot='ascii'
    integer :: iEnergy,iIono
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------
    allocate(Coord_DII(nDim,nEnergy,nIono),PlotState_IIV(nEnergy,nIono,nVar))


       PlotState_IIV = 0.0
       Coord_DII     = 0.0
       
       !Set Coordinates along field line and PA
       do iEnergy=1,nEnergy
          do iIono=1,nIono
             Coord_DII(E_,iEnergy,iIono) = EnergyGrid_I(iEnergy)             
             Coord_DII(S_,iEnergy,iIono) = FieldLineGrid_IC(iLine,iIono)/1e5
             PlotState_IIV(iEnergy,iIono,spec1_)  = &
                  ePhotoProdSpec1_IIC(iLine,iEnergy,iIono)&
                  /4.0/cPi/DeltaE_I(iEnergy)
             PlotState_IIV(iEnergy,iIono,spec2_)  = &
                  ePhotoProdSpec2_IIC(iLine,iEnergy,iIono)&
                  /4.0/cPi/DeltaE_I(iEnergy)
          enddo
       enddo

       ! set name for plotfile
       write(NamePlot,"(a,i4.4,a)") 'ephotoprod_iLine',iLine,'.out'
       
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
  subroutine plot_ionization_rate(iLine,nStep,time)
    use ModSeGrid,     ONLY: FieldLineGrid_IC,nIono,nEnergy, nPoint, &
         DeltaE_I,EnergyGrid_I
    use ModIoUnit,     ONLY: UnitTmp_
    use ModPlotFile,   ONLY: save_plot_file
    use ModNumConst,   ONLY: cRadToDeg,cPi
    use ModPlanetConst, only: Planet_, NamePlanet_I

    integer, intent(in) :: iLine, nStep
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
    integer :: iIon,iIono
    character(len=100) :: NamePlot
    logical,save :: IsFirstCall =.true.
    !--------------------------------------------------------------------------

    nVar=nIons!+2
    allocate(Coord_I(nIono),PlotState_IV(nIono,nVar))
    PlotState_IV = 0.0
    Coord_I     = 0.0
    
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       NamePlotVar=NamePlotVarEarth
    case('JUPITER')
              NamePlotVar=NamePlotVarJupiter
    end select
       

    !Set Coordinates along field line and PA
    do iIono=1,nIono
       Coord_I(iIono) = FieldLineGrid_IC(iLine,iIono)/1e5

       select case(NamePlanet_I(Planet_))
       case('EARTH')
          PlotState_IV(iIono,Oplus_)  = &
               PhotoIonRate1_IIC(iLine,Oplus_,iIono)
       case('JUPITER')
          PlotState_IV(iIono,H2plus_)  = &
               PhotoIonRate1_IIC(iLine,H2plus_,iIono)
          PlotState_IV(iIono,Heplus_)  = &
               PhotoIonRate1_IIC(iLine,Heplus_,iIono)
          PlotState_IV(iIono,Hplus_)  = &
               PhotoIonRate1_IIC(iLine,Hplus_,iIono)
          PlotState_IV(iIono,CH4plus_)  = &
               PhotoIonRate1_IIC(iLine,CH4plus_,iIono)
          PlotState_IV(iIono,CH3plus_)  = &
               PhotoIonRate1_IIC(iLine,CH3plus_,iIono)
          PlotState_IV(iIono,CH2plus_)  = &
               PhotoIonRate1_IIC(iLine,CH2plus_,iIono)
          PlotState_IV(iIono,CHplus_)  = &
               PhotoIonRate1_IIC(iLine,CHplus_,iIono)
!          PlotState_IV(iIono,8)  = &
!               NeutralDens1_IIC(iLine,H2_,iIono)
!          PlotState_IV(iIono,9)  = &
!               NeutralDens1_IIC(iLine,H_,iIono)
       end select
       
       
    enddo

    ! set name for plotfile
    write(NamePlot,"(a,i4.4,a)") 'PhotoIonization_iLine',iLine,'.out'
    
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
    use ModSeGrid,     ONLY: nLine, nPoint, nIono, nEnergy
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

    if(.not.allocated(eThermalDensity_IC)) &
         allocate(eThermalDensity_IC(nLine,nPoint))
    if(.not.allocated(eThermalTemp_IC)) &
         allocate(eThermalTemp_IC(nLine,nPoint))
    if(.not.allocated(mLat_I)) &
         allocate(mLat_I(nLine))
    if(.not.allocated(mLon_I)) &
         allocate(mLon_I(nLine))
    if(.not.allocated(gLat1_I)) &
         allocate(gLat1_I(nLine))
    if(.not.allocated(gLon1_I)) &
         allocate(gLon1_I(nLine))
    if(.not.allocated(gLat2_I)) &
         allocate(gLat2_I(nLine))
    if(.not.allocated(gLon2_I)) &
         allocate(gLon2_I(nLine))

    if(.not.allocated(NeutralDens1_IIC)) &
         allocate(NeutralDens1_IIC(nLine,nNeutralSpecies,nIono))
    if(.not.allocated(NeutralDens2_IIC)) &
         allocate(NeutralDens2_IIC(nLine,nNeutralSpecies,nIono))

    if(.not.allocated(NeutralTemp1_IC)) &
         allocate(NeutralTemp1_IC(nLine,nIono))
    if(.not.allocated(NeutralTemp2_IC)) &
         allocate(NeutralTemp2_IC(nLine,nIono))

    if(.not.allocated(ePhotoProdSpec1_IIC)) &
         allocate(ePhotoProdSpec1_IIC(nLine,nEnergy,nIono))
    if(.not.allocated(ePhotoProdSpec2_IIC)) &
         allocate(ePhotoProdSpec2_IIC(nLine,nEnergy,nIono))

    if(.not.allocated(PhotoIonRate1_IIC)) &
         allocate(PhotoIonRate1_IIC(nLine,nIons,nIono))
    if(.not.allocated(PhotoIonRate2_IIC)) &
         allocate(PhotoIonRate2_IIC(nLine,nIons,nIono))


  end subroutine allocate_background_arrays
  
  subroutine get_jupiter_atmos(DatafileName,nSpecies,iLine,NeutralDens_IC,NeutralTemp_C)
    use ModInterpolate, ONLY: linear
    use ModSeGrid, ONLY : FieldLineGrid_IC, nIono
    use ModIoUnit, ONLY : UnitTmp_
    character (len=100), intent(in) :: DatafileName ! 'JGITM-1D-atmos.dat'
    integer,            intent(in) :: nSpecies     ! 4
    integer,            intent(in) :: iLine
    real,               intent(out) :: NeutralDens_IC(nSpecies,nIono)
    real,               intent(out) :: NeutralTemp_C(nIono)
    
    integer, parameter :: nAltGrid= 10000
    integer            :: iIono
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
    do iIono=1,nIono
       NeutralDens_IC(1,iIono) = linear(AtmosArray(5,:), &        ! H2
            1,nAltGrid,FieldLineGrid_IC(iLine,iIono)/1e5,AtmosArray(1,:))
       NeutralDens_IC(2,iIono) = linear(AtmosArray(6,:), &        ! He
            1,nAltGrid,FieldLineGrid_IC(iLine,iIono)/1e5,AtmosArray(1,:))
       NeutralDens_IC(3,iIono) = linear(AtmosArray(7,:), &        ! H
            1,nAltGrid,FieldLineGrid_IC(iLine,iIono)/1e5,AtmosArray(1,:))
       NeutralDens_IC(4,iIono) = linear(AtmosArray(8,:), &        ! CH4
            1,nAltGrid,FieldLineGrid_IC(iLine,iIono)/1e5,AtmosArray(1,:))
       NeutralTemp_C(iIono) = linear(AtmosArray(2,:), &        ! Temp
            1,nAltGrid,FieldLineGrid_IC(iLine,iIono)/1e5,AtmosArray(1,:))
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
