Module ModCouplePWOMtoSE
  implicit none
  
  private !except
  
  integer :: nAltPw
  real,allocatable :: AltPw_C(:)
  real :: dAltPw
  integer :: nIonPW
  !named parameters for photoelectron coupling back to PWOM
  !earth
  integer :: OplusPW_=1
  !jupiter
  integer :: H3plusPW_=1,HplusPW_=2,H2plusPW_=3
  !public methods
  public :: init_pwom_se_coupling
  public :: get_se_for_pwom
  
  integer,allocatable :: iLineGlobal_I(:)
  integer::nLine
contains
  !=============================================================================
  ! set up coupling between pwom and stet for all lines. 
  ! thermal e density and temperature for each line and interpolate to STET 
  ! grid 
  subroutine init_pwom_se_coupling(IsVerbosePw,nAltPwIn,nLinePw,&
       iLineGlobalPw_I,&
       AltPwIn_C,PrecipEminPwIn,PrecipEmaxPwIn,PrecipEmeanPwIn,PrecipEfluxPwIn,&
       PolarRainEminPwIn,PolarRainEmaxPwIn,PolarRainEmeanPwIn,&
       PolarRainEfluxPwIn, OvationEminPwIn,OvationEmaxPwIn)
    use ModSeGrid, only: iLineGlobal,IsVerbose,AltPwUpper,set_altgrid,set_egrid
    use ModSeBackground,only: allocate_background_arrays,DoAlignDipoleRot,&
                              DoUsePWOM
    use ModElecTrans,only: PrecipEmin, PrecipEmax, &
         PrecipEmean,PrecipEflux,UsePrecipitation,&
         PolarRainEmin, PolarRainEmax, &
         PolarRainEmean,PolarRainEflux,UsePolarRain,&
         OvationEmin, OvationEmax, UseOvation
    use ModPlanetConst, only: Planet_, NamePlanet_I
    logical, intent(in) :: IsVerbosePw
    integer, intent(in) :: nAltPwIn, nLinePw,iLineGlobalPw_I(nLinePw)
    real,    intent(in) :: AltPwIn_C(nAltPwIn)
    real, optional, intent(in)::PrecipEminPwIn,PrecipEmaxPwIn, &
                                PrecipEmeanPwIn,PrecipEfluxPwIn
    real, optional, intent(in)::PolarRainEminPwIn,PolarRainEmaxPwIn, &
                                PolarRainEmeanPwIn,PolarRainEfluxPwIn
    real, optional, intent(in)::OvationEminPwIn,OvationEmaxPwIn
    !---------------------------------------------------------------------------
    DoUsePWOM = .true.    
    
    !save global line info
    nLine=nLinePw
    if (.not.allocated(iLineGlobal_I))allocate(iLineGlobal_I(nLine))
    iLineGlobal_I=iLineGlobalPw_I

    ! Set verbose base on PWOM input
    IsVerbose = IsVerbosePw

    ! Set up the grid parameters, actual grid will be set later by each line
    
    ! Allocate the background arrays
    if(IsVerbose) write(*,*) 'allocating background arrays'
    call allocate_background_arrays
    
    ! align dipole and rotation
    DoAlignDipoleRot = .false.


    ! store PWOM altitude grid in module for use by interpolation routines
    ! for passing info between STET and PWOM
    nAltPw=nAltPwIn
    if (.not.allocated(AltPw_C)) allocate(AltPw_C(nAltPw))
    AltPw_C=AltPwIn_C
    
    !set upper bound of PW in grid module
    AltPwUpper=AltPw_C(nAltPw)    

    ! save the grid spacing in PWOM
    dAltPw = AltPw_C(2)-AltPw_C(1)

    ! Set the incomming precipitation
    if (present(PrecipEminPwIn).and.present(PrecipEmaxPwIn) &
         .and.present(PrecipEmeanPwIn).and.present(PrecipEfluxPwIn)) then
       UsePrecipitation = .true.
       PrecipEmax=PrecipEmaxPwIn
       PrecipEmin=PrecipEminPwIn
       PrecipEmean=PrecipEmeanPwIn
       PrecipEflux=PrecipEfluxPwIn
    elseif(present(PrecipEminPwIn).or.present(PrecipEmaxPwIn) &
         .or.present(PrecipEmeanPwIn).or.present(PrecipEfluxPwIn)) then
       call con_stop('PW_error: SE precip update values incomplete')
    endif

    ! Set the incomming polar rain
    if (present(PolarRainEminPwIn).and.present(PolarRainEmaxPwIn) &
         .and.present(PolarRainEmeanPwIn).and.present(PolarRainEfluxPwIn)) then
       UsePolarRain = .true.
       PolarRainEmax=PolarRainEmaxPwIn
       PolarRainEmin=PolarRainEminPwIn
       PolarRainEmean=PolarRainEmeanPwIn
       PolarRainEflux=PolarRainEfluxPwIn
    elseif(present(PolarRainEminPwIn).or.present(PolarRainEmaxPwIn) &
         .or.present(PolarRainEmeanPwIn).or.present(PolarRainEfluxPwIn)) then
       call con_stop('PW_error: SE PolarRain update values incomplete')
    endif

    ! Set the incomming Ovation parameters (note that these will change over 
    ! time and be updated by get_stet_for_pwom)
    if (present(OvationEminPwIn).and.present(OvationEmaxPwIn))then
       UseOvation= .true.
       OvationEmax=OvationEmaxPwIn
       OvationEmin=OvationEminPwIn
    elseif(present(OvationEminPwIn).or.present(OvationEmaxPwIn))then
       call con_stop('PW_error: SE Ovation values incomplete')
    endif

    !set number of ions based on planet
    select case(NamePlanet_I(Planet_))
    case('EARTH')
       nIonPW=1
    case('JUPITER')
       nIonPW=3
    end select

    call set_egrid
    call set_altgrid

       
  end subroutine init_pwom_se_coupling
  
  !=============================================================================
  ! input the pwom grid, thermal e density, and Efield. run stet for iLine
  subroutine get_se_for_pwom(TimePw,UtPw,iLine,Coord_D,CoordG_D,CoordG2_D,&
       eDensPW_C,eTempPW_C,EfieldPW_C,Ap_I,F107,F107A,IYD,&
       SeDensPW_C, SeFluxPW_C, SeHeatPW_C, IonRatePW_C, PhotoIonRatePW_IC, &
       SecIonRatePW_IC,EMeanDiffPW,EFluxDiffPW,EMeanWavePW,&
       EFluxWavePW,EMeanMonoPW,EFluxMonoPW)
    use ModSeGrid, only: Efield_C,iLineGlobal,IsVerbose,calc_potential
    use ModSeBackground,only: mLat,mLon, gLat,gLon,&
         Idate, UT,fill_thermal_plasma_empirical,&
         plot_background,plot_ephoto_prod,get_neutrals_and_pe_spectrum
    use ModElecTrans,only: Time, EMeanDiff,EFluxDiff,EMeanWave,&
         EFluxWave,EMeanMono,EFluxMono,etrans

    implicit none
    ! Incomming time from PWOM
    real, intent(in) :: TimePw
    ! Universal Time for background
    real, intent(in) :: UtPW
    ! index of line we are working on
    integer, intent(in) :: iLine
    ! magnetic and geo Lat and Lon coord for base of line in degrees
    real,  intent(in) :: Coord_D(2),CoordG_D(2),CoordG2_D(2)
    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(in) :: eDensPW_C(nAltPw), eTempPW_C(nAltPw),EfieldPW_C(nAltPw)
    ! Ap and F107 values from PWOM to set thermosphere in STET
    real,  intent(in) :: Ap_I(7), F107, F107A
    integer,intent(in):: IYD !same as idate, but set in PWOM
    ! SE dens, flux, and heat from SE (interpolated to PWOM grid)
    real,  intent(out)::SeDensPW_C(nAltPw),SeFluxPW_C(nAltPw),SeHeatPW_C(nAltPw)
    ! Ionization rate from STET(interpolated to PWOM grid)
    real,optional,  intent(out):: IonRatePW_C(nAltPw)
    ! Photo Ionization rate from STET for each ion(interpolated to PWOM grid)
    real,optional,  intent(out):: PhotoIonRatePW_IC(nIonPW,nAltPw)
    ! Secondary Ionization rate from SE for each ion(interpolated to PWOM grid)
    real,optional,  intent(out):: SecIonRatePW_IC(nIonPW,nAltPw)
    
    !OVATION precipitation parameters
    real,optional,  intent(in) :: EmeanDiffPW,EFluxDiffPW,EMeanWavePW,&
         EFluxWavePW, EMeanMonoPW,EFluxMonoPW
    
    ! named parameters for coordinates
    integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID
    ! Is line open, for now always assume yes, but this could be passed
    logical :: IsOpen=.true.
    ! Conversion between statV/cm in PWOM to Volts/m needed for STET
    real, parameter :: cSTATVperCMtoVperM=30000.0
    !---------------------------------------------------------------------------
    if(IsVerbose) then
       write(*,*) '!!!!!!!!!!!!!!!!!!!'
       write(*,*) 'CALLING get_se_for_pwom at time=',TimePw
       write(*,*) 'Working on iLine (local,global)=',iLine,iLineGlobal_I(iLine)
       write(*,*) '!!!!!!!!!!!!!!!!!!!'
    endif
    ! Set the Time in SE to match the time in PWOM
    Time=TimePw
    
    !set global line
    iLineGlobal=iLineGlobal_I(iLine)

    ! set idate needed for MSIS call
    Idate=IYD
    
    ! set UT in ModBackground to match UT from PWOM
    UT=UtPw

    ! Set the Lat and Lon coordinates and Lshell for field line
    mLat=Coord_D(Lat_) 
    mLon=Coord_D(Lon_) 
    gLat=CoordG_D(Lat_) 
    gLon=CoordG_D(Lon_) 

    ! Set Efield to zero everywhere. So only PWOM efield is used.
    Efield_C(:)=0.0
    
    !set Ovation precip 
    if (present(EMeanDiffPW).or.present(EFluxDiffPW) .or. &
         present(EMeanWavePW).or.present(EFluxWavePW).or. &
         present(EMeanMonoPW).or.present(EFluxMonoPW))then
       EMeanDiff=EMeanDiffPW
       EFluxDiff=EFluxDiffPW
       EMeanWave=EMeanWavePW
       EFluxWave=EFluxWavePW
       EMeanMono=EMeanMonoPW
       EFluxMono=EFluxMonoPW
    endif

    ! Get the neutral atmosphere and photo e production spectrum
    call get_neutrals_and_pe_spectrum(F107,F107A,AP_I)

    ! Fill the background arrays
    call fill_thermal_plasma_empirical(F107,F107A,Time)
        

    ! Interpolate the PW values onto SE. This overwrites the thermal plasma 
    ! in the overlap region of the grids
    call interpolate_pwom_to_se(eDensPW_C,eTempPW_C,&
         EfieldPW_C*cSTATVperCMtoVperM)

    if(IsVerbose) then
       write(*,*) 'maxval(Efield_C(:))',maxval(Efield_C(:))
       write(*,*) 'maxval(EfieldPW_C(:))',&
            maxval(EfieldPW_C(:)*cSTATVperCMtoVperM)
    endif
    ! plot background after interp
!    call plot_background(iLine,1,time)
!    call con_stop('')     

    ! get the potential for the mapping above the 2-stream solution
    call calc_potential
    
    ! Get a new steady state solution 
    call etrans
    
    ! Interpolate the output back to PWOM grid
    if (present(IonRatePW_C) .and. present(PhotoIonRatePW_IC)) then
       call interpolate_se_to_PWOM(SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
            IonRatePW_C=IonRatePW_C,PhotoIonRatePW_IC=PhotoIonRatePW_IC, &
            SecIonRatePW_IC=SecIonRatePW_IC)
    elseif(present(IonRatePW_C).and. .not.present(PhotoIonRatePW_IC)) then
       call interpolate_se_to_PWOM(SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
            IonRatePW_C=IonRatePW_C)
    elseif(.not.present(IonRatePW_C) .and. present(PhotoIonRatePW_IC)) then
       call interpolate_se_to_PWOM(SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
            PhotoIonRatePW_IC=PhotoIonRatePW_IC,SecIonRatePW_IC=SecIonRatePW_IC)
    else
       call interpolate_se_to_PWOM(SeDensPW_C,SeFluxPW_C,SeHeatPW_C)
    endif
    
    if(IsVerbose) then
       write(*,*) '!!!!!!!!!!!!!!!!!!!'
       write(*,*) 'FINISH CALLING get_se_for_pwom at time=',TimePw
       write(*,*) '!!!!!!!!!!!!!!!!!!!'
    endif
!    call con_stop('')     
  end subroutine get_se_for_pwom
  
  
  !=============================================================================

  subroutine interpolate_se_to_PWOM(SeDensPW_C,SeFluxPW_C,SeHeatPW_C,&
       IonRatePW_C, PhotoIonRatePW_IC, SecIonRatePW_IC)
    use ModSeGrid,      only: AltExtended_C,nAltExtended,Alt_C,nAlt
    use ModElecTrans,   only: NumberDens_C,NumberFlux_C,HeatingRate_C, &
                              TotalIonizationRate_C, SecondaryIonRate_IC
    use ModPlanetConst, only: Planet_, NamePlanet_I
    use ModSeBackground,  only: PhotoIonRate_IC
    use ModInterpolate, only: linear
    implicit none

    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(out)::SeDensPW_C(nAltPw),SeFluxPW_C(nAltPw),SeHeatPW_C(nAltPw)
    real,  optional, intent(out)::IonRatePW_C(nAltPw)
    real,  optional, intent(out)::PhotoIonRatePW_IC(nIonPw,nAltPw)
    real,  optional, intent(out)::SecIonRatePW_IC(nIonPw,nAltPw)
    
    real :: Coord ! coordinate for interpolation
    integer :: iAlt
    real, parameter :: cEVtoErg=1.60217657e-12
    !named parameters for referencing SE ionization arrays
    integer, parameter :: H2plus_=1,Heplus_=2,Hplus_=3,CH4plus_=4,&
         CH3plus_=5,CH2plus_=6,CHplus_=7
    !--------------------------------------------------------------------------
    
    !loop over each point on PWOM grid, SE outputs to PWOM grid
    ALONG_LINE: do iAlt=1,nAltPw
       Coord= AltPw_C(iAlt)
       SeDensPW_C(iAlt) = &
            linear(NumberDens_C(:),1,nAltExtended,Coord,&
            AltExtended_C)
       SeFluxPW_C(iAlt) = &
            linear(NumberFlux_C(:),1,nAltExtended,Coord,&
            AltExtended_C)
       SeHeatPW_C(iAlt) = &
            linear(HeatingRate_C(:),1,nAltExtended,Coord,&
            AltExtended_C)*cEVtoErg
       if (present(IonRatePW_C)) then
          IonRatePW_C(iAlt) = &
               linear(TotalIonizationRate_C(:),1,nAltExtended,Coord,&
               AltExtended_C)
       endif
       if (present(PhotoIonRatePW_IC)) then
          select case(NamePlanet_I(Planet_))
          case('EARTH')
             !not working for Earth yet
             PhotoIonRatePW_IC(OplusPW_,iAlt) = 0
          case('JUPITER')
             If (Coord>Alt_C(nAlt)) then
                ! when pw above iono set photoionization rate to zero
                PhotoIonRatePW_IC(:,iAlt) = 0.0
             else
                PhotoIonRatePW_IC(H2plusPW_,iAlt) = &
                     linear(PhotoIonRate_IC(H2plus_,:),1,nAlt,&
                     Coord,Alt_C)
                PhotoIonRatePW_IC(HplusPW_,iAlt) = &
                     linear(PhotoIonRate_IC(Hplus_, :),1,nAlt, &
                     Coord,Alt_C)
                PhotoIonRatePW_IC(H3plusPW_,iAlt) = 0.0
             end if
          end select
       endif
       if (present(SecIonRatePW_IC)) then
          select case(NamePlanet_I(Planet_))
          case('EARTH')
             !not working for Earth yet
             SecIonRatePW_IC(OplusPW_,iAlt) = 0
          case('JUPITER')
             If (Coord>Alt_C(nAlt)) then
                ! when pw above iono set secondary ionization rate to zero
                SecIonRatePW_IC(:,iAlt) = 0.0
             else
                SecIonRatePW_IC(H2plusPW_,iAlt) = &
                     linear(SecondaryIonRate_IC(H2plus_,:),1,nAlt,&
                     Coord,Alt_C)
                SecIonRatePW_IC(HplusPW_,iAlt) = &
                     linear(SecondaryIonRate_IC(Hplus_, :),1,nAlt, &
                     Coord,Alt_C)
                SecIonRatePW_IC(H3plusPW_,iAlt) = 0.0
             end if
          end select
       endif
       !write(*,*) AltPw_C(iAlt)/1e5,SeDensPW_C(iAlt),SeFluxPW_C(iAlt),SeHeatPW_C(iAlt)
    enddo ALONG_LINE
  end subroutine interpolate_se_to_PWOM
  !=============================================================================
  ! interpolate PWOM thermal electron density and temperature onto SE grid
  
  subroutine interpolate_pwom_to_se(eDensPW_C,eTempPW_C,EfieldPW_C)
    use ModSeGrid,      only:  AltExtended_C,nAltExtended,Alt_C,nAlt,Efield_C
    use ModSeBackground,only: eThermalDensity_C,eThermalTemp_C
    use ModInterpolate, only: linear
    implicit none

    ! thermal e dens [/cc], temp [k] and E|| [V/m] from PWOM
    real,  intent(in) :: eDensPW_C(nAltPw), eTempPW_C(nAltPw),EfieldPW_C(nAltPw)
    
    real :: NormCoord ! normalized coordinate for interpolation
    integer :: iPoint
    real, parameter :: cKelvinToEV=8.621738e-5
    !--------------------------------------------------------------------------
    
    !assume Efield is zero outside of overlap region
    Efield_C(:)=0.0        

    !loop over each point on line and interpolate PWOM values to SE
    ALONG_LINE: do iPoint=1,nAltExtended
       if (AltExtended_C(iPoint) < AltPw_C(nAltPw)) then
          if (AltExtended_C(iPoint) > AltPw_C(1))then
             if (iPoint <=nAlt) then
                NormCoord=&
                     (Alt_C(iPoint)-AltPw_C(1))/dAltPw+1
                
                eThermalDensity_C(iPoint)=&
                     linear(eDensPW_C,1,nAltPw,NormCoord)
                eThermalTemp_C(iPoint)=&
                     linear(eTempPW_C,1,nAltPw,NormCoord)*cKelvinToEV
             else
                ! special treatment for Efield since we do not set field in iono
                Efield_C(iPoint)=&
                     linear(EfieldPW_C,1,nAltPw,NormCoord)
             end if
          endif
       else
          !once we are above the PW grid stop the do loop
          exit ALONG_LINE
       endif
    enddo ALONG_LINE
    
  end subroutine interpolate_pwom_to_se
  
end Module ModCouplePWOMtoSE
