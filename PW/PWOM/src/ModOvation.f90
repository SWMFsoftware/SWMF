module ModOvation
implicit none

private

logical,public:: UseOvation = .false.
logical,public:: DoPlotOvation = .false.
real,public :: StartTime
real,public :: OvationEmin,OvationEmax
! Time Cadence at which ovation data is provided
integer,parameter :: DtReadOvation = 300.0

! parameters for ovation grid
integer,parameter :: nMlt=96,nLat=80
!integer,parameter :: nMlt=80,nLat=96
real ,parameter   :: LatMax=89.5,LatMin=50
integer :: Mlt_=1,Lat_=2

! hold ovation data
real :: efluxDiff_G(0:nMlt+1,0:nLat+1)
real :: nfluxDiff_G(0:nMlt+1,0:nLat+1)
real :: CoordDiff_DG(2,0:nMlt+1,0:nLat+1)

real :: efluxMono_G(0:nMlt+1,0:nLat+1)
real :: nfluxMono_G(0:nMlt+1,0:nLat+1)
real :: CoordMono_DG(2,0:nMlt+1,0:nLat+1)

real :: efluxWave_G(0:nMlt+1,0:nLat+1)
real :: nfluxWave_G(0:nMlt+1,0:nLat+1)
real :: CoordWave_DG(2,0:nMlt+1,0:nLat+1)
! current ovation data file
Character(len=200) :: NameFile

integer :: iTimeRead_I(7)

character(len=3)::NameMonth

! public methods
public :: unit_test_ovation
public :: get_ovation_point
public :: read_ovation_all
public :: plot_ovation_polar    

contains
  !=============================================================================
  ! get Eflux[ergs/cm2/s] and Emean[eV] for input SmLat and SmLon for 
  ! each type of precip
  subroutine get_ovation_point(SmLat,SmLon,EMeanDiff,EFluxDiff,&
       EMeanWave,EFluxWave,EMeanMono,EFluxMono)
    use ModInterpolate, ONLY: bilinear
    real, intent(in) :: SmLat,SmLon !input in degrees
    real, intent(out):: EMeanDiff,EFluxDiff,&
         EMeanWave,EFluxWave,EMeanMono,EFluxMono
    real,parameter :: cLonToMlt=0.0666667 !lon (deg) to mlt
    real,parameter :: cErgToeV=6.242e11
    real :: nFluxDiff, nFluxMono, nFluxWave,Mlt
    !---------------------------------------------------------------------------
    
    if (SmLat>CoordDiff_DG(Lat_,1,1)) then
       Mlt = modulo(SmLon*cLonToMlt+12.0,24.0)
       EfluxDiff= &
            bilinear(efluxDiff_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordDiff_DG(Mlt_,:,1),CoordDiff_DG(Lat_,1,:),DoExtrapolate=.true.)
       nfluxDiff= &
            bilinear(nfluxDiff_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordDiff_DG(Mlt_,:,1),CoordDiff_DG(Lat_,1,:),DoExtrapolate=.true.)
       EMeanDiff=EfluxDiff/nFluxDiff*cErgToeV

       EfluxMono= &
            bilinear(efluxMono_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordMono_DG(Mlt_,:,1),CoordMono_DG(Lat_,1,:),DoExtrapolate=.true.)
       nfluxMono= &
            bilinear(nfluxMono_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordMono_DG(Mlt_,:,1),CoordMono_DG(Lat_,1,:),DoExtrapolate=.true.)
       EMeanMono=EfluxMono/nFluxMono*cErgToeV
       
       EfluxWave= &
            bilinear(efluxWave_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordWave_DG(Mlt_,:,1),CoordWave_DG(Lat_,1,:),DoExtrapolate=.true.)
       nfluxWave= &
            bilinear(nfluxWave_G,0,nMlt+1,0,nLat+1,(/Mlt,SmLat/),&
            CoordWave_DG(Mlt_,:,1),CoordWave_DG(Lat_,1,:),DoExtrapolate=.true.)
       EMeanWave=EfluxWave/nFluxWave*cErgToeV
    else
       ! do not interpolate below min latitude
       EfluxDiff=0.0
       EMeanDiff=0.0
       EfluxMono=0.0
       EMeanMono=0.0
       EfluxWave=0.0
       EMeanWave=0.0
    endif
    
    
  end subroutine get_ovation_point
  !=============================================================================
  subroutine read_ovation_all(tSimulation)
    real, intent(in) :: tSimulation
    
    call read_ovation(tSimulation,'diff')
    call read_ovation(tSimulation,'wave')
    call read_ovation(tSimulation,'mono')
  end subroutine read_ovation_all

  !=============================================================================
  subroutine read_ovation(tSimulation,TypePrecip)
    use ModTimeConvert, ONLY: time_real_to_int
    use ModIoUnit, ONLY: io_unit_new,UnitTmp_
    !use ModPwTime,ONLY: StartTime,
    real, intent(in) :: tSimulation
    character(len=4) :: TypePrecip
    real :: CurrentTime,TimeRead
    real :: eflux_G(0:nMlt+1,0:nLat+1)
    real :: nflux_G(0:nMlt+1,0:nLat+1)
    real :: Coord_DG(2,0:nMlt+1,0:nLat+1)
    integer, parameter :: Year_=1,Month_=2,Day_=3,Hour_=4,Minute_=5,Second_=6
    
    !file names
    Character(len=100) :: NameFileEflux,NameFileNumberFlux
    Character(len=200) :: header
    
    integer :: iMlt, iLat
    
    !---------------------------------------------------------------------------
    CurrentTime=StartTime+tSimulation
    ! get read time closest to simulation time without going over
    TimeRead = (floor(CurrentTime/DtReadOvation) * DtReadOvation)
    call time_real_to_int(TimeRead,iTimeRead_I)
    
    ! make sure seconds are cleaned up
    iTimeRead_I(Second_) = 0
     
    ! Get 3 character month string that corresponds to numerical month
    if (iTimeRead_I(Month_)==1) then
       NameMonth = 'Jan'
    elseif (iTimeRead_I(Month_)==2) then
       NameMonth = 'Feb'
    elseif (iTimeRead_I(Month_)==3) then
       NameMonth = 'Mar'
    elseif (iTimeRead_I(Month_)==4) then
       NameMonth = 'Apr'    
    elseif (iTimeRead_I(Month_)==5) then
       NameMonth = 'May'
    elseif (iTimeRead_I(Month_)==6) then
       NameMonth = 'Jun'
    elseif (iTimeRead_I(Month_)==7) then
       NameMonth = 'Jul'
    elseif (iTimeRead_I(Month_)==8) then
       NameMonth = 'Aug'
    elseif (iTimeRead_I(Month_)==9) then
       NameMonth = 'Sep'
    elseif (iTimeRead_I(Month_)==10) then
       NameMonth = 'Oct'
    elseif (iTimeRead_I(Month_)==11) then
       NameMonth = 'Nov'
    elseif (iTimeRead_I(Month_)==12) then
       NameMonth = 'Dec'
    endif
    
    !sample name: epoch_2013_energy_flux_for_1998-Sep-25_2355.txt
    !reconstruct the file name
    write(NameFile,"(a,a,a,i4.4,a,a,a,i2.2,a,i2.2,i2.2,a,a,a)") &
         'PW/OVATION/',TypePrecip,'/energyflux/epoch_2013_energy_flux_for_',&
         iTimeRead_I(Year_),'-',NameMonth,'-',iTimeRead_I(Day_),'_',&
         iTimeRead_I(Hour_),iTimeRead_I(Minute_),'_',TypePrecip,'.txt'   
    write(*,*) NameFile

    !open file for reading
    open(UnitTmp_,file=NameFile,status="old")
    !discard header
    read(UnitTmp_,*) header

    do iMlt=1,nMlt
       do iLat=1,nLat
          read(UnitTmp_,*) Coord_DG(Mlt_,iMLT,iLat),Coord_DG(Lat_,iMLT,iLat),&
               eflux_G(iMlt,iLat)
       enddo
    enddo
    close(UnitTmp_)

    ! Now get filename for number flux and read
    write(NameFile,"(a,a,a,i4.4,a,a,a,i2.2,a,i2.2,i2.2,a,a,a)") &
         'PW/OVATION/',TypePrecip,'/numberflux/epoch_2013_number_flux_for_',&
         iTimeRead_I(Year_),'-',NameMonth,'-',iTimeRead_I(Day_),'_',&
         iTimeRead_I(Hour_),iTimeRead_I(Minute_),'_',TypePrecip,'.txt'   

    !open file for reading
    open(UnitTmp_,file=NameFile,status="old")
    !discard header
    read(UnitTmp_,*) header

    do iMlt=1,nMlt
       do iLat=1,nLat
          read(UnitTmp_,*) Coord_DG(Mlt_,iMLT,iLat),Coord_DG(Lat_,iMLT,iLat),&
               nflux_G(iMlt,iLat)
       enddo
    enddo
    close(UnitTmp_)
    

    ! fill ghost cells
    Coord_DG(Mlt_,nMlt+1,:)=24.0
    Coord_DG(Mlt_,0,:)=-0.25
    Coord_DG(Lat_,nMlt+1,:)=Coord_DG(Lat_,1,:)
    Coord_DG(Lat_,0,:)=Coord_DG(Lat_,nMlt,:)

    Coord_DG(Lat_,:,nLat+1)=90.0
    Coord_DG(Lat_,:,0)=49.5
    Coord_DG(Mlt_,:,0)=Coord_DG(Mlt_,:,1)
    Coord_DG(Mlt_,:,nLat+1)=Coord_DG(Mlt_,:,nLat)

    eflux_G(nMlt+1,:)= eflux_G(1,:)
    eflux_G(0,:)= eflux_G(nMlt,:)

    eflux_G(:,0)= eflux_G(:,1)
    eflux_G(:,nLat+1)= eflux_G(:,nLat)
    
    nflux_G(nMlt+1,:)= nflux_G(1,:)
    nflux_G(0,:)= nflux_G(nMlt,:)

    nflux_G(:,0)= nflux_G(:,1)
    nflux_G(:,nLat+1)= nflux_G(:,nLat)
      
    ! put into appropriate array
    select case(TypePrecip)
    case('diff')
       CoordDiff_DG=Coord_DG
       efluxDiff_G= eflux_G
       nfluxDiff_G= nflux_G

    case('mono')
       CoordMono_DG=Coord_DG
       efluxMono_G= eflux_G
       nfluxMono_G= nflux_G

    case('wave')
       CoordWave_DG=Coord_DG
       efluxWave_G= eflux_G
       nfluxWave_G= nflux_G
       
    end select

  end subroutine read_ovation
  !===========================================================================
  subroutine plot_ovation
    use ModIoUnit,    ONLY: UnitTmp_
    use ModNumConst, ONLY: cDegToRad, cRadToDeg
    integer :: iLat, iMlt
    Character(len=200) :: NameFileOut
    integer, parameter :: Year_=1,Month_=2,Day_=3,Hour_=4,Minute_=5,Second_=6
    !filename for plotting output
    write(NameFileOut,"(a,i4.4,a,a,a,i2.2,a,i2.2,i2.2,a)") &
         'PW/OVATION/plots/precip_for_',&
         iTimeRead_I(Year_),'-',NameMonth,'-',iTimeRead_I(Day_),'_',&
         iTimeRead_I(Hour_),iTimeRead_I(Minute_),'.dat'   

    
    open(UnitTmp_,FILE=NameFileOut)
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Lat", "MLT", "efluxDiff" "nfluxDiff", "efluxMono", "nfluxMono", "efluxWave" "nfluxWave"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLat, &
         ', J=', nMLT+1,', DATAPACKING=POINT'
    
    do iMlt = 1, nMlt+1
       do iLat = 1, nLat
          write(UnitTmp_,"(100es18.10)") CoordDiff_DG(Lat_,iMLT,iLat),&
               CoordDiff_DG(Mlt_,iMLT,iLat),efluxDiff_G(iMlt,iLat),&
               nfluxDiff_G(iMlt,iLat),efluxMono_G(iMlt,iLat),&
               nfluxMono_G(iMlt,iLat),efluxWave_G(iMlt,iLat),&
               nfluxWave_G(iMlt,iLat)
       enddo
    enddo
    close(UnitTmp_)
  end subroutine plot_ovation
  !===========================================================================
  subroutine plot_ovation_polar
    use ModIoUnit,    ONLY: UnitTmp_
    use ModNumConst, ONLY: cDegToRad, cRadToDeg,cPi
    integer :: iLat, iMlt
    real :: theta,phi,xpos,ypos
    Character(len=200) :: NameFileOut
    integer, parameter :: Year_=1,Month_=2,Day_=3,Hour_=4,Minute_=5,Second_=6
    !---------------------------------------------------------------------------

    !filename for plotting output
    write(NameFileOut,"(a,i4.4,a,a,a,i2.2,a,i2.2,i2.2,a)") &
         'PW/OVATION/plots/polar_precip_for_',&
         iTimeRead_I(Year_),'-',NameMonth,'-',iTimeRead_I(Day_),'_',&
         iTimeRead_I(Hour_),iTimeRead_I(Minute_),'.dat'   

    open(UnitTmp_,FILE=NameFileOut)
    write(UnitTmp_,'(a)') &
         'VARIABLES = "Y", "X", "Lat", "Mlt", "efluxDiff" "nfluxDiff", "efluxMono", "nfluxMono", "efluxWave" "nfluxWave"'
    write(UnitTmp_,'(a,i3,a,i3,a)') 'Zone I=', nLat, &
         ', J=', nMLT+1,', DATAPACKING=POINT'
    
    do iMlt = 1, nMlt+1
       do iLat = 1, nLat
          phi=modulo((CoordDiff_DG(Mlt_,iMLT,iLat)-12.0),24.0)*cPi/12
          theta = (90.0-CoordDiff_DG(Lat_,iMLT,iLat))*cDegToRad
          xpos=sin(theta)*cos(phi)
          ypos=sin(theta)*sin(phi)
          write(UnitTmp_,"(100es18.10)") ypos,xpos,&
               CoordDiff_DG(Lat_,iMLT,iLat),CoordDiff_DG(Mlt_,iMLT,iLat),&
               efluxDiff_G(iMlt,iLat),&
               nfluxDiff_G(iMlt,iLat),efluxMono_G(iMlt,iLat),&
               nfluxMono_G(iMlt,iLat),efluxWave_G(iMlt,iLat),&
               nfluxWave_G(iMlt,iLat)
       enddo
    enddo
    close(UnitTmp_)
  end subroutine plot_ovation_polar

  !============================================================================
  subroutine unit_test_ovation
    use ModTimeConvert, ONLY: time_int_to_real
    integer :: iStartTime_I(7)=(/1998,9,25,23,55,0,0/)
    real :: time=0.0
!    real :: SmLat=77.0,SmLon=0.0
    real :: SmLat=66.0,SmLon=180.0
    real ::EMeanDiff,EFluxDiff,EMeanWave,EFluxWave,EMeanMono,EFluxMono
    !---------------------------------------------------------------------------
    write(*,*) 'testing for time'
    write(*,*) 'year:',iStartTime_I(1)
    write(*,*) 'month:',iStartTime_I(2)
    write(*,*) 'day:',iStartTime_I(3)
    write(*,*) 'hour:',iStartTime_I(4)
    write(*,*) 'min:',iStartTime_I(5)
    write(*,*) 'sec:',iStartTime_I(6)
    write(*,*) 'simulation time:', Time

    call time_int_to_real(iStartTime_I,StartTime)

    call read_ovation_all(time)
    
    call plot_ovation_polar
    call plot_ovation

    call get_ovation_point(SmLat,SmLon,EMeanDiff,EFluxDiff,&
         EMeanWave,EFluxWave,EMeanMono,EFluxMono)

    write(*,*) 'Eflux[ergs/cm2/s] and Emean[eV] for Lat, Lon:',SmLat,SmLon
    write(*,*) 'Diff:', EFluxDiff,EMeanDiff
    write(*,*) 'Mono:', EFluxMono,EMeanMono
    write(*,*) 'Wave:', EFluxWave,EMeanWave
    
  end subroutine unit_test_ovation


end module ModOvation
