! Initializes the stet code, 

subroutine stet_init(nLineIn,Coord_ID,Ap_I,F107,F107A, TimeIn)
  use ModSeGrid, only:update_grid,nLine,Efield_IC,DoIncludePotential,IsVerbose
  use ModSeBackground,only: allocate_background_arrays,mLat_I,mLon_I, &
       set_footpoint_locations,fill_thermal_plasma_empirical,plot_background,&
       plot_ephoto_prod,DoAlignDipoleRot,get_neutrals_and_pe_spectrum
  use ModNumConst,    ONLY: cDegToRad
  use ModSeState,     ONLY: Time,allocate_state_arrays
  implicit none 
  
  integer, intent(in):: nLineIn ! number of lines to be worked on
  real, intent(in) :: Coord_ID(nLineIn,2) ! Lat and Lon in degrees for each line
  
  !Thermospheric inputs
  real, intent(in) :: Ap_I(7),F107,F107A
  
  !Time input
  real, intent(in) :: TimeIn
  
  integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_ID
  
  integer :: iLine !loop variable for line
  !-----------------------------------------------------------------------------
  
  ! Set nLine to nLineIn
  nLine=nLineIn
  
  ! Set the Time to TimeIn
  Time=TimeIn
  
  !\
  ! Set up the grid
  !/
  if(IsVerbose) write(*,*) 'setting grid dimensions'
  call set_grid_dimensions_default
  
  do iLine=1,nLine
     call update_grid(iLine,Coord_ID(iLine,:),DoOnlySpatial=.true.)
  enddo

  !define the electric field 
  Efield_IC(:,:)=0.0

  ! after electric field definition we can finish setting up grid
  if (DoIncludePotential) then
     do iLine=1,nLine
        call update_grid(iLine,Coord_ID(iLine,:),DoOnlySpatial=.true.)
     enddo
  endif

  !\
  ! Set the background arrays, sources, and locations
  !/
  ! Allocate the background right
  if(IsVerbose) write(*,*) 'allocating background arrays'
  call allocate_background_arrays
  
  !align dipole and rotation
  DoAlignDipoleRot = .true.
  
  do iLine=1,nLine
     ! set location of field line from inputs
     mLat_I(iLine)=Coord_ID(iLine,Lat_) 
     mLon_I(iLine)=Coord_ID(iLine,Lon_) 
     
     !set glat and glon coords
     call set_footpoint_locations(iLine)
     if(IsVerbose) write(*,*) 'finished setting footpoints'
     ! Get the neutral atmosphere and photo e production spectrum
     call get_neutrals_and_pe_spectrum(iLine,F107,F107A,AP_I)
     if(IsVerbose) write(*,*) 'finished getting neutrals and pe spec'
     ! Fill the background arrays
     if(IsVerbose) write(*,*) 'filling background arrays'
     call fill_thermal_plasma_empirical(iLine,F107,F107A,Time)
  end do
  
  call allocate_state_arrays
  

end subroutine stet_init

!===============================================================================
! set the default grid dimensions and parameters. 
!This should only be done once at the  start of the run
subroutine set_grid_dimensions_default
  use ModSeGrid, only: allocate_grid_arrays, nIono,&
       DrIono1,nIono1,DrIono2,nIono2,DrIono3,nIono3,DrIono4,nIono4,TypeGridE,&
       nEnergy,EnergyMax,DeltaE, nTheta_II
  use ModPlanetConst, ONLY: Planet_, NamePlanet_I
  !-----------------------------------------------------------------------------
 ! For now this is the same as the unit test grid.
  write(*,*) 'creating grid'
  !  call create_se_test_grid
  
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
  !nEnergy=100
  !nEnergy=99
  !    nEnergy=94

  ! Energy grid for only photoelectrons

  nEnergy=299
  EnergyMax=300.5


  ! Energy grid for precipitation
!  nEnergy=999  
!  EnergyMax=1000.5
  DeltaE = 1.0
  
  ! Allocated the grid arrays and populate the bfield, sgrid, and PA grid  
  write(*,*) 'allocating arrays'
  call allocate_grid_arrays
  
  ! default values for theta grid. 
  nTheta_II(:,1)=5
  nTheta_II(:,2)=20
  nTheta_II(:,3)=90
  nTheta_II(:,4)=20
  nAngle = 135

end subroutine set_grid_dimensions_default

