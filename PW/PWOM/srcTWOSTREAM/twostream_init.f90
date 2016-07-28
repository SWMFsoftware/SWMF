! Initializes the stet code, 

subroutine twostream_init(Coord_D,Ap_I,F107,F107A, TimeIn)
  use ModSeGrid, only:set_egrid,set_altgrid,IsVerbose,Time
  use ModSeBackground,only: allocate_background_arrays,mLat,mLon, &
       set_footpoint_locations,fill_thermal_plasma_empirical,plot_background,&
       plot_ephoto_prod,DoAlignDipoleRot,get_neutrals_and_pe_spectrum
  use ModNumConst,    ONLY: cDegToRad

  implicit none 
  
  real, intent(in) :: Coord_D(2) ! Lat and Lon in degrees for each line
  
  !Thermospheric inputs
  real, intent(in) :: Ap_I(7),F107,F107A
  
  !Time input
  real, intent(in) :: TimeIn
  
  integer,parameter :: Lat_=1 ,Lon_=2 !named parameters for Coord_D
  
  !-----------------------------------------------------------------------------
  
  ! Set the Time to TimeIn
  Time=TimeIn
  
  !\
  ! Set up the grid
  !/
  if(IsVerbose) write(*,*) 'setting grids'
  call set_egrid
  call set_altgrid
  

!  !define the electric field 
!  Efield_IC(:,:)=0.0
!

  !\
  ! Set the background arrays, sources, and locations
  !/
  ! Allocate the background right
  if(IsVerbose) write(*,*) 'allocating background arrays'
  call allocate_background_arrays
  
  !align dipole and rotation
  DoAlignDipoleRot = .true.
  
  ! set location of field line from inputs
  mLat=Coord_D(Lat_) 
  mLon=Coord_D(Lon_) 
  
  !set glat and glon coords
  call set_footpoint_locations
  if(IsVerbose) write(*,*) 'finished setting footpoints'
  ! Get the neutral atmosphere and photo e production spectrum
  call get_neutrals_and_pe_spectrum(F107,F107A,AP_I)
  if(IsVerbose) write(*,*) 'finished getting neutrals and pe spec'
  ! Fill the background arrays
  if(IsVerbose) write(*,*) 'filling background arrays'
  call fill_thermal_plasma_empirical(F107,F107A,Time)
  
  
  !  call allocate_state_arrays
  
  
end subroutine twostream_init

