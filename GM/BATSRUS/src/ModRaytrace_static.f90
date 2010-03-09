!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
module ModRaytrace

  use ModSize
  use ModKind
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  save

  ! Logical parameter indicating static vs. dynamic allocation
  logical, parameter :: IsDynamicRaytrace = .false.

  ! Select between fast less accurate and slower but more accurate algorithms
  logical :: UseAccurateTrace    = .false. 
  logical :: UseAccurateIntegral = .true.

  ! Possible tasks
  logical :: DoTraceRay     = .true.  ! trace rays from all cell centers
  logical :: DoMapRay       = .false. ! map rays down to the ionosphere
  logical :: DoExtractRay   = .false. ! extract info along the rays into arrays
  logical :: DoIntegrateRay = .false. ! integrate some functions along the rays

  ! Use old IJK based logic for Cartesian tracing
  logical :: UseOldMethodOfRayTrace = .true.

  ! Number of rays per dimension on the starting grid
  ! This is needed for DoExtractRay = .true. only
  integer :: nRay_D(4) = (/0, 0, 0, 0/)

  ! The vector field to trace: B/U/J
  character         :: NameVectorField = 'B'

  ! How often shall we synchronize PE-s for the accurate algorithms
  real         :: DtExchangeRay = 0.1

  ! The minimum number of time steps between two ray traces on the same grid
  integer      :: DnRaytrace = 1

  ! Named parameters for ray status (must be less than east_=1)
  integer, parameter :: &
       ray_iono_ = 0, &
       ray_block_=-1, &
       ray_open_ =-2, &
       ray_loop_ =-3, &
       ray_body_ =-4, &
       ray_out_  =-5

  ! Ray and rayface contain the x,y,z coordinates for the foot point of a given
  ! field line for both directions, eg. 
  ! ray(2,1,i,j,k,iBLK) is the y coord for direction 1

  ! ray is for cell centers; rayface is for block surfaces with 
  ! a -0.5,-0.5,-0.5 shift in block normalized coordinates

  real, dimension(3,2,nI+1,nJ+1,nK+1,nBLK) :: ray, rayface

  ! Stored face and cell indices of the 2 rays starting from a face of a block
  integer :: rayend_ind(3,2,nI+1,nJ+1,nK+1,nBLK)

  ! Stored weights for the 2 rays starting from a face of a block
  real    :: rayend_pos(4,2,nI+1,nJ+1,nK+1,nBLK)

  ! Radius where ray tracing with numerical B stops and 
  ! radius and radius squared of ionosphere

  real :: R_raytrace=1., R2_raytrace=1.

!!! These could be allocatable arrays ???
  
  ! Node interpolated magnetic field components without B0
  real, dimension(1:nI+1,1:nJ+1,1:nK+1,nBLK):: bb_x,bb_y,bb_z

  ! Total magnetic field with second order ghost cells
  real, dimension(3,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: Bxyz_DGB

  ! Prefer open and closed field lines in interpolation ?!
  logical :: UsePreferredInterpolation

  ! Maximum length of ray
  real :: RayLengthMax

  ! Testing
  logical :: oktest_ray=.false.

  ! Constants to distinguish various ray types
  real, parameter :: rIonosphere = 1.0, rIonosphere2 = rIonosphere**2
  real, parameter :: &
       CLOSEDRAY= -(rIonosphere + 0.05), &
       OPENRAY  = -(rIonosphere + 0.1), &
       BODYRAY  = -(rIonosphere + 0.2), &
       LOOPRAY  = -(rIonosphere + 0.3), &
       NORAY    = -(rIonosphere + 100.0), &
       OUTRAY   = -(rIonosphere + 200.0)

  ! Base time for timed exchanges between rays
  real(Real8_) :: CpuTimeStartRay

  ! Number of rays found to be open based on the neighbors
  integer      :: nOpen

  ! Logical for raytracing in IE coupling
  logical :: DoTraceIE = .false.

  ! ----------- Variables for mapping rays to the ionosphere ------------
  real, allocatable:: RayMapLocal_DSII(:,:,:,:), RayMap_DSII(:,:,:,:)

  ! ----------- Variables for extracting variables along the ray --------
  logical :: DoExtractState = .false., DoExtractUnitSi = .false.

  ! ----------- Variables for integrals along the ray -------------------
  ! Named indexes
  integer, parameter :: &
       InvB_=1, Z0x_=2, Z0y_=3, Z0b_=4, RhoInvB_=5, pInvB_=6, &
       HpRhoInvB_ = 7, OpRhoInvB_ = 8, HppInvB_ =9, OppInvB_=10, &
       xEnd_=11, yEnd_=12, zEnd_=13, Length_=14

  ! Number of integrals
  integer, parameter :: nRayIntegral = 14

  ! Flow variables to be integrated (rho and P) other than the magnetic field
  real, dimension(2,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: Extra_VGB
  real, dimension(4,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: ExtraMulti_VGB

  ! Integrals for a local ray segment
  real :: RayIntegral_V(InvB_:OppInvB_)

  ! Integrals added up for all the local ray segments
  ! The fist index corresponds to the variables (index 0 shows closed vs. open)
  ! The second and third indexes correspond to the latitude and longitude of
  ! the IM/RCM grid
  real, allocatable :: RayIntegral_VII(:,:,:)
  real, allocatable :: RayResult_VII(:,:,:)

  ! Temporary array for extracting b.grad(B1) info
  real, allocatable :: bGradB1_DGB(:,:,:,:,:)
  logical:: DoExtractBGradB1 = .false.

  ! Conversion matrix between SM and GM coordinates 
  ! (to be safe initialized to unit matrix)
  real :: GmSm_DD(3,3) = reshape( (/ &
       1.,0.,0., &
       0.,1.,0., &
       0.,0.,1. /), (/3,3/) )

  integer :: iLatTest = 1, iLonTest = 1

contains

  !============================================================================

  subroutine xyz_to_latlon(Pos_D)

    use ModNumConst, ONLY: cTiny, cRadToDeg

    ! Convert xyz coordinates to latitude and longitude (in degrees)
    ! Put the latitude and longitude into the 1st and 2nd elements
    real, intent(inout) :: Pos_D(3)

    real :: x, y, z

    !-------------------------------------------------------------------------

    ! Check if this direction is closed or not
    if(Pos_D(1) > CLOSEDRAY)then

       ! Convert GM position into IM position
       Pos_D = matmul(Pos_D, GmSm_DD)

       ! Store input coordinates
       x = Pos_D(1); y = Pos_D(2); z = Pos_D(3)

       ! Make sure that asin will work, -1<= z <=1
       z = max(-1.0+cTiny, z)
       z = min( 1.0-cTiny, z)

       ! Calculate  -90 < latitude = asin(z)  <  90
       Pos_D(1) = cRadToDeg * asin(z)

       ! Calculate -180 < longitude = atan2(y,x) < 180
       if(abs(x) < cTiny .and. abs(y) < cTiny) x = 1.0
       Pos_D(2) = cRadToDeg * atan2(y,x)

       ! Get rid of negative longitude angles
       if(Pos_D(2) < 0.0) Pos_D(2) = Pos_D(2) + 360.0

    else
       ! Impossible values
       Pos_D(1) = -100.
       Pos_D(2) = -200.
    endif

  end subroutine xyz_to_latlon

  !============================================================================

  subroutine xyz_to_latlonstatus(Ray_DI)

    real, intent(inout) :: Ray_DI(3,2)

    integer :: iRay
    !-------------------------------------------------------------------------

    ! Convert 1st and 2nd elements into latitude and longitude
    do iRay=1,2
       call xyz_to_latlon(Ray_DI(:,iRay))
    end do

    ! Convert 3rd element into a status variable

    if(Ray_DI(3,1)>CLOSEDRAY .and. Ray_DI(3,2)>CLOSEDRAY)then
       Ray_DI(3,:)=3.      ! Fully closed
    elseif(Ray_DI(3,1)>CLOSEDRAY .and. Ray_DI(3,2)==OPENRAY)then
       Ray_DI(3,:)=2.      ! Half closed in positive direction
    elseif(Ray_DI(3,2)>CLOSEDRAY .and. Ray_DI(3,1)==OPENRAY)then
       Ray_DI(3,:)=1.      ! Half closed in negative direction
    elseif(Ray_DI(3,1)==OPENRAY .and. Ray_DI(3,2)==OPENRAY) then
       Ray_DI(3,:)=0.      ! Fully open
    elseif(Ray_DI(3,1)==BODYRAY)then
       Ray_DI(3,:)=-1.     ! Cells inside body
    elseif(Ray_DI(3,1)==LOOPRAY .and.  Ray_DI(3,2)==LOOPRAY) then
       Ray_DI(3,:)=-2.     ! Loop ray within block
    else
       Ray_DI(3,:)=-3.     ! Strange status
    end if

  end subroutine xyz_to_latlonstatus

  !============================================================================

  subroutine init_mod_raytrace

    ! True if ray array is still to be initialized
    logical :: DoInitRay = .true.

    ! Initialize ray array (write_logfile may use it before first ray tracing)
    if(DoInitRay)then
       ray       = 0.0
       DoInitRay = .false.
    end if

    if(IsDynamicRaytrace .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_raytrace allocated arrays'
    end if

  end subroutine init_mod_raytrace

  !============================================================================

  subroutine clean_mod_raytrace

    if(IsDynamicRaytrace .and. iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_raytrace deallocated arrays'
    end if

  end subroutine clean_mod_raytrace

  !============================================================================

end module ModRaytrace
