!^CFG COPYRIGHT UM
!^CFG FILE RAYTRACE
module ModRaytrace

  use ModSize
  use ModKind
  implicit none
  save

  ! Select between fast less accurate and slower but more accurate algorithms
  logical :: UseAccurateTrace = .false. 

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

  ! Testing
  logical :: oktest_ray=.false.

  real, parameter :: rIonosphere = 1.0, rIonosphere2 = rIonosphere**2
  real, parameter :: &
       CLOSEDRAY= -(rIonosphere + 0.05), &
       OPENRAY  = -(rIonosphere + 0.1), &
       BODYRAY  = -(rIonosphere + 0.2), &
       LOOPRAY  = -(rIonosphere + 0.3), &
       NORAY    = -(rIonosphere + 100.0), &
       OUTRAY   = -(rIonosphere + 200.0)

  real(Real8_) :: CpuTimeStartRay

  real         :: DtExchangeRay = 0.1

  integer      :: nOpen

  ! ----------- Variables for integrals along the ray -------------------
  ! True if the ray integrals are done
  logical :: DoIntegrate = .false.  

  ! Name indexes
  integer, parameter :: &
       ClosedRay_ = 0, InvB_=1, Z0x_=2, Z0y_=3, Z0b_=4, RhoInvB_=5, pInvB_=6 

  ! Number of integrals
  integer, parameter :: nRayIntegral = 6

  ! Flow variables to be integrated (rho and P) other than the magnetic field
  real, dimension(2,-1:nI+2,-1:nJ+2,-1:nK+2,nBLK) :: Extra_VGB

  ! Integrals for a local ray segment
  real :: RayIntegral_V(nRayIntegral)

  ! Integrals added up for all the local ray segments
  ! The fist index corresponds to the variables (index 0 shows closed vs. open)
  ! The second and third indexes correspond to the latitude and longitude of
  ! the IM/RCM grid
  real, allocatable :: RayIntegral_VII(:,:,:)

end module ModRaytrace
