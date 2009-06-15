
Module ModTides

  use ModSizeGITM, only: nLons, nLats, nBlocksMax

  implicit none

  real, allocatable :: u_gswm(:,:,:,:)
  real, allocatable :: v_gswm(:,:,:,:)
  real, allocatable :: T_gswm(:,:,:,:)
  real, allocatable :: lon_gswm(:)
  real, allocatable :: lat_gswm(:)

  integer :: nLatsGSWM, nLonsGSWM, nAltsGswm

  character(len=40), dimension(4) :: GSWM_name           ! EY: 18May09
  character(len=40), dimension(4) :: GSWM_file_name      ! EY: 18May09

  real, dimension(-1:nLons+2, -1:nLats+2, 2, nBlocksMax) :: &
       TidesNorth, TidesEast, TidesTemp

  real :: dLonGswm, dLatGswm

end Module ModTides


