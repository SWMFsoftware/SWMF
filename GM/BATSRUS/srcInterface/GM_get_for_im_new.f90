!===========================================================================

subroutine GM_get_for_im_new

  use CON_coupler, ONLY: Grid_C, IM_
  use ModNumConst, ONLY: cRadToDeg

  implicit none
  integer :: nLat, nLon
  !---------------------------------------------------------------------------

  ! Transformation matrix from IM to GM
  !GmIm_DD = transform_matrix(tSimulation,&
  !     Grid_C(IM_) % TypeCoord, Grid_C(GM_) % TypeCoord)

  ! Grid size of IM grid
  nLat = Grid_C(IM_) % nCoord_D(1)
  nLon = Grid_C(IM_) % nCoord_D(2)

  call integrate_ray_accurate(nLat,nLon,&
       90.0 - cRadToDeg*Grid_C(IM_) % Coord1_I,&
       cRadToDeg*Grid_C(IM_) % Coord2_I)

end subroutine GM_get_for_im_new

