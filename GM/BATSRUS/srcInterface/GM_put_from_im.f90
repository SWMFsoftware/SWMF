!^CFG COPYRIGHT UM
!^CMP FILE IM
!==========================================================================
subroutine GM_put_from_im(Buffer_II,iSizeIn,jSizeIn,NameVar)

  !call stop_mpi('RCM is OFF') !^CFG UNCOMMENT IF NOT RCM
  !^CFG IF RCM BEGIN
  use CON_coupler
  use ModImPressure                              ! Storage for IM pressure
  use ModNumConst
  implicit none
  character(len=*), parameter :: NameSub='GM_put_from_im'

  integer, intent(in) :: iSizeIn,jSizeIn
  real, intent(in) :: Buffer_II(iSizeIn,jSizeIn)
  character(len=*), intent(in) :: NameVar
  integer :: nCells_D(2)
  !---------------------------------------------------------------------------
  if(NameVar /= 'p') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  nCells_D=ncells_decomposition_d(IM_)
  if( iSizeIn /= nCells_D(1) .or. jSizeIn /= nCells_D(2) ) then

     write(*,*)NameSub//' grid sizes do not agree iSize,jSize,nCells=',&
          iSizeIn,jSizeIn,nCells_D(1:2)
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  if(.not.allocated(RCM_lat))then
     ! Allocate RCM_lat, RCM_lon, RCM_p
     call im_pressure_init(iSizeIn, jSizeIn)
     ! Convert colat, lon to lat-lon in degrees and store
     RCM_lat = (cHalfPi - Grid_C(IM_) % Coord1_I) * cRadToDeg
     RCM_lon = Grid_C(IM_)% Coord2_I              * cRadToDeg
  end if

  ! Store IM variable for internal use
  RCM_p = Buffer_II

  ! Make sure that ray tracing has been done, so the pressure can be applied
  call ray_trace

  !^CFG END RCM
end subroutine GM_put_from_im
