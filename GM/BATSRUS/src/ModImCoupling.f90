!^CFG COPYRIGHT UM
!^CFG FILE RCM
!==========================================================================
module ModImPressure

  integer :: isize,jsize

  real, dimension(:), allocatable :: RCM_lat, RCM_lon
  real, dimension(:,:), allocatable :: RCM_p

contains

  subroutine im_pressure_init(iSizeIn,jSizeIn)
    integer :: iSizeIn, jSizeIn, iError

    iSize = iSizeIn
    jSize = jSizeIn
    allocate(RCM_lat(iSize), RCM_lon(jSize), RCM_p(iSize,jSize), &
         stat=iError)

    RCM_p = -1.0

    if(iError/=0) call stop_mpi('GM:im_pressure_init could not allocate')

  end subroutine im_pressure_init

end module ModImPressure
!==========================================================================
subroutine get_im_pressure(iBLK,pIM)

  use ModImPressure
  use ModMain, ONLY : nI,nJ,nK,gcn
  use ModRaytrace, ONLY : ray
  use ModPhysics, ONLY : unitSI_p
  implicit none

  integer, intent(in) :: iBLK
  real, dimension(1-gcn:nI+gcn, 1-gcn:nJ+gcn, 1-gcn:nK+gcn), intent(out) &
       :: pIM

  integer :: i,j,k, n, i1,i2

  real :: rLat,rLon
  !--------------------------------------------------------------------------

  pIM = -1.

  ! If not allocated then there is no pressure obtained from IM yet
  ! Return with negative pressure, which will be ignored
  if(.not.allocated(RCM_p)) RETURN

  !\
  ! Check to see if cell centers are on closed fieldline
  !/
  if(any(ray(3,1,1:nI,1:nJ,1:nK,iBLK)==3.))then
     do i=1,nI; do j=1,nJ; do k=1,nK
        if(ray(3,1,i,j,k,iBLK)==3.) then
           rLat=ray(1,1,i,j,k,iBLK)
           rLon=ray(2,1,i,j,k,iBLK)

           !NOTE: RCM_lat in decending order
           i1=1
           do n=2,isize
              if(rLat<0.5*(RCM_lat(n-1)+RCM_lat(n)))then
                 i1=n
              end if
           end do

           !NOTE: RCM_lon in accending order
           i2=1
           do n=2,jsize
              if(rLon>0.5*(RCM_lon(n-1)+RCM_lon(n)))then
                 i2=n
              end if
           end do
           if(rLon>0.5*(RCM_lon(jsize)+RCM_lon(1)+360.)) i2=1

           pIM(i,j,k) = RCM_p(i1,i2)/unitSI_p
        end if
     end do; end do; end do
  end if

end subroutine get_im_pressure
!==========================================================================
