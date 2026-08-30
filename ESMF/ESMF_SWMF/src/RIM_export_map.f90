module RIM_export_map

  use, intrinsic :: iso_fortran_env, only: real64

  implicit none

  private
  public :: fill_rim_hemisphere

contains

  subroutine fill_rim_hemisphere(Source_II, Data_VII, iVar, IsNorth, &
       MinLon, MaxLon, MinLat, MaxLat, nLon, nLat)
    integer, intent(in) :: iVar, MinLon, MaxLon, MinLat, MaxLat
    integer, intent(in) :: nLon, nLat
    logical, intent(in) :: IsNorth
    real, intent(in) :: Source_II(nLat,nLon)
    real(real64), intent(inout) :: Data_VII(:,MinLon:,MinLat:)

    integer :: i, j, iPsi, iTheta, jFirst, jLast

    if(IsNorth) then
       jFirst = max(MinLat, nLat)
       jLast = MaxLat
    else
       jFirst = MinLat
       jLast = min(MaxLat, nLat - 1)
    end if

    do j = jFirst, jLast
       if(IsNorth) then
          iTheta = 2*nLat - j
       else
          iTheta = nLat - j + 1
       end if
       do i = MinLon, MaxLon
          iPsi = modulo(i + nLon/2 - 1, nLon - 1) + 1
          Data_VII(iVar,i,j) = Source_II(iTheta,iPsi)
       end do
    end do
  end subroutine fill_rim_hemisphere

end module RIM_export_map
