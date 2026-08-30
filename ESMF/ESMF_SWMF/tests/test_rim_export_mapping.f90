program test_rim_export_mapping
  use, intrinsic :: iso_fortran_env, only: real64
  use RIM_export_map, only: fill_rim_hemisphere

  implicit none

  integer, parameter :: nLat = 4, nLon = 7, nVar = 1
  real(real64), parameter :: Unset = -huge(1.0_real64)
  real :: North_II(nLat,nLon), South_II(nLat,nLon)
  real(real64) :: Data_VII(nVar,1:nLon,1:2*nLat-1)
  integer :: iLat, iLon, iPsi, iTheta, nError

  nError = 0
  do iPsi = 1, nLon
     do iTheta = 1, nLat
        North_II(iTheta,iPsi) = 1000 + 100*iTheta + iPsi
        South_II(iTheta,iPsi) = -1000 - 100*iTheta - iPsi
     end do
  end do
  Data_VII = Unset

  call fill_rim_hemisphere(North_II, Data_VII, 1, .true., &
       1, nLon, 1, 2*nLat-1, nLon, nLat)

  call check_value('north equator', Data_VII(1,4,4), 1401.0_real64)
  call check_value('north pole', Data_VII(1,4,7), 1101.0_real64)
  call check_value('north interior', Data_VII(1,4,5), 1301.0_real64)
  call check_value('north -180 seam', Data_VII(1,1,5), 1304.0_real64)
  call check_value('north +180 seam', Data_VII(1,7,5), 1304.0_real64)
  call check_value('south untouched', Data_VII(1,4,3), Unset)

  call fill_rim_hemisphere(South_II, Data_VII, 1, .false., &
       1, nLon, 1, 2*nLat-1, nLon, nLat)

  call check_value('equator remains north', Data_VII(1,4,4), 1401.0_real64)
  call check_value('south pole', Data_VII(1,4,1), -1401.0_real64)
  call check_value('south interior', Data_VII(1,4,2), -1301.0_real64)
  call check_value('south -180 seam', Data_VII(1,1,2), -1304.0_real64)
  call check_value('south +180 seam', Data_VII(1,7,2), -1304.0_real64)

  North_II = North_II + 10000.0
  South_II = South_II - 10000.0
  call fill_rim_hemisphere(North_II, Data_VII, 1, .true., &
       1, nLon, 1, 2*nLat-1, nLon, nLat)
  call fill_rim_hemisphere(South_II, Data_VII, 1, .false., &
       1, nLon, 1, 2*nLat-1, nLon, nLat)

  call check_value('second north fill', Data_VII(1,4,4), 11401.0_real64)
  call check_value('second south fill', Data_VII(1,4,1), -11401.0_real64)

  do iLat = 1, 2*nLat - 1
     if(iLat >= nLat) then
        iTheta = 2*nLat - iLat
     else
        iTheta = nLat - iLat + 1
     end if
     do iLon = 1, nLon
        iPsi = modulo(iLon + nLon/2 - 1, nLon - 1) + 1
        if(iLat >= nLat) then
           call check_value('full north mapping', Data_VII(1,iLon,iLat), &
                real(North_II(iTheta,iPsi), real64))
        else
           call check_value('full south mapping', Data_VII(1,iLon,iLat), &
                real(South_II(iTheta,iPsi), real64))
        end if
     end do
  end do

  if(nError > 0) error stop 1
  write(*,*) 'PASS: RIM export hemisphere mapping and refresh'

contains

  subroutine check_value(Name, Actual, Expected)
    character(len=*), intent(in) :: Name
    real(real64), intent(in) :: Actual, Expected

    if(Actual /= Expected) then
       nError = nError + 1
       write(*,*) 'FAIL: ', trim(Name), ' actual/expected=', Actual, Expected
    end if
  end subroutine check_value

end program test_rim_export_mapping
