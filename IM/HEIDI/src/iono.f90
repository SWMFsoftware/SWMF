!=====================================
! IONO: Stand-Alone Ionosphere Model |
!=====================================
program IONO

  use ModMain
  use ModIonosphere
  implicit none
  include "numv.h"
  
  real, dimension(1:NR+3)  :: Latfac
  real, dimension(NR+3,NT) :: Jfac
  real, dimension(NR+3,NT) :: FPOT

  integer :: Ifac, iteration, si, ei
  real    :: Julian_Day

  logical :: use_amie, use_mike, use_mhd

  COMMON/FACPOT/Ifac,Latfac,Jfac,FPOT

  use_mhd  = .false.
  use_amie = .false.
  use_mike = .true.

  si = 1
  ei = 96

  ionosphere_type        = 5               ! 0=const cond. 3=Full cond.
  Iono_Radius            = (6372.0+110.0)*1000.0   ! in meters
  PolarCapPedConductance  = 0.5
  StarLightPedConductance = 1.0   !Aaron's default value: 0.5

!
! June 91
!

  f107_flux              = 230.0           ! solar min = 60; solar max = 300
  Start_Time_Array(1)    = 1991
  Start_Time_Array(2)    = 6
  Start_Time_Array(3)    = 4
  Start_Time_Array(4)    = 0
  Start_Time_Array(5)    = 0
  Start_Time_Array(6)    = 0
  Start_Time_Array(7)    = 0

!
! May 97
!

  f107_flux              = 75.0           ! solar min = 60; solar max = 300
  Start_Time_Array(1)    = 1997
  Start_Time_Array(2)    = 5
  Start_Time_Array(3)    = 14
  Start_Time_Array(4)    = 0
  Start_Time_Array(5)    = 0
  Start_Time_Array(6)    = 0
  Start_Time_Array(7)    = 0

!
! Oct 98
!

  f107_flux              = 117           ! solar min = 60; solar max = 300
  Start_Time_Array(1)    = 1998
  Start_Time_Array(2)    = 10
  Start_Time_Array(3)    = 18
  Start_Time_Array(4)    = 0
  Start_Time_Array(5)    = 0
  Start_Time_Array(6)    = 0
  Start_Time_Array(7)    = 0

!
! Sep 98
!

!  f107_flux              = 140           ! solar min = 60; solar max = 300
!  Start_Time_Array(1)    = 1998
!  Start_Time_Array(2)    = 09
!  Start_Time_Array(3)    = 24
!  Start_Time_Array(4)    = 0
!  Start_Time_Array(5)    = 0
!  Start_Time_Array(6)    = 0
!  Start_Time_Array(7)    = 0

  call increment_real_world_time(3600.0 * (si-1))

  do iteration=si,ei

     call epencalc(iteration, use_mhd, use_amie, use_mike)

     call increment_real_world_time(3600.0)

  enddo

end program IONO
