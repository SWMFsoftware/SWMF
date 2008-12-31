
subroutine set_test

  use ModParamRIM
  use ModRIM
  use ModNumConst

  implicit none

  ! The problem with the tests is that we are driving them from the coupler,
  ! so the tests have to be in the coupler coordinate system, which is 
  ! Latitude, longitude, with latitude starting at the north pole and
  ! longitude starting at 12 MLT.

  if (TestName == "ideal") then

     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "               IDEAL TEST"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     OuterMagJrAll = &
          2.0e-6*exp(-abs((abs(LatitudeAll)-25.0*cDegToRad)) / &
          (3.0*cDegToRad)) * sin(LongitudeAll) &
          +2.0e-6*exp(-abs((abs(LatitudeAll)-155.0*cDegToRad)) / &
          (3.0*cDegToRad)) * sin(LongitudeAll) ! &
!          + &
!          3.0e-6*exp(-abs((abs(LatitudeAll)-0.0*cDegToRad)) / &
!          (5.0*cDegToRad)) &
!          -3.0e-6*exp(-abs((abs(LatitudeAll)-180.0*cDegToRad)) / &
!          (5.0*cDegToRad)) 

     OuterMagRhoAll = sin(LatitudeAll)**3 * 2.0e-19

     OuterMagPAll = &
          1.7e-7*exp(-abs((abs(LatitudeAll)-35.0*cDegToRad)) / &
          (2.0*cDegToRad)) * abs(sin(LongitudeAll/2 + cPi/4)) &
          +1.7e-7*exp(-abs((abs(LatitudeAll)-155.0*cDegToRad)) / &
          (2.0*cDegToRad)) * abs(sin(LongitudeAll/2 + cPi/4))

     OuterMagInvBAll = 1.0
     where(LatitudeAll <  10.0*cDegToRad) OuterMagInvBAll = -1.0
     where(LatitudeAll <  10.0*cDegToRad) OuterMagPAll = 0.0
     where(LatitudeAll <  10.0*cDegToRad) OuterMagRhoAll = 0.0
     where(LatitudeAll > 170.0*cDegToRad) OuterMagInvBAll = -1.0
     where(LatitudeAll > 170.0*cDegToRad) OuterMagPAll = 0.0
     where(LatitudeAll > 170.0*cDegToRad) OuterMagRhoAll = 0.0

     where(LatitudeAll > 40.0*cDegToRad .and. &
           LatitudeAll < 150.0*cDegToRad) OuterMagRhoAll = 0.0
     where(LatitudeAll > 40.0*cDegToRad .and. &
           LatitudeAll < 150.0*cDegToRad) OuterMagPAll = 0.0
     where(LatitudeAll > 40.0*cDegToRad .and. &
           LatitudeAll < 150.0*cDegToRad) OuterMagInvBAll = -1.0

  endif

end subroutine set_test

