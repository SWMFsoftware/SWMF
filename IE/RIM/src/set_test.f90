
subroutine set_test

  use ModParamRIM
  use ModRIM

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
          1.0e-6*exp(-abs((abs(LatitudeAll)-15.0*cDegToRad)) / &
          (5.0*cDegToRad)) * sin(LongitudeAll) &
          +1.0e-6*exp(-abs((abs(LatitudeAll)-165.0*cDegToRad)) / &
          (5.0*cDegToRad)) * sin(LongitudeAll) &
          + &
          3.0e-6*exp(-abs((abs(LatitudeAll)-0.0*cDegToRad)) / &
          (5.0*cDegToRad)) &
          -3.0e-6*exp(-abs((abs(LatitudeAll)-180.0*cDegToRad)) / &
          (5.0*cDegToRad)) 

     OuterMagRhoAll = sin(LatitudeAll) * 1.0e-20

     OuterMagPAll = &
          5.0e-9*exp(-abs((abs(LatitudeAll)-20.0*cDegToRad)) / &
          (5.0*cDegToRad)) * sin(LongitudeAll/2) &
          +5.0e-9*exp(-abs((abs(LatitudeAll)-160.0*cDegToRad)) / &
          (5.0*cDegToRad)) * sin(LongitudeAll/2)

     OuterMagInvBAll = 1.0
     where(LatitudeAll <  10.0*cDegToRad) OuterMagInvBAll = -1.0
     where(LatitudeAll <  10.0*cDegToRad) OuterMagPAll = 0.0
     where(LatitudeAll > 170.0*cDegToRad) OuterMagInvBAll = -1.0
     where(LatitudeAll > 170.0*cDegToRad) OuterMagPAll = 0.0

  endif

end subroutine set_test

