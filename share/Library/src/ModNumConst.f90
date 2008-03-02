!^CFG COPYRIGHT UM
Module ModNumConst
  use ModKind
  implicit none

  real, parameter:: &
       cZero      =     0.0,     &
       cHalf      =     0.5,     &
       cThird     =     1.0/3.0, &
       cQuarter   =     0.25,    &
       cEighth    =     0.125,   &
       cHundredth =     0.01,    &
       cTiny      =     0.000001,&
       cOne       =     1.0,     &
       cTwo       =     2.0,     &
       cThree     =     3.0,     &
       cFour      =     4.0,     &
       cNine	  =     9.0,     &
       cTen	  =    10.0,     &
       cTwelve	  =    12.0,     &
       cE1        =    10.0,     &
       cHundred   =   100.0,     &
       cE2        =   100.0,     &
       cE3        =  1000.0,     &
       cE6        =  1E6,        &
       cE9        =  1E9,        &
       cE12       =  1E12,       &
       cE15       =  1E15,       &
       cE18       =  1E18,       &
       cE21       =  1E21,       &
       cE24       =  1E24,       &
       cE27       =  1E27,       &
       cE30       =  1E30,       &
       cThousand  =  cE3,                                   &
       cMillion   =  cE6,                                   &
       cHuge      =  cE18,                                  &
       cSqrtTwo   =  1.4142135623730951,                    &
       cSqrtHalf  =  0.5*cSqrtTwo,                          &
       cPi        =  3.1415926535897932384626433832795,     &
       cTwoPi     =  2*cPi,                                 &
       cHalfPi    =  0.5*cPi,                               &
       cRadToDeg  =  180.0/cPi,                             &
       cDegToRad  =  cPi/180.0,                             &
       cTolerance =  0.0000000001

  real(Real8_), parameter:: &
       cTiny8     =  0.0000000001,                       &
       cZero8     =  0.0,                                &
       cPi8       =  3.1415926535897932384626433832795,  &
       cTwoPi8    =  2*cPi8

  real, parameter, dimension(3,3) :: cUnit_DD = reshape( &
       (/1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0/),&
       (/3,3/))

end module ModNumConst
