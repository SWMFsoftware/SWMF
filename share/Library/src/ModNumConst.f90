!^CFG COPYRIGHT UM
Module ModNumConst
  use ModKind
  implicit none
  real, parameter:: &
       cZero      =     0.0000000000000000000000000000000,    &
       cHalf      =     0.5000000000000000000000000000000,    &
       cThird     =     0.3333333333333333333333333333333,    &
       cQuarter   =     0.2500000000000000000000000000000,    &
       cEighth    =     0.1250000000000000000000000000000,    &
       cHundredth =     0.0100000000000000000000000000000,    &
       cTiny      =     0.0000010000000000000000000000000,    &
       cOne       =     1.0000000000000000000000000000000,    &
       cTwo       =     2.0000000000000000000000000000000,    &
       cThree     =     3.0000000000000000000000000000000,    &
       cFour      =     4.0000000000000000000000000000000,    &
       cNine	  =     9.0000000000000000000000000000000,    &
       cTen	  =    10.0000000000000000000000000000000,    &
       cTwelve	  =    12.0000000000000000000000000000000,    &
       cE1        =  cTen                                ,    &
       cHundred   =   100.0000000000000000000000000000000,    &
       cE2        =  cHundred                            ,    &
       cE3        =  1000.0000000000000000000000000000000,    &
       cE6        =  cE3*cE3,                                 &
       cE9        =  cE3*cE6,                                 &
       cE12       =  cE3*cE9,                                 &
       cE15       =  cE3*cE12,                                &
       cE18       =  cE3*cE15,                                &
       cE21       =  cE3*cE18,                                &
       cE24       =  cE3*cE21,                                &
       cE27       =  cE3*cE24,                                &
       cE30       =  cE3*cE27,                                &
       cThousand  =  cE3,                                     &
       cMillion   =  cE6,                                     &
       cHuge      =  cE18,                                    &
       cSqrtTwo   =  1.4142135623730951,                      &
       cSqrtHalf  =  cHalf*cSqrtTwo,                          &
       cPi        =  3.1415926535897932384626433832795,       &
       cTwoPi     =  cPi+cPi,                                 &
       cHalfPi    =  cHalf*cPi,                               &
       cRadToDeg  =  180.00000000000000000000/cPi,            &
       cDegToRad  =  cPi/180.00000000000000000000,            &
       cTolerance =  0.0000000001000000000000000000000000000000

  real(Real8_), parameter:: &
       cTiny8     =  0.000000000100000,                       &
       cZero8     =  0.0000000000000000000000000000000,       &
       cPi8       =  3.1415926535897932384626433832795,       &
       cTwoPi8    =  cPi8+cPi8

  real, parameter, dimension(3,3) :: cUnit_DD = reshape( &
       (/cOne,cZero,cZero, cZero,cOne,cZero, cZero,cZero,cOne/),&
       (/3,3/))

end module ModNumConst
