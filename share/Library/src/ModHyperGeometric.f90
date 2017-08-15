!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: ModHyperGeometric - calculates hypergeometric series
!INTERFACE:
module ModHyperGeometric
  use ModMpi, ONLY: &
       iRealPrec !1, if the code is compiled with double precision
  use ModNumConst, ONLY: cSqrtTwo, cPi
  implicit none
  real, parameter:: cTolerance_I(0:1) = (/1.0e-7, 1.0e-15/)
  real,parameter:: pK_LIMIT1 = 0.5*cSqrtTwo
contains
  !\
  ! Hypergeometric series
  real function hypergeom(a, b, c, z)
    !\
    ! Input parameters and argument
    !/
    real, intent(in) :: a, b, c, z
    !\
    ! Loop variable
    !/
    integer:: i
    !\
    ! Misc
    !/
    real :: aPlusI, bPlusI, cPlusI, rMember
    character(LEN=*), parameter:: NameSub = 'hypergeom'
    !----------------------------
    hypergeom = 1.0
    rMember = 1.0
    aPlusI = a -1.0
    bPlusI = b -1.0
    cPlusI = c -1.0
    do i = 1, 1000
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = cPlusI + 1.0
       rMember = rMember*aPlusI*bPlusI*z/(cPlusI*i)
       hypergeom = hypergeom + rMember 
       if(abs(rMember) < cTolerance_I(iRealPrec))RETURN
    end do
    call CON_stop(&
         'In '//NameSub//' convegence with 1000 terms is not achieved')
  end function hypergeom
  !=====================
  !\
  ! Compute the complete elliptic integral of 1st kind from the series
  ! representations given in Grandstein and Ryzhik  see formulae 8.113.1 
  ! (for 0<k<0.701) and 8.113(3) (for 0.701=<k<1) therein 
  !/ 
  subroutine calc_elliptic_int_1kind(ArgK,KElliptic)
    real, intent(in):: ArgK
    real, intent(out):: KElliptic
    !---------------------------
    !\
    !Loop variable
    !/
    integer:: i
    !\
    ! Misc
    !/
    real :: aPlusI,  rMember, LogFactor, ArgKPrime2
    character(LEN=*), parameter:: NameSub = 'calc_elliptic_int_1kind'
    !----------------------------


    ! Compute the CEI of 1st kind::

    if (abs(ArgK).lt.pK_LIMIT1) then
       KElliptic = 0.50*cPi*hypergeom( 0.50, 0.50, 1.0, ArgK**2)
       RETURN
    end if
    ! Initialize some variables::
    ArgKPrime2   = 1.0 -ArgK**2
    LogFactor    = 0.5*log(16.0/ArgKPrime2)
    KElliptic    = LogFactor
    rMember      = 1.0
    aPlusI       = -0.50
    do i = 1, 1000
       aPlusI = aPlusI + 1.0
       rMember = rMember*(aPlusI/i)**2*ArgKPrime2
       LogFactor = LogFactor - 0.50/(i*(i - 0.50))
       KElliptic = KElliptic + rMember*LogFactor 
       if(abs(rMember) < cTolerance_I(iRealPrec))RETURN
    end do
    call CON_stop(&
         'In '//NameSub//' convegence with 1000 terms is not achieved')
  end subroutine calc_elliptic_int_1kind
  !====================================================================
    !------------------------------------------------------------------
    ! Compute the complete elliptic integral of 2nd kind from the series
    ! representations given in Gradstein abd Ryzhik,
    ! see formulae 8.114.1 (for 0<k<0.701) and 8.114(3)(for 0.701=<k<1)
    ! therein::
  subroutine calc_elliptic_int_2kind(ArgK,EElliptic)

    real, intent(in):: ArgK
    real, intent(out):: EElliptic
    !\
    !Loop variable
    !/
    integer:: i
    !\
    ! Misc
    !/
    real :: aPlusI, bPlusI, rMember, LogFactor, ArgKPrime2 
    real :: OneOver2n2nMinus1
    character(LEN=*), parameter:: NameSub = 'calc_elliptic_int_2kind'
    !----------------------------
    ! Compute the CEI of second kind.
    if (abs(ArgK).lt.pK_LIMIT1) then
       EElliptic = 0.50*cPi*hypergeom(-0.50, 0.50, 1.0, ArgK**2)
       RETURN
    end if
    ! Initialize some variables::
    ArgKPrime2   = 1.0 -ArgK**2
    OneOver2n2nMinus1 = 0.50
    LogFactor    = 0.5*log(16.0/ArgKPrime2) - OneOver2n2nMinus1
    rMember      = 0.50*ArgKPrime2
    EElliptic    = 1 + LogFactor*rMember
    aPlusI       = -0.50
    bPlusI       =  0.50
    do i = 2, 1000
       aPlusI  = aPlusI + 1.0
       bPlusI  = bPlusI + 1.0
       rMember = rMember*aPlusI*bPlusI/(i*(i - 1.0))*ArgKPrime2
       LogFactor         = LogFactor - OneOver2n2nMinus1
       OneOver2n2nMinus1 = 0.250/(i*(i - 0.50))
       LogFactor         = LogFactor - OneOver2n2nMinus1
       EElliptic         = EElliptic + rMember*LogFactor 
       if(abs(rMember) < cTolerance_I(iRealPrec))RETURN
    end do
    call CON_stop(&
         'In '//NameSub//' convegence with 1000 terms is not achieved')
  end subroutine calc_elliptic_int_2kind
  !=================
end module ModHyperGeometric
