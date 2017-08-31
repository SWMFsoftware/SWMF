!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: ModHyperGeometric - calculates hypergeometric series
!INTERFACE:
module ModHyperGeometric
  use ModMpi, ONLY: &
       iRealPrec !1, if the code is compiled with double precision
  use ModNumConst, ONLY: cPi
  implicit none
  real, parameter:: cTolerance_I(0:1) = (/1.0e-7, 1.0e-15/)
  real, parameter:: cEiler    = 0.57721566490
contains
  real function psi_semi(n)
    integer, intent(in) :: n
    !\
    ! Definition of psi function: psi(z) = d\log(\Gamma(z))/dz
    ! For semiinteger z:
    ! psi(0.5 + n) = -C - log 4 +\sum_{k=1}^n{1/(k - 0.5)}
    integer :: k
    !--------
    psi_semi = -cEiler -log(4.0)
    do k=1, abs(n) 
       psi_semi = psi_semi +1.0/(real(k) - 0.50) 
    end do
  end function psi_semi
  !=================
  real function psi_int(n)
    integer, intent(in) :: n
    !\
    ! Definition of psi function: psi(z) = d\log(\Gamma(z))/dz
    ! For integer z:
    ! psi(1 + n) = -C +\sum_{k=1}^{n-1}{1.0/k}
    integer :: k
    !--------
    psi_int = -cEiler
    do k = 1, n -1 
       psi_int = psi_int +1.0/real(k) 
    end do
  end function psi_int
  !=================
  real function factorial(n)
    integer, intent(in) :: n
    ! n! = \Gamma(n+1)
    integer :: k
    !----------
    factorial = 1.0
    do k=2, n
       factorial = factorial*real(k)
    end do
  end function factorial
  !=================
  real function gamma_semi(n)
    integer, intent(in) :: n
    ! \Gamma(0.5 + n)
    integer :: k
    !------------
    gamma_semi = sqrt(cPi)
    do k = 1, n
       gamma_semi = gamma_semi * (real(k) - 0.50)
    end do
  end function gamma_semi
  !=================  
    
  !\
  ! Hypergeometric series
  real function hypergeom(A, B, C, Z)
    !\
    ! Input parameters and argument
    !/
    real, intent(in) :: A, B, C, Z
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
    i = 0
    do while(abs(rMember).ge.cTolerance_I(iRealPrec))
       i = i +1
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = cPlusI + 1.0
       rMember = rMember*aPlusI*bPlusI*z/(cPlusI*i)
       hypergeom = hypergeom + rMember 
    end do
  end function hypergeom
  !=====================
  real function hyper_semi_semi_int(nA, nB, nC, Z)
    integer, intent(in) :: nA, nB, nC 
    real,    intent(in) :: Z
    !\
    ! Calculate hypergeometric series F(a, b, c, z), if
    ! semiinteger a = 0.5 + nA, nA = 0, 1, 2
    ! semiinteger b = 0.5 + nB, nB = 0, 1, 2
    ! integer     c = nC
    real :: A, B, C
    integer :: n ! Discriminator= c - a -b
    !\
    !Loop variable
    !/
    integer :: i
    !\
    ! Misc
    !/
    real :: aPlusI, bPlusI, rMember, LogFactor, OneMinusZ
    character(LEN=*), parameter:: NameSub = 'hyper_semi_semi_int'
    !-----------
    !Real arguments of the hypergeometric function:
    A =  0.50 + real(nA); B = 0.50 + real(nB); C = real(nC)
    if (abs(z).lt.0.50) then
       !\
       !Direct summation of the hypergeometric series, well withing the
       !convergence radius
       !/
       hyper_semi_semi_int = hypergeom(&
            A=A,        &
            B=B,        &
            C=C,        &
            Z=Z)
       RETURN
    end if
    !\
    ! Use the analytic extension to the singular point z=1
    OneMinusZ = 1.0 - z
    ! The difference C - (A+B) is integer. Calculate this.
    !/
    n = nC - (1 + nA + nB)
    !\
    !The formulae for the "logarithmic case" (integer n) 
    !http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/
    !strongly depend on the sign of n. Consider case-by-case
    if(n==0)then
       LogFactor    = -log(OneMinusZ) + 2.0*psi_int(1) - &
            (psi_semi(nA) + psi_semi(nB))
       rMember      = factorial(nA + nB)/(gamma_semi(nA)*gamma_semi(nB))
       hyper_semi_semi_int = rMember*LogFactor
       aPlusI       = A - 1.0
       bPlusI       = B - 1.0
       i = 0
       do while(abs(rMember) .ge. cTolerance_I(iRealPrec))
          i = i + 1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/i**2*OneMinusZ
          LogFactor = LogFactor + 2.0/real(i) - 1.0/aPlusI - 1.0/bPlusI
          hyper_semi_semi_int = hyper_semi_semi_int + rMember*LogFactor
       end do
    else
       call CON_stop('In '//NameSub//' only n=0 case is implemented')
    end if
  end function hyper_semi_semi_int
  !=====================
  subroutine calc_elliptic_int_1kind(Z, KElliptic)
    real, intent(in):: Z
    real, intent(out):: KElliptic
    character(LEN=*), parameter:: NameSub = 'calc_elliptic_int_1kind'
    !----------------------------
    ! Calculate 2F1(0.5 +0, 0.5 + 0; 1; Z)
    KElliptic = 0.50*cPi*hyper_semi_semi_int(0, 0, 1, Z**2)
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
    if (ArgK**2.lt.0.50) then
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
