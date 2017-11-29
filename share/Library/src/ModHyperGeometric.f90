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
  use ModUtilities, ONLY: CON_stop

  implicit none
  real, parameter:: cTolerance_I(0:1) = (/1.0e-7, 1.0e-15/)
  real, parameter:: cEiler    = 0.57721566490
  real, parameter:: cSqrtPi   = 1.7724538509055159
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
    gamma_semi = cSqrtPi
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
    rMember   = 1.0
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
  real function hyper_semi_semi_int(nA, nB, nC, ZIn, OneMinusZIn)
    integer,        intent(in) :: nA, nB, nC 
    real, optional, intent(in) :: ZIn, OneMinusZIn
    real    :: Z, OneMinusZ
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
    real :: aPlusI, bPlusI, cPlusI, rMember, LogFactor
    character(LEN=*), parameter:: NameSub = 'hyper_semi_semi_int'
    !-----------
    if(present(ZIn))then
       Z = ZIn; OneMinusZ = 1.0 - Z
    elseif(present(OneMinusZIn))then
       OneMinusZ = OneMinusZIn; Z = 1.0 - OneMinusZ
    else
       call CON_stop(&
            NameSub//': ZIn or OneMinusZIn should be present')
    end if
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
    ! OneMinusZ = 1.0 - z
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
    elseif(n<0)then
       n = -n
       rMember      = factorial(nC-1)*factorial(n-1)/(OneMinusZ**n*&
            gamma_semi(nA)*gamma_semi(nB))
       hyper_semi_semi_int = rMember
       aPlusI       = A - (n + 1.0)
       bPlusI       = B - (n + 1.0)
       do i = 1, n-1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*real(n -i))*(-OneMinusZ)
          hyper_semi_semi_int = hyper_semi_semi_int + rMember
       end do
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = real(n)
       rMember = rMember*aPlusI*bPlusI/cPlusI*(-OneMinusZ)
       LogFactor    = -log(OneMinusZ) + psi_int(1) + psi_int(1+n) - &
            (psi_semi(nA) + psi_semi(nB))
       hyper_semi_semi_int = hyper_semi_semi_int + rMember*LogFactor
       i = 0
       do while(abs(rMember) .ge. cTolerance_I(iRealPrec))
          i = i + 1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          cPlusI = cPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*cPlusI)*OneMinusZ
          LogFactor = LogFactor + 1.0/real(i) + 1.0/cPlusI - &
               1.0/aPlusI - 1.0/bPlusI
          hyper_semi_semi_int = hyper_semi_semi_int + rMember*LogFactor
       end do
    else
       rMember      = factorial(nC-1)*factorial(n-1)/(&
            gamma_semi(nA + n)*gamma_semi(nB + n))
       hyper_semi_semi_int = rMember
       aPlusI       = A - 1.0
       bPlusI       = B - 1.0
       do i = 1, n-1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*real(n - i))*(-OneMinusZ)
          hyper_semi_semi_int = hyper_semi_semi_int + rMember
       end do
       aPlusI = aPlusI + 1.0
       bPlusI = bPlusI + 1.0
       cPlusI = real(n)
       rMember = rMember*aPlusI*bPlusI/cPlusI*(-OneMinusZ)
       LogFactor    = -log(OneMinusZ) + psi_int(1) + psi_int(1 + n) - &
            (psi_semi(nA + n) + psi_semi(nB + n))
       hyper_semi_semi_int = hyper_semi_semi_int + rMember*LogFactor
       i = 0
       do while(abs(rMember) .ge. cTolerance_I(iRealPrec))
          i = i + 1
          aPlusI = aPlusI + 1.0
          bPlusI = bPlusI + 1.0
          cPlusI = cPlusI + 1.0
          rMember = rMember*aPlusI*bPlusI/(real(i)*cPlusI)*OneMinusZ
          LogFactor = LogFactor + 1.0/real(i) + 1.0/cPlusI - &
               1.0/aPlusI - 1.0/bPlusI
          hyper_semi_semi_int = hyper_semi_semi_int + rMember*LogFactor
       end do
      
    end if
  end function hyper_semi_semi_int
  !=====================Toroid functions=========
  !\
  ! Calculate functions
  ! sqrt(2\sinh u)P^{-1}_{n - 1/2}(\cosh u)/(k^3(k^\prime)^n)
  ! sqrt(2\sinh u)Q^{-1}_{n - 1/2}(\cosh u)/(k^3(k^\prime)^n)
  ! The multiplier, k^3, is present in a definition 
  ! of the toroid functions and while calculated detivatives
  ! the multiplier is taken into account, however, it is not
  ! calculated, to avoid dividing zero by zero at k=0
  !/
  ! Herewith cosh u = (2 - k^2)/(2k^\prime);
  ! sinh u = k^2/(2k^\prime)
  ! k^\prime=sqrt(1 - k^2)=exp(-u)
  real function toroid_p(n, Kappa2In, KappaPrime2In)
    integer, intent(in):: n  !.ge.0
    real, optional, intent(in) :: Kappa2In, KappaPrime2In
    real :: Kappa2, KappaPrime2
    character(LEN=*), parameter:: NameSub = 'toroid_p'
    !-----------
    if(n < 0)call CON_stop(&
         NameSub//': argument n should be non-negative')
    if(present(Kappa2In))then
       toroid_p = 0.250*hyper_semi_semi_int(nA=1, nB=1+n, nC=3, &
            ZIn = Kappa2In)
    elseif(present(KappaPrime2In))then
       toroid_p = 0.250*hyper_semi_semi_int(nA=1, nB=1+n, nC=3, &
            OneMinusZIn = KappaPrime2In)
    else
       call CON_stop(&
            NameSub//': Kappa2In or KappaPrime2InIn should be present')
    end if
  end function toroid_p
  !====================
  real function toroid_q(n, KappaPrime2In, Kappa2In)
    integer, intent(in):: n  !.ge.1
    real, optional, intent(in) :: KappaPrime2In, Kappa2In
    character(LEN=*), parameter:: NameSub = 'toroid_q'
    !-----------
    if(n < 1)call CON_stop(&
         NameSub//': argument n should be non-negative')
    if(present(KappaPrime2In))then
       toroid_q = hyper_semi_semi_int(nA=1, nB=1+n, nC=n+1, &
            ZIn = KappaPrime2In)
    elseif(present(Kappa2In))then
       toroid_q = hyper_semi_semi_int(nA=1, nB=1+n, nC=n+1, &
            OneMinusZIn = Kappa2In)
    else
       call CON_stop(&
            NameSub//': Kappa2In or KappaPrime2InIn should be present')
    end if
    toroid_q = -toroid_q*cSqrtPi*gamma_semi(n-1)/factorial(n)
  end function toroid_q
 !====================
  real function toroid_q0(KappaPrime2In, Kappa2In)
    real, optional, intent(in) :: KappaPrime2In, Kappa2In
    character(LEN=*), parameter:: NameSub = 'toroid_q0'
    !-----------
    if(present(KappaPrime2In))then
       toroid_q0 = 2.0*cPi*hyper_semi_semi_int(nA=1, nB=1, nC=1, &
            ZIn = KappaPrime2In)
    elseif(present(Kappa2In))then
       toroid_q0 = 2.0*cPi*hyper_semi_semi_int(nA=1, nB=1, nC=1, &
            OneMinusZIn = Kappa2In)
    else
       call CON_stop(&
            NameSub//': Kappa2In or KappaPrime2InIn should be present')
    end if
  end function toroid_q0
  !====================
  real function scr_inductance(KappaPrime2)
    real, intent(in):: KappaPrime2
    !\
    ! for a superconducting ring calculate the ratio of inductance
    ! to \mu_0R_\infty as a function of (k^prime)^2. For a given R0 and a,
    ! R_\infty = sqrt(R_0^2-a^2) and k^\prime=a/(R_0+R_\infty)
    !\
    ! Inverse of inductances: from a geven harmonic harmonic and total
    real:: InvInductance, InvInductanceTotal
    real:: Tolerance
    !\
    ! Loop variable
    integer::n 
    !---------------
    InvInductance = toroid_q0(KappaPrime2In=KappaPrime2)/&
         (0.25*toroid_p(0, KappaPrime2In=KappaPrime2))
    Tolerance = InvInductance*cTolerance_I(iRealPrec)
    InvInductanceTotal = InvInductance 
    n = 0
    do while(InvInductance > Tolerance)
       n = n +1
       InvInductance = toroid_q(n, KappaPrime2In=KappaPrime2)/&
            ((0.25 - n*n)*toroid_p(n, KappaPrime2In=KappaPrime2))
       InvInductanceTotal = InvInductanceTotal + 2.0*InvInductance
    end do
    !\
    ! Invert and accunf for the common multiplier
    scr_inductance = 2.0*cPi**2/InvInductanceTotal
  end function scr_inductance
  !====================
  subroutine calc_elliptic_int_1kind(Z, KElliptic)
    real, intent(in):: Z
    real, intent(out):: KElliptic
    character(LEN=*), parameter:: NameSub = 'calc_elliptic_int_1kind'
    !----------------------------
    ! Calculate 2F1(0.5 +0, 0.5 + 0; 1; Z)
    KElliptic = 0.50*cPi*hyper_semi_semi_int(nA=0, nB=0, nC=1, ZIn=Z**2)
  end subroutine calc_elliptic_int_1kind
  !==================================================================
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
    EElliptic = 0.50*cPi*hyper_semi_semi_int(nA=-1, nB=0, nC=1, ZIn=ArgK**2)
  end subroutine calc_elliptic_int_2kind
  !===========The toroid functions=================================
  
end module ModHyperGeometric
