program taylor

  implicit none

  ! A fast way to calculate 
  !
  ! f(B0) = \int \sqrt( B0 - B(l) ) dl 
  !
  ! for a large number of B0 values
  !
  ! The trick is to expand the square root around 0
  !
  ! \sqrt(1 - x) = \sum_{n=0}^\infty a_n x^n
  !              = 1 - x/2 - \sum_{n=2}^\infty x^n*\frac{2n-3 !!}{2n !!}
  !              
  !
  ! where n!! = n*(n-2)*(n-4)*...
  ! Once the expansion is done the primitive integrals 
  !
  ! I_n(L) = \int_0^L B^n(l) dl
  !
  ! is calculated for n = 0...N where N is a not too large number (e.g. 10)
  ! Then f(B0) is approximated as
  !
  ! f(B0) = \sqrt(B0) * \sum_{n=0}^N a_n * [I_n(L2) - I_n(L1)] / B0^n
  !
  ! where L1 and L2 are the starting and ending points with
  ! B(L1) = B(L2) = B0 and B(l) < B0 for all L1 < l < L2.

  integer, parameter :: n        =  1000 ! number of points along field lines
  integer, parameter :: nIntegral=   100 ! number of integrals to do per line
  integer, parameter :: nTaylor  =    10 ! number of terms in the Taylor series

  real,    parameter :: Length = 100.0 ! length of field line

  real :: b_I(n)                       ! magnetic field values
  real :: x_I(n)                       ! locations
  real :: Dx_I(n)                      ! cell size (for uniform grid)

  real :: bMax, bMin                   ! Maximum and minimum values of B

  real :: bStart_I(nIntegral)          ! starting B values of the integrals


  real :: IntegralSimple_I(nIntegral)  ! results with simple integration
  real :: IntegralTaylor_I(nIntegral)  ! results with Taylor expansion

  integer :: i

  integer, parameter :: nSlow=10000, nFast=10000
  real  :: Time0, Time1, Time2
  !---------------------------------------------------------------------------

  write(*,*)'Test Taylor integral starting'
  write(*,*)'Number of field line points:      ', n
  write(*,*)'Number of starting values:        ', nIntegral
  write(*,*)'Number of terms in Taylor series: ', nTaylor
  call initialize
  write(*,*)'initialized values'
  Time0 = cputime()
  write(*,*)'Starting ',nSlow,' simple integrals'
  do i=1,nSlow
     call simple_integral
  end do
  Time1 = cputime()
  write(*,*)'Starting ',nFast,' Taylor series integrals'
  do i=1,nFast
     call taylor_integral
  end do
  Time2 = cputime()

  write(*,'(i8,a,f5.2,a)') nSlow,' simple integrals: ',Time1-Time0,' s'
  write(*,'(i8,a,f5.2,a)') nFast,' Taylor integrals: ',Time2-Time1,' s'

  write(*,*)'Compare results:'
  write(*,'(a3,4a10)') 'i', 'B0', 'Simple', 'Taylor','Diff%'
  do i = 1, nIntegral, max(1, nIntegral/10)
     write(*,'(i3,4f10.3)') &
          i, bStart_I(i), IntegralSimple_I(i), IntegralTaylor_I(i), &
          100*abs(IntegralTaylor_I(i)/IntegralSimple_I(i) - 1)
  end do

contains

  subroutine initialize
    integer :: i
    !------------------------------------------------------------------------
    do i = 1, n
       ! Uniform node based grid
       !Dx_I(i) = Length/(n-1.0)
       !x_I(i)  = (i-1.0)*Length/(n-1.0)

       ! Uniform cell centered grid
       Dx_I(i) = Length/n
       x_I(i)  = (i-0.5)*Length/n

       ! Constant
       ! b_I(i) = 2.0
       ! Parabola
       b_I(i)  = (x_I(i) - Length/2)**2
       
    end do

    bMin = minval(b_I)
    bMax = maxval(b_I)

    write(*,*)'bMin, bMax=', bMin, bMax

    do i = 1, nIntegral
       ! Constant
       !bStart_I(i) = 3.0
       ! Linear from bMin to bMax
       bStart_I(i) = i*(bMax - bMin)/real(nIntegral) + bMin
    end do

    !write(*,*)'bStart_I=',bStart_I

  end subroutine initialize
  !===========================================================================
  subroutine simple_integral

    integer :: i, iIntegral, iStart
    real    :: b0, Integral
    !-------------------------------------------------------------------------

    do iIntegral = 1, nIntegral

       B0 = bStart_I(iIntegral)

       ! find starting point
       do iStart = 1, n
          if( b_I(iStart) <= B0) EXIT
       end do

       !write(*,*)'iIntegral, b0, iStart, b(iStart)=',&
       !     iIntegral, B0, iStart, b_I(iStart)

       ! do integral

       Integral = 0.0
       do i = iStart, n

          if(b_I(i) > B0) EXIT

          Integral = Integral + sqrt(B0 - b_I(i)) * Dx_I(i)

       end do

       IntegralSimple_I(iIntegral) = Integral

       !write(*,*)'iIntegral, b0, iEnd,   b(iEnd)  =',&
       !     iIntegral, B0, i, b_I(i)

    end do

  end subroutine simple_integral
  !===========================================================================
  subroutine taylor_integral

    real :: a_I(0:nTaylor)              ! coefficients in the Taylor expansion
    real :: Primitive_II(0:n,0:nTaylor) ! primitive functions for powers of B

    integer :: i, iTaylor, iStart, iEnd, iIntegral
    real    :: b, bPower, B0, InvB0, InvB0Power, Sum
    !-------------------------------------------------------------------------

    ! Calculate the Taylor coefficients for sqrt(1-x)
    a_I(0) = 1.0
    do iTaylor = 1, nTaylor
       a_I(iTaylor) = (2*iTaylor-3.)/(2*iTaylor)*a_I(iTaylor-1)
    end do

    !write(*,*)'Taylor coefficients=',a_I

    ! Calculate the primitive functions for the powers of B
    ! Initialize i = 0 
    Primitive_II(0,:) = 0.0

    do i = 1, n
       ! The zeroth power (essentially same as x_I)
       Primitive_II(i,0) = Primitive_II(i-1,0) + Dx_I(i)

       ! The first power of B
       b = b_I(i)
       bPower = b
       Primitive_II(i,1) = Primitive_II(i-1,1) + bPower*Dx_I(i)
       
       ! The higher powers of B
       do iTaylor = 2, nTaylor
          ! Get the next power
          bPower = bPower*b
          Primitive_II(i,iTaylor) = Primitive_II(i-1,iTaylor) + bPower*Dx_I(i)
       end do

    end do

    ! Now calculate the definite integrals
    do iIntegral = 1, nIntegral

       B0    = bStart_I(iIntegral)
       InvB0 = 1.0/B0

       ! find starting point
       do iStart = 1, n
          if( b_I(iStart) <= B0) EXIT
       end do

       ! find end point
       do iEnd = n, 1, -1
          if( b_I(iEnd) <= B0) EXIT
       end do
       !write(*,*)'iIntegral, b0, iStart, b(iStart)=',&
       !     iIntegral, B0, iStart, b_I(iStart)
       !write(*,*)'iIntegral, b0, iEnd,   b(iEnd)  =',&
       !     iIntegral, B0, iEnd,   b_I(iEnd)

       ! estimate integral from Taylor series

       Sum = 0.0
       InvB0Power  = 1.0
       do iTaylor = 0, nTaylor

          Sum = Sum + a_I(iTaylor)*InvB0Power* &
               (Primitive_II(iEnd,iTaylor) - Primitive_II(iStart,iTaylor))

          !write(*,*)'iTaylor, a, B0^-n =', iTaylor, a_I(iTaylor), InvB0Power
          !write(*,*)'iStart, Primitive =', iStart,Primitive_II(iStart,iTaylor)
          !write(*,*)'iEnd,   Primitive =', iEnd,  Primitive_II(iEnd,iTaylor)
          !write(*,*)'Sum =',Sum

          InvB0Power = InvB0Power*InvB0
       end do

       IntegralTaylor_I(iIntegral) = Sum*sqrt(B0)

    end do

  end subroutine taylor_integral
  !===========================================================================
  double precision function cputime()

    ! Return cputime in seconds as a double precision number.
    integer:: clock,clockrate,count_max 
    !----------------------------------------------------------------------
    call system_clock(clock,clockrate,count_max) 
    cputime=clock/dble(clockrate)               

  end function cputime
  !========================================================================

end program taylor
