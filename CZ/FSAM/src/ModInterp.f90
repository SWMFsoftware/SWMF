module ModInterp

  implicit none

  private
  public :: lint
  public :: xtvd
  public :: xtvd_a
  public :: xint1d
  public :: xzc1d
  public :: xdel

contains

  !===================================================================================

  real function fq(x)
    implicit none

    real, intent(in) :: x
    real, parameter :: w = 1.D-5
    !-------------------------------------------------------------------------------

    if(x/w .gt. 10.D0) then
       fq = 1.D0
    else
       if(x/w .lt. -10.D0) then
          fq = -1.D0
       else
          fq = ( dexp(x/w) - dexp(-x/w) ) &
               / ( dexp(x/w) + dexp(-x/w) )
       endif
    endif

  end function fq

  !===================================================================================

  subroutine lint(xa,ya,n,x,y)
    ! linear interpolation at x
    implicit none

    integer, intent(in) :: n
    real, intent(in) :: xa(n),ya(n)
    real, intent(in) :: x
    real, intent(out) :: y
    real  :: h,a,b
    integer :: k, klo, khi
    !--------------------------------------------------------------------------------

    klo=1
    khi=n

    if((x.lt.xa(1)).or.(x.gt.xa(n))) then
       call con_stop('lint: x out of range!')
    endif

    do while ((khi-klo).gt.1)
       k = (khi+klo)/2
       if(xa(k).gt.x) then
          khi = k
       else
          klo = k
       endif
    enddo

    h = xa(khi) - xa(klo)
    if(h.eq.0.d0) then
       call con_stop( 'xa not strictly increasing')
    endif
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi)

  end subroutine lint

  !===================================================================================

  subroutine xtvd(iMax, ijkMax, dxa, dxbi, qint, dq, qlr)
    implicit none

    integer, intent(in) :: iMax, ijkMax
    real, intent(in)    :: dxa(iMax), dxbi(iMax)
    real, intent(in)    :: qint(ijkMax)
    real, intent(out)   :: dq(ijkMax),qlr(ijkMax)

    integer :: i, j, k
    real    :: dqi(iMax), dqim, dqip, ql, qr
    !-------------------------------------------------------------------------------

    do i=2,iMax-2
       dqim   = (qint(i  ) - qint(i-1) )*dxbi(i  )
       dqip   = (qint(i+1) - qint(i  ) )*dxbi(i+1)
       dqi(i) = sign(0.5d0, dqim)*max(0.d0, min(abs(dqim), sign(1.d0,dqim)*dqip))
    enddo

    do i=3,iMax-2
       ql     = qint(i-1) + dxa(i-1) * dqi(i-1)
       qr     = qint(i  ) - dxa(i  ) * dqi(i  )
       dq(i)  = qr - ql
       qlr(i) = 0.5d0*(qr + ql)
    enddo

  end subroutine xtvd

  !==================================================================================

  subroutine xtvd_a(iMax, ijkMax, dxa, dxai, qint, dq, qlr)
    implicit none

    integer, intent(in) :: iMax, ijkMax
    real, intent(in)    :: dxa(iMax), dxai(iMax)
    real, intent(in)    :: qint(ijkMax)
    real, intent(out)   :: dq(ijkMax), qlr(ijkMax)

    integer :: i, j, k
    real    :: dqi(iMax), dqim, dqip, ql, qr
    !--------------------------------------------------------------------------------

    do i=2,iMax-1
       dqim   = (qint(i  ) - qint(i-1))*dxai(i-1)
       dqip   = (qint(i+1) - qint(i  ))*dxai(i  )
       dqi(i) = sign(0.5d0, dqim)*max(0.d0, min(abs(dqim), sign(1.d0,dqim)*dqip))
    enddo

    do i=2,iMax-2
       ql     = qint(i  ) + dxa(i  )*dqi(i  )
       qr     = qint(i+1) - dxa(i  )*dqi(i+1)
       dq(i)  = qr - ql
       qlr(i) = 0.5d0*(qr + ql)
    enddo

  end subroutine xtvd_a

  !==================================================================================

  subroutine xdel(iMax, bint, dbint, dxa, dxbi)
    use ModPar,   ONLY: tiny, ijkn
    implicit none

    integer, intent(in)  :: iMax
    real, intent(in)  :: bint(1:ijkn)
    real, intent(out) :: dbint(1:ijkn)
    real, intent(in)  :: dxa(1:iMax), dxbi(1:iMax)

    integer :: i
    real :: dqm, dqp, q1, q2, dqi(1:iMax)
    !------------------------------------------------------------------------------ 

    do i = 2, iMax-2
       dqm    = (bint(i  ) - bint(i-1))*dxbi(i)
       dqp    = (bint(i+1) - bint(i  ))*dxbi(i+1)
       dqi(i) = max(dqm*dqp, 0.d0)*sign(1.d0, dqm + dqp)/max(abs(dqm + dqp), tiny)
    enddo

    do i = 3,iMax-2
       q1 = bint(i-1) + dxa(i-1)*dqi(i-1)
       q2 = bint(i  ) - dxa(i  )*dqi(i  )
       dbint(i) = q2 - q1
    enddo

  end subroutine xdel

  !==================================================================================

  subroutine xint1d(iMax, b, v, vfl, gfact, dxbi, dxa,  bint, vint)
    use ModPar,    ONLY: tiny, ijkn
    implicit none

    integer, intent(in) :: iMax
    real, intent(in),  dimension(1:ijkn) :: b, v, vfl
    real, intent(out), dimension(1:ijkn) :: bint, vint
    real, intent(in) :: gfact, dxbi(iMax), dxa(iMax)

    integer :: i
    real    ::  q1, q2, xi, dqm, dqp, dvi(iMax), dbi(iMax)
    !-------------------------------------------------------------------------------

    do  i = 2, iMax-2
       dqm   = (b(i  ) - b(i-1))*dxbi(i  )
       dqp   = (b(i+1) - b(i  ))*dxbi(i+1)
       dbi(i) = max(dqm*dqp, 0.d0)*sign(1.d0, dqm + dqp)/max(abs(dqm + dqp), tiny)
       dqm   = (v(i  ) - v(i-1))*dxbi(i  )
       dqp   = (v(i+1) - v(i  ))*dxbi(i+1)
       dvi(i) = max(dqm*dqp, 0.d0)*sign(1.d0, dqm + dqp)/max(abs(dqm + dqp), tiny)
    enddo

    do  i = 3, iMax-2
       q1      = b(i-1) + dxa(i-1)*dbi(i-1)
       q2      = b(i  ) - dxa(i  )*dbi(i  )
       xi      = vfl(i)*gfact
       bint(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xi*dxbi(i))*(q2 - q1)
       q1      = v(i-1) + dxa(i-1)*dvi(i-1)
       q2      = v(i  ) - dxa(i  )*dvi(i  )
       vint(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xi*dxbi(i))*(q2 - q1)
    enddo

    return

  end subroutine xint1d

  !================================================================================== 

  subroutine xzc1d(iMax, b, v, vp, vm, gfact, dxbi, dxa, bpch, bmch, vpch, vmch, DoBc)
    use ModPar,       ONLY: ijkn, tiny, myid1, nproc1
    use ModBoundary,  ONLY: niib, noib
    implicit none

    integer, intent(in) :: iMax
    real, intent(in),  dimension(1:ijkn) :: b, v, vp, vm
    real, intent(in) :: gfact, dxbi(iMax), dxa(iMax)
    real, intent(out), dimension(1:ijkn) :: bpch, bmch, vpch, vmch
    logical, intent(in), optional :: DoBc

    integer :: i
    real    :: q1, q2, xip, xim, dqm, dqp, dbi(iMax), dvi(iMax)
    !-------------------------------------------------------------------------------  

    do i = 2, iMax-2
       dqm   = (b(i  ) - b(i-1))*dxbi(i  )
       dqp   = (b(i+1) - b(i  ))*dxbi(i+1)
       dbi(i) = max(dqm*dqp, 0.d0)*sign(1.d0, dqm + dqp)/max(abs(dqm + dqp), tiny)
       dqm   = (v(i  ) - v(i-1))*dxbi(i  )
       dqp   = (v(i+1) - v(i  ))*dxbi(i+1)
       dvi(i) = max(dqm*dqp, 0.d0)*sign(1.d0, dqm + dqp)/max(abs(dqm + dqp), tiny)
    enddo

    do i = 3, iMax-2
       q1      = b(i-1) + dxa(i-1)*dbi(i-1)
       q2      = b(i  ) - dxa(i  )*dbi(i  )
       xip     = vp(i)*gfact
       bpch(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xip*dxbi(i))*(q2 - q1)
       xim     = vm(i)*gfact
       bmch(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xim*dxbi(i))*(q2 - q1)
       q1      = v(i-1) + dxa(i-1)*dvi(i-1)
       q2      = v(i  ) - dxa(i  )*dvi(i  )
       vpch(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xip*dxbi(i))*(q2 - q1)
       vmch(i) = 0.5d0*(q1 + q2) - 0.5d0*fq(xim*dxbi(i))*(q2 - q1)
    enddo

    if(present(DoBc))then
       if(DoBc)then
          if((myid1 == 0).and.(niib == 2)) then
             bmch(3) = bpch(3)
             vmch(3) = vpch(3)
          endif
          if((myid1 == nproc1-1).and.(noib == 2)) then
             bpch(iMax-2) = bmch(iMax-2)
             vpch(iMax-2) = vmch(iMax-2)
          endif
       endif
    end if

  end subroutine xzc1d

end module ModInterp
