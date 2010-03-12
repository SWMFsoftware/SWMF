!^CFG COPYRIGHT UM
module ModUser
  
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics
  
  include 'user_module.h' !list of public methods
  
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Stretched Dipole'
  
  real, parameter :: alpha = 1.1, beta = 1./alpha
  
contains
  
  !=============================================================================
  subroutine user_set_ics

    use ModMain,       ONLY: nI, nJ, nK, globalBLK, ProcTest
    use ModProcMH,     ONLY: iProc
    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK
    use ModVarIndexes, ONLY: Bx_,By_,Bz_

    integer :: i, j, k, iBLK 
    logical :: oktest, oktest_me
    real    :: x, y, z, a, b

    !--------------------------------------------------------------------------

    if(iProc==ProcTest)then
       write(*,*)'Initializing Stretched Dipole '
       write(*,*)'Parameters:'
       write(*,*)'alpha=', alpha
       write(*,*)'beta =', beta

       call set_oktest('user_set_ics',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    iBLK = globalBLK

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       x = x_BLK(i,j,k,iBLK)
       y = y_BLK(i,j,k,iBLK)
       z = z_BLK(i,j,k,iBLK)

       a = (x**2 + y**2 + z**2)**(5./2.)
       b = (x**2 + (y*beta)**2 + (alpha*z)**2)**(5./2.)

       State_VGB(Bx_,i,j,k,iBLK) = (-3. * z * x * alpha)/b + (3. * z * x)/a
       State_VGB(By_,i,j,k,iBLK) = (-3. * z * y * alpha)/b + (3. * z * y)/a
       State_VGB(Bz_,i,j,k,iBLK) = (1./alpha *b) *      &
            (-2. * (alpha*z)**2 + x**2 + (beta*y)**2) + &
            (2. * z**2 - x**2 -y**2)/a

    end do; end do; end do

  end subroutine user_set_ics
  
end module ModUser
