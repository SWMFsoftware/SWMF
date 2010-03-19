!^CFG COPYRIGHT UM
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_ics

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Stretched Dipole'

  real, parameter :: alpha = 2.0, beta = 3.0 ! Stretching factors

contains

  !=============================================================================
  subroutine user_set_ics

    use ModMain,       ONLY: nI, nJ, nK, globalBLK, ProcTest
    use ModProcMH,     ONLY: iProc
    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK
    use ModVarIndexes, ONLY: Bx_,By_,Bz_
    use ModPhysics,    ONLY: Bdp
    use ModB0,         ONLY: B0_DGB

    integer :: i, j, k, iBlock 
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

    iBlock = globalBLK

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)
       z = z_BLK(i,j,k,iBlock)
       !\
       ! The dipole field is stretched by alpha factor in the z direction and
       ! beta in the y direction.
       !/

       a = (sqrt(x**2 + y**2 + z**2))**5
       b = (sqrt(x**2 + (y*beta)**2 + (alpha*z)**2))**5

        !write(*,*) 'Bx', B0_DGB(1,i,j,k,iBlock), Bdp*(3*z * x)/a
        !write(*,*) 'By', B0_DGB(2,i,j,k,iBlock), Bdp*(3*z * y)/a
        !write(*,*) 'Bz', B0_DGB(3,i,j,k,iBlock), Bdp*(2*z**2  - x**2 - y**2)/a

       State_VGB(Bx_,i,j,k,iBlock) = Bdp*((3. * z * x * alpha)/b - (3. * z * x)/a)
       State_VGB(By_,i,j,k,iBlock) = Bdp*((3. * z * y * alpha)/b - (3. * z * y)/a)
       State_VGB(Bz_,i,j,k,iBlock) = &
            Bdp*((2. * (alpha*z)**2 - x**2 - (beta*y)**2)/(alpha*b) - &
            (  2.      *       z**2 - x**2 -        y**2)/a)

    end do; end do; end do

  end subroutine user_set_ics

end module ModUser
