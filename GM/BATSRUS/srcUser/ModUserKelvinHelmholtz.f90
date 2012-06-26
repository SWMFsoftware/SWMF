!^CFG COPYRIGHT UM
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_initial_perturbation

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'KELVIN-HELMHOLTZ INSTABILITY, G. Toth'

  real, parameter :: &
       xWidthUy=0.05, AmplUy=0.645, &
       xWidthUx=0.2, AmplUx=0.01, &
       yWaveUx=1.0, zWaveUx=0.0

contains

  subroutine user_initial_perturbation

    use ModMain,     ONLY: nBlock, Unused_B, ProcTest
    use ModProcMH,   ONLY: iProc
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUy_
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, y1, y2, z1, z2
    use ModNumConst, ONLY: cTwoPi
    use ModEnergy,   ONLY: calc_energy_ghost

    integer :: iBlock
    logical :: oktest, oktest_me

    !--------------------------------------------------------------------------

    if(iProc==PROCtest)then
       write(*,*)'Initializing Kelvin-Helmholtz problem'
       write(*,*)'Parameters:'
       write(*,*)'xWidthUy=',xWidthUy,' AmplUy =',AmplUy
       write(*,*)'xWidthUx=',xWidthUx,' AmplUx =',AmplUx
       write(*,*)'yWaveUx =',yWaveUx, ' zWaveUx=',zWaveUx

       call set_oktest('user_initial_perturbation',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if

    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE

       !Perturbation in Ux = 
       !    Ux0 * exp(-(x/xWidthUx)**2) * cos(ky*y) * cos(kz*z)

       where(abs(x_BLK(:,:,:,iBlock))<xWidthUx)            &
            State_VGB(RhoUx_,:,:,:,iBlock)=                   &
            AmplUx*exp(-(x_BLK(:,:,:,iBlock)/xWidthUx)**2)   &
            *cos(yWaveUx*cTwoPi/(y2-y1)*y_BLK(:,:,:,iBlock)) &
            *cos(zWaveUx*cTwoPi/(z2-z1)*z_BLK(:,:,:,iBlock)) &
            *State_VGB(Rho_,:,:,:,iBlock)

       ! Shear flow in Uy= Uy0 * tanh(x/xWidthUy)
       State_VGB(RhoUy_,:,:,:,iBlock) = &
            AmplUy*tanh(x_BLK(:,:,:,iBlock)/xWidthUy) &
            * State_VGB(Rho_,:,:,:,iBlock)

       call calc_energy_ghost(iBlock)
    end do

  end subroutine user_initial_perturbation

end module ModUser
