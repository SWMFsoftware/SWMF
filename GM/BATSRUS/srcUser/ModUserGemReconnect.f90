!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var

  use ModNumConst, ONLY: cTwoPi, cPi

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'GEM reconnection'

  ! GEM challenge parameters
  real:: Tp=0.01           ! plasma temperature
  real:: B0=0.0014         ! Background field
  real:: Lambda0=0.5       ! Width of current sheet
  real:: Apert = 0.2       ! amplitude of perturbation
  real:: GaussXInv = 2.0   ! X size of Gaussian perturbation in Az
  real:: GaussYInv = 2.0   ! Y size of Gaussian perturbation in Az
  real:: Kx = cTwoPi/25.6  ! X wave number of perturbation
  real:: Ky = cTwoPi/12.8  ! Y wave number of perturbation

  logical:: UseDoubleCurrentSheet = .false.
  logical:: UseSymmetric          = .false.
  logical:: UseUniformPressure    = .false.

contains

  subroutine user_read_inputs

    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModMain,      ONLY: UseUserUpdateStates, UseUserIcs, UseUserLogFiles

    real:: WaveLengthX, WaveLengthY, GaussX, GaussY

    character(len=100) :: NameCommand
    !-------------------------------------------------------------------------
    UseUserLogFiles = .true.
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#GEM')
          UseUserIcs = .true.
          call read_var('Amplitude', Apert)
          
       case('#GEMPARAM')
          call read_var('B0', B0)
          call read_var('Tp', Tp)
          call read_var('CurrentSheetWidth', Lambda0)

       case('#GEMDOUBLE')
          call read_var('UseDoubleCurrentSheet', UseDoubleCurrentSheet)
          call read_var('UseSymmetric', UseSymmetric)

       case('#GEMPRESSURE')
          call read_var('UseUniformPressure', UseUniformPressure)

       case('#GEMPERTURB')
          call read_var('GaussWidthX', GaussX)
          call read_var('GaussWidthY', GaussY)
          call read_var('WaveLengthX', WaveLengthX)
          call read_var('WaveLengthY', WaveLengthY)

          if(GaussX <= 0)then
             GaussXInv = 0.0
          else
             GaussXInv = 1.0/GaussX
          end if

          if(GaussY <= 0)then
             GaussYInv = 0.0
          else
             GaussYInv = 1.0/GaussY
          end if

          if(WaveLengthX <= 0.0)then
             Kx = 0.0
          else
             Kx = cTwoPi/WaveLengthX
          end if
          if(WaveLengthY <= 0.0)then
             Ky = 0.0
          else
             Ky = cTwoPi/WaveLengthY
          end if
       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs
  !============================================================================
  subroutine user_set_ics(iBlock)

    use ModGeometry, ONLY: Xyz_DGB

    use ModPhysics,  ONLY: ShockLeftState_V

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, rho_, Ppar_, p_, Pe_, &
         UseElectronPressure, UseAnisoPressure,Bz_
    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: iBlock

    real                :: x, y, a, Ly,Lx
    integer             :: i, j, k

    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    Ly = 12
    Lx = 24

    if (UseDoubleCurrentSheet) then
       ! Use double current sheets in a Harris equilibrium
       State_VGB(Bx_,:,:,:,iBlock) = &
            +B0*tanh((Xyz_DGB(y_,:,:,:,iBlock) + 0.25*Ly)/Lambda0) &
            -B0*tanh((Xyz_DGB(y_,:,:,:,iBlock) - 0.25*Ly)/Lambda0) &
            -B0
    else
       ! Single Harris current sheet
       State_VGB(Bx_,:,:,:,iBlock) = B0*tanh((Xyz_DGB(y_,:,:,:,iBlock))/Lambda0)
    end if

    if(UseUniformPressure)then
       ! Bz field set to make B^2 constant as specified by Ohia et. al
       State_VGB(Bz_,:,:,:,iBlock) = sqrt( &
            B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2 + ShockLeftState_V(Bz_)**2)
    else
       ! Modify thermal pressure(s) to balance magnetic pressure
       if(UseElectronPressure) then
          ! Distribute the correction proportionally between electrons and ions
          State_VGB(Pe_,:,:,:,iBlock) = ShockLeftState_V(Pe_)*(1.0 &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))
          State_VGB(p_,:,:,:,iBlock) = ShockLeftState_V(p_)*(1.0   &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))
       else
          State_VGB(p_,:,:,:,iBlock)  = ShockLeftState_V(p_) &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
       end if
       if(UseAnisoPressure) &
            ! parallel pressure
            State_VGB(Ppar_,:,:,:,iBlock) = ShockLeftState_V(Ppar_)*(1.0 &
            + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)               &
            /ShockLeftState_V(p_))
    end if

    
    ! Get density from the uniform temperature assumption
    State_VGB(rho_,:,:,:,iBlock) = State_VGB(p_,:,:,:,iBlock)/Tp
    
    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       
       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       
       if (UseDoubleCurrentSheet) then
          ! Double current sheets admit two modes; one with aligned X lines
          ! and the other with X lines separated by 180 degrees.
          if (UseSymmetric) then
             ! set intial perturbation Az = exp(-x^2/GaussX^2-y^2/Gauss^2)*cos(Kx*x)*sin(Ky*y)
             ! apply perturbation to reconnection sites only
             a = Apert*B0*exp(-x**2*GaussXInv**2 - (Ly/4 - abs(y))**2*GaussYInv**2)
             ! Bx = dAz/dy with x-component of perturbation shifted to center of grid
             State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) + &
                  a*(-2*(Ly/4 - abs(y))*GaussYInv**2*cos(Kx*x)*sin(Ky*y) + Ky*cos(Kx*x + cPi)*cos(Ky*y))
             ! By = -dAz/dx with x-component of perturbation shifted to center of grid
             State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
                  a*( 2*(Lx - abs(x))*GaussYInv**2*cos(Kx*x)*sin(Ky*y) + Kx*sin(Kx*x + cPi)*sin(Ky*y))  
          else
             ! set initial perturbation Az = exp(-x^2/GaussX^2-y^2/Gauss^2)*cos(Kx*x)
             ! apply perturbation to reconnection sites only by varying By
             if (y > 0)  then             
                a = Apert*B0*exp(-(x-Lx)**2*GaussXInv**2 - (y-Ly/4)**2*GaussYInv**2)
             else 
                a = Apert*B0*exp(-x**2*GaussXInv**2 - (y+Ly/4)**2*GaussYInv**2)
             end if
             ! By = -dAz/dx
             State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
                  a*(2*x*GaussXInv**2*cos(Kx*x)*cos(Ky*y) + Kx*sin(Kx*x))
          end if
          
       else
          ! set intial perturbation Az = exp(-x^2/GaussX^2-y^2/Gauss^2)*cos(Kx*x)*cos(Ky*y) 
          a = Apert*B0*exp(-x**2*GaussXInv**2 - y**2*GaussYInv**2)
          !  Bx = dAz/dy
          State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) + &
               a*(-2*y*GaussYInv**2*cos(Kx*x)*cos(Ky*y) - Ky*cos(Kx*x)*sin(Ky*y))
          ! By = -dAz/dx
          State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
               a*(2*x*GaussXInv**2*cos(Kx*x)*cos(Ky*y) + Kx*sin(Kx*x)*cos(Ky*y))
       end if
       
    end do; end do; end do
    
    
  end subroutine user_set_ics

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,     ONLY: nI, nJ, nK, nBlock, Unused_B
    use ModAdvance,  ONLY: By_, State_VGB
    use ModGeometry, ONLY: z2, z1
    use BATL_lib,    ONLY: CellFace_DB, CellSize_DB, Xyz_DGB

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    character (len=*), parameter :: Name='user_get_log_var'

    integer :: k1, k2, iBlock
    real:: y1, y2, dy1, dy2, HalfInvWidth, Flux
    !-------------------------------------------------------------------
    HalfInvWidth = 0.5/(z2-z1)
    VarValue=0.0
    select case(TypeVar)
    case('byflux')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          y1 = Xyz_DGB(y_,1,0,1,iBlock)
          y2 = Xyz_DGB(y_,1,nJ+1,1,iBlock)

          if(y1*y2 > 0) CYCLE
          k1 = -y1/CellSize_DB(y_,iBlock)
          k2 = k1 + 1
          dy1 = abs(Xyz_DGB(y_,1,k1,1,iBlock))/CellSize_DB(y_,iBlock)
          dy2 = 1.0 - dy1
          Flux = CellFace_DB(2,iBlock)*HalfInvWidth* &
               ( dy2*sum(abs(State_VGB(By_,1:nI,k1,1:nK,iBlock))) &
               + dy1*sum(abs(State_VGB(By_,1:nI,k2,1:nK,iBlock))))
          if(k1==0 .or. k2==nJ+1) Flux = 0.5*Flux
          VarValue = VarValue + Flux
       end do
    case default
       call stop_mpi('Unknown user logvar='//TypeVar)
    end select

  end subroutine user_get_log_var

end module ModUser
