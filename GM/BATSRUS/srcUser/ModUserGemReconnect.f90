!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var

  use ModNumConst, ONLY: cTwoPi

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

    use ModGeometry, ONLY: Xyz_DGB, y1, y2

    use ModPhysics,  ONLY: ShockLeftState_V

    use ModAdvance,  ONLY: State_VGB, Bx_, By_, rho_, Ppar_, p_, Pe_, &
         UseElectronPressure, UseAnisoPressure,Bz_
    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: iBlock

    real                :: x, y, a
    integer             :: i, j, k

    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    if (UseDoubleCurrentSheet) then
       ! Use double current sheets in a Harris equilibrium
       State_VGB(Bx_,:,:,:,iBlock) = &
            +B0*tanh((Xyz_DGB(y_,:,:,:,iBlock) + 0.25*(y2-y1))/Lambda0) &
            -B0*tanh((Xyz_DGB(y_,:,:,:,iBlock) - 0.25*(y2-y1))/Lambda0) &
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
          ! apply perturbation to reconnection sites only by varying By
          ! set initial perturbation Az = Apert*B0*cos(Kx*x)
          ! By = -dAz/dx
          State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + Apert*B0*Kx*sin(Kx*x)
          
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

    ! For TypeVar = byflux: 
    ! Integrate abs(By) along the current sheet at a fixed Y value
    ! Divide result by two to be compatible with GEM papers.

    use ModMain,     ONLY: nI, nJ, nK, nBlock, Unused_B
    use ModAdvance,  ONLY: By_, State_VGB
    use BATL_lib,    ONLY: CellFace_DB, CellSize_DB, Xyz_DGB, &
         CoordMin_D, CoordMax_D

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    character (len=*), parameter :: Name='user_get_log_var'

    integer :: j1, j2, iBlock
    real:: ySheet, y1, y2, Dy, dy1, dy2, HalfInvWidth, Flux
    !-------------------------------------------------------------------
    
    if(UseDoubleCurrentSheet)then
       ! Current sheet is at y=yMax/2
       ySheet = 0.5*CoordMax_D(y_)
    else
       ! Current sheet is at y=0
       ySheet = 0.0
    end if
    
    ! Width in Z direction should be ignored (it is one for 2D)
    ! The 0.5 is there to be compatible with GEM papers
    ! that did the integral for half of the domain x > 0.
    HalfInvWidth = 0.5/(CoordMax_D(z_) - CoordMin_D(z_))

    ! initialize log variable
    VarValue=0.0
    select case(TypeVar)
    case('byflux')
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          ! Check if the current sheet is in this block
          y1 = Xyz_DGB(y_,1,0,1,iBlock)
          y2 = Xyz_DGB(y_,1,nJ+1,1,iBlock)
          if( (y1 - ySheet)*(y2 - ySheet) > 0 ) CYCLE

          ! Get interpolation cells and distances
          Dy = CellSize_DB(y_,iBlock)
          j1 = (ySheet - y1)/Dy
          j2 = j1 + 1
          Dy1 = (ySheet - Xyz_DGB(y_,1,j1,1,iBlock))/Dy
          Dy2 = 1.0 - dy1

          ! Interpolate in Y, integrate in X and Z
          Flux = CellFace_DB(2,iBlock)*HalfInvWidth* &
               ( Dy2*sum(abs(State_VGB(By_,1:nI,j1,1:nK,iBlock))) &
               + Dy1*sum(abs(State_VGB(By_,1:nI,j2,1:nK,iBlock))))

          ! The flux is added up twice if ySheet is between two blocks
          if(j1==0 .or. j2==nJ+1) Flux = 0.5*Flux

          ! Add up total By flux
          VarValue = VarValue + Flux
       end do
    case default
       call stop_mpi('Unknown user logvar='//TypeVar)
    end select

  end subroutine user_get_log_var

end module ModUser
