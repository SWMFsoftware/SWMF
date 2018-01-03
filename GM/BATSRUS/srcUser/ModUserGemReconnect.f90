!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var

  use ModNumConst, ONLY: cTwoPi, cPi

  include 'user_module.h' ! list of public methods

  !\
  ! Here you must define a user routine Version number and a
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'GEM reconnection'

  ! GEM challenge parameters
  real:: Tp=0.01           ! plasma temperature in normalized unit
  real:: TpIO=0.01         ! plasma temperature in IO unit
  real:: B0=0.0014         ! Background field in normalized unit
  real:: B0IO=0.0014       ! Background field in IO unit
  real:: Lambda0=0.5       ! Width of current sheet
  real:: Apert = 0.2       ! amplitude of perturbation
  real:: GaussXInv = 2.0   ! X size of Gaussian perturbation in Az
  real:: GaussYInv = 2.0   ! Y size of Gaussian perturbation in Az
  real:: Kx = cTwoPi/25.6  ! X wave number of perturbation
  real:: Ky = cTwoPi/12.8  ! Y wave number of perturbation
  real:: ySheet = 0.0      ! Position of the current sheet

  real:: xB, Xt            ! Perturbation centers for double current sheet. 

  logical:: UseDoubleCurrentSheet = .false.
  logical:: UseUniformPressure    = .false.

  logical:: UseUniformIonPressure = .true.
  logical:: UseStandardGem        = .true.
  real   :: n0                    = 1.0

contains
  !============================================================================


  subroutine user_read_inputs

    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModMain,      ONLY: UseUserUpdateStates, UseUserIcs, UseUserLogFiles

    real:: WaveLengthX       ! X wavelength of perturbation
    real:: WaveLengthY       ! Y wavelength of perturbation
    real:: GaussX, GaussY    ! Width of Gaussian perturbation

    character(len=100) :: NameCommand
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    UseUserLogFiles = .true.
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#GEM')
          UseUserIcs = .true.
          call read_var('Amplitude', Apert)

       case('#GEMPARAM')
          call read_var('B0', B0IO)
          call read_var('Tp', TpIO)
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

       case('#UNIFORMIONPRESSURE')
          call read_var('UseUniformIonPressure', UseUniformIonPressure)

       case('#STANDARDGEM')
          call read_var('UseStandardGem', UseStandardGem)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do

    ! The two current sheets are at +ySheet and -ySheet
    ySheet = 0.0
    if(UseDoubleCurrentSheet) then
       ySheet = 0.25*WaveLengthY
       xB = -0.25*WaveLengthX
       xT =  0.25*WaveLengthX

    endif

    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================
  subroutine user_set_ics(iBlock)

    use BATL_lib,     ONLY: iTest, jTest, kTest, iProcTest, iBlockTest
    use ModGeometry, ONLY: Xyz_DGB, x1, x2, y1, y2
    use ModPhysics,  ONLY: ShockLeftState_V, ElectronCharge
    use ModAdvance,  ONLY: State_VGB, Bx_, By_, rho_, Ppar_, p_, Pe_, &
         UseElectronPressure, UseAnisoPressure, Bz_, RhoUx_, RhoUy_, RhoUz_, &
         UseEfield
    use ModSize,     ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK
    use ModMultiFluid
    use ModProcMH,         ONLY: iProc
    use ModPhysics,   ONLY: Io2No_V, UnitB_, UnitRho_, UnitP_

    integer, intent(in) :: iBlock

    real                :: x, y, a, a1, a2
    integer             :: i, j, k

    real :: nElec, nIon, uIon_D(3), uElec_D(3), Lx, Ly
    real :: Current_D(3), tmp

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    B0 = B0IO*Io2No_V(UnitB_)
    Tp = TpIO*Io2No_V(UnitP_)/Io2No_V(UnitRho_)

    if (UseDoubleCurrentSheet) then
       ! Use double current sheets in a Harris equilibrium
       State_VGB(Bx_,:,:,:,iBlock) = B0* &
            ( tanh((Xyz_DGB(y_,:,:,:,iBlock) + ySheet)/Lambda0)    &
            - tanh((Xyz_DGB(y_,:,:,:,iBlock) - ySheet)/Lambda0) - 1)
    else
       ! Single Harris current sheet
       State_VGB(Bx_,:,:,:,iBlock) = B0* &
            tanh(Xyz_DGB(y_,:,:,:,iBlock)/Lambda0)       
    end if

    if(UseUniformPressure)then
       ! Bz field set to make B^2 constant as specified by Ohia et. al
       State_VGB(Bz_,:,:,:,iBlock) = sqrt( &
            B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2 + ShockLeftState_V(Bz_)**2)
    else
       ! Modify thermal pressure(s) to balance magnetic pressure
       if(UseElectronPressure .and. .not. UseUniformIonPressure) then
          ! Distribute the correction proportionally between electrons and ions
          State_VGB(Pe_,:,:,:,iBlock) = ShockLeftState_V(Pe_)*(1.0 &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))
          State_VGB(p_,:,:,:,iBlock) = ShockLeftState_V(p_)*(1.0   &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)      &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))
       else if(UseElectronPressure .and. UseUniformIonPressure) then
          State_VGB(Pe_,:,:,:,iBlock) = ShockLeftState_V(Pe_)      &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
       else if (.not. UseEfield) then
          State_VGB(p_,:,:,:,iBlock)  = ShockLeftState_V(p_) &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
       end if
       if(UseAnisoPressure .and. .not. UseEfield) &
                                ! parallel pressure
            State_VGB(Ppar_,:,:,:,iBlock) = ShockLeftState_V(Ppar_)*(1.0 &
            + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)               &
            /ShockLeftState_V(p_))

       if (UseEfield) then
          if(UseUniformIonPressure) then
             State_VGB(iPIon_I(ElectronFirst_),:,:,:,iBlock) =   &
                  ShockLeftState_V(iPIon_I(ElectronFirst_))      &
                  + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
             if (UseAnisoPressure) then
                State_VGB(iPIon_I(ElectronFirst_),:,:,:,iBlock) =         &
                     ShockLeftState_V(iPIon_I(ElectronFirst_))            &
                     + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
                State_VGB(iPparIon_I(ElectronFirst_),:,:,:,iBlock) =      &
                     State_VGB(iPIon_I(ElectronFirst_),:,:,:,iBlock)
             end if
          else
             State_VGB(iPIon_I(ElectronFirst_),:,:,:,iBlock) =        &
                  ShockLeftState_V(iPIon_I(ElectronFirst_))*(1.0      &
                  + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2) &
                  /(ShockLeftState_V(iPIon_I(ElectronFirst_)) +       &
                  ShockLeftState_V(p_)))

             State_VGB(p_,:,:,:,iBlock) = ShockLeftState_V(p_)*(1.0 &
                  + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)    &
                  /(ShockLeftState_V(iPIon_I(ElectronFirst_))            &
                  + ShockLeftState_V(p_)))
          end if
       end if
    end if

    ! Get density from the uniform temperature assumption
    State_VGB(Rho_,:,:,:,iBlock) = State_VGB(p_,:,:,:,iBlock)/Tp

    if (UseEfield) then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          y = Xyz_DGB(y_,i,j,k,iBlock)

          nIon   = State_VGB(Rho_,i,j,k,iBlock)/MassFluid_I(IonFirst_)
          uIon_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)   &
               /State_VGB(Rho_,i,j,k,iBlock)

          Current_D     = 0

          if(UseDoubleCurrentSheet) then
             Current_D(z_) = -B0/Lambda0*(  &
                  1.0/(cosh((y+ySheet)/Lambda0))**2 -  &
                  1.0/(cosh((y-ySheet)/Lambda0))**2)
          else
             Current_D(z_) = -B0/Lambda0/(cosh(y/Lambda0))**2
          endif

          nElec   = nIon
          uElec_D = uIon_D - Current_D/nElec/ElectronCharge

          State_VGB(Ex_:Ez_,i,j,k,iBlock) = 0.0

          State_VGB(iRhoIon_I(ElectronFirst_),i,j,k,iBlock)   =  &
               nIon*MassIon_I(ElectronFirst_)
          State_VGB(iRhoUxIon_I(ElectronFirst_),i,j,k,iBlock) =  &
               State_VGB(iRhoIon_I(ElectronFirst_),i,j,k,iBlock) &
               *uElec_D(x_)
          State_VGB(iRhoUyIon_I(ElectronFirst_),i,j,k,iBlock) =  &
               State_VGB(iRhoIon_I(ElectronFirst_),i,j,k,iBlock) &
               *uElec_D(y_)
          State_VGB(iRhoUzIon_I(ElectronFirst_),i,j,k,iBlock) =  &
               State_VGB(iRhoIon_I(ElectronFirst_),i,j,k,iBlock) &
               *uElec_D(z_)

          if (i == iTest .and. j == jTest .and. k == kTest .and. &
               iBlock == iBlockTest .and. iProc == iProcTest) then
             write(*,*) 'Xyz_DGB        =', Xyz_DGB(:,i,j,k,iBlock)
             write(*,*) 'uIon_D         =', uIon_D
             write(*,*) 'uElec_D        =', uElec_D
             write(*,*) 'ElectronCharge =', ElectronCharge
          end if
       end do; end do; end do
    end if

    ! Size of the box
    Lx = x2 - x1
    Ly = y2 - y1

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)

       if (UseDoubleCurrentSheet) then
          ! Az = -exp(-(x-xT)^2/GaussX^2-(y-ysheet)^2/Gauss^2)*
          !      cos(Kx*(x-xT))*cos(Ky*(y-ysheet))  +
          !      exp(-(x-xB)^2/GaussX^2-(y+ysheet)^2/Gauss^2)*
          !      cos(Kx*(x-xB))*cos(Ky*(y+ysheet))
          a1 = -1*(Apert*B0*exp(-(x-xT)**2*GaussXInv**2 - (y-ySheet)**2*GaussYInv**2))
          a2 = (Apert*B0*exp(-(x-xB)**2*GaussXInv**2 - (y+ySheet)**2*GaussYInv**2))
          !  Bx = dAz/dy
          State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) + &
               a1*(-2*(y-ySheet)*GaussYInv**2*cos(Kx*(x-xT))*cos(Ky*(y-ySheet)) &
               - Ky*cos(Kx*(x-xT))*sin(Ky*(y-ySheet))) + &
               a2*(-2*(y+ySheet)*GaussYInv**2*cos(Kx*(x-xB))*cos(Ky*(y+ySheet)) &
               - Ky*cos(Kx*(x-xB))*sin(Ky*(y+ySheet)))

          ! By = -dAz/dx
          State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
               a1*(2*(x-xT)*GaussXInv**2*cos(Kx*(x-xT))*cos(Ky*(y-ySheet)) &
               + Kx*sin(Kx*(x-xT))*cos(Ky*(y-ySheet))) + &
               a2*(2*(x-xB)*GaussXInv**2*cos(Kx*(x-xB))*cos(Ky*(y+ySheet)) &
               + Kx*sin(Kx*(x-xB))*cos(Ky*(y+ySheet)))


       else if (UseStandardGEM) then
          State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) &
               - Apert * B0 * cPi/Ly *cos(cTwoPi*x/Lx) * sin(cPi*y/Ly)
          State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) &
               + Apert * B0 * cTwoPi/Lx * sin(cTwoPi*x/Lx) * cos(cPi*y/Ly)
       else
          ! Az = exp(-x^2/GaussX^2-y^2/Gauss^2)*cos(Kx*x)*cos(Ky*y) 
          a = Apert*B0*exp(-x**2*GaussXInv**2 - y**2*GaussYInv**2)
          !  Bx = dAz/dy
          State_VGB(Bx_,i,j,k,iBlock) = State_VGB(Bx_,i,j,k,iBlock) + &
               a*(-2*y*GaussYInv**2*cos(Kx*x)*cos(Ky*y) &
               - Ky*cos(Kx*x)*sin(Ky*y))
          ! By = -dAz/dx
          State_VGB(By_,i,j,k,iBlock) = State_VGB(By_,i,j,k,iBlock) + &
               a*(2*x*GaussXInv**2*cos(Kx*x)*cos(Ky*y) &
               + Kx*sin(Kx*x)*cos(Ky*y))
       end if

    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================

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
    real:: y1, y2, Dy, dy1, dy2, HalfInvWidth, Flux
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_get_log_var'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
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

    call test_stop(NameSub, DoTest)
  end subroutine user_get_log_var
  !============================================================================

end module ModUser
!==============================================================================
