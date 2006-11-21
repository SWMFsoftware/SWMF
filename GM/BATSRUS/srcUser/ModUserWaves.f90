!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! This is the default user module which contains empty methods defined
  ! in ModUserEmpty.f90
  !
  ! Please see the documentation, and the files ModUserEmpty.f90 and 
  ! srcUser/ModUserExamples.f90 for information about what the different user
  ! subroutines do and how to implement them for your specific problem.

  use ModUserEmpty:               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var

  use ModVarIndexes, ONLY: nVar

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'Hall MHD test, Yingjuan Ma'

  character (len=20) :: UserProblem='wave'

  real :: Width, Amplitude, Phase, Lamdax, Lamday, Lamdaz
  real,dimension(nVar):: Width_I=0.0, Ampl_I=0.0, Phase_I=0.0, &
       KxWave_I=0.0, KyWave_I=0.0,KzWave_I=0.0
  integer :: iVar             
  logical:: DoInitialize=.true.
  real :: Lx=25.6, Lz=12.8, lamda0=0.5, Ay=0.1, Tp=0.5 , B0=1.0  
  
contains

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModNumConst, ONLY : cTwoPi,cDegToRad
    implicit none

    character (len=100) :: NameCommand
    !-------------------------------------------------------------------------
    
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#USERPROBLEM')
          call read_var('UserProblem',UserProblem)
       case('#WAVE')
          call read_var('iVar',iVar)
          call read_var('Width',Width)
          call read_var('Amplitude',Amplitude)
          call read_var('Lamdax',Lamdax)          
          call read_var('Lamday',Lamday)
          call read_var('Lamdaz',Lamdaz)
          call read_var('Phase',Phase)
          Width_I(iVar)=Width
          Ampl_I(iVar)=Amplitude
          Phase_I(iVar)=Phase*cDegToRad
          !if the wavelength is smaller than 0.0, 
          !then the wave number is set to0
          KxWave_I(iVar) = max(0.0, cTwoPi/Lamdax)          
          KyWave_I(iVar) = max(0.0, cTwoPi/Lamday)          
          KzWave_I(iVar) = max(0.0, cTwoPi/Lamdaz)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  !=============================================================================

  subroutine user_set_ics
    use ModMain,     ONLY: globalBLK
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModAdvance,  ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, rho_, p_
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: ShockSlope
    use ModNumconst, ONLY: cOne,cPi, cTwoPi
    implicit none

    real,dimension(nVar):: state_I,KxTemp_I,KyTemp_I
    real :: SinSlope, CosSlope
    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    select case(UserProblem)
    case('wave')
       if(ShockSlope.ne.0.0.and.DoInitialize)then
          SinSlope=ShockSlope/sqrt(cOne+ShockSlope**2)
          CosSlope=      cOne/sqrt(cOne+ShockSlope**2)
          state_I(:)=Ampl_I(:)
          Ampl_I(rhoUx_) = &
               (CosSlope*state_I(rhoUx_)-SinSlope*state_I(rhoUy_))
          Ampl_I(rhoUy_) =  &
               (SinSlope*state_I(rhoUx_)+CosSlope*state_I(rhoUy_))
          Ampl_I(Bx_) = &
               CosSlope*state_I(Bx_)-SinSlope*state_I(By_)
          Ampl_I(By_) = &
               SinSlope*state_I(Bx_)+CosSlope*state_I(By_)

          KxTemp_I= KxWave_I
          KyTemp_I= KyWave_I
          KxWave_I= CosSlope*KxTemp_I-SinSlope*KyTemp_I
          KyWave_I= SinSlope*KxTemp_I+CosSlope*KyTemp_I

          !       if(UseHallResist) call init_hall_resist

          DoInitialize=.false.

          !       write(*,*)'KxWave_I(Bx_:Bz_),KyWave_I(Bx_:Bz_),KzWave_I(Bx_:Bz_)=',&
          !            KxWave_I(Bx_:Bz_),KyWave_I(Bx_:Bz_),KzWave_I(Bx_:Bz_)
          !       write(*,*)'       Ampl_I(Bx_:Bz_) =',       Ampl_I(Bx_:Bz_) 
          !       write(*,*)'      Phase_I(Bx_:Bz_) =',       Phase_I(Bx_:Bz_)

       end if

       do iVar=1,nVar
          where(abs(x_BLK(:,:,:,iBlock))<Width_I(iVar))   &          
               State_VGB(iVar,:,:,:,iBlock)=              &
               State_VGB(iVar,:,:,:,iBlock)               &
               + Ampl_I(iVar)*cos(Phase_I(iVar)           &
               + KxWave_I(iVar)*x_BLK(:,:,:,iBlock)       &
               + KyWave_I(iVar)*y_BLK(:,:,:,iBlock)       &
               + KzWave_I(iVar)*z_BLK(:,:,:,iBlock))
       end do

    case('GEM')
!       write(*,*)'GEM problem set up'
       State_VGB(Bx_,:,:,:,iBlock) = B0*tanh(z_BLK(:,:,:,iBlock)/lamda0)
       State_VGB(p_,:,:,:,iBlock)= State_VGB(p_,:,:,:,iBlock) &
            +0.5*(B0**2-State_VGB(Bx_,:,:,:,iBlock)**2)
       State_VGB(rho_,:,:,:,iBlock)= State_VGB(p_,:,:,:,iBlock)/Tp
       !!!set intial perturbation
       State_VGB(Bx_,:,:,:,iBlock) = State_VGB(Bx_,:,:,:,iBlock)-  Ay* cPi/ Lz &
            *cos(cTwoPi*x_BLK(:,:,:,iBlock)/Lx)*sin(cPi*z_BLK(:,:,:,iBlock)/Lz)
       State_VGB(Bz_,:,:,:,iBlock) = State_VGB(Bz_,:,:,:,iBlock)+ Ay* cTwoPi/ Lx &
            *sin(cTwoPi*x_BLK(:,:,:,iBlock)/Lx)*cos(cPi*z_BLK(:,:,:,iBlock)/Lz)
       
    case default
       if(iProc==0) call stop_mpi( &
            'user_set_ics: undefined user problem='//UserProblem)
       
    end select
  end subroutine user_set_ics

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar)

    use ModMain,     ONLY: nI, nJ, nK, nBlock, UnusedBlk
    use ModAdvance,  ONLY: Bz_, State_VGB
    use ModGeometry, ONLY: y2, y1, dx_BLK, dz_BLK, faz_BLK, z_BLK

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar

    character (len=*), parameter :: Name='user_get_log_var'

    integer :: k1, k2, iBlock
    real:: z1, z2, dz1, dz2, HalfInvWidth, Flux
    !-------------------------------------------------------------------
    HalfInvWidth = 0.5/(y2-y1)
    VarValue=0.0
    select case(TypeVar)
    case('bzflux')
       do iBlock = 1, nBlock
          if(unusedBlk(iBlock)) CYCLE
          z1 = z_BLK(1,1,0,iBlock)
          z2 = z_BLK(1,1,nK+1,iBlock)

          if(z1*z2 > 0) CYCLE
          k1 = -z1/dz_BLK(iBlock)
          k2 = k1 + 1
          dz1 = abs(z_BLK(1,1,k1,iBlock))/dz_BLK(iBlock)
          dz2 = 1.0 - dz1
          Flux = faz_BLK(iBlock)*HalfInvWidth* &
               ( dz2*sum(abs(State_VGB(Bz_,1:nI,1:nJ,k1,iBlock))) &
               + dz1*sum(abs(State_VGB(Bz_,1:nI,1:nJ,k2,iBlock))))
          if(k1==0 .or. k2==nK+1) Flux = 0.5*Flux
          VarValue = VarValue + Flux
       end do
    case default
       call stop_mpi('Unknown user logvar='//TypeVar)
    end select
  end subroutine user_get_log_var


end module ModUser
