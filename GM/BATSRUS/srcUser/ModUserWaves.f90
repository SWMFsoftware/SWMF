!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var,                &
       IMPLEMENTED4 => user_get_b0,                     &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_set_outerbcs

  use ModVarIndexes, ONLY: nVar

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: NameUserModule = &
       'Waves and GEM, Yingjuan Ma'

  character (len=20) :: UserProblem='wave'

  real :: Width, Amplitude, Phase, LambdaX, LambdaY, LambdaZ
  real,dimension(nVar):: Width_V=0.0, Ampl_V=0.0, Phase_V=0.0, &
       KxWave_V=0.0, KyWave_V=0.0,KzWave_V=0.0
  integer :: iVar             
  logical:: DoInitialize=.true.
  real :: Lx=25.6, Lz=12.8, Lambda0=0.5, Ay=0.1, Tp=0.5 , B0=1.0  

  ! The (rotated) unperturbed initial state with primitive variables
  real :: PrimInit_V(nVar)

  ! Velocity of wave (default is set for right going whistler wave test)
  real :: Velocity = 169.344
  
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
       case('#GEM')
          call read_var('Amplitude',Ay)
       case('#WAVESPEED')
          call read_var('Velocity',Velocity)
       case('#WAVE')
          call read_var('iVar',iVar)
          call read_var('Width',Width)
          call read_var('Amplitude',Amplitude)
          call read_var('LambdaX',LambdaX)          
          call read_var('LambdaY',LambdaY)
          call read_var('LambdaZ',LambdaZ)
          call read_var('Phase',Phase)
          Width_V(iVar)=Width
          Ampl_V(iVar)=Amplitude
          Phase_V(iVar)=Phase*cDegToRad
          !if the wavelength is smaller than 0.0, 
          !then the wave number is set to0
          KxWave_V(iVar) = max(0.0, cTwoPi/LambdaX)          
          KyWave_V(iVar) = max(0.0, cTwoPi/LambdaY)          
          KzWave_V(iVar) = max(0.0, cTwoPi/LambdaZ)

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

  subroutine user_set_ics
    use ModMain,     ONLY: globalBLK
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
    use ModAdvance,  ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_, Ux_, Uy_, &
         Bx_, By_, Bz_, rho_, p_
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: ShockSlope, ShockLeftState_V
    use ModNumconst, ONLY: cOne,cPi, cTwoPi
    implicit none

    real,dimension(nVar):: state_V,KxTemp_V,KyTemp_V
    real :: SinSlope, CosSlope
    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    select case(UserProblem)
    case('wave')

       if(DoInitialize)then

          DoInitialize=.false.

          PrimInit_V = ShockLeftState_V

          if(ShockSlope /= 0.0)then
             CosSlope = 1.0/sqrt(1+ShockSlope**2)
             SinSlope = ShockSlope*CosSlope

             State_V = Ampl_V
             Ampl_V(RhoUx_) = CosSlope*State_V(RhoUx_)-SinSlope*State_V(RhoUy_)
             Ampl_V(RhoUy_) = SinSlope*State_V(RhoUx_)+CosSlope*State_V(RhoUy_)
             Ampl_V(Bx_)    = CosSlope*State_V(Bx_)   - SinSlope*State_V(By_)
             Ampl_V(By_)    = SinSlope*State_V(Bx_)   + CosSlope*State_V(By_)

             KxTemp_V= KxWave_V
             KyTemp_V= KyWave_V
             KxWave_V= CosSlope*KxTemp_V - SinSlope*KyTemp_V
             KyWave_V= SinSlope*KxTemp_V + CosSlope*KyTemp_V

             State_V = ShockLeftState_V
             PrimInit_V(Ux_) = CosSlope*State_V(Ux_) - SinSlope*State_V(Uy_)
             PrimInit_V(Uy_) = SinSlope*State_V(Ux_) + CosSlope*State_V(Uy_)
             PrimInit_V(Bx_) = CosSlope*State_V(Bx_) - SinSlope*State_V(By_)
             PrimInit_V(By_) = SinSlope*State_V(Bx_) + CosSlope*State_V(By_)

             !write(*,*) &
             !    'KxWave_V(Bx_:Bz_),KyWave_V(Bx_:Bz_),KzWave_V(Bx_:Bz_)=',&
             !     KxWave_V(Bx_:Bz_),KyWave_V(Bx_:Bz_),KzWave_V(Bx_:Bz_)
             !write(*,*)'       Ampl_V(Bx_:Bz_) =',       Ampl_V(Bx_:Bz_) 
             !write(*,*)'      Phase_V(Bx_:Bz_) =',       Phase_V(Bx_:Bz_)

          end if
       end if

       do iVar=1,nVar
          where(abs( x_BLK(:,:,:,iBlock) + ShockSlope*y_BLK(:,:,:,iBlock) ) &
               < Width_V(iVar) )   &
               State_VGB(iVar,:,:,:,iBlock)=              &
               State_VGB(iVar,:,:,:,iBlock)               &
               + Ampl_V(iVar)*cos(Phase_V(iVar)           &
               + KxWave_V(iVar)*x_BLK(:,:,:,iBlock)       &
               + KyWave_V(iVar)*y_BLK(:,:,:,iBlock)       &
               + KzWave_V(iVar)*z_BLK(:,:,:,iBlock))
       end do

    case('GEM')
       ! write(*,*)'GEM problem set up'
       State_VGB(Bx_,:,:,:,iBlock) = B0*tanh(z_BLK(:,:,:,iBlock)/Lambda0)
       State_VGB(p_,:,:,:,iBlock)= State_VGB(p_,:,:,:,iBlock) &
            +0.5*(B0**2-State_VGB(Bx_,:,:,:,iBlock)**2)
       State_VGB(rho_,:,:,:,iBlock)= State_VGB(p_,:,:,:,iBlock)/Tp
       !!!set intial perturbation
       State_VGB(Bx_,:,:,:,iBlock) = State_VGB(Bx_,:,:,:,iBlock) &
            - Ay* cPi/ Lz &
            *cos(cTwoPi*x_BLK(:,:,:,iBlock)/Lx)*sin(cPi*z_BLK(:,:,:,iBlock)/Lz)
       State_VGB(Bz_,:,:,:,iBlock) = State_VGB(Bz_,:,:,:,iBlock) &
            + Ay* cTwoPi/ Lx &
            *sin(cTwoPi*x_BLK(:,:,:,iBlock)/Lx)*cos(cPi*z_BLK(:,:,:,iBlock)/Lz)
       
    case default
       if(iProc==0) call stop_mpi( &
            'user_set_ics: undefined user problem='//UserProblem)
       
    end select
  end subroutine user_set_ics

  !=====================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,     ONLY: nI, nJ, nK, nBlock, UnusedBlk
    use ModAdvance,  ONLY: Bz_, State_VGB
    use ModGeometry, ONLY: y2, y1, dx_BLK, dz_BLK, faz_BLK, z_BLK

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

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

  !=====================================================================
  subroutine user_get_b0(x, y, z, B0_D)

    real, intent(in) :: x, y, z
    real, intent(out):: B0_D(3)

    B0_D = (/0.2, 0.3, 0.4/)

  end subroutine user_get_b0

  !=====================================================================

  subroutine user_face_bcs(VarsGhostFace_V)

    use ModMain,    ONLY: x_, y_, z_, iTest, jTest, kTest, BlkTest
    use ModAdvance, ONLY: Ux_, Uy_, Uz_, By_, Bz_, State_VGB
    use ModFaceBc,  ONLY: FaceCoords_D, TimeBc, &
         VarsTrueFace_V, iFace, jFace, kFace, iBlockBc, iSide

    real, intent(out):: VarsGhostFace_V(nVar)

    integer :: iVar
    real :: Dx
    logical :: DoTest = .false.
    !-------------------------------------------------------------------------

    DoTest = iBlockBc == BlkTest
!DoTest = iFace == iTest .and. jFace == jTest .and. kFace == kTest .and. DoTest

!     if(DoTest)write(*,*)'face: iFace,jFace,kFace,iSide=',&
!          iFace,jFace,kFace,iSide
!     DoTest = .false.

    Dx = Velocity*TimeBc

!    if(DoTest) write(*,*)'Velocity, TimeBc, tSim, Dx=',&
!         Velocity, TimeBc, Dx

    do iVar = 1, nVar
       ! Both of these are primitive variables
       VarsGhostFace_V(iVar) = PrimInit_V(iVar)         &
            + Ampl_V(iVar)*cos(Phase_V(iVar)            &
            + KxWave_V(iVar)*(FaceCoords_D(x_) - Dx)    &
            + KyWave_V(iVar)*FaceCoords_D(y_)           &
            + KzWave_V(iVar)*FaceCoords_D(z_))

!       if(DoTest)write(*,*)'iVar, True, Ghost=',&
!            iVar, VarsTrueFace_V(iVar), VarsGhostFace_V(iVar)

    end do

  end subroutine user_face_bcs

  !=====================================================================

  subroutine user_set_outerbcs(iBlock,iSide, TypeBc, IsFound)

    use ModSize,     ONLY: nI, nJ, nK
    use ModPhysics,  ONLY: ShockSlope
    use ModMain,     ONLY: Time_Simulation, iTest, jTest, kTest, BlkTest
    use ModAdvance,  ONLY: nVar, Rho_, Ux_, Uz_, RhoUx_, RhoUz_, State_VGB
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, x1, x2, y1, y2, z1, z2, &
         r_BLK, XyzMin_D, XyzMax_D, TypeGeometry
    use ModSetOuterBC

    integer, intent(in) :: iBlock, iSide
    logical, intent(out) :: IsFound
    character (len=20),intent(in) :: TypeBc

    character (len=*), parameter :: Name='user_set_outerbcs'

    integer :: i,j,k,iVar
    real    :: Dx, x, y, z,r, rMin, rMax
!    logical :: DoTest = .false.
    !-------------------------------------------------------------------------

!    DoTest = iBlock == BlkTest

!    if(DoTest)then
!       write(*,*)'outer: iSide=',iSide
!       write(*,*)'x1,x2,y1,y2,z1,z2=',x1,x2,y1,y2,z1,z2
!       write(*,*)'XyzMin=',XyzMin_D
!       write(*,*)'XyzMax=',XyzMax_D
!    end if
    IsFound = .true.

    Dx = Velocity*Time_Simulation 

!Cartesian only code
!    do i = imin1g,imax2g,sign(1,imax2g-imin1g)
!       do j = jmin1g,jmax2g,sign(1,jmax2g-jmin1g)
!          do k = kmin1g,kmax2g,sign(1,kmax2g-kmin1g)

    if(TypeGeometry=='spherical_lnr')then
       rMin = exp(XyzMin_D(1)); rMax = exp(XyzMax_D(1));
    else
       rMin = XyzMin_D(1); rMax = XyzMax_D(1);
    end if

    do i=-1,nI+2
       do j=-1,nJ+2
          do k=-1,nK+2
             x = x_BLK(i,j,k,iBlk)
             y = y_BLK(i,j,k,iBlk)
             z = z_BLK(i,j,k,iBlk)
             r = r_BLK(i,j,k,iBLK)
             r = alog(r)

             if( x1<x .and. x<x2 .and. y1<y .and. y<y2 .and. z1<z .and. z<z2 &
                  .and. r > rMin .and. r < rMax) CYCLE

!             if(DoTest)write(*,*)'i,j,k,x,y,z,r=',i,j,k,x,y,z,r

             do iVar = 1, nVar

                ! Both of these are primitive variables
                State_VGB(iVar,i,j,k,iBlk) = PrimInit_V(iVar) &
                     + Ampl_V(iVar)*cos(Phase_V(iVar)               &
                     + KxWave_V(iVar)*(x - Dx)                      &
                     + KyWave_V(iVar)*y                             &
                     + KzWave_V(iVar)*z)
             end do
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlk) = &
                  State_VGB(Ux_:Uz_,i,j,k,iBlk)*State_VGB(Rho_,i,j,k,iBlk)
          end do
       end do
    end do

  end subroutine user_set_outerbcs

end module ModUser
