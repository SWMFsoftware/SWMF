!^CFG COPYRIGHT UM
!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_set_ics,                    &
       IMPLEMENTED3 => user_get_log_var,                &
       IMPLEMENTED4 => user_get_b0,                     &
       IMPLEMENTED5 => user_face_bcs,                   &
       IMPLEMENTED6 => user_set_outerbcs,               &
       IMPLEMENTED7 => user_amr_criteria,               &
       IMPLEMENTED8 => user_set_plot_var

  use ModVarIndexes, ONLY: nVar

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.2
  character (len=*), parameter :: NameUserModule = &
       'Waves and GEM, Yingjuan Ma'

  character (len=20)  :: UserProblem='wave'

  real                :: Width, Amplitude, Phase, LambdaX, LambdaY, LambdaZ
  real,dimension(nVar):: Width_V=0.0, Ampl_V=0.0, Phase_V=0.0, &
       KxWave_V=0.0, KyWave_V=0.0,KzWave_V=0.0
  integer   :: iPower_V(nVar)=1
  integer   :: iVar             
  logical   :: DoInitialize=.true.
  real      :: Lx=25.6, Lz=12.8, Lambda0=0.5, Ay=0.1, Tp=0.5 , B0=1.0  

  ! The (rotated) unperturbed initial state with primitive variables
  real      :: PrimInit_V(nVar)

  ! Velocity of wave (default is set for right going whistler wave test)
  real      :: Velocity = 169.344
  
  ! Variables used by the user problem AdvectSphere                                     
  real      :: pBackgrndIo, uBackgrndIo, FlowAngle ! in XY plane         
  real      :: NumDensBackgrndIo, NumDensMaxIo
  real      :: rSphere, rSphereRs, RhoBackgrndNo, RhoMaxNo, UxNo, UyNo
  logical   :: DoCalcAnalytic = .false., DoInitSphere = .false.

contains

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
   ! use ModPhysics,  ONLY: Si2No_V, Io2Si_V,Io2No_V,&
   !      UnitRho_, UnitU_, UnitP_, UnitN_, UnitX_
    use ModNumConst,  ONLY: cTwoPi,cDegToRad
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
       case('#WAVE','#WAVE2')
          call read_var('iVar',iVar)
          call read_var('Width',Width)
          call read_var('Amplitude',Amplitude)
          call read_var('LambdaX',LambdaX)          
          call read_var('LambdaY',LambdaY)
          call read_var('LambdaZ',LambdaZ)
          call read_var('Phase',Phase)
          Width_V(iVar) = Width
          Ampl_V(iVar)  = Amplitude
          Phase_V(iVar) = Phase*cDegToRad

          if(NameCommand == '#WAVE2')then
             iPower_V(iVar) = 2
          else
             iPower_V(iVar) = 1
          end if

          !if the wavelength is smaller than 0.0, 
          !then the wave number is set to0
          KxWave_V(iVar) = max(0.0, cTwoPi/LambdaX)          
          KyWave_V(iVar) = max(0.0, cTwoPi/LambdaY)          
          KzWave_V(iVar) = max(0.0, cTwoPi/LambdaZ)

       case('#ADVECTSPHERE')
          call read_var('DoInitSphere',      DoInitSphere     )
          call read_var('NumDensBackgrndIo', NumDensBackgrndIo)
          call read_var('pBackgrndIo',       pBackgrndIo      )
          call read_var('uBackgrndIo',       uBackgrndIo      )
          call read_var('FlowAngle',         FlowAngle        )
          call read_var('rSphereRs',         rSphereRs        )
          call read_var('NumDensMaxIo',      NumDensMaxIo     )
          
       case('#ANALYTIC')
          call read_var('DoCalcAnalytic', DoCalcAnalytic)
          
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

    use ModMain,     ONLY: globalBLK, unusedBLK, TypeCoordSystem
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, r_BLK
    use ModAdvance,  ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_, Ux_, Uy_, &
         Bx_, By_, Bz_, rho_, p_, Pe_, UseElectronPressure
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: ShockSlope, ShockLeftState_V, Si2No_V, Io2Si_V,&
                           Io2No_V, UnitRho_, UnitU_, UnitP_,UnitX_, UnitN_
    use ModNumconst, ONLY: cOne,cPi, cTwoPi, cDegToRad
    use ModSize,     ONLY: nI, nJ, nK, gcn
    use ModConst,    ONLY: cProtonMass, rSun

    implicit none

    real,dimension(nVar):: state_V,KxTemp_V,KyTemp_V
    real                :: SinSlope, CosSlope, rCell
    integer             :: i, j, k, iBlock
  
    character(len=*), parameter :: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    select case(UserProblem)

    case('AdvectSphere')
       ! This case describes an IC with uniform 1D flow of plasma in the XY
       ! plane, with no density or pressure gradients and no magnetic field.
       ! A sphere with higher density is embedded in the flow, positioned at
       ! the origin. The density profiles within the sphere is given by:
       ! rho = (RhoMax-RhoBackgrnd)* cos^2(pi*r/2 rSphere) + RhoBackgrnd

       ! Convert to normalized units                                  
       ! Flow angle is measured from the x axis                                      
       UxNo = uBackgrndIo*cos(cDegToRad*FlowAngle)*Io2No_V(UnitU_)
       UyNo = uBackgrndIo*sin(cDegToRad*FlowAngle)*Io2No_V(UnitU_)
       RhoBackgrndNo = NumDensBackgrndIo*Io2Si_V(UnitN_)* &
            cProtonMass*Si2No_V(UnitRho_) 
       RhoMaxNo       = NumDensMaxIo*Io2Si_V(UnitN_)* &
            cProtonMass*Si2No_V(UnitRho_)
       ! Convert rSphereRs to normalized units (needed for OH)
       rSphere = rSphereRs*rSun*Si2No_V(UnitX_)
       !\                                                                               
       ! Start filling in cells (including ghost cells)                                 
       !/                          
       if (DoInitSphere) then
          do k= 1-gcn,nK+gcn ; do j= 1-gcn,nJ+gcn ; do i=1-gcn,nI+gcn
             rCell = R_BLK(i,j,k,iBlock)
             if (rCell .le. rSphere)then
                ! inside the sphere                                                    
                State_VGB(rho_,i,j,k,iBlock) =RhoBackgrndNo + &
                     (RhoMaxNo - RhoBackgrndNo)*(cos(0.5*cPi*rCell/rSphere))**2
             else
                ! in background flow                                                   
                State_VGB(rho_,i,j,k,iBlock) = RhoBackgrndNo
             end if
          end do; end do ; end do
       else
          State_VGB(rho_,:,:,:,iBlock) = RhoBackgrndNo
       end if
       ! velocity                                                                       
       State_VGB(RhoUx_,:,:,:,iBlock) = UxNo*State_VGB(rho_,:,:,:,iBlock)
       State_VGB(RhoUy_,:,:,:,iBlock) = UyNo*State_VGB(rho_,:,:,:,iBlock)
       State_VGB(RhoUz_, :,:,:,iBlock) = 0.0
       State_VGB(Bx_:Bz_,:,:,:,iBlock) = 0.0
       State_VGB(p_,     :,:,:,iBlock) = pBackgrndIo*Io2No_V(UnitP_)

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
               State_VGB(iVar,:,:,:,iBlock) =        &
               State_VGB(iVar,:,:,:,iBlock)          &
               + Ampl_V(iVar)*cos(Phase_V(iVar)      &
               + KxWave_V(iVar)*x_BLK(:,:,:,iBlock)  &
               + KyWave_V(iVar)*y_BLK(:,:,:,iBlock)  &
               + KzWave_V(iVar)*z_BLK(:,:,:,iBlock))**iPower_V(iVar)
       end do

    case('GEM')
       ! write(*,*)'GEM problem set up'
       State_VGB(Bx_,:,:,:,iBlock) = B0*tanh(z_BLK(:,:,:,iBlock)/Lambda0)

       ! Modify pressure(s) to balance magnetic pressure
       if(UseElectronPressure) then
          ! Distribute the correction proportionally between electrons and ions
          State_VGB(Pe_,:,:,:,iBlock) = ShockLeftState_V(Pe_)*(1.0 &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2) &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))

          State_VGB(p_,:,:,:,iBlock) = ShockLeftState_V(p_)*(1.0 &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2) &
               /(ShockLeftState_V(Pe_) + ShockLeftState_V(p_)))

       else
          State_VGB(p_,:,:,:,iBlock)  = State_VGB(p_,:,:,:,iBlock) &
               + 0.5*(B0**2 - State_VGB(Bx_,:,:,:,iBlock)**2)
       end if

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
  subroutine user_set_plot_var(iBlock,NameVar,IsDimensional,&
       PlotVar_G, PlotVarBody, UsePlotVarBody,&
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModMain,    ONLY: nI, nJ, nK
    use ModPhysics, ONLY: NameTecUnit_V, NameIdlUnit_V, UnitRho_,No2Io_V

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    real,dimension(-1:nI+2, -1:nJ+2, -1:nK+2):: RhoExact_G, RhoError_G

    character (len=*), parameter :: NameSub = 'user_set_plot_var'
    !-------------------------------------------------------------------
    IsFound = .true.

    if(DoCalcAnalytic .and. UserProblem == 'AdvectSphere' ) & 
         call calc_analytic_sln_sphere(iBlock,RhoExact_G,RhoError_G)

    select case(NameVar)
    case('rhoexact')
       if (.not. DoCalcAnalytic) then
          write(*,*) NameSub,': cannot plot rhoexact when #ANALYTIC is not set to T'
          call CON_stop('Correct PARAM.in file')
       end if
       PlotVar_G = RhoExact_G*No2Io_V(UnitRho_)
       NameTecVar = 'RhoExact'
       NameTecUnit = NameTecUnit_V(UnitRho_)
       NameIdlUnit = NameIdlUnit_V(UnitRho_)
    case('rhoerr')
       if (.not. DoCalcAnalytic) then
          write(*,*) NameSub,': cannot plot rhoerror when #ANALYTIC is not set to T'
          call CON_stop('Correct PARAM.in file')
       end if
       PlotVar_G = RhoError_G*No2Io_V(UnitRho_)
       NameTecVar = 'RhoError'
       NameTecUnit = NameTecUnit_V(UnitRho_)
       NameIdlUnit = NameIdlUnit_V(UnitRho_)
    case default
       IsFound = .false.
    end select

  contains
    subroutine calc_analytic_sln_sphere(iBlock,RhoExact_G,RhoError_G)

      use ModMain,       ONLY: time_simulation
      use ModGeometry,   ONLY: x_BLK, y_BLK, z_BLK
      use ModNumConst,   ONLY: cPi
      use ModAdvance,    ONLY: State_VGB
      use ModVarIndexes, ONLY: Rho_
      use ModPhysics,    ONLY: Si2No_V, No2Si_V, UnitX_, UnitU_

      integer,intent(in)  :: iBlock
      real,dimension(-1:nI+2,-1:nJ+2,-1:nK+2),intent(out)::RhoExact_G,&
                                                           RhoError_G
      real    :: x, y, t
      real    :: rFromCenter, xSphereCenter, ySphereCenter
      integer :: i, j, k
      !real,dimension(1:3) :: r_D, rSphereCenter_D

      character(len=*),parameter  :: NameSub = 'calc_analytic_sln_sphere'
      !-----------------------------------------------------------------
      t = time_simulation
      do k=-1,nK+2 ; do j= -1,nJ+2 ; do i= -1,nI+2
         
         x = x_BLK(i,j,k,iBlock)
         y = y_BLK(i,j,k,iBlock)
    
         ! Find current location of sphere center
         xSphereCenter = UxNo*No2Si_V(UnitU_)*t*Si2No_V(UnitX_)
         ySphereCenter = UyNo*No2Si_V(UnitU_)*t*Si2No_V(UnitX_)
         
         ! Chcek if this cell is inside the sphere
         rFromCenter = sqrt((x - xSphereCenter)**2 + (y- ySphereCenter)**2)
         if (rFromCenter .le. rSphere) then
            RhoExact_G(i,j,k) = RhoBackgrndNo + (RhoMaxNo - RhoBackgrndNo)* &
                 cos(0.5*cPi*rFromCenter/rSphere)**2
         else
            RhoExact_G(i,j,k) = RhoBackgrndNo
         end if
      end do; end do ; end do
      
      RhoError_G = RhoExact_G - State_VGB(Rho_,:,:,:,iBlock)

    end subroutine calc_analytic_sln_sphere
 
  end subroutine user_set_plot_var
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

  !============================================================================

  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    use ModSize,     ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: State_VGB, Rho_
    use ModAMR,      ONLY: RefineCritMin_I, CoarsenCritMax

    ! Variables required by this user subroutine
    integer, intent(in)          :: iBlock
    real, intent(out)            :: UserCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical ,intent(inout)       :: IsFound

    real, parameter:: RhoMin = 2.0

    integer:: i, j, k
    !------------------------------------------------------------------
    ! These settings make sure that sorted refinement works
    ! such that all blocks with criteria 1 are refined,
    ! and all blocks with criteria 0 are coarsened if possible.
    RefineCritMin_I = 0.5
    CoarsenCritMax  = 0.5

    IsFound = .true.

    ! If density exceeds RhoMin, refine
    UserCriteria = 1.0
    do k = 1, nK; do j= 1, nJ; do i = 1, nI
       if(State_VGB(Rho_,i,j,k,iBlock) > RhoMin) RETURN
    end do; end do; end do

    ! No need to refine
    UserCriteria = 0.0

  end subroutine user_amr_criteria

end module ModUser
