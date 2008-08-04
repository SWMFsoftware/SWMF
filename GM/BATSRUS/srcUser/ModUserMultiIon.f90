!^CFG COPYRIGHT UM
!========================================================================
module ModUser
  ! do electro-magnetic force term in multi-ion point-implicitly

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_calc_sources,               &
       IMPLEMENTED3 => user_init_point_implicit,        &
       IMPLEMENTED4 => user_set_ics,                    &
       IMPLEMENTED5 => user_face_bcs

  use ModMultiFluid

  use ModNumConst, ONLY: cSqrtHalf, cDegToRad

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'POINT IMPLICIT MULTI-ION SOURCE, Toth and Ma'


  character (len=20) :: UserProblem='wave'
  real    :: LatitudeDim   = 45.0, SinLatitude = cSqrtHalf
  real    :: uOutDim = 2.0
  real    :: CollisionCoefDim = 1.0, CollisionCoef
  logical :: IsAnalytic = .true.

  real :: Width, Amplitude, Phase, Lamdax, Lamday, Lamdaz
  real,dimension(nVar):: Width_I=0.0, Ampl_I=0.0, Phase_I=0.0, &
       KxWave_I=0.0, KyWave_I=0.0,KzWave_I=0.0
  integer :: iVar             
  logical:: DoInitialize=.true.
  real :: Lx=25.6, Lz=12.8, lamda0=0.5, Ay=0.1, Tp=0.5 , B0=1.0  
  
contains

  !==========================================================================
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

       case("#FRICTION")
          call read_var('CollisionCoefDim',CollisionCoefDim)
       case('#DSDU')
          call read_var('IsAnalytic',IsAnalytic)
       case('#OUTFLOW')
          call read_var('uOutDim',uOutDim)
          call read_var('Latitude',LatitudeDim)
          SinLatitude = sin(LatitudeDim*cDegToRad)
       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do

  end subroutine user_read_inputs
  !==========================================================================
  subroutine user_calc_sources

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModProcMH,  ONLY: iProc
    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, Test_String, BlkTest, ProcTest
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModAdvance, ONLY: B0_DGB
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB
    use ModPhysics, ONLY: ElectronCharge, gm1, inv_gm1, &
         Si2No_V, No2Si_V, UnitTemperature_, UnitT_
    use ModMain,    ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: cross_product

    ! Variables for multi-ion MHD
    real    :: InvCharge, NumDens, InvNumDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, uIon2_D, u_D, uPlus_D, uPlusHallU_D
    real, dimension(3) :: Current_D, Force_D
    real, dimension(nIonFluid) :: NumDens_I, InvRho_I, Ux_I, Uy_I, Uz_I, Temp_I

    integer :: iBlock, i, j, k, jFluid, iFirstIons
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, TemperatureCoef, Heating
    real :: CollisionRate_II(nIonFluid, nIonFluid), CollisionRate

    character (len=*), parameter :: NameSub = 'user_calc_sources'
    logical :: DoTest, DoTestMe
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    ! Add source term n_s*(- u_+ - w_H + u_s )xB for multi-ions
    ! where u_+ is the number density weighted average ion velocity,
    ! and w_H = -J/(e n_e) is the Hall velocity. Here
    ! e is the electron charge and n_e is the electron number density.

    InvCharge = 1.0/ElectronCharge

    ! Rate = n*CoefDim / T^1.5 with T [K], n [/cc] and Rate [1/s]
    CollisionCoef = CollisionCoefDim &
         /No2Si_V(UnitTemperature_)**1.5/Si2No_V(UnitT_)

    do jFluid = 1, nIonFluid
       do iFluid = 1, nIonFluid
          CollisionRate_II(iFluid, jFluid) = CollisionCoef* &
               MassFluid_I(iFluid)*MassFluid_I(jFluid) &
               /(MassFluid_I(iFluid)+MassFluid_I(jFluid))
       end do
    end do

    ! Do not add
    iFirstIons = 1
    if(TypeFluid_I(1) == 'ion')iFirstIons = 2

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,globalBLK)

       if(TypeFluid_I(1) == 'ion')then
          ! Get first fluid quantities
          State_V(Rho_) = State_V(Rho_) &
               - sum(State_V(iRhoIon_I(2:nIonFluid)))
          State_V(RhoUx_) = State_V(RhoUx_) &
               - sum(State_V(iRhoUxIon_I(2:nIonFluid)))
          State_V(RhoUy_) = State_V(RhoUy_) &
               - sum(State_V(iRhoUyIon_I(2:nIonFluid)))
          State_V(RhoUz_) = State_V(RhoUz_) &
               - sum(State_V(iRhoUzIon_I(2:nIonFluid)))
          State_V(P_) = State_V(P_) &
               - sum(State_V(iPIon_I(2:nIonFluid)))
       end if

       ! Total magnetic field
       FullB_D = State_V(Bx_:Bz_) + B0_DGB(:,i,j,k,globalBLK)
       ! calculate number densities
       NumDens_I  = State_V(iRhoIon_I) / MassFluid_I(1:nIonFluid)
       NumDens    = sum(NumDens_I)
       InvNumDens = 1.0/NumDens

       Temp_I     = State_V(iPIon_I)/NumDens_I
       AverageTemp= sum(State_V(iPIon_I))*InvNumDens

       InvRho_I = 1.0/State_V(iRhoIon_I)
       Ux_I  = InvRho_I*State_V(iUxIon_I)
       Uy_I  = InvRho_I*State_V(iUyIon_I)
       Uz_I  = InvRho_I*State_V(iUzIon_I)

       ! calculate the average positive charge velocity
       uPlus_D(x_) = InvNumDens* sum(NumDens_I*Ux_I)
       uPlus_D(y_) = InvNumDens* sum(NumDens_I*Uy_I)
       uPlus_D(z_) = InvNumDens* sum(NumDens_I*Uz_I)

       ! Add the Hall velocity -J/(e n)
       if(index(Test_String,'newj') > 0)then
          Current_D = vInv_CB(i,j,k,globalBLK)*&
               ( bCrossArea_DX(:,i+1,j,k) - bCrossArea_DX(:,i,j,k) &
               + bCrossArea_DY(:,i,j+1,k) - bCrossArea_DY(:,i,j,k) &
               + bCrossArea_DZ(:,i,j,k+1) - bCrossArea_DZ(:,i,j,k))
       else
          call get_current(i,j,k,GlobalBlk,Current_D)
       end if
       uPlusHallU_D = uPlus_D - InvNumDens*InvCharge*Current_D

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       CollisionRate = CollisionCoef

       ! Calculate the source term for all the ion fluids
       do iFluid = iFirstIons, nIonFluid
          call select_fluid
          uIon_D = (/ Ux_I(iFLuid),  Uy_I(iFluid), Uz_I(iFluid) /)
          u_D    = uIon_D - uPlusHallU_D

          Force_D = &
               ElectronCharge*NumDens_I(iFluid)*cross_product(u_D, FullB_D) 

          Heating = 0.0

          do jFluid = 1, nIonFluid
             if(jFluid == iFluid) CYCLE

             ! Add collisional term
             uIon2_D = (/ Ux_I(jFLuid),  Uy_I(jFluid), Uz_I(jFluid) /)

             if(  InvNumDens*NumDens_I(iFluid) < 0.01 .or. &
                  InvNumDens*NumDens_I(jFluid) < 0.01) then
                CollisionRate = 100.0 * NumDens_I(iFluid) * NumDens_I(jFluid)
             else
                CollisionRate = CollisionRate_II(iFluid, jFluid) * &
                     NumDens_I(iFluid) * NumDens_I(jFluid) &
                     * TemperatureCoef
             end if

             Force_D = Force_D + CollisionRate*(uIon2_D - uIon_D)

             Heating = Heating + CollisionRate* &
                  ( 2*(Temp_I(jFluid) - Temp_I(iFluid)) &
                  + gm1*sum((uIon2_D - uIon_D)**2) )
          end do

          Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) = &
               Source_VC(iRhoUx_I(iFluid):iRhoUz_I(iFluid),i,j,k) + Force_D

          Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) + Heating

          Source_VC(Energy_-1+iFluid,i,j,k) = &
               Source_VC(Energy_-1+iFluid,i,j,k) + sum(Force_D*uIon_D) &
               + inv_gm1*Heating

          if(.not.IsAnalytic) CYCLE

          Coef   = ElectronCharge/MassFluid_I(iFluid) &
               *(1.0 - InvNumDens*NumDens_I(iFluid))
          CoefBx = Coef*FullB_D(x_)
          CoefBy = Coef*FullB_D(y_)
          CoefBz = Coef*FullB_D(z_)
          DsDu_VVC(iRhoUx, iRhoUy, i,j,k) =   CoefBz
          DsDu_VVC(iRhoUx, iRhoUz, i,j,k) = - CoefBy
          DsDu_VVC(iRhoUy, iRhoUz, i,j,k) =   CoefBx
          DsDu_VVC(iRhoUy, iRhoUx, i,j,k) = - CoefBz
          DsDu_VVC(iRhoUz, iRhoUx, i,j,k) =   CoefBy
          DsDu_VVC(iRhoUz, iRhoUy, i,j,k) = - CoefBx

          Coef = 0.5*CollisionRate*(NumDens - NumDens_I(iFluid))
          DsDu_VVC(iP,    iP,     i,j,k) = -2*Coef
          DsDu_VVC(iRhoUx,iRhoUx, i,j,k) = -Coef
          DsDu_VVC(iRhoUy,iRhoUy, i,j,k) = -Coef
          DsDu_VVC(iRhoUz,iRhoUz, i,j,k) = -Coef

       end do
       if( .not. IsAnalytic) CYCLE
       if(iFirstIons == 1)then
          Coef   = -ElectronCharge/MassFluid_I(1) &
               *InvNumDens*NumDens_I(1)
          CoefBx = Coef*FullB_D(x_)
          CoefBy = Coef*FullB_D(y_)
          CoefBz = Coef*FullB_D(z_)
          DsDu_VVC(RhoUx_, OpRhoUy_, i,j,k) =   CoefBz
          DsDu_VVC(RhoUx_, OpRhoUz_, i,j,k) = - CoefBy
          DsDu_VVC(RhoUy_, OpRhoUz_, i,j,k) =   CoefBx
          DsDu_VVC(RhoUy_, OpRhoUx_, i,j,k) = - CoefBz
          DsDu_VVC(RhoUz_, OpRhoUx_, i,j,k) =   CoefBy
          DsDu_VVC(RhoUz_, OpRhoUy_, i,j,k) = - CoefBx

          Coef = -0.5*CollisionRate*NumDens_I(1)
          DsDu_VVC(P_,     OpP_,     i, j, k) = -2*Coef
          DsDu_VVC(RhoUx_, OpRhoUx_, i, j, k) = -Coef
          DsDu_VVC(RhoUy_, OpRhoUy_, i, j, k) = -Coef
          DsDu_VVC(RhoUz_, OpRhoUz_, i, j, k) = -Coef

          Coef = -0.5*CollisionRate*NumDens_I(2)
          DsDu_VVC(OpP_,     P_,     i, j, k) = -2*Coef
          DsDu_VVC(OpRhoUx_, RhoUx_, i, j, k) = -Coef
          DsDu_VVC(OpRhoUy_, RhoUy_, i, j, k) = -Coef
          DsDu_VVC(OpRhoUz_, RhoUz_, i, j, k) = -Coef
       end if

       Coef   = -ElectronCharge/MassFluid_I(2) &
            *InvNumDens*NumDens_I(2)
       CoefBx = Coef*FullB_D(x_)
       CoefBy = Coef*FullB_D(y_)
       CoefBz = Coef*FullB_D(z_)
       DsDu_VVC(OpRhoUx_, RhoUy_, i,j,k) =   CoefBz
       DsDu_VVC(OpRhoUx_, RhoUz_, i,j,k) = - CoefBy
       DsDu_VVC(OpRhoUy_, RhoUz_, i,j,k) =   CoefBx
       DsDu_VVC(OpRhoUy_, RhoUx_, i,j,k) = - CoefBz
       DsDu_VVC(OpRhoUz_, RhoUx_, i,j,k) =   CoefBy
       DsDu_VVC(OpRhoUz_, RhoUy_, i,j,k) = - CoefBx

    end do; end do; end do

    if(DoTestMe)then
       write(*,*)NameSub,' CollisionCoef=',CollisionCoef
       write(*,*)NameSub,' CollisionRate=',CollisionRate
       write(*,*)NameSub,' AverageTemp  =',AverageTemp
       write(*,*)NameSub,' AverageTempDim=', &
            AverageTemp*No2Si_V(UnitTemperature_)
    end if

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !------------------------------------------------------------------------

    ! All ion momenta are implicit
    if(TypeFluid_I(1) == 'ions')then
       allocate(iVarPointImpl_I(4*nIonFluid))

       do iFluid = 1, nIonFluid
          iVarPointImpl_I(4*iFluid-3) = iRhoUx_I(iFluid)
          iVarPointImpl_I(4*iFluid-2) = iRhoUy_I(iFluid)
          iVarPointImpl_I(4*iFluid-1) = iRhoUz_I(iFluid)
          iVarPointImpl_I(4*iFluid)   = iP_I(iFluid)
       end do
    else
       allocate(iVarPointImpl_I(4*(nIonFluid-1)))
       do iFluid = 1, nIonFluid-1
          iVarPointImpl_I(4*iFluid-3) = iRhoUx_I(iFluid+1)
          iVarPointImpl_I(4*iFluid-2) = iRhoUy_I(iFluid+1)
          iVarPointImpl_I(4*iFluid-1) = iRhoUz_I(iFluid+1)
          iVarPointImpl_I(4*iFluid)   = iP_I(iFluid+1)
       end do
    end if

    IsPointImplMatrixSet = IsAnalytic

  end subroutine user_init_point_implicit

  !=====================================================================
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
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModMain, ONLY: x_, y_, z_
    use ModVarIndexes
    use ModPhysics, ONLY: BodyRho_I, BodyP_I, Io2No_V, UnitU_
    use ModFaceBc,  ONLY: FaceCoords_D, B0Face_D, VarsTrueFace_V

    real, intent(out):: VarsGhostFace_V(nVar)

    real :: zMin, uPerB0
    !--------------------------------------------------------------------------

    zMin = sqrt(sum(FaceCoords_D**2))*SinLatitude
    if( abs(FaceCoords_D(z_)) > zMin )then

       ! Outflow parallel with the magnetic field 
       uPerB0 = uOutDim*Io2No_V(UnitU_)/sqrt(sum(B0Face_D**2))

       ! Change sign for northward hemisphere so B0_D/B0 points outward
       if(FaceCoords_D(3) > 0) uPerB0 = -uPerB0

       VarsGhostFace_V(iRhoUx_I) = B0Face_D(x_)*uPerB0
       VarsGhostFace_V(iRhoUy_I) = B0Face_D(y_)*uPerB0
       VarsGhostFace_V(iRhoUz_I) = B0Face_D(z_)*uPerB0

       ! Apply body densities and pressures
       VarsGhostFace_V(iRho_I) = BodyRho_I
       VarsGhostFace_V(iP_I)   = BodyP_I

    else
       ! Apply ionosphere-like boundary conditions at low latitudes

       VarsGhostFace_V(Rho_)     = BodyRho_I(1)
       VarsGhostFace_V(OpRho_)   = VarsTrueFace_V(OpRho_)
       VarsGhostFace_V(iRhoUx_I) = -VarsTrueFace_V(iRhoUx_I)
       VarsGhostFace_V(iRhoUy_I) = -VarsTrueFace_V(iRhoUy_I)
       VarsGhostFace_V(iRhoUz_I) = -VarsTrueFace_V(iRhoUz_I)
       VarsGhostFace_V(iP_I)     = VarsTrueFace_V(iP_I)

    end if

    VarsGhostFace_V(Bx_:Bz_)  = VarsTrueFace_V(Bx_:Bz_)

  end subroutine user_face_bcs


end module ModUser
