module ModMultiIon

!!! "resistive terms" = ion-electron collisions to be added
!!! RZ geometry terms are missing

  ! Calculate source terms for multi-ion MHD. 
  ! For sake of numerical stability this should be done with 
  ! a point-implicit scheme. 
  ! Allow for extra point-implicit sources in ModUser

  use ModMultiFluid
  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: UseUserSource, &
       iTest, jTest, kTest, VarTest, BlkTest, ProcTest
  use ModUser, ONLY: user_calc_sources, user_init_point_implicit
  use ModSize, ONLY: nI, nJ, nK

  implicit none

  private

  ! Public methods and variables
  public:: multi_ion_set_parameters
  public:: multi_ion_set_restrict
  public:: multi_ion_source_expl
  public:: multi_ion_source_impl
  public:: multi_ion_init_point_impl
  public:: multi_ion_update

  logical, public              :: DoRestrictMultiIon = .false.

  ! does the cell have significant amount of multiple ion fluids?
  logical, public, allocatable :: IsMultiIon_CB(:,:,:,:)

  ! Local variables

  ! parameters for selecting the single ion region
  real :: MachNumberMultiIon = 0.0
  real :: ParabolaWidthMultiIon = 0.0

  ! electron pressure on the faces for grad Pe
  real, public:: &
       Pe_X(0:nI+2,0:nJ+1,0:nK+1), &
       Pe_Y(0:nI+1,0:nJ+2,0:nK+1), &
       Pe_Z(0:nI+1,0:nJ+1,0:nK+2)

  ! collision coefficient
  real :: CollisionCoefDim = -1.0
  real :: CollisionCoef = -1.0

  ! artificial friction parameters
  real    :: uCutOffDim = 0.0     ! cut-off velocity
  real    :: TauCutOffDim = -1.0  ! cut-off time scale
  integer :: nPowerCutOff = 0     ! cut-off exponent

  ! calculate analytic Jacobian for point-implicit scheme
  logical, parameter:: IsAnalyticJacobian = .true.

  ! how to reconcile ions with total fluid
  logical :: DoAddRho  = .false.
  logical :: DoAddRhoU = .true.

  ! Minimum pressure ratio for a minor fluid (so it remains positive)
  real:: LowPressureRatio = 1e-10

contains

  !===========================================================================
  subroutine multi_ion_set_parameters(NameCommand)

    use ModSize, ONLY: nI, nJ, nK, MaxBlock
    use ModPhysics, ONLY: LowDensityRatio
    use ModReadParam, ONLY: read_var
    use ModPointImplicit, ONLY: IsPointImplMatrixSet

    character(len=*), intent(in):: NameCommand
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#MHDIONS")
       call read_var('DoAddRho',  DoAddRho)
       call read_var('DoAddRhoU', DoAddRhoU)
    case("#MULTIION")
       call read_var('LowDensityRatio',    LowDensityRatio)
       call read_var('LowPressureRatio',   LowPressureRatio)
       call read_var('DoRestrictMultiIon', DoRestrictMultiIon)
       if(DoRestrictMultiIon)then
          call read_var('MachNumberMultiIon',    MachNumberMultiIon)
          call read_var('ParabolaWidthMultiIon', ParabolaWidthMultiIon)
       end if
       IsPointImplMatrixSet = IsAnalyticJacobian
       
    case("#COLLISION")
       call read_var('CollisionCoefDim', CollisionCoefDim)
       call read_var('TauCutOff', TauCutOffDim)
       if(TauCutOffDim > 0.0)then
          call read_var('uCutOffDim', uCutOffDim)
          call read_var('nPowerCutOff', nPowerCutOff)
       end if
    end select

    if(DoRestrictMultiIon .and. .not.allocated(IsMultiIon_CB)) &
         allocate(IsMultiIon_CB(nI,nJ,nK,MaxBlock))

  end subroutine multi_ion_set_parameters

  !===========================================================================

  subroutine multi_ion_set_restrict(iBlock)

    ! Identify regions where only one ion fluid is present.

    use ModSize,     ONLY: nI, nJ, nK, MaxBlock
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, p_
    use ModPhysics,  ONLY: g
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK

    integer, intent(in) :: iBlock

    real    :: Rho, p, RhoUx
    integer :: i, j, k
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub = 'multi_ion_set_restrict'
    !----------------------------------------------------------------------

    if (iBlock == BlkTest .and. iProc == ProcTest) then 
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false. ; DoTestMe = .false.
    end if
    
    do k=1,nK; do j=1,nJ; do i=1,nI
       ! Check if we are in the solar wind
       Rho   = State_VGB(Rho_,i,j,k,iBlock)
       p     = State_VGB(p_,i,j,k,iBlock)
       RhoUx = State_VGB(RhoUx_,i,j,k,iBlock)
          
       IsMultiIon_CB(i,j,k,iBlock) = .not. &
            (RhoUx < 0.0 .and. RhoUx**2 > MachNumberMultiIon**2*g*p*Rho &
            .and. y_BLK(i,j,k,iBlock)**2 + z_BLK(i,j,k,iBlock)**2 > &
            -ParabolaWidthMultiIon * x_BLK(i,j,k,iBlock))

       if(DoTestMe .and. i == iTest .and. j == jTest .and. k == kTest) then
          write(*,*) NameSub,'Rho, p, RhoUx =',Rho, p, RhoUx
          write(*,*) NameSub,'RhoUx**2, MachNumberMultiIon*g*p*Rho=', &
               RhoUx**2, MachNumberMultiIon*g*p*Rho
          write(*,*) NameSub,'y**2, z**2, -ParabolaWidthMultiIon*x=', &
               y_BLK(i,j,k,iBlock)**2 , z_BLK(i,j,k,iBlock)**2, &
                -ParabolaWidthMultiIon * x_BLK(i,j,k,iBlock)
          write(*,*) NameSub, ' IsMultiIon_CB=',  IsMultiIon_CB(i,j,k,iBlock) 
       end if

    end do; end do; end do
  
  end subroutine multi_ion_set_restrict

  !===========================================================================

  subroutine multi_ion_source_expl(iBlock)

    ! Add non-stiff source terms specific to multi-ion MHD
    !
    ! 1. d(rho u_s)/dt +=     (n_s/n_e)*grad p_e
    !    d(e_s)/dt     += u_s.(n_s/n_e)*grad p_e
    !
    !    where s is the index of the ion fluid,
    !    p_e is the electron pressure and
    !    n_e is the electron number density.
    !    The electron pressure may be solved for (UseElectronPressure is true)
    !    or can be a fixed fraction (ElectronTemperatureRatio) of the total 
    !    pressure.

    use ModMain,    ONLY: nDim, nI, nJ, nK, x_, y_, z_, &
         UseB0, UseBoris => boris_correction, UseBorisSimple
    use ModAdvance, ONLY: State_VGB, Source_VC, B0_DGB, &
         bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ, UseElectronPressure
    use ModPhysics, ONLY: InvClight2 => Inv_C2light, ElectronTemperatureRatio
    use ModGeometry,ONLY: UseCovariant, dx_BLK, dy_BLK, dz_BLK, vInv_CB, &
         FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB
    use ModCoordTransform, ONLY: cross_product

    integer, intent(in) :: iBlock

    ! For multi-ion MHD the gradient of electron pressure appears in 
    ! all the individual ion momentum equations as -n_i/n_e * grad Pe

    real :: State_V(nVar)
    real, dimension(nDim)     :: Current_D, FullB_D, Force_D
    real, dimension(nIonFluid):: ForceX_I, ForceY_I, ForceZ_I, ChargeDens_I
    real :: InvElectronDens
    real :: vInv, GradXPe, GradYPe, GradZPe

    ! Alfven Lorentz factor for Boris correction
    real :: Ga2

    integer :: i, j, k

    character(len=*), parameter:: NameSub = 'multi_ion_source_expl'
    logical :: DoTest, DoTestMe, DoTestCell
    !----------------------------------------------------------------------
    
    if(iProc == ProcTest .and. iBlock==BLkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false. ; DoTestMe = .false.
    end if

    if(DoTestMe)then
       write(*,*)NameSub,': initial Source_VC=', Source_VC(VarTest,iTest,jTest,kTest)
       write(*,*)NameSub,':         State_VGB=',State_VGB(:,iTest,jTest,kTest,iBlock)
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell = DoTestMe .and. iTest==i .and. jTest==j .and. kTest==k 
       
       vInv = vInv_CB(i,j,k,iBlock)

       State_V = State_VGB(:,i,j,k,iBlock)

       ChargeDens_I = ChargeIon_I * State_V(iRhoIon_I) / MassIon_I
       InvElectronDens = 1.0/sum(ChargeDens_I)
       
       if(InvElectronDens < 0.0) &
            call stop_mpi('negative electron denisty')

       ! Calculate Lorentz force = J x B
       Current_D = vInv* &
            ( bCrossArea_DX(:,i+1,j,k) - bCrossArea_DX(:,i,j,k) &
            + bCrossArea_DY(:,i,j+1,k) - bCrossArea_DY(:,i,j,k) &
            + bCrossArea_DZ(:,i,j,k+1) - bCrossArea_DZ(:,i,j,k))

       FullB_D = State_V(Bx_:Bz_)
       if(UseB0) FullB_D =  FullB_D + B0_DGB(:,i,j,k,iBlock)

       ! Lorentz force: J x B
       Force_D = cross_product(Current_D, FullB_D)
       
       if(DoTestCell)write(*,*)NameSub,':Force_D=', Force_D

       ! Subtract electron pressure gradient force
       if(UseMultiIon .and. &
            (UseElectronPressure .or. ElectronTemperatureRatio > 0.0))then

          if(UseCovariant)then
             ! grad Pe = (1/Volume)*Integral P_e dAreaVector over cell surface

             Force_D = Force_D - vInv* &
                  ( Pe_X(i+1,j,k)*FaceAreaI_DFB(:,i+1,j,k,iBlock) &
                  - Pe_X(i  ,j,k)*FaceAreaI_DFB(:,i  ,j,k,iBlock) &
                  + Pe_Y(i,j+1,k)*FaceAreaJ_DFB(:,i,j+1,k,iBlock) &
                  - Pe_Y(i,j  ,k)*FaceAreaJ_DFB(:,i,j  ,k,iBlock) &
                  + Pe_Z(i,j,k+1)*FaceAreaK_DFB(:,i,j,k+1,iBlock) &
                  - Pe_Z(i,j,k  )*FaceAreaK_DFB(:,i,j,k  ,iBlock) )
          else
             ! Gradient of Pe in Cartesian case
             Force_D(x_) = Force_D(x_) &
                  - (Pe_X(i+1,j,k) - Pe_X(i,j,k))/dx_BLK(iBlock)
             Force_D(y_) = Force_D(y_) &
                  - (Pe_Y(i,j+1,k) - Pe_Y(i,j,k))/dy_BLK(iBlock)
             Force_D(z_) = Force_D(z_) &
                  - (Pe_Z(i,j,k+1) - Pe_Z(i,j,k))/dz_BLK(iBlock)
          end if
          if(DoTestCell)write(*,*)NameSub,': after grad Pe, Force_D=', Force_D

       end if

       if(UseBoris .or. UseBorisSimple)then
          ! Simplified Boris correction
          ! (see the ASTRONUM 2009 proceedings paper by Toth et al.)
          ! Divide the number density by
          !
          ! 1 + V_s^2/c^2 = 1 + B^2/(c^2*n*M_s/q_s)
          !
          ! where we used V_s^2 = (q_s*M/M_s)*B^2/rho = B^2/(n*M_s/q_s)
          ! and q_s is the charge of species s in units of electron charge.

          Ga2 = sum(FullB_D**2)*InvClight2*InvElectronDens
          ChargeDens_I = ChargeDens_I/(1 + Ga2*ChargeIon_I/MassIon_I)
       end if

       ! Multiply by n_s/n_e for all ion fluids
       ForceX_I = ChargeDens_I*InvElectronDens*Force_D(x_)
       ForceY_I = ChargeDens_I*InvElectronDens*Force_D(y_)
       ForceZ_I = ChargeDens_I*InvElectronDens*Force_D(z_)

       if(DoTestCell)then
          write(*,*)Namesub,':ForceX_I, ForceY_I, ForceZ_I=', &
               ForceX_I, ForceY_I, ForceZ_I
          write(*,*)Namesub,': InvElectronDens, ChargeDens_I=', &
               InvElectronDens, ChargeDens_I
       end if

       ! Store ion momentum sources
       Source_VC(iRhoUxIon_I,i,j,k) = Source_VC(iRhoUxIon_I,i,j,k) + ForceX_I
       Source_VC(iRhoUyIon_I,i,j,k) = Source_VC(iRhoUyIon_I,i,j,k) + ForceY_I
       Source_VC(iRhoUzIon_I,i,j,k) = Source_VC(iRhoUzIon_I,i,j,k) + ForceZ_I

       ! Calculate ion energy sources = u_s.Force_s
       Source_VC(nVar+IonFirst_:nVar+IonLast_,i,j,k) = &
            Source_VC(nVar+IonFirst_:nVar+IonLast_,i,j,k) + &
            ( State_V(iRhoUxIon_I)*ForceX_I &
            + State_V(iRhoUyIon_I)*ForceY_I &
            + State_V(iRhoUzIon_I)*ForceZ_I &
            ) / State_V(iRhoIon_I)

    end do; end do; end do

    if(DoTestMe)write(*,*)NameSub,': final Source_VC=',&
         Source_VC(VarTest,iTest,jTest,kTest)

  end subroutine multi_ion_source_expl

  !===========================================================================

  subroutine multi_ion_source_impl

    ! Add 'stiff' source terms specific to multi-ion MHD:
    !
    ! 1. d(rho u_s)/dt +=      n_s*(- u_+ - w_H + u_s )xB 
    !    d(e_s)/dt     += u_s.[n_s*(- u_+ - w_H + u_s )xB]
    !    where s is the index for the ion fluid, 
    !    u_+ is the charge density weighted average ion velocity,
    !    w_H = -J/(e n_e) is the Hall velocity, 
    !    n_e is the electron number density, and
    !    e is the electron charge.
    !
    ! 2. ion-ion collisions if required.
    !
    ! 3. artificial friction term if required.
    !
    ! 4. user source terms if required.

    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplPerturbed, IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, UseB0,&
                          UseBoris => boris_correction, UseBorisSimple
    use ModAdvance, ONLY: State_VGB, Source_VC, B0_DGB
    use ModAdvance, ONLY: bCrossArea_DX, bCrossArea_DY, bCrossArea_DZ
    use ModGeometry,ONLY: vInv_CB, x_BLK, y_BLK, z_BLK
    use ModPhysics, ONLY: ElectronCharge, gm1, inv_gm1, &
         InvClight2 => Inv_C2light, Si2No_V, No2Si_V, Io2No_V, &
         UnitTemperature_, UnitT_, UnitU_
         
    use ModMain,    ONLY: x_, y_, z_
    use ModCoordTransform, ONLY: cross_product
    use ModNumConst,       ONLY: iLeviCivita_III
    use ModSize,           ONLY: nDim

    ! Variables for multi-ion MHD
    real:: InvElectronDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, uIon2_D, u_D, uPlus_D
    real, dimension(3) :: Force_D
    real, dimension(nIonFluid) :: &
         NumDens_I, ChargeDens_I, ChargeDensBoris_I, &
         Rho_I, InvRho_I, Ux_I, Uy_I, Uz_I, Temp_I

    ! Alfven Lorentz factor for Boris correction
    real :: Ga2

    integer :: iBlock, i, j, k, iIon, jIon, iRhoUx, iRhoUz, iP, iEnergy
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, TemperatureCoef, Heating
    real :: CollisionRate_II(nIonFluid, nIonFluid), CollisionRate

    ! Artificial friction
    real :: InvuCutOff2, InvTauCutOff

    logical :: DoTest, DoTestMe, DoTestCell

    ! Variables for analytic Jacobian
    integer :: iDim, jDim, kDim, iUi, iUk
    real    :: SignedB, ForceCoeff, Coeff, CoefJacobian, Du2
    real    :: Du_D(3)

    character(len=*), parameter :: NameSub = 'multi_ion_source_impl'
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if
    DoTestCell = .false.

    ! Add user defined point implicit source terms here
    ! Explicit user sources are added in calc_sources
    if(UsePointImplicit .and. UseUserSource) call user_calc_sources

    ! Do not evaluate multi-ion sources in the numerical Jacobian calculation
    ! (needed for the user source terms) 
    if(IsPointImplPerturbed .and. IsAnalyticJacobian) RETURN

    if(CollisionCoefDim > 0.0)then
       ! Rate = n*CoefDim / T^1.5 with T [K], n [/cc] and Rate [1/s]
       CollisionCoef = CollisionCoefDim &
            /No2Si_V(UnitTemperature_)**1.5/Si2No_V(UnitT_)

       do jIon = 1, nIonFluid
          do iIon = 1, nIonFluid
             CollisionRate_II(iIon, jIon) = CollisionCoef* &
                  MassIon_I(iIon)*MassIon_I(jIon) &
                  /(MassIon_I(iIon)+MassIon_I(jIon))
          end do
       end do
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       if(DoTestCell)write(*,*)NameSub, ' initial source = ',&
            Source_VC(VarTest,i,j,k)

       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,iBlock)

       ! Total magnetic field
       FullB_D = State_V(Bx_:Bz_) 
       if(UseB0) FullB_D =  FullB_D + B0_DGB(:,i,j,k,iBlock)

       ! calculate number densities
       NumDens_I     = State_V(iRhoIon_I) / MassIon_I
       ChargeDens_I  = NumDens_I * ChargeIon_I
       InvElectronDens = 1.0/sum(ChargeDens_I)

       if(InvElectronDens < 0.0) &
            call stop_mpi(NameSub//': negative electron denisty')

       ChargeDensBoris_I = ChargeDens_I
       
       if(UseBoris .or. UseBorisSimple)then
          ! See the ASTRONUM 2009 proceedings paper by Toth et al.
          !
          ! Boris correction: divide the number density by
          !
          ! 1 + V_s^2/c^2 = 1 + B^2/(c^2*n*M_s/q_s)
          !
          ! where we used V_s^2 = (q_s*M/M_s)*B^2/rho = B^2/(n*M_s/q_s)

          Ga2 = sum(FullB_D**2)*InvClight2*InvElectronDens
          ChargeDensBoris_I = ChargeDens_I/(1 + Ga2*ChargeIon_I/MassIon_I)
       end if

       Temp_I     = State_V(iPIon_I)/NumDens_I
       AverageTemp= sum(State_V(iPIon_I))/sum(NumDens_I)

       if(AverageTemp <= 0.0)then
          write(*,*)'ERROR: AverageTemp =',AverageTemp
          write(*,*)'i,j,k,iBlock,iProc =',i,j,k,iBlock,iProc
          write(*,*)'x,y,z              =', &
               x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
          write(*,*)'iRhoIon_I          =',iRhoIon_I
          write(*,*)'RhoIon_I           =',State_V(iRhoIon_I)
          write(*,*)'MassIon_I          =',MassIon_I
          write(*,*)'ChargeIon_I        =',ChargeIon_I
          write(*,*)'NumDens_I          =',NumDens_I
          write(*,*)'ChargeDens_I       =',ChargeDens_I
          write(*,*)'iPIon_I            =',iPIon_I
          write(*,*)'PIon_I             =',State_V(iPIon_I)
          write(*,*)'Temp_I             =',Temp_I
          call stop_mpi(NameSub//': non-positive average temperature')
       end if
       
       Rho_I    = State_V(iRhoIon_I)
       InvRho_I = 1.0/Rho_I
       Ux_I     = InvRho_I*State_V(iUxIon_I)
       Uy_I     = InvRho_I*State_V(iUyIon_I)
       Uz_I     = InvRho_I*State_V(iUzIon_I)

       ! calculate the average positive charge velocity
       uPlus_D(x_) = InvElectronDens * sum(ChargeDens_I*Ux_I)
       uPlus_D(y_) = InvElectronDens * sum(ChargeDens_I*Uy_I)
       uPlus_D(z_) = InvElectronDens * sum(ChargeDens_I*Uz_I)

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       if(TauCutOffDim > 0)then
          InvTauCutOff = 1.0/(Io2No_V(UnitT_)*TauCutOffDim)
          InvuCutOff2  = 1.0/(Io2No_V(UnitU_)*uCutOffDim)**2
       end if
       if(DoTestCell)then
          if(UseBoris)write(*,*) NameSub,'Ga2=',Ga2
          write(*,*) NameSub,' FullB_D  =', FullB_D
          write(*,*) NameSub,' uPlus_D  =', uPlus_D
       end if

       ! Calculate the source term for all the ion fluids
       do iIon = 1, nIonFluid
          uIon_D = (/ Ux_I(iIon),  Uy_I(iIon), Uz_I(iIon) /)
          u_D    = uIon_D - uPlus_D
          ForceCoeff = ElectronCharge*ChargeDensBoris_I(iIon)
          Force_D    = ForceCoeff * cross_product(u_D, FullB_D) 

          if(DoTestCell)then
             write(*,*) NameSub,' iIon =', iIon
             write(*,*) NameSub,' uIon_D   =', uIon_D
             write(*,*) NameSub,' u_D      =', u_D
             write(*,*) NameSub,' Force_D  =', Force_D
          end if

          ! Set corresponding matrix element
          if (IsAnalyticJacobian .and. UsePointImplicit) then
             do kDim = 1,nDim
                iUk = iUxIon_I(iIon) + kDim - 1
                do iDim = 1,nDim
                   if(kDim == iDim) CYCLE
                   jDim = 6 - kDim - iDim
                   SignedB = iLeviCivita_III(iDim, jDim, kDim)*FullB_D(jDim)

                   ! This Jacobian term occurs with respect to the same fluid
                   iUi = iUxIon_I(iIon) + iDim - 1
                   DsDu_VVC(iUk, iUi, i, j, k) = DsDu_VVC(iUk, iUi, i, j, k) & 
                        + ForceCoeff*SignedB*InvRho_I(iIon)

                   Coeff = ForceCoeff*SignedB*InvElectronDens
                   ! This term is with respect to any fluid
                   do jIon = 1, nIonFluid
                      iUi = iUxIon_I(jIon) + iDim - 1
                      DsDu_VVC(iUk,iUi,i,j,k) = &
                           DsDu_VVC(iUk,iUi,i,j,k) - &
                           Coeff*ChargeIon_I(jIon)/MassIon_I(jIon)
                   end do
                end do
             end do
          end if
          Heating = 0.0

          if(CollisionCoefDim > 0.0 .or. TauCutOffDim > 0.0)then
             do jIon = 1, nIonFluid
                if(jIon == iIon) CYCLE
                
                ! Add collisional terms
                uIon2_D = (/ Ux_I(jIon),  Uy_I(jIon), Uz_I(jIon) /)

                ! Physical collision
                if(CollisionCoefDim > 0.0)then
                   CollisionRate = CollisionRate_II(iIon, jIon) * &
                        NumDens_I(iIon) * NumDens_I(jIon) &
                        * TemperatureCoef
                else
                   CollisionRate = 0.0
                end if

                ! Artificial friction to keep the velocity difference in check
                ! We take the smaller of the two densities so that the 
                ! acceleration is independent of the density of 
                ! the minor species and the restriction works in all regions.
                ! The min function is symmetric, so momentum is conserved.
                ! u_0 is the cut-off velocity, Tau gives the time rate, and
                ! the power determines how sharp the cut-off is.
                if(TauCutOffDim > 0.0)then
                   ! CollisionRate = 
                   !  1/tau * min(rho^iIon, rho^jIon) * (du2/u_0^2)^n

                   Du2 = sum( (uIon2_D - uIon_D)**2 )
                   CollisionRate = CollisionRate + &
                        InvTauCutOff * min(Rho_I(iIon), Rho_I(jIon)) &
                        * ( InvUCutOff2 * Du2 ) ** nPowerCutOff
                end if

                Force_D = Force_D + CollisionRate * (uIon2_D - uIon_D)

!!! No heating for now
                !Heating = Heating + CollisionRate* &
                !     ( 2*(Temp_I(jIon) - Temp_I(iIon)) &
                !     + gm1*sum((uIon2_D - uIon_D)**2) )
             
                ! Calculate corresponding matrix elements
                if (TauCutOffDim > 0.0 .and. IsAnalyticJacobian) then

                   ! du = u^iIon - u^jIon
                   Du_D = uIon_D - uIon2_D

                   ! Common coefficient: CoefJacobian = 
                   !  1/tau * min(rho^iIon, rho^jIon) * (1/u_0)^2n * (du^2)^n-1
                   CoefJacobian = InvTauCutOff &
                        * min(Rho_I(iIon), Rho_I(jIon)) &
                        * InvUCutOff2 ** nPowerCutOff &
                        * Du2 ** (nPowerCutOff - 1)

                   ! Add dFriction/d(RhoU) elements to the Jacobian
                   do kDim = 1,nDim
                      ! k component of RhoU^iIon
                      iUk = iUxIon_I(iIon) + kDim - 1
                      do iDim = 1,nDim

                         ! dFriction^iIon_k/d(RhoU^iIon_i) = -CoefJacobian
                         !  *(2*n*du_i*du_k/rho^iIon + delta_ik*du^2/rho^iIon)

                         iUi = iUxIon_I(iIon) + iDim - 1
                         DsDu_VVC(iUk, iUi, i, j, k) = &
                              DsDu_VVC(iUk, iUi, i, j, k) &
                              - 2.0 * nPowerCutOff * InvRho_I(iIon) &
                              * Du_D(iDim) * Du_D(kDim) & 
                              *  CoefJacobian
                         if (iDim == kDim) DsDu_VVC(iUk, iUi, i, j, k) = &
                              DsDu_VVC(iUk, iUi, i, j, k) & 
                              - CoefJacobian * Du2 *InvRho_I(iIon)

                         ! dFriction^iIon_k/d(RhoU^jIon_i) = +CoefJacobian
                         !  *(2*n*du_i*du_k/rho^jIon + delta_ik*du^2/rho^jIon)

                         iUi = iUxIon_I(jIon) + iDim - 1
                         DsDu_VVC(iUk, iUi, i, j, k) = &
                              DsDu_VVC(iUk, iUi, i, j, k) &
                              + 2.0 * nPowerCutOff *InvRho_I(jIon) &
                              * Du_D(iDim) * Du_D(kDim) &
                              * CoefJacobian
                         if (iDim == kDim)DsDu_VVC(iUk, iUi, i, j, k)  = &
                              DsDu_VVC(iUk, iUi, i, j, k) & 
                              + CoefJacobian * Du2 *InvRho_I(jIon)
                      end do
                   end do
                end if
             end do

             iP = iPIon_I(iIon)
             Source_VC(iP,i,j,k) = Source_VC(iP,i,j,k) + Heating

          end if

          iRhoUx = iRhoUxIon_I(iIon); iRhoUz = iRhoUzIon_I(iIon)
          Source_VC(iRhoUx:iRhoUz,i,j,k) = Source_VC(iRhoUx:iRhoUz,i,j,k) &
               + Force_D

          iEnergy = Energy_-2+iIon+IonFirst_
          Source_VC(iEnergy,i,j,k) = Source_VC(iEnergy,i,j,k) &
               + sum(Force_D*uIon_D) &
               + inv_gm1*Heating

       end do

       if(DoTestCell)write(*,*)NameSub, ' final source = ',&
            Source_VC(VarTest,i,j,k)

    end do; end do; end do

    if(DoTestMe)then
       write(*,*)NameSub,' CollisionCoef=',CollisionCoef
       write(*,*)NameSub,' AverageTemp  =',AverageTemp
       write(*,*)NameSub,' AverageTempDim=', &
            AverageTemp*No2Si_V(UnitTemperature_)
    end if

  end subroutine multi_ion_source_impl
  !===========================================================================
  subroutine multi_ion_init_point_impl

    ! Select variables for point implicit evaluation. This is the union
    ! of the ion momenta and the variables selected (if any) in 
    ! ModUser::user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet

    logical :: IsPointImpl_V(nVar)
    integer :: iVar, iPointImplVar, nPointImplVar,n
    !------------------------------------------------------------------------

    IsPointImpl_V = .false.
    IsPointImplMatrixSet = IsAnalyticJacobian

    if(UseUserSource)then
       call user_init_point_implicit
       if(allocated(iVarPointImpl_I)) then
          IsPointImpl_V(iVarPointImpl_I) = .true.
          deallocate(iVarPointImpl_I)
       end if
    end if

    ! All ion momenta and pressures are implicit
    IsPointImpl_V(iRhoUxIon_I) = .true.
    IsPointImpl_V(iRhoUyIon_I) = .true.
    IsPointImpl_V(iRhoUzIon_I) = .true.
    ! IsPointImpl_V(iPIon_I)   = .true. !!! No heating in artificial friction

    nPointImplVar = count(IsPointImpl_V)

    allocate(iVarPointImpl_I(nPointImplVar))

    iPointImplVar = 0
    do iVar = 1, nVar
       if(.not. IsPointImpl_V(iVar)) CYCLE
       iPointImplVar = iPointImplVar + 1
       iVarPointImpl_I(iPointImplVar) = iVar
    end do

  end subroutine multi_ion_init_point_impl

  !==========================================================================

  subroutine multi_ion_update(iBlock, IsFinal)

    ! Resolve the update of total fluid vs. ion fluids:
    !   - take care of minor fluids with very small densities 
    !   - take care of conservation of total density and energy

    use ModEnergy, ONLY: calc_energy
    use ModAdvance, ONLY: State_VGB, Energy_GBI, &
         Rho_, p_, RhoUx_, RhoUy_, RhoUz_, UseElectronPressure
    use ModPhysics, ONLY: ElectronTemperatureRatio, LowDensityRatio

    integer, intent(in) :: iBlock
    logical, intent(in) :: IsFinal  ! true for the final update

    integer :: i, j, k
    real    :: State_V(nVar), Rho, InvRho, p, IonSum, InvSum
    real    :: TeRatio1, InvTeRatio1
    logical :: IsMultiIon
    !-----------------------------------------------------------------------

    TeRatio1    = 1 + ElectronTemperatureRatio
    InvTeRatio1 = 1 / TeRatio1

    do k=1,nK; do j=1,nJ; do i=1,nI

       State_V = State_VGB(:,i,j,k,iBlock)

       ! Total density 
       Rho    = State_V(Rho_)
       InvRho = 1/Rho

       ! Total pressure
       p      = sum(State_V(iPIon_I))

       ! Keep pressures above LowPressureRatio*pTotal
       State_VGB(iPIon_I,i,j,k,iBlock) = &
            max( State_V(iPIon_I), LowPressureRatio*p )

       if(.not.IsFinal)then

          ! Ion pressures are always added up into total pressure
          if(UseElectronPressure) then
             State_VGB(p_,i,j,k,iBlock) = p
          else
             State_VGB(p_,i,j,k,iBlock) = p * TeRatio1
          end if

          if(DoAddRho)then
             State_VGB(Rho_,i,j,k,iBlock) = sum(State_V(iRhoIon_I))
          else
             ! Distribute total density among ions
             InvSum = 1/sum(State_V(iRhoIon_I))
             State_VGB(iRhoIon_I,i,j,k,iBlock) = State_V(iRhoIon_I)*Rho*InvSum
          end if

          IonSum = sum(State_V(iRhoUxIon_I))
          if(DoAddRhoU .or. &
               maxval(State_V(iRhoUxIon_I))* &
               minval(State_V(iRhoUxIon_I)) <= 0)then
             State_VGB(RhoUx_,i,j,k,iBlock) = IonSum
          else
             State_VGB(iRhoUxIon_I,i,j,k,iBlock) = State_V(iRhoUxIon_I) &
                  *State_V(RhoUx_)/IonSum
          endif

          IonSum = sum(State_V(iRhoUyIon_I))
          if(DoAddRhoU .or. &
               maxval(State_V(iRhoUyIon_I)) &
               *minval(State_V(iRhoUyIon_I)) <= 0)then
             State_VGB(RhoUy_,i,j,k,iBlock) = IonSum
          else
             State_VGB(iRhoUyIon_I,i,j,k,iBlock) = State_V(iRhoUyIon_I) &
                  *State_V(RhoUy_)/IonSum
          endif

          IonSum = sum(State_V(iRhoUzIon_I))
          if(DoAddRhoU .or. &
               maxval(State_V(iRhoUzIon_I)) &
               *minval(State_V(iRhoUzIon_I)) <= 0)then
             State_VGB(RhoUz_,i,j,k,iBlock) = IonSum
          else
             InvSum = 1/IonSum
             State_VGB(iRhoUzIon_I,i,j,k,iBlock) = State_V(iRhoUzIon_I) &
                  *State_V(RhoUz_)/IonSum
          endif

          ! Nothing more to do if not final update
          CYCLE
       end if

       ! Check if we are in a region with multiple ions or not
       ! Note that IsMultiIon_CB is not necessarily allocated
       IsMultiIon = .true.
       if(DoRestrictMultiIon)IsMultiIon = IsMultiIon_CB(i,j,k,iBlock)

       if(IsMultiIon)then
          ! Add up ion fluids to total fluid
          ! This is necessary when point-implicit source terms are evaluated
          ! for the ion fluids only and their sum is not 0 due to user sources

          State_VGB(Rho_,  i,j,k,iBlock) = sum(State_V(iRhoIon_I))
          State_VGB(RhoUx_,i,j,k,iBlock) = sum(State_V(iRhoUxIon_I))
          State_VGB(RhoUy_,i,j,k,iBlock) = sum(State_V(iRhoUyIon_I))
          State_VGB(RhoUz_,i,j,k,iBlock) = sum(State_V(iRhoUzIon_I))
          if(UseElectronPressure)then
             State_VGB(p_,i,j,k,iBlock)  = p
          else
             State_VGB(p_,i,j,k,iBlock)  = p * TeRatio1
          end if
       else
          ! Put most of the stuff into the first ion fluid
          State_VGB(iRhoIon_I(1),i,j,k,iBlock) = &
               Rho*(1.0 - LowDensityRatio*(IonLast_ - IonFirst_))
          State_VGB(iRhoIon_I(2:nIonFluid),i,j,k,iBlock) = &
               Rho*LowDensityRatio

          ! Set ion velocities to be equal with the total
          State_VGB(iRhoUxIon_I,i,j,k,iBlock) = &
               State_VGB(iRhoIon_I,i,j,k,iBlock) &
               *InvRho*State_VGB(RhoUx_,i,j,k,iBlock)
          State_VGB(iRhoUyIon_I,i,j,k,iBlock) = &
               State_VGB(iRhoIon_I,i,j,k,iBlock) &
               *InvRho*State_VGB(RhoUy_,i,j,k,iBlock)
          State_VGB(iRhoUzIon_I,i,j,k,iBlock) = &
               State_VGB(iRhoIon_I,i,j,k,iBlock) &
               *InvRho*State_VGB(RhoUz_,i,j,k,iBlock)

          ! Set ion temperatures to be equal with the total
          if(UseElectronPressure)then
             p = State_V(p_)
          else
             p = State_V(p_) * InvTeRatio1
          end if
          State_VGB(iPIon_I,i,j,k,iBlock) = p*InvRho * &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*MassIon_I(1)/MassIon_I
       end if

    end do; end do; end do

    ! Reset total and ion energies
    call calc_energy(1, nI, 1, nJ, 1, nK, iBlock, 1, IonLast_)

  end subroutine multi_ion_update

  !==========================================================================

end module ModMultiIon
