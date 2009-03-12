module ModMultiIon

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

  ! calculate analytic Jacobian for point-implicit scheme?
  logical:: IsAnalyticJacobian = .false. !!! should be always true eventually

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
    case("#MULTIION")
       call read_var('LowDensityRatio',    LowDensityRatio)
       call read_var('DoRestrictMultiIon', DoRestrictMultiIon)
       if(DoRestrictMultiIon)then
          call read_var('MachNumberMultiIon',    MachNumberMultiIon)
          call read_var('ParabolaWidthMultiIon', ParabolaWidthMultiIon)
       end if
!!! Temporary. Should be done analytically all the time eventually
       call read_var('IsAnalyticJacobian', IsAnalyticJacobian)
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

    use ModMain,    ONLY: nI, nJ, nK, x_, y_, z_
    use ModAdvance, ONLY: State_VGB, Source_VC
    use ModGeometry,ONLY: UseCovariant, dx_BLK, dy_BLK, dz_BLK, vInv_CB, &
         FaceAreaI_DFB, FaceAreaJ_DFB, FaceAreaK_DFB

    integer, intent(in) :: iBlock
    
    ! For multi-ion MHD the gradient of electron pressure appears in 
    ! all the individual ion momentum equations as -n_i/n_e * grad Pe

    real :: NumDens_I(nIonFluid), InvNumDens
    real :: vInv, GradXPe, GradYPe, GradZPe, GradPe_D(3)
    integer :: i, j, k

    character(len=*), parameter:: NameSub = 'multi_ion_source_expl'
    !----------------------------------------------------------------------
    do k=1,nK; do j=1,nJ; do i=1,nI

       NumDens_I  = State_VGB(iRhoIon_I,i,j,k,iBlock) / MassIon_I
       InvNumDens = 1.0/sum(NumDens_I)

       if(UseCovariant)then
          vInv = vInv_CB(i,j,k,iBlock)

          GradPe_D = vInv* &
               ( Pe_X(i+1,j,k)*FaceAreaI_DFB(:,i+1,j,k,iBlock) &
               - Pe_X(i  ,j,k)*FaceAreaI_DFB(:,i  ,j,k,iBlock) &
               + Pe_Y(i,j+1,k)*FaceAreaJ_DFB(:,i,j+1,k,iBlock) &
               - Pe_Y(i,j  ,k)*FaceAreaJ_DFB(:,i,j  ,k,iBlock) &
               + Pe_Z(i,j,k+1)*FaceAreaK_DFB(:,i,j,k+1,iBlock) &
               - Pe_Z(i,j,k  )*FaceAreaK_DFB(:,i,j,k  ,iBlock) )
          GradXPe = GradPe_D(x_)
          GradYPe = GradPe_D(y_)
          GradZPe = GradPe_D(z_)
       else
          ! Simple implementation for Cartesian 
          GradXPe = (Pe_X(i+1,j,k) - Pe_X(i,j,k))/dx_BLK(iBlock)
          GradYPe = (Pe_Y(i,j+1,k) - Pe_Y(i,j,k))/dy_BLK(iBlock)
          GradZPe = (Pe_Z(i,j,k+1) - Pe_Z(i,j,k))/dz_BLK(iBlock)
       end if

       Source_VC(iRhoUxIon_I,i,j,k) = Source_VC(iRhoUxIon_I,i,j,k) &
            - NumDens_I*InvNumDens*GradXPe
       Source_VC(iRhoUyIon_I,i,j,k) = Source_VC(iRhoUyIon_I,i,j,k) &
            - NumDens_I*InvNumDens*GradYPe
       Source_VC(iRhoUzIon_I,i,j,k) = Source_VC(iRhoUzIon_I,i,j,k) &
            - NumDens_I*InvNumDens*GradZPe

     end do; end do; end do

  end subroutine multi_ion_source_expl

  !===========================================================================

  subroutine multi_ion_source_impl

    ! Evaluate the explicit or implicit or both source terms.
    ! If there is no explicit source term, the subroutine user_expl_source 
    ! and the corresponding calls can be removed.

    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplPerturbed, IsPointImplMatrixSet, DsDu_VVC
    use ModMain,    ONLY: GlobalBlk, nI, nJ, nK, UseB0,&
                          UseBoris => boris_correction
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
    real    :: InvCharge, NumDens, InvNumDens, pAverage, State_V(nVar)
    real, dimension(3) :: FullB_D, uIon_D, uIon2_D, u_D, uPlus_D, uPlusHallU_D
    real, dimension(3) :: Current_D, Force_D
    real, dimension(nIonFluid) :: &
         NumDens_I, Rho_I, InvRho_I, Ux_I, Uy_I, Uz_I, Temp_I

    ! Alfven Lorentz factor for Boris correction
    real :: Ga2

    integer :: iBlock, i, j, k, iIon, jIon, iRhoUx, iRhoUz, iP, iEnergy
    real :: CoefBx, CoefBy, CoefBz, Coef, AverageTemp, TemperatureCoef, Heating
    real :: CollisionRate_II(nIonFluid, nIonFluid), CollisionRate

    ! Artificial friction
    real :: InvuCutOff2, InvTauCutOff

    character (len=*), parameter :: NameSub = 'multi_ion_source_impl'
    logical :: DoTest, DoTestMe, DoTestCell

    ! Variables for analytic Jacobian
    integer :: iDim, jDim, kDim, iUi, iUk
    real    :: SignedB, ForceCoeff, Coeff, CoefJacobian, Du2
    real    :: Du_D(3)
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

    ! Add source term n_s*(- u_+ - w_H + u_s )xB for multi-ions
    ! where u_+ is the number density weighted average ion velocity,
    ! and w_H = -J/(e n_e) is the Hall velocity. Here
    ! e is the electron charge and n_e is the electron number density.

    InvCharge = 1.0/ElectronCharge

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
       NumDens_I  = State_V(iRhoIon_I) / MassIon_I
       NumDens    = sum(NumDens_I)
       InvNumDens = 1.0/NumDens

       Temp_I     = State_V(iPIon_I)/NumDens_I
       AverageTemp= sum(State_V(iPIon_I))*InvNumDens

       if(AverageTemp <= 0.0)then
          write(*,*)'ERROR: AverageTemp =',AverageTemp
          write(*,*)'i,j,k,iBlock,iProc =',i,j,k,iBlock,iProc
          write(*,*)'x,y,z              =', &
               x_BLK(i,j,k,iBlock), y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)
          write(*,*)'iRhoIon_I          =',iRhoIon_I
          write(*,*)'RhoIon_I           =',State_V(iRhoIon_I)
          write(*,*)'MassIon_I          =',MassIon_I
          write(*,*)'NumDens_I          =',NumDens_I
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
       uPlus_D(x_) = InvNumDens* sum(NumDens_I*Ux_I)
       uPlus_D(y_) = InvNumDens* sum(NumDens_I*Uy_I)
       uPlus_D(z_) = InvNumDens* sum(NumDens_I*Uz_I)

       ! Add the Hall velocity -J/(e n)
       ! old version: call get_current(i,j,k,iBlock,Current_D)
       Current_D = vInv_CB(i,j,k,iBlock)*&
            ( bCrossArea_DX(:,i+1,j,k) - bCrossArea_DX(:,i,j,k) &
            + bCrossArea_DY(:,i,j+1,k) - bCrossArea_DY(:,i,j,k) &
            + bCrossArea_DZ(:,i,j,k+1) - bCrossArea_DZ(:,i,j,k))

       uPlusHallU_D = uPlus_D - InvNumDens*InvCharge*Current_D

       TemperatureCoef = 1.0/(AverageTemp*sqrt(AverageTemp))

       Ga2 = 1.0
       if(UseBoris)then                       !^CFG IF BORISCORR BEGIN
          ! Reduce the J x B like terms by a factor of 
          ! Ga2 = 1/(1+v_A^2/c^2) = rho/(rho+B^2/c^2)
          Ga2 = State_V(Rho_)/(State_V(Rho_) + InvClight2*sum(FullB_D**2)) 
       end if                                 !^CFG END BORISCORR


       if(TauCutOffDim > 0)then
          InvTauCutOff = 1.0/(Io2No_V(UnitT_)*TauCutOffDim)
          InvuCutOff2  = 1.0/(Io2No_V(UnitU_)*uCutOffDim)**2
       end if
       if(DoTestCell)then
          if(UseBoris)write(*,*) NameSub,'Ga2=',Ga2
          write(*,*) NameSub,' FullB_D  =', FullB_D
          write(*,*) NameSub,' Current_D=', Current_D
          write(*,*) NameSub,' uPlus_D  =', uPlus_D
          write(*,*) NameSub,' uPlusHall=', uPlusHallU_D
       end if

       ! Calculate the source term for all the ion fluids
       do iIon = 1, nIonFluid
          ! call select_fluid
          uIon_D = (/ Ux_I(iIon),  Uy_I(iIon), Uz_I(iIon) /)
          u_D    = uIon_D - uPlusHallU_D
          ForceCoeff = Ga2 * ElectronCharge*NumDens_I(iIon)
          Force_D    = ForceCoeff * cross_product(u_D, FullB_D) 

          if(DoTestCell)then
             write(*,*) NameSub,' iIon =', iIon
             write(*,*) NameSub,' uIon_D   =', uIon_D
             write(*,*) NameSub,' u_D      =', u_D
             write(*,*) NameSub,' Force_D  =', Force_D
          end if

          ! Set corresponding matrix element
          if (IsAnalyticJacobian) then
             do kDim = 1,nDim
                iUk = iUxIon_I(iIon) + kDim - 1
                do iDim = 1,nDim
                   if(kDim == iDim) CYCLE
                   jDim = 6 - kDim - iDim
                   SignedB = iLeviCivita_III(iDim, jDim, kDim)*FullB_D(jDim)

                   ! This Jacobian term occurs with respect with the same fluid
                   iUi = iUxIon_I(iIon) + iDim - 1
                   DsDu_VVC(iUk, iUi, i, j, k) = DsDu_VVC(iUk, iUi, i, j, k) & 
                        + ForceCoeff*SignedB*InvRho_I(iIon)

                   Coeff = ForceCoeff*SignedB*InvNumDens
                   ! This term is with respect to any fluid
                   do jIon = 1, nIonFluid
                      iUi = iUxIon_I(jIon) + iDim - 1
                      DsDu_VVC(iUk,iUi,i,j,k) = &
                           DsDu_VVC(iUk,iUi,i,j,k) - Coeff/MassIon_I(jIon)
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

    use ModEnergy, ONLY: calc_energy
    use ModAdvance, ONLY: State_VGB, Energy_GBI, &
         Rho_, p_, RhoUx_, RhoUy_, RhoUz_, UseElectronPressure, Pe_
    use ModPhysics, ONLY: ElectronTemperatureRatio, LowDensityRatio

    ! Resolve the update of total fluid vs. ion fluids:
    !   - take care of minor fluids with very small densities 
    !   - take care of conservation of total density and energy

    integer, intent(in) :: iBlock
    logical, intent(in) :: IsFinal  ! true for the final update

    integer :: i, j, k
    real    :: Rho, InvRho, p, IonSum, InvSum
    logical :: IsMultiIon
    !-----------------------------------------------------------------------

    do k=1,nK; do j=1,nJ; do i=1,nI

       ! Total density and pressure
       Rho    = State_VGB(Rho_,i,j,k,iBlock)
       InvRho = 1/Rho
       p      = State_VGB(p_,i,j,k,iBlock)

       ! Remove electron pressure from the total pressure
       if(UseElectronPressure)then
          p = p -  State_VGB(Pe_,i,j,k,iBlock)
       else
          p = p / (1 + ElectronTemperatureRatio)
       end if

       if(.not.IsFinal)then
          InvSum = 1/sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
          State_VGB(iRhoIon_I,i,j,k,iBlock) = &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*Rho*InvSum

          ! Distribute total pressure among ions
          IonSum = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          State_VGB(iPIon_I,i,j,k,iBlock) = &
               State_VGB(iPIon_I,i,j,k,iBlock)*p*InvSum

          CYCLE
       end if

       ! Check if we are in a region with multiple ions or not
       ! Note that IsMultiIon_CB is not necessarily allocated
       IsMultiIon = .true.
       if(DoRestrictMultiIon)IsMultiIon = IsMultiIon_CB(i,j,k,iBlock)

       if(.not. IsMultiIon)then
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
               *InvRho*State_VGB(RhoUz_,i,j,k,iBlock)/Rho

          ! Set ion temperatures to be equal with the total
          State_VGB(iPIon_I,i,j,k,iBlock) = p*InvRho * &
               State_VGB(iRhoIon_I,i,j,k,iBlock) &
               *MassIon_I(1)/MassIon_I
       else
          ! Distribute total density among ions
          InvSum = 1/sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
          State_VGB(iRhoIon_I,i,j,k,iBlock) = &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*Rho*InvSum

          ! Distribute total pressure among ions
          InvSum = 1/sum(State_VGB(iPIon_I,i,j,k,iBlock))
          State_VGB(iPIon_I,i,j,k,iBlock) = &
               State_VGB(iPIon_I,i,j,k,iBlock)*p*InvSum

          ! Overwrite total momentum and energy with sum of ions
          IonSum = sum(State_VGB(iRhoUxIon_I,i,j,k,iBlock))
          Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) &
               + 0.5/Rho*(IonSum**2 - State_VGB(RhoUx_,i,j,k,iBlock)**2)
          State_VGB(RhoUx_,i,j,k,iBlock) = IonSum

          IonSum = sum(State_VGB(iRhoUyIon_I,i,j,k,iBlock))
          Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) + &
               0.5/Rho*(IonSum**2 - State_VGB(RhoUy_,i,j,k,iBlock)**2)
          State_VGB(RhoUy_,i,j,k,iBlock) = IonSum

          IonSum = sum(State_VGB(iRhoUzIon_I,i,j,k,iBlock))
          Energy_GBI(i,j,k,iBlock,1) = Energy_GBI(i,j,k,iBlock,1) + &
               0.5/Rho*(IonSum**2 - State_VGB(RhoUz_,i,j,k,iBlock)**2)
          State_VGB(RhoUz_,i,j,k,iBlock) = IonSum

       end if

    end do; end do; end do

    ! Set ion energies: no attempt to conserve the energy per ion fluid
    call calc_energy(1, nI, 1, nJ, 1, nK, iBlock, IonFirst_, IonLast_)

  end subroutine multi_ion_update

  !==========================================================================

end module ModMultiIon
