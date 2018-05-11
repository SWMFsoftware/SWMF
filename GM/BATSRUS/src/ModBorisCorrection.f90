module ModBorisCorrection

  use ModVarIndexes, ONLY: nVar, Rho_, RhoUx_, RhoUz_, Bx_, By_, Bz_, &
       Ux_, Uy_, Uz_, IsMhd
  use ModConst, ONLY: cLightSpeed
  use ModPhysics, ONLY: c2Light, InvClight2, ClightFactor, Si2No_V, UnitU_
  use ModCoordTransform, ONLY: cross_product
  use ModMain,    ONLY: UseB0, UseHalfStep, nStage
  use ModB0,      ONLY: B0_DGB, B0_DX, B0_DY, B0_DZ
  use ModAdvance, ONLY: State_VGB, StateOld_VGB, Source_VC, &
       Energy_GBI, EnergyOld_CBI, &
       LeftState_VX, RightState_VX, &
       LeftState_VY, RightState_VY, LeftState_VZ, RightState_VZ
  use BATL_lib,   ONLY: CellVolume_GB, Used_GB, nI, nJ, nK, nDim, MaxDim, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nINode, nJNode, nKNode, &
       x_, y_, z_, &
       get_region_indexes, block_inside_regions, test_start, test_stop

  implicit none
  private ! except

  public:: init_mod_boris_correction  ! initialize module
  public:: clean_mod_boris_correction ! clean module
  public:: read_boris_param    ! read parameters for (simple) Boris correction
  public:: mhd_to_boris        ! from classical to (simple) Boris variables
  public:: boris_to_mhd        ! from (simple) Boris to classical variables
  public:: add_boris_source    ! add source term proportional to div(E) 
  public:: boris_to_mhd_x, boris_to_mhd_y, boris_to_mhd_z ! convert faces
  public:: set_clight_cell     ! calculate speed of light at cell centers
  public:: set_clight_face     ! calculate speed of light at cell faces

  logical, public:: UseBorisRegion = .false.
  logical, public:: UseBorisCorrection = .false.
  logical, public:: UseBorisSimple = .false.

  ! Electric field . area vector for div(E) in Boris correction
  real, allocatable, public:: &
       EDotFA_X(:,:,:), EDotFA_Y(:,:,:), EDotFA_Z(:,:,:)

  !$omp threadprivate(EDotFA_X, EDotFA_Y, EDotFA_Z)

  ! Speed of light at cell centers and cell faces
  real, allocatable, public:: Clight_G(:,:,:), Clight_DF(:,:,:,:)
  
  ! Local variables ---------------

  ! Description of the region where Boris correction is used
  character(len=200):: StringBorisRegion = 'none'

  ! Indexes of regions defined with the #REGION commands
  integer, allocatable:: iRegionBoris_I(:)

contains
  !============================================================================
  subroutine init_mod_boris_correction

    use ModMain, ONLY: &
         iMinFace, iMaxFace, jMinFace, jMaxFace, kMinFace, kMaxFace
    !--------------------------------------------------------------------------
    ! Get signed indexes for Boris region(s)
    call get_region_indexes(StringBorisRegion, iRegionBoris_I)
    UseBorisRegion = allocated(iRegionBoris_I)
    
    !$omp parallel
    if(.not.allocated(EDotFA_X))then
       allocate(EDotFA_X(nI+1,jMinFace:jMaxFace,kMinFace:kMaxFace))
       allocate(EDotFA_Y(iMinFace:iMaxFace,nJ+1,kMinFace:kMaxFace))
       allocate(EDotFA_Z(iMinFace:iMaxFace,jMinFace:jMaxFace,nK+1))
    end if
    !$omp end parallel

  end subroutine init_mod_boris_correction

  !============================================================================
  subroutine clean_mod_boris_correction

    if(allocated(EDotFA_X)) deallocate(EDotFA_X)
    if(allocated(EDotFA_Y)) deallocate(EDotFA_Y)
    if(allocated(EDotFA_Z)) deallocate(EDotFA_Z)

  end subroutine clean_mod_boris_correction

  !============================================================================
  subroutine read_boris_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModMultiFluid, ONLY: nIonFluid

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_boris_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#BORIS")
       call read_var('UseBorisCorrection', UseBorisCorrection)
       if(UseBorisCorrection) then
          call read_var('ClightFactor', CLightFactor)
          if(IsMhd .and. nIonFluid == 1)then
             UseBorisSimple = .false.
          else
             ! For non-MHD equations only simplified Boris correction
             ! is possible
             UseBorisSimple     = .true.
             UseBorisCorrection = .false.
          end if
       else
          ClightFactor = 1.0
       end if

    case("#BORISSIMPLE")
       call read_var('UseBorisSimple', UseBorisSimple)
       if(UseBorisSimple) then
          call read_var('ClightFactor', ClightFactor)
          UseBorisCorrection=.false.
       else
          ClightFactor = 1.0
       end if

    case("#BORISREGION")
       call read_var('StringBorisRegion', StringBorisRegion)

    case default
       call stop_mpi(NameSub//': unknown NameCommand='//NameCommand)
    end select

  end subroutine read_boris_param
  !===========================================================================
  subroutine mhd_to_boris(iBlock)

    ! Convert StateOld_VGB and State_VGB from classical MHD variables
    ! to (simple) Boris variable

    integer, intent(in):: iBlock
    integer:: i, j, k

    real:: InvClight2Cell
    
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'mhd_to_boris'
    !------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseBorisRegion)then
       call set_clight_cell(iBlock)
    else
       InvClight2Cell = InvClight2
    end if
    
    if(UseBorisCorrection)then
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE

          if(UseBorisRegion) InvClight2Cell = 1/Clight_G(i,j,k)**2
          
          call mhd_to_boris_cell(StateOld_VGB(:,i,j,k,iBlock), &
               EnergyOld_CBI(i,j,k,iBlock,1))
          ! State_VGB is not used in 1-stage and HalfStep schemes
          if(.not.UseHalfStep .and. nStage > 1) &
               call mhd_to_boris_cell(State_VGB(:,i,j,k,iBlock), &
               Energy_GBI(i,j,k,iBlock,1))

       end do; end do; end do
    elseif(UseBorisSimple)then
       ! Update using simplified Boris correction, i.e. update
       !
       !    RhoUBorisSimple = (1+B^2/(rho*c^2)) * RhoU
       !
       ! instead of RhoU. See Gombosi et al JCP 2002, 177, 176 (eq. 38-39)

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE

          if(UseBorisRegion) InvClight2Cell = 1/Clight_G(i,j,k)**2
          
          call mhd_to_boris_simple_cell(StateOld_VGB(:,i,j,k,iBlock))
          ! State_VGB is not used in 1-stage and HalfStep schemes
          if(.not.UseHalfStep .and. nStage > 1) &
               call mhd_to_boris_simple_cell(State_VGB(:,i,j,k,iBlock))

       end do; end do; end do
    end if

    call test_stop(NameSub, DoTest, iBlock)
    !==========================================================================
  contains

    subroutine mhd_to_boris_cell(State_V, Energy)

      real, intent(inout):: State_V(nVar)
      real, intent(inout):: Energy

      real:: Rho, b_D(3), u_D(3)
      !------------------------------------------------------------------------
      b_D = State_V(Bx_:Bz_)
      if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

      Rho = State_V(Rho_)
      u_D = State_V(RhoUx_:RhoUz_)/Rho

      ! Gombosi et al. 2001, eq(12) with vA^2 = B^2/Rho
      !
      ! RhoUBoris = RhoU + (RhoU B^2 - B RhoU.B)/(Rho c^2)
      !           = U*(Rho + B^2/c^2 - B U.B/c^2
      State_V(RhoUx_:RhoUz_) = u_D*(Rho + sum(b_D**2)*InvClight2Cell) &
           - b_D*sum(u_D*b_D)*InvClight2Cell

      ! e_Boris = e + (UxB)^2/(2 c^2)   eq 92
      Energy = Energy + 0.5*sum(cross_product(u_D, b_D)**2)*InvClight2Cell

    end subroutine mhd_to_boris_cell
    !==========================================================================
    subroutine mhd_to_boris_simple_cell(State_V)

      real, intent(inout):: State_V(nVar)

      ! Replace MHD momentum with (1+B^2/(rho*c^2)) * RhoU
      ! Use B0=B0_D in the total magnetic field if present.

      real:: b_D(3), Factor
      character(len=*), parameter:: NameSub = 'mhd_to_boris_simple'
      !------------------------------------------------------------------------
      b_D = State_V(Bx_:Bz_)
      if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

      ! Gombosi et al. 2001, eq(38) with vA^2 = B^2/Rho
      Factor = 1 + sum(b_D**2)*InvClight2Cell/State_V(Rho_)
      State_V(RhoUx_:RhoUz_) = Factor*State_V(RhoUx_:RhoUz_)

    end subroutine mhd_to_boris_simple_cell

  end subroutine mhd_to_boris

  !============================================================================
  subroutine boris_to_mhd(iBlock)

    integer, intent(in):: iBlock

    integer:: i, j, k

    real:: Clight2Cell, InvClight2Cell
    
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'boris_to_mhd'
    !-------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(UseBorisRegion)then
       call set_clight_cell(iBlock)
    else
       Clight2Cell    = C2light
       InvClight2Cell = InvClight2
    end if

    if(UseBorisCorrection)then
       ! Convert StateOld_VGB and State_VGB back from
       ! semi-relativistic to classical MHD variables
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE
          
          if(UseBorisRegion)then
             Clight2Cell = Clight_G(i,j,k)**2
             InvClight2Cell = 1/Clight2Cell
          end if

          call boris_to_mhd_cell(StateOld_VGB(:,i,j,k,iBlock), &
               EnergyOld_CBI(i,j,k,iBlock,1))
          call boris_to_mhd_cell(State_VGB(:,i,j,k,iBlock), &
               Energy_GBI(i,j,k,iBlock,1))

       end do; end do; end do
    elseif(UseBorisSimple)then
       ! Convert mometum in StateOld_VGB and State_VGB back from
       ! enhanced momentum
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.Used_GB(i,j,k,iBlock)) CYCLE
          
          if(UseBorisRegion) Clight2Cell = Clight_G(i,j,k)**2

          call boris_simple_to_mhd_cell(StateOld_VGB(:,i,j,k,iBlock))
          call boris_simple_to_mhd_cell(State_VGB(:,i,j,k,iBlock))

       end do; end do; end do
    end if

    call test_stop(NameSub, DoTest, iBlock)

  contains
    !==========================================================================
    subroutine boris_to_mhd_cell(State_V, Energy)

      real, intent(inout):: State_V(nVar)
      real, intent(inout):: Energy

      ! Replace semi-relativistic momentum with classical momentum in State_V.
      ! Replace semi-relativistic energy density Energy with classical value.
      ! Use B0=B0_D in the total magnetic field if present.

      real:: RhoC2, b_D(3), RhoUBoris_D(3), u_D(3)
      !------------------------------------------------------------------------
      b_D = State_V(Bx_:Bz_)
      if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

      RhoC2       = State_V(Rho_)*cLight2Cell
      RhoUBoris_D = State_V(RhoUx_:RhoUz_)

      ! Gombosi et al. 2001, eq(16) with vA^2 = B^2/Rho, gA^2=1/(1+vA^2/c^2)
      !
      ! RhoU = [RhoUBoris + B*B.RhoUBoris/(Rho c^2)]/(1+B^2/Rho c^2)
      !      = (RhoUBoris*Rho*c^2 + B*B.RhoUBoris)/(Rho*c^2 + B^2)

      State_V(RhoUx_:RhoUz_) =(RhoC2*RhoUBoris_D + b_D*sum(b_D*RhoUBoris_D)) &
           /(RhoC2 + sum(b_D**2))

      ! e = e_boris - (U x B)^2/(2 c^2)   eq 92
      u_D = State_V(RhoUx_:RhoUz_)/State_V(Rho_)
      Energy = Energy - 0.5*sum(cross_product(u_D, b_D)**2)*InvClight2Cell

    end subroutine boris_to_mhd_cell
    !=========================================================================
    subroutine boris_simple_to_mhd_cell(State_V)

      real, intent(inout)          :: State_V(nVar)

      ! Replace (1+B^2/(rho*c^2)) * RhoU with RhoU
      ! Use B0=B0_D in the total magnetic field if present.

      real:: b_D(3), Factor
      !------------------------------------------------------------------------
      b_D = State_V(Bx_:Bz_)
      if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

      ! Gombosi et al. 2001, eq(38) with vA^2 = B^2/Rho
      Factor = 1 + sum(b_D**2)/(State_V(Rho_)*cLight2Cell)
      State_V(RhoUx_:RhoUz_) = State_V(RhoUx_:RhoUz_)/Factor

    end subroutine boris_simple_to_mhd_cell

  end subroutine boris_to_mhd
  !============================================================================
  subroutine add_boris_source(iBlock)

    ! Add E div(E)*(1/c0^2 - 1/c^2) source term to the momentum equation
    ! See eq 28 in Gombosi et al. 2001 JCP, doi:10.1006/jcph.2002.7009

    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: Coef, InvClightOrig2
    real :: FullB_D(MaxDim), E_D(MaxDim), DivE

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'add_boris_source'
    !-------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    
    if(UseBorisRegion) call set_clight_cell(iBlock)
    InvClightOrig2 = ClightFactor**2*InvClight2
    Coef = InvClightOrig2 - InvClight2
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.Used_GB(i,j,k,iBlock)) CYCLE

       if(UseBorisRegion)then
          ! 1/c0^2 - 1/c'^2
          Coef = InvClightOrig2 - 1/Clight_G(i,j,k)**2
          if(abs(Coef) < 1e-20) CYCLE
       end if

       FullB_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       if(UseB0)FullB_D = FullB_D + B0_DGB(:,i,j,k,iBlock)

       E_D = cross_product(FullB_D,&
            State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock))/&
            State_VGB(Rho_,i,j,k,iBlock)

       ! Calculate divergence of electric field
       DivE =                     EDotFA_X(i+1,j,k) - EDotFA_X(i,j,k)
       if(nDim > 1) DivE = DivE + EDotFA_Y(i,j+1,k) - EDotFA_Y(i,j,k)
       if(nDim > 2) DivE = DivE + EDotFA_Z(i,j,k+1) - EDotFA_Z(i,j,k)
       DivE = DivE/CellVolume_GB(i,j,k,iBlock)

       Source_VC(RhoUx_:RhoUz_,i,j,k) = Source_VC(RhoUx_:RhoUz_,i,j,k) &
            + Coef*DivE*E_D

    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)

  end subroutine add_boris_source

  !==========================================================================
  subroutine boris_to_mhd_x(iMin,iMax,jMin,jMax,kMin,kMax)

    ! Convert face centered Boris momenta to MHD velocities

    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    integer:: i, j, k
    real:: InvClight2Face
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       if(UseBorisRegion)then
          InvClight2Face = 1/Clight_DF(1,i,j,k)**2
       else
          InvClight2Face = InvClight2
       end if

       ! Left face values
       RhoInv = 1/LeftState_VX(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DX(x_,i,j,k) + LeftState_VX(Bx_,i,j,k)
          ByFull = B0_DX(y_,i,j,k) + LeftState_VX(By_,i,j,k)
          BzFull = B0_DX(z_,i,j,k) + LeftState_VX(Bz_,i,j,k)
       else
          BxFull = LeftState_VX(Bx_,i,j,k)
          ByFull = LeftState_VX(By_,i,j,k)
          BzFull = LeftState_VX(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  = InvClight2Face*RhoInv
       LeftState_VX(Ux_,i,j,k)=LeftState_VX(Ux_,i,j,k)*RhoInv
       LeftState_VX(Uy_,i,j,k)=LeftState_VX(Uy_,i,j,k)*RhoInv
       LeftState_VX(Uz_,i,j,k)=LeftState_VX(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VX(Ux_,i,j,k)*BxFull + &
            LeftState_VX(Uy_,i,j,k)*ByFull + &
            LeftState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris= 1/(1 + B2Full*RhoC2Inv)

       LeftState_VX(Ux_,i,j,k) = &
            Ga2Boris * (LeftState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
       LeftState_VX(Uy_,i,j,k) = &
            Ga2Boris * (LeftState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
       LeftState_VX(Uz_,i,j,k) = &
            Ga2Boris * (LeftState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv = 1/RightState_VX(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DX(x_,i,j,k) + RightState_VX(Bx_,i,j,k)
          ByFull = B0_DX(y_,i,j,k) + RightState_VX(By_,i,j,k)
          BzFull = B0_DX(z_,i,j,k) + RightState_VX(Bz_,i,j,k)
       else
          BxFull = RightState_VX(Bx_,i,j,k)
          ByFull = RightState_VX(By_,i,j,k)
          BzFull = RightState_VX(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =InvClight2Face*RhoInv
       RightState_VX(Ux_,i,j,k)=RightState_VX(Ux_,i,j,k)*RhoInv
       RightState_VX(Uy_,i,j,k)=RightState_VX(Uy_,i,j,k)*RhoInv
       RightState_VX(Uz_,i,j,k)=RightState_VX(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VX(Ux_,i,j,k)*BxFull + &
            RightState_VX(Uy_,i,j,k)*ByFull + &
            RightState_VX(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

       RightState_VX(Ux_,i,j,k) = &
            Ga2Boris * (RightState_VX(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VX(Uy_,i,j,k) = &
            Ga2Boris * (RightState_VX(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VX(Uz_,i,j,k) = &
            Ga2Boris * (RightState_VX(Uz_,i,j,k)+uBC2Inv*BzFull)

    end do; end do; end do

  end subroutine boris_to_mhd_x
  !==========================================================================
  subroutine boris_to_mhd_y(iMin,iMax,jMin,jMax,kMin,kMax)
    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    integer:: i, j, k
    real:: InvClight2Face
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       if(UseBorisRegion)then
          InvClight2Face = 1/Clight_DF(2,i,j,k)**2
       else
          InvClight2Face = InvClight2
       end if

       ! Left face values
       RhoInv = 1/LeftState_VY(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DY(x_,i,j,k) + LeftState_VY(Bx_,i,j,k)
          ByFull = B0_DY(y_,i,j,k) + LeftState_VY(By_,i,j,k)
          BzFull = B0_DY(z_,i,j,k) + LeftState_VY(Bz_,i,j,k)
       else
          BxFull = LeftState_VY(Bx_,i,j,k)
          ByFull = LeftState_VY(By_,i,j,k)
          BzFull = LeftState_VY(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =InvClight2Face*RhoInv
       LeftState_VY(Ux_,i,j,k)=LeftState_VY(Ux_,i,j,k)*RhoInv
       LeftState_VY(Uy_,i,j,k)=LeftState_VY(Uy_,i,j,k)*RhoInv
       LeftState_VY(Uz_,i,j,k)=LeftState_VY(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VY(Ux_,i,j,k)*BxFull + &
            LeftState_VY(Uy_,i,j,k)*ByFull + &
            LeftState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

       LeftState_VY(Ux_,i,j,k) = &
            Ga2Boris * (LeftState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
       LeftState_VY(Uy_,i,j,k) = &
            Ga2Boris * (LeftState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
       LeftState_VY(Uz_,i,j,k) = &
            Ga2Boris * (LeftState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv = 1/RightState_VY(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DY(x_,i,j,k) + RightState_VY(Bx_,i,j,k)
          ByFull = B0_DY(y_,i,j,k) + RightState_VY(By_,i,j,k)
          BzFull = B0_DY(z_,i,j,k) + RightState_VY(Bz_,i,j,k)
       else
          BxFull = RightState_VY(Bx_,i,j,k)
          ByFull = RightState_VY(By_,i,j,k)
          BzFull = RightState_VY(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv=InvClight2Face*RhoInv
       RightState_VY(Ux_,i,j,k)=RightState_VY(Ux_,i,j,k)*RhoInv
       RightState_VY(Uy_,i,j,k)=RightState_VY(Uy_,i,j,k)*RhoInv
       RightState_VY(Uz_,i,j,k)=RightState_VY(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VY(Ux_,i,j,k)*BxFull + &
            RightState_VY(Uy_,i,j,k)*ByFull + &
            RightState_VY(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

       RightState_VY(Ux_,i,j,k) = &
            Ga2Boris * (RightState_VY(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VY(Uy_,i,j,k) = &
            Ga2Boris * (RightState_VY(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VY(Uz_,i,j,k) = &
            Ga2Boris * (RightState_VY(Uz_,i,j,k)+uBC2Inv*BzFull)
    end do; end do; end do

  end subroutine boris_to_mhd_y
  !==========================================================================
  subroutine boris_to_mhd_z(iMin,iMax,jMin,jMax,kMin,kMax)

    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    ! Convert face centered Boris momenta/rho to MHD velocities

    integer:: i, j, k
    real:: InvClight2Face
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

       if(UseBorisRegion)then
          InvClight2Face = 1/Clight_DF(3,i,j,k)**2
       else
          InvClight2Face = InvClight2
       end if

       ! Left face values
       RhoInv = 1/LeftState_VZ(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DZ(x_,i,j,k) + LeftState_VZ(Bx_,i,j,k)
          ByFull = B0_DZ(y_,i,j,k) + LeftState_VZ(By_,i,j,k)
          BzFull = B0_DZ(z_,i,j,k) + LeftState_VZ(Bz_,i,j,k)
       else
          BxFull = LeftState_VZ(Bx_,i,j,k)
          ByFull = LeftState_VZ(By_,i,j,k)
          BzFull = LeftState_VZ(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =InvClight2Face*RhoInv
       LeftState_VZ(Ux_,i,j,k)=LeftState_VZ(Ux_,i,j,k)*RhoInv
       LeftState_VZ(Uy_,i,j,k)=LeftState_VZ(Uy_,i,j,k)*RhoInv
       LeftState_VZ(Uz_,i,j,k)=LeftState_VZ(Uz_,i,j,k)*RhoInv
       uBC2Inv= (LeftState_VZ(Ux_,i,j,k)*BxFull + &
            LeftState_VZ(Uy_,i,j,k)*ByFull + &
            LeftState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

       LeftState_VZ(Ux_,i,j,k) = &
            Ga2Boris * (LeftState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
       LeftState_VZ(Uy_,i,j,k) = &
            Ga2Boris * (LeftState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
       LeftState_VZ(Uz_,i,j,k) = &
            Ga2Boris * (LeftState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)

       ! Right face values
       RhoInv = 1/RightState_VZ(rho_,i,j,k)
       if(UseB0)then
          BxFull = B0_DZ(x_,i,j,k) + RightState_VZ(Bx_,i,j,k)
          ByFull = B0_DZ(y_,i,j,k) + RightState_VZ(By_,i,j,k)
          BzFull = B0_DZ(z_,i,j,k) + RightState_VZ(Bz_,i,j,k)
       else
          BxFull = RightState_VZ(Bx_,i,j,k)
          ByFull = RightState_VZ(By_,i,j,k)
          BzFull = RightState_VZ(Bz_,i,j,k)
       end if
       B2Full = BxFull**2 + ByFull**2 + BzFull**2
       RhoC2Inv  =InvClight2Face*RhoInv
       RightState_VZ(Ux_,i,j,k)=RightState_VZ(Ux_,i,j,k)*RhoInv
       RightState_VZ(Uy_,i,j,k)=RightState_VZ(Uy_,i,j,k)*RhoInv
       RightState_VZ(Uz_,i,j,k)=RightState_VZ(Uz_,i,j,k)*RhoInv
       uBC2Inv= (RightState_VZ(Ux_,i,j,k)*BxFull + &
            RightState_VZ(Uy_,i,j,k)*ByFull + &
            RightState_VZ(Uz_,i,j,k)*BzFull)*RhoC2Inv

       ! gammaA^2 = 1/[1+BB/(rho c^2)]
       Ga2Boris = 1/(1 + B2Full*RhoC2Inv)

       RightState_VZ(Ux_,i,j,k) = &
            Ga2Boris * (RightState_VZ(Ux_,i,j,k)+uBC2Inv*BxFull)
       RightState_VZ(Uy_,i,j,k) = &
            Ga2Boris * (RightState_VZ(Uy_,i,j,k)+uBC2Inv*ByFull)
       RightState_VZ(Uz_,i,j,k) = &
            Ga2Boris * (RightState_VZ(Uz_,i,j,k)+uBC2Inv*BzFull)
    end do; end do; end do

  end subroutine boris_to_mhd_z
  !===========================================================================
  subroutine set_clight_cell(iBlock)

    ! Set Clight_G, the reduced speed of light at cell centers
    ! including ghost cells

    integer, intent(in):: iBlock
    !-------------------------------------------------------------------------
    if(.not.allocated(Clight_G)) &
         allocate(Clight_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

    if(UseBorisRegion)then
       ! Get the Boris factor into Clight_G
       call block_inside_regions(iRegionBoris_I, iBlock, size(Clight_G), &
            'ghost', Value_I=Clight_G)
       ! Convert Boris factor to speed of light
       Clight_G = cLightSpeed * Si2No_V(UnitU_) * &
            (1 - (1 - ClightFactor)*Clight_G)
    else
       Clight_G = cLightSpeed * Si2No_V(UnitU_)
    end if

  end subroutine set_clight_cell
  !===========================================================================
  subroutine set_clight_face(iBlock)

    ! Set Clight_DF, the reduced speed of light at cell faces (only nDim sides)

    integer, intent(in):: iBlock
    !-------------------------------------------------------------------------
    if(.not.allocated(Clight_DF)) &
         allocate(Clight_DF(nDim,nINode,nJNode,nKNode))

    if(UseBorisRegion)then
       ! Get the Boris Factor into Clight_G
       call block_inside_regions(iRegionBoris_I, iBlock, size(Clight_DF), &
            'face', Value_I=Clight_DF)
       ! Convert Boris factor to speed of light
       Clight_DF = cLightSpeed * Si2No_V(UnitU_) * &
            (1 - (1 - ClightFactor)*Clight_DF)
    else
       Clight_DF = cLightSpeed * Si2No_V(UnitU_)
    end if

  end subroutine set_clight_face
  !===========================================================================

end module ModBorisCorrection
!==============================================================================
