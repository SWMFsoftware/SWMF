module ModBorisCorrection

  use ModVarIndexes,     ONLY: nVar, Rho_, RhoUx_, RhoUz_, Bx_, By_, Bz_, &
       Ux_, Uy_, Uz_, IsMhd
  use ModPhysics,        ONLY: c2Light, InvClight2, ClightFactor
  use ModCoordTransform, ONLY: cross_product
  use ModMain,    ONLY: UseB0
  use ModB0,      ONLY: B0_DGB, B0_DX, B0_DY, B0_DZ
  use ModAdvance, ONLY: State_VGB, Source_VC, LeftState_VX, RightState_VX, &
       LeftState_VY, RightState_VY, LeftState_VZ, RightState_VZ
  use BATL_lib,   ONLY: CellVolume_GB, Used_GB, nI, nJ, nK, nDim, MaxDim, &
       x_, y_, z_, get_region_indexes

  implicit none
  private ! except

  public:: init_mod_boris_correction  ! initialize module
  public:: clean_mod_boris_correction ! clean module
  public:: read_boris_param    ! read parameters for (simple) Boris correction
  public:: mhd_to_boris        ! from classical to semi-relativistic variables
  public:: boris_to_mhd        ! from semi-relativistic to classical variables
  public:: mhd_to_boris_simple ! from RhoU to (1+B^2/(Rho*c^2))*RhoU
  public:: boris_simple_to_mhd ! from (1+B^2/(Rho*c^2))*RhoU to RhoU
  public:: add_boris_source    ! add source term proportional to div(E) 
  public:: boris_to_mhd_x, boris_to_mhd_y, boris_to_mhd_z ! convert faces

  logical, public:: UseBorisCorrection = .false.
  logical, public:: UseBorisSimple = .false.

  ! Electric field . area vector for div(E) in Boris correction
  real, allocatable, public:: &
       EDotFA_X(:,:,:), EDotFA_Y(:,:,:), EDotFA_Z(:,:,:)

  !$omp threadprivate(EDotFA_X, EDotFA_Y, EDotFA_Z)

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

  !============================================================================
  subroutine mhd_to_boris(State_V, Energy, B0_D)

    real, intent(inout)          :: State_V(nVar)
    real, intent(inout), optional:: Energy
    real, intent(in),    optional:: B0_D(3)

    ! Replace MHD momentum with semi-relativistic momentum in State_V.
    ! Replace MHD energy density Energy with semi-relativistic value.
    ! Use B0=B0_D in the total magnetic field if present.

    real:: Rho, b_D(3), u_D(3)
    character(len=*), parameter:: NameSub = 'mhd_to_boris'
    !--------------------------------------------------------------------------
    b_D = State_V(Bx_:Bz_)
    if(present(b0_D)) b_D = b_D + B0_D

    Rho = State_V(Rho_)
    u_D = State_V(RhoUx_:RhoUz_)/Rho

    ! Gombosi et al. 2001, eq(12) with vA^2 = B^2/Rho
    !
    ! RhoUBoris = RhoU + (RhoU B^2 - B RhoU.B)/(Rho c^2)
    !           = U*(Rho + B^2/c^2 - B U.B/c^2
    State_V(RhoUx_:RhoUz_) = u_D*(Rho + sum(b_D**2)*InvClight2) &
         - b_D*sum(u_D*b_D)*InvClight2

    ! e_Boris = e + (UxB)^2/(2 c^2)   eq 92
    if(present(Energy)) &
         Energy = Energy + 0.5*sum(cross_product(u_D, b_D)**2)*InvClight2

  end subroutine mhd_to_boris
  !============================================================================
  subroutine boris_to_mhd(State_V, Energy, B0_D)

    real, intent(inout)          :: State_V(nVar)
    real, intent(inout), optional:: Energy
    real, intent(in),    optional:: B0_D(3)

    ! Replace semi-relativistic momentum with classical momentum in State_V.
    ! Replace semi-relativistic energy density Energy with classical value.
    ! Use B0=B0_D in the total magnetic field if present.

    real:: RhoC2, b_D(3), RhoUBoris_D(3), u_D(3)
    character(len=*), parameter:: NameSub = 'boris_to_mhd'
    !--------------------------------------------------------------------------
    b_D = State_V(Bx_:Bz_)
    if(present(b0_D)) b_D = b_D + B0_D

    RhoC2       = State_V(Rho_)*c2Light
    RhoUBoris_D = State_V(RhoUx_:RhoUz_)

    ! Gombosi et al. 2001, eq(16) with vA^2 = B^2/Rho, gA^2=1/(1+vA^2/c^2)
    !
    ! RhoU = [RhoUBoris + B*B.RhoUBoris/(Rho c^2)]/(1+B^2/Rho c^2)
    !      = (RhoUBoris*Rho*c^2 + B*B.RhoUBoris)/(Rho*c^2 + B^2)

    State_V(RhoUx_:RhoUz_) =(RhoC2*RhoUBoris_D + b_D*sum(b_D*RhoUBoris_D)) &
         /(RhoC2 + sum(b_D**2))

    if(present(Energy))then
       ! e = e_boris - (U x B)^2/(2 c^2)   eq 92
       u_D = State_V(RhoUx_:RhoUz_)/State_V(Rho_)
       Energy = Energy - 0.5*sum(cross_product(u_D, b_D)**2)*InvClight2
    end if

  end subroutine boris_to_mhd
  !============================================================================
  subroutine mhd_to_boris_simple(State_V, B0_D)

    real, intent(inout)          :: State_V(nVar)
    real, intent(in),    optional:: B0_D(3)

    ! Replace MHD momentum with (1+B^2/(rho*c^2)) * RhoU
    ! Use B0=B0_D in the total magnetic field if present.

    real:: b_D(3), Factor
    character(len=*), parameter:: NameSub = 'mhd_to_boris_simple'
    !--------------------------------------------------------------------------
    b_D = State_V(Bx_:Bz_)
    if(present(b0_D)) b_D = b_D + B0_D

    ! Gombosi et al. 2001, eq(38) with vA^2 = B^2/Rho
    Factor = 1 + sum(b_D**2)/(State_V(Rho_)*c2light)
    State_V(RhoUx_:RhoUz_) = Factor*State_V(RhoUx_:RhoUz_)

  end subroutine mhd_to_boris_simple
  !============================================================================
  subroutine boris_simple_to_mhd(State_V, B0_D)

    real, intent(inout)          :: State_V(nVar)
    real, intent(in),    optional:: B0_D(3)

    ! Replace (1+B^2/(rho*c^2)) * RhoU with RhoU
    ! Use B0=B0_D in the total magnetic field if present.

    real:: b_D(3), Factor
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'boris_simple_to_mhd'
    !--------------------------------------------------------------------------
    b_D = State_V(Bx_:Bz_)
    if(present(b0_D)) b_D = b_D + B0_D

    ! Gombosi et al. 2001, eq(38) with vA^2 = B^2/Rho
    Factor = 1 + sum(b_D**2)/(State_V(Rho_)*c2Light)
    State_V(RhoUx_:RhoUz_) = State_V(RhoUx_:RhoUz_)/Factor

  end subroutine boris_simple_to_mhd
  !============================================================================
  subroutine add_boris_source(iBlock)

    ! Add E div(E)*(1/c0^2 - 1/c^2) source term to the momentum equation
    ! See eq 28 in Gombosi et al. 2001 JCP, doi:10.1006/jcph.2002.7009

    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: Coef
    real :: FullB_D(MaxDim), E_D(MaxDim), DivE

    character(len=*), parameter:: NameSub = 'add_boris_source'
    !-------------------------------------------------------------------------
    Coef = (ClightFactor**2 - 1.0)*InvClight2
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.Used_GB(i,j,k,iBlock)) CYCLE

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

  end subroutine add_boris_source

  !==========================================================================
  subroutine boris_to_mhd_x(iMin,iMax,jMin,jMax,kMin,kMax)

    ! Convert face centered Boris momenta to MHD velocities

    integer, intent(in) :: iMin,iMax,jMin,jMax,kMin,kMax

    integer:: i, j, k
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

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
       RhoC2Inv  = InvClight2*RhoInv
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
       RhoC2Inv  =InvClight2*RhoInv
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
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

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
       RhoC2Inv  =InvClight2*RhoInv
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
       RhoC2Inv=InvClight2*RhoInv
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
    real:: RhoInv, RhoC2Inv, BxFull, ByFull, BzFull, B2Full, uBC2Inv, Ga2Boris
    !------------------------------------------------------------------------
    ! U_Boris=rhoU_Boris/rho
    ! U = 1/[1+BB/(rho c^2)]* (U_Boris + (UBorisdotB/(rho c^2) * B)

    do k=kMin, kMax; do j=jMin, jMax; do i=iMin,iMax

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
       RhoC2Inv  =InvClight2*RhoInv
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
       RhoC2Inv  =InvClight2*RhoInv
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

end module ModBorisCorrection
!==============================================================================
