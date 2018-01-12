module ModBorisCorrection

  use ModVarIndexes,     ONLY: nVar, Rho_, RhoUx_, RhoUz_, Bx_, Bz_, IsMhd
  use ModPhysics,        ONLY: c2Light, InvClight2, ClightFactor
  use ModCoordTransform, ONLY: cross_product

  implicit none
  private ! except

  public:: init_mod_boris_correction  ! initialize module
  public:: clean_mod_boris_correction ! clean module
  public:: read_boris_param    ! read parameters for (simple) Boris correction
  public:: mhd_to_boris        ! from classical to semi-relativistic variables
  public:: boris_to_mhd        ! from semi-relativistic to classical variables
  public:: mhd_to_boris_simple ! from RhoU to (1+B^2/(Rho*c^2))*RhoU
  public:: boris_simple_to_mhd ! from (1+B^2/(Rho*c^2))*RhoU to RhoU

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
    use BATL_lib, ONLY: nI, nJ, nK, get_region_indexes
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

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_boris_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#BORIS")
       call read_var('UseBorisCorrection', UseBorisCorrection)
       if(UseBorisCorrection) then
          call read_var('ClightFactor', CLightFactor)
          if(IsMhd)then
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

end module ModBorisCorrection
!==============================================================================
