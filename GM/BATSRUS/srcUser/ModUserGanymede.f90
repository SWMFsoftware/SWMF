!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop
  ! User module for Ganymede. Inherited from ModUserMercury.f90
  ! Must compile with MHDHypPe equation set.
  use ModUserEmpty,                          &
       IMPLEMENTED1 => user_init_session,    &
       IMPLEMENTED2 => user_set_ics,         &
       IMPLEMENTED3 => user_set_resistivity, &
       IMPLEMENTED4 => user_read_inputs,     &
       IMPLEMENTED5 => user_set_cell_boundary,&
       IMPLEMENTED6 => user_set_face_boundary

  include 'user_module.h' ! list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = 'Ganymede, Hongyang Zhou'

  real :: PlanetDensity=-1., PlanetPressure=-1., PlanetRadius=-1.

  real :: PlanetDensitySi=-1., PlanetPressureSi=-1., PlanetRadiusSi=-1.

  logical :: UsePlanetDensity=.false., UsePlanetPressure=.false.

  integer :: nLayer =0 ! Number of points in planet resistivity profile
  real, allocatable :: PlanetRadiusSi_I(:),PlanetRadius_I(:),&
       ResistivitySi_I(:),Resistivity_I(:)
  real, allocatable :: ResistivityRate(:)

  integer :: iLayer

contains
  !============================================================================

  subroutine user_init_session

    use CON_planet,     ONLY: RadiusPlanet, MassPlanet
    use ModNumConst,    ONLY: cPi
    use ModPhysics,     ONLY: Si2No_V,No2Si_V,UnitRho_, &
         UnitP_, UnitX_, BodyP_I, IonFirst_
    use ModIO,          ONLY: write_myname
    use ModProcMH,      ONLY: iProc
    use ModResistivity, ONLY: Si2NoEta
    use ModGeometry,    ONLY: TypeGeometry
    CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A22,E10.3)"
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_init_session'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

!!! Introduce PlanetDensitySi etc., read those and convert
!!! from that. Initialize these to -1. PlanetRadius should be
!!! in dimensional units.

    if (TypeGeometry /= 'spherical_lnr') &
         call stop_mpi('ERROR: Correct PARAM.in, need spherical grid.')

    if(PlanetDensitySi < 0.0) &
         PlanetDensitySi  = 3.0*MassPlanet/(4.0*cPi*RadiusPlanet**3)

    if(PlanetRadius < 0.0) &
         PlanetRadius = RadiusPlanet*Si2No_V(UnitX_)

    if(PlanetPressureSi < 0.0) &
         ! PlanetPressureSi = 1.0e-8*No2Si_V(UnitP_)
         PlanetPressureSi = BodyP_I(IonFirst_)*No2Si_V(UnitP_)

    PlanetDensity           = PlanetDensitySi*Si2No_V(UnitRho_)
    PlanetPressure          = PlanetPressureSi*Si2No_V(UnitP_)
    PlanetRadiusSi          = PlanetRadius*No2Si_V(UnitX_)

    if(nLayer > 1) then
       PlanetRadiusSi_I = PlanetRadius_I*No2Si_V(UnitX_)
       Resistivity_I = ResistivitySi_I*Si2NoEta
       do iLayer=2,nLayer
          ResistivityRate(iLayer-1) = &
               (Resistivity_I(iLayer) - Resistivity_I(iLayer-1))/&
               (PlanetRadius_I(iLayer) - PlanetRadius_I(iLayer-1))
       end do
    end if

    if(iProc==0) then
       call write_myname
       write(*,*) ''
       write(*,*) '   Resistive Planet Model'
       write(*,*) '   ----------------------'
       write(*,*) ''
       write(*,"(A29,E10.3)") '  Planet density  [kg/m^3] = ',PlanetDensitySi
       write(*,"(A29,E10.3)") '  Planet pressure [Pa]     = ',PlanetPressureSi
       write(*,"(A29,E10.3)") '  Planet radius   [m]      = ',PlanetRadiusSi
       if(nLayer > 0 ) then
          write(*,*) ''
          write(*,*) '   |-------- Planet Resistivity Profile -----|'
          write(*,*) '       Radius(SI)            Resistivity(SI)'
          do iLayer =1,nLayer
             write(*,"(A7,E10.3,A15,E10.3)") " ",PlanetRadiusSi_I(iLayer)," ",&
                  ResistivitySi_I(iLayer)
          end do
       else
          write(*,*) 'Conducting Planet (eta =0)'
       end if
       write(*,*) ''
       write(*,*) ''
    end if
    call test_stop(NameSub, DoTest)
  end subroutine user_init_session
  !============================================================================

  subroutine user_read_inputs

    use ModMain
    use ModProcMH,      ONLY: iProc
    use ModReadParam
    character(len=100) :: NameCommand
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case('#PLANETDENSITY')
          call read_var('UsePlanetDensity',UsePlanetDensity)
       case('#PLANETPRESSURE')
          call read_var('UsePlanetPressure',UsePlanetPressure)
       case("#RESISTIVEPLANET")
          UseResistivePlanet = .true.
          call read_var('PlanetDensitySi'       , PlanetDensitySi)
          call read_var('PlanetPressureSi'      , PlanetPressureSi)
          call read_var('PlanetRadius'        , PlanetRadius)
          call read_var('nResistivPoints', nLayer)
          if(nLayer == 1) then
             write(*,*) ' We need minimum 2 points for including resistivity profile'
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if

          if(nLayer > 1) then
             allocate(ResistivityRate(nLayer-1),&
                  PlanetRadiusSi_I(nLayer),&
                  PlanetRadius_I(nLayer), &
                  ResistivitySi_I(nLayer),&
                  Resistivity_I(nLayer))

             do iLayer=1,nLayer
                call read_var('Radius',PlanetRadius_I(iLayer))
                call read_var('Resistivity', ResistivitySi_I(iLayer))
             end do

             ! Check values
             do iLayer=2,nLayer
                if(PlanetRadius_I(iLayer-1) < &
                     PlanetRadius_I(iLayer)) then
                   write(*,*) 'ERROR: Shoud be decreasing Radius.'
                   call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
                end if
             end do
          end if

       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) then
             write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) ' *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do

    call test_stop(NameSub, DoTest)
  end subroutine user_read_inputs
  !============================================================================

  subroutine user_set_ics(iBlock)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: R_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Bx_, Bz_, Rho_, P_, Pe_
    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP, &
         iRho, iRhoUx, iRhoUz
    use ModProcMH,     ONLY: iProc
    use ModMain,       ONLY: UseSolidState
    use ModPhysics,    ONLY: FaceState_VI, solidBc_

    integer, intent(in) :: iBlock

    integer :: i,j,k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(R_BLK(1,1,1,iBlock) > PlanetRadius) RETURN

    if(UseSolidState) then
       ! hyzhou: the new logic do not require 2 layers of cells above
       ! I may even not need this: just use set_ICs would be fine.
       ! I need to look inside the logic
       ! Instead of PlaneRadius, I may use rsolid to be consistent.
       do iFluid = 1, nFluid
          call select_fluid
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(R_BLK(i,j,k,iBlock) > PlanetRadius) CYCLE
             State_VGB(iRho,i,j,k,iBlock) = FaceState_VI(Rho_,solidBc_)
             State_VGB(iP,i,j,k,iBlock)   = FaceState_VI(P_,solidBc_)
             ! Test for MhdHypPe
             State_VGB(Pe_,i,j,k,iBlock)  = FaceState_VI(Pe_,solidBc_)
             ! hyzhou: test for timestep
             ! State_VGB(iRho,i,j,k,iBlock) = PlanetDensity
             ! State_VGB(iP,i,j,k,iBlock)   = PlanetPressure
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
          end do; end do; end do
       end do
    else
       ! old logic as Mercury
       do iFluid = 1, nFluid
          call select_fluid
          do k=1,nK; do j=1,nJ; do i=1,nI
             if(R_BLK(i+nG,j,k,iBlock) > PlanetRadius) CYCLE
             State_VGB(iRho,i,j,k,iBlock) = PlanetDensity
             State_VGB(iP,i,j,k,iBlock)   = PlanetPressure
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
          end do; end do; end do
       end do
    end if

    ! hyzhou note: may be this is problematic: why should it be 0 in the
    ! mantle? I think maybe set it to Jovian wind is more reasonable.
    do k=1,nK; do j=1,nJ; do i=1,nI
       if(R_BLK(i,j,k,iBlock) > PlanetRadius) CYCLE
       State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModGeometry,   ONLY: R_BLK, Rmin_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, p_, Bx_, Bz_, &
                             iRho_I,iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I
    use ModEnergy,     ONLY: calc_energy_cell
    use BATL_lib,      ONLY: Xyz_DGB
    use ModMain,       ONLY: UseResistivePlanet, UseSolidState
    use ModProcMH,     ONLY: iProc
    use ModPhysics,    ONLY: BodyRho_I,BodyP_I,BodyNDim_I
    use ModB0,         ONLY: B0_DGB
    use ModMultiFluid, ONLY: nFluid

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    real :: r_D(3), dRhoUr_D(3), RhoUr, u_D(3), b_D(3)
    ! hyzhou: actually I can call iFluid, iRhoUx, iRhoUz from modmultifluid
    integer :: i, iG, j, k, iFluid, iRhoUx, iRhoUz
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_cell_boundary'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(.not. UseResistivePlanet .or. TypeBc /= 'ResistivePlanet') RETURN

    ! hyzhou: do not apply this cell boundary condition if UseSolidState=.true.
    if(UseSolidState) RETURN

    if(Rmin_BLK(iBlock) <= PlanetRadius) then
       do i = 1, nI
          if(R_BLK(i+nG,1,1,iBlock) >= PlanetRadius) CYCLE
          do k = 1, nK; do j = 1, nJ;
             ! Set density, pressure and momentum inside the planet
             ! and the nG ghost cells to fixed values.
             State_VGB(Rho_,i,j,k,iBlock) = PlanetDensity
             State_VGB(P_,i,j,k,iBlock)   = PlanetPressure
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = 0.0
          end do; end do
       end do

       if(r_BLK(MaxI,1,1,iBlock) >= PlanetRadius)then
          do i = MaxI, 1, -1
             ! Find the i index just outside the planet radius
             if(R_BLK(i-1,1,1,iBlock) < PlanetRadius ) EXIT
          end do

          do k = MinK, MaxK; do j = MinJ, MaxJ
             ! Get radial velocity
             r_D = Xyz_DGB(x_:z_,i,j,k,iBlock)/ r_BLK(i,j,k,iBlock)

             RhoUr = dot_product(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock),r_D)

             ! if(RhoUr > 0.0) then
             !   ! If flow is out of the planet, remove the radial component
             !   ! of the momentum so that the flow is tangential
             !   dRhoUr_D = -r_D*RhoUr
             ! else
             !   ! If flow is into the planet the flow is absorbed
             !   dRhoUr_D = 0.0
             ! end if

             ! Set nG cells inside the planet
             ! with zero gradient boundary condition
             ! if(RhoUr > 0.0) then
             !   do iG = i-nG, i-1
             !      State_VGB(RhoUx_:RhoUz_,iG,j,k,iBlock) = &
             !           State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) + dRhoUr_D
             !      u_D = State_VGB(RhoUx_:RhoUz_,iG,j,k,iBlock)/ &
             !           State_VGB(Rho_,iG,j,k,iBlock)

                   ! Based on my (Yuxi) experience, the time-accurate
                   ! part-implicit run will crash if fixed density and
                   ! pressure are used.

                   ! State_VGB(Rho_,iG,j,k,iBlock) = 1.0
                   ! State_VGB(P_,iG,j,k,iBlock) = PlanetPressure

                    ! float BC for Pressure & density
              !     State_VGB(Rho_,iG,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
              !     State_VGB(P_,iG,j,k,iBlock)   = State_VGB(P_,i,j,k,iBlock)
                   ! reflect for Rho*Ur
              !     State_VGB(RhoUx_:RhoUz_,iG,j,k,iBlock) = &
              !           u_D*State_VGB(Rho_,iG,j,k,iBlock)
              !  end do
             ! else
                ! Float BC
             !   do iG = i-nG, i-1
             !      State_VGB(Rho_,iG,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock)
             !      State_VGB(P_,iG,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock)
             !      State_VGB(RhoUx_:RhoUz_,iG,j,k,iBlock) = &
             !           State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)
             !   end do
             ! endif

             ! Fixed density if UsePlanetDensity=.true.
             if(UsePlanetDensity) then
                ! Velocity above the surface
                u_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)/ &
                     State_VGB(Rho_,i,j,k,iBlock)
                ! Copy velocity but use the planet density for momentum
                do iG = i-nG, i-1
                   State_VGB(iRho_I,iG,j,k,iBlock) = BodyRho_I
                   ! Modify the momentum
                   State_VGB(iRhoUx_I,iG,j,k,iBlock) = &
                         u_D(x_)*State_VGB(iRho_I,iG,j,k,iBlock)
                   State_VGB(iRhoUy_I,iG,j,k,iBlock) = &
                         u_D(y_)*State_VGB(iRho_I,iG,j,k,iBlock)
                   State_VGB(iRhoUz_I,iG,j,k,iBlock) = &
                         u_D(z_)*State_VGB(iRho_I,iG,j,k,iBlock)
                end do
             else
                ! Float density
                State_VGB(iRho_I,iG,j,k,iBlock) = &
                     State_VGB(iRho_I,i,j,k,iBlock)
                ! Float momentum
                State_VGB(iRhoUx_I,iG,j,k,iBlock) = &
                     State_VGB(iRhoUx_I,i,j,k,iBlock)
                State_VGB(iRhoUy_I,iG,j,k,iBlock) = &
                     State_VGB(iRhoUy_I,i,j,k,iBlock)
                State_VGB(iRhoUz_I,iG,j,k,iBlock) = &
                     State_VGB(iRhoUz_I,i,j,k,iBlock)
             end if

             ! Fixed pressure if UsePlanePressure=.true.
             if(UsePlanetPressure) then
                do iG = i-nG, i-1
                   State_VGB(iP_I,iG,j,k,iBlock) = BodyP_I
                end do
             else
                State_VGB(iP_I,iG,j,k,iBlock) = State_VGB(iP_I,i,j,k,iBlock)
             end if

             ! Velocity perpendicular to B is float, parallel to B is 0
             do iG=i-nG,i-1
                b_D = B0_DGB(:,i,j,k,iBlock) + State_VGB(Bx_:Bz_,i,j,k,iBlock)
                b_D = b_D/max(1e-30, sqrt(sum(b_D**2)))
                do iFluid = 1, nFluid
                   iRhoUx = iRhoUx_I(iFluid); iRhoUz = iRhoUz_I(iFluid)
                   State_VGB(iRhoUx:iRhoUz,iG,j,k,iBlock) = &
                        State_VGB(iRhoUx:iRhoUz,iG,j,k,iBlock) - &
                        sum(b_D*State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock))*b_D
                end do
             end do
          end do; end do

       end if

    end if

    call calc_energy_cell(iBlock)

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_cell_boundary
  !============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModVarIndexes
    use ModFaceBoundary, ONLY: VarsTrueFace_V, B0Face_D, iBoundary
    use ModPhysics, ONLY: FaceState_VI
    use ModMultiFluid, ONLY: iUx, iUz, iUx_I, iUz_I, iFluid, nFluid

    real, intent(out):: VarsGhostFace_V(nVar)
    real ::  bUnit_D(3)
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_face_boundary'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! If I am trying to use float BC for rho and p, no need for this
    VarsGhostFace_V = FaceState_VI(:,iBoundary)

    ! VarsGhostFace_V(Hyp_) = 0

    VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_)

    ! Float density and pressure
    ! VarsGhostFace_V(Rho_) = VarsTrueFace_V(Rho_)
    ! VarsGhostFace_V(P_) = VarsTrueFace_V(P_)

    ! First use B0Face_D + VarsTrueFace_V
    ! then try B0Face_D + State_VGB(Bx_:Bz_,iGhost,jGhost,kGhost,iBlock)

    bUnit_D = B0Face_D + VarsTrueFace_V(Bx_:Bz_)
    bUnit_D = bUnit_D/max(1e-30, sqrt(sum(bUnit_D**2)))
    do iFluid = 1, nFluid
       iUx = iUx_I(iFluid); iUz = iUz_I(iFluid)
       VarsGhostFace_V(iUx:iUz) = VarsTrueFace_V(iUx:iUz) - &
            sum(bUnit_D*VarsTrueFace_V(iUx:iUz))*bUnit_D
    end do

    call test_stop(NameSub, DoTest)
  end subroutine user_set_face_boundary
  !============================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)
    use ModEnergy
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModGeometry,   ONLY: R_BLK, Rmin_BLK
    use ModResistivity, ONLY: Eta0
    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)

    integer ::i,j,k,iLayer
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_resistivity'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    Eta_G = Eta0

    if(nLayer <2 ) RETURN
    if(Rmin_BLK(iBlock) > PlanetRadius_I(1)) RETURN

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do iLayer=nLayer-1,1,-1
          if(R_BLK(i,j,k,iBlock) < PlanetRadius_I(iLayer+1) ) CYCLE
          if(R_BLK(i,j,k,iBlock) > PlanetRadius_I(iLayer) ) CYCLE
          ! to avoid eta jumps adding Eta_G
          Eta_G(i,j,k) = Eta_G(i,j,k) + Resistivity_I(iLayer+1)+&
               (R_BLK(i,j,k,iBlock)-PlanetRadius_I(iLayer+1))* &
               ResistivityRate(iLayer)
       end do
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_resistivity
  !============================================================================

end module ModUser
!==============================================================================
