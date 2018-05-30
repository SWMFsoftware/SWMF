!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iBlockTest
  ! User module for Ganymede.
  ! Must compile with MHDHypPe equation set.
  ! Version 2.1 can work with both Cartesian and spherical grid.
  use ModUserEmpty,                          &
       IMPLEMENTED1 => user_init_session,    &
       IMPLEMENTED2 => user_set_ics,         &
       IMPLEMENTED3 => user_set_resistivity, &
       IMPLEMENTED4 => user_read_inputs,     &
       IMPLEMENTED5 => user_set_face_boundary,&
       IMPLEMENTED6 => user_set_cell_boundary

  include 'user_module.h' ! list of public methods

  real,              parameter :: VersionUserModule = 2.1
  character (len=*), parameter :: NameUserModule = 'Ganymede, Hongyang Zhou'

  real :: PlanetRadius=-1.
  real :: PlanetRadiusSi=-1.
  character (len=10) :: InitType = 'B1U1'

  integer :: iLayer
  integer :: nLayer =0 ! Number of points in planet resistivity profile
  real, allocatable :: PlanetRadiusSi_I(:),PlanetRadius_I(:),&
       ResistivitySi_I(:),Resistivity_I(:)
  real, allocatable :: ResistivityRate(:)

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

    if(iProc==0)then
       select case(TypeGeometry)
       case('spherical_lnr','spherical')
          write(*,*) 'Using spherical grid...'
       case('cartesian','rotatedcartesian')
          write(*,*) 'Using Cartesian grid...'
       case default
          call stop_mpi('ERROR: need spherical/ Cartesian grid!')
       end select
    end if

    if(PlanetRadius < 0.0) &
         PlanetRadius = RadiusPlanet*Si2No_V(UnitX_)

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
          write(*,*) 'Conducting Planet (eta = 0)'
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

       case('#INITTYPE')
          call read_var('InitType', InitType)

       case("#RESISTIVEPLANET")
          UseResistivePlanet = .true.
          call read_var('PlanetRadius'   , PlanetRadius)
          call read_var('nResistivPoints', nLayer)
          if(nLayer == 1) then
             write(*,*) ' We need minimum 2 points for including resistivity &
                  &profile'
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
                   write(*,*) 'ERROR: Should be decreasing Radius.'
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
    use ModGeometry,   ONLY: Rmin_BLK, R_BLK
    use ModSize,       ONLY: nI, nJ, nK, nG
    use ModVarIndexes, ONLY: Bx_, Bz_, Rho_, P_, Pe_
    use ModMultiFluid, ONLY: select_fluid, iFluid, nFluid, iP, &
         iRho, iRhoUx, iRhoUz
    use ModMain,       ONLY: UseSolidState,Coord1MinBc_,Coord1MaxBc_,solidBc_
    use ModPhysics,    ONLY: FaceState_VI, CellState_VI
    use BATL_lib,      ONLY: IsCartesian

    integer, intent(in) :: iBlock

    integer :: i,j,k
    integer :: Bc_
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_ics'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    if(IsCartesian)then
       Bc_ = Coord1MinBc_
    else
       Bc_ = Coord1MaxBc_
    end if

    select case(InitType)
    case('B1U1')
       ! Initialize mantle region
       if(Rmin_BLK(iBlock) <= PlanetRadius) then
          do iFluid = 1, nFluid
             call select_fluid
             do k=1,nK; do j=1,nJ; do i=1,nI
                if(R_BLK(i,j,k,iBlock) > PlanetRadius) CYCLE
                State_VGB(iRho,i,j,k,iBlock) = FaceState_VI(Rho_,solidBc_)
                State_VGB(iP,i,j,k,iBlock)   = FaceState_VI(P_,solidBc_)
                State_VGB(Pe_,i,j,k,iBlock)  = FaceState_VI(Pe_,solidBc_)
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             end do; end do; end do
          end do
       end if

       ! Magnetic field inside the mantle is dipole
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i,j,k,iBlock) > PlanetRadius)then
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  CellState_VI(Bx_:Bz_,Bc_)
          else
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
          end if
       end do; end do; end do

       ! Init for velocity
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid
             if(R_BLK(i,j,k,iBlock) < PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             else
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_)
             end if
          end do
       end do; end do; end do

       ! Uniform density         
       State_VGB(Rho_,:,:,:,iBlock) = CellState_VI(Rho_,Bc_)
       ! Uniform pressure                      
       State_VGB(p_ ,:,:,:,iBlock)  = CellState_VI(p_,Bc_)
       State_VGB(pe_,:,:,:,iBlock)  = CellState_VI(pe_,Bc_)

    case('B1U01')
       ! Initialize mantle region 
       if(Rmin_BLK(iBlock) <= PlanetRadius) then
          do iFluid = 1, nFluid
             call select_fluid
             do k=1,nK; do j=1,nJ; do i=1,nI
                if(R_BLK(i,j,k,iBlock) > PlanetRadius) CYCLE
                State_VGB(iRho,i,j,k,iBlock) = FaceState_VI(Rho_,solidBc_)
                State_VGB(iP,i,j,k,iBlock)   = FaceState_VI(P_,solidBc_)
                State_VGB(Pe_,i,j,k,iBlock)  = FaceState_VI(Pe_,solidBc_)
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             end do; end do; end do
          end do
       end if
             
       ! Magnetic field inside the mantle is dipole
       do k=1,nK; do j=1,nJ; do i=1,nI
          if(R_BLK(i,j,k,iBlock) > PlanetRadius)then
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  CellState_VI(Bx_:Bz_,Bc_)
          else
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
          end if
       end do; end do; end do
       
       ! Smooth init for velocity
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid
             if(R_BLK(i,j,k,iBlock) < 2.5*PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             elseif(R_BLK(i,j,k,iBlock) < 5.0*PlanetRadius)then
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_) *&
                     (R_BLK(i,j,k,iBlock)/2.5 - 1)
             else
                State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
                     CellState_VI(iRhoUx:iRhoUz,Bc_)
             end if
          end do
       end do; end do; end do
       
       ! Uniform density
       State_VGB(Rho_,:,:,:,iBlock) = CellState_VI(Rho_,Bc_)
       ! Uniform pressure
       State_VGB(p_ ,:,:,:,iBlock)  = CellState_VI(p_,Bc_)
       State_VGB(pe_,:,:,:,iBlock)  = CellState_VI(pe_,Bc_)

    case('B0U0') ! Upstream propagation 
       ! Magnetic field starts from dipole
       do k=1,nK; do j=1,nJ; do i=1,nI
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = 0.0
       end do; end do; end do
       
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iFluid=1,nFluid
             call select_fluid
             ! Velocity starts from zeros 
             State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = 0.0
             ! Uniform density
             State_VGB(iRho,i,j,k,iBlock) = CellState_VI(rho_,Bc_)
             ! Uniform pressure
             State_VGB(P_,i,j,k,iBlock) = CellState_VI(P_,Bc_)
             State_VGB(Pe_,i,j,k,iBlock) = CellState_VI(Pe_,Bc_)
          end do
       end do; end do; end do
       
    case default
       call stop_MPI('ERROR: unknown initialization type!')
    end select

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_ics
  !============================================================================

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    use ModAdvance,      ONLY: State_VGB
    use ModCellBoundary, ONLY: iMin, iMax, jMin, jMax, kMin, kMax
    use BATL_lib,        ONLY: Xyz_DGB
    use ModPhysics,      ONLY: CellState_VI
    use ModSize,         ONLY: nI
    use ModEnergy,       ONLY: calc_energy_cell

    integer,          intent(in)  :: iBlock, iSide
    character(len=*), intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer :: i,j,k

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_set_cell_boundary'
    !--------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)    

    ! This is only used for spherical grid with no outer face BC applied.
    ! If outer face BC is used, set the second parameter in command
    ! #OUTERBOUNDARY to none to avoid calling this function.

    ! Outer boundary condition                                             
    if(iSide == 2)then
       do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax
          if(Xyz_DGB(1,i,j,k,iBlock) < 0)then
             ! Upstream fixed                            
             State_VGB(:,i,j,k,iBlock) = CellState_VI(:,iSide)
          else
             ! Downstream float                                       
             State_VGB(:,i,j,k,iBlock) = State_VGB(:,nI,j,k,iBlock)
          end if
       end do; end do; end do

       IsFound = .true.
       RETURN
    end if
    
    ! What is this actually?
    call calc_energy_cell(iBlock)
    
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_cell_boundary

  !============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModVarIndexes
    use ModFaceBoundary, ONLY: VarsTrueFace_V, B0Face_D, iBoundary, &
         iFace, jFace, kFace, iside, iBlock => iBlockBc
    use ModPhysics,      ONLY: FaceState_VI
    use ModMultiFluid,   ONLY: iUx, iUz, iUx_I, iUz_I, iFluid, nFluid, iRho, iP
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: UseHyperbolicDivb

    real, intent(out):: VarsGhostFace_V(nVar)

    !integer:: iTrue, jTrue, kTrue, iBody, jBody, kBody
    integer:: i, j, k
    real ::  bUnit_D(3)
    character(len=*), parameter:: NameSub = 'user_set_face_boundary'
    !--------------------------------------------------------------------------

    ! Copy everything from physical face to ghost face
    ! rho, u, b, p, hyp
    VarsGhostFace_V = VarsTrueFace_V    

    !if(UseHyperbolicDivb) &
    !     VarsGhostFace_V(Hyp_) = 0

    ! For 1st order floating BC of magnetic field, the true face value is 
    ! the same as the first layer of physical cell center values.
    ! For 2nd order floating BC of magnetic field, extrapolation is needed,
    ! which requires two layers of physical cells' center values.

    ! Get ghost cell indexes from face indexes
    !i = iFace; j = jFace; k = kFace;
    !select case(iSide)
    !case(1)
    !   i = iFace - 1
    !case(3)
    !   j = jFace - 1
    !case(5)
    !   k = kFace - 1
    !end select

    ! copy the internal cell center value to the face
    !VarsGhostFace_V(Bx_:Bz_) = State_VGB(Bx_:Bz_,iFace,j,k,iBlock)

    ! Float for B
    VarsGhostFace_V(Bx_:Bz_) = VarsTrueFace_V(Bx_:Bz_)
    
    ! Fixed rho and p
    VarsGhostFace_V(rho_) = FaceState_VI(rho_,iBoundary)
    VarsGhostFace_V(p_)   = FaceState_VI(p_,iBoundary)
    VarsGhostFace_V(pe_)  = FaceState_VI(pe_,iBoundary)

    ! First use B0Face_D + VarsTrueFace_V
    ! then try B0Face_D + State_VGB(Bx_:Bz_,iGhost,jGhost,kGhost,iBlock)
    bUnit_D = B0Face_D + VarsTrueFace_V(Bx_:Bz_)
    bUnit_D = bUnit_D/max(1e-30, sqrt(sum(bUnit_D**2)))
    do iFluid = 1, nFluid
       iUx = iUx_I(iFluid); iUz = iUz_I(iFluid)
       VarsGhostFace_V(iUx:iUz) = VarsTrueFace_V(iUx:iUz) - &
            sum(bUnit_D*VarsTrueFace_V(iUx:iUz))*bUnit_D
    end do

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

    if(nLayer < 2) RETURN
    if(Rmin_BLK(iBlock) > PlanetRadius_I(1)) RETURN

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do iLayer=nLayer-1,1,-1
          if(R_BLK(i,j,k,iBlock) < PlanetRadius_I(iLayer+1) ) CYCLE
          if(R_BLK(i,j,k,iBlock) > PlanetRadius_I(iLayer) ) CYCLE
          ! Avoid eta jumps adding Eta_G
          Eta_G(i,j,k) = Eta_G(i,j,k) + Resistivity_I(iLayer+1)+&
               (R_BLK(i,j,k,iBlock)-PlanetRadius_I(iLayer+1))* &
               ResistivityRate(iLayer)
       end do
    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_set_resistivity
  !============================================================================

!  subroutine user_set_cell_boundary(iBlock, iSide, TypeBC, IsFound)
!
!    integer,          intent(in)  :: iBlock, iSide
!    character(len=*), intent(in)  :: TypeBc
!    logical,          intent(out) :: IsFound
!    !--------------------------------------------------------------------------
!    return
!
!  end subroutine user_set_cell_boundary

end module ModUser
!==============================================================================
