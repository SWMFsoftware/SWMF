module ModIonElectron

  use BATL_lib, ONLY: &
       test_start, test_stop, iTest, jTest, kTest, iVarTest

  ! methods related to the ion-electron fluid closures following Uri Schumlak

  use ModVarIndexes
  use ModMain, ONLY:  UseUserSource
  use ModAdvance, ONLY: State_VGB, Source_VC
  use ModPhysics, ONLY: C2light
  use ModGeometry, ONLY: true_cell
  use ModB0,       ONLY: UseB0, B0_DGB
  use ModMultiFluid, ONLY: nIonFluid, nTrueIon, ElectronFirst_, &
       iRhoIon_I, iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I,   &
       iPIon_I, ChargePerMass_I, MassIon_I, ChargeIon_I
  use BATL_lib,    ONLY: nI, nJ, nK, x_, y_, z_

  implicit none

  private ! except

  public:: ion_electron_source_impl
  public:: ion_electron_init_point_impl
  public:: read_ion_electron_param
  public:: correct_electronfluid_efield

  real, public:: HypEDecay = 0.1

  logical,public :: DoCorrectElectronFluid = .false.
  logical,public :: DoCorrectEfield        = .false.

  ! calculate analytic Jacobian for point-implicit scheme
  logical, parameter:: IsAnalyticJacobian = .true.

  integer, allocatable, public :: iVarUseCmax_I(:)

contains
  !============================================================================
  subroutine read_ion_electron_param(NameCommand)

    use ModMain,      ONLY: NameVarLower_V
    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: split_string

    character(len=*), intent(in):: NameCommand

    integer :: iVar, jVar, nVarUseCmax
    character(len=200) :: StringVarUseCmax
    character(len=20)  :: NameVarUseCmax_I(nVar)

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'read_ion_electron_param'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    select case(NameCommand)
    case("#HYPERBOLICDIVE")
       call read_var('HypEDecay', HypEDecay)
    case("#CORRECTELECTRONFLUID")
       call read_var('DoCorrectElectronFluid', DoCorrectElectronFluid)
    case("#CORRECTEFIELD")
       call read_var('DoCorrectEfield', DoCorrectEfield)
    case("#CMAXDIFFUSION")
       if(allocated(iVarUseCmax_I)) deallocate(iVarUseCmax_I)
       call read_var('StringVarUseCmax', StringVarUseCmax, IsLowerCase=.true.)
       if (StringVarUseCmax == 'all') then
          allocate(iVarUseCmax_I(nVar))
          iVarUseCmax_I = (/(iVar, iVar = 1, nVar)/)
       else if (StringVarUseCmax == 'default') then
          allocate(iVarUseCmax_I(4))
          iVarUseCmax_I = (/Ex_, Ey_, Ez_, HypE_/)
       else
          call split_string(StringVarUseCmax, NameVarUseCmax_I, nVarUseCmax)
          allocate(iVarUseCmax_I(nVarUseCmax))
          do iVar = 1, nVarUseCmax
             call get_iVar(NameVarUseCmax_I(iVar), iVarUseCmax_I(iVar))
          end do
       end if
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

    call test_stop(NameSub, DoTest)
  end subroutine read_ion_electron_param
  !============================================================================
  subroutine ion_electron_source_impl(iBlock)

    ! Add electron, ion and electric field source terms to Source_VC
    ! Calculate dS/dU Jacobian

    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplPerturbed, DsDu_VVC

    integer, intent(in)          :: iBlock

    real:: State_V(nVar), b_D(3)
    integer:: i, j, k, iVar
    integer:: iIon, iRhoUx, iRhoUy, iRhoUz
    real:: ChargePerMass

    logical :: DotestCell
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'ion_electron_source_impl'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    DoTestCell = .false.

    if(UsePointImplicit .and. .not. IsPointImplSource) then
       if (DoTest) write(*,*) NameSub, ' initial RETURN'
       RETURN
    end if

    ! Add user defined point implicit source terms here
    ! Explicit user sources are added in calc_sources
    if(UsePointImplicit .and. UseUserSource) call user_calc_sources(iBlock)

    ! Do not evaluate multi-ion sources in the numerical Jacobian calculation
    ! (needed for the user source terms)
    if(IsPointImplPerturbed .and. IsAnalyticJacobian) then
       if (DoTest) write(*,*) &
            NameSub, ' no evaluation in the numerical Jacobian calculation'
       RETURN
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       DoTestCell = DoTest .and. i==iTest .and. j==jTest .and. k==kTest

       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       State_V = State_VGB(:,i,j,k,iBlock)

       b_D = State_V(Bx_:Bz_)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       if (DoTestCell) then
          write(*,'(1x,2a,es20.12)') NameSub,' initial Source(testvar) =', &
               Source_VC(iVarTest,i,j,k)
       end if

       ! d(rhou)/dt += q/m*(rho*E + rhou x B)
       Source_VC(iRhoUxIon_I,i,j,k) = Source_VC(iRhoUxIon_I,i,j,k) &
            + ChargePerMass_I *               &
            ( State_V(iRhoIon_I)*State_V(Ex_) &
            + State_V(iRhoUyIon_I)*b_D(z_)    &
            - State_V(iRhoUzIon_I)*b_D(y_) )

       Source_VC(iRhoUyIon_I,i,j,k) = Source_VC(iRhoUyIon_I,i,j,k) &
            + ChargePerMass_I *               &
            ( State_V(iRhoIon_I)*State_V(Ey_) &
            + State_V(iRhoUzIon_I)*b_D(x_)    &
            - State_V(iRhoUxIon_I)*b_D(z_) )

       Source_VC(iRhoUzIon_I,i,j,k) = Source_VC(iRhoUzIon_I,i,j,k) &
            + ChargePerMass_I *               &
            ( State_V(iRhoIon_I)*State_V(Ez_) &
            + State_V(iRhoUxIon_I)*b_D(y_)    &
            - State_V(iRhoUyIon_I)*b_D(x_) )

       if (DoTestCell) then
          do iIon = 1, nIonFluid
             write(*,'(1x,2a,i3)')        NameSub,' iIon   =', iIon
             write(*,'(1x,2a,10es20.12)') NameSub,' uIon_D =', &
                  State_V(iRhoUxIon_I(iIOn))/State_V(iRhoIon_I(iIon)), &
                  State_V(iRhoUyIon_I(iIOn))/State_V(iRhoIon_I(iIon)), &
                  State_V(iRhoUzIon_I(iIOn))/State_V(iRhoIon_I(iIon))
             write(*,'(1x,2a,10es20.12)') NameSub,' E_D    =', State_V(Ex_:Ez_)
             write(*,'(1x,2a,10es20.12)') NameSub,' b_D    =', b_D
          end do
       end if

       ! dE/dt += -c^2*J = -c^2*sum(q/m*rho*u)
       Source_VC(Ex_,i,j,k) = &
            -C2light*sum(ChargePerMass_I * State_V(iRhoUxIon_I))
       Source_VC(Ey_,i,j,k) = &
            -C2light*sum(ChargePerMass_I * State_V(iRhoUyIon_I))
       Source_VC(Ez_,i,j,k) = &
            -C2light*sum(ChargePerMass_I * State_V(iRhoUzIon_I))

       ! Set corresponding matrix elements
       if (IsAnalyticJacobian .and. UsePointImplicit) then
          DsDu_VVC(iRhoUxIon_I,Ex_,i,j,k) = DsDu_VVC(iRhoUxIon_I,Ex_,i,j,k) &
               + ChargePerMass_I * State_V(iRhoIon_I)
          DsDu_VVC(iRhoUyIon_I,Ey_,i,j,k) = DsDu_VVC(iRhoUyIon_I,Ey_,i,j,k) &
               + ChargePerMass_I * State_V(iRhoIon_I)
          DsDu_VVC(iRhoUzIon_I,Ez_,i,j,k) = DsDu_VVC(iRhoUzIon_I,Ez_,i,j,k) &
               + ChargePerMass_I * State_V(iRhoIon_I)

          DsDu_VVC(Ex_,iRhoUxIon_I,i,j,k) = DsDu_VVC(Ex_,iRhoUxIon_I,i,j,k) &
               - C2light*ChargePerMass_I
          DsDu_VVC(Ey_,iRhoUyIon_I,i,j,k) = DsDu_VVC(Ey_,iRhoUyIon_I,i,j,k) &
               - C2light*ChargePerMass_I
          DsDu_VVC(Ez_,iRhoUzIon_I,i,j,k) = DsDu_VVC(Ez_,iRhoUzIon_I,i,j,k) &
               - C2light*ChargePerMass_I

          do iIon = 1, nIonFluid
             iRhoUx = iRhoUxIon_I(iIon)
             iRhoUy = iRhoUyIon_I(iIon)
             iRhoUz = iRhoUzIon_I(iIon)
             ChargePerMass = ChargePerMass_I(iIon)
             DsDu_VVC(iRhoUx,iRhoUy,i,j,k) = DsDu_VVC(iRhoUx,iRhoUy,i,j,k) &
                  + ChargePerMass*b_D(z_)
             DsDu_VVC(iRhoUx,iRhoUz,i,j,k) = DsDu_VVC(iRhoUx,iRhoUz,i,j,k) &
                  - ChargePerMass*b_D(y_)
             DsDu_VVC(iRhoUy,iRhoUz,i,j,k) = DsDu_VVC(iRhoUy,iRhoUz,i,j,k) &
                  + ChargePerMass*b_D(x_)
             DsDu_VVC(iRhoUy,iRhoUx,i,j,k) = DsDu_VVC(iRhoUy,iRhoUx,i,j,k) &
                  - ChargePerMass*b_D(z_)
             DsDu_VVC(iRhoUz,iRhoUx,i,j,k) = DsDu_VVC(iRhoUz,iRhoUx,i,j,k) &
                  + ChargePerMass*b_D(y_)
             DsDu_VVC(iRhoUz,iRhoUy,i,j,k) = DsDu_VVC(iRhoUz,iRhoUy,i,j,k) &
                  - ChargePerMass*b_D(x_)
          end do
       end if

       if (DoTestCell) then
          write(*,'(1x,2a,10es20.12)') &
               NameSub, ' ChargePerMass_I, net charge =',&
               ChargePerMass_I, sum(ChargePerMass_I*State_V(iRhoIon_I))
          write(*,*) NameSub, ' Source_VC='
          do iVar = 1, nVar
             write(*,'(2x,a,100es20.12)') &
                  NameVar_V(iVar), Source_VC(iVar,iTest,jTest,kTest)
          end do
       end if

    end do; end do; end do

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine ion_electron_source_impl
  !============================================================================
  subroutine ion_electron_init_point_impl

    ! Select variables for point implicit evaluation. This is the union
    ! of the ion and electron momenta, the electric field
    ! and the variables selected (if any) in ModUser::user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, &
         init_point_implicit_num
    use ModUserInterface ! user_init_point_implicit

    logical :: IsPointImpl_V(nVar)
    integer :: iVar, iPointImplVar, nPointImplVar

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'ion_electron_init_point_impl'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)
    IsPointImpl_V = .false.
    IsPointImplMatrixSet = IsAnalyticJacobian

    if(UseUserSource)then
       call user_init_point_implicit
       if(allocated(iVarPointImpl_I)) then
          IsPointImpl_V(iVarPointImpl_I) = .true.
          ! Set index array for variables to be perturbed
          if(.not.IsPointImplMatrixSet) call init_point_implicit_num
          deallocate(iVarPointImpl_I)
       end if
    end if

    ! All electron and ion momenta and the electric fields are implicit
    IsPointImpl_V(iRhoUxIon_I) = .true.
    IsPointImpl_V(iRhoUyIon_I) = .true.
    IsPointImpl_V(iRhoUzIon_I) = .true.
    IsPointImpl_V(Ex_:Ez_)     = .true.

    nPointImplVar = count(IsPointImpl_V)

    allocate(iVarPointImpl_I(nPointImplVar))

    iPointImplVar = 0
    do iVar = 1, nVar
       if(.not. IsPointImpl_V(iVar)) CYCLE
       iPointImplVar = iPointImplVar + 1
       iVarPointImpl_I(iPointImplVar) = iVar
    end do

    call test_stop(NameSub, DoTest)
  end subroutine ion_electron_init_point_impl
  !============================================================================

  subroutine correct_electronfluid_efield(State_VG, iMin, iMax, jMin, jMax, &
       kMin, kMax, iBlock, DoHallCurrentIn, DoGradPeIn, DoCorrectPeIn,      &
       DoCorrectEfieldIn)

    ! The subroutine will correct the electron fluid and electric field
    ! based on the Hall MHD case.

    use ModSize,       ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nG
    use ModMultiFluid, ONLY: iPIon_I
    use ModB0,         ONLY: UseB0, B0_DGB
    use ModCellGradient,   ONLY: calc_gradient
    use ModCurrent,        ONLY: get_current

    real,    intent(inout) :: State_VG(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    integer, intent(in)    :: iMin, iMax, jMin, jMax, kMin, kMax
    integer, intent(in)    :: iBlock
    logical, optional, intent(in) :: DoHallCurrentIn, DoGradPeIn
    logical, optional, intent(in) :: DoCorrectPeIn
    logical, optional, intent(in) :: DoCorrectEfieldIn

    integer :: i,j,k
    real    :: Current_D(3), GradPe_D(3), B0_D(3)
    real, allocatable :: GradPe_DG(:,:,:,:)
    logical :: DoHallCurrent, DoGradPe, DoCorrectEfield, DoCorrectPe

    logical :: DoTestCell

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'correct_electronfluid_efield'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! only single electron fluid is supported at this moment
    if(nElectronFluid /= 1) &
         call stop_mpi(NameSub//' only nElectronFluid = 1 is supported')

    DoHallCurrent = .false.
    if (present(DoHallCurrentIn)) DoHallCurrent = DoHallCurrentIn

    DoGradPe = .false.
    if (present(DoGradPeIn)) DoGradPe = DoGradPeIn

    DoCorrectPe = .false.
    if (present(DoCorrectPeIn)) DoCorrectPe = DoCorrectPeIn

    DoCorrectEfield = .false.
    if (present(DoCorrectEfieldIn)) DoCorrectEfield = DoCorrectEfieldIn

    if (DoGradPe) then
       if(.not. allocated(GradPe_DG)) &
            allocate(GradPe_DG(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK))
       call calc_gradient(iBlock, State_VG(iPIon_I(ElectronFirst_),:,:,:),nG, &
            GradPe_DG)
    else
       GradPe_D = 0.0
    end if

    do k = kMin, kMax; do j = jMin, jMax; do i = iMin, iMax
       DoTestCell = DoTest .and. i == iTest .and. j == jTest .and. k == kTest

       if (DoHallCurrent .and. true_cell(i,j,k,iBlock)) then
          call get_current(i,j,k,iBlock,Current_D)
       else
          Current_D = 0.0
       end if
       if (DoGradPe) GradPe_D = GradPe_DG(:,i,j,k)

       B0_D = 0.0
       if (UseB0) B0_D = B0_DGB(:,i,j,k,iBlock)

       call correct_electronfluid_efield_cell(State_VG(:,i,j,k), B0_D,     &
            Current_D, GradPe_D, DoCorrectPe, DoCorrectEfield, DoTestCell)
    end do; end do; end do

    if (allocated(GradPe_DG)) deallocate(GradPe_DG)

    call test_stop(NameSub, DoTest)
  end subroutine correct_electronfluid_efield
  !============================================================================

  subroutine correct_electronfluid_efield_cell(State_V, B0_D, Current_D,  &
       GradPe_D, DoCorrectPe, DoCorrectEfield, DoTest)

    use ModPhysics,        ONLY: ElectronCharge
    use ModCoordTransform, ONLY: cross_product

    real,    intent(inout) :: State_V(nVar)
    real,    intent(in)    :: B0_D(3)
    real,    intent(in)    :: Current_D(3)
    real,    intent(in)    :: GradPe_D(3)
    logical, intent(in)    :: DoCorrectPe
    logical, intent(in)    :: DoCorrectEfield
    logical, intent(in)    :: DoTest

    integer :: iVar
    real    :: nElec, InvElectronDens, uElec_D(3), Bfield_D(3), Efield_D(3)
    real    :: uPlus_D(3), Te
    real    :: StateOld_V(nVar)

    integer, parameter :: ElecUx_ = iRhoUxIon_I(ElectronFirst_),  &
         ElecUz_ = iRhoUzIon_I(ElectronFirst_)

    character(len=*), parameter:: NameSub = 'correct_electronfluid_efield_cell'
    !--------------------------------------------------------------------------

    ! save the old state_V for testing
    StateOld_V = State_V

    ! original electron temperature
    Te = State_V(iPIon_I(ElectronFirst_))*MassIon_I(ElectronFirst_) &
         /State_V(iRhoIon_I(ElectronFirst_))

    ! inv of electron charge density = 1/(e*ne)
    InvElectronDens = 1.0/sum( &
         ChargePerMass_I(1:nTrueIon)*State_V(iRhoIon_I(1:nTrueIon)) )

    ! electron number density
    nElec = ElectronCharge/InvElectronDens

    ! charge average ion velocity
    uPlus_D(x_) = InvElectronDens*sum(ChargePerMass_I(1:nTrueIon) &
         *State_V(iRhoUxIon_I(1:nTrueIon)))
    uPlus_D(y_) = InvElectronDens*sum(ChargePerMass_I(1:nTrueIon) &
         *State_V(iRhoUyIon_I(1:nTrueIon)))
    uPlus_D(z_) = InvElectronDens*sum(ChargePerMass_I(1:nTrueIon) &
         *State_V(iRhoUzIon_I(1:nTrueIon)))

    ! electron velocity = uPlus + u_H, uH = -j/ne/e
    uElec_D  = uPlus_D - Current_D*InvElectronDens

    ! calculate the electric field E = -ue x B - gradpe/ne/e
    Bfield_D = State_V(Bx_:Bz_) + B0_D
    Efield_D = -cross_product(uElec_D,Bfield_D) - GradPe_D*InvElectronDens

    ! correct the electron mass density and velocity
    State_V(iRhoIon_I(ElectronFirst_))  = nElec*MassFluid_I(ElectronFirst_)
    State_V(ElecUx_:ElecUz_)       = State_V(iRhoIon_I(ElectronFirst_))*uElec_D

    ! correct the electron pressure based on the original Te if needed
    if (DoCorrectPe) State_V(iPIon_I(ElectronFirst_)) = nElec*Te

    ! correct the electric field
    if (DoCorrectEfield) State_V(Ex_:Ez_) = Efield_D

    if (DoTest) then
       write(*,*) NameSub, ' DoCorrectEfield =', DoCorrectEfield
       write(*,'(1x,2a,10es13.5)')  NameSub, ' MassIon_I   =', MassIon_I
       write(*,'(1x,2a,10es13.5)')  NameSub, ' ChargeIon_I =', ChargeIon_I
       write(*,'(1x,2a,3es20.12)')  NameSub, ' uPlus_D     =', uPlus_D
       write(*,'(1x,2a,3es20.12)')  NameSub, ' Current_D   =', Current_D
       write(*,'(1x,2a,3es20.12)')  NameSub, ' uElec_D     =', uElec_D
       write(*,'(1x,2a,3es20.12)')  NameSub, ' Bfield_D    =', Bfield_D
       write(*,'(1x,2a,3es20.12)')  NameSub, ' Efield_D    =', Efield_D
       write(*,*) NameSub, ' NameVar, State_V (old/new) ='
       do iVar = 1, nVar
          write(*,'(1x,1a,3es20.12)') &
               NameVar_V(iVar), StateOld_V(iVar), State_V(iVar)
       end do
    end if
  end subroutine correct_electronfluid_efield_cell
  !============================================================================

end module ModIonElectron
!==============================================================================
