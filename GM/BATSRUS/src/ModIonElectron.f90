module ModIonElectron

  ! methods related to the ion-electron fluid closures following Uri Schumlak

  use ModVarIndexes
  use ModProcMH, ONLY: iProc
  use ModMain, ONLY:  UseUserSource, &
       iTest, jTest, kTest, VarTest, BlkTest, ProcTest
  use ModAdvance, ONLY: StateOld_VCB, State_VGB, Source_VC, time_BLK
  use ModPhysics, ONLY: C2light
  use ModGeometry, ONLY: true_cell
  use ModB0,       ONLY: UseB0, B0_DGB
  use ModMultiFluid, ONLY: nIonFluid, &
       iRhoIon_I, iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I, ChargePerMass_I
  use BATL_lib,    ONLY: nI, nJ, nK, x_, y_, z_

  implicit none

  private ! except

  public:: ion_electron_source_impl
  public:: ion_electron_init_point_impl
  public:: read_ion_electron_param
  
  ! calculate analytic Jacobian for point-implicit scheme
  logical, parameter:: IsAnalyticJacobian = .true.

contains
  !===========================================================================
  subroutine read_ion_electron_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_ion_electron_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_ion_electron_param
  !===========================================================================
  subroutine ion_electron_source_impl(iBlock)

    ! Add electron, ion and electric field source terms to Source_VC
    ! Calculate dS/dU Jacobian

    use ModPointImplicit, ONLY:  UsePointImplicit, IsPointImplSource, &
         IsPointImplPerturbed, DsDu_VVC

    integer, intent(in)          :: iBlock

    real:: State_V(nVar), b_D(3)
    integer:: i, j, k
    integer:: iIon, iRhoUx, iRhoUy, iRhoUz
    real:: ChargePerMass
    logical :: DoTest, DoTestMe, DotestCell
    character(len=*), parameter :: NameSub = 'ion_electron_source_impl'
    !-----------------------------------------------------------------------
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if
    DoTestCell = .false.

    ! Add user defined point implicit source terms here
    ! Explicit user sources are added in calc_sources
    if(UsePointImplicit .and. UseUserSource) call user_calc_sources(iBlock)

    ! Do not evaluate multi-ion sources in the numerical Jacobian calculation
    ! (needed for the user source terms) 
    if(IsPointImplPerturbed .and. IsAnalyticJacobian) RETURN

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       State_V = State_VGB(:,i,j,k,iBlock)

       b_D = State_V(Bx_:Bz_)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       ! d(rhou)/dt += q/m*(rho*E + rhou x B)
       Source_VC(iRhoUxIon_I,i,j,k) = Source_VC(iRhoUxIon_I,i,j,k) &
            + ChargePerMass_I *                &
            ( State_V(iRhoIon_I)*State_V(Ex_)  &
            + State_V(iRhoUyIon_I)*b_D(z_)     &
            - State_V(iRhoUzIon_I)*b_D(y_) )

       Source_VC(iRhoUyIon_I,i,j,k) = Source_VC(iRhoUyIon_I,i,j,k) &
            + ChargePerMass_I *                &
            ( State_V(iRhoIon_I)*State_V(Ey_)  &
            + State_V(iRhoUzIon_I)*b_D(x_)     &
            - State_V(iRhoUxIon_I)*b_D(z_) )

       Source_VC(iRhoUzIon_I,i,j,k) = Source_VC(iRhoUzIon_I,i,j,k) &
            + ChargePerMass_I *                &
            ( State_V(iRhoIon_I)*State_V(Ez_)  &
            + State_V(iRhoUxIon_I)*b_D(y_)     &
            - State_V(iRhoUyIon_I)*b_D(x_) )

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

    end do; end do; end do

  end subroutine ion_electron_source_impl
  !===========================================================================
  subroutine ion_electron_init_point_impl

    ! Select variables for point implicit evaluation. This is the union
    ! of the ion and electron momenta, the electric field
    ! and the variables selected (if any) in ModUser::user_init_point_implicit
    
    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    use ModUserInterface ! user_init_point_implicit

    logical :: IsPointImpl_V(nVar)
    integer :: iVar, iPointImplVar, nPointImplVar
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

  end subroutine ion_electron_init_point_impl

end module ModIonElectron
