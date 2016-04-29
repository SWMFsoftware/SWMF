module ModIonElectron

  use ModVarIndexes
  use ModAdvance, ONLY: StateOld_VCB, State_VGB, Source_VC, time_BLK
  use ModPhysics, ONLY: C2light
  use ModGeometry, ONLY: true_cell
  use ModB0,       ONLY: UseB0, B0_DGB
  use ModMultiFluid, ONLY: nIonFluid, &
       iRhoIon_I, iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I, ChargePerMass_I
  use BATL_lib,    ONLY: nI, nJ, nK, x_, y_, z_
  
  implicit none

  real:: BetaImpl = 1.0

contains
  !===========================================================================
  subroutine read_ion_electron_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'read_ion_electron_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case('#IMPLIONELECTRON')
       call read_var('BetaImpl', BetaImpl)
    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_ion_electron_param
  !===========================================================================
  subroutine add_ion_electron_source(iBlock, UseOld)

    ! Add electron, ion and electric field source terms to Source_VC
    ! If UseOld is present, use StateOld_VCB, otherwise State_VGB

    integer, intent(in)          :: iBlock
    logical, intent(in), optional:: UseOld

    real:: State_V(nVar), b_D(3)

    integer:: i, j, k
    !-----------------------------------------------------------------------

    ! For point-implicit scheme don't add explicit source at all
    if(BetaImpl > 0.0 .and. .not.present(UseOld)) RETURN

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       if(BetaImpl > 0.0)then
          State_V = StateOld_VCB(:,i,j,k,iBlock)
       else
          State_V = State_VGB(:,i,j,k,iBlock)
       end if

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

       ! de_s/dt += j_s.E = (q/m)_s*rhou_s.E 
       Source_VC(Energy_:Energy_+nIonFluid-1,i,j,k) =              &
            Source_VC(Energy_:Energy_+nIonFluid-1,i,j,k)           &
            + ChargePerMass_I * (State_V(iRhoUxIon_I)*State_V(Ex_) &
            +                    State_V(iRhoUyIon_I)*State_V(Ey_) &
            +                    State_V(iRhoUzIon_I)*State_V(Ez_) )

    end do; end do; end do

  end subroutine add_ion_electron_source

  !===========================================================================
  subroutine update_impl_ion_electron(iBlock)

    use ModCoordTransform, ONLY: cross_product
    use ModMain, ONLY: Cfl
    use ModEnergy, ONLY: calc_energy

    integer, intent(in):: iBlock

    ! Solve for ion and electron momenta and electric field update implicitly

    integer:: i, j, k, iIon, iRho, iRhoUx, iRhoUz

    ! Fluid density and full magnetic field (time centered)
    real:: Rho, b_D(3)

    ! Time step for the cell
    real:: Dt

    ! Coefficients of linear equation a*dM + BetaB x dM = Rhs
    real:: BetaDtQPerM, a, BetaB_D(3), Rhs_D(3)  

    ! Change in momentum and electric field
    real:: dRhoU_D(3), dE_D(3)
    !------------------------------------------------------------------------
    if(BetaImpl <= 0.0) RETURN

    ! Set the Source as the explicit multi-stage RHS
    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE
       Source_VC(1:nVar,i,j,k) = &
            (State_VGB(:,i,j,k,iBlock) - StateOld_VCB(:,i,j,k,iBlock)) &
            /(Cfl*time_BLK(i,j,k,iBlock))
    end do; end do; end do

    ! Add the stiff sources at old time level
    call add_ion_electron_source(iBlock, UseOld=.true.)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       ! Time step for this cell
       Dt = Cfl*time_BLK(i,j,k,iBlock)

       ! Time centered full magnetic field
       b_D = 0.5*(State_VGB(Bx_:Bz_,i,j,k,iBlock) +StateOld_VCB(Bx_:Bz_,i,j,k,iBlock))
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       ! RHS for electric field
       dE_D = Dt*Source_VC(Ex_:Ez_,i,j,k)

       ! do point-implicit update for each fluid one-by-one (including electrons)
       do iIon = 1, nIonFluid

          ! Set indexes
          iRho   = iRhoIon_I(iIon)
          iRhoUx = iRhoUxIon_I(iIon)
          iRhoUz = iRhoUzIon_I(iIon)

          BetaDtQPerM = BetaImpl*Dt*ChargePerMass_I(iIon)

          ! time centered ion/electron density at time level n
          Rho = 0.5*(State_VGB(iRho,i,j,k,iBlock) + StateOld_VCB(iRho,i,j,k,iBlock))

          ! Rhs_D = Dt*Rhs(momentum) + Beta*Dt*Q/M*Dt*Rhs(Efield)
          Rhs_D = Dt*Source_VC(iRhoUx:iRhoUz,i,j,k) + BetaDtQPerM*Rho*dE_D

          ! Magnetic field multiplied with Beta*DtQeMe
          BetaB_D = BetaDtQPerM*b_D

          ! Scalar coefficient in the vector equation below
          a = 1 + BetaDtQPerM**2 * C2light * Rho

          ! Analytic solution of A*dRhoU + BetaB x dRhoU = Rhs 
          dRhoU_D = (a*Rhs_D - cross_product(BetaB_D, Rhs_D) &
               + BetaB_D*sum(BetaB_D*Rhs_D)/a)&
               /(a**2 + sum(BetaB_D**2))

          ! Change in electric field due to implicit electron current
          dE_D = dE_D - C2light*BetaDtQPerM*dRhoU_D


          ! Update ion/electron momentum
          State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) = &
               State_VGB(iRhoUx:iRhoUz,i,j,k,iBlock) + dRhoU_D

       end do ! Ion loop

       ! Update electric field with all the fluid contributions
       State_VGB(Ex_:Ez_,i,j,k,iBlock) = &
            State_VGB(Ex_:Ez_,i,j,k,iBlock)+ dE_D

    end do; end do; end do

    ! Update energy since momentum changed
    call calc_energy(1,nI,1,nJ,1,nK,iBlock,IonFirst_,IonLast_)
    
  end subroutine update_impl_ion_electron
  !===========================================================================
end module ModIonElectron
