module ModIonElectron

  use ModVarIndexes
  use ModAdvance, ONLY: State_VGB, Source_VC
  use ModPhysics, ONLY: C2light
  use ModGeometry, ONLY: true_cell
  use ModB0,       ONLY: UseB0, B0_DGB
  use ModMultiFluid, ONLY: iRhoIon_I, iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I,&
       ChargePerMass_I
  use BATL_lib,    ONLY: nI, nJ, nK, x_, y_, z_
  
  implicit none

  real:: BetaImpl = 1.0

  ! Index of electron fluid
  integer, parameter:: El_ = IonLast_
  ! Index of the last TRUE (=heavy) ion
  integer, parameter:: i_ = El_ - 1

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
  subroutine add_ion_electron_source_expl(iBlock)

    integer, intent(in):: iBlock

    real:: State_V(nVar), b_D(3)

    integer:: i, j, k
    !-----------------------------------------------------------------------

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

       ! de_s/dt += j_s.E = (q/m)_s*rhou_s.E 
       Source_VC(Energy_:,i,j,k) = Source_VC(Energy_:,i,j,k) &
            + ChargePerMass_I * (State_V(iRhoUxIon_I)*State_V(Ex_) &
            +                    State_V(iRhoUyIon_I)*State_V(Ey_) &
            +                    State_V(iRhoUzIon_I)*State_V(Ez_) )

    end do; end do; end do

  end subroutine add_ion_electron_source_expl

  !===========================================================================
  subroutine add_ion_electron_source(iBlock)

    integer, intent(in):: iBlock

    real:: State_V(nVar), b_D(3)

    integer:: i, j, k
    !-----------------------------------------------------------------------

    if(BetaImpl <= 0.0) then
       call add_ion_electron_source_expl(iBlock)
       RETURN
    end if

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       State_V = State_VGB(:,i,j,k,iBlock)

       b_D = State_V(Bx_:Bz_)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       ! d(rhou)/dt += q/m*(rho*E + rhou x B)
       Source_VC(iRhoUxIon_I(1:i_),i,j,k) =          &
            Source_VC(iRhoUxIon_I(1:i_),i,j,k)       &
            + ChargePerMass_I(1:i_) *                &
            ( State_V(iRhoIon_I(1:i_))*State_V(Ex_)  &
            + State_V(iRhoUyIon_I(1:i_))*b_D(z_)     &
            - State_V(iRhoUzIon_I(1:i_))*b_D(y_) )

       Source_VC(iRhoUyIon_I(1:i_),i,j,k) =          &
            Source_VC(iRhoUyIon_I(1:i_),i,j,k)       &
            + ChargePerMass_I(1:i_) *                &
            ( State_V(iRhoIon_I(1:i_))*State_V(Ey_)  &
            + State_V(iRhoUzIon_I(1:i_))*b_D(x_)     &
            - State_V(iRhoUxIon_I(1:i_))*b_D(z_) )

       Source_VC(iRhoUzIon_I(1:i_),i,j,k) =          &
            Source_VC(iRhoUzIon_I(1:i_),i,j,k)       &
            + ChargePerMass_I(1:i_) *                &
            ( State_V(iRhoIon_I(1:i_))*State_V(Ez_)  &
            + State_V(iRhoUxIon_I(1:i_))*b_D(y_)     &
            - State_V(iRhoUyIon_I(1:i_))*b_D(x_) )

       ! electric field is updated later ?!
       Source_VC(Ex_,i,j,k) = &
            -C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUxIon_I(1:i_)))
       Source_VC(Ey_,i,j,k) = &
            -C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUyIon_I(1:i_)))
       Source_VC(Ez_,i,j,k) = &
            -C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUzIon_I(1:i_)))

       ! The energy update is explicit even for the electrons
       ! since there is no feed back from energy to E or ElRhoU
       ! de_s/dt += j_s.E = (q/m)_s*rhou_s.E 
       Source_VC(Energy_:,i,j,k) = Source_VC(Energy_:,i,j,k) &
            + ChargePerMass_I * (State_V(iRhoUxIon_I)*State_V(Ex_) &
            +                    State_V(iRhoUyIon_I)*State_V(Ey_) &
            +                    State_V(iRhoUzIon_I)*State_V(Ez_) )

    end do; end do; end do

  end subroutine add_ion_electron_source
  !===========================================================================
  subroutine add_impl_electron_source(iBlock)

    use ModCoordTransform, ONLY: cross_product
    use ModMain, ONLY: Dt

    integer, intent(in):: iBlock


    ! Solve for electron momentum and electric field sources implicitly

    integer:: i, j, k, iRhoEl, iRhoUxEl, iRhoUzEl

    real:: DtQeMe, a
    real:: Rho, RhoU_D(3), b_D(3), e_D(3), Rhs_D(3)
    real:: dRhoU_D(3), dE_D(3)
    !------------------------------------------------------------------------

    if(BetaImpl <= 0.0) RETURN

    DtQeMe = dt*ChargePerMass_I(El_)

    iRhoEl   = iRhoIon_I(El_)
    iRhoUxEl = iRhoUxIon_I(El_)
    iRhoUzEl = iRhoUzIon_I(El_)

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       ! Electron density and momentum
       Rho    = State_VGB(iRhoEl,i,j,k,iBlock)
       RhoU_D = State_VGB(iRhoUxEl:iRhoUzEl,i,j,k,iBlock)

       ! Full magnetic field
       b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       ! Electric field
       e_D = State_VGB(Ex_:Ez_,i,j,k,iBlock)

       ! Add ion current contribution explicitly
!       e_D(x_) = e_D(x_) &
!            -Dt*C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUxIon_I(1:i_)))
!       e_D(y_) = e_D(y_) &
!            -Dt*C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUyIon_I(1:i_)))
!       e_D(z_) = e_D(z_) &
!            -Dt*C2light*sum(ChargePerMass_I(1:i_) * State_V(iRhoUzIon_I(1:i_)))

       ! Rhs_D = dt*(Explicit source - extra term)
       Rhs_D = DtQeMe*( Rho*e_D + cross_product(RhoU_D, b_D)) &
            - BetaImpl*DtQeMe**2*Rho*C2light*RhoU_D

       ! Multiply the B field with Beta*DtQeMe as it only occurs together
       b_D = b_D*BetaImpl*DtQeMe

       ! Scalar coefficient in the vector equation below
       a = 1 + (BetaImpl*DtQeMe)**2 * C2light * Rho

       ! Analytic solution of A*dRhoU + B x dRhoU = Rhs 
       dRhoU_D = (a*Rhs_D - cross_product(b_D, Rhs_D) + b_D*sum(b_D*Rhs_D)/a)&
            /(a**2 + sum(b_D**2))

       ! Change in electric field due to implicit electron current
       dE_D = -C2light*DtQeMe*(BetaImpl*dRhoU_D + RhoU_D)

       ! Update electron momentum
       State_VGB(iRhoUxEl:iRhoUzEl,i,j,k,iBlock) = RhoU_D + dRhoU_D

       ! Update electric field
       State_VGB(Ex_:Ez_,i,j,k,iBlock) = e_D + dE_D

    end do; end do; end do
    
  end subroutine add_impl_electron_source
  !===========================================================================
end module ModIonElectron
