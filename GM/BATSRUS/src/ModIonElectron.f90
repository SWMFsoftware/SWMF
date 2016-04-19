module ModIonElectron

  use ModVarIndexes
  use ModAdvance, ONLY: State_VGB, Source_VC
  use ModGeometry, ONLY: true_cell
  use ModB0,       ONLY: UseB0, B0_DGB
  use ModMultiFluid, ONLY: iRhoIon_I, iRhoUxIon_I, iRhoUyIon_I, iRhoUzIon_I,&
       ChargePerMass_I
  use BATL_lib,    ONLY: nI, nJ, nK, x_, y_, z_

  implicit none

contains
  subroutine add_ion_electron_source(iBlock)

    integer, intent(in):: iBlock

    real:: State_V(nVar), b_D(3)

    integer:: i, j, k
    !-----------------------------------------------------------------------


    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not.true_cell(i,j,k,iBlock)) CYCLE

       State_V = State_VGB(:,i,j,k,iBlock)

       b_D = State_V(Bx_:Bz_)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

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

    end do; end do; end do

  end subroutine add_ion_electron_source

end module ModIonElectron
