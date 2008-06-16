!^CFG COPYRIGHT UM
!==============================================================================
module EEE_ModGetStateBC
  implicit none
  save

  private

  public :: EEE_get_state_BC

contains

  !============================================================================

  subroutine EEE_get_state_BC(x_D,Rho,U_D,B_D,p,Time,n_step,iteration_number)
    use EEE_ModTD99
    use EEE_ModShearFlow
    implicit none

    real, intent(in) :: x_D(3),Time
    real, intent(out) :: Rho,U_D(3),B_D(3),p
    integer, intent(in) :: n_step,iteration_number

    real :: Rho1,U1_D(3),B1_D(3),p1
    !--------------------------------------------------------------------------

    ! initialize perturbed state variables
    Rho=0.0; U_D=0.0; B_D=0.0; p=0.0

    if (DoTD99FluxRope.or.DoBqField) then
       call get_transformed_TD99fluxrope(x_D,B1_D,&
            U1_D,n_step,Iteration_Number,Rho1,Time)

       if(.not.DoBqField) U1_D=0.0

       Rho=Rho+Rho1; U_D=U_D+U1_D; B_D=B_D+B1_D
    end if

    if(UseShearFlow)then
       call get_shearflow(x_D,Time,U1_D,iteration_number)

       U_D = U_D + U1_D
    end if

  end subroutine EEE_get_state_BC

end module EEE_ModGetStateBC
