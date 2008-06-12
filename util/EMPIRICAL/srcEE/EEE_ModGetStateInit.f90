!^CFG COPYRIGHT UM
!==============================================================================
module EEE_ModGetStateInit
  implicit none
  save

  private

  public :: EEE_get_state_init

contains

  !============================================================================

  subroutine EEE_get_state_init(x_D,Rho,B_D,p,n_step,iteration_number)
    use EEE_ModGL98
    use EEE_ModTD99
    implicit none

    real, intent(in) :: x_D(3)
    real, intent(out) :: Rho,B_D(3),p
    integer, intent(in) :: n_step,iteration_number

    real :: U_D(3)
    real :: Rho1,U1_D(3),B1_D(3),p1
    !--------------------------------------------------------------------------

    ! initialize perturbed state variables
    Rho=0.0; U_D=0.0; B_D=0.0; p=0.0

    if(UseTD99Perturbation)then
       !\
       ! Add Titov & Demoulin (TD99) flux rope
       !/
       if(UseVariedCurrent)then
          Rho1=0.0; U1_D=0.0; B1_D=0.0
       else
          call get_transformed_TD99fluxrope(x_D,B1_D,&
               U1_D,n_step,iteration_number,Rho1)
       end if

       Rho=Rho+Rho1; U_D=U_D+U1_D; B_D=B_D+B1_D
    end if

    if(UseFluxRope)then
       !\
       ! Add Gibson & Low (GL98) flux rope
       !/
       call get_GL98_fluxrope(x_D,Rho1,p1,B1_D)

       call adjust_GL98_fluxrope(Rho1,p1)

       Rho=Rho+Rho1; B_D=B_D+B1_D; p=p+p1
    end if

  end subroutine EEE_get_state_init

end module EEE_ModGetStateInit
