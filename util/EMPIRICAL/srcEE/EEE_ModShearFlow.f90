!^CFG COPYRIGHT UM
!==============================================================================
module EEE_ModShearFlow
  use EEE_ModCommonVariables
  implicit none
  save

  private

  public :: set_parameters_shearflow,get_shearflow

  logical, public :: UseShearFlow=.false.
  real :: Longitude,Latitude

contains

  !============================================================================

  subroutine set_parameters_shearflow
    use ModReadParam, ONLY: read_var
    implicit none
    !--------------------------------------------------------------------------
    call read_var('UseShearFlow', UseShearFlow)
    call read_var('Longitude',    Longitude)
    call read_var('Latitude',     Latitude)

  end subroutine set_parameters_shearflow

  !============================================================================

  subroutine get_shearflow(x_D,U_D)
    !\
    ! Boundary shear flow that conserves the radial magnetic flux
    !/
    use EEE_ModCommonVariables
    implicit none

    real, intent(in) :: x_D(3)
    real, intent(out) :: U_D(3)
    !--------------------------------------------------------------------------
    U_D = 0.0

    U_D = U_D*No2Si_V(UnitU_)

  end subroutine get_shearflow

end module EEE_ModShearFlow
