module ModAbsorption
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes
  use ModPhysics
  use ModConst
  implicit none
contains
  subroutine calc_absorption(iCell, jCell, kCell, iBlock, Omega, Absorption)
    !The subroutine calculates the absorption coefficient, Absorption [m-1],
    ! at the circular frequency, omega.
    integer, intent(in)::iCell, jCell, kCell, iBlock !Cell and block number
    real, intent(in):: Omega     !Circular frequency, [rad/s]
    real, intent(out):: Absorption  !The absorption coefficient, [m-1]
    real:: Value_V(nVar)  !State vector
    reaL:: Rho  !Mass density, SI
    !---------------------
    Value_V(1:nVar) = State_VGB(1:nVar, iCell, jCell, kCell, iBlock)
    Rho = Value_V(rho_) * No2SI_V(UnitRho_)
  end subroutine calc_absorption
end module ModAbsorption
