module ModAbsorption
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes
  use ModPhysics
  use ModConst
  implicit none
contains
  subroutine calc_absorption(State_V, Omega, Absorption)
    use CRASH_ModTransport
    use ModUser,ONLY: user_material_properties
    !The subroutine calculates the absorption coefficient, Absorption [m-1],
    ! at the circular frequency, omega.
    real,intent(in):: State_V(nVar)
    real, intent(in):: Omega     !Circular frequency, [rad/s]
    real, intent(out):: Absorption  !The absorption coefficient, [m-1]
    real:: NAtomicSI, ZAverage, TeSI, TeEV
    reaL:: Rho  !Mass density, SI
    !---------------------
    
    call user_material_properties(State_V=State_V, NAtomicOut=NAtomicSI, TeOut=TeSI,&
         AverageIonChargeOut= ZAverage)
  
  end subroutine calc_absorption
end module ModAbsorption
