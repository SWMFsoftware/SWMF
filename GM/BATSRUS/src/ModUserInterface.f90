!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUserInterface

  ! Provide the interface information for user routines in user_interface.f90 
  ! This is required if the routine has optional arguments.
  ! The external can be called instead of calling the routines directly 
  ! from ModUser to avoid circular dependencies. This module can be used
  ! instead of ModUser to check/provide the subroutine interface.
  ! To avoid circular dependencies at the level of files, 
  ! this module cannot be in the same file as user_interface.f90.

  interface
     subroutine user_material_properties(State_V, i, j, k, iBlock, iDir, &
          EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
          EinternalOut, TeOut, PressureOut, &
          CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
          OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W)

       ! The State_V vector is in normalized units, all other physical
       ! quantities are in SI.
       !
       ! If the electron energy is used, then EinternalIn, EinternalOut,
       ! PressureOut, CvOut refer to the electron internal energies,
       ! electron pressure, and electron specific heat, respectively.
       ! Otherwise they refer to the total (electron + ion) internal energies,
       ! total (electron + ion) pressure, and the total specific heat.

       use ModAdvance,    ONLY: nWave
       use ModVarIndexes, ONLY: nVar

       real, intent(in) :: State_V(nVar)
       integer, optional, intent(in):: i, j, k, iBlock, iDir  ! cell/face index
       real, optional, intent(in)  :: EinternalIn             ! [J/m^3]
       real, optional, intent(in)  :: TeIn                    ! [K]
       real, optional, intent(out) :: NatomicOut              ! [1/m^3]
       real, optional, intent(out) :: AverageIonChargeOut     ! dimensionless
       real, optional, intent(out) :: EinternalOut            ! [J/m^3]
       real, optional, intent(out) :: TeOut                   ! [K]
       real, optional, intent(out) :: PressureOut             ! [Pa]
       real, optional, intent(out) :: CvOut                   ! [J/(K*m^3)]
       real, optional, intent(out) :: GammaOut                ! dimensionless
       real, optional, intent(out) :: HeatCondOut             ! [J/(m*K*s)]
       real, optional, intent(out) :: IonHeatCondOut          ! [J/(m*K*s)]
       real, optional, intent(out) :: TeTiRelaxOut            ! [1/s]
       real, optional, intent(out) :: &
            OpacityPlanckOut_W(nWave)                         ! [1/m]
       real, optional, intent(out) :: &
            OpacityRosselandOut_W(nWave)                      ! [1/m]

       ! Multi-group specific interface. The variables are respectively:
       !  Group Planckian spectral energy density
       real, optional, intent(out) :: PlanckOut_W(nWave)      ! [J/m^3]

     end subroutine user_material_properties
  end interface

end module ModUserInterface

