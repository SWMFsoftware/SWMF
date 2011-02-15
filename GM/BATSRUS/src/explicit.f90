!^CFG COPYRIGHT UM
!===========================================================================

subroutine advance_expl(DoCalcTimestep, iStageMax)

  use ModMain
  use ModFaceFlux,  ONLY: calc_face_flux
  use ModFaceValue, ONLY: calc_face_value
  use ModAdvance,   ONLY: UseUpdateCheck, DoFixAxis,set_b0_face
  use ModParallel,  ONLY: neiLev
  use ModGeometry,  ONLY: Body_BLK
  use ModCovariant, ONLY: is_axial_geometry 
  use ModBlockData, ONLY: set_block_data
  use ModImplicit,  ONLY: UsePartImplicit           !^CFG IF IMPLICIT
  use ModPhysics,   ONLY: No2Si_V, UnitT_
  use ModConserveFlux, ONLY: save_cons_flux, apply_cons_flux, &
       nCorrectedFaceValues, CorrectedFlux_VXB, &
       CorrectedFlux_VYB, CorrectedFlux_VZB
  use ModCoronalHeating, ONLY: get_coronal_heat_factor, UseUnsignedFluxModel

  use BATL_lib, ONLY: message_pass_face
  implicit none

  logical, intent(in) :: DoCalcTimestep
  integer, intent(in) :: iStageMax ! advance only part way
  integer :: iStage, iBlock

  character (len=*), parameter :: NameSub = 'advance_expl'

  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  !\
  ! Perform multi-stage update of solution for this time (iteration) step
  !/
  if(UsePartImplicit)call timing_start('advance_expl') !^CFG IF IMPLICIT

  if(UseBody2 .and. UseOrbit) call update_secondbody  !^CFG IF SECONDBODY

  STAGELOOP: do iStage = 1, nStage

     if(DoTestMe) write(*,*)NameSub,' starting stage=',iStage

     ! If heating model depends on global B topology, update here
     if(UseUnsignedFluxModel) call get_coronal_heat_factor

     call barrier_mpi2('expl1')

     do globalBLK = 1, nBlock
        if (unusedBLK(globalBLK)) CYCLE
        if(all(neiLev(:,globalBLK)/=1)) CYCLE
        ! Calculate interface values for L/R states of each 
        ! fine grid cell face at block edges with resolution changes
        !   and apply BCs for interface states as needed.
        call set_b0_face(globalBLK)
        call timing_start('calc_face_bfo')
        call calc_face_value(.true.,GlobalBlk)
        call timing_stop('calc_face_bfo')

        if(body_BLK(globalBLK))call set_BCs(Time_Simulation, .true.)

        ! Compute interface fluxes for each fine grid cell face at
        ! block edges with resolution changes.

        call timing_start('calc_fluxes_bfo')
        call calc_face_flux(.true., GlobalBlk)
        call timing_stop('calc_fluxes_bfo')

        ! Save conservative flux correction for this solution block
        call save_cons_flux(GlobalBlk)
        
     end do

     if(DoTestMe)write(*,*)NameSub,' done res change only'

     ! Message pass conservative flux corrections.
     call timing_start('send_cons_flux')
     if(UseBatl)then
        call message_pass_face(nCorrectedFaceValues, CorrectedFlux_VXB, &
             CorrectedFlux_VYB, CorrectedFlux_VZB, DoSubtractIn=.false.)
     else
        call message_pass_faces_9conserve
     end if
     call timing_stop('send_cons_flux')

     if(DoTestMe)write(*,*)NameSub,' done message pass'

     ! Multi-block solution update.
     do globalBLK = 1, nBlock

        if (unusedBLK(globalBLK)) CYCLE

        ! Calculate interface values for L/R states of each face
        !   and apply BCs for interface states as needed.
        call set_b0_face(globalBLK)
        call timing_start('calc_facevalues')
        call calc_face_value(.false., GlobalBlk)
        call timing_stop('calc_facevalues')

        if(body_BLK(globalBLK))call set_BCs(Time_Simulation, .false.)

        ! Compute interface fluxes for each cell.
        call timing_start('calc_fluxes')
        call calc_face_flux(.false., GlobalBlk)
        call timing_stop('calc_fluxes')

        ! Enforce flux conservation by applying corrected fluxes
        ! to each coarse grid cell face at block edges with 
        ! resolution changes.
        call apply_cons_flux(GlobalBlk)

        ! Compute source terms for each cell.
        call timing_start('calc_sources')
        call calc_sources
        call timing_stop('calc_sources')

        ! Calculate time step (both local and global
        ! for the block) used in multi-stage update
        ! for steady state calculations.
        if (.not.time_accurate .and. iStage == 1 &
             .and. DoCalcTimestep) call calc_timestep

        ! Update solution state in each cell.
        call timing_start('update_states')
        call update_states(iStage,globalBLK)
        call timing_stop('update_states')

!!!        if(iStage == nStage)call calc_electric_field(globalBLK)

        if(UseConstrainB .and. iStage==nStage)then    !^CFG IF CONSTRAINB BEGIN
           call timing_start('constrain_B')
           call get_VxB(GlobalBlk)
           call bound_VxB(GlobalBlk)
           call timing_stop('constrain_B')
        end if                                        !^CFG END CONSTRAINB

        ! Calculate time step (both local and global
        ! for the block) used in multi-stage update
        ! for time accurate calculations.
        if (time_accurate .and. iStage == nStage .and. DoCalcTimestep) &
             call calc_timestep

     end do ! Multi-block solution update loop.

     if(DoTestMe)write(*,*)NameSub,' done update blocks'

     call barrier_mpi2('expl2')

     if(is_axial_geometry() .and. DoFixAxis) call fix_axis_cells

     ! Check for allowable update percentage change.
     if(UseUpdateCheck)then
        call timing_start('update_check')
        call update_check(iStage)
        call timing_stop('update_check')

        if(DoTestMe)write(*,*)NameSub,' done update check'
     end if

     if(UseConstrainB .and. iStage==nStage)then    !^CFG IF CONSTRAINB BEGIN
        call timing_start('constrain_B')
        ! Correct for consistency at resolution changes
        call correct_VxB

        ! Update face centered and cell centered magnetic fields
        do iBlock = 1, nBlock
           if(unusedBLK(iBlock))CYCLE
           call constrain_B(iBlock)
           call Bface2Bcenter(iBlock)
        end do
        call timing_stop('constrain_B')
        if(DoTestMe)write(*,*)NameSub,' done constrain B'
     end if                                        !^CFG END CONSTRAINB

     if(DoCalcTimeStep) &
          Time_Simulation = Time_Simulation + Dt*No2Si_V(UnitT_)/nStage

     if(iStage<nStage)call exchange_messages

     if(DoTestMe)write(*,*)NameSub,' finished stage=',istage

     ! exit if exceeded requested maximum stage limit
     if ((iStageMax >= 0) .and. (iStage >= iStageMax)) EXIT STAGELOOP

  end do STAGELOOP  ! Multi-stage solution update loop.

  do iBlock = 1, nBlock
     if(.not.UnusedBlk(iBlock)) call set_block_data(iBlock)
  end do

  if(UsePartImplicit)call timing_stop('advance_expl') !^CFG IF IMPLICIT

  if(DoTestMe)write(*,*)NameSub,' finished'

end subroutine advance_expl

!^CFG IF SECONDBODY BEGIN
!===========================================================================

subroutine update_secondbody
  use ModMain,     ONLY: time_simulation,globalBLK,nBlock
  use ModConst,    ONLY: cTwoPi
  use ModPhysics,  ONLY: xBody2,yBody2,OrbitPeriod,PhaseBody2,DistanceBody2

  implicit none
  !-------------------------------------------------------------------------

  ! Update second body coordinates
  xBody2 = DistanceBody2*cos(cTwoPi*Time_Simulation/OrbitPeriod+PhaseBody2)
  yBody2 = DistanceBody2*sin(cTwoPi*Time_Simulation/OrbitPeriod+PhaseBody2)

  do globalBLK = 1, nBlock
     call set_boundary_cells(globalBLK)
     call fix_block_geometry(globalBLK)
  end do
  
  call set_body_flag
  call exchange_messages
  
end subroutine update_secondbody
!^CFG END SECONDBODY


!ADJOINT SPECIFIC BEGIN
!===========================================================================

subroutine advance_expl_adjoint(DoCalcTimestep)

  use ModMain
  use ModFaceFlux,  ONLY: calc_face_flux, calc_face_flux_adjoint
  use ModFaceValue, ONLY: calc_face_value, calc_face_value_adjoint
  use ModAdvance,   ONLY: UseUpdateCheck, DoFixAxis,set_b0_face, &
       State_VGB, Energy_GBI, StateOld_VCB, EnergyOld_CBI
  use ModGeometry,  ONLY: Body_BLK
  use ModCovariant, ONLY: is_axial_geometry 
  use ModBlockData, ONLY: set_block_data
  use ModPhysics,   ONLY: No2Si_V, UnitT_
  use ModConserveFlux, ONLY: save_cons_flux, apply_cons_flux, &
       nCorrectedFaceValues, CorrectedFlux_VXB, &
       CorrectedFlux_VYB, CorrectedFlux_VZB

  use BATL_lib, ONLY: message_pass_face
  use ModAdjoint

  implicit none

  logical, intent(in) :: DoCalcTimestep
  integer :: iStage, iBlock

  character (len=*), parameter :: NameSub = 'advance_expl_adjoint'

  real :: TimeOrig
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  !\
  ! Multi-stage update of adjoint for one explicit iteration
  !/

  ! Zero out adjoint ghost values
  AdjointPrev_VGB(1:nVar,1:nI,1:nJ,1:nK,:) = Adjoint_VGB(1:nVar,1:nI,1:nJ,1:nK,:)
  AdjEnergyPrev_GBI(1:nI,1:nJ,1:nK,:,:) = AdjEnergy_GBI(1:nI,1:nJ,1:nK,:,:)
  Adjoint_VGB = 0.0
  AdjEnergy_GBI = 0.0
  Adjoint_VGB(1:nVar,1:nI,1:nJ,1:nK,:) = AdjointPrev_VGB(1:nVar,1:nI,1:nJ,1:nK,:)
  AdjEnergy_GBI(1:nI,1:nJ,1:nK,:,:) = AdjEnergyPrev_GBI(1:nI,1:nJ,1:nK,:,:)
  

  ! Zero out AdjointPrev
  AdjointPrev_VGB   = 0.0
  AdjEnergyPrev_GBI = 0.0

  ! Should be time at start of forward iteration
  TimeOrig = time_simulation

  STAGELOOP: do iStage = nStage, 1, -1

     if(DoTestMe) write(*,*)NameSub,' starting adjoint stage=',iStage

     call barrier_mpi2('expl1_adjoint')

     ! advance forward state
     if (nStage > 1)then
        ! recall state from Old storage if already advanced before
        if (iStage /= nStage)then
           State_VGB(1:nVar,1:nI,1:nJ,1:nK,:) = & 
                StateOld_VCB(1:nVar,1:nI,1:nJ,1:nK,:)
           EnergyOld_CBI(1:nI,1:nJ,1:nK,:,:) = &
                Energy_GBI(1:nI,1:nJ,1:nK,:,:)
           call exchange_messages ! Old did not store ghost values
        end if

        ! advance state to iStage-1
        if (iStage > 1) call advance_expl(.true., iStage-1)
     end if

     ! Do not yet support this for adjoint
     if(is_axial_geometry() .and. DoFixAxis) &
          call stop_mpi(NameSub // ' Not yet supported')

     ! amr/res-change flux corrections not supported yet

     ! Multi-block adjoint solution update
     do globalBLK = 1, nBlock

        if (unusedBLK(globalBLK)) CYCLE

        ! Forward solution calculations to fill face values, fluxes, etc.
        call set_b0_face(globalBLK)
        call calc_face_value(.false., GlobalBlk)
        if(body_BLK(globalBLK))call set_BCs(Time_Simulation, .false.)
        call calc_face_flux(.false., GlobalBlk)
        call apply_cons_flux(GlobalBlk)
        call calc_sources
        call update_states(iStage,globalBLK) 

        ! Update adjoint solution state in each cell.
        ! AdjSource_VC initialized/set here
        ! AdjFlux_V* initialized/set here
        ! AdjointPrev (answer) is updated here
        ! State on block is recalled to stage-beginning value
        call timing_start('update_states_adjoint')
        call update_states_adjoint(iStage,globalBLK)
        call timing_stop('update_states_adjoint')

        ! Compute adjoint source terms for each cell.
        ! AdjuDotArea_*I initialized/set here
        call timing_start('calc_sources_adjoint')
        call calc_sources_adjoint
        call timing_stop('calc_sources_adjoint')

        ! Compute adjoint interface fluxes for each cell.
        call timing_start('calc_fluxes_adjoint')
        call calc_face_flux_adjoint(.false., GlobalBlk)
        call timing_stop('calc_fluxes_adjoint')

        ! body BCs not yet supported for adjoint
        if(body_BLK(globalBLK))call stop_mpi(NameSub // ' Not yet supported')

        ! TODO: Calculate interface values for L/R states of each face
        call timing_start('calc_facevalues_adjoint')
        call calc_face_value_adjoint(.false., GlobalBlk)
        call timing_stop('calc_facevalues_adjoint')

        ! apply BCs for interface states as needed -- B0 not supported for adjoint
        !call set_b0_face_adjoint(globalBLK)

     end do ! Multi-block solution update loop.

     if(DoTestMe)write(*,*)NameSub,' done update blocks adjoint'

     call barrier_mpi2('expl2_adjoint')

     ! TODO: Take care of ghosts in Adjoint_VGB, Energy_GBI
     ! should be additive, then zero out ghost values
     call exchange_messages_adjoint

     if(DoTestMe)write(*,*)NameSub,' finished adjoint stage=',istage

  end do STAGELOOP  ! Multi-stage solution update loop.

  ! add Prev to current adjoint, result is the final answer
  Adjoint_VGB(1:nVar,1:nI,1:nJ,1:nK,:) = &
       Adjoint_VGB(1:nVar,1:nI,1:nJ,1:nK,:) + AdjointPrev_VGB(1:nVar,1:nI,1:nJ,1:nK,:) 


  if(DoTestMe)write(*,*)NameSub,' finished'

end subroutine advance_expl_adjoint

!ADJOINT SPECIFIC END
