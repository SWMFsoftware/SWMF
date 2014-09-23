!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================================================================
module ModLocalTimeStep

  use ModSize, ONLY: nDim, nBlock, MaxBlock, nI, nJ, nK, nG, x_, y_, z_

  implicit none

  private ! except

  public:: advance_localstep  ! time accurate mode with subcycling
  logical, public :: UseLocalTimeStep = .false.

  real,    save, allocatable:: Time_B(:), TimeOld_B(:)
  integer, save, allocatable:: iStage_B(:)

  real:: DtMin, DtMax

contains
  !======================================================================
  subroutine advance_localstep(TimeSimulationLimit)

    use ModMain,       ONLY: Time_Simulation, Dt, Dt_BLK, Cfl
    use ModFaceFlux,   ONLY: calc_face_flux
    use ModFaceValue,  ONLY: calc_face_value
    use ModAdvance,    ONLY: State_VGB, Energy_GBI, &
         Flux_VX, Flux_VY, Flux_VZ
    use ModGeometry,   ONLY: Body_Blk, far_field_BCs_BLK
    use ModFaceBoundary, ONLY: set_face_boundary
    use ModCellBoundary, ONLY: set_cell_boundary
    use ModEnergy,     ONLY: calc_energy_ghost
    use ModPhysics,    ONLY: No2Si_V, UnitT_
    use ModCalcSource, ONLY: calc_source
    use ModTimeStepControl, ONLY: calc_timestep
    use ModVarIndexes, ONLY: nFluid, nVar
    use ModBlockData,  ONLY: set_block_data
    use ModCoronalHeating, ONLY: get_coronal_heat_factor, UseUnsignedFluxModel
    use ModResistivity, ONLY: set_resistivity, UseResistivity
    use BATL_lib,      ONLY: Unused_B, &
         message_pass_cell, store_face_flux, apply_flux_correction

    real,    intent(in) :: TimeSimulationLimit

    integer, parameter:: nFlux = nVar + nFluid

    real, allocatable:: Flux_VFD(:,:,:,:,:), &
         Flux_VXB(:,:,:,:,:), Flux_VYB(:,:,:,:,:), Flux_VZB(:,:,:,:,:)
    real:: TimeStage
    integer :: iStage, iBlock, nTimeStage, iStageBlock
    
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'advance_localstep'
    !-------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub,' starting'

    ! Time step (dt) for the first two steps are zeros. It is not 
    ! necessary and should be improved in the future. 
    call set_local_time_step(TimeSimulationLimit)

    if(.not. allocated(iStage_B)) then
       allocate(iStage_B(MaxBlock), Time_B(MaxBlock), TimeOld_B(MaxBlock))
       iStage_B  = 1
    endif

    if(.not. allocated(Flux_VFD)) allocate(&
         Flux_VFD(nFlux,1:nI+1,1:nJ+1,1:nK+1,nDim), &
         Flux_VXB(nFlux,nJ,nK,2,MaxBlock), & ! Two sides of the block
         Flux_VYB(nFlux,nI,nK,2,MaxBlock), &
         Flux_VZB(nFlux,nI,nJ,2,MaxBlock)) 

    Flux_VFD = 0.0
    Flux_VXB = 0.0
    Flux_VYB = 0.0
    Flux_VZB = 0.0

    Dt = DtMax
    nTimeStage = nint(DtMax/DtMin)
    TimeStage  = Time_Simulation
    Time_B    = Time_Simulation
    TimeOld_B = Time_Simulation


    ! Use 2nd order time step. Corresponding to the 'UseHalfStep' in 
    ! update_states_MHD.f90. 
    do iStage = 1, 2*nTimeStage
       call timing_start('message_pass_cell')
       call message_pass_cell(nVar, State_VGB, &
            DoSendCornerIn=.true., nProlongOrderIn=2, &
            TimeOld_B=TimeOld_B, Time_B=Time_B)
       call timing_stop('message_pass_cell')

       ! Calculate coronal heat factor
       ! This routine involves MPI_allreduce, so it cannot be done per block
       if(UseUnsignedFluxModel) call get_coronal_heat_factor

       ! Calculate resistivity
       ! This routine is usually local, but may involve message passing
       ! if DoMessagePassResistivity is true.
       ! A block-wise resistivity calculation should be implemented.
       if(UseResistivity) call set_resistivity

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE

          ! Skip block that already reached the current time of the stage
          if(Time_B(iBlock) > TimeStage*(1+1e-10)) CYCLE

          ! stage for block in the 2-stage update
          iStageBlock = iStage_B(iBlock)

          ! Update cell boundaries
          if (far_field_BCs_BLK(iBlock)) call set_cell_boundary( &
               nG, iBlock, nVar, State_VGB(:,:,:,:,iBlock))

          ! Make sure energy is set in the ghost cells
          call calc_energy_ghost(iBlock)

          call timing_start('calc_facevalues')
          call calc_face_value(.false., iBlock)
          call timing_stop('calc_facevalues')

          ! Update face boundaries
          if(body_BLK(iBlock)) &
               call set_face_boundary(iBlock, Time_B(iBlock), .false.)

          ! Calculate fluxes
          call timing_start('calc_fluxes')
          call calc_face_flux(.false., iBlock)
          call timing_stop('calc_fluxes')

          ! Calculate sources
          call calc_source(iBlock)

          ! Update state
          call timing_start('update_states')
          call update_states(iStageBlock, iBlock)
          call timing_stop('update_states')

          ! Need something like update_check in the future.

          ! Store fluxes for flux correction
          Flux_VFD(1:nFlux,1:nI+1,1:nJ,1:nK,x_) = &
               Flux_VX(1:nFlux,1:nI+1,1:nJ,1:nK)
          if(nJ>1) & ! 2D
               Flux_VFD(1:nFlux,1:nI,1:nJ+1,1:nK,y_) = &
               Flux_VY(1:nFlux,1:nI,1:nJ+1,1:nK)
          if(nK>1) & ! 3D
               Flux_VFD(1:nFlux,1:nI,1:nJ,1:nK+1,z_) = &
               Flux_VZ(1:nFlux,1:nI,1:nJ,1:nK+1)

          if(iStageBlock == 2) call store_face_flux(iBlock, nFlux, &
               Flux_VFD, Flux_VXB, Flux_VYB, Flux_VZB, &
               DtIn = Dt_BLK(iBlock)*Cfl, DoStoreCoarseFluxIn = .true.)

          ! Update block time
          TimeOld_B(iBlock) = Time_B(iBlock)
          Time_B(iBlock)    = Time_B(iBlock) &
               + Cfl*Dt_BLK(iBlock)/2*No2Si_V(UnitT_)

          ! Swap between iStage = 1 and 2 for this block
          iStage_B(iBlock) = 3 - iStageBlock

          ! Calculate time step limit it the block has reached the end
          ! of the time step
          if(abs(Time_B(iBlock) - Time_Simulation - Dt*No2Si_V(UnitT_)) &
               < 1.e-6) &
               call calc_timestep(iBlock)

          ! At this point the user has surely set all "block data"
          ! NOTE: The user has the option of calling set_block_data directly.
          call set_block_data(iBlock)

       enddo ! iBlock

       ! Apply flux correction at the end of full steps
       if(modulo(iStage, 4) == 0) call apply_flux_correction( &
            nVar, nFluid, State_VGB, Energy_GBI, &
            Flux_VXB, Flux_VYB, Flux_VZB, iStageIn=iStage/2)

       ! Update time for the stage
       TimeStage = TimeStage + DtMin/2*No2Si_V(UnitT_)

    enddo ! iStage

    Time_Simulation = Time_Simulation + Dt*No2Si_V(UnitT_)

    if(DoTestMe)write(*,*)NameSub,' finished'

  end subroutine advance_localstep

  !======================================================================
  subroutine set_local_time_step(TimeSimulationLimit)

    use ModGeometry,   ONLY: true_cell, true_BLK, MinDxValue, MaxDxValue
    use ModMain,       ONLY: Time_Simulation, Cfl, Dt_Blk
    use ModAdvance,    ONLY: time_BLK
    use BATL_lib,      ONLY: CellSize_DB, Unused_B
    use ModPhysics,    ONLY: No2Si_V, UnitT_, Si2No_V
    use ModProcMH,     ONLY: iComm
    use ModMpi

    real, intent(in) :: TimeSimulationLimit

    real   :: DtDxMinPe, DtDxMin 
    integer:: iError
    integer:: iBlock

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'set_local_time_step'
    !----------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub,' starting with TimeSimulationLimit=', &
         TimeSimulationLimit

    ! Find the smallest value of Dt_stable/Dx in the whole domain
    DtDxMinPe = minval(Dt_BLK(1:nBlock)/CellSize_DB(1,1:nBlock), &
         MASK=.not.Unused_B(1:nBlock))

    call MPI_allreduce(DtDxMinPe, DtDxMin, 1, MPI_REAL, MPI_MIN, iComm, iError)

    ! Smallest and largest time steps in the domain (includes Cfl number)
    DtMin = Cfl*MinDxValue*DtDxMin
    DtMax = Cfl*MaxDxValue*DtDxMin

    if(DoTestMe)write(*,*) NameSub,': original DtDxMin, DtMin, DtMax=', &
          DtDxMin, DtMin, DtMax

    ! Check if we reached the final time
    if(Time_Simulation + DtMax*No2Si_V(UnitT_) > TimeSimulationLimit) then
       ! For the last few time steps use uniform time steps
       DtMin = min(DtMin, &
            (TimeSimulationLimit - Time_Simulation)*Si2No_V(UnitT_) )
       DtMax = DtMin
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          Dt_BLK(iBlock) = DtMin/Cfl
       enddo
    else
       ! Set the block time step (without Cfl) proportional to refinement level
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          Dt_BLK(iBlock) = CellSize_DB(1,iBlock)*DtDxMin
       enddo
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       time_BLK(:,:,:,iBlock) = Dt_BLK(iBlock)

       ! Reset time step to zero inside body.
       if(.not.true_BLK(iBlock))then
          where(.not.true_cell(1:nI,1:nJ,1:nK,iBlock)) &
               time_BLK(:,:,:,iBlock) = 0.0
       end if
    enddo

    if(DoTestMe)write(*,*) NameSub,' finished with DtMin, DtMax=', &
         DtMin, DtMax

  end subroutine set_local_time_step

end module ModLocalTimeStep
