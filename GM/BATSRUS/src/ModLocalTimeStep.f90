!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================================================================
module ModLocalTimeStep

  use ModSize, ONLY: nDim, nBlock, MaxBlock, nI, nJ, nK, nG, x_, y_, z_
  use ModMain, ONLY: UseLocalTimeStep, UseLocalTimeStepNew, DtLimitDim
  use ModTimeStepControl, ONLY: UseMaxTimeStep
  use BATL_lib, ONLY: UseTimeLevel

  implicit none

  private ! except

  public:: advance_localstep    ! time accurate mode with subcycling
  public:: read_localstep_param ! read parameters for algorithm

  ! Local variables
  real:: DtMinSi, DtMaxSi

contains
  !======================================================================
  subroutine read_localstep_param(NameCommand, iSession)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand
    integer, intent(in):: iSession

    character(len=*), parameter:: NameSub = 'read_localstep_param'
    !-------------------------------------------------------------------
    select case(NameCommand)
    case('#LOCALTIMESTEP', '#SUBCYCLING')
       ! Check if we had it on already
       UseLocalTimeStepNew = .not.UseLocalTimeStep .and. iSession > 1

       call read_var('UseSubcycling', UseLocalTimeStep)
       if(UseLocalTimeStep)then
          call read_var('UseMaxTimeStep',   UseMaxTimeStep)
          UseTimeLevel = UseMaxTimeStep
          call read_var('DtLimitDim',       DtLimitDim)
       end if

       ! Check if the local time stepping was just switched on
       UseLocalTimeStepNew = UseLocalTimeStepNew .and. UseLocalTimeStep

    case default
       call stop_mpi(NameSub//': unknown command='//NameCommand)
    end select

  end subroutine read_localstep_param
  !===========================================================================
  subroutine advance_localstep(TimeSimulationLimit)

    use ModMain,       ONLY: Time_Simulation, Dt, Dt_BLK, Cfl, iStage, nStage,&
         nOrder, nOrderProlong, VarTest, iTest, jTest, kTest, BlkTest
    use ModFaceFlux,   ONLY: calc_face_flux
    use ModFaceValue,  ONLY: calc_face_value
    use ModAdvance,    ONLY: nFluid, nVar, State_VGB, Energy_GBI, &
         Flux_VX, Flux_VY, Flux_VZ
    use ModB0,         ONLY: set_b0_face
    use ModConserveFlux, ONLY: DoConserveFlux
    use ModGeometry,     ONLY: Body_Blk, far_field_BCs_BLK
    use ModFaceBoundary, ONLY: set_face_boundary
    use ModCellBoundary, ONLY: set_cell_boundary
    use ModEnergy,     ONLY: calc_energy_ghost
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitT_
    use ModCalcSource, ONLY: calc_source
    use ModTimeStepControl, ONLY: calc_timestep, DtMin, DtMax
    use ModBlockData,  ONLY: set_block_data
    use ModCoronalHeating, ONLY: get_coronal_heat_factor, UseUnsignedFluxModel
    use ModResistivity, ONLY: set_resistivity, UseResistivity
    use ModCoarseAxis, ONLY: UseCoarseAxis, coarsen_axis_cells
    use BATL_lib,      ONLY: Unused_B, min_tree_level,  &
         message_pass_cell, store_face_flux, apply_flux_correction

    real,    intent(in) :: TimeSimulationLimit

    ! Advance the solution by one time step. After the time step
    ! the simulation time should not exceed TimeSimulationLimit.
    ! Each block is advanced with local time steps proportional
    ! to the block level.
    ! For now a 2-stage scheme is applied per local time step,
    ! so it is a temporally 2nd order accurate scheme.

    ! Number of fluxes including the energy variables
    integer, parameter:: nFlux = nVar + nFluid
    real, allocatable:: Flux_VFD(:,:,:,:,:), &
         Flux_VXB(:,:,:,:,:), Flux_VYB(:,:,:,:,:), Flux_VZB(:,:,:,:,:)

    ! time and stage (1 or 2) info per block
    real   :: Time_B(MaxBlock), TimeOld_B(MaxBlock)
    integer:: iStage_B(MaxBlock)

    real   :: TimeStage, TimeEnd, DtSiTiny
    integer:: nTimeStage, iStageLocal, iLevelMin, iBlock


    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'advance_localstep'
    !-------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, &
         ' starting with TimeSimulationLimit=', TimeSimulationLimit

    if(UseMaxTimeStep)then
       ! DtMin and DtMax are set by ModTimeStepControl::set_global_timestep
       DtMinSi = DtMin*No2Si_V(UnitT_)*Cfl
       DtMaxSi = DtMax*No2Si_V(UnitT_)*Cfl

       if(DoTestMe)write(*,*) NameSub,' DtMinSi, DtMaxSi=', DtMinSi, DtMaxSi

    else
       ! set global and local time steps
       call set_local_time_step(TimeSimulationLimit)
       ! The total time step is the largest
       Dt = DtMaxSi*Si2No_V(UnitT_)
    end if

    if(DoConserveFlux)then
       if(.not. allocated(Flux_VFD)) allocate(&
            Flux_VFD(nFlux,nI+1,nJ+1,nK+1,nDim), &
            Flux_VXB(nFlux,nJ,nK,2,MaxBlock), & ! Two sides of the block
            Flux_VYB(nFlux,nI,nK,2,MaxBlock), &
            Flux_VZB(nFlux,nI,nJ,2,MaxBlock)) 

       ! Initialize fluxes for flux correction scheme
       Flux_VFD = 0.0
       Flux_VXB = 0.0
       Flux_VYB = 0.0
       Flux_VZB = 0.0


       if(DoTestMe)write(*,*) NameSub,' initialized conservative flux arrays'

    end if

    ! Initialize time and stage info for all blocks
    Time_B(1:nBlock)    = Time_Simulation
    TimeOld_B(1:nBlock) = Time_Simulation
    iStage_B(1:nBlock)  = 1

    ! A small fraction of the smallest time step in SI units
    DtSiTiny = 1e-6*DtMinSi

    ! Initialize time of current stage
    TimeStage = Time_Simulation

    ! Final simulation time
    TimeEnd = Time_Simulation + DtMaxSi

    ! Number of stages 
    nTimeStage = nint(DtMaxSi/DtMinSi)

    ! Loop over stages (nStage stages per local time step)
    do iStageLocal = 1, nStage*nTimeStage

       ! Number of grid levels involved in this stage
       iLevelMin = min_tree_level(iStageLocal)
       if(DoTestMe)write(*,*)'iStageLocal, iLevelMin=', iStageLocal, iLevelMin

       !write(*,*)'!!! left TimeOld, Time=', TimeOld_B(13), Time_B(13)
       !write(*,*)'!!! rite TimeOld, Time=', TimeOld_B(14), Time_B(14)
       !write(*,*)'!!! Bef left Rho( 3:6)=', State_VGB(3, 3:6,1,1,13)
       !write(*,*)'!!! Bef rite Rho(-1:2)=', State_VGB(3,-1:2,1,1,14)

       call timing_start('message_pass_cell')
       if(nOrder > 1 .and. nOrderProlong == 2)then
          call message_pass_cell(nVar, State_VGB, &
               TimeOld_B=TimeOld_B, Time_B=Time_B, iLevelMin=iLevelMin)
       else
          call message_pass_cell(nVar, State_VGB, &
               DoSendCornerIn=.true., &
               nProlongOrderIn=1, &
               DoRestrictFaceIn=.true., & 
               TimeOld_B=TimeOld_B, Time_B=Time_B, iLevelMin=iLevelMin)
       end if
       call timing_stop('message_pass_cell')

       !write(*,*)'!!! Aft left Rho( 3:6)=', State_VGB(3, 3:6,1,1,13)
       !write(*,*)'!!! Aft rite Rho(-1:2)=', State_VGB(3,-1:2,1,1,14)

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

          ! Calculate interface values for L/R states of each face
          !   and apply BCs for interface states as needed.
          call set_b0_face(iBlock)

          ! Skip block that already exceeded the current time of the stage
          if(Time_B(iBlock) > TimeStage + DtSiTiny) CYCLE

          ! stage for block in the nStage stage update
          iStage = iStage_B(iBlock)

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
          call update_states(iBlock)
          call timing_stop('update_states')

          ! Need something like update_check in the future.

          if(DoConserveFlux .and. iStage == nStage)then
             ! Store fluxes for flux correction
             Flux_VFD(1:nFlux,1:nI+1,1:nJ,1:nK,x_) = &
                  Flux_VX(1:nFlux,1:nI+1,1:nJ,1:nK)
             if(nJ>1) & ! 2D
                  Flux_VFD(1:nFlux,1:nI,1:nJ+1,1:nK,y_) = &
                  Flux_VY(1:nFlux,1:nI,1:nJ+1,1:nK)
             if(nK>1) & ! 3D
                  Flux_VFD(1:nFlux,1:nI,1:nJ,1:nK+1,z_) = &
                  Flux_VZ(1:nFlux,1:nI,1:nJ,1:nK+1)
             
             call store_face_flux(iBlock, nFlux, &
                  Flux_VFD, Flux_VXB, Flux_VYB, Flux_VZB, &
                  DtIn = Dt_BLK(iBlock)*Cfl, DoStoreCoarseFluxIn = .true.,&
                  DoReschangeOnlyIn = .not.UseMaxTimeStep)

             if(DoTestMe .and. iBlock == BlkTest) &
                  write(*,*) NameSub,' stored conservative flux'

          end if

          ! Update block time
          TimeOld_B(iBlock) = Time_B(iBlock)
          Time_B(iBlock)    = Time_B(iBlock) &
               + Cfl*Dt_BLK(iBlock)*No2Si_V(UnitT_)/nStage

          ! Swap between iStage = 1 and 2 for this block
          if(nStage==2) iStage_B(iBlock) = 3 - iStage

          ! Calculate time step limit for next time step
          ! if the block has reached the end of the time step
          if(Time_B(iBlock) >= TimeEnd - DtSiTiny) call calc_timestep(iBlock)

          ! STABILITY !!!
          !!!if(iStage == nStage)then
          !!!   DtBlkOld = Dt_Blk(iBlock)
          !!!   call calc_timestep(iBlock)
          !!!   if(iBlock==127) then
          !!!      write(*,*) &
          !!!           '!!! iBlock, DtBlkOld, max(time), min(time), DtBlk/DtMin=', &
          !!!           iBlock, DtBlkOld, maxval(time_BLK(:,:,:,iBlock)), &
          !!!           minval(time_BLK(:,:,:,iBlock)), Dt_BLK(iBlock)/DtMin
          !!!      write(*,*)'log2(Ratio)=', log(Dt_BLK(iBlock)/DtMin)/log(2.0)
          !!!      write(*,*)'floor(log2(Ratio))=', &
          !!!           floor(log(Dt_BLK(iBlock)/DtMin)/log(2.0))
          !!!      write(*,*)'2^floor(log2(Ratio))=', &
          !!!           2.0**floor(log(Dt_BLK(iBlock)/DtMin)/log(2.0))
          !!!   end if
          !!!   Dt_BLK(iBlock) = &
          !!!        DtMin*2.0**floor(log(Dt_BLK(iBlock)/DtMin)/log(2.0))
          !!!   if(iBlock==127) write(*,*)'!!! DtMin, new DtBlk=', &
          !!!        DtMin, Dt_BLK(iBlock)
          !!!
          !!!
          !!!   if(DtBlkOld < Dt_Blk(iBlock)) Dt_BLK(iBlock) = DtBlkOld
          !!!   time_BLK(:,:,:,iBlock) = Dt_BLK(iBlock)
          !!!   
          !!!   if(iBlock==127)write(*,*)'!!! iBlock, DtBlkOld, DtBlk=', &
          !!!        iBlock, DtBlkOld, Dt_BLK(iBlock)
          !!!end if

          ! At this point the user has surely set all "block data"
          ! NOTE: The user has the option of calling set_block_data directly.
          call set_block_data(iBlock)

       enddo ! iBlock

       ! Apply flux correction at the end of full steps
       if(DoConserveFlux .and. &
            (modulo(iStageLocal,2*nStage) == 0 .or. iStageLocal == nStage))then

          if(DoTestMe)write(*,*) NameSub, &
               ' before  conservative flux, test var=',&
               State_VGB(VarTest,iTest,jTest,kTest,BlkTest)

          call apply_flux_correction( &
               nVar, nFluid, State_VGB, Energy_GBI, &
               Flux_VXB, Flux_VYB, Flux_VZB, iStageIn=iStageLocal/nStage, &
               DoReschangeOnlyIn=.not.UseMaxTimeStep, DoTestIn=DoTestMe)

          if(DoTestMe)write(*,*) NameSub, &
               ' applied conservative flux, test var=',&
               State_VGB(VarTest,iTest,jTest,kTest,BlkTest)
       end if

       ! Coarsen the axis cells if requested
       if(UseCoarseAxis)call coarsen_axis_cells

       ! Update time for the stage
       TimeStage = TimeStage + DtMinSi/nStage

    enddo ! iStageLocal

    Time_Simulation = Time_Simulation + DtMaxSi

    if(DoTestMe)write(*,*)NameSub,' finished'

  end subroutine advance_localstep

  !======================================================================
  subroutine set_local_time_step(TimeSimulationLimit)

    use ModGeometry,   ONLY: true_cell, true_BLK, CellSize1Min, CellSize1Max
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
    DtMinSi = Cfl*CellSize1Min*DtDxMin*No2Si_V(UnitT_)
    DtMaxSi = Cfl*CellSize1Max*DtDxMin*No2Si_V(UnitT_)

    if(DoTestMe)write(*,*) NameSub,': original DtDxMin, DtMinSi, DtMaxSi=', &
          DtDxMin, DtMinSi, DtMaxSi

    ! Check if we reached the final time
    if(Time_Simulation + DtMaxSi > TimeSimulationLimit) then
       ! For the last few time steps use uniform time steps
       DtMinSi = min(DtMinSi, TimeSimulationLimit - Time_Simulation)
       DtMaxSi = DtMinSi
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          Dt_BLK(iBlock) = DtMinSi/Cfl * Si2No_V(UnitT_)
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

    if(DoTestMe)write(*,*) NameSub,' finished with DtMinSi, DtMaxSi=', &
         DtMinSi, DtMaxSi

  end subroutine set_local_time_step

end module ModLocalTimeStep
