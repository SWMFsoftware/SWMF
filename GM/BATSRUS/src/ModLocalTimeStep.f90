!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!  This code is a copyright protected software (c) 2002- University of Michigan
!===========================================================================
module ModLocalTimeStep
  use ModMain
  use ModFaceFlux,   ONLY: calc_face_flux
  use ModFaceValue,  ONLY: calc_face_value
  use ModAdvance,    ONLY: State_VGB,Energy_GBI, time_BLK,&
       Flux_VX, Flux_VY, Flux_VZ
  use ModGeometry,   ONLY: true_cell, true_BLK
  use ModPhysics,    ONLY: No2Si_V, UnitT_, Si2No_V
  use ModCalcSource, ONLY: calc_source
  use ModTimeStepControl, ONLY: calc_timestep
  use BATL_lib, ONLY: message_pass_cell, &
       store_face_flux, apply_flux_correction, &
       CellSize_DB
  use ModMpi
  use ModSize
  use ModProcMH

  implicit none

  logical:: UseLocalTimeStep = .false.
  real:: DtMin, DtMax
  real, allocatable:: iStage_B(:), Time_B(:), TimeOld_B(:)
contains
  !======================================================================
  subroutine advance_localstep(TimeSimulationLimit)
    implicit none
    real,    intent(in) :: TimeSimulationLimit
    integer :: iStage, iBlock, nTimeStage, iStageBlock, nFlux=nVar+nFluid
    character (len=*), parameter :: NameSub = 'advance_localstep'
    real:: TimeStage

    real, allocatable:: Flux_VFD(:,:,:,:,:), Flux_VXB(:,:,:,:,:), &
         Flux_VYB(:,:,:,:,:), Flux_VZB(:,:,:,:,:)

    integer:: i,j,k
    !-------------------------------------------------------------------------

    if(iteration_number == 1) then
       ! Calculate velocity for time step calculation. 
       call advance_expl(.true.,-1) 
       RETURN
    endif

    ! Time step (dt) for the first two steps are zeros. It is not 
    ! necessary and should be improved in the future. 
    call set_local_time_step(TimeSimulationLimit)

    if(.not. allocated(iStage_B)) then
       allocate(iStage_B(MaxBlock), Time_B(MaxBlock), &
            TimeOld_B(MaxBlock))
       iStage_B  = 1
    endif

    if(.not. allocated(Flux_VFD)) allocate(&
         Flux_VFD(nFlux,1:nI+1,1:nJ+1,1:nK+1,nDim), &
         Flux_VXB(nFlux,nJ,nK,2,MaxBlock), & ! Two layer ghost cells.
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

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          if(Time_B(iBlock) > TimeStage*(1+1e-10)) CYCLE
          iStageBlock = iStage_B(iBlock)

          call timing_start('calc_facevalues')
          call calc_face_value(.false., iBlock)
          call timing_stop('calc_facevalues')

          call timing_start('calc_fluxes')
          call calc_face_flux(.false., iBlock)
          call timing_stop('calc_fluxes')

          ! calc_source may not work now. 
          call calc_source(iBlock)

          call timing_start('update_states')
          call update_states(iStageBlock, iBlock)
          call timing_stop('update_states')

          ! Need something like update_check in the future.

          Flux_VFD(1:nFlux,1:nI+1,1:nJ,1:nK,x_) = Flux_VX(1:nFlux,1:nI+1,1:nJ,1:nK)
          if(nJ>1) then ! 2D
             Flux_VFD(1:nFlux,1:nI,1:nJ+1,1:nK,y_) = Flux_VY(1:nFlux,1:nI,1:nJ+1,1:nK)
             if(nK>1) Flux_VFD(1:nFlux,1:nI,1:nJ,1:nK+1,z_) = Flux_VZ(1:nFlux,1:nI,1:nJ,1:nK+1)              ! 3D
          endif

          if(iStageBlock == 2)  then 
             call store_face_flux(iBlock, nFlux, &
                  Flux_VFD, Flux_VXB, Flux_VYB, Flux_VZB, &
                  DtIn = Dt_BLK(iBlock)*Cfl, DoStoreCoarseFluxIn = .true.)
          endif

          TimeOld_B(iBlock) = Time_B(iBlock)
          Time_B(iBlock)    = Time_B(iBlock) + Cfl*Dt_BLK(iBlock)/2*No2Si_V(UnitT_)
          iStage_B(iBlock) = 3 - iStageBlock

          if(abs(Time_B(iBlock) - Time_Simulation - Dt*Cfl*No2Si_V(UnitT_)) < 1.e-6) &
               call calc_timestep(iBlock)
       enddo ! iBlock

       if(modulo(iStage,4) == 0) then
          call apply_flux_correction(nVar, nFluid, State_VGB, Energy_GBI, &
          Flux_VXB, Flux_VYB, Flux_VZB, iStageIn = iStage/2)
       endif
       TimeStage = TimeStage + Cfl*DtMin/2*No2Si_V(UnitT_)

    enddo ! iStage

    Time_Simulation = Time_Simulation + Cfl*Dt*No2Si_V(UnitT_)

  end subroutine advance_localstep

  !======================================================================
  subroutine set_local_time_step(TimeSimulationLimit)
    real, intent(in) :: TimeSimulationLimit
    real             :: DtMinPe, DtMaxPe, MaxVelocityPE, MaxVelocity 
    integer          :: iError
    integer          :: nTimeStage, iBlock
    !----------------------------------------------------------------------

    ! MaxVelocityPE may not exactly the physical max velocity. It is 
    ! just used in this subroutine to make sure the time step is proportional 
    ! to cell size. 
    MaxVelocityPE = maxval(CellSize_DB(1,1:nBlock)/Dt_BLK(1:nBlock), &
         MASK=.not.Unused_B(1:nBlock))
    call MPI_allreduce(MaxVelocityPE, MaxVelocity, 1, &
         MPI_REAL, MPI_MAX, iComm, iError)

    do iBlock = 1, nBlock
       Dt_BLK(iBlock) = CellSize_DB(1,iBlock)/MaxVelocity
    enddo

    DtMinPE = minval(Dt_BLK(1:nBlock), MASK=.not.Unused_B(1:nBlock))
    DtMaxPE = maxval(Dt_BLK(1:nBlock), MASK=.not.Unused_B(1:nBlock))

    call MPI_allreduce(DtMinPE, DtMin, 1, MPI_REAL, MPI_MIN, iComm, iError)
    call MPI_allreduce(DtMaxPE, DtMax, 1, MPI_REAL, MPI_MAX, iComm, iError)


    DtMax = floor( (DtMax + 1.e-10)/DtMin ) *DtMin
    Dt_BLK(1:nBlock) = floor( (Dt_BLK(1:nBlock) + 1.e-10)/DtMin ) *DtMin

    ! For the last few time steps, use the same time step.
    if(Time_Simulation+Cfl*DtMax*No2Si_V(UnitT_)*(1.+1.e-6) > &
         TimeSimulationLimit) Then
       if(Time_Simulation+Cfl*DtMin*No2Si_V(UnitT_)*(1.+1.e-6) > &
            TimeSimulationLimit) Then
          DtMin = (1+1.e-6)*(TimeSimulationLimit-Time_Simulation)*Si2No_V(UnitT_)/cfl
          Dt_BLK = DtMin
          DtMax = DtMin
       else
          Dt_BLK = DtMin
          DtMax = DtMin
       endif
    endif

    do iBlock = 1, nBlock
       time_BLK(:,:,:,iBlock) = Dt_BLK(iBlock)

       ! Reset time step to zero inside body.
       if(.not.true_BLK(iBlock))then
          where(.not.true_cell(1:nI,1:nJ,1:nK,iBlock)) &
               time_BLK(:,:,:,iBlock) = 0.0
       end if
    enddo

  end subroutine set_local_time_step

end module ModLocalTimeStep
