module ModTimeStepControl

  implicit none

  private ! except

  public:: read_time_step_control_param
  public:: control_time_step
  logical, public:: UseTimeStepControl = .false.

  ! Local variables
  integer             :: nVarControl
  integer, allocatable:: iVarControl_I(:)
  real,    allocatable:: VarRatio_I(:)

  real :: &
       RejectStepLevel1   = 0.3, &
       RejectStepLevel2   = 3.0, &
       RejectStepFactor   = 0.50, &
       ReduceStepLevel1   = 0.6, &
       ReduceStepLevel2   = 1.5, &
       ReduceStepFactor   = 0.95, &
       IncreaseStepLevel1 = 0.8, &
       IncreaseStepLevel2 = 1.2, &
       IncreaseStepFactor = 1.05

contains

  !===========================================================================
  subroutine read_time_step_control_param(NameCommand)

    use ModReadParam
    use ModVarIndexes, ONLY: nVar, NameVar_V
    use ModUtilities,  ONLY: lower_case, split_string

    character(len=*), intent(in):: NameCommand

    integer :: iControl, iVar
    character(len=100):: NameVar, NameVarControl
    character(len=20) :: NameVarControl_I(nVar)

    character(len=*), parameter:: NameSub='read_time_step_control_param'
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#CONTROLTIMESTEP", "#TIMESTEPCONTROL")
       call read_var('UseTimeStepControl', UseTimeStepControl)
    case("#CONTROLDECREASE")
       call read_var('RejectStepLevel1' ,  RejectStepLevel1)
       call read_var('ReduceStepLevel1' ,  ReduceStepLevel1)
       call read_var('IncreaseStepLevel1', IncreaseStepLevel1)
    case("#CONTROLINCREASE")
       call read_var('RejectStepLevel2' ,  RejectStepLevel2)
       call read_var('ReduceStepLevel2' ,  ReduceStepLevel2)
       call read_var('IncreaseStepLevel2', IncreaseStepLevel2)
    case("#CONTROLFACTOR")
       call read_var('RejectStepFactor',   RejectStepFactor)
       call read_var('ReduceStepFactor',   ReduceStepFactor)
       call read_var('IncreaseStepFactor', IncreaseStepFactor)
    case("#CONTROLVAR")
       call read_var('NameVarControl', NameVarControl, IsLowerCase=.true.)
       ! Split list of variables and match with NameVar_V
       call split_string(NameVarControl, nVar, NameVarControl_I, nVarControl)
       if(allocated(iVarControl_I)) deallocate(iVarControl_I, VarRatio_I)
       allocate( iVarControl_I(nVarControl), VarRatio_I(nVarControl) )
       do iControl=1, nVarControl
          do iVar = 1, nVar
             NameVar = NameVar_V(iVar)
             call lower_case(NameVar)
             if( NameVar == NameVarControl_I(iControl) ) EXIT
          end do
          if(iVar > nVar)call stop_mpi(NameSub//' invalid NameVarControl='//&
               NameVarControl)
          iVarControl_I(iControl) = iVar
       end do
    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_time_step_control_param

  !===========================================================================

  subroutine control_time_step

    use ModMain,    ONLY: nBlock, nI, nJ, nK, UnusedBlk, Dt, Cfl, &
         DtFixed, DtFixedOrig, UseDtFixed, Time_Simulation
    use ModAdvance, ONLY: Rho_, p_, &
         State_VGB, StateOld_VCB, Energy_GBI, EnergyOld_CBI, time_BLK
    use ModPhysics, ONLY: No2Si_V, UnitT_
    use ModProcMH,  ONLY: iProc, nProc, iComm
    use ModMpi

    integer:: iBlock, i, j, k, iError
    real   :: RelativeChangeMin,  RelativeChangeMax, Tmp, Factor
    real   :: CflOrig = -1.0

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='control_time_step'
    !-------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Use density and pressure as control variables by default
    if(nVarControl < 0)then
       nVarControl = 2
       allocate( iVarControl_I(nVarControl), VarRatio_I(nVarControl) )
       iVarControl_I = (/Rho_, p_/)
    end if

    ! Initialize CflOrig
    if(CflOrig < 0.0) CflOrig = Cfl

    ! Calculate the largest relative drop in the control variables
    RelativeChangeMin = 1e+6
    RelativeChangeMax = 1e-6

    do iBlock = 1, nBlock
       if(UnusedBlk(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          VarRatio_I = State_VGB(iVarControl_I,i,j,k,iBlock) &
               /    StateOld_VCB(iVarControl_I,i,j,k,iBlock)
          RelativeChangeMin = min(RelativeChangeMin, minval(VarRatio_I))
          RelativeChangeMax = max(RelativeChangeMin, maxval(VarRatio_I))
       end do; end do; end do
    end do

    if(nProc > 1)then
       ! Take minimum and maximum over all the processors
       Tmp = RelativeChangeMin
       call MPI_allreduce(Tmp, RelativeChangeMin, 1, MPI_REAL, MPI_MIN, &
            iComm, iError)
       Tmp = RelativeChangeMax
       call MPI_allreduce(Tmp, RelativeChangeMax, 1, MPI_REAL, MPI_MAX, &
            iComm, iError)
    end if

    ! Figure out time step reduction factor
    if(       RelativeChangeMin < RejectStepLevel1  &
         .or. RelativeChangeMax > RejectStepLevel2 )then
       ! Redo step if change is outside [RejectStepLevel1, RejectStepLevel2]
       Time_Simulation = Time_Simulation - Dt*No2Si_V(UnitT_)
       Dt = 0.0
       ! Do not use previous step in BDF2 scheme
       !!! n_prev = -1
       ! Reset the state variable, the energy and set time_BLK variable to 0
       do iBlock = 1, nBlock
          if(UnusedBlk(iBlock)) CYCLE
          State_VGB(:,1:nI,1:nJ,1:nK,iBlock)  = StateOld_VCB(:,:,:,:,iBlock)
          Energy_GBI(1:nI,1:nJ,1:nK,iBlock,:) = EnergyOld_CBI(:,:,:,iBlock,:)
          time_BLK(1:nI,1:nJ,1:nK,iBlock)     = 0.0
       end do
       ! Reduce next time step
       Factor = RejectStepFactor

    elseif(   RelativeChangeMin < ReduceStepLevel1 &
         .or. RelativeChangeMax > ReduceStepLevel2 )then
       ! Reduce next time step if pressure is reduced below ReduceStepLevel
       Factor = ReduceStepFactor

    elseif(RelativeChangeMin    > IncreaseStepLevel1  &
         .or. RelativeChangeMax > IncreaseStepLevel2 )then
       ! Increase next time step if change remained above IncreaseStepLevel
       ! and the last step was taken with DtFixed. Do not exceed DtFixedOrig
       Factor = IncreaseStepFactor
    else
       Factor = 1.0
    end if

    if(UseDtFixed)then
       DtFixed = min(DtFixedOrig, DtFixed*Factor)
    else
       Cfl     = min(CflOrig, Cfl*Factor)
    end if

    if(DoTestMe)then
       write(*,*) NameSub,': RelativeChangeMin,Max,Factor=', &
            RelativeChangeMin, RelativeChangeMax, Factor
       write(*,*) NameSub,': Dt, DtFixed, Cfl=',&
            Dt*No2Si_V(UnitT_), DtFixed*No2Si_V(UnitT_), Cfl
    end if

  end subroutine control_time_step

end module ModTimeStepControl
