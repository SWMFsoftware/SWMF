!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module SP_ModMain

  use SP_ModSize, ONLY: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_
  
  use SP_ModWrite, ONLY: &
       set_write_param, write_output, NamePlotDir

  use SP_ModGrid, ONLY: &
       nVar, &
       X_, Y_, Z_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, BOld_, RhoOld_,&
       iComm, iProc, nProc, nBlock, &
       Proc_, Block_, Begin_, End_, &
       LatMin, LatMax, LonMin, LonMax, &
       RMin, RBufferMin, RBufferMax, RMax, ROrigin, &
       iGridLocal_IB, iGridGlobal_IA, iNode_II, iNode_B, State_VIB, &
       CoordMin_DI, TypeCoordSystem,&
       set_grid_param, init_grid, get_node_indexes, fix_grid_consistency
  
  use SP_ModAdvance, ONLY: &
       TimeGlobal, iIterGlobal, DoTraceShock, UseDiffusion, &
       advance, set_injection_param, init_advance_const

  implicit none

  SAVE

  private ! except
  real :: DataInputTime
  ! Methods and variables from this module 
  public:: &
       read_param, initialize, run, finalize, check,&
       TimeGlobal, iIterGlobal, DataInputTime

  ! Methods and variables from ModSize
  public:: &
       nDim, nLat, nLon, nNode, &
       iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_

  ! Methods and variables from ModGrid
  public:: &
       nVar, &
       X_, Y_, Z_, Rho_, Bx_,By_,Bz_,B_, Ux_,Uy_,Uz_, T_, RhoOld_, BOld_,&
       iComm, iProc, nProc, nBlock, &
       Proc_, Block_, Begin_, End_,&
       LatMin, LatMax, LonMin, LonMax, &
       RMin, RBufferMin, RBufferMax, RMax, ROrigin,&
       iGridLocal_IB, iGridGlobal_IA, iNode_II, iNode_B, State_VIB, &
       CoordMin_DI, TypeCoordSystem,& 
       get_node_indexes

  ! Methods and variables from ModWrite

  ! Methods and variables from ModAdvance

  !\
  ! Logicals for actions
  !----------------------------------------------------------------------------
  ! run the component
  logical:: DoRun = .true.
  ! restart the run 
  logical:: DoRestart = .false.
  ! perform initialization
  logical:: DoInit = .true.
  !/


contains

  subroutine read_param(TypeAction)
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var, read_line, read_command
    character (len=*), intent(in)     :: TypeAction ! What to do  

    ! The name of the command
    character (len=100) :: NameCommand
    character (len=*), parameter :: NameSub='SP:read_param'
    !--------------------------------------------------------------------------
    if(DoInit)then
       DoInit=.false.
    end if

    ! Read the corresponding section of input file
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          DoRestart=.true.
       case('#GRID')
          call set_grid_param
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#SAVEPLOT')
          call set_write_param
       case('#COORDSYSTEM',"#COORDINATESYSTEM")
          call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
       case('#INJECTION')
          call set_injection_param
       case('#TEST')
          ! various test modes that allow to disable certain features
          call read_var('DoTraceShock', DoTraceShock)
          call read_var('UseDiffusion', UseDiffusion)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  end subroutine read_param

  !============================================================================

  subroutine initialize(TimeStart)
    ! initialize the model
    real, intent(in):: TimeStart
    character(LEN=*),parameter:: NameSub='SP:initialize'
    !--------------------------------------------------------------------------
    iIterGlobal = 0
    TimeGlobal = TimeStart
    call init_advance_const
    call init_grid
  end subroutine initialize

  !============================================================================

  subroutine run(TimeLimit)
    ! advance the solution in time
    real, intent(in):: TimeLimit
    !------------------------------
    if(DataInputTime <= TimeGlobal)&
         RETURN
    call fix_grid_consistency
    if(DoRun) &
         ! run the model
         call advance(min(DataInputTime,TimeLimit))
    ! update time & iteration counters
    iIterGlobal = iIterGlobal + 1
    TimeGlobal = min(DataInputTime,TimeLimit)
    
    call write_output(TimeGlobal, iIterGlobal)

  end subroutine run

  !============================================================================

  subroutine check
    use ModUtilities, ONLY: make_dir
    character(LEN=*),parameter:: NameSub='SP:check'
    !--------------------------------------------------------------------------
    ! Make output and check input directories
    if(iProc==0) call make_dir(NamePlotDir)
  end subroutine check

  !============================================================================

  subroutine finalize
    character(LEN=*),parameter:: NameSub='SP:finalize'
    !--------------------------------------------------------------------------
  end subroutine finalize

end module SP_ModMain
