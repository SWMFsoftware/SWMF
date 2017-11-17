!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPartSteady

  use BATL_lib, ONLY: &
       test_start, test_stop, lVerbose

  ! The Partially Steady State (PSS) scheme is based on the idea that
  ! one can evolve only the grid blocks that change in time, while
  ! the grid blocks containing a steady state can be ignored.
  ! Since information can travel by at most one cell/time step in an explicit
  ! scheme, it is sufficient to switch on blocks which neighbor the
  ! already evolving blocks.

  use ModProcMH
  use ModVarIndexes, ONLY: nVar
  use ModSize,       ONLY: MaxBlock, nI, nJ, nK, x_
  use ModMain,       ONLY: iNewDecomposition, nBlock, nBlockMax, &
       time_accurate, n_step
  use ModGeometry,   ONLY: CellSize_DB, CellSize1Min
  use ModParallel,   ONLY: NOBLK, NeiLev, NeiPe, NeiBlk
  use ModAdvance,    ONLY: iTypeAdvance_B, iTypeAdvance_BP, &
       SkippedBlock_, SteadyBlock_, SteadyBoundBlock_, ExplBlock_, &
       State_VGB, StateOld_VGB
  use BATL_lib,      ONLY: Unused_B, Unused_BP

  use ModMpi

  implicit none

  save

  private ! except

  logical, public :: UsePartSteady = .false.     ! True if the scheme is used
  logical, public :: IsSteadyState = .false.     ! True for full steady state
  logical, public :: IsNewSteadySelect = .false. ! True if selection changed
  integer, public :: MinCheckVar = 1     ! First variable to check for change
  integer, public :: MaxCheckVar = nVar  ! Last  variable to check for change

  real, public :: RelativeEps_V(nVar) = 0.001  ! Relative change per variable
  real, public :: AbsoluteEps_V(nVar) = 0.0001 ! Absolute change per variable

  public part_steady_switch ! Switch to and from the evolving/bounday blocks
  public part_steady_select ! Select the evolving/bounday blocks

  logical, parameter :: DoDebug = .false.

contains
  !============================================================================
  subroutine part_steady_switch(IsOn)

    ! If IsOn=.true., set used blocks to the evolving and boundary blocks
    ! If IsOn=.false., return to the original used blocks

    logical, intent(in) :: IsOn
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'part_steady_switch'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoDebug)write(*,*)'part_steady_switch called with IsOn=',IsOn

    if(IsOn)then
       ! Make block unused if it is skipped or steady
       Unused_BP(1:nBlockMax,:) = &
            iTypeAdvance_BP(1:nBlockMax,:) <= SteadyBlock_

       ! Change the decomposition index
       iNewDecomposition=mod(iNewDecomposition+1, 10000)
    else
       ! Restore the original unused blocks
       Unused_BP(1:nBlockMax,:) = &
            iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_

       ! Restore the decomposition index
       iNewDecomposition = mod(iNewDecomposition-1, 10000)
    end if
    ! Update local Unused_B array
    Unused_B(1:nBlockMax) = Unused_BP(1:nBlockMax,iProc)

    if(DoDebug)write(*,*)'part_steady_switch nBlockMax, nUnused=',&
         nBlockMax, count(Unused_B(1:nBlockMax))

    call test_stop(NameSub, DoTest)
  end subroutine part_steady_switch
  !============================================================================

  subroutine part_steady_select

    ! Select the blocks which are evolving and in the steady boundary
    ! Set IsNewSteadySelect = .true. if there is any change

    logical :: DoPreserveExpl = .false. ! Preserve explicit blocks
    logical :: IsChanged
    integer :: i, j, k, iVar, iBlock
    integer :: jBlock, jProc, iFace, nSubFace, iSubFace
    integer :: iError

    real :: dState_V(nVar), dStateLimit
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'part_steady_select'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    if(DoDebug)write(*,*)'part_steady_select'

    ! Check the advanced blocks (ExplBlock_ and SteadyBoundBlock_) for change
    IsChanged = .false.
    do iBlock = 1, nBlock
       ! Skip steady blocks
       if(iTypeAdvance_B(iBlock) <= SteadyBlock_) CYCLE

       ! Skip explicit blocks if Expl->Steady change is not allowed
       if(DoPreserveExpl .and. iTypeAdvance_B(iBlock) == ExplBlock_) CYCLE

       ! Calculate the maximum change in the block for each variable
       dState_V = 0.0
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iVar = MinCheckVar, MaxCheckVar
             ! Normalize maximum change to relative and absolute limits
             dState_V(iVar) = max(dState_V(iVar), &
                  + abs(State_VGB(iVar,i,j,k,iBlock) &
                  -     StateOld_VGB(iVar,i,j,k,iBlock)) &
                  / (RelativeEps_V(iVar) * abs(State_VGB(iVar,i,j,k,iBlock)) &
                  +  AbsoluteEps_V(iVar)))
          end do
       end do; end do; end do

       if(time_accurate)then
          ! Take into account the cell size difference between blocks
          ! so that the same FLUX has the same effect (dU/dt=dF/dx)
          dStateLimit = CellSize1Min / CellSize_DB(x_,iBlock)
       else
          ! The local time step takes care of most of the cell size differences
          ! Ideally one should use Limit = min(dx)/dx*dtcell/min(dtcell)
          ! but that would require extra memory and communication.
          dStateLimit = 1.0
       end if

       ! Check if the change is significant and modify block type if necessary
       if(any(dState_V(MinCheckVar:MaxCheckVar) > dStateLimit))then
          if(iTypeAdvance_B(iBlock) == SteadyBoundBlock_)then
             iTypeAdvance_B(iBlock) = ExplBlock_
             IsChanged = .true.
          end if
       elseif( .not. DoPreserveExpl) then
          if(iTypeAdvance_B(iBlock) ==  ExplBlock_) then
             iTypeAdvance_B(iBlock) = SteadyBoundBlock_
             IsChanged = .true.
          end if
       end if

    end do

    call MPI_allreduce(IsChanged, IsNewSteadySelect, 1, MPI_LOGICAL, MPI_LOR, &
          iComm, iError)

    if(DoDebug)write(*,*)'iProc, IsChanged, IsNewSteadySelect=',&
         iProc, IsChanged, IsNewSteadySelect

    ! In time accurate runs the explicit blocks should not be modified
    ! to steady state, because the blocks may change slowly for long time.
    ! In steady state runs, the blocks become all steady state in the end..
    DoPreserveExpl = time_accurate

    ! If no new evolving blocks were found, simply return
    if(.not. IsNewSteadySelect) RETURN

    ! Update the global information about evolving blocks
    call MPI_allgather(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
         iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)

    ! Find the blocks surrounding the evolving blocks
    ! First set all non-explicit blocks to steady
    where(iTypeAdvance_B == SteadyBoundBlock_) iTypeAdvance_B = SteadyBlock_

    ! Check the steady blocks if they are now part of the steady boundary
    BLOCKS: do iBlock = 1, nBlock

       ! Skip all other blocks
       if(iTypeAdvance_B(iBlock) /= SteadyBlock_) CYCLE

       if(DoDebug)write(*,*)'part_steady checking iBlock=',iBlock

       ! Check all faces and subfaces
       FACES: do iFace = 1, 6

          if(NeiLev(iFace,iBlock) == NOBLK) CYCLE FACES
          if(NeiLev(iFace,iBlock) == -1)then
             nSubFace = 4
          else
             nSubFace = 1
          end if

          SUBFACES: do iSubFace = 1, nSubFace
             jProc  = NeiPE( iSubFace,iFace,iBlock)
             jBlock = NeiBLK(iSubFace,iFace,iBlock)
             if(iTypeAdvance_BP(jBlock, jProc) >= ExplBlock_) then
                iTypeAdvance_B(iBlock) = SteadyBoundBlock_
                CYCLE BLOCKS
             end if
          end do SUBFACES
       end do FACES

    end do BLOCKS

    ! Update the global information about block types
    call MPI_allgather(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
         iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)

    ! Check for full steady state
    if(.not.time_accurate) &
         IsSteadyState = all(iTypeAdvance_BP(1:nBlockMax,:) /= ExplBlock_)

    if(iProc==0 .and. lVerbose>0) &
         write(*,*)'part_steady finished:',&
         ' nStep,nSkipped,Steady,Bound,ExplALL=',n_step, &
         count(iTypeAdvance_BP == SkippedBlock_), &
         count(iTypeAdvance_BP == SteadyBlock_), &
         count(iTypeAdvance_BP == SteadyBoundBlock_), &
         count(iTypeAdvance_BP == ExplBlock_)

    call test_stop(NameSub, DoTest)
  end subroutine part_steady_select
  !============================================================================

end module ModPartSteady
!==============================================================================
