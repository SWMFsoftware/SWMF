module ModPartSteady

  ! The Partially Steady State (PSS) scheme is based on the idea that
  ! one can evolve only the grid blocks that change in time, while
  ! the grid blocks containing a steady state can be ignored.
  ! Since information can travel by at most one cell/time step in an explicit
  ! scheme, it is sufficient to switch on blocks which neighbor the
  ! already evolving blocks. 

  use ModVarIndexes, ONLY: nVar
  use ModSize,       ONLY: MaxBlock
  use ModProcMH,     ONLY: nProc
  use ModMain,       ONLY: nBlock
  use ModAdvance,    ONLY: iTypeAdvance_B, iTypeAdvance_BP, &
       SkippedBlock_, SteadyBlock_, SteadyBoundBlock_, ExplBlock_

  implicit none

  save

  private ! except

  logical, public :: UsePartSteady = .false.     ! True if the scheme is used
  logical, public :: IsNewSteadySelect = .false. ! True if selection changed
  integer, public :: MinCheckVar = 1     ! First variable to check for change
  integer, public :: MaxCheckVar = nVar  ! Last  variable to check for change

  real, public :: RelativeEps_V(nVar) = 0.001  ! Relative change per variable
  real, public :: AbsoluteEps_V(nVar) = 0.0001 ! Absolute change per variable

  public part_steady_switch ! Switch to and from the evolving/bounday blocks
  public part_steady_select ! Select the evolving/bounday blocks

  logical, parameter :: DoDebug = .false.

contains
  !===========================================================================
  subroutine part_steady_switch(IsOn)

    use ModMain, ONLY: iNewDecomposition, nBlockMax, UnusedBLK
    use ModAMR,  ONLY: UnusedBlock_BP

    ! If IsOn=.true., set used blocks to the evolving and boundary blocks
    ! If IsOn=.false., return to the original used blocks

    logical, intent(in) :: IsOn
    !------------------------------------------------------------------------

    if(DoDebug)write(*,*)'part_steady_switch called with IsOn=',IsOn

    if(IsOn)then
       ! Store unused blocks
       where(UnusedBLK(1:nBlockMax)) &
            iTypeAdvance_B(1:nBlockMax) = SkippedBlock_ 

       where(UnusedBlock_BP(1:nBlockMax,:)) &
            iTypeAdvance_BP(1:nBlockMax,:) = SkippedBlock_

       ! Select UnusedBLK to be skipped or steady
       UnusedBLK(1:nBlockMax) = &
            iTypeAdvance_B(1:nBlockMax)    <= SteadyBlock_
       UnusedBlock_BP(1:nBlockMax,:) = &
            iTypeAdvance_BP(1:nBlockMax,:) <= SteadyBlock_

       ! Change the decomposition index
       iNewDecomposition=mod(iNewDecomposition+1, 10000)

    else
       ! Restore the original unused blocks
       UnusedBLK(1:nBlockMax)        = &
            iTypeAdvance_B(1:nBlockMax)    == SkippedBlock_
       UnusedBlock_BP(1:nBlockMax,:) = &
            iTypeAdvance_BP(1:nBlockMax,:) == SkippedBlock_

       ! Restore the decomposition index
       iNewDecomposition = mod(iNewDecomposition-1, 10000)
    end if

    if(DoDebug)write(*,*)'part_steady_switch nBlockMax, nUnused=',&
         nBlockMax, count(UnusedBLK(1:nBlockMax))

  end subroutine part_steady_switch
  !===========================================================================

  subroutine part_steady_select

    ! Select the blocks which are evolving and in the steady boundary
 !!! Set IsNewSteady = .true. if there is any change, .false. otherwise ???

    use ModMain,     ONLY: East_, Top_, n_step, lVerbose
    use ModProcMH,   ONLY: iProc, iComm
    use ModAdvance,  ONLY: State_VGB, StateOld_VCB, nI, nJ, nK, nIJK
    use ModParallel, ONLY: NOBLK, NeiLev, NeiPe, NeiBlk
    use ModMpi

    logical :: DoPreserveExpl = .false. ! Preserve explicit blocks
    logical :: IsChanged
    integer :: i, j, k, iVar, iBlock
    integer :: jBlock, jProc, iFace, nSubFace, iSubFace
    integer :: iError

    real :: dState_V(nVar), Norm_V(nVar)
    !--------------------------------------------------------------------------

    if(DoDebug)write(*,*)'part_steady_select'

!!! Possibly check cell by cell ???

    ! Check the advanced blocks (ExplBlock_ and SteadyBoundBlock_) for change
    IsChanged = .false.
    do iBlock = 1, nBlock
       ! Skip steady blocks
       if(iTypeAdvance_B(iBlock) <= SteadyBlock_) CYCLE

       ! Skip explicit blocks if Expl->Steady change is not allowed
       if(DoPreserveExpl .and. iTypeAdvance_B(iBlock) == ExplBlock_) CYCLE
       
       ! Calculate the change in the state
       dState_V = 0.0
       Norm_V  = 0.0
       do k=1,nK; do j=1,nJ; do i=1,nI
          do iVar = MinCheckVar, MaxCheckVar
             Norm_V(iVar)   = Norm_V(iVar) + &
                  abs(StateOld_VCB(iVar,i,j,k,iBlock))
             dState_V(iVar) = dState_V(iVar) &
                  + abs(State_VGB(iVar,i,j,k,iBlock) &
                  -  StateOld_VCB(iVar,i,j,k,iBlock))
          end do
       end do; end do; end do

!       if(DoDebug)write(*,*)'iBlock, dState_V, Norm_V=',&
!            iBlock, dState_V, Norm_V

       ! Normalize change in all variables
       dState_V = dState_V / (RelativeEps_V * Norm_V + nIJK*AbsoluteEps_V)

       ! Check if the change is significant and modify block type if necessary
       if(any(dState_V > 1.0))then
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

!!! For now preserve explicit blocks after the first selection
    DoPreserveExpl = .true.

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
       FACES: do iFace = East_, Top_

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

    if(iProc==0 .and. lVerbose>0) &
         write(*,*)'part_steady finished:',&
         ' nStep,nSkipped,Steady,Bound,ExplALL=',n_step, &
         count(iTypeAdvance_BP == SkippedBlock_), & 
         count(iTypeAdvance_BP == SteadyBlock_), & 
         count(iTypeAdvance_BP == SteadyBoundBlock_), &
         count(iTypeAdvance_BP == ExplBlock_)

  end subroutine part_steady_select

end module ModPartSteady
