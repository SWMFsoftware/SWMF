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

  logical, public :: UsePartSteady       ! True if the part steady scheme used
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

  subroutine part_steady_select(IsNew)

    ! Select the blocks which are evolving and in the steady boundary
    ! Return IsNew = .true. if there is any change, .false. otherwise.

    use ModMain,     ONLY: East_, Top_, Time_Accurate, Time_Simulation, &
         lVerbose
    use ModProcMH,   ONLY: iProc, iComm
    use ModAdvance,  ONLY: State_VGB, StateOld_VCB, nI, nJ, nK, nIJK
    use ModParallel, ONLY: NOBLK, NeiLev, NeiPe, NeiBlk
    use ModMpi

    logical, intent(out) :: IsNew

    logical :: IsFirstCall = .true., IsFirstSelect=.true.
    integer :: nEvolving, nEvolvingALL, nEvolvingLastALL = -1
    integer :: i, j, k, iVar, iBlock
    integer :: jBlock, jProc, iFace, nSubFace, iSubFace
    integer :: iError

    real :: dState_V(nVar), Norm_V(nVar)
    !--------------------------------------------------------------------------

    if(DoDebug)write(*,*)'part_steady_select'

    if(IsFirstCall)then
       IsFirstCall = .false.

       ! Assign all blocks to be in the steady boundary
       where(iTypeAdvance_BP /= SkippedBlock_) &
            iTypeAdvance_BP = SteadyBoundBlock_
       iTypeAdvance_B = iTypeAdvance_BP(:,iProc)

       ! There is nothing else to do
       RETURN
    end if

    ! The first time accurate time step with zero progress cannot be used
    if(Time_Accurate .and. Time_Simulation == 0) RETURN

    ! Check the steady boundary blocks for change
    do iBlock = 1, nBlock

       if(iTypeAdvance_B(iBlock) /= SteadyBoundBlock_) CYCLE

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

       ! Check if the change is significant and modify block type
       if(any(dState_V > 1.0)) &
            iTypeAdvance_B(iBlock) = ExplBlock_

    end do

    ! Count the number of evolving blocks
    nEvolving = count(iTypeAdvance_B(1:nBlock) == ExplBlock_)

    call MPI_allreduce(nEvolving, nEvolvingALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

    ! If no new evolving blocks were found, simply return
    IsNew = nEvolvingALL /= nEvolvingLastALL

    if(DoDebug)write(*,*)'nEvolvingALL,IsNew=',nEvolvingALL,IsNew

    if(.not. IsNew) RETURN

    ! Store number of evolving blocks
    nEvolvingLastALL = nEvolvingALL

    ! Update the global information about evolving blocks
    call MPI_allgather(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
         iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)


    if(IsFirstSelect)then
       ! The original setting is that all blocks are part of the 
       ! steady boundary, so that they are all checked for change.
       ! Here the block type needs to be changed to steady.
       where(iTypeAdvance_B == SteadyBoundBlock_) iTypeAdvance_B = SteadyBlock_
       IsFirstSelect = .false.
    endif

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
         write(*,*)'part_steady finished with nSkipped,Steady,Bound,ExplALL=',&
         count(iTypeAdvance_BP == SkippedBlock_), & 
         count(iTypeAdvance_BP == SteadyBlock_), & 
         count(iTypeAdvance_BP == SteadyBoundBlock_), &
         count(iTypeAdvance_BP == ExplBlock_)

  end subroutine part_steady_select

end module ModPartSteady
