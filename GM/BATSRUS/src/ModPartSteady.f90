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

  implicit none

  save

  private ! except

  logical, public :: UsePartSteady ! True if the part steady scheme used

  public IsEvolving_B       ! Logical array for evolving blocks
  public IsSteadyBoundary_B ! Logical array for blocks around evolving blocks
  public part_steady_init   ! Allocate and initialize variables
  public part_steady_clean  ! Deallocate variables
  public part_steady_switch ! Switch to and from the evolving/bounday blocks
  public part_steady_select ! Select the evolving/bounday blocks


  logical, parameter :: DoDebug = .true.

  ! Logicals for the local PE (*_B) and for all PE-s (*_BP) for
  ! 1. Blocks which were unused originally.
  ! 2. Blocks which have changed already are 'evolving'.
  ! 3. Blocks surrounding the evolving blocks from the 'steady boundary'.
  logical, allocatable :: &
       IsUnusedOrig_B(:),     IsUnusedOrig_BP(:,:), &
       IsEvolving_B(:),       IsEvolving_BP(:,:), &
       IsSteadyBoundary_B(:), IsSteadyBoundary_BP(:,:)


  ! The variables to be checked are indexed from MinCheckVar to MaxCheckVar
  integer :: MinCheckVar = 1, MaxCheckVar = nVar

  ! The significant relative and absolute changes for all variables
  real :: RelativeEps_V(nVar) = 0.001, AbsoluteEps_V(nVar) = 1.0e-4

contains
  !===========================================================================
  subroutine part_steady_init

    if(allocated(IsEvolving_B)) RETURN

    allocate( &
         IsUnusedOrig_B( MaxBlock),              &
         IsUnusedOrig_BP(MaxBlock,0:nProc-1),    &
         IsEvolving_B( MaxBlock),                &
         IsEvolving_BP(MaxBlock,0:nProc-1),      &
         IsSteadyBoundary_B( MaxBlock),          &
         IsSteadyBoundary_BP(MaxBlock,0:nProc-1) &
         )

    if(DoDebug)write(*,*)'part_steady_init: allocated variables'

  end subroutine part_steady_init

  !===========================================================================
  subroutine part_steady_clean

    if(.not.allocated(IsEvolving_B)) RETURN

    deallocate(IsEvolving_B, IsEvolving_BP, IsSteadyBoundary_B)

    if(DoDebug)write(*,*)'part_steady_clean: deallocated variables'

  end subroutine part_steady_clean
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
       IsUnusedOrig_B( 1:nBlockMax)   = UnusedBLK(1:nBlockMax)
       IsUnusedOrig_BP(1:nBlockMax,:) = UnusedBlock_BP(1:nBlockMax,:)

       ! Select UnusedBLK = IsUnusedOrig or not(evolving or steady_boundary)
       UnusedBLK(1:nBlockMax) = IsUnusedOrig_B(1:nBlockMax) .or. &
            .not. (IsEvolving_B(1:nBlockMax) &
            .or.   IsSteadyBoundary_B(1:nBlockMax) )

       UnusedBlock_BP(1:nBlockMax,:) = IsUnusedOrig_BP(1:nBlockMax,:) .or. &
            .not. (IsEvolving_BP(1:nBlockMax,:) &
            .or.   IsSteadyBoundary_BP(1:nBlockMax,:) )

       ! Change the decomposition index
       iNewDecomposition=mod(iNewDecomposition+1, 10000)

    else
       ! Restore the original unused blocks
       UnusedBLK(1:nBlockMax)        = IsUnusedOrig_B(1:nBlockMax)
       UnusedBlock_BP(1:nBlockMax,:) = IsUnusedOrig_BP(1:nBlockMax,:)

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

    use ModMain,     ONLY: East_, Top_, Time_Accurate, Time_Simulation
    use ModProcMH,   ONLY: iComm
    use ModAdvance,  ONLY: State_VGB, StateOld_VCB, nI, nJ, nK, nIJK
    use ModParallel, ONLY: NOBLK, NeiLev, NeiPe, NeiBlk
    use ModMpi

    logical, intent(out) :: IsNew

    logical :: IsFirstCall = .true.
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
       IsSteadyBoundary_B  = .true.
       IsSteadyBoundary_BP = .true.

       ! Assign all blocks to be non-evolving for now
       IsEvolving_B  = .false.
       IsEvolving_BP = .false.

       ! There is nothing else to do
       RETURN
    end if

    ! The first time accurate time step with zero progress cannot be used
    if(Time_Accurate .and. Time_Simulation == 0) RETURN

    ! Check the steady boundary blocks for change
    do iBlock = 1, nBlock

       if(.not.IsSteadyBoundary_B(iBlock)) CYCLE

       ! Calculate the change the state
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

       ! Check if the change is significant
       if(any(dState_V > 1.0)) then
          IsEvolving_B(iBlock)       = .true.
          IsSteadyBoundary_B(iBlock) = .false.
       end if

    end do

    ! Count the number of evolving blocks
    nEvolving = count(IsEvolving_B(1:nBlock))
    call MPI_allreduce(nEvolving, nEvolvingALL, 1, MPI_INTEGER, MPI_SUM, &
          iComm, iError)

    ! If no new evolving blocks were found, simply return
    IsNew = nEvolvingALL /= nEvolvingLastALL

    if(DoDebug)write(*,*)'nEvolvingALL,IsNew=',nEvolvingALL,IsNew

    if(.not. IsNew) RETURN

    ! Store number of evolving blocks
    nEvolvingLastALL = nEvolvingALL

    ! Update the global information about evolving blocks
    call MPI_allgather(IsEvolving_B, MaxBlock, MPI_LOGICAL, &
         IsEvolving_BP, MaxBlock, MPI_LOGICAL, iComm, iError)

    ! Check the unevolving blocks if they became part of the steady boundary
    IsSteadyBoundary_B = .false.
    BLOCKS: do iBlock = 1, nBlock

       ! Skip all other blocks
       if(  IsUnusedOrig_B(iBlock) .or. &
            IsEvolving_B(iBlock)  .or. &
            IsSteadyBoundary_B(iBlock)) CYCLE

       write(*,*)'checking iBlock=',iBlock

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
             if(IsEvolving_BP(jBlock, jProc)) then
                IsSteadyBoundary_B(iBlock) = .true.
                CYCLE BLOCKS
             end if
          end do SUBFACES
       end do FACES

    end do BLOCKS

    ! Update the global information about steady boundary blocks
    call MPI_allgather(IsSteadyBoundary_B, MaxBlock, MPI_LOGICAL, &
         IsSteadyBoundary_BP, MaxBlock, MPI_LOGICAL, iComm, iError)

    write(*,*)'nSteadyBoundary=',count(IsSteadyBoundary_BP)

  end subroutine part_steady_select

end module ModPartSteady
