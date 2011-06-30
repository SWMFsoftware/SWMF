!^CFG COPYRIGHT UM

module BATL_amr_criteria


  use BATL_tree, ONLY: nDesiredRefine,nNodeRefine,nDesiredCoarsen, &
       nNodeCoarsen, DoStrictAmr, iRank_A, Rank_A, nNodeSort


  implicit none

  SAVE

  private ! except

  public set_amr_criteria
  public clean_amr_criteria
  public read_amr_criteria
  public test_amr_criteria


  ! Choosing with blocks we want to refine is based a list of criteria 
  ! and a set of upper (refine)  and lower (coarsen) limits. The criteria 
  ! can be can be external or be calculated internally by estimating the 
  ! numerical errors (calc_error_criteria) based on the state variables.
  ! 
  integer, public            :: nAmrCrit = 0
  integer                    :: nAmrCritUsed = 0
  real, public, allocatable  :: AmrCrit_IB(:,:)

  !--- Moved to BATL_tree
  ! We can also specify the percentage of blocks we want to refine. For doing
  ! this we need to sort them into a priority list iRank_A. This priority list
  ! can also be used for refining/coarsening blocks to the point where we 
  ! to the point where we have no more blocks available, and not stopping 
  ! the program ( used with BATL_tree )
  !integer,public, allocatable :: iRank_A(:) !--- Moved to BATL_tree

  ! DoSortAmrCrit : true/false for generating the priority list iRank_A
  ! DoSoftAmrCrit : Standard behavior (false) is that criteria has to be 
  ! refined/coarsened or we will stop the program. If true we will 
  ! refine/coarsen them as much as we can, this is the default behavior for
  ! amr by percentage.
  logical,public :: &
       DoSortAmrCrit = .true.,&
       DoSoftAmrCrit = .false. ,&
       DoAutoAmr     = .false.
  
       !DoStrictAmr = .true. ,& !--- Moved to BATL_tree

  ! Try to make geometric dependence for refinement/coarsening
  logical,public :: &
       DoGeometryAmr = .false. ,&
       DoCritAmr     = .false.

  !--- Moved to BATL_tree
  ! nDesiredRefine and nDesiredCoarsen set the number of blocks that we
  ! want refined/coarsen by percentage of with DoSoftAmrCrit = .true.
  ! nNodeRefine and nNodeCoarsen set the number of blocks that will be
  ! refined/coarsen or the program with stop. Usually associated with 
  ! criteria
  !integer,public :: nDesiredRefine,nNodeRefine, &
  !     nDesiredCoarsen, nNodeCoarsen

  
  ! Local variables

  ! Percentage with want refined and coarsen
  real ::  PercentRefine=0.0, PercentCoarsen=0.0

  ! Threshold limits for refine or unrefined the grid (length nCrit) 
  real, allocatable, dimension(:)    :: CoarsenCrit_I, RefineCrit_I,&
       CoarsenCritAll_I, RefineCritAll_I

  integer, allocatable:: iVarCrit_I(:) ! Index to variables

  ! Parameters used by calc_error_criteria to estimate the errors
  real :: cAmrWavefilter = 1.0e-2
  real, parameter :: cEpsilon = 1.0d-16 ! avoid zero in denominator

  integer :: nExtCrit = 4 ! Number of External criteria, from BATSRUS
  integer :: nExtCritUsed = 0 ! Number of External criteria actually used, 
                              !4 :symmetry (BATSRUS)
  integer :: nIntCrit = 0 ! Number of internal criteria, 2nd order err estimate
  integer :: nAmrCritOld = 0, nBlockOld = 0

  ! How large the relative discrepancy in the Criteria for 
  ! symmetry points can be and still be refined in the same 
  ! way and the same for geometric variation
  real, parameter :: DeltaSymCrit = 1.0e-2
  real, parameter :: DeltaSymGeo  = 1.0e-6

  ! The max differens in the criteia to say thay have the same value,
  ! taking care of numerical fluctations
  real :: DeltaCritera = 0.0

  ! Converting index form running block and proc number to node numbers
  integer, allocatable :: iNode_I(:)

contains
  !============================================================================

  subroutine set_amr_criteria(nVar, State_VGB, nInCritExtUsed, CritExt_IB, &
       CoarsenCritExt_I, RefineCritExt_I)

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
         MaxBlock,nBlock
    use BATL_tree, ONLY: Unused_B
    use BATL_mpi, ONLY: iProc

    integer,  intent(in)  :: nVar
    real,    intent(in):: &                 ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    integer, intent(in), optional :: nInCritExtUsed
    real, intent(in), optional :: CritExt_IB(nExtCrit, nBlock)
    real, intent(in), optional :: CoarsenCritExt_I(nExtCrit),&
         RefineCritExt_I(nExtCrit)

    integer :: iCrit, iBlock
    !-----------------------------------------------------------------------
    if(DoGeometryAmr) &
         call CON_stop("set_amr_criteria :: DoGeometryAmr not implemented")

    !-------------- Setting number of criteria we are working with ----------
    nExtCritUsed = 0
    if(present(nInCritExtUsed)) &
         nExtCritUsed = nInCritExtUsed

    nAmrCritUsed = nIntCrit + nExtCritUsed
    nAmrCrit = nIntCrit + nExtCrit

!!$    if(iProc == 0) &
!!$         write(*,'(a50,4(i5))') &
!!$         " nAmrCrit, nIntCrit, nExtCrit, nAmrCritUsed =", &
!!$         nAmrCrit, nIntCrit, nExtCrit, nAmrCritUsed

    !------------ Collect all criteria external and internals -------------


    ! rescale arrays if necessary based on criteria and blocks
    if(nAmrCrit /= nAmrCritOld) then
       if(allocated(CoarsenCritAll_I)) deallocate( &
            CoarsenCritAll_I, RefineCritAll_I)
       allocate(CoarsenCritAll_I(nAmrCrit), RefineCritAll_I(nAmrCrit))
    end if

    if(nAmrCrit /= nAmrCritOld .or. nBlock /= nBlockOld) then
       if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
       allocate(AmrCrit_IB(nAmrCrit,nBlock))
       nAmrCritOld=nAmrCrit
       nBlockOld  = nBlock
       AmrCrit_IB = 0.0
    end if

    AmrCrit_IB = 0.0

    ! add external criteria into the list of all criteria
    if(present(CritExt_IB)) then
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          do iCrit = 1, nExtCrit
             AmrCrit_IB(nIntCrit+iCrit,iBlock) =  CritExt_IB(iCrit,iBlock)
          end do
       end do

       ! add external refinement and coarsening thresholds to internals
       do iCrit=1,nExtCrit
          RefineCritAll_I(nIntCrit+iCrit)  = RefineCritExt_I(iCrit)
          CoarsenCritAll_I(nIntCrit+iCrit) = CoarsenCritExt_I(iCrit)
          !print *,"iCrit ext :", iCrit, size(RefineCritExt_I)
          !print *," RefineCritExt_I(iCrit) : ",RefineCritExt_I(iCrit)
       end do
    end if

    ! Estimation of the numerical error
    if(nIntCrit > 0) &
         call calc_error_criteria(nVar, State_VGB )

    if(DoSortAmrCrit .or. .not.DoStrictAmr) then
       ! we make a amr priority list
       call sort_amr_criteria
    else
       ! refine only based on criteria
       call apply_unsorted_criteria
    end if


  end subroutine set_amr_criteria
  !===========================================================================
  subroutine sort_amr_criteria

    ! the routine will find the candidate for refinement or coursing based
    ! upon the percentage of refinement and the percentage wanted refined and
    ! and coarsened. It will ONLY modify the criteria the values in the 
    ! criteria list itself. The refinement criteria will overrule the 
    ! percentage refinement.

    use BATL_mpi, ONLY: iComm, iProc, nProc
    use ModMpi
    use BATL_size, ONLY: nBlock
    use ModSort, ONLY: sort_quick
    use BATL_tree, ONLY: Unused_BP,iStatusNew_A, Refine_, Coarsen_, &
         iNode_B, nNode, nChild, nNodeUsed, Unset_, diffRange

    ! Array containing all criteria for all blocks gathered
    real, allocatable :: AllCrit_II(:,:)

    ! CritSort_II contains all criteria for all blocks
    ! iIdxSort_II gives the new indexes after sorting
    ! and reused for other functions.
    integer, allocatable :: iIdxSort_II(:,:)
    real, allocatable :: CritSort_II(:,:)
    integer, allocatable :: iRank_I(:)
    real, allocatable :: Rank_I(:)
    integer :: iCrit, iBlock, iProces, k
    integer :: iCritSort, iSort
    integer :: iError
    real    :: Diff, Coeff, Crit

    integer :: nTotBlocks, nBlockMax, iTotalCrit, iHelpCrit
    integer, allocatable :: nBlock_P(:), nReciveCont_P(:), nRecivDisp_P(:)
    !-----------------------------------------------------------------------


    !if(iProc == 0) write(*,*) "sort_amr_criteria"


    iTotalCrit = 2
    iHelpCrit  = 3

    ! COMMENTS FOR FUTURE : we can use the Unused_B to only pick and send
    ! the data that is actively used, reducing communication


    !------------ Get information on all the blocks on all process -----------


    ! collect nBlock from each process 
    allocate(nBlock_P(nProc),nReciveCont_P(nProc),nRecivDisp_P(nProc))
    call MPI_Allgather(nBlock, 1,MPI_INTEGER, nBlock_P, 1, &
         MPI_INTEGER, iComm, iError)
    nBlockMax = maxval(nBlock_P)
    nTotBlocks = sum(nBlock_P)

    !------------ Collect all criteria for sorting and ranking ------------

    ! collect all criteria for the total grid
    allocate(AllCrit_II(nAmrCrit,nTotBlocks))
    AllCrit_II = -77

    ! store the Node indexing iNode_B for all processes
    if(allocated(iNode_I)) deallocate(iNode_I)
    allocate(iNode_I(nBlockMax*nProc))
    iNode_I = 0

    ! Set up the displacement and size for collecting the data for 
    ! AllCrit_II
    nReciveCont_P = nBlock_P*nAmrCrit
    nRecivDisp_P(1) = 0
    do iProces=2,nProc 
       nRecivDisp_P(iProces) = &
            nRecivDisp_P(iProces-1)+nBlock_P(iProces-1)*nAmrCrit
    end do

    ! gathering the criteria 
    call MPI_allgatherv(AmrCrit_IB, nAmrCrit*nBlock, MPI_REAL, &
         AllCrit_II, nReciveCont_P, nRecivDisp_P, &
         MPI_REAL, iComm, iError)


    ! Set up the displacement and size for collecting the data for 
    ! iNode_I
    nReciveCont_P = nBlock_P
    nRecivDisp_P(1) = 0
    do iProces=2,nProc 
       nRecivDisp_P(iProces) = nRecivDisp_P(iProces-1)+nBlock_P(iProces-1)
    end do

    ! gathering the node indexes
    call MPI_allgatherv(iNode_B, nBlock,MPI_INTEGER, &
         iNode_I, nReciveCont_P, nRecivDisp_P, &
         MPI_INTEGER, iComm, iError)

    ! Setting up arrays for sorting and index handling
    allocate( &
         CritSort_II(nTotBlocks,1), &
         iIdxSort_II(nTotBlocks,iHelpCrit)  )

    ! iRank_A is public, so deallocate it if necessary
    if(allocated(iRank_A)) deallocate(iRank_A,Rank_A)
    allocate(iRank_I(nTotBlocks), Rank_I(nTotBlocks),&
         iRank_A(nTotBlocks), Rank_A(nTotBlocks))

    iRank_A    = 0
    Rank_A     = 0.0

    iRank_I = 0
    Rank_I  = 0.0
    iIdxSort_II  = 0
    CritSort_II  = 0.0

    Coeff = 1./real(nAmrCritUsed*nTotBlocks) ! normalization for secondary rank
    iRank_I = nNode+1 ! larger than any possible rank
    iIdxSort_II(:,iHelpCrit) = 0 ! store info about if a block is already 
    !                               marked for refinement(+1) or coarsening(-1)
    nNodeRefine  = 0
    nNodeCoarsen = 0

    do iCrit=1,nAmrCritUsed


       ! copy criteria data into the sorting arrays
       iSort = 0
       do iProces = 1, nProc
          do iBlock = 1, nBlock_P(iProces)
             if(Unused_BP(iBlock,iProces-1)) CYCLE
             iSort = iSort +1
             CritSort_II(iSort,1) = &
                  AllCrit_II(iCrit,iBlock+nRecivDisp_P(iProces))
             iIdxSort_II(iSort, 1) = iSort
          end do
       end do
       nNodeSort = iSort

       ! Sort each criteria. Pass in first index of array arguments
       call sort_quick(nNodeSort, &
            CritSort_II(1:nNodeSort,1), iIdxSort_II(1:nNodeSort,1))


       iIdxSort_II(:,iTotalCrit) = iIdxSort_II(:,1)
       k = 1
       iSort = iIdxSort_II(1,iTotalCrit)
       !Value of first element (min)
       iIdxSort_II(iSort,1) = k
       Crit = CritSort_II(iSort,1)


       ! group together criteass which has a diffrence of 
       ! less then DeltaCritera by giving them the same
       ! sort ranking
       do iCritSort = 2,nNodeSort
          iSort = iIdxSort_II(iCritSort,iTotalCrit)
          if((CritSort_II(iSort,1)- Crit) < DeltaCritera) then
             iIdxSort_II(iSort,1) = k
          else
             k = k + 1
             iIdxSort_II(iSort,1) = k
             Crit = CritSort_II(iSort,1)
          end if
       end do

       do iSort = 1,nNodeSort
          iCritSort = iIdxSort_II(iSort,iTotalCrit)

          ! Use the minimum position (=rank) of the node in sorted criteria 
          ! array to set up the new ranking order for each node when doing amr
          iRank_I(iCritSort) = min(iRank_I(iCritSort),iIdxSort_II(iCritSort,1))

          ! if some of the nodes have the same minimal rank we will use the 
          ! average rank of all criteria as a secondary priority.
          Rank_I(iCritSort) = Rank_I(iCritSort)+Coeff*iIdxSort_II(iCritSort,1)
       end do



       ! To make sure that blocks that are marked for refinement/coarsening 
       ! by thresholds will come at the top/bottom of the ranking list 
       ! we will shift the values with the total number of blocks used.
       ! COMMENTS FOR FUTURE : the loop may be split for one that goes 
       ! form the top and one that goes form bottom so it do not need to go 
       ! though the whole list of elements

       !print *," iCrit =", iCrit
       !print *," RefineCritAll_I(iCrit) : ", RefineCritAll_I(iCrit),size(RefineCritAll_I)
       !print *," CoarsenCritAll_I(iCrit): ",CoarsenCritAll_I(iCrit),size(CoarsenCritAll_I) 
       ! Only if Criterias is in use
       if(RefineCritAll_I(iCrit) > -0.5 .and. CoarsenCritAll_I(iCrit) > -0.5) then
          do iSort = nNodeSort,1,-1
             iCritSort = iIdxSort_II(iSort,iTotalCrit)
             ! Block has already ben marked
             !if(iIdxSort_II(iCritSort,iHelpCrit) /= 0 ) CYCLE

             if(CritSort_II(iCritSort,1) > RefineCritAll_I(iCrit) &
                  .and. iIdxSort_II(iCritSort,iHelpCrit) /= 1) then
                ! Satisfies refinement threshold and not yet marked for refinement

                ! Shift up the rank above sorted values                
                Rank_I(iCritSort) = Rank_I(iCritSort) + nTotBlocks+1

                ! If node was originally marked for coarsening, 
                ! now there is one fewer node to be coarsened
                if(iIdxSort_II(iCritSort,iHelpCrit) == -1) &
                     nNodeCoarsen = nNodeCoarsen - 1

                ! Noe node is marked fore refinement
                nNodeRefine = nNodeRefine +1
                iIdxSort_II(iCritSort,iHelpCrit) = 1

             else if(CritSort_II(iCritSort,1) < CoarsenCritAll_I(iCrit) &
                  .and. iIdxSort_II(iCritSort,iHelpCrit) == 0) then
                ! Satisfies coarsening threshold and not yet marked at all !

                ! Shift down the rank below sorted values
                Rank_I(iCritSort) = Rank_I(iCritSort)-nTotBlocks-1

                ! Now node is marked for coarsening
                nNodeCoarsen = nNodeCoarsen +1
                iIdxSort_II(iCritSort,iHelpCrit) = -1 
             end if
          end do
       end if

    end do


    ! the final rank is a sum of the primary minimum rank
    ! and secondary average rank:
    Rank_I(1:nNodeSort) = iRank_I(1:nNodeSort) + Rank_I(1:nNodeSort)

    ! Get the number of elements to refine by percentage
    ! taking into account that each refinements generate nChild-1 
    ! new blocks
    nDesiredRefine  = floor(PercentRefine*nNodeUsed/(100.0*(nChild-1)))
    nDesiredCoarsen = floor(PercentCoarsen*nNodeUsed/100.0)

    ! make sure thay do not take more then thay can
    nDesiredRefine  = min(nNodeSort,nDesiredRefine)
    nDesiredCoarsen = min(nNodeSort,nDesiredCoarsen)

    ! nNodeRefine and nNodeCoarsen are based on thresholds.
    ! If DoSoftAmrCrit is false, these have to be refined/coarsened,
    ! or else the code stops.
    ! If DoSoftAmrCrit is true, we change the threshold based 
    ! values into "desired" values, so the code behaves as if 
    ! the criteria were based on percentages.
    if(DoSoftAmrCrit .or. .not.DoStrictAmr) then
       nDesiredRefine  = max(nDesiredRefine,  nNodeRefine)
       nDesiredCoarsen = max(nDesiredCoarsen, nNodeCoarsen)
       nNodeRefine  = 0
       nNodeCoarsen = 0
    end if

    !if(iProc == 0) then
    !   write(*,'(2(a30, i10))') " BATL nDesiredRefine  = ", &
    !        nDesiredRefine,",    nNodeRefine  = ",  nNodeRefine
    !   write(*,'(2(a30, i10))') " BATL nDesiredCoarsen = ", &
    !        nDesiredCoarsen,",    nNodeCoarsen = ", nNodeCoarsen
    !end if


    ! store the sorted indexing list for Rank_A
    iIdxSort_II(:,iTotalCrit) = 0 

    ! The quick sort algorithm assume a continuous indexing order,
    ! this will not necessary be true so we need to make a 
    ! mapping between the block numbers and the node numbers when
    ! there are unused blocks. This is stored in iIdxSort_II(iSort,1)
    k = 0
    iSort = 0
    do iProces = 1,nProc
       do iBlock = 1, nBlock_P(iProces)
          k = k +1
          if(Unused_BP(iBlock,iProces-1)) CYCLE
          iSort = iSort+1
          iIdxSort_II(iSort,1) = iNode_I(k)
       end do
    end do

    ! Finding final rank ordering
    call sort_quick(nNodeSort, Rank_I, iIdxSort_II(1,iTotalCrit))

    ! map local indexing in this subroutine to the node index
    do iSort = 1, nNodeSort
       iCritSort = iIdxSort_II(iSort,iTotalCrit)
       iRank_A(iSort) =iIdxSort_II(iCritSort,1)
       Rank_A(iSort) = Rank_I(iCritSort)
    end do



    ! making sure that the refinment and coarsning are starting from
    ! a position with have a transition in its criteria larger then
    ! diffRange
    if(nDesiredRefine > nNodeRefine .and.  nDesiredRefine < nNodeSort )then
       k=-1
       do iSort = nNodeSort-nDesiredRefine, nNodeSort-nNodeRefine
          diff =  abs(Rank_A(iSort) - Rank_A(nNodeSort- nDesiredRefine+1)) 
          if( diff > diffRange) EXIT
          k=k+1
       end do
       nDesiredRefine = nDesiredRefine - max(0,k)
    end if

    if(nDesiredCoarsen > nNodeCoarsen .and. nDesiredCoarsen < nNodeSort )then
       k=-1
       do iSort= nDesiredCoarsen+1, max(nNodeCoarsen,1), -1
          diff = abs(Rank_A(iSort) - Rank_A(nDesiredCoarsen))
          if(diff >diffRange) EXIT
          k=k+1
       end do
       nDesiredCoarsen = nDesiredCoarsen - max(0,k)
    end if

    ! Give iStatusNew_A based on the iRank_A and number of blocks we want for
    ! refinment and coarsning
    do iSort = nNodeSort, nNodeSort-max(nDesiredRefine,nNodeRefine)+1, -1
       iStatusNew_A(iRank_A(iSort)) = Refine_
    end do

    do iSort = 1, max(nDesiredCoarsen, nNodeCoarsen)
       if(iStatusNew_A(iRank_A(iSort)) /= Refine_) &
            iStatusNew_A(iRank_A(iSort)) =  Coarsen_
    end do

    deallocate(AllCrit_II)
    deallocate(CritSort_II, iIdxSort_II)
    deallocate(nBlock_P, nReciveCont_P, nRecivDisp_P)
    deallocate(iRank_I, Rank_I)

  end subroutine sort_amr_criteria


  !============================================================================

  ! The AMR criteria method used here are adapted from L{\"o}hner (1987).
  ! The flash3 user guide have also a description used.
  ! The method use the relative error ,E (Factor), of the approximation to 
  ! the differential equation:
  !
  !  dU      DU         d^2 U
  ! ----  = ----  + DX  ------
  !  dx      DX          dx^2
  ! 
  ! Giving a relative error
  !
  !             d^2 U
  !            -------
  !              dx^2
  !    E  =  -------------
  !           1     dU
  !          ---  ------
  !          DX     dx
  !
  ! To avoid problems with small linear gradient fluctuations in the 
  ! denominator we add it with a fraction of its average value.
  ! And to avoid zero in the denominator we add it with cEpsilon.
  ! The final expression is:
  !
  !                               d^2 U
  !                              -------
  !                               dx^2 
  !    E  =  ----------------------------------------------
  !           1     dU        1
  !          ---  ------  + ------ cAmrWavefilter*<U> + cEpsilon 
  !          DX     dx       DX^2
  !
  ! We use a wide stencil to do the refinement in time and 
  ! neglected the diagonal terms. 
  ! 
  ! The error factor for each variable will be normalized before
  ! threshold evaluation and also the sum of all the factors over
  ! all variables will be evaluated for deciding the refinement.
  !
  ! The subroutine only set iStatusNew_A array.
  !
  ! OPTIONAL:
  ! The user can sent in additional factors  used for deciding refinement
  ! together with there refinements criteria. This will be then used
  ! in addition to the error estimate described earlier.

  subroutine calc_error_criteria(nVar, State_VGB )!, & 
      ! nCritExt, CritExt_IB, CoarsenCritExt_I, RefineCritExt_I)

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
	 nI, nJ, nK, MaxDim, nDim, MaxBlock, nBlock
    use BATL_tree, ONLY: Unused_B
    use BATL_mpi, ONLY: iProc

    integer,  intent(in)   :: nVar
    real,    intent(in) :: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    !integer, intent(in), optional :: nCritExt ! num of external criteria
    !real, intent(inout), optional :: CritExt_IB(:,:)
    !real, intent(in),    optional :: CoarsenCritExt_I(:), RefineCritExt_I(:)

    real :: Crit, Crit_D(MaxDim)
    real :: Numerator, Denominator
    integer:: iBlock, iCrit, i, j, k, iVar
    ! number of blocks set out for refining and coarsning
    !------------------------------------------------------------------------
    Crit      = 1.0
    Crit_D    = 0.0
    Numerator   = 0.0
    Denominator = 0.0

    !if(iproc == 0) write(*,*) "calc_error_criteria"

    if(nIntCrit > nVar) &
         call CON_stop("calc_error_criteria :: More criteria then variables")
    

    nAmrCrit = nIntCrit + nExtCrit
    if(nAmrCrit /= nAmrCritOld) then
       if(allocated(CoarsenCritAll_I)) deallocate( &
            CoarsenCritAll_I, RefineCritAll_I)
       allocate(CoarsenCritAll_I(nAmrCrit), RefineCritAll_I(nAmrCrit))
    end if

    if(nAmrCrit /= nAmrCritOld .or. nBlock /= nBlockOld) then
       if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
       allocate(AmrCrit_IB(nAmrCrit,nBlock))
       AmrCrit_IB = 0.0
       nBlockOld=nBlock
    end if
    
    nAmrCritOld=nAmrCrit

    ! Calculate error estimates (1..nIntCrit)
    BLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Check refinement first
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          ! Only variables indexed in iVarCrit_I will decide the refinement
          do iCrit = 1, nIntCrit 
             iVar = iVarCrit_I(iCrit)

             Numerator = abs( &
                  State_VGB(iVar,i-2,j,k,iBlock)   - &
                  2.0*State_VGB(iVar,i,j,k,iBlock) + &
                  State_VGB(iVar,i+2,j,k,iBlock) )
             Denominator = (&
                  abs(State_VGB(iVar,i+2,j,k,iBlock) - &
                  State_VGB(iVar,i,j,k,iBlock)) + &
                  abs(State_VGB(iVar,i,j,k,iBlock) - &
                  State_VGB(iVar,i-2,j,k,iBlock))) + &
                  cAmrWavefilter * (&
                  abs(State_VGB(iVar,i+2,j,k,iBlock)) + &
                  abs(2.0*State_VGB(iVar,i,j,k,iBlock)) + &
                  abs(State_VGB(iVar,i-2,j,k,iBlock))) 
             Crit_D(1) = (Numerator/max(Denominator,cEpsilon)) 

             if(nDim >= 2) then
                Numerator =  abs( &
                     State_VGB(iVar,i,j-2,k,iBlock)   - &
                     2.0*State_VGB(iVar,i,j,k,iBlock) + &
                     State_VGB(iVar,i,j+2,k,iBlock) )
                Denominator = (&
                     abs(State_VGB(iVar,i,j+2,k,iBlock) - &
                     State_VGB(iVar,i,j,k,iBlock)) + &
                     abs(State_VGB(iVar,i,j,k,iBlock) - &
                     State_VGB(iVar,i,j-2,k,iBlock))) + &
                     cAmrWavefilter * (&
                     abs(State_VGB(iVar,i,j+2,k,iBlock)) + &
                     abs(2.0*State_VGB(iVar,i,j,k,iBlock)) + &
                     abs(State_VGB(iVar,i,j-2,k,iBlock)))
                Crit_D(2) = (Numerator/max(Denominator,cEpsilon))
             end if

             if(nDim >= 3) then
                Numerator = abs( &
                     State_VGB(iVar,i,j,k-2,iBlock)   - &
                     2.0*State_VGB(iVar,i,j,k,iBlock) + &
                     State_VGB(iVar,i,j,k+2,iBlock))
                Denominator = (&
                     abs(State_VGB(iVar,i,j,k+2,iBlock) - &
                     State_VGB(iVar,i,j,k,iBlock)) + &
                     abs(State_VGB(iVar,i,j,k,iBlock) - &
                     State_VGB(iVar,i,j,k-2,iBlock)) + &
                     cAmrWavefilter * (&
                     abs(State_VGB(iVar,i,j,k+2,iBlock)) + &
                     abs(2.0*State_VGB(iVar,i,j,k,iBlock)) + &
                     abs(State_VGB(iVar,i,j,k-2,iBlock))))
                Crit_D(3) = (Numerator/max(Denominator,cEpsilon))
             end if

             Crit = sqrt(sum(Crit_D**2))/nDim
             if( Crit >  AmrCrit_IB(iCrit,iBlock)) &
                  AmrCrit_IB(iCrit,iBlock) = Crit

          end do !end nVar
       end do; end do; end do

    end do BLOCK

    ! Fill in thresholds
    do iCrit=1,nIntCrit
       RefineCritAll_I(iCrit)  = RefineCrit_I(iCrit)
       CoarsenCritAll_I(iCrit) = CoarsenCrit_I(iCrit)
    end do

  end subroutine calc_error_criteria

  !============================================================================

  subroutine apply_unsorted_criteria

    use BATL_mpi,  ONLY: iProc 
    use BATL_size, ONLY: nBlock
    use BATL_tree, ONLY: iStatusNew_A, Refine_, Coarsen_, &
         Unused_B, iNode_B

    integer:: iBlock, iCrit
    ! number of blocks set out for refining and coarsning
    logical :: DoCoarsen
    !----------------------------------------------------------------------

    !if(iProc == 0) write(*,*) "apply_unsorted_criteria"

    BLOCK2:do iBlock = 1, nBlock

       if(Unused_B(iBlock)) CYCLE
       DoCoarsen = .true.

       !if(iproc==0)&
       !     write(*,'(2(f16.6),a8,4(f16.6))') &
       !     AmrCrit_IB(1:2,iBlock), " :: ", &
       !     RefineCritAll_I(1:2), CoarsenCritAll_I(1:2)
   
       do iCrit = 1, nAmrCritUsed

          ! Decide refinement based on normalized error factor for each
          ! chosen variable
          ! If one cell in the block needs refinement the block will
          ! be refined, But only if all cells in the block want to be
          ! coarsen the block will be flagged for coarsening

          if( AmrCrit_IB(iCrit,iBlock) > RefineCritAll_I(iCrit)) then
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             DoCoarsen = .false.
             CYCLE BLOCK2
          else if(AmrCrit_IB(iCrit,iBlock)  > CoarsenCritAll_I(iCrit)) then
             DoCoarsen = .false.
          end if
       end do

       if(DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_ 

    end do BLOCK2

  end subroutine apply_unsorted_criteria

  !============================================================================
  subroutine  read_amr_criteria(NameCommand)

    use ModReadParam, ONLY: read_var
    use BATL_tree, ONLY: MaxTotalBlock

    character(len=*), intent(in) :: NameCommand
    character (len=20) :: TypeAmr
    integer :: iCrit
    character(len=*), parameter:: NameSub='BATL_tree::read_amr_criteria'
    !-------------------------------------------------------------------------
    !nCrit          = 0
    !nAmrCrit       = 0
    !nAmrCritOld    = 0
    !cAmrWavefilter = 1.0e-2
    !nBlockOld      = 0

    select case(NameCommand)
    case("#AMRERRORCRIT") 
       call read_var('AmrWavefilter',cAmrWavefilter)  
       call read_var('nCrit', nIntCrit)
       nAmrCrit = nIntCrit
       nAmrCritOld = nIntCrit
       if(allocated(CoarsenCrit_I)) then
          deallocate(CoarsenCritAll_I, RefineCritAll_I)
          deallocate(CoarsenCrit_I, &
               RefineCrit_I, iVarCrit_I)
       end if
       allocate(CoarsenCrit_I(nIntCrit), &
            RefineCrit_I(nIntCrit),iVarCrit_I(nIntCrit))
       allocate(CoarsenCritAll_I(nIntCrit), RefineCritAll_I(nIntCrit))
       do iCrit = 1, nIntCrit
          call read_var('iVar',iVarCrit_I(iCrit))
          call read_var('CoarsenCrit',CoarsenCrit_I(iCrit))
          call read_var('RefineCrit',RefineCrit_I(iCrit))
       end do
       DoCritAmr = .true.
       DoAutoAmr = .true.
    case("#AMR") ! compatibilety with old system with BATSRUS amr options
       call read_var('PercentCoarsen', PercentCoarsen)
       call read_var('PercentRefine' , PercentRefine)
       call read_var('MaxTotalBlock',  MaxTotalBlock) 
       DoSortAmrCrit = PercentCoarsen > 0.0 .or. PercentRefine > 0.0       
       !case("#AMRTYPE")
       !   call read_var('DoAutoAmr', DoAutoAmr)
       !   call read_var('TypeAmr',TypeAmr)
       !   select case(TypeAmr)
       !   case("geometry")
       !      DoGeometryAmr = .true.
       !   case("criteria")
       !      DoGeometryAmr = .false.
       !   case("combined")
       !      DoGeometryAmr = .true.
       !   case default
       !      call CON_stop(NameSub//': unknown TypeAmr='//TypeAmr)
       !!   end select
       !   call read_var('IsStrictAmr'  ,DoStrictAmr)
       !call read_var('IsSoftAmrCrit',DoSoftAmrCrit)
    case("#AMRLIMIT")
       call read_var('PercentCoarsen', PercentCoarsen)
       call read_var('PercentRefine' , PercentRefine)
       call read_var('MaxTotalBlock',  MaxTotalBlock) 
       call read_var('DiffCriteriaLevel',  DeltaCritera)
       DoSortAmrCrit = PercentCoarsen > 0.0 .or. PercentRefine > 0.0
    case default
       call stop_mpi(NameSub//'incorect PARAM.in!')
    end select
  end subroutine read_amr_criteria

  !============================================================================

  subroutine clean_amr_criteria

    if(.not.allocated(CoarsenCrit_I)) RETURN
    deallocate(CoarsenCrit_I, RefineCrit_I, iVarCrit_I)
    deallocate(CoarsenCritAll_I, RefineCritAll_I)
    if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
    if(allocated(iNode_I)) deallocate(iNode_I)
    nAmrCrit     = 0
    nAmrCritUsed = 0
    nAmrCritOld  = 0
    nIntCrit     = 0
    nBlockOld = 0

  end subroutine clean_amr_criteria

  !============================================================================
  subroutine test_amr_criteria()
    use BATL_size, ONLY : MaxBlock,nBlock, iRatio, jRatio, kRatio
    use BATL_mpi, ONLY: iProc, nProc
    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
         refine_tree_node, distribute_tree, clean_tree, Unused_B, iNode_B, &
         nNode, show_tree, iStatusNew_A, &
         Coarsen_, Unset_, adapt_tree, move_tree, distribute_tree
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid
    use BATL_geometry, ONLY: init_geometry
    use BATL_amr, ONLY: do_amr, init_amr
    use BATL_size, ONLY: MaxDim, nDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
    ! For Random generation
    integer :: jSeed
    logical :: IsFirst
    integer, parameter :: iMPLIER=16807, &
         iMODLUS=2147483647, &
         iMOBYMP=127773, &
         iMOMDMP=2836

    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .false.
    real, parameter:: DomainMin_D(MaxDim) = (/ -24.0, -24.0, -24.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/ 24.0, 24.0, 24.0 /)
    integer :: iNode, iBlock

    real, allocatable :: Criterias_IB(:,:),AllCriterias_IBP(:,:,:)
    real, allocatable :: PreCriterias_IB(:,:)
    real, allocatable :: RefineLevel_I(:), CoursenLevel_I(:)

    integer :: nCritExt = 4
    integer, allocatable :: iA_I(:)
    real, allocatable :: TestState_VGB(:,:,:,:,:)
    integer :: nVar =3
    integer :: i,j,k,iVar
    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_amr_criteria'
    !-----------------------------------------------------------------------
    DoTestMe = iProc == 0

    if(DoTestMe) write(*,*) 'Starting ',NameSub

    call init_tree(MaxBlockTest)
    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
    call set_tree_root( nRootTest_D(1:nDim))

    call find_tree_node( (/0.5,0.5,0.5/), iNode)
    call refine_tree_node(iNode)
    call distribute_tree(.true.)
    call create_grid
    call init_amr

    call srand(123456789+iProc)

    allocate(Criterias_IB(nCritExt,nBlock), &
         AllCriterias_IBP(nCritExt,nBlock,nProc),&
         PreCriterias_IB(nCritExt,nBlock))
    allocate(RefineLevel_I(nCritExt),CoursenLevel_I(nCritExt))


    !===================== Begin Test  =======================

    !--------------------- internal --------------------------
    nIntCrit = 1

    allocate(CoarsenCrit_I(nVar), &
         RefineCrit_I(nVar),iVarCrit_I(nVar))

    allocate(TestState_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    CoarsenCrit_I = -1.0
    RefineCrit_I  =  1.0
    TestState_VGB = 1.0
    iVarCrit_I(nIntCrit)=1

    CoarsenCrit_I = -1.0
    RefineCrit_I  = 1.0

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK
          do j = MinJ, MaxJ
             do i = MinI,MaxI
                do iVar=1,nVar
                   TestState_VGB(iVar,i,j,k,iBlock) = &
                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
                end do
             end do
          end do
       end do
    end do


    call set_amr_criteria(nVar,TestState_VGB)

!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       write(*,'(i6,f16.12)') iNode_B(iBlock),AmrCrit_IB(1,iBlock)
!!$    end do
!!$
!!$    do iNode = 1, nNode-1
!!$       print *,iNode, iRank_A(iNode)
!!$    end do

    do iNode = 2, nNode-1
       if(iRank_A(iNode-1) > iRank_A(iNode)) then
          write(*,*) " ERROR in ",NameSub, " Internal test"
       end if
    end do

    !-------------------- external --------------------------
    Criterias_IB = 0.0
    RefineLevel_I =  1.0
    CoursenLevel_I = -1.0
    nExtCrit = 4 
    nCritExt = 1
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) then
          Criterias_IB(1,iBlock) = 10.0
       else
          Criterias_IB(1,iBlock) = AmrCrit_IB(1,iBlock)
       end if
    end do

    call set_amr_criteria(nVar,TestState_VGB, &
         nCritExt, Criterias_IB, CoursenLevel_I, RefineLevel_I)

    do iNode = 2, nNode-1
       if(iRank_A(iNode-1)>iRank_A(iNode)) then
          write(*,*) " ERROR in ",NameSub, "Externa=Intenal test"
       end if
    end do


    !-------------------- internal x 2 -------------------

    if(iRatio > 1 .and. jRatio > 1 .and. kRatio >1 ) then

       nCritExt = 0
       nExtCrit = 4 
       nIntCrit = 2
       iVarCrit_I(1:nIntCrit)=(/ 1,2 /)
       do iBlock=1,nBlock
          if(Unused_B(iBlock)) CYCLE
          do k = MinK, MaxK
             do j = MinJ, MaxJ
                do i = MinI,MaxI
                   TestState_VGB(1,i,j,k,iBlock) = &
                        dexp(0.1*(i*(35-iNode_B(iBlock)+1)))

                   TestState_VGB(2,i,j,k,iBlock) = &
                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
                end do
             end do
          end do
       end do

       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) then
             Criterias_IB(1,iBlock) = 10.0
          else
             Criterias_IB(1,iBlock) = AmrCrit_IB(1,nBlock-iBlock+1)
          end if

       end do


       call set_amr_criteria(nVar,TestState_VGB, &
            nCritExt, Criterias_IB, CoursenLevel_I, RefineLevel_I)

       allocate(iA_I(nNode-1))

       iA_I =(/ 1,35, 2, 34, 33,  3,  4, 32, 31,  5,  6, 30, 29,  7, &
            8, 28, 27,  9, 26, 10, 25, 11, 12, 24, 13, 23, 22, 15, 21, &
            16, 17, 20, 19, 18 /)

!!$       do iBlock = 1, nBlock
!!$          if(Unused_B(iBlock)) CYCLE
!!$          print *,"Criterias_IB(1,iBlock) = ", &
!!$               AmrCrit_IB(1:2,iBlock), iNode_B(iBlock)
!!$       end do
!!$
!!$       do iNode = 1, nNode-1
!!$          print *,iRank_A(iNode)," :: ", iA_I(iNode) 
!!$       end do

       do iNode = 1, nNode-1, 2
          if(iRank_A(iNode) /= iA_I(iNode)) &
               write(*,*) " ERROR in ",NameSub, "2 x Intenal test"
       end do

       deallocate(iA_I)
    end if
    !-------------------- testing levels ---------------------
    !intenals
    CoarsenCrit_I = -1.0
    RefineCrit_I  =  1.0
    !externals
    RefineLevel_I =  0.030
    CoursenLevel_I = 0.020

    PercentRefine  = 30.0
    PercentCoarsen = 10.0
    nExtCrit = 4 
    nCritExt = 1
    nIntCrit = 1

    ! init Internals
    iVarCrit_I(nIntCrit) = 1
    do iBlock=1,nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK
          do j = MinJ, MaxJ
             do i = MinI,MaxI
                !do iVar=1,nVar
                !print *,rand()
                TestState_VGB(2,i,j,k,iBlock) = &
                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
                !end do
             end do
          end do
       end do
    end do

    !init Externals
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       Criterias_IB(1,iBlock) = 0.001*(nNode - iNode_B(iBlock)+1)
    end do

    call set_amr_criteria(nVar,TestState_VGB, &
         nCritExt, Criterias_IB, CoursenLevel_I, RefineLevel_I)

    do k = nNodeCoarsen+1, nNode-1-nNodeRefine
       if( iRank_A(k) < nNodeRefine .and. & 
            iRank_A(k) > nNode-1 - nNodeCoarsen ) &
            write(*,*) "Error in seting refinment by criteria"
    end do

    ! test unused block
    Criterias_IB = 0.0
    RefineLevel_I =  1.0
    CoursenLevel_I = -1.0
    nExtCrit = 4 
    nCritExt = 0
    nIntCrit = 1
    iStatusNew_A = Unset_

    do iBlock=1,nBlock
       if(iNode_B(iBlock) > 3**nDim) then 
          iStatusNew_A(iNode_B(iBlock)) =  Coarsen_
       end if
    end do
    !call show_tree(NameSub,.true.)
    call adapt_tree
    call distribute_tree(DoMove=.false.)
    call do_amr(nVar,TestState_VGB)
    call move_tree
    !call show_tree(NameSub,.true.)

    do iBlock=1,nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK
          do j = MinJ, MaxJ
             do i = MinI,MaxI
                do iVar=1,nVar
                   TestState_VGB(iVar,i,j,k,iBlock) = &
                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
                end do
             end do
          end do
       end do
    end do

    call set_amr_criteria(nVar,TestState_VGB)

    do iNode = 2, nNode
       if(iRank_A(iNode-1)>iRank_A(iNode)) then
          write(*,*) " ERROR in ",NameSub, " unused  test"
       end if
    end do

    !===================== End Test  =======================
    deallocate(CoarsenCrit_I,RefineCrit_I,iVarCrit_I)
    deallocate(TestState_VGB)

    deallocate(Criterias_IB,PreCriterias_IB,AllCriterias_IBP)
    deallocate(RefineLevel_I,CoursenLevel_I)
    call clean_grid
    call clean_tree

  contains

    ! The saudo random number generator is for testing performense in
    ! parallel sorting
    subroutine srand(iSeed)
      integer, intent(in) :: iSeed
      jSeed = iSeed
      IsFirst = .true.
    end subroutine srand

    real function rand()
      !  A pseudo-random number generator implemented to make sure that 
      !  all platform reproduce the same sequence for testing and compering.
      !  The algorithm is based on "Integer Version 2" given in :
      !
      !       Park, Steven K. and Miller, Keith W., "Random Number Generators: 
      !       Good Ones are Hard to Find", Communications of the ACM, 
      !       October, 1988.

      integer :: nHvalue,nLvalue,nTestv
      integer, save :: nExtn

      if(IsFirst) then
         nExtn=jSeed
         IsFirst = .false.
      end if

      nHvalue = nExtn/iMOBYMP
      nLvalue = mod(nExtn,iMOBYMP)
      nTestv = iMPLIER*nLvalue - iMOMDMP*nHvalue
      if(nTestv > 0) then
         nExtn = nTestv
      else
         nExtn = nTestv + iMODLUS
      end if

      rand = real(nExtn)/real(iMODLUS)

    end function rand

  end subroutine test_amr_criteria

  end module BATL_amr_criteria
