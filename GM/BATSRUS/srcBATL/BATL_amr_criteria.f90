!^CFG COPYRIGHT UM

module BATL_amr_criteria

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

  implicit none

  SAVE

  private ! except

  public set_amr_criteria
  public clean_amr_criteria
  public read_amr_criteria_param
  public test_histogram_amr_criteria

  integer, public            :: nAmrCrit
  real, public, allocatable  :: AmrCrit_IB(:,:)

  ! Threshold limits for refine or unrefined the grid (length nCrit) 
  real, allocatable, dimension(:)    :: CoarsenCrit_I, RefineCrit_I, &
       CoarsenAmrCrit_I, RefineCritAll_I, GlobalCritMaxAll_I, CritMaxAll_I 

  integer, allocatable:: iVarCrit_I(:) ! Index to variables
  ! used in criteria
  real :: cAmrWavefilter
  real, parameter :: cEpsilon = 1.0d-16 ! avoid zero in denominator
  integer :: nCrit ! Number of variables used for amr internal decisions
  integer :: nAmrCritOld, nBlockOld

  ! How large the relativ descrepensy in the Critera for symerti points can be
  ! and stil be refined in the same way and the same for geometric variation
  real, parameter :: DeltaSymCrit = 1.0e-2
  real, parameter :: DeltaSymGeo  = 1.0e-6

  real ::  precentRefine=0.0, percentCoarsen=0.0
contains
  !============================================================================

  subroutine histogram_amr_criteria(nCrit,Crit_IB,CoarsenCrit_I,RefineCrit_I)

    ! the rutine will find the candidate for refinment or coarsing based
    ! upon the prosented of refinment and the prosented wanted refined and
    ! and coarsend. It will ONLY modify the criteria the values in the criteia 
    ! list itself. The refinemnet ciriteria will overrule the prosented refinment.


    use BATL_mpi, ONLY: iComm, iProc, nProc
    use ModMpi, ONLY: MPI_Allgather, MPI_REAL, MPI_INTEGER, MPI_SUM
    use BATL_size, ONLY: nBlock
    use ModSort, ONLY: sort_quick
    use BATL_tree, ONLY: Unused_B
    ! the last criteria is used for finding symetic
    ! points 
    integer, intent(in) :: nCrit
    real, intent(inout) :: Crit_IB(:,:)
    real, intent(in) :: CoarsenCrit_I(:), RefineCrit_I(:)

    integer :: iError
    real, allocatable :: AllCrit_IBP(:,:,:)
    ! SortC contaons all Crieteas for all blocks
    ! SortIdx_II gives the new indexes after sorting
    integer, allocatable :: SortIdx_II(:,:)
    ! convert from running index to Proc (1) and Block(2) index
    integer, allocatable :: Idx_II(:,:)
    real, allocatable :: SortCrit_II(:,:)
    real, allocatable :: MaxRestCrit(:,:)
    integer iCrit, iBlock, iProces, k,l
    ! idxSC : index of in sorted order of the Criteia
    ! idxSG : index of in sorted order of the Gemoetric (nCrit)
    integer :: nSort, idxSC, idxSG, nCont
    integer :: nDesired, nRefined, nCoarsened, currentBlocks
    real :: MinDif, diff, minValue, LocMaxCrit
    integer :: iMinCrit, iMaxCrit, LocMaxIdx, TotBlocks, MaxNBlock
    integer, allocatable :: AllBlocks(:),RecivCont(:),RecivDisp(:)
    integer, parameter :: one_=1
    !-----------------------------------------------------------------------
    allocate(AllBlocks(nProc),RecivCont(nProc),RecivDisp(nProc))
    AllBlocks(iProc+1) = nBlock
    call MPI_Allgather(nBlock, one_,MPI_INTEGER, AllBlocks,one_,MPI_INTEGER,IComm,iError)
 !   print *,"AllBlocks = ",AllBlocks
    MaxNBlock = maxval(AllBlocks)
    TotBlocks = sum(AllBlocks)

    allocate(AllCrit_IBP(nCrit,MaxNBlock,nProc))
    RecivCont(:) = AllBlocks(:)*nCrit
    RecivDisp(1) = 0
    do iProces=2,nProc 
       RecivDisp(iProces) = RecivDisp(iProces-1)+MaxNBlock*nCrit
    end do

!    print *, "RecivDisp = ",RecivDisp

    AllCrit_IBP = -77
    call MPI_allgatherv(Crit_IB(:,:),nCrit*nBlock,MPI_REAL, &
         AllCrit_IBP(:,:,:),RecivCont(:),RecivDisp(:),MPI_REAL,iComm,iError)

    !    print *, " AllCrit_IBP = ", AllCrit_IBP 

!!$    allocate(SortCrit_I((1+nBlock*nProc*nCrit)),SortIdx_I((1+nBlock*nProc*nCrit)))
!!$    allocate(Idx_II(nBlock*nProc*nCrit,3))
!!$
!!$    print *," SortIdx_II ", (SortIdx_II)
!!$
!!$    Idx_II = 0
!!$    SortIdx_I  = 0
!!$    SortCrit_I  = 0.0
!!$
!!$    call MPI_allgather(Crit_IB(:,:),nCrit*nBlock,MPI_REAL, &
!!$         AllCrit_IBP(:,:,:),nCrit*nBlock,MPI_REAL,iComm,iError)
!!$
!!$    k = 0
!!$    do iProces = 1,nProc
!!$       do iBlock = 1, nBlock
!!$          if(Unused_B(iBlock)) CYCLE
!!$          do iCrit=1,nCrit-1
!!$             k = k +1
!!$             SortCrit_I(k) = AllCrit_IBP(iCrit,iBlock,iProces)
!!$             SortIdx_I(k) = k
!!$             Idx_I(k,1) = iCrit
!!$             Idx_I(k,1) = iProces
!!$             Idx_I(k,2) = iBlock
!!$          end do
!!$       end do
!!$    end do
!!$    nSort = k
!!$
!!$    call sort_quick(nSort, SortCrit_I(1:nSort), SortIdx_I(1:nSort))
!!$
!!$    ! store the largetst value in the array
!!$    do iCrit=1,nCrit-1
!!$       SortCrit_II(nSort+1,iCrit) = SortCrit_II(SortIdx_II(nSort,iCrit),iCrit)
!!$       SortIdx_II(nSort+1,iCrit) = SortIdx_II(SortIdx_II(nSort,iCrit),iCrit)
!!$    end do
!!$
!!$    ! Test if the Criteias is normeixed, if not normize them
!!$    do iCrit=1,nCrit-1
!!$       if(SortCrit_II(SortIdx_II(1,iCrit),iCrit) /= 0.0 .or. SortCrit_II(SortIdx_II(nSort,iCrit),iCrit) /= 1.0) then
!!$          minValue = SortCrit_II(SortIdx_II(1,iCrit),iCrit)
!!$          diff = SortCrit_II(SortIdx_II(nSort,iCrit),iCrit) - minValue
!!$          if(diff /= 0.0) then
!!$             SortCrit_II(:,iCrit) = (SortCrit_II(:,iCrit) - minValue)/diff
!!$             AllCrit_IBP(iCrit,:,:) = (AllCrit_IBP(iCrit,:,:) -  minValue)/diff
!!$          end if
!!$       end if
!!$    end do
!!$
!!$    do iCrit=1,nCrit
!!$       print *,""
!!$       print *,"Sorted Criterias[",iCrit,"]",SortCrit_II(1:nSort,iCrit) 
!!$       print *,"Sorted Indexes[",iCrit,"]",SortIdx_II(1:nSort,iCrit)
!!$       print *,""
!!$    end do


    !===========================================================
    !=====
    !=====   THE OLD WAY
    !=====
    !===========================================================


    allocate(SortCrit_II((1+TotBlocks),(1+nCrit)),SortIdx_II((1+TotBlocks),(1+nCrit)))
    allocate(Idx_II(TotBlocks ,2))
    allocate(MaxRestCrit(2,nCrit-1))

!!$    call MPI_allgather(Crit_IB(:,:),nCrit*nBlock,MPI_REAL, &
!!$         AllCrit_IBP(:,:,:),nCrit*nBlock,MPI_REAL,iComm,iError)

    Idx_II = 0
    SortIdx_II  = 0
    SortCrit_II  = 0.0

    k = 0
    do iProces = 1,nProc
       do iBlock = 1, AllBlocks(iProces)
          !if(Unused_B(iBlock)) CYCLE
          k = k +1
          do iCrit=1,nCrit-1
             SortCrit_II(k,iCrit) = AllCrit_IBP(iCrit,iBlock,iProces)
             SortIdx_II(k, iCrit) = k
          end do
          Idx_II(k,1) = iProces
          Idx_II(k,2) = iBlock
       end do
    end do
    nSort = k
    ! Store Max value

    do iCrit=1,nCrit
       call sort_quick(nSort, SortCrit_II(1:nSort,iCrit), SortIdx_II(1:nSort,iCrit))
    end do

    ! store the largetst value in the array
    do iCrit=1,nCrit-1
       SortCrit_II(nSort+1,iCrit) = SortCrit_II(SortIdx_II(nSort,iCrit),iCrit)
       SortIdx_II(nSort+1,iCrit) = SortIdx_II(SortIdx_II(nSort,iCrit),iCrit)
    end do

    ! Test if the Criteias is normeixed, if not normize them
    do iCrit=1,nCrit-1
       if(SortCrit_II(SortIdx_II(1,iCrit),iCrit) /= 0.0 .or. SortCrit_II(SortIdx_II(nSort,iCrit),iCrit) /= 1.0) then
          minValue = SortCrit_II(SortIdx_II(1,iCrit),iCrit)
          diff = SortCrit_II(SortIdx_II(nSort,iCrit),iCrit) - minValue
          if(diff /= 0.0) then
             SortCrit_II(:,iCrit) = (SortCrit_II(:,iCrit) - minValue)/diff
             AllCrit_IBP(iCrit,:,:) = (AllCrit_IBP(iCrit,:,:) -  minValue)/diff
          end if
       end if
    end do

    !Crit_IB(:,:) = AllCrit_IBP(:,1:nBlock,iProc+1)
    !   do iCrit=1,nCrit
    !      print *,""
    !      print *,"Sorted Criterias[",iCrit,"]",SortCrit_II(1:nSort,iCrit) 
    !      print *,"Sorted Indexes[",iCrit,"]",SortIdx_II(1:nSort,iCrit)
    !      print *,""
    !   end do

    !-------------------------- REFINNG --------------------------
    nRefined =0
    do iCrit=1,nCrit-1
       do k = nSort,1,-1
          idxSC = SortIdx_II(k,iCrit)
          !print *, "k = ",k, SortCrit_II(idxSC,iCrit), RefineCrit_I(iCrit)
          if(SortCrit_II(idxSC,iCrit) < RefineCrit_I(iCrit)) EXIT
          SortIdx_II(idxSC,nCrit+1) = 1.0
          AllCrit_IBP(1:nCrit-1,Idx_II(idxSC,2),Idx_II(idxSC,1)) = &
               RefineCrit_I(1:nCrit-1)*(1.0+DeltaSymCrit)
       end do
       !nRefined = max(nRefined,nSort-k)
    end do

    nRefined = sum(SortIdx_II(1:nsort,nCrit+1))
    nDesired = int(precentRefine*sum(AllBlocks)/100.0)

    !    print *, "Number  Refinment by criteria : ",nRefined, " and we want : ",nDesired

    ! we want to refine more
    if(nDesired > nRefined) then
       !find whitch next criterea is closest to the refinement criteria

       ! make a list of the Criteria values that was next in the list for
       ! refinement for each of the criterias
       do iCrit=1,nCrit-1
          do k = nSort,1,-1
             idxSC = SortIdx_II(k,iCrit)
             !print *, "k = ",k, SortCrit_II(idxSC,iCrit), RefineCrit_I(iCrit)
             if(SortIdx_II(idxSC,nCrit+1) == 1.0 ) CYCLE
             MaxRestCrit(1,iCrit) = SortCrit_II(idxSC,iCrit)
             MaxRestCrit(2,iCrit) = k
             EXIT
          end do
       end do


       do k=1,nDesired-nRefined
          LocMaxCrit = MaxRestCrit(1,1)
          LocMaxIdx  = MaxRestCrit(2,1)
          iMaxCrit   = 1
          !Find the index and value of the next block to refine
          ! by findec the one with max criteria value
          do iCrit=1, nCrit-1
             if(LocMaxCrit < MaxRestCrit(1,iCrit)) then
                LocMaxCrit = MaxRestCrit(1,iCrit)
                LocMaxIdx = int(MaxRestCrit(2,iCrit))
                iMaxCrit = iCrit
             end if
          end do
          !          print *,"local max", MaxRestCrit(1,:), iMaxCrit, LocMaxIdx,SortIdx_II(LocMaxIdx,iMaxCrit)
          ! Mark the block for refinment by seting its value 
          ! larger then the refinement criteria
          IdxSC = SortIdx_II(int(MaxRestCrit(2,iMaxCrit)),iMaxCrit)
          SortIdx_II(idxSC,nCrit+1) = 1.0
          AllCrit_IBP(1:nCrit-1,Idx_II(idxSC,2),Idx_II(idxSC,1)) = &
               RefineCrit_I(1:nCrit-1)*(1.0+DeltaSymCrit)

          ! Pick the next value from the sorted criteria list 
          ! for the one block that was refined
          do iCrit=1, nCrit-1
             IdxSC = SortIdx_II(int(MaxRestCrit(2,iCrit)),iCrit)
             if(SortIdx_II(idxSC,nCrit+1) /= 1.0 ) CYCLE
             MaxRestCrit(2,iCrit) = MaxRestCrit(2,iCrit) -1
          end do
          IdxSC = SortIdx_II(int(MaxRestCrit(2,iMaxCrit)),iMaxCrit)
          MaxRestCrit(1,iMaxCrit) = SortCrit_II(IdxSC,iMaxCrit)
       end do
    end if

    !-------------------------- COURSNING --------------------------
    k = 0
    do iProces = 1,nProc
       do iBlock = 1,  AllBlocks(iProces)
!          !if(Unused_B(iBlock)) CYCLE
          k = k +1
          SortIdx_II(k,:) = k
          SortCrit_II(k,nCrit+1) = 0.0
          do iCrit=1,nCrit-1
             SortCrit_II(k,iCrit) = AllCrit_IBP(iCrit,iBlock,iProces)
             SortCrit_II(k,nCrit+1) = SortCrit_II(k,nCrit+1) + &
                  AllCrit_IBP(iCrit,iBlock,iProces)/abs(SortCrit_II(nSort+1,iCrit))
          end do
       end do
    end do

    SortCrit_II(:,nCrit+1)= SortCrit_II(:,nCrit+1)/(nCrit-1)
    call sort_quick(nSort, SortCrit_II(1:nSort,nCrit+1), SortIdx_II(1:nSort,nCrit+1))

    !    print *,'Sorted coarsning' , SortCrit_II(SortIdx_II(1:nSort,nCrit+1),nCrit+1)
!!$    print *,""
!!$    print *,"Sorted for coarsning Criterias[",iCrit,"]",SortCrit_II(iCrit,1:nSort) 
!!$    print *,"Sorted for coarsning Indexes[",iCrit,"]",SortIdx_II(iCrit,1:nSort)
!!$    print *,""

    nDesired = int(percentCoarsen*sum(AllBlocks)/100.0)


    ! set all criteria values that is smaller then CoarsenCrit_I 
    ! to a value larger then the limit so it will not be set for 
    ! coursning
    do iCrit=1,nCrit-1
       do k = 1,nSort
          idxSC = SortIdx_II(k,iCrit)
          if(SortCrit_II(idxSC,iCrit) >= CoarsenCrit_I(iCrit)) EXIT
          AllCrit_IBP(iCrit,Idx_II(idxSC,2),Idx_II(idxSC,1)) = &
               CoarsenCrit_I(iCrit)*(1.0+DeltaSymCrit)
       end do
    end do

    do k = 1,nSort
       idxSC = SortIdx_II(k,nCrit+1)
       if(SortCrit_II(idxSC,nCrit+1) >= CoarsenCrit_I(1)) EXIT
       do iCrit=1,nCrit-1
          AllCrit_IBP(iCrit,Idx_II(idxSC,2),Idx_II(idxSC,1)) = &
               CoarsenCrit_I(iCrit)*(1.0-DeltaSymCrit)
       end do
    end do



    ! set criteria vlue of all  nDesired blocks with smaalest values given by
    ! SortIdx_II(nCrit+1,k) to a vlaue under there CoarsenCrit_I 
    do k =1, nDesired
       idxSC = SortIdx_II(k,nCrit+1)
       do iCrit=1,nCrit-1
          AllCrit_IBP(iCrit,Idx_II(idxSC,2),Idx_II(idxSC,1)) = &
               CoarsenCrit_I(iCrit)*(1.0-DeltaSymCrit)
          !print *,"Corsening k = ",k, SortCrit_II(iCrit,idxSC),CoarsenCrit_I(iCrit)
       end do
    end do

    !   print *, " Number Coarsend by criteria : ",nCoarsened, " and we want : ",nDesired

    Crit_IB(:,:) = AllCrit_IBP(:,1:nBlock,iProc+1)

    deallocate(AllCrit_IBP)
    deallocate(SortCrit_II,SortIdx_II)
    deallocate(Idx_II)
    deallocate(MaxRestCrit)
    deallocate(AllBlocks,RecivCont,RecivDisp)
!!$  contains
!!$
!!$    subroutine set_criteria(crit, blk, proc, value) 
!!$      if(iproc == proc) &
!!$           Crit_IB(cit,blk) = value
!!$    end subroutine set_criteria

  end subroutine histogram_amr_criteria

  !============================================================================

  subroutine test_histogram_amr_criteria()

    use BATL_size, ONLY : MaxBlock,nBlock
    use BATL_mpi, ONLY: iComm, iProc, nProc
    use ModMpi, ONLY: MPI_allreduce, MPI_INTEGER, MPI_SUM

    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
         refine_tree_node, distribute_tree, clean_tree, Unused_B, iNode_B
    use BATL_grid, ONLY: init_grid, create_grid, clean_grid
    use BATL_geometry, ONLY: init_geometry

    use BATL_size, ONLY: MaxDim,  nDim

    ! For Random generation
    integer :: jseed, seed
    logical :: isFirst
    integer, parameter :: MPLIER=16807, &
         MODLUS=2147483647, &
         MOBYMP=127773, &
         MOMDMP=2836

    real, parameter :: pi=3.14159


    integer, parameter:: MaxBlockTest            = 50
    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
    logical, parameter:: IsPeriodicTest_D(MaxDim)= .false.
    real, parameter:: DomainMin_D(MaxDim) = (/ -24.0, -24.0, -24.0 /)
    real, parameter:: DomainMax_D(MaxDim) = (/ 24.0, 24.0, 24.0 /)
    integer :: iNode, iBlock, iCrit, n

    real, allocatable :: Criterias_IB(:,:),AllCriterias_IBP(:,:,:)
    real, allocatable :: PreCriterias_IB(:,:)
    real, allocatable :: RefineLevel_I(:), CoursenLevel_I(:)
    ! As if time of writing 4 April 2011 this is ment as reproduing
    ! the behavior of BATSRUS amr_physics so we only have 3 proper 
    ! criterias and the fourth is radial behaviour for symetri risens
    integer, parameter :: nCrit = 4, one_ = 1
    integer :: iProces, iError, TotBlocks ! Total number of used blocks
    integer :: TotRefine, nRefine, TotCoarsen, nCoarsen
    integer :: CoarsRefine(2), TotCoarsRefine(2)
    real :: Criteia
    integer, allocatable :: AllBlocks(:)
    logical:: DoTestMe
    character(len=*), parameter :: NameSub = 'test_histogram_amr_criteria'
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

    call srand(123456789+iProc)

    allocate(Criterias_IB(nCrit,nBlock), &
         AllCriterias_IBP(nCrit,nBlock,nProc),&
         PreCriterias_IB(nCrit,nBlock))
    allocate(RefineLevel_I(nCrit),CoursenLevel_I(nCrit))


    allocate(AllBlocks(nProc))
    AllBlocks(iProc+1) = nBlock
    call MPI_Allgather(nBlock, one_,MPI_INTEGER, AllBlocks,one_,MPI_INTEGER,IComm,iError)

!!$    RefineLevel_I  =  0.99
!!$    CoursenLevel_I =  0.05
!!$    precentRefine=35.0
!!$    percentCoarsen=35.0

    TotBlocks = 0
    do iProces = 1,nProc 
       do iBlock = 1, AllBlocks(iProces)
          if(Unused_B(iBlock)) CYCLE
          TotBlocks = TotBlocks +1
       end do
    end do

!!$    print *,""
!!$    print *," Before sorting"  
!!$    do iBlock = 1, nBlock
!!$       n = (iNode_B(iBlock)-1)
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do iCrit = 1, nCrit
!!$          select case(iCrit)
!!$          case(1)
!!$             Criterias_IB(iCrit,iBlock) = 1.0 - 1.0*n/(nProc*nBlock)
!!$          case(2)
!!$             Criterias_IB(iCrit,iBlock) = sin(2.0*pi*n/(nProc*nBlock))
!!$          case default
!!$             Criterias_IB(iCrit,iBlock) = rand()
!!$          end select
!!$       end do
!!$       n=n+1
!!$       print *, Criterias_IB(:,iBlock)
!!$    end do

    !===================== Test unchenged  =======================

    Criterias_IB = 0.0
    PreCriterias_IB = 0.0
    !    print *,""
    !    print *," Before sorting"  
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do iCrit=1,nCrit
          Criterias_IB(iCrit,iBlock) = rand()
       end do
       !print *, Criterias_IB(:,iBlock)
    end do
    !make sure the data set is normelized
    if(nProc == 1) then
       Criterias_IB(:,1) = 0.0
       Criterias_IB(:,2) = 1.0
    else if (iProc == 0) then
       Criterias_IB(:,1) = 0.0
    else if (iProc == 1) then
       Criterias_IB(:,1) = 1.0
    end if

    PreCriterias_IB(:,:) = Criterias_IB(:,:)

    RefineLevel_I  =  2.0
    CoursenLevel_I =  -2.0
    precentRefine  = 0.0
    percentCoarsen = 0.0

    call histogram_amr_criteria(nCrit,Criterias_IB,CoursenLevel_I,RefineLevel_I)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do iCrit=1,nCrit
          if (Criterias_IB(iCrit,iBlock) - PreCriterias_IB(iCrit,iBlock) /= 0.0) &
               write(*,*) NameSub, " Error, Criterias_IB has chenged, diff",&
               Criterias_IB(iCrit,iBlock),  PreCriterias_IB(iCrit,iBlock)
       end do
     !  print *, Criterias_IB(:,iBlock), " : ", PreCriterias_IB(:,iBlock)
    end do

    !===================== Test Procentage =======================
    Criterias_IB = 0.0
    !    print *,""
    !    print *," Before sorting"  
    do iBlock = 1, nBlock
       n = (iNode_B(iBlock)-1)
       if(Unused_B(iBlock)) CYCLE
       Criterias_IB(1:nCrit-1,iBlock) = 1.0 - 1.0*n/TotBlocks
       !       print *, Criterias_IB(:,iBlock)
    end do

    RefineLevel_I  =  2
    CoursenLevel_I =  -1
    precentRefine=35.0
    percentCoarsen=35.0

    call histogram_amr_criteria(nCrit,Criterias_IB,CoursenLevel_I,RefineLevel_I)

    nRefine = 0
    nCoarsen = 0
    do iBlock=1,nBlock
       if(Criterias_IB(1,iBlock) > 1.0) &
            nRefine = nRefine+1
       if(Criterias_IB(1,iBlock) < 0.0) &
            nCoarsen = nCoarsen+1
    end do
    CoarsRefine = (/nCoarsen,nRefine /) 
    call MPI_allreduce(CoarsRefine(1),TotCoarsRefine(1),2, MPI_INTEGER,&
         MPI_SUM,iComm, iError)

    if(TotCoarsRefine(2) /= int(sum(AllBlocks)*0.01*precentRefine)) &
         write(*,*) NameSub, " Error, diff in procent refine : ", &
         TotCoarsRefine(2) ,"  should be ", int(TotBlocks*0.01*precentRefine)
    if(TotCoarsRefine(1) /= int(sum(AllBlocks)*0.01*percentCoarsen)) &
         write(*,*) NameSub, " Error, diff in procent coarsen : ", &
         TotCoarsRefine(1),"  should be ", int(TotBlocks*0.01*percentCoarsen)

!!$    print *,""
!!$    print *," After sorting"  
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       print *, Criterias_IB(:,iBlock)
!!$    end do

    !===================== Test Level  =======================
    Criterias_IB = 0.0
    !    print *,""
    !    print *," Before sorting"
    TotBlocks =  sum(AllBlocks(:))
    !print *,"sum(AllBlocks(1:iProc))", sum(AllBlocks(1:iProc))
    do iBlock = 1, nBlock
       n = sum(AllBlocks(1:iProc)) + iBlock ! (iNode_B(iBlock)-1)
       if(Unused_B(iBlock)) CYCLE
       Criterias_IB(1:nCrit-1,iBlock) = 1.0 - 1.0*(n-1)/(TotBlocks-1)
       !print *, Criterias_IB(:,iBlock)
    end do


    !print *,Criterias_IB(1,:)

    RefineLevel_I  = 0.9
    CoursenLevel_I = 0.1
    precentRefine  = 0.0
    percentCoarsen = 0.0

    call histogram_amr_criteria(nCrit,Criterias_IB,CoursenLevel_I,RefineLevel_I)

    do iBlock = 1, nBlock
        n = sum(AllBlocks(1:iProc)) + iBlock!n = (iNode_B(iBlock)-1)
       if(Unused_B(iBlock)) CYCLE

       Criteia = 1.0 - 1.0*(n-1)/(TotBlocks-1)

       if(Criteia  >RefineLevel_I(1)) then
          if(Criterias_IB(1,iBlock) < RefineLevel_I(1)) then
             if(abs(Criterias_IB(1,iBlock)-RefineLevel_I(1)) > 1e-15)&
                  write(*,*) NameSub, " Error in refining for level refinment",&
                  "Criteria = ",Criterias_IB(1,iBlock) ,"Refinig level ", RefineLevel_I(1)
          end if
       end if

       if(Criteia < CoursenLevel_I(1) ) then
          if( Criterias_IB(1,iBlock) > CoursenLevel_I(1) ) then
             if(abs(Criterias_IB(1,iBlock)-CoursenLevel_I(1))> 1e-15) &
                  write(*,*) NameSub, " Error in coarsning for level refinment" ,&
                  "Criteria = ",Criterias_IB(1,iBlock) ," Coarsen level ",  CoursenLevel_I(1),&
                  Criterias_IB(1,iBlock)-CoursenLevel_I(1), iBlock, iProc
          end if
       end if
    end do

!!$    print *,""
!!$    print *," After sorting"  
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       print *, Criterias_IB(:,iBlock)
!!$    end do



    deallocate(Criterias_IB,PreCriterias_IB,AllCriterias_IBP)
    deallocate(RefineLevel_I,CoursenLevel_I)
    call clean_grid
    call clean_tree

  contains

    subroutine srand(iseed)
      integer, intent(in) :: iseed
      jseed = iseed
      isFirst = .true.
    end subroutine srand

    real function rand()
      !  A pseudo-random number generatur impelemented to make sure that 
      !  all platform reproduse the same sequnse for testing and comersions.
      !  The algorithm is based on "Integer Version 2" given in :
      !
      !       Park, Steven K. and Miller, Keith W., "Random Number Generators: 
      !       Good Ones are Hard to Find", Communications of the ACM, 
      !       October, 1988.

      integer hvlue,lvlue,testv
      integer, save :: nextn

      if(isFirst) then
         nextn=jseed
         isFirst = .false.
      end if

      hvlue = nextn/MOBYMP
      lvlue = mod(nextn,MOBYMP)
      testv = MPLIER*lvlue - MOMDMP*hvlue
      if(testv > 0) then
         nextn = testv
      else
         nextn = testv + MODLUS
      end if

      rand = real(nextn)/real(MODLUS)

    end function rand

  end subroutine test_histogram_amr_criteria




  subroutine set_amr_criteria(nVar, State_VGB, & 
       nCritExt, CritExt_IB, CoarsenCritExt_I, RefineCritExt_I)

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
	 nI, nJ, nK, MaxDim, nDim, MaxBlock, nBlock
    use BATL_tree, ONLY: iStatusNew_A, Refine_, Coarsen_, &
         Unused_B, iNode_B
    use BATL_mpi, ONLY: iComm, nProc
    use ModMpi, ONLY: MPI_allreduce, MPI_REAL, MPI_MAX

    integer,  intent(in)   :: nVar
    real,    intent(inout) :: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    integer, intent(in), optional :: nCritExt ! num of external criteria
    real, intent(inout), optional :: CritExt_IB(:,:)
    real, intent(in),    optional :: CoarsenCritExt_I(:), RefineCritExt_I(:)

    real :: Crit, Crit_D(MaxDim)
    real :: Numerator, Denominator
    integer:: iBlock, iCrit, i, j, k, iVar
    logical :: DoCoarsen
    integer :: iError
    !------------------------------------------------------------------------
    Crit      = 1.0
    Crit_D    = 0.0
    Numerator   = 0.0
    Denominator = 0.0

    if(nCrit > nVar) &
         call CON_stop("set_amr_criteria :: More criteria then variables")

    if(present(nCritExt)) then
       nAmrCrit = nCrit + nCritExt
       if(nAmrCrit /= nAmrCritOld) then
          nAmrCritOld=nAmrCrit
          if(allocated(GlobalCritMaxAll_I)) deallocate( &
               CritMaxAll_I, &
               GlobalCritMaxAll_I, &
               CoarsenAmrCrit_I, RefineCritAll_I)
          allocate(CritMaxAll_I(nAmrCrit), GlobalCritMaxAll_I(nAmrCrit), &
               CoarsenAmrCrit_I(nAmrCrit), RefineCritAll_I(nAmrCrit))
          CritMaxAll_I = 0.0
          GlobalCritMaxAll_I =0.0
       end if
    else
       nAmrCrit    = nCrit
       nAmrCritOld = nCrit
    end if

    if(nAmrCrit /= nAmrCritOld .or. nBlock /= nBlockOld) then
       if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
       allocate(AmrCrit_IB(nAmrCrit,nBlock))
       nBlockOld=nBlock
    end if

    ! Initialize all AMR criteria
    AmrCrit_IB = 0.0

    ! Calculate error estimates (1..nCrit)
    BLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       DoCoarsen = .true.
       ! Check refinement first
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          ! Only variables indexed in iVarCrit_I will decide the refinement
          do iCrit = 1, nCrit 
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
             ! Crit = sum(Crit_D)/nDim !!! cheaper
             if( Crit >  AmrCrit_IB(iCrit,iBlock)) &
                  AmrCrit_IB(iCrit,iBlock) = Crit

             !find max error Crit
             !if(present(CritExt_I) ) then 
             !   if( Crit > AmrCrit_IB(iCrit,iBlock)) &
             !        AmrCrit_IB(iCrit,iBlock) = Crit
             !else
             !   ! If one cell in the block needs refinement the block will
             !   ! be refined, But only if all cells in the block want to be
             !   ! coarsen the block will be flagged for coarsening
             !
             !   if( Crit > RefineCrit_I(iCrit)) then
             !      iStatusNew_A(iNode_B(iBlock)) = Refine_
             !      DoCoarsen = .false.
             !      CYCLE BLOCK
             !   else if(Crit > CoarsenCrit_I(iCrit)) then
             !      DoCoarsen = .false.
             !   end if
             !end if

          end do !end nVar
       end do; end do; end do

       !if( DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_

    end do BLOCK

    ! Fill in thresholds
    do iCrit=1,nCrit
       RefineCritAll_I(iCrit)  = RefineCrit_I(iCrit)
       CoarsenAmrCrit_I(iCrit) = CoarsenCrit_I(iCrit)
    end do

    ! Add external criteria and thresholds
    if(present(nCritExt)) then

       do iBlock = 1, nBlock
          do iCrit = 1, nCritExt
             AmrCrit_IB(nCrit+iCrit,iBlock) = CritExt_IB(iCrit,iBlock)
          end do
       end do

       do iCrit=1,nCritExt
          RefineCritAll_I(nCrit+iCrit)  = RefineCritExt_I(iCrit)
          CoarsenAmrCrit_I(nCrit+iCrit) = CoarsenCritExt_I(iCrit)
       end do

    end if

    ! Find the max criteria values on all processors
    do iCrit = 1, nAmrCrit
       CritMaxAll_I(iCrit) = maxval(AmrCrit_IB(iCrit,:))
    end do
    if(nProc > 1)then
       call MPI_allreduce(CritMaxAll_I, GlobalCritMaxAll_I,nAmrCrit, MPI_REAL,&
            MPI_MAX,iComm, iError)
    else
       GlobalCritMaxAll_I = CritMaxAll_I
    end if

    ! Normalize the external criteria only
    do iCrit = nCrit+1, nAmrCrit
       if(GlobalCritMaxAll_I(iCrit) /= 0.0) &
            AmrCrit_IB(iCrit,:) = AmrCrit_IB(iCrit,:)/GlobalCritMaxAll_I(iCrit)
    end do

    ! Set refinement and coarsening flags in iStatuNew_A
    BLOCK2:do iBlock = 1, nBlock

       if(Unused_B(iBlock)) CYCLE
       DoCoarsen = .true.
       do iCrit = 1, nAmrCrit

          ! Decide refinement based on normalized error factor for each
          ! chosen variable
          ! If one cell in the block needs refinement the block will
          ! be refined, But only if all cells in the block want to be
          ! coarsen the block will be flagged for coarsening
          if( AmrCrit_IB(iCrit,iBlock) > RefineCritAll_I(iCrit)) then
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             DoCoarsen = .false.
             CYCLE BLOCK2
          else if(AmrCrit_IB(iCrit,iBlock)  > CoarsenAmrCrit_I(iCrit)) then
             DoCoarsen = .false.
          end if
       end do

       if(DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_

    end do BLOCK2

  end subroutine set_amr_criteria

  !============================================================================
  subroutine  read_amr_criteria_param
    use ModReadParam, ONLY: read_var
    integer :: iCrit
    !-------------------------------------------------------------------------
    nCrit          = 0
    nAmrCrit       = 0
    nAmrCritOld    = 0
    cAmrWavefilter = 1.0e-2
    nBlockOld      = 0

    call read_var('AmrWavefilter',cAmrWavefilter)  
    call read_var('nCrit', nCrit)
    nAmrCrit = nCrit
    nAmrCritOld = nCrit
    if(allocated(CoarsenCrit_I)) then
       deallocate(CoarsenAmrCrit_I, RefineCritAll_I)
       deallocate(CoarsenCrit_I, RefineCrit_I, iVarCrit_I,GlobalCritMaxAll_I)
    end if
    allocate(CoarsenCrit_I(nCrit), RefineCrit_I(nCrit),iVarCrit_I(nCrit))
    allocate(CoarsenAmrCrit_I(nCrit), RefineCritAll_I(nCrit))
    allocate(GlobalCritMaxAll_I(nCrit))
    allocate(CritMaxAll_I(nCrit))
    do iCrit = 1, nCrit
       call read_var('iVar',iVarCrit_I(iCrit))
       call read_var('CoarsenCrit',CoarsenCrit_I(iCrit))
       call read_var('RefineCrit',RefineCrit_I(iCrit))
    end do

  end subroutine read_amr_criteria_param
  !============================================================================
  subroutine clean_amr_criteria

    if(.not.allocated(CoarsenCrit_I)) RETURN
    deallocate(CoarsenCrit_I, RefineCrit_I, iVarCrit_I)
    deallocate(CoarsenAmrCrit_I, RefineCritAll_I)
    deallocate(GlobalCritMaxAll_I)
    deallocate(CritMaxAll_I)
    nCrit    = 0
    nAmrCrit = 0
    nAmrCritOld = 0

  end subroutine clean_amr_criteria

end module BATL_amr_criteria
