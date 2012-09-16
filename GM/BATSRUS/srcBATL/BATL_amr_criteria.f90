!^CFG COPYRIGHT UM

! All amr criteria has the following form:
! * Refinemnet criteia, RefineCrit..._I
! * Coursening criteria,  CoarsenCrit...._I 
! * Geometric arae where it will be used, Area_I 
! * Quontity to compere to,  AmrCrit_IB ( grid resolution dx of chenage in Stat variable)
! * Maximum refinment level/resolution, MaxLevelCrit_I

! Overview of varables used for AMR, array dim is a indication of how they are used ONLY
! nAmrCritUsed                       : Total number of criterias used
! nAmrCrit                           : Total number of unique criterias used
! nExtCrit                           : Number of criteria coming form outside BATL, values in AmrCirt_IB
! nIntCirt                           : Number of criteria handeld inside BATL, values in AmrCirt_IB
! nPhysCritUsed                      : Number of criterias given by #AMRCRITERIA...... 
! nCritDxLevel                       : Number of geometic criterias containd in nPhysCritUsed
! nCritGeoBackword                   : Number of geometic criterias from #GRIDLEVEL/RESOLUTION, backward compatiblei
! nGeoCrit                           : =2, dx(1) and level(2)
! RefineCritAll_I[nAmrCritUsed]      : Limit for criteria to do refinment
! CoarsenCritAll_I[nAmrCritUsed]     : Limit for criteria for when to coursen
! ResolutionLimit_I[nAmrCritUsed]    : A criteria will not be applaid for block with a better resolution then indicated
! iResolutionLimit_I[nAmrCritUsed]   : Index if ResolutionLimit_I are applied to dx or level resolution criteria
! iMapToUiqCrit_I[nAmrCritUsed]      : Map form used criterias to unique criterias
! AreaAll_I[nmaxarea]                : List of all ares used by the amr criterias
! AmrCirt_IB[nAmrCrit,nBlock]        : Store the criteria values we compare to
! iAreaIdx_II[nareamax,nAmrCritUsed] : Give the index in AreaAll_I for each criteras areas
! nAreaOPerCrit_I[nAmrCritUsed]      : Number of areas appied for a criteria



module BATL_amr_criteria

  use BATL_mpi,  ONLY: iProc

  use BATL_tree, ONLY: nDesiredRefine,nNodeRefine,nDesiredCoarsen, &
       nNodeCoarsen, DoStrictAmr, iRank_A, Rank_A, nNodeSort

  use BATL_amr_geometry, ONLY: nGeoCrit,clean_amr_geometry,&
       isNewGeoParam, UseCrit_IB, UseCrit_B,AreaType,&
       apply_amr_geometry,nAmrCrit,nAmrCritUsed,AmrCrit_IB,nGeoCrit,&
       nCritGeoUsed,IsBatsrusAmr,iResolutionLimit_I,&
       CoarsenCritAll_I,RefineCritAll_I,ResolutionLimit_I,iVarCritAll_I,&
       AreaGeo_I,iAreaIdx_II,nAreaPerCritAll_I,&
       MaxArea,lNameArea,nCritGeoBackword,nCritDxLevel
  

  implicit none

  SAVE

  private ! except

  public set_amr_criteria
  public clean_amr_criteria
  public read_amr_criteria
  public test_amr_criteria
  public calc_error_amr_criteria
  public set_amr_geometry
  public init_amr_criteria
  public masked_amr_criteria

  ! DoSortAmrCrit : true/false for generating the priority list iRank_A
  ! DoSoftAmrCrit : standard behavior (false) is that criteria has to be 
  ! refined/coarsened or we will stop the program. If true we will 
  ! refine/coarsen them as much as we can, this is the default behavior for
  ! amr by percentage.
  logical, public :: &
       DoSortAmrCrit = .false., &
       DoSoftAmrCrit = .false., &
       DoAutoAmr     = .false.

  ! Need to init if PARAM.in was read (again)
  logical :: isNewPhysParam = .false.


  ! Used for combination of sorted and not sorted criteria
  logical :: DidUnsortedAMR = .false.

  ! Try to make geometric dependence for refinement/coarsening
  logical,public :: DoCritAmr     = .false.

  ! Percentage with want refined and coarsen
  real ::  PercentRefine=0.0, PercentCoarsen=0.0

  ! Threshold limits for refine or unrefined the grid (length nCrit) 
  real, allocatable, dimension(:)    :: &
       CoarsenCritPhys_I, RefineCritPhys_I

  
  ! Give max level (-) or resolution dependent on which index 
  ! idxMaxGeoCritPhys_I points to
  real, allocatable   :: MaxLevelCritPhys_I(:)
  integer,allocatable :: idxMaxGeoCritPhys_I(:)

  ! iMapToUiqCrit_I map indexes from all indexes (nAmrCritUsed) to 
  ! number of uniqe indexes (nAmrCrit) used by AmrCirt_IB
  integer, allocatable :: iMapToUiqCrit_I(:)

  ! Map the internal error criteria index to the state varable index
  integer, allocatable :: iMapToStateVar_I(:) 

  ! Parameters used by calc_error_amr_criteria to estimate the errors
  real :: cAmrWavefilter = 1.0e-2
  real, parameter :: cEpsilon = 1.0d-8 ! avoid zero in denominator

  integer :: nExtCrit = 0 ! Number of External criteria + sorting criteia from BATSRUS
  integer :: nIntCrit = 0 ! Number of internal criteria, 2nd order err estimate
  integer :: nPhysCritUsed = 0 ! total Number criterias read by read_amr_criteria
  integer :: iStartCrit =0, iEndCrit=0 !start of ccriteria segement we want to work on
  
  ! Storing names of areas for each criteria given by #AMRCRITERIA.....
  integer, public, allocatable :: nAreaPerCritPhys_I(:)

  ! How large the relative discrepancy in the Criteria for 
  ! symmetry points can be and still be refined in the same 
  ! way and the same for geometric variation
  real, parameter :: DeltaSymCrit = 1.0e-2
  real, parameter :: DeltaSymGeo  = 1.0e-6

  ! The max differens in the criteia to say thay have the same value,
  ! taking care of numerical fluctations
  real :: DeltaCritera = 1.0e-8

  ! Number of characters on TypeCriteria
  integer, parameter :: nChar=200
  ! store area type for each criteria [nArea,nCrit]
  character(len=nChar), allocatable :: AreaNamesPhys_II(:,:)

  ! Converting index form running block and proc number to node numbers
  integer, allocatable :: iNode_I(:)

contains
  !============================================================================

  subroutine init_amr_criteria(user_amr_geometry)

    use BATL_size, ONLY: MaxBlock,nBlock,MinI, MaxI,MinJ, MaxJ, MinK, MaxK, nDim
    use BATL_tree, ONLY: Unused_B
    use BATL_grid, ONLY: CoordMin_D, CoordMax_D
    use BATL_amr_geometry, ONLY: init_amr_geometry
    use BATL_mpi,  ONLY: iProc
    use ModNumConst,  ONLY: cUnit_DD

    interface
       subroutine user_amr_geometry(iBlock, iArea, DoRefine)
         integer, intent(in) :: iBlock, iArea
         logical,intent(out) :: DoRefine
       end subroutine user_amr_geometry
    end interface
    optional :: user_amr_geometry

    integer :: iBlock,iCrit,nCrit,iGeo,iArea,AreaSign
    logical :: DoTestMe = .false.
    character(len=lNameArea) :: regionname

    character(len=18), parameter :: NameSub = 'init_amr_criteria'
    !-------------------------------------------------------------------------

    if(.not.(isNewGeoParam .or. isNewPhysParam)) RETURN
    
    AreaGeo_I(0)%NameRegion = "ALL"
    AreaGeo_I(0)%Name ="all"
    AreaGeo_I(0)%Resolution = 0.0
    AreaGeo_I(0)%Center_D   = 0.0
    AreaGeo_I(0)%Size_D     = 0.0
    AreaGeo_I(0)%Radius1    = 0.0
    AreaGeo_I(0)%DoRotate   = .false.
    AreaGeo_I(0)%Rotate_DD  = cUnit_DD

    nAmrCrit =  nIntCrit + nExtCrit + nGeoCrit

    if(allocated(AmrCrit_IB)) deallocate(AmrCrit_IB)
    allocate(AmrCrit_IB(nAmrCrit,MaxBlock))
    AmrCrit_IB = -1.0

    ! Merging Data form Geometical and Physical refinment parameters
    nCrit = nAmrCritUsed + nCritGeoBackword

    if(allocated(RefineCritAll_I))then
         deallocate(RefineCritAll_I,&
                    CoarsenCritAll_I,&
                    iVarCritAll_I,&
                    ResolutionLimit_I,&
                    iResolutionLimit_I,&
                    iAreaIdx_II,&
                    nAreaPerCritAll_I)
    endif
    allocate(RefineCritAll_I(nCrit),&
             CoarsenCritAll_I(nCrit),&
             iVarCritAll_I(nCrit),&
             ResolutionLimit_I(nCrit),&
             iResolutionLimit_I(nCrit),&
             iAreaIdx_II(MaxArea,nCrit),&
             nAreaPerCritAll_I(nCrit))

    iVarCritAll_I = 0

    if(nAmrCritUsed > 0 ) then
       !Copy physics based criterias, from read_amr_criteia
       RefineCritAll_I(1:nAmrCritUsed)   = RefineCritPhys_I(1:nAmrCritUsed)
       CoarsenCritAll_I(1:nAmrCritUsed)  = CoarsenCritPhys_I(1:nAmrCritUsed)
       ResolutionLimit_I(1:nAmrCritUsed) = MaxLevelCritPhys_I(1:nAmrCritUsed)
       iVarCritAll_I(1:nAmrCritUsed)     = iMapToUiqCrit_I(1:nAmrCritUsed)
       nAreaPerCritAll_I(1:nAmrCritUsed) = nAreaPerCritPhys_I(1:nAmrCritUsed)

       ! fixing indexing as "physics" criteria comes first
       do iCrit=1,nAmrCritUsed
          iResolutionLimit_I(iCrit)   = max(0,maxval(iVarCritAll_I(1:nAmrCritUsed-nCritDxLevel))) +&
               idxMaxGeoCritPhys_I(iCrit)
       end do
    end if

    ! Copy over Geometry based criterias
    if(nCritGeoUsed > 0) call init_amr_geometry
    if(nAmrCritUsed > 0 ) then
       do iCrit=1,nAmrCritUsed
          BLOCKAREA: do iArea=1,nAreaPerCritAll_I(iCrit)
             regionname = trim(adjustl(AreaNamesPhys_II(iArea,iCrit)))

             if(regionname(1:1) == "-") then
                areasign = -1
                regionname =regionname(2:LEN_TRIM(regionname))
             else if (regionname(1:1) == "+") then
                areasign = +1
                regionname = regionname(2:LEN_TRIM(regionname))
             else
                areasign = +1
                regionname = regionname
             end if

             do iGeo=0,nCritGeoUsed
                if(trim(adjustl(regionname)) == trim(adjustl(AreaGeo_I(iGeo)%NameRegion))) then
                   iAreaIdx_II(iArea,iCrit) = areasign*iGeo
                   CYCLE BLOCKAREA
                end if
             end do

             if(iGeo > nCritGeoUsed ) &
                  call CON_stop(NameSub //' Can not find area name :'//trim(adjustl(regionname)))
          end do BLOCKAREA
       end do
    end if
    ! Total number of refine/coursening limmits
    nAmrCritUsed = nCrit

    if(iProc == 0) then
    end if

    ! Allocate Criteria mask
    if(.not. allocated(UseCrit_B)) allocate(UseCrit_B(MaxBlock))
    if(allocated(UseCrit_IB)) deallocate(UseCrit_IB)
    allocate(UseCrit_IB(nAmrCritUsed,MaxBlock))

    UseCrit_IB = .false.
    UseCrit_B  = .false.

    do iBlock=1,nBlock
       if(Unused_B(iBlock)) CYCLE
       call set_amr_geometry(iBlock,user_amr_geometry=user_amr_geometry)
    end do

    if(DoTestMe) then
       if(iProc == 0) then
          write(*,"(A17,100(F10.3))") "RefineCritAll_I   : ", RefineCritAll_I
          write(*,"(A17,100(F10.3))") "CoarsenCritAll_I  : ", CoarsenCritAll_I
          write(*,"(A17,100(F10.3))") "ResolutionLimit_I : ", ResolutionLimit_I
          write(*,"(A17,100(I10))") "iResolutionLimit_I     : ",iResolutionLimit_I
          write(*,"(A17,100(I10))") "iVarCritAll_I       : ", iVarCritAll_I
          write(*,"(A17,100(A10))") "AreaGeo_I%Name      : ", AreaGeo_I%Name
          write(*,"(A17,100(I10))") "iVarCritAll_I       : ", iVarCritAll_I
          write(*,"(A17,100(I10))") "AmrCrit_IB          : ", shape(AmrCrit_IB) 
          write(*,*) ""
          do iCrit=1,nAmrCritUsed
             write(*,*) "iAreaIdx_II( :,iCrit) ",iAreaIdx_II( 1:nAreaPerCritAll_I(iCrit),iCrit)
          end do
          write(*,*) ""
          write(*,*) " 1 : nAmrCritUsed = ", 1," : ",nAmrCritUsed
          write(*,*) " Geo range        = ", nPhysCritUsed +1," : ",  nAmrCritUsed
          write(*,*) " Phys range       = ",1," : ", nPhysCritUsed
          write(*,*) " Error Phys       = ",1," : ", nIntCrit
          write(*,*) " Ext calc range   = ",nIntCrit+1," : ",nIntCrit+nExtCrit
          write(*,*) " shape AmrCrit_IB = ",shape(AmrCrit_IB)
          write(*,*) ""
          write(*,*) ""
       end if
    end if

    ! All nesesary things are updated
    isNewPhysParam = .false.
    isNewGeoParam  = .false.

    ! cleaning arrays for read_amr_criteria
    if(nPhysCritUsed > 0 ) &
         deallocate(RefineCritPhys_I,&
         CoarsenCritPhys_I,&
         MaxLevelCritPhys_I, &
         iMapToUiqCrit_I, &
         idxMaxGeoCritPhys_I,&
         nAreaPerCritPhys_I,&
         AreaNamesPhys_II)    


  end subroutine init_amr_criteria
  !============================================================================

  subroutine set_amr_criteria(nVar, State_VGB, nInCritExtUsed, &
       CritExt_IB, Used_GB,TypeAmrIn,user_amr_geometry)

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
         MaxBlock,nBlock
    use BATL_tree, ONLY: Unused_B!,iStatusNew_A,Refine_

    interface
       subroutine user_amr_geometry(iBlock, iArea, DoRefine)
         integer, intent(in) :: iBlock, iArea
         logical,intent(out) :: DoRefine
       end subroutine user_amr_geometry
    end interface
    optional  :: user_amr_geometry

    integer,  intent(in)  :: nVar
    real,    intent(in),optional  :: & ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    integer, intent(in), optional :: nInCritExtUsed
    real, intent(in), optional :: CritExt_IB(nExtCrit, nBlock)
    logical, intent(in), optional:: &
         Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    character(3), intent(in), optional:: TypeAmrIn



      character(len=*), parameter :: NameSub = 'set_amr_criteria'
      character(3) :: TypeAmr
      integer :: iCrit, iBlock
      !-----------------------------------------------------------------------

      !nExtCritUsed = 0
      if(present(nInCritExtUsed)) then
         if(nInCritExtUsed /= nExtCrit) then
            write(*,*) nInCritExtUsed," /= ",nExtCrit
            write(*,*) ""
            call CON_stop("nExtCrit can only change by changes in PARAM.in")  
         end if
      end if

      TypeAmr = 'all'
      if(present(TypeAmrIn)) TypeAmr = TypeAmrIn 

      ! Compatible with old BATSRUS amr behavior
      if(IsBatsrusAmr .and.  TypeAmr=='all') then
         if(nPhysCritUsed > 0) then
            TypeAmr = 'phy'
         else
            TypeAmr = 'geo' 
         end if
      end if

      ! add external criteria into the list of all criteria
      if(present(CritExt_IB)) then
         do iBlock = 1, nBlock
            if(Unused_B(iBlock)) CYCLE
            do iCrit = 1, nExtCrit
               AmrCrit_IB(nIntCrit+iCrit,iBlock) = CritExt_IB(iCrit,iBlock)
            end do
         end do
      end if

      !-------- set up index Ranges and find blocks to refine/coursen -------
      select case(TypeAmr)
      case('all')
         iStartCrit = 1
         iEndCrit   = nAmrCritUsed
      case('geo')
         iStartCrit = nPhysCritUsed-nCritDxLevel +1
         iEndCrit   = nAmrCritUsed
         call apply_unsorted_criteria
         RETURN
      case('phy')
         iStartCrit = 1
         iEndCrit  = nPhysCritUsed-nCritDxLevel
      case default
         call CON_stop(NameSub // &
              ' ERROR: Unknown TypeAmr = '//TypeAmr)
      end select

      ! Estimation of the numerical error
      if(nIntCrit >0) &
           call calc_error_amr_criteria(nVar, State_VGB, Used_GB=Used_GB)

      if(DoSortAmrCrit .or. .not.DoStrictAmr) then
         ! we make an AMR priority list
         ! Only sort Physics based AMR
         ! Geomtry criteria is all or none 
         if(TypeAmr /= 'phy') then            !!! ALL !!!!
            iStartCrit = nPhysCritUsed-nCritDxLevel +1
            iEndCrit   = nAmrCritUsed
            call apply_unsorted_criteria
            DidUnsortedAMR = .true.
         end if
         ! Physical criterias
         iStartCrit = 1
         iEndCrit  = nPhysCritUsed
         call sort_amr_criteria
      else
         ! refine only based on criteria
         call apply_unsorted_criteria
      end if

      DidUnsortedAMR = .false.

    end subroutine set_amr_criteria
  !===========================================================================
  subroutine sort_amr_criteria

    ! the routine will find the candidate for refinement or coursing based
    ! upon the percentage of refinement and the percentage wanted refined and
    ! and coarsened. It will ONLY modify the criteria the values in the 
    ! criteria list itself. The refinement criteria will overrule the 
    ! percentage refinement.

    use BATL_mpi, ONLY: iComm, nProc
    use ModMpi
    use BATL_size, ONLY: nBlock
    use ModSort, ONLY: sort_quick
    use BATL_tree, ONLY: Unused_BP,iStatusNew_A, Refine_, Coarsen_, &
         Unset_, iNode_B, nNode, nNodeUsed, diffRange,Unused_B

    ! Array containing all criteria for all blocks gathered
    real, allocatable :: AllCrit_II(:,:)

    ! CritSort_II contains all criteria for all blocks
    ! iIdxSort_II gives the new indexes after sorting
    ! and reused for other functions.
    integer, allocatable :: iIdxSort_II(:,:)
    real, allocatable :: CritSort_II(:,:)
    integer, allocatable :: iRank_I(:)
    real, allocatable :: Rank_I(:),nRank_I(:)
    logical, allocatable ::  AllUSeCrit_II(:,:)
    ! AllRank_II(Nodes given by iNode_I,2) [rank,N ranks for node]
    real, allocatable :: AllRank_II(:,:)
    ! tranlsate sorting continous index to iNode_I index
    integer, allocatable :: iSortToNode_I(:)
    ! will have value 0: unset, 1:refine and -1:coursen to count
    ! the number of nodes to refine and/or coursen by criteria
    integer, allocatable ::nAmrTag(:)
    integer :: iCrit, iBlock, iProces, k
    integer :: iCritSort, iSort,iVarCrit,nNodeSortMax,iNodeSort
    integer :: iError
    real    :: Diff, Crit

    integer :: nTotBlocks, nBlockMax, iTotalCrit, iHelpCrit
    integer, allocatable :: nBlock_P(:), nReciveCont_P(:), nRecivDisp_P(:)
    integer, allocatable :: nUseReciveCont_P(:), nUseRecivDisp_P(:)

    ! indexes for AmrCrit_IB as we can have multiple crit using the same
    ! AmrCrit_IB criteria
    integer :: nVarStart, nVarEnd
    !-----------------------------------------------------------------------

    ! need extra array space to store temperal data
    iTotalCrit = 2
    iHelpCrit  = 3

    nVarStart = iVarCritAll_I(iStartCrit)
    nVarEnd   = iVarCritAll_I(iEndCrit)

    ! COMMENTS FOR FUTURE : we can use the Unused_B to only pick and send
    ! the data that is actively used, reducing communication


    !------------ Get information on all the blocks on all process -----------

    ! collect nBlock from each process 
    allocate(nBlock_P(nProc),nReciveCont_P(nProc),nRecivDisp_P(nProc))
    allocate(nUseReciveCont_P(nProc),nUseRecivDisp_P(nProc))
    call MPI_Allgather(nBlock, 1,MPI_INTEGER, nBlock_P, 1, &
         MPI_INTEGER, iComm, iError)
    nBlockMax = maxval(nBlock_P)
    nTotBlocks = sum(nBlock_P)

    !------------ Collect all criteria for sorting and ranking ------------

    ! collect all criteria for the total grid
    allocate(AllCrit_II(nAmrCrit,nTotBlocks))
    allocate(AllUSeCrit_II(nAmrCritUsed,nTotBlocks))
    AllCrit_II = -77

    ! store the Node indexing iNode_B for all processes
    if(allocated(iNode_I)) deallocate(iNode_I)
    allocate(iNode_I(nBlockMax*nProc))
    allocate(AllRank_II(nBlockMax*nProc,2))
    allocate(iSortToNode_I(nBlockMax*nProc))
    allocate(nAmrTag(nBlockMax*nProc))
    
    iNode_I = 0
    iSortToNode_I = -1
    AllRank_II =0
    nAmrTag = 0

    ! Set up the displacement and size for collecting the data for 
    ! AllCrit_II
    nReciveCont_P = nBlock_P*nAmrCrit
    nRecivDisp_P(1) = 0
    do iProces=2,nProc 
       nRecivDisp_P(iProces) = &
            nRecivDisp_P(iProces-1)+nBlock_P(iProces-1)*nAmrCrit
    end do

    ! Set up the displacement and size for collecting the data for 
    ! AllUseCrit_II
    nUseReciveCont_P = nBlock_P*nAmrCritUsed
    nUseRecivDisp_P(1) = 0
    do iProces=2,nProc 
       nUseRecivDisp_P(iProces) = &
            nUseRecivDisp_P(iProces-1)+nBlock_P(iProces-1)*nAmrCritUsed
    end do

    ! Get all criteria maskes, so we can exlude them
    call MPI_allgatherv(UseCrit_IB, nAmrCritUsed*nBlock, MPI_LOGICAL, &
         AllUseCrit_II, nUseReciveCont_P, nUseRecivDisp_P, &
         MPI_LOGICAL, iComm, iError)

       do iBlock=1,nBlock
          if(Unused_B(iBlock)) CYCLE
          !print *,iProc," inode (iBlock) :: ",iBlock,  iNode_B(iBlock)
          if(iStatusNew_A(iNode_B(iBlock))== Refine_) then
             AmrCrit_IB(nVarStart:nVarEnd,iBlock) = -777.0
             !print *,"WHAT Refine_", nVarStart,nVarEnd,iBlock,iNode_B(iBlock)
          end if
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
         iRank_A(nTotBlocks), Rank_A(nTotBlocks),nRank_I(nTotBlocks))

    iRank_A      = 0
    Rank_A       = 0.0
    iRank_I      = 0
    Rank_I       = 0.0
    iIdxSort_II  = 0
    CritSort_II  = 0.0
    nRank_I      = -1.0
    iRank_I = nNode+1 ! larger than any possible rank
    iIdxSort_II(:,iHelpCrit) = 0 ! store info about if a block is already 
    !                               marked for refinement(+1) or coarsening(-1)
    nNodeRefine  = 0
    nNodeCoarsen = 0
    nNodeSortMax = 0

    ! Working on 1 criteria at the time
    do iCrit=iStartCrit,iEndCrit


       iVarCrit = iVarCritAll_I(iCrit)
       ! copy criteria data into the sorting arrays
       iSort = 0
       do iProces = 1, nProc
          do iBlock = 1, nBlock_P(iProces)
             if(Unused_BP(iBlock,iProces-1)) CYCLE
             if(.not.AllUseCrit_II(iCrit,iBlock+nRecivDisp_P(iProces))) then
                CYCLE
             end if
             if(AllCrit_II(iVarCrit,iBlock+nRecivDisp_P(iProces)) == -777.0) then
                CYCLE
             end if
             iSort = iSort +1
             CritSort_II(iSort,1) = &
                  AllCrit_II(iVarCrit,iBlock+nRecivDisp_P(iProces))
             iIdxSort_II(iSort, 1) = iSort
             iSortToNode_I(iSort) = iBlock+nRecivDisp_P(iProces)
             ! nRank_I are reacalulated each time, but at the moment we have noe
             ! constant indexing so to get the corect match we to it each time
             nRank_I(iSortToNode_I(iSort)) = &
                  1.0/real(nTotBlocks*count(AllUseCrit_II(iStartCrit:iEndCrit,iBlock+nRecivDisp_P(iProces))))
          end do
       end do

       if(isort == 0) CYCLE

       nNodeSort = iSort
       nNodeSortMax = max(nNodeSort,nNodeSortMax)

        ! Sort each criteria. Pass in first index of array arguments
       call sort_quick(nNodeSort, &
            CritSort_II(1:nNodeSort,1), iIdxSort_II(1:nNodeSort,1))

       ! Store the sorting index in secound colum (iTotalCrit)
       iIdxSort_II(:,iTotalCrit) = iIdxSort_II(:,1)
       k = 1
       iSort = iIdxSort_II(1,iTotalCrit)
       !Value of first element (min)
       iIdxSort_II(iSort,1) = k
       Crit = CritSort_II(iSort,1)


       ! group together criteria which has a diffrence of 
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
          iNodeSort = iSortToNode_I(iCritSort)
          ! Use the minimum position (=rank) of the node in sorted criteria 
          ! array to set up the new ranking order for each node when doing amr
          iRank_I(iNodeSort) = min(iRank_I(iNodeSort),iIdxSort_II(iCritSort,1))

          ! if some of the nodes have the same minimal rank we will use the 
          ! average rank of all criteria as a secondary priority.
          Rank_I(iNodeSort) = Rank_I(iNodeSort)+nRank_I(iNodeSort)*iIdxSort_II(iCritSort,1)
       end do


       ! To make sure that blocks that are marked for refinement/coarsening 
       ! by thresholds will come at the top/bottom of the ranking list 
       ! we will shift the values with the total number of blocks used.
       ! COMMENTS FOR FUTURE : the loop may be split for one that goes 
       ! form the top and one that goes form bottom so it do not need to go 
       ! though the whole list of elements

       ! Only if Criterias is in use
       if(RefineCritAll_I(iCrit) > -0.5 .and. CoarsenCritAll_I(iCrit) > -0.5) then
          do iSort = nNodeSort,1,-1
             iCritSort = iIdxSort_II(iSort,iTotalCrit)
             iNodeSort = iSortToNode_I(iCritSort)
             ! Block has already ben marked
             !if(iIdxSort_II(iCritSort,iHelpCrit) /= 0 ) CYCLE

             if(CritSort_II(iCritSort,1) > RefineCritAll_I(iCrit) &
                  .and. iIdxSort_II(iCritSort,iHelpCrit) /= 1) then
                ! Satisfies refinement threshold and not yet marked for refinement

                ! tag for refinment
                nAmrTag(iNodeSort) = 1 

                ! Shift up the rank above sorted values                
                Rank_I(iNodeSort) = Rank_I(iNodeSort) + nTotBlocks+1


!                ! If node was originally marked for coarsening, 
!                ! now there is one fewer node to be coarsened
!                if(iIdxSort_II(iCritSort,iHelpCrit) == -1) &
!                     nNodeCoarsen = nNodeCoarsen - 1

                ! Noe node is marked fore refinement
                !nNodeRefine = nNodeRefine +1
                iIdxSort_II(iCritSort,iHelpCrit) = 1

             else if(CritSort_II(iCritSort,1) < CoarsenCritAll_I(iCrit) &
                  .and. iIdxSort_II(iCritSort,iHelpCrit) == 0) then
                ! Satisfies coarsening threshold and not yet marked at all !

                ! Shift down the rank below sorted values
                Rank_I(iNodeSort) = -Rank_I(iNodeSort)-nTotBlocks-1

                ! Now node is marked for coarsening
                !nNodeCoarsen = nNodeCoarsen +1
                nAmrTag(iNodeSort) = -1 
                iIdxSort_II(iCritSort,iHelpCrit) = -1 
             end if
          end do
       end if

      nNodeCoarsen = count(nAmrTag == -1)
      nNodeRefine  = count(nAmrTag == 1)

    end do

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
          if(.not.any(AllUseCrit_II(iStartCrit:iEndCrit,k))) CYCLE
          if(all(AllCrit_II(nVarStart:nVarEnd,k) == -777.0)) CYCLE
          iSort = iSort+1
          Rank_I(iSort) = iRank_I(k) + Rank_I(k)
          iIdxSort_II(iSort,1) = iNode_I(k)
          CritSort_II(iSort,1) = &
                AllCrit_II(iResolutionLimit_I(iStartCrit),iBlock+nRecivDisp_P(iProces))
       end do
    end do
    nNodeSort = iSort

    ! Get the number of elements to refine by percentage
    ! taking into account that each refinements generate nChild-1 
    ! new blocks
    nDesiredRefine  = floor(PercentRefine*nNodeUsed/(100.0))
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
       if( maxval(ResolutionLimit_I(iStartCrit:iEndCrit)) >= CritSort_II(iIdxSort_II(iSort,iTotalCrit),1)) CYCLE
       iStatusNew_A(iRank_A(iSort)) = Refine_
    end do
   
    if(DidUnsortedAMR) then
       ! Do not coarsen blocks not in the sorted list
       do iSort = max(nDesiredCoarsen, nNodeCoarsen),nNodeSort
         if(iStatusNew_A(iRank_A(iSort)) ==  Coarsen_)&
             iStatusNew_A(iRank_A(iSort)) = Unset_
       end do
   else
       do iSort = 1, max(nDesiredCoarsen, nNodeCoarsen)
         if(iStatusNew_A(iRank_A(iSort)) /= Refine_)&
             iStatusNew_A(iRank_A(iSort)) =  Coarsen_
       end do
   end if

    deallocate(iSortToNode_I)
    deallocate(AllRank_II)
    deallocate(AllUseCrit_II)
    deallocate(AllCrit_II)
    deallocate(CritSort_II, iIdxSort_II)
    deallocate(nBlock_P, nReciveCont_P, nRecivDisp_P)
    deallocate(iRank_I, Rank_I)
    deallocate(nRank_I)
    deallocate(nAmrTag)


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

  subroutine calc_error_amr_criteria(nVar, State_VGB,Used_GB )

    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
	 nI, nJ, nK, MaxDim, nDim, MaxBlock, nBlock
    use BATL_tree, ONLY: Unused_B

    integer,  intent(in)   :: nVar
    real,    intent(in) :: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    logical, intent(in), optional:: &
         Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)  ! used cells

    real :: Crit, Crit_D(MaxDim)
    real, parameter :: InvNdim = 1.0/nDim
    real :: Numerator, Denominator
    integer:: iBlock, iCrit, i, j, k, iVar
    ! number of blocks set out for refining and coarsning
    !------------------------------------------------------------------------
    Crit        = 1.0
    Crit_D      = 0.0
    Numerator   = 0.0
    Denominator = 0.0

    if(nIntCrit > nVar) &
         call CON_stop("calc_error_amr_criteria :: More criteria then variables")
    
    ! Calculate error estimates (1..nIntCrit)
    BLOCK: do iBlock = 1, nBlock

       ! Check refinement first
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          ! Only variables indexed in iVarCritAll_I will decide the refinement

          if(present(Used_GB)) then
             if(nDim == 1) then
                if(.not.all(Used_GB((i-2):(i+2),j, k,iBlock))) CYCLE
             end if
             if(nDim == 2) then
                if(.not.all(Used_GB((i-2):(i+2),(j-2):(j+2),k,iBlock))) CYCLE
             end if
             if(nDim == 3) then
                if(.not.all(Used_GB((i-2):(i+2),(j-2):(j+2),(k-2):(k+2),iBlock))) CYCLE
             end if
          end if

          do iCrit = 1, nIntCrit 
             if(.not.UseCrit_IB(iCrit,iBlock)) CYCLE
             
             
             iVar = iMapToStateVar_I(iCrit)

             Crit_D =0.0

             Numerator = abs( &
                  State_VGB(iVar,i-2,j,k,iBlock)   - &
                  2.0*State_VGB(iVar,i,j,k,iBlock) + &
                  State_VGB(iVar,i+2,j,k,iBlock) )
             Denominator =  (&
                  abs(State_VGB(iVar,i+2,j,k,iBlock) - &
                  State_VGB(iVar,i,j,k,iBlock)) + &
                  abs(State_VGB(iVar,i,j,k,iBlock) - &
                  State_VGB(iVar,i-2,j,k,iBlock))) + &
                  cAmrWavefilter * (&
                  abs(State_VGB(iVar,i+2,j,k,iBlock)) + &
                  abs(2.0*State_VGB(iVar,i,j,k,iBlock)) + &
                  abs(State_VGB(iVar,i-2,j,k,iBlock))) 
             Crit_D(1) = (Numerator/max(Denominator,cEpsilon)) 
             Crit = Crit_D(1)**2

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
                Crit = Crit + Crit_D(2)**2
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
                Crit = Crit + Crit_D(3)**2
             end if

             !Crit = (sum(Crit_D(1:nDim)**2))*InvNdim
             Crit = Crit *InvNdim
             if( Crit >  AmrCrit_IB(iCrit,iBlock)) &
                  AmrCrit_IB(iCrit,iBlock) = Crit

          end do !end nVar
       end do; end do; end do

    end do BLOCK

  end subroutine calc_error_amr_criteria

  !============================================================================

  subroutine apply_unsorted_criteria

    use BATL_size, ONLY: nBlock
    use BATL_tree, ONLY: iStatusNew_A, Refine_, Coarsen_, &
         Unused_B, iNode_B, iNode_B

    integer:: iBlock, iCrit, iVarCrit
    ! number of blocks set out for refining and coarsning
    logical :: DoCoarsen
    !----------------------------------------------------------------------

    BLOCK3:do iBlock = 1, nBlock

       if(Unused_B(iBlock)) CYCLE

       if(iStatusNew_A(iNode_B(iBlock)) == Refine_) CYCLE

       DoCoarsen = .true.


       do iCrit = iStartCrit, iEndCrit

          ! Criteria used in this block 
          if(.not.UseCrit_IB(iCrit,iBlock)) CYCLE

          iVarCrit = iVarCritAll_I(iCrit)

          if(AmrCrit_IB(iVarCrit,iBlock) > RefineCritAll_I(iCrit) .and. &
               AmrCrit_IB(iResolutionLimit_I(iCrit),iBlock) > ResolutionLimit_I(iCrit))then
             iStatusNew_A(iNode_B(iBlock)) = Refine_
             CYCLE BLOCK3
          else if(AmrCrit_IB(iVarCrit,iBlock) > CoarsenCritAll_I(iCrit))then
             DoCoarsen = .false.
          end if

       end do

       if(DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_ 

    end do BLOCK3

  end subroutine apply_unsorted_criteria

  !============================================================================

  subroutine read_amr_criteria(NameCommand, nCritInOut, NameCritOut_I,&
       NameStatVarIn_V, nStateVarIn, ReadExtraOut)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: lower_case
    use BATL_tree,    ONLY: MaxTotalBlock, iTree_IA,MaxLevel_,&
         MaxLevel
    use BATL_size,    ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
         MaxBlock, MaxDim, nDim
    use BATL_mpi,     ONLY: iProc
    use ModUtilities, ONLY: split_string

    character(len=*), intent(in) :: NameCommand
    integer :: iCrit,iAmrBox,iCritName,nCrit,iIntCrit,iStatVar

    ! Number and name of criteria to be used by BATSRUS 
    integer,           optional, intent(inout) :: nCritInOut
    character(len=20), optional, intent(out)   :: NameCritOut_I(:)
    logical,           optional, intent(out)   :: ReadExtraOut
    logical :: IsUniqueCritName, UseErrorCrit, ReadExtra, DoTestMe=.false.

    integer,           optional, intent(in):: nStateVarIn
    character(len=*),  optional, dimension(:), intent(in):: NameStatVarIn_V

    ! Max number of componets in TypeCriteia
    integer, parameter :: nMaxComponents=20

    character(len=nChar) :: CritName_I(nMaxComponents)
    integer :: nCritArgs,iCritPhy
    character(len=nChar) :: CritName
    character(len=20) :: NameStatVar
    logical :: IsLevel,IsRes
    logical :: DoAmr
    integer :: DnAmr
    real    :: DtAmr

    character(len=*), parameter:: NameSub='BATL_amr_criteria::read_amr_criteria'
    !-------------------------------------------------------------------------
    ReadExtra = .false.
    isNewPhysParam = .true.

    DoSortAmrCrit = .not. DoStrictAmr
    select case(NameCommand)
    case("#AMRCRITERIA", "#AMRCRITERIALEVEL", "#AMRCRITERIARESOLUTION")

       nCrit       = nCritInOut
       nAmrCrit    = nCrit
       iCritPhy = 0
       ! 'dx' and 'level' will be stord as the last indexes the array
       ! nCritDxLevel will give the number of geometric AMR agruments at the end
       nCritDxLevel  =  0 
       if(nCrit == 0) RETURN 

       ! deallocate,if they are already allocated
       if(allocated(CoarsenCritPhys_I)) &
            deallocate(CoarsenCritPhys_I, RefineCritPhys_I, &
            iMapToUiqCrit_I)
       if(allocated(MaxLevelCritPhys_I)) &
            deallocate(MaxLevelCritPhys_I,idxMaxGeoCritPhys_I)

       ! allocate all arrays
       allocate(CoarsenCritPhys_I(nCrit), RefineCritPhys_I(nCrit), &
            iMapToUiqCrit_I(nCrit),MaxLevelCritPhys_I(nCrit),&
            idxMaxGeoCritPhys_I(nCrit),nAreaPerCritPhys_I(nCrit),&
            AreaNamesPhys_II(nMaxComponents,nCrit))

       
       nAreaPerCritPhys_I = 0
       AreaNamesPhys_II = "" 

       if(present(nStateVarIn) .and. .not. allocated(iMapToStateVar_I))&
            allocate(iMapToStateVar_I(nCrit))

       ! Can only have #AMRCRITERIALEVEL or #AMRCRITERIARESOLUTION 
       ! in one session
       ! we always check for max Level
       idxMaxGeoCritPhys_I = 2
       IsLevel = .false.
       IsRes   = .false.
       if(NameCommand == "#AMRCRITERIALEVEL") then
          IsLevel = .true.
       else if(NameCommand == "#AMRCRITERIARESOLUTION") then
          IsRes = .true.
          idxMaxGeoCritPhys_I = 1
       end if

       ! Turn off old BATSRUS behavior when using 
       ! #AMRCRITERIALEVELg
       IsBatsrusAmr = .not.( IsLevel .or. IsRes)

       nIntCrit = 1
       nCritInOut = 0
       NameCritOut_I = "NULL"
       iMapToStateVar_I = -777
       UseErrorCrit = .false.

       do iCrit = 1, nCrit

          IsUniqueCritName = .true.
          !find index of the criteria from its name
          call read_var('TypeCriteria', CritName, IsLowerCase=.true.)

          call split_string(CritName,nMaxComponents, CritName_I, nCritArgs)

          CritName = trim(adjustl(CritName_I(1)))
          select case(CritName)
          case('error')
             iCritPhy = iCritPhy +1
             ! if(CritName(1:5) == 'error') then
             UseErrorCrit = .true.

             if(.not. present(nStateVarIn) &
                  .and. .not. present(NameStatVarIn_V))&
                  call CON_stop(NameCommand//' ERROR: Need a name table')

             ! find index associated with name => iStatVar
             do iStatVar = 1, nStateVarIn
                NameStatVar = NameStatVarIn_V(iStatVar)
                call lower_case(NameStatVar)
                if(NameStatVar == trim(adjustl(CritName_I(2))) ) EXIT
             end do
             !find if index iStatVar exsit in list before
             do iIntCrit = 1, nIntCrit
                if(iMapToStateVar_I(iIntCrit) == iStatVar) then
                   IsUniqueCritName = .false.
                end if
             end do
             ! if unique it will be added to list
             if(IsUniqueCritName .and. iStatVar <= nStateVarIn) then
                iMapToStateVar_I(nIntCrit) = iStatVar
                iMapToUiqCrit_I(iCritPhy) = -nIntCrit
                nIntCrit = nIntCrit+1
             end if

             call SetCritArea(3,nCritArgs,CritName_I,iCritPhy)
          case('transient')
             iCritPhy = iCritPhy +1
             !elseif(CritName(1:9) == 'transient') then
             if(.not. present(ReadExtraOut)) call CON_stop(NameCommand//&
                  ' ERROR: BATSRUS need flag to read more data')
             ReadExtra = .true.

             ! Find out it the name has bin used before
             do iCritName = 1, nCrit
                if(trim(adjustl(CritName_I(2))) == NameCritOut_I(iCritName)) then
                   iMapToUiqCrit_I(iCritPhy) = iCritName
                   IsUniqueCritName = .false.
                   EXIT
                end if
             end do

             ! Add it to the list if its unique
             if(IsUniqueCritName) then
                do iCritName = 1, nCritInOut+1
                   if(NameCritOut_I(iCritName) == "NULL" )then
                      iMapToUiqCrit_I(iCritPhy) = iCritName
                      NameCritOut_I(iCritName) = trim(adjustl(CritName_I(2)))
                      EXIT
                   end if
                end do
                nCritInOut = iMapToUiqCrit_I(iCritPhy)
             end if
             call SetCritArea(3,nCritArgs,CritName_I,iCritPhy)
          case('dx')
             ! at this time we do not know corect index for dx
             iMapToUiqCrit_I(nCrit-nCritDxLevel) = 10001
             call SetCritArea(2,nCritArgs,CritName_I,nCrit-nCritDxLevel)
          case('level') 
             ! at this time we do not know corect index for level
             iMapToUiqCrit_I(nCrit-nCritDxLevel) = 10002
             call SetCritArea(2,nCritArgs,CritName_I,nCrit-nCritDxLevel)
          case default
             iCritPhy = iCritPhy +1
             !else
             if(.not. present(NameCritOut_I) &
                  .and. .not. present(nCritInOut))&
                  call CON_stop(NameCommand//' ERROR: Need a name table')

             ! Find out if the name has bin used before
             do iCritName = 1, nCrit
                if(trim(CritName) == trim(NameCritOut_I(iCritName))) then
                   iMapToUiqCrit_I(iCritPhy) = iCritName
                   IsUniqueCritName = .false.
                   EXIT
                end if
             end do

             ! Add it to the list if its unique
             if(IsUniqueCritName) then
                do iCritName = 1, nCritInOut+1
                   if(NameCritOut_I(iCritName) == "NULL" )then
                      iMapToUiqCrit_I(iCritPhy) = iCritName
                      NameCritOut_I(iCritName) = CritName 
                      EXIT
                   end if
                end do
                nCritInOut = iMapToUiqCrit_I(iCritPhy)
             end if
             call SetCritArea(2,nCritArgs,CritName_I,iCritPhy)
             !end if
          end select
          select case(CritName)
          case('level')
             call read_var('RefineTo',RefineCritPhys_I(nCrit-nCritDxLevel))
             call read_var('CoarsenTo',CoarsenCritPhys_I(nCrit-nCritDxLevel))
             RefineCritPhys_I(nCrit-nCritDxLevel) = -RefineCritPhys_I(nCrit-nCritDxLevel)
             CoarsenCritPhys_I(nCrit-nCritDxLevel) = -(CoarsenCritPhys_I(nCrit-nCritDxLevel)+1)
             MaxLevelCritPhys_I(nCrit-nCritDxLevel) = -MaxLevel 
             if(IsRes) MaxLevelCritPhys_I(nCrit-nCritDxLevel) = 0.0 
             nCritDxLevel = nCritDxLevel+1
          case('dx')
             call read_var('RefineTo',RefineCritPhys_I(nCrit-nCritDxLevel))
             call read_var('CoarsenFrom',CoarsenCritPhys_I(nCrit-nCritDxLevel))
             MaxLevelCritPhys_I(nCrit-nCritDxLevel) = -MaxLevel
             if(IsRes)  MaxLevelCritPhys_I(nCrit-nCritDxLevel) = 0.0
             nCritDxLevel = nCritDxLevel+1
          case default   
             call read_var('CoarsenLimit',CoarsenCritPhys_I(iCritPhy))
             call read_var('RefineLimit',RefineCritPhys_I(iCritPhy))
             if(IsLevel .or. IsRes) call read_var('MaxResolution',MaxLevelCritPhys_I(iCritPhy))
             if(IsLevel) MaxLevelCritPhys_I(iCritPhy) = -1.0*MaxLevelCritPhys_I(iCritPhy)
          end select
       end do

       if(UseErrorCrit)&
            call read_var('AmrWavefilter',cAmrWavefilter)
       if(ReadExtra) &
            call read_var('ReadExtraOut',ReadExtraOut)

       DoCritAmr = .true.
       DoAutoAmr = .true.      
       nExtCrit     = nCritInOut
       nIntCrit = nIntCrit-1
       nAmrCritUsed = nCrit
       nPhysCritUsed = nCrit
       ! makeing the error criteia come before the physics criteria
       where(iMapToUiqCrit_I > 10000)
          iMapToUiqCrit_I = iMapToUiqCrit_I-10000 +nCritInOut
       end where
       where(iMapToUiqCrit_I > 0)
          iMapToUiqCrit_I = iMapToUiqCrit_I + nIntCrit
       end where
       iMapToUiqCrit_I = abs(iMapToUiqCrit_I)

       ! Extra stuff for #AMRCRITERIA
       if(.not. (IsLevel .or. IsRes))then
          if(allocated(iTree_IA)) then
             MaxLevelCritPhys_I(1:nAmrCritUsed) = -maxval(iTree_IA(MaxLevel_,:))
          else
             ! if the grid is not initzilised at this point, we set it to max
             ! value and lett propper threading take care of it
             MaxLevelCritPhys_I(1:nAmrCritUsed) = -MaxLevel 
          end if
       end if

       if(DoTestMe .and. iProc == 0) then
          write(*,*) " nCritInOut          = ", nCritInOut
          write(*,*) " nIntCrit            = ", nIntCrit
          write(*,*) " nAmrCritUsed        = ", nAmrCritUsed
          write(*,*) " NameCritOut_I       = ", NameCritOut_I
          write(*,*) " nPhysCritUsed       = ", nPhysCritUsed
          write(*,*) " iMapToStateVar_I     = ", iMapToStateVar_I
          write(*,*) " iMapToUiqCrit_I      = ", iMapToUiqCrit_I
          write(*,*) " CoarsenCritPhys_I   = ", CoarsenCritPhys_I
          write(*,*) " RefineCritPhys_I    = ", RefineCritPhys_I
          write(*,*) " nCritDxLevel        = ", nCritDxLevel
          write(*,*) " nAreaPerCritPhys_I  = ", nAreaPerCritPhys_I
          if(IsLevel .or. IsRes) write(*,*) " MaxLevelCritPhys_I  = ", MaxLevelCritPhys_I
          write(*,*) " AmrWavefilter       = ",cAmrWavefilter
       end if
    case("#AMR") ! compatibilety with old system with BATSRUS amr options
       call read_var('PercentCoarsen', PercentCoarsen)
       call read_var('PercentRefine' , PercentRefine)
       call read_var('MaxTotalBlock',  MaxTotalBlock) 
       DoSortAmrCrit = PercentCoarsen > 0.0 .or. PercentRefine > 0.0       
    case("#AMRLIMIT")
       !if(UseAmrMask) &
       !     call CON_stop(NameCommand//' ERROR: AMR masking not supported')
       call read_var('PercentCoarsen', PercentCoarsen)
       call read_var('PercentRefine' , PercentRefine)
       call read_var('MaxTotalBlock',  MaxTotalBlock) 
       call read_var('DiffCriteriaLevel',  DeltaCritera)
       DoSortAmrCrit = PercentCoarsen > 0.0 .or. PercentRefine > 0.0
    case("#DOAMR") 
       call read_var('DoAmr',DoAmr)
       if(DoAmr) then
          call read_var('DnAmr',DnAmr)
          call read_var('DtAmr',DtAmr)
          call read_var('IsStrictAmr'  ,DoStrictAmr)
       end if
       if(.not. DoStrictAmr) DoSortAmrCrit = .true.
    case default
       call CON_stop(NameSub//'incorect PARAM.in!')
    end select
   
    DoStrictAmr = .not. DoSortAmrCrit 
    !=================================================================
  contains 
    subroutine SetCritArea(iStartElmt, nElmt,Name_I,iCrit)

      integer, intent(in) :: iStartElmt, nElmt, iCrit
      character(len=nChar), intent(in):: Name_I(nMaxComponents)

      integer i,iElmt
      !---------------------------------------------------------------
      if(iStartElmt<= nElmt) then
         i=0
         do iElmt=iStartElmt, nElmt
            i=i+1
            AreaNamesPhys_II(i,iCrit) = trim(adjustl(Name_I(iElmt))) 
         end do
      else
         i=1
         AreaNamesPhys_II(1,iCrit) = "ALL"
      end if
      nAreaPerCritPhys_I(iCrit) = i
    end subroutine SetCritArea


  end subroutine read_amr_criteria

  !============================================================================
  subroutine set_amr_geometry(iBlock, user_amr_geometry)

    use BATL_grid, ONLY: Xyz_DGB
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,nDim

    interface
       subroutine user_amr_geometry(iBlock, iArea, DoRefine)
         integer, intent(in) :: iBlock, iArea
         logical,intent(out) :: DoRefine
       end subroutine user_amr_geometry
    end interface
    optional  :: user_amr_geometry

    integer, intent(in) :: iBlock

    integer :: i,j,k,iAmrBox,iCrit,iArea,idx
    real    :: Xyz_D(nDim) 
    logical :: UseBlock
    !--------------------------------------------------------

    if(nAmrCritUsed < 1) RETURN

    ! Find if Criteria should be used in block
    do iCrit=1,nAmrCritUsed
       UseCrit_IB(iCrit,iBlock) = .true.
       do iArea =1,nAreaPerCritAll_I(iCrit)
          idx = iAreaIdx_II(iArea,iCrit)

          call apply_amr_geometry(iBlock, AreaGeo_I(abs(idx)),&
               UseBlock,DoCalcCritIn= (iCrit==1 .and. iArea==1),&
               user_amr_geometry=user_amr_geometry)

          if(idx < 0) UseBlock = .not. UseBlock
          UseCrit_IB(iCrit,iBlock) = UseCrit_IB(iCrit,iBlock) .and. UseBlock

       end do
    end do

  end subroutine set_amr_geometry
  !============================================================================
  logical function masked_amr_criteria(iBlock,iCritExtIn)

    integer, intent(in) :: iBlock
    integer, optional, intent(in) :: iCritExtIn
    integer :: iCrit
    !--------------------------------------------------------

    if(present(iCritExtIn)) then
      masked_amr_criteria = .false.
      ! loop over physics criterias
      do iCrit=nIntCrit+1,nIntCrit+nPhysCritUsed-nCritDxLevel
        if(iVarCritAll_I(iCrit) == nIntCrit+iCritExtIn) &
           masked_amr_criteria = masked_amr_criteria .or. UseCrit_IB(iCrit,iBlock)
      end do
      masked_amr_criteria = .not. masked_amr_criteria
    else
      masked_amr_criteria = .not. any(UseCrit_IB(nIntCrit+1:nIntCrit+nPhysCritUsed-nCritDxLevel,iBlock))
    end if 

  end function masked_amr_criteria
  !============================================================================

  subroutine clean_amr_criteria

    if(allocated(MaxLevelCritPhys_I)) deallocate(MaxLevelCritPhys_I)
    if(allocated(iMapToUiqCrit_I)) deallocate(iMapToUiqCrit_I)
    if(allocated(iMapToStateVar_I)) deallocate(iMapToStateVar_I)
    if(allocated(iNode_I)) deallocate(iNode_I)

    nExtCrit = 1
    nIntCrit = 0
    DeltaCritera = 1.0e-8
    isNewPhysParam = .true.

    call clean_amr_geometry

  end subroutine clean_amr_criteria

  !============================================================================
  subroutine test_amr_criteria

!!$    use BATL_size, ONLY : MaxBlock,nBlock, iRatio, jRatio, kRatio
!!$    use BATL_mpi, ONLY: iProc, nProc
!!$    use BATL_tree, ONLY: init_tree, set_tree_root, find_tree_node, &
!!$         refine_tree_node, distribute_tree, clean_tree, Unused_B, iNode_B, &
!!$         nNode, show_tree, iStatusNew_A, &
!!$         Coarsen_, Unset_, adapt_tree, move_tree, distribute_tree
!!$    use BATL_grid, ONLY: init_grid, create_grid, clean_grid, CellSize_DB
!!$    use BATL_geometry, ONLY: init_geometry
!!$    use BATL_amr, ONLY: do_amr, init_amr
!!$    use BATL_size, ONLY: MaxDim, nDim, MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
!!$         nI, nJ, nK 
!!$    ! For Random generation
!!$    integer :: jSeed
!!$    logical :: IsFirst
!!$    integer, parameter :: iMPLIER=16807, &
!!$         iMODLUS=2147483647, &
!!$         iMOBYMP=127773, &
!!$         iMOMDMP=2836
!!$
!!$    integer, parameter:: MaxBlockTest            = 50
!!$    integer, parameter:: nRootTest_D(MaxDim)     = (/3,3,3/)
!!$    logical, parameter:: IsPeriodicTest_D(MaxDim)= .false.
!!$    real, parameter:: DomainMin_D(MaxDim) = (/ -24.0, -24.0, -24.0 /)
!!$    real, parameter:: DomainMax_D(MaxDim) = (/ 24.0, 24.0, 24.0 /)
!!$    integer :: iNode, iBlock
!!$
!!$    real, allocatable :: Criterias_IB(:,:),AllCriterias_IBP(:,:,:)
!!$    real, allocatable :: PreCriterias_IB(:,:)
!!$    real, allocatable :: RefineLevel_I(:), CoursenLevel_I(:)
!!$    logical, allocatable :: Used_GB(:,:,:,:)
!!$
!!$    integer :: nCritExt = 4
!!$    integer, allocatable :: iA_I(:)
!!$    real, allocatable :: TestState_VGB(:,:,:,:,:)
!!$    integer :: nVar =3
!!$    integer :: i,j,k,iVar
!!$    logical:: DoTestMe
!!$    character(len=*), parameter :: NameSub = 'test_amr_criteria'
!!$    !-----------------------------------------------------------------------
!!$    DoTestMe = iProc == 0
!!$
!!$    write(*,*) " Temural not testing :: test_amr_criteria"
!!$    RETURN
!!$
!!$    if(DoTestMe) write(*,*) 'Starting ',NameSub
!!$
!!$    call init_tree(MaxBlockTest)
!!$    call init_geometry( IsPeriodicIn_D = IsPeriodicTest_D(1:nDim) )
!!$    call init_grid( DomainMin_D(1:nDim), DomainMax_D(1:nDim) )
!!$    call set_tree_root( nRootTest_D(1:nDim))
!!$
!!$    call find_tree_node( (/0.5,0.5,0.5/), iNode)
!!$    call refine_tree_node(iNode)
!!$    call distribute_tree(.true.)
!!$    call create_grid
!!$    call init_amr
!!$
!!$    call srand(123456789+iProc)
!!$
!!$    allocate(Criterias_IB(nCritExt,nBlock), &
!!$         AllCriterias_IBP(nCritExt,nBlock,nProc),&
!!$         PreCriterias_IB(nCritExt,nBlock))
!!$    allocate(RefineLevel_I(nCritExt),CoursenLevel_I(nCritExt))
!!$
!!$
!!$    !===================== Begin Test  =======================
!!$    DoSortAmrCrit = .true.
!!$    !--------------------- internal --------------------------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    allocate(CoarsenCritAll_I(nVar), &
!!$         RefineCritAll_I(nVar),iVarCritAll_I(nVar),iMapToStateVar_I(nVar))
!!$
!!$    allocate(TestState_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$    nAmrCritUsed = 1
!!$    nAmrCrit = nIntCrit + nExtCritUsed
!!$    !allocate(AmrCrit_IB(nIntCrit+nExtCritUsed,nBlock))
!!$    !AmrCrit_IB = 0.0
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    !do iBlock = 1, nBlock
!!$    !   if(Unused_B(iBlock)) CYCLE
!!$    !   write(*,'(i6,f16.12)') iNode_B(iBlock),AmrCrit_IB(1,iBlock)
!!$    !end do
!!$    ! 
!!$    !do iNode = 1, nNode-1
!!$    !   print *,iNode, iRank_A(iNode)
!!$    !end do
!!$
!!$    do iNode = 2, nNode-1
!!$       if(iRank_A(iNode-1) > iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " Internal test"
!!$       end if
!!$    end do
!!$
!!$    !deallocate(AmrCrit_IB)
!!$
!!$
!!$
!!$    !-------------------- external --------------------------
!!$    Criterias_IB = 0.0
!!$    RefineLevel_I =  1.0
!!$    CoursenLevel_I = -1.0
!!$    ! external criteria copyed into the global criteia
!!$    CoarsenCritAll_I = CoursenLevel_I 
!!$    RefineCritAll_I  = RefineLevel_I
!!$
!!$    nExtCritUsed = 1
!!$    nIntCrit = 0
!!$    nCritExt = nExtCritUsed
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) then
!!$          Criterias_IB(1,iBlock) = 10.0
!!$       else
!!$          Criterias_IB(1,iBlock) = AmrCrit_IB(1,iBlock)
!!$       end if
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB, &
!!$         nCritExt, Criterias_IB)
!!$
!!$    do iNode = 2, nNode-1
!!$       if(iRank_A(iNode-1)>iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " External=Intenal test"
!!$       end if
!!$    end do
!!$
!!$    !--------------------- internal with masked cells --------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    UseAmrMask = .true.
!!$    nAmrBox = 1
!!$    allocate(AmrBox_DII(3,2,nAmrBox))
!!$    AmrBox_DII(1,1,1) = DomainMin_D(1) 
!!$    AmrBox_DII(1,2,1) = DomainMin_D(1) + CellSize_DB(1,1)*nI
!!$    AmrBox_DII(2,1,1) = DomainMin_D(2) 
!!$    AmrBox_DII(2,2,1) = DomainMin_D(2) + CellSize_DB(2,1)*nJ
!!$    AmrBox_DII(3,1,1) = DomainMin_D(3) 
!!$    AmrBox_DII(3,2,1) = DomainMin_D(3) + CellSize_DB(3,1)*nK
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       if( iNode_B(iBlock) == 1) then
!!$          if( AmrCrit_IB(1,iBlock)  == 0.0) &
!!$               write (*,*) " ERROR in ",NameSub, &
!!$               " in  Internal test masked cells", &
!!$               " AmrCrit_IB of Node == 1 shoud be none zero"
!!$       else
!!$          if( AmrCrit_IB(1,iBlock)  /= 0.0) &
!!$               write (*,*) " ERROR in ",NameSub, &
!!$               " in  Internal test masked cells", &
!!$               " AmrCrit_IB of Node /= 1 shoud be zero"
!!$       end if
!!$    end do
!!$
!!$   ! Using any becouse of the ghost cells
!!$   ! do iBlock = 1, nBlock
!!$   !    if(Unused_B(iBlock)) CYCLE
!!$   !    write(*,*) iNode_B(iBlock), &
!!$   !         any(DoAmr_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,iBlock))
!!$   ! end do
!!$
!!$
!!$    UseAmrMask = .false.
!!$    deallocate(AmrBox_DII)
!!$    deallocate(DoAmr_GB)
!!$    !--------------------- internal with masked cells --------
!!$
!!$    !--------------------- internal with masked body -------
!!$    nIntCrit = 1
!!$    nExtCritUsed = 0
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$
!!$    allocate(Used_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
!!$    Used_GB = .true.
!!$
!!$    CoarsenCritAll_I = -1.0
!!$    RefineCritAll_I  =  1.0
!!$    TestState_VGB = 1.0
!!$    iVarCritAll_I(nIntCrit)=1
!!$    iMapToStateVar_I(nIntCrit)=1
!!$
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    do iBlock = 1, nBlock
!!$       if(iNode_B(iBlock) /= nNode)  CYCLE
!!$
!!$       Used_GB(1,1,1,iBlock) = .false.
!!$
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = dexp(0.1*(i*0.5))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       TestState_VGB(:,1,1,1,iBlock) = 1.0e18
!!$    end do
!!$
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB,Used_GB=Used_GB)
!!$
!!$    if(iRank_A(1) /= nNode) &
!!$         write(*,*) " ERROR in ",NameSub, " Internal test masked body"
!!$
!!$    do iNode = 3, nNode-1
!!$       if(iRank_A(iNode-1) > iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " Internal test masked body"
!!$       end if
!!$    end do
!!$
!!$
!!$
!!$    deallocate(Used_GB)
!!$    !-------------------- internal x 2 -------------------
!!$
!!$    if(iRatio > 1 .and. jRatio > 1 .and. kRatio >1 ) then
!!$
!!$       nCritExt = 0
!!$       nIntCrit = 2
!!$       nExtCritUsed = nCritExt
!!$       nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$       iVarCritAll_I(1:nIntCrit)=(/ 1,2 /)
!!$       iMapToStateVar_I =  iVarCritAll_I
!!$       ! external criteria copyed into the global criteia
!!$       CoarsenCritAll_I = CoursenLevel_I 
!!$       RefineCritAll_I  = RefineLevel_I
!!$
!!$       do iBlock=1,nBlock
!!$          if(Unused_B(iBlock)) CYCLE
!!$          do k = MinK, MaxK
!!$             do j = MinJ, MaxJ
!!$                do i = MinI,MaxI
!!$                   TestState_VGB(1,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(35-iNode_B(iBlock)+1)))
!!$
!!$                   TestState_VGB(2,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$
!!$       do iBlock = 1, nBlock
!!$          if(Unused_B(iBlock)) then
!!$             Criterias_IB(1,iBlock) = 10.0
!!$          else
!!$             Criterias_IB(1,iBlock) = AmrCrit_IB(1,nBlock-iBlock+1)
!!$          end if
!!$
!!$       end do
!!$
!!$
!!$       call set_amr_criteria(nVar,TestState_VGB, &
!!$            nCritExt, Criterias_IB)
!!$
!!$       allocate(iA_I(nNode-1))
!!$
!!$       iA_I =(/ 1,35, 2, 34, 33,  3,  4, 32, 31,  5,  6, 30, 29,  7, &
!!$            8, 28, 27,  9, 26, 10, 25, 11, 12, 24, 13, 23, 22, 15, 21, &
!!$            16, 17, 20, 19, 18 /)
!!$
!!$       !do iBlock = 1, nBlock
!!$       !   if(Unused_B(iBlock)) CYCLE
!!$       !   print *,"Criterias_IB(1,iBlock) = ", &
!!$       !        AmrCrit_IB(1:2,iBlock), iNode_B(iBlock)
!!$       !end do
!!$
!!$       
!!$       !do iNode = 1, nNode-1
!!$       !   print *,iRank_A(iNode)," :: ", iA_I(iNode) 
!!$       !end do
!!$
!!$       do iNode = 1, nNode-1, 2
!!$          if(iRank_A(iNode) /= iA_I(iNode)) &
!!$               write(*,*) " ERROR in ",NameSub, "2 x Intenal test"
!!$       end do
!!$
!!$       deallocate(iA_I)
!!$    end if
!!$    !-------------------- testing levels ---------------------
!!$    !intenals
!!$    CoarsenCritAll_I(1) = -1.0
!!$    RefineCritAll_I(1)  =  1.0
!!$    !externals
!!$    RefineLevel_I(1) =  0.030
!!$    CoursenLevel_I(1) = 0.020
!!$
!!$
!!$    ! external criteria copyed into the global criteia
!!$    CoarsenCritAll_I(2) = CoursenLevel_I(1) 
!!$    RefineCritAll_I(2)  = RefineLevel_I(1)
!!$
!!$    PercentRefine  = 30.0
!!$    PercentCoarsen = 10.0
!!$    nCritExt = 1
!!$    nIntCrit = 1
!!$
!!$    nExtCritUsed = nCritExt
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$
!!$    ! init Internals
!!$    iVarCritAll_I(nIntCrit) = 1
!!$    iMapToStateVar_I(nIntCrit) = 1
!!$    do iBlock=1,nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                !do iVar=1,nVar
!!$                !print *,nrand()
!!$                TestState_VGB(2,i,j,k,iBlock) = &
!!$                     dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                !end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    !init Externals
!!$    do iBlock = 1, nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       Criterias_IB(1,iBlock) = 0.001*(nNode - iNode_B(iBlock)+1)
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB, &
!!$         nCritExt, Criterias_IB)
!!$
!!$    do k = nNodeCoarsen+1, nNode-1-nNodeRefine
!!$       if( iRank_A(k) < nNodeRefine .and. & 
!!$            iRank_A(k) > nNode-1 - nNodeCoarsen ) &
!!$            write(*,*) "Error in seting refinment by criteria"
!!$    end do
!!$
!!$    ! test unused block
!!$    Criterias_IB = 0.0
!!$    RefineLevel_I =  1.0
!!$    CoursenLevel_I = -1.0
!!$    nCritExt = 0
!!$    nIntCrit = 1
!!$
!!$    iStatusNew_A = Unset_
!!$    nExtCritUsed = nCritExt
!!$    nAmrCritUsed = nIntCrit + nExtCritUsed 
!!$
!!$    do iBlock=1,nBlock
!!$       if(iNode_B(iBlock) > 3**nDim) then 
!!$          iStatusNew_A(iNode_B(iBlock)) =  Coarsen_
!!$       end if
!!$    end do
!!$    !call show_tree(NameSub,.true.)
!!$    call adapt_tree
!!$    call distribute_tree(DoMove=.false.)
!!$    call do_amr(nVar,TestState_VGB)
!!$    call move_tree
!!$    !call show_tree(NameSub,.true.)
!!$
!!$    do iBlock=1,nBlock
!!$       if(Unused_B(iBlock)) CYCLE
!!$       do k = MinK, MaxK
!!$          do j = MinJ, MaxJ
!!$             do i = MinI,MaxI
!!$                do iVar=1,nVar
!!$                   TestState_VGB(iVar,i,j,k,iBlock) = &
!!$                        dexp(0.1*(i*(iNode_B(iBlock)+1)))
!!$                end do
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    call set_amr_criteria(nVar,TestState_VGB)
!!$
!!$    do iNode = 2, nNode
!!$       if(iRank_A(iNode-1)>iRank_A(iNode)) then
!!$          write(*,*) " ERROR in ",NameSub, " unused  test"
!!$       end if
!!$    end do
!!$
!!$
!!$    !===================== End Test  =======================
!!$    deallocate(CoarsenCritAll_I,RefineCritAll_I,iVarCritAll_I,iMapToStateVar_I)
!!$    deallocate(TestState_VGB)
!!$
!!$    deallocate(Criterias_IB,PreCriterias_IB,AllCriterias_IBP)
!!$    deallocate(RefineLevel_I,CoursenLevel_I)
!!$    call clean_grid
!!$    call clean_tree
!!$
!!$  contains
!!$
!!$    ! The saudo random number generator is for testing performense in
!!$    ! parallel sorting
!!$    subroutine srand(iSeed)
!!$      integer, intent(in) :: iSeed
!!$      jSeed = iSeed
!!$      IsFirst = .true.
!!$    end subroutine srand
!!$
!!$    real function rand()
!!$      !  A pseudo-random number generator implemented to make sure that 
!!$      !  all platform reproduce the same sequence for testing and compering.
!!$      !  The algorithm is based on "Integer Version 2" given in :
!!$      !
!!$      !       Park, Steven K. and Miller, Keith W., "Random Number Generators: 
!!$      !       Good Ones are Hard to Find", Communications of the ACM, 
!!$      !       October, 1988.
!!$
!!$      integer :: nHvalue,nLvalue,nTestv
!!$      integer, save :: nExtn
!!$
!!$      if(IsFirst) then
!!$         nExtn=jSeed
!!$         IsFirst = .false.
!!$      end if
!!$
!!$      nHvalue = nExtn/iMOBYMP
!!$      nLvalue = mod(nExtn,iMOBYMP)
!!$      nTestv = iMPLIER*nLvalue - iMOMDMP*nHvalue
!!$      if(nTestv > 0) then
!!$         nExtn = nTestv
!!$      else
!!$         nExtn = nTestv + iMODLUS
!!$      end if
!!$
!!$      rand = real(nExtn)/real(iMODLUS)
!!$
!!$    end function rand

  end subroutine test_amr_criteria

end module BATL_amr_criteria
