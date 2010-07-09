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

  integer, public            :: nAmrCrit
  real, public, allocatable  :: AmrCrit_IB(:,:)

  ! Threshold limits for refine or unrefined the grid (length nCrit) 
  real                               :: SumCoarsenCrit, SumRefineCrit
  real, allocatable, dimension(:)    :: CoarsenCrit_I, RefineCrit_I, &
       CoarsenAmrCrit_I, RefineCritAll_I, GlobalMaxCritAll_I, MaxCritAll_I 

  integer, allocatable, dimension(:) :: CritVar_I ! Index to variables
  ! used in criteria
  real :: cAmrWavefilter
  real, parameter :: cEpsilon = 1.0d-16 ! avoid zero in denominator
  integer :: nCrit ! Number of variables used for amr internal decisions
  integer :: nAmrCritOld, nBlockOld

contains
  !============================================================================

  subroutine set_amr_criteria(nVar, State_VGB, & 
       nCritExt, CritExt_IB, CoarsenCritExt_I, RefineCritExt_I)

    use BATL_grid, ONLY: CellSize_DB
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
    real, intent(inout), dimension(:,:), optional ::CritExt_IB
    real, intent(in), dimension(:), optional :: CoarsenCritExt_I, RefineCritExt_I

    real, dimension(MaxDim) :: Crit_D

    real :: Crit, maxCrit
    real :: Numerator, Denominator
    integer:: iBlock, iDim, iCrit, i, j, k, iVar
    logical :: DoCoarsen
    integer :: iError
    integer :: nAllreduce
    !------------------------------------------------------------------------
    Crit      = 1.0
    Crit_D    = 0.0
    Numerator   = 0.0
    Denominator = 0.0
    MaxCritAll_I = 0.0
    GlobalMaxCritAll_I =0.0

    if(nCrit > nVar) &
         call CON_stop("set_amr_criteria :: More criteria then variables")

    if(present(nCritExt)) then
       nAmrCrit = nCrit + nCritExt
       if(nAmrCrit /= nAmrCritOld) then
          nAmrCritOld=nAmrCrit
          if(allocated(GlobalMaxCritAll_I)) &
               deallocate(MaxCritAll_I,GlobalMaxCritAll_I, CoarsenAmrCrit_I, RefineCritAll_I)
          allocate(MaxCritAll_I(nAmrCrit), GlobalMaxCritAll_I(nAmrCrit), &
               CoarsenAmrCrit_I(nAmrCrit), RefineCritAll_I(nAmrCrit))
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
          do iCrit = 1, nCrit ! Only variables indexed in CritVar_I will decide the refinement
             iVar = CritVar_I(iCrit)

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
             !             if(present(CritExt_I) ) then 
             !                if( Crit > AmrCrit_IB(iCrit,iBlock)) &
             !                     AmrCrit_IB(iCrit,iBlock) = Crit
             !             else
             !                ! If one cell in the block needs refinement the block will
             !                ! be refined, But only if all cells in the block want to be
             !                ! coarsen the block will be flagged for coarsening
             !
             !                if( Crit > RefineCrit_I(iCrit)) then
             !                   iStatusNew_A(iNode_B(iBlock)) = Refine_
             !                   DoCoarsen = .false.
             !                   CYCLE BLOCK
             !                else if(Crit > CoarsenCrit_I(iCrit)) then
             !                   DoCoarsen = .false.
             !                end if
             !             end if


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
       MaxCritAll_I(iCrit) = maxval(AmrCrit_IB(iCrit,:))
    end do
    if(nProc > 1)then
       call MPI_allreduce(MaxCritAll_I, GlobalMaxCritAll_I,nAmrCrit, MPI_REAL, &
            MPI_MAX,iComm, iError)
    else
       GlobalMaxCritAll_I = MaxCritAll_I
    end if

    ! Normalize the external criteria only
    do iCrit = nCrit+1, nAmrCrit
       if(GlobalMaxCritAll_I(iCrit) /= 0.0) &
            AmrCrit_IB(iCrit,:) = AmrCrit_IB(iCrit,:)/GlobalMaxCritAll_I(iCrit)
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

       ! Decide refinement based mean of normalized error factor  
       Crit = sum(AmrCrit_IB(:,iBlock))
       if( Crit > SumRefineCrit )then
          iStatusNew_A(iNode_B(iBlock)) = Refine_
          DoCoarsen = .false.
          CYCLE BLOCK2
       else if(Crit > SumCoarsenCrit) then
          DoCoarsen = .false.
       end if

       if(DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_

    end do BLOCK2

  end subroutine set_amr_criteria

  !============================================================================
  subroutine  read_amr_criteria_param
    use ModReadParam, ONLY: read_var
    use BATL_size, ONLY: MaxBlock    
    integer :: iCrit
    !-------------------------------------------------------------------------
    nCrit          = 0
    nAmrCrit       = 0
    nAmrCritOld    = 0
    cAmrWavefilter = 1.0e-2
    nBlockOld      = 0
    SumCoarsenCrit = 9999999.9
    SumRefineCrit  = 9999999.9

    call read_var('AmrWavefilter',cAmrWavefilter)  
    call read_var('nCrit', nCrit)
    nAmrCrit = nCrit
    nAmrCritOld = nCrit
    if(allocated(CoarsenCrit_I)) then
       deallocate(CoarsenAmrCrit_I, RefineCritAll_I)
       deallocate(CoarsenCrit_I, RefineCrit_I, CritVar_I,GlobalMaxCritAll_I)
    end if
    allocate(CoarsenCrit_I(nCrit), RefineCrit_I(nCrit),CritVar_I(nCrit))
    allocate(CoarsenAmrCrit_I(nCrit), RefineCritAll_I(nCrit))
    allocate(GlobalMaxCritAll_I(nCrit))
    allocate(MaxCritAll_I(nCrit))
    do iCrit = 1, nCrit
       call read_var('iVar',CritVar_I(iCrit))
       call read_var('CoarsenCrit',CoarsenCrit_I(iCrit))
       call read_var('RefineCrit',RefineCrit_I(iCrit))
    end do
    call read_var('SumCoarsenCrit', SumCoarsenCrit)
    call read_var('SumRefineCrit',  SumRefineCrit)

  end subroutine read_amr_criteria_param
  !============================================================================
  subroutine clean_amr_criteria

    if(.not.allocated(CoarsenCrit_I)) RETURN
    deallocate(CoarsenCrit_I, RefineCrit_I, CritVar_I)
    deallocate(CoarsenAmrCrit_I, RefineCritAll_I)
    deallocate(GlobalMaxCritAll_I)
    deallocate(MaxCritAll_I)
    nCrit    = 0
    nAmrCrit = 0
    nAmrCritOld = 0

  end subroutine clean_amr_criteria

end module BATL_amr_criteria
