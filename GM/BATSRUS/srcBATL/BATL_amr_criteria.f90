!^CFG COPYRIGHT UM

module BATL_amr_criteria

  implicit none

  SAVE

  private ! except

  public set_amr_criteria
  public clean_amr_criteria
  public read_amr_criteria_param

  ! Threshold limits for refine or unrefined the grid (length nCrit) 
  real, allocatable, dimension(:) :: CoarsenCrit_I, RefineCrit_I
  integer, allocatable, dimension(:) :: CritVar_I ! Index to variables
  ! used in criteria
  real :: cAmrWavefilter
  real, parameter :: cEpsilon = 1.0d-16 ! avoid zero in denominator
  integer :: nCrit ! Number of variables used for amr decisions

contains
  !============================================================================

  subroutine set_amr_criteria(nVar,State_VGB, UseLocalStep)

    use BATL_grid, ONLY: CellSize_DB
    use BATL_size, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK,&
	 nI, nJ, nK, MaxDim, nDim, MaxBlock, nBlock
    use BATL_tree, ONLY: iStatusNew_A, Refine_, Coarsen_, &
         Unused_B, iNode_B


    integer,  intent(in)   :: nVar
    real,    intent(inout):: &                            ! state variables
         State_VGB(nVar,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock)
    logical, intent(in)   :: UseLocalStep

    real, dimension(MaxDim) :: Factor_D

    real :: Factor
    real :: Numerator, Denominator
    integer:: iBlock, iDim, iCrit, i, j, k, iVar
    logical :: DoCoarsen
    !------------------------------------------------------------------------
    Factor      = 1.0
    Factor_D    = 0.0
    Numerator   = 0.0
    Denominator = 0.0

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
    ! The subroutine only set iStatusNew_A array.   

    BLOCK: do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       DoCoarsen = .true.
       ! Check refinement first
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          do iCrit = 1, nCrit ! Only variables indexed in CritVar_I will deside the refinment
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
             Factor_D(1) = (Numerator/max(Denominator,cEpsilon)) 

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
                Factor_D(2) = (Numerator/max(Denominator,cEpsilon))
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
                Factor_D(3) = (Numerator/max(Denominator,cEpsilon))
             end if


       
             Factor = sqrt(sum(Factor_D**2))


             ! If one cell in the block needs refinement the block will
             ! be refined, But only if all cells in the block want to be
             ! coarsen the block will be flagged for coarsening
             
             if( Factor > RefineCrit_I(iCrit)) then
                iStatusNew_A(iNode_B(iBlock)) = Refine_
                DoCoarsen = .false.
                CYCLE BLOCK
             else if(Factor > CoarsenCrit_I(iCrit)) then
                DoCoarsen = .false.
             end if

          end do !end nVAr
       end do; end do; end do

       if( DoCoarsen) iStatusNew_A(iNode_B(iBlock)) =  Coarsen_

    end do BLOCK

  end subroutine set_amr_criteria

  !============================================================================
  subroutine  read_amr_criteria_param
    use ModReadParam, ONLY: read_var
    integer :: iCrit
    !-------------------------------------------------------------------------
    call read_var('AmrWavefilter',cAmrWavefilter)  
    call read_var('nCrit', nCrit)
    if(allocated(CoarsenCrit_I)) &	
         deallocate(CoarsenCrit_I, RefineCrit_I, CritVar_I)
    allocate(CoarsenCrit_I(nCrit), RefineCrit_I(nCrit),CritVar_I(nCrit))
    do iCrit = 1, nCrit
       call read_var('iVar',CritVar_I(iCrit))
       call read_var('CoarsenCrit',CoarsenCrit_I(iCrit))
       call read_var('RefineCrit',RefineCrit_I(iCrit))
    end do

  end subroutine read_amr_criteria_param
  !============================================================================
  subroutine clean_amr_criteria

    if(.not.allocated(CoarsenCrit_I)) RETURN
    deallocate(CoarsenCrit_I, RefineCrit_I, CritVar_I)
    nCrit = 0

  end subroutine clean_amr_criteria

end module BATL_amr_criteria
