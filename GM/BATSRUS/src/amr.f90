!^CFG COPYRIGHT UM
subroutine amr(DoFullMessagePass,TypeAmr)
  use ModProcMH
  use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,&
       unusedBLK, lVerbose, UseB, UseB0, Dt_BLK, nTrueCellsALL, &
       iNewGrid, iNewDecomposition
  use ModGeometry, ONLY : minDXvalue,maxDXvalue,true_cell
  use ModAMR, ONLY : &
       nAmrCriteria, DoProfileAmr, AmrCriteria_IB
  use ModAdvance, ONLY : DivB1_GB, iTypeAdvance_B, iTypeAdvance_BP, &
       nVar, State_VGB, &
       SkippedBlock_ !!!
  use ModBlockData, ONLY: clean_block_data
  use ModIO, ONLY : write_prefix, iUnitOut
  use ModMpi

  use BATL_lib,         ONLY: regrid_batl, set_amr_criteria, &
       MaxNode, nNode, iTree_IA, Status_, Used_, Proc_, Block_ !!!

  use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
  use ModMessagePass,   ONLY: exchange_messages
  use ModPartSteady,    ONLY: UsePartSteady
  use ModUser,          ONLY: user_specify_refinement

  implicit none

  logical, intent(in) :: DoFullMessagePass
  character(3), intent(in) :: TypeAmr

  integer :: iBlock

  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'amr'

  ! Check if we have the same grid as before, store old grid id
  integer, save :: iLastGrid=-1, iLastDecomposition=-1

  integer:: iNode !!!
  integer, allocatable:: iTypeAdvance_A(:) !!!

  !----------------------------------------------------------------------------

  call set_oktest(NameSub, DoTest, DoTestMe)

  ! Do message passing with second order accurate ghost cells

  if(DoTestMe)write(*,*)NameSub,' starts 2nd order accurate message passing'

  if(DoProfileAmr) call timing_start('amr::exchange_true')
  call exchange_messages(UseOrder2In=.true., &
       DoResChangeOnlyIn=.not.DoFullMessagePass)
  if(DoProfileAmr) call timing_stop('amr::exchange_true')

  if(UsePartSteady)then
     ! Convert iTypeAdvance_BP to _A !!! should use _A all the time
     allocate(iTypeAdvance_A(MaxNode))
     do iNode = 1, nNode
        if(iTree_IA(Status_,iNode) /= Used_) CYCLE
        iTypeAdvance_A(iNode) = &
             iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode))
     end do
  end if
  if(nAmrCriteria > 0)then
     AmrCriteria_IB(:,1:nBlockMax) = 0.0
     if(DoProfileAmr) call timing_start('amr::amr_criteria')
     call amr_criteria(AmrCriteria_IB,TypeAmr)
     if(DoProfileAmr) call timing_stop('amr::amr_criteria')
     if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
     call set_amr_criteria(nVar, State_VGB,&
          nAmrCriteria,AmrCriteria_IB,TypeAmrIn=TypeAmr,&
          user_amr_geometry=user_specify_refinement)
     if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')

  else
     if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
     call set_amr_criteria(nVar, State_VGB,TypeAmrIn=TypeAmr,&
          user_amr_geometry=user_specify_refinement)
     if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')
  end if

  if(DoProfileAmr) call timing_start('amr::regrid_batl')
  if(UsePartSteady)then
     call regrid_batl(nVar, State_VGB, Dt_BLK, DoTestIn=DoTestMe, &
          Used_GB=true_cell, iTypeNode_A=iTypeAdvance_A)
  else
     call regrid_batl(nVar, State_VGB, Dt_BLK, DoTestIn=DoTestMe, &
          Used_GB=true_cell)
  end if
  if(DoProfileAmr) call timing_stop('amr::regrid_batl')


  ! This should be eliminated by using iTypeAdvance_A everywhere !!!
  if(UsePartSteady)then
     ! restore iTypeAdvance_B and _BP
     iTypeAdvance_BP = SkippedBlock_
     iTypeAdvance_B  = SkippedBlock_
     do iNode = 1, nNode
        if(iTree_IA(Status_,iNode) /= Used_) CYCLE
        iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode)) = &
             iTypeAdvance_A(iNode)
        if(iTree_IA(Proc_,iNode) == iProc) &
             iTypeAdvance_B(iTree_IA(Block_,iNode)) = iTypeAdvance_A(iNode)
     end do
     deallocate(iTypeAdvance_A)
  end if

  if(DoProfileAmr) call timing_start('amr::set_batsrus_grid')
  call set_batsrus_grid
  if(DoProfileAmr) call timing_stop('amr::set_batsrus_grid')

  ! If the grid has not changed only the message passing has to be redone
  ! to reset ghost cells at resolution changes
  if(iNewGrid==iLastGrid .and. iNewDecomposition==iLastDecomposition) then
     if(DoProfileAmr) call timing_start('amr::exchange_noamr')
     call exchange_messages(DoResChangeOnlyIn=.true., UseOrder2In=.false.)
     if(DoProfileAmr) call timing_stop('amr::exchange_noamr')
     RETURN
  end if

  iLastGrid          = iNewGrid
  iLastDecomposition = iNewDecomposition

  if(DoProfileAmr) call timing_start('amr::count_true_cells')
  call count_true_cells
  if(DoProfileAmr) call timing_stop('amr::count_true_cells')

  ! Clean all dynamically stored block data
  call clean_block_data

  if(iProc==0 .and. lVerbose>0)then
     ! Write block/cell summary after AMR
     call write_prefix; write(iUnitOut,*) '|'
     call write_prefix; write(iUnitOut,*) &
          '|  AMR:  nBlockMax = ',nBlockMax,' nBLK = ',nBLK
     call write_prefix; write(iUnitOut,*) &
          '|  AMR:  Total number of blocks used = ', nBlockALL
     call write_prefix; write(iUnitOut,*) &
          '|  AMR:  Total number of cells = ', nBlockALL*nIJK
     call write_prefix; write(iUnitOut,*) &
          '|  AMR:  Total number of true cells = ', nTrueCellsALL
     call write_prefix; write(iUnitOut,*) &
          '|  Smallest cell dx: ', minDXvalue, &
          '  Largest cell dx: ',   maxDXvalue
     call write_prefix; write(iUnitOut,*) '|'
  end if

  if(DoProfileAmr) call timing_start('amr::set_batsrus_state')
  ! Fix energy and other variables in moved/refined/coarsened blocks
  call set_batsrus_state
  if(DoProfileAmr) call timing_stop('amr::set_batsrus_state')

  ! Update iTypeAdvance, and redo load balancing if necessary
  ! Load balance: move coords, data, and there are new blocks
  if(DoProfileAmr) call timing_start('amr::load_balance')
  call load_balance(.true.,.true.,.true.)
  if(DoProfileAmr) call timing_stop('amr::load_balance')
  ! redo message passing
  if(DoProfileAmr) call timing_start('amr::exchange_false')
  call exchange_messages(UseOrder2In=.false.)
  if(DoProfileAmr) call timing_stop('amr::exchange_false')

  if(UseB0)then
     ! Correct B0 face at newly created and removed resolution changes
     if(DoProfileAmr) call timing_start('amr::set_b0_source')
     do iBlock = 1, nBlock
        if (unusedBLK(iBlock)) CYCLE
        call set_b0_source(iBlock)
     end do
     if(DoProfileAmr) call timing_stop('amr::set_b0_source')
  end if

  ! Reset divb (it is undefined in newly created/moved blocks)
  if(UseB)then
     if(DoProfileAmr) call timing_start('amr::set_divb')
     DivB1_GB(:,:,:,1:nBlock) = -7.70
     if(DoProfileAmr) call timing_stop('amr::set_divb')
  end if

end subroutine amr
