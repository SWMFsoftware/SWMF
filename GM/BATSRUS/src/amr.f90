!^CFG COPYRIGHT UM
subroutine amr(DoMessagePass)
  use ModProcMH
  use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,MaxBlock,&
       unusedBLK,lVerbose,UseB,UseB0, UseBatl, Dt_BLK, nTrueCellsALL
  use ModGeometry, ONLY : minDXvalue,maxDXvalue,dx_BLK,true_cell
  use ModAMR, ONLY : automatic_refinement, RefineLimit_I, CoarsenLimit_I, &
       nRefineCrit
  use ModAdvance, ONLY : DivB1_GB, iTypeAdvance_B, iTypeAdvance_BP, &
       nVar, State_VGB
  use ModBlockData, ONLY: clean_block_data
  use ModIO, ONLY : write_prefix, iUnitOut
  use ModMpi

  use ModParallel,      ONLY: UsePlotMessageOptions
  use BATL_lib,         ONLY: regrid_batl, set_amr_criteria, Unused_B, iNode_B, &
       iStatusNew_A, Refine_, Coarsen_
  use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
  use ModUser,          ONLY: user_amr_criteria
  use ModBatlInterface, ONLY: useBatlTest
  implicit none

  logical, intent(in) :: DoMessagePass

  logical:: IsFound
  real :: UserCriteria
  integer :: iBlock, iError
  real :: minDX, maxDX    
  logical :: DoRefine_B(nBLK)

  logical :: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'amr'

  real :: refine_criteria(4, nBLK)
  !----------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)
  if(UseBatl)then
     ! Do message passing with second order accurate ghost cells

     if(DoTestMe)write(*,*)NameSub,' starts 2nd order accurate message passing'

     if(DoMessagePass)then
        if(.not.useBatlTest) UsePlotMessageOptions = .true.
        call exchange_messages
     end if
     if(automatic_refinement) then

        if(nRefineCrit > 0)then
           refine_criteria = 0.0
           call amr_criteria(refine_criteria)
           call set_amr_criteria(nVar, State_VGB,&
                nRefineCrit,refine_criteria,CoarsenLimit_I, RefineLimit_I)
        else
           call set_amr_criteria(nVar, State_VGB)
        end if

        call regrid_batl(nVar, State_VGB, Dt_BLK, DoTestIn=DoTestMe,Used_GB=true_cell)

        call set_batsrus_grid

        call count_true_cells
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

     else
        call specify_refinement(DoRefine_B)
        call regrid_batl(nVar, State_VGB, Dt_BLK, DoRefine_B, DoTestIn=DoTestMe,Used_GB=true_cell)
        call set_batsrus_grid

        if(iProc == 0 .and. lVerbose>0) then
           call write_prefix; write(iUnitOut,*) &
                'specify_refinement,  new number of blocks = ', nBlockALL
        end if
     end if

     ! Fix energy and other variables in moved/refined/coarsened blocks
     call set_batsrus_state

     if(DoMessagePass)then
        ! Update iTypeAdvance, and redo load balancing if necessary
        ! Load balance: move coords, data, and there are new blocks
        call load_balance(.true.,.true.,.true.)

        ! redo message passing
        UsePlotMessageOptions = .false.
        call exchange_messages
     end if
     RETURN !!! TODO: iTypeAdvance, B0, ModBlockData...
  end if

  ! Ensure ghostcells are up to date.
  ! UsePlotMessageOptions = .true. !!! this would be useful
  call exchange_messages

  if (automatic_refinement) then
     ! Physics based refinement.
     call amr_physics
  else
     ! Prespecified refinement.
     call specify_refinement(DoRefine_B)
     call refine_grid(DoRefine_B)
  end if

  ! Find new min and max dx
  minDX=minval(dx_BLK, MASK=(.not.unusedBLK))
  maxDX=maxval(dx_BLK, MASK=(.not.unusedBLK))
  call MPI_allreduce(minDX,minDXvalue,1,MPI_REAL,MPI_MIN,iComm,iError)
  call MPI_allreduce(maxDX,maxDXvalue,1,MPI_REAL,MPI_MAX,iComm,iError)

  ! Renumber blocks
  call number_soln_blocks

  ! Update the global advance info
  call MPI_allgather(iTypeAdvance_B, MaxBlock, MPI_INTEGER, &
       iTypeAdvance_BP, MaxBlock, MPI_INTEGER, iComm, iError)

  ! Clean all dynamically stored block data
  call clean_block_data

  ! Load balance: move coords, data, and there are new blocks
  call load_balance(.true.,.true.,.true.)
  call count_true_cells
  if(iProc==0.and.lVerbose>0)then
     ! Write block/cell summary after AMR
     call write_prefix
     write(iUnitOut,*) '|'
     call write_prefix 
     write(iUnitOut,*) '|  AMR:  nBlockMax = ',nBlockMax,' nBLK = ',nBLK
     call write_prefix 
     write(iUnitOut,*) '|  AMR:  Total number of blocks used = ', nBlockALL
     call write_prefix 
     write(iUnitOut,*) '|  AMR:  Total number of cells = ', nBlockALL*nIJK
     call write_prefix
     write(iUnitOut,*) '|  AMR:  Total number of true cells = ', nTrueCellsALL
     call write_prefix 
     write(iUnitOut,*) '|  Smallest cell dx: ',minDXvalue,&
          '  Largest cell dx: ',maxDXvalue
     call write_prefix
     write(iUnitOut,*) '|'
  end if

  ! Update ghost cells
  ! UsePlotMessageOptions = .false. !!! this would be useful (see above)
  call exchange_messages
  if(UseB0)then
     ! Correct B0 face at newly created and removed resolution changes
     do iBlock=1,nBlock
        if (unusedBLK(iBlock)) CYCLE
        call set_b0_source(iBlock)
     end do
  end if
  ! Reset divb (it is undefined in newly created/moved blocks)
  if(UseB)DivB1_GB=-7.70

end subroutine amr
