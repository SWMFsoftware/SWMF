!^CFG COPYRIGHT UM
subroutine amr(idepth)
  use ModProcMH
  use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,MaxBlock,&
       unusedBLK,lVerbose,UseB,UseB0
  use ModGeometry, ONLY : minDXvalue,maxDXvalue,dx_BLK
  use ModAMR, ONLY : automatic_refinement
  use ModAdvance, ONLY : DivB1_GB, iTypeAdvance_B, iTypeAdvance_BP
  use ModBlockData, ONLY: clean_block_data
  use ModIO, ONLY : write_prefix, iUnitOut
  use ModMpi
  implicit none

  integer, intent(in) :: idepth

  integer :: iBlock, iError
  real :: minDX,maxDX    
  logical :: local_refine(nBLK)

  !--------------------------------------------------------------------
  ! Ensure ghostcells are up to date.
  !
  call exchange_messages

  if (automatic_refinement) then
     !-----------------------------------------------------------------
     ! Physics based refinement.
     !
     call amr_physics

  else
     !-----------------------------------------------------------------
     ! Prespecified refinement.
     !
     call specify_refinement(local_refine)
     call refine_grid(local_refine)

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
  if(iProc==0.and.lVerbose>0)then
     ! Write block/cell summary after AMR
     call write_prefix; write(iUnitOut,*) '|'
     call write_prefix; 
     write(iUnitOut,*) '|  AMR:  nBlockMax = ',nBlockMax,' nBLK = ',nBLK
     call write_prefix; 
     write(iUnitOut,*) '|  AMR:  Total number of blocks used = ', nBlockALL
     call write_prefix; 
     write(iUnitOut,*) '|  AMR:  Total number of cells = ', nBlockALL*nIJK
     call write_prefix; 
     write(iUnitOut,*) '|  Smallest cell dx: ',minDXvalue,&
          '  Largest cell dx: ',maxDXvalue
     call write_prefix; write(iUnitOut,*) '|'
  end if

  ! Update ghost cells
  call exchange_messages
  if(UseB0)then
     ! Correct B0 face at newly created and removed resolution changes
     do iBlock=1,nBlock
        if (unusedBLK(iBlock)) CYCLE
        call set_b0_face(iBlock)
     end do
  end if
  ! Reset divb (it is undefined in newly created/moved blocks)
  if(UseB)DivB1_GB=-7.70

end subroutine amr
