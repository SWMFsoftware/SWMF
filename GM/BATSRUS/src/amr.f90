!^CFG COPYRIGHT UM
subroutine amr(idepth)
  use ModProcMH
  use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,unusedBLK,lVerbose
  use ModGeometry, ONLY : minDXvalue,maxDXvalue,dx_BLK
  use ModAMR, ONLY : automatic_refinement
  use ModAdvance, ONLY : DivB1_GB
  use ModIO, ONLY : write_prefix, iUnitOut
  use ModMpi
  implicit none

  integer, intent(in) :: idepth

  integer :: iBlock, nBlockMoved, iError
  real :: minDX,maxDX    
  logical :: local_refine(nBLK)

  !--------------------------------------------------------------------
  ! Ensure ghostcells are up to date.
  !
  call exchange_messages

  if (automatic_refinement) then              !^CFG IF NOT SIMPLE
     !-----------------------------------------------------------------
     ! Physics based refinement.
     !
     call amr_physics

  else                                        !^CFG IF NOT SIMPLE BEGIN
     !-----------------------------------------------------------------
     ! Prespecified refinement.
     !
     call specify_initial_refinement(local_refine, idepth)
     call refine_grid(local_refine)

  end if                                      !^CFG END SIMPLE
  call fixRefinementLevels

  ! Find new min and max dx
  minDX=minval(dx_BLK, MASK=(.not.unusedBLK))
  maxDX=maxval(dx_BLK, MASK=(.not.unusedBLK))
  call MPI_allreduce(minDX,minDXvalue,1,MPI_REAL,MPI_MIN,iComm,iError)
  call MPI_allreduce(maxDX,maxDXvalue,1,MPI_REAL,MPI_MAX,iComm,iError)


  !--------------------------------------------------------------------
  ! Fix up other variables due to refinement
  !
  call number_soln_blocks
  call load_balance(.true.,.true.,nBlockMoved)
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

  call find_neighbors
  call exchange_messages
  ! Correct B0 face at newly created and removed resolution changes
  do iBlock=1,nBlock
     if (unusedBLK(iBlock)) CYCLE
     call set_b0_face(iBlock)
  end do
  DivB1_GB=-7.70
end subroutine amr
