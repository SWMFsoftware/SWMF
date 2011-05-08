!^CFG COPYRIGHT UM
subroutine allocate_vars
  use ModProcMH
  use ModMain, ONLY : unusedBLK, UseBatl
  use ModAMR
  use ModParallel
  use ModMpi
  implicit none

  integer :: iPE, iBLK, ierror
  integer :: MaxBlockALL

  if(.not. UseBatl)then
     !\  
     ! Allocate and initialize list of available solution blocks.
     !/
     allocate( availableBLKs(0:nBLK,0:nProc-1), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: availableBLKs"
     availableBLKs(0,:)=1
     do iBLK=1,nBLK; do iPE=0,nProc-1
        availableBLKs(iBLK,iPE)=iBLK
     end do; end do

     !\
     ! Allocate and initialize global refinement/coarsening criteria and
     ! flagging variables.
     !/
     allocate( refine_criteria_list(4,nBLK,nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: refine_criteria_list"
     refine_criteria_list = 0.00

     allocate( refine_list(nBLK,nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: refine_list"
     refine_list = .false.

     allocate( coarsen_list(nBLK,nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: coarsen_list"
     coarsen_list = .false.

     allocate( ListToCoarsen(2,nBLK*nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: ListToCoarsen"
     ListToCoarsen = -1

     allocate( ListToRefine(2,nBLK*nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: ListToRefine"
     ListToRefine = -1

     allocate( SortB(3,nBLK*nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: SortB"
     SortB = -1

     allocate( SortP(3,nBLK*nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: SortP"
     SortP = -1

     allocate( SortC(3,nBLK*nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: SortC"
     SortC = 0.00

     !\
     ! Allocate and initialize other global variables.
     !/
     allocate( PE_time_list(nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: PE_time_list"
     PE_time_list = 0.00

     allocate( PE_float_list(nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: PE_float_list"
     PE_float_list = 0.00

     allocate( PE_mflop_list(nProc), stat=ierror )
     if (ierror > 0) write(*,*) "allocate_vars: PE = ",iProc, &
          " allocation error: PE_mflop_list"
     PE_mflop_list = 0.00
  end if

  !\
  ! allocate index arrays
  !/

  MaxBlockALL = MaxBlock * nProc
  allocate(iBlock_A(MaxBlockALL), iProc_A(MaxBlockALL), &
       iBlockRestartALL_A(MaxBlockALL))

  !\
  ! allocate unusedBlock_BP and set the default false value
  !/
  allocate(unusedBlock_BP(MaxBlock, 0:nProc-1))

  !\  
  ! Initialize unused block indicators.
  !/  
  unusedBLK     = .true.
  unusedBlock_BP= .false.

end subroutine allocate_vars
