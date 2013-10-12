!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan

module ModCT

  use ModSize
  use ModIO,         ONLY: iUnitOut, write_prefix
  use ModProcMH,     ONLY: iProc

  implicit none
  SAVE

  !\
  ! Variables for Constrained Transport
  !/

  ! Face centered magnetic field components, Bxface is centered on face X etc.

  real, allocatable :: Bxface_BLK(:,:,:,:)
  real, allocatable :: Byface_BLK(:,:,:,:)
  real, allocatable :: Bzface_BLK(:,:,:,:)

  ! Face centered normal magnetic field from the 4 finer neighbors
  ! on the two sides of the block (1=east/south/bot, 2=west/north/top)

  real, allocatable :: BxFaceFine_XQSB(:,:,:,:,:)
  real, allocatable :: ByFaceFine_YQSB(:,:,:,:,:)
  real, allocatable :: BzFaceFine_ZQSB(:,:,:,:,:)

  ! VxB stored at edges

  real, allocatable :: VxB_x(:,:,:,:)
  real, allocatable :: VxB_y(:,:,:,:)
  real, allocatable :: VxB_z(:,:,:,:)

  logical :: DoInitConstrainB = .true.

contains
  !============================================================================
  subroutine init_mod_ct


    if(allocated(Bxface_BLK)) return
    allocate(Bxface_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    allocate(Byface_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    allocate(Bzface_BLK(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    allocate(BxFaceFine_XQSB(nJ,nK,4,2,nBLK))
    allocate(ByFaceFine_YQSB(nI,nK,4,2,nBLK))
    allocate(BzFaceFine_ZQSB(nI,nJ,4,2,nBLK))
    allocate(VxB_x(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    allocate(VxB_y(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    allocate(VxB_z(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nBLK))
    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'init_mod_ct allocated arrays'
    end if

  end subroutine init_mod_ct
  !============================================================================
  subroutine clean_mod_ct

    if(.not.allocated(Bxface_BLK)) return
    deallocate(Bxface_BLK)
    deallocate(Byface_BLK)
    deallocate(Bzface_BLK)
    deallocate(BxFaceFine_XQSB)
    deallocate(ByFaceFine_YQSB)
    deallocate(BzFaceFine_ZQSB)
    deallocate(VxB_x)
    deallocate(VxB_y)
    deallocate(VxB_z)

    if(iProc==0)then
       call write_prefix
       write(iUnitOut,'(a)') 'clean_mod_ct deallocated arrays'
    end if

  end subroutine clean_mod_ct

end module ModCT
