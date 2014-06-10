!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHeidiDrifts
  !\
  ! Drift and diffusion variable definition module for the HEIDI program.
  ! Mike Liemohn, March 2006
  !/
  
  use ModHeidiSize, ONLY: nR, nT, NE, nPa, nS
  
  ! Define the atmospheric loss and charge exchange variables
  integer :: J6 = 0
  integer :: J18 = 0
  
  ! Define the advection drift rate variables
  real, allocatable :: VR(:,:,:,:)
  real :: P1(NR,NT)
  real, allocatable :: P2(:,:,:,:)
  real, allocatable :: ACHARGE(:,:,:,:,:)
  real, allocatable :: ATLOS(:,:,:,:,:)
  real, allocatable :: EDOT(:,:,:,:)
  real, allocatable :: MUDOT(:,:,:,:)
  real, allocatable :: VrConv(:,:,:,:)
  
  ! Define the pitch angle diffusion coefficient variables
  
  real, allocatable :: COULE(:,:,:,:,:)
  real, allocatable :: COULI(:,:,:,:,:)
  real, allocatable :: ATAE(:,:,:,:,:)
  real, allocatable :: BTAE(:,:,:,:,:)
  real, allocatable :: GTAE(:,:,:,:,:)
  real, allocatable :: ATAI(:,:,:,:,:)
  real, allocatable :: BTAI(:,:,:,:,:)
  real, allocatable :: GTAI(:,:,:,:,:)


contains

  subroutine init_mod_heidi_drifts


    if(allocated(VR)) RETURN
    allocate(VR(NR,NT,NE,NPA))
    VR = 0.0
    allocate(EDOT(NR,NT,NE,NPA))
    allocate(MUDOT(NR,NT,NE,NPA))
    allocate(VrConv(NR,NT,NE,NPA))
    allocate(COULE(NR,NT,NE,NPA,NS))
    allocate(COULI(NR,NT,NE,NPA,NS))
    allocate(ATAE(NR,NT,NE,NPA,NS))
    allocate(BTAE(NR,NT,NE,NPA,NS))
    allocate(GTAE(NR,NT,NE,NPA+1,NS))
    allocate(ATAI(NR,NT,NE,NPA,NS))
    allocate(BTAI(NR,NT,NE,NPA,NS))
    allocate(GTAI(NR,NT,NE,NPA+1,NS))

    allocate(P2(NR,NT,NE,NPA))    
    P2 = 0.0
    allocate(ACHARGE(NR,NT,NE,NPA,NS))
    ACHARGE = 0.0
    allocate(ATLOS(NR,NT,NE,NPA,NS))
    ATLOS = 0.0

  end subroutine init_mod_heidi_drifts

  subroutine clean_mod_heidi_drifts


    if(.not.allocated(VR)) RETURN
    deallocate(VR)
    deallocate(EDOT)
    deallocate(MUDOT)
    deallocate(VrConv)
    deallocate(COULE)
    deallocate(COULI)
    deallocate(ATAE)
    deallocate(BTAE)
    deallocate(GTAE)
    deallocate(ATAI)
    deallocate(BTAI)
    deallocate(GTAI)

    deallocate(P2)
    deallocate(ACHARGE)
    deallocate(ATLOS)
    
  end subroutine clean_mod_heidi_drifts

end module ModHeidiDrifts
