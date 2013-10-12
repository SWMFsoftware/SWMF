!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModIEInterface

  use CON_comp_param, only: MaxComp
  use ModKind

  implicit none

  type iecouplertype

     integer :: nLons
     integer :: nLats
     real(Real8_) :: Time

     real, allocatable :: Longitudes(:,:)
     real, allocatable :: Latitudes(:,:)

     integer, allocatable :: iProcLoc(:,:)
     integer, allocatable :: iLonLoc(:,:)
     integer, allocatable :: iLatLoc(:,:)

     real, allocatable :: rLatInter(:,:)
     real, allocatable :: rLonInter(:,:)

  end type iecouplertype

  type(iecouplertype) :: iecouplerinfo(MaxComp)

end module ModIEInterface
