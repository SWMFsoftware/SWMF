!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================!
module SP_ModUnit
  use SP_ModGrid, ONLY: nVar, LagrID_, FluxMax_
  implicit none
  !\
  ! unit of SEP energy is also applicable for ion temperature
  character(len=*), parameter :: NameEUnit = 'kev'
  real                        :: UnitEnergy
  ! simulated particles, used for converting momentum to energy
  ! is back. For different sorts of ions the sqrt(A) factor needs
  ! to be used in these relations.
  character(len=*), parameter :: NameParticle = 'proton'
  !/
  character(len=6), public, parameter:: NameVarUnit_V(LagrID_:FluxMax_) = (/&
       'none  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'amu/m3', &
       NameEUnit//'   ', &
       'm/s   ', &
       'm/s   ', &
       'm/s   ', &
       'T     ', &
       'T     ', &
       'T     ', &
       'J/m3  ', &
       'J/m3  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'm/s   ', &
       'T     ', &
       'none  ', &
       'amu/m3', &
       'T     ', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       'p.f.u.', &
       '??????'  /)
end module SP_ModUnit
