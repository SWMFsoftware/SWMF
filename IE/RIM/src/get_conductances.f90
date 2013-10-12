!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine get_conductance

  use ModParamRIM

  implicit none

  if (UseUAConductances) then

     write(*,*) "Sorry.... UA is not coupled right now.... "

     ! Fill in conductances here
     ! return

  endif

  


end subroutine get_conductance
