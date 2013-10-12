!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

subroutine UA_save_restart(TimeSimulation)

  use ModInputs

  implicit none

  real, intent(in) :: TimeSimulation

  call write_restart("UA/restartOUT/")

end subroutine UA_save_restart
