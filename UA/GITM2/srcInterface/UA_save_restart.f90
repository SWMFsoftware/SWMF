
subroutine UA_save_restart(TimeSimulation)

  use ModInputs

  implicit none

  real, intent(in) :: TimeSimulation

  call write_restart("UA/restartOUT/")

end subroutine UA_save_restart
