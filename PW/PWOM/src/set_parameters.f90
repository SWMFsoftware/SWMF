subroutine PW_set_parameters(NameAction)
  use Mod_PW, only: iUnitOut
  implicit none
  character (len=*), intent(in) :: NameAction
  character (len=*), parameter :: NameSub = 'PW_set_parameters'
  
  write(iUnitOut,*) NameSub,': called with action=',NameAction


end subroutine PW_set_parameters
