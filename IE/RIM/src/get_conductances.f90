
subroutine get_conductance

  use ModParamRIM

  implicit none

  if (UseUAConductances) then

     write(*,*) "Sorry.... UA is not coupled right now.... "

     ! Fill in conductances here
     ! return

  endif

  


end subroutine get_conductance
