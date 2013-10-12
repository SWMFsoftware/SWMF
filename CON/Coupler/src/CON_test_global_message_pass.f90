!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
Module CON_test_global_message_pass
  use CON_global_message_pass
  implicit none
  private !Except
  public::test_global_message_pass
  interface test_global_message_pass
     module procedure test_global_message_pass_dd
     module procedure test_global_message_pass_id
  end interface
contains
  subroutine test_global_message_pass_id(GridID_)
    integer,intent(in)::GridID_
    if(is_proc0())write(*,*)&
         'CON: routine test_global_message_pass is used by developers only'
  end subroutine test_global_message_pass_id
  subroutine test_global_message_pass_dd(DomainDecomposition)
    type(DomainDecompositionType),intent(in)::DomainDecomposition
    call test_global_message_pass_id(1)
  end subroutine test_global_message_pass_dd
end module CON_test_global_message_pass
