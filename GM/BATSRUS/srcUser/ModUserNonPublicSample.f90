!#NOTPUBLIC  email:darrens@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!==============================================================================
module ModUser

  ! This file contains a sample of a protected ModUser file.
  ! The top line above must follow the format shown.

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'NOTPUBLIC Protected Sample'

contains

  !============================================================================
  subroutine user_read_inputs

    use ModReadParam, ONLY: read_line, read_command, read_var
    character(len=100):: NameCommand
    real:: MyVariable

    character(len=*), parameter :: NameSub = "ModUser::user_read_inputs"
    !--------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case('#MYCOMMAND')
          call read_var('MyVariable', MyVariable)
       case('#USERINPUTEND')
          EXIT
       case default
          call stop_mpi(NameSub//': unknown command name='//trim(NameCommand))
       end select
    end do
  end subroutine user_read_inputs

  !============================================================================

end module ModUser
