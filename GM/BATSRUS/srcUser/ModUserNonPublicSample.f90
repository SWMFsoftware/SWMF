!#NOTPUBLIC  email:darrens@umich.edu  expires:12/31/2099
!^CFG COPYRIGHT UM
!==============================================================================
module ModUser

  ! This file contains a sample of a protected ModUser file.
  ! The top line above must follow the format shown.
  ! If a user file is included in a test of the code, it will be included
  !   regardless of the inclusion of the NONPUBLIC line.

  use ModUserEmpty,                                     &
       IMPLEMENTED1 => user_read_inputs

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'NOTPUBLIC Protected Sample'

contains

  !============================================================================
  subroutine user_read_inputs

    use ModReadParam,   ONLY: read_line, read_command, read_var
    character (len=100) :: NameCommand
    character(len=*), parameter :: NameSub = "ModUser::user_read_inputs"
    !--------------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE

       select case(NameCommand)
       case default
          call stop_mpi(NameSub//': unknown command name='//trim(NameCommand))
       end select
    end do
  end subroutine user_read_inputs

  !============================================================================

end module ModUser

