!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_ModTime
  use ModKind,    ONLY: Real8_ 
  implicit none
  SAVE
  ! Global interation and time
  real         :: SPTime  = 0.0
  integer      :: iIter   = 0
  ! StartTime converted to real
  real(Real8_) :: StartTime
  ! Time of the last data output
  real         :: DataInputTime
contains
  subroutine read_param(NameCommand)
    use ModReadParam, ONLY: read_var
    use ModTimeConvert, ONLY: time_int_to_real
    character(len=*), intent(in):: NameCommand ! From PARAM.in
    ! This is the same default value as in the SWMF
    integer      :: iStartTime_I(7) = (/2000,3,21,10,45,0,0/)  
    character(len=*), parameter :: NameSub='SP:read_param_time'
    !---------------
    select case(NameCommand) 
       case('#NSTEP')
          call read_var('nStep',iIter)
       case('#TIMESIMULATION')
          call read_var('tSimulation',SPTime)
          !The last time input occurred the same time
          DataInputTime = SPTime
       case("#STARTTIME", "#SETREALTIME")
          call read_var('iYear'  ,iStartTime_I(1))
          call read_var('iMonth' ,iStartTime_I(2))
          call read_var('iDay'   ,iStartTime_I(3))
          call read_var('iHour'  ,iStartTime_I(4))
          call read_var('iMinute',iStartTime_I(5))
          call read_var('iSecond',iStartTime_I(6))
          iStartTime_I(7) = 0
          call time_int_to_real(iStartTime_I, StartTime)
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
end module SP_ModTime
