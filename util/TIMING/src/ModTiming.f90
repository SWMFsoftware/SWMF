!^CFG COPYRIGHT UM
module ModTiming
  implicit none

  integer, parameter :: Real8_ = selected_real_kind(12)

  logical :: UseTiming=.false.

  integer, parameter :: maxtiming=100, maxclock=3

  character (LEN=40), dimension(maxtiming) :: sa_name
  character (LEN=40), parameter :: spaces = repeat(' ',40)
  character (LEN=79), parameter :: sepline = repeat('-',79)
  real(Real8_), dimension(maxtiming):: da_start=0.0, da_last=0.0, da_sum_other
  real(Real8_), dimension(maxtiming,maxclock) :: da_sum=0.0

  integer, dimension(maxtiming) :: ia_step=-1, ia_depth=0, ia_parent=1
  integer, dimension(maxtiming,maxclock)   :: ia_call=0,  ia_iter=0

  logical, dimension(maxtiming) :: la_active=.false.

  integer :: ntiming=0, i_last=1, step_reset(maxclock)=-1, &
             current_depth=0, max_depth=-1

  character (LEN=10) :: report_style='cumm'

  integer :: step=-1

  integer :: lVerbose=-1

  integer :: iProc=0

  character (LEN=2) :: NameComp='  '

  interface
     function timing_cpu()
       real(selected_real_kind(12)) :: timing_cpu
     end function timing_cpu
  end interface

end module ModTiming
