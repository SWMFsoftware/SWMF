!^CFG COPYRIGHT UM
!=============================================================================
! Timing module for general timing
!=============================================================================

subroutine timing_version(on,name,number)

  logical, intent(out)            :: on
  character (len=40), intent(out) :: name
  real, intent(out)               :: number

  on    =.true.
  name  ='TIMING by G. Toth (2001)'
  number=1.2

end subroutine timing_version

!==============================================================================
subroutine timing_active(value)

  use ModTiming, ONLY: UseTiming
  implicit none

  logical, intent(in) :: value
  !----------------------------------------------------------------------------

  UseTiming=value

end subroutine timing_active

!==============================================================================
subroutine timing_step(value)

  use ModTiming, ONLY: step
  implicit none
  integer, intent(in) :: value
  !----------------------------------------------------------------------------

  step = value

end subroutine timing_step

!==============================================================================
subroutine timing_comp_proc(value1,value2)

  use ModTiming, ONLY: NameComp, iProc
  implicit none
  character (len=*), intent(in) :: value1
  integer, intent(in) :: value2
  
  !----------------------------------------------------------------------------

  NameComp = value1
  iProc    = value2

end subroutine timing_comp_proc

!==============================================================================
subroutine timing_depth(value)

  use ModTiming, ONLY: max_depth
  implicit none
  integer, intent(in) :: value
  !----------------------------------------------------------------------------

  max_depth = value

end subroutine timing_depth

!==============================================================================
subroutine timing_report_style(value)

  use ModTiming, ONLY: report_style
  implicit none
  character (LEN=*), intent(in) :: value
  !----------------------------------------------------------------------------

  report_style = value

end subroutine timing_report_style

!==============================================================================
subroutine timing_param_put_i(name,value,error)

  use ModTiming
  implicit none

  character (LEN=*), intent(in) :: name
  integer, intent(in) :: value
  integer, intent(out):: error
  !----------------------------------------------------------------------------

  error=0
  select case(name)
  case('step')
     step=value
  case('depth')
     max_depth=value
  case('verbose')
     lVerbose=value
  case default
     error=-1
  end select

end subroutine timing_param_put_i

!==============================================================================

function timing_func_d(func_name,iclock,name,parent_name)

  use ModTiming
  implicit none

  real(Real8_) :: timing_func_d
  character (LEN=*), intent(in):: func_name, name, parent_name
  integer, intent(in) :: iclock

  integer     :: i, qclock, qiter, qcall
  real(Real8_) :: qsum
  !----------------------------------------------------------------------------

  ! write(*,*)'func_name=',func_name,' name=',name,' parent_name=',parent_name

  if(.not.UseTiming)then
     timing_func_d=-1.0
     RETURN
  end if

  do i=1,ntiming
     if(sa_name(i)==name .and. sa_name(ia_parent(i))==parent_name)EXIT
  end do

  if(i>ntiming)then
     timing_func_d=-1.0
     RETURN
  end if

  qclock = min( max(iclock,1), maxclock)

  if(la_active(i))then
     qsum = da_sum(i,qclock) + timing_cpu() - da_start(i)
  else
     qsum = da_sum(i,qclock)
  end if

  select case(func_name)
  case('sum')
     timing_func_d=qsum
  case('sum/iter')
     qiter=ia_iter(i,qclock)
     if(qiter<1)qiter=-1
     timing_func_d=qsum/qiter
  case('sum/call')
     qcall=ia_call(i,qclock)
     if(qcall<1)qcall=-1
     timing_func_d=qsum/qcall
  case default
     timing_func_d=-1.0
  end select

end function timing_func_d

!==============================================================================

subroutine timing_start(name)

  use ModTiming
  implicit none

  character (LEN=*), intent(in):: name

  integer :: i
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  current_depth = current_depth + 1

  if(max_depth >= 0 .and. current_depth > max_depth) RETURN

  if(lVerbose>2)write(*,*)'timing_start for ',name

  ! Search for previous timings of the same name at the same depth
  do i = i_last+1, ntiming
     if(ia_depth(i)==current_depth .and. &
          sa_name(i)==name) goto 100
  end do
  ! New name
  if(ntiming==maxtiming-1)write(*,*) &
       'WARNING: number of timings has reached maxtiming in ModTiming'

  if(ntiming==maxtiming) RETURN ! Cannot add more timing

  ntiming=ntiming+1
  i=ntiming
  sa_name(i)    = name
  ia_step(i)    = -1

100  continue

  ia_call(i,2:maxclock) = ia_call(i,2:maxclock)+1
  if(ia_step(i) < step) &
       ia_iter(i,2:maxclock) = ia_iter(i,2:maxclock)+1
  ia_step(i)     = step
  ia_depth(i)    = current_depth
  ia_parent(i)   = i_last
  la_active(i)   = .true.
  da_start(i)    = timing_cpu()
  i_last         = i

  if(lVerbose>2)write(*,*)'index, start:',i,da_start(i)
  
end subroutine timing_start

!==============================================================================
subroutine timing_stop(name)

  use ModTiming
  implicit none

  character (LEN=*), intent(in):: name

  integer     :: i
  real(Real8_) :: qnow
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  current_depth = current_depth - 1

  if(max_depth >= 0 .and. current_depth >= max_depth) RETURN

  if(lVerbose>2)write(*,*)'timing_stop  for ',name

  ! Write a warning for non-matching name
  if(sa_name(i_last)/=name.and.ntiming<maxtiming)write(*,*) &
       'WARNING in timing: unexpected STOP requested for ',name

  i=i_last
  qnow   = timing_cpu()
  da_sum(i,1)          = qnow-da_start(i)
  da_sum(i,2:maxclock) = da_sum(i,2:maxclock) + da_sum(i,1)
  la_active(i)         = .false.
  i_last               = ia_parent(i)

  if(lVerbose>2)then
     write(*,*)'index, stop :',i,qnow
     write(*,*)'clocks:',da_sum(i,:)
  end if

end subroutine timing_stop

!==============================================================================
subroutine timing_reset_all
  call timing_reset('#all',2)
end subroutine timing_reset_all

!==============================================================================
subroutine timing_reset(name,nclock)

  ! reset clocks 1 to nclock for name

  use ModTiming
  implicit none

  character (LEN=*), intent(in):: name
  integer, intent(in) :: nclock

  real(Real8_) :: qnow
  integer     :: i, qclock
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  if(lVerbose>1)write(*,*)'timing_reset 1..nclock=',nclock

  qclock = min(max(nclock,1),maxclock)

  qnow = timing_cpu()

  do i=1,ntiming
     if(name/='#all'.and.sa_name(i)/=name)CYCLE

     if(la_active(i))then
        da_sum(i,nclock+1:maxclock) = da_sum(i,nclock+1:maxclock) &
             + qnow-da_start(i)
        da_start(i) = qnow
        ia_step(i)  = step
        ia_iter(i,2:nclock) = 1
        ia_call(i,2:nclock) = 1
        da_sum(i,1:nclock)  = 0.0
     else
        ia_step(i) = -1
        ia_iter(i,1:nclock) = 0
        ia_call(i,1:nclock) = 0
        da_sum(i,1:nclock)  = 0.0
     end if
  end do

  if(name=='#all')step_reset(2:nclock)=step

end subroutine timing_reset

!==============================================================================

subroutine timing_show(name,iclock)

  ! Report timing results by clock iclock for name. 

  use ModTiming
  implicit none

  character (LEN=*), intent(in):: name
  integer,           intent(in):: iclock

  integer     :: i, i_parent, qclock, qiter, qcall
  real(Real8_) :: qnow, qsum, qsumparent
  character (LEN=40) :: s_parent
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  qclock = min(max(1,iclock),maxclock)

  do i=ntiming,1,-1
     if(sa_name(i)/=name)CYCLE

     i_parent=ia_parent(i)
     s_parent=sa_name(i_parent)
     qnow=timing_cpu()
     if(la_active(i))then
        qsum = da_sum(i,qclock) + qnow-da_start(i)
     else
        qsum = da_sum(i,qclock)
     end if

     if(qclock==1)then
        write(*,'(5a,f8.2,a)')                            &
                'Timing for last ',name,                  &
                ' (',s_parent(1:len_trim(s_parent)),'):', &
                qsum,' sec'
        RETURN
     end if

     if(la_active(i_parent))then
        qsumparent = da_sum(i_parent,qclock) + qnow-da_start(i_parent)
     else
        qsumparent = da_sum(i_parent,qclock)
     end if
     if(qsumparent<=0.0)qsumparent=-1.0

     qiter=ia_iter(i,qclock); if(qiter<1)qiter=-1
     qcall=ia_call(i,qclock); if(qcall<1)qcall=-1

     write(*,'(2a)',ADVANCE='NO')'Timing for ',name
     if(qclock==maxclock .or. step_reset(qclock)==-1)then
        write(*,'(a,i8,a)') ' at step',step,' :'
     else
        write(*,'(a,i8,a,i8,a)')&
             ' from step',step_reset(qclock),' to',step,' :'
     end if
     write(*,'(f8.2,a,f8.3,a,f8.3,a,f8.2,2a)')           &
                qsum,' sec, ',                            &
                qsum/qiter,' s/iter',                     &
                qsum/qcall,' s/call',                     &
                100.0*qsum/qsumparent,' % of ',s_parent(1:len_trim(s_parent))
  end do

end subroutine timing_show

!==============================================================================
subroutine timing_report
  use ModTiming, ONLY: report_style

  select case(report_style)
  case('tree')
     call timing_tree(2,-1)
  case('list')
     call timing_sort(2,-1,.false.)
  case default
     call timing_sort(2,-1,.true.)
  end select

end subroutine timing_report

!==============================================================================
subroutine timing_report_total

  use ModTiming, ONLY: report_style

  select case(report_style)
  case('tree')
     call timing_tree(3,-1)
  case('list')
     call timing_sort(3,-1,.false.)
  case default
     call timing_sort(3,-1,.true.)
  end select

end subroutine timing_report_total

!==============================================================================
subroutine timing_tree(iclock,show_depth)

  ! Produce a timing report that indicates calling sequence with order
  !    and nesting of calls with indentation
  ! Use clock iclock
  ! If show_depth is positive, show only the timings which are 
  !    at a nesting level not exceeding show_depth

  use ModTiming
  implicit none

  integer, intent(in):: iclock, show_depth

  ! Indices
  integer :: i, i_parent

  ! Temporary variables
  integer     :: qclock, qdepth, i_depth, qiter, qcall, indent
  real(Real8_) :: qsum, qsumparent, qnow

  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  qclock = min(max(iclock,2),maxclock)

  write(*,'(a79)')'----------------------------------------'// &
       '------------------------------------------'
  write(*,'(a)',ADVANCE='NO')'TIMING TREE'
  if(show_depth>0)&
       write(*,'(a,i2)',ADVANCE='NO')' of depth',show_depth
  if(step_reset(qclock)>=0)then
     write(*,'(a,i8,a,i8)',ADVANCE='NO') &
          ' from step',step_reset(qclock),' to',step
  else
     write(*,'(a,i8)',ADVANCE='NO')' at step',step
  end if
  write(*,'(a,i4)')' '//NameComp//' on PE ',iProc

  write(*,'(a20,a7,a8,a9,a9,a9,a9)') &
       'name'//spaces,'#iter','#calls','sec','s/iter','s/call','percent'
  write(*,'(a79)')'----------------------------------------'// &
       '------------------------------------------'
  qdepth = 0
  qnow   = timing_cpu()
  do i=1,ntiming
     
     if(show_depth>0 .and. ia_depth(i)>show_depth)CYCLE

     if(la_active(i))then
        da_sum(i,:) = da_sum(i,:) + qnow - da_start(i)
        da_start(i) = qnow
     endif
     qsum=da_sum(i,qclock)
     if(qsum.lt.0.0005)CYCLE

     qsumparent=da_sum(ia_parent(i),qclock)
     if(qsumparent.lt.qsum)qsumparent=-1.
     qiter=ia_iter(i,qclock); if(qiter<1)qiter=-1
     qcall=ia_call(i,qclock); if(qcall<1)qcall=-1
     ! Negative indent results in error because repeat(' ',-2) fails
     ! and spaces(1:-2) (although correct F90) also fails due to a
     ! PGF90 compiler bug
     indent=max(0,ia_depth(i)*2-4)
     write(*,'(a20,i7,i8,a,f9.2,f9.3,f9.3,f9.2)')                    &
          repeat(' ',indent)//sa_name(i),                            &
          qiter,                                                     &
          qcall,                                                     &
          repeat(' ',indent),                                        &
          qsum,                                                      &
          qsum/qiter,                                                &
          qsum/qcall,                                                &
          100*qsum/qsumparent

     if(ia_depth(i)==1) write(*,'(a79)')sepline
     
     ! Add up times for this depth and report missing part

     i_depth=ia_depth(i)
     if(i_depth > qdepth)then
        da_sum_other(i_depth) = qsum
     else
        da_sum_other(i_depth) = da_sum_other(i_depth) + qsum
     end if

     i_parent=i
     do qdepth = i_depth, ia_depth(i+1) + 1, -1

        ! Find parent (of parent (of parent...))
        i_parent=ia_parent(i_parent)
        
        ! Unmeasured timing = parent timing - sum of same level timing:
        qsum       = da_sum(i_parent,qclock)-da_sum_other(qdepth)

        ! reset sum for this depth
        da_sum_other(qdepth) = 0.0

        if(qsum<0.001) CYCLE
        qsumparent = da_sum(i_parent,qclock)
        if(qsumparent.lt.qsum)qsumparent=-1.
        qiter      = ia_iter(i_parent,qclock)
        if(qiter<1)qiter=-1

        indent = max(0,qdepth*2-4)
        write(*,'(a,a35,f9.2,f9.3,f18.2)')    &
             repeat(' ',indent),              &
             '#others'//spaces,               &
             qsum,                            &
             qsum/qiter,                      &
             100*qsum/qsumparent

     end do
     qdepth = i_depth
  end do
  write(*,'(a79)')sepline

end subroutine timing_tree

!==============================================================================
subroutine timing_sort(iclock,show_length,unique)

  ! Make a sorted report of timings made with clock iclock
  ! If show_length is positive, only that many items are listed.
  ! If unique is true, all calls to the same name are added up
  
  use ModTiming
  implicit none
  integer, intent(in) :: iclock, show_length
  logical, intent(in) :: unique

  integer :: ia_sort(maxtiming) ! indirect index array for sorting
  integer :: i, j, k, qclock
  real(Real8_)  :: qsum, qsummax, qnow
  character (LEN=40) :: s_name, s_parent

  integer     :: qntiming, showntiming
  real(Real8_) :: da_qsum(maxtiming)
  integer     :: ia_qcall(maxtiming), ia_qiter(maxtiming)
  character (LEN=40) :: sa_qname(maxtiming)
  
  !----------------------------------------------------------------------------
  if(.not.UseTiming)RETURN

  qclock = min(max(iclock,1),maxclock)

  qnow=timing_cpu()
  if(unique)then
     ! Add up timing and number of iterations and calls for identical names
     da_qsum(1:ntiming) =0.0
     ia_qiter(1:ntiming)=0
     ia_qcall(1:ntiming)=0
     qntiming=0
     do i=1,ntiming
        ! Check if the same name occured before or not
        do j=1,qntiming
           if(sa_name(i)==sa_name(j))EXIT
        end do
        sa_qname(j) = sa_name(i)
        if(la_active(i))then
           qsum = da_sum(i,qclock) + qnow - da_start(i)
        else
           qsum = da_sum(i,qclock)
        end if
        da_qsum(j)  = da_qsum(j)  + qsum
        ia_qiter(j) = ia_qiter(j) + ia_iter(i,qclock)
        ia_qcall(j) = ia_qcall(j) + ia_call(i,qclock)
        if (j > qntiming) qntiming = j
     end do
  else
     qntiming = ntiming
     sa_qname(1:ntiming) = sa_name(1:ntiming)
     where(la_active(1:ntiming))
        da_qsum(1:ntiming) = da_sum(1:ntiming,qclock) &
                           + qnow - da_start(1:ntiming)
     elsewhere
        da_qsum(1:ntiming) = da_sum(1:ntiming,qclock)
     end where
     ia_qiter(1:ntiming) = ia_iter(1:ntiming,qclock)
     ia_qcall(1:ntiming) = ia_call(1:ntiming,qclock)
  end if

  forall(i=1:qntiming)ia_sort(i)=i

  do i=1,qntiming-1
     do j=i+1,qntiming
        if(da_qsum(ia_sort(i)) < da_qsum(ia_sort(j)))then
           k=ia_sort(i)
           ia_sort(i)=ia_sort(j)
           ia_sort(j)=k
        end if
     end do
  end do

  qsummax=da_qsum(ia_sort(1))

  if(qsummax<=0.0)then
     write(*,*)'WARNING in timing_sort: Maximum timing is <= 0 !!!'
     RETURN
  end if

  write(*,'(a79)')sepline

  write(*,'(a,i8,a,i2)',ADVANCE='NO')'SORTED TIMING'
  if(show_length>0)&
       write(*,'(a,i3)',ADVANCE='NO')' of length',show_length
  if(qclock>1 .and. step_reset(qclock)>=0)then
     write(*,'(a,i8,a,i8)',ADVANCE='NO') &
          ' from step',step_reset(qclock),' to',step
  else
     write(*,'(a,i8)',ADVANCE='NO')' at step',step
  end if
  write(*,'(a,i4)')' '//NameComp//' on PE ',iProc

  write(*,'(a20)',ADVANCE='NO')                'name'//spaces
  if(.not.unique)write(*,'(a20)',ADVANCE='NO') '(parent)'//spaces
  write(*,'(a10,a10)',ADVANCE='NO')              'sec','percent'
  if(qclock>1)write(*,'(a10,a10)',ADVANCE='NO')  '#iter','#calls'
  write(*,*)
  write(*,'(a79)')sepline

  if(show_length>0)then
     showntiming=min(show_length,qntiming)
  else
     showntiming=qntiming
  endif

  do k=1,showntiming
     i=ia_sort(k)
     s_name  =sa_qname(i)
     qsum    =da_qsum(i)

     if(qsum < 0.001) CYCLE

     write(*,'(a20)',ADVANCE='NO')       s_name

     if(.not.unique)then
        s_parent=sa_name(ia_parent(i))
        write(*,'(a20)',ADVANCE='NO') &
          '('//s_parent(1:len_trim(s_parent))//')'//spaces
     end if

     write(*,'(f10.2,f10.2)',ADVANCE='NO') qsum, 100.0*qsum/qsummax

     if(qclock>1)write(*,'(i10,i10)',ADVANCE='NO') ia_qiter(i), ia_qcall(i)
     write(*,*)
     if(ia_depth(i)==1) write(*,'(a79)')sepline

  end do

  qsum=0
  do k=showntiming+1,qntiming
     i=ia_sort(k)
     qsum=qsum+da_qsum(i)
  end do
  if(qsum>0.0)then
     write(*,'(a20)',ADVANCE='NO')'#others'//spaces
     if(.not.unique)write(*,'(a20)',ADVANCE='NO')' '
     write(*,'(f8.2,f8.2)')qsum, 100.0*qsum/qsummax
  end if

  write(*,'(a79)')'----------------------------------------'// &
       '------------------------------------------'
end subroutine timing_sort


