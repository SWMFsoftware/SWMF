!
! This module mimics the coupling between real physics components.
! It uses some MPI communication to synchronize the components.
!
module CON_couple_all

  use CON_comp_param
  use CON_world, ONLY: use_comp, is_proc, i_proc, is_proc0, i_proc0, i_comm, &
       get_comp_info
  use ModMpi
  use CON_time, ONLY: nStep
  use CON_wrapper, ONLY: IsNewInputIe, TimeNewInputSp, DtCpu_C
  use ModUtilities, ONLY: sleep

  implicit none

  private   ! except

  public :: couple_all_init ! initialize all couplers
  public :: couple_two_comp     ! couple 2 components based on their IDs

  ! revision history:
  ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code

  character(len=*), parameter :: NameMod='CON_couple_all'

  integer :: iStatus_I(MPI_STATUS_SIZE), iError, iProc0Source

contains
  !============================================================================

  subroutine couple_all_init
    !--------------------------------------------------------------------------
    write(*,*)'couple_all_init was called'
  end subroutine couple_all_init
  !============================================================================

  subroutine couple_two_comp(iCompSource, iCompTarget, TimeSimulation)

    integer,  intent(in) :: iCompSource, iCompTarget ! component IDs
    real,     intent(in) :: TimeSimulation           ! coupling simulation time

    ! Couple two components given with their IDs. The simulation time
    ! is shared at the time of coupling. Call the appropriate coupling.
    ! Stop with an error for invalid component pairs.

    ! revision history:
    ! 27Aug03 - G. Toth <gtoth@umich.edu> initial prototype/prolog/code

    integer :: iUnitOut
    logical :: DoTest,DoTestMe

    character(len=*), parameter:: NameSub = 'couple_two_comp'
    !--------------------------------------------------------------------------
    call check_i_comp(iCompSource,NameSub//': source')
    call check_i_comp(iCompTarget,NameSub//': target')

    ! Return if any component is not used or if the PE is
    ! used by neither components.
    if(.not.(use_comp(iCompSource))) RETURN
    if(.not.(use_comp(iCompTarget))) RETURN
    if(.not.(is_proc(iCompSource).or.is_proc(iCompTarget))) RETURN

    write(*,*)NameSub,': coupling iProc=',i_proc(),' ',&
         NameComp_I(iCompSource),' --> ',NameComp_I(iCompTarget),&
         ' at nStep,tSim=',nStep,TimeSimulation

    iProc0Source = i_proc()

    ! SP is special as it can only solve up to the last coupling time
    if(iCompTarget == SP_)TimeNewInputSp = TimeSimulation

    ! IE is special as it does not have time.
    ! It solves when data is requested and new info is given.
    if(iCompTarget == IE_)IsNewInputIe = .true.

    if(iCompSource == IE_ .and. is_proc(IE_) .and. IsNewInputIe)then
       call get_comp_info(IE_,iUnitOut=iUnitOut)
       write(iUnitOut,*)'IE:',NameSub,' solve for ',NameComp_I(iCompTarget),&
            ' iProc,nStep,tSim=',i_proc(),nStep,TimeSimulation
       call sleep(DtCpu_C(IE_))
       IsNewInputIe = .false.
    end if

    ! if(is_proc0(iCompSource))write(*,*)NameSub,' sending iProc0Source=',&
    !     iProc0Source

    if(is_proc(iCompSource))call MPI_bcast(iProc0Source,1,MPI_INTEGER,&
         0,i_comm(iCompSource),iError)

    ! write(*,*)'bcast source done',i_proc()

    if(is_proc0(iCompSource)) call MPI_send(iProc0Source,1,MPI_INTEGER,&
         i_proc0(iCompTarget),1,i_comm(),iError)

    ! write(*,*)'send source done',i_proc()

    if(is_proc0(iCompTarget)) call MPI_recv(iProc0Source,1,MPI_INTEGER,&
         i_proc0(iCompSource),1,i_comm(),iStatus_I,iError)

    ! write(*,*)'recv target done',i_proc()

    if(is_proc0(iCompTarget)) call MPI_bcast(iProc0Source,1,MPI_INTEGER,&
         0,i_comm(iCompTarget),iError)

    ! write(*,*)'bcast target done',i_proc()

    ! if(is_proc0(iCompTarget))write(*,*)NameSub,' received iProc0Source=',&
    !     iProc0Source,' on iProc0Target=',i_proc()

  end subroutine couple_two_comp
  !============================================================================

end module CON_couple_all
!==============================================================================
