!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_transfer_data

  ! Transfer scalar or array of integers, reals or logicals 
  ! between two components.
  !
  ! Examples of usage:
  !
  ! Transfer a scalar of type integer that is
  ! available from all source processors and needed on all target processors:
  !
  !    if(is_proc(Source_) call get_data_from_source_comp(iData)
  !    call transfer_integer(Source_, Target_, iData, 
  !               UseSourceRootOnly=.false., UseTargetRootOnly=.false.)
  !    if(is_proc(Target_) call put_data_into_target_comp(iData)
  !
  ! Transfer a 2D array of type real that is obtained as the sum of the data 
  ! returned by the source processors. Result is needed on root target only:
  !
  !    allocate(Data_II(iSize,jSize))
  !    if(is_proc(Source_) call get_data_from_source_comp(Data_II)
  !    call transfer_real_array(Source_, Target_, iSize*jSize, Data_II, 
  !            UseSourceSum=.true., UseTargetRootOnly=.true.)
  !    if(is_proc0(Target_) call put_data_into_target_comp(Data_II)
  !    deallocate(Data_II)
  !
  ! The optional UseSourceSum argument indicates that the data obtained
  ! from the source processors should be added up either on to the root
  ! processor (default) or added up on all processors (If it is combined with
  ! UseSourceRootOnly = .false.).
  !
  ! The optional UseSourceRootOnly and UseTargetOnly parameters indicate
  ! if the data is availale/needed only on the Source/Target root processor.
  ! The default is the worst case scenario, i.e.:
  !
  !   UseSourceRootOnly = .true.  (data is only avaialable on the source root) 
  !   UseTargetRootOnly = .false. (data is needed on all target processors)
  !
  ! Setting these logicals to non-default values can reduce the communication.

  use ModMpi
  use CON_comp_param, ONLY: NameComp_I
  use CON_world

  implicit none

  private ! except
  
  public:: transfer_real
  public:: transfer_real_array
  public:: transfer_integer
  public:: transfer_integer_array
  public:: transfer_string
  public:: transfer_string_array

  ! Local variables
  integer:: iStatus_I(MPI_STATUS_SIZE)
  integer:: iError

  logical:: DoSourceSum       ! sum up data on source processors
  logical:: DoRootTransfer    ! transfer data from Source Root to Target Root?
  logical:: DoTargetBroadcast ! broadcast data from Target Root to Target procs

  logical:: UseAllReduce      ! use an MPI_allreduce or MPI_reduce

contains

  !===========================================================================

  subroutine transfer_data_action(DoTest, iCompSource, iCompTarget, &
       UseSourceSum, UseSourceRootOnlyIn, UseTargetRootOnlyIn)

    ! Set DoRootTransfer and DoTargetBroadcast logicals as needed

    logical, intent(in):: DoTest
    integer, intent(in):: iCompSource, iCompTarget
    logical, optional, intent(in):: &
         UseSourceSum, UseSourceRootOnlyIn, UseTargetRootOnlyIn

    integer:: iProc0Source, iProc0Target
    logical:: UseSourceRootOnly, UseTargetRootOnly

    character(len=*), parameter:: NameSub = 'transfer_data_action'
    !-------------------------------------------------------------------------
    ! Do source addition on source processors if needed. Default is false.
    DoSourceSum = .false.
    if(present(UseSourceSum)) DoSourceSum = UseSourceSum &
         .and. is_proc(iCompSource) .and. n_proc(iCompSource) > 0

    ! Default values for using root only represent the worst case:
    UseSourceRootOnly = .true.
    UseTargetRootOnly = .false.
    if(present(UseSourceRootOnlyIn)) UseSourceRootOnly = UseSourceRootOnlyIn
    if(present(UseTargetRootOnlyIn)) UseTargetRootOnly = UseTargetRootOnlyIn

    ! Should we use an MPI_all_reduce on the source processors
    UseAllReduce = DoSourceSum .and. .not. UseSourceRootOnly

    iProc0Source = i_proc0(iCompSource)
    iProc0Target = i_proc0(iCompTarget)

    if(DoTest)then
       write(*,*) NameSub,' starting for source, target=', &
            NameComp_I(iCompSource),', ',NameComp_I(iCompTarget)
       write(*,*) NameSub, ': DoSourceSum, UseAllReduce=', &
            DoSourceSum, UseAllReduce
       write(*,*) NameSub, &
            ': UseSourceRootOnly, TargetRootOnly, iProc0Source, Target=',&
            UseSourceRootOnly, UseTargetRootOnly, iProc0Source, iProc0Target
    end if

    ! transfer is ususally needed if source and target roots are different.
    DoRootTransfer = iProc0Source /= iProc0Target

    if(DoTest)write(*,*) NameSub,': initial DoRootTransfer=', DoRootTransfer

    ! If data is available on all source processors then
    ! check if target root is among them
    if(.not.UseSourceRootOnly .and. iProc0Target > iProc0Source)then
       DoRootTransfer = iProc0Target > i_proc_last(iCompSource) .or. &
            modulo(iProc0Target-iProc0Source, i_proc_stride(iCompSource)) /= 0

       if(DoTest)then
          write(*,*) NameSub,': i_proc_last(Source)=', i_proc_last(iCompSource)
          write(*,*) NameSub,': modulo(DnProc0,StrideSource)=', &
               modulo(iProc0Target-iProc0Source, i_proc_stride(iCompSource))
          write(*,*) NameSub,': final DoRootTransfer=', DoRootTransfer
       end if
    end if

    ! Check if broadcast on target is needed
    DoTargetBroadcast = n_proc(iCompTarget) > 1 .and. .not.UseTargetRootOnly

    if(DoTest)write(*,*) NameSub, &
         ': initial DoTargetBroadcast=', DoTargetBroadcast

    ! If there is no need to broadcast, we are done
    if(.not.DoTargetBroadcast) RETURN

    ! Only target processors may take part in the MPI broadcast
    if(.not.is_proc(iCompTarget))then
       if(DoTest)write(*,*) NameSub, &
            ': not a target proc, set DoTargetBroadcast=F'
       DoTargetBroadcast = .false.
       RETURN
    end if

    ! If only the source root has the data, broadcast is unavoidable
    if(UseSourceRootOnly) RETURN

    ! If the target root was not covered by source, it is unlikely 
    ! that the other target processors are all covered by source
    if(DoRootTransfer) RETURN

    ! Since DoRootTransfer is false, the target root is covered by source.
    ! Check if the last target processor is within the range of source procs,
    ! and the stride of target is an integer multiple of the stride of source.

    DoTargetBroadcast = &
         i_proc_last(iCompTarget) > i_proc_last(iCompSource) .or. &
         modulo(i_proc_stride(iCompTarget), i_proc_stride(iCompSource)) /= 0

    if(DoTest)then
       write(*,*)NameSub,': i_proc_last(Source)=', i_proc_last(iCompSource)
       write(*,*)NameSub,': i_proc_last(Target)=', i_proc_last(iCompTarget)
       write(*,*)NameSub,': modulo(stride(Target),stride(Source))=', &
            modulo(i_proc_stride(iCompTarget), i_proc_stride(iCompSource))
       write(*,*)NameSub,': final DoTargetBroadcast=', DoTargetBroadcast
       
    end if

  end subroutine transfer_data_action

  !===========================================================================
  subroutine transfer_integer_array(iCompSource, iCompTarget, nData, iData_I, &
       UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(in):: nData
    integer, intent(inout):: iData_I(:)

    logical, optional, intent(in):: &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly

    integer, parameter:: iTag = 1001

    integer, allocatable:: iDataLocal_I(:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_integer'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call transfer_data_action(DoTestMe, iCompSource, iCompTarget, &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    if(DoTestMe)then
       write(*,*) NameSub,': DoSourceSum, DoRootTransfer, DoTargetBroadcast=',&
            DoSourceSum, DoRootTransfer, DoTargetBroadcast
       call CON_stop('DEBUG')
    end if

    if(DoSourceSum)then
       allocate(iDataLocal_I(nData))
       iDataLocal_I = iData_I
       if(UseAllReduce)then
          call MPI_allreduce(iDataLocal_I, iData_I, nData, MPI_INTEGER, &
               MPI_SUM, i_comm(iCompSource), iError)
       else
          call MPI_reduce(iDataLocal_I, iData_I, nData, MPI_INTEGER, &
               MPI_SUM, 0, i_comm(iCompSource), iError)
       end if
       deallocate(iDataLocal_I)
    end if

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            iData_I, nData, MPI_INTEGER, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            iData_I, nData, MPI_INTEGER, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         iData_I, nData, MPI_INTEGER, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_integer_array

  !===========================================================================
  subroutine transfer_integer(iCompSource, iCompTarget, iData, &
       UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(inout):: iData

    logical, optional, intent(in):: &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly

    integer, parameter:: iTag = 1002

    integer:: iDataLocal

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_integer'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call transfer_data_action(DoTestMe, iCompSource, iCompTarget, &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    if(DoSourceSum)then
       iDataLocal = iData
       if(UseAllReduce)then
          call MPI_allreduce(iDataLocal, iData, 1, MPI_INTEGER, &
               MPI_SUM, i_comm(iCompSource), iError)
       else
          call MPI_reduce(iDataLocal, iData, 1, MPI_INTEGER, &
               MPI_SUM, 0, i_comm(iCompSource), iError)
       end if
    end if

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            iData, 1, MPI_INTEGER, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            iData, 1, MPI_INTEGER, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         iData, 1, MPI_INTEGER, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_integer

  !===========================================================================

  subroutine transfer_real_array(iCompSource, iCompTarget, nData, Data_I, &
       UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(in):: nData
    real, intent(inout):: Data_I(nData)

    logical, optional, intent(in):: &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly

    real, allocatable:: DataLocal_I(:)

    integer, parameter:: iTag = 1003
    !-------------------------------------------------------------------------
    call transfer_data_action(.false., iCompSource, iCompTarget, &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    if(DoSourceSum)then
       allocate(DataLocal_I(nData))
       DataLocal_I = Data_I
       if(UseAllReduce)then
          call MPI_allreduce(DataLocal_I, Data_I, nData, MPI_REAL, MPI_SUM, &
               i_comm(iCompSource), iError)
       else
          call MPI_reduce(DataLocal_I, Data_I, nData, MPI_REAL, MPI_SUM, &
               0, i_comm(iCompSource), iError)
       end if
       deallocate(DataLocal_I)
    end if

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            Data_I, nData, MPI_REAL, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            Data_I, nData, MPI_REAL, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         Data_I, nData, MPI_REAL, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_real_array

  !===========================================================================

  subroutine transfer_real(iCompSource, iCompTarget, Data, &
       UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    real, intent(inout):: Data

    logical, optional, intent(in):: &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly

    real:: DataLocal

    integer, parameter:: iTag = 1004
    !-------------------------------------------------------------------------
    call transfer_data_action(.false., iCompSource, iCompTarget, &
         UseSourceSum, UseSourceRootOnly, UseTargetRootOnly)

    if(DoSourceSum)then
       DataLocal = Data
       if(UseAllReduce)then
          call MPI_allreduce(DataLocal, Data, 1, MPI_REAL, MPI_SUM, &
               i_comm(iCompSource), iError)
       else
          call MPI_reduce(DataLocal, Data, 1, MPI_REAL, MPI_SUM, 0, &
               i_comm(iCompSource), iError)
       end if
    end if

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            Data, 1, MPI_REAL, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            Data, 1, MPI_REAL, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         Data, 1, MPI_REAL, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_real

  !===========================================================================

  subroutine transfer_string_array(&
       iCompSource, iCompTarget, nString, String_I,&
       UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(in):: nString
    character(len=*), intent(inout):: String_I(nString)

    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    integer:: nData

    integer, parameter:: iTag = 1005
    !-------------------------------------------------------------------------
    call transfer_data_action(.false., iCompSource, iCompTarget, &
         UseSourceRootOnly, UseTargetRootOnly)

    nData = nString*len(String_I(1))

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            String_I, nData, MPI_CHARACTER, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            String_I, nData, MPI_CHARACTER, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         String_I, nData, MPI_CHARACTER, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_string_array

  !===========================================================================

  subroutine transfer_string(iCompSource, iCompTarget, String,&
       UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    character(len=*), intent(inout):: String

    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    integer:: nData

    integer, parameter:: iTag = 1006
    !-------------------------------------------------------------------------
    call transfer_data_action(.false., iCompSource, iCompTarget, &
         UseSourceRootOnly, UseTargetRootOnly)

    nData = len(String)

    if(DoRootTransfer)then
       if(is_proc0(iCompSource)) call MPI_send( &
            String, nData, MPI_CHARACTER, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            String, nData, MPI_CHARACTER, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(DoTargetBroadcast) call MPI_bcast( &
         String, nData, MPI_CHARACTER, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_string

  !===========================================================================

end module CON_transfer_data
