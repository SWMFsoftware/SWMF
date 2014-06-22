module CON_transfer_data

  ! Transfer scalar or array of integers, reals or logicals 
  ! between two components.
  !
  ! Usage:
  !
  ! Transfer data that is available (and the same) on all source processors
  !
  !    if(is_proc(Source_) call get_data_from_source_comp(Data)
  !    call transfer_datatype(Source_, Target_, Data, 
  !                       UseSourceRootOnly=..., UseTargetRootOnly)
  !    if(is_proc(Target_) call put_data_into_target_comp(Data)
  !
  ! The data may be of type integer, real, or logical; scalar or array.
  !
  ! The optional UseSourceRootOnly and UseTargetOnly parameters indicate
  ! if the Data is availale/needed only on the Source/Target root processor.
  ! The default is the worst case scenario, i.e.:
  !
  !   UseSourceRootOnly = .true.  (data is only avaialable on the source root) 
  !   UseTargetRootOnly = .false. (data is needed on all target processors)
  !
  ! Setting these logicals to different vacan reduce the required communication

  use ModMpi
  use CON_world

  implicit none

  private ! except
  
  public:: transfer_integer
  public:: transfer_real
  public:: transfer_integer_array
  public:: transfer_real_array

  ! Local variables
  integer:: iStatus_I(MPI_STATUS_SIZE)
  integer:: iError

contains

  !===========================================================================
  subroutine transfer_integer_array(&
       iCompSource, iCompTarget, iData_I, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(inout):: iData_I(:)

    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    integer:: nData
    logical:: DoBcast

    integer, parameter:: iTag = 1001
    !-------------------------------------------------------------------------
    nData = size(iData_I)

    ! Default is UseSourceRootOnly = .true., so we check source root

    ! If the root processors are different send the data to target
    if(i_proc0(iCompSource) /= i_proc0(iCompTarget))then
       if(is_proc0(iCompSource)) call MPI_send( &
            iData_I, nData, MPI_INTEGER, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            iData_I, nData, MPI_INTEGER, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(n_proc(iCompTarget) == 1) RETURN

    ! Default is UseTargetRootOnly = .false., so we do a broadcast
    if(present(UseTargetRootOnly))then
       if(UseTargetRootOnly) RETURN
    end if

    if(is_proc(iCompTarget)) call MPI_bcast( &
         iData_I, nData, MPI_INTEGER, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_integer_array


  !===========================================================================

  subroutine transfer_real_array(&
       iCompSource, iCompTarget, Data_I, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    real, intent(inout):: Data_I(:)

    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    integer:: nData
    logical:: DoBcast

    integer, parameter:: iTag = 1002
    !-------------------------------------------------------------------------
    nData = size(Data_I)

    ! Default is UseSourceRootOnly = .true., so we check source root

    ! If the root processors are different send the data to target
    if(i_proc0(iCompSource) /= i_proc0(iCompTarget))then
       if(is_proc0(iCompSource)) call MPI_send( &
            Data_I, nData, MPI_REAL, i_proc0(iCompTarget),&
            iTag, i_comm(), iError)
       if(is_proc0(iCompTarget)) call MPI_recv( &
            Data_I, nData, MPI_REAL, i_proc0(iCompSource),&
            iTag, i_comm(), iStatus_I, iError)
    end if

    if(n_proc(iCompTarget) == 1) RETURN

    ! Default is UseTargetRootOnly = .false., so we do a broadcast
    if(present(UseTargetRootOnly))then
       if(UseTargetRootOnly) RETURN
    end if

    if(is_proc(iCompTarget)) call MPI_bcast( &
         Data_I, nData, MPI_REAL, 0, i_comm(iCompTarget), iError)

  end subroutine transfer_real_array

  !===========================================================================

  subroutine transfer_integer( &
       iCompSource, iCompTarget, iData, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    integer, intent(inout):: iData
    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    integer:: iData_I(1)
    !------------------------------------------------------------------------
    if(is_proc(iCompSource)) iData_I(1) = iData
    call transfer_integer_array( &
       iCompSource, iCompTarget, iData_I, UseSourceRootOnly, UseTargetRootOnly)
    if(is_proc(iCompTarget)) iData = iData_I(1)

  end subroutine transfer_integer

  !===========================================================================

  subroutine transfer_real( &
       iCompSource, iCompTarget, Data, UseSourceRootOnly, UseTargetRootOnly)

    integer, intent(in):: iCompSource, iCompTarget
    real, intent(inout):: Data
    logical, optional, intent(in):: UseSourceRootOnly, UseTargetRootOnly

    real:: Data_I(1)
    !------------------------------------------------------------------------
    if(is_proc(iCompSource)) Data_I(1) = Data
    call transfer_real_array( &
       iCompSource, iCompTarget, Data_I, UseSourceRootOnly, UseTargetRootOnly)
    if(is_proc(iCompTarget)) Data = Data_I(1)

  end subroutine transfer_real

  !===========================================================================

end module CON_transfer_data
