module ModBufferQueue
  
  ! Module contains implementation of a new data structure: buffer queue.
  ! Buffers in a queue have the same size, each contains a pointer to the next.

  ! APPLICATION: storage for data which size is not known in advance,
  ! the structure manages memory more efficiently and faster,
  ! new storage is allocated in chunks when and if the previous are exhausted.

  ! IMPORTANT NOTE: in the current implementation queue may only grow,
  ! buffers are not deallocated when queue is reset

  implicit none

  save

  private ! except

  public:: TypeBufferQueue
  public:: init_queue, reset_queue, add_to_queue, reset_peek_queue, peek_queue

  ! linked buffers: 
  ! each buffer has a predefined length;
  ! once a buffer is fully filled, the next one is allocated
  type TypeBufferQueue
     !\
     ! Data buffer itself
     !/
     integer, allocatable:: Buffer_I(:)
     !\
     ! Length of a single data record in Buffer_I
     !/
     integer:: nRecordLength
     !\
     ! Max number of data records in Buffer_I
     !/
     integer:: nRecordMax
     !\
     ! Current number of records in Buffer_I
     !/
     integer:: nRecord
     !\
     ! Current number of records in the WHOLE queue,
     ! used only at head of a queue, -1 everywhere else
     !/
     integer:: nRecordAll
     !\
     ! Pointer to the next buffer
     !/
     type(TypeBufferQueue), pointer:: PtrNext
     !\
     ! Pointer to the currently being filled buffer in the queue;
     ! used only at head of a queue, NULL everywhere else
     !/
     type(TypeBufferQueue), pointer:: PtrCurr
     !\
     ! Peeking buffer in the queue;
     ! used only at head of a queue, NULL everywhere else
     !/
     type(TypeBufferQueue), pointer:: PtrPeek
     !\
     ! Peeking location in the peeking buffer;
     ! used only at head of a queue, -1 everywhere else
     !/
     integer:: iRecordPeek
  end type TypeBufferQueue

contains

  function is_initialized(Head) Result(IsInit)
    ! detemine whether a queue has been initialized
    type(TypeBufferQueue), intent(in) :: Head
    logical:: IsInit
    !--------------------------------------------
    IsInit = allocated(Head % Buffer_I)
  end function is_initialized

  !===========================================================================

  function is_head(Buffer) Result(IsHead)
    ! detemine whether the provided buffer is a head of a queue
    type(TypeBufferQueue), intent(in) :: Buffer
    logical:: IsHead
    !--------------------------------------------
    IsHead = &
         associated(Buffer % PtrCurr) .and. (Buffer % nRecordAll  >= 0) .and. &
         associated(Buffer % PtrPeek) .and. (Buffer % iRecordPeek >= 0)
  end function is_head

  !===========================================================================

  function is_last(Buffer) Result(IsLast)
    ! detemine whether the provided buffer is the last in its queue
    type(TypeBufferQueue), intent(in) :: Buffer
    logical:: IsLast
    !--------------------------------------------
    IsLast = .not.associated(Buffer % PtrNext)
  end function is_last

  !===========================================================================

  subroutine init_queue(Head, nRecordMax, nRecordLength)
    ! set the characteristic sizes of a buffer queue
    type(TypeBufferQueue), target, intent(inout):: Head
    integer, intent(in)   :: nRecordMax
    integer, intent(in)   :: nRecordLength
    
    character(len=*),parameter:: NameSub='SP:ModBufferQueue:init_queue'
    !------------------------------------------------------------
    ! check whether the queue has aready been initialized
    if(is_initialized(Head)) call CON_stop(NameSub//&
         ': the buffer queue has already been initialized')
    ! init the queue
    Head % nRecordLength = nRecordLength
    Head % nRecordMax    = nRecordMax
    Head % nRecord       = 0
    Head % nRecordAll    = 0
    Head % iRecordPeek   = 1
    allocate(Head % Buffer_I(Head%nRecordMax * Head%nRecordLength))
    nullify( Head % PtrNext)
    Head % PtrPeek => Head
    Head % PtrCurr => Head
  end subroutine init_queue

  !===========================================================================

  subroutine init_next(Buffer)
    ! allocate the next buffer in the queue
    type(TypeBufferQueue), intent(inout):: Buffer
    character(len=*), parameter:: NameSub = 'SP:ModBufferQueue:init_next'
    !-------------------------------------------------------------
    ! check correctness
    if(.not.is_last(Buffer))&
         call CON_stop(NameSub//': trying to extend queue from the middle')
    ! allocate the pointer to the next buffer
    allocate(Buffer % PtrNext)
    allocate(Buffer % PtrNext % &
         Buffer_I(Buffer%nRecordMax * Buffer%nRecordLength))
    ! initialize its members
    Buffer % PtrNext % nRecordLength = Buffer % nRecordLength
    Buffer % PtrNext % nRecordMax    = Buffer % nRecordMax
    Buffer % PtrNext % nRecord       = 0
    Buffer % PtrNext % nRecordAll    =-1
    Buffer % PtrNext % iRecordPeek   =-1
    nullify(Buffer % PtrNext % PtrNext)
    nullify(Buffer % PtrNext % PtrCurr)
    nullify(Buffer % PtrNext % PtrPeek)
  end subroutine init_next

  !============================================================================

  subroutine add_to_queue(Head, iRecord_I)
    ! add new record to the buffer queue
    type(TypeBufferQueue), intent(inout):: Head
    integer,               intent(in)   :: iRecord_I(Head%nRecordLength)
    ! position of new record in the buffer, added for readability
    integer:: iRecordBegin, iRecordEnd
    character(len=*), parameter:: NameSub = 'SP:ModBufferQueue:add_to_queue'
    !--------------------------------------------------------------------------
    ! check correctness
    if(.not. is_head(Head)) call CON_stop(NameSub // &
         ': must provide the head of a queue as input')
    ! check if the last buffer in the queue is full
    if(Head % PtrCurr % nRecord < Head % nRecordMax)then
       ! the current buffer in the queue is NOT full, add new record to it;
       ! position of this record in the buffer
       iRecordBegin = 1 + Head % PtrCurr % nRecord * Head % nRecordLength 
       iRecordEnd   =(1 + Head % PtrCurr % nRecord)* Head % nRecordLength 
    else
       ! the current buffer is full
       ! if next buffer doesn't exist yet, then create a new one
       if(is_last(Head % PtrCurr))  call init_next(Head % PtrCurr)
       ! update pointer to the current buffer in the queue
       Head % PtrCurr => Head % PtrCurr % PtrNext
       ! position of this record in the buffer
       iRecordBegin = 1
       iRecordEnd   = Head % nRecordLength 
    end if
    ! put the record into buffer
    Head % PtrCurr % Buffer_I(iRecordBegin : iRecordEnd) = &
         iRecord_I(1 : Head % nRecordLength)
    ! update record counter
    Head % PtrCurr % nRecord = Head % PtrCurr % nRecord + 1
    Head % nRecordAll        = Head % nRecordAll + 1
  end subroutine add_to_queue

  !============================================================================
  
  subroutine reset_queue(Head)
    ! reset the whole queue
    type(TypeBufferQueue), target, intent(inout):: Head
    type(TypeBufferQueue), pointer:: Buffer
    !--------------------------------------------------------------------------
    Head % nRecordAll  = 0
    Head % iRecordPeek = 1
    Head % PtrCurr    =>Head
    Head % PtrPeek    =>Head
    Buffer => Head
    do 
       Buffer % nRecord = 0
       if(.not.associated(Buffer % PtrNext)) EXIT
       Buffer = Buffer % PtrNext
    end do
  end subroutine reset_queue

  !============================================================================

  subroutine reset_peek_queue(Head)
    ! reset peeking procedure in queue: move peeking location to front
    type(TypeBufferQueue), target, intent(inout):: Head

    character(len=*),parameter:: NameSub = 'SP:ModBufferQueue:reset_peek_queue'
    !-------------------------------------------------------------------------
    ! check correctness
    if(.not. is_head(Head)) call CON_stop(NameSub // &
         ': must provide the head of a queue as input')
    ! move peeking location to front
    Head % PtrPeek     => Head
    Head % iRecordPeek = 1
  end subroutine reset_peek_queue

  !============================================================================

  subroutine peek_queue(Head, iRecord_I)
    ! cycle throught records in the queue;
    ! NOTE: 
    type(TypeBufferQueue), intent(inout):: Head
    integer,               intent(out)  :: iRecord_I(Head % nRecordLength)
    
    integer:: iRecordBegin, iRecordEnd

    character(len=*), parameter:: NameSub = 'SP:ModBufferQueue:peek_queue'
    !------------------------------------------------------
    ! check correctness
    if(.not. is_head(Head)) call CON_stop(NameSub // &
         ': must provide the head of a queue as input')
    if(.not.associated(Head % PtrPeek)) call CON_stop(NameSub//&
         ': reached the end of queue, cannot peek further')

    ! peek at the current record
    iRecordBegin=(Head % iRecordPeek - 1) * Head % nRecordLength + 1
    iRecordEnd  = Head % iRecordPeek      * Head % nRecordLength
    iRecord_I   = Head % PtrPeek % Buffer_I(iRecordBegin : iRecordEnd)
  
    ! check if the end of peeking buffer has been reached
    if(Head % iRecordPeek == Head % nRecordMax)then
       ! move to next buffer
       Head % PtrPeek     => Head % PtrPeek % PtrNext
       Head % iRecordPeek =  1
    else
       ! move peek location within the same buffer
       Head % iRecordPeek = Head % iRecordPeek + 1
    end if
  end subroutine peek_queue

  !============================================================================
end module ModBufferQueue

