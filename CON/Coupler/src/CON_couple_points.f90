!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!BOP
!MODULE: CON_couple_pointwise - couple two components in a pointwise manner
!
!DESCRIPTION:
! Couple two components in a pointwise manner
!
!INTERFACE:
module CON_couple_points

  use CON_coupler

  !PUBLIC MEMBER FUNCTIONS:
  public :: couple_points_init
  public :: couple_points
  public :: couple_points_finalize

  type CouplePointsType

     ! Array of variable names
     character(len=lNameVar) :: NameVar

     ! id's of source and target models
     integer :: iCompSource, iCompTarget

     ! Number of dimensions and number of variables
     integer :: nDim=0, nVar=0, nData=0

     ! Router communicator info
     integer :: nProcSourceTarget=0, iCommSourceTarget=-1, iProcSourceTarget=-1, iProcWorld=-1

     ! number of processors of the OTHER component to communicate with
     integer :: nCoupleSource, nCoupleTarget

     ! processors of the OTHER component to communicate with
     integer, allocatable :: iCoupleProcSource_I(:), iCoupleProcTarget_I(:)

     ! Check if the Source and Target processors coincide
     logical:: IsSameLayout

     integer :: iDecompTarget=-1, iGridSource=-1, iGridTarget=-1

     ! number of entries received/sent by a processor during rendezvous
     integer, allocatable :: nCouplePointSource_I(:), nCouplePointTarget_I(:)

     ! Point positions local on a GM processor
     real, allocatable :: PosSource_DI(:,:)

     ! Original PT point positions
     real, pointer :: PosTarget_DI(:,:)

     integer, allocatable :: nPointSource_P(:), nPointTarget_I(:), nPointTarget_P(:), iPointTarget_I(:)
     integer :: nPointSource, nPointTarget
     real, allocatable :: DataSource_VI(:,:), DataTarget_VI(:,:)
  end type CouplePointsType

contains

  subroutine couple_points_init(Coupler)

    type(CouplePointsType), intent(inout) :: coupler
    integer :: iProc0Source, iProc0Target
    logical       :: UseMe = .true.
    integer:: iError, iStatus_I(MPI_STATUS_SIZE)
    character(len=*), parameter:: NameSub = 'couple_points_init'
    
    write(*,*)NameSub,': Initializing coupler'

    Coupler%iProcWorld = i_proc()

    ! Get the union communicator, the UseMe logical, 
    ! and the root proc indexes of the two components in the union comm.
    iProc0Source = i_proc0(Coupler%iCompSource)
    iProc0Target = i_proc0(Coupler%iCompTarget)
    call set_router_comm(Coupler%iCompSource, Coupler%iCompTarget, Coupler%iCommSourceTarget, UseMe, iProc0Source, iProc0Target)

    call MPI_comm_size(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, iError)
    call MPI_comm_rank(Coupler%iCommSourceTarget, Coupler%iProcSourceTarget, iError)

    ! create group on intersection of Coupler%iCompSource and Coupler%iCompTarget
    nProcCommon = n_proc(Coupler%iCompSource) + n_proc(Coupler%iCompTarget) - Coupler%nProcSourceTarget

    Coupler%IsSameLayout = n_proc(Coupler%iCompSource)      == n_proc(Coupler%iCompTarget) &
         .and.     i_proc0(Coupler%iCompSource)     == i_proc0(Coupler%iCompTarget) &
         .and.     i_proc_last(Coupler%iCompSource) == i_proc_last(Coupler%iCompTarget)

    !! Check if nDim in GM is the same as nDimPt
    !! nDim info could be stored in Grid_C!
    !if(is_proc(Coupler%iCompSource)) call GM_get_grid_info(nDim, iGridGm, iDecompGm)
    !if(.not.IsSameLayout) call MPI_bcast(&
    !     nDim, 1, MPI_INTEGER, iProc0Pt, Coupler%iCommSourceTarget, iError)
    !if(is_proc0(Coupler%iCompTarget))then
    !   call PT_get_grid_info(nDimPt, iGridPt, iDecompPt)
    !   if(nDim /= nDimPt)then
    !      write(*,*) NameSub,' ERROR: nDim, nDimPt=', nDim, nDimPt
    !      call CON_stop(NameSub// &
    !          ': GM and PT have different number of spatial dimensions')
    !   endif
    !end if

    write(*,*)NameSub,': Done initializing coupler'

  end subroutine couple_points_init

  subroutine couple_points(Coupler,Source_get, Target_Put, Source_get_grid_info, Target_get_grid_info)

     type(CouplePointsType), intent(inout) :: Coupler

    interface 
        subroutine Target_put(NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

          implicit none
          ! List of variables
          character(len=*), intent(inout):: NameVar

          ! Number of variables in Data_VI
          integer,          intent(inout):: nVar    

          ! Number of points in Pos_DI
          integer,          intent(inout):: nPoint  

          ! Position vectors
          real, pointer:: Pos_DI(:,:)               

          ! Recv data array
          real,    intent(in), optional:: Data_VI(:,:)

          ! Order of data
          integer, intent(in), optional:: iPoint_I(nPoint)

        end subroutine Target_put

        subroutine Source_get(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
             Data_VI)
          logical,          intent(in):: IsNew   ! true for new point array
          character(len=*), intent(in):: NameVar ! List of variables
          integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
          integer,          intent(in):: nDimIn  ! Dimensionality of positions
          integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

          real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
          real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array
        end subroutine Source_get

        subroutine Source_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine Source_get_grid_info

        subroutine Target_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine Target_get_grid_info

     end interface

    ! Is there a need to recalculate the data transfer route?
    logical:: IsNewRoute
    character(len=*), parameter:: NameSub = 'couple_points'

    ! Get new grid info
    ! Check parameters for changes
    ! If needed:
    !   (re)allocate router arrays
    !   (re)do mapping
    !   (re)do routing
    ! Get data from source model
    ! Send data to receiver model

    ! Get the data from source

    write(*,*) NameSub,': Allocating DatSource_VI'
    allocate(Coupler%DataSource_VI(Coupler%nVar,Coupler%nPointSource))
    if(is_proc(Coupler%iCompSource)) call Source_get( IsNewRoute, &
         Coupler%NameVar, Coupler%nVar, Coupler%nDim, Coupler%nPointSource, Coupler%PosSource_DI, Coupler%DataSource_VI)

    ! Transfer data to the correct processor
    allocate(Coupler%DataTarget_VI(Coupler%nVar, Coupler%nData))
    write(*,*) NameSub,': nVar=',Coupler%nVar
    if(Coupler%IsSameLayout)then
       write(*,*) NameSub, ': Transferring buffer (same layout)'
       call transfer_buffer(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, Coupler%nVar, &
            Coupler%nPointSource, Coupler%nPointSource_P, Coupler%DataSource_VI, &
            Coupler%nData,    Coupler%nPointTarget_P, Coupler%DataTarget_VI) 
    else
       write(*,*) NameSub, ': Transferring buffer (different layout)'
       call transfer_buffer_real(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, &
            Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I,           &
            Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I,           &
            Coupler%nVar, Coupler%nPointSource, Coupler%DataSource_VI, Coupler%nData, Coupler%DataTarget_VI)
    end if
    deallocate(Coupler%DataSource_VI)

    write(*,*) NameSub,': Passing data to target'
    ! Give the data to target
    if(is_proc(Coupler%iCompTarget)) call Target_put( &
         Coupler%NameVar, Coupler%nVar, Coupler%nPointTarget, Coupler%PosTarget_DI, Coupler%DataTarget_VI, Coupler%iPointTarget_I)

    deallocate(Coupler%DataTarget_VI)

  end subroutine couple_points

  subroutine couple_points_finalize(Coupler)

    type(CouplePointsType), intent(inout) :: coupler

    ! Deallocate arrays

  end subroutine couple_points_finalize

  !===========================================================================

  subroutine transfer_buffer(iComm, nProc, iProc, nData, &
       nBufferS, nBufferS_P, BufferS_I, &
       nBufferR, nBufferR_P, BufferR_I)

    use ModMpi

    integer, intent(in):: iComm                    ! MPI communicator
    integer, intent(in):: nProc                    ! number of processors
    integer, intent(in):: iProc                    ! local proc index
    integer, intent(in):: nData                    ! number of reals per point
    integer:: nBufferS, nBufferR                   ! total buffer sizes

    integer, intent(in):: nBufferS_P(0:nProc-1)        ! send buffer chunks 
    integer, intent(in):: nBufferR_P(0:nProc-1)        ! recv buffer chunks

    real, intent(in) :: BufferS_I(nData*nBufferS)  ! send buffer
    real, intent(out):: BufferR_I(nData*nBufferR)  ! recv buffer

    ! index of recv and send processors
    integer:: iProcR, iProcS

    integer:: iBufferS, iBufferR

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iRequestR, iRequestS, iError
    integer, allocatable:: iRequestS_I(:), iRequestR_I(:), iStatus_II(:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData

    allocate(iRequestS_I(nProc), iRequestR_I(nProc), &
         iStatus_II(MPI_STATUS_SIZE,nProc))

    ! Possibly optimize for local copies?

    ! Send the appropriate parts of the send buffer
    iRequestS = 0
    iBufferS  = 1
    do iProcR = 0, nProc - 1
       if(nBufferS_P(iProcR) == 0) CYCLE        ! Skip empty chunks

       iRequestS = iRequestS + 1                ! Count send requests

       ! Post send
       call MPI_isend(BufferS_I(iBufferS), nData*nBufferS_P(iProcR), &
            MPI_REAL, iProcR, iTag, iComm, iRequestS_I(iRequestS), iError)

       ! Jump to the starting point of the next chunk
       iBufferS  = iBufferS + nData*nBufferS_P(iProcR) 
    end do

    ! Recv the appropriate parts from the send processors
    iRequestR = 0
    iBufferR  = 1
    do iProcS = 0, nProc - 1
       if(nBufferR_P(iProcS) == 0) CYCLE            ! Skip empty chunks

       iRequestR = iRequestR + 1                    ! Count recv requests

       ! Post recv
       call MPI_irecv(BufferR_I(iBufferR), nData*nBufferR_P(iProcS), &
            MPI_REAL, iProcS, iTag, iComm, iRequestR_I(iRequestR), iError)

       ! Jump to the starting point of the next chunk
       iBufferR  = iBufferR + nData*nBufferR_P(iProcS)  

    end do

    ! wait for all receives to be completed
    if(iRequestR > 0) &
         call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

    ! wait for all sends to be completed
    if(iRequestS > 0) &
         call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)

    deallocate(iRequestS_I, iRequestR_I, iStatus_II)

  end subroutine transfer_buffer

  !===========================================================================

  subroutine transfer_buffer_real(&
       iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I,&
       nData, &
       nBufferS, BufferS_I, nBufferR, BufferR_I)

    ! This subroutine transfers array of reals
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present 

    use ModMpi


    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index

    integer, intent(in):: nCoupleS                 ! #   of procs to send to
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points  to send

    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iCoupleProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nCouplePointR_I(nCoupleR)! # of points  to recv

    integer, intent(in)   :: nData                 ! number of items per point
    integer, intent(in)   :: nBufferS                  ! send buffer size
    real,    intent(inout):: BufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                  ! recv buffer size
    real,    intent(inout):: BufferR_I(nData*nBufferR) ! recv buffer

    integer:: iProcR, iProcS ! R/S is proc to Recv from/Send to
    integer:: iCouple, iBuffer

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:), iStatus_II(:,:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer_direct'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData

    if(nCoupleS + nCoupleR == 0) RETURN

    allocate(iRequest_I(nCoupleS+nCoupleR))
    allocate(iStatus_II(MPI_STATUS_SIZE,nCoupleS+nCoupleR))

    !\
    ! Send buffer
    !/
    if(nCoupleS > 0)then
       !\
       ! Send array of reals
       !/
       iBuffer = 1
       do iCouple = 1, nCoupleS
          iProcS = iCoupleProcS_I(iCouple)
          call MPI_isend(&
               BufferS_I(iBuffer), nData*nCouplePointS_I(iCouple), &
               MPI_REAL, iProcS, iTag, iComm, iRequest_I(iCouple),iError)
          iBuffer = iBuffer + nData*nCouplePointS_I(iCouple)
       end do
    end if
    if(nCoupleR > 0)then
       !\
       ! Receive array of reals
       !/
       iBuffer = 1
       do iCouple = 1, nCoupleR
          iProcR = iCoupleProcR_I(iCouple)
          call MPI_irecv(&
               BufferR_I(iBuffer), nData*nCouplePointR_I(iCouple),&
               MPI_REAL, iProcR, iTag, iComm, &
               iRequest_I(nCoupleS+iCouple),iError)
          iBuffer = iBuffer + nData*nCouplePointR_I(iCouple)
       end do
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I, iStatus_II, iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_real

end module CON_couple_points
