!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module CON_couple_points

  ! Couple two components by sending and receiving positions of
  ! target points and values of the source component at these locations

  use CON_coupler

  implicit none

  public :: couple_points_init
  public :: couple_points

  ! Debugging
  logical, parameter:: DoDebug = .false.

  type CouplePointsType

     ! String of variable names
     character(len=lNameVar) :: NameVar

     ! Indexes of source and target components
     integer :: iCompSource, iCompTarget

     ! Number of dimensions for point coordinates,
     ! number of variables per point, number of data points found
     integer :: nDim=0, nVar=0, nData=0

     ! Router communicator info
     integer :: nProcUnion=-1, iCommUnion=-1
     integer :: iProcUnion=-1, iProcWorld=-1, nProcCommon=-1

     ! True if the Source and Target processors coincide
     logical:: IsSameLayout

     ! Number of source and target procs
     integer :: nProcSource=-1, nProcTarget=-1

     ! Root processor indexes with respect to the union communicator
     integer:: iProc0Source, iProc0Target

     ! Processor indexes with respect to the union communicator
     integer, allocatable:: iProcSource_P(:), iProcTarget_P(:)

     ! number of processors of the OTHER component to communicate with
     integer :: nCoupleSource, nCoupleTarget

     ! processors of the OTHER component to communicate with
     integer, allocatable :: iProcSourceCouple_I(:), iProcTargetCouple_I(:)

     ! number of entries received/sent by a processor during rendezvous
     integer, allocatable :: nPointSourceCouple_I(:), nPointTargetCouple_I(:)

     ! Point positions local on a source processor
     real, allocatable :: PosSource_DI(:,:)

     ! Original target point positions
     real, allocatable :: PosTarget_DI(:,:)

     ! Total number of points sent from this source processor
     integer :: nPointSource

     ! Total number of points received by this target processor
     integer:: nPointTarget

     ! Index of the i-th data point in the sorted data array received by target
     integer, allocatable :: iPointTarget_I(:)

     real, allocatable :: DataSource_VI(:,:), DataTarget_VI(:,:)

     integer:: iDecompLastSource = -1, iDecompLastTarget = -1

  end type CouplePointsType

contains
  !============================================================================
  subroutine couple_points_init(Coupler)

    type(CouplePointsType), intent(inout) :: Coupler

    integer:: iGroup, iProc, iError

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_points_init'
    !--------------------------------------------------------------------------

    if(Coupler%nProcUnion > 0) RETURN ! Already initialized coupler

    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe) &
         write(*,*) NameSub, ' starting for Source, Target=', &
         NameComp_I(Coupler%iCompSource),', ',NameComp_I(Coupler%iCompTarget)

    Coupler%iProcWorld = i_proc()

    ! Number of source and target processors
    Coupler%nProcSource = n_proc(Coupler%iCompSource)
    Coupler%nProcTarget = n_proc(Coupler%iCompTarget)

    ! Get the union communicator and the root proc indexes of the two
    ! components in the union communicator
    Coupler%iProc0Source = i_proc0(Coupler%iCompSource)
    Coupler%iProc0Target = i_proc0(Coupler%iCompTarget)
    call set_router_comm(Coupler%iCompSource, Coupler%iCompTarget, &
         Coupler%iCommUnion, &
         iProc=Coupler%iProc0Source, jProc=Coupler%iProc0Target)

    call MPI_comm_size(Coupler%iCommUnion, Coupler%nProcUnion, iError)
    call MPI_comm_rank(Coupler%iCommUnion, Coupler%iProcUnion, iError)

    ! Number of overlapping processors
    Coupler%nProcCommon = n_proc(Coupler%iCompSource) &
         + n_proc(Coupler%iCompTarget) - Coupler%nProcUnion

    ! The layouts agree if the union and the intersection coincide
    Coupler%IsSameLayout = Coupler%nProcCommon == Coupler%nProcUnion

    allocate( &
         Coupler%iProcSource_P(0:Coupler%nProcSource-1), &
         Coupler%iProcTarget_P(0:Coupler%nProcTarget-1)  )

    if(Coupler%IsSameLayout)then
       ! simple identity
       Coupler%iProcSource_P = [(iProc, iProc=0, Coupler%nProcSource-1)]
       Coupler%iProcTarget_P = [(iProc, iProc=0, Coupler%nProcTarget-1)]
    else
       ! Set up local proc index to union proc index arrays for both
       ! target and source components
       call MPI_Comm_Group(Coupler%iCommUnion, iGroup, iError)

       call MPI_group_translate_ranks( &
            i_group(Coupler%iCompSource), Coupler%nProcSource, &
            [(iProc, iProc=0, Coupler%nProcSource-1)], iGroup, &
            Coupler%iProcSource_P, iError)

       call MPI_group_translate_ranks( &
            i_group(Coupler%iCompTarget), Coupler%nProcTarget, &
            [(iProc, iProc=0, Coupler%nProcTarget-1)], iGroup, &
            Coupler%iProcTarget_P, iError)

       call MPI_group_free(iGroup, iError)
    end if
    if(DoTestMe)then
       write(*,*) NameSub,' nProcUnion, nProcCommon, IsSameLayout=', &
            Coupler%nProcUnion, Coupler%nProcCommon, Coupler%IsSameLayout
       write(*,*) NameSub,' nProcSource,  nProcTarget = ', &
            Coupler%nProcSource, Coupler%nProcTarget
       write(*,*) NameSub,' iProc0Source, iProc0Target = ', &
            Coupler%iProc0Source, Coupler%iProc0Target
       if(.not.Coupler%IsSameLayout)then
          write(*,*) NameSub,' iProcSource_P=', Coupler%iProcSource_P
          write(*,*) NameSub,' iProcTarget_P=', Coupler%iProcTarget_P
       end if
    end if

  end subroutine couple_points_init
  !============================================================================
  subroutine couple_points(Coupler, &
       source_get_grid_info, source_find_points, source_get, &
       target_get_grid_info, target_put)

    ! Transfer data from source to target at the locations defined by target
    ! The algorithm consists of the following steps:
    !
    ! Get new grid info
    ! If first time coupling or the grid or decomposition changed then
    !   (re)allocate router arrays
    !   (re)do mapping
    !   (re)do routing
    ! Get data from source component
    ! Send data to target component

    type(CouplePointsType), intent(inout) :: Coupler

    interface
       subroutine target_put(NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

         ! First call: Provides point positions in Pos_DI
         ! Afterwards: Get data for the points from Data_VI

         implicit none
         character(len=*), intent(inout):: NameVar ! List of variables
         integer,          intent(inout):: nVar    ! Number of vars in Data_VI
         integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
         real,    intent(in), optional:: Data_VI(:,:)     ! Recv data array
         integer, intent(in), optional:: iPoint_I(nPoint) ! Order of data
         real,   intent(out), optional, allocatable:: Pos_DI(:,:) ! Positions

       end subroutine target_put

       subroutine source_get(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
            Data_VI)
         implicit none
         logical,          intent(in):: IsNew   ! true for new point array
         character(len=*), intent(in):: NameVar ! List of variables
         integer,          intent(in):: nVarIn  ! Number of vars in Data_VI
         integer,          intent(in):: nDimIn  ! Dimensionality of positions
         integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

         real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
         real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array
       end subroutine source_get

       subroutine source_get_grid_info(nDim, iGrid, iDecomp)
         implicit none
         integer, intent(out) :: nDim, iGrid, iDecomp
       end subroutine source_get_grid_info

       subroutine target_get_grid_info(nDim, iGrid, iDecomp)
         implicit none
         integer, intent(out) :: nDim, iGrid, iDecomp
       end subroutine target_get_grid_info

       subroutine source_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

         implicit none

         integer, intent(in) :: nDimIn    ! dimension of position vectors
         integer, intent(in) :: nPoint    ! number of positions
         real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
         integer, intent(out):: iProc_I(nPoint)  ! processor owning position

       end subroutine source_find_points

    end interface

    ! Is there a need to recalculate the data transfer route?
    logical:: IsNewRoute

    ! Number of source and target processors
    integer:: nProcSource, nProcTarget

    ! Number of points sent to each target processor
    integer, allocatable :: nPointSource_P(:)

    ! Number of points received from each source processor
    integer, allocatable :: nPointTarget_P(:)

    ! source processor index for target points
    integer, allocatable:: iProcTarget_I(:)
    integer, allocatable:: iProcSource_I(:)

    ! Grid index
    integer :: iDecompSource, iDecompTarget
    integer :: iGridSource,   iGridTarget

    integer:: iError, iPoint

    ! Target positions sorted according to the corresponding source processors
    real, allocatable:: PosSortTarget_DI(:,:)

    ! Debugging
    character(len=4):: NameSourceTarget
    integer:: iUnitTest
    logical:: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'couple_points'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest .or. DoDebug)then
       NameSourceTarget = NameComp_I(Coupler%iCompSource) &
            //            NameComp_I(Coupler%iCompTarget)
       if(DoDebug) iUnitTest = i_proc() + 100
       if(DoTestMe) write(*,*) NameSub,' doing ',NameSourceTarget

       if(DoTest) call timing_start('pnt_cpl_'//NameSourceTarget)
    end if

    nProcSource = Coupler%nProcSource
    nProcTarget = Coupler%nProcTarget

    if(DoTestMe) write(*,*) NameSub, &
         ' nProcSource, nProcTarget=', nProcSource, nProcTarget

    ! Check if there is a need to redo the mapping
    if(is_proc(Coupler%iCompSource)) &
         call source_get_grid_info(Coupler%nDim, iGridSource, iDecompSource)

    if(is_proc(Coupler%iCompTarget)) &
         call target_get_grid_info(Coupler%nDim, iGridTarget, iDecompTarget)

    if(.not.Coupler%IsSameLayout)then

       call MPI_bcast(iDecompSource, 1, MPI_INTEGER, &
            Coupler%iProc0Source, Coupler%iCommUnion, iError)

       call MPI_bcast(iDecompTarget, 1, MPI_INTEGER, &
            Coupler%iProc0Target, Coupler%iCommUnion, iError)

    endif

    IsNewRoute = iDecompSource /= Coupler%iDecompLastSource &
         .or. Coupler%iDecompLastTarget /= iDecompTarget

    if(DoTestMe) write(*,*) NameSub,&
         ' iDecompSource, iDecompTarget, IsNewRoute, IsSameLayout=', &
         iDecompSource, iDecompTarget, IsNewRoute, Coupler%IsSameLayout

    Coupler%iDecompLastSource = iDecompSource
    Coupler%iDecompLastTarget = iDecompTarget

    if(IsNewRoute)then
       if(DoTest)call timing_start('pnt_route_'//NameSourceTarget)

       Coupler%nPointTarget = 0
       Coupler%nPointSource = 0
       Coupler%nData        = 0

       ! Get positions where info is needed from target.
       ! Target will allocate array.
       if(is_proc(Coupler%iCompTarget)) then
          call target_put(Coupler%NameVar, Coupler%nVar, &
               Coupler%nPointTarget, Pos_DI = Coupler%PosTarget_DI)
       else
          ! Array has to be allocated on other processors too
          allocate(Coupler%PosTarget_DI(Coupler%nDim,0))
       end if

       allocate(nPointSource_P(nProcTarget), nPointTarget_P(nProcSource))

       if(allocated(Coupler%iPointTarget_I)) deallocate(Coupler%iPointTarget_I)

       if(Coupler%IsSameLayout) then

          ! Find processors that own the target positions in source
          allocate(iProcTarget_I(Coupler%nPointTarget))

          call source_find_points(Coupler%nDim, Coupler%nPointTarget, &
               Coupler%PosTarget_DI, iProcTarget_I)

          ! Order points according to the owner processor indexes
          allocate(Coupler%iPointTarget_I(Coupler%nPointTarget))

          call get_buffer_order(Coupler%nProcUnion, &
               Coupler%nPointTarget, iProcTarget_I, &
               nPointTarget_P, Coupler%iPointTarget_I, Coupler%nData)

          deallocate(iProcTarget_I)

          ! Rearrange coordinate array according to processor order
          allocate(PosSortTarget_DI(Coupler%nDim,Coupler%nData))
          do iPoint = 1, Coupler%nPointTarget
             if(Coupler%iPointTarget_I(iPoint) < 0)CYCLE
             PosSortTarget_DI(:,Coupler%iPointTarget_I(iPoint)) = &
                  Coupler%PosTarget_DI(:,iPoint)
          end do
          deallocate(Coupler%PosTarget_DI)

          ! Set number of points to be received by the source component
          call set_recv_info(Coupler%iCommUnion, &
               Coupler%nProcUnion, Coupler%iProcUnion, &
               nPointTarget_P, Coupler%nPointSource, &
               nPointSource_P)

          if(DoTestMe) write(*,*) NameSub,' finished router for same layout'

       else
          ! Layouts are different
          if(allocated(Coupler%PosSource_DI)) &
               deallocate(Coupler%PosSource_DI)

          ! Setup communication pattern for finding points on Source
          ! This allocates and sets iProc*Couple_I, nPoint*Couple_I arrays
          call set_inquiry(&
               Coupler%iProcUnion, Coupler%iCompSource, &
               Coupler%iCompTarget, Coupler%nProcCommon, &
               Coupler%nCoupleTarget, Coupler%iProcTargetCouple_I, &
               Coupler%nPointTargetCouple_I, Coupler%nCoupleSource, &
               Coupler%iProcSourceCouple_I, Coupler%nPointSourceCouple_I)

          if(is_proc(Coupler%iCompTarget))then
             ! Number of points to send from target to the selected source
             ! processors
             Coupler%nPointTargetCouple_I(1:Coupler%nCoupleTarget-1) = &
                  Coupler%nPointTarget / Coupler%nCoupleTarget

             ! Last processor gets the rest of points
             Coupler%nPointTargetCouple_I(Coupler%nCoupleTarget) = &
                  Coupler%nPointTarget &
                  - (Coupler%nCoupleTarget - 1)&
                  * (Coupler%nPointTarget/Coupler%nCoupleTarget)
          end if

          ! send number of points from target to source that will be sent to
          ! "source_find_points." result is in Coupler%nPointSourceCouple_I
          ! (number of points to be recieved)
          call get_recv_buffer_size(Coupler%iCommUnion, &
               Coupler%nProcUnion, Coupler%iProcUnion, &
               Coupler%nCoupleTarget, Coupler%iProcTargetCouple_I, &
               Coupler%nPointTargetCouple_I, Coupler%nCoupleSource, &
               Coupler%iProcSourceCouple_I, Coupler%nPointSourceCouple_I)

          ! Allocate buffer for positions on source (zero size on target)
          if(is_proc(Coupler%iCompSource)) &
               Coupler%nPointSource = sum(Coupler%nPointSourceCouple_I)
          allocate(Coupler%PosSource_DI(Coupler%nDim,Coupler%nPointSource))

          ! Send positions from target to source, so source can find the owners
          call transfer_buffer_real( &
               Coupler%iCommUnion, Coupler%nProcUnion, Coupler%iProcUnion, &
               Coupler%nCoupleTarget, Coupler%iProcTargetCouple_I, &
               Coupler%nPointTargetCouple_I, &
               Coupler%nCoupleSource, Coupler%iProcSourceCouple_I, &
               Coupler%nPointSourceCouple_I, &
               Coupler%nDim, Coupler%nPointTarget, Coupler%PosTarget_DI, &
               Coupler%nPointSource, Coupler%PosSource_DI)

          ! Find processors that own the target positions in source
          allocate(iProcSource_I(Coupler%nPointSource))
          if(is_proc(Coupler%iCompSource))&
               call source_find_points(Coupler%nDim, Coupler%nPointSource, &
               Coupler%PosSource_DI, iProcSource_I)
          deallocate(Coupler%PosSource_DI)

          ! send owner processor indexes from source to target
          allocate(iProcTarget_I(Coupler%nPointTarget))
          call transfer_buffer_int( &
               Coupler%iCommUnion, Coupler%nProcUnion, Coupler%iProcUnion, &
               Coupler%nCoupleSource, Coupler%iProcSourceCouple_I, &
               Coupler%nPointSourceCouple_I, Coupler%nCoupleTarget, &
               Coupler%iProcTargetCouple_I, Coupler%nPointTargetCouple_I, &
               1, Coupler%nPointSource, iProcSource_I, &
               Coupler%nPointTarget, iProcTarget_I)

          if(is_proc(Coupler%iCompTarget))then

             allocate(Coupler%iPointTarget_I(Coupler%nPointTarget))

             call get_buffer_order(nProcSource, &
                  Coupler%nPointTarget, iProcTarget_I, &
                  nPointTarget_P, Coupler%iPointTarget_I, &
                  Coupler%nData)

             allocate(PosSortTarget_DI(Coupler%nDim,Coupler%nData))

             ! Rearrange coordinate array according to processor order
             do iPoint = 1, Coupler%nPointTarget
                if(Coupler%iPointTarget_I(iPoint) < 0) CYCLE
                PosSortTarget_DI(:,Coupler%iPointTarget_I(iPoint)) = &
                     Coupler%PosTarget_DI(:,iPoint)
             end do
          else
             allocate(PosSortTarget_DI(Coupler%nDim, Coupler%nData))
          end if

          deallocate(iProcTarget_I, Coupler%PosTarget_DI, iProcSource_I)

          ! Set number of points to be received by the source component
          call set_recv_info2( &
               Coupler%iCommUnion, Coupler%iProcUnion, &
               Coupler%iProc0Target, Coupler%iProcSource_P, &
               Coupler%iCompTarget, Coupler%iCompSource, &
               nPointTarget_P, Coupler%nPointSource, &
               nPointSource_P)

          if(DoTestMe) write(*,*) NameSub,' finished router for diff layout'

       end if ! not IsSameLayout

       ! Compact the routers to the non-zero size couplings
       if(is_proc(Coupler%iCompTarget))then

          call compact_router(nProcSource, &
               Coupler%iProcSource_P, nPointTarget_P, &
               Coupler%nCoupleTarget, &
               Coupler%iProcTargetCouple_I, Coupler%nPointTargetCouple_I)

          if(DoDebug)then
             write(iUnitTest,*) NameSourceTarget,&
                  ' nPointTarget_P=', nPointTarget_P
             write(iUnitTest,*) NameSourceTarget, &
                  ' new nCoupleTarget=', Coupler%nCoupleTarget
             write(iUnitTest,*) NameSourceTarget, &
                  ' new iProcTargetCouple_I =', Coupler%iProcTargetCouple_I
             write(iUnitTest,*) NameSourceTarget, &
                  ' new nPointTargetCouple_I=', Coupler%nPointTargetCouple_I
             flush(iUnitTest)
          end if
       end if
       if(is_proc(Coupler%iCompSource)) then
          call compact_router(nProcTarget, &
               Coupler%iProcTarget_P, nPointSource_P, &
               Coupler%nCoupleSource, &
               Coupler%iProcSourceCouple_I, Coupler%nPointSourceCouple_I)
          if(DoDebug)then
             write(iUnitTest,*) NameSourceTarget, &
                  ' nPointSource_P=', nPointSource_P
             write(iUnitTest,*) NameSourceTarget, &
                  ' new nCoupleSource=', Coupler%nCoupleSource
             write(iUnitTest,*) NameSourceTarget, &
                  ' new iProcSourceCouple_I =', Coupler%iProcSourceCouple_I
             write(iUnitTest,*) NameSourceTarget, &
                  ' new nPointSourceCouple_I=', Coupler%nPointSourceCouple_I
             flush(iUnitTest)
          end if
       end if

       ! Router is already compacted, information is stored
       deallocate(nPointSource_P, nPointTarget_P)

       ! Allocate buffer for positions on source (zero size on target)
       if(is_proc(Coupler%iCompSource)) Coupler%nPointSource = &
            sum(Coupler%nPointSourceCouple_I)
       if(allocated(Coupler%PosSource_DI)) deallocate(Coupler%PosSource_DI)
       allocate(Coupler%PosSource_DI(Coupler%nDim, Coupler%nPointSource))

       ! Send target point positions to the owner source processors
       call transfer_buffer_real(Coupler%iCommUnion, &
            Coupler%nProcUnion, Coupler%iProcUnion, &
            Coupler%nCoupleTarget, Coupler%iProcTargetCouple_I, &
            Coupler%nPointTargetCouple_I, Coupler%nCoupleSource, &
            Coupler%iProcSourceCouple_I, Coupler%nPointSourceCouple_I, &
            Coupler%nDim, Coupler%nData, PosSortTarget_DI, &
            Coupler%nPointSource, Coupler%PosSource_DI)

       deallocate(PosSortTarget_DI)

       if(DoTest)call timing_stop('pnt_route_'//NameSourceTarget)

    end if    ! IsNewRoute

    ! Transfer data from Source to Target

    ! Get the data from source
    if(DoTest)call timing_start('pnt_alloc_'//NameSourceTarget)
    allocate(Coupler%DataSource_VI(Coupler%nVar,Coupler%nPointSource))
    if(DoTest)call timing_stop('pnt_alloc_'//NameSourceTarget)

    if(DoTest)call timing_start('pnt_get_'//NameSourceTarget)
    if(is_proc(Coupler%iCompSource)) call source_get( IsNewRoute, &
         Coupler%NameVar, Coupler%nVar, Coupler%nDim, Coupler%nPointSource, &
         Coupler%PosSource_DI, Coupler%DataSource_VI)
    if(DoTest)call timing_stop('pnt_get_'//NameSourceTarget)

    ! Transfer data to the target processor
    if(DoTest)call timing_start('pnt_alloc_'//NameSourceTarget)
    allocate(Coupler%DataTarget_VI(Coupler%nVar, Coupler%nData))
    if(DoTest)call timing_stop('pnt_alloc_'//NameSourceTarget)

    if(DoTestMe) write(*,*) NameSub,' start transfer buffer'

    if(DoTest)call timing_start('pnt_transfer_'//NameSourceTarget)
    call transfer_buffer_real(Coupler%iCommUnion, &
         Coupler%nProcUnion, Coupler%iProcUnion, &
         Coupler%nCoupleSource, Coupler%iProcSourceCouple_I, &
         Coupler%nPointSourceCouple_I, Coupler%nCoupleTarget, &
         Coupler%iProcTargetCouple_I, Coupler%nPointTargetCouple_I, &
         Coupler%nVar, Coupler%nPointSource, Coupler%DataSource_VI, &
         Coupler%nData, Coupler%DataTarget_VI)

    if(DoTest)call timing_stop('pnt_transfer_'//NameSourceTarget)
    deallocate(Coupler%DataSource_VI)

    if(DoTestMe) write(*,*) NameSub,' start target_put'

    ! Give the data to target
    if(DoTest)call timing_start('pnt_put_'//NameSourceTarget)
    if(is_proc(Coupler%iCompTarget)) call target_put( &
         Coupler%NameVar, Coupler%nVar, Coupler%nPointTarget, &
         Coupler%DataTarget_VI, Coupler%iPointTarget_I)
    if(DoTest)call timing_stop('pnt_put_'//NameSourceTarget)

    deallocate(Coupler%DataTarget_VI)

    if(DoTest) call timing_stop('pnt_cpl_'//NameSourceTarget)
    if(DoTestMe) write(*,*) NameSub,'_',NameSourceTarget,' finished'

  end subroutine couple_points
  !============================================================================
  subroutine transfer_buffer_real(&
       iComm, nProc, iProc, &
       nCoupleS, iProcS_I, nPointS_I,&
       nCoupleR, iProcR_I, nPointR_I,&
       nData, &
       nBufferS, BufferS_I, nBufferR, BufferR_I)

    ! This subroutine transfers array of reals
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present

    use ModMpi

    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index

    integer, intent(in):: nCoupleS            ! #   of procs to send to
    integer, intent(in):: iProcS_I( nCoupleS) ! ids of procs to send to
    integer, intent(in):: nPointS_I(nCoupleS) ! # of points  to send

    integer, intent(in):: nCoupleR            ! #   of procs to recv from
    integer, intent(in):: iProcR_I( nCoupleR) ! ids of procs to recv from
    integer, intent(in):: nPointR_I(nCoupleR) ! # of points  to recv

    integer, intent(in)   :: nData            ! number of items per point
    integer, intent(in)   :: nBufferS                  ! send buffer size
    real,    intent(inout):: BufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                  ! recv buffer size
    real,    intent(inout):: BufferR_I(nData*nBufferR) ! recv buffer

    integer:: iProcR, iProcS ! R/S is proc to Recv from/Send to
    integer:: iCouple, iBuffer, nCoupleAll

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:)
    integer:: iRequest
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer_real'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData

    nCoupleAll = nCoupleS + nCoupleR

    if(nCoupleAll == 0) RETURN

    allocate(iRequest_I(nCoupleAll))

    ! Counter for sends and receives
    iRequest = 0

    ! Send buffer
    if(nCoupleS > 0)then

       ! Send array of reals
       iBuffer = 1
       do iCouple = 1, nCoupleS
          if(nPointS_I(iCouple) == 0) CYCLE
          iProcS = iProcS_I(iCouple)
          iRequest = iRequest + 1
          call MPI_isend(&
               BufferS_I(iBuffer), nData*nPointS_I(iCouple), &
               MPI_REAL, iProcS, iTag, iComm, iRequest_I(iRequest),iError)
          iBuffer = iBuffer + nData*nPointS_I(iCouple)
       end do
    end if
    if(nCoupleR > 0)then

       ! Receive array of reals
       iBuffer = 1
       do iCouple = 1, nCoupleR
          if(nPointR_I(iCouple) == 0) CYCLE
          iProcR = iProcR_I(iCouple)
          iRequest = iRequest + 1
          call MPI_irecv(&
               BufferR_I(iBuffer), nData*nPointR_I(iCouple),&
               MPI_REAL, iProcR, iTag, iComm, &
               iRequest_I(iRequest),iError)
          iBuffer = iBuffer + nData*nPointR_I(iCouple)
       end do
    end if

    ! Finalize transfer
    call MPI_waitall(iRequest, iRequest_I, MPI_STATUSES_IGNORE, iError)
    deallocate(iRequest_I)

  end subroutine transfer_buffer_real
  !============================================================================
  subroutine transfer_buffer_int(&
       iComm, nProc, iProc,&
       nCoupleS, iProcS_I, nPointS_I,&
       nCoupleR, iProcR_I, nPointR_I,&
       nData, &
       nBufferS, iBufferS_I, nBufferR, iBufferR_I)

    ! This subroutine transfers array of integers
    ! A processor can both send and recv
    ! Correctness: only buffers of the same type should be present

    use ModMpi

    integer, intent(in):: iComm  ! MPI communicator
    integer, intent(in):: nProc  ! number of processors
    integer, intent(in):: iProc  ! local proc index

    integer, intent(in):: nCoupleS                 ! #   of procs to send to
    integer, intent(in):: iProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nPointS_I(nCoupleS)! # of points  to send

    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nPointR_I(nCoupleR)! # of points  to recv

    integer, intent(in   ):: nData                 ! number of items per point
    integer, intent(in   ):: nBufferS                   ! send buffer size
    integer, intent(inout):: iBufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                   ! recv buffer size
    integer, intent(inout):: iBufferR_I(nData*nBufferR) ! recv buffer

    integer:: iProcR, iProcS ! R/S is proc to Recv from/Send to
    integer:: iCouple, iBuffer, nCoupleAll
    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:)
    integer:: iRequest
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'transfer_buffer_int'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub, 'called with nProc, iProc, nData=', &
         nProc, iProc, nData

    nCoupleAll = nCoupleS + nCoupleR
    if(nCoupleAll == 0) RETURN

    allocate(iRequest_I(nCoupleAll))
    iRequest = 0

    ! Send buffer
    if(nCoupleS > 0)then

       ! Send array of integers
       iBuffer = 1
       do iCouple = 1, nCoupleS
          if(nPointS_I(iCouple) == 0)CYCLE
          iProcS = iProcS_I(iCouple)
          iRequest = iRequest + 1
          call MPI_Isend(&
               iBufferS_I(iBuffer), nData*nPointS_I(iCouple), &
               MPI_INTEGER, iProcS, iTag, iComm, &
               iRequest_I(iRequest),iError)
          iBuffer = iBuffer + nData*nPointS_I(iCouple)
       end do
    end if
    if(nCoupleR > 0)then

       ! Receive array of integers
       iBuffer = 1
       do iCouple = 1, nCoupleR
          if(nPointR_I(iCouple) == 0)CYCLE
          iProcR = iProcR_I(iCouple)
          iRequest = iRequest + 1
          call MPI_Irecv(&
               iBufferR_I(iBuffer), nData*nPointR_I(iCouple),&
               MPI_INTEGER, iProcR, iTag, iComm, &
               iRequest_I(iRequest),iError)
          iBuffer = iBuffer + nData*nPointR_I(iCouple)
       end do
    end if

    ! Finalize transfer
    call MPI_waitall(iRequest, iRequest_I, MPI_STATUSES_IGNORE, iError)
    deallocate(iRequest_I)

  end subroutine transfer_buffer_int
  !============================================================================
  subroutine get_buffer_order(&
       nProc, nBuffer, iProc_I, nBuffer_P, iBuffer_I, nData)

    ! nBuffer_P returns the number of points in the buffer belong
    !    to each processor based on the iProc_I information.
    ! iBuffer_I returns the order of the points with increasing
    !    processor index.

    integer, intent(in)::  nProc, nBuffer    ! number of procs and points
    integer, intent(in)::  iProc_I(nBuffer)  ! proc index for each point
    integer, intent(out):: nBuffer_P(0:nProc-1)  ! number of points per procs
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering
    integer, intent(out):: nData ! number of points found on Source

    integer:: iBuffer, iProc
    integer, allocatable:: iBuffer_P(:)
    character(len=*), parameter:: NameSub = 'get_buffer_order'
    !--------------------------------------------------------------------------
    nData = nBuffer
    ! Buffer chunk sizes for each processor
    nBuffer_P = 0
    do iBuffer = 1, nBuffer
       iProc = iProc_I(iBuffer)
       if(iProc < 0) CYCLE
       nBuffer_P(iProc) = nBuffer_P(iProc) + 1
    end do

    ! Set starting point of chunks in the reordered buffer
    allocate(iBuffer_P(0:nProc-1))
    iBuffer_P(0) = 1
    do iProc = 1, nProc-1
       iBuffer_P(iProc) = iBuffer_P(iProc-1) + nBuffer_P(iProc-1)
    end do

    ! Get the order that is consecutive per processor index
    do iBuffer = 1, nBuffer
       ! This buffer belongs to iProc
       iProc = iProc_I(iBuffer)
       if(iProc < 0)then
          ! the point hasn't been found
          iBuffer_I(iBuffer) = - 1
          nData = nData - 1
       else
          ! It will occupy iBufferSort position
          iBuffer_I(iBuffer) = iBuffer_P(iProc)
          ! Jump to next position
          iBuffer_P(iProc) = iBuffer_P(iProc) + 1
       end if
    end do
    deallocate(iBuffer_P)

  end subroutine get_buffer_order
  !============================================================================
  subroutine set_recv_info(iComm, nProc, iProc, &
       nBufferS_P, nBufferR, nBufferR_P)

    ! Based on nBufferS_P (number of points sent to each recv proc)
    ! set nBufferR and nBufferR_P (number of points recv from each send proc)
    ! nBufferR = sum(nBufferR_P)

    integer, intent(in)::  iComm, nProc, iProc ! Router communicator
    integer, intent(in)::  nBufferS_P(nProc)   ! Send buffer chunks
    integer, intent(out):: nBufferR            ! Recv buffer size
    integer, intent(out):: nBufferR_P(nProc)   ! Recv buffer chunks

    integer, parameter:: iTag = 76

    integer, allocatable:: nPointRS_PP(:,:), nPoint_P(:)
    integer:: iProcR, iError
    ! Gather the nBufferS_P information onto the root.
    ! The resulting nPointRS_PP table contains the
    ! the number of points owned by iProcR and iProcS.
    !--------------------------------------------------------------------------
    if(iProc == 0)then
       allocate(nPointRS_PP(0:nProc-1,0:nProc-1), nPoint_P(0:nProc-1))
    else
       allocate(nPointRS_PP(1,1))
    end if

    call MPI_gather(nBufferS_P, nProc, MPI_INTEGER, &
         nPointRS_PP, nProc, MPI_INTEGER, 0, iComm, iError)

    ! Local copy for root
    if(iProc==0) nBufferR_P = nPointRS_PP(0,:)

    ! Now send the columns of the table to the recv processors
    do iProcR = 1, nProc-1
       if(iProc==0)then
          nPoint_P = nPointRS_PP(iProcR,:)
          call MPI_send(nPoint_P, nProc, MPI_INTEGER, iProcR, &
               iTag, iComm, iError)
       end if
       if(iProc == iProcR)then
          call MPI_recv(nBufferR_P, nProc, MPI_integer, 0, &
               iTag, iComm, MPI_STATUS_IGNORE, iError)
       end if
    end do
    if(iProc==0) deallocate(nPoint_P)
    deallocate(nPointRS_PP)

    ! Total number of points owned by this recv proc
    nBufferR = sum(nBufferR_P)

  end subroutine set_recv_info
  !============================================================================
  subroutine set_recv_info2( &
       iComm, iProc, iProc0S, iProcR_P, iCompS, iCompR, &
       nBufferS_P, nBufferR, nBufferR_P)

    ! Same as set_recv_info but for different layouts
    ! Based on nBufferS_P (number of points sent to each recv proc)
    ! set nBufferR and nBufferR_P (number of points recv from each send proc)
    ! nBufferR = sum(nBufferR_P)

    integer, intent(in)::  iComm, iProc    ! Router communicator
    integer, intent(in)::  iProc0S         ! Root of sender in router
    integer, intent(in)::  iProcR_P(0:)    ! Recv proc index in router
    integer, intent(in)::  iCompR, iCompS  ! Recv and send components
    integer, intent(in)::  nBufferS_P(:)   ! Send buffer chunks
    integer, intent(out):: nBufferR        ! Recv buffer size
    integer, intent(out):: nBufferR_P(:)   ! Recv buffer chunks

    integer, parameter:: iTag = 76

    integer:: iCommS, iProcS, nProcS, jProcR, nProcR

    integer, allocatable:: nPointRS_PP(:,:), nPointS_P(:)
    integer:: iError

    character(len=*), parameter:: NameSub = 'set_recv_info2'
    !--------------------------------------------------------------------------
    iCommS = i_comm(iCompS)
    iProcS = i_proc(iCompS)
    nProcS = n_proc(iCompS)
    nProcR = n_proc(iCompR)

    nBufferR_P = 0 ! initialize on all processors

    if(iProcS >= 0)then
       ! Gather the nBufferS_P information onto the root of the sender.
       ! The resulting nPointRS_PP table contains the
       ! the number of points transfered from iProcS to iProcR.
       if(iProcS == 0)then
          allocate(nPointRS_PP(0:nProcR-1,0:nProcS-1), nPointS_P(0:nProcS-1))
       else
          allocate(nPointRS_PP(1,1))
       end if

       call MPI_gather(nBufferS_P, nProcR, MPI_INTEGER, &
            nPointRS_PP, nProcR, MPI_INTEGER, 0, iCommS, iError)

       if(DoDebug) write(*,*) NameSub, ' nPointRS_PP=', nPointRS_PP
    end if

    ! Now send the nProcR columns of the table from iCompS root
    ! to the recv processors
    do jProcR = 0, nProcR - 1
       if(iProcS == 0)then
          if(iProc == iProcR_P(jProcR))then
             ! local copy (avoids MPI hanging)
             nBufferR_P = nPointRS_PP(jProcR,:)
          else
             ! Extract column so it is continuous in memory
             nPointS_P = nPointRS_PP(jProcR,:)
             ! Translate local rank jProcR to union rank
             call MPI_send(nPointS_P, nProcS, MPI_INTEGER, &
                  iProcR_P(jProcR), iTag, iComm, iError)
          end if
       end if
       if(iProc == iProcR_P(jProcR) .and. iProcS /= 0)then
          ! Recv from the root of iCommS
          call MPI_recv(nBufferR_P, nProcS, MPI_integer, iProc0S, &
               iTag, iComm, MPI_STATUS_IGNORE, iError)
       end if
    end do

    if(iProcS == 0) deallocate(nPointS_P)
    if(iProcS >= 0) deallocate(nPointRS_PP) !!! keep array to save time ?!?

    ! Total number of points owned by this recv proc
    nBufferR = sum(nBufferR_P)

  end subroutine set_recv_info2
  !============================================================================
  subroutine set_inquiry(iProc,                                  &
       iCompSource, iCompTarget, nProcCommon,                    &
       nCoupleTarget, iProcTargetCouple_I, nPointTargetCouple_I, &
       nCoupleSource, iProcSourceCouple_I, nPointSourceCouple_I)

    ! the subroutine returns number of processors to be communicated with
    !     nCoupleSource, nCoupleTarget
    ! and stores indexes of these processors in
    !     iProcSourceCouple_I, iProcTargetCouple_I
    ! iCompSource is the component which Data is stored on
    ! iCompTarget is the component which Data is to be sent to
    ! in the Union communicator procs on Source are assumed to come first

    integer,              intent(in ):: iProc
    integer,              intent(in ):: iCompSource, iCompTarget
    integer,              intent(in ):: nProcCommon

    integer,              intent(out):: nCoupleTarget
    integer, allocatable, intent(out):: iProcTargetCouple_I(:)
    integer, allocatable, intent(out):: nPointTargetCouple_I(:)

    integer,              intent(out):: nCoupleSource
    integer, allocatable, intent(out):: iProcSourceCouple_I(:)
    integer, allocatable, intent(out):: nPointSourceCouple_I(:)

    integer:: DnProc ! remainder of division (nCompSource/nCompTarget)**(+/-1)
    integer:: iProcLocal, iCouple
    integer:: nProcTarget, nProcSource
    integer:: nCouple, nCoupleOther, nProc0Other
    integer, allocatable:: iProc_I(:)
    !--------------------------------------------------------------------------
    nProcSource = n_proc(iCompSource)

    ! Here shared procs are considered to be on Source and not on Target
    ! they are processed separately in order to send inquiry to itself
    nProcTarget = n_proc(iCompTarget) - nProcCommon

    if(nProcTarget < 1)then
       ! there are no procs on Target only
       if(is_proc(iCompTarget))then
          ! Target and source procs coincide, communicates with itself
          nCoupleTarget = 1
          nCoupleSource = 1
       elseif(is_proc(iCompSource))then
          ! A source proc that is not involved
          nCoupleTarget = 0
          nCoupleSource = 0
       end if
       allocate(iProcTargetCouple_I(nCoupleTarget))
       allocate(iProcSourceCouple_I(nCoupleSource))

       if(nCoupleSource*nCoupleTarget > 0)then
          iProcTargetCouple_I(nCoupleTarget) = iProc
          iProcSourceCouple_I(nCoupleSource) = iProc
       end if

       allocate(nPointTargetCouple_I(nCoupleTarget))
       allocate(nPointSourceCouple_I(nCoupleSource))

       RETURN
    end if

    ! if nProcSource/nProcTarget is larger and not an exact multiple of the
    ! other, there is a discrepancy of DnProc processors
    if(nProcSource > nProcTarget)then
       DnProc = mod(nProcSource, nProcTarget)
    else
       DnProc = mod(nProcTarget, nProcSource)
    end if

    ! in the communicator iComm processor are ordered as follows
    ! [ SOURCE | TARGET ]
    ! find the local index of iProc on its component,
    ! index of 0-processor of the other component,
    ! # of procs to couple with: nCouple,
    ! # of procs other component is coupled with: nCoupleOther
    if(is_proc(iCompSource))then
       nCouple      = nProcTarget / nProcSource
       nCoupleOther = nProcSource / nProcTarget
       iProcLocal   = iProc
       nProc0Other  = nProcSource
    elseif(is_proc(iCompTarget))then
       nCouple      = nProcSource / nProcTarget
       nCoupleOther = nProcTarget / nProcSource
       iProcLocal   = iProc - nProcSource
       nProc0Other  = 0
    end if

    ! account for discrepancy
    if(iProcLocal < DnProc*(nCoupleOther + 1))then
       DnProc = 0
       nCoupleOther = nCoupleOther + 1
       nCouple = nCouple + 1
    elseif(nCouple == 0)then
       DnProc =  - DnProc
    end if

    nCouple      = max(1, nCouple)
    nCoupleOther = max(1, nCoupleOther)

    allocate(iProc_I(nCouple))

    do iCouple = 1, nCouple
       iProc_I(iCouple) = &
            nProc0Other + (iCouple - 1) + &
            (nCouple * iProcLocal + DnProc)/nCoupleOther
    end do

    if(  is_proc(iCompSource).and.&
         is_proc(iCompTarget)) then
       nCoupleTarget = 1
       nCoupleSource = nCouple + 1
       allocate(iProcTargetCouple_I(nCoupleTarget))
       allocate(iProcSourceCouple_I(nCoupleSource))
       iProcTargetCouple_I(1)                 = iProc
       iProcSourceCouple_I(nCoupleSource)     = iProc
       iProcSourceCouple_I(1:nCoupleSource-1) = iProc_I
    elseif(is_proc(iCompSource))then
       nCoupleTarget = 0
       nCoupleSource = nCouple
       allocate(iProcTargetCouple_I(nCoupleTarget))
       allocate(iProcSourceCouple_I(nCoupleSource))
       iProcSourceCouple_I = iProc_I
    elseif(is_proc(iCompTarget))then
       nCoupleTarget = nCouple
       nCoupleSource = 0
       allocate(iProcTargetCouple_I(nCoupleTarget))
       allocate(iProcSourceCouple_I(nCoupleSource))
       iProcTargetCouple_I = iProc_I
    end if

    deallocate(iProc_I)

    allocate(nPointTargetCouple_I(nCoupleTarget))
    allocate(nPointSourceCouple_I(nCoupleSource))

  end subroutine set_inquiry
  !============================================================================
  subroutine compact_router( &
       nProc, iProc_P, nPoint_P, nCouple, iProc_I, nPoint_I)

    ! Compact the nPoint_P(nProc) array to nPoint_I(nCouple)
    ! containng the non-zero elements only.

    integer, intent(in):: nProc
    integer, intent(in):: iProc_P(0:nProc-1) ! union index of component procs
    integer, intent(in):: nPoint_P(0:nProc-1)

    integer, intent(out):: nCouple                    ! number of non-0 entries
    integer, allocatable, intent(inout):: iProc_I(:)  ! proc index for coupler
    integer, allocatable, intent(inout):: nPoint_I(:) ! size of coupler

    integer:: iProc, iCouple

    ! Number of procs with data
    !--------------------------------------------------------------------------
    nCouple = count(nPoint_P /= 0)

    if(allocated(iProc_I))  deallocate(iProc_I);  allocate(iProc_I(nCouple))
    if(allocated(nPoint_I)) deallocate(nPoint_I); allocate(nPoint_I(nCouple))

    if(nCouple == 0) RETURN

    ! Fill in compacted arrays
    iCouple = 0
    do iProc = 0, nProc-1
       ! Skip processors that have no data to send.
       if(nPoint_P(iProc) == 0) CYCLE
       iCouple = iCouple + 1
       iProc_I(iCouple)  = iProc_P(iProc)
       nPoint_I(iCouple) = nPoint_P(iProc)
    end do

  end subroutine compact_router
  !============================================================================
  subroutine get_recv_buffer_size(iComm, nProc, iProc,&
       nCoupleS, iProcS_I, nPointS_I,&
       nCoupleR, iProcR_I, nPointR_I)

    ! This subroutine determines how large is the buffer
    ! to be transfered to recv component processor

    integer, intent(in)   :: iComm                  ! MPI communicator
    integer, intent(in)   :: nProc                  ! number of processors
    integer, intent(in)   :: iProc                  ! local proc index

    integer, intent(in):: nCoupleS                 !
    integer, intent(in):: iProcS_I( nCoupleS)! Couple procs
    integer, intent(in):: nPointS_I(nCoupleS)! # of points to send/recv

    integer, intent(in ):: nCoupleR                 !
    integer, intent(in ):: iProcR_I( nCoupleR)!
    integer, intent(out):: nPointR_I(nCoupleR)!

    ! index of recv and send processors
    integer:: iProcR, iProcS
    integer:: iCouple

    ! MPI stuff
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'get_recv_buffer_size'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub, 'called with nProc, iProc=', &
         nProc, iProc
    if(nCoupleS + nCoupleR == 0) RETURN

    allocate(iRequest_I(nCoupleS+nCoupleR))
    if(nCoupleS > 0)then

       ! Send buffer size
       do iCouple = 1, nCoupleS
          iProcS = iProcS_I(iCouple)
          call MPI_isend(&
               nPointS_I(iCouple), 1, MPI_INTEGER, &
               iProcS, iTag, iComm, iRequest_I(iCouple), iError)
       end do
    end if
    if(nCoupleR > 0)then

       ! Recv buffer size
       do iCouple = 1, nCoupleR
          iProcR = iProcR_I(iCouple)
          call MPI_irecv(&
               nPointR_I(iCouple), 1, MPI_INTEGER, &
               iProcR, iTag, iComm, iRequest_I(nCoupleS+iCouple), iError)
       end do
    end if

    ! Finalize transfer
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I, &
         MPI_STATUSES_IGNORE,iError)
    deallocate(iRequest_I)

  end subroutine get_recv_buffer_size
  !============================================================================

end module CON_couple_points
!==============================================================================
