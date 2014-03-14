!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PT

!BOP
!MODULE: CON_couple_gm_pt - couple GM and PT components
!
!DESCRIPTION:
! Couple GM and PT components both ways. 
!
!INTERFACE:
module CON_couple_gm_pt

  !USES:
  use CON_coupler

  use CON_couple_points

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pt_init ! initialize both couplings
  public :: couple_gm_pt      ! couple GM to PT

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe = .true., IsInitialized = .false.

  ! Check if the GM and PT processors coincide
  logical:: IsSameLayout

  ! Proc index inside SWMF
  integer:: iProcWorld

  ! Router communicator info
  integer:: iCommGmPt, nProcGmPt, iProcGmPt, iProc0Gm, iProc0Pt, nProcCommon

  type(CouplePointsType) :: Coupler
  
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt_init - initialize GM-PT couplings
  !INTERFACE:
  subroutine couple_gm_pt_init

    integer:: iError
    !integer:: nDimPt

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_pt_init'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    Coupler%iCompTarget = PT_
    Coupler%iCompSource = GM_
    
    if(IsInitialized) RETURN
    IsInitialized = .true.

    call couple_points_init(Coupler)

  end subroutine couple_gm_pt_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_pt - couple GM to PT
  !INTERFACE:
  subroutine couple_gm_pt(tSimulation)
    
    interface

       subroutine PT_put_from_gm( &
            NameVar, nVar, nPoint, Pos_DI, Data_VI, iPoint_I)

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
       end subroutine PT_put_from_gm
       subroutine GM_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

         implicit none

         ! dimension of position vectors
         integer, intent(in) :: nDimIn             

         ! number of positions
         integer, intent(in) :: nPoint                

         ! positions
         real,    intent(in) :: Xyz_DI(nDimIn,nPoint) 

         ! processor owning position
         integer, intent(out):: iProc_I(nPoint)       

       end subroutine GM_find_points


        subroutine GM_get_for_pt(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
             Data_VI)
          logical,          intent(in):: IsNew   ! true for new point array
          character(len=*), intent(in):: NameVar ! List of variables
          integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
          integer,          intent(in):: nDimIn  ! Dimensionality of positions
          integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

          real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
          real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array
        end subroutine GM_get_for_pt

        subroutine GM_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine GM_get_grid_info

        subroutine PT_get_grid_info(nDim, iGrid, iDecomp)
          integer, intent(out) :: nDim, iGrid, iDecomp
        end subroutine PT_get_grid_info
     end interface
    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere       (GM) source\\
    !    Particle Tracker           (PT) target
    !
    ! Send information from GM to PT. 
    !EOP

    ! Stored variables for GM->PT coupler
    !------------------------------------

    ! List of variables to pass
    character(len=lNameVar):: NameVar

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1
    integer:: iGridGm, iDecompGm, iGridPt, iDecompPt

    ! Number of local points for PT and GM cores
    integer:: nPointPt = 0, nPointGm = 0
    
    ! Number of data points found on GM component
    integer:: nData = 0

    ! Number of points that belong to a given processor of the OTHER component
    integer, allocatable, save:: nPointGm_P(:), nPointPt_P(:)

    ! Permutation of PT points after data is returned
    integer, allocatable, save:: iPointPt_I(:)

    ! Point positions local on a GM processor
    real, allocatable, save:: PosGm_DI(:,:)
    !---------------------------------------------------------------
    ! Temporary variables

    ! Is there a need to recalculate the data transfer route?
    logical:: IsNewRoute

    ! Storage for original PT point positions. 
    real, pointer:: PosPt_DI(:,:)

    ! Positions sorted according to the correspongin GM processors
    real, allocatable:: PosSortPt_DI(:,:)

    ! GM processor index for PT points
    integer, allocatable:: iProcPt_I(:)
    integer, allocatable:: iProcGm_I(:)

    ! Buffers for data on GM and PT
    real, allocatable:: DataGm_VI(:,:), DataPt_VI(:,:)

    integer:: iPoint, i, iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pt'

    ! Variables for general-case coupling

    ! number of processors of the OTHER component to communicate with
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    ! GM sends all its variables to PT. Take information from Grid_C
    Coupler%NameVar = Grid_C(Coupler%iCompSource)%NameVar
    Coupler%nVar    = Grid_C(Coupler%iCompSource)%nVar

    if(DoTest)write(*,*)NameSub,' starting iProc=',Coupler%iProcWorld

    write(*,*)NameSub,' starting iProc=',Coupler%iProcWorld

    ! After everything is initialized exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',Coupler%iProcWorld

    ! Allocate arrays for router
    if(.not.allocated(Coupler%nPointSource_P) .and. Coupler%IsSameLayout) then
       allocate(Coupler%nPointSource_P(Coupler%nProcSourceTarget), Coupler%nPointTarget_P(Coupler%nProcSourceTarget))
    elseif(.not.allocated(Coupler%nPointSource_P)) then
       allocate(Coupler%nPointSource_P(n_proc(Coupler%iCompTarget)), Coupler%nPointTarget_P(n_proc(Coupler%iCompSource)))
    end if

    ! Check if there is a need to redo the mapping
    call GM_get_grid_info(Coupler%nDim, iGridGm, iDecompGm)
    call PT_get_grid_info(Coupler%nDim, iGridPt, iDecompPt)

    if(.not.Coupler%IsSameLayout)then
       call MPI_bcast(&
            iDecompGm, 1, MPI_INTEGER, i_proc0(Coupler%iCompSource), Coupler%iCommSourceTarget, iError)
       call MPI_bcast(&
            iDecompPt, 1, MPI_INTEGER, i_proc0(Coupler%iCompTarget), Coupler%iCommSourceTarget, iError)
    endif

    IsNewRoute = iDecompGm /= iDecompLastGm .or. iDecompLastPt /= iDecompPt

    iDecompLastGm = iDecompGm
    iDecompLastPt = iDecompPt

    if(IsNewRoute)then
       Coupler%nPointTarget = 0
       Coupler%nPointSource = 0
       Coupler%nData    = 0
       ! Get positions where info is needed from PT. PT will allocate array.
       if(is_proc(Coupler%iCompTarget))then
          call PT_put_from_gm(Coupler%NameVar, Coupler%nVar, Coupler%nPointTarget, Coupler%PosTarget_DI)
       else
          Coupler%nPointTarget = 0
          allocate(Coupler%PosTarget_DI(Coupler%nDim,0))
       end if

       if(Coupler%IsSameLayout)then

          ! Find processors that own the PT positions in GM
          allocate(iProcPt_I(Coupler%nPointTarget))

          call GM_find_points(Coupler%nDim, Coupler%nPointTarget, Coupler%PosTarget_DI, iProcPt_I)

          ! Order points according to the owner processor indexes
          if(allocated(Coupler%iPointTarget_I)) deallocate(Coupler%iPointTarget_I)
          allocate(Coupler%iPointTarget_I(Coupler%nPointTarget))
          call get_buffer_order(Coupler%nProcSourceTarget, Coupler%nPointTarget, iProcPt_I, &
               Coupler%nPointTarget_P, Coupler%iPointTarget_I, Coupler%nData)

          deallocate(iProcPt_I)

          ! Rearrange coordinate array according to processor order
          allocate(PosSortPt_DI(Coupler%nDim,Coupler%nData))
          do iPoint = 1, Coupler%nPointTarget
             if(Coupler%iPointTarget_I(iPoint) < 0)CYCLE
             PosSortPt_DI(:,Coupler%iPointTarget_I(iPoint)) = Coupler%PosTarget_DI(:,iPoint)
          end do
          deallocate(Coupler%PosTarget_DI)

          ! Set number of points to be received on the GM component
          call set_recv_info(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, &
               Coupler%nPointTarget_P, Coupler%nPointSource, Coupler%nPointSource_P)

          ! Transfer PT positions to the GM processors that own them
          if(allocated(Coupler%PosSource_DI)) deallocate(Coupler%PosSource_DI)
          allocate(Coupler%PosSource_DI(Coupler%nDim,Coupler%nPointSource))
          call transfer_buffer(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, Coupler%nDim, &
               Coupler%nData,    Coupler%nPointTarget_P, PosSortPt_DI, &
               Coupler%nPointSource, Coupler%nPointSource_P, Coupler%PosSource_DI)

          deallocate(PosSortPt_DI)
       else
          !\
          ! Layouts ovelap
          !/
          if(allocated(Coupler%iCoupleProcSource_I ))deallocate(Coupler%iCoupleProcSource_I )
          if(allocated(Coupler%iCoupleProcTarget_I ))deallocate(Coupler%iCoupleProcTarget_I )
          if(allocated(Coupler%nCouplePointSource_I))deallocate(Coupler%nCouplePointSource_I)
          if(allocated(Coupler%nCouplePointTarget_I))deallocate(Coupler%nCouplePointTarget_I)
          if(allocated(Coupler%iPointTarget_I      ))deallocate(Coupler%iPointTarget_I      )
          if(allocated(Coupler%PosSource_DI        ))deallocate(Coupler%PosSource_DI        )

          ! Setup communication pattern for finding points on GM
          ! This allocates and sets iCoupleProc*_I, nCouplePoint*_I arrays
          call set_inquiry(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget,  &
               Coupler%iCompSource, Coupler%iCompTarget, nProcCommon,                        &
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I, &
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I)

          if(is_proc(Coupler%iCompTarget))then
             ! Number of points to send from PT to the selected GM processors
             Coupler%nCouplePointTarget_I(1:Coupler%nCoupleTarget-1) = Coupler%nPointTarget / Coupler%nCoupleTarget 
             ! Last processor gets the rest of points
             Coupler%nCouplePointTarget_I(Coupler%nCoupleTarget) = &
                  Coupler%nPointTarget - sum(Coupler%nCouplePointTarget_I(1:Coupler%nCoupleTarget-1))
          end if

          ! send number of points from PT to GM that will be sent to 
          ! "GM_find_points." result is in Coupler%nCouplePointSource_I (number of 
          ! points to be recieved)
          call get_recv_buffer_size(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, &
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I,         &
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I)

          ! Allocate buffer for positions on GM (zero size on PT)
          if(is_proc(Coupler%iCompSource)) Coupler%nPointSource = sum(Coupler%nCouplePointSource_I)
          allocate(Coupler%PosSource_DI(Coupler%nDim,Coupler%nPointSource))

          ! Send positions from PT to GM, so GM can find the owners
          call transfer_buffer_real(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget,&
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I,          &
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I,          &
               Coupler%nDim,                                                  &
               Coupler%nPointTarget, Coupler%PosTarget_DI,&
               Coupler%nPointSource, Coupler%PosSource_DI)

          ! Find processors that own the PT positions in GM
          allocate(iProcGm_I(Coupler%nPointSource))
          if(is_proc(Coupler%iCompSource))&
               call GM_find_points(Coupler%nDim, Coupler%nPointSource, Coupler%PosSource_DI, iProcGm_I)
          deallocate(Coupler%PosSource_DI)

          ! send owner processor indexes from GM to PT
          allocate(iProcPt_I(Coupler%nPointTarget))
          call transfer_buffer_int(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget,&
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I,&
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I,&
               1, &
               Coupler%nPointSource, iProcGm_I,&
               Coupler%nPointTarget, iProcPt_I)

          ! based on owner information set up the final communication pattern
          call  set_data_transfer(Coupler%iCommSourceTarget,                 &
               Coupler%iCompSource, Coupler%iCompTarget,                                     &
               Coupler%nPointSource, iProcGm_I,                          &
               Coupler%nPointTarget, iProcPt_I,                          &
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I, &
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I)


          if(is_proc(Coupler%iCompTarget))then
             ! Order points according to the owner processor indexes
             ! sort iProcPt_I according to order of procs in Coupler%iCoupleProcTarget_I
             allocate(Coupler%iPointTarget_I(Coupler%nPointTarget))
             call get_transfer_buffer_order(Coupler%iCompSource, Coupler%iCompTarget,           &
                  Coupler%nPointTarget, iProcPt_I,                          &
                  Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I, &
                  Coupler%iPointTarget_I, Coupler%nData)

             allocate(PosSortPt_DI(Coupler%nDim, Coupler%nData))

             ! Rearrange coordinate array according to processor order
             do iPoint = 1, Coupler%nPointTarget
                if(Coupler%iPointTarget_I(iPoint) < 0) CYCLE
                PosSortPt_DI(:,Coupler%iPointTarget_I(iPoint)) = Coupler%PosTarget_DI(:,iPoint)
             end do
          else
             allocate(PosSortPt_DI(Coupler%nDim, Coupler%nData))
          end if

          deallocate(iProcPt_I, Coupler%PosTarget_DI, iProcGm_I)

          ! Allocate buffer for positions on GM (zero size on PT) 
          if(is_proc(Coupler%iCompSource)) Coupler%nPointSource = sum(Coupler%nCouplePointSource_I)
          if(allocated(Coupler%PosSource_DI)) deallocate(Coupler%PosSource_DI)
          allocate(Coupler%PosSource_DI(Coupler%nDim, Coupler%nPointSource))

          ! Send PT point positions to the owner GM processors
          call transfer_buffer_real(Coupler%iCommSourceTarget, Coupler%nProcSourceTarget, Coupler%iProcSourceTarget, &
               Coupler%nCoupleTarget, Coupler%iCoupleProcTarget_I, Coupler%nCouplePointTarget_I,           &
               Coupler%nCoupleSource, Coupler%iCoupleProcSource_I, Coupler%nCouplePointSource_I,           &
               Coupler%nDim,                                                   &
               Coupler%nData,    PosSortPt_DI,              &
               Coupler%nPointSource, Coupler%PosSource_DI)

          deallocate(PosSortPt_DI)

       end if
    end if

    call couple_points(Coupler,GM_get_for_pt, PT_put_from_gm, PT_get_grid_info, GM_get_grid_info)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',Coupler%iProcWorld
  end subroutine couple_gm_pt

  !==========================================================================

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
    !----------------------------------------------------------------------
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

  !==========================================================================

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
    integer:: iProcR
    integer:: iError, iStatus_I(MPI_STATUS_SIZE)
    !-----------------------------------------------------------------------

    ! Gather the nBufferS_P information onto the root.
    ! The resulting nPointRS_PP table contains the
    ! the number of points owned by iProcR and iProcS.
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
               iTag, iComm, iStatus_I, iError)
       end if
    end do
    if(iProc==0) deallocate(nPoint_P)
    deallocate(nPointRS_PP)

    ! Total number of points owned by this recv proc
    nBufferR = sum(nBufferR_P)

  end subroutine set_recv_info


  !===========================================================================

  subroutine set_inquiry(iComm, nProc, iProc,                    &
       iCompSource, iCompTarget, nProcCommon,                    &
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I, &
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I)

    ! the subroutine returns number of processors to be communicated with
    !     nCoupleSource, nCoupleTarget 
    ! and stores indexes of these processors in 
    !     iCoupleProcSource_I, iCoupleProcTarget_I
    ! iCompSource is the component which Data is stored on
    ! iCompTarget is the component which Data is to be sent to
    ! in the Union communicator procs on Source are assumed to come first

    integer,              intent(in ):: iComm, nProc, iProc
    integer,              intent(in ):: iCompSource, iCompTarget
    integer,              intent(in ):: nProcCommon

    integer,              intent(out):: nCoupleTarget
    integer, allocatable, intent(out):: iCoupleProcTarget_I(:)
    integer, allocatable, intent(out):: nCouplePointTarget_I(:)

    integer,              intent(out):: nCoupleSource
    integer, allocatable, intent(out):: iCoupleProcSource_I(:)
    integer, allocatable, intent(out):: nCouplePointSource_I(:)

    integer:: dProc ! remainder of division (nCompSource/nCompTarget)**(+/-1)
    integer:: iProcLocal, iCouple
    integer:: nProcTarget, nProcSource
    integer:: nCouple, nCoupleOther, nProc0Other
    integer, allocatable:: iCoupleProc_I(:)
    integer:: iError
    !------------------------------------------------------------------
    ! Here shared procs are considered to be on Source and not on Target
    ! they are processed separately in order to send inquiry to itself
    nProcSource = n_proc(iCompSource)
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
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))

       if(nCoupleSource*nCoupleTarget > 0)then
          iCoupleProcTarget_I(nCoupleTarget) = iProc
          iCoupleProcSource_I(nCoupleSource) = iProc
       end if

       allocate(nCouplePointTarget_I(nCoupleTarget))
       allocate(nCouplePointSource_I(nCoupleSource))

       RETURN
    end if


    ! if nProcSource/nProcTarget is larger and not an exact multiple of the 
    ! other, there is a discrepancy of dProc processors
    if(nProcSource > nProcTarget)then
       dProc = MOD(nProcSource, nProcTarget)
    else
       dProc = MOD(nProcTarget, nProcSource)
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
    if(iProcLocal < dProc * (nCoupleOther+1))then
       dProc = 0
       nCoupleOther = nCoupleOther + 1
       nCouple = nCouple + 1
    elseif(nCouple == 0)then
       dProc =  - dProc
    end if
    
    nCouple      = MAX(1, nCouple)
    nCoupleOther = MAX(1, nCoupleOther)

    allocate(iCoupleProc_I( nCouple))

    do iCouple = 1, nCouple
       iCoupleProc_I(iCouple) = &
            nProc0Other + (iCouple - 1) + &
            (nCouple * iProcLocal + dProc)/nCoupleOther
    end do

    if(  is_proc(iCompSource).and.&
         is_proc(iCompTarget))then
       nCoupleTarget = 1
       nCoupleSource = nCouple + 1
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcTarget_I(1            ) = iProc
       iCoupleProcSource_I(nCoupleSource) = iProc
       iCoupleProcSource_I(1:nCoupleSource-1) = iCoupleProc_I(:)
    elseif(is_proc(iCompSource))then
       nCoupleTarget = 0
       nCoupleSource = nCouple
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcSource_I(:) = iCoupleProc_I(:)
    elseif(is_proc(iCompTarget))then
       nCoupleTarget = nCouple
       nCoupleSource = 0
       allocate(iCoupleProcTarget_I(nCoupleTarget))
       allocate(iCoupleProcSource_I(nCoupleSource))
       iCoupleProcTarget_I(:) = iCoupleProc_I(:)
    end if

    deallocate(iCoupleProc_I)

    allocate(nCouplePointTarget_I(nCoupleTarget))
    allocate(nCouplePointSource_I(nCoupleSource))

  end subroutine set_inquiry

  !===========================================================================

  subroutine get_recv_buffer_size(iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I)
    ! This subroutine determines how large is the buffer
    ! to be transfered to recv component processor
    integer, intent(in)   :: iComm                  ! MPI communicator
    integer, intent(in)   :: nProc                  ! number of processors
    integer, intent(in)   :: iProc                  ! local proc index
    !-----------------------------------------------------------------------
    integer, intent(in):: nCoupleS                 ! 
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! Couple procs
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points to send/recv
    !-----------------------------------------------------------------------
    integer, intent(in ):: nCoupleR                 !
    integer, intent(in ):: iCoupleProcR_I( nCoupleR)!
    integer, intent(out):: nCouplePointR_I(nCoupleR)!
    !-----------------------------------------------------------------------
    ! index of recv and send processors
    integer:: iProcR, iProcS 
    integer:: iCouple

    ! MPI stuff
    integer:: iStatus(MPI_STATUS_SIZE)
    integer, parameter:: iTag = 77
    integer:: iError
    integer, allocatable:: iRequest_I(:), iStatus_II(:,:)
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'get_recv_buffer_size'
    !-----------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*)NameSub, 'called with nProc, iProc=', &
         nProc, iProc
    if(nCoupleS + nCoupleR == 0) RETURN

    allocate(iRequest_I(nCoupleS+nCoupleR))
    allocate(iStatus_II(MPI_STATUS_SIZE,nCoupleS+nCoupleR))
    if(nCoupleS > 0)then
       !\
       ! Send buffer size
       !/
       do iCouple = 1, nCoupleS
          iProcS = iCoupleProcS_I(iCouple)
          call MPI_Isend(&
               nCouplePointS_I(iCouple), 1, MPI_INTEGER, &
               iProcS, iTag, iComm, iRequest_I(iCouple), iError)
       end do
    end if
    if(nCoupleR > 0)then
       !\
       ! Recv buffer size
       !/
       do iCouple = 1, nCoupleR
          iProcR = iCoupleProcR_I(iCouple)
          call MPI_Irecv(&
               nCouplePointR_I(iCouple), 1, MPI_INTEGER, &
               iProcR, iTag, iComm, iRequest_I(nCoupleS+iCouple), iError)
       end do
    end if
    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I,iStatus_II,iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine get_recv_buffer_size

  !==========================================================================

  subroutine set_data_transfer(iComm, iCompSource, iCompTarget,&
       nBufferSource, iProcSource_I, &
       nBufferTarget, iProcTarget_I, &
       nCoupleSource, iCoupleProcSource_I, nCouplePointSource_I,&
       nCoupleTarget, iCoupleProcTarget_I, nCouplePointTarget_I)

    ! 

    ! subroutine defines transfer based on processors 
    ! of the OTHER components, which own data items
    ! iProc_I

    ! Union communicator
    integer, intent(in):: iComm                        

    ! Source and target component indexes
    integer, intent(in):: iCompSource, iCompTarget

    ! Number of point positions asked from this source proc
    integer, intent(in):: nBufferSource
    ! Index of source processors owning these points
    integer, intent(in):: iProcSource_I(nBufferSource)

    ! Number of points required by the target processor
    integer, intent(in):: nBufferTarget
    ! Index of source processor owning the points of this target proc.
    integer, intent(in):: iProcTarget_I(nBufferTarget)

    ! Number of couplings performed by source on this processor
    integer, intent(out)               :: nCoupleSource
    ! Processor indexes in iComm that source is coupled with
    integer, intent(inout), allocatable:: iCoupleProcSource_I( :)
    ! Number of points per processor source is coupled with
    integer, intent(inout), allocatable:: nCouplePointSource_I(:)

    ! Number of couplings performed by target on this processor
    integer, intent(out)               :: nCoupleTarget
    ! Processor indexes in iComm that target is coupled with
    integer, intent(inout), allocatable:: iCoupleProcTarget_I( :)
    ! Number of points per processor target is coupled with
    integer, intent(inout), allocatable:: nCouplePointTarget_I(:)

    integer:: nProcSource, nProcTarget, nProcIntersect
    integer:: iBuffer, iProc, iCouple, iCoupleBuffer 
    integer:: iProcSourceLocal, iProcTargetLocal
    integer:: iError
    integer, allocatable:: iProcSourceLocal_P(:), iProcTargetLocal_P(:)
    integer, allocatable:: iProcSourceUnion_P(:), iProcTargetUnion_P(:)
    integer, allocatable:: nPoint_PP(:,:), nPointRecv_P(:), nPointSend_P(:) 
    integer, allocatable:: iCoupleProc_P(:), nCouplePoint_P(:)
    character(len=*), parameter:: NameSub = 'set_data_transfer'
    integer:: iTmp
    !----------------------------------------------------------------------
    nProcTarget    = n_proc(iCompTarget)
    nProcSource    = n_proc(iCompSource)

    if(is_proc(iCompTarget))then
       if(allocated(iCoupleProcTarget_I )) deallocate(iCoupleProcTarget_I)
       if(allocated(nCouplePointTarget_I)) deallocate(nCouplePointTarget_I)
       !\
       ! get nCouple and number of points per proc
       !/

       ! Allocate temporary processor index and buffer size arrays.
       ! nProcSource is the maximum possible number of processors to 
       ! communicate with. The actual size will be reduced later.
       allocate(iCoupleProc_P(0:nProcSource-1), &
            nCouplePoint_P(0:nProcSource-1))
       nCouplePoint_P = 0

       ! count points received from each source processors
       do iBuffer = 1, nBufferTarget
          iProc = iProcTarget_I(iBuffer)
          if(iProc < 0) CYCLE
          nCouplePoint_P(iProc) = nCouplePoint_P(iProc) + 1
       end do

       ! Skip source processors that have no data to send.
       ! Count the source processors that have data (nCoupleTarget).
       ! Compact the arrays to the source processors that have data.
       nCoupleTarget = 0
       do iProc = 0, nProcSource-1
          if(nCouplePoint_P(iProc) == 0) CYCLE
          ! One more source proc with data
          nCoupleTarget = nCoupleTarget + 1
          ! Move up the proc index and number of points 
          ! from iProc to nCoupleTarget-1
          iCoupleProc_P( nCoupleTarget-1) = iProc
          nCouplePoint_P(nCoupleTarget-1) = nCouplePoint_P(iProc)
       end do

       ! Copy compacted arrays into the output arguments
       allocate(&
            iCoupleProcTarget_I( nCoupleTarget), &
            nCouplePointTarget_I(nCoupleTarget))

       iCoupleProcTarget_I  = iCoupleProc_P(0:nCoupleTarget-1)
       nCouplePointTarget_I = nCouplePoint_P(0:nCoupleTarget-1)
       deallocate(iCoupleProc_P, nCouplePoint_P)
    end if

    if(is_proc(iCompSource))then
       allocate(iProcSourceUnion_P(0:nProcSource-1))
       allocate(iProcTargetUnion_P(0:nProcTarget-1))
       allocate(iProcSourceLocal_P(0:nProcSource+nProcTarget-1))
       allocate(iProcTargetLocal_P(0:nProcSource+nProcTarget-1))
       call get_rank_translation(&
            iComm, iCompSource, iCompTarget, nProcSource, nProcTarget, &
            iProcSourceLocal_P, iProcTargetLocal_P,                    &
            iProcSourceUnion_P, iProcTargetUnion_P)

       ! nPoint_PP is number of points to sent from Source to Target
       allocate(nPoint_PP(0:nProcSource-1, 0:nProcTarget-1))

       ! build nPoint_PP
       ! each proc on Source has only part of the info, store it into nPoint_PP
       nPoint_PP = 0
       if(nCoupleSource > 0)then
          iCouple = 1
          iProcTargetLocal = iProcTargetLocal_P(iCoupleProcSource_I(iCouple))
          iCoupleBuffer = nCouplePointSource_I(iCouple)
          do iBuffer = 1, nBufferSource
             iProc = iProcSource_I(iBuffer)
             if(iProc < 0) CYCLE
             if(iBuffer > iCoupleBuffer)then
                iCouple = iCouple + 1
                iCoupleBuffer = iCoupleBuffer + nCouplePointSource_I(iCouple)
                ! iCoupleProcSource_I(iCouple) is the index of target proc for
                ! union communicator that sent the request to this source 
                ! processor. iProcTargetLocal is the index of taget processor 
                ! on the target component
                iProcTargetLocal = iProcTargetLocal_P(&
                     iCoupleProcSource_I(iCouple))
             end if
             ! iProcSource_I is the source processor index on source component
             iProcSourceLocal = iProcSourceLocal_P(iProc)
             
             ! Count this point in the communication matrix
             nPoint_PP(iProcSourceLocal,iProcTargetLocal) = &
                  nPoint_PP(iProcSourceLocal,iProcTargetLocal) + 1
          end do
       end if

       deallocate(iCoupleProcSource_I, nCouplePointSource_I)

       !\
       ! Reduce the row nPoint_PP(iProcSourceLocal,:)
       !/
       allocate(nPointSend_P(0:nProcTarget-1), nPointRecv_P(0:nProcTarget-1))

       ! Each source processor collects the row of the matrix that it needs.
       ! Each row is added up from all source processors, because they
       ! were asked about point positions randomly.
       do iProcSourceLocal = 0, nProcSource-1
          ! Extract the particular row (partial info)
          nPointSend_P = nPoint_PP(iProcSourceLocal,:)
          ! Add up partial rows and collect it on iProcSourceLocal
          call MPI_Reduce(&
               nPointSend_P, nPointRecv_P, nProcTarget, MPI_INTEGER, MPI_SUM, &
               iProcSourceLocal, i_comm(iCompSource), iError)
       end do

       deallocate(nPoint_PP, nPointSend_P)

       ! Skip target processors that have no data to receive.
       ! Count the target processors that have data (nCoupleSource).
       ! Compact the arrays to the target processors that have data.
       ! 
       allocate(iCoupleProc_P(0:nProcTarget-1), &
            nCouplePoint_P(0:nProcTarget-1))
       nCoupleSource = 0
       do iProcTargetLocal = 0, nProcTarget-1
          if(nPointRecv_P(iProcTargetLocal) == 0) CYCLE
          nCoupleSource = nCoupleSource + 1

          ! Translate local rank iProcTargetLocal to union rank
          iCoupleProc_P(nCOupleSource-1) = iProcTargetUnion_P(iProcTargetLocal)
          nCouplePoint_P(nCoupleSource-1) = nPointRecv_P(iProcTargetLocal)
       end do

       ! Copy compacted arrays into the output arguments
       allocate(iCoupleProcSource_I( nCoupleSource))
       allocate(nCouplePointSource_I(nCoupleSource))
       iCoupleProcSource_I  = iCoupleProc_P(0:nCoupleSource-1)
       nCouplePointSource_I = nCouplePoint_P(0:nCoupleSource-1)

       deallocate(iCoupleProc_P, nCouplePoint_P, nPointRecv_P)
       deallocate(iProcSourceLocal_P, iProcSourceUnion_P)
       deallocate(iProcTargetLocal_P, iProcTargetUnion_P)
    end if

  end subroutine set_data_transfer

  !==========================================================================

  subroutine get_transfer_buffer_order(&
       iCompSource, iCompTarget, &
       nBuffer, iProc_I, &
       nCouple, iCoupleProc_I, nCouplePoint_I,&
       iBuffer_I, nData)
    ! iBuffer_I returns the order of the points with processor index
    ! corresponding to iCoupleProc_I
    ! nData is the number of points found on Source
    integer, intent(in) :: iCompSource, iCompTarget ! send, recv components
    integer, intent(in) :: nBuffer    ! number of points
    integer, intent(in) :: iProc_I(nBuffer)  ! proc index for each point
    integer, intent(in) :: nCouple
    integer, intent(in) :: iCoupleProc_I(nCouple), nCouplePoint_I(nCouple)
    integer, intent(out):: iBuffer_I(nBuffer)! index for reordering
    integer, intent(out):: nData ! number of points found on Source

    integer:: iBuffer, iCouple, iProc, nProc, nPointNotFound
    integer, allocatable:: iOrder_I(:)
    integer, allocatable:: iCoupleProcInv_P(:)
    character(len=*), parameter:: NameSub = 'get_transfer_buffer_order'
    !----------------------------------------------------------------------
    nProc = n_proc(iCompSource)
    nData = nBuffer
    ! Set starting point of chunks in the reordered buffer
    allocate(iOrder_I(nCouple))
    iOrder_I(1) = 1
    do iCouple = 2, nCouple
       iOrder_I(iCouple) = iOrder_I(iCouple-1) + nCouplePoint_I(iCouple-1)
    end do

    ! inversion of iCoupleProc_I, i.e. iCouple corresponding to iProc
    allocate(iCoupleProcInv_P(0:nProc-1))
    iCoupleProcInv_P = -1
    do iCouple = 1, nCouple
       iCoupleProcInv_P(iCoupleProc_I(iCouple)) = iCouple
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
          ! It will occupy iBuffer_I position
          iBuffer_I(iBuffer) = iOrder_I(iCoupleProcInv_P(iProc))
          ! Jump to next position
          iOrder_I(iCoupleProcInv_P(iProc)) = &
               iOrder_I(iCoupleProcInv_P(iProc)) + 1
       end if
    end do

    deallocate(iCoupleProcInv_P, iOrder_I)

  end subroutine get_transfer_buffer_order

  !===========================================================================

  subroutine transfer_buffer_int(&
       iComm, nProc, iProc,&
       nCoupleS, iCoupleProcS_I, nCouplePointS_I,&
       nCoupleR, iCoupleProcR_I, nCouplePointR_I,&
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
    integer, intent(in):: iCoupleProcS_I( nCoupleS)! ids of procs to send to
    integer, intent(in):: nCouplePointS_I(nCoupleS)! # of points  to send

    integer, intent(in):: nCoupleR                 ! #   of procs to recv from
    integer, intent(in):: iCoupleProcR_I( nCoupleR)! ids of procs to recv from
    integer, intent(in):: nCouplePointR_I(nCoupleR)! # of points  to recv

    integer, intent(in   ):: nData                 ! number of items per point
    integer, intent(in   ):: nBufferS                   ! send buffer size
    integer, intent(inout):: iBufferS_I(nData*nBufferS) ! send buffer
    integer, intent(in)   :: nBufferR                   ! recv buffer size
    integer, intent(inout):: iBufferR_I(nData*nBufferR) ! recv buffer

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
       ! Send array of integers
       !/
       iBuffer = 1
       do iCouple = 1, nCoupleS
          iProcS = iCoupleProcS_I(iCouple)
          call MPI_Isend(&
               iBufferS_I(iBuffer), nData*nCouplePointS_I(iCouple), &
               MPI_INTEGER, iProcS, iTag, iComm, &
               iRequest_I(iCouple),iError)       
          iBuffer = iBuffer + nData*nCouplePointS_I(iCouple)
       end do
    end if
    if(nCoupleR > 0)then
       !\
       ! Receive array of integers
       !/
       iBuffer = 1
       do iCouple = 1, nCoupleR
          iProcR = iCoupleProcR_I(iCouple)
          call MPI_Irecv(&
               iBufferR_I(iBuffer), nData*nCouplePointR_I(iCouple),&
               MPI_INTEGER, iProcR, iTag, iComm, &
               iRequest_I(nCoupleS+iCouple),iError)
          iBuffer = iBuffer + nData*nCouplePointR_I(iCouple)
       end do
    end if

    !\
    ! Finalize transfer
    !/
    call MPI_waitall(nCoupleS + nCoupleR, iRequest_I, iStatus_II, iError)
    deallocate(iRequest_I, iStatus_II)
  end subroutine transfer_buffer_int

  !===========================================================================

  subroutine get_rank_translation(                                      &
       iCommUnion, iCompSource, iCompTarget, nProcSource, nProcTarget,  &
       iProcSourceLocal_I, iProcTargetLocal_I,                          &
       iProcSourceUnion_I, iProcTargetUnion_I)
    ! for processor with rank iProc in the union communicator iCommUnion
    ! returns rank iProcLocal in the group associated with component iComp
    integer, intent(in) :: iCommUnion, iCompSource, iCompTarget
    integer, intent(in) :: nProcSource, nProcTarget
    integer, intent(out):: iProcSourceLocal_I(nProcSource+nProcTarget)
    integer, intent(out):: iProcTargetLocal_I(nProcSource+nProcTarget)
    integer, intent(out):: iProcSourceUnion_I(nProcSource)
    integer, intent(out):: iProcTargetUnion_I(nProcTarget)

    integer:: iGroupUnion
    integer:: iProc
    integer, allocatable:: iProc_I(:)
    integer:: iError
    !----------------------------------------------------------------------
    ! Create group associated with iCommUnion
    call MPI_Comm_Group(iCommUnion, iGroupUnion, iError)

    ! Translate ranks for Source Local to Union
    allocate(iProc_I(nProcSource))
    do iProc = 0, nProcSource-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         i_group(iCompSource), nProcSource, iProc_I, &
         iGroupUnion, iProcSourceUnion_I, iError)
    deallocate(iProc_I)

    ! Translate ranks for Target Local to Union
    allocate(iProc_I(nProcTarget))
    do iProc = 0, nProcTarget-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         i_group(iCompTarget), nProcTarget, iProc_I, &
         iGroupUnion, iProcTargetUnion_I, iError)
    deallocate(iProc_I)

    ! Translate ranks Union to Local
    allocate(iProc_I(nProcSource+nProcTarget))
    do iProc = 0, nProcSource+nProcTarget-1
       iProc_I(iProc+1) = iProc
    end do
    call MPI_Group_translate_ranks(&
         iGroupUnion, nProcSource+nProcTarget, iProc_I, &
         i_group(iCompSource), iProcSourceLocal_I, iError)
    call MPI_Group_translate_ranks(&
         iGroupUnion, nProcSource+nProcTarget, iProc_I, &
         i_group(iCompTarget), iProcTargetLocal_I, iError)
    deallocate(iProc_I)

    ! free union group
    call MPI_Group_free(iGroupUnion, iError)
  end subroutine get_rank_translation

end module CON_couple_gm_pt
