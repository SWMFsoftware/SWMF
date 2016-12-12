!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_pass_face_field
  use BATL_size, ONLY: MaxBlock, nBlock, nI, nJ, nK, nIjk_D, &
       MaxDim, nDim, jDim_, kDim_

  use BATL_mpi, ONLY: iComm, nProc, iProc, barrier_mpi
  
  use BATL_tree, ONLY: &
       iNodeNei_IIIB, DiLevelNei_IIIB, Unset_, iNode_B, &
       iTree_IA, Proc_, Block_, Unused_B
  implicit none
  
  SAVE
  
  private ! except
  !\
  ! PUBLIC MEMBERS
  !/
  !\
  ! Fills in ghost face values for the face-centered staggered
  ! (electric) field (in Particle-In-Cell code)
  !/
  public message_pass_field
  !\
  ! Add ghost face values of face-centered staggered vector 
  ! variables (such as the electric current components in PIC, which
  ! are created by the particles close to the block boundary)
  ! to all physical faces 
  !/
  public add_ghost_face_field
  ! Fast lookup tables for index ranges per dimension
  integer, parameter:: Min_=1, Max_=2
  integer:: iS_DIID(MaxDim,-1:1,Min_:Max_,MaxDim)
  integer:: iR_DIID(MaxDim,-1:1,Min_:Max_,MaxDim)
  
  ! Variables related to recv and send buffers
  integer, allocatable:: iBufferS_P(:), nBufferS_P(:), nBufferR_P(:)
  
  real,    allocatable:: BufferR_I(:), BufferS_I(:)

  !\
  ! Counter of currents through the physical faces
  !/
  real, allocatable :: Counter_FDB(:,:,:,:,:) 
  integer:: nBlockMax = -1
  
  integer, allocatable:: iRequestR_I(:), iRequestS_I(:), &
       iStatus_II(:,:)
contains
  subroutine message_pass_field(nG, Field_FDB, nWidthIn)
    use ModMpi

    ! Arguments
    integer, intent(in) :: nG    ! number of ghost cells for 1..nDim
    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1. By x the physical cells are marked below
    !  Left corner        Right corner n+2                 Ez (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_n+2     Not used: Ex (nX+2,:,:)
    !       _!_!x!x               _!_!_!_                  Ey (:,nY+2,:)
    !       _!_!_!                x!x!_!_        Phys Face Values: 
    !    -2  ! ! !                x!x!_!_        Ex(0,1,1),Ex(nI,1,1)...
    !       -2                             Ghost values: Ex(-1,1,1),Ey(0,1,1)..
    real, intent(inout) :: Field_FDB(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
         1-nG*kDim_:nK+nG*kDim_,MaxDim,MaxBlock)

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    
    ! nWidthIn is the number of ghost face layers to be set. Default is nG.


    ! Local variables

    logical, parameter :: UseRSend = .false.

    ! local variables corresponding to optional arguments
    integer :: nWidth

    ! Various indexes 
    integer :: iCountOnly     ! index for 2 stage scheme for count, sendrecv
    logical :: DoCountOnly    ! logical for count vs. sendrecv stages

    integer :: iDir, jDir, kDir
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

    integer :: iBufferS, iBufferR
    integer :: MaxBufferS = -1, MaxBufferR = -1
    integer:: iRequestR, iRequestS, iError
 
    character(len=*), parameter:: NameSub = 'BATL_pass_face_field'
    !--------------------------------------------------------------------------


    call timing_start('batl_pass')

    call timing_start('init_pass')

    ! Set values or defaults for optional arguments
    nWidth = nG
    if(present(nWidthIn)) nWidth = nWidthIn


    if(nWidth < 1 .or. nWidth > nG) call CON_stop(NameSub// &
         ' nWidth do not contain the ghost cells or too many')
 
    ! Set index ranges based on arguments
    call set_range

    if(nProc > 1)then
       ! Allocate fixed size communication arrays.
       if(.not.allocated(iBufferS_P))then
          allocate(iBufferS_P(0:nProc-1), nBufferS_P(0:nProc-1), &
               nBufferR_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
    end if


    call timing_stop('init_pass')


    do iCountOnly = 1, 2
       DoCountOnly = iCountOnly == 1

       ! No need to count data for send/recv in serial runs
       if(nProc == 1 .and. DoCountOnly) CYCLE

       call timing_start('local_pass')

       if(nProc>1)then
          if(DoCountOnly)then
             ! initialize buffer size
             nBufferR_P = 0
             nBufferS_P = 0
          else
             ! Make buffers large enough
             if(sum(nBufferR_P) > MaxBufferR) then
                if(allocated(BufferR_I)) deallocate(BufferR_I)
                MaxBufferR = sum(nBufferR_P) 
                allocate(BufferR_I(MaxBufferR))
             end if

             if(sum(nBufferS_P) > MaxBufferS) then
                if(allocated(BufferS_I)) deallocate(BufferS_I)
                MaxBufferS = sum(nBufferS_P)
                allocate(BufferS_I(MaxBufferS))
             end if

             ! Initialize buffer indexes for storing data into BufferS_I
             iBufferS = 0
             do iProcRecv = 0, nProc-1
                iBufferS_P(iProcRecv) = iBufferS
                iBufferS = iBufferS + nBufferS_P(iProcRecv)
             end do
          end if
       end if
       ! Loop through all blocks that may send a message
       do iBlockSend = 1, nBlock
          if(Unused_B(iBlockSend)) CYCLE
          do kDir = -1, 1
             ! Do not message pass in ignored dimensions
             if(nDim < 3 .and. kDir /= 0) CYCLE
             do jDir = -1, 1
                if(nDim < 2 .and. jDir /= 0) CYCLE
                do iDir = -1,1
                   ! Ignore inner parts of the sending block
                   if(iDir == 0 .and. jDir == 0 .and. kDir == 0) CYCLE
                   if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)/= Unset_)&
                        call do_equal
                end do ! iDir
             end do ! jDir
          end do ! kDir
       end do ! iBlockSend

       call timing_stop('local_pass')

    end do ! iCountOnly

    if(nProc==1)then
       call timing_stop('batl_pass')
       RETURN
    end if

    call timing_start('recv_pass')

    ! post requests
    iRequestR = 0
    iBufferR  = 1
    do iProcSend = 0, nProc - 1
       if(nBufferR_P(iProcSend) == 0) CYCLE
       iRequestR = iRequestR + 1

       call MPI_irecv(BufferR_I(iBufferR), nBufferR_P(iProcSend), &
            MPI_REAL, iProcSend, 10, iComm, iRequestR_I(iRequestR), &
            iError)

       iBufferR  = iBufferR  + nBufferR_P(iProcSend)
    end do

    call timing_stop('recv_pass')

    if(UseRSend) then
       call timing_start('barrier_pass')
       call barrier_mpi
       call timing_stop('barrier_pass')
    end if

    call timing_start('send_pass')

    ! post sends
    iRequestS = 0
    iBufferS  = 1
    do iProcRecv = 0, nProc-1
       if(nBufferS_P(iProcRecv) == 0) CYCLE
       iRequestS = iRequestS + 1

       if(UseRSend)then
          call MPI_rsend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iError)
       else
          call MPI_isend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iRequestS_I(iRequestS), &
               iError)
       end if

       iBufferS  = iBufferS  + nBufferS_P(iProcRecv)
    end do
    call timing_stop('send_pass')

    call timing_start('wait_pass')
    ! wait for all requests to be completed
    if(iRequestR > 0) &
         call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

    ! wait for all sends to be completed
    if(.not.UseRSend .and. iRequestS > 0) &
         call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)
    call timing_stop('wait_pass')

    call timing_start('buffer_to_state')
    call buffer_to_state
    call timing_stop('buffer_to_state')


    call timing_stop('batl_pass')

  contains

    !======================================================================
    subroutine buffer_to_state

      ! Copy buffer into recv block of Field_FDB

      integer:: iBufferR, i, j, k, iDim
      !------------------------------------------------------------------------
      jRMin = 1; jRMax = 1
      kRMin = 1; kRMax = 1

      iBufferR = 0
      do iProcSend = 0, nProc - 1
         if(nBufferR_P(iProcSend) == 0) CYCLE

         do
            iBufferR = iBufferR + 1
            iBlockRecv = nint(BufferR_I(iBufferR))
            do iDim = 1, MaxDim
               iRMin      = nint(BufferR_I(iBufferR+1))
               iRMax      = nint(BufferR_I(iBufferR+2))
               if(nDim > 1) jRMin = nint(BufferR_I(iBufferR+3))
               if(nDim > 1) jRMax = nint(BufferR_I(iBufferR+4))
               if(nDim > 2) kRMin = nint(BufferR_I(iBufferR+5))
               if(nDim > 2) kRMax = nint(BufferR_I(iBufferR+6))
               
               iBufferR = iBufferR + 2*nDim
               
               
               do k=kRMin,kRmax; do j=jRMin,jRMax; do i=iRMin,iRmax 
                  iBufferR = iBufferR + 1
                  Field_FDB(i,j,k,iDim,iBlockRecv) = &
                       BufferR_I(iBufferR)
               end do; end do; end do
            end do
            if(iBufferR >= sum(nBufferR_P(0:iProcSend))) EXIT
         end do
       end do
    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: nSize, iDim, iBufferS, i, j, k, iSend, jSend, kSend, iNodeRecv
      !------------------------------------------------------------------------

      iSend = (3*iDir + 3)/2
      jSend = (3*jDir + 3)/2
      kSend = (3*kDir + 3)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      
      if(DoCountOnly)then
         ! No need to count data for local copy
         if(iProc == iProcRecv) RETURN
         !\
         ! Number of integers to be send/received
         !/
         nSize = 1 + 2*nDim*MaxDim
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)

 
            ! Number of reals to send to and received from the other processor
            nSize = (iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) + nSize
         end do
         nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize
         nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
         RETURN
      end if
      if(iProc == iProcRecv)then
            ! Local copy
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)
            iSMin = iS_DIID(1,iDir,Min_,iDim)
            iSMax = iS_DIID(1,iDir,Max_,iDim)
            jSMin = iS_DIID(2,jDir,Min_,iDim)
            jSMax = iS_DIID(2,jDir,Max_,iDim)
            kSMin = iS_DIID(3,kDir,Min_,iDim)
            kSMax = iS_DIID(3,kDir,Max_,iDim)
            if((iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)==0)CYCLE
         
            Field_FDB(iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iDim, &
                 iBlockRecv)= Field_FDB(iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,&
                 iDim,iBlockSend)
         end do
      else
            ! Put data into the send buffer
            iBufferS = iBufferS_P(iProcRecv) +1

            BufferS_I(            iBufferS) = iBlockRecv
            do iDim = 1,MaxDim
               iRMin = iR_DIID(1,iDir,Min_,iDim)
               iRMax = iR_DIID(1,iDir,Max_,iDim)
               jRMin = iR_DIID(2,jDir,Min_,iDim)
               jRMax = iR_DIID(2,jDir,Max_,iDim)
               kRMin = iR_DIID(3,kDir,Min_,iDim)
               kRMax = iR_DIID(3,kDir,Max_,iDim)
               iSMin = iS_DIID(1,iDir,Min_,iDim)
               iSMax = iS_DIID(1,iDir,Max_,iDim)
               jSMin = iS_DIID(2,jDir,Min_,iDim)
               jSMax = iS_DIID(2,jDir,Max_,iDim)
               kSMin = iS_DIID(3,kDir,Min_,iDim)
               kSMax = iS_DIID(3,kDir,Max_,iDim)
               BufferS_I(            iBufferS+1) = iRMin
               BufferS_I(            iBufferS+2) = iRMax
               if(nDim > 1)BufferS_I(iBufferS+3) = jRMin
               if(nDim > 1)BufferS_I(iBufferS+4) = jRMax
               if(nDim > 2)BufferS_I(iBufferS+5) = kRMin
               if(nDim > 2)BufferS_I(iBufferS+6) = kRMax

               iBufferS = iBufferS + 2*nDim

               do k = kSMin,kSmax; do j = jSMin,jSMax; do i = iSMin,iSmax
                  iBufferS = iBufferS + 1
                  BufferS_I(iBufferS) = Field_FDB(i,j,k,iDim,iBlockSend)
               end do; end do; end do
            end do
            iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_equal

    !==========================================================================

    subroutine set_range

      integer:: iDim
      !Array index is the coordinate of the gridpoint with +1/2 being
      !approximated as 1. By x the physical cells are marked below
      !  Left corner        Right corner n+2                 Ez (:,:,nZ+2)    
      !       _!_!x!x               _!_!_!_n+2     Not used: Ex (nX+2,:,:)
      !       _!_!x!x               _!_!_!_                  Ey (:,nY+2,:)
      !       _!_!_!                x!x!_!_        Phys Face Values: 
      !    -2  ! ! !                x!x!_!_        Ex(0,1,1),Ex(nI,1,1)...
      !       -2                           Ghost values: Ex(-1,1,1),Ey(0,1,1)..
      !------------------------------------------------------------------------
      !\
      ! For x_, y_, z_ direction (the first index), the receiving block is at
      ! negative displcement with respect to the sending one
      !/
      !\
      ! 1. Send faces 1:nWidth
      !/
      iS_DIID(:,-1,Min_,:) = 1
      iS_DIID(:,-1,Max_,:) = nWidth
      do iDim = 1,MaxDim
         !\
         ! 2. Receive faces nIjk_D + 1:nIjk_D + nWidth
         !/
         iR_DIID(:,-1,Min_,iDim) = nIjk_D + 1
         iR_DIID(:,-1,Max_,iDim) = nIjk_D + nWidth
         !\
         ! !EXCEPT!!
         ! nIJK+nWidth face at which the iDim component is assigned
         ! is not used. Specifiaclly for nWidth=1 this face is not sent. 
         iS_DIID(iDim,-1,Max_,iDim) = iS_DIID(iDim,-1,Max_,iDim) - 1
         iR_DIID(iDim,-1,Max_,iDim) = iR_DIID(iDim,-1,Max_,iDim) - 1
      end do
      !\
      ! For x_, y_, z_ direction (the first index), the receiving block is at
      ! zero displcement with respect to the sending one
      !/
      !\
      ! 1. Send and receive faces are from 1 to nIJK_D
      !/
      iS_DIID(:, 0,Min_,:) = 1
      iR_DIID(:, 0,Min_,:) = 1
      do iDim = 1,MaxDim
         iS_DIID(:, 0,Max_,iDim) = nIjk_D
         iR_DIID(:, 0,Max_,iDim) = nIjk_D
         !\
         ! !EXCEPT!!! For the iDim component of the field
         ! Send and receive faces are from 0 to nIJK_D
         !/
         iS_DIID(iDim,0,Min_,iDim) = 0
         iR_DIID(iDim,0,Min_,iDim) = 0
      end do
      !\
      ! For x_, y_, z_ direction (the first index), the receiving block is at
      ! positive displcement with respect to the sending one
      !/
      !\
      ! 1. Recv  faces are from 1-nWidth to 0
      !/
      iR_DIID(:, 1,Min_,:) = 1 - nWidth
      iR_DIID(:, 1,Max_,:) = 0
      do iDim = 1,MaxDim
         !\
         ! Send faces are from nIjk_D + 1 - nWidth to nIjk_D
         !/
         iS_DIID(:, 1,Min_,iDim) = nIjk_D + 1 - nWidth
         iS_DIID(:, 1,Max_,iDim) = nIjk_D
         !\
         ! !EXCEPT!!! For the iDim component of the field
         ! Recv faces are up to -1 (0 face is physical)
         !/
         iS_DIID(iDim,1,Max_,iDim) = nIjk_D(iDim) - 1
         iR_DIID(iDim,1,Max_,iDim) =  - 1
      end do


    end subroutine set_range

  end subroutine message_pass_field
  !===========================
  subroutine add_ghost_face_field(nG, Current_FDB, nWidthIn)
    use ModMpi

    ! Arguments
    integer, intent(in) :: nG    ! number of ghost cells for 1..nDim
    real, intent(inout), dimension(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
         1-nG*kDim_:nK+nG*kDim_,MaxDim,MaxBlock):: Current_FDB
    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1. By x the physical cells are marked below
    !  Left corner        Right corner n+2                 Jz (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_n+2     Not used: Jx (nX+2,:,:)
    !       _!_!x!x               _!_!_!_                  Jy (:,nY+2,:)
    !       _!_!_!                x!x!_!_        Phys Face Values: 
    !    -2  ! ! !                x!x!_!_        Jx(0,1,1),Jx(nI,1,1)...
    !       -2                           Ghost values: Jx(-1,1,1),Jy(0,1,1)..
    !------------------------------------------------------------------------
    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    
    ! nWidthIn is the number of ghost face layers to be set. Default is all.


    ! Local variables

    logical, parameter :: UseRSend = .false.

    ! local variables corresponding to optional arguments
    integer :: nWidth

    ! Various indexes 
    integer :: iCountOnly     ! index for 2 stage scheme for count, sendrecv
    logical :: DoCountOnly    ! logical for count vs. sendrecv stages

    integer :: iDir, jDir, kDir
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

    integer :: iBufferS, iBufferR
    integer :: MaxBufferS = -1, MaxBufferR = -1
    integer:: iRequestR, iRequestS, iError
 
    character(len=*), parameter:: NameSub = 'BATL_pass::add_ghost_face_field'
    !--------------------------------------------------------------------------
    call timing_start('batl_pass')

    call timing_start('init_pass')
    !\
    ! Allocate block-based counter for currents through hysical phases
    !/
    if(nBlock>nBlockMax)then
       if(allocated(Counter_FDB))deallocate(Counter_FDB)
       allocate(Counter_FDB(0:nI, 1-jDim_:nJ, 1-kDim_:nK, MaxDim, nBlock))
       nBlockMax = nBlock
    end if
    !\
    ! Nullify counter
    !/
    Counter_FDB = 0.0
    !\
    ! Fill in physical faces of the Counter:
    !/
    do iBlockSend = 1, nBlock
       if(Unused_B(iBlockSend)) CYCLE
       Counter_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend) = &
            Current_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend)
       Counter_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend) = &
            Current_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend)
       Counter_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend) = &
            Current_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend)
    end do
    ! Set values or defaults for optional arguments
    nWidth = nG
    if(present(nWidthIn)) nWidth = nWidthIn

    if(nWidth < 1 .or. nWidth > nG) call CON_stop(NameSub// &
         ' nWidth do not contain the ghost cells or too many')

    ! Set index ranges based on arguments
    call set_range

    if(nProc > 1)then
       ! Allocate fixed size communication arrays.
       if(.not.allocated(iBufferS_P))then
          allocate(iBufferS_P(0:nProc-1), nBufferS_P(0:nProc-1), &
               nBufferR_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
    end if


    call timing_stop('init_pass')


    do iCountOnly = 1, 2
       DoCountOnly = iCountOnly == 1

       ! No need to count data for send/recv in serial runs
       if(nProc == 1 .and. DoCountOnly) CYCLE

       call timing_start('local_pass')

       ! Second order prolongation needs two stages: 
       ! first stage fills in equal and coarser ghost cells
       ! second stage uses these to prolong and fill in finer ghost cells

       if(nProc>1)then
          if(DoCountOnly)then
             ! initialize buffer size
             nBufferR_P = 0
             nBufferS_P = 0
          else
             ! Make buffers large enough
             if(sum(nBufferR_P) > MaxBufferR) then
                if(allocated(BufferR_I)) deallocate(BufferR_I)
                MaxBufferR = sum(nBufferR_P) 
                allocate(BufferR_I(MaxBufferR))
             end if

             if(sum(nBufferS_P) > MaxBufferS) then
                if(allocated(BufferS_I)) deallocate(BufferS_I)
                MaxBufferS = sum(nBufferS_P)
                allocate(BufferS_I(MaxBufferS))
             end if

             ! Initialize buffer indexes for storing data into BufferS_I
             iBufferS = 0
             do iProcRecv = 0, nProc-1
                iBufferS_P(iProcRecv) = iBufferS
                iBufferS = iBufferS + nBufferS_P(iProcRecv)
             end do
          end if
       end if
       ! Loop through all blocks that may send a message
       do iBlockSend = 1, nBlock
          if(Unused_B(iBlockSend)) CYCLE
          do kDir = -1, 1
             ! Do not message pass in ignored dimensions
             if(nDim < 3 .and. kDir /= 0) CYCLE
             do jDir = -1, 1
                if(nDim < 2 .and. jDir /= 0) CYCLE
                do iDir = -1,1
                   ! Ignore inner parts of the sending block
                   if(iDir == 0 .and. jDir == 0 .and. kDir == 0) CYCLE
                   if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)/= Unset_)&
                        call do_equal
                end do ! iDir
             end do ! jDir
          end do ! kDir
       end do ! iBlockSend

       call timing_stop('local_pass')

    end do ! iCountOnly

    if(nProc==1)then
       do iBlockSend = 1, nBlock
          if(Unused_B(iBlockSend)) CYCLE
          Current_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend) = &
               Counter_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend)
          Current_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend) = &
               Counter_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend)
          Current_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend) = &
               Counter_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend)
       end do
       call timing_stop('batl_pass')
       RETURN
    end if

    call timing_start('recv_pass')

    ! post requests
    iRequestR = 0
    iBufferR  = 1
    do iProcSend = 0, nProc - 1
       if(nBufferR_P(iProcSend) == 0) CYCLE
       iRequestR = iRequestR + 1

       call MPI_irecv(BufferR_I(iBufferR), nBufferR_P(iProcSend), &
            MPI_REAL, iProcSend, 10, iComm, iRequestR_I(iRequestR), &
            iError)

       iBufferR  = iBufferR  + nBufferR_P(iProcSend)
    end do

    call timing_stop('recv_pass')

    if(UseRSend) then
       call timing_start('barrier_pass')
       call barrier_mpi
       call timing_stop('barrier_pass')
    end if

    call timing_start('send_pass')

    ! post sends
    iRequestS = 0
    iBufferS  = 1
    do iProcRecv = 0, nProc-1
       if(nBufferS_P(iProcRecv) == 0) CYCLE
       iRequestS = iRequestS + 1

       if(UseRSend)then
          call MPI_rsend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iError)
       else
          call MPI_isend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iRequestS_I(iRequestS), &
               iError)
       end if

       iBufferS  = iBufferS  + nBufferS_P(iProcRecv)
    end do
    call timing_stop('send_pass')

    call timing_start('wait_pass')
    ! wait for all requests to be completed
    if(iRequestR > 0) &
         call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

    ! wait for all sends to be completed
    if(.not.UseRSend .and. iRequestS > 0) &
         call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)
    call timing_stop('wait_pass')

    call timing_start('buffer_to_state')
    call buffer_to_state
    call timing_stop('buffer_to_state')
    do iBlockSend = 1, nBlock
       if(Unused_B(iBlockSend)) CYCLE
       Current_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend) = &
            Counter_FDB(0:nI, 1:nJ, 1:nK, 1, iBlockSend)
       Current_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend) = &
            Counter_FDB(1:nI, 1-jDim_:nJ, 1:nK, 2, iBlockSend)
       Current_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend) = &
            Counter_FDB(1:nI, 1:nJ, 1-kDim_:nK, 3, iBlockSend)
    end do

    call timing_stop('batl_pass')

  contains

    !======================================================================
    subroutine buffer_to_state

      ! Copy buffer into recv block of Current_FDB

      integer:: iBufferR, i, j, k, iDim
      !------------------------------------------------------------------------
      jRMin = 1; jRMax = 1
      kRMin = 1; kRMax = 1

      iBufferR = 0
      do iProcSend = 0, nProc - 1
         if(nBufferR_P(iProcSend) == 0) CYCLE

         do
            iBufferR = iBufferR + 1
            iBlockRecv = nint(BufferR_I(iBufferR))
            do iDim = 1, MaxDim
               iRMin      = nint(BufferR_I(iBufferR+1))
               iRMax      = nint(BufferR_I(iBufferR+2))
               if(nDim > 1) jRMin = nint(BufferR_I(iBufferR+3))
               if(nDim > 1) jRMax = nint(BufferR_I(iBufferR+4))
               if(nDim > 2) kRMin = nint(BufferR_I(iBufferR+5))
               if(nDim > 2) kRMax = nint(BufferR_I(iBufferR+6))
               
               iBufferR = iBufferR + 2*nDim
               
               
               do k=kRMin,kRmax; do j=jRMin,jRMax; do i=iRMin,iRmax 
                  iBufferR = iBufferR + 1
                  Counter_FDB(i,j,k,iDim,iBlockRecv) =      &
                       Counter_FDB(i,j,k,iDim,iBlockRecv) + &
                       BufferR_I(iBufferR)
               end do; end do; end do
            end do
            if(iBufferR >= sum(nBufferR_P(0:iProcSend))) EXIT
         end do
       end do
    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: nSize, iDim, iBufferS, i, j, k, iSend, jSend, kSend, iNodeRecv
      !------------------------------------------------------------------------

      iSend = (3*iDir + 3)/2
      jSend = (3*jDir + 3)/2
      kSend = (3*kDir + 3)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      
      if(DoCountOnly)then
         ! No need to count data for local copy
         if(iProc == iProcRecv) RETURN
         !\
         ! Number of integers to be send/received
         !/
         nSize = 1 + 2*nDim*MaxDim
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)

 
            ! Number of reals to send to and received from the other processor
            nSize = (iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) + nSize
         end do
         nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize
         nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
         RETURN
      end if
      if(iProc == iProcRecv)then
         ! Local copy
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)
            iSMin = iS_DIID(1,iDir,Min_,iDim)
            iSMax = iS_DIID(1,iDir,Max_,iDim)
            jSMin = iS_DIID(2,jDir,Min_,iDim)
            jSMax = iS_DIID(2,jDir,Max_,iDim)
            kSMin = iS_DIID(3,kDir,Min_,iDim)
            kSMax = iS_DIID(3,kDir,Max_,iDim)
            if((iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)==0)CYCLE
         
            Counter_FDB(iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iDim,iBlockRecv) =    &
                 Counter_FDB(iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iDim,iBlockRecv) &
                 + Current_FDB(iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iDim,iBlockSend)
         end do
      else
            ! Put data into the send buffer
            iBufferS = iBufferS_P(iProcRecv) +1

            BufferS_I(            iBufferS) = iBlockRecv
            do iDim = 1,MaxDim
               iRMin = iR_DIID(1,iDir,Min_,iDim)
               iRMax = iR_DIID(1,iDir,Max_,iDim)
               jRMin = iR_DIID(2,jDir,Min_,iDim)
               jRMax = iR_DIID(2,jDir,Max_,iDim)
               kRMin = iR_DIID(3,kDir,Min_,iDim)
               kRMax = iR_DIID(3,kDir,Max_,iDim)
               iSMin = iS_DIID(1,iDir,Min_,iDim)
               iSMax = iS_DIID(1,iDir,Max_,iDim)
               jSMin = iS_DIID(2,jDir,Min_,iDim)
               jSMax = iS_DIID(2,jDir,Max_,iDim)
               kSMin = iS_DIID(3,kDir,Min_,iDim)
               kSMax = iS_DIID(3,kDir,Max_,iDim)
               BufferS_I(            iBufferS+1) = iRMin
               BufferS_I(            iBufferS+2) = iRMax
               if(nDim > 1)BufferS_I(iBufferS+3) = jRMin
               if(nDim > 1)BufferS_I(iBufferS+4) = jRMax
               if(nDim > 2)BufferS_I(iBufferS+5) = kRMin
               if(nDim > 2)BufferS_I(iBufferS+6) = kRMax

               iBufferS = iBufferS + 2*nDim

               do k = kSMin,kSmax; do j = jSMin,jSMax; do i = iSMin,iSmax
                  iBufferS = iBufferS + 1
                  BufferS_I(iBufferS) = Current_FDB(i,j,k,iDim,iBlockSend)
               end do; end do; end do
            end do
            iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_equal

    !==========================================================================

    subroutine set_range

      integer:: iDim
      !------------------------------------------------------------------------

            !Array index is the coordinate of the gridpoint with +1/2 being
      !approximated as 1. By x the physical cells are marked below
      !  Left corner        Right corner n+2                 Ez (:,:,nZ+2)    
      !       _!_!x!x               _!_!_!_n+2     Not used: Ex (nX+2,:,:)
      !       _!_!x!x               _!_!_!_                  Ey (:,nY+2,:)
      !       _!_!_!                x!x!_!_        Phys Face Values: 
      !    -2  ! ! !                x!x!_!_        Ex(0,1,1),Ex(nI,1,1)...
      !       -2                           Ghost values: Ex(-1,1,1),Ey(0,1,1)..
      !------------------------------------------------------------------------
      !\
      ! For x_, y_, z_ direction (the first index), the receiving block is at
      ! negative displcement with respect to the sending one
      iS_DIID(:,-1,Min_,:) = 1 - nWidth
      iS_DIID(:,-1,Max_,:) = 0
      do iDim = 1,MaxDim
         iR_DIID(:,-1,Min_,iDim) = nIjk_D + 1 -nWidth
         iR_DIID(:,-1,Max_,iDim) = nIjk_D
      end do
      !\
      ! Zero displacement
      !/
      iS_DIID(:, 0,Min_,:) = 1
      iR_DIID(:, 0,Min_,:) = 1
      do iDim = 1,MaxDim
         iS_DIID(:, 0,Max_,iDim) = nIjk_D
         iR_DIID(:, 0,Max_,iDim) = nIjk_D
         iS_DIID(iDim,0,Min_,iDim) = 0
         iR_DIID(iDim,0,Min_,iDim) = 0
      end do

      iR_DIID(:, 1,Min_,:) = 1
      iR_DIID(:, 1,Max_,:) = nWidth
  
      do iDim = 1,MaxDim
         iS_DIID(:, 1,Min_,iDim) = nIjk_D + 1
         iS_DIID(:, 1,Max_,iDIm) = nIjk_D + nWidth
         iS_DIID(iDim,1,Min_:Max_,iDim) = iS_DIID(iDim,1,Min_:Max_,iDim) - 1
         iR_DIID(iDim,1,Min_:Max_,iDim) = iR_DIID(iDim,1,Min_:Max_,iDim) - 1
      end do


    end subroutine set_range

  end subroutine add_ghost_face_field
  !=========================
  subroutine add_ghost_cell_field(nG, nVar, State_VGB, nWidthIn)
    use ModMpi

    ! Arguments
    integer, intent(in) :: nG    ! number of ghost cells for 1..nDim
    integer, intent(in) :: nVar
    real, intent(inout), dimension(1:nVar,1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
         1-nG*kDim_:nK+nG*kDim_,MaxBlock):: State_VGB
    !Array index is the coordinate of the gridpoint with +1/2 being
    !approximated as 1. By x the physical cells are marked below
    !  Left corner        Right corner n+2                 Jz (:,:,nZ+2)    
    !       _!_!x!x               _!_!_!_n+2     Not used: Jx (nX+2,:,:)
    !       _!_!x!x               _!_!_!_                  Jy (:,nY+2,:)
    !       _!_!_!                x!x!_!_        Phys Face Values: 
    !    -2  ! ! !                x!x!_!_        Jx(0,1,1),Jx(nI,1,1)...
    !       -2                           Ghost values: Jx(-1,1,1),Jy(0,1,1)..
    !------------------------------------------------------------------------

    ! Optional arguments
    integer, optional, intent(in) :: nWidthIn
    
    logical :: DoSendCorner

    ! Fill in the nVar variables in the ghost cells of Current_FDB.
    !
    ! nWidthIn is the number of ghost cell layers to be set. Default is all.


    ! Local variables

    logical, parameter :: UseRSend = .false.

    ! local variables corresponding to optional arguments
    integer :: nWidth

    ! Various indexes 
    integer :: iCountOnly     ! index for 2 stage scheme for count, sendrecv
    logical :: DoCountOnly    ! logical for count vs. sendrecv stages

    integer :: iDir, jDir, kDir
    integer :: iBlockRecv, iProcRecv, iBlockSend, iProcSend

    ! Index range for recv and send segments of the blocks
    integer :: iRMin, iRMax, jRMin, jRMax, kRMin, kRMax
    integer :: iSMin, iSMax, jSMin, jSMax, kSMin, kSMax

    integer :: iBufferS, iBufferR
    integer :: MaxBufferS = -1, MaxBufferR = -1
    integer:: iRequestR, iRequestS, iError
 
    character(len=*), parameter:: NameSub = 'BATL_pass::add_ghost_face_field'

    integer:: iBlock


    !--------------------------------------------------------------------------


    call timing_start('batl_pass')

    call timing_start('init_pass')

    ! Set values or defaults for optional arguments
    nWidth = nG
    if(present(nWidthIn)) nWidth = nWidthIn


    if(nWidth < 1 .or. nWidth > nG) call CON_stop(NameSub// &
         ' nWidth do not contain the ghost cells or too many')
    DoSendCorner = nWidth>1


    ! Set index ranges based on arguments
    call set_range

    if(nProc > 1)then
       ! Allocate fixed size communication arrays.
       if(.not.allocated(iBufferS_P))then
          allocate(iBufferS_P(0:nProc-1), nBufferS_P(0:nProc-1), &
               nBufferR_P(0:nProc-1))
          allocate(iRequestR_I(nProc-1), iRequestS_I(nProc-1))
          allocate(iStatus_II(MPI_STATUS_SIZE,nProc-1))
       end if
    end if


    call timing_stop('init_pass')


    do iCountOnly = 1, 2
       DoCountOnly = iCountOnly == 1

       ! No need to count data for send/recv in serial runs
       if(nProc == 1 .and. DoCountOnly) CYCLE

       call timing_start('local_pass')

       ! Second order prolongation needs two stages: 
       ! first stage fills in equal and coarser ghost cells
       ! second stage uses these to prolong and fill in finer ghost cells

       if(nProc>1)then
          if(DoCountOnly)then
             ! initialize buffer size
             nBufferR_P = 0
             nBufferS_P = 0
          else
             ! Make buffers large enough
             if(sum(nBufferR_P) > MaxBufferR) then
                if(allocated(BufferR_I)) deallocate(BufferR_I)
                MaxBufferR = sum(nBufferR_P) 
                allocate(BufferR_I(MaxBufferR))
             end if

             if(sum(nBufferS_P) > MaxBufferS) then
                if(allocated(BufferS_I)) deallocate(BufferS_I)
                MaxBufferS = sum(nBufferS_P)
                allocate(BufferS_I(MaxBufferS))
             end if

             ! Initialize buffer indexes for storing data into BufferS_I
             iBufferS = 0
             do iProcRecv = 0, nProc-1
                iBufferS_P(iProcRecv) = iBufferS
                iBufferS = iBufferS + nBufferS_P(iProcRecv)
             end do
          end if
       end if
       ! Loop through all blocks that may send a message
       do iBlockSend = 1, nBlock
          if(Unused_B(iBlockSend)) CYCLE
          do kDir = -1, 1
             ! Do not message pass in ignored dimensions
             if(nDim < 3 .and. kDir /= 0) CYCLE
             do jDir = -1, 1
                if(nDim < 2 .and. jDir /= 0) CYCLE
                ! Skip edges
                if(.not.DoSendCorner .and. jDir /= 0 .and. kDir /= 0) &
                     CYCLE
                do iDir = -1,1
                   ! Ignore inner parts of the sending block
                   if(iDir == 0 .and. jDir == 0 .and. kDir == 0) CYCLE
                   if(DiLevelNei_IIIB(iDir,jDir,kDir,iBlockSend)/= Unset_)&
                        call do_equal
                end do ! iDir
             end do ! jDir
          end do ! kDir
       end do ! iBlockSend

       call timing_stop('local_pass')

    end do ! iCountOnly

    call timing_start('recv_pass')

    ! post requests
    iRequestR = 0
    iBufferR  = 1
    do iProcSend = 0, nProc - 1
       if(nBufferR_P(iProcSend) == 0) CYCLE
       iRequestR = iRequestR + 1

       call MPI_irecv(BufferR_I(iBufferR), nBufferR_P(iProcSend), &
            MPI_REAL, iProcSend, 10, iComm, iRequestR_I(iRequestR), &
            iError)

       iBufferR  = iBufferR  + nBufferR_P(iProcSend)
    end do

    call timing_stop('recv_pass')

    if(UseRSend) then
       call timing_start('barrier_pass')
       call barrier_mpi
       call timing_stop('barrier_pass')
    end if

    call timing_start('send_pass')

    ! post sends
    iRequestS = 0
    iBufferS  = 1
    do iProcRecv = 0, nProc-1
       if(nBufferS_P(iProcRecv) == 0) CYCLE
       iRequestS = iRequestS + 1

       if(UseRSend)then
          call MPI_rsend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iError)
       else
          call MPI_isend(BufferS_I(iBufferS), nBufferS_P(iProcRecv), &
               MPI_REAL, iProcRecv, 10, iComm, iRequestS_I(iRequestS), &
               iError)
       end if

       iBufferS  = iBufferS  + nBufferS_P(iProcRecv)
    end do
    call timing_stop('send_pass')

    call timing_start('wait_pass')
    ! wait for all requests to be completed
    if(iRequestR > 0) &
         call MPI_waitall(iRequestR, iRequestR_I, iStatus_II, iError)

    ! wait for all sends to be completed
    if(.not.UseRSend .and. iRequestS > 0) &
         call MPI_waitall(iRequestS, iRequestS_I, iStatus_II, iError)
    call timing_stop('wait_pass')

    call timing_start('buffer_to_state')
    call buffer_to_state
    call timing_stop('buffer_to_state')


    call timing_stop('batl_pass')

  contains

    !======================================================================
    subroutine buffer_to_state

      ! Copy buffer into recv block of Current_FDB

      integer:: iBufferR, i, j, k, iDim
      !------------------------------------------------------------------------
      jRMin = 1; jRMax = 1
      kRMin = 1; kRMax = 1

      iBufferR = 0
      do iProcSend = 0, nProc - 1
         if(nBufferR_P(iProcSend) == 0) CYCLE

         do
            iBufferR = iBufferR + 1
            iBlockRecv = nint(BufferR_I(iBufferR))
            do iDim = 1, MaxDim
               iRMin      = nint(BufferR_I(iBufferR+1))
               iRMax      = nint(BufferR_I(iBufferR+2))
               if(nDim > 1) jRMin = nint(BufferR_I(iBufferR+3))
               if(nDim > 1) jRMax = nint(BufferR_I(iBufferR+4))
               if(nDim > 2) kRMin = nint(BufferR_I(iBufferR+5))
               if(nDim > 2) kRMax = nint(BufferR_I(iBufferR+6))
               
               iBufferR = iBufferR + 2*nDim
               
               
               do k=kRMin,kRmax; do j=jRMin,jRMax; do i=iRMin,iRmax 
                  State_VGB(:,i,j,k,iBlockRecv) =      &
                       State_VGB(:,i,j,k,iBlockRecv) + &
                       BufferR_I(iBufferR+1:iBufferR+nVar)
                  iBufferR = iBufferR + nVar
               end do; end do; end do
            end do
            if(iBufferR >= sum(nBufferR_P(0:iProcSend))) EXIT
         end do
       end do
    end subroutine buffer_to_state

    !==========================================================================

    subroutine do_equal

      integer :: nSize, iDim, iBufferS, i, j, k, iSend, jSend, kSend, iNodeRecv
      !------------------------------------------------------------------------

      iSend = (3*iDir + 3)/2
      jSend = (3*jDir + 3)/2
      kSend = (3*kDir + 3)/2

      iNodeRecv  = iNodeNei_IIIB(iSend,jSend,kSend,iBlockSend)
      iProcRecv  = iTree_IA(Proc_,iNodeRecv)
      iBlockRecv = iTree_IA(Block_,iNodeRecv)

      
      if(DoCountOnly)then
         ! No need to count data for local copy
         if(iProc == iProcRecv) RETURN
         !\
         ! Number of integers to be send/received
         !/
         nSize = 1 + 2*nDim*MaxDim
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)

 
            ! Number of reals to send to and received from the other processor
            nSize = (iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1) + nSize
         end do
         nBufferR_P(iProcRecv) = nBufferR_P(iProcRecv) + nSize
         nBufferS_P(iProcRecv) = nBufferS_P(iProcRecv) + nSize
         RETURN
      end if
      if(iProc == iProcRecv)then
            ! Local copy
         do iDim = 1,MaxDim
            iRMin = iR_DIID(1,iDir,Min_,iDim)
            iRMax = iR_DIID(1,iDir,Max_,iDim)
            jRMin = iR_DIID(2,jDir,Min_,iDim)
            jRMax = iR_DIID(2,jDir,Max_,iDim)
            kRMin = iR_DIID(3,kDir,Min_,iDim)
            kRMax = iR_DIID(3,kDir,Max_,iDim)
            iSMin = iS_DIID(1,iDir,Min_,iDim)
            iSMax = iS_DIID(1,iDir,Max_,iDim)
            jSMin = iS_DIID(2,jDir,Min_,iDim)
            jSMax = iS_DIID(2,jDir,Max_,iDim)
            kSMin = iS_DIID(3,kDir,Min_,iDim)
            kSMax = iS_DIID(3,kDir,Max_,iDim)
            if((iRMax-iRMin+1)*(jRMax-jRMin+1)*(kRMax-kRMin+1)==0)CYCLE
         
            State_VGB(:,iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iBlockRecv) =    &
                 State_VGB(:,iRMin:iRMax,jRMin:jRMax,kRMin:kRMax,iBlockRecv) &
                 + State_VGB(:,iSMin:iSMax,jSMin:jSMax,kSMin:kSMax,iBlockSend)
         end do
      else
            ! Put data into the send buffer
            iBufferS = iBufferS_P(iProcRecv) +1

            BufferS_I(            iBufferS) = iBlockRecv
            do iDim = 1,MaxDim
               iRMin = iR_DIID(1,iDir,Min_,iDim)
               iRMax = iR_DIID(1,iDir,Max_,iDim)
               jRMin = iR_DIID(2,jDir,Min_,iDim)
               jRMax = iR_DIID(2,jDir,Max_,iDim)
               kRMin = iR_DIID(3,kDir,Min_,iDim)
               kRMax = iR_DIID(3,kDir,Max_,iDim)
               iSMin = iS_DIID(1,iDir,Min_,iDim)
               iSMax = iS_DIID(1,iDir,Max_,iDim)
               jSMin = iS_DIID(2,jDir,Min_,iDim)
               jSMax = iS_DIID(2,jDir,Max_,iDim)
               kSMin = iS_DIID(3,kDir,Min_,iDim)
               kSMax = iS_DIID(3,kDir,Max_,iDim)
               BufferS_I(            iBufferS+1) = iRMin
               BufferS_I(            iBufferS+2) = iRMax
               if(nDim > 1)BufferS_I(iBufferS+3) = jRMin
               if(nDim > 1)BufferS_I(iBufferS+4) = jRMax
               if(nDim > 2)BufferS_I(iBufferS+5) = kRMin
               if(nDim > 2)BufferS_I(iBufferS+6) = kRMax

               iBufferS = iBufferS + 2*nDim

               do k = kSMin,kSmax; do j = jSMin,jSMax; do i = iSMin,iSmax
                  BufferS_I(iBufferS+1:iBufferS+nVar) = &
                       State_VGB(:,i,j,k,iBlockSend)
                  iBufferS = iBufferS + nVar
               end do; end do; end do
            end do
            iBufferS_P(iProcRecv) = iBufferS

      end if

    end subroutine do_equal

    !==========================================================================

    subroutine set_range

      integer:: iDim
      !------------------------------------------------------------------------

            !Array index is the coordinate of the gridpoint with +1/2 being
      !approximated as 1. By x the physical cells are marked below
      !  Left corner        Right corner n+2                 Ez (:,:,nZ+2)    
      !       _!_!x!x               _!_!_!_n+2     Not used: Ex (nX+2,:,:)
      !       _!_!x!x               _!_!_!_                  Ey (:,nY+2,:)
      !       _!_!_!                x!x!_!_        Phys Face Values: 
      !    -2  ! ! !                x!x!_!_        Ex(0,1,1),Ex(nI,1,1)...
      !       -2                           Ghost values: Ex(-1,1,1),Ey(0,1,1)..
      !------------------------------------------------------------------------
      !\
      ! For x_, y_, z_ direction (the first index), the receiving block is at
      ! negative displcement with respect to the sending one
      iS_DIID(:,-1,Min_,:) = 1
      iS_DIID(:,-1,Max_,:) = nWidth
      do iDim = 1,MaxDim
         iR_DIID(:,-1,Min_,iDim) = nIjk_D + 1
         iR_DIID(:,-1,Max_,iDim) = nIjk_D + nWidth
         iS_DIID(iDim,-1,Max_,iDim) = iS_DIID(iDim,-1,Max_,iDim) - 1
         iR_DIID(iDim,-1,Max_,iDim) = iR_DIID(iDim,-1,Max_,iDim) - 1
      end do
      iS_DIID(:, 0,Min_,:) = 1
      iR_DIID(:, 0,Min_,:) = 1
      do iDim = 1,MaxDim
         iS_DIID(:, 0,Max_,iDim) = nIjk_D
         iR_DIID(:, 0,Max_,iDim) = nIjk_D
         iS_DIID(iDim,0,Min_,iDim) = 0
         iR_DIID(iDim,0,Min_,iDim) = 0
      end do

      iR_DIID(:, 1,Min_,:) = 1 - nWidth
      iR_DIID(:, 1,Max_,:) = 0
      do iDim = 1,MaxDim
         iS_DIID(:, 1,Min_,iDim) = nIjk_D + 1 - nWidth
         iS_DIID(:, 1,Max_,iDim) = nIjk_D
         iS_DIID(iDim,1,Max_,iDim) = iS_DIID(iDim,1,Max_,iDim) - 1
         iR_DIID(iDim,1,Max_,iDim) = iR_DIID(iDim,1,Max_,iDim) - 1
      end do


    end subroutine set_range

  end subroutine add_ghost_cell_field
  !==================
end module BATL_pass_face_field
