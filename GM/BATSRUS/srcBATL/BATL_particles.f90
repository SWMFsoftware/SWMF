!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_particles

  use BATL_size, ONLY: nDim, MaxDim
  use BATL_tree, ONLY: Unset_
  use BATL_grid, ONLY: check_interpolate => check_interpolate_amr_gc

  implicit none
  SAVE
  logical  ::UseParticles = .false.
  !\
  ! Number of different sorts of particles
  !/ 
  integer, parameter:: nPType = 1
  !\
  ! Maximum particle number per processor per sort
  !/
  integer, parameter:: nParticleMax = 10000000

  type particles
     !\
     ! The number of parameters which characterize the 'state' of
     ! the particles
     !/
     integer:: nVar 
     integer:: nParticle
     !\
     ! nVar*nParticleMax array. The second index numerates 'particles'.
     !First nDim components are the cartesian coordinates the other 
     !can be are velocities, the momentum components, gyrokinetic
     !parameters or the physical particles, or, alternatively,
     !the magnetic field line ID and the number of the Lagrangian
     !fluid particle along the magnetic field line in application to
     !M-FLAMPA
     real,    pointer  :: State_VI(:,:)
     !\
     ! Array of the length of nParticleMax with the index enumerating
     ! particles. It stores the number of block which posesses an 
     ! information, sufficient to interpolate (with the possible use 
     ! of ghost cells) an information 
     integer,  pointer :: iBlock_I(:)
  end type particles
  type(particles),dimension(nPType) ::Of_I
  
contains
  !================
  subroutine set_pointer_to_particles(&
       iSortParticle, State_VI, iBlock_I, nVar, nParticle)
    integer,          intent(in)    :: iSortParticle
    real,    pointer, intent(inout) :: State_VI(:,:)
    integer, pointer, intent(inout) :: iBlock_I(:)
    integer,          intent(out)   :: nVar 
    integer,          intent(out) :: nParticle
    !------------
    nullify(State_VI); nullify(iBlock_I)
    State_VI  => Of_I(iSortParticle)%State_VI
    iBlock_I  => Of_I(iSortParticle)%iBlock_I
    nVar      =  Of_I(iSortParticle)%nVar
    nParticle =  Of_I(iSortParticle)%nParticle
  end subroutine set_pointer_to_particles
  !=================
  subroutine message_pass_particles
    use ModMpi
    use BATL_mpi, ONLY: iProc, nProc, iComm
    ! this subroutine passes particles between processors
    ! based on whether it is possible to interpolate background data
    ! to the current particle location
    !--------------------------------------------------------------------------
    integer          :: iSortParticle  ! loop variable
    real,    pointer :: State_VI(:,:)  ! state vec for particles of this sort
    integer, pointer :: iBlock_I(:)    ! blocks having particles of this sort 
    integer          :: nVar           ! # of variables including coordinates
    integer          :: nParticle      ! # of particles of this sort on proc

    ! number of particles to send to other procs
    integer:: nSend_P(nProc) 
    ! number of particles to recv by particle sort from other procs
    integer:: nRecv_P(nProc)
    ! offset for data to be sent/recv'd by procs in the BufferSend_I
    ! NOTE: starts with 0 (ZERO)
    integer:: iSendOffset_P(nProc)
    integer:: iRecvOffset_P(nProc)

    ! max number of variables 
    integer:: nVarMax         
    ! offset for particle data in the send BufferSend_I
    ! NOTE: offest values start with 0 (ZERO)
    integer, allocatable:: iSendOffset_I(:)  
    integer, allocatable:: iPeSend_I(:)! proc id to which send this particle

    real, allocatable:: BufferSend_I(:)! buffer of data to be sent
    integer          :: nBufferSend    ! size of BufferSend_I
    real, allocatable:: BufferRecv_I(:)! buffer of data to be recv'd
    integer          :: nBufferRecv    ! size of BufferRecv_I

    !-------------------------------------------------------------------------
    ! in order to reuse the same allocatable buffers for sending
    ! different sorts of particles, find the sort with MAX number of variables,
    ! nParticleMax as number of particles and allocate buffers accordingly
    nVarMax = 0
    do iSortParticle = 1, nPType
       if(nVarMax < Of_I(iSortParticle)%nVar)&
            nVarMax =  Of_I(iSortParticle)%nVar
    end do
    ! size of the buffer is (nParticles)*(nVar+1) with last 1 for block id
    nBufferSend = nParticleMax * (nVarMax+1)
    nBufferRecv = nParticleMax * (nVarMax+1)

    ! allocate buffers for send/recving data
    allocate(BufferSend_I(nBufferSend))
    allocate(BufferRecv_I(nBufferRecv))
    allocate(iSendOffset_I(nParticleMax))
    allocate(iPeSend_I(nParticleMax))


    ! now buffers are allocated, perform pass for all sorts
    do iSortParticle = 1, nPType
       call set_pointer_to_particles(&
            iSortParticle, State_VI, &
            iBlock_I, nVar,nParticle)
       call pass_this_sort
       Of_I(iSortParticle)%nParticle = nParticle
    end do
    ! deallocate buffer
    deallocate(BufferSend_I, iSendOffset_I, iPeSend_I, BufferRecv_I)


  contains

    subroutine pass_this_sort
      integer:: iParticle    ! loop variable
      real   :: Xyz_D(MaxDim)  ! particle coordinates
      logical:: IsPossible   ! can interpolate to Xyz_D on current block
      integer:: iBlock       ! current block containing the particle
      integer:: iBlockOut    ! block to be used to interpolate to Xyz_D
      integer:: iPeOut       ! proc that has iBlockOut
      integer:: iPe          ! loop variable
      integer:: iSendOffset  ! start position in BufferSend for particle data
      integer:: iRecvOffset  ! start position in BufferRecv for particle data
      integer:: nParticleStay! # of particles that stay on this proc
      logical:: IsOut        ! particle is out of domain
      integer:: iTag, iError, iRequest, iRequest_I(2*nProc)
      integer:: iStatus_II(MPI_STATUS_SIZE, 2*nProc)
      !------------------------------------------------------------------------
      ! reset parameters of the message_pass for this sort of particles
      nSend_P       = 0
      iSendOffset_I =-1
      iPeSend_I     =-1

      ! cycle through particles & find which should be passed to other procs
      do iParticle = 1, nParticle
         Xyz_D = 0
         Xyz_D(1:nDim)  = State_VI(1:nDim, iParticle)
         iBlock = iBlock_I(iParticle)
         call check_interpolate(Xyz_D, iBlock, &
              IsPossible, iPeOut, iBlockOut, IsOut)
         
         if(IsPossible) CYCLE ! don't need to pass this particle
         
         if(IsOut)then ! particle is out of the domain, don't pass it
            iBlock_I(iParticle) = Unset_
            CYCLE
         end if

         ! change the block 
         iBlock_I(iParticle) = iBlockOut
         
         if(iPeOut == iProc) CYCLE ! particle stays on this proc

         ! prepare this particle to be passed
         iPeSend_I(    iParticle) = iPeOut
         iSendOffset_I(iParticle) = nSend_P(iPeOut) * (nVar + 1)
         nSend_P(      iPeOut)    = nSend_P(iPeOut) + 1
      end do

      ! send size of messages
      iRequest = 0
      do iPe = 1, nProc
         if(iPe == iProc) CYCLE ! skip this proc
         iRequest = iRequest + 1
         call MPI_Isend(&
              nSend_P(iPe), 1, MPI_INTEGER, iPe, iTag, iComm, &
              iRequest_I(iRequest), iError)
         iRequest = iRequest + 1
         call MPI_Irecv(&
              nRecv_P(iPe), 1, MPI_INTEGER, iPe, iTag, iComm, &
              iRequest_I(iRequest), iError)
      end do
      ! finalize transfer
      call MPI_waitall(iRequest, iRequest_I, iStatus_II, iError)


      ! change total number of particles of this sort
      nParticle = nParticle - sum(nSend_P) + sum(nRecv_P)
      if(nParticle > nParticleMax)&
           call CON_stop("Exceeded allowed number of particles per sort")

      ! get offsets for procs in BufferSend_I & BufferRecv_I
      iSendOffset_P(1) = 0
      iRecvOffset_P(1) = 0
      do iPe = 2, nProc
         iSendOffset_P(iPe) = iSendOffset_P(iPe-1) + nSend_P(iPe-1) * (nVar+1)
         iRecvOffset_P(iPe) = iRecvOffset_P(iPe-1) + nRecv_P(iPe-1) * (nVar+1)
      end do

      ! fill the BufferSend and rearrange particle data
      nParticleStay = 0
      do iParticle = 1, nParticle
         if(iPeSend_I(iParticle) == -1)then
            ! particle stays on this proc
            nParticleStay = nParticleStay + 1
            iBlock_I(  nParticleStay) = iBlock_I(  iParticle)
            State_VI(:,nParticleStay) = State_VI(:,iParticle)
            CYCLE 
         end if
         iSendOffset = &
              iSendOffset_P(iPeSend_I(iParticle)) + iSendOffset_I(iParticle)
         BufferSend_I(iSendOffset+1:iSendOffset+nVar) = State_VI(:,iParticle)
         BufferSend_I(iSendOffset+nVar+1)             = iBlock_I(  iParticle)
      end do

      ! transfer data
      iRequest = 0
      do iPe = 1, nProc
         if(iPe == iProc) CYCLE ! skip this proc
         iRequest = iRequest + 1
         call MPI_Isend(&
              BufferSend_I(iSendOffset_P(iPe)+1), nSend_P(iPe)*(nVar+1), &
              MPI_REAL, &
              iPe, iTag, iComm, iRequest_I(iRequest), iError)
         iRequest = iRequest + 1
         call MPI_Irecv(&
              BufferRecv_I(iRecvOffset_P(iPe)+1), nRecv_P(iPe)*(nVar+1), &
              MPI_REAL, &
              iPe, iTag, iComm, iRequest_I(iRequest), iError)
      end do
      ! finalize transfer
      call MPI_waitall(iRequest, iRequest_I, iStatus_II, iError)

      ! finally, put particles from buffer to storage
      iRecvOffset = 0
      do iParticle = 1, sum(nRecv_P)
         State_VI(:, nParticleStay + iParticle) = &
              BufferRecv_I(iRecvOffset+1 : iRecvOffset+nVar)
         iBlock_I(   nParticleStay + iParticle) = &
              BufferRecv_I(iRecvOffset+nVar+1)
         iRecvOffset = iRecvOffset + nVar + 1
      end do
    end subroutine pass_this_sort
  end subroutine message_pass_particles
end module BATL_particles
