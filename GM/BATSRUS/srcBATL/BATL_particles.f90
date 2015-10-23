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
  integer, parameter:: nSortParticle = 1

  !\
  ! Use the pair RSend + IRecv or ISend + IRecv
  !/
  logical, parameter, private:: DoRSend = .true.

  type ParticleType
     !\
     ! The number of parameters which characterize the 'state' of
     ! the particles
     !/
     integer:: nVar 
     !\
     ! The current number of particles at a given processor and
     ! the maximum allowed number of particles
     integer:: nParticle, nParticleMax
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
  end type ParticleType

  type(ParticleType):: Particle_I(nSortParticle)

  ! offset for particle data in the send BufferSend_I
  ! NOTE: offest values start with 0
  integer, allocatable:: iSendOffset_I(:)  
  integer, allocatable:: iProcTo_I(:)! proc id to which send this particle

  real,    allocatable:: BufferSend_I(:)! buffer of data to be sent
  real,    allocatable:: BufferRecv_I(:)! buffer of data to be recv'd

contains

  !===========================================================================
  subroutine allocate_particles
    integer:: iParticleSort, nVar, nParticleMax
    !------------------------------------------------------------------------
    do iParticleSort = 1, nSortParticle
       nVar = Particle_I(iParticleSort)%nVar
       nParticleMax = Particle_I(iParticleSort)%nParticleMax
       nullify(Particle_I(iParticleSort)%State_VI)
       allocate(Particle_I(iParticleSort)%State_VI(1:nVar,1:nParticleMax))
       Particle_I(iParticleSort)%State_VI(1:nVar,1:nParticleMax) = 0.0
       nullify(Particle_I(iParticleSort)%iBlock_I)
       allocate(Particle_I(iParticleSort)%iBlock_I(1:nParticleMax))
       Particle_I(iParticleSort)%iBlock_I(1:nParticleMax) = -1
    end do
  end subroutine allocate_particles
  !===========================================================================
  subroutine allocate_buffers
    integer          :: iSortParticle  ! loop variable
    ! max number of particles 
    integer          :: nParticleMax  
    ! max size of buffer       
    integer          :: nBuffer    ! size of BufferSend_I,BufferRecv_I
    !-----------------------------------------------------------------
    ! in order to reuse the same allocatable buffers for sending
    ! different sorts of particles, find the sort with MAX number of variables,
    ! nParticleMax as number of particles and allocate buffers accordingly
    nBuffer = 0; nParticleMax = 0
    do iSortParticle = 1, nSortParticle
       ! size of the buffer is (nParticles)*(nVar+1) with last 1 for block id
       nBuffer = max(nBuffer, (Particle_I(iSortParticle)%nVar + 1)&
            *Particle_I(iSortParticle)%nParticleMax)

       nParticleMax = max(nParticleMax, Particle_I(iSortParticle)%nParticleMax)
    end do
    ! allocate buffers for send/recving data
    allocate(BufferSend_I(nBuffer))
    allocate(BufferRecv_I(nBuffer))
    allocate(iSendOffset_I(nParticleMax))
    allocate(iProcTo_I(nParticleMax))
    
  end subroutine allocate_buffers
  !===========================================================================
  subroutine set_pointer_to_particles(&
       iSortParticle, State_VI, iBlock_I, nVar, nParticle, nParticleMax)
    integer,          intent(in)    :: iSortParticle
    real,    pointer, intent(inout) :: State_VI(:,:)
    integer, pointer, intent(inout) :: iBlock_I(:)
    integer,          intent(out)   :: nVar 
    integer,          intent(out)   :: nParticle
    integer,          intent(out)   :: nParticleMax
    !-----------------------------------------------------------------------
    nullify(State_VI); nullify(iBlock_I)
    State_VI     => Particle_I(iSortParticle)%State_VI
    iBlock_I     => Particle_I(iSortParticle)%iBlock_I
    nVar         =  Particle_I(iSortParticle)%nVar
    nParticle    =  Particle_I(iSortParticle)%nParticle
    nParticleMax =  Particle_I(iSortParticle)%nParticleMax
  end subroutine set_pointer_to_particles
  !===========================================================================
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
    integer          :: nParticleMax   ! max # of particles of this sort on PE
    ! number of particles to send to other procs
    integer:: nSend_P(0:nProc-1) 
    ! number of particles to recv by particle sort from other procs
    integer:: nRecv_P(0:nProc-1)
    ! offset for data to be sent/recv'd by procs in the BufferSend_I
    ! NOTE: starts with 0 (ZERO)
    integer:: iSendOffset_P(0:nProc-1)
    integer:: iRecvOffset_P(0:nProc-1)
    !--------------
    if(.not.allocated(BufferSend_I))call allocate_buffers
    ! now buffers are allocated, perform pass for all sorts
    do iSortParticle = 1, nSortParticle
       call set_pointer_to_particles(&
            iSortParticle, State_VI, &
            iBlock_I, nVar, nParticle, nParticleMax)
       call pass_this_sort
       Particle_I(iSortParticle)%nParticle = nParticle
    end do
    ! deallocate buffer
    ! deallocate(BufferSend_I, iSendOffset_I, iProcTo_I, BufferRecv_I)


  contains
    !==========================================================================
    subroutine pass_this_sort
      integer:: iParticle    ! loop variable
      real   :: Xyz_D(MaxDim)! particle coordinates
      logical:: IsPossible   ! can interpolate to Xyz_D on current block
      integer:: iBlock       ! current block containing the particle
      integer:: iBlockOut    ! block to be used to interpolate to Xyz_D
      integer:: iProcOut     ! proc that has iBlockOut
      integer:: iProcTo      ! loop variable
      integer:: iProcFrom    ! loop variable
      integer:: iSendOffset  ! start position in BufferSend for particle data
      integer:: iRecvOffset  ! start position in BufferRecv for particle data
      integer:: nParticleStay! # of particles that stay on this proc
      logical:: IsOut        ! particle is out of domain
      integer:: iTag, iError, iRequest, iRequest_I(2*nProc)
      integer:: iStatus_II(MPI_STATUS_SIZE, 2*nProc)
      !-----------------------------------------------------------------------
      ! reset parameters of the message_pass for this sort of particles
      nSend_P       = 0; nRecv_P = 0
      iSendOffset_I =-1
      iProcTo_I     =-1

      ! cycle through particles & find which should be passed to other procs
      do iParticle = 1, nParticle
         Xyz_D = 0
         Xyz_D(1:nDim)  = State_VI(1:nDim, iParticle)
         iBlock = iBlock_I(iParticle)
         call check_interpolate(Xyz_D, iBlock, &
              IsPossible, iProcOut, iBlockOut, IsOut)
         
         if(IsPossible) CYCLE ! don't need to pass this particle
         
         if(IsOut)then ! particle is out of the domain, don't pass it
            iBlock_I(iParticle) = Unset_
            CYCLE
         end if

         ! change the block 
         iBlock_I(iParticle) = iBlockOut
         
         if(iProcOut == iProc) CYCLE ! particle stays on this proc

         ! prepare this particle to be passed
         iProcTo_I(    iParticle) = iProcOut
         iSendOffset_I(iParticle) = nSend_P(iProcOut) * (nVar + 1)
         nSend_P(      iProcOut)  = nSend_P(iProcOut) + 1
      end do

      ! send size of messages
      iRequest = 0
      do iProcFrom = 0, nProc - 1
         if(iProcFrom == iProc) CYCLE ! skip this proc
         iTag = iProc
         iRequest = iRequest + 1
         call MPI_Irecv(&
              nRecv_P(iProcFrom), 1, MPI_INTEGER, iProcFrom, iTag, iComm, &
              iRequest_I(iRequest), iError)
      end do

      if(DoRSend)then
         ! barrier: to guarantee that all recv's have been posted BEFORE Rsend
         call MPI_Barrier(iComm, iError)
         
         do iProcTo = 0, nProc - 1
            if(iProcTo == iProc) CYCLE ! skip this proc
            iTag = iProcTo
            call MPI_Rsend(&
                 nSend_P(iProcTo), 1, MPI_INTEGER, &
                 iProcTo, iTag, iComm, iError)
         end do
      else
         do iProcTo = 0, nProc - 1
            iTag = iProcTo
            iRequest = iRequest + 1
            call MPI_Isend(&
                 nSend_P(iProcTo), 1, MPI_INTEGER, &
                 iProcTo, iTag, iComm, iRequest_I(iRequest), iError)
         end do
      end if

      ! finalize transfer
      call MPI_waitall(iRequest, iRequest_I, iStatus_II, iError)

      ! get offsets for procs in BufferSend_I & BufferRecv_I
      iSendOffset_P(0) = 0
      do iProcTo = 1, nProc - 1
         iSendOffset_P(iProcTo) = &
              iSendOffset_P(iProcTo-1) + nSend_P(iProcTo-1) * (nVar+1)
      end do
      iRecvOffset_P(0) = 0
      do iProcFrom = 1, nProc - 1
         iRecvOffset_P(iProcFrom) = &
              iRecvOffset_P(iProcFrom-1) + nRecv_P(iProcFrom-1) * (nVar+1)
      end do

      ! fill the BufferSend and rearrange particle data
      nParticleStay = 0
      do iParticle = 1, nParticle
         if(iProcTo_I(iParticle) == -1)then
            ! particle stays on this proc
            nParticleStay = nParticleStay + 1
            iBlock_I(  nParticleStay) = iBlock_I(  iParticle)
            State_VI(:,nParticleStay) = State_VI(:,iParticle)
            CYCLE 
         end if
         iSendOffset = &
              iSendOffset_P(iProcTo_I(iParticle)) + iSendOffset_I(iParticle)
         BufferSend_I(iSendOffset+1:iSendOffset+nVar) = State_VI(:,iParticle)
         BufferSend_I(iSendOffset+nVar+1)             = iBlock_I(  iParticle)
      end do

      ! transfer data
      iRequest = 0
      do iProcFrom = 0, nProc - 1
         if(iProcFrom == iProc) CYCLE ! skip this proc
         if(nRecv_P(iProcFrom) > 0)then
            iRequest = iRequest + 1
            iTag = iProc
            call MPI_Irecv(&
                 BufferRecv_I(iRecvOffset_P(iProcFrom)+1), &
                 nRecv_P(iProcFrom)*(nVar+1), MPI_REAL, &
                 iProcFrom, iTag, iComm, iRequest_I(iRequest), iError)
         end if
      end do

      if(DoRSend)then
         ! barrier: to guarantee that all recv's have been posted BEFORE Rsend
         call MPI_Barrier(iComm, iError)
         
         do iProcTo = 0, nProc - 1
            if(iProcTo == iProc) CYCLE ! skip this proc
            if(nSend_P(iProcTo) > 0)then
               iTag = iProcTo
               call MPI_Rsend(&
                    BufferSend_I(iSendOffset_P(iProcTo)+1), &
                    nSend_P(iProcTo)*(nVar+1), MPI_REAL, &
                    iProcTo, iTag, iComm, iError)
            end if
         end do
      else
         do iProcTo = 0, nProc - 1
            if(nSend_P(iProcTo) > 0)then
               iTag = iProcTo
               iRequest = iRequest + 1
               call MPI_Isend(&
                    BufferSend_I(iSendOffset_P(iProcTo)+1), &
                    nSend_P(iProcTo)*(nVar+1), MPI_REAL, &
                    iProcTo, iTag, iComm, iRequest_I(iRequest), iError)
            end if
         end do
      end if
      ! finalize transfer
      call MPI_waitall(iRequest, iRequest_I, iStatus_II, iError)

      ! change total number of particles of this sort
      nParticle = nParticle - sum(nSend_P) + sum(nRecv_P)
      if(nParticle > nParticleMax)&
           call CON_stop("Exceeded allowed number of particles per sort=",&
           iSortParticle)

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
