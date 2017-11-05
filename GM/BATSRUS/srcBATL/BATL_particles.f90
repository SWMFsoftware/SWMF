!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module BATL_particles

  use BATL_size, ONLY: nDim, MaxDim, nKindParticle
  use BATL_grid, ONLY: check_interpolate => check_interpolate_amr_gc, &
       CoordMin_D, CoordMax_D
  use BATL_geometry, ONLY: IsCartesian, IsPeriodic_D
  implicit none

  private ! except

  ! Public methods and variables of this module
  public:: Particle_I
  public:: allocate_particles
  public:: set_pointer_to_particles
  public:: message_pass_particles
  public:: mark_undefined
  public:: remove_undefined_particles
  public:: get_particles
  public:: trace_particles
  SAVE
  logical  ::UseParticles = .false.
  !\
  ! Use the pair RSend + IRecv or ISend + IRecv
  !/
  logical, parameter, private:: DoRSend = .true.

  type ParticleType
     !\
     ! The number of parameters which characterize the 'state' of
     ! the particles
     !/
     integer:: nVar   ! # of real    parameters
     integer:: nIndex ! # of integer parameters
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
     ! (nIndex+1)*nParticleMax array with the indices enumerating
     ! particles. It ALWAYS stores the number of block (0th index) 
     ! which posesses an information, sufficient to interpolate 
     ! (with the possible use of ghost cells) an information 
     integer,  pointer :: iIndex_II(:,:)
  end type ParticleType

  type(ParticleType):: Particle_I(nKindParticle)

  ! offset for particle data in the send BufferSend_I
  ! NOTE: offest values start with 0
  integer, allocatable:: iSendOffset_I(:)  
  integer, allocatable:: iProcTo_I(:)! proc id to which send this particle

  real,    allocatable:: BufferSend_I(:)! buffer of data to be sent
  real,    allocatable:: BufferRecv_I(:)! buffer of data to be recv'd
  character(len=*),parameter  :: NameMod='BATL_particles'
  character(len=120), private :: StringError
  integer:: nU_I(2)
contains

  !===========================================================================
  subroutine allocate_particles
    integer:: iKindParticle, nVar, nIndex, nParticleMax
    !------------------------------------------------------------------------
    do iKindParticle = 1, nKindParticle
       nVar         = Particle_I(iKindParticle)%nVar
       nIndex       = Particle_I(iKindParticle)%nIndex
       nParticleMax = Particle_I(iKindParticle)%nParticleMax
       nullify( Particle_I(iKindParticle)%State_VI)
       allocate(Particle_I(&
               iKindParticle)%State_VI(1:nVar,1:nParticleMax))

       nullify( Particle_I(iKindParticle)%iIndex_II)
       allocate(Particle_I(&
            iKindParticle)%iIndex_II(0:nIndex,1:nParticleMax))
       call clean_particle_arr(iKindParticle, 1, nParticleMax)
    end do
  end subroutine allocate_particles
  !================================
  subroutine clean_particle_arr(iKindParticle, iParticleMin, iParticleMax)
    integer, intent(in) :: iKindParticle, iParticleMin, iParticleMax
    if(iParticleMin > iParticleMax) RETURN
    Particle_I(iKindParticle)%State_VI(:,iParticleMin:iParticleMax) = 0.0
    Particle_I(iKindParticle)%iIndex_II(:,iParticleMin:iParticleMax) = 0 
  end subroutine clean_particle_arr
  !================================
  subroutine allocate_buffers
    integer          :: iKindParticle  ! loop variable
    ! max number of particles 
    integer          :: nParticleMax  
    ! max size of buffer       
    integer          :: nBuffer    ! size of BufferSend_I,BufferRecv_I
    !-----------------------------------------------------------------
    ! in order to reuse the same allocatable buffers for sending
    ! different kinds of particles, find the kind with MAX number of variables,
    ! nParticleMax as number of particles and allocate buffers accordingly
    nBuffer = 0; nParticleMax = 0
    do iKindParticle = 1, nKindParticle
       ! size of the buffer is (nParticles)*(nVar+1) with last 1 for block id
       nBuffer = max(nBuffer, &
            (Particle_I(iKindParticle)%nVar+Particle_I(iKindParticle)%nIndex+&
            1)*Particle_I(iKindParticle)%nParticleMax)

       nParticleMax = max(nParticleMax, Particle_I(iKindParticle)%nParticleMax)
    end do
    ! allocate buffers for send/recving data
    allocate(BufferSend_I(nBuffer))
    allocate(BufferRecv_I(nBuffer))
    allocate(iSendOffset_I(nParticleMax))
    allocate(iProcTo_I(nParticleMax))
    
  end subroutine allocate_buffers
  !===========================================================================
  subroutine set_pointer_to_particles(&
       iKindParticle, State_VI, iIndex_II, &
       nVar, nIndex, nParticle, nParticleMax)
    integer,          intent(in)    :: iKindParticle
    real,    pointer, intent(inout) :: State_VI(:,:)
    integer, pointer, optional, intent(inout) :: iIndex_II(:,:)
    integer, optional,intent(out)   :: nVar 
    integer, optional,intent(out)   :: nIndex 
    integer, optional,intent(out)   :: nParticle
    integer, optional,intent(out)   :: nParticleMax
    !-----------------------------------------------------------------------
    nullify(State_VI)
    State_VI     => Particle_I(iKindParticle)%State_VI
    if(present(iIndex_II))then
       nullify(iIndex_II)
       iIndex_II => Particle_I(iKindParticle)%iIndex_II
    end if
    if(present(nVar))&
         nVar         =  Particle_I(iKindParticle)%nVar
    if(present(nIndex))&
         nIndex       =  Particle_I(iKindParticle)%nIndex
    if(present(nParticle))&
         nParticle    =  Particle_I(iKindParticle)%nParticle
    if(present(nParticleMax))&
         nParticleMax =  Particle_I(iKindParticle)%nParticleMax
  end subroutine set_pointer_to_particles
  !===========================================================================
  subroutine mark_undefined(iKind, iParticle)
    ! mark the particle of the given kind as undefined
    integer, intent(in):: iKind, iParticle
    !------------------------------------------------
    ! particle is considered undefined if its block number is negative,
    ! the absolute value is kept the same to retain information on its position
    Particle_I(iKind)%iIndex_II(0,iParticle) = &
         - abs(Particle_I(iKind)%iIndex_II(0,iParticle))
    ! abs() is used to prevent making an undefined particle a defined one
  end subroutine mark_undefined
  !===========================================================================
  subroutine remove_undefined_particles(iKindParticle)
    ! remove all particles with undefined block: iBlock_I<0
    integer, intent(in) :: iKindParticle
    integer :: iVar, iIndex ! loop variables
    real,    pointer :: State_VI(:,:)  ! state vec for particles of this kind
    integer, pointer :: iIndex_II(:,:) ! indices   for particles of this kind 
    integer          :: nVar           ! # of variables including coordinates
    integer          :: nIndex         ! # of indices including block number
    integer          :: nParticle      ! # of particles of this kind on proc
    integer          :: nParticleMax   ! max # of particles of this kind on PE
    integer          :: nUnset         ! # of particles with undefined block
    !-------------------------------------------------------------------------
    !\
    !Return, if no particles on the given Proc
    if(Particle_I(iKindParticle)%nParticle < 1)RETURN
   
    call set_pointer_to_particles(iKindParticle, &
         State_VI, iIndex_II, nVar, nIndex, nParticle, nParticleMax)
 
    nUnset = count(iIndex_II(0,1:nParticle)<0)
    !\
    ! Return, if there is no underfined particle
    if(nUnset==0) RETURN
    !\
    ! Set nParticle to 0 and return if all particles are undefined
    if(nParticle == nUNset)then
       call clean_particle_arr(iKindParticle, 1, nParticle)
       Particle_I(iKindParticle)%nParticle = 0
       RETURN
    end if
    do iVar = 1, nVar
       State_VI(iVar, 1:(nParticle-nUnset)) = PACK(&
            State_VI(iVar, 1:nParticle), &
            iIndex_II(0,   1:nParticle)>0 )
    end do
    do iIndex = nIndex, 0, -1
       iIndex_II(iIndex, 1:(nParticle-nUnset)) = PACK(&
            iIndex_II(iIndex, 1:nParticle), &
            iIndex_II(0,      1:nParticle)>0 )
    end do
    call clean_particle_arr(iKindParticle, nParticle - nUnset + 1, nParticle)
    Particle_I(iKindParticle)%nParticle = nParticle - nUnset
  end subroutine remove_undefined_particles
  !===========================================================================
  subroutine message_pass_particles(iKindParticleIn)
    use BATL_tree, ONLY: Unset_
    use ModMpi
    use BATL_mpi, ONLY: iProc, nProc, iComm
    ! this subroutine passes particles between processors
    ! based on whether it is possible to interpolate background data
    ! to the current particle location
    !--------------------------------------------------------------------------
    ! if present => pass only this kind
    integer, optional, intent(in):: iKindParticleIn

    integer          :: iKindParticle  ! loop variable
    real,    pointer :: State_VI(:,:)  ! state vec for particles of this kind
    integer, pointer :: iIndex_II(:,:) ! blocks having particles of this kind 
    integer          :: nVar           ! # of variables including coordinates
    integer          :: nIndex         ! # of indices including block number
    integer          :: nParticle      ! # of particles of this kind on proc
    integer          :: nParticleMax   ! max # of particles of this kind on PE
    ! number of particles to send to other procs
    integer:: nSend_P(0:nProc-1), nSendThis, nSendAll
    ! number of particles to recv by particle kind from other procs
    integer:: nRecv_P(0:nProc-1)
    ! offset for data to be sent/recv'd by procs in the BufferSend_I
    ! NOTE: starts with 0 (ZERO)
    integer:: iSendOffset_P(0:nProc-1)
    integer:: iRecvOffset_P(0:nProc-1)
    !--------------
    if(.not.allocated(BufferSend_I))call allocate_buffers
    ! now buffers are allocated, perform pass for all kinds
    do iKindParticle = 1, nKindParticle
       if(present(iKindParticleIn))then
          ! if kind is indicated => skip others
          if(iKindParticle /= iKindParticleIn) CYCLE
       end if
       call set_pointer_to_particles(&
            iKindParticle, State_VI, iIndex_II,&
            nVar, nIndex, nParticle, nParticleMax)
       call pass_this_kind
       Particle_I(iKindParticle)%nParticle = nParticle
    end do
  contains
    !==========================================================================
    subroutine pass_this_kind
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
      integer:: nParticleNew ! # of particles after message pass
      logical:: IsOut        ! particle is out of domain
      integer:: iTag, iError, iRequest, iRequest_I(2*nProc)
      integer:: iStatus_II(MPI_STATUS_SIZE, 2*nProc)
      !-----------------------------------------------------------------------
      ! reset parameters of the message_pass for this kind of particles
      nSend_P       = 0; nRecv_P = 0
      iSendOffset_I(1:nParticle) =-1
      iProcTo_I(    1:nParticle) = iProc

      ! cycle through particles & find which should be passed to other procs
      do iParticle = 1, nParticle
         Xyz_D = 0
         if(IsCartesian)then
            where(IsPeriodic_D(1:nDim))State_VI(1:nDim, iParticle) = &
                 CoordMin_D(1:nDIm) + &
                 modulo(State_VI(1:nDim, iParticle) - CoordMin_D(1:nDim), &
                 CoordMax_D(1:nDim) - CoordMin_D(1:nDim))
         end if
         Xyz_D(1:nDim)  = State_VI(1:nDim, iParticle)
         iBlock = iIndex_II(0, iParticle)
         call check_interpolate(Xyz_D, iBlock, iProcOut, iBlockOut)
         
         if(iProc == iProcOut .and. iBlock==iBlockOut)&
              CYCLE

         if(iBlockOut == Unset_)then ! particle is out of domain, don't pass it
            call mark_undefined(iKindParticle, iParticle)
            CYCLE
         end if

         ! change the block 
         iIndex_II(0, iParticle) = iBlockOut
         
         if(iProcOut == iProc) CYCLE ! particle stays on this proc

         ! prepare this particle to be passed
         iProcTo_I(    iParticle) = iProcOut
         iSendOffset_I(iParticle) = nSend_P(iProcOut) * (nVar + nIndex + 1)
         nSend_P(      iProcOut)  = nSend_P(iProcOut) + 1
      end do
      if(nProc==1)RETURN
      
      ! check if need to pass any particles in the communicator
      nSendThis = sum(nSend_P)
      call MPI_Allreduce(&
           nSendThis, nSendAll, 1, MPI_INTEGER, MPI_SUM, iComm, iError)
      if(nSendAll==0)RETURN

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
              iSendOffset_P(iProcTo-1) + nSend_P(iProcTo-1)*(nVar+nIndex+1)
      end do
      iRecvOffset_P(0) = 0
      do iProcFrom = 1, nProc - 1
         iRecvOffset_P(iProcFrom) = &
              iRecvOffset_P(iProcFrom-1) + nRecv_P(iProcFrom-1)*(nVar+nIndex+1)
      end do

      ! fill the BufferSend and rearrange particle data
      nParticleStay = 0
      do iParticle = 1, nParticle
         if(iProcTo_I(iParticle) == iProc)then
            ! particle stays on this proc
            nParticleStay = nParticleStay + 1
            iIndex_II(:, nParticleStay) = iIndex_II(:, iParticle)
            State_VI( :, nParticleStay) = State_VI( :, iParticle)
            CYCLE 
         end if
         iSendOffset = &
              iSendOffset_P(iProcTo_I(iParticle)) + iSendOffset_I(iParticle)
         BufferSend_I(iSendOffset+1:iSendOffset+nVar) = &
              State_VI( :, iParticle)
         BufferSend_I(iSendOffset+nVar+1:iSendOffset+nVar+nIndex+1) = &
              iIndex_II(:, iParticle)
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
                 nRecv_P(iProcFrom)*(nVar+nIndex+1), MPI_REAL, &
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
                    nSend_P(iProcTo)*(nVar+nIndex+1), MPI_REAL, &
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
                    nSend_P(iProcTo)*(nVar+nIndex+1), MPI_REAL, &
                    iProcTo, iTag, iComm, iRequest_I(iRequest), iError)
            end if
         end do
      end if
      ! finalize transfer
      call MPI_waitall(iRequest, iRequest_I, iStatus_II, iError)

      ! change total number of particles of this kind
      nParticleNew = nParticle - sum(nSend_P) + sum(nRecv_P)
      if(nParticleNew > nParticleMax)then
         write(StringError,'(a,i4)')&
              "Exceeded allowed number of particles per kind=",&
              iKindParticle
         call CON_stop(StringError)
      end if
      call clean_particle_arr(iKindParticle, nParticleNew +1, nParticle)
      nParticle = nParticleNew
      ! finally, put particles from buffer to storage
      iRecvOffset = 0
      do iParticle = 1, sum(nRecv_P)
         State_VI( :, nParticleStay + iParticle) = &
              BufferRecv_I(iRecvOffset+1 : iRecvOffset+nVar)
         iIndex_II(:, nParticleStay + iParticle) = &
              nint(BufferRecv_I(iRecvOffset+nVar+1:iRecvOffset+nVar+nIndex+1))
         iRecvOffset = iRecvOffset + nVar + nIndex + 1
      end do
    end subroutine pass_this_kind
  end subroutine message_pass_particles
  !==========================
  subroutine get_particles(iKindParticle, StateIn_VI, iLastIdIn, &
       iIndexIn_II, UseInputInGenCoord, DoReplace, nParticlePE)
    use BATL_mpi,      ONLY: iProc
    use BATL_geometry, ONLY: coord_to_xyz
    integer,         intent(in)  :: iKindParticle
    real,            intent(in)  :: StateIn_VI(:,:)
    integer,optional,intent(in)  :: iLastIdIn
    integer,optional,intent(in)  :: iIndexIn_II(:,:)
    logical, optional,intent(in) :: UseInputInGenCoord, DoReplace 
    integer,optional,intent(out) :: nParticlePE
    !\
    ! Data pointers for particles of a given sort
    !/
    real,                pointer :: State_VI(:,:)
    integer,             pointer :: iIndex_II(:,:)
    !\
    ! Size of data pointers
    !/
    integer  :: nVar, nIndex, nParticleOld, nParticleMax
    !\
    ! Size of input arrays
    !/
    integer  :: nVarIn, nParticleIn, nIndexIn, nU_I(2)
    !\
    ! Used if there is no input index array. In this case the particle
    ! Id, if desired, is formed as the order # of particle in the input 
    ! array(s) + iLastId 
    !/ 
    integer  :: iLastId 
    !Output parameters for check_interpolate routine:
    integer :: iProcOut, iBlockOut
    !Coordinates, may be modified by the check_interpolate routine
    real :: Xyz_D(MaxDim), Coord_D(MaxDim)
    !Loop Variables
    integer  :: iParticleIn, iVar
    !Misc
    integer  :: nParticle, iParticleNew
    logical  :: DoTransform
    character(len=*), parameter:: NameSub=NameMod//'::get_particles'
    !--------------------
    call set_pointer_to_particles(&
         iKindParticle, State_VI, iIndex_II, &
         nVar, nIndex, nParticleOld, nParticleMax)
    if(present(DoReplace))then
       if(DoReplace)then
          call clean_particle_arr(iKindParticle, 1, nParticleOld)
          nParticleOld = 0
       end if
    end if
    nU_I = ubound(StateIn_VI)
    nVarIn      = nU_I(1) 
    nParticleIn = nU_I(2)
    iLastId     = 0
    if(present(iLastIdIn)) iLastId = iLastIdIn
    if(present(iIndexIn_II))then
       nU_I = ubound(iIndexIn_II)
       nIndexIn = nU_I(1) 
    end if
    DoTransform = .false.
    if(present(UseInputInGenCoord))DoTransform = UseInputInGenCoord
    nParticle = 0
    do iParticleIn = 1, nParticleIn
       ! find block and processor suitable for interpolation
       if(DoTransform)then
          Coord_D = 0 
          Coord_D(1:nDim) = StateIn_VI(1:nDim,iParticleIn)
          call coord_to_xyz(Coord_D, Xyz_D)
       else
          Xyz_D = 0; Xyz_D(1:nDim) = StateIn_VI(1:nDim,iParticleIn)
       end if
       call check_interpolate(Xyz_D, &
            1, & ! input block ID doesn't matter
            iProcOut, iBlockOut)

       ! check whether point is outside of the domain
       if(iProcOut < 0)then
          write(StringError,'(a,es15.6,a, es15.6, a, es15.6)') &
               "Start point for a field line is outside of the domain: X = ",&
               StateIn_VI(1,iParticleIn), " Y = ", &
               StateIn_VI(2,iParticleIn), " Z = ", &
               StateIn_VI(3,iParticleIn) 
          call CON_stop(NameSub//':'//StringError)
       end if
       ! Assign particle to an appropriate processor
       if(iProc /= iProcOut) CYCLE
       nParticle = nParticle + 1
       iParticleNew = nParticleOld + nParticle
       if(iParticleNew > nParticleMax)then
          write(StringError,'(a,i4)')&
               "Exceeded allowed number of particles per kind=",&
               iKindParticle
          call CON_stop(NameSub//':'//StringError)
       end if
       State_VI(1:MaxDim,iParticleNew) = Xyz_D
       do iVar = MaxDim + 1, nVarIn
          State_VI(iVar,iParticleNew) = StateIn_VI(iVar,iParticleIn)
       end do
       if(present(iIndexIn_II))then
          !copy index array from inputs
          iIndex_II(1:nIndexIn,iParticleNew) = iIndexIn_II(:,iParticleIn)
       else
          !Create particle Id as 
          iIndex_II(min(1,nIndex),iParticleNew) = iLastId + iParticleIn 
       end if
       iIndex_II(0,iParticleNew) = iBlockOut
    end do
    Particle_I(iKindParticle)%nParticle = nParticleOld + nParticle
    if(present(nParticlePE))nParticlePE = nParticle
  end subroutine get_particles
  !==========================
  !Tracing the trajectories of particles (end points) of a given sort, 
  !which are displaced according to the function, displace_particle. 
  !For example, if the if the end points of are displaced by an 
  !elementray step size in the direction of the magnetic field, the 
  !resulting trajectories are the magnetic field lines. The displacement 
  !stops when the particle leaves the computational domain, or reaches 
  !a body, or goes beyond some "soft" boundary. The integration has been 
  !completed once all particles reach some boundary. 
  subroutine trace_particles(iKindParticle, displace_particle)
    use ModMpi
    use BATL_mpi, ONLY: iComm, nProc
    integer, intent(in) :: iKindParticle
    interface
       subroutine displace_particle(iParticle, IsEndOfSegment)
         implicit none
         integer, intent(in) :: iParticle
         logical, intent(out):: IsEndOfSegment
         !---------------
         !\
         ! IsEndOfSegment should be set to .true. if one of the
         ! following is true:
         ! (1) check_interpolate in the displaced location of the
         ! particle shows that the particle left the computational 
         ! domain. Such particle MUST be marked as "undefined":
         !\
         ! call mark_undefined(iKind, iParticle)
         !/
         ! (2) check_interpolate in the displaced location of the
         ! particle shows that the particle location can no longer 
         ! be interpolated at the given PE
         ! (3) other criteria are satisfied for ending the 
         ! trajectory of a given partiale (it reaches the internal 
         ! or "soft" boundary). Such particle MUST be 
         ! marked as "undefined":
         !\
         ! call mark_undefined(iKind, iParticle) 
         !/
       end subroutine displace_particle
    end interface
    !\
    ! Conditions for exiting a loop or the whole routine
    logical :: IsNotAnyParticle, DoStop, IsEndOfSegment
    !\
    !Loop variable
    integer :: iParticle
    !\
    !Misc
    integer :: iError
    !-----------
    !\
    !CYCLE till all particles leave the domain
    !/
    do 
       !\
       !For all particles at this PE
       !/ 
       do iParticle = 1, Particle_I(iKindParticle)%nParticle
          !Displace while the particle is at the PE domain
          SEGMENT:do 
             call displace_particle(iParticle,IsEndOfSegment)
             !\
             ! The end of segment is achieved if the particle
             ! (1) reaches the computational domain boundary.
             !     These particles should be marked using 
             !     mark_undefined procedure
             ! (2) passes to block assigned to different PE  
             if(IsEndOfSegment)EXIT SEGMENT
          end do SEGMENT
       end do
       !\
       ! Particles of the (1) kind are removed
       call remove_undefined_particles(iKindParticle)
       IsNotAnyParticle = Particle_I(iKindParticle)%nParticle == 0
       if(nProc>1)then
          call MPI_Allreduce(IsNotAnyParticle,&
               DoStop, 1, MPI_LOGICAL, MPI_LAND, iComm, iError)
       else
          DoStop = IsNotAnyParticle
       end if
       if(DoStop)RETURN
       !\
       ! Particles of (2) kind are sent to proper processor
       !/
       call message_pass_particles(iKindParticle)
    end do
  end subroutine trace_particles
  !=======================
end module BATL_particles
