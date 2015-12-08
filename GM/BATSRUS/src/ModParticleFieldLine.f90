!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModParticleFieldLine
  ! the module contains subroutine for extracting magnetic field lines 
  ! for passing to other codes;
  ! field line intergration is performed with use of BATL library including
  ! - continuous AMR interpolation
  ! - particle methods
  !============================================================================

  use BATL_lib, ONLY: &
       iProc,&
       MaxDim, nDim, &
       interpolate_grid_amr_gc, check_interpolate_amr_gc, &
       Particle_I, CellSize_DB, &
       message_pass_particles, remove_undefined_particles, allocate_particles
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes, ONLY: Bx_, By_, Bz_
  use ModB0, ONLY: get_b0
  use ModMain, ONLY: Body1
  use ModPhysics, ONLY: rBody

  implicit none
  SAVE

  private !except

  public:: init_particle_line
  public:: extract_particle_line
  public:: sort_particle_line
  public:: get_particle_data


  ! sorts of particles used to generate a magnetic field line
  integer, parameter:: &
       KindEnd_ = 1, &
       KindReg_ = 2

  ! variable in the state vector of a particle
  integer, parameter:: &
       ! coordinates of a particle
       x_    = 1, y_    = 2, z_    = 3, & 
       ! auxilary position, e.g. middle step in Runge-Kutta 2 method
       AuxX_ = 4, AuxY_ = 5, AuxZ_ = 6, & 
       ! auxilary field, e.g. stepsize
       Aux_   = 7, & 
       ! field line this particle lays on
       fl_    = 8, &
       ! index of the particle along this field line
       Index_ = 9

  ! number of variables in the state vector
  integer, parameter:: nVarParticle = 9

  ! maximum allowed number of field lines
  integer, parameter :: nFieldLineMax = 1000

  ! number of active field lines
  integer:: nFieldLine = 0

  ! keep track of number of particles per field line
  integer :: nParticleFieldLine_I(nFieldLineMax) = 0
  integer :: iParticleFieldLineOffset_I(nFieldLineMax) = 0

contains

  subroutine init_particle_line
    ! allocate containers for particles
    !------------------------------------------------------------------------
    Particle_I(KindReg_)%nParticleMax = 1000
    Particle_I(KindEnd_)%nParticleMax = 10
    Particle_I(KindReg_)%nVar = nVarParticle
    Particle_I(KindEnd_)%nVar = nVarParticle
    call allocate_particles
  end subroutine init_particle_line

  !==========================================================================

  subroutine extract_particle_line(nFieldLineIn, XyzInit_DI)
    ! extract nFieldLineIn magnetic field lines starting at XyzInit_DI;
    ! the whole field lines are extracted, i.e. they are traced forward
    ! and backward up until it reaches boundaries of the domain
    !------------------------------------------------------------------------
    integer, intent(in):: nFieldLineIn
    real,    intent(in):: XyzInit_DI(MaxDim, nFieldLineIn)

    integer :: nFieldLineNew ! number of new field lines initialized
    integer :: iFieldLine ! loop variable
    integer :: iBlock
    integer :: nParticleRegOld ! number of already existing regular particles
    integer :: iParticle ! loop variable
    logical :: IsCompleted_I(nFieldLineIn)

    ! direction of the magnetic field
    real:: Dir_D(MaxDim)

    ! direction of tracing: -1 -> backward, +1 -> forward
    integer:: iDirTrace


    ! parameters of end particles
    real,    pointer:: StateEnd_VI(:,:)
    integer, pointer:: iBlockEnd_I(:)
    ! parameters of regular particles
    real,    pointer:: StateReg_VI(:,:)
    integer, pointer:: iBlockReg_I(:)

    ! TEMPORARY: STORE WHETHER PARTICLES HAVE BEEN INITIALIZED
    logical, save:: IsInitialized = .false.
    !------------------------------------------------------------------------
    if(.not. IsInitialized)then
       call init_particle_line
       IsInitialized = .true.
    end if
    ! set a pointers to parameters of end particles
    nullify(StateEnd_VI); nullify(iBlockEnd_I)
    StateEnd_VI => Particle_I(KindEnd_)%State_VI
    iBlockEnd_I => Particle_I(KindEnd_)%iBlock_I

    ! set a pointers to parameters of regular particles
    nullify(StateReg_VI); nullify(iBlockReg_I)

    StateReg_VI => Particle_I(KindReg_)%State_VI
    iBlockReg_I => Particle_I(KindReg_)%iBlock_I

    !\
    ! Trace field lines
    !/
    ! initialize field lines
    Particle_I(KindEnd_)%nParticle = 0
    nFieldLineNew = 0
    nParticleRegOld = Particle_I(KindReg_)%nParticle
    do iFieldLine = 1, nFieldLineIn
       call start_line(XyzInit_DI(:, iFieldLine))
    end do
    call copy_end_to_regular

    ! trace in both directions
    do iDirTrace = -1, 1, 2
       TRACE: do
          ! copy last known coordinates to Auxilary 
          StateEnd_VI(AuxX_:AuxZ_,1:Particle_I(KindEnd_)%nParticle) = &
               StateEnd_VI(x_:z_, 1:Particle_I(KindEnd_)%nParticle)
          !\
          ! predictor step
          !/
          do iParticle = 1, Particle_I(KindEnd_)%nParticle
             ! get the direction of the magnetic field at original location
             call get_b_dir(StateEnd_VI(x_:z_, iParticle), Dir_D)
             ! find the step size 
             StateEnd_VI(Aux_, iParticle) = &
                  0.1*SQRT(sum(CellSize_DB(1:nDim,iBlockEnd_I(iParticle))**2))
             ! get middle location
             StateEnd_VI(x_:z_, iParticle) = StateEnd_VI(x_:z_, iParticle) + &
                  0.5 * iDirTrace * StateEnd_VI(Aux_, iParticle) * Dir_D
          end do
          !\
          ! Message pass: some particles may have moved to different procs
          !/
          call message_pass_particles(KindEnd_)

          !\
          ! remove particles that went outside of the domain
          !/
          if(Body1)then ! check if there is an inner body
             call check_inner_boundary_particle_line(KindEnd_)
          end if
          call remove_undefined_particles(KindEnd_)

          ! check if all field lines have been completed
          if(Particle_I(KindEnd_)%nParticle == 0) EXIT TRACE

          !\
          ! Corrector step
          !/
          do iParticle = 1, Particle_I(KindEnd_)%nParticle
             ! get the direction of the magnetic field in the middle
             call get_b_dir(StateEnd_VI(x_:z_, iParticle), Dir_D)
             ! get final location
             StateEnd_VI(x_:z_,iParticle)=StateEnd_VI(AuxX_:AuxZ_,iParticle)+&
                  iDirTrace * StateEnd_VI(Aux_, iParticle) * Dir_D
          end do
          !\
          ! Message pass: some particles may have moved to different procs
          !/
          call message_pass_particles(KindEnd_)

          !\
          ! remove particles that went outside of the domain
          !/
          if(Body1)then ! check if there is an inner body
             call check_inner_boundary_particle_line(KindEnd_)
          end if
          call remove_undefined_particles(KindEnd_)

          ! check if all field lines have been completed
          if(Particle_I(KindEnd_)%nParticle == 0) EXIT TRACE

          ! increase particle index & copy to regular
          StateEnd_VI(Index_,1:Particle_I(KindEnd_)%nParticle) = &
               StateEnd_VI(Index_,1:Particle_I(KindEnd_)%nParticle) + &
               iDirTrace
          call copy_end_to_regular

       end do TRACE

       ! if just finished backward tracing => return to initial particles
       if(iDirTrace < 0)then
          ! first fix indices of particles along field lines
          do iParticle = nParticleRegOld+1, Particle_I(KindReg_)%nParticle
             iFieldLine = nint(StateReg_VI(fl_, iParticle))
             StateReg_VI(Index_, iParticle) = &
                  StateReg_VI(Index_, iParticle) + &
                  nParticleFieldLine_I(iFieldLine)
          end do

          ! now copy initial points to KindEnd_:
          ! the initial particles are currently right after the particles,
          ! that were in the list before the start of this subroutine,
          ! i.e. occupy positions from (nParticleRegOld+1)
          StateEnd_VI(:, 1:nFieldLineNew) = &
               StateReg_VI(:, nParticleRegOld+1:nParticleRegOld+nFieldLineNew)
          iBlockEnd_I(   1:nFieldLineNew) = &
               iBlockReg_I(   nParticleRegOld+1:nParticleRegOld+nFieldLineNew)
          Particle_I(KindEnd_)%nParticle = nFieldLineNew
       end if
    end do

    ! finally, sort particle list with respect to indices along field lines
    call sort_particle_line

  contains

    subroutine start_line(XyzStart_D)
      real,    intent(in) :: XyzStart_D(MaxDim)

      real   :: Coord_D(MaxDim) ! generalized coordinates
      integer:: iProcOut, iBlockOut
      ! variables to call check_interpolate_amr_gc
      logical:: IsPossible, IsBoundary
      !----------------------------------------------------------------------
      ! find block and processor suitable for interpolation
      call check_interpolate_amr_gc(XyzStart_D, &
           1, IsPossible, & ! input block ID & output IsPossible don't matter
           iProcOut, iBlockOut, IsBoundary)

      ! check whether point is outside of the domain
      if(IsBoundary) RETURN

      !\
      ! Assign particle to an appropriate processor
      !/
      if(iProc /= iProcOut) RETURN

      Particle_I(KindEnd_)%nParticle = &
           Particle_I(KindEnd_)%nParticle + 1         
      nFieldLineNew = nFieldLineNew + 1
      nFieldLine    = nFieldLine    + 1

      StateEnd_VI(x_:z_, Particle_I(KindEnd_)%nParticle) = XyzStart_D
      StateEnd_VI(fl_,   Particle_I(KindEnd_)%nParticle) = nFieldLine

      ! set index for initial particle to be 0 for now, it's fixed later
      StateEnd_VI(Index_,Particle_I(KindEnd_)%nParticle) = 0
      iBlockEnd_I(Particle_I(KindEnd_)%nParticle) = iBlockOut
    end subroutine start_line
    !========================================================================
    subroutine copy_end_to_regular
      ! copies indicated variables of known end particles to regular particles
      integer, parameter:: iVarCopy_I(5) = (/x_, y_, z_, fl_, Index_/)
      !----------------------------------------------------------------------
      StateReg_VI(iVarCopy_I,&
           Particle_I(KindReg_)%nParticle+1:&
           Particle_I(KindReg_)%nParticle+Particle_I(KindEnd_)%nParticle) =&
           StateEnd_VI(iVarCopy_I, 1:Particle_I(KindEnd_)%nParticle)
      iBlockReg_I(&
           Particle_I(KindReg_)%nParticle+1:&
           Particle_I(KindReg_)%nParticle+Particle_I(KindEnd_)%nParticle) =&
           iBlockEnd_I(1:Particle_I(KindEnd_)%nParticle)
      Particle_I(KindReg_)%nParticle = &
           Particle_I(KindReg_)%nParticle + &
           Particle_I(KindEnd_)%nParticle
      ! also update the counter of particles per field line
      do iParticle = 1, Particle_I(KindEnd_)%nParticle
         iFieldLine = nint(StateEnd_VI(fl_,iParticle))
         nParticleFieldLine_I(iFieldLine) = &
              nParticleFieldLine_I(iFieldLine) + 1
      end do
    end subroutine copy_end_to_regular
    !========================================================================
    subroutine get_b_dir(Xyz_D, Dir_D)
      ! returns the direction of magnetic field 
      ! as well as the block used for interpolation
      real,    intent(in) :: Xyz_D(MaxDim)
      real,    intent(out):: Dir_D(MaxDim)

      ! magnetic field
      real   :: B_D(MaxDim) = 0.0
      ! inteprolation data: number of cells, cell indices, weights
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim) = 1
      integer:: iBlock
      !----------------------------------------------------------------------
      ! get potential part of the magnetic field at the current location
      call get_b0(Xyz_D, B_D)
      ! get the remaining part of the maagnetic field
      call interpolate_grid_amr_gc(Xyz_D, nCell, iCell_II, Weight_I)
      ! interpolate magnetic field value
      do iCell = 1, nCell
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         iBlock      = iCell_II(0, iCell)
         B_D = B_D + &
              State_VGB(Bx_:Bz_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
      end do
      ! normalize vector to unity
      Dir_D(1:nDim) = B_D(1:nDim) / sum(B_D(1:nDim)**2)**0.5
    end subroutine get_b_dir

  end subroutine extract_particle_line

  !========================================================================

  subroutine check_inner_boundary_particle_line(iKind)
    use BATL_lib, ONLY: Unset_
    ! check whether particles have crossed an inner boundary
    ! if so => set particle's block to Unset_
    integer, intent(in):: iKind

    integer:: iParticle ! loop variable
    !----------------------------------------------------------------------
    do iParticle = 1, Particle_I(iKind)%nParticle
       if(sum(Particle_I(iKind)%State_VI(x_:z_,iParticle)**2) < rBody*rBody)&
            Particle_I(iKind)%iBlock_I(iParticle) = Unset_
    end do
  end subroutine check_inner_boundary_particle_line

  !========================================================================

  subroutine sort_particle_line
    ! the subroutine sorts particles in the increasing order 
    integer, allocatable:: iOrder_I(:)
    integer:: iFieldLine, iParticle
    !----------------------------------------------------------------------
    allocate(iOrder_I(Particle_I(KindReg_)%nParticle))
    ! compute offsets (cumulutive sum) to different field lines
    iParticleFieldLineOffset_I(1) = 0
    do iFieldLine = 2, nFieldLine
       iParticleFieldLineOffset_I(iFieldLine) = &
            iParticleFieldLineOffset_I(iFieldLine-1) + &
            nParticleFieldLine_I(iFieldLine-1)
    end do
    ! find the order of particles
    do iParticle = 1, Particle_I(KindReg_)%nParticle
       iFieldLine = nint(Particle_I(KindReg_)%State_VI(fl_, iParticle))
       iOrder_I(iParticle) = &
            iParticleFieldLineOffset_I(iFieldLine) + &
            nint(Particle_I(KindReg_)%State_VI(Index_, iParticle))
    end do
    ! rearrange particles
    Particle_I(KindReg_)%State_VI(:, iOrder_I) = &
         Particle_I(KindReg_)%State_VI(:, 1:Particle_I(KindReg_)%nParticle)
    Particle_I(KindReg_)%iBlock_I(   iOrder_I) = &
         Particle_I(KindReg_)%iBlock_I(   1:Particle_I(KindReg_)%nParticle)

    deallocate(iOrder_I)

  end subroutine sort_particle_line

  !========================================================================
  
  subroutine get_particle_data(NameVar, DataOut_VI, nDataVar, nParticle)
    ! the subroutine gets variables specified in the string StringVar
    ! and writes them into DataOut_VI
    character(len=*),   intent(in) :: NameVar
    real,  pointer,     intent(out):: DataOut_VI(:,:)
    integer,            intent(out):: nDataVar
    integer,            intent(out):: nParticle

    ! mask for returning variables
    logical:: DoReturnVar_V(nVarParticle)
    ! string used for processing NameVar
    character(len=100):: NameVarProc
    integer:: nFullLen
    integer:: iParticle ! loop variable
    !----------------------------------------------------------------------
    ! first, determine which variables will be returned
    DoReturnVar_V = .false.
    NameVarProc =  trim(NameVar)
    nFullLen = len(NameVarProc)

    if(is_var('xx'))&
         DoReturnVar_V(x_) = .true.
    if(is_var('yy'))&
         DoReturnVar_V(y_) = .true.
    if(is_var('zz'))&
         DoReturnVar_V(z_) = .true.
    if(is_var('fl'))&
         DoReturnVar_V(fl_) = .true.
    if(is_var('id'))&
         DoReturnVar_V(Index_) = .true.

    ! check if there is an invalid variable requested
    if(len_trim(NameVarProc) > 0)&
         call CON_stop("Requiesting an invalid variable for particle line")

    ! return data
    if(associated(DataOut_VI)) deallocate(DataOut_VI)
    nDataVar  = count(DoReturnVar_V)
    nParticle = Particle_I(KindReg_)%nParticle
    allocate( DataOut_VI(nDataVar, nParticle) )
    do iParticle = 1, nParticle
       DataOut_VI(:, iParticle) = PACK(&
            Particle_I(KindReg_)%State_VI(:,iParticle),&
            MASK = DoReturnVar_V)
    end do

  contains
    
    function is_var(Name) result(IsPresent)
      ! finds if var Name is present
      ! if yes => remove it from NameVarProc
      character(len=*), intent(in):: Name
      logical:: IsPresent
      integer:: i, nLen
      !--------------------------------------------------------------------
      i = index(NameVarProc, Name)
      if(i==0)then
         IsPresent = .false.
         RETURN
      end if
      nLen = len(Name)
      NameVarProc = trim(NameVarProc(1:i-1) // NameVarProc(i+nLen+1:nFullLen))
      IsPresent = .true.
    end function is_var

  end subroutine get_particle_data

end module ModParticleFieldLine
