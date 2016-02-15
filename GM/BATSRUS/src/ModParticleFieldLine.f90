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
  use ModMain, ONLY: Body1
  use ModPhysics, ONLY: rBody

  implicit none
  SAVE

  private !except

  public:: read_particle_line_param
  public:: init_particle_line
  public:: extract_particle_line
  public:: get_particle_data


  ! kinds of particles used to generate a magnetic field line
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
       Aux_   = 7
  
  ! indices of a particle
  integer, parameter:: &
       ! field line this particle lays on
       fl_    = 1, &
       ! index of the particle along this field line
       Index_ = 2

  ! mode of spatial step at the field line extraction
  real   :: SpaceStepMin = 0.0
  real   :: SpaceStepMax = HUGE(SpaceStepMin)

  ! number of variables in the state vector
  integer, parameter:: nVarParticle = 7

  ! number of indices
  integer, parameter:: nIndexParticle = 2

  ! maximum allowed number of field lines
  integer, parameter :: nFieldLineMax = 1000

  ! number of active field lines
  integer:: nFieldLine = 0

  ! keep track of number of particles per field line
  integer :: nParticleFieldLine_I(nFieldLineMax) = 0
!  integer :: iParticleFieldLineOffset_I(nFieldLineMax) = 0

  !\
  ! initialization related info
  !/
  integer:: nLineInit
  real, allocatable:: XyzLineInit_DI(:,:)

contains

  subroutine read_particle_line_param(NameCommand)

    use ModMain,      ONLY: UseParticles
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: &
         NameSub = 'ModParticleFieldLine::read_particle_line_param'

    character(len=100) :: StringInitMode
    integer:: iLine, iDim ! loop variables
    !------------------------------------------------------------------------
    select case(NameCommand)
    case("#PARTICLELINE")
       call read_var('UseParticles', UseParticles)
       if(UseParticles)then
          ! read min and max values for space step
          ! based on their values space step may be
          ! - both negative => automatic (defined by grid resolution)
          ! - otherwise     => restricted by min and/or max
          !                    whichever is positive
          call read_var('SpaceStepMin', SpaceStepMin)
          call read_var('SpaceStepMax', SpaceStepMax)
          ! negative value means "undefined"
          if(SpaceStepMin < 0) SpaceStepMin = 0
          if(SpaceStepMax < 0) SpaceStepMax = HUGE(SpaceStepMin)
          !--------------------------------------------------------------
          call read_var('InitMode', StringInitMode, IsLowerCase=.true.)
          ! Initialization modes:
          ! - preset: starting points are set from PARAM.in file
          if(index(StringInitMode, 'preset') > 0)then
             call read_var('nLineInit', nLineInit)
             if(nLineInit <= 0)&
                  call stop_mpi(NameSub // &
                  ": invalid number of initialized particle lines")
             allocate(XyzLineInit_DI(MaxDim, nLineInit))
             do iLine = 1, nLineInit; do iDim = 1, MaxDim
                call read_var('XyzLineInit_DI', XyzLineInit_DI(iDim, iLine))
             end do; end do
          elseif(index(StringInitMode, 'import') > 0)then
             ! do nothing: particles are imported from elsewhere
          else
             call stop_mpi(NameSub //": unknown initialization mode")
          end if
       end if
    end select

  end subroutine read_particle_line_param

  !==========================================================================

  subroutine init_particle_line
    ! allocate containers for particles
    integer:: iLine ! loop variable
    !------------------------------------------------------------------------
    Particle_I(KindReg_)%nParticleMax = 10000 * nFieldLineMax
    Particle_I(KindEnd_)%nParticleMax = nFieldLineMax
    Particle_I(KindReg_)%nVar   = nVarParticle
    Particle_I(KindEnd_)%nVar   = nVarParticle
    Particle_I(KindReg_)%nIndex = nIndexParticle
    Particle_I(KindEnd_)%nIndex = nIndexParticle
    call allocate_particles
    if(nLineInit > 0)&
         ! extract initial field lines
         call extract_particle_line(nLineInit, XyzLineInit_DI)
  end subroutine init_particle_line

  !==========================================================================

  subroutine extract_particle_line(nFieldLineIn, XyzStart_DI)
    ! extract nFieldLineIn magnetic field lines starting at XyzStart_DI;
    ! the whole field lines are extracted, i.e. they are traced forward
    ! and backward up until it reaches boundaries of the domain
    !------------------------------------------------------------------------
    integer, intent(in):: nFieldLineIn
    real,    intent(in):: XyzStart_DI(MaxDim, nFieldLineIn)

    integer :: nLineThisProc ! number of new field lines initialized
    integer :: iFieldLine    ! loop variable
    integer :: iBlock
    integer :: nParticleOld  ! number of already existing regular particles
    integer :: iParticle     ! loop variable

    ! direction of the magnetic field
    real:: Dir_D(MaxDim)

    ! direction of tracing: -1 -> backward, +1 -> forward
    integer:: iDirTrace

    ! parameters of end particles
    real,    pointer:: StateEnd_VI(:,:)
    integer, pointer:: iIndexEnd_II(:,:)
    ! parameters of regular particles
    real,    pointer:: StateReg_VI(:,:)
    integer, pointer:: iIndexReg_II(:,:)
    !------------------------------------------------------------------------
    ! set a pointers to parameters of end particles
    nullify(StateEnd_VI); nullify(iIndexEnd_II)
    StateEnd_VI  => Particle_I(KindEnd_)%State_VI
    iIndexEnd_II => Particle_I(KindEnd_)%iIndex_II

    ! set a pointers to parameters of regular particles
    nullify(StateReg_VI); nullify(iIndexReg_II)
    StateReg_VI  => Particle_I(KindReg_)%State_VI
    iIndexReg_II => Particle_I(KindReg_)%iIndex_II

    !\
    ! Trace field lines
    !/
    ! initialize field lines
    Particle_I(KindEnd_)%nParticle = 0
    nLineThisProc = 0
    nParticleOld  = Particle_I(KindReg_)%nParticle
    nFieldLine    = nFieldLine + nFieldLineIn
    if(nFieldLine > nFieldLineMax)&
         call CON_stop('Limit for number of particle field lines exceeded')
    do iFieldLine = 1, nFieldLineIn
       call start_line(XyzStart_DI(:, iFieldLine), iFieldLine)
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
             call get_b_dir(Xyz_D = StateEnd_VI(x_:z_, iParticle),&
                  iBlock=iIndexEnd_II(0,iParticle),&
                  Dir_D = Dir_D)
             ! find the step size
             StateEnd_VI(Aux_, iParticle) = &
                  MIN(SpaceStepMax, MAX(SpaceStepMin,&
                  0.1*SQRT(&
                  sum(CellSize_DB(1:nDim,iIndexEnd_II(0,iParticle))**2)&
                  )))

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
          if(is_complete()) EXIT TRACE

          !\
          ! Corrector step
          !/
          do iParticle = 1, Particle_I(KindEnd_)%nParticle
             ! get the direction of the magnetic field in the middle
             call get_b_dir(Xyz_D=StateEnd_VI(x_:z_, iParticle),&
                  iBlock=iIndexEnd_II(0,iParticle),&
                  Dir_D=Dir_D)
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
          if(is_complete()) EXIT TRACE

          ! increase particle index & copy to regular
          iIndexEnd_II(Index_,1:Particle_I(KindEnd_)%nParticle) = &
               iIndexEnd_II(Index_,1:Particle_I(KindEnd_)%nParticle) + &
               iDirTrace
          call copy_end_to_regular

       end do TRACE

       ! if just finished backward tracing => return to initial particles
       if(iDirTrace < 0)then
          ! copy initial points to KindEnd_:
          ! the initial particles are currently right after the particles,
          ! that were in the list before the start of this subroutine,
          ! i.e. occupy positions from (nParticleOld+1)
          StateEnd_VI( :, 1:nLineThisProc) = &
               StateReg_VI( :, nParticleOld+1:nParticleOld+nLineThisProc)
          iIndexEnd_II(:, 1:nLineThisProc) = &
               iIndexReg_II(:, nParticleOld+1:nParticleOld+nLineThisProc)
          Particle_I(KindEnd_)%nParticle = nLineThisProc
       end if
    end do

  contains

    subroutine start_line(XyzStart_D, iFieldLineIndex)
      real,    intent(in) :: XyzStart_D(MaxDim)
      integer, intent(in) :: iFieldLineIndex

      real   :: Coord_D(MaxDim) ! generalized coordinates
      integer:: iProcOut, iBlockOut
      ! container for error message 
      character(len=100):: StringError

      character(len=*), parameter :: &
           NameSub = 'ModParticleFieldLine::extract_particle_line::start_line'
      !----------------------------------------------------------------------
      ! find block and processor suitable for interpolation
      call check_interpolate_amr_gc(XyzStart_D, &
           1, & ! input block ID doesn't matter
           iProcOut, iBlockOut)

      ! check whether point is outside of the domain
      if(iProcOut < 0)then
         write(StringError,'(a,es15.6,a, es15.6, a, es15.6)') &
              "Start point for a field line is outside of the domain: X = ",&
              XyzStart_D(1), " Y = ", XyzStart_D(2), " Z = ", XyzStart_D(3) 
         call stop_mpi(NameSub //": "//StringError)
      end if

      !\
      ! Assign particle to an appropriate processor
      !/
      if(iProc /= iProcOut) RETURN

      Particle_I(KindEnd_)%nParticle = &
           Particle_I(KindEnd_)%nParticle + 1         
      nLineThisProc = nLineThisProc + 1

      StateEnd_VI(x_:z_, Particle_I(KindEnd_)%nParticle) = XyzStart_D
      iIndexEnd_II(fl_,  Particle_I(KindEnd_)%nParticle) = iFieldLineIndex

      ! set index for initial particle to be 0
      iIndexEnd_II(Index_,Particle_I(KindEnd_)%nParticle) = 0
      iIndexEnd_II(0,     Particle_I(KindEnd_)%nParticle) = iBlockOut
    end subroutine start_line
    !========================================================================
    function is_complete() result(IsCompleteOut)
      use ModMpi
      use BATL_mpi, ONLY: iComm
      ! returns true if tracing is complete for all field line on all procs
      logical:: IsComplete, IsCompleteOut
      integer:: iError
      !----------------------------------------------------------------------
      if(Particle_I(KindEnd_)%nParticle < 0)STOP
      IsComplete = Particle_I(KindEnd_)%nParticle == 0
      ! reduce IsComplete variable from all processors
      call MPI_Allreduce(IsComplete, IsCompleteOut, 1, &
           MPI_LOGICAL, MPI_LAND, iComm, iError)               
    end function is_complete
    !========================================================================
    subroutine copy_end_to_regular
      ! copies indicated variables of known end particles to regular particles
      integer, parameter:: iVarCopy_I(3)   = (/x_, y_, z_/)
      integer, parameter:: iIndexCopy_I(3) = (/0, fl_, Index_/)
      !----------------------------------------------------------------------
      StateReg_VI(iVarCopy_I,&
           Particle_I(KindReg_)%nParticle+1:&
           Particle_I(KindReg_)%nParticle+Particle_I(KindEnd_)%nParticle) =&
           StateEnd_VI(iVarCopy_I, 1:Particle_I(KindEnd_)%nParticle)
      iIndexReg_II(iIndexCopy_I,&
           Particle_I(KindReg_)%nParticle+1:&
           Particle_I(KindReg_)%nParticle+Particle_I(KindEnd_)%nParticle) =&
           iIndexEnd_II(iIndexCopy_I, 1:Particle_I(KindEnd_)%nParticle)
      Particle_I(KindReg_)%nParticle = &
           Particle_I(KindReg_)%nParticle + &
           Particle_I(KindEnd_)%nParticle
      ! also update the counter of particles per field line
      do iParticle = 1, Particle_I(KindEnd_)%nParticle
         iFieldLine = iIndexEnd_II(fl_,iParticle)
         nParticleFieldLine_I(iFieldLine) = &
              nParticleFieldLine_I(iFieldLine) + 1
      end do
    end subroutine copy_end_to_regular
    !========================================================================
    subroutine get_b_dir(Xyz_D, iBlock, Dir_D)
      use ModMain, ONLY: UseB0
      use ModB0, ONLY: get_b0
      ! returns the direction of magnetic field 
      ! as well as the block used for interpolation
      real,    intent(in) :: Xyz_D(MaxDim)
      integer, intent(in) :: iBlock
      real,    intent(out):: Dir_D(MaxDim)

      ! magnetic field
      real   :: B_D(MaxDim) = 0.0
      ! inteprolation data: number of cells, cell indices, weights
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim) = 1
      !----------------------------------------------------------------------
      Dir_D = 0; B_D = 0
      ! get potential part of the magnetic field at the current location
      if(UseB0)call get_b0(Xyz_D, B_D)
      ! get the remaining part of the magnetic field
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
      ! interpolate magnetic field value
      do iCell = 1, nCell
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
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
            Particle_I(iKind)%iIndex_II(0, iParticle) = Unset_
    end do
  end subroutine check_inner_boundary_particle_line

  !========================================================================

  subroutine get_particle_data(NameVar, DataOut_VI,   nDataOut, nParticle)
    ! the subroutine gets variables specified in the string StringVar
    ! and writes them into DataOut_VI
    character(len=*),   intent(in) :: NameVar
    real,  pointer,     intent(out):: DataOut_VI(:,:)
    integer,            intent(out):: nDataOut
    integer,            intent(out):: nParticle

    ! mask for returning variables
    logical:: DoReturnVar_V(nVarParticle)
    logical:: DoReturnIndex_I(0:nIndexParticle)
    ! string used for processing NameVar
    character(len=100):: NameVarProc
    integer:: nFullLen
    integer:: nVarOut, nIndexOut
    integer:: iParticle ! loop variable
    !----------------------------------------------------------------------
    ! first, determine which variables will be returned
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.
    NameVarProc =  trim(NameVar)
    nFullLen = len(NameVarProc)

    if(is_var('xx'))&
         DoReturnVar_V(x_) = .true.
    if(is_var('yy'))&
         DoReturnVar_V(y_) = .true.
    if(is_var('zz'))&
         DoReturnVar_V(z_) = .true.
    if(is_var('fl'))&
         DoReturnIndex_I(fl_) = .true.
    if(is_var('id'))&
         DoReturnIndex_I(Index_) = .true.

    ! check if there is an invalid variable requested
    if(len_trim(NameVarProc) > 0)&
         call CON_stop("Requiesting an invalid variable for particle line")

    ! return data
    if(associated(DataOut_VI)) deallocate(DataOut_VI)
    nVarOut   = count(DoReturnVar_V)
    nIndexOut = count(DoReturnIndex_I)
    nDataOut  = nIndexOut + nVarOut
    nParticle = Particle_I(KindReg_)%nParticle
    allocate( DataOut_VI(nDataOut, nParticle) )
    do iParticle = 1, nParticle
       DataOut_VI(1:nVarOut,         iParticle) = PACK(&
            Particle_I(KindReg_)%State_VI( :,iParticle),&
            MASK = DoReturnVar_V)
       DataOut_VI(nVarOut+1:nDataOut,iParticle) = PACK(&
            Particle_I(KindReg_)%iIndex_II(:,iParticle),&
            MASK = DoReturnIndex_I)
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
