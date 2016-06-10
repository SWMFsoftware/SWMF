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
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_
  use ModMain, ONLY: Body1, NameThisComp
  use ModPhysics, ONLY: rBody

  implicit none
  SAVE

  private !except

  public:: read_particle_line_param
  public:: init_particle_line
  public:: extract_particle_line
  public:: advect_particle_line
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
       Aux_  = 7
  
  ! indices of a particle
  integer, parameter:: &
       ! field line this particle lays on
       fl_    = 1, &
       ! index of the particle along this field line
       Index_ = 2

  ! spatial step limits at the field line extraction
  real   :: SpaceStepMin = 0.0
  real   :: SpaceStepMax = HUGE(SpaceStepMin)

  ! ordering mode of particles along field lines
  integer:: iOrderMode = -1
  integer, parameter:: Field_ = 0, Radius_  = 1

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
    character(len=100) :: StringOrderMode
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
          !--------------------------------------------------------------
          call read_var('OrderMode', StringOrderMode, IsLowerCase=.true.)
          ! Ordering modes, particle index increases with:
          ! - radius:  distance from the center
          ! - field: along the direction of the magnetic field
          if(    index(StringOrderMode, 'field') > 0)then
             iOrderMode = Field_
          elseif(index(StringOrderMode, 'radius' ) > 0)then
             iOrderMode = Radius_
          else
             call CON_stop(NameThisComp//':'//NameSub //": unknown ordering mode")
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

  subroutine extract_particle_line(nFieldLineIn, XyzStart_DI, iTraceModeIn)
    ! extract nFieldLineIn magnetic field lines starting at XyzStart_DI;
    ! the whole field lines are extracted, i.e. they are traced forward
    ! and backward up until it reaches boundaries of the domain
    !------------------------------------------------------------------------
    integer,           intent(in):: nFieldLineIn
    real,              intent(in):: XyzStart_DI(MaxDim, nFieldLineIn)
    integer, optional, intent(in):: iTraceModeIn

    integer :: nLineThisProc ! number of new field lines initialized
    integer :: iFieldLine    ! loop variable
    integer :: iBlock
    integer :: nParticleOld  ! number of already existing regular particles
    integer :: iParticle     ! loop variable

    ! direction of the magnetic field
    real:: Dir1_D(MaxDim), Dir2_D(MaxDim)

    ! direction of tracing: -1 -> backward, +1 -> forward
    integer:: iDirTrace

    ! alignment of particle numbering with direction of B field:
    ! ( -1 -> reversed, +1 -> aligned)
    ! for radial ordering when magnetic field general direction may be inward
    integer:: iAlignment_I(nFieldLine+1: nFieldLine+nFieldLineIn)

    ! mode of tracing (see description below)
    integer:: iTraceMode

    ! starting point of a field line 
    ! (to satisfy intent INOUT for check_interpolate_amr_gc)
    real:: XyzStart_D(MaxDim)

    ! parameters of end particles
    real,    pointer:: StateEnd_VI(:,:)
    integer, pointer:: iIndexEnd_II(:,:)
    ! parameters of regular particles
    real,    pointer:: StateReg_VI(:,:)
    integer, pointer:: iIndexReg_II(:,:)

    ! constant in Ralston's method
    real, parameter:: cTwoThird = 2.0 / 3.0

    character(len=*), parameter:: NameSub='extract_particle_line'
    !------------------------------------------------------------------------
    ! set a pointers to parameters of end particles
    StateEnd_VI  => Particle_I(KindEnd_)%State_VI
    iIndexEnd_II => Particle_I(KindEnd_)%iIndex_II

    ! set a pointers to parameters of regular particles
    StateReg_VI  => Particle_I(KindReg_)%State_VI
    iIndexReg_II => Particle_I(KindReg_)%iIndex_II

    !\
    ! Trace field lines
    !/
    ! initialize field lines
    Particle_I(KindEnd_)%nParticle = 0
    nLineThisProc = 0
    nParticleOld  = Particle_I(KindReg_)%nParticle
    do iFieldLine = 1, nFieldLineIn
       XyzStart_D = XyzStart_DI(:, iFieldLine) 
       call start_line(XyzStart_D, nFieldLine + iFieldLine)
    end do
    call get_alignment()

    nFieldLine    = nFieldLine + nFieldLineIn
    if(nFieldLine > nFieldLineMax)&
         call CON_stop(NameThisComp//':'//NameSub//&
         ': Limit for number of particle field lines exceeded')
    call copy_end_to_regular

    ! check if trace mode is specified
    if(present(iTraceModeIn))then
       if(abs(iTraceModeIn) > 1)&
            call CON_stop(NameThisComp//':'//NameSub//': incorrect tracing mode')
       iTraceMode = iTraceModeIn
    else
       iTraceMode = 0
    end if

    ! tracing modes are the following:,
    ! -1 -> trace lines backward only
    !       do iDirTrace = -1, -1, 2
    !  0 -> trace lines in both directions
    !       do iDirTrace = -1,  1, 2
    ! +1 -> trace lines forward only
    !       do iDirTracs =  1,  1, 2
    do iDirTrace = 2*max(iTraceMode,0)-1, 2*min(iTraceMode,0)+1, 2
       TRACE: do
          ! copy last known coordinates to Auxilary 
          StateEnd_VI(AuxX_:AuxZ_,1:Particle_I(KindEnd_)%nParticle) = &
               StateEnd_VI(x_:z_, 1:Particle_I(KindEnd_)%nParticle)
          !\
          ! First stage of Ralston's method
          !/
          do iParticle = 1, Particle_I(KindEnd_)%nParticle
             ! get the direction of the magnetic field at original location
             call get_b_dir(&
                  Xyz_D = StateEnd_VI(x_:z_, iParticle),&
                  iBlock=iIndexEnd_II(0,iParticle),&
                  Dir_D = Dir1_D)
             ! find the step size
             StateEnd_VI(Aux_, iParticle) = &
                  MIN(SpaceStepMax, MAX(SpaceStepMin,&
                  0.1*SQRT(&
                  sum(CellSize_DB(1:nDim,iIndexEnd_II(0,iParticle))**2)&
                  )))

             ! get middle location
             StateEnd_VI(x_:z_, iParticle) = StateEnd_VI(x_:z_, iParticle) + &
                  iDirTrace * StateEnd_VI(Aux_, iParticle) * &
                  cTwoThird * Dir1_D * &
                  iAlignment_I(iIndexEnd_II(fl_, iParticle))
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
          ! Second stage of Ralston's method
          !/
          do iParticle = 1, Particle_I(KindEnd_)%nParticle
             ! get the direction of the magnetic field in the middle
             call get_b_dir(&
                  Xyz_D = StateEnd_VI(x_:z_, iParticle),&
                  iBlock=iIndexEnd_II(0,iParticle),&
                  Dir_D = Dir2_D)
             ! get final location
             StateEnd_VI(x_:z_,iParticle)=StateEnd_VI(AuxX_:AuxZ_,iParticle)+&
                  iDirTrace * StateEnd_VI(Aux_, iParticle) * &
                  (0.25 * Dir1_D + 0.75 * Dir2_D) * &
                  iAlignment_I(iIndexEnd_II(fl_, iParticle))
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

       ! if just finished backward tracing and need to trace in both dirs
       ! => return to initial particles
       if(iDirTrace < 0 .and. iTraceMode == 0)then
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
      real,    intent(inout):: XyzStart_D(MaxDim)
      integer, intent(in)   :: iFieldLineIndex

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
    subroutine get_alignment()
      use ModMpi
      use BATL_mpi, ONLY: iComm
      ! determine alignment of particle indexing with direction 
      ! of the magnetic field
      integer:: iParticle, iError
      real:: Dir_D(MaxDim)
      !------------------------------------------------------------------------
      if(iOrderMode == Field_)then
         iAlignment_I = 1
         RETURN
      end if

      iAlignment_I = 0
      do iParticle = 1, Particle_I(KindEnd_)%nParticle
         call get_b_dir(&
              Xyz_D = StateEnd_VI(x_:z_, iParticle),&
              iBlock=iIndexEnd_II(0,iParticle),&
              Dir_D = Dir_D)
         iAlignment_I(iIndexEnd_II(fl_, iParticle)) = &
              nint( SIGN(1.0, sum(Dir_D*StateEnd_VI(x_:z_,iParticle))) )
      end do
      call MPI_Allreduce(MPI_IN_PLACE, iAlignment_I, nFieldLineIn, &
           MPI_INTEGER, MPI_SUM, iComm, iError)
    end subroutine get_alignment
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
      character(len=200):: StringError
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
      if(all(B_D==0))then
         write(StringError,'(a,es15.6,a,es15.6,a,es15.6)') &
              NameThisComp//':'//NameSub//&
              ': trying to extract line at region with zero magnetic field'//&
              ' at location X = ', &
              Xyz_D(1), ' Y = ', Xyz_D(2), ' Z = ', Xyz_D(3) 
         call CON_stop(StringError)
      end if
      ! normalize vector to unity
      Dir_D(1:nDim) = B_D(1:nDim) / sum(B_D(1:nDim)**2)**0.5
    end subroutine get_b_dir

  end subroutine extract_particle_line

  !========================================================================

  subroutine advect_particle_line
    ! advect particles with the local plasma velocity
    use ModMain, ONLY: time_accurate, Dt, nStage
    ! parameters of particles
    real,    pointer:: State_VI(:,:)
    integer, pointer:: iIndex_II(:,:)

    ! local plasma velocity
    real:: VLocal_D(MaxDim)
    !loop variable
    integer:: iParticle

    character(len=*),parameter::NameSub='advect_particle_line'
    !----------------------------------------------------------------------
    ! particles can be advected only in the time accurate run
    if(.not.time_accurate) &
         RETURN

    ! set a pointers to parameters of particles
    State_VI  => Particle_I(KindReg_)%State_VI
    iIndex_II => Particle_I(KindReg_)%iIndex_II
    !\
    ! go over the list of particles and advect them
    !/
    do iParticle = 1, Particle_I(KindReg_)%nParticle
       !\
       ! simple advection (Euler integration)
       !/ 
       ! get the local velocity
       call get_v(&
            Xyz_D = State_VI(x_:z_, iParticle),&
            iBlock=iIndex_II(0,iParticle),&
            V_D = VLocal_D)
       ! update particle location
       State_VI(x_:z_, iParticle) = State_VI(x_:z_, iParticle) + &
            Dt * VLocal_D
    end do
    !\
    ! Message pass: some particles may have moved to different procs
    !/
    call message_pass_particles(KindReg_)
    call remove_undefined_particles(KindReg_)

  contains

    subroutine get_v(Xyz_D, iBlock, V_D)
      ! get local plasma velocity
      real,   intent(in) :: Xyz_D(MaxDim)
      integer,intent(in) :: iBlock
      real,   intent(out):: V_D(MaxDim)


      ! variables for AMR interpolation
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim) = 1
      !------------------------------------
      ! reset the interpoalted values
      V_D = 0!; Rho = 0; M_D = 0
      ! get the velocity
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
      ! interpolate the local density and momentum
      do iCell = 1, nCell
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         if(State_VGB(Rho_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell) <= 0)&
              call CON_stop(NameSub//": zero or negative plasma density")
         ! convert momentum to velocity
         V_D = V_D + &
              State_VGB(RhoUx_:RhoUz_,i_D(1),i_D(2),i_D(3),iBlock) / &
              State_VGB(Rho_,         i_D(1),i_D(2),i_D(3),iBlock) * &
              Weight_I(iCell)
      end do
      ! check of there is plasma at the location
!      if(Rho==0)&
!           call CON_stop(NameSub//&
!           ': trying to advect particle line at region with no plasma')
!      V_D = M_D / Rho
    end subroutine get_v
  end subroutine advect_particle_line

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

  subroutine get_particle_data(nData, NameVar, DataOut_VI, nParticle)
    ! the subroutine gets variables specified in the string StringVar
    ! and writes them into DataOut_VI
    integer,            intent(in) :: nData
    character(len=*),   intent(in) :: NameVar
    real, allocatable,  intent(out):: DataOut_VI(:,:)
    integer,            intent(out):: nParticle

    ! mask for returning variables
    logical:: DoReturnVar_V(nVarParticle)
    logical:: DoReturnIndex_I(0:nIndexParticle)
    ! number of variables/indices found in the request string
    integer:: nVarOut, nIndexOut
    ! loop variables
    integer:: iParticle, iData
    ! permutation of between data order of the request and order in which
    ! data are checked to be present in the request
    integer:: iDataPermutation_I(nData)

    ! data that can be return
    integer, parameter:: nVarAvail = 3, nIndexAvail = 2
    character(len=2), parameter:: NameVarAvail_V(nVarAvail) = &
         (/'xx', 'yy', 'zz'/)
    character(len=2), parameter:: NameIndexAvail_I(nIndexAvail) = &
         (/'fl', 'id'/)

    character(len=*), parameter:: NameSub = 'get_particle_data'
    !----------------------------------------------------------------------
    !\
    ! first, determine which data are requested
    !/
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.
    ! determine which variables are requested
    nVarOut = 0
    do iData = 1, nVarAvail
       DoReturnVar_V(iData)   = is_var(NameVarAvail_V(iData))
    end do
    nVarOut = count(DoReturnVar_V)
    ! determine which indices are requested
    nIndexOut = 0
    do iData = 1, nIndexAvail
       DoReturnIndex_I(iData) = is_var(NameIndexAvail_I(iData))
    end do
    nIndexOut = count(DoReturnIndex_I)
    ! check that number of data found is the same as requested
    if(nData /= nVarOut + nIndexOut)&
         call CON_stop(NameSub//': incorrect number of data is requested')

    !\
    ! return data
    !/
    if(allocated(DataOut_VI)) deallocate(DataOut_VI)
    nParticle = Particle_I(KindReg_)%nParticle
    allocate( DataOut_VI(nData, nParticle) )
    do iParticle = 1, nParticle
       ! store variables
       DataOut_VI(1:nVarOut,      iParticle) = PACK(&
            Particle_I(KindReg_)%State_VI( :,iParticle),&
            MASK = DoReturnVar_V)
       ! store indexes
       DataOut_VI(nVarOut+1:nData,iParticle) = PACK(&
            Particle_I(KindReg_)%iIndex_II(:,iParticle),&
            MASK = DoReturnIndex_I)
    end do
    ! permute data order so it corresponds to the order of the request
    DataOut_VI(:, :) = DataOut_VI(iDataPermutation_I, :)

  contains

    function is_var(Name) result(IsPresent)
      ! finds if datum Name is present
      character(len=*), intent(in):: Name
      logical:: IsPresent
      integer:: i
      !--------------------------------------------------------------------
      i = index(NameVar, Name)
      if(i==0)then
         IsPresent = .false.
         RETURN
      end if
      IsPresent = .true.
      ! save datum's order in the original request
      iDataPermutation_I(i/3 + 1) = nVarOut + iData
    end function is_var

  end subroutine get_particle_data

end module ModParticleFieldLine
