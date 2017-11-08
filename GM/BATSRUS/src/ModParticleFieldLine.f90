!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModParticleFieldLine

  ! the module contains subroutine for extracting magnetic field lines 
  ! for passing to other codes;
  ! field line intergration is performed with use of BATL library including
  ! - continuous AMR interpolation
  ! - particle methods

  use BATL_lib, ONLY: &
       iProc,&
       MaxDim, nDim, &
       check_interpolate_amr_gc, &
       Particle_I, Xyz_DGB, CellVolume_GB, &
       MaxBlock, nI, nJ, nK, nBlock, &
       allocate_particles, Unused_B, &
       message_pass_particles, remove_undefined_particles, &
       mark_undefined, check_particle_location, get_particles
       
  use ModBatlInterface, ONLY: interpolate_grid_amr_gc
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, B_, Bx_, Bz_
  use ModMain, ONLY: NameThisComp
  use ModPhysics, ONLY: rBody
  use ModCellGradient, ONLY: calc_gradient_ghost, &
       GradCosBR_DGB => GradVar_DGB

  implicit none
  SAVE

  private !except

  public:: n_particle_reg_max
  public:: n_particle_reg
  public:: read_particle_line_param
  public:: init_particle_line
  public:: set_soft_boundary
  public:: extract_particle_line
  public:: add_to_particle_line
  public:: advect_particle_line
  public:: get_particle_data
  public:: write_plot_particle

  interface get_particle_data
     module procedure get_particle_data_all
     module procedure get_particle_data_i_particle
  end interface

  ! Local variables ----------------------

  ! kinds of particles used to generate a magnetic field line
  integer, parameter:: &
       KindEnd_ = 1, &
       KindReg_ = 2

  ! variable in the state vector of a particle
  integer, parameter:: &
       ! coordinates of a particle
       x_    = 1, y_    = 2, z_    = 3, & 
       ! storage for particle's position at the beginning of the step
       XOld_ = 4, YOld_ = 5, ZOld_ = 6, & 
       ! stepsize during line extraction
       Ds_   = 7, &
       !grad CosBR at the beginning of extraction iteration
       GradX_= 8, GradZ_ = 10, &
       ! value of CosBR at the beginning of extraction iteration
       CosBR_    = 11, &
       ! previous direction
       DirX_ = 12, DirZ_ = 14

  ! indices of a particle
  integer, parameter:: &
       ! field line this particle lays on
       fl_        = 1, &
       ! index of the particle along this field line
       id_        = 2, &
       ! alignment of particle numbering with direction of B field:
       ! ( -1 -> reversed, +1 -> aligned); important for radial ordering 
       ! when magnetic field general direction may be inward
       Alignment_ = 3, &
       ! indicator for message passing in trace_particle_line
       Pass_ = 4

  integer, parameter::    &
       Normal_ = 0       ,&
       HalfStep_ = 1     ,&
       DoneFromScratch_ = 2
  

  ! data that can be requested and returned to outside this module
  ! ORDER MUST CORRESPOND TO INDICES ABOVE
  integer, parameter:: nVarAvail = 3, nIndexAvail = 2
  character(len=2), parameter:: NameVarAvail_V(nVarAvail) = &
       (/'xx', 'yy', 'zz'/)
  character(len=2), parameter:: NameIndexAvail_I(0:nIndexAvail) = &
       (/'bk','fl', 'id'/)


  ! spatial step limits at the field line extraction
  real   :: SpaceStepMin = 0.0
  real   :: SpaceStepMax = HUGE(SpaceStepMin)

  ! ordering mode of particles along field lines
  integer:: iOrderMode = -1
  integer, parameter:: Field_ = 0, Radius_  = 1

  ! number of variables in the state vector
  integer, parameter:: nVarParticleReg = 6
  integer, parameter:: nVarParticleEnd = 14

  ! number of indices
  integer, parameter:: nIndexParticleReg = 2
  integer, parameter:: nIndexParticleEnd = 4

  ! maximum allowed number of field lines
  integer, parameter :: nFieldLineMax = 1000

  ! number of active field lines
  integer:: nFieldLine = 0

  ! soft radial boundary that may be set only once during the run;
  ! "soft" means that particles are allowed beyond this boundary BUT:
  ! - during extraction first particle that crosses it is kept but extraction 
  !   of this line is stopped
  real:: RBoundarySoft = -1.0

  ! may need to apply corrections during line tracing
  logical:: UseCorrection=.false.
  real:: RCorrectionMin = 0.0, RCorrectionMax = Huge(1.0)

  ! initialization related info
  integer:: nLineInit
  real, allocatable:: XyzLineInit_DI(:,:)
  !\
  ! Shared variavles
  !/
  !\
  ! direction of tracing: -1 -> backward, +1 -> forward
  integer:: iDirTrace

  ! parameters of end particles
  real,    pointer:: StateEnd_VI(:,:)
  integer, pointer:: iIndexEnd_II(:,:)
contains
  !==========================================================================
  integer function n_particle_reg_max()
    n_particle_reg_max = Particle_I(KindReg_)%nParticleMax
  end function n_particle_reg_max
  !==========================================================================
  integer function n_particle_reg()
    n_particle_reg = Particle_I(KindReg_)%nParticle
  end function n_particle_reg
  !==========================================================================
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
             call CON_stop(NameThisComp//':'//NameSub //&
                  ": unknown ordering mode")
          end if
          !--------------------------------------------------------------
          ! field lines may need to be corrected during tracing
          call read_var('UseCorrection', UseCorrection)
          if(UseCorrection)then
             ! use may specify the boundaries where the correction is applied
             call read_var('RCorrectionMin', RCorrectionMin)
             if(RCorrectionMin < 0.0) RCorrectionMin = 0.0
             call read_var('RCorrectionMax', RCorrectionMax)
             if(RCorrectionMax < 0.0) RCorrectionMax = Huge(1.0)
          end if
       end if
    end select
  end subroutine read_particle_line_param

  !==========================================================================

  subroutine init_particle_line
    ! allocate containers for particles
    !------------------------------------------------------------------------
    Particle_I(KindReg_)%nParticleMax = 10000 * nFieldLineMax
    Particle_I(KindEnd_)%nParticleMax = nFieldLineMax
    Particle_I(KindReg_)%nVar   = nVarParticleReg
    Particle_I(KindEnd_)%nVar   = nVarParticleEnd
    Particle_I(KindReg_)%nIndex = nIndexParticleReg
    Particle_I(KindEnd_)%nIndex = nIndexParticleEnd
    call allocate_particles
    ! set pointers to parameters of end particles
    StateEnd_VI  => Particle_I(KindEnd_)%State_VI
    iIndexEnd_II => Particle_I(KindEnd_)%iIndex_II
    if(nLineInit > 0)&
         ! extract initial field lines
         call extract_particle_line(nLineInit, XyzLineInit_DI)
  end subroutine init_particle_line
  !==========================================================================
  subroutine set_soft_boundary(RBoundarySoftIn)
    ! set additional boundary imposed by user
    real, intent(in):: RBoundarySoftIn
    character(len=*), parameter:: NameSub = 'set_soft_boundary'
    !---------------------------------------------
    ! soft boundary can be set only once!
    if(RBoundarySoft < 0.0)then
       RBoundarySoft = RBoundarySoftIn
    else
       call CON_stop(NameThisComp//':'//NameSub//&
            ': soft boundary may be set only once and may not be changed!')
    end if
  end subroutine set_soft_boundary
  !=========================================================================
  subroutine extract_particle_line(nFieldLineIn, XyzStart_DI, iTraceModeIn, &
       iIndexStart_II,&
       UseInputInGenCoord)
    use BATL_lib, ONLY: coord_to_xyz
    use BATL_particles, ONLY: trace_particles
    ! extract nFieldLineIn magnetic field lines starting at XyzStart_DI;
    ! the whole field lines are extracted, i.e. they are traced forward
    ! and backward up until it reaches boundaries of the domain;
    ! requested coordinates may be different for different processor,
    ! if a certain line can't be started on a given processor, it is
    ! ignored, thus duplicates are avoided
    ! NOTE: different sets of lines may be request on different processors!
    !-----------------------------------------------------------------------
    integer,         intent(in)::nFieldLineIn
    real,            intent(in)::XyzStart_DI(MaxDim, nFieldLineIn)
    ! mode of tracing (forward, backward or both ways)
    integer,optional,intent(in)::iTraceModeIn
    ! initial particle indices for starting particles
    integer,optional,intent(in) :: &
         iIndexStart_II(nIndexParticleReg,nFieldLineIn)

    ! An input can be in generalized coordinates
    logical, optional,intent(in) :: UseInputInGenCoord 
    integer :: nLineThisProc! number of new field lines initialized locally
    integer :: nLineAllProc ! number of new field lines initialized on all PEs
    integer :: nParticleOld ! number of already existing regular particles

    ! mode of tracing (see description below)
    integer:: iTraceMode

    ! Cosine between direction of the field and radial direction: 
    ! for correcting the line's direction to prevent it from closing
    real, allocatable:: CosBR_CB(:,:,:,:)
    character(len=*), parameter:: NameSub='extract_particle_line'
    integer, pointer :: iIndexReg_II(:,:)
    !-----------------------------
    iIndexReg_II=> Particle_I(KindReg_)%iIndex_II
    ! initialize field lines
    call  get_particles(&
         iKindParticle      = KindEnd_,          &
         StateIn_VI         = XyzStart_DI,       &
         iLastIdIn          = nFieldLine,        &
         iIndexIn_II        = iIndexStart_II,    & 
         UseInputInGenCoord = UseInputInGenCoord,&
         DoReplace          = .true.,            &    
         nParticlePE        = nLineThisProc)

    nParticleOld  = Particle_I(KindReg_)%nParticle

    ! how many lines have been started on all processors
    call count_new_lines()
    nFieldLine    = nFieldLine + nLineAllProc
    if(nFieldLine > nFieldLineMax)&
         call CON_stop(NameThisComp//':'//NameSub//&
         ': Limit for number of particle field lines exceeded')

    call copy_end_to_regular
    ! allocate containers for cosine of angle between B and radial direction
    allocate(CosBR_CB(1:nI, 1:nJ, 1:nK, MaxBlock))
    call compute_cosbr
    call calc_gradient_ghost(CosBR_CB)
    ! free space
    deallocate(CosBR_CB)


    ! Trace field lines

    ! check if trace mode is specified
    if(present(iTraceModeIn))then
       if(abs(iTraceModeIn) > 1)&
            call CON_stop(&
            NameThisComp//':'//NameSub//': incorrect tracing mode')
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

       ! fix alignment of particle indexing with B field direction
       call get_alignment()

       call trace_particles(KindEnd_, particle_line)

       ! if just finished backward tracing and need to trace in both dirs
       ! => return to initial particles
       if(iDirTrace < 0 .and. iTraceMode == 0)then
          ! copy initial points to KindEnd_:
          ! the initial particles are currently right after the particles,
          ! that were in the list before the start of this subroutine,
          ! i.e. occupy positions from (nParticleOld+1)
          StateEnd_VI(x_:z_, 1:nLineThisProc) = Particle_I(KindReg_&
               )%State_VI(x_:z_, nParticleOld+1:nParticleOld+nLineThisProc)
          iIndexEnd_II(0:id_,1:nLineThisProc) = iIndexReg_II(&
               0:id_,nParticleOld+1:nParticleOld+nLineThisProc)
          Particle_I(KindEnd_)%nParticle = nLineThisProc
       end if
    end do
    ! Offset in id_
    call offset_id(nFieldLine - nLineAllProc + 1, nFieldLine)
  contains
    !========================================================================
    subroutine compute_cosbr()
      use ModMain, ONLY: UseB0
      use ModB0, ONLY: B0_DGB
      use ModGeometry, ONLY: R_BLK

      integer :: iBlock, i, j, k ! loop variables
      real    :: XyzCell_D(nDim), BCell_D(nDim)
      !------------------------------------------------------------------------
      ! precompute CosBR on the grid for all active blocks
      do iBlock = 1, nBlock
         if(Unused_B(iBlock))CYCLE
         do k = 1, nK; do j = 1, nJ; do i = 1, nI
            XyzCell_D = Xyz_DGB(1:nDim,i,j,k,iBlock)
            Bcell_D = State_VGB(Bx_:B_+nDim,i,j,k,iBlock)
            if(UseB0)BCell_D = Bcell_D + B0_DGB(1:nDim,i,j,k,iBlock)
            if(any(BCell_D /= 0.0) .and. R_BLK(i,j,k,iBlock) > 0.0)then
               CosBR_CB(i,j,k,iBlock) = sum(Bcell_D*XyzCell_D) / &
                    (sqrt(sum(BCell_D**2))*R_BLK(i,j,k,iBlock))
            else
               CosBR_CB(i,j,k,iBlock) = 0.0
            end if
         end do; end do; end do
      end do
    end subroutine compute_cosbr
    !=============
    subroutine  offset_id(iLineStart,iLineEnd)
      use ModMpi
      use BATL_mpi, ONLY: iComm, nProc
      integer, intent(in) :: iLineStart, iLineEnd
      integer:: iParticle, nParticle, iLine, iError
      integer, dimension(1:nFieldLineMax):: iOffsetLocal_I, iOffset_I
      !------------
      iOffsetLocal_I = 0; iOffset_I = 0  

      ! set a pointer to parameters of regular particles
      nParticle = Particle_I(KindReg_)%nParticle
      do iParticle = 1, nParticle
         iLine = iIndexReg_II(fl_,iParticle)

         ! Reset offset, if id_ is less than one
         iOffsetLocal_I(iLine) = max(iOffsetLocal_I(iLine),&
              1 - iIndexReg_II(id_,iParticle))
      end do
      ! find maximim offset from all processors
      if(nProc > 1) then
         call MPI_Allreduce(iOffsetLocal_I(iLineStart), &
              iOffset_I(iLineStart), 1 + iLineEnd - iLineStart, &
              MPI_INTEGER, MPI_MAX, iComm, iError)
      else
         iOffset_I(iLineStart:iLineEnd) = &
              iOffsetLocal_I(iLineStart:iLineEnd)
      end if
      do iParticle = 1, nParticle
         iLine = iIndexReg_II(fl_,iParticle)

         ! Apply offset
         iIndexReg_II(id_,iParticle) = iIndexReg_II(id_,iParticle) + &
              iOffset_I(iLine)
      end do
    end subroutine offset_id
    !========================================================================
    subroutine count_new_lines()
      ! gather information from all processors on how many new field lines
      ! have been started
      use ModMpi
      use BATL_mpi, ONLY: iComm, nProc
      integer:: iError
      !-----------------------------------------------------------------------
      if(nProc>1)then
         call MPI_Allreduce(nLineThisProc, nLineAllProc, 1, &
              MPI_INTEGER, MPI_SUM, iComm, iError)
      else
         nLineAllProc = nLineThisProc
      end if
    end subroutine count_new_lines
    !========================================================================
    subroutine get_alignment()
      ! determine alignment of particle indexing with direction 
      ! of the magnetic field
      integer:: iParticle
      real   :: DirB_D(MaxDim)
      !---------------------------------------------------------------------
      iIndexEnd_II(Pass_, 1:Particle_I(KindEnd_)%nParticle) = DoneFromScratch_
      if(iOrderMode == Field_)then
         iIndexEnd_II(Alignment_, 1:Particle_I(KindEnd_)%nParticle) = 1
         RETURN
      end if

      do iParticle = 1, Particle_I(KindEnd_)%nParticle
         call get_b_dir(&
              Xyz_D = StateEnd_VI(x_:z_, iParticle),&
              iBlock=iIndexEnd_II(0,iParticle),&
              Dir_D = DirB_D)
         iIndexEnd_II(Alignment_, iParticle) = &
              nint( SIGN(1.0, sum(DirB_D*StateEnd_VI(x_:z_, iParticle)) ))
      end do
    end subroutine get_alignment
    !=================
  end subroutine extract_particle_line
  !=========================================================================
  subroutine particle_line(iParticle, IsEndOfSegment)
    ! the subroutine is the tracing loop;
    ! tracing starts at KindEnd_ particles and proceeds in a given direction
    !-----------------------------------------------------------------------
    integer, intent(in) :: iParticle     ! loop variable
    logical, intent(out):: IsEndOfSegment
    !---------------
    ! direction of the magnetic field
    real, dimension(MaxDim) :: DirBCurr_D, DirBNext_D
    ! direction of the tangent to the line
    real, dimension(MaxDim) :: DirCurr_D, DirNext_D

    ! radii and vector of radial direction
    real :: ROld, RNew, DirR_D(MaxDim)
    ! cosine of angle between direction of field and radial direction
    real:: CosBR

    ! whether particle's location isn't covered by true cells
    logical:: IsBody

    ! whether to schedule a particle for message pass
    logical:: DoMove, IsGone
    
    character(len=*),parameter:: NameSub = "particle_line"
    !-------------------------------------------------------------
    select case(iIndexEnd_II(Pass_,iParticle))
    case(HalfStep_)
       ! change Pass_ field for particle passed after first stage
       iIndexEnd_II(Pass_,iParticle) = Normal_
       ! do not call the first stage, proceed to the second one
    case(DoneFromScratch_)
       ! Initialize the radial direction that corresponds 
       ! to the previous step
       StateEnd_VI(DirX_:DirZ_,  iParticle) = iDirTrace * &
            StateEnd_VI(x_:z_,  iParticle) / &
            sqrt(sum(StateEnd_VI(x_:z_, iParticle)**2))
       iIndexEnd_II(Pass_, iParticle) = Normal_
       call stage1(iParticle, IsEndOfSegment)
       if(IsEndOfSegment) RETURN
    case(Normal_)
       ! increase particle index & copy to regular
       iIndexEnd_II(id_,iParticle) = &
            iIndexEnd_II(id_,iParticle) + iDirTrace
       call copy_end_to_regular(iParticle)
       call stage1(iParticle, IsEndOfSegment)
       if(IsEndOfSegment) RETURN
    case default
       call stop_mpi(NameSub//': Unknown stage ID')
    end select
    call stage2(iParticle, IsEndOfSegment)
 contains
    !========================
    subroutine stage1(iParticle, IsEndOfSegment)
      integer, intent(in) :: iParticle
      Logical, intent(out):: IsEndOfSegment
      !-------------------
      IsEndOfSegment = .false.
      if(RBoundarySoft > 0.0)then
         if(sum(Particle_I(KindEnd_)%State_VI(1:nDim,iParticle)**2) &
              > rBoundarySoft**2)then
            call mark_undefined(KindEnd_, iParticle)
            IsEndOfSegment = .true.
            RETURN
         end if
      end if
      !\
      ! First stage of RK2 method
      !/
      ! copy last known coordinates to old coords
      StateEnd_VI(XOld_:ZOld_,iParticle) = &
           StateEnd_VI(x_:z_, iParticle)
      !\
      ! get the direction of the magnetic field at original location 
      ! and gradient of cosine of angle between field and
      ! radial direction
      !/
      call get_b_dir(&
           Xyz_D   = StateEnd_VI(x_:z_, iParticle),&
           iBlock  =iIndexEnd_II(0,iParticle),&
           Dir_D   = DirBCurr_D,&
           Grad_D  = StateEnd_VI(GradX_:GradZ_,iParticle),&
           StepSize= StateEnd_VI(Ds_,iParticle),&
           IsBody  = IsBody)
      
      ! particle's location isn't covered by mostly true cells, 
      ! i.e. not well outside central body
      if(IsBody)then
         call mark_undefined(KindEnd_,iParticle)
         IsEndOfSegment = .true.
         RETURN
      end if
      
      !\
      ! Limit the interpolated time step
      StateEnd_VI(Ds_, iParticle) = &
           MIN(SpaceStepMax, MAX(SpaceStepMin,&
           StateEnd_VI(Ds_, iParticle)))
      
      ! save direction for the second stage
      DirCurr_D = DirBCurr_D * &
           iDirTrace * iIndexEnd_II(Alignment_, iParticle)
      !\
      ! cosine of angle between field and radial direction
      !
      CosBR = iDirTrace * iIndexEnd_II(Alignment_, iParticle) * &
           sum(StateEnd_VI(x_:z_, iParticle)*DirCurr_D) / &
           sqrt(sum(StateEnd_VI(x_:z_, iParticle)**2))
      StateEnd_VI(CosBR_, iParticle) = CosBR
      
      ! correct the direction to prevent the line from closing
      if(UseCorrection) &
           call correct(DirCurr_D)

      ! get middle location
      StateEnd_VI(x_:z_, iParticle) = &
           StateEnd_VI(x_:z_, iParticle) + &
           0.50*StateEnd_VI(Ds_, iParticle) *  DirCurr_D
      
      ! update location and schedule for message pass
      call check_particle_location(  &
           iKindParticle = KindEnd_ ,&
           iParticle     = iParticle,&
           DoMove        = DoMove   ,&
           IsGone        = IsGone    )
      if(IsGone)then
         IsEndOfSegment = .true.
      elseif(DoMove) then
         iIndexEnd_II(Pass_, iParticle) = HalfStep_
         IsEndOfSegment = .true.
      end if
    end subroutine stage1
    !====================
    subroutine stage2(iParticle, IsEndOfSegment)
      integer, intent(in) :: iParticle
      Logical, intent(out):: IsEndOfSegment
      !-------------------
      IsEndOfSegment = .false.  
      ! get the direction of the magnetic field in the middle
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBNext_D,&
           IsBody  = IsBody)
      
      ! particle's location isn't covered by mostly true cells, 
      ! i.e. not well outside central body
      if(IsBody)then
         call mark_undefined(KindEnd_,iParticle)
         IsEndOfSegment = .true.
         RETURN
      end if
      
      ! direction at the 1st stage 
      DirNext_D = DirBNext_D *&
           iDirTrace * iIndexEnd_II(Alignment_, iParticle)
      CosBR = StateEnd_VI(CosBR_, iParticle)  
      
      ! correct the direction to prevent the line from closing
      if(UseCorrection) &
           call correct(DirNext_D)
      
      ! get final location
      StateEnd_VI(x_:z_,iParticle) = &
           StateEnd_VI(XOld_:ZOld_,iParticle) + &
           StateEnd_VI(Ds_, iParticle)*DirNext_D
      !\
      ! Achieve that the change in the radial distance
      ! equals StateEnd_VI(Ds_, iParticle)*sum(DirNext_D*DirR_D)
      !/
      ROld = sqrt(sum(StateEnd_VI(XOld_:ZOld_,iParticle)**2))
      DirR_D = StateEnd_VI(XOld_:ZOld_,iParticle)/ROld
      
      RNew = sqrt(sum(StateEnd_VI(x_:z_,iParticle)**2))
      StateEnd_VI(x_:z_,iParticle) = StateEnd_VI(x_:z_,iParticle)/RNew*&
           ( ROld + StateEnd_VI(Ds_, iParticle)*sum(DirNext_D*DirR_D) )
      ! update the direction of the previous step
      StateEnd_VI(DirX_:DirZ_,  iParticle) = DirNext_D
      
      ! update location and schedule for message pass
      call check_particle_location(  &
           iKindParticle = KindEnd_ ,&
           iParticle     = iParticle,&
           DoMove        = DoMove   ,&
           IsGone        = IsGone    )
      IsEndOfSegment = IsGone.or.DoMove
    end subroutine stage2
    !========================
    subroutine correct(Dir_D)
      ! corrects the direction to prevent the line from closing:
      ! if CosBR is far enough from the value of the alginment (+/-1),
      ! Dir_D is changed to a vector with max radial component 
      ! \perp \nabla CosBR to avoid reaching region with opposite sign of CosBR
      !-----------------------------------
      ! current direction
      real, intent(inout):: Dir_D(MaxDim)
      ! aux vectors
      real, dimension(MaxDim):: C_D, R_D, A_D
      ! deviation from radial direction
      real:: DCosBR
      ! threshold for deviation from radial direction to apply correction
      real:: DCosBRMax = 0.15
      ! scalars to define parameters of the correction
      real:: Misc, Rad
      real:: CRad, CTan
      real:: KappaTh, Kappa
      !-----------------------------------
      ! the deviation of the field line from radial direction
      DCosBR = abs(CosBR-iIndexEnd_II(Alignment_,iParticle))
      ! radius
      Rad = sqrt(sum(StateEnd_VI(x_:z_,iParticle)**2))

      ! if the line deviates to0 much -> correct its direction
      ! to go along surface B_R = 0 
      if(DCosBR > DCosBRMax .and. & 
           Rad > RCorrectionMin .and. Rad < RCorrectionMax)then
         ! unit vector along the direction of gradient
         ! with alignment accounted for
         C_D = iIndexEnd_II(Alignment_, iParticle) * &
              StateEnd_VI(GradX_:GradZ_, iParticle)/ &
              sqrt(sum(StateEnd_VI(GradX_:GradZ_, iParticle)**2))
         ! unit vector in radial direction
         R_D = StateEnd_VI(x_:z_, iParticle) / Rad
              
         ! vector \perp C_D with max possible radial component
         CRad = sum(C_D*R_D)
         CTan = sqrt(1.0 - CRad**2)
         A_D = iDirTrace * (CTan * R_D  - (C_D - CRad*R_D) * CRad / CTan)

         ! change the direction of the line
         Dir_D = A_D
         RETURN
      end if
      
      ! the line doesn't deviate much, but it may just have exited the region
      ! where the correction above has been applied; 
      ! to prevent sharp changes in the direction -> limit line's curvature

      ! cut-off curvature is based on the local grid size
      KappaTh = 0.1 / StateEnd_VI(Ds_, iParticle)
      ! curvature based on uncorrected direction
      Kappa  = &
           sqrt(2*(1-sum(Dir_D*StateEnd_VI(DirX_:DirZ_, iParticle)))) / &
           StateEnd_VI(Ds_, iParticle)
      ! if curvature isn't large -> don't apply the correction
      if(Kappa < KappaTh) RETURN

      ! remove scale the component of the current direction parallel
      ! to the one on the previous step;
      ! components \perp to it are scaled accordingly to keep ||Dir_D|| = 1
      Dir_D = Dir_D - &
           sum(Dir_D*StateEnd_VI(DirX_:DirZ_,iParticle)) * &
           StateEnd_VI(DirX_:DirZ_,iParticle)
      Dir_D = Dir_D / sqrt(sum(Dir_D**2))
      Misc  = 1 - 0.5 * (KappaTh*StateEnd_VI(Ds_, iParticle))**2
      Dir_D = Misc * StateEnd_VI(DirX_:DirZ_,iParticle) +&
           sqrt(1-Misc**2) * Dir_D
    end subroutine correct
    !====================
  end subroutine particle_line
  !========================================================================
  subroutine copy_end_to_regular(iParticleIn)
    integer, optional, intent(in) :: iParticleIn
    ! copies indicated variables of known end particles to regular particles
    integer, parameter:: iVarCopy_I(3)   = (/x_, y_, z_/)
    integer, parameter:: iIndexCopy_I(3) = (/0, fl_, id_/)
    integer:: iParticle, iParticleStart, iParticleEnd
    !----------------------------------------------------------------------
    if(present(iParticleIn))then
       iParticleStart = iParticleIn
       iParticleEnd   = iParticleIn
    else
       iParticleStart = 1
       iParticleEnd   = Particle_I(KindEnd_)%nParticle
    end if
    do iParticle = iParticleStart, iParticleEnd 
       Particle_I(KindReg_)%nParticle = Particle_I(KindReg_)%nParticle+1
       Particle_I(KindReg_)%State_VI(iVarCopy_I, &
            Particle_I(KindReg_)%nParticle) =&
            Particle_I(KindEnd_)%State_VI(iVarCopy_I, iParticle)
       Particle_I(KindReg_)%iIndex_II(iIndexCopy_I,&
            Particle_I(KindReg_)%nParticle) =&
            Particle_I(KindEnd_)%iIndex_II(iIndexCopy_I, iParticle)
    end do
  end subroutine copy_end_to_regular
  !================================
  subroutine get_b_dir(Xyz_D, iBlock, Dir_D, Grad_D, StepSize, IsBody)
    use ModMain, ONLY: UseB0
    use ModB0, ONLY: get_b0

    ! returns the direction of magnetic field 
    ! as well as the block used for interpolation

    real,    intent(in) :: Xyz_D(MaxDim)
    integer, intent(in) :: iBlock
    real,    intent(out):: Dir_D(MaxDim)
    real,    optional, intent(out):: Grad_D(MaxDim)
    real,    optional, intent(out):: StepSize
    logical, optional, intent(out):: IsBody

    ! magnetic field
    real   :: B_D(MaxDim) = 0.0
    ! interpolation data: number of cells, cell indices, weights
    integer:: nCell, iCell_II(0:nDim, 2**nDim)
    real   :: Weight_I(2**nDim)
    integer:: iCell ! loop variable
    integer:: i_D(MaxDim)
    character(len=200):: StringError
    character(len=*), parameter:: NameSub = "get_b_dir"
    !----------------------------------------------------------------------
    Dir_D = 0; B_D = 0
    ! get potential part of the magnetic field at the current location
    if(UseB0)call get_b0(Xyz_D, B_D)
    ! get the remaining part of the magnetic field
    call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I,&
         IsBody)

    ! interpolate magnetic field value
    do iCell = 1, nCell
       i_D = 1
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
    if(present(Grad_D))then
       Grad_D = 0
       ! interpolate grad(cos b.r)
       do iCell = 1, nCell
          i_D = 1
          i_D(1:nDim) = iCell_II(1:nDim, iCell)
          Grad_D(1:nDim) = Grad_D(1:nDim) + &
               GradCosBR_DGB(:,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
       end do
    end if
    if(present(StepSize))then
       StepSize = 0.0
       ! interpolate cell sizes. For non-cartesian grids 
       !the metic tensor is used
       do iCell = 1, nCell
          i_D = 1
          i_D(1:nDim) = iCell_II(1:nDim, iCell)
          StepSize = StepSize + &
               CellVolume_GB(i_D(1),i_D(2),i_D(3),iBlock) * Weight_I(iCell) 
       end do
       ! take a fraction of cubic root of inteprolated cell volume as step size
       StepSize = 0.1 * StepSize**(1.0/3.0)
    end if
  end subroutine get_b_dir
  !========================================================================
  subroutine add_to_particle_line(nParticleIn, XyzIn_DI, iIndexIn_II,&
       UseInputInGenCoord, DoReplace)

    use BATL_lib, ONLY: coord_to_xyz

    ! add particles with specified coordinates to the already existing lines
    integer, intent(in):: nParticleIn
    real,    intent(in):: XyzIn_DI(MaxDim, nParticleIn)
    integer, intent(in):: iIndexIn_II(nIndexParticleReg, nParticleIn)

    ! An input can be in generalized coordinates
    logical, optional,intent(in) :: UseInputInGenCoord 

    ! Whether to replace ALL old particles with the input
    logical, optional, intent(in):: DoReplace
    !\
    call get_particles(&
         iKindParticle      = KindReg_          ,&
         StateIn_VI         = XyzIn_DI          ,&
         iIndexIn_II        = iIndexIn_II       ,&
         UseInputInGenCoord = UseInputInGenCoord,&
         DoReplace          = DoReplace          )
  end subroutine add_to_particle_line

  !========================================================================

  subroutine advect_particle_line
    ! advect particles with the local plasma velocity
    use ModMain, ONLY: time_accurate, Dt
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

    ! set pointers to parameters of particles
    State_VI  => Particle_I(KindReg_)%State_VI
    iIndex_II => Particle_I(KindReg_)%iIndex_II

    ! go over the list of particles and advect them
    do iParticle = 1, Particle_I(KindReg_)%nParticle

       ! simple advection (Euler integration)
       ! get the local velocity
       call get_v(&
            Xyz_D = State_VI(x_:z_, iParticle),&
            iBlock=iIndex_II(0,iParticle),&
            V_D = VLocal_D)

       ! update particle location
       State_VI(x_:z_, iParticle) = State_VI(x_:z_, iParticle) + &
            Dt * VLocal_D
    end do

    ! Message pass: some particles may have moved to different procs
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
      integer:: i_D(MaxDim)
      !------------------------------------
      ! reset the interpoalted values
      V_D = 0 
      ! get the velocity
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
      ! interpolate the local density and momentum
      do iCell = 1, nCell
         i_D = 1
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         if(State_VGB(Rho_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell) <= 0)&
              call CON_stop(NameSub//": zero or negative plasma density")
         ! convert momentum to velocity
         V_D = V_D + &
              State_VGB(RhoUx_:RhoUz_,i_D(1),i_D(2),i_D(3),iBlock) / &
              State_VGB(Rho_,         i_D(1),i_D(2),i_D(3),iBlock) * &
              Weight_I(iCell)
      end do
    end subroutine get_v
  end subroutine advect_particle_line
  !=====================
  
  subroutine process_request(StringRequest, &
       nVarReturn, DoReturnVar_V, nIndexReturn, DoReturnIndex_I, iOrder_I)
    ! finds if variable Name is requested via StringRequest,
    ! this function is called in get_particle_data to identify the variable
    ! being requested from outside this module;
    ! the function also returns the position iPos of Name in StringRequest
    character(len=*),intent(in)  :: StringRequest
    integer,         intent(out) :: nVarReturn
    logical,         intent(out) :: DoReturnVar_V(nVarAvail)
    integer,         intent(out) :: nIndexReturn
    logical,         intent(out) :: DoReturnIndex_I(0:nIndexAvail)
    integer,         intent(out) :: iOrder_I(nVarAvail+nIndexAvail+1)

    integer:: iPos, iVar
    !--------------------------------------------------------------------
    do iVar = 1, nVarAvail
       iPos = index(StringRequest, NameVarAvail_V(iVar))
       DoReturnVar_V(iVar) = iPos > 0
       iOrder_I(iVar) = iPos
    end do
    nVarReturn = count(DoReturnVar_V)
    do iVar = 0, nIndexAvail
       iPos = index(StringRequest, NameIndexAvail_I(iVar))
       DoReturnIndex_I(iVar) = iPos > 0
       iOrder_I(nVarAvail + iVar + 1) = iPos
    end do
    nIndexReturn = count(DoReturnIndex_I)

    do iVar = 1, nVarReturn+nIndexReturn
       iPos = minloc(iOrder_I, 1, MASK = iOrder_I >= iVar)
       iOrder_I(iPos) = iVar
    end do
  end subroutine process_request
  
  !========================================================================

  subroutine get_particle_data_all(nData, NameVar, DataOut_VI)
    ! the subroutine gets variables specified in the string StringVar
    ! and writes them into DataOut_VI;
    ! data is for all particles otherwise
    ! NOTE: function n_particle_reg can be used to find size beforehand
    integer, intent(in) :: nData
    real,    intent(out):: DataOut_VI(nData,Particle_I(KindReg_)%nParticle)
    character(len=*),intent(in) :: NameVar

    ! mask for returning variables
    logical:: DoReturnVar_V(nVarParticleReg)
    logical:: DoReturnIndex_I(0:nIndexParticleReg)
    ! number of variables/indices found in the request string
    integer:: nVarOut, nIndexOut
    ! loop variables
    integer:: iParticle
    ! permutation of between data order of the request and order in which
    ! data are checked to be present in the request
    integer:: iOrder_I(nVarAvail+nIndexAvail+1)

    character(len=*), parameter:: NameSub = 'get_particle_data_all'
    !----------------------------------------------------------------------

    ! first, determine which data are requested
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.

    ! determine which variables are requested
    call process_request(NameVar, nVarOut, DoReturnVar_V(1:nVarAvail), &
         nIndexOut, DoReturnIndex_I(0:nIndexAvail), iOrder_I)
    ! check that number of data found is the same as requested
    if(nData /= nVarOut + nIndexOut)&
         call CON_stop(NameSub//': incorrect number of data is requested')

    ! return data
    do iParticle = 1, Particle_I(KindReg_)%nParticle
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
    DataOut_VI( pack(iOrder_I,iOrder_I>0), :) = DataOut_VI(:, :)
  end subroutine get_particle_data_all

  !========================================================================

  subroutine get_particle_data_i_particle(nData, NameVar, DataOut_V, iParticle)
    ! the subroutine gets variables specified in the string StringVar
    ! and writes them into DataOut_V;
    ! data is for a single particle specified as iParticle
    integer, intent(in) :: nData
    integer, intent(in) :: iParticle
    real,    intent(out):: DataOut_V(nData)
    character(len=*),intent(in) :: NameVar

    ! mask for returning variables
    logical:: DoReturnVar_V(nVarParticleReg)
    logical:: DoReturnIndex_I(0:nIndexParticleReg)
    ! number of variables/indices found in the request string
    integer:: nVarOut, nIndexOut
    ! permutation of between data order of the request and order in which
    ! data are checked to be present in the request
    integer:: iOrder_I(nVarAvail+nIndexAvail+1)

    character(len=*), parameter:: NameSub = 'get_particle_data_all'
    !----------------------------------------------------------------------
    ! first, determine which data are requested
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.

    ! determine which variables are requested
    call process_request(NameVar, nVarOut, DoReturnVar_V(1:nVarAvail), &
         nIndexOut, DoReturnIndex_I(0:nIndexAvail), iOrder_I)
    ! check that number of data found is the same as requested
    if(nData /= nVarOut + nIndexOut)&
         call CON_stop(NameSub//': incorrect number of data is requested')

    ! return data
    ! store variables
    DataOut_V(1:nVarOut) = PACK(&
         Particle_I(KindReg_)%State_VI( :,iParticle),&
         MASK = DoReturnVar_V)

    ! store indexes
    DataOut_V(nVarOut+1:nData) = PACK(&
         Particle_I(KindReg_)%iIndex_II(:,iParticle),&
         MASK = DoReturnIndex_I)

    ! permute data order so it corresponds to the order of the request
    DataOut_V( pack(iOrder_I,iOrder_I>0)) = DataOut_V(:)

  end subroutine get_particle_data_i_particle
  !==========================================================================

  subroutine write_plot_particle(iFile)

    ! Save particle data

    use ModProcMH,  ONLY: iProc, nProc
    use ModMain,    ONLY: n_step, time_accurate, time_simulation
    use ModIO,      ONLY: &
         StringDateOrTime, NamePlotDir, &
         TypeFile_I, plot_type, plot_form, Plot_
    use ModPlotFile,ONLY: save_plot_file

    integer, intent(in):: iFile

    character(len=100) :: NameFile, NameStart, NameVar
    character (len=80) :: NameProc
    ! container for data
    real, allocatable:: PlotVar_VI(:,:)

    logical:: IsIdl
    integer:: iPlotFile

    character(len=*), parameter :: NameSub = 'write_plot_particle'
    !------------------------------------------------------------------------

    iPlotFile = iFile - Plot_

    ! Set the name of the variables based on plot form
    select case(plot_form(iFile))
    case('tec')
       IsIdl = .false.
       NameVar = '"X", "Y", "Z"'
       NameVar = trim(NameVar)//', "FieldLine"'
    case default
       call CON_stop(NameSub//' ERROR invalid plot form='//plot_form(iFile))
    end select

    ! name of output files
    if(iPlotFile < 10)then
       write(NameStart,'(a,i1,a)') &
            trim(NamePlotDir)//trim(plot_type(iFile))//'_',iPlotFile
    else
       write(NameStart,'(a,i2,a)') &
            trim(NamePlotDir)//trim(plot_type(iFile))//'_',iPlotFile
    end if

    NameFile = NameStart

    if(time_accurate)then
       call get_time_string
       NameFile = trim(NameFile)// "_t"//StringDateOrTime
    end if
    write(NameFile,'(a,i7.7,a)') trim(NameFile) // '_n',n_step

    ! String containing the processor index
    if(nProc < 10000) then
       write(NameProc, '(a,i4.4,a)') "_pe", iProc
    elseif(nProc < 100000) then
       write(NameProc, '(a,i5.5,a)') "_pe", iProc
    else
       write(NameProc, '(a,i6.6,a)') "_pe", iProc
    end if

    NameFile = trim(NameFile) // trim(NameProc)

    if(IsIdl)then
       NameFile = trim(NameFile) // '.out'
    else
       NameFile = trim(NameFile) // '.dat'
    end if

    ! get the data on this processor
    if(allocated(PlotVar_VI)) deallocate(PlotVar_VI)
    allocate(PlotVar_VI(5,n_particle_reg()))
    call get_particle_data(5, 'xx yy zz fl id', PlotVar_VI)

    call save_plot_file(&
         NameFile, &
         TypeFileIn     = TypeFile_I(iFile), &
         StringHeaderIn = "TEMPORARY", &
         TimeIn         = time_simulation, &
         nStepIn        = n_step, &
         NameVarIn      = 'X Y Z FieldLine Index', &
         IsCartesianIn  = .false., &
         CoordIn_I      = PlotVar_VI(1, :), &
         VarIn_VI       = PlotVar_VI(2:,:))

    deallocate(PlotVar_VI)

  end subroutine write_plot_particle

end module ModParticleFieldLine
