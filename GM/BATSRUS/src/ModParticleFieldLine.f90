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
       Particle_I, CellSize_DB, Xyz_DGB, nI, nJ, nK, nBlock, &
       jDim_, kDim_, allocate_particles, Unused_B, &
       message_pass_particles, remove_undefined_particles, &
       mark_undefined
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, B_, Bx_, Bz_
  use ModMain, ONLY: Body1, NameThisComp
  use ModPhysics, ONLY: rBody

  implicit none
  SAVE

  private !except

  public:: read_particle_line_param
  public:: init_particle_line
  public:: extract_particle_line
  public:: add_to_particle_line
  public:: advect_particle_line
  public:: get_particle_data
  public:: set_soft_boundary
  public:: apply_soft_boundary
  public:: n_particle_reg_max
  public:: n_particle_reg


  interface get_particle_data
     module procedure get_particle_data_all
     module procedure get_particle_data_i_particle
  end interface

  

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
       Ds_   = 7, DsOld_ = 8, DsFactor_ = 9, &
       ! direction the tangent to the line
       DirX_    = 10, DirY_    = 11, DirZ_    = 12, &
       DirXOld_ = 13, DirYOld_ = 14, DirZOld_ = 15, &
       ! direction of the magnetic field
       DirBX_    = 16, DirBY_    = 17, DirBZ_    = 18, &
       DirBXOld_ = 19, DirBYOld_ = 20, DirBZOld_ = 21, &
       ! correction to prevent lines from closing
       CorrX_    = 22, CorrY_    = 23, CorrZ_    = 24, &
       ! value of CosBR at the beginning of extraction iteration
       CosBR_    = 25




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
       ! whether to copy the end particle to regular; 
       ! may need to skip copying if step size was changed
       DoCopy_    = 4, &
       ! method used for extraction of field line
       Method_    = 5

  ! extraction methods
  integer, parameter:: RK2_ = 1, GIRARD_ = 2

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
  integer, parameter:: nVarParticleEnd = 25

  ! number of indices
  integer, parameter:: nIndexParticleReg = 2
  integer, parameter:: nIndexParticleEnd = 5

  ! maximum allowed number of field lines
  integer, parameter :: nFieldLineMax = 1000

  ! number of active field lines
  integer:: nFieldLine = 0

  ! soft radial boundary that may be set only once during the run;
  ! "soft" means that particles are allowed beyond this boundary BUT:
  ! - during extraction first particle that crosses it is kept but extraction 
  !   of this line is stopped
  ! - all of the particles beyond the boundary are removed via
  !   public method apply_soft_boundary
  real:: RBoundarySoft = -1.0

  ! logical to switch between integration methods
  logical:: DoTraceGirard = .false.
  !\
  ! Cosine between direction of the field and radial direction: for correcting
  ! the line's direction to prevent it from closing
  !/
  integer, parameter :: nG = 1
  real, allocatable:: CosBR_GB(:,:,:,:), CosBR2_GB(:,:,:,:), &
       GradCosBR_DGB(:,:,:,:,:)

  !\
  ! initialization related info
  !/
  integer:: nLineInit
  real, allocatable:: XyzLineInit_DI(:,:)

contains
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

  !==========================================================================

  subroutine apply_soft_boundary
    ! check which particles went beyond soft boundary and remove them
    !----------------------------------------------------
    if(RBoundarySoft < 0.0) RETURN
    call check_soft_boundary(KindReg_)
    call remove_undefined_particles(KindReg_)
  end subroutine apply_soft_boundary

  !==========================================================================

  subroutine trace_particle_line(iDirTrace)
    use ModCoordTransform, ONLY: cross_product
    ! the subroutine is the tracing loop;
    ! tracing starts at KindEnd_ particles and proceeds in a given direction
    !------------------------------------------------------------------------
    integer, intent(in):: iDirTrace

    integer :: iParticle     ! loop variable

    ! direction of the magnetic field
    real, dimension(MaxDim):: DirBPrev_D, DirBCurr_D, DirBNext_D
    ! direction of the tangent to the line
    real, dimension(MaxDim):: DirPrev_D, DirCurr_D, DirNext_D
    ! correction to prevent line from closing
    real, dimension(MaxDim):: Corr_D
    ! vector used in the Girard's method
    real, dimension(MaxDim):: Omega_D

    ! cosine of angle between direction of field and radial direction
    real:: CosBR, DCosBR

    ! parameters of end particles
    real,    pointer:: StateEnd_VI(:,:)
    integer, pointer:: iIndexEnd_II(:,:)

    ! constants in RK2 method:
    ! dy/dt   = f(t,y)
    ! y_{n+1} = y_n + step * ( beta * k1 + (1-beta)*k2)
    ! k1      = f(t_n, y_n)
    ! k2      = f(t_n, y_n + step * alpha * k1)
    real:: AlphaRK2, BetaRK2
    real, parameter:: cTwoThird = 2.0 / 3.0! alpha in Ralston's method
    real, parameter:: cHalf     = 0.5      ! alpha in midpoint  method
    real, parameter:: cZero     = 0.0      ! beta  in Ralston's method
    real, parameter:: cQuarter  = 0.25     ! beta  in midpoint  method

    character(len=*),parameter:: NameSub = "trace_particle_line"
    !------------------------------------------------------------------------
    ! set pointers to parameters of end particles
    StateEnd_VI  => Particle_I(KindEnd_)%State_VI
    iIndexEnd_II => Particle_I(KindEnd_)%iIndex_II

    ! change const in RK2 based on which method is to be used
    if(DoTraceGirard)then
       AlphaRK2 = cHalf
       BetaRK2  = cZero
    else
       AlphaRK2 = cTwoThird   
       BetaRK2  = cQuarter
    end if

    ! initialize factor to apply to step size to 1
    ! as well as mark a particle for copying
    StateEnd_VI(DsFactor_,1:Particle_I(KindEnd_)%nParticle) = 1.0
    iIndexEnd_II(  DoCopy_,  1:Particle_I(KindEnd_)%nParticle) = 1
    iIndexEnd_II(Method_,  1:Particle_I(KindEnd_)%nParticle) = RK2_

    TRACE: do
       ! check if particles are beyond the soft boundary
       if(RBoundarySoft > 0.0)then
          call check_soft_boundary(KindEnd_)
          call remove_undefined_particles(KindEnd_)
       end if

       ! copy last known coordinates to Auxilary 
       StateEnd_VI(XOld_:ZOld_,1:Particle_I(KindEnd_)%nParticle) = &
            StateEnd_VI(x_:z_, 1:Particle_I(KindEnd_)%nParticle)
       !\
       ! First stage of RK2 method
       !/
       do iParticle = 1, Particle_I(KindEnd_)%nParticle
          select case(iIndexEnd_II(Method_, iParticle))
          case(RK2_)
             call stage1_rk2
          case(GIRARD_)
             call stage1_girard
          end select
       end do
       !\
       ! Message pass: some particles may have moved to different procs
       !/
       call message_pass_particles(KindEnd_)

       !\
       ! remove particles that went outside of the domain
       !/
       if(Body1)then ! check if there is an inner body
          call check_inner_boundary(KindEnd_)
       end if
       call remove_undefined_particles(KindEnd_)

       ! check if all field lines have been completed
       if(is_complete()) RETURN

       !\
       ! Second stage of RK2 method
       !/
       do iParticle = 1, Particle_I(KindEnd_)%nParticle
          select case(iIndexEnd_II(Method_, iParticle))
          case(RK2_)
             call stage2_rk2
          case(GIRARD_)
             call stage2_girard
          end select
       end do
       !\
       ! Message pass: some particles may have moved to different procs
       !/
       call message_pass_particles(KindEnd_)

       !\
       ! remove particles that went outside of the domain
       !/
       if(Body1)then ! check if there is an inner body
          call check_inner_boundary(KindEnd_)
       end if
       call remove_undefined_particles(KindEnd_)


       !\
       ! Stage 3:
       !/
       do iParticle = 1, Particle_I(KindEnd_)%nParticle
          select case(iIndexEnd_II(Method_, iParticle))
          case(RK2_)
             call stage3_rk2
          case(GIRARD_)
             call stage3_girard
          end select
       end do
       !\
       ! Message pass: some particles may have moved to different procs
       !/
       call message_pass_particles(KindEnd_)

       !\
       ! remove particles that went outside of the domain
       !/
       if(Body1)then ! check if there is an inner body
          call check_inner_boundary(KindEnd_)
       end if
       call remove_undefined_particles(KindEnd_)

       ! check if all field lines have been completed
       if(is_complete()) RETURN

       ! increase particle index & copy to regular
       where(iIndexEnd_II(DoCopy_,  1:Particle_I(KindEnd_)%nParticle)==1)
          iIndexEnd_II(id_,1:Particle_I(KindEnd_)%nParticle) = &
               iIndexEnd_II(id_,1:Particle_I(KindEnd_)%nParticle) + &
               iDirTrace
       end where
       call copy_end_to_regular

    end do TRACE

  contains
    !========================================================================
    subroutine stage1_rk2
      ! get the direction of the magnetic field at original location
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBCurr_D)
      ! find the step size
      StateEnd_VI(Ds_, iParticle) = &
           MIN(SpaceStepMax, MAX(SpaceStepMin,&
           0.1*SQRT(&
           sum(CellSize_DB(1:nDim,iIndexEnd_II(0,iParticle))**2)&
           ))) * &
           StateEnd_VI(DsFactor_, iParticle)

      ! save direction for the second stage
      StateEnd_VI(DirBXOld_:DirBZOld_, iParticle) = DirBCurr_D
      DirCurr_D = DirBCurr_D * &
           iDirTrace * iIndexEnd_II(Alignment_, iParticle)

      ! cosine of angle between field and radial direction at beginning of step
      CosBR = &
           sum(StateEnd_VI(x_:z_,iParticle) * DirBCurr_D) / &
           sum(StateEnd_VI(x_:z_,iParticle)**2)**0.5
      StateEnd_VI(CosBR_, iParticle) = CosBR

      ! get middle location
      StateEnd_VI(x_:z_, iParticle) = StateEnd_VI(x_:z_, iParticle) + &
           StateEnd_VI(Ds_, iParticle) * &
           AlphaRK2 * DirCurr_D
    end subroutine stage1_rk2
    !========================================================================
    subroutine stage2_rk2
      ! get the direction of the magnetic field in the middle
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBNext_D)
      ! direction at the 1st stage
      DirBCurr_D = StateEnd_VI(DirBXOld_:DirBZOld_, iParticle)

      DirNext_D = (BetaRK2 * DirBCurr_D + (1.0 - BetaRK2) * DirBNext_D) *&
           iDirTrace * iIndexEnd_II(Alignment_, iParticle)
      ! get final location
      StateEnd_VI(x_:z_,iParticle)=StateEnd_VI(XOld_:ZOld_,iParticle)+&
           StateEnd_VI(Ds_, iParticle) * &
           DirNext_D

      if(DoTraceGirard)then
         ! for the extraction using Girard's method
         ! save direction of the field and curve's tangent at half-step
         StateEnd_VI(DirBX_:DirBZ_, iParticle) = DirBNext_D
         StateEnd_VI(DirX_ :DirZ_,  iParticle) = DirNext_D
      end if
    end subroutine stage2_rk2
    !========================================================================
    subroutine stage3_rk2
      ! get the direction of the magnetic field
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBCurr_D)

      DCosBR = abs(&
           sum(StateEnd_VI(x_:z_, iParticle) * DirBCurr_D) /&
           sum(StateEnd_VI(x_:z_, iParticle)**2)**0.5- &
           StateEnd_VI(CosBR_,iParticle))
      if(DCosBR > 0.05)then
         iIndexEnd_II(DoCopy_, iParticle) = 0
         StateEnd_VI(DsFactor_, iParticle) = &
              0.5 * StateEnd_VI(DsFactor_, iParticle) * 0.05 / DCosBR  
         StateEnd_VI(x_:z_, iParticle) = &
              StateEnd_VI(XOld_:ZOld_, iParticle)
         RETURN
      end if

      iIndexEnd_II(DoCopy_, iParticle) = 1
      StateEnd_VI(DsFactor_, iParticle) = 1.0

      if(DoTraceGirard)then
         ! for the extraction using Girard's method
         ! save direction of the field and curve's tangent at half-step
         StateEnd_VI(DirBXOld_:DirBZOld_, iParticle) = &
              StateEnd_VI(DirBX_:DirBZ_, iParticle)
         StateEnd_VI(DirXOld_ :DirZOld_,  iParticle) = &
              StateEnd_VI(DirX_ :DirZ_,  iParticle)
         StateEnd_VI(DsOld_, iParticle) = StateEnd_VI(Ds_, iParticle)
         iIndexEnd_II(Method_, iParticle) = GIRARD_
      end if
    end subroutine stage3_rk2
    !========================================================================
    subroutine stage1_girard
      ! get the direction of the magnetic field at original location
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBCurr_D)

      ! find the step size
      StateEnd_VI(Ds_, iParticle) = &
           MIN(SpaceStepMax, MAX(SpaceStepMin,&
           0.1*SQRT(&
           sum(CellSize_DB(1:nDim,iIndexEnd_II(0,iParticle))**2)&
           ))) * &
           StateEnd_VI(DsFactor_, iParticle)

      ! direction of tangent and magnetic field at previous half-step
      DirBPrev_D = StateEnd_VI(DirBXOld_:DirBZOld_,iParticle)
      DirPrev_D = StateEnd_VI(DirXOld_ :DirZOld_, iParticle)

      ! get the correction if line starts closing
      CosBR = &
           sum(StateEnd_VI(x_:z_,iParticle) * DirBCurr_D) / &
           sum(StateEnd_VI(x_:z_,iParticle)**2)**0.5
      StateEnd_VI(CosBR_, iParticle) = CosBR

      if(abs(CosBR) < 0.5)then
         call get_grad_cosbr(&
              Xyz_D = StateEnd_VI(x_:z_, iParticle),&
              iBlock=iIndexEnd_II(0,iParticle),&
              Grad_D = Corr_D)
         Corr_D = Corr_D * &
              (1.0/CosBR - 2.0 * sign(1.0, CosBR))
      else
         Corr_D = 0.0
      end if
      StateEnd_VI(CorrX_:CorrZ_, iParticle) = Corr_D


      ! get tangent at the current location
      Omega_D = cross_product(&
           Corr_D * cHalf * StateEnd_VI(DsOld_,iParticle) + &
           (DirBCurr_D-DirBPrev_D), &
           DirPrev_D) * iDirTrace * iIndexEnd_II(Alignment_, iParticle)
      DirCurr_D = DirPrev_D + cross_product(DirPrev_D, Omega_D)
      StateEnd_VI(DirX_:DirZ_, iParticle) = DirCurr_D

      ! get next half-step location
      StateEnd_VI(x_:z_, iParticle) = StateEnd_VI(x_:z_, iParticle) + &
           StateEnd_VI(Ds_, iParticle) * &
           cHalf * DirCurr_D
    end subroutine stage1_girard
    !========================================================================
    subroutine stage2_girard
      ! get the direction of the magnetic field in the middle
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBNext_D)

      ! directions at previous half-step
      DirPrev_D = StateEnd_VI(DirXOld_:DirZOld_, iParticle)
      DirBPrev_D = StateEnd_VI(DirBXOld_:DirBZOld_, iParticle)

      ! tangent to line at the beginning of time-step
      DirCurr_D = StateEnd_VI(DirX_:DirZ_, iParticle)

      ! the correction to prevent line from closing
      Corr_D = StateEnd_VI(CorrX_:CorrZ_, iParticle)

      ! tangent at next half-step
      Omega_D = cross_product(&
           Corr_D*cHalf*( StateEnd_VI(Ds_,iParticle)+StateEnd_VI(DsOld_,iParticle) ) +&
           (DirBNext_D - DirBPrev_D), &
           DirCurr_D) * iDirTrace * iIndexEnd_II(Alignment_, iParticle)
      call rotate_girard(DirPrev_D, Omega_D, DirNext_D)


      StateEnd_VI(DirX_ :DirZ_, iParticle) = DirNext_D
      StateEnd_VI(DirBX_:DirBZ_,iParticle) = DirBNext_D

      ! get final location
      StateEnd_VI(x_:z_,iParticle) = StateEnd_VI(XOld_:ZOld_,iParticle)+&
           StateEnd_VI(Ds_, iParticle) * &
           DirNext_D
    end subroutine stage2_girard
    !========================================================================
    subroutine stage3_girard
      ! get the direction of the magnetic field
      call get_b_dir(&
           Xyz_D = StateEnd_VI(x_:z_, iParticle),&
           iBlock=iIndexEnd_II(0,iParticle),&
           Dir_D = DirBCurr_D)

      CosBR = StateEnd_VI(CosBR_,iParticle)
      if(sum(StateEnd_VI(x_:z_, iParticle) * DirBCurr_D) * CosBR < 0.0)then
         iIndexEnd_II(DoCopy_, iParticle) = 0
         StateEnd_VI(DsFactor_, iParticle) = &
              0.5 * StateEnd_VI(DsFactor_, iParticle)
         StateEnd_VI(x_:z_, iParticle) = &
              StateEnd_VI(XOld_:ZOld_, iParticle)
         RETURN
      end if

      DCosBR = abs(&
           sum(StateEnd_VI(x_:z_, iParticle) * DirBCurr_D) /&
           sum(StateEnd_VI(x_:z_, iParticle)**2)**0.5- &
           CosBR)
      if(DCosBR > 0.05)then
         iIndexEnd_II(DoCopy_, iParticle) = 0
         StateEnd_VI(DsFactor_, iParticle) = &
              0.5 * StateEnd_VI(DsFactor_, iParticle) * 0.05 / DCosBR  
         StateEnd_VI(x_:z_, iParticle) = &
              StateEnd_VI(XOld_:ZOld_, iParticle)
         RETURN
      end if

      iIndexEnd_II(DoCopy_, iParticle) = 1
      StateEnd_VI(DsFactor_, iParticle) = 1.0

      ! store for the next iteration
      StateEnd_VI(DsOld_, iParticle) = StateEnd_VI(Ds_, iParticle)
      StateEnd_VI(DirXOld_:DirZOld_, iParticle) = &
           StateEnd_VI(DirX_:DirZ_, iParticle)
      StateEnd_VI(DirBXOld_:DirBZOld_, iParticle) = &
           StateEnd_VI(DirBX_:DirBZ_, iParticle)
    end subroutine stage3_girard
    !========================================================================
    function is_complete() result(IsCompleteOut)
      use ModMpi
      use BATL_mpi, ONLY: iComm
      ! returns true if tracing is complete for all field line on all procs
      logical:: IsComplete, IsCompleteOut
      integer:: iError
      !----------------------------------------------------------------------
      if(Particle_I(KindEnd_)%nParticle < 0)&
           call CON_stop(NameThisComp//":"//NameSub//&
           ": negative number of particles of KindEnd_")
      IsComplete = Particle_I(KindEnd_)%nParticle == 0
      ! reduce IsComplete variable from all processors
      call MPI_Allreduce(IsComplete, IsCompleteOut, 1, &
           MPI_LOGICAL, MPI_LAND, iComm, iError)               
    end function is_complete
    !========================================================================
    subroutine rotate_girard(DirIn_D, Omega_D, DirOut_D)
      ! apply Girard's method for integrating
      !----------------------------------------------
      real, intent(in) :: DirIn_D(MaxDim)
      real, intent(in) :: Omega_D(MaxDim)
      real, intent(out):: DirOut_D(MaxDim)

      real:: DirAux_D(MaxDim)
      !----------------------------------------------------------------------
      DirAux_D = DirIn_D + 0.50*cross_product(DirIn_D, Omega_D)
      DirOut_D = DirIn_D + &
           cross_product(DirAux_D, Omega_D)/ (1 + 0.25*sum(Omega_D**2))
    end subroutine rotate_girard
  end subroutine trace_particle_line
  
  !==========================================================================

  subroutine extract_particle_line(nFieldLineIn, XyzStart_DI, iTraceModeIn, &
       iIndexStart_II,&
       UseInputInGenCoord)
    use BATL_geometry, ONLY: coord_to_xyz
    ! extract nFieldLineIn magnetic field lines starting at XyzStart_DI;
    ! the whole field lines are extracted, i.e. they are traced forward
    ! and backward up until it reaches boundaries of the domain;
    ! requested coordinates may be different for different processor,
    ! if a certain line can't be started on a given processor, it is
    ! ignored, thus duplicates are avoided
    !------------------------------------------------------------------------
    integer,         intent(in)::nFieldLineIn
    real,            intent(in)::XyzStart_DI(MaxDim, nFieldLineIn)
    ! mode of tracing (forward, backward or both ways)
    integer,optional,intent(in)::iTraceModeIn
    ! initial particle indices for starting particles
    integer,optional,intent(in)::iIndexStart_II(nIndexParticleReg,nFieldLineIn)
    !\
    ! An input can be in generalized coordinates
    !/
    logical, optional,intent(in) :: UseInputInGenCoord 
    logical :: DoTransformCoordToXyz
    integer :: nLineThisProc ! number of new field lines initialized locally
    integer :: nLineAllProc  ! number of new field lines initialized globally
    integer :: iFieldLine    ! loop variable
    integer :: nParticleOld  ! number of already existing regular particles

    ! direction of tracing: -1 -> backward, +1 -> forward
    integer:: iDirTrace

    ! mode of tracing (see description below)
    integer:: iTraceMode

    ! starting point of a field line 
    ! (to satisfy intent INOUT for check_interpolate_amr_gc)
    real:: XyzStart_D(MaxDim)
    ! its index along the field line
    integer:: iIndexStart_I(nIndexParticleReg)

    ! parameters of end particles
    real,    pointer:: StateEnd_VI(:,:)
    integer, pointer:: iIndexEnd_II(:,:)
    ! parameters of regular particles
    real,    pointer:: StateReg_VI(:,:)
    integer, pointer:: iIndexReg_II(:,:)

    character(len=*), parameter:: NameSub='extract_particle_line'
    !------------------------------------------------------------------------
    ! set pointers to parameters of end particles
    StateEnd_VI  => Particle_I(KindEnd_)%State_VI
    iIndexEnd_II => Particle_I(KindEnd_)%iIndex_II

    ! set pointers to parameters of regular particles
    StateReg_VI  => Particle_I(KindReg_)%State_VI
    iIndexReg_II => Particle_I(KindReg_)%iIndex_II
    DoTransformCoordToXyz = .false.
    if(present(UseInputInGenCoord))&
         DoTransformCoordToXyz = UseInputInGenCoord

    !\
    ! allocate containers for cosine of angle between B and radial direction
    !/
    if(DoTraceGirard)then
       allocate(CosBR_GB(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
            1-nG*kDim_:nK+nG*kDim_,nBlock))
       CosBR_GB = 0.0
       allocate(CosBR2_GB(1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
            1-nG*kDim_:nK+nG*kDim_,nBlock))
       CosBR2_GB = 0.0
       allocate(GradCosBR_DGB(nDim,1-nG:nI+nG,1-nG*jDim_:nJ+nG*jDim_,&
            1-nG*kDim_:nK+nG*kDim_,nBlock))
       GradCosBR_DGB = 0.0
       call compute_grad_cosbr
    end if

    !\
    ! Trace field lines
    !/
    ! initialize field lines
    Particle_I(KindEnd_)%nParticle = 0
    nLineThisProc = 0
    nParticleOld  = Particle_I(KindReg_)%nParticle
    do iFieldLine = 1, nFieldLineIn
       if(DoTransformCoordToXyz)then
          call coord_to_xyz(XyzStart_DI(:, iFieldLine), XyzStart_D)
       else
          XyzStart_D = XyzStart_DI(:, iFieldLine) 
       end if
       if(present(iIndexStart_II)) then
          iIndexStart_I = iIndexStart_II(:,iFieldLine)
       else
          iIndexStart_I = (/nFieldLine + iFieldLine, 0/)
       end if
       call start_line(XyzStart_D, iIndexStart_I)
    end do

    ! how many lines have been started on all processors
    call count_new_lines()

    nFieldLine    = nFieldLine + nLineAllProc
    if(nFieldLine > nFieldLineMax)&
         call CON_stop(NameThisComp//':'//NameSub//&
         ': Limit for number of particle field lines exceeded')
    call copy_end_to_regular

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

       call trace_particle_line(iDirTrace)

       ! if just finished backward tracing and need to trace in both dirs
       ! => return to initial particles
       if(iDirTrace < 0 .and. iTraceMode == 0)then
          ! copy initial points to KindEnd_:
          ! the initial particles are currently right after the particles,
          ! that were in the list before the start of this subroutine,
          ! i.e. occupy positions from (nParticleOld+1)
          StateEnd_VI(x_:z_, 1:nLineThisProc) = &
               StateReg_VI(x_:z_, nParticleOld+1:nParticleOld+nLineThisProc)
          iIndexEnd_II(0:id_,1:nLineThisProc) = &
               iIndexReg_II(0:id_,nParticleOld+1:nParticleOld+nLineThisProc)
          Particle_I(KindEnd_)%nParticle = nLineThisProc
       end if
    end do
    !\
    ! Offset in id_
    !/
    call offset_id(nFieldLine - nLineAllProc+1,nFieldLine)

    !\
    ! free space
    if(DoTraceGirard)&
       deallocate(CosBR_GB, CosBR2_GB, GradCosBR_DGB)

  contains

    subroutine start_line(XyzStart_D, iIndexStart_I)
      real,    intent(inout):: XyzStart_D(MaxDim)
      integer, intent(in)   :: iIndexStart_I(nIndexParticleReg)

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

      ! Exclude particles that are in the zero magnetic field region
      if(do_exclude(XyzStart_D, iBlockOut))&
           RETURN

      Particle_I(KindEnd_)%nParticle = &
           Particle_I(KindEnd_)%nParticle + 1         
      nLineThisProc = nLineThisProc + 1

      StateEnd_VI(x_:z_,Particle_I(KindEnd_)%nParticle) = XyzStart_D
      iIndexEnd_II(fl_, Particle_I(KindEnd_)%nParticle) = iIndexStart_I(fl_)
      iIndexEnd_II(id_, Particle_I(KindEnd_)%nParticle) = iIndexStart_I(id_)
      iIndexEnd_II(0,   Particle_I(KindEnd_)%nParticle) = iBlockOut
    end subroutine start_line
    !=============
    subroutine  offset_id(iLineStartIn,iLineEndIn)
      use ModMpi
      use BATL_mpi, ONLY: iComm, nProc
      integer, optional, intent(in) :: iLineStartIn, iLineEndIn
      integer:: iLineStart, iLineEnd, iParticle, nParticle, iLine, iError
      integer:: iOffsetLocal_I(1:nFieldLineMax), iOffset_I(1:nFieldLineMax)
      integer, pointer:: iIndex_II(:,:)
      !-----------------------------
      iLineStart = 1 
      if(present(iLineStartIn))iLineStart = iLineStartIn
      iLineEnd = nFieldLineMax
      if(present(iLineEndIn))iLineEnd = iLineEndIn
      iOffsetLocal_I = 0; iOffset_I = 0  
      !\
      ! set a pointer to parameters of regular particles
      !/
      iIndex_II => Particle_I(KindReg_)%iIndex_II
      nParticle = Particle_I(KindReg_)%nParticle
      do iParticle = 1, nParticle
         iLine = iIndex_II(fl_,iParticle)
         !\
         ! Reset offset, if id_ is less than one
         !/
         iOffsetLocal_I(iLine) = max(iOffsetLocal_I(iLine),&
              1 - iIndex_II(id_,iParticle))
      end do
      ! find maximim offset from all processors
      if(nProc > 1) then
         call MPI_Allreduce(iOffsetLocal_I(iLineStart), &
              iOffset_I(iLineStart), 1 + iLineEnd - iLineStart, &
              MPI_INTEGER, MPI_MAX, iComm, iError)
      else
         iOffset_I(iLineStart:iLineEnd) = &
              iOffset_I(iLineStart:iLineEnd)
      end if
      do iParticle = 1, nParticle
         iLine = iIndex_II(fl_,iParticle)
         !\
         ! Apply offset
         !/
         iIndex_II(id_,iParticle) = iIndex_II(id_,iParticle) + &
         iOffsetLocal_I(iLine)
      end do
    end subroutine offset_id
    !========================================================================
    function do_exclude(Xyz_D, iBlock) result(DoExcludeOut)
      use ModMain, ONLY: UseB0
      use ModB0, ONLY: get_b0
      ! determine if the starting point and therefore the whole field line
      ! should be excluded; criterion: zero magnetic field
      real,    intent(in):: Xyz_D(MaxDim)
      integer, intent(in):: iBlock
      logical:: DoExcludeOut

      ! interpolated magnetic field
      real   :: B_D(MaxDim) = 0.0
      ! interpolation data: number of cells, cell indices, weights
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim)
      !----------------------------------------------------------------------
      B_D = 0
      ! get potential part of the magnetic field at the current location
      if(UseB0)call get_b0(Xyz_D, B_D)
      ! get the remaining part of the magnetic field
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
      ! interpolate magnetic field value
      do iCell = 1, nCell
         i_D = 1
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         B_D = B_D + &
              State_VGB(Bx_:Bz_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
      end do
      DoExcludeOut = all(B_D==0)
    end function do_exclude
    !========================================================================
    subroutine count_new_lines()
      ! gather information from all processors on how many new field lines
      ! have been started
      use ModMpi
      use BATL_mpi, ONLY: iComm
      integer:: iError
      !------------------------------------------------------------------------
      call MPI_Allreduce(nLineThisProc, nLineAllProc, 1, &
           MPI_INTEGER, MPI_SUM, iComm, iError)
    end subroutine count_new_lines
    !========================================================================
    subroutine get_alignment()
      ! determine alignment of particle indexing with direction 
      ! of the magnetic field
      integer:: iParticle
      real:: Dir_D(MaxDim)
      !------------------------------------------------------------------------
      if(iOrderMode == Field_)then
         iIndexEnd_II(Alignment_, 1:Particle_I(KindEnd_)%nParticle) = 1
         RETURN
      end if

      do iParticle = 1, Particle_I(KindEnd_)%nParticle
         call get_b_dir(&
              Xyz_D = StateEnd_VI(x_:z_, iParticle),&
              iBlock=iIndexEnd_II(0,iParticle),&
              Dir_D = Dir_D)
         iIndexEnd_II(Alignment_, iParticle) = &
              nint( SIGN(1.0, sum(Dir_D*StateEnd_VI(x_:z_,iParticle))) )
      end do
    end subroutine get_alignment
    !========================================================================
    subroutine compute_grad_cosbr()
      use BATL_pass_cell, ONLY: message_pass_cell
      use ModMain, ONLY: UseB0
      use ModB0, ONLY: B0_DGB
      use ModGeometry, ONLY: R_BLK
      use ModCellGradient, ONLY: calc_gradient
      use ModFaceGradient, ONLY: set_block_field2
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
            CosBR_GB(i,j,k,iBlock) = sum(Bcell_D*XyzCell_D) / &
                 (sqrt(sum(BCell_D**2))*R_BLK(i,j,k,iBlock))
         end do; end do; end do
      end do
      ! fill ghost celss via message pass
      call message_pass_cell(&
           nG              = 1, &
           State_GB        = CosBR_GB, &
           nProlongOrderIn = 1)
      
      do iBlock = 1, nBlock
         if(Unused_B(iBlock))CYCLE
         ! correct ghost cells
         call set_block_field2(&
              iBlock    = iBlock, &
              nVar      = 1, &
              Field1_VG = CosBR2_GB, &
              Field_VG  = CosBR_GB)
         ! compute gradient
         call calc_gradient(iBlock, CosBR_GB(:,:,:,iBlock), nG, &
              GradCosBR_DGB(:,:,:,:,iBlock))
      end do
      ! fill ghost cells for gradient via message pass
      call message_pass_cell(&
           nVar            = nDim,&
           nG              = 1, &
           State_VGB       = GradCosBR_DGB, &
           nProlongOrderIn = 1)
    end subroutine compute_grad_cosbr
  end subroutine extract_particle_line

  !========================================================================
  subroutine copy_end_to_regular
    ! copies indicated variables of known end particles to regular particles
    integer, parameter:: iVarCopy_I(3)   = (/x_, y_, z_/)
    integer, parameter:: iIndexCopy_I(3) = (/0, fl_, id_/)
    integer:: iParticle
    !----------------------------------------------------------------------
    do iParticle = 1, Particle_I(KindEnd_)%nParticle
       if(Particle_I(KindEnd_)%iIndex_II(DoCopy_,iParticle)==0)&
            CYCLE
       Particle_I(KindReg_)%nParticle = Particle_I(KindReg_)%nParticle+1
       Particle_I(KindReg_)%State_VI(iVarCopy_I, &
            Particle_I(KindReg_)%nParticle) =&
            Particle_I(KindEnd_)%State_VI(iVarCopy_I, iParticle)
       Particle_I(KindReg_)%iIndex_II(iIndexCopy_I,&
            Particle_I(KindReg_)%nParticle) =&
            Particle_I(KindEnd_)%iIndex_II(iIndexCopy_I, iParticle)
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
    call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
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
  end subroutine get_b_dir
  !========================================================================
  subroutine get_grad_cosbr(Xyz_D, iBlock, Grad_D)
    use ModMain, ONLY: UseB0
    use ModB0, ONLY: get_b0
    ! returns the direction of magnetic field 
    ! as well as the block used for interpolation
    real,          intent(in) :: Xyz_D(MaxDim)
    integer,       intent(in) :: iBlock
    real,          intent(out):: Grad_D(MaxDim)

    ! interpolation data: number of cells, cell indices, weights
    integer:: nCell, iCell_II(0:nDim, 2**nDim)
    real   :: Weight_I(2**nDim)
    integer:: iCell ! loop variable
    integer:: i_D(MaxDim)
    character(len=200):: StringError
    character(len=*), parameter:: NameSub = "get_grad_cosbr"
    !----------------------------------------------------------------------
    Grad_D = 0
    ! get the remaining part of the magnetic field
    call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I)
    ! interpolate magnetic field value
    do iCell = 1, nCell
       i_D = 1
       i_D(1:nDim) = iCell_II(1:nDim, iCell)
       Grad_D(1:nDim) = Grad_D(1:nDim) + &
            GradCosBR_DGB(:,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
    end do
  end subroutine get_grad_cosbr
  !========================================================================
  subroutine add_to_particle_line(nParticleIn, XyzIn_DI, iIndexIn_II,&
       UseInputInGenCoord, DoReplace)
    !\
    ! Very strange, but this is the first use of the BATL procedure to 
    ! transform the generalized to cartesian coords in BATSRUS
    !/
    use BATL_geometry, ONLY: coord_to_xyz
    ! add particles with specified coordinates to the already existing lines
    integer, intent(in):: nParticleIn
    real,    intent(in):: XyzIn_DI(MaxDim, nParticleIn)
    integer, intent(in):: iIndexIn_II(nIndexParticleReg, nParticleIn)
    !\
    ! An input can be in generalized coordinates
    !/
    logical, optional,intent(in) :: UseInputInGenCoord 
    !\
    ! Whether to replace ALL old particles with the input
    !/
    logical, optional, intent(in):: DoReplace
    logical :: DoTransformCoordToXyz
    real   :: Xyz_D(MaxDim)
    integer:: iIndex_I(nIndexParticleReg)
    integer:: iParticle
    integer:: iProcOut, iBlockOut
    ! parameters of regular particles
    real,    pointer:: StateReg_VI(:,:)
    integer, pointer:: iIndexReg_II(:,:)

    character(len=100):: StringError
    character(len=*), parameter:: NameSub = 'add_to_particle_line'
    !----------------------------------------------------------------------
    StateReg_VI => Particle_I(KindReg_)%State_VI
    iIndexReg_II=> Particle_I(KindReg_)%iIndex_II
    DoTransformCoordToXyz = .false.
    if(present(UseInputInGenCoord))&
         DoTransformCoordToXyz = UseInputInGenCoord
    if(present(DoReplace))then
       if(DoReplace) Particle_I(KindReg_)%nParticle = 0
    end if
    do iParticle = 1, nParticleIn
       if(DoTransformCoordToXyz)then
          call coord_to_xyz(XyzIn_DI(:,   iParticle), Xyz_D)
       else
          Xyz_D   = XyzIn_DI(:,   iParticle)
       end if
       iIndex_I= iIndexIn_II(:,iParticle)

       ! find block and processor suitable for interpolation
       call check_interpolate_amr_gc(Xyz_D, &
            1, & ! input block ID doesn't matter
            iProcOut, iBlockOut)

       ! check whether point is outside of the domain
       if(iProcOut < 0)then
          write(StringError,'(a,es15.6,a, es15.6, a, es15.6)') &
               "Point for a field line is outside of the domain: X = ",&
               Xyz_D(1), " Y = ", Xyz_D(2), " Z = ", Xyz_D(3) 
          call stop_mpi(NameThisComp//':'//NameSub //": "//StringError)
       end if

       ! Assign particle to an appropriate processor
       if(iProc /= iProcOut) RETURN

       Particle_I(KindReg_)%nParticle = Particle_I(KindReg_)%nParticle + 1

       StateReg_VI(x_:z_,Particle_I(KindReg_)%nParticle) = Xyz_D
       iIndexReg_II(fl_, Particle_I(KindReg_)%nParticle) = iIndex_I(fl_)
       iIndexReg_II(id_, Particle_I(KindReg_)%nParticle) = iIndex_I(id_)
       iIndexReg_II(0,   Particle_I(KindReg_)%nParticle) = iBlockOut
    end do
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
      integer:: i_D(MaxDim)
      !------------------------------------
      ! reset the interpoalted values
      V_D = 0!; Rho = 0; M_D = 0
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
      ! check of there is plasma at the location
      !      if(Rho==0)&
      !           call CON_stop(NameSub//&
      !           ': trying to advect particle line at region with no plasma')
      !      V_D = M_D / Rho
    end subroutine get_v
  end subroutine advect_particle_line

  !========================================================================

  subroutine check_inner_boundary(iKind)
    ! check whether particles have crossed an inner boundary
    ! if so => set particle's block to negative value
    integer, intent(in):: iKind

    integer:: iParticle ! loop variable
    !----------------------------------------------------------------------
    do iParticle = 1, Particle_I(iKind)%nParticle
       if(sum(Particle_I(iKind)%State_VI(x_:z_,iParticle)**2) < rBody*rBody)&
            call mark_undefined(iKind, iParticle)
    end do
  end subroutine check_inner_boundary

  !========================================================================

  subroutine check_soft_boundary(iKind)
    ! check whether particles have crossed an inner boundary
    ! if so => set particle's block to negative value
    integer, intent(in):: iKind

    integer:: iParticle ! loop variable
    !----------------------------------------------------------------------
    do iParticle = 1, Particle_I(iKind)%nParticle
       if(sum(Particle_I(iKind)%State_VI(x_:z_,iParticle)**2) > rBoundarySoft**2)&
            call mark_undefined(iKind, iParticle)
    end do
  end subroutine check_soft_boundary

  !========================================================================
  
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
    !\
    ! first, determine which data are requested
    !/
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.
    ! determine which variables are requested
    call process_request(NameVar, nVarOut, DoReturnVar_V(1:nVarAvail), &
         nIndexOut, DoReturnIndex_I(0:nIndexAvail), iOrder_I)
    ! check that number of data found is the same as requested
    if(nData /= nVarOut + nIndexOut)&
         call CON_stop(NameSub//': incorrect number of data is requested')

    !\
    ! return data
    !/
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
    !\
    ! first, determine which data are requested
    !/
    DoReturnVar_V   = .false.
    DoReturnIndex_I = .false.
    ! determine which variables are requested
    call process_request(NameVar, nVarOut, DoReturnVar_V(1:nVarAvail), &
         nIndexOut, DoReturnIndex_I(0:nIndexAvail), iOrder_I)
    ! check that number of data found is the same as requested
    if(nData /= nVarOut + nIndexOut)&
         call CON_stop(NameSub//': incorrect number of data is requested')

    !\
    ! return data
    !/
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

end module ModParticleFieldLine
