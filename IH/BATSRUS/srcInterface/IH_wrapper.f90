!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module IH_wrapper

  ! Wrapper for IH_BATSRUS Inner Heliosphere (IH) component

  implicit none

  private ! except

  ! CON wrapper
  public:: IH_set_param
  public:: IH_init_session
  public:: IH_run
  public:: IH_save_restart
  public:: IH_finalize

  ! Global buffer coupling
  public:: IH_get_for_global_buffer

  ! Coupling toolkit
  public:: IH_synchronize_refinement
  public:: IH_get_for_mh
  public:: IH_get_for_mh_with_xyz
  public:: IH_put_from_mh

  ! Coupling with SC
  public:: IH_set_buffer_grid
  public:: IH_set_buffer_grid_get_info
  public:: IH_save_global_buffer
  public:: IH_match_ibc

  ! Coupling with SP
  public:: IH_get_for_sp
  public:: IH_get_a_line_point

  ! Coupling with GM
  public:: IH_get_for_gm

contains
  !==========================================================================

  subroutine IH_init_session(iSession, TimeSimulation)

    use IH_ModMain,     ONLY: Time_Simulation
    use CON_physics, ONLY: get_time

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_init_session'

    logical :: IsUninitialized = .true.
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)

    if(IsUninitialized)then

       call get_time(tSimulationOut=Time_Simulation)

       call IH_BATS_setup
       IsUninitialized = .false.
    end if
    call IH_BATS_init_session

    if(DoTest)write(*,*)NameSub,' finished for session ',iSession

  end subroutine IH_init_session

  !==========================================================================
  subroutine IH_set_param(CompInfo, TypeAction)

    use CON_comp_info
    use IH_ModProcMH
    use IH_ModIO, ONLY: iUnitOut, StringPrefix, STDOUT_, NamePlotDir
    use IH_ModRestartFile, ONLY: NameRestartInDir, NameRestartOutDir
    use IH_ModMain, ONLY : CodeVersion, NameThisComp, &
         time_accurate, StartTime, iStartTime_I
    use CON_physics, ONLY: get_time
    use ModTimeConvert, ONLY: time_real_to_int

    character (len=*), parameter :: NameSub='IH_set_param'

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do

    logical :: DoTest,DoTestMe
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with TypeAction,iProc=',&
         TypeAction,iProc

    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true.,                        &
            NameVersion='IH_BATSRUS (Univ. of Michigan)', &
            Version    =CodeVersion)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc,&
            Name=NameThisComp)

       NamePlotDir(1:2)      = NameThisComp
       NameRestartInDir(1:2) = NameThisComp
       NameRestartOutDir(1:2)= NameThisComp
    case('READ','CHECK')
       call get_time( &
            DoTimeAccurateOut = time_accurate, &
            tStartOut         = StartTime)
       call time_real_to_int(StartTime,iStartTime_I)

       call IH_set_parameters(TypeAction)
    case('STDOUT')
       iUnitOut=STDOUT_
       if(iProc==0)then
          StringPrefix = NameThisComp//':'
       else
          write(StringPrefix,'(a,i4.4,a)')NameThisComp,iProc,':'
       end if
    case('FILEOUT')
       call get(CompInfo,iUnitOut=iUnitOut)
       StringPrefix=''
    case('GRID')
       call IH_set_grid
    case default
       call CON_stop(NameSub//' SWMF_ERROR: invalid TypeAction='//TypeAction)
    end select
  end subroutine IH_set_param

  !============================================================================

  subroutine IH_finalize(TimeSimulation)

    use IH_ModMain, ONLY: time_loop

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_finalize'

    integer :: iError
    !--------------------------------------------------------------------------
    ! We are not advancing in time any longer
    time_loop = .false.

    call IH_BATS_save_files('FINAL')

    call IH_error_report('PRINT',0.,iError,.true.)

  end subroutine IH_finalize

  !============================================================================

  subroutine IH_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='IH_save_restart'

    call IH_BATS_save_files('RESTART')

  end subroutine IH_save_restart

  !============================================================================

  subroutine IH_run(TimeSimulation,TimeSimulationLimit)

    use IH_ModProcMH, ONLY: iProc
    use IH_ModMain, ONLY: Time_Simulation

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    character(len=*), parameter :: NameSub='IH_run'

    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' called with tSim, tSimLimit, iProc=',&
         TimeSimulation, TimeSimulationLimit, iProc

    if(abs(Time_Simulation-TimeSimulation)>0.0001) then
       write(*,*)NameSub, &
            ' IH time=',Time_Simulation,' SWMF time=',TimeSimulation
       call CON_stop(NameSub//': IH and SWMF simulation times differ')
    end if

    call IH_BATS_advance(TimeSimulationLimit)

    ! Return time after the time step
    TimeSimulation = Time_Simulation

  end subroutine IH_run

  !============================================================================
  subroutine IH_set_buffer_grid(DD)
    use IH_ModBuffer,ONLY:&
         set_spher_buffer_grid,set_buffer_name,&
         DomainDecompositionType,&
         LocalBufferDD
    use CON_coupler,ONLY:IH_,is_proc

    type(DomainDecompositionType),&
         intent(out)::DD

    call set_spher_buffer_grid(&
         DD,IH_,IsLocal=.false.)
    if(.not.is_proc(IH_))return

    call set_spher_buffer_grid(&
         LocalBufferDD,IH_,IsLocal=.true.)
    call set_buffer_name('IH_from_sc')

  end subroutine IH_set_buffer_grid
  !============================================================================
  subroutine IH_set_grid
    use CON_comp_param
    use IH_domain_decomposition
    use CON_coupler
    use IH_ModMain, ONLY: TypeCoordSystem, nVar, NameVarCouple
    use IH_ModPhysics,ONLY:No2Si_V, UnitX_
    use IH_ModGeometry, ONLY : TypeGeometry, LogRGen_I
    !--------------------------------------------------------------------------

    ! Here we should set the IH (MH) grid descriptor
    if(done_dd_init(IH_))return
    call init_decomposition(&
         GridID_=IH_,&
         CompID_=IH_,&
         nDim=3,     &
         IsTreeDecomposition=.true.)

    ! Allocate array on non_IH processors (???why is this needed???)
    if(.not.allocated(LogRGen_I))then
       allocate(LogRGen_I(1))
       LogRGen_I = 0.0
    end if

    call set_coord_system(&
         GridID_=IH_,&
         TypeCoord=TypeCoordSystem,&
         UnitX=No2Si_V(UnitX_),&
         TypeGeometry=TypeGeometry,&
         Coord1_I=LogRGen_I, &
         nVar = nVar, &
         NameVar = NameVarCouple)

    if(is_Proc(IH_))then
       !Initialize the local grid

       call init_decomposition(&
            DomainDecomposition=MH_DomainDecomposition,&
            CompID_=IH_,&
            nDim=3,&
            IsTreeDecomposition=.true.)

       !Get the octree root array
       call MH_get_root_decomposition(MH_DomainDecomposition)

       !Get the whole octree after the initial refinement
       call MH_update_local_decomposition(MH_DomainDecomposition)

       MH_DomainDecomposition%IsLocal=.true.
    end if

    !Repeat the initialization at the global grid level:
    !Octree root array:
    if(is_proc0(IH_))call MH_get_root_decomposition(IH_)

    !Broadcast root array:
    call bcast_decomposition(IH_)

    !Synchronize global and local grids:
    call synchronize_refinement(&
         GridID_=IH_,&
         localDD=MH_domaindecomposition)

  end subroutine IH_set_grid

  !===============================================================

  subroutine IH_synchronize_refinement(iProc0,iCommUnion)

    use IH_domain_decomposition
    use CON_comp_param

    integer, intent(in) ::iProc0,iCommUnion

    !Synchronize the local grid decomposition to accomodate the
    !grid change.

    if(is_proc(IH_)) &
         call MH_update_local_decomposition(MH_DomainDecomposition)

    call synchronize_refinement(&
         GridID_=IH_,&
         LocalDD=MH_domaindecomposition,&
         iProcUnion=iProc0,&
         iCommUnion=iCommUnion)

  end subroutine IH_synchronize_refinement

  !============================================================================

  subroutine IH_set_buffer_grid_get_info(CompID_, &
       nR, nPhi, nTheta, BufferMinMax_DI)

    use IH_domain_decomposition, ONLY: is_proc
    use IH_ModMain,              ONLY: BuffR_, nPhiBuff, nThetaBuff,&
         nRBuff, BufferMin_D, BufferMax_D, dSphBuff_D

    integer, intent(in)     :: CompID_
    integer, intent(out)    :: nR, nPhi, nTheta
    real, intent(out)       :: BufferMinMax_DI(3,2)

    integer  :: nCell_D(3)
    logical :: DoTest, DoTestMe

    character(len=*), parameter :: NameSub = 'IH_set_buffer_grid_get_info'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest, DoTestMe)

    ! Make sure only a coupling target component is executing this 
    if(.not. is_proc(CompID_)) RETURN

    ! Return buffer size and limits to SWMF calling routine
    BufferMinMax_DI(:,1) = BufferMin_D
    BufferMinMax_DI(:,2) = BufferMax_D

    nR     = nRBuff 
    nPhi   = nPhiBuff
    nTheta = nThetaBuff

    ! Calculate grid spacing and save in IH_BATSRUS
    nCell_D = (/nR, nPhi, nTheta/)
    dSphBuff_D = (BufferMax_D - BufferMin_D)/real(nCell_D)
    dSphBuff_D(BuffR_) = (BufferMax_D(BuffR_) - BufferMin_D(BuffR_)) &
         /real(nCell_D(BuffR_) - 1)

    if(DoTest) then
       write(*,*) NameSub,': with nR, nPhi, nTheta = ',nCell_D
       write(*,*) 'BufferMin_D: ',BufferMin_D
       write(*,*) 'BufferMax_D: ',BufferMax_D
       write(*,*) 'dSph_D: ',dSphBuff_D
    end if

  end subroutine IH_set_buffer_grid_get_info

  !============================================================================

  subroutine IH_save_global_buffer(nVar, nR, nPhi, nTheta,BufferIn_VG)

    use IH_ModMain, ONLY: BufferState_VG

    integer,intent(in) :: nVar, nR, nPhi, nTheta
    real,intent(in)    :: BufferIn_VG(nVar,nR,0:nPhi+1,0:nTheta+1)

    character(len=*), parameter :: NameSub = 'IH_save_global_buffer'
    !-------------------------------------------------------------
    if(.not. allocated(BufferState_VG))&
         allocate(BufferState_VG(nVar, nR, 0:nPhi+1, 0:nTheta+1))
    BufferState_VG = BufferIn_VG

  end subroutine IH_save_global_buffer

  !===========================================================================
  subroutine IH_match_ibc

    use IH_ModMessagePass, ONLY: exchange_messages, fill_in_from_buffer
    use IH_ModGeometry,ONLY:R_BLK
    use IH_BATL_lib,  ONLY: Xyz_DGB
    use IH_ModMain,   ONLY:&
         nI,nJ,nK, BufferMax_D, MaxDim,nBlock, Unused_B
    use IH_ModAdvance,ONLY:nVar,State_VGB,rho_,rhoUx_,rhoUz_,Ux_,Uz_
    use IH_ModProcMH, ONLY:iProc
    use IH_ModIO,     ONLY:IsRestartCoupler

    character(len=*), parameter :: NameSub='IH_match_ibc'
    character(len=*), parameter :: StringTest ='IH_fill_buffer_only'

    integer  :: iBlock
    integer  :: i,j,k
    real     :: x_D(MaxDim), rBuffMax
    logical  :: DoTest,DoTestMe
    ! ------------------------------------------------------------------------
    if(IsRestartCoupler) RETURN

    rBuffMax = BufferMax_D(1)

    call CON_set_do_test(StringTest, DoTest, DoTestMe)
    if(DoTest .and. iProc == 0)  write(*,*) &
         NameSub,' in test mode: no filling of cells outside the buffer grid.'

    ! Fill all spatial domain with values depend on the BC
    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE

       ! Fill in the cells, covered by the bufer grid, including ghost cells
       call fill_in_from_buffer(iBlock)

       ! Fill in the physical cells, which are outside the buffer grid
       ! When testing, do not fill cells outside the buffer
       if(.not. DoTest) then 
          do k = 1, nK; do j = 1 , nJ; do i = 1, nI
             if(R_BLK(i,j,k,iBlock) < rBuffMax)CYCLE

             ! For each grid point, get the values at the base (buffer) 
             x_D = Xyz_DGB(:,i,j,k,iBlock)*rBuffMax/R_BLK(i,j,k,iBlock)

             ! The grid point values are extracted from the base values
             call IH_get_from_spher_buffer_grid(&
                  x_D, nVar, State_VGB(:,i,j,k,iBlock))

             !Transform primitive variables to conservative ones:
             State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)=&
                  State_VGB(Ux_:Uz_,i,j,k,iBlock)*&
                  State_VGB(rho_,i,j,k,iBlock)

             !Scale as (r/R)^2:
             State_VGB(:,i,j,k,iBlock)=&
                  State_VGB(:,i,j,k,iBlock)*&
                  (rBuffMax/R_BLK(i,j,k,iBlock))**2

          end do; end do; end do
       end if
    end do

    ! Fill in the ghostcells, calculate energy
    call exchange_messages 

  end subroutine IH_match_ibc

  !============================================================================

  subroutine IH_get_for_global_buffer(&
       nR, nPhi,nTheta, BufferMinMax_DI, &
       TimeCoupling, iCompSource, iCompTarget, Buffer_VG)

    ! DESCRIPTION

    ! This subroutines fills a buffer grid by interpolating from a source 
    ! IH_BATSRUS grid using second-order trilinear interpolation.

    ! The buffer grid can be a spherical shell, or a segment of such a shell.

    ! All state variables in the source grid are interpolated, but only those
    ! needed for coupling (as determined by CON_coupler) are actually passed. 

    ! The filled buffer state vector is converted to SI units and vector 
    ! quantities are rotated to the target component coordinate system.

    ! INPUT:

    ! nR, nPhi, nTheta: grid spacing for the buffer grid
    ! BufferMinMAx_DI : Buffer grid minimum and maximum coordinates, in all
    ! dimensions.

    ! OUTPUT:

    ! Buffer_VG : defined for all coupling variables and all buffer grid points
    ! (including buffer ghost cells).

    ! REVISION HISTORY
    ! 30Dec2011 R. Oran   - initial version

    !USES:
    use IH_ModProcMH,         ONLY: iProc
    use IH_ModSize,           ONLY: nI, nJ, nK
    use IH_ModMain,           ONLY: UseB0, BuffR_, BuffPhi_, BuffTheta_
    use IH_ModAdvance,        ONLY: &
         State_VGB, UseElectronPressure, UseAnisoPressure
    use IH_ModB0, ONLY: B0_DGB
    use IH_ModPhysics,        ONLY: &
         No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, &
         UnitX_, UnitEnergyDens_,  inv_g
    use IH_ModVarIndexes,     ONLY: &
         Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, Pe_, &
         Ppar_, WaveFirst_, WaveLast_, Ehot_, nVar
    use CON_coupler,       ONLY: &
         iVar_V, DoCoupleVar_V, RhoCouple_, RhoUxCouple_,&
         RhoUzCouple_, PCouple_, BxCouple_, BzCouple_,  &
         PeCouple_, PparCouple_, WaveFirstCouple_,  &
         WaveLastCouple_, Bfield_, Wave_, EhotCouple_, &
         AnisoPressure_, ElectronPressure_, MultiFluid_,&
         MultiSpecie_, CollisionlessHeatFlux_, nVarCouple, Grid_C
    use ModCoordTransform, ONLY: sph_to_xyz
    use CON_axes,          ONLY: transform_matrix, transform_velocity
    use ModInterpolate,    ONLY: trilinear
    use IH_BATL_lib,          ONLY: xyz_to_coord, CoordMin_DB, CellSize_DB, nDim

    !INPUT ARGUMENTS:
    ! Buffer size and limits
    integer,intent(in) :: nR, nPhi, nTheta
    real, intent(in)   :: TimeCoupling
    real, intent(in)   :: BufferMinMax_DI(nDim,2)
    integer,intent(in) :: iCompSource, iCompTarget

    ! OUTPUT ARGUMENTS
    ! State variables to be fiiled in all buffer grid points
    real,dimension(nVarCouple,nR,0:nPhi+1,0:nTheta+1),intent(out):: Buffer_VG

    ! variables for defining the buffer grid

    integer :: nCell_D(3)
    real    :: SphMin_D(3), SphMax_D(3), dSph_D(3)

    ! Variables for interpolating from a IH_BATSRUS block to a buffer grid point

    ! Store complete interpolated state vector
    real :: StateInPoint_V(nVar)

    ! Store interpolated state variables needed for coupling
    real :: Buffer_V(nVarCouple), B0_D(3)

    ! Buffer grid cell center coordinates
    real :: CoordBuffer_D(3), XyzBuffer_D(3)

    ! Buffer grid cell center position  normalized by grid spacing
    ! (in IH_BATSRUS grid generalized coordinates)
    real :: BufferNorm_D(3)

    real :: SourceToTarget_DD(3,3)
    ! variable indices in buffer
    integer   :: &
         iRhoCouple,       &
         iRhoUxCouple,     &
         iRhoUzCouple,     &   
         iPCouple,         &       
         iPeCouple,        &      
         iPparCouple,      &    
         iBxCouple,        &      
         iBzCouple,        &      
         iWaveFirstCouple, &
         iWaveLastCouple,  &
         iEhotCouple


    integer   :: iPhiNew, iBlock, iPe, iR, iPhi, iTheta, iDim
    integer   :: iFound, jFound, kFound
    real      :: x, y, z, r, theta, phi

    ! Variables for testing 
    integer :: i, j ,k
    real    :: State_VG(nVar,-1:nI+2,-1:nJ+2,-1:nK+2), Btot_D(3)
    logical :: DoTest, DoTestMe
    character (len=*), parameter :: StringTest='IH_impose_par_flow_buffer'

    character (len=*), parameter :: NameSub='IH_get_for_buffer_grid'
    !--------------------------------------------------------------------------
    call CON_set_do_test(StringTest,DoTest,DoTestMe)
    if (DoTest .and. iProc == 0) write(*,*) &
         NameSub, ' in test mode: Imposing parallel flow inside buffer grid.'

    Buffer_VG = 0.0

    ! get variable indices in buffer
    iRhoCouple       = iVar_V(RhoCouple_)
    iRhoUxCouple     = iVar_V(RhoUxCouple_)
    iRhoUzCouple     = iVar_V(RhoUzCouple_)
    iPCouple         = iVar_V(PCouple_)
    iPeCouple        = iVar_V(PeCouple_)
    iPparCouple      = iVar_V(PparCouple_)
    iBxCouple        = iVar_V(BxCouple_)
    iBzCouple        = iVar_V(BzCouple_)
    iWaveFirstCouple = iVar_V(WaveFirstCouple_)
    iWaveLastCouple  = iVar_V(WaveLastCouple_)
    iEhotCouple      = iVar_V(EhotCouple_)

    ! Calculate buffer grid spacing
    nCell_D  = (/nR, nPhi, nTheta/)
    SphMin_D = BufferMinMax_DI(:,1)
    SphMax_D = BufferMinMax_DI(:,2)

    dSph_D     = (SphMax_D - SphMin_D)/real(nCell_D)
    dSph_D(BuffR_) = (SphMax_D(BuffR_) - SphMin_D(BuffR_))/(nCell_D(BuffR_)-1)

    ! Loop over buffer grid points
    do iR = 1, nR ; do iPhi = 1, nPhi ; do iTheta = 1, nTheta

       ! Find the coordinates of the current buffer grid point, 
       r     =  SphMin_D(BuffR_)     + (iR - 1)*dSph_D(BuffR_)
       Phi   =  SphMin_D(BuffPhi_)   + (real(iPhi)-0.5)*dSph_D(BuffPhi_)
       Theta =  SphMin_D(BuffTheta_) + (real(iTheta)-0.5)*dSph_D(BuffTheta_)

       ! Convert to xyz
       call sph_to_xyz(r, Theta, Phi, x, y, z)

       ! Find the block and PE in the IH_BATSRUS grid
       call IH_xyz_to_peblk(x,y,z,iPe,iBlock, .FALSE., iFound, jFound, kFound)

       ! Check if this block belongs to this processor
       if (iProc /= iPe) CYCLE

       ! Convert buffer grid point coordinate to IH_BATSRUS generalized coords
       XyzBuffer_D = (/x,y,z/)
       call xyz_to_coord(XyzBuffer_D, CoordBuffer_D)

       ! Buffer grid point position normalized by the grid spacing
       BufferNorm_D = (CoordBuffer_D - CoordMin_DB(:,iBlock)) &
            / CellSize_DB(:,iBlock) + 0.5

       if(DoTest) then
          ! Impose U||B prior to interpolation
          do i=-1,nI+2 ; do j=-1,nJ+2 ; do k=-1, nK+2
             State_VG(:,i,j,k) = State_VGB(:,i,j,k,iBlock)
             Btot_D = State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock)

             State_VG(RhoUx_:RhoUz_,i,j,k) = Btot_D*                 &
                  sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)*Btot_D)/ & 
                  (sum(Btot_D**2)+1e-40)

          end do; end do; end do

          ! Interpolate from the modified state in the block 
          ! to the buffer grid point
          StateInPoint_V = &
               trilinear(State_VG, nVar, -1, nI+2, -1, nJ+2, -1, nK+2, &
               BufferNorm_D) !, DoExtrapolate = .TRUE.)

       else

          ! Interpolate from the true solution block to the buffer grid point
          StateInPoint_V = &
               trilinear(State_VGB(:,:,:,:,iBlock), &
               nVar, -1, nI+2, -1, nJ+2, -1, nK+2, &
               BufferNorm_D) !, DoExtrapolate = .TRUE.)
       end if

       ! Fill in the coupled state variables

       Buffer_V(iRhoCouple)= StateInPoint_V(rho_)
       Buffer_V(iRhoUxCouple:iRhoUzCouple) = &
            StateInPoint_V(rhoUx_:rhoUz_)

       if(DoCoupleVar_V(Bfield_)) then
          if(UseB0)then
             B0_D = &
                  trilinear(B0_DGB(:,:,:,:,iBlock), &
                  3, -1, nI+2, -1, nJ+2, -1, nK+2, &
                  BufferNorm_D, DoExtrapolate = .TRUE.)
             Buffer_V(iBxCouple:iBzCouple) = &
                  StateInPoint_V(Bx_:Bz_) + B0_D
          else
             Buffer_V(iBxCouple:iBzCouple) = &
                  StateInPoint_V(Bx_:Bz_)
          end if
       end if

       if(DoCoupleVar_V(Wave_)) &
            Buffer_V(iWaveFirstCouple:iWaveLastCouple) = &
            StateInPoint_V(WaveFirst_:WaveLast_)

       Buffer_V(iPCouple)  = StateInPoint_V(p_) 

       if(DoCoupleVar_V(ElectronPressure_))then
          Buffer_V(iPeCouple) = StateInPoint_V(Pe_)
       else if(UseElectronPressure)then
          Buffer_V(iPCouple) = Buffer_V(iPCouple) + StateInPoint_V(Pe_)
       end if

       if(DoCoupleVar_V(AnisoPressure_)) Buffer_V(iPparCouple) = &
            StateInPoint_V(Ppar_)

       if(DoCoupleVar_V(CollisionlessHeatFlux_)) Buffer_V(iEhotCouple) = &
            StateInPoint_V(Ehot_)

       ! Convert to SI units
       Buffer_V(iRhoCouple) = &
            Buffer_V(iRhoCouple) * No2Si_V(UnitRho_)
       Buffer_V(iRhoUxCouple:iRhoUzCouple)= &
            Buffer_V(iRhoUxCouple:iRhoUzCouple) *No2Si_V(UnitRhoU_)
       Buffer_V(iPCouple) = Buffer_V(iPCouple) * No2Si_V(UnitP_)

       if(DoCoupleVar_V(Bfield_)) Buffer_V(iBxCouple:iBzCouple) = &
            Buffer_V(iBxCouple:iBzCouple)*No2Si_V(UnitB_)

       if(DoCoupleVar_V(Wave_)) &
            Buffer_V(iWaveFirstCouple:iWaveLastCouple) = &
            Buffer_V(iWaveFirstCouple:iWaveLastCouple) &
            * No2Si_V(UnitEnergyDens_)

       if(DoCoupleVar_V(ElectronPressure_)) Buffer_V(iPeCouple) = &
            Buffer_V(iPeCouple)*No2Si_V(UnitP_)

       if(DoCoupleVar_V(AnisoPressure_))Buffer_V(iPparCouple) = &
            Buffer_V(iPparCouple)*No2Si_V(UnitP_)

       if(DoCoupleVar_V(CollisionlessHeatFlux_))Buffer_V(iEhotCouple) = &
            Buffer_V(iEhotCouple)*No2Si_V(UnitEnergyDens_)

       ! ------------------------------------------------------
       !! Perform vector transformations if necessary
       !!        The followinf can be usefull if the source in an inertial
       !!        frame and the target is in a rotating frame.
       !!        WARNING: If you uncomment these lines make sure to disable
       !!                 any transformations done when buffer is read by target (e.g. IH_ModBuffer)
       !! START:
       !
       !if (Grid_C(iCompSource)%TypeCoord /=  &
       !     Grid_C(iCompTarget)%TypeCoord) then
       !   !Transform velocity
       !   ! NOTE: This transformation is only valid for a single fluid
       !   Buffer_V(iRhoUxCouple:iRhoUzCouple)=Buffer_V(iRhoCouple)*&
       !        transform_velocity(TimeCoupling,&
       !        Buffer_V(iRhoUxCouple:iRhoUzCouple)/Buffer_V(iRhoCouple),&
       !        No2Si_V(UnitX_)*XyzBuffer_D,&
       !        Grid_C(iCompSource)%TypeCoord,&
       !        Grid_C(iCompTarget)%TypeCoord)
       !
       !   ! Transform magnetic field
       !   SourceToTarget_DD = transform_matrix(TimeCoupling, &
       !        Grid_C(iCompSource)%TypeCoord, Grid_C(iCompTarget)%TypeCoord)

       !   Buffer_V(iBxCouple:iBzCouple) = &
       !        matmul(SourceToTarget_DD,Buffer_V(iBxCouple:iBzCouple))
       !end if
       !! END vector transformation

       ! DONE - fill the buffer grid
       Buffer_VG(:,iR, iPhi,iTheta) = Buffer_V

    end do; end do; end do

    ! Fill buffer grid ghost cells
    do iPhi = 1, nPhi 
       iPhiNew = iPhi + nPhi/2
       if (iPhiNew > nPhi) iPhiNew = iPhiNew - nPhi
       Buffer_VG(:,:,iPhi, 0) = Buffer_VG(:,:,iPhiNew, 1)
       Buffer_VG(:,:,iPhi,nTheta+1) = Buffer_VG(:,:,iPhiNew, nTheta)
    end do
    Buffer_VG(:,:,0,:) = Buffer_VG(:,:,nPhi,:)
    Buffer_VG(:,:,nPhi+1,:) = Buffer_VG(:,:,1,:)
  end subroutine IH_get_for_global_buffer

  !============================================================================

  subroutine IH_get_for_mh(nPartial,iGetStart,Get,W,State_V,nVar)

    !USES:
    use IH_ModAdvance, ONLY: State_VGB, UseElectronPressure, UseAnisoPressure
    use IH_ModB0,      ONLY: B0_DGB
    use IH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, inv_g
    use IH_ModPhysics, ONLY: UnitEnergyDens_
    use IH_ModAdvance, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, WaveFirst_, &
         WaveLast_, Pe_, Ppar_, Ehot_
    use IH_ModMain,    ONLY: UseB0, nDim

    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use CON_coupler,ONLY: OH_, iVar_V, DoCoupleVar_V, &
         RhoCouple_, RhoUxCouple_, &
         RhoUzCouple_, PCouple_, BxCouple_, BzCouple_, PeCouple_, PparCouple_,&
         WaveFirstCouple_, WaveLastCouple_, Bfield_, Wave_, AnisoPressure_, &
         ElectronPressure_, MultiFluid_, MultiSpecie_, EhotCouple_, &
         CollisionlessHeatFlux_

    !INPUT ARGUMENTS:
    integer,intent(in)              ::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)   ::Get
    type(WeightPtrType),intent(in)  ::W
    real,dimension(nVar),intent(out)::State_V

    integer   :: iGet, i, j, k, iBlock
    real      :: Weight,X_D(nDim)
    integer   :: &
         iRhoCouple,       &
         iRhoUxCouple,     &
         iRhoUzCouple,     &   
         iPCouple,         &       
         iPeCouple,        &      
         iPparCouple,      &    
         iBxCouple,        &      
         iBzCouple,        &      
         iWaveFirstCouple, &
         iWaveLastCouple,  &
         iEhotCouple

    character (len=*), parameter :: NameSub='IH_get_for_mh'

    ! 'Safety' parameter, to keep the coupler toolkit unchanged
    ! Useless, to my mind.
    integer, parameter:: Get4_ = OH_
    !--------------------------------------------------------------------------
    ! get variable indices in buffer
    iRhoCouple       = iVar_V(RhoCouple_)
    iRhoUxCouple     = iVar_V(RhoUxCouple_)
    iRhoUzCouple     = iVar_V(RhoUzCouple_)
    iPCouple         = iVar_V(PCouple_)
    iPeCouple        = iVar_V(PeCouple_)
    iPparCouple      = iVar_V(PparCouple_)
    iBxCouple        = iVar_V(BxCouple_)
    iBzCouple        = iVar_V(BzCouple_)
    iWaveFirstCouple = iVar_V(WaveFirstCouple_)
    iWaveLastCouple  = iVar_V(WaveLastCouple_)
    iEhotCouple      = iVar_V(EhotCouple_)

    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)

    State_V(iRhoCouple)= State_VGB(rho_,i,j,k,iBlock)*Weight
    State_V(iRhoUxCouple:iRhoUzCouple) = &
         State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)*Weight

    if(DoCoupleVar_V(Bfield_)) then
       if(UseB0)then
          State_V(iBxCouple:iBzCouple) = &
               (State_VGB(Bx_:Bz_,i,j,k,iBlock)+ B0_DGB(:,i,j,k,iBlock))*Weight
       else
          State_V(iBxCouple:iBzCouple) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
       end if
    end if

    if(DoCoupleVar_V(Wave_)) &
         State_V(iWaveFirstCouple:iWaveLastCouple) = &
         State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)*Weight

    State_V(iPCouple)  = State_VGB(p_,i,j,k,iBlock) *Weight

    if(DoCoupleVar_V(ElectronPressure_))then
       State_V(iPeCouple) = &
            State_VGB(Pe_,i,j,k,iBlock)*Weight
    else if(UseElectronPressure)then
       State_V(iPCouple) = &
            State_V(iPCouple) + State_VGB(Pe_,i,j,k,iBlock)*Weight
    end if

    if(DoCoupleVar_V(AnisoPressure_)) State_V(iPparCouple) = &
         State_VGB(Ppar_,i,j,k,iBlock)*Weight

    if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(iEhotCouple) = &
         State_VGB(Ehot_,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)

       State_V(iRhoCouple) = &
            State_V(iRhoCouple) + &
            State_VGB(rho_,i,j,k,iBlock)*Weight 
       State_V(iRhoUxCouple:iRhoUzCouple) = &
            State_V(iRhoUxCouple:iRhoUzCouple) + &
            State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) *Weight
       if(DoCoupleVar_V(Bfield_)) then
          if(UseB0)then
             State_V(iBxCouple:iBzCouple) = &
                  State_V(iBxCouple:iBzCouple) + &
                  (State_VGB(Bx_:Bz_,i,j,k,iBlock) + &
                  B0_DGB(:,i,j,k,iBlock))*Weight
          else 
             State_V(iBxCouple:iBzCouple) = &
                  State_V(iBxCouple:iBzCouple) + &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
          end if
       end if

       if(DoCoupleVar_V(Wave_)) &
            State_V(iWaveFirstCouple:iWaveLastCouple) = &
            State_V(iWaveFirstCouple:iWaveLastCouple) + &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(AnisoPressure_)) State_V(iPparCouple) = &
            State_V(iPparCouple) + &
            State_VGB(Ppar_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(iEhotCouple) = &
            State_V(iEhotCouple) + &
            State_VGB(Ehot_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(ElectronPressure_))then
          State_V(iPeCouple) = State_V(iPeCouple) + &
               State_VGB(Pe_,i,j,k,iBlock)*Weight
          State_V(iPCouple) = State_V(iPCouple) + &
               State_VGB(p_,i,j,k,iBlock) *Weight

       else if(UseElectronPressure)then
          State_V(iPCouple) = State_V(iPCouple) &
               + (State_VGB(p_,i,j,k,iBlock) + &
               State_VGB(Pe_,i,j,k,iBlock))*Weight
       else
          State_V(iPCouple) = State_V(iPCouple) + &
               State_VGB(p_,i,j,k,iBlock) *Weight
       end if
    end do

    ! Convert to SI units
    State_V(iRhoCouple) = &
         State_V(iRhoCouple) * No2Si_V(UnitRho_)
    State_V(iRhoUxCouple:iRhoUzCouple)= &
         State_V(iRhoUxCouple:iRhoUzCouple) *No2Si_V(UnitRhoU_)
    State_V(iPCouple) = State_V(iPCouple) * No2Si_V(UnitP_)

    if(DoCoupleVar_V(Bfield_)) State_V(iBxCouple:iBzCouple) = &
         State_V(iBxCouple:iBzCouple)*No2Si_V(UnitB_)

    if(DoCoupleVar_V(Wave_)) &
         State_V(iWaveFirstCouple:iWaveLastCouple) = &
         State_V(iWaveFirstCouple:iWaveLastCouple) &
         * No2Si_V(UnitEnergyDens_)

    if(DoCoupleVar_V(ElectronPressure_)) State_V(iPeCouple) = &
         State_V(iPeCouple)*No2Si_V(UnitP_)

    if(DoCoupleVar_V(AnisoPressure_))State_V(iPparCouple) = &
         State_V(iPparCouple)*No2Si_V(UnitP_)

    if(DoCoupleVar_V(CollisionlessHeatFlux_))State_V(iEhotCouple) = &
         State_V(iEhotCouple)*No2Si_V(UnitEnergyDens_)

  end subroutine IH_get_for_mh
  
  !===========================================================================

  subroutine IH_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,W,State_V,nVar)

    !USES:
    use IH_ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use IH_ModB0,         ONLY: B0_DGB
    use IH_ModPhysics,    ONLY: &
         No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, UnitX_, inv_g, &
         UnitEnergyDens_
    use IH_ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, &
         WaveFirst_, WaveLast_, Pe_, Ppar_, Ehot_
    use IH_ModMain,       ONLY: UseB0
    use IH_BATL_lib,      ONLY: Xyz_DGB
    use CON_coupler,   ONLY: SC_, IH_, iVar_V, DoCoupleVar_V,&
         RhoCouple_, RhoUxCouple_, RhoUzCouple_, &
         PCouple_, BxCouple_, BzCouple_, PeCouple_, EhotCouple_, &
         PparCouple_, WaveFirstCouple_, WaveLastCouple_, &
         Bfield_, Wave_, AnisoPressure_, CollisionlessHeatFlux_, &
         ElectronPressure_, MultiFluid_, MultiSpecie_
    use CON_router, ONLY: IndexPtrType, WeightPtrType
    use CON_coupler,ONLY: SC_

    !INPUT ARGUMENTS:
    integer,intent(in)               :: nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)    :: Get
    type(WeightPtrType),intent(in)   :: W
    real,dimension(nVar),intent(out) ::State_V

    integer::iGet, i, j, k, iBlock
    real :: Weight

    character (len=*), parameter :: NameSub='IH_get_for_mh_with_xyz'

    !The meaning of state intdex in buffer and in model can be 
    !different.
    integer   :: &
         iRhoCouple,       &
         iRhoUxCouple,     &
         iRhoUzCouple,     &
         iPCouple,         &
         iPeCouple,        &
         iPparCouple,      &
         iBxCouple,        &
         iBzCouple,        &
         iWaveFirstCouple, &
         iWaveLastCouple,  &
         iEhotCouple,      &
         BuffX_,           &
         BuffZ_

    integer, parameter:: Get4_ = SC_
    !--------------------------------------------------------------------------
    ! get variable indices in buffer
    iRhoCouple       = iVar_V(RhoCouple_)
    iRhoUxCouple     = iVar_V(RhoUxCouple_)
    iRhoUzCouple     = iVar_V(RhoUzCouple_)
    iPCouple         = iVar_V(PCouple_)
    iPeCouple        = iVar_V(PeCouple_)
    iPparCouple      = iVar_V(PparCouple_)
    iBxCouple        = iVar_V(BxCouple_)
    iBzCouple        = iVar_V(BzCouple_)
    iWaveFirstCouple = iVar_V(WaveFirstCouple_)
    iWaveLastCouple  = iVar_V(WaveLastCouple_)
    iEhotCouple      = iVar_V(EhotCouple_)

    BuffX_ = nVar - 2
    BuffZ_ = nVar

    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)

    State_V(iRhoCouple)= State_VGB(rho_,i,j,k,iBlock)*Weight
    State_V(iRhoUxCouple:iRhoUzCouple) = &
         State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)*Weight

    if(DoCoupleVar_V(Bfield_)) then
       if(UseB0)then
          State_V(iBxCouple:iBzCouple) = &
               (State_VGB(Bx_:Bz_,i,j,k,iBlock)+ B0_DGB(:,i,j,k,iBlock))*Weight
       else
          State_V(iBxCouple:iBzCouple) = &
               State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
       end if
    end if

    if(DoCoupleVar_V(Wave_)) &
         State_V(iWaveFirstCouple:iWaveLastCouple) = &
         State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)*Weight

    State_V(iPCouple)  = State_VGB(p_,i,j,k,iBlock)*Weight

    if(DoCoupleVar_V(ElectronPressure_))then
       State_V(iPeCouple) = &
            State_VGB(Pe_,i,j,k,iBlock)*Weight
    else if(UseElectronPressure)then
       State_V(iPCouple) = &
            State_V(iPCouple) + State_VGB(Pe_,i,j,k,iBlock)*Weight
    end if

    if(DoCoupleVar_V(AnisoPressure_)) State_V(iPparCouple) = &
         State_VGB(Ppar_,i,j,k,iBlock)*Weight

    if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(iEhotCouple) = &
         State_VGB(Ehot_,i,j,k,iBlock)*Weight

    State_V(BuffX_:BuffZ_) = &
         Xyz_DGB(:,i,j,k,iBlock)*State_VGB(rho_,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)

       State_V(iRhoCouple) = &
            State_V(iRhoCouple) + &
            State_VGB(rho_,i,j,k,iBlock)*Weight
       State_V(iRhoUxCouple:iRhoUzCouple) = &
            State_V(iRhoUxCouple:iRhoUzCouple) + &
            State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) *Weight
       if(DoCoupleVar_V(Bfield_)) then
          if(UseB0)then
             State_V(iBxCouple:iBzCouple) = &
                  State_V(iBxCouple:iBzCouple) + &
                  (State_VGB(Bx_:Bz_,i,j,k,iBlock) + &
                  B0_DGB(:,i,j,k,iBlock))*Weight
          else
             State_V(iBxCouple:iBzCouple) = &
                  State_V(iBxCouple:iBzCouple) + &
                  State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
          end if
       end if

       if(DoCoupleVar_V(Wave_)) &
            State_V(iWaveFirstCouple:iWaveLastCouple) = &
            State_V(iWaveFirstCouple:iWaveLastCouple) + &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(AnisoPressure_)) State_V(iPparCouple) = &
            State_V(iPparCouple) + &
            State_VGB(Ppar_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(iEhotCouple) = &
            State_V(iEhotCouple) + &
            State_VGB(Ehot_,i,j,k,iBlock)*Weight

       if(DoCoupleVar_V(ElectronPressure_))then
          State_V(iPeCouple) = State_V(iPeCouple) + &
               State_VGB(Pe_,i,j,k,iBlock)*Weight
          State_V(iPCouple) = State_V(iPCouple) + &
               State_VGB(p_,i,j,k,iBlock) *Weight

       else if(UseElectronPressure)then
          State_V(iPCouple) = State_V(iPCouple) &
               + (State_VGB(p_,i,j,k,iBlock) + &
               State_VGB(Pe_,i,j,k,iBlock))*Weight
       else
          State_V(iPCouple) = State_V(iPCouple) + &
               State_VGB(p_,i,j,k,iBlock) *Weight
       end if

       State_V(BuffX_:BuffZ_) = State_V(BuffX_:BuffZ_) + &
            Xyz_DGB(:,i,j,k,iBlock)*State_VGB(rho_,i,j,k,iBlock)*Weight
    end do

    ! Convert to SI units
    State_V(iRhoCouple) = &
         State_V(iRhoCouple) * No2Si_V(UnitRho_)
    State_V(iRhoUxCouple:iRhoUzCouple)= &
         State_V(iRhoUxCouple:iRhoUzCouple) *No2Si_V(UnitRhoU_)
    State_V(iPCouple) = State_V(iPCouple) * No2Si_V(UnitP_)

    if(DoCoupleVar_V(Bfield_)) State_V(iBxCouple:iBzCouple) = &
         State_V(iBxCouple:iBzCouple)*No2Si_V(UnitB_)

    if(DoCoupleVar_V(Wave_)) &
         State_V(iWaveFirstCouple:iWaveLastCouple) = &
         State_V(iWaveFirstCouple:iWaveLastCouple) &
         * No2Si_V(UnitEnergyDens_)

    if(DoCoupleVar_V(ElectronPressure_)) State_V(iPeCouple) = &
         State_V(iPeCouple)*No2Si_V(UnitP_)

    if(DoCoupleVar_V(AnisoPressure_))State_V(iPparCouple) = &
         State_V(iPparCouple)*No2Si_V(UnitP_)

    if(DoCoupleVar_V(CollisionlessHeatFlux_))State_V(iEhotCouple) = &
         State_V(iEhotCouple)*No2Si_V(UnitEnergyDens_)


  end subroutine IH_get_for_mh_with_xyz

  !============================================================================

  subroutine IH_get_a_line_point(nPartial,iGetStart,Get,W,State_V,nVar)
    !USES:
    use IH_ModAdvance,ONLY: State_VGB, Bx_, Bz_
    use IH_ModB0,     ONLY: B0_DGB
    use IH_BATL_lib,  ONLY: CellSize_DB
    use IH_ModMain,   ONLY: UseB0
    use CON_router

    !INPUT ARGUMENTS:
    integer,             intent(in) ::nPartial, iGetStart, nVar
    type(IndexPtrType),  intent(in) ::Get
    type(WeightPtrType), intent(in) ::W
    real,                intent(out)::State_V(nVar)

    integer::iGet, i, j, k, iBlock
    real :: Weight
    !--------------------------------------------------------------------------

    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)
    if(UseB0)then
       State_V(1:3)= &
            Weight*(State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))
    else
       State_V(1:3)= &
            Weight*State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    State_V(4:6)= CellSize_DB(:,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       if(UseB0)then
          State_V(1:3) = State_V(1:3)+ &
               Weight*(State_VGB(Bx_:Bz_,i,j,k,iBlock) &
               + B0_DGB(:,i,j,k,iBlock))
       else
          State_V(1:3) = State_V(1:3) + Weight*State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end if
       State_V(4:6) = State_V(4:6) + Weight*CellSize_DB(:,iBlock)
    end do
  end subroutine IH_get_a_line_point

  !============================================================================
  !BOP
  !ROUTINE: IH_put_from_mh - transform and put the data got from MH
  !INTERFACE:
  subroutine IH_put_from_mh(nPartial,&
       iPutStart,&
       Put,& 
       Weight,&
       DoAdd,&
       StateSI_V,&
       nVar)
    !USES:
    use CON_router,    ONLY: IndexPtrType, WeightPtrType
    use CON_coupler,   ONLY: OH_, iVar_V, DoCoupleVar_V, &
         RhoCouple_, RhoUxCouple_, RhoUzCouple_, &
         PCouple_, BxCouple_, BzCouple_, PeCouple_, EhotCouple_, &
         PparCouple_, WaveFirstCouple_, WaveLastCouple_, &
         Bfield_, Wave_, ElectronPressure_, AnisoPressure_, &
         CollisionlessHeatFlux_      
    use IH_ModAdvance,    ONLY: State_VGB, UseElectronPressure, &
         UseAnisoPressure
    use IH_ModB0,         ONLY: B0_DGB
    use IH_ModPhysics,    ONLY: Si2No_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, &
         UnitEnergyDens_
    use IH_ModMain,       ONLY: UseB0
    use IH_ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_, &
         WaveFirst_, WaveLast_, Pe_, Ppar_, Ehot_

    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::Weight
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::StateSI_V

    !REVISION HISTORY:
    !18JUL03     I.Sokolov <igorsok@umich.edu> - intial prototype/code
    !23AUG03                                     prolog
    !03SEP03     G.Toth    <gtoth@umich.edu>   - simplified
    !05APR11     R. Oran   <oran@umich.edu>    - Use non-fixed coupling indices
    !                                          derived by the coupler according
    !                                           to actual variable names
    !                                           (see use CON_coupler).
    !                                           Handle anisotropic pressure.
    !                                 
    !EOP

    character (len=*), parameter :: NameSub='IH_put_from_mh'

    real,dimension(nVar)::State_V
    integer             ::iPut, i, j, k, iBlock
    integer   :: &
         iRhoCouple,       &
         iRhoUxCouple,     &
         iRhoUzCouple,     &
         iPCouple,         &
         iPeCouple,        &
         iPparCouple,      &
         iBxCouple,        &
         iBzCouple,        &
         iWaveFirstCouple, &
         iWaveLastCouple,  &
         iEhotCouple

    integer, parameter:: PutFrom_ = OH_
    !--------------------------------------------------------------------------
    ! get variable indices in buffer
    iRhoCouple       = iVar_V(RhoCouple_)
    iRhoUxCouple     = iVar_V(RhoUxCouple_)
    iRhoUzCouple     = iVar_V(RhoUzCouple_)
    iPCouple         = iVar_V(PCouple_)
    iPeCouple        = iVar_V(PeCouple_)
    iPparCouple      = iVar_V(PparCouple_)
    iBxCouple        = iVar_V(BxCouple_)
    iBzCouple        = iVar_V(BzCouple_)
    iWaveFirstCouple = iVar_V(WaveFirstCouple_)
    iWaveLastCouple  = iVar_V(WaveLastCouple_)
    iEhotCouple      = iVar_V(EhotCouple_)

    ! Convert state variable in buffer to nirmalized units.
    State_V(iRhoCouple) = StateSI_V(iRhoCouple) * Si2No_V(UnitRho_)

    State_V(iRhoUxCouple:iRhoUzCouple) = &
         StateSI_V(iRhoUxCouple:iRhoUzCouple) * Si2No_V(UnitRhoU_)

    State_V(iPCouple) = StateSI_V(iPCouple) * Si2No_V(UnitP_)

    if(DoCoupleVar_V(Bfield_)) State_V(iBxCouple:iBzCouple) = &
         StateSI_V(iBxCouple:iBzCouple)* Si2No_V(UnitB_)

    if(DoCoupleVar_V(Wave_)) &
         State_V(iWaveFirstCouple:iWaveLastCouple) = &
         StateSI_V(iWaveFirstCouple:iWaveLastCouple) &
         * Si2No_V(UnitEnergyDens_)

    if(DoCoupleVar_V(ElectronPressure_))State_V(iPeCouple) = &
         StateSI_V(iPeCouple)*Si2No_V(UnitP_)

    if(DoCoupleVar_V(AnisoPressure_)) State_V(iPparCouple) = &
         StateSI_V(iPparCouple)*Si2No_V(UnitP_)

    if(DoCoupleVar_V(CollisionlessHeatFlux_)) State_V(iEhotCouple) = &
         StateSI_V(iEhotCouple)*Si2No_V(UnitEnergyDens_)

    i      = Put%iCB_II(1,iPutStart)
    j      = Put%iCB_II(2,iPutStart)
    k      = Put%iCB_II(3,iPutStart)
    iBlock = Put%iCB_II(4,iPutStart)

    if(DoAdd)then
       State_VGB(rho_,i,j,k,iBlock) = &
            State_VGB(rho_,i,j,k,iBlock) + &
            State_V(iRhoCouple)

       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = &
            State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) + &
            State_V(iRhoUxCouple:iRhoUzCouple)

       if (DoCoupleVar_V(Bfield_)) State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
            State_VGB(Bx_:Bz_,i,j,k,iBlock) + &
            State_V(iBxCouple:iBzCouple)

       if(DoCoupleVar_V(Wave_)) &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = &
            State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) + &
            State_V(iWaveFirstCouple:iWaveLastCouple)

       if(DoCoupleVar_V(ElectronPressure_))then
          State_VGB(Pe_,i,j,k,iBlock) = &
               State_VGB(Pe_,i,j,k,iBlock) + State_V(iPeCouple)
          State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
               + State_V(iPCouple)
       else if(UseElectronPressure)then
          State_VGB(Pe_,i,j,k,iBlock) = State_VGB(Pe_,i,j,k,iBlock) &
               + 0.5*State_V(iPCouple)
          ! correct pressure state variable
          State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
               + 0.5*State_V(iPCouple)
       else
          State_VGB(p_,i,j,k,iBlock) = State_VGB(p_,i,j,k,iBlock) &
               +State_V(iPCouple)
       end if

       if(DoCoupleVar_V(AnisoPressure_))then
          State_VGB(Ppar_,i,j,k,iBlock) = &
               State_VGB(Ppar_,i,j,k,iBlock) + State_V(iPparCouple)
       else if(UseAnisoPressure)then
          State_VGB(Ppar_,i,j,k,iBlock) = State_VGB(Ppar_,i,j,k,iBlock) &
               + State_V(iPCouple)
       end if

       if(DoCoupleVar_V(CollisionlessHeatFlux_))then
          State_VGB(Ehot_,i,j,k,iBlock) = &
               State_VGB(Ehot_,i,j,k,iBlock) + State_V(iEhotCouple)
       endif

    else

       State_VGB(rho_,i,j,k,iBlock)= State_V(iRhoCouple)
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) = &
            State_V(iRhoUxCouple:iRhoUzCouple)

       if(DoCoupleVar_V(Bfield_)) then
          if(UseB0)then
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  State_V(iBxCouple:iBzCouple) - &
                  B0_DGB(:,i,j,k,iBlock)
          else
             State_VGB(Bx_:Bz_,i,j,k,iBlock) = &
                  State_V(iBxCouple:iBzCouple)
          end if
       end if

       if(DoCoupleVar_V(Wave_))State_VGB(WaveFirst_:WaveLast_,i,j,k,iBlock) = &
            State_V(iWaveFirstCouple:iWaveLastCouple)

       State_VGB(p_,i,j,k,iBlock) = State_V(iPCouple)
       if(DoCoupleVar_V(AnisoPressure_))then
          State_VGB(Ppar_,i,j,k,iBlock) = State_V(iPparCouple)
       else if(UseAnisoPressure)then
          State_VGB(Ppar_,i,j,k,iBlock) = State_V(iPCouple)
       end if

       if(DoCoupleVar_V(CollisionlessHeatFlux_))then
          State_VGB(Ehot_,i,j,k,iBlock) = State_V(iEhotCouple)
       endif

       if(DoCoupleVar_V(ElectronPressure_))then
          State_VGB(Pe_,i,j,k,iBlock) = State_V(iPeCouple)
       else if(UseElectronPressure)then
          State_VGB(Pe_,i,j,k,iBlock) = 0.5*State_V(iPCouple)
          State_VGB(p_,i,j,k,iBlock) = 0.5*State_V(iPCouple)
       end if
    end if

  end subroutine IH_put_from_mh

  !============================================================================

  subroutine IH_get_for_sp(nPartial,iGetStart,Get,W,State_V,nVar)
    !USES:
    use IH_ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_, P_
    use IH_ModB0,      ONLY: B0_DGB
    use IH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitU_, UnitB_, UnitX_
    use CON_router
    use IH_BATL_lib, ONLY: Xyz_DGB
    use IH_ModMain,  ONLY: UseB0

    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    integer::iGet, i, j, k, iBlock
    real :: Weight
    !The meaning of state intdex in buffer and in model can be 
    !different. Below are the conventions for buffer:
    integer,parameter::&
         BuffRho_ = 1, &
         BuffUx_  = 2, &
         BuffUz_  = 4, &
         BuffBx_  = 5, &
         BuffBy_  = 6, &
         BuffBz_  = 7, &
         BuffP_   = 8, &
         BuffX_   = 9, &
         BuffZ_   =11
    !--------------------------------------------------------------------------

    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)

    State_V(BuffRho_)= Weight*State_VGB(rho_,i,j,k,iBlock)
    State_V(BuffUx_:BuffUz_)= Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)/&
         State_VGB(rho_,i,j,k,iBlock)
    if(UseB0)then
       State_V(BuffBx_:BuffBz_)= Weight*(State_VGB(Bx_:Bz_,i,j,k,iBlock)+ &
            B0_DGB(:,i,j,k,iBlock))
    else
       State_V(BuffBx_:BuffBz_)= Weight*State_VGB(Bx_:Bz_,i,j,k,iBlock)
    end if
    State_V(BuffP_)= Weight*State_VGB(P_,i,j,k,iBlock)
    State_V(BuffX_:BuffZ_) = Xyz_DGB(:,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       State_V(1) = State_V(1) + &
            Weight*State_VGB(rho_,i,j,k,iBlock)
       State_V(BuffUx_:BuffUz_) =  State_V(BuffUx_:BuffUz_) + &
            Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)/&
            State_VGB(rho_,i,j,k,iBlock)
       if(UseB0)then
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) + &
               Weight*(State_VGB(Bx_:Bz_,i,j,k,iBlock)+ B0_DGB(:,i,j,k,iBlock))
       else
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) + &
               Weight*State_VGB(Bx_:Bz_,i,j,k,iBlock)
       end if
       State_V(BuffP_)= State_V(BuffP_) + &
            Weight*State_VGB(P_,i,j,k,iBlock)
       State_V(BuffX_:BuffZ_) = State_V(BuffX_:BuffZ_) + &
            Xyz_DGB(:,i,j,k,iBlock)*Weight     
    end do
    ! Convert momentum to velocity and convert everything to SI units
    State_V(BuffUx_:BuffUz_) = State_V(BuffUx_:BuffUz_)  *No2Si_V(UnitU_)
    State_V(1)               = State_V(1)                *No2Si_V(UnitRho_)
    State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_)  *No2Si_V(UnitB_)
    State_V(BuffP_)          = State_V(BuffP_)           *No2Si_V(UnitP_)
    State_V(BuffX_:BuffZ_)   = State_V(BuffX_:BuffZ_)    *No2Si_V(UnitX_) 

  end subroutine IH_get_for_sp

  !============================================================================

  subroutine IH_get_for_gm(&
       nPartial,iGetStart,Get,W,State_V,nVar,TimeCoupling)

    !USES:
    use IH_ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_,P_
    use IH_ModB0,      ONLY: B0_DGB
    use IH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_
    use IH_ModMain,    ONLY: UseRotatingFrame,UseB0
    use CON_router

    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V
    real,intent(in)::TimeCoupling

    integer::iGet, i, j, k, iBlock
    real :: Weight, Momentum_D(3),Density

    character (len=*), parameter :: NameSub='IH_get_for_gm'
    !The meaning of state intdex in buffer and in model can be 
    !different. Below are the conventions for buffer:
    integer,parameter::&
         BuffRho_  =1,&
         BuffRhoUx_=2,&
         BuffRhoUz_=4,&
         BuffBx_   =5,&
         BuffBz_   =7,&
         BuffP_    =8


    !----------------------------------------------------------

    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)

    Density= State_VGB(rho_,         i,j,k,iBlock)
    State_V(BuffRho_)          = Density*Weight

    Momentum_D=State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)
    if(UseRotatingFrame)call add_density_omega_cross_r

    State_V(BuffRhoUx_:BuffRhoUz_) = Momentum_D*Weight
    if(UseB0)then
       State_V(BuffBx_:BuffBz_) = &
            (State_VGB(Bx_:Bz_,i,j,k,iBlock) + B0_DGB(:,i,j,k,iBlock))*Weight
    else
       State_V(BuffBx_:BuffBz_) = &
            State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
    end if

    State_V(BuffP_)  = State_VGB(P_,i,j,k,iBlock)*Weight

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)

       Density = State_VGB(rho_,i,j,k,iBlock) 
       State_V(BuffRho_)=State_V(BuffRho_) + Density*Weight

       Momentum_D = State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)

       if(UseRotatingFrame)call add_density_omega_cross_r

       State_V(BuffRhoUx_:BuffRhoUz_) = State_V(BuffRhoUx_:BuffRhoUz_) &
            + Momentum_D*Weight
       if(UseB0)then
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) &
               + (State_VGB(Bx_:Bz_,i,j,k,iBlock) &
               + B0_DGB(:,i,j,k,iBlock))*Weight
       else
          State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_) &
               + State_VGB(Bx_:Bz_,i,j,k,iBlock)*Weight
       end if
       State_V(BuffP_) = State_V(BuffP_) &
            + State_VGB(P_,i,j,k,iBlock)*Weight
    end do

    ! Convert to SI units
    State_V(BuffRho_)             = State_V(BuffRho_)       *No2Si_V(UnitRho_)
    State_V(BuffRhoUx_:BuffRhoUz_)= &
         State_V(BuffRhoUx_:BuffRhoUz_)                     *No2Si_V(UnitRhoU_)
    State_V(BuffBx_:BuffBz_)      = State_V(BuffBx_:BuffBz_)*No2Si_V(UnitB_)
    State_V(BuffP_)               = State_V(BuffP_)         *No2Si_V(UnitP_)

  contains
    !==========================================================================
    subroutine add_density_omega_cross_r
      ! Add Omega x R term. For IH Omega_D = (0,0,OmegaBody)
      use IH_BATL_lib,    ONLY: Xyz_DGB, x_, y_
      use IH_ModPhysics,  ONLY: OmegaBody
      !------------------------------------------------------------------------
      Momentum_D(x_) = Momentum_D(x_) &
           - Density*OmegaBody*Xyz_DGB(y_,i,j,k,iBlock)
      Momentum_D(y_)= Momentum_D(y_) &
           + Density*OmegaBody*Xyz_DGB(x_,i,j,k,iBlock)
    end subroutine add_density_omega_cross_r

  end subroutine IH_get_for_gm

end module IH_wrapper
