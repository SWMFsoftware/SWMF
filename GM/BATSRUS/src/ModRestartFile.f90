!^CFG COPYRIGHT UM

module ModRestartFile

  use ModProcMH,     ONLY: iProc
  use ModIO,         ONLY: Unit_Tmp, nFile, Dt_Output, Dn_Output, Restart_, &
       restart, save_restart_file, write_prefix, iUnitOut
  use ModMain,       ONLY: GlobalBlk, Global_Block_Number, nI, nJ, nK, Gcn, &
       nBlock, UnusedBlk, ProcTest, BlkTest, iTest, jTest, kTest, &
       n_step, Time_Simulation, dt_BLK, Cfl, CodeVersion, nByteReal, &
       NameThisComp
  use ModVarIndexes, ONLY: nVar, DefaultState_V
  use ModAdvance,    ONLY: State_VGB
  use ModGeometry,   ONLY: dx_BLK, dy_BLK, dz_BLK, xyzStart_BLK
  use ModParallel,   ONLY: iBlockRestartALL_A
  use ModIO,         ONLY: Restart_Bface                    !^CFG IF CONSTRAINB
  use ModCT,         ONLY: BxFace_BLK,ByFace_BLK,BzFace_BLK !^CFG IF CONSTRAINB
  use ModMain,       ONLY: UseConstrainB                    !^CFG IF CONSTRAINB
  use ModImplicit,   ONLY: n_prev, w_prev, dt_prev          !^CFG IF IMPLICIT
  use ModKind,       ONLY: Real4_, Real8_

  implicit none

  private ! except

  public read_restart_parameters
  public write_restart_files 
  public read_restart_files
  public read_octree_file
  public init_mod_restart_file

  ! Directories for input and output restart files
  character(len=100), public :: NameRestartInDir ="GM/restartIN/"
  character(len=100), public :: NameRestartOutDir="GM/restartOUT/"

  ! Local variables
  character(len=*), parameter :: StringRestartExt=".rst"
  character(len=*), parameter :: NameHeaderFile  ="restart.H"
  character(len=*), parameter :: NameDataFile    ="data.rst"
  character(len=*), parameter :: NameBlkFile     ="blk"

  logical :: RestartBlockLevels=.false. ! Load LEVmin,LEVmax in octree restart
  integer :: nByteRealRead = 8     ! Real precision in restart files

  ! One can use 'block' or 'one' format for input and output restart files
  ! The input format is set to 'block' for backwards compatibility
  character (len=10)  :: TypeRestartInFile ='block'

  ! The output format is platform dependent and it is reset in restart_init
  character (len=10)  :: TypeRestartOutFile='one'  

  character(len=100) :: NameFile

  ! Temporaray variables to read arbitrary precision data files
  real (Real8_) :: Dt8, Time8, Dxyz8_D(3), Xyz8_D(3)
  real (Real4_) :: Dt4, Time4, Dxyz4_D(3), Xyz4_D(3)
  real (Real8_) :: State8_CV(nI,nJ,nK,nVar), State8_VC(nVar,nI,nJ,nK)
  real (Real4_) :: State4_CV(nI,nJ,nK,nVar), State4_VC(nVar,nI,nJ,nK)
  !^CFG IF CONSTRAINB BEGIN
  real (Real8_) :: B8_X(nI+1,nJ,nK), B8_Y(nI,nJ+1,nK), B8_Z(nI,nJ,nK+1)
  real (Real4_) :: B4_X(nI+1,nJ,nK), B4_Y(nI,nJ+1,nK), B4_Z(nI,nJ,nK+1)
  !^CFG END CONSTRAINB

contains

  subroutine init_mod_restart_file
    use ModMPI, ONLY: MPI_HEADER_FILE

    if(MPI_HEADER_FILE == 'mpif90_Linux_openmpi.h')then
       ! There are some issues when the 'one' restart file is used 
       ! on grendel and nyx with openmpi (although it may have nothing
       ! to do with openmpi itself.
       TypeRestartOutFile = 'block'
    else
       ! On other machines the format 'one' works correctly.
       TypeRestartOutFile = 'one'
    end if
    if(iProc==0)then
       call write_prefix;
       write(iUnitOut,*) &
            'init_mod_restart_file: setting TypeRestartOutFile = ',&
            trim(TypeRestartOutFile)
    end if

    NameRestartInDir(1:2)  = NameThisComp
    NameRestartOutDir(1:2) = NameThisComp

  end subroutine init_mod_restart_file

  !============================================================================

  subroutine read_restart_parameters(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: fix_dir_name, check_dir
    use ModMain,      ONLY: UseStrict

    character(len=*), intent(in) :: NameCommand
    character(len=*), parameter:: NameSub = 'read_restart_parameters'
    !--------------------------------------------------------------------------

    select case(NameCommand)
    case("#SAVERESTART")
       call read_var('DoSaveRestart',save_restart_file)
       if(save_restart_file)then
          if(iProc==0)call check_dir(NameRestartOutDir)
          call read_var('DnSaveRestart',dn_output(restart_))
          call read_var('DtSaveRestart',dt_output(restart_))
          nfile=max(nfile,restart_)
       end if
    case("#NEWRESTART")
       restart=.true.
       call read_var('DoRestartBFace',restart_Bface) !^CFG IF CONSTRAINB
    case("#BLOCKLEVELSRELOADED")
       ! Sets logical for upgrade of restart files 
       ! to include LEVmin and LEVmax
       RestartBlockLevels=.true.
    case("#PRECISION")
       call read_var('nByteReal',nByteRealRead)
       if(nByteReal/=nByteRealRead)then
          if(iProc==0) write(*,'(a,i1,a,i1)') NameSub// &
               ' WARNING: BATSRUS was compiled with ',nByteReal,&
               ' byte reals, requested precision is ',nByteRealRead
          if(UseStrict)call stop_mpi(NameSub// &
               ' ERROR: differing precisions for reals')
       end if
    case("#RESTARTINDIR")
       call read_var("NameRestartInDir",NameRestartInDir)
       call fix_dir_name(NameRestartInDir)
       if (iProc==0) call check_dir(NameRestartInDir)
    case("#RESTARTOUTDIR")
       call read_var("NameRestartOutDir",NameRestartOutDir)
       call fix_dir_name(NameRestartOutDir)
       if (iProc==0) call check_dir(NameRestartOutDir)
    case("#RESTARTINFILE")
       call read_var('TypeRestartInFile',TypeRestartInFile)
    case("#RESTARTOUTFILE")
       call read_var('TypeRestartOutFile',TypeRestartOutFile)
    case default
       call stop_mpi(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine read_restart_parameters

  !============================================================================

  subroutine write_restart_files

    integer :: iBlock
    character(len=*), parameter :: NameSub='write_restart_files'
    !------------------------------------------------------------------------
    call timing_start(NameSub)

    call write_octree_file
    if(iProc==0)call write_restart_header
    select case(TypeRestartOutFile)
    case('block')
       do iBlock = 1,nBlock
          if (.not.unusedBLK(iBlock)) call write_restart_file(iBlock)
       end do
    case('direct','one')
       call write_one_restart_file
    case default
       call stop_mpi('Unknown TypeRestartOutFile='//TypeRestartOutFile)
    end select
    if(iProc==0)call save_advected_points

    call timing_stop(NameSub)

  end subroutine write_restart_files

  !===========================================================================

  subroutine read_restart_files

    integer :: iBlock
    character(len=*), parameter :: NameSub='read_restart_files'
    !------------------------------------------------------------------------
    call timing_start(NameSub)
    select case(TypeRestartInFile)
    case('block')
       do iBlock = 1, nBlock
          if (.not.unusedBLK(iBlock)) call read_restart_file(iBlock)
       end do
    case('direct','one')
       call read_one_restart_file
    case default
       call stop_mpi('Unknown TypeRestartInFile='//TypeRestartinFile)
    end select

    do iBlock = 1, nBlock
       if (.not.unusedBLK(iBlock)) call fix_block_geometry(iBlock)
    end do
    call set_body_flag
    call timing_stop(NameSub)

  end subroutine read_restart_files

  !===========================================================================

  subroutine write_restart_header

    use ModMain,       ONLY: Dt, NameThisComp, TypeCoordSystem,&
         nBlockAll, Body1, Time_Accurate, iStartTime_I, IsStandAlone
    use ModMain,       ONLY: UseBody2                     !^CFG IF SECONDBODY
    use ModVarIndexes, ONLY: NameEquation, nVar, nFluid
    use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2
    use ModGeometry, ONLY: XyzMin_D, XyzMax_D, &             
         TypeGeometry, UseCovariant, UseVertexBasedGrid      
    use ModParallel, ONLY: proc_dims
    use ModUser,     ONLY: NameUserModule, VersionUserModule
    use ModPhysics
    use CON_planet,  ONLY: NamePlanet
    use ModReadParam,ONLY: i_line_command
    use ModIO,       ONLY: NameMaxTimeUnit

    implicit none

    integer :: iFluid
    !--------------------------------------------------------------------------

    if (iProc/=0) RETURN

    open(unit_tmp,file=trim(NameRestartOutDir)//NameHeaderFile)

    write(unit_tmp,'(a)')'#CODEVERSION'
    write(unit_tmp,'(f5.2,a35)')CodeVersion,'CodeVersion'
    write(unit_tmp,*)

    write(unit_tmp,'(a)')'#USERMODULE'
    write(unit_tmp,'(a)')       NameUserModule
    write(unit_tmp,'(f5.2,a35)')VersionUserModule,'VersionUserModule'

    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#COMPONENT'
    write(unit_tmp,'(a2,a38)')NameThisComp,'NameComp'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#PRECISION'
    write(unit_tmp,'(i1,a39)')nByteReal,'nByteReal'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#EQUATION'
    write(unit_tmp,'(a,a32)')NameEquation,'NameEquation'
    write(unit_tmp,'(i8,a32)')nVar,'nVar'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#CHECKGRIDSIZE'
    write(unit_tmp,'(i8,a32)') nI,'nI'
    write(unit_tmp,'(i8,a32)') nJ,'nJ'
    write(unit_tmp,'(i8,a32)') nK,'nK'
    write(unit_tmp,'(i8,a32)') nBlockALL,'MinBlockALL'
    if (IsStandAlone .and. NameThisComp == 'GM') then
       write(unit_tmp,*)
       write(unit_tmp,'(a)')'#PLANET'
       write(unit_tmp,'(a,a32)') NamePlanet,'NamePlanet'
       if(i_line_command("#IDEALAXES", iSessionIn=1) > 0)then
          write(unit_tmp,*)
          write(unit_tmp,'(a)')'#IDEALAXES'
       end if
    end if
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#NEWRESTART'
    write(unit_tmp,'(l1,a39)')UseConstrainB,'DoRestartBFace'!^CFG IF CONSTRAINB
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#RESTARTINFILE'
    ! Note that the output file format is saved as the input for next restart
    write(unit_tmp,'(a,a30)')TypeRestartOutFile,'TypeRestartInFile'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#BLOCKLEVELSRELOADED'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#NSTEP'
    write(unit_tmp,'(i8,a32)')n_step,'nStep'
    write(unit_tmp,*)
    if(n_prev == n_step)then                            !^CFG IF IMPLICIT BEGIN
       write(unit_tmp,'(a)')'#NPREVIOUS'
       write(unit_tmp,'(i8,a32)')n_prev,'nPrev'
       write(unit_tmp,'(1pe20.12,a20)')dt_prev,'DtPrev'
       write(unit_tmp,*)
    end if                                              !^CFG END IMPLICIT
    write(unit_tmp,'(a)')'#STARTTIME'
    write(unit_tmp,'(i8,a32)')iStartTime_I(1),'iYear'
    write(unit_tmp,'(i8,a32)')iStartTime_I(2),'iMonth'
    write(unit_tmp,'(i8,a32)')iStartTime_I(3),'iDay'
    write(unit_tmp,'(i8,a32)')iStartTime_I(4),'iHour'
    write(unit_tmp,'(i8,a32)')iStartTime_I(5),'iMinute'
    write(unit_tmp,'(i8,a32)')iStartTime_I(6),'iSecond'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#TIMESIMULATION'
    write(unit_tmp,'(es15.8,a25)')time_simulation,'tSimulation'
    write(unit_tmp,*)
    if(UseCovariant)then                        
       write(unit_tmp,'(a)')'#GRIDGEOMETRY'
       write(unit_tmp,'(a)')trim(TypeGeometry)
       write(unit_tmp,*)
       write(unit_tmp,'(a)')'#VERTEXBASEDGRID'
       write(unit_tmp,'(l1,a39)') UseVertexBasedGrid,'UseVertexBasedGrid'
       write(unit_tmp,*)
       write(unit_tmp,'(a)')'#LIMITGENCOORD1'                   
       write(unit_tmp,'(1pe13.5,a27)')XyzMin_D(1),'XyzMin_D(1)' 
       write(unit_tmp,'(1pe13.5,a27)')XyzMax_D(1),'XyzMax_D(1)' 
       write(unit_tmp,*)
    end if                                      
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#GRID'
    write(unit_tmp,'(i8,a32)')proc_dims(1),'nRootBlockX'
    write(unit_tmp,'(i8,a32)')proc_dims(2),'nRootBlockY'
    write(unit_tmp,'(i8,a32)')proc_dims(3),'nRootBlockZ'
    write(unit_tmp,'(1pe13.5,a27)')x1,'xMin'
    write(unit_tmp,'(1pe13.5,a27)')x2,'xMax'
    write(unit_tmp,'(1pe13.5,a27)')y1,'yMin'
    write(unit_tmp,'(1pe13.5,a27)')y2,'yMax'
    write(unit_tmp,'(1pe13.5,a27)')z1,'zMin'
    write(unit_tmp,'(1pe13.5,a27)')z2,'zMax'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#COORDSYSTEM'
    write(unit_tmp,'(a3,a37)') TypeCoordSystem,'TypeCoordSystem'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#SOLARWIND'
    write(unit_tmp,'(1pe15.7,a25)')SW_n_dim,  'SwNDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_T_dim,  'SwTDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_Ux_dim, 'SwUxDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_Uy_dim, 'SwUyDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_Uz_dim, 'SwUzDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_Bx_dim, 'SwBxDdim'
    write(unit_tmp,'(1pe15.7,a25)')SW_By_dim, 'SwByDim'
    write(unit_tmp,'(1pe15.7,a25)')SW_Bz_dim, 'SwBzDim'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#IOUNITS'
    write(unit_tmp,'(a20,a20)')TypeIoUnit,'TypeIoUnit'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#NORMALIZATION'
    write(unit_tmp,'(a)')'READ'
    write(unit_tmp,'(1pe23.15,a17)')No2Si_V(UnitX_),   'No2SiUnitX'
    write(unit_tmp,'(1pe23.15,a17)')No2Si_V(UnitU_),   'No2SiUnitU'
    write(unit_tmp,'(1pe23.15,a17)')No2Si_V(UnitRho_), 'No2SiUnitRho'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'#PLOTFILENAME'
    write(unit_tmp,'(a10,a30)') NameMaxTimeUnit, 'NameMaxTimeUnit'
    write(unit_tmp,*)

    if(body1)then
       write(unit_tmp,'(a)')'#BODY'
       write(unit_tmp,'(l1,a39)')      .true., 'UseBody'
       write(unit_tmp,'(1pe13.5,a27)') rBody, 'rBody'
       if(NameThisComp=='GM') &
            write(unit_tmp,'(1pe13.5,a27)') rCurrents, 'rCurrents'
       do iFluid = IonFirst_, nFluid
          write(unit_tmp,'(1pe13.5,a27)') BodyNDim_I(iFluid), 'BodyNDim'
          write(unit_tmp,'(1pe13.5,a27)') BodyTDim_I(iFluid), 'BodyTDim'
       end do
       write(unit_tmp,*)
    end if
    !^CFG IF SECONDBODY BEGIN
    if(UseBody2)then
       write(unit_tmp,'(a)')'#SECONDBODY'
       write(unit_tmp,'(l1,a39)')     UseBody2,      'UseBody2'
       write(unit_tmp,'(1pe13.5,a27)')Rbody2,        'rBody2'
       write(unit_tmp,'(1pe13.5,a27)')xbody2,        'xBody2'
       write(unit_tmp,'(1pe13.5,a27)')ybody2,        'yBody2'
       write(unit_tmp,'(1pe13.5,a27)')zbody2,        'zBody2'
       write(unit_tmp,'(1pe13.5,a27)')rCurrentsBody2,'rCurrentsBody2'
       write(unit_tmp,'(1pe13.5,a27)')RhoDimBody2,   'RhoDimBody2'
       write(unit_tmp,'(1pe13.5,a27)')tDimBody2,     'tDimBody2'
       write(unit_tmp,*)
    end if
    !^CFG END SECONDBODY
    write(unit_tmp,'(a)')'#END'
    write(unit_tmp,*)
    write(unit_tmp,'(a)')'Additional info'
    write(unit_tmp,*)
    write(unit_tmp,'(l8,a)') time_accurate,   ' time_accurate'
    write(unit_tmp,*)
    if(time_accurate)write(unit_tmp,'(2(1pe13.5),a)')&
         time_simulation, dt,                 ' time_simulation, dt'

    write(unit_tmp,'(a)')'Io2Si_V='
    write(unit_tmp,'(100es13.5)') Io2Si_V
    write(unit_tmp,'(a)')'No2Io_V='
    write(unit_tmp,'(100es13.5)') No2Io_V

    close(unit_tmp)

  end subroutine write_restart_header

  !============================================================================

  subroutine read_restart_file(iBlock)

    use ModEnergy, ONLY: calc_energy_cell

    integer, intent(in) :: iBlock

    integer   :: iVar, i, j, k, iError, iBlockRestart
    character :: StringDigit
    real      :: tSimulationRead

    character (len=*), parameter :: NameSub='read_restart_file'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------
    if(iProc==PROCtest.and.iBlock==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    iBlockRestart = iBlockRestartALL_A(global_block_number(iBlock))
    write(StringDigit,'(i1)') max(5,1+int(alog10(real(iBlockRestart))))

    write(NameFile,'(a,i'//StringDigit//'.'//StringDigit//',a)') &
         trim(NameRestartInDir)//NameBlkFile,iBlockRestart,StringRestartExt

    open(unit_tmp, file=NameFile, status='old', form='UNFORMATTED',&
         iostat = iError)

    if(iError /= 0) call stop_mpi(NameSub// &
         ' read_restart_file could not open: '//trim(NameFile))

    ! Fill in ghost cells
    do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
       State_VGB(1:nVar, i, j, k, iBlock) = DefaultState_V(1:nVar)
    end do;end do;end do

    ! Do not overwrite time_simulation which is read from header file
    if(nByteRealRead == 8)then
       read(unit_tmp, iostat = iError) Dt8, Time8
       dt_BLK(iBlock) = Dt8
       tSimulationRead   = Time8

       read(unit_tmp, iostat = iError) Dxyz8_D, Xyz8_D
       Dx_BLK(iBlock) = Dxyz8_D(1)
       Dy_BLK(iBlock) = Dxyz8_D(2)
       Dz_BLK(iBlock) = Dxyz8_D(3)
       XyzStart_BLK(:,iBlock) = Xyz8_D

       read(Unit_tmp, iostat = iError) State8_CV
       do iVar = 1, nVar
          State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock) = State8_CV(:,:,:,iVar)
       end do

       if(Restart_Bface)then                          !^CFG IF CONSTRAINB BEGIN
          read(Unit_tmp, iostat = iError) b8_X, b8_Y, b8_Z               
          BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock) = b8_X
          ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock) = b8_Y
          BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock) = b8_Z
       end if                                         !^CFG END CONSTRAINB
       if(n_prev==n_step) then                        !^CFG IF IMPLICIT BEGIN
          read(Unit_tmp, iostat = iError) State8_CV
          w_prev(:,:,:,:,iBlock) = State8_CV
       end if                                         !^CFG END IMPLICIT
    else
       read(unit_tmp, iostat = iError) Dt4, Time4
       dt_BLK(iBlock) = Dt4
       tSimulationRead   = Time4

       read(unit_tmp, iostat = iError) Dxyz4_D, Xyz4_D
       Dx_BLK(iBlock) = Dxyz4_D(1)
       Dy_BLK(iBlock) = Dxyz4_D(2)
       Dz_BLK(iBlock) = Dxyz4_D(3)
       XyzStart_BLK(:,iBlock) = Xyz4_D

       read(Unit_tmp, iostat = iError) State4_CV
       do iVar = 1, nVar
          State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock) = State4_CV(:,:,:,iVar)
       end do

       if(Restart_Bface)then                          !^CFG IF CONSTRAINB BEGIN
          read(Unit_tmp, iostat = iError) b4_X, b4_Y, b4_Z               
          BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock) = b4_X
          ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock) = b4_Y
          BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock) = b4_Z
       end if                                         !^CFG END CONSTRAINB
       if(n_prev==n_step) then                        !^CFG IF IMPLICIT BEGIN
          read(Unit_tmp, iostat = iError) State4_CV
          w_prev(:,:,:,:,iBlock) = State4_CV
       end if                                         !^CFG END IMPLICIT
    endif

    if(iError /= 0) call stop_mpi(NameSub// &
         ' could not read data from '//trim(NameFile))

    close(unit_tmp)

    if(CodeVersion>5.60 .and. CodeVersion <7.00) &
         dt_BLK(iBlock)=dt_BLK(iBlock)/cfl

    call calc_energy_cell(iBlock)

    if(Dx_BLK(iBlock) < 0 .or. Dy_BLK(iBlock) < 0 .or. Dz_BLK(iBlock) < 0 &
         .or. Dt_BLK(iBlock) < 0 .or. tSimulationRead < 0)then
       write(*,*)NameSub,': corrupt restart data!!!'
       write(*,*)'iBlock  =', iBlock
       write(*,*)'Dxyz    =', Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock)
       write(*,*)'Dt,tSim =', Dt_BLK(iBlock), tSimulationRead
       write(*,*)'XyzStart=', XyzStart_BLK(:,iBlock)
       write(*,*)'State111=', State_VGB(1:nVar,1,1,1,iBlock)
       call stop_mpi(NameSub//': corrupt restart data!!!')
    end if

    if(DoTestMe)then
       write(*,*)NameSub,': iProc, iBlock =',iProc, iBlock
       write(*,*)NameSub,': dt,tSimRead =',dt_BLK(iBlock),tSimulationRead
       write(*,*)NameSub,': dx,dy,dz_BLK=',&
            dx_BLK(iBlock),dy_BLK(iBlock),dz_BLK(iBlock)
       write(*,*)NameSub,': xyzStart_BLK=',xyzStart_BLK(:,iBlock)
       write(*,*)NameSub,': State_VGB   =', &
            State_VGB(:,Itest,Jtest,Ktest,iBlock)
       write(*,*)NameSub,' finished'
    end if

  end subroutine read_restart_file

  !===========================================================================

  subroutine write_restart_file(iBlock)

    integer, intent(in) :: iBlock

    character (len=*), parameter :: NameSub='write_restart_file'
    integer:: iVar, iBlockRestart
    character:: StringDigit
    !--------------------------------------------------------------------

    iBlockRestart = global_block_number(iBlock)
    write(StringDigit,'(i1)') max(5,int(1+alog10(real(iBlockRestart))))

    write(NameFile,'(a,i'//StringDigit//'.'//StringDigit//',a)') &
         trim(NameRestartOutDir)//NameBlkFile,iBlockRestart,StringRestartExt

    open(unit_tmp, file=NameFile, status="replace", form='UNFORMATTED')

    write(Unit_tmp)  dt_BLK(iBlock),time_Simulation
    write(Unit_tmp) &
         dx_BLK(iBlock),dy_BLK(iBlock),dz_BLK(iBlock),&
         xyzStart_BLK(:,iBlock)
    write(Unit_tmp) &
         ( State_VGB(iVar,1:nI,1:nJ,1:nK,iBlock), iVar=1,nVar)
    if(UseConstrainB)then                            !^CFG iF CONSTRAINB BEGIN
       write(Unit_tmp) &
            BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock),&
            ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock),&
            BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock)
    end if                                           !^CFG END CONSTRAINB
    if(n_prev==n_step) write(Unit_tmp) &             !^CFG IF IMPLICIT
         w_prev(:,:,:,:,iBlock)                   !^CFG IF IMPLICIT
    close(unit_tmp)

  end subroutine write_restart_file

  !============================================================================

  subroutine open_one_restart_file(DoRead)

    logical, intent(in) :: DoRead

    integer :: lRecord, l, iError
    character(len=*), parameter :: NameSub='open_one_restart_file'
    logical :: DoTest, DoTestme
    !-------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)
    if(DoTestMe)write(*,*) NameSub,' starting with DoRead=',DoRead

    ! Calculate the record length for the first block
    inquire (IOLENGTH = lRecord ) &
         Dt_BLK(1), Dx_BLK(1), Dy_BLK(1), Dz_BLK(1),&
         XyzStart_BLK(:,1), &
         State_VGB(1:nVar,1:nI,1:nJ,1:nK,1)

    if(DoRead .and. Restart_Bface .or. &         !^CFG iF CONSTRAINB BEGIN
         .not.DoRead .and. UseConstrainB)then       
       inquire (IOLENGTH = l) &
            BxFace_BLK(1:nI+1,1:nJ,1:nK,1),&
            ByFace_BLK(1:nI,1:nJ+1,1:nK,1),&
            BzFace_BLK(1:nI,1:nJ,1:nK+1,1)
       lRecord = lRecord + l
    end if                                       !^CFG END CONSTRAINB
    if(n_prev==n_step)then                       !^CFG IF IMPLICIT BEGIN
       inquire (IOLENGTH = l) w_prev(:,:,:,:,1)
       lRecord = lRecord + l
    end if                                       !^CFG END IMPLICIT

    if(DoTestMe)write(*,*) NameSub,' nByteReal, nByteRealRead, lRecord=',&
          nByteReal, nByteRealRead, lRecord

    if(DoRead)then
       if(nByteReal /= nByteRealRead) &
            lRecord = (lRecord * nByteRealRead)/nByteReal

       NameFile = trim(NameRestartInDir)//NameDataFile

       open(Unit_Tmp, file=NameFile, &
            RECL = lRecord, ACCESS = 'direct', FORM = 'unformatted', &
            status = 'old', iostat=iError)
    else
       NameFile = trim(NameRestartOutDir)//NameDataFile
       ! Delete and open file from proc 0
       if(iProc==0) open(Unit_Tmp, file=NameFile, &
               RECL = lRecord, ACCESS = 'direct', FORM = 'unformatted', &
               status = 'replace', iostat=iError)

       ! Make sure that all processors wait until the file is re-opened
       call barrier_mpi

       if(iProc > 0)open(Unit_Tmp, file=NameFile, &
            RECL = lRecord, ACCESS = 'direct', FORM = 'unformatted', &
            status = 'old', iostat=iError)
    end if
    if(iError /= 0)then
       write(*,*) NameSub,': ERROR for DoRead=',DoRead
       call stop_mpi(NameSub//': could not open file='//NameFile)
    end if

  end subroutine open_one_restart_file

  !============================================================================

  subroutine read_one_restart_file

    character (len=*), parameter :: NameSub='read_one_restart_file'
    integer :: i, j, k, iBlock, iRec
    logical :: IsRead, DoTest, DoTestMe
    !-------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    call open_one_restart_file(DoRead = .true.)

    if(DoTestMe)write(*,*) NameSub,' starting with nBlock=', nBlock

    do iBlock = 1, nBlock

       if(UnusedBlk(iBlock)) CYCLE
       ! Use the global block index as the record number
       iRec = iBlockRestartALL_A(global_block_number(iBlock))

       if(DoTestMe) write(*,*) NameSub,' iBlock,iBlockGlobal, iRec=',&
            iBlock, global_block_number(iBlock), iRec

       ! Fill in ghost cells
       do k=1-gcn,nK+gcn; do j=1-gcn,nK+gcn; do i=1-gcn,nK+gcn
          State_VGB(1:nVar, i, j, k, iBlock) = DefaultState_V(1:nVar)
       end do; end do; end do

       IsRead = .false.
       if(nByteRealRead == 4)then
          if(Restart_Bface)then                       !^CFG IF CONSTRAINB BEGIN
             ! Read with face centered magnetic field for constrained transport
             read(Unit_Tmp, rec=iRec) Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  B4_X, B4_Y, B4_Z
             BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock) = B4_X
             ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock) = B4_Y
             BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock) = B4_Z
             IsRead = .true.
          endif                                       !^CFG END CONSTRAINB
          if(n_prev==n_step)then                      !^CFG IF IMPLICIT BEGIN
             ! Read with previous state for sake of implicit BDF2 scheme
             read(Unit_Tmp, rec=iRec) Dt4, Dxyz4_D, Xyz4_D, State4_VC, &
                  State4_CV
             w_prev(:,:,:,:,iBlock) = State4_CV
             IsRead = .true.
          end if                                       !^CFG END IMPLICIT
          if(.not.IsRead) &
               read(Unit_Tmp, rec=iRec) Dt4, Dxyz4_D, Xyz4_D, State4_VC

          Dt_BLK(iBlock) = Dt4
          Dx_BLK(iBlock) = Dxyz4_D(1)
          Dy_BLK(iBlock) = Dxyz4_D(2)
          Dz_BLK(iBlock) = Dxyz4_D(3)
          XyzStart_BLK(:,iBlock) = Xyz4_D
          State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock) = State4_VC

       else
          if(Restart_Bface)then                       !^CFG IF CONSTRAINB BEGIN
             ! Read with face centered magnetic field for constrained transport
             read(Unit_Tmp, rec=iRec) Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  B8_X, B8_Y, B8_Z
             BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock) = B8_X
             ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock) = B8_Y
             BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock) = B8_Z
             IsRead = .true.
          endif                                       !^CFG END CONSTRAINB
          if(n_prev==n_step)then                      !^CFG IF IMPLICIT BEGIN
             ! Read with previous state for sake of implicit BDF2 scheme
             read(Unit_Tmp, rec=iRec) Dt8, Dxyz8_D, Xyz8_D, State8_VC, &
                  State8_CV
             w_prev(:,:,:,:,iBlock) = State8_CV
             IsRead = .true.
          end if                                       !^CFG END IMPLICIT
          if(.not.IsRead) &
               read(Unit_Tmp, rec=iRec) Dt8, Dxyz8_D, Xyz8_D, State8_VC

          Dt_BLK(iBlock) = Dt8
          Dx_BLK(iBlock) = Dxyz8_D(1)
          Dy_BLK(iBlock) = Dxyz8_D(2)
          Dz_BLK(iBlock) = Dxyz8_D(3)
          XyzStart_BLK(:,iBlock) = Xyz8_D
          State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock) = State8_VC
       end if

       if(Dx_BLK(iBlock) < 0 .or. Dy_BLK(iBlock) < 0 .or. Dz_BLK(iBlock) < 0 &
            .or. Dt_BLK(iBlock) < 0)then
          write(*,*)NameSub,': corrupt restart data!!!'
          write(*,*)'iBlock  =',iBlock
          write(*,*)'Dxyz    =',Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock)
          write(*,*)'Dt      =', Dt_BLK(iBlock)
          write(*,*)'XyzStart=',XyzStart_BLK(:,iBlock)
          write(*,*)'State111=',State_VGB(1:nVar,1,1,1,iBlock)
          call stop_mpi(NameSub//': corrupt restart data!!!')
       end if
    end do

    close(Unit_Tmp)

  end subroutine read_one_restart_file

  !============================================================================

  subroutine write_one_restart_file

    character (len=*), parameter :: NameSub='write_one_restart_file'
    integer :: iBlock, iRec
    !--------------------------------------------------------------------

    call open_one_restart_file(DoRead = .false.)

    do iBlock = 1, nBlock

       if(UnusedBlk(iBlock)) CYCLE
       ! Use the global block index as the record number
       iRec = global_block_number(iBlock)

       if(UseConstrainB)then                          !^CFG IF CONSTRAINB BEGIN
          ! Save face centered magnetic field 
          write(Unit_tmp, rec=iRec)  Dt_BLK(iBlock),&
               Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock),&
               XyzStart_BLK(:,iBlock), &
               State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock), &
               BxFace_BLK(1:nI+1,1:nJ,1:nK,iBlock),&
               ByFace_BLK(1:nI,1:nJ+1,1:nK,iBlock),&
               BzFace_BLK(1:nI,1:nJ,1:nK+1,iBlock)
          CYCLE
       endif                                          !^CFG END CONSTRAINB
       if(n_prev==n_step)then                         !^CFG IF IMPLICIT BEGIN
          ! Save previous time step for sake of BDF2 scheme
          write(Unit_tmp, rec=iRec) &
               Dt_BLK(iBlock), &
               Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock), &
               XyzStart_BLK(:,iBlock), &
               State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock), &
               w_prev(:,:,:,:,iBlock)
          CYCLE
       endif                                          !^CFG END IMPLICIT

       write(Unit_tmp, rec=iRec) &
            Dt_BLK(iBlock), &
            Dx_BLK(iBlock), Dy_BLK(iBlock), Dz_BLK(iBlock), &
            XyzStart_BLK(:,iBlock), &
            State_VGB(1:nVar,1:nI,1:nJ,1:nK,iBlock)
    end do

    close(Unit_Tmp)

  end subroutine write_one_restart_file

  !============================================================================

  subroutine read_octree_file
    use ModProcMH
    use ModParallel, ONLY : nBLK,proc_dims
    use ModOctree
    use ModIoUnit, Only : UNITTMP_
    use ModIO, ONLY : iUnitOut, write_prefix

    implicit none

    integer :: i,j,k, total_number_of_blocks_needed, BlksPerPE, iError,nError
    integer, dimension(3) :: r_proc_dims
    character (len=4), Parameter :: octree_ext=".rst"
    logical :: isRoot
    type (adaptive_block_ptr) :: octree

    open(UNITTMP_, file=trim(NameRestartInDir)//"octree"//octree_ext, &
         status="old", form="UNFORMATTED")
    read(UNITTMP_) r_proc_dims(1),r_proc_dims(2),r_proc_dims(3)
    read(UNITTMP_) total_number_of_blocks_needed

    if ( (r_proc_dims(1) /= proc_dims(1)) .or. &
         (r_proc_dims(2) /= proc_dims(2)) .or. &
         (r_proc_dims(3) /= proc_dims(3)) ) then
       write(*,*) "read_octree_file: PE = ",iProc, &
            " Initial processor outlay incorrect, ", &
            & "proc_dims = ",proc_dims," r_proc_dims = ", &
            r_proc_dims
       call stop_mpi('ERROR in read_octree_file')
    end if

    if (total_number_of_blocks_needed > nProc*nBLK) then
       write(*,*) "read_octree_file: PE = ",iProc, &
            " Error, insufficient number of solution blocks, ", &
            "total_number_of_blocks_needed = ",total_number_of_blocks_needed
       call stop_mpi('ERROR in read_octree_file')
    end if

    BlksPerPE = ((total_number_of_blocks_needed-1)/nProc)+1
    if(iProc==0)then
       call write_prefix; write(iUnitOut,*) &
            'Reading restart files with ',total_number_of_blocks_needed, &
            ' blocks (',BlksPerPE,' per PE)'
    end if

    do k = 1, proc_dims(3)
       do j = 1, proc_dims(2)
          do i = 1, proc_dims(1)
             isRoot = .true.
             octree % ptr => octree_roots(i, j, k) % ptr
             call read_octree_soln_block(octree, BlksPerPE, isRoot)
          end do
       end do
    end do

    close(UNITTMP_)

  end subroutine read_octree_file
  !============================================================================

  recursive subroutine read_octree_soln_block(octree, BlksPerPE, isRoot)
    use ModProcMH
    use ModMain, ONLY : nBlockMax, unusedBLK
    use ModGeometry, ONLY : xyzStart, dxyz
    use ModIoUnit, ONLY : UNITTMP_
    use ModAMR, ONLY : local_cube,local_cubeBLK,availableBLKs
    use ModOctree
    ! use ModMpi
    implicit none

    type (adaptive_block_ptr) :: octree
    integer, intent(inout) :: BlksPerPE
    integer :: iPE, iBLK, iLEV, childNumber, numberBLK, ii, minPE, iError
    integer :: iLEVmin, iLEVmax
    logical, intent(inout) :: isRoot
    logical :: sol_blk_used
    type (adaptive_block_ptr) :: child
    real, dimension(4,3):: xyzends

    integer :: iChild
    !-------------------------------------------------------------------------

    if (associated(octree % ptr)) then
       read(UNITTMP_) numberBLK, childNumber, iPE, iBLK, iLEV, sol_blk_used


       octree % ptr % number  = numberBLK
       octree % ptr % child_number = childNumber
       octree % ptr % used    = sol_blk_used
       octree % ptr % refine  = .false.
       octree % ptr % coarsen = .false.
       octree % ptr % body    = .false.
       octree % ptr % IsExtraBoundaryOrPole = .false.  
       octree % ptr % IsOuterBoundary = .false.  
       iPE = octree % ptr % PE
       iBLK = octree % ptr % BLK

       if(RestartBlockLevels)then
          !        if (iProc == 0) &
          read(UNITTMP_) iLEVmin,iLEVmax

          !        call MPI_BCAST(iLEVmin,    1, MPI_INTEGER, 0, iComm, iError)
          !        call MPI_BCAST(iLEVmax,    1, MPI_INTEGER, 0, iComm, iError)

          octree % ptr % LEVmin = iLEVmin
          octree % ptr % LEVmax = iLEVmax
       end if

       if (octree % ptr % used) then
          global_block_ptrs(iBLK, iPE+1) % ptr => octree % ptr
          nBlockMax =  max(nBlockMax, iBLK)

          if (iProc == iPE) unusedBLK(iBLK) = .false.

          if(isRoot) then
             isRoot = .false.
          else
             availableBLKs(0,iPE)=availableBLKs(0,iPE)+1
          end if

          if(availableBLKs(0,iPE)<=BlksPerPE) then
             ! grab next from this processor
             local_cube = iPE
             local_cubeBLK = availableBLKs(0,iPE)
          else
             ! grab next from PE with minimum blocks
             minPE = iPE
             do ii=0,nProc-1
                if(availableBLKs(0,ii)<availableBLKs(0,minPE)) minPE=ii
             end do
             local_cube = minPE
             local_cubeBLK = availableBLKs(0,minPE)
          end if
       else
          local_cube = iPE
          local_cubeBLK = iBLK

          call refine_octree_block(octree, &
               local_cube, local_cubeBLK, -1,-1)

          do iChild = 1, 8
             child % ptr => octree % ptr % child(iChild)%ptr
             child % ptr % PE  = local_cube(1)
             child % ptr % BLK = local_cubeBLK(1)
             call read_octree_soln_block(child, BlksPerPE, isRoot)
          end do
       end if
    end if

  end subroutine read_octree_soln_block

  !============================================================================
  subroutine write_octree_file
    use ModProcMH
    use ModMain, ONLY : nBlockALL
    use ModParallel, ONLY : proc_dims
    use ModOctree
    use ModIoUnit, ONLY : UNITTMP_
    use ModIO, ONLY : write_prefix, iUnitOut
    implicit none

    integer :: i, j, k

    type (adaptive_block_ptr) :: octree

    character (len=4), Parameter :: octree_ext=".rst"

    if (iProc /= 0) return

    call write_prefix; write(iUnitOut,*) '=> Writing restart files ...'
    open(UNITTMP_, file=trim(NameRestartOutDir)//"octree"//octree_ext, &
         status="replace", form="UNFORMATTED")
    write(UNITTMP_) proc_dims(1),proc_dims(2),proc_dims(3)
    write(UNITTMP_) nBlockALL

    do k = 1, proc_dims(3)
       do j = 1, proc_dims(2)
          do i = 1, proc_dims(1)
             octree % ptr => octree_roots(i, j, k) % ptr
             call write_octree_soln_block(octree)
          end do
       end do
    end do

    close(UNITTMP_)

  end subroutine write_octree_file

  !===========================================================================
  recursive subroutine write_octree_soln_block(octree)
    use ModProcMH
    use ModOctree
    use ModIoUnit, ONLY : UNITTMP_
    implicit none

    type (adaptive_block_ptr) :: octree

    integer :: iPE, iBLK, iLEV, childNumber, numberBLK, icube
    integer :: iLEVmin, iLEVmax
    logical :: sol_blk_used

    type (adaptive_block_ptr) :: child

    if (associated(octree % ptr)) then
       numberBLK = octree % ptr % number
       childNumber = octree % ptr % child_number
       iPE = octree % ptr % PE
       iBLK = octree % ptr % BLK
       iLEV = octree % ptr % LEV
       sol_blk_used = octree % ptr % used
       iLEVmin = octree % ptr % LEVmin
       iLEVmax = octree % ptr % LEVmax

       write(UNITTMP_) numberBLK, childNumber, iPE, iBLK, iLEV, sol_blk_used
       write(UNITTMP_) iLEVmin,iLEVmax
       if (.not. octree % ptr % used) then
          do icube = 1, 8
             child % ptr => octree % ptr % child(iCube)%ptr
             call write_octree_soln_block(child)
          end do
       end if
    end if

  end subroutine write_octree_soln_block

end module ModRestartFile
