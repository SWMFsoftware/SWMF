!^CFG COPYRIGHT UM
!==============================================================================
subroutine write_restart_header
  use ModProcMH
  use ModMain
  use ModVarIndexes, ONLY: NameEquation, nVar
  use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2
  use ModGeometry, ONLY: XyzMin_D, XyzMax_D, &             !^CFG IF COVARIANT
       TypeGeometry, UseCovariant, UseVertexBasedGrid      !^CFG IF COVARIANT
  use ModParallel, ONLY: proc_dims
  use ModImplicit, ONLY: n_prev, dt_prev                   !^CFG IF IMPLICIT
  use ModUser,     ONLY: NameUserModule, VersionUserModule
  use ModPhysics
  use ModIO

  implicit none

  integer:: iProblemType=11
  !--------------------------------------------------------------------------

  if (iProc/=0) RETURN
  iProblemType = problem_type

  open(unit_tmp,file=trim(NameRestartOutDir)//'restart.H')

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
  write(unit_tmp,*)
  write(unit_tmp,'(a)')'#PROBLEMTYPE'
  write(unit_tmp,'(i8,a32)') iProblemType,'iProblemType'
  write(unit_tmp,*)
  write(unit_tmp,'(a)')'#NEWRESTART'
  write(unit_tmp,'(l1,a39)') UseConstrainB,'DoRestartBFace' !^CFG IF CONSTRAINB
  write(unit_tmp,*)
  write(unit_tmp,'(a)')'#BLOCKLEVELSRELOADED'
  write(unit_tmp,*)
  write(unit_tmp,'(a)')'#NSTEP'
  write(unit_tmp,'(i8,a32)')n_step,'nStep'
  write(unit_tmp,*)
  if(n_prev == n_step)then                              !^CFG IF IMPLICIT BEGIN
     write(unit_tmp,'(a)')'#NPREVIOUS'
     write(unit_tmp,'(i8,a32)')n_prev,'nPrev'
     write(unit_tmp,'(1pe20.12,a20)')dt_prev,'DtPrev'
     write(unit_tmp,*)
  end if                                                !^CFG END IMPLICIT
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
  if(UseCovariant)then                        !^CFG IF COVARIANT BEGIN
     write(unit_tmp,'(a)')'#COVARIANTGEOMETRY'
     write(unit_tmp,'(a)')trim(TypeGeometry)
     write(unit_tmp,*)
     write(unit_tmp,'(a)')'#VERTEXBASEDGRID'
     write(unit_tmp,'(l1,a39)') UseVertexBasedGrid,'UseVertexBasedGrid'
     write(unit_tmp,*)
     write(unit_tmp,'(a)')'#LIMITGENCOORD1'                   
     write(unit_tmp,'(1pe13.5,a27)')XyzMin_D(1),'XyzMin_D(1)' 
     write(unit_tmp,'(1pe13.5,a27)')XyzMax_D(1),'XyzMax_D(1)' 
     write(unit_tmp,*)
  end if                                      !^CFG END COVARIANT

  if(NameThisComp=='GM')then
     write(unit_tmp,'(a)')'#SOLARWIND'
     write(unit_tmp,'(1pe15.7,a25)')SW_rho_dim,'SwRhoDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_T_dim,  'SwTDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_Ux_dim, 'SwUxDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_Uy_dim, 'SwUyDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_Uz_dim, 'SwUzDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_Bx_dim, 'SwBxDdim'
     write(unit_tmp,'(1pe15.7,a25)')SW_By_dim, 'SwByDim'
     write(unit_tmp,'(1pe15.7,a25)')SW_Bz_dim, 'SwBzDim'
     write(unit_tmp,*)
  end if
  if(body1)then
     write(unit_tmp,'(a)')'#BODY'
     write(unit_tmp,'(l1,a39)')body1,'UseBody'
     write(unit_tmp,'(1pe13.5,a27)')Rbody,'rBody'
     if(NameThisComp=='GM')then
        write(unit_tmp,'(1pe13.5,a27)')Rcurrents,   'rCurrents'
        write(unit_tmp,'(1pe13.5,a27)')Body_rho_dim,'BodyRhoDim'
        write(unit_tmp,'(1pe13.5,a27)')Body_T_dim,  'BodyTDim'
     end if
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
  write(unit_tmp,'(6(1pe13.5),a)') &
       unitUSER_x, unitUSER_t, unitUSER_rho,&
       unitUSER_rhoU,unitUSER_B, unitUSER_P,' unitUSER_ x,t,rho,rhoU,B,p'
  write(unit_tmp,'(6(1pe13.5),a)') &
       unitSI_x, unitSI_t, unitSI_rho, &
       unitSI_rhoU, unitSI_B, unitSI_P,     ' unitSI_ x,t,rho,rhoU,B,p'

  close(unit_tmp)

end subroutine write_restart_header
!==============================================================================

subroutine read_restart_file
  use ModProcMH,     ONLY: iProc
  use ModIO,         ONLY: NameRestartInDir, FileName, Unit_Tmp, nByteRealRead
  use ModMain,       ONLY: GlobalBlk, Global_Block_Number, nI, nJ, nK, Gcn, &
       ProcTest, BlkTest, iTest, jTest, kTest, &
       n_step, dt_BLK, Cfl, CodeVersion
  use ModVarIndexes, ONLY: nVar, DefaultState_V
  use ModAdvance,    ONLY: State_VGB
  use ModGeometry,   ONLY: dx_BLK, dy_BLK, dz_BLK, xyzStart_BLK
  use ModParallel,   ONLY: iBlockRestartALL_A
  use ModCT,         ONLY: BxFace_BLK,ByFace_BLK,BzFace_BLK !^CFG IF CONSTRAINB
  use ModIO,         ONLY: Restart_BFace                    !^CFG IF CONSTRAINB
  use ModImplicit,   ONLY: n_prev, w_prev                   !^CFG IF IMPLICIT
  use ModKind,       ONLY: Real4_, Real8_
  implicit none

  integer :: iVar, i, j, k, iError
  real    :: tSimulationRead
  character (len=4), Parameter :: restart_ext=".rst"

  ! Temporaray variables to read arbitrary precision data files
  real (Real8_) :: Dt8, Time8, Dxyz8_D(3), Xyz8_D(3)
  real (Real4_) :: Dt4, Time4, Dxyz4_D(3), Xyz4_D(3)
  real (Real8_) :: State8_CV(1:nI,1:nJ,1:nK,nVar)
  real (Real4_) :: State4_CV(1:nI,1:nJ,1:nK,nVar)
  !^CFG IF CONSTRAINB BEGIN
  real (Real8_) :: b8_X(1:nI+1,1:nJ,1:nK), b8_Y(1:nI,1:nJ+1,1:nK), &
       b8_Z(1:nI,1:nJ,1:nK+1)
  real (Real4_) :: b4_X(1:nI+1,1:nJ,1:nK), b4_Y(1:nI,1:nJ+1,1:nK), &
       b4_Z(1:nI,1:nJ,1:nK+1)
  !^CFG END CONSTRAINB

  character (len=*), parameter :: NameSub='read_restart_file'
  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------
  if(iProc==PROCtest.and.globalBLK==BLKtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  write(filename,'(a,i5.5,a)') &
       trim(NameRestartInDir)//'blk',&
       iBlockRestartALL_A(global_block_number(globalBLK)),&
       restart_ext

  open(unit_tmp, file=filename, status='old', form='UNFORMATTED',&
       iostat = iError)

  if(iError /= 0) call stop_mpi(NameSub// &
       ' read_restart_file could not open: '//trim(filename))

  ! Fill in ghost cells
  do k=1-gcn,nK+gcn; do j=1-gcn,nK+gcn; do i=1-gcn,nK+gcn
     State_VGB(1:nVar, i, j, k, globalBLK) = DefaultState_V
  end do;end do;end do

  ! Do not overwrite time_simulation which is read from restart.H
  if(nByteRealRead == 8)then
     read(unit_tmp, iostat = iError) Dt8, Time8
     dt_BLK(globalBLK) = Dt8
     tSimulationRead   = Time8

     read(unit_tmp, iostat = iError) Dxyz8_D, Xyz8_D
     Dx_BLK(globalBLK) = Dxyz8_D(1)
     Dy_BLK(globalBLK) = Dxyz8_D(2)
     Dz_BLK(globalBLK) = Dxyz8_D(3)
     XyzStart_BLK(:,globalBLK) = Xyz8_D

     read(Unit_tmp, iostat = iError) State8_CV
     do iVar = 1, nVar
        State_VGB(iVar,1:nI,1:nJ,1:nK,globalBLK) = State8_CV(:,:,:,iVar)
     end do

     if(Restart_Bface)then                            !^CFG IF CONSTRAINB BEGIN
        read(Unit_tmp, iostat = iError) b8_X, b8_Y, b8_Z               
        BxFace_BLK(1:nI+1,1:nJ,1:nK,globalBLK) = b8_X
        ByFace_BLK(1:nI,1:nJ+1,1:nK,globalBLK) = b8_Y
        BzFace_BLK(1:nI,1:nJ,1:nK+1,globalBLK) = b8_Z
     end if                                           !^CFG END CONSTRAINB
     if(n_prev==n_step) then                          !^CFG IF IMPLICIT BEGIN
        read(Unit_tmp, iostat = iError) State8_CV
        w_prev(:,:,:,:,globalBLK) = State8_CV
     end if                                           !^CFG END IMPLICIT
  else
     read(unit_tmp, iostat = iError) Dt4, Time4
     dt_BLK(globalBLK) = Dt4
     tSimulationRead   = Time4

     read(unit_tmp, iostat = iError) Dxyz4_D, Xyz4_D
     Dx_BLK(globalBLK) = Dxyz4_D(1)
     Dy_BLK(globalBLK) = Dxyz4_D(2)
     Dz_BLK(globalBLK) = Dxyz4_D(3)
     XyzStart_BLK(:,globalBLK) = Xyz4_D

     read(Unit_tmp, iostat = iError) State4_CV
     do iVar = 1, nVar
        State_VGB(iVar,1:nI,1:nJ,1:nK,globalBLK) = State4_CV(:,:,:,iVar)
     end do

     if(Restart_Bface)then                            !^CFG IF CONSTRAINB BEGIN
        read(Unit_tmp, iostat = iError) b4_X, b4_Y, b4_Z               
        BxFace_BLK(1:nI+1,1:nJ,1:nK,globalBLK) = b4_X
        ByFace_BLK(1:nI,1:nJ+1,1:nK,globalBLK) = b4_Y
        BzFace_BLK(1:nI,1:nJ,1:nK+1,globalBLK) = b4_Z
     end if                                           !^CFG END CONSTRAINB
     if(n_prev==n_step) then                          !^CFG IF IMPLICIT BEGIN
        read(Unit_tmp, iostat = iError) State4_CV
        w_prev(:,:,:,:,globalBLK) = State4_CV
     end if                                           !^CFG END IMPLICIT
  endif

  if(iError /= 0) call stop_mpi(NameSub// &
       ' could not read data from '//trim(filename))

  close(unit_tmp)

  if(CodeVersion>5.60 .and. CodeVersion <7.00) &
       dt_BLK(globalBLK)=dt_BLK(globalBLK)/cfl

  call correctE

  if(DoTestMe)then
     write(*,*)NameSub,': iProc, globalBLK =',iProc, globalBLK
     write(*,*)NameSub,': dt,tSimRead =',dt_BLK(globalBLK),tSimulationRead
     write(*,*)NameSub,': dx,dy,dz_BLK=',dx_BLK(globalBLK),dy_BLK(globalBLK),&
          dz_BLK(globalBLK)
     write(*,*)NameSub,': xyzStart_BLK=',xyzStart_BLK(:,globalBLK)
     write(*,*)NameSub,': State_VGB   =', &
          State_VGB(:,Itest,Jtest,Ktest,globalBLK)
     write(*,*)NameSub,' finished'
  end if

end subroutine read_restart_file

!=============================================================================

subroutine write_restart_file
  use ModProcMH
  use ModMain
  use ModVarIndexes
  use ModAdvance, ONLY : State_VGB
  use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, xyzStart_BLK
  use ModCT, ONLY : BxFace_BLK,ByFace_BLK,BzFace_BLK       !^CFG IF CONSTRAINB
  use ModImplicit, ONLY: n_prev, w_prev                    !^CFG IF IMPLICIT
  use ModIO
  implicit none

  character (len=4), Parameter :: restart_ext=".rst"
  integer::iVar
  !--------------------------------------------------------------------

  if (global_block_number(globalBLK) <= 0 ) then
     write(filename,'(a,i4.4,a,i3.3,a)') &
          trim(NameRestartOutDir)//'pe',iProc,'_blk',globalBLK,restart_ext
  else
     write(filename,'(a,i5.5,a)') &
          trim(NameRestartOutDir)//'blk',&
          global_block_number(globalBLK),restart_ext
  end if

  open(unit_tmp, file=filename, status="replace", form='UNFORMATTED')

  write(Unit_tmp)  dt_BLK(globalBLK),time_Simulation
  write(Unit_tmp) &
       dx_BLK(globalBLK),dy_BLK(globalBLK),dz_BLK(globalBLK),&
       xyzStart_BLK(:,globalBLK)
  write(Unit_tmp) &
       ( State_VGB(iVar,1:nI,1:nJ,1:nK,globalBLK), iVar=1,nVar)
  if(UseConstrainB)then                            !^CFG iF CONSTRAINB BEGIN
     write(Unit_tmp) &
          BxFace_BLK(1:nI+1,1:nJ,1:nK,globalBLK),&
          ByFace_BLK(1:nI,1:nJ+1,1:nK,globalBLK),&
          BzFace_BLK(1:nI,1:nJ,1:nK+1,globalBLK)
  end if                                           !^CFG END CONSTRAINB
  if(n_prev==n_step) write(Unit_tmp) &             !^CFG IF IMPLICIT
       w_prev(:,:,:,:,globalBLK)                   !^CFG IF IMPLICIT
  close(unit_tmp)

end subroutine write_restart_file
