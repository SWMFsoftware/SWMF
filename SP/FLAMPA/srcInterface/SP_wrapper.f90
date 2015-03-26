!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==================================================================
module SP_wrapper

  use ModUtilities,ONLY:check_allocate
  use SP_ModMain,ONLY: prefix, iStdOut, DoWriteAll, DoRun,&
       SaveMhData, DoReadMhData, XyzLine_D, RBoundSC, RBoundIH, nSmooth
  implicit none

  save

  private ! except

  public:: SP_set_param
  public:: SP_init_session
  public:: SP_run
  public:: SP_save_restart
  public:: SP_finalize

  ! coupling with MHD components
  public:: SP_get_line_param
  public:: SP_put_input_time
  public:: SP_put_from_mh


  ! Local variables



contains

  !=============================================================
  function SP_xyz_i(iPoint)
    use CON_coupler

    real,dimension(3)::SP_xyz_i
    integer,intent(in)::iPoint
    SP_xyz_i= point_state_v('SP_XyzSP',3,iPoint)
  end function SP_xyz_i
  !=============================================================
  subroutine SP_get_line_param(DsOut,&
       XyzOut_D,&
       RSCOut,&     !^CMP IF SC
       RIHOut &     !^CMP IF IH
       )
    use SP_ModMain, ONLY: DsResolution

    real,intent(out)::DsOut,XyzOut_D(3)
    real,intent(out)::RSCOut         !^CMP IF SC
    real,intent(out)::RIHOut         !^CMP IF IH
    DsOut=DsResolution
    XyzOut_D=XyzLine_D
    RSCOut=RBoundSC                  !^CMP IF SC
    RIHOut=RBoundIH                  !^CMP IF IH
  end subroutine SP_get_line_param
  !==================================================================
  subroutine SP_put_input_time(TimeIn)
    use SP_ModMain, ONLY: DInner_I, DataInputTime, nX, SP_allocate
    use CON_coupler

    real,intent(in)::TimeIn
    integer,dimension(2)::nU_I
    DataInputTime=TimeIn
    if(.not.allocated(DInner_I))then
       nU_I=ubound_vector('SP_XyzSP')
       nX=nU_I(2)
       call SP_allocate
    end if
    nX=0
  end subroutine SP_put_input_time
  !==================================================================
  subroutine SP_set_param(CompInfo,TypeAction)
    use CON_comp_info
    use CON_coupler
    use SP_ModMain, ONLY:iProc, nProc, iComm, DoInit, DoRestart, &
         iDataSet, SP_Time
    use ModConst, ONLY: rSun
    use ModIOUnit,ONLY:STDOUT_
    use SP_ModSetParam, ONLY:i_session_read, SP_set_parameters
    !-----------------------
    type(CompInfoType),intent(inout)       :: CompInfo
    character(len=*), intent(in)           :: TypeAction
    character (len=*), parameter :: NameSub='SP_set_param'
    !-----------------------------------------------------------!
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        =.true.,                        &
            NameVersion=&
            'SEP FLAMPA, Sokolov&Roussev ', &
            Version    =1.1)
    case('STDOUT')
       iStdOut=STDOUT_
       prefix='SP: '
       DoWriteAll=.false.
    case('FILEOUT')
       call get(CompInfo,iUnitOut=iStdOut)
       prefix=''
    case('CHECK')
       write(iStdOut,*)NameSub//': CHECK iSession =',i_session_read()
       if(DoInit)then
          DoInit=.false.
       end if
       if(DoRestart)then
          write(iStdOut,*)prefix,'Restart from nStep=',iDataSet,&
               ', tSimulation=',SP_Time
          DoRestart=.false.
       end if
    case('READ')
       call SP_set_parameters(TypeAction)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       if(nProc/=1)call CON_stop(&
            'The present SEP version can not use more than 1 PE')
    case('GRID')
       !Initialize grid
       call init_decomposition(SP_,SP_,1)
       call set_coord_system(SP_,'HGI',rSun)
    case default
    end select
  end subroutine SP_set_param
  !===========================================================================
  subroutine SP_put_from_mh(nPartial,&
       iPutStart,&
       Put,&
       W,&
       DoAdd,&
       Buff_I,nVar)
    use SP_ModMain, ONLY: nX, State_VI, X_DI
    use CON_router

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer::iCell
    real:: Weight
    iCell=Put%iCB_II(1,iPutStart)
    Weight=W%Weight_I(iPutStart)
    nX=max(nX,iCell)
    if(DoAdd)then
       State_VI(:,iCell)=State_VI(:,iCell)+Buff_I(:)*Weight
    else
       !Get the lagrangian mesh coordinates from global vector
       X_DI(:,iCell)=&
            SP_xyz_i(iCell)
       State_VI(:,iCell)=Buff_I(:)*Weight
    end if
  end subroutine SP_put_from_mh
  !======================================================================
  subroutine SP_save_mhdata
    use SP_ModMain, ONLY: iDataSet, nX, DataInputTime, X_DI, State_VI
    use ModIoUnit
    use ModConst 

    character(LEN=*),parameter::IO_dir='./SP/MHDATA/'
    character(LEN=50)::NameFile
    integer::iFile,iPoint,i
    write(NameFile,'(a,i4.4,a)')trim(IO_dir)//'mhdata_',iDataSet,'.dat'
    iFile=io_unit_new()
    open(iFile,file=trim(NameFile),status='unknown',&
         form='formatted')
    write(iFile,'(a,f13.6,a,i6)')'Time_Simulation=  ',DataInputTime,'         ,  n_Step=',iDataSet
    do iPoint=1,nX
       write(iFile,'(11(e13.6,2x))')X_DI(:,iPoint)/rSun,State_VI(1:8,iPoint)
    end do
    close(iFile)
  end subroutine SP_save_mhdata
  !==========================================================================

  subroutine sp_smooth_data
    use SP_ModMain, ONLY:State_VI, nX
    use ModNumConst


    integer::i,iVar
    do iVar=1,8
       SMOOTH:do i=1,nSmooth
          if(.not.any(&
               2.0*&
               abs(State_VI(iVar,1:nX-2)-State_VI(iVar,2:nX-1))<&
               abs(State_VI(iVar,3:nX  )-State_VI(iVar,2:nX-1)).or.&
               abs(State_VI(iVar,1:nX-2)-State_VI(iVar,2:nX-1))>&
               2.0*&
               abs(State_VI(iVar,3:nX  )-State_VI(iVar,2:nX-1))))&
               EXIT SMOOTH
          where(2.0*&
               abs(State_VI(iVar,1:nX-2)-State_VI(iVar,2:nX-1))<&
               abs(State_VI(iVar,3:nX  )-State_VI(iVar,2:nX-1)).or.&
               abs(State_VI(iVar,1:nX-2)-State_VI(iVar,2:nX-1))>&
               2.0*&
               abs(State_VI(iVar,3:nX  )-State_VI(iVar,2:nX-1)))&
               State_VI(iVar,2:nX-1)=0.50*(State_VI(iVar,1:nX-2)+&
               State_VI(iVar,3:nX))
       end do SMOOTH
    end do
  end subroutine sp_smooth_data
  !==================================================================
  subroutine SP_init_session(iSession,TimeIn)
    use ModNumConst, ONLY: cZero
    use SP_ModMain,ONLY: SP_Time, SP_diffusive_shock
    integer,  intent(in) :: iSession    ! session number (starting from 1)
    real,     intent(in) :: TimeIn      ! seconds from start time
    logical,save::DoInit=.true.
    if(.not.DoInit)return
    DoInit=.false.
    if(SP_Time==cZero)SP_Time=TimeIn
    call SP_diffusive_shock(TypeAction='INIT')
  end subroutine SP_init_session
  !==================================================================
  subroutine SP_run(tInOut,tLimit)
    use SP_ModReadMhData,ONLY:read_ihdata_for_sp,mh_transform_for_flampa
    use ModNumConst
    use SP_ModMain, ONLY: iDataSet, DataInputTime, SP_diffusive_shock

    real,intent(inout)::tInOut
    real,intent(in):: tLimit
    character(LEN=15):: NameFile
    if(.not.DoReadMhData)then
       if(nSmooth>0)call SP_smooth_data
       call mh_transform_for_flampa
       tInOut=max(tInOut,tLimit)
    else
       call read_ihdata_for_sp(1,0)
       IF(nSmooth>0)call SP_smooth_data
       tInOut=DataInputTime
    end if
    if(SaveMhData)call SP_save_mhdata
    if(DoRun)call  SP_diffusive_shock("RUN",DataInputTime)
    iDataSet=iDataSet+1
  end subroutine SP_run
  !==================================================================
  subroutine SP_finalize(TimeSimulation)
    use SP_ModMain, ONLY: SP_diffusive_shock

    real,intent(in)::TimeSimulation
    call SP_diffusive_shock("FINALIZE")
  end subroutine SP_finalize
  !==================================================================
  subroutine SP_save_restart(TimeSimulation)
    use SP_ModMain, ONLY: SP_Time, iDataSet
    use ModIOUnit,ONLY:io_unit_new
    use CON_coupler

    real,     intent(in) :: TimeSimulation
    integer::iFile
    !-----------------------------------------------
    iFile=io_unit_new()
    open(iFile,FILE='SP/restartOUT/restart.H',&
         STATUS='replace',ERR=10)
    write(iFile,*)
    write(iFile,'(a)')'#RESTART'
    write(iFile,*)
    write(iFile,'(a)')'#TSIMULATION'
    write(iFile,*)SP_Time
    write(iFile,*)
    write(iFile,'(a)')'#NSTEP'
    write(iFile,*)iDataSet
    write(iFile,*)
    close(iFile)
    !  call SP_save_f
    !  if(.not.SaveMhData)call SP_save_mh_data
    RETURN
10  call CON_stop('SP/restartOUT directory is not available')
  end subroutine SP_save_restart

end module SP_wrapper
