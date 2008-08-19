!=============================================================!
module SP_ModProc
  implicit none
  integer::iProc=-1,nProc=-1,iComm=-1
end module SP_ModProc
!=============================================================!
module SP_ModIhData
  use ModUtilities,ONLY:check_allocate
  use SP_ModMain
  implicit none
  save
  real::XyzLine_D(3)=(/0.0,0.0,0.0/)
  real::RBoundSC=1.2      !^CMP IF SC
  real::RBoundIH=21.0     !^CMP IF IH
  logical::DoRun=.true.,SaveMhData=.false.,DoReadMhData=.false.
  logical::DoRestart=.false.
  integer::nSmooth=0
contains
  function SP_xyz_i(iPoint)
    use CON_coupler
    implicit none
    real,dimension(3)::SP_xyz_i
    integer,intent(in)::iPoint
    SP_xyz_i= point_state_v('SP_XyzSP',3,iPoint)
  end function SP_xyz_i 
end module SP_ModIhData

!=============================================================!
subroutine SP_get_line_param(DsOut,&
                             XyzOut_D,&
                             RSCOut,&     !^CMP IF SC
                             RIHOut &     !^CMP IF IH
                             )
  use SP_ModIhData
  implicit none
  real,intent(out)::DsOut,XyzOut_D(3)
  real,intent(out)::RSCOut         !^CMP IF SC
  real,intent(out)::RIHOut         !^CMP IF IH
  DsOut=DsResolution
  XyzOut_D=XyzLine_D
  RSCOut=RBoundSC                  !^CMP IF SC
  RIHOut=RBoundIH                  !^CMP IF IH
end subroutine SP_get_line_param 
!=============================================================!
subroutine SP_put_input_time(TimeIn)
  use SP_ModIhData
  use CON_coupler
  implicit none
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
!=============================================================!
subroutine SP_set_param(CompInfo,TypeAction)
  use CON_comp_info
  use CON_coupler
  use ModReadParam
  use ModIOUnit,ONLY:STDOUT_
  use SP_ModIhData
  use SP_ModProc
  implicit none
  type(CompInfoType),intent(inout)       :: CompInfo
  character(len=*), intent(in)           :: TypeAction
  character (len=*), parameter :: NameSub='SP_set_param'
  ! The name of the command
  character (len=100) :: NameCommand
  logical,save::DoInit=.true.
  integer::iVerbose=0
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
     write(iStdOut,*)NameSub//': CHECK iSession =',i_session_read()
     if(DoInit)then
        DoInit=.false.
     end if
     do
        if(.not.read_line() ) EXIT
        if(.not.read_command(NameCommand)) CYCLE
        select case(NameCommand)
        case('#RESTART')
           DoRestart=.true.
        case('#LINE')
           call read_var('XyzLine_D(x_)',XyzLine_D(1))
           call read_var('XyzLine_D(y_)',XyzLine_D(2))
           call read_var('XyzLine_D(z_)',XyzLine_D(3))
           call read_var('RBoundSC'     ,RBoundSC    )!^CMP IF SC 
           call read_var('RBoundIH'     ,RBoundIH    )!^CMP IF IH
        case('#RTRANSIENT')
!           call read_var('rTransient',rTransient)
        case('#DORUN')
           call read_var('DoRun',DoRun)
        case('#SAVEMHDATA')
           call read_var('SaveMhData',SaveMhData)
        case('#DOREADMHDATA')
           call read_var('DoReadMhData',DoReadMhData)
        case('#NSTEP')
           call read_var('nStep',iDataSet)
        case('#TSIMULATION')
           call read_var('tSimulation',SP_Time)
        case('#PLOT')
!           call read_var('DnPlot',kfriss)
        case('#VERBOSE')
           call read_var('iVerbose',iVerbose)
           if(iVerbose>0)DoWriteAll=.true.
           if(DoWriteAll.and.iProc==0)&
                write(*,*)prefix,' Verbose everything'
        case('#NSMOOTH')
           call read_var('nSmooth',nSmooth)
        case default
           call CON_stop(NameSub//&
                ': Unknown command '&
                //NameCommand)
        end select
     end do
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
!======================================================================!
subroutine SP_put_from_mh(nPartial,&
     iPutStart,&
     Put,&
     W,&
     DoAdd,&
     Buff_I,nVar)
  use CON_router
  use SP_ModIhData
  implicit none
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
  use SP_ModIhData
  use ModIoUnit
  use ModConst 
  implicit none
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
  use ModNumConst
  use SP_ModIhData
  implicit none
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
!========================================================================!
subroutine SP_init_session(iSession,TimeIn)
  use SP_ModIhData
  implicit none
  integer,  intent(in) :: iSession    ! session number (starting from 1)
  real,     intent(in) :: TimeIn      ! seconds from start time
  logical,save::DoInit=.true.
  if(.not.DoInit)return
  DoInit=.false.
  if(SP_Time==cZero)SP_Time=TimeIn
  call SP_diffusive_shock(TypeAction='INIT')
end subroutine SP_init_session
!=============================================================!
subroutine SP_run(tInOut,tLimit)
  use SP_ModIhData
  use ModNumConst
  implicit none
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
!=============================================================!
subroutine SP_finalize(TimeSimulation)
  use SP_ModMain
  implicit none
  real,intent(in)::TimeSimulation
  call SP_diffusive_shock("FINALIZE")
end subroutine SP_finalize
!=============================================================!
subroutine SP_save_restart(TimeSimulation)
  use ModIOUnit,ONLY:io_unit_new
  use SP_ModIhData
  use CON_coupler
  implicit none
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
  return
10 call CON_stop('SP/restartOUT directory is not available')
end subroutine SP_save_restart
 
