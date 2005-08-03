!=============================================================!
module SP_ModProc
  implicit none
  integer::iProc=-1,nProc=-1,iComm=-1
end module SP_ModProc
!=============================================================!
module SP_ModIhData
  use ModNumConst
  use ModUtilities,ONLY:check_allocate
  implicit none
  include 'stdout.h'
  real,allocatable,dimension(:,:),save::State_VI
  real,allocatable,dimension(:,:),save::Xyz_DI
  integer::nPoint
  logical::DoInitArr=.true.
  integer,parameter::nResolution=10
  real,parameter::DsResolution=cOne/nResolution
  real::XyzLine_D(3)=(/cZero,cZero,cZero/)
  real::RBoundSC=1.1      !^CMP IF SC
  real::RBoundIH=21.0     !^CMP IF IH
  real::tSimulation=cZero,DataInputTime=cZero
  logical::DoRun=.true.,SaveMhData=.false.,DoReadMhData=.false.
  logical::DoRestart
  integer::nSmooth=0
  integer::nStep=0 
end module SP_ModIhData
!==============================================================================
subroutine sp_set_ihdata(nPointIn,XyzIn_DI)
  use SP_ModIhData
  implicit none
  integer,intent(in)::nPointIn
  real,dimension(3,nPointIn),intent(in)::XyzIn_DI
  integer::iError
  nPoint=nPointIn
  
  if(DoInitArr)then
     allocate(Xyz_DI(3,nPoint),stat=iError)
     call check_allocate(iError,&
            'Xyz_DI in sp_set_ihdata')
     allocate(State_VI(8,nPoint),stat=iError)
     call check_allocate(iError,&
            'State_VI in sp_set_ihdata')
     DoInitArr=.false.
  end if
  Xyz_DI(:,1:nPoint)=XyzIn_DI(:,1:nPoint)
end subroutine sp_set_ihdata
!=============================================================!
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
  use CON_axes
  implicit none
  real,intent(in)::TimeIn
  DataInputTime=TimeIn
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
        write(iStdOut,*)prefix,'Restart from nStep=',nStep,&
             ', tSimulation=',tSimulation
!        call SP_initial
!        call SP_read_mh_data
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
!           call read_var('DoRun',DoRun)
        case('#SAVEMHDATA')
           call read_var('SaveMhData',SaveMhData)
        case('#DOREADMHDATA')
           call read_var('DoReadMhData',DoReadMhData)
        case('#NSTEP')
           call read_var('nStep',nStep)
        case('#TSIMULATION')
           call read_var('tSimulation',tSimulation)
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
     call set_coord_system(SP_,'HGI',cOne)
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
  if(DoAdd)then
     State_VI(:,iCell)=State_VI(:,iCell)+Buff_I(:)*Weight
  else
     State_VI(:,iCell)=Buff_I(:)*Weight
  end if
end subroutine SP_put_from_mh
!======================================================================
subroutine write_ihdata
  use SP_ModIhData
  use ModIoUnit
  implicit none
  character(LEN=*),parameter::IO_dir='./SP/'
  character(LEN=50)::NameFile
  integer::iFile,iPoint,i
  write(NameFile,'(a,i4.4,a)')trim(IO_dir)//'ihdata_',nStep,'.dat'
  iFile=io_unit_new()
  open(iFile,file=trim(NameFile),status='unknown',&
          form='formatted')
  write(iFile,*)'tSimulation=',tSimulation,',  nStep=',nStep
  do iPoint=1,nPoint
     write(iFile,'(11(e13.6,2x))')Xyz_DI(:,iPoint),State_VI(1:8,iPoint)
  end do
  close(iFile)
end subroutine write_ihdata
!==========================================================================

subroutine sp_smooth_ihdata
  use ModNumConst
  use SP_ModIhData
  implicit none
  integer::i,iVar
  do iVar=1,8
     SMOOTH:do i=1,int(cHalf/DsResolution)
        if(.not.any(&
             cTwo*&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))<&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1)).or.&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))>&
             cTwo*&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1))))&
             EXIT SMOOTH
        where(cTwo*&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))<&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1)).or.&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))>&
             cTwo*&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1)))&
             State_VI(iVar,2:nPoint-1)=cHalf*(State_VI(iVar,1:nPoint-2)+&
             State_VI(iVar,3:nPoint))
     end do SMOOTH
  end do
end subroutine sp_smooth_ihdata      
!========================================================================!
subroutine SP_init_session(iSession,TimeIn)
  use SP_ModIhData
  implicit none
  integer,  intent(in) :: iSession    ! session number (starting from 1)
  real,     intent(in) :: TimeIn      ! seconds from start time
  logical,save::DoInit=.true.
  integer::      jnext,jsep,jstep,iStep
  common /SP_spmain/jnext,jsep,jstep,istep
  if(.not.DoInit)return
  DoInit=.false.
  if(tSimulation==cZero)tSimulation=TimeIn
end subroutine SP_init_session
!=============================================================!
subroutine SP_run(tInOut,tLimit)
  use SP_ModIhData,ONLY:tSimulation,DataInputTime,nStep,nSmooth
!  use SP_ModIhData,ONLY:DoRun,SaveMhData,DoReadMhData
  use ModNumConst
  implicit none
  real,intent(inout)::tInOut
  real,intent(in):: tLimit
  nStep=nStep+1
!  if(.not.DoReadMhData)then
 !    if(nSmooth>0)call SP_smooth_data
     tInOut=max(tInOut,tLimit)
!  else
 !    call SP_read_mh_data
     tInOut=DataInputTime
!  end if
 ! if(SaveMhData)call SP_save_mh_data
 ! if(DoRun)call SP_sharpen_and_run(tSimulation,DataInputTime)
end subroutine SP_run
!=============================================================!
subroutine SP_finalize(TimeSimulation)
  implicit none
  real,intent(in)::TimeSimulation
!  call SP_closetime
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
  write(iFile,*)tSimulation
  write(iFile,*)
  write(iFile,'(a)')'#NSTEP'
  write(iFile,*)nStep
  write(iFile,*)
  close(iFile)
!  call SP_save_f
!  if(.not.SaveMhData)call SP_save_mh_data
  return
10 call CON_stop('SP/restartOUT directory is not available')
end subroutine SP_save_restart
 
