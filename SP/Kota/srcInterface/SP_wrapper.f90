!==============================================================================
module SP_ModIhData
  implicit none
  real,allocatable,dimension(:,:),save::State_VI
  real,allocatable,dimension(:,:),save::Xyz_DI
  integer::nPoint
  logical::DoInit=.true.
  real::DsResolution=0.1,XyzLine_D(3)
  real::RBoundSC=1.0,RBoundIH=1.0
  real,save::Time_Simulation,DataInputTime
  integer::n_Step
end module SP_ModIhData
!==============================================================================
subroutine SP_put_input_time(TimeIn)
  use SP_ModIhData,ONLY:DataInputTime
  implicit none
  real,intent(in)::TimeIn
  DataInputTime=TimeIn
end subroutine SP_put_input_time
!==============================================================================
subroutine sp_set_ihdata(nPointIn,XyzIn_DI)
  use SP_ModIhData
  use ModUtilities
  implicit none
  integer,intent(in)::nPointIn
  real,dimension(3,nPoint),intent(in)::XyzIn_DI
  integer::iError
  nPoint=nPointIn
  
  if(DoInit)then
     allocate(Xyz_DI(3,nPoint),stat=iError)
     call check_allocate(iError,&
            'Xyz_DI in sp_set_ihdata')
     allocate(State_VI(8,nPoint),stat=iError)
     call check_allocate(iError,&
            'State_VI in sp_set_ihdata')
     DoInit=.false.
  end if
  Xyz_DI(:,1:nPoint)=XyzIn_DI(:,1:nPoint)
end subroutine sp_set_ihdata
!==================================================
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

subroutine SP_get_line_param(DsOut,XyzOut_D,RSCOut,RIHOut)
  use SP_ModIhData
  implicit none
  real,intent(out)::DsOut,XyzOut_D(3),RSCOut,RIHOut
  DsOut=DsResolution
  XyzOut_D=XyzLine_D
  RSCOut=RBoundSC
  RIHOut=RBoundIH
end subroutine SP_get_line_param


!======================================================================
subroutine write_ihdata
  use SP_ModIhData
  use ModIoUnit
!  use ModIO,  ONLY: IO_dir
!  use ModMain,ONLY:iCoupleSP,n_step,Time_Simulation
  implicit none
  character(LEN=50)::NameFile
  integer::iFile,iPoint,i
!  write(NameFile,'(a,i4.4,a)')trim(IO_dir)//'ihdata_',iCoupleSP,'.dat'
  iFile=io_unit_new()
  open(iFile,file=trim(NameFile),status='unknown',&
          form='formatted')
  write(iFile,*)'Time_Simulation=',Time_Simulation,',  n_step=',n_step
  do iPoint=1,nPoint
     write(iFile,'(11(e13.6,2x))')Xyz_DI(:,iPoint),State_VI(1:8,iPoint)
  end do
  close(iFile)
end subroutine write_ihdata
!==========================================================================

subroutine sp_smooth_ihdata
  use ModNumConst
  use SP_ModIhData
!  use ModMain,ONLY:DsResolution
  implicit none
  integer::i,iVar
  do i=1,int(cHalf/DsResolution)
     do iVar=1,8
        where(cTwo*&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))<&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1)).or.&
             abs(State_VI(iVar,1:nPoint-2)-State_VI(iVar,2:nPoint-1))>&
             cTwo*&
             abs(State_VI(iVar,3:nPoint  )-State_VI(iVar,2:nPoint-1)))&
        State_VI(iVar,2:nPoint-1)=cHalf*(State_VI(iVar,1:nPoint-2)+&
             State_VI(iVar,3:nPoint))
     end do
  end do
end subroutine sp_smooth_ihdata      
!======================================================================
subroutine SP_init_session(iSession,TimeSimulation)
  use SP_ModIhData,ONLY:Time_Simulation
  implicit none
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  Time_Simulation=Time_Simulation
end subroutine SP_init_session
!======================================================================
subroutine SP_finalize(TimeSimulation)
  implicit none
  real,intent(in)::TimeSimulation
  call closetime
end subroutine SP_finalize
!======================================================================
subroutine SP_set_param(CompInfo,TypeAction)
  use CON_comp_info
  use ModIOUnit,ONLY:STDOUT_
  implicit none
  include 'stdout.h'
  type(CompInfoType),intent(inout)       :: CompInfo
  character(len=*), intent(in)           :: TypeAction
   select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true.,                        &
          NameVersion=&
          'SEP acceleration with pitch-angle dependance, Kota(2004)', &
          Version    =1.1)
  case('STDOUT')
     iStdOut=STDOUT_
     prefix='SP'
  case('FILEOUT')
     call get(CompInfo,iUnitOut=iStdOut)
     prefix=''
  case default
  end select
end subroutine SP_set_param
!=========================================================
subroutine SP_save_restart(TimeSimulation) 
  implicit none
  real,     intent(in) :: TimeSimulation 
end subroutine SP_save_restart
