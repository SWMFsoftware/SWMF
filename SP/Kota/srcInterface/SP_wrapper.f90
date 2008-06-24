!^CFG COPYRIGHT UM
!=============================================================!
module SP_Mod
  use ModConst
  implicit none
  include 'coupler.h'
  include 'stdout.h'
  save
  logical:: UseSelfSimilarity,UseRefresh
  common/SP_log/UseSelfSimilarity,UseRefresh
  real,parameter::DsResolution=cOne/nResolution
  real::XyzLine_D(3)=(/cZero,cZero,cZero/)
  real::RBoundSC=1.1      !^CMP IF SC
  real::RBoundIH=21.0     !^CMP IF IH
  real::tSimulation=cZero,DataInputTime=cZero
  integer::nStep=0
  logical::DoRun=.true.,SaveMhData=.false.,DoReadMhData=.false.
  integer::iProc=-1,nProc=-1,iComm=-1
  integer::nSmooth=0
contains
  function SP_xyz_i(iPoint)
    use CON_coupler
    implicit none
    real,dimension(3)::SP_xyz_i
    integer,intent(in)::iPoint
    SP_xyz_i= point_state_v('SP_XyzSP',3,iPoint)
  end function SP_xyz_i
end module SP_Mod
!=============================================================!
subroutine SP_set_param(CompInfo,TypeAction)
  use CON_comp_info
  use CON_coupler
  use ModReadParam
  use ModIOUnit,ONLY:STDOUT_
  use SP_Mod
  implicit none
  type(CompInfoType),intent(inout)       :: CompInfo
  character(len=*), intent(in)           :: TypeAction
  character (len=*), parameter :: NameSub='SP_set_param'
  ! The name of the command
  character (len=100) :: NameCommand
  logical,save::DoInit=.true.
  integer::iVerbose=0
  integer::kfriss,kacc
  real:: time,tmax,dlnt0,dt1,dta
  common /SP_times/  time,tmax,dlnt0,dt1,dta,kfriss,kacc
  !-----------------------------------------------------------!
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true.,                        &
          NameVersion=&
          'SEP with pitch-angle scattering (Kota)', &
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
        call sp_set_defaults
        DoInit=.false.
     end if
     if(DoRestart)then
        write(iStdOut,*)prefix,'Restart from nStep=',nStep,&
             ', tSimulation=',tSimulation
        call SP_initial
        call SP_read_mh_data
        DoRestart=.false.
     end if
  case('READ')
     write(iStdOut,*)NameSub//': CHECK iSession =',i_session_read()
     if(DoInit)then
        call sp_set_defaults
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
           call read_var('rTransient',rTransient)
        case('#DORUN')
           call read_var('DoRun',DoRun)
        case('#SAVEMHDATA')
           call read_var('SaveMhData',SaveMhData)
        case('#DOREADMHDATA')
           call read_var('DoReadMhData',DoReadMhData)
        case('#NSTEP')
           call read_var('nStep',nStep)
        case('#TSIMULATION')
           call read_var('tSimulation',tSimulation)
        case('#PLOT')
           call read_var('DnPlot',kfriss)
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

     call set_coord_system(SP_,'HGI',cAU)
  case default
  end select
end subroutine SP_set_param
!======================================================================!
subroutine SP_init_session(iSession,TimeIn)
  use SP_Mod
  implicit none
  integer,  intent(in) :: iSession    ! session number (starting from 1)
  real,     intent(in) :: TimeIn      ! seconds from start time
  logical,save::DoInit=.true.
  integer::      jnext,jsep,jstep,iStep
  common /SP_spmain/jnext,jsep,jstep,istep
  if(.not.DoInit)return
  DoInit=.false.
  UseSelfSimilarity=.false.
  UseRefresh=.false.
  if(tSimulation==cZero)tSimulation=TimeIn
  DoRestart=.false.
  call SP_init
  call SP_clean_coupler
  iMax=0
  iShock=1
  if(jSep/=jStep)then
     jSep=jStep
     if(iProc==0)&
          write(iStdOut,*)prefix, 'SWMF resets jsep=',jSep
  end if
end subroutine SP_init_session
!=============================================================!
subroutine SP_get_line_param(DsOut,&
                             XyzOut_D,&
                             RSCOut,&     !^CMP IF SC
                             RIHOut &     !^CMP IF IH
                             )
  use SP_Mod
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
  use SP_Mod
  use CON_coupler
  use CON_axes
  implicit none
  real,intent(in)::TimeIn
  DataInputTime=TimeIn
end subroutine SP_put_input_time
!=============================================================!
subroutine SP_put_from_mh(nPartial,&
     iPutStart,&
     Put,&
     W,&
     DoAdd,&
     Buff_I,nVar)
  use SP_Mod
  use CON_router
  implicit none
  integer,intent(in)::nPartial,iPutStart,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::Buff_I
!==========================================================
!============INTERFACE TO F77 CODE=========================
!The data which are accepted by SEP module from the corona or
!inner heliosphere
!with the identifier old: storage to save the data got in the
!course of previous coupling
!
!rx,ry.rz - the coordinates of the Lagrangian mesh, in au

  real,parameter::RsPerAu=RSun/cAU

!
!vx,vy,vz - three components of the velocity in
!the Lagrangian point, in au/hour

  real,parameter::AuPerHourInv=cSecondPerHour/cAU

!
!bx,by,bz - three components of the magnetic field in 
!the Lagrangian point, in nT

!
!dd - density in nucleons/cm^3

  real,parameter::NucleonPerCm3Inv = 1E-6/cProtonMass

!
!pp - pressure in erg/cm^3

  real, parameter :: ErgPerCm3Inv = 10.0

!      real rx,ry,rz,vx,vy,vz,bx,by,bz,dd,pp
!      integer  iMax,nMax
!      PARAMETER(nMax=2500)
!      common /ihcoord/ rx(nMax),ry(nMax),rz(nMax)
!      common /ihvel/   vx(nMax), vy(nMax), vz(nMax)
!      common /ihmagf/ bx(nMax), by(nMax), bz(nMax)
!      common /ihpdi/     dd(nMax), pp(nMax), iMax
!      real vxOld,vyOld,vzOld
!      common /oldvel/     vxOld(nMax), vyOld(nMax), vzOld(nMax)
!      integer Old_,New_
!      parameter(Old_=1,New_=2)
!      integer iShock,iShockOld
!      common/ishock/iShock,iShockOld
!      real Smooth_VII
!      common/smooth/Smooth_VII(11,nMax,Old_:New_)
!      integer nResolution
!      PARAMETER(nResolution=10)

  integer,parameter::&
       BuffRho_=1,&
       BuffUx_=2, &
       BuffUz_=4, &
       BuffBx_=5, &
       BuffBz_=7, &
       BuffP_=8,  &
       x_  =1,    &
       z_  =3,    &   
       rho_=4,    &
       Ux_ =5,    &
       Uz_ =7,    &
       Bx_ =8,    &
       Bz_ =10,   &
       P_  =11



  integer::iCell
  real:: Weight
  iCell=Put%iCB_II(1,iPutStart)
  if(iCell>nMax)return
  iMax=max(iCell,iMax)
  Weight=W%Weight_I(iPutStart)
  if(DoAdd)then
     Smooth_VII(rho_,iCell,New_)=&
          Smooth_VII(rho_,iCell,New_) + Buff_I(BuffRho_)*&
                                 Weight*NucleonPerCm3Inv

     Smooth_VII(Ux_:Uz_,iCell,New_)=&
          Smooth_VII(Ux_:Uz_,iCell,New_) + Buff_I(BuffUx_:BuffUz_)*&
                                 Weight*AuPerHourInv 

     Smooth_VII(Bx_:Bz_,iCell,New_)=&
         Smooth_VII(Bx_:Bz_,iCell,New_) + Buff_I(BuffBx_:BuffBz_)*&
                                 Weight*1E9    !Tesla to nanotesla

     Smooth_VII(P_,iCell,New_)    =&
          Smooth_VII(P_,iCell,New_) + Buff_I(BuffP_)*&
                                 Weight*ErgPerCm3Inv
  else
     !Save data from previous coupling
     Smooth_VII(:,iCell,Old_)=Smooth_VII(:,iCell,New_)
     !Get the lagrangian mesh coordinates from global vector
     Smooth_VII(x_:z_,iCell,New_)=&
         SP_xyz_i(iCell)*RsPerAu
     
     Smooth_VII(rho_,iCell,New_)=Buff_I(BuffRho_)*&
                                 Weight*NucleonPerCm3Inv

     Smooth_VII(Ux_:Uz_,iCell,New_)=Buff_I(BuffUx_:BuffUz_)*&
                                 Weight*AuPerHourInv 

     Smooth_VII(Bx_:Bz_,iCell,New_)=Buff_I(BuffBx_:BuffBz_)*&
                                 Weight*1E9  !Tesla to nanotesla

     Smooth_VII(P_,iCell,New_)    =Buff_I(BuffP_)*&
                                 Weight*ErgPerCm3Inv
    
  end if
end subroutine SP_put_from_mh
!=============================================================!
subroutine SP_run(tInOut,tLimit)
  use SP_Mod,ONLY:tSimulation,DataInputTime,nStep,nSmooth
  use SP_Mod,ONLY:DoRun,SaveMhData,DoReadMhData
  use ModNumConst
  implicit none
  real,intent(inout)::tInOut
  real,intent(in):: tLimit
  nStep=nStep+1
  if(.not.DoReadMhData)then
     if(nSmooth>0)call SP_smooth_data
     tInOut=max(tInOut,tLimit)
  else
     call SP_read_mh_data
     tInOut=DataInputTime
  end if
  if(SaveMhData)call SP_save_mh_data
  if(DoRun)call SP_sharpen_and_run(tSimulation,DataInputTime)
end subroutine SP_run
!=============================================================!
subroutine SP_smooth_data
  use SP_Mod
  implicit none
  integer::iVar,iSmooth
  integer,parameter::iVarMinSmooth=4,iVarMaxSmooth=11
  logical:: DoSmooth_I(2:iMax-1)
  do iVar=iVarMinSmooth,iVarMaxSmooth
     SMOOTH:do iSmooth=1,nSmooth
        DoSmooth_I=abs(Smooth_VII(iVar,1:iMax-2,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_))>   &
                       2*abs(Smooth_VII(iVar,3:iMax  ,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_)).or.&
                       2*abs(Smooth_VII(iVar,1:iMax-2,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_))<   &
                       abs(Smooth_VII(iVar,3:iMax  ,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_)).or.&
                      (Smooth_VII(iVar,1:iMax-2,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_))*   &
                      (Smooth_VII(iVar,3:iMax  ,New_)-    &
                       Smooth_VII(iVar,2:iMax-1,New_))>cZero
        if(.not.any(DoSmooth_I))CYCLE SMOOTH
        where(DoSmooth_I)
           Smooth_VII(iVar,2:iMax-1,New_)=&
                (Smooth_VII(iVar,1:iMax-2,New_)+    &
                 Smooth_VII(iVar,3:iMax  ,New_))*cHalf
        end where
     end do SMOOTH
  end do
end subroutine SP_smooth_data
!==============================================================
subroutine SP_save_mh_data
  use SP_Mod
  use ModIoUnit
  implicit none
  integer::iFile,iLine
  character(LEN=20)::NameFile
  write(NameFile,'(a,i5.5,a)')'./SP/mh_',nStep,'.dat'
  if(iProc==0)write(iStdOut,*)prefix,'Save file '//NameFile
  iFile=io_unit_new()
  open(iFile,file=NameFile,status='replace')
  write(iFile,'(1es13.5,a)') DataInputTime,' = DataInputTime'
  write(iFile,'(i8,a)') iMax,' = iMax'
  write(iFile,'(a)')&
       'x[AU] y[AU] z[AU] rho[cm-3] ux[AU/hour] uy[AU/hour] '//&
       'uz[AU/hour] Bx[nT] By[nT] Bz[nT] p[erg/cm3]'
  do iLine=1,iMax
     write(iFile,'(100es13.5)')Smooth_VII(:,iLine,New_)
  end do
  close(iFile)
end subroutine SP_save_mh_data
!=============================================================!
subroutine SP_read_mh_data
  use SP_Mod
  use ModIoUnit
  implicit none
  integer::iFile,iLine
  character(LEN=20)::NameFile
  write(NameFile,'(a,i5.5,a)')'./SP/mh_',nStep,'.dat'
  if(iProc==0)write(iStdOut,*)prefix,'Read file '//NameFile
  iFile=io_unit_new()
  open(iFile,FILE=NameFile,STATUS='old',ERR=10)
  read(iFile,*)DataInputTime
  read(iFile,*)iMax
  read(iFile,*)
  do iLine=1,iMax
     Smooth_VII(:,iLine,Old_)=Smooth_VII(:,iLine,New_)
     read(iFile,*)Smooth_VII(:,iLine,New_)
  end do
  close(iFile)
  return
10 call CON_stop(NameFile//' is not available')  
end subroutine SP_read_mh_data
!=============================================================!
subroutine SP_finalize(TimeSimulation)
  implicit none
  real,intent(in)::TimeSimulation
  call SP_closetime
end subroutine SP_finalize
!=============================================================!
subroutine SP_save_restart(TimeSimulation)
  use ModIOUnit,ONLY:io_unit_new
  use SP_Mod
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
  call SP_save_f
  if(.not.SaveMhData)call SP_save_mh_data
  return
10 call CON_stop('SP/restartOUT directory is not available')
end subroutine SP_save_restart
 
