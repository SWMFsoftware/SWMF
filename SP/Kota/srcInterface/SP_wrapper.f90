!=============================================================!
subroutine get_io_unit_new(iNewUnit)
  use ModIoUnit,ONLY:io_unit_new
  integer,intent(out)::iNewUnit
  iNewUnit=io_unit_new()
end subroutine get_io_unit_new
!=============================================================!
module SP_Mod
  use ModConst
  implicit none
  include 'coupler.h'
  include 'stdout.h'
  save
  logical UseSelfSimilarity,UseRefresh
  common/log/UseSelfSimilarity,UseRefresh
  real,parameter::DsResolution=cOne/nResolution
  real::XyzLine_D(3)=(/cZero,cZero,cZero/)
  real::RBoundSC=1.1      !^CMP IF SC
  real::RBoundIH=21.0     !^CMP IF IH
  real::tSimulation,DataInputTime
end module SP_Mod
!=============================================================!
subroutine SP_set_param(CompInfo,TypeAction)
  use CON_comp_info
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
  !-----------------------------------------------------------!
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true.,                        &
          NameVersion=&
          'SEP acceleration with pitch-angle scattering, Kota(2004)', &
          Version    =1.1)
  case('STDOUT')
     iStdOut=STDOUT_
     prefix='SP'
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
        case('#LINE')
           call read_var('XyzLine_D(x_)',XyzLine_D(1))
           call read_var('XyzLine_D(y_)',XyzLine_D(2))
           call read_var('XyzLine_D(z_)',XyzLine_D(3))
           call read_var('RBoundSC'     ,RBoundSC    )!^CMP IF SC 
           call read_var('RBoundIH'     ,RBoundIH    )!^CMP IF IH  
        case default
           call CON_stop(NameSub//&
                ': Unknown command '&
                //NameCommand)
        end select
     end do
  case default
  end select
end subroutine SP_set_param
!============================================================!
subroutine SP_init_session(iSession,TimeIn)
  use SP_Mod,ONLY:tSimulation,UseSelfSimilarity,UseRefresh,iMax
  implicit none
  integer,  intent(in) :: iSession    ! session number (starting from 1)
  real,     intent(in) :: TimeIn      ! seconds from start time
  logical,save::DoInit=.true.
  if(.not.DoInit)return
  DoInit=.false.
  UseSelfSimilarity=.false.
  UseRefresh=.false.
  tSimulation=TimeIn
  call SP_init
  call SP_clean_coupler
  iMax=0
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
  use SP_Mod,ONLY:DataInputTime
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
  use CON_router
  use SP_Mod
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

  real,parameter::NucleonPerCm3Inv=cOne/(cProtonMass*cE6)

!
!pp - pressure in erg/cm^3

  real,parameter::ErgPerCm3Inv=cE1

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
                                 Weight*cE9    !Tesla to nanotesla

     Smooth_VII(P_,iCell,New_)    =&
          Smooth_VII(P_,iCell,New_) + Buff_I(BuffP_)*&
                                 Weight*ErgPerCm3Inv
  else
     !Save data from previous coupling
     Smooth_VII(:,iCell,Old_)=Smooth_VII(:,iCell,New_)
     !Get the lagrangian mesh coordinates from global vector
     Smooth_VII(x_:z_,iCell,New_)=&
          point_state_v('SP_Xyz_DI',3,iCell)*RsPerAu
     
     Smooth_VII(rho_,iCell,New_)=Buff_I(BuffRho_)*&
                                 Weight*NucleonPerCm3Inv

     Smooth_VII(Ux_:Uz_,iCell,New_)=Buff_I(BuffUx_:BuffUz_)*&
                                 Weight*AuPerHourInv 

     Smooth_VII(Bx_:Bz_,iCell,New_)=Buff_I(BuffBx_:BuffBz_)*&
                                 Weight*cE9  !Tesla to nanotesla

     Smooth_VII(P_,iCell,New_)    =Buff_I(BuffP_)*&
                                 Weight*ErgPerCm3Inv
    
  end if
end subroutine SP_put_from_mh
!=============================================================!
subroutine SP_run(tInOut,tLimit)
  use SP_Mod,ONLY:tSimulation,DataInputTime
  implicit none
  real,intent(inout)::tInOut
  real,intent(in):: tLimit
  tInOut=max(tInOut,tLimit)
  call SP_sharpen_and_run(tSimulation,DataInputTime)
end subroutine SP_run
!======================================================================
subroutine SP_finalize(TimeSimulation)
  implicit none
  real,intent(in)::TimeSimulation
  call closetime
end subroutine SP_finalize
!======================================================================
subroutine SP_save_restart(TimeSimulation) 
  implicit none
  real,     intent(in) :: TimeSimulation 
end subroutine SP_save_restart
 
