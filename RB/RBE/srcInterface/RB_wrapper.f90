!^CFG COPYRIGHT UM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!BOI
subroutine RB_set_param(CompInfo, TypeAction)

  !USES:
  use CON_comp_info
  use ModUtilities
  use ModReadParam

  implicit none

  character (len=*), intent(in)     :: TypeAction ! which action to perform
  type(CompInfoType), intent(inout) :: CompInfo   ! component information

  character (len=*), parameter :: NameSub='RB_set_param'
  integer :: iError
  character (len=100) :: NameCommand
  logical             :: UseStrict=.true.

  integer :: iProc, nProc, iComm

  !------------------------------------------------------------------------
  !if(iProc>=0)then
  !   call RB_write_prefix;  
  !   write(iUnitOut,*) NameSub,' TypeAction= ',TypeAction, &
  !        '    iProc=',iProc
  !end if
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                                     &
          Use=.true.,                                       &
          NameVersion='Radiation Belt Environment, M. Fok', &
          Version=1.0)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
     if(nProc>1)call CON_stop('RB_ERROR this version can run on 1 PE only!')
  case('STDOUT')
     !iUnitOut=STDOUT_
     !StringPrefix='RB:'
  case('FILEOUT')
     !call get(CompInfo,iUnitOut=iUnitOut)
     !StringPrefix=''
  case('READ')
     call RB_set_parameters('READ')
     call readinputdata

  case('CHECK')
     ! We should check and correct parameters here
     !if(iProc==0)then
     !   call RB_write_prefix;  write(iUnitOut,*)&
     !        NameSub,': CHECK iSession =',i_session_read()
     !end if
  case('GRID')
     call RB_set_grid
  case default
     call CON_stop(NameSub//' RB_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine RB_set_param
!============================================================================
subroutine RB_set_grid
  use RBE_grid, ONLY: ir, ip
  use RBE_cgrid, ONLY: xlati, phi
  use ModNumConst
  use CON_coupler, ONLY: set_grid_descriptor, is_proc, RB_

  implicit none

  character (len=*), parameter :: NameSub='RB_set_grid'
  integer :: iSize,jSize
  real, dimension (:,:), allocatable :: gridLat,gridLT
  real :: Radius_I(1)
  logical :: IsInitialized=.false.
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)'RB_set_grid_descriptor called, IsInitialized=',&
       IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  Radius_I(1) = (6375.0+120.0)*1000.0 ! radial size of the ionosphere in meters

  ! RB grid size in generalized coordinates
  call set_grid_descriptor( RB_,                 & ! component index
       nDim=2,                                   & ! dimensionality
       nRootBlock_D=(/1,1/),                     & ! single block
       nCell_D=(/ir+2, ip/),                     & ! size of cell based grid
       XyzMin_D=(/cHalf, cHalf/),                & ! min gen.coords for cells
       XyzMax_D=(/ir+1.5,ip-0.5/),               & ! max gen.coords for cells
       TypeCoord='SMG',                          & ! solar magnetic coord
       Coord1_I=cRadToDeg*xlati,                 & ! latitude in degrees
       Coord2_I=cRadToDeg*phi,                   & ! longitude in degrees
       Coord3_I=Radius_I,                        & ! radial size in meters
       IsPeriodic_D=(/.false.,.true./))            ! periodic in longitude

  if(DoTest)then
     write(*,*)NameSub,' ir,ip=',ir,ip
     write(*,*)NameSub,' size(xlati)=',size(xlati),' size(phi)=',size(phi)
     write(*,*)NameSub,' xlati=',xlati
     write(*,*)NameSub,' phi=',phi
  end if

end subroutine RB_set_grid
!==============================================================================

subroutine RB_init_session(iSession, TimeSimulation)

  implicit none

  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='RB_init_session'
  !------------------------------------------------------------------------
  call rbe_init

end subroutine RB_init_session
!==============================================================================

subroutine RB_run(TimeSimulation,TimeSimulationLimit)

  use rbe_time, ONLY: t, dt
  use rbe_cread2, ONLY: dtmax

  implicit none

  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
  real, intent(inout) :: TimeSimulation   ! current time of component
  
  character(len=*), parameter :: NameSub='RB_run'

  !------------------------------------------------------------------------

  dt = min(dtmax, 0.5*(TimeSimulationLimit - TimeSimulation))
  call rbe_run

  ! return time at the end of the time step to CON
  TimeSimulation   = t
  
end subroutine RB_run
!===========================================================================

subroutine RB_finalize(TimeSimulation)

  !USES:
  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='RB_finalize'

  !-------------------------------------------------------------------------

  !call RB_write_prefix; write(iUnitOut,*) &
  !     NameSub,' at TimeSimulation=',TimeSimulation

end subroutine RB_finalize
!===========================================================================

subroutine RB_save_restart(TimeSimulation)

  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='RB_save_restart'

  !-------------------------------------------------------------------------
  call rbe_save_restart

end subroutine RB_save_restart
!===========================================================================

subroutine RB_put_from_gm(Integral_IIV,iSizeIn,jSizeIn,nIntegralIn,&
            BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)

!  use ModIoUnit, ONLY: UnitTmp_
!  use ModRiceRB
  
  implicit none

  integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
  real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
  integer, intent(in) :: nVarLine, nPointLine
  real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)

  character (len=*),intent(in) :: NameVar
  real, intent(in) :: tSimulation

  real, parameter :: noValue=-99999.

  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='RB_put_from_gm'
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(nVarLine /= 4) then
     write(*,*)'nVarLine=',nVarLine
     call CON_stop(NameSub//' invalid nVarLine (should be 4)')
  end if

  if(DoTestMe)then
     write(*,*)NameSub,' iSizeIn,jSizeIn,nIntegralIn=',&
          iSizeIn,jSizeIn,nIntegralIn
     write(*,*)NameSub,' nVarLine,nPointLine=',nVarLine,nPointLine
     write(*,*)NameSub,' Integral_IIV(1,23,:)=',Integral_IIV(1,23,:)
     write(*,*)NameSub,' BufferLine_VI(:,1) =',BufferLine_VI(:,1)
     write(*,*)NameSub,' BufferLine_VI(:,2) =',BufferLine_VI(:,2)
     write(*,*)NameSub,' IMF=',Integral_IIV(1:8,1,4)
  end if

  !n = 0
  !do iLat = 1, nLat
  !   do iLon = 1, nLon
  !      n = n+1
  !      
  !   end do
  !end do

end subroutine RB_put_from_gm
!============================================================================
