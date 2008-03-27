!^CFG COPYRIGHT UM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!BOI
subroutine IM_set_param(CompInfo, TypeAction)

  !USES:
  use CON_comp_info
  use ModUtilities
  use ModReadParam
  use CON_coupler, ONLY: Couple_CC, GM_, IM_, IE_
  use rbe_cread2,  ONLY: UseGm,UseIE
  implicit none

  character (len=*), intent(in)     :: TypeAction ! which action to perform
  type(CompInfoType), intent(inout) :: CompInfo   ! component information

  character (len=*), parameter :: NameSub='IM_set_param'
  integer :: iError
  character (len=100) :: NameCommand
  logical             :: UseStrict=.true.

  integer :: iProc, nProc, iComm

  !------------------------------------------------------------------------

  UseGm = Couple_CC(GM_, IM_) % DoThis
  UseIE = Couple_CC(IE_, IM_) % DoThis
  
  !if(iProc>=0)then
  !   call IM_write_prefix;  
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
     if(nProc>1)call CON_stop('IM_ERROR this version can run on 1 PE only!')
  case('STDOUT')
     !iUnitOut=STDOUT_
     !StringPrefix='IM:'
  case('FILEOUT')
     !call get(CompInfo,iUnitOut=iUnitOut)
     !StringPrefix=''
  case('READ')
     call IM_set_parameters('READ')
     call readinputdata

  case('CHECK')
     ! We should check and correct parameters here
     !if(iProc==0)then
     !   call IM_write_prefix;  write(iUnitOut,*)&
     !        NameSub,': CHECK iSession =',i_session_read()
     !end if
  case('GRID')
     call IM_set_grid
  case default
     call CON_stop(NameSub//' IM_ERROR: invalid TypeAction='//TypeAction)
  end select

end subroutine IM_set_param
!============================================================================
subroutine IM_set_grid
  use IME_grid, ONLY: ir, ip
  use IME_cgrid, ONLY: xlati, phi
  use ModNumConst
  use CON_coupler, ONLY: set_grid_descriptor, is_proc, IM_

  implicit none

  character (len=*), parameter :: NameSub='IM_set_grid'
  integer :: iSize,jSize
  real, dimension (:,:), allocatable :: gridLat,gridLT
  real :: Radius_I(1)
  logical :: IsInitialized=.false.
  logical :: DoTest, DoTestMe
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTest)write(*,*)'IM_set_grid_descriptor called, IsInitialized=',&
       IsInitialized
  if(IsInitialized) return
  IsInitialized=.true.

  Radius_I(1) = (6375.0+120.0)*1000.0 ! radial size of the ionosphere in meters

  ! IM grid size in generalized coordinates
  call set_grid_descriptor( IM_,                 & ! component index
       nDim=2,                                   & ! dimensionality
       nRootBlock_D=(/1,1/),                     & ! single block
       nCell_D=(/ir+2, ip/),                     & ! size of cell based grid
       XyzMin_D=(/cHalf, cHalf/),                & ! min gen.coords for cells
       XyzMax_D=(/ir+2.5,ip-0.5/),               & ! max gen.coords for cells
       TypeCoord='SMG',                          & ! solar magnetic coord
       Coord1_I=cRadToDeg*xlati(0:ir+1),         & ! latitude in degrees
       Coord2_I=mod(cRadToDeg*phi+180.0,360.0),  & ! longitude in degrees
       Coord3_I=Radius_I,                        & ! radial size in meters
       IsPeriodic_D=(/.false.,.true./))            ! periodic in longitude

  if(DoTest)then
     write(*,*)NameSub,' ir,ip=',ir,ip
     write(*,*)NameSub,' size(xlati)=',size(xlati),' size(phi)=',size(phi)
     write(*,*)NameSub,' xlati=',xlati
     write(*,*)NameSub,' phi=',phi
  end if

end subroutine IM_set_grid
!==============================================================================

subroutine IM_init_session(iSession, TimeSimulation)
  use rbe_constant
  use rbe_cread2
  use rbe_cgrid

  implicit none

  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='IM_init_session'
  !------------------------------------------------------------------------
  ! GM info needed before initialization just set up latitude/longitude grid

  call grids(re,rc,xme,xmp,q,c,js)
  
end subroutine IM_init_session
!==============================================================================

subroutine IM_run(TimeSimulation,TimeSimulationLimit)

  use rbe_time, ONLY: t, dt
  use rbe_cread2, ONLY: dtmax

  implicit none

  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
  real, intent(inout) :: TimeSimulation   ! current time of component
  
  Logical, save :: IsInitiallized = .false.
  character(len=*), parameter :: NameSub='IM_run'

  !------------------------------------------------------------------------
  
  if (.not. IsInitiallized) then
     call rbe_init
     IsInitiallized = .true.
  endif

  dt = min(dtmax, 0.5*(TimeSimulationLimit - TimeSimulation))
  call rbe_run

  ! return time at the end of the time step to CON
  TimeSimulation   = t
  
end subroutine IM_run
!===========================================================================

subroutine IM_finalize(TimeSimulation)

  !USES:
  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='IM_finalize'

  !-------------------------------------------------------------------------

  !call IM_write_prefix; write(iUnitOut,*) &
  !     NameSub,' at TimeSimulation=',TimeSimulation

end subroutine IM_finalize
!===========================================================================

subroutine IM_save_restart(TimeSimulation)

  implicit none

  real,     intent(in) :: TimeSimulation   ! seconds from start time
  character(len=*), parameter :: NameSub='IM_save_restart'

  !-------------------------------------------------------------------------
  call rbe_save_restart

end subroutine IM_save_restart
!===========================================================================

subroutine IM_put_from_gm(Integral_IIV,iSizeIn,jSizeIn,nIntegralIn,&
            BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)
  use ModGmRb
  use rbe_grid,    ONLY: nLat => ir, nLon => ip
  use rbe_constant,ONLY: rEarth => re
  use rbe_cread2,  ONLY: xnswa,vswa,bxw,byw,bzw,nsw,iyear,iday
  implicit none

  integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
  real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
  integer, intent(in) :: nVarLine, nPointLine
  real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)

  character (len=*),intent(in) :: NameVar
  real, intent(in) :: tSimulation

  real, parameter :: noValue=-99999.
  integer :: n,iLat,iLon
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='IM_put_from_gm'
  logical,save :: IsFirstCall = .true.
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
     !write(*,*)NameSub,' Integral_IIV(21,1,:)=',Integral_IIV(21,1,:)
     write(*,*)NameSub,' BufferLine_VI(:,1) =',BufferLine_VI(:,1)
     write(*,*)NameSub,' BufferLine_VI(:,2) =',BufferLine_VI(:,2)
     write(*,*)NameSub,' IMF=',Integral_IIV(1:8,1,4)
  end if
  
  if (allocated(StateLine_VI)) then
     deallocate(StateLine_VI,StateIntegral_IIV)
  endif
  
  if (.not.allocated(StateLine_VI)) then
     allocate(StateLine_VI(nVarLine,nPointLine),&
          StateIntegral_IIV(0:iSizeIn-1,jSizeIn,nIntegralIn))
  endif
  
  StateLine_VI      = BufferLine_VI
  StateIntegral_IIV = Integral_IIV
  nPoint=nPointLine
  
  !Convert Units
  StateLine_VI(2,:) = StateLine_VI(2,:) / rEarth ! m --> Earth Radii
  StateLine_VI(3,:) = StateLine_VI(3,:) / rEarth ! m --> Earth Radii

  !Solar wind values
  xnswa(1) = Integral_IIV(1,1,4)*1.0e-6                   !m^-3 -->/cc
  vswa (1) = sqrt(sum(Integral_IIV(2:4,1,4)**2.0))*1.0e-3 !m/s-->km/s
  
  bxw(1) = Integral_IIV(5,1,4)*1.0e9      !T --> nT
  byw(1) = Integral_IIV(6,1,4)*1.0e9      !T --> nT
  bzw(1) = Integral_IIV(7,1,4)*1.0e9      !T --> nT

  nsw = 1
  
  iyear=2002
  iday=1
    

  ! create an index array on the first call
  if (IsFirstCall) then
     n = 0
     do iLon = 1, nLon
        do iLat = 0, nLat+1
           n = n+1
           iLineIndex_II(iLon,iLat) = n
        end do
     end do
     IsFirstCall = .false.
  endif

end subroutine IM_put_from_gm
!============================================================================

subroutine IM_put_from_ie(Buffer_IIV, iSize, jSize, nVarIn, &
                 Name_V, iBlock)
  
  use rbe_grid,    ONLY: nLat => ir, nLon => ip
  use rbe_cgrid,   ONLY: Lat_I => xlati, Lon_I => phi
  use rbe_convect, ONLY: Potential_II => potent
  use CON_coupler, ONLY: Grid_C, IE_
  use ModNumConst, ONLY: cTwoPi,cPi,cHalfPi
  use ModInterpolate, ONLY: bilinear
  implicit none

  character(len=*), parameter :: NameSub='IM_put_from_ie'

  !INPUT ARGUMENTS:
  integer, intent(in):: iSize, jSize, nVarIn, iBlock
  real, intent(in) :: Buffer_IIV(iSize, jSize, nVarIn)
  character(len=*), intent(in) :: Name_V(nVarIn)

  integer, parameter :: nVar = 2
  integer, parameter :: South_ = 1, North_ = 2

  logical :: IsPotFound, IsJrFound
  real    :: dPhiIono,dThetaIono
  integer :: iLat, iLon, iVar, nThetaIono, nPhiIono
  real,dimension(:,:),allocatable  :: x,y,z
  Character(len=100) :: NameElectrodynamics
  !----------------------------------------------------------------------------
  if(iBlock /= north_) RETURN
  
  nThetaIono = Grid_C(IE_) % nCoord_D(1)
  nPhiIono   = Grid_C(IE_) % nCoord_D(2)
  if(nThetaIono /= 2*iSize - 1 .or. nPhiIono /= jSize)then
     write(*,*)NameSub,': Grid_C(IE_)%nCoord_D(1:2)=',&
          Grid_C(IE_) % nCoord_D(1:2)
     write(*,*)NameSub,': iSize,2*iSize-1,jSize=',iSize,2*iSize-1,jSize
     call CON_stop(NameSub//' ERROR: Inconsistent IE grid sizes')
  endif
  ! Get ionospheric grid spacing for use in interpolation
  dPhiIono   = cTwoPi / (nPhiIono-1)
  dThetaIono = maxval( Grid_C(IE_) % Coord1_I(:) ) / (nThetaIono -1)

  IsPotFound = .false.
  do iVar = 1, nVarIn
     select case(Name_V(iVar))
     case('Pot')
        IsPotFound = .true.
        do iLon=1,nLon
           do iLat=0,nLat+1
              !Interpolate IE potential onto IM grid
              ! Note that the IM grid is 180 degrees rotated relative to 
              ! the usual SM coordinates used by IE
              Potential_II(iLat,iLon) = &
                   bilinear (Buffer_IIV(:, :, iVar),1,iSize,1,jSize,&
                   (/ (Lat_I(iLat))/dThetaIono+1,&
                      modulo(Lon_I(iLon)+cPi,cTwoPi)/dPhiIono+1 /) )
           end do
        end do
     end select
  end do
  

  if(.not.IsPotFound)then
     write(*,*)NameSub,': Name_V=',Name_V
     call CON_stop(NameSub//' could not find Pot')
  end if
  
end subroutine IM_put_from_ie
!==============================================================================
