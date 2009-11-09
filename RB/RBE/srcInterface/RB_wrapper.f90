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
  use CON_coupler, ONLY: Couple_CC, GM_, RB_, IE_
  use rbe_cread2,  ONLY: UseGm,UseIE
  implicit none

  character (len=*), intent(in)     :: TypeAction ! which action to perform
  type(CompInfoType), intent(inout) :: CompInfo   ! component information

  character (len=*), parameter :: NameSub='RB_set_param'
  integer :: iError
  character (len=100) :: NameCommand
  logical             :: UseStrict=.true.

  integer :: iProc, nProc, iComm

  !------------------------------------------------------------------------
  UseGm = Couple_CC(GM_, RB_) % DoThis
  UseIE = Couple_CC(IE_, RB_) % DoThis
  
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
       nCell_D=(/ir, ip/),                     & ! size of cell based grid
       XyzMin_D=(/cHalf, cHalf/),                & ! min gen.coords for cells
       XyzMax_D=(/ir-0.5,ip-0.5/),               & ! max gen.coords for cells
       TypeCoord='SMG',                          & ! solar magnetic coord
       Coord1_I=cRadToDeg*xlati(1:ir),         & ! latitude in degrees
       Coord2_I=mod(cRadToDeg*phi+180.0,360.0),  & ! longitude in degrees
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
  use rbe_constant
  use rbe_cread2
  use rbe_cgrid

  implicit none

  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='RB_init_session'
  !------------------------------------------------------------------------
  ! GM info needed before initialization just set up latitude/longitude grid

  call grids(re,rc,xme,xmp,q,c,js)
  
end subroutine RB_init_session
!==============================================================================

subroutine RB_run(TimeSimulation,TimeSimulationLimit)

  use rbe_time, ONLY: t, dt
  use rbe_cread2, ONLY: dtmax

  implicit none

  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded
  real, intent(inout) :: TimeSimulation   ! current time of component
  
  Logical, save :: IsInitiallized = .false.
  character(len=*), parameter :: NameSub='RB_run'

  !------------------------------------------------------------------------
  
  if (.not. IsInitiallized) then
     call rbe_init
     IsInitiallized = .true.
  endif

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
  call rbe_save_result(.true., .false.)

end subroutine RB_save_restart
!===========================================================================

subroutine RB_put_from_gm(Integral_IIV,iSizeIn,jSizeIn,nIntegralIn,&
            BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)
  use ModGmRb
  use rbe_grid,    ONLY: nLat => ir, nLon => ip
  use rbe_constant,ONLY: rEarth => re
  use rbe_cread2,  ONLY: xnswa,vswa,bxw,byw,bzw,nsw,iyear,iday,UseSmooth
  use ModPrerunField,ONLY: DoWritePrerun, save_prerun
  implicit none

  integer, intent(in) :: iSizeIn, jSizeIn, nIntegralIn
  real,    intent(in) :: Integral_IIV(iSizeIn,jSizeIn,nIntegralIn)
  integer, intent(in) :: nVarLine, nPointLine
  real,    intent(in) :: BufferLine_VI(nVarLine, nPointLine)

  character (len=*),intent(in) :: NameVar
  real, intent(in) :: tSimulation

  real, parameter :: noValue=-99999.
  real :: SwDensMax, SwVelMax, SwDensMin, SwVelMin
  integer :: n,iLat,iLon
  logical :: DoTest, DoTestMe
  character(len=*), parameter :: NameSub='RB_put_from_gm'
  logical,save :: IsFirstCall = .true.
  !-------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  if(NameVar /= 'Z0x:Z0y:Z0b:I_I:S_I:R_I:B_I:IMF') &
       call CON_stop(NameSub//' invalid NameVar='//NameVar)

  if(nVarLine /= nVar) then
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
     write(*,*)NameSub,' IMF: Density  = ',Integral_IIV(1,1,6)
     write(*,*)NameSub,' IMF: Velocity = ',Integral_IIV(2,1,6)
     write(*,*)NameSub,' IMF: Bx       = ',Integral_IIV(5,1,6)
     write(*,*)NameSub,' IMF: By       = ',Integral_IIV(6,1,6)
     write(*,*)NameSub,' IMF: Bz       = ',Integral_IIV(7,1,6)
  end if
  
  if (allocated(StateLine_VI)) then
     deallocate(StateLine_VI,StateIntegral_IIV)
  endif
  
  if (.not.allocated(StateLine_VI)) then
     allocate(StateLine_VI(nVarLine,nPointLine),&
          StateIntegral_IIV(iSizeIn,jSizeIn,nIntegralIn))
  endif
  
  StateLine_VI      = BufferLine_VI
  StateIntegral_IIV = Integral_IIV
  nPoint    = nPointLine
  nIntegral = nIntegralIn
  !Convert Units
  StateLine_VI(2,:) = StateLine_VI(2,:) / rEarth ! m --> Earth Radii
  StateLine_VI(3,:) = StateLine_VI(3,:) / rEarth ! m --> Earth Radii

  !Solar wind values
  if(IsFirstCall .or. (.not. UseSmooth)) then
     xnswa(1) = Integral_IIV(1,1,6)*1.0e-6                   !m^-3 -->/cc
     vswa (1) = sqrt(sum(Integral_IIV(2:4,1,6)**2.0))*1.0e-3 !m/s-->km/s
  else
     ! Update Solar wind value, but do not let them change more than 5 percent 
     ! per update
     SwDensMax = 1.05*xnswa(1)
     SwDensMin = 0.95*xnswa(1)
     SwVelMax  = 1.05*vswa(1)
     SwVelMin  = 0.95*vswa(1)
     xnswa(1) = min(SwDensMax,Integral_IIV(1,1,6)*1.0e-6)
     xnswa(1) = max(SwDensMin,xnswa(1))
     vswa(1)  = min(SwVelMax,sqrt(sum(Integral_IIV(2:4,1,6)**2.0))*1.0e-3)
     vswa(1)  = max(SwVelMin,vswa(1))
  endif
  bxw(1) = Integral_IIV(5,1,6)*1.0e9      !T --> nT
  byw(1) = Integral_IIV(6,1,6)*1.0e9      !T --> nT
  bzw(1) = Integral_IIV(7,1,6)*1.0e9      !T --> nT

  nsw = 1
  
  iyear=2002
  iday=1
    

  ! create an index array on the first call
  if (IsFirstCall) then
     n = 0
     do iLon = 1, nLon
        do iLat = 1, nLat
           n = n+1
           iLineIndex_II(iLon,iLat) = n
        end do
     end do
     IsFirstCall = .false.
  endif

  if (DoWritePrerun) call save_prerun(tSimulation)
end subroutine RB_put_from_gm
!============================================================================

subroutine RB_put_from_ie(Buffer_IIV, iSize, jSize, nVarIn, &
                 Name_V, iBlock)
  
  use rbe_grid,    ONLY: nLat => ir, nLon => ip
  use rbe_cgrid,   ONLY: Lat_I => xlati, Lon_I => phi
  use rbe_convect, ONLY: Potential_II => potent
  use CON_coupler, ONLY: Grid_C, IE_
  use ModNumConst, ONLY: cTwoPi,cPi,cHalfPi
  use ModInterpolate, ONLY: bilinear
  use ModPrerunField,ONLY: DoWritePrerun, save_prerun_IE
  use rbe_time,    ONLY: tSimulation => t
  implicit none

  character(len=*), parameter :: NameSub='RB_put_from_ie'

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
           do iLat=1,nLat
              !Interpolate IE potential onto RB grid
              ! Note that the RB grid is 180 degrees rotated relative to 
              ! the usual SM coordinates used by IE
              Potential_II(iLat,iLon) = &
                   bilinear (Buffer_IIV(:, :, iVar),1,iSize,1,jSize,&
                   (/ (Lat_I(iLat))/dThetaIono+1,&
                      modulo(Lon_I(iLon)+cPi,cTwoPi)/dPhiIono+1 /) )
           end do
        end do
     end select
  end do
  
  if (DoWritePrerun) call save_prerun_IE(tSimulation)
  if(.not.IsPotFound)then
     write(*,*)NameSub,': Name_V=',Name_V
     call CON_stop(NameSub//' could not find Pot')
  end if
  
end subroutine RB_put_from_ie
!==============================================================================

subroutine RB_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
  ! Puts satellite locations and names from GM into RB.
  use ModRbSat, ONLY: nRbSats, DoWriteSats, NameSat_I, SatLoc_3I
  use ModNumConst,   ONLY: cDegToRad
  
  implicit none
  character (len=*),parameter :: NameSub='RB_put_sat_from_gm'

  ! Arguments
  integer, intent(in)            :: nSats
  real, intent(in)               :: Buffer_III(4,2,nSats)
  character(len=100), intent(in) :: Buffer_I(nSats)

  ! Internal variables
  integer :: iError, iSat, l1, l2
  !--------------------------------------------------------------------------- 
  ! Activate satellite writing in RCM
  DoWriteSats = .true.
  nRbSats = nSats

  ! Check allocation of sat tracing variables
  if(allocated(SatLoc_3I)) deallocate(SatLoc_3I)
  if(allocated(NameSat_I)) deallocate(NameSat_I)

  allocate(SatLoc_3I(4,2,nRbSats), stat=iError)
  allocate(NameSat_I(nRbSats),     stat=iError)

  ! Assign incoming values, remove path and extension from name.
  SatLoc_3I = Buffer_III
  SatLoc_3I(4,2,:)=SatLoc_3I(4,2,:)
  do iSat=1, nSats
     l1 = index(Buffer_I(iSat), '/', back=.true.) + 1
     l2 = index(Buffer_I(iSat), '.') - 1
     if (l1-1<=0) l1=1
     if (l2+1<=0) l2=len_trim(Buffer_I(iSat))
     NameSat_I(iSat) = Buffer_I(iSat)(l1:l2)
  end do

  ! Change to correct units (degrees to radians)
!  SatLoc_3I(1,2,:) = (90. - SatLoc_3I(1,2,:)) * cDegToRad
!  SatLoc_3I(2,2,:) =        SatLoc_3I(2,2,:)  * cDegToRad

end subroutine RB_put_sat_from_gm

!==============================================================================
