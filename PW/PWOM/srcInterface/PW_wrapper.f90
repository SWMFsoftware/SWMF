!^CFG COPYRIGHT UM
! Wrapper for the empty PWOM (PW) component
!==========================================================================
subroutine PW_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use CON_coupler
  use ModIoUnit, only: STDOUT_
  use ModPWOM, only: iUnitOut, iProc, nProc, iComm, StringPrefix, &
       nTotalLine, nAlt
  use ModCommonVariables, only: Altd

  implicit none

  character (len=*), parameter :: NameSub='PW_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do

  integer :: i
  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.true., &
          NameVersion='PWOM (A. Glocer et al.)', &
          Version    =1.0)
  case('MPI')
     call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)

  case('READ')
     call PW_set_parameters('READ')

  case('CHECK')
     ! call PW_set_parameters('CHECK')

  case('STDOUT')

     iUnitOut=STDOUT_
     if(nProc==1)then
        StringPrefix='PW:'
     else
        write(StringPrefix,'(a,i3.3,a)')'PW',iProc,':'
     end if

  case('FILEOUT')

     call get(CompInfo,iUnitOut=iUnitOut)
     StringPrefix=''

  case('GRID')
     ! We pretend to have an nRadius*nLine 2D grid, 
     ! because it is unstructured in theta,phi
     call set_grid_descriptor( &
          PW_,                             &! component index
          nDim=2,                          &! dimensionality
          nRootBlock_D=(/1,nProc/),        &! distributed in the second dimension
          nCell_D =(/ nAlt, nTotalLine /),          &! size of the grid
          XyzMin_D=(/Altd(1), 1.0/),                &! min altitude and index
          XyzMax_D=(/Altd(nAlt), real(nTotalLine)/),&! max altitude and index
          Coord1_I = Altd(1:nAlt),                  &! altitudes
          Coord2_I= (/ (real(i), i=1,nTotalLine) /),   &! indexes
          TypeCoord='SMG')                           ! 

  case default
     call CON_stop(NameSub//': PW_ERROR: empty version cannot be used!')
  end select

end subroutine PW_set_param

!==============================================================================

subroutine PW_init_session(iSession, TimeSimulation)
  use ModPWOM, ONLY: iProc, UseIE, Time
  use ModPwTime
  use ModTimeConvert, ONLY: time_real_to_int
  use CON_coupler, ONLY: Couple_CC, IE_, PW_
  use CON_physics, ONLY: get_time
  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  
  character(len=*), parameter :: NameSub='PW_init_session'

  logical :: DoInitialize = .true.
  !----------------------------------------------------------------------------
  UseIE = Couple_CC(IE_, PW_) % DoThis
  
  if(DoInitialize) then
     ! Initialize PW
     call PW_initialize
     ! Set time related variables for UA 
     call get_time(tStartOut = StartTime)
     CurrentTime = StartTime + TimeSimulation
     call time_real_to_int(StartTime, iStartTime)
     if(Time /= TimeSimulation .and. iProc==0) &
          write(*,*)'WARNING ',NameSub,': PWOM Time=',Time, &
          ' differs from SWMF TimeSimulation=',TimeSimulation 
     Time = TimeSimulation
  endif
  DoInitialize = .false.

end subroutine PW_init_session

!==============================================================================

subroutine PW_finalize(TimeSimulation)

  use ModPWOM, ONLY: iLine, nLine, iUnitGraphics, iUnitOutput,nLog,iLineGlobal,&
       r_C, State_CVI, GeoMagLat_I,GeoMagLon_I,     &
       ThetaLine_I, PhiLine_I, xLine_I, yLine_I, zLine_I, &
       xLineOld_I, yLineOld_I, zLineOld_I, UthetaLine_I,  &
       UphiLine_I, UxLine_I, UyLine_I, UzLine_I,          &
       OmegaLine_I, JrLine_I, iThetaLine_I,iPhiLine_I,    &
       NameRestartIn, NameRestart, NameGraphics,          &
       NameOutput,  iUnitRestart, iUnitRestartIn

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_finalize'
  !-------------------------------------------------------------------------
  if (nLog == -1) then
     do iLine=1,nLine
        close(iUnitOutput(iLine))
     enddo
  elseif(nLog ==0) then
     !do nothing in this case
  elseif(nLog==iLineGlobal(iLine)) then
     close(iUnitOutput(iLine))
  else
  end if

  ! Deallocate variables needed for simulation 
  deallocate(r_C, State_CVI, GeoMagLat_I,GeoMagLon_I,     &
       ThetaLine_I, PhiLine_I, xLine_I, yLine_I, zLine_I, &
       xLineOld_I, yLineOld_I, zLineOld_I, UthetaLine_I,  &
       UphiLine_I, UxLine_I, UyLine_I, UzLine_I,          &
       OmegaLine_I, JrLine_I, iThetaLine_I,iPhiLine_I,    &
       NameRestartIn, NameRestart, NameGraphics,          &
       NameOutput,  iUnitRestart, iUnitRestartIn,         &
       iUnitGraphics,iUnitOutput, iLineGlobal)


end subroutine PW_finalize

!==============================================================================

subroutine PW_save_restart(TimeSimulation)
  
  Use ModPWOM, only: &
       nAlt,r_C,GeoMagLat_I,GeoMagLon_I,DtVertical,&
       nStep,NameRestart, &
       State_CVI,nLine
  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  integer :: iLine
  character(len=*), parameter :: NameSub='PW_save_restart'
  !---------------------------------------------------------------------------
  
  do iLine=1,nLine
     call PW_write_restart(&
      nAlt,r_C,GeoMagLat_I(iLine),GeoMagLon_I(iLine),TimeSimulation,DtVertical,&
      nStep,NameRestart(iLine), &
      State_CVI(:,:,iLine))
  enddo



end subroutine PW_save_restart

!==============================================================================

subroutine PW_run(TimeSimulation,TimeSimulationLimit)

  use ModPWOM, ONLY: iLine, nLine, Time, nStep, DtHorizontalOrig, &
       DtHorizontal, DtOutput, DoPlotElectrodynamics, DtPlotElectrodynamics, &
       Tmax

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)
  if(DoTestMe)write(*,*) NameSub,': Time,TimeSimulation,TimeSimulationLimit=',&
       Time,TimeSimulation,TimeSimulationLimit

  DtHorizontal = min(DtHorizontalOrig, TimeSimulationLimit - Time)

  Tmax = TimeSimulationLimit

  if(DoTestMe)write(*,*) NameSub,': DtHorizontalOrig,DtHorizontal=',&
       DtHorizontalOrig,DtHorizontal

  ! Avoid taking tiny little steps
  if(DtHorizontal < 1e-6)then
     Time           = TimeSimulationLimit
     TimeSimulation = TimeSimulationLimit
     RETURN
  end if


  do iLine=1,nLine
     call move_line
     call PW_advance_line
  end do
  !Output the electrodynamics info
  if (DoPlotElectrodynamics) then
     if (floor(Time/DtPlotElectrodynamics) &
          /= floor((Time-DtHorizontal)/DtPlotElectrodynamics) ) &
          call PW_print_electrodynamics
  endif
  TimeSimulation = Time

  if(DoTestMe)write(*,*) NameSub,': final nStep, TimeSimulation=',&
       nStep, TimeSimulation

end subroutine PW_run

!==============================================================================

subroutine PW_put_from_ie(Buffer_IIV, iSize, jSize, nVarIn, &
                 Name_V, iBlock)

  use ModPWOM, ONLY: allocate_ie_variables, Phi_G, Theta_G, Potential_G, Jr_G
  use CON_coupler, ONLY: Grid_C, IE_
  implicit none

  character(len=*), parameter :: NameSub='PW_put_from_ie'

  !INPUT ARGUMENTS:
  integer, intent(in):: iSize, jSize, nVarIn, iBlock
  real, intent(in) :: Buffer_IIV(iSize, jSize, nVarIn)
  character(len=*), intent(in) :: Name_V(nVarIn)

  integer, parameter :: nVar = 2
  integer, parameter :: South_ = 1, North_ = 2

  logical :: IsPotFound, IsJrFound

  integer :: i, j, iVar, nThetaIono, nPhiIono
  !----------------------------------------------------------------------------
  if(iBlock /= north_) RETURN

  if(.not.allocated(Phi_G))then
     nThetaIono = Grid_C(IE_) % nCoord_D(1)
     nPhiIono   = Grid_C(IE_) % nCoord_D(2)
     if(nThetaIono /= 2*iSize - 1 .or. nPhiIono /= jSize)then
        write(*,*)NameSub,': Grid_C(IE_)%nCoord_D(1:2)=',&
             Grid_C(IE_) % nCoord_D(1:2)
        write(*,*)NameSub,': iSize,2*iSize-1,jSize=',iSize,2*iSize-1,jSize
        call CON_stop(NameSub//' ERROR: Inconsistent IE grid sizes')
     endif

     call allocate_ie_variables(jSize, iSize)

     do i = 1, iSize
        Theta_G(1:jSize,i) = Grid_C(IE_) % Coord1_I(i)
     end do
     do j = 1, jSize
        Phi_G( j,1:iSize) = Grid_C(IE_) % Coord2_I(j)
     end do

     ! This has to be initialized (will be replaced below)
     Potential_G = 0.0 
     call PW_get_electrodynamics
     call initial_line_location
  end if

  IsPotFound = .false.
  IsJrFound  = .false.
  do iVar = 1, nVarIn
     select case(Name_V(iVar))
     case('Pot')
        IsPotFound = .true.
        do i=1,iSize
           do j=1,jSize
              Potential_G(j,i) = Buffer_IIV(iSize+1-i, j, iVar)
           end do
        end do
        case('Jr')
        IsJrFound = .true.
        do i=1,iSize
           do j=1,jSize
              Jr_G(j,i) = Buffer_IIV(iSize+1-i, j, iVar)
           end do
        end do
     end select

  end do

  if(.not.IsPotFound .or. .not.IsJrFound)then
     write(*,*)NameSub,': Name_V=',Name_V
     call CON_stop(NameSub//' could not find Pot or Jr')
  end if

  call PW_get_electrodynamics

end subroutine PW_put_from_ie
!==============================================================================

subroutine PW_get_for_gm(Buffer_VI, nVar, nLineTotal, Name_V, tSimulation)

  use ModPWOM, only : iComm,nProc,&
                      ThetaLine_I,PhiLine_I, &
                      State_CVI,&
                      nLine,nAlt,&
                      nLine_P, nLineBefore_P
  use ModCommonPlanet, only: RhoO_,RhoH_,RhoHe_,uO_,uH_,uHe_
  use ModMpi

  implicit none
  character (len=*),parameter :: NameSub='PW_get_for_gm'

  integer, intent(in)           :: nVar, nLineTotal
  real, intent(out)             :: Buffer_VI(nVar, nLineTotal)
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar,i,iError
  integer :: iSendCount

  integer, allocatable :: iDisplacement_P(:), iRecieveCount_P(:)
  real,    allocatable :: SendBuffer_VI(:,:)

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! The sizes and displacements of MPI messages
  allocate(iRecieveCount_P(nProc), iDisplacement_P(nProc))
  allocate(SendBuffer_VI(nVar, nLine))
  iRecieveCount_P = nVar*nLine_P
  iDisplacement_P = nVar*nLineBefore_P

  ! Fill buffer using SI units for sending on each proc
  do iVar=1,nVar
     select case (Name_V(iVar))
     case('CoLat')
        SendBuffer_VI(iVar,:)=ThetaLine_I(1:nLine)
     case('Longitude')
        SendBuffer_VI(iVar,:)=PhiLine_I(1:nLine)
     case('Density1')
        ! g/cm^3 = 1000 * kg/m^3
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,RhoH_ ,1:nLine)*1000.0
     case('Density2')
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,RhoO_ ,1:nLine)*1000.0
     case('Density3')
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,RhoHe_,1:nLine)*1000.0
     case('Velocity1')
        ! cm/s = 0.01*m/s
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,uH_ ,1:nLine)*0.01
     case('Velocity2')
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,uO_ ,1:nLine)*0.01
     case('Velocity3')
        SendBuffer_VI(iVar,:)=State_CVI(nAlt,uHe_,1:nLine)*0.01
     case default
        call CON_stop(NameSub//': unknown variable name='//Name_V(iVar))
     end select
  enddo

  ! The recieve buffer is Buffer_VI allocated in CON_couple_pw_gm
  iSendCount=nVar*nLine


  ! Gather all data to the root processor 
  call MPI_GATHERV(SendBuffer_VI, iSendCount, MPI_REAL, &
       Buffer_VI, iRecieveCount_P, iDisplacement_P, MPI_REAL, &
       0, iComm, iError)

  if(DoTestMe)then
     write(*,*)NameSub,' nVar, nLineTotal=',nVar, nLineTotal
     write(*,*)NameSub,' SendBuffer_VI=',SendBuffer_VI
  end if

  deallocate(iRecieveCount_P, iDisplacement_P, SendBuffer_VI)

end subroutine PW_get_for_gm

!==============================================================================

subroutine PW_put_from_gm(nTotalLine,Buffer_I)
  use ModPWOM, ONLY: nLine,iLineGlobal
  use ModGmPressure
  use CON_coupler, ONLY: Couple_CC, GM_, PW_
  implicit none
  
  integer,intent(in) :: nTotalLine
  real, intent(in)   :: Buffer_I(nTotalLine)
  integer            :: iLine
  !----------------------------------------------------------------------------
  
  !Make sure pressure array is allocated, if it isn't then allocate
  if(.not.allocated(p_I)) then
     allocate(p_I(nLine))
  end if
  
  !Fill pressure array with GM pressure
  do iLine=1,nLine
     p_I(iLine) = Buffer_I(iLineGlobal(iLine))
  enddo
  
  !Convert from SI to CGS
  p_I = p_I * 10.0  ! N/m^2 --> dynes/cm^2
  
  UseGmToPw = Couple_CC(GM_, PW_) % DoThis  
  
end subroutine PW_put_from_gm
