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
  use ModPWOM, ONLY: UseIE
  use CON_coupler, ONLY: Couple_CC, IE_, PW_
  
  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  
  character(len=*), parameter :: NameSub='PW_init_session'

  logical :: DoInitialize = .true.
  !----------------------------------------------------------------------------
  UseIE = Couple_CC(IE_, PW_) % DoThis
  
  if(DoInitialize) call PW_initialize
  DoInitialize = .false.

end subroutine PW_init_session

!==============================================================================

subroutine PW_finalize(TimeSimulation)

  use ModPWOM, ONLY: iLine, nLine, iUnitGraphics, iUnitOutput

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_finalize'
  !-------------------------------------------------------------------------
  do iLine=1,nLine
     close(UNIT=iUnitGraphics(iLine))
  enddo
  close(UNIT=iUnitOutput)

end subroutine PW_finalize

!==============================================================================

subroutine PW_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='PW_save_restart'

  call CON_stop(NameSub//': PW_ERROR: not yet implemented!')

end subroutine PW_save_restart

!==============================================================================

subroutine PW_run(TimeSimulation,TimeSimulationLimit)

  use ModPWOM, ONLY: iLine, nLine, Time, DtHorizontal,DToutput

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'
  !---------------------------------------------------------------------------
  DtHorizontal = min(DtHorizontal, TimeSimulationLimit - Time)
  
  do iLine=1,nLine
     call move_line
     call PW_advance_line
  end do
  if (floor(Time/DToutput) .ne. floor((Time-DtHorizontal)/DToutput) ) &
       call PW_print_electrodynamics
  TimeSimulation = Time

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
                      dOxyg_CI,dHyd_CI,dHel_CI,uOxyg_CI,uHyd_CI,uHel_CI,&
                      nLine,nAlt,&
                      nLine_P, nLineBefore_P

  use ModMpi

  implicit none
  character (len=*),parameter :: NameSub='PW_get_for_gm'

  integer, intent(in)           :: nVar, nLineTotal
  real, intent(out)             :: Buffer_VI(nVar, nLineTotal)
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar,i,iError
  integer :: iSendCount

  integer, save, allocatable :: iDisplacement_P(:), iRecieveCount_P(:)
  real, save, allocatable :: SendBuffer_VI(:,:)

  logical :: DoTest, DoTestMe
  !--------------------------------------------------------------------------
  call CON_set_do_test(NameSub, DoTest, DoTestMe)

  ! The sizes and displacements of MPI messages
  if(.not.allocated(iRecieveCount_P)) then
     allocate(iRecieveCount_P(nProc), iDisplacement_P(nProc))
     allocate(SendBuffer_VI(nVar, nLine))
     iRecieveCount_P = nVar*nLine_P
     iDisplacement_P = nVar*nLineBefore_P
  end if

  ! Fill buffer using SI units for sending on each proc
  do iVar=1,nVar
     select case (Name_V(iVar))
     case('CoLat')
        SendBuffer_VI(iVar,:)=ThetaLine_I(1:nLine)
     case('Longitude')
        SendBuffer_VI(iVar,:)=PhiLine_I(1:nLine)
     case('Density1')
        ! g/cm^3 = 1000 * kg/m^3
        SendBuffer_VI(iVar,:)=dOxyg_CI(nAlt,1:nLine)*1000.0
     case('Density2')
        SendBuffer_VI(iVar,:)=dHyd_CI(nAlt,1:nLine)*1000.0
     case('Density3')
        SendBuffer_VI(iVar,:)=dHel_CI(nAlt,1:nLine)*1000.0
     case('Velocity1')
        ! cm/s = 0.01*m/s
        SendBuffer_VI(iVar,:)=uOxyg_CI(nAlt,1:nLine)*0.01
     case('Velocity2')
        SendBuffer_VI(iVar,:)=uHyd_CI(nAlt,1:nLine)*0.01
     case('Velocity3')
        SendBuffer_VI(iVar,:)=uHel_CI(nAlt,1:nLine)*0.01
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

end subroutine PW_get_for_gm
