!^CFG COPYRIGHT UM
! Wrapper for the empty PWOM (PW) component
!==========================================================================
subroutine PW_set_param(CompInfo, TypeAction)

  use CON_comp_info
  use ModIoUnit, only: STDOUT_
  use ModPWOM, only: iUnitOut, iProc, nProc, iComm, StringPrefix

  implicit none

  character (len=*), parameter :: NameSub='PW_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
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
     ! Do nothing
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

  use ModPWOM, ONLY: iLine, nLine, Time, DtMax, Dt,DToutput

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'
  !---------------------------------------------------------------------------
  Dt = min(DtMax, TimeSimulationLimit - Time)
  
  do iLine=1,nLine
     call MoveFluxTube
     call PW_advance_line
  end do
  if (floor(Time/DToutput) .ne. floor((Time-Dt)/DToutput) ) &
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

  call PW_get_electrodynamic

end subroutine PW_put_from_ie
!==============================================================================

subroutine PW_get_for_gm(Buffer_IIV, nFieldLines, nVar, Name_V, &
     tSimulation)

  use ModPWOM, only : icomm,errcode,nProc,nTotalLine,nVar,&
                      FieldLineTheta,FieldLinePhi, &
                      dOxyg,dHyd,dHel,uOxyg,uHyd,uHel

  implicit none
  character (len=*),parameter :: NameSub='PW_get_for_gm'

  integer, intent(in)           :: nTotalLine, nVar
  real, intent(out)             :: Buffer_IIV(iSize,jSize,nVar)
  character (len=*),intent(in)  :: Name_V(nVar)
  real,             intent(in)  :: tSimulation

  integer :: iVar
  real    :: tSimulationTmp,
  real    :: SendBuffer(nLine,nVar)
  integer :: iDisplacement_V(nProc),iRecieveCount_V(nProc)
  !--------------------------------------------------------------------------
  if(iSize /= IONO_nTheta .or. jSize /= IONO_nPsi)then
     write(*,*)NameSub//' incorrect buffer size=',iSize,jSize,&
          ' IONO_nTheta,IONO_nPsi=',IONO_nTheta, IONO_nPsi
     call CON_stop(NameSub//' SWMF_ERROR')
  end if

  ! Make sure that the most recent result is provided
  tSimulationTmp = tSimulation
  call PW_run(tSimulationTmp,tSimulation)

  ! Prepare buffer for sending on each proc
  do iVar=1,nVar
     select case (Name_V(iVar))
     case('CoLat    ')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=FieldLineTheta(iLine)
        enddo
        
     case('Longitude')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=FieldLinePhi(iLine)
        enddo
        
     case('Density1 ')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=dOxyg(nDim,iLine)
        enddo
        
     case('Density2 ')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=dHyd(nDim,iLine)
        enddo
        
     case('Density3 ')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=dHel(nDim,iLine)
        enddo
        
     case('Velocity1')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=uOxyg(nDim,iLine)
        enddo
        
     case('Velocity2')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=uHyd(nDim,iLine)
        enddo
        
     case('Velocity3')
        do iLine=1,nLine
           SendBuffer(iLine,iVar)=uHel(nDim,iLine)
        enddo
        
     end select
     
  enddo

  ! The recieve buffer is Buffer_IIV allocated in CON_couple_pw_gm
  
  ! create the displacement array and recieve count array for MPI_GATHERV
  do i=0,nProc-1
     if (i .lt. mod(nTotalLine,nProc)) then
        iDisplacement_V(i)=&
             i*ceiling(real(nTotalLine)/real(nProc))+1
        
        iRecieveCount_V(i)=&
             ceiling(real(nTotalLine)/real(nProc))
     else
        iDisplacement_V(i)=&
             (mod(nTotalLine,nProc))*ceiling(real(nTotalLine)/real(nProc)) &
             + ((i)-mod(nTotalLine,nProc))                        &
             *floor(real(nTotalLine)/real(nProc))+1
        
        iDisplacement_V(i)=&
             floor(real(nTotalLine)/real(nProc))
     endif
  enddo
  iSendCount=nLine
  
  ! Gather all data to be passed on the root processor (0)
  
    do iVar=1,nVar
       
       call MPI_GATHERV(SendBuffer(1,iVar), iSendCount, mpi_real, &
            Buffer_IIV(1,iVar), iRecieveCount,iDisplacement_V, mpi_real,&
            0, icomm, errcode)
    
    enddo
  
  

end subroutine PW_get_for_gm
