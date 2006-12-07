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
  use ModPWOM, ONLY: Theta_G,Phi_G,SigmaH_G,SigmaP_G,Jr_G,Potential_G,&
       NamePhiNorth,iPhi,nPhi,iTheta,nTheta, allocate_ie_variables
  use ModIoUnit, ONLY: UnitTmp_
  
  implicit none
  

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
  logical  :: UseIE=.false.
  character(len=*), parameter :: NameSub='PW_init_session'

  logical :: DoInitialize = .true.
  !----------------------------------------------------------------------------
  if(DoInitialize) call PW_initialize
  DoInitialize = .false.
  
  if (.not. UseIE) then
     open(UnitTmp_, FILE=NamePhiNorth)
     call allocate_ie_variables(257, 65)
     do iPhi=1,nPhi
        do iTheta=1,nTheta
           read(unit=UnitTmp_,fmt='(6(1PE13.5))') &
               Theta_G(iPhi,iTheta),Phi_G(iPhi,iTheta),SigmaH_G(iPhi,iTheta),&
               SigmaP_G(iPhi,iTheta),Jr_G(iPhi,iTheta),Potential_G(iPhi,iTheta)
           
        enddo
     enddo
     close(UnitTmp_)
  endif
  

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

  use ModPWOM, ONLY: iLine, nLine, Time, DtMax, Dt

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='PW_run'
  logical :: IsFirstCall=.true.
  !---------------------------------------------------------------------------

  Dt = min(DtMax, TimeSimulationLimit - Time)
  !call PW_get_electrodynamic
  if (IsFirstCall) then
     call PW_get_electrodynamic
     call Initial_Line_Location
     IsFirstCall=.false.
  endif
  
  do iLine=1,nLine
     call MoveFluxTube
     call PW_advance_line
  end do
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

  end if


  IsPotFound = .false.
  IsJrFound  = .false.
  do iVar = 1, nVarIn
     select case(Name_V(iVar))
     case('Pot')
        IsPotFound = .true.
        Potential_G(1:jSize, 1:iSize) = transpose(Buffer_IIV(:, :, iVar))
     case('Jr')
        IsJrFound = .true.
        Jr_G(1:jSize, 1:iSize) = transpose(Buffer_IIV(:, :, iVar))
     end select
  end do
  if(.not.IsPotFound .or. .not.IsJrFound)then
     write(*,*)NameSub,': Name_V=',Name_V
     call CON_stop(NameSub//' could not find Pot or Jr')
  end if

end subroutine PW_put_from_ie
