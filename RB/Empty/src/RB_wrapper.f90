!^CFG COPYRIGHT UM

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!               Space Weather Modeling Framework (SWMF)                !
!    Center for Space Environment Modeling, The University of Michigan !
!-----------------------------------------------------------------------
!BOI
! !TITLE: Wrapper for the "empty" Radiation Belts Component 
! !AUTHORS: Ovsei Volberg
! !AFFILIATION: CSEM, The University of Michigan
! !DATE: April 13, 2004 - the inital version was written
! !INTRODUCTION: This wrapper provides the "empty" interface
!                for the  Radiation Belts Component             
!EOI
!-------------------------------------------------------------------------

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  RB_set_param
!
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_set_param(CompInfo, TypeAction)
!
! !USES:
!
  use CON_comp_info
 
  implicit none
!
! !INPUT PARAMETERS: 
!
  character (len=*), intent(in)     :: TypeAction ! which action to perform
!
! !INPUT/OUTPUT PARAMETERS: 
!
  type(CompInfoType), intent(inout) :: CompInfo   ! component information

! !OUTPUT PARAMETERS:
!
! !FILES USED:  CON_comp_info.f90
!
! !REVISION HISTORY: 
!
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:

  character (len=*), parameter :: NameSub='RB_set_param'

  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,                         &
          Use=.false.,                           &
          NameVersion='Empty', &
          Version=0.0)
  case default
     call CON_stop(NameSub//' RB_ERROR: empty version cannot be used!')
  end select

end subroutine RB_set_param
!============================================================================

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  RB_set_grid
!
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_set_grid
!
! !USES:
!
  implicit none
!
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
! !OUTPUT PARAMETERS:
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
! 
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:
  character (len=*), parameter :: NameSub='RB_set_grid'

 
    call CON_stop(NameSub//' RB_ERROR: empty version cannot be used!')
 
end subroutine RB_set_grid
!==============================================================================

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  RB_init_session
!
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_init_session(iSession, TimeSimulation)
!
! !USES:
!
  implicit none
!
! !INPUT PARAMETERS: 
!
 
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time
!
! !INPUT/OUTPUT PARAMETERS: 
!
! !OUTPUT PARAMETERS:
!
! !FILES USED:  
!
! !REVISION HISTORY: 
! 
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:

  character(len=*), parameter :: NameSub='RB_init_session'

  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')

end subroutine RB_init_session
!============================================================================

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  RB_run
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_run(TimeSimulation,TimeSimulationLimit)

! !USES:
  implicit none

! !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

! !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

! !OUTPUT PARAMETERS:
!
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
!  
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:

  character(len=*), parameter :: NameSub='RB_run'
  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
end subroutine RB_run
!===========================================================================

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  
!
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_finalize(TimeSimulation)

! !USES:
  implicit none

! !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time
! !OUTPUT PARAMETERS:
!
! !BUGS:  
!
! !SEE ALSO: 
!
! !SYSTEM ROUTINES: 
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
!  
! 
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:
  character(len=*), parameter :: NameSub='RB_finalize'

  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
end subroutine RB_finalize
!===========================================================================
!==============================================================================

!BOP -------------------------------------------------------------------------
!
! !ROUTINE:  RB_save_restart
!
! !DESCRIPTION: 
! 
! !INTERFACE:
!
subroutine RB_save_restart(TimeSimulation)
!
! !USES:
!
  implicit none
!
! !INPUT PARAMETERS: 
!
  real,     intent(in) :: TimeSimulation   ! seconds from start time
!
! !INPUT/OUTPUT PARAMETERS: 
!
! !OUTPUT PARAMETERS:
!
! !FILES USED:  
!
! !REVISION HISTORY: 
!
!  
! 
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:

  character(len=*), parameter :: NameSub='RB_save_restart'

  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
end subroutine RB_save_restart
!==============================================================================

!BOP -------------------------------------------------------------------------
!
!ROUTINE:  RB_put_from_gm
!
!DESCRIPTION: 
! 
!INTERFACE:
!
subroutine RB_put_from_gm(Buffer_IIV,iSizeIn,jSizeIn,nVarIn,NameVar)
  !
  !USES:
  !
  implicit none
  !
  !INPUT PARAMETERS: 
  !
  integer, intent(in) :: iSizeIn,jSizeIn,nVarIn
  real, dimension(iSizeIn,jSizeIn,nVarIn), intent(in) :: Buffer_IIV
  character (len=*),intent(in) :: NameVar
  !
  !INPUT/OUTPUT PARAMETERS: 
  !
  !OUTPUT PARAMETERS:
  !
  !FILES USED:  
  !
  !REVISION HISTORY: 
  !
  !
  !
!EOP -------------------------------------------------------------------------

!LOCAL VARIABLES:

  character(len=*), parameter :: NameSub='RB_put_from_gm'

  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
end subroutine RB_put_from_gm
!============================================================================
subroutine RB_put_from_ie(Buffer_IIV, iSize, jSize, nVarIn, &
                 Name_V, iBlock)
 
  implicit none

  character(len=*), parameter :: NameSub='RB_put_from_ie'

  !INPUT ARGUMENTS:
  integer, intent(in):: iSize, jSize, nVarIn, iBlock
  real, intent(in) :: Buffer_IIV(iSize, jSize, nVarIn)
  character(len=*), intent(in) :: Name_V(nVarIn)
  
  call CON_stop(NameSub//': RB_ERROR: empty version cannot be used!')
end subroutine RB_put_from_ie
!==============================================================================
subroutine RB_put_sat_from_gm(nSats, Buffer_I, Buffer_III)
  implicit none

  integer, intent(in)            :: nSats
  real, intent(in)               :: Buffer_III(4,2,nSats)
  character(len=100), intent(in) :: Buffer_I(nSats)
end subroutine RB_put_sat_from_gm
