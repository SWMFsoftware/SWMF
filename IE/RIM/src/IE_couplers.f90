!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

!------------------------------------------------------------------------
! Routines for setting the number of grid cells
!------------------------------------------------------------------------

subroutine IE_setnMlts(iComponent, nMLTsIn, iError)

  use CON_comp_param, only: MaxComp
  use ModIEInterface

  implicit none

  integer, intent(in)  :: iComponent, nMLTsIn
  integer, intent(out) :: iError

  iError = 0

  if (iComponent > MaxComp .or. iComponent <= 0) then
     write(*,*) "Error in IE_setmlts. Component not correct!"
     write(*,*) "Valid Range for component : ",1,MaxComp
     write(*,*) "Value sent in : ",iComponent
     iError = 1
  else
     iecouplerinfo(iComponent)%nLons = nMLTsIn
  endif

end subroutine IE_setnMlts

!--------------------------------------------------------------

subroutine IE_setnLats(iComponent, nLatsIn, iError)

  use CON_comp_param, only: MaxComp
  use ModIEInterface

  implicit none

  integer, intent(in)  :: iComponent, nLatsIn
  integer, intent(out) :: iError

  iError = 0

  if (iComponent > MaxComp .or. iComponent <= 0) then
     write(*,*) "Error in IE_setmlts. Component not correct!"
     write(*,*) "Valid Range for component : ",1,MaxComp
     write(*,*) "Value sent in : ",iComponent
     iError = 1
  else
     iecouplerinfo(iComponent)%nLats = nLatsIn
  endif

end subroutine IE_setnLats

!--------------------------------------------------------------

subroutine IE_setgrid(iComponent, MLTsIn, LatsIn, iError)

  use CON_comp_param
  use ModIEInterface

  integer, intent(in) :: iComponent
  real, dimension(iecouplerinfo(iComponent)%nLons,&
                  iecouplerinfo(iComponent)%nLats), intent(in) :: MLTsIn,LatsIn
  integer, intent(out) :: iError

  integer :: nLatsIE, nLonsIE
  real, allocatable :: LatitudeIE(:,:), LongitudeIE(:,:)

  logical :: IsIEGridRefined = .true.

  iError = 0

  if (.not.IsIEGridRefined) then
     ! Check to see if we need to redo grid - should be NO!
  endif

  if (IsIEGridRefined) then
     ! get IE grid (i.e. message pass from iProc = 0 on IE to here)
     ! make new interpolation indices (stored in allocable arrays above)
     IsIEGridRefined = .false.
  else

  end if

end subroutine IE_setgrid
