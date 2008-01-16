
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
