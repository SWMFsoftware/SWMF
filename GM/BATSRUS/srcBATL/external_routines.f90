! Collection of external routines called by BATL
subroutine user_specify_region(iArea, iBlock, nValue, NameLocation, &
     IsInside, IsInside_I, Value_I)

  implicit none

  integer,   intent(in):: iArea        ! area index in BATL_region
  integer,   intent(in):: iBlock       ! block index
  integer,   intent(in):: nValue       ! number of output values
  character, intent(in):: NameLocation ! c, g, x, y, z, or n

  logical, optional, intent(out) :: IsInside
  logical, optional, intent(out) :: IsInside_I(nValue)
  real,    optional, intent(out) :: Value_I(nValue)

  character(len=*), parameter:: NameSub = 'user_specify_region'
  !----------------------------------------------------------------------------
  call CON_stop(NameSub//' is not implemented')

end subroutine user_specify_region
!==============================================================================
