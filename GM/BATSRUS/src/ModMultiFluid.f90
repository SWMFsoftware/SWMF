module ModMultiFluid

  use ModVarIndexes

  implicit none
  save

  integer :: iRho, iRhoUx, iRhoUy, iRhoUz, iP, iEnergy, iUx, iUy, iUz

contains

  subroutine select_fluid(iFluid)
    integer, intent(in) :: iFluid
    integer :: i

    if(iFluid == 1)then
       iRho   = Rho_
       iRhoUx = RhoUx_
       iRhoUy = RhoUy_
       iRhoUz = RhoUz_
       iP     = P_
    else
       i = iVarFluid_I(iFluid)
       iRho   = i + 1
       iRhoUx = i + 2
       iRhoUy = i + 3
       iRhoUz = i + 4
       iP     = i + 5
    end if
    iEnergy= nVar + iFluid

    iUx    = iRhoUx
    iUy    = iRhoUy
    iUz    = iRhoUz

  end subroutine select_fluid

end module ModMultiFluid
