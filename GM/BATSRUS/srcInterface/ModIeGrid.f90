!^CFG COPYRIGHT UM
!^CMP FILE IE
module ModIeGrid

  ! This module contains easily accessible information about the IE grid

  implicit none
  save

  integer           :: nThetaIono = -1, nPhiIono = -1
  real, allocatable :: ThetaIono_I(:), PhiIono_I(:)
  real              :: rIonosphere

  real, private     :: dThetaIono, dPhiIono

contains

  subroutine init_mod_ie_grid(iSize, jSize)
    ! The ionosphere works on two hemispheres with a node based grid
    ! iSize is the number of latitude nodes from the pole to the equator.
    ! jSize is the number of longitude nodes (including a periodic overlap)

    use ModNumConst, ONLY: cPi, cTwoPi
    use ModPhysics,  ONLY: Si2No_V, UnitX_
    use CON_coupler, ONLY: Grid_C, IE_

    integer, intent(in) :: iSize, jSize
    real :: rPlanet, IonoHeight
    character(len=*), parameter :: NameSub='init_mod_ie_grid'
    !-------------------------------------------------------------------------

    if(nThetaIono > 0) RETURN

    nThetaIono = Grid_C(IE_) % nCoord_D(1)
    nPhiIono   = Grid_C(IE_) % nCoord_D(2)

    if(nThetaIono /= iSize .or. nPhiIono /= jSize)then
       write(*,*)NameSub,': Grid_C(IE_)%nCoord_D(1:2)=',&
            Grid_C(IE_) % nCoord_D(1:2)
       write(*,*)NameSub,': iSize,2*iSize-1,jSize=',iSize,2*iSize-1,jSize
       call stop_mpi(NameSub//' ERROR: Inconsistent IE grid sizes')
    endif

    allocate(ThetaIono_I(nThetaIono), PhiIono_I(nPhiIono))
    ThetaIono_I = Grid_C(IE_) % Coord1_I
    PhiIono_I   = Grid_C(IE_) % Coord2_I
    rIonosphere = Grid_C(IE_) % Coord3_I(1) * Si2No_V(UnitX_)

    dThetaIono = cPi    / (nThetaIono - 1)
    dPhiIono   = cTwoPi / (nPhiIono - 1)

  end subroutine init_mod_ie_grid
  !===========================================================================
  subroutine get_ie_grid_index(Theta, Phi, ThetaNorm, PhiNorm)
    real, intent(in) :: Theta, Phi
    real, intent(out):: ThetaNorm, PhiNorm

    ThetaNorm = Theta / dThetaIono
    PhiNorm   = Phi   / dPhiIono
  end subroutine get_ie_grid_index

end module ModIeGrid
