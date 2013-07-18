!This code is a copyright protected software (c) 2002- University of Michigan

!==============================================================================
module ModHeatFluxCollisionless

  implicit none
  save

  private ! except

  ! Public methods
  public :: read_heatflux_collisionless_param
  public :: get_gamma_collisionless
  public :: update_heatflux_collisionless

  ! Parameters for heat flux region
  logical, public :: UseHeatFluxRegion = .false.
  real, public :: rCollisional, rCollisionless

  ! Parameters for collisionless heat conduction
  logical, public :: UseHeatFluxCollisionless = .false.
  real :: CollisionlessAlpha = 1.05

contains

  !============================================================================

  subroutine read_heatflux_collisionless_param(NameCommand)

    use ModReadParam, ONLY: read_var

    character(len=*), intent(in) :: NameCommand

    character(len=*), parameter :: &
         NameSub = 'ModHeatFluxCollision::read_heatflux_collisionless_param'
    !--------------------------------------------------------------------------

    select case(NameCommand)

    case("#HEATFLUXREGION")
       call read_var('UseHeatFluxRegion', UseHeatFluxRegion)
       if(UseHeatFluxRegion)then
          call read_var('rCollisional', rCollisional)
          call read_var('rCollisionless', rCollisionless)
       end if

    case("#HEATFLUXCOLLISIONLESS")
       call read_var('UseHeatFluxCollisionless', UseHeatFluxCollisionless)
       if(UseHeatFluxCollisionless)then
          call read_var('CollisionlessAlpha', CollisionlessAlpha)
       endif

    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

  end subroutine read_heatflux_collisionless_param

  !============================================================================

  subroutine update_heatflux_collisionless(iBlock)

    use BATL_lib,      ONLY: Xyz_DGB
    use ModVarIndexes, ONLY: Pe_, P_, Ehot_
    use ModSize
    use ModAdvance,    ONLY: State_VGB, UseElectronPressure
    use ModPhysics,    ONLY: inv_gm1
    use ModEnergy,     ONLY: calc_energy_cell

    integer, intent(in) :: iBlock

    integer:: i, j, k, iP
    real:: Gamma
    !--------------------------------------------------------------------------

    ! We use a varying gamma for the electrons to parameterize the
    ! collisionless heat flux of Hollweg (1976).

    iP = p_
    if(UseElectronPressure) iP = Pe_

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       call get_gamma_collisionless(Xyz_DGB(:,i,j,k,iBlock), Gamma)

       State_VGB(iP,i,j,k,iBlock) = (Gamma - 1) &
            *(inv_gm1*State_VGB(iP,i,j,k,iBlock) &
            + State_VGB(Ehot_,i,j,k,iBlock))
       State_VGB(Ehot_,i,j,k,iBlock) = State_VGB(iP,i,j,k,iBlock) &
            *(1.0/(Gamma - 1) - inv_gm1)
    end do; end do; end do

    call calc_energy_cell(iBlock)

  end subroutine update_heatflux_collisionless

  !============================================================================

  subroutine get_gamma_collisionless(x_D, Gamma)

    use BATL_lib,   ONLY: MaxDim
    use ModPhysics, ONLY: g, inv_gm1
    use ModAdvance, ONLY: UseElectronPressure
    use ModMain,    ONLY: UseHeatConduction

    real, intent(in) :: x_D(MaxDim)
    real, intent(out) :: Gamma

    real :: r
    real :: GammaCollisionless
    !--------------------------------------------------------------------------

    r = sqrt(sum(x_D**2))

    if(UseElectronPressure)then
       GammaCollisionless = (inv_gm1*g + 1.5*CollisionlessAlpha) &
            /(inv_gm1 + 1.5*CollisionlessAlpha)
    else
       GammaCollisionless = (inv_gm1*g + 0.75*CollisionlessAlpha) &
            /(inv_gm1 + 0.75*CollisionlessAlpha)
    end if

    if(.not.UseHeatConduction)then
       Gamma = GammaCollisionless
    elseif(rCollisionless < 0.0)then
       Gamma = GammaCollisionless &
            + (g-GammaCollisionless)/((r/rCollisional)**2 + 1)
    elseif(r <= rCollisional)then
       Gamma = g
    else
       Gamma = GammaCollisionless + (g-GammaCollisionless)* &
            exp(-((r-rCollisional)/(rCollisionless-rCollisional))**2)
    end if

  end subroutine get_gamma_collisionless

end module ModHeatFluxCollisionless
