module ModIonoVelocity

  implicit none
  SAVE
  private ! except


  public :: read_iono_velocity_param
  public :: apply_iono_velocity

  ! Local variables 
  logical :: UseIonoVelocity = .false.
  real    :: rCoupleUiono    = 3.5
  real    :: TauCoupleUiono  = 20.0

contains

  subroutine read_iono_velocity_param

    use ModReadParam, ONLY: read_var

    call read_var('UseIonoVelocity', UseIonoVelocity)
    call read_var('rCoupleUiono'   , rCoupleUiono)
    call read_var('TauCoupleUiono',  TauCoupleUiono)

  end subroutine read_iono_velocity_param

  !===========================================================================

  subroutine apply_iono_velocity

    use ModMain,    ONLY: nI, nJ, nK, nBlock, time_accurate, time_simulation, &
         Dt, UnusedBlk, UseB0, UseRotatingBc
    use ModAdvance, ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, Bx_, Bz_
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, r_BLK, Rmin_BLK
    use ModB0,      ONLY: B0_DGB
    use ModPhysics, ONLY: Si2No_V, UnitT_, rBody
    use ModEnergy,  ONLY: calc_energy_cell

    real :: Factor, RhoUdotB
    real, dimension(3) :: Xyz_D, Uiono_D, Urot_D, RhoUiono_D, RhoU_D, b_D

    integer :: i,j,k,iBlock
    character (len=*), parameter :: NameSub='apply_iono_velocity'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    if(.not.UseIonoVelocity) RETURN

    call set_oktest(NameSub, DoTest, DoTestMe)

    if(time_accurate)then
       ! Ramp up is based on physical time: u' = u + dt/tau * (uIE - u)
       ! A typical value might be 10 sec, to get close to the IE velocity
       ! in 20 seconds

       Factor = min(1.0, Dt/(TauCoupleUiono*Si2No_V(UnitT_)))

    else
       ! Ramp up is based on number of iterations: u' = u + (uIE-u)/(1+nTau)
       ! A typical value might be 10, to get close to the E x B velocity
       ! in 20 iterations

       Factor = 1.0/(1.0 + TauCoupleUiono)

    end if

    do iBlock = 1, nBlock
       if(unusedBLK(iBlock)) CYCLE
       if(rMin_BLK(iBlock) > rCoupleUiono) CYCLE

       do k = 1, nK; do j = 1, nJ; do i = 1, nI

          if( r_BLK(i,j,k,iBlock) > rCoupleUiono) CYCLE
          if( r_BLK(i,j,k,iBlock) < rBody) CYCLE

          Xyz_D = (/ x_BLK(i,j,k,iBlock), &
               y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock) /)

          b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
          
          if(UseB0)b_D = b_D + b0_DGB(:,i,j,k,iBlock)

          ! Calculate E x B velocity
          call calc_inner_bc_velocity(Time_Simulation, Xyz_D, b_D, uIono_D)

          ! Add rotational velocity if required
          if (UseRotatingBc) then
             call calc_corotation_velocities(Xyz_D, uRot_D)
             uIono_D = uIono_D + uRot_D
          end if

          ! Convert to momentum
          RhoUIono_D = State_VGB(Rho_,i,j,k,iBlock)*uIono_D

          ! The original momentum
          RhoU_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)

          ! Store field aligned momentum component
          RhoUdotB = sum(RhoU_D*b_D)

          ! Push momenta towards the Rho*(E x B + Omega x r) value
          RhoU_D = RhoU_D + Factor*(RhoUIono_D - RhoU_D)

          ! Restore field aligned momentum component
          RhoU_D = RhoU_D + b_D*(RhoUdotB - sum(RhoU_D*b_D))/sum(b_D**2)
          
          ! Store result
          State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = RhoU_D

       end do; end do; end do

       ! Recalculate the energy
       call calc_energy_cell(iBlock)

    end do

  end subroutine apply_iono_velocity

end module ModIonoVelocity
