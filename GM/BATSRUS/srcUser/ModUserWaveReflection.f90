!^CFG FILE USERFILES
module ModUser

  use ModUserEmpty, ONLY:               &
       user_read_inputs,                &
!!!       user_set_ics,                    &
       user_initial_perturbation,       &
       user_set_boundary_cells,        &
       user_face_bcs,                   &
       user_set_outerbcs,               &
       user_specify_initial_refinement, &
!!!       user_amr_criteria,               &
       user_write_progress,             &
       user_get_log_var,                &
       user_calc_sources,               &
       user_heat_source,                &
       user_get_b0,                     &
       user_update_states

  include 'user_module.h'

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'USER ROUTINES FOR IMPLICIT PAPER, G. Toth'

  real :: xLeft, xRight

contains

  subroutine user_set_ics
    use ModMain,     ONLY: nI, nJ, nK, globalBLK, nBlock
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, P_, Bx_, By_
    use ModPhysics,  ONLY: ShockSlope, g

    real :: Potential_G(-1:nI+2,-1:nJ+2)
    real, parameter :: pPerturb = 1.1
    real :: SinSlope, CosSlope, pTotal
    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    !  if(iBlock==nBlock)write(*,*)' !!! Setting Rho=1000-2x-y !!!'

    !  State_VGB(Rho_,:,:,:,iBlock) = 1000.0 &
    !  -2*x_BLK(:,:,:,iBlock) -y_BLK(:,:,:,iBlock) !!! + 3*z_BLK(:,:,:,iBlock)

    !  RETURN

    xLeft  = 25.0
    xRight = 30.0

    if(ShockSlope == 0.0)then
       ! Perturb pressure and return
       where(     x_BLK(:,:,:,iBlock) >= xLeft &
            .and. x_BLK(:,:,:,iBlock) <= xRight) &
            State_VGB(P_,:,:,:,iBlock)=pPerturb*State_VGB(P_,:,:,:,iBlock)
       RETURN
    endif

    ! Store total pressure
    pTotal = State_VGB(P_,1,1,1,iBlock) + &
         0.5*(State_VGB(Bx_,1,1,1,iBlock)**2 &
         +    State_VGB(By_,1,1,1,iBlock)**2 )

    ! Calculate rotated magnetic field
    ! Original vector potential: A_z = 0.1*y - 100*min(0,x)
    ! Rotated  vector potential: A_z = 0.1*(Cos*y-Sin*x)-100*min(0,Cos*x+Sin*y)

    CosSlope = cos(atan(ShockSlope))
    SinSlope = sin(atan(ShockSlope))

    Potential_G = &
         0.1*(CosSlope*Y_BLK(:,:,1,iBlock)-SinSlope*x_BLK(:,:,1,iBlock)) &
         -100*min(0., &
         CosSlope*x_BLK(:,:,1,iBlock) + SinSlope*Y_BLK(:,:,1,iBlock))

    ! B = curl A so Bx = dA_z/dy and By = -dAz/dx
    do j=1,nJ; do i=1,nI
       State_VGB(Bx_,i,j,:,iBlock) = &
            +(Potential_G(i,j+1)-Potential_G(i,j-1)) / (2*Dy_BLK(iBlock))
       State_VGB(By_,i,j,:,iBlock) = &
            -(Potential_G(i+1,j)-Potential_G(i-1,j)) / (2*Dx_BLK(iBlock))
    end do; end do

    ! Recalculate pressure
    State_VGB(P_,1:nI,1:nJ,1:nK,iBlock) = pTotal - &
         0.5*(State_VGB(Bx_,1:nI,1:nJ,1:nK,iBlock)**2 &
         +    State_VGB(By_,1:nI,1:nJ,1:nK,iBlock)**2 )

    ! Recalculate momentum and density
    do k=1,nK; do j=1,nJ; do i=1,nI
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
            State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) * &
            g * State_VGB(P_,i,j,k,iBlock) / State_VGB(Rho_,i,j,k,iBlock)

       State_VGB(Rho_,i,j,k,iBlock) = g * State_VGB(P_,i,j,k,iBlock)
    end do; end do; end do

    ! Rotate pressure perturbation parameters
    CosSlope = cos(atan(ShockSlope))
    xLeft  = xLeft /CosSlope
    xRight = xRight/CosSlope

    ! Perturb pressure
    where(     x_BLK(:,:,:,iBlock) >= xLeft -ShockSlope*Y_BLK(:,:,:,iBlock) &
         .and. x_BLK(:,:,:,iBlock) <= xRight-ShockSlope*Y_BLK(:,:,:,iBlock)) &
         State_VGB(P_,:,:,:,iBlock) = pPerturb*State_VGB(P_,:,:,:,iBlock)

  end subroutine user_set_ICs

  !========================================================================

  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    use ModSize, ONLY     : nI, nJ, nK
    use ModGeometry, ONLY : x_BLK, y_BLK, dx_BLK
    use ModPhysics,  ONLY : ShockSlope
    use ModMain, ONLY     : time_simulation

    ! Variables required by this user subroutine
    character (len=*),intent(in) :: TypeCriteria
    integer, intent(in)          :: iBlock
    real, intent(out)            :: UserCriteria
    logical ,intent(inout)       :: IsFound

    real :: xCenter, yCenter, xShifted, x_I(5), Distance
    !------------------------------------------------------------------
    xCenter = 0.5*(x_BLK(nI,nJ,nK,iBlock)+x_BLK(1,1,1,iBlock))
    yCenter = 0.5*(y_BLK(nI,nJ,nK,iBlock)+y_BLK(1,1,1,iBlock))
    xShifted = xCenter + ShockSlope*yCenter

    x_I(1) = xLeft  - 0.9*time_simulation
    x_I(2) = xLeft  + 1.1*time_simulation
    x_I(3) = xRight - 0.9*time_simulation
    x_I(4) = xRight + 1.1*time_simulation
    x_I(5) = -0.1*time_simulation

    Distance = minval(abs(x_I - xShifted)) / dx_BLK(iBlock)

    if(Distance <= nI/2+2)then
       ! Block cuts one of the discontinuities
       ! The larger the cell size the more important to refine it
       UserCriteria = 10000.0*dx_BLK(iBlock)
    else
       ! The block does not cut through. The closer it is the more
       ! useful a refinement is.
       UserCriteria = 1./Distance
    endif

    IsFound = .true.

  end subroutine user_amr_criteria

end module ModUser
