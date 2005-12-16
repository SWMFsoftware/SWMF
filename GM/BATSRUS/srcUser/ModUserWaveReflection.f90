module ModUser

  use ModUserEmpty, ONLY:               &
!!!       user_read_inputs,                &
!!!       user_init_session,               &
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
       'WAVE REFLECTION, G. Toth'

  real :: xLeft0 = 25.0, xRight0 = 30.0, cSoundX0 = 1.0
  real :: DistanceMin = 1.0
  real :: xLeft, xRight, cSoundX, Ux

contains

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam

    character (len=100) :: NameCommand
    !-------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#DISTANCE")
          call read_var('DistanceMin',DistanceMin)
       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do
  end subroutine user_read_inputs

  !=====================================================================

  subroutine user_init_session

    use ModProcMH,  ONLY: iProc
    use ModAdvance, ONLY: Ux_
    use ModPhysics, ONLY: ShockSlope, Shock_Lstate

    real :: CosSlope

    ! Rotate pressure perturbation parameters
    CosSlope = cos(atan(ShockSlope))
    xLeft   = xLeft0 /CosSlope
    xRight  = xRight0/CosSlope

    ! Also fix the sound speed projected to the X axis
    cSoundX = cSoundX0/CosSlope

    Ux = Shock_Lstate(Ux_)/CosSlope

    if(iProc==0)then
       write(*,*)'user_init_session: ShockSlope  =',ShockSlope
       write(*,*)'user_init_session: xLeft,xRight=',xLeft,xRight
       write(*,*)'user_init_session: cSoundX,  Ux=',cSoundX,Ux
    endif

  end subroutine user_init_session

  !=====================================================================

  subroutine user_set_ics
    use ModMain,     ONLY: nI, nJ, nK, globalBLK, nBlock
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK
    use ModAdvance,  ONLY: State_VGB, Rho_, RhoUx_, RhoUz_, P_, Bx_, By_
    use ModPhysics,  ONLY: ShockSlope, Shock_Lstate, g

    real :: Potential_G(-1:nI+2,-1:nJ+2)
    real, parameter :: pPerturb = 1.1
    real :: SinSlope, CosSlope, pTotal
    integer :: i, j, k, iBlock
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    if(ShockSlope == 0.0)then
       ! Perturb pressure and return
       where(     x_BLK(:,:,:,iBlock) >= xLeft0 &
            .and. x_BLK(:,:,:,iBlock) <= xRight0) &
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
         Shock_Lstate(Bx_)* &
         (CosSlope*Y_BLK(:,:,1,iBlock)-SinSlope*x_BLK(:,:,1,iBlock)) &
         -Shock_Lstate(By_)*min(0., &
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

    ! Perturb pressure
    where(     x_BLK(:,:,:,iBlock) >= xLeft -ShockSlope*Y_BLK(:,:,:,iBlock) &
         .and. x_BLK(:,:,:,iBlock) <= xRight-ShockSlope*Y_BLK(:,:,:,iBlock)) &
         State_VGB(P_,:,:,:,iBlock) = pPerturb*State_VGB(P_,:,:,:,iBlock)

  end subroutine user_set_ics

  !========================================================================

  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    use ModSize,     ONLY : nI, nJ, nK
    use ModGeometry, ONLY : x_BLK, y_BLK, dx_BLK
    use ModPhysics,  ONLY : ShockSlope
    use ModMain,     ONLY : time_simulation
    use ModAMR,      ONLY : RefineCritMin_I, CoarsenCritMax

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

    !write(*,*)'xLeft,xRight,Ux,cSoundX=',xLeft,xRight,Ux,cSoundX
    !call stop_mpi('Debug')

    ! Location of sound wave edges and the tangential discontinuity
    x_I(1) = xLeft  + (Ux-cSoundX)*time_simulation
    x_I(2) = xLeft  + (Ux+cSoundx)*time_simulation
    x_I(3) = xRight + (Ux-cSoundX)*time_simulation
    x_I(4) = xRight + (Ux+cSoundX)*time_simulation
    x_I(5) =           Ux         *time_simulation

    ! Reflect left going sound wave edges
    if(x_I(1) < x_I(5)) x_I(1) = 2*x_I(5) - x_I(1)
    if(x_I(3) < x_I(5)) x_I(3) = 2*x_I(5) - x_I(3)

    Distance = minval(abs(x_I - xShifted))

    if(Distance <= nI/2*dx_BLK(iBlock) + DistanceMin)then
       UserCriteria = 1.0
    else
       UserCriteria = 0.0
    endif

    ! Do not refine blocks far from discontinuity (crit=0.0)
    ! Do not coarsen blocks near discontinuity    (crit=1.0)
    RefineCritMin_I = 0.5
    CoarsenCritMax  = 0.5

    IsFound = .true.

  end subroutine user_amr_criteria

end module ModUser
