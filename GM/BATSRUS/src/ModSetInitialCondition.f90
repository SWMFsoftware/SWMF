!^CFG COPYRIGHT UM
subroutine set_ICs
  use ModMain
  use ModAdvance
  use ModGeometry, ONLY : x2,y2,z2,x_BLK,y_BLK,z_BLK,R_BLK,true_cell
  use ModIO, ONLY : restart
  use ModImplicit,ONLY: UsePartImplicit             !^CFG IF IMPLICIT
  use ModPhysics
  use ModNumConst
  implicit none

  real :: Rmax, SinSlope, CosSlope
  real :: B4, dB4dx, zeta4, q4, epsi4, plobe, &
       XFace, YFace, ZFace
  integer::i,j,k,iVar

  B0xCell_BLK(:,:,:,globalBLK) = 0.00
  B0yCell_BLK(:,:,:,globalBLK) = 0.00
  B0zCell_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_x_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_y_BLK(:,:,:,globalBLK) = 0.00
  B0xFace_z_BLK(:,:,:,globalBLK) = 0.00
  B0yFace_z_BLK(:,:,:,globalBLK) = 0.00
  B0zFace_z_BLK(:,:,:,globalBLK) = 0.00

  time_BLK(:,:,:,globalBLK) = 0.00

  fbody_x_BLK(:,:,:,globalBLK) = 0.00
  fbody_y_BLK(:,:,:,globalBLK) = 0.00
  fbody_z_BLK(:,:,:,globalBLK) = 0.00
  qheat_BLK(:,:,:,globalBLK) = 0.00

  Flux_VX = cZero
  Flux_VY = cZero
  Flux_VZ = cZero




  if(UsePartImplicit)&                              !^CFG IF IMPLICIT
       call init_conservative_facefluxes(globalBLK) !^CFG IF IMPLICIT

  if(unusedBLK(globalBLK))then  
     do iVar=1,nVar
        State_VGB(iVar,:,:,:,globalBLK)   = DefaultState_V(iVar)
     end do
  else


     !\
     ! If used, initialize solution variables and parameters.
     !/
     select case (problem_type)                !^CFG IF NOT SIMPLE BEGIN
     case(problem_dissipation)                    !^CFG IF DISSFLUX
        if(.not.restart) call set_ICs_diss_test   !^CFG IF DISSFLUX
     case (problem_uniform)
        !\
        ! Initialize solution variables for
        ! a uniform flow problem.
        !/
        if(.not.restart)then
           State_VGB(Bz_,:,:,:,globalBLK) = -1.0
           State_VGB(P_,:,:,:,globalBLK)  =  0.6
        endif

     case (problem_shocktube)

        !\
        ! Initialize solution variables for
        ! the (rotated) shock tube problem.
        !/
        if(.not.restart)then
           SinSlope=ShockSlope/sqrt(cOne+ShockSlope**2)
           CosSlope=      cOne/sqrt(cOne+ShockSlope**2)
           where (X_BLK(:,:,:,globalBLK) < -ShockSlope*Y_BLK(:,:,:,globalBLK))
              State_VGB(rho_,:,:,:,globalBLK)   =  shock_Lstate(rho_)
              State_VGB(rhoUx_,:,:,:,globalBLK) =  shock_Lstate(rho_)* &
                   (CosSlope*shock_Lstate(ux_)-SinSlope*shock_Lstate(uy_))
              State_VGB(rhoUy_,:,:,:,globalBLK) =  shock_Lstate(rho_)* &
                   (SinSlope*shock_Lstate(ux_)+CosSlope*shock_Lstate(uy_))
              State_VGB(rhoUz_,:,:,:,globalBLK) = shock_Lstate(rho_)*shock_Lstate(uz_)
              State_VGB(Bx_,:,:,:,globalBLK) = &
                   CosSlope*shock_Lstate(bx_)-SinSlope*shock_Lstate(by_)
              State_VGB(By_,:,:,:,globalBLK) = &
                   SinSlope*shock_Lstate(bx_)+CosSlope*shock_Lstate(by_)
              State_VGB(Bz_,:,:,:,globalBLK) =  shock_Lstate(bz_)
              State_VGB(P_,:,:,:,globalBLK)  =  shock_Lstate(P_)
           elsewhere
              State_VGB(rho_,:,:,:,globalBLK)   =  shock_Rstate(rho_)
              State_VGB(rhoUx_,:,:,:,globalBLK) =  shock_Rstate(rho_)* &
                   (CosSlope*shock_Rstate(ux_)-SinSlope*shock_Rstate(uy_))
              State_VGB(rhoUy_,:,:,:,globalBLK) =  shock_Rstate(rho_)* &
                   (SinSlope*shock_Rstate(ux_)+CosSlope*shock_Rstate(uy_))
              State_VGB(rhoUz_,:,:,:,globalBLK) = shock_Rstate(rho_)*shock_Rstate(uz_)
              State_VGB(Bx_,:,:,:,globalBLK) = &
                   CosSlope*shock_Rstate(bx_)-SinSlope*shock_Rstate(by_)
              State_VGB(By_,:,:,:,globalBLK) = &
                   SinSlope*shock_Rstate(bx_)+CosSlope*shock_Rstate(by_)
              State_VGB(Bz_,:,:,:,globalBLK) =  shock_Rstate(bz_)
              State_VGB(P_,:,:,:,globalBLK)  =  shock_Rstate(P_)
           end where
        end if ! not restart

        if(index(test_string,'testBsplit')>0)then
           write(*,*)'Splitting B between B0 and B1 to test the code!!!'
           State_VGB(Bx_,:,:,:,globalBLK)        = State_VGB(Bx_,:,:,:,globalBLK) &
                -shock_Lstate(bx_)*0.7
           B0xCell_BLK(:,:,:,globalBLK)   = shock_Lstate(bx_)*0.7
           B0xFace_x_BLK(:,:,:,globalBLK) = shock_Lstate(bx_)*0.7
           B0xFace_y_BLK(:,:,:,globalBLK) = shock_Lstate(bx_)*0.7
           B0xFace_z_BLK(:,:,:,globalBLK) = shock_Lstate(bx_)*0.7


           State_VGB(By_,:,:,:,globalBLK)        = State_VGB(By_,:,:,:,globalBLK) &
                -shock_Lstate(by_)*0.4
           B0yCell_BLK(:,:,:,globalBLK)   = shock_Lstate(by_)*0.4
           B0yFace_x_BLK(:,:,:,globalBLK) = shock_Lstate(by_)*0.4
           B0yFace_y_BLK(:,:,:,globalBLK) = shock_Lstate(by_)*0.4
           B0yFace_z_BLK(:,:,:,globalBLK) = shock_Lstate(by_)*0.4

           State_VGB(Bz_,:,:,:,globalBLK)        = State_VGB(Bz_,:,:,:,globalBLK) &
                -shock_Lstate(bz_)*0.2
           B0zCell_BLK(:,:,:,globalBLK)   = shock_Lstate(bz_)*0.2
           B0zFace_x_BLK(:,:,:,globalBLK) = shock_Lstate(bz_)*0.2
           B0zFace_y_BLK(:,:,:,globalBLK) = shock_Lstate(bz_)*0.2
           B0zFace_z_BLK(:,:,:,globalBLK) = shock_Lstate(bz_)*0.2
        endif
     case (problem_heliosphere)
        !\
        ! Initialize solution variables for
        ! a heliospheric flow problem.
        !/

        !\
        ! Calculate intrinsic field values (Heliosphere).
        !/
        call set_b0(globalBLK)

        !\
        ! Calculate the volume averaged body forces due
        ! to gravivational acceleration (Heliosphere).
        !/
        call body_force_averages

        !\
        ! Calculate the volume averaged heating sources
        ! (Heliosphere).
        !/
        if(UseUserHeating)     call user_heat_source           !^CFG IF USERFILES
        if(.not.restart)then
           !\
           ! Initialize solution quantities (Heliosphere).
           !/
           Rmax = max(21.00, sqrt(x2**2+y2**2+z2**2))
           if (UseInertial) then 
              where (R_BLK(:,:,:,globalBLK) > rBody)
                 State_VGB(rho_,:,:,:,globalBLK) = 1.00/R_BLK(:,:,:,globalBLK)**3
                 State_VGB(P_,:,:,:,globalBLK) = inv_g*State_VGB(rho_,:,:,:,globalBLK)
                 State_VGB(rhoUx_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* &
                      (x_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)) 
                 State_VGB(rhoUy_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* &
                      (y_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)) 
                 State_VGB(rhoUz_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* &
                      (z_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK))
                 State_VGB(Bx_,:,:,:,globalBLK)=0.00
                 State_VGB(By_,:,:,:,globalBLK)=0.00
                 State_VGB(Bz_,:,:,:,globalBLK)=0.00
              elsewhere
                 State_VGB(rho_,:,:,:,globalBLK) = 1.00
                 State_VGB(P_,:,:,:,globalBLK) = inv_g
                 State_VGB(Bx_,:,:,:,globalBLK)=0.00
                 State_VGB(By_,:,:,:,globalBLK)=0.00
                 State_VGB(Bz_,:,:,:,globalBLK)=0.00
              end where
           else
              where (R_BLK(:,:,:,globalBLK) > rBody)
                 State_VGB(rho_,:,:,:,globalBLK) = 1.00/R_BLK(:,:,:,globalBLK)**3
                 State_VGB(P_,:,:,:,globalBLK) = inv_g*State_VGB(rho_,:,:,:,globalBLK)
                 State_VGB(rhoUx_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      (40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* &
                      (x_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)) + &
                      OMEGAbody*y_BLK(:,:,:,globalBLK))
                 State_VGB(rhoUy_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      (40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* &
                      (y_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)) - &
                      OMEGAbody*x_BLK(:,:,:,globalBLK))
                 State_VGB(rhoUz_,:,:,:,globalBLK) = State_VGB(rho_,:,:,:,globalBLK)* &
                      40.0*((R_BLK(:,:,:,globalBLK)-1.)/(Rmax-1.))* & 
                      (z_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK))
                 State_VGB(Bx_,:,:,:,globalBLK)=0.00
                 State_VGB(By_,:,:,:,globalBLK)=0.00
                 State_VGB(Bz_,:,:,:,globalBLK)=0.00
              elsewhere
                 State_VGB(rho_,:,:,:,globalBLK) = 1.00
                 State_VGB(P_,:,:,:,globalBLK) = inv_g
                 State_VGB(rhoUx_,:,:,:,globalBLK) = 0.00
                 State_VGB(rhoUy_,:,:,:,globalBLK) = 0.00
                 State_VGB(rhoUz_,:,:,:,globalBLK) = 0.00
                 State_VGB(Bx_,:,:,:,globalBLK)=0.00
                 State_VGB(By_,:,:,:,globalBLK)=0.00
                 State_VGB(Bz_,:,:,:,globalBLK)=0.00
              end where
           endif
        end if
     case (problem_arcade)

        if(.not.restart)call set_ICs_arc

        call body_force_averages

     case (problem_cme)

        if(.not.restart)call set_ICs_cme

        call body_force_averages

     case (problem_comet)
        !\
        ! Initialize solution variables for
        ! comet problem.
        !/

        if(.not.restart)then
           do iVar=1,nVar
              State_VGB(iVar,:,:,:,globalBLK)   = CellState_VI(iVar,1)
           end do
        end if

     case (problem_rotation)
        !\
        ! Initialize solution variables for
        ! Rotating Dipole Test Case
        !/
        !\
        ! Calculate intrinsic field values (Magnetosphere).
        !/
        call set_b0(globalBLK)

        !\
        ! Calculate the volume averaged body forces due
        ! to gravivational acceleration (Magnetosphere).
        !/
        call body_force_averages

        if(.not.restart)then


           where (Z_BLK(:,:,:,globalBLK)<=.666667*X_BLK(:,:,:,globalBLK) .and. &
                Z_BLK(:,:,:,globalBLK)>=.333332*X_BLK(:,:,:,globalBLK) .and. &
                Y_BLK(:,:,:,globalBLK)<=.666667*X_BLK(:,:,:,globalBLK) .and. &
                Y_BLK(:,:,:,globalBLK)>=.333332*X_BLK(:,:,:,globalBLK) .and. &
                X_BLK(:,:,:,globalBLK)>=0.0)
              State_VGB(rho_,:,:,:,globalBLK)   = SW_rho
              State_VGB(P_,:,:,:,globalBLK)     = SW_p
              State_VGB(rhoUx_,:,:,:,globalBLK) = SW_rho*SW_Ux* &
                   X_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(rhoUy_,:,:,:,globalBLK) = SW_rho*SW_Ux* &
                   Y_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(rhoUz_,:,:,:,globalBLK) = SW_rho*SW_Ux* &
                   Z_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(Bx_,:,:,:,globalBLK)    = SW_Bx
              State_VGB(By_,:,:,:,globalBLK)    = SW_By
              State_VGB(Bz_,:,:,:,globalBLK)    = SW_Bz
           elsewhere 
              State_VGB(rho_,:,:,:,globalBLK)   = 0.000000000001*SW_rho
              State_VGB(P_,:,:,:,globalBLK)     = SW_p
              State_VGB(rhoUx_,:,:,:,globalBLK) = 0.000000000001*SW_rho*SW_Ux* &
                   X_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(rhoUy_,:,:,:,globalBLK) = 0.000000000001*SW_rho*SW_Ux* &
                   Y_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(rhoUz_,:,:,:,globalBLK) = 0.000000000001*SW_rho*SW_Ux* &
                   Z_BLK(:,:,:,globalBLK)/R_BLK(:,:,:,globalBLK)
              State_VGB(Bx_,:,:,:,globalBLK)    = SW_Bx
              State_VGB(By_,:,:,:,globalBLK)    = SW_By
              State_VGB(Bz_,:,:,:,globalBLK)    = SW_Bz
           end where
        end if

     case (problem_diffusion)
        !\
        ! Initialize solution variables for
        ! a magnetic diffusion test case problem.
        !/
        if(.not.restart)then
           State_VGB(rho_,:,:,:,globalBLK)   =  SW_rho
           State_VGB(rhoUx_,:,:,:,globalBLK) =  0.0
           State_VGB(rhoUy_,:,:,:,globalBLK) =  0.0
           State_VGB(rhoUz_,:,:,:,globalBLK) =  SW_rho*SW_Uz
           State_VGB(By_,:,:,:,globalBLK) =  0.0
           State_VGB(Bz_,:,:,:,globalBLK) =  0.0
           State_VGB(P_,:,:,:,globalBLK)  =  inv_g

           where (Z_BLK(:,:,:,globalBLK) > 0.0)
              State_VGB(Bx_,:,:,:,globalBLK) =  SW_Bx
           elsewhere
              State_VGB(Bx_,:,:,:,globalBLK) =  -SW_Bx
           end where
        end if

     case (problem_earth)                       !^CFG END SIMPLE
        !\
        ! Initialize solution variables for
        ! Earth Magnetosphere flow problem.
        !/

        !\
        ! Calculate intrinsic field values (Magnetosphere).
        !/
        call set_b0(globalBLK)
        call body_force_averages

        if(.not.restart)then
           !\
           ! Initialize solution quantities (Magnetosphere).
           !/
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   = CellState_VI(:,body1_)
              end if
           end do;end do;end do
        end if

     case (problem_saturn,problem_jupiter)    !^CFG IF NOT SIMPLE BEGIN
        !\
        ! Initialize solution variables for
        ! Saturn and Jupiter Magnetosphere flow problem.
        !/
        !\
        ! Calculate intrinsic field values (Magnetosphere).
        !/
        call set_b0(globalBLK)

        !\
        ! Calculate the volume averaged body forces due
        ! to gravivational acceleration (Magnetosphere).
        !/
        call body_force_averages

        if(.not.restart)then
           !\
           ! Initialize solution quantities (Magnetosphere).
           !/
           !\
           ! Initialize solution quantities (Magnetosphere).
           !/
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
              end if
           end do;end do;end do
        end if

     case (problem_venus)
        !\
        ! Initialize solution variables for
        ! Venus Ionosphere flow problem.
        !/
        !\
        ! Calculate intrinsic field values.
        !/
        call set_b0(globalBLK)

        !\
        ! Calculate the volume averaged body forces due
        ! to gravivational acceleration.
        !/
        call body_force_averages

        if(.not.restart)then
           !\
           ! Initialize solution quantities.
           !/
          
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
              end if
           end do;end do;end do
        endif

     case (problem_mars)
        !\
        ! Initialize solution variables for
        ! Mars Ionosphere flow problem.
        !/
        !\
        ! Calculate intrinsic field values.
        !/
        call set_b0(globalBLK)

        !\
        ! Calculate the volume averaged body forces due
        ! to gravivational acceleration.
        !/
        call body_force_averages

        if(.not.restart)then
           !\
           ! Initialize solution quantities.
           !/
          
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
              end if
           end do;end do;end do
        endif
     case (problem_cylinder)
        !\
        ! Initialize solution variables for
        ! 2 Dimensional mhd shock problem.
        !/

        if(.not.restart)then
           !\
           ! Initialize solution quantities.
           !/
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
              end if
           end do;end do;end do
        end if

     case (problem_sphere)
        !\
        ! Initialize solution variables for
        ! 3 Dimensional mhd shock problem.
        !/

        if(.not.restart)then
           !\
           ! Initialize solution quantities.
           !/
           do k=1-gcn,nK+gcn;do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
              if(true_cell(i,j,k,globalBLK))then
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,1)
              else
                 State_VGB(:,i,j,k,globalBLK)   =  CellState_VI(:,body1_)
              end if
           end do;end do;end do
        end if

     end select             !^CFG END SIMPLE             
     if(UseUserICs) call user_set_ICs       !^CFG IF USERFILES
  end if ! unusedBLK
!^CFG IF NOT SIMPLE BEGIN
  ! Eliminate B1 around body if necessary
  if(.not.restart)then
     if(UseConstrainB)then  !^CFG IF CONSTRAINB BEGIN
        call constrain_ICs
     else                   !^CFG END CONSTRAINB
        if(index(test_string,'DriftIn')>0)then
           where(x_BLK(:,:,:,globalBLK)<10.)
              State_VGB(Bx_,:,:,:,globalBLK)=0.0
              State_VGB(By_,:,:,:,globalBLK)=0.0
              State_VGB(Bz_,:,:,:,globalBLK)=0.0
              State_VGB(P_,:,:,:,globalBLK)=&
                   State_VGB(P_,:,:,:,globalBLK)+ &
                   0.5*(SW_Bx**2+SW_By**2+SW_Bz**2)
           elsewhere
              State_VGB(Bx_,:,:,:,globalBLK)=SW_Bx
              State_VGB(By_,:,:,:,globalBLK)=SW_By
              State_VGB(Bz_,:,:,:,globalBLK)=SW_Bz
           end where
        end if
     end if                 !^CFG IF CONSTRAINB
  end if
!^CFG END SIMPLE 
  !\
  ! Compute energy from set values above.
  !/
  call correctE
end subroutine set_ICs


!^CFG IF NOT SIMPLE BEGIN
subroutine coronal_hole_boundary(RadiusSun, sin2Theta_coronal_hole)
  implicit none

  real, intent(in) :: RadiusSun
  real, intent(out) :: sin2Theta_coronal_hole 

  ! Determine coronal hole boundary.

  !====== Milestone 2 Coronal hole parameters ==========
  !  sin2Theta_coronal_hole = 0.250000    ! 30 degrees
  ! sin2Theta_coronal_hole = 0.1786062   ! 25 degrees
  !  sin2Theta_coronal_hole = 0.116980    ! 20 degrees
  !  sin2Theta_coronal_hole = 1.000000    ! 90 degrees
  !=====================================================

  !====== Milestone 3 Coronal hole parameters ==========
  sin2Theta_coronal_hole = 0.090423978 ! 17.5 degrees
  !=====================================================

  !  sin2Theta_coronal_hole = 1.000000    ! 90 degrees
  !  sin2Theta_coronal_hole = 0.32898993  ! 35 degrees
  !  sin2Theta_coronal_hole = 0.250000    ! 30 degrees
  !  sin2Theta_coronal_hole = 0.1786062   ! 25 degrees
  !  sin2Theta_coronal_hole = 0.14644661  ! 22.5 degrees
  !  sin2Theta_coronal_hole = 0.116980    ! 20 degrees
  !  sin2Theta_coronal_hole = 0.090423978 ! 17.5 degrees
  !  sin2Theta_coronal_hole = 0.066987    ! 15 degrees
  !  sin2Theta_coronal_hole = 0.030154    ! 10 degrees

  if (RadiusSun > 1.00 .and. RadiusSun <= 7.00) then
     sin2Theta_coronal_hole = min(sin2Theta_coronal_hole + &
          (1.00-sin2Theta_coronal_hole)* &
          (RadiusSun-1.00)/(9.00-1.00), &
          1.00)
  else if (RadiusSun > 7.00) then
     sin2Theta_coronal_hole = min(sin2Theta_coronal_hole + &
          (1.00-sin2Theta_coronal_hole)* &
          (7.00-1.00)/(9.00-1.00), &
          1.00)
     sin2Theta_coronal_hole = min(sin2Theta_coronal_hole + &
          (1.00-sin2Theta_coronal_hole)* &
          (RadiusSun-7.00)/(47.00-7.00), &
          1.00)
  end if

end subroutine coronal_hole_boundary
!^CFG END SIMPLE 
