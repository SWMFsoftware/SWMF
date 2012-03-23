!^CFG COPYRIGHT UM
subroutine amr_criteria(Crit_IB)

  use ModSize,     ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       x_, y_, z_, MaxBlock
  use ModMain,     ONLY: nBlock, UseBatl, UseB0, UseUserAmr, UnusedBlk
  use ModGeometry, ONLY: &
       x_BLK, y_BLK, z_BLK, r_BLK, dx_BLK, dy_BLK, dz_BLK, true_cell
  use ModAdvance,  ONLY: State_VGB, StateOld_VCB, B0_DGB, &
       Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, P_
  use ModAMR,      ONLY: nRefineCrit,RefineCrit
  use ModPhysics,  ONLY: rCurrents
  use ModPhysics,  ONLY: UseSunEarth
  use ModUser,     ONLY: user_amr_criteria
  use ModCurrent,  ONLY: get_current
  use BATL_lib,    ONLY: DoAmr_B, UseAmrMask
  use ModNumConst, ONLY: cSqrtTwo, cTiny
  implicit none

  real, intent(out) :: Crit_IB(4,MaxBlock)

  real :: UserCriteria

  logical :: UseSwitchAMR, IsFound
  integer :: iBlock, iCrit, iVar, i, j, k
  real :: xxx,yyy,zzz,RR, RcritAMR,AMRsort_1,AMRsort_2

  real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: &
       Var_G, Rho_G, RhoUx_G, RhoUy_G, RhoUz_G, Bx_G, By_G, Bz_G, P_G

  ! X, Y and Z derivatives for vectors and scalars
  real, dimension(1:nI,1:nJ,1:nK) :: &
       GradXVarX_C, GradXVarY_C, GradXVarZ_C, GradX_C, &
       GradYVarX_C, GradYVarY_C, GradYVarZ_C, GradY_C, &
       GradZVarX_C, GradZVarY_C, GradZVarZ_C, GradZ_C

  real:: Current_D(3)

  character(len=*), parameter:: NameSub='amr_criteria'
  !--------------------------------------------------------------------------

  ! initialize all criteria to zero
  Crit_IB = 0.0
  do iBlock = 1, nBlock
     if (unusedBLK(iBlock)) CYCLE

     if(UseBatl) then
        if(UseAmrMask) then
           if(.not.DoAmr_B(iBlock)) CYCLE
        end if
     end if

     ! set 4th criteria to block radius, used in amr_physics to preserve 
     ! symmetry (not used by BATL !!!)
     Crit_IB(4,iBlock) = maxval(R_BLK(1:nI,1:nJ,1:nK,iBlock))

     ! Initialize values to use below for criteria
     if (UseSunEarth) then
        RcritAMR = 1 + &
             1.5*cSqrtTwo*max(dx_BLK(iBlock), dy_BLK(iBlock), dz_BLK(iBlock))
     else
        RcritAMR = 0.0
     end if

     do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
        Rho_G(i,j,k)  = State_VGB(Rho_,i,j,k,iBlock)
        RhoUx_G(i,j,k)= State_VGB(RhoUx_,i,j,k,iBlock)
        RhoUy_G(i,j,k)= State_VGB(RhoUy_,i,j,k,iBlock)
        RhoUz_G(i,j,k)= State_VGB(RhoUz_,i,j,k,iBlock)
        Bx_G(i,j,k)   = State_VGB(Bx_,i,j,k,iBlock)
        By_G(i,j,k)   = State_VGB(By_,i,j,k,iBlock)
        Bz_G(i,j,k)   = State_VGB(Bz_,i,j,k,iBlock)
        P_G(i,j,k)    = State_VGB(P_,i,j,k,iBlock)
     end do; end do; end do

     do iCrit = 1, nRefineCrit
        select case(RefineCrit(iCrit))
        case('gradt')
           ! Temperature gradient.
           Var_G = P_G/Rho_G
           call calc_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
           Crit_IB(iCrit,iBlock) = &
                sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2))
        case('gradlogrho')
           ! Log of density gradient.
           Var_G = log10(Rho_G)
           call calc_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
           Crit_IB(iCrit,iBlock) = &
                sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2))
        case('gradlogp')
           ! Log of pressure gradient
           Var_G = log10(P_G)
           call calc_gradient( iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
           Crit_IB(iCrit,iBlock) = &
                sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2))

        case('gradp')
           ! Pressure gradient 2.
           call calc_gradient(iBlock, P_G, GradX_C,GradY_C,GradZ_C)
           Crit_IB(iCrit,iBlock) = &
                sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2))

        case('grade')
           ! Electric field gradient.
           if(UseB0)then
              Var_G = sqrt( &
                   ( -((RhoUy_G/Rho_G)* &
                   (Bz_G+B0_DGB(z_,:,:,:,iBlock)) - &
                   (RhoUz_G/Rho_G)* &
                   (By_G+B0_DGB(y_,:,:,:,iBlock))) )**2 + &
                   ( -((RhoUz_G/Rho_G)* &
                   (Bx_G+B0_DGB(x_,:,:,:,iBlock)) - &
                   (RhoUx_G/Rho_G)* &
                   (Bz_G+B0_DGB(z_,:,:,:,iBlock))) )**2 + &
                   ( -((RhoUx_G/Rho_G)* &
                   (By_G+B0_DGB(y_,:,:,:,iBlock)) - &
                   (RhoUy_G/Rho_G)* &
                   (Bx_G+B0_DGB(x_,:,:,:,iBlock))) )**2 )
           else
              Var_G = sqrt( &
                   (RhoUy_G*Bz_G - RhoUz_G*By_G)**2 + &
                   (RhoUz_G*Bx_G - RhoUx_G*Bz_G)**2 + &
                   (RhoUx_G*By_G - RhoUy_G*Bx_G)**2 )/Rho_G
           end if
           call calc_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
           Crit_IB(iCrit,iBlock) = &
                sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2))

        case('curlv', 'curlu')
           ! Curl of velocity
           call calc_gradient(iBlock, RhoUx_G/Rho_G, &
                GradXVarX_C, GradYVarX_C, GradZVarX_C)
           call calc_gradient(iBlock, RhoUy_G/Rho_G, &
                GradXVarY_C, GradYVarY_C, GradZVarY_C)
           call calc_gradient(iBlock, RhoUz_G/Rho_G, &
                GradXVarZ_C, GradYVarZ_C, GradZVarZ_C)
           Crit_IB(iCrit,iBlock) = sqrt(maxval( &
                (GradYVarZ_C - GradZVarY_C)**2 + &
                (GradZVarX_C - GradXVarZ_C)**2 + &
                (GradXVarY_C - GradYVarX_C)**2))

        case('curlb')
           ! Curl of magnetic field (current)
           call calc_gradient(iBlock, Bx_G, &
                GradXVarX_C, GradYVarX_C, GradZVarX_C)
           call calc_gradient(iBlock, By_G, &
                GradXVarY_C, GradYVarY_C, GradZVarY_C)
           call calc_gradient(iBlock, Bz_G, &
                GradXVarZ_C, GradYVarZ_C, GradZVarZ_C)
           Crit_IB(iCrit,iBlock) = sqrt(maxval( &
                (GradYVarZ_C - GradZVarY_C)**2 + &
                (GradZVarX_C - GradXVarZ_C)**2 + &
                (GradXVarY_C - GradYVarX_C)**2))

        case('j2')
           Crit_IB(iCrit,iBlock) = 0.0
           do k=1,nK; do j=1,nJ; do i=1,nI
              call  get_current(i, j, k, iBlock, Current_D)
              Crit_IB(iCrit,iBlock) = max(Crit_IB(iCrit,iBlock),&
                   sum(Current_D**2))
           end do;end do;end do

        case('divu','divv')
           ! Divergence of velocity (this is REALLY INEFFICIENT !!!)
           call calc_gradient( iBlock, RhoUx_G/Rho_G, &
                GradXVarX_C, GradYVarX_C, GradZVarX_C)
           call calc_gradient( iBlock, RhoUy_G/Rho_G, &
                GradXVarY_C,GradYVarY_C,GradZVarY_C)
           call calc_gradient( iBlock, RhoUz_G/Rho_G, &
                GradXVarZ_C,GradYVarZ_C,GradZVarZ_C)

           Crit_IB(iCrit,iBlock) = &
                maxval(abs(GradXVarX_C + GradYVarY_C + GradZVarZ_C))

        case('rcurrents')	
           ! Inverse distance from Rcurrents, squared
           Var_G(1:nI,1:nJ,1:nK) = 1.0/((max(cTiny, &
                abs(Rcurrents - R_BLK(1:nI,1:nJ,1:nK,iBlock))))**2)
           Crit_IB(iCrit,iBlock) = maxval(Var_G(1:nI,1:nJ,1:nK),&
                MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))

        case default
           ! WARNING if we do not find the criteria in the abou list we 
           ! will search for it among 'transient' criteria

           if (UseUserAMR .or. index(RefineCrit(iCrit),'user') > 0 ) then
              IsFound=.false.
              call user_amr_criteria(iBlock, &
                   UserCriteria, RefineCrit(iCrit), IsFound)
              if (IsFound) then
                 Crit_IB(iCrit,iBlock) = userCriteria
              else            
                 write(*,*) NameSub, &
                      ' User refinement criteria=', trim(RefineCrit(iCrit)), &
                      ' not found in user_amr_criteria in the user module!'
                 call stop_mpi('Fix user_amr_criteria or PARAM.in!')
              end if
           else
              xxx = 0.5*(x_BLK(nI,nJ,nK,iBlock)+x_BLK(1,1,1,iBlock))
              yyy = 0.5*(y_BLK(nI,nJ,nK,iBlock)+y_BLK(1,1,1,iBlock))
              zzz = 0.5*(z_BLK(nI,nJ,nK,iBlock)+z_BLK(1,1,1,iBlock))
              RR = sqrt(xxx**2+yyy**2+zzz**2 )
              if (UseSunEarth) then
                 UseSwitchAMR = RR > RcritAMR
              else
                 UseSwitchAMR = RR > RcritAMR .and. abs(zzz) <= dz_BLK(iBlock)
              end if
              if (UseSwitchAMR) then
                 !\
                 ! Use dynamic refinement if there is a transient event 
                 !/
                 call trace_transient(RefineCrit(iCrit),iCrit,iBlock,AMRsort_1)
                 Crit_IB(iCrit,iBlock) = AMRsort_1
                 !\              
                 ! Restrict the refinement to the particular ray Sun-Earth only
                 ! Only if UseSunEarth == .true.
                 !/
                 if (UseSunEarth) then
                    call refine_sun_earth_cyl(iBlock,xxx,yyy,zzz,AMRsort_2)
                 else
                    AMRsort_2 = 1.0
                 end if
                 Crit_IB(iCrit,iBlock) = AMRsort_2*Crit_IB(iCrit,iBlock)
              end if
           end if
        end select
     end do ! iCrit
  end do ! iBlock

contains
  !============================================================================
  subroutine calc_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)

    ! This is an interface to cartesian or covariant_gradient.

    use ModGeometry,ONLY: UseCovariant

    integer, intent(in):: iBlock
    real,    intent(in):: Var_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real, intent(out), dimension(1:nI,1:nJ,1:nK):: GradX_C, GradY_C, GradZ_C

    !--------------------------------------------------------------------------
    if(UseCovariant)then                               
       call covariant_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)  
    else                                                
       call cartesian_gradient(iBlock, Var_G, GradX_C, GradY_C, GradZ_C)
    end if
  end subroutine calc_gradient

  !============================================================================
  subroutine cartesian_gradient(iBlock, Var_G, GradX_C,GradY_C,GradZ_C)

    use ModGeometry, ONLY: &
         body_blk, true_cell, fAX_BLK, fAY_BLK, fAZ_BLK, vInv_CB

    integer,intent(in) :: iBlock

    real, intent(in) :: Var_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real, intent(out), dimension(1:nI,1:nJ,1:nK):: GradX_C,GradY_C,GradZ_C

    real, dimension(0:nI+1, 0:nJ+1, 0:nK+1) :: OneTrue_G

    integer :: i, j, k
    !---------------------------------------------------

    if(.not.body_blk(iBlock)) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          GradX_C(i,j,k) = 0.5*fAX_BLK(iBlock)*&
               (Var_G(i+1,j,k) - Var_G(i-1,j,k))*vInv_CB(i,j,k,iBlock)
          GradY_C(i,j,k) = 0.5*fAY_BLK(iBlock)*&
               (Var_G(i,j+1,k) - Var_G(i,j-1,k))*vInv_CB(i,j,k,iBlock)
          GradZ_C(i,j,k) = 0.5*fAZ_BLK(iBlock)*&
               (Var_G(i,j,k+1) - Var_G(i,j,k-1))*vInv_CB(i,j,k,iBlock)
       end do; end do; end do
    else
       where(true_cell(0:nI+1, 0:nJ+1, 0:nK+1,iBlock)) 
          OneTrue_G = 1.0
       elsewhere
          OneTrue_G = 0.0
       end where
       !
       !\
       ! Where .not.true_cell, all the gradients are zero
       ! In true_cell the input to gradient from the face neighbor
       ! is ignored, if the face neighbor is .not.true_cell, the input
       ! from the opposite cell is doubled in this case
       !/
       !
       do k=1,nK; do j=1,nJ; do i=1,nI
          GradX_C(i,j,k) = 0.5*fAX_BLK(iBlock)*&
               OneTrue_G(i,j,k)*(&
               (Var_G(i+1,j,k) - Var_G(i,j,k))*&
               OneTrue_G(i+1,j,k)*&
               (2.0 - OneTrue_G(i-1,j,k)) + &
               (Var_G(i,j,k)-Var_G(i-1,j,k))*&
               OneTrue_G(i-1,j,k)*&
               (2.0 - OneTrue_G(i+1,j,k)) )*vInv_CB(i,j,k,iBlock)

          GradY_C(i,j,k) = 0.5*fAY_BLK(iBlock)*&
               OneTrue_G(i,j,k)*(&
               (Var_G(i,j+1,k) - Var_G(i,j,k))*&
               OneTrue_G(i,j+1,k)*&
               (2.0 - OneTrue_G(i,j-1,k))+&
               (Var_G(i,j,k)-Var_G(i,j-1,k))*&
               OneTrue_G(i,j-1,k)*&
               (2.0 - OneTrue_G(i,j+1,k)) )*vInv_CB(i,j,k,iBlock)

          GradZ_C(i,j,k) = 0.5*fAZ_BLK(iBlock)*&
               OneTrue_G(i,j,k)*(&
               (Var_G(i,j,k+1)-Var_G(i,j,k))*&
               OneTrue_G(i,j,k+1)*&
               (2.0 - OneTrue_G(i,j,k-1))+&
               (Var_G(i,j,k)-Var_G(i,j,k-1))*&
               OneTrue_G(i,j,k-1)*&
               (2.0 - OneTrue_G(i,j,k+1)) )*vInv_CB(i,j,k,iBlock)
       end do; end do; end do
    end if
  end subroutine cartesian_gradient

  !===========================================================================

  subroutine trace_transient(NameCrit,iCrit,iBlock,refine_crit)
    use ModAMR,      ONLY:nRefineCrit,RefineCrit, nRefineLevelIC

    character(len=*), intent(in) :: NameCrit

    integer, intent(in) :: iBlock,iCrit
    real, intent(out) :: refine_crit
    real :: AMRsort
    real, dimension(1:nI,1:nJ,1:nK) :: scrARR

    real, dimension(1:nI, 1:nJ, 1:nK) :: RhoOld_C, RhoUxOld_C, &
         RhoUyOld_C, RhoUzOld_C, BxOld_C, ByOld_C, BzOld_C, POld_C    

    ! if we do a amr pysics refinmnet at startup we do not have any
    ! old walue. Can be inporved in the future
    if(nRefineLevelIC>0 ) then
       do k=1,nK; do j=1,nJ; do i=1,nI
          RhoOld_C(i,j,k)  = State_VGB(rho_,i,j,k,iBlock)
          RhoUxOld_C(i,j,k)= State_VGB(rhoUx_,i,j,k,iBlock)
          RhoUyOld_C(i,j,k)= State_VGB(rhoUy_,i,j,k,iBlock)
          RhoUzOld_C(i,j,k)= State_VGB(rhoUz_,i,j,k,iBlock)
          BxOld_C(i,j,k)   = State_VGB(Bx_,i,j,k,iBlock)
          ByOld_C(i,j,k)   = State_VGB(By_,i,j,k,iBlock)
          BzOld_C(i,j,k)   = State_VGB(Bz_,i,j,k,iBlock)
          POld_C(i,j,k)    = State_VGB(P_,i,j,k,iBlock)
       end do; end do; end do
    else
       do k=1,nK; do j=1,nJ; do i=1,nI
          RhoOld_C(i,j,k)  = StateOld_VCB(rho_,i,j,k,iBlock)
          RhoUxOld_C(i,j,k)= StateOld_VCB(rhoUx_,i,j,k,iBlock)
          RhoUyOld_C(i,j,k)= StateOld_VCB(rhoUy_,i,j,k,iBlock)
          RhoUzOld_C(i,j,k)= StateOld_VCB(rhoUz_,i,j,k,iBlock)
          BxOld_C(i,j,k)   = StateOld_VCB(Bx_,i,j,k,iBlock)
          ByOld_C(i,j,k)   = StateOld_VCB(By_,i,j,k,iBlock)
          BzOld_C(i,j,k)   = StateOld_VCB(Bz_,i,j,k,iBlock)
          POld_C(i,j,k)    = StateOld_VCB(P_,i,j,k,iBlock)
       end do; end do; end do
    end if

    select case(NameCrit)
    case('P_dot','p_dot')
       !\
       ! refine_crit = abs(|p|-|p|_o)/max(|p|,|p|_o,cTiny)
       ! over all the cells of block iBlock
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(P_G(1:nI,1:nJ,1:nK) - POld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,P_G(1:nI,1:nJ,1:nK), &
            POld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('T_dot','t_dot')
       !\
       ! refine_crit = abs(|T|-|T|_o)/max(|T|,|T|_o,cTiny)
       ! over all the cells of block iBlock
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(P_G(1:nI,1:nJ,1:nK)/Rho_G(1:nI,1:nJ,1:nK)  - &
            POld_C(1:nI,1:nJ,1:nK)/RhoOld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,P_G(1:nI,1:nJ,1:nK)/ &
            Rho_G(1:nI,1:nJ,1:nK),POld_C(1:nI,1:nJ,1:nK)/ &
            RhoOld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('Rho_dot','rho_dot')
       !\
       ! refine_crit = abs(|rho|-|rho|_o)/max(|rho|,|rho|_o,cTiny)
       ! over all the cells of block iBlock
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(Rho_G(1:nI,1:nJ,1:nK) - RhoOld_C(1:nI,1:nJ,1:nK))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny,Rho_G(1:nI,1:nJ,1:nK), &
            RhoOld_C(1:nI,1:nJ,1:nK))
       refine_crit = maxval(scrARR)
    case('RhoU_dot','rhoU_dot','rhou_dot')
       !\
       ! refine_crit = abs(|rhoU|-|rhoU|_o)/max(|rhoU|,|rhoU|_o,cTiny)
       ! over all the cells of block iBlock
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2              + &
            RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                                    + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2  + &
            RhoUz_G(1:nI,1:nJ,1:nK)**2),sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2 + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = maxval(scrARR)
    case('B_dot','b_dot')
       !\
       ! refine_crit = abs(|B|-|B|_o)/max(|B|,|B|_o,cTiny)
       ! over all the cells of block iBlock
       !/
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2           + &
            By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2                                 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
            Bz_G(1:nI,1:nJ,1:nK)**2),sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = maxval(scrARR)
    case('meanUB') 
       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2              + &
            RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                                    + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2  + &
            RhoUz_G(1:nI,1:nJ,1:nK)**2),sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2 + &
            RhoUyOld_C(1:nI,1:nJ,1:nK)**2 + RhoUzOld_C(1:nI,1:nJ,1:nK)**2))
       AMRsort = maxval(scrARR)

       scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2           + &
            By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2)      - &
            sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2                                 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
            sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
            Bz_G(1:nI,1:nJ,1:nK)**2),sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2 + &
            ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
       refine_crit = AMRsort*maxval(scrARR)
    case('Rho_2nd_1')
       scrARR(1:nI,1:nJ,1:nK) = ( & 
            abs(Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                                + &
            abs(Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                                + &
            abs(Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
            2 * Rho_G(1:nI,1:nJ,1:nK)))                               / &
            Rho_G(1:nI,1:nJ,1:nK)
       refine_crit = maxval(scrARR)
    case('Rho_2nd_2')
       scrARR(1:nI,1:nJ,1:nK) = abs( & 
            (Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                             + &
            (Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
            2 * Rho_G(1:nI,1:nJ,1:nK))                             + &
            (Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
            2 * Rho_G(1:nI,1:nJ,1:nK)))                            / &
            Rho_G(1:nI,1:nJ,1:nK)
       refine_crit = maxval(scrARR)
    case default
       call stop_mpi('Unknown RefineCrit='//NameCrit)
    end select

  end subroutine trace_transient
end subroutine amr_criteria
!==============================================================================
subroutine refine_sun_earth_cone(iBlock,xBLK,yBLK,zBLK,refine_profile)

  !!! The code below is WAY overcomplicated. To be removed.

  use ModMain,     ONLY: BLKtest
  use ModProcMH,   ONLY: iProc
  use ModPhysics,  ONLY: Rbody,xEarth,yEarth,zEarth,InvD2Ray
  use ModNumConst, ONLY: cRadToDeg

  ! This subroutine aims to restrict the refinement mainly along the ray 
  ! Sun-Earth, in a cone with a user-defined opening angle.

  implicit none
  
  integer, intent(in) :: iBlock
  real, intent(in) :: xBLK,yBLK,zBLK
  real, intent(out) :: refine_profile

  real :: rBLK, xxx,yyy,zzz, cutFACT
  real :: signY,cosPHI,cosTHETA
  real :: signY_BLK,cosPHI_BLK,cosTHETA_BLK
  !---------------------------------------------------------------------------

  cutFACT = InvD2Ray*2*cRadToDeg
  !\
  ! For InvD2Ray = 1. ==> refine_profile = 0.174587 at angle 5deg around the ray.
  ! For InvD2Ray = 2. ==> refine_profile = 0.030481 at angle 5deg around the ray.
  !/
  xxx = xEarth
  yyy = yEarth
  zzz = zEarth
  
  if (yyy == 0.0) then
     signY = 1.0
  else
     signY = abs(yyy)/yyy
  end if
  cosTHETA = zzz/sqrt(xxx**2+yyy**2+zzz**2)
  cosPHI   = xxx/sqrt(xxx**2+yyy**2)
  if ((iProc==0).and.(iBlock==BLKtest)) then
     write(*,*) ''
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<'
     write(*,*) '                 Position of the Earth'
     write(*,*) '' 
     write(*,*) 'cosPHI   =',cosPHI
     write(*,*) 'PHI      =',acos(cosPHI)*cRadToDeg
     write(*,*) 'cosTHETA =',cosTHETA
     write(*,*) 'THETA    =',acos(cosTHETA)*cRadToDeg
     write(*,*) '' 
     write(*,*) '>>>>>>>>>>>>>>>>>>>                  <<<<<<<<<<<<<<<<<<<<<'
     write(*,*) '' 
  end if
  rBLK = sqrt(xBLK**2+yBLK**2+zBLK**2)
  if (rBLK.gt.Rbody) then
     if (yBLK == 0.0) then
        signY_BLK = 1.0
     else
        signY_BLK = abs(yBLK)/yBLK
     end if
     cosTHETA_BLK = zBLK/sqrt(xBLK**2+yBLK**2+zBLK**2)
     cosPHI_BLK   = xBLK/sqrt(xBLK**2+yBLK**2)
     
     refine_profile = abs(0.5*(signY+signY_BLK)* &
          exp(-cutFACT*(acos(cosPHI_BLK)-acos(cosPHI))**2)* &
          exp(-cutFACT*(acos(cosTHETA_BLK)-acos(cosTHETA))**2))
  else
     refine_profile = 0.0
  end if
  
end subroutine refine_sun_earth_cone
!==============================================================================
subroutine refine_sun_earth_cyl(iBlock,xBLK,yBLK,zBLK,refine_profile)

  !!! The code below is WAY overcomplicated. To be removed.

  use ModPhysics, ONLY: rBody, xEarth, yEarth, zEarth, InvD2Ray

  implicit none
  
  integer, intent(in) :: iBlock
  real, intent(in) :: xBLK,yBLK,zBLK
  real, intent(out) :: refine_profile
  real :: rBLK, xxx,yyy,zzz, cutFact
  real :: cosPHI,sinPHI,cosTHETA,sinTHETA
  real :: dist2BLK,yPrimeBLK

  ! This subroutine aims to restrict the refinement mainly along the ray 
  ! Sun-Earth, in a cylinder with user-defined profile across.
  !----------------------------------------------------------------------------
  cutFact = (InvD2Ray*10)**2
  !\
  ! For InvD2Ray = 1. ==> refine_profile = 0.3679 at distance 0.1*Rsun from the ray
  ! For InvD2Ray = 2. ==> refine_profile = 0.0183 at distance 0.1*Rsun from the ray
  ! For InvD2Ray = 3. ==> refine_profile = 0.0001 at distance 0.1*Rsun from the ray
  !/


  xxx = xEarth
  yyy = yEarth
  zzz = zEarth
  cosTHETA = zzz/ &
             sqrt(xxx**2+yyy**2+zzz**2)
  sinTHETA = sqrt(xxx**2+yyy**2)/ &
             sqrt(xxx**2+yyy**2+zzz**2)
  cosPHI   = xxx/sqrt(xxx**2+yyy**2)
  sinPHI   = yyy/sqrt(xxx**2+yyy**2)

  dist2BLK  = (xBLK*sinPHI-yBLK*cosPHI)**2                + &
              (-xBLK*cosTHETA*cosPHI-yBLK*cosTHETA*sinPHI + &
              zBLK*sinTHETA)**2
  yPrimeBLK = (xBLK*sinTHETA*cosPHI+yBLK*sinTHETA*sinPHI  + &
              zBLK*cosTHETA)
  
  rBLK = sqrt(xBLK**2+yBLK**2+zBLK**2)
  if ((rBLK.gt.Rbody).and.(yPrimeBLK >= 0.0)) then
     if (cutFact*dist2BLK <= 150.0) then
        refine_profile = exp(-cutFact*dist2BLK)
     else
        refine_profile = 0.0
     end if
  else
     refine_profile = 0.0
  end if
  
end subroutine refine_sun_earth_cyl
