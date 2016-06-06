!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModAMR

  use ModCellGradient, ONLY: calc_gradient

  implicit none
  save

  private  ! except

  public:: init_mod_amr
  public:: clean_mod_amr
  public:: do_amr

  ! Local and global refinement/coarsening and neighbor parameters.

  ! Refinement criteria parameters
  integer, public  :: nAmrCriteria
  character(len=20), public, allocatable:: RefineCrit(:)

  real, allocatable:: AmrCriteria_IB(:,:)

  ! Refinement parameters.
  integer, public:: nRefineLevelIC
  logical, public:: DoAmr = .false.
  integer, public:: DnAmr = -1
  real,    public:: DtAmr = -1.0
  logical, public :: automatic_refinement

  integer, public:: min_block_level, max_block_level
  real,    public:: min_cell_dx, max_cell_dx
  logical, public:: DoProfileAmr = .false.

contains
  !============================================================================
  subroutine init_mod_amr(nCrit)

    use ModMain, ONLY: MaxBlock

    integer, intent(in) :: nCrit
    !-----------------------------------------------------------------------

    ! clean for each time we have new refinment criteia
    call clean_mod_amr

    allocate(RefineCrit(nCrit),              &
         AmrCriteria_IB(nCrit,MaxBlock))

  end subroutine init_mod_amr
  !============================================================================
  subroutine clean_mod_amr

    if(allocated(RefineCrit)) deallocate( &
         RefineCrit,&
         AmrCriteria_IB)

  end subroutine clean_mod_amr
  !============================================================================
  subroutine do_amr(DoFullMessagePass,TypeAmr)

    use ModProcMH
    use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,&
         lVerbose, UseB, Dt_BLK, nTrueCellsALL, &
         iNewGrid, iNewDecomposition, UseHighOrderAMR
    use ModGeometry, ONLY: CellSizeMin, CellSizeMax, true_cell
    use ModAdvance,  ONLY: DivB1_GB, iTypeAdvance_B, iTypeAdvance_BP, &
         nVar, State_VGB, &
         SkippedBlock_ !!!
    use ModLoadBalance, ONLY: load_balance
    use ModRaytrace, ONLY: ray
    use ModBlockData, ONLY: clean_block_data
    use ModIO, ONLY : write_prefix, iUnitOut
    use ModMpi

    use BATL_lib,         ONLY: regrid_batl, set_amr_criteria, &
         MaxNode, nNode, iTree_IA, nLevelMin, nLevelMax, &
         IsLogRadius, IsGenRadius, Status_, Used_, Proc_, Block_ !!!

    use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
    use ModMessagePass,   ONLY: exchange_messages
    use ModPartSteady,    ONLY: UsePartSteady
    use ModLocalTimeStep, ONLY: UseLocalTimeStep
    use ModVarIndexes, ONLY: DefaultState_V

    logical, intent(in) :: DoFullMessagePass
    character(3), intent(in) :: TypeAmr

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'do_amr'

    ! Check if we have the same grid as before, store old grid id
    integer, save :: iLastGrid=-1, iLastDecomposition=-1

    integer:: iNode !!!
    integer, allocatable:: iTypeAdvance_A(:) !!!

    !--------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Do message passing with second order accurate ghost cells

    if(DoTestMe)write(*,*)NameSub,' starts 2nd order accurate message passing'

    if(DoProfileAmr) call timing_start('amr::exchange_true')
    call exchange_messages(UseOrder2In=.true., &
         DoResChangeOnlyIn=.not.DoFullMessagePass)
    if(DoProfileAmr) call timing_stop('amr::exchange_true')

    if(UsePartSteady)then
       ! Convert iTypeAdvance_BP to _A !!! should use _A all the time
       allocate(iTypeAdvance_A(MaxNode))
       do iNode = 1, nNode
          if(iTree_IA(Status_,iNode) /= Used_) CYCLE
          iTypeAdvance_A(iNode) = &
               iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode))
       end do
    end if
    if(nAmrCriteria > 0)then
       AmrCriteria_IB(:,1:nBlockMax) = 0.0
       if(DoProfileAmr) call timing_start('amr::amr_criteria')
       call amr_criteria(AmrCriteria_IB,TypeAmr)
       if(DoProfileAmr) call timing_stop('amr::amr_criteria')
       if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
       call set_amr_criteria(nVar, State_VGB,&
            nAmrCriteria,AmrCriteria_IB,TypeAmrIn=TypeAmr)
       if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')

    else
       if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
       call set_amr_criteria(nVar, State_VGB,TypeAmrIn=TypeAmr)
       if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')
    end if

    if(DoProfileAmr) call timing_start('amr::regrid_batl')
    if(UsePartSteady)then
       call regrid_batl(nVar, State_VGB, Dt_BLK, DoTestIn=DoTestMe, &
            Used_GB=true_cell, iTypeNode_A=iTypeAdvance_A)
    else
       call regrid_batl(nVar, State_VGB, Dt_BLK, &
            DoBalanceEachLevelIn=UseLocalTimeStep, DoTestIn=DoTestMe, &
            Used_GB=true_cell, UseHighOrderAMRIn=UseHighOrderAMR, &
            DefaultStateIn_V=DefaultState_V)
    end if
    if(DoProfileAmr) call timing_stop('amr::regrid_batl')


    ! This should be eliminated by using iTypeAdvance_A everywhere !!!
    if(UsePartSteady)then
       ! restore iTypeAdvance_B and _BP
       iTypeAdvance_BP = SkippedBlock_
       iTypeAdvance_B  = SkippedBlock_
       do iNode = 1, nNode
          if(iTree_IA(Status_,iNode) /= Used_) CYCLE
          iTypeAdvance_BP(iTree_IA(Block_,iNode),iTree_IA(Proc_,iNode)) = &
               iTypeAdvance_A(iNode)
          if(iTree_IA(Proc_,iNode) == iProc) &
               iTypeAdvance_B(iTree_IA(Block_,iNode)) = iTypeAdvance_A(iNode)
       end do
       deallocate(iTypeAdvance_A)
    end if

    if(DoProfileAmr) call timing_start('amr::set_batsrus_grid')
    call set_batsrus_grid
    if(DoProfileAmr) call timing_stop('amr::set_batsrus_grid')

    ! If the grid has not changed only the message passing has to be redone
    ! to reset ghost cells at resolution changes
    if(iNewGrid==iLastGrid .and. iNewDecomposition==iLastDecomposition) then
       if(DoProfileAmr) call timing_start('amr::exchange_noamr')
       call exchange_messages(DoResChangeOnlyIn=.true., UseOrder2In=.false.)
       if(DoProfileAmr) call timing_stop('amr::exchange_noamr')
       RETURN
    end if

    iLastGrid          = iNewGrid
    iLastDecomposition = iNewDecomposition

    if(DoProfileAmr) call timing_start('amr::count_true_cells')
    call count_true_cells
    if(DoProfileAmr) call timing_stop('amr::count_true_cells')

    ! Clean all dynamically stored block data
    call clean_block_data

    if(iProc==0 .and. lVerbose>0)then
       ! Write block/cell summary after AMR
       call write_prefix; write(iUnitOut,*) '|'
       call write_prefix; write(iUnitOut,*) &
            '|  AMR:  nBlockMax = ',nBlockMax,' nBLK = ',nBLK
       call write_prefix; write(iUnitOut,*) &
            '|  AMR:  Total number of blocks used = ', nBlockALL
       call write_prefix; write(iUnitOut,*) &
            '|  AMR:  Total number of cells       = ', nBlockALL*nIJK
       call write_prefix; write(iUnitOut,*) &
            '|  AMR:  Total number of true cells  = ', nTrueCellsALL
       call write_prefix; write(iUnitOut,*) &
            '  AMR:   Min and max AMR levels      = ', nLevelMin, nLevelMax
       if(IsLogRadius .or. IsGenRadius)then
          call write_prefix; write(iUnitOut,*) &
               '  AMR:   Min and max cell size in Phi= ', &
               CellSizeMin, CellSizeMax
       else
          call write_prefix; write(iUnitOut,*) &
               '  AMR:   Min and max cell size in x/r= ', &
               CellSizeMin, CellSizeMax
       endif
       call write_prefix; write(iUnitOut,*) '|'
    end if

    if(DoProfileAmr) call timing_start('amr::set_batsrus_state')
    ! Fix energy and other variables in moved/refined/coarsened blocks
    call set_batsrus_state
    if(DoProfileAmr) call timing_stop('amr::set_batsrus_state')

    ! Update iTypeAdvance, and redo load balancing if necessary
    if(DoProfileAmr) call timing_start('amr::load_balance')
    call load_balance(DoMoveCoord=.true., DoMoveData=.true., IsNewBlock=.true.)
    if(DoProfileAmr) call timing_stop('amr::load_balance')
    ! redo message passing
    if(DoProfileAmr) call timing_start('amr::exchange_false')
    call exchange_messages(UseOrder2In=.false.)
    if(DoProfileAmr) call timing_stop('amr::exchange_false')

    ! Reset divb (it is undefined in newly created/moved blocks)
    if(UseB)then
       if(DoProfileAmr) call timing_start('amr::set_divb')
       DivB1_GB(:,:,:,1:nBlock) = -7.70

       ! write_log_file may use ray array before another ray tracing
       if(allocated(ray)) ray(:,:,:,:,:,1:nBlock) = 0.0

       if(DoProfileAmr) call timing_stop('amr::set_divb')
    end if

  end subroutine do_amr
  !============================================================================
  subroutine amr_criteria(Crit_IB,TypeAmr)

    use ModSize,       ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         x_, y_, z_, MaxBlock
    use ModMain,       ONLY: nBlock, UseB0, UseUserAmr, Unused_B,&
         DoThinCurrentSheet
    use ModGeometry,   ONLY: r_BLK, true_cell
    use ModAdvance,    ONLY: State_VGB, StateOld_VCB, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, P_
    use ModB0,         ONLY: B0_DGB
    use ModPhysics,    ONLY: rCurrents
    use ModPhysics,    ONLY: UseSunEarth
    use ModCurrent,    ONLY: get_current
    use BATL_lib,      ONLY: Xyz_DGB, CellSize_DB, masked_amr_criteria
    use ModNumConst,   ONLY: cSqrtTwo, cTiny
    use ModVarIndexes, ONLY: SignB_
    use ModUserInterface ! user_amr_criteria

    real, intent(out) :: Crit_IB(nAmrCriteria,maxBlock)
    character(len=3), intent(in) :: TypeAmr
    real :: userCriteria

    logical :: UseSwitchAMR, IsFound
    integer :: iBlock, iCrit, i, j, k
    real :: Xyz_D(3), RR, RcritAMR,AMRsort_1,AMRsort_2

    real, allocatable, save, dimension(:,:,:):: &
         Var_G, Rho_G, RhoUx_G, RhoUy_G, RhoUz_G, Bx_G, By_G, Bz_G, P_G

    ! X, Y and Z derivatives for vectors and scalars
    real, dimension(1:nI,1:nJ,1:nK) :: &
         GradXVarX_C, GradXVarY_C, GradXVarZ_C, GradX_C, &
         GradYVarX_C, GradYVarY_C, GradYVarZ_C, GradY_C, &
         GradZVarX_C, GradZVarY_C, GradZVarZ_C, GradZ_C

    real:: Current_D(3)

    character(len=*), parameter:: NameSub='amr_criteria'

    ! Needed for the 'currentsheet'
    real :: rDotB_G(nI,nJ,0:nK+1)

    !--------------------------------------------------------------------------

    if(.not.allocated(Rho_G)) allocate( &
         Var_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         Rho_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         RhoUx_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         RhoUy_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         RhoUz_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         Bx_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK),    &
         By_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK),    &
         Bz_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK),    &
         p_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK))

    ! initialize all criteria to zero
    Crit_IB = 0.0
    do iBlock = 1, nBlock
       if (Unused_B(iBlock)) CYCLE
       if (masked_amr_criteria(iBlock)) CYCLE

       ! Initialize values to use below for criteria
       if (UseSunEarth) then
          RcritAMR = 1 + 1.5*cSqrtTwo*maxval(CellSize_DB(:,iBlock))
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

       do iCrit = 1,nAmrCriteria

          if (masked_amr_criteria(iBlock,iCritExtIn=iCrit)) CYCLE
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

          case('currentsheet')
             if(SignB_>1 .and. DoThinCurrentSheet)then

                if(maxval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))>0.0.and. &
                     minval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))<0.0)then
                   Crit_IB(iCrit,iBlock) = 1.0
                else
                   Crit_IB(iCrit,iBlock) = 0.0
                end if

             else
                ! Calculate BdotR including the ghost cells in 3rd dimension only
                if(UseB0)then
                   do k=0, nK+1; do j=1, nJ; do i=1, nI
                      rDotB_G(i,j,k) = sum( Xyz_DGB(:,i,j,k,iBlock) &
                           * (B0_DGB(:,i,j,k,iBlock) &
                           +  State_VGB(Bx_:Bz_,i,j,k,iBlock)))
                   end do; end do; end do
                else
                   do k=0, nK+1; do j=1, nJ; do i=1, nI
                      rDotB_G(i,j,k) = sum(Xyz_DGB(:,i,j,k,iBlock)   &
                           * State_VGB(Bx_:Bz_,i,j,k,iBlock))
                   end do; end do; end do
                end if

                if(maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny) then
                   Crit_IB(iCrit,iBlock) = 1.0
                else
                   Crit_IB(iCrit,iBlock) = 0.0
                end if

             end if
          case default
             ! WARNING if we do not find the criteria in the above list we 
             ! will search for it among 'transient' criteria
             if (UseUserAMR .or. RefineCrit(iCrit) == 'user') then
                IsFound=.false.
                call user_amr_criteria(iBlock, &
                     UserCriteria, RefineCrit(iCrit), IsFound)
                Crit_IB(iCrit,iBlock) = userCriteria
             else
                Xyz_D = 0.5*(Xyz_DGB(:,nI,nJ,nK,iBlock) &
                     +       Xyz_DGB(:,1,1,1,iBlock))
                RR = sqrt(sum(Xyz_D**2))
                if (UseSunEarth) then
                   UseSwitchAMR = RR > RcritAMR
                else
                   UseSwitchAMR = RR > RcritAMR .and. &
                        abs(Xyz_D(z_)) <= CellSize_DB(z_,iBlock)
                end if
                if (UseSwitchAMR) then
                   ! Use dynamic refinement if there is a transient event 
                   call trace_transient(RefineCrit(iCrit),iCrit,iBlock,&
                        AMRsort_1)
                   Crit_IB(iCrit,iBlock) = AMRsort_1
                   ! Restrict the refinement to the Sun-Earth ray only
                   ! Only if UseSunEarth == .true.
                   if (UseSunEarth) then
                      call refine_sun_earth_cyl(iBlock, &
                           Xyz_D(x_), Xyz_D(y_), Xyz_D(z_), AMRsort_2)
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
    !==========================================================================

    subroutine trace_transient(NameCrit,iCrit,iBlock,refine_crit)

      character(len=*), intent(in) :: NameCrit

      integer, intent(in) :: iBlock,iCrit
      real, intent(out) :: refine_crit
      real :: AMRsort
      real, dimension(nI,nJ,nK) :: scrARR

      real, dimension(nI,nJ,nK) :: RhoOld_C, RhoUxOld_C, &
           RhoUyOld_C, RhoUzOld_C, BxOld_C, ByOld_C, BzOld_C, POld_C    
      !-----------------------------------------------------------------------

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
         scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2       + &
              RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
              sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                            + &
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
         scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2         + &
              By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2)           - &
              sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2                              + &
              ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
         scrARR(1:nI,1:nJ,1:nK) = scrARR(1:nI,1:nJ,1:nK) / max(cTiny, &
              sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
              Bz_G(1:nI,1:nJ,1:nK)**2),sqrt(BxOld_C(1:nI,1:nJ,1:nK)**2 + &
              ByOld_C(1:nI,1:nJ,1:nK)**2 + BzOld_C(1:nI,1:nJ,1:nK)**2))
         refine_crit = maxval(scrARR)
      case('meanUB') 
         scrARR(1:nI,1:nJ,1:nK) = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2        + &
              RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)       - &
              sqrt(RhoUxOld_C(1:nI,1:nJ,1:nK)**2                             + &
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
  !============================================================================
  subroutine refine_sun_earth_cone(iBlock,xBLK,yBLK,zBLK,refine_profile)

!!! The code below is WAY overcomplicated. To be removed.

    use ModMain,     ONLY: BLKtest
    use ModProcMH,   ONLY: iProc
    use ModPhysics,  ONLY: Rbody,xEarth,yEarth,zEarth,InvD2Ray
    use ModNumConst, ONLY: cRadToDeg

    ! This subroutine aims to restrict the refinement mainly along the ray 
    ! Sun-Earth, in a cone with a user-defined opening angle.

    integer, intent(in) :: iBlock
    real, intent(in) :: xBLK,yBLK,zBLK
    real, intent(out) :: refine_profile

    real :: rBLK, xxx,yyy,zzz, cutFACT
    real :: signY,cosPHI,cosTHETA
    real :: signY_BLK,cosPHI_BLK,cosTHETA_BLK
    !--------------------------------------------------------------------------

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
  !============================================================================
  subroutine refine_sun_earth_cyl(iBlock,xBLK,yBLK,zBLK,refine_profile)

!!! The code below is WAY overcomplicated. To be removed.

    use ModPhysics, ONLY: rBody, xEarth, yEarth, zEarth, InvD2Ray

    integer, intent(in) :: iBlock
    real, intent(in) :: xBLK,yBLK,zBLK
    real, intent(out) :: refine_profile
    real :: rBLK, xxx,yyy,zzz, cutFact
    real :: cosPHI,sinPHI,cosTHETA,sinTHETA
    real :: dist2BLK,yPrimeBLK

    ! This subroutine aims to restrict the refinement mainly along the ray 
    ! Sun-Earth, in a cylinder with user-defined profile across.
    !--------------------------------------------------------------------------
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

end module ModAMR
