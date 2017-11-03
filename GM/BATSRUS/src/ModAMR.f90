!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModAMR

  use ModCellGradient, ONLY: calc_gradient

  implicit none
  save

  private  ! except

  public:: init_mod_amr    ! initialize module
  public:: clean_mod_amr   ! deallocate variables
  public:: prepare_amr     ! message pass and set AMR criteria
  public:: do_amr          ! perform AMR
  public:: set_levels      ! copy min/max_block_level into iTree_IA

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

  ! Probably iTypeAdvance_BP should be replaced by this?!
  ! This array is allocated by prepare_amr and deallocated by do_amr
  integer, allocatable:: iTypeAdvance_A(:)

contains
  !============================================================================
  subroutine init_mod_amr(nCrit)

    use ModSize, ONLY: MaxBlock

    integer, intent(in) :: nCrit
    !-----------------------------------------------------------------------

    ! clean for each time we have new refinment criteia
    call clean_mod_amr

    allocate(RefineCrit(nCrit), AmrCriteria_IB(nCrit,MaxBlock))

  end subroutine init_mod_amr
  !============================================================================
  subroutine clean_mod_amr

    if(allocated(RefineCrit)) deallocate(RefineCrit, AmrCriteria_IB)

  end subroutine clean_mod_amr
  !============================================================================
  subroutine prepare_amr(DoFullMessagePass, TypeAmr)

    use ModMain,     ONLY: nBlockMax
    use ModAdvance,  ONLY: iTypeAdvance_BP, nVar, State_VGB
    use BATL_lib,    ONLY: set_amr_criteria, &
         MaxNode, nNode, iTree_IA, &
         Status_, Used_, Proc_, Block_ !!!

    use ModMessagePass,   ONLY: exchange_messages
    use ModPartSteady,    ONLY: UsePartSteady

    logical, intent(in) :: DoFullMessagePass
    character(3), intent(in) :: TypeAmr

    integer:: iNode

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'prepare_amr'
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
       call amr_criteria(AmrCriteria_IB)
       if(DoProfileAmr) call timing_stop('amr::amr_criteria')
       if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
       call set_amr_criteria(nVar, State_VGB,&
            nAmrCriteria, AmrCriteria_IB, TypeAmrIn=TypeAmr)
       if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')
    else
       if(DoProfileAmr) call timing_start('amr::set_amr_criteria')
       call set_amr_criteria(nVar, State_VGB, TypeAmrIn=TypeAmr)
       if(DoProfileAmr) call timing_stop('amr::set_amr_criteria')
    end if

  end subroutine prepare_amr

  !============================================================================
  subroutine do_amr

    use ModProcMH
    use ModMain, ONLY : nIJK,nBLK,nBlock,nBlockMax,nBlockALL,&
         lVerbose, UseB, Dt_BLK, &
         iNewGrid, iNewDecomposition, UseHighOrderAMR, &
         UseLocalTimeStep
    use ModGeometry, ONLY: CellSizeMin, CellSizeMax, true_cell, nTrueCells, &
         count_true_cells
    use ModAdvance,  ONLY: DivB1_GB, iTypeAdvance_B, iTypeAdvance_BP, &
         nVar, State_VGB, &
         SkippedBlock_ !!!
    use ModLoadBalance, ONLY: load_balance
    use ModFieldTrace, ONLY: ray
    use ModBlockData, ONLY: clean_block_data
    use ModIO, ONLY : write_prefix, iUnitOut
    use ModMpi

    use BATL_lib,         ONLY: regrid_batl, &
         nNode, iTree_IA, nLevelMin, nLevelMax, &
         IsLogRadius, IsGenRadius, Status_, Used_, Proc_, Block_ !!!

    use ModBatlInterface, ONLY: set_batsrus_grid, set_batsrus_state
    use ModMessagePass,   ONLY: exchange_messages
    use ModPartSteady,    ONLY: UsePartSteady
    use ModVarIndexes, ONLY: DefaultState_V

    ! Check if we have the same grid as before, store old grid id
    integer:: iLastGrid=-1, iLastDecomposition=-1

    integer:: iNode

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'do_amr'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

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
            '|  AMR:  Total number of true cells  = ', nTrueCells
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
  subroutine amr_criteria(Crit_IB)

    use ModSize,       ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
         x_, y_, z_, MaxBlock
    use ModMain,       ONLY: nBlock, UseB0, UseUserAmr, Unused_B, &
         DoThinCurrentSheet
    use ModGeometry,   ONLY: r_BLK, true_cell
    use ModAdvance,    ONLY: State_VGB, StateOld_VGB, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Bx_, By_, Bz_, P_
    use ModB0,         ONLY: B0_DGB
    use ModPhysics,    ONLY: No2Io_V, UnitU_, UnitJ_, UnitP_, &
         UnitTemperature_, UnitElectric_, rCurrents
    use ModCurrent,    ONLY: get_current
    use BATL_lib,      ONLY: Xyz_DGB, masked_amr_criteria
    use ModNumConst,   ONLY: cTiny
    use ModVarIndexes, ONLY: SignB_
    use ModUserInterface ! user_amr_criteria

    real, intent(out):: Crit_IB(nAmrCriteria,MaxBlock)

    real :: UserCriteria

    logical :: IsFound
    integer :: iBlock, iCrit, i, j, k

    real, allocatable, save, dimension(:,:,:):: &
         Var_G, Rho_G, RhoUx_G, RhoUy_G, RhoUz_G, Bx_G, By_G, Bz_G, P_G

    ! X, Y and Z derivatives for vectors and scalars
    real, dimension(1:nI,1:nJ,1:nK) :: &
         GradXVarX_C, GradXVarY_C, GradXVarZ_C, GradX_C, &
         GradYVarX_C, GradYVarY_C, GradYVarZ_C, GradY_C, &
         GradZVarX_C, GradZVarY_C, GradZVarZ_C, GradZ_C

    real:: Current_D(3)

    ! Needed for the 'currentsheet'
    real :: rDotB_G(nI,nJ,0:nK+1)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='amr_criteria'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

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
                  sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2)) &
                  *No2Io_V(UnitTemperature_)
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
                  sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2)) &
                  * No2Io_V(UnitP_)

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
                  sqrt(maxval(GradX_C**2 + GradY_C**2 + GradZ_C**2)) &
                  * No2Io_V(UnitElectric_)

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
                  (GradXVarY_C - GradYVarX_C)**2)) &
                  * No2Io_V(UnitU_)

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
                  (GradXVarY_C - GradYVarX_C)**2)) &
                  *No2Io_V(UnitJ_)

          case('j2')
             Crit_IB(iCrit,iBlock) = 0.0
             do k=1,nK; do j=1,nJ; do i=1,nI
                call  get_current(i, j, k, iBlock, Current_D)
                Crit_IB(iCrit,iBlock) = max(Crit_IB(iCrit,iBlock),&
                     sum(Current_D**2))
             end do;end do;end do
             Crit_IB(iCrit,iBlock) = Crit_IB(iCrit,iBlock) &
                  *No2Io_V(UnitJ_)**2

          case('divu','divv')
             ! Divergence of velocity (this is REALLY INEFFICIENT !!!)
             call calc_gradient( iBlock, RhoUx_G/Rho_G, &
                  GradXVarX_C, GradYVarX_C, GradZVarX_C)
             call calc_gradient( iBlock, RhoUy_G/Rho_G, &
                  GradXVarY_C,GradYVarY_C,GradZVarY_C)
             call calc_gradient( iBlock, RhoUz_G/Rho_G, &
                  GradXVarZ_C,GradYVarZ_C,GradZVarZ_C)

             Crit_IB(iCrit,iBlock) = &
                  maxval(abs(GradXVarX_C + GradYVarY_C + GradZVarZ_C)) &
                  *No2Io_V(UnitU_)

          case('rcurrents')	
             ! Inverse distance from Rcurrents, squared
             ! The new geometric methods are better, but this is kept 
             ! so that the tests keep running
             Var_G(1:nI,1:nJ,1:nK) = 1.0/((max(cTiny, &
                  abs(Rcurrents - r_BLK(1:nI,1:nJ,1:nK,iBlock))))**2)
             Crit_IB(iCrit,iBlock) = maxval(Var_G(1:nI,1:nJ,1:nK),&
                  MASK=true_cell(1:nI,1:nJ,1:nK,iBlock))

          case('currentsheet')
             if(SignB_>1 .and. DoThinCurrentSheet)then
                if(maxval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))>0 .and. &
                     minval(State_VGB(SignB_,1:nI,1:nJ,0:nK+1,iBlock))<0)then
                   Crit_IB(iCrit,iBlock) = 1.0
                else
                   Crit_IB(iCrit,iBlock) = 0.0
                end if

             else
                ! Calculate R.B including the ghost cells in 3rd dimension only
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
                ! Use dynamic refinement if there is a transient event 
                call trace_transient(RefineCrit(iCrit), Crit_IB(iCrit,iBlock))
             end if
          end select
       end do ! iCrit
    end do ! iBlock

  contains
    !==========================================================================

    subroutine trace_transient(NameCrit, RefineCrit)

      character(len=*), intent(in) :: NameCrit

      real, intent(out) :: RefineCrit

      real, dimension(nI,nJ,nK) :: Tmp_C, RhoOld_C, RhoUxOld_C, &
           RhoUyOld_C, RhoUzOld_C, BxOld_C, ByOld_C, BzOld_C, POld_C    
      !-----------------------------------------------------------------------
      do k=1,nK; do j=1,nJ; do i=1,nI
         RhoOld_C(i,j,k)  = StateOld_VGB(rho_,i,j,k,iBlock)
         RhoUxOld_C(i,j,k)= StateOld_VGB(rhoUx_,i,j,k,iBlock)
         RhoUyOld_C(i,j,k)= StateOld_VGB(rhoUy_,i,j,k,iBlock)
         RhoUzOld_C(i,j,k)= StateOld_VGB(rhoUz_,i,j,k,iBlock)
         BxOld_C(i,j,k)   = StateOld_VGB(Bx_,i,j,k,iBlock)
         ByOld_C(i,j,k)   = StateOld_VGB(By_,i,j,k,iBlock)
         BzOld_C(i,j,k)   = StateOld_VGB(Bz_,i,j,k,iBlock)
         POld_C(i,j,k)    = StateOld_VGB(P_,i,j,k,iBlock)
      end do; end do; end do

      select case(NameCrit)
      case('P_dot','p_dot')
         !\
         ! RefineCrit = abs(|p|-|p|_o)/max(|p|,|p|_o,cTiny)
         ! over all the cells of block iBlock
         !/
         Tmp_C = abs(P_G(1:nI,1:nJ,1:nK) - POld_C)
         Tmp_C = Tmp_C / max(cTiny,P_G(1:nI,1:nJ,1:nK), &
              POld_C)
         RefineCrit = maxval(Tmp_C)
      case('T_dot','t_dot')
         !\
         ! RefineCrit = abs(|T|-|T|_o)/max(|T|,|T|_o,cTiny)
         ! over all the cells of block iBlock
         !/
         Tmp_C = abs(P_G(1:nI,1:nJ,1:nK)/Rho_G(1:nI,1:nJ,1:nK)  - &
              POld_C/RhoOld_C)
         Tmp_C = Tmp_C / max(cTiny,P_G(1:nI,1:nJ,1:nK)/ &
              Rho_G(1:nI,1:nJ,1:nK),POld_C/ &
              RhoOld_C)
         RefineCrit = maxval(Tmp_C)
      case('Rho_dot','rho_dot')
         !\
         ! RefineCrit = abs(|rho|-|rho|_o)/max(|rho|,|rho|_o,cTiny)
         ! over all the cells of block iBlock
         !/
         Tmp_C = abs(Rho_G(1:nI,1:nJ,1:nK) - RhoOld_C)
         Tmp_C = Tmp_C / max(cTiny,Rho_G(1:nI,1:nJ,1:nK), &
              RhoOld_C)
         RefineCrit = maxval(Tmp_C)
      case('RhoU_dot','rhoU_dot','rhou_dot')
         !\
         ! RefineCrit = abs(|rhoU|-|rhoU|_o)/max(|rhoU|,|rhoU|_o,cTiny)
         ! over all the cells of block iBlock
         !/
         Tmp_C = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2       + &
              RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2)      - &
              sqrt(RhoUxOld_C**2 + RhoUyOld_C**2 + RhoUzOld_C**2))
         Tmp_C = Tmp_C / max(cTiny,                                           &
              sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2  + &
              RhoUz_G(1:nI,1:nJ,1:nK)**2), &
              sqrt(RhoUxOld_C**2 + RhoUyOld_C**2 + RhoUzOld_C**2))
         RefineCrit = maxval(Tmp_C)
      case('B_dot','b_dot')
         !\
         ! RefineCrit = abs(|B|-|B|_o)/max(|B|,|B|_o,cTiny)
         ! over all the cells of block iBlock
         !/
         Tmp_C = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2                + &
              By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2) - &
              sqrt(BxOld_C**2 + ByOld_C**2 + BzOld_C**2))
         Tmp_C = Tmp_C / max(cTiny, &
              sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
              Bz_G(1:nI,1:nJ,1:nK)**2),                                 &
              sqrt(BxOld_C**2 + ByOld_C**2 + BzOld_C**2))
         RefineCrit = maxval(Tmp_C)
      case('meanUB')
         !  (d|rhoU|/dt)/|rhoU| * (d|B|/dt)/|B|
         Tmp_C = abs(sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2                   + &
              RhoUy_G(1:nI,1:nJ,1:nK)**2 + RhoUz_G(1:nI,1:nJ,1:nK)**2) - &
              sqrt(RhoUxOld_C**2 + RhoUyOld_C**2 + RhoUzOld_C**2))
         Tmp_C = Tmp_C / max(cTiny,                                          &
              sqrt(RhoUx_G(1:nI,1:nJ,1:nK)**2 + RhoUy_G(1:nI,1:nJ,1:nK)**2 + &
              RhoUz_G(1:nI,1:nJ,1:nK)**2),                                   &
              sqrt(RhoUxOld_C**2 + RhoUyOld_C**2 + RhoUzOld_C**2))
         RefineCrit = maxval(Tmp_C)

         Tmp_C = abs(sqrt(Bx_G(1:nI,1:nJ,1:nK)**2                + &
              By_G(1:nI,1:nJ,1:nK)**2 + Bz_G(1:nI,1:nJ,1:nK)**2) - &
              sqrt(BxOld_C**2 + ByOld_C**2 + BzOld_C**2))
         Tmp_C = Tmp_C / max(cTiny,                                     &
              sqrt(Bx_G(1:nI,1:nJ,1:nK)**2 + By_G(1:nI,1:nJ,1:nK)**2  + &
              Bz_G(1:nI,1:nJ,1:nK)**2),                                 &
              sqrt(BxOld_C**2 + ByOld_C**2 + BzOld_C**2))
         RefineCrit = RefineCrit*maxval(Tmp_C)
      case('Rho_2nd_1')
         ! (|d2Rho/dx2| + |d2Rho/dy2| + |d2Rho/dz2|)/rho
         Tmp_C = ( & 
              abs(Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
              2 * Rho_G(1:nI,1:nJ,1:nK))                            + &
              abs(Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
              2 * Rho_G(1:nI,1:nJ,1:nK))                            + &
              abs(Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
              2 * Rho_G(1:nI,1:nJ,1:nK)))                           / &
              Rho_G(1:nI,1:nJ,1:nK)
         RefineCrit = maxval(Tmp_C)
      case('Rho_2nd_2')
         ! (|d2Rho/dx2  +  d2Rho/dy2  +  d2Rho/dz2|)/rho
         Tmp_C = abs( & 
              (Rho_G(0:nI-1,1:nJ,1:nK) + Rho_G(2:nI+1,1:nJ,1:nK) - &
              2 * Rho_G(1:nI,1:nJ,1:nK))                         + &
              (Rho_G(1:nI,0:nJ-1,1:nK) + Rho_G(1:nI,2:nJ+1,1:nK) - &
              2 * Rho_G(1:nI,1:nJ,1:nK))                         + &
              (Rho_G(1:nI,1:nJ,0:nK-1) + Rho_G(1:nI,1:nJ,2:nK+1) - &
              2 * Rho_G(1:nI,1:nJ,1:nK)))                        / &
              Rho_G(1:nI,1:nJ,1:nK)
         RefineCrit = maxval(Tmp_C)
      case default
         call stop_mpi('Unknown RefineCrit='//NameCrit)
      end select

    end subroutine trace_transient

  end subroutine amr_criteria
  !============================================================================
  subroutine set_levels

    use BATL_lib, ONLY: iTree_IA, MinLevel_, MaxLevel_

    iTree_IA(MinLevel_,:) = min_block_level
    iTree_IA(MaxLevel_,:) = max_block_level

  end subroutine set_levels
  !=======================================================================

end module ModAMR
