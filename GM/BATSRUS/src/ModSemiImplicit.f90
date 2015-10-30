!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModSemiImplicit

  use ModSemiImplVar
  use ModProcMH,   ONLY: iProc
  use ModImplicit, ONLY: nStencil
  use ModLinearSolver, ONLY: LinearSolverParamType

  implicit none
  save

  private ! except

  public:: read_semi_impl_param  ! Read parameters
  public:: init_mod_semi_impl    ! Initialize variables
  public:: clean_mod_semi_impl   ! Deallocate variables
  public:: advance_semi_impl     ! Advance semi implicit operator

  ! The index of the variable currently being solved for
  integer, public:: iVarSemi

  ! These arrays with *Semi_* and *SemiAll_* are indexed by compact iBlockSemi
  ! and have a single ghost cell at most.
  ! The SemiAll_ variables are indexed from 1..nVarSemiAll
  real, allocatable:: DconsDsemiAll_VCB(:,:,:,:,:) ! dCons/dSemi derivatives
  real, allocatable:: SemiAll_VCB(:,:,:,:,:)       ! Semi-implicit vars
  real, allocatable:: NewSemiAll_VCB(:,:,:,:,:)    ! Updated semi-impl vars
  real, allocatable:: ResSemi_VCB(:,:,:,:,:)       ! Result of Matrix(Semi)
  real, allocatable:: JacSemi_VVCIB(:,:,:,:,:,:,:) ! Jacobian/preconditioner

  ! Store Difference of U^* (after explicit update)  and U^n. 
  real, allocatable:: DeltaSemiAll_VCB(:,:,:,:,:) 

  ! Linear arrays for RHS, unknowns, pointwise Jacobi preconditioner
  real, allocatable:: Rhs_I(:), x_I(:), JacobiPrec_I(:)

  ! Default parameters for the semi-implicit linear solver
  type(LinearSolverParamType):: SemiParam = LinearSolverParamType( &
       .true.,      &! DoPrecond
       'left',      &! TypePrecondSide
       'BILU',      &! TypePrecond
       0.0,         &! PrecondParam (Gustafsson for MBILU)
       'GMRES',     &! TypeKrylov
       'rel',       &! TypeStop
       1e-3,        &! ErrorMax
       100,         &! MaxMatvec
       100,         &! nKrylovVector
       .false.,     &! UseInitialGuess
       -1.0, -1, -1) ! Error, nMatvec, iError (return values)

  ! Index of the test block
  integer:: iBlockSemiTest = 1

contains
  !============================================================================
  subroutine read_semi_impl_param(NameCommand)

    use ModReadParam,     ONLY: read_var
    use ModLinearSolver,  ONLY: &
         Jacobi_, BlockJacobi_, GaussSeidel_, Dilu_, Bilu_

    character(len=*), intent(in) :: NameCommand

    integer :: i
    character(len=*), parameter:: NameSub = 'read_semi_impl_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#SEMIIMPLICIT', '#SEMIIMPL')
       call read_var('UseSemiImplicit', UseSemiImplicit)
       if(UseSemiImplicit)then
          call read_var('TypeSemiImplicit', TypeSemiImplicit)
          i = index(TypeSemiImplicit,'split')
          UseSplitSemiImplicit = i > 0
          if(UseSplitSemiImplicit) &
               TypeSemiImplicit = TypeSemiImplicit(1:i-1)
       end if


    case("#SEMIIMPLICITSTABLE")
       call read_var('UseStableImplicit',UseStableImplicit)
       
    case("#SEMICOEFF", "#SEMIIMPLCOEFF", "#SEMIIMPLICITCOEFF")
       call read_var('SemiImplCoeff', SemiImplCoeff)

    case("#SEMIPRECONDITIONER", "#SEMIPRECOND")
       call read_var('DoPrecond',   SemiParam%DoPrecond)
       call read_var('TypePrecond', SemiParam%TypePrecond, IsUpperCase=.true.)
       select case(SemiParam%TypePrecond)
       case('HYPRE')
       case('JACOBI')
          SemiParam%PrecondParam = Jacobi_
       case('BLOCKJACOBI')
          SemiParam%PrecondParam = BlockJacobi_
       case('GS')
          SemiParam%PrecondParam = GaussSeidel_
       case('DILU')
          SemiParam%PrecondParam = Dilu_
       case('BILU')
          SemiParam%PrecondParam = Bilu_
       case('MBILU')
          call read_var('GustafssonPar', SemiParam%PrecondParam)
          SemiParam%PrecondParam = -SemiParam%PrecondParam
       case default
          call stop_mpi(NameSub//' invalid TypePrecond='// &
               SemiParam%TypePrecond)
       end select
    case("#SEMIKRYLOV")
       call read_var('TypeKrylov', SemiParam%TypeKrylov, IsUpperCase=.true.)
       call read_var('ErrorMax', SemiParam%ErrorMax)
       call read_var('MaxMatvec', SemiParam%MaxMatvec)
       SemiParam%nKrylovVector = SemiParam%MaxMatvec

    case('#SEMIKRYLOVSIZE')
       call read_var('nKrylovVector', SemiParam%nKrylovVector)
    case default
       call stop_mpi(NameSub//' invalid NameCommand='//NameCommand)
    end select

    if(UseSemiImplicit .and. UseStableImplicit .and. &
         (TypeSemiImplicit /= 'cond' .and. TypeSemiImplicit /='parcond')) &
         call stop_mpi(NameSub//' StableImplicit is not available for '&
         //TypeSemiImplicit)
    
  end subroutine read_semi_impl_param
  !============================================================================
  subroutine init_mod_semi_impl

    use ModLinearSolver, ONLY: bilu_
    use BATL_lib, ONLY: MaxDim, nI, nJ, nK, nIJK, &
         MinI, MaxI, MinJ, MaxJ, MinK, MaxK, MaxBlock
    use ModVarIndexes, ONLY: nWave
    use ModResistivity, ONLY: BxImpl_

    character(len=*), parameter:: NameSub = 'init_mod_semi_impl'
    !------------------------------------------------------------------------
    if(allocated(SemiAll_VCB)) RETURN

    allocate(iBlockFromSemi_B(MaxBlock))

    nVectorSemi = 0
    if(allocated(iVectorSemi_I)) deallocate(iVectorSemi_I)

    select case(TypeSemiImplicit)
    case('radiation')
       if(nWave == 1 .or. UseSplitSemiImplicit)then
          ! Radiative transfer with point implicit temperature: I_w
          nVarSemiAll = nWave
       else
          ! Radiative transfer: (electron) temperature and waves
          nVarSemiAll = 1 + nWave
       end if
    case('cond', 'parcond')
       ! Heat conduction: T
       nVarSemiAll = 1
    case('radcond')
       ! Radiative transfer and heat conduction: I_w and T
       nVarSemiAll = nWave + 1
    case('resistivity')
       ! (Hall) resistivity: magnetic field
       nVarSemiAll = MaxDim
       nVectorSemi = 1
       allocate(iVectorSemi_I(nVectorSemi))
       iVectorSemi_I(1) = BxImpl_
    case default
       call stop_mpi(NameSub//': nVarSemiAll unknown for'//TypeSemiImplicit)
    end select

    if(UseSplitSemiImplicit)then
       ! Split semi-implicit scheme solves 1 implicit variable at a time
       nVarSemi = 1
    else
       ! Unsplit (semi-)implicit scheme solves all impl. variables together
       nVarSemi = nVarSemiAll
    end if
    ! Set range of (semi-)implicit variables for unsplit scheme.
    ! For split scheme this will be overwritten.
    iVarSemiMin = 1
    iVarSemiMax = nVarSemi

    ! Fix preconditioner type from DILU to BILU for a single variable
    ! because BILU is optimized for scalar equation.
    if(nVarSemi == 1 .and. SemiParam%TypePrecond == 'DILU')then
       SemiParam%TypePrecond  = 'BILU'
       SemiParam%PrecondParam = Bilu_
    end if

    if(nVarSemi > 1 .and. SemiParam%TypePrecond == 'HYPRE' .and. iProc==0) &
         call stop_mpi( &
         'HYPRE preconditioner requires split semi-implicit scheme!')

    allocate(SemiAll_VCB(nVarSemiAll,nI,nJ,nK,MaxBlock))
    allocate(NewSemiAll_VCB(nVarSemiAll,nI,nJ,nK,MaxBlock))
    allocate(DconsDsemiAll_VCB(nVarSemiAll,nI,nJ,nK,MaxBlock))

    allocate(SemiState_VGB(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    allocate(ResSemi_VCB(nVarSemi,nI,nJ,nK,MaxBlock))

    allocate(JacSemi_VVCIB(nVarSemi,nVarSemi,nI,nJ,nK,nStencil,MaxBlock))

    allocate(DeltaSemiAll_VCB(nVarSemiAll,nI,nJ,nK,MaxBlock))

    ! Variables for the flux correction
    if(  TypeSemiImplicit == 'parcond' .or. &
         TypeSemiImplicit == 'resistivity' .or. UseAccurateRadiation)then
       allocate( &
            FluxImpl_VXB(nVarSemi,nJ,nK,2,MaxBlock), &
            FluxImpl_VYB(nVarSemi,nI,nK,2,MaxBlock), &
            FluxImpl_VZB(nVarSemi,nI,nJ,2,MaxBlock) )
       FluxImpl_VXB = 0.0
       FluxImpl_VYB = 0.0
       FluxImpl_VZB = 0.0
    end if

    ! Linear arrays
    allocate(Rhs_I(nVarSemi*nIJK*MaxBlock))
    allocate(x_I(nVarSemi*nIJK*MaxBlock))
    if(SemiParam%TypePrecond == 'JACOBI')then
       allocate(JacobiPrec_I(nVarSemi*nIJK*MaxBlock))
    else
       allocate(JacobiPrec_I(1))
    end if

  end subroutine init_mod_semi_impl
  !============================================================================
  subroutine clean_mod_semi_impl

    ! Deallocate all variables

    if(allocated(iBlockFromSemi_B))  deallocate(iBlockFromSemi_B)
    if(allocated(SemiAll_VCB))       deallocate(SemiAll_VCB)
    if(allocated(NewSemiAll_VCB))    deallocate(NewSemiAll_VCB)
    if(allocated(DeltaSemiAll_VCB))  deallocate(DeltaSemiAll_VCB)
    if(allocated(DconsDsemiAll_VCB)) deallocate(DconsDsemiAll_VCB)
    if(allocated(SemiState_VGB))     deallocate(SemiState_VGB)
    if(allocated(ResSemi_VCB))       deallocate(ResSemi_VCB)
    if(allocated(JacSemi_VVCIB))     deallocate(JacSemi_VVCIB)
    if(allocated(FluxImpl_VXB))      deallocate(FluxImpl_VXB)
    if(allocated(FluxImpl_VYB))      deallocate(FluxImpl_VYB)
    if(allocated(FluxImpl_VZB))      deallocate(FluxImpl_VZB)
    if(allocated(Rhs_I))             deallocate(Rhs_I)
    if(allocated(x_I))               deallocate(x_I)
    if(allocated(JacobiPrec_I))      deallocate(JacobiPrec_I)
    if(allocated(iVectorSemi_I))     deallocate(iVectorSemi_I)

  end subroutine clean_mod_semi_impl
  !============================================================================
  subroutine advance_semi_impl

    ! Advance semi-implicit terms

    use ModProcMH, ONLY: iComm
    use ModMain, ONLY: NameThisComp, time_accurate, BlkTest, ProcTest
    use ModAdvance, ONLY: DoFixAxis
    use ModCoarseAxis, ONLY: UseCoarseAxis, coarsen_axis_cells
    use ModGeometry, ONLY: true_cell
    use ModImplHypre, ONLY: hypre_initialize, hypre_preconditioner
    use ModLinearSolver, ONLY: solve_linear_multiblock
    use ModMessagePass, ONLY: exchange_messages
    use BATL_lib, ONLY: nDim, nI, nJ, nK, nBlock, Unused_B

    use ModRadDiffusion,   ONLY: &
         get_impl_rad_diff_state, set_rad_diff_range, update_impl_rad_diff
    use ModHeatConduction, ONLY: &
         get_impl_heat_cond_state, update_impl_heat_cond
    use ModResistivity,    ONLY: &
         get_impl_resistivity_state, update_impl_resistivity
    use ModFieldLineThread, ONLY: UseFieldLineThreads, advance_threads, Heat_

    integer :: iBlockSemi, iBlock, iError1, i, j, k, iVar, n

    logical :: DoTest, DoTestMe, DoTestKrylov

    character(len=20) :: NameSub = 'MH_advance_semi_impl'
    !--------------------------------------------------------------------------
    NameSub(1:2) = NameThisComp
    call set_oktest('krylov', DoTest, DoTestKrylov)
    call set_oktest(NameSub, DoTest, DoTestMe) 
    if(DoTestMe) write(*,*)NameSub,' starting'

    if(TypeSemiImplicit /= 'resistivity')then
       ! All used blocks are solved for with the semi-implicit scheme
       ! except for (Hall) resistivity
       nBlockSemi = 0
       do iBlock = 1, nBlock
          if(Unused_B(iBlock)) CYCLE
          nBlockSemi = nBlockSemi + 1
          iBlockFromSemi_B(nBlockSemi) = iBlock

          ! Set the test block
          if(iProc == ProcTest .and. iBlock == BlkTest) &
               iBlockSemiTest = nBlockSemi
       end do
       DconsDsemiAll_VCB(:,:,:,:,1:nBlockSemi) = 1.0 
    else
       ! For (Hall) resistivity the number of semi-implicit blocks will be 
       ! set in get_impl_resistivity_state, so initialize up to nBlock
       DconsDsemiAll_VCB(:,:,:,:,1:nBlock) = 1.0
    end if

    ! Get current state and dCons/dSemi derivatives
    select case(TypeSemiImplicit)
    case('radiation', 'radcond', 'cond')
          call get_impl_rad_diff_state(SemiAll_VCB, DconsDsemiAll_VCB, &
               DeltaSemiAll_VCB=DeltaSemiAll_VCB, &
               DoCalcDeltaIn=UseStableImplicit)
    case('parcond')
       call get_impl_heat_cond_state(SemiAll_VCB, DconsDsemiAll_VCB, &
            DeltaSemiAll_VCB=DeltaSemiAll_VCB, &
            DoCalcDeltaIn=UseStableImplicit)
    case('resistivity')
       call get_impl_resistivity_state(SemiAll_VCB)
    case default
       call stop_mpi(NameSub//': no get_impl_state implemented for' &
            //TypeSemiImplicit)
    end select

    ! Re-initialize HYPRE preconditioner
    if(SemiParam%TypePrecond == 'HYPRE') call hypre_initialize

    ! For nVarSemi = 1,  loop through all semi-implicit variables one-by-one
    ! For nVarSemi = nVarSemiAll, do all (semi-)implicit variables together
    do iVarSemi = 1, nVarSemiAll, nVarSemi

       if(UseSplitSemiImplicit)then
          iVarSemiMin = iVarSemi
          iVarSemiMax = iVarSemi
          call set_rad_diff_range(iVarSemi)
       end if

       ! Set right hand side
       call get_semi_impl_rhs(SemiAll_VCB, Rhs_I)

       ! Calculate Jacobian matrix if required
       if(SemiParam%DoPrecond)then
          call timing_start('impl_jacobian')
          call get_semi_impl_jacobian
          call timing_stop('impl_jacobian')
       endif

       ! solve implicit system
       call solve_linear_multiblock( SemiParam, &
            nVarSemi, nDim, nI, nJ, nK, nBlockSemi, iComm, &
            semi_impl_matvec, Rhs_I, x_I, DoTestKrylov, &
            JacSemi_VVCIB, JacobiPrec_I, cg_precond, hypre_preconditioner)

       if(SemiParam%iError /= 0 .and. iProc == 0 .and. time_accurate) &
            call error_report(NameSub//': Krylov solver failed, Krylov error',&
            SemiParam%Error, iError1, .true.)

       ! NewSemiAll_VCB = SemiAll_VCB + Solution
       n=0
       do iBlockSemi = 1, nBlockSemi; do k=1,nK; do j=1,nJ; do i=1,nI
          do iVar = iVarSemiMin, iVarSemiMax
             n = n + 1
             if(true_cell(i,j,k,iBlockFromSemi_B(iBlockSemi)))then
                NewSemiAll_VCB(iVar,i,j,k,iBlockSemi) = &
                     SemiAll_VCB(iVar,i,j,k,iBlockSemi) + x_I(n)
             else
                NewSemiAll_VCB(iVar,i,j,k,iBlockSemi) = &
                     SemiAll_VCB(iVar,i,j,k,iBlockSemi)
             end if
          enddo
       enddo; enddo; enddo; enddo
    end do ! Splitting

    ! Put back semi-implicit result into the explicit code
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)
       select case(TypeSemiImplicit)
       case('radiation', 'radcond', 'cond')
          call update_impl_rad_diff(iBlock, iBlockSemi, &
               NewSemiAll_VCB(:,:,:,:,iBlockSemi), &
               SemiAll_VCB(:,:,:,:,iBlockSemi), &
               DconsDsemiAll_VCB(:,:,:,:,iBlockSemi))
       case('parcond')
          call update_impl_heat_cond(iBlock, iBlockSemi, &
               NewSemiAll_VCB(:,:,:,:,iBlockSemi), &
               SemiAll_VCB(:,:,:,:,iBlockSemi), &
               DconsDsemiAll_VCB(:,:,:,:,iBlockSemi))
       case('resistivity')
          call update_impl_resistivity(iBlock, &
               NewSemiAll_VCB(:,:,:,:,iBlockSemi))
       case default
          call stop_mpi(NameSub//': no update_impl implemented for' &
               //TypeSemiImplicit)
       end select
    end do



    if(DoFixAxis)call fix_axis_cells
    if(UseCoarseAxis)call coarsen_axis_cells
    ! Exchange messages, so ghost cells of all blocks are updated
    if(UseFieldLineThreads) &
         call advance_threads(Heat_)
    call exchange_messages

  end subroutine advance_semi_impl

  !============================================================================
  subroutine get_semi_impl_rhs_block(iBlock, SemiState_VG, RhsSemi_VC, &
       IsLinear)

    use BATL_lib,          ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK
    use ModLinearSolver,   ONLY: UsePDotADotP
    use ModRadDiffusion,   ONLY: get_rad_diffusion_rhs
    use ModHeatConduction, ONLY: get_heat_conduction_rhs
    use ModResistivity,    ONLY: get_resistivity_rhs

    integer, intent(in):: iBlock
    real, intent(inout):: SemiState_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK)
    real, intent(out)  :: RhsSemi_VC(nVarSemi,nI,nJ,nK)
    logical, intent(in):: IsLinear

    character(len=*), parameter:: NameSub = 'get_semi_impl_rhs_block'
    !------------------------------------------------------------------------
    select case(TypeSemiImplicit)
    case('radiation', 'radcond', 'cond')
       if(.not.IsLinear) UsePDotADotP = .false.
       call get_rad_diffusion_rhs(iBlock, SemiState_VG, &
            RhsSemi_VC, IsLinear=IsLinear)
    case('parcond')
       call get_heat_conduction_rhs(iBlock, SemiState_VG, &
            RhsSemi_VC, IsLinear=IsLinear)
    case('resistivity')
       call get_resistivity_rhs(iBlock, SemiState_VG, &
            RhsSemi_VC, IsLinear=IsLinear)
    case default
       call stop_mpi(NameSub//': no get_rhs implemented for' &
            //TypeSemiImplicit)
    end select


    
  end subroutine get_semi_impl_rhs_block
  !============================================================================
  subroutine get_semi_impl_rhs(SemiAll_VCB, RhsSemi_VCB)

    use ModAdvance,        ONLY: time_BLK
    use ModMain,           ONLY: time_accurate, dt, Cfl
    use ModGeometry,       ONLY: far_field_BCs_BLK
    use ModSize,           ONLY: nI, nJ, nK
    use ModCellBoundary,   ONLY: set_cell_boundary
    use BATL_lib,          ONLY: message_pass_cell, message_pass_face, &
         apply_flux_correction_block, CellVolume_GB

    real, intent(in) :: SemiAll_VCB(nVarSemiAll,nI,nJ,nK,nBlockSemi)
    real, intent(out):: RhsSemi_VCB(nVarSemi,nI,nJ,nK,nBlockSemi)

    integer :: iBlockSemi, iBlock, i, j, k

    integer :: iVarSemi_

    real :: DtLocal

    character(len=*), parameter:: NameSub = 'get_semi_impl_rhs'
    !------------------------------------------------------------------------
    ! Fill in SemiState so it can be message passed
    SemiState_VGB = 0.0
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          SemiState_VGB(:,i,j,k,iBlock) = &
               SemiAll_VCB(iVarSemiMin:iVarSemiMax,i,j,k,iBlockSemi)
       end do; end do; end do
    end do

    ! Message pass to fill in ghost cells 
    select case(TypeSemiImplicit)
    case('radiation', 'radcond', 'cond')
       if(UseAccurateRadiation)then
          call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=2, &
               nProlongOrderIn=1, nCoarseLayerIn=2, DoRestrictFaceIn = .true.)
       else
          call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=1, &
               nProlongOrderIn=1, DoSendCornerIn=.false., &
               DoRestrictFaceIn=.true.)
       end if
    case('parcond','resistivity')
       call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=2, &
            nProlongOrderIn=1, nCoarseLayerIn=2, DoRestrictFaceIn = .true.)
    case default
       call stop_mpi(NameSub//': no get_rhs message_pass implemented for' &
            //TypeSemiImplicit)
    end select

    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)

       ! Apply boundary conditions (1 layer of outer ghost cells)
       if(far_field_BCs_BLK(iBlock))&
            call set_cell_boundary(1, iBlock, nVarSemi, &
            SemiState_VGB(:,:,:,:,iBlock), iBlockSemi, IsLinear=.false.)

       call get_semi_impl_rhs_block(iBlock, SemiState_VGB(:,:,:,:,iBlock), &
            RhsSemi_VCB(:,:,:,:,iBlockSemi), IsLinear=.false.)
    end do

    if(TypeSemiImplicit == 'parcond' .or. TypeSemiImplicit == 'resistivity' &
         .or. UseAccurateRadiation)then
       call message_pass_face(nVarSemi, &
            FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)

       do iBlockSemi = 1, nBlockSemi
          iBlock = iBlockFromSemi_B(iBlockSemi)

          ! there are no ghost cells for RhsSemi_VCB
          call apply_flux_correction_block(iBlock, nVarSemi, 0, 0, &
               RhsSemi_VCB(:,:,:,:,iBlockSemi), &
               Flux_VXB=FluxImpl_VXB, Flux_VYB=FluxImpl_VYB, &
               Flux_VZB=FluxImpl_VZB)
       end do
    end if

    ! Multiply with cell volume (makes matrix symmetric)
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          RhsSemi_VCB(:,i,j,k,iBlockSemi) = RhsSemi_VCB(:,i,j,k,iBlockSemi) &
               *CellVolume_GB(i,j,k,iBlock)
       end do; end do; end do
    end do
    
    if(UseStableImplicit) then
       DtLocal = dt
       do iBlockSemi = 1, nBlockSemi
          iBlock = iBlockFromSemi_B(iBlockSemi)
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not. time_accurate) &
                  DtLocal = max(1.0e-30,Cfl*time_BLK(i,j,k,iBlock))
             RhsSemi_VCB(:,i,j,k,iBlockSemi) = &
                  RhsSemi_VCB(:,i,j,k,iBlockSemi) &
                  + DeltaSemiAll_VCB(:,i,j,k,iBlockSemi) &
                  *DconsDsemiAll_VCB(:,i,j,k,iBlockSemi) &
                  *CellVolume_GB(i,j,k,iBlock)/DtLocal                  
          end do; end do; end do
       enddo
    endif
    
  end subroutine get_semi_impl_rhs
  !============================================================================
  subroutine semi_impl_matvec(x_I, y_I, MaxN)

    ! Calculate y_I = A.x_I where A is the linearized sem-implicit operator

    use ModAdvance,  ONLY: time_BLK
    use ModGeometry, ONLY: far_field_BCs_BLK
    use ModMain, ONLY: dt, time_accurate, Cfl
    use ModSize, ONLY: nI, nJ, nK
    use ModLinearSolver,   ONLY: UsePDotADotP, pDotADotPPe
    use ModCellBoundary,   ONLY: set_cell_boundary
    use BATL_lib, ONLY: message_pass_cell, message_pass_face, &
         apply_flux_correction_block, CellVolume_GB

    integer, intent(in) :: MaxN
    real,    intent(in) :: x_I(MaxN)
    real,    intent(out):: y_I(MaxN)

    integer :: iBlockSemi, iBlock, i, j, k, iVar, n
    real :: Volume, DtLocal

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'semi_impl_matvec'
    !------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)
    call timing_start(NameSub)

    ! Fill in StateSemi so it can be message passed
    n = 0
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)
       do k = 1, nK; do j = 1, nJ; do i = 1, nI; do iVar = 1, nVarSemi
          n = n + 1
          SemiState_VGB(iVar,i,j,k,iBlock) = x_I(n)
       end do; end do; end do; end do
    end do

    ! Message pass to fill in ghost cells 
    select case(TypeSemiImplicit)
    case('radiation', 'radcond', 'cond')

       if(UseAccurateRadiation)then
          UsePDotADotP = .false.

          call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=2, &
               nProlongOrderIn=1, nCoarseLayerIn=2, DoRestrictFaceIn = .true.)
       else
          ! Initialize the computation of (p . A . P) form
          UsePDotADotP = SemiParam%TypeKrylov == 'CG'

          pDotADotPPe = 0.0

          call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=1, &
               nProlongOrderIn=1, DoSendCornerIn=.false., &
               DoRestrictFaceIn=.true.)
       end if
    case('parcond','resistivity')
       call message_pass_cell(nVarSemi, SemiState_VGB, nWidthIn=2, &
            nProlongOrderIn=1, nCoarseLayerIn=2, DoRestrictFaceIn = .true.)
    case default
       call stop_mpi(NameSub//': no get_rhs message_pass implemented for' &
            //TypeSemiImplicit)
    end select

    n = 0
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)

       if(far_field_BCs_BLK(iBlock))call set_cell_boundary( &
            1, iBlock, nVarSemi, &
            SemiState_VGB(:,:,:,:,iBlock), iBlockSemi, IsLinear=.true.)

       call get_semi_impl_rhs_block(iBlock, SemiState_VGB(:,:,:,:,iBlock), &
            ResSemi_VCB(:,:,:,:,iBlockSemi), IsLinear=.true.)

       if(UsePDotADotP)then
          DtLocal = Dt
          if(UseSplitSemiImplicit)then
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.time_accurate) &
                     DtLocal = max(1.0e-30,Cfl*time_BLK(i,j,k,iBlock))
                Volume = CellVolume_GB(i,j,k,iBlock)/DtLocal
                n = n + 1
                pDotADotPPe = pDotADotPPe +  &
                     Volume*x_I(n)**2 &
                     *DconsDsemiAll_VCB(iVarSemi,i,j,k,iBlockSemi) &
                     /(DtLocal * SemiImplCoeff)
             end do; enddo; enddo
          else
             do k = 1, nK; do j = 1, nJ; do i = 1, nI
                if(.not.time_accurate) &
                     DtLocal = max(1.0e-30,Cfl*time_BLK(i,j,k,iBlock))
                Volume = CellVolume_GB(i,j,k,iBlock) !!! /DtLocal ???
                do iVar = 1, nVarSemi
                   n = n + 1
                   pDotADotPPe = pDotADotPPe +  &
                        Volume*x_I(n)**2 &
                        *DconsDsemiAll_VCB(iVar,i,j,k,iBlockSemi) &
                        /(DtLocal * SemiImplCoeff)
                enddo
             enddo; enddo; enddo
          end if
       end if

    end do

    if(TypeSemiImplicit == 'parcond' .or. TypeSemiImplicit == 'resistivity' &
         .or. UseAccurateRadiation)then
       call message_pass_face(nVarSemi, &
            FluxImpl_VXB, FluxImpl_VYB, FluxImpl_VZB)

       do iBlockSemi = 1, nBlockSemi
          iBlock = iBlockFromSemi_B(iBlockSemi)

          ! zero ghost cells for ResSemi_VCB
          call apply_flux_correction_block(iBlock, nVarSemi,0, 0, &
               ResSemi_VCB(:,:,:,:,iBlockSemi), &
               Flux_VXB=FluxImpl_VXB, Flux_VYB=FluxImpl_VYB, &
               Flux_VZB=FluxImpl_VZB)
       end do
    end if

    n=0
    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)

       DtLocal = dt
       if(UseSplitSemiImplicit)then
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.time_accurate) &
                  DtLocal = max(1.0e-30,Cfl*time_BLK(i,j,k,iBlock))
             Volume = CellVolume_GB(i,j,k,iBlock)
             n = n + 1
             y_I(n) = Volume* &
                  (x_I(n)*DconsDsemiAll_VCB(iVarSemi,i,j,k,iBlockSemi)/DtLocal&
                  - SemiImplCoeff * ResSemi_VCB(1,i,j,k,iBlockSemi))
          end do; enddo; enddo
       else
          do k = 1, nK; do j = 1, nJ; do i = 1, nI
             if(.not.time_accurate) &
                  DtLocal = max(1.0e-30,Cfl*time_BLK(i,j,k,iBlock))
             Volume = CellVolume_GB(i,j,k,iBlock)
             do iVar = 1, nVarSemi
                n = n + 1
                y_I(n) = Volume* &
                     (x_I(n)*DconsDsemiAll_VCB(iVar,i,j,k,iBlockSemi)/DtLocal &
                     - SemiImplCoeff * ResSemi_VCB(iVar,i,j,k,iBlockSemi))
             enddo
          enddo; enddo; enddo
       end if
    end do


    if (UsePDotADotP)then
       pDotADotPPe = pDotADotPPe * SemiImplCoeff
    else
       pDotADotPPe = 0.0
    end if

    ! Apply left preconditioning if required: y --> P_L.y
    if(SemiParam%DoPrecond .and. SemiParam%TypeKrylov /= 'CG') &
         call semi_precond(MaxN, y_I)

    call timing_stop(NameSub)

  end subroutine semi_impl_matvec
  !============================================================================
  subroutine cg_precond(Vec_I, PrecVec_I, MaxN)

    ! Set PrecVec = Prec.Vec where Prec is the preconditioner matrix. 
    ! This routine is used by the Preconditioned Conjugate Gradient method

    integer, intent(in) :: MaxN
    real,    intent(in) :: Vec_I(MaxN)
    real,    intent(out):: PrecVec_I(MaxN)
    !-------------------------------------------------------------------------
    PrecVec_I = Vec_I
    call semi_precond(MaxN, PrecVec_I)

  end subroutine cg_precond
  !============================================================================
  subroutine semi_precond(MaxN, Vec_I)

    ! Multiply Vec with the preconditioner matrix. 
    ! This routine is used by the Preconditioned Conjugate Gradient method

    use BATL_size,       ONLY: nDim, nI, nJ, nK
    use ModLinearSolver, ONLY: precond_left_multiblock
    use ModImplHypre,    ONLY: hypre_preconditioner

    integer, intent(in)   :: MaxN
    real,    intent(inout):: Vec_I(MaxN)
    !-------------------------------------------------------------------------
    select case(SemiParam%TypePrecond)
    case('HYPRE')
       call hypre_preconditioner(MaxN, Vec_I)
    case('JACOBI')
       Vec_I = JacobiPrec_I(1:MaxN)*Vec_I
    case default
       call precond_left_multiblock(SemiParam,           &
            nVarSemi, nDim, nI, nJ, nK, nBlockSemi, JacSemi_VVCIB, Vec_I)
    end select

  end subroutine semi_precond
  !============================================================================

  subroutine test_semi_impl_jacobian

    ! Calculate the Jacobian Jac_VVI = d(RHS)/d(Var) for the test cell 
    ! using numerical derivatives of the RHS obtained with 
    ! get_semi_impl_rhs_block

    use BATL_lib, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK, nI, nJ, nK, &
         IsRotatedCartesian, rot_to_cart

    use ModMain, ONLY: iTest, jTest, kTest, BlkTest

    ! Local variables
    real, parameter:: RelEps = 1e-6, AbsEps = 1e-8

    real, allocatable:: SemiState_VG(:,:,:,:), SemiPert_VG(:,:,:,:), &
         Rhs_VC(:,:,:,:), RhsPert_VC(:,:,:,:), &
         JacAna_VVCI(:,:,:,:,:,:), JacNum_VVI(:,:,:), &
         JacAna_VV(:,:), JacNum_VV(:,:)

    real:: JacAna, JacNum

    integer:: i, j, k, iPert, jPert, kPert, iBlock, iStencil, iVar, jVar

    character(len=*), parameter:: NameSub = 'test_semi_impl_jacobian'
    !--------------------------------------------------------------------------

    allocate( &
         SemiState_VG(nVarSemi, MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         SemiPert_VG(nVarSemi,MinI:MaxI,MinJ:MaxJ,MinK:MaxK), &
         Rhs_VC(nVarSemi,nI,nJ,nK), RhsPert_VC(nVarSemi,nI,nJ,nK), &
         JacAna_VVCI(nVarSemi,nVarSemi,nI,nJ,nK,nStencil), &
         JacNum_VVI(nVarSemi,nVarSemi,nStencil), &
         JacAna_VV(nVarSemi,nVarSemi), JacNum_VV(nVarSemi,nVarSemi))

    ! For sake of simpler code
    iBlock = BlkTest

    ! Get analytic Jacobian
    call get_semi_impl_jacobian_block(iBlock, JacAna_VVCI)

    ! Set semi-implicit state with nG ghost cells. 
    ! The ghost cells are set but not used.
    SemiState_VG = 1.0
    SemiState_VG(:,1:nI,1:nJ,1:nK) = &
         SemiAll_VCB(iVarSemiMin:iVarSemiMax,:,:,:,iBlockSemiTest)

    ! Initialize perturbed state
    SemiPert_VG = SemiState_VG

    ! Get original RHS (use SemiPert to keep compiler happy)
    call get_semi_impl_rhs_block(iBlock, SemiPert_VG, Rhs_VC, IsLinear=.false.)

    ! Select test cell
    i = iTest; j = jTest; k = kTest

    ! Initialize numerical Jacobian
    JacNum_VVI = 0.0

    ! Loop over stencil
    do iStencil = 1, nStencil

       ! Set indexes for perturbed cell
       iPert = i; jPert = j; kPert = k
       select case(iStencil)
       case(2)
          iPert = i - 1
          if(iPert == 0) CYCLE
       case(3)
          iPert = i + 1
          if(iPert > nI) CYCLE
       case(4)
          jPert = j - 1
          if(jPert == 0) CYCLE
       case(5)
          jPert = j + 1
          if(jPert > nJ) CYCLE
       case(6)
          kPert = k - 1
          if(kPert == 0) CYCLE
       case(7)
          kPert = k + 1
          if(kPert > nK) CYCLE
       end select

       ! Loop over variables to be perturbed
       do iVar = 1, nVarSemi

          ! Perturb the variable
          SemiPert_VG(iVar,iPert,jPert,kPert) &
               = SemiState_VG(iVar,iPert,jPert,kPert)*(1 + RelEps) + AbsEps

          ! Calculate perturbed RHS
          call get_semi_impl_rhs_block(iBlock, SemiPert_VG, RhsPert_VC, &
               IsLinear=.false.)

          ! Calculate Jacobian elements
          JacNum_VVI(:,iVar,iStencil) = &
               (RhsPert_VC(:,i,j,k) - Rhs_VC(:,i,j,k)) / &
               (  SemiPert_VG(iVar,iPert,jPert,kPert)    &
               - SemiState_VG(iVar,iPert,jPert,kPert) )

          ! Reset the variable to the original value
          SemiPert_VG(iVar,iPert,jPert,kPert) = &
               SemiState_VG(iVar,iPert,jPert,kPert)
       end do
    end do

    do iStencil = 1, nStencil

       JacAna_VV = JacAna_VVCI(:,:,i,j,k,iStencil)
       JacNum_VV = JacNum_VVI(:,:,iStencil)

       if(IsRotatedCartesian .and. TypeSemiImplicit == 'resistivity')then
          JacAna_VV = rot_to_cart(JacAna_VV)
          JacNum_VV = rot_to_cart(JacNum_VV)
       end if

       do jVar = 1, nVarSemi; do iVar = 1, nVarSemi
          JacAna = JacAna_VV(iVar,jVar)
          JacNum = JacNum_VV(iVar,jVar)
          write(*,'(a,3i2,a,4es13.5)') &
               NameSub//': iVar,jVar,iStencil=', iVar, jVar, iStencil, &
               ' JacAna, JacNum, AbsDiff, RelDiff=', &
               JacAna, JacNum, abs(JacAna - JacNum), &
               2*abs(JacAna - JacNum)/(abs(JacAna) + abs(JacNum) + 1e-30)
       end do; end do
    end do

    deallocate(SemiState_VG, SemiPert_VG, Rhs_VC, RhsPert_VC, &
         JacAna_VVCI, JacNum_VVI, JacAna_VV, JacNum_VV)

  end subroutine test_semi_impl_jacobian

  !============================================================================
  subroutine get_semi_impl_jacobian_block(iBlock, JacSemi_VVCI)

    use BATL_lib,          ONLY: nI, nJ, nK
    use ModRadDiffusion,   ONLY: add_jacobian_rad_diff
    use ModHeatConduction, ONLY: add_jacobian_heat_cond
    use ModResistivity,    ONLY: add_jacobian_resistivity

    integer, intent(in) :: iBlock
    real,    intent(out):: JacSemi_VVCI(nVarSemi,nVarSemi,nI,nJ,nK,nStencil)

    character(len=*), parameter:: NameSub = 'get_semi_impl_jacobian_block'
    !------------------------------------------------------------------------
    ! All elements have to be set
    JacSemi_VVCI = 0.0

    select case(TypeSemiImplicit)
    case('radiation', 'radcond', 'cond')
       call add_jacobian_rad_diff(iBlock, nVarSemi, JacSemi_VVCI)
    case('parcond')
       call add_jacobian_heat_cond(iBlock, nVarSemi, JacSemi_VVCI)
    case('resistivity')
       call add_jacobian_resistivity(iBlock, nVarSemi, JacSemi_VVCI)
    case default
       call stop_mpi(NameSub//': no add_jacobian implemented for' &
            //TypeSemiImplicit)
    end select

  end subroutine get_semi_impl_jacobian_block

  !============================================================================
  subroutine get_semi_impl_jacobian

    use ModAdvance, ONLY: time_BLK
    use ModMain,    ONLY: nI, nJ, nK, Dt, time_accurate, Cfl
    use ModGeometry, ONLY: true_cell
    use ModImplicit, ONLY: UseNoOverlap
    use ModImplHypre, ONLY: hypre_set_matrix_block, hypre_set_matrix
    use BATL_lib, ONLY: CellVolume_GB

    integer :: iBlockSemi, iBlock, i, j, k, iStencil, iVar
    real    :: Coeff, DtLocal

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'get_semi_impl_jacobian'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe) call test_semi_impl_jacobian

    ! The HYPRE AMG preconditioner requires the overlap between blocks
    ! For all other preconditioners it is better to avoid the overpap
    if(SemiParam%TypePrecond=='HYPRE') UseNoOverlap = .false.

    do iBlockSemi = 1, nBlockSemi
       iBlock = iBlockFromSemi_B(iBlockSemi)

       ! Get dR/dU
       call get_semi_impl_jacobian_block(iBlock, &
            JacSemi_VVCIB(:,:,:,:,:,:,iBlockSemi))

       ! Form A = Volume*(1/dt - SemiImplCoeff*dR/dU) 
       !    symmetrized for sake of CG       
       do iStencil = 1, nStencil; do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          Coeff = -SemiImplCoeff*CellVolume_GB(i,j,k,iBlock)
          JacSemi_VVCIB(:, :, i, j, k, iStencil, iBlockSemi) = &
               Coeff * JacSemi_VVCIB(:, :, i, j, k, iStencil, iBlockSemi)
       end do; end do; end do; end do
       DtLocal = dt
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          if(.not.time_accurate) &
               DtLocal = max(1.0e-30, Cfl*time_BLK(i,j,k,iBlock))
          Coeff = CellVolume_GB(i,j,k,iBlock)/DtLocal
          if(UseSplitSemiImplicit)then
             JacSemi_VVCIB(1,1,i,j,k,1,iBlockSemi) = &
                  Coeff*DconsDsemiAll_VCB(iVarSemi,i,j,k,iBlockSemi) &
                  + JacSemi_VVCIB(1,1,i,j,k,1,iBlockSemi)
          else
             do iVar = 1, nVarSemi
                JacSemi_VVCIB(iVar,iVar,i,j,k,1,iBlockSemi) = &             
                     Coeff*DconsDsemiAll_VCB(iVar,i,j,k,iBlockSemi) &
                     + JacSemi_VVCIB(iVar,iVar,i,j,k,1,iBlockSemi) 
             end do
          end if          
       end do; end do; end do

       if(SemiParam%TypePrecond == 'HYPRE') &
            call hypre_set_matrix_block(iBlockSemi, &
            JacSemi_VVCIB(1,1,1,1,1,1,iBlockSemi))
    end do

    if(SemiParam%TypePrecond == 'HYPRE') call hypre_set_matrix(.true.)

    UseNoOverlap = .true.

  end subroutine get_semi_impl_jacobian

end module ModSemiImplicit

