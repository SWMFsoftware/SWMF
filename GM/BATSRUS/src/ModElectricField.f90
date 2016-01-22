!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModElectricField

  ! Calculate electric field Efield_DGB from MHD quantities.
  !
  ! Calculate potential and inductive part of the electric field: Epot and Eind 
  ! Efield = Epot + Eind
  !
  ! where curl(Epot) = 0 and div(Eind)=0.
  !
  ! This is similar to the projection scheme. 
  ! NOTE: For convenience the sign of the potential is reversed.
  !
  ! Epot = grad(Potential), and Potential satisfies the Poisson equation
  !
  ! Laplace(Potential) = div(E)

  use BATL_lib, ONLY: nDim, nI, nJ, nK, nIJK, nG, &
       MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       nBlock, MaxBlock, Unused_B, &
       iProc, iComm, message_pass_cell

  use ModGeometry,     ONLY: far_field_bcs_blk, true_cell
  use ModCellBoundary, ONLY: set_cell_boundary
  use ModCellGradient, ONLY: calc_gradient, calc_divergence
  use ModLinearSolver, ONLY: LinearSolverParamType, solve_linear_multiblock

  implicit none
  SAVE

  private ! except

  public:: calc_inductive_e

  ! Total, potential and inductive electric fields
  real, public, allocatable:: Efield_DGB(:,:,:,:,:)
  real, public, allocatable:: Epot_DGB(:,:,:,:,:)
  real, public, allocatable:: Eind_DGB(:,:,:,:,:)

  ! Electric potential
  real, public, allocatable:: Potential_GB(:,:,:,:)

  ! local variables
  real, allocatable:: DivE_CB(:,:,:,:)
  real, allocatable:: Laplace_C(:,:,:)

  ! Default parameters for the linear solver
  type(LinearSolverParamType):: SolverParam = LinearSolverParamType( &
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

  ! number of blocks used
  integer:: nBlockUsed
  
  ! number of unknowns to solve for
  integer:: nVarAll

  ! indicate if we are in the linear stage or not
  logical:: IsLinear = .false.

contains
  !=========================================================================
  subroutine get_electric_field

    ! Fill in physical cells of Efield_DGB with electric field
    integer:: iBlock
    !----------------------------------------------------------------------
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call get_electric_field_block(iBlock)
    end do

  end subroutine get_electric_field

  !=========================================================================
  subroutine get_electric_field_block(iBlock)

    ! Fill in physical cells of Efield_DGB with electric field for block iBlock

    use ModVarIndexes, ONLY: Rho_, RhoUx_, RhoUz_, Bx_, Bz_
    use ModAdvance, ONLY: State_VGB
    use ModB0, ONLY: UseB0, B0_DGB
    use ModCoordTransform, ONLY: cross_product

    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: b_D(3), u_D(3)
    !----------------------------------------------------------------------
    if(.not.allocated(Efield_DGB)) &
         allocate(Efield_DGB(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    do k = 1, nK; do j = 1, nJ; do i = 1, nI
       if(.not. true_cell(i,j,k,iBlock))then
          Efield_DGB(:,i,j,k,iBlock) = 0.0
          CYCLE
       end if

       ! Ideal MHD case
       b_D = State_VGB(Bx_:Bz_,i,j,k,iBlock)
       if(UseB0) b_D = b_D + B0_DGB(:,i,j,k,iBlock)

       ! Single ion fluid case: E = - u x B
       u_D = State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) &
            /State_VGB(Rho_,i,j,k,iBlock)

       Efield_DGB(:,i,j,k,iBlock) = -cross_product(u_D, b_D)
    end do; end do; end do

  end subroutine get_electric_field_block

  !=========================================================================
  subroutine calc_inductive_e

    ! Calculate the inductive part of the electric field Eind

    integer:: iBlock
    !-----------------------------------------------------------------------

    nBlockUsed = nBlock - count(Unused_B(1:nBlock))
    nVarAll    = nBlockUsed*nIJK

    call get_electric_field

    call calc_div_e

    ! Calculate Potential_GB from the Poisson equation 
    ! fill ghost cells at the end
    call calc_potential

    ! Calculate Epot and Eind
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       ! Epot = grad(Potential)
       call calc_gradient(iBlock, Potential_GB(:,:,:,iBlock), nG, &
            Epot_DGB(:,:,:,:,iBlock), UseBodyCell=.true.)

       ! Eind = E - Epot
       Eind_DGB(:,:,:,:,iBlock) = Efield_DGB(:,:,:,:,iBlock) - &
            Epot_DGB(:,:,:,:,iBlock)
    end do

  end subroutine calc_inductive_e
  !=========================================================================
  subroutine calc_div_e

    integer:: iBlock
    !----------------------------------------------------------------------
    if(.not.allocated(DivE_CB)) allocate(DivE_CB(nI,nJ,nK,nBlock))

    ! Fill in ghost cells
    call message_pass_cell(3, Efield_DGB)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Fill outer ghost cells with floating values
       if(far_field_bcs_blk(iBlock)) &
            call set_cell_boundary(nG, iBlock, 3, Efield_DGB(:,:,:,:,iBlock),&
            TypeBcIn = 'float')

       ! Calculate DivE for this block (it has no ghost cells: nG=0)
       call calc_divergence(iBlock, Efield_DGB(:,:,:,:,iBlock), &
            0, DivE_CB(:,:,:,iBlock), UseBodyCell=.true.)
    end do

  end subroutine calc_div_e
  !=========================================================================
  subroutine calc_potential

    real, allocatable:: Rhs_I(:), Potential_I(:)

    integer:: iBlock, i, j, k, n

    logical:: DoTest, DoTestMe

    character(len=*), parameter:: NameSub = 'calc_potential'
    !---------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(.not.allocated(Potential_GB)) &
         allocate(Potential_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    ! Make potential zero in unused cells
    Potential_GB(:,:,:,1:nBlock) = 0.0

    ! Make Epot = Efield in outer boundary ghost cells
    Epot_DGB(:,:,:,:,1:nBlock) = Efield_DGB(:,:,:,:,1:nBlock)

    allocate(Rhs_I(nVarAll), Potential_I(nVarAll))

    ! Substitute 0 into the Laplace operator WITH boundary conditions
    ! to get the RHS
    Potential_I = 0.0
    IsLinear = .false.
    call matvec(Potential_I, Rhs_I, nVarAll)
    Rhs_I = -Rhs_I

    ! Add div(E) to the RHS
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          Rhs_I(n) = Rhs_I(n) + DivE_CB(i,j,k,iBlock)
       end do; end do; end do
    end do

    ! Start the linear solve stage
    IsLinear = .true.

    ! Make potential and electric field 0 in unused and boundary cells
    Potential_GB(:,:,:,1:nBlock) = 0.0
    Epot_DGB(:,:,:,:,1:nBlock)   = 0.0

    ! solve Poisson equation Laplace(Potential) = div(E)
    call solve_linear_multiblock( SolverParam, &
         1, nDim, nI, nJ, nK, nBlockUsed, iComm, &
         matvec_inductive_e, Rhs_I, Potential_I, DoTestMe)

    if(SolverParam%iError /= 0 .and. iProc == 0) &
         write(*,*) NameSub,' failed in ModInductiveE'

    ! Put solution x_I into Potential_GB
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          if(true_cell(i,j,k,iBlock)) &
               Potential_GB(i,j,k,iBlock) = Potential_I(n)
       end do; end do; end do
    end do

    ! Fill in ghost cells and boundary cells using true BCs
    IsLinear = .false.
    call bound_phi

    deallocate(Rhs_I, Potential_I)

  end subroutine calc_potential
  !=========================================================================
  subroutine bound_potential

    use ModGeometry, ONLY: body_blk, r_BLK
    use ModPhysics,  ONLY: rBody
    use BATL_lib, ONLY: Xyz_DGB

    ! Fill in ghost cells and boundary cells for the potential

    integer:: iBlock, i, j, k
    real:: rInside
    character(len=10):: TypeBc
    !----------------------------------------------------------------------
    ! Fill in ghost cells for Potential_GB
    call message_pass_cell(Potential_GB)

    if(IsLinear)then
       ! In the linear stage use float BC that corresponds to E_normal=0
       TypeBc = 'float'
    else
       ! Apply gradient = Enormal BC
       TypeBc = 'gradpot'
    end if

    ! Fill outer ghost cells 
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(far_field_bcs_blk(iBlock)) &
            call set_cell_boundary(nG, iBlock, 1, Potential_GB(:,:,:,iBlock),&
            TypeBcIn = TypeBc)
    end do

    if(IsLinear) RETURN

    ! Fill in inner boundary cells with ionosphere potential
    rInside = 1.01
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(.not.body_BLK(iBLOCK)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
          if(r_BLK(i,j,k,iBlock) > rBody) CYCLE
          if(r_BLK(i,j,k,iBlock) < rInside) CYCLE
          call get_ie_potential(Xyz_DGB(:,i,j,k,iBlock), &
               Potential_GB(i,j,k,iBlock))
       end do; end do; end do
    end do

  end subroutine bound_potential

  !=========================================================================
  subroutine matvec_inductive_e(x_I, y_I, MaxN)

    ! Calculate y = Laplace(x)

    integer, intent(in):: MaxN
    real,    intent(in) :: x_I(MaxN)
    real,    intent(out):: y_I(MaxN)

    integer:: iBlock, i, j, k, n
    !---------------------------------------------------------------------
    if(.not.allocated(Laplace_C)) allocate(Laplace_C(nI,nJ,nK))

    ! x_I -> Potential_GB
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          if(true_cell(i,j,k,iBlock)) Potential_GB(i,j,k,iBlock) = x_I(n)
       end do; end do; end do
    end do

    ! Boundary cells are filled in with zero
    call bound_potential

    ! Calculate gradient of potential
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       call calc_gradient(iBlock, Potential_GB(:,:,:,iBlock), nG, &
            Epot_DGB(:,:,:,:,iBlock), UseBodyCell=.true.)
    end do

    ! Fill in ghost cells
    call message_pass_cell(3,Epot_DGB)

    ! Calculate Laplace_C and copy it into y_I
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Fill outer ghost cells with floating values
       if(.not.IsLinear .and. far_field_bcs_blk(iBlock)) &
            call set_cell_boundary(nG, iBlock, 3, Epot_DGB(:,:,:,:,iBlock),&
            TypeBcIn = 'efield')

       call calc_divergence(iBlock, Epot_DGB(:,:,:,:,iBlock), 0, Laplace_C, &
            UseBodyCell=.true.)

       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          y_I(n) = Laplace_C(i,j,k)
       end do; end do; end do
    end do

  end subroutine matvec_inductive_e

end module ModElectricField

