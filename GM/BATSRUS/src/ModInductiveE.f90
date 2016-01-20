!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModInductiveE

  ! Calculate inductive electric field Eind 
  ! from the total electric field Efield
  ! minus the potential electric field Epot.
  ! Epot = -grad Potential, and Potential satisfies the Poisson equation
  !
  ! Laplace Potential = - div(E)

  implicit none
  SAVE

  private ! except

  public:: calc_inductive_e

  ! Total and inductive electric fields
  real, public, allocatable:: Efield_DGB(:,:,:,:,:)
  real, public, allocatable:: Eind_DGB(:,:,:,:,:)

  ! Electric potential
  real, public, allocatable:: Potential_GB(:,:,:,:)

  ! local variables
  real, allocatable:: DivE_CB(:,:,:,:)

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

contains
  !=========================================================================
  subroutine calc_inductive_e

    if(.not.allocated(Efield_DGB)) &
         allocate(Efield_DGB(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))
    if(.not.allocated(DivE_CB)) &
         allocate(DivE_CB(nI,nJ,nK,MaxBlock))
    if(.not.allocated(Potential_GB)) &
         allocate(Potential_GB(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock))

    nBlockUsed = nBlock - count(Unused_B(1:nBlock))
    nVarAll    = nBlockUsed*nIJK

    call get_electric_field(Efield_DGB)

    call calc_div_e

    ! Calculate Potential_GB from the Poisson equation 
    ! fill ghost cells at the end
    call calc_potential

    ! Calculate the inductive E field
    ! fill ghost cells at the end
    call calc_inductive_field

  end subroutine calc_inductive_e
  !=========================================================================
  subroutine calc_div_e

    integer:: iBlock
    !----------------------------------------------------------------------
    ! Fill in ghost cells
    call message_pass_cell(3, Efield_DGB)

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       ! Fill outer ghost cells with floating values
       if(far_field_bcs_blk(iBlock)) &
            call set_cell_boundary(nG, iBlock, 3, Efield_DGB(:,:,:,:,iBlock),&
            TypeBcIn = 'float')

       ! Calculate DivE for this block
       call calc_divergence(iBlock, Efield_DGB(:,:,:,:,iBlock), &
            nG=0, DivE_CB(:,:,:,iBlock))
    end do

  end subroutine calc_div_e
  !=========================================================================
  subroutine calc_potential

    character(len=*), parameter:: calc_potential

    allocate(Rhs_I(nVarAll), Potential_I(nVarAll))

    ! Substitute 0 into the Laplace operator WITH boundary conditions
    ! to get the RHS
    Potential_I = 0.0
    UseIonoPotential = .true.
    call matvec(Potential_I, Rhs_I, nVarAll)
    Rhs_I = -Rhs_I

    ! Add -div(E) to the RHS
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          if(.not.true_cell(i,j,k,iBlock)) CYCLE
          Rhs_I(n) = Rhs_I(n) - DivE_CB(i,j,k,iBlock)
       end do; end do; end do
    end do
    UseIonoPotential = .false.

    ! solve implicit system
    call solve_linear_multiblock( SolverParam, &
         1, nDim, nI, nJ, nK, nBlockUsed, iComm, &
         inductive_e_matvec, Rhs_I, x_I, DoTestKrylov)

    if(SemiParam%iError /= 0 .and. iProc == 0) &
         write(*,*) NameSub,' failed in ModInductiveE'

    ! Put solution into Potential_GB
    n = 0
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       do k = 1, nK; do j = 1, nJ; do i = 1, nI
          n = n + 1
          if(true_cell(i,j,k,iBlock))then
             Potential_GB(i,j,k,iBlock) = x_I(n)
          else
             Potential_GB(i,j,k,iBlock) = 0.0
          end if
       end do; end do; end do
    end do

  end subroutine calc_potential
  !=========================================================================
  subroutine bound_potential


    ! Fill outer ghost cells with floating values
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE

       if(far_field_bcs_blk(iBlock)) &
            call set_cell_boundary(nG, iBlock, 1, Potential_GB(:,:,:,iBlock),&
            TypeBcIn = 'float')
    end do

    if(.not.UseIonoPotential) CYCLE

    ! Fill in inner boundary cells with ionosphere potential
    rInside = rBody - 0.4
    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(body_BLK(iBLOCK))then
          do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
             if(r_BLK(i,j,k,iBlock) > rBody) CYCLE
             if(r_BLK(i,j,k,iBlock) < rInside) CYCLE
             call get_ie_potential(Xyz_DGB(:,i,j,k,iBlock), &
                  Potential_GB(i,j,k,iBlock))
          end do; end do; end do
       end do
    end do

  end subroutine bound_potential

  !=========================================================================
  subroutine matvec_inductive_e(x_I, y_I, MaxN)

    ! x_I -> Potential_GB


    call message_pass_cell(Potential_GB)
    call bound_potential
    call get_gradient
    call message_pass_cell(Potential_GB)
    call get_divergence

    ! Laplace_CB -> y_I

  end subroutine matvec_inductive_e

end module ModInductiveE

