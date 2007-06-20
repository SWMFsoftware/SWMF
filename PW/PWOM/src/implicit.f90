subroutine PW_implicit(ImplPar, DtExpl, DtImpl, nCell, nVar, State_GV, &
     calc_residual, update_boundary)

  use ModLinearSolver, ONLY: prehepta, Lhepta, Uhepta

  implicit none

  real,    intent(in) :: ImplPar, DtExpl, DtImpl
  integer, intent(in) :: nCell, nVar
  real, intent(inout) :: State_GV(-1:nCell+2, nVar)

  interface
     subroutine calc_residual(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)
       implicit none
       integer, intent(in) :: nOrder, nCell, nVar
       real,    intent(in) :: Dt
       real,    intent(in) :: State_GV(-1:nCell+2, nVar)
       real,    intent(out):: Resid_CV(nCell, nVar)
     end subroutine calc_residual

     subroutine update_boundary(nCell, nVar, State_GV)
       implicit none
       integer, intent(in)    :: nCell, nVar
       real,    intent(inout) :: State_GV(-1:nCell+2, nVar)
     end subroutine update_boundary
  end interface

  real, parameter :: Eps = 1.e-6
  real, allocatable, dimension(:,:) :: &
       RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV

  integer, parameter :: nDiag = 3
  real, allocatable :: Matrix_VVCI(:, :, :, :)
       
  real, allocatable :: Norm_V(:)   ! second norm of variables
  real, allocatable :: x_I(:)      ! linear vector of right hand side/unknowns
  real :: Coeff

  integer :: iVar, jVar, i, iStencil, iDiag, iX

  !-----------------------------------------------------------------------

  allocate(StateEps_GV(-1:nCell+2, nVar), &
       RightHand_CV(nCell, nVar), ResidOrig_CV(nCell, nVar),  &
       ResidEps_CV(nCell, nVar), Norm_V(nVar), x_I(nCell*nVar), &
       Matrix_VVCI(nVar, nVar, nCell, nDiag))

  ! Make sure that ghost cells are up-to-date
  call update_boundary(nCell, nVar, State_GV)

  ! calculate right hand side
  call calc_residual(2, DtImpl, nCell, nVar, State_GV, RightHand_CV)

  ! Calculate the unperturbed residual with first order scheme
  call calc_residual(1, DtExpl, nCell, nVar, State_GV, ResidOrig_CV)

  ! Calculate the norm for the variables
  do iVar = 1, nVar
     Norm_V(iVar) = sqrt(sum(State_GV(1:nCell,iVar)**2)/nCell)
  end do

  ! Calculate the dR/dU matrix
  do jVar = 1, nVar
     do iStencil = 1, 3
        ! Get perturbed state
        StateEps_GV(1:nCell,:) = State_GV
        do i = iStencil, nCell, 3
           StateEps_GV(i, jVar) = StateEps_GV(i, jVar) + Eps*Norm_V(jVar)
        end do
        call update_boundary(nCell, nVar, StateEps_GV)
        call calc_residual(1, DtExpl, nCell, nVar, StateEps_GV, ResidEps_CV)

        ! Jacobian is multiplied with -ImplPar*DtImpl
        Coeff= -Implpar*DtImpl/(Eps*Norm_V(jVar)*DtExpl)
        do i = 1, nCell
           iDiag = modulo(i-iStencil,3)+1
           do iVar = 1, nVar
              Matrix_VVCI(iVar,jVar,i,iDiag) = &
                   Coeff*(ResidEps_CV(i,iVar) - ResidOrig_CV(i,iVar))
           end do
        end do

     end do
  end do

  ! Add the diagonal part J = I - delta t*dR/dU
  do i=1, nCell
     do iVar = 1, nVar
        Matrix_VVCI(iVar,iVar,i,1) = Matrix_VVCI(iVar,iVar,i,1) + 1.0
     end do
  end do
 
  ! L-U decomposition
  call prehepta(nCell,nVar,nCell,nCell,0.0,&
       Matrix_VVCI(:,:,:,1), Matrix_VVCI(:,:,:,2), Matrix_VVCI(:,:,:,3))

  ! Put right hand side into a linear vector
  iX = 0
  do i = 1, nCell
     do iVar = 1, nVar
        iX=iX + 1
        x_I(iX) = RightHand_CV(i, iVar)
     end do
  end do
  ! x --> L^{-1}.rhs
  call Lhepta(nCell, nVar, nCell, nCell, x_I,&
       Matrix_VVCI(:,:,:,1), Matrix_VVCI(:,:,:,2))

  ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
  call Uhepta(.true.,nCell,nVar,nCell,nCell, x_I, Matrix_VVCI(:,:,:,3))

  ! Update the solution (x = U^n+1 - U^n)
  iX = 0
  do i = 1, nCell
     do iVar = 1, nVar
        iX=iX + 1
        State_GV(i, iVar) = State_GV(i, iVar) + x_I(iX)
     end do
  end do

  call update_boundary(nCell, nVar, State_GV)

  deallocate(RightHand_CV, ResidOrig_CV, StateEps_GV, ResidEps_CV, &
       Norm_V, x_I, Matrix_VVCI)

end subroutine PW_implicit

!=======================================================================
subroutine calc_residual(nOrder, Dt, nCell, nVar, State_GV, Resid_CV)

  implicit none

  integer, intent(in) :: nOrder, nCell, nVar
  real,    intent(in) :: Dt
  real,    intent(in) :: State_GV(-1:nCell+2, nVar)
  real,    intent(out):: Resid_CV(nCell, nVar)

  call CON_stop('calc_residual to be implemented')

end subroutine calc_residual

!=======================================================================
subroutine update_boundary(nCell, nVar, State_GV)

  implicit none

  integer, intent(in)    :: nCell, nVar
  real,    intent(inout) :: State_GV(-1:nCell+2, nVar)

  call CON_stop('update_boundary to be implemented')

end subroutine update_boundary
