subroutine PW_implicit(State_VC)

  use ModCommonVariables, ONLY: nDim
  use ModLinearSolver, ONLY: prehepta, Lhepta, Uhepta

  implicit none

  ! These should be in external modules !!!
  integer, parameter :: nVar = 10
  real :: DtExpl, DtImpl

  real, intent(inout) :: State_VC(nVar, nDim)


  real, parameter :: ImplPar = 1.0, Eps = 1.e-6

  real, dimension(nVar, nDim) :: StateOld_VC, RightHand_VC, ResidOrig_VC, &
       StateEps_VC, ResidEps_VC

  real :: Matrix_VVCI(nVar, nVar, nDim, 3)
       
  real :: Norm_V(nVar)         ! second norm of variables
  real :: x_I(nDim*nVar)       ! linear vector of right hand side/unknowns
  real :: Coeff

  integer :: iVar, jVar, i, iStencil, iDiag, iX

  !-----------------------------------------------------------------------

  ! calculate right hand side

  StateOld_VC = State_VC

  ! Note !! use explicit heat conduction
  call calc_residual(2, DtImpl, State_VC, RightHand_VC)

  ! Calculate the unperturbed residual with first order scheme
  call calc_residual(1, DtExpl, State_VC, ResidOrig_VC)

  ! Calculate the norm for the variables
  do iVar = 1, nVar
     Norm_V(iVar) = sqrt(sum(State_VC(iVar,:)**2)/nDim)
  end do

  ! Calculate the dR/dU matrix
  do jVar = 1, nVar
     do iStencil = 1, 3
        ! Get perturbed state
        StateEps_VC = StateOld_VC
        do i = iStencil, nDim, 3
           StateEps_VC(jVar, i) = StateEps_VC(jVar, i) + Eps*Norm_V(jVar)
        end do
        call update_boundary(StateEps_VC)
        call calc_residual(1, DtExpl, StateEps_VC, ResidEps_VC)

        Coeff= -Implpar*DtImpl/(Eps*Norm_V(jVar)*DtExpl)
        do i = 1, nDim
           iDiag = modulo(i-iStencil,3)+1
           do iVar = 1, nVar
              Matrix_VVCI(iVar,jVar,i,iDiag) = &
                   Coeff*(ResidEps_VC(iVar,i) - ResidOrig_VC(iVar,i))
           end do
        end do

     end do
  end do

  ! Add the diagonal part J = I - delta t*dR/dU
  do i=1, nDim
     do iVar = 1, nVar
        Matrix_VVCI(iVar,iVar,i,1) = Matrix_VVCI(iVar,iVar,i,1) + 1.0
     end do
  end do
 
  ! L-U decomposition
  call prehepta(nDim,nVar,nDim,nDim,0.0,Matrix_VVCI(:,:,:,1),&
       Matrix_VVCI(:,:,:,2), Matrix_VVCI(:,:,:,3))

  ! Put right hand side into a linear vector
  iX = 0
  do i = 1, nDim
     do iVar = 1, nVar
        iX=iX + 1
        x_I(iX) = RightHand_VC(iVar,i)
     end do
  end do
  ! x --> L^{-1}.rhs
  call Lhepta(nDim,nVar,nDim,nDim,x_I,Matrix_VVCI(:,:,:,1), Matrix_VVCI(:,:,:,2))

  ! x --> U^{-1}.L^{-1}.rhs = A^{-1}.rhs
  call Uhepta(.true.,nDim,nVar,nDim,nDim,x_I,Matrix_VVCI(:,:,:,3))

  ! Update the solution (x = U^n+1 - U^n)
  iX = 0
  do i = 1, nDim
     do iVar = 1, nVar
        iX=iX + 1
        State_VC(iVar, i) = StateOld_VC(iVar, i) + x_I(iX)
     end do
  end do

end subroutine PW_implicit

!=======================================================================
subroutine calc_residual(nOrder, DtIn, State_VC, Resid_VC)

  use ModCommonVariables, ONLY: nDim
  implicit none

  integer, parameter :: nVar = 10

  integer, intent(in) :: nOrder
  real,    intent(in) :: DtIn
  real,    intent(in) :: State_VC(nVar, nDim)
  real,    intent(out):: Resid_VC(nVar, nDim)

  call CON_stop('calc_residual to be implemented')

end subroutine calc_residual

!=======================================================================
subroutine update_boundary(State_VC)
  use ModCommonVariables, ONLY: nDim
  implicit none

  integer, parameter :: nVar = 10
  real, intent(inout) :: State_VC(nVar, nDim)

  call CON_stop('update_boundary to be implemented')

end subroutine update_boundary
