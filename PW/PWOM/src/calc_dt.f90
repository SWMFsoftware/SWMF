subroutine calc_dt
  use ModCommonVariables, ONLY: State_GV,StateOld_GV,nDim,nIon,iRho_I,iP_I,Dt
  use ModNumConst,        ONLY: cHalf
  use ModPWOM,            ONLY: DtVertical
  implicit none
  real, parameter   :: RatioMin = 0.8, RatioMax = 1.2, cIncrease = 1.2
  real, parameter   :: DtMin = 1.0e-1
  real, allocatable :: pRatioState_CV(:,:),RhoRatioState_CV(:,:)
  real :: RhoRatioMin, RhoRatioMax, pRatioMin, pRatioMax
  integer :: iIon
  !----------------------------------------------------------------------------
  
  if (.not.allocated(pRatioState_CV)) &
       allocate(pRatioState_CV(nDim,nIon),RhoRatioState_CV(nDim,nIon))

  ! Calculate the largest relative change in p and rho state variables
  do iIon=1,nIon
     pRatioState_CV(1:nDim,iIon) = &
          State_GV(1:nDim,iP_I(iIon)) &
          / StateOld_GV(1:nDim,iP_I(iIon))
     RhoRatioState_CV(1:nDim,iIon) = &
          State_GV(1:nDim,iRho_I(iIon)) &
          / StateOld_GV(1:nDim,iRho_I(iIon))
  enddo
  
  ! Find the Min and Max of the changes
  pRatioMin   = minval(pRatioState_CV)  ;pRatioMax   = maxval(pRatioState_CV)
  RhoRatioMin = minval(RhoRatioState_CV);RhoRatioMax = maxval(RhoRatioState_CV)

!  write(*,*) 'pRatioMin,pRatioMax',pRatioMin,pRatioMax
!  write(*,*) 'RhoRatioMin,RhoRatioMax',RhoRatioMin,RhoRatioMax

  ! Change Dt based on change from timestep to timestep
  if (pRatioMin < RatioMin .or. pRatioMax > RatioMax &
       .or. RhoRatioMin < RatioMin .or. RhoRatioMax > RatioMax) then
!     write(*,*) 'reducing orig timestep, Dt = ',Dt 
     Dt = Dt * cHalf
     Dt = max(Dt,DtMin)
     Dt = min(Dt,DtVertical)
!     write(*,*) 'reducing timestep, Dt = ',Dt
  else
!     write(*,*) 'increasing orig timestep, Dt = ',Dt 
     Dt = Dt * cIncrease
     Dt = min(Dt,DtVertical)
     Dt = max(Dt,DtMin)
!     write(*,*) 'increasing timestep, Dt = ',Dt
  endif
  deallocate(RhoRatioState_CV,pRatioState_CV)
end subroutine calc_dt
