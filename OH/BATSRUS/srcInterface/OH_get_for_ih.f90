!^CFG COPYRIGHT UM
subroutine OH_get_for_ih(&
     nPartial,iGetStart,Get,W,State_V,nVar)

  !USES:
  use OH_ModAdvance,ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, &
       rho_, rhoUx_, rhoUy_, rhoUz_, Bx_, By_, Bz_,P_
       
  use OH_ModPhysics, ONLY: No2Si_V, UnitRho_, UnitP_, UnitRhoU_, UnitB_, UnitX_
  use CON_router

  use OH_ModGeometry, ONLY: x_BLK, y_BLK, z_BLK

  implicit none

  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V

  integer::iGet, i, j, k, iBlock
  real :: Weight

  character (len=*), parameter :: NameSub='OH_get_for_ih'
  !The meaning of state intdex in buffer and in model can be 
  !different. Below are the conventions for buffer:
  integer,parameter::&
       BuffRho_  =1,&
       BuffRhoUx_=2, BuffUx_=BuffRhoUx_, &
       BuffRhoUz_=4, BuffUz_=BuffRhoUz_, &
       BuffBx_   =5,&
       BuffBy_   =6,&
       BuffBz_   =7,&
       BuffP_    =8,&
       BuffX_    =9,BuffZ_=11
 

  !----------------------------------------------------------
 
  i      = Get%iCB_II(1,iGetStart)
  j      = Get%iCB_II(2,iGetStart)
  k      = Get%iCB_II(3,iGetStart)
  iBlock = Get%iCB_II(4,iGetStart)
  Weight = W%Weight_I(iGetStart)
  
  State_V(BuffRho_)          = &
       State_VGB(rho_,         i,j,k,iBlock) *Weight
  State_V(BuffRhoUx_:BuffRhoUz_) = &
       State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock) *Weight
  State_V(BuffBx_)           = &
       (State_VGB(Bx_,          i,j,k,iBlock) + &
       B0xCell_BLK(i,j,k,iBlock))*Weight
  State_V(BuffBy_)           = &
       (State_VGB(By_,          i,j,k,iBlock) + &
       B0yCell_BLK(i,j,k,iBlock))*Weight
  State_V(BuffBz_)           = &
       (State_VGB(Bz_,          i,j,k,iBlock) + &
       B0zCell_BLK(i,j,k,iBlock))*Weight
  State_V(BuffP_)            = &
       State_VGB(P_,       i,j,k,iBlock) *Weight
  State_V(BuffX_:BuffZ_)     = &
       (/x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)/)&
      *State_VGB(rho_,         i,j,k,iBlock)*Weight

  do iGet=iGetStart+1,iGetStart+nPartial-1
     i      = Get%iCB_II(1,iGet)
     j      = Get%iCB_II(2,iGet)
     k      = Get%iCB_II(3,iGet)
     iBlock = Get%iCB_II(4,iGet)
     Weight = W%Weight_I(iGet)
     State_V(BuffRho_)             =State_V(BuffRho_)             +&
          State_VGB(Rho_,i,j,k,iBlock) *Weight 
     State_V(BuffRhoUx_:BuffRhoUz_)=State_V(BuffRhoUx_:BuffRhoUz_)+&
          State_VGB(RhoUx_:rhoUz_,i,j,k,iBlock) *Weight
     State_V(BuffBx_)              =State_V(BuffBx_)              +&
          (State_VGB(Bx_,i,j,k,iBlock) + &
          B0xCell_BLK(i,j,k,iBlock))*Weight
     State_V(BuffBy_)              =State_V(BuffBy_)              +&
          (State_VGB(By_,i,j,k,iBlock) + &
          B0yCell_BLK(i,j,k,iBlock))*Weight
     State_V(BuffBz_)              =State_V(BuffBz_)              +&
          (State_VGB(Bz_,i,j,k,iBlock) + &
          B0zCell_BLK(i,j,k,iBlock))*Weight
     State_V(BuffP_)               =State_V(BuffP_)               +&
          State_VGB(P_,i,j,k,iBlock) *Weight
     State_V(BuffX_:BuffZ_)        = State_V(BuffX_:BuffZ_)       +&
          (/x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)/)&
          *State_VGB(rho_,i,j,k,iBlock)*Weight
  end do
  
  ! Convert to SI units
  State_V(BuffRho_)             = State_V(BuffRho_)       *No2Si_V(UnitRho_)
  State_V(BuffRhoUx_:BuffRhoUz_)= &
       State_V(BuffRhoUx_:BuffRhoUz_)                     *No2Si_V(UnitRhoU_)
  State_V(BuffBx_:BuffBz_)      = State_V(BuffBx_:BuffBz_)*No2Si_V(UnitB_)
  State_V(BuffP_)               = State_V(BuffP_)         *No2Si_V(UnitP_)
  State_V(BuffX_:BuffZ_)        = &
       State_V(BuffX_:BuffZ_)          * No2Si_V(UnitRho_)*No2Si_V(UnitX_) 


end subroutine OH_get_for_ih
