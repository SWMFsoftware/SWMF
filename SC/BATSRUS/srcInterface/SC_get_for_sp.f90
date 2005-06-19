!==============================================================================
subroutine SC_get_for_sp(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  !USES:
  use SC_ModAdvance,ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, &
       rho_, rhoUx_, rhoUy_, rhoUz_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_

  use SC_ModPhysics,ONLY:UnitSI_rho,UnitSI_p,UnitSI_U,UnitSI_B,unitSI_x
  use CON_router
  use SC_ModGeometry, ONLY: x_BLK, y_BLK, z_BLK
  implicit none

  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V

  integer::iGet, i, j, k, iBlock
  real :: Weight
  !The meaning of state intdex in buffer and in model can be 
  !different. Below are the conventions for buffer:
  integer,parameter::&
       BuffRho_  =1,&
       BuffUx_=2,&
       BuffUz_=4,&
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
  
  State_V(BuffRho_)= Weight*State_VGB(rho_,i,j,k,iBlock)
  State_V(BuffUx_:BuffUz_)= Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)/&
       State_VGB(rho_,i,j,k,iBlock)
  State_V(BuffBx_)= Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
  State_V(BuffBy_)= Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
  State_V(BuffBz_)= Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
  State_V(BuffP_)= Weight*State_VGB(P_,i,j,k,iBlock)
  State_V(BuffX_:BuffZ_)     = &
       (/x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)/)&
      *Weight
  do iGet=iGetStart+1,iGetStart+nPartial-1
     i      = Get%iCB_II(1,iGet)
     j      = Get%iCB_II(2,iGet)
     k      = Get%iCB_II(3,iGet)
     iBlock = Get%iCB_II(4,iGet)
     Weight = W%Weight_I(iGet)
     State_V(1) = State_V(1) + &
          Weight*State_VGB(rho_,i,j,k,iBlock)
     State_V(BuffUx_:BuffUz_) =  State_V(BuffUx_:BuffUz_) + &
          Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)/&
          State_VGB(rho_,i,j,k,iBlock)
     State_V(BuffBx_) = State_V(BuffBx_) + &
          Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
     State_V(BuffBy_) = State_V(BuffBy_) + &
          Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
     State_V(BuffBz_) = State_V(BuffBz_) + &
          Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
     State_V(BuffP_)= State_V(BuffP_) + &
          Weight*State_VGB(P_,i,j,k,iBlock)
     State_V(BuffX_:BuffZ_)     =State_V(BuffX_:BuffZ_)+ &
       (/x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), z_BLK(i,j,k,iBlock)/)&
      *Weight     
  end do
  ! Convert momentum to velocity and convert everything to SI units
  State_V(BuffUx_:BuffUz_) = State_V(BuffUx_:BuffUz_)*UnitSI_U
  State_V(1)   = State_V(1)*UnitSI_rho
  State_V(BuffBx_:BuffBz_) = State_V(BuffBx_:BuffBz_)*UnitSI_B
  State_V(BuffP_)   = State_V(BuffP_)*UnitSI_p
  State_V(BuffX_:BuffZ_) = State_V(BuffX_:BuffZ_)    *UnitSI_x 
end subroutine SC_get_for_sp
!==================================================================
subroutine SC_get_a_line_point(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  !USES:
  use SC_ModAdvance,ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, &
       Bx_, By_, Bz_
  use SC_ModGeometry,ONLY:dx_BLK,dy_BLK,dz_BLK
  use CON_router
  
  implicit none
  
  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iGetStart,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
  real,dimension(nVar),intent(out)::State_V
  
  integer::iGet, i, j, k, iBlock
  real :: Weight
  !----------------------------------------------------------
  
  i      = Get%iCB_II(1,iGetStart)
  j      = Get%iCB_II(2,iGetStart)
  k      = Get%iCB_II(3,iGetStart)
  iBlock = Get%iCB_II(4,iGetStart)
  Weight = W%Weight_I(iGetStart)

  State_V(1)= Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
  State_V(2)= Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
  State_V(3)= Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
  State_V(4)= dx_BLK(iBlock)*Weight
  State_V(5)= dy_BLK(iBlock)*Weight
  State_V(6)= dz_BLK(iBlock)*Weight
  do iGet=iGetStart+1,iGetStart+nPartial-1
     i      = Get%iCB_II(1,iGet)
     j      = Get%iCB_II(2,iGet)
     k      = Get%iCB_II(3,iGet)
     iBlock = Get%iCB_II(4,iGet)
     Weight = W%Weight_I(iGet)
     State_V(1)=  State_V(1)+ Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
     State_V(2)=  State_V(2)+ Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
     State_V(3)=  State_V(3)+ Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
     State_V(4)=  State_V(4)+ dx_BLK(iBlock)*Weight
     State_V(5)=  State_V(5)+ dy_BLK(iBlock)*Weight
     State_V(6)=  State_V(6)+ dz_BLK(iBlock)*Weight
  end do
end subroutine SC_get_a_line_point


