module SP_ModIhData
  implicit none
  real,allocatable,dimension(:,:),save::State_VI
  real,allocatable,dimension(:,:),save::Xyz_DI
  integer::nPoint
  logical::DoInit=.true.
end module SP_ModIhData
!==============================================================================
subroutine sp_set_ihdata(nPointIn,XyzIn_DI)
  use SP_ModIhData
  use CON_world,ONLY:check_allocate
  implicit none
  integer,intent(in)::nPointIn
  real,dimension(3,nPoint),intent(in)::XyzIn_DI
  integer::iError
  nPoint=nPointIn
  
  if(DoInit)then
     allocate(Xyz_DI(3,nPoint),stat=iError)
     call check_allocate(iError,&
            'Xyz_DI in sp_set_ihdata')
     allocate(State_VI(8,nPoint),stat=iError)
     call check_allocate(iError,&
            'State_VI in sp_set_ihdata')
     DoInit=.false.
  end if
  Xyz_DI(:,1:nPoint)=XyzIn_DI(:,1:nPoint)
end subroutine sp_set_ihdata
!==============================================================================
subroutine SP_put_from_ih(nPartial,&
     iPutStart,&
     Put,&
     W,&
     DoAdd,&
     Buff_I,nVar)
  use CON_router
  use SP_ModIhData
  implicit none
  integer,intent(in)::nPartial,iPutStart,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::Buff_I
  integer::iCell
  real:: Weight
  iCell=Put%iCB_II(1,iPutStart)
  Weight=W%Weight_I(iPutStart)
  if(DoAdd)then
     State_VI(:,iCell)=State_VI(:,iCell)+Buff_I(:)*Weight
  else
     State_VI(:,iCell)=Buff_I(:)*Weight
  end if
end subroutine SP_put_from_ih
!==============================================================================

subroutine IH_get_for_sp(&
     nPartial,iGetStart,Get,W,State_V,nVar)
  !USES:
  use ModAdvance,ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK, &
       rho_, rhoUx_, rhoUy_, rhoUz_, Ux_, Uy_, Uz_, Bx_, By_, Bz_, P_

  use ModPhysics,ONLY:UnitSI_rho,UnitSI_p,UnitSI_U,UnitSI_B
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
  
  State_V(1)= Weight*State_VGB(rho_,i,j,k,iBlock)
  State_V(2:4)= Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)
  State_V(5)= Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
  State_V(6)= Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
  State_V(7)= Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
  State_V(8)= Weight*State_VGB(P_,i,j,k,iBlock)
  do iGet=iGetStart+1,iGetStart+nPartial-1
     i      = Get%iCB_II(1,iGet)
     j      = Get%iCB_II(2,iGet)
     k      = Get%iCB_II(3,iGet)
     iBlock = Get%iCB_II(4,iGet)
     Weight = W%Weight_I(iGet)
     State_V(1) = State_V(1) + &
          Weight*State_VGB(rho_,i,j,k,iBlock)
     State_V(2:4) =  State_V(2:4) + &
          Weight*State_VGB(rhoUx_:rhoUz_,i,j,k,iBlock)
     State_V(5) = State_V(5) + &
          Weight*(State_VGB(Bx_,i,j,k,iBlock)+ B0xCell_BLK(i,j,k,iBlock))
     State_V(6) = State_V(6) + &
          Weight*(State_VGB(By_,i,j,k,iBlock)+ B0yCell_BLK(i,j,k,iBlock))
     State_V(7) = State_V(7) + &
          Weight*(State_VGB(Bz_,i,j,k,iBlock)+ B0zCell_BLK(i,j,k,iBlock))
     State_V(8)= State_V(8) + &
          Weight*State_VGB(P_,i,j,k,iBlock)     
  end do
  ! Convert momentum to velocity and convert everything to SI units
  State_V(2:4) = State_V(2:4)/State_V(1)*UnitSI_U
  State_V(1)   = State_V(1)*UnitSI_rho
  State_V(5:7) = State_V(5:7)*UnitSI_B
  State_V(8)   = State_V(8)*UnitSI_p
end subroutine IH_get_for_sp
subroutine write_ihdata
  use SP_ModIhData
  use ModIoUnit
  use ModMain,ONLY:iCoupleSP,n_step,Time_Simulation
  implicit none
  character(LEN=14)::NameFile
  integer::iFile,iPoint,i
  write(NameFile,'(a,i4.4,a)')'ihdata',iCoupleSP,'.dat'
  iFile=io_unit_new()
  open(iFile,file=NameFile,status='unknown',&
          form='formatted')
  write(iFile,*)'Time_Simulation=',Time_Simulation,',  n_step=',n_step
  do iPoint=1,nPoint
     write(iFile,'(11(e13.6,2x))')Xyz_DI(:,iPoint),State_VI(1:8,iPoint)
  end do
  close(iFile)
end subroutine write_ihdata
