!==============================================================================
!BOP
!ROUTINE: GM_put_from_ih - tranfrom and put the data got from IH_
!INTERFACE:
subroutine GM_put_from_ih(nPartial,&
     iPutStart,&
     Put,& 
     Weight,&
     DoAdd,&
     StateSI_V,&
     nVar)
  !USES:
  use CON_coupler, ONLY: IndexPtrType, WeightPtrType
  use ModAdvance, ONLY: State_VGB,rho_,rhoUx_,rhoUy_,rhoUz_,Bx_,By_,Bz_,P_,&
       B0xCell_BLK, B0yCell_BLK, B0zCell_BLK

  use ModPhysics, nVarGm => nVar

  implicit none

  !INPUT ARGUMENTS:
  integer,intent(in)::nPartial,iPutStart,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::Weight
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::StateSI_V

  !REVISION HISTORY:
  !18JUL03     I.Sokolov <igorsok@umich.edu> - intial prototype/code
  !23AUG03                                     prolog
  !03SEP03     G.Toth    <gtoth@umich.edu>   - simplified
  !EOP

  character (len=*), parameter :: NameSub='GM_put_from_ih.f90'

  real,dimension(nVar)::State_V
  integer::iPut, i, j, k, iBlock
  !-----------------------------------------------------------------------
  if(nVar /= nVarGm)then
     write(*,*)NameSub//' from IH nVar=',nVar,' nVarGm=',nVarGm
     call CON_stop(NameSub//' ERROR: nVar /= nVarGm')
  end if

  State_V(rho_)          = StateSI_V(rho_)         / UnitSI_rho
  State_V(rhoUx_:rhoUz_) = StateSI_V(rhoUx_:rhoUz_)/(UnitSI_rho*UnitSI_U)
  State_V(Bx_:Bz_)       = StateSI_V(Bx_:Bz_)      / UnitSI_B
  State_V(P_)            = StateSI_V(P_)           / UnitSI_p

  i      = Put%iCB_II(1,iPutStart)
  j      = Put%iCB_II(2,iPutStart)
  k      = Put%iCB_II(3,iPutStart)
  iBlock = Put%iCB_II(4,iPutStart)

  if(DoAdd)then
     State_VGB(:,i,j,k,iBlock) = State_VGB(:,i,j,k,iBlock) + State_V
  else
     State_VGB(:,i,j,k,iBlock)   = State_V
     State_VGB(Bx_,i,j,k,iBlock) = State_V(Bx_) - B0xCell_BLK(i,j,k,iBlock)
     State_VGB(By_,i,j,k,iBlock) = State_V(By_) - B0yCell_BLK(i,j,k,iBlock)
     State_VGB(Bz_,i,j,k,iBlock) = State_V(Bz_) - B0zCell_BLK(i,j,k,iBlock)
  end if
end subroutine GM_put_from_ih
!==============================================================================
