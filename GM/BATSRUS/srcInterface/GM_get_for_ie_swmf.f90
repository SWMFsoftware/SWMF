!^CFG COPYRIGHT UM
!^CMP FILE IE
!==================================================================!
!BOP
!ROUTINE: GM_get_for_ie - gets FACurrents for coupling via SWMF
!INTERFACE:
subroutine GM_get_for_ie_swmf(&
     nPartial,iGet,Get,W,State_V,nVar)
!USES:
  use ModMain, ONLY    : TypeCoordSystem, Time_Simulation
  use CON_physics, ONLY: get_planet_field, transform_matrix
  use CON_coupler,ONLY : IndexPtrType,WeightPtrType,cOne
  use ModGeometry,ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics,ONLY  : UnitSI_J,UnitSI_X, RBody
  use ModAdvance,ONLY  : State_VGB, Bx_, By_, Bz_
  use ModMappingParam,ONLY:rMap,DipoleSign
  implicit none
!INPUT ARGUMENTS:  
  integer,intent(in)::nPartial,iGet,nVar
  type(IndexPtrType),intent(in)::Get
  type(WeightPtrType),intent(in)::W
!OUTPUT ARGUMENTS:
  real,dimension(nVar),intent(out)::State_V
!REVISION HISTORY:
!14AUG03 - Gabor Toth <gtoth@umich,edu> and Aaron Ridley - 
!          initial prototype/code for coupling via MPI
!21Aug03 - I.Sokolov <igorsok@umich.edu> - major revision for
!          coupling via SWMF/prolog
!EOP
  real,dimension(3)::XyzGm_D,XyzIe_D,XyzSmg_D
  real,dimension(3)::FieldGm_D,FieldIe_D
  real::dir
  integer::i,j,k,iBLK
  
  i=Get%iCB_II(1,iGet)
  j=Get%iCB_II(2,iGet)
  k=Get%iCB_II(3,iGet)
  iBLK=Get%iCB_II(4,iGet)
  
  State_V(1) = 0.50* &              
       ((State_VGB(Bz_,i,j+1,k,iBLK) - &
       State_VGB(Bz_,i,j-1,k,iBLK))/dy_BLK(iBLK) - &
       (State_VGB(By_,i,j,k+1,iBLK) - &
       State_VGB(By_,i,j,k-1,iBLK))/dz_BLK(iBLK))
  
  State_V(2) =  0.50* &
       ((State_VGB(Bx_,i,j,k+1,iBLK) - &
       State_VGB(Bx_,i,j,k-1,iBLK))/dz_BLK(iBLK) - &
       (State_VGB(Bz_,i+1,j,k,iBLK) - &
       State_VGB(Bz_,i-1,j,k,iBLK))/dx_BLK(iBLK))
  
  State_V(3) = 0.50* &
       ((State_VGB(By_,i+1,j,k,iBLK) - &
       State_VGB(By_,i-1,j,k,iBLK))/dx_BLK(iBLK) - &
       (State_VGB(Bx_,i,j+1,k,iBLK) - &
       State_VGB(Bx_,i,j-1,k,iBLK))/dy_BLK(iBLK))
  
  !Transformation to the SI system:
  State_V(1:3)=State_V(1:3)*UnitSI_J
  if(nVar==3)return
  
  XyzGm_D(1)=x_BLK(i,j,k,iBLK)
  XyzGm_D(2)=y_BLK(i,j,k,iBLK)
  XyzGm_D(3)=z_BLK(i,j,k,iBLK)
  
  ! Transform input coordinates to SMG system (used by IE)
  XyzSmg_D = matmul(transform_matrix(&
       Time_Simulation,TypeCoordSystem,'SMG'),&
       XyzGm_D)
  ! Check if the initial point is on the northern hemisphere
  If(XyzSmg_D(3)>=0)then
     call Get_Mapping_Point(XyzGm_D, 1.0, XyzIe_D)
     dir=DipoleSign
  else
     call Get_Mapping_Point(XyzGm_D, -1.0, XyzIe_D)
     dir=-DipoleSign
  end If
  call Get_MagneticField_Orient(XyzGm_D,State_V(4:6))
  call get_planet_field(Time_Simulation,XyzGm_D,&
       TypeCoordSystem//' NORM',FieldGm_D)
  call get_planet_field(Time_Simulation,XyzIe_D,&
       'SMG NORM', FieldIe_D)
  ! Get the ratio between the magnitude of the ionosphere B 
  ! and the magnetosphere B
  State_V(7)=dir*sqrt(dot_product(FieldIe_D,FieldIe_D)/&
       dot_product(FieldGm_D,FieldGm_D) )
  
  ! Get the absolute value of the cos of the angle
  !  between the magnetic field and radial direction
  State_V(8) = &
       abs(sum(FieldIe_D*XyzIe_D))/&
       (rMap*sqrt(dot_product(FieldIe_D,FieldIe_D)))
  if(nVar==11)State_V(9:11)=XyzIe_D
end subroutine GM_get_for_ie_swmf
!=====================================================================
