!^CFG COPYRIGHT UM
!BOP
!MODULE: ModLagrangianGrid - advects the fragment of SP grid in IH or SC
!INTERFACE:  
module ModLagrangianGrid
  !USES:
  use MH_domain_decomposition
  use CON_global_message_pass
  use ModMain, ONLY: NameThisComp
  use CON_comp_param,ONLY:NameComp_I
  !DESCRIPTION:
  !This file is an instantiation of the general advance_vector routine
  !as applied for advection of lagrangian points
  !NameVector here can be 'IM_Xyz_DI' or
  !'SP_Xyz_DI', the coordinates of the points along the chosen 
  !magnetic field lines. They are advected according to the equation
  !dX/dt=u(X,t), the velocity field is obtained from the MHD solution
  !EOP
  !BOP
  !REVISION HISTORY:
  !I.V.Sokolov<igorsok@umich.edu>- 6/30/04-prototype,code,prolog
  !EOP
  implicit none
  private !Except
  logical,save::DoSkip=.false.
  character(LEN=10),save::NameVector,NameMask
  type(RouterType),save::Router
  type(GridDescriptorType),save::LineGrid,MhGrid
  type(DomainDecompositionType),save::LineDD
  !BOP
  !PUBLIC MEMBERS:
  public:: advance_lagrangian_grid
  !EOP
contains
!=======================================================================
  subroutine advance_lagrangian_grid(tStart,tFinal)
    use CON_integrator,ONLY:&
         init_router_for_vector,&
         check_if_can_integrate,&
         advance_vector
    real,intent(in)::tStart,tFinal
    integer::iError,iPoint,iCompLine
    logical,save::DoInit=.true.
    !------------------------------------
    if(DoSkip)return
    call timing_start('lagrangian_grid')
    if(DoInit)then
       DoInit=.false.
       DoSkip=.true.
       
       !Figure out if the component needs to advect 
       !lagrangian points
       do iCompLine=1,MaxComp
          if(used_mask(&
               NameComp_I(iCompLine)&
               //'_IsIn'//&
               NameThisComp)&
               )then
             !Initialize
             NameVector=NameComp_I(iCompLine)&
                  //'_Xyz_DI'
             NameMask=NameComp_I(iCompLine)&
                  //'_IsIn'//&
                  NameThisComp
             call set_standard_grid_descriptor(&
                  MH_DomainDecomposition,GridDescriptor=MhGrid)
             call init_router_for_vector(&
                  NameVector=NameVector,&
                  SourceGD=MhGrid,&
                  LineDD=LineDD,&
                  LineGD=LineGrid,&
                  Router=Router)
             call check_if_can_integrate(NameVector)
             DoSkip=.false.
             !End of initialization
          end if
       end do
       if(DoSkip)then
          call timing_stop('lagrangian_grid')
          return !Nothing to advect in this component
       end if
    end if
    call  advance_vector(tStart,&
         tFinal,&
         NameVector,&
         NameMask,&
         d_xyz=get_u,&
         SourceGD=MhGrid,&
         LineGD=LineGrid,&
         Router=Router)
    call timing_stop('lagrangian_grid')
  end subroutine advance_lagrangian_grid
  subroutine get_u(&
       nPartial,iGetStart,Get,W,State_V,nVar)
    !USES:
    use ModPhysics,ONLY: unitSI_t
    use ModAdvance,ONLY: State_VGB,StateOld_VCB, &
         rho_,rhoUx_,rhoUz_
    use ModMain,ONLY:Time_Simulation,dt
    use CON_router
    use CON_integrator,ONLY:tNow
    
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
    
    State_V(1:3)= Weight*(&
         (Time_Simulation-tNow)*&
         StateOld_VCB(RhoUx_:rhoUz_,i,j,k,iBlock)/&
         StateOld_VCB(rho_,i,j,k,iBlock)+&
         (tNow+dt*unitSI_t-Time_Simulation)*&
         State_VGB(RhoUx_:rhoUz_,i,j,k,iBlock)/&
         State_VGB(rho_,i,j,k,iBlock))/(dt*unitSI_t**2)
    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       State_V(1:3)=  State_V(1:3)+Weight*(&
            (Time_Simulation-tNow)*&
            StateOld_VCB(RhoUx_:rhoUz_,i,j,k,iBlock)/&
            StateOld_VCB(rho_,i,j,k,iBlock)+&
            (tNow+dt*unitSI_t-Time_Simulation)*&
            State_VGB(RhoUx_:rhoUz_,i,j,k,iBlock)/&
            State_VGB(rho_,i,j,k,iBlock))/(dt*unitSI_t**2)
    end do
  end subroutine get_u
end module ModLagrangianGrid
!============interface to BATS_methods====================
subroutine update_lagrangian_grid(tStart,tFinal)
  use ModLagrangianGrid,ONLY:advance_lagrangian_grid
  implicit none
  real,intent(in)::tStart,tFinal
  call advance_lagrangian_grid(tStart,tFinal)
end subroutine update_lagrangian_grid
