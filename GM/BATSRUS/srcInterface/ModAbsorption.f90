!^CFG COPYRIGHT UM
!BOP
!MODULE: ModDensityAndGradient - provide the density and gradients at SC_Xyz_DI
!INTERFACE:
module ModAbsorption
!USES:
  use MH_domain_decomposition
  use CON_global_message_pass
  use CON_integrator
  use ModMain, ONLY: nDim
  use ModProcMH, ONLY: iProc,iComm
  use ModMpi
  use ModAdvance, ONLY: State_VGB
  use ModVarIndexes
  use ModPhysics
  use ModConst
  !DESCRIPTION:
  !This file is an instantiation of the general advance_vector routine
  !as applied for extracting the data for the plasma density, its gradient
  !as well as the absorption coefficient at the circular frequency, Omega,
  !at the points which coordinates are in the g;pbal vector 'MH_Pos_DI',
  !where MH=MH,SC, or IH
  !EOP
  !BOP
  !REVISION HISTORY:
  !
  !EOP
  implicit none
  private !Except

  logical,save::DoInit=.true.

  character(LEN=10),save::NameVector,NameMask

  type(RouterType),save::Router
  type(GridDescriptorType),save::LineGrid,MhGrid
  type(DomainDecompositionType),save::LineDD

  real,allocatable,save,dimension(:)::Density_I,DeltaSNew_I
  real,allocatable,save,dimension(:,:)::GradDensity_DI
  real,allocatable,save,dimension(:)::AbsorptionCoeff_I

  !PUBLIC MEMBERS:
  public:: get_density_and_absorption, NameVector
  public:: GradDensity_DI, Density_I, DeltaSNew_I, AbsorptionCoeff_I
  real,public:: Omega  !Circular frequency, [rad/s]
  !EOP

contains
  subroutine calc_absorption(State_V, Absorption)
    !use CRASH_ModTransport
    use ModUser,ONLY: user_material_properties
    !The subroutine calculates the absorption coefficient, Absorption [m-1],
    ! at the circular frequency, omega.
    real,intent(in):: State_V(nVar)

    real, intent(out):: Absorption  !The absorption coefficient, [m-1]
    real:: NAtomicSI, ZAverage, TeSI, TeEV
    reaL:: Rho  !Mass density, SI
    !---------------------
    
    call user_material_properties(State_V=State_V, NAtomicOut=NAtomicSI, TeOut=TeSI,&
         AverageIonChargeOut= ZAverage)
  
  end subroutine calc_absorption
  !===============================
  subroutine get_density_and_absorption(nRay)
    !
    ! Calculates plasma density, Density_I, and its gradient, 
    ! GradDensity_DI(3,nRay), at specified locations Position_DI(3,nRay)
    ! Also, it provides appropriate step, DeltaSNew_I, conceivably dependent
    ! on the numeric grid size
    !
    integer,intent(in)::nRay
    integer::iError
 
    call timing_start('get_plasma_density')
    if(DoInit)then
       DoInit=.false.
       call set_standard_grid_descriptor(&
            MH_DomainDecomposition,GridDescriptor=MhGrid)
   
       call init_router_for_vector(&
            NameVector=NameVector,&
            SourceGD=MhGrid,&
            LineDD=LineDD,&
            LineGD=LineGrid,&
            Router=Router)
       call check_if_can_integrate(NameVector)
    end if
 
    call set_router(& 
         GridDescriptorSource=MhGrid,&
         GridDescriptorTarget=LineGrid,&
         Router=Router,&
         NameMappingVector=NameVector,&
         interpolate=interpolation_fix_reschange)
   
    call global_message_pass(Router=Router,&
         nVar=nDim+3,&    !Gradient of density (nDim) + 3 (density, DeltaS, absorption
         fill_buffer=get_density_local,&
         apply_buffer=put_density_value)

    call MPI_BCAST(GradDensity_DI(1,1),3*nRay,MPI_REAL,0,iComm,iError)
    call MPI_BCAST(Density_I(1),         nRay,MPI_REAL,0,iComm,iError)
    call MPI_BCAST(DeltaSNew_I(1),       nRay,MPI_REAL,0,iComm,iError)
    call MPI_BCAST(AbsorptionCoeff_I(1),       nRay,MPI_REAL,0,iComm,iError)
    !call bcast_global_vector(&
    !      NameVector,&
    !      0,&
    !      Router%iComm)
  
    call timing_stop('get_plasma_density')
  end subroutine get_density_and_absorption
  !================================================
  subroutine get_density_local(&
       nPartial,iGetStart,Get,W,State_V,nVar)
    !USES:
    use ModAdvance,ONLY: State_VGB,StateOld_VCB, &
         rho_
    use ModGeometry,ONLY:Dx_BLK,Dy_BLK,Dz_BLK
    use CON_router
    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    integer::iGet, i, j, k, iBlock
    real :: Weight, Absorption
    !----------------------------------------------------------
    ! call stop_MPI('+++stop_MPI: in get_density_local, '// &
    !     'before the first statement')
 
    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)
    State_V(1)= Weight*(&
         State_VGB(rho_,i+1,j,k,iBlock)-State_VGB(rho_,i-1,j,k,iBlock))&
         /(2*Dx_BLK(iBlock))
    State_V(2)= Weight*(&
         State_VGB(rho_,i,j+1,k,iBlock)-State_VGB(rho_,i,j-1,k,iBlock))&
         /(2*Dy_BLK(iBlock))
    State_V(3)= Weight*(&
         State_VGB(rho_,i,j,k+1,iBlock)-State_VGB(rho_,i,j,k-1,iBlock))&
         /(2*Dz_BLK(iBlock))
    State_V(nDim+1)= Weight*&
         State_VGB(rho_,i,j,k,iBlock)
    State_V(nDim+1+1)=Weight*&
         min(Dx_BLK(iBlock),Dy_BLK(iBlock),Dz_BLK(iBlock))
    !Save Absorption coeff
    call calc_absorption(State_VGB(:, i, j, k, iBlock), Absorption)
    State_V(nDim+3)=Weight*Absorption

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       State_V(1)= State_V(1)+Weight*(&
            State_VGB(rho_,i+1,j,k,iBlock)-State_VGB(rho_,i-1,j,k,iBlock))&
            /Dx_BLK(iBlock)
       State_V(2)= State_V(2)+Weight*(&
            State_VGB(rho_,i,j+1,k,iBlock)-State_VGB(rho_,i,j-1,k,iBlock))&
            /Dy_BLK(iBlock)
       State_V(3)= State_V(3)+Weight*(&
            State_VGB(rho_,i,j,k+1,iBlock)-State_VGB(rho_,i,j,k-1,iBlock))&
            /Dz_BLK(iBlock)
       State_V(nDim+1)  = State_V(nDim+1)+Weight*&
            State_VGB(rho_,i,j,k,iBlock)
       State_V(nDim+1+1)= State_V(nDim+1+1)+Weight*&
            min(Dx_BLK(iBlock),Dy_BLK(iBlock),Dz_BLK(iBlock))
       !Add Absorption coeff
       call calc_absorption(State_VGB(:, i, j, k, iBlock), Absorption)
       State_V(nDim+3) = State_V(nDim+3) + Weight*Absorption

    end do
  end subroutine get_density_local

  !====================================================================

  subroutine put_density_value(nPartial,&
       iPutStart,&
       Put,&
       W,&
       DoAdd,&
       Buff_I,nVar)
    implicit none
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer::iCell

    ! call stop_MPI('+++stop_MPI: in put_density_value, '// &
    !     'before the first statement')


    iCell=Put%iCB_II(1,iPutStart)
    if(DoAdd)then
       GradDensity_DI(:,iCell)= GradDensity_DI(:,iCell)+&
            Buff_I(1:nDim)
       Density_I(iCell)= Density_I(iCell)+&
            Buff_I(nDim+1)
       DeltaSNew_I(iCell) = DeltaSNew_I(iCell)+&
            Buff_I(nDim+1+1)
       !Add absorption coeff
       AbsorptionCoeff_I(iCell) =  AbsorptionCoeff_I(iCell) + &
            Buff_I(nDim+3)
    else
       GradDensity_DI(:,iCell)= &
            Buff_I(1:nDim)
       Density_I(iCell)= &
            Buff_I(nDim+1)
       DeltaSNew_I(iCell) = &
            Buff_I(nDim+1+1)
       
       !Put absorption coeff
       AbsorptionCoeff_I(iCell) = Buff_I(nDim+3)
       
    end if
  end subroutine put_density_value
end module ModAbsorption
