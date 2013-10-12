!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_integrator
  use CON_global_message_pass
  implicit none
  private !Except
  real::tNow,Dt
  real,dimension(:,:),pointer::Xyz_DI
  real,dimension(:,:),allocatable,save::XyzStored_DI
  public::tNow
  public::init_router_for_vector
  public::check_if_can_integrate
  public::advance_vector
  integer::nDim
  integer::nU_I(2) 
contains
  subroutine init_router_for_vector(&
       NameVector,&   !Global vector name, intent(in)
       SourceGD,& !GridDesctriptor for the source field 
       LineDD,&   !Domain Decomposition, intent(out)
       LineGD,&   !GirdDescriptor,save,intent(out)
       Router)    !resulting router, intent(out)
    character(LEN=*),intent(in)::NameVector
    type(GridDescriptorType),intent(in)::SourceGD
    type(DomainDecompositionType),intent(out)::LineDD
    type(GridDescriptorType),intent(out)::LineGD
    type(RouterType),intent(out)::Router
    !------------------------------------------------
    if(.not.is_local_grid(SourceGD%DD%Ptr))&
         call CON_stop(&
         'For a while init_router_for_vector is only used with local GDs')
    nU_I=ubound_vector(NameVector)
    call init_decomposition(LineDD,&
            compid_grid(SourceGD%DD%Ptr),1,&
            IsLocal=.true.)
    call get_root_decomposition(&
         LineDD,&                 !GridDescroptor to be constructed
         iRootMapDim_D=(/1/),&!The block amount, along each direction(D)
         XyzMin_D=(/cHalf/),&      !Minimal gen. coordinates, along each D 
         XyzMax_D=(/cHalf+nU_I(2)/),& !Maximal gen. coordinates, along each D
         nCells_D=(/nU_I(2)/))
    call set_standard_grid_descriptor(LineDD,&
         GridDescriptor=LineGD)
    call init_router(&
         SourceGD,& !GridDesctriptor for the source field (in) 
         LineGD,&   !GirdDescriptor,save,intent(out)
         Router,&   !resulting router, intent(out)
         nIndexesTarget=1)
  end subroutine init_router_for_vector
  !=============================================================!
  subroutine check_if_can_integrate(NameVector)
    character(LEN=*),intent(in)::NameVector
    integer:: nUNow_I(2),iError
    nU_I=ubound_vector(NameVector)
    if(allocated(XyzStored_DI))then
       nUNow_I=ubound(XyzStored_DI)
       if(all(nUNow_I>=nU_I))return
       nU_I=max(nU_I,nUNow_I)
       deallocate(XyzStored_DI)
    end if
    allocate(XyzStored_DI(nU_I(1),nU_I(2)),stat=iError)
    call check_allocate(iError,'XyzStored_DI')
  end subroutine check_if_can_integrate
  !=======================================================
  subroutine advance_vector(tStart,&
       tFinal,&
       NameVector,&
       NameMask,&
       d_xyz,&
       SourceGD,&
       LineGD,&
       Router)
    real,intent(in)::tStart,tFinal
    character(LEN=*),intent(in)::NameVector
    character(LEN=*),intent(in)::NameMask
    optional::NameMask
    interface
       subroutine d_xyz(nPartial,&
            iGetStart,&
            Get,&
            Weight,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::Weight
         real,dimension(nVar),intent(out)::Buff_I
       end subroutine d_xyz
    end interface
    type(GridDescriptorType),intent(in)::SourceGD,LineGD
    type(RouterType)::Router
    
    if(tFinal<=tStart)return
    call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
    nU_I=ubound(Xyz_DI)
    nDim=SourceGD%nDim
    Dt=(tFinal-tStart)*cHalf
    !  Dt=cTiny
    
    !First stage of the Runge-Kutta integration by time    
    call set_router(& 
         GridDescriptorSource=SourceGD,&
         GridDescriptorTarget=LineGD,&
         Router=Router,&
         NameMappingVector=NameVector,&
         NameMask=NameMask,&
         interpolate=interpolation_fix_reschange)
    XyzStored_DI(1:nDim,1:nU_I(2))=&
         Xyz_DI(1:nDim,1:nU_I(2))! Store Xyz_DI
    tNow=tStart
    call global_message_pass(Router,&
         nVar=nDim,&
         fill_buffer=d_xyz,&
         apply_buffer=put_d_xyz)
    call bcast_global_vector(&
         NameVector,&
         0,&
         Router%iComm,&
         NameMask)
    tNow=tNow+dt
    dt=dt+dt
    !Second stage of the Runge-Kutta integration by time    
    call set_router(& 
         GridDescriptorSource=SourceGD,&
         GridDescriptorTarget=LineGD,&
         Router=Router,&
         NameMappingVector=NameVector,&
         NameMask=NameMask,&
         interpolate=interpolation_fix_reschange)
    call global_message_pass(Router,&
         nVar=nDim,&
         fill_buffer=d_xyz,&
         apply_buffer=put_d_xyz)
    call bcast_global_vector(&
         NameVector,&
         0,&
         Router%iComm,&
         NameMask)
  end subroutine advance_vector
!=======================================================
  subroutine put_d_xyz(nPartial,&
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
    iCell=Put%iCB_II(1,iPutStart)
    if(DoAdd)then
       Xyz_DI(1:nDim,iCell)=&
            Xyz_DI(1:nDim,iCell)+dt*Buff_I(:)
    else
       Xyz_DI(1:nDim,iCell)=&
            XyzStored_DI(1:nDim,iCell)+dt*Buff_I(:)
    end if
  end subroutine put_d_xyz
  !================================================
end module CON_integrator
