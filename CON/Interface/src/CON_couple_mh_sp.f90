!^CFG COPYRIGHT UM
Module CON_couple_mh_to_sp
  use CON_comp_param
  use CON_global_message_pass
  use ModConst
  implicit none
  private !Except
  public::init_couple_mh_to_sp
  public::couple_mh_to_sp
  type(RouterType),save,private::Router
  type(GridDescriptorType),save::IH_GridDescriptor !Source
  type(GridDescriptorType),save::SP_GridDescriptor !Target
  real,allocatable,dimension(:,:)::XyzTemp_DI
  real,dimension(:,:),pointer,save::Xyz_DI
  logical,dimension(:),pointer,save::Used_I
  integer,parameter::nPointMax=5000
  integer::nPoint
  public::nPoint
  integer::iPoint
  integer::iError
  real::BAndDXyz_I(1:6)!The interpolated values of full B and DXyz
  real::SignBDotR
  real::DsResolution,XyzLine_DI(3)
contains
  subroutine init_couple_mh_to_sp
    interface
       subroutine IH_get_a_line_point(nPartial,&
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
       end subroutine IH_get_a_line_point
       subroutine SC_get_a_line_point(nPartial,&
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
       end subroutine SC_get_a_line_point
    end interface
    !======================================================================


    call set_standard_grid_descriptor(IH_,GridDescriptor=&
         IH_GridDescriptor)
    call init_decomposition(SP_,SP_,1)
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         SP_,&                             !GridDescroptor to be constructed
         iRootMapDim_D=(/1/),&             !The block amount, along each direction(D)
         XyzMin_D=(/cHalf/),&              !Minimal gen. coordinates, along each D 
         XyzMax_D=(/cHalf+cOne/),& !Maximal gen. coordinates, along each D
         nCells_D=(/1/),&
         PE_I=(/0/))                     !PEs layout, throughout all the blocks

    call set_standard_grid_descriptor(SP_,GridDescriptor=&
         SP_GridDescriptor)

    call init_router(IH_GridDescriptor,SP_GridDescriptor,&
         Router)

    call IH_synchronize_refinement

    !Set the initial point for the line

    allocate(XyzTemp_DI(3,nPointMax),stat=iError)
    call check_allocate(iError,&
         'XyzTemp_DI in CON_couple_ih_to_sp')
    iPoint=1
    XyzTemp_DI(:,iPoint)=XyzLine_DI(:)
    !   XyzTemp_DI(1,iPoint)=15.
    !   XyzTemp_DI(2,iPoint)=0.
    !   XyzTemp_DI(3,iPoint)=0.2
    do while (dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint))>=cOne)
       if(iPoint>1)then
          if(dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint))>=&
               dot_product(XyzTemp_DI(:,iPoint-1),XyzTemp_DI(:,iPoint-1)))then
             nPoint=iPoint
             EXIT
          end if
       end if
       if(iPoint==nPointMax)call CON_stop(&
            'Insufficient number of points in couple_ih_sp =',iPoint)

       call set_router(& 
            GridDescriptorSource=IH_GridDescriptor,&
            GridDescriptorTarget=SP_GridDescriptor,&
            Router=Router,&
            mapping=xyz_point_mapping,&
            interpolate=interpolation_fix_reschange)     
       call global_message_pass(Router,&
            nVar=6,&
            fill_buffer=IH_get_a_line_point,&
            apply_buffer=SP_put_a_line_point)
       if(is_proc0(SP_))then
          if(iPoint==1)then
             SignBDotR=cOne
             XyzTemp_DI(:,iPoint+1)=XyzTemp_DI(:,iPoint)-DsResolution*&
                  cHalf*BAndDXyz_I(1:3)/sqrt(dot_product(&
                  BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
                  sqrt(dot_product(&
                  BAndDXyz_I(4:6),BAndDXyz_I(4:6)))          
             if(dot_product(XyzTemp_DI(:,iPoint+1),XyzTemp_DI(:,iPoint+1))>&
                  dot_product(XyzTemp_DI(:,iPoint),XyzTemp_DI(:,iPoint)))&
                  SignBDotR=-cOne
          end if
          XyzTemp_DI(:,iPoint+1)=XyzTemp_DI(:,iPoint)-DsResolution*&
               cHalf*BAndDXyz_I(1:3)/sqrt(dot_product(&
               BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
               sqrt(dot_product(&
               BAndDXyz_I(4:6),BAndDXyz_I(4:6)))*&
               SignBDotR
       end if

       nPoint=iPoint
       iPoint=iPoint+1
       call MPI_bcast(&
            XyzTemp_DI(1,iPoint),3,MPI_REAL,i_proc0(SP_),i_comm(),iError)
       call set_router(& 
            GridDescriptorSource=IH_GridDescriptor,&
            GridDescriptorTarget=SP_GridDescriptor,&
            Router=Router,&
            mapping=xyz_point_mapping,&
            interpolate=interpolation_fix_reschange)     
       call global_message_pass(Router,&
            nVar=6,&
            fill_buffer=IH_get_a_line_point,&
            apply_buffer=SP_put_a_line_point)
       if(is_proc0(SP_))then     
          XyzTemp_DI(:,iPoint)=XyzTemp_DI(:,iPoint-1)-DsResolution*&
               BAndDXyz_I(1:3)/sqrt(dot_product(&
               BAndDXyz_I(1:3),BAndDXyz_I(1:3)))*&
               sqrt(dot_product(&
               BAndDXyz_I(4:6),BAndDXyz_I(4:6)))*&
               SignBDotR
       end if
       call MPI_bcast(&
            XyzTemp_DI(1,iPoint),3,MPI_REAL,i_proc0(SP_),i_comm(),iError)
    end do
    call allocate_vector('SP_Xyz_DI',3,nPoint)
    call associate_with_global_vector(Xyz_DI,'SP_Xyz_DI')
    call allocate_mask('SP_IsInIH','SP_Xyz_DI')
    call associate_with_global_mask(Used_I,'SP_IsInIH')
    do iPoint=1,nPoint
       Xyz_DI(:,iPoint)= XyzTemp_DI(:,nPoint-iPoint+1)
       Used_I(iPoint)=.true.
    end do

    deallocate(XyzTemp_DI)
    if(is_proc0(SP_))&
         call get_root_decomposition(&
         SP_,&                            !GridDescroptor to be constructed
         iRootMapDim_D=(/1/),&            !The block amount, along each direction(D)
         XyzMin_D=(/cHalf/),&             !Minimal gen. coordinates, along each D 
         XyzMax_D=(/cHalf+real(nPoint)/),&!Maximal gen. coordinates, along each D
         nCells_D=(/nPoint/),&
         PE_I=(/0/))                    !PEs layout, throughout all the blocks
    SP_GridDescriptor%iGridPointMax_D=(/nPoint/)
    if(is_proc(SP_))call sp_set_ihdata(nPoint,Xyz_DI)
  end subroutine init_couple_mh_to_sp
  !========================================================================
  subroutine xyz_point_mapping(&
       SP_nDim,SP_Xyz_D,IH_nDim,IH_Xyz_D,IsInterfacePoint)
    use ModConst
    integer,intent(in)::SP_nDim,IH_nDim
    real,dimension(SP_nDim),intent(in)::SP_Xyz_D
    real,dimension(IH_nDim),intent(out)::IH_Xyz_D
    logical,intent(out)::IsInterfacePoint
    IH_Xyz_D=XyzTemp_DI(:,iPoint)
    IsInterfacePoint=.true.
  end subroutine xyz_point_mapping

  subroutine SP_put_a_line_point(nPartial,&
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
    real:: Weight
    Weight=W%Weight_I(iPutStart)
    if(DoAdd)then
       BAndDXyz_I=BAndDXyz_I+Buff_I(:)*Weight
    else
       BAndDXyz_I=Buff_I(:)*Weight
    end if
  end subroutine SP_put_a_line_point
 !==================================================================
  subroutine couple_mh_to_sp
    interface
       subroutine IH_get_for_sp(nPartial,&
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
       end subroutine IH_get_for_sp
    end interface
    interface
       subroutine SP_put_from_ih(nPartial,&
            iPutStart,&
            Put,&
            Weight,&
            DoAdd,&
            Buff_I,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iPutStart,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::Weight
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::Buff_I
       end subroutine SP_put_from_ih
    end interface
    if(.not.Router%IsProc)return
    call IH_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call bcast_global_vector('SP_Xyz_DI',&
         Router%iProc0Source,&
         Router%iCommUnion,&
         'SP_IsInIH')
    call set_router(& 
         GridDescriptorSource=IH_GridDescriptor,&
         GridDescriptorTarget=SP_GridDescriptor,&
         Router=Router,&
         NameMappingVector='SP_Xyz_DI',&
         NameMask='SP_IsInIH',&
         interpolate=interpolation_fix_reschange)     
    call global_message_pass(Router,&
         nVar=8,&
         fill_buffer=IH_get_for_sp,&
         apply_buffer=SP_put_from_ih)
    if(is_proc(SP_))then
       call sp_set_ihdata(nPoint,Xyz_DI(:,1:nPoint))
       call sp_smooth_ihdata
       if(is_proc0(SP_))call write_ihdata
    end if
    !=========================================================================
  end subroutine couple_mh_to_sp
  !===========================================================================
end Module CON_couple_mh_to_sp
