!^CMP COPYRIGHT UM
!^CMP FILE IH
!^CMP FILE SC
!BOP
!MODULE: CON_couple_ih_sc - couple IH to SC one way
!INTERFACE:
module CON_couple_ih_sc
!DESCRIPTION:
!This coupler is  based on the global message  pass.
!The toolkit for coupling the different data sets within a
!single component code, or within a framework as it is now 
!
!The features of the present coupler: the domain decompositions
!both for a target and for a source are with AMR.
!IH is a source, SC is a target

  !USES:
  use CON_coupler
  use ModConst

  implicit none
  private !except

  !
  !PUBLIC MEMBER FUNCTIONS:
  public:: couple_ih_sc_init
  public:: couple_ih_sc,couple_sc_ih

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V.<igorsok@umich.edu> - prototype for ih-gm
  ! 7/04/04                                 - version for ih-sc
  ! 7/20/04                                 - version for sc-buffer
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: IH_iGridRealization=-2
  integer :: SC_iGridInScIh=-2
  integer :: SC_iGridInIhSc=-2
  
 
  type(RouterType),save             :: RouterIhSc
  type(RouterType),save             :: RouterScBuff
  type(GridDescriptorType),save     :: SC_SourceGrid
  type(GridDescriptorType),save     :: SC_TargetGrid
  type(GridDescriptorType),save     :: IH_Grid
  type(GridDescriptorType),save     :: BuffGD
  type(DomainDecompositionType),&
       save,target                  :: BuffDD
  logical :: DoInitialize=.true., DoTest, DoTestMe

  character(len=*), parameter :: NameMod='couple_ih_sc'


contains

  !===============================================================!
  subroutine couple_ih_sc_init
    interface
       subroutine IH_set_buffer_grid(DD)
         use CON_domain_decomposition
         implicit none
         type(DomainDecompositionType),&
              intent(out)::DD
       end subroutine IH_set_buffer_grid
    end interface

    if(.not.DoInitialize)return
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
 !   if(DoTest)write(*,*)'couple_ih_sc_init iProc=',i_proc()

    call init_coupler(              &    
       iCompSource=IH_,             & ! component index for source
       iCompTarget=SC_,             & ! component index for target
       nGhostPointTarget=2,         & ! number of halo points in target
       GridDescriptorSource=IH_Grid,& ! OUT!\
       GridDescriptorTarget=SC_TargetGrid,& !-General coupler variables 
       Router=RouterIhSc)             ! OUT!/

    IH_iGridRealization=-1
    SC_iGridInIhSc     =-1
    SC_iGridInScIh     =-1
    
    call IH_set_buffer_grid(BuffDD)
    call set_standard_grid_descriptor(&
         BuffDD,          &
         Standard_=Nodes_,&
         nGhostGridPoints=1,  &
         GridDescriptor=BuffGD)
    call set_standard_grid_descriptor(&
         SC_,GridDescriptor=SC_SourceGrid)
    call init_buffer_grid_couple(&
         SourceGD=SC_SourceGrid,&
         TargetGD=BuffGD, &
         RouterToBuffer=RouterScBuff,    &
         nVar=8,&
         NameBuffer='IH_from_sc')  !Version for the first order in time
  end subroutine couple_ih_sc_init
  !===============================================================!
!BOP
!IROUTINE: couple_ih_sc - interpolate and get MHD state at SC outer ghostpoints
!INTERFACE:
  subroutine couple_ih_sc(TimeCoupling)
    !INPUT ARGUMENTS:
    interface
       subroutine SC_put_from_ih(nPartial,&
            iPutStart,&
            Put,& 
            Weight,&
            DoAdd,&
            StateSI_V,&
            nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iPutStart,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::Weight
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::StateSI_V
       end subroutine SC_put_from_ih
       subroutine IH_get_for_sc(&
            nPartial,iGetStart,Get,W,State_V,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::W
         real,dimension(nVar),intent(out)::State_V
       end subroutine IH_get_for_sc
    end interface

    real,intent(in)::TimeCoupling
!EOP
    if(.not.RouterIhSc%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)
!    if(DoTest)write(*,*)'couple_ih_sc iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(RouterIhSc%iProc0Source,&
                                   RouterIhSc%iCommUnion)
    call SC_synchronize_refinement(RouterIhSc%iProc0Target,&
                                   RouterIhSc%iCommUnion)


    if(IH_iGridRealization/=i_realization(IH_).or.&     
         SC_iGridInIhSc/=i_realization(SC_))then  

!       if(DoTest)write(*,*)'couple_ih_sc call set_router iProc=',i_proc()

       call set_router(&
            GridDescriptorSource=IH_Grid,&
            GridDescriptorTarget=SC_TargetGrid,&
            Router=RouterIhSc,&
            is_interface_block=boundary_block,&
            interface_point_coords=outer_cells, &
            interpolate=interpolation_fix_reschange)
       IH_iGridRealization=i_realization(IH_)
       SC_iGridInIhSc       =i_realization(SC_)
 !      if(DoTest)write(*,*)'couple_ih_sc passed set_router iProc=',i_proc()
    end if

    call couple_comp(&
         RouterIhSc,&
         nVar=8,&
         fill_buffer=IH_get_for_sc,&
         apply_buffer=SC_put_from_ih)

  end subroutine couple_ih_sc
  !======================================================!
  logical function boundary_block(lGlobalTreeNode)

    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsBoundary_D

    IsBoundary_D=is_right_boundary_d(&
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode).or.&
         is_left_boundary_d(&
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    boundary_block=any(IsBoundary_D)

  end function Boundary_block
  !===============================================================!
  subroutine outer_cells(&
       GridDescriptor,&
       lGlobalTreeNode,&
       Xyz_D,&
       nIndexes,&
       Index_I,&
       IsInterfacePoint)

    type(GridDescriptorType),intent(in):: GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    logical,intent(out)::IsInterfacePoint
    real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D
    integer, dimension(nIndexes),intent(inout)::Index_I
    logical,dimension(3)::&
         IsLeftFace_D,IsRightFace_D
    integer,parameter::x_=1,y_=2,z_=3

    IsLeftFace_D=Index_I(x_:z_)<1.and.is_left_boundary_d(&
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    IsRightFace_D=Index_I(x_:z_)>&
         ncells_decomposition_d(SC_TargetGrid%DD%Ptr).and.&
         is_right_boundary_d(SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    IsInterfacePoint=any(IsRightFace_D.or.IsLeftFace_D)
  end subroutine outer_cells 
  !========================================================!
!===============================================================!
  !BOP
!IROUTINE: couple_sc_ih - interpolate and get MHD state at the IH buffer grid 
!INTERFACE:
  subroutine couple_sc_ih(TimeCoupling)
    use ModIoUnit
    !INPUT ARGUMENTS:
    interface
       subroutine SC_get_for_ih(&
            nPartial,iGetStart,Get,W,State_V,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::W
         real,dimension(nVar),intent(out)::State_V
       end subroutine SC_get_for_ih
    end interface
    integer::iPoint,nU_I(2)
    real,intent(in)::TimeCoupling
    integer,save::iCoupling=0
    integer::iFile
    character(LEN=21)::NameFile
!EOP
    if(.not.RouterScBuff%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)
!    if(DoTest)write(*,*)'couple_sc_ih iProc=',i_proc()

    if(SC_iGridInScIh/=SC_iGridInIhSc)then  

!       if(DoTest)write(*,*)'couple_sc_ih call set_router iProc=',i_proc()

       call set_router(&
            GridDescriptorSource=SC_SourceGrid,&
            GridDescriptorTarget=BuffGD,&
            Router=RouterScBuff,&
            mapping=buffer_grid_point,&
            interpolate=interpolation_fix_reschange)
!       if(DoTest)write(*,*)'couple_sc_ih passed set_router iProc=',i_proc()
       SC_iGridInScIh= SC_iGridInIhSc
    end if
    call couple_buffer_grid(&
         RouterScBuff,&
         nVar=8,&
         fill_buffer=SC_get_for_ih,&
         NameBuffer='IH_from_sc',&
         TargetID_=IH_)
!    if(DoTest)write(*,*)'couple_sc_ih: passed messagepass iProc=',&
!         i_proc()
    if(DoTest.and.is_proc0(compid_grid(BuffGD%DD%Ptr)))then
       nU_I=ubound_vector('IH_from_sc')
       iCoupling=iCoupling+1
       iFile=io_unit_new()
       write(NameFile,'(a,i4.4,a)')'./IH/from_sc_',iCoupling,'.dat'
       open(iFile,FILE=NameFile,STATUS='unknown')
       do iPoint=1,nU_I(2)
          write(iFile,*)point_state_v('IH_from_sc',8,iPoint)
       end do
       close(iFile)
    end if
    if(DoTest)write(*,*)'Couple passed at PE=',i_proc()
  end subroutine couple_sc_ih
  !======================================================!
  subroutine buffer_grid_point(&
       nDimFrom,Sph_D,nDimTo,Xyz_D,IsInterfacePoint)         
    integer,intent(in)::nDimFrom,nDimTo       
    real,dimension(nDimFrom),intent(in)::Sph_D
    real,dimension(nDimTo),intent(out)::Xyz_D
    logical,intent(out)::IsInterfacePoint
    !The order of spherical indexes as in BATSRUS
    integer,parameter::R_=1,Psi_=2,Theta_=3,x_=1,y_=2,z_=3
    real::RSinTheta
    IsInterfacePoint=.true.
    !Transform to cartesian 
    RSinTheta=Sph_D(R_)*sin(Sph_D(Theta_))!To be modified
    Xyz_D(x_)=RSinTheta*cos(Sph_D(Psi_))  !\because SC grid
    Xyz_D(y_)=RSinTheta*sin(Sph_D(Psi_))  !/may be spherical
    Xyz_D(z_)=Sph_D(R_)*cos(Sph_D(Theta_))!To be modified
    !Xyz_D should be transformed to Carrington coordinates,
    !which have not been implemented yet. 
  end subroutine buffer_grid_point
end module CON_couple_ih_sc

