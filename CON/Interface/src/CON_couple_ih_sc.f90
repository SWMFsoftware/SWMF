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
  public:: couple_ih_sc

  real,public::IH_SC_CouplingTime

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V.<igorsok@umich.edu> - prototype for ih-gm
  ! 7/04/04                                 - version for ih-sc
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: IH_iGridRealization=-2
  integer :: SC_iGridRealization=-2
  
 
  type(RouterType),save         :: Router
  type(GridDescriptorType),save :: SC_Grid
  type(GridDescriptorType),save :: IH_Grid

  logical :: DoInitialize=.true., DoTest, DoTestMe

  character(len=*), parameter :: NameMod='CON_couple_sc_gm'


contains

  !===============================================================!
  subroutine couple_ih_sc_init

    if(.not.DoInitialize)return
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm_init iProc=',i_proc()

    call init_coupler(              &    
       iCompSource=IH_,             & ! component index for source
       iCompTarget=SC_,             & ! component index for target
       nGhostPointTarget=2,         & ! number of halo points in target
       GridDescriptorSource=IH_Grid,& ! OUT!\
       GridDescriptorTarget=SC_Grid,& ! OUT!-General coupler variables 
       Router=Router)                 ! OUT!/

    IH_iGridRealization=-1
    SC_iGridRealization=-1


  end subroutine couple_ih_sc_init
  !===============================================================!
!BOP
!IROUTINE: couple_ih_gm - interpolate and get  MHD state at GM west ghostpoints
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

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call SC_synchronize_refinement(Router%iProc0Target,Router%iCommUnion)

    IH_SC_CouplingTime=TimeCoupling


    if(IH_iGridRealization/=i_realization(IH_).or.&     
         SC_iGridRealization/=i_realization(SC_))then  

       if(DoTest)write(*,*)'couple_ih_gm call set_router iProc=',i_proc()

       call set_router(&
            GridDescriptorSource=IH_Grid,&
            GridDescriptorTarget=SC_Grid,&
            Router=Router,&
            is_interface_block=boundary_block,&
            interface_point_coords=outer_cells, &
            interpolate=interpolation_fix_reschange)
       IH_iGridRealization=i_realization(IH_)
       SC_iGridRealization=i_realization(SC_)
    end if

    call couple_comp(&
         Router,&
         nVar=8,&
         fill_buffer=IH_get_for_sc,&
         apply_buffer=SC_put_from_ih)

  end subroutine couple_ih_sc
  !===============================================================!
  logical function boundary_block(lGlobalTreeNode)

    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsBoundary_D

    IsBoundary_D=is_right_boundary_d(&
         SC_Grid%DD%Ptr,lGlobalTreeNode).or.&
         is_left_boundary_d(&
         SC_Grid%DD%Ptr,lGlobalTreeNode)
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
         SC_Grid%DD%Ptr,lGlobalTreeNode)
    IsRightFace_D=Index_I(x_:z_)>&
         ncells_decomposition_d(SC_Grid%DD%Ptr).and.&
         is_right_boundary_d(SC_Grid%DD%Ptr,lGlobalTreeNode)
    IsInterfacePoint=any(IsRightFace_D.or.IsLeftFace_D)
  end subroutine outer_cells 
end module CON_couple_ih_sc

