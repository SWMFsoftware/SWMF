!^CMP COPYRIGHT UM
!^CMP FILE IH
!^CMP FILE GM
!BOP
!MODULE: CON_couple_ih_gm - couple IH to GM one way
!INTERFACE:
module CON_couple_ih_gm
!DESCRIPTION:
! This coupler uses the SWMF coupling toolkit.
! Both the IH and GM grids use AMR.
! IH is a source, GM is a target.
! It is assumed that IH uses the solar radius rSun as distance unit, while
! GM uses the Earth's radius as distance unit.

! Typically the mapping is charcterized by some real matrices
! and/or vectors. In the particular case of the cartesian
! generalized coordinates this is the position of the origin of
! the GM frame of reference in the IH frame of refernce  and the
! matrix which characterizes the turn of the IH frame of
! reference with respect to GM :
!
! $$
! {\bf V}_{IH}={\bf A}_{IH,GM}\cdot{\bf V}_{GM},\\
! {\bf V}_{GM}={\bf A}_{GM,IH}\cdot{\bf V}_{IH}
! $$
  !USES:
  use CON_coupler
  use CON_time,ONLY:TimeStart
  use ModConst
  use CON_axes, ONLY: transform_matrix, vPlanetHgi_D, XyzPlanetHgi_D

  implicit none

  save

  private !except

  !
  !PUBLIC MEMBER FUNCTIONS:
  public:: couple_ih_gm_init
  public:: couple_ih_gm

  real, public:: CouplingTimeIhGm

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V. <igorsok@umich.edu> - initial prototype
  ! 9/02/03 G.Toth <gtoth@umich.edu> - minor changes
  ! 6/18/05 Sokolov - generalized to arbitrary coordinate systems
  ! 6/22/05 G.Toth  - redo mapping and router to allow for relative motion
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: IH_iGridRealization=-2
  integer :: GM_iGridRealization=-2

  real, dimension(3)   :: XyzPlanetIh_D, vPlanetIh_D
  real, dimension(3,3) :: GmToIh_DD, IhToGm_DD, HgiToIh_DD

  ! Maximum time difference [s] without remap
  ! The 600 s corresponds to about 0.1 degree rotation between IH and GM
  real :: dTimeMappingMax = 600.0

  type(RouterType)         :: Router
  type(GridDescriptorType) :: GM_Grid
  type(GridDescriptorType) :: IH_Grid

  logical :: DoInitialize=.true., DoTest, DoTestMe

  character(len=*), parameter :: NameMod='CON_couple_ih_gm'

contains

  !===============================================================!
  subroutine couple_ih_gm_init

    if(.not.DoInitialize)return
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm_init iProc=',i_proc()

    call init_coupler(              &    
       iCompSource=IH_,             & ! component index for source
       iCompTarget=GM_,             & ! component index for target
       nGhostPointTarget=2,         & ! number of halo points in target
       GridDescriptorSource=IH_Grid,& ! OUT!\
       GridDescriptorTarget=GM_Grid,& ! OUT!-General coupler variables 
       Router=Router)                 ! OUT!/

    IH_iGridRealization=-1
    GM_iGridRealization=-1

    ! Initialize the coordinate transformation

  end subroutine couple_ih_gm_init
  !===============================================================!
  !BOP
  !IROUTINE: couple_ih_gm - get IH solution in the GM inflow boundary ghost cells
  !INTERFACE:
  subroutine couple_ih_gm(TimeCoupling)
    !INPUT ARGUMENTS:
    interface
       subroutine GM_put_from_ih(nPartial,&
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
       end subroutine GM_put_from_ih
    end interface
    real,intent(in)::TimeCoupling
    !EOP

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call GM_synchronize_refinement(Router%iProc0Target,Router%iCommUnion)

    CouplingTimeIhGm=TimeCoupling

    if(IH_iGridRealization/=i_realization(IH_).or.&     
         GM_iGridRealization/=i_realization(GM_).or.&
         TimeCoupling - TimeCouplingLast > dTimeMappingMax)then  

       ! Set the transformation matrices, compute the position and
       ! orbital veloctiy of the Earth in IH at the coupling time.

       GmToIh_DD  = transform_matrix(TimeCoupling , &
            Grid_C(GM_) % TypeCoord, Grid_C(IH_) % TypeCoord)
       IhToGm_DD  = transpose(GmToIh_DD)

       HgiToIh_DD = transform_matrix(TimeCoupling, &
            'HGI', Grid_C(IH_) % TypeCoord)

       ! Transform the Earth position in HGI to IH coordinates 
       ! and change distance units to that used in IH
 
       XyzPlanetIh_D = matmul(HgiToIh_DD, XyzPlanetHgi_D)/Grid_C(IH_)%UnitX

       ! Transform the Earth velocity from Hgi to IH coordinates 

       vPlanetIh_D = matmul(HgiToIh_DD, vPlanetHgi_D)

       if(DoTest)write(*,*)'couple_ih_gm call set_router iProc=',i_proc()
       
       if(DoTestMe)write(*,*)'couple_ih_gm_init XyzPlanetIh_D=',&
            XyzPlanetIh_D

       call set_router(&
            GridDescriptorSource=IH_Grid,&
            GridDescriptorTarget=GM_Grid,&
            Router=Router,&
            is_interface_block=GM_is_west_block,&
            interface_point_coords=GM_west_cells, &
            mapping=map_gm_ih,&
            interpolate=interpolation_fix_reschange)

       IH_iGridRealization = i_realization(IH_)
       GM_iGridRealization = i_realization(GM_)
       TimeCouplingLast    = TimeCoupling
    end if

    call couple_comp(&
         Router,&
         nVar=8,&
         fill_buffer=IH_get_for_gm_and_transform,&
         apply_buffer=GM_put_from_ih)

  end subroutine couple_ih_gm
  !===============================================================!
  logical function GM_is_west_block(lGlobalTreeNode)

    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsRightBoundary_D
    integer,parameter::x_=1

    IsRightBoundary_D=is_right_boundary_d(&
         GM_Grid%DD%Ptr,lGlobalTreeNode)
    GM_is_west_block=IsRightBoundary_D(x_)

  end function GM_is_west_block
  !===============================================================!
  subroutine GM_west_cells(&
       GridDescriptor,&
       lGlobalTreeNode,&
       nDim,&
       Xyz_D,&
       nIndexes,&
       i_D,&
       IsInterfacePoint)

    type(GridDescriptorType),intent(in):: GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    integer,intent(in)::nDim
    real,intent(inout)::Xyz_D(nDim)
    integer,intent(inout)::i_D(nIndexes)
    logical,intent(out)::IsInterfacePoint

    logical,dimension(3)::IsLeftFace_D,IsRightFace_D
    integer,parameter::x_=1,y_=2,z_=3

    IsLeftFace_D=i_D(x_:z_)<1
    IsRightFace_D=i_D(x_:z_)>&
         ncells_decomposition_d(GridDescriptor%DD%Ptr)
    IsInterfacePoint=IsRightFace_D(x_).and..not.&
         (any(IsLeftFace_D(y_:z_)).or.any(IsRightFace_D(y_:z_)))

  end subroutine GM_west_cells

  !===============================================================!

  subroutine map_gm_ih(&
       GM_nDim,GM_Xyz_D,IH_nDim,IH_Xyz_D,IsInterfacePoint)

    integer,intent(in)::GM_nDim,IH_nDim
    real,dimension(GM_nDim),intent(in)::GM_Xyz_D
    real,dimension(IH_nDim),intent(out)::IH_Xyz_D
    logical,intent(out)::IsInterfacePoint
    !In each mapping the corrdinates of the TARGET grid point (GM)
    !shoud be be transformed to the SOURCE (IH) generalized coords.

    IH_Xyz_D = XyzPlanetIh_D + matmul(GmToIh_DD, GM_Xyz_D)*&
         Grid_C(GM_)%UnitX/Grid_C(IH_)%UnitX
    IsInterfacePoint=.true.

  end subroutine map_gm_ih

  !===============================================================!

  subroutine IH_get_for_gm_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V

    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7
    !------------------------------------------------------------
    call IH_get_for_gm(&
       nPartial,iGetStart,Get,w,State_V,nVar,CouplingTimeIhGm)

    State_V(RhoUx_:RhoUz_)=&
         matmul(IhToGm_DD,State_V(RhoUx_:RhoUz_) &
         - State_V(Rho_)*vPlanetIh_D )
    State_V(Bx_:Bz_)=matmul(IhToGm_DD,State_V(Bx_:Bz_))

  end subroutine IH_get_for_gm_and_transform

end module CON_couple_ih_gm

