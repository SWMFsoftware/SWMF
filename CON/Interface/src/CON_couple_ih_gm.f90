!^CMP COPYRIGHT UM
!^CMP FILE IH
!^CMP FILE GM
!BOP
!MODULE: CON_couple_ih_gm - couple IH to GM one way
!INTERFACE:
module CON_couple_ih_gm
!DESCRIPTION:
!This coupler is  based on the global message  pass.
!The toolkit for coupling the different data sets within a
!single component code, or within a framework as it is now 
!
!The features of the present coupler: the domain decompositions
!both for a target and for a source are with AMR.
!IH is a source, GM is a target

!Typically the mapping is charcterized by some real matrices
!and/or vectors. In the particular case of the cartesian
!generalized coordinates this is the position of the origin of
!the GM frame of reference in the IH frame of refernce  and the
!matrix which characterizes the turn of the IH frame of
!reference with respect to GM :
!
!$$
!{\bf V}_{IH}={\bf A}_{IH,GM}\cdot{\bf V}_{GM},\\
!{\bf V}_{GM}={\bf A}_{GM,IH}\cdot{\bf V}_{IH}
!$$
  !USES:
  use CON_coupler
  use CON_time,ONLY:TimeStart
  use ModConst
  use CON_geopack
  use CON_axes, ONLY: transform_matrix, vPlanetHgi_D

  implicit none
  private !except

  !
  !PUBLIC MEMBER FUNCTIONS:
  public:: couple_ih_gm_init
  public:: couple_ih_gm

  real,public::IH_GM_CouplingTime

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V. <igorsok@umich.edu> - initial prototype
  ! 9/02/03 G.Toth <gtoth@umich.edu> - minor changes
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: IH_iGridRealization=-2
  integer :: GM_iGridRealization=-2
  integer :: IH_GM_iMapping=-1

  real,dimension(3)   :: IH_EarthPosition_D,vPlanetIh_D
  real,dimension(3,3) :: IH_a_GM_DD,GM_a_IH_DD

  type(RouterType),save         :: Router
  type(GridDescriptorType),save :: GM_Grid
  type(GridDescriptorType),save :: IH_Grid

  logical :: DoInitialize=.true., DoTest, DoTestMe

  character(len=*), parameter :: NameMod='CON_couple_ih_gm'

  ! The integer time array for GeoPack, 
  ! only first 6 elements are used by GeoPack
  integer :: TimeGeoPack(7)

contains

  !===============================================================!
  subroutine couple_ih_gm_init

    if(.not.DoInitialize)return
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm_init iProc=',i_proc()

    !Initialize the transformation matrices, compute the position 
    !of the Earth and its orbital velocity in IH, at the time zero  

    IH_a_GM_DD         = transform_matrix(cZero, &
         Grid_C(GM_) % TypeCoord, Grid_C(IH_) % TypeCoord)

    GM_a_IH_DD         = transpose(IH_a_GM_DD)

    !Compute the Earth position in Hgi and transform 
    !to IH coordinates
 
    IH_EarthPosition_D = -cAU/rSun*SunEMBDistance*matmul(&
         transform_matrix(cZero, 'HGI', &
         Grid_C(IH_) % TypeCoord),HgiGse_DD(:,1))

    !Compute the Earth velocity in Hgi and transform 
    !to IH coordinates 

    vPlanetIh_D = matmul(transform_matrix(cZero,'HGI',&
            Grid_C(IH_) % TypeCoord),vPlanetHgi_D)

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
!IROUTINE: couple_ih_gm - interpolate and get  MHD state at GM west ghostpoints
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
    integer::iMappingNew

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call GM_synchronize_refinement(Router%iProc0Target,Router%iCommUnion)

    IH_GM_CouplingTime=TimeCoupling

    call update_mapping_parameters(iMappingNew)

    if(IH_iGridRealization/=i_realization(IH_).or.&     
         GM_iGridRealization/=i_realization(GM_).or.&
         IH_GM_iMapping/=iMappingNew)then  

       if(DoTest)write(*,*)'couple_ih_gm call set_router iProc=',i_proc()
       
       if(DoTestMe)write(*,*)'couple_ih_gm_init IH_EarthPosition_D=',&
            IH_EarthPosition_D
       call set_router(&
            GridDescriptorSource=IH_Grid,&
            GridDescriptorTarget=GM_Grid,&
            Router=Router,&
            is_interface_block=GM_west_block,&
            interface_point_coords=GM_west_cells, &
            mapping=GM_IH_mapping,&
            interpolate=interpolation_fix_reschange)
       IH_iGridRealization=i_realization(IH_)
       GM_iGridRealization=i_realization(GM_)
       IH_GM_iMapping=iMappingNew      
    end if

    call couple_comp(&
         Router,&
         nVar=8,&
         fill_buffer=IH_get_for_gm_and_transform,&
         apply_buffer=GM_put_from_ih)

  end subroutine couple_ih_gm
  !===============================================================!
  logical function GM_west_block(lGlobalTreeNode)

    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsRightBoundary_D
    integer,parameter::x_=1

    IsRightBoundary_D=is_right_boundary_d(&
         GM_Grid%DD%Ptr,lGlobalTreeNode)
    GM_west_block=IsRightBoundary_D(x_)

  end function GM_west_block
  !===============================================================!
  subroutine GM_west_cells(&
       GridDescriptor,&
       lGlobalTreeNode,&
       nDim,&
       Xyz_D,&
       nIndexes,&
       Index_I,&
       IsInterfacePoint)

    type(GridDescriptorType),intent(in):: GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    integer,intent(in)::nDim
    real,intent(inout)::Xyz_D(nDim)
    integer,intent(inout)::Index_I(nIndexes)
    logical,intent(out)::IsInterfacePoint

    logical,dimension(3)::IsLeftFace_D,IsRightFace_D
    integer,parameter::x_=1,y_=2,z_=3

    IsLeftFace_D=Index_I(x_:z_)<1
    IsRightFace_D=Index_I(x_:z_)>&
         ncells_decomposition_d(GridDescriptor%DD%Ptr)
    IsInterfacePoint=IsRightFace_D(x_).and..not.&
         (any(IsLeftFace_D(y_:z_)).or.any(IsRightFace_D(y_:z_)))

  end subroutine GM_west_cells

  !===============================================================!

  subroutine GM_IH_mapping(&
       GM_nDim,GM_Xyz_D,IH_nDim,IH_Xyz_D,IsInterfacePoint)

    integer,intent(in)::GM_nDim,IH_nDim
    real,dimension(GM_nDim),intent(in)::GM_Xyz_D
    real,dimension(IH_nDim),intent(out)::IH_Xyz_D
    logical,intent(out)::IsInterfacePoint

    IH_Xyz_D=IH_EarthPosition_D+&
         matmul(IH_a_GM_DD,GM_Xyz_D)*REarth/RSun
    IsInterfacePoint=.true.

  end subroutine GM_IH_Mapping

  !===============================================================!

  subroutine update_mapping_parameters(iMappingNew)

    use ModKind
    use CON_physics, ONLY: get_time, time_real_to_int
    use CON_geopack
    ! The time at the beginning of the simulation
    real(Real8_), save :: tStart
    integer,intent(out)::iMappingNew

    ! For sake of simplicity the initial conversion matrix is used.
    iMappingNew=1
    RETURN

    !Initialize the transformation matrices, compute the position 
    !of the Earth and its orbital velocity in IH, at the coupling time  

    !This matrix is needed to map the points of the target (GM)
    !to the source grid (IH), for interpolation
    IH_a_GM_DD         = transform_matrix(IH_GM_CouplingTime, &
            Grid_C(GM_) % TypeCoord, Grid_C(IH_) % TypeCoord)

    !The transposed matrix is needed to convert vectors to be sent
    !from the source (IH) to target (GM)

    GM_a_IH_DD         = transpose(IH_a_GM_DD)

    !Compute the Earth position in Hgi and transform 
    !to IH coordinates

    IH_EarthPosition_D = -cAU/rSun*SunEMBDistance*&
         matmul(transform_matrix(IH_GM_CouplingTime,'HGI',&
         Grid_C(IH_) % TypeCoord),HgiGse_DD(:,1))

    !Compute the Earth velocity in Hgi and transform 
    !to IH coordinates

    vPlanetIh_D = matmul(transform_matrix(IH_GM_CouplingTime,'HGI',&
         Grid_C(IH_) % TypeCoord),vPlanetHgi_D)

  end subroutine update_mapping_parameters

  !===============================================================!

  subroutine IH_get_for_gm_and_transform(&
       nPartial,iGetStart,Get,W,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    integer, parameter :: rho_=1, rhoUx_=2, rhoUz_=4, Bx_=5, Bz_=7
    !------------------------------------------------------------
    call IH_get_for_gm(&
       nPartial,iGetStart,Get,W,State_V,nVar,IH_GM_CouplingTime)

    State_V(rhoUx_:rhoUz_)=&
         matmul(GM_a_IH_DD,State_V(rhoUx_:rhoUz_) &
         - State_V(rho_)*vPlanetIh_D )
    State_V(Bx_:Bz_)=matmul(GM_a_IH_DD,State_V(Bx_:Bz_))
  end subroutine IH_get_for_gm_and_transform

end module CON_couple_ih_gm

