!^CMP COPYRIGHT UM
!^CMP FILE GM
!^CMP FILE IE
!BOP
!MODULE:CON_couple_gm_ie_swmf - two-way GM_IE coupling using SWMF 
!INTERFACE:
module CON_couple_gm_ie_swmf
!DESCRIPTION:
!Couple between two components:
!    Global Magnetosphere       (GM) source
!    Ionosphere Electrodynamics (IE) target
!Coupler features: 
!    At least one of the grids (GM) has an AMR
!    Router is constructed from source (GM)
!EOP
!BOP
!USES:
  use CON_coupler

  implicit none

  private !Except

!PUBLIC MEMBER FUNCTIONS:
  public::init_couple_gm_ie_swmf
  public::couple_gm_ie_swmf
  public::couple_ie_gm_swmf
!EOP
  !\
  ! General coupling variables
  !/

  ! Name of this interface
  character (len=*), parameter :: NameSub='couple_gm_ie'
  type(GridDescriptorType),save:: GridGm
  type(GridDescriptorType),save:: GridIe
  type(GridDescriptorType),save:: GridGmIe
  type(RouterType),save::RouterGmIe
  type(RouterType),save::RouterIeGm

  integer,parameter::GmIeGrid_=MaxComp+1

  real::rCurrent2,rIonosphere
  real::AxisMag_D(3)=(/cZero,cZero,cOne/),GmIe_DD(3,3)
  ! Size of the 2D spherical structured 
  !(possibly non-uniform) IE grid
  integer:: nPoint_B(2)
  !Limit value for colatitude
  real::ColatLim_B(2)

  integer::iGridRealizationGm=-2
  integer::iMappingGmIe=-2
  integer::iMappingIeGm=-2
  integer::iMappingNew=-2
  
  logical :: DoTest, DoTestMe
  integer :: iProcWorld
  logical :: DoInit=.true.

  logical:: IsNewFacPoint
  real::TimeCoupling=cZero
contains
!BOP
!IROUTINE: init_couple_gm_ie -  initialize coupling
!INTERFACE:
  subroutine init_couple_gm_ie_swmf
    !INPUT ARGUMENTS:
    interface
       subroutine GMIE_set_grid
         implicit none
       end subroutine GMIE_set_grid
       subroutine GM_get_mapping_param_for_ie(& 
            rCurrentGm,rIonosphere)    !srcGM/map_ionosphere_bc.f90
         implicit none
         real,intent(out)::rCurrentGm,rIonosphere
       end subroutine GM_get_mapping_param_for_ie
    end interface
    integer::iError
    real::rCurrentsGm

    !EOP
    if(.not.DoInit)return
    DoInit=.false.
    call init_coupler(                          &
         iCompSource=GM_,                       &            
         iCompTarget=IE_,                       &  
         GridDescriptorSource=GridGm,           &
         GridDescriptorTarget=GridIe,           &
         nGhostPointTarget=1,                   &
         StandardTarget_=Nodes_,                &
         Router=RouterGmIe)
    if(.not.done_dd_init(GmIeGrid_))call GMIE_set_grid
    call set_standard_grid_descriptor(&
         GridID_=GmIeGrid_,&
         nGhostGridPoints=1,&
         Standard_=Nodes_,&
         GridDescriptor=GridGmIe)
    call init_router(&
         GridDescriptorSource=GridIe,&
         GridDescriptorTarget=GridGmIe,&
         Router=RouterIeGm)

    if(.not.RouterGmIe%IsProc)return
    iGridRealizationGm=-1
    iMappingGmIe=-1
    iMappingIeGm=-1
    if(is_proc0(GM_))then
       call GM_get_mapping_param_for_ie(rCurrentsGm,rIonosphere)
       rCurrent2=rCurrentsGm**2
    end if
    call MPI_BCAST(rCurrent2,1,MPI_REAL,&
         RouterGmIe%iProc0Source,RouterGmIe%iCommUnion,iError)
    call MPI_BCAST(rIonosphere,1,MPI_REAL,&
         RouterGmIe%iProc0Source,RouterGmIe%iCommUnion,iError)
    
  end subroutine init_couple_gm_ie_swmf
!========================COUPLER VIA SWMF=================!
!BOP
!IROUTINE: couple_gm_ie - coupling GM and IE via a router.
!INTERFACE: 
  subroutine couple_gm_ie_swmf(tSimulation)
!INPUT ARGUMENTS:
    interface
       subroutine GM_get_for_ie_swmf(&      !srcGM/GM_get_for_ie_swmf.f90
            nPartial,iGet,Get,W,State_V,nVar)
         use CON_router
         implicit none
         !INPUT ARGUMENTS  
         integer,intent(in)::nPartial,iGet,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::W
         !OUTPUT ARGUMENTS
         real,dimension(nVar),intent(out)::State_V
       end subroutine GM_get_for_ie_swmf

       subroutine IE_check_allocation_south(nPoint)!srcIE/IE_wrapper.f90
         implicit none
         integer,intent(in)::nPoint
       end subroutine IE_check_allocation_south

       subroutine IE_check_allocation_north(nPoint)!srcIE/IE_wrapper.f90
         implicit none
         integer,intent(in)::nPoint
       end subroutine IE_check_allocation_north

       subroutine GM_synchronize_refinement(&     !srcGM/GM_wrapper.f90
            iProc0,iCommUnion)
         implicit none
         integer,intent(in)::iProc0,iCommUnion
       end subroutine GM_synchronize_refinement

       subroutine IE_interpolate(iBlock,ColatLim) !srcIE/IE_interpolate.f90
         implicit none
         integer,intent(in)::iBlock
         real,intent(in)::ColatLim
       end subroutine IE_interpolate
    
       subroutine ionosphere_fac(iBlock) !srcIE/iono_coupling.f90
         implicit none
         integer,intent(in)::iBlock
       end subroutine ionosphere_fac
    end interface
    real,intent(in)::tSimulation
!REVISION HISTORY:
!15AUG03     - Toth&Ridley,initial prototype for coupler via MPI 
!25AUG03     - I.Sokolov<igorsok@umich.edu> - code/prolog
!EOP

    ! Logical for new mapping points
    logical :: IsNewFacPoint

    integer::iPE,iPoint,iBlock,iMappingNew,iError
    !--------------------------------------------------------------

! After everything is initialized exclude PEs which are not involved
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(iGridRealizationGm==-2)call init_couple_gm_ie_swmf
    if(.not.RouterGmIe%IsProc)return
    TimeCoupling=tSimulation
    call MPI_BCAST(TimeCoupling,1,MPI_REAL,&
         RouterGmIe%iProc0Source, RouterGmIe%iCommUnion,iError)
    call update_mapping(iMappingNew)
    call GM_synchronize_refinement(&
         RouterGmIe%iProc0Source,&    !The GlobalGrid is a source
         RouterGmIe%iCommUnion)

    IsNewFacPoint=iGridRealizationGm/=i_realization(GM_).or.&
         iMappingGmIe/=iMappingNew
    if(DoTest)write(*,*)NameSub//': iMapping,iGrid,IsNewFacPoint',&
         iMappingNew,i_realization(GM_),IsNewFacPoint
    if(IsNewFacPoint)then
!BOP:
!DESCRIPTION:
!Constructor fot the router
!\begin{verbatim}
       call construct_router_from_source(&
            GridGm,&
            GridIe,&
            RouterGmIe,&
            is_interface_block=is_block_near_r_current,&
            interface_point_coords=near_r_current,&
            mapping=sort_hemispheres)
!\end{verbatim}
!EOP
       !Allocate arrays in IE_ to get data
       nPoint_B=0
       do iPE=0,RouterGmIe%nProc-1
          if(RouterGmIe%nPut_P(iPE)==0)CYCLE
          do iPoint=1,RouterGmIe%nPut_P(iPE)
             nPoint_B(RouterGmIe%iPut_P(iPE)%iCB_II(3,iPoint))=&
                  nPoint_B(RouterGmIe%iPut_P(iPE)%iCB_II(3,iPoint))+1
          end do
       end do
       if(nPoint_B(1)>0)then
          call IE_check_allocation_north(nPoint_B(1))
          write(*,*) &
             'IE: Northern Magnetospheric FAC Mapping Points = ', &
             nPoint_B(1)
       end if
       if(nPoint_B(2)>0)then
          call IE_check_allocation_south(nPoint_B(2))
          write(*,*) &
               'IE: Southern Magnetospheric FAC Mapping Points = ', &
               nPoint_B(2)
          write(*,*)'IE:'
       end if
       !Done allocation

       !Coupling: send-receive FAC(3 values), 
       !          magnetic field data (5 values),
       !          location in the ionosphere(3 values)
       nPoint_B=0
       ColatLim_B(1)=cZero; ColatLim_B(2)=cPi
       call couple_comp(RouterGmIe,&
            11,&
            GM_get_for_ie_swmf,&
            interface_to_ie_put_from_gm)
       iGridRealizationGm=i_realization(GM_)
       iMappingGmIe=iMappingNew
       if(is_proc(IE_))then
          do iBlock=1,2
             if(nPoint_B(iBlock)>0)then
                call IE_interpolate(iBlock,ColatLim_B(iBlock))
                if(nPoint_B(iBlock)>8)call ionosphere_fac(iBlock)
             end if
          end do
          if(DoTest)write(*,*)NameSub//':',i_proc(),': ColatLimB=',&
               ColatLim_B,' nPoint_B:',nPoint_B
       end if
    else
       nPoint_B=0
       ColatLim_B(1)=cZero; ColatLim_B(2)=cPi
       call couple_comp(RouterGmIe,&
            3,&
            GM_get_for_ie_swmf,&
            interface_to_ie_put_from_gm)
       if(is_proc(IE_))then
          do iBlock=1,2
             if(nPoint_B(iBlock)>8)call ionosphere_fac(iBlock)
          end do
       if(DoTest)write(*,*)NameSub//':',i_proc(),': ColatLimB=',&
            ColatLim_B,' nPoint_B:',nPoint_B
       end if
    end if
  end subroutine couple_gm_ie_swmf
  !
!BOP
!INTERFACE:
  subroutine interface_to_ie_put_from_gm(&
       nPartial,iPutStart,Put,W,DoAdd,State_V,nVar)
!INPUT ARGUMENTS:
    interface
       subroutine IE_put_from_gm_swmf(&       !srcIE/IE_wrapper.f90
            State_V,nVar,iBlock,nPoint,ColatLim)
         implicit none
         !INPUT ARGUMENTS
         integer,intent(in)::nVar,iBlock,nPoint
         real,dimension(nVar),intent(in):: State_V
         real,intent(inout)::ColatLim
       end subroutine IE_put_from_gm_swmf
    end interface

    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd !Is not used
    real,dimension(nVar),intent(in)::State_V
!EOP
    integer::iBlock
    iBlock=Put%iCB_II(3,iPutStart)
    nPoint_B(iBlock)=nPoint_B(iBlock)+1
    call IE_put_from_gm_swmf(&
         State_V,&
         nVar,&
         iBlock,&
         nPoint_B(iBlock),&
         ColatLim_B(iBlock))
  end subroutine interface_to_ie_put_from_gm
!BOP
!IROUTINE: near_r_current - GM_ points near r=rCurrent surface to get FAC 
!INTERFACE: 
  subroutine near_r_current(&
       GridDescriptor,&
       lGlobalTreeNode,&
       Xyz_D,&
       nIndexes,&
       Index_I,&
       IsInterfacePoint)
!INPUT ARGUMENTS:
    type(GridDescriptorType),intent(in)::GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    real,dimension(GridDescriptor%nDim),intent(inout)::Xyz_D
    integer, dimension(nIndexes),intent(inout)::Index_I
!OUTPUT ARGUMENTS:
    logical,intent(out)::IsInterfacePoint
!EOP
    real,dimension(3)::XyzNei_D,DXyz_D

    IsInterfacePoint=.false.
    if(is_inside_r_current(Xyz_D))return
    DXyz_D=d_xyz_cell_d(GM_,lGlobalTreeNode)
    XyzNei_D(2:3)=Xyz_D(2:3);XyzNei_D(1)=Xyz_D(1)+DXyz_D(1)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
    XyzNei_D(2:3)=Xyz_D(2:3);XyzNei_D(1)=Xyz_D(1)-DXyz_D(1)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
    XyzNei_D=Xyz_D;XyzNei_D(2)=Xyz_D(2)+DXyz_D(2)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
    XyzNei_D=Xyz_D;XyzNei_D(2)=Xyz_D(2)-DXyz_D(2)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
    XyzNei_D=Xyz_D;XyzNei_D(3)=Xyz_D(3)+DXyz_D(3)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
    XyzNei_D=Xyz_D;XyzNei_D(3)=Xyz_D(3)-DXyz_D(3)
    if(is_inside_r_current(XyzNei_D))then
       IsInterfacePoint=.true.
       return
    end if
  end subroutine near_r_current
  !================================================================!
!BOP
!IROUTINE: is_block_near_r_current - block involving GM-points to get FAC 
!INTERFACE: 
  logical function is_block_near_r_current(iBlock)
!INPUT PARAMETERS:
    integer,intent(in)::iBlock
!EOP

    real,dimension(3)::Xyz_D,DXyz_D
    integer::i,j,k
    integer,dimension(3),parameter::i_D=(/1,0,0/)
    integer,dimension(3),parameter::j_D=(/0,1,0/)
    integer,dimension(3),parameter::k_D=(/0,0,1/)
    Xyz_D=xyz_block_d(GM_,iBlock)
    DXyz_D=d_xyz_block_d(GM_,iBlock)
    is_block_near_r_current=all(Xyz_D*(Xyz_D+DXyz_D)<=cZero)
    if(is_block_near_r_current)return
    do k=0,1; do j=0,1; do i=0,1
       is_block_near_r_current=&
            is_inside_r_current(Xyz_D+DXyz_D*(i*i_D+j*j_D+k*k_D))
       if(is_block_near_r_current)return
    end do;end do; end do
    is_block_near_r_current=.false.
  end function is_block_near_r_current

  logical function is_inside_r_current(Xyz_D)
    real,dimension(3),intent(in)::Xyz_D
    is_inside_r_current=dot_product(Xyz_D,Xyz_D)<rCurrent2
  end function is_inside_r_current
  !================================================================!
!BOP
!IROUTINE:sort_hemispherees - mapping for the router
!INTERFACE:
  subroutine sort_hemispheres(&
       nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
!INPUT ARGUMENTS:
    integer,intent(in)::nDimFrom,nDimTo       
    real,dimension(nDimFrom),intent(in)::XyzFrom_D
!OUTPUT ARGUMENTS:
    real,dimension(nDimTo),intent(out)::XyzTo_D
    logical,intent(out)::IsInterfacePoint
!DESCRIPTION:
!The mapping used to construct the router is very simple in this
!case, because one should only find out, into which of the two
!hemispheres (North, South) the point XyzFrom\_D is projected along 
!the field line. Exact position of the point XyzTo\_D is not needed 
!EOP
    IsInterfacePoint=.true.
    if(dot_product(XyzFrom_D,AxisMag_D)<0)then
       XyzTo_D=xyz_grid_d(GridIe,2,(/1,1/))
    else
       XyzTo_D=xyz_grid_d(GridIe,1,(/1,1/))
    end if
  end subroutine sort_hemispheres
  !================================================================!
  !==================IE TO GM COUPLER==============================!
!BOP
!IROUTINE: couple_gm_ie - coupling GM and IE via a router.
!INTERFACE: 
  subroutine couple_ie_gm_swmf(tSimulation)
!INPUT ARGUMENTS:
    interface
       subroutine IE_get_for_gm_swmf(nPartial,iGetStart,Get,W,State_V,nVar)
         use CON_router                          !See srcIE/IE_wrapper
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::W
         real,dimension(nVar),intent(out)::State_V
       end subroutine IE_get_for_gm_swmf
       subroutine GM_put_from_ie_swmf(nPartial,iPut,Put,W,DoAdd,State_V,nVar)
         use CON_router                           !See srcGM
         implicit none
         integer,intent(in)::nPartial,iPut,nVar
         type(IndexPtrType),intent(in)::Put
         type(WeightPtrType),intent(in)::W
         logical,intent(in)::DoAdd
         real,dimension(nVar),intent(in)::State_V
       end subroutine GM_put_from_ie_swmf
       subroutine transform_phi_bc_to_u_bc
         implicit none
       end subroutine transform_phi_bc_to_u_bc
       subroutine IE_run(tStart,tEnd)
         implicit none
         real::tStart,tEnd
       end subroutine IE_run
    end interface
    real,intent(in)::tSimulation
!REVISION HISTORY:
!15AUG03     - Toth&Ridley,initial prototype for coupler via MPI 
!05SEP03     - I.Sokolov<igorsok@umich.edu> - code/prolog
!EOP
    integer::iError
    real    :: tSimulationTmp
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(.not.RouterIeGm%IsProc)return
    TimeCoupling=tSimulation
    call MPI_BCAST(TimeCoupling,1,MPI_REAL,&
         RouterIeGm%iProc0Target, RouterIeGm%iCommUnion,iError)
  ! Make sure that the most recent result is provided
    if(is_proc(IE_))then
       tSimulationTmp = TimeCoupling
       call IE_run(tSimulationTmp,tSimulationTmp)
    end if

    call update_mapping(iMappingNew)
    if(iMappingNew/=iMappingIeGm)then
!BOP:
!DESCRIPTION:
!Constructor fot the router
!\begin{verbatim}
       call set_router(&
            GridDescriptorSource=GridIe,&
            GridDescriptorTarget=GridGmIe,&
            Router=RouterIeGm,&
            mapping=r_bound_to_r_iono,&
            interpolate=bilinear_interpolation)
!\end{verbatim}
!EOP
       iMappingIeGm=iMappingNew
    end if
    call couple_comp(RouterIeGm,    1,&
         IE_get_for_gm_swmf,GM_put_from_ie_swmf)
    if(DoTest.and.is_proc0(GM_)) call GM_print_variables('IE_swmf')
    if(is_proc(GM_))call transform_phi_bc_to_u_bc
  end subroutine couple_ie_gm_swmf


  !================================================================!
!BOP
!IROUTINE:sort_hemispherees - mapping for the router
!INTERFACE:
  subroutine r_bound_to_r_iono(&
       nDimFrom,XyzFrom_D,nDimTo,XyzTo_D,IsInterfacePoint)
    !USES:
    use CON_physics, ONLY: map_planet_field
!INPUT ARGUMENTS:
    integer,intent(in)::nDimFrom,nDimTo       
    real,dimension(nDimFrom),intent(in)::XyzFrom_D
!OUTPUT ARGUMENTS:
    real,dimension(nDimTo),intent(out)::XyzTo_D
    logical,intent(out)::IsInterfacePoint
!DESCRIPTION:
!The mapping used to construct the router is very simple in this
!case, because one should only find out, into which of the two
!hemispheres (North, South) the point XyzFrom\_D is projected along 
!the field line. Exact position of the point XyzTo\_D is not needed 
!EOP
    real::XyzSpherGm_D(2),XyzGm_D(3),XyzIe_D(3),XyzSmg_D(3)
    real::rBoundary,r2
    integer::iHemisphere
    real :: Sph_D(2)
!----------------------------------------------------------------------------------------

    rBoundary=Grid_C(GmIeGrid_)%Coord3_I(1)
    IsInterfacePoint=.true.
    call gen_to_stretched(XyzFrom_D,XyzSpherGm_D,2,GmIeGrid_)
    XyzGm_D(1)=sin(XyzSpherGm_D(1))*cos(XyzSpherGm_D(2))*rBoundary
    XyzGm_D(2)=sin(XyzSpherGm_D(1))*sin(XyzSpherGm_D(2))*rBoundary
    XyzGm_D(3)=cos(XyzSpherGm_D(1))*rBoundary
    XyzSmg_D=matmul(GmIe_DD, XyzGm_D)
    if(rBoundary>=1.05)then
       ! Get the mapping point
       call map_planet_field(TimeCoupling,XyzSmg_D,'SMG NORM',rIonosphere, &
            XyzIe_D,iHemisphere)
    else
       XyzIe_D= XyzSmg_D
    end if
    
    !Neither atan2 nor atan2_check are sufficiently stable for the 
    !application below:
    if(all(abs(XyzIe_D(1:2))<cTolerance))then
       Sph_D(1)=cPi*(cHalf-sign(cHalf,XyzIe_D(3)))
       Sph_D(2)=cZero
    else
       Sph_D(1)=atan2(&
            sqrt(dot_product(XyzIe_D(1:2),XyzIe_D(1:2))),&
            XyzIe_D(3))
       Sph_D(2)=atan2(XyzIe_D(2),XyzIe_D(1))
       if(Sph_D(2)<cZero)Sph_D(2)=Sph_D(2)+cTwoPi
    end if
    call stretched_to_gen(Sph_D,XyzTo_D,2,GmIeGrid_)
  end subroutine r_bound_to_r_iono
!BOP
!INTERFACE:
  subroutine update_mapping(iMappingNew)
    use CON_physics,ONLY:transform_matrix
!OUTPUT ARGUMENTS:
    integer,intent(out)::iMappingNew
!EOP
    real::AxisMagNew_D(3),GmIeNew_DD(3,3)
    GmIeNew_DD=transform_matrix(TimeCoupling,Grid_C(GM_)%TypeCoord,'SMG')
    AxisMagNew_D=matmul((/cZero,cZero,cOne/),GmIeNew_DD)
    if(iMappingGmIe==-1.or.&
         dot_product(AxisMagNew_D,AxisMag_D)<cOne-cTiny)then
       iMappingNew=iMappingGmIe+1
       AxisMag_D=AxisMagNew_D
       GmIe_DD=GmIeNew_DD
    else
       iMappingNew=iMappingGmIe
    end if
    if(DoTest.and.is_proc0(GM_))write(*,*)'Axis and transformation matrix:',&
         AxisMag_D, GmIe_DD
  end subroutine update_mapping
  !=====================================================================!
end module CON_couple_gm_ie_swmf
