!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IH
!^CMP FILE GM

module CON_couple_gm_ih

  ! This coupler uses the SWMF coupling toolkit.
  ! Both the IH and GM grids use AMR.
  ! IH is a source, GM is a target.
  ! Typically IH uses the solar radius rSun as distance unit,
  ! while GM uses the Earth's radius as distance unit.
  !
  ! The mapping is characterized by rotation matrices and/or
  ! translation vectors.
  ! In the particular case this is the position of the origin of
  ! the GM frame of reference in the IH frame of refernce  and the
  ! matrix which characterizes the rotation of the IH frame of
  ! reference with respect to GM:
  !
  ! $$
  ! {\bf V}_{IH}={\bf A}_{IH,GM}\cdot{\bf V}_{GM},
  ! {\bf V}_{GM}={\bf A}_{GM,IH}\cdot{\bf V}_{IH}
  ! $$
  use CON_coupler
  use CON_time, ONLY:TimeStart
  use ModConst
  use CON_axes, ONLY: transform_matrix, vPlanetHgi_D, XyzPlanetHgi_D

  use IH_wrapper, ONLY: IH_synchronize_refinement, &
       IH_get_for_gm

  use GM_wrapper, ONLY: GM_synchronize_refinement, &
       GM_put_from_mh, GM_is_right_boundary_d

  implicit none

  save

  private ! except

  public:: couple_ih_gm_init
  public:: couple_ih_gm

  real, public:: CouplingTimeIhGm

  ! revision history:
  ! 7/23/03 Sokolov I.V. <igorsok@umich.edu> - initial prototype
  ! 9/02/03 G.Toth <gtoth@umich.edu> - minor changes
  ! 6/18/05 Sokolov - generalized to arbitrary coordinate systems
  ! 6/22/05 G.Toth  - redo mapping and router to allow for relative motion

  ! To trace the possible changes in the grids and/or mapping
  integer :: IH_iGridRealization=-2
  integer :: GM_iGridRealization=-2

  real, dimension(3)   :: XyzPlanetIh_D, vPlanetIh_D
  real, dimension(3,3) :: GmToIh_DD, IhToGm_DD, HgiToIh_DD

  ! Maximum time difference [s] without remap
  ! The 600 s corresponds to about 0.1 degree rotation between IH and GM
  real :: dTimeMappingMax = 600.0

  type(RouterType)         :: Router
  type(GridType) :: GM_Grid
  type(GridType) :: IH_Grid
  type(LocalGridType)        :: GM_LocalGrid

  logical :: DoInitialize=.true., DoTest, DoTestMe

  character(len=*), parameter :: NameMod='CON_couple_gm_ih'

contains
  !============================================================================
  subroutine couple_ih_gm_init
    !--------------------------------------------------------------------------
    if(.not.DoInitialize)RETURN
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm_init iProc=',i_proc()

    call init_coupler(              &
         iCompSource=IH_,             & ! component index for source
         iCompTarget=GM_,             & ! component index for target
         nGhostPointTarget=2,         & ! number of halo points in target
         GridSource=IH_Grid,& ! OUT
         GridTarget=GM_Grid,& ! OUT!-General coupler variables
         LocalGridTarget  =GM_LocalGrid,  & ! OUT!-optional
         Router=Router)                 ! OUT
    IH_iGridRealization=-1
    GM_iGridRealization=-1
    call set_couple_var_info(iCompSource=IH_, iCompTarget=GM_)

  end subroutine couple_ih_gm_init
  !============================================================================
  subroutine couple_ih_gm(TimeCoupling)
    real,intent(in)::TimeCoupling

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_ih_gm iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call GM_synchronize_refinement(Router%iProc0Target,Router%iCommUnion)

    CouplingTimeIhGm=TimeCoupling
    call set_couple_var_info(IH_, GM_)
    if(IH_iGridRealization/=i_realization(IH_).or.&
         GM_iGridRealization/=i_realization(GM_).or.&
         TimeCoupling - TimeCouplingLast > dTimeMappingMax)then
       if(DoTestMe)then
          write(*,*)'GM_iGridRealization=',GM_iGridRealization
          write(*,*)'i_realization(GM_)=',i_realization(GM_)
          write(*,*)'IH_iGridRealization=',IH_iGridRealization
          write(*,*)'i_realization(IH_)=',i_realization(IH_)
          write(*,*)'TimeCoupling, TimeCouplingLast=',&
               TimeCoupling, TimeCouplingLast
       end if
       if(GM_iGridRealization/=i_realization(GM_))then

          ! reset local Grid
          call clean_gd(GM_LocalGrid)
          call set_local_gd(i_proc(), GM_grid, GM_LocalGrid)
       end if
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
            GridSource=IH_Grid,&
            GridTarget=GM_LocalGrid,&
            Router=Router,&
            is_interface_block=GM_is_west_block,&
            interface_point_coords=GM_west_cells, &
            mapping=map_gm_ih,&
            interpolate=interpolation_amr_gc)

       IH_iGridRealization = i_realization(IH_)
       GM_iGridRealization = i_realization(GM_)
       TimeCouplingLast    = TimeCoupling
    end if

    call couple_comp(&
         Router,&
         nVar=nVarCouple,&
         fill_buffer=IH_get_for_gm_and_transform,&
         apply_buffer=GM_put_from_mh)

  end subroutine couple_ih_gm
  !============================================================================
  logical function GM_is_west_block(iBlockLocal)

    integer,intent(in):: iBlockLocal

    logical:: IsRightBoundary_D(GM_grid%nDim)
    integer, parameter:: x_=1
    !--------------------------------------------------------------------------
    IsRightBoundary_D = GM_is_right_boundary_d(iBlockLocal)
    GM_is_west_block=IsRightBoundary_D(x_)

  end function GM_is_west_block
  !============================================================================
  subroutine GM_west_cells(nDim, Xyz_D, nIndex, i_D, IsInterfacePoint)
    integer, intent(in):: nDim, nIndex
    real,    intent(inout):: Xyz_D(nDim)
    integer, intent(inout):: i_D(nIndex)
    logical, intent(out)  :: IsInterfacePoint

    logical:: IsLeftFace_D(3), IsRightFace_D(3)
    integer, parameter:: x_=1, y_=2, z_=3
    !--------------------------------------------------------------------------
    IsLeftFace_D=i_D(x_:z_)  < 1
    IsRightFace_D=i_D(x_:z_) > GM_Grid%Domain%Ptr%nCell_D
    IsInterfacePoint=IsRightFace_D(x_).and..not.&
         (any(IsLeftFace_D(y_:z_)).or.any(IsRightFace_D(y_:z_)))

  end subroutine GM_west_cells
  !============================================================================
  subroutine map_gm_ih( &
       GM_nDim, XyzGm_D, IH_nDim, XyzIh_D, IsInterfacePoint)

    integer, intent(in) :: GM_nDim, IH_nDim
    real,    intent(in) :: XyzGm_D(GM_nDim)
    real,    intent(out):: XyzIh_D(IH_nDim)
    logical, intent(out):: IsInterfacePoint

    ! In each mapping the corrdinates of the TARGET grid point (GM)
    ! shoud be be transformed to the SOURCE (IH) generalized coords.
    !--------------------------------------------------------------------------
    XyzIh_D = XyzPlanetIh_D + matmul(GmToIh_DD, XyzGm_D)*&
         Grid_C(GM_)%UnitX/Grid_C(IH_)%UnitX
    IsInterfacePoint=.true.

  end subroutine map_gm_ih
  !============================================================================
  subroutine IH_get_for_gm_and_transform( &
       nPartial, iGetStart, Get, w, State_V, nVar)

    integer, intent(in)::nPartial, iGetStart, nVar
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType), intent(in):: w
    real, intent(out):: State_V(nVar)

    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7
    !--------------------------------------------------------------------------
    call IH_get_for_gm(&
         nPartial,iGetStart,Get,w,State_V,nVar,CouplingTimeIhGm)

    State_V(RhoUx_:RhoUz_)=&
         matmul(IhToGm_DD,State_V(RhoUx_:RhoUz_) &
         - State_V(Rho_)*vPlanetIh_D )
    State_V(Bx_:Bz_)=matmul(IhToGm_DD,State_V(Bx_:Bz_))

  end subroutine IH_get_for_gm_and_transform
  !============================================================================
end module CON_couple_gm_ih
!==============================================================================
