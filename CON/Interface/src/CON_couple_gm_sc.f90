!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE SC
!^CMP FILE GM

module CON_couple_gm_sc

  ! This coupler uses the SWMF coupling toolkit.
  ! Both the SC and GM grids use AMR.
  ! SC is a source, GM is a target.
  ! Typically SC uses the solar radius rSun as distance unit,
  ! while GM uses the Earth's radius as distance unit.
  !
  ! The mapping is characterized by rotation matrices and/or
  ! translation vectors.
  ! In the particular case this is the position of the origin of
  ! the GM frame of reference in the SC frame of refernce  and the
  ! matrix which characterizes the rotation of the SC frame of
  ! reference with respect to GM:
  !
  ! $$
  ! {\bf V}_{SC}={\bf A}_{SC,GM}\cdot{\bf V}_{GM},
  ! {\bf V}_{GM}={\bf A}_{GM,SC}\cdot{\bf V}_{SC}
  ! $$

  use CON_coupler
  use CON_transfer_data, ONLY: transfer_real_array, transfer_integer
  use CON_time, ONLY: TimeStart
  use ModConst
  use CON_axes, ONLY: transform_matrix, vPlanetHgi_D, XyzPlanetHgi_D

  use SC_wrapper, ONLY: SC_synchronize_refinement, SC_get_for_gm, &
       SC_save_global_buffer, SC_set_buffer_grid_get_info, SC_xyz_to_coord

  use GM_wrapper, ONLY: GM_synchronize_refinement, GM_put_from_mh, &
       GM_is_right_boundary_d, GM_get_for_global_buffer

  implicit none

  save

  private ! except

  public:: couple_sc_gm_init
  public:: couple_sc_gm
  public:: couple_gm_sc

  real, public:: CouplingTimeScGm

  ! revision history:
  ! 7/23/03 Sokolov I.V. <igorsok@umich.edu> - initial prototype
  ! 9/02/03 G.Toth <gtoth@umich.edu> - minor changes
  ! 6/18/05 Sokolov - generalized to arbitrary coordinate systems
  ! 6/22/05 G.Toth  - redo mapping and router to allow for relative motion

  ! To trace the possible changes in the grids and/or mapping
  integer :: SC_iGridRealization=-2
  integer :: GM_iGridRealization=-2

  ! Three-dimensional coordinate vector
  integer, parameter :: nDim = 3
  real, dimension(nDim)   :: XyzPlanetSc_D, vPlanetSc_D
  real, dimension(nDim,nDim) :: GmToSc_DD, ScToGm_DD, HgiToSc_DD

  ! Maximum time difference [s] without remap
  ! The 600 s corresponds to about 0.1 degree rotation between SC and GM
  real :: dTimeMappingMax = 600.0

  type(RouterType)         :: Router
  type(GridType) :: GM_Grid
  type(GridType) :: SC_Grid
  type(LocalGridType)        :: GM_LocalGrid

  logical :: DoInitialize=.true., DoTest, DoTestMe

  ! Size and limits of the 3D spherical buffer grid
  integer :: iSize, jSize, kSize
  real    :: BufferMinMaxGm_DI(nDim,2)

  character(len=*), parameter :: NameMod='CON_couple_sc_gm'

contains
  !============================================================================
  subroutine couple_sc_gm_init
    !--------------------------------------------------------------------------
    if(.not.DoInitialize)RETURN
    DoInitialize=.false.

    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_sc_gm_init iProc=',i_proc()

    call init_coupler(              &
         iCompSource=SC_,             & ! component index for source
         iCompTarget=GM_,             & ! component index for target
         nGhostPointTarget=2,         & ! number of halo points in target
         GridSource=SC_Grid,& ! OUT
         GridTarget=GM_Grid,& ! OUT!-General coupler variables
         LocalGridTarget  =GM_LocalGrid,  & ! OUT!-optional
         Router=Router)                 ! OUT
    SC_iGridRealization=-1
    GM_iGridRealization=-1
    call set_couple_var_info(iCompSource=GM_, iCompTarget=SC_)
    ! Set buffer grid location and size in SC, and retrieve them for coupler
    if(is_proc(SC_)) then
       call SC_set_buffer_grid_get_info(iSize, jSize, kSize, BufferMinMaxGm_DI)

       ! Convert units for radial coordinate  before passing to Gm
       BufferMinMaxGm_DI(1,:) = BufferMinMaxGm_DI(1,:) &
            *(Grid_C(SC_)%UnitX/Grid_C(Gm_)%UnitX)
    end if

    ! Pass buffer size
    call transfer_integer(SC_, GM_, iSize, jSize, kSize, &
         UseSourceRootOnly = .false.)

    ! Pass buffer boundary info
    call transfer_real_array(SC_, GM_, 6, BufferMinMaxGm_DI, &
         UseSourceRootOnly = .false.)

  end subroutine couple_sc_gm_init
  !============================================================================
  subroutine couple_sc_gm(TimeCoupling)
    real,intent(in)::TimeCoupling

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameMod,DoTest,DoTestMe)
    if(DoTest)write(*,*)'couple_sc_gm iProc=',i_proc()

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call SC_synchronize_refinement(Router%iProc0Source,Router%iCommUnion)
    call GM_synchronize_refinement(Router%iProc0Target,Router%iCommUnion)

    CouplingTimeScGm=TimeCoupling
    call set_couple_var_info(SC_, GM_)
    if(SC_iGridRealization/=i_realization(SC_).or.&
         GM_iGridRealization/=i_realization(GM_).or.&
         TimeCoupling - TimeCouplingLast > dTimeMappingMax)then
       if(GM_iGridRealization/=i_realization(GM_))then

          ! reset local Grid
          call clean_gd(GM_LocalGrid)
          call set_local_gd(i_proc(), GM_grid, GM_LocalGrid)
       end if
       ! Set the transformation matrices, compute the position and
       ! orbital veloctiy of the Earth in SC at the coupling time.

       GmToSc_DD  = transform_matrix(TimeCoupling , &
            Grid_C(GM_) % TypeCoord, Grid_C(SC_) % TypeCoord)
       ScToGm_DD  = transpose(GmToSc_DD)

       HgiToSc_DD = transform_matrix(TimeCoupling, &
            'HGI', Grid_C(SC_) % TypeCoord)

       ! Transform the Earth position in HGI to SC coordinates
       ! and change distance units to that used in SC

       XyzPlanetSc_D = matmul(HgiToSc_DD, XyzPlanetHgi_D)/Grid_C(SC_)%UnitX

       ! Transform the Earth velocity from Hgi to SC coordinates

       vPlanetSc_D = matmul(HgiToSc_DD, vPlanetHgi_D)

       if(DoTest)write(*,*)'couple_sc_gm call set_router iProc=',i_proc()

       if(DoTestMe)write(*,*)'couple_sc_gm_init XyzPlanetSc_D=',&
            XyzPlanetSc_D
       call set_router(&
            GridSource=SC_Grid,&
            GridTarget=GM_LocalGrid,&
            Router=Router,&
            is_interface_block=GM_is_west_block,&
            interface_point_coords=GM_west_cells, &
            mapping=map_gm_sc,&
            interpolate=interpolation_amr_gc)

       SC_iGridRealization = i_realization(SC_)
       GM_iGridRealization = i_realization(GM_)
       TimeCouplingLast    = TimeCoupling
    end if

    call couple_comp(&
         Router,&
         nVar=nVarCouple,&
         fill_buffer=SC_get_for_gm_and_transform,&
         apply_buffer=GM_put_from_mh)

  end subroutine couple_sc_gm
  !============================================================================
  logical function GM_is_west_block(iBlockLocal)

    integer,intent(in)  :: iBlockLocal
    logical             :: IsRightBoundary_D(GM_grid%nDim)
    integer, parameter:: x_=1
    !--------------------------------------------------------------------------
    IsRightBoundary_D = GM_is_right_boundary_d(iBlockLocal)
    GM_is_west_block=IsRightBoundary_D(x_)

  end function GM_is_west_block
  !============================================================================
  subroutine GM_west_cells(nDim, Xyz_D, nIndex, i_D, IsInterfacePoint)

    integer, intent(in):: nDim, nIndex
    real,    intent(inout):: Xyz_D(nDim)
    integer, intent(inout)       :: i_D(nIndex)
    logical, intent(out)         :: IsInterfacePoint

    logical:: IsLeftFace_D(nDim), IsRightFace_D(nDim)
    integer, parameter:: x_=1, y_=2, z_=3
    !--------------------------------------------------------------------------
    IsLeftFace_D = i_D(x_:z_)  < 1
    IsRightFace_D = i_D(x_:z_) > GM_Grid%Domain%Ptr%nCell_D
    IsInterfacePoint = IsRightFace_D(x_) .and. &
         .not.(any(IsLeftFace_D(y_:z_)) .or. any(IsRightFace_D(y_:z_)))

  end subroutine GM_west_cells
  !============================================================================
  subroutine map_gm_sc( &
       GM_nDim, XyzGm_D, SC_nDim, CoordSc_D, IsInterfacePoint)

    integer,intent(in) :: GM_nDim, SC_nDim
    real,   intent(in) :: XyzGm_D(GM_nDim)
    real,   intent(out):: CoordSc_D(SC_nDim)
    logical,intent(out)::IsInterfacePoint

    ! In each mapping the corrdinates of the TARGET grid point (GM)
    ! shoud be be transformed to the SOURCE (SC) generalized coords.

    real :: XyzSc_D(nDim)
    !--------------------------------------------------------------------------
    XyzSc_D = XyzPlanetSc_D + matmul(GmToSc_DD, XyzGm_D)*&
         Grid_C(GM_)%UnitX/Grid_C(SC_)%UnitX
    call SC_xyz_to_coord(XyzSc_D, CoordSc_D)
    IsInterfacePoint=.true.

  end subroutine map_gm_sc
  !============================================================================
  subroutine SC_get_for_gm_and_transform( &
       nPartial, iGetStart, Get, w, State_V, nVar)

    integer, intent(in):: nPartial, iGetStart, nVar
    type(IndexPtrType), intent(in):: Get
    type(WeightPtrType), intent(in):: w
    real, intent(out)::State_V(nVar)

    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7
    !--------------------------------------------------------------------------
    call SC_get_for_gm(&
         nPartial,iGetStart,Get,w,State_V,nVar,CouplingTimeScGm)

    State_V(RhoUx_:RhoUz_)=&
         matmul(ScToGm_DD,State_V(RhoUx_:RhoUz_) &
         - State_V(Rho_)*vPlanetSc_D )
    State_V(Bx_:Bz_)=matmul(ScToGm_DD,State_V(Bx_:Bz_))

  end subroutine SC_get_for_gm_and_transform
  !============================================================================
  subroutine couple_gm_sc(TimeCoupling)

    real, intent(in) :: TimeCoupling     ! simulation time at coupling

    ! Couple between two components:
    !    General Magnetosphere (GM)  source
    !    Solar Corona           (SC)  target
    !
    ! The GM component sends the state variables to a buffer grid.
    ! SC uses the buffer grid to calculate the boundary condition on Body2.

    ! Array to store state vector on all buffer grid points
    real, allocatable :: Buffer_VIII(:,:,:,:)

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_gm_sc'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(DoTest.and.is_proc0(GM_))&
         write(*,'(a,es12.5)')NameSub//': starting, Time=', TimeCoupling
    ! Transfer buffer grid from GM to SC to be used for inner boundary
    allocate(Buffer_VIII(nVarCouple,iSize,jSize,kSize))
    if(is_proc(GM_)) call GM_get_for_global_buffer(iSize, jSize, kSize, &
         BufferMinMaxGm_DI, Buffer_VIII)

    ! Add up Buffer on SC processors and transfer to IH
    call transfer_real_array(GM_, SC_, size(Buffer_VIII), Buffer_VIII, &
         UseSourceSum=.true.)

    if(is_proc(SC_)) call SC_save_global_buffer( &
         nVarCouple, iSize, jSize, kSize, Buffer_VIII)
    deallocate(Buffer_VIII)

    if(DoTest.and.is_proc0(GM_))&
         write(*,'(a,es12.5)')NameSub//': finished, Time=', TimeCoupling

  end subroutine couple_gm_sc
  !============================================================================
end module CON_couple_gm_sc
!==============================================================================
