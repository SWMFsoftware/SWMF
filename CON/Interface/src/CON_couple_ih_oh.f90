!^CMP COPYRIGHT UM
!^CMP FILE OH
!^CMP FILE IH
!BOP
!MODULE: CON_couple_ih_oh - couple OH to IH outer boundary (one way)
!INTERFACE:
module CON_couple_ih_oh

  !DESCRIPTION:
  ! This coupler uses the SWMF parallel coupling toolkit.
  ! The IH grid is coupled to a buffer grid in OH. The buffer grid
  ! uses the same coordinate system as OH, so the transformation is
  ! done in the IH wrapper.
  !
  ! The OH grid is coupled to the outer ghost cells of the IH grid directly.
  ! Both OH and IH use AMR grids, the buffer is a simple spherical grid.
  
  !USES:
  use CON_coupler
  use CON_axes, ONLY: transform_matrix,transform_velocity
  use ModConst

  implicit none
  private !except
  !
  !PUBLIC MEMBER FUNCTIONS:
  public:: couple_oh_ih_init
  public:: couple_oh_ih,couple_ih_oh

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V.<igorsok@umich.edu> - prototype for ih-gm
  ! 7/04/04                                 - version for ih-sc
  ! 7/20/04                                 - version for sc-buffer
  ! 6/03/08 Oran R. <oran@umich.edu>        - version for ih-oh
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: OH_iGridRealization=-2
  integer :: IH_iGridInIhOh=-2
  integer :: IH_iGridInOhIh=-2

  ! OH <-> IH conversion matrices
  real :: OhToIh_DD(3,3),IhToOh_DD(3,3)

  ! Maximum time difference [s] without remap 
  ! The 600 s corresponds to about 0.1 degree rotation between IH and OH
  real :: dTimeMappingMax = 600.0
 
  type(RouterType),save             :: RouterOhIh
  type(RouterType),save             :: RouterIhBuff
  type(GridDescriptorType),save     :: IH_SourceGrid
  type(GridDescriptorType),save     :: IH_TargetGrid
  type(GridDescriptorType),save     :: OH_Grid
  type(GridDescriptorType),save     :: BuffGD
  type(DomainDecompositionType),&
       save,target                  :: BuffDD
  logical :: DoInitialize=.true., DoTest, DoTestMe
  real :: tNow
  character(len=*), parameter :: NameMod='couple_oh_ih'
  logical::IsSphericalIh=.false.
  integer::iError
contains

  !===============================================================!
  subroutine couple_oh_ih_init
    interface
       subroutine OH_set_buffer_grid(Dd)
         use CON_domain_decomposition
         implicit none
         type(DomainDecompositionType),&
              intent(out)::Dd
       end subroutine OH_set_buffer_grid
    end interface

    if(.not.DoInitialize)return
    DoInitialize=.false.
    
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    call init_coupler(              &    
       iCompSource=OH_,             & ! component index for source
       iCompTarget=IH_,             & ! component index for target
       nGhostPointTarget=2,         & ! number of halo points in target
       GridDescriptorSource=OH_Grid,& ! OUT!\
       GridDescriptorTarget=IH_TargetGrid,& !-General coupler variables 
       Router=RouterOhIh)             ! OUT!/
    
   
    OH_iGridRealization=-1
    IH_iGridInOhIh     =-1
    IH_iGridInIhOh     =-1
    
    call OH_set_buffer_grid(BuffDD)
    call set_standard_grid_descriptor(&
         BuffDD,          &
         Standard_=Nodes_,&
         nGhostGridPoints=1,  &
         GridDescriptor=BuffGD)
    call set_standard_grid_descriptor(&
         IH_,GridDescriptor=IH_SourceGrid)
    call init_buffer_grid_couple(&
         SourceGD=IH_SourceGrid,&
         TargetGD=BuffGD, &
         RouterToBuffer=RouterIhBuff,    &
         nVar=8,&
         NameBuffer='OH_from_ih')  !Version for the first order in time
  end subroutine couple_oh_ih_init
  !===============================================================!
  !BOP
  !IROUTINE: couple_oh_ih - get OH solution at IH outer ghostpoints
  !INTERFACE:
  subroutine couple_oh_ih(TimeCoupling)
    !INPUT ARGUMENTS:
    interface
       subroutine IH_put_from_mh(nPartial,&
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
       end subroutine IH_put_from_mh
    end interface

    real,intent(in)::TimeCoupling
    !EOP

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !-------------------------------------------------------------------------

    if(.not.RouterOhIh%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call OH_synchronize_refinement(RouterOhIh%iProc0Source,&
                                   RouterOhIh%iCommUnion)
    call IH_synchronize_refinement(RouterOhIh%iProc0Target,&
                                   RouterOhIh%iCommUnion)

    ! Redo the router if any of the grids changed or 
    ! the two coordinate systems have rotated away too much
    if(OH_iGridRealization/=i_realization(OH_).or.&     
         IH_iGridInOhIh/=i_realization(IH_) .or. &
         (Grid_C(OH_) % TypeCoord /= Grid_C(IH_) % TypeCoord &
         .and. TimeCoupling - TimeCouplingLast > dTimeMappingMax)) then  
       ! Recalculate the IH to OH transformation matrix used in mapping
       ! a target point (IH) to the source (OH)
       !(it is time dependent in general)

       IhToOh_DD = transform_matrix(TimeCoupling, &
            Grid_C(IH_) % TypeCoord, Grid_C(OH_) % TypeCoord)

       !Recalculate the tramsposed matrix, used in transforming vectors
       !of velocity and magnetic field to be sent from OH to IH
       OhToIh_DD = transpose(IhToOh_DD)

       call set_router(&
            GridDescriptorSource=OH_Grid,&
            GridDescriptorTarget=IH_TargetGrid,&
            Router=RouterOhIh,&
            is_interface_block=is_boundary_block,&
            interface_point_coords=outer_cells, &
            mapping=map_ih_oh, &
            interpolate=interpolation_fix_reschange)

       OH_iGridRealization = i_realization(OH_)
       IH_iGridInOhIh      = i_realization(IH_)
       TimeCouplingLast    = TimeCoupling
       tNow=TimeCoupling
    end if
    call couple_comp(&
         RouterOhIh,&
         nVar=8,&
         fill_buffer=OH_get_for_ih_and_transform,&
         apply_buffer=IH_put_from_mh)

  end subroutine couple_oh_ih
  !======================================================!
  logical function is_boundary_block(lGlobalTreeNode)
    integer,parameter::R_=1
    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsBoundary_D

    IsBoundary_D=is_right_boundary_d(&
         IH_TargetGrid%DD%Ptr,lGlobalTreeNode)
    !For spherical domain
    is_boundary_block=IsBoundary_D(R_).and.IsSphericalIh

    !For cartesian box   
     IsBoundary_D= IsBoundary_D.or.&
         is_left_boundary_d(&
         IH_TargetGrid%DD%Ptr,lGlobalTreeNode)

    is_boundary_block=is_boundary_block.or. &
         (any(IsBoundary_D).and.(.not.IsSphericalIh))
   
    
  end function is_boundary_block
  !===============================================================!
  subroutine outer_cells(&
       GridDescriptor,&
       lGlobalTreeNode,&
       nDim,&
       Xyz_D,&
       nIndexes,&
       i_D,&
       IsInterfacePoint)

    type(GridDescriptorType),intent(in):: GridDescriptor
    integer,intent(in)::lGlobalTreeNode,nIndexes
    logical,intent(out)::IsInterfacePoint
    integer,intent(in)::nDim
    real,intent(inout)::Xyz_D(nDim)
    integer,intent(inout)::i_D(nIndexes)

    logical,dimension(3)::IsLeftFace_D,IsRightFace_D
    integer,parameter::x_=1,y_=2,z_=3
    integer,parameter::R_=1
   
    IsLeftFace_D=i_D(x_:z_)<1.and.is_left_boundary_d(&
         IH_TargetGrid%DD%Ptr,lGlobalTreeNode)
    IsRightFace_D=i_D(x_:z_)>&
         ncells_decomposition_d(IH_TargetGrid%DD%Ptr).and.&
         is_right_boundary_d(IH_TargetGrid%DD%Ptr,lGlobalTreeNode)

    !For spherical grid
    IsInterfacePoint=IsRightFace_D(R_).and.IsSphericalIh

    !For Cartesian grid:

    IsInterfacePoint=IsInterfacePoint.or.&
         (any(IsRightFace_D.or.IsLeftFace_D).and.(.not.IsSphericalIh))
  end subroutine outer_cells 
  !========================================================!
  subroutine map_ih_oh(&
       IH_nDim,IH_Xyz_D,OH_nDim,OH_Xyz_D,IsInterfacePoint)

    integer,intent(in)::OH_nDim,IH_nDim
    real,dimension(IH_nDim),intent(in)::IH_Xyz_D
    real,dimension(OH_nDim),intent(out)::OH_Xyz_D
    logical,intent(out)::IsInterfacePoint

    OH_Xyz_D = matmul(IhToOh_DD, IH_Xyz_D)*&
         Grid_C(IH_)%UnitX/Grid_C(OH_)%UnitX

    IsInterfacePoint=.true.

  end subroutine map_ih_oh
  !=================================================================!

  subroutine OH_get_for_ih_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7,&
         BuffX_=9,BuffZ_=11
    !------------------------------------------------------------
    call OH_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)

    !Transform velocity
    State_V(RhoUx_:RhoUz_)=State_V(Rho_)*&
         transform_velocity(tNow,&
         State_V(RhoUx_:RhoUz_)/State_V(Rho_),&
         State3_V(BuffX_:BuffZ_)/State_V(Rho_),&
         Grid_C(OH_)%TypeCoord,Grid_C(IH_)%TypeCoord)
    
    State_V(Bx_:Bz_)=matmul(OhToIh_DD,State_V(Bx_:Bz_))

  end subroutine OH_get_for_ih_and_transform

!===============================================================!
!===============================================================!
!===============================================================!
!===============================================================
!BOP
!IROUTINE: couple_ih_oh - interpolate and get MHD state at the OH buffer grid 
!INTERFACE:
  subroutine couple_ih_oh(TimeCoupling)
    use ModIoUnit
    !INPUT ARGUMENTS:
    interface
       subroutine IH_get_for_mh(&
            nPartial,iGetStart,Get,w,State_V,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::w
         real,dimension(nVar),intent(out)::State_V
       end subroutine IH_get_for_mh
    end interface
    integer::iPoint,nU_I(2)
    real,intent(in)::TimeCoupling
    integer,save::iCoupling=0
    integer::iFile
    character(LEN=21)::NameFile
    logical::DoneMatchIBC=.false.
!EOP
    if(.not.RouterIhBuff%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    call IH_synchronize_refinement(RouterIhBuff%iProc0Source,&
                                   RouterIhBuff%iCommUnion)


    if(IH_iGridInIhOh/=i_realization(IH_))then  
       call set_router(&
            GridDescriptorSource=IH_SourceGrid,&
            GridDescriptorTarget=BuffGD,&
            Router=RouterIhBuff,&
            mapping=buffer_grid_point,&
            interpolate=interpolation_fix_reschange)
       IH_iGridInIhOh= i_realization(IH_)
    end if

    call couple_buffer_grid(&
         RouterIhBuff,&
         nVar=8,&
         fill_buffer=IH_get_for_mh,&
         NameBuffer='OH_from_ih',&
         TargetID_=OH_)
    if(.not.DoneMatchIBC)then
       DoneMatchIBC=.true.
       if(is_proc(OH_))call OH_match_ibc
    end if
    if(DoTest.and.is_proc0(compid_grid(BuffGD%DD%Ptr)))then
       nU_I=ubound_vector('OH_from_ih')
       iCoupling=iCoupling+1
       iFile=io_unit_new()
       write(NameFile,'(a,i4.4,a)')'./OH/from_ih_',iCoupling,'.dat'
       open(iFile,FILE=NameFile,STATUS='unknown')
       do iPoint=1,nU_I(2)
          write(iFile,*)point_state_v('OH_from_ih',8,iPoint)
       end do
       close(iFile)
    end if
    if(DoTest)write(*,*)'Couple passed at PE=',i_proc()
  end subroutine couple_ih_oh
  !======================================================!
  subroutine buffer_grid_point(&
       nDimFrom,Sph_D,nDimTo,IH_Xyz_D,IsInterfacePoint)         

    ! Transform from the spherical buffer grid to the Cartesian IH grid

    integer,intent(in)                    :: nDimFrom, nDimTo       
    real,dimension(nDimFrom), intent(in)  :: Sph_D
    real,dimension(nDimTo),   intent(out) :: IH_Xyz_D
    logical,intent(out)::IsInterfacePoint

    ! The order of spherical indexes as in BATSRUS
    ! This is a left handed system !!! should be replaced !!!
    integer,parameter::r_=1,Psi_=2,Theta_=3,x_=1,y_=2,z_=3
    real :: rSinTheta
    !-----------------------------------------------------------------------
    RSinTheta    = Sph_D(r_)*sin(Sph_D(Theta_))!To be modified
    IH_Xyz_D(x_) = RSinTheta*cos(Sph_D(Psi_))  !\if IH grid
    IH_Xyz_D(y_) = RSinTheta*sin(Sph_D(Psi_))  !/is spherical
    IH_Xyz_D(z_) = Sph_D(r_)*cos(Sph_D(Theta_))!To be modified

    IsInterfacePoint=.true.
  end subroutine buffer_grid_point
  
end module CON_couple_ih_oh

