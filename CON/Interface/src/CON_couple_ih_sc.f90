!^CMP COPYRIGHT UM
!^CMP FILE IH
!^CMP FILE SC
!BOP
!MODULE: CON_couple_ih_sc - couple IH and SC both ways
!INTERFACE:
module CON_couple_ih_sc

  !DESCRIPTION:
  ! This coupler uses the SWMF parallel coupling toolkit.
  ! The SC grid is coupled to a buffer grid in IH. The buffer grid
  ! uses the same coordinate system as SC, so the transformation is
  ! done in the IH wrapper.
  !
  ! The IH grid is coupled to the outer ghost cells of the SC grid directly.
  ! Both SC and IH use AMR grids, the buffer is a simple spherical grid.
  
  !USES:
  use CON_coupler
  use CON_axes, ONLY: transform_matrix,transform_velocity
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

  ! IH <-> SC conversion matrices
  real :: IhToSc_DD(3,3),ScToIh_DD(3,3)

  ! Maximum time difference [s] without remap 
  ! The 600 s corresponds to about 0.1 degree rotation between SC and IH
  real :: dTimeMappingMax = 600.0
 
  type(RouterType),save             :: RouterIhSc
  type(RouterType),save             :: RouterScBuff
  type(GridDescriptorType),save     :: SC_SourceGrid
  type(GridDescriptorType),save     :: SC_TargetGrid
  type(GridDescriptorType),save     :: IH_Grid
  type(GridDescriptorType),save     :: BuffGD
  type(DomainDecompositionType),&
       save,target                  :: BuffDD
  logical :: DoInitialize=.true., DoTest, DoTestMe
  real :: tNow
  character(len=*), parameter :: NameMod='couple_ih_sc'
  logical::IsSphericalSc=.false. , UseGenRSc = .false., UseLogRSc = .false.
  integer::iError
  
  !Parameters of the stretched grid, if needed
  integer :: nGenRGridSc
  real    :: DeltaGen

  ! Minimum number of variables in the different models
  integer :: nVarIhSc

contains
  !===============================================================!
  subroutine couple_ih_sc_init

    interface
       subroutine IH_set_buffer_grid(Dd)
         use CON_domain_decomposition
         implicit none
         type(DomainDecompositionType),&
              intent(out)::Dd
       end subroutine IH_set_buffer_grid
    end interface

    if(.not.DoInitialize)return
    DoInitialize=.false.
    
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    nVarIhSc = min(Grid_C(IH_)%nVar, Grid_C(SC_)%nVar)

    IsSphericalSc = index(Grid_C(SC_) % TypeGeometry,'spherical') > 0 
    UseLogRSc     = index(Grid_C(SC_) % TypeGeometry,'lnr'      ) > 0
    UseGenRSc     = index(Grid_C(SC_) % TypeGeometry,'genr'     ) > 0
    if(UseGenRSc) then
       nGenRGridSc = size(Grid_C(SC_) % Coord1_I,1)
       if(nGenRGridSc==1) &
            call CON_stop('Stretched grid in SC is not properly initialized')
       DeltaGen = 1.0/(nGenRGridSC - 1)
    end if

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
         nVar = nVarIhSc, &
         NameBuffer='IH_from_sc')  !Version for the first order in time

  end subroutine couple_ih_sc_init
  !===============================================================!
  !BOP
  !IROUTINE: couple_ih_sc - get IH solution at SC outer ghostpoints
  !INTERFACE:
  subroutine couple_ih_sc(TimeCoupling)
    !INPUT ARGUMENTS:
    interface
       subroutine SC_put_from_mh(nPartial,&
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
       end subroutine SC_put_from_mh
    end interface

    real,intent(in)::TimeCoupling
    !EOP

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !-------------------------------------------------------------------------

    if(.not.RouterIhSc%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call IH_synchronize_refinement(RouterIhSc%iProc0Source,&
                                   RouterIhSc%iCommUnion)
    call SC_synchronize_refinement(RouterIhSc%iProc0Target,&
                                   RouterIhSc%iCommUnion)

    ! Redo the router if any of the grids changed or 
    ! the two coordinate systems have rotated away too much
    if(IH_iGridRealization/=i_realization(IH_).or.&     
         SC_iGridInIhSc/=i_realization(SC_) .or. &
         (Grid_C(IH_) % TypeCoord /= Grid_C(SC_) % TypeCoord &
         .and. TimeCoupling - TimeCouplingLast > dTimeMappingMax)) then  
       ! Recalculate the SC to IH transformation matrix used in mapping
       ! a target point (SC) to the source (IH)
       !(it is time dependent in general)

       ScToIh_DD = transform_matrix(TimeCoupling, &
            Grid_C(SC_) % TypeCoord, Grid_C(IH_) % TypeCoord)

       !Recalculate the tramsposed matrix, used in transforming vectors
       !of velocity and magnetic field to be sent from IH to SC
       IhToSc_DD = transpose(ScToIh_DD)

       call set_router(&
            GridDescriptorSource=IH_Grid,&
            GridDescriptorTarget=SC_TargetGrid,&
            Router=RouterIhSc,&
            is_interface_block=is_boundary_block,&
            interface_point_coords=outer_cells, &
            mapping=map_sc_ih, &
            interpolate=interpolation_fix_reschange)

       IH_iGridRealization = i_realization(IH_)
       SC_iGridInIhSc      = i_realization(SC_)
       TimeCouplingLast    = TimeCoupling
       tNow=TimeCoupling
    end if

    call couple_comp(&
         RouterIhSc,&
         nVar = nVarIhSc, &
         fill_buffer=IH_get_for_sc_and_transform,&
         apply_buffer=SC_put_from_mh)

  end subroutine couple_ih_sc
  !======================================================!
  logical function is_boundary_block(lGlobalTreeNode)
    integer,parameter::R_=1
    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsBoundary_D

    IsBoundary_D=is_right_boundary_d(&
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    !For spherical domain
    is_boundary_block=IsBoundary_D(R_).and.IsSphericalSc
    
    !Now if IsSphericalSc ==.true. then is_boundary_block is .true.
    !only if the right block boundary along the radial coordinate
    !is the SC domain boundary.
    !Else is_boundary_block = .false.

    !For cartesian box   
     IsBoundary_D= IsBoundary_D.or.&
         is_left_boundary_d(&
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
     !Now IsBoundary_D(iDim) is true if any of the boundaries along 
     !the direction iDim is the SC boundary

    is_boundary_block=is_boundary_block.or. &
         (any(IsBoundary_D).and.(.not.IsSphericalSc))
    !Now if IsSphericalSc = .true. the value of is_boundary_block does
    !not change. Otherwise is_boundary_block is true if any of the
    !block boundaries is the boundary of the SC domain
   
    
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
         SC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    IsRightFace_D=i_D(x_:z_)>&
         ncells_decomposition_d(SC_TargetGrid%DD%Ptr).and.&
         is_right_boundary_d(SC_TargetGrid%DD%Ptr,lGlobalTreeNode)

    !For spherical grid
    IsInterfacePoint=IsRightFace_D(R_).and.IsSphericalSc

    !For Cartesian grid:

    IsInterfacePoint=IsInterfacePoint.or.&
         (any(IsRightFace_D.or.IsLeftFace_D).and.(.not.IsSphericalSc))
  end subroutine outer_cells 
  !========================================================!
  subroutine map_sc_ih(&
       SC_nDim,SC_XyzIn_D,IH_nDim,IH_Xyz_D,IsInterfacePoint)

    integer,intent(in)::IH_nDim,SC_nDim
    real,dimension(SC_nDim),intent(in)::SC_XyzIn_D
    real,dimension(IH_nDim),intent(out)::IH_Xyz_D
    logical,intent(out)::IsInterfacePoint
    
    real, dimension(SC_nDim) :: SC_Xyz_D
    integer, parameter :: R_ = 1, Phi_=2, Theta_ = 3, x_ = 1, y_ = 2, z_ = 3
    real :: R, Gen, Phi, Theta, rSinTheta
    !--------------------------------------- 
    !In each mapping the corrdinates of the TARGET grid point (SC)
    !shoud be be transformed to the SOURCE (IH) generalized coords.
    if(.not.IsSphericalSc)then
       SC_Xyz_D = SC_XyzIn_D
    else
       !transform to dimensionless cartesian Xyz
       R = SC_Xyz_D(R_)
       if(UseLogRSc)then
          R = exp(R)
       elseif(UseGenRSc)then
       end if
      
       Phi = SC_Xyz_D(Phi_) 
       Theta = SC_Xyz_D(Theta_)
       rSinTheta = R *sin(Theta)
       
       SC_Xyz_D(x_) = rSinTheta * cos(Phi) 
       SC_Xyz_D(y_) = rSinTheta * sin(Phi) 
       SC_Xyz_D(z_) = R * cos(Theta)
    end if
    

    IH_Xyz_D = matmul(ScToIh_DD, SC_Xyz_D)*&
         Grid_C(SC_)%UnitX/Grid_C(IH_)%UnitX
    IsInterfacePoint=.true.

  end subroutine map_sc_ih
  !=================================================================!

  subroutine IH_get_for_sc_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7
    integer::BuffX_,BuffZ_
    !------------------------------------------------------------

    BuffX_ = nVarIhSc + 1
    BuffZ_ = nVarIhSc + 3

    call IH_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)

    !Transform velocity
    State_V(RhoUx_:RhoUz_)=State_V(Rho_)*&
         transform_velocity(tNow,&
         State_V(RhoUx_:RhoUz_)/State_V(Rho_),&
         State3_V(BuffX_:BuffZ_)/State_V(Rho_),&
         Grid_C(IH_)%TypeCoord,Grid_C(SC_)%TypeCoord)
    
    State_V(Bx_:Bz_)=matmul(IhToSc_DD,State_V(Bx_:Bz_))

  end subroutine IH_get_for_sc_and_transform

!===============================================================!
!===============================================================!
!===============================================================!
!===============================================================
!BOP
!IROUTINE: couple_sc_ih - interpolate and get MHD state at the IH buffer grid 
!INTERFACE:
  subroutine couple_sc_ih(TimeCoupling)
    use ModIoUnit
    !INPUT ARGUMENTS:
    interface
       subroutine SC_get_for_mh(&
            nPartial,iGetStart,Get,w,State_V,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::w
         real,dimension(nVar),intent(out)::State_V
       end subroutine SC_get_for_mh
    end interface
    integer::iPoint,nU_I(2)
    real,intent(in)::TimeCoupling
    integer,save::iCoupling=0
    integer::iFile
    character(LEN=21)::NameFile
    logical::DoneMatchIBC=.false.
!EOP
    if(.not.RouterScBuff%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    call SC_synchronize_refinement(RouterScBuff%iProc0Source,&
                                   RouterScBuff%iCommUnion)


    if(SC_iGridInScIh/=i_realization(SC_))then  
       call set_router(&
            GridDescriptorSource=SC_SourceGrid,&
            GridDescriptorTarget=BuffGD,&
            Router=RouterScBuff,&
            mapping=buffer_grid_point,&
            interpolate=interpolation_fix_reschange)
       SC_iGridInScIh= i_realization(SC_)
    end if

    call couple_buffer_grid(&
         RouterScBuff,&
         nVar = nVarIhSc, &
         fill_buffer=SC_get_for_mh,&
         NameBuffer='IH_from_sc',&
         TargetID_=IH_)

    if(.not.DoneMatchIBC)then
       DoneMatchIBC=.true.
       if(is_proc(IH_))call IH_match_ibc
    end if
    if(DoTest.and.is_proc0(compid_grid(BuffGD%DD%Ptr)))then
       nU_I=ubound_vector('IH_from_sc')
       iCoupling=iCoupling+1
       iFile=io_unit_new()
       write(NameFile,'(a,i4.4,a)')'./IH/from_sc_',iCoupling,'.dat'
       open(iFile,FILE=NameFile,STATUS='unknown')
       do iPoint=1,nU_I(2)
          write(iFile,*)point_state_v('IH_from_sc', nVarIhSc, iPoint)
       end do
       close(iFile)
    end if
    if(DoTest)write(*,*)'Couple passed at PE=',i_proc()
  end subroutine couple_sc_ih
  !======================================================!
  subroutine buffer_grid_point(&
       nDimFrom,Sph_D,nDimTo,SC_Coord_D,IsInterfacePoint)         

    ! Transform from the spherical buffer grid to the Cartesian SC grid

    integer,intent(in)                    :: nDimFrom, nDimTo       
    real,dimension(nDimFrom), intent(in)  :: Sph_D
    real,dimension(nDimTo),   intent(out) :: SC_Coord_D
    logical,intent(out)::IsInterfacePoint

    ! The order of spherical indexes as in BATSRUS
    ! This is a left handed system 
    integer,parameter::r_=1,Psi_=2,Theta_=3,x_=1,y_=2,z_=3
    
    real :: rSinTheta, BuffXyz_D(x_:z_)
    !-----------------------------------------------------------------------
    !In each mapping the corrdinates of the TARGET grid point (Buffer)
    !shoud be be transformed to the SOURCE (SC) generalized coords.
    
    if(.not.IsSphericalSc)then
   
       RSinTheta    = Sph_D(r_)*sin(Sph_D(Theta_))!To be modified
    
       BuffXyz_D(x_) = RSinTheta*cos(Sph_D(Psi_))  
       BuffXyz_D(y_) = RSinTheta*sin(Sph_D(Psi_))  
       BuffXyz_D(z_) = Sph_D(r_)*cos(Sph_D(Theta_))
    
       !\
       ! The buffer grid coordinates are normalized by the unit of length
       ! of IH. Therefore,
       !/
       SC_Coord_D = BuffXyz_D *&
            Grid_C(IH_)%UnitX/Grid_C(SC_)%UnitX
    else
       !\
       ! The buffer grid coordinates are normalized by the unit of length
       ! of IH. Therefore,
       !/
       Sc_Coord_D(R_) =  Sph_D(r_) *&
            Grid_C(IH_)%UnitX/Grid_C(SC_)%UnitX
       
       if(UseLogRSc) then
          Sc_Coord_D(R_) = log(Sc_Coord_D(R_))
       end if
       Sc_Coord_D(Psi_:Theta_) = Sph_D(Psi_:Theta_) 
    
    end if

    IsInterfacePoint=.true.
  end subroutine buffer_grid_point
  
end module CON_couple_ih_sc

