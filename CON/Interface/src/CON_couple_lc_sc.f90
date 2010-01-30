!^CMP COPYRIGHT UM
!^CMP FILE SC
!^CMP FILE LC
!BOP
!MODULE: CON_couple_sc_lc - couple SC and LC both ways
!INTERFACE:
module CON_couple_sc_lc

  !DESCRIPTION:
  ! This coupler uses the SWMF parallel coupling toolkit.
  ! The LC grid is coupled to a buffer grid in SC. The buffer grid
  ! uses the same coordinate system as LC, so the transformation is
  ! done in the SC wrapper.
  !
  ! The SC grid is coupled to the outer ghost cells of the LC grid directly.
  ! Both LC and SC use AMR grids, the buffer is a simple spherical grid.
  
  !USES:
  use CON_coupler
  use CON_axes, ONLY: transform_matrix,transform_velocity
  use ModConst

  implicit none
  private !except
  !
  !PUBLIC MEMBER FUNCTIONS:
  public:: couple_sc_lc_init
  public:: couple_sc_lc,couple_lc_sc

  !REVISION HISTORY:
  ! 7/23/03 Sokolov I.V.<igorsok@umich.edu> - prototype for sc-gm
  ! 7/04/04                                 - version for sc-lc
  ! 7/20/04                                 - version for lc-buffer
  !EOP

  !To trace the possible changes in the grids and/or mapping
  integer :: SC_iGridRealization=-2
  integer :: LC_iGridInLcSc=-2
  integer :: LC_iGridInScLc=-2

  ! SC <-> LC conversion matrices
  real :: ScToLc_DD(3,3),LcToSc_DD(3,3)

  ! Maximum time difference [s] without remap 
  ! The 600 s corresponds to about 0.1 degree rotation between LC and SC
  real :: dTimeMappingMax = 600.0
 
  type(RouterType),save             :: RouterScLc
  type(RouterType),save             :: RouterLcBuff
  type(GridDescriptorType),save     :: LC_SourceGrid
  type(GridDescriptorType),save     :: LC_TargetGrid
  type(GridDescriptorType),save     :: SC_Grid
  type(GridDescriptorType),save     :: BuffGD
  type(DomainDecompositionType),&
       save,target                  :: BuffDD
  logical :: DoInitialize=.true., DoTest, DoTestMe
  real :: tNow
  character(len=*), parameter :: NameMod='couple_sc_lc'
  logical::IsSphericalLc=.false. , UseGenRLc = .false., UseLogRLc = .false.
  integer::iError
  
  !Parameters of the stretched grid, if needed
  integer :: nGenRGridLc
  real    :: DeltaGen
   
contains
  !===============================================================!
  subroutine couple_sc_lc_init
    interface
       subroutine SC_set_buffer_grid(Dd)
         use CON_domain_decomposition
         implicit none
         type(DomainDecompositionType),&
              intent(out)::Dd
       end subroutine SC_set_buffer_grid
    end interface

    if(.not.DoInitialize)return
    DoInitialize=.false.
    
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    IsSphericalLc = index(Grid_C(LC_) % TypeGeometry,'spherical') > 0 
    UseLogRLc     = index(Grid_C(LC_) % TypeGeometry,'lnr'      ) > 0
    UseGenRLc     = index(Grid_C(LC_) % TypeGeometry,'genr'     ) > 0
    if(UseGenRLc) then
       nGenRGridLc = size(Grid_C(LC_) % Coord1_I,1)
       if(nGenRGridLc==1) &
            call CON_stop('Stretched grid in LC is not properly initialized')
       DeltaGen = 1.0/(nGenRGridLC - 1)
    end if

    call init_coupler(              &    
       iCompSource=SC_,             & ! component index for source
       iCompTarget=LC_,             & ! component index for target
       nGhostPointTarget=2,         & ! number of halo points in target
       GridDescriptorSource=SC_Grid,& ! OUT!\
       GridDescriptorTarget=LC_TargetGrid,& !-General coupler variables 
       Router=RouterScLc)             ! OUT!/
    
   
    SC_iGridRealization=-1
    LC_iGridInScLc     =-1
    LC_iGridInLcSc     =-1
    
    call SC_set_buffer_grid(BuffDD)
    call set_standard_grid_descriptor(&
         BuffDD,          &
         Standard_=Nodes_,&
         nGhostGridPoints=1,  &
         GridDescriptor=BuffGD)
    call set_standard_grid_descriptor(&
         LC_,GridDescriptor=LC_SourceGrid)
    call init_buffer_grid_couple(&
         SourceGD=LC_SourceGrid,&
         TargetGD=BuffGD, &
         RouterToBuffer=RouterLcBuff,    &
         nVar=8,&
         NameBuffer='SC_from_lc')  !Version for the first order in time
  end subroutine couple_sc_lc_init
  !===============================================================!
  !BOP
  !IROUTINE: couple_sc_lc - get SC solution at LC outer ghostpoints
  !INTERFACE:
  subroutine couple_sc_lc(TimeCoupling)
    !INPUT ARGUMENTS:
    interface
       subroutine LC_put_from_mh(nPartial,&
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
       end subroutine LC_put_from_mh
    end interface

    real,intent(in)::TimeCoupling
    !EOP

    ! Last coupling time
    real :: TimeCouplingLast = -1.0
    !-------------------------------------------------------------------------

    if(.not.RouterScLc%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    ! Synchronize and broadcast domain decompostion (AMR may have changed it)
    call SC_synchronize_refinement(RouterScLc%iProc0Source,&
                                   RouterScLc%iCommUnion)
    call LC_synchronize_refinement(RouterScLc%iProc0Target,&
                                   RouterScLc%iCommUnion)

    ! Redo the router if any of the grids changed or 
    ! the two coordinate systems have rotated away too much
    if(SC_iGridRealization/=i_realization(SC_).or.&     
         LC_iGridInScLc/=i_realization(LC_) .or. &
         (Grid_C(SC_) % TypeCoord /= Grid_C(LC_) % TypeCoord &
         .and. TimeCoupling - TimeCouplingLast > dTimeMappingMax)) then  
       ! Recalculate the LC to SC transformation matrix used in mapping
       ! a target point (LC) to the source (SC)
       !(it is time dependent in general)

       LcToSc_DD = transform_matrix(TimeCoupling, &
            Grid_C(LC_) % TypeCoord, Grid_C(SC_) % TypeCoord)

       !Recalculate the tramsposed matrix, used in transforming vectors
       !of velocity and magnetic field to be sent from SC to LC
       ScToLc_DD = transpose(LcToSc_DD)

       call set_router(&
            GridDescriptorSource=SC_Grid,&
            GridDescriptorTarget=LC_TargetGrid,&
            Router=RouterScLc,&
            is_interface_block=is_boundary_block,&
            interface_point_coords=outer_cells, &
            mapping=map_lc_sc, &
            interpolate=interpolation_fix_reschange)

       SC_iGridRealization = i_realization(SC_)
       LC_iGridInScLc      = i_realization(LC_)
       TimeCouplingLast    = TimeCoupling
       tNow=TimeCoupling
    end if
    call couple_comp(&
         RouterScLc,&
         nVar=8,&
         fill_buffer=SC_get_for_lc_and_transform,&
         apply_buffer=LC_put_from_mh)

  end subroutine couple_sc_lc
  !======================================================!
  logical function is_boundary_block(lGlobalTreeNode)
    integer,parameter::R_=1
    integer,intent(in)::lGlobalTreeNode
    logical,dimension(3)::IsBoundary_D

    IsBoundary_D=is_right_boundary_d(&
         LC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    !For spherical domain
    is_boundary_block=IsBoundary_D(R_).and.IsSphericalLc
    
    !Now if IsSphericalLc ==.true. then is_boundary_block is .true.
    !only if the right block boundary along the radial coordinate
    !is the LC domain boundary.
    !Else is_boundary_block = .false.

    !For cartesian box   
     IsBoundary_D= IsBoundary_D.or.&
         is_left_boundary_d(&
         LC_TargetGrid%DD%Ptr,lGlobalTreeNode)
     !Now IsBoundary_D(iDim) is true if any of the boundaries along 
     !the direction iDim is the LC boundary

    is_boundary_block=is_boundary_block.or. &
         (any(IsBoundary_D).and.(.not.IsSphericalLc))
    !Now if IsSphericalLc = .true. the value of is_boundary_block does
    !not change. Otherwise is_boundary_block is true if any of the
    !block boundaries is the boundary of the LC domain
   
    
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
         LC_TargetGrid%DD%Ptr,lGlobalTreeNode)
    IsRightFace_D=i_D(x_:z_)>&
         ncells_decomposition_d(LC_TargetGrid%DD%Ptr).and.&
         is_right_boundary_d(LC_TargetGrid%DD%Ptr,lGlobalTreeNode)

    !For spherical grid
    IsInterfacePoint=IsRightFace_D(R_).and.IsSphericalLc

    !For Cartesian grid:

    IsInterfacePoint=IsInterfacePoint.or.&
         (any(IsRightFace_D.or.IsLeftFace_D).and.(.not.IsSphericalLc))
  end subroutine outer_cells 
  !========================================================!
  subroutine map_lc_sc(&
       LC_nDim,LC_XyzIn_D,SC_nDim,SC_Xyz_D,IsInterfacePoint)

    integer,intent(in)::SC_nDim,LC_nDim
    real,dimension(LC_nDim),intent(in)::LC_XyzIn_D
    real,dimension(SC_nDim),intent(out)::SC_Xyz_D
    logical,intent(out)::IsInterfacePoint
    
    real, dimension(LC_nDim) :: LC_Xyz_D
    integer, parameter :: R_ = 1, Phi_=2, Theta_ = 3, x_ = 1, y_ = 2, z_ = 3
    real :: R, Gen, Phi, Theta, rSinTheta
    !--------------------------------------- 
    !In each mapping the corrdinates of the TARGET grid point (LC)
    !shoud be be transformed to the SOURCE (SC) generalized coords.
    if(.not.IsSphericalLc)then
       LC_Xyz_D = LC_XyzIn_D
    else
       !transform to dimensionless cartesian Xyz
       R = LC_Xyz_D(R_)
       if(UseLogRLc)then
          R = exp(R)
       elseif(UseGenRLc)then
       end if
      
       Phi = LC_Xyz_D(Phi_) 
       Theta = LC_Xyz_D(Theta_)
       rSinTheta = R *sin(Theta)
       
       LC_Xyz_D(x_) = rSinTheta * cos(Phi) 
       LC_Xyz_D(y_) = rSinTheta * sin(Phi) 
       LC_Xyz_D(z_) = R * cos(Theta)
    end if
    

    SC_Xyz_D = matmul(LcToSc_DD, LC_Xyz_D)*&
         Grid_C(LC_)%UnitX/Grid_C(SC_)%UnitX
    IsInterfacePoint=.true.

  end subroutine map_lc_sc
  !=================================================================!

  subroutine SC_get_for_lc_and_transform(&
       nPartial,iGetStart,Get,w,State_V,nVar)

    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::w
    real,dimension(nVar),intent(out)::State_V
    real,dimension(nVar+3)::State3_V
    integer, parameter :: Rho_=1, RhoUx_=2, RhoUz_=4, Bx_=5, Bz_=7,&
         BuffX_=9,BuffZ_=11
    !------------------------------------------------------------
    call SC_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)

    !Transform velocity
    State_V(RhoUx_:RhoUz_)=State_V(Rho_)*&
         transform_velocity(tNow,&
         State_V(RhoUx_:RhoUz_)/State_V(Rho_),&
         State3_V(BuffX_:BuffZ_)/State_V(Rho_),&
         Grid_C(SC_)%TypeCoord,Grid_C(LC_)%TypeCoord)
    
    State_V(Bx_:Bz_)=matmul(ScToLc_DD,State_V(Bx_:Bz_))

  end subroutine SC_get_for_lc_and_transform

!===============================================================!
!===============================================================!
!===============================================================!
!===============================================================
!BOP
!IROUTINE: couple_lc_sc - interpolate and get MHD state at the SC buffer grid 
!INTERFACE:
  subroutine couple_lc_sc(TimeCoupling)
    use ModIoUnit
    !INPUT ARGUMENTS:
    interface
       subroutine LC_get_for_mh(&
            nPartial,iGetStart,Get,w,State_V,nVar)
         use CON_router
         implicit none
         integer,intent(in)::nPartial,iGetStart,nVar
         type(IndexPtrType),intent(in)::Get
         type(WeightPtrType),intent(in)::w
         real,dimension(nVar),intent(out)::State_V
       end subroutine LC_get_for_mh
    end interface
    integer::iPoint,nU_I(2)
    real,intent(in)::TimeCoupling
    integer,save::iCoupling=0
    integer::iFile
    character(LEN=21)::NameFile
    logical::DoneMatchIBC=.false.
!EOP
    if(.not.RouterLcBuff%IsProc)return
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    call LC_synchronize_refinement(RouterLcBuff%iProc0Source,&
                                   RouterLcBuff%iCommUnion)


    if(LC_iGridInLcSc/=i_realization(LC_))then  
       call set_router(&
            GridDescriptorSource=LC_SourceGrid,&
            GridDescriptorTarget=BuffGD,&
            Router=RouterLcBuff,&
            mapping=buffer_grid_point,&
            interpolate=interpolation_fix_reschange)
       LC_iGridInLcSc= i_realization(LC_)
    end if

    call couple_buffer_grid(&
         RouterLcBuff,&
         nVar=8,&
         fill_buffer=LC_get_for_mh,&
         NameBuffer='SC_from_lc',&
         TargetID_=SC_)
    if(.not.DoneMatchIBC)then
       DoneMatchIBC=.true.
       if(is_proc(SC_))call SC_match_ibc
    end if
    if(DoTest.and.is_proc0(compid_grid(BuffGD%DD%Ptr)))then
       nU_I=ubound_vector('SC_from_lc')
       iCoupling=iCoupling+1
       iFile=io_unit_new()
       write(NameFile,'(a,i4.4,a)')'./SC/from_lc_',iCoupling,'.dat'
       open(iFile,FILE=NameFile,STATUS='unknown')
       do iPoint=1,nU_I(2)
          write(iFile,*)point_state_v('SC_from_lc',8,iPoint)
       end do
       close(iFile)
    end if
    if(DoTest)write(*,*)'Couple passed at PE=',i_proc()
  end subroutine couple_lc_sc
  !======================================================!
  subroutine buffer_grid_point(&
       nDimFrom,Sph_D,nDimTo,LC_Coord_D,IsInterfacePoint)         

    ! Transform from the spherical buffer grid to the Cartesian LC grid

    integer,intent(in)                    :: nDimFrom, nDimTo       
    real,dimension(nDimFrom), intent(in)  :: Sph_D
    real,dimension(nDimTo),   intent(out) :: LC_Coord_D
    logical,intent(out)::IsInterfacePoint

    ! The order of spherical indexes as in BATSRUS
    ! This is a left handed system 
    integer,parameter::r_=1,Psi_=2,Theta_=3,x_=1,y_=2,z_=3
    
    real :: rSinTheta, BuffXyz_D(x_:z_)
    !-----------------------------------------------------------------------
    !In each mapping the corrdinates of the TARGET grid point (Buffer)
    !shoud be be transformed to the SOURCE (LC) generalized coords.
    
    if(.not.IsSphericalLc)then
   
       RSinTheta    = Sph_D(r_)*sin(Sph_D(Theta_))!To be modified
    
       BuffXyz_D(x_) = RSinTheta*cos(Sph_D(Psi_))  
       BuffXyz_D(y_) = RSinTheta*sin(Sph_D(Psi_))  
       BuffXyz_D(z_) = Sph_D(r_)*cos(Sph_D(Theta_))
    
       !\
       ! The buffer grid coordinates are normalized by the unit of length
       ! of SC. Therefore,
       !/
       LC_Coord_D = BuffXyz_D *&
            Grid_C(SC_)%UnitX/Grid_C(LC_)%UnitX
    else
       !\
       ! The buffer grid coordinates are normalized by the unit of length
       ! of SC. Therefore,
       !/
       Lc_Coord_D(R_) =  Sph_D(r_) *&
            Grid_C(SC_)%UnitX/Grid_C(LC_)%UnitX
       
       if(UseLogRLc) then
          Lc_Coord_D(R_) = log(Lc_Coord_D(R_))
       end if
       Lc_Coord_D(Psi_:Theta_) = Sph_D(Psi_:Theta_) 
    
    end if

    IsInterfacePoint=.true.
  end subroutine buffer_grid_point
  
end module CON_couple_sc_lc

