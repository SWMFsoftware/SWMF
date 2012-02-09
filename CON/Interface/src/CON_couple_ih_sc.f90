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

  ! Variables used by global MPI coupler
  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe=.true., IsInitialized=.false., DoMatchIBC = .true.

  ! Size and limits of the 3D spherical buffer grid
  integer, save :: iSize, jSize, kSize, nCell_D(3)
  real, save    :: BufferMinMaxSc_DI(3,2), BufferMinMaxIh_DI(3,2)


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
    !--------------------------------------------------------------------------
    ! REDIRECT to couple_ih_sc_init_global if needed
    if (UseGlobalMpiCoupler_CC(SC_,IH_) .or. UseGlobalMpiCoupler_CC(IH_,SC_) .or. &
         index(Grid_C(SC_) % TypeGeometry,'spherical') > 0)then
       UseGlobalMpiCoupler_CC(SC_,IH_) = .TRUE.
       UseGlobalMpiCoupler_CC(IH_,SC_) = .TRUE.

       call couple_ih_sc_init_global
       RETURN
    end if


    if(.not.DoInitialize)return
    DoInitialize=.false.
    
    call CON_set_do_test(NameMod,DoTest,DoTestMe)

    call set_couple_var_info(SC_,IH_)

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
         nVar = nVarCouple, &
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
    ! REDIRECT to couple_ih_sc_global if needed
    if(UseGlobalMpiCoupler) then
       call couple_ih_sc_global(TimeCoupling)
       RETURN
    end if

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
         nVar = nVarCouple, &
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

    integer,intent(in)              :: nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)   :: Get
    type(WeightPtrType),intent(in)  :: w
    real,dimension(nVar),intent(out):: State_V
    real,dimension(nVar+3)          :: State3_V

    ! variable indices in buffer
    integer   ::  &
         iRhoCouple,   &
         iRhoUxCouple, &
         iRhoUzCouple, &
         iBxCouple,    &
         iBzCouple

    integer :: BuffX_,BuffZ_
    !------------------------------------------------------------
    ! get variable indices in buffer
    iRhoCouple       = iVar_V(RhoCouple_)
    iRhoUxCouple     = iVar_V(RhoUxCouple_)
    iRhoUzCouple     = iVar_V(RhoUzCouple_)
    iBxCouple        = iVar_V(BxCouple_)
    iBzCouple        = iVar_V(BzCouple_)

    BuffX_ = nVarCouple + 1
    BuffZ_ = nVarCouple + 3

    call IH_get_for_mh_with_xyz(&
       nPartial,iGetStart,Get,w,State3_V,nVar+3)
    State_V=State3_V(1:nVar)

    !Transform velocity
    ! NOTE: This transformation is only valid for a single fluid
    State_V(iRhoUxCouple:iRhoUzCouple)=State_V(iRhoCouple)*&
         transform_velocity(tNow,&
         State_V(iRhoUxCouple:iRhoUzCouple)/State_V(iRhoCouple),&
         State3_V(BuffX_:BuffZ_)/State_V(iRhoCouple),&
         Grid_C(IH_)%TypeCoord,Grid_C(SC_)%TypeCoord)
    
    ! Transform magnetic field
    State_V(iBxCouple:iBzCouple) = &
         matmul(IhToSc_DD,State_V(iBxCouple:iBzCouple))

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
    ! --------------------------------------------------------
    ! REDIRECT to couple_sc_ih_global if needed
    if(UseGlobalMpiCoupler) then
       call couple_sc_ih_global(TimeCoupling)
       RETURN
    end if

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
         nVar = nVarCouple, &
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
          write(iFile,*)point_state_v('IH_from_sc', nVarCouple, iPoint)
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

  ! =========================================================================
  !                   SUBROUTINES FOR GLOBAL MPI COUPLER
  ! =========================================================================
  !DESCRIPTION:
  
  ! Couple SC and IH components via a global buffer grid
  ! The subroutines:
  !                CON_couple_ih_sc_init
  !                CON_couple_sc_ih
  !                CON_couple_ih_sc
  
  ! redirect to the corresponding subroutines below, with the suffix "_global"
  ! in case UseGlobalMpiCoupler =.TRUE.
  !

  !REVISION HISTORY:
  ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version as external
  !                                       subroutines
  ! 08/27/2003 G.Toth - combined into a module
  ! 12/12/2011 R.Oran <oran@umich.edu> version for two BATSRUS components using
  ! a global buffer grid

  ! =======================================================================
  !IROUTINE: couple_ih_sc_init_global - initialize SC-IH couplings
  !INTERFACE:
  subroutine couple_ih_sc_init_global
 

    logical :: DoTest, DoTestMe
    integer :: iCommWorld, iError, iStatus_I(MPI_STATUS_SIZE)
    real    :: IHToScUnitX
    character(len=*), parameter :: NameSub='couple_ih_sc_init_global'

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    ! Share buffer grid info (set in IH) with SC.
    ! Set buffer grid in IH.
    ! Calculate union communicator

    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(DoTest) write(*,*) NameSub, ' started'

    UseGlobalMpiCoupler = .true.
    if(i_proc() == 0) write(*,*) 'Using global MPI coupler'
    ! Determine which state variables should be coupled
    call set_couple_var_info(SC_,IH_)

    ! Set UseMe to .true. for the participating PE-s
    UseMe = is_proc(SC_) .or. is_proc(IH_)

    ! Set buffer grid location and size in IH, and retrieve them for coupler
    if(is_proc(IH_)) then
       call IH_set_buffer_grid_get_info(IH_,iSize, jSize, kSize, BufferMinMaxIh_DI)
       if(DoTest .and. is_proc0(IH_)) &
            write(*,*) 'In IH, rBuffMinMax: ',BufferMinMaxIh_DI(1,:)

       ! Convert units for radial coordinate  before passing to SC
       BufferMinMaxSc_DI = BufferMinMaxIh_DI
       IhToScUnitX = (Grid_C(IH_)%UnitX/Grid_C(SC_)%UnitX)
       BufferMinMaxSc_DI(1,:) = BufferMinMaxIh_DI(1,:)*IhToScUnitX

       ! Package info for passing via MPI
       nCell_D = (/iSize,jSize,kSize/)
    end if
 
    ! Share with SC head node.
    if(i_proc0(IH_) /= i_proc0(SC_))then
       ! MPI share buffer grid size
       if(is_proc0(IH_)) then
          call MPI_send(nCell_D, 3, MPI_INTEGER, i_proc0(SC_),&
               1, iCommWorld, iError)
      end if
       if(is_proc0(SC_)) then
          call MPI_recv(nCell_D, 3, MPI_INTEGER, i_proc0(IH_),&
               1, iCommWorld, iStatus_I, iError)
       end if

       ! MPI share buffer grid limits
       if(is_proc0(IH_)) &
            call MPI_send(BufferMinMaxSc_DI, 6, MPI_REAL, i_proc0(SC_),&
            1, iCommWorld, iError)
       if(is_proc0(SC_))&
            call MPI_recv(BufferMinMaxSc_DI, 6, MPI_REAL, i_proc0(IH_),&
            1, iCommWorld, iStatus_I, iError)
    end if

    ! Broadcast to all SC nodes.
    if(n_proc(SC_)>1 .and. is_proc(SC_)) then
         call MPI_bcast(nCell_D, 3, MPI_INTEGER, 0, i_comm(SC_), iError)
         call MPI_bcast(BufferMinMaxSc_DI, 6, MPI_REAL, 0, i_comm(SC_), iError)
      end if

    if(is_proc0(SC_) .and. DoTest) &
         write(*,*) 'In SC, rBuffer Min/Max: ',BufferMinMaxSc_DI(1,:)

  end subroutine couple_ih_sc_init_global

  !BOP =======================================================================
  !IROUTINE: couple_sc_ih_global - couple SC component to IH component

  !INTERFACE:
  subroutine couple_sc_ih_global(TimeCoupling)
    
    use ModMpi,    ONLY: MPI_reduce

    !INPUT ARGUMENTS:
    real, intent(in) :: TimeCoupling     ! simulation time at coupling


    !DESCRIPTION:
    ! Couple between two components:
    !    Solar Corona      (SC)  source
    !    Inner Heliosphere (IH)  target
    !
    ! The SC component sends the state variables to a buffer grid.
    ! IH uses the buffer grid to calculate the inner boundary conditions.

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    ! Array to store state vector on all buffer grid points
    real, dimension(:,:,:,:), allocatable :: Buffer_VIII, BufferGlobal_VIII

    ! MPI related variables
    integer :: iError, nSize, iStatus_I(MPI_STATUS_SIZE)

    ! Variable for output file
    integer :: i,j,k, iFile
    integer, save :: iCoupling = 0

    character (len=*), parameter :: NameSub='couple_sc_ih_global'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()
    ! Exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, SCi_iProc0, IHi_iProc0=', &
         iProcWorld,i_proc0(SC_),i_proc0(IH_)

    ! Allocate and intialize local buffer on SC processors, including ghost cells
    if(is_proc(SC_)) then
       allocate(Buffer_VIII(nVarCouple,iSize,0:jSize+1,0:kSize+1), stat=iError)
       Buffer_VIII = 0.0
    end if

    ! Allocate global buffer on all processors
   allocate(BufferGlobal_VIII(nVarCouple,iSize,0:jSize+1,0:kSize+1), stat=iError)
    BufferGlobal_VIII = 0.0

    nSize = iSize*(jSize+2)*(kSize+2)*nVarCouple

    if(DoTest)write(*,*)NameSub,', finished allocating global buffer'
    ! Fill in the coupled state variables
    if(is_proc(SC_)) then
       call SC_get_for_global_buffer(iSize,jSize,kSize, &
            BufferMinMaxIh_DI,TimeCoupling,SC_, IH_, Buffer_VIII)
       ! Collect to the SC root PE
       call MPI_reduce(Buffer_VIII, BufferGlobal_VIII, nSize, MPI_REAL,MPI_SUM,&
            0, i_comm(SC_), iError)
    end if

    if(DoTest)write(*,*)NameSub,', SC finished filling global buffer'
    if(i_proc0(SC_) /= i_proc0(IH_))then
       ! Pass filled buffer to IH root PE
       if(is_proc0(SC_)) &
            call MPI_send(BufferGlobal_VIII,nSize,MPI_REAL,i_Proc0(IH_),&
            1,i_comm(),iError)
       if(is_proc0(IH_)) &
            call MPI_recv(BufferGlobal_VIII,nSize,MPI_REAL,i_Proc0(SC_),&
            1,i_comm(),iStatus_I,iError)
    end if

    ! Broadcast buffer to all IH processors
    if(n_proc(IH_)>1 .and. is_proc(IH_)) &
         call MPI_bcast(BufferGlobal_VIII,nSize,MPI_REAL,0,i_comm(IH_),iError)

    if(DoTest)write(*,*)NameSub,', variables transferred',&
         ', iProc:',iProcWorld

    !\
    ! Save buffer in IH for later use
    !/
    if(is_proc(IH_))then
       call IH_save_global_buffer(nVarCouple, iSize, jSize, kSize,BufferGlobal_VIII)
       if(DoTest) &
            write(*,*)NameSub//' iProc, Buffer(:,1,1,1)=',&
            iProcWorld,Buffer_VIII(:,1,1,1)
    end if

    if(DoMatchIBC) then
       DoMatchIBC = .false.
       if(is_proc(IH_))call IH_match_ibc
    end if
    !\
    ! Deallocate buffer to save memory
    !/
    if(is_proc(SC_)) deallocate(Buffer_VIII)
    deallocate(BufferGlobal_VIII)

    if(DoTest)write(*,*)NameSub,', variables deallocated',&
         ', iProc:',iProcWorld

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_sc_ih_global

  !BOP =======================================================================
  !IROUTINE: couple_ih_sc_global - couple IH to SC
  !INTERFACE:
  subroutine couple_ih_sc_global(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    !DESCRIPTION:
    ! Couple between two components:
    !    Inner Heliosphere (IH) source
    !    Solar Corona      (SC) target
    !
    ! Send state variable from IH to outer cells in SC.
    !EOP

    logical :: DoTest, DoTestMe
    integer :: iProcWorld

    ! Buffer for state variable to fill outer cells of SC

    real, dimension(:,:,:,:), allocatable :: Buffer_VIII

    ! MPI related variables
    integer :: iError, nSize, iStatus_I(MPI_STATUS_SIZE)

    integer :: iBlock, iProcTo
    character (len=*), parameter :: NameSub='couple_ih_sc_global'
    !-------------------------------------------------------------------------
    call CON_stop(NameSub//' is not yet implemented. Correct #COUPLERTYPE command.')
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iProcWorld = i_proc()

    if(DoTest)write(*,*)NameSub,' starting iProc=',iProcWorld

    ! Exclude PEs which are not involved
    if(.not.UseMe) RETURN

    if(DoTest)write(*,*)NameSub,' starting, iProc=',iProcWorld
    if(DoTest)write(*,*)NameSub,', iProc, IHi_iProc0, i_proc0(SC_)=', &
         iProcWorld,i_proc0(IH_),i_proc0(SC_)

    !\
    ! Allocate buffers for the variables both in SC and IH
    !/
    allocate(Buffer_VIII(nVarCouple,iSize, jSize, kSize), stat=iError)
    !\
    ! Get IH state variable for SC
    !/
    if(is_proc(IH_)) &
         call IH_get_for_mh(Buffer_VIII, iSize, jSize, kSize, nVarCouple)

    !\                                    
    ! Transfer variables from IH to SC
    !/
    if(i_proc0(SC_) /= i_proc0(IH_))then
       nSize = iSize*jSize*kSize*nVarCouple
       if(is_proc0(IH_)) &
            call MPI_send(Buffer_VIII, nSize, MPI_REAL, i_proc0(SC_),&
            1, i_comm(), iError)
       if(is_proc0(SC_)) &
            call MPI_recv(Buffer_VIII, nSize, MPI_REAL, i_proc0(IH_),&
            1,i_comm(), iStatus_I, iError)
    end if
    !\
    ! Put variables into IH
    !/
    if(is_proc(SC_)) &
         !call SC_put_from_mh(Buffer_VIII, iSize+1, jSize, kSize, nVarCouple)

    !\
    ! Deallocate buffer to save memory
    !/
    deallocate(Buffer_VIII)

    if(DoTest)write(*,*)NameSub,' finished, iProc=',iProcWorld

  end subroutine couple_ih_sc_global

  
end module CON_couple_ih_sc

