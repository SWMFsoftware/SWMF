!===============================================================!
!^CFG COPYRIGHT UM
!BOP
!MODULE: CON_coupler - methods for general coupler variables
!DESCRIPTION:
! This class can be used in the couplers and wrappers.
! It should contain everything needed for the coupling.
! Provides simplified methods to access Igor Sokolov's
! general coupler library.
!INTERFACE:
module CON_coupler

  !USES:
  use CON_comp_param
  use CON_world
  use CON_buffer_grid
  use CON_global_message_pass, &
       couple_comp => global_message_pass
  use CON_router
  use CON_time, ONLY: FreqType
  use CON_test_global_message_pass
  use ModUtilities, ONLY: check_allocate
  implicit none

  !PUBLIC TYPES:
  public :: CoordSystemType ! Extend grid descriptor type with coordinate info
  integer, parameter :: lTypeCoord=3

  integer, parameter :: lTypeGeometry = 15

  type CoordSystemType    
     real, dimension(:),pointer    :: Coord1_I, Coord2_I, Coord3_I
     integer, dimension(3)         :: nCoord_D
     character (len=lTypeCoord)    :: TypeCoord
     character (len=lTypeGeometry) :: TypeGeometry
     real::UnitX
  end type CoordSystemType

  !PUBLIC DATA MEMBERS:
  public :: Grid_C          ! Store the coordinate information for components

  type(CoordSystemType), target, save :: Grid_C(MaxComp+3)

  public :: Couple_CC           ! Frequency of couplings

  type(FreqType) :: Couple_CC(MaxComp,MaxComp) = &
       FreqType(.false.,-1,-1.0,0,0.0)

  public :: MaxCouple           ! Maximum number of couplings

  integer, parameter :: MaxCouple = 23

  public :: nCouple             ! Actual number of couplings

  integer :: nCouple = MaxCouple

  public :: iCompCoupleOrder_II ! The order of couplings

  integer :: iCompCoupleOrder_II(2,MaxCouple) = reshape ( (/&
       ! This is the default order based on the propagation of information
       ! from component to component
       SC_, IH_, & 
       IH_, SC_, & 
       SC_, SP_, & ! The order of these two couplings is mandatory
       IH_, SP_, & ! Do not modify them, please.
       IH_, GM_, &
       GM_, IE_, &
       GM_, IM_, &
       GM_, PW_, &
       GM_, RB_, &
       UA_, IE_, &
       UA_, LA_, &
       LA_, UA_, &
       IE_, IM_, &
       IM_, GM_, &
       IM_, IE_, &
       PW_, GM_, &
       IE_, UA_, &
       IE_, GM_, &
       IE_, PW_, &
       IE_, PS_, &
       IE_, RB_, &
       IH_, OH_, &
       OH_, IH_  &
       /), (/2, MaxCouple/) )

  public :: DoCoupleOnTime_C ! Should do coupling limit the time step

  logical :: DoCoupleOnTime_C(MaxComp) = .true.

  !PUBLIC MEMBER FUNCTIONS:
  public :: set_coord_system      !Sets coordinate information for a component
  public :: gen_to_stretched      !\ Transform generalized coordinates to
  public :: stretched_to_gen      !  stretched ones and vice versa for 
                                  !/  non-uniform structured grid

  public :: init_coupler          !Initializes grids and router
  public :: set_router            !Presets the router indexes
  public :: couple_comp           !Couple two components via global message pass

  public :: check_couple_symm     !Check if coupling is symmetric

  !REVISION HISTORY:
  ! 07/22/03 Gabor Toth <gtoth@umich.edu> - initial prototype 
  ! 08/22/03 I.Sokolov  <igorsok@umich.edu> - modifications of
  !              the procedure interfaces for new coupling toolkit
  !              and init_coupler is added.
  ! 09/02/03 G.Toth - improved Protex description. 
  !                   Use renaming for couple_comp and set_router
  ! 09/07/03 I.Sokolov - bug fixes in init_coord_system_all and
  !                 in set_coord_system. Add gen_to_stretched and 
  !                 streched_to_gen
  ! 08/10/04 I.Sokolov - to avoid an inproper use of init_coord_system_all
  !                      which destroys the Grid_C structure. 
  !EOP

  character(len=*), parameter, private :: NameMod='CON_coupler'
contains
  !===============================================================!
  subroutine set_coord_system( &
       GridID_,       &! Grid ID
       TypeCoord,     &! Coordinate system type (MAG,GEO,..)
       UnitX,         &! Unit of length, in SI (meters)
       TypeGeometry,  &! Geometry type (cartesian, spherical_lnr, etc) 
       Coord1_I,      &! Non-uniform coords in 1st dim (optional)
       Coord2_I,      &! Non-uniform coords in 2nd dim (optional)
       Coord3_I,      &! Non-uniform coords in 3rd dim (optional)
       iProc0In,      &
       iCommIn) 

    character(len=*), parameter :: NameSub='set_coord_system'

    integer, intent(in) :: GridID_
    character(len=lTypeCoord), intent(in)  :: TypeCoord

    ! Optional parameters

    character(len=lTypeGeometry), optional, intent(in)  :: &
         TypeGeometry

    real,    intent(in), optional :: &
         Coord1_I(:), Coord2_I(:), Coord3_I(:)
    integer,intent(in),optional:: iProc0In,iCommIn
    real,intent(in),optional::UnitX

    character(LEN=lTypeGeometry), parameter :: &
         TypeGeometryBlanck = '               '  ! 15 spaces

    integer :: iProc0, iComm, iError
    logical :: IsRoot
    type(CoordSystemType), pointer :: ThisGrid
    logical,save::DoInit=.true.
    !-------------------------------------------------------------!
    if(DoInit)call init_coord_system_all
    if(present(iProc0In).and.present(iCommIn))then
       iProc0=iProc0In
       iComm=iCommIn
    else
       iProc0 = i_proc0(compid_grid(GridID_))
       iComm=i_comm()
    end if

    ThisGrid => Grid_C(GridID_)

    IsRoot=is_proc0(compid_grid(GridID_))

    if(IsRoot)Thisgrid%TypeCoord=TypeCoord

    ! Broadcast the coordinate type
    call MPI_bcast(ThisGrid%TypeCoord,lTypeCoord,MPI_CHARACTER,&
         iProc0,iComm,iError)
    if(present(UnitX))then
       if(IsRoot)Thisgrid%UnitX=UnitX
       ! Broadcast the unit of length
       call MPI_bcast(ThisGrid%UnitX,1,MPI_REAL,&
            iProc0,iComm,iError)
    end if

    if(IsRoot)then
       if(present(TypeGeometry))then
          Thisgrid%TypeGeometry = TypeGeometryBlanck
          Thisgrid%TypeGeometry(1:len_trim(TypeGeometry)) = &
               trim(TypeGeometry)
       else
          Thisgrid%TypeGeometry = 'cartesian      '
       end if
    end if

    ! Broadcast the geometry type
    call MPI_bcast(ThisGrid%TypeGeometry,lTypeGeometry,MPI_CHARACTER,&
         iProc0,iComm,iError)

    ! Get the size of the coordinate arrays
    if(IsRoot)then
       ThisGrid%nCoord_D = 0
       if(present(Coord1_I)) ThisGrid%nCoord_D(1)=size(Coord1_I)
       if(present(Coord2_I)) ThisGrid%nCoord_D(2)=size(Coord2_I)
       if(present(Coord3_I)) ThisGrid%nCoord_D(3)=size(Coord3_I)
    end if

    call MPI_bcast(ThisGrid%nCoord_D,3,MPI_INTEGER,iProc0,&
         iComm,iError)

    if(ThisGrid%nCoord_D(1)>0)then 
       if(associated(ThisGrid%Coord1_I))deallocate(ThisGrid%Coord1_I)
       allocate(ThisGrid%Coord1_I(&
            ThisGrid%nCoord_D(1)),stat=iError)
       call check_allocate(iError,NameSub//'Coord1_I')
       if(IsRoot)then
          ThisGrid%Coord1_I = Coord1_I
       end if
       call MPI_bcast(ThisGrid%Coord1_I,&
            ThisGrid%nCoord_D(1),MPI_REAL,&
            iProc0,iComm,iError)
    end if

    if(ThisGrid%nCoord_D(2)>0) then
       if(associated(ThisGrid%Coord2_I))deallocate(ThisGrid%Coord2_I)
       allocate(ThisGrid%Coord2_I(&
            ThisGrid%nCoord_D(2)),stat=iError)
       call check_allocate(iError,NameSub//'Coord2_I')
       if(IsRoot)then
          ThisGrid%Coord2_I = Coord2_I
       end if
       call MPI_bcast(ThisGrid%Coord2_I,ThisGrid%nCoord_D(2),&
            MPI_REAL,iProc0,iComm,iError)
    end if

    if(ThisGrid%nCoord_D(3)>0)then
       if(associated(ThisGrid%Coord3_I))deallocate(ThisGrid%Coord3_I)
       allocate(ThisGrid%Coord3_I(&
            ThisGrid%nCoord_D(3)),stat=iError)
       call check_allocate(iError,NameSub//'Coord3_I')
       if(IsRoot)then
          ThisGrid%Coord3_I = Coord3_I
       end if
       call MPI_bcast(ThisGrid%Coord3_I,ThisGrid%nCoord_D(3),&
            MPI_REAL,iProc0,iComm,iError)
    end if
  contains
    subroutine init_coord_system_all
      integer :: iComp
      DoInit=.false.
      do iComp=1,MaxComp+3
         nullify(&
              Grid_C(iComp)%Coord1_I, &
              Grid_C(iComp)%Coord2_I, &
              Grid_C(iComp)%Coord3_I)
         Grid_C(iComp)%nCoord_D = 0
         Grid_C(iComp)%UnitX=cOne
      end do
    end subroutine init_coord_system_all
  end subroutine set_coord_system
  !=============================================================!
  !BOP
  !IROUTINE: gen_to_stretched - transform coords for structured non-uniform grid 
  !INTERFACE:
  subroutine gen_to_stretched(XyzGen_D,XyzStretched_D,nDim,GridID_,DoExtrapolate)
    !INPUT ARGUMENTS:
    integer,intent(in)::nDim,GridID_
    real,dimension(nDim),intent(in)::XyzGen_D

    !Note, that the PRESENCE of this parameter means to do extrapolation means 
    !to do extrapolation, while the VALUE of it, if present, is meaningless
    logical, OPTIONAL, intent(in):: DoExtrapolate 
    !OUTPUT ARGUMENTS:
    real,dimension(nDim),intent(out)::XyzStretched_D
    !DESCRIPTION:
    !   Trasforms generalized coordinates (which for the stretched grids are usually
    !   nothing but the grid point index) to streched coordinates
    !EOP
    real:: OneIfExtrapolate = 1.0
    !------------------------------!
    if(present(DoExtrapolate))then
       OneIfExtrapolate = 1.0
    else
       OneIfExtrapolate = 0.0
    end if

    XyzStretched_D=XyzGen_D
    if(associated(Grid_C(GridID_)%Coord1_I))&
         call stretch(1,Grid_C(GridID_)%Coord1_I)
    if(associated(Grid_C(GridID_)%Coord2_I).and.nDim>1)&
         call stretch(2,Grid_C(GridID_)%Coord2_I)
    if(associated(Grid_C(GridID_)%Coord3_I).and.nDim>2)&
         call stretch(3,Grid_C(GridID_)%Coord3_I)
  contains
    subroutine stretch(iDim,Coord_I)
      integer,intent(in)::iDim
      real,dimension(:),intent(in)::Coord_I
      integer::iL,iU,Number
      real::Fraction
      !--------------!

      iL=lbound(Coord_I,1)
      iU=ubound(Coord_I,1)

      if(XyzStretched_D(iDim)<=real(iL))then
         XyzStretched_D(iDim)=Coord_I(iL) &
              + OneIfExtrapolate * (XyzStretched_D(iDim) - real(iL)) *&
              (Coord_I(iL+1) - Coord_I(iL))
      elseif(XyzStretched_D(iDim)>=real(iU))then
         XyzStretched_D(iDim)=Coord_I(iU)  &
              + OneIfExtrapolate * (XyzStretched_D(iDim) - real(iU)) *&
              (Coord_I(iU) - Coord_I(iU-1))
      else
         Number=floor(XyzStretched_D(iDim))
         Fraction=XyzStretched_D(iDim)-real(Number)
         XyzStretched_D(iDim)=Coord_I(Number)*(cOne-Fraction)+&
              Coord_I(Number+1)*Fraction
      end if
    end subroutine stretch
  end subroutine gen_to_stretched
  !=============================================================!
  !BOP
  !IROUTINE: stretched_to_gen - transform coords for structured non-uniform grid 
  !INTERFACE:
  subroutine stretched_to_gen(XyzStretched_D,XyzGen_D,nDim,GridID_,DoExtrapolate)
    !INPUT ARGUMENTS:
    integer,intent(in)::nDim,GridID_
    real,dimension(nDim),intent(in)::XyzStretched_D

    !Note, that the PRESENCE of this parameter means to do extrapolation means 
    !to do extrapolation, while the VALUE of it, if present, is meaningless
    logical, OPTIONAL, intent(in):: DoExtrapolate 
    !OUTPUT ARGUMENTS:
    real,dimension(nDim),intent(out)::XyzGen_D
    !DESCRIPTION:
    !   Trasforms stretched coordinates  to  generalized coordinates 
    !   (which for the stretched grids are usually
    !   nothing but the grid point index)
    !EOP
    
    real:: OneIfExtrapolate = 1.0
    !------------------------------!
    if(present(DoExtrapolate))then
       OneIfExtrapolate = 1.0
    else
       OneIfExtrapolate = 0.0
    end if

    XyzGen_D=XyzStretched_D
    if(associated(Grid_C(GridID_)%Coord1_I))&
         call gen(1,Grid_C(GridID_)%Coord1_I)
    if(associated(Grid_C(GridID_)%Coord2_I).and.nDim>1)&
         call gen(2,Grid_C(GridID_)%Coord2_I)
    if(associated(Grid_C(GridID_)%Coord3_I).and.nDim>2)&
         call gen(3,Grid_C(GridID_)%Coord3_I)
  contains
    subroutine gen(iDim,Coord_I)
      integer,intent(in)::iDim
      real,dimension(:),intent(in)::Coord_I
      integer::iL,iU,Number
      iL=lbound(Coord_I,1)
      iU=ubound(Coord_I,1)
      if(XyzGen_D(iDim)<=Coord_I(iL))then
         XyzGen_D(iDim)=real(iL) +&
              OneIfExtrapolate * (XyzGen_D(iDim) - Coord_I(iL))/&
              (Coord_I(iL+1) - Coord_I(iL))
      elseif(XyzGen_D(iDim)>=Coord_I(iU))then
         XyzGen_D(iDim)=real(iU) +&
              OneIfExtrapolate * (XyzGen_D(iDim) - Coord_I(iU))/&
              (Coord_I(iU) - Coord_I(iU-1))
      else
         Number=maxloc(Coord_I,1,MASK=Coord_I<=XyzGen_D(iDim))
         XyzGen_D(iDim)=real(Number)+(XyzGen_D(iDim)-Coord_I(Number))/&
              (Coord_I(Number+1)-Coord_I(Number))
      end if
    end subroutine gen
  end subroutine stretched_to_gen
  !=================================================================
  !=============================================================!
  ! simplified interfaces for Igor's coupler
  !BOP
  !INTERFACE:       
  subroutine set_grid_descriptor( &
       iComp,        & ! Component ID
       nDim,         & ! Dimensionality of the grid
       TypeCoord,    & ! Coordinate system type (GSM,GSE,..)
       Coord1_I,     & ! Non-uniform coordinates in 1st dim (op)
       Coord2_I,     & ! Non-uniform coordinates in 2nd dim (op)
       Coord3_I,     & ! Non-uniform coordinates in 3rd dim (op)
       nRootBlock_D, & ! Number of root blocks in all dimensions
       XyzMin_D,     & ! Minimum generalized coordinates
       XyzMax_D,     & ! Maximum generalized coordinates
       nCell_D,      & ! Number of cells in a block in all dims
       iProc_A,      & ! Processor index for all the blocks (op)
       iBlock_A,     & ! Block index for all the blocks (op)
       IsPeriodic_D)   ! Periodicity for all dimesnsions (op)
    !"op"=optional
    !INPUT ARGUMENTS:
    integer, intent(in) :: iComp, nDim
    integer, intent(in) :: nRootBlock_D(nDim), nCell_D(nDim)
    real,    intent(in) :: XyzMin_D(nDim), XyzMax_D(nDim)
    character(len=lTypeCoord), intent(in)  :: TypeCoord

    ! Optional parameters
    real,    intent(in), optional :: Coord1_I(:)
    real,    intent(in), optional :: Coord2_I(:)
    real,    intent(in), optional :: Coord3_I(:)
    integer, intent(in), optional :: iProc_A(:), iBlock_A(:)
    logical, intent(in), optional :: IsPeriodic_D(:)
    !DESCRIPTION: 
    ! Describe and broadcast non-octree grids
    !EOP
    !-----------------------------------------------------------!
    integer:: GridID
    logical:: DoTest,DoTestMe
    !-----------------------------------------------------------!
    call CON_set_do_test('test_grids',DoTest,DoTestMe)
    GridID=iComp
    if(done_dd_init(iComp))return
    call init_decomposition(&
         GridID,                             &!Decomposition ID_
         iComp,                              &! component index
         nDim)                                ! dimensionality
    call set_coord_system(&
         GridID,                             &!Decomposition ID_
         TypeCoord,                          &
         Coord1_I=Coord1_I,                  &
         Coord2_I=Coord2_I,                  &
         Coord3_I=Coord3_I)                           
    if(is_proc0(iComp))call get_root_decomposition(&
         GridID,                             &!Decomposition ID_
         nRootBlock_D ,                      &
         XyzMin_D,                           &
         XyzMax_D,                           &
         nCell_D,                            &
         iProc_A,                            &
         iBlock_A,                           &
         IsPeriodic_D)
    call bcast_decomposition(GridID)
    if(DoTest)&
         call test_global_message_pass(GridID)
  end subroutine set_grid_descriptor


  !=======================================================================
  !IROUTINE: init_coupler - initializes coupler between two components
  subroutine init_coupler(  &    
       iCompSource,         &! component index for source
       nGhostPointSource,   &! number of halo points in Source 
       StandardSource_,     &! CellCentered_ or Nodes_
       nIndexSource,        &! number of indexes for source grid
       iCompTarget,         &! component index for target
       nGhostPointTarget,   &! number of halo points in target 
       StandardTarget_,     &! CellCentered_ or Nodes_
       nIndexTarget,        &! number of indexes for target grid
       GridDescriptorSource,& ! OUT!\
       GridDescriptorTarget,& ! OUT!-General coupler variables 
       Router)                ! OUT!/

    !INPUT ARGUMENTS:
    integer,intent(in)          :: iCompSource, iCompTarget
    integer,intent(in),optional :: nIndexSource, nIndexTarget
    integer,intent(in),optional :: &
         StandardSource_, StandardTarget_
    integer,intent(in),optional :: &
         nGhostPointSource, nGhostPointTarget

    !OUTPUT ARGUMENTS:
    type(GridDescriptorType),intent(out)::GridDescriptorSource
    type(GridDescriptorType),intent(out)::GridDescriptorTarget

    !EOP
    type(RouterType)::Router
    integer::StandardSourceHere_
    integer::StandardTargetHere_
    integer::nGhostPointSourceHere
    integer::nGhostPointTargetHere
    integer::nIndexSourceHere, nIndexTargetHere

    StandardSourceHere_=CellCentered_
    StandardTargetHere_=CellCentered_
    nGhostPointSourceHere=0
    nGhostPointTargetHere=0   
    if(present(nGhostPointSource))&
         nGhostPointSourceHere=nGhostPointSource
    if(present(StandardSource_))&
         StandardSourceHere_=StandardSource_
    call set_standard_grid_descriptor(iCompSource,        &
         nGhostPointSource,  &
         StandardSourceHere_,&
         GridDescriptorSource)

    if(present(nGhostPointTarget))&
         nGhostPointTargetHere=nGhostPointTarget
    if(present(StandardTarget_))&
         StandardTargetHere_=StandardTarget_
    call set_standard_grid_descriptor(iCompTarget,        &
         nGhostPointTarget,  &
         StandardTargetHere_,&
         GridDescriptorTarget)

    if(present(nIndexSource))then
       nIndexSourceHere=nIndexSource
    else
       nIndexSourceHere=ndim_grid(iCompSource)+1
    end if

    if(present(nIndexTarget))then
       nIndexTargetHere=nIndexTarget
    else
       nIndexTargetHere=ndim_grid(iCompTarget)+1
    end if
    call init_router(&
         GridDescriptorSource,&
         GridDescriptorTarget,&
         Router,&
         nIndexSource,&
         nIndexTarget)
  end subroutine init_coupler

  !=============================================================!

  subroutine check_couple_symm(iComp1,iComp2,NameCaller)
    integer, intent(in)           :: iComp1,iComp2
    character (len=*), intent(in) :: NameCaller

    call check_i_comp(iComp1,NameCaller//' iComp1 ')
    call check_i_comp(iComp2,NameCaller//' iComp2 ')

    if(Couple_CC(iComp1,iComp2)%Dn /= Couple_CC(iComp2,iComp1)%Dn ) &
         call CON_stop(NameCaller//' SWMF_ERROR '//&
         'Couple_CC % Dn is not symmetric for '//&
         NameComp_I(iComp1)//' and '//NameComp_I(iComp2))
    if(Couple_CC(iComp1,iComp2) % Dt /= Couple_CC(iComp2,iComp1) % Dt)then
       call CON_stop(NameCaller//' SWMF_ERROR '//&
            'Couple_CC % Dt is not symmetric for '//&
            NameComp_I(iComp1)//' and '//NameComp_I(iComp2))
    end if

  end subroutine check_couple_symm

  !============================================================================

  subroutine set_router_comm(iComp1,iComp2,iCommRouter,UseMe,iProc)

    ! Return the union communicator iCommRouter for two components indexed by
    ! iComp1 and iComp2
    ! Also return UseMe=.true. if the processor is part of the union
    ! and UseMe=.false. if the processor is not part of the union
    ! If the optional iProc argument is present, it is transformed from
    ! the value valid for iCommWorld to the value valid for iCommRouter.

    character(len=*), parameter :: NameSub=NameMod//'::set_router_comm'
    integer, intent(in) :: iComp1,iComp2
    integer, intent(out):: iCommRouter
    logical, intent(out), optional   :: UseMe
    integer, intent(inout), optional :: iProc
    integer :: iGroup1, iGroup2, iGroupUnion, iProcUnion, iError

    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,': starting for iComp1,iComp2,iProc=',&
         iComp1,iComp2,i_proc()

    iGroup1 = i_group(iComp1)
    iGroup2 = i_group(iComp2)

    !write(*,*)NameSub,': call MPI_group_union for iGroup1,iGroup2 =',&
    !     iGroup1,iGroup2
    call MPI_group_union(iGroup1,iGroup2, iGroupUnion, iError)

    !write(*,*)NameSub,': call MPI_comm_create, iProc=',i_proc()
    call MPI_comm_create(i_comm(),iGroupUnion,iCommRouter,iError)

    call MPI_group_rank(iGroupUnion,iProcUnion,iError)
    !write(*,*)NameSub,': iGroupUnion,iProcUnion=',iGroupUnion,iProcUnion

    if(present(UseMe)) UseMe = iProcUnion /= MPI_UNDEFINED

    if(present(iProc))then
       if(DoTest)write(*,*)NameSub,': bcast iProcUnion from iProc=',iProc
       call MPI_bcast(iProcUnion,1,MPI_INTEGER,iProc,i_comm(),iError)
       iProc=iProcUnion
       if(DoTest)write(*,*)NameSub,': iProc modified to ',iProc
    end if

    call MPI_group_free(iGroupUnion,iError)
    if(DoTest)write(*,*)NameSub,': finished, iProc=',i_proc()

  end subroutine set_router_comm

end module CON_coupler
