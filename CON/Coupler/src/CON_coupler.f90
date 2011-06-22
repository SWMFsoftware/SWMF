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

  SAVE

  !PUBLIC TYPES:
  public :: CoordSystemType ! Extend grid descriptor type with coordinate info
  integer, parameter :: lTypeCoord    = 3
  integer, parameter :: lTypeGeometry = 15
  integer, parameter :: lNameVar      = 200

  type CoordSystemType    
     real, dimension(:),pointer    :: Coord1_I, Coord2_I, Coord3_I
     integer, dimension(3)         :: nCoord_D
     character (len=lTypeCoord)    :: TypeCoord
     character (len=lTypeGeometry) :: TypeGeometry
     integer                       :: nVar
     character (len=lNameVar)      :: NameVar
     real                          :: UnitX
  end type CoordSystemType

  !PUBLIC DATA MEMBERS:
  public :: Grid_C          ! Store the coordinate information for components

  type(CoordSystemType), target :: Grid_C(MaxComp+3)

  public :: Couple_CC           ! Frequency of couplings

  type(FreqType) :: Couple_CC(MaxComp,MaxComp) = &
       FreqType(.false.,-1,-1.0,0,0.0)

  public :: MaxCouple           ! Maximum number of couplings

  integer, parameter :: MaxCouple = 25

  public :: nCouple             ! Actual number of couplings

  integer :: nCouple = MaxCouple

  public :: iCompCoupleOrder_II ! The order of couplings

  ! This is the default order based on the propagation of information
  ! from component to component
  integer :: iCompCoupleOrder_II(2,MaxCouple) = reshape ( (/&
       LC_, SC_, &
       SC_, LC_, &
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
       IM_, GM_, &
       IM_, IE_, &
       PW_, GM_, &
       IE_, GM_, &
       IE_, IM_, &
       IE_, PW_, &
       IE_, PS_, &
       IE_, RB_, &
       IE_, UA_, &
       IH_, OH_, &
       OH_, IH_  &
       /), (/2, MaxCouple/) )

  public :: DoCoupleOnTime_C ! Should do coupling limit the time step

  logical :: DoCoupleOnTime_C(MaxComp) = .true.

  integer, public :: iCompSourceCouple, iCompTargetCouple

  ! no. of variables to be coupled by any pair of source and target components
  integer, public :: nVarCouple, nVarCouple_CC(MaxComp,MaxComp)

  ! no. of state variable groups for which coupling is implemented
  integer,parameter,public     :: nCoupleVarGroup = 7

  ! named indices for variable groups for coupling
  integer,parameter,public :: &
       Bfield_               = 1, &
       AnisoPressure_     = 2, &
       ElectronPressure_  = 3, &
       Wave_              = 4, &
       MultiFluid_        = 5, &
       MultiSpecie_       = 6, &
       Material_          = 7

  logical, public :: &
       DoCoupleVar_V(nCoupleVarGroup) = .false. , &
       DoCoupleVar_VCC(nCoupleVarGroup,MaxComp,MaxComp) = .false.

  ! no. of variable indices known to the coupler
  integer,parameter, public  :: nVarIndexCouple = 12

  ! Fixed indices for mapping actual variable indices
  integer,parameter,public ::     &
       RhoCouple_           = 1,  &
       RhoUxCouple_         = 2,  &
       RhoUzCouple_         = 3,  &
       BxCouple_            = 4,  &
       BzCouple_            = 5,  &
       PCouple_             = 6,  &
       PeCouple_            = 7,  &
       PparCouple_          = 8,  &
       WaveFirstCouple_     = 9,  &
       WaveLastCouple_      = 10, &
       MaterialFirstCouple_ = 11, &
       MaterialLastCouple_  = 12

  ! vector storing the actual values of variable indices inside a
  ! coupled component 
  integer, public  :: &
       iVar_V(nVarIndexCouple) = 0, &
       iVar_VCC(nVarIndexCouple,MaxComp, MaxComp) = 0

  !PUBLIC MEMBER FUNCTIONS:
  public :: set_coord_system    ! Sets coordinate information for a component
  public :: gen_to_stretched    ! Transform generalized coordinates to
  public :: stretched_to_gen    ! stretched ones and vice versa for 
  !                               non-uniform structured grid

  public :: init_coupler        ! Initializes grids and router
  public :: set_router          ! Presets the router indexes
  public :: couple_comp         ! Couple two components via global message pass

  public :: check_couple_symm   ! Check if coupling is symmetric
  public :: set_couple_var_info ! Determine which variables to couple and
  !                               find their indices in each component.

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
  ! 01/30/10 G. Toth - added nVar and NameVar to Grid_C
  ! 04/07/11 R. Oran - added subroutine set_couple_var_info for determining
  !                    which variables should be coupled and finding their
  !                    indices in the source and target components.
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
       nVar,          &! number of variables per cell/node
       NameVar,       &! variable names
       iProc0In,      &
       iCommIn) 

    use ModUtilities, ONLY: lower_case

    character(len=*), parameter :: NameSub='set_coord_system'

    integer, intent(in) :: GridID_
    character(len=*), intent(in)  :: TypeCoord

    ! Optional parameters
    character(len=*), intent(in), optional:: TypeGeometry

    real,             intent(in), optional :: &
         Coord1_I(:), Coord2_I(:), Coord3_I(:)

    real,             intent(in), optional:: UnitX

    integer,          intent(in), optional:: nVar
    character(len=*), intent(in), optional:: NameVar

    integer,          intent(in), optional:: iProc0In, iCommIn

    integer :: iProc0, iComm, iError
    logical :: IsRoot
    type(CoordSystemType), pointer :: ThisGrid
    logical :: DoInit=.true.
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

    ! Broadcast the coordinate type
    if(IsRoot)Thisgrid%TypeCoord=TypeCoord
    call MPI_bcast(ThisGrid%TypeCoord,lTypeCoord,MPI_CHARACTER,&
         iProc0,iComm,iError)

    if(present(TypeGeometry))then
       ! Broadcast the geometry type
       if(IsRoot) Thisgrid%TypeGeometry = TypeGeometry
       call MPI_bcast(ThisGrid%TypeGeometry, lTypeGeometry, MPI_CHARACTER,&
            iProc0, iComm, iError)
    end if

    if(present(nVar))then
       ! Broadcast number of variables per cell
       if(IsRoot) ThisGrid%nVar = nVar
       call MPI_bcast(ThisGrid%nVar, 1, MPI_INTEGER, iProc0, iComm, iError)
    end if

    if(present(NameVar))then
       ! Broadcast list of variable names
       if(IsRoot) then
          ThisGrid%NameVar = NameVar
          call lower_case(ThisGrid%NameVar)
       end if
       call MPI_bcast(ThisGrid%NameVar, lNameVar, MPI_CHARACTER, &
            iProc0, iComm, iError)
    end if

    if(present(UnitX))then
       ! Broadcast the unit of length
       if(IsRoot)Thisgrid%UnitX=UnitX
       call MPI_bcast(ThisGrid%UnitX,1,MPI_REAL,&
            iProc0,iComm,iError)
    end if

    ! Get the size of the coordinate arrays and broadcast it
    if(IsRoot)then
       ThisGrid%nCoord_D = 0
       if(present(Coord1_I)) ThisGrid%nCoord_D(1)=size(Coord1_I)
       if(present(Coord2_I)) ThisGrid%nCoord_D(2)=size(Coord2_I)
       if(present(Coord3_I)) ThisGrid%nCoord_D(3)=size(Coord3_I)
    end if
    call MPI_bcast(ThisGrid%nCoord_D, 3, MPI_INTEGER, iProc0, iComm, iError)

    ! Allocate and broadcast coordinate arrays
    if(ThisGrid%nCoord_D(1)>0)then 
       if(associated(ThisGrid%Coord1_I))deallocate(ThisGrid%Coord1_I)
       allocate(ThisGrid%Coord1_I(ThisGrid%nCoord_D(1)))
       if(IsRoot) ThisGrid%Coord1_I = Coord1_I
       call MPI_bcast(ThisGrid%Coord1_I, ThisGrid%nCoord_D(1), MPI_REAL,&
            iProc0, iComm, iError)
    end if

    if(ThisGrid%nCoord_D(2)>0) then
       if(associated(ThisGrid%Coord2_I))deallocate(ThisGrid%Coord2_I)
       allocate(ThisGrid%Coord2_I(ThisGrid%nCoord_D(2)))
       if(IsRoot) ThisGrid%Coord2_I = Coord2_I
       call MPI_bcast(ThisGrid%Coord2_I, ThisGrid%nCoord_D(2), MPI_REAL, &
            iProc0, iComm, iError)
    end if

    if(ThisGrid%nCoord_D(3)>0)then
       if(associated(ThisGrid%Coord3_I))deallocate(ThisGrid%Coord3_I)
       allocate(ThisGrid%Coord3_I(ThisGrid%nCoord_D(3)))
       if(IsRoot) ThisGrid%Coord3_I = Coord3_I
       call MPI_bcast(ThisGrid%Coord3_I, ThisGrid%nCoord_D(3), MPI_REAL, &
            iProc0, iComm, iError)
    end if

  contains

    subroutine init_coord_system_all
      integer :: iComp
      DoInit=.false.
      do iComp=1, MaxComp+3
         nullify(&
              Grid_C(iComp)%Coord1_I, &
              Grid_C(iComp)%Coord2_I, &
              Grid_C(iComp)%Coord3_I)
         Grid_C(iComp)%nCoord_D    = 0
         Grid_C(iComp)%TypeCoord   = '???'
         Grid_C(iComp)%TypeGeometry = 'cartesian'
         Grid_C(iComp)%nVar         = 0
         Grid_C(iComp)%NameVar      = '???'
         Grid_C(iComp)%UnitX        = 1.0
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
      if (Coord_I(iL) < Coord_I(iU)) then
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
      else
         iU=lbound(Coord_I,1)
         iL=ubound(Coord_I,1)
         if(XyzGen_D(iDim)<=Coord_I(iL))then
            XyzGen_D(iDim)=real(iL) +&
                 OneIfExtrapolate * (XyzGen_D(iDim) - Coord_I(iL))/&
                 (Coord_I(iL) - Coord_I(iL-1))
         elseif(XyzGen_D(iDim)>=Coord_I(iU))then
            XyzGen_D(iDim)=real(iU) +&
                 OneIfExtrapolate * (XyzGen_D(iDim) - Coord_I(iU))/&
                 (Coord_I(iU) - Coord_I(iU+1))
         else
            Number=maxloc(Coord_I,1,MASK=Coord_I<=XyzGen_D(iDim))
            XyzGen_D(iDim)=real(Number)+(XyzGen_D(iDim)-Coord_I(Number))/&
                 (Coord_I(Number)-Coord_I(Number-1))
         end if
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
       IsPeriodic_D, & ! Periodicity for all dimesnsions (op)
       nVar,         & ! Number of variables per grid cell
       NameVar       ) ! Variable names 

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

    integer, intent(in), optional :: nVar
    character(len=lNameVar), intent(in), optional:: NameVar

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
         nVar=nVar,                          &
         NameVar=NameVar,                    &
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
  !========================================================================== 
  subroutine set_couple_var_info(iCompSource, iCompTarget)

    ! Determine coupling flags and variable indices used for accesing
    ! the buffer grid data. 

    ! Handle coupling of any two sets of variables  
    ! with a minimal no. of assumptions about variable indices. 
    ! Allows coupling two components with/without: 
    ! Ppar, Pe, multi/single fluid/specie, waves. 

    ! REVISION HISTORY:    
    ! Feb 2011 - R. Oran - initial version.

    use ModProcessVarName,  ONLY: process_var_name, nVarMax
    use ModUtilities,       ONLY: split_string, join_string

    integer,intent(in)             :: iCompSource, iCompTarget
    character(len=1500)            :: NameVarSource, NameVarTarget
    character(len=15),allocatable  :: NameVarSource_V(:), NameVarTarget_V(:)
    character(len=15)              :: NameList_V(nVarMax)
    logical,allocatable            :: IsFoundVarSource_V(:)
    integer      :: nVarSource, nVarTarget, iVarSource, iVarTarget
    integer      :: nDensitySource, nDensityTarget, nDensityCouple
    integer      :: nSpeedSource, nSpeedTarget, nSpeedCouple
    integer      :: nPSource, nPTarget, nPCouple
    integer      :: nPparSource, nPparTarget, nPparCouple
    integer      :: nWaveSource, nWaveTarget, nWaveCouple 
    integer      :: nMaterialSource, nMaterialTarget, nMaterialCouple      
    logical      :: DoTest, DoTestMe

    character(len=*), parameter    :: NameSub = NameMod//'::set_couple_var_info'
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if (DoTest .and. i_proc() == 0 ) write(*,*) NameSub, ' started.'

    iCompSourceCouple = iCompSource
    iCompTargetCouple = iCompTarget

    ! The following coupling flags are set to true if both source and target
    ! have the relevant state variables, see below. 
    ! NOTE: If both components have multiple densities, all species/fluids state
    ! variables should be coupled, hence a further check is made to ensure that
    ! the fluids/species in both components are identical.        

    DoCoupleVar_V = .false.
    ! The elements of DoCoupleVar_V will be set to true if:
    ! Bfield_               :  Both have a magnetic field.
    ! AnisoPressure_     :  Both use anisotropic pressure.
    ! ElectronPressure_  :  Both use electron pressure.
    ! Wave_              :  Both use the same # of waves.
    ! MultiFluid_        :  Both use the same neutral fluids
    ! MultiSpecie_       :  Both use the same species     
    ! Material_          :  Both use the same # of materials

    nDensitySource  = 0 ; nDensityTarget  = 0 ; nDensityCouple  = 0
    nSpeedSource    = 0 ; nSpeedTarget    = 0 ; nSpeedCouple    = 0
    nPSource        = 0 ; nPTarget        = 0 ; nPCouple        = 0
    nPparSource     = 0 ; nPpartarget     = 0 ; nPparCouple     = 0
    nWaveSource     = 0 ; nWaveTarget     = 0 ; nWaveCouple     = 0
    nMaterialSource = 0 ; nMaterialTarget = 0 ; nMaterialCouple = 0

    !\        
    ! process variable names
    !/                                                  
    ! Here the external subroutine process_var_names is called, which will cast  
    ! the variable names in the source and target components into a            
    ! standardized form.                                         
    ! Variable names are defined in equation modules  and often different names
    ! are given to the same physical quantity. This is especially relevant to
    ! multi fluid/specie models.                                              
    ! The subroutine also returns the number of distinct densities and speeds,
    ! pressure (total and parallel), number of waves and material.
    ! See description inside the subroutine for more details.            

    NameVarSource = ' '//Grid_C(iCompSource)%NameVar
    NameVarTarget = ' '//Grid_C(iCompTarget)%NameVar

    ! Separate NameVar(Source/Target) into a string arrays of variable names. 
    call split_string(NameVarSource, nVarMax, NameList_V, nVarSource)
    allocate(NameVarSource_V(nVarSource))
    NameVarSource_V(1:nVarSource) = NameList_V(1:nVarSource)
    NameList_V = ''

    call split_string(NameVarTarget, nVarMax, NameList_V, nVarTarget)
    allocate(NameVarTarget_V(nVarTarget))
    NameVarTarget_V(1:nVarTarget) = NameList_V(1:nVarTarget)

    ! process source variable names                
    if (DoTest) then
       write(*,*) NameComp_I(iCompSource), ':  variable names for I/O:'
       write(*,*) trim(NameVarSource)
    end if
    call process_var_name(nVarSource, NameVarSource_V, &
         nDensitySource, nSpeedSource, nPSource, nPparSource, &
         nWaveSource, nMaterialSource)
    call join_string(nVarSource, NameVarSource_V,NameVarSource)
    if(DoTest) then
       write(*,*) NameComp_I(iCompSource), ':  variable names used by coupler:'
       write(*,*) trim(NameVarSource)
       write(*,*) ' '
    end if

    ! process target variable names       
    if(DoTest) then
       write(*,*) NameComp_I(iCompTarget), ': variable names for I/O:'
       write(*,*) trim(NameVarTarget)
    end if
    call process_var_name(nVarTarget, NameVarTarget_V, &
         nDensityTarget, nSpeedTarget, nPTarget, nPparTarget, &
         nWaveTarget, nMaterialTarget)
    call join_string(nVarTarget, NameVarTarget_V,NameVarTarget)
    if (DoTest) then
       write(*,*) NameComp_I(iCompTarget), ': variabale names used by coupler:'
       write(*,*) trim(NameVarTarget)
       write(*,*) ' '
    end if

    !\                                                          
    ! Get info required for coupling, depending on which variables are present 
    !/                                                 
    ! Logical array to report whether a variable name in the source is treated
    ! by any of the cases below and thus it can be handeled by the coupler.
    ! For brevity, all elements are initialized to .true.
    ! An element will be set to .false. if the variable name is not found.
    allocate(IsFoundVarSource_V(nVarSource))
    IsFoundVarSource_V = .true.

    do iVarSource = 1, nVarSource

       select case(NameVarSource_V(iVarSource))

       case('Rho', 'P', 'Ew','Eint', 'Hyp', 'My', 'Mz', 'By', 'Bz')
          ! Do nothing.
          ! Rho, P assumed to always be present
          ! ew, EInt, hyp : internal variables, not to be coupled.      
          ! By, Bz, My, Mz already covered by other cases.       

       case('Mx')
          ! Check that My and Mz immediately follow Mx 
          if ( NameVarSource_V(iVarSource + 1) /= 'My' .and. &
               NameVarSource_V(iVarSource + 2) /= 'Mz') then
             write(*,*) 'Error in ModEquation for ', NameComp_I(iCompSource)
             write(*,*) 'X, Y and Z Momenta must have consecutive indices!'
             call CON_stop(NameSub)
          end if

       case('Bx')
          if(index(NameVarTarget,' Bx ') > 0) then
             DoCoupleVar_V(Bfield_) = .true.
             ! Check that By and Bz immediately follow Bx 
             if ( NameVarSource_V(iVarSource + 1) /= 'By' .and. &
                  NameVarSource_V(iVarSource + 2) /= 'Bz') then
                write(*,*) 'Error in ModEquation for ', NameComp_I(iCompSource)
                write(*,*) 'Bx, By and Bz must have consecutive indices!'
                call CON_stop('ERROR in'//NameSub)
             end if
          end if

       case('Pe')
          DoCoupleVar_V(ElectronPressure_) = index(NameVarTarget, ' Pe ') > 0

       case('Ppar')
          DoCoupleVar_V(AnisoPressure_) = index(NameVarTarget,' Ppar ') > 0

       case('i01')
          ! Enumerated names for waves                                     
          ! Coupling components with and without waves is allowed.             
          ! If waves are present in both source and target, they should have 
          ! the same number. Otherwise a CON_stop is issued.                

          if(nWaveSource == nWaveTarget)then
             nWaveCouple  = nWaveSource
             DoCoupleVar_V(Wave_) = .true.
          else if (nWaveTarget == 0) then
             if(i_proc()==0) &
                  write(*,*) 'Coupling components with and without waves!!!'
          else
             write(*,*) 'SWMF error found by ',NameSub
             write(*,*) 'Cannot couple components with different nWave>0!'
             call CON_stop(NameSub//': change nWave (use Config.pl).')
          end if

       case('m01')
          ! Verify that target uses the same # of materials.
          ! Otherwise stop with error.                                       
          if(nMaterialSource == nMaterialTarget)then
             DoCoupleVar_V(Material_) = .true.
             nMaterialCouple = nMaterialSource
          else
             write(*,*) 'SWMF error found by ',NameSub
             write(*,*) 'Cannot couple components with different nMaterial!'
             call CON_stop(NameSub//': change nMaterial (use Config.pl).')
          end if

       case default
          ! The only variable names that are left are associated with either:
          !   - a neutral/ionized fluid other than the main (M)HD component.
          !   - additional waves.                 
          !   - additional materials.         

          if(nDensitySource == 1 .and. &
               nWaveSource < 1 .and. nMaterialSource < 1) then
             ! Source variable name is unknown, stop.                        
             write(*,*) 'SWMF error found in ', NameSub
             write(*,*) 'Coupling of variable ', NameVarSource_V(iVarSource)
             write(*,*) ' used by ', NameComp_I(iCompSource), &
                  ' component is undefined!'
             call CON_stop(NameSub//' check variable names!')
          end if

          ! Report status of this variable as not found.   
          ! This will be corrected for waves and materials after looping over
          ! all names is complete.                            
          IsFoundVarSource_V(iVarSource) = .false.
       end select
    end do

    ! Correct "found" status for waves and materials  
    if (nWaveSource >= 1) &
         IsFoundVarSource_V(WaveFirstCouple_:WaveLastCouple_) = .true.
    if (nMaterialSource >= 1) &
         IsFoundVarSource_V(MaterialFirstCouple_:MaterialLastCouple_) = .true.

    !\                                                    
    ! Check multi fluid/specie variables.
    !/
    if ( nDensitySource == nDensityTarget .and. &
         nSpeedSource == nSpeedTarget .and. nDensitySource >1) then

       ! Verify that the same specie/fluid variables are present in both
       ! source and target.                               
       SOURCELOOP: do iVarSource = 1,nVarSource

          ! Only check varaiables that were not yet accounted for. 
          if(IsFoundVarSource_V(iVarSource)) CYCLE

          ! Look up source variable name in the target
          do iVarTarget = 1, nVarTarget
             if (NameVarSource_V(iVarSource) == NameVarTarget_V(iVarTarget)) then
                IsFoundVarSource_V(iVarSource) = .true.
                CYCLE SOURCELOOP
             end if
          end do

          if(.not. IsFoundVarSource_V(iVarSource)) then
             ! At least one fluids/specie name does not match 
             write(*,*) 'SWMF error found in ', NameSub
             write(*,*) 'Coupled components with unmatching densities/speeds:'
             write(*,*) 'Component ', NameComp_I(iCompSource),' uses '// &
                  NameVarSource_V(iVarSource)
             write(*,*) 'Component ', NameComp_I(iCompTarget), 'does not.'
             call CON_stop(NameSub//': check ModEquation and recompile!')
          end if
       end do SOURCELOOP

       if(nSpeedSource > 1) then
          DoCoupleVar_V(MultiFluid_) = .true.
       else
          DoCoupleVar_V(MultiSpecie_) = .true.
       end if
    end if
    if(nDensitySource /= nDensityTarget .or. &
         nSpeedSource /= nSpeedTarget) then
       if( (nDensitySource == 1 .and. nSpeedSource == 1) .or. &
            (nDensityTarget == 1 .and. nSpeedTarget == 1) ) then
          ! Coupling multi<-> single fluid is allowed.        
          ! In this case the additional fluid variables are not coupled, and
          ! their boundary conditions should be implemented separately. 
          DoCoupleVar_V(MultiFluid_) = .false.
          DoCoupleVar_V(MultiSpecie_)   = .false.
       else
          ! Both components have differing # of multiple densities or speeds 
          write(*,*) 'SWMF error found in ', NameSub
          write(*,*) 'Coupled SWMF components use different # of fluids/species!'
          call CON_stop(NameSub//': check ModEquation and recompile!')
       end if
    end if

    ! calculate  indices and # of variables transfered to buffer grid             
    nDensityCouple = min(nDensitySource, nDensityTarget)
    nSpeedCouple   = min(nSpeedSource,   nSpeedTarget)
    nPCouple       = min(nPSource,       nPTarget)
    nPparCouple    = min(nPparSource,    nPparTarget)

    iVar_V(RhoCouple_) = 1 
    iVar_V(RhoUxCouple_) = 2
    iVar_V(RhoUzCouple_) = 4
    nVarCouple = 4

    iVar_V(PCouple_) = nVarCouple + 1
    nVarCouple = nVarCouple + 1

    if (DoCoupleVar_V(Bfield_)) then
       iVar_V(BxCouple_) = nVarCouple + 1
       iVar_V(BzCouple_) = nVarCouple + 3
       nVarCouple = nVarCouple + 3
    end if

    if (DoCoupleVar_V(AnisoPressure_)) then      
       nVarCouple = nVarCouple + 1
       iVar_V(PparCouple_) = nVarCouple
    end if

    if (DoCoupleVar_V(ElectronPressure_)) then
       nVarCouple = nVarCouple + 1
       iVar_V(PeCouple_) = nVarCouple
    end if

    if (DoCoupleVar_V(Wave_)) then
       iVar_V(WaveFirstCouple_) = nVarCouple + 1
       iVar_V(WaveLastCouple_)  = iVar_V(WaveFirstCouple_) + nWaveSource
       nVarCouple = iVar_V(WaveLastCouple_)
    end if

    if (DoCoupleVar_V( Material_)) then
       iVar_V(MaterialFirstCouple_) = nVarCouple + 1
       iVar_V(MaterialLastCouple_)  = &
            iVar_V(MaterialFirstCouple_) + nMaterialSource
       nVarCouple = iVar_V(MaterialLastCouple_)
    end if

    if (nVarCouple >  nVarSource) then
       write(*,*) 'SWMF Error: # of coupled variables exceeds nVarSource'
       call CON_stop(NameSub//' error in calculating nVarCouple')
    end if

    ! Store coupling info to avoid recalculation at next coupling time
    DoCoupleVar_VCC(:,iCompSource, iCompTarget) = DoCoupleVar_V
    DoCoupleVar_VCC(:,iCompTarget, iCompSource) = DoCoupleVar_V

    iVar_VCC(:,iCompSource, iCompTarget) = iVar_V
    iVar_VCC(:,iCompTarget, iCompSource) = iVar_V

    nVarCouple_CC(iCompSource,iCompTarget) = nVarCouple
    nVarCouple_CC(iCompTarget,iCompSource) = nVarCouple


    if( i_proc() ==0) then !<<<<<<<< change to if(DoTest) ??                  
       write(*,*) '---------------------------------------------'
       write(*,*) NameSub,':'
       write(*,*) 'Source: ', NameComp_I(iCompSource)
       write(*,*) 'Target: ', NameComp_I(iCompTarget)
       write(*,*) 'nDensity(Source/Target/Couple):'
       write(*,*) nDensitySource, nDensityTarget,nDensityCouple
       write(*,*) 'nSpeed(Source/Target/Couple):'
       write(*,*) nSpeedSource, nSpeedTarget,nSpeedCouple
       write(*,*) 'nP(Source/Target/Couple):'
       write(*,*) nPSource, nPTarget,nPCouple
       write(*,*) 'nPpar(Source/Target/Couple):'
       write(*,*) nPparSource, nPparTarget,nPparCouple
       write(*,*) '---------------------------------------------'
       write(*,*) NameSub,' for:'
       write(*,*) 'Source: ',         NameComp_I(iCompSource)
       write(*,*) 'Target: ',         NameComp_I(iCompTarget)
       write(*,*) 'Coupling flags:'
       write(*,*) 'Magnetic field: ', DoCoupleVar_V(Bfield_)
       write(*,*) 'Pe: ',             DoCoupleVar_V(ElectronPressure_)
       write(*,*) 'Ppar: ',           DoCoupleVar_V(AnisoPressure_)
       write(*,*) 'Waves: ',          DoCoupleVar_V(Wave_)
       write(*,*) 'Neutrals: ',       DoCoupleVar_V(MultiFluid_)
       write(*,*) 'Ions: ',           DoCoupleVar_V(MultiSpecie_)
       write(*,*) '---------------------------------------------'
       write(*,*) 'Coupled variable indices in source component:'
       write(*,*) 'Rho: ',      iVar_V(RhoCouple_)
       write(*,*) 'RhoUx: ',    iVar_V(RhoUxCouple_)
       write(*,*) 'RhoUz: ',    iVar_V(RhoUzCouple_)
       write(*,*) 'Bx: ',       iVar_V(BxCouple_)
       write(*,*) 'Bz: ',       iVar_V(BzCouple_)
       write(*,*) 'P: ',        iVar_V(PCouple_)
       write(*,*) 'Pe: ',       iVar_V(PeCouple_)
       write(*,*) 'Ppar: ',     iVar_V(PparCouple_)
       write(*,*) 'WaveFirst: ',iVar_V(WaveFirstCouple_)
       write(*,*) 'WaveLast: ', iVar_V(WaveLastCouple_)
       write(*,*) 'nVarCouple: ',nVarCouple
       write(*,*) '---------------------------------------------'

    end if

    deallocate(NameVarSource_V, NameVarTarget_V, IsFoundVarSource_V)

  end subroutine set_couple_var_info
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
