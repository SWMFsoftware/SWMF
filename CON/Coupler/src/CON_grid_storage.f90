!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module CON_grid_storage

  ! revision history:
  ! 09SEP03              I.Sokolov, igorsok@umich.edu - initial prototype/code
  ! 12SEP03              version for any operating system
  ! 16JAN05              G.Toth removed the obsolete GmIe_grid

  use CON_comp_param, ONLY: EE_, GM_, IE_, IH_, IM_, OH_, PC_, PS_, &
       PT_, PW_, RB_, SC_, SP_, UA_, CZ_
  use CON_world
  use CON_domain_decomposition, ONLY: DomainType, DomainPointerType, &
       init_decomposition_dd, get_root_decomposition_dd,  &
       bcast_decomposition_dd, synchronize_refinement_dd, &
       associate_dd_pointer_dd
  use ModUtilities, ONLY: CON_stop

  implicit none

  integer, public, parameter:: MaxGrid = MaxComp

  type(DomainPointerType), private :: Domain_I(MaxGrid)
  logical   :: DoneDomainInit_C(MaxGrid)= .false.
  ! Introduced to bypass the HALEM compiler restrictions
  integer :: nDim_C(MaxGrid) = 0
  interface init_decomposition
     module procedure init_decomposition_id
     module procedure init_decomposition_dd
  end interface

  interface get_root_decomposition
     module procedure get_root_decomposition_id
     module procedure get_root_decomposition_dd
  end interface

  interface bcast_decomposition
     module procedure bcast_decomposition_id
     module procedure bcast_decomposition_dd
  end interface

  interface synchronize_refinement
     module procedure synchronize_refinement_id
     module procedure synchronize_refinement_dd
  end interface

  interface associate_dd_pointer
     module procedure associate_dd_pointer_dd
     module procedure associate_dd_pointer_id
  end interface associate_dd_pointer

  type(DomainType), save, target :: &
       EeGrid, GmGrid, IeGrid, IhGrid, ImGrid, OhGrid, PcGrid, PsGrid, &
       PtGrid, PwGrid, RbGrid, ScGrid, SpGrid, UaGrid, CzGrid

contains
  !============================================================================
  subroutine init_grid_storage(Domain_I, GridID_)

    integer,                 intent(in )   :: GridID_
    type(DomainPointerType), intent(inout) :: Domain_I(MaxGrid)
    ! information for the global grids is stored at all PEs so it is
    ! important to reduce the memory requirements. This short procedure
    ! describes how the memory is allocated for the domain decomposition
    ! structure. This solution satisfies the most picky SGI compiler,
    ! but requires to add manually an identifier for the domain
    ! decomposition while adding a new component to the framework.

    character(len=*), parameter:: NameSub = 'init_grid_storage'
    !--------------------------------------------------------------------------
    select case(GridID_)
    case(EE_)
       Domain_I(GridID_)%Ptr => EeGrid
    case(GM_)
       Domain_I(GridID_)%Ptr => GmGrid
    case(IE_)
       Domain_I(GridID_)%Ptr => IeGrid
    case(IM_)
       Domain_I(GridID_)%Ptr => ImGrid
    case(IH_)
       Domain_I(GridID_)%Ptr => IhGrid
    case(OH_)
       Domain_I(GridID_)%Ptr => OhGrid
    case(PC_)
       Domain_I(GridID_)%Ptr => PcGrid
    case(PS_)
       Domain_I(GridID_)%Ptr => PsGrid
    case(PT_)
       Domain_I(GridID_)%Ptr => PtGrid
    case(PW_)
       Domain_I(GridID_)%Ptr => PwGrid
    case(RB_)
       Domain_I(GridID_)%Ptr => RbGrid
    case(SP_)
       Domain_I(GridID_)%Ptr => SpGrid
    case(SC_)
       Domain_I(GridID_)%Ptr => ScGrid
    case(UA_)
       Domain_I(GridID_)%Ptr => UaGrid
    case(CZ_)
       Domain_I(GridID_)%Ptr => CzGrid
    case default
       write(*,*)'ERROR in ModInitGridStorage: GridID = ',GridID_
       call CON_stop(NameSub//': not implemented grid ID')
    end select

  end subroutine init_grid_storage
  !============================================================================
  logical function done_dd_init(GridID_)
    integer,intent(in)::GridID_
    ! Returns .true. if the domain decomposition is initialaized
    !--------------------------------------------------------------------------
    done_dd_init=.false.
    if(GridID_>0 .and. GridID_<=MaxGrid)&
         done_dd_init = DoneDomainInit_C(GridID_)
  end function done_dd_init
  !============================================================================
  subroutine init_decomposition_id(&
       GridID_, CompID_, nDim, IsTreeDD)

    ! Initialization for the decomposition.
    ! Note that if the decomposition is a tree decomposition and
    ! another order of children is accepted in the component,
    ! different from that assumed to be standard one
    ! the array iShift_DI should be reassigned in set_root.

    integer,           intent(in) :: GridID_
    integer,           intent(in) :: CompID_,nDim
    logical, optional, intent(in) :: IsTreeDD
    !--------------------------------------------------------------------------
    if(GridID_ > MaxGrid .or. GridID_ <=0 )call CON_stop(&
         'Prrohibited value for GridID_')
    if(done_dd_init(GridID_)) RETURN
    call init_grid_storage(Domain_I, GridID_)
    call init_decomposition_dd(Domain_I(GridID_)%Ptr, CompID_, nDim, IsTreeDD)
    DoneDomainInit_C(GridID_) = .true.
    nDim_C(GridID_) = nDim

  end subroutine init_decomposition_id
  !============================================================================
  subroutine get_root_decomposition_id(&
       GridID_,      &! ID for Decomposition to be constructed
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,   &! As in DomainType
       CoordMax_D,   &! As in DomainType
       nCell_D,      &! As in DomainType
       PE_I,         &! PE layout
       iBlock_I,     &! Local Block Number layout
       IsPeriodic_D, &! As in DomainType
       iShift_DI,    &! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)

    ! To get a decomposition domain, even the tree one, the root
    ! decomposition should be first constructed. PE here are the ranks in the
    ! LOCAL communicator for the component

    integer, intent(in) :: GridID_
    integer, intent(in) :: iRootMapDim_D(nDim_C(GridID_))
    real,    intent(in) :: CoordMin_D(nDim_C(GridID_))
    real,    intent(in) :: CoordMax_D(nDim_C(GridID_))
    integer, intent(in) :: nCell_D(nDim_C(GridID_))
    integer, optional, intent(in) :: PE_I(:)
    integer, optional, intent(in) :: iBlock_I(:)
    logical, optional, intent(in) :: IsPeriodic_D(nDim_C(GridID_))
    integer, optional, intent(in) :: iShift_DI(:,:)
    logical, optional, intent(in) :: DoGlueMargins
    integer, optional, intent(in) :: iDirMinusGlue
    integer, optional, intent(in) :: iDirPlusGlue
    integer, optional, intent(in) :: iDirCycle
    !--------------------------------------------------------------------------
    call get_root_decomposition_dd(Domain_I(GridID_)%Ptr, &
       iRootMapDim_D,&! As in DomainType
       CoordMin_D,   &! As in DomainType
       CoordMax_D,   &! As in DomainType
       nCell_D,      &! As in DomainType
       PE_I,         &! PE layout
       iBlock_I,     &! Local Block Number layout
       IsPeriodic_D, &! As in DomainType
       iShift_DI,    &! As in DomainType
       DoGlueMargins,&! As in DomainType
       iDirMinusGlue,&! As in DomainType
       iDirPlusGlue,& ! As in DomainType
       iDirCycle)     ! As in DomainType)

  end subroutine get_root_decomposition_id
  !============================================================================
  subroutine bcast_decomposition_id(GridID_)

    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    call bcast_decomposition_dd(Domain_I(GridID_)%Ptr)

  end subroutine bcast_decomposition_id
  !============================================================================
  subroutine synchronize_refinement_id(&
       GridID_, LocalDomain, iProcUnion, iCommUnion)

    ! If any of the optional parameters is not present, the global
    ! decomposition at all the PEs of the global communicator is
    ! synchronized with the local one at root processor of the
    ! component. Otherwise the global
    ! decomposition at all the PEs of the communicator iProcUnion is
    ! synchronized with the local one at the PE having the rank
    ! iProcUnion in the communicator iCommUnion.
    ! Recalculate local PE ranks (of the local grid) to their values
    ! in the global communicator (i_comm()).

    integer,           intent(in) :: GridID_      ! Global grid ID
    type(DomainType),  intent(in) :: LocalDomain  ! Local Grid with which the
    integer, optional, intent(in) :: iProcUnion,iCommUnion
    !--------------------------------------------------------------------------
    call synchronize_refinement_dd(&
            Domain_I(GridID_)%Ptr, LocalDomain, iProcUnion, iCommUnion)

  end subroutine synchronize_refinement_id
  !============================================================================
  integer function i_realization(GridID_)

    integer, intent(in) :: GridID_
    !--------------------------------------------------------------------------
    i_realization = Domain_I(GridID_)%Ptr%iRealization

  end function i_realization
  !============================================================================
  subroutine associate_dd_pointer_id(GridID_,DomainPointer)

    integer,                 intent(in) :: GridID_
    type(DomainPointerType), intent(out):: DomainPointer
    !--------------------------------------------------------------------------
    nullify(DomainPointer%Ptr)
    DomainPointer%Ptr => Domain_I(GridID_)%Ptr

  end subroutine associate_dd_pointer_id
  !============================================================================
  integer function ndim_id(GridID_)

    integer, intent(in) :: GridID_
    !--------------------------------------------------------------------------
    ndim_id = Domain_I(GridID_)%Ptr%nDim

  end function ndim_id
  !============================================================================
  function ncell_id(GridID_)

    integer,intent(in) :: GridID_
    integer :: ncell_id(nDim_C(GridID_))
    !--------------------------------------------------------------------------
    ncell_id = Domain_I(GridID_)%Ptr%nCell_D

  end function ncell_id
  !============================================================================
  integer function compid(GridID_)

    integer,intent(in)::GridID_
    !--------------------------------------------------------------------------
    compid = Domain_I(GridID_)%Ptr%CompID_

  end function compid
  !============================================================================
end module CON_grid_storage
!==============================================================================

