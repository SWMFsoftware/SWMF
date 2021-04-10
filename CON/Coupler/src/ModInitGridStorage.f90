!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModInitGridStorage
  use CON_world, ONLY: MaxComp
  use CON_comp_param, ONLY: EE_, GM_, IE_, IH_, IM_, OH_, PC_, PS_, &
       PT_, PW_, RB_, SC_, SP_, UA_, CZ_
  use CON_domain_decomposition, ONLY: DomainPointerType, DomainType
  implicit none
  PRIVATE ! Except
  integer, public, parameter:: MaxGrid = MaxComp
  type(DomainType), save, target :: &
       EeGrid, GmGrid, IeGrid, IhGrid, ImGrid, OhGrid, PcGrid, PsGrid, &
       PtGrid, PwGrid, RbGrid, ScGrid, SpGrid, UaGrid, CzGrid
  ! Public member:
  ! Sets pointer for the component on target DD
  public :: init_grid_storage
contains
  !============================================================================
  ! revision history:
  ! 09SEP03              I.Sokolov<igorsok@umich.edu - initial prototype/code
  ! 12SEP03              version for any operating system
  ! 16JAN05              G.Toth removed the obsolete GmIe_grid
  subroutine init_grid_storage(Domain_I, GridID_)
    integer,                 intent(in )   :: GridID_
    type(DomainPointerType), intent(inout) :: Domain_I(MaxGrid)
    ! information for the global grids is stored at all PEs so it is
    ! important to reduce the memory requirements. This short procedure
    ! describes how the memory is allocated for the domain decomposition
    ! structure. This solution satisfies the most picky SGI compiler,
    ! but requires to add manually an identifier for the domain
    ! decomposition while adding a new component to the framework.
    !--------------------------------------------------------------------------
    select case(GridID_)
    case(EE_)
       Domain_I(GridID_)%Ptr=>EeGrid
    case(GM_)
       Domain_I(GridID_)%Ptr=>GmGrid
    case(IE_)
       Domain_I(GridID_)%Ptr=>IeGrid
    case(IM_)
       Domain_I(GridID_)%Ptr=>ImGrid
    case(IH_)
       Domain_I(GridID_)%Ptr=>IhGrid
    case(OH_)
       Domain_I(GridID_)%Ptr=>OhGrid
    case(PC_)
       Domain_I(GridID_)%Ptr=>PcGrid
    case(PS_)
       Domain_I(GridID_)%Ptr=>PsGrid
    case(PT_)
       Domain_I(GridID_)%Ptr=>PtGrid
    case(PW_)
       Domain_I(GridID_)%Ptr=>PwGrid
    case(RB_)
       Domain_I(GridID_)%Ptr=>RbGrid
    case(SP_)
       Domain_I(GridID_)%Ptr=>SpGrid
    case(SC_)
       Domain_I(GridID_)%Ptr=>ScGrid
    case(UA_)
       Domain_I(GridID_)%Ptr=>UaGrid
    case(CZ_)
       Domain_I(GridID_)%Ptr=>CzGrid
    case default
       write(*,*)'ERROR in ModInitGridStorage: GridID = ',GridID_
       call CON_stop('ERRORin ModInitGridStorage: not implemented grid ID')
    end select
  end subroutine init_grid_storage
  !============================================================================
end module ModInitGridStorage
!==============================================================================
