!BOP
!MODULE: ModInitGridStorage - set the number of grids and optimize the memory to store them
!INTERFACE:
module ModInitGridStorage
  !USES:
  use CON_world,ONLY:MaxComp
  use CON_comp_param
  use CON_domain_decomposition,ONLY:DDPointerType,DomainDecompositionType
  implicit none
  !EOP
  integer,parameter::MaxGrid=MaxComp+3
  integer,parameter::GmIeGrid_=MaxComp+1
  type(DomainDecompositionType),private,save,target::&
       GmIeGrid,GmGrid,IhGrid,UaGrid,IeGrid,ImGrid,RbGrid,SpGrid,ScGrid
contains
  !BOP
  !REVISION HISTORY:
  !09SEP03              I.Sokolov<igorsok@umich.edu - initial prototype/code
  !12SEP03              version for any operational system:  code/prolog
  !BOP
  !IROUTINE: init_grid_storage - initialize a storage for data relating to the component domain decomposition
  !INTERFACE:
  subroutine init_grid_storage(DD_I,GridID_)
    !USES:
    !INPUT ARGUMENTS:
    integer,intent(in)::GridID_
    type(DDPointerType),dimension(MaxGrid),intent(inout)::DD_I
    !DESCRIPTION: 
    !             information for the global grids is stored at each of the PEs so it is
    !             important to reduce the memory requirements. This short procedure describes
    !             how the memory is allocated for the domain decomposition structure. 
    !             This solution satisfies the most picky SCI compiler, but requires to add manually
    !             an identifier for the domain decomposition while addung a new component to the framework. 
    !EOP
    select case(GridID_)
    case(GmIeGrid_)
       DD_I(GridID_)%Ptr=>GmIeGrid
    case(GM_)
       DD_I(GridID_)%Ptr=>GmGrid
    case(IE_)
       DD_I(GridID_)%Ptr=>IeGrid
    case(IM_)
       DD_I(GridID_)%Ptr=>ImGrid
    case(RB_)
       DD_I(GridID_)%Ptr=>RbGrid
    case(UA_)
       DD_I(GridID_)%Ptr=>UaGrid
    case(IH_)
       DD_I(GridID_)%Ptr=>IhGrid
    case(SP_)
       DD_I(GridID_)%Ptr=>SpGrid
    case(SC_)
       DD_I(GridID_)%Ptr=>ScGrid
    case default
       call CON_stop('Not implemented grid #',GridID_)
    end select
  end subroutine init_grid_storage
end module ModInitGridStorage
