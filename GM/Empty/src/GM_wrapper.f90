!^CFG COPYRIGHT UM
! Wrapper for Global Magnetosphere (GM) component
!==========================================================================
subroutine GM_set_param(CompInfo, TypeAction)

  use CON_comp_info

  implicit none

  character (len=*), parameter :: NameSub='GM_set_param'

  ! Arguments
  type(CompInfoType), intent(inout) :: CompInfo   ! Information for this comp.
  character (len=*), intent(in)     :: TypeAction ! What to do
  !-------------------------------------------------------------------------
  select case(TypeAction)
  case('VERSION')
     call put(CompInfo,&
          Use        =.false., &
          NameVersion='Empty', &
          Version    =0.0)

  case default
     call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
  end select

end subroutine GM_set_param

!==============================================================================

subroutine GM_init_session(iSession, TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  integer,  intent(in) :: iSession         ! session number (starting from 1)
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_init_session'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_init_session

!==============================================================================

subroutine GM_finalize(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_finalize'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_finalize

!==============================================================================

subroutine GM_save_restart(TimeSimulation)

  implicit none

  !INPUT PARAMETERS:
  real,     intent(in) :: TimeSimulation   ! seconds from start time

  character(len=*), parameter :: NameSub='GM_save_restart'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_save_restart

!==============================================================================

subroutine GM_run(TimeSimulation,TimeSimulationLimit)

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  real, intent(inout) :: TimeSimulation   ! current time of component

  !INPUT ARGUMENTS:
  real, intent(in) :: TimeSimulationLimit ! simulation time not to be exceeded

  character(len=*), parameter :: NameSub='GM_run'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_run

!==============================================================================

subroutine GM_get_for_ie_swmf(nPartial,iGet,Get,W,State_V,nVar)

  character(len=*), parameter :: NameSub='GM_get_for_ie_swmf'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_get_for_ie_swmf

!==============================================================================
subroutine GMIE_set_grid
  implicit none
  call CON_stop('GMIE_set_grid: GM_ERROR: empty version cannot be used!')
end subroutine GMIE_set_grid
!==============================================================================
subroutine GM_get_mapping_param_for_ie(& 
     rCurrentGm,rIonosphere)    !srcGM/map_ionosphere_bc.f90
  implicit none
  real,intent(out)::rCurrentGm,rIonosphere
  call CON_stop(&
       'GM_get_mapping_param_for_ie: GM_ERROR: empty version cannot be used!')
end subroutine GM_get_mapping_param_for_ie
!==============================================================================
subroutine GM_get_for_im(Buffer_IIV,iSize,jSize,nVar,NameVar)
  implicit none

  integer, intent(in) :: iSize,jSize,nVar
  real, intent(out), dimension(iSize,jSize,nVar) :: Buffer_IIV
  character (len=*), intent(in) :: NameVar

  character (len=*), parameter :: NameSub='GM_get_for_im'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
end subroutine GM_get_for_im
!==============================================================================
subroutine GM_get_for_rb(Buffer_IIV,iSize,jSize,nVar,NameVar)
  implicit none

  integer, intent(in) :: iSize,jSize,nVar
  real, intent(out), dimension(iSize,jSize,nVar) :: Buffer_IIV
  character (len=*), intent(in) :: NameVar

  character (len=*), parameter :: NameSub='GM_get_for_rb'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
end subroutine GM_get_for_rb
!==============================================================================
subroutine GM_get_for_ie(Buffer_IV,nPoint,nVar,NameVar)
  implicit none

  integer, intent(in) :: nPoint,nVar
  real, intent(out), dimension(nPoint,nVar) :: Buffer_IV
  character (len=*), intent(in) :: NameVar

  character (len=*), parameter :: NameSub='GM_get_for_ie'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
end subroutine GM_get_for_ie
!==============================================================================
subroutine GM_put_from_im(Buffer_II,iSizeIn,jSizeIn,NameVar)
  implicit none

  integer, intent(in) :: iSizeIn,jSizeIn
  real, intent(in) :: Buffer_II(iSizeIn,jSizeIn)
  character(len=*), intent(in) :: NameVar

  character(len=*), parameter :: NameSub='GM_put_from_im'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')
end subroutine GM_put_from_im
!==============================================================================

!==========================================================================
subroutine GM_put_from_ie(Buffer_II,iSize,jSize,NameVar)

  implicit none

  character(len=*), parameter :: NameSub='GM_put_from_ie'

  integer, intent(in) :: iSize,jSize
  real, intent(in) :: Buffer_II(iSize,jSize)
  character(len=*), intent(in) :: NameVar

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_put_from_ie
!==============================================================================

!==============================================================================

subroutine GM_calc_fac(IsNewFacPoint, nPoint_I, ColatLim_B)

  implicit none

  character(len=*), parameter :: NameSub='GM_calc_fac'

  ! Logical for new positions and
  ! the number of mapped FAC points for North and South hemispheres
  logical, intent(out) :: IsNewFacPoint
  integer, intent(out) :: nPoint_I(2)
  real, intent(out)    :: ColatLim_B(2)

  

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_calc_fac

!==============================================================================
subroutine GM_put_from_ie_swmf(nPartial,iPut,Put,W,DoAdd,State_V,nVar)
  use CON_router                           !See srcGM
  implicit none
  integer,intent(in)::nPartial,iPut,nVar
  type(IndexPtrType),intent(in)::Put
  type(WeightPtrType),intent(in)::W
  logical,intent(in)::DoAdd
  real,dimension(nVar),intent(in)::State_V
  call CON_stop(&
       'GM_put_from_ie_swmf: GM_ERROR: empty version cannot be used!')
end subroutine GM_put_from_ie_swmf
!==============================================================================
!subroutine ionosphere_fac(iBlock) !srcIE/iono_coupling.f90
!  implicit none
!  integer,intent(in)::iBlock
!  call CON_stop('ionosphere_fac: GM_ERROR: empty version cannot be used!')
!end subroutine ionosphere_fac
!=============================================================================
subroutine transform_phi_bc_to_u_bc
  implicit none
  call CON_stop(&
       ' transform_phi_bc_to_u_bc: GM_ERROR: empty version cannot be used!')
end subroutine transform_phi_bc_to_u_bc
!=============================================================================
real function logvar_ionosphere()
  call CON_stop('logvar_ionosphere: GM_ERROR: empty version cannot be used!')
  logvar_ionosphere = -777.77
end function logvar_ionosphere
!=============================================================================
subroutine calc_inner_bc_velocities
  call CON_stop('calc_inner_bc_velocities: '// &
       'GM_ERROR: empty version cannot be used!')
end subroutine calc_inner_bc_velocities
!=============================================================================
subroutine GM_calc_iono_bcs

  implicit none

  character(len=*), parameter :: NameSub='GM_calc_iono_bcs'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_calc_iono_bcs

!==============================================================================

subroutine GM_print_variables(NameSource)

  implicit none

  character(len=*), parameter :: NameSub='GM_print_variables'

  character(len=*), intent(in) :: NameSource

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_print_variables

!==============================================================================

subroutine GM_synchronize_refinement(iProc0,iCommUnion)

  implicit none
  integer,intent(in) :: iProc0,iCommUnion
  character(len=*), parameter :: NameSub='GM_synchronize_refinement'
  
  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_synchronize_refinement

!==============================================================================

subroutine GM_put_from_ih(nPartial,iPutStart,Put,Weight,DoAdd,StateSI_V,&
     nVar)

  ! Derived type arguments, it is easier not to declare them
  character(len=*), parameter :: NameSub='GM_put_from_ih'

  call CON_stop(NameSub//': GM_ERROR: empty version cannot be used!')

end subroutine GM_put_from_ih

