!^CMP COPYRIGHT UM
!^CMP FILE IE
!^CMP FILE IM

!BOP
!MODULE: CON_couple_ie_im - couple IE and IM components
!
!DESCRIPTION:
! Couple between two components:\\
!    Ionosphere Electrodynamics (IE) Source\\
!    Inner Magnetosphere (IM)        Target
!INTERFACE:
module CON_couple_ie_im

  !USES:
  use CON_coupler

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_im_init ! initialize coupling
  public :: couple_ie_im      ! couple IE to IM

  !REVISION HISTORY:
  ! 06/26/2003 G.Toth <gtoth@umich.edu> - initial version as external
  !                                       subroutines
  !                  implemented MPI, MCT, and SWMF style couplings
  !                  finally SWMF coupling is retained
  ! 08/27/2003 G.Toth - formed a module out of the external subroutine
  ! 09/02/2003 G.Toth - merged with I.Sokolov's version
  ! 09/15/2003 I.Sokolov - non-uniform ionosphere is added
  !EOP

  type(GridDescriptorType)::IE_Grid           !Source!!
  type(GridDescriptorType)::IM_Grid           !Target!!
  type(RouterType),save:: RouterIeIm
  logical :: IsInitialized = .false.
 
  logical :: DoTest, DoTestMe

  ! Name of this interface
  character (len=*), parameter :: NameSub='couple_ie_im'
contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_im_init - initialize IE-UA coupling
  !INTERFACE:
  subroutine couple_ie_im_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. Since both IE and IM grids are
    ! static, the router is formed here for the whole run.
    !EOP
    !------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.
    
    ! Initialize the coupler including communicator for this router
    ! This must be called by ALL processors due to MPI restrictions!

    call init_coupler(                          &
         iCompSource=IE_,                       &
         nGhostPointSource=1,                   &      
         StandardSource_=Nodes_,                & ! from IE nodes
         iCompTarget=IM_,                       &
         nIndexTarget=2,                        & ! number of indexes for IM: iColat,iLon 
         GridDescriptorSource=IE_Grid,          &
         GridDescriptorTarget=IM_Grid,          &
         Router=RouterIeIm)

    ! It is time to leave for non-involved PEs
    if(.not.RouterIeIm%IsProc) RETURN
    
    ! Both grids are static, it is sufficient to set the router once
    call set_router(&
         IE_Grid,                     &
         IM_Grid,                     &
         RouterIeIm,& ! all blocks (just 1)and all cells on IM 
         mapping=map_im_to_ie,    & ! mapping between IM and IE coords
         interpolate=bilinear_interpolation) ! from IE nodes

  end subroutine couple_ie_im_init

  !BOP =======================================================================
  !IROUTINE: couple_ie_im - couple IE to IM component
  !INTERFACE:
  subroutine couple_ie_im(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ionosphere Electrodynamics (IE) Source\\
    !    Inner Magnetosphere (IM)        Target
    !
    ! Send electrostatic potential from IE to IM.
    !EOP

    external IE_get_for_im,IM_put_from_ie
    integer, parameter :: nVarIeIm=4
    real :: tSimulationTmp
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    
    ! After everything is initialized exclude PEs which are not involved
    if(.not.RouterIeIm%IsProc) RETURN

    ! Make sure that IE provides the most recent results
    if(is_proc(IE_)) then
       ! Use temporary variable for the intent(inout) argument of IE_run.
       tSimulationTmp = tSimulation
       call IE_run(tSimulationTmp, tSimulation)
    endif

    call couple_comp(& 
         RouterIeIm, nVarIeIm,&
         fill_buffer =IE_get_for_im,&
         apply_buffer=IM_put_from_ie)

    if(is_proc(IM_)) call IM_put_from_ie_complete

    if(DoTest.and.is_proc0(IM_)) call IM_print_variables('IE')

  end subroutine couple_ie_im

  !============================================================================

  subroutine map_im_to_ie(&
       IMi_nDim,  &
       IMr1_Xyz_D,&
       IEi_nDim,&
       IEr1_Xyz_D,&
       IsInterfacePoint)
    
    ! Map IM generalized coordinates into IE generalized coordinates
    integer, intent(in) :: IMi_nDim,IEi_nDim
    real, intent(in)    :: IMr1_Xyz_D(IMi_nDim)
    real, intent(out)   :: IEr1_Xyz_D(IEi_nDim)
    logical,intent(out) :: IsInterfacePoint
 
    real :: ColetLon_D(2)
    integer :: iColat, iLon
    !------------------------------------------------------------------------
    IsInterfacePoint=.true.

    iColat = nint(IMr1_Xyz_D(1))
    iLon = nint(IMr1_Xyz_D(2))
    
    if(  iColat<1 .or. iColat>Grid_C(IM_)% nCoord_D(1) .or. &
         iLon<1 .or. iLon>Grid_C(IM_)% nCoord_D(2)       )then
       write(*,*)'map_im_to_ie: IMr1_Xyz_D=',IMr1_Xyz_D
       write(*,*)'map_im_to_ie: iColat,iLon,nCoord_D=',&
            iColat,iLon,Grid_C(IM_) % nCoord_D
       call CON_stop('map_im_to_ie SWMF_ERROR: index out of range!')
    end if
    
    !For structured but non-uniform IM grid:
    call gen_to_stretched(IMr1_Xyz_D, &!in:generalized IM coords(indexes) 
                          ColetLon_D, &!out:stretched coords (radians) 
                          2,          &!IM_grid dimension
                          IM_)         !IM_grid ID

    !For structured but non-uniform ionosphere grid
                     
    call stretched_to_gen(ColetLon_D,&!in:stretched coords (radians)
                          IEr1_Xyz_D,&!out:generalized IE cords(indexes)
                          2,         &!IE_grid dimension
                          IE_)        !IE_grid ID
                            
  end subroutine map_im_to_ie

end module CON_couple_ie_im
