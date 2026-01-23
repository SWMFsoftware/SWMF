!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE IM

!
! Couple between two components:
!    Ionosphere Electrodynamics (IE) Source
!    Inner Magnetosphere (IM)        Target
module CON_couple_ie_im

  use CON_coupler
  use CON_world, ONLY: use_comp
  use CON_comp_param, ONLY: UA_
  use CON_transfer_data, ONLY: transfer_real_array

  use IE_wrapper, ONLY: IE_get_for_im, IE_put_from_im, IE_put_from_im_complete,&
                        IE_get_info_for_im

  use IE_wrapper, ONLY: IE_get_for_gm ! this is used by RAM-SCB coupler ???
  use IE_wrapper, ONLY: IE_run        ! forces IE and IM run sequentially???

  use IM_wrapper, ONLY: IM_get_info_for_ie, IM_get_for_ie, IM_put_from_ie_mpi, &
       IM_put_from_ie, IM_put_from_ie_complete

  implicit none

  private ! except

  public :: couple_ie_im_init ! initialize coupling
  public :: couple_ie_im      ! couple IE to IM
  public :: couple_im_ie      ! couple IE to IM

  ! revision history:
  ! 06/26/2003 G.Toth <gtoth@umich.edu> - initial version as external
  !                                       subroutines
  !                  implemented MPI, MCT, and SWMF style couplings
  !                  finally SWMF coupling is retained
  ! 08/27/2003 G.Toth - formed a module out of the external subroutine
  ! 09/02/2003 G.Toth - merged with I.Sokolov's version
  ! 09/15/2003 I.Sokolov - non-uniform ionosphere is added

  character(len=lNameVersion):: NameVersionIm
  logical :: IsInitialized = .false.

  ! Variables for coupler with coupling toolkit
  type(GridType)::IE_Grid           ! Source
  type(GridType)::IM_Grid           ! Target
  type(RouterType),save:: RouterIeIm, RouterImIe
  type(LocalGridType) :: IE_LocalGrid, IM_LocalGrid
  logical :: DoTest, DoTestMe

  ! Variables for the simple coupler
  integer, save :: nTheta, nPhi
  integer, save :: nVarImIe=3


  ! Name of this interface
  character (len=*), parameter :: NameMod='CON_couple_ie_im'

contains
  !============================================================================

  subroutine couple_ie_im_init

    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. Since both IE and IM grids are
    ! static, the router is formed here for the whole run.

    use CON_world, ONLY: get_comp_info

    integer :: nEngIM
    !--------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! This coupler does not work for RAM, because RAM grid is not on the
    ! ionosphere.
    call get_comp_info(IM_,NameVersion=NameVersionIm)
    if(NameVersionIm(1:3) == 'RAM')then
       ! IE-IM/RAM coupling uses MPI
       nTheta = size(Grid_C(IE_) % Coord1_I)
       nPhi   = size(Grid_C(IE_) % Coord2_I)
    else
       ! Get extra info from IE when using CIMI
       if(NameVersionIm(1:3) == 'CIM')then
         ! IE-IM/CIMI coupling is being updated
         call IM_get_info_for_ie(nEngIM)
         call IE_get_info_for_im(use_comp(UA_), nEngIM, nVarImIe)
      end if

       ! IE-IM/RCM coupling uses the coupling toolkit
       call init_coupler(                 &
            iCompSource=IE_,              &
            nGhostPointSource=1,          &
            StandardSource_=Nodes_,       & ! from IE nodes
            iCompTarget=IM_,              &
            nIndexTarget=2,               & ! IM grid size: iColat,iLon
            GridSource=IE_Grid, &
            GridTarget=IM_Grid, &
            LocalGridTarget=IM_LocalGrid, &
            Router=RouterIeIm)

       ! It is time to leave for non-involved PEs
       if(RouterIeIm%IsProc) call set_router( &
            IE_Grid,             &
            IM_LocalGrid,                &
            RouterIeIm,             &
            mapping=map_im_to_ie,   &
            interpolate=bilinear_interpolation)

       call init_coupler(                 &
            iCompSource=IM_,              &
            nGhostPointSource=1,          &
            StandardSource_=Nodes_,       & ! from IM nodes
            iCompTarget=IE_,              &
            nIndexTarget=2,               & ! IE grid size: iColat,iLon
            GridSource=IM_Grid, &
            GridTarget=IE_Grid, &
            LocalGridTarget=IE_LocalGrid, &
            Router=RouterImIe)
       ! Both grids are static, it is sufficient to set the router once
       if(RouterImIe%IsProc) call set_router(  &
            IM_Grid,                &
            IE_LocalGrid,           &
            RouterImIe,             &
            mapping=map_ie_to_im,   &
            interpolate=bilinear_interpolation)
    end if

  end subroutine couple_ie_im_init
  !============================================================================
  subroutine couple_im_ie(tSimulation)

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Inner Magnetosphere (IM)        Source
    !    Ionosphere Electrodynamics (IE) Target
    !
    ! Send field-align current from IM to IE.

    character(len=*), parameter:: NameSub = 'couple_im_ie'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! After everything is initialized exclude PEs which are not involved
    if(.not.RouterImIe%IsProc) RETURN

    call couple_comp(&
         RouterImIe, nVarImIe, &
         fill_buffer =IM_get_for_ie,&
         apply_buffer=IE_put_from_im)

    if(is_proc(IE_)) call IE_put_from_im_complete

  end subroutine couple_im_ie
  !============================================================================

  subroutine couple_ie_im(tSimulation)

    real, intent(in) :: tSimulation     ! simulation time at coupling

    ! Couple between two components:
    !    Ionosphere Electrodynamics (IE) Source
    !    Inner Magnetosphere (IM)        Target
    !
    ! Send electrostatic potential from IE to IM.

    integer, parameter :: nVarIeIm=4
    real :: tSimulationTmp
    character(len=*), parameter:: NameSub = 'couple_ie_im'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(NameVersionIm(1:3) == 'RAM')then
       call couple_mpi
    else
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
    end if

  contains
    !==========================================================================

    subroutine couple_mpi

      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Variable to pass potential on the 2D IE grid
      real, allocatable :: Buffer_IIV(:,:,:)

      integer, parameter:: nVar = 1
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(DoTest)write(*,*)NameSubSub,', iProc=', i_Proc()

      ! Allocate buffer both on IM and IE processors. IE sends two variables.
      allocate(Buffer_IIV(nTheta,nPhi,nVar))
      ! Need to call ALL IE processors here
      if(is_proc(IE_)) call IE_get_for_gm( &
           Buffer_IIV, nTheta, nPhi, nVar, ['potential'], tSimulation)

      ! IM/RAM wants only potential. The result is on IE root processor only.
      call transfer_real_array(IE_, IM_, nTheta*nPhi, Buffer_IIV)

      if(is_proc(IM_)) &
           call IM_put_from_ie_mpi(nTheta, nPhi, Buffer_IIV)
      deallocate(Buffer_IIV)

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',i_proc()

    end subroutine couple_mpi
    !==========================================================================

  end subroutine couple_ie_im
  !============================================================================

  subroutine map_ie_to_im(&
       IEi_nDim,  &
       IEr1_Xyz_D,&
       IMi_nDim,&
       IMr1_Xyz_D,&
       IsInterfacePoint)

    ! Map IE generalized coordinates into IM generalized coordinates
    integer, intent(in) :: IEi_nDim,IMi_nDim
    real, intent(in)    :: IEr1_Xyz_D(IEi_nDim)
    real, intent(out)   :: IMr1_Xyz_D(IMi_nDim)
    logical,intent(out) :: IsInterfacePoint

    real :: ColatLon_D(2)
    integer :: iColat, iLon

    character(len=*), parameter:: NameSub = 'map_ie_to_im'
    !--------------------------------------------------------------------------
    IsInterfacePoint=.true.

    iCoLat = nint(IEr1_Xyz_D(1))
    iLon   = nint(IEr1_Xyz_D(2))

    if ( iLon<1 .or. iLon>Grid_C(IE_)% nCoord_D(2) .or. &
         iCoLat<1 .or. iCoLat>Grid_C(IE_)% nCoord_D(1)) then
       write(*,*)'map_ie_to_im: iColat,Grid_C(IE_)% nCoord_D(1)=',&
            iColat, Grid_C(IE_)% nCoord_D(1)
       write(*,*)'map_ie_to_im: iLon,Grid_C(IE_)% nCoord_D(2)=',&
            iLon,Grid_C(IE_)% nCoord_D(2)
       call CON_stop(NameSub//' SWMF_ERROR: index out of range!')
    end if

    ! For structured but non-uniform IM grid:
    call gen_to_stretched(IEr1_Xyz_D, &! in:generalized IE coords(indexes)
                          ColatLon_D, &! out:stretched coords (radians)
                          2,          &! IE_grid dimension
                          IE_)         ! IE_grid ID

    ! For structured but non-uniform ionosphere grid

    call stretched_to_gen(ColatLon_D,&! in:stretched coords (radians)
                          IMr1_Xyz_D,&! out:generalized IM cords(indexes)
                          2,         &! IM_grid dimension
                          IM_)        ! IM_grid ID

  end subroutine map_ie_to_im
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

    real :: ColatLon_D(2)
    integer :: iColat, iLon
    character(len=*), parameter:: NameSub = 'map_im_to_ie'
    !--------------------------------------------------------------------------
    IsInterfacePoint=.true.

    iColat = nint(IMr1_Xyz_D(1))
    iLon   = nint(IMr1_Xyz_D(2))
    if (iLon == 0) iLon = Grid_C(IM_)% nCoord_D(2)

    if(  iColat<1 .or. iColat>Grid_C(IM_)% nCoord_D(1) .or. &
         iLon<1 .or. iLon>Grid_C(IM_)% nCoord_D(2)       )then
       write(*,*)'map_im_to_ie: IMr1_Xyz_D=',IMr1_Xyz_D
       write(*,*)'map_im_to_ie: iColat,iLon,nCoord_D=',&
            iColat,iLon,Grid_C(IM_) % nCoord_D
       call CON_stop(NameSub//' SWMF_ERROR: index out of range!')
    end if

    ! For structured but non-uniform IM grid:
    call gen_to_stretched(IMr1_Xyz_D, &! in:generalized IM coords(indexes)
         ColatLon_D, &! out:stretched coords (radians)
         2,          &! IM_grid dimension
         IM_)         ! IM_grid ID

    ! For structured but non-uniform ionosphere grid

    call stretched_to_gen(ColatLon_D,&! in:stretched coords (radians)
                          IEr1_Xyz_D,&! out:generalized IE cords(indexes)
                          2,         &! IE_grid dimension
                          IE_)        ! IE_grid ID

  end subroutine map_im_to_ie
  !============================================================================

end module CON_couple_ie_im
!==============================================================================
