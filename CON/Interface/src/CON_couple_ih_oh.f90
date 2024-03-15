!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IH
!^CMP FILE OH
module CON_couple_ih_oh

  ! This coupler uses the SWMF parallel coupling toolkit.
  ! The IH grid is coupled to a buffer grid in OH. The buffer grid
  ! uses the same coordinate system as IH, so the transformation is
  ! done in the OH wrapper.
  !
  ! The OH grid is coupled to the outer ghost cells of the IH grid directly.
  ! Both IH and OH use AMR grids, the buffer is a simple spherical grid.

  use CON_coupler
  use CON_transfer_data, ONLY: transfer_real_array, transfer_integer

  use IH_wrapper, ONLY: &
       IH_get_for_global_buffer

  use OH_wrapper, ONLY: &
       OH_match_ibc, &
       OH_set_buffer_grid_get_info, &
       OH_save_global_buffer

  implicit none
  private ! except

  public:: couple_ih_oh_init
  public:: couple_ih_oh
  public:: couple_oh_ih

  ! revision history:
  ! 7/23/03 Sokolov I.V.<igorsok@umich.edu> - prototype for ih-gm
  ! 7/04/04                                 - version for ih-sc
  ! 7/20/04                                 - version for sc-buffer

  logical       :: IsInitialized=.false., DoMatchIBC = .true.

  ! Size and limits of the 3D spherical buffer grid
  integer, save :: iSize, jSize, kSize
  real, save    :: BufferMinMaxIh_DI(3,2)

  character(len=*), parameter :: NameMod='couple_ih_oh'

contains
  !============================================================================

  subroutine couple_ih_oh_init

    ! Couple IH and OH components via a buffer grid
    ! The subroutines:
    !                CON_couple_ih_oh_init
    !                CON_couple_ih_oh
    !                CON_couple_oh_ih (not implemented)

    ! revision history:
    ! 07/25/2003 G.Toth <gtoth@umich.edu> - initial version as external
    !                                       subroutines
    ! 08/27/2003 G.Toth - combined into a module
    ! 12/01/2004 G.Toth - the GM->IE coupling is rewritten for Jr(iSize,jSize)

    ! 12/12/2011 R.Oran <oran@umich.edu> version for two BATSRUS components
    !                                    using a buffer grid

    logical :: DoTest, DoTestMe

    ! This subroutine should be called from all PE-s
    ! Share buffer grid info (set in OH) with IH.

    character(len=*), parameter:: NameSub = 'couple_ih_oh_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(.not.is_proc(IH_) .and. .not.is_proc(OH_)) RETURN

    if(DoTest) write(*,*) NameSub, ' started'

    ! Determine which state variables should be coupled,
    ! pass this info to IH and OH
    call set_couple_var_info(IH_, OH_)

    ! Set buffer grid location and size in OH, and retrieve them for coupler
    if(is_proc(OH_))then
       call OH_set_buffer_grid_get_info( &
            iSize, jSize, kSize, BufferMinMaxIh_DI)

       ! Convert units for radial coordinate
       BufferMinMaxIh_DI(1,:) = BufferMinMaxIh_DI(1,:) &
            *(Grid_C(OH_)%UnitX/Grid_C(IH_)%UnitX)
    end if

    ! Pass buffer size
    call transfer_integer(OH_, IH_, iSize, jSize, kSize, &
         UseSourceRootOnly = .false.)

    ! Pass buffer boundary info
    call transfer_real_array(OH_, IH_, 6, BufferMinMaxIh_DI, &
         UseSourceRootOnly = .false.)

  end subroutine couple_ih_oh_init
  !============================================================================

  subroutine couple_ih_oh(TimeCoupling)

    real, intent(in) :: TimeCoupling     ! simulation time at coupling

    ! Couple between two components:
    !    Inner Heliosphere (IH)  source
    !    Outer Heliosphere (OH)  target
    !
    ! The IH component sends the state variables to a buffer grid.
    ! OH uses the buffer grid to calculate the inner boundary conditions.

    ! Array to store state vector on all buffer grid points
    real, allocatable:: Buffer_VIII(:,:,:,:)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_ih_oh'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    if(DoTest.and.is_proc0(OH_))&
         write(*,'(a,es12.5)')NameSub//': starting, Time=', TimeCoupling

    ! Transfer buffer grid from IH to OH to be used for inner boundary
    allocate(Buffer_VIII(nVarCouple,iSize,jSize,kSize))
    if(is_proc(IH_)) call IH_get_for_global_buffer(iSize, jSize, kSize, &
         BufferMinMaxIh_DI, Buffer_VIII)

    ! Add up Buffer on IH processors and transfer to OH
    call transfer_real_array(IH_, OH_, size(Buffer_VIII), Buffer_VIII, &
         UseSourceSum=.true.)

    if(is_proc(OH_)) call OH_save_global_buffer( &
         nVarCouple, iSize, jSize, kSize, Buffer_VIII)
    deallocate(Buffer_VIII)

    ! Apply initial boundary condition in OH
    if(DoMatchIBC) then
       DoMatchIBC = .false.
       if(is_proc(OH_)) call OH_match_IBC
    end if

    if(DoTest.and.is_proc0(OH_))&
         write(*,'(a,es12.5)')NameSub//': finished, Time=', TimeCoupling

  end subroutine couple_ih_oh
  !============================================================================

  subroutine couple_oh_ih(TimeCoupling)

    real, intent(in) :: TimeCoupling

    ! Couple between two components:
    !    Outer Heliosphere (OH) source
    !    Inner Heliosphere (IH) target
    !
    ! Send state variable from OH to outer cells in IH.

    ! Buffer for state variable to fill outer cells of IH

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_oh_ih'
    !--------------------------------------------------------------------------
    call CON_stop(NameSub// &
         ' is not yet implemented. Correct #COUPLERTYPE command in PARAM.in')
  end subroutine couple_oh_ih
  !============================================================================

end module CON_couple_ih_oh
!==============================================================================

