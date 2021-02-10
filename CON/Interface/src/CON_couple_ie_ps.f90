!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE PS

! BOP
! MODULE: CON_couple_ie_ps - couple IE and PS components
!
!DESCRIPTION:
! Couple between two components:\\
!    Ionosphere Electrodynamics (IE) Source\\
!    Plasmasphere (PS)        Target
!INTERFACE:
module CON_couple_ie_ps

  !USES:
  use CON_coupler

  use IE_wrapper, ONLY: IE_get_for_ps
  use PS_wrapper, ONLY: PS_put_from_ie

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_ie_ps_init ! initialize coupling
  public :: couple_ie_ps      ! couple IE to PS
  public :: couple_ps_ie      ! couple PS to IE

  !REVISION HISTORY:
  ! 01/31/2012 A.Dodger <adodger@umich.edu> - Updated to allow IE-PS coupling
  ! EOP

  character(len=lNameVersion):: NameVersionPs
  logical :: IsInitialized = .false.

  logical :: DoTest, DoTestMe

  ! Variables for the simple coupler
  logical, save :: UseMe
  integer, save :: iProc0Ie, iProc0Ps, iCommWorld

  ! Name of this interface
  character (len=*), parameter :: NameMod='CON_couple_ie_ps'

contains
  !============================================================================

  ! BOP =======================================================================
  ! IROUTINE: couple_ie_ps_init - initialize IE-PS coupling
  !INTERFACE:
  subroutine couple_ie_ps_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. Since both IE and PS grids are
    ! static, the router is formed here for the whole run.
    ! EOP

    use CON_world, ONLY: get_comp_info
    !--------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! IE-PS coupling uses MPI
    UseMe = is_proc(IE_) .or. is_proc(PS_)

    ! Get useful node information:
    iProc0Ps   = i_proc0(PS_)
    iProc0Ie   = i_proc0(IE_)
    iCommWorld = i_comm()

  end subroutine couple_ie_ps_init
  !============================================================================
  ! BOP =======================================================================
  ! IROUTINE: couple_ps_ie - couple PS to IE component
  !INTERFACE:
  subroutine couple_ps_ie(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Plasmasphere (PS)        Source\\
    !    Ionosphere Electrodynamics (IE) Target
    !
    ! Send field-align current from PS to IE.
    ! EOP

    integer, parameter :: nVarPsIe=3
    real :: tSimulationTmp
    character(len=*), parameter:: NameSub = 'couple_ps_ie'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    write(*,*) "There is no PS->IE coupling!"

  end subroutine couple_ps_ie
  !============================================================================

  ! BOP =======================================================================
  ! IROUTINE: couple_ie_ps - couple IE to PS component
  !INTERFACE:
  subroutine couple_ie_ps(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Ionosphere Electrodynamics (IE) Source\\
    !    Plasmasphere (PS)        Target
    !
    ! Send electrostatic potential from IE to PS.
    ! EOP

    integer, parameter :: nVarIePs=4
    real :: tSimulationTmp
    integer :: iProcWorld

    character(len=*), parameter:: NameSub = 'couple_ie_ps'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)
    iProcWorld = i_proc()
    call couple_mpi

    if(DoTest)write(*,*)NameSub,': finished iProc=',iProcWorld

  contains
    !==========================================================================

    subroutine couple_mpi

      use CON_coupler
      character (len=*), parameter :: NameSubSub=NameSub//'.couple_mpi'

      ! Variable to pass is potential on the 2D IE grid
      real, dimension(:,:), allocatable ::Potential_Out
      ! MPI status variable
      integer :: iStatus_I(MPI_STATUS_SIZE)
      integer :: iSize, jSize
      real :: Buffer_II(Grid_C(IE_) % nCoord_D(1),Grid_C(IE_) % nCoord_D(2))

      ! General error code
      integer :: iError
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,', iProc, iProc0Ie, iProc0Ps=', &
           iProcWorld, iProc0Ie, iProc0Ps

      ! Ensure that grids are compatable.  Dynamic transformations forthcoming.
      if(Grid_C(IE_)%TypeCoord /= Grid_C(PS_)%TypeCoord) call CON_stop( &
           NameSub//' ERROR: PS/IE use different coordinate systems.')

      ! Size of IE grid:
      iSize = Grid_C(IE_) % nCoord_D(1)
      jSize = Grid_C(IE_) % nCoord_D(2)

      ! Get Potential from IE:
      if(is_proc0(IE_)) &
           call IE_get_for_ps(Buffer_II, iSize, jSize, tSimulation)

      ! Transfer variables from IE to PS if on different nodes:
      if(iProc0Ie /= iProc0Ps)then
         ! Share electric potential (Buffer_II) from IE to PS:
         if(is_proc0(IE_)) &
              call MPI_send(Buffer_II, size(Buffer_II), &
              MPI_REAL, iProc0Ie, 1, iCommWorld, iError)
         if(is_proc0(PS_)) &
              call MPI_recv(Buffer_II, size(Buffer_II), &
              MPI_REAL, iProc0Ps, 1, iCommWorld, iStatus_I, iError)
      end if

      if(DoTest) write(*,* )NameSubSub,', variables transferred iProc:', &
           iProcWorld

      ! Put variables into PS
      if(is_proc0(PS_)) &
           call PS_put_from_ie(iSize, jSize, Buffer_II)

      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi
    !==========================================================================

  end subroutine couple_ie_ps
  !============================================================================

end module CON_couple_ie_ps
