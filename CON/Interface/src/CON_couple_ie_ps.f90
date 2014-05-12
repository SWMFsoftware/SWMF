!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE IE
!^CMP FILE PS

!BOP
!MODULE: CON_couple_ie_ps - couple IE and PS components
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
  !EOP

  character(len=lNameVersion):: NameVersionPs
  logical :: IsInitialized = .false.

  ! Variables for coupler with coupling toolkit
  type(GridDescriptorType)::IE_Grid           ! Source
  type(GridDescriptorType)::PS_Grid           ! Target
  type(RouterType),save:: RouterIePs, RouterPsIe

  logical :: DoTest, DoTestMe

  ! Variables for the simple coupler
  logical, save :: UseMe
  integer, save :: nTheta, nPhi
  integer, save :: iProc0Ie, iProc0Ps, iCommWorld

  ! Name of this interface
  character (len=*), parameter :: NameMod='CON_couple_ie_ps'

contains

  !BOP =======================================================================
  !IROUTINE: couple_ie_ps_init - initialize IE-PS coupling
  !INTERFACE:
  subroutine couple_ie_ps_init

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s so that
    ! a union group can be formed. Since both IE and PS grids are
    ! static, the router is formed here for the whole run.
    !EOP

    use CON_world, ONLY: get_comp_info
    !------------------------------------------------------------------------

    if(IsInitialized) RETURN
    IsInitialized = .true.

    ! This coupler does not work for RAM, because RAM grid is not on the
    ! ionosphere. 
    call get_comp_info(PS_,NameVersion=NameVersionPs)

    !    if(NameVersionPs(1:3) == 'DGC')then
    ! IE-PS coupling uses MPI
    UseMe = is_proc(IE_) .or. is_proc(PS_)

    nTheta = size(Grid_C(IE_) % Coord1_I)
    nPhi   = size(Grid_C(IE_) % Coord2_I)

    iProc0Ps   = i_proc0(PS_)
    iProc0Ie   = i_proc0(IE_)
    iCommWorld = i_comm()

  end subroutine couple_ie_ps_init
  !BOP =======================================================================
  !IROUTINE: couple_ps_ie - couple PS to IE component
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
    !EOP

    integer, parameter :: nVarPsIe=3
    real :: tSimulationTmp
    character(len=*), parameter:: NameSub = NameMod//'::couple_ps_ie'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    write(*,*) "There is no PS->IE coupling!"

  end subroutine couple_ps_ie

  !BOP =======================================================================
  !IROUTINE: couple_ie_ps - couple IE to PS component
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
    !EOP

    integer, parameter :: nVarIePs=4
    real :: tSimulationTmp
    integer :: iProcWorld
    character(len=*), parameter:: NameSub = NameMod//'::couple_ie_ps'
    !-------------------------------------------------------------------------
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
      integer :: i,j
      integer :: iPS_Size, jPS_Size, iSize, jSize
      real :: iPS(Grid_C(PS_) % nCoord_D(1))
      real :: jPS(Grid_C(PS_) % nCoord_D(2))
      real :: Buffer_II(Grid_C(IE_) % nCoord_D(1),Grid_C(IE_) % nCoord_D(2))

      ! General error code
      integer :: iError
      !------------------------------------------------------------------------

      ! After everything is initialized exclude PEs which are not involved
      if(.not.UseMe) RETURN

      if(DoTest)write(*,*)NameSubSub,', iProc, iProc0Ie, iProc0Ps=', &
           iProcWorld, iProc0Ie, iProc0Ps

      ! Get Plasmasphere grid size
      iPs_Size = Grid_C(PS_) % nCoord_D(1)
      jPs_Size = Grid_C(PS_) % nCoord_D(2)

      ! Get Plasmasphere Grid
      iPS = 90.0 - (Grid_C(PS_) % Coord1_I)
      jPS = (Grid_C(PS_) % Coord2_I)-180.    

      do j=1, jPS_Size
         if (jPS(j) < 0) then
            jPS(j) = jPS(j) + 360
         else 
            jPS(j) = jPS(j)
         endif
      enddo

      ! Allocate buffers both on PS and IE root processors
      allocate(Potential_out(iPS_Size,jPS_Size))

      ! Get Potential from IE, then Bilinear Interpolate
      if(is_proc0(IE_)) &
           call IE_get_for_ps(Buffer_II, iSize, jSize, tSimulation)

      ! Transfer variables from IE to PS
      if(iProc0Ie /= iProc0Ps)then
         if(is_proc0(IE_)) &
              call MPI_send(Potential_out, size(Potential_out), &
              MPI_REAL, iProc0Ps, 1, iCommWorld, iError)
         if(is_proc0(PS_)) &
              call MPI_recv(Potential_out, size(Potential_out), &
              MPI_REAL, iProc0Ie, 1, iCommWorld, iStatus_I, iError)
      end if

      if(DoTest) write(*,* )NameSubSub,', variables transferred iProc:', &
           iProcWorld

      ! Put variables into PS
      if(is_proc0(PS_)) &
           call PS_put_from_ie(iPs_Size, jPs_Size, Potential_out)

      !\
      ! Deallocate buffer to save memory
      !/
      deallocate(Potential_Out)

      if(DoTest)write(*,*)NameSubSub,', variables deallocated',&
           ', iProc:',iProcWorld
      if(DoTest)write(*,*)NameSubSub,' finished, iProc=',iProcWorld

    end subroutine couple_mpi

  end subroutine couple_ie_ps

end module CON_couple_ie_ps
