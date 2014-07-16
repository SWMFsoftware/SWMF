!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE RB

!BOP
!MODULE: CON_couple_gm_rb - couple GM and RB components
!
!DESCRIPTION:
! Couple GM and RB components one way for now. 
!
!INTERFACE:
module CON_couple_gm_rb

  !USES:
  use CON_coupler
  use CON_transfer_data, ONLY: &
       transfer_integer, transfer_real_array, transfer_string_array

  use GM_wrapper, ONLY: &
       GM_get_for_rb, GM_get_for_rb_trace, GM_get_sat_for_rb, GM_satinit_for_rb
  use RB_wrapper, ONLY: &
       RB_put_from_gm, RB_put_sat_from_gm

  implicit none

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_rb_init ! initialize both couplings
  public :: couple_gm_rb      ! couple GM to RB

  !REVISION HISTORY:
  ! 05/21/2004 O.Volberg - initial version
  !
  !EOP

  ! Communicator and logicals to simplify message passing and execution
  logical :: IsInitialized = .false.

  ! Size of the 2D spherical structured (possibly non-uniform) RB grid
  integer, save :: iSize, jSize

  ! Number of satellites in GM that will also be traced in RB
  integer, save :: nShareSats
contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_rb_init - initialize GM-RB coupling
  !INTERFACE:
  subroutine couple_gm_rb_init

    !DESCRIPTION:
    ! Store RB grid size. Transfer number of shared satellites from GM to RB.
    !EOP
    !------------------------------------------------------------------------
    if(IsInitialized) RETURN
    IsInitialized = .true.

    if(.not.(is_proc(RB_) .or. is_proc(GM_))) RETURN

    iSize = Grid_C(RB_) % nCoord_D(1)
    jSize = Grid_C(RB_) % nCoord_D(2)

    ! Set number of satellites shared between GM and RB for tracing.
    if(is_proc(GM_))call GM_satinit_for_rb(nShareSats)
    call transfer_integer(GM_, RB_, nShareSats, UseSourceRootOnly=.false.)

  end subroutine couple_gm_rb_init

  !BOP =======================================================================
  !IROUTINE: couple_gm_rb - couple GM to RB component
  !INTERFACE: couple_gm_rb(tSimulation)
  subroutine couple_gm_rb(tSimulation)

    !INPUT ARGUMENTS:
    real, intent(in) :: tSimulation     ! simulation time at coupling

    !DESCRIPTION:
    ! Couple between two components:\\
    !    Global Magnetosphere (GM) source\\
    !    Radiation Belt       (RB) target
    !
    ! Send field line volumes, average density and pressure and
    ! geometrical information.
    !EOP

    !\
    ! Coupling variables
    !/

    ! Number of integrals to pass
    integer, parameter :: nIntegral=6

    ! Names of variables to pass
    character (len=*), parameter:: NameVar='x:y:bmin:I_I:S_I:R_I:B_I:rho:p'

    ! Number of variables and points saved into the line data
    ! Probably nVarLine should be a parameter (=4)
    integer :: nVarLine, nPointLine

    ! Buffer for the variables on the 2D RB grid and line data
    real, allocatable :: Integral_IIV(:,:,:), BufferLine_VI(:,:)

    ! Buffer for satellite locations   
    real, allocatable :: SatPos_DII(:,:,:)
    
    ! Buffer for satellite names   
    character (len=100), allocatable:: NameSat_I(:)
    
    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='couple_gm_rb'
    !-------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting, iProc=', i_proc()

    if (nShareSats > 0) then
       ! If RB sat tracing is enabled, get sat locations from GM
       allocate(SatPos_DII(4,2,nShareSats), NameSat_I(nShareSats))
       if(is_proc(GM_)) &
            call GM_get_sat_for_rb(SatPos_DII, NameSat_I, nShareSats)

       ! The satellite trace is only on the root Source
       call transfer_real_array(GM_, RB_, 4*2*nShareSats, SatPos_DII)

       ! The satellite names are known all source procs
       call transfer_string_array(GM_, RB_, nShareSats, NameSat_I, &
            UseSourceRootOnly = .false.)

       if(is_proc(RB_)) &
            call RB_put_sat_from_gm(nShareSats, NameSat_I, SatPos_DII)
       deallocate(SatPos_DII, NameSat_I)
    end if

    ! Allocate buffers both in GM and RB
    allocate(Integral_IIV(iSize,jSize,nIntegral))

    ! Get field line integrals from GM. Only GM root returns the results.
    if(is_proc(GM_)) &
         call GM_get_for_rb_trace(iSize, jSize, NameVar, nVarLine, nPointLine)

    ! Send over size information
    call transfer_integer(GM_, RB_, nVarLine)
    call transfer_integer(GM_, RB_, nPointLine)

    if(is_proc0(GM_) .or. is_proc(RB_))then
       allocate(BufferLine_VI(nVarLine, nPointLine))
    else
       allocate(BufferLine_VI(1,1))
    end if

    if(is_proc0(GM_))call GM_get_for_rb(Integral_IIV, iSize, jSize, nIntegral, &
            BufferLine_VI, nVarLine, nPointLine, NameVar)

    call transfer_real_array(GM_, RB_, size(BufferLine_VI), BufferLine_VI)
    call transfer_real_array(GM_, RB_, size(Integral_IIV), Integral_IIV)
       
    if(is_proc(RB_)) call RB_put_from_gm(Integral_IIV,iSize,jSize,nIntegral,&
         BufferLine_VI,nVarLine,nPointLine,NameVar,tSimulation)

    deallocate(Integral_IIV, BufferLine_VI)

    if(DoTest)write(*,*)NameSub,' finished, iProc=', i_proc()

  end subroutine couple_gm_rb

end module CON_couple_gm_rb
