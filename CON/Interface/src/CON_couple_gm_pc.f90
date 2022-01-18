!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PC

!
! Couple GM and PC components both ways.
!
module CON_couple_gm_pc

  use CON_coupler

  use CON_couple_points

  use GM_wrapper
  use PC_wrapper

  implicit none
  save

  private ! except

  public :: couple_gm_pc_init ! initialize both couplings
  public :: couple_gm_pc      ! couple GM to PC
  public :: couple_pc_gm      ! couple PC to GM
  public :: couple_gm_pc_grid_info ! send grid information from GM to PC

  ! revision history:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version

  ! Router communicator info
  type(CouplePointsType) :: CouplerGMtoPC,  CouplerPCtoGM

contains
  !============================================================================

  subroutine couple_gm_pc_init

    use CON_transfer_data, ONLY: &
         transfer_integer, transfer_integer_array, transfer_real_array

    ! integer:: nDimPt

    logical :: DoTest, DoTestMe

    ! GM sends PC a number of integers and reals
    integer, allocatable :: iParam_I(:)
    real, allocatable :: Param_I(:)
    integer :: nParamInt, nParamReal

    ! This subroutine should be called from all PE-s
    character(len=*), parameter:: NameSub = 'couple_gm_pc_init'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    ! Initialize point couplers

    CouplerGMtoPC%iCompSource = GM_
    CouplerGMtoPC%iCompTarget = PC_
    CouplerGMtoPC%NameVar     = Grid_C(GM_)%NameVar
    CouplerGMtoPC%nVar        = Grid_C(GM_)%nVar

    call couple_points_init(CouplerGMtoPC)

    CouplerPCtoGM%iCompSource = PC_
    CouplerPCtoGM%iCompTarget = GM_
    CouplerPCtoGM%NameVar     = Grid_C(GM_)%NameVar
    CouplerPCtoGM%nVar        = Grid_C(GM_)%nVar

    call couple_points_init(CouplerPCtoGM)

    if(.not.(is_proc(GM_) .or. is_proc(PC_))) RETURN

    ! Get the number of integer and real parameters to pass
    if(is_proc(GM_))call GM_get_for_pc_init(nParamInt, nParamReal)
     call transfer_integer(GM_, PC_, nParamInt, nParamReal, &
         UseSourceRootOnly=.false., UseTargetRootOnly=.false.)
    allocate(iParam_I(nParamInt), Param_I(nParamReal))

    ! Transfer integer and real parameters from GM to PC
    if(is_proc(GM_)) &
         call GM_get_for_pc_init(nParamInt, nParamReal, iParam_I, Param_I)

    call transfer_integer_array(GM_, PC_, nParamInt, iParam_I, &
         UseSourceRootOnly=.false., UseTargetRootOnly=.false.)

    call transfer_real_array(GM_, PC_, nParamReal, Param_I, &
         UseSourceRootOnly=.false., UseTargetRootOnly=.false.)

    if(is_proc(PC_)) &
         call PC_put_from_gm_init(nParamInt, nParamReal, iParam_I, Param_I, &
         CouplerGMtoPC%NameVar)

    deallocate(iParam_I, Param_I)

  end subroutine couple_gm_pc_init
  !============================================================================

  subroutine couple_gm_pc_grid_info()
    use CON_transfer_data, ONLY: &
         transfer_integer, transfer_integer_array

    logical :: DoTest, DoTestMe

    ! GM sends PC a number of integers and reals
    integer, allocatable :: Int_I(:), nSize_I(:)
    integer :: nInt, nPicGrid

    ! This subroutine should be called from all PE-s
    character(len=*), parameter:: NameSub = 'couple_gm_pc_grid_info'
    !--------------------------------------------------------------------------

    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    if(.not.(is_proc(GM_) .or. is_proc(PC_))) RETURN

    ! Get the number of integers to pass
    if(is_proc(GM_))call GM_get_for_pc_grid_info(nInt, nPicGrid)
    call transfer_integer(GM_, PC_, nInt, UseSourceRootOnly=.false., &
         UseTargetRootOnly=.false.)
    call transfer_integer(GM_, PC_, nPicGrid, UseSourceRootOnly=.false., &
         UseTargetRootOnly=.false.)

    if(nInt>0) then
       allocate(Int_I(nInt))
       allocate(nSize_I(nPicGrid))

       if(is_proc(GM_)) &
            call GM_get_for_pc_grid_info(nInt, nPicGrid, &
            nSize_I, Int_I)

       ! Transfer integers from GM to PC
       call transfer_integer_array(GM_, PC_, nInt, Int_I, &
            UseSourceRootOnly=.false., UseTargetRootOnly=.false.)
       call transfer_integer_array(GM_, PC_, nPicGrid, nSize_I, &
            UseSourceRootOnly=.false., UseTargetRootOnly=.false.)

       if(is_proc(PC_)) &
            call PC_put_from_gm_grid_info(nInt, nPicGrid, nSize_I, Int_I)

       deallocate(Int_I)
       deallocate(nSize_I)
    endif

  end subroutine couple_gm_pc_grid_info
  !============================================================================

  subroutine couple_gm_pc(tSimulation)

    use CON_transfer_data, ONLY: transfer_real

    real, intent(in) :: tSimulation

    ! BATSRUS time step in SI units. Keeps changing.
    real:: DtSi

    logical :: DoTest, DoTestMe

    ! Name of this interface

    character(len=*), parameter:: NameSub = 'couple_gm_pc'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    call couple_gm_pc_grid_info

    if(DoTest)write(*,*)NameSub,' starting iProc=', i_proc()

    if(IsTightCouple_CC(GM_,PC_)) then
       if(is_proc(GM_)) call GM_get_for_pc_dt(DtSi)
       call transfer_real(GM_,PC_,DtSi)
       if(is_proc(PC_)) call PC_put_from_gm_dt(DtSi)
    else
       if(is_proc(PC_)) call PC_put_from_gm_dt(Couple_CC(GM_,PC_)%Dt)
    endif

    call couple_points(CouplerGMtoPC, GM_get_grid_info, GM_find_points, &
         GM_get_for_pc, PC_get_grid_info, PC_put_from_gm)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_gm_pc
  !============================================================================
  subroutine couple_pc_gm(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'couple_pc_gm'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerPCtoGM%iProcWorld

    call couple_points(CouplerPCtoGM, PC_get_grid_info,  PC_find_points , &
         PC_get_for_gm, GM_get_grid_info, GM_put_from_pc)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_pc_gm
  !============================================================================

end module CON_couple_gm_pc
!==============================================================================
