!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!^CMP FILE GM
!^CMP FILE PC

!BOP
!MODULE: CON_couple_gm_pc - couple GM and PC components
!
!DESCRIPTION:
! Couple GM and PC components both ways. 
!
!INTERFACE:
module CON_couple_gm_pc

  !USES:
  use CON_coupler

  use CON_couple_points

  use GM_wrapper
  use PC_wrapper

  use ModMpi

  implicit none
  save

  private ! except

  !PUBLIC MEMBER FUNCTIONS:

  public :: couple_gm_pc_init ! initialize both couplings
  public :: couple_gm_pc      ! couple GM to PC
  public :: couple_pc_gm      ! couple PC to GM

  !REVISION HISTORY:
  ! 09/24/2013 G.Toth <gtoth@umich.edu> - initial version
  !EOP

  ! Router communicator info
  type(CouplePointsType) :: CouplerGMtoPC,  CouplerPCtoGM

contains

  !BOP =======================================================================
  !IROUTINE: couple_gm_pc_init - initialize GM-PC couplings
  !INTERFACE:
  subroutine couple_gm_pc_init

    integer:: iError, iCommWorld, iStatus_I(MPI_STATUS_SIZE)
    !integer:: nDimPt

    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='couple_gm_pc_init'

    ! iParam_I gives the information to allocate ParamReal_I
    ! needed for seting up the grid and particle constants
    integer :: iParam_I(4)
    real, pointer, dimension(:) :: ParamReal_I

    integer :: n

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iCommWorld = i_comm()

    ! Initialize point couplers

    CouplerGMtoPC%iCompSource = GM_
    CouplerGMtoPC%iCompTarget = PC_
    CouplerGMtoPC%NameVar     = Grid_C(GM_)%NameVar
    CouplerGMtoPC%nVar        = Grid_C(GM_)%nVar
    CouplerGMtoPC%nDim        = iParam_I(3)

    call couple_points_init(CouplerGMtoPC)

    CouplerPCtoGM%iCompSource = PC_
    CouplerPCtoGM%iCompTarget = GM_
    CouplerPCtoGM%NameVar     = Grid_C(GM_)%NameVar
    CouplerPCtoGM%nVar        = Grid_C(GM_)%nVar
    CouplerPCtoGM%nDim        = iParam_I(3)

    call couple_points_init(CouplerPCtoGM)

    ! Get the integer parameters. The 0 is the size of real params for now.
    if(is_proc(GM_)) call GM_get_for_pc_init(iParam_I, 0)

    if(.not. CouplerGMtoPC%IsSameLayout)then
       if(i_proc0(GM_) /= i_proc0(PC_)) then
          if(is_proc0(GM_)) call MPI_send( &
               iParam_I, 4, MPI_INTEGER, i_proc0(PC_),&
               1001, iCommWorld, iError)
          if(is_proc0(PC_)) call MPI_recv( &
               iParam_I, 4, MPI_INTEGER, i_proc0(GM_),&
               1001, iCommWorld, iStatus_I, iError)
       end if
       if(n_proc(PC_) > 1 .and. is_proc(PC_)) call MPI_bcast( &
            iParam_I, 4, MPI_INTEGER, 0, i_comm(PC_),iError)
    end if

    ! n = number of species *3 (mass, charge, temrature ratio)
    !     + numerb of dimentions * 9 ( xmin, xmax, dx ) for y and z also
    !     + 3 ( Lnorm, Unorm, Mnorm) 
    n = iParam_I(1)*3 + iParam_I(2)*9 + 3
    allocate(ParamReal_I(n))

    ! Transfer real parameters from GM to PC
    if(is_proc(GM_)) &
         call GM_get_for_pc_init(iParam_I, n, ParamReal_I)

    if(.not. CouplerGMtoPC%IsSameLayout)then
       if(i_proc0(GM_) /= i_proc0(PC_)) then
          if(is_proc0(GM_)) call MPI_send( &
               ParamReal_I, n, MPI_REAL, i_proc0(PC_),&
               1002, iCommWorld, iError)

          if(is_proc0(PC_)) call MPI_recv( &
               ParamReal_I, n, MPI_REAL, i_proc0(GM_),&
               1002, iCommWorld, iStatus_I, iError)

       end if
       if(n_proc(PC_) > 1 .and. is_proc(PC_)) call MPI_bcast( &
            ParamReal_I, n, MPI_REAL, 0, i_comm(PC_),iError)
    end if

    if(is_proc(PC_)) &
         call PC_put_from_gm_init(iParam_I, ParamReal_I, n)

  end subroutine couple_gm_pc_init

  !=======================================================================
  subroutine couple_gm_pc(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    ! Grid index

    integer:: iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pc'

    integer :: iStatus_I(MPI_STATUS_SIZE)

    ! Time step in SI units. Keeps changing 
    real    :: DtSi
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=', i_proc()

    call couple_points(CouplerGMtoPC, GM_get_grid_info, GM_find_points, &
         GM_get_for_pc, PC_get_grid_info, PC_put_from_gm)

    ! Transfer time step from GM to PC
    if(CouplerGMtoPC%IsSameLayout)then
       call GM_get_for_pc_dt(DtSi)
       call PC_put_from_gm_dt(DtSi)
    else
       if(is_proc0(GM_)) then
          call GM_get_for_pc_dt(DtSi)

          if(i_proc0(GM_) /= i_proc0(PC_))then
             if(is_proc0(GM_)) call MPI_send( &
                  DtSi, 1, MPI_REAL, i_proc0(PC_), &
                  1003, i_comm(), iError)
             if(is_proc0(PC_)) call MPI_recv( &
                  DtSi, 1, MPI_REAL, i_proc0(GM_), &
                  1003, i_comm(), iStatus_I, iError)
          end if

       end if

       if(is_proc(PC_))then
          if(n_proc(PC_) > 1) call MPI_bcast( &
               DtSi, 1, MPI_REAL, 0, i_comm(PC_),iError)
          call PC_put_from_gm_dt(DtSi)
       end if
    end if

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_gm_pc
  !=======================================================================
  subroutine couple_pc_gm(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    ! The first time there is no back coupling
    ! This is actually wrong for restarts!!!
    logical:: IsFirstTime = .true.

    logical :: DoTest, DoTestMe
    character (len=*), parameter :: NameSub='couple_pc_gm'
    !-------------------------------------------------------------------------

    if (IsFirstTime)  then
       IsFirstTime = .false.
       ! Finnishing the setup of the IPIC3D solver
       ! after GM -> PC coupling
       call  PC_finilize_init_session
    end if

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerPCtoGM%iProcWorld

    call couple_points(CouplerPCtoGM, PC_get_grid_info,  PC_find_points , &
         PC_get_for_gm, GM_get_grid_info, GM_put_from_pc)

    if(DoTest) write(*,*) NameSub,' finished, iProc=', i_proc()

  end subroutine couple_pc_gm

end module CON_couple_gm_pc
