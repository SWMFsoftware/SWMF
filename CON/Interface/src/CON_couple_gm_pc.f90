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

  ! Communicator and logicals to simplify message passing and execution
  logical       :: UseMe = .true.

  ! Check if the GM and PT processors coincide
  logical:: IsSameLayout

  ! Proc index inside SWMF
  integer:: iProcWorld

  ! Router communicator info
  integer:: iCommGmPt, nProcGmPt, iProcGmPt, iProc0Gm, iProc0Pt, nProcCommon

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

    ! ParamInt_I gives the information to allocate ParamReal_I
    ! needed for seting up the grid and particle constants
    integer :: ParamInt_I(4)
    real, pointer, dimension(:) :: ParamReal_I

    integer :: n

    !DESCRIPTION:
    ! This subroutine should be called from all PE-s
    !EOP
    !------------------------------------------------------------------------
    call CON_set_do_test(NameSub,DoTest,DoTestMe)

    iCommWorld = i_comm()

    ParamInt_I = -1   
    n = 0
 
    if(is_proc0(GM_)) then
       call GM_get_for_pc_init(ParamInt_I,n)!,ParamReal_I,-1)
       call MPI_send(ParamInt_I, 4, MPI_INTEGER, i_proc0(PC_),&
            1001, iCommWorld, iError)

       !n = nSpecis + nRegions*9 + normalization(3)
       n = ParamInt_I(1)*3 + ParamInt_I(2)*9 + 3
       allocate(ParamReal_I(n))
       !if(.not. allocated(ParamReal_I) ) allocate(ParamReal_I(n))
       call GM_get_for_pc_init(ParamInt_I, n, ParamReal_I)
       call MPI_send(ParamReal_I, n, MPI_DOUBLE, i_proc0(PC_),&
            1002, iCommWorld, iError)
    end if

    if(is_proc(PC_)) then
       if (is_proc0(PC_)) then
          call MPI_recv(ParamInt_I, 4, MPI_INTEGER, i_proc0(GM_),&
               1001, iCommWorld, iStatus_I, iError)
       end if

       call MPI_bcast(ParamInt_I, 4, MPI_INTEGER, 0, i_comm(PC_),iError)


       !n = nSpecis + nRegions * 6
       n = ParamInt_I(1)*3 + ParamInt_I(2)*9 + 3
        allocate(ParamReal_I(n))

       if (is_proc0(PC_)) then
          call MPI_recv(ParamReal_I, n, MPI_DOUBLE, i_proc0(GM_),&
               1002, iCommWorld, iStatus_I, iError)
       end if

       call MPI_bcast(ParamReal_I, n, MPI_DOUBLE, 0, i_comm(PC_),iError)

       call PC_put_from_gm_init(ParamInt_I, ParamReal_I, n)

    end if

    CouplerGMtoPC%iCompTarget = PC_
    CouplerGMtoPC%iCompSource = GM_
    CouplerGMtoPC%nDim        = ParamInt_I(3)

    call couple_points_init(CouplerGMtoPC)

    CouplerPCtoGM%iCompTarget = GM_
    CouplerPCtoGM%iCompSource = PC_
    CouplerPCtoGM%nDim        = ParamInt_I(3)

    call couple_points_init(CouplerPCtoGM)

  end subroutine couple_gm_pc_init

!=======================================================================
  subroutine couple_gm_pc(tSimulation)

    !INPUT ARGUMENT:
    real, intent(in) :: tSimulation

    ! List of variables to pass
    character(len=lNameVar):: NameVar

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1

    integer:: iPoint, i, iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_gm_pc'

    ! Variables for general-case coupling

    ! number of processors of the OTHER component to communicate with
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)

    integer :: iStatus_I(MPI_STATUS_SIZE)
    real    :: SIDt
    !-------------------------------------------------------------------------

    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    CouplerGMtoPC%NameVar = Grid_C(CouplerGMtoPC%iCompSource)%NameVar
    CouplerGMtoPC%nVar    = Grid_C(CouplerGMtoPC%iCompSource)%nVar

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerGMtoPC%iProcWorld



    call couple_points(CouplerGMtoPC, GM_get_grid_info,  GM_find_points, &
         GM_get_for_pc, PC_get_grid_info, PC_put_from_gm)

    ! old argument list
    !call couple_points(CouplerGMtoPC,GM_get_for_pc, PC_put_from_gm, &
    !      PC_get_grid_info, GM_get_grid_info,  GM_find_points)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',CouplerGMtoPC%iProcWorld
 
    if(is_proc0(GM_)) then
      call GM_get_for_pc_dt(SIDt)
      call MPI_send(SIDt, 1, MPI_DOUBLE, i_proc0(PC_),&
                    1003, i_comm(), iError)
    end if

   if(is_proc(PC_)) then
       if (is_proc0(PC_)) then
          call MPI_recv(SIDt, 1, MPI_DOUBLE, i_proc0(GM_),&
               1003, i_comm(), iStatus_I, iError)
       end if
       call MPI_bcast(SIDt, 1, MPI_DOUBLE, 0, i_comm(PC_),iError)
       call PC_put_from_gm_dt(SIDt)
   end if

  end subroutine couple_gm_pc
!=======================================================================
  subroutine couple_pc_gm(tSimulation)

    ! List of variables to pass
    real, intent(in) :: tSimulation

    character(len=lNameVar):: NameVar

    ! Grid index
    integer:: iDecompLastGm = -1, iDecompLastPt = -1

    integer:: iPoint, i, iError

    logical :: DoTest, DoTestMe

    ! Name of this interface
    character (len=*), parameter :: NameSub='couple_pc_gm'

    ! Variables for general-case coupling

    ! number of processors of the OTHER component to communicate with
    integer, save:: nCoupleGm, nCouplePt

    ! processors of the OTHER component to communicate with
    integer, allocatable, save:: iCoupleProcGm_I(:), iCoupleProcPt_I(:)

    ! number of entries received/sent by a processor during rendezvous
    integer, allocatable, save:: nCouplePointGm_I(:), nCouplePointPt_I(:)

    integer :: iStatus_I(MPI_STATUS_SIZE)
    real    :: SIDt
    logical, save :: isFirstTime = .true.

    !-------------------------------------------------------------------------

    CouplerPCtoGM%NameVar = Grid_C(CouplerPCtoGM%iCompSource)%NameVar
    CouplerPCtoGM%nVar    = Grid_C(CouplerPCtoGM%iCompSource)%nVar

    if (isFirstTime)  then
      isFirstTime = .false.
      RETURN
    end if


    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTest)write(*,*)NameSub,' starting iProc=',CouplerPCtoGM%iProcWorld

    call couple_points(CouplerPCtoGM, PC_get_grid_info,  PC_find_points , &
                    PC_get_for_gm, GM_get_grid_info, GM_put_from_pc)

    ! old arguments list
    !call couple_points(CouplerPCtoGM,PC_get_for_gm, GM_put_from_pc, &
    !          PC_get_grid_info, GM_get_grid_info,  PC_find_points)

    if(DoTest) write(*,*) NameSub,' finished, iProc=',CouplerPCtoGM%iProcWorld
 
  end subroutine couple_pc_gm

end module CON_couple_gm_pc
