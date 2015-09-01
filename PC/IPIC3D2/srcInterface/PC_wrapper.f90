!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module PC_wrapper

  ! Wrapper for the IPIC3D (PC) component

  implicit none

  private ! except

  public:: PC_set_param
  public:: PC_init_session
  public:: PC_finilize_init_session
  public:: PC_run
  public:: PC_save_restart
  public:: PC_finalize

  ! CON_coupler_points
  public:: PC_get_grid_info
  public:: PC_find_points

  ! GM coupler
  public:: PC_put_from_gm_dt
  public:: PC_put_from_gm_init
  public:: PC_put_from_gm
  public:: PC_get_for_gm

  ! Local variables
  integer :: nDim
  integer, save :: iProc
  real :: DtSi
contains

  !==========================================================================
  subroutine PC_set_param(CompInfo, TypeAction)

    use CON_coupler, ONLY: PC_, set_coord_system, &
         init_decomposition, done_dd_init
    use CON_comp_info
    use ModReadParam

    ! Arguments
    type(CompInfoType), intent(inout):: CompInfo   ! Information for this comp.
    character (len=*), intent(in)    :: TypeAction ! What to do

    ! The PARAM.in segment for IPIC3D
    character(len=lStringLine), allocatable :: StringLineF_I(:) 

    character (len=*), parameter :: NameSub='PC_set_param'
    integer, save :: iComm, nProc

    !-------------------------------------------------------------------------
    select case(TypeAction)
    case('VERSION')
       call put(CompInfo,&
            Use        = .true., &
            NameVersion= 'IPIC3D', &
            Version    = 1.0)
    case('MPI')
       call get(CompInfo, iComm=iComm, iProc=iProc, nProc=nProc)
       call pc_wrapper_c_init_mpi(iComm, iProc, nProc)
    case('READ', 'CHECK')

       ! get section of PARAM.in that contains the PC module
       allocate(StringLineF_I(i_line_read()+1:n_line_read()))
       call read_text(StringLineF_I)

if(iProc==0)then
write(*,*)NameSub,': TypeAction: ',TypeAction
end if

       if(n_line_read()-i_line_read()+1 > 0) then
          call pc_wrapper_c_read_param(StringLineF_I, &
               n_line_read()-i_line_read()+1, lStringLine,iProc)
       end if
    case('STDOUT')

    case('FILEOUT')

    case('GRID')
       ! Grid info depends on IPIC3D (BATSRUS)
       ! Using this constrution to only communicate NameVar and nVar
       if(.not.done_dd_init(PC_))then
          call init_decomposition(PC_,PC_,3,.true.)
          call set_coord_system(PC_,'GSM', 1.0, &
               NameVar = "Variables given by coupled code!!!!", nVar=-1)
       endif
    case default
       call CON_stop(NameSub//': PC_ERROR: unknown Action='//TypeAction)
    end select

  end subroutine PC_set_param


  !============================================================================

  subroutine PC_init_session(iSession, TimeSimulation)

    !INPUT PARAMETERS:
    integer,  intent(in) :: iSession         ! session number (starting from 1)
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_init_session'
    !--------------------------------------------------------------------------

    call pc_wrapper_c_init(TimeSimulation)

  end subroutine PC_init_session

  !============================================================================

  subroutine PC_finilize_init_session

    character(len=*), parameter :: NameSub='PC_finilize_init_session'
    !--------------------------------------------------------------------------

    call pc_wrapper_c_finilize_init

  end subroutine PC_finilize_init_session
  !============================================================================

  subroutine PC_finalize(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_finalize'
    !-------------------------------------------------------------------------

    ! call pc_wrapper_c_end()

  end subroutine PC_finalize

  !============================================================================

  subroutine PC_save_restart(TimeSimulation)

    !INPUT PARAMETERS:
    real,     intent(in) :: TimeSimulation   ! seconds from start time

    character(len=*), parameter :: NameSub='PC_save_restart'
    !--------------------------------------------------------------------------

    call pc_wrapper_c_save_restart()

  end subroutine PC_save_restart

  !============================================================================

  subroutine PC_run(TimeSimulation, TimeSimulationLimit)

    !INPUT/OUTPUT ARGUMENTS:
    real, intent(inout):: TimeSimulation   ! current time of component

    !INPUT ARGUMENTS:
    real, intent(in):: TimeSimulationLimit ! simulation time not to be exceeded

    real:: Dt, Time

    logical:: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='PC_run'
    !--------------------------------------------------------------------------
    call CON_set_do_test(NameSub, DoTest, DoTestMe)

    if(DoTestMe)write(*,*) NameSub, &
         ' starting with TimeSimulation, TimeSimulationLimit=', &
         TimeSimulation, TimeSimulationLimit

    !! Sync timestep with the wrapper
    Dt = TimeSimulationLimit - TimeSimulation
    if(DtSi > Dt .and. Dt >= 0.0 ) then 
        DtSi =  Dt 
        ! set the right time step
        call pc_wrapper_c_set_dt(DtSi)
    end if

    call pc_wrapper_c_run(Time)
    TimeSimulation = TimeSimulation + DtSi    

    if(DoTestMe)write(*,*) NameSub, &
         ' finishing with TimeSimulation, DtSi=', &
         TimeSimulation, DtSi

  end subroutine PC_run

  !============================================================================

  subroutine PC_get_grid_info(nDimOut, iGridOut, iDecompOut)

    integer, intent(out):: nDimOut    ! grid dimensionality
    integer, intent(out):: iGridOut   ! grid index
    integer, intent(out):: iDecompOut ! decomposition index

    logical:: IsFirstTime = .true.

    character(len=*), parameter :: NameSub = 'PC_get_grid_info'
    !--------------------------------------------------------------------------

    !call pc_wrapper_c_get_ndim(nDim)
    nDimOut    = nDim
    if(IsFirstTime)then
       iGridOut   = 1
       iDecompOut = 1
    else
       iGridOut   = 2
       iDecompOut = 2
    end if
    IsFirstTime = .false.

  end subroutine PC_get_grid_info
  !============================================================================
  subroutine PC_put_from_gm_dt(INDtSi)

    real,    intent(in) :: INDtSi

    character(len=*), parameter :: NameSub = 'PC_put_from_gm_dt'
    !--------------------------------------------------------------------------

    ! Store the time step, set it when we do PC_run
    DtSi = INDtSi

    call pc_wrapper_c_set_dt(DtSi)

  end subroutine PC_put_from_gm_dt
  !============================================================================
  subroutine PC_put_from_gm_init(ParamInt_I, ParamReal_I, n, NameVar, nNameVar)

    integer, intent(in) :: n
    integer, intent(in) :: ParamInt_I(3)
    real,    intent(in) :: ParamReal_I(n)
    character(len=*), intent(inout):: NameVar ! List of variables
    integer, intent(in) :: nNameVar

    character(len=*), parameter :: NameSub = 'PC_put_from_gm_init'
    !--------------------------------------------------------------------------
    nDim = ParamInt_I(3)
    call pc_wrapper_c_set_ndim(ParamInt_I(3))
    call pc_wrapper_c_from_gm_init(ParamInt_I, ParamReal_I, n, NameVar, &
         nNameVar) 

  end subroutine PC_put_from_gm_init
  !============================================================================
  subroutine PC_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

    integer, intent(in) :: nDimIn                ! dimension of position vector
    integer, intent(in) :: nPoint                ! number of positions
    real,    intent(in) :: Xyz_DI(nDimIn,nPoint) ! positions
    integer, intent(out):: iProc_I(nPoint)       ! processor owning position

    ! Find array of points and return processor indexes owning them
    ! Could be generalized to return multiple processors...

    character(len=*), parameter:: NameSub = 'PC_find_points'
    !--------------------------------------------------------------------------

    call pc_wrapper_c_find_points(nDimIn, nPoint, Xyz_DI, iProc_I)

  end subroutine PC_find_points
  !============================================================================
  subroutine PC_put_from_gm( &
       NameVar, nVar, nPoint, Data_VI, iPoint_I, Pos_DI)

    character(len=*), intent(inout):: NameVar ! List of variables
    integer,          intent(inout):: nVar    ! Number of variables in Data_VI
    integer,          intent(inout):: nPoint  ! Number of points in Pos_DI
    real,    intent(in), optional  :: Data_VI(:,:)    ! Recv data array
    integer, intent(in), optional  :: iPoint_I(nPoint)! Order of data
    real, intent(out), allocatable, optional:: Pos_DI(:,:)  ! Position vectors

    character(len=*), parameter :: NameSub='PC_put_from_gm'
    !--------------------------------------------------------------------------

    if(.not. present(Data_VI))then
       call pc_wrapper_c_get_ngridpoints(nPoint)

       if(allocated(Pos_DI)) deallocate(Pos_DI)
       allocate(Pos_DI(nDim,nPoint))

       call pc_wrapper_c_get_grid(Pos_DI,nDim*nPoint)
       RETURN
    end if

    call pc_wrapper_c_set_state_var(NameVar//char(0), len(NameVar//char(0)), &
         nVar, Data_VI, iPoint_I)

  end subroutine PC_put_from_gm
  !============================================================================
  subroutine PC_get_for_gm(IsNew, NameVar, nVarIn, nDimIn, nPoint, Xyz_DI, &
       Data_VI)

    !use CON_coupler, ONLY: i_proc, PC_, n_proc

    logical,          intent(in):: IsNew   ! true for new point array
    character(len=*), intent(in):: NameVar ! List of variables
    integer,          intent(in):: nVarIn  ! Number of variables in Data_VI
    integer,          intent(in):: nDimIn  ! Dimensionality of positions
    integer,          intent(in):: nPoint  ! Number of points in Xyz_DI

    real, intent(in) :: Xyz_DI(nDimIn,nPoint)  ! Position vectors
    real, intent(out):: Data_VI(nVarIn,nPoint) ! Data array

    !integer :: i

    character(len=*), parameter :: NameSub='GM_get_for_pc'
    !--------------------------------------------------------------------------

    call pc_wrapper_c_get_state_var( &
         nDimIn, nPoint, Xyz_DI, Data_VI, nVarIn, &
         NameVar//char(0), len(NameVar//char(0)));

!!$    if(i_proc(PC_) == 0) then
!!$       OPEN(UNIT=9, FILE="stat_for_pos_0.txt", ACTION="write") 
!!$       do i=1,nPoint
!!$          write(9,*) Xyz_DI(1,i),", ",Xyz_DI(2,i),", ",Data_VI(1,i),", ",i_proc(PC_) 
!!$       end do
!!$       CLOSE(9)
!!$    end if

  end subroutine PC_get_for_gm

end module PC_wrapper
