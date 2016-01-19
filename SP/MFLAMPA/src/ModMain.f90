!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMain

  implicit none

  SAVE

  private ! except

  public:: read_param, initialize, run, finalize
  public:: iComm, iProc, nProc

  !\
  ! Logicals for actions
  !----------------------------------------------------------------------------
  ! run the component
  logical:: DoRun = .true.
  ! restart the run 
  logical:: DoRestart = .false.
  ! perform initialization
  logical:: DoInit = .true.

  !\
  ! MPI information
  !----------------------------------------------------------------------------
  integer:: iComm = -1
  integer:: iProc = -1
  integer:: nProc = -1

  !\
  ! Containers for coordinates and data
  !----------------------------------------------------------------------------
  real, allocatable:: State_VCB(:,:,:,:,:)

contains

  subroutine read_param(TypeAction)
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var, read_line, read_command
    
    character (len=*), intent(in)     :: TypeAction ! What to do  

    ! The name of the command
    character (len=100) :: NameCommand
    character (len=*), parameter :: NameSub='SP_set_parameters'
    !--------------------------------------------------------------------------
!    write(iStdOut,*)NameSub//': CHECK iSession =',i_session_read()
    if(DoInit)then
       DoInit=.false.
    end if

    ! Read the corresponding section of input file
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case('#RESTART')
          DoRestart=.true.
       case('#LINE')
          !          call read_var('XyzLine_D(x_)',XyzLine_D(1))
          !          call read_var('XyzLine_D(y_)',XyzLine_D(2))
          !          call read_var('XyzLine_D(z_)',XyzLine_D(3))
          !          call read_var('RBoundSC'     ,RBoundSC    )
          !          call read_var('RBoundIH'     ,RBoundIH    )

       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#SAVEMHDATA')
          !          call read_var('SaveMhData',SaveMhData)
       case('#DOREADMHDATA')
          !          call read_var('DoReadMhData',DoReadMhData)
       case('#NSTEP')
!          call read_var('nStep',iDataSet)
       case('#TSIMULATION')
!          call read_var('tSimulation',SP_Time)
       case('#PLOT')
          !           call read_var('DnPlot',kfriss)
       case('#VERBOSE')
!          call read_var('iVerbose',iVerbose)
!          if(iVerbose>0)DoWriteAll=.true.
!          if(DoWriteAll.and.iProc==0)&
!               write(*,*)prefix,' Verbose everything'
       case('#NSMOOTH')
!          call read_var('nSmooth',nSmooth)
       case('#TEST')
!          call read_var('DoTest', DoTest)
       case default
          call CON_stop(NameSub//&
               ': Unknown command '&
               //NameCommand)
       end select
    end do
  end subroutine read_param

  !============================================================================

  subroutine initialize
    ! allocate arrays used in this model
    integer:: iError
    character(LEN=*),parameter:: NameSub='initialize'
    !--------------------------------------------------------------------------
    if(.not.allocated(State_VCB))&
         allocate(State_VCB(3,1,1,1,1),stat=iError)
    !    call check_allocate(iError,NameSub//'State_VCB')
    State_VCB(:,1,1,1,1) = (/2.5,0.0,0.0/)
  end subroutine initialize

  !============================================================================

  subroutine run
  end subroutine run

  !============================================================================

  subroutine finalize
  end subroutine finalize

end module ModMain
