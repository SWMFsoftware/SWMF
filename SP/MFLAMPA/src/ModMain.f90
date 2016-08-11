!  Copyright (C) 2002 Regents of the University of Michigan
! portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModMain

  use ModSize, ONLY: &
       nDim, nVar, nLat, nLon, nNode, &
       RMin, iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_, R_, Lat_, Lon_, Bx_, By_, Bz_
  
  use ModWrite, ONLY: &
       set_write_param, write_output, NamePlotDir

  use ModGrid, ONLY: &
       iComm, iProc, nProc, nBlock, &
       Proc_, Block_, Begin_, End_,&
       LatMin, LatMax, LonMin, LonMax, &
       iGrid_IA, iNode_II, iNode_B, State_VIB, CoordOrigin_DA, &
       set_grid_param, init_grid, get_node_indexes
  
  implicit none

  SAVE

  private ! except

  ! Methods and variables from this module 
  public:: &
       read_param, initialize, advance, finalize, check,&
       TimeGlobal, iIterGlobal

  ! Methods and variables from ModSize
  public:: &
       nDim, nVar, nLat, nLon, nNode, &
       RMin, iParticleMin, iParticleMax, nParticle,&
       Particle_, OriginLat_, OriginLon_, R_, Lat_, Lon_, Bx_, By_, Bz_

  ! Methods and variables from ModGrid
  public:: &
       iComm, iProc, nProc, nBlock, &
       Proc_, Block_, Begin_, End_,&
       TypeCoordSystem, LatMin, LatMax, LonMin, LonMax, &
       iGrid_IA, iNode_II, iNode_B, State_VIB, CoordOrigin_DA, &
       get_node_indexes

  ! Methods and variables from ModWrite

  ! Coordinate system
  character(len=3) :: TypeCoordSystem = 'HGI'

  !\
  ! Logicals for actions
  !----------------------------------------------------------------------------
  ! run the component
  logical:: DoRun = .true.
  ! restart the run 
  logical:: DoRestart = .false.
  ! perform initialization
  logical:: DoInit = .true.
  !/

  !\
  ! Global interation and time
  !-----------------------------
  real   :: TimeGlobal  = -1.0
  integer:: iIterGlobal = -1
  !/

contains

  subroutine read_param(TypeAction)
    ! Read input parameters for SP component
    use ModReadParam, ONLY: read_var, read_line, read_command
    character (len=*), intent(in)     :: TypeAction ! What to do  

    ! The name of the command
    character (len=100) :: NameCommand
    character (len=100) :: StringPlot
    character (len=*), parameter :: NameSub='SP:read_param'
    !--------------------------------------------------------------------------
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
       case('#GRID')
          call set_grid_param
       case('#DORUN')
          call read_var('DoRun',DoRun)
       case('#SAVEPLOT')
          call read_var('StringPlot', StringPlot)
          call set_write_param(StringPlot)
       case("#COORDSYSTEM")
          call read_var('TypeCoordSystem',TypeCoordSystem,IsUpperCase=.true.)
       case default
          call CON_stop(NameSub//': Unknown command '//NameCommand)
       end select
    end do
  end subroutine read_param

  !============================================================================

  subroutine initialize
    ! initialize the model
    integer:: iError
    !
    integer:: iLat, iLon, iNode, iBlock, iProcNode
    character(LEN=*),parameter:: NameSub='SP:initialize'
    !--------------------------------------------------------------------------
    iIterGlobal = 0
    call init_grid
  end subroutine initialize

  subroutine advance
    iIterGlobal = iIterGlobal + 1
    call write_output(TimeGlobal, iIterGlobal)
  end subroutine advance

  !============================================================================

  subroutine check
    use ModUtilities, ONLY: make_dir
    character(LEN=*),parameter:: NameSub='SP:check'
    !--------------------------------------------------------------------------
    ! Make output and check input directories
    if(iProc==0) call make_dir(NamePlotDir)
  end subroutine check

  !============================================================================

  subroutine finalize
    character(LEN=*),parameter:: NameSub='SP:finalize'
    !--------------------------------------------------------------------------
  end subroutine finalize

end module ModMain
