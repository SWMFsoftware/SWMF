 !  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
module ModPlotBox

  use ModIO, ONLY: plot_dx, plot_range, plot_normal, TypeCoordPlot_I, plot_form, &
       TypeFile_I
  use ModIoUnit, ONLY: UnitTmp_, UnitTmp2_
  use ModNumConst, ONLY: cRadtoDeg, cDegToRad, cPi

  implicit none
  SAVE

  private   ! except
  public:: init_plot_box
  public:: set_plot_box
  public:: write_plot_box

  ! Size of current plot:
  integer :: nX, nY, nZ

  ! Resolution and range of current plot:
  real :: dX, dY, dZ
  real :: xMin, xMax, yMin, yMax, zMin, zMax

  ! Rotation angles around axis, center coordinates
  real :: xAngle, yAngle, zAngle
  real :: X0_D(3)
  ! Local results container:
  ! Array of values written to file:
  real, allocatable :: PlotVar_VIII(:,:,:,:)

  ! Coordinate conversion matrix
  real:: PlotToGm_DD(3,3)

contains
  !==========================================================================
  subroutine init_plot_box(iFile, nPlotVar)

    ! Set up the box grid for this plot file

    use ModMain,           ONLY: time_simulation, TypeCoordSystem
    use CON_axes,          ONLY: transform_matrix
    use ModCoordTransform, ONLY: show_rot_matrix

    integer, intent(in):: iFile, nPlotVar
    
    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'init_plot_box'
    !-----------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! Allocate results array and set up the box
    if(allocated(PlotVar_VIII)) RETURN

    ! Get plot area info from ModIO arrays:
    dX   = abs(plot_dx(1, iFile))
    dY   = abs(plot_dx(2, iFile))
    dZ   = abs(plot_dx(3, iFile))
    xMin = plot_range(1, iFile)
    xMax = plot_range(2, iFile)
    yMin = plot_range(3, iFile)
    yMax = plot_range(4, iFile)
    zMin = plot_range(5, iFile)
    zMax = plot_range(6, iFile)

    ! Get box orientation and center from ModIO arrays:
    X0_D = (/ (xMax+xMin)/2.0,(yMax+yMin)/2.0,(yMax+yMin)/2.0 /)
    xAngle = plot_normal(1,iFile) * cDegtoRad
    yAngle = plot_normal(1,iFile) * cDegtoRad
    zAngle = plot_normal(1,iFile) * cDegtoRad

    ! Set number of points:
    nX = nint((xMax - xMin)/dX) + 1
    nY = nint((yMax - yMin)/dY) + 1
    nZ = nint((zMax - zMin)/dZ) + 1

    ! Ensure dX, dY and dZ are compatible with the ranges
    dX = (xMax - xMin)/max(1, nX - 1)
    dY = (yMax - yMin)/max(1, nY - 1)
    dZ = (zMax - zMin)/max(1, nZ - 1)

    ! The 0 element is to count the number of blocks that
    ! contribute to a plot variable.
    allocate(PlotVar_VIII(0:nPlotVar,nX,nY,nZ))
    PlotVar_VIII = 0.0

    ! Get coordinate transformation matrix:
    PlotToGm_DD = transform_matrix(Time_Simulation, &
         TypeCoordPlot_I(iFile), TypeCoordSystem)

    if (DoTestMe) then
       write(*,*) NameSub//' iFile, nPlotVar=      ', iFile, nPlotVar
       write(*,*) NameSub//' Raw plot_dx=          ', plot_dx(:,iFile)
       write(*,*) NameSub//' Raw plot_range=       ', plot_range(:,iFile)
       write(*,*) NameSub//' dX, dY, dZ =      ', dX, dY, dZ
       write(*,*) NameSub//' x, y, z range = ',  &
            xMin, xMax, yMin,yMax,zMin,zMax
       write(*,*) NameSub,' nX, nY, nZ = ', nX, nY, nZ
       write(*,*) NameSub,' PlotToGm_DD:'
       call show_rot_matrix(PlotToGm_DD)
    end if

  end subroutine init_plot_box

  !============================================================================
  subroutine set_plot_box(iBlock, nPlotvar, Plotvar_GV)

    ! Interpolate the plot variables for block iBlock
    ! onto the spherical shell of the plot area.

    use ModGeometry,    ONLY: rMin_BLK
    use ModMain,        ONLY: BlkTest
    use ModInterpolate, ONLY: trilinear
    use BATL_lib,       ONLY: CoordMin_DB, nIjk_D, CellSize_DB, xyz_to_coord, MinI, MaxI,&
         MinJ, MaxJ, Mink, MaxK

    use ModCoordTransform, ONLY: rlonlat_to_xyz

    ! Arguments
    integer, intent(in) :: iBlock
    integer, intent(in) :: nPlotvar
    real,    intent(in) :: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)

    ! Local variables
    integer :: i, j, k, iVar

    real :: x, y, z, r, Lon, Lat
    real :: Xyz_D(3)
    real :: Coord_D(3), CoordNorm_D(3)

    character(len=*), parameter :: NameSub='set_plot_box'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    ! Check testing for block
    if(iBlock == BlkTest)then
       call CON_set_do_test(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if (DoTestMe) write(*,*) NameSub//' Called for iBlock=      ', iBlock

    ! Find box points and interpolate PlotVar
    
    
    
    


    
  end subroutine set_plot_box

  !============================================================================
  subroutine write_plot_box(iFile, nPlotVar, NameVar_V, &
       NameUnit, NameFile)

    ! Collect results from all blocks and write to single output file.
    use ModMpi
    use ModMain,     ONLY: time_simulation, n_step
    use ModProcMH,   ONLY: iProc, nProc, iComm
    use ModPlotFile, ONLY: save_plot_file

    integer,          intent(in) :: iFile, nPlotvar
    character(len=*), intent(in) :: NameFile, NameVar_V(nPlotVar), NameUnit

    integer :: iVar, iError, i, j, k
    character(len=500) :: NameVar

    logical :: DoTest,DoTestMe
    character(len=*), parameter :: NameSub='write_plot_box'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! This subroutine does not support HDF output.
    if(plot_form(iFile) == 'hgf') call CON_stop(NameSub// &
         ': HDF file type is not supported for BOX output.')

    ! Collect results to head node
    if(nProc > 1) call MPI_reduce_real_array(PlotVar_VIII, size(PlotVar_VIII),&
         MPI_SUM, 0, iComm, iError)

    ! Save results to disk
    if(iProc==0) then
       ! Build a single-line list of variable names.
       NameVar = 'x y z'
       do iVar=1, nPlotVar
          NameVar = trim(NameVar)  // ' ' // trim(NameVar_V(iVar))
       end do
       do k = 1, nZ
          do j = 1, nY; 
             do i = 1, nX; 
                if(PlotVar_VIII(0,i,j,k) <= 1.0) CYCLE
                PlotVar_VIII(1:nPlotVar,i,j,k) = &
                     PlotVar_VIII(1:nPlotVar,i,j,k) &
                     / PlotVar_VIII(0,i,j,k)
             end do; 
          end do; 
       end do

       ! Call save_plot_file to write data to disk.
       call save_plot_file(NameFile, &
            TypeFileIn=TypeFile_I(iFile), &
            StringHeaderIn=NameUnit, &
            nStepIn=n_step, &
            TimeIn=time_simulation, &
            NameVarIn = NameVar, &
            CoordMinIn_D = (/xMin, yMin, zMin/), &
            CoordMaxIn_D = (/xMax, yMax, zMax/), &
            VarIn_VIII = PlotVar_VIII(1:,:,:,:))
    end if

    ! Deallocate results arrays:.
    deallocate(PlotVar_VIII)

  end subroutine write_plot_box

  !============================================================================

end module ModPlotBox
