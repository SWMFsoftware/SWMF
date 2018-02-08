!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPlotBox

  use BATL_lib, ONLY: &
       test_start, test_stop, iBlockTest

  use ModIO, ONLY: plot_dx, plot_range, plot_normal, TypeCoordPlot_I, &
       plot_form, TypeFile_I, plot_type, ObsPos_DI, IsObsBox_I

  use ModNumConst,       ONLY  : cDegToRad, cTwoPi

  use ModCoordTransform, ONLY  : xyz_to_lonlat

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
  real :: xMin, xMax, yMin, yMax, zMin, zMax, xLen, yLen, zLen

  ! Rotation angles around axis, center coordinates
  real :: xAngle, yAngle, zAngle
  real :: Xyz0_D(3), Xyz0Hgi_D(3)

  ! Local results container:
  ! Array of values written to file:
  real, allocatable :: PlotVar_VIII(:,:,:,:)

  ! Coordinate conversion matrix
  real:: PlotToGm_DD(3,3)

contains
  !============================================================================
  subroutine init_plot_box(iFile, nPlotVar)

    ! Set up the box grid for this plot file

    use ModMain,           ONLY: time_simulation, TypeCoordSystem
    use CON_axes,          ONLY: transform_matrix
    use ModCoordTransform, ONLY: show_rot_matrix, cross_product

    integer, intent(in)        :: iFile, nPlotVar

    real :: Los_D(3), aUnit_D(3), bUnit_D(3), ObsPos_D(3)

    logical                    :: DoTest
    character(len=*), parameter:: NameSub = 'init_plot_box'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Allocate results array and set up the box
    if(allocated(PlotVar_VIII)) RETURN

    ! Get box resolution, center and size from ModIO arrays:
    dX   = abs(plot_dx(1, iFile))
    dY   = abs(plot_dx(2, iFile))
    dZ   = abs(plot_dx(3, iFile))
    Xyz0_D = plot_range(1:3, iFile)
    xLen = plot_range(4, iFile)
    yLen = plot_range(5, iFile)
    zLen = plot_range(6, iFile)

    xMin = plot_range(1, iFile) - xLen/2
    xMax = plot_range(1, iFile) + xLen/2
    yMin = plot_range(2, iFile) - yLen/2
    yMax = plot_range(2, iFile) + yLen/2
    zMin = plot_range(3, iFile) - zLen/2
    zMax = plot_range(3, iFile) + zLen/2

    ! Get coordinate transformation matrix:
    PlotToGm_DD = transform_matrix(Time_Simulation, &
         TypeCoordPlot_I(iFile), TypeCoordSystem)

    if(IsObsBox_I(iFile))then
       ! This is the tilt (roll)
       xAngle = plot_normal(1,iFile) * cDegtoRad

       ! Translate image center from LOS coordinates to TypeCoordPlot_I(iFile)
       ObsPos_D = ObsPos_DI(:,iFile)
       Los_D = ObsPos_D/sqrt(sum(ObsPos_D**2))
       aUnit_D = cross_product((/0.,0.,1./), Los_D)
       aUnit_D = aUnit_D/sqrt(sum(aUnit_D**2))
       bUnit_D = cross_product(Los_D, aUnit_D)
       bUnit_D = bUnit_D/sqrt(sum(bUnit_D**2))
       Xyz0Hgi_D =  Xyz0_D(1)*Los_D + Xyz0_D(2)*aUnit_D + Xyz0_D(3)*bUnit_D

       ! Observer position is with respect to center of box.
       ! Convert observer location to longitude and latitude.
       ! Coordinate system is TypeCoordPlot_I(iFile)
       call xyz_to_lonlat(ObsPos_DI(:,iFile) - Xyz0Hgi_D, zAngle, yAngle)
       zAngle = cTwoPi - zAngle
    else
       ! Get box orientation from ModIO arrays:
       xAngle = plot_normal(1,iFile) * cDegtoRad
       yAngle = plot_normal(2,iFile) * cDegtoRad
       zAngle = plot_normal(3,iFile) * cDegtoRad
    end if

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

    if (DoTest) then
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

    call test_stop(NameSub, DoTest)
  end subroutine init_plot_box
  !============================================================================

  subroutine set_plot_box(iBlock, nPlotvar, Plotvar_GV)

    ! Interpolate the plot variables for block iBlock
    ! onto the spherical shell of the plot area.

    use ModInterpolate, ONLY: trilinear
    use BATL_lib,       ONLY: CoordMin_DB, nIjk_D, CellSize_DB, xyz_to_coord, &
         MinI, MaxI, MinJ, MaxJ, Mink, MaxK
    use ModPhysics,     ONLY: rBody

    use ModCoordTransform, ONLY: rot_matrix_x, rot_matrix_y, rot_matrix_z

    ! Arguments
    integer, intent(in) :: iBlock
    integer, intent(in) :: nPlotvar
    real,    intent(in) :: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)

    ! Local variables
    integer :: i, j, k, iVar

    real :: Xyz_D(3), XyzRot_D(3), XyzGm_D(3)
    real :: Coord_D(3), CoordNorm_D(3)

    ! Check testing for block
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_plot_box'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)

    ! Shift box elements with origin of box             
    do k = 1, nZ
       Xyz_D(3) = zMin + (k-1)*dZ - Xyz0_D(3)

       do j = 1, nY
          Xyz_D(2) = yMin + (j-1)*dY - Xyz0_D(2)

          do i = 1, nX
             Xyz_D(1) = xMin + (i-1)*dX - Xyz0_D(1)

             ! Rotate box
             XyzRot_D = matmul(rot_matrix_z(-zAngle), &
                  matmul(rot_matrix_y(-yAngle), &
                  matmul(rot_matrix_x(-xAngle), Xyz_D)))

             ! Shift box back and Get Gm coordinates (i.e., TypeCoordSystem)
             XyzGm_D = matmul(PlotToGm_DD, XyzRot_D + Xyz0Hgi_D)

             ! When inside Body keep default plot values
             if(sqrt(sum(XyzGm_D**2)) < rBody)CYCLE

             ! Get generalized coordinates
             call xyz_to_coord(XyzGm_D, Coord_D)

             ! Normalize to block coordinates
             CoordNorm_D = &
                  (Coord_D - CoordMin_DB(:,iBlock))/CellSize_DB(:,iBlock) + 0.5

             ! Check if point is inside block
             if(any(CoordNorm_D < 0.5)) CYCLE
             if(any(CoordNorm_D > nIjk_D + 0.5)) CYCLE

             ! Compute the interpolated values at the current location
             PlotVar_VIII(0,i,j,k) = 1.0
             do iVar=1, nPlotVar
                ! Interpolate up to ghost cells.
                PlotVar_VIII(iVar,i,j,k) = trilinear(PlotVar_GV(:,:,:,iVar),&
                     MinI, MaxI, MinJ, MaxJ, MinK, MaxK, CoordNorm_D)
             end do
          end do
       end do
    end do

    call test_stop(NameSub, DoTest, iBlock)
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

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_plot_box'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! This subroutine does not support HDF output.
    if(plot_form(iFile) == 'hdf') call stop_mpi(NameSub// &
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
       NameVar = trim(NameVar) // ' xAngle yAngle zAngle'

       ! Correct for double counting in MPI_reduce_real_array
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
            ParamIn_I = (/ xAngle, yAngle, zAngle /), &
            NameVarIn = NameVar, &
            CoordMinIn_D = (/ xMin, yMin, zMin /), &
            CoordMaxIn_D = (/ xMax, yMax, zMax /), &
            VarIn_VIII = PlotVar_VIII(1:,:,:,:))
    end if

    ! Deallocate results arrays:.
    deallocate(PlotVar_VIII)

    call test_stop(NameSub, DoTest)
  end subroutine write_plot_box
  !============================================================================

end module ModPlotBox
!==============================================================================
