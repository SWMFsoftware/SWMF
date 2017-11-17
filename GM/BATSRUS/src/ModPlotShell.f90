!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModPlotShell

  use BATL_lib, ONLY: &
       test_start, test_stop, iBlockTest

  use ModIO
  use ModNumConst, ONLY: cRadtoDeg, cDegToRad

  implicit none

  SAVE

  private   ! except
  public:: init_plot_shell
  public:: set_plot_shell
  public:: write_plot_shell

  ! Size of current plot:
  integer :: nR, nLon, nLat

  ! Ranges for current plot:
  real :: dR, dLat, dLon
  real :: rMin, rMax, LonMin, LonMax, LatMin, LatMax

  ! Local results container:
  ! Array of values written to file:
  real, allocatable :: PlotVar_VIII(:,:,:,:)

  ! Coordinate conversion matrix
  real:: PlotToGm_DD(3,3)

contains
  !============================================================================
  subroutine init_plot_shell(iFile, nPlotVar)

    ! set up the shell grid for this plot file

    use ModMain,           ONLY: time_simulation, TypeCoordSystem
    use CON_axes,          ONLY: transform_matrix
    use ModCoordTransform, ONLY: show_rot_matrix

    integer, intent(in):: iFile, nPlotVar

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_plot_shell'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! Allocate results array and set up all the spherical shell
    if(allocated(PlotVar_VIII)) RETURN

    ! Get plot area info from ModIO arrays:
    dR     = abs(plot_dx(1, iFile))
    dLon   = plot_dx(2, iFile) * cDegtoRad
    dLat   = plot_dx(3, iFile) * cDegtoRad
    rMin   = plot_range(1, iFile)
    rMax   = plot_range(2, iFile)
    LonMin = plot_range(3, iFile) * cDegtoRad
    LonMax = plot_range(4, iFile) * cDegtoRad
    LatMin = plot_range(5, iFile) * cDegtoRad
    LatMax = plot_range(6, iFile) * cDegtoRad

    ! Set number of points:
    nR   = nint((rMax - rMin)/dR)       + 1
    nLon = nint((LonMax - LonMin)/dLon) + 1
    nLat = nint((LatMax - LatMin)/dLat) + 1

    ! Ensure dR, dLon and dLat are compatible with the ranges
    dR   = (rMax - rMin)/max(1, nR - 1)
    dLon = (LonMax - LonMin)/max(1, nLon - 1)
    dLat = (LatMax - LatMin)/max(1, nLat - 1)

    ! The 0 element is to count the number of blocks that
    ! contribute to a plot variable.
    allocate(PlotVar_VIII(0:nPlotVar,nR,nLon,nLat))
    PlotVar_VIII = 0.0

    ! Get coordinate transformation matrix:
    PlotToGm_DD = transform_matrix(Time_Simulation, &
         TypeCoordPlot_I(iFile), TypeCoordSystem)

    if (DoTest) then
       write(*,*) NameSub//' iFile, nPlotVar=      ', iFile, nPlotVar
       write(*,*) NameSub//' Raw plot_dx=          ', plot_dx(:,iFile)
       write(*,*) NameSub//' Raw plot_range=       ', plot_range(:,iFile)
       write(*,*) NameSub//' dR, dLon, dLat =      ', dR, dLon, dLat
       write(*,*) NameSub//' r, Lon, Lat range = ',  &
            rMin, rMax, LonMin,LonMax,LatMin,LatMax
       write(*,*) NameSub,' nR, nLon, nLat, dR, dLon, dLat = ', &
            nR, nLon, nLat, dR, dLon, dLat
       write(*,*) NameSub,' PlotToGm_DD:'
       call show_rot_matrix(PlotToGm_DD)
    end if

    call test_stop(NameSub, DoTest)
  end subroutine init_plot_shell
  !============================================================================

  subroutine set_plot_shell(iBlock, nPlotvar, Plotvar_GV)

    ! Interpolate the plot variables for block iBlock
    ! onto the spherical shell of the plot area.

    use ModGeometry,    ONLY: rMin_BLK
    use ModInterpolate, ONLY: trilinear
    use BATL_lib,       ONLY: CoordMin_DB, nIjk_D, CellSize_DB, xyz_to_coord
    use ModCoordTransform, ONLY: rlonlat_to_xyz

    ! Arguments
    integer, intent(in) :: iBlock
    integer, intent(in) :: nPlotvar
    real,    intent(in) :: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)

    ! Local variables
    integer :: i, j, k, iVar

    real :: r, Lon, Lat
    real :: XyzPlot_D(3), XyzGm_D(3)
    real :: Coord_D(3), CoordNorm_D(3)

    ! Check testing for block
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'set_plot_shell'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    if(iBlock == iBlockTest)then
       call CON_set_do_test(NameSub, DoTest, DoTest)
    else
       DoTest = .false.; DoTest = .false.
    end if

    if (DoTest) write(*,*) NameSub//' Called for iBlock=      ', iBlock

    ! Loop through shell points and interpolate PlotVar
    do i = 1, nR
       r = rMin + (i-1)*dR

       if(r < rMin_BLK(iBlock)) CYCLE
       ! if(r > rMax_B(iBlock)) CYCLE

       do j = 1, nLon
          Lon = LonMin + (j-1)*dLon
          do k = 1, nLat
             Lat = LatMin + (k-1)*dLat

             ! Convert to Cartesian coordinates
             call rlonlat_to_xyz(r,Lon,Lat,XyzPlot_D)

             ! Convert from plot coordinates to BATSRUS grid:
             XyzGm_D = matmul(PlotToGm_DD, XyzPlot_D)

             ! Convert to generalized coordinates
             call xyz_to_coord(XyzGm_D, Coord_D)

             ! Normalize to block coordinates
             CoordNorm_D = &
                  (Coord_D - CoordMin_DB(:,iBlock))/CellSize_DB(:,iBlock) + 0.5

             ! Check if point is inside block
             if(any(CoordNorm_D < 0.5)) CYCLE
             if(any(CoordNorm_D > nIjk_D + 0.5)) CYCLE

             ! compute the interpolated values at the current location
             PlotVar_VIII(0,i,j,k) = 1.0
             do iVar=1, nPlotVar
                ! Interpolate up to ghost cells.
                PlotVar_VIII(iVar,i,j,k) = trilinear(PlotVar_GV(:,:,:,iVar),&
                     MinI, MaxI, MinJ, MaxJ, MinK, MaxK, CoordNorm_D)
             end do

          end do ! lon loop
       end do    ! lat loop
    end do       ! r loop

    call test_stop(NameSub, DoTest, iBlock)
  end subroutine set_plot_shell
  !============================================================================

  subroutine write_plot_shell(iFile, nPlotVar, NameVar_V, &
       NameUnit, NameFile)

    ! Collect results from all blocks and write to single output file.
    use ModMpi
    use ModMain,     ONLY: time_simulation, n_step
    use ModProcMH,   ONLY: iProc, nProc, iComm
    use ModPlotFile, ONLY: save_plot_file

    integer,          intent(in) :: iFile, nPlotvar
    character(len=*), intent(in) :: NameFile, NameVar_V(nPlotVar), NameUnit

    integer :: iVar, iR, iLon, iLat, iError
    character(len=500) :: NameVar

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_plot_shell'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest)

    ! This subroutine does not support HDF output.
    if(plot_form(iFile) == 'hdf') call CON_stop(NameSub// &
         ': HDF file type not supported for Geo Sphere output.')

    ! Collect results to head node
    if(nProc > 1) call MPI_reduce_real_array(PlotVar_VIII, size(PlotVar_VIII),&
         MPI_SUM, 0, iComm, iError)

    ! Save results to disk
    if(iProc==0) then
       ! Build a single-line list of variable names.
       NameVar = 'r lon lat'
       do iVar=1, nPlotVar
          NameVar = trim(NameVar)  // ' ' // trim(NameVar_V(iVar))
       end do

       do iR = 1, nR; do iLon = 1, nLon; do iLat = 1, nLat
          if(PlotVar_VIII(0,iR,iLon,iLat) <= 1.0) CYCLE
          PlotVar_VIII(1:nPlotVar,iR,iLon,iLat) = &
               PlotVar_VIII(1:nPlotVar,iR,iLon,iLat) &
               / PlotVar_VIII(0,iR,iLon,iLat)
       end do; end do; end do

       ! Call save_plot_file to write data to disk.
       call save_plot_file(NameFile, &
            TypeFileIn=TypeFile_I(iFile), &
            StringHeaderIn=NameUnit, &
            nStepIn=n_step, &
            TimeIn=time_simulation, &
            NameVarIn = NameVar, &
            CoordMinIn_D = (/rMin, cRadtoDeg*LonMin, cRadtoDeg*LatMin/), &
            CoordMaxIn_D = (/rMax, cRadtoDeg*LonMax, cRadtoDeg*LatMax/), &
            VarIn_VIII = PlotVar_VIII(1:,:,:,:))
    end if

    ! Deallocate results arrays:.
    deallocate(PlotVar_VIII)

    call test_stop(NameSub, DoTest)
  end subroutine write_plot_shell
  !============================================================================

end module ModPlotShell
!==============================================================================
