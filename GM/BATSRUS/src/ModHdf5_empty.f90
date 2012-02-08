module ModHdf5

  ! This empty module is used when HDF5 plotting is not enabled

  implicit none

  integer, parameter:: lNameVar = 10

contains
  !============================================================================
  subroutine hdf5_setup

    ! do nothing

  end subroutine hdf5_setup
  !=========================================================================
  subroutine write_plot_hdf5(filename, PlotVarNames, PlotVarUnits, nPlotVar, &
       xmin, xmax, ymin, ymax, zmin, zmax, nBLKcells)

    use ModIO, ONLY : nplotvarmax

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar), intent(in):: plotVarUnits(nPlotVar)
    real,                    intent(in):: xMin, xMax, yMin, yMax, zMin, zMax
    integer,                 intent(in):: nBLKcells
    !----------------------------------------------------------------------
    write (*,*) "ERROR: HDF5 plotting is not enabled!"

  end subroutine write_plot_hdf5
  !=========================================================================
  subroutine write_var_hdf5(PlotVar, nPlotVar, H5Index, iBLK)

    use ModMain, ONLY : nI, nJ, nK, nBlockMax
    use ModIO, ONLY: nPlotVarMax

    integer, intent(in) :: nPlotVar, H5Index, iBLK
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVarMax)
    !----------------------------------------------------------------------
    write (*,*) "ERROR: HDF5 plotting not enabled!"

  end subroutine write_var_hdf5

end module ModHdf5
