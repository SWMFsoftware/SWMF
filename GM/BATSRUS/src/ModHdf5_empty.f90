!============================================================================
! This empty module is used when HDF5 support is not enabled
!============================================================================

module ModHdf5

contains

  subroutine hdf5_setup()

  end subroutine hdf5_setup

  subroutine write_plot_hdf5(filename, plotVarNames, nPlotVar, ifile)

    use ModIO, ONLY : nplotvarmax

    character (len=80), intent(in)  :: filename
    character (len=10), intent(in)  :: plotVarNames(nplotvarmax)
    integer, intent(in)  :: nPlotVar, ifile

    write (*,*) "ERROR: HDF5 support not enabled!"

  end subroutine write_plot_hdf5

  subroutine write_var_hdf5(PlotVar, nPlotVar, iBLK)

    use ModMain, ONLY : nI, nJ, nK, nBlockMax
    use ModIO, ONLY: nPlotVarMax

    integer, intent(in) :: nPlotVar, iBLK
    real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVarMax)

    write (*,*) "ERROR: HDF5 support not enabled!"

  end subroutine write_var_hdf5

end module ModHdf5
