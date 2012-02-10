module ModHdf5

  ! This empty module is used when HDF5 plotting is not enabled

  implicit none

  integer, parameter:: lNameVar = 10

contains
  !=========================================================================
  subroutine write_plot_hdf5(NameFile, NameVar_V, NameUnit_V, &
       nPlotVar, xMin, xMax, yMin, yMax, zMin, zMax, nCell)

    character(len=80),       intent(in):: NameFile
    integer,                 intent(in):: nPlotVar
    character(len=lNameVar), intent(in):: NameVar_V(nPlotVar)
    character(len=lNameVar), intent(in):: NameUnit_V(nPlotVar)
    real,                    intent(in):: xMin, xMax, yMin, yMax, zMin, zMax
    integer,                 intent(in):: nCell
    !----------------------------------------------------------------------
    write (*,*) "ERROR: HDF5 plotting is not enabled!"

  end subroutine write_plot_hdf5
  !=========================================================================
  subroutine write_var_hdf5(PlotVar_GI, nPlotVar, iBlockH5, iBlock)

    use ModMain, ONLY : nI, nJ, nK

    integer, intent(in) :: nPlotVar, iBlockH5, iBlock
    real, intent(in) :: PlotVar_GI(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
    !----------------------------------------------------------------------
    write (*,*) "ERROR: HDF5 plotting not enabled!"

  end subroutine write_var_hdf5

end module ModHdf5
