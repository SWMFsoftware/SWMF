module ModHdf5

  ! This empty module is used when HDF5 plotting is not enabled

  implicit none

  integer, parameter:: lNameVar = 10

contains
  !=========================================================================
  subroutine write_plot_hdf5(filename, plot_type,plotVarNames, plotVarUnits,&
       nPlotVar,isCutFile, nonCartesian,plot_dimensional, xmin, xmax, &
       ymin, ymax, zmin, zmax)

    use ModProcMH, ONLY: iProc

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=3),        intent(in) :: plot_type
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    logical, intent(in) :: isCutFile, nonCartesian,plot_dimensional
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
 
    !----------------------------------------------------------------------
    if(iProc==0) write (*,*) "ERROR: HDF5 plotting is not enabled!"

  end subroutine write_plot_hdf5
  !=========================================================================
  subroutine write_var_hdf5(PlotVar_GI, nPlotVar, iBlockH5, iBlock, &
       IsNonCartesian)

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    integer, intent(in) :: nPlotVar, iBlockH5, iBlock
    real, intent(in) :: PlotVar_GI(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
    logical, intent(in) :: IsNonCartesian
    !----------------------------------------------------------------------

  end subroutine write_var_hdf5

  subroutine write_cut_var_hdf5(iFile, plotType, iBlock, H5Index,nPlotVar,PlotVar, &
     xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
     isNonCartesian,nCell,H5Advance)

      use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK
      character (len=3), intent(in) :: plotType
      integer, intent(in)   :: iFile, iBlock, H5Index
      integer, intent(in)   :: nPlotVar
      real,    intent(in)   :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
      real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
      logical, intent(in) :: isNonCartesian
      real,    intent(inout):: DxBlock,DyBlock,DzBlock
      integer, intent(out)  :: nCell 
      logical, intent(out) :: H5advance

  end subroutine write_cut_var_hdf5


end module ModHdf5
