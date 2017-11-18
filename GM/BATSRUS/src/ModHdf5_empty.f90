!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModHdf5

  ! This empty module is used when HDF5 plotting is not enabled

  implicit none

  integer, parameter:: lNameVar = 10

contains
  !============================================================================
  subroutine init_hdf5_plot(iFile, plotType, nPlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian, NotACut)

    ! Arguments

    integer, intent(in) :: iFile
    integer, intent(in) :: nPlotVar
    real,    intent(in) :: xMin, xMax, yMin, yMax, zMin, zMax
    logical, intent(in) :: isNonCartesian, NotACut
    character (len=*), intent(in) :: plotType
    real,    intent(inout):: DxBlock, DyBlock, DzBlock

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'init_hdf5_plot'
    !--------------------------------------------------------------------------
  end subroutine init_hdf5_plot
  !============================================================================

  subroutine write_plot_hdf5(filename, plotType, plotVarNames, plotVarUnits,&
       nPlotVar, NotACut, nonCartesian, IsSphPlot, plot_dimensional, &
       xmin, xmax, ymin, ymax, zmin, zmax)

    use ModProcMH, ONLY: iProc

    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lnamevar), intent(in):: plotvarnames(nplotvar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    character (len=*), intent(in) :: plotType
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
    logical, intent(in) :: NotACut, plot_dimensional, nonCartesian, IsSphPlot
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_plot_hdf5'
    !--------------------------------------------------------------------------
    if(iProc==0) write (*,*) "ERROR: HDF5 plotting is not enabled!"

  end subroutine write_plot_hdf5
  !============================================================================
  subroutine write_var_hdf5(iFile, plotType, iBlock, H5Index,nPlotVar,PlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian, NotACut, nCell,H5Advance)

    ! Save all cells within plotting range, for each processor

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK

    ! Arguments

    integer, intent(in)   :: iFile, iBlock, H5Index
    integer, intent(in)   :: nPlotVar
    real,    intent(in)   :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
    real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
    logical, intent(in) :: isNonCartesian, NotACut
    character (len=*), intent(in) :: plotType
    real,    intent(inout):: DxBlock,DyBlock,DzBlock
    integer, intent(out)  :: nCell
    logical, intent(out) :: H5advance

    logical:: DoTest
    character(len=*), parameter:: NameSub = 'write_var_hdf5'
    !--------------------------------------------------------------------------
  end subroutine write_var_hdf5
  !============================================================================

end module ModHdf5
!==============================================================================
