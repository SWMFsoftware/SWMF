module ModHdf5

  ! This empty module is used when HDF5 plotting is not enabled

  implicit none

  integer, parameter:: lNameVar = 10

contains
  !=========================================================================

  subroutine write_plot_hdf5(filename, plotType, plotVarNames, plotVarUnits,&
       nPlotVar, NotACut, nonCartesian, IsSphPlot, plot_dimensional, xmin, xmax, &
       ymin, ymax, zmin, zmax)

    use BATL_tree, only: iNode_B, iTree_IA, Coord0_,&
         Level_,iMortonNode_A
    use BATL_lib, only : CoordMin_DB, CoordMax_DB
    use ModNumConst
    use ModMpi
    use BATL_lib, ONLY : CoordMin_DB, MaxDim
    use BATL_mpi, only : barrier_mpi

    use ModProcMH, ONLY: iProc
    
    integer,                 intent(in):: nPlotVar
    character(len=80),       intent(in):: filename
    character(len=lnamevar), intent(in):: plotvarnames(nplotvar)
    character(len=lNameVar),  intent(in):: plotVarUnits(nPlotVar)
    character (len=*), intent(in) :: plotType
    real, intent(in)  ::  xMin, xMax, yMin, yMax, zMin, zMax
    logical, intent(in) :: NotACut, plot_dimensional, nonCartesian, IsSphPlot
    !----------------------------------------------------------------------
    
    if(iProc==0) write (*,*) "ERROR: HDF5 plotting is not enabled!"

  end subroutine write_plot_hdf5
  !=========================================================================
  subroutine write_var_hdf5(iFile, plotType, iBlock, H5Index,nPlotVar,PlotVar, &
       xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock,&
       isNonCartesian, NotACut, nCell,H5Advance)

    ! Save all cells within plotting range, for each processor

    use ModSize, ONLY: MinI, MaxI, MinJ, MaxJ, MinK, MaxK



    implicit none

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

  end subroutine write_var_hdf5

  !=========================================================================
  subroutine close_sph_hdf5_plot(nPlotVar)
    integer, intent(in) :: nPlotVar
  end subroutine close_sph_hdf5_plot
  !=========================================================================

  subroutine write_sph_var_hdf5(iCell, jCell, H5Index, rplot, theta_out,phi_out,&
   dtheta_plot, dphi_plot, nTheta, nPhi, nPlotVar, PointVar, plotVarNames)

    integer, intent(in) :: iCell, jCell, nPlotVar, nTheta, nPhi, H5Index
    real, intent(in) :: theta_out, phi_out, PointVar(nPlotVar)
    real, intent(in) :: rplot, dtheta_plot, dphi_plot
    character(len=*), intent(in):: plotVarNames(nPlotVar)
    end subroutine write_sph_var_hdf5
  !=========================================================================
 !==================================================================================
  subroutine init_sph_hdf5_plot(nPlotVar, filename, plotVarNames, PlotVarUnits, nTheta,&
  nPhi, rplot) 
  
    integer, intent(in) :: nPlotVar, nTheta, nPhi
    character(len=80),       intent(in):: filename
    character(len=*), intent(in):: plotVarNames(nPlotVar)
    character(len=*),  intent(in):: plotVarUnits(nPlotVar)
    real, intent(in) :: rplot
  end subroutine init_sph_hdf5_plot
end module ModHdf5
