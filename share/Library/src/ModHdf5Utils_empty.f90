module ModHdf5Utils
  
  implicit none

  contains
   subroutine save_hdf5_file(FileName,TypePosition, TypeStatus, StringHeader,&
            nStep, NumberOfBlocksUsed, Time, nDimOut, nParam, nVar,&
            n_D, NameVar, NameUnits, MinimumBlockIjk, XYZMinMax, PlotVarBlk,&
            iComm, CoordMin, CoordMax)
    use ModUtilities, only: split_string

    integer, intent(in) :: nDimOut, nParam, nVar, n_D(3),iComm, nStep,NumberOfBlocksUsed
    character (len=*), intent(in) :: FileName
    character (len=*), intent(in) :: TypePosition, NameVar(nVar)
    character (len=10), intent(in) :: TypeStatus
    character (len=500), intent(in) :: StringHeader
    integer, intent(in) :: MinimumBlockIjk(:,:)
    real, intent(in) :: Time, PlotVarBlk(:,:,:,:,:), XYZMinMax(:,:,:)
    real, optional, intent(in) :: CoordMin(:), CoordMax(:)
    character (len=*), intent(in) ::  NameUnits
   end subroutine save_hdf5_file

  !=====================================================================
end module

