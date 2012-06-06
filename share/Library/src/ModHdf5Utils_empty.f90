module ModHdf5Utils
  
  use ModMpiOrig

  implicit none

  contains
  subroutine losHdf5PlotWrite(fileName, TypePosition, TypeStatus, StringHeader,&
    nStep, nBlocksUsed,nBlkUsedGlobal, offset,Time, nDimOut, nParam, nVar, n_D,  NameVar, PlotVarNameList,&
    unitList,coordMin, coordMax,iCoordMin, bBox, plotVarBlk, codeVersion, mpiComm, iProc, nProc)
    
    
    character (len=*), intent(in) :: fileName
    character (len=*), intent(in) :: TypePosition
    character (len=*), intent(in) :: PlotVarNameList(:)
    character (len=*), intent(in) :: unitList(:)
    character (len=10), intent(in) :: TypeStatus
    character (len=500), intent(in) :: StringHeader
    character (len=500), intent(in) :: NameVar
    character (len=501) :: HeaderString

    integer, intent(in) :: nDimOut, nParam, nVar, n_D(3), mpiComm, nStep, nBlocksUsed,nBlkUsedGlobal, offset
    integer, intent(in) :: iCoordMin(:,:), iProc, nProc
    real, intent(in) :: Time, coordMin(:), coordMax(:), codeVersion
    real, intent(in) ::plotVarBlk(:,:,:,:,:), bbox(:,:,:)


    write(*,*) "WARNING: HDF5 plotting is not enabled!"
    

  end subroutine losHdf5PlotWrite 
end module

