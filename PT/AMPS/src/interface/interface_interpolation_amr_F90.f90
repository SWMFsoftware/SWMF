!Module for interfacing AMR interpolation procedure written in FORTRAN
module ModInterfaceInterpolationAMR

  !inclduing AMR interpolation module itself
  USE ModInterpolateAMR, ONLY: interpolate_amr

  implicit none

contains

  subroutine find(&
       nDim, Xyz_D, iProc, iBlock, &
       XyzCorner_D, Dxyz_D, IsOut)
    !\
    ! Block AMR Grid characteristic:
    ! the search routine, which returns, for a given point,
    ! the block to which this points belong and the processor, at which
    ! this block is allocated, as well as the block parameters:
    ! coordinates of the left corner (the point in the block with the
    ! minvalue of coordinates) and (/Dx, Dy, Dz/)
    !/
    integer, intent(in) :: nDim
    !\
    ! "In"- the coordinates of the point, "out" the coordinates of the
    ! point with respect to the block corner. In the most cases
    ! XyzOut_D = XyzIn_D - XyzCorner_D, the important distinction,
    ! however, is the periodic boundary, near which the jump in the
    ! stencil coordinates might occur. To handle the latter problem,
    ! we added the "out" intent. The coordinates for the stencil
    ! and input point are calculated and recalculated below with
    ! respect to the block corner.
    !/
    real,  intent(inout):: Xyz_D(nDim)
    integer, intent(out):: iProc, iBlock !processor and block number
    !\
    ! Block left corner coordinates and the grid size:
    !/
    real,    intent(out):: XyzCorner_D(nDim), Dxyz_D(nDim)
    logical, intent(out):: IsOut !Point is out of the domain. 
    !-----------------------------------------------------------------
    ! call interface find function written in C
    call INTERFACE__INTERPOLATION_AMR__find(&
         nDim, Xyz_D, iProc, iBlock, XyzCorner_D, Dxyz_D, IsOut)
  end subroutine find

  !===============================================================

  subroutine INTERFACE__INTERPOLATION_AMR__interpolate_amr(&
       nDim, XyzIn_D, nIndexes, nCell_D, nGridOut, Weight_I,&
       iIndexes_II, IsSecondOrder, UseGhostCell)
    !\
    ! INPUT PARAMETERS
    !/
    !\
    ! Number of dimensions
    !/
    integer, intent(in) :: nDim
    !\
    ! Number of indexes. Usually, nIndexes = nDim + 1, for three cell indexes
    ! and a block number. If the latter information is not needed,
    ! nIndexes=nDim should be added
    !/
    integer, intent(in) :: nIndexes
    !\
    ! Point coordinates
    !/
    real,    intent(in) :: XyzIn_D(nDim)
    !\
    ! Block AMR Grid characteristic: number of cells
    !/
    integer, intent(in) :: nCell_D(nDim)
    !\
    ! Yet another Block AMR Grid characteristic:
    ! the search routine, which returns, for a given point,
    ! the block to which this points belong and the processor, at which
    ! this block is allocated, as well as the block parameters:
    ! coordinates of the left corner (the point in the block with the
    ! minvalue of coordinates) and (/Dx, Dy, Dz/)
    !/
    !\
    ! Do or do not use ghost cells
    !/
    logical, intent(in):: UseGhostCell
    !\
    !OUTPUT PARAMETERS
    !/
    !\
    !Number of grid points involved into interpolation
    !/
    integer, intent(out):: nGridOut
    !\
    ! Interpolation weights (only the first nGridOut values are meaningful
    !/
    real,    intent(out):: Weight_I(2**nDim)
    !\
    ! Cell(+block) indexes and processor number for grid points to be
    ! invilved into interpolation. iProc numbers are stored in
    ! iIndexes_II(0,:)
    !/
    integer, intent(out):: iIndexes_II(0:nIndexes,2**nDim)
    !\
    !The following is true if stencil does not employ
    !the out-of-grid points
    !/
    logical, intent(out):: IsSecondOrder  
    !--------------------------------------------------------------------
    ! call the interpolation prcodure itself
    call interpolate_amr(&
         nDim=nDim, &
         XyzIn_D=XyzIn_D, &
         nIndexes=nDim+1,&
         find=find, &
         nCell_D=nCell_D(1:nDim),&
         nGridOut=nGridOut,&
         Weight_I=Weight_I,&
         iIndexes_II=iIndexes_II,&
         IsSecondOrder=IsSecondOrder,&
         UseGhostCell=UseGhostCell)
  end subroutine INTERFACE__INTERPOLATION_AMR__interpolate_amr
  
end module ModInterfaceInterpolationAMR
