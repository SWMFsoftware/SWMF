!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!===========================TESTS============================================
module ModTestInterpolateAMR
  use ModInterpolateCellAMR, ONLY: interpolate=>interpolate_amr_cell, &
       order_connectivity_list
  implicit none
  !\
  ! Shift of the iGrid point in the stencil with respect to the
  ! first one
  !/

contains
  !==================================================================
  subroutine interpolate_cell_amr(nDim, Xyz_D, XyzCell_D, DXyz_D,   &
       nNeighbor, XyzNeighbor_DI, DiLevelNei_I, IsOut_I,&
       WeightOut_I, IsSecondOrder)
    integer, intent(in) :: nDim, nNeighbor
    real,    intent(in) :: Xyz_D(nDim), XyzCell_D(nDim), DXyz_D(nDim)
    real,    intent(in) :: XyzNeighbor_DI(nDim, nNeighbor)
    integer, intent(in) :: DiLevelNei_I(nNeighbor)
    logical, intent(in) :: IsOut_I(nNeighbor)
    real,    intent(out):: WeightOut_I(0:nNeighbor)
    logical, intent(out):: IsSecondOrder
    integer,parameter   :: iCellId_I(1) = 0
    integer             :: iCellId_II(1,nNeighbor)
    real                :: DXyzInv_D(nDim)
    integer             :: iCellId_IIII(1,2**nDim,2**nDim,2**nDim)
    real                :: XyzGrid_DIII(nDim,0:2**nDim,2**nDim,2**nDim)
    real                :: Weight_I(2**nDim)
    integer             :: iLevel_II(2**nDim,2**nDim), iIndexes_II(1,2**nDim)
    logical:: IsOut_II( 2**nDim,2**nDim) 
    integer:: nGridOut, iNeighbor
    !-------------------
    do iNeighbor = 1, nNeighbor; iCellId_II(1,iNeighbor) = iNeighbor; end do

    call order_connectivity_list(nDim, XyzCell_D, DXyz_D, 1, iCellId_I,&
       nNeighbor, XyzNeighbor_DI, DiLevelNei_I, iCellId_II, IsOut_I,   &
       DXyzInv_D, iCellId_IIII, XyzGrid_DIII, iLevel_II, IsOut_II)

    call interpolate(nDim, Xyz_D, XyzCell_D, DXyzInv_D,       &
       1, iCellId_IIII, XyzGrid_DIII, iLevel_II, IsOut_II,             &
       nGridOut, Weight_I, iIndexes_II,                                &
       IsSecondOrder)
    WeightOut_I = 0 
    WeightOut_I(iIndexes_II(1,1:nGridOut)) = Weight_I(1:nGridOut) 
  end subroutine interpolate_cell_amr
end module ModTestInterpolateAMR

program test_interpolate_cell_amr

  use ModTestInterpolateAMR, ONLY: interpolate_cell_amr
  implicit none
  !Test grid
  !        ___________
  !       |   |   |   |
  !       |___|___|___|
  !       |   |   |_|  
  !    ___|___|___|_|_
  !   |       |   |   |
  !   |       |___|___|
  !   |       |
  !   |_______|
  !
  real,parameter:: XyzNeighbor_DI(2,9) = reshape( (/&
       -1.00, -1.00,        &!
        0.50, -0.50,        &!
        1.50, -0.50,        &!
       -0.50,  0.50,        &!
        1.25,  0.25,        &!
        1.25,  0.75,        &!
       -0.50,  1.50,        &!
        0.50,  1.50,        &!
        1.50,  1.50  /),(/2,9/) )
  real, parameter:: XyzCell_D(2) = 0.50, DXyz_D(2) = 1.0
  integer,parameter:: DiLevelNei_I(9) = (/-1, 0, 0, 0, 1, 1, 0, 0, 0/) 
  real:: Xyz_D(2), XyzInterpolated_D(2), WeightOut_I(0:9)
  logical :: IsSecondOrder
  logical,parameter:: IsOut_I(9) = .false.
  integer:: iPoint
  !-----------------
  Xyz_D =(/0.70, 0.30 /)
  call interpolate_cell_amr(2, Xyz_D, XyzCell_D, DXyz_D,   &
       9, XyzNeighbor_DI, DiLevelNei_I, IsOut_I,&
       WeightOut_I, IsSecondOrder)
  XyzInterpolated_D = XyzCell_D*WeightOut_I(0)
  do iPoint = 1, 9
     XyzInterpolated_D = XyzInterpolated_D + &
          XyzNeighbor_DI(:, iPoint)*WeightOut_I(iPoint)
  end do
  write(*,*)'Xyz_D=',Xyz_D
  write(*,*)'WeightOut_I=',WeightOut_I
  write(*,*)'XyzInterpolated_D=',XyzInterpolated_D
end program test_interpolate_cell_amr

subroutine CON_stop(StringError)

  implicit none
  character (len=*), intent(in) :: StringError
  !----------------------------------------------------------------------------

  write(*,'(a)')StringError
  write(*,'(a)')'!!! SWMF_ABORT !!!'
  stop

end subroutine CON_stop

