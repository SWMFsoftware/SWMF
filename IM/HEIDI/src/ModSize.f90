module ModSize

  !\
  ! Block parameters.
  !/
  integer, parameter :: nBLK=2

  !\
  ! Cell array parameters.
  !   Minimum is 2x2x2 for fixed grid, 4x4x4 with AMR.
  !   Use even numbers only.
  !/
  integer, parameter :: nCellsI=6,  &
                        nCellsJ=6,  &
                        nCellsK=6
  integer, parameter :: nI=nCellsI, nJ=nCellsJ, nK=nCellsK
  integer, parameter, dimension(3) :: nCells=(/nCellsI,nCellsJ,nCellsK/)
  
  !\
  ! number of ghostcells. KEEP FIXED AT 2!
  !/
  integer, parameter :: gcn=2
  
end module ModSize
