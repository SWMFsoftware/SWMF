!^CFG COPYRIGHT UM


!========================================================================

module ModUser
  use ModNumConst, ONLY: cHalf,cTwo,cThree,&
       cFour,cE1,cHundred,cHundredth,cZero,&
       cOne
  use ModMain,     ONLY: UseUserB0,UseUserHeating
  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK


  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Manchester, Roussev'

contains


  subroutine user_set_boundary_cells(iBLK)

    ! Set the boundary cell information IsBoundaryCell_GI(:,:,:,ExtraBc_) 
    ! for a sphere of radius 1 around the origin.
    ! Allow resolution change.

    use ModGeometry
    use ModBoundaryCells,ONLY:SaveBoundaryCells
    use ModPhysics,ONLY:rBody
    implicit none
    integer, intent(in):: iBLK

    IsBoundaryCell_GI(:,:,:,ExtraBc_) = R_BLK(:,:,:,iBLK)<rBody
    if(SaveBoundaryCells)return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')
  end subroutine user_set_boundary_cells
end module ModUser

