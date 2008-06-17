!^CFG COPYRIGHT UM


!========================================================================

module ModUser
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_set_boundary_cells,&
       IMPLEMENTED2 => user_specify_initial_refinement

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Sokolov'
  !This version allows to refine the current sheet without increasing the
  !number of blocks intersected by the inner boundary
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
  !----------------------------------------------------------
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
    use ModMain,ONLY:time_loop,nI,nJ,nK
    use ModAMR,ONLY:InitialRefineType
    use ModNumConst
    use ModAdvance,ONLY:State_VGB,Bx_,By_,Bz_
    use ModGeometry
    use ModPhysics,ONLY:rBody
    implicit none
    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real::BDotRMin,BDotRMax
    integer::i,j,k
    !-------------------------------------------------------------------
    select case (InitialRefineType)
    case('coupledhelio')
       if(.not.time_loop)then
          !refine to have resolution not worse 4.0 and
          !refine the body intersecting blocks
          refineBlock=minR<=rBody.or.dx_BLK(iBLK)>4.01
       elseif(dx_BLK(iBLK)<0.40.or.far_field_BCs_BLK(iBLK))then
          refineBlock=.false. !Do not refine body or outer boundary
       else
          !refine heliosheath
          BDotRMin=cZero
          do k=0,nK+1;do j=1,nJ
             BDotRMin=min( BDotRMin,minval(&
                  State_VGB(Bx_,1:nI,j,k,iBLK)*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK)*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK)*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          BDotRMax=cZero
          do k=0,nK+1;do j=1,nJ
             BDotRMax=max( BDotRMax,maxval(&
                  State_VGB(Bx_,1:nI,j,k,iBLK)*&
                  x_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(By_,1:nI,j,k,iBLK)*&
                  y_BLK(1:nI,j,k,iBLK)+&
                  State_VGB(Bz_,1:nI,j,k,iBLK)*&
                  z_BLK(1:nI,j,k,iBLK)))
          end do;end do
          refineBlock=BDotRMin<-cTiny.and.&
               BDotRMax>cTiny
       end if
       found=.true.
    end select
  end subroutine user_specify_initial_refinement
end module ModUser

