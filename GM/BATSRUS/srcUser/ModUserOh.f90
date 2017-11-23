!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModUser

  use BATL_lib, ONLY: &
       test_start, test_stop
  use ModUserEmpty,               &
       IMPLEMENTED2 => user_specify_initial_refinement

  include 'user_module.h' ! list of public methods

  real, parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'HELIOSPHERE, Sokolov'
  ! This version allows to refine the current sheet without increasing the
  ! number of blocks intersected by the inner boundary
contains
  !============================================================================



  subroutine user_specify_initial_refinement(iBlock,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
    use ModMain, ONLY:time_loop,nI,nJ,nK
    use ModAMR, ONLY:InitialRefineType
    use ModNumConst
    use ModAdvance, ONLY:State_VGB,Bx_,By_,Bz_
    use ModGeometry
    use ModPhysics, ONLY:rBody
    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBlock

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real::BDotRMin,BDotRMax
    integer::i,j,k
    logical:: DoTest
    character(len=*), parameter:: NameSub = 'user_specify_initial_refinement'
    !--------------------------------------------------------------------------
    call test_start(NameSub, DoTest, iBlock)
    select case (InitialRefineType)
    case('coupledhelio')
       if(.not.time_loop)then
          ! refine to have resolution not worse 4.0 and
          ! refine the body intersecting blocks
          refineBlock=minR<=rBody.or.CellSize_DB(x_,iBlock)>4.01
       elseif(CellSize_DB(x_,iBlock)<0.40.or.far_field_BCs_BLK(iBlock))then
          refineBlock=.false. ! Do not refine body or outer boundary
       else
          ! refine heliosheath
          BDotRMin=0.0
          do k=0,nK+1;do j=1,nJ
             BDotRMin=min( BDotRMin,minval(&
                  State_VGB(Bx_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(x_,1:nI,j,k,iBlock)+&
                  State_VGB(By_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(y_,1:nI,j,k,iBlock)+&
                  State_VGB(Bz_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(z_,1:nI,j,k,iBlock)))
          end do;end do
          BDotRMax=0.0
          do k=0,nK+1;do j=1,nJ
             BDotRMax=max( BDotRMax,maxval(&
                  State_VGB(Bx_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(x_,1:nI,j,k,iBlock)+&
                  State_VGB(By_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(y_,1:nI,j,k,iBlock)+&
                  State_VGB(Bz_,1:nI,j,k,iBlock)*&
                  Xyz_DGB(z_,1:nI,j,k,iBlock)))
          end do;end do
          refineBlock=BDotRMin<-1e-6.and.&
               BDotRMax>1e-6
       end if
       found=.true.
    end select
    call test_stop(NameSub, DoTest, iBlock)
  end subroutine user_specify_initial_refinement
  !============================================================================
end module ModUser
!==============================================================================

