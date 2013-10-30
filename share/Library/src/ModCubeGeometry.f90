!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModCubeGeometry
  implicit none
  SAVE
  !=================ARRAYS FOR A CUBE====================================!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! z-axis|
  !       |  7----------8
  !       | /|         /|
  !       |/ |        / |
  !       5----------6  |
  !       |  |   _   |  |
  !       |  |   /|y-axi$
  !       |  |  /    |  |
  !       |  | /     |  |
  !       |  |/      |  |
  !       |  3----------4
  !       | /        | /
  !       |/         |/
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. The functions used in more than one place,
  ! therefore we delegate them here.   

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir+1
  !(3) is connected with the given one by the edge of direction iDir+2
  !(4) is connected to the given one by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       1, 3, 5, 7,   1, 5, 2, 6,   1, 2, 3, 4, & 
       2, 4, 6, 8,   2, 6, 1, 5,   2, 1, 4, 3, &
       3, 1, 7, 5,   3, 7, 4, 8,   3, 4, 1, 2, &
       4, 2, 8, 6,   4, 8, 3, 7,   4, 3, 2, 1, &
       5, 7, 1, 3,   5, 1, 6, 2,   5, 6, 7, 8, &
       6, 8, 2, 4,   6, 2, 5, 1,   6, 5, 8, 7, &
       7, 5, 3, 1,   7, 3, 8, 4,   7, 8, 5, 6, &
       8, 6, 4, 2,   8, 4, 7, 3,   8, 7, 6, 5  &
       /), (/4,3,8/))

  !Vertexes (enumerated by the first undex), which form 
  !the face of direction iDir (second index) and does not 
  !include the given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) is connected to the given one by the edge of direction iDir
  !(2) is connected to (1) by the edge of direction iDir+1
  !(3) is connected to (1) by the edge of direction iDir+2
  !(4) is connected to  by the face diagonal of direction iDir
  integer, dimension(4,3,8), parameter:: iOppositeFace_IDI = &
       reshape((/&   ! yz face   ! xz face      ! xy face ! 
       2, 4, 6, 8,   3, 7, 4, 8,   5, 6, 7, 8, & 
       1, 3, 5, 7,   4, 8, 3, 7,   6, 5, 8, 7, &
       4, 2, 8, 6,   1, 5, 2, 6,   7, 8, 5, 6, &
       3, 1, 7, 5,   2, 6, 1, 5,   8, 7, 6, 5, &
       6, 8, 2, 4,   7, 3, 8, 4,   1, 2, 3, 4, &
       5, 7, 1, 3,   8, 4, 7, 3,   2, 1, 4, 3, &
       8, 6, 4, 2,   5, 1, 6, 2,   3, 4, 1, 2, &
       7, 5, 3, 1,   6, 2, 5, 1,   4, 3, 2, 1  &
       /), (/4,3,8/))
    !Number of the vertex connected by 
    !the edge of direction iDir (second index) 
    !with the given vertex (first index)
    integer, dimension(8,3), parameter:: iEdge_ID = &
         reshape((/&   !Number of the connected vertex
         2 ,1, 4, 3, 6, 5, 8, 7,        & !Edge x
         3, 4, 1, 2, 7, 8, 5, 6,        & !Edge y
         5, 6, 7, 8, 1, 2, 3, 4         & !Edge z
         /),(/8,3/))
  !=================ARRAYS FOR A RECTANGLE==============!
  !For a cubic stencil enumerated as follows:
  !
  !       ^
  ! y-axis|
  !       |   
  !       3----------4  
  !       |          |  
  !       |          !
  !       |          |  
  !       |          |  
  !       |          |  
  !       |          !
  !       |          |
  !       |          |
  !       1----------2------> x-axis
  !       
  ! we provide several functions characterizing its geometry in terms of 
  ! the grid point numbers. The functions used in more than one place,
  ! therefore we delegate them here. 
  

  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) including the 
  !given vertex (the third index). 
  !When the first index equals 1,2,3,4, the vertex, accordingly: 
  !(1) coincides with the given one
  !(2) is connected with the given one by the edge of direction iDir
  integer, dimension(2,2,4), parameter:: iSide_IDI = &
       reshape((/&   ! x side   ! y side ! 
       1, 2,    1, 3, & !Vertex 1
       2, 1,    2, 4, & !Vertex 2
       3, 4,    3, 1, & !Vertex 3
       4, 3,    4, 2  & !Vertex 4
       /), (/2,2,4/))

  !Vertexes (enumerated by the first undex), which form 
  !the side of direction iDir (second index) and does not include the
  !given vertex (the third index). 
  !When the first index equals 1,2 the vertex, accordingly: 
  !(1) is connected to the given one by the side of direction 1+mod(iDir,2)
  !(2) is connected to (1) by the edge of direction iDir
  integer, dimension(2,2,4), parameter:: iOppositeSide_IDI = &
       reshape((/&   ! x side   ! y side ! 
       3, 4,    2, 4, & !Vertex 1
       4, 3,    1, 3, & !Vertex 2
       1, 2,    4, 2, & !Vertex 3
       2, 1,    3, 1  & !Vertex 4
       /), (/2,2,4/))
  !\
  !Array used to sort corners: 
  !(1) assign a type to it (Case_) 
  !(2) assign a grid point of the stencil chosen as the first one (Grid_)
  !(3) assign a direction, if the configuration can be characterized (Dir_)
  !by direction

  logical:: DoInit = .true.
  integer, parameter:: Grid_=1, Dir_=2, Case_=3
  integer:: iSortStencil3_II(Grid_:Case_,0:258) = 0
  integer:: iSortStencil2_II(Grid_:Case_, 0:15) = 0
  
  !\
  ! Different cases of stencil (will be stored in Case_ column
  !/
  integer, parameter:: Uniform_ = 0, Face_ = 1, Edge_ = 2
  integer, parameter:: FiveTetrahedra_ = 3
  integer, parameter:: OneFine_        = 4, OneCoarse_      = 5
  integer, parameter:: FineMainDiag_   = 6, FineFaceDiag_   = 7
  integer, parameter:: CoarseMainDiag_ = 8, CoarseFaceDiag_ = 9
  integer, parameter:: FineEdgePlusOne_ = 10, ThreeFineOnFace_ =11
  integer, parameter:: CoarseEdgePlusOne_ = 12, ThreeCoarseOnFace_=13
  integer, parameter:: ThreeCoarseOnFacePlusOne_ = 14, CoarseChain_   = 15
  integer, parameter:: Transition2Edge_ = 256, Transition2Corner_ = 257
  integer, parameter:: TransitionJunction_ = 258
  !\
  ! Analogous for 2D
  !/
  integer, parameter:: Trapezoid_ = 1
  integer, parameter:: Rhombus_ = 3
  integer, parameter:: &
        Coarse_  = 0,               &
        Fine_    = 1

  !Parameters to enumerate faces or face diagonal
  integer, parameter:: Xy_ = 3, Xz_ = 2, Yz_ = 1
  
  !Number of the vertex connected by 
  !the face diagonal across the face of direction iDir (second index) 
  !with the given vertex (first index)
  integer, dimension(8,3), parameter:: iFaceDiag_ID = &
       reshape((/&   !Number of the connected vertex
       7 ,8, 5, 6, 3, 4, 1, 2,        & !Face yz
       6, 5, 8, 7, 2, 1, 4, 3,        & !Face xz
       4, 3, 2, 1, 8, 7, 6, 5         & !Face xy
       /),(/8,3/))
  
  !Number of the vertex connected by 
  !the main diagonal with the given vertex (index)
  integer, dimension(8), parameter:: iMainDiag_I = &
       (/ 8, 7, 6, 5, 4, 3, 2, 1 /)
contains
  !==========================
  subroutine init_sort_stencil
    integer:: iCase, iGrid, iLoc, iDir, iMisc
    integer:: iLevel_I(8)
    integer :: nDim = 3, nGrid = 8
    integer, parameter:: x_ = 1, y_ = 2, z_ = 3

    !------------------------
    nDim = 3; nGrid = 8
    iSortStencil3_II = 0

    DoInit = .false. !Do this only once
    !Zero and 255 cases are both uniform. !2 cases, 2 totally, 254 left
    iSortStencil3_II(Grid_:Dir_,0:255:255) = 1
    iSortStencil3_II(Case_,0:255:255) = Uniform_
    CASE:do iCase = 1, 254
       !\
       ! Generate 'binary' iLevel_I from iCase
       !/
       iLevel_I=0
       iGrid = 1
       iMisc = iCase
       do while(iMisc > 0)
          iLevel_I(iGrid) = mod(iMisc,2)
          iMisc = (iMisc -  iLevel_I(iGrid))/2
          iGrid = iGrid +1
       end do
       !\
       ! Sort out faces 6 cases, 8 totally, 248 left
       !/
       do iDir = 1, 3
          if(all(iLevel_I(iFace_IDI(:,iDir,1))==Fine_).and.&
               all(iLevel_I(iOppositeFace_IDI(:,iDir,1))==Coarse_))then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir !  2*iDir -1
             iSortStencil3_II(Case_,iCase) = Face_
             CYCLE CASE
          elseif(all(iLevel_I(iFace_IDI(:,iDir,1))==Coarse_).and.&
               all(iLevel_I(iOppositeFace_IDI(:,iDir,1))==Fine_))then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir !2*iDir 
             iSortStencil3_II(Case_,iCase) = Face_
             CYCLE CASE
          end if
       end do
       !\
       ! Sort out edges 30 cases, 38 totally, 218 left
       !/
       do iDir = 1, 3
          if(all(iLevel_I(iFace_IDI(:,iDir,1))==&
               iLevel_I(iOppositeFace_IDI(:,iDir,1)) ) )then
             iSortStencil3_II(Grid_,iCase) = 1
             iSortStencil3_II(Dir_, iCase) = iDir
             iSortStencil3_II(Case_,iCase) = Edge_
             CYCLE CASE
          end if
       end do
       !\
       ! Find configurations which split into five tetrahedra
       !   26 cases 64 totally 192 left
       !/
       do iGrid = 1,nGrid
          if(  all(iLevel_I(iEdge_ID(iGrid,:))==Coarse_).and.&
               all(iLevel_I(iFaceDiag_ID(iGrid,:))==Fine_))then
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iGrid
             iSortStencil3_II(Case_,iCase) = FiveTetrahedra_
             CYCLE CASE
          end if
       end do

       select case( count(iLevel_I==Fine_))
       case(1)                          ! 8 cases 72 totally 184 left 
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil3_II(Dir_,iCase) = 1
          iSortStencil3_II(Grid_,iCase) = iLoc
          iSortStencil3_II(Case_,iCase) = OneFine_
          CYCLE CASE 
       case(7)                          ! 8 cases 80 totally 176 left
          iLoc = minloc(iLevel_I,DIM=1)
          iSortStencil3_II(Dir_,iCase)  = 1
          iSortStencil3_II(Grid_,iCase) = iLoc
          iSortStencil3_II(Case_,iCase) = OneCoarse_
          CYCLE CASE
       case(2)
          iLoc = maxloc(iLevel_I,DIM=1)
          if(iLevel_I(iMainDiag_I(iLoc))==Fine_)then 

             ! 4 cases   84 totally 172 left
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iLoc
             iSortStencil3_II(Case_,iCase) = FineMainDiag_
             CYCLE CASE
          else                      ! 12 cases  96 totally 160 left 
             do iDir = Yz_,Xy_
                if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Fine_)then
                   iSortStencil3_II(Dir_,iCase)  = iDir
                   iSortStencil3_II(Grid_,iCase) = iLoc
                   iSortStencil3_II(Case_,iCase) = FineFaceDiag_
                   CYCLE CASE
                end if
             end do
          end if
       case(6)                       ! 4 cases 100 totally 156 left
          iLoc = minloc(iLevel_I,DIM=1)
          if(iLevel_I(iMainDiag_I(iLoc))==Coarse_)then  
             iSortStencil3_II(Dir_,iCase)  = 1
             iSortStencil3_II(Grid_,iCase) = iLoc
             iSortStencil3_II(Case_,iCase) = CoarseMainDiag_
             CYCLE CASE
          else                      ! 12 cases 112 totally 144 left
             do iDir = Yz_,Xy_
                if(iLevel_I(iFaceDiag_ID(iLoc, iDir))==Coarse_)then
                   iSortStencil3_II(Dir_,iCase) = iDir
                   iSortStencil3_II(Grid_,iCase) = iLoc
                   iSortStencil3_II(Case_,iCase) = CoarseFaceDiag_
                   CYCLE CASE
                end if
             end do
          end if
       case(3)          !         24 cases  136 totally 120 left
          do iGrid = 1, nGrid
             if(iLevel_I(iGrid)==Fine_)then
                if(iLevel_I(iMainDiag_I(iGrid))==Fine_)then
                   do iDir = Yz_,Xy_
                      if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Fine_)then
                         iSortStencil3_II(Dir_,iCase)  = iDir
                         iSortStencil3_II(Grid_,iCase) = iGrid
                         iSortStencil3_II(Case_,iCase) = FineEdgePlusOne_
                         CYCLE CASE
                      end if
                   end do
                else            ! 24 cases 160 totally 96 left
                   do iDir = 1,nDim
                      if(all(iLevel_I(iFace_IDI(1:3,iDir,iGrid))==Fine_))then 
                         iSortStencil3_II(Dir_,iCase)  = iDir
                         iSortStencil3_II(Grid_,iCase) = iGrid
                         iSortStencil3_II(Case_,iCase) = ThreeFineOnFace_
                         CYCLE CASE
                      end if
                   end do
                end if
             end if
          end do
       case(5)          !         24 cases 184 totally 72 left 
          do iGrid = 1, nGrid
             if(iLevel_I(iGrid)/=Coarse_)CYCLE
             if(iLevel_I(iMainDiag_I(iGrid))==Coarse_)then
                do iDir = Yz_,Xy_
                   if(iLevel_I(iFaceDiag_ID(iGrid, iDir))==Coarse_)then 
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = CoarseEdgePlusOne_
                      CYCLE CASE          
                   end if
                end do
             end if
             !                     24 cases 208 totally 48 left
             do iDir = 1,nDim
                if(all(iLevel_I(iFace_IDI(1:3, iDir, iGrid))==Coarse_))then
                   iSortStencil3_II(Dir_,iCase)  = iDir
                   iSortStencil3_II(Grid_,iCase) = iGrid
                   iSortStencil3_II(Case_,iCase) = ThreeCoarseOnFace_
                   CYCLE CASE
                end if
             end do
          end do
       case(4)                 
          do iGrid = 1,nGrid
             if(iLevel_I(iGrid)==Coarse_)then
                do iDir = 1,nDim    !   24 cases 232 totally 24 left 
                   if(all(&
                        iLevel_I(iOppositeFace_IDI(2:4,iDir,iGrid))==Coarse_&
                        ))then
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = &
                           ThreeCoarseOnFacePlusOne_
                      CYCLE CASE
                   end if
                   !
                   !\
                   ! 24 cases 256 totally 0 left
                   !/
                   if(iLevel_I(iEdge_ID(iGrid,iDir))==Coarse_.and.&
                        iLevel_I(iEdge_ID(iGrid,1 + mod(iDir,3)))==Coarse_&
                        .and.iLevel_I(iEdge_ID(iEdge_ID(iGrid,iDir),&
                        1 + mod(1 + iDir,3)))==Coarse_)then
                      iSortStencil3_II(Dir_,iCase)  = iDir
                      iSortStencil3_II(Grid_,iCase) = iGrid
                      iSortStencil3_II(Case_,iCase) = CoarseChain_
                      CYCLE CASE
                   end if
                end do
             end if
          end do
       end select
    end do CASE
    iSortStencil3_II(Grid_:Dir_,Transition2Edge_:TransitionJunction_)  = 1
    iSortStencil3_II(Case_,Transition2Edge_:TransitionJunction_)  = &
         (/Transition2Edge_, Transition2Corner_, TransitionJunction_/)
    !=================2 dimensional case=======
    nDim = 2; nGrid = 4
   
    !Zero and 15 cases are both uniform. !2 cases, 2 totally, 14 left
    iSortStencil2_II(Grid_:Dir_,0:15:15) = 1
    iSortStencil2_II(Case_,0:15:15) = Uniform_
    CASE2:do iCase = 1, 14
       !\
       ! Generate 'binary' iLevel_I from iCase
       !/
       iLevel_I=0
       iGrid = 1
       iMisc = iCase
       do while(iMisc > 0)
          iLevel_I(iGrid) = mod(iMisc,2)
          iMisc = (iMisc -  iLevel_I(iGrid))/2
          iGrid = iGrid +1
       end do
       !\
       ! Sort out faces 4 cases, 6 totally, 10 left
       !/
       do iDir = 1, nDim
          if(all(iLevel_I(iSide_IDI(:,iDir,1))==Fine_).and.&
               all(iLevel_I(iOppositeSide_IDI(:,iDir,1))==Coarse_))then
             iSortStencil2_II(Grid_,iCase) = 1
             iSortStencil2_II(Dir_, iCase) = iDir !  2*iDir -1
             iSortStencil2_II(Case_,iCase) = Trapezoid_
             CYCLE CASE2
          elseif(all(iLevel_I(iSide_IDI(:,iDir,1))==Coarse_).and.&
               all(iLevel_I(iOppositeSide_IDI(:,iDir,1))==Fine_))then
             iSortStencil2_II(Grid_,iCase) = 1
             iSortStencil2_II(Dir_, iCase) = iDir !2*iDir 
             iSortStencil2_II(Case_,iCase) = Trapezoid_
             CYCLE CASE2
          end if
       end do
       select case( count(iLevel_I(1:4)==Fine_))
       case(1)                          ! 4 cases 10 totally 6 left 
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil2_II(Dir_,iCase) = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = OneFine_
          CYCLE CASE2
       case(3)                          ! 4 cases 14 totally 2 left
          iLoc = minloc(iLevel_I(1:4),DIM=1)
          iSortStencil2_II(Dir_,iCase)  = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = OneCoarse_
          CYCLE CASE2
       case(2)                          !          
          iLoc = maxloc(iLevel_I,DIM=1)
          iSortStencil2_II(Dir_,iCase) = 1
          iSortStencil2_II(Grid_,iCase) = iLoc
          iSortStencil2_II(Case_,iCase) = Rhombus_
       end select
    end do CASE2
  end subroutine init_sort_stencil
end module ModCubeGeometry
