!^CFG COPYRIGHT UM
subroutine specify_initial_refinement(refb, lev)
  use ModSize
  use ModMain, ONLY : &
       UseUserSpecifyRefinement,&
       body1,UseRotatingBc,unusedBLK
  use ModGeometry, ONLY : XyzMin_D,XyzMax_D,XyzStart_BLK,&
       x_BLK, y_BLK, z_BLK, dx_BLK, dy_BLK,dz_BLK, &
       TypeGeometry,x1,x2,far_field_BCs_BLK
  use ModPhysics, ONLY : Rbody,Rcurrents
  use ModAMR, ONLY : InitialRefineType
  use ModNumConst
  use ModUser, ONLY: user_specify_initial_refinement
  implicit none

  logical, intent(out) :: refb(nBLK)
  integer, intent(in) :: lev

  integer :: iBLK
  real :: xxx,yyy,zzz,RR, xxPoint,yyPoint,zzPoint
  real :: minRblk,maxRblk,xx1,xx2,yy1,yy2,zz1,zz2, tmpminRblk,tmpmaxRblk
  real :: critx, critDX
  real :: x_off,x_off1,x_off2,curve_rad,curve_rad1,curve_rad2
  real :: mach,frac,bound1,bound2
  real :: lscale,rcyl, minx, miny, minz, maxx, maxy, maxz

  real,parameter::cRefinedTailCutoff=(cOne-cQuarter)*(cOne-cEighth)
  logical::IsFound
  character(len=*), parameter :: NameSub='specify_initial_refinement'
  logical :: oktest, oktest_me

  !----------------------------------------------------------------------------
  call set_oktest('initial_refinement',oktest,oktest_me)

  refb = .false.

  do iBLK = 1,nBLK
     if (.not. unusedBLK(iBLK)) then
        ! Block min, max radius values
        xx1 = cHalf*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
        xx2 = cHalf*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
        yy1 = cHalf*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
        yy2 = cHalf*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
        zz1 = cHalf*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
        zz2 = cHalf*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))
        
        select case(TypeGeometry)                              
        case('cartesian')                                      
           ! Block center coordinates
           xxx = cHalf*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK)) 
           yyy = cHalf*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
           zzz = cHalf*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))         
           RR = sqrt( xxx*xxx + yyy*yyy + zzz*zzz )
           minRblk = sqrt(&
             minmod(xx1,xx2)**2 + minmod(yy1,yy2)**2 + minmod(zz1,zz2)**2)
           maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
                (max(abs(yy1),abs(yy2)))**2 + &
                (max(abs(zz1),abs(zz2)))**2)
           if(body1.and.maxRblk<rBody)CYCLE
        case('spherical')                                          
           minRblk = XyzStart_BLK(1,iBLK)-cHalf*dx_BLK(iBLK)            
           maxRblk = minRblk+nI*dx_BLK(iBLK)                            
           RR=cHalf*(minRblk+maxRblk)                                   
           xxx=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                cos(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           yyy=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                sin(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           zzz=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))
        case('spherical_lnr')                                          
           minRblk =exp( XyzStart_BLK(1,iBLK)-cHalf*dx_BLK(iBLK))            
           maxRblk =exp(nI*dx_BLK(iBLK))*minRblk                            
           RR=cHalf*(minRblk+maxRblk)                                   
           xxx=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                cos(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           yyy=RR*cos(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))*& 
                sin(XyzStart_BLK(2,iBLK)+cHalf*(nJ-1)*dy_BLK(iBLK))       
           zzz=RR*sin(XyzStart_BLK(3,iBLK)+cHalf*(nK-1)*dz_BLK(iBLK))
        case default
           minRblk=-cOne
           maxRblk=-cOne
           RR=-cOne
           xxx=-cOne
           yyy=-cOne
           zzz=-cOne
           !In case of an arbitrary geometry the variables above
           !are initialized but they are meaningless
        end select                                                

        select case (InitialRefineType)
        case ('none')
           ! Refine no blocks
        case ('all')
           ! Refine all used blocks
           refb(iBLK) = .true.
        case('outerboundary')
           refb(iBLK) =  far_field_BCs_BLK(iBLK)

        case ('3Dbodyfocus')
           ! Refine, focusing on body
           if (maxRblk > Rbody) then
              if (lev <= 2) then
                 ! Refine all blocks first two times through
                 refb(iBLK) = .true.
              else
                 ! Refine blocks intersecting body
                 if (minRblk < 1.5*Rbody) refb(iBLK) = .true.
              end if
           end if

        case('xplanefocus')
           ! concentrate refinement around x=0 plane
           if (lev <= 2) then
              ! Refine all blocks first two times through
              refb(iBLK) = .true.
           else
              ! Refine blocks intersecting x=0 plane
              if(xx1<=0.and.0<=xx2)refb(iBLK) = .true.
           end if

        case ('spherefocus')
           ! Refine, focusing spherically - like
           ! for a body, but can be used when there is
           ! no body
           if (lev <= 2) then
              ! Refine all blocks first two times through
              refb(iBLK) = .true.
           else
              ! Refine blocks intersecting R=4.5 sphere
              if (minRblk < 4.5) refb(iBLK) = .true.
           end if

        case default
           IsFound=.false.
           if (UseUserSpecifyRefinement) &
                call user_specify_initial_refinement(iBLK,refb(iBLK),lev,dx_BLK(iBLK), &
                xxx,yyy,zzz,RR,minx,miny,minz,minRblk,maxx,maxy,maxz,maxRblk,IsFound)
           if(.not.IsFound) &
           call stop_mpi(NameSub//' ERROR: unknown InitialRefineType='// &
                trim(InitialRefineType)//'!!!')
        end select
     end if
  end do

  call specify_area_refinement(refb)

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

end subroutine specify_initial_refinement

!==============================================================================

subroutine specify_area_refinement(DoRefine_B)

  !DESCRIPTION:
  ! Set DoRefine_B to .true. for blocks touching the predefined areas
  ! if the area has a finer resolution than the block

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: MaxBlock, nBlock, nBlockMax, nI, nJ, nK, UnusedBlk, &
       BlkTest, Test_String
  use ModAMR,      ONLY: nArea, Area_I
  use ModGeometry, ONLY: dx_BLK, UseCovariant
  use ModNodes,    ONLY: NodeX_NB, NodeY_NB, NodeZ_NB

  implicit none

  !INPUT/OUTPUT ARGUMENTS:
  logical, intent(inout) :: DoRefine_B(MaxBlock)

  !LOCAL VARIABLES:
  logical :: DoRefine
  integer, parameter :: nCorner = 8
  real     :: CornerOrig_DI(3, nCorner), Corner_DI(3, nCorner)
  real     :: DistMin_D(3), DistMax_D(3), Radius1Sqr
  integer  :: i, j, k, iDim, iArea, iBlock, iCorner
  real     :: CurrentResolution

  ! These variables are needed for generalized coordinates
  real :: Xyz_D(3)
  real, dimension(nI+1,nJ+1,nK+1):: x_N, y_N, z_N, R2_N

  character(len=*), parameter :: NameSub = 'specify_area_refinement'

  logical :: DoTest, DoTestMe, DoTestBlock
  !---------------------------------------------------------------------------
  if(nArea <= 0) RETURN


  call set_oktest(NameSub,DoTest,DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' nArea, nBlock, nBlockMax, MaxBlock=',&
       nArea, nBlock, nBlockMax, MaxBlock

  BLOCK: do iBlock = 1, nBlockMax

     if( UnusedBlk(iBlock) ) CYCLE BLOCK

     DoTestBlock = DoTestMe .and. iBlock == BlkTest

     ! No need to check block if it is to be refined already
     if(DoRefine_B(iBlock)) CYCLE BLOCK

     CurrentResolution = dx_BLK(iBlock)

     if(DoTestBlock)write(*,*)NameSub,' CurrentResolution=',CurrentResolution

     iCorner = 0
     do k=1,nK+1,nK; do j=1,nJ+1,nJ; do i=1,nI+1,nI
        iCorner = iCorner+1
        CornerOrig_DI(1,iCorner) = NodeX_NB(i,j,k,iBlock)
        CornerOrig_DI(2,iCorner) = NodeY_NB(i,j,k,iBlock)
        CornerOrig_DI(3,iCorner) = NodeZ_NB(i,j,k,iBlock)
     end do; end do; end do

     AREA: do iArea = 1, nArea

        if(DoTestBlock)write(*,*)NameSub,' iArea,Name,Resolution=',&
             iArea, ' ',trim(Area_I(iArea) % Name),' ',&
             Area_I(iArea) % Resolution

        ! No need to check area if block is finer than area resolution
        if(Area_I(iArea) % Resolution >= CurrentResolution) CYCLE AREA

        ! Check if area refines the whole domain
        if(Area_I(iArea) % Name == 'all')then
           DoRefine_B(iBlock) = .true.
           CYCLE AREA
        endif

        ! Shift corner coordinates to the center of area
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = &
                CornerOrig_DI(:,iCorner) - Area_I(iArea) % Center_D
        end do

        ! Rotate corners into the orientation of the area if required
        if(Area_I(iArea) % DoRotate) &
             Corner_DI = matmul(Area_I(iArea) % Rotate_DD, Corner_DI)

        ! Normalize coordinates to the size of the area in all 3 directions
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = Corner_DI(:,iCorner) / Area_I(iArea) % Size_D
        end do

        ! Calculate maximum and minimum distances in all 3 directions
        do iDim = 1, 3
           DistMax_D(iDim) = maxval(abs(Corner_DI(iDim,:)))

           if( maxval(Corner_DI(iDim,:))*minval(Corner_DI(iDim,:)) <= 0.0)then
              ! The block covers the center point in this dimension
              DistMin_D(iDim) = 0.0
           else
              ! Select the point that is closer in this dimension
              DistMin_D(iDim) = minval(abs(Corner_DI(iDim,:)))
           end if
        end do

        ! This occurs multiple times
        Radius1Sqr = Area_I(iArea) % Radius1

        if(DoTestBlock)then
           write(*,*)NameSub,' DistMin_D=',DistMin_D
           write(*,*)NameSub,' DistMax_D=',DistMax_D
        end if

        ! Check if this area is intersecting with the block
        select case( Area_I(iArea) % Name)
        case('brick')
           DoRefine = all( DistMin_D < 1.0 )
        case('sphere')
           DoRefine = sum(DistMin_D**2) < 1.0
        case('shell')
           ! Check if block intersects with the enclosing sphere
           ! but it is not fully inside the inner sphere
           DoRefine = sum(DistMin_D**2)<1.0 .and. sum(DistMax_D**2)>Radius1Sqr
        case('cylinderx')
           DoRefine = DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0
        case('cylindery')
           DoRefine = DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0
        case('cylinderz')
           DoRefine = DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0
        case('ringx')
           ! Check if block intersects with the enclosing cylinder
           ! but it is not fully inside the inner cylinder
           DoRefine = DistMin_D(1) < 1.0 .and. sum(DistMin_D(2:3)**2) < 1.0 &
                .and. sum(DistMax_D(2:3)**2) > Radius1Sqr 
        case('ringy')
           DoRefine = DistMin_D(2) < 1.0 .and. sum(DistMin_D(1:3:2)**2) < 1.0 &
                .and. sum(DistMax_D(1:3:2)**2) > Radius1Sqr
        case('ringz')
           DoRefine = DistMin_D(3) < 1.0 .and. sum(DistMin_D(1:2)**2) < 1.0 &
                .and. sum(DistMax_D(1:2)**2) > Radius1Sqr 

        case default
           call stop_mpi(NameSub // &
                ' ERROR: Unknown NameArea = '//Area_I(iArea) % Name)

        end select

        if(DoTestBlock)write(*,*) NameSub,' DoRefine (from corners)=',DoRefine

        if(DoRefine) DoRefine_B(iBlock) = .true.

        if(.not.UseCovariant)then
           if(DoRefine_B(iBlock)) EXIT AREA
           CYCLE
        end if

        ! Covariant case
        if(.not.DoRefine_B(iBlock)) CYCLE AREA

        ! Check if the covariant block is really inside the area
        ! Check all nodes of the block
        do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
           Xyz_D(1) = NodeX_NB(i,j,k,iBlock)
           Xyz_D(2) = NodeY_NB(i,j,k,iBlock)
           Xyz_D(3) = NodeZ_NB(i,j,k,iBlock)

           ! Shift to area center
           Xyz_D = Xyz_D - Area_I(iArea) % Center_D

           ! Rotate into area coordinates
           if(Area_I(iArea) % DoRotate) &
                Xyz_D = matmul(Area_I(iArea) % Rotate_DD, Xyz_D)

           ! Rescale coordinates to the size of the area in all directions
           Xyz_D = Xyz_D / Area_I(iArea) % Size_D

           ! We only need the absolute values of the coordinates
           x_N(i,j,k) = abs(Xyz_D(1))
           y_N(i,j,k) = abs(Xyz_D(2))
           z_N(i,j,k) = abs(Xyz_D(3))
        end do; end do; end do

        select case( Area_I(iArea) % Name)
        case('brick')
           if( any( x_N<1.0 .and. y_N<1.0 .and. z_N<1.0) ) EXIT AREA
        case('sphere')
           if( any(x_N**2 + y_N**2 + z_N**2 < 1.0) ) EXIT AREA
        case('shell')
           R2_N = x_N**2 + y_N**2 + z_N**2
           if( any(R2_N < 1.0 .and. R2_N > Radius1Sqr )) EXIT AREA
        case('cylinderx')
           if( any(x_N < 1.0 .and. y_N**2+z_N**2 < 1.0 ) ) EXIT AREA
        case('cylindery')
           if( any(y_N < 1.0 .and. x_N**2+z_N**2 < 1.0 ) ) EXIT AREA
        case('cylinderz')
           if( any(z_N < 1.0 .and. x_N**2+y_N**2 < 1.0 ) ) EXIT AREA
        case('ringx')
           R2_N = y_N**2+z_N**2
           if( any(x_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr) ) EXIT AREA
        case('ringy')
           R2_N = x_N**2+z_N**2
           if( any(y_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr) ) EXIT AREA
        case('ringz')
           R2_N = x_N**2+y_N**2
           if( any(z_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr) ) EXIT AREA
        case default
           call stop_mpi(NameSub // &
                ' ERROR: Unknown NameArea = '//Area_I(iArea) % Name)

        end select

        ! The block is not inside
        DoRefine_B(iBlock) = .false.

        if(DoTestBlock)write(*,*)NameSub,' DoRefine=false for iArea=',iArea

     end do AREA

     if(DoTestBlock)write(*,*)NameSub,' DoRefine final=',DoRefine_B(iBlock)

  end do BLOCK

  if(DoTest)write(*,*)NameSub,' on iProc=',iProc,&
       ' number of selected blocks=',count(DoRefine_B)

end subroutine specify_area_refinement
