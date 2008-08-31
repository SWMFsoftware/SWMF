!^CFG COPYRIGHT UM
subroutine specify_refinement(DoRefine_B)

  !DESCRIPTION:
  ! Set DoRefine_B to .true. for blocks touching the predefined areas
  ! if the area has a finer resolution than the block

  use ModProcMH,   ONLY: iProc
  use ModMain,     ONLY: MaxBlock, nBlock, nBlockMax, nI, nJ, nK, UnusedBlk, &
       BlkTest, Test_String, r_, Phi_, Theta_, x_, y_, z_,UseB0
  use ModAMR,      ONLY: nArea, AreaType, Area_I, lNameArea
  use ModGeometry, ONLY: dx_BLK, UseCovariant, x_BLK, y_BLK, z_BLK, &
       TypeGeometry
  use ModNodes,    ONLY: NodeX_NB, NodeY_NB, NodeZ_NB
  use ModUser,     ONLY: user_specify_refinement

  ! Needed for the 'currentsheet' area type only
  use ModAdvance,  ONLY: State_VGB, Bx_, By_, Bz_, B0_DGB
  use ModGeometry, ONLY: far_field_BCs_BLK
  use ModNumConst, ONLY: cTiny, cRadToDeg

  implicit none

  !OUTPUT ARGUMENTS:
  logical, intent(out) :: DoRefine_B(MaxBlock)

  !LOCAL VARIABLES:
  type(AreaType) :: Area
  character(len=lNameArea) :: NameArea

  logical :: DoRefine, IsSpecialArea
  integer, parameter :: nCorner = 8
  real     :: CornerOrig_DI(3, nCorner), Corner_DI(3, nCorner)
  real     :: DistMin_D(3), DistMax_D(3), Radius1Sqr
  integer  :: i, j, k, iDim, iArea, iBlock, iCorner
  real     :: CurrentResolution

  ! These variables are needed for generalized coordinates
  real :: Xyz_D(3)
  real, dimension(nI+1,nJ+1,nK+1):: x_N, y_N, z_N, R2_N

  ! Needed for the 'currentsheet'
  real :: rDotB_G(nI,nJ,0:nK+1)

  character(len=*), parameter :: NameSub = 'specify_refinement'

  logical :: DoTest, DoTestMe, DoTestBlock
  !---------------------------------------------------------------------------
  if(nArea <= 0) RETURN

  call set_oktest(NameSub,DoTest,DoTestMe)

  if(DoTestMe)write(*,*)NameSub,' nArea, nBlock, nBlockMax, MaxBlock=',&
       nArea, nBlock, nBlockMax, MaxBlock

  DoRefine_B = .false.

  BLOCKLOOP: do iBlock = 1, nBlockMax

     if( UnusedBlk(iBlock) ) CYCLE BLOCKLOOP

     DoTestBlock = DoTestMe .and. iBlock == BlkTest

     CurrentResolution = dx_BLK(iBlock)

     if(DoTestBlock)write(*,*)NameSub,' CurrentResolution=',CurrentResolution

     iCorner = 0
     do k=1,nK+1,nK; do j=1,nJ+1,nJ; do i=1,nI+1,nI
        iCorner = iCorner+1
        CornerOrig_DI(1,iCorner) = NodeX_NB(i,j,k,iBlock)
        CornerOrig_DI(2,iCorner) = NodeY_NB(i,j,k,iBlock)
        CornerOrig_DI(3,iCorner) = NodeZ_NB(i,j,k,iBlock)

        if(DoTestBlock)write(*,*)NameSub,' iCorner, CornerOrig_D=',&
             iCorner, CornerOrig_DI(:,iCorner)
     end do; end do; end do

     AREALOOP: do iArea = 1, nArea
        
        Area = Area_I(iArea)
        NameArea = Area % Name

        if(DoTestBlock)write(*,*)NameSub,' iArea,Name,Resolution=',&
             iArea, ' ',trim(NameArea),' ',Area % Resolution

        ! No need to check area if block is finer than area resolution
        if(Area % Resolution >= CurrentResolution) CYCLE AREALOOP

        ! Treat special cases first
        IsSpecialArea = .true.
        select case(NameArea)
        case('all')
           DoRefine_B(iBlock) = .true.
        case('user')
           call user_specify_refinement(iBlock, iArea, DoRefine_B(iBlock))
        case('currentsheet')

           ! Calculate BdotR including ghost cells in all directions
           if(UseB0)then
              do k=0, nK+1; do j=1, nJ; do i=1, nI
                 rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
                      * (B0_DGB(x_,i,j,k,iBlock) + State_VGB(Bx_,i,j,k,iBlock)) &
                      +              y_BLK(i,j,k,iBlock)   &
                      * (B0_DGB(y_,i,j,k,iBlock) + State_VGB(By_,i,j,k,iBlock)) &
                      +              z_BLK(i,j,k,iBlock)   &
                      * (B0_DGB(z_,i,j,k,iBlock) + State_VGB(Bz_,i,j,k,iBlock))
              end do; end do; end do
           else
              do k=0, nK+1; do j=1, nJ; do i=1, nI
                 rDotB_G(i,j,k) = x_BLK(i,j,k,iBlock)   &
                      * State_VGB(Bx_,i,j,k,iBlock) &
                      +           y_BLK(i,j,k,iBlock)   &
                      * State_VGB(By_,i,j,k,iBlock) &
                      +           z_BLK(i,j,k,iBlock)   &
                      * State_VGB(Bz_,i,j,k,iBlock)
              end do; end do; end do
           end if
           DoRefine_B(iBlock) = &
                maxval(rDotB_G) > cTiny .and. minval(rDotB_G) < -cTiny

        case default
           IsSpecialArea = .false.
        end select

        if(IsSpecialArea)then
           if(DoRefine_B(iBlock)) EXIT AREALOOP
           EXIT AREALOOP
        end if

        ! Check if it is a brick in the generalized coordinates
        if(NameArea == "brick_gen" .and. TypeGeometry /= 'cartesian')then
           ! Convert corners to generalized coordinates
           do iCorner = 1, nCorner
              call xyz_to_gen(CornerOrig_DI(:,iCorner), Corner_DI(:,iCorner))
           end do

           ! Convert angles to degrees and ln(r) to r
           ! Make sure that phi=360 is not 0 but really 360
           select case(TypeGeometry)
           case('spherical')
              Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
              where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
              where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
              Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg

           case('spherical_lnr')
              Corner_DI(r_,:)     = exp(Corner_DI(r_,:))
              Corner_DI(Phi_,:)   = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
              where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
              where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
              Corner_DI(Theta_,:) = Corner_DI(Theta_,:)*cRadToDeg

           case('cylindrical', 'axial_torus')
              Corner_DI(Phi_,:) = modulo(Corner_DI(Phi_,:)*cRadToDeg, 360.0)
              where(Corner_DI(Phi_,3:4) < cTiny) Corner_DI(Phi_,3:4) = 360.0
              where(Corner_DI(Phi_,7:8) < cTiny) Corner_DI(Phi_,7:8) = 360.0
           end select
        else
           Corner_DI = CornerOrig_DI
        end if

        ! Shift corner coordinates to the center of area
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = Corner_DI(:,iCorner) - Area % Center_D
        end do

        ! Rotate corners into the orientation of the area if required
        if(Area % DoRotate) Corner_DI = matmul(Area % Rotate_DD, Corner_DI)

        ! Normalize coordinates to the size of the area in all 3 directions
        do iCorner = 1, nCorner
           Corner_DI(:,iCorner) = Corner_DI(:,iCorner) / Area % Size_D
        end do

        ! Calculate maximum and minimum distances in all 3 directions
        ! Avoid rounding errors if possible
        do iDim = 1, 3
           DistMax_D(iDim) = (1-cTiny)*maxval(abs(Corner_DI(iDim,:)))

           if( maxval(Corner_DI(iDim,:))*minval(Corner_DI(iDim,:)) <= 0.0)then
              ! The block covers the center point in this dimension
              DistMin_D(iDim) = 0.0
           else
              ! Select the point that is closer in this dimension
              DistMin_D(iDim) = (1+cTiny)*minval(abs(Corner_DI(iDim,:)))
           end if
        end do

        ! This occurs multiple times
        Radius1Sqr = (Area % Radius1)**2

        if(DoTestBlock)then
           write(*,*)NameSub,' DistMin_D=',DistMin_D
           write(*,*)NameSub,' DistMax_D=',DistMax_D
        end if

        ! Check if this area is intersecting with the block
        select case( NameArea)
        case('brick', 'brick_gen')
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
                ' ERROR: Unknown NameArea = '//NameArea)

        end select

        if(DoTestBlock)write(*,*) NameSub,' DoRefine (from corners)=',DoRefine

        if(DoRefine) DoRefine_B(iBlock) = .true.

        if(NameArea == 'brick_gen' .or. .not. UseCovariant )then
           if(DoRefine_B(iBlock)) EXIT AREALOOP
           CYCLE AREALOOP
        end if

        ! Covariant case
        if(.not.DoRefine_B(iBlock)) CYCLE AREALOOP

        ! Check if the covariant block is really inside the area
        ! Check all nodes of the block
        do k=1,nK+1; do j=1,nJ+1; do i=1,nI+1
           Xyz_D(1) = NodeX_NB(i,j,k,iBlock)
           Xyz_D(2) = NodeY_NB(i,j,k,iBlock)
           Xyz_D(3) = NodeZ_NB(i,j,k,iBlock)

           ! Shift to area center
           Xyz_D = Xyz_D - Area % Center_D

           ! Rotate into area coordinates
           if(Area % DoRotate) &
                Xyz_D = matmul(Area % Rotate_DD, Xyz_D)

           ! Rescale coordinates to the size of the area in all directions
           Xyz_D = Xyz_D / Area % Size_D

           ! We only need the absolute values of the coordinates
           x_N(i,j,k) = abs(Xyz_D(1))
           y_N(i,j,k) = abs(Xyz_D(2))
           z_N(i,j,k) = abs(Xyz_D(3))
        end do; end do; end do

        select case( NameArea)
        case('brick')
           if( any( x_N<1.0 .and. y_N<1.0 .and. z_N<1.0) ) EXIT AREALOOP
        case('sphere')
           if( any(x_N**2 + y_N**2 + z_N**2 < 1.0) ) EXIT AREALOOP
        case('shell')
           R2_N = x_N**2 + y_N**2 + z_N**2
           if( any(R2_N < 1.0 .and. R2_N > Radius1Sqr )) EXIT AREALOOP
        case('cylinderx')
           if( any(x_N < 1.0 .and. y_N**2+z_N**2 < 1.0 ) ) EXIT AREALOOP
        case('cylindery')
           if( any(y_N < 1.0 .and. x_N**2+z_N**2 < 1.0 ) ) EXIT AREALOOP
        case('cylinderz')
           if( any(z_N < 1.0 .and. x_N**2+y_N**2 < 1.0 ) ) EXIT AREALOOP
        case('ringx')
           R2_N = y_N**2+z_N**2
           if(any(x_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) EXIT AREALOOP
        case('ringy')
           R2_N = x_N**2+z_N**2
           if(any(y_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) EXIT AREALOOP
        case('ringz')
           R2_N = x_N**2+y_N**2
           if(any(z_N<1.0 .and. R2_N<1.0 .and. R2_N>Radius1Sqr)) EXIT AREALOOP
        case default
           call stop_mpi(NameSub //' ERROR: Unknown NameArea = '//NameArea)
        end select

        ! The block is not inside
        DoRefine_B(iBlock) = .false.

        if(DoTestBlock)write(*,*)NameSub,' DoRefine=false for iArea=',iArea

     end do AREALOOP

     if(DoTestBlock)write(*,*)NameSub,' DoRefine final=',DoRefine_B(iBlock)

  end do BLOCKLOOP

  if(DoTest)write(*,*)NameSub,' on iProc=',iProc,&
       ' number of selected blocks=',count(DoRefine_B)

end subroutine specify_refinement
