subroutine fix_axis_cells

  use ModProcMH, ONLY: iComm
  use ModMain, ONLY: nI, nJ, nK, nBlock, UnusedBlk, iTest, jTest, kTest, &
       BlkTest
  use ModAdvance, ONLY: nVar, State_VGB, Energy_GBI, rFixAxis, r2FixAxis
  use ModGeometry, ONLY: TypeGeometry, XyzMin_D, XyzMax_D, MinDxValue, &
       x_Blk, y_Blk, r_BLK, rMin_BLK, far_field_bcs_blk, vInv_CB,&
       r_to_gen
  use ModEnergy, ONLY: calc_energy_point
  use ModParallel, ONLY: NeiLBot, NeiLTop, NOBLK
  use ModMpi
  implicit none

  real, allocatable:: Buffer_VIII(:,:,:,:), SumBuffer_VIII(:,:,:,:)

  integer, parameter :: Sum_=1, SumXLeft_=2, SumXRight_=3, &
       SumYLeft_=4, SumYRight_=5, Geom_=6
  integer, parameter :: Volume_=1, SumX_=2, SumX2_=3
  integer, parameter :: North_=1, South_=2
  integer :: i, j, k, kMin, kMax, kOut, iBlock, iHemisphere, iR, nR, iError
  integer :: iVar, nAxisCell
  real :: r, x, y, Volume, InvVolume, SumX, SumXAvg, InvSumX2, dLeft, dRight
  real :: State_V(nVar), dStateDx_V(nVar), dStateDy_V(nVar)

  logical:: DoTest, DoTestMe
  character(len=*), parameter:: NameSub = 'fix_axis_cells'
  !--------------------------------------------------------------------------
  call set_oktest(NameSub, DoTest, DoTestMe)

  if(DoTestMe)then
     if(.not.UnusedBlk(BlkTest)) &
          write(*,*) NameSub,' initial state, energy=', &
          State_VGB(:,iTest,jTest,kTest,BlkTest), &
          Energy_GBI(iTest,jTest,kTest,BlkTest,:)
  end if

  if(TypeGeometry == 'cylindrical')then
     call fix_axis_cells_cyl
     RETURN
  end if

  ! Maximum number of cells along the axis
  nR = nint((XyzMax_D(1)-XyzMin_D(1))/MinDxValue)
  allocate(Buffer_VIII(nVar, nR, 1:Geom_, North_:South_), &
       SumBuffer_VIII(nVar, nR, 1:Geom_, North_:South_))

  Buffer_VIII    = 0.0
  SumBuffer_VIII = 0.0

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

     if(rMin_BLK(iBlock) < r2FixAxis) then 
        nAxisCell = 2
     elseif(rMin_BLK(iBlock) < rFixAxis) then 
        nAxisCell = 1
     else
        CYCLE
     end if

     ! Determine hemisphere
     if( NeiLTop(iBlock) == NOBLK )then
        iHemisphere = North_; kMin = nK+1-nAxisCell; kMax = nK; kOut = kMin-1
     elseif( NeiLBot(iBlock) == NOBLK) then
        iHemisphere = South_; kMin = 1; kMax = nAxisCell; kOut = kMax+1
     else
        CYCLE
     endif

     do i=1,nI
        r = r_Blk(i,1,1,iBlock)
        if(TypeGeometry == 'spherical_lnr') r = alog(r)
        if(TypeGeometry == 'spherical_genr') r = r_to_gen(r)
        iR = ceiling( (r - XyzMin_D(1))/MinDxValue + 0.1)

        ! Average small cells in the kMin:kMax ring(s)
        do k=kMin, kMax; 
           Volume = 1.0/vInv_CB(i,1,k,iBlock)
           Buffer_VIII(Volume_,iR,Geom_,iHemisphere) = &
                Buffer_VIII(Volume_,iR,Geom_,iHemisphere) + nJ*Volume
           do j=1,nJ
              Buffer_VIII(:,iR,Sum_,iHemisphere) = &
                   Buffer_VIII(:,iR,Sum_,iHemisphere) &
                   + Volume*State_VGB(:,i,j,k,iBlock)
           end do
        end do

        ! Calculate moments of the values of the kOut ring
        do j = 1, nJ
           x = x_BLK(i,j,kOut,iBlock)
           y = y_BLK(i,j,kOut,iBlock)
           State_V = State_VGB(:,i,j,kOut,iBlock)

           Buffer_VIII(:,iR,SumXLeft_,iHemisphere) = &
                Buffer_VIII(:,iR,SumXLeft_,iHemisphere)  + State_V*max(0.,-x)
           Buffer_VIII(:,iR,SumXRight_,iHemisphere) = &
                Buffer_VIII(:,iR,SumXRight_,iHemisphere) + State_V*max(0.,x)

           Buffer_VIII(:,iR,SumYLeft_,iHemisphere) = &
                Buffer_VIII(:,iR,SumYLeft_,iHemisphere)  + State_V*max(0.,-y)
           Buffer_VIII(:,iR,SumYRight_,iHemisphere) = &
                Buffer_VIII(:,iR,SumYRight_,iHemisphere) + State_V*max(0.,y)

           ! Add up abs(x) and x**2
           Buffer_VIII(SumX_,iR,Geom_,iHemisphere) = &
                Buffer_VIII(SumX_,iR,Geom_,iHemisphere) + abs(x)
           Buffer_VIII(SumX2_,iR,Geom_,iHemisphere) = &
                Buffer_VIII(SumX2_,iR,Geom_,iHemisphere) + x**2
        end do
     end do
  end do
     
  ! Collect all contributions around the axis
  call MPI_allreduce(Buffer_VIII, SumBuffer_VIII, nVar*nR*Geom_*2, MPI_REAL, &
       MPI_SUM, iComm, iError)

  ! Overwrite cells around the axis with a linear slope
  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

     if(rMin_BLK(iBlock) < r2FixAxis) then 
        nAxisCell = 2
     elseif(rMin_BLK(iBlock) < rFixAxis) then 
        nAxisCell = 1
     else
        CYCLE
     end if

     if( NeiLTop(iBlock) == NOBLK )then
        iHemisphere = North_; kMax = nK + 1; kMin = nK + 1 - nAxisCell
     elseif( NeiLBot(iBlock) == NOBLK) then
        iHemisphere = South_; kMin = 0; kMax = nAxisCell
     else
        CYCLE
     endif

     do i = 1, nI
        r = r_Blk(i,1,1,iBlock)
        if(TypeGeometry == 'spherical_lnr') r = alog(r)
        if(TypeGeometry == 'spherical_genr') r = r_to_gen(r)
        iR = ceiling( (r - XyzMin_D(1))/MinDxValue + 0.1)

        InvVolume= 1.0/SumBuffer_VIII(Volume_,iR,Geom_,iHemisphere)
        SumX     = 0.5*SumBuffer_VIII(SumX_  ,iR,Geom_,iHemisphere)
        InvSumX2 = 0.5/SumBuffer_VIII(SumX2_ ,iR,Geom_,iHemisphere)

        State_V = InvVolume*SumBuffer_VIII(:,iR,Sum_,iHemisphere)

        ! Limit the slope for each variable
        do iVar = 1, nVar
           SumXAvg = SumX*State_V(iVar)

           dLeft  = SumXAvg - SumBuffer_VIII(iVar,iR,SumXLeft_,iHemisphere)
           dRight = SumBuffer_VIII(iVar,iR,SumXRight_,iHemisphere) - SumXAvg
           dStateDx_V(iVar) = (sign(0.5,dLeft)+sign(0.5,dRight))*InvSumX2 &
               *min(1.5*abs(dLeft),1.5*abs(dRight),0.5*abs(dLeft+dRight))

           dLeft  = SumXAvg - SumBuffer_VIII(iVar,iR,SumYLeft_,iHemisphere)
           dRight = SumBuffer_VIII(iVar,iR,SumYRight_,iHemisphere) - SumXAvg
           dStateDy_V(iVar) = InvSumX2 * (sign(0.5,dLeft)+sign(0.5,dRight)) &
                *min(1.5*abs(dLeft),1.5*abs(dRight),0.5*abs(dLeft+dRight))
        end do

        ! Apply fit to each cell within the supercell
        do k = kMin, kMax; do j = 1, nJ

           State_VGB(:,i,j,k ,iBlock) = State_V &
                + dStateDx_V*x_BLK(i,j,k,iBlock) &
                + dStateDy_V*y_BLK(i,j,k,iBlock)

           call calc_energy_point(i,j,k ,iBlock)
        end do; end do
     end do

  end do

  deallocate(Buffer_VIII, SumBuffer_VIII)

  if(DoTestMe)then
     if(.not.UnusedBlk(BlkTest)) &
          write(*,*) NameSub,' final state, energy=', &
          State_VGB(:,iTest,jTest,kTest,BlkTest), &
          Energy_GBI(iTest,jTest,kTest,BlkTest,:)
  end if

end subroutine fix_axis_cells
!=============================================================================
subroutine fix_axis_cells_cyl

  use ModProcMH, ONLY: iComm
  use ModMain, ONLY: nJ, nK, nBlock, UnusedBlk
  use ModAdvance, ONLY: nVar, State_VGB, r2FixAxis
  use ModGeometry, ONLY: XyzMin_D, XyzMax_D, MinDxValue, &
       x_Blk, y_Blk, z_Blk, Dx_Blk, Dz_Blk, vInv_CB
  use ModEnergy, ONLY: calc_energy_point
  use ModParallel, ONLY: NeiLEast, NOBLK
  use ModMpi
  implicit none

  real, allocatable:: Buffer_VII(:,:,:), SumBuffer_VII(:,:,:)

  integer, parameter :: Sum_=1, SumXLeft_=2, SumXRight_=3, &
       SumYLeft_=4, SumYRight_=5, Geom_=6
  integer, parameter :: Volume_=1, SumX_=2, SumX2_=3
  integer :: i, j, k, iBlock, iZ, nZ, iError
  integer :: iVar, nAxisCell
  real :: MinDzValue
  real :: x, y, z, Volume, InvVolume, SumX, SumXAvg, InvSumX2, dLeft, dRight
  real :: State_V(nVar), dStateDx_V(nVar), dStateDy_V(nVar)

  !--------------------------------------------------------------------------

  ! Maximum number of cells along the axis
  MinDzValue = MinDxValue*dz_Blk(1)/dx_Blk(1)
  nZ = nint((XyzMax_D(3)-XyzMin_D(3))/MinDzValue)
  allocate(Buffer_VII(nVar, nZ, 1:Geom_), SumBuffer_VII(nVar, nZ, 1:Geom_))

  Buffer_VII    = 0.0
  SumBuffer_VII = 0.0

  if (r2FixAxis > 0.1) then 
     nAxisCell = 2
  else
     nAxisCell = 1
  end if

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. NeiLeast(iBlock) /= NOBLK) CYCLE

     do k=1,nK
        z = z_Blk(1,1,k,iBlock)
        iZ = ceiling( (z - XyzMin_D(3))/MinDzValue + 0.1)

        do i=1, nAxisCell
           Volume = 1.0/vInv_CB(i,1,k,iBlock)
           Buffer_VII(Volume_,iZ,Geom_) = Buffer_VII(Volume_,iZ,Geom_) &
                + nJ*Volume
           do j=1,nJ
              Buffer_VII(:,iZ,Sum_) = Buffer_VII(:,iZ,Sum_) &
                   + Volume*State_VGB(:,i,j,k,iBlock)
           end do
        end do
        do j=1,nJ
           x       = x_BLK(nAxisCell+1,j,k,iBlock)
           y       = y_BLK(nAxisCell+1,j,k,iBlock)
           State_V = State_VGB(:,nAxisCell+1,j,k,iBlock)

           Buffer_VII(:,iZ,SumXLeft_) = Buffer_VII(:,iZ,SumXLeft_) &
                + State_V*max(0.,-x)
           Buffer_VII(:,iZ,SumXRight_) = Buffer_VII(:,iZ,SumXRight_) &
                + State_V*max(0.,x)

           Buffer_VII(:,iZ,SumYLeft_) = Buffer_VII(:,iZ,SumYLeft_) &
                + State_V*max(0.,-y)
           Buffer_VII(:,iZ,SumYRight_) = Buffer_VII(:,iZ,SumYRight_) &
                + State_V*max(0.,y)

           ! Add up abs(x) and x**2
           Buffer_VII(SumX_,iZ,Geom_)  = Buffer_VII(SumX_,iZ,Geom_) + abs(x)
           Buffer_VII(SumX2_,iZ,Geom_) = Buffer_VII(SumX2_,iZ,Geom_) + x**2
        end do
     end do
  end do
     
  call MPI_allreduce(Buffer_VII, SumBuffer_VII, nVar*nZ*Geom_, MPI_REAL, &
       MPI_SUM, iComm, iError)

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. NeiLeast(iBlock) /= NOBLK) CYCLE

     do k=1,nK
        z = z_Blk(1,1,k,iBlock)
        iZ = ceiling( (z - XyzMin_D(3))/MinDzValue + 0.1)

        InvVolume = 1.0/SumBuffer_VII(Volume_,iZ,Geom_)
        SumX      = 0.5*SumBuffer_VII(SumX_  ,iZ,Geom_)
        InvSumX2  = 0.5/SumBuffer_VII(SumX2_ ,iZ,Geom_)

        State_V = InvVolume*SumBuffer_VII(:,iZ,Sum_)

        ! Limit the slope for each variable
        do iVar = 1, nVar
           SumXAvg = SumX*State_V(iVar)

           dLeft  = SumXAvg - SumBuffer_VII(iVar,iZ,SumXLeft_)
           dRight = SumBuffer_VII(iVar,iZ,SumXRight_) - SumXAvg
           dStateDx_V(iVar) = (sign(0.5,dLeft)+sign(0.5,dRight))*InvSumX2 &
               *min(1.5*abs(dLeft),1.5*abs(dRight),0.5*abs(dLeft+dRight))

           dLeft  = SumXAvg - SumBuffer_VII(iVar,iZ,SumYLeft_)
           dRight = SumBuffer_VII(iVar,iZ,SumYRight_) - SumXAvg
           dStateDy_V(iVar) = InvSumX2 * (sign(0.5,dLeft)+sign(0.5,dRight)) &
                *min(1.5*abs(dLeft),1.5*abs(dRight),0.5*abs(dLeft+dRight))
        end do

        ! Apply fit to each cell within the supercell
        do j=1, nJ; do i=1, nAxisCell

           State_VGB(:,i,j,k,iBlock) = State_V &
                + dStateDx_V*x_BLK(i,j,k,iBlock) &
                + dStateDy_V*y_BLK(i,j,k,iBlock)

           call calc_energy_point(i,j,k,iBlock)
        end do; end do
     end do

  end do

  deallocate(Buffer_VII, SumBuffer_VII)

end subroutine fix_axis_cells_cyl
