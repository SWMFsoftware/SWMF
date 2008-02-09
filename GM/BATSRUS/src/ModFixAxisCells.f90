subroutine fix_axis_cells

  use ModProcMH, ONLY: iComm
  use ModMain, ONLY: nI, nJ, nK, nBlock, UnusedBlk
  use ModAdvance, ONLY: nVar, State_VGB, rFixAxis, r2FixAxis
  use ModGeometry, ONLY: TypeGeometry, XyzMin_D, XyzMax_D, MinDxValue, &
       x_Blk, y_Blk, z_Blk, r_BLK, rMin_BLK, Dy_Blk, far_field_bcs_blk
  use ModConst, ONLY: cTwoPi
  use ModEnergy, ONLY: calc_energy_point
  use ModParallel, ONLY: NeiLBot, NeiLTop, NOBLK
  use ModMpi
  implicit none

  real, allocatable:: Buffer_VIII(:,:,:,:), SumBuffer_VIII(:,:,:,:)

  integer, parameter :: Sum_=1, SumXLeft_=2, SumXRight_=3, &
       SumYLeft_=4, SumYRight_=5, SumX_=6
  integer, parameter :: North_=1, South_=2
  integer :: i, j, k, kMin, kMax, kOut, iBlock, iHemisphere, iR, nR, iError
  integer :: iVar, nAxisCell
  real :: r, x, y, InvNCell, SumX, SumXAvg, InvSumX2, dLeft, dRight
  real :: State_V(nVar), dStateDx_V(nVar), dStateDy_V(nVar)

  !--------------------------------------------------------------------------

  ! Maximum number of cells in radial direction
  nR = nint((XyzMax_D(1)-XyzMin_D(1))/MinDxValue)
  if(.not.allocated(Buffer_VIII)) &
       allocate(Buffer_VIII(nVar, nR, 1:SumX_, North_:South_), &
                SumBuffer_VIII(nVar, nR, 1:SumX_, North_:South_))

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
        iR = ceiling( (r - XyzMin_D(1))/MinDxValue + 0.1)

        do k=kMin, kMax; do j=1,nJ
           Buffer_VIII(:,iR,Sum_,iHemisphere) = &
                Buffer_VIII(:,iR,Sum_,iHemisphere) &
                + State_VGB(:,i,j,k,iBlock)
        end do; end do
        do j=1,nJ
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

           ! Use first two variables of buffer to add up x and x**2
           Buffer_VIII(1,iR,SumX_,iHemisphere) = &
                Buffer_VIII(1,iR,SumX_,iHemisphere) + abs(x)
           Buffer_VIII(2,iR,SumX_,iHemisphere) = &
                Buffer_VIII(2,iR,SumX_,iHemisphere) + x**2
        end do
     end do
  end do
     
  call MPI_allreduce(Buffer_VIII, SumBuffer_VIII, nVar*nR*SumX_*2, MPI_REAL, &
       MPI_SUM, iComm, iError)

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

     if(rMin_BLK(iBlock) < r2FixAxis) then 
        nAxisCell = 2
     elseif(rMin_BLK(iBlock) < rFixAxis) then 
        nAxisCell = 1
     else
        CYCLE
     end if

     InvNCell = 1.0/(nAxisCell*nint(cTwoPi/Dy_BLK(iBlock)))

     if( NeiLTop(iBlock) == NOBLK )then
        iHemisphere = North_; kMax = nK + 1; kMin = nK + 1 - nAxisCell
     elseif( NeiLBot(iBlock) == NOBLK) then
        iHemisphere = South_; kMin = 0; kMax = nAxisCell
     else
        CYCLE
     endif

     do i=1,nI
        r = r_Blk(i,1,1,iBlock)
        if(TypeGeometry == 'spherical_lnr') r = alog(r)
        iR = ceiling( (r - XyzMin_D(1))/MinDxValue + 0.1)
        State_V = SumBuffer_VIII(:,iR,Sum_,iHemisphere)*InvNCell

        SumX     = 0.5*SumBuffer_VIII(1,iR,SumX_,iHemisphere)
        InvSumX2 = 0.5/max(SumBuffer_VIII(2,iR,SumX_,iHemisphere),1e-10)

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
        do k=kMin, kMax; do j=1,nJ

           State_VGB(:,i,j,k ,iBlock) = State_V &
                + dStateDx_V*x_BLK(i,j,k,iBlock) &
                + dStateDy_V*y_BLK(i,j,k,iBlock)

           call calc_energy_point(i,j,k ,iBlock)
        end do; end do
     end do

  end do

  deallocate(Buffer_VIII, SumBuffer_VIII)

end subroutine fix_axis_cells
