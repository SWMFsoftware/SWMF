!^CFG FILE COVARIANT

subroutine fix_axis_cells

  use ModProcMH, ONLY: iComm
  use ModMain, ONLY: nI, nJ, nK, nBlock, UnusedBlk
  use ModAdvance, ONLY: nVar, State_VGB
  use ModGeometry, ONLY: TypeGeometry, XyzMin_D, XyzMax_D, MinDxValue, &
       x_Blk, y_Blk, z_Blk, r_BLK, Dy_Blk, far_field_bcs_blk
  use ModConst, ONLY: cTwoPi
  use ModEnergy, ONLY: calc_energy_point
  use ModParallel, ONLY: NeiLBot, NeiLTop, NOBLK
  use ModMpi
  implicit none

  real, allocatable:: Buffer_VIII(:,:,:,:), SumBuffer_VIII(:,:,:,:)

  integer, parameter :: Sum_=1, SumX_=2, SumY_=3, SumXX_=4
  integer, parameter :: North_=1, South_=2
  integer :: i, j, k, k1, iBlock, iHemisphere, iR, nR, nPhi, iError
  real :: r, StateAvg_V(nVar), dState_V(nVar)

  !--------------------------------------------------------------------------

  ! Maximum number of cells in radial direction
  nR = ceiling((XyzMax_D(1)-XyzMin_D(1))/MinDxValue)
  if(.not.allocated(Buffer_VIII)) &
       allocate(Buffer_VIII(nVar, nR, Sum_:SumXX_, North_:South_), &
                SumBuffer_VIII(nVar, nR, Sum_:SumXX_, North_:South_))

  Buffer_VIII    = 0.0
  SumBuffer_VIII = 0.0

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

     ! Determine hemisphere
     if( NeiLTop(iBlock) == NOBLK )then
        iHemisphere = North_
        k = nK
     elseif( NeiLBot(iBlock) == NOBLK) then
        iHemisphere = South_
        k = 1
     else
        CYCLE
     endif
     
     do i=1,nI
        r = r_Blk(i,1,k,iBlock)
        if(TypeGeometry == 'spherical_lnr') r = alog(r)
        iR = ceiling(r/MinDxValue+0.1)
        do j=1,nJ
           Buffer_VIII(:,iR,Sum_,iHemisphere) = &
                Buffer_VIII(:,iR,Sum_,iHemisphere) &
                + State_VGB(:,i,j,k,iBlock)
           Buffer_VIII(:,iR,SumX_,iHemisphere) = &
                Buffer_VIII(:,iR,SumX_,iHemisphere) &
                + State_VGB(:,i,j,k,iBlock)*x_BLK(i,j,k,iBlock)
           Buffer_VIII(:,iR,SumY_,iHemisphere) = &
                Buffer_VIII(:,iR,SumY_,iHemisphere) &
                + State_VGB(:,i,j,k,iBlock)*y_BLK(i,j,k,iBlock)
           Buffer_VIII(:,iR,SumXX_,iHemisphere) = &
                Buffer_VIII(:,iR,SumXX_,iHemisphere) &
                + x_BLK(i,j,k,iBlock)**2
        end do
     end do
  end do
     
  call MPI_allreduce(Buffer_VIII, SumBuffer_VIII, nVar*nR*8, MPI_REAL, &
       MPI_SUM, iComm, iError)

  do iBlock = 1, nBlock
     if(unusedBlk(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

     nPhi = nint(cTwoPi/Dy_BLK(iBlock))

     if( NeiLTop(iBlock) == NOBLK )then
        iHemisphere = North_
        k = nK
        k1= nK+1
     elseif( NeiLBot(iBlock) == NOBLK) then
        iHemisphere = South_
        k = 1
        k1= 0
     else
        CYCLE
     endif

     do i=1,nI
        r = r_Blk(i,1,k,iBlock)
        if(TypeGeometry == 'spherical_lnr') r = alog(r)
        iR = ceiling(r/MinDxValue+0.1)
        do j=1,nJ
           StateAvg_V = SumBuffer_VIII(:,iR,Sum_,iHemisphere)/nPhi
           dState_V = &
                ( SumBuffer_VIII(:,iR,SumX_,iHemisphere)*x_BLK(i,j,k,iBlock) &
                + SumBuffer_VIII(:,iR,SumY_,iHemisphere)*y_BLK(i,j,k,iBlock) &
                )/SumBuffer_VIII(:,iR,SumXX_,iHemisphere)

           State_VGB(:,i,j,k ,iBlock) = StateAvg_V + dState_V
           State_VGB(:,i,j,k1,iBlock) = StateAvg_V - dState_V

           call calc_energy_point(i,j,k ,iBlock)
           call calc_energy_point(i,j,k1,iBlock)
        end do
     end do

  end do

  deallocate(Buffer_VIII, SumBuffer_VIII)

end subroutine fix_axis_cells
