!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf

module ModFixAxisCells

  implicit none

  private ! except

  public:: fix_axis_cells

contains

  subroutine fix_axis_cells

    use ModProcMH, ONLY: iComm
    use ModMain, ONLY: nI, nJ, nK, nBlock, Unused_B, iTest, jTest, kTest, &
         BlkTest, x_, y_
    use ModAdvance, ONLY: nVar, State_VGB, Energy_GBI, rFixAxis, r2FixAxis
    use ModGeometry, ONLY: TypeGeometry, XyzMin_D, XyzMax_D, CellSize1Min, &
         r_BLK, rMin_BLK, far_field_bcs_blk
    use ModEnergy, ONLY: calc_energy_point
    use BATL_lib, ONLY: CoordMin_DB, CoordMax_DB, Lat_, &
         IsCylindrical, IsRLonLat, IsLogRadius, IsGenRadius, radius_to_gen, &
         Xyz_DGB, CellVolume_GB
    use ModNumConst, ONLY: cHalfPi
    use ModMpi

    real, allocatable:: Buffer_VIII(:,:,:,:)

    integer, parameter :: Sum_=1, SumXLeft_=2, SumXRight_=3, &
         SumYLeft_=4, SumYRight_=5, Geom_=6
    integer, parameter :: Volume_=1, SumX_=2, SumX2_=3
    integer :: i, j, k, kMin, kMax, kOut, iBlock, iHemisphere, iR, nR, iError
    integer :: iVar, nAxisCell
    real :: Beta
    real :: r, x, y, Volume, InvVolume, SumX, SumXAvg, InvSumX2, dLeft, dRight
    real :: State_V(nVar), dStateDx_V(nVar), dStateDy_V(nVar)

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'fix_axis_cells'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    if(DoTestMe)then
       if(.not.Unused_B(BlkTest)) &
            write(*,*) NameSub,' initial state, energy=', &
            State_VGB(:,iTest,jTest,kTest,BlkTest), &
            Energy_GBI(iTest,jTest,kTest,BlkTest,:)
    end if

    if(IsCylindrical)then
       call fix_axis_cells_cyl
       RETURN
    end if

    if(.not.IsRLonLat) &
         call stop_mpi(NameSub//': invalid geometry='//TypeGeometry)

    ! Maximum number of cells along the axis
    nR = nint((XyzMax_D(1)-XyzMin_D(1))/CellSize1Min)
    allocate(Buffer_VIII(nVar, nR, 1:Geom_, 1:2))

    Buffer_VIII    = 0.0

    do iBlock = 1, nBlock
       if(Unused_B(iBlock))CYCLE

       if(rMin_BLK(iBlock) < r2FixAxis) then 
          nAxisCell = 2
       elseif(rMin_BLK(iBlock) < rFixAxis) then 
          nAxisCell = 1
       else
          CYCLE
       end if

       ! Determine hemisphere
       if(    CoordMax_DB(Lat_,iBlock) > cHalfPi-1e-8)then
          iHemisphere = 1; kMin = nK+1-nAxisCell; kMax = nK; kOut = kMin-1
       elseif(CoordMin_DB(Lat_,iBlock) < -cHalfPi+1e-8)then
          iHemisphere = 2; kMin = 1; kMax = nAxisCell; kOut = kMax+1
       else
          CYCLE
       end if

       do i = 1, nI
          r = r_Blk(i,1,1,iBlock)
          if(IsLogRadius) r = alog(r)
          if(IsGenRadius) call radius_to_gen(r)
          iR = ceiling( (r - XyzMin_D(1))/CellSize1Min + 0.1)

          ! Average small cells in the kMin:kMax ring(s)
          do k=kMin, kMax; 
             Volume = CellVolume_GB(i,1,k,iBlock)
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
             x = Xyz_DGB(x_,i,j,kOut,iBlock)
             y = Xyz_DGB(y_,i,j,kOut,iBlock)
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
    call MPI_allreduce(MPI_IN_PLACE, Buffer_VIII, nVar*nR*Geom_*2, MPI_REAL, &
         MPI_SUM, iComm, iError)

    ! Overwrite cells around the axis with a linear slope
    do iBlock = 1, nBlock
       if(Unused_B(iBlock) .or. .not. far_field_BCs_BLK(iBlock)) CYCLE

       if(rMin_BLK(iBlock) < r2FixAxis) then 
          nAxisCell = 2
          Beta = 0.8
       elseif(rMin_BLK(iBlock) < rFixAxis) then 
          nAxisCell = 1
          Beta = 1.5
       else
          CYCLE
       end if

       if(    CoordMax_DB(Lat_,iBlock) > cHalfPi-1e-8)then
          iHemisphere = 1; kMax = nK + 1; kMin = nK + 1 - nAxisCell
       elseif(CoordMin_DB(Lat_,iBlock) < -cHalfPi+1e-8)then
          iHemisphere = 2; kMin = 0; kMax = nAxisCell
       else
          CYCLE
       end if

       do i = 1, nI
          r = r_Blk(i,1,1,iBlock)
          if(IsLogRadius) r = alog(r)
          if(IsGenRadius) call radius_to_gen(r)
          iR = ceiling( (r - XyzMin_D(1))/CellSize1Min + 0.1)

          InvVolume= 1.0/Buffer_VIII(Volume_,iR,Geom_,iHemisphere)
          SumX     = 0.5*Buffer_VIII(SumX_  ,iR,Geom_,iHemisphere)
          InvSumX2 = 2.0/Buffer_VIII(SumX2_ ,iR,Geom_,iHemisphere)

          State_V = InvVolume*Buffer_VIII(:,iR,Sum_,iHemisphere)

          ! Limit the slope for each variable
          do iVar = 1, nVar
             SumXAvg = SumX*State_V(iVar)

             dLeft  = SumXAvg - Buffer_VIII(iVar,iR,SumXLeft_,iHemisphere)
             dRight = Buffer_VIII(iVar,iR,SumXRight_,iHemisphere) - SumXAvg
             dStateDx_V(iVar) = (sign(0.5,dLeft)+sign(0.5,dRight))*InvSumX2 &
                  *min(Beta*abs(dLeft),Beta*abs(dRight),0.5*abs(dLeft+dRight))

             dLeft  = SumXAvg - Buffer_VIII(iVar,iR,SumYLeft_,iHemisphere)
             dRight = Buffer_VIII(iVar,iR,SumYRight_,iHemisphere) - SumXAvg
             dStateDy_V(iVar) = InvSumX2 * (sign(0.5,dLeft)+sign(0.5,dRight)) &
                  *min(Beta*abs(dLeft),Beta*abs(dRight),0.5*abs(dLeft+dRight))
          end do

          ! Apply fit to each cell within the supercell
          do k = kMin, kMax; do j = 1, nJ

             State_VGB(:,i,j,k ,iBlock) = State_V &
                  + dStateDx_V*Xyz_DGB(x_,i,j,k,iBlock) &
                  + dStateDy_V*Xyz_DGB(y_,i,j,k,iBlock)

             call calc_energy_point(i,j,k ,iBlock)
          end do; end do
       end do

    end do

    deallocate(Buffer_VIII)

    if(DoTestMe)then
       if(.not.Unused_B(BlkTest)) &
            write(*,*) NameSub,' final state, energy=', &
            State_VGB(:,iTest,jTest,kTest,BlkTest), &
            Energy_GBI(iTest,jTest,kTest,BlkTest,:)
    end if

  end subroutine fix_axis_cells
  !============================================================================
  subroutine fix_axis_cells_cyl

    use ModProcMH, ONLY: iComm
    use ModMain, ONLY: nJ, nK, nBlock, Unused_B, x_, y_, z_
    use ModAdvance, ONLY: nVar, State_VGB, r2FixAxis
    use ModGeometry, ONLY: XyzMin_D, XyzMax_D, CellSize1Min
    use ModEnergy, ONLY: calc_energy_point
    use ModMpi
    use BATL_lib, ONLY: Xyz_DGB, CellSize_DB, CellVolume_GB, CoordMin_DB, r_

    real, allocatable:: Buffer_VII(:,:,:)

    integer, parameter :: Sum_=1, SumXLeft_=2, SumXRight_=3, &
         SumYLeft_=4, SumYRight_=5, Geom_=6
    integer, parameter :: Volume_=1, SumX_=2, SumX2_=3
    integer :: i, j, k, iBlock, iZ, nZ, iError
    integer :: iVar, nAxisCell
    real :: MinDzValue, Beta
    real :: x, y, z, Volume, InvVolume, SumX, SumXAvg, InvSumX2, dLeft, dRight
    real :: State_V(nVar), dStateDx_V(nVar), dStateDy_V(nVar)

    character(len=*), parameter:: NameSub = 'fix_axis_cells_cyl'
    !--------------------------------------------------------------------------

    ! Maximum number of cells along the axis
    if(nK == 1)then
       nZ = 1
    else
       MinDzValue = CellSize1Min*CellSize_DB(z_,1)/CellSize_DB(x_,1)
       nZ = nint((XyzMax_D(3)-XyzMin_D(3))/MinDzValue)
    end if
    allocate(Buffer_VII(nVar,nZ,1:Geom_))

    Buffer_VII    = 0.0

    if (r2FixAxis > 0.1) then 
       nAxisCell = 2
       Beta = 0.8
    else
       nAxisCell = 1
       Beta = 1.5
    end if

    do iBlock = 1, nBlock
       if(Unused_B(iBlock)) CYCLE
       if(CoordMin_DB(r_,iBlock) > 0.0) CYCLE

       do k = 1, nK
          if(nK == 1)then
             iZ = 1
          else
             z = Xyz_DGB(z_,1,1,k,iBlock)
             iZ = ceiling( (z - XyzMin_D(3))/MinDzValue + 0.1)
          end if

          do i=1, nAxisCell
             Volume = CellVolume_GB(i,1,k,iBlock)
             Buffer_VII(Volume_,iZ,Geom_) = Buffer_VII(Volume_,iZ,Geom_) &
                  + nJ*Volume
             do j=1,nJ
                Buffer_VII(:,iZ,Sum_) = Buffer_VII(:,iZ,Sum_) &
                     + Volume*State_VGB(:,i,j,k,iBlock)
             end do
          end do
          do j=1,nJ
             x       = Xyz_DGB(x_,nAxisCell+1,j,k,iBlock)
             y       = Xyz_DGB(y_,nAxisCell+1,j,k,iBlock)
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

    call MPI_allreduce(MPI_IN_PLACE, Buffer_VII, nVar*nZ*Geom_, MPI_REAL, &
         MPI_SUM, iComm, iError)

    do iBlock = 1, nBlock

       if(Unused_B(iBlock)) CYCLE
       if(CoordMin_DB(r_,iBlock) > 0.0) CYCLE

       do k = 1, nK
          if(nK == 1)then
             iZ = 1
          else
             z = Xyz_DGB(z_,1,1,k,iBlock)
             iZ = ceiling( (z - XyzMin_D(3))/MinDzValue + 0.1)
          end if

          InvVolume = 1.0/Buffer_VII(Volume_,iZ,Geom_)
          SumX      = 0.5*Buffer_VII(SumX_  ,iZ,Geom_)
          InvSumX2  = 2.0/Buffer_VII(SumX2_ ,iZ,Geom_)

          State_V = InvVolume*Buffer_VII(:,iZ,Sum_)

          ! Limit the slope for each variable
          do iVar = 1, nVar
             SumXAvg = SumX*State_V(iVar)

             dLeft  = SumXAvg - Buffer_VII(iVar,iZ,SumXLeft_)
             dRight = Buffer_VII(iVar,iZ,SumXRight_) - SumXAvg
             dStateDx_V(iVar) = (sign(0.5,dLeft)+sign(0.5,dRight))*InvSumX2 &
                  *min(Beta*abs(dLeft),Beta*abs(dRight),0.5*abs(dLeft+dRight))

             !if(iVar==1)then
             !   write(*,*)'!!! iBlock=',iBlock
             !   write(*,*)'!!! SumX, InvSumX2 =', SumX, InvSumX2
             !   write(*,*)'!!! Avg, SumXAvg   =', State_V(iVar), SumXAvg 
             !   write(*,*)'!!! SumXLeft, Right=', &
             !        Buffer_VII(iVar,iZ,SumXLeft_), &
             !        Buffer_VII(iVar,iZ,SumXRight_)
             !   write(*,*)'!!! dLeft, dRight  =',  dLeft, dRight
             !   write(*,*)'!!! dStateDx_V     =', dStateDx_V(iVar)
             !end if

             dLeft  = SumXAvg - Buffer_VII(iVar,iZ,SumYLeft_)
             dRight = Buffer_VII(iVar,iZ,SumYRight_) - SumXAvg
             dStateDy_V(iVar) = InvSumX2 * (sign(0.5,dLeft)+sign(0.5,dRight)) &
                  *min(Beta*abs(dLeft),Beta*abs(dRight),0.5*abs(dLeft+dRight))

          end do

          !write(*,*)'!!! x     =', Xyz_DGB(x_,nAxisCell,1:nJ,1,iBlock)
          !write(*,*)'!!! Before=',State_VGB(1,nAxisCell,1:nJ,1,iBlock)


          ! Apply fit to each cell within the supercell
          do j=1, nJ; do i=1, nAxisCell

             State_VGB(:,i,j,k,iBlock) = State_V &
                  + dStateDx_V*Xyz_DGB(x_,i,j,k,iBlock) &
                  + dStateDy_V*Xyz_DGB(y_,i,j,k,iBlock)

             call calc_energy_point(i,j,k,iBlock)
          end do; end do

          !write(*,*)'!!! After =',State_VGB(1,nAxisCell,1:nJ,1,iBlock)

       end do

    end do

    deallocate(Buffer_VII)

  end subroutine fix_axis_cells_cyl

end module ModFixAxisCells
