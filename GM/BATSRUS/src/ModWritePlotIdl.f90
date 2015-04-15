!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
subroutine write_plot_idl(iFile, iBlock, nPlotVar, PlotVar, &
     xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock, nCell)

  ! Save all cells within plotting range, for each processor

  use ModProcMH
  use ModMain,     ONLY: PROCtest, BLKtest, test_string
  use ModGeometry, ONLY: x1, x2, y1, y2, z1, z2, XyzStart_BLK
  use ModIO,       ONLY: save_binary, plot_type1, plot_dx, plot_range, &
       nPlotVarMax
  use ModNumConst, ONLY: cPi, cTwoPi
  use ModKind,     ONLY: nByteReal
  use ModIoUnit,   ONLY: UnitTmp_
  use BATL_size,   ONLY: nGI, nGJ, nGK
  use BATL_lib,    ONLY: IsRLonLat, IsCylindrical, IsCartesianGrid, &
       CoordMin_D, CoordMax_D, CoordMin_DB, CellSize_DB, &
       nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK, &
       Xyz_DGB, x_, y_, z_, Phi_

  implicit none

  ! Arguments

  integer, intent(in)   :: iFile, iBlock
  integer, intent(in)   :: nPlotVar
  real,    intent(in)   :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
  real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
  real,    intent(inout):: DxBlock,DyBlock,DzBlock
  integer, intent(out)  :: nCell

  ! Local variables
  ! Indices and coordinates
  integer :: iVar, i, j, k, i2, j2, k2, iMin, iMax, jMin, jMax, kMin, kMax
  integer :: nRestrict, nRestrictX, nRestrictY, nRestrictZ
  real :: Coord_D(3), x, y, z, ySqueezed, Dx, Restrict
  real :: xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
  real :: Plot_V(nPlotVarMax)
  logical:: IsBinary, DoSaveGenCoord

  real:: cHalfMinusTiny

  character(len=*), parameter :: NameSub = 'write_plot_idl'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

  if(iProc==PROCtest .and. iBlock==BLKtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if

  if(nByteReal == 8)then
     cHalfMinusTiny = 0.5*(1.0 - 1e-9)
  else
     cHalfMinusTiny = 0.5*(1.0 - 1e-6)
  end if

  IsBinary = save_binary .and. plot_type1 /= 'cut_pic'

  ! Save generalized coordinates for cuts out of non-Cartesian grids
  DoSaveGenCoord = plot_type1(1:3) == 'cut' .and. .not. IsCartesianGrid

  if(index(test_string,'SAVEPLOTALL')>0)then

     ! Save all cells of block including ghost cells
     DxBlock = CellSize_DB(x_,iBlock)
     DyBlock = CellSize_DB(y_,iBlock) 
     DzBlock = CellSize_DB(z_,iBlock)

     plot_range(1,iFile) = CoordMin_D(1) - nGI*DxBlock
     plot_range(2,iFile) = CoordMax_D(1) + nGI*DxBlock
     plot_range(3,iFile) = CoordMin_D(2) - nGJ*DyBlock
     plot_range(4,iFile) = CoordMax_D(2) + nGJ*DyBlock
     plot_range(5,iFile) = CoordMin_D(3) - nGK*DzBlock
     plot_range(6,iFile) = CoordMax_D(3) + nGK*DzBlock
     plot_Dx(1,iFile) = DxBlock
     plot_Dx(2,iFile) = DyBlock
     plot_Dx(3,iFile) = DzBlock

     nCell=0
     do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
        nCell=nCell+1
        if(DoSaveGenCoord)then
           Coord_D = CoordMin_DB(:,iBlock) &
                + ((/i,j,k/)-0.5)*CellSize_DB(:,iBlock)
        else
           Coord_D = Xyz_DGB(:,i,j,k,iBlock)
        end if
        if(IsBinary)then
           write(UnitTmp_) DxBlock, Coord_D, PlotVar(i,j,k,1:nPlotVar)
        else
           do iVar=1,nPlotVar
              Plot_V(iVar) = PlotVar(i,j,k,iVar)
              if(abs(Plot_V(iVar))<1.0d-99)Plot_V(iVar)=0.0
           end do
           write(UnitTmp_,'(50(1pe13.5))') DxBlock, Coord_D, Plot_V(1:nPlotVar)
        endif
     end do; end do; end do

     RETURN

  end if

  ! The range for the cell centers is Dx/2 wider
  xMin1 = xMin - cHalfMinusTiny*CellSize_DB(x_,iBlock)
  xMax1 = xMax + cHalfMinusTiny*CellSize_DB(x_,iBlock)
  yMin1 = yMin - cHalfMinusTiny*CellSize_DB(y_,iBlock)
  yMax1 = yMax + cHalfMinusTiny*CellSize_DB(y_,iBlock)
  zMin1 = zMin - cHalfMinusTiny*CellSize_DB(z_,iBlock)
  zMax1 = zMax + cHalfMinusTiny*CellSize_DB(z_,iBlock)

  if(IsRLonLat .or. IsCylindrical .and. .not.DoSaveGenCoord)then
     ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
     ySqueezed = mod(xyzStart_BLK(Phi_,iBlock),cPi)
     ! Make sure that small angles are moved to Pi degrees for y=0 cut
     if(ySqueezed < 0.25*cPi .and. &
          abs(yMin+yMax-cTwoPi) < 1e-6 .and. yMax-yMin < 0.01) &
          ySqueezed = ySqueezed + cPi
  else
     ySqueezed = xyzStart_BLK(y_,iBlock)
  end if                                         

  if(DoTestMe)then
     write(*,*) NameSub, 'xMin1,xMax1,yMin1,yMax1,zMin1,zMax1=',&
          xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
     write(*,*) NameSub, 'xyzStart_BLK=',iBlock,xyzStart_BLK(:,iBlock)
     write(*,*) NameSub, 'ySqueezed =',ySqueezed
     write(*,*) NameSub, 'xyzEnd=', &
          xyzStart_BLK(x_,iBlock)+(nI-1)*CellSize_DB(x_,iBlock),&
          ySqueezed + (nJ-1)*CellSize_DB(y_,iBlock),&
          xyzStart_BLK(z_,iBlock)+(nK-1)*CellSize_DB(z_,iBlock)
  end if

  ! If block is fully outside of cut then cycle
  if(  xyzStart_BLK(x_,iBlock) > xMax1.or.&
       xyzStart_BLK(x_,iBlock)+(nI-1)*CellSize_DB(x_,iBlock) < xMin1.or.&
       ySqueezed > yMax1.or.&
       ySqueezed+(nJ-1)*CellSize_DB(y_,iBlock) < yMin1.or.&  
       xyzStart_BLK(z_,iBlock) > zMax1.or.&
       xyzStart_BLK(z_,iBlock)+(nK-1)*CellSize_DB(z_,iBlock) < zMin1)&
       RETURN

  Dx = plot_Dx(1,iFile)
  DxBlock = CellSize_DB(x_,iBlock)
  DyBlock = CellSize_DB(y_,iBlock)
  DzBlock = CellSize_DB(z_,iBlock)

  ! Calculate index limits of cells inside cut
  iMin = max(1 ,floor((xMin1-xyzStart_BLK(x_,iBlock))/DxBlock)+2)
  iMax = min(nI,floor((xMax1-xyzStart_BLK(x_,iBlock))/DxBlock)+1)

  jMin = max(1 ,floor((yMin1-ySqueezed)/DyBlock)+2)
  jMax = min(nJ,floor((yMax1-ySqueezed)/DyBlock)+1)

  kMin = max(1 ,floor((zMin1-xyzStart_BLK(z_,iBlock))/DzBlock)+2)
  kMax = min(nK,floor((zMax1-xyzStart_BLK(z_,iBlock))/DzBlock)+1)

  if(DoTestMe)then
     write(*,*) NameSub, 'iMin,iMax,jMin,jMax,kMin,kMax=',&
          iMin,iMax,jMin,jMax,kMin,kMax
     write(*,*) NameSub, 'DxBlock,x1,y1,z1',DxBlock,xyzStart_BLK(:,iBlock)
     write(*,*) NameSub, 'ySqueezed  =',ySqueezed
     write(*,*) NameSub, 'xMin1,xMax1=',xMin1,xMax1
     write(*,*) NameSub, 'yMin1,yMax1=',yMin1,yMax1
     write(*,*) NameSub, 'zMin1,zMax1=',zMin1,zMax1
  end if

  ! Count number of cells saved into the plot
  nCell = 0
  if(DxBlock >= Dx)then
     ! Cell is equal or coarser than Dx, save all cells in cut
     do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax
        x = Xyz_DGB(x_,i,j,k,iBlock)
        y = Xyz_DGB(y_,i,j,k,iBlock)
        z = Xyz_DGB(z_,i,j,k,iBlock)

        ! Check if we are inside the Cartesian box
        if(x<x1 .or. x>x2 .or. y<y1 .or. y>y2 .or. z<z1 .or. z>z2) CYCLE

        if(DoSaveGenCoord)then
           Coord_D = CoordMin_DB(:,iBlock) &
                + ((/i,j,k/)-0.5)*CellSize_DB(:,iBlock)
        else
           Coord_D = Xyz_DGB(:,i,j,k,iBlock)
        end if

        if(IsBinary)then
           write(UnitTmp_) DxBlock, Coord_D, PlotVar(i,j,k,1:nPlotVar)
        else
           do iVar=1, nPlotVar
              Plot_V(iVar) = PlotVar(i,j,k,iVar)
              if(abs(Plot_V(iVar)) < 1.0d-99) Plot_V(iVar) = 0.0
           end do
           write(UnitTmp_,'(50es13.5)') DxBlock, Coord_D, Plot_V(1:nPlotVar)
        endif
        nCell = nCell+1
     end do; end do; end do
  else
     ! Block is finer then required resolution
     ! Calculate restriction factor
     nRestrict = min(nI, nint(Dx/DxBlock))

     ! Calclulate restricted cell size
     DxBlock    = nRestrict*DxBlock
     DyBlock    = nRestrict*DyBlock
     DzBlock    = nRestrict*DzBlock

     ! Restriction is limited by the width of the plotting region
     nRestrictX = min(iMax-iMin+1, nRestrict)
     nRestrictY = min(jMax-jMin+1, nRestrict)
     nRestrictZ = min(kMax-kMin+1, nRestrict)

     ! Factor for taking the average
     Restrict = 1./(nRestrictX*nRestrictY*nRestrictZ)

     if(DoTestMe) write(*,*) NameSub, 'nRestrict, X, Y, Z,Restrict=',&
          nRestrict, nRestrictX, nRestrictY, nRestrictZ, Restrict

     ! Loop for the nRestrictX*nRestrictY*nRestrictZ bricks inside the cut
     do k = kMin, kMax, nRestrictZ
        k2 = k + nRestrictZ - 1
        do j = jMin, jMax, nRestrictY
           j2 = j + nRestrictY - 1
           do i = iMin, iMax, nRestrictX
              i2 = i + nRestrictX - 1
              x = 0.5*(Xyz_DGB(x_,i,j,k,iBlock) + Xyz_DGB(x_,i2,j2,k2,iBlock))
              y = 0.5*(Xyz_DGB(y_,i,j,k,iBlock) + Xyz_DGB(y_,i2,j2,k2,iBlock))
              z = 0.5*(Xyz_DGB(z_,i,j,k,iBlock) + Xyz_DGB(z_,i2,j2,k2,iBlock))

              if(x<x1 .or. x>x2 .or. y<y1 .or. y>y2 .or. z<z1 .or. z>z2) &
                   CYCLE

              if(DoSaveGenCoord)then
                 Coord_D = CoordMin_DB(:,iBlock) &
                      + (0.5*(/i+i2,j+j2,k+k2/) - 0.5)*CellSize_DB(:,iBlock)
              else
                 Coord_D = (/x, y, z/)
              end if

              do iVar=1,nPlotVar
                 Plot_V(iVar) = Restrict*sum(PlotVar(i:i2,j:j2,k:k2,iVar))
              end do
              if(IsBinary)then
                 write(UnitTmp_)DxBlock, Coord_D, Plot_V(1:nPlotVar)
              else
                 do iVar = 1, nPlotVar
                    if(abs(Plot_V(iVar)) < 1.0d-99)Plot_V(iVar)=0.0
                 end do
                 write(UnitTmp_,'(50es13.5)')DxBlock, Coord_D,&
                      Plot_V(1:nPlotVar)
              endif
              nCell=nCell+1
           end do
        end do
     end do
  end if

end subroutine write_plot_idl
