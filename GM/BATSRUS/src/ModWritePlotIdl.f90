!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_plot_idl(iFile, iBlock, nPlotVar, PlotVar, &
     xMin, xMax, yMin, yMax, zMin, zMax, DxBlock, DyBlock, DzBlock, nCell)

  ! Save all cells within plotting range, for each processor

  use ModProcMH
  use ModMain, ONLY: nI, nJ, nK, PROCtest, BLKtest, test_string, &
       x_, y_, z_, Phi_
  use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, Dx_BLK, Dy_BLK, Dz_BLK,&
       x1, x2, y1, y2, z1, z2, XyzStart_BLK, XyzMin_D, XyzMax_D
  use ModCovariant, ONLY: is_axial_geometry          
  use ModPhysics, ONLY : No2Io_V, UnitX_
  use ModIO
  use ModNumConst

  implicit none

  ! Arguments

  integer, intent(in)   :: iFile, iBlock
  integer, intent(in)   :: nPlotVar
  real,    intent(in)   :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar)
  real,    intent(in)   :: xMin,xMax,yMin,yMax,zMin,zMax
  real,    intent(inout):: DxBlock,DyBlock,DzBlock
  integer, intent(out)  :: nCell

  ! Local variables
  ! Indices and coordinates
  integer :: iVar, i, j, k, i2, j2, k2, iMin, iMax, jMin, jMax, kMin, kMax
  integer :: nRestrict, nRestrictX, nRestrictY, nRestrictZ
  real :: x,y,z,Dx,Restrict
  real :: xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
  real :: DxBlockOut
  real :: Plot_V(nPlotVarMax)

  real :: ySqueezed

  real, parameter:: cHalfMinusTiny=cHalf*(cOne-cTiny)

  character(len=*), parameter :: NameSub = 'write_plot_idl'
  logical :: DoTest, DoTestMe
  !---------------------------------------------------------------------------

  if(iProc==PROCtest .and. iBlock==BLKtest)then
     call set_oktest(NameSub, DoTest, DoTestMe)
  else
     DoTest=.false.; DoTestMe=.false.
  end if


  if(index(test_string,'SAVEPLOTALL')>0)then
     ! Save all cells of block including ghost cells

     nCell=0
     DxBlock=Dx_BLK(iBlock); DyBlock=Dy_BLK(iBlock); DzBlock=Dz_BLK(iBlock)
     DxBlockOut = DxBlock
     if (plot_dimensional(iFile)) DxBlockOut = DxBlock*No2Io_V(UnitX_)

     plot_range(1,iFile)=XyzMin_D(1)-2*DxBlock
     plot_range(2,iFile)=XyzMax_D(1)+2*DxBlock
     plot_range(3,iFile)=XyzMin_D(2)-2*DyBlock
     plot_range(4,iFile)=XyzMax_D(2)+2*DyBlock
     plot_range(5,iFile)=XyzMin_D(3)-2*DzBlock
     plot_range(6,iFile)=XyzMax_D(3)+2*DzBlock
     plot_Dx(1,iFile)=DxBlock
     plot_Dx(2,iFile)=DyBlock
     plot_Dx(3,iFile)=DzBlock

     do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
        nCell=nCell+1
        if (plot_dimensional(iFile)) then
           x = x_BLK(i,j,k,iBlock)*No2Io_V(UnitX_)
           y = y_BLK(i,j,k,iBlock)*No2Io_V(UnitX_)
           z = z_BLK(i,j,k,iBlock)*No2Io_V(UnitX_)
        else
           x = x_BLK(i,j,k,iBlock)
           y = y_BLK(i,j,k,iBlock)
           z = z_BLK(i,j,k,iBlock)
        end if
        if(save_binary)then
           write(unit_tmp)DxBlockOut,x,y,z,PlotVar(i,j,k,1:nPlotVar)
        else
           do iVar=1,nPlotVar
              Plot_V(iVar) = PlotVar(i,j,k,iVar)
              if(abs(Plot_V(iVar))<1.0e-99)Plot_V(iVar)=0.0
           end do
           write(unit_tmp,'(50(1pe13.5))') &
                DxBlockOut, x, y, z, Plot_V(1:nPlotVar)
        endif
     end do; end do; end do

     RETURN

  end if

  ! The range for the cell centers is Dx/2 wider
  xMin1 = xMin - cHalfMinusTiny*Dx_BLK(iBlock)
  xMax1 = xMax + cHalfMinusTiny*Dx_BLK(iBlock)
  yMin1 = yMin - cHalfMinusTiny*Dy_BLK(iBlock)
  yMax1 = yMax + cHalfMinusTiny*Dy_BLK(iBlock)
  zMin1 = zMin - cHalfMinusTiny*Dz_BLK(iBlock)
  zMax1 = zMax + cHalfMinusTiny*Dz_BLK(iBlock)

  nCell = 0
  if(is_axial_geometry())then                 
     ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
     ySqueezed = mod(xyzStart_BLK(Phi_,iBlock),cPi)
     ! Make sure that small angles are moved to Pi degrees for y=0 cut
     if(ySqueezed < 0.25*cPi .and. &
          abs(yMin+yMax-cTwoPi) < cTiny .and. yMax-yMin < 0.01) &
          ySqueezed = ySqueezed + cPi
  else                                          
     ySqueezed = xyzStart_BLK(y_,iBlock)
  end if                                         

  if(DoTestMe)then
     write(*,*) NameSub, 'xMin1,xMax1,yMin1,yMax1,zMin1,zMax1=',&
          xMin1,xMax1,yMin1,yMax1,zMin1,zMax1
     write(*,*) NameSub, 'xyzStart_BLK=',iBlock,xyzStart_BLK(:,iBlock)
     write(*,*) NameSub, 'ySqueezed =',ySqueezed
     write(*,*) NameSub, 'xyzEnd=',xyzStart_BLK(x_,iBlock)+(nI-1)*Dx_BLK(iBlock),&
          ySqueezed+(nJ-1)*Dy_BLK(iBlock),&
          xyzStart_BLK(z_,iBlock)+(nK-1)*Dz_BLK(iBlock)
  end if

  ! If block is fully outside of cut then cycle
  if(  xyzStart_BLK(x_,iBlock) > xMax1.or.&
       xyzStart_BLK(x_,iBlock)+(nI-1)*Dx_BLK(iBlock) < xMin1.or.&
       ySqueezed > yMax1.or.&
       ySqueezed+(nJ-1)*Dy_BLK(iBlock) < yMin1.or.&  
       xyzStart_BLK(z_,iBlock) > zMax1.or.&
       xyzStart_BLK(z_,iBlock)+(nK-1)*Dz_BLK(iBlock) < zMin1)&
       RETURN

  Dx = plot_Dx(1,iFile)
  DxBlock=Dx_BLK(iBlock); DyBlock=Dy_BLK(iBlock); DzBlock=Dz_BLK(iBlock)

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

  if(DxBlock>=Dx)then
     ! Cell is equal or coarser than Dx, save all cells in cut
     DxBlockOut = DxBlock
     if (plot_dimensional(iFile))DxBlockOut = DxBlockOut*No2Io_V(UnitX_)
     do k=kMin,kMax; do j=jMin,jMax; do i=iMin,iMax
        x = x_BLK(i,j,k,iBlock)
        y = y_BLK(i,j,k,iBlock)
        z = z_BLK(i,j,k,iBlock)

        if(x<x1 .or. x>x2 .or. y<y1 .or. y>y2 .or. z<z1 .or. z>z2) CYCLE

        if (plot_dimensional(iFile)) then
           x = x*No2Io_V(UnitX_)
           y = y*No2Io_V(UnitX_)
           z = z*No2Io_V(UnitX_)
        end if
        if(save_binary)then
           write(unit_tmp) DxBlockOut, x, y, z, PlotVar(i,j,k,1:nPlotVar)
        else
           do iVar=1,nPlotVar
              Plot_V(iVar) = PlotVar(i,j,k,iVar)
              if(abs(Plot_V(iVar)) < 1.0e-99) Plot_V(iVar) = 0.0
           end do
           write(unit_tmp,'(50(1pe13.5))') &
                DxBlockOut, x, y, z, Plot_V(1:nPlotVar)
        endif
        nCell = nCell+1
     end do; end do; end do
  else
     ! Block is finer then required resolution
     ! Calculate restriction factor
     nRestrict = min(nI,nint(Dx/DxBlock))

     ! Calclulate restricted cell size
     DxBlock    = nRestrict*DxBlock
     DyBlock    = nRestrict*DyBlock
     DzBlock    = nRestrict*DzBlock
     DxBlockOut = DxBlock
     if (plot_dimensional(iFile))DxBlockOut = DxBlockOut*No2Io_V(UnitX_)

     ! Restriction is limited by the width of the plotting region
     nRestrictX = min(iMax-iMin+1, nRestrict)
     nRestrictY = min(jMax-jMin+1, nRestrict)
     nRestrictZ = min(kMax-kMin+1, nRestrict)

     ! Factor for taking the average
     Restrict=1./(nRestrictX*nRestrictY*nRestrictZ)

     if(DoTestMe) write(*,*) NameSub, 'nRestrict, X, Y, Z,Restrict=',&
          nRestrict, nRestrictX, nRestrictY, nRestrictZ, Restrict

     ! Loop for the nRestrictX*nRestrictY*nRestrictZ bricks inside the cut
     do k=kMin,kMax,nRestrictZ
        k2=k+nRestrictZ-1
        do j=jMin,jMax,nRestrictY
           j2=j+nRestrictY-1
           do i=iMin,iMax,nRestrictX
              i2=i+nRestrictX-1
              x = 0.5*(x_BLK(i,j,k,iBlock)+x_BLK(i2,j2,k2,iBlock))
              y = 0.5*(y_BLK(i,j,k,iBlock)+y_BLK(i2,j2,k2,iBlock))
              z = 0.5*(z_BLK(i,j,k,iBlock)+z_BLK(i2,j2,k2,iBlock))

              if(x<x1 .or. x>x2 .or. y<y1 .or. y>y2 .or. z<z1 .or. z>z2) &
                   CYCLE

              if(plot_dimensional(iFile))then
                 x = x*No2Io_V(UnitX_)
                 y = y*No2Io_V(UnitX_)
                 z = z*No2Io_V(UnitX_)
              end if
              do iVar=1,nPlotVar
                 Plot_V(iVar) = Restrict*sum(PlotVar(i:i2,j:j2,k:k2,iVar))
              end do
              if(save_binary)then
                 write(unit_tmp)DxBlockOut,x,y,z,Plot_V(1:nPlotVar)
              else
                 do iVar=1,nPlotVar
                    if(abs(Plot_V(iVar)) < 1.0E-99)Plot_V(iVar)=0.0
                 end do
                 write(unit_tmp,'(50es13.5)')DxBlockOut,x,y,z,&
                      Plot_V(1:nPlotVar)
              endif
              nCell=nCell+1
           end do
        end do
     end do
  end if

end subroutine write_plot_idl
