!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
!=============================================================================
module ModPlotShell

  use ModIO
  use ModNumConst, ONLY: cRadtoDeg, cDegToRad, cPi, cZero, cOne

  implicit none
  save

  private   ! except
  public:: write_plot_sph
  public:: set_plot_shell
  public:: write_plot_shell
  
  ! Size of current plot:
  integer :: nRad, nLon, nLat

  ! Ranges for current plot:
  real :: dR, dLat, dLon
  real :: radMin, radMax, lonMin, lonMax, latMin, latMax

  ! Local results container:
  ! Array of values written to file:
    real, allocatable :: PlotVarLocal_VIII(:,:,:,:)

contains
  !=============================================================================
  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

  !=============================================================================
  subroutine set_plot_shell(iFile, iBLK, nPlotvar, Plotvar)
    ! Save all cells within plotting range, for each processor, for a GeoSphere
    ! plot type.  iFile is the identifying integer for the file and is used to
    ! obtain the plot range and resolution.  iBlk is the number of the current
    ! block from which values are obtained.  nPlotVar PlotVar give the number
    ! of variables to save and an IxJxKxnPlotVar array of those values within 
    ! the block. 
    
    use ModMain,        ONLY: time_simulation, TypeCoordSystem
    use CON_axes,       ONLY: transform_matrix
    use ModGeometry,    ONLY: nI,nJ,nK, Xyz_DGB, CellSize_DB 
    use ModMain,        ONLY: BlkTest
    use ModInterpolate, ONLY: trilinear
    implicit none

    ! Arguments
    integer, intent(in) :: iFile, iBLK
    integer, intent(in) :: nPlotvar
    real,    intent(in) :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
    
    ! Local variables
    integer :: i, j, k, iVar
    
    real :: rNow, lonNow, latNow
    real :: x,y,z,xNorm, yNorm, zNorm, dxblk
    real :: xBlkMin,xBlkMax,yBlkMin,yBlkMax,zBlkMin,zBlkMax, minRblk,maxRblk
    real :: XyzGeo_D(3), XyzGm_D(3), GmToGeo_DD(3,3)
    
    character(len=200) :: NameVar, StringHeader
    
    character(len=*), parameter :: NameSub='set_plot_shell'
    logical :: DoTest, DoTestMe, oktest, oktest_me
    !---------------------------------------------------------------------------
    ! Check testing for subroutine:
    call CON_set_do_test(NameSub, DoTest, DoTestMe)
    
    ! Check testing for blocks:
    if(iBLK == BlkTest)then
       call set_oktest(NameSub,oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if
    
    ! Get plot area info from ModIO arrays:
    dR     = abs(plot_dx(1, iFile))
    dLon   = plot_dx(2, iFile) * cDegtoRad
    dLat   = plot_dx(3, iFile) * cDegtoRad
    radMin = plot_range(1, iFile)
    radMax = plot_range(2, iFile)
    lonMin = plot_range(3, iFile) * cDegtoRad
    lonMax = plot_range(4, iFile) * cDegtoRad
    latMin = plot_range(5, iFile) * cDegtoRad
    latMax = plot_range(6, iFile) * cDegtoRad
    
    if (DoTestMe) then
       write(*,*) NameSub//' Called for iBlk=      ', iBLK
       write(*,*) NameSub//' Raw plot_dx=          ', plot_dx(:,iFile)
       write(*,*) NameSub//' Raw plot_range=       ', plot_range(:,iFile)
       write(*,*) NameSub//' dR, dLon, dLat =      ', dR, dLon, dLat
       write(*,*) NameSub//' rad, Lon, Lat range = ',  &
            radMin, radMax, lonMin,lonMax,latMin,latMax
    end if
    
    ! Ensure good choices of dR, dLon, and dLat.
    if(mod(radMax-radMin, dR)  >1E-5) &
         call CON_stop(NameSub//': radius range not evenly divisible by dR')
    if(mod(lonMax-lonMin, dLon)>1E-5) &
         call CON_stop(NameSub//': Longitude range not evenly divisible by dLon')
    if(mod(latMax-latMin, dLat)>1E-5) &
         call CON_stop(NameSub//': Latitude range not evenly divisible by dLat')
    
    ! Set number of points:
    nRad = (radMax-radMin) / dR   + 1
    nLon = (lonMax-lonMin) / dLon + 1
    nLat = (latMax-latMin) / dLat + 1
    
    if(DoTestMe) write(*,*) NameSub//' nRad, nLon, nLat = ', nRad, nLon, nLat
    
    ! Allocate local & global results array:
    if(.not.allocated(PlotVarLocal_VIII)) then
       allocate(PlotVarLocal_VIII(nPlotVar, nRad, nLon, nLat))
       PlotVarLocal_VIII = 0.0
    end if
    
    ! Get coordinate transformation matrix:
    GmToGeo_DD = transform_matrix(Time_Simulation, TypeCoordSystem, 'GEO')
    
    ! get the max and min radial distance for this block so that we can check
    ! whether or not this block can to contibute to the plot.
    xBlkMin = 0.50*(Xyz_DGB(x_, 0, 0, 0,iBLK)+Xyz_DGB(x_,   1,   1  , 1,iBLK))
    xBlkMax = 0.50*(Xyz_DGB(x_,nI,nJ,nK,iBLK)+Xyz_DGB(x_,nI+1,nJ+1,nK+1,iBLK))
    yBlkMin = 0.50*(Xyz_DGB(y_, 0, 0, 0,iBLK)+Xyz_DGB(y_,   1,   1,   1,iBLK))
    yBlkMax = 0.50*(Xyz_DGB(y_,nI,nJ,nK,iBLK)+Xyz_DGB(y_,nI+1,nJ+1,nK+1,iBLK))
    zBlkMin = 0.50*(Xyz_DGB(z_, 0, 0, 0,iBLK)+Xyz_DGB(z_,   1,   1,   1,iBLK))
    zBlkMax = 0.50*(Xyz_DGB(z_,nI,nJ,nK,iBLK)+Xyz_DGB(z_,nI+1,nJ+1,nK+1,iBLK))

    if(DoTestMe) then
       write(*,*)NameSub//' limits of current block:'
       write(*,*)'      x:', xBlkMin, xBlkMax
       write(*,*)'      y:', xBlkMin, xBlkMax
       write(*,*)'      z:', xBlkMin, xBlkMax
       write(*,*)'cell sizes = ', CellSize_DB(x_:z_, iBlk)
    end if

    do i=1, nRad
       ! Set current radius:
       rNow = radMin+(i-1)*dR
       
       minRblk = sqrt(minmod(xBlkMin,xBlkMax)**2 + &
            minmod(yBlkMin,yBlkMax)**2 + &
            minmod(zBlkMin,zBlkMax)**2)
       maxRblk = sqrt((max(abs(xBlkMin),abs(xBlkMax)))**2 + &
            (max(abs(yBlkMin),abs(yBlkMax)))**2 + &
            (max(abs(zBlkMin),abs(zBlkMax)))**2)
       
       ! If the block does not intersect the sphere that we are plotting, cycle.
       ! This reduces unnecessary calculations.
       if(minRblk>rNow .or. maxRblk<rNow) cycle
       
       do j=1,nlon
          ! Current lat value:
          lonNow = (j-1)*dlon
          
          do k=1,nlat
             ! Current lon value:
             latNow = (k-1)*dlat
             
             if(DoTestMe) write(*,*)NameSub//': i,j,k = ',i,j,k
             
             ! get the cartesian coordinate from the spherical coordinates
             XyzGeo_D = (/&
                  rNow*cos(latNow)*cos(lonNow), &
                  rNow*cos(latNow)*sin(lonNow), &
                  rNow*sin(latNow)/)
             
             ! Convert from geographic coordinates to BATSRUS grid:
             XyzGm_D = matmul(GmToGeo_DD, XyzGeo_D)
             x = XyzGm_D(1)
             y = XyzGm_D(2)
             z = XyzGm_D(3)
             
             ! If point outside of this block, cycle.
             if ( x<xBlkMin .or. x>=xBlkMax .or. &
                  y<yBlkMin .or. y>=yBlkMax .or. &
                  z<zBlkMin .or. z>=zBlkMax) cycle
             
             if(DoTestMe) then
                write(*,*) 'rlatlon = ', rNow, lonNow, latNow
                write(*,*) 'geoxyz  = ', XyzGeo_D
                write(*,*) 'BATSxyz = ', x, y, z
             end if
             
             ! Get point in normalized coordinates, for speed:
             xNorm = (x-Xyz_DGB(x_,1,1,1,iBLK)) / CellSize_DB(x_,iBLK) +1.
             yNorm = (y-Xyz_DGB(y_,1,1,1,iBLK)) / CellSize_DB(y_,iBLK) +1.
             zNorm = (z-Xyz_DGB(z_,1,1,1,iBLK)) / CellSize_DB(z_,iBLK) +1.
             
             if(DoTestMe) write(*,*) 'NORM VARIABLES: ', xNorm, yNorm, zNorm

             !compute the interpolated values at the current location.
             do iVar=1, nPlotVar
                ! Interpolate up to ghost cells.
                PlotVarLocal_VIII(iVar, i,j,k) = trilinear(PlotVar(:,:,:,iVar),&
                     0,nI+1, 0,nJ+1, 0,nK+1, (/xNorm, yNorm, zNorm/))
             end do
             
          end do !lon loop
       end do    !lat loop
    end do       !rad loop

  end subroutine set_plot_shell

  !=============================================================================
  subroutine write_plot_shell(iFileIn, nPlotVarIn, NameVarIn_V, &
       NameUnitIn, NameFileIn)
    ! Collect results from all blocks and write to single output file.
    use ModMpi
    use ModMain,     ONLY: time_simulation, n_step
    use ModProcMH,   ONLY: iProc, nProc, iComm
    use ModPlotFile, ONLY: save_plot_file

    implicit none

    integer,          intent(in) :: iFileIn, nPlotvarIn
    character(len=*), intent(in) :: NameFileIn,NameUnitIn,NameVarIn_V(nPlotVarIn)

    ! Array of values written to file:
    real, allocatable :: PlotVar_VIII(:,:,:,:)

    integer :: iVar, iError
    character(len=200) :: NameVar

    logical :: DoTest,DoTestMe
    character(len=*), parameter :: NameSub='write_plot_shell'
    !---------------------------------------------------------------------------
    ! This subroutine does not support HDF output.
    if(plot_form(iFileIn) == 'hdf') call CON_stop(NameSub// &
         ': HDF file type not supported for Geo Sphere output.')
    
    ! Allocate final result array:
    if(iProc==0) allocate(PlotVar_VIII(nPlotVarIn, nRad, nLon, nLat))

    ! Collapse results on to head node:
    if(nProc>1)then
       call MPI_reduce(PlotVarLocal_VIII, PlotVar_VIII, &
            nPlotVarIn * nRad * nLon * nLat,            &
            MPI_REAL, MPI_SUM, 0, iComm, iError)
    else
       PlotVar_VIII = PlotVarLocal_VIII
    end if
    
    ! Save results to disk:
    if(iProc==0) then
       ! Build a single-line list of variable names.
       NameVar = 'r lon lat'
       do iVar=1, nPlotVarIn
          NameVar = trim(NameVar)  // ' ' // trim(NameVarIn_V(iVar))
       end do

       ! Call save_plot_file to write data to disk.
       call save_plot_file(NameFileIn, &
            TypeFileIn=TypeFile_I(iFileIn), &
            StringHeaderIn=NameUnitIn, &
            nStepIn=n_step, &
            TimeIn=time_simulation, &
            NameVarIn = NameVar, &
            CoordMinIn_D = (/radMin, cRadtoDeg*lonMin, cRadtoDeg*latMin/), &
            CoordMaxIn_D = (/radMax, cRadtoDeg*lonMax, cRadtoDeg*latMax/), &
            VarIn_VIII = PlotVar_VIII)
    end if
    
    ! Deallocate results arrays:.
    deallocate(PlotVarLocal_VIII)
    if(iProc==0) deallocate(PlotVar_VIII)
    
  end subroutine write_plot_shell

  !=============================================================================
  subroutine write_plot_sph(iFile,iBLK,nPlotvar,Plotvar, &
       nTheta,nPhi,rPlot,plotvarnames,H5Index,nBlkCellsN,nBlkCellsS)
    ! CANDIDATE FOR REMOVAL.  SHELL PLOTTING TOOLS PREFERRED.
    ! Save all cells within plotting range, for each processor
    
    use ModPhysics, ONLY : No2Io_V, UnitX_
    use ModGeometry, ONLY: nI,nJ,nK, Xyz_DGB, CellSize_DB 
    use ModMain,     ONLY: BlkTest
    use ModHdf5, ONLY: write_sph_var_hdf5
    implicit none
    
    integer, parameter:: lNameVar = 10
    ! Arguments
    integer, intent(in)  :: iFile, iBLK,nTheta,nPhi
    integer, intent(in)  :: nPlotvar
    integer, intent(inout) :: H5Index
    integer, intent(out) :: nBlkCellsN,nBlkCellsS
    real, intent(in)     :: PlotVar(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nplotvar)
    real, intent(in)     :: rPlot
    character(len=lNameVar), intent(in):: plotVarNames(nPlotVar)
    ! Local variables
    ! Indices and coordinates
    integer :: i,j,iVar
    integer :: i1,j1,k1,i2,j2,k2
    
    real :: x,y,z,xx,yy,zz
    real :: dx1,dy1,dz1,dx2,dy2,dz2
    real :: xx1,xx2,yy1,yy2,zz1,zz2,minRblk, maxRblk
    real :: rplot_out,theta_plot, phi_plot, dphi_plot, dtheta_plot, dxblk
    real :: theta_out, phi_out
    real :: PointVar(nPlotvarMax)

    logical :: oktest,oktest_me
    !---------------------------------------------------------------------------
    
    if(iBLK == BlkTest)then
       call set_oktest('write_plot_sph',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    end if
    
    nBLKcellsN = 0
    nBLKcellsS = 0
    
    ! get the max and min radial distance for this block so that we can check
    ! whether or not this block contibutes to the plot.
    xx1 = 0.50*(Xyz_DGB(x_, 0, 0, 0,iBLK)+Xyz_DGB(x_,   1,   1  , 1,iBLK))
    xx2 = 0.50*(Xyz_DGB(x_,nI,nJ,nK,iBLK)+Xyz_DGB(x_,nI+1,nJ+1,nK+1,iBLK))
    yy1 = 0.50*(Xyz_DGB(y_, 0, 0, 0,iBLK)+Xyz_DGB(y_,   1,   1,   1,iBLK))
    yy2 = 0.50*(Xyz_DGB(y_,nI,nJ,nK,iBLK)+Xyz_DGB(y_,nI+1,nJ+1,nK+1,iBLK))
    zz1 = 0.50*(Xyz_DGB(z_, 0, 0, 0,iBLK)+Xyz_DGB(z_,   1,   1,   1,iBLK))
    zz2 = 0.50*(Xyz_DGB(z_,nI,nJ,nK,iBLK)+Xyz_DGB(z_,nI+1,nJ+1,nK+1,iBLK))
    minRblk = sqrt(&
         minmod(xx1,xx2)**2 + minmod(yy1,yy2)**2 + minmod(zz1,zz2)**2)
    maxRblk = sqrt((max(abs(xx1),abs(xx2)))**2 + &
         (max(abs(yy1),abs(yy2)))**2 + &
         (max(abs(zz1),abs(zz2)))**2)
    
    ! If the block does not intersect the sphere that we are plotting then
    ! cycle
    if(minRblk>rplot .or. maxRblk<rplot)  return
    if (plot_form(ifile) == "hdf") then
       dtheta_plot = cPi/real(ntheta)
    else
       dtheta_plot = cPi/real(ntheta-1)
    end if
    dphi_plot   = 2.0*cPi/real(nphi)
    
    do i=1,ntheta
       if (plot_form(ifile) == "hdf") then !batl hdf plots use cell centered values
          theta_plot = (i-.5)*dtheta_plot 
       else
          theta_plot = (i-1)*dtheta_plot
       end if
       do j=1,nphi
          if (plot_form(ifile) == 'hdf') then
             phi_plot = (j-.5)*dphi_plot
          else
             phi_plot = (j-1)*dphi_plot
          end if
          
          
          !!SAVE  The following code was originally used to write out the spherical plot
          !!SAVE  files in a tilted frame which is the same as the ionosphere frame.
          !!SAVE  In order to allow spherical and x,y,z cuts and 3d files to exist in
          !!SAVE  the same coordiate system, this is changed so that the spherical plots
          !!SAVE  are written in the BATSRUS cartesian system.
          !!SAVE
          !!SAVE        ! get the coordinates in the tilted frame (same as the ionosphere)
          !!SAVE        x = rplot*sin(theta_plot)*cos(phi_plot)
          !!SAVE        y = rplot*sin(theta_plot)*sin(phi_plot)
          !!SAVE        z = rplot*cos(theta_plot)
          !!SAVE        xyz(1) = x; xyz(2) = y; xyz(3) = z
          !!SAVE        ! rotate these into the correct magnetospheric coordinates
          !!SAVE        call Rotate_into_Magneto(xyz, ThetaTilt)
          !!SAVE        x = xyz(1); y = xyz(2); z = xyz(3)
          
          ! get the cartesian coordinate from the spherical coordinates
          x = rplot*sin(theta_plot)*cos(phi_plot)
          y = rplot*sin(theta_plot)*sin(phi_plot)
          z = rplot*cos(theta_plot)
          
          ! check to see if this point is inside the block - if so print it out
          if (x >= xx1 .and. x < xx2 .and. &
               y >= yy1 .and. y < yy2 .and. &
               z >= zz1 .and. z < zz2 ) then
             
             !compute the interpolated values at the current location.
             
             PointVar=0.0
             
             ! Convert to normalized coordinates (index and position are the same)
             xx=(x-Xyz_DGB(x_,1,1,1,iBLK))/CellSize_DB(x_,iBLK)+1.
             yy=(y-Xyz_DGB(y_,1,1,1,iBLK))/CellSize_DB(y_,iBLK)+1.
             zz=(z-Xyz_DGB(z_,1,1,1,iBLK))/CellSize_DB(z_,iBLK)+1.
             
             ! Determine cell indices corresponding to location qx
             i1=floor(xx); i2=i1+1
             j1=floor(yy); j2=j1+1
             k1=floor(zz); k2=k1+1
             
             ! Distance relative to the cell centers
             dx1=xx-i1; dx2=1.-dx1
             dy1=yy-j1; dy2=1.-dy1
             dz1=zz-k1; dz2=1.-dz1
             
             ! Bilinear interpolation in 3D
             do iVar=1,nPlotVar
                PointVar(iVar) = &
                     dx1*(   dy1*(   dz1*PlotVar(i2,j2,k2,iVar)+&
                     dz2*PlotVar(i2,j2,k1,iVar))+&
                     dy2*(   dz1*PlotVar(i2,j1,k2,iVar)+&
                     dz2*PlotVar(i2,j1,k1,iVar)))+&
                     dx2*(   dy1*(   dz1*PlotVar(i1,j2,k2,iVar)+&
                     dz2*PlotVar(i1,j2,k1,iVar))+&
                     dy2*(   dz1*PlotVar(i1,j1,k2,iVar)+&
                     dz2*PlotVar(i1,j1,k1,iVar)))
             end do
             
             rplot_out = rplot
             if (plot_dimensional(ifile)) rplot_out = rplot_out*No2Io_V(UnitX_)
             phi_out = phi_plot*180.0/cPi        !output in degrees
             theta_out = theta_plot*180.0/cPi    !output in degrees
             
             if (theta_plot <= (cPi/2.0+.1*dtheta_plot)) then
                nBLKcellsN=nBLKcellsN+1
                select case(plot_form(ifile))
                case('tec')
                   write(unit_tmp,'(i7,i7, 30(E14.6))')&
                        i,j,x,y,z,theta_out,phi_out,PointVar(1:nplotvar)
                case('idl')
                   dxblk = 1.0  ! just chosen to get the scaling right.  
                   ! since this is R and R is constant.
                   if(save_binary)then
                      write(unit_tmp)&
                           dxblk,rplot_out,theta_out,phi_out,PointVar(1:nplotvar)
                   else
                      write(unit_tmp,'(20(1pe13.5))')&
                           dxblk,rplot_out,theta_out,phi_out,PointVar(1:nplotvar)
                   endif
                case('hdf')
                   call write_sph_var_hdf5(i, j, H5Index, rplot_out, theta_out,phi_out,&
                        dtheta_plot, dphi_plot, nTheta, nPhi, nPlotVar, PointVar, &
                        plotVarNames)
                   
                end select
             end if
             if (theta_plot >= (cPi/2.0-.1*dtheta_plot)) then
                theta_out = 180.- theta_out      !output in degrees
                nBLKcellsS=nBLKcellsS+1
                select case(plot_form(ifile))
                case('tec')
                   write(unit_tmp2,'(i7,i7, 30(E14.6))')&
                        ntheta-i+1,j,x,y,z,theta_out,phi_out,PointVar(1:nplotvar)
                case('idl')
                   dxblk = 1.0  ! just chosen to get the scaling right.  
                   ! since this is R and R is constant.
                   if(save_binary)then
                      write(unit_tmp2)&
                           dxblk,rplot_out,theta_out,phi_out,PointVar(1:nplotvar)
                   else
                      write(unit_tmp2,'(20(1pe13.5))')&
                           dxblk,rplot_out,theta_out,phi_out,PointVar(1:nplotvar)
                   endif
                case('hdf')
                   call write_sph_var_hdf5(i, j, H5Index, rplot_out, theta_out,phi_out,&
                        dtheta_plot, dphi_plot, nTheta, nPhi, nPlotVar, PointVar, &
                        plotVarNames)
                   
                end select
             end if
             H5Index = H5Index + 1
             
          end if
          
       end do
    end do
end subroutine write_plot_sph

!=============================================================================
end module ModPlotShell
