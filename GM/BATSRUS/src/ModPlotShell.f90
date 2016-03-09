!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!=============================================================================
module ModPlotShell

  use ModIO
  use ModNumConst, ONLY: cRadtoDeg, cDegToRad, cPi

  implicit none
  save

  private   ! except
  public:: write_plot_sph
  public:: set_plot_shell
  public:: write_plot_shell
  
  ! Size of current plot:
  integer :: nR, nLon, nLat

  ! Ranges for current plot:
  real :: dR, dLat, dLon
  real :: rMin, rMax, LonMin, LonMax, LatMin, LatMax

  ! Local results container:
  ! Array of values written to file:
  real, allocatable :: PlotVar_VIII(:,:,:,:)

  ! Coordinate conversion matrix
  real:: PlotToGm_DD(3,3)

contains
  !============================================================================
  real function minmod(x, y)
    real, intent(in) :: x, y
    minmod = max(0.0, min(abs(x), sign(1.0,x)*y))
  end function minmod

  !============================================================================
  subroutine set_plot_shell(iFile, iBlock, nPlotvar, Plotvar_GV)

    ! Interpolate the plot variables for block iBlock
    ! onto the spherical shell of the plot area of plot file iFile.
    
    use ModMain,        ONLY: time_simulation, TypeCoordSystem
    use CON_axes,       ONLY: transform_matrix
    use ModGeometry,    ONLY: rMin_BLK
    use ModMain,        ONLY: BlkTest
    use ModInterpolate, ONLY: trilinear
    use BATL_lib,       ONLY: CoordMin_DB, nIjk_D, CellSize_DB, xyz_to_coord
    use ModCoordTransform, ONLY: rlonlat_to_xyz

    ! Arguments
    integer, intent(in) :: iFile, iBlock
    integer, intent(in) :: nPlotvar
    real,    intent(in) :: PlotVar_GV(MinI:MaxI,MinJ:MaxJ,MinK:MaxK,nPlotVar)
    
    ! Local variables
    integer :: i, j, k, iVar
    
    real :: r, Lon, Lat
    real :: XyzPlot_D(3), XyzGm_D(3)
    real :: Coord_D(3), CoordNorm_D(3)
    
    character(len=*), parameter :: NameSub='set_plot_shell'
    logical :: DoTest, DoTestMe
    !--------------------------------------------------------------------------
    ! Check testing for block
    if(iBlock == BlkTest)then
       call CON_set_do_test(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if
    
    ! Allocate results array and set up all the spherical shell
    if(.not.allocated(PlotVar_VIII)) then

       ! Get plot area info from ModIO arrays:
       dR     = abs(plot_dx(1, iFile))
       dLon   = plot_dx(2, iFile) * cDegtoRad
       dLat   = plot_dx(3, iFile) * cDegtoRad
       rMin   = plot_range(1, iFile)
       rMax   = plot_range(2, iFile)
       LonMin = plot_range(3, iFile) * cDegtoRad
       LonMax = plot_range(4, iFile) * cDegtoRad
       LatMin = plot_range(5, iFile) * cDegtoRad
       LatMax = plot_range(6, iFile) * cDegtoRad
    
       ! Set number of points:
       nR   = nint((rMax - rMin)/dR)       + 1
       nLon = nint((LonMax - LonMin)/dLon) + 1
       nLat = nint((LatMax - LatMin)/dLat) + 1
    
       ! Ensure dR, dLon and dLat are compatible with the ranges
       dR   = (rMax - rMin)/max(1, nR - 1)
       dLon = (LonMax - LonMin)/max(1, nLon - 1)
       dLat = (LatMax - LatMin)/max(1, nLat - 1)

       ! The 0 element is to count the number of blocks that
       ! contribute to a plot variable.
       allocate(PlotVar_VIII(0:nPlotVar,nR,nLon,nLat))
       PlotVar_VIII = 0.0

       ! Get coordinate transformation matrix:
       PlotToGm_DD = transform_matrix(Time_Simulation, &
            TypeCoordPlot_I(iFile), TypeCoordSystem)
    end if
     
    if (DoTestMe) then
       write(*,*) NameSub//' Called for iBlock=      ', iBlock
       write(*,*) NameSub//' Raw plot_dx=          ', plot_dx(:,iFile)
       write(*,*) NameSub//' Raw plot_range=       ', plot_range(:,iFile)
       write(*,*) NameSub//' dR, dLon, dLat =      ', dR, dLon, dLat
       write(*,*) NameSub//' r, Lon, Lat range = ',  &
            rMin, rMax, LonMin,LonMax,LatMin,LatMax
       write(*,*) NameSub,' nR, nLon, nLat, dR, dLon, dLat = ', &
            nR, nLon, nLat, dR, dLon, dLat
    end if

    ! Loop through shell points and interpolate PlotVar
    do i = 1, nR
       r = rMin + (i-1)*dR

       if(r < rMin_BLK(iBlock)) CYCLE
       ! if(r > rMax_B(iBlock)) CYCLE

       do j = 1, nLon
          Lon = LonMin + (j-1)*dLon
          do k = 1, nLat
             Lat = LatMin + (k-1)*dLat

             ! Convert to Cartesian coordinates
             call rlonlat_to_xyz(r,Lon,Lat,XyzPlot_D)

             ! Convert from plot coordinates to BATSRUS grid:
             XyzGm_D = matmul(PlotToGm_DD, XyzPlot_D)

             ! Convert to generalized coordinates
             call xyz_to_coord(XyzGm_D, Coord_D)

             ! Normalize to block coordinates
             CoordNorm_D = &
                  (Coord_D - CoordMin_DB(:,iBlock))/CellSize_DB(:,iBlock) + 0.5

             ! Check if point is inside block
             if(any(CoordNorm_D < 0.5)) CYCLE
             if(any(CoordNorm_D > nIjk_D + 0.5)) CYCLE

             ! compute the interpolated values at the current location
             PlotVar_VIII(0,i,j,k) = 1.0
             do iVar=1, nPlotVar
                ! Interpolate up to ghost cells.
                PlotVar_VIII(iVar,i,j,k) = trilinear(PlotVar_GV(:,:,:,iVar),&
                     MinI, MaxI, MinJ, MaxJ, MinK, MaxK, CoordNorm_D)
             end do
             
          end do ! lon loop
       end do    ! lat loop
    end do       ! r loop

  end subroutine set_plot_shell

  !============================================================================
  subroutine write_plot_shell(iFile, nPlotVar, NameVar_V, &
       NameUnit, NameFile)

    ! Collect results from all blocks and write to single output file.
    use ModMpi
    use ModMain,     ONLY: time_simulation, n_step
    use ModProcMH,   ONLY: iProc, nProc, iComm
    use ModPlotFile, ONLY: save_plot_file

    integer,          intent(in) :: iFile, nPlotvar
    character(len=*), intent(in) :: NameFile, NameVar_V(nPlotVar), NameUnit

    integer :: iVar, iR, iLon, iLat, iError
    character(len=500) :: NameVar

    logical :: DoTest,DoTestMe
    character(len=*), parameter :: NameSub='write_plot_shell'
    !--------------------------------------------------------------------------
    call set_oktest(NameSub, DoTest, DoTestMe)

    ! This subroutine does not support HDF output.
    if(plot_form(iFile) == 'hdf') call CON_stop(NameSub// &
         ': HDF file type not supported for Geo Sphere output.')
    
    ! Collect results to head node
    if(nProc > 1) call MPI_reduce(MPI_IN_PLACE, PlotVar_VIII, &
         size(PlotVar_VIII), MPI_REAL, MPI_SUM, 0, iComm, iError)

    ! Save results to disk
    if(iProc==0) then
       ! Build a single-line list of variable names.
       NameVar = 'r lon lat'
       do iVar=1, nPlotVar
          NameVar = trim(NameVar)  // ' ' // trim(NameVar_V(iVar))
       end do

       do iR = 1, nR; do iLon = 1, nLon; do iLat = 1, nLat
          if(PlotVar_VIII(0,iR,iLon,iLat) <= 1.0) CYCLE
          PlotVar_VIII(1:nPlotVar,iR,iLon,iLat) = &
               PlotVar_VIII(1:nPlotVar,iR,iLon,iLat) &
               / PlotVar_VIII(0,iR,iLon,iLat)
       end do; end do; end do

       ! Call save_plot_file to write data to disk.
       call save_plot_file(NameFile, &
            TypeFileIn=TypeFile_I(iFile), &
            StringHeaderIn=NameUnit, &
            nStepIn=n_step, &
            TimeIn=time_simulation, &
            NameVarIn = NameVar, &
            CoordMinIn_D = (/rMin, cRadtoDeg*LonMin, cRadtoDeg*LatMin/), &
            CoordMaxIn_D = (/rMax, cRadtoDeg*LonMax, cRadtoDeg*LatMax/), &
            VarIn_VIII = PlotVar_VIII(1:,:,:,:))
    end if
    
    ! Deallocate results arrays:.
    deallocate(PlotVar_VIII)
    
  end subroutine write_plot_shell

  !============================================================================
  subroutine write_plot_sph(iFile,iBLK,nPlotvar,Plotvar, &
       nTheta,nPhi,rPlot,plotvarnames,H5Index,nBlkCellsN,nBlkCellsS)
    ! CANDIDATE FOR REMOVAL.  SHELL PLOTTING TOOLS PREFERRED.
    ! Save all cells within plotting range, for each processor
    
    use ModPhysics, ONLY : No2Io_V, UnitX_
    use ModGeometry, ONLY: nI,nJ,nK, Xyz_DGB, CellSize_DB 
    use ModMain,     ONLY: BlkTest
    use ModHdf5, ONLY: write_sph_var_hdf5
    
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
    !--------------------------------------------------------------------------
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
