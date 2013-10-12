!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!This code is a copyright protected software (c) 2002- University of Michigan
!=============================================================================
subroutine write_plot_sph(iFile,iBLK,nPlotvar,Plotvar, &
     nTheta,nPhi,rPlot,plotvarnames,H5Index,nBlkCellsN,nBlkCellsS)

  ! Save all cells within plotting range, for each processor

  use ModPhysics, ONLY : No2Io_V, UnitX_
  use ModGeometry, ONLY: nI,nJ,nK, Xyz_DGB, CellSize_DB 
  use ModMain,     ONLY: BlkTest
  use ModNumConst
  use ModIO
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

contains

  real function minmod(x,y)
    real, intent(in) :: x,y
    minmod = max(cZero,min(abs(x),sign(cOne,x)*y))
  end function minmod

end subroutine write_plot_sph
