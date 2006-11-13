!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_plot_idl(ifile,iBLK,nplotvar,plotvar, &
     xmin,xmax,ymin,ymax,zmin,zmax, &
     dxblk,dyblk,dzblk,nBLKcells)

  ! Save all cells within plotting range, for each processor


  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,PROCtest,BLKtest,test_string,x_,y_,z_,Phi_
  use ModGeometry,ONLY:x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK,&
       is_axial_geometry,x1,x2,y1,y2,z1,z2,&          !^CFG IF COVARIANT
       XyzStart_BLK,XyzMin_D,XyzMax_D
  use ModPhysics, ONLY : unitUSER_x
  use ModIO
  use ModNumConst
  implicit none

  ! Arguments

  integer, intent(in) :: ifile, iBLK
  integer, intent(in) :: nplotvar
  real, intent(in) :: PlotVar(-1:nI+2,-1:nJ+2,-1:nK+2,nplotvar)
  real, intent(in)     :: xmin,xmax,ymin,ymax,zmin,zmax
  real, intent(inout)  :: dxblk,dyblk,dzblk
  integer, intent(out) :: nBLKcells

  ! Local variables
  ! Indices and coordinates
  integer :: i,j,k,i2,j2,k2,imin,imax,jmin,jmax,kmin,kmax,r,rx,ry,rz,iVar
  real :: x,y,z,dx,r_3
  real :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1
  real :: dxblk_out

  real :: ySqueezed

  real,parameter:: cHalfMinusTiny=cHalf*(cOne-cTiny)

  logical :: oktest,oktest_me
  !---------------------------------------------------------------------------

  if(iProc==PROCtest .and. iBLK==BLKtest)then
     call set_oktest('write_plot_idl',oktest,oktest_me)
  else
     oktest=.false.; oktest_me=.false.
  end if


  if(index(test_string,'SAVEPLOTALL')>0)then
     ! Save all cells of block including ghost cells

     nBLKcells=0
     dxblk=dx_BLK(iBLK); dyblk=dy_BLK(iBLK); dzblk=dz_BLK(iBLK)
     dxblk_out = dxblk
     if (plot_dimensional(ifile)) dxblk_out = dxblk*unitUSER_x

     plot_range(1,ifile)=XyzMin_D(1)-2*dxblk
     plot_range(2,ifile)=XyzMax_D(1)+2*dxblk
     plot_range(3,ifile)=XyzMin_D(2)-2*dyblk
     plot_range(4,ifile)=XyzMax_D(2)+2*dyblk
     plot_range(5,ifile)=XyzMin_D(3)-2*dzblk
     plot_range(6,ifile)=XyzMax_D(3)+2*dzblk
     plot_dx(1,ifile)=dxblk
     plot_dx(2,ifile)=dyblk
     plot_dx(3,ifile)=dzblk

     do k=-1,nK+2; do j=-1,nJ+2; do i=-1,nI+2
        nBLKcells=nBLKcells+1
        if (plot_dimensional(ifile)) then
           x = x_BLK(i,j,k,iBLK)*unitUSER_x
           y = y_BLK(i,j,k,iBLK)*unitUSER_x
           z = z_BLK(i,j,k,iBLK)*unitUSER_x
        else
           x = x_BLK(i,j,k,iBLK)
           y = y_BLK(i,j,k,iBLK)
           z = z_BLK(i,j,k,iBLK)
        end if
        if(save_binary)then
           write(unit_tmp)dxblk_out,x,y,z,PlotVar(i,j,k,1:nplotvar)
        else
           write(unit_tmp,'(20(1pe13.5))')dxblk_out,x,y,z,&
                PlotVar(i,j,k,1:nplotvar)
        endif
     end do; end do; end do

     return

  end if

  ! The range for the cell centers is dx/2 wider
  xmin1=xmin-cHalfMinusTiny*dx_BLK(iBLK)
  xmax1=xmax+cHalfMinusTiny*dx_BLK(iBLK)
  ymin1=ymin-cHalfMinusTiny*dy_BLK(iBLK)
  ymax1=ymax+cHalfMinusTiny*dy_BLK(iBLK)
  zmin1=zmin-cHalfMinusTiny*dz_BLK(iBLK)
  zmax1=zmax+cHalfMinusTiny*dz_BLK(iBLK)

  nBLKcells = 0
  if(is_axial_geometry())then                  !^CFG IF COVARIANT BEGIN
     ! Make sure that angles around 3Pi/2 are moved to Pi/2 for x=0 cut
     ySqueezed = mod(xyzStart_BLK(Phi_,iBLK),cPi)
     ! Make sure that small angles are moved to Pi degrees for y=0 cut
     if(ySqueezed < 0.25*cPi) ySqueezed = ySqueezed + cPi
  else                                          !^CFG END COVARIANT
     ySqueezed = xyzStart_BLK(y_,iBLK)
  end if                                        !^CFG IF COVARIANT 

  ! If block is fully outside of cut then cycle
  if(xyzStart_BLK(x_,iBLK)>xmax1.or.&
       xyzStart_BLK(x_,iBLK)+(nI-1)*dx_BLK(iBLK)<xmin1.or.&
       ySqueezed>ymax1.or.&
       ySqueezed+(nJ-1)*dy_BLK(iBLK)<ymin1.or.&  
       xyzStart_BLK(z_,iBLK)>zmax1.or.&
       xyzStart_BLK(z_,iBLK)+(nK-1)*dz_BLK(iBLK)<zmin1)&
       return

  dx = plot_dx(1,ifile)
  dxblk=dx_BLK(iBLK); dyblk=dy_BLK(iBLK); dzblk=dz_BLK(iBLK)

  ! Calculate index limits of cells inside cut
  imin=max(1 ,floor((xmin1-xyzStart_BLK(x_,iBLK))/dxblk)+2)
  imax=min(nI,floor((xmax1-xyzStart_BLK(x_,iBLK))/dxblk)+1)

  jmin=max(1 ,floor((ymin1-ySqueezed)/dyblk)+2)
  jmax=min(nJ,floor((ymax1-ySqueezed)/dyblk)+1)

  kmin=max(1 ,floor((zmin1-xyzStart_BLK(z_,iBLK))/dzblk)+2)
  kmax=min(nK,floor((zmax1-xyzStart_BLK(z_,iBLK))/dzblk)+1)


  if(oktest_me.and.iBLK==BLKtest)then
     write(*,*)'imin,imax,jmin,jmax,kmin,kmax=',&
          imin,imax,jmin,jmax,kmin,kmax
     write(*,*)'dxblk,x1,y1,z1',dxblk,xyzStart_BLK(:,iBLK)
     write(*,*)'ySqueezed  =',ySqueezed
     write(*,*)'xmin1,xmax1=',xmin1,xmax1
     write(*,*)'ymin1,ymax1=',ymin1,ymax1
     write(*,*)'zmin1,zmax1=',zmin1,zmax1
  end if

  if(dxblk>=dx)then
     ! Cell is equal or coarser than dx, save all cells in cut
     dxblk_out = dxblk
     if (plot_dimensional(ifile))dxblk_out = dxblk_out*unitUSER_x
     do k=kmin,kmax
        do j=jmin,jmax
           do i=imin,imax
              x = x_BLK(i,j,k,iBLK)
              y = y_BLK(i,j,k,iBLK)
              z = z_BLK(i,j,k,iBLK)

              if(x<x1.or.x>x2.or.y<y1.or.y>y2.or.&                !^CFG IF COVARIANT
                   z<z1.or.z>z2) CYCLE                            !^CFG IF COVARIANT

              if (plot_dimensional(ifile)) then
		 x = x*unitUSER_x
		 y = y*unitUSER_x
		 z = z*unitUSER_x
              end if
              if(save_binary)then
                 write(unit_tmp)dxblk_out,x,y,z,&
   		      PlotVar(i,j,k,1:nplotvar)
              else
                 write(unit_tmp,'(20(1pe13.5))')dxblk_out,x,y,z,& 
                      PlotVar(i,j,k,1:nplotvar)
              endif
              nBLKcells=nBLKcells+1
           end do
        end do
     end do
  else
     ! Block is finer then required resolution
     ! Calculate restriction factor
     r=min(nI,nint(dx/dxblk))

     ! Calclulate restricted cell size
     dxblk=r*dxblk; dyblk=r*dyblk; dzblk=r*dzblk
     dxblk_out = dxblk
     if (plot_dimensional(ifile))dxblk_out = dxblk_out*unitUSER_x

     ! Restriction is limited by the width of the plotting region
     rx=min(imax-imin+1,r)
     ry=min(jmax-jmin+1,r)
     rz=min(kmax-kmin+1,r)

     ! Factor for taking the average
     r_3=1./(rx*ry*rz)

     if(oktest_me.and.iBLK==BLKtest)then
        write(*,*)'r,rx,ry,rz,r_3=',r,rx,ry,rz,r_3
     end if

     ! Loop for the rx*ry*rz bricks inside the cut
     do k=kmin,kmax,rz
        k2=k+rz-1
        do j=jmin,jmax,ry
           j2=j+ry-1
           do i=imin,imax,rx
              i2=i+rx-1
              x=0.5*(x_BLK(i,j,k,iBLK)+x_BLK(i2,j2,k2,iBLK))
              y=0.5*(y_BLK(i,j,k,iBLK)+y_BLK(i2,j2,k2,iBLK))
              z=0.5*(z_BLK(i,j,k,iBLK)+z_BLK(i2,j2,k2,iBLK))

              if(x<x1.or.x>x2.or.y<y1.or.y>y2.or.&                !^CFG IF COVARIANT
                   z<z1.or.z>z2) CYCLE                            !^CFG IF COVARIANT 

              if(plot_dimensional(ifile))then
                 x=x*unitUSER_x
                 y=y*unitUSER_x
                 z=z*unitUSER_x
              end if
              if(save_binary)then
                 write(unit_tmp)dxblk_out,x,y,z,&
                      (r_3*sum(PlotVar(i:i2,j:j2,k:k2,iVar)),iVar=1,nplotvar)
              else
                 write(unit_tmp,'(20(1pe13.5))')dxblk_out,x,y,z,&
                      (r_3*sum(PlotVar(i:i2,j:j2,k:k2,iVar)),iVar=1,nplotvar)
              endif
              nBLKcells=nBLKcells+1
           end do
        end do
     end do
  end if

end subroutine write_plot_idl
