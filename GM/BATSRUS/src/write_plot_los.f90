!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_plot_los(ifile)
  ! Purpose:  creates a synthic coronagraph image of Thomson scattered 
  !           white light by inegrating light scattered in the line of sight.
  !           This is accomplished by looping over all blocks per processor 
  !           and looping over lines of sight and then integrating on those
  !           lines.
  !           Written by Chip Manchester, KC Hansen
  !                   some improvements by Gabor Toth
  !           July     2001
  !           January  2002 modified for improved image plane
  !           December 2003 fixed sign error in scattering coefficient b_los
  !           January  2004 fix to accept 0 in LOS vector
  !           January  2004 fixed declaration for norm_los(3), r_Pix(3) 
  !           February 2004 fix integration and make 2nd order accurate
  !                         fix save_file in mainf.90 update ghost cells
  !                         include forgotten plot_pars1=plot_pars(ifile)
  !                         flags for PB and LW and limb darkening parameter
  !                         improved block-line distance calculation
  !                         fix IDL output and plot filenames
  !                         use dynamic arrays and contained subroutines
  !                         exclude whole block if outside LOS image
  !                         simplify block-line distance calculation
  !                         simplify body-line distance calculation
  !                         moved plot variable loop inside line integration
  !                         dimensionalize some of the plot variables

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,time_simulation,unusedBLK, &
       time_accurate,TimeH4,TimeM2,TimeS2
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics, ONLY : unitUSER_x
  use ModIO
  use ModAdvance, ONLY : rho_,rhoUx_,State_VGB
  use ModNumConst, ONLY : cTiny
  use ModMpi
  implicit none

  ! Arguments

  integer, intent(in) :: ifile

  ! Local variables

  integer :: iError

  ! Plot variables
  integer, parameter :: neqparmax=10
  real, allocatable :: PlotVar(:,:,:), PlotBLK(:,:,:), los_image(:,:,:)

  real ::     eqpar(neqparmax)
  character (len=10) :: eqparnames(neqparmax)
  character (len=10) :: plotvarnames(nplotvarlosmax)

  integer :: neqpar, nplotvar
  integer :: iPix, jPix
  real    :: dx_Pix, x_Pix, dy_Pix, y_Pix 
  real    :: norm_los(3), pos(3), a_Pix(3), b_Pix(3), r_Pix(3)
  real    :: x_los_blk, y_los_blk, z_los_blk, los_to_blk
  real    :: x_los_org, y_los_org, z_los_org, los_to_org 
  real    :: radius_of_block
  real    :: a_mag, b_mag

  real    :: XyzBlockCenter_D(3), xLosBlock, yLosBlock

  character (LEN=79) :: allnames, StringHeadLine
  character (LEN=500) :: unitstr_TEC, unitstr_IDL
  character (LEN=4) :: file_extension
  character (LEN=40) :: file_format

  ! block and variable Indices
  integer :: iBLK,iVar

  logical :: oktest,oktest_me,DoTiming,DoTimingMe
  !---------------------------------------------------------------------------

  ! Initialize stuff
  call set_oktest('write_plot_los',oktest,oktest_me)
  call set_oktest('los_timing',DoTiming,DoTimingMe)

  call timing_start('write_plot_los')

  where(los_vector(:,ifile)==0.0) los_vector(:,ifile)=cTiny
  norm_los = los_vector(:,ifile)/sqrt(sum(los_vector(:,ifile)**2))

  dx_Pix = x_size_image/real(n_pix_X -1) 
  dy_Pix = y_size_image/real(n_pix_Y -1)

  if(oktest .and. iProc==0) then
     write(*,*) 'vector',los_vector(:,ifile)
     write(*,*) 'norm_los', norm_los
     write(*,*) 'x_size_image ',x_size_image
     write(*,*) 'y_size_image ',y_size_image
     write(*,*) 'xoffset,yoffset',xoffset,yoffset 
     write(*,*) 'pos',pos
     write(*,*) 'dx,dy',dx_Pix,dy_Pix
     write(*,*) 'nx,ny',n_pix_X,n_pix_Y
  end if

  unitstr_TEC = ''
  unitstr_IDL = ''

  plot_type1=plot_type(ifile)
  plot_vars1=plot_vars(ifile)
  plot_pars1=plot_pars(ifile)

  if(oktest_me)write(*,*)'ifile=',ifile,' plot_type=',plot_type1, &
       ' form = ',plot_form(ifile)

  call split_str(plot_vars1,nplotvarlosmax,plotvarnames,nplotvar)
  call split_str(plot_pars1,neqparmax,eqparnames,neqpar)
  call set_eqpar(ifile-plot_,neqpar,eqparnames,eqpar)

  allnames='x y '//trim(plot_vars1)//' '//plot_pars1

  if(oktest_me) then
     write(*,*) plot_vars1
     write(*,*) nplotvar,plotvarnames(1:nplotvar)
  end if

  ! Get the headers that contain variables names and units
  select case(plot_form(ifile))
  case('tec')
     call get_TEC_los_variables(ifile,nplotvar,plotvarnames,unitstr_TEC)
     if(oktest .and. iProc==0) write(*,*)unitstr_TEC
  case('idl')
     call get_IDL_los_units(ifile,nplotvar,plotvarnames,unitstr_IDL)
     if(oktest .and. iProc==0) write(*,*)unitstr_IDL
  end select
  a_mag = sqrt(norm_los(1)**2 + norm_los(2)**2)
  a_Pix(1) =  norm_los(2)/a_mag
  a_Pix(2) = -norm_los(1)/a_mag
  a_Pix(3) =  0.0
  b_mag = sqrt((norm_los(1)*norm_los(3))**2 + (norm_los(2)*norm_los(3))**2 &
       + (norm_los(1)**2 + norm_los(2)**2)**2)
  b_Pix(1) = -norm_los(1)*norm_los(3)/b_mag
  b_Pix(2) = -norm_los(2)*norm_los(3)/b_mag
  b_Pix(3) = (norm_los(1)**2 + norm_los(2)**2)/b_mag

  pos(1) = -xoffset*a_Pix(1) -yoffset*b_Pix(1) 
  pos(2) = -xoffset*a_Pix(2) -yoffset*b_Pix(2)
  pos(3) = -xoffset*a_Pix(3) -yoffset*b_Pix(3)

  allocate( &
       PlotVar(n_pix_X,n_pix_Y,nplotvar), &
       PlotBLK(n_pix_X,n_pix_Y,nplotvar), &
       los_image(n_pix_X,n_pix_Y,nplotvar))

  PlotVar = 0.0

  if(DoTiming)call timing_start('los_block_loop')

  ! loop over blocks
  do iBLK = 1,nBLK

     if (unusedBLK(iBLK)) CYCLE

     radius_of_block = 0.5*sqrt(&
          ((nI+1)*dx_BLK(iBLK))**2 + &
          ((nJ+1)*dy_BLK(iBLK))**2 + &
          ((nK+1)*dz_BLK(iBLK))**2)

     !position of the block center
     XyzBlockCenter_D(1) = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(2) = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(3) = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))

     ! calculate position of block center on the los image
     XyzBlockCenter_D = XyzBlockCenter_D - Pos
     xLosBlock = dot_product(XyzBlockCenter_D,a_pix)
     yLosBlock = dot_product(XyzBlockCenter_D,b_pix)

     ! Check if block is inside the LOS image
     if(xLosBlock < -radius_of_block) CYCLE
     if(xLosBlock > +radius_of_block + x_size_image) CYCLE
     if(yLosBlock < -radius_of_block) CYCLE
     if(yLosBlock > +radius_of_block + y_size_image) CYCLE

     ! Initialize plot variable for this block
     PlotBLK = 0.0

     ! Loop over pixels
     do jPix=1,n_pix_Y
        ! Y position of the pixel on the image plane
        y_Pix = real(jPix -1)*dy_Pix

        ! Check if block can intersect this pixel
        if(abs(y_Pix-yLosBlock)>radius_of_block)CYCLE

        do iPix=1,n_pix_X

           ! X position of the pixel on the image plane
           x_Pix = real(iPix -1)*dx_Pix

           ! Check if block can intersects this pixel
           if( (x_Pix-xLosBlock)**2 + (y_Pix-yLosBlock)**2 > &
                radius_of_block**2 ) CYCLE

           ! Check if pixel is within occultation radius
           if((x_Pix - xOffset)**2 + (y_Pix - yOffset)**2 <= radius_occult**2)&
                CYCLE

           ! Calculate contribution of this block to this pixel
           call set_plotvar_los(ifile-plot_)

        end do ! jPix loop
     end do ! iPix loop

     ! sum over blocks on a pe, then sum over pe's
     PlotVar = PlotVar + PlotBLK
  end do       !iBLK

  if(DoTiming)call timing_stop('los_block_loop')

  if (plot_dimensional(ifile)) call dimensionalize_plotvar_los(ifile-plot_)

  ! collect the pixels on one node and then write out the file 
  if(nProc>1)then 
     call MPI_REDUCE(PlotVar,los_image,n_pix_x*n_pix_y*nplotvar, &
          MPI_REAL,MPI_SUM,0,iComm,iError)
  else
     los_image=PlotVar
  end if

  if(DoTiming)call timing_start('los_save_plot')
  if (iProc==0) then

     select case(plot_form(ifile))
     case('tec')
        file_extension='.dat'
     case('idl')
        file_extension='.out'
     end select

     if (ifile-plot_ > 9) then
        file_format='("' // trim(NamePlotDir) // '",a,i2,a,i7.7,a)'
     else
        file_format='("' // trim(NamePlotDir) // '",a,i1,a,i7.7,a)'
     end if

     if(time_accurate)then
        call gettimestring
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_t"//TimeH4//TimeM2//TimeS2//"_n",n_step,&
             file_extension
     else
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_n",n_step,file_extension
     end if
     open(unit_tmp,file=filename,status="unknown",err=999)

     ! write header file
     select case(plot_form(ifile))
     case('tec')
        write(unit_tmp,*) 'TITLE="BATSRUS: Synthetic Image"'
        write(unit_tmp,'(a)')trim(unitstr_TEC)
        write(unit_tmp,*) 'ZONE T="LOS Image"', &
             ', I=',n_pix_X,', J=',n_Pix_Y,', K=1, F=POINT'
        ! Write point values

        do iPix=1,n_pix_X
           do jPix=1,n_pix_Y
              x_Pix = real(iPix-1)*dx_Pix
              y_Pix = real(jPix-1)*dy_Pix

              if (plot_dimensional(ifile)) then
                 write(unit_tmp,fmt="(30(E14.6))") x_Pix*unitUSER_x, &
                      y_Pix*unitUSER_x,los_image(iPix,jPix,1:nPlotVar)
              else
                 write(unit_tmp,fmt="(30(E14.6))") x_Pix,y_Pix, &
                      los_image(iPix,jPix,1:nPlotVar)
              end if

           end do
        end do

     case('idl')
        ! description of file contains units, physics and dimension
        StringHeadLine = 'LOS integrals_var22'
        write(unit_tmp,"(a)") StringHeadLine

        ! 2 in the next line means 2 dimensional plot, 1 in the next line
        !  is for a 2D cut, in other words one dimension is left out)
        write(unit_tmp,"(i7,1pe13.5,3i3)") &
             n_step,time_simulation,2,neqpar,nplotvar

        ! Grid size
        write(unit_tmp,"(2i4)") n_pix_X, n_pix_Y

        ! Equation parameters
        write(unit_tmp,"(100(1pe13.5))") eqpar(1:neqpar)

        ! Coordinate, variable and equation parameter names
        write(unit_tmp,"(a)")allnames

        ! Data
        do jPix=1,n_pix_Y
           y_Pix = real(jPix-1)*dy_Pix

           do iPix=1,n_pix_X
              x_Pix = real(iPix-1)*dx_Pix

              if (plot_dimensional(ifile)) then
                 x_Pix = x_Pix * unitUSER_x
                 y_Pix = y_Pix * unitUSER_x
              end if

              write(unit_tmp,fmt="(30(1pe13.5))") &
                   x_Pix, y_Pix, los_image(iPix,jPix,1:nPlotVar)

           end do
        end do

     end select
     close(unit_tmp)
  end if  !iProc ==0
  if(DoTiming)call timing_stop('los_save_plot')

  deallocate(PlotVar, PlotBLK, los_image)

  if(oktest_me)write(*,*)'write_plot_los finished'

  call timing_stop('write_plot_los')

  return

999 continue

  call stop_mpi("Error in opening or writing file in write_plot_los")

contains
  !===========================================================================

  subroutine set_plotvar_los(iplotfile)

    use ModProcMH
    use ModIO
    use ModGeometry
    use ModIO, ONLY : mu_los, TypeLosImage
    implicit none

    integer, intent(in) :: iPlotFile

    ! Local variables
    integer :: iVar
    integer :: i, j, k, i_los, nline_seg, counter 
    real :: x_los, y_los, z_los, r_los  
    real :: x_q, y_q, z_q, q_mag
    real :: a_los, b_los, c_los, d_los
    real :: sin_omega, cos_omega, cos_theta
    real :: rho_los, s_los_sqrd, ds_los, direction
    real :: point_in(3), point_1(3), point_2(3), point_los(3) 
    real :: intrsct(2,3,3), face_location(2,3) 
    real :: xx1, xx2, yy1, yy2, zz1, zz2 

    !-------------------------------------------------------------------------
    !if(DoTiming)call timing_start('los_set_plotvar')

    ! Get the 3D location of the pixel
    r_Pix = pos + x_Pix*a_Pix + y_Pix*b_Pix

    !x_los, y_los, z_los, r_los give the position of the point on the los
    !mu_los parameter related to the limb darkening
    !face_location give the locations of the faces of the block
    !face_location(2,3) = x1, y1, z1---x2, y2, z2

    nline_seg = (nI + nJ + nK)

    !Determine the location of the block faces
    xx1 = 0.50*(x_BLK( 0, 0, 0,iBLK)+x_BLK(   1,   1  , 1,iBLK))
    xx2 = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(nI+1,nJ+1,nK+1,iBLK))
    yy1 = 0.50*(y_BLK( 0, 0, 0,iBLK)+y_BLK(   1,   1,   1,iBLK))
    yy2 = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(nI+1,nJ+1,nK+1,iBLK))
    zz1 = 0.50*(z_BLK( 0, 0, 0,iBLK)+z_BLK(   1,   1,   1,iBLK))
    zz2 = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(nI+1,nJ+1,nK+1,iBLK))

    face_location(1,1) = xx1
    face_location(1,2) = yy1
    face_location(1,3) = zz1
    face_location(2,1) = xx2
    face_location(2,2) = yy2
    face_location(2,3) = zz2

    !Determine where the line of sight enters and exits the block
    !loop over the number of block face pairs, face directions and coordinates
    do i=1,2       !face loop
       intrsct(i,1,1) = face_location(i,1)
       intrsct(i,2,2) = face_location(i,2)
       intrsct(i,3,3) = face_location(i,3)

       do j=1,3     !direction loop
          do k=1,3   !coordinate loop
             if (j /= k) then  
                intrsct(i,j,k) = r_Pix(k) + &
                     (norm_los(k)/norm_los(j))*(face_location(i,j) - r_Pix(j))
             end if
          end do
       end do
    end do

    !which of the 6 points are on the block?
    counter = 0
    do i=1,2 
       do j=1,3 
          if( (intrsct(i,j,1) >= xx1) .and. (intrsct(i,j,1) <= xx2)) then
             if( (intrsct(i,j,2) >= yy1) .and. (intrsct(i,j,2) <= yy2)) then
                if( (intrsct(i,j,3) >= zz1) .and. (intrsct(i,j,3) <= zz2)) then
                   counter = counter + 1
                   do k=1,3  !fixed..do loop was misplaced
                      if(i == 1) point_1(k) = intrsct(i,j,k)
                      if(i == 2) point_2(k) = intrsct(i,j,k)
                   end do
                end if
             end if
          end if
       end do
    end do

    ! Check if the los cuts through the block 
    if(counter == 2) then           !integrate
       s_los_sqrd    = 0.0
       direction     = 0.0

       do i=1,3 
          direction  = direction + (point_1(i) - point_2(i))*norm_los(i)
          s_los_sqrd = s_los_sqrd + (point_1(i) - point_2(i))**2 
       end do

       if(direction > 0) then 
          point_in = point_2
       else 
          point_in = point_1
       endif

       ds_los = sqrt(s_los_sqrd)/real(nline_seg)

       !if(DoTiming)call timing_start('los_integral')
       do i_los=1,nline_seg
          x_los = point_in(1) + (i_los-0.5)*ds_los*norm_los(1)
          y_los = point_in(2) + (i_los-0.5)*ds_los*norm_los(2)
          z_los = point_in(3) + (i_los-0.5)*ds_los*norm_los(3)
          r_los = sqrt(x_los**2 + y_los**2 + z_los**2)

          point_los(1) = x_los
          point_los(2) = y_los
          point_los(3) = z_los

          sin_omega = 1.0/r_los
          cos_omega = sqrt(1.0 - sin_omega**2)

          !omega and functions of omega are unique to a given line of sight
          a_los = cos_omega*sin_omega**2
          b_los = -0.125*( 1.0 - 3.0*sin_omega**2 - (cos_omega**2/sin_omega)* &
               (1.0 + 3.0*sin_omega**2)*log((1.0 + sin_omega)/cos_omega) )
          c_los = 4.0/3.0 - cos_omega - (1.0/3.0)*cos_omega**3
          d_los = 0.125*( 5.0 + sin_omega**2 - (cos_omega**2/sin_omega) * &
               (5.0 - sin_omega**2)*log((1.0 + sin_omega)/cos_omega) )

          z_q =   (norm_los(1)**2 + norm_los(2)**2)*z_los            &
               - norm_los(3)*(norm_los(1)*x_los + norm_los(2)*y_los)  
          x_q = x_los + (norm_los(1)/norm_los(3)) * (z_q - z_los)
          y_q = y_los + (norm_los(2)/norm_los(3)) * (z_q - z_los)
          q_mag = sqrt(x_q**2 + y_q**2 + z_q**2)

          cos_theta = q_mag/r_los

          ! interpolate density at this point with bilinear interpolation
          rho_los = point_value_los(point_los,iBLK)


          do iVar=1,nplotvar
             TypeLosImage = plotvarnames(iVar)
             select case(TypeLosImage)
             case ('len')
                ! Integrate the length of the integration lines
                PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + ds_los
                     
             !case('vlos','Vlos','ulos','Ulos')
             !   ! Integrate the velocity
             !   PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
             !        (rhoU_los/rho_los)*ds_los

             case('WL','wl')
                ! White light with limb darkening
                PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                     rho_los*( &
                     (1.0 - mu_los)*(2.0*c_los - a_los*cos_theta**2) &
                     + mu_los*(2.0*d_los - b_los*cos_theta**2) )*ds_los

             case('PB','pb')
                ! Polarization brightness
                PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                  rho_los*( (1.0 - mu_los)*a_los + mu_los*b_los) &
                  *cos_theta**2*ds_los

             case('rho','RHO')
                ! Simple density integral
                PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                     rho_los*ds_los

             case default
                PlotBLK(iPix,jPix,iVar)=-7777.
             end select
          end do ! iVar

       end do !line segment interation loop
       !if(DoTiming)call timing_stop('los_integral')

    end if !if counter = 2

    !if(DoTiming)call timing_stop('los_set_plotvar')

  end subroutine set_plotvar_los

  !============================================================================

  function point_value_los(Xpnt,iBLK)
    !Find the value at an arbitrary location X,Y,Z,iBLK

    ! Arguments
    integer, intent(in) :: iBLK
    real, intent(in) :: Xpnt(3)

    ! Local variables:
    real    :: point_value_los
    real    :: x(3)
    real    :: dx1,dy1,dz1,dx2,dy2,dz2
    integer :: i1,j1,k1,i2,j2,k2
    !--------------------------------------------------------------------------
    point_value_los = 0.0

    ! Convert to normalized coordinates (index and position are the same)
    x(1)=(Xpnt(1)-x_BLK(1,1,1,iBLK))/dx_BLK(iBLK)+1.
    x(2)=(Xpnt(2)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1.
    x(3)=(Xpnt(3)-z_BLK(1,1,1,iBLK))/dz_BLK(iBLK)+1.

    ! Determine cell indices corresponding to location qx
    i1=floor(x(1)); i2=i1+1
    j1=floor(x(2)); j2=j1+1
    k1=floor(x(3)); k2=k1+1

    ! Distance relative to the cell centers
    dx1=x(1)-real(i1); dx2=1.-dx1
    dy1=x(2)-real(j1); dy2=1.-dy1
    dz1=x(3)-real(k1); dz2=1.-dz1

    ! Bilinear interpolation in 3D
    point_value_los = &
         dx1*(   dy1*(   dz1*State_VGB(rho_,i2,j2,k2,iBLK)+&
         dz2*State_VGB(rho_,i2,j2,k1,iBLK))+&
         dy2*(   dz1*State_VGB(rho_,i2,j1,k2,iBLK)+&
         dz2*State_VGB(rho_,i2,j1,k1,iBLK)))+&
         dx2*(   dy1*(   dz1*State_VGB(rho_,i1,j2,k2,iBLK)+&
         dz2*State_VGB(rho_,i1,j2,k1,iBLK))+&
         dy2*(   dz1*State_VGB(rho_,i1,j1,k2,iBLK)+&
         dz2*State_VGB(rho_,i1,j1,k1,iBLK)))

  end function point_value_los

  !==========================================================================

  subroutine dimensionalize_plotvar_los(iplotfile)

    use ModPhysics, ONLY : unitUSER_x,unitUSER_U,unitSI_x,unitSI_rho
    use ModIO
    implicit none

    integer, intent(in) :: iPlotFile

    character (len=10) :: s

    integer :: iVar,i,j,k
    !--------------------------------------------------------------------------

    do iVar=1,nPlotVar
       s=plotvarnames(iVar)

       select case(s)
       case ('len')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*unitUSER_x
       case('rho')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*unitSI_rho*unitSI_x
       case('vlos','Vlos','ulos','Ulos')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*unitUSER_U
       case default
          ! no normalization
       end select

    end do ! iVar
  end subroutine dimensionalize_plotvar_los

end subroutine write_plot_los

!==============================================================================

subroutine get_TEC_los_variables(iFile,nplotvar,plotvarnames,unitstr_TEC)

  use ModPhysics, ONLY : unitstr_TEC_x,unitstr_TEC_U
  use ModIO, ONLY: plot_dimensional
  implicit none

  ! Arguments

  integer, intent(in) :: Nplotvar,iFile
  character (LEN=10), intent(in) :: plotvarnames(Nplotvar)
  character (len=500), intent(out) :: unitstr_TEC 
  character (len=10) :: s

  integer :: iVar, len


  !\
  ! This routine takes the plot_var information and loads the header file with
  ! the appropriate string of variable names and units
  !/

  if (plot_dimensional(ifile)) then
     write(unitstr_TEC,'(a)') 'VARIABLES = '
     write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'"X '//&
          trim(unitstr_TEC_x)
     write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'", "Y '//&
          trim(unitstr_TEC_x)
  else
     write(unitstr_TEC,'(a)') 'VARIABLES = "X", "Y'
  end if

  do iVar = 1, nplotvar

     write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'", "'

     s=plotvarnames(iVar)

     if (plot_dimensional(ifile)) then

        select case(s)
        case ('len')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'len'//' '//&
                trim(unitstr_TEC_x)
        case('rho')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'`r [m^-^2]'
        case('vlos','Vlos','ulos','Ulos')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'u.s'//' '//&
                trim(unitstr_TEC_U)
        case('wl')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'`wl [m^-^2]'//' '
        case('pb')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'`pb [m^-^2]'//' '

           ! DEFAULT FOR A BAD SELECTION
        case default
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'Default'

        end select

     else

        select case(s)
        case ('len')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'len'
        case('rho')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'`r'
        case('vlos','Vlos','ulos','Ulos')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'u.s'
        case('wl')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'wl'
        case('pb')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'pb'

           ! DEFAULT FOR A BAD SELECTION
        case default
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'Default'

        end select

     end if

  end do

  write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'"'

end subroutine get_TEC_los_variables

!======================================================================
subroutine get_IDL_los_units(ifile,nplotvar,plotvarnames,unitstr_IDL)

  use ModPhysics, ONLY : unitstr_IDL_x,unitstr_IDL_U
  use ModIO, ONLY : plot_dimensional

  implicit none

  ! Arguments

  integer, intent(in) :: iFile,Nplotvar
  character (LEN=10), intent(in) :: plotvarnames(Nplotvar)
  character (len=79), intent(out) :: unitstr_IDL 
  character (len=10) :: s

  integer :: iVar, len


  !\
  ! This routine takes the plot_var information and loads the header file with
  ! the appropriate string of unit values
  !/

  if (plot_dimensional(ifile)) then
     write(unitstr_IDL,'(a)') trim(unitstr_IDL_x)//' '//&
          trim(unitstr_IDL_x)//' '//&
          trim(unitstr_IDL_x)
  else
     write(unitstr_IDL,'(a)') 'normalized variables'
  end if

  if (plot_dimensional(ifile)) then

     do iVar = 1, nplotvar

        s=plotvarnames(iVar)

        select case(s)
        case ('len')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//&
                trim(unitstr_IDL_x)
        case('rho')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//'[m^-^2]'
        case('vlos','Vlos','ulos','Ulos')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//&
                trim(unitstr_IDL_U)
        case('wl')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//'[m^-^2]'
        case('pb')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//'[m^-^2]'
           ! DEFAULT FOR A BAD SELECTION
        case default
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//'" Dflt"'

        end select

     end do

  end if

end subroutine get_IDL_los_units
!==============================================================================
