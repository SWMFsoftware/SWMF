!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_plot_los(iFile)
  ! Purpose:  creates a synthic coronagraph image of Thomson scattered 
  !           white light by inegrating light scattered in the line of sight.
  !           This is accomplished by looping over all blocks per processor 
  !           and looping over lines of sight and then integrating on those
  !           lines.
  !           Written by Chip Manchester, KC Hansen
  !                   some improvements by Gabor Toth
  !                   some changes by Noe Lugaz
  !           July     2001
  !           January  2002 modified for improved image plane
  !           December 2003 fixed sign error in scattering coefficient b_los
  !           January  2004 fix to accept 0 in LOS vector
  !           January  2004 fixed declaration for norm_los(3), r_Pix(3) 
  !           February 2004 fix integration and make 2nd order accurate
  !                         fix save_file in main.f90 update ghost cells
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
  !           January 2006  compatibility with framework
  !                         rotation of the coordinates to HGI
  !                         wide-angle line-of-sight (for STEREO)
  !                         change in the parameters: satellite_position, 
  !                         offset_angle   
  !                         Cartesian grid and circular image centered
  !                         at the Sun (no offset angle)

  use ModProcMH
  use ModMain, ONLY : nI,nJ,nK,n_step,time_simulation,unusedBLK, &
       time_accurate,nBlock, NameThisComp,rBuffMax,TypeCoordSystem
  use ModGeometry, ONLY : x_BLK,y_BLK,z_BLK,dx_BLK,dy_BLK,dz_BLK
  use ModPhysics, ONLY : No2Io_V, UnitX_
  use ModIO
  use ModAdvance, ONLY : rho_, State_VGB
  use ModNumConst, ONLY : cTiny, cUnit_DD, cTolerance
  use ModMpi
  use CON_axes, ONLY : transform_matrix
  use ModCoordTransform, ONLY : rot_matrix_z
  use ModUtilities, ONLY: lower_case
  implicit none

  ! Arguments

  integer, intent(in) :: iFile

  ! Local variables

  integer :: iError

  ! File specific parameters
  integer :: nPix
  real    :: aOffset, bOffset, rSizeImage, rSizeImage2, rOccult, rOccult2,&
       OffsetAngle


  ! Plot variables
  integer, parameter :: neqparmax=10
  real, allocatable :: PlotVar(:,:,:), PlotBLK(:,:,:), los_image(:,:,:)

  real ::     eqpar(neqparmax)
  character (len=10) :: eqparnames(neqparmax)
  character (len=10) :: plotvarnames(nplotvarlosmax)
  character (len=10) :: NameVar


  integer :: nEqpar, nPlotvar
  integer :: iPix, jPix
  real    ::  x_Pix, y_Pix 
  real    :: ImageCenter_D(3), a_Pix(3), b_Pix(3), r_Pix(3),LosPix_D(3)
  real    :: rBlockSize, rBlockCenter
  real    :: a_mag, b_mag
  real    :: SizePix, r2Pix    
  real    :: BlockDistance, ObsDistance, Ratio
  real    :: XyzBlockCenter_D(3), xLosBlock, yLosBlock

  real, dimension(3,3) :: FromHgi_DD
  real, dimension(3) :: Los_D, ObsPos_D

  ! rInner in IH and rOuter in SC should be identical!
  ! rInner in SC should be smaller than the occulting radius
  ! rInner in IH should be larger than the inner boundary radius
  ! rOuter in SC should be smaller than the size of the domain
  real :: rInner, rInner2, rOuter, rOuter2

  character (LEN=79) :: allnames, StringHeadLine
  character (LEN=500) :: unitstr_TEC, unitstr_IDL
  character (LEN=4) :: file_extension
  character (LEN=40) :: file_format

  ! block and variable Indices
  integer :: iBLK, iVar

  logical :: oktest,oktest_me,DoTiming,DoTimingMe, DoCheckBlock
  logical :: UseScattering
  !---------------------------------------------------------------------------

  ! Initialize stuff
  call set_oktest('write_plot_los',oktest,oktest_me)
  call set_oktest('los_timing',DoTiming,DoTimingMe)

  call timing_start('write_plot_los')

  ! Set rInner and rOuter depending on component
  select case(NameThisComp)
  case('SC')
     rInner = 0.5 
     rOuter = rBuffMax
  case('IH')
     rInner = rBuffMax
     rOuter = 1000.0
  case('GM')
     rInner = 0.0 ! needed for comet applications
     rOuter = 1000.0
  end select
  rInner2 = rInner**2
  rOuter2 = rOuter**2

  if(NameThisComp == 'GM')then
     ! Do not convert to HGI
     FromHgi_DD = cUnit_DD
  else
     ! Convert to HGI
     FromHgi_DD = transform_matrix(Time_Simulation,'HGI', TypeCoordSystem)
  end if

  ! Set file specific parameters
  nPix       = n_pix_r(iFile)
  aOffset    = xOffset(iFile)
  bOffset    = yOffset(iFile)
  rSizeImage = r_size_image(iFile)
  rSizeImage2= rSizeImage**2
  rOccult    = radius_occult(iFile)
  rOccult2   = rOccult**2
  OffsetAngle= offset_angle(iFile)

  ! Rotate observation point from HGI system to the current coordinate system
  ObsPos_D    = matmul(FromHgi_DD, ObsPos_DI(:,ifile))
  ObsDistance = sqrt(sum(ObsPos_D**2))
  ! Normalize line of sight vector pointing towards the origin
  Los_D       = -ObsPos_D/ObsDistance
  ! Rotation with offset angle
  Los_D =matmul( rot_matrix_z(OffsetAngle), Los_D)
  ! Observer distance from image plane
  ObsDistance = abs(sum(ObsPos_D*Los_D))

  ! Make zero components slightly different from zero
  where(Los_D == 0.0) Los_D = cTiny

  SizePix = 2*rSizeImage/(nPix - 1)

  if(oktest .and. iProc==0) then
     write(*,*) 'ObsPos         =',ObsPos_DI(:,ifile)
     write(*,*) 'Los_D          =', Los_D
     write(*,*) 'rSizeImage     =',rSizeImage
     write(*,*) 'aOffset,bOffset=', aOffset, bOffset
     write(*,*) 'ImageCenter_D  =',ImageCenter_D
     write(*,*) 'SizePix        =',SizePix
     write(*,*) 'nPix           =',nPix
  end if

  unitstr_TEC = ''
  unitstr_IDL = ''

  plot_type1=plot_type(ifile)
  plot_vars1=plot_vars(ifile)
  plot_pars1=plot_pars(ifile)

  if(oktest_me)write(*,*)'ifile=',ifile,' plot_type=',plot_type1, &
       ' form = ',plot_form(ifile)

  call lower_case(plot_vars1)
  call split_str(plot_vars1,nPlotvarLosMax,plotvarnames,nplotvar)
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
  a_mag = sqrt(Los_D(1)**2 + Los_D(2)**2)
  a_Pix(1) =  Los_D(2)/a_mag
  a_Pix(2) = -Los_D(1)/a_mag
  a_Pix(3) =  0.0
  b_mag = sqrt((Los_D(1)*Los_D(3))**2 + (Los_D(2)*Los_D(3))**2 &
       + (Los_D(1)**2 + Los_D(2)**2)**2)
  b_Pix(1) = -Los_D(1)*Los_D(3)/b_mag
  b_Pix(2) = -Los_D(2)*Los_D(3)/b_mag
  b_Pix(3) = (Los_D(1)**2 + Los_D(2)**2)/b_mag

  ! 3D vector pointing from the Sun to the image center
  ImageCenter_D = ObsPos_D + ObsDistance*Los_D + aOffset*a_Pix + bOffset*b_Pix

  ! Make offset to be relative to the Sun (and not the projected observer)
  aOffset = dot_product(ImageCenter_D, a_Pix)
  bOffset = dot_product(ImageCenter_D, b_Pix)

!!!aOffset = aOffset + dot_product(ObsPos_D, a_Pix)

  allocate( &
       PlotVar(nPix,nPix,nplotvar), &
       PlotBLK(nPix,nPix,nplotvar), &
       los_image(nPix,nPix,nplotvar))

  PlotVar = 0.0

  UseScattering = any(plotvarnames(1:nPlotVar) == 'wl') &
       .or.       any(plotvarnames(1:nPlotVar) == 'pb')


  if(DoTiming)call timing_start('los_block_loop')

  ! loop over blocks
  do iBLK = 1,nBlock

     if (unusedBLK(iBLK)) CYCLE

     rBlockSize = 0.5*sqrt(&
          ((nI+1)*dx_BLK(iBLK))**2 + &
          ((nJ+1)*dy_BLK(iBLK))**2 + &
          ((nK+1)*dz_BLK(iBLK))**2)

     !position of the block center
     XyzBlockCenter_D(1) = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(2) = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(3) = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))
     rBlockCenter = sqrt(sum(XyzBlockCenter_D**2))

     if(rBlockCenter < rInner-rBlockSize) CYCLE

     if(rBlockCenter > rOuter+rBlockSize) CYCLE

     ! calculate position of block center on the los image
     BlockDistance = dot_product(Los_D, XyzBlockCenter_D - ObsPos_D)

     DoCheckBlock = BlockDistance > 0

     if(DoCheckBlock)then
        Ratio = ObsDistance/BlockDistance
        ! 3D vector from the image center to the projected block center
        XyzBlockCenter_D = Ratio*(XyzBlockCenter_D -  ObsPos_D) + ObsPos_D &
             - ImageCenter_D
        xLosBlock = dot_product(XyzBlockCenter_D,a_pix)
        yLosBlock = dot_product(XyzBlockCenter_D,b_pix)

        ! Project block size
        rBlockSize = rBlockSize*Ratio

        ! Check if block is inside the LOS image
        if(rSizeImage < sqrt(xLosBlock**2+yLosBlock**2) - rBlockSize) CYCLE
     end if
     ! Initialize plot variable for this block
     PlotBLK = 0.0

     ! Loop over pixels
     do jPix=1,nPix

        ! Y position of the pixel on the image plane
        y_Pix = (jPix -1) * SizePix -rSizeImage

        ! Check if block can intersect this pixel
        if(DoCheckBlock)then
           if(abs(y_Pix-yLosBlock)>rBlockSize)CYCLE
        end if

        do iPix=1,nPix

           ! X position of the pixel on the image plane
           x_Pix = (iPix - 1) * SizePix - rSizeImage

           ! Check if block can intersect this pixel
           if(DoCheckBlock)then
              if( (x_Pix-xLosBlock)**2 + (y_Pix-yLosBlock)**2 > &
                   rBlockSize**2 ) CYCLE 
           end if

           r2Pix = (x_Pix + aOffset)**2 + (y_Pix + bOffset)**2
           ! Check if pixel is within occultation radius
           if( r2Pix  <= rOccult2 ) CYCLE

           r2Pix = x_Pix**2 + y_Pix**2
           ! Check if pixel is outside the circular region
           if( r2Pix > rSizeImage2 ) CYCLE 

           ! Calculate contribution of this block to this pixel
           call set_plotvar_los

        end do ! jPix loop
     end do ! iPix loop

     ! sum over blocks on a pe, then sum over pe's
     PlotVar = PlotVar + PlotBLK
  end do       !iBLK

  if(DoTiming)call timing_stop('los_block_loop')

  if (plot_dimensional(ifile)) call dimensionalize_plotvar_los

  ! collect the pixels on one node and then write out the file 
  if(nProc>1)then
     call MPI_REDUCE(PlotVar,los_image,nPix*nPix*nplotvar, &
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
        call get_time_string
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_t"//StringDateOrTime//"_n",n_step,&
             file_extension
     else
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_n",n_step,file_extension
     end if
     open(unit_tmp,file=filename,status="replace",err=999)

     ! write header file
     select case(plot_form(ifile))
     case('tec')
        write(unit_tmp,*) 'TITLE="BATSRUS: Synthetic Image"'
        write(unit_tmp,'(a)')trim(unitstr_TEC)
        write(unit_tmp,*) 'ZONE T="LOS Image"', &
             ', I=',nPix,', J=',nPix,', K=1, F=POINT'
        ! Write point values
        do iPix=1,nPix
           x_Pix = (iPix - 1) * SizePix - rSizeImage
           do jPix=1,nPix
              y_Pix = (jPix - 1) * SizePix - rSizeImage

              if (plot_dimensional(ifile)) then
                 write(unit_tmp,fmt="(30(E14.6))") x_Pix*No2Io_V(UnitX_), &
                      y_Pix*No2Io_V(UnitX_),los_image(iPix,jPix,1:nPlotVar)
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
        write(unit_tmp,"(2i4)") nPix, nPix

        ! Equation parameters
        write(unit_tmp,"(100es13.5)") eqpar(1:neqpar)

        ! Coordinate, variable and equation parameter names
        write(unit_tmp,"(a)")allnames

        ! Data
        do jPix=1,nPix
           y_Pix = (jPix - 1) * SizePix - rSizeImage
           do iPix=1,nPix
              x_Pix = (iPix - 1) * SizePix - rSizeImage

              if (plot_dimensional(ifile)) then
                 x_Pix = x_Pix * No2Io_V(UnitX_)
                 y_Pix = y_Pix * No2Io_V(UnitX_)
              end if

              write(unit_tmp,fmt="(30es18.10)") &
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

  subroutine set_plotvar_los

    ! Local variables
    integer :: i, j, k, counter
    real :: intrsct(2,3,3), face_location(2,3)
    real :: xx1, xx2, yy1, yy2, zz1, zz2
    real :: point_1(3), point_2(3)
    real :: R2Point1, R2Point2,rLine_D(3),rLine2
    real :: coeff1,coeff2,coeff3
    real :: Discr
    real :: Solution1, Solution1_D(3), Solution2, Solution2_D(3)
    logical :: IsOuter, IsGoodSolution1, IsGoodSolution2  

    !-------------------------------------------------------------------------
    !if(DoTiming)call timing_start('los_set_plotvar')

    ! Get the 3D location of the pixel
    r_Pix = ImageCenter_D + x_Pix*a_Pix + y_Pix*b_Pix

    !x_los, y_los, z_los, r_los give the position of the point on the los
    !mu_los parameter related to the limb darkening
    !face_location give the locations of the faces of the block
    !face_location(2,3) = x1, y1, z1---x2, y2, z2 

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
    LosPix_D = ObsPos_D - r_Pix
    LosPix_D = LosPix_D/sqrt(sum(LosPix_D**2))
    where(LosPix_D ==0.0) LosPix_D = cTiny

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
                     (LosPix_D(k)/LosPix_D(j))*(face_location(i,j) - r_Pix(j))
             end if
          end do
       end do
    end do

    !which of the 6 points are on the block?
    counter = 0
    CHECK: do i=1,2 
       do j=1,3 
          if( (intrsct(i,j,1) >= xx1) .and. (intrsct(i,j,1) <= xx2)) then
             if( (intrsct(i,j,2) >= yy1) .and. (intrsct(i,j,2) <= yy2)) then
                if( (intrsct(i,j,3) >= zz1) .and. (intrsct(i,j,3) <= zz2)) then
                   counter = counter + 1
                   if(counter == 1) point_1 = intrsct(i,j,:)
                   if(counter == 2) then
                      point_2 = intrsct(i,j,:)
                      ! If point 2 is different from point 1, we are done
                      if(sum(abs(point_1 - point_2)) > cTolerance) EXIT CHECK
                      ! Ignore the second point, keep checking
                      counter = 1
                   end if
                end if
             end if
          end if
       end do
    end do CHECK

    ! Check if the los cuts through the block 
    if(counter /= 2) RETURN 

    R2Point1 = sum(point_1**2)
    R2Point2 = sum(point_2**2)

    ! Check if the whole segment is inside rInner
    if( R2Point1 <= rInner2 .and. R2Point2 <= rInner2) RETURN

    ! Check if the whole segment is outside rOuter
    rLine_D = r_Pix - LosPix_D*dot_product(LosPix_D, r_Pix)
    rLine2  = sum(rLine_D**2)
    if( rLine2 > rOuter2 ) RETURN

    ! Check if there is a need to calculate an intersection

    ! Do we intersect the outer sphere
    IsOuter = R2Point1 > rOuter2 .or. R2Point2 > rOuter2

    ! Do we intersect the inner or outer spheres
    if( IsOuter .or. &
         (rLine2 < rInner2 .and. rBlockCenter < rInner+rBlockSize) ) then

       coeff1 = sum((point_2 - point_1)**2)
       coeff2 = 2*dot_product(point_1, point_2 - point_1)

       if( IsOuter ) then
          coeff3 = R2Point1 - rOuter2
       else
          coeff3 = R2Point1 - rInner2
       end if

       Discr = coeff2**2-4*coeff1*coeff3

       if(Discr < 0.0)then
          write(*,*)'Warning: Discr=',Discr
          !   call stop_mpi("Negative discriminant")
          RETURN
       end if

       ! Line of sight tangent to the outer sphere
       if(IsOuter.AND.Discr==0.0)RETURN

       ! Find the two intersections (distance from point1 towards point2)
       Discr = sqrt(Discr)
       Solution1 = (-coeff2-Discr)/(2*coeff1)
       Solution2 = (-coeff2+Discr)/(2*coeff1)

       Solution1_D = point_1 + (point_2 - point_1) * Solution1
       Solution2_D = point_1 + (point_2 - point_1) * Solution2


       ! Check if the solutions are within the segment
       IsGoodSolution1 = (Solution1 >= 0.0 .and. Solution1 <= 1.0)
       IsGoodSolution2 = (Solution2 >= 0.0 .and. Solution2 <= 1.0)

       if(IsOuter)then
          ! For outer sphere replace
          ! outlying point1 with solution1 and
          ! outlying point2 with solution2
          if(R2Point1 > rOuter2) then
             if(IsGoodSolution1)then
                point_1 = Solution1_D
             else
                RETURN
             end if
          end if
          if(R2Point2 > rOuter2) then
             if(IsGoodSolution2)then
                point_2 = Solution2_D
             else
                RETURN
             end if
          end if
       else
          ! For inner sphere replace 
          ! internal point1 with solution2 and 
          ! internal point2 with solution1
          if(R2Point1 < rInner2) point_1 = Solution2_D
          if(R2Point2 < rInner2) point_2 = Solution1_D
          ! Weird case: the segment cuts the inner sphere
          if(IsGoodSolution1 .and. IsGoodSolution2)then
             ! Need to do two integrals:
             ! from point1 to solution1 and
             ! from point2 to solution2
             if(Discr > 0.0)then
                if(Solution1 > cTiny) &
                     call integrate_segment(point_1, Solution1_D)
                if(solution2< 1 - cTiny) &
                     call integrate_segment(point_2, Solution2_D)
                RETURN
             end if
          end if

       end if
    end if

    call integrate_segment(point_1, point_2)

  end subroutine set_plotvar_los

  !===========================================================================

  subroutine integrate_segment(Point_1, Point_2)

    real, intent(in) :: Point_1(3), Point_2(3)
    integer, parameter ::  nline_seg = nI + nJ + nK

    real :: Direction, s_los_sqrd
    integer :: i_los
    real :: x_los, y_los, z_los, r_los
    real :: x_q, y_q, z_q, q_mag
    real :: a_los, b_los, c_los, d_los
    real :: sin_omega, cos_omega, cos_theta,Sin2Omega,Cos2Omega, Logarithm
    real :: rho_los, ds_los
    real :: point_in(3), point_los(3)
    !------------------------------------------------------------------------

    Direction  = dot_product((Point_1 - Point_2), LosPix_D)
    s_los_sqrd = sum((Point_2 - Point_1)**2)

    if(direction > 0) then 
       point_in = point_2
    else 
       point_in = point_1
    endif
    ds_los = sqrt(s_los_sqrd) / nline_seg

    !if(DoTiming)call timing_start('los_integral')
    do i_los=1,nline_seg
       x_los = point_in(1) + (i_los-0.5)*ds_los*LosPix_D(1)
       y_los = point_in(2) + (i_los-0.5)*ds_los*LosPix_D(2)
       z_los = point_in(3) + (i_los-0.5)*ds_los*LosPix_D(3)
       r_los = sqrt(x_los**2 + y_los**2 + z_los**2)
       point_los(1) = x_los
       point_los(2) = y_los
       point_los(3) = z_los

       if(UseScattering .and. r_los > 1.0)then
          ! This calculation is useful for light scattering in SC and IH
          ! as it assumes that the radiation comes from a central 
          ! body with radius 1. Normally setting rOccult > 1 ensures r_los > 1.
          sin_omega = 1.0/r_los
          Sin2Omega = sin_omega**2
          Cos2Omega = 1 - Sin2Omega
          cos_omega = sqrt(Cos2Omega)
          Logarithm = log((1.0 + sin_omega)/cos_omega)  

          !omega and functions of omega are unique to a given line of sight
          a_los = cos_omega*Sin2Omega
          b_los = -0.125*( 1.0 - 3.0*Sin2Omega - (Cos2Omega/sin_omega)* &
               (1.0 + 3.0*Sin2Omega)*Logarithm )
          c_los = 4.0/3.0 - cos_omega - (1.0/3.0)*cos_omega*Cos2Omega
          d_los = 0.125*( 5.0 + sin_omega**2 - (Cos2omega/sin_omega) * &
               (5.0 - Sin2Omega)*Logarithm )

          z_q =   (LosPix_D(1)**2 + LosPix_D(2)**2)*z_los            &
               - LosPix_D(3)*(LosPix_D(1)*x_los + LosPix_D(2)*y_los)
          x_q = x_los + (LosPix_D(1)/LosPix_D(3)) * (z_q - z_los)
          y_q = y_los + (LosPix_D(2)/LosPix_D(3)) * (z_q - z_los)
          q_mag = sqrt(x_q**2 + y_q**2 + z_q**2)

          cos_theta = q_mag/r_los       
       end if

       ! interpolate density at this point with bilinear interpolation
       rho_los = point_value_los(point_los, iBLK)

       do iVar = 1, nPlotvar
          NameVar = plotvarnames(iVar)
          select case(NameVar)
          case ('len')
             ! Integrate the length of the integration lines
             PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + ds_los

             !case('vlos','Vlos','ulos','Ulos')
             !   ! Integrate the velocity
             !   PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
             !        (rhoU_los/rho_los)*ds_los

          case('wl')
             ! White light with limb darkening
             if(r_los > 1.0) &
                  PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                  rho_los*( &
                  (1.0 - mu_los)*(2.0*c_los - a_los*cos_theta**2) &
                  + mu_los*(2.0*d_los - b_los*cos_theta**2) )*ds_los

          case('pb')
             ! Polarization brightness
             if(r_los > 1.0) &
                  PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                  rho_los*( (1.0 - mu_los)*a_los + mu_los*b_los) &
                  *cos_theta**2*ds_los

          case('rho')
             ! Simple density integral
             PlotBLK(iPix,jPix,iVar) = PlotBLK(iPix,jPix,iVar) + &
                  rho_los*ds_los

          case default
             PlotBLK(iPix,jPix,iVar)=-7777.
          end select
       end do ! iVar

    end do !line segment interation loop 

  end subroutine integrate_segment

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
    x(1)=(Xpnt(1)-x_BLK(1,1,1,iBLK))/dx_BLK(iBLK)+1./2.
    x(2)=(Xpnt(2)-y_BLK(1,1,1,iBLK))/dy_BLK(iBLK)+1./2.
    x(3)=(Xpnt(3)-z_BLK(1,1,1,iBLK))/dz_BLK(iBLK)+1./2.

    ! Determine cell indices corresponding to location qx
    i1=floor(x(1)); i2=i1+1
    j1=floor(x(2)); j2=j1+1
    k1=floor(x(3)); k2=k1+1

    ! Distance relative to the cell centers
    dx1=x(1)-i1; dx2=1.-dx1
    dy1=x(2)-j1; dy2=1.-dy1
    dz1=x(3)-k1; dz2=1.-dz1

    ! Trilinear interpolation of density
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

  subroutine dimensionalize_plotvar_los

    use ModPhysics, ONLY : No2Io_V, No2Si_V, UnitX_, UnitU_, UnitRho_
    !--------------------------------------------------------------------------

    do iVar = 1, nPlotVar
       NameVar = plotvarnames(iVar)

       select case(NameVar)
       case ('len')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*No2Si_V(UnitX_)
       case('rho')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*No2Si_V(UnitRho_)*No2Si_V(UnitX_)
       case('vlos','Vlos','ulos','Ulos')
          PlotVar(:,:,iVar)=PlotVar(:,:,iVar)*No2Io_V(UnitU_)
       case default
          ! no normalization
       end select

    end do ! iVar
  end subroutine dimensionalize_plotvar_los

end subroutine write_plot_los

!==============================================================================

subroutine get_TEC_los_variables(iFile,nplotvar,plotvarnames,unitstr_TEC)

  use ModPhysics, ONLY : NameTecUnit_V, UnitX_, UnitU_
  use ModIO, ONLY: plot_dimensional
  implicit none

  ! Arguments

  integer, intent(in) :: Nplotvar,iFile
  character (LEN=10), intent(in) :: plotvarnames(Nplotvar)
  character (len=500), intent(out) :: unitstr_TEC 
  character (len=10) :: s

  integer :: iVar
  !--------------------------------------------------------------------------

  !\
  ! This routine takes the plot_var information and loads the header file with
  ! the appropriate string of variable names and units
  !/

  if (plot_dimensional(ifile)) then
     write(unitstr_TEC,'(a)') 'VARIABLES = '
     write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'"X '//&
          trim(NameTecUnit_V(UnitX_))
     write(unitstr_TEC,'(a)') trim(unitstr_TEC)//'", "Y '//&
          trim(NameTecUnit_V(UnitX_))
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
                trim(NameTecUnit_V(UnitX_))
        case('rho')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'`r [m^-^2]'
        case('vlos','Vlos','ulos','Ulos')
           write(unitstr_TEC,'(a)') & 
                trim(unitstr_TEC)//'u.s'//' '//&
                trim(NameTecUnit_V(UnitU_))
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

!==============================================================================
subroutine get_IDL_los_units(ifile,nplotvar,plotvarnames,unitstr_IDL)

  use ModPhysics, ONLY : NameIdlUnit_V, UnitX_, UnitU_
  use ModIO, ONLY : plot_dimensional

  implicit none

  ! Arguments

  integer, intent(in) :: iFile,Nplotvar
  character (LEN=10), intent(in) :: plotvarnames(Nplotvar)
  character (len=79), intent(out) :: unitstr_IDL 
  character (len=10) :: s

  integer :: iVar
  !----------------------------------------------------------------------------

  !\
  ! This routine takes the plot_var information and loads the header file with
  ! the appropriate string of unit values
  !/

  if (plot_dimensional(ifile)) then
     write(unitstr_IDL,'(a)') trim(NameIdlUnit_V(UnitX_))//' '//&
          trim(NameIdlUnit_V(UnitX_))//' '//&
          trim(NameIdlUnit_V(UnitX_))
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
                trim(NameIdlUnit_V(UnitX_))
        case('rho')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//'[m^-^2]'
        case('vlos','Vlos','ulos','Ulos')
           write(unitstr_IDL,'(a)') & 
                trim(unitstr_IDL)//' '//&
                trim(NameIdlUnit_V(UnitU_))
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
