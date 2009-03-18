!^CFG COPYRIGHT UM
!=============================================================================
subroutine write_plot_los(iFile)

  ! Purpose:  Integrate some quantities along several lines of sight and
  !           create a 2D image of the integrated quantities.
  !           The viewing point can be inifinitely far or at a finite distance.
  !           Applications include integrating density, 
  !           creating a synthetic coronagraph image of Thomson scattered 
  !           white light by integrating light scattered in the line of sight.
  !
  !           The algorithm loops over all blocks per processor (in parallel)
  !           and over lines of sight and then the results obtained on the
  !           processors are added up.
  !
  !           Written by Chip Manchester, KC Hansen
  !                   some improvements by Gabor Toth
  !                   some changes by Noe Lugaz
  !
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
  !           March 2009    Allow integration of functions different from rho

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
  use ModCoordTransform, ONLY : rot_matrix_z, cross_product
  use ModUtilities, ONLY: lower_case
  use ModPlotFile, ONLY: save_plot_file
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
  real, allocatable :: ImagePe_VII(:,:,:), Image_VII(:,:,:)

  real ::     eqpar(neqparmax)
  character (len=10) :: eqparnames(neqparmax)
  character (len=10) :: plotvarnames(nplotvarlosmax)
  character (len=10) :: NameVar


  integer :: nEqpar, nPlotvar
  integer :: iPix, jPix
  real    ::  x_Pix, y_Pix 
  real    :: ImageCenter_D(3), aPix_D(3), bPix_D(3), r_Pix(3),LosPix_D(3)
  real    :: rBlockSize, rBlockCenter
  real    :: SizePix, r2Pix    
  real    :: BlockDistance, ObsDistance, Ratio
  real    :: XyzBlockCenter_D(3), CellSize_D(3), xLosBlock, yLosBlock

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
  logical :: UseScattering, UseRho

  character(len=*), parameter :: NameSub = 'write_plot_los'
  !---------------------------------------------------------------------------

  ! Initialize stuff
  call set_oktest(NameSub, oktest, oktest_me)
  call set_oktest('los_timing', DoTiming, DoTimingMe)

  call timing_start(NameSub)

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
     rOuter = 1e30
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
  ObsPos_D    = matmul(FromHgi_DD, ObsPos_DI(:,iFile))
  ObsDistance = sqrt(sum(ObsPos_D**2))
  ! Normalize line of sight vector pointing towards the origin
  Los_D       = -ObsPos_D/ObsDistance
  ! Rotation with offset angle
  Los_D =matmul( rot_matrix_z(OffsetAngle), Los_D)
  ! Observer distance from image plane
  ObsDistance = abs(sum(ObsPos_D*Los_D))

  ! Make zero components slightly different from zero
  where(Los_D == 0.0) Los_D = cTiny

  ! Pixel size for node based pixel grid
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
     write(*,*) 'plot variables, UseRho=', plot_vars1, UseRho
     write(*,*) 'nPlotVar, PlotVarNames_V=', &
          nPlotVar,plotvarnames(1:nplotvar)
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

  ! Create unit vectors aPix_D and bPix_D orthogonal to the 
  ! central line of sight to setup the coordinate system in the viewing plane
  if(abs(Los_D(3)) < maxval(abs(Los_D(1:2))))then
     aPix_D = cross_product(Los_D, (/0.,0.,1./))
  else
     ! Viewing along the Z axis more or less
     aPix_D = cross_product(Los_D, (/0.,1.,0./))
  end if
  aPix_D = aPix_D/sqrt(sum(aPix_D**2))
  bPix_D = cross_product(Los_D, aPix_D)
  bPix_D = bPix_D/sqrt(sum(bPix_D**2))

  ! 3D vector pointing from the origin to the image center
  ImageCenter_D = ObsPos_D + ObsDistance*Los_D &
       + aOffset*aPix_D + bOffset*bPix_D

  ! Make offset to be relative to the Sun (and not the projected observer)
  aOffset = dot_product(ImageCenter_D, aPix_D)
  bOffset = dot_product(ImageCenter_D, bPix_D)

!!!aOffset = aOffset + dot_product(ObsPos_D, aPix_D)

  allocate( &
       ImagePe_VII(nPlotVar,nPix,nPix), &
       Image_VII(nPlotVar,nPix,nPix))

  ImagePe_VII = 0.0

  ! Do we need to apply scattering
  UseScattering = any(plotvarnames(1:nPlotVar) == 'wl') &
       .or.       any(plotvarnames(1:nPlotVar) == 'pb')

  ! Do we need to calculate density (also for white light and polarization)
  UseRho = UseScattering .or. any(plotvarnames(1:nPlotVar) == 'rho')

  if(DoTiming)call timing_start('los_block_loop')

  ! loop over blocks
  do iBLK = 1, nBlock

     if (unusedBLK(iBLK)) CYCLE

     CellSize_D = (/ dx_BLK(iBlk), dy_BLK(iBlk), dz_BLK(iBlk) /)

     rBlockSize = 0.5*sqrt(&
          ((nI+1)*dx_BLK(iBLK))**2 + &
          ((nJ+1)*dy_BLK(iBLK))**2 + &
          ((nK+1)*dz_BLK(iBLK))**2)

     !position of the block center
     XyzBlockCenter_D(1) = 0.50*(x_BLK(nI,nJ,nK,iBLK)+x_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(2) = 0.50*(y_BLK(nI,nJ,nK,iBLK)+y_BLK(1,1,1,iBLK))
     XyzBlockCenter_D(3) = 0.50*(z_BLK(nI,nJ,nK,iBLK)+z_BLK(1,1,1,iBLK))
     rBlockCenter = sqrt(sum(XyzBlockCenter_D**2))

     if(rBlockCenter < rInner - rBlockSize) CYCLE

     if(rBlockCenter > rOuter + rBlockSize) CYCLE

     ! calculate position of block center on the los image
     BlockDistance = dot_product(Los_D, XyzBlockCenter_D - ObsPos_D)

     DoCheckBlock = BlockDistance > 0

     if(DoCheckBlock)then
        Ratio = ObsDistance/BlockDistance
        ! 3D vector from the image center to the projected block center
        XyzBlockCenter_D = Ratio*(XyzBlockCenter_D -  ObsPos_D) + ObsPos_D &
             - ImageCenter_D
        xLosBlock = dot_product(XyzBlockCenter_D,aPix_D)
        yLosBlock = dot_product(XyzBlockCenter_D,bPix_D)

        ! Project block size
        rBlockSize = rBlockSize*Ratio

        ! Check if block is inside the LOS image
        if(rSizeImage < sqrt(xLosBlock**2+yLosBlock**2) - rBlockSize) CYCLE
     end if

     ! Loop over pixels
     do jPix = 1, nPix

        ! Y position of the pixel on the image plane
        y_Pix = (jPix -1) * SizePix - rSizeImage

        ! Check if block can intersect this pixel
        if(DoCheckBlock)then
           if(abs(y_Pix - yLosBlock) > rBlockSize)CYCLE
        end if

        do iPix = 1, nPix

           ! X position of the pixel on the image plane
           x_Pix = (iPix - 1) * SizePix - rSizeImage

           ! Check if block can intersect this pixel
           if(DoCheckBlock)then
              if( (x_Pix - xLosBlock)**2 + (y_Pix-yLosBlock)**2 > &
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
  end do       !iBLK

  if(DoTiming)call timing_stop('los_block_loop')

  ! add up the pixels from all PE-s to the root node and save LOS file 
  if(nProc > 1)then
     call MPI_REDUCE(ImagePe_VII, Image_VII, nPix*nPix*nPlotvar, &
          MPI_REAL, MPI_SUM, 0, iComm, iError)
  else
     Image_VII = ImagePe_VII
  end if

  if(iProc==0 .and. plot_dimensional(iFile)) call dimensionalize_plotvar_los

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
             ifile-plot_,"_t"//trim(StringDateOrTime)//"_n",n_step,&
             file_extension
     else
        write(filename,file_format) &
             trim(plot_type1)//"_",&
             ifile-plot_,"_n",n_step,file_extension
     end if

     ! write header file
     select case(plot_form(ifile))
     case('tec')
        open(unit_tmp,file=filename,status="replace",IOSTAT = iError)
        if(iError /= 0)call stop_mpi(NameSub//" ERROR opening "//filename)

        write(unit_tmp,*) 'TITLE="BATSRUS: Synthetic Image"'
        write(unit_tmp,'(a)')trim(unitstr_TEC)
        write(unit_tmp,*) 'ZONE T="LOS Image"', &
             ', I=',nPix,', J=',nPix,', K=1, F=POINT'
        ! Write point values
        do iPix = 1, nPix
           x_Pix = (iPix - 1) * SizePix - rSizeImage
           do jPix = 1, nPix
              y_Pix = (jPix - 1) * SizePix - rSizeImage

              if (plot_dimensional(ifile)) then
                 write(unit_tmp,fmt="(30(E14.6))") x_Pix*No2Io_V(UnitX_), &
                      y_Pix*No2Io_V(UnitX_), Image_VII(1:nPlotVar,iPix,jPix)
              else
                 write(unit_tmp,fmt="(30(E14.6))") x_Pix, y_Pix, &
                      Image_VII(1:nPlotVar,iPix,jPix)
              end if

           end do
        end do
        close(unit_tmp)
     case('idl')
        ! description of file contains units, physics and dimension
        StringHeadLine = 'LOS integrals_var22'
        ! set the size of plot image
        x_Pix = rSizeImage 
        if (plot_dimensional(ifile)) x_Pix = x_pix * No2Io_V(UnitX_)
        ! Write
        call save_plot_file(filename, &
             TypeFileIn = TypeIdlFile_I(iFile), &
             StringHeaderIn = StringHeadLine, &
             nStepIn = n_step, &
             TimeIn = time_simulation, &
             ParamIn_I = eqpar(1:neqpar), &
             NameVarIn = allnames, &
             nDimIn = 2, & 
             CoordMinIn_D = (/-x_Pix, -x_Pix/), &
             CoordMaxIn_D = (/+x_Pix, +x_Pix/), &
             VarIn_VII = Image_VII)
        
     end select
     
  end if  !iProc ==0
  if(DoTiming)call timing_stop('los_save_plot')

  deallocate(ImagePe_VII, Image_VII)

  if(oktest_me)write(*,*) NameSub,' finished'

  call timing_stop(NameSub)

contains

  !===========================================================================

  subroutine set_plotvar_los

    ! Local variables
    integer :: i, j, k, counter
    real :: intrsct(2,3,3), face_location(2,3)
    real :: xx1, xx2, yy1, yy2, zz1, zz2
    real :: Point1_D(3), Point2_D(3)
    real :: R2Point1, R2Point2,rLine_D(3),rLine2
    real :: coeff1,coeff2,coeff3
    real :: Discr
    real :: Solution1, Solution1_D(3), Solution2, Solution2_D(3)
    logical :: IsOuter, IsGoodSolution1, IsGoodSolution2  

    !-------------------------------------------------------------------------
    !if(DoTiming)call timing_start('los_set_plotvar')

    ! Get the 3D location of the pixel
    r_Pix = ImageCenter_D + x_Pix*aPix_D + y_Pix*bPix_D

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
                   if(counter == 1) Point1_D = intrsct(i,j,:)
                   if(counter == 2) then
                      Point2_D = intrsct(i,j,:)
                      ! If point 2 is different from point 1, we are done
                      if(sum(abs(Point1_D - Point2_D)) > cTolerance) EXIT CHECK
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

    R2Point1 = sum(Point1_D**2)
    R2Point2 = sum(Point2_D**2)

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

       coeff1 = sum((Point2_D - Point1_D)**2)
       coeff2 = 2*dot_product(Point1_D, Point2_D - Point1_D)

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

       Solution1_D = Point1_D + (Point2_D - Point1_D) * Solution1
       Solution2_D = Point1_D + (Point2_D - Point1_D) * Solution2


       ! Check if the solutions are within the segment
       IsGoodSolution1 = (Solution1 >= 0.0 .and. Solution1 <= 1.0)
       IsGoodSolution2 = (Solution2 >= 0.0 .and. Solution2 <= 1.0)

       if(IsOuter)then
          ! For outer sphere replace
          ! outlying point1 with solution1 and
          ! outlying point2 with solution2
          if(R2Point1 > rOuter2) then
             if(IsGoodSolution1)then
                Point1_D = Solution1_D
             else
                RETURN
             end if
          end if
          if(R2Point2 > rOuter2) then
             if(IsGoodSolution2)then
                Point2_D = Solution2_D
             else
                RETURN
             end if
          end if
       else
          ! For inner sphere replace 
          ! internal point1 with solution2 and 
          ! internal point2 with solution1
          if(R2Point1 < rInner2) Point1_D = Solution2_D
          if(R2Point2 < rInner2) Point2_D = Solution1_D
          ! Weird case: the segment cuts the inner sphere
          if(IsGoodSolution1 .and. IsGoodSolution2)then
             ! Need to do two integrals:
             ! from point1 to solution1 and
             ! from point2 to solution2
             if(Discr > 0.0)then
                if(Solution1 > cTiny) &
                     call integrate_segment(Point1_D, Solution1_D)
                if(solution2< 1 - cTiny) &
                     call integrate_segment(Point2_D, Solution2_D)
                RETURN
             end if
          end if

       end if
    end if

    call integrate_segment(Point1_D, Point2_D)

  end subroutine set_plotvar_los

  !===========================================================================

  subroutine integrate_segment(XyzStart_D, XyzEnd_D)

    ! Integrate variables from XyzStart_D to XyzEnd_D
    ! The line is split into nSegment segments of length Ds

    use ModGeometry,    ONLY: XyzStart_BLK
    use ModInterpolate, ONLY: trilinear
    use ModUser,        ONLY: user_set_plot_var

    real, intent(in) :: XyzStart_D(3), XyzEnd_D(3)
    integer, parameter ::  nSegment = nI + nJ + nK

    integer :: iSegment
    real :: x_q, y_q, z_q, q_mag
    real :: a_los, b_los, c_los, d_los
    real :: SinOmega, CosOmega, Cos2Theta, Sin2Omega, Cos2Omega, Logarithm

    real :: Ds             ! Length of line segment
    real :: XyzLos_D(3)    ! Coordinate of center of line segment    
    real :: xLos, yLos, zLos, rLos ! Coordinates and radial distance
    real :: CoordNorm_D(3) ! Normalized coordinates of current point
    real :: Rho            ! Density at the point
    real :: Value          ! Value of the LOS variable at the point

    ! Variables for user defined LOS variables
    integer :: iBlockLast = -1, iVarLast = -1
    logical :: IsFound, UseBody
    character(len=1):: NameTecVar, NameTecUnit, NameIdlUnit
    real    :: ValueBody
    real, allocatable, save:: PlotVar_GV(:,:,:,:)
    !------------------------------------------------------------------------

    !if(DoTiming)call timing_start('los_integral')

    ! Length of a segment
    Ds = sqrt(sum((XyzEnd_D - XyzStart_D)**2)) / nSegment

    do iSegment = 1, nSegment
       XyzLos_D = XyzStart_D &
            + (iSegment - 0.5)/nSegment*(XyzEnd_D - XyzStart_D)
       rLos = sqrt(sum(XyzLos_D**2))
       xLos = XyzLos_D(1)
       yLos = XyzLos_D(2)
       zLos = XyzLos_D(3)

       if(UseScattering .and. rLos > 1.0)then
          ! This calculation is useful for light scattering in SC and IH
          ! as it assumes that the radiation comes from a central 
          ! body with radius 1. Normally setting rOccult > 1 ensures rLos > 1.
          SinOmega = 1.0/rLos
          Sin2Omega = SinOmega**2
          Cos2Omega = 1 - Sin2Omega
          CosOmega = sqrt(Cos2Omega)
          Logarithm = log((1.0 + SinOmega)/CosOmega)  

          !omega and functions of omega are unique to a given line of sight
          a_los = CosOmega*Sin2Omega
          b_los = -0.125*( 1.0 - 3.0*Sin2Omega - (Cos2Omega/SinOmega)* &
               (1.0 + 3.0*Sin2Omega)*Logarithm )
          c_los = 4.0/3.0 - CosOmega - (1.0/3.0)*CosOmega*Cos2Omega
          d_los = 0.125*( 5.0 + Sin2Omega - (Cos2omega/SinOmega) * &
               (5.0 - Sin2Omega)*Logarithm )

          z_q =   (LosPix_D(1)**2 + LosPix_D(2)**2)*zLos            &
               - LosPix_D(3)*(LosPix_D(1)*xLos + LosPix_D(2)*yLos)
          x_q = xLos + (LosPix_D(1)/LosPix_D(3)) * (z_q - zLos)
          y_q = yLos + (LosPix_D(2)/LosPix_D(3)) * (z_q - zLos)
          q_mag = sqrt(x_q**2 + y_q**2 + z_q**2)

          Cos2Theta = (q_mag/rLos)**2
       end if

       ! Calculate normalized position
       ! XyzStart contains the coordinates of cell 1,1,1, hence add 1
       CoordNorm_D = (XyzLos_D - XyzStart_BLK(:,iBlk))/CellSize_D + 1

       ! interpolate density if it is needed by any of the plot variables
       if(UseRho) Rho = trilinear(State_VGB(Rho_,:,:,:,iBlk), &
            -1, nI+2, -1, nJ+2, -1, nK+2, CoordNorm_D)

       do iVar = 1, nPlotvar
          Value = 0.0 ! initialize to 0 so that if statements below work right
          NameVar = plotvarnames(iVar)
          select case(NameVar)
          case ('len')
             ! Integrate the length of the integration lines
             Value = 1.0

          case('wl')
             ! White light with limb darkening
             if(rLos > 1.0) Value = Rho*( &
                  (1 - mu_los)*(2*c_los - a_los*Cos2Theta) &
                  + mu_los*(2*d_los - b_los*Cos2Theta) )

          case('pb')
             ! Polarization brightness
             if(rLos > 1.0) Value = &
                  Rho*( (1.0 - mu_los)*a_los + mu_los*b_los)*Cos2Theta

          case('rho')
             ! Simple density integral
             Value = Rho

          case('sphere10')
             Value = max(0.0, 100.0 - rLos**2)
             
          case('cube10')
             Value = product( 0.5 + sign(0.5, 10.0 - abs(XyzLos_D)) )

          case default
             ! Obtain user defined plot function for the whole block
             if(iBlk /= iBlockLast .or. iVar > iVarLast)then
                iBlockLast = iBlk
                iVarLast   = iVar
                if(.not.allocated(PlotVar_GV)) &
                     allocate(PlotVar_GV(-1:nI+2,-1:nJ+2,-1:nK+2,nPlotVar))
                call user_set_plot_var(iBlk, NameVar, &
                     plot_dimensional(iFile), &
                     PlotVar_GV(:,:,:,iVar), &
                     ValueBody, UseBody, NameTecVar, NameTecUnit, NameIdlUnit,&
                     IsFound)
                if(.not. IsFound)then
                   PlotVar_GV(:,:,:,iVar)=-7777.
                   if(iProc==0.and.iBLK==1)write(*,*) &
                        NameSub, ' WARNING: unknown plotvarname=', NameVar
                end if
             end if
             ! Interpolate value
             Value = trilinear(PlotVar_GV(:,:,:,iVar), &
                  -1, nI+2, -1, nJ+2, -1, nK+2, CoordNorm_D)
          end select

          ImagePe_VII(iVar,iPix,jPix) = ImagePe_VII(iVar,iPix,jPix) + Value*Ds

       end do ! iVar

    end do !line segment interation loop 

  end subroutine integrate_segment

  !==========================================================================

  subroutine dimensionalize_plotvar_los

    use ModPhysics, ONLY : No2Io_V, No2Si_V, UnitX_, UnitRho_
    !--------------------------------------------------------------------------

    do iVar = 1, nPlotVar
       NameVar = plotvarnames(iVar)

       select case(NameVar)
       case ('len')
          Image_VII(iVar,:,:) = Image_VII(iVar,:,:)*No2Si_V(UnitX_)
       case('rho')
          Image_VII(iVar,:,:) = Image_VII(iVar,:,:) &
               *No2Si_V(UnitRho_)*No2Si_V(UnitX_)
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
