module ModRadioWaveImage

  use ModNumConst
  use ModConst
  use ModCoordTransform
  use ModRadioWaveRaytracing
  use CON_global_vector, ONLY: allocate_vector, associate_with_global_vector

  implicit none


  public :: ray_bunch_intensity

contains !=========================================================

  subroutine ray_bunch_intensity(XyzObserver_D, RadioFrequency, ImageRange_I, &
       rIntegration, nXPixel, nYPixel, Intensity_II)
    !
    !   The subroutine ray_bunch_intensity() is an interface to the radiowave
    ! raytracer and emissivity integrator ray_path() from the 
    ! ModRadioWaveRaytracing module. It calculates the initial slopes of the
    ! rays pointing at the pixels of the image plane containing the sun at 
    ! its center. It calculates an appropriate radius of the integration
    ! sphere. It then calculates the points at which the (straight) rays 
    ! penetrate the integration sphere as well as their slopes. At this point 
    ! it makes a series of calls to the ray_path() routine that calculates
    ! both ray trajectories and emissivity integrals for each ray.
    !   The ray refraction and plasma emissivity calculation is based solely
    ! on the plasma density. The density is obtained from the MHD simulation.
    ! The density data are located in the memories of several processors, and
    ! the get_plasma_density() program from the ModDensityAndGradient module
    ! makes the density acquisition. 
    !
    ! The ray_bunch_intensity() parameters are
    ! as follows:
    ! 
    ! XyzObserver_D:   position of observer; actually it is the location 
    !     of the radiotelescope.
    ! RadioFrequency:  frequency in Hertz at wich the radiotelescope works.
    ! ImageRange_I:    two pairs of coordinates of the "lower left" and "upper
    !     right" corners of the rectangular image plane. The image plane 
    !     center coincides with the sun center (usually, the numeric domain 
    !     origin. 
    ! rIntegration:    radius of the spherical surface around the sun within
    !     which the emissivity integration is performed. The plasma 
    !     emissivity is considered insignificant beyond the integration sphere.
    !     Also, the plasma density outside the integration sphere is small
    !     enough to assume no refraction. Therefore the radiowave rays from
    !     the radiotelescope to the integration sphere surface are treated as
    !     straight lines.
    ! nXPixel, nYPixel: the image resolution in horizontal and vertical pixel
    !     numbers. 
    ! Intensity_II: radiotelescope image of the size nXPixel by nYPixel. Each
    !     pixel presents the radiation intensity at the specified frequency. 
    !     The intencity of an individual pixel is calculated as the integral
    !     of the plasma emissivity along the corresponding ray trajectory. 
    ! Written by Leonid Benkevitch.
    !
    !

    use ModMain,ONLY:NameThisComp
    use ModDensityAndGradient,ONLY: get_plasma_density, NameVector, &
         Density_I, GradDensity_DI, DeltaSNew_I

    implicit none

    real, intent(in) :: XyzObserver_D(3)          
    ! Observer's position
    real, intent(in) :: RadioFrequency
    real, intent(in) :: ImageRange_I(4)        
    ! (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
    real, intent(in) :: rIntegration           
    ! Radius of "integration sphere"
    integer, intent(in) ::  nXPixel, nYPixel   
    ! Dimensions of the raster in pixels
    real, dimension(nYPixel,nXPixel), intent(out) :: Intensity_II          
    ! The result of emissivity integration
    logical, dimension(nYPixel*nXPixel) :: ExcludeRay_I   
    ! A ray is excluded from processing if it is .true.
    real, dimension(nYPixel,nXPixel) :: XPixel_II, YPixel_II   
    ! Pixel coordinates INSIDE of the image plane
    real :: Normal_D(3)                                          
    ! Unity vector normal to image plane
    real :: ZAxisOrt_D(3) = (/0, 0, 1/)
    real :: Tau_D(3), Xi_D(3)                                   
    ! Image plane inner Cartesian orts
    real :: XyzObservLen
    real, dimension(3,nYPixel,nXPixel) :: XyzPixel_DII          
    ! SGI pixel coordinates
    real, dimension(3,nYPixel,nXPixel) :: Slope_DII  
    ! SGI unity slope vectors for all the line-of-sights pointing at
    ! the pixels 
    real, dimension(nYPixel,nXPixel) :: ObservToIntSphereDist_II  
    ! Distance from the radiotelescope to the integration sphere
    real, dimension(nYPixel,nXPixel) :: ObservToIntSphere2_II  
    ! Squared distance from the radiotelescope to the integration sphere
    real, dimension(3,nYPixel,nXPixel) :: Position_DII        ! 
    real, dimension(nYPixel,nXPixel) :: XPosition_II, YPosition_II
    real, dimension(nYPixel,nXPixel) :: ZPosition_II, SolarDistSqr_II
    real, dimension(nYPixel*nXPixel) :: SolarDistSqr_I
    real :: XPixelDel, YPixelDel    
    real :: SlopeUnscaled_D(3)
    real :: XLower, XUpper, YLower, YUpper
    real, dimension(3,nXPixel*nYPixel) :: Slope_DI
    real, pointer, dimension(:,:) :: Position_DI
    real, dimension(nXPixel*nYPixel) :: Intensity_I, RayPath_I,DeltaS_I
    logical, dimension(nXPixel*nYPixel) :: RayFlag_I         
    !.true. if a ray is OK; .false. otherwise
    logical, save :: NewEntry = .true.
    ! Must be set to .true. before a call to ray_path() with new value of nRay
    real :: OneAU = 215.0, Tolerance = 0.01, DeltaS = 1.0
    real :: MaxRayPath = 60.
    real :: DensityCr 
    real :: XPixel, YPixel, SolarDistMin, MinRayPath
    real ::  PercentRayLeft, rIntegrationSqr
    integer :: nRay, nIteration, i, j, iRay, nRayInsideIntSphere
    integer, dimension(nXPixel*nYPixel) :: RayInsideIntSphere_I
    logical :: deb = .false.
    real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe
    logical, save :: DoAllocate = .true.

    !
    ! Total number of rays and pixels in the raster
    !
    nRay = nXPixel*nYPixel  

    !
    ! Calculate the critical density from the frequency
    !
    DensityCr = cPi*cProtonMass*cElectronMass*1e6 &
         *(RadioFrequency/ProtonChargeSGSe)**2

    !
    ! Determine the image plane inner coordinates of pixel centers
    !
    XLower = ImageRange_I(1)
    YLower = ImageRange_I(2)
    XUpper = ImageRange_I(3)
    YUpper = ImageRange_I(4)
    XPixelDel = (XUpper - XLower)/nXPixel
    YPixelDel = (YUpper - YLower)/nYPixel
    XPixel_II(1,:) = (/ (XLower + (real(j)-0.5)*XPixelDel, j = 1, nXPixel) /)
    do i = 2, nXPixel
       XPixel_II(i,:) = XPixel_II(1,:)
    end do
    YPixel_II(:,1) = (/ (YLower + (real(j)-0.5)*YPixelDel, j = 1, nYPixel) /)
    do i = 2, nYPixel
       YPixel_II(:,i) = YPixel_II(:,1)
    end do

    !
    ! Determune the orts, Tau and Xi, of the inner coordinate system
    ! of the image plane
    !
    XyzObservLen = sqrt(sum(XyzObserver_D**2))
    Normal_D = XyzObserver_D/XyzObservLen
    Tau_D = cross_product(ZAxisOrt_D, Normal_D)
    Xi_D = cross_product(Normal_D, Tau_D)

    !
    ! Calculate coordinates of all the pixels in the SGI
    !
    do i = 1, nYPixel
       do j = 1, nXPixel
          XyzPixel_DII(:,i,j) = XPixel_II(i,j)*Tau_D + YPixel_II(i,j)*Xi_D
          SlopeUnscaled_D = XyzPixel_DII(:,i,j) - XyzObserver_D
          Slope_DII(:,i,j) = SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2)) ! v
       end do
    end do

    !
    ! Find the points on the integration sphere where it intersects
    ! with the straight "rays" 
    !
    do i = 1, nYPixel
       do j = 1, nXPixel
          ObservToIntSphereDist_II(i,j) = -sum(XyzObserver_D*Slope_DII(:,i,j))&
               - sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) &
               - sum(cross_product(Slope_DII(:,i,j), XyzObserver_D)**2))
          ObservToIntSphere2_II(i,j) = -sum(XyzObserver_D*Slope_DII(:,i,j)) + &
               sqrt(sum((Slope_DII(:,i,j)*rIntegration)**2) &
               - sum(cross_product(Slope_DII(:,i,j), XyzObserver_D)**2))
          Position_DII(:,i,j) = XyzObserver_D &
               + Slope_DII(:,i,j)*ObservToIntSphereDist_II(i,j)
       end do
    end do

    !
    ! Initial allocation of the global vector of ray positions
    ! and other dynamic arrays
    !
    if(DoAllocate)then
       DoAllocate=.false.
       ! Symbolic name of the global vector:
       NameVector = NameThisComp//'_Pos_DI'
       call allocate_vector(NameVector, 3, nRay)
       call associate_with_global_vector(Position_DI, NameVector)
       allocate(Density_I(nRay), GradDensity_DI(3,nRay),DeltaSNew_I(nRay))
    end if

    !
    ! Do emissivity integration inside of the integration sphere 
    !
    Position_DI = reshape(Position_DII, (/3, nRay/))
    Slope_DI = reshape(Slope_DII, (/3, nRay/))
    Intensity_I = 0
    RayFlag_I = .false.
    !NewEntry = .true.
    ExcludeRay_I = .false.
    RayInsideIntSphere_I = 1
    DeltaS_I = DeltaS
    RayPath_I = 0
    nIteration = 0
    rIntegrationSqr = rIntegration**2 + 0.01

    MinRayPath = 0.0
    nRayInsideIntSphere = nRay

    do while(nRayInsideIntSphere .gt. 0)
       nIteration = nIteration + 1

       !write(*,*) 'ray_bunch_intensity(): nIteration = ', nIteration

       SolarDistSqr_I = sum(Position_DI**2,1)
       where(SolarDistSqr_I .gt. rIntegrationSqr) 
          RayInsideIntSphere_I = 0
          ExcludeRay_I = .true.
       end where

       nRayInsideIntSphere = sum(RayInsideIntSphere_I)

       call ray_path(get_plasma_density, nRay, ExcludeRay_I, Slope_DI, &
            DeltaS_I, Tolerance, DensityCr, Intensity_I, RayFlag_I)
       RayPath_I = RayPath_I + DeltaS_I
       MinRayPath = minval(RayPath_I)

    end do

    Intensity_II = reshape(Intensity_I, (/nYPixel,nXPixel/))

  end subroutine ray_bunch_intensity

end module ModRadioWaveImage

!================================================

subroutine get_ray_bunch_intensity(XyzObserver_D, RadioFrequency, &
     ImageRange_I, rIntegration, nXPixel, nYPixel, Intensity_II)

  use ModRadioWaveImage, ONLY: ray_bunch_intensity

  implicit none

  real, intent(in) :: XyzObserver_D(3)          
  ! Observer position
  real, intent(in) :: RadioFrequency  
  real, intent(in) :: ImageRange_I(4)         
  ! (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
  real, intent(in) :: rIntegration           
  ! Radius of "integration sphere"
  integer, intent(in) ::  nXPixel, nYPixel   
  ! Dimensions of the raster in pixels
  real, dimension(nYPixel,nXPixel), intent(out) :: Intensity_II          
  ! The result from emissivity integration

  call ray_bunch_intensity(XyzObserver_D, RadioFrequency, ImageRange_I, &
       rIntegration, nXPixel, nYPixel, Intensity_II)

end subroutine get_ray_bunch_intensity
