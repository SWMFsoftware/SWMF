!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
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
    !INPUTS
    ! Observer's position
    real, intent(in) :: XyzObserver_D(3)          
 
    real, intent(in) :: RadioFrequency
    ! (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
    real, intent(in) :: ImageRange_I(4)        
    ! Radius of "integration sphere"
    real, intent(in) :: rIntegration           
    ! Dimensions of the raster in pixels
    integer, intent(in) ::  nXPixel, nYPixel   
    ! The result of emissivity integration
    real, dimension(nXPixel,nYPixel), intent(out) :: Intensity_II          
    ! A ray is excluded from processing if it is .true.
    logical, dimension(nYPixel*nXPixel) :: ExcludeRay_I   
    ! Pixel coordinates INSIDE of the image plane
    real, dimension(nXPixel,nYPixel) :: XPixel_II, YPixel_II  
    ! Unity vector normal to image plane 
    real :: Normal_D(3)                                          
    real :: ZAxisOrt_D(3) = (/0, 0, 1/)
    ! Image plane inner Cartesian orts
    real :: Tau_D(3), Xi_D(3)                                   
    real :: XyzObservLen
    ! SGI pixel coordinates
    real, dimension(3,nYPixel,nXPixel) :: XyzPixel_DII
    ! SGI unity slope vectors for all the line-of-sights pointing at 
    ! the pixels          
    real, dimension(3,nYPixel,nXPixel) :: Slope_DII
    ! Distance from the radiotelescope to the integration sphere  
    real :: ObservToIntSphereDist
   
    real, dimension(3,nXPixel,nYPixel) :: Position_DII        ! 
    real, dimension(nXPixel*nYPixel) :: SolarDistSqr_I
    real :: XPixelDel, YPixelDel    
    real :: SlopeUnscaled_D(3)
    real :: XLower, XUpper, YLower, YUpper
    real, dimension(3,nXPixel*nYPixel) :: Slope_DI
    real, pointer, dimension(:,:) :: Position_DI
    real, dimension(nXPixel*nYPixel) :: Intensity_I, DeltaS_I
    !.true. if a ray is OK; .false. otherwise
    ! Must be set to .true. before a call to ray_path() with new value of nRay
    logical, dimension(nXPixel*nYPixel) :: RayFlag_I         
    real,parameter :: Tolerance = 0.01, DeltaS = 1.0
    real :: DensityCr 
    real :: rIntegrationSqr
    integer :: nRay, nIteration, i, j, nRayInsideIntSphere
    logical :: IsRayInsideIntSphere_I(nXPixel*nYPixel)
    real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe
    logical, save :: DoAllocate = .true.
    !------------------------------------

    !
    ! Total number of rays and pixels in the raster
    !
    nRay = nXPixel*nYPixel  

    !
    ! Calculate the critical density from the frequency, in CGS
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
    XPixel_II(:,1) = (/ (XLower + (real(i) - 0.5)*XPixelDel, i = 1, nXPixel) /)
    do j = 2, nYPixel
       XPixel_II(:,j) = XPixel_II(:,1)
    end do
    YPixel_II(1,:) = (/ (YLower + (real(j) - 0.5)*YPixelDel, j = 1, nYPixel) /)
    do i = 2, nXPixel
       YPixel_II(i,:) = YPixel_II(1,:)
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
    do j = 1, nYPixel
       do i = 1, nXPixel
          XyzPixel_DII(:,i,j) = XPixel_II(i,j)*Tau_D + YPixel_II(i,j)*Xi_D
          SlopeUnscaled_D = XyzPixel_DII(:,i,j) - XyzObserver_D
          Slope_DII(:,i,j) = SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2)) ! v
       end do
    end do

    !
    ! Find the points on the integration sphere where it intersects
    ! with the straight "rays" 
    !
    do j = 1, nYPixel
       do i = 1, nXPixel
          !\
          ! Solve a duadratic equation,
          !|| XyzObs_D + ObservToIntSphereDist*Slope_DI || = rIntegration
          !or  ObservToIntSphereDist**2 + 
          !  + 2*ObservToIntSphereDist*XyzObs_D
          !  + XyzObservLen**2 - rIntegration**2 = 0
          !/  
          ObservToIntSphereDist = -sum(XyzObserver_D*Slope_DII(:,i,j))&
               - sqrt(rIntegration**2 - XyzObservLen**2 &
               + sum(Slope_DII(:,i,j)*XyzObserver_D)**2)
          Position_DII(:,i,j) = XyzObserver_D &
               + Slope_DII(:,i,j)*ObservToIntSphereDist
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
    ExcludeRay_I = .false.
    IsRayInsideIntSphere_I = .true.
    DeltaS_I = DeltaS
    rIntegrationSqr = rIntegration**2 + 0.01
    nRayInsideIntSphere = nRay

    do while(nRayInsideIntSphere .gt. 0)
       SolarDistSqr_I = sum(Position_DI**2,1)
       ExcludeRay_I = SolarDistSqr_I .gt. rIntegrationSqr 
          IsRayInsideIntSphere_I = .not.ExcludeRay_I

       nRayInsideIntSphere = count(IsRayInsideIntSphere_I)

       call ray_path(get_plasma_density, nRay, ExcludeRay_I, Slope_DI, &
            DeltaS_I, Tolerance, DensityCr, Intensity_I, RayFlag_I)
    end do

    Intensity_II = reshape(Intensity_I, (/nXPixel,nYPixel/))

  end subroutine ray_bunch_intensity

end module ModRadioWaveImage

!================================================

subroutine get_ray_bunch_intensity(XyzObserver_D, RadioFrequency, &
     ImageRange_I, rIntegration, nXPixel, nYPixel, Intensity_II)

  use ModRadioWaveImage, ONLY: ray_bunch_intensity

  implicit none
  ! Observer position
  real, intent(in) :: XyzObserver_D(3)          
  real, intent(in) :: RadioFrequency  
  ! (x0, y0, x1, y1), i.e. (XLower, YLower, XUpper, YUpper)
  real, intent(in) :: ImageRange_I(4)         
  ! Radius of "integration sphere"
  real, intent(in) :: rIntegration           
 ! Dimensions of the raster in pixels
  integer, intent(in) ::  nXPixel, nYPixel  
  ! The result from emissivity integration 
  real, dimension(nYPixel,nXPixel), intent(out) :: Intensity_II          
  !--------------------------------
  call ray_bunch_intensity(XyzObserver_D, RadioFrequency, ImageRange_I, &
       rIntegration, nXPixel, nYPixel, Intensity_II)
end subroutine get_ray_bunch_intensity
