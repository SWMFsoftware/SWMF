!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModRadioWaveImage
  use ModConst, ONLY: cPi, cElectronMass, cProtonMass
  use ModCoordTransform, ONLY: cross_product
  use CON_global_vector, ONLY: allocate_vector, associate_with_global_vector
  use BATL_size, ONLY    : MaxDim
  implicit none
  SAVE

  public :: ray_bunch_intensity
  real, pointer, dimension(:,:) :: Coord_DI
  !\
  !Frequency related:  in GGSE
  real :: DensityCr = -1.0, DensityCrInv = -1.0
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
    use BATL_geometry, ONLY: IsCartesianGrid, xyz_to_coord, coord_to_xyz
    use ModMain,ONLY:NameThisComp
    use ModDensityAndGradient,ONLY: NameVector, &
         Density_I, GradDensity_DI, DeltaSNew_I

    implicit none
    !INPUTS
    ! Observer's position
    real, intent(in) :: XyzObserver_D(MaxDim)          
 
    real, intent(in) :: RadioFrequency
    ! (XLower, YLower, XUpper, YUpper)
    real, intent(in) :: ImageRange_I(4)        
    ! Radius of "integration sphere"
    real, intent(in) :: rIntegration           
    ! Dimensions of the raster in pixels
    integer, intent(in) ::  nXPixel, nYPixel
    !\
    ! OUTPUTS  
    ! The result of emissivity integration
    real, dimension(nXPixel,nYPixel), intent(out) :: Intensity_II          
    !/
                            
    !\
    !Misc functions of inputs 
    !/     
    real :: XyzObservLen !=sqrt(sum(XyzObs_D**2))
    !\
    ! Same as ImageRange_I(4)
    real :: XLower, XUpper, YLower, YUpper
    !\
    ! rIntegration related vars:
    real :: rIntegrationSqr !=rIntegration**2 + Eps
    !\
    !nPixel related:

    integer :: nRay !=nXPixelX*nYPixel

    !\
    !Pixel grid
    !/ 
    real :: DxPixel, DyPixel 
    ! Pixel coordinates INSIDE the image plane
    real, dimension(nXPixel,nYPixel) :: XPixel_II, YPixel_II  
    ! Unity vector normal to the image plane 
    real :: Normal_D(MaxDim)                                          
    real :: ZAxisOrt_D(MaxDim) = (/0, 0, 1/)
    ! Image plane inner Cartesian orts
    real :: Tau_D(MaxDim), Xi_D(MaxDim) 
    !\
    ! Pixel 3D Cartesian coordinates
    real, dimension(MaxDim) :: XyzPixel_D

    !\
    ! Initial directions of rays
    !/
    ! Unity slope vectors for all the line-of-sights pointing at 
    ! the pixels          
    real :: SlopeUnscaled_D(MaxDim), Slope_DI(MaxDim,nXPixel*nYPixel)
    !\
    ! Distance from the radiotelescope to the integration sphere  
    real :: ObservToIntSphereDist
    !\
    ! Intersection point of the line of sight with the integration sphere
    !/
    real :: XyzAtIntSphere_D(MaxDim)

    !\
    ! Global vector: generalized coordinates of the ray endpoints
    !/
    

    real :: SolarDistSqr, Xyz_D(MaxDim)   
  
    !\
    ! A ray is excluded from processing if it is .true.
    logical, dimension(nYPixel*nXPixel) :: UnusedRay_I
    real, dimension(nXPixel*nYPixel) :: Intensity_I, DeltaS_I
    !\
    ! Parameters
    real,parameter :: Tolerance = 0.01, DeltaS = 1.0
    real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe
    logical, save :: DoAllocate = .true.
    !\
    !Loop variables
    integer ::  i, j, iRay 
    !------------------------------------

    !
    ! Total number of rays and pixels in the raster
    !
    nRay = nXPixel*nYPixel  
    !
    ! Initial allocation of the global vector of ray positions
    ! and other dynamic arrays
    !
    if(DoAllocate)then
       DoAllocate=.false.
       ! Symbolic name of the global vector:
       NameVector = NameThisComp//'_Pos_DI'
       call allocate_vector(NameVector, 3, nRay)
       call associate_with_global_vector(Coord_DI, NameVector)
       allocate(Density_I(nRay), GradDensity_DI(MaxDim,nRay),DeltaSNew_I(nRay))
    end if
    !
    ! Calculate the critical density from the frequency, in CGS
    !
    DensityCr = cPi*cProtonMass*cElectronMass*1e6 &
         *(RadioFrequency/ProtonChargeSGSe)**2
    DensityCrInv = 1.0/DensityCr

    !
    ! Determine the image plane inner coordinates of pixel centers
    !
    XLower = ImageRange_I(1)
    YLower = ImageRange_I(2)
    XUpper = ImageRange_I(3)
    YUpper = ImageRange_I(4)
    DxPixel = (XUpper - XLower)/nXPixel
    DyPixel = (YUpper - YLower)/nYPixel
    XPixel_II(:,1) = (/ (XLower + (real(i) - 0.5)*DxPixel, i = 1, nXPixel) /)
    do j = 2, nYPixel
       XPixel_II(:,j) = XPixel_II(:,1)
    end do
    YPixel_II(1,:) = (/ (YLower + (real(j) - 0.5)*DyPixel, j = 1, nYPixel) /)
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
    Tau_D = Tau_D/sqrt(sum(Tau_D**2))
    Xi_D = cross_product(Normal_D, Tau_D)

    !
    ! Calculate coordinates of all the pixels in the SGI
    !
    iRay = 0
    do j = 1, nYPixel
       do i = 1, nXPixel
          XyzPixel_D = XPixel_II(i,j)*Tau_D + YPixel_II(i,j)*Xi_D
          SlopeUnscaled_D = XyzPixel_D - XyzObserver_D
          iRay = iRay + 1
          Slope_DI(:,iRay) = SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2)) 

          !\
          ! Find the points on the integration sphere where it intersects
          ! with the straight "rays" 
          !/
    
          !\
          ! Solve a duadratic equation,
          !|| XyzObs_D + ObservToIntSphereDist*Slope_DI || = rIntegration
          !or  ObservToIntSphereDist**2 + 
          !  + 2*ObservToIntSphereDist*XyzObs_D
          !  + XyzObservLen**2 - rIntegration**2 = 0
          !/  
          ObservToIntSphereDist = -sum(XyzObserver_D*Slope_DI(:,iRay))&
               - sqrt(rIntegration**2 - XyzObservLen**2 &
               + sum(Slope_DI(:,iRay)*XyzObserver_D)**2)
          XyzAtIntSphere_D = XyzObserver_D &
               + Slope_DI(:,iRay)*ObservToIntSphereDist
          if(IsCartesianGrid)then
             Coord_DI(:,iRay) = XyzAtIntSphere_D 
          else
             call xyz_to_coord(XyzAtIntSphere_D, Coord_DI(:,iRay))
          end if
       end do
    end do
    !
    ! Do emissivity integration inside of the integration sphere 
    !
    Intensity_I = 0
    UnusedRay_I = .false.
    DeltaS_I = DeltaS
    rIntegrationSqr = rIntegration**2 + 0.01
    do while(.not.all(UnusedRay_I))
       !\
       ! Advance rays through one step.
       call ray_path(nRay, UnusedRay_I, Slope_DI, &
            DeltaS_I, Tolerance, Intensity_I)
       !Exclude rays which are out the integration sphere
       RAYS:do iRay = 1, nRay
          if(UnusedRay_I(iRay))CYCLE RAYS
          if(IsCartesianGrid)then
             SolarDistSqr = sum(Coord_DI(:,iRay)**2)
          else
             call coord_to_xyz(Coord_DI(:,iRay),Xyz_D)
             SolarDistSqr = sum(Xyz_D**2)
          end if
          UnusedRay_I(iRay) = SolarDistSqr .gt. rIntegrationSqr
       end do RAYS
    end do

    Intensity_II = reshape(Intensity_I, (/nXPixel,nYPixel/))

  end subroutine ray_bunch_intensity
  !=========
  subroutine ray_path(nRay, UnusedRay_I, Slope_DI, &
       DeltaS_I, ToleranceInit, Intensity_I)
    !use ModCoordTransform
    !use CON_global_vector, ONLY: associate_with_global_vector
    use ModDensityAndGradient, ONLY: NameVector, get_plasma_density, &
         GradDensity_DI, Density_I, DeltaSNew_I
    use ModPhysics, ONLY   : No2Si_V, UnitRho_, UnitX_
    use ModProcMH, ONLY    : iProc
    use BATL_geometry, ONLY: CoordMin_D, CoordMax_D, &
         IsCartesianGrid, xyz_to_coord, coord_to_xyz
    use ModMain,     ONLY  : TypeCellBc_I
    !   The subroutine ray_path() makes raytracing and emissivity integration
    ! along ray paths. It works on a group of rays (a beam). Each ray is
    ! specified by its Cartesian Xyz_DI and its direction cosines 
    ! Slope_DI. The subroutine calculates new Xyz_DI, which is 
    ! DeltaS_I away from the initial position. It calculates the intensity
    ! along the step as emissivity multiplied by DeltaS and adds this to the
    ! Intensity_I. Thus, after series of calls to ray_path(), the Intensity_I
    ! contains the result of integration along the paths of every ray. 
    !   The electromagnetic rays are refracted in plasma due to its non-zero
    ! refractive index, which is the square root of the dielectric
    ! permittivity \epsilon. The \epsilon, in turn, is a function of the 
    ! plasma density. The external subprogram, get_plasma_density(),
    ! provides the plasma Density_I along with its gradient, GradDensity_DI,
    ! at the Xyz_DI. The \epsilon can be calculated as
    !         1 - Density_I/DensityCr, 
    ! where DensityCr is the "critical" plasma density at which the dielectric
    ! permittivity \epsilon falls down to zero. The value of DensityCr is
    ! proportional to the square of the wave frequency. For example,
    ! for the 42.3 MHz radiowaves the critical density is
    ! ~3.71x10^(-17) g/cm^3. A surface where the plasma density achieves the
    ! critical value acts like a a mirror. No radiation penetrates the
    ! critical surface. The radiowaves can only travel in the regions
    ! with lower density.
    !   The emissivity w of plasma at a selected plasma frequency is
    ! calculated as a polynomial function
    !         w = (Dens2DensCr)^2*(0.5 - Dens2DensCr)^2,
    ! where Dens2DensCr is the quotient Density_I/DensityCr. So, here the
    ! plasma density is used for calculation of both emissivity
    ! and dielectric permittyvity.
    !
    !   The ray refraction and plasma emissivity calculation is based solely
    ! on the plasma density. The density is obtained from the MHD simulation.
    ! The density data are located in the memories of several processors, and
    ! the get_plasma_density() program from the ModDensityAndGradient module
    ! makes the density acquisition. 
    ! 
    !   The parameters of the ray_path() are briefly described below.
    !
    ! get_plasma_density():  external subroutine that returns the plasma
    !     Density_I and its gradient GradDensity_DI. It also provides the
    !     recommended step size DeltaSNew_I and asserts the UnusedRay_i (sets
    !     it to .true.) for a ray, should it (illegally) penetrate into the
    !     region with "negative" dielectric permittivity.
    ! nRay:           number of rays being processed.
    ! UnusedRay_I:  the caller program can use this logical array to stop
    !     processing of any individual ray when it already finished travelling
    !     inside of a specified part of space. Set the corresponding element
    !     of UnusedRay_I to .true. to leave it unprocessed during the
    !     subsequent calls to ray_path(). Before the first call to ray_path()
    !     all the elements of UnusedRay_I should be set to .false.; then all
    !     the rays will be processed.
    ! Xyz_DI:    Cartesian position vectors of the rays.
    ! Slope_DI:       Direction cosines of the rays.
    ! DeltaS_I:       Current step values for the rays. Set the elements of
    !     Delta_I to some reasonable value, say, 1.0. The DeltaS_I elements
    !     are individually modified by ray_path() to satisfy the precision
    !     requirements set by ToleranceInit.
    ! ToleranceInit:   determines the precision of ray paths calculation.
    !     ToleranceInit is the inverse of the minimum number of ray
    !     trajectory points per one radian of the ray curvature. If this
    !     requirement is not satisfied, the corresponding element of DeltaS_I
    !     is decreased. Do not set the ToleranceInit to any value greater than
    !     0.1 (which means 10 points per curvature radian): it will be
    !     internally reduced to 0.1 anyway.   
    ! DensityCr: the plasma density at which its dielectric permittivity
    !     becomes zero for chosen wave frequency.
    ! Intensity_I:    the intensities of each ray calculated as the integral
    !     of emissivity along the ray path. During each call to ray_path(),
    !     each element of the Intensity_I is incremented by the integral along
    !     the step DeltaS_I. Set all the elements of Intensity_I to 0.0 before
    !     the first call to the ray_path().
    ! 
    ! Written by Leonid Benkevitch.
    !
    !
    integer, intent(in) :: nRay   ! # of pixels in the raster
    real,    intent(inout), dimension(MaxDim,nRay) :: Slope_DI
    real,    intent(inout), dimension(nRay) :: Intensity_I, DeltaS_I
    real,    intent(in) ::  ToleranceInit

    ! a ray is excluded from processing if it is .true. 
    logical, intent(inout), dimension(nRay) :: UnusedRay_I   

    logical, save :: IsNewEntry = .true.
   
    integer, parameter :: nSplitDeltaS = 2

    real,    save, dimension(MaxDim)       :: Slope1_D, Omega_D
    real,    save, dimension(MaxDim)       :: ProjSlopeOnMinusGradEps_D
    real,    save, dimension(MaxDim)       :: StepX_D, StepY_D, RelGradRefrInx_D 
    real,    save, dimension(MaxDim)       :: GradDielPerm_D, PositionHalfBack_D

    real, save, allocatable             :: Xyz_DI(:,:)
    real, save, pointer, dimension(:)   :: DistanceToCritSurf_I

    ! IsOKRay_I: .true. for shallow rays; is set to .false. for steep rays.
    logical, save, pointer, dimension(:)   :: IsOKRay_I       
    

    real :: HalfDeltaS                        ! DeltaS halved
    real :: DielPerm, DielPermHalfBack, Dens2DensCr
    real :: Coef, Curv, Curv1
    real :: LCosAl ! L*cos(Alpha),   
    ! where L is inverse grad of \epsilon, Alpha is the incidence angle
    real :: GradDielPermSqr, GradEpsDotSlope
    real :: ParabLen
    real, save :: &
         Tolerance=0.1, ToleranceSqr=1.0e-2, StepMin = 1.0
    
    integer:: iRay, iDim

    character(LEN=20),save:: TypeBoundaryDown_D(MaxDim), TypeBoundaryUp_D(MaxDim)
    !---------------

    if (IsNewEntry) then
       IsNewEntry = .false.
       allocate(Xyz_DI(MaxDim,nRay))
       Xyz_DI(MaxDim,nRay) = 0.0
       allocate(DistanceToCritSurf_I(nRay))
       DistanceToCritSurf_I = 0.0
       allocate(IsOKRay_I(nRay))
       IsOKRay_I = .true.

       !  minimum ten points between a vacuum and a critical surface and
       !  minimum 10 points over 1 rad of the curvature
       Tolerance = min(ToleranceInit,0.1)  

       
       ToleranceSqr = Tolerance**2 
       
       ! One (ten-thousandth) hundredth of average step 
       StepMin = 1e-2*sum(DeltaS_I)/nRay
       do iDim=1,MaxDim
          TypeBoundaryDown_D(iDim) = trim(TypeCellBc_I(2*iDim-1))
          TypeBoundaryUp_D(iDim)   = trim(TypeCellBc_I(2*iDim))
       end do
       if(iProc==0)then
          write(*,*)'TypeBoundaryDown_D=',TypeBoundaryDown_D
          write(*,*)'TypeBoundaryUp_D=',TypeBoundaryUp_D
          write(*,*)'StepMin=',StepMin
       end if
          
    end if

    !Start the predictor step of the GIRARD scheme: 

    do iRay = 1, nRay
       if (UnusedRay_I(iRay)) CYCLE  ! Do not process the rays that are done
       if(IsCartesianGrid)then
          Xyz_DI(:, iRay)= Coord_DI(:,iRay)
       else
          call coord_to_xyz(Coord_DI(:,iRay), Xyz_DI(:, iRay))
       end if
       Xyz_DI(:,iRay) = Xyz_DI(:,iRay) + &
            Slope_DI(:,iRay)*0.50*DeltaS_I(iRay)
       if(IsCartesianGrid)then
          Coord_DI(:,iRay) = Xyz_DI(:, iRay)
       else
          call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
       end if
       ! Now Xyz_DI moved by 1/2 DeltaS !!!
       call check_bc(iRay)
       if(UnusedRay_I(iRay))then
          !\
          !Ray end is behind the boundary of the computational domain
          !Push the end point backward
          !/
          Xyz_DI(:,iRay) = Xyz_DI(:,iRay) - &
               Slope_DI(:,iRay)*0.50*DeltaS_I(iRay)
          if(IsCartesianGrid)then
             Coord_DI(:,iRay) = Xyz_DI(:, iRay)
          else
             call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
          end if
       end if
    end do

    !
    ! The following routine collects the density values at the Xyz_DI
    ! from all relevant processors
    !

    call get_plasma_density(nRay)
    
    
    !In making radio images this is a bad pixel which should be further processed 
    UnusedRay_I = UnusedRay_I .or.  Density_I >= DensityCr ! "bad rays" are done
    
    do iRay = 1, nRay

       if (UnusedRay_I(iRay)) CYCLE  ! Do not process the rays that are done
       if(IsCartesianGrid)then
          Xyz_DI(:, iRay)= Coord_DI(:,iRay)
       else
          call coord_to_xyz(Coord_DI(:,iRay), Xyz_DI(:, iRay))
       end if

       HalfDeltaS = 0.50*DeltaS_I(iRay)
       ! Original Position (at an integer point):
       PositionHalfBack_D = Xyz_DI(:,iRay) - Slope_DI(:,iRay)*HalfDeltaS
       Dens2DensCr = Density_I(iRay)*DensityCrInv

       DielPerm = 1.0 - Dens2DensCr
       GradDielPerm_D = -GradDensity_DI(:,iRay)*DensityCrInv
       GradDielPermSqr = sum(GradDielPerm_D**2)

       GradEpsDotSlope = sum(GradDielPerm_D*Slope_DI(:,iRay))

       DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS

       !Check positivity:
       if( DielPermHalfBack<=0.0)then
          UnusedRay_I(iRay) = .true.
          CYCLE
       end if

       Curv = (0.50*HalfDeltaS/DielPermHalfBack)**2* &
            (GradDielPermSqr - GradEpsDotSlope**2)
       !The curvature squared characterizes the magnitude of
       !the tranverse gradient of the dielectric permittivity

       if (IsOKRay_I(iRay)) then

          !
          ! Check if the trajectory curvature is too sharp
          ! to meet the Tolerance. If so, reduce the DeltaS 
          ! step for the iRay-th ray and leave it until
          ! the next call.
          !
          !

          if (Curv .ge. ToleranceSqr) then
             DeltaS_I(iRay) = 0.99 * &
                  DeltaS_I(iRay)/(2*sqrt(Curv/ToleranceSqr))
             if(IsCartesianGrid)then
                Coord_DI(:,iRay) = PositionHalfBack_D
             else
                call xyz_to_coord(PositionHalfBack_D, Coord_DI(:,iRay))
             end if
             CYCLE
          end if

          ! Test if some of the next points can get into the prohibited
          ! part of space with "negative" dielectric permittivity
          !
          ! Here we check the magnitude of the longitudinal gradient
          ! of the dielectric permittivity

          if (3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack) then
             !
             ! Mark the ray as steep; 
             ! memorize the distance to the critical surface;
             ! reduce step
             !
             IsOKRay_I(iRay) = .false.
             DistanceToCritSurf_I(iRay) = DielPermHalfBack  &
                  /sqrt(GradDielPermSqr)
             DeltaS_I(iRay) =  max(&
                  0.50*Tolerance*DistanceToCritSurf_I(iRay),&
                  StepMin)
             if(IsCartesianGrid)then
                Coord_DI(:,iRay) = PositionHalfBack_D
             else
                call xyz_to_coord(PositionHalfBack_D, Coord_DI(:,iRay))
             end if
             CYCLE
          end if

       end if ! IsOKRay_I

       !
       ! Either switch to opposite pranch of parabola
       ! or make a Boris step
       !

       if ((3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack)&
            .or.(DeltaS_I(iRay)<StepMin.and.&
            GradEpsDotSlope<0)) then

          ! Switch to the opposite branch of parabolic trajectory
          !
          ! When a next step can drive the ray into the area with the
          ! plasma density greater then its critical value, then a special 
          ! technique of "parabolic ray reflection" is employed. 
          ! It can be shown that a ray trajectory in the medium with constant
          ! gradient of dielectric permittivity is a parabola. If the step
          ! DeltaS_I is small enough we can assume the grad \epsilon constant
          ! and hence assume that the ray approaches the critical surface 
          ! in a parabolic path.
          ! We calculate the parameters of the parabola, namely:
          ! StepX_D -- vector along the -grad \epsilon, with the length equal 
          !     to the distance from the PositionHalfBack_D to the parabola 
          !     extremum;  
          ! StepY_D -- vector perpendicular to StepX_D, lying inside of the
          !     parabola plane, pointing at the opposite parabola branch, and
          !     with the length equal the distance from PositionHalfBack_D to
          !     the opposite branch of parabola.
          ! The parabolic reflection is just a replacement of the Xyz_DI
          ! with the symmetric point at the opposite branch of parabola, and 
          ! changing the "incident" direction Slope_DI to the "departing" ray
          ! direction according to Snell law.
          !

          LCosAl = -GradEpsDotSlope/GradDielPermSqr
          ProjSlopeOnMinusGradEps_D = -LCosAl*GradDielPerm_D  
          ! Here v_proj

          StepY_D = Slope_DI(:,iRay) - ProjSlopeOnMinusGradEps_D 
          ! Here c; |c| = sin \alpha

          Slope_DI(:,iRay) = Slope_DI(:,iRay) - 2*ProjSlopeOnMinusGradEps_D

          !
          ! We need to have 
          !   |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*DielPermHalfBack*L
          ! Now the direction of Step_Y is right, the length of it is equal
          ! to  sin(\alpha).
          ! Multiply it by L*Cos(\alpha)*DielPermHalfBack
          !

          StepY_D = 4*StepY_D*LCosAl*DielPermHalfBack 

          Xyz_DI(:,iRay) = PositionHalfBack_D + StepY_D
          if(IsCartesianGrid)then
             Coord_DI(:,iRay) = Xyz_DI(:, iRay)
          else
             call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
          end if
          call check_bc(iRay)
          if(UnusedRay_I(iRay))then
             !\
             !Ray end is behind the boundary of the computational domain
             !Push the end point backward
             !/
             Xyz_DI(:,iRay) = PositionHalfBack_D
             if(IsCartesianGrid)then
                Coord_DI(:,iRay) = Xyz_DI(:, iRay)
             else
                call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
             end if
          end if
          !
          ! Step_X is in the direction of - grad \epsilon,
          ! whose vector is of the length of 1/L.
          ! The length of Step_X is cos^2(\alpha)*L*DielPermHalfBack 
          ! Thus,
          !

          StepX_D = (DielPermHalfBack*LCosAl**2)*GradDielPerm_D


          ParabLen = sqrt(sum((2*StepX_D)**2) + sum(StepY_D**2))
          Intensity_I(iRay) = Intensity_I(iRay)  &
               + ParabLen*(Dens2DensCr**2)*(0.50 - Dens2DensCr)**2
       else 

          ! Make a step using Boris' algorithm

          Coef = 0.50*HalfDeltaS/(1.0 - Dens2DensCr)
          RelGradRefrInx_D = Coef*GradDielPerm_D    
          ! grad(n) = grad(eps(i+1/2))/(2*eps(i+1/2))

          Omega_D = cross_product(RelGradRefrInx_D, Slope_DI(:,iRay))
          Slope1_D = Slope_DI(:,iRay) &
               + cross_product(Slope_DI(:,iRay), Omega_D)

          Omega_D = cross_product(RelGradRefrInx_D, Slope1_D)
          Slope1_D = Slope_DI(:,iRay) &
               + cross_product(Slope_DI(:,iRay), Omega_D)

          Curv1 = sum(Omega_D**2)
          Coef = 2/(1 + Curv1)
          Slope_DI(:,iRay) = Slope_DI(:,iRay) &
               + Coef*cross_product(Slope1_D, Omega_D)

          Xyz_DI(:,iRay) = Xyz_DI(:,iRay) &
               + Slope_DI(:,iRay)*HalfDeltaS
          if(IsCartesianGrid)then
             Coord_DI(:,iRay) = Xyz_DI(:, iRay)
          else
             call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
          end if
          call check_bc(iRay)
          if(UnusedRay_I(iRay))then
             !\
             !Ray end is behind the boundary of the computational domain
             !Push the end point backward
             !/
             Xyz_DI(:,iRay) = Xyz_DI(:,iRay)&
               - Slope_DI(:,iRay)*HalfDeltaS
             if(IsCartesianGrid)then
                Coord_DI(:,iRay) = Xyz_DI(:, iRay)
             else
                call xyz_to_coord(Xyz_DI(:, iRay), Coord_DI(:,iRay))
             end if
          end if
          Intensity_I(iRay) = Intensity_I(iRay) &
               + DeltaS_I(iRay)*(Dens2DensCr**2)*(0.50 - Dens2DensCr)**2
       end if
       if(UnusedRay_I(iRay))CYCLE
       !
       !   The code below makes gradual increases of the DeltaS up to the value
       ! specified in DeltaSNew. The smooth step increase is required so as
       ! not to get into the space behind the critical surface, stepping with
       ! DeltaS that instantly changes from a very little size to the normal
       ! DeltaSNew length. DeltaS is usually reduced in a close vicinity of 
       ! the critical surface, where the ray is travelling along a very sharp
       ! curve with high curvature. For many rays it means fractioning of the
       ! DeltaS down several orders of magnitude, therefore the new step trial
       ! should start from a bigger step of the same order of magnitude.
       !   This problem is solved using a non-linear difference equation:
       !           Y(i+1) = [2 - Y(i)/X(i)]*Y(i),
       ! where X(i) is the desired final DeltaS length from DeltaSNew, and
       ! Y(i) is the actual DeltaS length. A brief analysis of the equation
       ! shows that, when Y is small compared to X, the next value of Y will
       ! be almost 2*X, so the DeltaS would grow in a geometrical progression.
       ! However, as Y approaches X, its growth rate becomes slower. However,
       ! Y always reaches X in several steps. One can check that for Y = X the
       ! next value of Y is always that of X.
       !     The behaviour of Y can also be characterized as exponential
       ! growth when Y is close to zero, and exponential saturation when Y is
       ! close to X. 
       ! 
       if (IsOKRay_I(iRay)) then
          !
          ! For shallow rays the DeltaS is increased unconditionally
          !
          DeltaS_I(iRay) = max(2 - DeltaS_I(iRay)/DeltaSNew_I(iRay),1.0)*&
               min(DeltaSNew_I(iRay), DeltaS_I(iRay))
       else 
          !
          ! If the iRay-th ray is marked as steep (i.e. "not gentle" or "not
          ! shallow") then the code below increases its DeltaS only if the
          ! current distance to the critical surface, calculated as
          !      \epsilon / grad \epsilon,
          ! is greater than this distance value saved along with marking the
          ! ray as steep in the DistanceToCritSurf_I. 
          !   This can happen in two cases: 
          ! (1) either the ray was "parabola reflected" and
          ! after several steps it went away from the surface by the distance
          ! where the parabola switching occurred; 
          ! (2)or the ray is not steep any more because 
          ! the current DeltaS is so small, that the next step does not
          ! penetrate the critical surface.
          !   The ray is reverted to "gentle" or "shallow"
          !
          if (DielPermHalfBack .gt. &
               DistanceToCritSurf_I(iRay)*sqrt(GradDielPermSqr)) then
             IsOKRay_I(iRay) = .true.
             DeltaS_I(iRay) = max(2 - DeltaS_I(iRay)/DeltaSNew_I(iRay), 1.0)*&
                  min(DeltaSNew_I(iRay), DeltaS_I(iRay))
          end if
       end if
    end do
  contains
    !==================
    subroutine check_bc(iRay)
      integer, intent(in)::iRay
      !------------
      UnusedRay_I(iRay) = any(Coord_DI(:,iRay) < CoordMin_D).or.&
           any(Coord_DI(:,iRay) > CoordMax_D)
    end subroutine check_bc
    !====================
  end subroutine ray_path
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
  real, dimension(nXPixel,nYPixel), intent(out) :: Intensity_II          
  !--------------------------------
  call ray_bunch_intensity(XyzObserver_D, RadioFrequency, ImageRange_I, &
       rIntegration, nXPixel, nYPixel, Intensity_II)
end subroutine get_ray_bunch_intensity
