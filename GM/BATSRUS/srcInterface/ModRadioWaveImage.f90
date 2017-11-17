!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!BOP
!MODULE: ModDensityAndGradient - provide the density and gradients at SC_Xyz_DI
!INTERFACE:
module ModDensityAndGradient
  !USES:
  use MH_domain_decomposition
  use CON_global_message_pass
  use ModMain, ONLY: MaxDim
  use ModProcMH, ONLY:iComm
  use ModMpi
  implicit none
  private !Except
  logical,save :: DoInit=.true.
  character(LEN=10), save :: NameVector, NameMask
  type(RouterType),save::Router
  type(GridDescriptorType),save::LineGrid,MhGrid
  type(DomainDecompositionType),save::LineDD
  real,allocatable,save,dimension(:)::Density_I,DeltaSNew_I
  real,allocatable,save,dimension(:,:)::GradDensity_DI

  !PUBLIC MEMBERS:
  public::get_plasma_density, NameVector, NameMask
  public::GradDensity_DI,Density_I,DeltaSNew_I
  !EOP

contains
  !==========================
  subroutine get_plasma_density(nRay)
    !
    ! Calculates plasma density, Density_I, and its gradient, 
    ! GradDensity_DI(3,nRay), at specified locations Position_DI(3,nRay)
    ! Also, it provides appropriate step, DeltaSNew_I, conceivably dependent
    ! on the numeric grid size
    !
    integer,intent(in)::nRay
    !\
    !Misc
    integer :: nU_I(2), iError

    call timing_start('get_plasma_density')
    if(DoInit)then
       DoInit=.false.
       call set_standard_grid_descriptor(&
            MH_DomainDecomposition,GridDescriptor=MhGrid)

       nU_I=ubound_vector(NameVector)
       call init_decomposition(LineDD,&
            compid_grid(MhGrid%DD%Ptr),1,&
            IsLocal=.true.)
       call get_root_decomposition(&
            LineDD,&                 !GridDescroptor to be constructed
            iRootMapDim_D=(/1/),&!The block amount, along each direction(D)
            XyzMin_D=(/cHalf/),&      !Minimal gen. coordinates, along each D 
            XyzMax_D=(/cHalf+nU_I(2)/),& !Maximal gen. coordinates, along each D
            nCells_D=(/nU_I(2)/))
       call set_standard_grid_descriptor(LineDD,&
            GridDescriptor=LineGrid)
       call init_router(&
            MhGrid,& !GridDesctriptor for the source field (in) 
            LineGrid,&   !GirdDescriptor,save,intent(out)
            Router,&   !resulting router, intent(out)
            nIndexTarget=1)
       DoInit=.false.
    end if

    call set_router(& 
         GridDescriptorSource=MhGrid,&
         GridDescriptorTarget=LineGrid,&
         Router=Router,&
         NameMappingVector=NameVector,&
         NameMask=NameMask, &
         interpolate=interpolation_amr_gc)

    call global_message_pass(Router=Router,&
         nVar=MaxDim+1+1,&
         fill_buffer=get_density_local,&
         apply_buffer=put_density_value)

    call MPI_BCAST(GradDensity_DI(1,1),3*nRay,MPI_REAL,0,iComm,iError)
    call MPI_BCAST(Density_I(1),         nRay,MPI_REAL,0,iComm,iError)
    call MPI_BCAST(DeltaSNew_I(1),       nRay,MPI_REAL,0,iComm,iError)
    call timing_stop('get_plasma_density')
  end subroutine get_plasma_density
  !================================================
  subroutine get_density_local(&
       nPartial,iGetStart,Get,W,State_V,nVar)
    !USES:
    use ModAdvance, ONLY: State_VGB
    use ModVarIndexes, ONLY: Rho_
    use ModCellGradient, ONLY: GradDensity_DGB => GradVar_DGB
    use ModGeometry,ONLY: CellSize_DB, x_, y_, z_
    use ModPhysics, ONLY: No2Si_V, UnitRho_
    use CON_router
    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::State_V

    integer::iGet, i, j, k, iBlock
    real :: Weight
    !----------------------------------------------------------
    i      = Get%iCB_II(1,iGetStart)
    j      = Get%iCB_II(2,iGetStart)
    k      = Get%iCB_II(3,iGetStart)
    iBlock = Get%iCB_II(4,iGetStart)
    Weight = W%Weight_I(iGetStart)
    State_V(1:MaxDim)= Weight*GradDensity_DGB(1:MaxDim,i,j,k,iBlock)
    State_V(MaxDim+1)= Weight*&
         State_VGB(rho_,i,j,k,iBlock)
    State_V(MaxDim+1+1)=Weight*&
         minval(CellSize_DB(:,iBlock))

    do iGet=iGetStart+1,iGetStart+nPartial-1
       i      = Get%iCB_II(1,iGet)
       j      = Get%iCB_II(2,iGet)
       k      = Get%iCB_II(3,iGet)
       iBlock = Get%iCB_II(4,iGet)
       Weight = W%Weight_I(iGet)
       State_V(1:MaxDim)= State_V(1:MaxDim) + Weight*&
            GradDensity_DGB(1:MaxDim,i,j,k,iBlock)
       State_V(MaxDim+1)  = State_V(MaxDim+1)+Weight*&
            State_VGB(rho_,i,j,k,iBlock)
       State_V(MaxDim+1+1)= State_V(MaxDim+1+1)+Weight*&
            minval(CellSize_DB(:,iBlock))

    end do
    !Convert density to SI
    State_V(1:MaxDim+1) = State_V(1:MaxDim+1) * No2Si_V(UnitRho_)
  end subroutine get_density_local

  !====================================================================

  subroutine put_density_value(nPartial,&
       iPutStart,&
       Put,&
       W,&
       DoAdd,&
       Buff_I,nVar)
    implicit none
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_I
    integer::iCell

    iCell=Put%iCB_II(1,iPutStart)
    !Convert densities from kg/m3 to g/cm3, 
    !the transformation coefficient is 1.0e-3
    if(DoAdd)then
       GradDensity_DI(:,iCell)= GradDensity_DI(:,iCell)+&
            Buff_I(1:MaxDim) * 1.0e-3
       Density_I(iCell)= Density_I(iCell)+&
            Buff_I(MaxDim+1) * 1.0e-3
       DeltaSNew_I(iCell) = DeltaSNew_I(iCell)+&
            Buff_I(MaxDim+1+1)
    else
       GradDensity_DI(:,iCell)= &
            Buff_I(1:MaxDim)  * 1.0e-3
       Density_I(iCell)= &
            Buff_I(MaxDim+1)  * 1.0e-3
       DeltaSNew_I(iCell) = &
            Buff_I(MaxDim+1+1)
    end if
  end subroutine put_density_value
end module ModDensityAndGradient
!===============================================
module ModRadioWaveImage
  use ModConst, ONLY: cPi, cElectronMass, cProtonMass
  use ModCoordTransform, ONLY: cross_product
  use CON_global_vector, ONLY: allocate_vector, associate_with_global_vector,&
       allocate_mask, associate_with_global_mask
  use BATL_lib, ONLY: iProc, MaxDim, nDim, check_interpolate_amr_gc, &
       Particle_I, MaxBlock, nI, nJ, nK, nBlock, &
       Unused_B, &
       message_pass_particles, remove_undefined_particles, &
       mark_undefined, Particle_I, check_particle_location
  use ModParticles, ONLY: allocate_particles
  use ModBatlInterface, ONLY:  interpolate_grid_amr_gc
  implicit none
  SAVE

  public :: ray_bunch_intensity
  real, pointer    :: Coord_DI(:,:), State_VI(:,:)
  integer, pointer :: iIndex_II(:,:)
  logical, pointer :: IsMask_I(:)
  !\
  !Frequency related:  in GGSE
  real :: DensityCr = -1.0, DensityCrInv = -1.0
  !\
  !Tolerance parameters::
  !  minimum 100 points between a vacuum and a critical surface and
  !  minimum 100 points over 1 rad of the curvature
  real, parameter :: Tolerance =  1.0e-2, Tolerance2   =  1.0e-4
 
  !\
  ! Initial step size
  real,parameter  :: DeltaS = 1.0
  ! One  hundredth of average step        
  real, parameter :: StepMin = 1.0e-2
  !\
  ! Radius squared of a "soft boundary". If the ray goes beyond 
  ! this boundary, the integration for this ray stops.
  real :: rIntegration2 !=rIntegration**2 + Eps
  !\
  ! Kind of particles used to integrate emissivity over rays 
  integer :: KindRay_ = -1
  !\
  ! Components of the ray vector
  !\
  ! Cartesian coordinates
  integer, parameter :: x_ = 1, y_ = 2, z_ = 3
  ! Unity slope vectors for all the line-of-sights pointing at 
  ! the pixels
  integer, parameter :: SlopeX_ = 4, SlopeY_ = 5, SlopeZ_ = 6 
  !\
  ! Integration step
  integer, parameter :: Ds_ = 7
  !\
  ! Distance to a critical surface, at which Density = DensityCr
  integer, parameter :: Dist2Cr_ = 8
  integer, parameter :: nVar = Dist2Cr_
  !\
  !Indices
  integer, parameter :: Ray_ = 1, Steepness_ = 2, OK_ = 0, Bad_ = 1
  integer, parameter :: Status_ = 3, Predictor_ = 1
  integer, parameter :: nIndex = Steepness_
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
    use ModDensityAndGradient,ONLY: NameVector, NameMask, &
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
    !nPixel related:

    integer :: nRay !=nXPixelX*nYPixel

    !\
    !Pixel grid
    !/ 
    real :: DxPixel, DyPixel 
    ! Pixel coordinates INSIDE the image plane
    real :: XPixel, YPixel  
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
             
    real :: SlopeUnscaled_D(MaxDim)
    !\
    ! Distance from the radiotelescope to the integration sphere  
    real :: ObservToIntSphereDist
    !\
    ! Intersection point of the line of sight with the integration sphere
    !/
    real :: XyzAtIntSphere_D(MaxDim)

    !\
    ! A ray is excluded from processing if it is .true.
    logical, dimension(nYPixel*nXPixel) :: UnusedRay_I
    real, dimension(nXPixel*nYPixel) :: Intensity_I
    !\
    ! Parameters
    real, parameter :: ProtonChargeSGSe = 4.8e-10 !SGSe
    logical, save   :: DoAllocate = .true.
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
       NameMask   = NameThisComp//'_Mask_I' 
       call allocate_vector(NameVector, 3, nRay)
       NameMask   = NameThisComp//'_Mask_I' 
       call allocate_mask(NameMask, nRay)
       call associate_with_global_vector(Coord_DI, NameVector)
       call associate_with_global_mask(IsMask_I, NameMask)
       call allocate_particles(&
            iKindParticle = KindRay_, &
            nVar          = nVar    , &
            nIndex        = nIndex  , &
            nParticleMax  = nRay      )
       nullify(State_VI);   State_VI => Particle_I(KindRay_)%State_VI
       nullify(iIndex_II); iIndex_II => Particle_I(KindRay_)%iIndex_II
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
       YPixel = YLower + (real(j) - 0.5)*DyPixel
       do i = 1, nXPixel
          XPixel = XLower + (real(i) - 0.5)*DxPixel
          XyzPixel_D = XPixel*Tau_D + YPixel*Xi_D
          SlopeUnscaled_D = XyzPixel_D - XyzObserver_D
          iRay = iRay + 1
          State_VI(SlopeX_:SlopeZ_,iRay) = &
               SlopeUnscaled_D/sqrt(sum(SlopeUnscaled_D**2)) 

          !\
          ! Find the points on the integration sphere where it intersects
          ! with the straight "rays" 
          !/
    
          !\
          ! Solve a quadratic equation,
          !   XyzObs_D + ObservToIntSphereDist*Slope_DI || = rIntegration
          !or  ObservToIntSphereDist**2 + 
          !  + 2*ObservToIntSphereDist*XyzObs_D
          !  + XyzObservLen**2 - rIntegration**2 = 0
          !/  
          ObservToIntSphereDist = -sum(XyzObserver_D*&
               State_VI(SlopeX_:SlopeZ_,iRay))&
               - sqrt(rIntegration**2 - XyzObservLen**2 &
               + sum(State_VI(SlopeX_:SlopeZ_,iRay)&
               *XyzObserver_D)**2)
          State_VI(x_:z_,iRay) = XyzObserver_D &
               + State_VI(SlopeX_:SlopeZ_,iRay)&
               *ObservToIntSphereDist
          iIndex_II(Ray_,iRay) = iRay
       end do
    end do
    !
    ! Do emissivity integration inside of the integration sphere 
    !
    Intensity_I = 0
    UnusedRay_I = .false.
    State_VI(Ds_,:) = DeltaS
    rIntegration2 = rIntegration**2 + 0.01
    do while(.not.all(UnusedRay_I))
       !\
       ! Advance rays through one step.
       call ray_path(nRay, UnusedRay_I, Intensity_I)
    end do

    Intensity_II = reshape(Intensity_I, (/nXPixel,nYPixel/))

  end subroutine ray_bunch_intensity
  !=========
  subroutine ray_path(nRay, UnusedRay_I, Intensity_I)
    use ModDensityAndGradient, ONLY: NameVector, get_plasma_density, &
         GradDensity_DI, Density_I, DeltaSNew_I
    use ModPhysics, ONLY   : No2Si_V, UnitRho_, UnitX_
    use ModProcMH, ONLY    : iProc
    use BATL_geometry, ONLY: IsCartesianGrid, xyz_to_coord
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
    real,    intent(inout), dimension(nRay) :: Intensity_I
  
    ! a ray is excluded from processing if it is .true. 
    logical, intent(inout), dimension(nRay) :: UnusedRay_I   
    !\
    ! Predictor
    real :: HalfDeltaS                        
    !\
    ! Coordinates at the begining of stage 
    !/
    real  ::  PositionHalfBack_D(MaxDim)
    !\
    ! Density/Density_cr
    !/
    real  ::  Dens2DensCr
    !\
    ! Dielectric permeability, 1 - Dens2DensCr
    !/
    real  :: DielPerm
    !\
    ! Gradient of dielectric permeability
    !/
    real  ::  GradDielPerm_D(MaxDim)
    !\
    ! Gradient of dielectric permeability squared
    real  ::  GradDielPerm2
 
    real :: DielPermHalfBack
    real :: Curv  
    real :: GradEpsDotSlope
    !\
    ! Loop variable
    !/
    integer:: iParticle
    !\
    ! Output from check_step program. If true, the integration step
    ! should be repeated with reduced timestep
    !/
    logical :: DoReturn = .false.
    !\
    ! Heliocentric distance squared, to compare with 
    ! rIntegration squared
    real ::  SolarDist2  
    !----------------------------

    !Start the predictor step of the GIRARD scheme: 

    do iParticle = 1, nRay
       !\
       ! Do not process the rays that are done
       !/
       if (UnusedRay_I(iParticle)) CYCLE
       !\
       ! For the slope calculated at the previous
       ! stage we advance the ray at half step  
       State_VI(x_:z_,iParticle) = &
            State_VI(x_:z_,iParticle) + &
            State_VI(SlopeX_:SlopeZ_,iParticle)&
            *0.50*State_VI(Ds_,iParticle)
       if(IsCartesianGrid)then
          Coord_DI(:,iParticle) = State_VI(x_:z_, iParticle)
       else
          call xyz_to_coord(State_VI(x_:z_, iParticle), &
               Coord_DI(:,iParticle))
       end if
       ! Now Xyz_DI moved by 1/2 DeltaS !!!
       !\
       ! Check if at half step the ray escaped the domain
       call check_particle_location(KindRay_, iParticle,&
            IsGone = UnusedRay_I(iParticle))
    end do

    !
    ! The following routine collects the density values at the Xyz_DI
    ! from all relevant processors
    !
    IsMask_I(1:nRay) = .not.UnusedRay_I(1:nRay)
    call get_plasma_density(nRay)
    
    
    !In making radio images the ray accidentally reachig the overdense plasma 
    !region  should not be further processed 
    UnusedRay_I = UnusedRay_I .or.  Density_I >= DensityCr ! "bad rays" are done
    
    do iParticle = 1, nRay
       !\
       ! Do not process the rays that are done
       !/
       if (UnusedRay_I(iParticle)) CYCLE  
 
       HalfDeltaS = 0.50*State_VI(Ds_,iParticle)
       ! Original Position (at an integer point):
       PositionHalfBack_D = State_VI(x_:z_,iParticle) - &
            State_VI(SlopeX_:SlopeZ_,iParticle)*HalfDeltaS
       Dens2DensCr = Density_I(iParticle)*DensityCrInv

       DielPerm       = 1.0 - Dens2DensCr
       GradDielPerm_D = -GradDensity_DI(:,iParticle)*DensityCrInv
       GradDielPerm2  = sum(GradDielPerm_D**2)

       GradEpsDotSlope = sum(GradDielPerm_D*&
            State_VI(SlopeX_:SlopeZ_,iParticle))

       DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS

       !Check positivity:
       if( DielPermHalfBack<=0.0)then
          UnusedRay_I(iParticle) = .true.
          CYCLE
       end if


       if (iIndex_II(Steepness_,iParticle)==OK_) then
          call check_step(iParticle, DoReturn)
          if(DoReturn) CYCLE
       end if 

       !
       ! Either switch to opposite pranch of parabola
       ! or make a Boris step
       call advance_ray(iParticle)
       call check_particle_location(KindRay_, iParticle,&
            IsGone = UnusedRay_I(iParticle))
       !Exclude rays which are out the integration sphere
       SolarDist2 = sum(State_VI(x_:z_,iParticle)**2)
       UnusedRay_I(iParticle) = UnusedRay_I(iParticle)&
            .or. SolarDist2 > rIntegration2
       if(UnusedRay_I(iParticle))CYCLE
       !\
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
       !/
       !\
       ! For shallow rays the DeltaS is increased unconditionally
       !/
       !\
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
       !/
       if (iIndex_II(Steepness_,iParticle)==OK_.or.DielPermHalfBack .gt. &
               State_VI(Dist2Cr_,iParticle)*sqrt(GradDielPerm2)) then
          iIndex_II(Steepness_,iParticle) = OK_
          State_VI(Ds_,iParticle) = &
               max(2 - State_VI(Ds_,iParticle)/&
               DeltaSNew_I(iParticle), 1.0)*&
               min(DeltaSNew_I(iParticle), &
               State_VI(Ds_,iParticle))
       end if
    end do
  contains
    !==================
    subroutine check_step(iParticle, DoReturn)
      integer, intent(in) :: iParticle
      logical, intent(out):: DoReturn
      !\
      ! Misc
      real :: Curv
      !-------------------
      DoReturn = .false.
      !\
      !The curvature squared characterizes the magnitude of
      !the tranverse gradient of the dielectric permittivity
      Curv = (0.50*HalfDeltaS/DielPermHalfBack)**2* &
           (GradDielPerm2 - GradEpsDotSlope**2)

      if (Curv .ge. Tolerance2) then
         !\
         ! Check if the trajectory curvature is too sharp
         ! to meet the Tolerance. If so, reduce the DeltaS 
         ! step for the iRay-th ray and leave it until
         ! the next call.
         !/
         State_VI(Ds_,iParticle) = 0.99 * &
              State_VI(Ds_,iParticle)/&
              (2*sqrt(Curv/Tolerance2))
         State_VI(x_:z_,iParticle) = &
              PositionHalfBack_D
         DoReturn = .true.
         RETURN
      end if

      ! Check if some of the next points can get into the prohibited
      ! part of space with "negative" dielectric permittivity
      !
      ! Here we check the magnitude of the longitudinal gradient
      ! of the dielectric permittivity

      if (3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack) then
         !
         ! Mark the ray as steep; 
         ! store the distance to the critical surface;
         ! reduce step
         !
         iIndex_II(Steepness_,iParticle) = Bad_
         State_VI(Dist2Cr_,iParticle) = DielPermHalfBack  &
              /sqrt(GradDielPerm2)
         State_VI(Ds_,iParticle) =  max(&
              0.50*Tolerance*State_VI(Dist2Cr_,iParticle),&
              StepMin)
         State_VI(x_:z_,iParticle) = PositionHalfBack_D
         DoReturn = .true.
      end if
    end subroutine check_step
    !==================
    subroutine advance_ray(iParticle)
      integer, intent(in):: iParticle
      !\
      ! Omega = \nabla n/n\times d\vec{x}/ds
      ! Analogous to the magnetic field times q/m in the
      ! equation for the particle gyration
      !/
      real  ::  Omega_D(MaxDim)
      !\
      ! Slope at the stage of predictor
      !/
      real  ::  Slope1_D(MaxDim)
      
      real,    dimension(MaxDim)       :: ProjSlopeOnMinusGradEps_D
      real,    dimension(MaxDim)       :: StepX_D, StepY_D
      ! Below, L is inverse grad of \epsilon, Alpha is the incidence angle
      real :: LCosAl ! L*cos(Alpha)
      ! Length of parabola
      real :: ParabLen
      !\
      ! grad(n)/n, with n=\sqrt(\varepsilon)
      !/
      real  ::  RelGradRefrInx_D(MaxDim)
      !\
      !Misc
      real :: Coef, Curv  
      !---------------
      if ((3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack)&
           .or.(State_VI(Ds_,iParticle)<StepMin.and.&
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

         LCosAl = -GradEpsDotSlope/GradDielPerm2
         ProjSlopeOnMinusGradEps_D = -LCosAl*GradDielPerm_D  
         ! Here v_proj

         StepY_D = State_VI(SlopeX_:SlopeZ_,iParticle) &
              - ProjSlopeOnMinusGradEps_D 
         ! Here c; |c| = sin \alpha

         State_VI(SlopeX_:SlopeZ_,iParticle) = &
              State_VI(SlopeX_:SlopeZ_,iParticle) &
              - 2*ProjSlopeOnMinusGradEps_D

         !
         ! We need to have 
         !   |Step_Y| = 4 * sin(\alpha)*cos(\alpha)*DielPermHalfBack*L
         ! Now the direction of Step_Y is right, the length of it is equal
         ! to  sin(\alpha).
         ! Multiply it by L*Cos(\alpha)*DielPermHalfBack
         !

         StepY_D = 4*StepY_D*LCosAl*DielPermHalfBack 

         State_VI(x_:z_,iParticle) = &
              PositionHalfBack_D + StepY_D
         !
         ! Step_X is in the direction of - grad \epsilon,
         ! whose vector is of the length of 1/L.
         ! The length of Step_X is cos^2(\alpha)*L*DielPermHalfBack 
         ! Thus,
         !

         StepX_D = (DielPermHalfBack*LCosAl**2)*GradDielPerm_D


         ParabLen = sqrt(sum((2*StepX_D)**2) + sum(StepY_D**2))
         Intensity_I(iParticle) = Intensity_I(iParticle)  &
              + ParabLen*(Dens2DensCr**2)*(0.50 - Dens2DensCr)**2
      else 

         ! Make a step using Boris' algorithm
         !\
         ! grad(n)/n = grad(eps(i+1/2))/(2*eps(i+1/2))
         ! RelGradRefrInx_D =  grad(n)/n*ds/2
         !/
         Coef = 0.50*HalfDeltaS/(1.0 - Dens2DensCr)
         RelGradRefrInx_D = Coef*GradDielPerm_D    

         Omega_D = cross_product(RelGradRefrInx_D, &
              State_VI(SlopeX_:SlopeZ_,iParticle))
         Slope1_D = State_VI(SlopeX_:SlopeZ_,iParticle) &
              + cross_product(&
              State_VI(SlopeX_:SlopeZ_,iParticle), Omega_D)

         Omega_D = cross_product(RelGradRefrInx_D, Slope1_D)
         Slope1_D = State_VI(SlopeX_:SlopeZ_,iParticle) &
              + cross_product(&
              State_VI(SlopeX_:SlopeZ_,iParticle), Omega_D)

         Curv = sum(Omega_D**2)
         Coef = 2/(1 + Curv)
         State_VI(SlopeX_:SlopeZ_,iParticle) = &
              State_VI(SlopeX_:SlopeZ_,iParticle) &
              + Coef*cross_product(Slope1_D, Omega_D)

         State_VI(x_:z_,iParticle) = &
              State_VI(x_:z_,iParticle) &
              + State_VI(SlopeX_:SlopeZ_,iParticle)&
              *HalfDeltaS
         Intensity_I(iParticle) = Intensity_I(iParticle) &
              + State_VI(Ds_,iParticle)*&
              (Dens2DensCr**2)*(0.50 - Dens2DensCr)**2
      end if
    end subroutine advance_ray
    !====================
    subroutine get_local_density(&
         iParticle, Density, GradDensity_D, StepSize, IsBody)
         !USES:
      use ModAdvance, ONLY: State_VGB
      use ModVarIndexes, ONLY: Rho_
      use ModCellGradient, ONLY: GradDensity_DGB => GradVar_DGB
      use ModGeometry,ONLY: CellSize_DB, x_, y_, z_
      use ModPhysics, ONLY: No2Si_V, UnitRho_
      ! returns the direction of magnetic field for a given particle
      integer, intent(in):: iParticle
      real,    intent(out):: Density
      real,    intent(out):: GradDensity_D(MaxDim)
      real,    intent(out):: StepSize
      logical, intent(out):: IsBody
      
      !Coordinates and block #
      real     :: Xyz_D(MaxDim)
      integer  :: iBlock
      ! interpolation data: number of cells, cell indices, weights
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim)
      character(len=200):: StringError
      character(len=*), parameter:: NameSub = "get_b_dir"
      !----------------------------------------------------------------------
      !Coordinates and block #
      Xyz_D   = State_VI(x_:z_, iParticle)
      iBlock  = iIndex_II(0,iParticle)
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I,&
           IsBody)
      if(IsBody)then
         call mark_undefined(KindRay_,iParticle);RETURN
      end if
      Density = 0
      do iCell = 1, nCell
         i_D = 1
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         Density = Density + &
              State_VGB(Rho_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
      end do
    end subroutine get_local_density
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
