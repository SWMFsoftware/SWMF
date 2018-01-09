!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModRadioWaveImage
  use ModConst, ONLY: cPi, cElectronMass, cElectronChargeSquaredJm
  use BATL_lib, ONLY: iProc, MaxDim, nDim, x_, y_, z_
  use ModParticles, ONLY: allocate_particles,&
       mark_undefined, Particle_I, check_particle_location,&
       put_particles, trace_particles, interpolate_grid_amr_gc
  implicit none
  SAVE
  PRIVATE!Except
  !Public members
  public :: ray_bunch_intensity, check_allocate
  ! Intensity_I:    the intensities of each ray calculated as the integral
  !     of emissivity along the ray path. During each call to ray_path(),
  !     each element of the Intensity_I is incremented by the integral along
  !     the step DeltaS_I. Set all the elements of Intensity_I to 0.0 before
  !     the first call to the ray_path().
  public:: Intensity_I, StateIn_VI
  real, allocatable :: Intensity_I(:)
  real, allocatable :: StateIn_VI(:,:)

  public:: nRay, SlopeX_, SlopeZ_
  integer :: nRay !=nXPixelX*nYPixel
  !\
  !Frequency related:  
  real :: NElectronSiCr = -1.0, NElectronSiCrInv = -1.0
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
  public:: rIntegration2
  !\
  ! Radius squared of a "soft boundary". If the ray goes beyond 
  ! this boundary, the integration for this ray stops.
  real :: rIntegration2 !=rIntegration**2 + Eps
  !\
  ! Kind of particles used to integrate emissivity over rays 
  integer :: KindRay_ = -1
  !\
  ! Pointers to the arrays
  real,    pointer  :: State_VI(:,:)
  integer, pointer  :: iIndex_II(:,:)
  !\
  ! Components of the ray vector
  !\
  ! Cartesian coordinates: x_ = 1, y_ = 2, z_ = 3
  ! Unity slope vectors for all the line-of-sights pointing at 
  ! the pixels
  integer, parameter :: SlopeX_ = 4, SlopeZ_ = 6 
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
  integer, parameter :: nIndex = Status_
contains
  !=========================================================
  subroutine check_allocate
    logical, save   :: DoAllocate = .true.
    !\
    ! Initial allocation of the global vector of ray positions
    ! and other dynamic arrays
    !/
    if(DoAllocate)then
       DoAllocate = .false.
       allocate(Intensity_I(nRay),StateIn_VI(x_:SlopeZ_,1:nRay))
       nullify(State_VI); nullify(iIndex_II)
    end if
    if(nRay > size(Intensity_I))then
       deallocate(Intensity_I,StateIn_VI)
       allocate(Intensity_I(nRay),StateIn_VI(x_:SlopeZ_,1:nRay))
    end if
    call allocate_particles(&
         iKindParticle = KindRay_, &
         nVar          = nVar    , &
         nIndex        = nIndex  , &
         nParticleMax  = nRay      )
    State_VI  => Particle_I(KindRay_)%State_VI
    iIndex_II => Particle_I(KindRay_)%iIndex_II
  end subroutine check_allocate
  !=========================================================
  subroutine ray_bunch_intensity(RadioFrequency)
    use ModProcMH, ONLY: nProc, iComm, iProc
    use ModMpi
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
    ! RadioFrequency:  frequency in Hertz at wich the radiotelescope works.
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
    real, intent(in) :: RadioFrequency
    !\
    integer :: iError
    !------------------------------------
    !
    ! Calculate the critical density from the frequency
    !
    NElectronSiCr = cPi*cElectronMass&
         *RadioFrequency**2/cElectronChargeSquaredJm
    NElectronSiCrInv = 1.0/NElectronSiCr
 

    call put_particles(KindRay_ ,         &
         StateIn_VI         = StateIn_VI, &
         DoReplace          = .true.)
    !\
    ! Do emissivity integration inside of the integration sphere 
    !/
    Intensity_I = 0
    State_VI(Ds_,:     )    = DeltaS
    State_VI(Dist2Cr_,:)    = 0.0
    iIndex_II(Steepness_:Status_,:) = OK_
    call trace_particles(KindRay_, ray_path)
    if(nProc > 1) call MPI_reduce_real_array(Intensity_I(1), nRay, MPI_SUM, &
         0, iComm, iError)
  end subroutine ray_bunch_intensity
  !=========
  subroutine ray_path(iParticle, IsEndOfSegment)
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
    ! get_local_density():  external subroutine that returns the plasma
    !     Density and its gradient GradDensity_D. It also provides the
    !     recommended step size DeltaSNew.
    ! ToleranceInit:   determines the precision of ray paths calculation.
    !     ToleranceInit is the inverse of the minimum number of ray
    !     trajectory points per one radian of the ray curvature. If this
    !     requirement is not satisfied, the corresponding element of DeltaS_I
    !     is decreased. Do not set the ToleranceInit to any value greater than
    !     0.1 (which means 10 points per curvature radian): it will be
    !     internally reduced to 0.1 anyway.   
    ! DensityCr: the plasma density at which its dielectric permittivity
    !     becomes zero for chosen wave frequency.
    ! Written by Leonid Benkevitch.
    !
    !
    integer, intent(in)   :: iParticle
    logical, intent(out)  :: IsEndOfSegment
    real ::    GradNElectronSi_D(MaxDim), NElectronSi, DeltaSNew  
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
    real  ::  DielPerm
    !\
    ! Gradient of dielectric permeability
    !/
    real  ::  GradDielPerm_D(MaxDim)
    !\
    ! Gradient of dielectric permeability squared
    real  ::  GradDielPerm2
    !\
    ! Misc
    real :: DielPermHalfBack  
    real :: GradEpsDotSlope
    !\
    ! Output from check_step program. If true, the integration step
    ! should be repeated with reduced timestep
    !/
    logical :: DoReturn, IsGone, DoMove
    !----------------------------
    IsEndOfSegment = .false.
    !\
    !Start the predictor step of the GIRARD scheme: 
    if(iIndex_II(Status_,iParticle)/=Predictor_)then
       !\
       ! For the slope calculated at the previous
       ! stage we advance the ray at half step  
       State_VI(x_:z_,iParticle) = State_VI(x_:z_,iParticle) + &
            State_VI(SlopeX_:SlopeZ_,iParticle)&
            *0.50*State_VI(Ds_,iParticle)
       ! Now Xyz_DI moved by 1/2 DeltaS !!!
       !\
       ! Check if at half step the ray escaped the domain
       call check_particle_location(KindRay_, iParticle,&
            DoMove=DoMove, IsGone = IsGone)
       if(IsGone)then
          IsEndOfSegment = .true.
          RETURN
       elseif(DoMove)then
          iIndex_II(Status_,iParticle) = Predictor_
          IsEndOfSegment = .true.
          RETURN
       end if
    else
       iIndex_II(Status_,iParticle) = OK_
    end if
    !\
    ! The following routine calculates density, gradient and 
    ! step size
    !/
    call get_local_density(NElectronSi, GradNElectronSi_D, DeltaSNew,&
         IsEndOfSegment)
    !In making radio images the ray accidentally reachig the overdense plasma 
    !region  should not be further processed 
    if(IsEndOfSegment.or.  NElectronSi >= NElectronSiCr)then
       IsEndOfSegment = .true.
       call mark_undefined(KindRay_,iParticle)
       RETURN
    end if
    !\
    ! Do not process the rays that are done
    !/ 
    HalfDeltaS = 0.50*State_VI(Ds_,iParticle)
    ! Original Position (at an integer point):
    PositionHalfBack_D = State_VI(x_:z_,iParticle) - &
         State_VI(SlopeX_:SlopeZ_,iParticle)*HalfDeltaS
    Dens2DensCr = NElectronSi*NElectronSiCrInv
    
    DielPerm       = 1.0 - Dens2DensCr
    GradDielPerm_D = -GradNElectronSi_D*NElectronSiCrInv
    GradDielPerm2  = sum(GradDielPerm_D**2)

    GradEpsDotSlope = sum(GradDielPerm_D*&
         State_VI(SlopeX_:SlopeZ_,iParticle))
    
    DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS
    
    !Check positivity:
    if( DielPermHalfBack<=0.0)then
       IsEndOfSegment = .true.
       call mark_undefined(KindRay_,iParticle)
       RETURN
    end if

    if (iIndex_II(Steepness_,iParticle)==OK_) then
       call check_step(DoReturn)
       if(DoReturn) RETURN
    end if
    !
    ! Either switch to opposite pranch of parabola
    ! or make a Boris step
    call advance_ray
    call check_particle_location(KindRay_, iParticle,&
         IsGone = IsEndOfSegment)
    !Exclude rays which are out the integration sphere
    if(IsEndOfSegment.or.&
         sum(State_VI(x_:z_,iParticle)**2) > rIntegration2)then
       IsEndOfSegment = .true.
       call mark_undefined(KindRay_,iParticle)
       RETURN
    end if
    call get_new_step_size(DeltaSNew)
  contains
    !==================
    subroutine check_step(DoReturn)
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
    subroutine advance_ray
      use ModCoordTransform, ONLY: cross_product
      use ModUser, ONLY: user_material_properties
      use ModVarIndexes, ONLY: nVar
      use ModWaves,      ONLY: nWave, WaveFirst_, WaveLast_, FrequencySi_W
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
      
      real :: ProjSlopeOnMinusGradEps_D(MaxDim)
      real :: StepX_D(MaxDim),  StepY_D(MaxDim)
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
      real :: Coef, Curv, State_V(nVar)
      !\
      ! Ray ID
      integer :: iRay
      !\
      ! Realistic emission
      !/
      real, dimension(nWave):: PlanckSpectrum_W, AbsorptionCoef_W
      !---------------
      iRay = iIndex_II(Ray_,iParticle)
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
         !\
         ! Get state vector at the mid point of the newly added ray segment
         !/
         State_V = -1
         !\
         ! Calculate absorption coefficient
         !/
         call user_material_properties(State_V, &
              OpacityEmissionOut_W = AbsorptionCoef_W,&
              PlanckOut_W = PlanckSpectrum_W)
         !\
         ! Realistic emission
         !/
         !Intensity_I(iRay) = Intensity_I(iRay)  &
         !+ ParabLen*AbsorptionCoef_W(1)*PlanckSpectrum_W(1)&
         !* FrequencySi_W(WaveLast_)
         Intensity_I(iRay) = Intensity_I(iRay)  &
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
         Intensity_I(iRay) = Intensity_I(iRay) &
              + State_VI(Ds_,iParticle)*&
              (Dens2DensCr**2)*(0.50 - Dens2DensCr)**2
      end if
    end subroutine advance_ray
    !====================
    subroutine get_new_step_size(StepSizeNew)
      real, intent(in) :: StepSizeNew
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
              max(2 - State_VI(Ds_,iParticle)/StepSizeNew, 1.0)*&
              min(StepSizeNew, State_VI(Ds_,iParticle))
      end if
    end subroutine get_new_step_size
    !====================
    subroutine get_local_density(&
         NElectronSi, GradNElectronSi_D, DeltaSNew, IsBody)
      !USES:
      use ModAdvance, ONLY: State_VGB
      use ModVarIndexes, ONLY: Rho_
      use ModCellGradient, ONLY: GradDensity_DGB => GradVar_DGB
      use ModGeometry,ONLY: CellSize_DB
      use ModPhysics, ONLY: No2Si_V, UnitN_
      real,    intent(out):: NElectronSi
      real,    intent(out):: GradNElectronSi_D(MaxDim)
      real,    intent(out):: DeltaSNew
      logical, optional, intent(out):: IsBody
      real :: Density
      real :: GradDensity_D(MaxDim)
      !Coordinates and block #
      real     :: Xyz_D(MaxDim)
      integer  :: iBlock
      ! interpolation data: number of cells, cell indices, weights
      integer:: nCell, iCell_II(0:nDim, 2**nDim)
      real   :: Weight_I(2**nDim)
      integer:: iCell ! loop variable
      integer:: i_D(MaxDim)
      character(len=*), parameter:: NameSub = 'get_local_density'
      !----------------------------------------------------------------------
      !Coordinates and block #
      Xyz_D   = State_VI(x_:z_, iParticle)
      iBlock  = iIndex_II(0,iParticle)
      call interpolate_grid_amr_gc(Xyz_D, iBlock, nCell, iCell_II, Weight_I,&
           IsBody)
      if(IsBody)then
         call mark_undefined(KindRay_,iParticle);RETURN
      end if
      Density = 0.0; GradDensity_D = 0.0; DeltaSNew = 0.0
      do iCell = 1, nCell
         i_D = 1
         i_D(1:nDim) = iCell_II(1:nDim, iCell)
         Density = Density + &
              State_VGB(Rho_,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
         GradDensity_D = GradDensity_D + &
              GradDensity_DGB(:,i_D(1),i_D(2),i_D(3),iBlock)*Weight_I(iCell)
         DeltaSNew = DeltaSNew + &
              Weight_I(iCell)*minval(CellSize_DB(:,iBlock))
      end do
      NElectronSi     = Density      *No2Si_V(UnitN_)
      GradNElectronSi_D = GradDensity_D*No2Si_V(UnitN_)
    end subroutine get_local_density
    !====================
  end subroutine ray_path
end module ModRadioWaveImage

