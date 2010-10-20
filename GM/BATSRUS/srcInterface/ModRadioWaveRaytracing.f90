!^CFG COPYRIGHT UM
module ModRadioWaveRaytracing

  use ModNumConst
  use ModCoordTransform
  use CON_global_vector, ONLY: associate_with_global_vector
  use ModDensityAndGradient, ONLY: NameVector, &
       GradDensity_DI, Density_I, DeltaSNew_I, EnergyDeposition_I
  use ModPhysics, ONLY : No2Si_V, UnitRho_, UnitX_
  use ModProcMH, ONLY: iProc
  use ModMain,ONLY: UseLaserPackage
  use ModAbsorption, ONLY: AbsorptionCoeff_I
  use ModGeometry, ONLY: XyzMin_D, XyzMax_D, nDim
  use ModMain,     ONLY: TypeBC_I

  implicit none

  public :: ray_path
  Logical:: UseParallel = .false.


contains !=========================================================

  subroutine ray_path(get_plasma_density, nRay, UnusedRay_I, Slope_DI, &
       DeltaS_I, ToleranceInit, DensityCr, Intensity_I, IsBehindCr_I)

    !
    !   The subroutine ray_path() makes raytracing and emissivity integration
    ! along ray paths. It works on a group of rays (a beam). Each ray is
    ! specified by its Cartesian Position_DI and its direction cosines 
    ! Slope_DI. The subroutine calculates new Position_DI, which is 
    ! DeltaS_I away from the initial position. It calculates the intensity
    ! along the step as emissivity multiplied by DeltaS and adds this to the
    ! Intensity_I. Thus, after series of calls to ray_path(), the Intensity_I
    ! contains the result of integration along the paths of every ray. 
    !   The electromagnetic rays are refracted in plasma due to its non-zero
    ! refractive index, which is the square root of the dielectric
    ! permittivity \epsilon. The \epsilon, in turn, is a function of the 
    ! plasma density. The external subprogram, get_plasma_density(),
    ! provides the plasma Density_I along with its gradient, GradDensity_DI,
    ! at the Position_DI. The \epsilon can be calculated as
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
    !     recommended step size DeltaSNew_I and asserts the IsBehindCr_I (sets
    !     it to .true.) for a ray should it (illegally) penetrate into the
    !     region with "negative" dielectric permittivity.
    ! nRay:           number of rays being processed.
    ! UnusedRay_I:  the caller program can use this logical array to stop
    !     processing of any individual ray when it already finished travelling
    !     inside of a specified part of space. Set the corresponding element
    !     of UnusedRay_I to .true. to leave it unprocessed during the
    !     subsequent calls to ray_path(). Before the first call to ray_path()
    !     all the elements of UnusedRay_I should be set to .false.; then all
    !     the rays will be processed.
    ! Position_DI:    Cartesian position vectors of the rays.
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
    ! IsBehindCr_I:      the .true. elements of this logical array indicate that
    !     the corresponding rays penetrated into the "prohibited" region of
    !     space with the plasma density above its critical value. Normally, it
    !     should never happen. However, in case the algorithm made such an
    !     error, the flagged rays should be considered as "bad" and thrown
    !     away from the resultant Intensity_I. Set all the elements of
    !     IsBehindCr_I to .false. before calling ray_path() for the first time.
    ! 
    ! Written by Leonid Benkevitch.
    !
    !

    implicit none

    interface
       subroutine get_plasma_density(nRay)
         implicit none
         integer,intent(in) :: nRay
       end subroutine get_plasma_density
    end interface

    integer, intent(in) :: nRay   ! # of pixels in the raster
    real,    intent(inout), dimension(3,nRay) :: Slope_DI
    real,    intent(inout), dimension(nRay) :: Intensity_I, DeltaS_I
    real,    intent(in) ::  ToleranceInit, DensityCr

    ! .true. if a ray is OK, .false. otherwise:
    logical, intent(inout), dimension(nRay):: IsBehindCr_I

    ! a ray is excluded from processing if it is .true. 
    logical, intent(inout), dimension(nRay) :: UnusedRay_I   

    logical, save:: IsNewEntry = .true.
   
    integer, parameter :: nSplitDeltaS = 2

    real,    save, dimension(3)       :: Slope1_D, Omega_D
    real,    save, dimension(3)       :: ProjSlopeOnMinusGradEps_D
    real,    save, dimension(3)       :: StepX_D, StepY_D, RelGradRefrInx_D 
    real,    save, dimension(3)       :: GradDielPerm_D, PositionHalfBack_D

    real, save, pointer, dimension(:,:) :: Position_DI
    real, save, pointer, dimension(:)   :: DistanceToCritSurf_I

    ! GentleRay_I: .true. for shallow rays; is set to .false. for steep rays.
    logical, save, pointer, dimension(:)   :: GentleRay_I       
    

    real :: HalfDeltaS                        ! DeltaS halved
    real :: DielPerm, DielPermHalfBack, Dens2DensCr, Dens2DensCr1, Dens2DensCr2
    real :: Coef, Curv, Curv1
    real :: LCosAl ! L*cos(Alpha),   
    ! where L is inverse grad of \epsilon, Alpha is the incidence angle
    real :: GradDielPermSqr, GradEpsDotSlope
    real :: ParabLen, GradDielPerm
    real, save :: &
         Tolerance=0.1, ToleranceSqr=1.0e-2, DensityCrInv=1, AbsoluteMinimumStep=1
    
    integer, save :: nCall=0
    integer:: iRay

    character(LEN=20):: TypeBoundaryDown_D(nDim), TypeBoundaryUp_D(nDim)
    !---------------

    if (IsNewEntry) then
       IsNewEntry = .false.
       call associate_with_global_vector(Position_DI,NameVector)
       allocate(DistanceToCritSurf_I(nRay))
       DistanceToCritSurf_I = cZero
       allocate(GentleRay_I(nRay))
       GentleRay_I = .true.
       DensityCrInv = cOne/DensityCr

       !  minimum ten points between a vacuum and a critical surface and
       !  minimum 10 points over 1 rad of the curvature
       Tolerance = min(ToleranceInit,0.1)  

       
       ToleranceSqr = Tolerance**2 
       
       ! One (ten-thousandth) hundredth of average step 
       AbsoluteMinimumStep = 1e-2*sum(DeltaS_I)/nRay 
       TypeBoundaryDown_D = TypeBc_I(1:2*nDim-1:2)
       TypeBoundaryUp_D   = TypeBc_I(2:2*nDim:2)
    end if

    nCall = nCall + 1
   

    !Start the predictor step of the GIRARD scheme: 
    do iRay = 1, nRay
       if (UnusedRay_I(iRay)) CYCLE  ! Do not process the rays that are done
       Position_DI(:,iRay) = Position_DI(:,iRay) + &
            Slope_DI(:,iRay)*cHalf*DeltaS_I(iRay) 
       ! Now Position_DI moved by 1/2 DeltaS !!!
       call check_bc
    end do

    !
    ! The following routine collects the density values at the Position_DI
    ! from all relevant processors
    !

    call get_plasma_density(nRay)
    
    
    !In making radio images this is a bad pixel which should be further processed 
    IsBehindCr_I = Density_I>= DensityCr.and..not.UnusedRay_I
    UnusedRay_I = UnusedRay_I .or. IsBehindCr_I ! "bad rays" are done
    
    !In solving the laser deposition energy the total remnant energy 
    !should be deposited if the ray pamatrates through the critical surface
    if(UseLaserPackage)then
       where( IsBehindCr_I)
          EnergyDeposition_I = Intensity_I
          Intensity_I = 0.0
       endwhere
    end if

  

    do iRay = 1, nRay

       if (UnusedRay_I(iRay)) CYCLE  ! Do not process the rays that are done


       HalfDeltaS = cHalf*DeltaS_I(iRay)
       ! Original Position (at an integer point):
       PositionHalfBack_D = Position_DI(:,iRay) - Slope_DI(:,iRay)*HalfDeltaS
       Dens2DensCr = Density_I(iRay)*DensityCrInv

       DielPerm = cOne - Dens2DensCr
       GradDielPerm_D = -GradDensity_DI(:,iRay)*DensityCrInv
       GradDielPermSqr = sum(GradDielPerm_D**2)

       GradEpsDotSlope = sum(GradDielPerm_D*Slope_DI(:,iRay))

       DielPermHalfBack = DielPerm  - GradEpsDotSlope*HalfDeltaS

       Curv = (cHalf*HalfDeltaS/DielPermHalfBack)**2* &
            (GradDielPermSqr - GradEpsDotSlope**2)

       if (GentleRay_I(iRay)) then

          !
          ! Check if the trajectory curvature is too sharp
          ! to meet the Tolerance. If so, reduce the DeltaS 
          ! step for the iRay-th ray and leave it until
          ! the next call.
          !
          !

          if (Curv .ge. ToleranceSqr) then
             DeltaS_I(iRay) = DeltaS_I(iRay)/(2*sqrt(Curv/ToleranceSqr))
             Position_DI(:,iRay) = PositionHalfBack_D
             if(UseLaserPackage)EnergyDeposition_I(iRay) = 0.0
             CYCLE
          end if

          !
          ! Test if some of the next points can get into the prohibited
          ! part of space with "negative" dielectric permittivity
          !

          if (3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack) then
             !
             ! Mark the ray as steep; 
             ! memorize the distance to the critical surface;
             ! reduce step
             !
             GentleRay_I(iRay) = .false.
             DistanceToCritSurf_I(iRay) = DielPermHalfBack  &
                  /sqrt(GradDielPermSqr)
             DeltaS_I(iRay) =  cHalf*Tolerance*DistanceToCritSurf_I(iRay)
             Position_DI(:,iRay) = PositionHalfBack_D
             if(UseLaserPackage)EnergyDeposition_I(iRay) = 0.0
             CYCLE
          end if

       end if ! GentleRay_I

       !
       ! Either switch to opposite pranch of parabola
       ! or make a Boris step
       !

       if ((3*GradEpsDotSlope*HalfDeltaS <= -DielPermHalfBack) &
            .or. (DeltaS_I(iRay) .lt. AbsoluteMinimumStep)) then

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
          ! The parabolic reflection is just a replacement of the Position_DI
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

          Position_DI(:,iRay) = PositionHalfBack_D + StepY_D

          !
          ! Step_X is in the direction of - grad \epsilon,
          ! whose vector is of the length of 1/L.
          ! The length of Step_X is cos^2(\alpha)*L*DielPermHalfBack 
          ! Thus,
          !

          StepX_D = (DielPermHalfBack*LCosAl**2)*GradDielPerm_D


          ParabLen = sqrt(sum((2*StepX_D)**2) + sum(StepY_D**2))
          if(.not.UseLaserPackage)then
             Intensity_I(iRay) = Intensity_I(iRay)  &
                  + ParabLen*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2
          else
             EnergyDeposition_I(iRay) = Intensity_I(iRay) * &
                  ParabLen * AbsorptionCoeff_I(iRay)
             Intensity_I(iRay) = Intensity_I(iRay) * &
                  (1 - ParabLen * AbsorptionCoeff_I(iRay)) 
          end if

       else 

          ! Make a step using Boris' algorithm

          Coef=cHalf*HalfDeltaS/(cOne - Dens2DensCr)
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

          Position_DI(:,iRay) = Position_DI(:,iRay) &
               + Slope_DI(:,iRay)*HalfDeltaS
          if(UseLaserPackage)then
             EnergyDeposition_I(iRay) = Intensity_I(iRay) *  &
                  DeltaS_I(iRay) * AbsorptionCoeff_I(iRay)
             Intensity_I(iRay) = Intensity_I(iRay) * &
                  (1 - DeltaS_I(iRay) * AbsorptionCoeff_I(iRay)) 
          else
             Intensity_I(iRay) = Intensity_I(iRay) &
                  + DeltaS_I(iRay)*(Dens2DensCr**2)*(cHalf - Dens2DensCr)**2
          end if
       end if

       call check_bc

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
       if (GentleRay_I(iRay)) then
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
             GentleRay_I(iRay) = .true.
             DeltaS_I(iRay) = max(2 - DeltaS_I(iRay)/DeltaSNew_I(iRay), 1.0)*&
                  min(DeltaSNew_I(iRay), DeltaS_I(iRay))
          end if
       end if
    end do
  contains
    !==================
    subroutine check_bc
      integer::iDim
      do iDim = 1, nDim
         if(Position_DI(iDim, iRay) < XyzMin_D(iDim))then
            if(TypeBoundaryDown_D(iDim)=='reflect')then
               Position_DI(iDim, iRay) = &
                    2 * XyzMin_D(iDim) -  Position_DI(iDim, iRay)
               Slope_DI(iDim, iRay) = - Slope_DI(iDim, iRay)
            else
               UnusedRay_I(iRay) = .true.
               return
            end if
         else if(Position_DI(iDim, iRay) > XyzMax_D(iDim))then
            if(TypeBoundaryUp_D(iDim)=='reflect')then
               Position_DI(iDim, iRay) = &
                    2 * XyzMax_D(iDim) -  Position_DI(iDim, iRay)
               Slope_DI(iDim, iRay) = - Slope_DI(iDim, iRay)
            else
               UnusedRay_I(iRay) = .true.
               return
            end if
         end if
      end do
    end subroutine check_bc
    !====================
  end subroutine ray_path

end module ModRadioWaveRaytracing


