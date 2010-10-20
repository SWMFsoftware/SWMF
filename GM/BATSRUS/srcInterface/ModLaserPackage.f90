!=================================THE LASER PACKAGE STARTS HERE=================!
!Our intent is to make ray_path routine re-usable,
!which is a reason to put the package here.
!===============================================================================!
module ModLaserPulse
  implicit none
  !The laser pulse is assumed to have the linear raise front
  !and linear decay with the constant irradiance within the
  !pulse,  at  t_raise < t < tPulse - tDecay
  real:: tRaise = 1.0e-10 ![s]
  real:: tDecay = 1.0e-10 ![s]
  real:: tPulse = 1.1e-9  ![s]
  real:: IrradianceSi = 3.8e12 ![J/s]
contains
  real function irradiance_t(TimeSi)
    real,intent(in)::TimeSi
    !----------------------
    irradiance_t = IrradianceSi *&
         max(0.0, min(1.0,       &
         TimeSi/tRaise,          &
         (tPulse -TimeSi)/tDecay))
         
  end function irradiance_t
end module ModLaserPulse
!===========================
module ModBeams
  use ModPhysics,  ONLY: Si2No_V, UnitX_
  implicit none
  SAVE
  !Beam geometry: 'rz', '2d', '3d'
  character(LEN=2):: TypeBeam='unknown'
  !\
  ! Geometry of the 'rz'-beam:
  !/
  !
  !\     |
  !\\    |
  !\\\   |
  ! \\\  |
  !  \\\ |
  !   \\\|
  !    \\|_ _ _ _ _
  !     \|    I
  !      | This height is yCr   
  !______|____V____ axis of symmetry
  !      |
  !      |
  !     /|
  !    //|    In the presented case nRayPerBeam = 1 
  !   ///|    (half of all rays except for the Central Ray (CR))
  !  /// |    The intensity falls as exp(-r^2/rBeam^2),
  ! ///  |    where r is the transverse distance from the CR
  !///   |
  !//slope (approximately +60 deg in the presented case))
  !/-----|
  !
  !< This coordinate is xPlane

  integer:: nRayPerBeam

  integer:: nBeam = 0
  real   :: rBeamMuM = 1.0, xPlaneMum = -60.0

  real   :: rBeam, xPlane

  real   :: BeamParam_II(3,192)

  !Named indexes:
  integer, parameter:: SlopeDeg_ = 1, yCrMuM_ = 2, AmplitudeRel_ = 3

  integer:: nRayTotal = -1
  real, allocatable, dimension(:,:):: XyzRay_DI, SlopeRay_DI 
  real, allocatable, dimension(:)  :: Amplitude_I
 
contains
  subroutine get_rays
    !Convert MuM = 1e-6 m to dimensionless:

    rBeam  =  rBeamMuM * 1.0e-6 * Si2No_V(UnitX_)
    xPlane = xPlaneMuM * 1.0e-6 * Si2No_V(UnitX_)

    select case(TypeBeam)
    case('rz','RZ')
       call rz_beam_rays
    case('2d','2D','3d', '3D')
       call stop_mpi('TypeBeam='//TypeBeam//' is not yet implemented')
    case default
       call stop_mpi('Unknown TypeBeam='//TypeBeam)
    end select
    
    !Normalize amplitudes:
    Amplitude_I(1:nRayTotal) = Amplitude_I(1:nRayTotal)/&
         sum(Amplitude_I(1:nRayTotal))
  end subroutine get_rays
  !======================
  subroutine rz_beam_rays
    use ModGeometry, ONLY: TypeGeometry, y1
    use ModConst,    ONLY: cDegToRad 
    
    real:: CosTheta, SinTheta,  yCrCentral, yPlaneCentral, yPlane
    real:: rDistance, BeamAmplitude
    integer:: iRay, iBeam


    if(TypeGeometry/='rz')call CON_stop(&
         'Dont use TypeBeam='//TypeBeam//' with TypeGeometry='//TypeGeometry)

    !Allocation is excessive, nRayTotal is not yet known:

    allocate(  XyzRay_DI(3, nBeam * (2 * nRayPerBeam +1)))
    allocate(SLopeRay_DI(3, nBeam * (2 * nRayPerBeam +1)))
    allocate(Amplitude_I(nBeam * (2 * nRayPerBeam +1)))

    nRayTotal = 0
    do iBeam = 1, nBeam
       cosTheta =  cos(cDegToRad * BeamParam_II(SlopeDeg_, iBeam))

       !Positive direction of the slope is taken for the beams
       !converging to the axis:

       sinTheta = -sin(cDegToRad * BeamParam_II(SlopeDeg_, iBeam))

       !Transform yCrMuM to dimensionless
       yCrCentral =  BeamParam_II(yCrMuM_, iBeam) &
            * 1.0e-6 * Si2No_V(UnitX_)

       yPlaneCentral = yCrCentral + xPlane*SinTheta/CosTheta

       BeamAmplitude = BeamParam_II(AmplitudeRel_,iBeam)
 
       do iRay = -nRayPerBeam, nRayPerBeam
          !We neglect exp(-2.25)\approx0.1 and chose the beam margin
          !to be at 1.5 rBeam from the central ray:

           rDistance = iRay * rBeam * 1.5 /nRayPerBeam 
           yPlane = yPlaneCentral + rDistance/CosTheta

           !Do not include rays with the starting point 
           !being otside the computational domain:
           if(abs(yPlane) >= y1)CYCLE
           nRayTotal = nRayTotal +1
           
           Amplitude_I(nRayTotal) = BeamAmplitude * &
                exp(-(rDistance/rBeam)**2) *   &
                abs( yCrCentral + rDistance/CosTheta)

           if(yPlane > 0)then

              XyzRay_DI(:, nRayTotal) = (/xPlane, yPlane, 0.0/)
              SlopeRay_DI(:, nRayTotal) = &
                   (/CosTheta, SinTheta, 0.0/)

           else

              XyzRay_DI(:, nRayTotal) = (/xPlane, -yPlane, 0.0/)
              SlopeRay_DI(:, nRayTotal) = &
                   (/CosTheta, -SinTheta, 0.0/)

           end if
       end do
    end do
  end subroutine rz_beam_rays
end module ModBeams

subroutine read_laser_pulse_param
  use ModReadparam
  use ModLaserPulse
  use ModBeams
  use ModMain, ONLY: UseLaserPackage
  implicit none

  character (len=100) :: NameCommand

  !-------------
  !Usage
  !#LASERPULSE
  !T                UseLaserPackage
  !3.8e12           IrradianceSI
  !1.1e-9           tPulse
  !1.0e-10          tRaise
  !1.0e-10          tDecay
  !
  !#RZBEAMS or #BEAMS2D or #BEAMS3D
  !30               nRayPerBeam
  !438.0            rBeamMuM            
  !-100.0           xPlaneMuM
  !parallel         TypeConvergence
  !
  !#BEAM
  !10.0              SlopeDeg
  !0.0              yCrMuM
  !1.0              AmplitudeRel
  !
  !#BEAM
  !20.0             SlopeDeg
  !200.0            yCrMuM
  !0.7              AmplitudeRel
  !
  !#BEAM
  !30.0             SlopeDeg
  !400.0            yCrMuM
  !0.7              AmplitudeRel
  !
  !#END_LASERPULSE                                
  !
  !-------------
  call read_var('UseLaserPackage', UseLaserPackage)
  if(.not.UseLaserPackage) return
  call read_var('IrradianceSI'   , IrradianceSI   )
  call read_var('tPulse'         , tPulse         )
  call read_var('tRaise'         , tRaise         )
  call read_var('tDecay'         , tDecay         )
  do
     if(.not.read_line() ) EXIT
     if(.not.read_command(NameCommand)) CYCLE
     select case(NameCommand)
     case('#END_LASERPULSE')
        return
     case('#RZBEAMS')
        TypeBeam = 'rz'
        call read_var('nRayPerBeam', nRayPerBeam)
        call read_var('rBeamMuM'   , rBeamMuM   )
        call read_var('xPlaneMuM'  , xPlaneMuM  )
     case('#BEAM')
        nBeam = nBeam +1
        call read_var('SlopeDeg', &
             BeamParam_II(SlopeDeg_, nBeam))
        call read_var('yCrMuM', &
             BeamParam_II(yCrMuM_, nBeam))
        call read_var('AmplitudeRel', &
             BeamParam_II(AmplitudeRel_, nBeam))
     case default
        call stop_mpi(&
             'ERROR in read_laser_pulse_param: unknown command='//NameCommand)
     end select
  end do
end subroutine read_laser_pulse_param
!===========================
module ModLaserPackage
  use CON_global_message_pass
  use ModAbsorption, ONLY: DensityCrSi, NameMask,&
       LineGrid, MhGrid, get_density_and_absorption
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)
  
  use ModBeams, ONLY: nRayTotal, XyzRay_DI, SlopeRay_DI, Amplitude_I
  use ModRadioWaveRaytracing, ONLY: ray_path
  use ModDensityAndGradient, ONLY: EnergyDeposition_I
  implicit none
  
  real, allocatable, save:: SourceE_CB(:,:,:,:)
  real, dimension(:,:), pointer:: Position_DI
  logical,dimension(:), pointer:: DoRay_I
  real, allocatable, save:: Slope_DI(:,:), Intensity_I(:), DeltaS_I(:)
  logical, allocatable,save:: Unused_I(:), IsBehindCr_I(:)
  logical, save:: DoInit = .true.
  integer, save:: nRay = -1
  
  real:: DeltaS = 1.0, Tolerance = 0.1

  type(RouterType),save::Router
  
  !------------
contains
  !================================
  subroutine init_laser_package
    use ModSize, ONLY: MaxBlock, nI, nJ, nK
    use ModDensityAndGradient, ONLY: Density_I, GradDensity_DI,DeltaSNew_I, &
         NameVector
    use ModAbsorption, ONLY: AbsorptionCoeff_I
    use ModMain,ONLY:NameThisComp
    use CON_global_vector, ONLY: allocate_vector, associate_with_global_vector
    use CON_global_vector, ONLY: allocate_mask, associate_with_global_mask
    use ModBeams, ONLY: get_rays
    !--------------------------------------
    DoInit = .false.
    call get_rays
    nRay = nRayTotal
    
    !Allocate arrays for rays
    allocate(SourceE_CB(nI,nJ,nK,MaxBlock))

    if(nRay >0)then
       allocate(Density_I(nRay), GradDensity_DI(3,nRay),DeltaSNew_I(nRay))
       allocate(AbsorptionCoeff_I(nRay))
       allocate(EnergyDeposition_I(nRay))
    end if

    NameVector=NameThisComp//'_Rays_DI'
    call allocate_vector(NameVector, 3, nRayTotal)
    call associate_with_global_vector(Position_DI, NameVector)
 

    NameMask = NameThisComp//'_UsedRay'
    call allocate_mask(NameMask, nRayTotal)
    call associate_with_global_mask(DoRay_I, NameMask)

    allocate(Slope_DI(3, nRayTotal))
    allocate(Intensity_I(nRayTotal), DeltaS_I(nRayTotal))
    allocate(Unused_I(nRayTotal), IsBehindCr_I(nRayTotal))

    !Initialize all arrays:
    Position_DI(:, 1:nRayTotal) = XyzRay_DI(:,   1:nRayTotal)
    Slope_DI(:,    1:nRayTotal) = SlopeRay_DI(:, 1:nRayTotal)
    Intensity_I(   1:nRayTotal) = Amplitude_I(   1:nRayTotal)
    DoRay_I = .true.

    Unused_I = .false.; IsBehindCr_I = .false.

    DeltaS_I = DeltaS
    
    call get_density_and_absorption(nRayTotal)
    
    call init_router(&
         GridDescriptorSource = LineGrid, &
         GridDescriptorTarget = MhGrid, &
         Router=Router,                 &
         nIndexesSource = 1)
    

  end subroutine init_laser_package
  !================================
  subroutine get_impl_energy_source
    use ModAbsorption, ONLY: NameMask
    use ModDensityAndGradient, ONLY: NameVector
        

    !--------------------------
    if(DoInit) call init_laser_package
    SourceE_CB(:,:,:,:) = 0.0

    !Initialize all arrays:
    Position_DI(:, 1:nRayTotal) = XyzRay_DI(:,   1:nRayTotal)
    Slope_DI(:,    1:nRayTotal) = SlopeRay_DI(:, 1:nRayTotal)
    Intensity_I(   1:nRayTotal) = Amplitude_I(   1:nRayTotal)
    DoRay_I = .true.

    Unused_I = .false.; IsBehindCr_I = .false.
    DeltaS_I = DeltaS

    do while(.not.all(Unused_I))
       EnergyDeposition_I=0.0
       !Propagate each of rays through the distance of DeltaS
       call ray_path(get_density_and_absorption, nRay, Unused_I, Slope_DI, &
            DeltaS_I, Tolerance, DensityCrSi, Intensity_I, IsBehindCr_I)
       DoRay_I = (.not.Unused_I).or.IsBehindCr_I

       !Save EnergyDeposition_I to SourceE_CB
       call construct_router_from_source(&
            GridDescriptorSource = LineGrid, &
            GridDescriptorTarget = MhGrid, &
            Router=Router,                 &
            NameMappingVector=NameVector,  &
            NameMask=NameMask,             &
            interpolate=interpolation_fix_reschange)
   
       call global_message_pass(Router=Router,&
            nVar=1,&    !Energy deposition only
            fill_buffer=get_energy_deposition,&
            apply_buffer=put_energy_deposition)
       DoRay_I = (.not.Unused_I)
    end do

  end subroutine get_impl_energy_source
  !=============================
  subroutine get_energy_deposition(&
       nPartial,iGetStart,Get,W,Buff_V,nVar)
    use CON_router
    !INPUT ARGUMENTS:
    integer,intent(in)::nPartial,iGetStart,nVar
    type(IndexPtrType),intent(in)::Get
    type(WeightPtrType),intent(in)::W
    real,dimension(nVar),intent(out)::Buff_V

    integer::iCell

    iCell=Get%iCB_II(1,iGetStart)
    Buff_V(1) = EnergyDeposition_I(iCell)
   
  end subroutine get_energy_deposition

  !====================================================================

  subroutine put_energy_deposition(nPartial,&
       iPutStart,&
       Put,&
       W,&
       DoAdd,&
       Buff_V,nVar)
    implicit none
    integer,intent(in)::nPartial,iPutStart,nVar
    type(IndexPtrType),intent(in)::Put
    type(WeightPtrType),intent(in)::W
    logical,intent(in)::DoAdd
    real,dimension(nVar),intent(in)::Buff_V
    integer::iPut, i, j, k, iBlock
    real :: Weight

    do iPut = iPutStart, iPutStart - 1 + nPartial
       i      = Put%iCB_II(1,iPut)
       j      = Put%iCB_II(2,iPut)
       k      = Put%iCB_II(3,iPut)
       iBlock = Put%iCB_II(4,iPut)
       Weight = W%Weight_I(iPut)
       
       SourceE_CB(i, j, k, iBlock) = Buff_V(1) * Weight
    end do

  end subroutine put_energy_deposition
end module ModLaserPackage

subroutine add_laser_energy_deposition
  !This routine should add the laser energy deposition to the "Explicit residual"
  !within the semi-implicit scheme. Effectively, in this way the heat conduction
  !equation, which otherwise is solved within the semi-implicit scheme WITH ZERO
  !right-hand-side, is solved with the right-hand-side including the laser energy
  !deposition
  use ModSize, ONLY: nI, nJ, nK
  use ModImplicit, ONLY: ResExpl_VCB, iTeImpl, nImplBlk, impl2iblk
  !Here ResExpl_VCB is the energy source array, to which the laser energy deposition 
  !should be added. This array is the dimensionless form of the energy source in each
  !control volume, divided by the time step, Delta t. The dimensional characteristic
  !of the laser emission with the meaning of the energy delivired to the plasma is
  !the irradiance (power). Therefore, the pointwise sources of energy calculated by
  !this subroutine are multipled by :
  !LaserIrradiance*Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  use ModPhysics, ONLY: Si2No_V, UnitEnergydens_, UnitX_, UnitT_
  use ModLaserPackage, ONLY: get_impl_energy_source, SourceE_CB
  use ModMain, ONLY: Time_Simulation
  use ModLaserPulse, ONLY: irradiance_t
  implicit none
  real:: Irradiance
  integer :: iImplBlock, iBlock, i, j, k, iVar

  character(len=*), parameter:: NameSub = 'add_laser_energy_deposition'
  !----------------
  Irradiance = irradiance_t(Time_Simulation) * &
       Si2No_V(UnitEnergydens_)*Si2No_V(UnitX_)**3/Si2No_V(UnitT_)

  call get_impl_energy_source

  do iImplBlock = 1, nImplBLK
     iBlock = impl2iBLK(iImplBlock)
     do k = 1, nK; do j = 1, nJ; do i = 1, nI
        ResExpl_VCB(iTeImpl,i,j,k,iImplBlock) = &
             ResExpl_VCB(iTeImpl,i,j,k,iImplBlock) + &
             SourceE_CB(i,j,k,iBlock) * Irradiance
     end do; end do; end do
  end do
  
end subroutine add_laser_energy_deposition
!==========================================
