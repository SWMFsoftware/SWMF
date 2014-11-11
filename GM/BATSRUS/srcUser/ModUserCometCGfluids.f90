!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!#NOTPUBLIC  email:zghuang@umich.edu  expires:12/31/2099
!========================================================================
module ModUser

  use ModUserEmpty,               &
       IMPLEMENTED1  => user_set_boundary_cells,         &
       IMPLEMENTED2  => user_read_inputs,                &
       IMPLEMENTED3  => user_init_session,               &
       IMPLEMENTED4  => user_set_face_boundary,          &
       IMPLEMENTED5  => user_calc_sources,               &
       IMPLEMENTED6  => user_update_states,              &
       IMPLEMENTED7  => user_set_resistivity,            &
       IMPLEMENTED8  => user_material_properties,        &
       IMPLEMENTED9  => user_init_point_implicit,        &
       IMPLEMENTED10 => user_set_plot_var,               &
       IMPLEMENTED11 => user_set_ICs,                    &
       IMPLEMENTED12 => user_get_log_var,                &
       IMPLEMENTED13 => user_set_cell_boundary,          &
       IMPLEMENTED14 => user_amr_criteria,               &
       IMPLEMENTED15 => user_initial_perturbation

  use ModSize
  use ModProcMH,    ONLY: iProc
  use ModNumConst,  ONLY: cPi
  use ModAdvance,   ONLY: Pe_, UseElectronPressure
  use ModMultiFluid

  include 'user_module.h' !list of public methods

  !\
  ! Here you must define a user routine Version number and a 
  ! descriptive string.
  !/
  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: NameUserModule = &
       'CG Comet, G. Toth & H. Zhenguang, 2014'

  character (len=100) :: NameShapeFile

  ! Use the CG shape or not. If not, then use a spherical body.
  logical :: DoUseCGShape = .true.
  real    :: rSphericalBodySi = 2.0e3
  real    :: rSphericalBody

  integer:: nTriangle
  real, allocatable:: XyzTriangle_DII(:,:,:), Normal_DI(:,:)
  real :: rMinShape = 0.0, rMaxShape = 0.0

  ! Position of the sun at the start time
  real :: LatSun=43.5, LonSun=0.0
  real :: NormalSun_D(3)

  ! Rotation of the comet (changes the direction of the Sun)
  real:: RotationCometHour = 12.0

  ! Angular velocity
  real:: OmegaComet

  ! Maximum change in longitude before updating the boundary conditions
  real:: AngleUpdateDeg = 10.0

  ! The time between updates
  real:: DtUpdateSi

  ! minimum and maximum temperature
  real :: TempCometMinDim, TempCometMaxDim, TempCometMin, TempCometMax

  ! Temperature at 75.5 degree to determin the temperature slope
  real :: TempComet75Dim, TempComet75

  ! Minimum and maximum production rate
  real :: ProductionRateMaxSi, ProductionRateMinSi
  real :: ProductionRateMax, ProductionRateMin

  ! Maximum solar zenith angle for dayside production rate
  real :: SolarAngleMaxDim, SolarAngleMax

  ! Parameters for y=ax+b to mimic the production rate and 
  ! temperature distribution
  real :: SlopeProduction, bProduction, SlopeTemp, bTemp

  ! Constant parameters to calculate uNormal and temperature from TempCometLocal
  real :: TempToUnormal
  real :: TempToPressure
  real :: NumdensToRho

  ! Inner boundary condition for ions
  character (len=10) :: TypeBodyBC = 'default'
  logical :: UseSwBC        = .false.
  logical :: UseReflectedBC = .false.

  ! FaceCoordsTest_D
  real :: FaceCoordsX=0.0, FaceCoordsY=0.0, FaceCoordsZ=0.0
  real :: FaceCoordsTest_D(3) = (/0.0, 0.0, 0.0/)

  ! Last step and time the inner boundary values were saved for each block
  integer :: nStepSave_B(MaxBlock) = -100
  real    :: TimeSimulationSave_B(MaxBlock) = -1e30
  integer :: nStepSaveCalcRates_B(MaxBlock) = -100
  integer :: nStepPritSetFace = -100

  !++++++++++++++++++++++++++++++++++++++++++++++
  ! From Martin & Xianzhe 2014
  !++++++++++++++++++++++++++++++++++++++++++++++
  integer, parameter, public :: nNeutral = 1
  integer, parameter :: Neu1_  =  1

  !! Ion species names
  integer, parameter :: SW_   =  1
  integer, parameter :: H2Op_ =  2

  real :: Qprod = 1e22
  real :: Tmin, rHelio, vHI, uHaser
  real, dimension(MaxI-MinI+1,MaxJ-MinJ+1,MaxK-MinK+1,nBLK) :: ne20eV_GB = 0.

  character (len=6), parameter, public :: NameNeutral_I(nNeutral) = &
       (/ 'Neu1  ' /)
contains
  !============================================================================
  subroutine user_read_inputs

    use ModMain, ONLY: UseUserInitSession, UseExtraBoundary
    use ModReadParam

    character (len=100) :: NameCommand

    character(len=*), parameter:: NameSub = 'user_read_inputs'
    !-------------------------------------------------------------------------

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#SHAPEFILE")
          call read_var('NameShapeFile' ,NameShapeFile)
       case("#USECGSHAPE")
          call read_var('DoUseCGShape',     DoUseCGShape)
          call read_var('rSphericalBodySi', rSphericalBodySi)
       case("#SUNDIRECTION")
          call read_var('LatSun', LatSun)
          call read_var('LonSun', LonSun)
       case("#COMETSTATE")
          call read_var('ProductionRateMinSi', ProductionRateMinSi)
          call read_var('ProductionRateMaxSi', ProductionRateMaxSi)
          call read_var('SolarAngleMaxDim',     SolarAngleMaxDim)
          call read_var('TempCometMinDim',      TempCometMinDim)
          call read_var('TempCometMaxDim',      TempCometMaxDim)
          call read_var('TempComet75',       TempComet75Dim)
       case("#IONIZATIONPARAM")
          call read_var('rHelio', rHelio)   !! Heliocentric distance [AU]
          call read_var('vHI', vHI)         !! Ionization frequency for cometary heavy ions (1/lifetime of cometary heavy neutrals)
          call read_var('Tmin', Tmin)       !! Minimum ion temperature (enforced in update states)
          call read_var('Qprod', Qprod)     !! Total production rate
          call read_var('uHaser', uHaser)   !! velocity for Haser model used in preset conditions
       case("#COMETROTATION")
          call read_var('RotationCometHour', RotationCometHour)
          call read_var('AngleUpdateDeg', AngleUpdateDeg)
       case("#BODYBC")
          call read_var('TypeBodyBC', TypeBodyBC)
          select case(TypeBodyBC)
             case('solarwind')
                UseSwBC        = .true.
                UseReflectedBC = .false.
             case('reflected')
                UseSwBC        = .false.
                UseReflectedBC = .true.
             case('default')
                UseSwBC        = .false.
                UseReflectedBC = .false.
             case default
                if(iProc==0) call stop_mpi( &
               NameSub//' invalid body type='//trim(NameCommand))
             end select
       case("#TESTFACECOORDS")
          call read_var('FaceCoordsX', FaceCoordsX)
          call read_var('FaceCoordsY', FaceCoordsY)
          call read_var('FaceCoordsZ', FaceCoordsZ)
       case('#USERINPUTEND')
          EXIT
       case default
          if(iProc==0) call stop_mpi( &
               NameSub//' invalid command='//trim(NameCommand))
       end select
    end do

    ! UseUserInitSession = .true.
    ! UseExtraBoundary   = .true.

  end subroutine user_read_inputs
  !===========================================================================
  subroutine user_init_session

    ! Read shape file and convert units

    use ModMain, ONLY: Time_Simulation, n_step
    use ModPhysics, ONLY: Io2No_V, Si2No_V, No2Si_V, &
         UnitU_, UnitTemperature_, UnitT_, &
         UnitN_, UnitX_, gm1
    use ModNumConst, ONLY: cTwoPi, cDegToRad
    use ModCoordTransform, ONLY: dir_to_xyz
    use ModConst, ONLY: cBoltzmann, cAtomicMass
    use ModVarIndexes, ONLY: MassFluid_I
    use ModBlockData, ONLY: MaxBlockData
    use ModIO, ONLY: restart

    character(len=*), parameter :: NameSub='user_init_session'

    !------------------------------------------------------------------------

    if (DoUseCGShape) then
       ! We need to have unit conversions before reading the shape file 
       ! which contains everything in SI units
       call read_shape_file
       if (iProc ==0) &
            write(*,*) NameSub, ': reading CG shape file.'
    else
       if (iProc ==0) &
            write(*,*) NameSub, ': using spherical body: rSphericalBodySi =', &
            rSphericalBodySi
    end if

    rSphericalBody = rSphericalBodySi*Si2No_V(UnitX_)
    ProductionRateMax = &
         ProductionRateMaxSi / Si2No_V(UnitX_)**2 / Io2No_V(UnitT_)
    ProductionRateMin = &
         ProductionRateMinSi / Si2No_V(UnitX_)**2 / Io2No_V(UnitT_)

    SolarAngleMax = SolarAngleMaxDim * cDegToRad
    TempCometMin = TempCometMinDim * Io2No_V(UnitTemperature_)
    TempCometMax = TempCometMaxDim * Io2No_V(UnitTemperature_)
    TempComet75  = TempComet75Dim  * Io2No_V(UnitTemperature_)

    ! From Huebner & Markiewicz 2000, eq 8:
    ! uNormal = sqrt( pi kT / 2*m ), so TempToUnormal = sqrt( pi k/ 2m )
    ! and also unit conversions of temperature to SI, and velocity from SI
    TempToUnormal = sqrt(cPi*cBoltzmann/(2*MassFluid_I(nFluid)*cAtomicMass) * &
         No2Si_V(UnitTemperature_))*Si2No_V(UnitU_)

    ! From Huebner & Markiewicz, 2000 eq. 13)
    ! T' = T*((8 + 2 f_rv - pi)/(2(f_rv + 3)) = T(2f + 2 - pi)/(2f) 
    !    = T[1 - (pi-2)/4*(gamma - 1)]
    ! so TempToPressure = [1 - (pi-2)/4*(gamma - 1)]/Mass
    TempToPressure = (1 - gm1*0.25*(cPi-2))/MassFluid_I(nFluid)

    ! Number density calculated as production rate/velocity
    ! which is in units of 1/(length cubed).
    ! This has to be multiplied with mass of molecule which is in amu,
    ! but mass density is measured in amu/cm^3 so the 
    ! normalized number density is in 1/cm^3.
    NumdensToRho = MassFluid_I(nFluid)/(No2Si_V(UnitN_)*No2Si_V(UnitX_)**3)

    ! Calculate the parameters for production rate (y = a*cos(theta)+b)
    SlopeProduction = &
         (ProductionRateMax - ProductionRateMin) / (1-cos(SolarAngleMax))
    bProduction     = &
         (ProductionRateMin - ProductionRateMax*cos(SolarAngleMax)) / &
         (1-cos(SolarAngleMax))

    ! Calculate the parameters for temperature (y = a/cos(theta)+b)
    SlopeTemp = (TempCometMax - TempComet75)/(1 - 1/cos(75.5*cDegToRad))
    bTemp     = TempCometMax - SlopeTemp

    ! Angular velocity of the comet
    OmegaComet = cTwoPi / (RotationCometHour*3600 * Si2No_V(UnitT_))

    ! Frequency of boundary condition updates
    DtUpdateSi = AngleUpdateDeg*cDegToRad / abs(OmegaComet) * No2Si_V(UnitT_)

    call dir_to_xyz((90-LatSun)*cDegToRad, LonSun*cDegToRad, NormalSun_D)

    ! Maximum amount of data to be stored in ModBlockData
    ! In practice this is a rather generous overestimate.
    ! The first 8 variables: 5 are needed to store the neu1 face B.C. and
    ! 3 are used to store the face normal. The last 1 is for the shading
    ! in the photoionization.
    MaxBlockData = 8*(nI+1)*(nJ+1)*(nK+1) + nI*nJ*nK

    if (restart) then
       nStepSave_B          = n_step
       TimeSimulationSave_B = Time_Simulation
       nStepSaveCalcRates_B = n_step
       nStepPritSetFace     = n_step
       FaceCoordsTest_D     = (/FaceCoordsX, FaceCoordsY, FaceCoordsZ/)
    end if

    if(iProc==0)then
       write(*,*) 'rSphericalBodySi, rSphericalBody =', &
            rSphericalBodySi, rSphericalBody
       write(*,*) 'ProductionRateMaxSi, ProductionRateMax =', &
            ProductionRateMaxSi, ProductionRateMax
       write(*,*) 'ProductionRateMinSi, ProductionRateMin =', &
            ProductionRateMinSi, ProductionRateMin
       write(*,*) 'SolarAngleMaxDim, SolarAngleMax =', SolarAngleMaxDim, SolarAngleMax
       write(*,*) 'TempCometMinDim, TempCometMin =', TempCometMinDim, TempCometMin
       write(*,*) 'TempCometMaxDim, TempCometMax =', TempCometMaxDim, TempCometMax
       write(*,*) 'TempComet75Dim,  TempComet75  =', TempComet75Dim,  TempComet75
       write(*,*) 'TempToUnormal, TempToPressure =', TempToUnormal, TempToPressure
       write(*,*) 'MassFluid_I =', MassFluid_I
       write(*,*) 'TempToUn*sqrt(TempCometMax)  =', TempToUnormal*sqrt(TempCometMax)
       write(*,*) 'SlopeProduction, bProduction =', SlopeProduction, bProduction
       write(*,*) 'SlopeTemp, bTemp             =', SlopeTemp/Io2No_V(UnitTemperature_),&
            bTemp/Io2No_V(UnitTemperature_)
       write(*,*)'RotationComet, Omega  =', RotationCometHour, OmegaComet
       write(*,*)'AngleUpdateDeg, DtUpdateSi =', AngleUpdateDeg, DtUpdateSi
       write(*,*)'LatSun, LonSun=', LatSun, LonSun
       write(*,*)'NormalSun_D   =', NormalSun_D
    end if
  end subroutine user_init_session
  !===========================================================================
  subroutine user_set_boundary_cells(iBlock)

    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, &
         r_BLK, true_cell, x2
    use ModMain, ONLY: ProcTest, BlkTest, iTest, jTest, kTest

    integer, intent(in):: iBlock

    integer :: i, j, k
    real    :: XyzInside_D(3)
    logical :: DoTest, DoTestMe
    character(len=*), parameter :: NameSub='user_set_boundary_cells'
    !------------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    if (DoUseCGShape) then
       ! Place a point inside rMinShape sphere with transcendent coordinates
       ! to reduce chances of hitting the edge or corner of triangles

       XyzInside_D = rMinShape*(/cPi/10,cPi**2/50,cPi**3/700/)

       do k = MinK, MaxK; do j = MinJ, MaxJ; do i=MinI, MaxI
          ! Check if we are close enough
          if(r_BLK(i,j,k,iBlock) > rMaxShape) then
             IsBoundaryCell_GI(i,j,k,ExtraBc_) = .false.
          elseif(r_BLK(i,j,k,iBlock) < rMinShape) then
             IsBoundaryCell_GI(i,j,k,ExtraBc_) = .true.
          else
             ! Connect cell center with a point inside.
             ! If the line segment does not intersect the shape or it intersects
             ! even times then the point is inside the shape.
             IsBoundaryCell_GI(i,j,k,ExtraBc_) = .not. is_segment_intersected( &
                  XyzInside_D, Xyz_DGB(:,i,j,k,iBlock), IsOddIn=.true.)
          end if
       end do; end do; end do

       if(DoTestMe .and. iBlock==BLKtest .and. iProc==PROCtest)then
          write(*,*)NameSub,': iProc, iBlock =',iProc, iBlock
          write(*,*)NameSub,': is_segment_intersected =', &
               is_segment_intersected(XyzInside_D, &
               Xyz_DGB(:,Itest,Jtest,Ktest,iBlock), &
               IsOddIn=.true.)
          write(*,*)NameSub,': true cell?     ',&
               true_cell(Itest,Jtest,Ktest,iBlock)
       end if
    else
       IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
            r_BLK(:,:,:,iBlock) <= rSphericalBody
    end if

  end subroutine user_set_boundary_cells

  !=========================================================================
  subroutine read_shape_file

    use ModPhysics, ONLY: Si2No_V, UnitX_
    use ModIoUnit, ONLY: UnitTmp_
    use ModCoordTransform, ONLY: cross_product
    use ModRandomNumber, ONLY: random_real

    logical:: DoReadShapeFile = .true.

    integer:: nPoint, i, j, iPoint, iTriangle, iPoint1, iPoint2, iPoint3
    integer:: iSeed = 7

    real, allocatable:: Xyz_DI(:,:)

    character(len=100):: String1, String2

    character(len=*), parameter:: NameSub = 'read_shape_file'
    !-----------------------------------------------------------------------
    if(.not.DoReadShapeFile) RETURN
    DoReadShapeFile = .false.

    if(iProc==0)write(*,*) NameSub,' reading shape file ',trim(NameShapeFile)

    open(UnitTmp_, file=NameShapeFile)

    read(UnitTmp_, '(a)') String1
    read(UnitTmp_, *) String1, nPoint, String2
    if(String2 /= 'POINTS')call stop_mpi(NameSub//&
         ' no POINTS in '//trim(String2)//' of '//trim(NameShapeFile))
    read(UnitTmp_, *) String1, nTriangle, String2
    if(String2 /= 'TRIANGLES')call stop_mpi(NameSub//&
         ' no TRIANGLES in '//trim(String2)//' of '//trim(NameShapeFile))

    if(iProc==0)write(*,*) NameSub,' nPoint=', nPoint,' nTriangle=',nTriangle

    allocate(Xyz_DI(3,nPoint), &
         XyzTriangle_DII(3,3,nTriangle), Normal_DI(3,nTriangle))

    read(UnitTmp_, '(a)') String1
    do iPoint = 1, nPoint
       read(UnitTmp_,*) String1, i, j, Xyz_DI(:,iPoint)

       ! Perturb vertices of all triangles to avoid the the situation that 
       ! a line segment is parallel to a triangle plane in 
       ! is_segment_intersected

       Xyz_DI(1, iPoint) = Xyz_DI(1, iPoint) + random_real(iSeed)*1e-5
       Xyz_DI(2, iPoint) = Xyz_DI(2, iPoint) + random_real(iSeed)*1e-5
       Xyz_DI(3, iPoint) = Xyz_DI(3, iPoint) + random_real(iSeed)*1e-5

       ! Convert from SI units to normalized unit
       Xyz_DI(:,iPoint) = Xyz_DI(:,iPoint) * Si2No_V(UnitX_) 
    end do
    do iTriangle = 1, nTriangle
       read(UnitTmp_,*) String1, i, j, iPoint1, iPoint2, iPoint3
       XyzTriangle_DII(:,1,iTriangle) = Xyz_DI(:,iPoint1)
       XyzTriangle_DII(:,2,iTriangle) = Xyz_DI(:,iPoint2)
       XyzTriangle_DII(:,3,iTriangle) = Xyz_DI(:,iPoint3)

       Normal_DI(:,iTriangle) = cross_product( &
            Xyz_DI(:,iPoint2) - Xyz_DI(:,iPoint1), &
            Xyz_DI(:,iPoint3) - Xyz_DI(:,iPoint1))
       Normal_DI(:,iTriangle) = Normal_DI(:,iTriangle) / &
            sqrt(sum(Normal_DI(:,iTriangle)**2))
    end do

    !write(*,*)'!!! XyzTriangle_DII(:,:,1)=',XyzTriangle_DII(:,:,1)
    !write(*,*)'!!! XyzTriangle_DII(:,:,n)=',XyzTriangle_DII(:,:,nTriangle)

    rMinShape = sqrt(minval(sum(Xyz_DI**2,DIM=1)))
    rMaxShape = sqrt(maxval(sum(Xyz_DI**2,DIM=1)))

    if(iProc==0)write(*,*) NameSub,' rMinShape, rMaxShape=', &
         rMinShape, rMaxShape

    deallocate(Xyz_DI)

    close(UnitTmp_)

  end subroutine read_shape_file

  !==============================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModMain, ONLY: n_step, time_simulation
    use ModVarIndexes,   ONLY: nVar, Rho_, p_
    use ModFaceBoundary, ONLY: iFace, jFace, kFace, FaceCoords_D, &
         iBoundary, VarsTrueFace_V, iSide, iBlock => iBlockBc, TimeBc
    use ModGeometry,    ONLY: ExtraBc_, Xyz_DGB
    use ModPhysics, ONLY: UnitRho_, BodyRho_I, BodyP_I, Si2No_V, &
         ElectronPressureRatio
    use ModSolarwind,    ONLY: get_solar_wind_point
    use ModBlockData, ONLY: use_block_data, clean_block_data, &
         get_block_data, put_block_data

    logical :: DoTestHere=.true., IsIlluminated = .false.

    real, intent(out):: VarsGhostFace_V(nVar)

    integer:: iTrue, jTrue, kTrue, iBody, jBody, kBody


    real :: XyzIntersect_D(3), XyzStart_D(3), XyzEnd_D(3)
    real :: XyzTrueCell_D(3), XyzBodyCell_D(3)
    real :: Normal_D(3), CosAngle
    real :: TempCometLocal, uNormal, ProductionRateLocal

    logical :: DoWriteOnce = .true.

    integer :: iIonFluid
    real    :: uNormalIon_I(nIonFluid)

    character(len=*), parameter:: NameSub = 'user_set_face_boundary'
    !------------------------------------------------------------------------

    if (iProc /= 0) DoTestHere = .false.

!!! Outer boundaries
    if(iBoundary >0) then
       call get_solar_wind_point(TimeBc, FaceCoords_D, VarsGhostFace_V)

       ! Float condition for neu1
       VarsGhostFace_V(Neu1Rho_:Neu1P_)    = VarsTrueFace_V(Neu1Rho_:Neu1P_)
       RETURN
    end if

    if (iBoundary /= ExtraBc_) call stop_mpi(NameSub//' bad iBoundary value')

!!! Body boundaries

    ! Default body boundary conditions
    if (UseSwBC) then
       ! Boundary condition with solar wind values
       if (DoWriteOnce) then
          write(*,*) NameSub, ': solar wind body conditions.'
          DoWriteOnce = .false.
       end if

       call get_solar_wind_point(TimeBc, FaceCoords_D, VarsGhostFace_V)
    else
       ! Floating boundary condition
       if (DoWriteOnce .and. .not. UseReflectedBC) then
          write(*,*) NameSub, ': floating body conditions.'
          DoWriteOnce = .false.
       end if

       VarsGhostFace_V = VarsTrueFace_V
    end if

    !! Neutral boundary conditions -----------------------------------------
    if (DoUseCGShape) then

       ! We can use the saved values if no AMR is done
       if(use_block_data(iBlock)) then

          call get_block_data(iBlock, 5, VarsGhostFace_V(Neu1Rho_:Neu1P_))
          call get_block_data(iBlock, 3, Normal_D)

          if (.not. UseSwBC) call set_ion_face_boundary

          if ((n_step <= nStepPritSetFace+2) .and. &
               sum(abs(FaceCoords_D - FaceCoordsTest_D)) < 1e-8 ) then
             write(*,*) '=============== n_step ', n_step, '===================='
             write(*,*) 'FaceCoords_D  =', FaceCoords_D
             write(*,*) 'Normal_D      =', Normal_D
             write(*,*) 'P             =', VarsGhostFace_V(P_)
             write(*,*) 'Pe            =', VarsGhostFace_V(Pe_)
             write(*,*) 'SwRho         =', VarsGhostFace_V(SwRho_)
             write(*,*) 'SwU_D         =', VarsGhostFace_V(SwRhoUx_:SwRhoUz_)
             write(*,*) 'SwP           =', VarsGhostFace_V(SwP_)
             write(*,*) 'H2OpRho       =', VarsGhostFace_V(H2OpRho_)
             write(*,*) 'H2OpU_D       =', VarsGhostFace_V(H2OpRhoUx_:H2OpRhoUz_)
             write(*,*) 'H2OpP         =', VarsGhostFace_V(H2OpP_)
             write(*,*) 'Neu1Rho       =', VarsGhostFace_V(Neu1Rho_)
             write(*,*) 'Neu1u_D       =', VarsGhostFace_V(Neu1Ux_:Neu1Uz_)
             write(*,*) 'Neu1p         =', VarsGhostFace_V(Neu1p_)
          end if
          RETURN
       end if

       ! Empty the block storage if we redo the calculation
       if(use_block_data(iBlock)) call clean_block_data(iBlock)

       ! Save step and simulation time info
       nStepSave_B(iBlock)          = n_step
       TimeSimulationSave_B(iBlock) = Time_Simulation

       ! Default indexes for the true and body cells
       iTrue = iFace; jTrue = jFace; kTrue = kFace
       iBody = iFace; jBody = jFace; kBody = kFace

       select case(iSide)
       case(1)
          iBody = iFace - 1
       case(2)
          iTrue = iFace - 1
       case(3)
          jBody = jFace - 1
       case(4)
          jTrue = jFace - 1
       case(5)
          kBody = kFace -1
       case(6)
          kTrue = kFace -1
       end select

       XyzBodyCell_D = Xyz_DGB(:,iBody,jBody,kBody,iBlock)
       XyzTrueCell_D = Xyz_DGB(:,iTrue,jTrue,kTrue,iBlock)

       ! Find the intersection point between the true cell and the body cell
       ! that is closest to the true cell
       if (.not. is_segment_intersected(XyzTrueCell_D, XyzBodyCell_D, IsOddIn = .true., &
            XyzIntersectOut_D=XyzIntersect_D, NormalOut_D = Normal_D))then
          write(*,*) 'XyzTrueCell_D =', XyzTrueCell_D
          write(*,*) 'XyzBodyCell_D =', XyzBodyCell_D
          write(*,*) NameSub,' error for face =', iFace, jFace, kFace
          write(*,*) NameSub,' error for iside, iBlock=', iSide, iBlock
          call stop_mpi(NameSub// &
               ': No intersection points are found between true and the body cells')
       end if

       ! Fix the normal direction if it is not pointing outward
       if (sum(Normal_D*(XyzTrueCell_D - XyzBodyCell_D)) < 0.0) &
            Normal_D = -Normal_D
    else
       ! Use a spherical body instead of a real CG shape
       ! Normal is in the r direction
       Normal_D = FaceCoords_D/sqrt(sum(FaceCoords_D*FaceCoords_D))
    end if

    ! Calculate the cos angle between the surface normal and the sun direction
    CosAngle = sum(Normal_D*NormalSun_D)

    ! Set local outflow parameters as default that may be overwritten if illuminated
    TempCometLocal      = TempCometMin
    ProductionRateLocal = ProductionRateMin

    if (CosAngle > 0.0) then
       if (DoUseCGShape) then
          ! See whether the intersection point is in the shade by going towards the Sun
          ! and checking for intersection with the shape
          XyzStart_D = XyzIntersect_D + 1e-9*rMaxShape*NormalSun_D
          XyzEnd_D   = XyzIntersect_D +    2*rMaxShape*NormalSun_D
          if(.not.is_segment_intersected(XyzStart_D, XyzEnd_D)) &
               IsIlluminated = .true.
       else
          IsIlluminated = .true.
       end if
    end if

    if (IsIlluminated) then
       ! Increase temperature of the face if it is illuminated
       TempCometLocal      = max( TempCometMin, &
            SlopeTemp / CosAngle + bTemp)
       ! Increase neutral production rate
       ProductionRateLocal = max( ProductionRateMin, &
            SlopeProduction * CosAngle + bProduction )
    end if

    ! Calculate the normal velocity
    uNormal = sqrt(TempCometLocal)*TempToUnormal

    VarsGhostFace_V(Neu1Ux_:Neu1Uz_) = Normal_D*uNormal
    VarsGhostFace_V(Neu1Rho_)    = ProductionRateLocal/uNormal*NumdensToRho
    VarsGhostFace_V(Neu1P_)      = &
         VarsGhostFace_V(Neu1Rho_)*TempCometLocal*TempToPressure

    if (.not. UseSwBC) call set_ion_face_boundary

    if (DoTestHere .and. IsIlluminated .and. CosAngle > 0.5) then
       FaceCoordsTest_D = FaceCoords_D

       write(*,*) 'FaceCoords_D: ', FaceCoords_D
       if (DoUseCGShape) then
          write(*,*) 'TestFace_D: ', (XyzBodyCell_D + XyzTrueCell_D)/2
          write(*,*) 'XyzTrueCell_D =', XyzTrueCell_D
          write(*,*) 'XyzBodyCell_D =', XyzBodyCell_D
          write(*,*) 'XyzIntersect_D=', XyzIntersect_D
          write(*,*) 'XyzStart_D    =', XyzStart_D
          write(*,*) 'XyzEnd_D      =', XyzEnd_D
       end if
       write(*,*) 'Normal_D      =', Normal_D
       write(*,*) 'CosAngle      =', CosAngle
       write(*,*) 'P             =', VarsGhostFace_V(P_)
       write(*,*) 'Pe            =', VarsGhostFace_V(Pe_)
       write(*,*) 'SwRho         =', VarsGhostFace_V(SwRho_)
       write(*,*) 'SwU_D         =', VarsGhostFace_V(SwRhoUx_:SwRhoUz_)
       write(*,*) 'SwP           =', VarsGhostFace_V(SwP_)
       write(*,*) 'H2OpRho       =', VarsGhostFace_V(H2OpRho_)
       write(*,*) 'H2OpU_D       =', VarsGhostFace_V(H2OpRhoUx_:H2OpRhoUz_)
       write(*,*) 'H2OpP         =', VarsGhostFace_V(H2OpP_)
       write(*,*) 'Neu1Rho       =', VarsGhostFace_V(Neu1Rho_)
       write(*,*) 'Neu1u_D       =', VarsGhostFace_V(Neu1Ux_:Neu1Uz_)
       write(*,*) 'Neu1uNormal   =', uNormal
       write(*,*) 'Neu1p         =', VarsGhostFace_V(Neu1p_)
       DoTestHere=.false.
       nStepPritSetFace = n_step
    end if

    IsIlluminated = .false.

    ! Store for future time steps
    if (DoUseCGShape)  then
       call put_block_data(iBlock, 5, VarsGhostFace_V(Neu1Rho_:Neu1P_))
       call put_block_data(iBlock, 3, Normal_D)
    end if

  contains
    !=====================================================================
    subroutine set_ion_face_boundary

      integer :: iUx, iUz

      ! Projection length of U_ on the local surface radius vector
      do iIonFluid=1,nIonFluid
         uNormalIon_I(iIonFluid) = sum(&
              VarsTrueFace_V(iRhoUxIon_I(iIonFluid):iRhoUzIon_I(iIonFluid))* &
              Normal_D)
      end do

      !  BdotR = dot_product(VarsTrueFace_V(Bx_:Bz_),FaceCoords_D)/ &
      !  dot_product(FaceCoords_D,FaceCoords_D)

      ! Projection vectors
      !  BRefl_D = BdotR*FaceCoords_D


      ! Bz component propagated through moon, Bx and By didn't
      !  VarsGhostFace_V(Bx_:By_) = 0.0
      !  VarsGhostFace_V(Bz_)     = SW_Bz

      do iIonFluid=1,nIonFluid
         if (UseReflectedBC) then
            ! Reflected raidal velocity uG = uT - 2*u_normal
            if (DoWriteOnce) then
               write(*,*) NameSub, ': reflected boundary condition'
               DoWriteOnce = .false.
            end if

            if (uNormalIon_I(iIonFluid) > 0.0) then
               iUx = iRhoUxIon_I(iIonFluid)
               iUz = iRhoUzIon_I(iIonFluid)
               VarsGhostFace_V(iUx:iUz) = VarsTrueFace_V(iUx:iUz) - &
                    2.0*uNormalIon_I(iIonFluid)*Normal_D
            end if
         else
            ! set outward flux body value (Comet's surface not considered as plasma source)
            ! leave inward flux untouched
            if (uNormalIon_I(iIonFluid) > 0.0) then
               iUx = iRhoUxIon_I(iIonFluid)
               iUz = iRhoUzIon_I(iIonFluid)
               !VarsGhostFace_V(iUx_I(iIonFluid):iUz_I(iIonFluid)) = 0.0
               VarsGhostFace_V(iUx:iUz) = 0.0
               VarsGhostFace_V(iRhoIon_I(iIonFluid)) = BodyRho_I(iIonFluid)
               VarsGhostFace_V(iPIon_I(iIonFluid))   = BodyP_I(iIonFluid)
            endif
         end if
      end do

      do iIonFluid=1,nIonFluid
         if (uNormalIon_I(iIonFluid) > 0.0 .and. &
              any(VarsGhostFace_V(iRhoIon_I) > 1e-2*Si2NO_V(UnitRho_))) then
            write(*,*) 'n_step, iIonFluid, FaceCoords_D =', &
                 n_step, iIonFluid, FaceCoords_D
            write(*,*) 'uNormalIon_I(iIonFluid)                =', &
                 uNormalIon_I(iIonFluid)
            write(*,*) 'BodyRho_I(iIonFluid), BodyP_I(iIonFluid) =', &
                 BodyRho_I(iIonFluid), BodyP_I(iIonFluid)
            call stop_mpi('Plasma source at the surface????????')
         end if
      end do

      VarsGhostFace_V(Rho_)   = sum(VarsGhostFace_V(iRhoIon_I))
      VarsGhostFace_V(RhoUx_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUxIon_I))/ &
           sum(VarsGhostFace_V(iRhoIon_I))
      VarsGhostFace_V(RhoUy_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUyIon_I))/ &
           sum(VarsGhostFace_V(iRhoIon_I))
      VarsGhostFace_V(RhoUz_) = sum(VarsGhostFace_V(iRhoIon_I)*VarsGhostFace_V(iRhoUzIon_I))/ &
           sum(VarsGhostFace_V(iRhoIon_I))

      if(UseElectronPressure) then
         VarsGhostFace_V(P_)  = sum(VarsGhostFace_V(iPIon_I))
         VarsGhostFace_V(Pe_) = VarsGhostFace_V(P_)*ElectronPressureRatio
      else
         VarsGhostFace_V(P_)  = sum(VarsGhostFace_V(iPIon_I)) &
              *(1.+ElectronPressureRatio)
      end if
    end subroutine set_ion_face_boundary
    !===========================================================
  end subroutine user_set_face_boundary

  !============================================================================
  logical function is_segment_intersected(Xyz1_D, Xyz2_D, &
       IsOddIn, XyzIntersectOut_D, NormalOut_D)

    ! Check if a line segment connecting Xyz1_D and Xyz2_D intersects
    ! the shape. If IsEvenIn is present and true, check if the number
    ! of intersections is odd or even.

    ! Using algorithm: http://geomalgorithms.com/a06-_intersect-2.html

    real, intent(in):: Xyz1_D(3) ! segment start coordinates
    real, intent(in):: Xyz2_D(3) ! segment end   coordinates
    logical, optional:: IsOddIn  ! check for odd number of intersects
    real,    optional :: XyzIntersectOut_D(3)
    real,    optional :: NormalOut_D(3)

    logical:: IsOdd

    integer:: iTriangle, nIntersect, iMinRatio

    integer, parameter:: MaxIntersect = 10
    real:: Ratio, Ratio1, Ratio2, Ratio_I(MaxIntersect)
    integer:: iTriangle_I(MaxIntersect)
    real, dimension(3):: v1_D, Xyz_D, u_D, v_D, w_D
    real:: nDotP2P1, nDotV1P1, u2, v2, uDotV, wDotU, wDotV, InvDenom

    character(len=*), parameter:: NameSub = 'is_segment_intersected'
    !------------------------------------------------------------------------
    ! Default is to check for first intersection only
    ! (for shadows and convex shapes)
    IsOdd = .false.
    if(present(IsOddIn))   IsOdd = IsOddIn

    nIntersect = 0
    do iTriangle = 1, nTriangle
       ! Find intersection of line segment with the plane of the triangle
       nDotP2P1 = sum(Normal_DI(:,iTriangle)*(Xyz2_D - Xyz1_D))

       ! Check if the segment is parallel to the plane
       if(abs(nDotP2P1) < 1e-12) then
          if (abs(sum(Normal_DI(:,iTriangle)*(Xyz1_D - XyzTriangle_DII(:,1,iTriangle)))) &
               < 1e-12) then
             write(*,*) 'segment lies in the same plane: iTriangle: ', Xyz2_D
             write(*,*) 'Test:', &
                  sum(Normal_DI(:,iTriangle)*(Xyz1_D - XyzTriangle_DII(:,1,iTriangle)))
             CYCLE
          else 
             CYCLE
          end if
       end if

       ! Vertex 1 of triangle
       v1_D = XyzTriangle_DII(:,1,iTriangle)

       nDotV1P1 = sum(Normal_DI(:,iTriangle)*(v1_D -  Xyz1_D))

       ! Intersection is at P1 + Ratio * (P2 - P1)
       Ratio = nDotV1P1 / nDotP2P1

       ! Check if point is inside the segment
       if(Ratio <    -1e-12) CYCLE
       if(Ratio > 1.0+1e-12) CYCLE

       ! Intersection point
       Xyz_D =  Xyz1_D + Ratio*(Xyz2_D -  Xyz1_D)

       ! Calculate the barycentric coordinates Ratio1 and Ratio2
       ! The intersection point is inside if the conditions
       ! 0 < Ratio1, Ratio2 and Ratio1 + Ratio2 < 1 both hold.

       ! Vectors relative to the first vertex
       u_D = XyzTriangle_DII(:,2,iTriangle) - v1_D 
       v_D = XyzTriangle_DII(:,3,iTriangle) - v1_D 
       w_D = Xyz_D                          - v1_D

       u2 = sum(u_D**2)
       v2 = sum(v_D**2)
       uDotV = sum(u_D*v_D)
       wDotU = sum(w_D*u_D)
       wDotV = sum(w_D*v_D)

       InvDenom = 1.0/(uDotV**2 - u2*v2)

       ! First barycentric coordinate
       Ratio1 = (uDotV*wDotV - v2*wDotU)*InvDenom
       if(Ratio1 <    -1e-12) CYCLE
       if(Ratio1 > 1.0+1e-12) CYCLE

       ! Second barycentric coordinate
       Ratio2 = (uDotV*wDotU - u2*wDotV)*InvDenom

       !       if (abs(Xyz2_D(1) -xTest) <= 1e-4 .and. &
       !            abs(Xyz2_D(2) - yTest) <= 1e-4 .and. &
       !            abs(Xyz2_D(3) - zTest) <= 1e-4) then
       !          if (abs(Ratio2) < 1e-10 .or. abs(Ratio2 - 1.0) < 1e-10) then
       !             write(*,*) 'iTriangle: ', iTriangle, 'Ratio1: ', Ratio1, 'Ratio2: ', Ratio2
       !          end if
       !       end if

       if(         Ratio2 <    -1e-12) CYCLE
       if(Ratio1 + Ratio2 > 1.0+1e-12) CYCLE

       if(.not. IsOdd)then
          ! The line segment intersected the triangle
          is_segment_intersected = .true.
          RETURN
       end if

       ! Check if this intersection is different from previous ones
       if(nIntersect > 0)then
          if( any(abs(Ratio - Ratio_I(1:nIntersect)) < 1e-12) ) CYCLE
       end if

       ! New intersection was found
       nIntersect =  nIntersect + 1

       if(nIntersect > MaxIntersect) call stop_mpi(NameSub// &
            ': too many intersections, increase MaxIntersect')

       ! Store the position along the segment into Ratio_I
       Ratio_I(nIntersect) = Ratio
       iTriangle_I(nIntersect) = iTriangle
    end do

    ! Only record the closeset intersection to the true cell
    if(present(NormalOut_D)) then
       iMinRatio = minloc(Ratio_I(1:nIntersect), 1)
       iTriangle = iTriangle_I(iMinRatio)
       NormalOut_D = Normal_DI(:,iTriangle)

       if (nIntersect > 1) then          
          !          write(*,*) 'nIntersect: ', nIntersect, 'Ratio_I: ', Ratio_I(1:nIntersect)
          !          write(*,*) 'Xyz1_D: ', Xyz1_D
          !          write(*,*) 'Xyz2_D: ', Xyz2_D
          !          do iIntersect = 1, nIntersect
          !             iTriangle = Triangle_I(iIntersect)
          !             write(*,*) XyzTriangle_DII(:,1,iTriangle)
          !             write(*,*) XyzTriangle_DII(:,2,iTriangle)
          !             write(*,*) XyzTriangle_DII(:,3,iTriangle)
          !          end do
          write(*,*)  'Ratio_I(iMinRatio)            = ', Ratio_I(iMinRatio) 
          write(*,*)  'minval(Ratio_I(1:nIntersect)) = ', minval(Ratio_I(1:nIntersect))
          write(*,*)  'Ratio_I(1:nIntersect)', Ratio_I(1:nIntersect)
       end if
    end if
    if(present(XyzIntersectOut_D)) then
       XyzIntersectOut_D = Xyz1_D + minval(Ratio_I(1:nIntersect))*(Xyz2_D -  Xyz1_D)
    end if

    if(.not. IsOdd)then
       ! The line segment was not intersected by any triangle
       is_segment_intersected = .false.
       RETURN
    end if

    ! We got to the other side of the surface
    ! if there were odd number of intersections
    is_segment_intersected = modulo(nIntersect, 2) == 1

  end function is_segment_intersected


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++ From Martin & Xianzhe ++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_electron_collision_rates(Te,nElec,i,j,k,iBlock,fen_I,fei_I)

    ! calculate all collision rates for electrons (fen, fei)
    ! (used for sources & resistivity)

    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: No2SI_V, UnitN_
    use ModMultiFluid, ONLY: MassIon_I, ChargeIon_I
    use ModGeometry,   ONLY: Xyz_DGB

    integer,intent(in) :: i,j,k,iBlock   
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(out)   :: fen_I(nNeutral)
    real,intent(out)   :: fei_I(nIonFluid)

    real :: sqrtTe

    logical :: DoTestElectronCollision = .false.
    character(len=*), parameter :: NameSub='calc_electron_collision_rates'
    !----------------------------------------------------------------------

    !! electron-neutral and electron-ion collision rates
    !! provide all rates in SI units

    if (Te < 0) then
       write(*,*) NameSub,': i,j,k,iBlock =', i,j,k,iBlock
       write(*,*) NameSub,': State_VGB    =', State_VGB(:,i,j,k,iBlock)
       write(*,*) NameSub,': Xyz_DGB      =', Xyz_DGB(:,i,j,k,iBlock)
       call stop_mpi(NameSub//': negative Te')
    end if

    !! reduced temperature ~ Te
    sqrtTe = sqrt(Te)

    !! initialize all collision rates with zero
    fei_I = 0. ; fen_I = 0.

    !! Electron - neutral collision rates
    !! e - H2O,  Itikawa, Planet. Space Sci., 1971 and Itikawa, Phys. Fluids 1983
    !    fen_I(H2O_) = 2.745E-5*NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)/1e6*Te**(-0.62)      !! rate in [1/s]
    fen_I(Neu1_) = 2.745E-5*State_VGB(Neu1Rho_, i, j, k, iBlock)/1e6*Te**(-0.62)/MassFluid_I(nFluid)*No2SI_V(UnitN_)

    !!Electron - ion collision rates
    !! e - H2Op, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fei_I(H2Op_) = 54.5*ChargeIon_I(H2Op_)**2*State_VGB(H2OpRho_,i,j,k,iBlock)/MassIon_I(H2Op_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)      !! rate in [1/s]
    !! e - Hp, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fei_I(SW_) = 54.5*ChargeIon_I(SW_)**2*State_VGB(SwRho_,i,j,k,iBlock)/MassIon_I(SW_)* &
         No2SI_V(UnitN_)/1E6/(Te*sqrtTe)      !! rate in [1/s]

    if (DoTestElectronCollision) then
       write(*,*) 'fen_I  =', fen_I
       write(*,*) 'fei_I  =', fei_I
    end if

  end subroutine calc_electron_collision_rates

  !========================================================================
  subroutine user_calc_rates(Ti_I,Te,i,j,k,iBlock,nElec,nIon_I,fin_II,fii_II,fie_I,alpha_I,kin_IIII,v_II,&
       ve_II,uElec_D,uIon_DI,Qexc_II,Qion_II, DoCalcShading, IsIntersectedShapeR)

    ! calculate all rates not involving electron collisions

    use ModConst,    ONLY: cElectronMass, cProtonMass
    use ModMain,     ONLY: iTest, jTest, kTest, BlkTest, &
         n_step, iTest, ProcTest
    use ModGeometry, ONLY: R_BLK, Xyz_DGB

    integer,intent(in) :: i,j,k,iBlock
    real,intent(in)    :: Ti_I(nIonFluid)
    real,intent(in)    :: Te
    real,intent(in)    :: nElec
    real,intent(in)    :: nIon_I(nIonFluid)
    real,intent(in)    :: uElec_D(3)
    real,intent(in)    :: uIon_DI(3,nIonFluid)
    logical,intent(in) :: DoCalcShading
    real,intent(out)   :: fin_II(nIonFluid,nNeutral)
    real,intent(out)   :: fii_II(nIonFluid,nIonFluid)
    real,intent(out)   :: fie_I(nIonFluid)
    real,intent(out)   :: alpha_I(nIonFluid)
    real,intent(out)   :: kin_IIII(nIonFluid,nNeutral,nNeutral,nIonFluid)
    real,intent(out)   :: v_II(nNeutral,nIonFluid)
    real,intent(out)   :: ve_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qexc_II(nNeutral,nIonFluid)
    real,intent(out)   :: Qion_II(nNeutral,nIonFluid)
    real,intent(inout) :: IsIntersectedShapeR

    real :: Tred, Mred
    real :: DistProjection2, CosAngleTmp, NCol, sigma, J3, log10Te, sqrtTe
    real,dimension(nNeutral,nIonFluid) :: sigma_e
    integer :: n
    real, save :: ElImpRate_I(nNeutral,61)!, ElCrossSect_I(61)


    real, save :: sigmaeh2o = 4.53E-21  !! Ionization cross section for 20 eV electrons [m^2]
    real, save :: ve = 2.65E6           !! Speed of 20 eV electrons [m/s]

    logical :: DoTest, DoTestMe
    logical :: IsIntersectedShape

    real    :: nTmp

    character(len=*), parameter:: NameSub = 'user_calc_rates'
    !-----------------------------------------------------------------
    if(iBlock==BlkTest .and. i==iTest .and. j==jTest .and. &
         k==kTest .and. iProc==ProcTest) then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.
       DoTestMe=.false.
    end if

    ! H2O and H electron impact rate depending on electron temperature (Cravens et al. 1987)
    ElImpRate_I(Neu1_,1:61) = (/ 0.00E+00, 1.14E-16, 2.03E-16, 3.04E-16, 4.37E-16, 6.34E-16, 9.07E-16, &
         1.28E-15, 1.79E-15, 2.34E-15, 3.15E-15, 4.35E-15, 5.54E-15, 6.90E-15, 8.47E-15, 1.05E-14, 1.25E-14, &
         1.51E-14, 1.80E-14, 2.09E-14, 2.41E-14, 2.74E-14, 3.09E-14, 3.46E-14, 3.90E-14, 4.34E-14, 4.80E-14, &
         5.24E-14, 5.70E-14, 6.15E-14, 6.63E-14, 7.15E-14, 7.61E-14, 8.03E-14, 8.47E-14, 8.90E-14, 9.30E-14, &
         9.72E-14, 1.01E-13, 1.05E-13, 1.08E-13, 1.11E-13, 1.14E-13, 1.16E-13, 1.18E-13, 1.21E-13, 1.23E-13, &
         1.24E-13, 1.25E-13, 1.25E-13, 1.25E-13, 1.23E-13, 1.23E-13, 1.23E-13, 1.21E-13, 1.20E-13, 1.17E-13, &
         1.15E-13, 1.12E-13, 1.10E-13, 1.07E-13 /)

    ! ! Hydrogen-electron impact ionization cross section depending on electron energy
    ! ElCrossSect_I(1:61) = (/ 0.000E-00, 0.114E-20, 0.189E-20, 0.257E-20, 0.320E-20, 0.377E-20, &
    !      0.429E-20, 0.473E-20, 0.510E-20, 0.542E-20, 0.568E-20, 0.587E-20, 0.600E-20, 0.609E-20, 0.613E-20, &
    !      0.613E-20, 0.608E-20, 0.600E-20, 0.588E-20, 0.575E-20, 0.559E-20, 0.541E-20, 0.521E-20, 0.501E-20, &
    !      0.480E-20, 0.457E-20, 0.435E-20, 0.414E-20, 0.392E-20, 0.369E-20, 0.348E-20, 0.328E-20, 0.307E-20, &
    !      0.288E-20, 0.270E-20, 0.252E-20, 0.235E-20, 0.219E-20, 0.204E-20, 0.190E-20, 0.176E-20, 0.163E-20, &
    !      0.152E-20, 0.141E-20, 0.130E-20, 0.120E-20, 0.111E-20, 0.103E-20, 0.095E-20, 0.088E-20, 0.081E-20, &
    !      0.075E-20, 0.069E-20, 0.063E-20, 0.059E-20, 0.054E-20, 0.050E-20, 0.045E-20, 0.042E-20, 0.039E-20, &
    !      0.036E-20 /)
    !----------------------------------------------------------------------

    !! provide all rates in SI units

    !! initialize all collision rates with zero
    fin_II = 0. ; fii_II = 0. ; fie_I = 0. ; kin_IIII = 0. ; alpha_I = 0. ; v_II = 0. ; ve_II = 0.
    sigma_e = 0. ; Qexc_II = 0. ; Qion_II = 0.

    !! Ionization rates (add opacity correction/shadowing when needed)
    v_II(Neu1_,H2Op_) = vHI                           !! from PARAM.in
    v_II(Neu1_,SW_) = vHI*1e-10                       !! to avoid too low densities

    !! Electron excess energies from ionization (increases electron pressure)
    Qexc_II(Neu1_,H2Op_) = 1.9226E-18 ! 12.0 eV, Huebner 1992

    !! Ionization potential for electron impact ionization (needs to be delivered by the ionizing electron)
    Qion_II(Neu1_,H2Op_) = 2.02e-18 ! 12.6 eV Joshipura et al. (2007)

    ! ! UV opacity
    J3 = 4.5E14/(rHelio**2) ! J3 = 4.5E14 [m^-2*s^-1]: lambda < 984A solar flux @ 1 AU, Marconi, 1982)
    sigma = (v_II(Neu1_,H2Op_))/J3
    ! Alternative:
    ! Cross section of J3 flux for ionization (lambda < 984A) [m^2], Marconi, 1982)
    ! sigma_13=2.4E-18 cm^2: H2O + hv -> OH + H
    ! sigma_23=1.0E-18 cm^2: H2O + hv -> H2 + O(1D)
    ! sigma_33=8.4E-18 cm^2: H2O + hv -> H2Op + e
    ! sigma = 1.18E-21 ! sum(sigma_i3)

    ! (+++++++++need to fixed++++++++++++)
    ! Dist distance from sun-comet line, only neutral H2O considered
    !    Dist = sqrt(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)*&
    !         rPlanetSI+0.1
    !    write(*,*)  'rBody =', rBody

    ! New Block, need to check whether the cell is in the shade
    if(DoCalcShading .and. DoUseCGShape) then

       if (i == 1 .and. j == 1 .and. k ==1 .and. iBlock ==1) then
          write(*,*) NameSub, ': doing calculations. n_step, iProc =', &
               n_step, iProc
       end if

       CosAngleTmp     = sum(Xyz_DGB(:,i,j,k,iBlock)*NormalSun_D)
       DistProjection2 = R_BLK(i,j,k,iBlock)**2 - CosAngleTmp**2

       if (DistProjection2 < rMinShape**2 .and. CosAngleTmp < 0) then
          IsIntersectedShapeR = 1.0
       else if (DistProjection2 < rMinShape**2 .and. CosAngleTmp > 0) then
          IsIntersectedShapeR = 0.0
       else if (DistProjection2 > rMaxShape**2) then
          IsIntersectedShapeR = 0.0
       else
          if (is_segment_intersected(Xyz_DGB(:,i,j,k,iBlock), &
               Xyz_DGB(:,i,j,k,iBlock)+5*rMaxShape*NormalSun_D)) then
             IsIntersectedShapeR = 1.0
          else
             IsIntersectedShapeR = 0.0
          end if
       end if
    end if

    if (.not.DoUseCGShape) then
       CosAngleTmp    = sum(Xyz_DGB(:,i,j,k,iBlock)*NormalSun_D)
       DistProjection2= R_BLK(i,j,k,iBlock)**2 - CosAngleTmp**2
       if (DistProjection2 < rSphericalBody**2 .and. CosAngleTmp < 0) then
          IsIntersectedShapeR = 1.0
       else
          IsIntersectedShapeR = 0.0
       end if
    end if

    if (IsIntersectedShapeR == 1.0) then
       IsIntersectedShape = .true.
    else if (IsIntersectedShapeR == 0.0) then
       IsIntersectedShape = .false.
    else
       write(*,*) 'iProc, iBlock =', iProc, iBlock
       write(*,*) 'IsIntersectedShapeR =',&
            IsIntersectedShapeR
       call stop_mpi('IsIntersectedShapeR /= 0.0 or 1.0')
    end if

    if (IsIntersectedShape) then
       v_II = v_II*1e-9
    else
       NCol = 0
       v_II = v_II*exp(-sigma*NCol) + v_II*1e-9
    end if

    !    if ( is_segment_intersected(Xyz_DGB(:,i,j,k,iBlock), Xyz_DGB(:,i,j,k,iBlock)+5*rMaxShape*NormalSun_D) ) then
    !       v_II = v_II*1e-9 ! Inside the body's shadow
    !    else
    !       ! N total number of water-type molecules in upstream column
    !       !       uNeutr = sqrt(UnxNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2+UnyNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2+&
    !       !            UnzNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)**2)
    !!       uNeutr = sqrt(sum((State_VGB(Neu1RhoUx_:Neu1RhoUz_,i,j,k,iBlock)/State_VGB(Neu1Rho_,i,j,k,iBlock))**2))
    !       !  (+++++++++need to fixed++++++++++++)
    !       !       NCol = Qprod/uNeutr/Dist/4.*(-atan(Xyz_DGB(x_,i,j,k,iBlock)*rPlanetSI/Dist)/cPi+0.5)
    !       NCol = 0
    !       v_II = v_II*exp(-sigma*NCol) + v_II*1e-9
    !       !if(i==iTest.and.j==jTest.and.k==kTest.and.iBlock==BlkTest) then
    !       !   write(*,*)'sigma      = ',sigma
    !       !   write(*,*)'Dist       = ',Dist
    !       !   write(*,*)'NCol       = ',NCol
    !       !   write(*,*)'Correction = ',exp(-sigma*NCol)
    !       !end if
    !    end if

    if (Te <= 0.0) then
       write(*,*) '!!!!!!!!!!!!!!!!!!!!!! Te <= 0: Te =', Te
       call stop_mpi('!!!!!!!!!!!!!!!!!!!!!! Te<0')
    end if

    log10Te = log10(max(Te,1.0))
    nTmp = (log10Te-4.45053516)/0.0415 + 1.0

    if (abs(nTmp) > 1e6) then
       write(*,*) 'iProc, Te, log10Te, nTmp   =', iProc, Te, log10Te, nTmp
    end if

    ! H2O electron impact ionization cross section after Cravens et al. 1987: H2O & e -> H2Op + 2e
    nTmp = min(60.0, max(1.0,  nTmp))
    n = int(nTmp)
    if (n > 60) n = 60
    if (n < 1) n = 1
    ve_II(Neu1_,H2Op_) = max(nElec*((log10Te-((n-1.0)*0.0415+4.45053516))/0.0415*(ElImpRate_I(Neu1_,n+1)-ElImpRate_I(Neu1_,n))+&
         ElImpRate_I(Neu1_,n)),0.0) ! linear interpolation
    v_II(Neu1_,H2Op_) = v_II(Neu1_,H2Op_) + ve_II(Neu1_,H2Op_) ! v_II is the total ionization rate, photons and electrons!
    ! TestArray(1,i,j,k,iBlock) = ve_II(H2O_,H2Op_)/v_II(H2O_,H2Op_)

    ! Number of energetic electrons (number of equivalent 20 eV electrons)
    ne20eV_GB(i,j,k,iBlock) = ve_II(Neu1_,H2Op_)/(sigmaeh2o*ve)

    !! ********** Ion-neutral collision/charge exchange rates ********** 
    !! Example(s)
    ! resonant H+ & O -> O+ & H  subtracts H+ and adds O+
    ! kin_IIII(Hp_,O_,Op_,H_) = 6.61E-11/1E6*sqrt(Ti_I(Hp_))*(1.0-0.047*log10(Ti_I(Hp)))**2    !! rate in [m^3/s]
    ! resonant O+ & H -> H+ & O  subtracts O+ and adds H+
    ! kin_IIII(Op_,H_,Hp_,O_) = 4.63E-12/1E6*sqrt(TnNeutral(H_,i-MinI+1,j-MinJ+1,k-MinK+1)+TOp_/16.)  !! rate in [m^3/s]


    !! H2Op & H2O -> H2Op & H2O    ! non-resonant
    !! fin_II(H2Op_,H2O_) = 0.
    !! H2Op & H2O -> H2O & H2Op    ! resonant
    kin_IIII(H2Op_,Neu1_,Neu1_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996)

    !! Hp & H2O -> H & H2Op    ! resonant, estimate to get the same drag on SW-protons, neutral H2O product is just placeholder for fast neutral hydrogen (unused here)
    !uSWBulk2 = sum(State_VGB(iRhoUxIon_I(SW_):iRhoUzIon_I(SW_),i,j,k,iBlock)**2) / &
    !     State_VGB(iRhoIon_I(SW_),i,j,k,iBlock)**2*No2SI_V(UnitU_)**2
    !uSWTherm2 = 8.*cBoltzmann*Ti_I(SW_)/(MassIon_I(SW_)*cProtonMass*cPi)
    !kin_IIII(SW_,H2O_,H2O_,H2Op_) = 1E-4*3.0E-15*sqrt(uSWBulk2 + uSWTherm2) ! sigma=3e-15 cm^2, Cometopause of comet Halley, Ip, ApJ, 343, 956-952, 1989

    !! Hp & H2O -> H2O & H2Op    ! resonant, estimate to get the same drag on SW-protons, neutral H2O product is just placeholder for fast neutral hydrogen (unused here)
    kin_IIII(SW_,Neu1_,Neu1_,H2Op_) = 1E-6*1.7E-9 !! Gombosi et al., J. Geophys. Res., (1996), estimated

    !! ********** Ion-ion collision rates ********** 
    ! SWp - SWp is left zero because the they do not result in a change in the source terms
    ! Hp - Hp is left zero because the they do not result in a change in the source terms
    ! H2Op - H2Op is left zero because the they do not result in a change in the source terms

    ! H2Op - SWp, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    Tred = (MassIon_I(H2Op_)*Ti_I(SW_)+MassIon_I(SW_)*Ti_I(H2Op_))/(MassIon_I(H2Op_)+MassIon_I(SW_)) ! reduced temp
    Mred = MassIon_I(H2Op_)*MassIon_I(SW_)/(MassIon_I(H2Op_)+MassIon_I(SW_)) ! reduced mass
    fii_II(H2Op_,SW_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SW_)**2/MassIon_I(H2Op_)*&
         sqrt(Mred)*1e-6*nIon_I(SW_)/(Tred*sqrt(Tred))
    ! SWp - H2Op, Coulomb collision, Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    fii_II(SW_,H2Op_) = 1.27*ChargeIon_I(H2Op_)**2*ChargeIon_I(SW_)**2/MassIon_I(SW_)*&
         sqrt(Mred)*1e-6*nIon_I(H2Op_)/(Tred*sqrt(Tred))

    !! Ion - electron collision rates, reduced mass=~me and reduced temperature=~Te
    sqrtTe = sqrt(Te)
    !! H2Op - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(H2Op_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(H2Op_)*ChargeIon_I(H2Op_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]
    !! Hp - e, Schunk and Nagy, Ionospheres, Cambridge University Press, 2000
    fie_I(SW_) = 1.27*sqrt(cElectronMass/cProtonMass)/MassIon_I(SW_)*ChargeIon_I(SW_)**2*nElec/ &
         1E6/(Te*sqrtTe)      !! rate in [1/s]

    !! ********** Ion-electron recombination rates ********** 

    !! Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    if (Te < 800.) then 
       alpha_I(H2Op_) = 1E-6*1.57E-5*Te**(-0.569) !! rate in [m^3/s]
    elseif (Te<4000) then
       alpha_I(H2Op_) = 1E-6*4.73E-5*Te**(-0.74)  !! rate in [m^3/s]
    else
       alpha_I(H2Op_) = 1E-6*1.03E-3*Te**(-1.111) !! rate in [m^3/s]
    end if

    !     if (Te < 200.) then 
    !        alpha_I(H2Op_) = 1E-6*7E-7*sqrt(300./Te) !! rate in [m^3/s]
    !     else
    !        alpha_I(H2Op_) = 2.342*1E-6*7E-7*Te**(0.2553-0.1633*log10(Te)) !! rate in [m^3/s]
    !     end if

    alpha_I(SW_)   = 1E-6*4.8E-12*(250/Te)**0.7  !! Schunk and Nagy, Ionospheres,Cambridge University Press, 2000
    !alpha_I(SW_)   = 1E-6*3.5E-12*(Te/300)**(-0.7)  !! Schmidt et al., Comput. Phys. Commun. (1988)


    if (DoTestMe) then
       write(*,*) NameSub
       write(*,*) ' fin_II   =', fin_II
       write(*,*) ' fii_II   =', fii_II
       write(*,*) ' fie_I    =', fie_I
       write(*,*) ' alpha_I  =', alpha_I
       write(*,*) ' kin_IIII =', kin_IIII
       write(*,*) ' v_II     =', v_II
       write(*,*) ' ve_II    =', ve_II
       write(*,*) ' Qexc_II  =', Qexc_II
       write(*,*) ' Qion_II  =', Qion_II
       write(*,*) ' IsIntersectedShapeR =', &
            IsIntersectedShapeR
       write(*,*) ' IsIntersectedShape  =', IsIntersectedShape
    end if

  end subroutine user_calc_rates

  !========================================================================

  subroutine user_calc_sources(iBlock)

    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, Dt_BLK
    use ModAdvance,    ONLY: State_VGB, Source_VC, Rho_, RhoUx_, RhoUy_, RhoUz_, &
         Bx_,By_,Bz_, P_
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cProtonMass
    use ModGeometry,   ONLY: r_BLK, Xyz_DGB
    use ModCurrent,    ONLY: get_current
    use ModProcMH,     ONLY: iProc
    use ModPhysics,    ONLY: SW_Ux, SW_Uy, SW_Uz, UnitN_, UnitRho_, UnitU_, UnitP_, UnitT_, UnitB_, &
         ElectronPressureRatio, ElectronCharge, Si2No_V, No2Si_V, UnitEnergyDens_, UnitJ_, UnitRhoU_, &
         UnitTemperature_
    use ModPointImplicit, ONLY: UsePointImplicit, IsPointImplSource
    use ModVarIndexes, ONLY: MassFluid_I
    use ModBlockData, ONLY: use_block_data, &
         get_block_data, put_block_data

    integer, intent(in) :: iBlock

    real, dimension(1:nI,1:nJ,1:nK) :: nElec_C, Te_C, SBx_C, SBy_C, SBz_C, SPe_C
    real, dimension(4,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoTerm_IIC
    real, dimension(5,1:nIonFluid,1:nI,1:nJ,1:nK) :: SRhoUxTerm_IIC, SRhoUyTerm_IIC, SRhoUzTerm_IIC
    real, dimension(8,1:nIonFluid,1:nI,1:nJ,1:nK) :: SPTerm_IIC
    real, dimension(8,1:nI,1:nJ,1:nK) :: SPeTerm_IC

    real, dimension(1:3,1:nI,1:nJ,1:nK) :: Current_DC, uIonMean_DC, uElec_DC
    real, dimension(1:3,1:nIonFluid,1:nI,1:nJ,1:nK) :: uIon_DIC
    real, dimension(1:3,1:nI,1:nJ,1:nK) :: uNeu1_DC
    real, dimension(1:nI,1:nJ,1:nK) :: TempNeu1_C, nNeu1_C
    real, dimension(1:nIonFluid,1:nNeutral,1:nI,1:nJ,1:nK) :: fin_IIC, uIonNeu2_IIC
    real, dimension(1:nNeutral,1:nI,1:nJ,1:nK) :: fen_IC, uNeuElec2_IC
    real, dimension(1:nIonFluid,1:nIonFluid,1:nI,1:nJ,1:nK) :: fii_IIC, uIonIon2_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: Ti_IC, uIonElec2_IC, fei_IC, fie_IC, &
         nIon_IC, SRho_IC, SRhoUx_IC, SRhoUy_IC, SRhoUz_IC, SP_IC
    real, dimension(1:nNeutral,1:nIonFluid) :: Qexc_II, Qion_II
    real, dimension(1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: v_IIC, ve_IIC
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: alpha_IC
    real, dimension(1:nIonFluid) :: fiiTot_I, finTot_I, vAdd_I, kinAdd_I, kinSub_I
    real, dimension(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,1:nI,1:nJ,1:nK) :: kin_IIIIC

    logical :: DoTest, DoTestMe
    real :: theta, fenTot, feiTot,logTe
    integer :: i,j,k,iNeutral,jNeutral,iIonFluid,jIonFluid,iTerm


    logical :: DoCalcShading = .false.
    integer, save :: iBlockLast = -100
    real,    save :: IsIntersectedShapeR_III(nI,nJ,nK) = -1.0

    character(len=*), parameter:: NameSub='user_calc_sources'
    !----------------------------------------------------------------------

    ! Do not evaluate any source terms explicitly when running pointimplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Evaluate source terms explicitly even when running pointimplicit
    !if(UsePointImplicit .and. IsPointImplSource) RETURN

    !! Limit region for evaluation for source term evaluation
    !! if(RMin_BLK(iBlock) > 2.) RETURN

    if(iBlock == BlkTest.and.iProc==ProcTest) then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    !    write(*,*) 'calc_sources: iBlock =', iBlock

    !! Set the source arrays for this block to zero
    SRho_IC        = 0.
    SRhoTerm_IIC   = 0.
    SRhoUx_IC      = 0.
    SRhoUxTerm_IIC = 0.
    SRhoUy_IC      = 0.
    SRhoUyTerm_IIC = 0.
    SRhoUz_IC      = 0.
    SRhoUzTerm_IIC = 0.
    SBx_C          = 0.
    SBy_C          = 0.
    SBz_C          = 0.
    SP_IC          = 0.
    SPTerm_IIC     = 0.
    SPe_C          = 0.
    SPeTerm_IC     = 0.


    ! nElec_C is the electron/ion density in SI units ( n_e=sum(n_i*Zi) )
    do k=1,nK; do j=1,nJ; do i=1,nI
       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    !! ion velocity components in SI
    uIon_DIC(1,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUxIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(2,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUyIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIon_DIC(3,1:nIonFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUzIon_I,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(iRhoIon_I,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uIonMean_DC(1:3,1:nI,1:nJ,1:nK) = 0.

    !! Neu1 velocity componet in SI
    uNeu1_DC(1, 1:nI,1:nJ,1:nK) = State_VGB(Neu1RhoUx_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Neu1Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uNeu1_DC(2, 1:nI,1:nJ,1:nK) = State_VGB(Neu1RhoUy_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Neu1Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)
    uNeu1_DC(3, 1:nI,1:nJ,1:nK) = State_VGB(Neu1RhoUz_,1:nI,1:nJ,1:nK,iBlock) / &
         State_VGB(Neu1Rho_,1:nI,1:nJ,1:nK,iBlock)*No2SI_V(UnitU_)

    ! Neu1 temperature in SI
    TempNeu1_C = State_VGB(Neu1P_,1:nI,1:nJ,1:nK,iBlock)* &
         MassFluid_I(nFluid)/State_VGB(Neu1Rho_,1:nI,1:nJ,1:nK,iBlock) * &
         No2SI_V(UnitTemperature_)

    ! Neu1 density in SI
    nNeu1_C  = State_VGB(Neu1Rho_,1:nI,1:nJ,1:nK,iBlock)/MassFluid_I(nFluid)*No2SI_V(UnitN_)

    do iIonFluid=1,nIonFluid
       uIonMean_DC(1,1:nI,1:nJ,1:nK) = uIonMean_DC(1,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(2,1:nI,1:nJ,1:nK) = uIonMean_DC(2,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid)
       uIonMean_DC(3,1:nI,1:nJ,1:nK) = uIonMean_DC(3,1:nI,1:nJ,1:nK)+nIon_IC(iIonFluid,1:nI,1:nJ,1:nK)* &
            uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)/nElec_C(1:nI,1:nJ,1:nK)*ChargeIon_I(iIonFluid) 
    end do

    !! (u_i-u_n)^2 in SI
    !    do iIonFluid=1,nIonFluid
    !       uIonNeu2_IIC(iIonFluid,iNeutral,1:nI,1:nJ,1:nK) = &
    !            (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)- &
    !            UnxNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
    !            (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)- &
    !            UnyNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2+&
    !            (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)- &
    !            UnzNeutral_IG(iNeutral,2-MinI:nI+1-MinI,2-MinJ:nJ+1-MinJ,2-MinK:nK+1-MinK))**2
    !    end do

    do iIonFluid=1,nIonFluid
       do iNeutral=1,nNeutral
          uIonNeu2_IIC(iIonFluid,iNeutral,1:nI,1:nJ,1:nK) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)-uNeu1_DC(1,1:nI,1:nJ,1:nK))**2 + &
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)-uNeu1_DC(2,1:nI,1:nJ,1:nK))**2 + &
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)-uNeu1_DC(3,1:nI,1:nJ,1:nK))**2
       end do
    end do

    !! (u_i1-u_i2)^2 in SI
    do iIonFluid=1,nIonFluid
       do jIonFluid=1,nIonFluid
          uIonIon2_IIC(iIonFluid,jIonFluid,:,:,:) = &
               (uIon_DIC(1,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(1,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(2,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(2,jIonFluid,1:nI,1:nJ,1:nK))**2+&
               (uIon_DIC(3,iIonFluid,1:nI,1:nJ,1:nK)-uIon_DIC(3,jIonFluid,1:nI,1:nJ,1:nK))**2
       end do
    end do

    if (UseElectronPressure) then
       ! Electron temperature calculated from electron pressure
       ! Ion temperature is calculated from ion pressure
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/&
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
               nElec_C(i,j,k))
       end do; end do; end do
    else
       ! Electron temperature calculated from pressure assuming Te_C=Ti_IC*ElectronTemperatureRatio:
       ! p=nkT with n_e=n_i*Z_i (quasi-neutrality), n=n_e+n_i and p=p_e+p_i=p_i*(1+ElectronPressureRatio)
       do k=1,nK; do j=1,nJ; do i=1,nI
          Ti_IC(1:nIonFluid,i,j,k) = State_VGB(iPIon_I,i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_))
          Te_C(i,j,k) = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
               NO2SI_V(UnitP_)/(cBoltzmann*nElec_C(i,j,k))
       end do; end do; end do
    end if

    if (DoUseCGShape) then
       ! Need to calculate the shading for the photoionization behind
       ! the body for the first source term of the point-implicit solver.
       if (iBlock /= iBlockLast) then
          iBlocklast = iBlock
          if (use_block_data(iBlock)) then
             call get_block_data(iBlock,nI,nJ,nK, IsIntersectedShapeR_III)
             DoCalcShading = .false.
          else
             DoCalcShading = .true.
          end if
       else
          DoCalcShading = .false.
       end if
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI
       ! No need to evaluate source terms for cells inside the body
       ! if((Body1).and.(R_BLK(i,j,k,iBlock)<rBody)) CYCLE

       call get_current(i,j,k,iBlock,Current_DC(:,i,j,k))

       ! calculate uElec_DC from Hall velocity -J/(e*n) [m/s]
       uElec_DC(1:3,i,j,k) = uIonMean_DC(1:3,i,j,k)-Current_DC(1:3,i,j,k)/(nElec_C(i,j,k)*Si2No_V(UnitN_)*&
            ElectronCharge)*No2SI_V(UnitU_)

       call calc_electron_collision_rates(Te_C(i,j,k),nElec_C(i,j,k),i,j,k,iBlock,fen_IC(1:nNeutral,i,j,k), &
            fei_IC(1:nIonFluid,i,j,k))
       call user_calc_rates(Ti_IC(1:nIonFluid,i,j,k),Te_C(i,j,k),i,j,k,iBlock,nElec_C(i,j,k),nIon_IC(1:nIonFluid,i,j,k),&
            fin_IIC(1:nIonFluid,1:nNeutral,i,j,k),fii_IIC(1:nIonFluid,1:nIonFluid,i,j,k),fie_IC(1:nIonFluid,i,j,k),&
            alpha_IC(1:nIonFluid,i,j,k),kin_IIIIC(1:nIonFluid,1:nNeutral,1:nNeutral,1:nIonFluid,i,j,k),&
            v_IIC(1:nNeutral,1:nIonFluid,i,j,k),ve_IIC(1:nNeutral,1:nIonFluid,i,j,k),uElec_DC(1:3,i,j,k),&
            uIon_DIC(1:3,1:nIonFluid,i,j,k),Qexc_II(1:nNeutral,1:nIonFluid),Qion_II(1:nNeutral,1:nIonFluid), &
            DoCalcShading, IsIntersectedShapeR_III(i,j,k))

       !! Zeroth moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinAdd_I = 0. ; kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        nNeu1_C(i,j,k)
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        nNeu1_C(i,j,k)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       !       do iNeutral=1,nNeutral
       !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
       !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
       !       end do
       do iNeutral=1,nNeutral
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               nNeu1_C(i,j,k)
       end do

       !! Sources divideded into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       SRhoTerm_IIC(1,1:nIonFluid,i,j,k) = vAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I              !! newly ionized neutrals
       SRhoTerm_IIC(2,1:nIonFluid,i,j,k) = kinAdd_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I            !! mass added through ion-neutral charge exchange
       SRhoTerm_IIC(3,1:nIonFluid,i,j,k) = -kinSub_I(1:nIonFluid)*Si2No_V(UnitN_)/Si2No_V(UnitT_)*MassIon_I           !! mass removed through ion-neutral charge exchange
       SRhoTerm_IIC(4,1:nIonFluid,i,j,k) = -alpha_IC(1:nIonFluid,i,j,k)*(nElec_C(i,j,k)* &                            !! loss due to recombination
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)/Si2No_V(UnitT_))*MassIon_I


       !! First moment, x component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(1,iIonFluid,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       !       do iNeutral=1,nNeutral
       !          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
       !               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
       !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
       !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
       !               (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,1:nIonFluid,i,j,k))
       !       end do
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (uNeu1_DC(1,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               nNeu1_C(i,j,k)*&
               (uNeu1_DC(1,i,j,k)-uIon_DIC(1,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                   !                        (UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(1,jIonFluid,i,j,k))
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*nNeu1_C(i,j,k)*&
                        (uNeu1_DC(1,i,j,k)-uIon_DIC(1,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUxTerm_IIC(1,iIonFluid,i,j,k) = SRhoUxTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(1,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(1,i,j,k)
       SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(1,1:nIonFluid,i,j,k)-uIonMean_DC(1,i,j,k))*Si2No_V(UnitU_)
       SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, y component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(2,iIonFluid,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          !          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
          !               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
          !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
          !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
          !               (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,1:nIonFluid,i,j,k))
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (uNeu1_DC(2,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               nNeu1_C(i,j,k)*&
               (uNeu1_DC(2,i,j,k)-uIon_DIC(2,1:nIonFluid,i,j,k))
       end do


       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*&
                   !                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                   !                        (UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(2,jIonFluid,i,j,k))
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*&
                        nNeu1_C(i,j,k)*&
                        (uNeu1_DC(2,i,j,k)-uIon_DIC(2,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                       !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUyTerm_IIC(1,iIonFluid,i,j,k) = SRhoUyTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(2,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &             !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(2,i,j,k)
       SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(2,1:nIonFluid,i,j,k)-uIonMean_DC(2,i,j,k))*Si2No_V(UnitU_)
       SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &    !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &    !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       !! First moment, z component
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
       fiiTot_I = 0. ; finTot_I = 0. ; vAdd_I = 0.
       do iIonFluid=1,nIonFluid                                                                                       !! momentum transfer by ion-ion collisions
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*&                        !! ion-ion collisions
               (uIon_DIC(3,iIonFluid,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do                                                                                                         !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          !          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
          !               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
          !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
          !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
          !               (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,1:nIonFluid,i,j,k))
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*&                         !! ion-neutral collisions
               (uNeu1_DC(3,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               nNeu1_C(i,j,k)*&
               (uNeu1_DC(3,i,j,k)-uIon_DIC(3,1:nIonFluid,i,j,k))
       end do

       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                   !                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*&
                   !                        (UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uIon_DIC(3,jIonFluid,i,j,k))
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        nNeu1_C(i,j,k)*&
                        (uNeu1_DC(3,i,j,k)-uIon_DIC(3,jIonFluid,i,j,k))
                end do
             end do
          end do
       end do

       SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k) = (vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                 !! newly photoionized neutrals
            kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*Si2No_V(UnitN_)*MassIon_I* &                                      !! new ions from charge exchange
            Si2No_V(UnitU_)
       ! Add u_s*drho_s/dt for d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt
       do iIonFluid=1,nIonFluid
          SRhoUzTerm_IIC(1,iIonFluid,i,j,k) = SRhoUzTerm_IIC(1,iIonFluid,i,j,k)+sum(SRhoTerm_IIC(1:4,iIonFluid,i,j,k))&
               *uIon_DIC(3,iIonFluid,i,j,k)*Si2No_V(UnitU_)
       end do
       SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)/ElectronCharge* &            !! current dissipation, ion-electron collisions
            MassIon_I*nIon_IC(1:nIonFluid,i,j,k)/nElec_C(i,j,k)*Current_DC(3,i,j,k)
       SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k) = -fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*MassIon_I*&
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitN_)*(uIon_DIC(3,1:nIonFluid,i,j,k)-uIonMean_DC(3,i,j,k))*Si2No_V(UnitU_)
       SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*fiiTot_I(1:nIonFluid)* &   !! ion-ion collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I
       SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k) = nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*finTot_I(1:nIonFluid)* &   !! ion neutral collisions
            Si2No_V(UnitU_)/Si2No_V(UnitT_)*cProtonMass*MassIon_I

       ! (u_n-u_e)^2 difference in neutral and electron speeds qubed [m^2/s^2]
       do iNeutral=1,nNeutral
          !          uNeuElec2_IC(iNeutral,i,j,k) = ((UnxNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(1,i,j,k))**2 &
          !               +(UnyNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(2,i,j,k))**2 &
          !               +(UnzNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-uElec_DC(3,i,j,k))**2)
          uNeuElec2_IC(iNeutral,i,j,k) =( (uNeu1_DC(1,i,j,k)-uElec_DC(1,i,j,k))**2 &
               +(uNeu1_DC(2,i,j,k)-uElec_DC(2,i,j,k))**2 &
               +(uNeu1_DC(3,i,j,k)-uElec_DC(3,i,j,k))**2 )
       end do

       ! (u_i-u_e)^2 difference in ion and electron speeds qubed [m^2/s^2]
       do iIonFluid=1,nIonFluid
          uIonElec2_IC(iIonFluid,i,j,k) = (uIon_DIC(1,iIonFluid,i,j,k)-uElec_DC(1,i,j,k))**2+&
               (uIon_DIC(2,iIonFluid,i,j,k)-uElec_DC(2,i,j,k))**2+&
               (uIon_DIC(3,iIonFluid,i,j,k)-uElec_DC(3,i,j,k))**2
       end do

       !! Second moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"       
       kinSub_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! subtraction to individual fluid from charge exchange [1/(m^3*s)]
                   !                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + &!!nIon_IC(iIonFluid,i,j,k)* &
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)*NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
                   kinSub_I(iIonFluid) = kinSub_I(iIonFluid) + &!!nIon_IC(iIonFluid,i,j,k)* &
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        nNeu1_C(i,j,k)
                end do
             end do
          end do
       end do

       vAdd_I = 0.
       do iNeutral=1,nNeutral
          !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)*&
          !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)*&
               nNeu1_C(i,j,k)
       end do

       SPTerm_IIC(1,1:nIonFluid,i,j,k) = -(kinSub_I(1:nIonFluid)/Si2No_V(UnitT_)+ &                                  !! lost ions through charge exchange and recombination
            alpha_IC(1:nIonFluid,i,j,k)*nElec_C(i,j,k)/Si2No_V(UnitT_))*State_VGB(iPIon_I,i,j,k,iBlock)

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               cBoltzmann*(Ti_IC(iIonFluid,i,j,k)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(2,1:nIonFluid,i,j,k) = 2.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)        

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          !          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
          !               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+NeutralMass_I(nFluid)/cProtonMass)*&
          !               cBoltzmann*(TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Ti_IC(1:nIonFluid,i,j,k))
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)/(MassIon_I(1:nIonFluid)+MassFluid_I(nFluid))*&
               cBoltzmann*(TempNeu1_C(i,j,k)-Ti_IC(1:nIonFluid,i,j,k))
       end do

       SPTerm_IIC(3,1:nIonFluid,i,j,k) = 2.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(4,1:nIonFluid,i,j,k) = 2.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*nIon_IC(1:nIonFluid,i,j,k)*&
            cBoltzmann*(Te_C(i,j,k)-Ti_IC(1:nIonFluid,i,j,k))*Si2No_V(UnitEnergyDens_)
       SPTerm_IIC(5,1:nIonFluid,i,j,k) = 2./3.*fie_IC(1:nIonFluid,i,j,k)/Si2No_V(UnitT_)*cElectronMass*&             !! ion-electron collisional exchange (due to Hall velocity)
            nIon_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitRho_)*uIonElec2_IC(1:nIonFluid,i,j,k)*Si2No_V(UnitU_)**2

       fiiTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-ion collisions
       do iIonFluid=1,nIonFluid
          fiiTot_I(1:nIonFluid) = fiiTot_I(1:nIonFluid)+fii_IIC(1:nIonFluid,iIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*MassIon_I(iIonFluid)/(MassIon_I(1:nIonFluid)+MassIon_I(iIonFluid))*&
               uIonIon2_IIC(1:nIonFluid,iIonFluid,i,j,k)
       end do

       finTot_I(1:nIonFluid) = 0.                                                                                    !! momentum transfer by ion-neutral collisions
       do iNeutral=1,nNeutral
          finTot_I(1:nIonFluid) = finTot_I(1:nIonFluid)+fin_IIC(1:nIonFluid,iNeutral,i,j,k)*nIon_IC(1:nIonFluid,i,j,k)*&  
               MassIon_I(1:nIonFluid)*MassFluid_I(nFluid)/(MassIon_I(1:nIonFluid)+MassFluid_I(nFluid))*&
               uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)
       end do

       SPTerm_IIC(6,1:nIonFluid,i,j,k) = 2./3.*fiiTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2
       SPTerm_IIC(7,1:nIonFluid,i,j,k) = 2./3.*finTot_I(1:nIonFluid)/Si2No_V(UnitT_)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2


       do iNeutral=1,nNeutral
          !          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
          !               NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)
          vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
               nNeu1_C(i,j,k)* &
               uIonNeu2_IIC(1:nIonFluid,iNeutral,i,j,k)
       end do
       kinAdd_I = 0.
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             do iNeutral=1,nNeutral
                do jNeutral=1,nNeutral
                   !! addition to individual fluid from charge exchange [1/(m*s^2)]
                   !                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)*&
                   !                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                   !                        NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uIonNeu2_IIC(jIonFluid,iNeutral,i,j,k)
                   kinAdd_I(jIonFluid) = kinAdd_I(jIonFluid) + nIon_IC(iIonFluid,i,j,k)*&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)* &
                        nNeu1_C(i,j,k)* &
                        uIonNeu2_IIC(jIonFluid,iNeutral,i,j,k)
                end do
             end do
          end do
       end do

       SPTerm_IIC(8,1:nIonFluid,i,j,k) = 1./3.*(vAdd_I(1:nIonFluid)/Si2No_V(UnitT_)+kinAdd_I(1:nIonFluid)/Si2No_V(UnitT_))*&
            MassIon_I(1:nIonFluid)*Si2No_V(UnitN_)*Si2No_V(UnitU_)**2

       if (UseElectronPressure) then
          SPeTerm_IC(1,i,j,k) = -sum(alpha_IC(1:nIonFluid,i,j,k)*nIon_IC(1:nIonFluid,i,j,k))/ &                           !! lost electrons through recombination
               Si2No_V(UnitT_)*State_VGB(Pe_,i,j,k,iBlock)

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             !             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
             !                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)*uNeuElec2_IC(iNeutral,i,j,k)
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+v_IIC(iNeutral,1:nIonFluid,i,j,k)* &
                  nNeu1_C(i,j,k)* &
                  uNeuElec2_IC(iNeutral,i,j,k)
          end do
          SPeTerm_IC(2,i,j,k) = 1./3.*cElectronMass*sum(vAdd_I)*Si2No_V(UnitRho_)/Si2No_V(UnitT_)*Si2No_V(UnitU_)**2      !! new electrons through photoionized neutrals

          feiTot = 0.
          do iIonFluid=1,nIonFluid
             feiTot = feiTot+fei_IC(iIonFluid,i,j,k)/MassIon_I(iIonFluid)*&
                  (Ti_IC(iIonFluid,i,j,k)-Te_C(i,j,k))
          end do
          SPeTerm_IC(3,i,j,k) = 2.*cElectronMass*Si2No_V(UnitRho_)/Si2No_V(UnitN_)*&                                      !! ion-electron collisional exchange (thermal motion)
               nElec_C(i,j,k)*cBoltzmann*Si2No_V(UnitEnergyDens_)*feiTot/Si2No_V(UnitT_)

          fenTot = 0.
          do iNeutral=1,nNeutral
             !             fenTot = fenTot+fen_IC(iNeutral,i,j,k)/MassFluid_I(nFluid)*&
             !                  (TnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)-Te_C(i,j,k))
             fenTot = fenTot+fen_IC(iNeutral,i,j,k)/MassFluid_I(nFluid)/cProtonMass*&
                  (TempNeu1_C(i,j,k)-Te_C(i,j,k))
          end do
          SPeTerm_IC(4,i,j,k) = 2.*cElectronMass*nElec_C(i,j,k)*cBoltzmann*&                                              !! electron-neutral collisional exchange (thermal motion)
               Si2No_V(UnitEnergyDens_)*fenTot/Si2No_V(UnitT_)

          SPeTerm_IC(5,i,j,k) = 2./3.*sum(fei_IC(1:nIonFluid,i,j,k)*uIonElec2_IC(1:nIonFluid,i,j,k))/ &                   !! ion-electron collisional exchange (due to Hall velocity)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          SPeTerm_IC(6,i,j,k) = 2./3.*sum(fen_IC(1:nNeutral,i,j,k)*uNeuElec2_IC(1:nNeutral,i,j,k))/&                      !! electron-neutral collisional exchange (bulk motion)
               Si2No_V(UnitT_)*cElectronMass*nElec_C(i,j,k)*Si2No_V(UnitRho_)*Si2No_V(UnitU_)**2

          vAdd_I(1:nIonFluid) = 0.
          do iNeutral=1,nNeutral
             !             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+((v_IIC(iNeutral,1:nIonFluid,i,j,k)-&
             !                  ve_IIC(iNeutral,1:nIonFluid,i,j,k))*Qexc_II(iNeutral,1:nIonFluid)- &
             !                  ve_IIC(iNeutral,1:nIonFluid,i,j,k)*Qion_II(iNeutral,1:nIonFluid))*ChargeIon_I(1:nIonFluid)* &          
             !                  NnNeutral_IG(iNeutral,i-MinI+1,j-MinJ+1,k-MinK+1)
             vAdd_I(1:nIonFluid) = vAdd_I(1:nIonFluid)+((v_IIC(iNeutral,1:nIonFluid,i,j,k)-&
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k))*Qexc_II(iNeutral,1:nIonFluid)- &
                  ve_IIC(iNeutral,1:nIonFluid,i,j,k)*Qion_II(iNeutral,1:nIonFluid))*ChargeIon_I(1:nIonFluid)* &          
                  nNeu1_C(i,j,k)
          end do
          SPeTerm_IC(7,i,j,k) = 2./3.*sum(vAdd_I)*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                                ! heating of electrons due to ionization excess energy

          logTe = log(Te_C(i,j,k))
          SPeTerm_IC(8,i,j,k) = exp(-188.4701+33.2547*logTe-2.0792*logTe**2+0.0425*logTe**3)                              !! electron cooling due to collisions w/ water vapor
          !          if(Te_C(i,j,k)<1.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)) then
          !             SPeTerm_IC(8,i,j,k)=4.5e-9/(0.5*TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))* &
          !                  (Te_C(i,j,k)-TnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1))
          !          else
          !             SPeTerm_IC(8,i,j,k)=SPeTerm_IC(8,i,j,k)+4.5e-9
          !          end if
          if(Te_C(i,j,k)<1.5*TempNeu1_C(i,j,k)) then
             SPeTerm_IC(8,i,j,k)=4.5e-9/(0.5*TempNeu1_C(i,j,k))* &
                  (Te_C(i,j,k)-TempNeu1_C(i,j,k))
          else
             SPeTerm_IC(8,i,j,k)=SPeTerm_IC(8,i,j,k)+4.5e-9
          end if

          !          SPeTerm_IC(8,i,j,k) = -2./3.*NnNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*nElec_C(i,j,k)* &
          !               SPeTerm_IC(8,i,j,k)/1e6*1.60217733e-19*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                           
          SPeTerm_IC(8,i,j,k) = -2./3.*nNeu1_C(i,j,k)*nElec_C(i,j,k)* &
               SPeTerm_IC(8,i,j,k)/1e6*1.60217733e-19*Si2No_V(UnitEnergyDens_)/Si2No_V(UnitT_)                           

       end if

       !! sum up individual terms
       do iTerm=1,4
          SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(4,1:nIonFluid,i,j,k)
       do iTerm=1,5
          SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(iTerm,1:nIonFluid,i,j,k)
          SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SRhoUx_IC(1:nIonFluid,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)+SRhoUxTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUy_IC(1:nIonFluid,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)+SRhoUyTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SRhoUz_IC(1:nIonFluid,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)+SRhoUzTerm_IIC(5,1:nIonFluid,i,j,k)
       do iTerm=1,8
          SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(iTerm,1:nIonFluid,i,j,k)
       end do
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(1,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(2,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(3,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(4,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(5,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(6,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(7,1:nIonFluid,i,j,k)
       ! SP_IC(1:nIonFluid,i,j,k) = SP_IC(1:nIonFluid,i,j,k)+SPTerm_IIC(8,1:nIonFluid,i,j,k)
       if(UseElectronPressure) then
          SPe_C(i,j,k) = sum(SPeTerm_IC(1:8,i,j,k))
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(1,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(2,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(3,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(4,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(5,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(6,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(7,i,j,k)
          ! SPe_C(i,j,k) = SPe_C(i,j,k) + SPeTerm_IC(8,i,j,k)
       end if

       Source_VC(iRhoIon_I   ,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)    + Source_VC(iRhoIon_I   ,i,j,k)
       Source_VC(iRhoUxIon_I ,i,j,k) = SRhoUx_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUxIon_I ,i,j,k)
       Source_VC(iRhoUyIon_I ,i,j,k) = SRhoUy_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUyIon_I ,i,j,k)
       Source_VC(iRhoUzIon_I ,i,j,k) = SRhoUz_IC(1:nIonFluid,i,j,k)  + Source_VC(iRhoUzIon_I ,i,j,k)
       Source_VC(iPIon_I     ,i,j,k) = SP_IC(1:nIonFluid,i,j,k)      + Source_VC(iPIon_I     ,i,j,k)

       Source_VC(Rho_   ,i,j,k) = sum(SRho_IC(1:nIonFluid,i,j,k))    + Source_VC(Rho_   ,i,j,k)
       Source_VC(rhoUx_ ,i,j,k) = sum(SRhoUx_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUx_ ,i,j,k)
       Source_VC(rhoUy_ ,i,j,k) = sum(SRhoUy_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUy_ ,i,j,k)
       Source_VC(rhoUz_ ,i,j,k) = sum(SRhoUz_IC(1:nIonFluid,i,j,k))  + Source_VC(rhoUz_ ,i,j,k)
       Source_VC(Bx_    ,i,j,k) = SBx_C(i,j,k)                       + Source_VC(Bx_    ,i,j,k)
       Source_VC(By_    ,i,j,k) = SBy_C(i,j,k)                       + Source_VC(By_    ,i,j,k)
       Source_VC(Bz_    ,i,j,k) = SBz_C(i,j,k)                       + Source_VC(Bz_    ,i,j,k)
       if(UseElectronPressure) then
          if (DoTestMe) then
             write(*,*) NameSub, ': Source_VC(Pe_) before applying user term =', Source_VC(Pe_    ,i,j,k)
             write(*,*) NameSub, ': User source_VC(Pe_)                      =', SPe_C(i,j,k)
          end if
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))   + Source_VC(P_     ,i,j,k)
          Source_VC(Pe_    ,i,j,k) = SPe_C(i,j,k)                    + Source_VC(Pe_    ,i,j,k)
       else
          Source_VC(P_     ,i,j,k) = sum(SP_IC(1:nIonFluid,i,j,k))*(1.+ElectronPressureRatio) + &
               Source_VC(P_     ,i,j,k)
       end if

    end do;  end do;  end do

    if (DoCalcShading .and. DoUseCGShape) then
       call put_block_data(iBlock,nI,nJ,nK,IsIntersectedShapeR_III)
    end if

    if(DoTestMe) then
       write(*,*) NameSub
       write(*,*)'Inputs: '
       i=iTest ; j=jTest ; k=kTest
       theta=acos((-SW_Ux*Xyz_DGB(x_,i,j,k,iBlock)-SW_Uy*Xyz_DGB(y_,i,j,k,iBlock)&
            -SW_Uz*Xyz_DGB(z_,i,j,k,iBlock))/R_BLK(i,j,k,iBlock)/&
            (SW_Ux**2+SW_Uy**2+SW_Uz**2)**0.5)
123    format (A13,ES25.16,A15,A3,F7.2,A3)
       write(*,123)'x         = ',Xyz_DGB(x_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'y         = ',Xyz_DGB(y_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'z         = ',Xyz_DGB(z_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'r         = ',R_BLK(i,j,k,iBlock)," [rPlanet]"
       write(*,123)'SW_Ux     = ',SW_Ux*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uy     = ',SW_Uy*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'SW_Uz     = ',SW_Uz*No2SI_V(UnitU_)," [m/s]"
       write(*,123)'Tmin      = ',Tmin," [K]"
       write(*,*)''
       write(*,*)'Neutrals:'
       do iNeutral=1,nNeutral
          write(*,124)'Neutral species #',iNeutral,': ', NameNeutral_I(iNeutral)," (",&
               MassFluid_I(nFluid)," amu)"
          write(*,123)'n_n       = ',nNeu1_C(i,j,k)," [m^-3]"
          write(*,123)'m_n       = ',MassFluid_I(nFluid)," [amu]"
          write(*,123)'unx       = ',uNeu1_DC(1,i,j,k)," [m/s]"
          write(*,123)'uny       = ',uNeu1_DC(2,i,j,k)," [m/s]"
          write(*,123)'unz       = ',uNeu1_DC(3,i,j,k)," [m/s]"
          write(*,123)'Tn        = ',TempNeu1_C(i,j,k)," [K]"
       end do
       write(*,*)''
       write(*,*)'Total plasma phase (e- and i+):'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'uRhox     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       write(*,123)'uRhoz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_)," [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'Ptot      = ',(State_VGB(P_,i,j,k,iBlock)+State_VGB(Pe_,i,j,k,iBlock))*&
               No2SI_V(UnitP_)," [kg/(m*s^2)]"
       else
          write(*,123)'Ptot      = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [kg/(m*s^2)]"
       end if
       write(*,123)'Bx        = ',State_VGB(Bx_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'By        = ',State_VGB(By_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'Bz        = ',State_VGB(Bz_,i,j,k,iBlock)*No2SI_V(UnitB_)," [T]"
       write(*,123)'uMeanx    = ',uIonMean_DC(1,i,j,k)," [m/s]"
       write(*,123)'uMeany    = ',uIonMean_DC(2,i,j,k)," [m/s]"
       write(*,123)'uMeanz    = ',uIonMean_DC(3,i,j,k)," [m/s]"       
       write(*,123)'jx        = ',Current_DC(1,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jy        = ',Current_DC(2,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,123)'jz        = ',Current_DC(3,i,j,k)*No2SI_V(UnitJ_)," [A/m^2]"
       write(*,*)''
       write(*,123)'SRho      = ',sum(SRho_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
            100.*sum(SRho_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(Rho_,i,j,k,iBlock)),"%)"
       if (State_VGB(RhoUx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUx_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUx    = ',sum(SRhoUx_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUy_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUy_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUy_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUy    = ',sum(SRhoUy_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(RhoUz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
               " (", 100.*sum(SRhoUz_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(RhoUz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SRhoUz    = ',sum(SRhoUz_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
       end if
       if (State_VGB(Bx_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBx_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bx_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBx       = ',SBx_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(By_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBy_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(By_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBy       = ',SBy_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if (State_VGB(Bz_,i,j,k,iBlock) /= 0.) then
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"," (", &
               100.*SBz_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Bz_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SBz       = ',SBz_C(i,j,k)*No2SI_V(UnitB_)/No2SI_V(UnitT_)," [T/s]"
       end if
       if(UseElectronPressure) then
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/(State_VGB(P_,i,j,k,iBlock)),"%)"
       else
          write(*,123)'SP        = ',sum(SP_IC(1:nIonFluid,i,j,k))*No2SI_V(UnitP_)/No2SI_V(UnitT_)*&
               (1+ElectronPressureRatio)," [Pa/s]"," (",100.*sum(SP_IC(1:nIonFluid,i,j,k))*Dt_BLK(iBlock)/ &
               (State_VGB(P_,i,j,k,iBlock))*(1+ElectronPressureRatio),"%)"
       end if
       write(*,*)''
       write(*,123)'dt        = ',Dt_BLK(iBlock)*No2SI_V(UnitT_)," [s]"
       write(*,*)''
       write(*,*)'Individual ion fluids:'
       do iIonFluid=1,nIonFluid
          write(*,124)'Ion species     #',iIonFluid,': ',NameFluid_I(iIonFluid+1)," (",&
               MassIon_I(iIonFluid)," amu/",ChargeIon_I(iIonFluid)," e)"
124       format (A17,I2,A3,A7,A3,F5.2,A5,F5.1,A3)
          write(*,123)'Ux        = ',uIon_DIC(1,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uy        = ',uIon_DIC(2,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'Uz        = ',uIon_DIC(3,iIonFluid,i,j,k)," [m/s]"
          write(*,123)'ni        = ',nIon_IC(iIonFluid,i,j,k)," [m^-3]"
          write(*,123)'Ti        = ',Ti_IC(iIonFluid,i,j,k)," [K]"
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'SRho      = ',SRho_IC(iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRho_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT1   = ',SRhoTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT2   = ',SRhoTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT3   = ',SRhoTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SRhoT4   = ',SRhoTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRho_)/No2SI_V(UnitT_)," [kg/(m^3*s)]"," (", &
               100.*SRhoTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          if (State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUx_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUx    = ',SRhoUx_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT1 = ',SRhoUxTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT2 = ',SRhoUxTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT3 = ',SRhoUxTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT4 = ',SRhoUxTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUxT5 = ',SRhoUxTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUy_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (", 100.*SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUy    = ',SRhoUy_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT1 = ',SRhoUyTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT2 = ',SRhoUyTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT3 = ',SRhoUyTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT4 = ',SRhoUyTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUyT5 = ',SRhoUyTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          if (State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) /= 0.) then
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUz_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]",&
                  " (",100.*SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          else
             write(*,123)'SRhoUz    = ',SRhoUz_IC(iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT1 = ',SRhoUzTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT2 = ',SRhoUzTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT3 = ',SRhoUzTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT4 = ',SRhoUzTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
             write(*,123)' SRhoUzT5 = ',SRhoUzTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
          end if
          write(*,123)'SP        = ',SP_IC(iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SP_IC(iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT1     = ',SPTerm_IIC(1,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(1,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT2     = ',SPTerm_IIC(2,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(2,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT3     = ',SPTerm_IIC(3,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(3,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT4     = ',SPTerm_IIC(4,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(4,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT5     = ',SPTerm_IIC(5,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(5,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT6     = ',SPTerm_IIC(6,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(6,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT7     = ',SPTerm_IIC(7,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(7,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
          write(*,123)' SPT8     = ',SPTerm_IIC(8,iIonFluid,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
               100.*SPTerm_IIC(8,iIonFluid,i,j,k)*Dt_BLK(iBlock)/(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)),"%)"
       end do
       write(*,*)''
       write(*,*)'Electrons:'
       write(*,123)'n_e       = ',nElec_C(i,j,k)," [m^-3]"
       if (UseElectronPressure) then
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if
       write(*,123)'Uex       = ',uElec_DC(1,i,j,k)," [m/s]"
       write(*,123)'Uey       = ',uElec_DC(2,i,j,k)," [m/s]"
       write(*,123)'Uez       = ',uElec_DC(3,i,j,k)," [m/s]"
       write(*,123)'Te        = ',Te_C(i,j,k)," [K]"
       if(UseElectronPressure) then
          if (State_VGB(Pe_,i,j,k,iBlock).gt.0.) then
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPe_C(i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(1,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(2,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(3,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(4,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(5,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(6,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(7,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
             write(*,123)' SPeT8    = ',SPeTerm_IC(8,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"," (", &
                  100.*SPeTerm_IC(8,i,j,k)*Dt_BLK(iBlock)/(State_VGB(Pe_,i,j,k,iBlock)),"%)"
          else
             write(*,123)'SPe       = ',SPe_C(i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT1    = ',SPeTerm_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT2    = ',SPeTerm_IC(2,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT3    = ',SPeTerm_IC(3,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT4    = ',SPeTerm_IC(4,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT5    = ',SPeTerm_IC(5,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT6    = ',SPeTerm_IC(6,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT7    = ',SPeTerm_IC(7,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
             write(*,123)' SPeT8    = ',SPeTerm_IC(8,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
          end if
       end if
       write(*,*)''
       write(*,*)'Ion-electron combinations:'
       do iIonFluid=1,nIonFluid
          write(*,*)NameFluid_I(iIonFluid+1), '&  e'
          write(*,123)'fei       = ',fei_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'fie       = ',fie_IC(iIonFluid,i,j,k)," [1/s]"
          write(*,123)'|u_ime|   = ',sqrt(uIonElec2_IC(iIonFluid,i,j,k))," [m/s]"
          write(*,123)'alpha     = ',alpha_IC(iIonFluid,i,j,k)," [m^3/s]"
       end do
       write(*,*)''
       write(*,*)'Ion-ion combinations:'
       do iIonFluid=1,nIonFluid
          do jIonFluid=1,nIonFluid
             if (iIonFluid.ne.jIonFluid) then
                write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameFluid_I(jIonFluid+1)
                write(*,123)' fii      = ',fii_IIC(iIonFluid,jIonFluid,i,j,k)," [1/s]"
                write(*,123)' |u_imi|  = ',sqrt(uIonIon2_IIC(iIonFluid,jIonFluid,i,j,k))," [m/s]"
             end if
          end do
       end do
       write(*,*)''
       write(*,*)'Ion-neutral combinations:'
       do iIonFluid=1,nIonFluid
          do iNeutral=1,nNeutral
             write(*,*)NameFluid_I(iIonFluid+1), '&  ',NameNeutral_I(iNeutral)
             write(*,123)' v_phio   = ',v_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' v_eio    = ',ve_IIC(iNeutral,iIonFluid,i,j,k)," [1/s]"
             write(*,123)' fin      = ',fin_IIC(iIonFluid,iNeutral,i,j,k)," [1/s]"
             write(*,123)' |u_imn|  = ',sqrt(uIonNeu2_IIC(iIonFluid,iNeutral,i,j,k))," [m/s]"
             write(*,*)' kin (Ion & Neutral-> Neutral & Ion):'
             do jIonFluid=1,nIonFluid
                do jNeutral=1,nNeutral
                   write(*,*)' ',NameFluid_I(iIonFluid+1),'&  ',NameNeutral_I(iNeutral),'->  ', &
                        NameNeutral_I(jNeutral),'&  ',NameFluid_I(jIonFluid+1),'=',&
                        kin_IIIIC(iIonFluid,iNeutral,jNeutral,jIonFluid,i,j,k)," [m^3/s]"  
                end do
             end do
          end do
       end do
       write(*,*)''
       write(*,*)'Electron-neutral combinations:'
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral)
          write(*,123)'fen      = ',fen_IC(iNeutral,i,j,k)," [1/s]"
          write(*,123)'|u_nme|  = ',sqrt(uNeuElec2_IC(iNeutral,i,j,k))," [m/s]"
       end do
       write(*,*)''
    end if

  end subroutine user_calc_sources

  !========================================================================

  subroutine user_update_states(iStage,iBlock)
    use ModAdvance,  ONLY: State_VGB, Energy_GBI
    use ModPhysics,  ONLY: SW_N, LowDensityRatio, cBoltzmann, ElectronPressureRatio, Si2No_V, &
         No2Si_V, UnitN_, UnitP_!, UnitB_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: true_cell
    use ModMain,     ONLY: ProcTest, BlkTest, iTest, jTest, kTest, VarTest

    integer,intent(in) :: iStage, iBlock

    integer :: i,j,k,iIonFluid
    logical :: DoTest, DoTestMe   
    real, dimension(1:nI,1:nJ,1:nK)             :: nElec_C
    real, dimension(1:nIonFluid,1:nI,1:nJ,1:nK) :: nIon_IC
    
    character(len=*), parameter :: NameSub='user_update_states'

    !----------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if
    
    call update_states_MHD(iStage,iBlock)

    ! Enforce minimum temperature (pressure), Tmin, if temperatures Ti_IC or Te_C are below

    do k=1,nK; do j=1,nJ; do i=1,nI
       if (.not. true_cell(i,j,k,iBlock)) CYCLE

       do iIonFluid=1,nIonFluid
          ! set minimum mass density (and in these locations Ti = Tmin and vi=vbulkplasma)
          if(State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) < SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2) then
             State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) = SW_n*MassIon_I(iIonFluid)*LowDensityRatio**2
             State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)/ &
                  MassIon_I(iIonFluid)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)

             !! Fix SW
             State_VGB(iRhoUxIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iPIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock)/ &
                  MassIon_I(SW_)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)

          end if
       end do

       ! ! fix solar wind inside cavity to minimum value
       ! if(sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)*No2Si_V(UnitB_)**2 < 1e-9**2) then
       ! if(R_BLK(i,j,k,iBlock)<2.e-5) then
       !    State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) = SW_n*MassIon_I(SW_)*LowDensityRatio**2
       !    State_VGB(iRhoUxIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iRhoUyIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iRhoUzIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
       !         State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       !    State_VGB(iPIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock)/ &
       !         MassIon_I(SW_)*No2SI_V(UnitN_)*cBoltzmann*Tmin*SI2No_V(UnitP_)          
       ! end if

       State_VGB(Rho_,i,j,k,iBlock) = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))

       nIon_IC(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
       nElec_C(i,j,k) = sum(nIon_IC(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))

       do iIonFluid=1,nIonFluid
          ! set minimum pressure
          if(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*NO2SI_V(UnitP_) < &
               nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin) then
             State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock) = &
                  nIon_IC(iIonFluid,i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
          if(State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*NO2SI_V(UnitP_) < &
               nIon_IC(iIonFluid,i,j,k)*cBoltzmann*1000.0) then
             !! Fix SW
             State_VGB(iRhoUxIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUyIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iRhoUzIon_I(SW_),i,j,k,iBlock) = State_VGB(iRhoIon_I(SW_),i,j,k,iBlock) * &
                  State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
             State_VGB(iPIon_I(SW_),i,j,k,iBlock) = &
                  nIon_IC(SW_,i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if

       end do

       if(UseElectronPressure) then
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          if (State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_) < nElec_C(i,j,k)*cBoltzmann*Tmin) then
             State_VGB(Pe_,i,j,k,iBlock) = nElec_C(i,j,k)*cBoltzmann*Tmin*SI2No_V(UnitP_)
          end if
       else
          State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
       end if


    end do; end do; end do

    call calc_energy_cell(iBlock)

    if (DoTestMe) &
         write(*,*) NameSub, ' after user term =', &
         State_VGB(VarTest, iTest, jTest, kTest, iBlock), &
         Energy_GBI(iTest,jTest,kTest,iBlock,:)

  end subroutine user_update_states

  !========================================================================

  subroutine derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I, nElec, EtaSi)
    use ModResistivity,  ONLY: Eta0SI
    use ModConst,        ONLY: cElectronMass, cElectronCharge, cMu
    use ModMain,         ONLY: iTest, jTest, kTest, BlkTest, ProcTest
    use ModProcMH,       ONLY: iProc

    integer, intent(in)  :: iBlock, i, j, k
    real,    intent(in)  :: TeSI
    real,    intent(in)  :: nIon_I(1:nIonFluid)
    real,    intent(in)  :: nElec
    real,    intent(out) :: EtaSi

    real :: EtaSiColl!, EtaSiSpitzer, lnL
    !    real, save :: SpitzerCoef, EtaPerpSpitzerSi
    logical :: DoTest, DoTestMe=.true.
    real :: eeSigma!, B0_D(3)
    real, dimension(nIonFluid) :: fei_I, eiSigma_I
    real, dimension(nNeutral)  :: fen_I, enSigma_I
    integer :: iIonFluid, iNeutral

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('derive_cell_diffusivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! Spitzer formulation from Stoecker "Taschenbuch der Physik", Verlag "Harri Deutsch"
    ! lnL = log(1e7*TeSI**1.5/sqrt(nElec))
    ! EtaSiSpitzer = cElectronCharge**2*lnL/(32.*sqrt(2*cPi/cElectronMass*(cBoltzmann*TeSI)**3)*cEps**2)/cMu


    !! Collisional type resisitivity/diffusivity
    call calc_electron_collision_rates(TeSI,nElec,i,j,k,iBlock,fen_I(1:nNeutral),fei_I(1:nIonFluid))
    eiSigma_I(1:nIonFluid) = cElectronCharge**2*nElec/((fei_I(1:nIonFluid)+1E-20)*cElectronMass) 
    enSigma_I(1:nNeutral) = cElectronCharge**2*nElec/((fen_I(1:nNeutral)+1E-20)*cElectronMass)
    !! Eta_G is calculated from both conductivities using Kirchhoff's rule:
    !! 1/sigma_tot = 1/eiSigma_I+1/enSigma_I
    !! The resulting conductivity is close to Spitzer conductivity far from the comet and
    !! decreases due to abundant electron-neutral collisions close to the nucleus
    !! EtaSiColl = 1/(sigma_tot*mu_0) magnetic diffusivity [m^2/s]
    EtaSiColl = (sum(1/eiSigma_I(1:nIonFluid))+sum(1/enSigma_I(1:nNeutral)))/cMu
    !! Total diffusivity [m^2/s]
    EtaSi = Eta0SI + EtaSiColl

    ! TestArray(1,i,j,k,iBlock) = EtaSiColl
    ! TestArray(2,i,j,k,iBlock) = EtaSiSpitzer
    ! TestArray(3,i,j,k,iBlock) = EtaSiSpitzer/EtaSiColl

    if(DoTestMe) then
       write(*,*)'derive_cell_diffusivity:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       do iIonFluid=1,nIonFluid
          write(*,*)'e & ',NameFluid_I(iIonFluid+1),':'
          write(*,*)'s_ei  = ',eiSigma_I(iIonFluid)," [1/(Ohm*m)]"
       end do
       do iNeutral=1,nNeutral
          write(*,*)'e & ',NameNeutral_I(iNeutral),':'
          write(*,*)'s_en  = ',enSigma_I(iNeutral)," [1/(Ohm*m)]"
       end do
       write(*,*)'e & e:'
       write(*,*)'s_ee  = ',eeSigma," [1/(Ohm*m)]"
       write(*,*)''
       write(*,*)'Eta0   = ',Eta0Si," [m^2/s]"
       write(*,*)'Eta_en = ',sum(1/enSigma_I(1:nNeutral))/cMu," [m^2/s]"
       write(*,*)'Eta_ei = ',sum(1/eiSigma_I(1:nIonFluid))/cMu," [m^2/s]"
       write(*,*)'Eta_ee = ',1/eeSigma/cMu," [m^2/s]"
       write(*,*)'Eta_eX = ',EtaSiColl," [m^2/s]"
       write(*,*)'EtaTot = ',EtaSi," [m^2/s]"
       write(*,*)''
    end if

  end subroutine derive_cell_diffusivity

  !========================================================================

  subroutine user_set_resistivity(iBlock, Eta_G)

    use ModPhysics,     ONLY: No2Si_V, Si2No_V, &
         UnitN_, UnitX_, UnitT_, UnitP_, ElectronPressureRatio
    use ModProcMH,      ONLY: iProc
    use ModMain,        ONLY: ProcTest, BlkTest, iTest, jTest, kTest
    use ModAdvance,     ONLY: State_VGB
    use ModVarIndexes,  ONLY: Pe_, P_
    use ModConst,       ONLY: cBoltzmann
    use ModMultiFluid,  ONLY: MassIon_I

    integer, intent(in) :: iBlock
    real, intent(out) :: Eta_G(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) 

    integer :: i, j, k
    logical :: DoTest, DoTestMe=.true.
    real, dimension(1:nIonFluid,MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: nIon_IG
    real, dimension(MinI:MaxI,MinJ:MaxJ,MinK:MaxK) :: Te_G, nElec_G
    real :: EtaSi

    !---------------------------------------------------------------------

    if(iProc==PROCtest .and. iBlock == BlkTest) then
       call set_oktest('user_set_resistivity',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    ! nElec_G is the electron/ion density in SI units (n_e=n_itot)
    nElec_G = 0.
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       nIon_IG(1:nIonFluid,i,j,k) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
       nElec_G(i,j,k) = sum(nIon_IG(1:nIonFluid,i,j,k)*ChargeIon_I(1:nIonFluid))
    end do; end do; end do

    if (UseElectronPressure) then
       Te_G = State_VGB(Pe_,:,:,:,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann* &
            nElec_G)
    else
       Te_G(:,:,:) = State_VGB(P_,:,:,:,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec_G(:,:,:))
    end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI

       call derive_cell_diffusivity(iBlock, i, j, k, Te_G(i,j,k), nIon_IG(1:nIonFluid,i,j,k), nElec_G(i,j,k), EtaSi)
       Eta_G(i,j,k) = EtaSi*SI2No_V(UnitX_)**2/SI2No_V(UnitT_)

    end do; end do; end do

    if(DoTestMe) then
       write(*,*)'user_set_resistivity:'
       write(*,*)'Te    = ',Te_G(iTest,jTest,kTest)," [K]"
       write(*,*)'n_e   = ',nElec_G(iTest,jTest,kTest)," [m^-3]"
       write(*,*)'Eta   = ',Eta_G(iTest,jTest,kTest)*No2SI_V(UnitX_)**2/No2SI_V(UnitT_)," [m^2/s]"
    end if

  end subroutine user_set_resistivity

  !========================================================================

  subroutine user_material_properties(State_V, i,j,k,iBlock,iDir, &
       EinternalIn, TeIn, NatomicOut, AverageIonChargeOut, &
       EinternalOut, TeOut, PressureOut,   &
       CvOut, GammaOut, HeatCondOut, IonHeatCondOut, TeTiRelaxOut, &
       OpacityPlanckOut_W, OpacityRosselandOut_W, PlanckOut_W, &
       EntropyOut)

    use ModPhysics,      ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitX_, ElectronPressureRatio, inv_gm1
    use ModVarIndexes,   ONLY: nVar, p_
    use ModConst,        ONLY: cElectronCharge, cBoltzmann, cMu, cElectronMass
    use ModAdvance,      ONLY: State_VGB
    use ModMain,         ONLY: ProcTest, iTest, jTest, kTest, BlkTest
    use ModProcMH,       ONLY: iProc
    use ModGeometry,     ONLY: Xyz_DGB

    !------------------------------------------------------------------------
    ! The State_V vector is in normalized units
    real, intent(in) :: State_V(nVar)
    integer, optional, intent(in) :: i, j, k, iBlock, iDir
    real, optional, intent(in)  :: EinternalIn                  ! [J/m^3]
    real, optional, intent(in)  :: TeIn                         ! [K]
    real, optional, intent(out) :: NatomicOut                   ! [1/m^3]
    real, optional, intent(out) :: AverageIonChargeOut          ! dimensionless
    real, optional, intent(out) :: EinternalOut                 ! [J/m^3]
    real, optional, intent(out) :: TeOut                        ! [K]
    real, optional, intent(out) :: PressureOut                  ! [Pa]   
    real, optional, intent(out) :: CvOut                        ! [J/(K*m^3)]  
    real, optional, intent(out) :: GammaOut                     ! dimensionless
    real, optional, intent(out) :: HeatCondOut                  ! [W/(m*K)]   
    real, optional, intent(out) :: IonHeatCondOut               ! [W/(m*K)]
    real, optional, intent(out) :: TeTiRelaxOut                 ! [1/s]  
    real, optional, intent(out) :: OpacityPlanckOut_W(nWave)    ! [1/m] 
    real, optional, intent(out) :: OpacityRosselandOut_W(nWave) ! [1/m] 
    real, optional, intent(out) :: PlanckOut_W(nWave)           ! [J/m^3] 
    real, optional, intent(out) :: EntropyOut

    real, save :: KappaCoeffSI = (cBoltzmann/cElectronCharge)**2/cMu
    real :: nElec, EtaSI, TeSI
    real, dimension(nIonFluid) :: nIon_I
    logical :: DoTest, DoTestMe=.true.

    real :: xmin, xmax, HeatCondFactor, widthmax, widthmin, xMaxyz, xMinyz

    !----------------------------------------------------------------------

    if(iBlock==BlkTest.and.i==iTest.and.j==jTest.and.k==kTest.and.iProc==ProcTest) then
       call set_oktest('user_material_properties',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*NO2SI_V(UnitN_)
    nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
    if (UseElectronPressure) then
       TeSI = State_VGB(Pe_,i,j,k,iBlock)*NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    else
       TeSI = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/(1.+ElectronPressureRatio)*&
            NO2SI_V(UnitP_)/(cBoltzmann*nElec)
    end if
    if(present(CvOut)) CvOut = cBoltzmann*nElec*inv_gm1
    if(present(TeOut)) TeOut = TeSI
    if(present(AverageIonChargeOut).or.present(NatomicOut)) AverageIonChargeOut = nElec/sum(nIon_I(1:nIonFluid))
    if(present(NatomicOut)) NatomicOut = nElec/AverageIonChargeOut

    if(present(HeatCondOut)) then
       !!write(*,*)'iBlock = ',iBlock,'iNeutralBlockLast = ',iNeutralBlockLast

       xmin =  75000e3*Si2No_V(UnitX_) ! cometopause 75'000 km
       xmax = 100000e3*Si2No_V(UnitX_) ! cometopause max
       widthmin = 2.*xmin              ! flaring ratio 2 (Cravens 1989)
       widthmax = 2.*xmax              ! flaring ratio 2 (Cravens 1989)

       xMaxyz = -(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)/widthmax**2*xmax+xmax

       if(Xyz_DGB(x_,i,j,k,iBlock) > xMaxyz) then
          HeatCondFactor = 0.0
          HeatCondOut = 0.0
       else
          call derive_cell_diffusivity(iBlock, i, j, k, TeSI, nIon_I(1:nIonFluid), nElec, EtaSi)

          xMinyz = -(Xyz_DGB(y_,i,j,k,iBlock)**2+Xyz_DGB(z_,i,j,k,iBlock)**2)/widthmin**2*xmin+xmin
          if (Xyz_DGB(x_,i,j,k,iBlock) < xMinyz) then
             HeatCondFactor = 1.0
          else
             HeatCondFactor=(xMaxyz-Xyz_DGB(x_,i,j,k,iBlock))/(xMaxyz-xMinyz)
          end if
          HeatCondOut = TeSI/EtaSI*KappaCoeffSI*HeatCondFactor
       end if

       ! TestArray(1,i,j,k,iBlock) = nElec*sqrt(cBoltzmann*TeSi/cElectronMass)*cBoltzmann*TeSI/HeatCondOut
       ! TestArray(1,i,j,k,iBlock) = HeatCondOut
       ! TestArray(1,i,j,k,iBlock) = HeatCondFactor

    end if

    if(DoTestMe) then
       write(*,*)'user_material_properties:'
       write(*,*)'n_e    = ',nElec," [1/m^3]"
       write(*,*)'Te     = ',TeSI," [K]"
       if(present(CvOut)) write(*,*)'Cv     = ',CvOut,' [J/(K*m^3)]'
       if(present(HeatCondOut)) then
          write(*,*)'Eta    = ',EtaSI," [m^2/s]"
          write(*,*)'Kappa  = ',HeatCondOut," [W/(m*K)]"
          write(*,*)'Ffree  = ',nElec*sqrt(cBoltzmann*TeSi/cElectronMass)*cBoltzmann*TeSI," [W/m^2]"
       end if
       write(*,*)''
    end if
  end subroutine user_material_properties

  !========================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !----------------------------------------------------------------------

    !! Source terms are evaluated explicitly!
    !RETURN

    ! All ion momenta are implicit
    if(UseElectronPressure)then
       allocate(iVarPointImpl_I(5*nIonFluid + 1))
       iVarPointImpl_I(5*nIonFluid + 1) = Pe_
    else
       allocate(iVarPointImpl_I(5*nIonFluid))
    end if

    do iFluid = 1, nIonFluid
       iVarPointImpl_I(5*iFluid-4) = iRhoIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-3) = iRhoUxIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-2) = iRhoUyIon_I(iFluid)
       iVarPointImpl_I(5*iFluid-1) = iRhoUzIon_I(iFluid)
       iVarPointImpl_I(5*iFluid)   = iPIon_I(iFluid)
    end do

    IsPointImplMatrixSet = .false.
    !IsAsymmetric= .false.

  end subroutine user_init_point_implicit

  !========================================================================

  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional,&
       PlotVar_G, PlotVarBody, UsePlotVarBody,&
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance,    ONLY: State_VGB
    use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitU_, UnitT_, &
         ElectronCharge, ElectronPressureRatio, UnitTemperature_
    use ModVarIndexes, ONLY: P_, Pe_
    use ModConst,      ONLY: cBoltzmann
    use ModCurrent,    ONLY: get_current
    use ModMultiFluid, ONLY: MassIon_I
    use ModMain,       ONLY: Dt_BLK

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(-1:nI+2, -1:nJ+2, -1:nK+2)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    integer :: i, j, k, iIonFluid
    real :: nElec
    real, dimension(3)           :: Current_I, uIonMean_I
    real, dimension(nIonFluid)   :: nIon_I
    real, dimension(3,nIonFluid) :: uIon_I


    !--------------------------------------------------------------------------

    IsFound = .true.

    select case(NameVar)
    case('nn1')
       NameIdlUnit = '1/cm^3'
       NameTecUnit = '[1/cm^3]'
       !do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(Neu1Rho_,i,j,k,iBlock)/MassFluid_I(nFluid)*No2Si_V(UnitN_)
       end do; end do; end do

    case('unx1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(Neu1RhoUx_,i,j,k,iBlock)/ &
               State_VGB(Neu1Rho_,i,j,k,iBlock) *No2Si_V(UnitU_)
       end do; end do; end do

    case('uny1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(Neu1RhoUy_,i,j,k,iBlock)/ &
               State_VGB(Neu1Rho_,i,j,k,iBlock) *No2Si_V(UnitU_)
       end do; end do; end do

    case('unz1')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(Neu1RhoUz_,i,j,k,iBlock)/ &
               State_VGB(Neu1Rho_,i,j,k,iBlock) *No2Si_V(UnitU_)
       end do; end do; end do

    case('tn1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(Neu1P_,i,j,k,iBlock)* &
               MassFluid_I(nFluid)/State_VGB(Neu1Rho_,i,j,k,iBlock) * &
               No2SI_V(UnitTemperature_)
       end do; end do; end do

    case('te')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          if(UseElectronPressure)then
             PlotVar_G(i,j,k) = State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)/&
                  (cBoltzmann*nElec)
          else
             PlotVar_G(i,j,k) = State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)*ElectronPressureRatio/&
                  (1.+ElectronPressureRatio)/(cBoltzmann*nElec)
          end if
       end do; end do; end do

    case('ti1')
       NameIdlUnit = 'K'
       NameTecUnit = '[K]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = State_VGB(iPIon_I(1),i,j,k,iBlock)*NO2SI_V(UnitP_)/ &
               (cBoltzmann*State_VGB(iRhoIon_I(1),i,j,k,iBlock)/MassIon_I(1)*NO2SI_V(UnitN_))
       end do; end do; end do

    case('uex')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(1,1:nIonFluid) = State_VGB(iRhoUxIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(1) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(1) = uIonMean_I(1)+nIon_I(iIonFluid)* &
                  uIon_I(1,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(1)-Current_I(1)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uey')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(2,1:nIonFluid) = State_VGB(iRhoUyIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(2) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(2) = uIonMean_I(2)+nIon_I(iIonFluid)* &
                  uIon_I(2,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(2)-Current_I(2)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('uez')
       NameIdlUnit = 'km/s'
       NameTecUnit = '[km/s]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          call get_current(i,j,k,iBlock,Current_I)
          nIon_I(1:nIonFluid) = State_VGB(iRhoIon_I,i,j,k,iBlock)/MassIon_I*No2SI_V(UnitN_)
          uIon_I(3,1:nIonFluid) = State_VGB(iRhoUzIon_I,i,j,k,iBlock) / &
               State_VGB(iRhoIon_I,i,j,k,iBlock)*No2SI_V(UnitU_)
          nElec = sum(nIon_I(1:nIonFluid)*ChargeIon_I(1:nIonFluid))
          uIonMean_I(3) = 0.
          do iIonFluid=1,nIonFluid
             uIonMean_I(3) = uIonMean_I(3)+nIon_I(iIonFluid)* &
                  uIon_I(3,iIonFluid)/nElec*ChargeIon_I(iIonFluid)
          end do
          PlotVar_G(i,j,k) = (uIonMean_I(3)-Current_I(3)/(nElec*Si2No_V(UnitN_)*&
               ElectronCharge)*No2SI_V(UnitU_))/1E3
       end do; end do; end do

    case('dt')
       NameIdlUnit = 's'
       NameTecUnit = '[s]'
       PlotVar_G(:,:,:) = Dt_BLK(iBlock)*No2SI_V(UnitT_)

    case('ns')
       NameIdlUnit = ' '   
       NameTecUnit = '[ ]'
       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
          PlotVar_G(i,j,k) = ne20eV_GB(i,j,k,iBlock)
       end do; end do; end do
       ! case('testarray1')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(1,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray2')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(2,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray3')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(3,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('testarray4')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=1,nK; do j=1,nJ; do i=1,nI ! only idl
       !       !       do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = TestArray(4,i,j,k,iBlock)
       !    end do; end do; end do
       ! case('fluxlim')
       !    NameIdlUnit = ' '   
       !    NameTecUnit = '[ ]'
       !    do k=0,nK+1; do j=0,nJ+1; do i=0,nI+1
       !       PlotVar_G(i,j,k) = FluxLimited_GB(i,j,k,iBlock)
       !    end do; end do; end do
    case default
       IsFound = .false.
    end select

    UsePlotVarBody = .false.
    PlotVarBody    = 0.0

  end subroutine user_set_plot_var



  !========================================================================

  subroutine user_preset_conditions(i,j,k,iBlock)
    ! This is applied as initial conditions and in the upstream boundary for the semi-implicit heat conduction
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: SW_rho, SW_Ux, SW_Uy, SW_Uz, SW_p, LowDensityRatio, &
         ElectronPressureRatio, Io2No_V, UnitTemperature_, &
         SW_Bx, SW_By, SW_Bz, No2Si_V, UnitX_, Si2No_V, UnitU_, UnitN_
    use ModGeometry, ONLY:Xyz_DGB, r_BLK
    use ModVarIndexes, ONLY: MassFluid_I

    integer,intent(in) :: i, j, k, iBlock
    real:: RhoSw, RhoNeu1

    !--------------------------
    RhoSw = SW_rho*(1.0 - LowDensityRatio*(IonLast_ - IonFirst_))
    State_VGB(SwRho_,i,j,k,iBlock)     = RhoSw
    State_VGB(SwRhoUx_,i,j,k,iBlock)   = RhoSw*SW_Ux
    State_VGB(SwRhoUy_,i,j,k,iBlock)   = RhoSw*SW_Uy
    State_VGB(SwRhoUz_,i,j,k,iBlock)   = RhoSw*SW_Uz
    State_VGB(SwP_,i,j,k,iBlock)       = SW_p*(1.0-LowDensityRatio*(IonLast_-IonFirst_))

    State_VGB(H2OpRho_,i,j,k,iBlock)   = SW_rho*LowDensityRatio
    State_VGB(H2OpRhoUx_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Ux
    State_VGB(H2OpRhoUy_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Uy
    State_VGB(H2OpRhoUz_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Uz
    State_VGB(H2OpP_,i,j,k,iBlock)     = SW_p*LowDensityRatio !*MassIon_I(1)/MassIon_I(2)

    ! Neutral values are pre-set by Haser model
    RhoNeu1 = Qprod/( 4.*cPi*(R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_))**2*uHaser ) * &
         exp(-vHI*R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_)/uHaser)*Si2No_V(UnitN_) * &
         MassFluid_I(nFluid)
    State_VGB(Neu1Rho_,i,j,k,iBlock)   = RhoNeu1
    State_VGB(Neu1Ux_:Neu1Uz_,i,j,k,iBlock) = RhoNeu1*uHaser*Si2No_V(UnitU_) * &
         Xyz_DGB(:,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)
    State_VGB(Neu1P_,i,j,k,iBlock)     = State_VGB(Neu1Rho_,i,j,k,iBlock)/MassFluid_I(nFluid) &
         *50*Io2No_V(UnitTemperature_)

    State_VGB(Rho_,i,j,k,iBlock)       = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
    State_VGB(RhoUx_,i,j,k,iBlock)     = sum(State_VGB(iRhoUxIon_I,i,j,k,iBlock))
    State_VGB(RhoUy_,i,j,k,iBlock)     = sum(State_VGB(iRhoUyIon_I,i,j,k,iBlock))
    State_VGB(RhoUz_,i,j,k,iBlock)     = sum(State_VGB(iRhoUzIon_I,i,j,k,iBlock))

    State_VGB(Bx_,i,j,k,iBlock) = SW_Bx
    State_VGB(By_,i,j,k,iBlock) = SW_By
    State_VGB(Bz_,i,j,k,iBlock) = SW_Bz

    if(UseElectronPressure) then
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))
       State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
    else
       State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))* &
            (1.+ElectronPressureRatio)
    end if

  end subroutine user_preset_conditions

  !========================================================================

  subroutine user_set_ICs(iBlock)
    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest, Body1_
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: ElectronPressureRatio, No2Si_V, UnitRho_, &
         UnitRhoU_, UnitP_, CellState_VI
    use ModGeometry, ONLY: true_cell

    integer, intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true.
    integer :: i, j, k, iIonFluid
    ! !-------------------------------------------------------------------------
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_ICs', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
       if (.not.true_cell(i,j,k,iBlock)) then
          State_VGB(:,i,j,k,iBlock) = CellState_VI(:,Body1_)
          ! State_VGB(iPIon_I,i,j,k,iBlock) = BodyNDim_I*cBoltzmann*BodyTDim_I*1e6*SI2No_V(UnitP_)
          ! if (UseElectronPressure) then
          !    State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          !    State_VGB(Pe_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock)*&
          !         ElectronPressureRatio
          ! else
          !    State_VGB(P_,i,j,k,iBlock) = &
          !         sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
          ! end if
       else
          call user_preset_conditions(i,j,k,iBlock)
       end if

    end do; end do ; end do


    if(DoTestMe) then
       i=iTest ; j=jTest ; k=kTest
123    format (A13,ES25.16,A15)
       do iIonFluid=1,nIonFluid
          write(*,*)'Ion species #',iIonFluid,': ',NameFluid_I(iIonFluid+1)
          write(*,123)'Rho       = ',State_VGB(iRhoIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
          write(*,123)'rhoUx     = ',State_VGB(iRhoUxIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUy     = ',State_VGB(iRhoUyIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'rhoUz     = ',State_VGB(iRhoUzIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
               " [kg/(m^2*s)]"
          write(*,123)'Pi        = ',State_VGB(iPIon_I(iIonFluid),i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       end do
       write(*,*)''
       write(*,*)'Total:'
       write(*,123)'Rho       = ',State_VGB(Rho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
       write(*,123)'rhoUx     = ',State_VGB(RhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUy     = ',State_VGB(RhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       write(*,123)'rhoUz     = ',State_VGB(RhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
            " [kg/(m^2*s)]"
       if (UseElectronPressure) then
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(Pe_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"
       else
          write(*,123)'PiTot     = ',State_VGB(P_,i,j,k,iBlock)/(1.+ElectronPressureRatio)* &
               No2SI_V(UnitP_)," [Pa]"
          write(*,123)'Pe        = ',State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio/&
               (1.+ElectronPressureRatio)*No2SI_V(UnitP_)," [Pa]"
       end if

    end if

  end subroutine user_set_ICs

  !========================================================================

  subroutine user_set_cell_boundary(iBlock,iSide, TypeBc, IsFound)

    use ModAdvance,  ONLY: State_VGB
    use ModImplicit, ONLY: StateSemi_VGB, iTeImpl
    use ModSize,     ONLY: nI, MaxI, MinJ, MaxJ, MinK, MaxK
    use ModPhysics,  ONLY: Si2No_V, UnitTemperature_

    integer,          intent(in)  :: iBlock, iSide
    character(len=*),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    integer:: i, j, k
    real :: TeSi

    character(len=*), parameter :: NameSub = 'user_set_cell_boundary'
    !-------------------------------------------------------------------
    IsFound = .true.

    if(TypeBc == 'usersemi')then
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i = nI+1, MaxI

          call user_preset_conditions(i,j,k,iBlock)
          call user_material_properties(State_VGB(:,i,j,k,iBlock), &
               i, j, k, iBlock, TeOut=TeSi)
          StateSemi_VGB(iTeImpl,i,j,k,iBlock) = TeSi*Si2No_V(UnitTemperature_)


       end do; end do; end do

       RETURN
    elseif(TypeBc == 'usersemilinear')then
       RETURN
    end if


  end subroutine user_set_cell_boundary

  !============================================================================

  subroutine user_get_log_var(VarValue, TypeVar, Radius)

    use ModMain,       ONLY: Dt_BLK, BLKtest
    use ModPhysics,    ONLY: No2SI_V, UnitT_

    real, intent(out)             :: VarValue
    character (len=*), intent(in) :: TypeVar
    real, intent(in), optional    :: Radius

    character (len=*), parameter  :: NameSub = 'user_get_log_var'
    !-------------------------------------------------------------------------
    select case(TypeVar)
    case('dtpnt')
       VarValue = Dt_BLK(BLKtest)*No2SI_V(UnitT_)
    case default
       VarValue = -7777.0
    end select

  end subroutine user_get_log_var

  !============================================================================

  subroutine user_amr_criteria(iBlock, UserCriteria, TypeCriteria, IsFound)

    ! Set UserCriteria = 1.0 for refinement, 0.0 for coarsening.                                                                                                                             
    use BATL_lib,    ONLY: iNode_B, iTree_IA, Level_
    use ModAdvance,  ONLY: State_VGB, H2OpRho_
    use ModPhysics,  ONLY: No2SI_V, UnitN_
    use ModGeometry, ONLY: Xyz_DGB

    ! Variables required by this user subroutine

    integer, intent(in)          :: iBlock
    real, intent(out)            :: UserCriteria
    character (len=*),intent(in) :: TypeCriteria
    logical ,intent(inout)       :: IsFound

    integer:: i, j, k, nLevel, iNode
    !------------------------------------------------------------------
    IsFound = .true.

    UserCriteria = 0.0 ! only refinement

    iNode = iNode_B(iBlock)
    nLevel = iTree_IA(Level_,iNode)

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI

       if ((Xyz_DGB(x_,i,j,k,iBlock) < 0.1).and. &
            (Xyz_DGB(x_,i,j,k,iBlock) > -0.4).and. &
            (Xyz_DGB(y_,i,j,k,iBlock) <  0.05).and. &
            (Xyz_DGB(y_,i,j,k,iBlock) > -0.05).and. &
            (Xyz_DGB(z_,i,j,k,iBlock) <  0.2).and. &
                                ! (Xyz_DGB(z_,i,j,k,iBlock) > -0.1).and. &
                                ! (R_BLK(i,j,k,iBlock) > 0.01)) then
            (Xyz_DGB(z_,i,j,k,iBlock) > -0.1)) then


          if (State_VGB(H2OpRho_,i,j,k,iBlock)/MassIon_I(H2Op_)*No2SI_V(UnitN_)>1.0) then
             UserCriteria = 1.0
             RETURN
          end if
          !if (State_VGB(H2OpRhoUx_,i,j,k,iBlock)/State_VGB(H2OpRho_,i,j,k,iBlock)*No2SI_V(UnitU_)<-410e3) then
          !   UserCriteria = 1.0
          !   RETURN
          !end if
       end if
    end do; end do; end do

  end subroutine user_amr_criteria


  subroutine user_initial_perturbation
    ! This is applied to reset all ions after the neutral background
    ! is fully developed

    use ModMain, ONLY: Unused_B
    use ModAdvance,  ONLY: P_, Pe_, State_VGB
    use ModPhysics,  ONLY: SW_rho, SW_Ux, SW_Uy, SW_Uz, SW_p, LowDensityRatio, &
         ElectronPressureRatio, &
         SW_Bx, SW_By, SW_Bz
    use ModGeometry, ONLY: true_cell
    use ModEnergy,    ONLY: calc_energy_cell

    integer :: i, j, k, iBlock
    real:: RhoSw

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    ! -------------------------------------------------------------------------

    do iBlock = 1, nBlock;

       if(Unused_B(iBlock)) CYCLE
       do k = MinK, MaxK; do j = MinJ, MaxJ; do i=MinI, MaxI
          if(.not. true_cell(i,j,k,iBlock)) CYCLE

          RhoSw = SW_rho*(1.0 - LowDensityRatio*(IonLast_ - IonFirst_))
          State_VGB(SwRho_,i,j,k,iBlock)     = RhoSw
          State_VGB(SwRhoUx_,i,j,k,iBlock)   = RhoSw*SW_Ux
          State_VGB(SwRhoUy_,i,j,k,iBlock)   = RhoSw*SW_Uy
          State_VGB(SwRhoUz_,i,j,k,iBlock)   = RhoSw*SW_Uz
          State_VGB(SwP_,i,j,k,iBlock)       = SW_p*(1.0-LowDensityRatio*(IonLast_-IonFirst_))

          State_VGB(H2OpRho_,i,j,k,iBlock)   = SW_rho*LowDensityRatio
          State_VGB(H2OpRhoUx_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Ux
          State_VGB(H2OpRhoUy_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Uy
          State_VGB(H2OpRhoUz_,i,j,k,iBlock) = SW_rho*LowDensityRatio*SW_Uz
          State_VGB(H2OpP_,i,j,k,iBlock)     = SW_p*LowDensityRatio !*MassIon_I(1)/MassIon_I(2)

          State_VGB(Rho_,i,j,k,iBlock)       = sum(State_VGB(iRhoIon_I,i,j,k,iBlock))
          State_VGB(RhoUx_,i,j,k,iBlock)     = sum(State_VGB(iRhoUxIon_I,i,j,k,iBlock))
          State_VGB(RhoUy_,i,j,k,iBlock)     = sum(State_VGB(iRhoUyIon_I,i,j,k,iBlock))
          State_VGB(RhoUz_,i,j,k,iBlock)     = sum(State_VGB(iRhoUzIon_I,i,j,k,iBlock))

          State_VGB(Bx_,i,j,k,iBlock) = SW_Bx
          State_VGB(By_,i,j,k,iBlock) = SW_By
          State_VGB(Bz_,i,j,k,iBlock) = SW_Bz

          if(UseElectronPressure) then
             State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))
             State_VGB(Pe_,i,j,k,iBlock)     = State_VGB(P_,i,j,k,iBlock)*ElectronPressureRatio
          else
             State_VGB(P_,i,j,k,iBlock)      = sum(State_VGB(iPIon_I,i,j,k,iBlock))* &
                  (1.+ElectronPressureRatio)
          end if

       end do; end do; end do

       call calc_energy_cell(iBlock)

    end do

  end subroutine user_initial_perturbation

end module ModUser
