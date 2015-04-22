!#NOTPUBLIC  email:shyinsi@umich.edu  expires:12/31/2099
!This code is a copyright protected software (c) 2002- University of Michigan
!========================================================================

module ModUser

  use ModSize      
  use ModVarIndexes 
  use ModMultiFluid
  use ModUserEmpty,                               &
       IMPLEMENTED1  => user_read_inputs,         &
       IMPLEMENTED2  => user_calc_sources,        &
       IMPLEMENTED3  => user_set_face_boundary,   &
       IMPLEMENTED4  => user_init_point_implicit, &
       IMPLEMANTED5  => user_set_ICs,             &
       IMPLEMENTED6  => user_set_boundary_cells,  &
       IMPLEMENTED7  => user_update_states,       &
       IMPLEMENTED8  => user_init_session,        &
       IMPLEMENTED9  => user_set_plot_var
 
  use ModProcMH, ONLY: iProc
  use ModMain, ONLY: xTest, yTest, zTest
  use ModVarIndexes, ONLY: nVar
  use ModNumConst, ONLY: cPi
  use ModPhysics, ONLY: Io2No_V, Si2No_V, No2Si_V, &
       UnitRho_, UnitU_, UnitTemperature_, UnitT_, &
       UnitP_, UnitN_, UnitX_, Gamma, GammaMinus1      


  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 0.1
  character (len=*), parameter :: NameUserModule = &
       'Yinsi Shou, 2 neutral fluid & 6 dust Comet module, Apr 2014'

  integer, parameter, public :: nNeutral = 2 !! number of neutral species
  !! Neutral species names
  integer, parameter, public :: nDust = 6 !! number of dust species

  character (len=6), parameter, public :: NameNeutral_I(nNeutral) = &
       (/ 'H2O  ','H    '/)

  character (len=6), parameter, public :: NameDust_I(nDust) = &
       (/ 'Dust1 ','Dust2 ','Dust3 ','Dust4 ','Dust5 ','Dust6 '/)

  integer, parameter :: &
        H2ORho_= 1,   &
        H2ORhoUx_= 2, H2OUx_= 2,  &
        H2ORhoUy_= 3, H2OUy_= 3,  &
        H2ORhoUz_= 4, H2OUz_= 4, &
        H2OP_= 5,     &
        HRho_= 6,   &
        HRhoUx_= 7, HUx_= 7,  &
        HRhoUy_= 8, HUy_= 8,  &
        HRhoUz_= 9, HUz_= 9, &
        HP_= 10


  integer, parameter :: H2O_  =  1,H_  =  2
  integer, parameter :: Dust1_=1,Dust2_=2,Dust3_=3,Dust4_=4,Dust5_=5,Dust6_=6 
 

  real, dimension(nNeutral) :: uNeutr_I,DestructRate_I 
   
  integer :: iNeutralBlockLast = -1
  real :: Qprod, Tmin, rHelio,Tn,unr_km,ion_rate, bodyRadius, dust2gas
  integer, parameter :: nNeuFluid = nNeutral  
  integer, parameter :: iRhoNeu_I(nNeuFluid)   = iRho_I(1:nNeutral)
  integer, parameter :: iRhoUxNeu_I(nNeuFluid) = iRhoUx_I(1:nNeutral)
  integer, parameter :: iRhoUyNeu_I(nNeuFluid) = iRhoUy_I(1:nNeutral)
  integer, parameter :: iRhoUzNeu_I(nNeuFluid) = iRhoUz_I(1:nNeutral)
  integer, parameter :: iPNeu_I(nNeuFluid)     = iP_I(1:nNeutral)
  real ::  MassNeu_I(nNeuFluid)   !in amu

  integer, parameter :: nDustFluid = nDust  
  integer, parameter :: iRhoDust_I(nDustFluid)   = iRho_I(nNeutral+1:nFluid)
  integer, parameter :: iRhoUxDust_I(nDustFluid) = iRhoUx_I(nNeutral+1:nFluid)
  integer, parameter :: iRhoUyDust_I(nDustFluid) = iRhoUy_I(nNeutral+1:nFluid)
  integer, parameter :: iRhoUzDust_I(nDustFluid) = iRhoUz_I(nNeutral+1:nFluid)
  integer, parameter :: iPDust_I(nDustFluid)     = iP_I(nNeutral+1:nFluid)
  real ::  MassDust_I(nDustFluid)   !in amu
  real, parameter :: dustradius_I(nDustFluid)= (/1e-7, 1e-6,1e-5,1e-4,1e-3,1e-2/) 
  real, parameter ::  dustNumDistr_I(nDustFluid) = dustradius_I**(-4)  ! dust number ~ radius^(-4)
  real, parameter ::  dustMassDistr_I(nDustFluid) = dustNumDistr_I*dustradius_I**(3)
  real ::  dustDensity = 1e3 ! kg/m^3
  ! assmue dust particles of all sizes have a constant density 
  
  character (len=100) :: NameGravityFile

  real:: Gravity_DCB(3,MinI:MaxI,MinJ:MaxJ,MinK:MaxK,MaxBlock) = 0.
  integer :: addGravity
  logical :: UseGravity(MaxBlock)= .false.
  ! for real  shape
  character (len=100) :: NameShapeFile

  integer:: nTriangle
  real, allocatable:: XyzTriangle_DII(:,:,:), Normal_DI(:,:)
  real :: rMinShape = 0.0, rMaxShape = 0.0

  ! Position of the sun at the start time
  real :: LatSun=43.5, LonSun=0.0

  ! Rotation of the comet (changes the direction of the Sun)
  real:: RotationCometHour = 12.0
  
  ! Angular velocity
  real:: OmegaCometSi

  ! Maximum change in longitude before updating the boundary conditions
  real:: AngleUpdateDeg = 10.0

  ! The time between updates
  real:: DtUpdateSi
  integer:: DnUpdate = 0, nStepStart


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
  real :: SlopeProduction, bProduction, SlopeTemp, bTemp, cos75

  ! Coefficients to calculate uNormal and temperature from TempCometLocal
  real :: TempToUnormal
  real :: TempToPressure

  ! Last step and time the inner boundary values were saved for each block
  integer:: nStepSave_B(MaxBlock) = -100
  real :: TimeSimulationSave_B(MaxBlock) = -1e30

contains


  !========================================================================
  subroutine user_read_inputs

    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    character (len=100) :: NameCommand

    !-----------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input Comet starts'
    endif

    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case("#COMET")
          call read_var('Qprod' , Qprod)            !! Neutral gas production rate [1/s]
          call read_var('rHelio', rHelio)           !! Heliocentric distance [AU]
          call read_var('Tmin' , Tmin)              !! Minimum ion temperature (enforced in update states)
          call read_var('Tn', Tn)                   !! Neutral temperature at the inner boundary (K)
          call read_var('unr',unr_km)               !! the neutral velocity at the inner boundary (km/s)
          call read_var('ionization_rate',ion_rate) !!
          call read_var('rbody',bodyRadius)         !! in km
          call read_var('dust2gas', dust2gas)       !! dust2gas mass ratio
          call read_var('dustdensity', dustDensity) !! dust density in kg/m^3
          call read_var('addGravity', addGravity)
          if (addGravity /=0) addGravity=1 
         bodyRadius = bodyRadius *1e3   ! in m, SI unit
         uNeutr_I(H2O_) = unr_km*1e3  ! [m/s] 
         DestructRate_I(H2O_)=ion_rate
!=========================================
       case("#GRAVITYFILE")
          call read_var('NameGravityFile' ,NameGravityFile)        
       case("#SHAPEFILE")
          call read_var('NameShapeFile' ,NameShapeFile)
       case("#SUNDIRECTION")
          call read_var('LatSun', LatSun)
          call read_var('LonSun', LonSun)

          ! If the Sun's position is changed, recalculate illumination
          nStepSave_B = -100
          TimeSimulationSave_B = -1e30

       case("#COMETSTATE")
          call read_var('ProductionRateMinSi', ProductionRateMinSi)
          call read_var('ProductionRateMaxSi', ProductionRateMaxSi)
          call read_var('SolarAngleMaxDim',    SolarAngleMaxDim)
          call read_var('TempCometMinDim',     TempCometMinDim)
          call read_var('TempCometMaxDim',     TempCometMaxDim)
          call read_var('TempComet75',         TempComet75Dim)
       case("#COMETROTATION")
          call read_var('RotationCometHour', RotationCometHour)
          call read_var('AngleUpdateDeg',    AngleUpdateDeg)
          call read_var('DnUpdate',          DnUpdate)
          nStepStart = n_step

!=========================================
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input Comet ends'
          endif
          EXIT
       case default
          if(iProc==0) then
             call write_myname; write(*,*) &
                  'ERROR: Invalid user defined #COMMAND in user_read_inputs. '
             write(*,*) '--Check user_read_inputs for errors'
             write(*,*) '--Check to make sure a #USERINPUTEND command was used'
             write(*,*) '  *Unrecognized command was: '//NameCommand
             call stop_mpi('ERROR: Correct PARAM.in or user_read_inputs!')
          end if
       end select
    end do
!===================================real comet
  UseUserInitSession = .true.
  UseExtraBoundary   = .true.
!===================================


  end subroutine user_read_inputs

!==============================================================================
 subroutine user_init_session

    ! Read shape file and convert units

    use ModMain, ONLY: Time_Simulation
    use ModNumConst, ONLY: cTwoPi, cDegToRad
    use ModCoordTransform, ONLY: dir_to_xyz
    use ModConst, ONLY: cBoltzmann, cAtomicMass
    use ModVarIndexes, ONLY: MassFluid_I
    use ModBlockData, ONLY: MaxBlockData
    !------------------------------------------------------------------------
    ! We need to have unit conversions before reading the shape file 
    ! which contains everything in SI units
    call read_shape_file


    ! Production rate = number density * velocity, the unit is in [m^-2 s^-1]
    ! But, applying 1 / Si2No_V(UnitX_)**2 / Io2No_V(UnitT_) is not correct
    ! because Si2No_V(UnitN_) /= 1/Si2No_V(UnitX_)**3
    ProductionRateMax = &
         ProductionRateMaxSi * Si2No_V(UnitN_) * Si2No_V(UnitU_)
    ProductionRateMin = &
         ProductionRateMinSi * Si2No_V(UnitN_) * Si2No_V(UnitU_)

    SolarAngleMax= SolarAngleMaxDim * cDegToRad
    TempCometMin = TempCometMinDim * Io2No_V(UnitTemperature_)
    TempCometMax = TempCometMaxDim * Io2No_V(UnitTemperature_)
    TempComet75  = TempComet75Dim  * Io2No_V(UnitTemperature_)

    ! From Gamma. Toth's derivations
    ! uNormal = sqrt( kT / m ), so TempToUnormal = sqrt( k/ m )
    ! and also unit conversions of temperature to SI, and velocity from SI
    TempToUnormal = sqrt(cBoltzmann/(MassFluid_I(1)*cAtomicMass) * &
         No2Si_V(UnitTemperature_))*Si2No_V(UnitU_)

    ! From Gamma. Toth's derivations
    ! T' = T/Gamma so and p = n*T' = rho*T/(Gamma*m) so TempToPressure = 1/(Gamma*m)
    TempToPressure = 1/(Gamma*MassFluid_I(1))

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
    OmegaCometSi = cTwoPi / (RotationCometHour*3600)

    ! Frequency of boundary condition updates
    DtUpdateSi = AngleUpdateDeg*cDegToRad / abs(OmegaCometSi)

    ! Maximum amount of data to be stored in ModBlockData
    ! This is for the inner boundary conditions. 
    ! In practice this is a rather generous overestimate.
    MaxBlockData = nVar*(nI+1)*(nJ+1)*(nK+1)

    if(iProc==0)then
       write(*,*) 'ProductionRateMaxSi, ProductionRateMax =', &
            ProductionRateMaxSi, ProductionRateMax
       write(*,*) 'ProductionRateMinSi, ProductionRateMin =', &
            ProductionRateMinSi, ProductionRateMin
       write(*,*) 'SolarAngleMaxDim, SolarAngleMax =', &
            SolarAngleMaxDim, SolarAngleMax
       write(*,*) 'TempCometMinDim, TempCometMin   =', &
            TempCometMinDim, TempCometMin
       write(*,*) 'TempCometMaxDim, TempCometMax   =', &
            TempCometMaxDim, TempCometMax
       write(*,*) 'TempComet75Dim,  TempComet75    =', &
            TempComet75Dim,  TempComet75
       write(*,*) 'TempToUnormal, TempToPressure   =', &
            TempToUnormal, TempToPressure
       write(*,*) 'MassFluid_I,                    =', MassFluid_I
       write(*,*) 'TempToUn*sqrt(TempCometMax)     =', &
            TempToUnormal*sqrt(TempCometMax)
       write(*,*) 'SlopeProduction, bProduction    =', &
            SlopeProduction, bProduction
       write(*,*) 'SlopeTemp, bTemp                =', &
            SlopeTemp/Io2No_V(UnitTemperature_),bTemp/Io2No_V(UnitTemperature_)
       write(*,*)'LatSun, LonSun initial           =', &
            LatSun, LonSun
       write(*,*)'RotationCometHour, OmegaCometSi  =', &
            RotationCometHour, OmegaCometSi
       write(*,*)'AngleUpdateDeg, DtUpdateSi       =', &
            AngleUpdateDeg, DtUpdateSi
    end if

  end subroutine user_init_session
!===========================================================================
subroutine user_set_boundary_cells(iBlock)

    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, r_BLK
    integer, intent(in):: iBlock

    integer:: i, j, k
    real:: XyzInside_D(3)
    !------------------------------------------------------------------------
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

  end subroutine user_set_boundary_cells
!=============================================================================
subroutine read_gravity_file

    use ModPlotFile, ONLY: read_plot_file
    use BATL_lib, ONLY: find_grid_block, nNodeUsed 
    use ModGeometry, ONLY: Xyz_DGB  
    logical:: DoReadGravityFile = .true.
    integer:: iPoint, nPoint, i, j, k, iBlock
    integer:: iCell, nCell, iProcFound
    integer:: iCell_D(3)
    real :: Dist_D(3)
    real, allocatable:: Gravity_DI(:,:), Xyz_DI(:,:)
    character(len=100):: String1, String2

    character(len=*), parameter:: NameSub = 'read_gravity_file'
  !-----------------------------------------------------------------------
   if(.not.DoReadGravityFile) RETURN
    DoReadGravityFile = .false.
  
    if(iProc==0)write(*,*) NameSub,' reading gravity file ',trim(NameGravityFile)

    call read_plot_file(NameGravityFile, n1Out = nCell)

    if(nCell /= nI*nj*nK*nNodeUsed) then
    write (*,*) 'ncell=', nCell
    write(*,*) 'ni,nj,nk,nNodeused:', nI,nj,nK,nNodeUsed
    call stop_mpi(NameSub//' Wrong grid ')
    end if

    allocate(Gravity_DI(3,nCell), Xyz_DI(3,nCell))
    call read_plot_file(NameGravityFile, CoordOut_DI=Xyz_DI, VarOut_VI=Gravity_DI)

  !Loop through the elements and put into a block based data:

   do iCell = 1, nCell
        call find_grid_block(Xyz_DI(:,iCell), iProcFound, iBlock, iCell_D, Dist_D)
        if(iProcFound /= iProc) CYCLE
        where(Dist_D > 0.5) iCell_D = iCell_D + 1
        i = iCell_D(1); j = iCell_D(2); k = iCell_D(3)
        Gravity_DCB(:,i,j,k,iBlock) = Gravity_DI(:,iCell)
    end do



   end subroutine read_gravity_file

!=============================================================================
subroutine read_shape_file

    use ModIoUnit, ONLY: UnitTmp_
    use ModCoordTransform, ONLY: cross_product
    use ModRandomNumber, ONLY: random_real
    
    logical:: DoReadShapeFile = .true.

    integer:: nPoint, i, j, iPoint, iTriangle, iPoint1, iPoint2, iPoint3
    integer:: iSeed = 7

    real, allocatable:: Xyz_DI(:,:)
    real :: RandomNumber

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
  subroutine user_calc_sources(iBlock)

    use ModMain,       ONLY: nI, nJ, nK, iTest, jTest, kTest, &
         BlkTest, PROCtest, iteration_number, Dt_BLK, n_Step
    use ModAdvance,    ONLY: State_VGB, Source_VC
    use ModConst,      ONLY: cBoltzmann, cElectronMass, cElectronCharge, cProtonMass
    use ModGeometry,   ONLY: Rmin_BLK, r_BLK, Xyz_DGB
    use ModProcMH,     ONLY: iProc
    use ModPhysics
    use ModPointImplicit, ONLY: UsePointImplicit_B, UsePointImplicit, IsPointImplSource

    integer, intent(in) :: iBlock

   
    real, dimension(1:nI,1:nJ,1:nK) :: uNeu2_C
    real, dimension(3,1:nNeuFluid,1:nI,1:nJ,1:nK) :: SRhoUxNeuTerm_IIC, SRhoUyNeuTerm_IIC, SRhoUzNeuTerm_IIC
    real, dimension(4,1:nNeuFluid,1:nI,1:nJ,1:nK) :: SPNeuTerm_IIC

    real, dimension(1:nNeuFluid,1:nI,1:nJ,1:nK) :: SRhoNeu_IC, &
         SRhoUxNeu_IC, SRhoUyNeu_IC, SRhoUzNeu_IC, SPNeu_IC
    real, dimension(1:nNeuFluid,1:nI,1:nJ,1:nK) :: NnNeutral_IC,UnxNeutral_IC, &
         UnyNeutral_IC,UnzNeutral_IC,Tn_IC
   
    real, dimension(1:nDustFluid,1:nI,1:nJ,1:nK) :: SRhoDust_IC, &
         SRhoUxDust_IC, SRhoUyDust_IC, SRhoUzDust_IC
    real, dimension(1:nDustFluid,1:nI,1:nJ,1:nK) :: NnDust_IC,UxDust_IC, &
         UyDust_IC,UzDust_IC
    real, dimension(1:nDustFluid,1:nI,1:nJ,1:nK) ::  &
         SRhoUxDustGravity_IC, SRhoUyDustGravity_IC, SRhoUzDustGravity_IC   

    real, dimension(1:nDustFluid) :: relMach_I,Trec_I,Tdust_I,Ch_I,Udiff_I
   

    real  :: vH2O_Destr, vH_prod
    real  :: temp1, temp2, excess_energy 
    logical :: DoTest, DoTestMe=.true.
    integer :: i,j,k,iNeutral,jNeutral,iTerm,iDim,iDust
    real  :: cDrag=2.0
    real  :: thermalU
    real :: cross_section  !between H2O and H




    !----------------------------------------------------------------------

    ! Do not evaluate any source terms explicitly when running pointimplicit
    ! if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    ! Evaluate source terms explicitly even when running pointimplicit
    ! if(UsePointImplicit .and. IsPointImplSource) RETURN

    !! Limit region for evaluation for source term evaluation
    !! if(RMin_BLK(iBlock) > 2.) RETURN
   
    if(iBlock == BlkTest) then
       call set_oktest('user_calc_sources',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

  
    MassNeu_I = MassFluid_I(1:nNeuFluid)
     
    MassDust_I = 4/3*cPi*dustradius_I**3*dustDensity/cProtonMass
   
    vH2O_Destr = 1e-6/9 !s^-1 in SI
    vH_prod = 1e-5/9 !s^-1 in SI
    excess_energy = 1.6e-19*1.9  ! 1.9 eV in SI    
    Tdust_I=(/280,280,280,280,280,280/)
      
     SRhoNeu_IC = 0.0
     SRhoUxNeu_IC=0.0 
     SRhoUyNeu_IC = 0.0
     SRhoUzNeu_IC=0.0 
     SPNeu_IC= 0.0

     SRhoDust_IC = 0.0
     SRhoUxDust_IC=0.0 
     SRhoUyDust_IC = 0.0
     SRhoUzDust_IC=0.0 
     
     SRhoUxNeuTerm_IIC = 0.0 
     SRhoUyNeuTerm_IIC = 0.0
     SRhoUzNeuTerm_IIC = 0.0
     SPNeuTerm_IIC = 0.0 
    
    if (addGravity==1 .and. .not.UseGravity(iBlock)) then 
    call read_gravity_file
    UseGravity(iBlock)= .true.
    end if

    do k=1,nK; do j=1,nJ; do i=1,nI
     ! Calculate normalized neutral density 
       NnNeutral_IC(1:nNeuFluid,i,j,k) = State_VGB(iRhoNeu_I,i,j,k,iBlock)/MassNeu_I*No2Si_V(UnitN_)
       Tn_IC(1:nNeuFluid,i,j,k)= State_VGB(iPNeu_I,i,j,k,iBlock)*No2Si_V(UnitP_)&
                                 /cBoltzmann/NnNeutral_IC(1:nNeuFluid,i,j,k)

     ! Calculate dust density
      NnDust_IC(1:nDustFluid,i,j,k) = State_VGB(iRhoDust_I,i,j,k,iBlock)/MassDust_I*No2Si_V(UnitN_)
    end do; end do; end do


     !!Calculate normalized neutral velocity 
    UnxNeutral_IC(1:nNeuFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUxNeu_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoNeu_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)
    UnyNeutral_IC(1:nNeuFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUyNeu_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoNeu_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)
    UnzNeutral_IC(1:nNeuFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUzNeu_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoNeu_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)

      !!Calculate dust velocity 
    UxDust_IC(1:nDustFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUxDust_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoDust_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)
    UyDust_IC(1:nDustFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUyDust_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoDust_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)
    UzDust_IC(1:nDustFluid,1:nI,1:nJ,1:nK)=State_VGB(iRhoUzDust_I,1:nI,1:nJ,1:nK, iBlock)/ &
         State_VGB(iRhoDust_I,1:nI,1:nJ,1:nK,iBlock)*No2Si_V(UnitU_)

     uNeu2_C(1:nI,1:nJ,1:nK) =  &
          (UnxNeutral_IC(H2O_,1:nI,1:nJ,1:nK)-UnxNeutral_IC(H_,1:nI,1:nJ,1:nK))**2.+&
          (UnyNeutral_IC(H2O_,1:nI,1:nJ,1:nK)-UnyNeutral_IC(H_,1:nI,1:nJ,1:nK))**2.+& 
          (UnzNeutral_IC(H2O_,1:nI,1:nJ,1:nK)-UnzNeutral_IC(H_,1:nI,1:nJ,1:nK))**2.


    do k=1,nK; do j=1,nJ; do i=1,nI
     
   thermalU = sqrt(2*cBoltzmann*Tn_IC(H2O_,i,j,k)/MassNeu_I(H2O_)/cProtonMass)
  ! Udiff_I=sqrt((UnxNeutral_IC(H2O_,i,j,k)-UxDust_IC(Dust1_:Dust6_,i,j,k))**2&
  !        +(UnyNeutral_IC(H2O_,i,j,k)-UyDust_IC(Dust1_:Dust6_,i,j,k))**2&
  !    +(UnzNeutral_IC(H2O_,i,j,k)-UzDust_IC(Dust1_:Dust6_,i,j,k))**2)
   Udiff_I=sqrt(UnxNeutral_IC(H2O_,i,j,k)**2.+UnyNeutral_IC(H2O_,i,j,k)**2+&
         UnzNeutral_IC(H2O_,i,j,k)**2.) - &
           sqrt(UxDust_IC(Dust1_:Dust6_,i,j,k)**2.+ &
         UyDust_IC(Dust1_:Dust6_,i,j,k)**2.+UzDust_IC(Dust1_:Dust6_,i,j,k)**2.)  

   relMach_I = Udiff_I/thermalU
       !! Zeroth moment
         
   SRhoNeu_IC(H2O_,i,j,k)=-vH2O_Destr*State_VGB(H2ORho_,i,j,k,iBlock)/Si2No_V(UnitT_)
 
   SRhoNeu_IC(H_,i,j,k) = vH_prod*State_VGB(H2ORho_,i,j,k,iBlock)/MassNeu_I(H2O_)*&
                          MassNeu_I(H_)/Si2No_V(UnitT_)
       !dusts have no src in continuity eqn. The srcs are zeroed already.     
    

       !! First moment, x component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
    temp1 = MassNeu_I(H_)*MassNeu_I(H2O_)/(MassNeu_I(H_)+MassNeu_I(H2O_))*NnNeutral_IC(H_,i,j,k)*&
              NnNeutral_IC(H2O_,i,j,k)*cross_section
    temp2 = (Tn_IC(H_,i,j,k)/MassNeu_I(H_)+Tn_IC(H2O_,i,j,k)/MassNeu_I(H2O_))/cProtonMass*8.*&
              cBoltzmann/cPi
    temp2 = sqrt(temp2+ uNeu2_C(i,j,k))      

    SRhoUxNeuTerm_IIC(2,H2O_,i,j,k)=temp1*temp2*(UnxNeutral_IC(H_,i,j,k)-UnxNeutral_IC(H2O_,i,j,k))&
    *Si2No_V(UnitN_)*Si2No_V(UnitU_)/Si2No_V(UnitT_)

    SRhoUxNeuTerm_IIC(1,H2O_,i,j,k)= SRhoNeu_IC(H2O_,i,j,k)*UnxNeutral_IC(H2O_,i,j,k)*&
                                         Si2No_V(UnitU_)


    SRhoUxDust_IC(Dust1_:Dust6_,i,j,k)=0.5*cDrag*cPi*(dustradius_I(1:6)**2)*&
        (UnxNeutral_IC(H2O_,i,j,k)-UxDust_IC(Dust1_:Dust6_,i,j,k))*abs(Udiff_I)* &
        NnNeutral_IC(H2O_,i,j,k)*NnDust_IC(Dust1_:Dust6_,i,j,k)*MassNeu_I(H2O_)&
        *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)  

    SRhoUxNeuTerm_IIC(3,H2O_,i,j,k)=-sum(SRhoUxDust_IC(Dust1_:Dust6_,i,j,k))
   
    SRhoUxDustGravity_IC(:,i,j,k) = MassDust_I*NnDust_IC(:,i,j,k)*Gravity_DCB(1,i,j,k,iBlock) &
                                    *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)
         
  
    SRhoUxNeuTerm_IIC(2,H_,i,j,k)=-SRhoUxNeuTerm_IIC(2,H2O_,i,j,k)
    SRhoUxNeuTerm_IIC(1,H_,i,j,k)= SRhoNeu_IC(H_,i,j,k)*UnxNeutral_IC(H_,i,j,k)*Si2No_V(UnitU_)
  

       !! First moment, y component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
     SRhoUyNeuTerm_IIC(2,H2O_,i,j,k)=temp1*temp2*(UnyNeutral_IC(H_,i,j,k)-UnyNeutral_IC(H2O_,i,j,k))&
    *Si2No_V(UnitN_)*Si2No_V(UnitU_)/Si2No_V(UnitT_)

     SRhoUyNeuTerm_IIC(1,H2O_,i,j,k)= SRhoNeu_IC(H2O_,i,j,k)*UnyNeutral_IC(H2O_,i,j,k)*&
                                         Si2No_V(UnitU_)


        SRhoUyDust_IC(Dust1_:Dust6_,i,j,k)=0.5*cDrag*cPi*(dustradius_I(1:6)**2)*&
        (UnyNeutral_IC(H2O_,i,j,k)-UyDust_IC(Dust1_:Dust6_,i,j,k))*abs(Udiff_I)* &
        NnNeutral_IC(H2O_,i,j,k)*NnDust_IC(Dust1_:Dust6_,i,j,k)*MassNeu_I(H2O_)&
        *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)  

    SRhoUyNeuTerm_IIC(3,H2O_,i,j,k)=-sum(SRhoUyDust_IC(Dust1_:Dust6_,i,j,k))
     
    SRhoUyDustGravity_IC(:,i,j,k) = MassDust_I*NnDust_IC(:,i,j,k)*Gravity_DCB(2,i,j,k,iBlock) &
                                    *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)
    

    SRhoUyNeuTerm_IIC(2,H_,i,j,k)=-SRhoUyNeuTerm_IIC(2,H2O_,i,j,k)
    SRhoUyNeuTerm_IIC(1,H_,i,j,k)= SRhoNeu_IC(H_,i,j,k)*UnyNeutral_IC(H_,i,j,k)*Si2No_V(UnitU_)
  
            

 !! First moment, z component
       !! d(rho_s*u_s)/dt = rho_s*du_s/dt + u_s*drho_s/dt combined from zeroth and first moment by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"
   SRhoUzNeuTerm_IIC(2,H2O_,i,j,k)=temp1*temp2*(UnzNeutral_IC(H_,i,j,k)-UnzNeutral_IC(H2O_,i,j,k))&
    *Si2No_V(UnitN_)*Si2No_V(UnitU_)/Si2No_V(UnitT_)

   SRhoUzNeuTerm_IIC(1,H2O_,i,j,k)= SRhoNeu_IC(H2O_,i,j,k)*UnzNeutral_IC(H2O_,i,j,k)*&
                                         Si2No_V(UnitU_)


        SRhoUzDust_IC(Dust1_:Dust6_,i,j,k)=0.5*cDrag*cPi*(dustradius_I(1:6)**2)*&
       (UnzNeutral_IC(H2O_,i,j,k)-UzDust_IC(Dust1_:Dust6_,i,j,k))*abs(Udiff_I)* &
        NnNeutral_IC(H2O_,i,j,k)*NnDust_IC(Dust1_:Dust6_,i,j,k)*MassNeu_I(H2O_)&
        *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)  

   SRhoUzNeuTerm_IIC(3,H2O_,i,j,k)=-sum(SRhoUzDust_IC(Dust1_:Dust6_,i,j,k))

   SRhoUzDustGravity_IC(:,i,j,k) = MassDust_I*NnDust_IC(:,i,j,k)*Gravity_DCB(3,i,j,k,iBlock) &
                                    *SI2No_V(UnitN_)*SI2No_V(UnitU_)/SI2No_V(UnitT_)


   SRhoUzNeuTerm_IIC(2,H_,i,j,k)=-SRhoUzNeuTerm_IIC(2,H2O_,i,j,k)
   SRhoUzNeuTerm_IIC(1,H_,i,j,k)= SRhoNeu_IC(H_,i,j,k)*UnzNeutral_IC(H_,i,j,k)*Si2No_V(UnitU_)
      

       !! Second moment
       !! Sources separated into the terms by Tamas' "Transport Equations for Multifluid Magnetized Plasmas"      
    SPNeuTerm_IIC(1,H2O_,i,j,k) = -vH2O_Destr*State_VGB(P_,i,j,k,iBlock)
    
 !   SPNeuTerm_IIC(2,H2O_,i,j,k) =GammaMinus1*vH_prod*excess_energy*NnNeutral_IC(H2O_,i,j,k)* &
 !                                MassNeu_I(H_)/(MassNeu_I(H_)+MassNeu_I(H2O_))&
 !                                *SI2No_V(UnitP_)/SI2No_V(UnitT_)  
     
   
   temp2=temp2*NnNeutral_IC(H_,i,j,k)*NnNeutral_IC(H2O_,i,j,k)*MassNeu_I(H2O_)&
          *cProtonMass*cross_section

   SPNeuTerm_IIC(3,H2O_,i,j,k) = 2.*cBoltzmann * (Tn_IC(H_,i,j,k)-Tn_IC(H2O_,i,j,k))/(MassNeu_I(H_)*cProtonMass)&
                            + 2./3.*uNeu2_C(i,j,k)

   SPNeuTerm_IIC(3,H2O_,i,j,k) = SPNeuTerm_IIC(3,H2O_,i,j,k)*temp2*MassNeu_I(H_)/&
            (MassNeu_I(H_)+MassNeu_I(H2O_))*SI2No_V(UnitP_)/SI2No_V(UnitT_) 

  !  where(relMach_I<1e-3) relMach_I=1e-3
  !  Trec_I = Tn_IC(H2O_,i,j,k)/gp1*(2*Gamma+2*GammaMinus1*relMach_I**2-GammaMinus1/&
  !         (0.5+relMach_I**2+relMach_I*cPi**(-0.5)*exp(-relMach_I**2)/erf(relMach_I)))
  !  Ch_I = gp1/GammaMinus1*cBoltzmann/(8*cProtonMass*MassNeu_I(H2O_))/relMach_I**2*&
  !         (relMach_I*cPi**(-0.5)*exp(-relMach_I**2)+(0.5+relMach_I**2)*erf(relMach_I))
  !  SPNeuTerm_IIC(4,H2O_,i,j,k)= GammaMinus1*4*cPi*State_VGB(Rho_,i,j,k,iBlock)* &
  !             sum(NnDust_IC(1:6,i,j,k)*dustradius_I**2*relMach_I*thermalU*(Tdust_I-Trec_I)*Ch_I)*&
  !             No2SI_V(UnitRho_)*SI2No_V(UnitP_)




    SPNeuTerm_IIC(1,H_,i,j,k) = GammaMinus1*vH_prod*excess_energy* NnNeutral_IC(H2O_,i,j,k)* &
                                MassNeu_I(H2O_)/(MassNeu_I(H_)+MassNeu_I(H2O_))&
                                 *SI2No_V(UnitP_)/SI2No_V(UnitT_)
      
    SPNeuTerm_IIC(2,H_,i,j,k) = 2.*cBoltzmann * (Tn_IC(H2O_,i,j,k)-Tn_IC(H_,i,j,k))/(MassNeu_I(H2O_)*cProtonMass)&
                            + 2./3.*uNeu2_C(i,j,k)

    SPNeuTerm_IIC(2,H_,i,j,k) = SPNeuTerm_IIC(2,H_,i,j,k)*temp2*MassNeu_I(H2O_)/(MassNeu_I(H_)+MassNeu_I(H2O_))*SI2No_V(UnitP_)/SI2No_V(UnitT_)

       !! sum up individual terms
    !   do iTerm=1,2
    !      SRho_IC(1:nIonFluid,i,j,k) = SRho_IC(1:nIonFluid,i,j,k)+SRhoTerm_IIC(iTerm,1:nIonFluid,i,j,k)
    !   end do
      
        do iTerm=1,3
           SRhoUxNeu_IC(1:nNeuFluid,i,j,k) = SRhoUxNeu_IC(1:nNeuFluid,i,j,k)+&
                                             SRhoUxNeuTerm_IIC(iTerm,1:nNeuFluid,i,j,k)
           SRhoUyNeu_IC(1:nNeuFluid,i,j,k) = SRhoUyNeu_IC(1:nNeuFluid,i,j,k)+& 
                                             SRhoUyNeuTerm_IIC(iTerm,1:nNeuFluid,i,j,k)
           SRhoUzNeu_IC(1:nNeuFluid,i,j,k) = SRhoUzNeu_IC(1:nNeuFluid,i,j,k)+&
                                             SRhoUzNeuTerm_IIC(iTerm,1:nNeuFluid,i,j,k)
         
        end do
       
       do iTerm=1,3
          SPNeu_IC(1:nNeuFluid,i,j,k) = SPNeu_IC(1:nNeuFluid,i,j,k)+&
                                        SPNeuTerm_IIC(iTerm,1:nNeuFluid,i,j,k)
         
       end do
      
      

          
  
      Source_VC(iRhoNeu_I   ,i,j,k) = SRhoNeu_IC(1:nNeuFluid,i,j,k) +  Source_VC(iRhoNeu_I   ,i,j,k)
      Source_VC(iRhoUxNeu_I ,i,j,k) = SRhoUxNeu_IC(1:nNeuFluid,i,j,k)+ Source_VC(iRhoUxNeu_I ,i,j,k)
      Source_VC(iRhoUyNeu_I ,i,j,k) = SRhoUyNeu_IC(1:nNeuFluid,i,j,k)+ Source_VC(iRhoUyNeu_I ,i,j,k)
      Source_VC(iRhoUzNeu_I ,i,j,k) = SRhoUzNeu_IC(1:nNeuFluid,i,j,k)+ Source_VC(iRhoUzNeu_I ,i,j,k)
      Source_VC(iPNeu_I,i,j,k) = SPNeu_IC(1:nNeuFluid,i,j,k) + Source_VC(iPNeu_I,i,j,k)

      
      Source_VC(iRhoDust_I   ,i,j,k) = SRhoDust_IC(1:nDustFluid,i,j,k)+ Source_VC(iRhoDust_I,i,j,k)
      Source_VC(iRhoUxDust_I ,i,j,k) = SRhoUxDust_IC(1:nDustFluid,i,j,k)+ &
                                       SRhoUxDustGravity_IC(:,i,j,k)+ Source_VC(iRhoUxDust_I ,i,j,k)
      Source_VC(iRhoUyDust_I ,i,j,k) = SRhoUyDust_IC(1:nDustFluid,i,j,k)+ &
                                       SRhoUyDustGravity_IC(:,i,j,k)+ Source_VC(iRhoUyDust_I ,i,j,k)
      Source_VC(iRhoUzDust_I ,i,j,k) = SRhoUzDust_IC(1:nDustFluid,i,j,k)+ & 
                                       SRhoUzDustGravity_IC(:,i,j,k)+ Source_VC(iRhoUzDust_I ,i,j,k)
  end do;  end do;  end do

    

    if(iBlock == BlkTest .and. iProc==PROCtest) then
       call set_oktest('user_calc_sources',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if
   
    ! write(*,*) 'DoTest, DoTestMe:', DoTest, DoTestMe
    if(DoTestMe) then
       write(*,*)'user_calc_sources:'
       write(*,*)'Inputs: '
       i=iTest ; j=jTest ; k=kTest; 
123    format (A13,ES25.16,A15,A3,F7.2,A3)
       write(*,123)'x         = ',Xyz_DGB(x_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'y         = ',Xyz_DGB(y_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'z         = ',Xyz_DGB(z_,i,j,k,iBlock)," [rPlanet]"
       write(*,123)'r         = ',R_BLK(i,j,k,iBlock)," [rPlanet]"

       write(*,*) 'Massdust_I=',MassDust_I
       write(*,*) 'MassFluid_I=',MassFluid_I
       write(*,*) 'dustMassDistr=',dustMassDistr_I(Dust1_:Dust6_)/sum(dustMassDistr_I)
       write(*,*)'Neutrals:'
       do iNeutral=1,nNeutral
          write(*,123)'n_nIC       = ',NnNeutral_IC(iNeutral,i,j,k)," [m^-3]"
          write(*,123)'unxIC       = ',UnxNeutral_IC(iNeutral,i,j,k)," [m/s]"
          write(*,123)'unyIC       = ',UnyNeutral_IC(iNeutral,i,j,k)," [m/s]"
          write(*,123)'unzIC       = ',UnzNeutral_IC(iNeutral,i,j,k)," [m/s]"
       end do

       do iDust=1,nDust
          write(*,123)'n_dustIC       = ',NnDust_IC(iDust,i,j,k)," [m^-3]"
          write(*,123)'uxIC_dust      = ',UxDust_IC(iDust,i,j,k)," [m/s]"
          write(*,123)'uyIC           = ',UyDust_IC(iDust,i,j,k)," [m/s]"
          write(*,123)'uzIC           = ',UzDust_IC(iDust,i,j,k)," [m/s]"
          write(*,123)'SRhoUxDust     = ',SRhoUxDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]" 
          write(*,123)'SRhoUyDust     = ',SRhoUyDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]"
          write(*,123)'SRhoUzDust     = ',SRhoUzDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]"
          write(*,123)'SRhoUxDustGravity= ',SRhoUxDustGravity_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]"
          write(*,123)'SRhoUyDustGravity= ',SRhoUyDustGravity_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]" 
          write(*,123)'SRhoUzDustGravity= ',SRhoUzDustGravity_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_),"[kg m^-2/s^2]"       
          write(*,123)'gx=',Gravity_DCB(1,i,j,k,iBlock),"[m/s^2]"
          write(*,123)'gy=',Gravity_DCB(2,i,j,k,iBlock),"[m/s^2]" 
          write(*,123)'gz=',Gravity_DCB(3,i,j,k,iBlock),"[m/s^2]" 

          write(*,123)'DeltaUx=',SRhoUxDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_)&
                       /(NnDust_IC(iDust,i,j,k)*MassDust_I(iDust)*cProtonMass),"[m/s^2]"
          write(*,123)'DeltaUy=',SRhoUyDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_)&
                       /(NnDust_IC(iDust,i,j,k)*MassDust_I(iDust)*cProtonMass),"[m/s^2]"
          write(*,123)'DeltaUz=',SRhoUzDust_IC(iDust,i,j,k)*No2SI_V(UnitRhoU_)&
                       /(NnDust_IC(iDust,i,j,k)*MassDust_I(iDust)*cProtonMass),"[m/s^2]" 
      end do

     
           write(*,*)'Neutral     '

         
!           write(*,123)'RhoNeu       = ',State_VGB(NeuRho_,i,j,k,iBlock)*No2SI_V(UnitRho_)," [kg/m^3]"
!           write(*,123)'rhoUxNeu     = ',State_VGB(NeuRhoUx_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
!                " [kg/(m^2*s)]"
!           write(*,123)'rhoUyNeu     = ',State_VGB(NeuRhoUy_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
!                " [kg/(m^2*s)]"
!           write(*,123)'rhoUzNeu     = ',State_VGB(NeuRhoUz_,i,j,k,iBlock)*No2SI_V(UnitRhoU_),&
!                " [kg/(m^2*s)]"
!           write(*,123)'PNeu       = ',State_VGB(NeuP_,i,j,k,iBlock)*No2SI_V(UnitP_)," [Pa]"


     write(*,123)'SRhoNeuH2O= ',SRhoNeu_IC(1,i,j,k)*No2SI_V(UnitRho_)&
                   /No2SI_V(UnitT_)," [kg/(m^3*s)]"

          
     write(*,123)'SRhoUxNeuH2O    = ',SRhoUxNeu_IC(1,i,j,k)&
                              *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
     
     write(*,123)'SRhoUxNeuT1H2O    = ',SRhoUxNeuTerm_IIC(1,1,i,j,k)&
                               *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"

     write(*,123)'SRhoUxNeuT2H2O    = ',SRhoUxNeuTerm_IIC(2,1,i,j,k)&
                              *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"


     write(*,123)'SRhoUyNeuH2O= ',SRhoUyNeu_IC(1,i,j,k)*&
                           No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
         
     write(*,123)'SRhoUyNeuT1H2O= ',SRhoUyNeuTerm_IIC(1,1,i,j,k)&
                          *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
  
     write(*,123)'SRhoUyNeuT2H2O= ',SRhoUyNeuTerm_IIC(2,1,i,j,k)&
                          *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"


     write(*,123)'SRhoUzNeuH2O= ',SRhoUzNeu_IC(1,i,j,k)*No2SI_V(UnitRhoU_)&
                          /No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
     
     write(*,123)'SRhoUzNeuT1H2O= ',SRhoUzNeuTerm_IIC(1,1,i,j,k)&
                        *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"
 
     write(*,123)'SRhoUzNeuT2H2O= ',SRhoUzNeuTerm_IIC(2,1,i,j,k)&
                        *No2SI_V(UnitRhoU_)/No2SI_V(UnitT_)," [kg/(m^2*s^2)]"

     write(*,123)'SPNeuH2O= ',SPNeu_IC(1,i,j,k)*No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
     write(*,123)'SPNeuT1H2O= ',SPNeuTerm_IIC(1,1,i,j,k)*&
                                      No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
     write(*,123)'SPNeuT2H2O      = ',SPNeuTerm_IIC(2,1,i,j,k)*&
                                    No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
     write(*,123)'SPNeuT3H2O      = ',SPNeuTerm_IIC(3,1,i,j,k)*&
                                    No2SI_V(UnitP_)/No2SI_V(UnitT_)," [Pa/s]"
    
 
   end if 


  end subroutine user_calc_sources

  !========================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !----------------------------------------------------------------------

    !! Source terms are evaluated explicitly!
    !RETURN
    if(.not. allocated(iVarPointImpl_I)) &
         allocate(iVarPointImpl_I(5*nFluid))

  !  Martin's code
  !  do iFluid = 1, nIonFluid
  !     iVarPointImpl_I(5*iFluid-4) = iRhoIon_I(iFluid)
  !     iVarPointImpl_I(5*iFluid-3) = iRhoUxIon_I(iFluid)
  !     iVarPointImpl_I(5*iFluid-2) = iRhoUyIon_I(iFluid)
  !     iVarPointImpl_I(5*iFluid-1) = iRhoUzIon_I(iFluid)
  !     iVarPointImpl_I(5*iFluid)   = iPIon_I(iFluid)
  !  end do

   do iFluid = 1, nFluid
       iVarPointImpl_I(5*iFluid-4) = iRho_I(iFluid)
       iVarPointImpl_I(5*iFluid-3) = iRhoUx_I(iFluid)
       iVarPointImpl_I(5*iFluid-2) = iRhoUy_I(iFluid)
       iVarPointImpl_I(5*iFluid-1) = iRhoUz_I(iFluid)
       iVarPointImpl_I(5*iFluid)   = iP_I(iFluid)
    end do

    IsPointImplMatrixSet = .false.
    !IsAsymmetric= .false.

  end subroutine user_init_point_implicit

  !========================================================================

    subroutine user_set_ICs(iBlock)
    use ModIO,       ONLY: restart
    use ModProcMH,   ONLY: iProc
    use ModMain,     ONLY: iTest, jTest, kTest, ProcTest, BlkTest, Body1_, Body1
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics
    use ModConst,    ONLY: cBoltzmann
    use ModGeometry, ONLY: R_BLK,Xyz_DGB

    integer, intent(in) :: iBlock

    logical :: DoTest, DoTestMe=.true.
    integer :: i, j, k, iIonFluid
    ! !-------------------------------------------------------------------------
    if(iProc==PROCtest .and. iBlock==BLKtest)then
       call set_oktest('user_set_ICs', DoTest, DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    
    MassNeu_I=MassFluid_I(1:nNeuFluid)
    MassFluid_I(nNeuFluid+1:nFluid)=4/3*cPi*dustradius_I**3*dustDensity/cProtonMass
 
   

    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
     !  if ((Body1) .and. (R_BLK(i,j,k,iBlock) < bodyRadius*SI2No_V(UnitX_))) then
          ! State_VGB(:,i,j,k,iBlock) = CellState_VI(:,Body1_)
          ! State_VGB(iUx_I,i,j,k,iBlock) = 0.
          ! State_VGB(iUy_I,i,j,k,iBlock) = 0.
          ! State_VGB(iUz_I,i,j,k,iBlock) = 0.
          ! State_VGB(iPIon_I,i,j,k,iBlock) = BodyNDim_I*cBoltzmann*BodyTDim_I*1e6*SI2No_V(UnitP_)
          ! if (UseElectronPressure) then
          !    State_VGB(P_,i,j,k,iBlock) = sum(State_VGB(iPIon_I,i,j,k,iBlock))
          !    State_VGB(Pe_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock)*&
          !         ElectronPressureRatio
          ! else
          !    State_VGB(P_,i,j,k,iBlock) = &
          !         sum(State_VGB(iPIon_I,i,j,k,iBlock))*(1.+ElectronPressureRatio)
          ! end if
      ! else

           State_VGB(H2ORho_,i,j,k,iBlock)= 1e11*Si2No_V(UnitN_)
        !   State_VGB(H2ORho_,i,j,k,iBlock)= Qprod/(4.*cPi*(R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_))**2&
        !       *uNeutr_I(H2O_))*exp(-DestructRate_I(H2O_)/(0.9*uNeutr_I(H2O_))*R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_))&
      !         *Si2No_V(UnitN_)*MassNeu_I(H2O_)

      !     State_VGB(H2ORhoUx_,i,j,k,iBlock)= State_VGB(H2ORho_,i,j,k,iBlock)*&       !dimensionless 
      !         0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(x_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_) 
      !     State_VGB(H2ORhoUy_,i,j,k,iBlock)= State_VGB(H2ORho_,i,j,k,iBlock)*&       !dimensionless 
      !         0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(y_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_)               
      !     State_VGB(H2ORhoUz_,i,j,k,iBlock)= State_VGB(H2ORho_,i,j,k,iBlock)*&       !dimensionless 
       !        0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(z_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_)
                  
           State_VGB(H2ORhoUx_,i,j,k,iBlock)= 0.0
                                
           State_VGB(H2ORhoUy_,i,j,k,iBlock)= 0.0
                                
           State_VGB(H2ORhoUz_,i,j,k,iBlock)= 0.0
                                
        State_VGB(H2OP_,i,j,k,iBlock)= State_VGB(H2ORho_,i,j,k,iBlock)/MassNeu_I(H2O_)*& 
                             30*Io2No_V(UnitTemperature_)

     !   State_VGB(HRho_,i,j,k,iBlock)=Qprod*1.205E-11*1e6*(1E2*R_BLK(i,j,k,iBlock)*No2Si_V(UnitX_))**(-1.6103)&
     !                     *MassNeu_I(H_)*SI2No_V(UnitN_)
        State_VGB(HRho_,i,j,k,iBlock) = 1e8*SI2No_V(UnitN_) 
        State_VGB(HRhoUx_,i,j,k,iBlock)= 0.0
        State_VGB(HRhoUy_,i,j,k,iBlock)= 0.0
        State_VGB(HRhoUz_,i,j,k,iBlock)= 0.0
        State_VGB(HP_,i,j,k,iBlock)= State_VGB(HRho_,i,j,k,iBlock)/MassNeu_I(H_)*& 
                             30*Io2No_V(UnitTemperature_)
    !    State_VGB(HRhoUx_,i,j,k,iBlock)= State_VGB(HRho_,i,j,k,iBlock)*&       !dimensionless                   
    !       0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(x_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_)             
    !    State_VGB(HRhoUy_,i,j,k,iBlock)= State_VGB(HRho_,i,j,k,iBlock)*&       !dimensionless                   
    !       0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(y_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_)
    !    State_VGB(HRhoUz_,i,j,k,iBlock)= State_VGB(HRho_,i,j,k,iBlock)*&       !dimensionless                   
    !       0.5*0.7*uNeutr_I(H2O_)*Xyz_DGB(z_,i,j,k,iBlock)/R_BLK(i,j,k,iBlock)*SI2No_V(UnitU_)
  



        State_VGB(iRhoDust_I,i,j,k,iBlock)= State_VGB(H2ORho_,i,j,k,iBlock)*dust2gas*&
                           dustMassDistr_I(Dust1_:Dust6_)/sum(dustMassDistr_I)*1e-3
        
     !   State_VGB(iRhoUxDust_I,i,j,k,iBlock)=State_VGB(iRhoDust_I,i,j,k,iBlock)*&   !dimensionless 
     !                      0.5*UnxNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*SI2No_V(UnitU_)
     !   State_VGB(iRhoUyDust_I,i,j,k,iBlock)=State_VGB(iRhoDust_I,i,j,k,iBlock)*&   !dimensionless 
     !                      0.5*UnyNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*SI2No_V(UnitU_)
     !   State_VGB(iRhoUzDust_I,i,j,k,iBlock)=State_VGB(iRhoDust_I,i,j,k,iBlock)*&   !dimensionless 
     !                      0.5*UnzNeutral_IG(H2O_,i-MinI+1,j-MinJ+1,k-MinK+1)*SI2No_V(UnitU_)
 
        State_VGB(iRhoUxDust_I,i,j,k,iBlock)= 0.0
        State_VGB(iRhoUyDust_I,i,j,k,iBlock)= 0.0
        State_VGB(iRhoUzDust_I,i,j,k,iBlock)= 0.0
        State_VGB(iPDust_I,i,j,k,iBlock)= 1e-20*SI2No_V(UnitP_)  
        
     !  end if
    end do; end do ; end do
     
  



    if(DoTestMe) then
       i=iTest ; j=jTest ; k=kTest
123    format (A13,ES25.16,A15)
  
     write(*,*) 'iRhoDust_I', iRhoDust_I
     write(*,*) 'dustmassdistr_I',dustMassDistr_I
     write(*,*) 'mass/sum',dustMassDistr_I/sum(dustMassDistr_I)

    end if

  end subroutine user_set_ICs

  !============================================================================

  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModSize,         ONLY: x_,y_,z_
    use ModVarIndexes
    use ModFaceBoundary, ONLY: TimeBc, iFace, jFace, kFace, FaceCoords_D, &
         iBoundary, VarsTrueFace_V, iSide, iBlock => iBlockBc
    use ModMain,         ONLY: body1_, n_step, time_simulation, time_accurate, Dt
    use ModPhysics,      ONLY: LowDensityRatio,UnitTemperature_, Io2No_V,  &
        NO2SI_V, UnitP_, UnitRho_, UnitRhoU_, UnitU_, UnitB_, UnitN_, Io2SI_V, &
        rBody,cPi,rPlanetSI,SI2No_V,UnitX_
    use ModGeometry,      ONLY: R_BLK, Xyz_DGB,ExtraBc_
    use BATL_lib, ONLY: CellSize_DB  
    use ModNumConst, ONLY : cDegToRad, cRadToDeg
    use ModCoordTransform, ONLY: dir_to_xyz
    use ModBlockData, ONLY: use_block_data, clean_block_data, &
         get_block_data, put_block_data
    use ModConst, ONLY: cProtonMass

    logical :: FirstCall = .true., DoTest, DoTestMe=.true.
    !real    :: UdotR(nIonFluid), URefl_D(1:3,nIonFluid)
    !real    :: BdotR, BRefl_D(1:3)

    real, intent(out):: VarsGhostFace_V(nVar) 
    !=======================================================================
    logical :: DoTestHere=.true., IsIlluminated = .false.  

    integer:: iTrue, jTrue, kTrue, iBody, jBody, kBody

    real :: XyzIntersect_D(3), XyzStart_D(3), XyzEnd_D(3)
    real :: TestFace_D(3)
    real :: XyzTrueCell_D(3), XyzBodyCell_D(3)
    real :: Normal_D(3), CosAngle
    real :: TempCometLocal, uNormal, ProductionRateLocal
    real :: LonSunNow
    real :: FaceCoordsTest_D(3) = 0.0
    real :: maxLiftableSize, gravityNormal, cDrag=2.0
    integer:: nStepLonSun = -1
    real, save :: NormalSun_D(3)

   
    character(len=*), parameter:: NameSub = 'user_set_face_boundary'

    !=======================================================================

    !------------------------------------------------------------------------

    if(DoTestMe.and.FirstCall) then
       call set_oktest('user_set_face_boundary',DoTest,DoTestMe)
    else
       DoTest=.false.; DoTestMe=.false.
    end if

    MassNeu_I=MassFluid_I(1:nNeuFluid)
  
!=====================================================================================
    !! Outer boundaries
    if(iBoundary > 0) then

       VarsGhostFace_V = VarsTrueFace_V
          
    !! Body boundaries
    else if (iBoundary == ExtraBc_) then
     
  !   write(*,*) 'extraBc Test1'   !Yinsi
!==================================================================================
    ! We can use the saved values if 
    ! not too much time or time step has passed since the last save
    if (use_block_data(iBlock))  then
       call get_block_data(iBlock, nVar, VarsGhostFace_V)
     !  VarsGhostFace_V(HRho_:HP_) = VarsTrueFace_V(HRho_:HP_)
       RETURN
    end if

   
   
   !  write(*,*) 'extraBc test2' ! Yinsi
   !=====================================================================================

    ! Recalculate LonSunNow if needed. Note that we use Time_Simulation
    ! that only changes after the full time step is done
    ! and not ModFaceBoundary::TimeBc that can change during subcycling.
  
       call set_oktest(NameSub, DoTest, DoTestMe)

     
       ! Use input direction
          LonSunNow = LonSun
        if (DoTest) write(*,*) NameSub, &
               ': iProc, LonSun, LonSunNow=', iProc, LonSun, LonSunNow


       ! Get the direction vector to the Sun
       call dir_to_xyz((90-LatSun)*cDegToRad, LonSunNow*cDegToRad, NormalSun_D)


    ! Save step and simulation time info
    nStepSave_B(iBlock)          = n_step
    TimeSimulationSave_B(iBlock) = Time_Simulation

    ! Floating boundary condition by default
    VarsGhostFace_V = VarsTrueFace_V

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
    if (.not. is_segment_intersected( &
         XyzTrueCell_D, XyzBodyCell_D, IsOddIn = .true., &
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

    ! Set local outflow parameters as default that may be overwritten if illuminated
    TempCometLocal      = TempCometMin
    ProductionRateLocal = ProductionRateMin

    ! Check if Sun light hits the shape
    CosAngle = sum(Normal_D*NormalSun_D)

    if (CosAngle > 0.0) then

       ! See whether the intersection point is in the shade by going towards the Sun
       ! and checking for intersection with the shape
       XyzStart_D = XyzIntersect_D + 1e-9*rMaxShape*NormalSun_D
       XyzEnd_D   = XyzIntersect_D +    2*rMaxShape*NormalSun_D
       if(.not.is_segment_intersected(XyzStart_D, XyzEnd_D)) then
          IsIlluminated = .true.

          ! Increase temperature of the face if it is illuminated
          TempCometLocal      = max( TempCometMin, &
               SlopeTemp / CosAngle + bTemp)
          ! Increase neutral production rate
          ProductionRateLocal = max( ProductionRateMin, &
               SlopeProduction * CosAngle + bProduction )
       end if
    end if

    ! Calculate the normal velocity
    uNormal = sqrt(TempCometLocal)*TempToUnormal
    
    VarsGhostFace_V(Ux_:Uz_) = Normal_D*uNormal
    VarsGhostFace_V(Rho_)    = ProductionRateLocal/uNormal*MassFluid_I(1)
    VarsGhostFace_V(P_)      = &
         VarsGhostFace_V(Rho_)*TempCometLocal*TempToPressure

    VarsGhostFace_V(HUx_:HUz_) = Normal_D*uNormal
    VarsGhostFace_V(HRho_)= VarsGhostFace_V(Rho_)/1e7    
    VarsGhostFace_V(HP_)      = &
         VarsGhostFace_V(HRho_)*TempCometLocal*TempToPressure*18 
! since temp2pressure is for h2o


    VarsGhostFace_V(iRhoDust_I) = VarsGhostFace_V(H2ORho_)*dust2gas*&
                               dustMassDistr_I(Dust1_:Dust6_)/sum(dustMassDistr_I)
    VarsGhostFace_V(iPDust_I) = 1e-20*SI2No_V(UnitP_)

    gravityNormal = -sum(Normal_D*gravity_DCB(:,iFace,jFace,kFace,iBlock))  ! positve direction for gnormal is outward
    maxLiftableSize = 3/8.0*cDrag*ProductionRateLocal*No2SI_V(UnitN_)*No2SI_V(UnitU_)* &
                     cProtonMass*MassNeu_I(H2O_)*uNormal*No2SI_V(UnitU_)/&
                     (dustDensity*gravityNormal)
              
    if (addGravity==1) then
       where (dustRadius_I > maxLiftableSize)
       VarsGhostFace_V(iRhoDust_I)= VarsGhostFace_V(iRhoDust_I)*1e-3;
       end where
    end if    
   ! VarsGhostFace_V(iRhoUxDust_I)=Normal_D(1)*uNormal
     VarsGhostFace_V(iRhoUxDust_I)= 0.0
                             
   ! VarsGhostFace_V(iRhoUyDust_I)=Normal_D(2)*uNormal
     VarsGhostFace_V(iRhoUyDust_I)= 0.0                     
   ! VarsGhostFace_V(iRhoUzDust_I)=Normal_D(3)*uNormal
     VarsGhostFace_V(iRhoUzDust_I)= 0.0  

    if (.false. .and. DoTestHere .and. IsIlluminated .and. CosAngle > 0.5) then
       FaceCoordsTest_D = FaceCoords_D
       
       write(*,*) 'FaceCoords_D  =', FaceCoords_D
       write(*,*) 'XyzTrueCell_D =', XyzTrueCell_D
       write(*,*) 'XyzBodyCell_D =', XyzBodyCell_D
       write(*,*) 'XyzIntersect_D=', XyzIntersect_D
       write(*,*) 'XyzStart_D    =', XyzStart_D 
       write(*,*) 'XyzEnd_D      =', XyzEnd_D 
       write(*,*) 'Normal_D      =', Normal_D
       write(*,*) 'CosAngle      =', CosAngle
       write(*,*) 'ProductionRate=', ProductionRateLocal
       write(*,*) 'TempLocal [K] =', No2Si_V(UnitTemperature_)*TempCometLocal
       write(*,*) 'Temperature[K]=', No2Si_V(UnitTemperature_)* &
            VarsGhostFace_V(p_)*MassFluid_I(1)/VarsGhostFace_V(Rho_)
       write(*,*) 'Rho           =', VarsGhostFace_V(Rho_)
       write(*,*) 'u_D           =', VarsGhostFace_V(Ux_:Uz_)
       write(*,*) 'uNormal       =', uNormal
       write(*,*) 'p             =', VarsGhostFace_V(p_)
       write(*,*) 'Mach number   =', &
            uNormal/sqrt(Gamma*VarsGhostFace_V(p_)/VarsGhostFace_V(Rho_))
       DoTestHere=.false.
    end if

    IsIlluminated = .false.

    ! Store for future time steps
    call put_block_data(iBlock, nVar, VarsGhostFace_V)


    end if


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

    integer:: iTriangle, nIntersect, iIntersect, iMinRatio

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

 !=================================================================================
subroutine user_update_states(iStage,iBlock)
    use ModAdvance,  ONLY: State_VGB
    use ModPhysics,  ONLY: cBoltzmann, Si2No_V, &
         No2Si_V, UnitN_, UnitP_
    use ModEnergy,   ONLY: calc_energy_cell
    use ModGeometry, ONLY: r_BLK

    integer,intent(in) :: iStage, iBlock


    !----------------------------------------------------------------------                                                                                

    call update_states_MHD(iStage,iBlock)

    ! Enforce a constant Ht pressure                                                                                                                       

   !  State_VGB(Dust1P_,:,:,:,iBlock) = maxval(State_VGB(Dust1P_,:,:,:,iBlock))*No2SI_V(UnitN_)*&
   !                                  50*cBoltzmann*SI2No_V(UnitP_)
    
       State_VGB(iPDust_I,:,:,:,iBlock) = 1e-20*SI2No_V(UnitP_)

    call calc_energy_cell(iBlock)

  end subroutine user_update_states


!================================================================================

subroutine user_set_plot_var(iBlock, NameVar, IsDimensional,&
      PlotVar_G, PlotVarBody, UsePlotVarBody,&
      NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

   use ModAdvance,    ONLY: State_VGB, RhoUx_, RhoUy_, RhoUz_
   use ModPhysics,    ONLY: No2Si_V, Si2No_V, UnitP_, UnitN_, UnitU_, UnitT_, &
        ElectronCharge, ElectronPressureRatio
   use ModVarIndexes, ONLY: Rho_, P_, Pe_
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

   integer :: i, j, k


  !--------------------------------------------------------------------------

   IsFound = .true.

 
select case(NameVar)
case('dt')
      NameIdlUnit = 's'
      NameTecUnit = '[s]'
      PlotVar_G(:,:,:) = Dt_BLK(iBlock)*No2SI_V(UnitT_)

case('gx')
      NameIdlUnit = 'm/s^2'
      NameTecUnit = '[m/s^2]'
     
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
      PlotVar_G(i,j,k)= Gravity_DCB(1,i,j,k,iBlock)
    end do; end do; end do 
    
case('gy')
      NameIdlUnit = 'm/s^2'
      NameTecUnit = '[m/s^2]'
     
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
      PlotVar_G(i,j,k)= Gravity_DCB(2,i,j,k,iBlock)
    end do; end do; end do

case('gz')
      NameIdlUnit = 'm/s^2'
      NameTecUnit = '[m/s^2]'
     
    do k=MinK,MaxK; do j=MinJ,MaxJ; do i=MinI,MaxI
      PlotVar_G(i,j,k)= Gravity_DCB(3,i,j,k,iBlock)
    end do; end do; end do

 


case default
      IsFound = .false.
end select

   UsePlotVarBody = .false.
   PlotVarBody    = 0.0

end subroutine user_set_plot_var


!=================================================================================

end module ModUser
