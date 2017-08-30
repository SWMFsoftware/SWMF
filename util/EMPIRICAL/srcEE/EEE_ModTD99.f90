!  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module EEE_ModTD99
  use EEE_ModCommonVariables
  use ModHyperGeometric
  implicit none
  save
  private

  public :: set_parameters_TD99, get_transformed_TD99fluxrope

  !----------------------------------------------------------------------
  ! Logical variables related to the magnetic field computation::
  logical, public :: DoBqField=.false.

  logical :: UseVariedCurrent=.false.

  logical :: DoEquilItube=.false.
  logical :: DoRevCurrent=.false.
  logical :: DoMaintainEpot=.false.
  real    :: CurrentRiseTime,CurrentStartTime

  !----------------------------------------------------------------------
  ! Variables related to the position of the flux rope::
  real, parameter :: Li=0.5

  !----------------------------------------------------------------------
  ! Variables related to the flux rope properties::
  real :: ITube = 0.0, RTube  = 0.0, aTube = 0.0
  real :: Depth = 0.0, aRatio = 0.0
  real :: Mass = 0.0
  real :: InvH0=0.0, Rho0=0.0
  real :: ItubeSaved=0.0

  !----------------------------------------------------------------------
  ! Variables related to the properties of the strapping field, Bq::
  integer :: nStepSaved=-1
  real, parameter :: AlphaRamp=9.52381E-04         !in [-]
  real, parameter :: VTransX=1.500E+03             !in [m/s]
  real, parameter :: VTransY=-2.900E+04            !in [m/s]
  real, parameter :: UVorCMax0=2.5                 !in units of 100 km/s
  real, parameter :: BqZMax0=3.768210E+01          !in [Gauss]
  real :: BqZMax=0.0, BqZMaxSaved=0.0, UVorCMax=0.0
  real :: q=0.0, L=0.0

  !----------------------------------------------------------------------
  ! Declare the rotational matrix of coordinate transformation::
  real, dimension(3,3) :: Rotate_DD
  !\
  ! direction of the magnetic moment of the configuration
  real :: UnitX_D(3) 
contains
  !============================================================================

  subroutine set_parameters_TD99(NameCommand)
    use ModReadParam, ONLY: read_var

    character(len=*), intent(in):: NameCommand

    character(len=*), parameter:: NameSub = 'set_parameters_TD99'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case("#TD99FLUXROPE")
       call read_var('UseCme',              UseCme)
       call read_var('UseVariedCurrent'   , UseVariedCurrent)
       call read_var('CurrentStartTime'   , CurrentStartTime)
       call read_var('CurrentRiseTime '   , CurrentRiseTime)
       call read_var('DoAddFluxRope'      , DoAddFluxRope)
       call read_var('DoEquilItube'       , DoEquilItube)
       call read_var('DoRevCurrent'       , DoRevCurrent)
       call read_var('aratio'             , aratio)
       call read_var('Itube'              , Itube)
       call read_var('Rtube'              , Rtube)
       call read_var('aTube'              , aTube)
       call read_var('Depth'              , Depth)
       call read_var('Mass'               , Mass)
       call read_var('LongitudeCme'       , LongitudeCme)
       call read_var('LatitudeCme'        , LatitudeCme)
       call read_var('OrientationCme'     , OrientationCme)
       call read_var('DoBqField'          , DoBqField)
       call read_var('q'             , q)
       call read_var('L'             , L)
    case("#CME")
       call read_var('Current',     Itube)
       call read_var('RadiusMajor', Rtube)
       call read_var('RadiusMinor', aTube)
       call read_var('Depth',       Depth)
       call read_var('Mass',        Mass)
    case default
       call CON_stop(NameSub//' unknown NameCommand='//NameCommand)
    end select

  end subroutine set_parameters_TD99

  !============================================================================

  real function varied_current(Time)
    real, intent(in):: Time
    !--------------------------------------------------------------------------
    varied_current = min(&
         max(Time-CurrentStartTime,0.0)/CurrentRiseTime,1.0)*&
         ItubeSaved 
  end function varied_current

  !============================================================================

  real function time_increment(Time)
    real, intent(in):: Time
    !--------------------------------------------------------------------------
    time_increment = min(&
         max(Time-CurrentStartTime,0.0)/CurrentRiseTime,1.0)
  end function time_increment

  !============================================================================

  subroutine get_transformed_TD99fluxrope(RFace_D,BFRope_D,UVorT_D,&
       n_step,Iteration_Number,RhoFRope,Time)

    real, dimension(3), intent(in) :: RFace_D
    real, dimension(3), intent(out) :: BFRope_D,UVorT_D
    integer, intent(in) :: n_step,Iteration_Number
    real, intent(out), optional :: RhoFRope
    real, intent(in), optional :: Time

    logical:: DoFirstCall=.true.
    real:: atemp,Itemp
    real:: UVorR
    real, dimension(3):: B1FRopeTemp_D
    real, dimension(3):: R1Face_D,B1FRope_D,B1qField_D
    real, dimension(3):: UVorC_D,U1VorC_D
    !--------------------------------------------------------------------------
    ! Initialize the TD99 model parameters once::

    if (DoFirstCall) then
       call init_TD99_parameters
       DoFirstCall=.false.
    endif
    if (present(Time).and.UseVariedCurrent) &
         Itube = varied_current(Time)

    ! Check if the potential electric field needs to be applied.
    ! Maintain Epot for (CurrentStartTime < Time < CurrentRiseTime)::

    if (present(Time).and.DoBqField) then
       DoMaintainEpot = (&
            (Time.gt.CurrentStartTime).and.&
            (Time.lt.CurrentRiseTime))
    else
       DoMaintainEpot = .false.
    endif

    ! Compute the flux rope, and transform coordinates and vectors to
    ! position the flux rope in the desired way::

    R1Face_D = matmul(Rotate_DD,RFace_D)
    if (DoAddFluxRope) then
       call compute_TD99_FluxRope(R1Face_D,B1FRope_D,RhoFRope)
       if (DoRevCurrent) then
          Itemp = Itube; Itube = -Itemp
          atemp = aTube; aTube = aratio*atemp
          call compute_TD99_FluxRope(R1Face_D,B1FRopeTemp_D)
          B1FRope_D = B1FRope_D+B1FRopeTemp_D
          Itube = Itemp
          aTube = atemp
       endif
    else
       B1FRope_D = 0.0
       RhoFRope  = 0.0
    endif
    U1VorC_D = 0.0
    if (DoBqField) then
       if (present(Time)) then
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot,n_step,Iteration_Number,Time)
       else
          call compute_TD99_BqField(R1Face_D,B1qField_D,&
               U1VorC_D,DoMaintainEpot,n_step,Iteration_Number)
       endif
       B1FRope_D = B1FRope_D+B1qField_D
    endif
    BFRope_D = matmul(B1FRope_D,Rotate_DD)
    UVorC_D  = matmul(U1VorC_D,Rotate_DD)

    ! Compute the tangential component of the velocity field, UVorT_D,
    ! associated with the potential electric field, Epot::

    UVorR   = dot_product(RFace_D,UVorC_D)
    UVorT_D = UVorC_D-RFace_D*UVorR

    BFRope_D = BFRope_D*No2Si_V(UnitB_)
    UVorT_D = UVorT_D*No2Si_V(UnitU_)
    RhoFRope = RhoFRope*No2Si_V(UnitRho_)

  end subroutine get_transformed_TD99fluxrope

  !============================================================================

  subroutine init_TD99_parameters

    use ModCoordTransform, ONLY: rot_matrix_x,rot_matrix_y,rot_matrix_z

    real:: AlphaRope,LInduct,WFRope,FootSepar,ItubeDim
    !--------------------------------------------------------------------------

    ! Compute the magnetic energy, WFRope, associated with the portion
    ! of the flux rope current that is above the solar surface::

    InvH0 = cGravitation*Msun/Rsun*Si2No_V(UnitU_)**2   ! in [-]
    AlphaRope  = 2.0*acos(Depth/Rtube)                 ! in [rad]
    FootSepar  = Rtube*sin(0.5*AlphaRope)/1.0e6         ! in [Mm]
    LInduct    = cMu*(0.5*AlphaRope/cPi)*Rtube*(log(8.0 &
         *(Rtube-Depth)/aTube) - 1.75)            ! in [H]
    WFRope     = 0.5*LInduct*Itube**2*1.0e7             ! in [ergs]

    ! Compute the average density inside the flux rope assuming that the
    ! total amount of prominence mass is Mass (=10^16g=10^13kg)::

    Rho0  = Mass/(AlphaRope*Rtube*cPi*aTube**2)
    ! in [kg/m^3]

    ! Define the normalized model parameters here::

    ! Flux rope::
    Rtube = Rtube*Si2No_V(UnitX_)
    aTube = aTube*Si2No_V(UnitX_)
    ItubeDim   = Itube
    Itube = ItubeDim*Si2No_V(UnitJ_)*Si2No_V(UnitX_)**2 ! in [A]
    Rho0  = Rho0*Si2No_V(UnitRho_)

    ! Save the maximum value of the current for possible use in
    ! varied_current case::

    ItubeSaved = Itube 

    ! Strapping field::

    Depth     = Depth*Si2No_V(UnitX_)
    L     = L*Si2No_V(UnitX_)
    q     = q*Si2No_V(UnitB_)*Si2No_V(UnitX_)**2

    ! Construct the rotational matrix, Rotate_DD, to position the
    ! flux rope in the desired way on the solar surface::

    Rotate_DD = matmul(rot_matrix_y(-0.5*cPi),&
         rot_matrix_x(-OrientationCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_y(LatitudeCme*cDegToRad))
    Rotate_DD = matmul(Rotate_DD,          &
         rot_matrix_z(-LongitudeCme*cDegToRad))
    UnitX_D = (/1.0, 0.0, 0.0/)

    if (iProc==0) then
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'    Twisted Flux Rope Model by Titov & Demoulin, 1999.     '
       write(*,*) prefix
       write(*,*) prefix,'>>>>>>>>>>>>>>>>>>>                   <<<<<<<<<<<<<<<<<<<<<'
       write(*,*) prefix
       write(*,*) prefix,'Depth      = ',Depth*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'Rtube  = ', &
            Rtube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'aTube  = ', &
            aTube*No2Si_V(UnitX_)/1.0E6,'[Mm]'
       write(*,*) prefix,'atube/Rtube = ',aTube/Rtube,'[-]'
       write(*,*) prefix,'Itube  = ',ItubeDim,'[A]'
       write(*,*) prefix,'aratio = ',aratio,'[-]'
       write(*,*) prefix,'Mass   = ',Mass*1.0e3,'[g] '
       write(*,*) prefix,'Rho0   = ',Rho0*No2Io_V(UnitRho_),'[g/cm^3]'
       write(*,*) prefix
       write(*,*) prefix,'q      = ', &
            q*No2Si_V(UnitB_)*No2Si_V(UnitX_)**2,'[T m^2]'
       write(*,*) prefix,'L      = ',L*No2Si_V(UnitX_)/1.0e6,'[Mm]'
       write(*,*) prefix
       write(*,*) prefix,'Free energy of flux rope is ',WFRope,'Ergs.'
       write(*,*) prefix,'Separation of flux rope ends is ',FootSepar,'Mm,'
       write(*,*) prefix,'   or ',cPi*FootSepar*1.0e6/(2.0*Rsun)*cRadToDeg,'deg.'
       write(*,*) prefix
       if (UseVariedCurrent) then
          write(*,*) prefix,'>>>>>       UseVariedCurrent is set to .true.!!!      <<<<<'
          write(*,*) prefix,'CurrentStartTime = ',CurrentStartTime,'[s]'
          write(*,*) prefix,'CurrentRiseTime  = ',CurrentRiseTime,'[s]'
          write(*,*) prefix
       endif
    endif
    if (DoEquilItube) then

       ! Compute the equilibrium toroidal current, Itube, based
       ! on the force balance in direction normal to the surface of
       ! the flux tube.

       Itube = 8.0*cPi*q*L*Rtube &
            *(L**2+Rtube**2)**(-1.5) &
            /(alog(8.0*Rtube/aTube) &
            -1.5+Li/2.0)                           ! in [-]
       WFRope    = 0.5*LInduct*(ItubeDim)**2*1.0e7      ! in [ergs]
    endif
    if (DoEquilItube.and.iProc==0) then
       write(*,*) prefix,'The strapping field, Bq, is added and the EQUILIBRIUM value'
       write(*,*) prefix,'of Itube is computed!!!'
       write(*,*) prefix
       write(*,*) prefix,'The value of Itube is reset to :: ',Itube
       write(*,*) prefix,'The free energy of the flux rope is :: ',WFRope,'Ergs.'
       write(*,*) prefix
    endif

  end subroutine init_TD99_parameters

  !=====================================================================!

  subroutine compute_TD99_BqField(RFace_D,BqField_D,UVorC_D,&
       DoMaintainEpot,n_step,Iteration_Number,TimeNow)

    real, intent(in), dimension(3) :: RFace_D
    real, intent(out), dimension(3) :: BqField_D,UVorC_D
    logical, intent(in) :: DoMaintainEpot
    integer, intent(in) :: n_step,Iteration_Number
    real, intent(in), optional :: TimeNow

    ! Variables related to coordinates::
    real:: R2Plus,R2Mins
    real, dimension(3):: RPlus_D,RMins_D

    ! Variables related to computations of potential electric field::
    real:: BqZOverBqZ0,BqZFunction
    real, dimension(3):: EpotC_D,UTranC_D
    real, dimension(3):: GradBqZ_D,GradPsiC_D
    !--------------------------------------------------------------------

    ! Compute the locations, RMins_D and RPlus_D, of the two magnetic
    ! charges, -/+q::

    if (present(TimeNow)) then
       RPlus_D(x_) = RFace_D(x_)-L &
            - VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(x_) = RFace_D(x_)+L &
            + VTransX*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RPlus_D(y_) = RFace_D(y_) &
            - VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
       RMins_D(y_) = RFace_D(y_) &
            + VTransY*(TimeNow-CurrentStartTime)*Si2No_V(UnitX_)
    else
       RPlus_D(x_) = RFace_D(x_)-L
       RMins_D(x_) = RFace_D(x_)+L
       RPlus_D(y_) = RFace_D(y_)
       RMins_D(y_) = RPlus_D(y_)
    endif
    RPlus_D(z_) = RFace_D(z_) + Depth -1.0
    RMins_D(z_) = RPlus_D(z_)
    R2Plus = sqrt(dot_product(RPlus_D,RPlus_D))
    R2Mins = sqrt(dot_product(RMins_D,RMins_D))

    ! Compute the field of the strapping magnetic field, BqField_D::

    BqField_D = q*(RPlus_D/R2Plus**3 - RMins_D/R2Mins**3)

    ! Update the values of BqZMax and BqZMaxSaved once at each iteration::

    if (n_step/=nStepSaved) then
       nStepSaved  = n_step
       BqZMaxSaved = BqZMax
       BqZMax  = 0.0
    end if
    if (iteration_number==1) &
         BqZMaxSaved = BqZMax0/No2Io_V(UnitB_)
    BqZMax = max(abs(BqField_D(z_)),BqZMax)

    ! Apply Epot only if DoMaintainEpot=.true.!!!

    if (DoMaintainEpot) then

       ! Compute the gradient of the z-component of the Bq field on the
       ! solar surface in Cartesian geometry -- GradBqZ_D(x_:z_)::

       GradBqZ_D = 3.0*q*RMins_D(z_)*&
            (RMins_D/R2Mins**5 - RPlus_D/R2Plus**5)
       GradBqZ_D(z_)    = 0.0
 
       ! Compute the gradient of the scalar potential in Cartesian
       ! geometry -- GradPsiC_D::

       BqZOverBqZ0 = min(1.0,abs(BqField_D(z_))/BqZMaxSaved)
       BqZFunction = max(0.0,1.0-BqZOverBqZ0**2)
       if (BqZOverBqZ0.gt.0.0) then
          GradPsiC_D = GradBqZ_D*2.0 &
               *AlphaRamp*BqZFunction*BqZOverBqZ0*exp(-BqZFunction)
       else
          GradPsiC_D = 0.0
       endif

       ! Compute the potential electric field on the solar surface to
       ! be applied at one of the spots -- EpotC_D(x_:z_).
       ! This is given by:
       ! EpotC_D(x_:z_) = sign(BqField_D(z_))*BqZOverBqZ0*GradPsiC_D::

       EpotC_D = sign(1.0,BqField_D(z_))*GradPsiC_D*BqZOverBqZ0

       ! Compute the plasma velocity, UVorC_D, associated with the static
       ! eletric field, EpotC_D.
       ! This is given by:
       ! UVorC_D(x_:z_) = sign(BqField_D(z_))*GradPsiC_D(x_:z_) X Ez,
       ! where Ez = (0,0,1)

       UVorC_D(x_)  =  GradPsiC_D(y_)*sign(1.0,BqField_D(z_))
       UVorC_D(y_)  = -GradPsiC_D(x_)*sign(1.0,BqField_D(z_))
       UVorC_D(z_)  =  0.0
       UVorC_D      = UVorC_D*UVorCMax0
       if (iteration_number.gt.2) &
            UVorC_D = UVorC_D*UVorCMax0/UVorCMax

       ! Compute the translational velocity, UTranC_D, at which the two
       ! q-sources are moved with respect to each other.

       if (present(TimeNow)) then
          UTranC_D(x_) = (VTransX*Si2No_V(UnitU_))*BqZOverBqZ0 &
               *sign(1.0,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(y_) = (VTransY*Si2No_V(UnitU_))*BqZOverBqZ0 &
               *sign(1.0,BqField_D(z_))*exp(-BqZFunction)
          UTranC_D(z_) = 0.0
       else
          UTranC_D = 0.0
       endif
       ! Add the translational velocity, UTranC_D, to UVorC_D::
       UVorC_D = UVorC_D+UTranC_D
    else
       EpotC_D = 0.0
       UVorC_D = 0.0
    endif

  end subroutine compute_TD99_BqField

  !=====================================================================!

  subroutine compute_TD99_FluxRope(RFace_D,BFRope_D,RhoFRope)
    use ModCoordTransform, ONLY: cross_product
    !\__                                                             __/!
    !    Twisted Magnetic Field Configuration by Titov & Demoulin '99   !
    !                                                                   !
    ! An instability that causes a CME eruption is expected to occur at !
    ! R > L*sqrt(2). For a detailed description of the initial state    !
    ! refer to A&A, 1999, v.351, pp.707-720                             !
    !                                                                   !
    ! ___  This module was written by Ilia Roussev on June 10, 2002 ___ !
    !/                                                                 \!

    real, intent(in), dimension(3):: RFace_D
    real, intent(out), dimension(3):: BFRope_D

    real, intent(out), optional:: RhoFRope

    !\
    ! Coordinates relative to the configuration center:
    !/
    real :: XyzRel_D(3)
    !\
    ! Distance to the configuration center squared:
    !/
    real :: R2 !=sum(XyzRel_D**2)

    real:: xxx, R2Face
    real:: RMinus, RPlus2, Rperp
    real:: ThetaUVy,ThetaUVz
    real:: Kappa,dKappadx,dKappadr
    real:: KappaA,dKappaAdr
    !Values of the hypergeometric functions, 
    !2F1(3/2, 3/2 ; 3; Kappa**2) and 2F1(1/2, 1/2; 1; Kappa**3)
    real:: F32323, F12121

    ! Complete elliptic integrals of related variables::
    real:: KElliptic, EElliptic
    ! Vector potential related variables::
    real:: Ak,dAkdk
    real:: AkA,dAkdkA,d2Akdk2A
    real:: AI,dAIdx,dAIdr
    ! Flux-rope related variables::
    real:: BIPhi_D(3), B0_D(3), B0
    !--------------------------------------------------------------------
    ! Assign X,Y,Z coordinates at which to compute the magnetic field::
    !\
    ! xxx - coordinate along the axis of symmetry, equal to 0 at the 
    ! plane of symmetry.
    ! zzz - is the heliocetric coordinate along the line passing through 
    ! the center of configuration, which is at the depth d below the 
    ! photosphere level.
    ! yyy coordinate is equal to zero at the axis of symmetry. 
    XyzRel_D = RFace_D - (/0.0, 0.0, 1 - Depth/)
    R2 = sum(XyzRel_D**2); xxx = sum(XyzRel_D*UnitX_D)
    R2Face = sqrt(dot_product(RFace_D,RFace_D))

    !\
    ! Field in the center of configuration
    !/
    B0 = 0.5*ITube/RTube; B0_D = B0*UnitX_D
    ! Compute Rperp and TubalDist::

    Rperp = sqrt(R2 - xxx**2)
    RMinus = sqrt(xxx**2 + (Rperp - Rtube)**2)
    RPlus2 = (Rperp + Rtube)**2 + xxx**2


    ! Define the model input, Kappa

    Kappa = 2.0*sqrt(Rperp*Rtube/RPlus2)
    if (RMinus.ge.aTube) then
       ! Compute the field and density outside the current torus   
       !\
       !Initialize redundant output
       !/
       if (present(RhoFRope))RhoFRope=0.0    
 
       ! Compute the vector potential, Ak, of the magnetic field 
       ! produced by the ring current Itube and its derivatives
       if(Kappa < 0.7)then       
          F32323 = hypergeom(1.50, 1.50, 3.0, Kappa**2)
          Ak     = 0.250*F32323
          F12121 = hypergeom(0.50, 0.50, 1.0, Kappa**2)
          dAkDk  = (F12121 - 0.1250*(2.0 - Kappa**2)*F32323)/&
               (1.0 - Kappa**2)
       else
          ! Truncate the value of Kappa:: 
          if (abs(1.0-Kappa).lt.cTiny/10.0) &
               Kappa = 1.0-cTiny/10.0

          ! Compute the vector potential in the internal, AIin, and
          ! external (outside the current torus), AIex, regions::   
          
          call calc_elliptic_int_1kind(Kappa,KElliptic)
          call calc_elliptic_int_2kind(Kappa,EElliptic)
          Ak  =(4.0/cPi)*((2.0-Kappa**2)*KElliptic - 2.0*EElliptic)/Kappa**4
          !\
          ! Calculate derivative of (Ak*k*3) over k using formulae:
          ! dK/dk = E/(k*(1-k^2)) - K/k, dE/dk = (E - K)/k
          ! Then, divide by k^2.
          !/
          dAkdk    = (4.0/cPi)*((2.0-Kappa**2)*EElliptic/(1.0-Kappa**2) &
               - 2.0*KElliptic)/Kappa**4
       end if        
          ! Obtain the BI field in the whole space from the corresponding
          ! vector potential, AI -->
          ! BI = curl(AI*ThetaUV) = BFRope_D(x_:z_)::
       
       BFRope_D =B0*(Rtube/sqrt(RPlus2))**3*&
            (dAkDk*(2*xxx*XyzRel_D + (RTube**2 - R2)*UnitX_D)/RPlus2 &
            +Ak*UnitX_D)
       !No toroidal field outside the filament
       BIPhi_D = 0.0
    else
       !\
       ! Compute the field and density inside the current torus
       !/
       ! 1.
       ! Add the prominence material inside the flux rope, assuming that the
       ! total amount mass is 10^13kg, and that the desnity scale-height is
       ! the same as the pressure scale-height, 1/InvH0 (i.e., iso-thermal
       ! atmoshpere)::
       
       if (present(RhoFRope))&
            RhoFRope = Rho0*exp(-10.0*(RMinus/aTube)**6) &
            *exp(-InvH0*abs(R2Face-1.0))
       !2.    
       ! Define the model input, KappaA. A given point is charakterized by
       ! two coordinates, say, RPepr and RMinus. In this case,
       ! if we denote Kappa=known_function(RPerp,RMinus), then 
       ! KappaA=known_function(RPepr,ATube)
       KappaA = 2.0*sqrt(Rperp*Rtube &
            /(4.0*Rperp*Rtube+aTube**2))
       call calc_elliptic_int_1kind(KappaA,KElliptic)
       call calc_elliptic_int_2kind(KappaA,EElliptic)
       !\
       ! This is the vector potential and its k derivative at the boundary
       !/
       Ak      = ((2.0-KappaA**2)*KElliptic &
            - 2.0*EElliptic)/KappaA
       dAkdk   = (2.0-KappaA**2)*EElliptic &
            /(KappaA**2*(1.0-KappaA**2)) &
            - 2.0*KElliptic/KappaA**2
       ! Compute the vector potential, AI by smpoothly continuing the 
       ! value from the boundary
       AI = Itube/(2.0*cPi)*sqrt(Rtube/Rperp)*&
            (Ak + dAkdk*(Kappa - KappaA))
       d2Akdk2A = ((7.0*KappaA**2-4.0 &
            - KappaA**4)*EElliptic/(1.0-KappaA**2) &
            + (4.0-5.0*KappaA**2)*KElliptic) &
            /(KappaA**3*(1.0-KappaA**2))
       ! Derive the BI field from the corresponding vector potential,
       ! AI (this involves the comp. of some nasty derivatives)::


       dKappaAdr = KappaA*aTube**2/(2.0*Rperp &
            *(4.0*Rperp*Rtube+aTube**2))
       dKappadx  = -xxx*Kappa/(RPerp*RPlus2)
       dKappadr  = Kappa*(Rtube**2 - R2) &
            /(2.0*Rperp*RPlus2)

       ! Derivative of AI with respect to `x` and `rperp`:: 

       dAIdx   = Itube/(2.0*cPi)*sqrt(Rtube/Rperp) &
            *(dAkdk*dKappadx)
       dAIdr   = Itube/(2.0*cPi)*sqrt(Rtube/Rperp) &
            *(dAkdk*dKappadr+d2Akdk2A*dKappaAdr &
            *(Kappa-KappaA))
       ! Obtain the BI field in the whole space from the corresponding
       ! vector potential, AI -->
       ! BI = curl(AI*ThetaUV) = BFRope_D(x_:z_)::

       BFRope_D = -dAIdx*XyzRel_D + (dAIdr+AI/(2.0*Rperp))*UnitX_D
       ! Compute the toroidal field (BIphix, BIphiy, BIphiz)
       ! produced by the azimuthal current Iphi. This is needed to ensure
       ! that the flux rope configuration is force free. 
       BIPhi_D = abs(Itube)/(2.0*cPi*RPerp*aTube**2) &
            *sqrt(2.0*(aTube**2-RMinus**2))*&
            cross_product(UnitX_D,XyzRel_D)
    end if
    ! Add the field of the azimuthal current, Iphi::
    ! Compute the field produced by the ring current, Itube, both
    ! inside and outside the torus, BI = BFRope_D(x_:z_)::
    BFRope_D = BFRope_D + BIPhi_D
  end subroutine compute_TD99_FluxRope

end module EEE_ModTD99
