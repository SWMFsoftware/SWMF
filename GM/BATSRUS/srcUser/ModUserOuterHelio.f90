!  Copyright (C) 2002 Regents of the University of Michigan, 
!  portions used with permission 
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!April 30, 2007 implementing MultiFluid
!June 01, 2007 correcting normalization
!June 08, 2007 correcting source terms
!August 18, 2007 Implementing 3-fluids
!October 23, 2007 more little corrections
!January 01-04, 2008 Source terms in point-implicit form - with help
! of Gabor Toth
!January 08 all source terms for PopI implicit
!January 25 comments added
!January 31-Feb 02
!April 28, 2008 Michigan visit
!May 02-setting the initial conditions to be closer to the 
!final solution
!May 30-31
!June30
!July 20-WORKED! 
!July 28-4 neutral fluids
!September 28 - source terms as McNutt
!March 25, 2015 - adapted to be solved with either multifluid neutrals or single
!ion MHD used for PT-OH coupling by A. Michael with help from G. Toth
!July 2015 - trying not to have body at 30AU for the neutrals
!September 2015 - A.Michael added time dependent code to user_update_states
! in subroutine calc_time_dep_sw, it reads in file solarwind.dat
!note TimeDep code was implemented without the Heliopause function
!September 2015 - Added user_init_session
!==============================================================================
module ModUser

  use ModSize,       ONLY: nI, nJ, nK, MinI, MaxI, MinJ, MaxJ, MinK, MaxK
  use ModMain,       ONLY: body1_, PROCtest, BLKtest, iTest, jTest, kTest, &
       nBlock, Unused_B, Time_Simulation
  use ModPhysics,    ONLY: Gamma, GammaMinus1, OmegaBody, &
       UnitX_, Io2Si_V, Si2Io_V, Si2No_V, No2Io_V, No2Si_V, Io2No_V, &
       NameTecUnit_V, NameIdlUnit_V, UnitAngle_, UnitDivB_, UnitEnergyDens_, &
       UnitJ_, UnitN_, UnitRho_, UnitU_, rBody, UnitB_, UnitP_, &
       UnitTemperature_, UnitT_, UnitRhoU_
  use ModNumConst,      ONLY: cRadToDeg, cTwoPi
  use ModConst,         ONLY: cBoltzmann, cProtonMass
  use ModAdvance,    ONLY: State_VGB, Source_VC, ExtraSource_ICB
  use ModGeometry,   ONLY: Xyz_DGB, r_BLK, true_cell
  use ModVarIndexes, ONLY: nVar, Rho_, Ux_, Uy_, Uz_, RhoUx_, RhoUy_, RhoUz_, &
       Bx_, By_, Bz_, p_, Energy_, iRho_I, iRhoUx_I, iRhoUy_I, iRhoUz_I, iP_I, &
       nFluid, NameVar_V
  use ModProcMH,     ONLY: iProc 
  use ModMultiFluid, ONLY: UseNeutralFluid, RhoNeutralsISW, RhoNeutralsISW_dim, &
       PNeutralsISW, PNeutralsISW_dim, UxNeutralsISW, UyNeutralsISW, UzNeutralsISW, &
       TNeutralsISW_dim, UxNeutralsISW_dim, UyNeutralsISW_dim, UzNeutralsISW_dim, &
       mProtonMass, iFluid, iRho, iRhoUx, iRhoUy, iRhoUz, iP, iEnergy, &
       RhoBcFactor_I, uBcFactor_I, select_fluid
  use ModUserEmpty,                                     &
       IMPLEMENTED1  => user_read_inputs,               &
       IMPLEMENTED2  => user_set_face_boundary,         &
       IMPLEMENTED3  => user_normalization,             &
       IMPLEMENTED4  => user_set_cell_boundary,         &
       IMPLEMENTED5  => user_set_ics,                   &
       IMPLEMENTED6  => user_initial_perturbation,      &
       IMPLEMENTED7  => user_update_states,             &
       IMPLEMENTED8  => user_action,                    & 
       IMPLEMENTED9  => user_io_units,                  &
       IMPLEMENTED10 => user_set_plot_var,              &
       IMPLEMENTED11 => user_calc_sources,              &
       IMPLEMENTED12 => user_init_point_implicit,       &
       IMPLEMENTED13 => user_init_session,              &
       IMPLEMENTED14 => user_set_boundary_cells

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 2.3
  character (len=*), parameter :: &
       NameUserModule = 'Outer Heliosphere with 4 neutrals, Opher & Toth'

  !these are the variable used for multiflud neutrals 
  !that are not defined in ModEquationMhd.f90 used for K-MHD
  integer, parameter :: &
       NeuRho_    = min(nVar, 9), &
       NeuRhoUx_  = min(nVar-2, 10), NeuUx_ = NeuRhoUx_, &
       NeuRhoUy_  = min(nVar-1, 11), NeuUy_ = NeuRhoUy_, &
       NeuRhoUz_  = min(nVar  , 12), NeuUz_ = NeuRhoUz_, &
       NeuP_      = min(nVar  , 13), &
       Ne2Rho_    = min(nVar  , 14), &
       Ne2RhoUx_  = min(nVar-2, 15), Ne2Ux_ = Ne2RhoUx_, &
       Ne2RhoUy_  = min(nVar-1, 16), Ne2Uy_ = Ne2RhoUy_, &
       Ne2RhoUz_  = min(nVar  , 17), Ne2Uz_ = Ne2RhoUz_, &
       Ne2P_      = min(nVar  , 18), &
       Ne3Rho_    = min(nVar  , 19), &
       Ne3RhoUx_  = min(nVar-2, 20), Ne3Ux_ = Ne3RhoUx_, &
       Ne3RhoUy_  = min(nVar-1, 21), Ne3Uy_ = Ne3RhoUy_, &
       Ne3RhoUz_  = min(nVar  , 22), Ne3Uz_ = Ne3RhoUz_, &
       Ne3P_      = min(nVar  , 23), &
       Ne4Rho_    = min(nVar  , 24), &
       Ne4RhoUx_  = min(nVar-2, 25), Ne4Ux_ = Ne4RhoUx_, &
       Ne4RhoUy_  = min(nVar-1, 26), Ne4Uy_ = Ne4RhoUy_, &
       Ne4RhoUz_  = min(nVar  , 27), Ne4Uz_ = Ne4RhoUz_, &
       Ne4P_      = min(nVar  , 28)

  ! Named indexes for fluids
  integer, parameter :: Ion_ = 1, Neu_ = min(nFluid,2), &
       Ne2_ = min(nFluid,3), Ne3_ = min(nFluid,4), Ne4_= min(nFluid,5) 

  logical :: UseSource_I(Ion_:Ne4_) = .true.

  real :: OmegaSun   = 0.0  ! normalized angular speed of Sun
  real :: ParkerTilt = 0.0  ! Bphi/Br at the equator at r=rBody

  integer :: iTableSolarWind = -1 ! initialization is needed

  ! SWH variables.
  !/
  real :: &
       SWH_rho=0.0, SWH_rho_dim=0.0, &
       SWH_p=0.0  , SWH_T_dim  =0.0, &
       SWH_Ux=0.0 , SWH_Ux_dim=0.0 , &
       SWH_Uy=0.0 , SWH_Uy_dim=0.0 , &
       SWH_Uz=0.0 , SWH_Uz_dim=0.0 , &
       SWH_Bx=0.0 , SWH_Bx_dim=0.0 , &
       SWH_By=0.0 , SWH_By_dim=0.0 , &
       SWH_Bz=0.0 , SWH_Bz_dim=0.0

  ! 
  ! VLISM variables.
  !/
  real ::      VLISW_T_dim=0.0  , &
       VLISW_a_dim=0.0  , &
       VLISW_rho=0.0 , VLISW_rho_dim=0.0, &
       VLISW_p=0.0  , VLISW_p_dim=0.0   , &
       VLISW_Ux=0.0 , VLISW_Ux_dim=0.0 , &
       VLISW_Uy=0.0  , VLISW_Uy_dim=0.0 , &
       VLISW_Uz=0.0 , VLISW_Uz_dim=0.0 , &
       VLISW_Bx=0.0 , VLISW_Bx_dim=0.0 , &
       VLISW_By=0.0 , VLISW_By_dim=0.0 , &
       VLISW_Bz=0.0 , VLISW_Bz_dim=0.0 , &
       VLISW_B_factor=0.0

  real :: VLISW_p_dim1=0.0, VLISW_p1=0.0

  ! Fast solar wind
  real:: &
       SWfast_rho=0.0, SWfast_rho_dim=0.0, &
       SWfast_p=0.0  , SWfast_p_dim=0.0  , &
       SWfast_Ux=0.0 , SWfast_Ux_dim=0.0 , &
       SWfast_Uy=0.0 , SWfast_Uy_dim=0.0 , &
       SWfast_Uz=0.0 , SWfast_Uz_dim=0.0

  ! mass of neutrals in 
  real :: mNeutrals

  ! Velocity, temperature, Mach number and radius limits for the populations
  real :: TempPop1LimitDim = 1e5    ! [K]
  real :: uPop1LimitDim    = 100.0  ! [km/s]
  real :: MachPop2Limit    = 0.9
  real :: MachPop3Limit    = 1.5
  real :: MachPop4Limit    = 2.0
  real :: rPop3Limit       = 50.0   ! [AU] it is all Pop3 out to rPop3Limit

  integer :: iFluidProduced_C(nI, nJ, nK)

  ! Variables for the reflective shape
  real:: rCylinder = -1.0, zCylinder = -1.0
  real:: rCrescent = -1.0, xCrescentCenter = -1.0

contains

  !=========================================================================

  subroutine user_read_inputs

    use ModReadParam

    character (len=100) :: NameCommand
    character (len=*), parameter :: NameSub = 'user_read_inputs'
    !-------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT
       case("#SOLARWINDH")
          call read_var('SWH_rho_dim',SWH_rho_dim)
          call read_var('SWH_T_dim'  ,SWH_T_dim)
          call read_var('SWH_Ux_dim' ,SWH_Ux_dim)
          call read_var('SWH_Uy_dim' ,SWH_Uy_dim)
          call read_var('SWH_Uz_dim' ,SWH_Uz_dim)
          call read_var('SWH_Bx_dim' ,SWH_Bx_dim)
          call read_var('SWH_By_dim' ,SWH_By_dim)
          call read_var('SWH_Bz_dim' ,SWH_Bz_dim)
       case("#SOLARWINDFAST")
          call read_var('SWfast_rho_dim',SWfast_rho_dim)
          call read_var('SWfast_Ux_dim' ,SWfast_Ux_dim)
          call read_var('SWfast_Uy_dim' ,SWfast_Uy_dim)
          call read_var('SWfast_Uz_dim' ,SWfast_Uz_dim)
       case("#GLOBALHELIOSPHERE")
          call read_var('Rbody',Rbody)
       case("#VLISW")
          call read_var('VLISW_rho_dim' ,VLISW_rho_dim)
          call read_var('VLISW_T_dim'  ,VLISW_T_dim)
          call read_var('VLISW_Ux_dim' ,VLISW_Ux_dim)
          call read_var('VLISW_Uy_dim' ,VLISW_Uy_dim)
          call read_var('VLISW_Uz_dim' ,VLISW_Uz_dim)
          call read_var('VLISW_Bx_dim' ,VLISW_Bx_dim)
          call read_var('VLISW_By_dim' ,VLISW_By_dim)
          call read_var('VLISW_Bz_dim' ,VLISW_Bz_dim)
       case("#SOURCES")
          call read_var('UseIonSource', UseSource_I(Ion_))
          call read_var('UseNeuSource', UseSource_I(Neu_))
          call read_var('UseNe2Source', UseSource_I(Ne2_))
          call read_var('UseNe3Source', UseSource_I(Ne3_))
          call read_var('UseNe4Source', UseSource_I(Ne4_))
       case("#REGIONS")
          call read_var('TempPop1LimitDim', TempPop1LimitDim)
          call read_var('uPop1LimitDim',    uPop1LimitDim)
          call read_var('MachPop2Limit',    MachPop2Limit)
          call read_var('MachPop3Limit',    MachPop3Limit)
          call read_var('rPop3Limit',       rPop3Limit)
          call read_var('MachPop4Limit',    MachPop4Limit)
          ! case("#FACTORS")
          !    call read_var('RhoNeuFactor', RhoNeuFactor)
          !    call read_var('uNeuFactor'  , uNeuFactor)
          !    call read_var('RhoNe2Factor', RhoNe2Factor)
          !    call read_var('uNe2Factor'  , uNe2Factor)
          !    call read_var('RhoNe3Factor', RhoNe3Factor)
          !    call read_var('uNe3Factor'  , uNe3Factor)
          !    call read_var('RhoNe4Factor', RhoNe4Factor)
          !    call read_var('uNe4Factor'  , uNe4Factor)
       case("#INNERCYLINDER")
          call read_var('rCylinder', rCylinder)
          call read_var('zCylinder', zCylinder)
       case("#INNERCRESCENT")
          call read_var('rCrescent', rCrescent)
          call read_var('xCrescentCenter', xCrescentCenter)
       case default
          if(iProc==0) call stop_mpi(NameSub// &
               ': unrecognized command: '//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !==========================================================================
  subroutine user_set_face_boundary(VarsGhostFace_V)

    use ModFaceBoundary, ONLY: iBoundary, FaceCoords_D, VarsTrueFace_V, &
         iFace, jFace, kFace, iBlockBc

    use ModCoordTransform, ONLY: rot_xyz_sph

    real, intent(out):: VarsGhostFace_V(nVar)

    ! local variables
    real:: xFace, yFace, zFace
    real:: SinTheta
    real:: Bsph_D(3), Vsph_D(3)

    real :: pSolarWind, Pmag, PmagEquator

    real :: XyzSph_DD(3,3) ! rotation matrix Xyz_D = matmul(XyzSph_DD, Sph_D)

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='user_set_face_boundary'
    !-------------------------------------------------------------------
    if(iBoundary /= body1_ .and. iBoundary /= 1)then
       write(*,*) NameSub,': iBoundary=', iBoundary
       call stop_mpi(NameSub//' is not implemented for this boundary!')
    end if

    if(iProc == ProcTest .and. iBlockBc == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(DoTestMe)write(*,*) NameSub,' starting with iBoundary=', iBoundary

    if(iBoundary == 1)then
       ! The LISW enters at the 1st boundary (negative x)

       ! Ion inflow
       VarsGhostFace_V(Rho_) = VLISW_rho
       VarsGhostFace_V(Ux_)  = VLISW_Ux
       VarsGhostFace_V(Uy_)  = VLISW_Uy
       VarsGhostFace_V(Uz_)  = VLISW_Uz
       VarsGhostFace_V(Bx_)  = VLISW_Bx
       VarsGhostFace_V(By_)  = VLISW_By
       VarsGhostFace_V(Bz_)  = VLISW_Bz
       VarsGhostFace_V(p_ )  = VLISW_p

       if(UseNeutralFluid)then
          ! Inflow for the neutral Pop IV

          ! Use supersonic outflow by for most fluids
          VarsGhostFace_V = VarsTrueFace_V

          VarsGhostFace_V(Ne4Rho_) = RhoNeutralsISW
          VarsGhostFace_V(Ne4Ux_ ) = UxNeutralsISW
          VarsGhostFace_V(Ne4Uy_ ) = UyNeutralsISW
          VarsGhostFace_V(Ne4Uz_ ) = UzNeutralsISW
          VarsGhostFace_V(Ne4P_  ) = PNeutralsISW
       endif

       RETURN
    end if

    ! Make sure that OmegaSun and ParkerTilt are set
    if(OmegaSun == 0.0) call set_omega_parker_tilt

    XyzSph_DD = rot_xyz_sph(FaceCoords_D)

    xFace = FaceCoords_D(1)
    yFace = FaceCoords_D(2)
    zFace = FaceCoords_D(3)

    SinTheta = sqrt((xFace**2 + YFace**2)/(xFace**2 + yFace**2 + zFace**2))

    ! calculating the spherical Parker field components
    ! SWH_Bx is the value of the field at the pole B0

    ! Note: use -zFace to invert polarity
    ! polarity for 1997    Bsph_D(1) =  sign(SWH_Bx, zFace)             ! Br
    Bsph_D(1) =  sign(SWH_Bx, -zFace)             ! Br  !good for 2005 field
    Bsph_D(2) =  0.0                             ! Btheta
    Bsph_D(3) = -Bsph_D(1)*SinTheta*ParkerTilt   ! Bphi

    Vsph_D    = (/ SWH_Ux, 0.0, 0.0 /)           ! Vr, Vtheta, Vphi

    ! This is wrong:
    ! VphiSolar = OMEGAbody*(6.96E5)*sinTheta/unitUSER_U
    ! Vphi=omega*r*sinTheta - 6.96E5 is the solar radii in km
    ! No2Si_V(UnitU_) is in m/s so its ok
    ! No2Io_V(UnitU_) is in km/s
    ! factor 1000 for test
    ! --------------
    ! Correct:
    ! From angular momentum conservation 
    !    vPhi(rFace)*rFace*SinTheta = vPhi(rSun)*rSun*SinTheta
    ! and
    !    vPhi(rSun) = rSun*SinTheta*Omega
    ! so
    !    VphiSolar = OmegaSun*SinTheta*rSun**2/rFace
    ! Approximately VphiSolar = 0.0 (about 3km/s at 10AU)

    !! Introducing the Fast Solar Wind
    !! sin2Theta_fast_wind = 0.250000    ! 30 degrees
    !!  sin2Theta_fast_wind = 0.1786062   ! 25 degrees
    !!  sin2Theta_fast_wind = 0.116980    ! 20 degrees
    !!  sin2Theta_fast_wind = 1.000000    ! 90 degrees
    !!      if (sinTheta*sinTheta > sin2Theta_fast_wind) then
    !!       else
    !!          ! FAST WIND
    !!          VrSolarWind     =  SWfast_Ux
    !!          VthetaSolarWind =  0.0
    !!          VphiSolarWind   =  VphiSolar
    !!          BrSolarWind     =  BrSolar
    !!          BthetaSolarWind =  BthetaSolar
    !!          BphiSolarWind   =  BphiSolar
    !!          RhoSolarWind    =  SWfast_rho
    !!      end if
    !! Latitude variating wind
    !! VrSolarWind = SWH_Ux+((SWfast_Ux - SWH_Ux)*cosTheta*cosTheta)
    !! RhoSolarWind = SWH_rho + ((SWfast_rho - SWH_rho)*cosTheta*cosTheta)
    !! pressure will vary for fast and slow wind         PSolarWind = SWH_p
    !/
    !!test      VrSolarWind = SWH_Ux

    ! Calculate pressure (equilibrium, but why? It is supersonic!)
    Pmag = sum(Bsph_D**2) / 2.0

    ! magnetic pressure at the equator (actually wrong, neglects Br=SWH_Bx)
    PmagEquator = (SWH_Bx*ParkerTilt)**2/2
    pSolarWind  = SWH_p + PmagEquator -  Pmag

    ! Apply boundary conditions for ions
    VarsGhostFace_V(Rho_)    = SWH_rho
    VarsGhostFace_V(p_)      = SWH_p !!! pSolarWind
    VarsGhostFace_V(Ux_:Uz_) = matmul(XyzSph_DD, Vsph_D)
    VarsGhostFace_V(Bx_:Bz_) = matmul(XyzSph_DD, Bsph_D)

    if(UseNeutralFluid)then

       ! NeuRho is PopI; NeuIIRho is PopII and NeuIIIRho is PopIII
       !
       ! Pop I is going through the inner BCs    

       ! soft boundary for Pop I-IV

       VarsGhostFace_V(NeuRho_:Ne4P_) = VarsTrueFace_V(NeuRho_:Ne4P_)

    endif

    if(DoTestMe)then
       write(*,*) NameSub,' FaceCoord=', FaceCoords_D
       write(*,*) NameSub,' i,j,kFace=', iFace, jFace, kFace
       write(*,*) NameSub,' SinTheta =', SinTheta
       write(*,*) NameSub,' Bsph_D   =', Bsph_D
       write(*,*) NameSub,' Vsph_D   =', Vsph_D
       write(*,*) NameSub,' Pmag, Peq=', Pmag, PmagEquator
       write(*,*) NameSub,' RhoGhost =', VarsGhostFace_V(Rho_)
       write(*,*) NameSub,' pGhost   =', VarsGhostFace_V(p_)
       write(*,*) NameSub,' Ughost   =', VarsGhostFace_V(Ux_:Uz_)
       write(*,*) NameSub,' Bghost   =', VarsGhostFace_V(Bx_:Bz_)
       if(UseNeutralFluid)then
          write(*,*) NameSub,'Pop1=',VarsGhostFace_V(NeuRho_:NeuP_)
          write(*,*) NameSub,'Pop2=',VarsGhostFace_V(Ne2Rho_:Ne2P_)
          write(*,*) NameSub,'Pop3=',VarsGhostFace_V(Ne3Rho_:Ne3P_)
       end if
    end if

  end subroutine user_set_face_boundary

  !-------------------------------------------------------------------
  subroutine user_normalization

    use ModConst, ONLY: cAU, cProtonMass
    use ModPhysics, ONLY: No2Si_V, UnitX_, UnitU_, UnitRho_

    character (len=*), parameter :: NameSub='user_normalization'
    logical :: DoTest, DoTestMe
    !-------------------------------------------------------------------

    call set_oktest(NameSub, DoTest, DoTestMe)

    No2Si_V(UnitX_)  = cAU                                          ! m
    No2Si_V(UnitU_)  = sqrt(Gamma*cBoltzmann*SWH_T_dim/cProtonMass) ! m/s
    No2Si_V(UnitRho_)= cProtonMass*SWH_rho_dim*1.0E+6               ! kg/m^3

    if(DoTestMe)then
       write(*,*)NameSub,' No2Si_V(UnitX_)  =',No2Si_V(UnitX_)
       write(*,*)NameSub,' No2Si_V(UnitU_)  =',No2Si_V(UnitU_)
       write(*,*)NameSub,' No2Si_V(UnitRho_)=',No2Si_V(UnitRho_)
    end if

  end subroutine user_normalization
  !-------------------------------------------------------------------

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    ! The ISM enters at the east boundary (negative x)

    integer,intent(in)::iBlock, iSide
    character (len=*),intent(in) :: TypeBc
    logical,intent(out) :: IsFound

    integer :: i, iVar
    !----------------------------------------------------------------------
    IsFound = .true.

    State_VGB(rho_,-1:2,:,:,iBlock)=VLISW_rho
    State_VGB(rhoUx_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Ux
    State_VGB(rhoUy_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Uy
    State_VGB(rhoUz_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Uz
    State_VGB(Bx_,-1:2,:,:,iBlock)=VLISW_Bx
    State_VGB(By_,-1:2,:,:,iBlock)=VLISW_By
    State_VGB(Bz_,-1:2,:,:,iBlock)=VLISW_Bz
    State_VGB(p_,-1:2,:,:,iBlock)=VLISW_p

    if(UseNeutralFluid)then
       !
       ! PopIV is the one one coming with the ISW
       ! The separation between Pop IV and Pop I is arbitrary so
       ! we took the separation as Vlad in x=-1500AU

       State_VGB(Ne4Rho_,-1:2,:,:,iBlock)   = RhoNeutralsISW
       State_VGB(Ne4RhoUx_,-1:2,:,:,iBlock) = RhoNeutralsISW*UxNeutralsISW
       State_VGB(Ne4RhoUy_,-1:2,:,:,iBlock) = RhoNeutralsISW*UyNeutralsISW
       State_VGB(Ne4RhoUz_,-1:2,:,:,iBlock) = RhoNeutralsISW*UzNeutralsISW
       State_VGB(Ne4P_,-1:2,:,:,iBlock)     = PNeutralsISW

       !
       ! In general you should specify as many values as many incoming 
       ! characteristic waves are present. For a neutral fluid this 
       ! is 0 for supersonic outflow, 1 for subsonic outflow, 
       ! 4 for subsonic inflow and 5 for supersonic inflow. 

       !ausg26 cond
       !!State_VGB(NeuRho_,-1:2,:,:,iBlock)   = 0.639*RhoNeutralsISW    
       !!State_VGB(NeuRhoUx_,-1:2,:,:,iBlock) = 0.639*RhoNeutralsISW*UxNeutralsISW    
       !!State_VGB(NeuRhoUy_,-1:2,:,:,iBlock) = 0.639*RhoNeutralsISW*UyNeutralsISW   
       !!State_VGB(NeuRhoUz_,-1:2,:,:,iBlock) = 0.639*RhoNeutralsISW*UzNeutralsISW    
       ! hydro run State_VGB(NeuP_,-1:2,:,:,iBlock)     = 0.6798*PNeutralsISW
       !!State_VGB(NeuP_,-1:2,:,:,iBlock)     = 0.639*PNeutralsISW

       !
       !\
       ! PopII and III supersonic outflow
       !/
       !!do iVar = Ne2Rho_, Ne3P_; do i = -1, 2
       !conditions of August 26
       do iVar = NeuRho_, Ne3P_; do i = -1, 2
          State_VGB(iVar,i,:,:,iBlock) = State_VGB(iVar,3,:,:,iBlock)
       end do; end do
    endif

  end subroutine user_set_cell_boundary

  !=====================================================================

  subroutine user_set_ics(iBlock)

    use ModPhysics,  ONLY: rBody
    use ModCoordTransform, ONLY: rot_xyz_sph, rot_matrix_y, rot_matrix_z

    integer, intent(in) :: iBlock

    integer :: i,j,k

    real :: x, y, z, r
    real :: b_D(3), v_D(3), bSph_D(3), vSph_D(3)
    real :: SinTheta, SignZ

    real :: XyzSph_DD(3,3) ! rotation matrix Xyz_D = matmul(XyzSph_DD, Sph_D)

    ! These variables are not used now
    ! real :: thetaN, sinthetaN, lambda, RhoSolarW
    ! real :: sin2Theta_fast_wind

    ! tilted declarations
    ! Magnetic field and velocity vectors in spherical coordinates
    real :: bParker_D(3), vSolar_D(3)

    ! Rotation matrix between solar (solar means tilted coordinates) and simulation coordinates
    real, save :: RotMatrix_DD(3,3)

    ! I included an extra matrix to rotate between tilted and intertial coordinates
    real, save :: RotMatrixInv_DD(3,3)   !teste

    ! Position vector in solar coordinates (merav) at the grid point 
    real, dimension(3) :: Xyz_D


    ! The last time the rotation matrix was calculated
    real :: TimeLast = -1.0

    ! Initialize only once
    logical :: IsUninitialized = .true.

    ! Tilt of the rotation axis in degrees
    real :: ThetaTiltDeg, THETAtilt, rot_period_dim

    ! added in June 2010
    real :: xt, yt, zt, rt

    ! end of tilt declarations

    character(len=*), parameter:: NameSub = 'user_set_ics'
    logical :: DoTest, DoTestMe, DoTestCell
    !--------------------------------------------------------------------------

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    ! Make sure that OmegaSun and ParkerTilt are set
    if(OmegaSun == 0.0)call set_omega_parker_tilt

    if(IsUninitialized)then
       ! Initialize some constants

       THETAtiltDeg = 0.0
       ! Convert to radians :
       THETAtilt = cTwoPi*THETAtiltDeg/360.00

       ! Defining the rotation components of the Sun
       rot_period_dim = 26.0*24.0                   ! rotation period in hours
       OmegaBody = cTwoPi/(rot_period_dim*3600.00)

       IsUninitialized = .false.
    endif
    if(Time_Simulation /= TimeLast)then
       !      Calculate the rotation matrix for the current time
       !      The order and the signs should be carefully checked and verified

       THETAtiltDeg = 0.0
       ! Convert to radians :
       THETAtilt = cTwoPi*THETAtiltDeg/360.00

       ! Defining the rotation components of the Sun
       rot_period_dim = 26.0*24.0                   ! rotation period in hours
       OmegaBody = cTwoPi/(rot_period_dim*3600.00)


       ! to perform the negative rotation from TILTED -> INERTIAL

       !\
       ! Merav's Original Code
       !/
       !RotMatrix_DD = matmul(&
       !     rot_matrix_y(THETAtilt), &        ! rotate around y axis by tilt
       !     rot_matrix_z(OMEGAbody*Time_Simulation))  ! rotate around z axis by rotation
       !
       !RotMatrixInv_DD = matmul(&
       !     rot_matrix_y(-THETAtilt), &        ! rotate around y axis by tilt
       !     rot_matrix_z(-OmegaBody*Time_Simulation))
       !\
       ! End Merav's Original Code
       !/

       ! Gabor's  Edits
       RotMatrix_DD = matmul(&
            rot_matrix_y(THETAtilt), &        ! rotate around y axis by tilt
            rot_matrix_z(OMEGAbody*Time_Simulation))  ! rotate around z axis by rotation

       RotMatrixInv_DD = matmul(&
            rot_matrix_z(-OmegaBody*Time_Simulation), &
            rot_matrix_y(-THETAtilt))        ! rotate around y axis by tilt
       ! end Gabor's edits

       TimeLast = Time_Simulation
    end if

    do i=MinI,MaxI; do j=MinJ,MaxJ; do k=MinK,MaxK

       if(.not. true_cell(i,j,k,iBlock)) CYCLE

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)
       r = r_BLK(i,j,k,iBlock)

       !june 2010
       Xyz_D = matmul( (/x, y, z/), RotMatrix_DD)

       xt = Xyz_D(1)
       yt = Xyz_D(2)
       zt = Xyz_D(3) 
       rt = sqrt(xt**2+yt**2+zt**2)

       XyzSph_DD = rot_xyz_sph(x,y,z)

       ! theta is the angle measured from the pole
       SinTheta = sqrt(xt**2+yt**2)/rt

       !2015
       Vsph_D = (/ SWH_Ux, 0., 0. /)

       vSolar_D(1) = vSph_D(1) !Radial velocity
       vSolar_D(2) = vSph_D(2)
       vSolar_D(3) = vSph_D(3) !VERIFY THAT!!

       !calculating the Parker B field spherical components Bsph_D

       !2015
       rBody = 30.
       OmegaSun   = cTwoPi/(26.0*24.0*3600.00*Si2No_V(UnitT_))
       ParkerTilt = OmegaSun*rBody/SWH_Ux

       ! good for polarity of 1997       SignZ = sign(1.0, z)
       SignZ = -sign(1.0,z)  ! good for 2005

       ! dipole configuration
       !Bsph_D(1) = SignZ*SWH_Bx*(rBody/rt)**2  ! Br
       !Bsph_D(2) = 0.0                        ! Btheta
       !Bsph_D(3) = -SignZ*SWH_Bx*SinTheta*ParkerTilt*(rBody/rt) !Bphi

       !monopole magnetic field
       Bsph_D(1) = SWH_Bx*(rBody/rt)**2  ! Br
       Bsph_D(2) = 0.0                        ! Btheta
       Bsph_D(3) = SWH_Bx*SinTheta*ParkerTilt*(rBody/rt) !Bphi

       ! tilted
       bParker_D(1) = Bsph_D(1) !Radial component
       bParker_D(2) = Bsph_D(2) !Theta component
       bParker_D(3) = Bsph_D(3) !Phi component

       bParker_D    = matmul(RotMatrixInv_DD, bParker_D)

       ! magnetic field components in cartesian coordinates
       b_D = matmul(XyzSph_DD, bParker_D)

       State_VGB(Bx_:Bz_,i,j,k,iBlock) = b_D

       !velocity components in cartesian coordinates
       v_D = matmul(XyzSph_DD, Vsph_D)

       ! density and pressure
       State_VGB(Rho_,i,j,k,iBlock) = SWH_rho * (rBody/r)**2
       State_VGB(P_,i,j,k,iBlock)   = SWH_p   * (rBody/r)**(2*Gamma)

       ! momentum
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock)*v_D

       if(UseNeutralFluid)then
          !\
          ! PopI
          !/
          State_VGB(NeuRho_,i,j,k,iBlock)  =  RhoNeutralsISW
          State_VGB(NeuP_,i,j,k,iBlock) =   PNeutralsISW
          State_VGB(NeuRhoUx_,i,j,k,iBlock)=RhoNeutralsISW*UxNeutralsISW
          State_VGB(NeuRhoUy_,i,j,k,iBlock)=RhoNeutralsISW*UyNeutralsISW
          State_VGB(NeuRhoUz_,i,j,k,iBlock)=RhoNeutralsISW*UzNeutralsISW


          !! State_VGB(NeuRho_,i,j,k,iBlock) = &
          !!    RhoNeutralsISW*Exp(-lambda*thetaN/((r**2)*SinThetaN)) &
          !!    *Exp(-lambda*thetaN/(r*SinThetaN))              
          !! State_VGB(NeuP_,i,j,k,iBlock) = &
          !!    PNeutralsISW*Exp(-lambda*thetaN/(r*SinThetaN))

          !\
          ! PopII - I set it to be about 100km/s radially. 
          ! The temperature is set to be 10^5K for this population.
          ! v_D is the plasma velocity, we take one quarter of that.
          !/
          State_VGB(Ne2Rho_,i,j,k,iBlock) = 1.e-3 
          State_VGB(Ne2P_,i,j,k,iBlock)   = 1.e-3 
          State_VGB(Ne2RhoUx_:Ne2RhoUz_,i,j,k,iBlock) = 1.e-3

          !\ 
          ! PopIII
          !/
          State_VGB(Ne3Rho_,i,j,k,iBlock) = 0.1 * State_VGB(Rho_,i,j,k,iBlock)
          State_VGB(Ne3P_,i,j,k,iBlock)   = 0.1 * State_VGB(p_,i,j,k,iBlock)
          State_VGB(Ne3RhoUx_:Ne3RhoUz_,i,j,k,iBlock) = &
               State_VGB(Ne3Rho_,i,j,k,iBlock)*v_D

          !\
          ! PopIV
          !/
          State_VGB(Ne4Rho_,i,j,k,iBlock)  =RhoNeutralsISW
          State_VGB(Ne4RhoUx_,i,j,k,iBlock)=RhoNeutralsISW*UxNeutralsISW 
          State_VGB(Ne4RhoUy_,i,j,k,iBlock)=RhoNeutralsISW*UyNeutralsISW 
          State_VGB(Ne4RhoUz_,i,j,k,iBlock)=RhoNeutralsISW*UzNeutralsISW
          State_VGB(Ne4P_,i,j,k,iBlock) = PNeutralsISW
       endif

       if(DoTestCell)then
          write(*,*)NameSub,' x, y, z, r             =', x, y, z, r
          write(*,*)NameSub,' SignZ, SWH_Bx, SinTheta=',SignZ, SWH_Bx, SinTheta
          write(*,*)NameSub,' OmegaSun, rBody, SWH_Ux=',OmegaSun,rBody,SWH_Ux
          write(*,*)NameSub,' Vsph_D                 =',Vsph_D
          write(*,*)NameSub,' v_D                    =',v_D
          write(*,*)NameSub,' Bsph_D                 =',Bsph_D
          write(*,*)NameSub,' b_D    =',State_VGB(Bx_:Bz_,i,j,k,iBlock)
          write(*,*)NameSub,' RhoU_D =',State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)
          write(*,*)NameSub,' rho,   =',State_VGB(Rho_,i,j,k,iBlock)
          write(*,*)NameSub,' p      =',State_VGB(P_,i,j,k,iBlock)
       end if

    end do; end do; end do

  end subroutine user_set_ics

  !=====================================================================
  subroutine user_initial_perturbation

    use ModEnergy, only: calc_energy

    ! Set Pop II to be the same state as the ions with smaller density
    integer:: iBlock

    character (len=*), parameter :: NameSub = 'user_initial_perturbation'
    !-------------------------------------------------------------------

    !This subroutine is not needed when not using the 4 neutral fluids
    if(.not.UseNeutralFluid) &
         call CON_stop(NameSub//':  no neutral fluids present')

    do iBlock = 1, nBlock

       if( Unused_B(iBlock) ) CYCLE

       State_VGB(Ne2Rho_,:,:,:,iBlock) = RhoBcFactor_I(Ne2_) * &
            State_VGB(Rho_,:,:,:,iBlock)

       State_VGB(Ne2RhoUx_:Ne2RhoUz_,:,:,:,iBlock) = &
            RhoBcFactor_I(Ne2_)*uBcFactor_I(Ne2_)* &
            State_VGB(RhoUx_:RhoUz_,:,:,:,iBlock)

       State_VGB(Ne2P_,:,:,:,iBlock) = RhoBcFactor_I(Ne2_) * &
            State_VGB(p_,:,:,:,iBlock)

       call calc_energy(MinI,MaxI, MinJ,MaxJ, MinK,MaxK, iBlock, Ne2_, Ne2_)
    end do

  end subroutine user_initial_perturbation

  !=====================================================================

  subroutine user_update_states(iBlock)

    use ModAdvance, ONLY: StateOld_VCB, State_VGB, EnergyOld_CBI, Energy_GBI
    use ModGeometry, ONLY: rMin_BLK

    integer,intent(in):: iBlock

    integer:: i, j, k

    !------------------------------------------------------------------
    call update_states_MHD(iBlock)

    ! No need to check blocks outside:
    if(rMin_BLK(iBlock) > rBody) RETURN

    do k=1, nK; do j=1, nJ; do i=1, nI
       if(r_BLK(i,j,k,iBlock)  > rBody) CYCLE

       if(iTableSolarwind > 0)then
          ! Calculate the time dependent solar wind
          call calc_time_dep_sw(i,j,k,iBlock)
       else
          ! Retain initial condition if time dependent solar wind is not used
          State_VGB(Rho_:p_,i,j,k,iBlock) = StateOld_VCB(Rho_:p_,i,j,k,iBlock) 
          Energy_GBI(i,j,k,iBlock,1) = EnergyOld_CBI(i,j,k,iBlock,1)
       endif

    end do; end do; end do

  contains
    !======================================================================
    subroutine calc_time_dep_sw(i,j,k,iBlock)

      use BATL_lib,       ONLY: Xyz_DGB
      use ModCoordTransform, ONLY: rot_xyz_sph
      use ModLookupTable,    ONLY: interpolate_lookup_table
      use ModEnergy,      ONLY: calc_energy

      integer,intent(in):: i, j, k, iBlock

      ! variables for Solar Cycle
      real :: Rho, Ur, Temp, p, x, y, z, r, Latitude
      real :: Bsph_D(3), Vsph_D(3)

      real :: XyzSph_DD(3,3) ! rotation matrix Xyz_D = matmul(XyzSph_DD,Sph_D)

      real, parameter:: LengthCycle = 662206313.647 ! length of solar cycle

      real :: TimeCycle ! holds current time of the simulation
      real :: Value_I(3)
      real :: SinTheta
      !---------------------------------------------------------------------

      ! calculate the latitude of the cell
      x = Xyz_DGB(1,i,j,k,iBlock)
      y = Xyz_DGB(2,i,j,k,iBlock)
      z = Xyz_DGB(3,i,j,k,iBlock)
      r = r_BLK(i,j,k,iBlock)

      XyzSph_DD = rot_xyz_sph(x,y,z)

      SinTheta = sqrt(x**2+y**2)/r

      ! calculating latitude of the cell
      Latitude = cRadToDeg*asin(z/r)

      ! calculating time relative to the solar cycle
      TimeCycle = modulo(Time_Simulation, LengthCycle) 

      ! interpolating the value of Rho, Vr, and Temp 
      ! at the cell from the lookup table
      call interpolate_lookup_table(iTableSolarwind, Latitude, TimeCycle, &
           Value_I)

      Ur  = Value_I(1)*Io2No_V(UnitU_)
      Rho = Value_I(2)*Io2No_V(UnitRho_)
      Temp= Value_I(3)*Io2No_V(UnitTemperature_)
      p = 2.0*Rho*Temp

      ! Spherical velocity, Vr, Vtheta, Vphi constant with  radial distance
      Vsph_D    = (/ Ur, 0.0, 0.0 /)         

      !\
      !monopole with By negative and a time varying B 
      ! time-dependent behavior of B taken from Michael et al. 2015
      !/  
      Bsph_D(1) = (SQRT(0.5)/rBody**2)*(9.27638+ &
           7.60832d-8*TimeCycle-1.91555*SIN(1.28737d-8*TimeCycle)+ &
           0.144184*SIN(2.22823d-8*TimeCycle)+ &
           47.7758*SIN(2.18788d-10*TimeCycle)+ &
           83.5522*SIN(-1.20266d-9*TimeCycle))*Io2No_V(UnitB_) ! Br
      Bsph_D(2) =  0.0                             ! Btheta
      Bsph_D(3) = Bsph_D(1)*SinTheta*ParkerTilt*SWH_Ux/Ur ! Bphi for vary B
      !Bsph_D(3) = Bsph_D(1)*SinTheta*ParkerTilt   ! Bphi for B independent of Vsw

      !\
      !for dipole with a time varying B
      !/
      !Bsph_D(1) = sign((SQRT(0.5)/rBody**2)*(9.27638+ &
      !     7.60832d-8*TimeCycle-1.91555*SIN(1.28737d-8*TimeCycle)+ &
      !     0.144184*SIN(2.22823d-8*TimeCycle)+ &
      !     47.7758*SIN(2.18788d-10*TimeCycle)+ &
      !     83.5522*SIN(-1.20266d-9*TimeCycle))*Io2No_V(UnitB_), z) ! Br
      !Bsph_D(2) =  0.0                             ! Btheta
      !Bsph_D(3) = -Bsph_D(1)*SinTheta*ParkerTilt   ! Bphi for B independent of Vsw
      !Bsph_D(3) = -Bsph_D(1)*SinTheta*ParkerTilt*SWH_Ux/Ur ! Bphi for vary B

      ! Scale density, pressure, and magnetic field with radial distance
      Rho = Rho*(rBody/r)**2
      p     = p*(rBody/r)**(2*Gamma)
      Bsph_D(1) = Bsph_D(1)*(rBody/r)**2
      Bsph_D(3) = Bsph_D(3)*(rBody/r)

      ! Setting the state variables
      State_VGB(Rho_,i,j,k,iBlock) = Rho

      ! Velocity converted to 3 components of momentum in Cartesian coords
      State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = matmul(XyzSph_DD, Rho*Vsph_D)

      ! Spherical magnetic field converted to Cartesian components
      State_VGB(Bx_:Bz_,i,j,k,iBlock) = matmul(XyzSph_DD, Bsph_D)

      ! Sets pressure and energy only for the ion fluid
      State_VGB(p_,i,j,k,iBlock) = p
      call calc_energy(i, i, j, j, k, k, iBlock, 1, 1)

    end subroutine calc_time_dep_sw

  end subroutine user_update_states

  !=====================================================================

  subroutine user_action(NameAction)

    character(len=*), intent(in):: NameAction

    character(len=*), parameter:: StringFormat = '(10X,A19,F15.6,A11,F15.6)'
    !-----------------------------------------------------------------------

    if(NameAction /= 'write progress') RETURN

    write(*,StringFormat) 'SWH_rho_dim [n/cc]:',SWH_rho_dim,'SWH_rho:',SWH_rho

    write(*,StringFormat) 'SWH_Ux_dim  [km/s]:',SWH_Ux_dim,'SWH_Ux:',SWH_Ux
    write(*,StringFormat) 'SWH_Uy_dim  [km/s]:',SWH_Uy_dim,'SWH_Uy:',SWH_Uy 
    write(*,StringFormat) 'SWH_Uz_dim  [km/s]:',SWH_Uz_dim,'SWH_Uz:',SWH_Uz
    write(*,StringFormat) 'SWH_T_dim   [   K]:',SWH_T_dim,'SWH_p:',SWH_p
    write(*,StringFormat) 'SWH_Bx_dim  [  nT]:',SWH_Bx_dim,'SWH_Bx:',SWH_Bx 
    write(*,StringFormat) 'SWH_By_dim  [  nT]:',SWH_By_dim,'SWH_By:',SWH_By
    write(*,StringFormat) 'SWH_Bz_dim  [  nT]:',SWH_Bz_dim,'SWH_Bz:',SWH_Bz
    write(*,'(10X,A19,F15.6)')           'SWH_T_dim   [   K]:',SWH_T_dim
    ! fast solar wind
    write(*,*)
    write(*,*)
    write(*,StringFormat) 'VLISW_rho_dim[n/cc]:',VLISW_rho_dim,'VLISW_rho:',VLISW_rho 
    write(*,StringFormat) 'VLISW_Ux_dim[km/s]: ',VLISW_Ux_dim,'VLISW_Ux:',VLISW_Ux
    write(*,StringFormat) 'VLISW_Uy_dim[km/s]: ',VLISW_Uy_dim,'VLISW_Uy:',VLISW_Uy
    write(*,StringFormat) 'VLISW_Uz_dim[km/s]: ',VLISW_Uz_dim,'VLISW_Uz:',VLISW_Uz 
    write(*,StringFormat) 'VLISW_p_dim [nPa]: ',VLISW_p_dim,'VLISW_p:',VLISW_p 
    write(*,StringFormat) 'VLISW_Bx_dim[nT]: ',VLISW_Bx_dim,'VLISW_Bx:',VLISW_Bx
    write(*,StringFormat) 'VLISW_By_dim[nT]:',VLISW_By_dim,'VLISW_By:',VLISW_By
    write(*,StringFormat) 'VLISW_Bz_dim[nT]:',VLISW_Bz_dim,'VLISW_Bz:',VLISW_Bz
    write(*,'(10X,A19,F15.6)') 'VLISW_a_dim[km/s]: ',VLISW_a_dim
    write(*,'(10X,A19,F15.6)') 'VLISW_T_dim[K]: ',VLISW_T_dim! 
    if(UseNeutralFluid)then
       !neutrals
       write(*,*)     
       write(*,StringFormat) 'RhoNeutralsISW_dim:',RhoNeutralsISW_dim ,'RhoNeutralsISW:',RhoNeutralsISW 
       write(*,StringFormat) 'UxNeutralsISW_dim:',UxNeutralsISW_dim,'UxNeutralsISW:',UxNeutralsISW 
       write(*,StringFormat) 'UyNeutralsISW_dim:',UyNeutralsISW_dim,'UyNeutralsISW:',UyNeutralsISW
       write(*,StringFormat) 'UzNeutralsISW_dim:',UzNeutralsISW_dim,'UzNeutralsISW:',UzNeutralsISW 
       write(*,StringFormat) 'PNeutralsISW_dim:',PNeutralsISW_dim,'PNeutralsISW:',PNeutralsISW
       write(*,'(10X,A19,F15.6)') 'TNeutralsISW_dim:',TNeutralsISW_dim     
    endif
    write(*,*)

  end subroutine user_action

  !=====================================================================
  subroutine user_io_units

    use ModConst, ONLY: cAU, cProtonMass
    character (len=*), parameter :: NameSub='user_io_units'
    !-------------------------------------------------------------------------
    ! set strings for writing Tecplot output
    !/
    NameTecUnit_V(UnitX_)            = 'AU'
    NameTecUnit_V(UnitRho_)          = 'amu/cm3'
    NameTecUnit_V(UnitN_)            = '#/cm3'
    NameTecUnit_V(UnitU_)            = 'km/s'
    NameTecUnit_V(UnitP_)            = 'dyne/cm^2'
    NameTecUnit_V(UnitB_)            = 'nT'
    NameTecUnit_V(UnitRhoU_)         = 'g/cm^2/s'
    NameTecUnit_V(UnitEnergydens_)   = 'erg/cm^3'
    NameTecUnit_V(UnitJ_)            = 'uA/m^2'
    NameTecUnit_V(UnitDivB_)         = 'G/cm'
    NameTecUnit_V(UnitAngle_)        = 'deg'

    NameIdlUnit_V(UnitX_)            = 'AU'
    NameIdlUnit_V(UnitRho_)          = 'amu/cc'
    NameIdlUnit_V(UnitN_)            = '/cc'
    NameIdlUnit_V(UnitU_)            = 'km/s'
    NameIdlUnit_V(UnitP_)            = 'dyne/cm^2'
    NameIdlUnit_V(UnitB_)            = 'nT'
    NameIdlUnit_V(UnitRhoU_)         = 'g/cm^2/s'
    NameIdlUnit_V(UnitEnergydens_)   = 'erg/cm^3'
    NameIdlUnit_V(UnitJ_)            = 'uA/m^2'
    NameIdlUnit_V(UnitDivB_)         = 'G/cm'
    NameIdlUnit_V(UnitAngle_)        = 'deg'
 
    Io2Si_V(UnitX_)           = cAU                       ! R  
    Io2Si_V(UnitRho_)         = 1.0E+6*cProtonMass        ! Mp/cm^3
    Io2Si_V(UnitN_)           = 1.0E+6                    ! #/cm^3
    Io2Si_V(UnitU_)           = 1.0E+3                    ! km/s
    Io2Si_V(UnitP_)           = 1.0E-1                    ! dyne/cm^2
    Io2Si_V(UnitB_)           = 1.0E-9                    ! nT
    Io2Si_V(UnitRhoU_)        = 1.0E+1                    ! g/cm^2/s
    Io2Si_V(UnitEnergydens_)  = 1.0E-1                    ! erg/cm^3
    Io2Si_V(UnitJ_)           = 1.0E-6                    ! uA/m^2
    Io2Si_V(UnitDivB_)        = 1.0E-2                    ! Gauss/cm
    Io2Si_V(UnitAngle_)       = cRadToDeg                 ! degrees

    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    !  normalization of SWH and VLISW and Neutrals

    VLISW_a_dim    = No2Io_V(UnitU_)*(VLISW_T_dim/SWH_T_dim)
    VLISW_p_dim1    = No2Io_V(UnitP_)/Gamma &
         *(VLISW_rho_dim/SWH_rho_dim)*(VLISW_T_dim/SWH_T_dim)
    ! 
    ! Pressure of plasma = 2*T_ion*rho_ion

    VLISW_B_factor = No2Io_V(UnitB_)*sqrt((VLISW_T_dim/SWH_T_dim) &
         *(VLISW_rho_dim/SWH_rho_dim))

    VLISW_rho = VLISW_rho_dim*Io2No_V(UnitRho_)
    VLISW_p1   = VLISW_p_dim1*Io2No_V(UnitP_)
    VLISW_p    = 2.*VLISW_T_dim*Io2No_V(UnitTemperature_)*VLISW_rho

    !merav
    !write(*,*) 'VLISW_p1',VLISW_p1
    !write(*,*) 'VLISW_p',VLISW_p
    !merav

    VLISW_Ux  = VLISW_Ux_dim*Io2No_V(UnitU_)
    VLISW_Uy  = VLISW_Uy_dim*Io2No_V(UnitU_)
    VLISW_Uz  = VLISW_Uz_dim*Io2No_V(UnitU_)
    VLISW_Bx  = VLISW_Bx_dim*Io2No_V(UnitB_)
    VLISW_By  = VLISW_By_dim*Io2No_V(UnitB_)
    VLISW_Bz  = VLISW_Bz_dim*Io2No_V(UnitB_)

    ! Latitude Dependent Wind
    SWfast_rho = SWfast_rho_dim*Io2No_V(UnitRho_)
    SWfast_p   = SWfast_p_dim*Io2No_V(UnitP_) 
    SWfast_Ux  = SWfast_Ux_dim*Io2No_V(UnitU_)
    SWfast_Uy  = SWfast_Uy_dim*Io2No_V(UnitU_)
    SWfast_Uz  = SWfast_Uz_dim*Io2No_V(UnitU_)

    SWH_rho = SWH_rho_dim*Io2No_V(UnitRho_)

    ! Pressure of plasma = 2*T_ion*rho_ion
    SWH_p   = 2.*SWH_T_dim*Io2No_V(UnitTemperature_)*SWH_rho

    !merav
    !write(*,*) 'SWH_p1',SWH_p1
    !write(*,*) 'SWH_p',SWH_p
    !merav


    SWH_Ux  = SWH_Ux_dim*Io2No_V(UnitU_)
    SWH_Uy  = SWH_Uy_dim*Io2No_V(UnitU_)
    SWH_Uz  = SWH_Uz_dim*Io2No_V(UnitU_)
    SWH_Bx  = SWH_Bx_dim*Io2No_V(UnitB_)
    SWH_By  = SWH_By_dim*Io2No_V(UnitB_)
    SWH_Bz  = SWH_Bz_dim*Io2No_V(UnitB_)

    SWfast_p_dim = No2Io_V(UnitP_)/Gamma*(SWfast_rho_dim/SWH_rho_dim)
    SWfast_p = SWfast_p_dim*Io2No_V(UnitP_)
    !
    !
    ! The units of rho_dim are n/cc and unitUSER_rho Gamma/cc

    if(UseNeutralFluid)then
       !/
       !merav june01    PNeutralsISW   = RhoNeutralsISW * 
       !TNeutralsISW_dim*Io2No_V(UnitTemperature_)
       !merav june01  SWH_p   = SWH_rho * SW_T_dim*Io2No_V(UnitTemperature_)

       RhoNeutralsISW = RhoNeutralsISW_dim*Io2No_V(UnitRho_)
       PNeutralsISW_dim = No2Io_V(UnitP_)/Gamma*(RhoNeutralsISW_dim/SWH_rho_dim)*(TNeutralsISW_dim /SWH_T_dim)

       !PNeutralsISW1 = PNeutralsISW_dim*Io2No_V(UnitP_)
       PNeutralsISW = TNeutralsISW_dim*Io2No_V(UnitTemperature_)*RhoNeutralsISW 
       UxNeutralsISW  = UxNeutralsISW_dim*Io2No_V(UnitU_)
       UyNeutralsISW  = UyNeutralsISW_dim*Io2No_V(UnitU_)
       UzNeutralsISW  = UzNeutralsISW_dim*Io2No_V(UnitU_)
       mNeutrals    = mProtonMass*cProtonMass

       !merav
       !write(*,*) 'PNeutralsISW',PNeutralsISW
       !write(*,*) 'PNeutralsISW1',PNeutralsISW1
       !merav
    endif

  end subroutine user_io_units
  !==============================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK

    integer,          intent(in)   :: iBlock
    character(len=*), intent(in)   :: NameVar
    logical,          intent(in)   :: IsDimensional
    real,             intent(out)  :: PlotVar_G(MinI:MaxI, MinJ:MaxJ, MinK:MaxK)
    real,             intent(out)  :: PlotVarBody
    logical,          intent(out)  :: UsePlotVarBody
    character(len=*), intent(inout):: NameTecVar
    character(len=*), intent(inout):: NameTecUnit
    character(len=*), intent(inout):: NameIdlUnit
    logical,          intent(out)  :: IsFound

    character (len=*), parameter :: NameSub='user_set_plot_var'

    !-------------------------------------------------------------------

    UsePlotVarBody = .true.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('srho')
       NameTecVar = 'Srho'
       PlotVar_G(1:nI,1:nJ,1:nK) = Source_VC(NeuRho_,:,:,:)
    case('fluid')
       if(UseNeutralFluid)then
          call select_region(iBlock)
          PlotVar_G(1:nI,1:nJ,1:nK) = iFluidProduced_C
       else
          !This is not needed when not using the 4 nutral fluids
          call CON_stop(NameSub//': no neutral fluids present')
       endif
    case('mach')
       PlotVar_G = &
            sqrt( sum(State_VGB(RhoUx_:RhoUz_,:,:,:,iBlock)**2, DIM=1)      &
            /    (Gamma *State_VGB(p_,:,:,:,iBlock)*State_VGB(Rho_,:,:,:,iBlock)) )
       !merav addition
    case('machalfven')
       PlotVar_G = & 
            sqrt( sum(State_VGB(RhoUx_:RhoUz_,:,:,:,iBlock)**2, DIM=1)      &
            /   ( sum(State_VGB(Bx_:Bz_,:,:,:,iBlock)**2, DIM=1) &
            * State_VGB(Rho_,:,:,:,iBlock)) )

       !end of merav addition         
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !====================================================================

  subroutine user_calc_sources(iBlock)
    !\
    ! Calculates the charge exchange cross sections for the neutrals
    ! This subroutine calculate the charge exchange between the ionized 
    ! component and the three population of neutrals, Neu, Ne2, Ne3 
    ! following the notation of McNutt (1998)
    ! The three population of neutrals are the neutrals created by 
    ! charge exchange in different regions of the outer heliosphere
    !
    ! Neu are the neutrals the comes from the interstellar medium  (PopI)
    ! Ne2 are created between the Termination Shock and Heliopause (PopII)
    ! Ne3 are the neutrals created inside the Termination Shock    (PopIII)
    !
    ! As an example for Neu the source terms inside the Termination Shock will 
    ! take into account their creations and destruction; 
    ! outside the Termination Shock they will be just destroyed 
    ! (see more details below). 
    !
    ! The _I(1) is the ionized fluid and _I(Neu_)-_I(Ne4_) are the neutral fluids
    !/
    ! 
    ! History of implementation:
    ! Written by Merav Opher April 21, 2002;
    ! Help From Gabor Toth
    !/
    ! *modified for the version 7.5.2 May, 2002*
    ! *revistied and modified on Nov, 2002*
    ! source terms for the plasma
    ! modified to take of the S=0 Feb 08 (initialization) 
    ! modified June 08, 2007 to take into account multi-fluids
    ! modified August 18;September10, 2007 to take into account 3-fluids
    ! october 27 checking units 
    ! January 01, 2008 implementing implicit source terms 
    ! January 08, 2008 all source terms for PopI as implicit
    ! January 24, 2008 comments included
    ! January-June BUGS fixed
    ! September - source terms written as McNutt (1998)
    ! 
    !-------------------------------------------------------------------
    use ModPointImplicit, ONLY:  UsePointImplicit, &
         IsPointImplSource, IsPointImplPerturbed
    use ModNumConst

    integer, intent(in) :: iBlock


    real :: cth
    real :: State_V(nVar)  
    real :: Ux, Uy, Uz, U2

    real, dimension(nFluid) :: &
         Ux_I, Uy_I, Uz_I, U2_I, Temp_I, &
         UThS_I, URelS_I, URelSdim_I, UStar_I, Sigma_I, Rate_I, &
         UStarM_I, SigmaN_I, RateN_I, &
         I0xp_I, I0px_I, I2xp_I, I2px_I, &
         JxpUx_I, JxpUy_I, JxpUz_I, JpxUx_I, JpxUy_I, JpxUz_I, &
         Kxp_I, Kpx_I, Qepx_I, QmpxUx_I, QmpxUy_I, QmpxUz_I


    integer :: i, j, k

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'user_calc_sources'
    !-----------------------------------------------------------------------

    if(iBlock == BLKtest .and. iProc == PROCtest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    endif

    !updating sources from AMPS in PT-OH coupling
    if(.not.UseNeutralFluid)then
       if(allocated(ExtraSource_ICB))then

          do k = 1, nK; do j = 1, nJ; do i = 1, nI

             !skipping the cells inside rBody
             if(.not.true_cell(i,j,k,iBlock))CYCLE

             !Extract conservative variables
             State_V = State_VGB(:,i,j,k,iBlock)
             
             Ux  = State_V(RhoUx_)/State_V(Rho_)
             Uy  = State_V(RhoUy_)/State_V(Rho_)
             Uz  = State_V(RhoUz_)/State_V(Rho_)

             U2  = Ux**2 + Uy**2 + Uz**2

             !updating the source terms
             Source_VC(Rho_,i,j,k) = Source_VC(Rho_,i,j,k) &
                  + ExtraSource_ICB(1,i,j,k,iBlock)
             Source_VC(RhoUx_,i,j,k) = Source_VC(RhoUx_,i,j,k) &
                  + ExtraSource_ICB(2,i,j,k,iBlock)
             Source_VC(RhoUy_,i,j,k) = Source_VC(RhoUy_,i,j,k) &
                  + ExtraSource_ICB(3,i,j,k,iBlock)
             Source_VC(RhoUz_,i,j,k) = Source_VC(RhoUz_,i,j,k) &
                  + ExtraSource_ICB(4,i,j,k,iBlock)
             Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) &
                  + ExtraSource_ICB(5,i,j,k,iBlock)
             Source_VC(p_,i,j,k) = Source_VC(p_,i,j,k) &
                  + GammaMinus1*( Source_VC(Energy_,i,j,k) &
                  - Ux*Source_VC(RhoUx_,i,j,k) &
                  - Uy*Source_VC(RhoUy_,i,j,k) &
                  - Uz*Source_VC(RhoUz_,i,j,k) &
                  + 0.5*U2*Source_VC(Rho_,i,j,k) )

          end do; end do; end do
       end if

       RETURN
    end if

    ! Do not provide explicit source term when point-implicit scheme is used
    ! IsPointImplSource is true only when called from ModPointImplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    !  calculating some constants cBoltzmann is J/K 
    cth = 2.0*cBoltzmann/mNeutrals

    ! Figure out which neutral population is produced at this point
    if(.not.IsPointImplPerturbed) call select_region(iBlock) 

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! Extract conservative variables
       State_V = State_VGB(:, i, j, k, iBlock)

!!! Experiment: zero temperature Neu
!!!if(DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest) &
!!!     State_V(NeuP_) = 1e-30 !!!

       Ux_I  = State_V(iRhoUx_I)/State_V(iRho_I)
       Uy_I  = State_V(iRhoUy_I)/State_V(iRho_I)
       Uz_I  = State_V(iRhoUz_I)/State_V(iRho_I)

       ! Velocity square for the ionized and three population of neutrals; 
       U2_I  = Ux_I**2 + Uy_I**2 + Uz_I**2

       ! Temperature for the ionized and three population of neutrals (K)
       ! T = p/rho 
       ! since P_plasma=2T_proton*rho_ion; T_proton=0.5P_plasma/rho_ion

       Temp_I       = (State_V(iP_I)/State_V(iRho_I))*No2Si_V(UnitTemperature_)
       Temp_I(Ion_) = 0.5*Temp_I(Ion_)

       ! Thermal speed (squared) for ionized and three populations of neutrals
       ! UThS units are (m/s)^2
       UThS_I = cth*Temp_I 

       ! Relative velocity between neutrals and ionized fluid squared
       !! URelS_I = (Ux_I - Ux_I(1))**2 &
       !!     +    (Uy_I - Uy_I(1))**2 &
       !!     +    (Uz_I - Uz_I(1))**2 

       URelS_I(1) =0.

       URelS_I(Neu_) = (Ux_I(Neu_) - Ux_I(1))**2 &
            +    (Uy_I(Neu_) - Uy_I(1))**2 &
            +    (Uz_I(Neu_) - Uz_I(1))**2

       URelS_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(1))**2 &
            +    (Uy_I(Ne2_) - Uy_I(1))**2 &
            +    (Uz_I(Ne2_) - Uz_I(1))**2

       URelS_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(1))**2 &
            +    (Uy_I(Ne3_) - Uy_I(1))**2 &
            +    (Uz_I(Ne3_) - Uz_I(1))**2

       URelS_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(1))**2 &
            +    (Uy_I(Ne4_) - Uy_I(1))**2 &
            +    (Uz_I(Ne4_) - Uz_I(1))**2

       ! Calculating Cross Section Sigma_I for the different neutrals
       !
       ! Incorporating units to calculate the charge exchange cross sections
       ! No2Si_V(UnitU_) has units of m/s like cstartT so UReldim and UStar 
       ! has units of m/s

       !! URelSdim_I  = URelS_I * No2Si_V(UnitU_)**2

       URelSdim_I(1)=0.
       URelSdim_I(Neu_)  = URelS_I(Neu_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne2_)  = URelS_I(Ne2_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne3_)  = URelS_I(Ne3_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne4_)  = URelS_I(Ne4_) * No2Si_V(UnitU_)**2

       ! UStar_I has units of m/s

       !!  UStar_I  = sqrt(URelSdim_I + (4./cPi)*(UThS_I +UThS_I(1)))

       UStar_I(1)=0. 
       UStar_I(Neu_) =  sqrt(URelSdim_I(Neu_) + (4./cPi)*(UThS_I(Neu_) +UThS_I(1)))
       UStar_I(Ne2_) =  sqrt(URelSdim_I(Ne2_) + (4./cPi)*(UThS_I(Ne2_) +UThS_I(1)))
       UStar_I(Ne3_) =  sqrt(URelSdim_I(Ne3_) + (4./cPi)*(UThS_I(Ne3_) +UThS_I(1)))
       UStar_I(Ne4_) =  sqrt(URelSdim_I(Ne4_) + (4./cPi)*(UThS_I(Ne4_) +UThS_I(1)))

       ! UStarM_I has units of m/s

       !!   UStarM_I  = sqrt(URelSdim_I + (64./(9.*cPi))*(UThS_I +UThS_I(1)))

       UStarM_I(1)=1.0
       UStarM_I(Neu_)  = sqrt(URelSdim_I(Neu_) + (64./(9.*cPi))*(UThS_I(Neu_) +UThS_I(1)))
       UStarM_I(Ne2_)  = sqrt(URelSdim_I(Ne2_) + (64./(9.*cPi))*(UThS_I(Ne2_) +UThS_I(1)))
       UStarM_I(Ne3_)  = sqrt(URelSdim_I(Ne3_) + (64./(9.*cPi))*(UThS_I(Ne3_) +UThS_I(1)))
       UStarM_I(Ne4_)  = sqrt(URelSdim_I(Ne4_) + (64./(9.*cPi))*(UThS_I(Ne4_) +UThS_I(1)))


       ! Maher and Tinsley cross section Sigma 
       ! UStar has to have units of cm/s so the factor 100 is to pass m to cm
       ! Sigma has units of units of m^2

       !! Sigma_I =((1.64E-7 - (6.95E-9)*log(UStarM_I*100.))**2)*(1.E-4)
       !! SigmaN_I =((1.64E-7 - (6.95E-9)*log(UStar_I*100.))**2)*(1.E-4)       

       Sigma_I(1)=0.
       SigmaN_I(1)=0.

       Sigma_I(Neu_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Neu_)*100.))**2)*(1.E-4)
       SigmaN_I(Neu_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Neu_)*100.))**2)*(1.E-4)
       Sigma_I(Ne2_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne2_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne2_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne2_)*100.))**2)*(1.E-4)
       Sigma_I(Ne3_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne3_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne3_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne3_)*100.))**2)*(1.E-4)
       Sigma_I(Ne4_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne4_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne4_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne4_)*100.))**2)*(1.E-4)


       ! New Cross Section from Lindsay and Stebbings, 2005
       !!Sigma_I =((2.2835E-7 - (1.062E-8)*log(UStarM_I*100.))**2)*(1.E-4)

       !!SigmaN_I =((2.2835E-7 - (1.062E-8)*log(UStar_I*100.))**2)*(1.E-4)

       !!Sigma_I(Neu_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Neu_)*100.))**2)*(1.E-4)
       !!SigmaN_I(Neu_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Neu_)*100.))**2)*(1.E-4)
       !!Sigma_I(Ne2_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne2_)*100.))**2)*(1.E-4)
       !!SigmaN_I(Ne2_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne2_)*100.))**2)*(1.E-4)
       !!Sigma_I(Ne3_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne3_)*100.))**2)*(1.E-4)
       !!SigmaN_I(Ne3_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne3_)*100.))**2)*(1.E-4)
       !!Sigma_I(Ne4_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne4_)*100.))**2)*(1.E-4)
       !!SigmaN_I(Ne4_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne4_)*100.))**2)*(1.E-4)

       ! Calculating Rate  = \nu * nH * mp where nH is the density of neutrals
       ! \nu = Sigma*np*u_star where np is the density of the ionized flow and 
       ! For each population of neutrals there will be another Rate
       ! The charge exhange cross section 100 to change ustar to cm/s
       ! Rate has no units (m^2*m/s*s*m-3 )

       !! Rate_I =Sigma_I*State_V(Rho_)*State_V(iRho_I)*UStarM_I  &
       !!      *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       Rate_I(1)=0.

       Rate_I(Neu_) =Sigma_I(Neu_)*State_V(Rho_)*State_V(iRho_I(Neu_))*UStarM_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne2_) =Sigma_I(Ne2_)*State_V(Rho_)*State_V(iRho_I(Ne2_))*UStarM_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne3_) =Sigma_I(Ne3_)*State_V(Rho_)*State_V(iRho_I(Ne3_))*UStarM_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne4_) =Sigma_I(Ne4_)*State_V(Rho_)*State_V(iRho_I(Ne4_))*UStarM_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)


       !!RateN_I =SigmaN_I*State_V(Rho_)*State_V(iRho_I)*UStar_I  &
       !!     *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       RateN_I(1)=0.      

       RateN_I(Neu_) =SigmaN_I(Neu_)*State_V(Rho_)*State_V(iRho_I(Neu_))*UStar_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne2_) =SigmaN_I(Ne2_)*State_V(Rho_)*State_V(iRho_I(Ne2_))*UStar_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne3_) =SigmaN_I(Ne3_)*State_V(Rho_)*State_V(iRho_I(Ne3_))*UStar_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne4_) =SigmaN_I(Ne4_)*State_V(Rho_)*State_V(iRho_I(Ne4_))*UStar_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       ! Calculating the terms that enter in the Source terms
       ! The expressions for I0, Jxp, Kxp, Qexp are taken from Zank et al. 1996
       ! For example for population I; Neu:
       !
       ! Source(density) = I0xpNe2 + I0xpNe3 in region 1 
       !                 = - I0pxNeu         otherwise
       !
       ! Source(momentum) = QmxpNeu + JxpNe2 + JxpNe3 in region 1
       !                  = - JpxNeu                  otherwise
       !
       ! Source(energy) = QexpNeu + KxpNe2 + KxpNe3 in region 1
       !                = -KpxNeu                   otherwise
       ! 
       ! 'xp' indicates neutrals->protons charge exchange rates and
       ! 'px' indicates protons->neutrals
       ! charge exchange
       ! For example: 
       ! I0xpNe2 is the term of creation of Neu by charge exchange p-Ne2
       ! I0xpNe3 is the term of creation of Neu by charge exchange p-Ne3

       I0xp_I =RateN_I 
       I0px_I =RateN_I 

       I2xp_I  = Rate_I*(UStar_I/UStarM_I)*UThS_I(1)/No2Si_V(UnitU_)**2
       I2px_I  = Rate_I*(UStar_I/UStarM_I)*UThS_I/No2Si_V(UnitU_)**2

       ! units are fine: (Uth2/ustar)*termxp is unitless as it should be

       JxpUx_I  = Ux_I(1)*Rate_I 
       JxpUy_I  = Uy_I(1)*Rate_I   
       JxpUz_I  = Uz_I(1)*Rate_I 

       JpxUx_I  = Ux_I*Rate_I  
       JpxUy_I  = Uy_I*Rate_I   
       JpxUz_I  = Uz_I*Rate_I  

       !QmpxUx_I = JpxUx_I - JxpUx_I
       !QmpxUy_I = JpxUy_I - JxpUy_I
       !QmpxUz_I = JpxUz_I - JxpUz_I

       QmpxUx_I(1)=0.
       QmpxUy_I(1)=0.
       QmpxUz_I(1)=0.

       QmpxUx_I(Neu_) = (Ux_I(Neu_) - Ux_I(1))*Rate_I(Neu_)
       QmpxUx_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(1))*Rate_I(Ne2_)
       QmpxUx_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(1))*Rate_I(Ne3_)
       QmpxUx_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(1))*Rate_I(Ne4_)

       QmpxUy_I(Neu_) = (Uy_I(Neu_) - Uy_I(1))*Rate_I(Neu_)
       QmpxUy_I(Ne2_) = (Uy_I(Ne2_) - Uy_I(1))*Rate_I(Ne2_)
       QmpxUy_I(Ne3_) = (Uy_I(Ne3_) - Uy_I(1))*Rate_I(Ne3_)
       QmpxUy_I(Ne4_) = (Uy_I(Ne4_) - Uy_I(1))*Rate_I(Ne4_)

       QmpxUz_I(Neu_) = (Uz_I(Neu_) - Uz_I(1))*Rate_I(Neu_)
       QmpxUz_I(Ne2_) = (Uz_I(Ne2_) - Uz_I(1))*Rate_I(Ne2_)
       QmpxUz_I(Ne3_) = (Uz_I(Ne3_) - Uz_I(1))*Rate_I(Ne3_)
       QmpxUz_I(Ne4_) = (Uz_I(Ne4_) - Uz_I(1))*Rate_I(Ne4_)

       Kxp_I = 0.5*U2_I(1)*Rate_I  + I2xp_I   

       Kpx_I = 0.5*U2_I*Rate_I  + I2px_I   

       Qepx_I = Kpx_I - Kxp_I

       ! Calculate the source terms for this cell
       call calc_source_cell

    end do; end do; end do

  contains

    !==========================================================================
    subroutine calc_source_cell

      use ModPhysics,   ONLY: GammaMinus1_I

      ! Calculate source temrs for one cell. The pressures source is
      ! S(p) = (gamma-1)[S(e) - u.S(rhou) + 0.5 u**2 S(rho)]

      real:: Source_V(nVar + nFluid)
      integer:: iVar
      !---------------------------------------------------------------------

      Source_V = 0.0

      do iFluid = Neu_, Ne4_
         if(.not.UseSource_I(iFluid)) CYCLE
         call select_fluid
         if (iFluid == iFluidProduced_C(i,j,k)) then
            Source_V(iRho)    = sum(I0xp_I(Neu_:Ne4_))  - I0xp_I(iFluid)
            Source_V(iRhoUx)  = sum(JxpUx_I(Neu_:Ne4_)) - JpxUx_I(iFluid)
            Source_V(iRhoUy)  = sum(JxpUy_I(Neu_:Ne4_)) - JpxUy_I(iFluid)
            Source_V(iRhoUz)  = sum(JxpUz_I(Neu_:Ne4_)) - JpxUz_I(iFluid)
            Source_V(iEnergy) = sum(Kxp_I(Neu_:Ne4_))   - Kpx_I(iFluid)
         else
            Source_V(iRho)    = - I0px_I(iFluid)
            Source_V(iRhoUx)  = - JpxUx_I(iFluid)
            Source_V(iRhoUy)  = - JpxUy_I(iFluid)
            Source_V(iRhoUz)  = - JpxUz_I(iFluid)
            Source_V(iEnergy) = - Kpx_I(iFluid)
         end if
         Source_V(iP) = GammaMinus1_I(iFluid)* ( Source_V(iEnergy) &
              - Ux_I(iFluid)*Source_V(iRhoUx) &
              - Uy_I(iFluid)*Source_V(iRhoUy) &
              - Uz_I(iFluid)*Source_V(iRhoUz) &
              + 0.5*U2_I(iFluid)*Source_V(iRho) )
      end do

      if(UseSource_I(Ion_))then
         !!Source_V(RhoUx_) = sum( QmpxUx_I(Neu_:Ne4_) )
         !!Source_V(RhoUy_) = sum( QmpxUy_I(Neu_:Ne4_) )
         !!Source_V(RhoUz_) = sum( QmpxUz_I(Neu_:Ne4_) )
         Source_V(RhoUx_) = QmpxUx_I(Neu_) + QmpxUx_I(Ne2_) &
              + QmpxUx_I(Ne3_) + QmpxUx_I(Ne4_)
         Source_V(RhoUy_) = QmpxUy_I(Neu_) + QmpxUy_I(Ne2_) &
              + QmpxUy_I(Ne3_) + QmpxUy_I(Ne4_)
         Source_V(RhoUz_) = QmpxUz_I(Neu_) + QmpxUz_I(Ne2_) &
              + QmpxUz_I(Ne3_) + QmpxUz_I(Ne4_)

         Source_V(Energy_)= sum( Qepx_I(Neu_:Ne4_) )

         Source_V(p_) = GammaMinus1* ( Source_V(Energy_) &
              - Ux_I(Ion_)*Source_V(RhoUx_) &
              - Uy_I(Ion_)*Source_V(RhoUy_) &
              - Uz_I(Ion_)*Source_V(RhoUz_) ) 
      end if

      Source_VC(:,i,j,k) = Source_VC(:,i,j,k) + Source_V

      if(DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest)then

         Source_V = Source_VC(:,i,j,k)

         write(*,*) NameSub, ' iFluidProduced=', iFluidProduced_C(i,j,k)
         do iVar = 1, nVar + nFLuid
            write(*,*) ' Source(',NameVar_V(iVar),')=',Source_V(iVar)
         end do
         write(*,*) ' Temp_I    ', Temp_I
         write(*,*) ' Rate_I    ', Rate_I
         write(*,*) ' Sigma_I   ', Sigma_I
         write(*,*) ' UStar_I   ', UStar_I
         write(*,*) ' UStarM_I  ', UStarM_I
         write(*,*) ' Ux_I      ', Ux_I
         write(*,*) ' U2_I      ', U2_I
         write(*,*) ' UTh_I     ', sqrt(UThS_I)
         write(*,*) ' URelDim_I ', sqrt(URelSdim_I)
         write(*,*) ' uDim_I    ', sqrt(U2_I)*No2Io_V(UnitU_)
         write(*,*) ' I2xp_I    ', I2xp_I
         write(*,*) ' I2px_I    ', I2px_I
         write(*,*) ' JxpUx_I   ', JxpUx_I
         write(*,*) ' JpxUx_I   ', JpxUx_I
         write(*,*) ' QmpxUx_I  ', QmpxUx_I
         write(*,*) ' Kxp_I     ', Kxp_I
         write(*,*) ' Kpx_I     ', Kpx_I
         write(*,*) ' Qepx_I    ', Qepx_I 
      end if

    end subroutine calc_source_cell

  end subroutine user_calc_sources

  !============================================================================

  subroutine set_omega_parker_tilt

    ! Calculate angular velocity in normalized units
    ! Note: the rotation period is 25.38 days in ModConst.f90:
    ! OmegaSun = cTwoPi/(RotationPeriodSun*Si2No_V(UnitT_))

    OmegaSun   = cTwoPi/(26.0*24.0*3600.00*Si2No_V(UnitT_))
    ParkerTilt = OmegaSun*rBody/SWH_Ux

  end subroutine set_omega_parker_tilt

  !============================================================================

  subroutine select_region(iBlock)

    ! set the global variabls iFluidProduced_C
    ! to select which neutral fluid is produced in each cell of the block

    integer, intent(in):: iBlock

    integer :: i, j, k
    real    :: InvRho, U2, p, Mach2, TempDim, U2Dim, B2, MachAlfven2
    real    :: MachMagneto2

    character(len=*), parameter:: NameSub = 'select_region'
    !------------------------------------------------------------------------

    !This subroutine is not needed when not using the 4 neutral fluids                   
    if(.not.UseNeutralFluid) call CON_stop(NameSub//': no neutral fluids present')

    ! Produce fluid3 at the inner boundary

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       if(r_BLK(i,j,k,iBlock) < rPop3Limit) then
          iFluidProduced_C(i,j,k) = Ne3_
          CYCLE
       end if

       InvRho = 1.0/State_VGB(Rho_,i,j,k,iBlock)
       p      = State_VGB(p_,i,j,k,iBlock)
       U2     = sum(State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)**2)*InvRho**2
       U2Dim  = U2*No2Io_V(UnitU_)**2

       !merav modifications
       B2 = sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)        

       ! Square of Alfven Mach Number
       MachAlfven2 = U2/(B2*InvRho + 1.E-30)
       MachMagneto2 = U2/(1.E-10 + Gamma*p*InvRho + B2*InvRho)

       !end of merav modifications
       ! Square of Mach number
       Mach2      = U2/(Gamma*p*InvRho)

       ! Temperature in Kelvins
       TempDim = InvRho*p*No2Si_V(UnitTemperature_)

       ! Apply full source except near the boundaries between regions
!!!october11       if (MachPop4Limit**2 < Mach2 .and. uPop1LimitDim**2 > U2Dim) then
!!!       if (MachPop4Limit**2 < MachMagneto2 .and. uPop1LimitDim**2 > U2Dim) then
!!!       if (MachPop4Limit**2 < MachAlfven2 .and. uPop1LimitDim**2 > U2Dim .and. MachPop3Limit**2 < Mach2) then
!!!july12 use sonic Mach number - Berci
       if (MachPop4Limit**2 < Mach2 .and. uPop1LimitDim**2 > U2Dim) then  
          !Outside the bow shock
          iFluidProduced_C(i,j,k) = Ne4_
       elseif( TempPop1LimitDim > TempDim .and. uPop1LimitDim**2 > U2Dim)then
          !Outside the heliopause
          iFluidProduced_C(i,j,k) = Neu_
       elseif( MachPop2Limit**2 > Mach2 )then
          ! Heliosheath
          iFluidProduced_C(i,j,k) = Ne2_
       elseif( Mach2 > MachPop3Limit**2 )then
          ! Inside termination shock
          iFluidProduced_C(i,j,k) = Ne3_
       else
          ! No neutrals are produced in this region (but they are destroyed)
          iFluidProduced_C(i,j,k) = 0
       end if

    end do; end do; end do

  end subroutine select_region

  !============================================================================

  subroutine user_init_point_implicit

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, EpsPointImpl_V

    character(len=*), parameter:: NameSub = 'user_init_point_implicit'
    !------------------------------------------------------------------------

    !This subroutine is not needed when not using the 4 neutral fluids                    
    if(.not.UseNeutralFluid) call CON_stop(NameSub//': no neutral fluids present')

    ! Allocate and set iVarPointImpl_I
    ! In this example there are 3 implicit variables

    ! All the neutrals momenta and plasma are implicit 
    ! (3 neutral fluid and 1 ion)

    allocate(iVarPointImpl_I(20))

    iVarPointImpl_I = (/Rho_,RhoUx_, RhoUy_, RhoUz_, P_, &
         NeuRho_, NeuRhoUx_, NeuRhoUy_, NeuRhoUz_, NeuP_, &
         Ne2Rho_, Ne2RhoUx_, Ne2RhoUy_, Ne2RhoUz_, Ne2P_, &
         Ne3Rho_, Ne3RhoUx_, Ne3RhoUy_, Ne3RhoUz_, Ne3P_, &
         Ne4Rho_, Ne4RhoUx_, Ne4RhoUy_, Ne4RhoUz_, Ne4P_/)

    ! Because the second and third populations of neutrals have initially 
    ! very small values I'm
    ! setting EpsPointImpl_V to be small for these variables ??? !!!
    EpsPointImpl_V(Ne2Rho_:Ne4P_) = 1.e-11

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set.

    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

  !=======================================================================

  subroutine user_init_session

    use ModLookupTable,   ONLY: i_lookup_table

    ! Get index for the solar wind table holding time-dep. solar cycle values
    iTableSolarwind = i_lookup_table('solarwind2d')

  end subroutine user_init_session

  !=======================================================================
  subroutine user_set_boundary_cells(iBlock)
    use ModGeometry, ONLY: ExtraBc_, IsBoundaryCell_GI, Xyz_DGB, x2

    integer, intent(in):: iBlock
    integer:: i, j, k
    real:: x, y, z, x0, z0, Ratio, d, r
    !---------------------------------------------------------------------
    if(rCylinder > 0.0)then
       IsBoundaryCell_GI(:,:,:,ExtraBc_) = &
            Xyz_DGB(x_,:,:,:,iBlock)**2 + &
            Xyz_DGB(y_,:,:,:,iBlock)**2 < rCylinder**2 .and. &
            abs(Xyz_DGB(z_,:,:,:,iBlock)) < zCylinder 
       RETURN
    end if
    if(rCrescent < 0.0) RETURN

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI
       x = Xyz_DGB(x_,i,j,k,iBlock) - xCrescentCenter
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)

       ! point at the centerline of the crescent
       if(x < 0)then
          ! Half circle
          Ratio = xCrescentCenter/sqrt(x**2 + z**2)
          x0 = x*Ratio
          z0 = z*Ratio
       else
          ! horizontal lines
          x0 = x
          z0 = sign(xCrescentCenter, z)
       end if
       ! Distance from the centerline
       d = sqrt( (x-x0)**2 + y**2 + (z-z0)**2 )
       
       ! The radius of the crescent is a function of x. 
       if(x < 0.5*xCrescentCenter)then
          r = rCrescent
       else
          r = rCrescent*max(0.0, 1.5 - x/xCrescentCenter)
       end if
       IsBoundaryCell_GI(i,j,k,ExtraBc_) = d < r
    end do; end do; end do

  end subroutine user_set_boundary_cells

end module ModUser
