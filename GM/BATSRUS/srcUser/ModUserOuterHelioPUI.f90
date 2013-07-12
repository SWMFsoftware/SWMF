!#NOTPUBLIC  email:mopher@bu.edu  expires:12/31/2099

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
!February 2011 - PUI fluid
!April 2011 - Sources for PUI 
!July 2012 - Use sonic Mach number to define Region 1 - Region 4 boundary 
!==============================================================================
module ModUser

  use ModSize,     ONLY: nI, nJ, nK, nBLK
  use ModMain
  use ModPhysics
!  use ModSetOuterBC
  use ModAdvance,  ONLY : State_VGB
  use ModGeometry, ONLY : Xyz_DGB, r_BLK, true_cell
  use ModVarIndexes
  use ModProcMH
  use ModMultiFluid
  use ModUserEmpty,                                     &
       IMPLEMENTED1  => user_read_inputs,               &
       IMPLEMENTED2  => user_set_face_boundary,         &
       IMPLEMENTED3  => user_normalization,             &
       IMPLEMENTED4  => user_set_cell_boundary,         &
       IMPLEMENTED5  => user_set_ics,                   &
       IMPLEMENTED6  => user_initial_perturbation,      &
       IMPLEMENTED8  => user_action,                    & 
       IMPLEMENTED9  => user_io_units,                  &
       IMPLEMENTED10 => user_set_plot_var,              &
       IMPLEMENTED11 => user_calc_sources,              &
       IMPLEMENTED12 => user_init_point_implicit


  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 1.0
  character (len=*), parameter :: &
       NameUserModule = 'Outer Heliosphere with 4 neutrals and 2 ion fluids, Opher & Toth'

  ! Named indexes for fluids
   integer, parameter :: ALL_ =1 , SWH_ = 2, Pu3_ = 3, Neu_ = 4, Ne2_ = 5, Ne3_ = 6, Ne4_= 7 ! defined in ModEquation

  logical :: UseSource_I(SWH_:Ne4_) = .true.

  real :: OmegaSun   = 0.0  ! normalized angular speed of Sun
  real :: ParkerTilt = 0.0  ! Bphi/Br at the equator at r=rBody

  ! SWH variables.
  !/
  real :: SWH_a_dim=0.0  , &
       SWH_rho=0.0, SWH_rho_dim=0.0, &
       SWH_p=0.0  , SWH_T_dim  =0.0, &
       SWH_Ux=0.0 , SWH_Ux_dim=0.0 , &
       SWH_Uy=0.0 , SWH_Uy_dim=0.0 , &
       SWH_Uz=0.0 , SWH_Uz_dim=0.0 , &
       SWH_Bx=0.0 , SWH_Bx_dim=0.0 , &
       SWH_By=0.0 , SWH_By_dim=0.0 , &
       SWH_Bz=0.0 , SWH_Bz_dim=0.0 , &
       SWH_B_factor=0.0

! Is this used?
  real, dimension(0:1) :: &
       SWH_rho_t,  &
       SWH_p_t  ,  &
       SWH_Ux_t ,  &
       SWH_Uy_t ,  &
       SWH_Uz_t ,  &
       SWH_Bx_t ,  &
       SWH_By_t ,  &
       SWH_Bz_t ,  &
       SWH_time_t
  ! 
  ! VLISM variables.
  !/
  real ::      SW_B_factor=0.0  
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
  
  real :: SWH_p1=0.0, PNeutralsISW1=0.0, PU3_p1=0.0
 
  real :: SWfast_T_dim=0.0, &
       SWfast_a_dim=0.0, &
       SWfast_rho=0.0, SWfast_rho_dim=0.0, &
       SWfast_p=0.0  , SWfast_p_dim=0.0  , &
       SWfast_Ux=0.0 , SWfast_Ux_dim=0.0 , &
       SWfast_Uy=0.0 , SWfast_Uy_dim=0.0 , &
       SWfast_Uz=0.0 , SWfast_Uz_dim=0.0,  &
       SWfast_Bx=0.0 , SWfast_Bx_dim=0.0 , &
       SWfast_By=0.0 , SWfast_By_dim=0.0 , &
       SWfast_Bz=0.0 , SWfast_Bz_dim=0.0 , &
       SWfast_B_factor=0.0

  ! B. Zieger commented out 05/17/2012
  ! specified by use statement in MH_set_parameters

  !real :: TNeutralsISW_dim=0.0, &
  !     RhoNeutralsISW=0.0, RhoNeutralsISW_dim=0.0 , &
  !     PNeutralsISW=0.0  , PNeutralsISW_dim=0.0  , &
  !     UxNeutralsISW=0.0 , UxNeutralsISW_dim=0.0 , &
  !     UyNeutralsISW=0.0 , UyNeutralsISW_dim=0.0 , &
  !     UzNeutralsISW=0.0 , UzNeutralsISW_dim=0.0 ,  &
  !     mNeutralsmp, mNeutrals

  real :: mNeutralsmp, mNeutrals

  real, dimension(0:1) :: &
       RhoNeutralsISW_t,  &
       PneutralsISW_t  ,  &
       UxNeutralsISW_t ,  &
       UyNeutralsISW_t ,  &
       UzNeutralsISW_t

  ! Velocity, temperature, Mach number and radius limits for the populations
  real :: TempPop1LimitDim = 1e5    ! [K]
  real :: uPop1LimitDim    = 100.0  ! [km/s]
  real :: MachPop2Limit    = 0.9
  real :: MachPop3Limit    = 1.5
  real :: MachPop4Limit    = 2.0
  real :: rPop3Limit       = 50.0   ! [AU] it is all Pop3 out to rPop3Limit

  ! Various factors in initial and boundary conditions
  real :: RhoNeuFactor = 1.0,   uNeuFactor = 1.0
  real :: RhoNe2Factor = 1.e-3, uNe2Factor = 1.0
  real :: RhoNe3Factor = 0.01,  uNe3Factor = 1.0
  real :: RhoNe4Factor = 1.0,  uNe4Factor = 1.0

  integer :: iFluidProduced_C(nI, nJ, nK)

  real :: Pu3_a_dim=0.0  , &
       Pu3_rho=0.0, Pu3_rho_dim=0.0, &
       Pu3_p=0.0  , Pu3_T_dim  =0.0, &
       Pu3_Ux=0.0 , Pu3_Ux_dim=0.0 , &
       Pu3_Uy=0.0 , Pu3_Uy_dim=0.0 , &
       Pu3_Uz=0.0 , Pu3_Uz_dim=0.0 
  ! Is this used?
  real, dimension(0:1) :: &
       Pu3_rho_t,  &
       Pu3_p_t  ,  &
       Pu3_Ux_t ,  &
       Pu3_Uy_t ,  &
       Pu3_Uz_t ,  &
       Pu3_time_t
 !some extra variables for constant pressure at inner boundary
  real :: pPUI_30 = 0.0 , &
          pSolarWind_30 = 0.0
contains

  !=========================================================================

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam

    character (len=100) :: NameCommand
    character (len=*), parameter :: Name='user_read_inputs'
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
       case("#NEUTRALS")
          call read_var('RhoNeutralsISW_dim' ,RhoNeutralsISW_dim)
          call read_var('TNeutralsISW_dim' ,TNeutralsISW_dim)
          call read_var('UxNeutralsISW_dim' ,UxNeutralsISW_dim)
          call read_var('UyNeutralsISW_dim' ,UyNeutralsISW_dim)
          call read_var('UzNeutralsISW_dim' ,UzNeutralsISW_dim)
          call read_var('mNeutralsmp',mNeutralsmp)
	case("#PICKUPION3")
          call read_var('Pu3_rho_dim',Pu3_rho_dim)
          call read_var('Pu3_T_dim'  ,Pu3_T_dim)
          call read_var('Pu3_Ux_dim' ,Pu3_Ux_dim)
          call read_var('Pu3_Uy_dim' ,Pu3_Uy_dim)
          call read_var('Pu3_Uz_dim' ,Pu3_Uz_dim)
 
          ! This is a flag to define how many populations of Neutrals to run
          case("#SOURCES")
          call read_var('UseSWHSource', UseSource_I(SWH_))
          call read_var('UsePu3Source', UseSource_I(Pu3_))
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
       case("#FACTORS")
          call read_var('RhoNeuFactor', RhoNeuFactor)
          call read_var('uNeuFactor'  , uNeuFactor)
          call read_var('RhoNe2Factor', RhoNe2Factor)
          call read_var('uNe2Factor'  , uNe2Factor)
          call read_var('RhoNe3Factor', RhoNe3Factor)
          call read_var('uNe3Factor'  , uNe3Factor)
          call read_var('RhoNe4Factor', RhoNe4Factor)
          call read_var('uNe4Factor'  , uNe4Factor)
       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
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
    !C.P. edited
    real:: Bsph_D(3), Vsph_D(3), VPUIsph_D(3)

    real :: pSolarWind,pPUI, p_frac, Ptot, Pmag, PmagEquator

    real :: XyzSph_DD(3,3) 

    logical :: DoTest, DoTestMe
    character(len=*), parameter:: NameSub='user_set_face_boundary'
    !-------------------------------------------------------------------

    if(iBoundary /= body1_) &
         call stop_mpi(NameSub//' only inner BC is implemented!')
    

    if(iProc == ProcTest .and. iBlockBc == BlkTest)then

       call set_oktest(NameSub, DoTest, DoTestMe)

    else
       DoTest = .false.; DoTestMe = .false.
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
    Bsph_D(1) =  sign(SWH_Bx, zFace)             ! Br  !good for 2005 field
    Bsph_D(2) =  0.0                             ! Btheta
    Bsph_D(3) = -Bsph_D(1)*SinTheta*ParkerTilt   ! Bphi

    Vsph_D    = (/ SWH_Ux, 0.0, 0.0 /)           ! Vr, Vtheta, Vphi
    VPUIsph_D = (/ Pu3_Ux, 0.0, 0.0 /)    
   
    ! Calculate pressure (equilibrium at a given inner boundary)
    Pmag = sum(Bsph_D**2) / 2.0

    ! magnetic pressure at the equator (actually wrong, neglects Br=SWH_Bx)
    PmagEquator = (SWH_Bx*ParkerTilt)**2/2
    
    !Merav says numerical reasons only and FIX FOR PUI

    !this method garruntees equal pressure over inner boundary.  Simplier splitting method likely pl
   p_frac = SWH_p/(SWH_p+PU3_p) 
   pSolarWind = SWH_p + p_frac*(PmagEquator - Pmag)
   pPUI = PU3_p +(1-p_frac)* (PmagEquator - Pmag)
  
  !new inner boundary pressures  
  pPUI_30 = pPUI
  PSolarWind_30 = pSolarWind

    ! Apply boundary conditions for ions
    VarsGhostFace_V(SWHRho_)    = SWH_rho
    VarsGhostFace_V(SWHP_)      = pSolarWind ! SWH_p 
    VarsGhostFace_V(SWHUx_:SWHUz_) = matmul(XyzSph_DD, Vsph_D)
    VarsGhostFace_V(Bx_:Bz_) = matmul(XyzSph_DD, Bsph_D)
 
    VarsGhostFace_V(Pu3Rho_)    = PU3_rho
    VarsGhostFace_V(Pu3p_)      = pPUI  !PU3_p
    VarsGhostFace_V(Pu3Ux_:Pu3Uz_) = matmul(XyzSph_DD, VPUIsph_D)
    !
    
    !Total ion fluid, should be less hard coded, use ion notation later 
    VarsGhostFace_V(Rho_)    = VarsGhostFace_V(SWHRho_) + VarsGhostFace_V(Pu3Rho_)
    VarsGhostFace_V(P_)      = VarsGhostFace_V(SWHp_) + VarsGhostFace_V(Pu3P_)
    VarsGhostFace_V(Ux_:Uz_) =  (VarsGhostFace_V(SWHRho_)*VarsGhostFace_V(SWHUx_:SWHUz_)&
         + VarsGhostFace_V(Pu3Rho_)*VarsGhostFace_V(Pu3Ux_:Pu3Uz_))/VarsGhostFace_V(Rho_)

    ! NeuRho is PopI
    ! Pop I is going through the inner BCs    

   ! soft boundary for Pop I-IV

       VarsGhostFace_V(NeuRho_:NeuP_) = VarsTrueFace_V(NeuRho_:NeuP_)

    ! PopII leaves the domain at a supersonic velocity 
    ! (50km/s while for their temperature 1.E5K their C_s=30km/s)
    ! For the transient case when it flows inward, we use a fraction of ion parameters

    if( sum(VarsTrueFace_V(Ne2Ux_:Ne2Uz_)*FaceCoords_D) > 0.0)then
       VarsGhostFace_V(Ne2Rho_)       = VarsGhostFace_V(Rho_)    *RhoNe2Factor
       VarsGhostFace_V(Ne2P_)         = VarsGhostFace_V(p_)      *RhoNe2Factor
       VarsGhostFace_V(Ne2Ux_:Ne2Uz_) = VarsGhostFace_V(Ux_:Uz_) *uNe2Factor
    else
       VarsGhostFace_V(Ne2Rho_:Ne2P_) = VarsTrueFace_V(Ne2Rho_:Ne2P_)
    end if

    ! Pop III has the velocity and temperature of the ions at inner boundary
    ! the density is taken to be a fraction of the ions

   if( sum(VarsTrueFace_V(Ne3Ux_:Ne3Uz_)*FaceCoords_D) > 0.0)then
       VarsGhostFace_V(Ne3Rho_)       = VarsGhostFace_V(Rho_)    *RhoNe3Factor
       VarsGhostFace_V(Ne3P_)         = VarsGhostFace_V(p_)      *RhoNe3Factor
       VarsGhostFace_V(Ne3Ux_:Ne3Uz_) = VarsGhostFace_V(Ux_:Uz_) *uNe3Factor
   else
       VarsGhostFace_V(Ne3Rho_:Ne3P_) = VarsTrueFace_V(Ne3Rho_:Ne3P_)
   end if

    ! Pop IV 
   if( sum(VarsTrueFace_V(Ne4Ux_:Ne4Uz_)*FaceCoords_D) > 0.0)then
       VarsGhostFace_V(Ne4Rho_)       = VarsGhostFace_V(Rho_)    *RhoNe4Factor
       VarsGhostFace_V(Ne4P_)         = VarsGhostFace_V(p_)      *RhoNe4Factor
       VarsGhostFace_V(Ne4Ux_:Ne4Uz_) = VarsGhostFace_V(Ux_:Uz_) *uNe4Factor
   else
       VarsGhostFace_V(Ne4Rho_:Ne4P_) = VarsTrueFace_V(Ne4Rho_:Ne4P_)
   end if

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
       write(*,*) NameSub,' Pop1     =', VarsGhostFace_V(NeuRho_:NeuP_)
       write(*,*) NameSub,' Pop2     =', VarsGhostFace_V(Ne2Rho_:Ne2P_)
       write(*,*) NameSub,' Pop3     =', VarsGhostFace_V(Ne3Rho_:Ne3P_)
      !added by C.P.
       write(*,*) NameSub,' Pop4     =', VarsGhostFace_V(Ne4Rho_:Ne4P_)
       write(*,*) NameSub,' VPUIsph_D=', VPUIsph_D
       write(*,*) NameSub,' Pui3     =', VarsGhostFace_V(Pu3Rho_:Pu3p_) 
       write(*,*) NameSub,' pPUI, PU3_p, pSolarwind, SWH_p =', pPUI, PU3_p, pSolarwind, SWH_p 
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


    No2Si_V(UnitX_)  = cAU                                      ! m
    No2Si_V(UnitU_)  = sqrt(g*cBoltzmann*SWH_T_dim/cProtonMass) ! m/s
    No2Si_V(UnitRho_)= cProtonMass*SWH_rho_dim*1.0E+6           ! kg/m^3

    if(DoTestMe)then
       write(*,*)NameSub,' No2Si_V(UnitX_)  =',No2Si_V(UnitX_)
       write(*,*)NameSub,' No2Si_V(UnitU_)  =',No2Si_V(UnitU_)
       write(*,*)NameSub,' No2Si_V(UnitRho_)=',No2Si_V(UnitRho_)
    end if

  end subroutine user_normalization
  !-------------------------------------------------------------------

  subroutine user_set_cell_boundary(iBlock, iSide, TypeBc, IsFound)

    ! The ISM enters at the east boundary (negative x)

    use ModMain
    use ModVarIndexes
    use ModProcMH
    use ModAdvance, ONLY : State_VGB
    use ModMultiFluid

    integer,intent(in)::iBlock, iSide
    character (len=*),intent(in) :: TypeBc
    logical,intent(out) :: IsFound

    integer :: i, iVar
    !----------------------------------------------------------------------
    IsFound = .true.
    State_VGB(SWHrho_,-1:2,:,:,iBlock)=VLISW_rho 
    State_VGB(SWHrhoUx_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Ux
    State_VGB(SWHrhoUy_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Uy
    State_VGB(SWHrhoUz_,-1:2,:,:,iBlock)=VLISW_rho*VLISW_Uz
    State_VGB(Bx_,-1:2,:,:,iBlock)=VLISW_Bx
    State_VGB(By_,-1:2,:,:,iBlock)=VLISW_By
    State_VGB(Bz_,-1:2,:,:,iBlock)=VLISW_Bz
    State_VGB(SWHp_,-1:2,:,:,iBlock)=VLISW_p
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

    !
    !\

    ! PopII and III supersonic outflow
    !/

    do iVar = NeuRho_, Ne3P_; do i = -1, 2
       State_VGB(iVar,i,:,:,iBlock) = State_VGB(iVar,3,:,:,iBlock)
    end do; end do

    !do iVar = Pu3Rho_, Pu3P_; do i = -1, 2
    !   State_VGB(iVar,i,:,:,iBlock) = State_VGB(iVar,3,:,:,iBlock)
    State_VGB(Pu3Rho_,-1:2,:,:,iBlock) = 10E-6*VLISW_rho
    State_VGB(Pu3RhoUx_,-1:2,:,:,iBlock) = 10E-6*VLISW_rho*VLISW_Ux
    State_VGB(Pu3RhoUy_,-1:2,:,:,iBlock) =10E-6*VLISW_rho*VLISW_Uy
    State_VGB(Pu3RhoUz_,-1:2,:,:,iBlock) =10E-6*VLISW_rho*VLISW_Uz
    State_VGB(Pu3P_,-1:2,:,:,iBlock) = 10E-6*VLISW_p
    ! end do; end do

    State_VGB(Rho_,-1:2,:,:,iBlock)  = State_VGB(SWHRho_,-1:2,:,:,iBlock) +  State_VGB(Pu3Rho_,-1:2,:,:,iBlock)
    State_VGB(P_,-1:2,:,:,iBlock)  = State_VGB(SWHP_,-1:2,:,:,iBlock) + State_VGB(Pu3P_,-1:2,:,:,iBlock)
    State_VGB(RhoUx_,-1:2,:,:,iBlock)  = (State_VGB(SWHRhoUx_,-1:2,:,:,iBlock) +  State_VGB(Pu3RhoUx_,-1:2,:,:,iBlock))
    State_VGB(RhoUy_,-1:2,:,:,iBlock)  = (State_VGB(SWHRhoUy_,-1:2,:,:,iBlock) +  State_VGB(Pu3RhoUy_,-1:2,:,:,iBlock))
    State_VGB(RhoUz_,-1:2,:,:,iBlock)  = (State_VGB(SWHRhoUz_,-1:2,:,:,iBlock) +  State_VGB(Pu3RhoUz_,-1:2,:,:,iBlock)) 

  end subroutine user_set_cell_boundary

  !=====================================================================

  subroutine user_set_ics(iBlock)

    use ModVarIndexes    
    use ModAdvance,  ONLY: State_VGB    
    !C.P. edited
    use ModPhysics,  ONLY: rBody, No2Si_V, UnitX_
    use ModCoordTransform, ONLY: rot_xyz_sph
    !C.P. added
    use ModConst, ONLY: cAU

    integer, intent(in) :: iBlock

    integer :: i,j,k
    real :: x, y, z, r, rho0
    real :: b_D(3), v_D(3), bSph_D(3), vSph_D(3), vPUI_D(3), vPUISph_D(3)
    real :: SinTheta, SignZ

    ! WARNING: ASSUMES 30 AU INNER BOUNDARY
    ! may wish to use more accuate 1 AU and sigma
    real :: r30  !cm
    real :: sigma !cm^2

    real :: XyzSph_DD(3,3) ! rotation matrix Xyz_D = matmul(XyzSph_DD, Sph_D)

    character(len=*), parameter:: NameSub = 'user_set_ics'
    logical :: DoTest, DoTestMe, DoTestCell
    !No longer needed? 
    !--------------------------------------------------------------------------
    r30 = 30.0 * cAU *100.0  ! in cm, cAU in m
    sigma = 2E-15 !cm^2
    !--------------------------------------------------------------------------
    if(iProc == ProcTest .and. iBlock == BlkTest)then

       call set_oktest(NameSub, DoTest, DoTestMe)

    else
       DoTest = .false.; DoTestMe = .false.
    end if

    !Make sure that OmegaSun and ParkerTilt are set
    if(OmegaSun == 0.0)call set_omega_parker_tilt

    do k = MinK, MaxK; do j = MinJ, MaxJ; do i = MinI, MaxI


       if(.not. true_cell(i,j,k,iBlock)) CYCLE

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       x = Xyz_DGB(x_,i,j,k,iBlock)
       y = Xyz_DGB(y_,i,j,k,iBlock)
       z = Xyz_DGB(z_,i,j,k,iBlock)
       r = r_BLK(i,j,k,iBlock)

       XyzSph_DD = rot_xyz_sph(x,y,z)

       ! theta is the angle measured from the pole
       SinTheta = sqrt(x**2+y**2)/r

       ! for the neutrals thetaN angle is relative to the X axis
       ! thetaN = atan(sqrt(y**2+z**2)/(x+cTiny))
       ! sinThetaN=sqrt(y**2+z**2)/r
       ! lambda = 4.0
       !so pra iniciar eu fixei lambda como 4.0 
       !(olha a expressao 3.53 da tese do Linde)

       !calculating the Parker B field spherical components Bsph_D

       ! good for polarity of 1997       SignZ = sign(1.0, z)
       SignZ = sign(1.0,z)  ! good for 2005

       Bsph_D(1) = SignZ*SWH_Bx*(rBody/r)**2  ! Br
       Bsph_D(2) = 0.0                        ! Btheta
       Bsph_D(3) = -SignZ*SWH_Bx*SinTheta*ParkerTilt*(rBody/r) !Bphi


       Vsph_D = (/ SWH_Ux, 0., 0. /)
       vPUISph_D = (/ Pu3_Ux, 0., 0. /)

       ! magnetic field components in cartesian coordinates
       b_D = matmul(XyzSph_DD, Bsph_D)

       State_VGB(Bx_:Bz_,i,j,k,iBlock) = b_D


       !velocity components in cartesian coordinates
       v_D = matmul(XyzSph_DD, Vsph_D)
       vPUI_D = matmul(XyzSph_DD, VPUIsph_D)
       ! density and pressure
       State_VGB(SWHRho_,i,j,k,iBlock) = SWH_rho * (rBody/r)**2
       State_VGB(SWHP_,i,j,k,iBlock)   = SWH_p   * (rBody/r)**(2*g) 

       ! momentum
       State_VGB(SWHRhoUx_:SWHRhoUz_,i,j,k,iBlock) = State_VGB(SWHrho_,i,j,k,iBlock)*v_D


       !\
       ! PopI
       !/
       State_VGB(NeuRho_,i,j,k,iBlock)  =  RhoNeutralsISW
       State_VGB(NeuP_,i,j,k,iBlock) =   PNeutralsISW
       State_VGB(NeuRhoUx_,i,j,k,iBlock)=RhoNeutralsISW*UxNeutralsISW
       State_VGB(NeuRhoUy_,i,j,k,iBlock)=RhoNeutralsISW*UyNeutralsISW
       State_VGB(NeuRhoUz_,i,j,k,iBlock)=RhoNeutralsISW*UzNeutralsISW

       !! This profile takes into account loss due to PUI
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
       State_VGB(Ne2Rho_,i,j,k,iBlock) = 1.e-3 * RhoNeutralsISW
       State_VGB(Ne2P_,i,j,k,iBlock)   = 1.e-3 * PNeutralsISW 
       State_VGB(Ne2RhoUx_:Ne2RhoUz_,i,j,k,iBlock) = &
            0.25*State_VGB(Ne2Rho_,i,j,k,iBlock)*v_D
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


       ! If filling in with realistic PUI profile must do after neutrals populated
       ! neglect photoionization at >30 AU, ~10% effect
       ! must be a-dimensional
       !start with PUI's only to 100 AU
       !if (r <= rPop3Limit) then

       ! No production yet
       State_VGB(Pu3Rho_,i,j,k,iBlock) = Pu3_rho * (rBody/r)**2
       State_VGB(Pu3P_,i,j,k,iBlock)   = Pu3_p   * (rBody/r)**(2*g)

       ! as if there was production
       !  State_VGB(Pu3Rho_,i,j,k,iBlock) = SWH_rho_dim * sigma * r30**2 / (r*Io2Si_V(UnitX_)*100.0)* (State_VGB(NeuRho_,i,j,k,iBlock) + &
       !      State_VGB(Ne2Rho_,i,j,k,iBlock) + State_VGB(Ne3Rho_,i,j,k,iBlock) + State_VGB(Ne4Rho_,i,j,k,iBlock))   ! placeholder PUI density equation
       !    State_VGB(Pu3P_,i,j,k,iBlock)   = State_VGB(Pu3Rho_,i,j,k,iBlock)*Pu3_T_dim *Io2No_V(UnitTemperature_) 
       !else       
       !set to small value
       !   State_VGB(Pu3Rho_,i,j,k,iBlock) = 10E-6*State_VGB(Rho_,i,j,k,iBlock) 
       !   State_VGB(Pu3P_,i,j,k,iBlock) = 10E-6*State_VGB(P_,i,j,k,iBlock)
       !end if

       !momentum
       State_VGB(Pu3RhoUx_:Pu3RhoUz_,i,j,k,iBlock) = &
            State_VGB(Pu3Rho_,i,j,k,iBlock)*vPUI_D
       !


       State_VGB(Rho_,i,j,k,iBlock)  = &
            State_VGB(SWHRho_,i,j,k,iBlock) +  State_VGB(Pu3Rho_,i,j,k,iBlock)
       State_VGB(P_,i,j,k,iBlock) = &
            State_VGB(SWHP_,i,j,k,iBlock) + State_VGB(Pu3P_,i,j,k,iBlock)
       !momentum
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock)  = &
            State_VGB(SWHRhoUx_:SWHRhoUz_,i,j,k,iBlock) &
            + State_VGB(Pu3RhoUx_:Pu3RhoUz_,i,j,k,iBlock)

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
	  write(*,*)NameSub,' VPUIsph_D                 =',VPUIsph_D
          write(*,*)NameSub,' vPUI_D                    =',vPUI_D
          write(*,*)NameSub,' Pu3Rho,   =',State_VGB(Pu3Rho_,i,j,k,iBlock)
          write(*,*)NameSub,' Pu3p      =',State_VGB(Pu3P_,i,j,k,iBlock)

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
    do iBlock = 1, nBlock

       if( Unused_B(iBlock) ) CYCLE

       State_VGB(Ne2Rho_,:,:,:,iBlock) = RhoNe2Factor * &
            State_VGB(Rho_,:,:,:,iBlock)

       State_VGB(Ne2RhoUx_:Ne2RhoUz_,:,:,:,iBlock) = RhoNe2Factor*uNe2Factor* &
            State_VGB(RhoUx_:RhoUz_,:,:,:,iBlock)

       State_VGB(Ne2P_,:,:,:,iBlock) = RhoNe2Factor * &
            State_VGB(p_,:,:,:,iBlock)

       call calc_energy(-1, nI+2, -1, nJ+2, -1, nK+2, iBlock, Ne2_, Ne2_)
    end do

  end subroutine user_initial_perturbation

  !=====================================================================

  subroutine user_action(NameAction)
    use ModMain
    use ModPhysics
    use ModMultiFluid

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
    !neutrals
    write(*,*)     
    write(*,StringFormat) 'RhoNeutralsISW_dim:',RhoNeutralsISW_dim ,'RhoNeutralsISW:',RhoNeutralsISW 
    write(*,StringFormat) 'UxNeutralsISW_dim:',UxNeutralsISW_dim,'UxNeutralsISW:',UxNeutralsISW 
    write(*,StringFormat) 'UyNeutralsISW_dim:',UyNeutralsISW_dim,'UyNeutralsISW:',UyNeutralsISW
    write(*,StringFormat) 'UzNeutralsISW_dim:',UzNeutralsISW_dim,'UzNeutralsISW:',UzNeutralsISW 
    write(*,StringFormat) 'PNeutralsISW_dim:',PNeutralsISW_dim,'PNeutralsISW:',PNeutralsISW
    write(*,'(10X,A19,F15.6)') 'TNeutralsISW_dim:',TNeutralsISW_dim     
    write(*,*)
    !PUIs
    !C.P. added
    write(*,*)     
    write(*,StringFormat) 'PU3_rho_dim [n/cc]:',PU3_rho_dim,'SWH_rho:',PU3_rho
    write(*,StringFormat) 'PU3_Ux_dim  [km/s]:',PU3_Ux_dim,'PU3_Ux:',PU3_Ux
    write(*,StringFormat) 'PU3_Uy_dim  [km/s]:',PU3_Uy_dim,'PU3_Uy:',PU3_Uy 
    write(*,StringFormat) 'PU3_Uz_dim  [km/s]:',PU3_Uz_dim,'PU3_Uz:',PU3_Uz
    write(*,StringFormat) 'PU3_T_dim   [   K]:',PU3_T_dim,'PU3_p:',PU3_p
    write(*,'(10X,A19,F15.6)')           'PU3_T_dim   [   K]:',PU3_T_dim
    write(*,*)  
    !------------------------------------------------------------------

  end subroutine user_action

  !=====================================================================
  subroutine user_io_units

    use ModPhysics
    use ModConst, ONLY: cAU, cProtonMass
    !
    character (len=*), parameter :: Name='user_io_units'
    !-------------------------------------------------------------------------
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
    VLISW_p_dim1    = No2Io_V(UnitP_)*inv_g &
         *(VLISW_rho_dim/SWH_rho_dim)*(VLISW_T_dim/SWH_T_dim)

    ! Pressure of plasma = 2*T_ion*rho_ion

    VLISW_B_factor = No2Io_V(UnitB_)*sqrt((VLISW_T_dim/SWH_T_dim) &
         *(VLISW_rho_dim/SWH_rho_dim))

    VLISW_rho = VLISW_rho_dim*Io2No_V(UnitRho_)
    VLISW_p1   = VLISW_p_dim1*Io2No_V(UnitP_)
    VLISW_p    = 2.*VLISW_T_dim*Io2No_V(UnitTemperature_)*VLISW_rho

    VLISW_Ux  = VLISW_Ux_dim*Io2No_V(UnitU_)
    VLISW_Uy  = VLISW_Uy_dim*Io2No_V(UnitU_)
    VLISW_Uz  = VLISW_Uz_dim*Io2No_V(UnitU_)
    VLISW_Bx  = VLISW_Bx_dim*Io2No_V(UnitB_)
    VLISW_By  = VLISW_By_dim*Io2No_V(UnitB_)
    VLISW_Bz  = VLISW_Bz_dim*Io2No_V(UnitB_)

    SWH_rho = SWH_rho_dim*Io2No_V(UnitRho_)
    ! Pressure of plasma = 2*T_ion*rho_ion

    SWH_p1   = SWH_T_dim*Io2No_V(UnitTemperature_)*SWH_rho
    SWH_p   = 2.*SWH_T_dim*Io2No_V(UnitTemperature_)*SWH_rho

    SWH_Ux  = SWH_Ux_dim*Io2No_V(UnitU_)
    SWH_Uy  = SWH_Uy_dim*Io2No_V(UnitU_)
    SWH_Uz  = SWH_Uz_dim*Io2No_V(UnitU_)
    SWH_Bx  = SWH_Bx_dim*Io2No_V(UnitB_)
    SWH_By  = SWH_By_dim*Io2No_V(UnitB_)
    SWH_Bz  = SWH_Bz_dim*Io2No_V(UnitB_)

    PU3_Ux  = PU3_Ux_dim*Io2No_V(UnitU_)
    PU3_Uy  = PU3_Uy_dim*Io2No_V(UnitU_)
    PU3_Uz  = PU3_Uz_dim*Io2No_V(UnitU_)
    PU3_rho = PU3_rho_dim*Io2No_V(UnitRho_)
    ! Pressure of sw plasma = 2*T_ion*rho_ion, assume electron pressure  = solar wind press (ok to same order), so no factor of 2 for PUI
    PU3_p1   = PU3_T_dim*Io2No_V(UnitTemperature_)*PU3_rho
    PU3_p   = PU3_T_dim*Io2No_V(UnitTemperature_)*PU3_rho

    ! The units of rho_dim are n/cc and unitUSER_rho g/cc
    !/

    RhoNeutralsISW = RhoNeutralsISW_dim*Io2No_V(UnitRho_)
    PNeutralsISW_dim = No2Io_V(UnitP_)*inv_g*(RhoNeutralsISW_dim/SWH_rho_dim)*(TNeutralsISW_dim /SWH_T_dim)

    !PNeutralsISW1 = PNeutralsISW_dim*Io2No_V(UnitP_)
    PNeutralsISW = TNeutralsISW_dim*Io2No_V(UnitTemperature_)*RhoNeutralsISW 
    UxNeutralsISW  = UxNeutralsISW_dim*Io2No_V(UnitU_)
    UyNeutralsISW  = UyNeutralsISW_dim*Io2No_V(UnitU_)
    UzNeutralsISW  = UzNeutralsISW_dim*Io2No_V(UnitU_)
    mNeutrals    = mNeutralsmp*cProtonMass

    ! set strings for writing Tecplot output
    !/
    NameTecUnit_V(UnitX_)            = 'AU'
    NameTecUnit_V(UnitRho_)          = '#/cm3'
    NameTecUnit_V(UnitU_)            = 'km/s'
    NameTecUnit_V(UnitP_)            = 'dyne/cm^2'
    NameTecUnit_V(UnitB_)            = 'nT'

  end subroutine user_io_units
  !==============================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModSize, ONLY: nI, nJ, nK
    use ModAdvance,  ONLY: Source_VC

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

    character (len=*), parameter :: Name='user_set_plot_var'

    !-------------------------------------------------------------------

    UsePlotVarBody = .true.
    PlotVarBody    = 0.0

    IsFound = .true.
    select case(NameVar)
    case('srho')
       NameTecVar = 'Srho'
       PlotVar_G(1:nI,1:nJ,1:nK) = Source_VC(NeuRho_,:,:,:)
    case('fluid')
       call select_region(iBlock)
       PlotVar_G(1:nI,1:nJ,1:nK) = iFluidProduced_C
    case('mach')
       !C.P. edited to add PUI
       PlotVar_G = &
            sqrt( sum((State_VGB(SWHRhoUx_:SWHRhoUz_,:,:,:,iBlock) &
            + State_VGB(Pu3RhoUx_:Pu3RhoUz_,:,:,:,iBlock))**2, DIM=1) &
            / (g *(State_VGB(SWHp_,:,:,:,iBlock) + State_VGB(Pu3p_,:,:,:,iBlock)) &
            * (State_VGB(SWHRho_,:,:,:,iBlock)+State_VGB(Pu3Rho_,:,:,:,iBlock))) )
       !merav addition
       !edited to add PUI, could make more generic from first ion to last ion
    case('machalfven')
       PlotVar_G = & 
            sqrt( sum((State_VGB(SWHRhoUx_:SWHRhoUz_,:,:,:,iBlock) &
            + State_VGB(Pu3RhoUx_:Pu3RhoUz_,:,:,:,iBlock))**2, DIM=1)      &
            /   (( sum(State_VGB(Bx_:Bz_,:,:,:,iBlock)**2) &
            * (State_VGB(SWHRho_,:,:,:,iBlock) + State_VGB(Pu3Rho_,:,:,:,iBlock)) ) ) ) 

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
    ! Ne4 are the neutrals created between bow shock and heliopause (PopIV)

    ! As an example for Neu the source terms inside the Termination Shock will 
    ! take into account their creations and destruction; 
    ! outside the Termination Shock they will be just destroyed 
    ! (see more details below). 
    !
    ! _I(1) is the SW ionized fluid, _I(2) is the Pick-up ion fluid and _I(3)-_I(6) are the neutral fluids
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
    use ModProcMH
    use ModPointImplicit, ONLY:  UsePointImplicit, UsePointImplicit_B, &
         IsPointImplSource, IsPointImplPerturbed
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModPhysics
    use ModNumConst

    integer, intent(in) :: iBlock

    character (len=*), parameter :: Name='user_calc_sources'

    real :: cth
    real :: State_V(nVar)  

    real, dimension(nFluid) :: &
         Ux_I, Uy_I, Uz_I, U2_I, Temp_I, &
         UThS_I, URelS_I, URelSdim_I, UStar_I, Sigma_I, Rate_I, &
         UStarM_I, SigmaN_I, RateN_I, &
         I0xp_I, I0px_I, I2xp_I, I2px_I, &
         JxpUx_I, JxpUy_I, JxpUz_I, JpxUx_I, JpxUy_I, JpxUz_I, &
         Kxp_I, Kpx_I, Qepx_I, QmpxUx_I, QmpxUy_I, QmpxUz_I
    !For PU3
    real, dimension(nFluid) :: &
         URelSPu3_I, URelSPu3dim_I, UStarPu3_I, SigmaPu3_I, RatePu3_I, &
         UStarMPu3_I, SigmaNPu3_I, RateNPu3_I, &
         I0xpu3_I, I0pu3x_I, I2xpu3_I, I2pu3x_I, &
         Jxpu3Ux_I, Jxpu3Uy_I, Jxpu3Uz_I, Jpu3xUx_I, Jpu3xUy_I, Jpu3xUz_I, &
         Kxpu3_I, Kpu3x_I, Qepu3x_I, Qmpu3xUx_I, Qmpu3xUy_I, Qmpu3xUz_I

    logical:: UseSourceNe2P = .true., UseEnergySource = .true.

    integer :: i, j, k

    logical:: DoTest, DoTestMe
    character(len=*), parameter:: NameSub = 'user_calc_sources'
    !-----------------------------------------------------------------------
    ! Do not provide explicit source term when point-implicit scheme is used
    ! IsPointImplSource is true only when called from ModPointImplicit

    !As per Gabor's recommendation for explicit user sources

    !if(UsePointImplicit .and. .not. IsPointImplSource) RETURN
    if(IsPointImplSource) RETURN

    if(iBlock == BLKtest .and. iProc == PROCtest)then

       call set_oktest(NameSub, DoTest, DoTestMe)

    else
       DoTest = .false.; DoTestMe = .false.
    endif

    !  calculating some constants cBoltzmann is J/K 

    cth = 2.0*cBoltzmann/mNeutrals

    !ALL FLUIDS, including total ion

    ! Figure out which neutral population is produced at this point
    if(.not.IsPointImplPerturbed) call select_region(iBlock) 

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! Extract conservative variables
       State_V = State_VGB(:, i, j, k, iBlock)
       !Production rates of neutrals through charge exchange between sw ions and neutrals, 0 rate for sw ions with other ions

       !This calculates the array for flow values for all populations
       Ux_I  = State_V(iRhoUx_I)/State_V(iRho_I)
       Uy_I  = State_V(iRhoUy_I)/State_V(iRho_I)
       Uz_I  = State_V(iRhoUz_I)/State_V(iRho_I)

       ! Velocity square for the two ionized and four population of neutrals; 
       U2_I  = Ux_I**2 + Uy_I**2 + Uz_I**2


       ! Temperature for the two ionized and four population of neutrals (K)
       ! T = p/rho 
       ! since P_plasma=2T_proton*rho_ion; T_proton=0.5P_plasma/rho_ion
       Temp_I       = (State_V(iP_I)/State_V(iRho_I))*No2Si_V(UnitTemperature_)
       Temp_I(SWH_) = 0.5*Temp_I(SWH_)


       ! Thermal speed (squared) for ionized and three populations of neutrals
       ! UThS units are (m/s)^2
       UThS_I = cth*Temp_I !array of all thermal speeds here

       ! Relative velocity between neutrals and ionized fluid squared
       !! URelS_I = (Ux_I - Ux_I(1))**2 &
       !!     +    (Uy_I - Uy_I(1))**2 &
       !!     +    (Uz_I - Uz_I(1))**2 

       !For SW
       URelS_I(SWH_) =0.
       URelS_I(Pu3_) =0.

       URelS_I(Neu_) = (Ux_I(Neu_) - Ux_I(SWH_))**2 &
            +    (Uy_I(Neu_) - Uy_I(SWH_))**2 &
            +    (Uz_I(Neu_) - Uz_I(SWH_))**2

       URelS_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(SWH_))**2 &
            +    (Uy_I(Ne2_) - Uy_I(SWH_))**2 &
            +    (Uz_I(Ne2_) - Uz_I(SWH_))**2

       URelS_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(SWH_))**2 &
            +    (Uy_I(Ne3_) - Uy_I(SWH_))**2 &
            +    (Uz_I(Ne3_) - Uz_I(SWH_))**2

       URelS_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(SWH_))**2 &
            +    (Uy_I(Ne4_) - Uy_I(SWH_))**2 &
            +    (Uz_I(Ne4_) - Uz_I(SWH_))**2

       URelsdim_I(SWH_)=0.
       URelSdim_I(Pu3_)=0.
       URelSdim_I(Neu_)  = URelS_I(Neu_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne2_)  = URelS_I(Ne2_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne3_)  = URelS_I(Ne3_) * No2Si_V(UnitU_)**2
       URelSdim_I(Ne4_)  = URelS_I(Ne4_) * No2Si_V(UnitU_)**2

       !for PU3
       URelSPu3_I(SWH_) =0.
       URelSPu3_I(Pu3_) =0.

       URelSPu3_I(Neu_) = (Ux_I(Neu_) - Ux_I(PU3_))**2 &
            +    (Uy_I(Neu_) - Uy_I(PU3_))**2 &
            +    (Uz_I(Neu_) - Uz_I(PU3_))**2

       URelSPu3_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(PU3_))**2 &
            +    (Uy_I(Ne2_) - Uy_I(PU3_))**2 &
            +    (Uz_I(Ne2_) - Uz_I(PU3_))**2

       URelSPu3_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(PU3_))**2 &
            +    (Uy_I(Ne3_) - Uy_I(PU3_))**2 &
            +    (Uz_I(Ne3_) - Uz_I(PU3_))**2

       URelSPu3_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(PU3_))**2 &
            +    (Uy_I(Ne4_) - Uy_I(PU3_))**2 &
            +    (Uz_I(Ne4_) - Uz_I(Pu3_))**2

       URelSPu3dim_I(SWH_)=0.
       URelSPu3dim_I(Pu3_)=0.
       URelSPu3dim_I(Neu_)  = URelSPu3_I(Neu_) * No2Si_V(UnitU_)**2
       URelSPu3dim_I(Ne2_)  = URelSPu3_I(Ne2_) * No2Si_V(UnitU_)**2
       URelSPu3dim_I(Ne3_)  = URelSPu3_I(Ne3_) * No2Si_V(UnitU_)**2
       URelSPu3dim_I(Ne4_)  = URelSPu3_I(Ne4_) * No2Si_V(UnitU_)**2

       ! Calculating Cross Section Sigma_I for the different neutrals
       !
       ! Incorporating units to calculate the charge exchange cross sections
       ! No2Si_V(UnitU_) has units of m/s like cstartT so UReldim and UStar 
       ! has units of m/s

       !! URelSdim_I  = URelS_I * No2Si_V(UnitU_)**2

       !For SW
       UStar_I(SWH_)=0.
       UStar_I(Pu3_)=0.
       UStar_I(Neu_) =  sqrt(URelSdim_I(Neu_) + (4./cPi)*(UThS_I(Neu_) +UThS_I(SWH_)))
       UStar_I(Ne2_) =  sqrt(URelSdim_I(Ne2_) + (4./cPi)*(UThS_I(Ne2_) +UThS_I(SWH_)))
       UStar_I(Ne3_) =  sqrt(URelSdim_I(Ne3_) + (4./cPi)*(UThS_I(Ne3_) +UThS_I(SWH_)))
       UStar_I(Ne4_) =  sqrt(URelSdim_I(Ne4_) + (4./cPi)*(UThS_I(Ne4_) +UThS_I(SWH_)))

       ! UStar_I has units of m/s

       !!  UStar_I  = sqrt(URelSdim_I + (4./cPi)*(UThS_I +UThS_I(SWH_)))

       !For PU3 
       UStarPu3_I(SWH_)=0. 
       UStarPu3_I(Pu3_)=0. 
       UStarPu3_I(Neu_) =  sqrt(URelSPu3dim_I(Neu_) + (4./cPi)*(UThS_I(Neu_) +UThS_I(PU3_)))
       UStarPu3_I(Ne2_) =  sqrt(URelSPu3dim_I(Ne2_) + (4./cPi)*(UThS_I(Ne2_) +UThS_I(PU3_)))
       UStarPu3_I(Ne3_) =  sqrt(URelSPu3dim_I(Ne3_) + (4./cPi)*(UThS_I(Ne3_) +UThS_I(PU3_)))
       UStarPu3_I(Ne4_) =  sqrt(URelSPu3dim_I(Ne4_) + (4./cPi)*(UThS_I(Ne4_) +UThS_I(PU3_)))


       !!  UStar_I  = sqrt(URelSdim_I + (4./cPi)*(UThS_I +UThS_I(SWH_)))

       ! UStarM_I has units of m/s

       !!   UStarM_I  = sqrt(URelSdim_I + (64./(9.*cPi))*(UThS_I +UThS_I(SWH_)))


       !missing factor of 4/pi on the UthS_I speed?
       UStarM_I(SWH_)=0.
       UStarM_I(Pu3_)=0.
       UStarM_I(Neu_)  = sqrt(URelSdim_I(Neu_) + (64./(9.*cPi))*(UThS_I(Neu_) +UThS_I(SWH_)))
       UStarM_I(Ne2_)  = sqrt(URelSdim_I(Ne2_) + (64./(9.*cPi))*(UThS_I(Ne2_) +UThS_I(SWH_)))
       UStarM_I(Ne3_)  = sqrt(URelSdim_I(Ne3_) + (64./(9.*cPi))*(UThS_I(Ne3_) +UThS_I(SWH_)))
       UStarM_I(Ne4_)  = sqrt(URelSdim_I(Ne4_) + (64./(9.*cPi))*(UThS_I(Ne4_) +UThS_I(SWH_)))

       !For PU3

       UStarMPu3_I(SWH_)=0.
       UStarMPu3_I(Pu3_)=0.
       UStarMPu3_I(Neu_)  = sqrt(URelSPu3dim_I(Neu_) + (64./(9.*cPi))*(UThS_I(Neu_) +UThS_I(PU3_)))
       UStarMPu3_I(Ne2_)  = sqrt(URelSPu3dim_I(Ne2_) + (64./(9.*cPi))*(UThS_I(Ne2_) +UThS_I(PU3_)))
       UStarMPu3_I(Ne3_)  = sqrt(URelSPu3dim_I(Ne3_) + (64./(9.*cPi))*(UThS_I(Ne3_) +UThS_I(PU3_)))
       UStarMPu3_I(Ne4_)  = sqrt(URelSPu3dim_I(Ne4_) + (64./(9.*cPi))*(UThS_I(Ne4_) +UThS_I(PU3_)))


       ! Maher and Tinsley cross section Sigma 
       ! UStar has to have units of cm/s so the factor 100 is to pass m to cm
       ! Sigma has units of units of m^2  

       !for SW
       Sigma_I(SWH_)=0.
       SigmaN_I(SWH_)=0.
       Sigma_I(Pu3_)=0.
       SigmaN_I(Pu3_)=0.

       !For PU3
       SigmaPu3_I(SWH_)=0.
       SigmaNPu3_I(SWH_)=0.
       SigmaPu3_I(Pu3_)=0.
       SigmaNPu3_I(Pu3_)=0.


       ! Sigma_I(Neu_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Neu_)*100.))**2)*(1.E-4)
       !SigmaN_I(Neu_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Neu_)*100.))**2)*(1.E-4)
       !Sigma_I(Ne2_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne2_)*100.))**2)*(1.E-4)
       !SigmaN_I(Ne2_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne2_)*100.))**2)*(1.E-4)
       !Sigma_I(Ne3_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne3_)*100.))**2)*(1.E-4)
       !SigmaN_I(Ne3_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne3_)*100.))**2)*(1.E-4)
       !Sigma_I(Ne4_) =((1.64E-7 - (6.95E-9)*log(UStarM_I(Ne4_)*100.))**2)*(1.E-4)
       !SigmaN_I(Ne4_) =((1.64E-7 - (6.95E-9)*log(UStar_I(Ne4_)*100.))**2)*(1.E-4)


       ! New Cross Section from Lindsay and Stebbings, 2005
       !!Sigma_I =((2.2835E-7 - (1.062E-8)*log(UStarM_I*100.))**2)*(1.E-4)
       !!SigmaN_I =((2.2835E-7 - (1.062E-8)*log(UStar_I*100.))**2)*(1.E-4)
       !(4.15-0.531*log(E))^2(1-e^*67.5/E)^4.5 in 10^-16 cm^2, E in keV
       ! (1-e^a3/E) term ~ 1 at E<10keV
       !  log(v^2) = 2log(v), but what about the 2.2935E-7 term?  Should be 4.15 ??
       !some unit/dimension thing?
       Sigma_I(Neu_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Neu_)*100.))**2)*(1.E-4)
       SigmaN_I(Neu_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Neu_)*100.))**2)*(1.E-4)
       Sigma_I(Ne2_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne2_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne2_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne2_)*100.))**2)*(1.E-4)
       Sigma_I(Ne3_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne3_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne3_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne3_)*100.))**2)*(1.E-4)
       Sigma_I(Ne4_) =((2.2835E-7 - (1.062E-8)*log(UStarM_I(Ne4_)*100.))**2)*(1.E-4)
       SigmaN_I(Ne4_)=((2.2835E-7 - (1.062E-8)*log(UStar_I(Ne4_)*100.))**2)*(1.E-4)

       !For PU3
       SigmaPu3_I(Neu_) =((2.2835E-7 - (1.062E-8)*log(UStarMPu3_I(Neu_)*100.))**2)*(1.E-4)
       SigmaNPu3_I(Neu_)=((2.2835E-7 - (1.062E-8)*log(UStarPu3_I(Neu_)*100.))**2)*(1.E-4)
       SigmaPu3_I(Ne2_) =((2.2835E-7 - (1.062E-8)*log(UStarMPu3_I(Ne2_)*100.))**2)*(1.E-4)
       SigmaNPu3_I(Ne2_)=((2.2835E-7 - (1.062E-8)*log(UStarPu3_I(Ne2_)*100.))**2)*(1.E-4)
       SigmaPu3_I(Ne3_) =((2.2835E-7 - (1.062E-8)*log(UStarMPu3_I(Ne3_)*100.))**2)*(1.E-4)
       SigmaNPu3_I(Ne3_)=((2.2835E-7 - (1.062E-8)*log(UStarPu3_I(Ne3_)*100.))**2)*(1.E-4)
       SigmaPu3_I(Ne4_) =((2.2835E-7 - (1.062E-8)*log(UStarMPu3_I(Ne4_)*100.))**2)*(1.E-4)
       SigmaNPu3_I(Ne4_)=((2.2835E-7 - (1.062E-8)*log(UStarPu3_I(Ne4_)*100.))**2)*(1.E-4)





       ! Calculating Rate  = \nu * nH * mp where nH is the density of neutrals
       ! \nu = Sigma*np*u_star where np is the density of the ionized flow and 
       ! For each population of neutrals there will be another Rate
       ! The charge exhange cross section 100 to change ustar to cm/s
       ! Rate has no units (m^2*m/s*s*m-3 )

       !! Rate_I =Sigma_I*State_V(SWHRho_)*State_V(iRho_I)*UStarM_I  &
       !!      *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       !For SW  
       Rate_I(SWH_)=0.
       Rate_I(Pu3_)=0.

       Rate_I(Neu_) =Sigma_I(Neu_)*State_V(SWHRho_)*State_V(iRho_I(Neu_))*UStarM_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne2_) =Sigma_I(Ne2_)*State_V(SWHRho_)*State_V(iRho_I(Ne2_))*UStarM_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne3_) =Sigma_I(Ne3_)*State_V(SWHRho_)*State_V(iRho_I(Ne3_))*UStarM_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(Ne4_) =Sigma_I(Ne4_)*State_V(SWHRho_)*State_V(iRho_I(Ne4_))*UStarM_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       !For PU3      
       RatePu3_I(SWH_)=0.
       RatePu3_I(Pu3_)=0.

       RatePu3_I(Neu_) =SigmaPu3_I(Neu_)*State_V(PU3Rho_)*State_V(iRho_I(Neu_))*UStarMPu3_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RatePu3_I(Ne2_) =SigmaPu3_I(Ne2_)*State_V(PU3Rho_)*State_V(iRho_I(Ne2_))*UStarMPu3_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RatePu3_I(Ne3_) =SigmaPu3_I(Ne3_)*State_V(PU3Rho_)*State_V(iRho_I(Ne3_))*UStarMPu3_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RatePu3_I(Ne4_) =SigmaPu3_I(Ne4_)*State_V(PU3Rho_)*State_V(iRho_I(Ne4_))*UStarMPu3_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)


       !!RateN_I =SigmaN_I*State_V(Rho_)*State_V(iRho_I)*UStar_I  &
       !!     *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       !for SW
       RateN_I(SWH_)=0. 
       RateN_I(Pu3_)=0.     

       RateN_I(Neu_) =SigmaN_I(Neu_)*State_V(SWHRho_)*State_V(iRho_I(Neu_))*UStar_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne2_) =SigmaN_I(Ne2_)*State_V(SWHRho_)*State_V(iRho_I(Ne2_))*UStar_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne3_) =SigmaN_I(Ne3_)*State_V(SWHRho_)*State_V(iRho_I(Ne3_))*UStar_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateN_I(Ne4_) =SigmaN_I(Ne4_)*State_V(SWHRho_)*State_V(iRho_I(Ne4_))*UStar_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)

       !For PU3
       RateNPu3_I(SWH_)=0.
       RateNPu3_I(Pu3_)=0.

       RateNPu3_I(Neu_) =SigmaNPu3_I(Neu_)*State_V(PU3Rho_)*State_V(iRho_I(Neu_))*UStarPu3_I(Neu_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateNPu3_I(Ne2_) =SigmaNPu3_I(Ne2_)*State_V(PU3Rho_)*State_V(iRho_I(Ne2_))*UStarPu3_I(Ne2_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateNPu3_I(Ne3_) =SigmaNPu3_I(Ne3_)*State_V(PU3Rho_)*State_V(iRho_I(Ne3_))*UStarPu3_I(Ne3_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       RateNPu3_I(Ne4_) =SigmaNPu3_I(Ne4_)*State_V(PU3Rho_)*State_V(iRho_I(Ne4_))*UStarPu3_I(Ne4_)  &
            *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)



       ! Calculating the terms that enter in the Source terms
       ! The expressions for I0, Jxp, Kxp, Qexp are taken from Zank, Pauls, Williams, and Hall 1996
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
       ! 'xp' indicates neutrals-> SW protons charge exchange rates and
       ! 'px' indicates SW protons->neutrals
       ! 'xpu3' indicates neutrals -> PUI 3
       ! 'pu3x' indicates PUI3 -> neutrals
       ! charge exchange
       ! For example: 
       ! I0xpNe2 is the term of creation of Neu by charge exchange p-Ne2
       ! I0xpNe3 is the term of creation of Neu by charge exchange p-Ne3

       !For SW
       I0xp_I =RateN_I 
       I0px_I =RateN_I 
       !On 8/10 replaced UThS_I(1) with UThS_I(SWH_)
       I2xp_I  = Rate_I*(UStar_I/UStarM_I)*UThS_I(SWH_)/No2Si_V(UnitU_)**2    
       I2px_I  = Rate_I*(UStar_I/UStarM_I)*UThS_I/No2Si_V(UnitU_)**2


       !For PUI's
       I0xpu3_I = RateNPu3_I! no neutral -> pu3 for now, and this notation doesn't make much sense, would be another pu pop
       I0pu3x_I = RateNPu3_I
       I2xpu3_I =  RatePu3_I*(UStarPu3_I/UStarMPu3_I)*UThS_I(PU3_)/No2Si_V(UnitU_)**2
       I2pu3x_I = RatePu3_I*(UStarPu3_I/UStarMPu3_I)*UThS_I/No2Si_V(UnitU_)**2
       ! units are fine: (Uth2/ustar)*termxp is unitless as it should be
       !For SW  
       JxpUx_I  = Ux_I(SWH_)*Rate_I 
       JxpUy_I  = Uy_I(SWH_)*Rate_I   
       JxpUz_I  = Uz_I(SWH_)*Rate_I

       JpxUx_I  = Ux_I*Rate_I  
       JpxUy_I  = Uy_I*Rate_I   
       JpxUz_I  = Uz_I*Rate_I 

       !For PU3
       Jxpu3Ux_I  = Ux_I(PU3_)*RatePu3_I
       Jxpu3Uy_I  = Uy_I(PU3_)*RatePu3_I
       Jxpu3Uz_I  = Uz_I(PU3_)*RatePu3_I

       Jpu3xUx_I  = Ux_I*RatePu3_I
       Jpu3xUy_I  = Uy_I*RatePu3_I
       Jpu3xUz_I  = Uz_I*RatePu3_I

       !QmpxUx_I = JpxUx_I - JxpUx_I ! this is for neutrals, which can be created or destroyed
       !QmpxUy_I = JpxUy_I - JxpUy_I
       !QmpxUz_I = JpxUz_I - JxpUz_I

       !For SW
       QmpxUx_I(SWH_)=0.
       QmpxUy_I(SWH_)=0.
       QmpxUz_I(SWH_)=0.
       QmpxUx_I(Pu3_)=0.
       QmpxUy_I(Pu3_)=0.
       QmpxUz_I(Pu3_)=0.

       QmpxUx_I(Neu_) = (Ux_I(Neu_) - Ux_I(SWH_))*Rate_I(Neu_)
       QmpxUx_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(SWH_))*Rate_I(Ne2_)
       QmpxUx_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(SWH_))*Rate_I(Ne3_)
       QmpxUx_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(SWH_))*Rate_I(Ne4_)

       QmpxUy_I(Neu_) = (Uy_I(Neu_) - Uy_I(SWH_))*Rate_I(Neu_)
       QmpxUy_I(Ne2_) = (Uy_I(Ne2_) - Uy_I(SWH_))*Rate_I(Ne2_)
       QmpxUy_I(Ne3_) = (Uy_I(Ne3_) - Uy_I(SWH_))*Rate_I(Ne3_)
       QmpxUy_I(Ne4_) = (Uy_I(Ne4_) - Uy_I(SWH_))*Rate_I(Ne4_)

       QmpxUz_I(Neu_) = (Uz_I(Neu_) - Uz_I(SWH_))*Rate_I(Neu_)
       QmpxUz_I(Ne2_) = (Uz_I(Ne2_) - Uz_I(SWH_))*Rate_I(Ne2_)
       QmpxUz_I(Ne3_) = (Uz_I(Ne3_) - Uz_I(SWH_))*Rate_I(Ne3_)
       QmpxUz_I(Ne4_) = (Uz_I(Ne4_) - Uz_I(SWH_))*Rate_I(Ne4_)

       ! For PUI
       Qmpu3xUx_I(SWH_)=0.
       Qmpu3xUy_I(SWH_)=0.
       Qmpu3xUz_I(SWH_)=0.
       Qmpu3xUx_I(Pu3_)=0.
       Qmpu3xUy_I(Pu3_)=0.
       Qmpu3xUz_I(Pu3_)=0.

       Qmpu3xUx_I(Neu_) = (Ux_I(Neu_) - Ux_I(PU3_))*RatePu3_I(Neu_)
       Qmpu3xUx_I(Ne2_) = (Ux_I(Ne2_) - Ux_I(PU3_))*RatePu3_I(Ne2_)
       Qmpu3xUx_I(Ne3_) = (Ux_I(Ne3_) - Ux_I(PU3_))*RatePu3_I(Ne3_)
       Qmpu3xUx_I(Ne4_) = (Ux_I(Ne4_) - Ux_I(PU3_))*RatePu3_I(Ne4_)

       Qmpu3xUy_I(Neu_) = (Uy_I(Neu_) - Uy_I(PU3_))*RatePu3_I(Neu_)
       Qmpu3xUy_I(Ne2_) = (Uy_I(Ne2_) - Uy_I(PU3_))*RatePu3_I(Ne2_)
       Qmpu3xUy_I(Ne3_) = (Uy_I(Ne3_) - Uy_I(PU3_))*RatePu3_I(Ne3_)
       Qmpu3xUy_I(Ne4_) = (Uy_I(Ne4_) - Uy_I(PU3_))*RatePu3_I(Ne4_)

       Qmpu3xUz_I(Neu_) = (Uz_I(Neu_) - Uz_I(PU3_))*RatePu3_I(Neu_)
       Qmpu3xUz_I(Ne2_) = (Uz_I(Ne2_) - Uz_I(PU3_))*RatePu3_I(Ne2_)
       Qmpu3xUz_I(Ne3_) = (Uz_I(Ne3_) - Uz_I(PU3_))*RatePu3_I(Ne3_)
       Qmpu3xUz_I(Ne4_) = (Uz_I(Ne4_) - Uz_I(PU3_))*RatePu3_I(Ne4_)

       !For SW
       Kxp_I = 0.5*U2_I(SWH_)*Rate_I  + I2xp_I   

       Kpx_I = 0.5*U2_I*Rate_I  + I2px_I   

       Qepx_I = Kpx_I - Kxp_I 

       !For PU3
       Kxpu3_I = 0.5*U2_I(Pu3_)*RatePu3_I  + I2xpu3_I

       Kpu3x_I = 0.5*U2_I*RatePu3_I  + I2pu3x_I

       Qepu3x_I = Kpu3x_I - Kxpu3_I
       ! Calculate the source terms for this cell
       call calc_source_cell

    end do; end do; end do

  contains

    !==========================================================================
    subroutine calc_source_cell

      !C.P. needs to work here 

      ! Calculate source temrs for one cell. The pressures source is
      ! S(p) = (gamma-1)[S(e) - u.S(rhou) + 0.5 u**2 S(rho)]

      real:: Source_V(nVar + nFluid)
      integer:: iVar
      !---------------------------------------------------------------------

      Source_V = 0.0

      !on 8/10 changed some I0xp to I0px, which reflects anayltical equations better, but is mathematically indentical
      do iFluid = Neu_, Ne4_
         if(.not.UseSource_I(iFluid)) CYCLE
         call select_fluid
         !have losses to PUI, need sources, I think just a sign change on SW - source in zone 3
         if (iFluid == iFluidProduced_C(i,j,k)) then
            !iRho etc = change in density etc  
            Source_V(iRho)    = sum(I0xp_I(Neu_:Ne4_))  - I0px_I(iFluid) - I0pu3x_I(iFluid) !loss from interaction with PU3 
            Source_V(iRhoUx)  = sum(JxpUx_I(Neu_:Ne4_)) - JpxUx_I(iFluid) - Jpu3xUx_I(iFluid)
            Source_V(iRhoUy)  = sum(JxpUy_I(Neu_:Ne4_)) - JpxUy_I(iFluid) - Jpu3xUy_I(iFluid)
            Source_V(iRhoUz)  = sum(JxpUz_I(Neu_:Ne4_)) - JpxUz_I(iFluid) - Jpu3xUz_I(iFluid)
            Source_V(iEnergy) = sum(Kxp_I(Neu_:Ne4_))   - Kpx_I(iFluid) - Kpu3x_I(iFluid)
         else
            Source_V(iRho)    = - I0px_I(iFluid) - I0pu3x_I(iFluid) 
            Source_V(iRhoUx)  = - JpxUx_I(iFluid) - Jpu3xUx_I(iFluid)
            Source_V(iRhoUy)  = - JpxUy_I(iFluid) - Jpu3xUy_I(iFluid) 
            Source_V(iRhoUz)  = - JpxUz_I(iFluid) - Jpu3xUz_I(iFluid)
            Source_V(iEnergy) = - Kpx_I(iFluid) - Kpu3x_I(iFluid)
         end if
         Source_V(iP) = (g-1)* ( Source_V(iEnergy) &
              - Ux_I(iFluid)*Source_V(iRhoUx) &
              - Uy_I(iFluid)*Source_V(iRhoUy) &
              - Uz_I(iFluid)*Source_V(iRhoUz) &
              + 0.5*U2_I(iFluid)*Source_V(iRho) )
      end do

      !only make Pu3 in region before TS, and no loss of pu3
      if (Ne3_ == iFluidProduced_C(i,j,k)) then
         if(UseSource_I(SWH_) .and. UseSource_I(Pu3_))then
            Source_V(SWHRho_) = -sum( I0px_I(Neu_:Ne4_))  ! loss to PUIm symetric, so technically I)xp_I
            Source_V(SWHRhoUx_) = -sum(JxpUx_I(Neu_:Ne4_)) !QmpxUx_I(Neu_) + QmpxUx_I(Ne2_) + QmpxUx_I(Ne3_) + QmpxUx_I(Ne4_)
            Source_V(SWHRhoUy_) = -sum(JxpUy_I(Neu_:Ne4_)) !QmpxUy_I(Neu_) + QmpxUy_I(Ne2_) + QmpxUy_I(Ne3_) + QmpxUy_I(Ne4_)
            Source_V(SWHRhoUz_) = -sum(JxpUz_I(Neu_:Ne4_)) !QmpxUz_I(Neu_) + QmpxUz_I(Ne2_) + QmpxUz_I(Ne3_) + QmpxUz_I(Ne4_)

            Source_V(SWHEnergy_)= -sum( Kxp_I(Neu_:Ne4_) )

            !found mistake here on 7/26 and corrected, and on 10/6
            !added to pressure source term due to density source terms on 8/10
            Source_V(SWHp_) = (g-1)* ( Source_V(SWHEnergy_) &
                 - Uy_I(SWH_)*Source_V(SWHRhoUx_) &
                 - Uy_I(SWH_)*Source_V(SWHRhoUy_) &
                 - Uz_I(SWH_)*Source_V(SWHRhoUz_) &
                 + 0.5*U2_I(SWH_)*Source_V(SWHRho_) )

            !C.P. found mistake here on 7/27, with wrong J/K term, should be Jxp, not Jpx
            !added to pressure source term due to density source terms on 8/10
            Source_V(Pu3Rho_) = sum(I0xp_I(Neu_:Ne4_)) - sum(I0xpu3_I(Neu_:Ne4_))
            Source_V(Pu3RhoUx_) = sum(QmpxUx_I(Neu_:Ne4_)) + sum(JxpUx_I(Neu_:Ne4_))+sum(-Jxpu3Ux_I(Neu_:Ne4_))
            Source_V(Pu3RhoUy_) = sum(QmpxUy_I(Neu_:Ne4_)) + sum(JxpUy_I(Neu_:Ne4_))+sum(-Jxpu3Uy_I(Neu_:Ne4_))
            Source_V(Pu3RhoUz_) = sum(QmpxUz_I(Neu_:Ne4_)) + sum(JxpUz_I(Neu_:Ne4_))+sum(-Jxpu3Uz_I(Neu_:Ne4_))
            Source_V(Pu3Energy_)= sum(Qepx_I(Neu_:Ne4_)) + sum(Kxp_I(Neu_:Ne4_) ) +sum(-Kxpu3_I(Neu_:Ne4_))
            Source_V(Pu3p_) = (g-1)* ( Source_V(Pu3Energy_) &
                 - Ux_I(Pu3_)*Source_V(Pu3RhoUx_) &
                 - Uy_I(Pu3_)*Source_V(Pu3RhoUy_) &
                 - Uz_I(Pu3_)*Source_V(Pu3RhoUz_) &
                 + 0.5*U2_I(Pu3_)*Source_V(Pu3Rho_) ) 

         else if(UseSource_I(SWH_) .and. UseSource_I(PU3_))  then
            !this case isn't perfectly valid anymore, as SWH is not total ion fluid
            Source_V(SWHRho_) = 0.0 
            Source_V(SWHRhoUx_) = sum(QmpxUx_I(Neu_:Ne4_))
            Source_V(SWHRhoUy_) = sum(QmpxUy_I(Neu_:Ne4_)) 
            Source_V(SWHRhoUz_) = sum(QmpxUz_I(Neu_:Ne4_))
            !flipped from xp and -px to px and -xp on Oct 4 2011
            Source_V(SWHEnergy_)=  sum( +Kpx_I(Neu_:Ne4_) )+ sum( -Kxp_I(Neu_:Ne4_) )

            Source_V(SWHp_) = (g-1)* ( Source_V(SWHEnergy_) &
                 - Ux_I(SWH_)*Source_V(SWHRhoUx_) &
                 - Uy_I(SWH_)*Source_V(SWHRhoUy_) &
                 - Uz_I(SWH_)*Source_V(SWHRhoUz_) )

            Source_V(Pu3Rho_)   = sum(-I0xpu3_I(Neu_:Ne4_))
            Source_V(Pu3RhoUx_) = sum(-Jxpu3Ux_I(Neu_:Ne4_))
            Source_V(Pu3RhoUy_) = sum(-Jxpu3Uy_I(Neu_:Ne4_))
            Source_V(Pu3RhoUz_) = sum(-Jxpu3Uz_I(Neu_:Ne4_))
            Source_V(Pu3Energy_)= sum(-Kxpu3_I(Neu_:Ne4_))
            Source_V(Pu3p_) = (g-1)* ( Source_V(Pu3Energy_) &
                 - Ux_I(Pu3_)*Source_V(Pu3RhoUx_) &
                 - Uy_I(Pu3_)*Source_V(Pu3RhoUy_) &
                 - Uz_I(Pu3_)*Source_V(Pu3RhoUz_) &
                 + 0.5*U2_I(Pu3_)*Source_V(Pu3Rho_) )

         endif
         !otherwise don't do anything in this case,sources not on
      else if(UseSource_I(SWH_) .and. UseSource_I(PU3_)) then
         !Should this assumption be changed, right now acting like only SW interacts with neutrals outside of Region 3, ok approx?
         Source_V(SWHRho_) = 0.0
         Source_V(SWHRhoUx_) = sum(QmpxUx_I(Neu_:Ne4_))
         Source_V(SWHRhoUy_) = sum(QmpxUy_I(Neu_:Ne4_))
         Source_V(SWHRhoUz_) = sum(QmpxUz_I(Neu_:Ne4_))
         !fixed px and xp on Oct 4,2011
         Source_V(SWHEnergy_)=  sum( +Kpx_I(Neu_:Ne4_) )+ sum( -Kxp_I(Neu_:Ne4_) )

         Source_V(SWHp_) = (g-1)* ( Source_V(SWHEnergy_) &
              - Ux_I(SWH_)*Source_V(SWHRhoUx_) &
              - Uy_I(SWH_)*Source_V(SWHRhoUy_) &
              - Uz_I(SWH_)*Source_V(SWHRhoUz_) )

         Source_V(Pu3Rho_)   = sum(-I0xpu3_I(Neu_:Ne4_))
         Source_V(Pu3RhoUx_) = sum(-Jxpu3Ux_I(Neu_:Ne4_))
         Source_V(Pu3RhoUy_) = sum(-Jxpu3Uy_I(Neu_:Ne4_))
         Source_V(Pu3RhoUz_) = sum(-Jxpu3Uz_I(Neu_:Ne4_))
         Source_V(Pu3Energy_)= sum(-Kxpu3_I(Neu_:Ne4_))
         Source_V(Pu3p_) = (g-1)* ( Source_V(Pu3Energy_) &
              - Ux_I(Pu3_)*Source_V(Pu3RhoUx_) &
              - Uy_I(Pu3_)*Source_V(Pu3RhoUy_) &
              - Uz_I(Pu3_)*Source_V(Pu3RhoUz_) &
              + 0.5*U2_I(Pu3_)*Source_V(Pu3Rho_)) 

         !otherwise do nothing, no ion source on, pui produced  
      endif

      !Sum up sources for total ion fluid
      Source_V(Rho_) = Source_V(SWHRho_) +  Source_V(Pu3Rho_)
      Source_V(RhoUx_) = Source_V(SWHRhoUx_) +  Source_V(Pu3RhoUx_)
      Source_V(RhoUy_) =Source_V(SWHRhoUy_) +  Source_V(Pu3RhoUy_)
      Source_V(RhoUz_) = Source_V(SWHRhoUz_) +  Source_V(Pu3RhoUz_)
      Source_V(P_) = Source_V(SWHp_) + Source_V(PU3p_)

      Source_VC(:,i,j,k) = Source_VC(:,i,j,k) + Source_V

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
    real    :: InvRho, U2,SWHU2, p, Mach2, TempDim,TempDimSW, U2Dim, B2, rho, MachAlfven2,T_ave
    real    :: MachMagneto2
    !------------------------------------------------------------------------

    ! Produce fluid3 at the inner boundary

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       if(r_BLK(i,j,k,iBlock) < rPop3Limit) then
          iFluidProduced_C(i,j,k) = Ne3_
          CYCLE
       end if

       InvRho = 1.0/(State_VGB(SWHRho_,i,j,k,iBlock)+ State_VGB(Pu3Rho_,i,j,k,iBlock)) !rho is sum of ions density

       !U2     = sum((State_VGB(SWHRhoUx_:SWHRhoUz_,i,j,k,iBlock)+State_VGB(Pu3RhoUx_:Pu3RhoUz_,i,j,k,iBlock))**2)*InvRho**2 !+ added by C.P.


       p      = State_VGB(SWHp_,i,j,k,iBlock) + State_VGB(Pu3p_,i,j,k,iBlock)
       !U2     = (State_VGB(SWHRhoUx_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock) + State_VGB(Pu3RhoUx_,i,j,k,iBlock)/State_VGB(PU3Rho_,i,j,k,iBlock))**2 + &
       !         (State_VGB(SWHRhoUy_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock) + State_VGB(Pu3RhoUy_,i,j,k,iBlock)/State_VGB(PU3Rho_,i,j,k,iBlock))**2 + &
       !         (State_VGB(SWHRhoUz_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock) + State_VGB(Pu3RhoUz_,i,j,k,iBlock)/State_VGB(PU3Rho_,i,j,k,iBlock))**2 

       ! U2     = (sum(State_VGB(SWHRhoUx_:SWHRhoUz_,i,j,k,iBlock)**2)/State_VGB(SWHRho_,i,j,k,iBlock)**2 + sum(State_VGB(Pu3RhoUx_:Pu3RhoUz_,i,j,k,iBlock)**2)/State_VGB(PU3Rho_,i,j,k,iBlock)**2)/2  !/ number of fluids 

       U2 = (State_VGB(RhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**2 &
            + (State_VGB(RhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**2 &
            + (State_VGB(RhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock))**2

       SWHU2 = (State_VGB(SWHRhoUx_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock))**2 &
            + (State_VGB(SWHRhoUy_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock))**2 &
            + (State_VGB(SWHRhoUz_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock))**2

       U2Dim  = U2*No2Io_V(UnitU_)**2

       T_ave = (State_VGB(SWHp_,i,j,k,iBlock) + State_VGB(Pu3p_,i,j,k,iBlock)) &
            /(State_VGB(SWHRho_,i,j,k,iBlock) + State_VGB(Pu3Rho_,i,j,k,iBlock)) 



       !merav modifications
       rho = State_VGB(SWHrho_,i,j,k,iBlock) + State_VGB(Pu3rho_,i,j,k,iBlock) 
       B2 = sum(State_VGB(Bx_:Bz_,i,j,k,iBlock)**2)        

       ! Square of Alfven Mach Number
       MachAlfven2 = U2*rho/(B2+1.E-30)
       !MachMagneto2 = U2/((1.E-10)+(g*p*InvRho)+(B2*InvRho))

       MachMagneto2 = SWHU2/((1.E-10)+(g*T_ave)+(B2*InvRho)) !U2 or U2DIM?!!



       !end of merav modifications
       ! Square of Mach number
       !Mach2      = U2/(g*p*InvRho)
       Mach2       = U2/(g*T_ave)


       ! Temperature in Kelvins
       TempDim = InvRho*p*No2Si_V(UnitTemperature_)
       !Account for e- pressure
       TempDimSW = 1./2.*State_VGB(SWHp_,i,j,k,iBlock)/State_VGB(SWHRho_,i,j,k,iBlock)*No2Si_V(UnitTemperature_)
       ! if (TempPop1LimitDim > TempDimSW) write(*,*) &
       ! TempDimSW,State_VGB(SWHp_,i,j,k,iBlock),State_VGB(SWHRho_,i,j,k,iBlock), &
       ! No2Si_V(UnitTemperature_),r_BLK(i,j,k,iBlock), Xyz_DGB(:,i,j,k,iBlock), TempDim

!!!july12 use sonic Mach number - Berci
       if (MachPop4Limit**2 < Mach2 .and. uPop1LimitDim**2 > U2Dim) then
!!!october11          if (MachPop4Limit**2 < MachMagneto2 .and. uPop1LimitDim**2 > U2Dim) then
!!!       if (MachPop4Limit**2 < MachAlfven2 .and. uPop1LimitDim**2 > U2Dim .and. MachPop3Limit**2 < Mach2) then
          !Outside the bow shock
          iFluidProduced_C(i,j,k) = Ne4_
       elseif( TempPop1LimitDim > TempDim .and. uPop1LimitDim**2 > U2Dim)then
          ! elseif( TempPop1LimitDim > TempDimSW .and. uPop1LimitDim**2 > U2Dim)then
          !add spatial and checking for region 2 qualifications, Oct 26, 2011 
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
       !write(*,*) MachMagneto2, T_ave, No2Si_V(UnitTemperature_), B2, InvRho, &
       ! U2, U2Dim, No2Io_V(UnitU_), Xyz_DGB(:,i,j,k,iBlock),r_BLK(i,j,k,iBlock), &
       ! iFluidProduced_C(i,j,k)

    end do; end do; end do

  end subroutine select_region

  !============================================================================

  subroutine user_init_point_implicit

    use ModVarIndexes, ONLY: &
         SWHRho_,SWHRhoUx_,SWHRhoUy_,SWHRhoUz_,SWHP_, &
         Pu3Rho_, Pu3RhoUx_, Pu3RhoUy_, Pu3RhoUz_, Pu3P_, & 
         NeuRho_,NeuRhoUx_,NeuRhoUy_, NeuRhoUz_,NeuP_,&
         Ne2Rho_,Ne2RhoUx_,Ne2RhoUy_,Ne2RhoUz_,Ne2P_,&
         Ne3Rho_,Ne3RhoUx_,Ne3RhoUy_,Ne3RhoUz_,Ne3P_,&
         Ne4Rho_,Ne4RhoUx_,Ne4RhoUy_,Ne4RhoUz_,Ne4P_


    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, EpsPointImpl_V
    !------------------------------------------------------------------------

    ! Allocate and set iVarPointImpl_I
    ! In this example there are 3 implicit variables


    ! All the neutrals momenta and plasma are implicit 
    ! (4 neutral fluid and 2 ion)

    !commented out to avoid user implicit terms
    !allocate(iVarPointImpl_I(20))

    ! iVarPointImpl_I = (/SWHRho_,SWHRhoUx_, SWHRhoUy_, SWHRhoUz_,SWHP_, &
    !     Pu3Rho_, Pu3RhoUx_, Pu3RhoUy_, Pu3RhoUz_, Pu3P_, &
    !     NeuRho_, NeuRhoUx_, NeuRhoUy_, NeuRhoUz_, NeuP_, &
    !     Ne2Rho_, Ne2RhoUx_, Ne2RhoUy_, Ne2RhoUz_, Ne2P_, &
    !     Ne3Rho_, Ne3RhoUx_, Ne3RhoUy_, Ne3RhoUz_, Ne3P_, &
    !     Ne4Rho_, Ne4RhoUx_, Ne4RhoUy_, Ne4RhoUz_, Ne4P_ /)


    ! Because the second and third populations of neutrals have initially 
    ! very small values I'm
    ! setting EpsPointImpl_V to be small for these variables ??? !!!
    !EpsPointImpl_V(Ne2Rho_:Ne4P_) = 1.e-11

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set.

    !IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

end module ModUser
