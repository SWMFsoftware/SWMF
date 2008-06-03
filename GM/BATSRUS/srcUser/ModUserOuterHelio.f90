!April 30, 2007 implementing MultiFluid
!June 01, 2007 correcting normalization
!June 08, 2007 correcting source terms
!August 18, 2007 Implementing 3-fluids
!October 23, 2007 more little corrections
!January 01, 2008 Source terms in point-implicit form - with help
! of Gabor Toth
!==============================================================================
module ModUser

  use ModSize,     ONLY: nI,nJ,nK,gcn,nBLK
  use ModMain
  use ModPhysics
  use ModSetOuterBC
  use ModAdvance,  ONLY : State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK, far_field_BCs_BLK, MaxBoundary
  use ModVarIndexes
  use ModConst, ONLY: cLightSpeed, cElectronCharge, cElectronMass
  use ModProcMH
  use ModMultiFluid
  use ModUserEmpty,               &
       IMPLEMENTED1 => user_read_inputs,                &
       IMPLEMENTED2 => user_face_bcs,                   &
       IMPLEMENTED3 => user_normalization,               &
       IMPLEMENTED4 => user_set_outerbcs,               &
       IMPLEMENTED5 => user_set_ics,                    &
       IMPLEMENTED8 => user_write_progress,            & 
       IMPLEMENTED9 => user_io_units,                  &
       IMPLEMENTED10 => user_set_plot_var,              &
       IMPLEMENTED11 => user_calc_sources,              &
       IMPLEMENTED12 => user_init_point_implicit


  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 2.0
  ! I am calling Version Module = 2.0 the version with 3-fluids
  character (len=*), parameter :: NameUserModule = 'Global Heliosphere'

  real :: rFriction = 1.0

  ! 
  ! SWH variables.
  !/
  real ::      SWH_T_dim=0.0  , &
       SWH_a_dim=0.0  , &
       SWH_rho=0.0 , SWH_rho_dim=0.0, &
       SWH_p=0.0  , SWH_p_dim=0.0   , &
       SWH_Ux=0.0 , SWH_Ux_dim=0.0 , &
       SWH_Uy=0.0  , SWH_Uy_dim=0.0 , &
       SWH_Uz=0.0 , SWH_Uz_dim=0.0 , &
       SWH_Bx=0.0 , SWH_Bx_dim=0.0 , &
       SWH_By=0.0 , SWH_By_dim=0.0 , &
       SWH_Bz=0.0 , SWH_Bz_dim=0.0 , &
       SWH_B_factor=0.0

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

  real, dimension(0:1) :: &
       VLISW_rho_t,  &
       VLISW_p_t  ,  &
       VLISW_Ux_t ,  &
       VLISW_Uy_t ,  &
       VLISW_Uz_t ,  &
       VLISW_Bx_t ,  &
       VLISW_By_t ,  &
       VLISW_Bz_t ,  &
       VLISW_time_t
  !
  ! FASTSW variables.
  !/
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

  real, dimension(0:1) :: &
       SWfast_rho_t,  &
       SWfast_p_t  ,  &
       SWfast_Ux_t ,  &
       SWfast_Uy_t ,  &
       SWfast_Uz_t ,  &
       SWfast_Bx_t,   &
       SWfast_By_t ,  &
       SWfast_Bz_t ,  &
       SWfast_time_t
  !
  ! neutrals variables
  !/

  real :: TNeutralsISW_dim=0.0, &
       RhoNeutralsISW=0.0, RhoNeutralsISW_dim=0.0 , &
       PNeutralsISW=0.0  , PNeutralsISW_dim=0.0  , &
       UxNeutralsISW=0.0 , UxNeutralsISW_dim=0.0 , &
       UyNeutralsISW=0.0 , UyNeutralsISW_dim=0.0 , &
       UzNeutralsISW=0.0 , UzNeutralsISW_dim=0.0 ,  &
       mNeutralsmp, mNeutrals

  real, dimension(0:1) :: &
       RhoNeutralsISW_t,  &
       PneutralsISW_t  ,  &
       UxNeutralsISW_t ,  &
       UyNeutralsISW_t ,  &
       UzNeutralsISW_t

contains

  !=========================================================================

  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam

    character (len=100) :: NameCommand
    character (len=*), parameter :: Name='user_read_inputs'

    integer:: i, j, k, n, m

    !-------------------------------------------------------------------
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)

       case('#USERINPUTEND')
          if(iProc==0) write(*,*)'USERINPUTEND'
          EXIT
       case("#SOLARWINDH")
          call read_var('SWH_rho_dim' ,SWH_rho_dim)
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
          ! Probably in this case i will define rFriction as larger than the domain
       case("#FRICTION")
          call read_var('rFriction',rFriction)

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !=========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,     ONLY: nDim,West_,North_,Top_
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModPhysics  
    use ModNumConst, ONLY: cZero,cOne,cTwo,cTolerance
    use ModSetOuterBC
    use ModProcMH
    use ModFaceBc, ONLY: iBoundary, FaceCoords_D, VarsTrueFace_V
    !
    real, intent(out):: VarsGhostFace_V(nVar)

    ! Variables required by this user subroutine
    !/
    real:: XFace,YFace,ZFace
    real:: VxFaceOutside,VyFaceOutside,VzFaceOutside
    real:: BxFaceOutside,ByFaceOutside,BzFaceOutside
    real:: VrFaceOutside,VthetaFaceOutside,VphiFaceOutside,&
         VrFaceInside,VthetaFaceInside,VphiFaceInside,     &
         BrFaceOutside,BthetaFaceOutside,BphiFaceOutside,  &
         BrFaceInside,BthetaFaceInside,BphiFaceInside
    real :: RhoFaceInside, PFaceInside
    real:: cosTheta,sinTheta,cosPhi,sinPhi,RFace
    real, dimension(1:3):: location,v_phi
    real:: cosThetaT,sinThetaT,cosPhiT,sinPhiT
    ! Arguments
    !
    ! Global Heliospheric Variables
    !/
    real :: VrSolarWind,VthetaSolarWind,VphiSolarWind,BrSolarWind
    real :: BthetaSolarWind,BphiSolarWind,RhoSolarWind,PSolarWind
    real :: sg,BrSolar,BthetaSolar,BphiSolar,B0mag,Btot,Pmag,VphiSolar
    real :: sin2Theta_fast_wind, PmagEquator
    real::  rot_period_dim,OMEGAbodyH
    real :: PSolarWindT
    !/
    !  integer,intent(in)::iBLK

    !-------------------------------------------------------------------

    No2Si_V(UnitX_)= 215.0*Rsun                                ! m
    No2Si_V(UnitU_)= sqrt(g*cBoltzmann*SWH_T_dim/cProtonMass)
    No2Si_V(UnitRho_)=cProtonMass*SWH_rho_dim*1.0E+6           ! kg/m^

    XFace = FaceCoords_D(1)
    YFace = FaceCoords_D(2)
    ZFace = FaceCoords_D(3)
    VxFaceOutside = VarsTrueFace_V(Ux_)
    VyFaceOutside = VarsTrueFace_V(Uy_)
    VzFaceOutside = VarsTrueFace_V(Uz_)
    BxFaceOutside = VarsTrueFace_V(Bx_)
    ByFaceOutside = VarsTrueFace_V(By_)
    BzFaceOutside = VarsTrueFace_V(Bz_)

    RFace = sqrt(XFace**2 + YFace**2 + ZFace**2)
    !
    ! Rotate to spherical coordinates
    !/

    !Apply boundary conditions

    select case(iBoundary)
    case(body1_)
       cosTheta = ZFace/RFace
       sinTheta = sqrt(XFace**2+YFace**2)/RFace
       cosPhi = XFace/sqrt(XFace**2+YFace**2+cTolerance**2)
       sinPhi = YFace/sqrt(XFace**2+YFace**2+cTolerance**2)
       VrFaceOutside = (VxFaceOutside*XFace + &
            VyFaceOutside*YFace + &
            VzFaceOutside*ZFace) / RFace
       VthetaFaceOutside = ((VxFaceOutside*XFace + &
            VyFaceOutside*YFace) * ZFace - &
            VzFaceOutside*(XFace**2+YFace**2)) / &
            (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
       VphiFaceOutside = (VyFaceOutside*XFace - &
            VxFaceOutside*YFace)*sinTheta / &
            ((XFace**2+YFace**2+cTolerance**2)/RFace)
       BrFaceOutside = (BxFaceOutside*XFace + &
            ByFaceOutside*YFace + &
            BzFaceOutside*ZFace) / RFace
       BthetaFaceOutside = ((BxFaceOutside*XFace + &
            ByFaceOutside*YFace) * ZFace - &
            BzFaceOutside*(XFace**2+YFace**2)) / &
            (sqrt(XFace**2+YFace**2+cTolerance**2)*RFace)
       BphiFaceOutside = (ByFaceOutside*XFace - &
            BxFaceOutside*YFace)*sinTheta / &
            ((XFace**2+YFace**2+cTolerance**2)/RFace)

       ! defining the rotation components of the Sun

       rot_period_dim = 26.0*24.0     ! rotation period in hours
       OMEGAbodyH = (2.00*cPi/(rot_period_dim*3600.00)) !in sec^-1
       Rbody = 30. !define it

       !
       ! calculating the parker field components, Br, Btheta and Bphi
       ! SW_Bx is the value of the field in the pole B0
       !/

       sg = sign(1.0,cosTheta)
       !     sg = -sign(1.0,cosTheta)
       ! should I invert the polarity for this cycle?

       B0mag       =  SWH_Bx
       BrSolar     =  sg*B0mag

       ! the factor 1.496E8 is to convert Rbody in AU to km; 

       BphiSolar   = -(sg*B0mag*OMEGAbodyH*Rbody*(1.496E8)*sinTheta)/(SWH_Ux_dim)
       BthetaSolar =  0.0

       !      VphiSolar = OMEGAbody*(6.96E5)*sinTheta/unitUSER_U
       ! Vphi=omega*r*sinTheta - 6.96E5 is the solar radii in km
       ! No2Si_V(UnitU_) is in m/s so its ok
       ! No2Io_V(UnitU_) is in km/s
       ! factor 1000 for test

       VphiSolar = OMEGAbodyH*(6.96E5)*sinTheta/No2Io_V(UnitU_)

       Btot = sqrt((BrSolar**2)+(BphiSolar**2))

       ! magnetic pressure Btot*unitUSER_B will be in nT  and the pressure in Pa
       ! No2Si_V(UnitB_) is in T

       Pmag = (((Btot*No2Io_V(UnitB_)*1.0E-9)**2)/(2.0*cMu))/(inv_g*0.1*No2Io_V(UnitP_))

       !magnetic pressure in the equator
       !
       !     PmagEquator = (((SWH_Bx*27.9*unitUSER_B*1.0E-
       !9)**2)/(2.0*cMu))/(inv_g*unitSI_p)

       ! The factor 27.9 comes from Rbody*OMEGAbody/SW_Ux
       !/
       !
       PmagEquator = (((SWH_Bx*27.9*No2Io_V(UnitB_)*1.0E-9)**2)/(2.0*cMu))/(inv_g*0.1*No2Io_V(UnitP_))

       ! Introducing the Fast Solar Wind
!!!    sin2Theta_fast_wind = 0.250000    ! 30 degrees
       !  sin2Theta_fast_wind = 0.1786062   ! 25 degrees
       !  sin2Theta_fast_wind = 0.116980    ! 20 degrees
       !  sin2Theta_fast_wind = 1.000000    ! 90 degrees

!!!      if (sinTheta*sinTheta > sin2Theta_fast_wind) then
       !SLOW WIND
       VrSolarWind     = SWH_Ux
       VthetaSolarWind = 0.0
       VphiSolarWind   = VphiSolar
       BrSolarWind     = BrSolar
       BthetaSolarWind = BthetaSolar
       BphiSolarWind   = BphiSolar
       RhoSolarWind    = SWH_rho
!!!       else
!!!          ! FAST WIND
!!!          VrSolarWind     =  SWfast_Ux

!!!          VthetaSolarWind =  0.0

!!!          VphiSolarWind   =  VphiSolar

!!!          BrSolarWind     =  BrSolar

!!!          BthetaSolarWind =  BthetaSolar

!!!          BphiSolarWind   =  BphiSolar

!!!          RhoSolarWind    =  SWfast_rho

!!!      end if

       ! Latitude variating wind


!!!         VrSolarWind = SWH_Ux+((SWfast_Ux - SWH_Ux)*cosTheta*cosTheta)

!!!         RhoSolarWind = SWH_rho + ((SWfast_rho - SWH_rho)*cosTheta*cosTheta)

       ! pressure will vary for fast and slow wind         PSolarWind = SWH_p

       !
       !
       !/
!!!test      VrSolarWind = SWH_Ux

       VthetaSolarWind = 0.0
       VphiSolarWind = VphiSolar
       BrSolarWind = BrSolar
       BthetaSolarWind = BthetaSolar
       BphiSolarWind = BphiSolar

       ! making sure that the total pressure Pt+Pb will be constant (PmagEquator is Pb 

       !at the equator)
       !/
!!!               PSolarWindT = SWH_p + (SWH_p - SWfast_p)*cosTheta*cosTheta

!!!               PsolarWind  = PSolarWindT + PmagEquator - Pmag

       PSolarWind = SWH_p + PmagEquator -  Pmag

       ! ISW
       !
       VrFaceInside      =   VrSolarWind
       VthetaFaceInside =    VthetaSolarWind
       VphiFaceInside   =    VphiSolarWind
       BrFaceInside      =   BrSolarWind
       BthetaFaceInside  =   BthetaSolarWind
       BphiFaceInside    =   BphiSolarWind
       RhoFaceInside     =   RhoSolarWind
       PFaceInside       =   PSolarWind

       VarsGhostFace_V(rho_) = RhoFaceInside
       VarsGhostFace_V(p_) = PFaceInside

       ! Rotate back to cartesian coordinates
       VarsGhostFace_V(Ux_) = VrFaceInside*XFace/RFace+&
            VthetaFaceInside*cosTheta*cosPhi          -&
            VphiFaceInside*sinPhi
       VarsGhostFace_V(Uy_) = VrFaceInside*YFace/RFace+&
            VthetaFaceInside*cosTheta*sinPhi          +&
            VphiFaceInside*cosPhi
       VarsGhostFace_V(Uz_) = VrFaceInside*ZFace/RFace-&
            VthetaFaceInside*sinTheta
       VarsGhostFace_V(Bx_) = BrFaceInside*XFace/RFace+&
            BthetaFaceInside*cosTheta*cosPhi          -&
            BphiFaceInside*sinPhi
       VarsGhostFace_V(By_) = BrFaceInside*YFace/RFace+&
            BthetaFaceInside*cosTheta*sinPhi          +&
            BphiFaceInside*cosPhi
       VarsGhostFace_V(Bz_) = BrFaceInside*ZFace/RFace-&
            BthetaFaceInside*sinTheta
       !attempt to deal with neutrals

       !!    VarsGhostFace_V(NeuRho_) = RhoNeutralsISW
       !!    VarsGhostFace_V(NeuRhoUx_)= 0.0
       !!    VarsGhostFace_V(NeuRhoUy_)= 0.0
       !!    VarsGhostFace_V(NeuRhoUz_)= 0.0
       !!    VarsGhostFace_V(NeuP_)= PNeutralsISW

       ! NeuRho is PopI; NeuIIRho is PopII and NeuIIIRho is PopIII
       !
       ! Pop I
       !/  
       VarsGhostFace_V(NeuRho_) = VarsTrueFace_V(NeuRho_)
       VarsGhostFace_V(NeuRhoUx_) = VarsTrueFace_V(NeuRhoUx_)
       VarsGhostFace_V(NeuRhoUy_) = VarsTrueFace_V(NeuRhoUy_)
       VarsGhostFace_V(NeuRhoUz_) = VarsTrueFace_V(NeuRhoUz_)
       VarsGhostFace_V(NeuP_) = VarsTrueFace_V(NeuP_)
       !
       ! Pop II
       !/  
       VarsGhostFace_V(Ne2Rho_) = VarsTrueFace_V(Ne2Rho_) 
       VarsGhostFace_V(Ne2RhoUx_) = VarsTrueFace_V(Ne2RhoUx_)
       VarsGhostFace_V(Ne2RhoUy_) = VarsTrueFace_V(Ne2RhoUy_)
       VarsGhostFace_V(Ne2RhoUz_) = VarsTrueFace_V(Ne2RhoUz_)
       VarsGhostFace_V(Ne2P_) = VarsTrueFace_V(Ne2P_)
       !
       ! Pop III
       !/  
       VarsGhostFace_V(Ne3Rho_) = VarsTrueFace_V(Ne3Rho_) 
       VarsGhostFace_V(Ne3RhoUx_) = VarsTrueFace_V(Ne3RhoUx_)
       VarsGhostFace_V(Ne3RhoUy_) = VarsTrueFace_V(Ne3RhoUy_)
       VarsGhostFace_V(Ne3RhoUz_) = VarsTrueFace_V(Ne3RhoUz_)
       VarsGhostFace_V(Ne3P_) = VarsTrueFace_V(Ne3P_)

    end select
    !-------------------------------------------------------------------

  end subroutine user_face_bcs

  !============================================================================
  subroutine user_normalization

    use ModProcMH, ONLY:iProc  
    use ModMain 
    use ModPhysics 
    use ModVarIndexes       
    character (len=*), parameter :: Name='user_normalization'
    !-------------------------------------------------------------------

    SWH_rho_dim = 7.8E-3 !n/cm^3
    SWH_T_dim=1.609E3   !K
    No2Si_V(UnitX_)= 215.0*Rsun                                ! m
    No2Si_V(UnitU_)= sqrt(g*cBoltzmann*SWH_T_dim/cProtonMass)
    No2Si_V(UnitRho_)=cProtonMass*SWH_rho_dim*1.0E+6           ! kg/m^3

  end subroutine user_normalization

  !============================================================================

  subroutine user_set_outerbcs(iBlock,iSide,TypeBc,found)

    use ModMain
    use ModVarIndexes
    use ModAdvance,   ONLY : State_VGB
    use ModPhysics,   ONLY : CellState_VI
    use ModSetOuterBC
    use ModProcMH
    use ModAdvance, ONLY : State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModMultiFluid

    integer,intent(in)::iBlock, iSide
    logical,intent(out) :: found
    character (len=20),intent(in) :: TypeBc
    real :: time_now

    found = .true.

    State_VGB(rho_,:,:,:,iBlock)=VLISW_rho
    State_VGB(rhoUx_,:,:,:,iBlock)=VLISW_rho*VLISW_Ux
    State_VGB(rhoUy_,:,:,:,iBlock)=VLISW_rho*VLISW_Uy
    State_VGB(rhoUz_,:,:,:,iBlock)=VLISW_rho*VLISW_Uz
    State_VGB(Bx_,:,:,:,iBlock)=VLISW_Bx
    State_VGB(By_,:,:,:,iBlock)=VLISW_By
    State_VGB(Bz_,:,:,:,iBlock)=VLISW_Bz
    State_VGB(p_,:,:,:,iBlock)=VLISW_p
    !
    ! PopI
    !/    
    State_VGB(NeuRho_,:,:,:,iBlock)=RhoNeutralsISW    
    State_VGB(NeuRhoUx_,:,:,:,iBlock)=RhoNeutralsISW*UxNeutralsISW    
    State_VGB(NeuRhoUy_,:,:,:,iBlock)=RhoNeutralsISW*UyNeutralsISW   
    State_VGB(NeuRhoUz_,:,:,:,iBlock)=RhoNeutralsISW*UzNeutralsISW    
    State_VGB(NeuP_,:,:,:,iBlock)=PNeutralsISW
    !
    ! PopII
    !/ 
    !we want for NeuIIRho outflow

    State_VGB(Ne2Rho_,:,:,:,iBlock)=RhoNeutralsISW*cTiny
    State_VGB(Ne2RhoUx_,:,:,:,iBlock)=RhoNeutralsISW*UxNeutralsISW*cTiny
    State_VGB(Ne2RhoUy_,:,:,:,iBlock)=0.0
    State_VGB(Ne2RhoUz_,:,:,:,iBlock)=0.0
    State_VGB(Ne2P_,:,:,:,iBlock)=PNeutralsISW*cTiny

    !
    !    call BC_cont(Ne2Rho_,Ne2p_)
    ! erro discover Nov 11, 2007
    !
    ! PopIII
    !/ 

    call BC_cont(Ne3Rho_,Ne3p_)

    !we want for NeuIII outflow
  end subroutine user_set_outerbcs

  !=====================================================================

  subroutine user_set_ics    
    use ModMain,  ONLY: globalBLK    
    use ModGeometry, ONLY: r_BLK    
    use ModVarIndexes    
    use ModAdvance,  ONLY: State_VGB    
    use ModPhysics,  ONLY: rBody, BodyRho_I, BodyP_I    

    implicit none    
    integer :: iBlock    

    integer :: i,j,k

    real :: x0,y0,z0,rho0
    real, dimension(3) :: B0, vel
    real :: rp,Br,Btheta,Bphi
    real :: B0mag,VrS,VthetaS,VphiS
    real :: cos_theta,sin_theta,cos_phi,sin_phi,sg
    real :: rot_period_dim, OMEGAbodyH
    real :: thetaN, sinthetaN, lambda, RhoSolarW
    real :: sin2Theta_fast_wind

    !--------------------------------------------------------------------------

    iBlock = globalBLK
    rot_period_dim = 26.0*24.0     ! rotation period in hours
    OMEGAbodyH = (2.00*cPi/(rot_period_dim*3600.00))
    Rbody = 30. !define it

    SWH_rho_dim = 7.8E-3 !n/cm^3      
    SWH_T_dim=1.609E3   !K      
    No2Si_V(UnitX_)= 215.0*Rsun    ! m      
    No2Si_V(UnitU_)= sqrt(g*cBoltzmann*SWH_T_dim/cProtonMass)      
    No2Si_V(UnitRho_)=cProtonMass*SWH_rho_dim*1.0E+6             ! kg/m^3 

    do i=1-gcn,nI+gcn
       do j=1-gcn,nJ+gcn
          do k=1-gcn,nK+gcn

             x0 = x_BLK(i,j,k,iBlock)
             y0 = y_BLK(i,j,k,iBlock)
             z0 = z_BLK(i,j,k,iBlock)

             ! defining the rotation components of the Sun

             ! transforming x0,y0,z0 to r, theta and phi
             ! theta is the angle measured from the pole

             rp=sqrt(x0**2+y0**2+z0**2)

             cos_theta=z0/rp
             sin_theta=sqrt(x0**2+y0**2)/rp
             cos_phi=x0/sqrt(x0**2+y0**2)
             sin_phi=y0/sqrt(x0**2+y0**2)

             ! for the neutrals
             thetaN = atan(sqrt(y0**2+z0**2)/(x0+cTiny))
             sinThetaN=sqrt(y0**2+z0**2)/rp
             lambda = 4.0
             !so pra iniciar eu fixei lambda como 4.0 (olha a expressao 3.53 da tese do Linde)


             !calculating the parker field components, Br, Btheta and Bphi


             sg=sign(1.0,cos_theta)
             !       sg=-sign(1.0,cos_theta)
             B0mag  =  SWH_Bx
             Br     =  sg*B0mag*((Rbody/rp)**2)
             Bphi   = -((sg*B0mag*OMEGAbodyH*Rbody*(1.496E8)*sin_theta)/(SWH_Ux_dim))*(Rbody/rp)
             Btheta =  0.0

             ! Vphi =omega*Rs*sintheta*(Rs/r)

             VphiS = OMEGAbodyH*(6.96E5)*sin_theta/No2Io_V(UnitU_)
             VrS = SWH_Ux
             VthetaS = 0.

             ! Introducing the Fast Solar Wind
             !!    sin2Theta_fast_wind = 0.250000    ! 30 degrees
             !  sin2Theta_fast_wind = 0.1786062   ! 25 degrees
             !  sin2Theta_fast_wind = 0.116980    ! 20 degrees
             !  sin2Theta_fast_wind = 1.000000    ! 90 degrees

             !!      if (sin_theta*sin_theta > sin2Theta_fast_wind) then
             !!          !SLOW WIND
             !!          VrS     = SWH_Ux
             !!          RhoSolarW    = SWH_rho
             !!       else
             !!          ! FAST WIND
             !!          VrS     =  SWfast_Ux
             !!          RhoSolarW    =  SWfast_rho
             !!      end if

             !!I still need to impemenet that the density vary from solar wind slow and fast

             !!
             if (rp>Rbody) then
                ! magnetic field components in cartesian coordinates    
                ! x component of magnetic field  
                ! magnetic field components in cartesian coordinates
                !x component of magnetic field
                B0(1)=Br    *sin_theta*cos_phi  + &
                     Btheta*cos_theta*cos_phi - &
                     Bphi  *sin_phi
                !y component of magnetic field
                B0(2)=Br    *sin_theta*sin_phi  + &
                     Btheta*cos_theta*sin_phi + &
                     Bphi  *cos_phi
                !z component of magnetic field
                B0(3)=Br    *cos_theta          - &
                     Btheta*sin_theta

                VphiS = 0.      !NECESSARY TO WORK 

                State_VGB(Bx_,i,j,k,iBlock) =  B0(1) 
                State_VGB(By_,i,j,k,iBlock) =  B0(2)
                State_VGB(Bz_,i,j,k,iBlock) =  B0(3)

                !velocity components in cartesian coordinates
                !x component of magnetic field
                vel(1)=VrS    *(x0/rp)  + &
                     VthetaS*cos_theta*cos_phi - &
                     VphiS  *sin_phi
                !y component of magnetic field
                vel(2)=VrS    *(y0/rp)  + &
                     VthetaS*cos_theta*sin_phi + &
                     VphiS  *cos_phi

                !z component of magnetic field
                vel(3)=VrS    *(z0/rp)  - &
                     VthetaS*sin_theta

                ! density and pressure
                State_VGB(rho_,i,j,k,iBlock) = SWH_rho*((Rbody/rp)**2)
                State_VGB(P_,i,j,k,iBlock) = SWH_p*((Rbody/rp)**(2*g))
                State_VGB(rhoUx_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock)*vel(1)
                State_VGB(rhoUy_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock)*vel(2)
                State_VGB(rhoUz_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock)*vel(3) 
                !
                ! PopI
                !/
                State_VGB(NeuRho_,i,j,k,iBlock) = RhoNeutralsISW*Exp(-lambda*thetaN/(rp*SinThetaN))
                State_VGB(NeuP_,i,j,k,iBlock) = PNeutralsISW*Exp(-lambda*thetaN/(rp*SinThetaN))
                State_VGB(NeuRhoUx_,i,j,k,iBlock) = UxNeutralsISW*RhoNeutralsISW*Exp(-lambda*thetaN/(rp*SinThetaN)) 

                !
                ! PopII - o que fazer com set_ICs popII?
                !/
                State_VGB(Ne2Rho_,i,j,k,iBlock) = RhoNeutralsISW*cTiny
                State_VGB(Ne2P_,i,j,k,iBlock) = PNeutralsISW*cTiny
                State_VGB(Ne2RhoUx_,i,j,k,iBlock) = UxNeutralsISW*RhoNeutralsISW*cTiny  

                !
                ! PopIII - o que fazer com set_ICs popIII?
                !/
                State_VGB(Ne3Rho_,i,j,k,iBlock) = RhoNeutralsISW*cTiny
                State_VGB(Ne3P_,i,j,k,iBlock) = PNeutralsISW*cTiny
                State_VGB(Ne3RhoUx_,i,j,k,iBlock) = 0.0

             end if
          end do
       end do
    end do

  end subroutine user_set_ics

  !=====================================================================

  subroutine user_write_progress
    use ModMain
    use ModPhysics
    use ModMultiFluid
    implicit none

    !
    !-------------------------------------------------------------------


    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_rho_dim [n/cc]:',SWH_rho_dim,'SWH_rho:',SWH_rho

    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_Ux_dim  [km/s]:',SWH_Ux_dim,'SWH_Ux:',SWH_Ux
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_Uy_dim  [km/s]:',SWH_Uy_dim,'SWH_Uy:',SWH_Uy 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_Uz_dim  [km/s]:',SWH_Uz_dim,'SWH_Uz:',SWH_Uz
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_p_dim   [ nPa]:',SWH_p_dim,'SWH_p:',SWH_p
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_Bx_dim  [  nT]:',SWH_Bx_dim,'SWH_Bx:',SWH_Bx 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_By_dim  [  nT]:',SWH_By_dim,'SWH_By:',SWH_By
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWH_Bz_dim  [  nT]:',SWH_Bz_dim,'SWH_Bz:',SWH_Bz
    write(*,'(10X,A19,F15.6)')           'SWH_T_dim   [   K]:',SWH_T_dim
    ! fast solar wind
    write(*,*)
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWfast_rho_dim:',SWfast_rho_dim,'SWfast_rho:',SWfast_rho
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWfast_Ux_dim[km/s]:',SWfast_Ux_dim,'SWfast_rho:',SWfast_rho
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWfast_Uy_dim[km/s]:',SWfast_Uy_dim,'SWfast_Uy:',SWfast_Uy 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWfast_Uz_dim[km/s]:',SWfast_Uz_dim,'SWfast_Uz:',SWfast_Uz 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'SWfast_p_dim[nPa]:',SWfast_p_dim,'SWfast_p:',SWfast_p  
    write(*,'(10X,A19,F15.6)')           'SWfast_a_dim[km/s]:',SWfast_a_dim
    write(*,'(10X,A19,F15.6)')           'SWfast_T_dim[K]:',SWfast_T_dim
    write(*,*)
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_rho_dim[n/cc]:',VLISW_rho_dim,'VLISW_rho:',VLISW_rho 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_Ux_dim[km/s]: ',VLISW_Ux_dim,'VLISW_Ux:',VLISW_Ux
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_Uy_dim[km/s]: ',VLISW_Uy_dim,'VLISW_Uy:',VLISW_Uy
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_Uz_dim[km/s]: ',VLISW_Uz_dim,'VLISW_Uz:',VLISW_Uz 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_p_dim [nPa]: ',VLISW_p_dim,'VLISW_p:',VLISW_p 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_Bx_dim[nT]: ',VLISW_Bx_dim,'VLISW_Bx:',VLISW_Bx
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_By_dim[nT]:',VLISW_By_dim,'VLISW_By:',VLISW_By
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'VLISW_Bz_dim[nT]:',VLISW_Bz_dim,'VLISW_Bz:',VLISW_Bz
    write(*,'(10X,A19,F15.6)')           'VLISW_a_dim[km/s]: ',VLISW_a_dim
    write(*,'(10X,A19,F15.6)')           'VLISW_T_dim[K]: ',VLISW_T_dim! 
    !neutrals
    write(*,*)     
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'RhoNeutralsISW_dim:',RhoNeutralsISW_dim ,'RhoNeutralsISW:',RhoNeutralsISW 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'UxNeutralsISW_dim:',UxNeutralsISW_dim,'UxNeutralsISW:',UxNeutralsISW 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'UyNeutralsISW_dim:',UyNeutralsISW_dim,'UyNeutralsISW:',UyNeutralsISW
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'UzNeutralsISW_dim:',UzNeutralsISW_dim,'UzNeutralsISW:',UzNeutralsISW 
    write(*,'(10X,A19,F15.6,A11,F15.6)') 'PNeutralsISW_dim:',PNeutralsISW_dim,'PNeutralsISW:',PNeutralsISW
    write(*,'(10X,A19,F15.6)')           'TNeutralsISW_dim:',TNeutralsISW_dim     
    write(*,*)

    !------------------------------------------------------------------

  end subroutine user_write_progress

  !=====================================================================

  subroutine user_io_units
    use ModPhysics 
    use ModProcMH, ONLY:iProc
    use ModMain
    use ModVarIndexes
    !
    character (len=*), parameter :: Name='user_io_units'

    Io2Si_V(UnitX_)           = 215*Rsun                  ! R  
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

    ! if(globalBLK==100) then
    !  write(*,*) 'Io2Si_V(UnitX_)',Io2Si_V(UnitX_)
    ! end if

    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    !  normalization of SWH and VLISW and Neutrals
    VLISW_a_dim    = No2Io_V(UnitU_)*(VLISW_T_dim/SWH_T_dim)
    VLISW_p_dim    = No2Io_V(UnitP_)*inv_g*(VLISW_rho_dim/SWH_rho_dim)*(VLISW_T_dim/SWH_T_dim)
    VLISW_B_factor = No2Io_V(UnitB_)*sqrt((VLISW_T_dim/SWH_T_dim)*(VLISW_rho_dim/SWH_rho_dim))

    VLISW_rho = VLISW_rho_dim*Io2No_V(UnitRho_)
    VLISW_p   = VLISW_p_dim*Io2No_V(UnitP_)
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

    SWH_p_dim = No2Io_V(UnitP_)*inv_g*(SWH_rho_dim/SWH_rho_dim)*(SWH_T_dim/SWH_T_dim)

    SWH_p   = SWH_p_dim*Io2No_V(UnitP_)
    SWH_Ux  = SWH_Ux_dim*Io2No_V(UnitU_)
    SWH_Uy  = SWH_Uy_dim*Io2No_V(UnitU_)
    SWH_Uz  = SWH_Uz_dim*Io2No_V(UnitU_)
    SWH_Bx  = SWH_Bx_dim*Io2No_V(UnitB_)
    SWH_By  = SWH_By_dim*Io2No_V(UnitB_)
    SWH_Bz  = SWH_Bz_dim*Io2No_V(UnitB_)

    SWfast_p_dim = No2Io_V(UnitP_)*inv_g*(SWfast_rho_dim/SWH_rho_dim)
    SWfast_p = SWfast_p_dim*Io2No_V(UnitP_)
    !
    !
    ! The units of rho_dim are n/cc and unitUSER_rho g/cc
    !/
    !merav june01    PNeutralsISW   = RhoNeutralsISW * 
    !TNeutralsISW_dim*Io2No_V(UnitTemperature_)
    !merav june01  SWH_p   = SWH_rho * SW_T_dim*Io2No_V(UnitTemperature_)

    RhoNeutralsISW = RhoNeutralsISW_dim*Io2No_V(UnitRho_)
    PNeutralsISW_dim = No2Io_V(UnitP_)*inv_g*(RhoNeutralsISW_dim/SWH_rho_dim)*(TNeutralsISW_dim /SWH_T_dim)

    PNeutralsISW = PNeutralsISW_dim*Io2No_V(UnitP_)
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

    integer :: i, j, k, iVar
    character (len=*), parameter :: Name='user_set_plot_var'
    !-------------------------------------------------------------------

    UsePlotVarBody = .true.
    PlotVarBody    = 0.0
    
    select case(NameVar)
    case('srho')
       NameTecVar = 'Srho'
       PlotVar_G(1:nI,1:nJ,1:nK) = Source_VC(NeuRho_,:,:,:)
       IsFound = .true.
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !====================================================================

  subroutine user_calc_sources
    ! Calculates the charge exchange cross sections for the neutrals
    ! 
    ! Written by Merav Opher April 21, 2002
    !/
    ! *modified for the version 7.5.2 May, 2002*
    ! *revistied and modified on Nov, 2002*
    ! source terms for the plasma
    ! modified to take of the S=0 Feb 08 (initialization) 
    ! modified June 08, 2007 to take into account multi-fluids
    ! modified August 18;September10, 2007 to take into account 3-fluids
    ! october 27 checking units 
    ! January 01, implementing implicit source terms 
    ! Help of Gabor Toth
    !-------------------------------------------------------------------
    use ModProcMH
    use ModPointImplicit, ONLY:  UsePointImplicit, UsePointImplicit_B, &
         IsPointImplSource, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC

    use ModMain
    use ModVarIndexes, ONLY: rho_, rhoUx_, rhoUy_, rhoUz_,p_,Bx_, By_, Bz_, nVar
    use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, R_BLK,&
         body_BLK, Rmin_BLK, vInv_CB
    use ModAdvance,  ONLY: State_VGB, Source_VC, &
         Rho_, RhoUx_, RhoUy_, RhoUz_, Energy_
    use ModPhysics
    use ModNumConst
    use ModCovariant, ONLY: UseCovariant 

    character (len=*), parameter :: Name='user_calc_sources'
    real ::  Ux, Uy, Uz
    real :: usq, usqdim, Tproton
    real :: cstar,cstarl,cstarln,cTpro,cth, term3, term4
    real :: test,test1,  testterm1, testterm2, testusq, testneutrals
    real :: testterm3, testterm4 
    real :: Entropy, Mask 
    real::UTh2
    real::x0,y0,z0,rp

    real, dimension(3)::I0xp,I0px,Tneutral,urel2,urel2dim,vNeutrals,UThN2
    real, dimension(3)::JpxUx,JpxUy,JpxUz,ureldim,termpx,termxp,termepx
    real, dimension(3)::JxpUx,JxpUy,JxpUz,QmxpUx,QmxpUy,QmxpUz,Kxp,Kpx,Qexp
    real, dimension(3)::ustar,sigma,n_H,nuc,nucNORM,UxN,UyN,UzN,Friction
    real, dimension(3)::QmpxUx,QmpxUy,QmpxUz,Qepx

    integer :: xN  !this is the index to indicate Pop 1, 2, 3 of neutrals

    logical :: OutsideHeliopause = .false.
    logical:: oktest=.false.,oktest_me

    integer :: i, j, k, iDim, iBlock
    real    :: Coef

    !-----------------------------------------------------------------------
    iBlock = GlobalBlk

    ! Only blocks within radius rFriction need to be point implicit
    UsePointImplicit_B(iBlock) = rMin_BLK(iBlock) < rFriction

    ! Tproton has units of kelvin
    ! It's the proton temperature that enters here. Tproton = 0.5*Tplasma
    ! see Liewer & Brackbill '97
    ! unitUSER_p is in cgs dyne/cm2. The factor of 1E3 bring the mp to g and
    ! the factor of 1.E7 bring kB to erg/K
    !
    if(oktest_me)then
       write(*,*)'Source(p,E)', Source_VC(P_:P_+1,iTest,jTest,kTest)
    end if

    !  calculating first some constants
    cTpro = 1.0/(cBoltzmann*1.E5)
    cstar = (128.0*cBoltzmann)/(9.0*cPi*mNeutrals)
    cstarl = 8.0*cBoltzmann/(cPi*mNeutrals)
    cth = 2.0*cBoltzmann/mNeutrals
    cstarln = (9.0*cPi*cPi/16.0)*cstarl
    !cstar eh J/kg
    !cstar*T eh (m/s)^2 
    !No2Si_V m/s
    !/
    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       x0 = x_BLK(i,j,k,iBlock)
       y0 = y_BLK(i,j,k,iBlock)
       z0 = z_BLK(i,j,k,iBlock)

       rp=sqrt(x0**2+y0**2+z0**2)

       Ux = State_VGB(rhoUx_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       Uy = State_VGB(rhoUy_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock)
       Uz = State_VGB(rhoUz_,i,j,k,iBlock)/State_VGB(Rho_,i,j,k,iBlock) 

       UxN(1)=State_VGB(NeuRhoUx_,i,j,k,iBlock)/State_VGB(NeuRho_,i,j,k,iBlock)
       UyN(1)=State_VGB(NeuRhoUy_,i,j,k,iBlock)/State_VGB(NeuRho_,i,j,k,iBlock)
       UzN(1)=State_VGB(NeuRhoUz_,i,j,k,iBlock)/State_VGB(NeuRho_,i,j,k,iBlock)

       UxN(2)=State_VGB(Ne2RhoUx_,i,j,k,iBlock)/State_VGB(Ne2Rho_,i,j,k,iBlock)
       UyN(2)=State_VGB(Ne2RhoUy_,i,j,k,iBlock)/State_VGB(Ne2Rho_,i,j,k,iBlock)
       UzN(2)=State_VGB(Ne2RhoUz_,i,j,k,iBlock)/State_VGB(Ne2Rho_,i,j,k,iBlock)

       UxN(3)=State_VGB(Ne3RhoUx_,i,j,k,iBlock)/State_VGB(Ne3Rho_,i,j,k,iBlock)
       UyN(3)=State_VGB(Ne3RhoUy_,i,j,k,iBlock)/State_VGB(Ne3Rho_,i,j,k,iBlock)
       UzN(3)=State_VGB(Ne3RhoUz_,i,j,k,iBlock)/State_VGB(Ne3Rho_,i,j,k,iBlock)

       usq = sqrt(Ux**2.0 + Uy**2.0 + Uz**2.0) 

       usqdim = usq*(No2Si_V(UnitU_))*1.E-3

       !units of usqdim is km/s

       Tproton = 0.5*cTpro*State_VGB(P_,i,j,k,iBlock)*No2Io_V(UnitP_)/ & 
            (State_VGB(Rho_,i,j,k,iBlock)*No2Io_V(UnitRho_))

       Tneutral(1) = 0.5*cTpro*State_VGB(NeuP_,i,j,k,iBlock)*No2Io_V(UnitP_)/ &
            (State_VGB(NeuRho_,i,j,k,iBlock)*No2Io_V(UnitRho_))

       Tneutral(2) = 0.5*cTpro*State_VGB(Ne2P_,i,j,k,iBlock)*No2Io_V(UnitP_)/ &
            (State_VGB(Ne2Rho_,i,j,k,iBlock)*No2Io_V(UnitRho_))

       Tneutral(3) = 0.5*cTpro*State_VGB(Ne3P_,i,j,k,iBlock)*No2Io_V(UnitP_)/ &
            (State_VGB(Ne3Rho_,i,j,k,iBlock)*No2Io_V(UnitRho_))

       ! Calculating the thermal speed for the three populations of neutrals
       ! units of UthN2 and Uth2 are (m/s)^2

       UThN2(1)=cth*Tneutral(1)
       UThN2(2)=cth*Tneutral(2)
       UThN2(3)=cth*Tneutral(3)
       UTh2=cth*Tproton

       ! Relative velocity between neutrals and plasma 
       ! the sqrt is to ensure that urel >0.

       vNeutrals(1) = sqrt(UxN(1)**2.0 + UyN(1)**2.0 + UzN(1)**2.0)
       !!               urel2(1) = (vNeutrals(1)-usq)**2.0
       !! i am not doing this way above because of regions where urel2 drops to zero or negative

       urel2(1) = (UxN(1) - Ux)**2 + (UyN(1) - Uy)**2 + (UzN(1) - Uz)**2 

       vNeutrals(2) = sqrt(UxN(2)**2.0 + UyN(2)**2.0 + UzN(2)**2.0)

       urel2(2) = (UxN(2) - Ux)**2 + (UyN(2) - Uy)**2 + (UzN(2) - Uz)**2
       !!              urel2(2) = (vNeutrals(2)-usq)**2.0

       vNeutrals(3) = sqrt(UxN(3)**2.0 + UyN(3)**2.0 + UzN(3)**2.0)
       urel2(3) = (UxN(3) - Ux)**2 + (UyN(3) - Uy)**2 + (UzN(3) - Uz)**2
       !!              urel2(3) = (vNeutrals(3)-usq)**2.0

       ! Incorporating units in variables so we can calculate the charge exchange cross 
       ! sections
       ! No2Si_V has units of m/s like cstartT so ureldim and ustar has units of m/s

       urel2dim(1) = urel2(1)*((No2Si_V(UnitU_))**2.0)
       ustar(1) = sqrt(urel2dim(1) + cstar*(Tproton + TNeutral(1)))   
       ureldim(1)=sqrt(urel2dim(1)) 

       ! here i put cstar  to be consistent with i did before but
       ! seem to me Zank et al. 1996 is using cstartl 
       ! also before i was using ureldim (uxN-ux)**2 etc instead of vneutrals-usq like i am 
       ! doing here
       !
       urel2dim(2) = urel2(2)*((No2Si_V(UnitU_))**2.0)
       ustar(2) = sqrt(urel2dim(2) + cstarl*(Tproton + TNeutral(2)))   
       ureldim(2)=sqrt(urel2dim(2))

       urel2dim(3) = urel2(3)*((No2Si_V(UnitU_))**2.0)
       ustar(3) = sqrt(urel2dim(3) + cstarl*(Tproton + TNeutral(3)))   
       ureldim(3)=sqrt(urel2dim(3))

       !Calculating some common terms
       ! units of tempx and termxp are (m/s)^-1

       termpx(1)=1.0/sqrt(4*(cstarl*Tproton+urel2dim(1))+cstarln*Tneutral(1))
       termpx(2)=1.0/sqrt(4*(cstarl*Tproton+urel2dim(2))+cstarln*Tneutral(2))
       termpx(3)=1.0/sqrt(4*(cstarl*Tproton+urel2dim(3))+cstarln*Tneutral(3))     

       termxp(1)=1.0/sqrt(4*(cstarl*Tneutral(1)+urel2dim(1))+cstarln*Tproton)     
       termxp(2)=1.0/sqrt(4*(cstarl*Tneutral(2)+urel2dim(2))+cstarln*Tproton)
       termxp(3)=1.0/sqrt(4*(cstarl*Tneutral(3)+urel2dim(3))+cstarln*Tproton) 

       ! units of termepx is m/s

       termepx(1)=sqrt(cstarl*Tproton+cstar*Tneutral(1)+urel2dim(1))
       termepx(2)=sqrt(cstarl*Tproton+cstar*Tneutral(2)+urel2dim(2))
       termepx(3)=sqrt(cstarl*Tproton+cstar*Tneutral(3)+urel2dim(3))

       ! Maher and Tinsley cross section - units of cm^2 (Timur's paper has a misprint 
       ! where he
       ! had + instead of a -
       ! sigma has ureldim look at Baranov et al. 1991              
       ! ureldim has to have units of cm/s so the factor 100 is to pass m to cm
       ! probken when ureldim=0             sigma = (1.64E-7 - (6.95E-9)*log(ureldim*100))**2.0 
       ! Linde uses ustar, page 53

       sigma(1) = (1.64E-7 - (6.95E-9)*log(ustar(1)*100.))**2.0
       sigma(2) = (1.64E-7 - (6.95E-9)*log(ustar(2)*100.))**2.0
       sigma(3) = (1.64E-7 - (6.95E-9)*log(ustar(3)*100.))**2.0

       ! Get the neutral density in units
       ! The charge exhange cross section 
       ! 100 to change ustar to cm/s

       nuc(1) = sigma(1)*State_VGB(NeuRho_,i,j,k,iBlock)*No2Io_V(UnitRho_)*ustar(1)*100.
       nuc(2) = sigma(2)*State_VGB(Ne2Rho_,i,j,k,iBlock)*No2Io_V(UnitRho_)*ustar(2)*100.
       nuc(3) = sigma(3)*State_VGB(Ne3Rho_,i,j,k,iBlock)*No2Io_V(UnitRho_)*ustar(3)*100.
       ! Now everything normalized

       nucNORM(1)=nuc(1)*No2Io_V(UnitT_)
       nucNORM(2)=nuc(2)*No2Io_V(UnitT_)
       nucNORM(3)=nuc(3)*No2Io_V(UnitT_)

       !
       ! source terms for neutrals
       !/
       ! we are ignoring the creation of very energetic neutrals, which are created by charge exchange
       ! see Liewer et al. 1996
       ! we probably want to put a flag based on entropy do identify the heliopause   

       Mask = 0.0

       if (rp > 40.0) then 
          Mask = 0.5*(1.0-sign(1.0,usqdim-24.8))+(1.0-sign(1.0,usqdim-400.0))+ &
               0.5*(1.0-sign(1.0,Tproton - 1.15E4))
       end if

       Friction(1) = nucNORM(1)*State_VGB(Rho_,i,j,k,iBlock)
       Friction(2) = nucNORM(2)*State_VGB(Rho_,i,j,k,iBlock)
       Friction(3) = nucNORM(3)*State_VGB(Rho_,i,j,k,iBlock)

       !
       ! Declaration of all the I',J,K !following my notes of Sep 10, 2007
       !/

       I0xp(1)=Friction(1)
       I0xp(2)=Friction(2)
       I0xp(3)=Friction(3)
       I0px(1)=Friction(1)
       I0px(2)=Friction(2)
       I0px(3)=Friction(3)

       ! units are fine: (Uth2/ustar)*termxp is unitless as it should be

       JxpUx(1)=Friction(1)*(Ux+((UxN(1)-Ux)*(UTh2/ustar(1))*termxp(1)))
       JxpUy(1)=Friction(1)*(Uy+((UyN(1)-Uy)*(UTh2/ustar(1))*termxp(1)))
       JxpUz(1)=Friction(1)*(Uz+((UzN(1)-Uz)*(UTh2/ustar(1))*termxp(1)))

       JxpUx(2)=Friction(2)*(Ux+((UxN(2)-Ux)*(UTh2/ustar(2))*termxp(2)))
       JxpUy(2)=Friction(2)*(Uy+((UyN(2)-Uy)*(UTh2/ustar(2))*termxp(2)))
       JxpUz(2)=Friction(2)*(Uz+((UzN(2)-Uz)*(UTh2/ustar(2))*termxp(2)))

       JxpUx(3)=Friction(3)*(Ux+((UxN(3)-Ux)*(UTh2/ustar(3))*termxp(3)))
       JxpUy(3)=Friction(3)*(Uy+((UyN(3)-Uy)*(UTh2/ustar(3))*termxp(3)))
       JxpUz(3)=Friction(3)*(Uz+((UzN(3)-Uz)*(UTh2/ustar(3))*termxp(3)))

       JpxUx(1)=Friction(1)*(UxN(1)+((Ux-UxN(1))*(UThN2(1)/ustar(1))*termpx(1)))
       JpxUy(1)=Friction(1)*(UyN(1)+((Uy-UyN(1))*(UThN2(1)/ustar(1))*termpx(1)))
       JpxUz(1)=Friction(1)*(UzN(1)+((Uz-UzN(1))*(UThN2(1)/ustar(1))*termpx(1)))

       JpxUx(2)=Friction(2)*(UxN(2)+((Ux-UxN(2))*(UThN2(2)/ustar(2))*termpx(2)))
       JpxUy(2)=Friction(2)*(UyN(2)+((Uy-UyN(2))*(UThN2(2)/ustar(2))*termpx(2)))
       JpxUz(2)=Friction(2)*(UzN(2)+((Uz-UzN(2))*(UThN2(2)/ustar(2))*termpx(2)))


       JpxUx(3)=Friction(3)*(UxN(3)+((Ux-UxN(3))*(UThN2(3)/ustar(3))*termpx(3)))
       JpxUy(3)=Friction(3)*(UyN(3)+((Uy-UyN(3))*(UThN2(3)/ustar(3))*termpx(3)))
       JpxUz(3)=Friction(3)*(UzN(3)+((Uz-UzN(3))*(UThN2(3)/ustar(3))*termpx(3)))

       QmxpUx(1)=JxpUx(1)-JpxUx(1)
       QmxpUy(1)=JxpUy(1)-JpxUy(1)
       QmxpUz(1)=JxpUz(1)-JpxUz(1)

       QmxpUx(2)=JxpUx(2)-JpxUx(2)
       QmxpUy(2)=JxpUy(2)-JpxUy(2)
       QmxpUz(2)=JxpUz(2)-JpxUz(2)

       QmxpUx(3)=JxpUx(3)-JpxUx(3)
       QmxpUy(3)=JxpUy(3)-JpxUy(3)
       QmxpUz(3)=JxpUz(3)-JpxUz(3)

       Kxp(1)=0.5*Friction(1)*(usq**2+(1.5*(UTh2/ustar(1))*termepx(1)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*usq*(vNeutrals(1)-usq)*(UThN2(1)/ustar(1))*termpx(1))

       Kxp(2)=0.5*Friction(2)*(usq**2+(1.5*(UTh2/ustar(2))*termepx(2)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*usq*(vNeutrals(2)-usq)*(UThN2(2)/ustar(2))*termpx(2))

       Kxp(3)=0.5*Friction(3)*(usq**2+(1.5*(UTh2/ustar(3))*termepx(3)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*usq*(vNeutrals(3)-usq)*(UThN2(3)/ustar(3))*termpx(3))

       Kpx(1)=0.5*Friction(1)*(vNeutrals(1)**2+(1.5*(UThN2(1)/ustar(1))*termepx(1)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*vNeutrals(1)*(usq-vNeutrals(1))*(UThN2(1)/ustar(1))*termpx(1))

       Kpx(2)=0.5*Friction(2)*(vNeutrals(2)**2+(1.5*(UThN2(2)/ustar(2))*termepx(2)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*vNeutrals(2)*(usq-vNeutrals(2))*(UThN2(2)/ustar(2))*termpx(2))

       Kpx(3)=0.5*Friction(3)*(vNeutrals(3)**2+(1.5*(UThN2(3)/ustar(3))*termepx(3)*(1/((No2Si_V(UnitU_))**2.0)))+ &
            2*vNeutrals(3)*(usq-vNeutrals(3))*(UThN2(3)/ustar(3))*termpx(3))

       Qexp(1)=Kxp(1)-Kpx(1)
       Qexp(2)=Kxp(2)-Kpx(2)
       Qexp(3)=Kxp(3)-Kpx(3)

       ! Check if point implicit scheme is on (part implicit may switch it off)
       ! Also check if this particular block is point implicit or not
       if(.not.(UsePointImplicit .and. UsePointImplicit_B(iBlock)))then
          ! Add all source terms if we do not use the point implicit scheme
          call user_expl_source
          call user_impl_source
       elseif(IsPointImplSource)then
          ! Add implicit sources only
          call user_impl_source
       else
          ! Add explicit sources only
          call user_expl_source
       end if
    end do; end do; end do

  contains

    !==========================================================================
    subroutine user_expl_source
      ! Add explicit source here

      ! Here come the explicit source terms
      !
      ! Pop I
      !/
      if (Mask >= 3.0) then
         Source_VC(NeuRho_,i,j,k) = Source_VC(NeuRho_,i,j,k) +I0xp(2) + I0xp(3)
         Source_VC(NeuEnergy_,i,j,k) = Source_VC(NeuEnergy_,i,j,k) +        &
              Qexp(1) + Kxp(2) + Kxp(3)
         ! writing the source P as SP = [SE -uxS(rhoux)-uyS(rhoUy) - uzS(rhoUz)]*gamma-1

      else
         Source_VC(NeuRho_,i,j,k) = Source_VC(NeuRho_,i,j,k) - I0px(1)
         Source_VC(NeuEnergy_,i,j,k) = Source_VC(NeuEnergy_,i,j,k) &
              - Kpx(1)
      end if
      !
      ! Pop II
      !/
      if (Mask < 3.0) then
         if (Mask > 1.0) then
            Source_VC(Ne2Rho_,i,j,k) = Source_VC(Ne2Rho_,i,j,k) + &
                 I0xp(1) + I0xp(3)
            Source_VC(Ne2Energy_,i,j,k) = Source_VC(Ne2Energy_,i,j,k) &
                 + Qexp(2) + Kxp(1) + Kxp(3)
         else
            Source_VC(Ne2Rho_,i,j,k) = Source_VC(Ne2Rho_,i,j,k) -I0px(2) 
            Source_VC(Ne2Energy_,i,j,k) = Source_VC(Ne2Energy_,i,j,k) &
                 - Kpx(2)
         end if
      end if
      !
      ! Pop III
      !/
      if (Mask <= 1.0) then
         Source_VC(Ne3Rho_,i,j,k) = Source_VC(Ne3Rho_,i,j,k) & 
              + I0xp(1) + I0xp(2)
         Source_VC(Ne3Energy_,i,j,k) = Source_VC(Ne3Energy_,i,j,k) &
              + Qexp(3) + Kxp(1) + Kxp(2)
      else
         Source_VC(Ne3Rho_,i,j,k) = Source_VC(Ne3Rho_,i,j,k) -I0px(3) 
         Source_VC(Ne3Energy_,i,j,k) = Source_VC(Ne3Energy_,i,j,k) &
              -Kpx(3)
      end if
      !
      ! Source terms for Plasma
      !/
      Source_VC(Rho_,i,j,k)   = Source_VC(Rho_,i,j,k)  + 0.0

      !correcting the source terms (see Gombosi 2000)

      Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k) +        &
           Qepx(1)+ Qepx(2) + Qepx(3)

    end subroutine user_expl_source
    !==========================================================================
    subroutine user_impl_source

      !  This is where the point implicit source terms go
      if(r_BLK(i,j,k,iBlock) > rFriction ) RETURN


      !   All the neutrals momenta are implicit

      ! Source terms for Plasma

      Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k) +    &
           (g-1)*(Qepx(1)+ Qepx(2) + Qepx(3) &
           - Ux*Source_VC(RhoUx_,i,j,k) &
           - Uy*Source_VC(RhoUy_,i,j,k) &
           - Uz*Source_VC(RhoUz_,i,j,k))


      !
      ! Pop I
      !/
      if (Mask >= 3.0) then

         Source_VC(NeuRhoUx_,i,j,k) = Source_VC(NeuRhoUx_,i,j,k) +  &
              QmxpUx(1) + JxpUx(2) + JxpUx(3)
         Source_VC(NeuRhoUy_,i,j,k) = Source_VC(NeuRhoUy_,i,j,k) +  &
              QmxpUy(1) + JxpUy(2) + JxpUy(3)
         Source_VC(NeuRhoUz_,i,j,k) = Source_VC(NeuRhoUz_,i,j,k) +  & 
              QmxpUz(1) + JxpUz(2) + JxpUz(3)

         Source_VC(NeuP_,i,j,k) = Source_VC(NeuP_,i,j,k) +    &
              (g-1)*(Qexp(1) + Kxp(2) + Kxp(3) &
              - UxN(1)*Source_VC(NeuRhoUx_,i,j,k) &
              - UyN(1)*Source_VC(NeuRhoUy_,i,j,k) &
              - UzN(1)*Source_VC(NeuRhoUz_,i,j,k)) 
      else
         Source_VC(NeuRhoUx_,i,j,k) = Source_VC(NeuRhoUx_,i,j,k) &
              -JpxUx(1)
         Source_VC(NeuRhoUy_,i,j,k) = Source_VC(NeuRhoUy_,i,j,k) &
              -JpxUy(1)
         Source_VC(NeuRhoUz_,i,j,k) = Source_VC(NeuRhoUz_,i,j,k) &
              -JpxUz(1)
         Source_VC(NeuP_,i,j,k) = Source_VC(NeuP_,i,j,k) +    &
              (g-1)*( -kpx(1) &
              - UxN(1)*Source_VC(NeuRhoUx_,i,j,k) &
              - UyN(1)*Source_VC(NeuRhoUy_,i,j,k) &
              - UzN(1)*Source_VC(NeuRhoUz_,i,j,k))
      end if
      !
      ! Pop II
      !/
      if (Mask < 3.0) then
         if (Mask > 1.0) then
            Source_VC(Ne2RhoUx_,i,j,k) = Source_VC(Ne2RhoUx_,i,j,k) +  &
                 QmxpUx(2) + JxpUx(1) + JxpUx(3)
            Source_VC(Ne2RhoUy_,i,j,k) = Source_VC(Ne2RhoUy_,i,j,k) +  &
                 QmxpUy(2) + JxpUy(1) + JxpUy(3)
            Source_VC(Ne2RhoUz_,i,j,k) = Source_VC(Ne2RhoUz_,i,j,k) +  & 
                 QmxpUz(2) + JxpUz(1) + JxpUz(3)
            Source_VC(Ne2P_,i,j,k) = Source_VC(Ne2P_,i,j,k) +    &
                 (g-1)*(Qexp(2) + Kxp(1) + Kxp(3) &
                 - UxN(2)*Source_VC(Ne2RhoUx_,i,j,k) &
                 - UyN(2)*Source_VC(Ne2RhoUy_,i,j,k) &
                 - UzN(2)*Source_VC(Ne2RhoUz_,i,j,k))
         else
            Source_VC(Ne2RhoUx_,i,j,k) = Source_VC(Ne2RhoUx_,i,j,k) &
                 -JpxUx(2)
            Source_VC(Ne2RhoUy_,i,j,k) = Source_VC(Ne2RhoUy_,i,j,k) &
                 -JpxUy(2)
            Source_VC(Ne2RhoUz_,i,j,k) = Source_VC(Ne2RhoUz_,i,j,k) &
                 -JpxUz(2)
            Source_VC(Ne2P_,i,j,k) = Source_VC(Ne2P_,i,j,k) +    &
                 (g-1)*(- Kpx(2) &
                 - UxN(2)*Source_VC(Ne2RhoUx_,i,j,k) &
                 - UyN(2)*Source_VC(Ne2RhoUy_,i,j,k) &
                 - UzN(2)*Source_VC(Ne2RhoUz_,i,j,k))

         end if
      end if
      !
      ! Pop III
      !/
      if (Mask <= 1.0) then
         Source_VC(Ne3RhoUx_,i,j,k) = Source_VC(Ne3RhoUx_,i,j,k) +  &
              QmxpUx(3) + JxpUx(2) + JxpUx(1)
         Source_VC(Ne3RhoUy_,i,j,k) = Source_VC(Ne3RhoUy_,i,j,k) +  &
              QmxpUy(3) + JxpUy(2) + JxpUy(1)
         Source_VC(Ne3RhoUz_,i,j,k) = Source_VC(Ne3RhoUz_,i,j,k) +  & 
              QmxpUz(3) + JxpUz(2) + JxpUz(1)
         Source_VC(Ne3P_,i,j,k) = Source_VC(Ne3P_,i,j,k) +    &
              (g-1)*(Qexp(3) + Kxp(1) + Kxp(2) &
              - UxN(3)*Source_VC(Ne3RhoUx_,i,j,k) &
              - UyN(3)*Source_VC(Ne3RhoUy_,i,j,k) &
              - UzN(3)*Source_VC(Ne3RhoUz_,i,j,k))
      else
         Source_VC(Ne3RhoUx_,i,j,k) = Source_VC(Ne3RhoUx_,i,j,k)  &
              -JpxUx(3)
         Source_VC(Ne3RhoUy_,i,j,k) = Source_VC(Ne3RhoUy_,i,j,k)  &
              -JpxUy(3)
         Source_VC(Ne3RhoUz_,i,j,k) = Source_VC(Ne3RhoUz_,i,j,k)  &
              -JpxUz(3)
         Source_VC(Ne3P_,i,j,k) = Source_VC(Ne3P_,i,j,k) +    &
              (g-1)*(-Kpx(3) &
              - UxN(3)*Source_VC(Ne3RhoUx_,i,j,k) &
              - UyN(3)*Source_VC(Ne3RhoUy_,i,j,k) &
              - UzN(3)*Source_VC(Ne3RhoUz_,i,j,k))
      end if
      !
      ! Source terms for Plasma
      !/
      Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k) +        &
           QmpxUx(1) + QmpxUx(2) + QmpxUx(3)
      Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k) +        &
           QmpxUy(1) + QmpxUy(2) + QmpxUy(3)
      Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k) +        & 
           QmpxUz(1) + QmpxUz(2) + QmpxUz(3)

    end subroutine user_impl_source

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModVarIndexes, ONLY: RhoUx_,RhoUy_,RhoUz_,P_,NeuRhoUx_,NeuRhoUy_,&
         NeuRhoUz_,NeuP_,Ne2RhoUx_,Ne2RhoUy_,Ne2RhoUz_,Ne2P_,&
         Ne3RhoUx_,Ne3RhoUy_,Ne3RhoUz_,Ne3P_

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet
    !------------------------------------------------------------------------

    ! Allocate and set iVarPointImpl_I
    ! In this example there are 3 implicit variables


    ! All the neutrals momenta and plasma are implicit 
    ! (3 neutral fluid and 1 ion)

    allocate(iVarPointImpl_I(16))

    iVarPointImpl_I = (/RhoUx_, RhoUy_, RhoUz_, P_, NeuRhoUx_ ,&
         NeuRhoUy_, NeuRhoUz_ , NeuP_, Ne2RhoUx_, Ne2RhoUy_, Ne2RhoUz_ ,&
         Ne2P_, Ne3RhoUx_, Ne3RhoUy_, Ne3RhoUz_, Ne3P_/)


    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.

    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

end module ModUser
