! ModUser Outer Helio written by Merav Opher
! Modified for coupling with IH (Museum Project): boundary at 10AU and ISM flowing from +x direction.
! for this project its a single fluid so i deleted the calc_sources and the multi-fluid parts
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
       IMPLEMENTED6 => user_amr_criteria,               &
       IMPLEMENTED7 => user_write_progress,            & 
       IMPLEMENTED8 => user_io_units,                  &
       IMPLEMENTED9 => user_set_plot_var,              &
       IMPLEMENTED10 => user_calc_sources              
      

  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 2.0
  ! I am calling Version Module = 2.0 the version with 3-fluids
  character (len=*), parameter :: NameUserModule = 'Global Heliosphere'

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

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !=================  !  SUBROUTINE USER_SET_INNER_BCS  !======================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModSize,     ONLY: nDim,West_,North_,Top_
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModPhysics  
    use ModNumConst, ONLY: cTolerance
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
       Rbody = 10. !define it

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

       ! VphiSolar = OMEGAbody*(6.96E5)*sinTheta/unitUSER_U
       ! Vphi=omega*r*sinTheta - 6.96E5 is the solar radii in km
       ! No2Si_V(UnitU_) is in m/s so its ok
       ! No2Io_V(UnitU_) is in km/s
       ! factor 1000 for test

       VphiSolar = OMEGAbodyH*(6.96E5)*sinTheta/No2Io_V(UnitU_)
       VphiSolar =0. !in 30AU it is negligble Merav June 04

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
 

    end select
    !-------------------------------------------------------------------

  end subroutine user_face_bcs

  !-------------------------------------------------------------------
  subroutine user_normalization

    use ModProcMH, ONLY:iProc  
    use ModMain 
    use ModPhysics 
    use ModVarIndexes       
    character (len=*), parameter :: Name='user_normalization'
    !-------------------------------------------------------------------

    SWH_rho_dim = 7.02E-2 !n/cm^3
    SWH_T_dim=6.962E3   !K
    No2Si_V(UnitX_)= 215.0*Rsun                                ! m
    No2Si_V(UnitU_)= sqrt(g*cBoltzmann*SWH_T_dim/cProtonMass)
    No2Si_V(UnitRho_)=cProtonMass*SWH_rho_dim*1.0E+6           ! kg/m^3

  end subroutine user_normalization
  !-------------------------------------------------------------------

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
    Rbody = 10. !define it

    SWH_rho_dim = 7.02E-2 !n/cm^3      
    SWH_T_dim=6.962E3   !K      
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

                

             end if
          end do
       end do
    end do

  end subroutine user_set_ics

  !=====================================================================

  subroutine user_amr_criteria(iBLK, UserCriteria, TypeCriteria, IsFound)

    use ModMain
    use ModAdvance
    use ModGeometry, ONLY:x_BLK,y_BLK,z_BLK,R_BLK,&
         dx_BLK,dy_BLK,dz_BLK,true_cell
    use ModPhysics
    use ModConst
    ! 
    ! Variables required by this user subroutine::
    !/
    integer, intent(in):: iBLK
    logical:: IsInRange
    logical, intent(out):: IsFound
    integer:: i,j,k

    real, intent(out):: userCriteria
    !
    ! Local variables::
    !/
    real:: dsMin,dsMax,dsTwo
    real:: XCell,YCell,ZCell,RCell,RCenter
    real:: B0xCell,B0yCell,B0zCell,MinBr,MaxBr
    real:: BIxCell,BIyCell,BIzCell
    real, dimension(1-gcn:nI+gcn,1-gcn:nJ+gcn,1-gcn:nK+gcn):: Br_D
    logical,dimension(3)::IsGhostCell_D

    character (len=20),intent(in):: TypeCriteria
    !-------------------------------------------------------------------

    !
    ! Find the radial location of the center of the block and

    ! the min/max cell size::
    !/
    !    RCenter = cEighth*&

  end subroutine user_amr_criteria

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
    
    case default
       IsFound = .false.
    end select

  end subroutine user_set_plot_var

  !====================================================================

  subroutine user_calc_sources
    !\
    ! Calculates the charge exchange cross sections for the neutrals
    !
    ! 
    !-------------------------------------------------------------------
    use ModProcMH
    use ModPointImplicit, ONLY:  UsePointImplicit, UsePointImplicit_B, &
         IsPointImplSource, iVarPointImpl_I, IsPointImplMatrixSet, DsDu_VVC
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModGeometry, ONLY : dx_BLK, dy_BLK, dz_BLK, R_BLK,&
         body_BLK, Rmin_BLK, vInv_CB
    use ModPhysics
    use ModNumConst

    character (len=*), parameter :: Name='user_calc_sources'
    
    real:: x0,y0,z0,rp
    real :: cTpro, cth, usqdim, Tproton, Mask
    real :: State_V(nVar)  

   
    logical:: oktest=.false.,oktest_me

    ! to help testing the calc sources for the different populations

    

    integer :: i, j, k, jFluid, iDim, iBlock

    !-----------------------------------------------------------------------
   
end subroutine user_calc_sources
end module ModUser
