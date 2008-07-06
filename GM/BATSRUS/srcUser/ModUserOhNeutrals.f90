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
!July 03
!==============================================================================
module ModUser

  use ModSize,     ONLY: nI, nJ, nK, gcn, nBLK
  use ModMain
  use ModPhysics
  use ModSetOuterBC
  use ModAdvance,  ONLY : State_VGB
  use ModGeometry, ONLY : x_BLK, y_BLK, z_BLK
  use ModVarIndexes
  use ModProcMH
  use ModMultiFluid
  use ModUserEmpty,                                     &
       IMPLEMENTED1  => user_read_inputs,               &
       IMPLEMENTED2  => user_face_bcs,                  &
       IMPLEMENTED3  => user_normalization,             &
       IMPLEMENTED4  => user_set_outerbcs,              &
       IMPLEMENTED5  => user_set_ics,                   &
       IMPLEMENTED8  => user_write_progress,            & 
       IMPLEMENTED9  => user_io_units,                  &
       IMPLEMENTED10 => user_set_plot_var,              &
       IMPLEMENTED11 => user_calc_sources,              &
       IMPLEMENTED12 => user_init_point_implicit


  include 'user_module.h' !list of public methods

  real,              parameter :: VersionUserModule = 2.1
  ! I am calling Version Module = 2.0 the version with 3-fluids
  character (len=*), parameter :: &
       NameUserModule = 'Outer Heliosphere with 3 neutrals, Opher & Toth'

  logical :: NeutralsI   = .false.
  logical :: NeutralsII  = .false.
  logical :: NeutralsIII = .false. 

  real :: OmegaSun   = 0.0  ! normalized angular speed of Sun
  real :: ParkerTilt = 0.0  ! Bphi/Br at the equator at r=rBody

  ! SWH variables.
  !/
  real ::      SWH_T_dim=0.0  , &
       SWH_a_dim=0.0  , &
       SWH_rho=0.0, SWH_rho_dim=0.0, &
       SWH_p=0.0  , SWH_p_dim=0.0   , &
       SWH_Ux=0.0 , SWH_Ux_dim=0.0 , &
       SWH_Uy=0.0 , SWH_Uy_dim=0.0 , &
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

  !real, dimension(0:1) :: &
  !     VLISW_rho_t,  &
  !     VLISW_p_t  ,  &
  !     VLISW_Ux_t ,  &
  !     VLISW_Uy_t ,  &
  !     VLISW_Uz_t ,  &
  !     VLISW_Bx_t ,  &
  !     VLISW_By_t ,  &
  !     VLISW_Bz_t   &
  !     VLISW_time_t
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

  !real, dimension(0:1) :: &
  !     SWfast_rho_t,  &
  !     SWfast_p_t  ,  &
  !     SWfast_Ux_t ,  &
  !     SWfast_Uy_t ,  &
  !     SWfast_Uz_t ,  &
  !     SWfast_Bx_t,   &
  !     SWfast_By_t ,  &
  !     SWfast_Bz_t ,  &
  !     SWfast_time_t
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
          ! This is a flag to define how many population of Neutrals to run
       case("#NEUPOPULATIONS")
          call read_var('NeutralsI' ,NeutralsI)
          call read_var('NeutralsII' ,NeutralsII)
          call read_var('NeutralsIII' ,NeutralsIII)

       case default
          if(iProc==0) call stop_mpi( &
               'read_inputs: unrecognized command: '//NameCommand)
       end select
    end do

  end subroutine user_read_inputs

  !==========================================================================
  subroutine user_face_bcs(VarsGhostFace_V)

    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModPhysics  
    use ModSetOuterBC
    use ModProcMH
    use ModFaceBc, ONLY: iBoundary, FaceCoords_D, VarsTrueFace_V, &
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
    character(len=*), parameter:: NameSub='user_face_bcs'
    !-------------------------------------------------------------------

    if(iBoundary /= body1_) &
         call stop_mpi(NameSub//' only inner BC is implemented!')

    if(iProc == ProcTest .and. iBlockBc == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    XyzSph_DD = rot_xyz_sph(FaceCoords_D)
    
    xFace = FaceCoords_D(1)
    yFace = FaceCoords_D(2)
    zFace = FaceCoords_D(3)

    SinTheta = sqrt((xFace**2 + YFace**2)/(xFace**2 + yFace**2 + zFace**2))

    ! calculating the spherical Parker field components
    ! SWH_Bx is the value of the field at the pole B0

    ! Note: use -zFace to invert polarity
    Bsph_D(1) =  sign(SWH_Bx, zFace)             ! Br
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
    VarsGhostFace_V(Rho_) = SWH_rho
    VarsGhostFace_V(p_)   = pSolarWind
    VarsGhostFace_V(Ux_:Uz_) = matmul(XyzSph_DD, Vsph_D)
    VarsGhostFace_V(Bx_:Bz_) = matmul(XyzSph_DD, Bsph_D)

    ! NeuRho is PopI; NeuIIRho is PopII and NeuIIIRho is PopIII
    !
    ! Pop I is going through the inner BCs

    VarsGhostFace_V(NeuRho_) = RhoNeutralsISW !*0.9
    VarsGhostFace_V(NeuUx_)  = UxNeutralsISW
    VarsGhostFace_V(NeuUy_)  = UyNeutralsISW
    VarsGhostFace_V(NeuUz_)  = UzNeutralsISW
    VarsGhostFace_V(NeuP_)   = PNeutralsISW   !*0.9

    ! PopII leaves the domain at a supersonic velocity 
    ! (50km/s while for their temperature 1.E5K their C_s=30km/s)
    ! For the transient case when it flows inward, use 0.001*ion density

    if( sum(VarsTrueFace_V(Ne2Ux_:Ne2Uz_)*FaceCoords_D) > 0.0)then
       VarsGhostFace_V(Ne2Rho_)       = 1.E-3*RhoNeutralsISW
       VarsGhostFace_V(Ne2Ux_:Ne2Uz_) = VarsGhostFace_V(Ux_:Uz_)
       VarsGhostFace_V(Ne2P_)         = 1.E-3*VarsGhostFace_V(p_)
    else
       VarsGhostFace_V(Ne2Rho_:Ne2P_) = VarsTrueFace_V(Ne2Rho_:Ne2P_)
    end if

    ! Pop III has the velocity and temperature of the ions at inner boundary
    ! the density is taken to be 10% of the ions

    VarsGhostFace_V(Ne3Rho_)       = 0.1*VarsGhostFace_V(rho_)
    VarsGhostFace_V(Ne3P_)         = 0.1*VarsGhostFace_V(p_)
    VarsGhostFace_V(Ne3Ux_:Ne3Uz_) = VarsGhostFace_V(Ux_:Uz_)

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
       write(*,*) NameSub,'Pop1=',VarsGhostFace_V(NeuRho_:NeuP_)
       write(*,*) NameSub,'Pop2=',VarsGhostFace_V(Ne2Rho_:Ne2P_)
       write(*,*) NameSub,'Pop3=',VarsGhostFace_V(Ne3Rho_:Ne3P_)
    end if

  end subroutine user_face_bcs

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

  subroutine user_set_outerbcs(iBlock, iSide, TypeBc, IsFound)

    ! The ISM enters at the east boundary (negative x)

    use ModMain
    use ModVarIndexes
    use ModSetOuterBC
    use ModProcMH
    use ModAdvance, ONLY : State_VGB
    use ModMultiFluid

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
    !
    ! PopI
    !
    State_VGB(NeuRho_,-1:2,:,:,iBlock)   = RhoNeutralsISW    
    State_VGB(NeuRhoUx_,-1:2,:,:,iBlock) = RhoNeutralsISW*UxNeutralsISW    
    State_VGB(NeuRhoUy_,-1:2,:,:,iBlock) = RhoNeutralsISW*UyNeutralsISW   
    State_VGB(NeuRhoUz_,-1:2,:,:,iBlock) = RhoNeutralsISW*UzNeutralsISW    
    State_VGB(NeuP_,-1:2,:,:,iBlock)     = PNeutralsISW
    !
    ! In general you should specify as many values as many incoming 
    ! characteristic waves are present. For a neutral fluid this 
    ! is 0 for supersonic outflow, 1 for subsonic outflow, 
    ! 4 for subsonic inflow and 5 for supersonic inflow. 
    !\
    ! PopII and III supersonic outflow
    !/
    do iVar = Ne2Rho_, Ne3P_; do i = -1, 2
       State_VGB(iVar,i,:,:,iBlock) = State_VGB(iVar,3,:,:,iBlock)
    end do; end do

  end subroutine user_set_outerbcs

  !=====================================================================

  subroutine user_set_ics

    use ModMain,  ONLY: globalBLK    
    use ModGeometry, ONLY: r_BLK    
    use ModVarIndexes    
    use ModAdvance,  ONLY: State_VGB    
    use ModPhysics,  ONLY: rBody
    use ModCoordTransform, ONLY: rot_xyz_sph

    implicit none    
    integer :: iBlock    

    integer :: i,j,k

    real :: x, y, z, r, rho0
    real :: b_D(3), v_D(3), bSph_D(3), vSph_D(3)
    real :: SinTheta, SignZ

    real :: XyzSph_DD(3,3) ! rotation matrix Xyz_D = matmul(XyzSph_DD, Sph_D)

    ! These variables are not used now
    ! real :: thetaN, sinthetaN, lambda, RhoSolarW
    ! real :: sin2Theta_fast_wind

    character(len=*), parameter:: NameSub = 'user_set_ics'
    logical :: DoTest, DoTestMe, DoTestCell
    !--------------------------------------------------------------------------
    iBlock = globalBLK

    if(iProc == ProcTest .and. iBlock == BlkTest)then
       call set_oktest(NameSub, DoTest, DoTestMe)
    else
       DoTest = .false.; DoTestMe = .false.
    end if

    if(OmegaSun == 0.0)then
       ! Calculate angular velocity in normalized units
       ! Note: the rotation period is 25.38 days in ModConst.f90:
       ! OmegaSun = cTwoPi/(RotationPeriodSun*Si2No_V(UnitT_))

       OmegaSun   = cTwoPi/(26.0*24.0*3600.00*Si2No_V(UnitT_))
       ParkerTilt = OmegaSun*rBody/SWH_Ux

       if(DoTestMe)then
          write(*,*)NameSub,' OmegaSun         =',OmegaSun
          write(*,*)NameSub,' ParkerTilt       =',ParkerTilt
       end if
    end if

    do i=1-gcn,nI+gcn; do j=1-gcn,nJ+gcn; do k=1-gcn,nK+gcn
       r  = r_BLK(i,j,k,iBlock)
       if(r < rBody) CYCLE

       DoTestCell = DoTestMe .and. i==iTest .and. j==jTest .and. k==kTest

       x = x_BLK(i,j,k,iBlock)
       y = y_BLK(i,j,k,iBlock)
       z = z_BLK(i,j,k,iBlock)

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

       SignZ = sign(1.0, z)
       !SignZ = -sign(1.0,z)

       Bsph_D(1) = SignZ*SWH_Bx*(rBody/r)**2  ! Br
       Bsph_D(2) = 0.0                        ! Btheta
       Bsph_D(3) = -SignZ*SWH_Bx*SinTheta*ParkerTilt*(rBody/r) !Bphi

       ! wrong:   Vphi = OmegaSun*(6.96E5)*sin_theta/No2Io_V(UnitU_)
       ! correct: Vphi = OmegaSun*sin_theta*rSun**2/r
       ! rSun must be in normalized units, of course
       ! Vphi is approximately 0. (30km/s at 1AU, 1km/s at 30AU)

       Vsph_D = (/ SWH_Ux, 0., 0. /)

       !! Introducing the Fast Solar Wind
       !!    sin2Theta_fast_wind = 0.250000    ! 30 degrees
       !!  sin2Theta_fast_wind = 0.1786062   ! 25 degrees
       !!  sin2Theta_fast_wind = 0.116980    ! 20 degrees
       !!  sin2Theta_fast_wind = 1.000000    ! 90 degrees
       !!      if (sin_theta*sin_theta > sin2Theta_fast_wind) then
       !!          !SLOW WIND
       !!          VrS     = SWH_Ux
       !!          RhoSolarW    = SWH_rho
       !!       else
       !!          ! FAST WIND
       !!          VrS     =  SWfast_Ux
       !!          RhoSolarW    =  SWfast_rho
       !!      end if
       !!
       !! I still need to impemenet that the density 
       !! varies between slow and fast solar wind
       !!

       ! magnetic field components in cartesian coordinates
       b_D = matmul(XyzSph_DD, Bsph_D)

       State_VGB(Bx_:Bz_,i,j,k,iBlock) = b_D

       !velocity components in cartesian coordinates
       v_D = matmul(XyzSph_DD, Vsph_D)

       ! density and pressure
       State_VGB(Rho_,i,j,k,iBlock) = SWH_rho * (rBody/r)**2
       State_VGB(P_,i,j,k,iBlock)   = SWH_p   * (rBody/r)**(2*g)

       ! momentum
       State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = State_VGB(rho_,i,j,k,iBlock)*v_D

       !\
       ! PopI
       !/
       State_VGB(NeuRho_,i,j,k,iBlock)  =RhoNeutralsISW
       State_VGB(NeuRhoUx_,i,j,k,iBlock)=RhoNeutralsISW*UxNeutralsISW 
       State_VGB(NeuRhoUy_,i,j,k,iBlock)=RhoNeutralsISW*UyNeutralsISW 
       State_VGB(NeuRhoUz_,i,j,k,iBlock)=RhoNeutralsISW*UzNeutralsISW
       State_VGB(NeuP_,i,j,k,iBlock) = PNeutralsISW

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

    ! if(globalBLK==100) then
    !  write(*,*) 'Io2Si_V(UnitX_)',Io2Si_V(UnitX_)
    ! end if

    Si2Io_V = 1/Io2Si_V
    No2Io_V = No2Si_V*Si2Io_V
    Io2No_V = 1/No2Io_V

    !  normalization of SWH and VLISW and Neutrals
    VLISW_a_dim    = No2Io_V(UnitU_)*(VLISW_T_dim/SWH_T_dim)
    VLISW_p_dim    = No2Io_V(UnitP_)*inv_g &
         *(VLISW_rho_dim/SWH_rho_dim)*(VLISW_T_dim/SWH_T_dim)
    VLISW_B_factor = No2Io_V(UnitB_)*sqrt((VLISW_T_dim/SWH_T_dim) &
         *(VLISW_rho_dim/SWH_rho_dim))

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
    !\
    ! Calculates the charge exchange cross sections for the neutrals
    ! This subroutine calculate the charge exchange between the ionized 
    ! component and the three population of neutrals, Neu, Ne2, Ne3 
    ! following the notation of Pauls et al. (1995) and Zank et al (1996)
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
    ! The _I(1) is the ionized fluid and _I(2)-_I(4) are the neutral fluids
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
    ! 
    !-------------------------------------------------------------------
    use ModProcMH
    use ModPointImplicit, ONLY:  UsePointImplicit, UsePointImplicit_B, &
         IsPointImplSource
    use ModMain
    use ModVarIndexes
    use ModAdvance
    use ModPhysics
    use ModNumConst

    character (len=*), parameter :: Name='user_calc_sources'
    
    real :: cth, usqdim, Tproton, Mask
    real :: State_V(nVar)  

    real, dimension(nFluid) :: &
         Ux_I, Uy_I, Uz_I, UTotal_I, UTotalDim_I, Temp_I, &
         UThS_I, URelS_I, URelSdim_I, UStar_I, Sigma_I, Rate_I, &
         Termpx_I, Termxp_I, Termexp_I, Termepx_I, &
         I0xp_I, I0px_I, I1xp_I, I1px_I, I2xp_I, I2px_I, &
         JxpUx_I, JxpUy_I, JxpUz_I, JpxUx_I, JpxUy_I, JpxUz_I, &
         QmxpUx_I, QmxpUy_I, QmxpUz_I, Kxp_I, Kpx_I, Qexp_I, Qepx_I, &
         QmpxUx_I, QmpxUy_I, QmpxUz_I

    logical:: oktest=.false., oktest_me

    integer :: i, j, k, iBlock, iRegion

    !-----------------------------------------------------------------------
    ! Do not provide explicit source term when point-implicit scheme is used
    ! IsPointImplSource is true only when called from ModPointImplicit
    if(UsePointImplicit .and. .not. IsPointImplSource) RETURN

    iBlock = GlobalBlk

    !  calculating some constants cBoltzmann is J/K 
    
    cth = 2.0*cBoltzmann/mNeutrals

    do k = 1, nK; do j = 1, nJ; do i = 1, nI

       ! Extract conservative variables
       State_V = State_VGB(:,i,j,k,globalBLK)

       Ux_I  = State_V(iRhoUx_I)/State_V(iRho_I)
       Uy_I  = State_V(iRhoUy_I)/State_V(iRho_I)
       Uz_I  = State_V(iRhoUz_I)/State_V(iRho_I)

       ! Total velocity of the ionized and three population of neutrals; 
       UTotal_I  = sqrt(Ux_I**2 + Uy_I**2 + Uz_I**2)

       ! units of usqdim is km/s
       UTotalDim_I  = UTotal_I*No2Si_V(UnitU_)*1.E-3    

       ! Temperature for the ionized and three population of neutrals (K)
       ! Temp_I(1) is the temperature of the ionized component
       ! It's the proton temperature that enters here. Tproton = 0.5*Tplasma
       ! see Liewer & Brackbill '97 
       ! T = p/rho 

       Temp_I    = (State_V(iP_I)/State_V(iRho_I))*No2Si_V(UnitTemperature_)
       Temp_I(1) = 0.5*Temp_I(1)

       ! Thermal speed (squared) for ionized and three populations of neutrals
       ! UThS units are (m/s)^2
       UThS_I = cth*Temp_I 

       ! Relative velocity between neutrals and ionized fluid squared
       URelS_I = (Ux_I - Ux_I(1))**2 &
            +    (Uy_I - Uy_I(1))**2 &
            +    (Uz_I - Uz_I(1))**2 

       ! Calculating Cross Section Sigma_I for the different neutrals
       !
       ! Incorporating units to calculate the charge exchange cross sections
       ! No2Si_V(UnitU_) has units of m/s like cstartT so UReldim and UStar 
       ! has units of m/s

       URelSdim_I  = URelS_I * No2Si_V(UnitU_)**2

       ! UStar_I was slightly different in Liewer et al.. 
       ! I am using Zank et al.'96 UStar has units of m/s

       UStar_I  = sqrt(URelSdim_I + (4./cPi)*(UThS_I +UThS_I(1)))


       ! Maher and Tinsley cross section Sigma 
       ! Linde et al. 1998 has a misprint where there is a + instead of a -
       ! Baranov et al. 1991 uses URelSdim instead of UStar; 
       ! Linde's, Thesis uses UStar
       ! Problem when URelSdim=0         
       ! UStar has to have units of cm/s so the factor 100 is to pass m to cm
       ! Sigma has units of units of m^2
       
       Sigma_I =((1.64E-7 - (6.95E-9)*log(UStar_I*100.))**2)*(1.E-4)

       ! Calculating Rate  = \nu * nH * mp where nH is the density of neutrals
       ! \nu = Sigma*np*u_star where np is the density of the ionized flow and 
       ! For each population of neutrals there will be another Rate
       ! The charge exhange cross section 100 to change ustar to cm/s
       ! Rate has no units (m^2*m/s*s*m-3 )

       Rate_I =Sigma_I*State_V(Rho_)*State_V(iRho_I)*UStar_I  &
              *No2Si_V(UnitRho_)*No2Si_V(UnitT_)*(1./cProtonMass)
       Rate_I(1)=0.
       ! Calculating some common terms
       ! Tempx and Termxp have units of (m/s)^-1

       Termpx_I=1.0/sqrt(4.*((4./cPi)*UThS_I(1)+URelSdim_I)+(9.*cPi/4.)*UThS_I)
       Termxp_I=1.0/sqrt(4.*((4./cPi)*UThS_I+URelSdim_I)+(9.*cPi/4.)*UThS_I(1))

       ! Termepx has units of m/s

       Termexp_I =sqrt((4./cPi)*UThS_I +(64./(9.*cPi))*UThS_I(1)+URelSdim_I) 
       Termepx_I =sqrt((4./cPi)*UThS_I(1)+(64./(9.*cPi))*UThS_I +URelSdim_I) 

       ! Mask is used to find in the domain to the regions separating where the different
       ! populations of Neutrals are created. In regions where the 
       ! different populations are not created they are just destroyed by
       ! charge exchange. 
       ! Mask = 1 inside the Termination Shock; = 2 between the 
       ! Termination Shock and the Heliopause; = 3 beyond the Heliopause

       usqdim= UTotalDim_I(1)
       Tproton=Temp_I(1)
       Mask = 0.0

       Mask = 0.5*(1.0-sign(1.0,usqdim-24.8))+(1.0-sign(1.0,usqdim-350.0)) &
               + 0.5*(1.0-sign(1.0,Tproton - 1.15E4))

       if(Mask <= 1.0)then 
            iRegion = 1       
       endif
       if (Mask >= 3.0)then 
            iRegion = 3
       endif
       if (Mask < 3.0 .and. Mask > 1.0) then 
            iRegion = 2
       endif 
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

       I0xp_I =Rate_I 
       I0px_I =Rate_I 

       I1xp_I  = Rate_I*(UThS_I(1)/UStar_I)*Termxp_I 
       I1px_I  = Rate_I*(UThS_I/UStar_I)*Termpx_I  

       I2xp_I  = Rate_I*(UThS_I(1)/UStar_I)*Termexp_I/No2Si_V(UnitU_)**2
       I2px_I  = Rate_I*(UThS_I/UStar_I)*Termepx_I/No2Si_V(UnitU_)**2

       ! units are fine: (Uth2/ustar)*termxp is unitless as it should be

       JxpUx_I  = Ux_I(1)*Rate_I  + (Ux_I -Ux_I(1))*I1xp_I 
       JxpUy_I  = Uy_I(1)*Rate_I  + (Uy_I -Uy_I(1))*I1xp_I 
       JxpUz_I  = Uz_I(1)*Rate_I  + (Uz_I -Uz_I(1))*I1xp_I 

       JpxUx_I  = Ux_I*Rate_I  + (Ux_I(1)-Ux_I)*I1px_I 
       JpxUy_I  = Uy_I*Rate_I  + (Uy_I(1)-Uy_I)*I1px_I 
       JpxUz_I  = Uz_I*Rate_I  + (Uz_I(1)-Uz_I)*I1px_I 

       QmxpUx_I =JxpUx_I -JpxUx_I 
       QmxpUy_I =JxpUy_I -JpxUy_I 
       QmxpUz_I =JxpUz_I -JpxUz_I 

       QmpxUx_I = -QmxpUx_I 
       QmpxUy_I = -QmxpUy_I 
       QmpxUz_I = -QmxpUz_I    

       Kxp_I = 0.5*(UTotal_I(1)**2)*Rate_I  + (3./4.)*I2xp_I   &
            + (Ux_I(1)*(Ux_I-Ux_I(1)) &
            +  Uy_I(1)*(Uy_I-Uy_I(1)) &
            +  Uz_I(1)*(Uz_I-Uz_I(1)))*I1xp_I
       
       Kpx_I = 0.5*(UTotal_I**2)*Rate_I  + (3./4.)*I2px_I   &
            + (Ux_I*(Ux_I(1)-Ux_I) &
            +  Uy_I*(Uy_I(1)-Uy_I) &
            +  Uz_I*(Uz_I(1)-Uz_I))*I1px_I 

       Qexp_I = Kxp_I -Kpx_I 
       Qepx_I = -Qexp_I  


       ! The energy source is used in the implicit for pressure source, 
       ! and there are no explicit sources other than that. 
       ! So call both subroutines no matter what.
       call user_expl_source 
       call user_impl_source 

    end do; end do; end do

    if(iBlock==BLKtest .and. iProc==PROCtest)then
       call set_oktest('user_calc_sources',oktest,oktest_me)
    else
       oktest=.false.; oktest_me=.false.
    endif

    if(oktest_me)then
       write(*,*) 'Source_VC(Rho_,TestCell)',Source_VC(Rho_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(RhoUx_,TestCell)',Source_VC(RhoUx_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(RhoUy_,TestCell)',Source_VC(RhoUy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(RhoUz_,TestCell)',Source_VC(RhoUz_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(P_,TestCell)',Source_VC(P_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Energy_,TestCell)',Source_VC(Energy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(NeuRho_,TestCell)',Source_VC(NeuRho_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(NeuRhoUx_,TestCell)',Source_VC(NeuRhoUx_,Itest,Jtest,Ktest) 
       write(*,*) 'Source_VC(NeuRhoUy_,TestCell)',Source_VC(NeuRhoUy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(NeuRhoUz_,TestCell)',Source_VC(NeuRhoUz_,Itest,Jtest,Ktest)        
       write(*,*) 'Source_VC(NeuP_,TestCell)',Source_VC(NeuP_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(NeuEnergy_,TestCell)',Source_VC(NeuEnergy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2Rho_,TestCell)',Source_VC(Ne2Rho_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2RhoUx_,TestCell)',Source_VC(Ne2RhoUx_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2RhoUy_,TestCell)',Source_VC(Ne2RhoUy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2RhoUz_,TestCell)',Source_VC(Ne2RhoUz_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2P_,TestCell)',Source_VC(Ne2P_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne2Energy_,TestCell)',Source_VC(Ne2Energy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne3Rho_,TestCell)',Source_VC(Ne3Rho_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne3RhoUx_,TestCell)',Source_VC(Ne3RhoUx_,Itest,Jtest,Ktest) 
       write(*,*) 'Source_VC(Ne3RhoUy_,TestCell)',Source_VC(Ne3RhoUy_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne3RhoUz_,TestCell)',Source_VC(Ne3RhoUz_,Itest,Jtest,Ktest)        
       write(*,*) 'Source_VC(Ne3P_,TestCell)',Source_VC(Ne3P_,Itest,Jtest,Ktest)
       write(*,*) 'Source_VC(Ne3Energy_,TestCell)',Source_VC(Ne3Energy_,Itest,Jtest,Ktest)   
       write(*,*) 'Mask',Mask
       write(*,*) 'Temp_I(1),Temp_I(2),Temp_(3),Temp_(4)',Temp_I(1),Temp_I(2),Temp_I(3),Temp_I(4)
       write(*,*) 'Rate_I(1),Rate_I(2),Rate_I(3),Rate_I(4)',Rate_I(1),Rate_I(2),Rate_I(3),Rate_I(4)
       write(*,*) 'UStar_I(1),UStar_I(2),UStar_I(3),UStar_I(4)',UStar_I(1),UStar_I(2),UStar_I(3),UStar_I(4)
       write(*,*) 'UThS_I',UThS_I
       write(*,*) 'URelSdim_I',URelSdim_I
       write(*,*) 'UTotal_I',UTotal_I
       write(*,*) 'I1xp_I',I1xp_I
       write(*,*) 'I1px_I',I1px_I
       write(*,*) 'I2xp_I',I2xp_I
       write(*,*) 'I2px_I',I2px_I
       write(*,*) 'JxpUx_I',JxpUx_I
       write(*,*) 'JpxUx_I',JpxUx_I
       write(*,*) 'QmxpUx_I',QmxpUx_I
       write(*,*) 'QmpxUx_I',QmpxUx_I
       write(*,*) 'Kxp_I',Kxp_I
       write(*,*) 'Kpx_I',Kpx_I
       write(*,*) 'Qexp_I',Qexp_I
       write(*,*) 'Qepx_I',Qepx_I 
   end if

  contains

    !==========================================================================
    subroutine user_expl_source
      ! Add explicit source here

      ! Here come the explicit source terms
      !
      if(NeutralsI) then      
         if (iRegion == 3) then
            Source_VC(NeuEnergy_,i,j,k) = Source_VC(NeuEnergy_,i,j,k)         &
                 + Qexp_I(2) + Kxp_I(3) + Kxp_I(4)
         else
            Source_VC(NeuEnergy_,i,j,k) = Source_VC(NeuEnergy_,i,j,k) &
                 - Kpx_I(2)
         end if
      end if
      ! Pop II
      !/
      if(NeutralsII) then
         if (iRegion == 2) then
               Source_VC(Ne2Energy_,i,j,k) = Source_VC(Ne2Energy_,i,j,k) &
                    + Qexp_I(3) + Kxp_I(2) + Kxp_I(4)
         else
               Source_VC(Ne2Energy_,i,j,k) = Source_VC(Ne2Energy_,i,j,k) &
                    - Kpx_I(3)
         end if      
      end if
      !
      ! Pop III
      !/
      if(NeutralsIII) then
         if (iRegion == 1) then
            Source_VC(Ne3Energy_,i,j,k) = Source_VC(Ne3Energy_,i,j,k) &
                 + Qexp_I(4) + Kxp_I(2) + Kxp_I(3)
         else
            Source_VC(Ne3Energy_,i,j,k) = Source_VC(Ne3Energy_,i,j,k) &
                 - Kpx_I(4)
         end if
      end if
      !
      ! Source terms for Plasma
      !/     
      !correcting the source terms (see Gombosi 2000)

      Source_VC(Energy_,i,j,k) = Source_VC(Energy_,i,j,k)         &
           + Qepx_I(2)+ Qepx_I(3) + Qepx_I(4)

    end subroutine user_expl_source
    !==========================================================================
    subroutine user_impl_source

      !\
      ! Source terms for Plasma
      ! writing the source P as SP = [SE -uxS(rhoux)-uyS(rhoUy) - uzS(rhoUz) +0.5u**2 S(rho)]*(gamma-1)
      !/
      ! Source terms for the pressure should be taken from the source terms of the energy

      Source_VC(rhoUx_,i,j,k) = Source_VC(rhoUx_,i,j,k)         &
           + QmpxUx_I(2) + QmpxUx_I(3) + QmpxUx_I(4)
      Source_VC(rhoUy_,i,j,k) = Source_VC(rhoUy_,i,j,k)         &
           + QmpxUy_I(2) + QmpxUy_I(3) + QmpxUy_I(4)
      Source_VC(rhoUz_,i,j,k) = Source_VC(rhoUz_,i,j,k)         & 
           + QmpxUz_I(2) + QmpxUz_I(3) + QmpxUz_I(4)
      Source_VC(P_,i,j,k) = Source_VC(P_,i,j,k)     &
!!$              +(g-1.)*(3./4.)*(I2px_I(2) - I2xp_I(2) &
!!$                           +  I2px_I(3) - I2xp_I(3) &
!!$                           +  I2px_I(4) - I2xp_I(4) )
           + (g-1.0)*(Source_VC(Energy_,i,j,k) &
           - Ux_I(1)*Source_VC(RhoUx_,i,j,k) &
           - Uy_I(1)*Source_VC(RhoUy_,i,j,k) &
           - Uz_I(1)*Source_VC(RhoUz_,i,j,k)) 

! lines for testing

      ! Pop I
      !/
      if(NeutralsI) then
         if (iRegion == 3) then
            Source_VC(NeuRho_,i,j,k) = Source_VC(NeuRho_,i,j,k) &
                 + I0xp_I(3) + I0xp_I(4)
            Source_VC(NeuRhoUx_,i,j,k) = Source_VC(NeuRhoUx_,i,j,k)   &
                 + QmxpUx_I(2) + JxpUx_I(3) + JxpUx_I(4)
            Source_VC(NeuRhoUy_,i,j,k) = Source_VC(NeuRhoUy_,i,j,k)   &
                 + QmxpUy_I(2) + JxpUy_I(3) + JxpUy_I(4)
            Source_VC(NeuRhoUz_,i,j,k) = Source_VC(NeuRhoUz_,i,j,k)   & 
                 + QmxpUz_I(2) + JxpUz_I(3) + JxpUz_I(4)
            Source_VC(NeuP_,i,j,k) = Source_VC(NeuP_,i,j,k)     &
!!$                 + (g-1.)*(3./4.)*( I2xp_I(2)+ I2xp_I(3)+ I2xp_I(4)- I2px_I(2) )           
                 + (g-1.)*(Source_VC(NeuEnergy_,i,j,k) &
                 - Ux_I(2)*Source_VC(NeuRhoUx_,i,j,k) &
                 - Uy_I(2)*Source_VC(NeuRhoUy_,i,j,k) &
                 - Uz_I(2)*Source_VC(NeuRhoUz_,i,j,k) &
                 + 0.5*(UTotal_I(2)**2)*Source_VC(NeuRho_,i,j,k))         
          else
            Source_VC(NeuRho_,i,j,k) = Source_VC(NeuRho_,i,j,k) &
                 - I0px_I(2)
            Source_VC(NeuRhoUx_,i,j,k) = Source_VC(NeuRhoUx_,i,j,k) &
                 - JpxUx_I(2)
            Source_VC(NeuRhoUy_,i,j,k) = Source_VC(NeuRhoUy_,i,j,k) &
                 - JpxUy_I(2)
            Source_VC(NeuRhoUz_,i,j,k) = Source_VC(NeuRhoUz_,i,j,k) &
                 - JpxUz_I(2)
            Source_VC(NeuP_,i,j,k) = Source_VC(NeuP_,i,j,k)     &
!!$                    + (g-1.)*(3./4.)*(-1.)* I2px_I(2)
                 + (g-1.)*(Source_VC(NeuEnergy_,i,j,k) &
                 - Ux_I(2)*Source_VC(NeuRhoUx_,i,j,k) &
                 - Uy_I(2)*Source_VC(NeuRhoUy_,i,j,k) &
                 - Uz_I(2)*Source_VC(NeuRhoUz_,i,j,k) &
                 + 0.5*(UTotal_I(2)**2)*Source_VC(NeuRho_,i,j,k))

         end if
      end if
      !
      ! Pop II
      !/
      if(NeutralsII) then
         if (iRegion == 2) then
               Source_VC(Ne2Rho_,i,j,k) = Source_VC(Ne2Rho_,i,j,k)  &
                    + I0xp_I(2) + I0xp_I(4)
               Source_VC(Ne2RhoUx_,i,j,k) = Source_VC(Ne2RhoUx_,i,j,k)   &
                    + QmxpUx_I(3) + JxpUx_I(2) + JxpUx_I(4)
               Source_VC(Ne2RhoUy_,i,j,k) = Source_VC(Ne2RhoUy_,i,j,k)   &
                    + QmxpUy_I(3) + JxpUy_I(2) + JxpUy_I(4)
               Source_VC(Ne2RhoUz_,i,j,k) = Source_VC(Ne2RhoUz_,i,j,k)   & 
                    + QmxpUz_I(3) + JxpUz_I(2) + JxpUz_I(4)
!!!               Source_VC(Ne2P_,i,j,k) = Source_VC(Ne2P_,i,j,k)     &
!!$                    + (g-1.)*(3./4.)*( I2xp_I(3)+ I2xp_I(2)+  I2xp_I(4)- I2px_I(3) )      
!!!                    + (g-1.)*(Source_VC(Ne2Energy_,i,j,k) &
!!!                    - Ux_I(3)*Source_VC(Ne2RhoUx_,i,j,k) &
!!!                    - Uy_I(3)*Source_VC(Ne2RhoUy_,i,j,k) &
!!!                    - Uz_I(3)*Source_VC(Ne2RhoUz_,i,j,k) &
!!!                    + 0.5*(UTotal_I(3)**2)*Source_VC(Ne2Rho_,i,j,k))
          else
               Source_VC(Ne2Rho_,i,j,k) = Source_VC(Ne2Rho_,i,j,k) &
                    - I0px_I(3) 
               Source_VC(Ne2RhoUx_,i,j,k) = Source_VC(Ne2RhoUx_,i,j,k) &
                    - JpxUx_I(3)
               Source_VC(Ne2RhoUy_,i,j,k) = Source_VC(Ne2RhoUy_,i,j,k) &
                    - JpxUy_I(3)
               Source_VC(Ne2RhoUz_,i,j,k) = Source_VC(Ne2RhoUz_,i,j,k) &
                    - JpxUz_I(3)
!!!               Source_VC(Ne2P_,i,j,k) = Source_VC(Ne2P_,i,j,k)     &
!!$                    + (g-1.)*(3./4.)*(-1.)* I2px_I(3)
!!!                    + (g-1.)*(Source_VC(Ne2Energy_,i,j,k) &
!!!                    - Ux_I(3)*Source_VC(Ne2RhoUx_,i,j,k) &
!!!                    - Uy_I(3)*Source_VC(Ne2RhoUy_,i,j,k) &
!!!                    - Uz_I(3)*Source_VC(Ne2RhoUz_,i,j,k) &
!!!                    + 0.5*(UTotal_I(3)**2)*Source_VC(Ne2Rho_,i,j,k))
          end if     
      end if
      !
      ! Pop III
      !/
      if(NeutralsIII) then
         if (iRegion == 1) then
            Source_VC(Ne3Rho_,i,j,k) = Source_VC(Ne3Rho_,i,j,k) & 
                 + I0xp_I(2) + I0xp_I(3)
            Source_VC(Ne3RhoUx_,i,j,k) = Source_VC(Ne3RhoUx_,i,j,k)   &
                 + QmxpUx_I(4) + JxpUx_I(3) + JxpUx_I(2)
            Source_VC(Ne3RhoUy_,i,j,k) = Source_VC(Ne3RhoUy_,i,j,k)   &
                 + QmxpUy_I(4) + JxpUy_I(3) + JxpUy_I(2)
            Source_VC(Ne3RhoUz_,i,j,k) = Source_VC(Ne3RhoUz_,i,j,k)   & 
                 + QmxpUz_I(4) + JxpUz_I(3) + JxpUz_I(2)
            Source_VC(Ne3P_,i,j,k) = Source_VC(Ne3P_,i,j,k)     &
!!$                 + (g-1.)*(3./4.)*( I2xp_I(4)+ I2xp_I(3)+ I2xp_I(2)- I2px_I(4) )
                 + (g-1.)*(Source_VC(Ne3Energy_,i,j,k) &
                 - Ux_I(4)*Source_VC(Ne3RhoUx_,i,j,k) &
                 - Uy_I(4)*Source_VC(Ne3RhoUy_,i,j,k) &
                 - Uz_I(4)*Source_VC(Ne3RhoUz_,i,j,k) &
                 + 0.5*(UTotal_I(4)**2)*Source_VC(Ne3Rho_,i,j,k))
         else
            Source_VC(Ne3Rho_,i,j,k) = Source_VC(Ne3Rho_,i,j,k) &
                 - I0px_I(4) 
            Source_VC(Ne3RhoUx_,i,j,k) = Source_VC(Ne3RhoUx_,i,j,k)  &
                 - JpxUx_I(4)
            Source_VC(Ne3RhoUy_,i,j,k) = Source_VC(Ne3RhoUy_,i,j,k)  &
                 - JpxUy_I(4)
            Source_VC(Ne3RhoUz_,i,j,k) = Source_VC(Ne3RhoUz_,i,j,k)  &
                 - JpxUz_I(4)
            Source_VC(Ne3P_,i,j,k) = Source_VC(Ne3P_,i,j,k)     &
!!$                 + (g-1.)*(3./4.)*(-1.)* I2px_I(4)
                 + (g-1.)*(Source_VC(Ne3Energy_,i,j,k) &
                 - Ux_I(4)*Source_VC(Ne3RhoUx_,i,j,k) &
                 - Uy_I(4)*Source_VC(Ne3RhoUy_,i,j,k) &
                 - Uz_I(4)*Source_VC(Ne3RhoUz_,i,j,k) &
                 + 0.5*(UTotal_I(4)**2)*Source_VC(Ne3Rho_,i,j,k))
         end if
      end if

    end subroutine user_impl_source

  end subroutine user_calc_sources

  !============================================================================

  subroutine user_init_point_implicit

    use ModVarIndexes, ONLY: Rho_,RhoUx_,RhoUy_,RhoUz_,P_,NeuRho_,NeuRhoUx_,NeuRhoUy_,&

         NeuRhoUz_,NeuP_,Ne2Rho_,Ne2RhoUx_,Ne2RhoUy_,Ne2RhoUz_,Ne2P_,&
         Ne3Rho_,Ne3RhoUx_,Ne3RhoUy_,Ne3RhoUz_,Ne3P_

    use ModPointImplicit, ONLY: iVarPointImpl_I, IsPointImplMatrixSet, EpsPointImpl_V
    !------------------------------------------------------------------------

    ! Allocate and set iVarPointImpl_I
    ! In this example there are 3 implicit variables


    ! All the neutrals momenta and plasma are implicit 
    ! (3 neutral fluid and 1 ion)

    allocate(iVarPointImpl_I(20))

    iVarPointImpl_I = (/Rho_,RhoUx_, RhoUy_, RhoUz_, P_, NeuRho_,NeuRhoUx_ ,&
         NeuRhoUy_, NeuRhoUz_ , NeuP_, Ne2Rho_,Ne2RhoUx_, Ne2RhoUy_, Ne2RhoUz_ ,&
         Ne2P_, Ne3Rho_,Ne3RhoUx_, Ne3RhoUy_, Ne3RhoUz_, Ne3P_/)

    !Because the second and third populations of neutrals have initially very small values I'm
    !setting EpsPointImpl_V to be small for their variables
    EpsPointImpl_V(Ne2Rho_:Ne3P_)=1.e-11

    ! Tell the point implicit scheme if dS/dU will be set analytically
    ! If this is set to true the DsDu_VVC matrix has to be set below.

    IsPointImplMatrixSet = .false.

  end subroutine user_init_point_implicit

end module ModUser
