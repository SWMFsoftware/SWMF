!^CFG COPYRIGHT UM
!==============================================================================

module ModUser
  use ModExpansionFactors
  use ModMagnetogram
  use ModUserEmpty,                                     &
       IMPLEMENTED1  => user_read_inputs,                &
       IMPLEMENTED2  => user_set_ics,                    &
       IMPLEMENTED3  => user_get_b0,                     &
       IMPLEMENTED4  => user_update_states,              &
       IMPLEMENTED5  => user_specify_initial_refinement, &
       IMPLEMENTED6  => user_set_outerbcs,               &
       IMPLEMENTED7  => user_set_boundary_cells,         &
       IMPLEMENTED8  => user_face_bcs,                   &
       IMPLEMENTED9  => user_set_plot_var,               &
       IMPLEMENTED10 => user_get_log_var

  include 'user_module.h' !list of public methods

  real, parameter :: VersionUserModule = 1.1
  character (len=*), parameter :: &
       NameUserModule = 'EMPIRICAL SC on spherical grids'

  ! module variables related to shear flows at the coronal base
  ! ShearTime should be in seconds
  logical :: DoShearFlow
  integer :: iVorticalFlow
  real :: MaxBnARDim, ShearAmplitude, ShearAngle, xShear_D(ndim), &
       RampUpTime, RampDownTime, StartTimeRampDown, ShearTime

  ! module variables for Bipolar Active Region
  logical :: DoBipolarAR
  integer, parameter :: nDipMax=3
  integer :: nDip
  real :: xDip_DI(nDim,nDipMax), MDipUnit_DI(nDim,nDipMax), DipoleDim

  ! module variables for additional refinement window around Active Regions
  logical :: DoARRefinement
  real :: MinDelThetaAR, MaxRadiusAR, MinLattitudeAR, MaxLattitudeAR, &
       MinLongitudeAR, MaxLongitudeAR

  ! module variables for second refinement window
  logical :: DoARRefinement2
  real :: MinDelThetaAR2, MaxRadiusAR2, MinLattitudeAR2, MaxLattitudeAR2, &
       MinLongitudeAR2, MaxLongitudeAR2

contains

  !============================================================================
  subroutine user_read_inputs
    use ModMain
    use ModProcMH,    ONLY: iProc
    use ModReadParam
    use ModIO,        ONLY: write_prefix, write_myname, iUnitOut

    integer:: i
    character (len=100) :: NameCommand
    character (len=lStringLine)   :: NameModel
    !--------------------------------------------------------------------------

    if(iProc==0.and.lVerbose > 0)then
       call write_prefix; write(iUnitOut,*)'User read_input HELIOSPHERE starts'
    endif
    do
       if(.not.read_line() ) EXIT
       if(.not.read_command(NameCommand)) CYCLE
       select case(NameCommand)
       case("#PFSSM")
          call read_var('UseUserB0'  ,UseUserB0)
          if (UseUserB0)then
             call read_magnetogram_file
             call read_var('dt_UpdateB0',dt_UpdateB0)
             DoUpdateB0 = dt_updateb0 > 0.0
          endif
       case("#EMPIRICALSW")
          call read_var('NameModel',NameModel)
          call set_empirical_model(trim(NameModel))
       case('#USERINPUTEND')
          if(iProc==0.and.lVerbose > 0)then
             call write_prefix;
             write(iUnitOut,*)'User read_input HELIOSPHERE ends'
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

       call init_AR_parameters
       call init_AR_refinement_window
       call init_shear_profile
    end do
  end subroutine user_read_inputs

  !============================================================================
  subroutine init_AR_parameters
    use ModNumConst, ONLY: cDegToRad
    implicit none

    integer :: iDip
    real :: rDip, PhiFirstDip, DelPhiDip, Theta, Alpha, Phi
    real :: CosPhi, SinPhi, CosTheta, SinTheta, CosAlpha, SinAlpha
    real :: x, y, z
    !--------------------------------------------------------------------------
    DoBipolarAR = .true.

    nDip = 3
    DipoleDim = -0.003
    rDip = 0.94
    PhiFirstDip = -3.0
    DelPhiDip = 3.0

    ! anti-clockwise angle Alpha in yz-plane between dipole moment
    ! and z unit vector
    Alpha = 0.0

    do iDip=1,nDip
       Theta = 90.0
       Phi = PhiFirstDip + (iDip-1)*DelPhiDip

       SinTheta = sin(Theta*cDegToRad)
       CosTheta = cos(Theta*cDegToRad)
       SinPhi = sin(Phi*cDegToRad)
       CosPhi = cos(Phi*cDegToRad)
       CosAlpha = cos(Alpha*cDegToRad)
       SinAlpha = sin(Alpha*cDegToRad)

       x = rDip*CosPhi*SinTheta
       y = rDip*SinPhi*SinTheta
       z = rDip*CosTheta

       xDip_DI(x_,iDip) = x
       xDip_DI(y_,iDip) = y*CosAlpha - z*SinAlpha
       xDip_DI(z_,iDip) = y*SinAlpha + z*CosAlpha

       MDipUnit_DI(x_,iDip) = 0.0
       MDipUnit_DI(y_,iDip) =-SinAlpha
       MDipUnit_DI(z_,iDip) = CosAlpha
    end do

  end subroutine init_AR_parameters

  !============================================================================
  subroutine init_shear_profile
    use ModNumConst, ONLY: cDegToRad
    implicit none

    real :: ShearLattitude, ShearLongitude
    real :: SinLattitude, CosLattitude, SinLongitude, CosLongitude
    !--------------------------------------------------------------------------
    DoShearFlow = .true.

    ! 1 = Hydrodynamical vortical flow
    ! 2 = Magnetic flux conserving voritical flow
    iVorticalFlow = 1

    MaxBnARDim = 20.0 ! [Gauss]
    ShearAmplitude = 1.088e-2   ! 2.0 percent U_A
    ShearAngle = 10.0*cDegToRad
    RampUpTime = 2400.0 ! [Seconds]
    RampDownTime = 2400.0 ! [Seconds]
    StartTimeRampDown = 20000.0 ! [Seconds]
    ShearTime = 10000.0 ! [Seconds]

    ShearLattitude = 0.0  ! [Degrees]
    ShearLongitude = 0.0  ! [Degrees]

    SinLattitude = sin(ShearLattitude*cDegToRad)
    CosLattitude = cos(ShearLattitude*cDegToRad)
    SinLongitude = sin(ShearLongitude*cDegToRad)
    CosLongitude = cos(ShearLongitude*cDegToRad)

    xShear_D(1) = CosLattitude*CosLongitude
    xShear_D(2) = CosLattitude*SinLongitude
    xShear_D(3) = SinLattitude
    
  end subroutine init_shear_profile

  !============================================================================
  subroutine init_AR_refinement_window
    use ModNumConst, ONLY: cDegToRad
    implicit none
    !--------------------------------------------------------------------------
    DoARRefinement = .true.

    MinDelThetaAR = 0.25
    MaxRadiusAR = 1.5
    MinLattitudeAR =-6.2
    MaxLattitudeAR = 6.2
    MinLongitudeAR =-14.6
    MaxLongitudeAR = 14.6

    MinDelThetaAR = MinDelThetaAR*cDegToRad
    MinLattitudeAR = MinLattitudeAR*cDegToRad
    MaxLattitudeAR = MaxLattitudeAR*cDegToRad
    MinLongitudeAR = MinLongitudeAR*cDegToRad
    MaxLongitudeAR = MaxLongitudeAR*cDegToRad

    DoARRefinement2 = .false.

    MinDelThetaAR2 = 0.125
    MaxRadiusAR2 = 1.3
    MinLattitudeAR2 =-3.4
    MaxLattitudeAR2 = 3.4
    MinLongitudeAR2 =-14.25
    MaxLongitudeAR2 = 14.25

    MinDelThetaAR2 = MinDelThetaAR2*cDegToRad
    MinLattitudeAR2 = MinLattitudeAR2*cDegToRad
    MaxLattitudeAR2 = MaxLattitudeAR2*cDegToRad
    MinLongitudeAR2 = MinLongitudeAR2*cDegToRad
    MaxLongitudeAR2 = MaxLongitudeAR2*cDegToRad

  end subroutine init_AR_refinement_window

  !============================================================================
  subroutine user_set_boundary_cells(iBLK)
    use ModGeometry,      ONLY: ExtraBc_, IsBoundaryCell_GI, r_Blk
    use ModBoundaryCells, ONLY: SaveBoundaryCells
    use ModPhysics,       ONLY: rBody
    implicit none

    integer, intent(in):: iBLK

    character (len=*), parameter :: Name='user_set_boundary_cells'
    !--------------------------------------------------------------------------
    IsBoundaryCell_GI(:,:,:,ExtraBc_) = r_Blk(:,:,:,iBLK) < rBody

    if(SaveBoundaryCells) return
    call stop_mpi('Set SaveBoundaryCells=.true. in PARAM.in file')

  end subroutine user_set_boundary_cells

  !============================================================================
  subroutine user_face_bcs(VarsGhostFace_V)
    use ModFaceBc,     ONLY: FaceCoords_D, VarsTrueFace_V, B0Face_D, TimeBc, &
         iBlockBc, jFace, kFace
    use ModMain,       ONLY: time_accurate,x_,y_,z_, UseRotatingFrame, &
         Time_Simulation
    use ModNumConst,   ONLY: cTolerance
    use ModPhysics,    ONLY: inv_gm1, OmegaBody
    use ModSize,       ONLY: nDim
    use ModVarIndexes, ONLY: nVar,Ew_,rho_,Ux_,Uy_,Uz_,Bx_,By_,Bz_,P_
    implicit none

    real, intent(out) :: VarsGhostFace_V(nVar)

    real :: Dens, Pres, Gamma, r, Un, B1n, FullBn, Bnorm
    real, dimension(nDim) :: UnitR_D, U_D, UparB_D, &
         B1_D, B1n_D, B1t_D, FullB_D, Bunit_D
    real :: UTheta, UPhi, Ux, Uy, Uz
    !--------------------------------------------------------------------------
    r = sqrt(sum(FaceCoords_D**2))
    UnitR_D  = FaceCoords_D/r

    U_D   = VarsTrueFace_V(Ux_:Uz_)
    Un    = dot_product(UnitR_D,U_D)
    B1_D  = VarsTrueFace_V(Bx_:Bz_)
    B1n_D = dot_product(UnitR_D,B1_D)*UnitR_D
    B1t_D = B1_D - B1n_D

    !\
    ! Update BCs for induction field
    !/
    VarsGhostFace_V(Bx_:Bz_) = B1t_D

    !\
    ! Update BCs for velocity field
    !/
    FullB_D = B0Face_D + VarsGhostFace_V(Bx_:Bz_)
    FullBn = dot_product(UnitR_D,FullB_D)

    if(abs(FullBn) < cTolerance .or. Un < 0.0)then
       ! No flow at polarity inversion lines and no backflow
       VarsGhostFace_V(Ux_:Uz_) = -Un*UnitR_D
    else
       Bnorm    = sqrt(dot_product(FullB_D,FullB_D))
       Bunit_D  = FullB_D/Bnorm
       UparB_D  = dot_product(VarsTrueFace_V(Ux_:Uz_),Bunit_D)*Bunit_D
       VarsGhostFace_V(Ux_:Uz_) = UparB_D
    end if

!!!    if(DoShearFlow .and. Time_Accurate .and. Time_Simulation>0.0 &
!!!        .and. Time_Simulation<StartTimeRampDown+RampDownTime)then
!!!    if(DoShearFlow .and. Time_Simulation<ShearTime .and. Time_Accurate)then
    if(DoShearFlow .and. Time_Accurate)then
       call get_shear_flow(UTheta,UPhi,Ux,Uy,Uz,jFace,kFace,iBlockBc,FullBn)
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) + Ux
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) + Uy
       VarsGhostFace_V(Uz_) = VarsGhostFace_V(Uz_) + Uz
    end if

    !\
    ! Update BCs for the mass density, pressure, Wave Energy
    !/
    call plasma_at_base(FaceCoords_D(x_),FaceCoords_D(y_), &
         FaceCoords_D(z_),r,Dens,Pres)
    call my_gamma_emp(FaceCoords_D(x_),FaceCoords_D(y_), &
         FaceCoords_D(z_),Gamma)

    VarsGhostFace_V(rho_) = 2.0*Dens - VarsTrueFace_V(rho_)
    VarsGhostFace_V(P_) = VarsGhostFace_V(rho_)*Pres/Dens
    VarsGhostFace_V(Ew_) = VarsGhostFace_V(rho_) &
         *Pres/Dens*(1.0/(Gamma-cOne)-inv_gm1)

    !\
    ! Apply corotation
    !/
    if (.not.UseRotatingFrame) then
       VarsGhostFace_V(Ux_) = VarsGhostFace_V(Ux_) &
            - 2.0*OmegaBody*FaceCoords_D(y_)
       VarsGhostFace_V(Uy_) = VarsGhostFace_V(Uy_) &
            + 2.0*OmegaBody*FaceCoords_D(x_)
    end if

  end subroutine user_face_bcs

  !============================================================================
  subroutine get_shear_flow(UTheta,UPhi,Ux,Uy,Uz,j,k,iBlock,FullBn)
    use ModGeometry, ONLY: Dy_Blk, Dz_Blk, XyzStart_Blk
    use ModNumConst, ONLY: cTolerance, cHalfPi
    use ModPhysics,  ONLY: rBody
    implicit none

    real, intent(out) :: UTheta, UPhi, Ux, Uy, Uz
    integer, intent(in) :: j, k, iBlock
    real, intent(in) :: FullBn

    real :: dTheta, dPhi, r, Lattitude, Theta, Phi
    real :: SinTheta, CosTheta, SinPhi, CosPhi
    real :: ShearProfileL, ShearProfileR, FullBnL, FullBnR
    !--------------------------------------------------------------------------

    dTheta = Dz_Blk(iBlock)
    dPhi = Dy_Blk(iBlock)

    r = rBody
    Lattitude = XyzStart_Blk(3,iBlock) + real(k-1)*dTheta
    ! Theta is the co-lattitude
    Theta = cHalfPi - Lattitude
    Phi = XyzStart_Blk(2,iBlock) + real(j-1)*dPhi

    SinTheta = sin(Theta)
    CosTheta = cos(Theta)
    SinPhi = sin(Phi)
    CosPhi = cos(Phi)

    select case(iVorticalFlow)
    case(1)
       ShearProfileL = shear_profile_HD(Phi-0.5*dPhi,Theta,iBlock)
       ShearProfileR = shear_profile_HD(Phi+0.5*dPhi,Theta,iBlock)
       UTheta = 1.0/(r*SinTheta) &
            *(ShearProfileR-ShearProfileL)/dPhi

       ShearProfileL = shear_profile_HD(Phi,Theta-0.5*dTheta,iBlock)
       ShearProfileR = shear_profile_HD(Phi,Theta+0.5*dTheta,iBlock)
       UPhi = -1.0/(r) &
            *(ShearProfileR-ShearProfileL)/dTheta
    case(2)
       if(abs(FullBn)<cTolerance)then
          ! No flow at polarity inversion lines
          UTheta = 0.0
          UPhi = 0.0
          Ux = 0.0
          Uy = 0.0
          Uz = 0.0
          return
       else
          ShearProfileL = shear_profile_FixBr(Phi-0.5*dPhi,Theta,iBlock,FullBnL)
          ShearProfileR = shear_profile_FixBr(Phi+0.5*dPhi,Theta,iBlock,FullBnR)
          if(FullBnL*FullBnR<0.0 .and. abs(FullBnL)>cTolerance &
               .and. abs(FullBnR)>cTolerance)then
             UTheta = 0.0
          else
             UTheta = 1.0/(r*SinTheta*FullBn) &
                  *(ShearProfileR-ShearProfileL)/dPhi
          end if

          ShearProfileL = shear_profile_FixBr(Phi,Theta-0.5*dTheta,iBlock,FullBnL)
          ShearProfileR = shear_profile_FixBr(Phi,Theta+0.5*dTheta,iBlock,FullBnR)
          if(FullBnL*FullBnR<0.0 .and. abs(FullBnL)>cTolerance &
               .and. abs(FullBnR)>cTolerance)then
             UPhi = 0.0
          else
             UPhi = -1.0/(r*FullBn) &
                  *(ShearProfileR-ShearProfileL)/dTheta
          end if
       end if
    end select

    Ux = UTheta*CosTheta*CosPhi - UPhi*SinPhi
    Uy = UTheta*CosTheta*SinPhi + UPhi*CosPhi
    Uz =-UTheta*SinTheta

  end subroutine get_shear_flow

  !============================================================================
  real function shear_profile_HD(Phi,Theta,iBlock)
    use ModMain,     ONLY: Time_Simulation
    use ModNumConst, ONLY: cHalfPi, cPi, cTwoPi, cDegToRad
    use ModPhysics,  ONLY: Io2No_V, UnitB_, rBody
    use ModVarIndexes
    implicit none

    real, intent(in) :: Phi, Theta
    integer, intent(in) :: iBlock

    real :: Lattitude, Longitude
    real :: DelTheta1, DelTheta2, DelTheta3, DelPhi1, DelPhi2
    real :: Mask, Ramp
    !--------------------------------------------------------------------------

    DelTheta1 = 2.5*cDegToRad
    DelTheta2 = 3.5*cDegToRad
    DelTheta3 = 13.5*cDegToRad
    Lattitude = cHalfPi - Theta

    if(abs(Lattitude)<DelTheta1)then
       Mask = sin(cHalfPi*Lattitude/DelTheta1)**2
    else if(abs(Lattitude)<DelTheta2)then
       Mask = 1.0
    else if(abs(Lattitude)<DelTheta3)then
       Mask = cos(cHalfPi*(abs(Lattitude)-DelTheta2)/(DelTheta3-DelTheta2))**2
    else
       shear_profile_HD = 0.0
       return
    end if

    DelPhi1 = 6.0*cDegToRad
    DelPhi2 = 12.0*cDegToRad

    if( (Phi<DelPhi2 .and. Phi>-DelPhi2) .or.&
        (Phi<-cTwoPi+DelPhi2 .and. Phi>-cTwoPi-DelPhi2) .or.&
        (Phi< cTwoPi+DelPhi2 .and. Phi> cTwoPi-DelPhi2))then

       Longitude = Phi
       if(Phi>-DelPhi2 .and. Phi>DelPhi2) Longitude = Phi - cTwoPi
       if(Phi<-DelPhi2 .and. Phi<DelPhi2) Longitude = Phi + cTwoPi

       if(abs(Longitude)>DelPhi1)then
          Mask = Mask*cos(cHalfPi*(abs(Longitude)-DelPhi1) &
               /(DelPhi2-DelPhi1))**2
       end if
    else
       shear_profile_HD = 0.0
       return
    end if

!!!    Ramp = 1.0

!!!    Ramp = sin(cPi*Time_Simulation/ShearTime)

!!!    if(Time_Simulation<StartTimeRampDown)then
       Ramp = min(Time_Simulation/RampUpTime,1.0)
!!!    else
!!!       Ramp = max(1.0-(Time_Simulation-StartTimeRampDown)/RampDownTime,0.0)
!!!    end if

    shear_profile_HD = ShearAmplitude*Mask*Ramp
       
  end function shear_profile_HD

  !============================================================================
  real function shear_profile_FixBr(Phi,Theta,iBlock,FullBn)
    use ModMain,     ONLY: Time_Simulation
    use ModNumConst, ONLY: cHalfPi, cPi, cTwoPi, cDegToRad
    use ModPhysics,  ONLY: Io2No_V, UnitB_, rBody
    use ModVarIndexes
    implicit none

    real, intent(in) :: Phi, Theta
    integer, intent(in) :: iBlock
    real, intent(out) :: FullBn

    real :: x, y, z, r, Angle
    real :: B0_D(nDim), UnitR_D(nDim), Del_D(nDim)
    real :: Mask, Ramp, MaxBnAR, BnRatio
    !--------------------------------------------------------------------------
    r = rBody

    x = r*sin(Theta)*cos(Phi)
    y = r*sin(Theta)*sin(Phi)
    z = r*cos(Theta)

    call user_get_b0(x,y,z,B0_D)

    UnitR_D(x_) = x/r
    UnitR_D(y_) = y/r
    UnitR_D(z_) = z/r

    FullBn = dot_product(UnitR_D,B0_D)

    Del_D = xShear_D - UnitR_D
    Angle = 2.0*asin(0.5*sqrt(dot_product(Del_D,Del_D)))

    Ramp = 1.0

!!!    Ramp = sin(cPi*Time_Simulation/ShearTime)

!!!    if(Time_Simulation<StartTimeRampDown)then
!!!       Ramp = min(Time_Simulation/RampUpTime,1.0)
!!!    else
!!!       Ramp = max(1.0-(Time_Simulation-StartTimeRampDown)/RampDownTime,0.0)
!!!    end if

    MaxBnAR = MaxBnARDim*Io2No_V(UnitB_)
    BnRatio = FullBn/MaxBnAR

    shear_profile_FixBr = ShearAmplitude*FullBn**3*exp(1.0-BnRatio**2) &
         *Mask*Ramp

  end function shear_profile_FixBr

  !============================================================================
  subroutine user_set_outerbcs(iBlock, iSide, TypeBc, IsFound)
    use ModAdvance,  ONLY: State_VGB,B0xFace_x_Blk,B0yFace_x_Blk, &
         B0zFace_x_Blk,Hyp_
    use ModMain,     ONLY: UseHyperbolicDivb, UseRotatingFrame, &
         Time_Simulation, Time_Accurate
    use ModNumConst, ONLY: cTolerance, cPi
    use ModGeometry, ONLY: x_Blk, y_Blk, z_Blk, r_Blk
    use ModPhysics,  ONLY: inv_gm1, OmegaBody
    use ModSize
    use ModVarIndexes
    use ModFaceValue, ONLY: UseLogRhoLimiter
    implicit none

    integer,          intent(in)  :: iBlock, iSide
    character(len=20),intent(in)  :: TypeBc
    logical,          intent(out) :: IsFound

    character (len=*), parameter :: Name='user_set_outerbcs'

    integer :: i, j, k, iMin, jMin, kMin, iMax, jMax, kMax
    real :: B1_D(nDim), B1n, B1t_D(nDim), Bnorm, Bunit_D(nDim), &
         FullB_D(nDim), FullBn
    real :: U_D(nDim), Un, UparB_D(nDim), Umesh_D(nDim)
    real, dimension(nDim) :: UnitR_D
    real :: x, y, z, r, Dens, Pres, Gamma
    real :: UTheta, UPhi, Ux, Uy, Uz
    !--------------------------------------------------------------------------
    IsFound = .true.

    if(iSide/=East_) call stop_mpi('Wrong iSide in user_set_outerBCs')

    iMin=1-gcn;iMax=0
    jMin=1-gcn;jMax=nJ+gcn
    kMin=1-gcn;kMax=nK+gcn

    do k=kMin,kMax; do j=jMin,jMax
       x = 0.5*sum(x_Blk(0:1,j,k,iBlock))
       y = 0.5*sum(y_Blk(0:1,j,k,iBlock))
       z = 0.5*sum(z_Blk(0:1,j,k,iBlock))
       r = 0.5*sum(r_Blk(0:1,j,k,iBlock))

       call plasma_at_base(x,y,z,r,Dens,Pres)
       call my_gamma_emp(x,y,z,Gamma)

       if(UseLogRhoLimiter)then
          State_VGB(Rho_,1,j,k,iBlock)=log(State_VGB(Rho_,1,j,k,iBlock))
          Dens = log(Dens)
       end if
       State_VGB(Rho_,iMax,j,k,iBlock) = 2.0*Dens &
            - State_VGB(Rho_,iMax+1,j,k,iBlock)
       do i=iMax-1,iMin,-1
          State_VGB(Rho_,i,j,k,iBlock) = 2.0*State_VGB(Rho_,i+1,j,k,iBlock) &
               - State_VGB(Rho_,i+2,j,k,iBlock)
       end do
       if(UseLogRhoLimiter)then
          State_VGB(Rho_,iMin:1,j,k,iBlock) = &
               exp(State_VGB(Rho_,iMin:1,j,k,iBlock))
          Dens = exp(Dens)
       end if

       do i=iMin,iMax
          State_VGB(P_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
               *Pres/Dens
          State_VGB(Ew_,i,j,k,iBlock) = State_VGB(Rho_,i,j,k,iBlock) &
               *Pres/Dens*(1.0/(Gamma-1.0)-inv_gm1)
       end do

       UnitR_D(x_)=x_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       UnitR_D(y_)=y_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)
       UnitR_D(z_)=z_Blk(1,j,k,iBlock)/r_Blk(1,j,k,iBlock)


       ! B1_D and U_D should be the reconstructed variables
       B1_D  = State_VGB(Bx_:Bz_,1,j,k,iBlock)
       B1n   = dot_product(UnitR_D,B1_D)
       B1t_D = B1_D - B1n*UnitR_D

       do i=iMin,iMax
          State_VGB(Bx_:Bz_,i,j,k,iBlock) = B1t_D
       end do

       U_D = State_VGB(RhoUx_:RhoUz_,1,j,k,iBlock)/State_VGB(Rho_,1,j,k,iBlock)
       Un = dot_product(UnitR_D,U_D)

       FullB_D(x_) = B0xFace_x_BLK(1,j,k,iBlock) + B1t_D(x_)
       FullB_D(y_) = B0yFace_x_BLK(1,j,k,iBlock) + B1t_D(y_)
       FullB_D(z_) = B0zFace_x_BLK(1,j,k,iBlock) + B1t_D(z_)
       FullBn = dot_product(UnitR_D,FullB_D)

       if(abs(FullBn) < cTolerance .or. Un < 0.0)then
          ! No flow at polarity inversion lines and no backflow
          do i=iMin,iMax
             Umesh_D = State_VGB(RhoUx_:RhoUz_,1-i,j,k,iBlock) &
                  /State_VGB(Rho_,1-i,j,k,iBlock)
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                  -State_VGB(Rho_,i,j,k,iBlock) &
                  *dot_product(UnitR_D,Umesh_D)*UnitR_D
          end do
       else
          Bnorm   = sqrt(dot_product(FullB_D,FullB_D))
          Bunit_D = FullB_D/Bnorm
          UparB_D = dot_product(U_D,Bunit_D)*Bunit_D
          do i=iMin,iMax
             State_VGB(RhoUx_:RhoUz_,i,j,k,iBlock) = &
                  State_VGB(Rho_,i,j,k,iBlock)*UparB_D
          end do
       end if

!!!       if(DoShearFlow .and. Time_Accurate .and. Time_Simulation>0.0 &
!!!           .and. Time_Simulation<StartTimeRampDown+RampDownTime)then
!!!       if(DoShearFlow .and. Time_Simulation<ShearTime .and. Time_Accurate)then
       if(DoShearFlow .and. Time_Accurate)then
          call get_shear_flow(UTheta,UPhi,Ux,Uy,Uz,j,k,iBlock,FullBn)
          do i=iMin,iMax
             State_VGB(RhoUx_,i,j,k,iBlock) = &
                  State_VGB(RhoUx_,i,j,k,iBlock) &
                  + State_VGB(Rho_,i,j,k,iBlock)*Ux
             State_VGB(RhoUy_,i,j,k,iBlock) = &
                  State_VGB(RhoUy_,i,j,k,iBlock) &
                  + State_VGB(Rho_,i,j,k,iBlock)*Uy
             State_VGB(RhoUz_,i,j,k,iBlock) = &
                  State_VGB(RhoUz_,i,j,k,iBlock) &
                  + State_VGB(Rho_,i,j,k,iBlock)*Uz
          end do
       end if

    end do; end do

    !\
    ! Apply corotation
    !/
    if(.not.UseRotatingFrame)then
       State_VGB(RhoUx_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = &
            State_VGB(RhoUx_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            -2.0*OmegaBody*y_Blk(iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            *State_VGB(Rho_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock)
       State_VGB(RhoUy_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = &
            State_VGB(RhoUy_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            -2.0*OmegaBody*x_Blk(iMin:iMax,jMin:jMax,kMin:kMax,iBlock) &
            *State_VGB(Rho_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock)
    end if

    if(UseHyperbolicDivb)then
       State_VGB(Hyp_,iMin:iMax,jMin:jMax,kMin:kMax,iBlock) = 0.0
    end if

  end subroutine user_set_outerbcs

  !============================================================================
  subroutine user_set_ics
    use ModMain,      ONLY: globalBLK, nI, nJ, nK, UseHyperbolicDivb
    use ModVarIndexes
    use ModAdvance,   ONLY: State_VGB, Hyp_
    use ModPhysics,   ONLY: inv_gm1, GBody, Si2No_V, UnitX_
    use ModGeometry,  ONLY: x_Blk, y_Blk, z_Blk, r_Blk
    implicit none

    integer :: iBLK, i, j, k
    integer :: IterCount
    logical :: oktest,oktest_me
    real :: x, y, z, r, Dens, Pres, Gamma
    real :: Ur, Ur0, Ur1, del, UBase, rTransonic, UEscape

    real, parameter :: Epsilon = 1.0e-6
    !--------------------------------------------------------------------------
    call set_oktest('user_set_ics',oktest,oktest_me)

    iBLK = globalBLK

    UEscape = sqrt(-GBody*2.0/(rSun*Si2No_V(UnitX_)))

    !\
    ! Initialize MHD wind with Parker's solution
    ! construct solution which obeys
    !   rho x u_r x r^2 = constant
    !/
    rTransonic = 0.25*UEscape**2
    if(.not.(rTransonic>exp(1.0))) call stop_mpi('sonic point inside Sun')

    UBase=rTransonic**2*exp(1.5-2.0*rTransonic)

    do k=1,nK;do j=1,nJ; do i=1,nI
       x = x_BLK(i,j,k,iBLK)
       y = y_BLK(i,j,k,iBLK)
       z = z_BLK(i,j,k,iBLK)
       r = r_BLK(i,j,k,iBLK)

       if(r > rTransonic)then
          !\
          ! Inside supersonic region
          !/
          Ur0 = 1.0
          IterCount = 0
          do
             IterCount = IterCount + 1
             Ur1 = sqrt(UEscape**2/r - 3.0 &
                  + 2.0*log(16.0*Ur0*r**2/UEscape**4))
             del = abs(Ur1 - Ur0)
             if(del < Epsilon)then
                Ur = Ur1
                EXIT
             elseif(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end do
       else
          !\
          ! Inside subsonic region
          !/
          Ur0 = 1.0
          IterCount = 0
          do
             IterCount = IterCount + 1
             Ur1 = (UEscape**2/(4.0*r))**2 &
                  *exp(0.5*(Ur0**2 + 3.0 - UEscape**2/r))
             del = abs(Ur1 - Ur0)
             if(del < Epsilon)then
                Ur = Ur1
                EXIT
             elseif(IterCount < 1000)then
                Ur0 = Ur1
                CYCLE
             else
                call stop_mpi('PARKER > 1000 it.')
             end if
          end do
       end if

       call plasma_at_base(x,y,z,r,Dens,Pres)
       State_VGB(Rho_,i,j,k,iBLK) = Dens*UBase/(r**2*Ur)
       State_VGB(RhoUx_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK)*Ur*x/r
       State_VGB(RhoUy_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK)*Ur*y/r
       State_VGB(RhoUz_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK)*Ur*z/r
       State_VGB(Bx_:Bz_,i,j,k,iBLK) = 0.0
       if(UseHyperbolicDivb)then
          State_VGB(Hyp_,i,j,k,iBLK) = 0.0
       end if
       State_VGB(P_,i,j,k,iBLK) = State_VGB(Rho_,i,j,k,iBLK)*Pres/Dens
       call my_gamma_emp(x,y,z,Gamma)
       State_VGB(Ew_,i,j,k,iBLK) = State_VGB(P_,i,j,k,iBLK) &
            *(1.0/(Gamma-1.0)-inv_gm1) 
    end do;end do; end do

  end subroutine user_set_ics

  !============================================================================
  subroutine plasma_at_base(x,y,z,r,Dens,Pres)
    use ModPhysics,  ONLY: BodyRho_I, Si2No_V, UnitTemperature_
    implicit none

    real, intent(in) :: x, y, z, r
    real, intent(out) :: Dens, Pres

    real :: UFinal, URatio
    !--------------------------------------------------------------------------
!    call get_bernoulli_integral(x/r,y/r,z/r,UFinal)
!    URatio = UFinal/UMin
!    Dens  = BodyRho_I(1)*(1.0/URatio)**2
!    Pres  = Dens*T0*Si2No_V(UnitTemperature_)/min(URatio,2.0)
    Dens  = BodyRho_I(1)
    Pres  = Dens*T0*Si2No_V(UnitTemperature_)

  end subroutine plasma_at_base

  !============================================================================
  subroutine user_get_b0(xInput,yInput,zInput,B0_D)
    use ModPhysics,  ONLY: Io2No_V,UnitB_
    implicit none

    real, intent(in):: xInput,yInput,zInput
    real, intent(out), dimension(3):: B0_D
    !--------------------------------------------------------------------------

    call get_magnetogram_field(xInput,yInput,zInput,B0_D)

    ! The additional Bipolar B0 is not used for determining the varying gamma
    if(DoBipolarAR) call get_AR_field(xInput,yInput,zInput,B0_D)

    B0_D = B0_D*Io2No_V(UnitB_)

  end subroutine user_get_b0

  !============================================================================
  subroutine get_AR_field(x,y,z,B0_D)
    use ModNumConst, ONLY: cTolerance
    implicit none

    real, intent(in) :: x, y, z
    real, intent(inout) :: B0_D(nDim)

    integer :: iDip
    real :: xShift_D(nDim), r, r_inv, r2_inv, r3_inv
    real :: Dp, Bdp_D(nDim)
    !--------------------------------------------------------------------------
    do iDip=1,nDip
       xShift_D(1) = x - xDip_DI(1,iDip)
       xShift_D(2) = y - xDip_DI(2,iDip)
       xShift_D(3) = z - xDip_DI(3,iDip)

       r = max(sqrt(dot_product(xShift_D,xShift_D)),cTolerance)
       r_inv = 1.0/r
       r2_inv = r_inv*r_inv
       r3_inv = r_inv*r2_inv

       Bdp_D = DipoleDim*MDipUnit_DI(:,iDip)

       Dp = dot_product(Bdp_D,xShift_D)*3.0*r2_inv

       B0_D = B0_D + (Dp*xShift_D-Bdp_D)*r3_inv
    end do

  end subroutine get_AR_field

  !============================================================================
  subroutine user_update_states(iStage,iBlock)
    use ModVarIndexes
    use ModSize
    use ModAdvance, ONLY: State_VGB
    use ModMain,    ONLY: nStage
    use ModPhysics, ONLY: inv_gm1
    use ModGeometry,ONLY: x_BLK, y_BLK, z_BLK, R_BLK
    use ModEnergy,  ONLY: calc_energy_cell
    use ModExpansionFactors, ONLY:  GammaSS
    implicit none

    integer, intent(in) :: iStage, iBlock

    integer :: i, j, k
    real    :: Gammacell
    !--------------------------------------------------------------------------
    call update_states_MHD(iStage,iBlock)
    !\
    ! Begin update of pressure and relaxation energy::
    !/
    !  if (iStage/=nStage) return
    do k=1,nK; do j=1,nJ; do i=1,nI
       call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
            z_BLK(i,j,k,iBlock),Gammacell)
       call correct_gamma(i,j,k,iBlock,Gammacell)

       State_VGB(P_,i,j,k,iBlock) = (Gammacell-1.0) &
            *(inv_gm1*State_VGB(P_,i,j,k,iBlock) + State_VGB(Ew_,i,j,k,iBlock))
       State_VGB(Ew_,i,j,k,iBlock) = State_VGB(P_,i,j,k,iBlock) &
            *(1.0/(Gammacell-1.0)-inv_gm1)
    end do; end do; end do
    call calc_energy_cell(iBlock)
    !\
    ! End update of pressure and relaxation energy::
    !/
  end subroutine user_update_states

  !============================================================================
  subroutine correct_gamma(i,j,k,iBlock,Gamma)
    use ModAdvance, ONLY: State_VGB, B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModGeometry,ONLY: R_BLK
    use ModVarIndexes
    implicit none

    integer, intent(in) :: i, j, k, iBlock
    real, intent(inout) :: Gamma
    !--------------------------------------------------------------------------
    !\
    ! Change the empirical gamma near the current sheet based on the
    ! plasma beta
    !/
    if(R_BLK(i,j,k,iBlock) > Rs_PFSSM) &
         Gamma = Gamma - (Gamma - GammaSS)*max(0.0, &
         -1.0 + 2.0*State_VGB(P_,i,j,k,iBlock)/(State_VGB(P_,i,j,k,iBlock) &
         +((State_VGB(Bx_,i,j,k,iBlock) + B0xCell_BLK(i,j,k,iBlock))**2 &
         + (State_VGB(By_,i,j,k,iBlock) + B0yCell_BLK(i,j,k,iBlock))**2 &
         + (State_VGB(Bz_,i,j,k,iBlock) + B0zCell_BLK(i,j,k,iBlock))**2) &
         *0.25*(R_BLK(i,j,k,iBlock)/Rs_PFSSM)**1.5))

  end subroutine correct_gamma

  !============================================================================
  subroutine user_specify_initial_refinement(iBLK,refineBlock,lev,DxBlock, &
       xCenter,yCenter,zCenter,rCenter,                        &
       minx,miny,minz,minR,maxx,maxy,maxz,maxR,found)
    use ModMain,ONLY:time_loop,nI,nJ,nK
    use ModAMR,ONLY:InitialRefineType, initial_refine_levels
    use ModNumConst, ONLY: cTwoPi, cDegToRad
    use ModAdvance,ONLY:&
         State_VGB,Bx_,By_,Bz_,B0xCell_BLK,B0yCell_BLK,B0zCell_BLK
    use ModGeometry
    use ModPhysics,ONLY:rBody
    use ModOctree
    implicit none

    logical,intent(out) :: refineBlock, found
    integer, intent(in) :: lev
    real, intent(in)    :: DxBlock
    real, intent(in)    :: xCenter,yCenter,zCenter,rCenter
    real, intent(in)    :: minx,miny,minz,minR
    real, intent(in)    :: maxx,maxy,maxz,maxR
    integer, intent(in) :: iBLK

    character (len=*), parameter :: Name='user_specify_initial_refinement'
    real :: BDotRMin,BDotRMax,CritR
    integer :: i,j,k, levmin, levbodyfocus
    logical :: IsSouth, IsNorth
    real :: MinThetaBlk, MaxThetaBlk, MinPhiBlk, MaxPhiBlk
    !--------------------------------------------------------------------------
    ! Do not use minx,miny,minz,maxx,maxy,maxz, since they are not set.

    select case (InitialRefineType)
    case ('helio_init')

       select case(TypeGeometry)
       case('spherical_lnr')
          ! assumes 1x2x1 blocks on root level
          levmin=4
          levbodyfocus=5
          if(.not.time_loop)then
             !refine to have resolution not worse than levmin
             if(lev<=levmin)then
                RefineBlock = .true.
             elseif(lev<=levbodyfocus)then
                ! Additional bodyfocus refinement
                CritR=rBody + (exp(XyzMax_D(1))-exp(XyzMin_D(1))) &
                     /real(2**(lev+2-levmin))

                ! But no additional refinement near poles
                IsSouth = XyzStart_BLK(3,iBlk) - Dz_Blk(iBlk) &
                     < XyzMin_D(3) + real(nK*(2**(lev-levmin)-2))*Dz_Blk(iBlk)
                IsNorth = XyzStart_BLK(3,iBlk) + real(nK)*Dz_Blk(iBlk) &
                     > XyzMax_D(3) - real(nK*(2**(lev-levmin)-2))*Dz_Blk(iBlk)

!!!                if( .not.(IsSouth.or.IsNorth) )then
                if( MinR < CritR .and. .not.(IsSouth.or.IsNorth) )then
                   RefineBlock = .true.
                else
                   RefineBlock = .false.
                endif
             else
                RefineBlock = .false.
             end if
          else
             MinPhiBlk = XyzStart_BLK(2,iBlk) - 0.5*Dy_Blk(iBlk)
             MaxPhiBlk = MinPhiBlk + real(nJ)*Dy_Blk(iBlk)
             MinThetaBlk = XyzStart_BLK(3,iBlk) - 0.5*Dz_Blk(iBlk)
             MaxThetaBlk = MinThetaBlk + real(nK)*Dz_Blk(iBlk)
             if(DoARRefinement &
                  .and. MinR < MaxRadiusAR &
                  .and. MaxThetaBlk > MinLattitudeAR &
                  .and. MinThetaBlk < MaxLattitudeAR &
                  .and. ((MaxPhiBlk > MinLongitudeAR &
                  .and.   MinPhiBlk < MaxLongitudeAR) &
                  .or.   (MaxPhiBlk > MinLongitudeAR-cTwoPi &
                  .and.   MinPhiBlk < MaxLongitudeAR-cTwoPi) &
                  .or.   (MaxPhiBlk > MinLongitudeAR+cTwoPi &
                  .and.   MinPhiBlk < MaxLongitudeAR+cTwoPi)) &
                  .and. Dz_Blk(iBlk) > MinDelThetaAR)then
                RefineBlock = .true.
             else if(DoARRefinement2 &
                  .and. MinR < MaxRadiusAR2 &
                  .and. MaxThetaBlk > MinLattitudeAR2 &
                  .and. MinThetaBlk < MaxLattitudeAR2 &
                  .and. ((MaxPhiBlk > MinLongitudeAR2 &
                  .and.   MinPhiBlk < MaxLongitudeAR2) &
                  .or.   (MaxPhiBlk > MinLongitudeAR2-cTwoPi &
                  .and.   MinPhiBlk < MaxLongitudeAR2-cTwoPi) &
                  .or.   (MaxPhiBlk > MinLongitudeAR2+cTwoPi &
                  .and.   MinPhiBlk < MaxLongitudeAR2+cTwoPi)) &
                  .and. Dz_Blk(iBlk) > MinDelThetaAR2)then
                RefineBlock = .true.
             elseif( Dz_Blk(iBlk) > 2.0*cDegToRad )then
                call refine_heliosheath
             else
                RefineBlock = .false.
             end if
          end if
       case default
          call stop_mpi('user_specify_initial_refinement is ' &
               //'not implemented for geometry= '//TypeGeometry)
       end select
       found=.true.
    end select

    contains
      subroutine refine_heliosheath
        !----------------------------------------------------------------------

        BDotRMin=0.0
        do k=0,nK+1;do j=1,nJ
           BDotRMin=min( BDotRMin,minval(&
                (B0xCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(Bx_,1:nI,j,k,iBLK))*&
                x_BLK(1:nI,j,k,iBLK)+&
                (B0yCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(By_,1:nI,j,k,iBLK))*&
                y_BLK(1:nI,j,k,iBLK)+&
                (B0zCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(Bz_,1:nI,j,k,iBLK))*&
                z_BLK(1:nI,j,k,iBLK)))
        end do;end do
        BDotRMax=0.0
        do k=0,nK+1;do j=1,nJ
           BDotRMax=max( BDotRMax,maxval(&
                (B0xCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(Bx_,1:nI,j,k,iBLK))*&
                x_BLK(1:nI,j,k,iBLK)+&
                (B0yCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(By_,1:nI,j,k,iBLK))*&
                y_BLK(1:nI,j,k,iBLK)+&
                (B0zCell_BLK(1:nI,j,k,iBLK)+&
                State_VGB(Bz_,1:nI,j,k,iBLK))*&
                z_BLK(1:nI,j,k,iBLK)))
        end do;end do
        refineBlock =BDotRMin<-cTiny.and.&
             BDotRMax>cTiny
      end subroutine refine_heliosheath

  end subroutine user_specify_initial_refinement

  !============================================================================
  subroutine user_set_plot_var(iBlock, NameVar, IsDimensional, &
       PlotVar_G, PlotVarBody, UsePlotVarBody, &
       NameTecVar, NameTecUnit, NameIdlUnit, IsFound)

    use ModAdvance, ONLY: State_VGB
    use ModSize, ONLY: nI, nJ, nK, gcn
    use ModGeometry, ONLY: x_BLK, y_BLK, z_BLK, r_BLK, Dy_Blk, Dz_Blk, &
         XyzStart_Blk
    use ModNumConst, ONLY: cHalfPi
    use ModPhysics, ONLY: No2Io_V, UnitU_, UnitB_, UnitX_, NameTecUnit_V, rBody
    use ModVarIndexes
    use ModExpansionFactors, ONLY:  GammaSS
    implicit none

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

    integer :: i,j,k
    real :: Gamma
    real :: x, y, z, r, Lattitude, Theta, Phi, UnitR_D(nDim)
    real :: B0_D(nDim), FullBn
    real :: UTheta, UPhi, Ux, Uy, Uz
    !--------------------------------------------------------------------------

    IsFound=.true.
    
    PlotVar_G = 0.0

    select case(NameVar)
    case('g0')
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock),Gamma)

          PlotVar_G(i,j,k) = Gamma
       end do; end do; end do
    case('g1')
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          call my_gamma_emp(x_BLK(i,j,k,iBlock),y_BLK(i,j,k,iBlock), &
               z_BLK(i,j,k,iBlock),Gamma)
          call correct_gamma(i,j,k,iBlock,Gamma)

          PlotVar_G(i,j,k) = Gamma
       end do; end do; end do
    case('utheta','uphi','b0r')
       r = rBody
       do k=1-gcn,nK+gcn; do j=1-gcn,nJ+gcn; do i=1-gcn,nI+gcn
          Lattitude = XyzStart_Blk(3,iBlock) + real(k-1)*Dz_Blk(iBlock)
          ! Theta is the co-lattitude
          Theta = cHalfPi - Lattitude
          Phi = XyzStart_Blk(2,iBlock) + real(j-1)*Dy_Blk(iBlock)
          x = r*sin(Theta)*cos(Phi)
          y = r*sin(Theta)*sin(Phi)
          z = r*cos(Theta)
          call user_get_b0(x,y,z,B0_D)
          UnitR_D(x_) = x/r
          UnitR_D(y_) = y/r
          UnitR_D(z_) = z/r
          FullBn = dot_product(UnitR_D,B0_D)
          select case(NameVar)
          case('utheta')
             call get_shear_flow(UTheta,UPhi,Ux,Uy,Uz,j,k,iBlock,FullBn)
             PlotVar_G(i,j,k) = UTheta
             if(IsDimensional) PlotVar_G(i,j,k) = &
                  PlotVar_G(i,j,k)*No2Io_V(UnitU_)
             NameTecUnit = NameTecUnit_V(UnitU_)
          case('uphi')
             call get_shear_flow(UTheta,UPhi,Ux,Uy,Uz,j,k,iBlock,FullBn)
             PlotVar_G(i,j,k) = UPhi
             if(IsDimensional) PlotVar_G(i,j,k) = &
                  PlotVar_G(i,j,k)*No2Io_V(UnitU_)
             NameTecUnit = NameTecUnit_V(UnitU_)
          case('b0r')
             PlotVar_G(i,j,k) = FullBn
             if(IsDimensional) PlotVar_G(i,j,k) = &
                  PlotVar_G(i,j,k)*No2Io_V(UnitB_)
             NameTecUnit = NameTecUnit_V(UnitB_)
          end select
       end do; end do; end do
       NameTecVar = NameVar
    end select

    UsePlotVarBody=.true.
    PlotVarBody=1.0

  end subroutine user_set_plot_var

  !============================================================================
  subroutine user_get_log_var(VarValue, TypeVar, Radius)
    use ModIO,         ONLY: write_myname
    use ModMain,       ONLY: unusedBLK, nBLK
    use ModVarIndexes, ONLY: Ew_,Bx_,By_,Bz_,rho_,rhoUx_,rhoUy_,rhoUz_,P_ 
    use ModAdvance,    ONLY: State_VGB, tmp1_BLK, &
         B0xCell_BLK, B0yCell_BLK, B0zCell_BLK
    use ModPhysics,    ONLY: inv_gm1, No2Io_V, UnitEnergyDens_, UnitX_
    implicit none

    real, intent(out)            :: VarValue
    character (len=*), intent(in):: TypeVar
    real, intent(in), optional :: Radius

    character (len=*), parameter :: Name='user_get_log_var'
    integer :: iBLK
    real :: UnitEnergy
    real, external :: integrate_BLK
    !--------------------------------------------------------------------------

    UnitEnergy = No2Io_V(UnitEnergyDens_)*No2Io_V(UnitX_)**3

    select case(TypeVar)
    case('em_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               (B0xcell_BLK(:,:,:,iBLK)+State_VGB(Bx_,:,:,:,iBLK))**2 &
               + (B0ycell_BLK(:,:,:,iBLK)+State_VGB(By_,:,:,:,iBLK))**2 &
               + (B0zcell_BLK(:,:,:,iBLK)+State_VGB(Bz_,:,:,:,iBLK))**2
       end do
       VarValue = UnitEnergy*0.5*integrate_BLK(1,tmp1_BLK)
    case('ek_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = &
               ( State_VGB(rhoUx_,:,:,:,iBLK)**2 &
               + State_VGB(rhoUy_,:,:,:,iBLK)**2 &
               + State_VGB(rhoUz_,:,:,:,iBLK)**2) &
               /State_VGB(rho_,:,:,:,iBLK)             
       end do
       VarValue = UnitEnergy*0.5*integrate_BLK(1,tmp1_BLK)
    case('et_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(P_,:,:,:,iBLK)
       end do
       VarValue = UnitEnergy*inv_gm1*integrate_BLK(1,tmp1_BLK)
    case('ew_t')
       do iBLK=1,nBLK
          if (unusedBLK(iBLK)) CYCLE
          tmp1_BLK(:,:,:,iBLK) = State_VGB(Ew_,:,:,:,iBLK)
       end do
       VarValue = UnitEnergy*integrate_BLK(1,tmp1_BLK)
    case default
       VarValue = -7777.
       call write_myname;
       write(*,*) 'Warning in set_user_logvar: unknown logvarname = ',TypeVar
    end select

  end subroutine user_get_log_var

  !============================================================================
  subroutine my_gamma_emp(xx,yy,zz,gammaOut)

    ! Subroutine my_gamma_emp
    ! Provides the distribution of the polytropic index, complying with
    ! the WSA or Fisk semi-empirical models

    use ModExpansionFactors
    use ModNumConst
    implicit none

    real, intent(in) :: xx,yy,zz
    real, intent(out)   :: gammaOut 
    real :: RR,Uf,BernoulliFactor
    real, parameter :: gammaIH=1.5
    real, parameter :: R1=2.50,R2=12.50
    integer,parameter::nPowerIndex=2
    !--------------------------------------------------------------------------
    !\
    ! Calculate cell-centered spherical coordinates::
    RR   = sqrt(xx**2+yy**2+zz**2)
    !\
    ! Avoid calculating inside a critical radius = 0.5*Rsun
    !/
    if (RR <max(Ro_PFSSM-dR*nRExt,0.90*Ro_PFSSM)) then 
       gammaOut= gammaSS
       RETURN
    end if

    ! Calculate gamma
    if(RR >= R2)then
       gammaOut=gammaIH
    else if(RR >= R1)then
       gammaOut=gammaSS+(RR-R1)*(gammaIH-gammaSS)/(R2-R1)
    else
       call get_bernoulli_integral(xx,yy,zz,Uf)
       BernoulliFactor=(cHalf*Uf**2+cSunGravitySI)/&
            (T0*cBoltzmann/cProtonMass)&
            *(R1-RR)*&
            & (Ro_PFSSM/RR)**nPowerIndex/ (R1-Ro_PFSSM)+ GammaSS&
            &/(GammaSS-cOne)*(cOne- (R1-RR)*(Ro_PFSSM/RR)&
            &**nPowerIndex/ (R1-Ro_PFSSM))
       gammaOut = BernoulliFactor/(BernoulliFactor-cOne)
    end if

  end subroutine my_gamma_emp

end module ModUser
